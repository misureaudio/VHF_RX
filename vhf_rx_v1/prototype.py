"""VHF NBFM receiver prototype, v1 (see docs/RSPdx_VHF_NBFM_Receiver_Design.md).

Front end: SDRplay RSPdx or RSP2/RSP2pro — both are supported and
co-exist in this codebase; the connected unit is probed at startup and
selected with --radio (default: auto). Per-device sample scale and LNA
limits come from sdrplay_api.PROFILES.

Examples
--------
  python prototype.py --band 2m --duration 10 --auto-ch 3
  python prototype.py --band 70cm --duration 10 --ch 434.600e6 433.500e6
  python prototype.py --radio rsp2pro --band 2m --ch 145.500e6 --play --ctcss 88.5 --duration 30
  python prototype.py --band 2m --duration 10 --gr 40 --lna 20   # quieter front end (RSPdx states)
"""
from __future__ import annotations

import argparse
import json
import queue
import sys
import threading
import time

import numpy as np

from nbfm import BandSpectrum, ChannelDemod, CTCSS_TONES, Squelch
from sdrplay_api import (AGC_50HZ, BW_5M, BW_8M, PROFILES, SDRPlay, SDRPlayError)

FS = 6.048e6          # max rate at full bit depth: 14-bit on RSPdx, 16-bit
                      # on RSP2/RSP2pro (see sdrplay_api.PROFILES)
CHUNK = 1 << 16       # 65536 samples (~10.8 ms) per DSP pass: amortizes the
                      # per-call overhead of the decimation chain ~60x vs
                      # the ~1 k-sample service blocks

BANDS = {
    "2m":   dict(rf=145.0e6, bw=BW_5M, fmin=144.0e6, fmax=146.0e6),
    "70cm": dict(rf=434.0e6, bw=BW_8M, fmin=432.0e6, fmax=436.0e6),
}


def band_lna_max(profile, band: str) -> int:
    """Highest valid LNAstate for the band on the selected radio."""
    return profile.lna_max_2m if band == "2m" else profile.lna_max_70cm


class RxWorker:
    """Consumes IQ blocks, runs the DSP chain, publishes state."""

    def __init__(self, args, band: dict, sample_scale: float, device_info: dict):
        self.args = args
        self.band = band
        self.fs = FS
        # int16 full scale of the connected radio (8192 RSPdx / 32768
        # RSP2-RSP2pro) — a wrong value shifts every dBFS quantity by 12 dB
        self.sample_scale = sample_scale
        self.device_info = device_info
        # unbounded: the service delivers in bursts (USB batching) that
        # far exceed any fixed depth; the worker's capacity (8.5 MSPS)
        # exceeds the average feed rate, so the queue stays bounded in
        # practice (a full burst ≈ 5 MB)
        self.q: "queue.Queue" = queue.Queue(maxsize=0)
        self.drops = 0
        self.n_samples = 0
        self.n_blocks = 0
        self.t_start = 0.0
        self.t_stop = 0.0
        self.total_gain_db = float("nan")
        self.overloads = 0

        # one FFT per window serves both the channel grid and the waterfall
        self.spectrum = BandSpectrum(self.fs, band["rf"], band["fmin"], band["fmax"])
        self.demods: dict[float, ChannelDemod] = {}
        self.squelches: dict[float, Squelch] = {}
        self.audio_buf: dict[float, list] = {}
        self.metrics: dict[float, dict] = {}
        self.active_now: list[float] = []
        self.player = None          # playback.AudioPlayer, set by main() with --play
        self._lock = threading.Lock()
        self._stop = threading.Event()

    # -- callback side (service thread) --------------------------------
    def on_iq(self, xi, xq, n, reset):
        # keep this as thin as possible: view the raw int16 samples and
        # queue owned copies; all conversion/DSP happens on the worker
        # thread. (as_array works on the ctypes pointer the service passes;
        # the view aliases the service buffer, so .copy() is mandatory.)
        try:
            i = np.ctypeslib.as_array(xi, shape=(n,)).copy()
            q = np.ctypeslib.as_array(xq, shape=(n,)).copy()
            self.q.put_nowait((i, q))
        except queue.Full:
            self.drops += 1
        except Exception:
            import traceback; traceback.print_exc()

    def on_event(self, event_id, tuner, gain_db):
        # event 0 = gain change (gain in dB is valid); other events
        # (e.g. 1 = power overload) carry no meaningful gain value
        if event_id == 0 and gain_db == gain_db:  # not nan
            self.total_gain_db = gain_db
        if event_id == 1:       # power overload
            self.overloads += 1

    # -- worker side ----------------------------------------------------
    def _ensure_demod(self, f: float, ctcss_hz: float | None = None):
        if f not in self.demods:
            d = ChannelDemod(f, self.fs)
            d.set_ctcss(ctcss_hz)
            self.demods[f] = d
            self.squelches[f] = Squelch()
            self.audio_buf[f] = []
            self.metrics[f] = {"carrier_db": float("-inf"), "ctcss_tone": None,
                               "squelch_open": False, "open_samples": 0}

    def _pick_channels(self, cp, chans, noise_db):
        explicit = self.args.channels
        if explicit:
            return [f for f in explicit if self.band["fmin"] <= f <= self.band["fmax"]]
        if cp is None:
            return self.active_now
        # the decimation chain only reaches channels within ~±300 kHz of
        # the LO: the /7 stage's anti-alias FIR cuts at 0.40·fs/7 =
        # 345.6 kHz @ 6.048 MSPS, so a channel further out is filtered
        # before demod (its 16 kHz signal is pure noise floor — verified
        # live: a +18 dB signal at 144.000 with LO 145.000 demodulates to
        # nothing). The waterfall still shows it; we just don't demod it.
        win = 300e3
        inwin = (chans >= self.band["rf"] - win) & (chans <= self.band["rf"] + win)
        thr = noise_db + self.args.squelch_floor_db
        db = 10.0 * np.log10(np.maximum(cp, 1e-30))
        idx = np.where((db >= thr) & inwin)[0]
        if len(idx) == 0:
            return []
        order = idx[np.argsort(db[idx])[::-1]]
        top = [float(chans[i]) for i in order[:self.args.auto_ch]]
        # hysteresis: keep currently active channels above a lower bar
        keep = [f for f in self.active_now
                if inwin[int(np.argmin(np.abs(chans - f)))]
                and (10.0 * np.log10(max(float(cp[int(np.argmin(np.abs(chans - f)))]), 1e-30))
                    >= noise_db + max(6.0, self.args.squelch_floor_db - 6.0))]
        merged = list(dict.fromkeys(top + keep))
        return merged[:max(self.args.auto_ch, len(keep))]

    def _process_chunk(self, chunk: np.ndarray):
        """Run the full DSP pass on one ~64 k-sample chunk."""
        cp, chans, noise_db = self.spectrum.update(chunk)
        active = self._pick_channels(cp, chans, noise_db)
        self.active_now = active
        for f in list(self.demods):
            if f not in active and f not in (self.args.channels or []):
                del self.demods[f]
                del self.squelches[f]
                del self.audio_buf[f]
                del self.metrics[f]
                if self.player is not None:
                    self.player.drop(f)
        for f in active:
            ctcss = self.args.ctcss if f in (self.args.channels or []) else None
            self._ensure_demod(f, ctcss)
            audio, m = self.demods[f].process(chunk)
            sq = self.squelches[f].update(m["carrier_db"], audio.size * 2)
            mt = self.metrics[f]
            mt["carrier_db"] = m["carrier_db"]
            mt["ctcss_tone"] = m["ctcss_tone"]
            mt["squelch_open"] = sq
            if sq:
                mt["open_samples"] += audio.size
            if audio.size:
                buf = self.audio_buf[f]
                buf.append(audio)
                if sum(a.size for a in buf) > self.args.duration * 8000 * 2:
                    del buf[:len(buf) // 2]
            # play tap: non-destructive copy (audio_buf above feeds the WAV
            # save; the old --play code popped from it and raced save())
            if self.player is not None and audio.size:
                self.player.push(f, audio, sq)
        if self.player is not None:
            self.player.flush()

    def run(self):
        self.t_start = time.monotonic()
        pending: list[np.ndarray] = []
        pending_n = 0
        t_wait = t_conv = t_proc = 0.0
        t_last = time.monotonic()
        while not self._stop.is_set():
            try:
                i16, q16 = self.q.get(timeout=0.5)
            except queue.Empty:
                t_wait += time.monotonic() - t_last
                t_last = time.monotonic()
                continue
            t_wait += time.monotonic() - t_last
            t0 = time.monotonic()
            self.n_blocks += 1
            iq = ((i16.astype(np.float32) + 1j * q16.astype(np.float32))
                  / self.sample_scale).astype(np.complex64)
            self.n_samples += iq.size
            t1 = time.monotonic()
            pending.append(iq)
            pending_n += iq.size
            while pending_n >= CHUNK:
                allc = np.concatenate(pending)
                self._process_chunk(allc[:CHUNK])
                rest = allc[CHUNK:]
                pending = [rest] if rest.size else []
                pending_n = rest.size
            t_proc += time.monotonic() - t1
            t_conv += t1 - t0
            t_last = time.monotonic()
            if t_last - self.t_start >= self.args.duration:
                break
        if pending:
            self._process_chunk(np.concatenate(pending))
        self.t_stop = time.monotonic()
        wall = max(self.t_stop - self.t_start, 1e-9)
        print(f"[worker-timing] wall {wall:.1f}s  "
              f"queue-wait {t_wait / wall * 100:5.1f}%  "
              f"convert {t_conv / wall * 100:5.1f}%  "
              f"process {t_proc / wall * 100:5.1f}%")

    def stop(self):
        self._stop.set()

    # -- report ----------------------------------------------------------
    def report(self) -> dict:
        rep = {
            "band": self.args.band,
            "rf_hz": self.band["rf"],
            "fs_hz": self.fs,
            **self.device_info,
            "lna": self.band["lna"],
            "gr_db": self.args.gr,
            "agc": "50Hz",
            "agc_setpoint_dbfs": self.args.agc_setpoint,
            "total_gain_db": self.total_gain_db,
            "duration_s": self.t_stop - self.t_start,
            "blocks": self.n_blocks,
            "samples": self.n_samples,
            "measured_rate_mspS": self.n_samples / max(self.t_stop - self.t_start, 1e-9) / 1e6,
            "dropped_blocks": self.drops,
            "overload_events": self.overloads,
            "noise_floor_db": self.spectrum.noise_db,
            "channels": {},
        }
        for f, mt in self.metrics.items():
            cp_idx = int(np.argmin(np.abs(self.spectrum.channels - f)))
            p_db = 10.0 * np.log10(max(float(self.spectrum.ch_power_last[cp_idx]), 1e-30)) \
                if getattr(self.spectrum, "ch_power_last", None) is not None else None
            rep["channels"][f"{f/1e6:.4f}"] = {
                "power_db": p_db,
                "snr_over_floor_db": (p_db - self.spectrum.noise_db) if p_db is not None else None,
                "carrier_db_last": mt["carrier_db"],
                "squelch_open": mt["squelch_open"],
                "open_fraction": mt["open_samples"] / max(self.args.duration * 8000, 1),
                "ctcss_tone": mt["ctcss_tone"],
            }
        return rep

    def save(self, outdir: str) -> list[str]:
        import os
        os.makedirs(outdir, exist_ok=True)
        saved = []
        tag = f"{self.args.band}_{int(self.band['rf']/1e6)}"
        wf = f"{outdir}/waterfall_{tag}.png"
        self.spectrum.to_png(wf)
        saved.append(wf)
        for f, buf in self.audio_buf.items():
            if not buf:
                continue
            audio = np.concatenate(buf)
            peak = float(np.max(np.abs(audio))) if audio.size else 0.0
            if peak > 0:
                audio = audio / peak * 0.9
            wav = f"{outdir}/audio_{f/1e6:.4f}.wav"
            _write_wav(wav, audio.astype(np.float32), 8000)
            saved.append(wav)
        rep_path = f"{outdir}/report_{tag}.json"
        with open(rep_path, "w") as fh:
            json.dump(self.report(), fh, indent=2)
        saved.append(rep_path)
        return saved


def _write_wav(path, audio, fs):
    from scipy.io import wavfile
    pcm = np.clip(audio * 32767.0, -32768, 32767).astype(np.int16)
    wavfile.write(path, fs, pcm)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--band", choices=list(BANDS), default="2m")
    ap.add_argument("--radio", default="auto", choices=["auto", "rspdx", "rsp2pro"],
                    help="SDRplay front end: auto = probe (sole device is used; "
                         "if both are connected the RSPdx is preferred), rspdx, rsp2pro")
    ap.add_argument("--freq", type=float, default=None, help="LO frequency Hz (overrides band default)")
    ap.add_argument("--fs", type=float, default=FS)
    ap.add_argument("--gr", type=int, default=50, help="gain reduction dB (20..59)")
    ap.add_argument("--lna", type=int, default=None,
                    help="LNA state (default: band max for the selected radio; "
                         "RSPdx 0-26/0-20, RSP2pro 0-8/0-5)")
    ap.add_argument("--agc-setpoint", type=int, default=-20, help="AGC setpoint dBFS")
    ap.add_argument("--ch", dest="channels", type=float, nargs="+", default=None,
                    help="fixed channel list (Hz); default: auto-select")
    ap.add_argument("--auto-ch", type=int, default=3, help="max auto-selected channels")
    ap.add_argument("--squelch-floor-dB", dest="squelch_floor_db", type=float, default=12.0)
    ap.add_argument("--ctcss", type=float, default=None,
                    help="CTCSS tone (Hz) for fixed channels; must be in the TS 103 236 table")
    ap.add_argument("--duration", type=float, default=10.0)
    ap.add_argument("--play", action="store_true",
                    help="play audio on the local output device (48 kHz mono, sounddevice)")
    ap.add_argument("--vol", type=float, default=1.0,
                    help="playback volume multiplier on the fixed 1/5000 gain (3 kHz dev = 0.6 FS at 1.0)")
    ap.add_argument("--device", default=None,
                    help="output device index or name substring (default: system default)")
    ap.add_argument("--out", default="output")
    args = ap.parse_args()

    band = dict(BANDS[args.band])
    if args.freq:
        band["rf"] = args.freq
    if args.ctcss is not None:
        near = min(CTCSS_TONES, key=lambda t: abs(t - args.ctcss))
        if abs(near - args.ctcss) > 0.01 * near:
            ap.error(f"CTCSS {args.ctcss} Hz not in ETSI TS 103 236 table (nearest {near})")
        args.ctcss = near
    for f in (args.channels or []):
        if abs(f - band["rf"]) > 300e3:
            print(f"[rx] WARNING: {f/1e6:.4f} MHz is >300 kHz from the LO "
                  f"({band['rf']/1e6:.3f} MHz) — outside the demod chain's reach "
                  f"(h1 anti-alias cutoff 345.6 kHz); only the noise floor will be "
                  f"demodulated. Use --freq near the channel.")

    rx = SDRPlay()
    dev = rx.find_device(args.radio)
    prof = rx.profile
    print(f"[rx] API v{rx.api_version:.2f}, device {rx._ser(dev)} (hwVer {dev.hwVer}) "
          f"-> {prof.label}, int16 full scale {prof.sample_scale:.0f}")
    if args.lna is None:
        args.lna = band_lna_max(prof, args.band)
    elif not (0 <= args.lna <= band_lna_max(prof, args.band)):
        ap.error(f"LNA state {args.lna} out of range for {prof.label} "
                 f"band {args.band} (0..{band_lna_max(prof, args.band)})")
    band["lna"] = args.lna
    device_info = {
        "radio": prof.key,
        "radio_label": prof.label,
        "serial": rx._ser(dev),
        "hwver": int(dev.hwVer),
        "sample_scale": prof.sample_scale,
    }
    rx.get_params()
    rx.configure(fs_hz=args.fs, rf_hz=band["rf"], bw_khz=band["bw"],
                 lna=args.lna, gr_db=args.gr, agc=AGC_50HZ,
                 agc_setpoint_dbfs=args.agc_setpoint)
    worker = RxWorker(args, band, prof.sample_scale, device_info)
    rx.start(worker.on_iq, worker.on_event)
    print(f"[rx] streaming: LO {band['rf']/1e6:.3f} MHz, fs {args.fs/1e6:.3f} MSPS, "
          f"IF BW {band['bw']/1e3:.0f} MHz, LNA {args.lna}, GR {args.gr} dB")

    th = threading.Thread(target=worker.run, daemon=True)
    th.start()

    player = None
    if args.play:
        from playback import AudioPlayer
        player = AudioPlayer(vol=args.vol)
        dev = (int(args.device) if args.device is not None
               and str(args.device).isdigit() else args.device)
        try:
            player.start(device=dev)
            import sounddevice as sd
            q = sd.query_devices(player.stream.device, "output")
            print(f"[rx] audio: 48 kHz mono -> {q['name']}")
        except Exception as e:
            print(f"[rx] WARNING: audio output unavailable "
                  f"({type(e).__name__}: {e}); recording only")
            player = None
    worker.player = player

    try:
        th.join(timeout=args.duration + 5)
    except KeyboardInterrupt:
        pass
    worker.stop()
    time.sleep(0.3)
    if player:
        player.stop()
    saved = worker.save(args.out)
    rx.stop()

    rep = worker.report()
    print("\n================ REPORT ================")
    print(f"band={rep['band']}  LO={rep['rf_hz']/1e6:.3f} MHz  fs={rep['fs_hz']/1e6:.3f} MSPS "
          f"(measured {rep['measured_rate_mspS']:.3f})  blocks={rep['blocks']} drops={rep['dropped_blocks']} "
          f"overloads={rep['overload_events']}")
    print(f"total gain={rep['total_gain_db']:.1f} dB  noise floor={rep['noise_floor_db']:.1f} dBFS")
    for f, c in rep["channels"].items():
        snr = f"{c['snr_over_floor_db']:.1f} dB" if c["snr_over_floor_db"] is not None else "n/a"
        print(f"  {f} MHz: P={c['power_db']:.1f} dBFS  SNR={snr}  "
              f"squelch={'OPEN' if c['squelch_open'] else 'closed'} "
              f"(open {c['open_fraction']*100:.0f}%)  ctcss={c['ctcss_tone']}")
    if player is not None:
        und_ms = player.underrun_frames / player.fs_out * 1000.0
        print(f"audio: {und_ms:.0f} ms underrun (starved while a channel "
              f"was open) over {rep['duration_s']:.0f} s "
              f"({'ok' if und_ms < 50 else 'check CPU load / SDR gaps'})")
    print("saved:", *saved, sep="\n  ")


if __name__ == "__main__":
    main()
