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
from sdrplay_api import (AGC_50HZ, BW_1_536M, BW_5M, BW_8M, PROFILES, SDRPlay,
                         SDRPlayError)

FS = 6.048e6          # max rate at full bit depth: 14-bit on RSPdx, 16-bit
                      # on RSP2/RSP2pro (see sdrplay_api.PROFILES)
CHUNK = 1 << 16       # 65536 samples (~10.8 ms) per DSP pass: amortizes the
                      # per-call overhead of the decimation chain ~60x vs
                      # the ~1 k-sample service blocks

BANDS = {
    "2m":   dict(rf=145.0e6, bw=BW_5M, fmin=144.0e6, fmax=146.0e6),
    "70cm": dict(rf=434.0e6, bw=BW_8M, fmin=432.0e6, fmax=436.0e6),
    # wideband FM broadcast (EU 87.5-108 MHz; user range 88-108 MHz) — Phase 2
    # band plumbing. rf = band centre, bw = IF BW, fmin/fmax = band edges,
    # spacing = channel grid (100 kHz EU), dev = design deviation (50 kHz EU),
    # tau = de-emphasis (50 us EU broadcast), audio_fs = final audio rate (48 kHz).
    # fs = 2.016 MSPS (RSPdx 144 kHz grid, = 48 kHz x 7 x 2 x 3): the LO is
    # retuned onto the acquired channel (Phase 3), so the demod only needs
    # ~+/-124 kHz around DC — the 6.048 MSPS rate (8 MHz IF) was pure overhead.
    # Measured on the dev machine: the DSP chain at 6.048 MSPS needs ~0.95 of
    # one core -> the worker runs at 3.6 MSPS (59.7 % of real time) and the
    # audio ring underruns 40 % of a 300 s run. At 2.016 MSPS the chain costs
    # ~1/3 (full-rate DDC, h1 FIR and BPF all scale with fs) -> no underruns.
    # bw = BW_1_536M: the largest IF BW that fits the 2.016 MSPS span
    # (+/-0.768 MHz); main() swaps to BW_8M if --fs is overridden higher.
    "wbfm": dict(rf=98.0e6, bw=BW_1_536M, fmin=88.0e6, fmax=108.0e6,
                 spacing=100e3, dev=50e3, tau=50e-6, audio_fs=48000,
                 wideband=True, fs=2.016e6, half_bw=50e3),
}


def band_lna_max(profile, band: str) -> int:
    """Highest valid LNAstate for the band on the selected radio."""
    if band in ("2m", "wbfm"):
        # 88-108 MHz broadcast sits in the VHF LNA bucket (below 144 MHz);
        # 144-146 MHz (2 m) uses the same bucket.
        return profile.lna_max_2m
    if band == "70cm":
        return profile.lna_max_70cm
    return profile.lna_max_70cm


class RxWorker:
    """Consumes IQ blocks, runs the DSP chain, publishes state."""

    def __init__(self, args, band: dict, sample_scale: float, device_info: dict,
                 retune=None):
        self.args = args
        self.band = band
        self.retune = retune      # callable(rf_hz) -> None : LO retune (Phase 3)
        # per-band rate: wideband runs at 2.016 MSPS (see BANDS), NBFM at FS
        self.fs = band.get("fs", FS)
        # int16 full scale of the connected radio (8192 RSPdx / 32768
        # RSP2-RSP2pro) — a wrong value shifts every dBFS quantity by 12 dB
        self.sample_scale = sample_scale
        self.device_info = device_info
        # unbounded: the service delivers in bursts (USB batching) that
        # far exceed any fixed depth; the worker's capacity must exceed the
        # feed rate or the queue grows without bound (measured: at 6.048 MSPS
        # wideband the worker only sustained ~3.6 MSPS -> 40 % underruns;
        # the 2.016 MSPS wideband default leaves ~3x headroom)
        self.q: "queue.Queue" = queue.Queue(maxsize=0)
        self.drops = 0
        self.n_samples = 0
        self.n_blocks = 0
        self.t_start = 0.0
        self.t_stop = 0.0
        self.total_gain_db = float("nan")
        self.overloads = 0

        # one FFT per window serves both the channel grid and the waterfall
        # band is clamped to the resolvable span (f_center +/- min(fs/2, bw/2)):
        # the nfft=8192 window only resolves +/- fs/2 around f_center, and the
        # IF filter only passes +/- bw/2 — a wider band (e.g. the 88-108 MHz
        # broadcast band) would put channel edges outside the FFT -> empty bin
        # range -> IndexError; at 2.016 MSPS wideband the IF BW (1.536 MHz)
        # is the tighter bound, so the waterfall spans +/-0.768 MHz.
        span_half = min(self.fs / 2, band["bw"] * 1e3 / 2)
        span_lo = band["rf"] - span_half
        span_hi = band["rf"] + span_half
        # half_bw: channel-power window. NBFM: +/-6.25 kHz (half the 12.5 kHz
        # channel). Wideband: +/-50 kHz (the dev span) — with the narrow NBFM
        # window a heavily modulated broadcast station (strong low-frequency
        # program content, e.g. a 74 Hz bass/hum) spreads its carrier power
        # into sidebands far outside +/-6.25 kHz (J0(beta) ~ 0 for beta >> 1),
        # the measured channel power dips >6 dB and auto-select flaps
        # select/drop every ~0.5 s (verified live at 103.0 MHz). The +/-50 kHz
        # window measures the occupied span, which is ~constant in total power.
        self.spectrum = BandSpectrum(self.fs, band["rf"],
                                     max(band["fmin"], span_lo), min(band["fmax"], span_hi),
                                     spacing=band.get("spacing", 12.5e3),
                                     half_bw=band.get("half_bw", 6.25e3))
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
            if str(self.args.band) == "wbfm":
                # wideband FM broadcast: re-scaled chain from Phase 1
                # (complex rate fs/14 — 432 kHz @ 6.048 MSPS, 144 kHz @ 2.016 MSPS;
                # BPF +/-(dev+fa+margin), LP 15 kHz, final audio 48 kHz via
                # /3 @ 2.016 or /9 @ 6.048, de-emph tau=50 us for EU broadcast).
                # rf_center=f: the LO is retuned onto the channel (Phase 3), so the
                # SDR centres the IQ on f and the carrier sits at DC -> the DDC
                # offset is (f - rf_center) = 0. Without this the DDC downconverts
                # by the absolute f (aliasing ~184 kHz off), moving the carrier out
                # of the BPF -> silent/crackling audio.
                # pll_f_n: 25 kHz at 6.048 MSPS (fs_ch = 432 kHz, w_n = 0.36
                # rad/smp): the NBFM loop (1.75 kHz) low-passes the audio at
                # its loop bandwidth (PLL demod output = deviation x L(jw), the
                # closed-loop phase transfer) -> broadcast audio above ~2 kHz
                # was destroyed. A 25 kHz loop tracks the full 0-15 kHz audio
                # band (|L(15 kHz)| = 0.94) and still rejects the out-of-band
                # FM noise (the 15 kHz audio LP then removes the pilot/subcarrier).
                # Rate-scaled for lower rates: at 2.016 MSPS (fs_ch = 144 kHz)
                # 25 kHz gives w_n = 1.09 rad/smp, which is UNSTABLE for the
                # 2nd-order arctan2 loop (limit cycle: 68 Hz peak larger than
                # the audio tones — verified in tests/test_wbfm_lowfs.py).
                # 0.14 x fs_ch keeps w_n <= 0.87 rad/smp (stable, |L(15k)| = 0.97).
                pll_f_n = min(25000.0, 0.14 * self.fs / 14)
                d = ChannelDemod(f, self.fs, fs_audio=self.band.get("audio_fs", 48000),
                                 wideband=True, bp_pass_hz=self.band.get("bp_pass_hz", 62000.0),
                                 bp_stop_hz=self.band.get("bp_stop_hz", 124000.0),
                                 bp_taps=self.band.get("bp_taps", 128),
                                 tau=self.band.get("tau", 50e-6),
                                 pll_f_n=pll_f_n, rf_center=f)
            else:
                # NBFM: the LO stays at the band centre (no retune), so the
                # carrier sits at (f - band_rf) in the IQ and the DDC must
                # downconvert by that offset, not the absolute f (which aliases
                # ~152 kHz off at 6.048 MSPS -> carrier out of the 7.5 kHz BPF).
                d = ChannelDemod(f, self.fs, rf_center=self.band["rf"])
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
        # the decimation chain only reaches channels within the /7 stage's
        # anti-alias FIR: h1 cuts at 0.40·fs/7 (345.6 kHz @ 6.048 MSPS,
        # 115.2 kHz @ 2.016 MSPS), so a channel further out is filtered
        # before demod (its signal is pure noise floor — verified
        # live: a +18 dB signal at 144.000 with LO 145.000 demodulates to
        # nothing). The waterfall still shows it; we just don't demod it.
        # 0.85 margin: at 2.016 MSPS this is 98.9 kHz, i.e. only the tuned
        # channel on the 100 kHz grid (the LO is retuned onto it anyway).
        win = 0.85 * 0.40 * self.fs / 7
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
        prev_active = self.active_now
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
        if (not active) and prev_active and str(self.args.band) == "wbfm" and self.retune is not None:
            # retune back to the band centre while the demod set is empty
            self.retune(self.band["rf"])
        for f in active:
            ctcss = self.args.ctcss if f in (self.args.channels or []) else None
            was = f in self.demods
            self._ensure_demod(f, ctcss)
            if (not was) and self.retune is not None and str(self.args.band) == "wbfm":
                # Phase 3 retune-on-detect: the RSPdx max IF BW (8 MHz) cannot span
                # 88-108 MHz in one tuning -> retune the LO onto the acquired
                # channel before demodulating it.
                self.retune(f)
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
                # keep ~2x the run length of final-rate audio (was hard-coded
                # 8 kHz: at 48 kHz wideband it discarded everything but the
                # last ~1/6 of the run)
                af = self.band.get("audio_fs", 8000)
                if sum(a.size for a in buf) > self.args.duration * af * 2:
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
            "spacing": self.band.get("spacing", 12.5e3),
            "dev_hz": self.band.get("dev"),
            "audio_fs_hz": self.band.get("audio_fs", 8000),
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
                "open_fraction": mt["open_samples"] / max(self.args.duration * self.band.get("audio_fs", 8000), 1),
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
            _write_wav(wav, audio.astype(np.float32),
                       int(self.band.get("audio_fs", 8000)))
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
    ap.add_argument("--antenna", default="a", choices=["a", "b"],
                    help="input port: a (default) or b (4.7 V bias-T available, "
                         "disabled unless --bias-t is given)")
    ap.add_argument("--bias-t", action="store_true",
                    help="enable the 4.7 V bias feed on port B (passive "
                         "antennas: leave OFF, the default)")
    ap.add_argument("--freq", type=float, default=None, help="LO frequency Hz (overrides band default)")
    ap.add_argument("--fs", type=float, default=None,
                    help="sample rate Hz (default: per-band — 2.016e6 wideband, "
                         "6.048e6 NBFM; RSPdx grid: multiples of 144 kHz)")
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
                    help="playback volume multiplier on the deviation-normalised gain "
                         "(0.6 / dev_max: 1/5000 for 3 kHz NBFM, ~1/83333 for 50 kHz WBFM)")
    ap.add_argument("--device", default=None,
                    help="output device index or name substring (default: system default)")
    ap.add_argument("--out", default="output")
    args = ap.parse_args()

    band = dict(BANDS[args.band])
    if args.freq:
        band["rf"] = args.freq
    # --fs default: per-band rate (wideband 2.016 MSPS, NBFM 6.048 MSPS);
    # the worker builds its chain from band["fs"], so keep them in lockstep.
    if args.fs is None:
        args.fs = band.get("fs", FS)
    band["fs"] = args.fs
    if args.band == "wbfm":
        # IF BW must fit the span: BW_1_536M (+/-0.768 MHz) at 2.016 MSPS;
        # BW_8M when the user overrides --fs higher (span-clipped to +/-fs/2).
        band["bw"] = BW_1_536M if args.fs <= 2.5e6 else BW_8M
    if args.ctcss is not None:
        near = min(CTCSS_TONES, key=lambda t: abs(t - args.ctcss))
        if abs(near - args.ctcss) > 0.01 * near:
            ap.error(f"CTCSS {args.ctcss} Hz not in ETSI TS 103 236 table (nearest {near})")
        args.ctcss = near
    h1_cut = 0.40 * args.fs / 7   # /7 stage anti-alias cutoff (demod reach)
    for f in (args.channels or []):
        if abs(f - band["rf"]) > h1_cut:
            print(f"[rx] WARNING: {f/1e6:.4f} MHz is >{h1_cut/1e3:.0f} kHz from the LO "
                  f"({band['rf']/1e6:.3f} MHz) — outside the demod chain's reach "
                  f"(h1 anti-alias cutoff {h1_cut/1e3:.1f} kHz @ {args.fs/1e6:.3f} MSPS); "
                  f"only the noise floor will be demodulated. Use --freq near the channel.")

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
    antenna_idx = 0 if args.antenna == "a" else 1
    device_info = {
        "radio": prof.key,
        "radio_label": prof.label,
        "serial": rx._ser(dev),
        "hwver": int(dev.hwVer),
        "sample_scale": prof.sample_scale,
        "antenna": args.antenna,
        "bias_t": bool(args.bias_t),
    }
    rx.get_params()
    rx.configure(fs_hz=args.fs, rf_hz=band["rf"], bw_khz=band["bw"],
                 lna=args.lna, gr_db=args.gr, agc=AGC_50HZ,
                 agc_setpoint_dbfs=args.agc_setpoint,
                 antenna=antenna_idx, bias_t=args.bias_t)
    worker = RxWorker(args, band, prof.sample_scale, device_info, retune=rx.retune)
    rx.start(worker.on_iq, worker.on_event)
    print(f"[rx] streaming: LO {band['rf']/1e6:.3f} MHz, fs {args.fs/1e6:.3f} MSPS, "
          f"IF BW {band['bw']/1e3:.1f} MHz, LNA {args.lna}, GR {args.gr} dB, "
          f"port {args.antenna.upper()} (bias-T {'ON' if args.bias_t else 'off'})")

    th = threading.Thread(target=worker.run, daemon=True)
    th.start()

    player = None
    if args.play:
        from playback import AudioPlayer
        # fs_in MUST be set at construction: __init__ derives the resample
        # ratio (_up) and the ramp length from it. Wideband final audio is
        # 48 kHz (fs_in == fs_out, no resample); NBFM is 8 kHz (6x upsample).
        # (Assigning player.fs_in after construction left _up=6 in place and
        # played the 48 kHz stream 6x too slow -> "8 kHz" band-limited sound.)
        player = AudioPlayer(fs_in=int(band.get("audio_fs", 8000)), vol=args.vol)
        # deviation-normalised playback gain: 0.6 / dev_max per band
        # (NBFM @ 3 kHz -> 0.6 FS; WBFM @ 50 kHz -> 0.6 FS)
        player.gain = (0.6 / band.get("dev", 3e3)) * args.vol
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
