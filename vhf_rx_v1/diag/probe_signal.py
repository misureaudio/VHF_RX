"""Live probe: what is on a given channel?

Streams the band, demodulates ONE channel, accumulates the 16 kHz complex
channel signal + 8 kHz audio, then classifies the signal:

  - CW/RTTY beacon: one narrow spectral line (>> neighbours)
  - wideband noise-like: flat channel spectrum (oscillator, digital, or
    out-of-band leak)
  - voice FM: power concentrated in the audio band, structured lines

Usage:
  python diag/probe_signal.py --band 2m --ch 144.0e6 --duration 15
  python diag/probe_signal.py --band 70cm --duration 15   # strongest channel
"""
from __future__ import annotations

import argparse
import os
import queue
import sys
import threading
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from nbfm import BandSpectrum, ChannelDemod           # noqa: E402
from prototype import BANDS, CHUNK, FS, band_lna_max  # noqa: E402
from sdrplay_api import AGC_50HZ, SDRPlay             # noqa: E402


class ProbeWorker:
    def __init__(self, band: dict, duration: float, channels: list[float],
                 sample_scale: float):
        self.band = band
        self.duration = duration
        self.channels = channels
        self.q: "queue.Queue" = queue.Queue(maxsize=0)
        self.fs = FS
        self.sample_scale = sample_scale
        self.spectrum = BandSpectrum(self.fs, band["rf"], band["fmin"], band["fmax"])
        self.demods: dict[float, ChannelDemod] = {}
        self.ch_buf: dict[float, list] = {f: [] for f in channels}
        self.au_buf: dict[float, list] = {f: [] for f in channels}
        self.n_blocks = 0
        self.drops = 0
        self.overloads = 0
        self.gain_db = float("nan")
        self._stop = threading.Event()
        self._t0 = 0.0

    def on_iq(self, xi, xq, n, reset):
        try:
            i = np.ctypeslib.as_array(xi, shape=(n,)).copy()
            q = np.ctypeslib.as_array(xq, shape=(n,)).copy()
            self.q.put_nowait((i, q))
        except queue.Full:
            self.drops += 1

    def on_event(self, event_id, tuner, gain_db):
        if event_id == 0 and gain_db == gain_db:
            self.gain_db = gain_db
        if event_id == 1:
            self.overloads += 1

    def _ensure(self, f: float):
        if f not in self.demods:
            self.demods[f] = ChannelDemod(f, self.fs)
            self.ch_buf[f] = []
            self.au_buf[f] = []

    def run(self):
        self._t0 = time.monotonic()
        pending, pending_n = [], 0
        active = list(self.channels)
        while not self._stop.is_set():
            try:
                i16, q16 = self.q.get(timeout=0.5)
            except queue.Empty:
                continue
            self.n_blocks += 1
            iq = ((i16.astype(np.float32) + 1j * q16.astype(np.float32))
                  / self.sample_scale).astype(np.complex64)
            pending.append(iq)
            pending_n += iq.size
            while pending_n >= CHUNK:
                allc = np.concatenate(pending)
                chunk, rest = allc[:CHUNK], allc[CHUNK:]
                pending = [rest] if rest.size else []
                pending_n = rest.size
                cp, chans, noise_db = self.spectrum.update(chunk)
                if not self.channels:                       # auto: strongest
                    if cp is not None:
                        thr = noise_db + 12.0
                        db = 10.0 * np.log10(np.maximum(cp, 1e-30))
                        idx = np.where(db >= thr)[0]
                        if len(idx):
                            top = float(chans[int(idx[np.argmax(db[idx])])])
                            if top not in active:
                                active.append(top)
                                self._ensure(top)
                for f in active:
                    self._ensure(f)
                    audio, _ = self.demods[f].process(chunk)
                    self.ch_buf[f].append(self.demods[f].x16_last.copy())
                    self.au_buf[f].append(audio)
            if time.monotonic() - self._t0 >= self.duration:
                break
        self._stop.set()

    def stop(self):
        self._stop.set()


def classify(f: float, x16: np.ndarray, audio: np.ndarray) -> None:
    fs16 = 16000.0
    x16 = x16[: x16.size // 2 * 2]      # even length for rfft
    n = x16.size
    w = np.hanning(n)
    # rfft is real-input-only; for the complex channel signal use the
    # positive-frequency half of fft (mathematically identical to rfft)
    X = np.fft.fft(x16 * w)[: n // 2 + 1]
    f16 = np.fft.rfftfreq(n, 1 / fs16)
    p = np.abs(X)
    k0 = int(np.argmax(p))
    pk_db = 20 * np.log10(p[k0] / max(p.sum() / len(p), 1e-30))
    # occupied bandwidth: 4*sigma of the power-weighted frequency (covers
    # 99% for a Gaussian-like spectrum; ~0 for a CW line)
    m = np.abs(f16) <= 6250.0
    fm, pm = f16[m], p[m] ** 2
    tot = float(pm.sum())
    if tot > 0:
        mu = float(np.sum(fm * pm) / tot)
        var = float(np.sum((fm - mu) ** 2 * pm) / tot)
        obw = 4.0 * float(np.sqrt(var))
    else:
        obw = 0.0
    # spectral flatness of the channel signal (0.3-3.4 kHz vs rest of band)
    mv = (np.abs(f16) >= 300) & (np.abs(f16) <= 3400)
    pv = p[mv] ** 2
    flat = float(np.exp(np.mean(np.log(pv + 1e-30))) / max(pv.mean(), 1e-30))
    # audio-band vs upper-band power on the 16k complex (carrier-detect metric)
    audio_sum = float(np.sum(p[(np.abs(f16) >= 31) & (np.abs(f16) <= 3031)] ** 2))
    ref_sum = float(np.sum(p[(np.abs(f16) >= 4000) & (np.abs(f16) <= 7970)] ** 2))
    carr_db = 10.0 * np.log10(max(audio_sum, 1e-30) / max(ref_sum, 1e-30))
    # top 5 lines
    print(f"=== {f/1e6:.4f} MHz: {n/fs16:.1f} s of 16 kHz channel signal")
    print(f"  peak line: {f16[k0]:+.1f} Hz, {pk_db:.1f} dB above mean bin")
    print(f"  99% occupied bandwidth: {obw/1000:.2f} kHz   "
          f"flatness(0.3-3.4k): {flat:.3f}   carrier-detect: {carr_db:+.1f} dB")
    top = np.argsort(p[m])[::-1][:5]
    for i in top:
        print(f"    line {f16[m][i]:+9.1f} Hz   {20*np.log10(p[m][i]/p[m].max()):+5.1f} dB")
    if audio.size:
        wa = np.hanning(audio.size)
        Ya = np.fft.rfft(audio * wa)
        fa = np.fft.rfftfreq(audio.size, 1 / 8000)
        amp = 2 * np.abs(Ya) / np.sum(wa)
        k1 = int(np.argmax(amp))
        topa = np.argsort(amp)[::-1][:4]
        print(f"  audio: rms {np.sqrt(np.mean(audio**2)):.4f}, top lines "
              + ", ".join(f"{fa[i]:.0f} Hz" for i in topa)
              + f" (strongest {amp[k1]:.4f} at {fa[k1]:.0f} Hz)")
    # classification
    if pk_db > 25 and obw < 500.0:
        print(f"  -> NARROW LINE at {f16[k0]:+.0f} Hz: CW/RTTY beacon or "
              f"single-tone artifact (no voice demod expected)")
    elif flat > 0.5:
        print("  -> WIDEBAND FLAT: noise-like / digital / out-of-band leak "
              "(not voice FM)")
    elif carr_db > 10:
        print("  -> audio-band concentrated: FM voice or voice-like signal")
    else:
        print("  -> mixed: see lines above")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--band", choices=list(BANDS), default="2m")
    ap.add_argument("--freq", type=float, default=None,
                    help="LO frequency in Hz (overrides the band default)")
    ap.add_argument("--ch", dest="channels", type=float, nargs="*", default=None,
                    help="channel(s) in Hz; default: auto-select strongest")
    ap.add_argument("--gr", type=int, default=50)
    ap.add_argument("--agc-setpoint", type=int, default=-20)
    ap.add_argument("--duration", type=float, default=15.0)
    args = ap.parse_args()

    band = dict(BANDS[args.band])
    if args.freq:
        band["rf"] = args.freq
    rx = SDRPlay()
    dev = rx.find_device()
    prof = rx.profile
    print(f"[probe] API v{rx.api_version:.2f}, device {rx._ser(dev)} "
          f"(hwVer {dev.hwVer}) -> {prof.label}")
    rx.get_params()
    rx.configure(fs_hz=FS, rf_hz=band["rf"], bw_khz=band["bw"],
                 lna=band_lna_max(prof, args.band),
                 gr_db=args.gr, agc=AGC_50HZ,
                 agc_setpoint_dbfs=args.agc_setpoint)
    w = ProbeWorker(band, args.duration, args.channels or [],
                    prof.sample_scale)
    rx.start(w.on_iq, w.on_event)
    print(f"[probe] LO {band['rf']/1e6:.3f} MHz, {args.duration:.0f} s, "
          f"channels: {args.channels or 'auto'}")
    th = threading.Thread(target=w.run, daemon=True)
    th.start()
    th.join(timeout=args.duration + 10)
    w.stop()
    time.sleep(0.3)
    rx.stop()
    print(f"[probe] blocks={w.n_blocks} drops={w.drops} overloads={w.overloads} "
          f"gain={w.gain_db:.1f} dB noise_floor={w.spectrum.noise_db:.1f} dB")
    cp = getattr(w.spectrum, "ch_power_last", None)
    if cp is not None:
        db = 10.0 * np.log10(np.maximum(cp, 1e-30))
        floor = w.spectrum.noise_db
        order = np.argsort(db)[::-1][:6]
        print("  strongest channels (full-band 12.5 kHz grid, absolute freq):")
        for i in order:
            print(f"    {w.spectrum.channels[i]/1e6:9.4f} MHz  "
                  f"{db[i]:6.1f} dB  ({db[i]-floor:+5.1f} dB over floor)")
    for f in w.demods:
        x16 = np.concatenate(w.ch_buf[f]) if w.ch_buf[f] else np.zeros(0)
        audio = np.concatenate(w.au_buf[f]) if w.au_buf[f] else np.zeros(0)
        if x16.size > 16000:
            classify(f, x16, audio)


if __name__ == "__main__":
    main()
