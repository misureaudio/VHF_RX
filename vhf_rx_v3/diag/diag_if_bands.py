"""Full-band energy distribution of received FM content at 103.0 MHz (diag).

Stages measured:
  raw IQ  (6.048 MSPS complex)  : amplitude + band energy
  baseband (432 kHz complex)    : amplitude bands (pilot check @19 kHz)
  IF     (instantaneous freq)   : band energy = the actual FM content
  audio  (48 kHz)               : band energy after PLL/LP/de-emph
"""
from __future__ import annotations

import argparse
import queue
import sys
import threading
import time

import numpy as np

sys.path.insert(0, "D:/Source/hermes-dir/vhf_rx_v3")

from nbfm import ChannelDemod                            # noqa: E402
from prototype import BANDS, FS, band_lna_max          # noqa: E402
from sdrplay_api import AGC_50HZ, SDRPlay             # noqa: E402


def bands(x, fs, edges, names):
    x = np.asarray(x)
    if np.iscomplexobj(x):
        f = np.fft.fftfreq(x.size, 1 / fs)
        X = np.abs(np.fft.fft(x * np.hanning(x.size))) ** 2
        m0 = (f >= 0) & (f < fs / 2)
        f, X = f[m0], X[m0]
    else:
        f = np.fft.rfftfreq(x.size, 1 / fs)
        X = np.abs(np.fft.rfft(x * np.hanning(x.size))) ** 2
    tot = X.sum()
    out = []
    for lo, hi, nm in zip(edges[:-1], edges[1:], names):
        m = (f >= lo) & (f < hi)
        out.append((nm, 100.0 * X[m].sum() / tot))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--freq", type=float, default=103.0e6)
    ap.add_argument("--duration", type=float, default=10.0)
    args = ap.parse_args()

    band = dict(BANDS["wbfm"])
    band["rf"] = args.freq
    rx = SDRPlay()
    dev = rx.find_device()
    prof = rx.profile
    rx.get_params()
    rx.configure(fs_hz=FS, rf_hz=band["rf"], bw_khz=band["bw"],
                 lna=band_lna_max(prof, "wbfm"), gr_db=50,
                 agc=AGC_50HZ, agc_setpoint_dbfs=-20, antenna=0, bias_t=False)
    q = queue.Queue(maxsize=0)

    def on_iq(xi, xq, n, reset):
        i = np.ctypeslib.as_array(xi, shape=(n,)).copy()
        qq = np.ctypeslib.as_array(xq, shape=(n,)).copy()
        q.put_nowait((i, qq, int(reset)))

    rx.start(on_iq, None)
    f = args.freq
    d = ChannelDemod(f, FS, fs_audio=48000, wideband=True,
                     bp_pass_hz=62000.0, bp_stop_hz=124000.0, bp_taps=128,
                     tau=50e-6, rf_center=f)
    t0 = time.monotonic()
    x16_all, au_all, iq_all, resets = [], [], [], 0
    # skip the first 1 s (LO/DC-corrector/AGC settling) for the records
    skip_until = t0 + 1.0
    while time.monotonic() - t0 < args.duration:
        try:
            i16, q16, rst = q.get(timeout=0.2)
        except Exception:
            continue
        resets += rst
        iq = ((i16.astype(np.float32) + 1j * q16.astype(np.float32))
              / prof.sample_scale).astype(np.complex64)
        if time.monotonic() >= skip_until:
            iq_all.append(iq)
            audio, _ = d.process(iq)
            x16_all.append(d.x16_last.copy())
            au_all.append(audio)
        else:
            d.process(iq)
    rx.stop()

    iq = np.concatenate(iq_all)
    x16 = np.concatenate(x16_all)
    au = np.concatenate(au_all)
    print(f"[dbg] captured {iq.size/FS:.2f} s (after 1 s skip), "
          f"reset events: {resets}")

    # raw IQ amplitude: where is the carrier + sideband extent?
    iq = iq[: iq.size // 2 * 2]
    w0 = np.hanning(iq.size)
    X0 = np.abs(np.fft.fft(iq * w0))
    f0 = np.fft.fftfreq(iq.size, 1 / FS)
    p0 = X0 ** 2
    tot0 = p0.sum()
    print("RAW IQ band energy (%):")
    for lo, hi in [(-3e6, 3e6), (-1e6, 1e6), (-300e3, 300e3),
                   (-100e3, 100e3), (-30e3, 30e3), (-10e3, 10e3),
                   (-1e3, 1e3), (-200, 200)]:
        m = (f0 >= lo) & (f0 < hi)
        print(f"  {lo/1e3:8.0f} .. {hi/1e3:8.0f} kHz : "
              f"{100*p0[m].sum()/tot0:6.2f}%")
    # carrier line vs sidebands: power in |f|<50 Hz vs 50 Hz..300 kHz
    mc = np.abs(f0) < 50.0
    ms = (np.abs(f0) >= 50.0) & (np.abs(f0) <= 300e3)
    print(f"  carrier line |f|<50 Hz: {100*p0[mc].sum()/tot0:6.2f}%   "
          f"sidebands 50 Hz-300 kHz: {100*p0[ms].sum()/tot0:6.2f}%")

    # baseband (432 kHz complex): pilot check + band energy
    x16 = x16[: x16.size // 2 * 2]
    print("BASEBAND (432 kHz complex) band energy (%):")
    for nm, v in bands(x16, 432000.0,
                       [0, 1e3, 10e3, 19e3, 20e3, 30e3, 50e3, 100e3, 216e3],
                       ["0-1k", "1k-10k", "19k-20k", "20k-30k",
                        "30k-50k", "50k-100k", "100k-216k"]):
        print(f"  {nm:10s} : {v:6.2f}%")
    # explicit pilot line at 19 kHz
    w1 = np.hanning(x16.size)
    X1 = np.fft.fft(x16 * w1)
    f1 = np.fft.fftfreq(x16.size, 1 / 432000.0)
    k19 = int(np.argmin(np.abs(f1 - 19000.0)))
    pk = 20 * np.log10(np.abs(X1[k19]) / (np.abs(X1).mean() + 1e-30))
    k0 = int(np.argmin(np.abs(f1)))
    rel = 20 * np.log10(np.abs(X1[k19]) / (np.abs(X1[k0]) + 1e-30))
    print(f"  19 kHz pilot line: {pk:6.1f} dB above mean bin "
          f"({rel:+6.1f} dB vs carrier)")

    # IF: the actual FM content
    ip = np.imag(x16[1:] * np.conj(x16[:-1]))
    ip = np.arctan2(ip, np.real(x16[1:] * np.conj(x16[:-1])))
    ip_hz = ip * (432000.0 / (2.0 * np.pi))
    ip_hz = ip_hz[: ip_hz.size // 2 * 2]
    print("IF (FM content, Hz deviation) band energy (%):")
    for nm, v in bands(ip_hz, 432000.0,
                       [0, 100, 300, 1e3, 3e3, 10e3, 19e3, 20e3,
                        50e3, 100e3, 216e3],
                       ["0-100", "100-300", "300-1k", "1k-3k", "3k-10k",
                        "10k-19k", "19k-20k", "20k-50k", "50k-100k",
                        "100k-216k"]):
        print(f"  {nm:10s} : {v:6.2f}%")

    # audio (48 kHz)
    au = au[: au.size // 2 * 2]
    print("AUDIO (48 kHz, after PLL/LP/de-emph) band energy (%):")
    for nm, v in bands(au.astype(np.float64), 48000.0,
                       [0, 100, 300, 1e3, 3e3, 6e3, 10e3, 15e3, 24e3],
                       ["0-100", "100-300", "300-1k", "1k-3k", "3k-6k",
                        "6k-10k", "10k-15k", "15k-24k"]):
        print(f"  {nm:10s} : {v:6.2f}%")
    print(f"  audio rms {np.sqrt(np.mean(au**2)):.4f}, "
          f"peak {np.abs(au).max():.4f}")


if __name__ == "__main__":
    main()
