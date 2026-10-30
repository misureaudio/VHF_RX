"""Offline: wideband chain at 2.016 MSPS (RSPdx 144 kHz grid).

The wideband band defaults to 2.016 MSPS (BANDS["wbfm"]["fs"] in
prototype.py): at 6.048 MSPS the chain needs ~0.95 of one core on the dev
machine, the worker sustains only ~3.6 MSPS (59.7 % of real time) and the
audio ring underruns 40 % of a 300 s run. At 2.016 MSPS the chain costs ~1/3.

Checks:
  1. Decimation grid: fs/14 = 144 kHz -> final /3 lands on 48 kHz.
  2. The production rate-scaled PLL loop (0.14 x fs_ch = 20.16 kHz,
     w_n = 0.87 rad/sample) is stable and recovers the tones.
  3. Reference: the 25 kHz loop (w_n = 1.09 rad/sample) is UNSTABLE at
     144 kHz (limit cycle) — printed for the record, not asserted.

Run: D:/Source/hermes-dir/.venv/Scripts/python.exe vhf_rx_v3/tests/test_wbfm_lowfs.py
"""
import sys

import numpy as np

sys.path.insert(0, r"D:/Source/hermes-dir/vhf_rx_v3")

from nbfm import ChannelDemod  # noqa: E402

FS = 2.016e6
F0 = 103.0e6
BLOCK = 16384
FAIL = 0


def fm_modulate(audio_hz, f_carrier, fs, amp=0.5):
    """audio_hz: instantaneous frequency deviation in Hz. Complex baseband."""
    phi = 2.0 * np.pi * np.cumsum(audio_hz) / fs
    t = np.arange(audio_hz.size) / fs
    return amp * np.exp(1j * (phi + 2.0 * np.pi * f_carrier * t))


def check(name, ok, detail):
    global FAIL
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")
    if not ok:
        FAIL += 1


def tone_freqs(x, fs):
    spec = np.abs(np.fft.rfft(x * np.hanning(x.size)))
    freqs = np.fft.rfftfreq(x.size, 1.0 / fs)
    out = []
    i = 1
    while i < spec.size - 1:
        if spec[i] > spec[i - 1] and spec[i] >= spec[i + 1] and spec[i] > 0.02 * spec.max():
            out.append(float(freqs[i]))
        i += 1
    return out


def run(fs, pll_f_n, seconds=1.5, label=""):
    n = int(seconds * fs)
    t = np.arange(n) / fs
    audio = (5000.0 * np.sin(2.0 * np.pi * 1000.0 * t)
             + 3000.0 * np.sin(2.0 * np.pi * 10000.0 * t))
    iq = fm_modulate(audio, 0.0, fs)  # carrier at DC (LO on the station)
    d = ChannelDemod(F0, fs, fs_audio=48000.0, wideband=True,
                     bp_pass_hz=62000.0, bp_stop_hz=124000.0, bp_taps=128,
                     tau=50e-6, pll_f_n=pll_f_n, rf_center=F0)
    outs = [d.process(iq[i:i + BLOCK])[0] for i in range(0, iq.size, BLOCK)]
    out = np.concatenate(outs[int(0.6 * fs / BLOCK):])
    lines = tone_freqs(out, 48000.0)
    f1 = min(lines, key=lambda f: abs(f - 1000.0))
    f10 = min(lines, key=lambda f: abs(f - 10000.0))
    spec = np.abs(np.fft.rfft(out * np.hanning(out.size)))
    k1, k10 = int(f1 * out.size / 48000.0), int(f10 * out.size / 48000.0)
    a1, a10 = spec[k1], spec[k10]
    # Hann window coherent gain 0.5 x DFT N/2 -> expected peak A*N/4:
    e1 = 5000.0 * 0.995 * out.size / 4.0   # de-emphasis(1k) ~ 0.995, loop gain ~ 1
    e10 = 3000.0 * 0.306 * out.size / 4.0  # de-emphasis(10k) ~ 0.306
    print(f"  {label} tones: 1k -> {f1:.1f} Hz (amp {a1/e1*100:.0f}% of expected "
          f"{e1:.0f}), 10k -> {f10:.1f} Hz (amp {a10/e10*100:.0f}% of expected {e10:.0f})")
    # top-5 peaks for inspection:
    top = np.argsort(spec)[::-1][:5]
    top_hz = np.fft.rfftfreq(out.size, 1.0 / 48000.0)[top]
    top_db = 20.0 * np.log10(spec[top] / spec.max() + 1e-30)
    print(f"  {label} top peaks: " + ", ".join(f"{h:.0f} Hz @ {db:.1f} dB" for h, db in zip(top_hz, top_db)))
    # slip/limit-cycle detection: residual power outside DC and the tone bands
    resid = spec.copy()
    resid[:10] = 0.0
    resid[k1 - 10:k1 + 11] = 0.0
    resid[k10 - 10:k10 + 11] = 0.0
    resid_db = 10.0 * np.log10((resid ** 2).sum() / (spec.max() ** 2) + 1e-30)
    return f1, f10, a1 / e1, a10 / e10, resid_db


def main():
    print(f"== wideband chain @ {FS/1e6:.3f} MSPS ==")
    print(f"   fs/14 = {FS/14:.1f} Hz, final factor = {int(round(FS/14/48000))}")

    # production config: prototype.py uses pll_f_n = min(25000, 0.14*fs/14)
    prod_loop = min(25000.0, 0.14 * FS / 14)
    f1, f10, r1, r10, resid = run(FS, prod_loop,
                                  label=f"production pll_f_n={prod_loop:.0f} "
                                        f"(w_n={2*np.pi*prod_loop/(FS/14):.2f} rad/smp)")
    check("production tone freqs", abs(f1 - 1000.0) < 5.0 and abs(f10 - 10000.0) < 8.0,
          f"1k={f1:.1f} Hz, 10k={f10:.1f} Hz")
    check("production tone amplitudes", 0.85 < r1 < 1.15 and 0.80 < r10 < 1.25,
          f"1k {r1*100:.0f}%, 10k {r10*100:.0f}%")
    check("production no slip/limit-cycle (residual < -30 dB)", resid < -30.0,
          f"residual {resid:.1f} dB rel. peak")

    print()
    f1, f10, r1, r10, resid = run(FS, 25000.0,
                                  label="25k loop @ 144 kHz (UNSTABLE reference, "
                                        "w_n=1.09 rad/smp — not asserted)")
    print(f"  (25 kHz reference: 1k={f1:.1f} Hz {r1*100:.0f}%, 10k={f10:.1f} Hz "
          f"{r10*100:.0f}%, residual {resid:.1f} dB — limit cycle expected)")

    print()
    print("RESULT:", "FAIL" if FAIL else "PASS", f"({FAIL} failures)")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
