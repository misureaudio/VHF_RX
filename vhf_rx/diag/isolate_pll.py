"""Isolate the PLL: feed a pure tone directly into _pll_block at 16k (no
DDC, no stages) and measure the output amplitude vs input deviation.
This localizes demod accuracy to the PLL itself (2026-10-25: the PLL alone
matches theory exactly — 1.028/1.220/1.533 at 300/1000/3000 Hz — so any
full-chain discrepancy is in the DDC/decimation/stages, not the loop).
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from nbfm import _pll_block  # noqa: E402

F_N, ZETA, FS = 1750.0, 0.707, 16000.0
w_n = 2 * np.pi * F_N / FS
Kp, Ki = 2 * ZETA * w_n, w_n ** 2
scale = FS / (2 * np.pi)
DEV = 1800.0
N = 16000 * 2  # 2 s
SKIP = 320     # 20 ms lock transient


def measure(f0):
    n = N
    t = np.arange(n) / FS
    # input: carrier at 0 Hz (baseband), deviation DEV*sin(2 pi f0 t)
    phase = 2 * np.pi * np.cumsum(DEV * np.sin(2 * np.pi * f0 * t)) / FS
    x = np.exp(1j * phase).astype(np.complex64)
    out, _, _ = _pll_block(x, 0.0, 0.0, Kp, Ki, scale)
    # skip 20 ms lock transient
    out = out[SKIP:]
    w = np.hanning(out.size)
    Y = np.fft.rfft(out * w)
    amp = 2 * np.abs(Y) / np.sum(w)
    k = int(round(f0 * out.size / FS))
    return float(amp[k - 1:k + 2].max())


print("PLL-only (no stages), 1.8 kHz dev tone:")
for f0 in (300.0, 1000.0, 3000.0):
    a = measure(f0)
    print(f"  {f0:5.0f} Hz: out {a:7.1f} Hz  / dev {DEV} = {a / DEV:.3f}")

# Also: what does the NCO frequency phi look like in time at 300 Hz?
# Check the DC level and the 300 Hz amplitude separately.
t = np.arange(N) / FS
phase = 2 * np.pi * np.cumsum(DEV * np.sin(2 * np.pi * 300 * t)) / FS
x = np.exp(1j * phase).astype(np.complex64)
out, _, _ = _pll_block(x, 0.0, 0.0, Kp, Ki, scale)
out = out[SKIP:]
print(f"\n300 Hz case: output DC (mean) = {np.mean(out):.3f} Hz, "
      f"rms = {np.sqrt(np.mean(out**2)):.3f} Hz")
