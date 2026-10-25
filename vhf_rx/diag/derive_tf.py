"""Derive the exact discrete PLL transfer function (input freq -> NCO freq)
and compare with the continuous-time approximation. Kp, Ki for
f_n=1750 Hz, zeta=0.707, fs=16 kHz.

Result (2026-10-25): the exact discrete |T_f(e^{jw})| =
|[Kp(1-u)+Ki] / [(1-u)^2 + u(Kp(1-u)+Ki)]| (u = e^{-jw}) matches the measured
PLL output to 0.1 % (300/1000/3000 Hz -> 1.028/1.220/1.533); the
continuous-time (Kp s + Ki)/(s^2 + Kp s + Ki) matches at low f but is wrong
near/above the loop BW (0.845 vs measured 1.44 at 3 kHz). Implemented in
nbfm.pll_tracking_response.
"""
import numpy as np

F_N, ZETA, FS = 1750.0, 0.707, 16000.0
w_n = 2 * np.pi * F_N / FS
Kp, Ki = 2 * ZETA * w_n, w_n ** 2
print(f"Kp={Kp:.4f}  Ki={Ki:.4f}  w_n={w_n:.4f} rad/sample")


def Tf_discrete(f):
    """|T_f(e^{jw})| = |[Kp(1-u)+Ki] / [(1-u)^2 + u(Kp(1-u)+Ki)]|, u=e^{-jw}."""
    u = np.exp(-2j * np.pi * f / FS)
    num = Kp * (1 - u) + Ki
    den = (1 - u) ** 2 + u * (Kp * (1 - u) + Ki)
    return abs(num / den)


def Tf_continuous(f):
    """|(Kp*s+Ki)/(s^2+Kp*s+Ki)|, s=jw, continuous-time 2nd-order loop."""
    w = 2 * np.pi * f / FS
    num = np.sqrt(Ki ** 2 + Kp ** 2 * w ** 2)
    den = np.sqrt((Ki - w ** 2) ** 2 + Kp ** 2 * w ** 2)
    return num / den


print(f"{'f(Hz)':>8} {'|Tf_disc|':>10} {'|Tf_cont|':>10}")
for f in (300.0, 1000.0, 3000.0):
    print(f"{f:8.0f} {Tf_discrete(f):10.4f} {Tf_continuous(f):10.4f}")
print("\nmeasured PLL (diag/isolate_pll.py): 300:1.028  1000:1.220  3000:1.533")
