"""Verify the open-sleeve physics: does a CLOSE-coupled 70 cm sleeve
(spacing ~10-25 mm, reactive near-field) match the 1.36 lambda driven element
on the 70 cm band, versus the 307 mm (0.45 lambda) Yagi-like spacing?

Driven: 936 mm (0.452 lambda @145 MHz) at x=0, centre-fed.
Sleeve: 316 mm (0.459 lambda @435 MHz).
Test single offset sleeves at several spacings, and a symmetric pair.
"""
import numpy as np
import antenna_mom as am

L2 = 936.0e-3
A2 = 3.0e-3
N2 = 201
L7 = 316.0e-3
A7 = 3.0e-3
N7 = 61
f7 = 435.0e6
lam7 = 299792458.0 / f7
print(f"70cm lambda = {lam7*1e3:.1f} mm; driven = {L2/lam7:.3f} lambda")

def solve_spacing(x7_list, label):
    wires = [am.Wire(0.0, 0.0, -L2/2, L2/2, A2, N2, feed=True)]
    for x7 in x7_list:
        wires.append(am.Wire(x7, 0.0, -L7/2, L7/2, A7, N7, feed=True))
    feed = N2 // 2
    Zin, I = am.solve(wires, f7, feed, V0=1.0, nq_off=24, nq_self=400)
    print(f"{label:28s} Zin={Zin.real:7.1f} {Zin.imag:+8.1f}j  SWR={am.swr(Zin):5.2f}")
    return Zin

print("\n--- Single 70cm sleeve, varying spacing (435 MHz) ---")
for mm in (307, 100, 50, 25, 20, 15, 10):
    x = mm * 1e-3
    solve_spacing([x], f"sleeve @ +{mm} mm ({x/lam7:.3f} lam)")

print("\n--- Symmetric pair of 70cm sleeves (435 MHz) ---")
for mm in (25, 20, 15, 10):
    x = mm * 1e-3
    solve_spacing([-x, x], f"pair @ +/-{mm} mm ({x/lam7:.3f} lam)")
