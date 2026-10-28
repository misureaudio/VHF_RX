"""Re-verify the 2-port Z-parameters for the optimised dual-band design at
435 MHz with the corrected z_parameters sign, and cross-check that the
single-port open Z11 equals solve()'s Zin (self-consistency)."""
import numpy as np
import antenna_mom as am

L2, L7, x7 = 936.0e-3, 320.0e-3, 280.0e-3
A2 = A7 = 3.0e-3
N2, N7 = 201, 61
w2 = am.Wire(0.0, 0.0, -L2/2, L2/2, A2, N2, feed=True)
w7 = am.Wire(x7, 0.0, -L7/2, L7/2, A7, N7, feed=True)
wo = [w2, w7]
f7 = 435.0e6

# Single-port: drive driven element, parasitic present (open).
feed = N2 // 2
Zin, I = am.solve(wo, f7, feed, V0=1.0, nq_off=24, nq_self=300)
print(f"solve()  Zin (driven, 70cm) = {Zin.real:.2f} {Zin.imag:+.2f}j  SWR={am.swr(Zin):.2f}")

# 2-port open-circuit Z.
feed_segs = [N2//2, N2 + N7//2]
Zp, _ = am.z_parameters(wo, f7, feed_segs, nq_off=24, nq_self=300)
print(f"2-port Z11 = {Zp[0,0].real:.2f} {Zp[0,0].imag:+.2f}j")
print(f"2-port Z12 = {Zp[0,1].real:.2f} {Zp[0,1].imag:+.2f}j   |Z12|={abs(Zp[0,1]):.2f}")
print(f"2-port Z22 = {Zp[1,1].real:.2f} {Zp[1,1].imag:+.2f}j")
print(f"loaded Zin (parasitic open, ZL=inf) = {am.loaded_zin(Zp, 1e9).real:.2f} {am.loaded_zin(Zp,1e9).imag:+.2f}j")
print(f"Z11 - Z12^2/Z22 = {Zp[0,0] - Zp[0,1]*Zp[1,0]/Zp[1,1]:.2f}")
print()
print(f"self-consistency: |Z11 - Zin_solve| = {abs(Zp[0,0]-Zin):.3f} ohm")
print(f"coupling coeff |Z12|/sqrt(|Z11 Z22|) = {abs(Zp[0,1])/np.sqrt(abs(Zp[0,0])*abs(Zp[1,1])):.3f}")
