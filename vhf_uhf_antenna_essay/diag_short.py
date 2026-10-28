"""Short-dipole Rr diagnostic: compare the solver's Re(Zin), the power-balance
Rr, and the analytical 20*pi^2*(L/lambda)^2 for a centre-fed short dipole.
Also print the current distribution to check it is ~triangular."""
import numpy as np
import antenna_mom as am

f = 144.0e6
lam = am.c0 / f
a = 3.0e-3

for Lf in (0.02, 0.05, 0.10):
    L = Lf * lam
    N = 51
    w = am.Wire(0, 0, -L/2, L/2, a, N, feed=True)
    Zin, I = am.solve([w], f, N//2)
    Iu = I / I[N//2]                 # unit feed current
    Prad = am.radiated_power([w], Iu, f)
    Rr_Zin = Zin.real
    Rr_Prad = 2 * Prad / abs(1.0)**2
    Rr_an = 20*np.pi**2*Lf**2
    # current distribution (normalised to peak)
    Ic = np.abs(Iu)
    print(f"\nL = {Lf*100:.0f}% lam  (L = {L*1e3:.1f} mm)")
    print(f"  Re(Zin)      = {Rr_Zin:.5f} ohm")
    print(f"  Rr (Prad)    = {Rr_Prad:.5f} ohm")
    print(f"  Rr (analytic)= {Rr_an:.5f} ohm   (20 pi^2 (L/lam)^2)")
    print(f"  ratio Re(Zin)/analytic = {Rr_Zin/Rr_an:.3f}")
    print(f"  Im(Zin)      = {Zin.imag:.1f} ohm")
    # show the current profile (every 5th segment)
    prof = ", ".join(f"{Ic[i]:.2f}" for i in range(0, N, 6))
    print(f"  |I| profile (feed centre): {prof}")
