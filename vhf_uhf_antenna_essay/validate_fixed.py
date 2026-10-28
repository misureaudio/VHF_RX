"""Validation battery for the corrected thin-wire MoM solver (antenna_mom.py).

Checks, in order:
  1. Half-wave dipole anchors (0.47lam/3mm -> 70-j10, 0.50lam/1mm -> 81+j27).
  2. Short-dipole radiation resistance Rr = 20 pi^2 (L/lambda)^2.
  3. Convergence in N (mesh refinement).
  4. Resonant length scan (X=0).
  5. Classic dipole impedance table (Balanis Table 5.1).
  6. Power balance (Pin vs Prad for a unit feed current).
  7. Gain / pattern / 3-dB beamwidth.
  8. Two-wire (driven + parasitic) smoke test.
"""
import time
import numpy as np
import antenna_mom as am

f0 = 144.0e6
lam = am.c0 / f0
a = 3.0e-3

def dipole(L, a=a, N=201, f=f0, nq_off=24, nq_self=400):
    w = am.Wire(0, 0, -L/2, L/2, a, N, feed=True)
    return am.solve([w], f, N//2, nq_off=nq_off, nq_self=nq_self)

print("=" * 72)
print("VALIDATION  (corrected thin-wire MoM, Pocklington EFIE)")
print("=" * 72)

# --- 1. half-wave anchors ------------------------------------------------
print("\n[1] Half-wave dipole anchors @ 144 MHz")
for tag, Lf, ar in [("0.47 lam, a=3mm", 0.47, 3e-3),
                    ("0.50 lam, a=1mm", 0.50, 1e-3)]:
    t = time.time()
    Zin, I = dipole(Lf*lam, a=ar)
    print(f"    {tag:16s}: Zin = {Zin.real:7.2f} {Zin.imag:+7.2f} j ohm   "
          f"SWR(50)={am.swr(Zin):5.2f}   ({time.time()-t:.1f}s)")
print("    [classic: 0.47lam/3mm ~ 70 - j10 ; 0.50lam/1mm ~ 81 + j27]")

# --- 2. short dipole Rr --------------------------------------------------
print("\n[2] Short-dipole radiation resistance  Rr = 20 pi^2 (L/lam)^2")
for Lf in (0.02, 0.05):
    Zin, I = dipole(Lf*lam, N=51)
    Rr_ex = 20*np.pi**2*Lf**2
    print(f"    L={Lf*100:4.1f}% lam: Zin = {Zin.real:8.4f} {Zin.imag:+10.2f} j   "
          f"(Rr exact = {Rr_ex:8.4f})")

# --- 3. convergence in N -------------------------------------------------
print("\n[3] Convergence in N  (0.47 lam dipole @144 MHz)")
Zref, _ = dipole(0.47*lam, N=401)
print(f"    N= 401 (ref): Z = {Zref.real:8.3f} {Zref.imag:+9.3f} j")
for N in (101, 201, 301):
    Z, _ = dipole(0.47*lam, N=N)
    print(f"    N= {N:4d}:    Z = {Z.real:8.3f} {Z.imag:+9.3f} j   "
          f"|dZ| = {abs(Z-Zref):.3e}")

# --- 4. resonant length scan --------------------------------------------
print("\n[4] Resonant length scan (X=0)")
prev = None
found = None
for Lm in np.arange(0.40, 0.56, 0.002):
    Z, _ = dipole(Lm*lam, N=151)
    if prev is not None and np.sign(prev.imag) != np.sign(Z.imag):
        found = Lm
        break
    prev = Z
print(f"    resonant L/lambda ~ {found}   [classic ~0.473]")

# --- 5. classic table ----------------------------------------------------
print("\n[5] Classic dipole impedance table (Balanis Table 5.1, a->0)")
for Lf, (Rc, Xc) in {0.25: (18.8, -108.8), 0.5: (73.1, 42.5),
                     0.75: (243.8, 104.5), 1.0: (152.3, -46.1)}.items():
    Z, _ = dipole(Lf*lam, N=201)
    print(f"    L/lambda={Lf:4.2f}: Z = {Z.real:8.1f} {Z.imag:+8.1f} j   "
          f"[classic {Rc:6.1f} {Xc:+7.1f}]")

# --- 6. power balance ----------------------------------------------------
print("\n[6] Power balance  (unit feed current, 0.47 lam dipole)")
w = am.Wire(0, 0, -0.47*lam/2, 0.47*lam/2, a, 201, feed=True)
Zin, I = am.solve([w], f0, 201//2)
# normalise to unit feed current
Iu = I / I[201//2]
Pin = 0.5 * Zin.real * abs(1.0)**2
Prad = am.radiated_power([w], Iu, f0)
print(f"    Pin = {Pin:.5f} W,  Prad = {Prad:.5f} W,  ratio = {Prad/Pin:.5f}")
print("    [ratio should be ~1; the small deficit is conductor loss]")

# --- 7. gain / pattern ---------------------------------------------------
print("\n[7] Gain / pattern  (0.47 lam dipole)")
m = am.pattern_metrics([w], Iu, f0)
print(f"    D = {m['D']:.3f} ({m['D_dBi']:.2f} dBi)   "
      f"G = {m['G']:.3f} ({m['G_dBi']:.2f} dBi)   eff = {m['eff']*100:.2f}%")
print(f"    E-plane 3-dB bw = {m['bw_E_deg']:.1f} deg   "
      f"H-plane 3-dB bw = {m['bw_H_deg']:.1f} deg   "
      f"F/B = {m['F/B_dB']:.1f} dB")
print("    [classic: D ~ 1.5 (2.15 dBi), E-plane bw ~ 78 deg]")

# --- 8. two-wire smoke test ---------------------------------------------
print("\n[8] Two-wire (driven + parasitic) smoke test")
w1 = am.Wire(0, 0, -0.47*lam/2, 0.47*lam/2, a, 201, feed=True)
w2 = am.Wire(0.30, 0, -0.15*lam/2, 0.15*lam/2, a, 61)   # parasitic
Zin2, I2 = am.solve([w1, w2], f0, 201//2)
print(f"    Zin (with parasitic) = {Zin2.real:7.2f} {Zin2.imag:+7.2f} j ohm")
print(f"    parasitic peak |I|/|I_feed| = "
      f"{np.max(np.abs(I2[201:]))/np.abs(I2[201//2]):.3f}")
Zp, _ = am.z_parameters([w1, w2], f0, [201//2, 201+61//2])
print(f"    open-circuit Z11 = {Zp[0,0].real:7.2f} {Zp[0,0].imag:+7.2f} j   "
      f"Z12 = {Zp[0,1].real:7.2f} {Zp[0,1].imag:+7.2f} j")
print(f"    loaded Zin (parasitic open, ZL=inf) = "
      f"{am.loaded_zin(Zp, 1e9).real:7.2f} "
      f"{am.loaded_zin(Zp, 1e9).imag:+7.2f} j")

print("\nAll validation checks complete.")
