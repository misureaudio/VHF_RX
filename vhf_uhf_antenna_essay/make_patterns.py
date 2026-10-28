"""Polar pattern plots for the corrected dual-band open-sleeve dipole
(symmetric 70 cm sleeve pair at x = +/-25 mm) at 145 MHz and 435 MHz.

PLANE DEFINITIONS (corrected):
  E-plane  = ELEVATION cut at phi=0 (x-z plane): theta swept 0..360.
             Figure-8, lobes broadside (horizontal), nulls end-fire (vertical).
  H-plane  = AZIMUTHAL cut at theta=90 (broadside, x-y plane): phi swept 0..360.
             For a z-directed dipole this is (nearly) OMNIDIRECTIONAL — a
             circle, never end-fire nulls.

Angle convention for both: polar angle 0 at the top, increasing clockwise.
E-plane: 0=+z (up), 90=+x (right).  H-plane: 0=+x (sleeve+), 180=-x (sleeve-).
Small red dots mark the sleeve directions on each plot.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import antenna_mom as am

L2, A2, N2 = 936.0e-3, 3.0e-3, 201
L7, x7 = 316.0e-3, 25.0e-3
N7 = 61
wires = [am.Wire(0.0, 0.0, -L2/2, L2/2, A2, N2, feed=True),
         am.Wire(x7, 0.0, -L7/2, L7/2, 3.0e-3, N7, feed=True),
         am.Wire(-x7, 0.0, -L7/2, L7/2, 3.0e-3, N7, feed=True)]
feed = N2 // 2
f0, f1 = 145.0e6, 435.0e6

def full_elevation(U, TH, PH):
    """Full 360 deg ELEVATION pattern (theta as the polar angle)."""
    th = TH[:, 0]
    i_r = int(np.argmin(np.abs(PH[0, :] - 0.0)))
    i_l = int(np.argmin(np.abs(PH[0, :] - np.pi)))
    Ur, Ul = U[:, i_r], U[:, i_l]
    n = 720
    alpha = np.linspace(0.0, 2.0*np.pi, n, endpoint=False)
    Uf = np.empty(n)
    for j, a in enumerate(alpha):
        Uf[j] = np.interp(a, th, Ur) if a <= np.pi else np.interp(2.0*np.pi - a, th, Ul)
    return alpha, Uf

def azimuthal(U, TH, PH):
    """AZIMUTHAL pattern at theta=90 (phi as the polar angle), 0..360."""
    th = TH[:, 0]
    ph = PH[0, :]
    ith = int(np.argmin(np.abs(th - np.pi / 2)))
    return ph, U[ith, :]

results = {}
for f in (f0, f1):
    Zin, I = am.solve(wires, f, feed, V0=1.0, nq_off=24, nq_self=400)
    _, U, TH, PH = am.far_field(wires, I, f, ntheta=721, nphi=720)
    aE, UE = full_elevation(U, TH, PH)
    aH, UH = azimuthal(U, TH, PH)
    results[f] = dict(Zin=Zin, aE=aE, UE=UE, aH=aH, UH=UH)
    print(f"{f/1e6:.0f} MHz: Zin={Zin.real:.1f}{Zin.imag:+.1f}j  SWR={am.swr(Zin):.2f}")

dB = lambda U: 10.0*np.log10(np.clip(U / U.max(), 1e-6, None))

fig = plt.figure(figsize=(12.0, 10.5))
titles = {f0: "145 MHz  (2 m)", f1: "435 MHz  (70 cm)"}
for row, f in enumerate((f0, f1)):
    r = results[f]
    for col, (a, U, pl, ylim, slv) in enumerate(
            ((r["aE"], r["UE"], "E-plane  (elevation, phi=0)", (-25, 1),
              (np.pi/2, 3*np.pi/2)),
             (r["aH"], r["UH"], "H-plane  (azimuth, theta=90)", (-6, 1),
              (0.0, np.pi)))):
        ax = fig.add_subplot(2, 2, row*2 + col + 1, projection="polar")
        d = dB(U)
        ax.plot(a, d, color="tab:blue", lw=2.0)
        ax.fill(a, d, color="tab:blue", alpha=0.22)
        ax.set_ylim(*ylim)
        ax.set_yticks([ylim[0], 0])
        ax.set_yticklabels([str(ylim[0]), "0 dB"], fontsize=8, color="0.4")
        ax.set_rlabel_position(60)
        ax.set_theta_offset(np.pi/2)
        ax.set_theta_direction(-1)
        # Sleeve directions (red dots just inside the outer ring).
        rr = ylim[1] - 0.35
        for ang in slv:
            ax.plot([ang], [rr], marker="o", color="tab:red", ms=6, zorder=5)
        ax.set_title(titles[f] + "   " + pl, fontsize=11, pad=22)
fig.suptitle("Dual-band VHF/UHF open-sleeve dipole — radiation patterns "
             "(normalised, 0 dB at max)\n"
             "driven 936 mm (centre-fed)  +  symmetric 70 cm sleeves at x = +/-25 mm   "
             "(red dots = sleeve directions)",
             fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.955])
out = "patterns_dualband.png"
fig.savefig(out, dpi=150)
print("wrote", out)
