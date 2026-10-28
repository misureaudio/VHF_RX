"""Dual-band VHF-UHF open-sleeve dipole design for the 144-148 / 430-440 MHz
amateur bands.

Design: a centre-fed 2 m driven dipole (resonant near 145 MHz) plus a 70 cm
parasitic element (resonant near 435 MHz) offset to one side.  The parasitic
is open-circuited; its current is induced.  The feedpoint is on the driven
element.

Steps:
  1. Baseline DK7ZB geometry -> SWR sweeps on both bands + patterns.
  2. 2-D geometry optimisation (parasitic length x position) to minimise the
     70 cm SWR at 435 MHz while keeping the 2 m SWR good.
  3. Optimised design -> full SWR sweeps + patterns + 2-port Z-parameters.
  4. Write results to dualband_results.json and dualband_sweeps.csv.
"""
import json
import time
import numpy as np
import antenna_mom as am

c = am.c0
lam2 = c / 145.0e6     # 2 m band wavelength ~2.07 m
lam7 = c / 435.0e6     # 70 cm band wavelength ~0.69 m

A2 = 3.0e-3           # driven element radius (6 mm rod)
A7 = 3.0e-3           # parasitic radius
N2 = 201
N7 = 61

def make_wires(L2, L7, x7):
    w1 = am.Wire(0, 0, -L2/2, L2/2, A2, N2, feed=True)
    w2 = am.Wire(x7, 0, -L7/2, L7/2, A7, N7)
    return [w1, w2]

def swr_at(wires, f, feed):
    Zin, I = am.solve(wires, f, feed, nq_off=24, nq_self=300)
    return Zin, am.swr(Zin), I

print("=" * 74)
print("DUAL-BAND VHF-UHF OPEN-SLEEVE DIPOLE  (144-148 / 430-440 MHz)")
print("=" * 74)

# --------------------------------------------------------------------------
# 1. Baseline DK7ZB geometry
# --------------------------------------------------------------------------
L2b, L7b, x7b = 0.936, 0.316, 0.307
print(f"\n[1] Baseline DK7ZB:  2m driven L={L2b*1e3:.0f}mm  "
      f"70cm para L={L7b*1e3:.0f}mm @ x={x7b*1e3:.0f}mm")
wb = make_wires(L2b, L7b, x7b)
feed = N2 // 2

for band, f0, f1 in [("2m", 144.0, 148.0), ("70cm", 430.0, 440.0)]:
    print(f"\n    SWR sweep {band}:")
    print(f"    {'f(MHz)':>8s} {'Re(Z)':>8s} {'Im(Z)':>9s} {'SWR':>6s}")
    for ff in np.arange(f0, f1 + 1e-9, 1.0):
        Zin, s, _ = swr_at(wb, ff*1e6, feed)
        print(f"    {ff:8.1f} {Zin.real:8.2f} {Zin.imag:+9.2f} {s:6.2f}")

# --------------------------------------------------------------------------
# 2. Optimise parasitic (L7, x7) for 70cm SWR at 435 MHz
# --------------------------------------------------------------------------
print("\n[2] Optimising 70cm parasitic (L7, x7) for min SWR @ 435 MHz")
best = None
L7s = np.arange(0.28, 0.36 + 1e-9, 0.02)
x7s = np.arange(0.24, 0.40 + 1e-9, 0.04)
grid = []
for L7 in L7s:
    for x7 in x7s:
        w = make_wires(L2b, L7, x7)
        Zin, s, _ = swr_at(w, 435.0e6, feed)
        grid.append((s, L7, x7, Zin))
        print(f"    L7={L7*1e3:5.0f}mm x7={x7*1e3:5.0f}mm  "
              f"SWR(435)={s:6.2f}  Z={Zin.real:7.1f}{Zin.imag:+7.1f}j")
grid.sort(key=lambda t: t[0])
s7, L7o, x7o, Z7o = grid[0]
print(f"\n    Best 70cm: L7={L7o*1e3:.0f}mm x7={x7o*1e3:.0f}mm  "
      f"SWR={s7:.2f}")

# check the 2m SWR at the optimised geometry
w2m = make_wires(L2b, L7o, x7o)
Z2m, s2m, _ = swr_at(w2m, 145.0e6, feed)
print(f"    2m SWR at optimised geometry (145 MHz) = {s2m:.2f}  "
      f"(Z={Z2m.real:.1f}{Z2m.imag:+.1f}j)")

# --------------------------------------------------------------------------
# 3. Optimised design: full sweeps + patterns + 2-port Z
# --------------------------------------------------------------------------
print("\n[3] Optimised design:  2m driven L=%.0fmm  70cm para L=%.0fmm @ x=%.0fmm"
      % (L2b*1e3, L7o*1e3, x7o*1e3))
wo = make_wires(L2b, L7o, x7o)

results = {"baseline": {"L2_mm": L2b*1e3, "L7_mm": L7b*1e3, "x7_mm": x7b*1e3},
           "optimised": {"L2_mm": L2b*1e3, "L7_mm": L7o*1e3, "x7_mm": x7o*1e3}}

for name, wires in [("baseline", wb), ("optimised", wo)]:
    results[name]["swr"] = {}
    for band, f0, f1 in [("2m", 144.0, 148.0), ("70cm", 430.0, 440.0)]:
        row = []
        for ff in np.arange(f0, f1 + 1e-9, 0.5):
            Zin, s, _ = swr_at(wires, ff*1e6, feed)
            row.append({"f_MHz": ff, "Re_Z": round(Zin.real, 2),
                        "Im_Z": round(Zin.imag, 2), "SWR": round(s, 3)})
        results[name]["swr"][band] = row
        swrs = [r["SWR"] for r in row]
        print(f"    {name:10s} {band:5s}: SWR min={min(swrs):.2f} max={max(swrs):.2f} "
              f"mean={np.mean(swrs):.2f}")

    for tag, ff in [("145MHz", 145.0e6), ("435MHz", 435.0e6)]:
        Zin, _, I = swr_at(wires, ff, feed)
        m = am.pattern_metrics(wires, I, ff, ntheta=721, nphi=720)
        results[name][tag] = {
            "Zin": [round(Zin.real, 2), round(Zin.imag, 2)],
            "SWR": round(am.swr(Zin), 3),
            "D_dBi": round(m["D_dBi"], 2), "G_dBi": round(m["G_dBi"], 2),
            "eff_pct": round(m["eff"]*100, 2),
            "bw_E_deg": round(m["bw_E_deg"], 1), "bw_H_deg": round(m["bw_H_deg"], 1),
            "F/B_dB": round(m["F/B_dB"], 1),
        }
        print(f"    {name:10s} {tag}: Zin={Zin.real:.1f}{Zin.imag:+.1f}j "
              f"SWR={am.swr(Zin):.2f}  G={m['G_dBi']:.2f} dBi  "
              f"bwE={m['bw_E_deg']:.0f}deg  F/B={m['F/B_dB']:.1f}dB")

# 2-port open-circuit Z-parameters at 435 MHz (optimised)
Zp, _ = am.z_parameters(wo, 435.0e6, [feed, N2 + N7//2], nq_off=24, nq_self=300)
results["optimised"]["Z2port_435"] = {
    "Z11": [round(Zp[0,0].real,2), round(Zp[0,0].imag,2)],
    "Z12": [round(Zp[0,1].real,2), round(Zp[0,1].imag,2)],
    "Z22": [round(Zp[1,1].real,2), round(Zp[1,1].imag,2)],
    "loaded_Zin": [round(am.loaded_zin(Zp, 1e9).real,2),
                   round(am.loaded_zin(Zp, 1e9).imag,2)],
}
print(f"\n    2-port Z @435MHz (optimised): Z11={Zp[0,0].real:.1f}{Zp[0,0].imag:+.1f}j "
      f"Z12={Zp[0,1].real:.1f}{Zp[0,1].imag:+.1f}j Z22={Zp[1,1].real:.1f}{Zp[1,1].imag:+.1f}j")
print(f"    loaded Zin (parasitic open) = "
      f"{am.loaded_zin(Zp,1e9).real:.1f}{am.loaded_zin(Zp,1e9).imag:+.1f}j")

with open("dualband_results.json", "w") as fh:
    json.dump(results, fh, indent=2)

# CSV of the sweeps
import csv
with open("dualband_sweeps.csv", "w", newline="") as fh:
    wtr = csv.writer(fh)
    wtr.writerow(["design", "band", "f_MHz", "Re_Z", "Im_Z", "SWR"])
    for name in ("baseline", "optimised"):
        for band in ("2m", "70cm"):
            for r in results[name]["swr"][band]:
                wtr.writerow([name, band, r["f_MHz"], r["Re_Z"], r["Im_Z"], r["SWR"]])

print("\nWrote dualband_results.json and dualband_sweeps.csv")
