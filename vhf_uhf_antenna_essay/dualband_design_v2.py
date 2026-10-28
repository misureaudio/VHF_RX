"""Corrected dual-band VHF/UHF OPEN-SLEEVE dipole design.

Driven: 936 mm (2 m) at x=0, centre-fed.
Sleeves: a SYMMETRIC PAIR of 70 cm sleeves at x = +/- x7, close-coupled
(reactive near-field, ~0.03 lambda) so they act as an impedance transformer,
NOT a weakly-coupled Yagi parasitic.

Steps:
 1. Small optimisation: sleeve length x spacing, minimise 70 cm SWR @435 MHz.
 2. Band sweeps (2 m and 70 cm) for the chosen geometry.
 3. 3-port open-circuit Z (driven, sleeve+, sleeve-) at 435 MHz + loaded Zin.
 4. Pattern metrics at 145 and 435 MHz (corrected H-plane).
Writes dualband_v2_results.json.
"""
import json
import numpy as np
import antenna_mom as am

C = 299792458.0
L2, A2, N2 = 936.0e-3, 3.0e-3, 201
f2c, f7c = 145.0e6, 435.0e6
NQ_OFF, NQ_SELF = 24, 400

def build(L7, x7):
    N7 = max(41, int(round(L7 / (5.2e-3))))
    wires = [am.Wire(0.0, 0.0, -L2/2, L2/2, A2, N2, feed=True),
             am.Wire(x7, 0.0, -L7/2, L7/2, 3.0e-3, N7, feed=True),
             am.Wire(-x7, 0.0, -L7/2, L7/2, 3.0e-3, N7, feed=True)]
    return wires, N7

def zin70(L7, x7):
    wires, N7 = build(L7, x7)
    Zin, I = am.solve(wires, f7c, N2//2, V0=1.0, nq_off=NQ_OFF, nq_self=NQ_SELF)
    return Zin, wires, N7

print("=== Optimisation: sleeve length x spacing (70 cm SWR @435 MHz) ===")
best = None
for L7 in (300.0e-3, 316.0e-3, 332.0e-3):
    for x7 in (20.0e-3, 25.0e-3, 30.0e-3):
        Zin, _, _ = zin70(L7, x7)
        s = am.swr(Zin)
        tag = f"L7={L7*1e3:.0f}mm x7={x7*1e3:.0f}mm"
        print(f"  {tag:24s} Zin={Zin.real:6.1f} {Zin.imag:+7.1f}j  SWR={s:5.2f}")
        if best is None or s < best[0]:
            best = (s, L7, x7, Zin)
swr_best, L7b, x7b, _ = best
print(f"BEST: L7={L7b*1e3:.0f} mm, x7={x7b*1e3:.0f} mm, SWR={swr_best:.2f}")

wires, N7 = build(L7b, x7b)
feed = N2 // 2
res = {"geometry": {"L2_mm": L2*1e3, "L7_mm": L7b*1e3, "x7_mm": x7b*1e3,
                    "a_mm": 3.0, "N2": N2, "N7": N7, "sleeves": "symmetric pair"},
       "best_70cm_swr": swr_best}

print("\n=== Band sweeps ===")
for label, band in (("2m", range(144, 149)), ("70cm", range(430, 441, 2))):
    rows = []
    for fMHz in band:
        f = fMHz * 1e6
        Zin, _ = am.solve(wires, f, feed, V0=1.0, nq_off=NQ_OFF, nq_self=NQ_SELF)
        rows.append({"f_MHz": fMHz, "Re": round(Zin.real, 1),
                     "Im": round(Zin.imag, 1), "SWR": round(am.swr(Zin), 2)})
        print(f"  {label} {fMHz} MHz: {Zin.real:6.1f} {Zin.imag:+7.1f}j  SWR={am.swr(Zin):5.2f}")
    res[label] = rows

print("\n=== 3-port open-circuit Z @435 MHz (driven, sleeve+, sleeve-) ===")
n = len(wires)
segs = [w.N for w in wires]
feed_segs = [segs[0]//2, segs[0]+segs[1]//2, segs[0]+segs[1]+segs[2]//2]
Zfull = am.build_Z(wires, f7c, nq_off=NQ_OFF, nq_self=NQ_SELF)
Z3 = np.zeros((3, 3), dtype=complex)
for p2 in range(3):
    A = Zfull.copy()
    b = np.zeros(len(Zfull), dtype=complex)
    for p1 in range(3):
        row = np.zeros(len(Zfull)); row[feed_segs[p1]] = 1.0
        A[feed_segs[p1], :] = row
        b[feed_segs[p1]] = 1.0 if p1 == p2 else 0.0
    I = np.linalg.solve(A, b)
    for p1 in range(3):
        Z3[p1, p2] = -Zfull[feed_segs[p1], :] @ I
print(f"  Z11={Z3[0,0].real:.1f} {Z3[0,0].imag:+.1f}j")
print(f"  Z12={Z3[0,1].real:.1f} {Z3[0,1].imag:+.1f}j   |Z12|={abs(Z3[0,1]):.1f}")
print(f"  Z13={Z3[0,2].real:.1f} {Z3[0,2].imag:+.1f}j   |Z13|={abs(Z3[0,2]):.1f}")
print(f"  Z22={Z3[1,1].real:.1f} {Z3[1,1].imag:+.1f}j")
print(f"  Z23={Z3[1,2].real:.1f} {Z3[1,2].imag:+.1f}j")
# Loaded Zin of driven (ports 2,3 open): Zin' = Z11 - Z1v Zvv^-1 Zv1
Z1v = Z3[0, 1:3]
Zvv = Z3[1:3, 1:3]
Zv1 = Z3[1:3, 0]
Zin_loaded = Z3[0, 0] - Z1v @ np.linalg.solve(Zvv, Zv1)
print(f"  loaded Zin (sleeves open) = {Zin_loaded.real:.1f} {Zin_loaded.imag:+.1f}j  SWR={am.swr(Zin_loaded):.2f}")
res["zport_435"] = {
    "Z11": [Z3[0,0].real, Z3[0,0].imag], "Z12": [Z3[0,1].real, Z3[0,1].imag],
    "Z13": [Z3[0,2].real, Z3[0,2].imag], "Z22": [Z3[1,1].real, Z3[1,1].imag],
    "Z23": [Z3[1,2].real, Z3[1,2].imag],
    "absZ12": abs(Z3[0,1]),
    "loaded_Zin": [Zin_loaded.real, Zin_loaded.imag],
    "loaded_SWR": am.swr(Zin_loaded)}

print("\n=== Pattern metrics (corrected H-plane) ===")
for fMHz in (145, 435):
    f = fMHz * 1e6
    _, I = am.solve(wires, f, feed, V0=1.0, nq_off=NQ_OFF, nq_self=NQ_SELF)
    pm = am.pattern_metrics(wires, I, f)
    print(f"  {fMHz} MHz: D={pm['D_dBi']:.2f} dBi G={pm['G_dBi']:.2f} dBi "
          f"eff={pm['eff']*100:.1f}% bwE={pm['bw_E_deg']:.0f} bwH={pm['bw_H_deg']:.0f} "
          f"F/B={pm['F/B_dB']:.1f} dB")
    res[f"pattern_{fMHz}"] = {k: (round(v, 3) if isinstance(v, float) else v)
                              for k, v in pm.items()
                              if k not in ("theta", "phi", "Ue", "Uh", "ang_e", "ang_h")}

with open("dualband_v2_results.json", "w") as fh:
    json.dump(res, fh, indent=2)
print("\nwrote dualband_v2_results.json")
