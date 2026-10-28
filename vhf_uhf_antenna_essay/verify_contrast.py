"""Verify the two contrast-table entries not in the design optimisation grid:
 (1) a SINGLE offset 70 cm sleeve at x=+25 mm (2-wire system);
 (2) a SYMMETRIC PAIR at x=+/-10 mm (3-wire system).
Both at 435 MHz. Report loaded Zin and SWR, plus |Z12| for the single case.
"""
import numpy as np
import antenna_mom as am

C = 299792458.0
L2, A2, N2 = 936.0e-3, 3.0e-3, 201
L7, A7 = 316.0e-3, 3.0e-3
f7c = 435.0e6
NQ_OFF, NQ_SELF = 24, 400

def loaded_zin(wires, feed, f):
    """Input impedance of `feed` port with all other ports open-circuited."""
    n = len(wires)
    segs = [w.N for w in wires]
    feed_segs = []
    s = 0
    for w in wires:
        feed_segs.append(s + w.N // 2)
        s += w.N
    Zfull = am.build_Z(wires, f, nq_off=NQ_OFF, nq_self=NQ_SELF)
    Zp = np.zeros((n, n), dtype=complex)
    for p2 in range(n):
        A = Zfull.copy()
        b = np.zeros(len(Zfull), dtype=complex)
        for p1 in range(n):
            row = np.zeros(len(Zfull)); row[feed_segs[p1]] = 1.0
            A[feed_segs[p1], :] = row
            b[feed_segs[p1]] = 1.0 if p1 == p2 else 0.0
        I = np.linalg.solve(A, b)
        for p1 in range(n):
            Zp[p1, p2] = -Zfull[feed_segs[p1], :] @ I
    Zin_loaded = Zp[0, 0] - Zp[0, 1:] @ np.linalg.solve(Zp[1:, 1:], Zp[1:, 0])
    return Zp, Zin_loaded

# (1) Single offset sleeve at x=+25 mm
N7 = 61
w1 = [am.Wire(0.0, 0.0, -L2/2, L2/2, A2, N2, feed=True),
      am.Wire(25.0e-3, 0.0, -L7/2, L7/2, A7, N7, feed=True)]
Zp1, Zin1 = loaded_zin(w1, N2//2, f7c)
print(f"single +25mm : Zin={Zin1.real:6.1f} {Zin1.imag:+7.1f}j  SWR={am.swr(Zin1):5.2f}   "
      f"|Z12|={abs(Zp1[0,1]):5.1f}  Z11={Zp1[0,0].real:.1f}{Zp1[0,0].imag:+.1f}j")

# (2) Symmetric pair at x=+/-10 mm
w2 = [am.Wire(0.0, 0.0, -L2/2, L2/2, A2, N2, feed=True),
      am.Wire(10.0e-3, 0.0, -L7/2, L7/2, A7, N7, feed=True),
      am.Wire(-10.0e-3, 0.0, -L7/2, L7/2, A7, N7, feed=True)]
Zp2, Zin2 = loaded_zin(w2, N2//2, f7c)
print(f"pair  +/-10mm: Zin={Zin2.real:6.1f} {Zin2.imag:+7.1f}j  SWR={am.swr(Zin2):5.2f}   "
      f"|Z12|={abs(Zp2[0,1]):5.1f}  |Z13|={abs(Zp2[0,2]):5.1f}")

# (3) Symmetric pair at x=+/-20 mm (for the |Z12| column of the contrast table)
w3 = [am.Wire(0.0, 0.0, -L2/2, L2/2, A2, N2, feed=True),
      am.Wire(20.0e-3, 0.0, -L7/2, L7/2, A7, N7, feed=True),
      am.Wire(-20.0e-3, 0.0, -L7/2, L7/2, A7, N7, feed=True)]
Zp3, Zin3 = loaded_zin(w3, N2//2, f7c)
print(f"pair  +/-20mm: Zin={Zin3.real:6.1f} {Zin3.imag:+7.1f}j  SWR={am.swr(Zin3):5.2f}   "
      f"|Z12|={abs(Zp3[0,1]):5.1f}  |Z13|={abs(Zp3[0,2]):5.1f}")

# (4) Single offset sleeve at x=+15 mm (the 'pulls the real part down' claim)
w4 = [am.Wire(0.0, 0.0, -L2/2, L2/2, A2, N2, feed=True),
      am.Wire(15.0e-3, 0.0, -L7/2, L7/2, A7, N7, feed=True)]
Zp4, Zin4 = loaded_zin(w4, N2//2, f7c)
print(f"single +15mm : Zin={Zin4.real:6.1f} {Zin4.imag:+7.1f}j  SWR={am.swr(Zin4):5.2f}   "
      f"|Z12|={abs(Zp4[0,1]):5.1f}")

# (5) Single offset sleeve at x=+307 mm (far-field / Yagi-like, the contrast anchor)
w5 = [am.Wire(0.0, 0.0, -L2/2, L2/2, A2, N2, feed=True),
      am.Wire(307.0e-3, 0.0, -L7/2, L7/2, A7, N7, feed=True)]
Zp5, Zin5 = loaded_zin(w5, N2//2, f7c)
print(f"single +307mm: Zin={Zin5.real:6.1f} {Zin5.imag:+7.1f}j  SWR={am.swr(Zin5):5.2f}   "
      f"|Z12|={abs(Zp5[0,1]):5.1f}")
