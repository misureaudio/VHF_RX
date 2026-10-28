"""
antenna_mom.py
==============
Thin-wire Method-of-Moments (Pocklington EFIE) solver for co-planar,
parallel, z-directed PEC wire antennas.  numpy only.

Formulation
-----------
Each wire is a straight perfect-electric-conductor (PEC) filament along z at
transverse position (x0, y0), radius a, spanning z1..z2, split into N
segments of length dl = (z2-z1)/N.  The unknown is the z-directed surface
current, expanded in piecewise-constant (PW-const) basis functions::

    I(z) = sum_n I_n * chi_n(z),   chi_n = 1 on segment n, 0 elsewhere.

The electric-field integral equation is enforced tangentially (z) at each
segment centre, the test point being offset one radius a from the axis (on
the wire surface).  For a PEC::

    E_z^scat(r_m) = - E_z^inc(r_m).

The scattered z-field is the Pocklington kernel integral.  The scalar (G) and
vector (d2G/dz2) potential contributions are combined BEFORE evaluation so
the 1/R^2 near-field blow-ups cancel analytically.  With
G = e^{-jkR}/(4 pi R), R = sqrt(rho^2 + (z_m - z')^2), rho = transverse
distance (a for the self term, wire spacing for mutual terms)::

    E_z^scat(r_m) = (j w mu0 / 4 pi) * int I(z') K(z_m, z') dz'

    K = e^{-jkR} [ (z_m-z')^2/R^3 - 2j (z_m-z')^2/(k R^4)
                   - 2 (z_m-z')^2/(k^2 R^5)
                   + j rho^2/(k R^4) + rho^2/(k^2 R^5) - 1/R ]

There is NO separate endpoint / charge term: the scalar-potential
contribution is absorbed into the kernel integral via the fundamental theorem
of calculus (g(z1)-g(z2) = -int dg/dz' dz').  Adding a separate endpoint term
double-counts it and wrecks the reactance.

The kernel peaks sharply (width ~ wire radius a) at the centre of the SELF
segment, so the diagonal is integrated with fine Gauss-Legendre quadrature
(nq_self) while the smooth off-diagonal terms use coarse quadrature (nq_off).
A single coarse quadrature under-integrates the diagonal peak and collapses
Zin toward zero.

Why collocation (point matching), not Galerkin with this kernel
---------------------------------------------------------------
The kernel above is the *Pocklington* kernel: it is derived so that the EFIE
holds POINTWISE in the observation coordinate z.  It is therefore paired with
collocation (test at a point).  If one instead multiplies by a test function
f_m(z) and integrates over z (Galerkin), the (1/k^2) d2G/dz2 term must be
integrated by parts, which (i) replaces d2G/dz2 by f_m''(z) G and (ii)
produces endpoint terms [f_m dG/dz - f_m' G].  Feeding the un-integrated
Pocklington kernel directly into a double (Galerkin) integral skips both
steps and is mathematically wrong: it is exactly the defect that made the
earlier PWM/Galerkin version return Z ~ 0 - j1500 ohm for a half-wave dipole
instead of ~73 + j42 ohm.  This module therefore uses collocation.

Excitation (delta-gap model)
----------------------------
A port voltage V0 is applied across an infinitesimal gap at the feed segment:
E_z^inc = V0 delta(z - z_feed).  In the point-matching form the right-hand
side is b[feed] = V0 (the delta dispersed as V0/dl over the feed segment is
absorbed into the row scaling of the matrix).  The input impedance is::

    Z_in = -V0 / I_feed,

where I_feed is the solved current on the feed segment.  The minus sign comes
from E_scat = -E_inc.  (Using +V0 flips the sign of the radiation resistance;
a NEGATIVE Rr is the obvious red flag.)

Multi-wire (driven + parasitic)
-------------------------------
For a driven wire plus open-circuited parasitics the EFIE is enforced on ALL
segments.  The source is nonzero only at the driven wire's feed segment; the
parasitic segments have E_inc = 0 (PEC, open circuit), so their currents are
induced.  The feedpoint impedance of the driven wire is Z_in = -V0/I_feed
from the full-system solve.  For a true two-port (both wires fed) the
open-circuit Z-parameters come from the unit-current method (see
``z_parameters``) and a loaded 2-port reduces to
Z11' = Z11 - Z12 Z21/(Z22 + ZL).

Far field
---------
Coherent sum of Hertzian elements (E_phi = 0 for z-directed wires)::

    E_theta(th, ph) = (j eta0 k e^{-jkr}/4 pi r) sin(th)
        * sum_n I_n dl exp(-jk (x0 sin th cos ph + y0 sin th sin ph
                              + z_n cos th))

Directivity D = 4 pi U_max / int int U sin th dth dph with U = |E|^2;
gain = D * efficiency, efficiency from the conductor surface resistance.

Validation anchors (reproduce before trusting ANY downstream number)
--------------------------------------------------------------------
- Centre-fed 0.47 lambda, a = 3 mm dipole : Zin ~ 70 - j10 ohm (SWR ~1.46 @50).
- Centre-fed 0.50 lambda, a = 1 mm dipole : Zin ~ 81 + j27 ohm.
- Short dipole (feed-normalised sinusoidal current): Rr = 20 pi^2 (L/lambda)^2.
  (80 pi^2 (L/lambda)^2 is the UNIFORM-current Hertzian idealisation - do NOT
  compare the solver against it; that mismatch makes a 28 %-high result look
  "3x low".)
- Half-wave dipole directivity ~ 1.5 (2.15 dBi), E-plane 3-dB width ~ 78 deg.
"""
import numpy as np

# ---- physical constants -------------------------------------------------
c0 = 299792458.0
mu0 = 4.0e-7 * np.pi
eps0 = 8.8541878128e-12
eta0 = np.sqrt(mu0 / eps0)          # ~376.73 ohm
sigma_cu = 5.8e7                    # S/m (copper)

_gauss_cache = {}
def _gauss(n):
    if n not in _gauss_cache:
        _gauss_cache[n] = np.polynomial.legendre.leggauss(n)
    return _gauss_cache[n]


class Wire:
    """Straight PEC wire along z at (x0, y0), radius a, from z1 to z2 (m),
    split into N segments."""
    def __init__(self, x0, y0, z1, z2, a, N, feed=False):
        self.x0, self.y0, self.a = float(x0), float(y0), float(a)
        self.z1, self.z2 = float(z1), float(z2)
        self.N = int(N)
        self.dl = (self.z2 - self.z1) / self.N
        self.zc = self.z1 + (np.arange(self.N) + 0.5) * self.dl   # centres
        self.zs = self.z1 + np.arange(self.N + 1) * self.dl       # node z's
        self.feed = bool(feed)

    @property
    def L(self):
        return self.z2 - self.z1


def _global_segments(wires):
    segs = []
    for wi, wr in enumerate(wires):
        for i in range(wr.N):
            segs.append((wi, i, wr.zc[i], wr.x0, wr.y0, wr.a, wr.dl))
    return segs


def build_Z(wires, f, nq_off=24, nq_self=400):
    """MoM impedance matrix [ohm] (Pocklington EFIE, point matching,
    PW-const basis).  Fine quadrature on the self (diagonal) element, coarse
    off-diagonal.  See module docstring for the kernel and the
    no-endpoint-term rule."""
    k = 2.0 * np.pi * f / c0
    w = 2.0 * np.pi * f
    segs = _global_segments(wires)
    n = len(segs)
    Z = np.zeros((n, n), dtype=complex)
    jwmu4p = 1j * w * mu0 / (4.0 * np.pi)
    for m in range(n):
        wim, im, zcm, xm, ym, am, dlm = segs[m]
        tx, ty = xm + am, ym                 # test point on the surface (+x)
        for nn in range(n):
            win, ino, zcn, xn, yn, an, _ = segs[nn]
            rho2 = (tx - xn) ** 2 + (ty - yn) ** 2
            z1, z2 = wires[win].zs[ino], wires[win].zs[ino + 1]
            nq = nq_self if m == nn else nq_off
            xg, wg = _gauss(nq)
            zq = 0.5 * (z2 - z1) * xg + 0.5 * (z1 + z2)
            wq = 0.5 * (z2 - z1) * wg
            dz = zcm - zq
            R = np.sqrt(rho2 + dz ** 2)
            eR = np.exp(-1j * k * R)
            K = eR * (dz ** 2 / R ** 3 - 2j * dz ** 2 / (k * R ** 4)
                      - 2 * dz ** 2 / (k ** 2 * R ** 5)
                      + 1j * rho2 / (k * R ** 4) + rho2 / (k ** 2 * R ** 5)
                      - 1.0 / R)
            Z[m, nn] = jwmu4p * np.sum(wq * K) * dlm
    return Z


def solve(wires, f, feed_seg, V0=1.0, nq_off=24, nq_self=400):
    """Return (Z_in, I) for a delta-gap V0 at global segment feed_seg.

    Sign: E_scat = -E_inc  =>  Z_in = -V0 / I_feed.
    """
    Z = build_Z(wires, f, nq_off=nq_off, nq_self=nq_self)
    V = np.zeros(len(Z), dtype=complex)
    V[feed_seg] = V0
    I = np.linalg.solve(Z, V)
    Zin = -V0 / I[feed_seg]
    return Zin, I


def z_parameters(wires, f, feed_segs, nq_off=24, nq_self=400):
    """Open-circuit Z-parameters for ports at the given global feed segments.

    feed_segs: list of global segment indices (one per port).
    Returns (Z, Idists): Z the n-port open-circuit Z-matrix [ohm]; Idists[p]
    the full current vector when port p is driven by a unit feed current with
    all other ports open.
    """
    n = len(feed_segs)
    Zfull = build_Z(wires, f, nq_off=nq_off, nq_self=nq_self)
    Z = np.zeros((n, n), dtype=complex)
    Idists = []
    for p2 in range(n):
        A = Zfull.copy()
        b = np.zeros(len(Zfull), dtype=complex)
        # Drive port p2 with a unit feed current; OPEN every other port
        # (I_feed = 0) so the column is a true open-circuit Z-column.
        for p1 in range(n):
            row = np.zeros(len(Zfull))
            row[feed_segs[p1]] = 1.0
            A[feed_segs[p1], :] = row
            b[feed_segs[p1]] = 1.0 if p1 == p2 else 0.0
        I = np.linalg.solve(A, b)
        Idists.append(I)
        for p1 in range(n):
            # Sign consistent with solve(): the port voltage is the NEGATIVE
            # of the EFIE residual (E_scat = -E_inc), so Z = -Zfull @ I.
            Z[p1, p2] = -Zfull[feed_segs[p1], :] @ I
    return Z, Idists


def loaded_zin(Z, ZL):
    """Loaded input impedance of a 2-port: Z11' = Z11 - Z12 Z21/(Z22 + ZL)."""
    return Z[0, 0] - Z[0, 1] * Z[1, 0] / (Z[1, 1] + ZL)


# ---- far field ----------------------------------------------------------
def far_field(wires, I, f, ntheta=361, nphi=181):
    """Far field for z-directed wires.  Returns (E_theta, U, TH, PH) where
    E_theta is the complex E-field (relative, r = 1 m) and U = |E|^2 the power
    pattern on the (theta, phi) grid.  E_phi = 0."""
    k = 2.0 * np.pi * f / c0
    th = np.linspace(0.0, np.pi, ntheta)
    ph = np.linspace(0.0, 2.0 * np.pi, nphi, endpoint=False)
    TH, PH = np.meshgrid(th, ph, indexing='ij')
    st, ct = np.sin(TH), np.cos(TH)
    S = np.zeros_like(TH, dtype=complex)
    for wi, wr in enumerate(wires):
        Iw = I[wi * wr.N:(wi + 1) * wr.N]
        for i in range(wr.N):
            amp = Iw[i] * wr.dl
            p = (wr.x0 * st * np.cos(PH) + wr.y0 * st * np.sin(PH)
                 + wr.zc[i] * ct)
            S += amp * np.exp(-1j * k * p)
    Et = (1j * eta0 * k / (4.0 * np.pi)) * st * S
    U = np.abs(Et) ** 2
    return Et, U, TH, PH


def radiated_power(wires, I, f, ntheta=361, nphi=181):
    """Total radiated power [W] for a peak-phasor current I (P = 1/2 Re ...)."""
    _, U, TH, PH = far_field(wires, I, f, ntheta=ntheta, nphi=nphi)
    dth = TH[1, 0] - TH[0, 0]
    dph = PH[0, 1] - PH[0, 0]
    P = np.sum(U * np.sin(TH)) * dth * dph
    return 0.5 * P / eta0


def directivity(U, TH, PH):
    """Directivity [-] from a power pattern U on a (theta, phi) grid."""
    dth = TH[1, 0] - TH[0, 0]
    dph = PH[0, 1] - PH[0, 0]
    total = np.sum(U * np.sin(TH)) * dth * dph
    umax = U.max()
    return 4.0 * np.pi * umax / total


def efficiency(wires, I, f):
    """Radiating efficiency from the conductor surface resistance (copper)."""
    Rs = np.sqrt(np.pi * f * mu0 / sigma_cu)
    Ploss = 0.0
    for wi, wr in enumerate(wires):
        Iw = I[wi * wr.N:(wi + 1) * wr.N]
        Ploss += 0.5 * (Rs / (2.0 * np.pi * wr.a)) * np.sum(np.abs(Iw) ** 2) * wr.dl
    Prad = radiated_power(wires, I, f)
    Pin = Prad + Ploss
    return Prad / Pin if Pin > 0 else 1.0


def swr(Zin, Z0=50.0):
    g = (Zin - Z0) / (Zin + Z0)
    m = abs(g)
    if m >= 1.0:
        return float('inf')
    return (1.0 + m) / (1.0 - m)


def s11(Zin, Z0=50.0):
    g = (Zin - Z0) / (Zin + Z0)
    return g, 20.0 * np.log10(abs(g) if abs(g) > 0 else 1e-30)


def pattern_metrics(wires, I, f, ntheta=721, nphi=720):
    """E-plane / H-plane patterns, directivity, gain, 3-dB beamwidth, F/B.

    E-plane: elevation cut at phi=0 (the x-z plane), theta swept 0..180 deg.
    H-plane: AZIMUTHAL cut at theta=90 deg (broadside, the x-y plane),
    phi swept 0..360 deg.  For a z-directed dipole the H-plane is therefore
    (nearly) omnidirectional and must NOT show end-fire nulls; taking a second
    elevation cut (phi=90) and calling it the H-plane is the classic blunder.
    """
    Et, U, TH, PH = far_field(wires, I, f, ntheta=ntheta, nphi=nphi)
    D = directivity(U, TH, PH)
    eff = efficiency(wires, I, f)
    G = D * eff
    th = TH[:, 0]
    ph = PH[0, :]
    iph_e = 0
    Ue = U[:, iph_e]                       # E-plane: theta sweep at phi=0
    ang_e = th
    ith = int(np.argmin(np.abs(th - np.pi / 2)))
    Uh = U[ith, :]                         # H-plane: phi sweep at theta=90
    ang_h = ph

    def bw3d(Uplane, ang):
        a = np.sqrt(Uplane)
        i0 = int(np.argmax(a))
        half = a[i0] / np.sqrt(2.0)
        il, ir = i0, i0
        while il > 0 and a[il] > half:
            il -= 1
        while ir < len(a) - 1 and a[ir] > half:
            ir += 1
        return np.degrees(ang[ir] - ang[il])

    # F/B: forward (theta=90, phi=0) vs backward (theta=90, phi=180).
    # (Dividing by the whole-plane minimum is WRONG: it picks the endfire null
    #  at theta=0/180 and returns a nonsense ~250 dB.)
    fwd = U[ith, iph_e]
    bwd = U[ith, PH.shape[1] // 2]
    FBr = 10.0 * np.log10(fwd / max(bwd, 1e-30))
    return {
        "theta": th, "phi": ph, "Ue": Ue, "Uh": Uh,
        "ang_e": ang_e, "ang_h": ang_h,
        "D": D, "eff": eff, "G": G,
        "D_dBi": 10 * np.log10(D), "G_dBi": 10 * np.log10(G),
        "bw_E_deg": bw3d(Ue, ang_e), "bw_H_deg": bw3d(Uh, ang_h),
        "F/B_dB": FBr,
    }
