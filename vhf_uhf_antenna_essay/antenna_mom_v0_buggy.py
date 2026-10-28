"""
Thin-wire Method-of-Moments (Pocklington EFIE) solver for co-planar,
parallel, z-directed wires, all centered at z = 0.

Conventions
-----------
Time dependence e^{+j w t} (engineering convention, cf. Balanis).
Free-space Green's function  G(R) = exp(-j k R) / (4 pi R),  k = w / c.

From A_z = (mu0/4 pi) Int I(z') G dz', the Lorenz gauge and
E = -jw A - grad Phi one obtains the scattered tangential field on the
wire surface (observation point on the surface, source on the axis of a
parallel wire at transverse distance d):

    E_z^s(z) = -(j w mu0 / 4 pi) Int I(z') [ G + (1/k^2) d2G/dz2 ] dz'

The EFIE for a perfect conductor,  E_z^i + E_z^s = 0,  becomes

    (j w mu0 / 4 pi) Int I(z') [ G + (1/k^2) d2G/dz2 ] dz' = E_z^i(z).

With R = sqrt(d^2 + u^2), u = z - z', d = wire radius (self term) or wire
center spacing (mutual term for parallel wires), the second derivative in
closed form (verified against finite differences):

    d2G/dz^2 = e^{-jkR}/(4 pi) * [ (2u^2 - d^2)/R^5
                                   + jk (2u^2 - d^2)/R^4
                                   + (jk)^2 u^2 / R^3 ].

Basis: piecewise sinusoidal (PWM), f_n(z) = sin(k (z - z_{n-1})) on the
n-th segment [z_{n-1}, z_n].  The PWM is continuous at the junctions up
to O(k Delta), which vanishes with mesh refinement.

Ports (delta-gap model).  The source field of a port voltage V is
E_z^i = V delta(z) at the wire center.  The MoM right-hand side is
b_m = V f_m(0) on the center segment m; the feed current is
I(0) = I_m f_m(0);  Z_in = V / I(0).  For open-circuit Z-parameters,
port q is replaced by the constraint I_q(0) = 1 (the center-segment row
is replaced by that linear constraint) and the port voltage is the
scattered field point value V_p = E_z^s(0) at the feed point, summed
over all wires.

Loaded 2-ports are handled algebraically from the open-circuit matrix:
Z11' = Z11 - Z12 Z21 / (Z22 + ZL2), and current distributions are
superposed from the unit-current solutions.
"""
import numpy as np

c0 = 299_792_458.0
mu0 = 4.0e-7 * np.pi
eps0 = 8.8541878128e-12
eta0 = np.sqrt(mu0 / eps0)


def G_and_d2G(k, u, d):
    """Green's function and d2G/dz2 for observation-source separation
    u (longitudinal) and d (transverse).  Vectorized."""
    u = np.asarray(u, dtype=float)
    d = np.asarray(d, dtype=float)
    R = np.sqrt(d * d + u * u)
    e = np.exp(-1j * k * R)
    G = e / (4.0 * np.pi * R)
    d2G = (e / (4.0 * np.pi)) * (
        (2.0 * u * u - d * d) / R ** 5
        + 1j * k * (2.0 * u * u - d * d) / R ** 4
        + (1j * k) ** 2 * u * u / R ** 3
    )
    return G, d2G


def kernel(k, z, zp, d):
    """EFIE kernel (j w mu0) [ G + (1/k^2) d2G/dz2 ] = (j k eta0) [...].
    G includes 1/(4 pi).  Verified symbolically from A_z, Phi (Lorenz).
    z: (nz,), zp: (nzp,)."""
    G, d2G = G_and_d2G(k, z[:, None] - zp[None, :], d)
    return (1j * k * eta0) * (G + d2G / k ** 2)


class ParallelWireAntenna:
    """wires: list of dicts {x (m), a (radius m), L (m), N (segments),
    feed (bool)}; f: frequency (Hz)."""

    def __init__(self, wires, f, nq=48):
        self.wires = wires
        self.f = f
        self.w = 2.0 * np.pi * f
        self.k = self.w / c0
        self.nq = nq
        x, wt = np.polynomial.legendre.leggauss(nq)
        self.xq, self.wq = 0.5 * x + 0.5, 0.5 * wt
        self.edges = [(-w["L"] / 2.0 + np.arange(w["N"] + 1) * (w["L"] / w["N"]))
                      for w in wires]
        self.nseg = [w["N"] for w in wires]
        self.Ntot = int(sum(self.nseg))
        self.offsets = []
        off = 0
        for w in wires:
            self.offsets.append(off)
            off += w["N"]
        self._build_quadrature()

    def _build_quadrature(self):
        self.zpts = []    # (nseg, nq) quadrature z on each wire
        self.wts = []     # (nseg, nq) length weights
        self.fdiag = []   # (nseg, nq): f_n evaluated on its own segment
        for i, w in enumerate(self.wires):
            e = self.edges[i]
            zp, wt = [], []
            for n in range(w["N"]):
                z = e[n] + (e[n + 1] - e[n]) * self.xq
                zp.append(z)
                wt.append((e[n + 1] - e[n]) * self.wq)
            self.zpts.append(np.array(zp))
            self.wts.append(np.array(wt))
            self.fdiag.append(np.array(
                [np.sin(self.k * (z - e[n])) for n, z in enumerate(self.zpts[i])]))

    def center_seg(self, i):
        N = self.nseg[i]
        return N // 2 - 1 if N % 2 == 0 else N // 2

    def f_center(self, i):
        """f_m(0) for the center segment m of wire i."""
        m = self.center_seg(i)
        return np.sin(self.k * (0.0 - self.edges[i][m]))

    def matrix(self):
        A = np.zeros((self.Ntot, self.Ntot), dtype=complex)
        for i in range(len(self.wires)):
            wi = self.wires[i]
            for j in range(len(self.wires)):
                wj = self.wires[j]
                d = abs(wi["x"] - wj["x"]) if i != j else wi["a"]
                Zi = self.zpts[i]
                oi, oj = self.offsets[i], self.offsets[j]
                for m in range(self.nseg[i]):
                    wm = self.wts[i][m] * self.fdiag[i][m]
                    for n in range(self.nseg[j]):
                        K = kernel(self.k, Zi[m], self.zpts[j][n], d)
                        A[oi + m, oj + n] = np.sum(
                            K * wm[:, None]
                            * (self.wts[j][n] * self.fdiag[j][n])[None, :])
        return A

    def scattered_Ez(self, I, wire_idx, zobs=0.0):
        """Point value E_z^s on the surface of wire wire_idx at z = zobs
        (self term included)."""
        k = self.k
        E = 0.0
        for j, wj in enumerate(self.wires):
            d = abs(self.wires[wire_idx]["x"] - wj["x"]) \
                if j != wire_idx else self.wires[wire_idx]["a"]
            for n in range(self.nseg[j]):
                Ij = I[self.offsets[j] + n]
                if Ij == 0:
                    continue
                K = kernel(k, np.array([zobs]), self.zpts[j][n], d)[0]
                E += Ij * np.sum(K * self.wts[j][n])
        return -E

    def unit_current(self, port_idx, A=None):
        """Current vector when port (wire index) port_idx is driven by a
        unit feed current (open circuit elsewhere)."""
        A = self.matrix() if A is None else A.copy()
        b = np.zeros(self.Ntot, dtype=complex)
        i2 = port_idx
        m2 = self.center_seg(i2)
        row = np.zeros(self.Ntot)
        row[self.offsets[i2] + m2] = self.f_center(i2)
        A[self.offsets[i2] + m2, :] = row
        b[self.offsets[i2] + m2] = 1.0
        return np.linalg.solve(A, b)

    def z_matrix(self, ports):
        """Open-circuit Z-parameters for wire indices in `ports`.
        Port voltage: delta-gap model, V_p = (A I)_m / f_m(0) with m the
        center segment of wire p (same quadrature as the matrix).
        Returns (Z, Idists) with Idists[p] the unit-current distribution
        of port p."""
        n = len(ports)
        Z = np.zeros((n, n), dtype=complex)
        Idists = []
        A = self.matrix()
        for p2 in range(n):
            I = self.unit_current(ports[p2], A=A)
            Idists.append(I)
            for p1 in range(n):
                i1 = ports[p1]
                m1 = self.center_seg(i1)
                Z[p1, p2] = np.dot(A[self.offsets[i1] + m1, :], I) \
                    / self.f_center(i1)
        return Z, Idists

    def z_in(self, port_idx, loads=None):
        """Input impedance at port (wire index) port_idx with the other
        fed wires terminated in loads {wire_idx: Z_load} (exact 2-port
        algebra from the open-circuit matrix)."""
        ports = [port_idx] + [i for i in sorted(loads or {}) if i != port_idx]
        Z, _ = self.z_matrix(ports)
        Zin = Z[0, 0]
        if len(ports) > 1:
            ZL = loads[ports[1]]
            Zin = Z[0, 0] - Z[0, 1] * Z[1, 0] / (Z[1, 1] + ZL)
        return Zin

    def loaded_currents(self, port_idx, V, loads):
        """Current vector for source V at port_idx with other fed wires
        loaded per loads {wire_idx: Z_load}."""
        ports = [port_idx] + [i for i in sorted(loads or {}) if i != port_idx]
        Z, Idists = self.z_matrix(ports)
        I1 = V / (Z[0, 0] - Z[0, 1] * Z[1, 0] / (Z[1, 1] + loads[ports[1]])) \
            if len(ports) > 1 else V / Z[0, 0]
        I2 = -Z[1, 0] * I1 / (Z[1, 1] + loads[ports[1]]) if len(ports) > 1 else 0.0
        return I1 * Idists[0] + I2 * Idists[1]


def far_field(ant, I, theta, phi):
    """E_theta at r = 1 m for z-directed wires (E_phi = 0)."""
    k = ant.k
    th = np.atleast_1d(np.asarray(theta, dtype=float))
    ph = np.atleast_1d(np.asarray(phi, dtype=float))
    st, ct = np.sin(th), np.cos(th)
    sp, cp = np.sin(ph), np.cos(ph)
    E = np.zeros(th.shape, dtype=complex)
    for i, w in enumerate(ant.wires):
        e = ant.edges[i]
        posph = np.exp(1j * k * w["x"] * st * cp)
        M = np.zeros(th.shape, dtype=complex)
        for n in range(ant.nseg[i]):
            In = I[ant.offsets[i] + n]
            if In == 0:
                continue
            zq = ant.zpts[i][n]
            base = np.sin(k * (zq - e[n])) * ant.wts[i][n]        # (nq,)
            M += In * (np.exp(1j * k * zq[None, :] * ct[:, None]).T @ base)
        E += posph * M
    return (1j * eta0 * k / (4.0 * np.pi)) * np.sin(th) * E


def radiated_power(ant, I, n_theta=361, n_phi=181):
    th = np.linspace(0.0, np.pi, n_theta)
    ph = np.linspace(0.0, 2.0 * np.pi, n_phi, endpoint=False)
    P = 0.0
    dth = th[1] - th[0]
    dph = ph[1] - ph[0]
    for p in ph:
        Et = far_field(ant, I, th, np.full_like(th, p))
        P += np.trapezoid(np.abs(Et) ** 2 * np.sin(th), th) * dph
    return 0.5 * P / eta0


def pattern(ant, I, n_theta=721):
    th = np.linspace(0.0, np.pi, n_theta)
    Ee = far_field(ant, I, th, np.zeros_like(th))            # x-z plane
    Eh = far_field(ant, I, th, np.full_like(th, np.pi / 2))  # y-z plane
    P = radiated_power(ant, I)
    Umax = max(np.max(np.abs(Ee) ** 2), np.max(np.abs(Eh) ** 2))
    gain = 4.0 * np.pi * Umax / P
    return {"theta": th, "Ee": Ee, "Eh": Eh, "P_rad": P,
            "gain_dBi": 10.0 * np.log10(gain)}


def bw_3dB(theta, E):
    a = np.abs(E)
    i0 = int(np.argmax(a))
    half = a[i0] / np.sqrt(2.0)
    il = i0
    while il > 0 and a[il] > half:
        il -= 1
    ir = i0
    while ir < len(a) - 1 and a[ir] > half:
        ir += 1
    return np.degrees(theta[ir] - theta[il])


def s11(Zin, Z0=50.0):
    g = (Zin - Z0) / (Zin + Z0)
    return g, 20.0 * np.log10(abs(g))
