"""Static-limit cross-check: EFIE capacitance operator vs direct electrostatics.

PWM current I_n = 1 on all segments (sawtooth).  Compares
  C_EFIE    from the static EFIE kernel  (2u^2-a^2)/R^5 / (4 pi)
  C_energy  from W = 1/2 int int rho rho' G_e dz dz',  rho = (j/w) I'
For a correct kernel these must agree (EFIE is the variational
electrostatics; the continuity equation supplies rho).
"""
import numpy as np
from antenna_mom import c0, eps0

f = 5e6
k = 2*np.pi*f/c0
w = 2*np.pi*f
L = 0.02*c0/144e6
a = 3e-3
N = 40
D = L/N
edges = -L/2 + np.arange(N+1)*D
zc = N//2
f_center0 = np.sin(k*D/2)


def Ks(u):
    R = np.sqrt(a*a + u*u)
    return (2*u*u - a*a)/(4*np.pi*R**5)


n = 800
S = 0.0
for nseg in range(N):
    zp = edges[nseg] + np.linspace(0, D, n)
    fp = np.sin(k*(zp - edges[nseg]))
    zm = edges[zc] + np.linspace(0, D, n)
    fm = np.sin(k*(zm - edges[zc]))
    U = zm[:, None] - zp[None, :]
    K = Ks(U)
    S += np.trapezoid(np.trapezoid(K*fm[:, None]*fp[None, :], zp, axis=1), zm)

V = 1j/(w*eps0)*(1/(4*np.pi))*S/f_center0
I0 = np.sin(k*D/2)
C_EFIE = I0/(1j*w*V)

# direct electrostatics on one half
n2 = 8000
zall = np.linspace(-L/2, L/2, n2)
Iz = np.zeros_like(zall)
Ip = np.zeros_like(zall)
for nseg in range(N):
    m1 = (zall >= edges[nseg])*(zall <= edges[nseg+1])
    m2 = (zall > edges[nseg])*(zall < edges[nseg+1])
    Iz[m1] = np.sin(k*(zall[m1] - edges[nseg]))
    Ip[m2] = k*np.cos(k*(zall[m2] - edges[nseg]))
rho = 1j/w*Ip
Rfull = np.sqrt(a*a + (zall[:, None]-zall[None, :])**2)
Wtot = 0.5*np.trapezoid(
    np.trapezoid(rho[:, None]*np.conj(rho[None, :])/(4*np.pi*Rfull),
                 zall, axis=1), zall).real
Q_half = 1j/w*(0.0 - I0)
C_dipole_energy = abs(Q_half)**2/(2*Wtot)

print(f"I0 = {I0:.3e}, S = {S:.4e}")
print(f"V = {V.real:.4e} + j{V.imag:.4e}")
print(f"C_EFIE        = {C_EFIE.real*1e12:.4f} pF   (imag {C_EFIE.imag*1e12:.2e})")
print(f"C_dipole (energy) = {C_dipole_energy*1e12:.4f} pF")
print(f"ratio C_EFIE/C_energy = {C_EFIE.real/C_dipole_energy:.3f}")
print(f"spheroid reference  = {2*np.pi*eps0*L/np.log(L/(2*a))*1e12:.4f} pF")
