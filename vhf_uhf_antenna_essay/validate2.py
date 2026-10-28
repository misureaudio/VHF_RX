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
# PWM current with I_n = 1 for all n:  I(z) = sin(k(z - edge_n)) on segment n
# EFIE static "voltage": V = j/(w eps0) * (1/4pi) * S / f_center0, S = sum_n int Ks fc fn
S = 0.0
Iz_total = 0.0
for nseg in range(N):
    zp = edges[nseg] + np.linspace(0, D, n)
    fp = np.sin(k*(zp - edges[nseg]))
    Iz_total += np.trapezoid(fp, zp)
    zm = edges[zc] + np.linspace(0, D, n)
    fm = np.sin(k*(zm - edges[zc]))
    U = zm[:,None] - zp[None,:]
    K = Ks(U)
    S += np.trapezoid(np.trapezoid(K*fm[:,None]*fp[None,:], zp, axis=1), zm)

V = 1j/(w*eps0)*(1/(4*np.pi))*S/f_center0
# Q = int rho dz, rho = (j/w) I'  ->  Q = (j/w) [I(L/2) - I(-L/2)] = 0 for a center-fed dipole!
# For a dipole the two halves carry opposite charge: Q_half = (j/w) int_0^{L/2} I' dz = (j/w)(I(L/2)-I(0)) = -(j/w) I(0)
I0 = np.sin(k*D/2)   # current at feed point (center of center segment)
Q_half = 1j/w*(0 - I0)
# energy: W = 1/2 sum over (half,half) pairs with |charge|: use full line with sign
# charge density rho(z) = (j/w) I'(z);  W = 1/2 int int rho(z) rho(z')/(4 pi R) dz dz'
zall = np.linspace(-L/2, L/2, 8000)
Iz = np.zeros_like(zall)
for nseg in range(N):
    m = (zall >= edges[nseg]) & (zall <= edges[nseg+1])
    Iz[m] = np.sin(k*(zall[m] - edges[nseg]))
# I' piecewise: k cos(k(z-edge)) on each segment (ignore junction deltas, O(kD) small)
Ip = np.zeros_like(zall)
for nseg in range(N):
    m = (zall > edges[nseg]) & (zall < edges[nseg+1])
    Ip[m] = k*np.cos(k*(zall[m] - edges[nseg]))
rho = 1j/w*Ip
Rm = np.sqrt(a*a + (zall[:,None]-zall[None,:])**2)
W = 0.5*np.trapezoid(
    np.trapezoid(np.real(rho*np.conj(rho))*0 + (rho*np.conj(rho).T*0), None) ) # placeholder
# W = 1/2 int dz int dz' rho(z) rho*(z') G_e
integrand = rho[:,None]*np.conj(rho[None,:])/(4*np.pi*Rm)
W = 0.5*np.trapezoid(np.trapezoid(integrand, zall, axis=1), zall).real
C_energy = abs(Q_half)**2/(W*2)  # W = Q^2/(2C) for each half (C_half); dipole C = 2*C_half
C_half_energy = abs(Q_half)**2/(2*W_half) if False else None
# compute W for one half only
mh = zall >= 0
zh = zall[mh]; rh = rho[mh]
Rh = np.sqrt(a*a + (zh[:,None]-zh[None,:])**2)
Wh = 0.5*np.trapezoid(np.trapezoid(rh[:,None]*np.conj(rh[None,:])/(4*np.pi*Rh), zh, axis=1), zh).real
C_half_energy = abs(Q_half)**2/(2*Wh)
C_dipole_energy = 2*C_half_energy

# EFIE capacitance: V relates to I0 via V = I0/(j w C_eff)  ->  C_eff = I0/(j w V)
C_EFIE = I0/(1j*w*V)
print(f"I0 = {I0:.3e}, S = {S:.4e}")
print(f"V = {V.real:.4e} + j{V.imag:.4e}")
print(f"C_EFIE  = {C_EFIE.real*1e12:.4f} pF   (imag {C_EFIE.imag*1e12:.2e})")
print(f"C_dipole (energy) = {C_dipole_energy*1e12:.4f} pF")
print(f"ratio C_EFIE/C_energy = {C_EFIE.real/C_dipole_energy:.3f}")