"""Diagnosis: compare (A) the validated reference approach (PW-const, point
matching, Pocklington K form), (B) the essay's current approach (PWM,
Galerkin, G+d2G/k^2 form), and (C) PWM + point matching + Pocklington K.
All for a centre-fed half-wave dipole at 144 MHz, a=3mm -> expect ~73+j42."""
import numpy as np

c = 299792458.0
mu0 = 4.0 * np.pi * 1e-7
eta0 = np.sqrt(mu0 / (8.8541878128e-12))

f = 144.0e6
k = 2*np.pi*f/c
w = 2*np.pi*f
lam = c/f
L = 0.5*lam
a = 3e-3
N = 100
D = L/N
edges = -L/2 + np.arange(N+1)*D
zc = -L/2 + (np.arange(N)+0.5)*D

def gauss(n):
    x, wt = np.polynomial.legendre.leggauss(n)
    return 0.5*x+0.5, 0.5*wt

def Kpock(dz, rho2, k):
    R = np.sqrt(rho2 + dz*dz)
    eR = np.exp(-1j*k*R)
    return eR*(dz*dz/R**3 - 2j*dz*dz/(k*R**4) - 2*dz*dz/(k**2*R**5)
               + 1j*rho2/(k*R**4) + rho2/(k**2*R**5) - 1.0/R)

# ---- (A) reference: PW-const, point matching, skill's exact formula -------
def ref_solve(N=N, a=a, L=L):
    D = L/N
    edges = -L/2+np.arange(N+1)*D
    zc = -L/2+(np.arange(N)+0.5)*D
    jwmu4p = 1j*w*mu0/(4*np.pi)
    Z = np.zeros((N,N), complex)
    for m in range(N):
        for n in range(N):
            rho2 = a*a
            z1, z2 = edges[n], edges[n+1]
            nq = 400 if m==n else 32
            xg, wq = gauss(nq)
            zq = 0.5*(z2-z1)*xg+0.5*(z1+z2)
            dz = zc[m]-zq
            K = Kpock(dz, rho2, k)
            Z[m,n] = jwmu4p*np.sum(wq*K)*D   # skill multiplies by dl_m
    V = np.zeros(N, complex); V[N//2] = 1.0
    I = np.linalg.solve(Z, V)
    return -1.0/I[N//2]

# ---- (B) essay current: PWM, Galerkin, G+d2G/k^2 --------------------------
def essay_solve(N=N, a=a, L=L, nq=48):
    D = L/N
    edges = -L/2+np.arange(N+1)*D
    xg, wq = gauss(nq)
    A = np.zeros((N,N), complex)
    for m in range(N):
        zm = edges[m]+(edges[m+1]-edges[m])*xg
        wmm = (edges[m+1]-edges[m])*wq*np.sin(k*(zm-edges[m]))
        for n in range(N):
            zn = edges[n]+(edges[n+1]-edges[n])*xg
            wnn = (edges[n+1]-edges[n])*wq*np.sin(k*(zn-edges[n]))
            u = zm[:,None]-zn[None,:]
            R = np.sqrt(a*a+u*u)
            e = np.exp(-1j*k*R)
            G = e/(4*np.pi*R)
            d2G = (e/(4*np.pi))*((2*u*u-a*a)/R**5 + 1j*k*(2*u*u-a*a)/R**4 + (1j*k)**2*u*u/R**3)
            K = (1j*k*eta0)*(G+d2G/k**2)
            A[m,n] = np.sum(K*wmm[:,None]*wnn[None,:])
    # unit feed current: replace feed row with f_center constraint
    m2 = N//2
    fc = np.sin(k*(0.0-edges[m2]))
    Ac = A.copy()
    row = np.zeros(N); row[m2] = fc
    Ac[m2,:] = row
    b = np.zeros(N, complex); b[m2] = 1.0
    I = np.linalg.solve(Ac, b)
    Zin = 1.0/(I[m2]*fc)   # V=1, I_feed = I_m f_m(0)
    return Zin

# ---- (C) PWM + point matching + Pocklington K -----------------------------
def pwm_point(N=N, a=a, L=L, nq_self=400, nq_off=32, sign=-1.0):
    D = L/N
    edges = -L/2+np.arange(N+1)*D
    zc = -L/2+(np.arange(N)+0.5)*D
    jwmu4p = 1j*w*mu0/(4*np.pi)
    A = np.zeros((N,N), complex)
    for m in range(N):
        for n in range(N):
            rho2 = a*a
            z1, z2 = edges[n], edges[n+1]
            nq = nq_self if m==n else nq_off
            xg, wq = gauss(nq)
            zq = 0.5*(z2-z1)*xg+0.5*(z1+z2)
            fn = np.sin(k*(zq-z1))
            dz = zc[m]-zq
            K = Kpock(dz, rho2, k)
            A[m,n] = jwmu4p*np.sum(wq*K*fn)   # point matching, PWM basis
    m2 = N//2
    fc = np.sin(k*(0.0-edges[m2]))
    Ac = A.copy()
    row = np.zeros(N); row[m2] = fc
    Ac[m2,:] = row
    b = np.zeros(N, complex); b[m2] = 1.0
    I = np.linalg.solve(Ac, b)
    Zin = sign*1.0/(I[m2]*fc)
    return Zin, I

print("half-wave dipole 144 MHz, a=3mm, N=100  (classic 73.1 + j42.5)")
try:
    ZA = ref_solve()
    print(f"(A) PW-const point-match (reference) : {ZA.real:8.2f} {ZA.imag:+8.2f} j")
except Exception as e:
    print("(A) failed:", repr(e))
try:
    ZB = essay_solve()
    print(f"(B) PWM Galerkin G+d2G (essay current): {ZB.real:8.2f} {ZB.imag:+8.2f} j")
except Exception as e:
    print("(B) failed:", repr(e))
for s in (-1.0, +1.0):
    ZC, I = pwm_point(sign=s)
    print(f"(C) PWM point-match K, sign={s:+.0f}   : {ZC.real:8.2f} {ZC.imag:+8.2f} j")
