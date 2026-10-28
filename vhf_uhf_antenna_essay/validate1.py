"""Validation battery 1: kernel closed form + classical dipole impedances."""
import sys, time
import numpy as np
sys.path.insert(0, r"C:\Users\MATTIA\source\hermes-dir\vhf_uhf_antenna_essay")
from antenna_mom import (G_and_d2G, ParallelWireAntenna, c0, radiated_power,
                         pattern, bw_3dB)

lam = lambda f: c0 / f

# --- 1. closed-form d2G vs finite differences ------------------------------
k = 3.0
u = np.array([0.0, 0.3e-3, 1.0e-3, 5.0e-3, 1.0e-2, 0.05, 0.2])
d = 3.0e-3
h = 1e-6
_, d2G_cf = G_and_d2G(k, u, d)

def Gfun(u_):
    R = np.sqrt(d * d + u_ * u_)
    return np.exp(-1j * k * R) / (4 * np.pi * R)

d2G_fd = (Gfun(u + h) - 2 * Gfun(u) + Gfun(u - h)) / h ** 2
rel = np.abs(d2G_cf - d2G_fd) / np.maximum(np.abs(d2G_fd), 1e-30)
print("d2G closed form vs FD: max rel err =", rel.max())
assert rel.max() < 1e-4, "kernel closed form FAILED"

# --- 2. lambda/2 dipole, 144 MHz -------------------------------------------
f = 144.0e6
L = 0.5 * lam(f)
a = 3.0e-3

def dipole_Z(f, L, a=3e-3, N=100, nq=48):
    ant = ParallelWireAntenna([{"x": 0.0, "a": a, "L": L, "N": N, "feed": True}], f, nq=nq)
    Z, _ = ant.z_matrix([0])
    return Z[0, 0]

t0 = time.time()
Z = dipole_Z(f, L, N=100)
print(f"lam/2 dipole @144 MHz (N=100): Z = {Z.real:.2f} + j{Z.imag:.2f}  "
      f"({time.time()-t0:.1f}s)   [classic 73.1 + j42.5]")

# --- 3. convergence in N ----------------------------------------------------
print("\nconvergence (lam/2 @144 MHz):")
Zref = dipole_Z(f, L, N=200, nq=48)
print(f"  N= 200: Z = {Zref.real:8.3f} + j{Zref.imag:8.3f}   (reference)")
for N in (40, 80, 120):
    Z = dipole_Z(f, L, N=N)
    print(f"  N={N:4d}: Z = {Z.real:8.3f} + j{Z.imag:8.3f}   "
          f"|dZ vs N=200| = {abs(Z - Zref):.2e}")

# --- 4. resonant length scan ------------------------------------------------
print("\nresonant length scan (X=0):")
prev = None
found = None
for Lm in np.arange(0.40, 0.52, 0.002):
    Z = dipole_Z(f, Lm * lam(f), N=80)
    if prev is not None and np.sign(prev.imag) != np.sign(Z.imag):
        found = Lm
        break
    prev = Z
print(f"  resonant L/lambda ~ {found}   [classic 0.473]")

# --- 5. other classic impedances --------------------------------------------
print("\nclassic table check (e^{+jwt} convention, Balanis Table 5.1):")
for frac, (Rc, Xc) in {0.25: (18.8, -108.8), 0.5: (73.1, 42.5),
                       0.75: (243.8, 104.5), 1.0: (152.3, -46.1)}.items():
    Z = dipole_Z(f, frac * lam(f), N=100)
    print(f"  L/lambda={frac}: Z = {Z.real:8.1f} + j{Z.imag:8.1f}   "
          f"[classic {Rc} + j{Xc}]")

# --- 6. power balance --------------------------------------------------------
ant = ParallelWireAntenna([{"x": 0.0, "a": a, "L": L, "N": 100, "feed": True}], f)
I = ant.unit_current(0)
Zin = dipole_Z(f, L, N=100)
Pin = 0.5 * np.real(Zin)
P = radiated_power(ant, I)
print(f"\npower balance: Pin = {Pin:.4f} W, Prad = {P:.4f} W, "
      f"ratio = {P/Pin:.4f}")

# --- 7. gain / pattern -------------------------------------------------------
pat = pattern(ant, I)
print(f"gain = {pat['gain_dBi']:.2f} dBi  [classic 2.15 dBi], "
      f"bw3dB E-plane = {bw_3dB(pat['theta'], pat['Ee']):.1f} deg "
      f"[classic ~78 deg]")
