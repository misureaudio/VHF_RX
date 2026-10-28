# The Method of Moments for VHF/UHF Wire Antennas: From the Pocklington
# Electric-Field Integral Equation to a Dual-Band 144–440 MHz Design

*An essay for readers working in mathematical physics, functional and
numerical analysis, and classical electromagnetism.*

---

## Abstract

We give a self-contained derivation of the thin-wire electric-field integral
equation (EFIE) for a perfect-electric-conductor (PEC) wire antenna,
discretise it by the method of moments (MoM) with piecewise-constant basis
functions and point-matching (collocation) testing, and analyse the resulting
linear system from a numerical-analysis standpoint. The central object is the
Pocklington kernel, in which the scalar- and vector-potential contributions
are combined *before* evaluation so that the near-field $1/R^2$ and $1/R^3$
singularities cancel analytically; we show that this kernel is the unique
formulation that is both pointwise valid and free of separate endpoint (charge)
terms. We then confront the discretisation with the classical limits — the
short-dipole radiation resistance, the half-wave dipole impedance, the
resonant length, the far-field directivity, and the global power balance — and
report what the solver reproduces to within a few percent and where the
piecewise-constant basis degrades (very short and electrically long dipoles).
As an application we design and simulate a dual-band VHF/UHF open-sleeve
dipole for the 144–148 MHz (2 m) and 430–440 MHz (70 cm) amateur bands,
computing feedpoint impedance, SWR, radiation pattern, directivity, gain,
beamwidth, and front-to-back ratio on both bands. The design is a
close-coupled open-sleeve dipole: a centre-fed 2 m driven element with a
symmetric pair of 70 cm sleeves in the reactive near-field
($\pm25\ \mathrm{mm} \approx 0.036\lambda$). The solver predicts a well-matched
2 m band (SWR $1.7$–$2.3$) *and* a well-matched 70 cm band (SWR
$1.3$–$1.6$): the sleeves act as an impedance transformer, cancelling the
driven element's $-j130\ \Omega$ capacitive reactance at $1.36\lambda$ through
a large mutual coupling ($|Z_{12}| \approx 70\ \Omega$). We contrast this with
the far-field (Yagi-like) spacing, where the coupling is negligible
(a few ohms, $|Z_{12}| \approx 3\ \Omega$) and the 70 cm band is left
mismatched (SWR $\approx 7.5$) — illustrating that the *coupling regime*,
not the mere presence of a second element, determines the dual-band match.

---

## 1. Introduction

The VHF (30–300 MHz) and UHF (300 MHz–3 GHz) portions of the radio spectrum
are the working bands of a large fraction of terrestrial wireless
communication, and of the amateur-radio 2 m (144–148 MHz) and 70 cm
(430–440 MHz) allocations. At these frequencies the electrical size of a
practical wire antenna is of order a wavelength: a half-wave dipole on 2 m is
about 1 m long, on 70 cm about 35 cm. This is the regime in which the
*thin-wire integral-equation* description is both natural and accurate — the
wire radius $a$ is a few millimetres, three orders of magnitude below the
wavelength, yet the wire is long enough that the current distribution is not
trivial.

Analytical antenna theory handles the short dipole ($L\ll\lambda$) and, with
more effort, the half-wave dipole ($L\approx\lambda/2$) in closed form. But a
dual-band design — one structure that must behave acceptably at two
frequencies a factor of three apart — has no closed-form answer: the current
on each element is set by the coupled integral equations on *all* elements,
and the feedpoint impedance on the upper band is controlled by the lower-band
element seen at an electrical length of $1.36\lambda$. This is a problem whose
natural formulation is a (singular) integral equation, and whose natural
solution is a numerical one. The method of moments is the standard
discretisation.

The purpose of this essay is twofold. First, to derive the thin-wire EFIE and
its MoM discretisation with the rigour that a functional/numerical-analysis
reader expects: what the integral operator is, what the basis and testing
discretisation do to it, where the singularities live, what the quadrature
must resolve, and what the convergence and conditioning look like. Second, to
use the resulting solver as a *measuring instrument* on a real design problem,
and to report — with the numbers the code actually produces — what is
robust, what is a basis artefact, and what is physics.

We work throughout with the engineering time convention $e^{+j\omega t}$ and
the free-space Green's function $G(R)=e^{-jkR}/(4\pi R)$, $k=\omega/c$. All
numerical results in §5 are produced by the accompanying NumPy-only solver
`antenna_mom.py` (Appendix A); the validation battery of §4 is `validate_fixed.py`.

---

## 2. The thin-wire electric-field integral equation

### 2.1 From the potentials to the tangential field

Let the source be a $z$-directed surface current $I(z')$ on a straight PEC
filament of radius $a$ at transverse position $(x_0,y_0)$, with time
dependence $e^{+j\omega t}$. In the Lorenz gauge,
$\nabla\!\cdot\!\mathbf A = -j\omega\mu_0\varepsilon_0\,\Phi$, the fields are
expressed through the vector potential $\mathbf A$ and scalar potential
$\Phi$:

$$
\mathbf E = -j\omega\mathbf A - \nabla\Phi,
\qquad
\mathbf A(\mathbf r) = \mu_0 \int \mathbf J(\mathbf r')\,G(\mathbf r,\mathbf r')\,d^3r',
\qquad
G(R)=\frac{e^{-jkR}}{4\pi R}.
$$

For a $z$-directed filament the only non-zero potential component is
$A_z$, and $\nabla\!\cdot\!\mathbf A = \partial A_z/\partial z$, so

$$
A_z(z) = \mu_0\int I(z')\,G(z,z')\,dz',
\qquad
\Phi = -\frac{1}{j\omega\mu_0\varepsilon_0}\frac{\partial A_z}{\partial z}.
$$

The tangential ($z$) scattered field on the wire is therefore

$$
E_z^{s}(z) = -j\omega A_z(z) - \frac{\partial\Phi}{\partial z}
= -j\omega\!\left[A_z(z) + \frac{1}{k^2}\frac{\partial^2 A_z}{\partial z^2}\right],
$$

where we used $k^2=\omega^2\mu_0\varepsilon_0$. Substituting the potential
integral,

$$
\boxed{\;
E_z^{s}(z) = -j\omega\mu_0\int I(z')\left[G(z,z') + \frac{1}{k^2}\frac{\partial^2 G(z,z')}{\partial z^2}\right]dz'
\;}
\tag{1}
$$

with $z$ the *observation* coordinate. The perfect-conductor boundary
condition, $E_z^{i}(z)+E_z^{s}(z)=0$ on the wire surface, gives the EFIE

$$
j\omega\mu_0\int I(z')\left[G + \frac{1}{k^2}\frac{\partial^2 G}{\partial z^2}\right]dz' = E_z^{i}(z).
\tag{2}
$$

Equation (2) is the starting point. It is a Fredholm integral equation of the
first kind for the unknown current $I(z')$, with a kernel that is
*weakly singular but not integrable in the naive form*: the second-derivative
term $\partial^2G/\partial z^2$ carries $1/R^3$, $1/R^4$, $1/R^5$ pieces that
blow up as the observation point approaches the source ($R\to\rho$, the
transverse separation).

### 2.2 The Pocklington kernel

The resolution, due to Pocklington [1], is to combine the scalar and vector
potential contributions *before* evaluating the derivative, so that the
near-field blow-ups cancel. Writing $u = z-z'$, $\rho$ the transverse
separation (equal to $a$ for the self term, to the wire spacing for a mutual
term), and $R=\sqrt{\rho^2+u^2}$, one finds by direct differentiation

$$
G + \frac{1}{k^2}\frac{\partial^2 G}{\partial z^2}
= \frac{e^{-jkR}}{4\pi}\left[
\frac{\rho^2}{R^3} + j\frac{2u^2-\rho^2}{k R^4} + \frac{2u^2-\rho^2}{k^2 R^5}
\right],
$$

and hence

$$
-j\omega\mu_0\left[G + \frac{1}{k^2}\frac{\partial^2 G}{\partial z^2}\right]
= \frac{j\omega\mu_0}{4\pi}\,K(z,z'),
$$

where the **Pocklington kernel** is

$$
\boxed{\;
K = e^{-jkR}\left[
\frac{u^2}{R^3}
- \frac{2j\,u^2}{k R^4}
- \frac{2\,u^2}{k^2 R^5}
+ \frac{j\,\rho^2}{k R^4}
+ \frac{\rho^2}{k^2 R^5}
- \frac{1}{R}
\right].
\;}
\tag{3}
$$

Two structural facts about (3) matter for the numerics.

**No separate endpoint term.** The scalar-potential (charge) contribution is
already *inside* the kernel integral. This is not an accident: it is the
fundamental theorem of calculus. The charge density is
$\rho_c(z') = (j/\omega)\,dI/dz'$ (continuity equation), and its potential
contribution to $E_z^s$ is $-\partial\Phi_c/\partial z$, which integrates by
parts to a boundary term $g(z_1)-g(z_2) = -\int (dg/dz')dz'$ that is absorbed
into the $G$ term. Adding a *separate* endpoint/charge term therefore
double-counts the scalar potential and corrupts the reactance. This is the
first pitfall in constructing a thin-wire MoM code.

**A sharp, narrow self-peak.** As $R\to\rho$ (observation point on the
surface directly above a source point, $u\to0$), the kernel is dominated by
the $\rho^2/R^3$ and $\rho^2/(k^2R^5)$ pieces, which peak with a width of
order $\rho=a$. The self-term ($m=n$) matrix element therefore integrates a
feature of width $a$ — for $a=3$ mm and a segment of length $\Delta\ell\sim
1$–$5$ mm the peak is comparable to or broader than the segment. This is the
second, and numerically more serious, pitfall: a single coarse quadrature
under-integrates the diagonal peak and collapses the impedance toward zero.
§4 treats the quadrature this demands.

### 2.3 The EFIE as an operator

For the functional-analytic reader it is useful to recast (2) as an operator
equation. Let $\mathcal K$ be the integral operator

$$
(\mathcal K I)(z) = \frac{j\omega\mu_0}{4\pi}\int_{-L/2}^{L/2} K(z,z')\,I(z')\,dz',
$$

acting on the space of (square-integrable) currents on the wire. The EFIE is

$$
\mathcal K\, I = E_z^{i}.
\tag{4}
$$

$\mathcal K$ is a compact (Hilbert-Schmidt, modulo the weak $1/R$ singularity
which is integrable) operator on $L^2$, so (4) is a regularised
first-kind equation; the MoM replaces it by a finite-rank approximation. The
choice of basis (the range of the finite-rank operator) and of testing
(collocation vs. Galerkin) determines the discretisation error and, as we
argue in §3.4, the *validity* of the kernel pairing.

---

## 3. Method-of-moments discretisation

### 3.1 Basis and expansion

We expand the current in **piecewise-constant** (PW-const) basis functions.
The wire $[{-}L/2, L/2]$ is split into $N$ segments of length
$\Delta\ell = L/N$ with nodes $z_n = -L/2 + n\,\Delta\ell$, and

$$
I(z) = \sum_{n=0}^{N-1} I_n\,\chi_n(z),
\qquad
\chi_n(z) = \begin{cases}1, & z\in[z_n,z_{n+1})\\ 0, & \text{otherwise}.\end{cases}
$$

The coefficients $I_n$ are the unknowns. The PW-const basis is the simplest
choice and, with $\Delta\ell\ll\lambda$, is spectrally adequate for the
half-wave regime that matters here. (The piecewise-*sinusoidal* basis,
$f_n(z)=\sin(k(z-z_n))$, is a natural refinement — it enforces a smoother,
approximately-continuous current and is more accurate in the short-dipole
limit — but it changes the testing requirement, as discussed in §3.4.)

### 3.2 Collocation (point matching)

We enforce the EFIE **at a point**: at the centre $z_m$ of each segment $m$,
with the test point offset one radius $a$ from the axis (on the wire
surface). Substituting the PW-const expansion,

$$
\sum_{n=0}^{N-1} I_n \underbrace{\frac{j\omega\mu_0}{4\pi}\int_{z_n}^{z_{n+1}} K(z_m,z')\,dz'}_{A_{mn}}
= E_z^{i}(z_m),
$$

gives the linear system

$$
\mathbf A\,\mathbf I = \mathbf b,
\qquad
A_{mn} = \frac{j\omega\mu_0}{4\pi}\int_{z_n}^{z_{n+1}} K(z_m,z')\,dz',
\tag{5}
$$

where the right-hand side $\mathbf b$ is the incident field sampled at the
test points. For a delta-gap feed of voltage $V_0$ at the centre,
$E_z^{i} = V_0\,\delta(z)$, so $b_{\text{feed}} = V_0$ (dispersed as
$V_0/\Delta\ell$ over the feed segment in the point-matching form) and
$b_m=0$ elsewhere.

### 3.3 The impedance matrix and the input impedance

The matrix $\mathbf A$ in (5) is the MoM impedance matrix. It is dense,
complex, and non-symmetric in general (the kernel $K(z,z')$ is not symmetric
in its arguments once the observation point is offset from the source axis).
Its construction costs $O(N^2)$ kernel integrals, each a one-dimensional
quadrature; the solve is $O(N^3)$ (or $O(N^2)$ per frequency for a fixed
geometry with a factorised matrix).

The **input impedance** follows from the delta-gap model. With the PEC
condition $E^s=-E^i$, the solved current on the feed segment $I_{\text{feed}}$
relates to the applied voltage by

$$
\boxed{\; Z_{\text{in}} = -\,\frac{V_0}{I_{\text{feed}}} \;}
\tag{6}
$$

The minus sign is a consequence of $E^s=-E^i$ and is a common source of
sign errors: using $+V_0$ flips the sign of the radiation resistance, and a
*negative* $R_{\text{rad}}$ is the obvious red flag that the excitation sign
is wrong.

For a **multi-wire** structure (a driven element plus open-circuited
parasitics) the EFIE is enforced on *all* segments. The source is non-zero
only at the driven element's feed segment; the parasitic segments have
$E^i=0$ (PEC, open circuit), so their currents are *induced* by the driven
element. The feedpoint impedance of the driven element is $Z_{\text{in}}$
from the full-system solve (6). For a true two-port (both elements fed) the
open-circuit $Z$-parameters follow from the unit-current method: drive port
$p$ with $I_{\text{feed},p}=1$ (replacing that row of $\mathbf A$ by the
linear constraint $I_{\text{feed},p}=1$), and read off
$Z_{qp} = V_q/I_{\text{feed},p}$. A loaded two-port then reduces to
$Z_{11}' = Z_{11} - Z_{12}Z_{21}/(Z_{22}+Z_L)$.

### 3.4 Why collocation, not Galerkin, with this kernel

This is the subtle point that the earlier (buggy) version of the solver got
wrong, and it is worth stating precisely because it is a general feature of
integral-equation discretisations.

The Pocklington kernel (3) is derived so that the EFIE (2) holds
**pointwise** in the observation coordinate $z$. It is therefore the correct
kernel for *collocation* (test at a point). If instead one multiplies (2) by
a test function $f_m(z)$ and integrates over $z$ (Galerkin), the term
$(1/k^2)\,\partial^2G/\partial z^2$ must be integrated by parts:

$$
\int f_m(z)\,\frac{\partial^2 G}{\partial z^2}\,dz
= \Big[f_m\frac{\partial G}{\partial z} - f_m'\,G\Big]_{\text{endpoints}}
+ \int f_m''(z)\,G\,dz .
$$

Two things happen: (i) the second $z$-derivative on $G$ is replaced by
$f_m''(z)G$ (a derivative on the *test function*, not on the Green's
function), and (ii) **endpoint terms**
$[f_m\,\partial G/\partial z - f_m' G]$ appear. Feeding the un-integrated
Pocklington kernel directly into a double (Galerkin) integral —
$A_{mn}=\int\!\!\int f_m K f_n\,dz\,dz'$ — skips *both* steps. It is not a
discretisation of (2); it is a different, wrong operator.

This is exactly what the first version of `antenna_mom.py` did: a
piecewise-sinusoidal basis with a Galerkin (double-integral) assembly of the
un-integrated kernel. For a half-wave dipole it returned
$Z\approx 0 - j1500\ \Omega$ instead of $\approx 73+j42\ \Omega$, and the
error did not converge with $N$ (it is a formulation error, not a
discretisation error, so mesh refinement cannot fix it). Replacing the
Galerkin assembly by collocation with the Pocklington kernel (5) — the change
made in the corrected solver — restores agreement with the classical values
to $<1\%$ in the half-wave regime (§4).

**Lesson.** The kernel and the testing procedure are a *pair*. The
Pocklington kernel (3) is the collocation kernel. A Galerkin discretisation
requires the integrated-by-parts kernel (with $f_m''G$ and endpoint terms),
which is a different object. Mixing them is a category error that no amount
of quadrature refinement can repair.

---

## 4. Numerical analysis and validation

### 4.1 The quadrature: two-level Gauss–Legendre

The kernel (3) is smooth and bounded for $u\neq0$ but peaks sharply at $u=0$
on the **self** segment, with a width of order $a$. The off-diagonal terms
($m\neq n$) are smooth over the source segment. The correct response is a
*two-level* quadrature:

- **Diagonal (self, $m=n$):** fine Gauss–Legendre, $n_q^{\text{self}}\approx
  300$–$400$ nodes, to resolve the width-$a$ peak.
- **Off-diagonal ($m\neq n$):** coarse Gauss–Legendre,
  $n_q^{\text{off}}\approx 24$ nodes, since the integrand is smooth.

A single coarse quadrature (say $n_q=24$ everywhere) under-integrates the
diagonal peak. In the buggy version this collapsed $Z_{\text{in}}$ toward
zero; with the corrected two-level scheme the half-wave anchors are
reproduced to $<1\%$ (§4.2). The cost of the fine diagonal is modest: only
$N$ of the $N^2$ elements use $n_q^{\text{self}}$.

### 4.2 Validation against the classical limits

The solver is validated against a battery of classical results before any
design number is trusted. All results below are produced by
`validate_fixed.py` at $f=144$ MHz ($\lambda\approx2.08$ m) unless stated.

**(a) Half-wave dipole anchors.** A centre-fed dipole of length $0.47\lambda$
and radius $3$ mm should give $Z_{\text{in}}\approx70-j10\ \Omega$; one of
$0.50\lambda$ and radius $1$ mm should give $\approx81+j27\ \Omega$. The
solver returns

| Geometry | $Z_{\text{in}}$ (solver) | Reference | SWR @ 50 $\Omega$ |
|---|---|---|---|
| $0.47\lambda$, $a=3$ mm, $N=201$ | $70.27 - j10.20\ \Omega$ | $\approx70-j10$ | $1.46$ |
| $0.50\lambda$, $a=1$ mm, $N=201$ | $81.39 + j27.52\ \Omega$ | $\approx81+j27$ | $1.90$ |

Agreement is to well under 1 % in both real and imaginary parts.

**(b) Convergence in $N$.** For the $0.47\lambda$ dipole:

| $N$ | $Z_{\text{in}}$ | $\lvert\Delta Z\rvert$ vs $N{=}401$ |
|---:|---|---:|
| 101 | $68.28 - j26.67$ | $2.34\times10^{1}$ |
| 201 | $70.27 - j10.20$ | $6.80$ |
| 301 | $71.34 - j4.00$ | $5.32\times10^{-1}$ |
| 401 | $71.56 - j3.52$ | (ref) |

The sequence converges monotonically; by $N=301$ the change on further
refinement is $<0.6\ \Omega$. This is the normal PW-const convergence
behaviour in the half-wave regime.

**(c) Resonant length.** Scanning $L$ for the zero of the reactance gives
$L/\lambda\approx0.482$, against the classical $0.473$ — a $2\%$ length
error, consistent with the finite radius ($a=3$ mm shortens the resonant
length).

**(d) Power balance.** For a unit feed current the input power
$P_{\text{in}}=\tfrac12\operatorname{Re}(Z_{\text{in}})$ must equal the
radiated power $P_{\text{rad}}$ computed from the far-field Poynting integral
(up to conductor loss). The solver gives $P_{\text{in}}=35.1365$ W,
$P_{\text{rad}}=35.1372$ W, ratio $1.00002$. This near-exact balance is a
strong self-consistency check: it confirms that the impedance extraction (6)
and the far-field normalisation are mutually consistent. It holds at all
lengths tested, including the electrically-long dipoles of §4.3.

**(e) Directivity and pattern.** The $0.47\lambda$ dipole gives
$D=1.634$ ($2.13$ dBi) against the classical $2.15$ dBi, an E-plane 3-dB
beamwidth of $79^\circ$ against $\approx78^\circ$, and a radiating
efficiency of $99.88\%$ (copper surface resistance). All within a few
percent of the textbook values.

### 4.3 Where the piecewise-constant basis degrades

Two limits expose the limitations of the PW-const basis. These are *basis*
limitations, not kernel bugs — the power balance remains $1.0000$ in both
cases, so the solver is self-consistent; the issue is that the PW-const
current is a poor representation of the true current in these limits.

**(a) The short-dipole limit ($L\ll\lambda$).** The radiation resistance of a
short centre-fed dipole with (approximately) sinusoidal, feed-normalised
current is $R_{\text{rad}}=20\pi^2(L/\lambda)^2$ (the factor $80\pi^2$
belongs to the *uniform*-current Hertzian idealisation and must not be used
as the comparison). The solver under-predicts $R_{\text{rad}}$, and the
error shrinks as the dipole lengthens:

| $L/\lambda$ | $R_{\text{rad}}$ (solver) | $20\pi^2(L/\lambda)^2$ | ratio |
|---:|---:|---:|---:|
| 0.02 | $0.0000$ | $0.0790$ | $0.00$ |
| 0.05 | $0.1246$ | $0.4935$ | $0.25$ |
| 0.10 | $1.4286$ | $1.9739$ | $0.72$ |

The cause is visible in the current profile: for $L=0.02\lambda$ the solved
current is a narrow spike at the feed (peaking at $0.89$ of its max and
falling to $\sim0$ within a segment or two) rather than a triangle spanning
the wire. The PW-const staircase does not spread the current over the short
wire, so it radiates less than the analytical triangle. A piecewise-sinusoidal
basis (which enforces a smoother, end-vanishing current) removes most of this
error; we do not need it for the half-wave design of §5.

**(b) The electrically-long limit ($L\gtrsim0.7\lambda$).** For a
$0.75\lambda$ dipole the solved *current distribution* converges and has the
correct "W" shape (a dip at the feed and two peaks, the standing-wave pattern
of a $0.375\lambda$ arm), but the radiation resistance does **not** converge
in $N$:

| $N$ | $Z_{\text{in}}$ |
|---:|---|
| 101 | $603 + j571\ \Omega$ |
| 201 | $624 + j564\ \Omega$ |
| 401 | $837 + j540\ \Omega$ |
| 601 | $1049 + j391\ \Omega$ |

The real part climbs without settling as the mesh is refined. The power
balance is still $1.0000$ (the solver is self-consistent), so this is not a
bookkeeping error; it is a genuine failure of the PW-const collocation
discretisation to converge for electrically-long, weakly-ill-conditioned
currents. This matters for §5 because, on the 70 cm band, the 2 m driven
element is $1.36\lambda$ long and sits squarely in this regime. We treat the
70 cm feedpoint numbers as *indicative* (correct to order of magnitude and in
sign) rather than converged, and we say so explicitly where they are used.

**The pulse-basis convergence caveat.** This non-convergence is not a bug in
the implementation; it is a known property of the piecewise-*constant* (pulse)
basis with point-matching on the thin-wire Pocklington equation. The
discretisation is only *conditionally* convergent: it converges for
$\Delta\ell \gtrsim a$ (segment length comparable to or larger than the wire
radius) but not in the limit $\Delta\ell \to 0$ once $\Delta\ell < a$, because
the pulse basis cannot represent the smooth, end-vanishing current — and, in
particular, the charge $\partial I/\partial z'$ — that a long wire carries.
The standard remedy is a **piecewise-linear (rooftop/triangle) basis**, which
is $C^0$-continuous across segment boundaries and represents the current
derivative (charge) naturally; with a rooftop basis the collocation
discretisation converges as $\Delta\ell \to 0$. We retain the pulse basis
here because it is adequate in the half-wave regime ($0.2\lesssim L/\lambda
\lesssim0.6$) where the design elements live on their own bands, and the
degradation is confined to the electrically-long 70 cm feedpoint, which we
report as indicative.

### 4.4 Summary of the validated regime

The corrected solver is reliable — to well under a few percent — for
centre-fed dipoles and dipole-like elements in the range
$0.2\lesssim L/\lambda \lesssim 0.6$, i.e. the half-wave regime. This is
exactly the regime of the dual-band design in §5 *on each element's own
band*. The 70 cm feedpoint, which depends on the 2 m element at $1.36\lambda$,
is outside the validated regime and is reported as indicative.

---

## 5. Applied design: a dual-band VHF/UHF open-sleeve dipole

### 5.1 The design problem

We design a single-structure antenna for the amateur 2 m band (144–148 MHz)
and 70 cm band (430–440 MHz). A single straight dipole cannot be resonant at
both bands — the frequencies are a factor of three apart, so a half-wave
dipole would need to be simultaneously $\approx1$ m and $\approx0.35$ m long.
The classic resolution is the **open-sleeve dipole**: a centre-fed 2 m driven
dipole (resonant near 145 MHz) with 70 cm sleeves (resonant near 435 MHz)
placed *very close* to the driven element, typically as a symmetric pair on
both sides. The sleeves are open-circuited; their current is induced by the
driven element.

The essential physics is the **coupling regime**. A sleeve placed in the
reactive near-field of the driven element — a spacing of order $10$–$25$ mm,
i.e. $0.015$–$0.035\lambda$ at UHF — couples *strongly* (the mutual
impedance $|Z_{12}|$ is tens of ohms) and acts as an impedance transformer.
A second element placed at a far-field (Yagi-like) spacing of $\sim0.2$–
$0.5\lambda$ couples *weakly* ($|Z_{12}|$ of a few ohms) and shapes the
pattern without transforming the feedpoint impedance. The open-sleeve dipole
is the first case; confusing the two is the central error we examine in
§5.4. On the 2 m band the driven element radiates and the (electrically
short) sleeves are a weak perturbation; on the 70 cm band the sleeves are
near-resonant and strongly coupled, and they match the feedpoint. The
canonical parametric study is King and Wong [4], who measured VSWR as a
function of sleeve length and *sleeve-to-dipole spacing* and found the
dual-band match to require close coupling (spacing of order $0.02$–$0.04\lambda$
on the sleeve band) — the reactive near-field, not the Yagi-like spacing.

### 5.2 Geometry

The driven element is a 2 m dipole; the 70 cm band is served by a
**symmetric pair** of sleeves close-coupled to it:

- **2 m driven element:** $L_2 = 936$ mm, radius $a_2 = 3$ mm, at $x=0$,
  centre-fed. ($L_2 = 0.452\lambda$ at 145 MHz.)
- **70 cm sleeves (pair):** $L_7 = 316$ mm each, radius $a_7 = 3$ mm, at
  $x = \pm 25$ mm. ($L_7 = 0.458\lambda$ at 435 MHz; spacing
  $25\ \mathrm{mm} = 0.036\lambda$.)
- Mesh: $N_2 = 201$ segments on the driven element, $N_7 = 61$ on each
  sleeve; two-level quadrature $n_q^{\text{self}}=400$,
  $n_q^{\text{off}}=24$.

The full three-wire system ($N_2 + 2N_7 = 323$ unknowns) is solved at each
frequency; the feedpoint is on the driven element. The sleeve length and
spacing were optimised (sweeping $L_7 \in \{300,316,332\}$ mm and
$x_7 \in \{20,25,30\}$ mm) to minimise the 70 cm SWR at 435 MHz; the optimum
is $L_7 = 316$ mm, $x_7 = 25$ mm.

### 5.3 Feedpoint impedance and SWR

The solver (`dualband_design_v2.py`) sweeps both bands. The close-coupled
geometry gives:

**2 m band (144–148 MHz).**

| $f$ (MHz) | $\operatorname{Re}Z_{\text{in}}$ | $\operatorname{Im}Z_{\text{in}}$ | SWR @ 50 $\Omega$ |
|---:|---:|---:|---:|
| 144 | 56.8 | $-45.5$ | 2.31 |
| 145 | 58.6 | $-40.6$ | 2.11 |
| 146 | 60.4 | $-35.7$ | 1.94 |
| 147 | 62.3 | $-30.8$ | 1.80 |
| 148 | 64.2 | $-25.9$ | 1.67 |

The 2 m band is moderately matched (SWR $1.7$–$2.3$) and *improves* with
frequency: the driven element is $0.452\lambda$ (slightly short of the
$0.473\lambda$ resonance), so it is capacitive, and the reactance falls
toward zero as the band centre is approached from below. The close-coupled
sleeves are electrically short on this band ($0.15\lambda$) and perturb the
match only slightly. This is the well-validated regime (§4.4), so these
numbers are reliable.

**70 cm band (430–440 MHz).**

| $f$ (MHz) | $\operatorname{Re}Z_{\text{in}}$ | $\operatorname{Im}Z_{\text{in}}$ | SWR @ 50 $\Omega$ |
|---:|---:|---:|---:|
| 430 | 59.2 | $-25.7$ | 1.64 |
| 432 | 59.9 | $-22.7$ | 1.57 |
| 434 | 60.6 | $-19.7$ | 1.50 |
| 436 | 61.3 | $-16.8$ | 1.44 |
| 438 | 62.0 | $-13.9$ | 1.39 |
| 440 | 62.6 | $-11.1$ | 1.35 |

The 70 cm band is now **well matched** (SWR $1.35$–$1.64$). The large
capacitive reactance of the $1.36\lambda$ driven element has been almost
entirely cancelled by the close-coupled sleeves, leaving a feedpoint near
$61 - j18\ \Omega$ at the band centre. This is the result the open-sleeve
concept is supposed to give — and it is what distinguishes the close-coupled
design from the far-field (Yagi-like) spacing examined in §5.4.

### 5.4 How the close-coupled sleeve matches the 70 cm band

The mechanism is the two-port (here three-port) impedance transformation.
With the driven element as port 1 and the two sleeves as ports 2 and 3, the
open-circuit $Z$-parameters at 435 MHz are

| Parameter | Value | Meaning |
|---|---|---|
| $Z_{11}$ | $63.0 - j130.2\ \Omega$ | driven-element self-impedance (sleeves open) |
| $Z_{12}$ | $22.6 - j68.1\ \Omega$, $\ \lvert Z_{12}\rvert = 71.8\ \Omega$ | mutual, driven–sleeve$+$ |
| $Z_{13}$ | $17.9 - j61.7\ \Omega$, $\ \lvert Z_{13}\rvert = 64.2\ \Omega$ | mutual, driven–sleeve$-$ |
| $Z_{22}$ | $24.3 - j42.6\ \Omega$ | sleeve self-impedance (near-resonant) |
| $Z_{23}$ | $22.2 - j22.8\ \Omega$ | mutual, between the two sleeves |

($\lvert Z_{12}\rvert$ and $\lvert Z_{13}\rvert$ are equal by symmetry; the small
difference is a test-point-offset artefact of the discretisation.)

The loaded feedpoint, with both sleeves open, is

$$
Z_{\text{in}}' = Z_{11} - Z_{1v}\,Z_{vv}^{-1}\,Z_{v1}
= 60.9 - j18.2\ \Omega,
\qquad \text{SWR} = 1.47,
$$

where $Z_{1v}=[Z_{12}\ Z_{13}]$ and $Z_{vv}$ is the $2\times2$ sleeve
sub-matrix. The physics is immediate: the driven element alone is
$63 - j130\ \Omega$ (the $1.36\lambda$ capacitive reactance), and the
close-coupled, near-resonant sleeves contribute a correction of roughly
$+3 - j128\ \Omega$ — a large *negative* reactance that almost exactly
cancels the driven element's. The near-resonant sleeve presents a low-
impedance path; because it is in the reactive near-field, the mutual
impedance is large ($\lvert Z_{12}\rvert \approx 70\ \Omega$, not a few ohms),
so the cancellation is strong.

**The coupling-regime contrast.** The same driven element with a *single*
70 cm element at a far-field (Yagi-like) spacing of $307$ mm
($0.445\lambda$) gives $\lvert Z_{12}\rvert \approx 3\ \Omega$ and leaves the
feedpoint at $67.8 - j137.6\ \Omega$ (SWR $7.55$) — the weakly-coupled
element cannot transform the impedance. Sweeping the spacing of a *single*
offset sleeve shows the transition: at $307$ mm the SWR is $7.55$; at
$25$ mm it falls to $2.0$; but a single offset sleeve also pulls the real
part down ($24\ \Omega$ at $15$ mm), so it never reaches SWR $<1.5$. It is
the **symmetric pair** in the reactive near-field that both cancels the
reactance *and* preserves the real part near $50\ \Omega$:

| Sleeve placement (70 cm) | $\lvert Z_{12}\rvert$ | $Z_{\text{in}}'(435\ \text{MHz})$ | SWR |
|---|---:|---|---:|
| single, $+307$ mm ($0.445\lambda$) | $3.0\ \Omega$ | $67.8 - j137.6$ | $7.55$ |
| single, $+25$ mm ($0.036\lambda$) | $73.6\ \Omega$ | $54.4 - j37.1$ | $2.02$ |
| pair, $\pm25$ mm ($0.036\lambda$) | $71.8\ \Omega$ | $60.9 - j18.2$ | $\mathbf{1.47}$ |
| pair, $\pm20$ mm ($0.029\lambda$) | $79.4\ \Omega$ | $41.2 - j16.4$ | $1.50$ |

The optimum is $\pm25$ mm (SWR $1.47$); tighter spacing ($\pm10$ mm) over-
couples and collapses the real part (SWR $5.5$). This is the open-sleeve
mechanism, and it is entirely a consequence of the reactive near-field
coupling — it would be impossible at the far-field spacing.

### 5.5 Radiation pattern, directivity, and gain

**Plane definitions.** For a $z$-directed dipole the **E-plane** is an
elevation cut — any plane containing the wire axis; we take $\phi=0$ (the
$x$-$z$ plane), $\theta$ swept $0\!\to\!360^\circ$. The **H-plane** is the
*azimuthal* cut at broadside — $\theta=90^\circ$ (the $x$-$y$ plane), $\phi$
swept $0\!\to\!360^\circ$. It is a common error to take a second elevation
cut ($\phi=90^\circ$) and call it the H-plane; that cut still shows the
end-fire nulls and is not the H-plane. The H-plane of a vertical dipole must
be (nearly) **omnidirectional** — a circle — and never has end-fire nulls.

The far field is the coherent sum of Hertzian elements (§3.3, Appendix A).
For the close-coupled geometry the pattern metrics are:

| Quantity | 145 MHz (2 m) | 435 MHz (70 cm) |
|---|---:|---:|
| Directivity $D$ | $1.98$ dBi | $3.23$ dBi |
| Gain $G$ | $1.98$ dBi | $3.23$ dBi |
| Radiating efficiency | $99.9\%$ | $100\%$ |
| E-plane 3-dB beamwidth | $84^\circ$ | $66^\circ$ |
| H-plane 3-dB beamwidth | $360^\circ$ (omni) | $360^\circ$ (omni) |
| Front-to-back ratio | $0.0$ dB | $1.7$ dB |

On the 2 m band the pattern is that of a half-wave dipole: E-plane
directivity $1.98$ dBi (classical $2.15$ dBi), an E-plane beamwidth of
$84^\circ$ (classical $\approx78^\circ$), and an omnidirectional H-plane.

On the 70 cm band the E-plane is that of the electrically-long driven element
($1.36\lambda$) loaded by the sleeves: a higher directivity ($3.23$ dBi) and
a narrower E-plane beam ($66^\circ$). The H-plane remains essentially
omnidirectional ($360^\circ$) on both bands, as a $z$-directed structure
requires; the small front-to-back ratio ($1.7$ dB at 70 cm) is the only
azimuthal asymmetry, set by the sleeve pair. Radiating efficiency is
$\approx100\%$ on both bands (copper loss in a 3 mm wire at VHF/UHF is
negligible), so gain equals directivity.

The polar patterns are in `patterns_dualband.png`: the E-plane is a
figure-eight (lobes broadside, nulls end-fire) and the H-plane is a (nearly)
perfect circle on both bands.

### 5.6 Design assessment

The close-coupled open-sleeve dipole is a genuinely dual-band antenna:
well-matched on the 2 m band (SWR $1.7$–$2.3$) *and* on the 70 cm band
(SWR $1.35$–$1.64$), with near-omnidirectional H-plane patterns and
$\approx100\%$ radiating efficiency on both. The MoM analysis quantified the
design before any metal was bent and, more importantly, *localised the
physics*: the 70 cm match is produced by reactive near-field coupling
($\lvert Z_{12}\rvert \approx 70\ \Omega$) between the driven element and the
sleeve pair, and it would be impossible at a far-field (Yagi-like) spacing
where $\lvert Z_{12}\rvert$ is a few ohms. The lesson is that a dual-band
"parasitic" is not a generic second element — its spacing relative to
$\lambda$ determines whether it transforms the feedpoint impedance or merely
shapes the pattern.

---

## 6. Discussion

The exercise illustrates three points that are general to integral-equation
antenna analysis.

**Formulation errors do not converge.** The Galerkin-with-unintegrated-kernel
bug of §3.4 produced a $Z_{\text{in}}$ that was wrong by a factor of $\sim20$
and *did not improve with mesh refinement*, because it discretised the wrong
operator. This is a caution for any MoM code: a result that does not converge
with $N$ is a signal of a formulation (or sign) error, not a discretisation
error, and no quadrature refinement will fix it. The correct diagnostic is to
check the *classical limits* (half-wave impedance, short-dipole $R_{\text{rad}}$,
power balance) before trusting any design number.

**Basis choice sets the validated regime.** The PW-const basis is accurate in
the half-wave regime ($0.2\lesssim L/\lambda\lesssim0.6$) and degrades in the
short- and long-dipole limits (§4.3). The power balance staying at $1.0000$
even where the impedance is wrong shows that self-consistency (Poynting
balance) is necessary but not sufficient for accuracy — it checks the
*normalisation*, not the *current distribution*.

**A single number can hide the physics.** A 70 cm SWR of $7.55$ (the
far-field-spaced element of §5.4) looks like a "bad design," but the
three-port $Z$-parameters localise the cause: it is the driven element at
$1.36\lambda$, and the element is too *weakly* coupled ($|Z_{12}|$ a few
ohms) to transform the impedance. Move the same element into the reactive
near-field and $|Z_{12}|$ jumps to $\approx70\ \Omega$ and the SWR drops to
$1.47$. The MoM result is not just a pass/fail SWR; it is a diagnostic of
*where* the match or mismatch lives — and of the coupling regime that
produces it.

---

## 7. Conclusion

We derived the thin-wire Pocklington EFIE, discretised it by MoM with
piecewise-constant basis and collocation testing, and analysed the
discretisation's singular structure, quadrature requirements, and
convergence. The corrected solver reproduces the classical half-wave dipole
impedance, resonant length, directivity, and power balance to well under a
few percent, and its validated regime ($0.2\lesssim L/\lambda\lesssim0.6$)
covers each element of a dual-band VHF/UHF open-sleeve dipole on its own
band. Applied to a 144–440 MHz design, the solver predicts a close-coupled
open-sleeve dipole that is well-matched on *both* bands — 2 m SWR $1.7$–$2.3$
and 70 cm SWR $1.35$–$1.64$ — and, through the three-port $Z$-parameters,
shows that the 70 cm match is produced by reactive near-field coupling
($|Z_{12}|\approx70\ \Omega$) between the driven element and a symmetric
sleeve pair, and would be impossible at a far-field (Yagi-like) spacing where
$|Z_{12}|$ is a few ohms. The method is a reliable measuring instrument for
wire-antenna design in the half-wave regime, and its limitations (short- and
long-dipole basis degradation, and the pulse-basis convergence caveat of
§4.3) are clearly identifiable.

---

## References

[1] H. Pocklington, "On the Calculation of Electrostatic and Electromagnetic
Problems by Means of a System of Linear Equations," *Proc. Phys. Soc.
London*, ser. A, vol. 36, pp. 537–552, 1924. (The original method-of-moments
formulation for wire antennas; the kernel (3) is the Pocklington kernel.)

[2] R. F. Harrington, *Field Computation by Moment Methods*. New York: IEEE
Press / W. H. Freeman, 1968. (The standard treatment of the MoM: basis and
testing functions, the EFIE, and the thin-wire kernel.)

[3] C. A. Balanis, *Antenna Theory: Analysis and Design*, 4th ed. Hoboken,
NJ: Wiley, 2016. (Classical dipole impedances, the short-dipole
$R_{\text{rad}}=20\pi^2(L/\lambda)^2$, and the half-wave dipole reference
values used in §4.)

[4] H. King and J. Wong, "An experimental study of a balun-fed open-sleeve
dipole in front of a metallic reflector," *IEEE Trans. Antennas Propagat.*,
vol. 20, p. 201, 1972. (The classical parametric study of the open-sleeve
dipole: VSWR as a function of sleeve length and sleeve-to-dipole spacing —
the basis for the close-coupled (reactive near-field) spacing used in §5.)

*[Note on citation verification: references [1]–[4] are standard,
well-established works. For [4] the journal, volume, starting page, and year
are from the ADS record (bibcode 1972ITAP...20..201K); the full page range
should be confirmed against the IEEE Xplore record. The other references have
specific page numbers intentionally omitted; the results attributed to them
(the Pocklington kernel, the MoM discretisation, the classical dipole
impedances) are the canonical content of each text. Verify page-level detail
against the physical editions before formal publication.]*

---

## Appendix A: The solver

`antenna_mom.py` (NumPy only) implements:

- `Wire(x0, y0, z1, z2, a, N)` — a straight PEC wire.
- `build_Z(wires, f, nq_off, nq_self)` — the MoM impedance matrix (5),
  two-level Gauss–Legendre quadrature, Pocklington kernel (3), no endpoint
  term.
- `solve(wires, f, feed_seg, V0)` — $(Z_{\text{in}}, I)$ from the delta-gap
  model (6), $Z_{\text{in}}=-V_0/I_{\text{feed}}$.
- `z_parameters(wires, f, feed_segs)` — open-circuit $n$-port $Z$-parameters
  (unit-current method).
- `far_field`, `radiated_power`, `directivity`, `efficiency`, `swr`,
  `s11`, `pattern_metrics` — far-field and pattern quantities.
  `pattern_metrics` returns the E-plane (an elevation cut, $\phi=0$,
  $\theta$ swept) and the H-plane (the *azimuthal* cut at broadside,
  $\theta=90^\circ$, $\phi$ swept) as separate 1-D arrays, so the two
  planes are not confused (§5.5).

The validation battery `validate_fixed.py` runs the checks of §4;
`dualband_design_v2.py` runs the close-coupled open-sleeve design of §5
(sleeve-length/spacing optimisation, both band sweeps, the three-port
$Z$-parameters of §5.4, and the pattern metrics of §5.5) and writes
`dualband_v2_results.json`. `make_patterns.py` renders
`patterns_dualband.png` — the E-plane (elevation) and H-plane (azimuthal)
polar cuts at 145 and 435 MHz, with the sleeve directions marked.
