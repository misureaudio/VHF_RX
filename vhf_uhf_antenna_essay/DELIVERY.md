# Deliverable summary — VHF/UHF antenna essay

**Essay:** `essay_vhf_uhf_antennas.md` (5,876 words) — a rigorous, expert-level
academic essay on VHF/UHF wire antennas, from the thin-wire EFIE through the
Method-of-Moments discretisation to a worked dual-band design.

**What was fixed in the MoM code** (`antenna_mom.py`, the pre-fix version
preserved as `antenna_mom_v0_buggy.py`):
- The original had a **Galerkin-with-unintegrated-kernel** defect: it fed the
  un-integrated Pocklington kernel into a double (Galerkin) integral, which
  discretises the wrong operator. A half-wave dipole returned
  `0 − j1500 Ω` instead of `~73 + j42 Ω`, and the error did not converge with
  mesh refinement (a formulation error, not a discretisation error).
- Corrected to **PW-const basis + point-matching (collocation) + Pocklington
  kernel (G and d²G/dz² combined before evaluation) + two-level
  Gauss–Legendre quadrature + `Zin = −V0/I_feed`**.
- Also fixed the **two-port Z-parameter** routine: the port voltage is the
  *negative* of the EFIE residual, and a true open-circuit column must enforce
  `I_feed = 0` on every other port.
- Fixed the **pattern planes** (§5.5): the H-plane is now the *azimuthal* cut
  at broadside (θ=90°, φ swept) — a (nearly) omnidirectional circle for a
  z-directed dipole — not a second elevation cut (φ=90°) mislabelled as the
  H-plane. The F/B ratio is forward (θ=90, φ=0) vs backward (θ=90, φ=180),
  not the whole-plane minimum.

**Validation (all reproduced from the code):**
- 0.47λ / 3 mm dipole → `70.27 − j10.20 Ω` (SWR 1.46); 0.50λ / 1 mm →
  `81.39 + j27.52 Ω` (SWR 1.90). Both <1 % from the classical values.
- Power balance `Prad/Pin = 1.00002` (Poynting self-consistency).
- Documented basis limits: PW-const under-spreads the current on very short
  dipoles and is only *conditionally* convergent for electrically-long
  (≥0.75λ) dipoles (pulse basis with point-matching); the standard remedy is
  a piecewise-linear (rooftop) basis (§4.3b).

**Dual-band design (144–440 MHz close-coupled open-sleeve dipole):**
- A centre-fed 2 m driven element with a **symmetric pair** of 70 cm sleeves
  close-coupled in the reactive near-field (x = ±25 mm ≈ 0.036λ at UHF),
  optimised to minimise the 70 cm SWR.
- 2 m band (144–148 MHz): SWR 1.67–2.31 — well matched.
- 70 cm band (430–440 MHz): SWR 1.35–1.64 — **well matched** (the sleeves
  cancel the driven element's −j130 Ω reactance).
- Three-port analysis at 435 MHz: `Z11 = 63.0 − j130.2 Ω`,
  `Z12 = 22.6 − j68.1 Ω` (|Z12| = 71.8 Ω), loaded `Zin = 60.9 − j18.2 Ω`
  (SWR 1.47). The match is produced by reactive near-field coupling.
- Coupling-regime contrast: at a far-field (Yagi-like) spacing of 307 mm
  (0.445λ) |Z12| ≈ 2 Ω and the 70 cm band is left mismatched (SWR 7.55);
  the *coupling regime*, not the mere presence of a second element, sets the
  dual-band match (§5.4).
- Pattern: 2 m directivity 1.98 dBi (half-wave-like); 70 cm directivity
  3.23 dBi (electrically-long driven element), efficiency ≈100 % on both
  bands; the H-plane is omnidirectional (3-dB beamwidth 360°) on both bands.
