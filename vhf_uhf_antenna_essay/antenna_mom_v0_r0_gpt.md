Your traceback shows a **NumPy array broadcasting error in `antenna_mom.py`, inside `far_field()`**, when `radiated_power()` tries to calculate the antenna's radiated power.

However, there are **two separate problems** in your output:

1. **The immediate exception:** incompatible array dimensions in the far-field calculation.
2. **A more fundamental numerical problem:** your calculated antenna impedances are dramatically different from the classical dipole values, and a `ComplexWarning` indicates that complex-valued quantities are being converted to real numbers.

Fixing the broadcasting error will allow the validation script to proceed, but it will not fix the impedance calculations.

## 1. Why the traceback occurs

The relevant error is:

```text
ValueError: operands could not be broadcast together with shapes (10,) (361,)
```

The traceback points to this expression in `antenna_mom.py`, line 226:

```python
np.exp(1j * k * zq * ct) * ant.wts[i][n]
```

Your validation script calls:

```python
P = __import__("antenna_mom").radiated_power(ant, I)
```

at line 72. validate1

The call chain is:

```text
validate1.py
    └── radiated_power(ant, I)
          └── far_field(ant, I, th, ...)
                └── np.exp(1j * k * zq * ct) * ant.wts[i][n]
```

The problem is that the two arrays have different shapes:

- `zq`: apparently has shape `(10,)`, representing 10 quadrature points along a wire segment.
- `ct`: has shape `(361,)`, representing 361 observation angles, presumably spanning 0° to 180°.
- `np.exp(1j * k * zq * ct)` attempts to multiply these arrays element by element.

NumPy broadcasting aligns array dimensions from the right. Since neither dimension is 1 and `10 != 361`, the multiplication fails.

### How to fix the broadcasting error

You need to calculate the phase factor for **every combination of quadrature point and observation angle**.

Use two-dimensional broadcasting:

```python
np.exp(1j * k * zq[:, None] * ct[None, :])
```

Here:
- `zq[:, None]` has shape `(10, 1)`.
- `ct[None, :]` has shape `(1, 361)`.
- Their product has shape `(10, 361)`.

This produces a phase factor for every quadrature-point/observation-angle pair.

**Important:** this is the correction to the failing expression, not necessarily the complete correction to `far_field()`. The surrounding summation must also be consistent with these dimensions.

A likely correction inside `far_field()` is to change the phase-factor calculation from:

```python
np.exp(1j * k * zq * ct)
```

to:

```python
np.exp(1j * k * zq[:, None] * ct[None, :])
```

Then inspect how the result is combined with `ant.wts[i][n]` and the current coefficients. If the quadrature weights have shape `(10,)`, they will generally need to be reshaped to `(10, 1)` as well. The final summation over quadrature points should leave an array of shape `(361,)`, one far-field value per observation angle.

---

## 2. The `ComplexWarning` is another important problem

Your output contains:

```text
antenna_mom.py:137: ComplexWarning: Casting complex values to real discards the imaginary part
  A[oi + m, oj + n] = float(
```

This warning means that a complex-valued result is being converted to a real number, discarding its imaginary component.

For an electromagnetic Method of Moments (MoM) calculation, the impedance matrix generally needs to retain complex-valued entries because the electromagnetic kernel contains a complex exponential.

For example, a Green's function may contain:

\[
G(R)=\frac{e^{-jkR}}{4\pi R}
\]

Since

\[
e^{-jkR}=\cos(kR)-j\sin(kR),
\]

the kernel has both real and imaginary components.

If your matrix assembly discards the imaginary part, the resulting system no longer represents the intended complex-valued electromagnetic interaction.

### What to inspect

Look at line 137 in `antenna_mom.py`. If it resembles:

```python
A[oi + m, oj + n] = float(expression)
```

and `expression` is supposed to represent a complex impedance-matrix element, the conversion is incorrect.

A possible correction is:

```python
A[oi + m, oj + n] = complex(expression)
```

However, **do not apply this blindly**. The correct fix depends on whether the expression is intended to represent a complex matrix element or a mathematically real quantity. You should also verify that the matrix `A` was initialized with a complex dtype, for example:

```python
A = np.zeros((size, size), dtype=complex)
```

The key requirement is that complex-valued electromagnetic interactions must not lose their imaginary components during matrix assembly.

---

## 3. Your calculated dipole impedances are not physically credible

Consider your half-wave dipole at 144 MHz:

```text
lam/2 dipole @144 MHz (N=100):
Z = -0.00 + j319630.78
[classic 73.1 + j42.5]
```

The classical reference value is approximately

\[
Z_{\mathrm{in}}\approx 73.1+j42.5\ \Omega.
\]

Your computed result is approximately

\[
Z_{\mathrm{in}}\approx j319631\ \Omega.
\]

This is not a small numerical discrepancy. It indicates a substantial problem in the electromagnetic formulation, numerical integration, matrix assembly, or feed-point impedance calculation.

Your convergence results reinforce this concern:

| Segments \(N\) | Calculated reactance |
|---:|---:|
| 20 | \(+j24{,}132\ \Omega\) |
| 40 | \(+j8{,}235\ \Omega\) |
| 80 | \(+j81{,}786\ \Omega\) |
| 160 | \(-j122{,}174{,}285\ \Omega\) |

The impedance does not converge as the number of segments increases. Instead, it becomes wildly unstable.

Increasing the segment count is not going to resolve a fundamental error of this magnitude.

I would investigate these areas in order:

1. **Complex impedance-matrix assembly:** verify that the imaginary components are preserved.
2. **Kernel formulation:** verify the signs and derivatives of the Green's function and confirm the appropriate self-interaction treatment.
3. **Basis and testing functions:** verify the discretization and any assumptions about the wire current.
4. **Feed-point impedance:** check how the excitation is applied and how the input voltage and current are calculated.
5. **Quadrature:** verify that the integration weights, segment dimensions, and singular or near-singular interactions are handled correctly.

The `ComplexWarning` makes matrix assembly a particularly important place to start.

---

## 4. There is also a bug in the power-balance test

Your validation script contains:

```python
for frac, (Rc, Xc) in {
    0.25: (18.8, -108.8),
    0.5: (73.1, 42.5),
    0.75: (243.8, 104.5),
    1.0: (152.3, -46.1)
}.items():
    Z = dipole_Z(f, frac * lam(f), N=80, nq=8)
```

This loop overwrites `Z` on every iteration. When it finishes, `Z` corresponds to the **full-wavelength dipole**, not the half-wave dipole.

Later, the power-balance calculation uses:

```python
ant = ParallelWireAntenna(
    [{"x": 0.0, "a": a, "L": L, "N": 100, "feed": True}],
    f
)
I = ant.unit_current(0)
Pin = 0.5 * np.real(Z * np.abs(1.0) ** 2)
P = __import__("antenna_mom").radiated_power(ant, I)
```

Here, `ant` represents the half-wave dipole, but `Z` still holds the impedance calculated for the full-wave dipole.

That makes the power-balance calculation inconsistent.

The relevant code is in the power-balance and classical-impedance sections of your uploaded script. validate1

### How to fix it

Calculate the half-wave dipole impedance independently and use that value for the power-balance test.

For example:

```python
# --- 6. power balance ---
ant = ParallelWireAntenna(
    [{"x": 0.0, "a": a, "L": L, "N": 100, "feed": True}],
    f
)

I = ant.unit_current(0)

Z_half = dipole_Z(f, L, a=a, N=100)

Pin = 0.5 * np.real(Z_half)
P = __import__("antenna_mom").radiated_power(ant, I)

print(
    f"\npower balance: Pin = {Pin:.4f} W, "
    f"Prad = {P:.4f} W, ratio = {P / Pin:.4f}"
)
```

This assumes the excitation voltage is 1 V RMS, consistent with the convention \(P_{\mathrm{in}}=\tfrac12\operatorname{Re}(Z)|I|^2\) when the phasor current is a peak-amplitude quantity. Verify the convention used by `unit_current()` before relying on the numerical power balance.

This correction fixes the validation-script bookkeeping error, not the underlying impedance problem.

---

## 5. Recommended debugging sequence

I would address the problems in this order:

| Priority | Problem | Required action |
|---|---|---|
| 1 | Complex-to-real conversion | Inspect line 137 and preserve complex values wherever required. |
| 2 | Far-field broadcasting | Reshape quadrature-point and observation-angle arrays to compatible dimensions. |
| 3 | Impedance convergence | Verify the MoM formulation, matrix assembly, and feed-point calculation. |
| 4 | Power-balance bookkeeping | Calculate and use the impedance corresponding to the actual antenna. |
| 5 | Classical validation | Re-run the impedance and radiation tests after correcting the preceding issues. |

### What I need to identify the precise fix

You uploaded `validate1.py`, which establishes where the exception originates, but the failing implementation is in `antenna_mom.py`.

**Please upload `antenna_mom.py` as well.** The most useful sections are:

- Lines 110–155: impedance-matrix assembly and the `ComplexWarning`.
- Lines 190–245: `far_field()` and `radiated_power()`.
- The implementations of `unit_current()` and `z_matrix()`.

With that file, I can trace the array dimensions through the actual implementation, determine whether the proposed broadcasting correction is sufficient, and investigate why the calculated dipole impedance differs so dramatically from the classical result.