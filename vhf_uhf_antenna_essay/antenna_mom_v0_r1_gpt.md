I've inspected `antenna_mom.py`. There is an important discrepancy between the file you uploaded and the code shown in your original traceback: **the uploaded file already contains a broadcasting correction, and its matrix assembly does not contain the `float(...)` conversion mentioned in the warning.** antenna_mom

This suggests that the file you uploaded may differ from the version Python executed, or that the source was modified after the failed run.

There are also issues in the electromagnetic formulation worth investigating. I'll separate the immediate runtime problem from the numerical problems so you can fix them in the right order.

## 1. The broadcasting error: your uploaded code already addresses it

Your original traceback showed:

```python
np.exp(1j * k * zq * ct) * ant.wts[i][n]
```

The uploaded `antenna_mom.py`, however, contains:

```python
zq = ant.zpts[i][n]
base = np.sin(k * (zq - e[n])) * ant.wts[i][n]
M += In * (base @ np.exp(1j * k * zq[None, :] * ct[:, None]))
```

This version is in `far_field()`. antenna_mom

The expression

```python
zq[None, :] * ct[:, None]
```

correctly produces a matrix with dimensions `(361, 10)` when there are 361 observation angles and 10 quadrature points.

The subsequent matrix-vector product

```python
base @ np.exp(1j * k * zq[None, :] * ct[:, None])
```

then produces one field value per observation angle.

**Therefore, don't replace this expression with another broadcasting fix.** It is already correct in the uploaded file.

### Check which file Python is actually importing

From your project directory, run:

```powershell
py -c "import antenna_mom; print(antenna_mom.__file__)"
```

Verify that the path points to the `antenna_mom.py` you uploaded.

Then inspect the actual source Python is importing:

```powershell
py -c "import antenna_mom, inspect; print(inspect.getsource(antenna_mom.far_field))"
```

If the output shows the old expression without `None` indexing, you have identified the immediate cause: Python is executing an older version of the module.

If it shows the uploaded implementation, rerun the validation script and examine the new traceback, if any.

## 2. The complex-value warning does not match this file

Your original output reported:

```text
antenna_mom.py:137: ComplexWarning: Casting complex values to real discards the imaginary part
A[oi + m, oj + n] = float(...)
```

But in the uploaded file, the corresponding assignment is:

```python
A[oi + m, oj + n] = np.sum(
    K * wm[:, None]
    * (self.wts[j][n] * self.fdiag[j][n])[None, :]
)
```

The matrix is initialized as:

```python
A = np.zeros((self.Ntot, self.Ntot), dtype=complex)
```

Both choices are appropriate for a complex-valued MoM impedance matrix. antenna_mom

Consequently, **the uploaded version does not exhibit the specific complex-to-real conversion reported in your traceback**.

Run the import-path checks above before making any changes to this section. If the warning persists with the uploaded implementation, we should investigate the actual file being executed and reproduce the warning against that version.

## 3. A normalization inconsistency exists in the electromagnetic kernel

This is a separate issue I found by examining the uploaded source.

Your documentation defines the Green's function as

\[
G(R)=\frac{e^{-jkR}}{4\pi R}.
\]

It then gives the kernel prefactor as

\[
\frac{j\omega\mu_0}{4\pi}
\left(G+\frac{1}{k^2}\frac{d^2G}{dz^2}\right).
\]

The implementation follows that prefactor:

```python
return (1j * k * eta0 / (4.0 * np.pi)) * (
    G + d2G / k ** 2
)
```

However, because `G` already includes the factor \(1/(4\pi)\), this introduces an additional \(1/(4\pi)\) compared with the standard formulation using the Green's function as defined in your file.

To see the inconsistency, start from the vector potential:

\[
A_z(z)=\mu_0\int I(z')G(z,z')\,dz',
\]

where \(G=e^{-jkR}/(4\pi R)\).

The resulting tangential electric field gives an EFIE kernel with prefactor

\[
j\omega\mu_0
\left(G+\frac{1}{k^2}\frac{d^2G}{dz^2}\right),
\]

up to the overall sign determined by the chosen field and excitation conventions.

Thus, if you retain the existing definition of \(G\), the corresponding kernel prefactor should be

```python
1j * k * eta0
```

rather than

```python
1j * k * eta0 / (4.0 * np.pi)
```

The relevant implementation is in `kernel()`. antenna_mom

**This is a genuine formulation inconsistency, but it is not sufficient to explain your entire impedance error.** Correcting this normalization alone would rescale the impedance; it would not explain the extreme instability as the segment count increases.

I would correct it only after confirming that the documented Green's function and field convention are the ones you intend to implement.

## 4. Why the impedance results still need investigation

The uploaded code has a more fundamental numerical challenge: the matrix assembly evaluates the second derivative of the Green's function for self-interactions, using the wire radius to regularize the distance.

The second derivative contains terms proportional to

\[
\frac{1}{R^5},\qquad \frac{1}{R^4},\qquad \frac{1}{R^3}.
\]

These terms can become very large when the observation and source points are close together. Your implementation evaluates them through quadrature on each segment. antenna_mom

For a thin-wire Pocklington solver, this makes the treatment of self-interaction and the discretization especially important. A correct closed-form derivative does not, by itself, guarantee accurate numerical integration of the resulting kernel.

Your reported impedance sequence,

\[
X_{20}\approx 24\,132\ \Omega,\qquad
X_{80}\approx 81\,786\ \Omega,\qquad
X_{160}\approx-122\,174\,285\ \Omega,
\]

is a strong indication that the assembled system is not behaving numerically as intended.

The next checks should be:

1. Correct the kernel normalization consistently.
2. Check the conditioning of the assembled matrix as \(N\) increases.
3. Inspect the self-segment matrix elements separately from mutual-segment elements.
4. Verify the feed-point constraint and the extraction of input impedance.
5. Compare the numerical current distribution with the expected approximately sinusoidal distribution of a center-fed half-wave dipole.

## 5. What I recommend doing now

Follow this sequence:

**First:** confirm the actual imported file using the two PowerShell commands above.

**Second:** rerun `py -m validate1` and check whether the original broadcasting error and `ComplexWarning` remain.

**Third:** correct the kernel normalization and investigate self-interaction integration before trusting the calculated impedances.

**Fourth:** rerun the convergence test. Do not use the classical impedance table as a pass/fail benchmark until the matrix formulation and feed-point extraction have been verified.

One additional issue in `validate1.py` remains: its power-balance calculation uses the impedance left over from the preceding length scan rather than the half-wave impedance. That should also be corrected, as discussed earlier.

**The bottom line:** your uploaded `antenna_mom.py` already fixes the broadcasting expression shown in the original traceback, but the normalization inconsistency and the unstable impedance convergence need separate attention. I would not regard the solver as validated until those numerical issues have been resolved.