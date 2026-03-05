# Lightning Laplace Solver for Corner Domains via Rational Functions

A MATLAB implementation of the **"Lightning" solver** for Laplace's equation on 2D polygonal domains with corner singularities, based on the rational function approximation method with exponentially clustered poles.

## Overview

Solving Laplace's equation on domains with corners is notoriously difficult using classical polynomial or spectral methods because corner singularities cause slow convergence. The **Lightning algorithm** overcomes this by:

- Placing rational poles exponentially clustered near each corner of the domain.
- Solving the resulting boundary-value problem via weighted least squares.
- Achieving **root-exponential convergence** in the number of degrees of freedom (DoF), far outperforming polynomial-only approaches.
- Supporting both **Dirichlet** and **Neumann** boundary conditions.

The method is described in detail in:

> **Gopal, A. & Trefethen, L. N.** (2019).  
> *Solving Laplace problems with corner singularities via rational functions.*  
> SIAM Journal on Numerical Analysis, 57(5), 2074–2094.

---

## Repository Contents

| File | Description |
|------|-------------|
| `Lightning Laplace Solver for Corner Domains via Rational Functions.m` | Main MATLAB script with 5 test examples and the full solver implementation |
| `Lightning Laplace Solver for Corner Domains via Rational Functions.mlx` | MATLAB Live Script version with embedded output |
| `Lightning Laplace Solver for Corner Domains via Rational Functions.pdf` | PDF export of the Live Script with results and plots |

---

## Requirements

- **MATLAB** (R2019b or later recommended)
- *(Optional)* **[Chebfun](https://www.chebfun.org/)** — required only for AAA rational compression (`'aaa'` option)

---

## Quick Start

Open MATLAB and run the main script to execute all 5 built-in example problems:

```matlab
run('Lightning Laplace Solver for Corner Domains via Rational Functions.m')
```

Or call the `laplace` function directly for your own domain:

```matlab
% Define polygon vertices as a vector of complex numbers
P = [1, 1i, -1, -1i];              % diamond-shaped domain

% Define boundary data as a function handle
h = @(z) -log(abs(z));             % logarithmic boundary condition

% Solve and plot
laplace(P, h, 'tol', 1e-8, 'rel');
```

---

## Function Reference

### `laplace(P, g, options...)`

Solves the Laplace boundary value problem on a polygonal domain and returns the solution function.

#### Inputs

| Parameter | Type | Description |
|-----------|------|-------------|
| `P` | vector / cell array / string | Domain specification (see below) |
| `g` | function handle / vector / cell array | Boundary condition data |

#### Domain Specification (`P`)

- **Complex vector**: vertices of a polygon, e.g. `P = [1, 1i, -1, -1i]`
- **Integer `n`**: random polygon with `n` vertices
- **Negative integer `-n`**: random domain with `n` circular arc sides
- **Cell array**: mixed straight/circular arc sides; each cell is `{vertex}` (straight) or `{vertex, radius}` (arc)
- **Named strings**: `'L'` (L-shaped), `'pent'` (pentagon), `'snow'` (snowflake), `'iso'` (isosceles shape)

#### Options (name-value pairs / flags)

| Option | Type | Description |
|--------|------|-------------|
| `'tol', value` | numeric | Convergence tolerance (default: `1e-6`) |
| `'rel'` | flag | Use relative error weighting near corners (auto-enabled for discontinuous data) |
| `'steps'` | flag | Show step-by-step convergence plots |
| `'noplots'` | flag | Suppress all plots |
| `'slow'` | flag | Distribute poles more evenly (slower but useful for debugging) |
| `'noarnoldi'` | flag | Disable Arnoldi orthogonalization (educational/diagnostic use) |
| `'aaa'` | flag | Compress solution with AAA rational approximation (requires Chebfun) |

#### Outputs

| Output | Description |
|--------|-------------|
| `u` | Solution function handle: `u(z)` evaluates the real part (harmonic function) |
| `maxerr` | Estimated maximum error on the boundary |
| `f` | Complex analytic function handle (`u = real(f)`) |
| `Z` | Boundary sample points used in the least-squares solve |
| `Zplot` | Boundary points for plotting |
| `A` | Least-squares matrix |
| `v` | Coefficient vector |

---

## Built-in Examples

### Example 1 — Pinched Diamond
**Non-convex** polygon with symmetric re-entrant corners. Boundary data: `-log|z|` (logarithmic singularity).

```matlab
z = 0.3 * exp(0.25i * pi);
P = [1, z, 1i, 1i*z, -1, -z, -1i, -1i*z];
laplace(P, @(z) -log(abs(z)), 'tol', 1e-8, 'rel');
```

### Example 2 — Wedding Cake
**Non-convex** staircase shape with multiple 90° and 270° corners. Boundary data: distance from a point above the domain (`|z - 5i|`).

```matlab
P = [4, 4+1i, 3+1i, 3+2i, 2+2i, 2+3i, 1+3i, 1+4i, ...
    -1+4i, -1+3i, -2+3i, -2+2i, -3+2i, -3+1i, -4+1i, -4];
laplace(P, @(z) abs(z - 5i), 'tol', 1e-8, 'rel');
```

### Example 3 — Asymmetric Starburst
**Non-convex** star with acute asymmetric corners. Boundary data: high-frequency oscillatory function (`cos(3x)·sin(3y)`).

```matlab
angles = [0, 45, 80, 110, 160, 210, 260, 290, 330] * (pi/180);
radii  = [1.0, 0.2, 1.3, 0.4, 0.9, 0.1, 1.1, 0.5, 0.8];
P = radii .* exp(1i * angles);
laplace(P, @(z) cos(real(z)*3) .* sin(imag(z)*3), 'tol', 1e-8, 'rel');
```

### Example 4 — Pinched Bowtie
**Non-convex** bowtie with a very narrow choke point. Boundary data: parabolic bowl (`|z|²`). Tests robustness under severe ill-conditioning.

```matlab
P = [-2-2i, 2-2i, 0.2-0.1i, 2+2i, -2+2i, -0.2+0.1i];
laplace(P, @(z) abs(z).^2, 'tol', 1e-8, 'rel');
```

### Example 5 — Diamond
**Convex** diamond with acute salient corners. Serves as a clean baseline. Boundary data: `-log|z|`.

```matlab
P = [1, 0.5i, -1, -0.5i];
laplace(P, @(z) -log(abs(z)), 'tol', 1e-10, 'rel');
```

---

## Algorithm Details

The solver builds a basis of **rational functions** of the form:

```
f(z) = Σ aₖ (z - wc)^k  +  Σ bⱼ dⱼ / (z - pⱼ)
```

where:
- `wc` is the centroid of the domain for translation/scale invariance.
- The polynomial part uses an **Arnoldi orthogonalized** basis for numerical stability.
- Poles `pⱼ` are placed **outside** the domain along the outward bisector of each corner, with spacings growing like `exp(4(√k − √n))` (exponential clustering).
- The unknown coefficients are determined by a **weighted least-squares** fit to the boundary data.

The number of poles per corner is **adaptively increased** in each iteration: corners where the local error is large relative to the global error receive extra poles, until the target tolerance is met or a maximum of 100 poles per corner is reached.

---

## Output Plots

Each call to `laplace` produces a two-panel figure:

| Left panel | Right panel |
|------------|-------------|
| Convergence curve: error vs. √(DoF) on a log scale | Contour plot of the harmonic solution `u(z)` on the domain, with poles shown as red dots |

---

## License

This project is released for educational and research purposes. See the repository for full details.
