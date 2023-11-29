
# Fortran toolkit

###### Borja Petit, 2023

---

**THIS SITE IS UNDER CONSTRUCTION**

<a name="inicio"></a>

This Fortran code provides a list of functions and subroutines that I typically use in my research, which can be classified in four categories:

- General purpose
- Statistics
- Linear algebra
- Optimization

All functions and subroutines assume double precision for real variables.

Efficiency of the algorithms is not guaranteed, and the code is distributed under the MIT license.

You can download the code from this [link](https://borjapetit.github.io/fortran_toolkit/toolkit.f90).

If you find any mistake, or have any suggestion, [contact me](mailto:bpetit@cunef.edu).

---

## General purpose

- [```grid```](docs/grid.md): generate a grid for a continuous variable.
- [```interpolation```](docs/interpolation.md): interpolate a value over a grid, returning position and distance.
- [```interpolate```](docs/interpolate.md): linearly interpolate a value over an n-dimensional grid (up to dimension 6).
- [```timing```](docs/timing.md): returns the number of seconds since 00:00h of the 1st day of the month (robust to parallelization).
- [```multiplo```](docs/multiplo.md): returns ```.TRUE.``` if an integer is a multiple of another user-provided integer.
- [```iseven```](docs/iseven.md): returns ```.TRUE.``` if a user-provided integer is even.
- [```error```](docs/error.md): print error message and interrupt execution.

## Statistics

- [```varmean```](docs/varmean.md): returns the average of a variable, and allows for weights and mask.
- [```varvar```](docs/varvar.md): returns the variance of a variable, and allows for weights and mask.
- [```varstd```](docs/varstd.md): returns the standard deviation of a variable, and allows for weights and mask.
- [```correlation```](docs/correlation.md): returns the correlation of two variables, and allows for weights and mask.
- [```percentile```](docs/percentile.md): returns the i-th percentile of a distribution, and allows for weights and mask.
- [```olsreg```](docs/olsreg.md): returns the OLS coefficients of a linear regression (with up to 8 explanatory variables), and allows for weights and mask.
- [```tauchen```](docs/tauchen.md): returns the transition matrix for a discretized AR(1) process.
- [```normaldist```](docs/normaldist.md): returns the distribution for a normal random variable given some mean and variance.
- [```randomnormal```](docs/randomnormal.md): returns a random draw (either a scalar or a vector) for a normal distribution.
- [```cdfn```](docs/cdfn.md): returns the cdf of a standard normal distribution.


## Linear algebra

- [```vect```](docs/vect.md): transform a matrix of into a vector (similar to ```reshape```).
- [```cumsum```](docs/cumsum.md): returns the vector with cumulative sum of a vector (as Matlab's ```cumsum``` function).
- [```diag```](docs/diag.md): returns a vector with the main diagonal of a matrix.
- [```inverse```](docs/inverse.md): returns the inverse of a squared matrix.


## Optimization

Algorithms for single-valued univariate equations:

- [```golden```](docs/golden.md): maximize a single-valued univariate equation using the Golden Search algorithm.
- [```brent```](docs/brent.md): find the root of a single-valued univariate equation using the Brent's method.

Algorithms for single-valued multivariate equations:

- [```simplex```](docs/simplex.md): minimize a single-valued multivariate equation using the Nelder-Mead algorithm.

Algorithms for systems of equations:

- [```lmmin```](docs/lmmin.md): minimize a multivariate system of equations using the Levenberg–Marquardt algorithm.

  The Levenberg–Marquardt algorithm uses the Jacobian of the system to find the minimum. When evaluating the objective function is time-costly, computing the Jacobian may take too long. In those cases, one potential way of speeding up the algorithm is to update the Jacobian matrix using the Broyden's method, that does not require further function evaluations.

  - [```broyden```](docs/broyden.md): updates a Jacobian matrix using the Broyden's method.


> **Constrained optimization**:<br>
None of the algorithms in this toolkit is explicitly written to allow for constraints, but one can transform a constrained optimization problem into an unconstrained one using the [```normalize``` and ```denormalize```](docs/normalize.md) subroutines.

<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>

