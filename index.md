
**THIS SITE IS UNDER CONSTRUCTION**

<a name="inicio"></a>

This Fortran code provides a list of functions and subroutines that I typically use in my research. All functions and subroutines assume double precision for real variables.

Efficiency of the algorithms is not guaranteed, and the code is distributed under the MIT license.

If you find any mistake, or have any suggestion, [contact me](mailto:bpetit@cunef.edu).

---

In the preambule, several parameters used across functions and subroutines are defined:

```fortran
  integer  , parameter :: dp    = kind(1.0d00)             ! double-precision real variables
  real(dp) , parameter :: cero  = dble(0.00000000000000)   ! zero
  real(dp) , parameter :: medio = dble(0.50000000000000)   ! one half
  real(dp) , parameter :: uno   = dble(1.00000000000000)   ! one
  real(dp) , parameter :: dos   = dble(2.00000000000000)   ! two
  real(dp) , parameter :: cien  = dble(100.000000000000)   ! one hundred
  real(dp) , parameter :: mil   = dble(1000.00000000000)   ! one thousand
  real(dp) , parameter :: tolvl = dble(0.00000000010000)   ! tolerance level
```

---

## Econ-related functions and subroutines

- [```crra```](docs/crra.md): a constant relative risk aversion utilit function.
- [```ces```](docs/ces.md): a constant elasticity of substitution aggregator.

## General purpose

- [```grid```](docs/grid.md): generate a grid for a continuous variable.
- [```interpolation```](docs/interpolation.md): interpolate a value over a grid, returning position and distance.
- [```interpolate```](docs/interpolate.md): linearly interpolate a value over an n-dimensional grid (up to dimension 6).
- [```timing```](docs/timing.md): timing function robust to parallelization.
- [```multiplo```](docs/multiplo.md): check if an integer is a multiple of another user-provided integer.
- [```iseven```](docs/iseven.md): check if a user-provided integer is even.
- [```error```](docs/error.md): print error message and pause execution.
- [```num2text```](docs/num2text.md): converts a real number or an integer into a string.

## Statistics

- [```varmean```](docs/varmean.md): returns the average of a variable (allows for weights and mask).
- [```varvar```](docs/varvar.md): returns the variance of a variable (allows for weights and mask).
- [```varstd```](docs/varstd.md): returns the standard deviation of a variable (allows for weights and mask).
- [```correlation```](docs/correlation.md): returns the correlation between two variables (allows for weights and mask).
- [```percentile```](docs/percentile.md): returns the i-th percentile of a distribution (allows for weights and mask).
- [```olsreg```](docs/olsreg.md): returns the OLS coefficients of a linear regression with up to 8 explanatory variables (allows for weights and mask).
- [```tauchen```](docs/tauchen.md): returns the transition matrix for a discretized AR(1) process.
- [```normaldist```](docs/normaldist.md): returns the distribution for a normal random variable given some mean and variance.
- [```randomnormal```](docs/randomnormal.md): returns a random draw (either a scalar or a vector) for a normal distribution.
- [```cdfn```](docs/cdfn.md): returns the cdf of a standard normal distribution.
- [```shuffle_vect```](docs/shuffle_vect.md): this subroutine returns a vector with the shuffled values of a used-provided vector using the Fisher-Yates algorithm.
- [```fisheryates```](docs/fisheryates.md): this subroutine returns a vector filled with integers from 1 to $n$ and then shuffled using the Fisher-Yates algorithm.


## Linear algebra

- [```vect```](docs/vect.md): transform a matrix of into a vector (similar to ```reshape```).
- [```cumsum```](docs/cumsum.md): returns the vector with cumulative sum of a vector (as Matlab's ```cumsum``` function).
- [```diag```](docs/diag.md): returns a vector with the main diagonal of a matrix.
- [```inverse```](docs/inverse.md): returns the inverse of a squared matrix.
- [```transmat```](docs/transmat.md): returns the transpose of a matrix.

## Optimization

*Algorithms for single-valued univariate equations:*

- [```golden```](docs/golden.md): maximize a single-valued univariate equation using the Golden Search algorithm.
- [```brent```](docs/brent.md): find the root of a single-valued univariate equation using the Brent's method.

*Algorithm for single-valued multivariate equations:*

- [```simplex```](docs/simplex.md): minimize a single-valued multivariate equation using the Nelder-Mead algorithm.

*Algorithm for systems of equations:*

- [```lmmin```](docs/lmmin.md): minimize a multivariate system of equations using the Levenberg–Marquardt algorithm.

  The Levenberg–Marquardt algorithm uses the Jacobian of the system to find the minimum. When evaluating the objective function is time-costly, computing the Jacobian may take too long. In those cases, one potential way of speeding up the algorithm is to update the Jacobian matrix using the Broyden's method, that does not require further function evaluations.

*Other functions:*

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

