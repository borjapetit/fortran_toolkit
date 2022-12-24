
# A Toolkit for Fortran90

## General propuse:

- [```grid```](grid.md): generate a grid for a continuous variable.
- [```interpolation```](interpolation.md): interpolate a value over a grid, returning position and distance.
- [```interpolate```](interpolate.md): linearly interpolate a value over an n-dimensional grid, with n $\leq$ 6.
- [```timing```](timing.md): returns the number of seconds since 00:00h of the 1st day of the month (_robust to parelalization_).
- [```multiplo```](multiplo.md): returns 1 if an integer is a multiple of another user-provided integer.
- [```iseven```](iseven.md): returns 1 if a user-provided integer is even.
- [```error```](error.md): print error message and interrupt execution.

## Statistics:

The toolkit contians functions to compute basic **summary statistics** (mean, variance, standard deviation and correlation). All of them allow for weigths. the user can compute them either unconditionally or conditionally, using a ```mask```.

- [```varmean```](varmean.md): returns the average of a variable, allowing for weigths and mask.
- [```varvar```](varvar.md): returns the variance of a variable, allowing for weigths and mask.
- [```varstd```](varstd.md): returns the standard deviation of a variable, allowing for weigths and mask.
- [```correlation```](correlation.md): returns the correlation of two variables, allowing for weigths and mask.
- [```percentile```](percentile.md): returns the i-th percentile of a variables, allowing for weigths and mask.

asdasdasd

- [```olsreg```](olsreg.md): returns the ols coefficients of a n-var regression (with n<=8), allowing for weigths and mask.

sdasdas

- [```tauchen```](tauchen.md): returns the transition matrix for a discretized ar(1) process.
- [```normaldist```](normaldist.md): returns the distribution for a nomral random variable.
- [```randomnormal```](randomnormal.md): returns a random draw for a nomal distribution
- [```cdfn```](cdfn.md): retutns the cdf of a nomabl distribution.

## Linear algebra

- [```vect```](vect.md): transform a matrix of NxM into a vector of NxM rows
- [```cumsum```](cumsum.md): returns the vector with cummulative sum of a vector (as Matlab's cumsum function)
- [```diag```](diag.md): returns the main diagonal of a matrix
- [```inverse```](inverse.md): returns the invesrse of a sqaured matrix

## Optimization

Algorithms for single-valued univariate equations:

- [```golden```](golden.md): maximize a single-valued univariate equation using the Golden search algorithm.
- [```brent```](brent.md): find the root a single-valued univariate equation using the Brent's method. The user must supply an interval containing the root.

Algorithms for single-valued multivariate equations:

- [```simplex```](simplex.md): minimize a single-valued multivariate equation using the Simplex algorithm.

Algorithms for systems of equations:

- [```lmmin```](lmmin.md): minimize a multivariate system of equations using the Levenberg–Marquardt algorithm.


The Levenberg–Marquardt algorithm uses the Jacobian of the system to find the minimum. When evaluating the objective function is time-costly, computing the Jacobian may take too long. In those cases, one potential way of speeding up the algorithm is to update the Jacobian matrix using the Broyden's method, that does not require further function evaluations.

- [```broyden```](broyden.md): updates a Jacobian matrix using the Broyden's method.


None of the algorithms in this toolkit is explicitly written to allow for inequality constraints, but one can transforme a contrained optimization problem into an uncontrained one using ```normalize``` and ```denormalize```.

- [```normalize```](normalize.md): transform a bounded variable into an unbounded one.
- [```denormalize```](denormalize.md): transform a undonded variable into an counded one.