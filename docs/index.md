
# A Toolkit for Fortran90

## General propuse:

- [```grid```](grid.md): generate a grid for a continuous variable.
- [```interpolation```](interpolation.md): interpolate a value over a grid, returning position and distance.
- [```interpolate```](interpolate.md): linearly interpolate a value over an n-dimensional grid, with n $\leq$ 6.
- [```timing```](timing.md): returns the number of seconds since 00:00h of the 1st day of the month (_robust to parelalization_).
- [```multiplo```](#multiplo): returns 1 if an integer is a multiple of another user-provided integer.
- [```iseven```](#iseven): returns 1 if a user-provided integer is even.
- [```error```](#error): print error message and interrupt execution.

## Statistics:

The toolkit contians functions to compute basic **summary statistics** (mean, variance, standard deviation and correlation). All of them allow for weigths. the user can compute them either unconditionally or conditionally, using a ```mask```.

- [```varmean```](#varmean): returns the average of a variable, allowing for weigths and mask.
- [```varvar```](#varvar): returns the variance of a variable, allowing for weigths and mask.
- [```varstd```](#varstd): returns the standard deviation of a variable, allowing for weigths and mask.
- [```correlation```](#correlation): returns the correlation of two variables, allowing for weigths and mask.
- [```percentile```](#percentile): returns the i-th percentile of a variables, allowing for weigths and mask.

asdasdasd

- [```olsreg```](#olsreg): returns the ols coefficients of a n-var regression (with n<=8), allowing for weigths and mask.

sdasdas

- [```tauchen```](#tauchen): returns the transition matrix for a discretized ar(1) process.
- [```normaldist```](#normaldist): returns the distribution for a nomral random variable.
- [```randomnormal```](#randomnormal): returns a random draw for a nomal distribution
- [```cdfn```](#cdfn): retutns the cdf of a nomabl distribution.

## Linear algebra:

- [```vect```](#vect): transform a matrix of NxM into a vector of NxM rows
- [```cumsum```](#cumsum): returns the vector with cummulative sum of a vector (as Matlab's cumsum function)
- [```diag```](#diag): returns the main diagonal of a matrix
- [```transmat```](#transmat): returns the transpose of a square matrix
- [```inverse```](#inverse): returns the invesrse of a sqaured matrix

## Optimization

This toolkit contains several optimization algorithms. None of these are explicitly writen to allow for inequality constraints, but one can transforme a contrained optimization problem into an uncontrained one using ```normalize``` and ```denormalize```.

- [```simplex```](#simplex): Simplex algorithm
- [```lmmin```](#lmmin): Levenberg–Marquardt algorithm
- [```golden```](#golden): Golden search algorithm
- [```brent```](#brent): Brent method
- [```normalize```](#normalize): transform a bounded variable into an unbounded one [for optimizaton]
- [```denormalize```](#denormalize): transform a undonded variable into an counded one [for optimizaton]
- [```broyden```](#broyden): updates a Jacobian matrix using the Broyden's method.
