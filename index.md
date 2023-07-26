
**UNDER CONSTRUCTION**

<a name="inicio"></a>

This Fortran code provides a list of functions and subroutines that I typically use in my research, which can be classified in four categories:

- General propuse
- Statistics
- Linear algebra
- Optimization

All functions/subroutines assume double precision for real variables.

Efficiency of the algorithms is not guaranteed, and the code is distributed under the MIT license.

You can download the code from this [link](https://borjapetit.github.io/fortran_toolkit/src/toolkit.f90).

If you find any mistake, or have any suggestion, [contact me](mailto:bpetit@cunef.edu).

---

### General propuse

- [grid](#grid): generate a grid for a continuous variable.
- [interpolation](#interpolation.md): interpolate a value over a grid, returning position and distance.
- [interpolate](docs/interpolate.md): linearly interpolate a value over an n-dimensional grid (up to dimension 6).
- [timing](docs/timing.md): returns the number of seconds since 00:00h of the 1st day of the month (robust to parelalization).
- [multiplo](docs/multiplo.md): returns ```.TRUE.``` if an integer is a multiple of another user-provided integer.
- [iseven](docs/iseven.md): returns ```.TRUE.``` if a user-provided integer is even.
- [error](docs/error.md): print error message and interrupt execution.


<a name="grid"></a>

#### grid 

```fortran
function grid(maxv,minv,n,s) result(v)
  implicit none
  real(kind=8)            :: v(n)  ! output: generated grid
  real(kind=8)            :: maxv  ! input: upper-bound of x
  real(kind=8)            :: minv  ! input: lower-bound of x
  integer                 :: n     ! input: number of points in teh grid
  real(kind=8) , optional :: s     ! input: (optional) curvatura parameter 
```

This function creates a grid of $\texttt{n}$ points between $\texttt{maxv}$ and $\texttt{minv}$ with a curvature of $\texttt{s}$:

- if $\texttt{s}=1$ or missing: linear grid (default)
- if $\texttt{s}>1$: more grids points around $\texttt{maxv}$
- if $\texttt{s}<1$: more grids points around $\texttt{minv}$

**Dependencies**: [`error`](error.md)

**Example**

```fortran
! create a linear grid with 100 points between 0 and 10
vector = grid( 10.d0 , 0.0d0 , 100 )  

! create a quadratic grid with 500 points between -1 and 1
vector = grid( 1.d0 , -1.0d0 , 500 , 2.0d0 )  
```

[(back to index)](#inicio)

---

<a name="interpolation"></a>

#### interpolation 

```fortran
subroutine interpolation(pos,wth,xnow,xgrid)
  implicit none
  real(kind=8), intent(in)  :: xnow,xgrid(:)
  real(kind=8), intent(out) :: wth
  integer     , intent(out) :: pos
```

This subroutine finds the closest point in a given grid and return its position (```pos```) and the relative distance (```wth```) such that:

$$\texttt{xnow} = \texttt{vec}(\texttt{pos})\cdot\texttt{wth} \ + \ \texttt{vec}(\texttt{pos}-1)\cdot(\texttt{1.0d0}-\texttt{wth})$$

- If ```xnow``` < min(```xgrid```) the subroutine takes min(```xgrid```) as the interpolated value of ```xnow```.
- If ```xnow``` > max(```xgrid```) the subroutine takes max(```xgrid```) as the interpolated value of ```xnow```.

The vector ```xvec``` should be monotonic, so that 

$$\texttt{xgrid}(i) > \texttt{xgrid(i-1)} \ \ \ \text{or} \ \ \  \texttt{xgrid}(i) < \texttt{xgrid(i-1)}  \ \ \ \ \ \forall i\in[1,n]$$

This subroutine is mainly used by the function [```interpolate```](interpolate.md).

**Dependencies**: [```error```](error.md)

[(back to index)](../index.md)

**Example**

```fortran
! vector over which we wan to interpolate
vec = [1,2,3]

! interpolate x = 2.3 
call interpolation(pos,wth,2.3,vec)

! output: pos = 3, wth = 0.3
! 2.3 = vec(pos)*wth + vec(pos-1)*(one-wth) = 3*0.3 + 2*0.7

! interpolate x = 0
call interpolation(pos,wth,0,vec)

! since 0 is smaller than min(vec), the subroutine returns the min of vec
! output: pos = 2, wth = 0.0
! 1 = vec(pos)*wth + vec(pos-1)*(one-wth) = 2*0.0 + 1*1.0
```






---

### Statistics

- [```varmean```](docs/varmean.md): returns the average of a variable, allowing for weigths and mask.
- [```varvar```](docs/varvar.md): returns the variance of a variable, allowing for weigths and mask.
- [```varstd```](docs/varstd.md): returns the standard deviation of a variable, allowing for weigths and mask.
- [```correlation```](docs/correlation.md): returns the correlation of two variables, allowing for weigths and mask.
- [```percentile```](docs/percentile.md): returns the i-th percentile of a distribution, allowing for weigths and mask.
- [```olsreg```](docs/olsreg.md): returns the OLS coefficients of a linear regression (with up to 8 explanatory variables), allowing for weigths and mask.
- [```tauchen```](docs/tauchen.md): returns the transition matrix for a discretized AR(1) process.
- [```normaldist```](docs/normaldist.md): returns the distribution for a normal random variable given some mean and variance.
- [```randomnormal```](docs/randomnormal.md): returns a random draw (either a scalar or a vector) for a normal distribution.
- [```cdfn```](docs/cdfn.md): returns the cdf of a standard normal distribution.

---

### Linear algebra

- [```vect```](docs/vect.md): transform a matrix of into a vector (similar to ```reshape```).
- [```cumsum```](docs/cumsum.md): returns the vector with cummulative sum of a vector (as Matlab's ```cumsum``` function).
- [```diag```](docs/diag.md): returns a vector with the main diagonal of a matrix.
- [```inverse```](docs/inverse.md): returns the inverse of a sqaured matrix.

---

### Optimization

Algorithms for single-valued univariate equations:

- [```golden```](docs/golden.md): maximize a single-valued univariate equation using the Golden Search algorithm.
- [```brent```](docs/brent.md): find the root of a single-valued univariate equation using the Brent's method.

Algorithms for single-valued multivariate equations:

- [```simplex```](docs/simplex.md): minimize a single-valued multivariate equation using the Nelder-Mead algorithm.

Algorithms for systems of equations:

- **[```lmmin```](docs/lmmin.md): minimize a multivariate system of equations using the Levenberg–Marquardt algorithm.

  The Levenberg–Marquardt algorithm uses the Jacobian of the system to find the minimum. When evaluating the objective function is time-costly, computing the Jacobian may take too long. In those cases, one potential way of speeding up the algorithm is to update the Jacobian matrix using the Broyden's method, that does not require further function evaluations.

  - **[```broyden```](docs/broyden.md): updates a Jacobian matrix using the Broyden's method.

**Constrained optimization**: none of the algorithms in this toolkit is explicitly written to allow for constraints, but one can transform a contrained optimization problem into an uncontrained one using the [```normalize``` and ```denormalize```](docs/normalize.md) subroutines.

---

