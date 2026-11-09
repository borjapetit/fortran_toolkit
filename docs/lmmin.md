
<span style="text-align:right;display:block;">
<a href="https://borjapetit.github.io/fortran_toolkit/">Back to index</a>
</span>

## ```lmmin```

```fortran
subroutine lmmin(func,x,y,iy,exitcode,x0,itermax,damp,tol,toleach,shock,usebro,iprint,tolrel,maxstall)
  implicit none
  external                             :: func        ! user-supplied function to be minimize
  real(kind=8) , intent(out)           :: x(:)        ! values of "x" at minimum
  real(kind=8) , intent(out)           :: y(:)        ! valuf of "func" at "x"
  integer      , intent(out)           :: iy          ! number of function evaluations
  integer      , intent(out)           :: exitcode    ! indicator of convergence
  real(kind=8) , intent(in)            :: x0(:)       ! initial guess
  real(kind=8) , intent(in) , optional :: shock       ! shock to parameter values (as %)
  real(kind=8) , intent(in) , optional :: damp        ! initial damping factor
  real(kind=8) , intent(in) , optional :: tol         ! tolerance level for the sum of squared errors
  real(kind=8) , intent(in) , optional :: toleach     ! tolerance level for each equation
  real(kind=8) , intent(in) , optional :: tolrel      ! relative improvement below which an iteration is unproductive
  integer      , intent(in) , optional :: itermax     ! max number of function evaluations
  integer      , intent(in) , optional :: iprint      ! indicator for printing behaviour
  integer      , intent(in) , optional :: usebro      ! indicator for the use of Broyden method to update Jacobian
  integer      , intent(in) , optional :: maxstall    ! max number of consecutive unproductive iterations
```

This subroutine applies the Levenberg–Marquardt algorithm (click [here](https://en.wikipedia.org/wiki/Levenberg–Marquardt_algorithm) for more information) which minimizes the sum of squared errors of a (possibly nonlinear) system of multivariate equations.

Given the current point, with residuals ```y``` and (numerical) Jacobian ```J```, the step ```dx``` solves ```(J'J + lambda D) dx = -J'y```, where ```D``` is a diagonal scaling matrix (the largest diagonal of ```J'J``` observed so far for each unknown, with a positive floor; Moré, 1978) and ```lambda``` is the damping factor. The linear system is solved by Cholesky decomposition. The damping factor is updated with the gain ratio (actual over predicted reduction in the error), as in Nielsen (1999).

The user should input the function ```func``` and the initial guess ```x0```. The function ```func``` must be of the form:

```fortran
  function func(xvec) result(resid)
    real(kind=8)               :: xvec(:)
    real(kind=8) , allocatable :: resid(:)
    ...
  end function func
```

The size of ```x``` must be equal to the size of ```x0``` (number of unknowns), and the size of ```y``` must be equal to the size of the vector returned by ```func``` (number of equations); the latter is not checked by the subroutine.

The subroutine returns the value(s) of ```x``` with the lowest sum of squared errors found, the residuals at that point (```y```), the number of function evaluations (```iy```), and an indicator, ```exitcode```:
- ```exitcode``` = 0: the algorithm found a root (sum of squared errors below ```tol```, or every residual below ```toleach``` in absolute value)
- ```exitcode``` = 1: the gradient ```J'y``` is close to zero (or the linear system cannot be solved)
- ```exitcode``` = 2: the step in ```x``` is close to zero
- ```exitcode``` = 3: no progress: ```maxstall``` consecutive unproductive iterations
- ```exitcode``` = 9: maximum number of function evaluations reached

Optionally, the user can also supply:
- ```shock```  (double precision): relative shock to parameter values to compute the numerical Jacobian, default = 0.05. If the residuals do not change, the shock is doubled (at most 5 times); if they still do not change, the corresponding column of the Jacobian is set to zero.
- ```damp```  (double precision): initial damping factor, default = 1.0
- ```tol``` (double precision): tolerance level for the sum of squared errors, default = 1.0d-8
- ```toleach``` (double precision): tolerance level for each equation, default = 1.0d-8
- ```tolrel``` (double precision): relative improvement in the error below which an iteration is considered unproductive, default = 1.0d-6
- ```maxstall``` (integer): number of consecutive unproductive iterations after which the algorithm stops, default = 10
- ```itermax``` (integer): max number of function evaluations, default = 500
- ```iprint``` (integer): indicator for printing behaviour
  - ```iprint``` = 0: don't print anything (default)
  - ```iprint``` = 1: print main results
  - ```iprint``` = 2: print main results and each iteration
- ```usebro``` (integer): indicator for the use of Broyden method to update Jacobian
  - ```usebro``` = 0: re-compute the Jacobian numerically after every successful step (default)
  - ```usebro``` = 1: update the existing Jacobian using the Broyden method; a new numerical Jacobian is computed when a Broyden-updated Jacobian fails to produce a better point

**Note**: for constrained optimization problems, one can make use of the [```normalize and denormalize```](normalize.md) subroutines.

**Internal dependencies**: none (it only uses the kind ```dp``` and the constants of ```toolkit_pars```)

**Updates**: November 2026. (1) The shock used to compute a column of the Jacobian is enlarged at most 5 times when the residuals do not change (it was enlarged without limit, until ```itermax```, and the test was on the sum of squared errors instead of the residuals); (2) scaling with a positive floor and Cholesky solve, instead of the inverse of ```(1+lambda) diag(J'J)```, which was singular with a zero column; (3) damping updated with the gain ratio; (4) stop after ```maxstall``` unproductive iterations (```exitcode``` = 3) and relative test on the step; (5) points evaluated to compute the Jacobian are used when they improve on the current point, and the best point evaluated is returned.
