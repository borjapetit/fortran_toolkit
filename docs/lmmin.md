
<span style="text-align:right;display:block;">
<a href="https://borjapetit.github.io/fortran_toolkit/">Back to index</a>
</span>

## ```lmmin```

```fortran
subroutine lmmin(func,x,y,iy,exitcode,x0,itermax,damp,tol,shock,usebro,iprint)
  implicit none
  external                             :: func        ! user-supplied function to be minimize
  real(kind=8) , intent(out)           :: x(:)        ! values of "x" at minimum
  real(kind=8) , intent(out)           :: y(:)        ! valuf of "func" at "x"
  integer      , intent(out)           :: iy          ! number of function evaluations
  integer      , intent(out)           :: exitcode    ! indicator of convergence
  real(kind=8) , intent(in)            :: x0(:)       ! initial guess
  real(kind=8) , intent(in) , optional :: shock       ! shock to parameter values (as %)
  real(kind=8) , intent(in) , optional :: damp        ! damping factor
  real(kind=8) , intent(in) , optional :: tol         ! tolerance level
  integer      , intent(in) , optional :: itermax     ! max number of functione valuations
  integer      , intent(in) , optional :: iprint      ! indicator for printing behaviour
  integer      , intent(in) , optional :: usebro      ! indicator for the use of Broyden method to update Jacobian
```

This subroutine applies the Levenberg–Marquardt algorithm (click [here](https://en.wikipedia.org/wiki/Levenberg–Marquardt_algorithm) for more information) which minimizes the sum of squared errors of a (possibly nonlinear) system of multivariate equations.

This subroutine applies the Nelder-Mead algorithm (click [here](https://en.wikipedia.org/wiki/Nelder–Mead_method) for more information).

The user should input the function ```func``` and the initial guess ```x0```. The function ```func``` must be of the form:

```fortran
  function func(xvec) result(resid)
    real(kind=8)               :: xvec(:)
    real(kind=8) , allocatable :: resid(:)
    ...
  end function func
```

The subroutine returns the value(s) of ```x``` that makes ```func``` smaller than ```tol``` (close enoughs to zero), the number of function evaluations (```iy```), and an indicator, ```exitcode```:
- ```exitcode``` = 0: the algorithm found a root
- ```exitcode``` = 1: the jacobian of the system is close to singular
- ```exitcode``` = 2: the step in ```x``` is closed to xero
- ```exitcode``` = 9: maximum number of iterations reached

Optionally, the user can also supply:
- ```shock```  (double precision): shock to parameter values (as %), default = 0.05
- ```damp```  (double precision): damping factor, default = 1.0
- ```tol``` (double precision): tolerance level, default = 1.0d-8
- ```itermax``` (integer): max number of functione valuations, default = 500
- ```iprint``` (integer): indicator for printing behaviour
  - ```iprint``` = 0: don't print anything (default)
  - ```iprint``` = 1: print main results
  - ```iprint``` = 2: print main results and each iteration
- ```usebro``` (integer): indicator for the use of Broyden method to update Jacobian
  - ```usebro``` = 0: re-compute the Jacobian numerically (default)
  - ```usebro``` = 1: update the existing Jacobian using the Broyden method


**Note**: for constrained optimization problems, one can make use of the [```normalize and denormalize```](normalize.md) subroutines.

**Internal dependencies**: [```inverse```](inverse.md), [```broyden```](broyden.md)


