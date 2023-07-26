### lmmin

```fortran
subroutine lmmin_states_both(func,x,y,iy,ind,x0,itermax,damp,tol,toleach,shock,usebro,iprint)
  implicit none
  external                             :: func        ! user-supplied function to be minimize
  real(kind=8) , intent(out)           :: x(:)        ! values of "x" at minimum
  real(kind=8) , intent(out)           :: y(:)        ! valuf of "func" at "x"
  integer      , intent(out)           :: iy          ! number of function evaluations
  integer      , intent(out)           :: ind         ! indicator of convergence
  real(kind=8) , intent(in)            :: x0(:)       ! initial guess
  real(kind=8) , intent(in) , optional :: shock       ! shock to parameter values (as %)
  real(kind=8) , intent(in) , optional :: damp        ! damping factor
  real(kind=8) , intent(in) , optional :: tol         ! tolerance level
  real(kind=8) , intent(in) , optional :: toleach     ! tolerance level for each function
  integer      , intent(in) , optional :: itermax     ! max number of functione valuations
  integer      , intent(in) , optional :: iprint      ! indicator for printing behaviour
  integer      , intent(in) , optional :: usebro      ! indicator for the use of Broyden method to update Jacobian
```

This subroutine applies the Levenberg–Marquardt algorithm (click [here](https://en.wikipedia.org/wiki/Levenberg–Marquardt_algorithm) for more information) which minimizes the sum of squared errors of a (possibly nonlinear) system of multivariate equations.


This subroutine applies the Nelder-Mead algorithm (click [here](https://en.wikipedia.org/wiki/Nelder–Mead_method) for more information).


The user should input the function ```func``` and the initial guess ```x0```. The function ```func``` must be of the form:

```fortran
function func(x) result(f)
  implicit none
  real(kind=8) :: x(:),f
  f = ... some function of (x1,x2,...) ...
end function func
```

The subroutine returns the value(s) of ```x``` that makes ```func``` smaller than ```tol``` (close enought to zero), the number of function evaluations (```numiter```), and an indicator, ```exitcode```:

- ```exitcode``` = 0: the algorithm found a root
- ```exitcode``` = 1: the root is not within the interval (```x0```, ```x1```)
- ```exitcode``` = 9: maximum number of function evaluations reached

Optionally, the user can also supply a maximun number of fucntion evaluations (```itermax```), the level of tolerance (```tol```). Finally, the user can also control what the subroutine prints by setting the corresponding value of ```iprint```:

- ```iprint``` = 0: don't print anything (default)
- ```iprint``` = 1: print main results
- ```iprint``` = 2: print main results and each iteration

**Note**: for constrained optimization problems, one can make use of the [```normalize and denormalize```](normalize.md) subroutines.

**Dependencies**: [```inverse```](inverse.md), [```broyden```](broyden.md)

[(back to index)](../index.md)

---

**Example**

```fortran
program example

  use toolkit , only : dp,simplex

  implicit none

  real(dp) :: y
  real(dp) :: x0(2)
  real(dp) :: x1(2)
  integer  :: errorcode
  integer  :: itercount

  ! set initial guess to x1 = 2 and x2 = 3
  x0 = (/ 2.0d0 , 3.0d0 /)

  ! find the root of the Matyas function
  call simplex(matyas,x1,y,itercount,errorcode,x0)

  ! find the root of the Matyas function
  !     set a tolerance of 0.5
  !     set a max number of function evaluations to 1000
  !     print every iteration
  call simplex(matyas,x1,y,itercount,errorcode,x0,itermax=1000,tol=0.5d0,iprint=1)

  return

  contains

  ! Matyas function
  function matyas(x) result(y)
    implicit none
    real(dp) :: x(:),y
    y = dble(0.26)*( x(1)*x(1) + x(2)*x(2)) - dble(0.48)*x(1)*x(2)
    return
  end function matyas

end program example
```
