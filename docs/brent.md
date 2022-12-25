
### brent

```fortran
subroutine brent(func,x,numiter,exitcode,x0,x1,itermax,tol,iprint)
  external                             :: func      ! user-supplied function
  real(kind=8) , intent(out)           :: x         ! root
  integer      , intent(out)           :: numiter   ! number of functions evaluations
  integer      , intent(out)           :: exitcode  ! exit code
  real(kind=8) , intent(in)            :: x0        ! min value of x
  real(kind=8) , intent(in)            :: x1        ! max value of x
  real(kind=8) , intent(in) , optional :: tol       ! detault = 1.0d-8
  integer      , intent(in) , optional :: itermax   ! detault = 500
  integer      , intent(in) , optional :: iprint    ! detault = 0
```

Find the root of a single-valued univariate equation. The user must supply a function ```func```, and a maximum and minimum value of ```x```. The root of ```func``` must be between ```x0``` and ```x1```. The function must be of the form:

```fortran
function func(x) result(f)
  implicit none
  real(kind=8) :: x,f
  f = ... some function of x ...
end function func
```

Optionally, the user can also supply a maximun number of fucntion evaluations (```itermax```), the level of tolerance (```tol```). Finally, the user can also control what the subroutien prints by setting the corresponding value of ```iprint```:

- ```iprint``` = 0: don't print anything (default)
- ```iprint``` = 1: print main results
- ```iprint``` = 2: print main results and each iteration

The subroutine returns the value of ```x``` that makes ```func``` smaller than ```tol``` in absolute value (close enought to zero), the number of function evaluations (```numiter```), and an indicator, ```exitcode```:

- ```exitcode``` = 0: the algorithm found a root
- ```exitcode``` = 1: the root is not within the interval (```x0```, ```x1```)
- ```exitcode``` = 2: maximum number of function evaluations reached

**Dependencies**: none

[(back to index)](../inicio.md)

---

**Example**

```fortran
program example

  use toolkit , only : dp,brent

  implicit none

  real(dp) :: xvar,xmax,xmin
  integer  :: numiter,exitcode

  ! max and min values for x
  xmax = 15.0 ; xmin = 0.0

  ! find the root
  call brent(func,xvar,numiter,exitcode,xmin,xmax)

  ! find the root, printing every iteration
  call brent(func,xvar,numiter,exitcode,xmin,xmax,iprint=2)

  ! set a tolerance of 0.5 and a maximum of function evaluations of 1000
  call brent(func,xvar,numiter,exitcode,xmin,xmax,itermax=1000,tol=0.5d0,iprint=2)

  return

  contains

  function func(x) result(y)
    implicit none
    real(dp) :: x,y
    y = 2.0d0*(x**2.0) - 3.5d0*x + 5.0d0
    return
  end function func

end program example
```
