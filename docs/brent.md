## brent

Find the root of a single-valued univariate equation. The user must supply a function $f(x)$, and a maximum and minimum value of $x$, such that $x^*\in[x_0,x_1]$, with $f(x^*)=0$. 

**Dependencies**: none

[(back to index)](inicio.md)

---

### Sintax

```fortran
subroutine brent(func,x,numiter,exitcode,x0,x1,itermax,tol,iprint)
```

### Inputs

```fortran
external                :: func    ! user-supplied function
real(kind=8) ,          :: x0      ! min value of x
real(kind=8) ,          :: x1      ! max value of x
real(kind=8) , optional :: tol     ! detault = 1.0d-8
integer      , optional :: itermax ! detault = 500
integer      , optional :: iprint  ! detault = 0

! Set the value of iprint to control what's printed:
!   iprint = 0: don't print anything (default)
!   iprint = 1: print main results
!   iprint = 2: print main results and each iteration
```

### Outputs

```fortran
real(kind=8) :: x         ! root
integer      :: numiter   ! number of functions evaluations
integer      :: exitcode  ! exit code

! exitcode = 0: the algorithm found a root
! exitcode = 1: the root is not within the interval (x0,x1)
! exitcode = 2: maximum number of function evaluations reached
```

### Example

```fortran
program example

  use toolkit , only : dp,brent

  implicit none

  real(kind=8) :: xvar,xmax,xmin
  integer      :: numiter,exitcode

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
    real(kind=8) :: x(:),y
    y = 2.0d0*(x**2.0) - 3.5d0*x + 5.0d0
    return
  end func

end program example
```
