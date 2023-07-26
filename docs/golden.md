### golden

```fortran
subroutine golden(func,x,y,xmax,xmin,itermax,tol)
  implicit none
  external                             :: func
  real(kind=8) , intent(out)           :: x
  real(kind=8) , intent(out)           :: y
  real(kind=8) , intent(in)            :: xmax
  real(kind=8) , intent(in)            :: xmin
  real(kind=8) , intent(in) , optional :: tol
  integer      , intent(in) , optional :: itermax
```

This subroutine finds the maximum of a user supplied single-valued function, ```func```, with one unknown using the Golden Search algorithm, where

- ```x``` is the $\arg\max$ ```func```.
- ```y``` is ```func```(```x```)

The function ```func``` must be f the form:

```fortran
function func(x) result(f)
  implicit none
  real(kind=8) :: x,f
  f = ... some function of x ...
end function func
```

The user must also supply a range for the variable $x$, ```xmin``` and ```xmax```. Optionally, the user can also supply a maximun number of fucntion evaluations (```itermax```, 500 by default), the level of tolerance (```tol```, 1.0d-8 by default).


**Dependencies**: [```error```](error.md)

[(back to index)](../index.md)

---

**Example**

```fortran
program example

  use toolkit , only : dp,golden

  implicit none

  real(dp) :: xvar,f
  real(dp) :: xmax,xmin

  ! max and min values for x
  xmax = 12.0 ; xmin = 4.0

  ! maximize the function func
  call golden(func,xvar,f,xmax,xmin)

  ! maximize the function func set a tolerance of 0.5 and 
  ! a maximum of function evaluations of 1000
  call golden(func,xvar,f,xmax,xmin,itermax=1000,tol=0.5d0)

  return

  contains

  function func(x) result(y)
    implicit none
    real(dp) :: x,y
    y = x**2.0d0 - 3.5d0*x
    return
  end function func

end program example
```