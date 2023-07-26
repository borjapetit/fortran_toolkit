## golden

```fortran
subroutine golden(func,x,y,numiter,xmax,xmin,itermax,tol,iprint)
  implicit none
  external                             :: func      ! user-supplied function
  real(kind=8) , intent(out)           :: x         ! output: arg max of func
  real(kind=8) , intent(out)           :: y         ! output: max of func
  integer      , intent(out)           :: numiter   ! output: number of function evaluations
  real(kind=8) , intent(in)            :: xmax      ! input: upper-bound of x
  real(kind=8) , intent(in)            :: xmin      ! input: lower-bound of x
  real(kind=8) , intent(in) , optional :: tol       ! input: (optional) level of tolerance [detault = 1.0d-8]
  integer      , intent(in) , optional :: itermax   ! input: (optional) maximum function evaluations [default = 500]
  integer      , intent(in) , optional :: iprint    ! input: (optional) control what's printed [detault = 0]
```

This subroutine finds the maximum of a user-supplied single-valued function, $\texttt{func}$, with one unknown using the Golden Search algorithm. The function $\texttt{func}$ must be f the form:

```fortran
function func(x) result(f)
  implicit none
  real(kind=8) :: x,f
  f = ... some function of x ...
end function func
```

The user must also supply an upper- and lower-bound range ($\texttt{xmin}$ and $\texttt{xmax}$) for the function argument. Optionally, the user can also supply a maximun number of fucntion evaluations ($\texttt{itermax}$, 500 by default), the level of tolerance ($\texttt{tol}$, 1.0d-8 by default). Finally, the user can also control what it is printing during execution by setting the corresponding value of $\texttt{iprint}$:

- $\texttt{iprint}$ = 0: don't print anything (default)
- $\texttt{iprint}$ = 1: print warnings
- $\texttt{iprint}$ = 2: print warnings and every iteration

The subroutine returns the value of $\texttt{x}$ that maximizes $\texttt{func}$, the value of $\texttt{func}$ at $\texttt{x}$ ($\texttt{y}$),  and the number of function evaluations ($\texttt{numiter}$)

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