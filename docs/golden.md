## golden

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

_Dependencies_: none

[(back to index)](inicio.md)