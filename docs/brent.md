## brent

```fortran
subroutine brent(func,x,iy,ind,x0,x1,itermax,tol)
  implicit none
  external                             :: func
  real(kind=8) , intent(out)           :: x
  integer      , intent(out)           :: iy
  integer      , intent(out)           :: ind
  real(kind=8) , intent(in)            :: x0
  real(kind=8) , intent(in)            :: x1
  real(kind=8) , intent(in) , optional :: tol
  integer      , intent(in) , optional :: itermax
```

_Dependencies_: none

[(back to index)](inicio.md)