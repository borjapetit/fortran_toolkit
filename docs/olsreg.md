## olsreg

```fortran
subroutine olsreg(coeffs,yvec,x1vec,x2vec,...,x8vec,w,mask,iprint)
  implicit none
  real(kind=8) , intent(out)           :: coeffs(:)
  real(kind=8) , intent(in)            :: yvec(:)
  real(kind=8) , intent(in)            :: x1vec(:)
  real(kind=8) , intent(in) , optional :: x2vec(:),x3vec(:),...,x8vec(:)
  real(kind=8) , intent(in) , optional :: w(:)     ! same length of "yvec"
  logical      , intent(in) , optional :: mask(:)  ! same length of "yvec"
  integer      , intent(in) , optional :: iprint

  ! Internal dependencies: varmean, varvar
```

[(back to index)](inicio.md)