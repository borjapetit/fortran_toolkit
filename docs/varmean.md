## varmean

```fortran
function varmean(var,w,mask) result(meanvar)
  implicit none
  real(kind=8)            :: meanvar
  real(kind=8)            :: var(:)
  real(kind=8) , optional :: w(:)     ! same length of "var"
  logical      , optional :: mask(:)  ! same length of "var"
```

This function returns the mean of a variable ```var``` given some (optional) weigths ```w```. The user can also supply as ```mask``` to compute the conditional mean. If supplied, the vector ```w``` should have the same size as ```var```. If not supplied, the program assums uniform weigthing.

_Example_: compute the mean of a vector

```fortran
! without weigths
xvar = (/ 1.0, 4.0, 4.0, 9.0 /)
mean = varmean(xvar)  ! mean = 4.5

! without weigths and conditional on xvar>0
xvar = (/ 1.0, 4.0, 4.0, 9.0 /)
mean = varmean(xvar, mask = xvar.gt.2.0d0 )  ! mean = 5.66

! with weigths
xvar = (/ 1.0, 4.0, 4.0, 9.0 /)
wvar = (/ 2.0, 4.0, 5.0, 2.0 /)
mean = varmean(xvar,w = wvar)  ! mean = 4.3076
```
[(back to index)](inicio.md)