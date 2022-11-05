## varstd

```fortran
function varstd(var,w,mask) result(stdvar)
  implicit none
  real(kind=8)            :: stdvar
  real(kind=8)            :: var(:)
  real(kind=8) , optional :: w(:)     ! same length as "var"
  logical      , optional :: mask(:)  ! same length as "var"

  ! Internal dependencies: varvar
```

This function returns the stadard deviation of a variable ```var``` given some (optional) weigths ```w```. The user can also supply as ```mask``` to compute the conditional stadard deviation. If supplied, the vector ```w``` should have the same size as ```var```. If not supplied, the program assums uniform weigthing.

[(back to index)](inicio.md)