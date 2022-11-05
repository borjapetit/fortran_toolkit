## varvar

```fortran
function varvar(var,w,mask) result(variance)
  implicit none
  real(kind=8)            :: variance
  real(kind=8)            :: var(:)
  real(kind=8) , optional :: w(:)     ! same length of "var"
  logical      , optional :: mask(:)  ! same length of "var"
```

This function returns the variance of a variable ```var``` given some (optional) weigths ```w```. The user can also supply as ```mask``` to compute the conditional mean. If supplied, the vector ```w``` should have the same size as ```var```. If not supplied, the program assums uniform weigthing.

_Example_: compute the variance of a vector

```fortran
! without weigths
xvar = (/ 1.0, 4.0, 4.0, 9.0 /)
mean = varmean(xvar)  ! mean = 4.5
var  = varvar(xvar)   ! var  = 11.0
```

**Dependencies**: [```varmean```](varmean.md)

[(back to index)](inicio.md)