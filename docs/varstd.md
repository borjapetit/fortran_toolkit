# varstd

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

_Example_: compute the standard deviation of a vector ```xvar = (/ 1.0, 4.0, 4.0, 9.0 /)```

```fortran
! without weigths
std = varstd(xvar)   ! var  = 3.316

! conditional on "x" being greater than 1
std = varsted(xvar , mask = xvar.gt.1.0d0 )   ! var = 2.887
```

**Dependencies**: [```varmean```](varmean.md),  [```varvar```](varvar.md)

[(back to index)](index.md)