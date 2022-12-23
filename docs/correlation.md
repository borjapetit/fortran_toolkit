## correlation

```fortran
function correlation(xvar1,xvar2,w,mask) result(corr)
  implicit none
  real(kind=8)            :: corr
  real(kind=8)            :: xvar1(:),xvar2(:)  ! both vectors should have the same length
  real(kind=8) , optional :: w(:)               ! same length of "xvar1" and "xvar2"
  logical      , optional :: mask               ! same length of "xvar1" and "xvar2"

  ! Internal dependencies: varmean, varstd
```

This function returns the correlation coefficient between two variables ```xvar1``` and ```xvar2``` given some (optional) weigths ```w```. The user can also supply a ```mask``` to compute the conditional stadard deviation. If supplied, the vector ```w``` should have the same size as ```var```. If not supplied, the program assumes uniform weigthing.

[(back to index)](inicio.md)