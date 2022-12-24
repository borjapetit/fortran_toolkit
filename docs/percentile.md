## percentile

```fortran
function percentile(xvec,pct,w,mask) result(cutoff)
  implicit none
  real(kind=8)            :: cutoff
  real(kind=8)            :: xvec(:)
  real(kind=8)            :: pct
  real(kind=8) , optional :: w(:)     ! same length of "xvec"
  logical      , optional :: mask(:)  ! same length of "xvec"

  ! Internal dependencies: none
```

This function returns the percentile ```pct``` for a distribution ```xvec```, given some (optional) weigths ```w```. The user can also supply as ```mask``` to compute the conditional correlation. If supplied, the vector ```w``` should have the same size as ```var```. If not supplied, the program assums uniform weigthing.

_Example_: given a vector ```xvec``` with a sample of a variable ```x```, find the 60th percentile:

```fortran
pc60 = percentile(xvec,60.0d0)
```

[(back to index)](inicio.md)