# percentile

```fortran
function percentile(xvec,pct,w,mask) result(cutoff)
  implicit none
  real(kind=8)            :: cutoff
  real(kind=8)            :: xvec(:)
  real(kind=8)            :: pct
  real(kind=8) , optional :: w(:)     ! same length of "xvec"
  logical      , optional :: mask(:)  ! same length of "xvec"
```

This function returns the percentile ```pct``` for a distribution ```xvec```, given some (optional) weigths ```w```. the percentile ```pct``` should be between 0 and 1. If supplied, the vector ```w``` should have the same size as ```var```. If not supplied, the program assums uniform weigthing.

The user can also supply a ```mask``` to compute the conditional percentile. The input ```mask``` is a logical array of the same size of ```xvec```. For example:

```fortran
pc60 = percentile(xvec,0.60d0,mask = var.gt.0.0d0 .and. var.lt.5.0d0)
```

This computes the 60th percentile of ```xvec``` conditional on ```xvec``` being between 0 and 5.

**Dependencies**: none

[(back to index)](index.md)
