
## interpolation

```fortran
subroutine interpolation(pos,wth,xnow,xgrid)
  implicit none
  real(kind=8), intent(in)  :: xnow,xgrid(:)
  real(kind=8), intent(out) :: wth
  integer     , intent(out) :: pos
```

This subroutine finds the closest point in a given grid and return its position and the relative distance.

_Example_: consider a vector ```vec = [1,2,3]```. We can use ```interpolation``` to find the closest point of ```xnow = 2.3```:

```fortran
call interpolation(pos,wth,2.3,vec)

! Results: pos = 3, wth = 0.3
! 2.3 = vec(pos)*wth + vec(pos-1)*(one-wth) = 3*0.3 + 2*0.7
```

This subroutine is mainly used by the function [```interpolate```](interpolate.md).

**Dependencies**: none

[(back to index)](index.md)
