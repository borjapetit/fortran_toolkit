## interpolation

```fortran
subroutine interpolation(pos,wth,xnow,xgrid)
  implicit none
  real(kind=8), intent(in)  :: xnow,xgrid(:)
  real(kind=8), intent(out) :: wth
  integer     , intent(out) :: pos
```

This subroutine finds the closest point in a given grid and return its position and the relative distance.

- If ```xnow```$<\min$(```xgrid```) the subroutine takes $\min($```xgrid```$)$ as the interpolated value of ```xnow```.
- If ```xnow```$>\max$(```xgrid```) the subroutine takes $\max($```xgrid```$)$ as the interpolated value of ```xnow```.

_Example_:

```fortran
! vector over which we wan to interpolate
vec = [1,2,3]

! interpolate x = 2.3 
call interpolation(pos,wth,2.3,vec)

! Results: pos = 3, wth = 0.3
! 2.3 = vec(pos)*wth + vec(pos-1)*(one-wth) = 3*0.3 + 2*0.7

! interpolate x = 0
call interpolation(pos,wth,0,vec)

! if the value is smaller than min(vec), it returns the min of vec
! Results: pos = 2, wth = 0.0
! 1 = vec(pos)*wth + vec(pos-1)*(one-wth) = 2*0.0 + 1*1.0
```

This subroutine is mainly used by the function [```interpolate```](interpolate.md).

**Dependencies**: none

[(back to index)](index.md)
