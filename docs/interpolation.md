
<span style="text-align:right;display:block;">
<a href="https://borjapetit.github.io/fortran_toolkit/">Back to index</a>
</span>

## ```interpolation```

```fortran
subroutine interpolation(pos,wth,xnow,xgrid)
  implicit none
  real(kind=8), intent(in)  :: xnow,xgrid(:)
  real(kind=8), intent(out) :: wth
  integer     , intent(out) :: pos
```

This subroutine finds the closest point in a given grid and return its position (```pos```) and the relative distance (```wth```) such that:

```
xnow = vec(pos)·wth + vec(pos-1)·(1-wth)
```

- If ```xnow``` < min(```xgrid```) the subroutine takes min(```xgrid```) as the interpolated value of ```xnow```.
- If ```xnow``` > max(```xgrid```) the subroutine takes max(```xgrid```) as the interpolated value of ```xnow```.

The vector ```xgrid``` should be monotonically increasing, so that ```xgrid(i) > xgrid(i-1)``` is satisfied for any `i`.

If ```xgrid``` has a single element, the subroutine returns ```pos``` = 1 and ```wth``` = 1. If ```xnow``` is NaN, it prints an error message and returns a valid position with ```wth``` = NaN, so that the interpolated value is NaN. If ```xgrid``` is empty, it prints an error message and returns ```pos``` = 0.

This subroutine is mainly used by the function [```interpolate```](interpolate.md).

**Internal dependencies**: [```error```](error.md)

---

**Example**

```fortran
! vector over which we wan to interpolate
vec = [1.0d0,2.0d0,3.0d0]

! interpolate x = 2.3 
call interpolation(pos,wth,2.3d0,vec)

! output: pos = 3, wth = 0.3
! 2.3 = vec(pos)*wth + vec(pos-1)*(one-wth) = 3*0.3 + 2*0.7

! interpolate x = 0
call interpolation(pos,wth,0.0d0,vec)

! since 0 is smaller than min(vec), the subroutine returns the min of vec
! output: pos = 2, wth = 0.0
! 1 = vec(pos)*wth + vec(pos-1)*(one-wth) = 2*0.0 + 1*1.0
```

