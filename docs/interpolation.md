---
exclude: true
---

```fortran
subroutine interpolation(pos,wth,xnow,xgrid)
  implicit none
  real(kind=8), intent(in)  :: xnow,xgrid(:)
  real(kind=8), intent(out) :: wth
  integer     , intent(out) :: pos
```

This subroutine finds the closest point in a given grid and return its position (```pos```) and the relative distance (```wth```) such that:

$$\texttt{xnow} = \texttt{vec}(\texttt{pos})\cdot\texttt{wth} \ + \ \texttt{vec}(\texttt{pos}-1)\cdot(\texttt{1.0d0}-\texttt{wth})$$

- If ```xnow``` < min(```xgrid```) the subroutine takes min(```xgrid```) as the interpolated value of ```xnow```.
- If ```xnow``` > max(```xgrid```) the subroutine takes max(```xgrid```) as the interpolated value of ```xnow```.

The vector ```xvec``` should be monotonic, so that 

$$\texttt{xgrid}(i) > \texttt{xgrid(i-1)} \ \ \ \text{or} \ \ \  \texttt{xgrid}(i) < \texttt{xgrid(i-1)}  \ \ \ \ \ \forall i\in[1,n]$$

This subroutine is mainly used by the function [```interpolate```](interpolate.md).

**Dependencies**: [```error```](error.md)

[(back to index)](../index.md)

---

**Example**

```fortran
! vector over which we wan to interpolate
vec = [1,2,3]

! interpolate x = 2.3 
call interpolation(pos,wth,2.3,vec)

! output: pos = 3, wth = 0.3
! 2.3 = vec(pos)*wth + vec(pos-1)*(one-wth) = 3*0.3 + 2*0.7

! interpolate x = 0
call interpolation(pos,wth,0,vec)

! since 0 is smaller than min(vec), the subroutine returns the min of vec
! output: pos = 2, wth = 0.0
! 1 = vec(pos)*wth + vec(pos-1)*(one-wth) = 2*0.0 + 1*1.0
```

