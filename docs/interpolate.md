## interpolate

```fortran
function interpolate(x1,x2,...,xn,y1,y2,...,yn,mat) result(xi)
  implicit none
  real(kind=8) :: xi
  real(kind=8) :: x1,x2,...,xn
  real(kind=8) :: y1(:),y2(:),...,yn(:)
  real(kind=8) :: mat(:,:,...,:)
```

This function returns the linearly interpolated value of an n-dimensional function. The variables ```x1```, ```x2```, ..., ```xn``` are the values of the variables to be interpolated over their coresponding grids ```y1```, ```y2```, ...., ```yn```, and ```mat``` is an n-dimensional array with the results.

- The array ```mat``` must have at most, dimension 6 (so ```n```$\leq6$).
- If ```xn```$<\min$(```yn```) the subroutine takes $\min($```yn```$)$ as the value of ```x```.
- If ```xn```$>\max$(```yn```) the subroutine takes $\max($```yn```$)$ as the value of ```x```.

_Example 1_: We have a 2-dimensional array ```mat``` whose $(i,j)$-element is the value of some function $f$ evaluated at the $i$-element of the vector $x$, and the $j$-element of the vector $y$. Then, we can interpolate the value of $f$ at $(x_0,y_0)$ as:

```fortran
xi = interpolate(x_0,y_0,x,y,mat)
```

_Example 2_: We have a 3-dimensional array ```mat``` whose $(i,j,k)$-element is the value of some function $f$ evaluated at the $i$-element of the vector $x$, the $j$-element of the vector $y$, and the $k$-element of the vecor $z$. Then, we can interpolate the value of $f$ at $(x_0,y_0,z_0)$ as:

```fortran
xi = interpolate(x_0,y_0,z_0,x,y,z,mat)
```

**Dependencies**: [```interpolation```](interpolation.md)

[(back to index)](index.md)