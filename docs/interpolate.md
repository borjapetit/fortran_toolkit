# interpolate

```fortran
function interpolate(x1,x2,...,xn,y1,y2,...,yn,mat) result(xi)
  implicit none
  real(kind=8) :: xi
  real(kind=8) :: x1,x2,...,xn
  real(kind=8) :: y1(:),y2(:),...,yn(:)
  real(kind=8) :: mat(:,:,...,:)
```

This function returns the linearly interpolated value of an n-dimensional function, with ```n```<7. The variables ```x1```, ```x2```, ..., ```xn``` are the values of the variables to be interpolated over their coresponding grids ```y1```, ```y2```, ...., ```yn```, and ```mat``` is an n-dimensional array with the results.

**Dependencies**: [```interpolation```](interpolation.md)

[(back to index)](../index.md)

---

**Example**

We have a 2-dimensional array ```mat``` whose (i,j)-element is the value of some function ```func``` evaluated at the i-element of the vector ```x```, and the j-element of the vector ```y```:

```fortran
do i=1,nx
  do j=1,ny
    mat(i,j) = func(x(i),y(j))
  end do
end do
```

Then, we can interpolate the value of ```func``` at ```(x_0,y_0)``` as:

```fortran
func0 = interpolate(x_0,y_0,x,y,mat)
```

If the array ```mat``` is 3-dimensional, we can interpolate the value of ```func``` at ```(x_0,y_0,z_0)``` as

```fortran
func0 = interpolate(x_0,y_0,z_0,x,y,z,mat)
```

