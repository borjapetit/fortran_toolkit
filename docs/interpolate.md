
### interpolate

```fortran
function interpolate(x1,x2,...,xn,y1,y2,...,yn,mat) result(xi)
  implicit none
  real(kind=8) :: x1    , x2    , ... , xn
  real(kind=8) :: y1(:) , y2(:) , ... , yn(:)   
  real(kind=8) :: mat(:,:,...,:)
```

This function returns the linearly interpolated value of an n-dimensional function, with ```n``` < 7. The variables ```x1```, ```x2```, ..., ```xn``` are the values of the variables to be interpolated over their coresponding grids ```y1```, ```y2```, ...., ```yn```, and ```mat``` is an n-dimensional array with the results.

**Dependencies**: [```error```](error.md), [```interpolation```](interpolation.md)

**Note**: ```interpolate``` is an interface that calls specific functions depending on the dimension of ```mat```. These functions are ```interpolate1d```, ```interpolate2d```, ..., ```interpolate6d```.

[(back to index)](../index.md)

---

**Example**

We have a 2-dimensional array ```mat``` whose (i,j)-element is the value of some function ```func``` evaluated at the i-element of the vector ```xgrid```, and the j-element of the vector ```ygrid```:

```fortran
do i=1,nx
  do j=1,ny
    mat(i,j) = func(xgrid(i),ygrid(j))
  end do
end do
```

Then, we can interpolate the value of ```func``` at ```(x0,y0)``` as:

```fortran
func0 = interpolate(x0,y0,xgrid,ygrid,mat)
```

If the array ```mat``` is 3-dimensional, we can interpolate the value of ```func``` at ```(x0,y0,z0)``` as

```fortran
func0 = interpolate(x0,y0,z0,xgrid,ygrid,zgrid,mat)
```

