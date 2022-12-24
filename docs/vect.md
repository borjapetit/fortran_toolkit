# vect

```fortran
function vec(mat) result(vec)
  implicit none
  real(kind=8)  :: mat(:,:,...,:)
  real(kind=8)  :: vec(:)
```

This function returns a 1-dimensional array ```vec``` with all the elements of a user-supplied ```n```-dimensional array ```mat```, where ```n```$\leq5$.

_Example_:

```fortran
mat(:,1) = (/ 1 , 2 /)
mat(:,2) = (/ 3 , 4 /)

vec = vect(mat)

print * , 'vec =', vec   ! vec =  1.00  3.00  2.00  4.00
```

**Dependencies**: none

[(back to index)](index.md)
