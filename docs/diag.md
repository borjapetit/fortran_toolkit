## diag

```fortran
function diag(mat) result(vec)
  implicit none
  real(kind=8) :: mat(:,:),vec(size(mat,dim=1))

  ! Dependencies: none
```

This function returns the main diagonal of a matric ```mat```.

_Example_:

```fortran
mat(:,1) = (/ 1.0, 2.0, 1.0 /)
mat(:,2) = (/ 3.0, 3.0, 4.0 /)
mat(:,3) = (/ 5.0, 1.0, 3.0 /)

vec = diag(mat)

print * , 'vec = ', vec   ! vec =  1.00  3.00  3.00
```

[(back to index)](inicio.md)