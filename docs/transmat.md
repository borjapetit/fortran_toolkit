## transmat

```fortran
function transmat(mat) result(matt)
  implicit none
  real(kind=8) , intent(in)  :: mat(:,:)
  real(kind=8)               :: matt(size(mat,2),size(mat,1))
```

_Dependencies_: none

This function returns the transpose of a matrix ```mat```.

[(back to index)](inicio.md)