
## ```inverse```

```fortran
function inverse(mat) result(imat)
  implicit none
  real(kind=8) :: mat(:,:)
  real(kind=8) :: imat(size(m,1),size(m,1))
```

This function returns the inverse of a squared matrix ```mat``` implementing a LU decomposition.

**Dependencies**: [```error```](error.md)

[(back to index)](../index.md)

---