
## ```trasnmat```

```fortran
function transmat(mat) result(matt)
  implicit none
  real(dp) :: mat(:,:)
  real(dp) :: matt(size(mat,2),size(mat,1))
```

This function returns the transpose of ```mat``` by applying:
```fortran
forall (i=1:size(mat,1)) matt(:,i) = mat(i,:)
```

**Dependencies**: none

[(back to index)](../index.md)
