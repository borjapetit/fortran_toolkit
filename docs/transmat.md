
<span style="text-align:right;display:block;">
<a href="https://borjapetit.github.io/fortran_toolkit/">Back to index</a>
</span>

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

**Internal dependencies**: none
