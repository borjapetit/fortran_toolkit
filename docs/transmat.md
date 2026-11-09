
<span style="text-align:right;display:block;">
<a href="https://borjapetit.github.io/fortran_toolkit/">Back to index</a>
</span>

## ```transmat```

```fortran
function transmat(mat) result(matt)
  implicit none
  real(kind=8) :: mat(:,:)
  real(kind=8) :: matt(size(mat,2),size(mat,1))
```

This function returns the transpose of ```mat``` (which does not need to be square), equivalent to the intrinsic ```transpose```, by applying:
```fortran
forall (i=1:size(mat,1)) matt(:,i) = mat(i,:)
```

**Internal dependencies**: none
