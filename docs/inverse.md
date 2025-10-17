
<span style="text-align:right;display:block;">
<a href="https://borjapetit.github.io/fortran_toolkit/">Back to index</a>
</span>

## ```inverse```

```fortran
function inverse(mat) result(imat)
  implicit none
  real(kind=8) :: mat(:,:)
  real(kind=8) :: imat(size(m,1),size(m,1))
```

This function returns the inverse of a squared matrix ```mat``` implementing a LU decomposition.

**Internal dependencies**: [```error```](error.md)

---