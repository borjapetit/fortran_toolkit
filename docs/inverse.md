
<span style="text-align:right;display:block;">
<a href="https://borjapetit.github.io/fortran_toolkit/">Back to index</a>
</span>

## ```inverse```

```fortran
function inverse(m) result(im)
  implicit none
  real(kind=8) , intent(in) :: m(:,:)
  real(kind=8)              :: im(size(m,1),size(m,1))
```

This function returns the inverse of a square matrix ```m``` using Gauss-Jordan elimination with partial pivoting. The input matrix is not modified.

If ```m``` is not square, or if it is (numerically) singular, the function prints an error message (see [```error```](error.md)) and returns a matrix of zeros.

**Internal dependencies**: [```error```](error.md)

---