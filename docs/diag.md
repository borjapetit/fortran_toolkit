
<span style="text-align:right;display:block;">
<a href="https://borjapetit.github.io/fortran_toolkit/">Back to index</a>
</span>

## ```diag```

```fortran
function diag(mat) result(vec)
  implicit none
  real(kind=8) :: mat(:,:)              ! input
  real(kind=8) :: vec(size(mat,dim=1))  ! output: vector with main diagonal of mat
```

This function returns the main diagonal of a square matrix `mat`. If `mat` is not square, the function prints an error message (see [`error`](error.md)).

**Internal dependencies**: [```error```](error.md)

---

**Example**

```fortran
mat(:,1) = (/ 1.0d0, 2.0d0, 1.0d0 /)
mat(:,2) = (/ 3.0d0, 3.0d0, 4.0d0 /)
mat(:,3) = (/ 5.0d0, 1.0d0, 3.0d0 /)

vec = diag(mat)

print * , 'vec = ', vec   ! vec =  1.00  3.00  3.00
```
