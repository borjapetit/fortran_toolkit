
<span style="text-align:right;display:block;">
<a href="https://borjapetit.github.io/fortran_toolkit/">Back to index</a>
</span>

## ```inverse```

```fortran
function inverse(m,method) result(im)
  implicit none
  real(kind=8)     , intent(in)            :: m(:,:)
  character(len=*) , intent(in) , optional :: method
  real(kind=8)                             :: im(size(m,1),size(m,1))
```

This function returns the inverse of a square matrix ```m```. The input matrix is not modified. The optional argument ```method``` selects the algorithm:

- ```method``` = 'gj' (default): Gauss-Jordan elimination with partial pivoting (```inverse_gj```).
- ```method``` = 'qr': QR decomposition with column pivoting (```inverse_qr```, which solves ```m x = e_i``` with ```qrsolve``` for each column). More robust for ill-conditioned matrices; returns a stable pseudo-inverse when the matrix is numerically rank-deficient.
- ```method``` = 'svd': one-sided Jacobi singular value decomposition (```inverse_svd```). The most robust option: returns the Moore-Penrose pseudo-inverse, ignoring singular values below ```tol*svmax*n```.

Note that 'qr' and 'svd' treat a nearly singular matrix as rank-deficient, so for very ill-conditioned matrices (e.g. a Hilbert matrix of order 8) they return a pseudo-inverse and ```m*im``` is not the identity.

With ```method``` = 'gj', if ```m``` is not square, or if it is (numerically) singular, the function prints an error message (see [```error```](error.md)) and returns a matrix of zeros.

**Internal dependencies**: [```error```](error.md)

---