

# Fortran toolkit

### Borja Petit

---

## ```broyden```

```fortran
subroutine broyden(j1,j0,x1,x0,f1,f0)
  implicit none
  real(kind=8) , intent(in)  :: x1(:)    ! input: initial value of x
  real(kind=8) , intent(in)  :: x0(:)    ! input: final value of x
  real(kind=8) , intent(in)  :: f1(:)    ! input: value of function at x = x0
  real(kind=8) , intent(in)  :: f0(:)    ! input: value of function at x = x1
  real(kind=8) , intent(in)  :: j0(:,:)  ! input: jacobian at x = x0
  real(kind=8) , intent(out) :: j1(:,:)  ! output: updated jacobian
```

This subroutine applies the Boryden's method to update a Jacobian matrix. You can learn more about this method in this [link](https://en.wikipedia.org/wiki/Broyden%27s_method). In particular, this subroutine returns an opproximation to the jacobian matrix around the point `x = x1`, without function evaluations.

The user must supply a pair of points `x0` and `x1`, the value of the function evaluated at those points `f0` and `f1`, and the jacobian matrix `j0` evaluated at `x0`.

*Note*: this function is used in [`lmmin`](lmmin.md)

**Dependencies**: none

[(back to index)](../index.md)

--
