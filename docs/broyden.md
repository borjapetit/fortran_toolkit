
<span style="text-align:right;display:block;">
<a href="https://borjapetit.github.io/fortran_toolkit/">Back to index</a>
</span>

## ```broyden```

```fortran
subroutine broyden(jac,x1,x0,f1,f0)
  implicit none
  real(kind=8) , intent(inout) :: jac(:,:)  ! input: jacobian at x = x0; output: updated jacobian at x = x1
  real(kind=8) , intent(in)    :: x1(:)     ! input: new value of x
  real(kind=8) , intent(in)    :: x0(:)     ! input: previous value of x
  real(kind=8) , intent(in)    :: f1(:)     ! input: value of function at x = x1
  real(kind=8) , intent(in)    :: f0(:)     ! input: value of function at x = x0
```

This subroutine applies the Broyden's method to update a Jacobian matrix. You can learn more about this method in this [link](https://en.wikipedia.org/wiki/Broyden%27s_method). In particular, this subroutine returns an approximation to the Jacobian matrix around the point `x = x1`, without function evaluations:

$$J_1 = J_0 + \frac{\left(\Delta f - J_0\,\Delta x\right)\Delta x'}{\Delta x'\Delta x}, \ \ \ \ \Delta x = x_1-x_0, \ \ \Delta f = f_1-f_0$$

The user must supply a pair of points `x0` and `x1`, the value of the function evaluated at those points `f0` and `f1`, and the Jacobian matrix `jac` evaluated at `x0`, which is overwritten with the updated Jacobian. The sizes must satisfy `size(jac,1) = size(f0) = size(f1)` and `size(jac,2) = size(x0) = size(x1)`.

If the dimensions do not match, or if the step `x1-x0` is (close to) zero, the subroutine prints an error message and leaves `jac` unchanged.

*Note*: [`lmmin`](lmmin.md) uses its own internal copy of this update.

**Internal dependencies**: [`error`](error.md)
