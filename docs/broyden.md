
<span style="text-align:right;display:block;">
<a href="https://borjapetit.github.io/fortran_toolkit/">Back to index</a>
</span>

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

**Internal dependencies**: [`error`](error.md)

---

#### Code:

```fortran
  subroutine broyden(j1,j0,x1,x0,f1,f0)
    implicit none
    real(dp) , intent(in)  :: x1(:),f1(:),x0(:),f0(:),j0(:,:)
    real(dp) , intent(out) :: j1(:,:)
    real(dp)               :: df(size(f0,1))
    real(dp)               :: dx(size(x1,1))
    real(dp)               :: den
    integer                :: i
    if (size(f1).ne.size(j0,1)) then
      call error(' error in broyden: size(y) /= size(j,1)')
      j1 = j0
      return
    end if
    if (size(x1).ne.size(j0,2)) then
      call error(' error in broyden: size(x) /= size(j,2)')
      j1 = j0
      return
    end if
    den = sum(dx(:)*dx(:))
    if (abs(den).lt.tolvl) then
      call error(' error in broyden: denominator close to 0')
      j1 = j0
      return
    end if
    dx(:) = x1(:)-x0(:)
    df(:) = f1(:)-f0(:)
    do i=1,size(df)
      j1(i,:) = j0(i,:) + dx(:)*( df(i)-sum(j0(i,:)*dx(:)) ) / den
    end do
    return
  end subroutine broyden
```

