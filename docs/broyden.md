## broyden

```fortran
subroutine broyden(j1,j0,x1,x0,f1,f0)
  implicit none
  real(kind=8) , intent(in)  :: x1(:),f1(:)
  real(kind=8) , intent(in)  :: x0(:),f0(:)
  real(kind=8) , intent(in)  :: j0(:,:)
  real(kind=8) , intent(out) :: j1(:,:)
```

_Dependencies_: none

This subroutine applies the Boryden's method to update a Jacobian matrix.

Imagine we have an $m$-dimensional function $f$ in $n$ unknows. We evaluate two points $x_0$ and $x_1$, $f_1 = f(x_1)$ and $f_0 = f(x_0)$, and we compute the numerical jacobian of the function $f$ around $x=x_0$. This subroutine returns an opproximation to the jacobian matrix arounf the point $x=x_1$. The user must supply a pair of points ```x0``` and ```x1```, the value of the function evaluated at thos epoints ```f0``` and ```f1```, and the jacobian matrix ```j0``` evaluated at ```x0```.

You can learn more about this method in this [link](https://en.wikipedia.org/wiki/Broyden%27s_method).

[(back to index)](inicio.md)