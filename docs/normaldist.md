## normaldist

```fortran
subroutine normaldist(xvec,mu,sigma,n,dist)
  implicit none
  integer     , intent(in)  :: n
  real(kind=8), intent(in)  :: mu,sigma,xvec(n)
  real(kind=8), intent(out) :: dist(n)

  ! Internal dependencies: cdfn
```

This function returns the distribution of a normal random variable with mean ```mu``` and standadrd deviation ```sigma```.

[(back to index)](inicio.md)