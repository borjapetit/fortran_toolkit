# normaldist

```fortran
subroutine normaldist(xvec,mu,sigma,n,dist)
  implicit none
  integer     , intent(in)  :: n         ! number of points in the grid for "x"
  real(kind=8), intent(in)  :: xvec(n)   ! grid of points
  real(kind=8), intent(in)  :: mu        ! mean
  real(kind=8), intent(in)  :: sigma     ! standard deviation
  real(kind=8), intent(out) :: dist(n)   ! vector with the distribution
```

This function returns the distribution of a normal random variable with mean ```mu``` and standadrd deviation ```sigma```, over a grid ```xvec``` of size ```n```.

**Dependencies**: [```cdfn```](cdfn.md)

[(back to index)](index.md)
