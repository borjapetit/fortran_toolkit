## tauchen

```fortran
subroutine tauchen(xvec,rho,mu,sigma,n,pmat)
  implicit none
  integer     , intent(in)  :: n
  real(kind=8), intent(in)  :: rho,mu,sigma,xvec(n)
  real(kind=8), intent(out) :: pmat(n,n)

  ! Internal dependencies: normaldist
```

This function returns the transition matrix for a discretized AR(1) process of the form:
$$x' = \mu + \rho x + \sigma \epsilon, \hspace{0.2cm} \epsilon \sim N(0,1)$$
The vector with values of $x$, ```xvec```, if of dimension ```n``` and does not need to be equally spaced.

[(back to index)](inicio.md)