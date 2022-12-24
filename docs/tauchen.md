# tauchen

```fortran
subroutine tauchen(xvec,rho,mu,sigma,n,pmat)
  implicit none
  integer     , intent(in)  :: n
  real(kind=8), intent(in)  :: rho,mu,sigma,xvec(n)
  real(kind=8), intent(out) :: pmat(n,n)
```

This function returns the transition matrix for a discretized AR(1) process of the form:
$$x' = \texttt{mu} + \text{rho}\cdot x + \text{sigma}\cdot\epsilon, \hspace{0.2cm} \epsilon \sim N(0,1)$$
The vector with the values of $x$, ```xvec```, is of dimension ```n``` and does not need to be equally spaced.

**Dependencies**: [```normaldist```](normaldist.md)

[(back to index)](index.md)
