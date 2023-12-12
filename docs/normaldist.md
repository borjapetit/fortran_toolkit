
# Fortran toolkit

###### Borja Petit

---

## ```normaldist```

```fortran
subroutine normaldist(xvec,mu,sigma,n,dist)
  implicit none
  integer     , intent(in)  :: n         ! number of points in the grid for "x"
  real(kind=8), intent(in)  :: xvec(n)   ! grid of points
  real(kind=8), intent(in)  :: mu        ! mean
  real(kind=8), intent(in)  :: sigma     ! standard deviation
  real(kind=8), intent(out) :: dist(n)   ! vector with the distribution
```

This subroutine discretized the distribution of a normal random variable with mean ```mu``` and standard deviation ```sigma```, over a grid ```xvec``` of size ```n```. The user should provide a vector with values of the variable $x$ over which to compute the distribution. For any $i\in[2,\texttt{n}-1]$:

$$\texttt{dist}(i) = \Phi\left( \frac{\texttt{z}(i+1)+\texttt{z}(i)}{2}\right) - \Phi \left( \frac{\texttt{z}(i)+\texttt{z}(i-1)}{2}\right)$$

where $\Phi(\cdot)$ is the cdf of a standard normal distribution and $\texttt{z} = (\texttt{xvec}-\texttt{mu})/\texttt{sigma}$ is a vector with the normalized values of ```xvec```. And for $i\in\{1,\texttt{n}\}$:

$$\texttt{dist}(1) = \Phi\left( \frac{\texttt{z}(2)+\texttt{z}(1)}{2}\right) \ \ \ \ \ \ \text{and} \ \ \ \ \ \ 
\texttt{dist}(n) = 1 - \Phi\left( \frac{\texttt{z}(\texttt{n})+\texttt{z}(\texttt{n}-1)}{2}\right) $$

**Dependencies**: [```cdfn```](cdfn.md)

[(back to index)](../index.md)

---