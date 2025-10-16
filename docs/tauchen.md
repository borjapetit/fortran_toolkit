

# Fortran toolkit

### Borja Petit
---

## ```tauchen```

```fortran
subroutine tauchen(xvec,rho,mu,sigma,n,pmat)
  implicit none
  integer     , intent(in)  :: n
  real(kind=8), intent(in)  :: rho,mu,sigma,xvec(n)
  real(kind=8), intent(out) :: pmat(n,n)
```

This function returns the transition matrix for a discretized AR(1) process of the form:

$$x' = \texttt{mu} + \texttt{rho} \cdot  x + \texttt{sigma} \cdot u ,  \ \ \ \ \ \text{with} u \sim N(0,1)$$

The vector with the values of $x$, ```xvec```, is of dimension ```n``` and does not need to be equally spaced.

**Dependencies**: [```normaldist```](normaldist.md)

[(back to index)](../index.md)

---

**Example**

Imagine a variable $x$ that follows an AR(1) procress with an autocorrelation of 0.8 and subject to normal shocks with stadard deviation 0.2.

```fortran
! create a 100-point equallly soaced grid for the variable "x"
xvec = grid( 3*0.20 , -3*0.20 , 100 )

! get the transition matrix for the discretized AR(1) process
call tauchen(xvec,0.80,0.00,0.20,100,pmat)
```
