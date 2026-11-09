
<span style="text-align:right;display:block;">
<a href="https://borjapetit.github.io/fortran_toolkit/">Back to index</a>
</span>

## ```tauchen```

```fortran
subroutine tauchen(xvec,rho,mu,sigma,n,pmat)
  implicit none
  integer     , intent(in)  :: n
  real(kind=8), intent(in)  :: rho,mu,sigma,xvec(n)
  real(kind=8), intent(out) :: pmat(n,n)
```

This subroutine returns the transition matrix for a discretized AR(1) process of the form:

$$x' = \texttt{mu} + \texttt{rho} \cdot  x + \texttt{sigma} \cdot u ,  \ \ \ \ \ \text{with} \ u \sim N(0,1)$$

The vector with the values of $x$, ```xvec```, is of dimension ```n``` (```n``` $\geq 2$), should be increasing, and does not need to be equally spaced. The element ```pmat(i,j)``` is the probability of moving from ```xvec(i)``` to ```xvec(j)```; each row is computed with [```normaldist```](normaldist.md) with mean ```mu + rho*xvec(i)``` and standard deviation ```sigma```, and sums to 1.

Note that ```mu``` is the intercept of the process, not its unconditional mean, which is ```mu/(1-rho)```.

**Internal dependencies**: [```normaldist```](normaldist.md)

---

**Example**

Imagine a variable $x$ that follows an AR(1) process with an autocorrelation of 0.8 and subject to normal shocks with standard deviation 0.2.

```fortran
! create a 100-point equally spaced grid for the variable "x"
! (+/- 3 unconditional standard deviations: 0.20/sqrt(1-0.80**2) = 0.333)
xvec = grid( 1.0d0 , -1.0d0 , 100 )

! get the transition matrix for the discretized AR(1) process
call tauchen(xvec,0.80d0,0.00d0,0.20d0,100,pmat)
```
