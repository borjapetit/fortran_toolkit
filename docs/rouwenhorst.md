
<span style="text-align:right;display:block;">
<a href="https://borjapetit.github.io/fortran_toolkit/">Back to index</a>
</span>

## ```rouwenhorst_grid``` and ```rouwenhorst_transition```

```fortran
subroutine rouwenhorst_grid(n,rho,sigma_e,zgrid)
  implicit none
  integer      , intent(in)  :: n
  real(kind=8) , intent(in)  :: rho,sigma_e
  real(kind=8) , intent(out) :: zgrid(n)

subroutine rouwenhorst_transition(n,rho,pmat)
  implicit none
  integer      , intent(in)  :: n
  real(kind=8) , intent(in)  :: rho
  real(kind=8) , intent(out) :: pmat(n,n)
```

These subroutines discretize an AR(1) process ```z' = rho*z + e```, with ```e ~ N(0,sigma_e^2)```, using the Rouwenhorst method, which is more accurate than [```tauchen```](tauchen.md) for highly persistent processes.

- ```rouwenhorst_grid``` returns a grid of ```n``` equally spaced points in ```[-psi,psi]```, with ```psi = sigma_z*sqrt(n-1)``` and ```sigma_z = sigma_e/sqrt(1-rho^2)```.
- ```rouwenhorst_transition``` returns the ```n x n``` transition matrix (rows sum to one).

**Internal dependencies**: [```error```](error.md)

---

**Example**

```fortran
  real(kind=8) :: zgrid(7),pmat(7,7)
  call rouwenhorst_grid(7,0.95d0,0.1d0,zgrid)
  call rouwenhorst_transition(7,0.95d0,pmat)
```
