# randomnormal

```fortran
subroutine randomnormal(shock,mu,std)
  implicit none
  real(kind=8) , intent(in)  :: mu,std
  real(kind=8) , intent(out) :: shock ! or shock(:)
```

This function returns a random number draw from a distribution $N(\mu,\sigma)$. The output, ```shock```, can either be a scalar or a vector of dimension-$n$.

_Note_: ```randomnormal``` is an interface that calls ```randomnormal_scalar``` or ```randomnormal_vec``` depending on whether ```shock``` is a scalar or a vector.

**Dependencies**: none

[(back to index)](index.md)
