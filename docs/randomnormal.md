
<span style="text-align:right;display:block;">
<a href="https://borjapetit.github.io/fortran_toolkit/">Back to index</a>
</span>

## ```randomnormal```

```fortran
subroutine randomnormal(shock,mu,std)
  implicit none
  real(kind=8) , intent(in)  :: mu,std
  real(kind=8) , intent(out) :: shock ! or shock(:)
```

This subroutine returns a random number/vector drawn from a normal distribution with mean ```mu``` and standard deviation ```std```, using the ratio-of-uniforms method of Leva (1992). The output, ```shock```, can either be a scalar or a vector of dimension-$n$ (it does not work with matrices). In the vector case, the elements are independent draws (the sample is not re-centered, so its sample mean is not exactly ```mu```).

The subroutine uses the intrinsic ```random_number```, so the user can fix the seed with the intrinsic ```random_seed```.

**Note**: ```randomnormal``` is an interface that calls ```randomnormal_scalar``` or ```randomnormal_vec``` depending on whether ```shock``` is a scalar or a vector.

**Internal dependencies**: none
