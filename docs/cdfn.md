## cdfn

```fortran
elemental function cdfn(x) result(f)
  implicit none
  real(kind=8), intent(in) :: x
  real(kind=8)             :: f

  ! Internal dependencies: none
```

This function returns the cdf of a standard normal distribution, ```f``` $=\Phi(x)$. This subroutine is defined as ```elemental```, which implies that it can be call for both scalars and arrays.

_Example 1_: Scalar output

```fortran
print * , 'Result =', cdf(0.0)    ! Result = 0.500
print * , 'Result =', cdf(-1.0)   ! Result = 0.158
```

_Example 2_: Vector otput

```fortran
vec = (/ 0.0 , -1.0 /)

print * , 'Result =', cdf(vec)   ! Result = 0.500 0.158
```

[(back to index)](inicio.md)