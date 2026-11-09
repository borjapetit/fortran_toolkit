
<span style="text-align:right;display:block;">
<a href="https://borjapetit.github.io/fortran_toolkit/">Back to index</a>
</span>

## ```cdfn```

```fortran
elemental function cdfn(x) result(f)
  implicit none
  real(kind=8), intent(in) :: x   ! input
  real(kind=8)             :: f   ! output: cdf of standard normal at x
```

This function returns the cdf of a standard normal distribution, $\Phi(x) = P(X\leq x)$ with $X\sim N(0,1)$, using the algorithm AS 66 (Hill, 1973). The function is defined as `elemental`, so `x` can be a scalar or an array (the output has the same shape as `x`).

**Internal dependencies**: none

---

**Example**

```fortran
print * , 'Result =', cdfn(0.0d0)                    ! Result =   0.500
print * , 'Result =', cdfn(-1.0d0)                   ! Result =   0.158
print * , 'Result =', cdfn( (/ 0.0d0 , -1.0d0 /) )   ! Result =   0.500   0.158
```