
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

This function returns the cdf of a standard normal distribution. This subroutine is defined as `elemental`.

**Internal dependencies**: none

---

**Example**

```fortran
print * , 'Result =', cdf(0.0)                  ! Result =   0.500
print * , 'Result =', cdf(-1.0)                 ! Result =   0.158
print * , 'Result =', cdf( (/ 0.0 , -1.0 /) )   ! Result =   0.500   0.158
```