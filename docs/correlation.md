
<span style="text-align:right;display:block;">
<a href="https://borjapetit.github.io/fortran_toolkit/">Back to index</a>
</span>

## ```correlation```

```fortran
function correlation(xvar1,xvar2,w,mask) result(corr)
  implicit none
  real(kind=8)            :: corr      ! output: correlation coefficient between xvar1 and xvar2
  real(kind=8)            :: xvar1(:)  ! input: first variable
  real(kind=8)            :: xvar2(:)  ! input: second variable
  real(kind=8) , optional :: w(:)      ! input: (optional) weights
  logical      , optional :: mask(:)   ! input: (optional) restrictions on observations
```

This function returns the correlation coefficient between two (same size) vectors `xvar1` and `xvar2` given some (optional) weights `w`:
- If supplied, the vector `w` should have the same size as `xvar1` and `xvar2`.
- If not supplied, the program assumes uniform weighting.

The user can also supply a `mask` to compute the conditional correlation. The variable `mask` is a logical array of the same size of `xvar1` and `xvar2` indicating which observations should be considered.

**Arrays**: ```xvar1``` and ```xvar2``` (and ```w``` and ```mask```) can also be 2-dimensional arrays, all of the same shape. In that case the function uses all their elements, as if they were vectorized with [```vect```](vect.md): ```correlation(vect(xvar1),vect(xvar2),w=vect(w),mask=vect(mask))```.

If either variable is constant (zero standard deviation), the function prints an error message and returns 0.

**Internal dependencies**: [`varmean`](varmean.md),  [`varstd`](varstd.md), [`error`](error.md), [`vect`](vect.md)

---

**Example**

```fortran
corr = correlation(x1,x2,mask = x1.gt.0.0d0 .and. x2.lt.5.0d0)
```

This computes the correlation between vectors `x1` and `x2` only for those pair of points that satify the condition `x1`$>0$ and `x2`$<5$.





