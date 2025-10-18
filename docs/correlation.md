
<span style="text-align:right;display:block;">
<a href="https://borjapetit.github.io/fortran_toolkit/">Back to index</a>
</span>

## ```correlation```

```fortran
function correlation(x1,x2,w,mask) result(corr)
  implicit none
  real(kind=8)            :: corr    ! output: correlation coefficient between x1 and x2
  real(kind=8)            :: x1(:)   ! input: first variable
  real(kind=8)            :: x2(:)   ! input: second variable
  real(kind=8) , optional :: w(:)    ! input: (optional) weights
  logical      , optional :: mask    ! input: (optional) retrictions on observations
```

This function returns the correlation coefficient between two (same size) vectors `x1` and `x2` given some (optional) weights `w`:
- If supplied, the vector `w` should have the same size as `x1` and `x2`.
- If not supplied, the program assumes uniform weigthing: 

The user can also supply a `mask` to compute the conditional correlation. The variable `mask` is a logical array of the same size of `x1` and `x2` indicating which observations should be considered.

**Internal dependencies**: [`varmean`](varmean.md),  [`varstd`](varstd.md), [`error`](error.md)

---

**Example**

```fortran
corr = correlation(x1,x2,mask = x1.gt.0.0d0 .and. x2.lt.5.0d0)
```

This computes the correlation between vectors `x1` and `x2` only for those pair of points that satify the condition `x1`$>0$ and `x2`$<5$.





