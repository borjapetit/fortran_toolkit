

# Fortran toolkit

### Borja Petit

---

## ```correlation```

```fortran
function correlation(xvar1,xvar2,w,mask) result(corr)
  implicit none
  real(kind=8)            :: corr      ! output: correlation coefficient between xvar1 and xvar2
  real(kind=8)            :: xvar1(:)  ! input: first variable
  real(kind=8)            :: xvar2(:)  ! input: second variable
  real(kind=8) , optional :: w(:)      ! input: (optional) weights
  logical      , optional :: mask      ! input: (optional) retrictions on observations
```

This function returns the correlation coefficient between two (same size) vectors `xvar1` and `xvar2` given some (optional) weights `w`. If supplied, the vector `w` should have the same size as `xvar1` and `xvar2`. If not supplied, the program assumes uniform weigthing.

The user can also supply a `mask` to compute the conditional correlation. The variable `mask` is a logical array of the same size of `xvar1` and `xvar2` indicating which observations should be considered.

**Dependencies**: [`varmean`](varmean.md),  [`varvar`](varvar.md)

[(back to index)](../index.md)

---

**Example**

```fortran
corr = correlation(xvar1,xvar2,mask = xvar1.gt.0.0d0 .and. xvar2.lt.5.0d0)
```

This computes the correlation between vectors `xvar1` and `xvar2` only for those pair of points that satify the condition `xvar1}>0$ and `xvar2}<5$.





