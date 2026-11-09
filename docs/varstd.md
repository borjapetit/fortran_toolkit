
<span style="text-align:right;display:block;">
<a href="https://borjapetit.github.io/fortran_toolkit/">Back to index</a>
</span>

## ```varstd```

```fortran
function varstd(var,w,mask) result(stdvar)
  implicit none
  real(kind=8)            :: stdvar
  real(kind=8)            :: var(:)
  real(kind=8) , optional :: w(:)   
  logical      , optional :: mask(:)
```

This function returns the standard deviation of a variable ```var``` given some (optional) weights ```w```. If supplied, the vector ```w``` should have the same size as ```var```. If not supplied, the program assumes uniform weighting. The standard deviation is the square root of [```varvar```](varvar.md), i.e., the population standard deviation (the weighted sum of squared deviations is divided by the sum of weights, without degrees-of-freedom correction).

The user can also supply a ```mask``` to compute the conditional standard deviation. The input ```mask``` is a logical array of the same size of ```var```.

**Arrays**: ```var``` (and ```w``` and ```mask```) can also be 2-dimensional arrays, all of the same shape. In that case the function uses all their elements, as if they were vectorized with [```vect```](vect.md): ```varstd(vect(var),w=vect(w),mask=vect(mask))```.

**Internal dependencies**: [```error```](error.md), [```varmean```](varmean.md),  [```varvar```](varvar.md), [```vect```](vect.md)

---

**Example**

```fortran
standdesv = varstd(var,mask = var.gt.0.0d0 .and. var.lt.5.0d0)
```

This computes the standard deviation of ```var``` conditional on ```var``` being between 0 and 5.
