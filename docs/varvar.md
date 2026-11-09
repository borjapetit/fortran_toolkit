
<span style="text-align:right;display:block;">
<a href="https://borjapetit.github.io/fortran_toolkit/">Back to index</a>
</span>

## ```varvar```

```fortran
function varvar(var,w,mask) result(variance)
  implicit none
  real(kind=8)            :: variance
  real(kind=8)            :: var(:)
  real(kind=8) , optional :: w(:)   
  logical      , optional :: mask(:)
```

This function returns the (population) variance of a variable ```var```, ```sum(w*(var-mean)**2)/sum(w)```, given some (optional) weights ```w```. If supplied, the vector ```w``` should have the same size as ```var```. If not supplied, the program assumes uniform weighting.

The user can also supply a ```mask``` to compute the conditional variance. The input ```mask``` is a logical array of the same size of ```var```. Observations with ```var``` = NaN are excluded from both the mean and the variance.

**Arrays**: ```var``` (and ```w``` and ```mask```) can also be 2-dimensional arrays, all of the same shape. In that case the function uses all their elements, as if they were vectorized with [```vect```](vect.md): ```varvar(vect(var),w=vect(w),mask=vect(mask))```.

**Internal dependencies**: [```error```](error.md), [```varmean```](varmean.md), [```vect```](vect.md)

---

**Example**

```fortran
variance = varvar(var,mask = var.gt.0.0d0 .and. var.lt.5.0d0)
```

This computes the variance of ```var``` conditional on ```var``` being between 0 and 5.
