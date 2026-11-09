
<span style="text-align:right;display:block;">
<a href="https://borjapetit.github.io/fortran_toolkit/">Back to index</a>
</span>

## ```varmean```

```fortran
function varmean(var,w,mask) result(meanvar)
  implicit none
  real(kind=8)            :: meanvar
  real(kind=8)            :: var(:)
  real(kind=8) , optional :: w(:)     
  logical      , optional :: mask(:)  
```

This function returns the (weighted) mean of a variable ```var```, ```sum(w*var)/sum(w)```, given some (optional) weights ```w```. If supplied, the vector ```w``` should have the same size as ```var```. If not supplied, the program assumes uniform weighting.

The user can also supply a ```mask``` to compute the conditional mean. The input ```mask``` is a logical array of the same size of ```var```.

**Arrays**: ```var``` (and ```w``` and ```mask```) can also be 2-dimensional arrays, all of the same shape. In that case the function uses all their elements, as if they were vectorized with [```vect```](vect.md): ```varmean(vect(var),w=vect(w),mask=vect(mask))```.

**Internal dependencies**: [```error```](error.md), [```vect```](vect.md)

---

**Example**

```fortran
mean = varmean(var,mask = var.gt.0.0d0 .and. var.lt.5.0d0)
```

This computes the mean of ```var``` conditional on ```var``` being between 0 and 5.
