
<span style="text-align:right;display:block;">
<a href="https://borjapetit.github.io/fortran_toolkit/">Back to index</a>
</span>

## ```percentile```

```fortran
function percentile(xvec,pct,w,mask) result(cutoff)
  implicit none
  real(kind=8)            :: cutoff
  real(kind=8)            :: xvec(:)
  real(kind=8)            :: pct
  real(kind=8) , optional :: w(:)    
  logical      , optional :: mask(:) 
```

This function returns the percentile ```pct``` for a distribution ```xvec```, given some (optional) weights ```w```. If supplied, the vector ```w``` should have the same size as ```xvec```. If not supplied, the program assumes uniform weighting.

The value of ```pct``` can be given either as a fraction (0 $\leq$ ```pct``` $\leq$ 1) or as a percentage (1 < ```pct``` $\leq$ 100). Note that ```pct``` = 1 is read as a fraction (i.e., the maximum): for the 1st percentile use ```pct``` = 0.01.

The subroutine finds the percentile ```pct``` using a Bisection method.

The user can also supply a ```mask``` to compute the conditional percentile. The input ```mask``` is a logical array of the same size of ```xvec```.

**Arrays**: ```xvec``` (and ```w``` and ```mask```) can also be 2-dimensional arrays, all of the same shape. In that case the function uses all their elements, as if they were vectorized with [```vect```](vect.md): ```percentile(vect(xvec),pct,w=vect(w),mask=vect(mask))```.

**Internal dependencies**: [```error```](error.md), [```vect```](vect.md)

---

**Example**

```fortran
pc60 = percentile(xvec,0.60d0,mask = xvec.gt.0.0d0 .and. xvec.lt.5.0d0)
```

This computes the 60th percentile of ```xvec``` conditional on ```xvec``` being between 0 and 5 (```pct``` = 60.0d0 gives the same result).
