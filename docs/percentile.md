
# Fortran toolkit

###### Borja Petit

---

## ```percentile```

```fortran
function percentile(var,pct,w,mask) result(cutoff)
  implicit none
  real(kind=8)            :: cutoff
  real(kind=8)            :: var(:)
  real(kind=8)            :: pct
  real(kind=8) , optional :: w(:)    
  logical      , optional :: mask(:) 
```

This function returns the percentile ```pct``` for a distribution ```var```, given some (optional) weights ```w```. The value of ```pct``` should be between 0 and 1. If supplied, the vector ```w``` should have the same size as ```var```. If not supplied, the program assums uniform weigthing.

The subroutine finds the percentile ```pct``` using a Bisection method.

The user can also supply a ```mask``` to compute the conditional percentile. The input ```mask``` is a logical array of the same size of ```var```.

**Dependencies**: [```error```](error.md)

[(back to index)](../index.md)

---

**Example**

```fortran
pc60 = percentile(var,0.60d0,mask = var.gt.0.0d0 .and. var.lt.5.0d0)
```

This computes the 60th percentile of ```var``` conditional on ```var``` being between 0 and 5.
