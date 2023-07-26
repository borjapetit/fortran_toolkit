### varmean

```fortran
function varmean(var,w,mask) result(meanvar)
  implicit none
  real(kind=8)            :: meanvar
  real(kind=8)            :: var(:)
  real(kind=8) , optional :: w(:)     
  logical      , optional :: mask(:)  
```

This function returns the mean of a variable ```var``` given some (optional) weigths ```w```. If supplied, the vector ```w``` should have the same size as ```var```. If not supplied, the program assums uniform weigthing.

The user can also supply a ```mask``` to compute the conditional mean. The input ```mask``` is a logical array of the same size of ```var```.

**Dependencies**: none

[(back to index)](../index.md)

---

**Example**

```fortran
mean = varmean(var,mask = var.gt.0.0d0 .and. var.lt.5.0d0)
```

This computes the mean of ```var``` conditional on ```var``` being between 0 and 5.
