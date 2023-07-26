---
title: varstd
exclude: true
---

```fortran
function varstd(var,w,mask) result(stdvar)
  implicit none
  real(kind=8)            :: stdvar
  real(kind=8)            :: var(:)
  real(kind=8) , optional :: w(:)   
  logical      , optional :: mask(:)
```

This function returns the stadard deviation of a variable ```var``` given some (optional) weigths ```w```. If supplied, the vector ```w``` should have the same size as ```var```. If not supplied, the program assums uniform weigthing.

The user can also supply a ```mask``` to compute the conditional standard deviation. The input ```mask``` is a logical array of the same size of ```var```.

**Dependencies**: [```error```](error.md), [```varmean```](varmean.md),  [```varvar```](varvar.md)

[(back to index)](../index.md)

---

**Example**

```fortran
standdesv = varstd(var,mask = var.gt.0.0d0 .and. var.lt.5.0d0)
```

This computes the standard deviation of ```var``` conditional on ```var``` being between 0 and 5.
