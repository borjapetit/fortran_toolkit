### correlation

```fortran
function correlation(xvar1,xvar2,w,mask) result(corr)
  implicit none
  real(kind=8)            :: corr
  real(kind=8)            :: xvar1(:),xvar2(:)  ! both vectors should have the same length
  real(kind=8) , optional :: w(:)               ! same length of "xvar1" and "xvar2"
  logical      , optional :: mask               ! same length of "xvar1" and "xvar2"
```

This function returns the correlation coefficient between two (same size) vectors ```xvar1``` and ```xvar2``` given some (optional) weigths ```w```. If supplied, the vector ```w``` should have the same size as ```xvar1``` and ```xvar2```. If not supplied, the program assumes uniform weigthing.

The user can also supply a ```mask``` to compute the conditional correlation. The variable ```mask``` is a logical array of the same size of ```xvar1``` and ```xvar2```. For example:

```fortran
corr = correlation(xvar1,xvar2,mask = xvar1.gt.0.0d0 .and. xvar2.lt.5.0d0)
```

This computes the correlation between vectors ```xvar1``` and ```xvar2``` only for those pair of points that satify the condition ```xvar1.gt.0.0d0``` and ```xvar2.lt.5.0d0```.

**Dependencies**: [```varmean```](varmean.md),  [```varvar```](varvar.md)

[(back to index)](../index.md)
