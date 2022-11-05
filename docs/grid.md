### grid

```fortran
function grid(maxv,minv,n,s) result(v)
  implicit none
  real(kind=8)            :: maxv,minv,v(n)
  real(kind=8) , optional :: s
  integer                 :: n
```

_Dependencies_: none

This function creates a grid of ```n``` points between ```maxv``` and ```minv``` with a curvature of ```s```

- if ```s```$=1$: linear grid _(default)_
- if ```s```$>1$: more grids points around ```maxv```
- if ```s```$<1$: more grids points around ```min```

_Example_ 1: create a linear grid with 100 points between 0 and 10:

```fortran
vector = grid( 10.d0 , 0.0d0 , 100 )  
```

_Example_ 2: create a quadratic grid with 500 points between -1 and 1:

```fortran
vector = grid( 1.d0 , -1.0d0 , 400 , 2.0d0 )  
```

[(back to index)](index.md)
