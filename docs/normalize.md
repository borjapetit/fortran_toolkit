## normalize

```fortran
subroutine normalize(y,x,xmax,xmin)
  implicit none
  real(kind=8) , intent(in)  :: xmax,xmin,x
  real(kind=8) , intent(out) :: y
```

_Dependencies_: none

This subroutine takes a bounded varibale  ```x```, contrained to be between ```xmin``` and ```xmax```, and applies the transformation 
$$ \texttt{y} = \log\left( \frac{ \texttt{x} - \texttt{xmin}}{\texttt{xmax} - \texttt{x}} \right) $$
to return an unbounded variable ```y```.

**Note**: This subroutine is useful to use stadard optimization algorithms for contrained problems. For isntance, one can use the Simplex method to minimize a function of ```x```, where ```x``` should be in the unit interval, by making use of ```normalize``` and ```denormalize```.

[(back to index)](inicio.md)