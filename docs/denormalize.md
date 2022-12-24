# denormalize

```fortran
subroutine denormalize(y,x,xmax,xmin)
  implicit none
  real(kind=8) , intent(in)  :: xmax,xmin,y
  real(kind=8) , intent(out) :: x
```

This subroutine takes an unbounded varibale ```y``` and applies the transformation
$$ \texttt{x} = \texttt{xmin} + \left(\frac{\exp(y)}{1+\exp(y)} \right)(\texttt{xmax}-\texttt{xmin}) $$
to return a bounded variable ```x```, contrained to be between ```xmin``` and ```xmax```.

```fortran
call denormalize ( betau , beta , 1.0d0 , 0.0d0 )

! if betau = 3 --> beta = 0.952
! if betau = -2 --> beta = 0.119
```

**Note**: This subroutine is useful to use stadard optimization algorithms for contrained problems. For isntance, one can use the Simplex method to minimize a function of ```x```, where ```x``` should be in the unit interval, by making use of ```normalize``` and ```denormalize```.

**Dependencies**: none

[(back to index)](index.md)
