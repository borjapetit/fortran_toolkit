## denormalize

```fortran
subroutine denormalize(y,x,xmax,xmin)
  implicit none
  real(kind=8) , intent(in)  :: y     ! input: unbounded variable
  real(kind=8) , intent(in)  :: xmax  ! input: upper-bound of the bounded variable
  real(kind=8) , intent(in)  :: xmin  ! input: lower-bound of the bounded variable
  real(kind=8) , intent(out) :: x     ! output: bounded variable corresponding to y
```

This subroutine takes an unbounded varibale $\texttt{y}$ and applies the transformation

$$ \texttt{x} = \texttt{xmin} + \left(\frac{\exp(y)}{1+\exp(y)} \right)(\texttt{xmax}-\texttt{xmin}) $$

to return a bounded variable $\texttt{x}$, contrained to be between $\texttt{xmin}$ and $\texttt{xmax}$.

_Note_: This subroutine, along with $\texttt{normalize}$, is useful to use stadard optimization algorithms for contrained problems.

**Dependencies**: none

[(back to index)](index.md)

---

**Example**

Given an unbounded variable $\texttt{betau}$, we want to get the corresponding bounded variable in the unit interval:
```fortran
call denormalize ( betau , beta , 1.0d0 , 0.0d0 )

! if betau = 3  --> beta = 0.952
! if betau = -2 --> beta = 0.119
```
