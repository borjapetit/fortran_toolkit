
## grid

```fortran
function grid(maxv,minv,n,s) result(v)
  implicit none
  real(kind=8)            :: v(n)  ! output: generated grid
  real(kind=8)            :: maxv  ! input: upper-bound of x
  real(kind=8)            :: minv  ! input: lower-bound of x
  integer                 :: n     ! input: number of points in teh grid
  real(kind=8) , optional :: s     ! input: (optional) curvatura parameter 
```

This function creates a grid of $\texttt{n}$ points between $\texttt{maxv}$ and $\texttt{minv}$ with a curvature of $\texttt{s}$:

- if $\texttt{s}=1$ or missing: linear grid (default)
- if $\texttt{s}>1$: more grids points around $\texttt{maxv}$
- if $\texttt{s}<1$: more grids points around $\texttt{minv}$

**Dependencies**: [`error`](error.md)

[(back to index)](../index.md)

---

**Example**

```fortran
! create a linear grid with 100 points between 0 and 10
vector = grid( 10.d0 , 0.0d0 , 100 )  

! create a quadratic grid with 500 points between -1 and 1
vector = grid( 1.d0 , -1.0d0 , 500 , 2.0d0 )  
```


