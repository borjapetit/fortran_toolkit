
## ```ces```

```fortran
function ces(x1,x2,a,b) result(u)
  implicit none
  real(kind=8) :: x1,x2,a,b
  real(kind=8) :: u
```

This function aggregates ```x1``` and ```x2```with a CES function:
$$ u = 
\begin{cases}
  \left( a\,x_1^b + (1-a)\,x_2^b\right)^{1/b} & \text{ if } b>0 \\ \\ 
  x_1^a \,x_2^{1-a} & \text{ if } b\approx 0
\end{cases}$$

**Dependencies**: none

[(back to index)](../index.md)
