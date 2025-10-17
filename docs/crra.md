
<span style="text-align:right;display:block;">
<a href="https://borjapetit.github.io/fortran_toolkit/">Back to index</a>
</span>

## ```crra```

```fortran
function crra(x,b) result(u)
  implicit none
  real(kind=8) :: x,b
  real(kind=8) :: u
```

This function evaluates a CRRA utility function:

$$ u = 
\begin{cases}
  \dfrac{x^b}{b} & \text{ if } b>0 \\ \\ 
  \log{x} & \text{ if } b\approx 0
\end{cases}
$$

where ```x``` and ```b```are supplied by the user.

**Internal dependencies**: none

