
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
  \dfrac{x^b}{b} & \text{ if } b\neq 0 \\ \\ 
  \log{x} & \text{ if } b\approx 0
\end{cases}
$$

where ```x``` and ```b``` are supplied by the user. Note that ```b``` is not the coefficient of relative risk aversion, $\sigma$, but $b = 1-\sigma$: the usual form $x^{1-\sigma}/(1-\sigma)$ is obtained with ```b``` $=1-\sigma$, and $\log x$ (the case $\sigma=1$) with ```b``` $=0$.

If ```x``` is (close to) zero, the function returns 0 when ```b``` $>0$, and ```-huge(x)``` when ```b``` $<0$ or ```b``` $\approx 0$.

**Internal dependencies**: none

