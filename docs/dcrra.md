
<span style="text-align:right;display:block;">
<a href="https://borjapetit.github.io/fortran_toolkit/">Back to index</a>
</span>

## ```dcrra```

```fortran
elemental function dcrra(x,b,tolfun) result(u)
  implicit none
  real(kind=8) , intent(in)            :: x,b
  real(kind=8) , intent(in) , optional :: tolfun
  real(kind=8)                         :: u
```

This function returns the derivative of the [```crra```](crra.md) function with respect to ```x``` (marginal utility): ```x**(b-1)``` if ```b``` is different from zero, and ```1/x``` (derivative of ```log(x)```) if ```b``` = 0. Below ```tolfun``` (default = ```tolvl```) the derivative is extended linearly, so that it is defined for any ```x```.

The function is elemental, so ```x``` and ```b``` can be arrays.

**Internal dependencies**: none

---
