
<span style="text-align:right;display:block;">
<a href="https://borjapetit.github.io/fortran_toolkit/">Back to index</a>
</span>

## ```golden```

```fortran
subroutine golden(func,x,y,xmax,xmin,itermax,tol,iprint)
  implicit none
  external                             :: func      ! user-supplied function
  real(kind=8) , intent(out)           :: x         ! output: arg max of func
  real(kind=8) , intent(out)           :: y         ! output: max of func
  real(kind=8) , intent(in)            :: xmax      ! input: upper-bound of x
  real(kind=8) , intent(in)            :: xmin      ! input: lower-bound of x
  real(kind=8) , intent(in) , optional :: tol       ! input: (optional) level of tolerance [default = 1.0d-8]
  integer      , intent(in) , optional :: itermax   ! input: (optional) maximum function evaluations [default = 500]
  integer      , intent(in) , optional :: iprint    ! input: (optional) control what's printed [default = 0]
```

This subroutine finds the maximum of a user-supplied single-valued function, `func`, with one unknown using the Golden Search algorithm. The function `func` must be of the form:

```fortran
function func(x) result(f)
  implicit none
  real(kind=8) :: x,f
  f = ... some function of x ...
end function func
```

The idea of this algorithm is to keep shrinking an interval inside which the optimal $x$ lies, keeping two interior points placed at the golden ratio, so that each iteration requires a single function evaluation and reduces the interval by a factor of 0.618. The user needs to provide the initial interval (`xmin`,`xmax`). The function should be unimodal in that interval; otherwise the algorithm returns a local maximum.

Optionally, the user can also supply a maximum number of iterations (`itermax`, 500 by default) and the level of tolerance (`tol`, 1.0d-8 by default), relative to the size of the interior points: the algorithm stops when `abs(x3-x0)` $\leq$ `tol*(abs(x1)+abs(x2))`. Finally, the user can also control what it is printed during execution by setting the corresponding value of `iprint`:

- `iprint` = 0: don't print anything (default)
- `iprint` = 1: print the final result
- `iprint` = 2: print the final result and every iteration

The subroutine returns the value of `x` that maximizes `func` and the value of `func` at `x` (`y`).

**Internal dependencies**: none

---

**Example**

We want to maximize the function

$$f(x) = 10 x^{0.3} - x$$

This function is maximized at $x=4.8040$, with $f(4.8040) = 11.2093$.

The following code maximizes the function $f(x)$ using the Golden Search algorithm:

```fortran
program example

  use toolkit , only : dp,golden

  implicit none

  real(dp) :: xvar,fx
  real(dp) :: xmax,xmin

  ! max and min values for x
  xmax = 10.0d0 ; xmin = 0.0d0

  ! maximize the function func
  call golden(func,xvar,fx,xmax,xmin,iprint=2)

  return

  contains

  function func(x) result(y)
    implicit none
    real(dp) :: x,y
    y = 10.0d0*x**0.3d0 - x
    return
  end function func

end program example
```

The output is:

```
  starting golden search algorithm 
                                   
  iteration =    1 | x0 =     0.0000 | x1 =     2.3607 | x2 =     3.8197 | x3 =     6.1803 | f(x1) =    10.5786 | f(x2) =    11.1292
  iteration =    2 | x0 =     2.3607 | x1 =     3.8197 | x2 =     4.7214 | x3 =     6.1803 | f(x1) =    11.1292 | f(x2) =    11.2088
  iteration =    3 | x0 =     3.8197 | x1 =     4.7214 | x2 =     5.2786 | x3 =     6.1803 | f(x1) =    11.2088 | f(x2) =    11.1937
  iteration =    4 | x0 =     3.8197 | x1 =     4.3769 | x2 =     4.7214 | x3 =     5.2786 | f(x1) =    11.1953 | f(x2) =    11.2088
  iteration =    5 | x0 =     4.3769 | x1 =     4.7214 | x2 =     4.9342 | x3 =     5.2786 | f(x1) =    11.2088 | f(x2) =    11.2081
  iteration =    6 | x0 =     4.3769 | x1 =     4.5898 | x2 =     4.7214 | x3 =     4.9342 | f(x1) =    11.2059 | f(x2) =    11.2088
  ...
  iteration =   37 | x0 =     4.8040 | x1 =     4.8040 | x2 =     4.8040 | x3 =     4.8040 | f(x1) =    11.2093 | f(x2) =    11.2093
  iteration =   38 | x0 =     4.8040 | x1 =     4.8040 | x2 =     4.8040 | x3 =     4.8040 | f(x1) =    11.2093 | f(x2) =    11.2093
  iteration =   39 | x0 =     4.8040 | x1 =     4.8040 | x2 =     4.8040 | x3 =     4.8040 | f(x1) =    11.2093 | f(x2) =    11.2093
    
  solved: iterations =   39 | x =     4.8040 | f(x) =    11.2093
```





