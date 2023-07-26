---
exclude: true
---

```fortran
subroutine golden(func,x,y,numiter,xmax,xmin,itermax,tol,iprint)
  implicit none
  external                             :: func      ! user-supplied function
  real(kind=8) , intent(out)           :: x         ! output: arg max of func
  real(kind=8) , intent(out)           :: y         ! output: max of func
  integer      , intent(out)           :: numiter   ! output: number of function evaluations
  real(kind=8) , intent(in)            :: xmax      ! input: upper-bound of x
  real(kind=8) , intent(in)            :: xmin      ! input: lower-bound of x
  real(kind=8) , intent(in) , optional :: tol       ! input: (optional) level of tolerance [detault = 1.0d-8]
  integer      , intent(in) , optional :: itermax   ! input: (optional) maximum function evaluations [default = 500]
  integer      , intent(in) , optional :: iprint    ! input: (optional) control what's printed [detault = 0]
```

This subroutine finds the maximum of a user-supplied single-valued function, $\texttt{func}$, with one unknown using the Golden Search algorithm. The function $\texttt{func}$ must be of the form:

```fortran
function func(x) result(f)
  implicit none
  real(kind=8) :: x,f
  f = ... some function of x ...
end function func
```

The idea of this algorithm is to keep shrinking an interval $(x_1,x_2)$ inside which the optimal $x$ lays. The user needs to provide the initial interval $(\texttt{xmin},\texttt{xmax})$.

Optionally, the user can also supply a maximun number of function evaluations ($\texttt{itermax}$, 500 by default), the level of tolerance ($\texttt{tol}$, 1.0d-8 by default). Finally, the user can also control what it is printed during execution by setting the corresponding value of $\texttt{iprint}$:

- $\texttt{iprint}$ = 0: don't print anything (default)
- $\texttt{iprint}$ = 1: print warnings
- $\texttt{iprint}$ = 2: print warnings and every iteration

The subroutine returns the value of $\texttt{x}$ that maximizes $\texttt{func}$, the value of $\texttt{func}$ at $\texttt{x}$ ($\texttt{y}$),  and the number of function evaluations ($\texttt{numiter}$)

**Dependencies**: [```error```](error.md)

[(back to index)](../index.md)

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

  real(dp) :: xvar,f
  real(dp) :: xmax,xmin

  ! max and min values for x
  xmax = 10.0 ; xmin = 0.0

  ! maximize the function func
  call golden(func,xvar,f,xmax,xmin,iprint=2)

  return

  contains

  function func(x) result(y)
    implicit none
    real(dp) :: x,y
    y = x**2.0d0 - 3.5d0*x
    return
  end function func

end program example
```

The output is:

```
 starting golden search algorithm

 Iteration =   19  | x1 =     2.3607  x2 =     3.8197  | f(x1) =    10.5786  f(x2) =    11.1292
 Iteration =   20  | x1 =     3.8197  x2 =     5.2786  | f(x1) =    11.1292  f(x2) =    11.1937
 Iteration =   21  | x1 =     5.2786  x2 =     6.1803  | f(x1) =    11.1937  f(x2) =    11.0901
 Iteration =   22  | x1 =     4.7214  x2 =     5.2786  | f(x1) =    11.2088  f(x2) =    11.1937
 Iteration =   23  | x1 =     4.3769  x2 =     4.7214  | f(x1) =    11.1953  f(x2) =    11.2088
 Iteration =   24  | x1 =     4.7214  x2 =     4.9342  | f(x1) =    11.2088  f(x2) =    11.2081
 Iteration =   25  | x1 =     4.5898  x2 =     4.7214  | f(x1) =    11.2059  f(x2) =    11.2088
 Iteration =   26  | x1 =     4.7214  x2 =     4.8027  | f(x1) =    11.2088  f(x2) =    11.2093
 Iteration =   27  | x1 =     4.8027  x2 =     4.8529  | f(x1) =    11.2093  f(x2) =    11.2091
 Iteration =   28  | x1 =     4.7716  x2 =     4.8027  | f(x1) =    11.2092  f(x2) =    11.2093
 Iteration =   29  | x1 =     4.8027  x2 =     4.8219  | f(x1) =    11.2093  f(x2) =    11.2093
 Iteration =   30  | x1 =     4.7908  x2 =     4.8027  | f(x1) =    11.2093  f(x2) =    11.2093
 Iteration =   31  | x1 =     4.8027  x2 =     4.8100  | f(x1) =    11.2093  f(x2) =    11.2093
 Iteration =   32  | x1 =     4.7981  x2 =     4.8027  | f(x1) =    11.2093  f(x2) =    11.2093
 Iteration =   33  | x1 =     4.8027  x2 =     4.8055  | f(x1) =    11.2093  f(x2) =    11.2093
 Iteration =   34  | x1 =     4.8009  x2 =     4.8027  | f(x1) =    11.2093  f(x2) =    11.2093
 Iteration =   35  | x1 =     4.8027  x2 =     4.8037  | f(x1) =    11.2093  f(x2) =    11.2093
 Iteration =   36  | x1 =     4.8037  x2 =     4.8044  | f(x1) =    11.2093  f(x2) =    11.2093
 Iteration =   37  | x1 =     4.8033  x2 =     4.8037  | f(x1) =    11.2093  f(x2) =    11.2093
 Iteration =   38  | x1 =     4.8037  x2 =     4.8040  | f(x1) =    11.2093  f(x2) =    11.2093
 Iteration =   39  | x1 =     4.8040  x2 =     4.8041  | f(x1) =    11.2093  f(x2) =    11.2093
 Iteration =   40  | x1 =     4.8039  x2 =     4.8040  | f(x1) =    11.2093  f(x2) =    11.2093
 Iteration =   41  | x1 =     4.8040  x2 =     4.8040  | f(x1) =    11.2093  f(x2) =    11.2093
 Iteration =   42  | x1 =     4.8040  x2 =     4.8040  | f(x1) =    11.2093  f(x2) =    11.2093
 Iteration =   43  | x1 =     4.8040  x2 =     4.8040  | f(x1) =    11.2093  f(x2) =    11.2093
 Iteration =   44  | x1 =     4.8040  x2 =     4.8040  | f(x1) =    11.2093  f(x2) =    11.2093
 Iteration =   45  | x1 =     4.8040  x2 =     4.8040  | f(x1) =    11.2093  f(x2) =    11.2093
 Iteration =   46  | x1 =     4.8040  x2 =     4.8040  | f(x1) =    11.2093  f(x2) =    11.2093
 Iteration =   47  | x1 =     4.8040  x2 =     4.8040  | f(x1) =    11.2093  f(x2) =    11.2093
 Iteration =   48  | x1 =     4.8040  x2 =     4.8040  | f(x1) =    11.2093  f(x2) =    11.2093
 Iteration =   49  | x1 =     4.8040  x2 =     4.8040  | f(x1) =    11.2093  f(x2) =    11.2093
 Iteration =   50  | x1 =     4.8040  x2 =     4.8040  | f(x1) =    11.2093  f(x2) =    11.2093
 Iteration =   51  | x1 =     4.8040  x2 =     4.8040  | f(x1) =    11.2093  f(x2) =    11.2093
 Iteration =   52  | x1 =     4.8040  x2 =     4.8040  | f(x1) =    11.2093  f(x2) =    11.2093
 Iteration =   53  | x1 =     4.8040  x2 =     4.8040  | f(x1) =    11.2093  f(x2) =    11.2093
 Iteration =   54  | x1 =     4.8040  x2 =     4.8040  | f(x1) =    11.2093  f(x2) =    11.2093
 Iteration =   55  | x1 =     4.8040  x2 =     4.8040  | f(x1) =    11.2093  f(x2) =    11.2093
 Iteration =   56  | x1 =     4.8040  x2 =     4.8040  | f(x1) =    11.2093  f(x2) =    11.2093
 Iteration =   57  | x1 =     4.8040  x2 =     4.8040  | f(x1) =    11.2093  f(x2) =    11.2093
 Iteration =   58  | x1 =     4.8040  x2 =     4.8040  | f(x1) =    11.2093  f(x2) =    11.2093
 
 Solved: Iterations =   58  | x =     4.8040  | f(x) =    11.2093
```





