
<span style="text-align:right;display:block;">
<a href="https://borjapetit.github.io/fortran_toolkit/">Back to index</a>
</span>

## ```simplex```

```fortran
subroutine simplex(func,x,y,numiter,exitcode,x0,itermax,tol,iprint,tolx,step,fzero,lb,ub)
  implicit none
  external                             :: func    
  real(kind=8) , intent(out)           :: x(:)    
  real(kind=8) , intent(out)           :: y       
  integer      , intent(out)           :: numiter   
  integer      , intent(out)           :: exitcode 
  real(kind=8) , intent(in)            :: x0(:)
  real(kind=8) , intent(in) , optional :: tol  
  integer      , intent(in) , optional :: itermax 
  integer      , intent(in) , optional :: iprint  
  real(kind=8) , intent(in) , optional :: tolx
  real(kind=8) , intent(in) , optional :: step
  logical      , intent(in) , optional :: fzero
  real(kind=8) , intent(in) , optional :: lb(:)
  real(kind=8) , intent(in) , optional :: ub(:)
```

This subroutine minimizes a function using the Nelder-Mead algorithm (click [here](https://en.wikipedia.org/wiki/Nelder–Mead_method) for more information), with the adaptive coefficients of Gao and Han (2012, *Computational Optimization and Applications* 51, 259–277), which perform better than the standard ones in high-dimensional problems: reflection 1, expansion 1+2/n, contraction 0.75−1/(2n), and shrink 1−1/n.

The user should input the function ```func``` and the initial guess ```x0```. The function ```func``` must be of the form:

```fortran
function func(x) result(f)
  implicit none
  real(kind=8) :: x(:),f
  f = ... some function of (x1,x2,...) ...
end function func
```

The subroutine returns the value(s) of ```x``` that minimize ```func```, the value of ```func``` at that point (```y```), the number of function evaluations (```numiter```), and an indicator, ```exitcode```:
- ```exitcode``` = 0: ```func``` is smaller than ```tol``` (only if ```fzero``` = .true.)
- ```exitcode``` = 1: the algorithm has converged: the spread of ```func``` in the simplex is smaller than ```tol``` and the size of the simplex is smaller than ```tolx``` (or the size of the simplex is smaller than ```tol```)
- ```exitcode``` = 9: maximum number of function evaluations reached

Optionally, the user can also supply:
- ```itermax``` (integer): max number of function evaluations, default = 500 (the check is done at the beginning of each iteration, so the final count can exceed ```itermax``` by up to ```n```+1 evaluations)
- ```tol``` (double precision): tolerance level, default = 1.0d-8
- ```tolx``` (double precision): tolerance on the size of the simplex, default = sqrt(```tol```)
- ```step``` (double precision): size of the initial simplex, default = 0.20
- ```fzero``` (logical): if .true. (default), the algorithm stops as soon as ```func``` is smaller than ```tol```, which is useful when ```func``` is a sum of squared residuals. Set ```fzero``` = .false. for general minimization problems (e.g., when ```func``` can be negative).
- ```lb```, ```ub``` (double precision vectors of dimension "n"): lower and upper bounds for ```x``` (box constraints, ```lb``` <= ```x``` <= ```ub```). Either or both can be supplied; use ```-huge(1.0d0)``` or ```huge(1.0d0)``` for an unknown without a lower or an upper bound. The bounds are imposed by an internal transformation: the simplex works on a free variable ```z```, and ```func``` is evaluated at ```x(z)```:
  - both bounds: ```x = lb + (ub-lb)*sin(z)**2```
  - lower bound only: ```x = lb + z**2```
  - upper bound only: ```x = ub - z**2```
  - no bounds: ```x = z```

  Hence ```func``` is never evaluated outside the box, the bounds can be attained, and the subroutine returns ```x``` in its original units. For bounded unknowns, ```step``` and ```tolx``` refer to ```z```, not to ```x```. If ```x0``` is outside the box, it is projected onto it.
- ```iprint``` (integer): indicator for printing behaviour
  - ```iprint``` = 0: don't print anything (default)
  - ```iprint``` = 1: print main results
  - ```iprint``` = 2: print main results and each iteration

**Note**: for constrained optimization problems, supply ```lb``` and/or ```ub``` (there is no need to transform the variables with [```normalize and denormalize```](normalize.md)).

**Example with bounds**

```fortran
  ! minimize func with x(1) >= 0 and 4 <= x(2) <= 7
  call simplex(func,x,y,numiter,exitcode,x0,lb=(/0.0d0,4.0d0/),ub=(/huge(1.0d0),7.0d0/),fzero=.false.)
```

**Internal dependencies**: none

**Updates**: November 2026. Corrected the expansion and contraction points (previously ```(1-g)*xm - g*xhigh``` and ```(1-b)*xm - b*xhigh```, instead of ```xm + g*(xr-xm)``` and ```xm + b*(xhigh-xm)```); added the outside contraction, the adaptive coefficients, the convergence criterion on the spread of ```func```, and the optional inputs ```tolx```, ```step``` and ```fzero```; optional bounds ```lb``` and ```ub```.

---

**Example**

```fortran
program example

  use toolkit , only : dp,simplex

  implicit none

  real(dp) :: y
  real(dp) :: x0(2)
  real(dp) :: x1(2)
  integer  :: errorcode
  integer  :: itercount

  ! set initial guess to x1 = 2 and x2 = 3
  x0 = (/ 2.0d0 , 3.0d0 /)

  ! find the minimum of the Matyas function
  call simplex(matyas,x1,y,itercount,errorcode,x0)

  ! find the minimum of the Matyas function
  !     set a tolerance of 0.5
  !     set a max number of function evaluations to 1000
  !     print every iteration
  call simplex(matyas,x1,y,itercount,errorcode,x0,itermax=1000,tol=0.5d0,iprint=2)

  ! minimize a function that can be negative: do not stop when func < tol
  call simplex(negmatyas,x1,y,itercount,errorcode,x0,fzero=.false.)

  return

  contains

  ! Matyas function
  function matyas(x) result(y)
    implicit none
    real(dp) :: x(:),y
    y = dble(0.26)*( x(1)*x(1) + x(2)*x(2)) - dble(0.48)*x(1)*x(2)
    return
  end function matyas

  ! Matyas function minus one (minimum = -1)
  function negmatyas(x) result(y)
    implicit none
    real(dp) :: x(:),y
    y = matyas(x) - dble(1.0)
    return
  end function negmatyas

end program example
```
