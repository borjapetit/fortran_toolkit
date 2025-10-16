
## ```brent```

```fortran
subroutine brent(func,x,numiter,exitcode,x0,x1,itermax,tol,iprint)
  external                             :: func      ! user-supplied function
  real(kind=8) , intent(out)           :: x         ! output: root
  integer      , intent(out)           :: numiter   ! output: number of functions evaluations
  integer      , intent(out)           :: exitcode  ! output: exit code
  real(kind=8) , intent(in)            :: x0        ! input: lower-bound of x
  real(kind=8) , intent(in)            :: x1        ! input: upper-bound of x
  real(kind=8) , intent(in) , optional :: tol       ! input: (optional) level of tolerance [default = 1.0d-8]
  integer      , intent(in) , optional :: itermax   ! input: (optional) maximum function evaluations [default = 500]
  integer      , intent(in) , optional :: iprint    ! input: (optional) control what's printed [default = 0]
```

Find the root of a single-valued univariate equation using the Brent's method (more info [here](https://en.wikipedia.org/wiki/Brent%27s_method)). The user must specify a function ```func```, and an two values of ```x```. If ```func```(```x0```) and ```func```(```x1```) have opposite signs the existence of a root is guaranteed (if not, the algorithm still runs). And if there are more than one root, the algorithm converges to one of them depending on the initial points.

The function ```func``` must be of the form:

```fortran
function func(x) result(f)
  implicit none
  real(kind=8) :: x,f
  f = ... some function of x ...
end function func
```

Optionally, the user can also supply a level of tolerance (```tol```) and maximum number of function evaluations (```itermax```). Finally, the user can also control what it is printing during execution by setting the corresponding value of ```iprint```:
- ```iprint``` = 0: don't print anything (default)
- ```iprint``` = 1: print warnings
- ```iprint``` = 2: print warnings and every iteration

The subroutine returns the value of ```x``` that makes ```func``` close enoughs to zero, the number of function evaluations (```numiter```), and an indicator, ```exitcode```:
- ```exitcode``` = 0: the algorithm found a root
- ```exitcode``` = 1: either ```x0``` or ```x1``` is a root of ```func```
- ```exitcode``` = 2: the root is not within the interval (```x0```, ```x1```)
- ```exitcode``` = 3: the points ```x0``` and ```x1``` are too close
- ```exitcode``` = 9: maximum number of function evaluations reached

**Dependencies**: none

[(back to index)](../index.md)

---

**Example**

The following example finds the root of the function

$$f(x) = x^2 - 3.5 x - 5$$

which has two roots at $x=-1.0895$ and $x=4.5895$.

```fortran
program example
    use toolkit , only : dp,brent
    implicit none
    real(dp) :: xvar
    integer  :: numiter
    integer  :: exitcode
      
    write(*,*) 'find the root: initial points 0 and 15, with f(0)<0 and f(15)>0'

    call brent(func,xvar,numiter,exitcode,0.0d0,15.0d0,iprint=2)

    write(*,*) 'find the root: initial points 5 and 10, with f(5)>0 and f(10)>0'

    call brent(func,xvar,numiter,exitcode,5.0d0,10.0d0,iprint=2)
  
    return
  
    contains
  
    function func(x) result(y)
      implicit none
      real(dp) :: x,y
      y = x**2.0 - 3.5d0*x - 5.0d0
      return
    end function func
  
  end program example
```

The output of this program is:

```
 find the root: initial points 0 and 15, with f(0)<0 and f(15)>0

 starting brent algorithm

   Iteration =    4 | xs =     0.4348 | f(xs) =    -6.3327
   Iteration =    6 | xs =    -1.5559 | f(xs) =     2.8668
   Iteration =    8 | xs =    -0.9356 | f(xs) =    -0.8501
   Iteration =   10 | xs =    -1.0775 | f(xs) =    -0.0679
   Iteration =   12 | xs =    -1.0895 | f(xs) =     0.0003
   Iteration =   14 | xs =    -1.0895 | f(xs) =    -0.0000
   Iteration =   16 | xs =    -1.0895 | f(xs) =    -0.0000

   solved: a root is found

   iter =   16
   x    =    -1.0895
   func =    -0.0000

 find the root: initial points 5 and 10, with f(5)>0 and f(10)>0

 starting brent algorithm

   root is not bracketed: existence is not guaranteed

   Iteration =    4 | xs =     4.7826 | f(xs) =     1.1342
   Iteration =    6 | xs =     4.6021 | f(xs) =     0.0718
   Iteration =    8 | xs =     4.5899 | f(xs) =     0.0024
   Iteration =   10 | xs =     4.5895 | f(xs) =     0.0000
   Iteration =   12 | xs =     4.5895 | f(xs) =     0.0000

   solved: a root is found

   iter =   12
   x    =     4.5895
   func =     0.0000
```
