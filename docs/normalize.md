### normalize & denormalize

```fortran
subroutine normalize(y,x,xmax,xmin)
  implicit none
  real(kind=8) , intent(in)  :: xmax,xmin,x
  real(kind=8) , intent(out) :: y
```

This subroutine takes a bounded varibale  ```x```, contrained to be between ```xmin``` and ```xmax```, and return an unbounded variable ```y``` using the following transformation;

$$y = \log\left(\frac{\texttt{x} - \texttt{xmin}}{\texttt{xmax} - \texttt{x}}\right)$$

```fortran
subroutine denormalize(y,x,xmax,xmin)
  implicit none
  real(kind=8) , intent(in)  :: xmax,xmin,y
  real(kind=8) , intent(out) :: x
```

This subroutine takes an unbounded varibale ```y``` and returns a bounded variable ```x```, contrained to be between ```xmin``` and ```xmax```, applying the following transformation

$$\texttt{x} = \texttt{xmin} + (\texttt{xmax} - \texttt{xmin})\cdot\left( \frac{\exp(\texttt{y})}{ 1 + \exp(\texttt{y})}\right)$$

**Note**: These subroutines are useful to solve contrained problems using standard uncontrained optimization algorithms.

**Dependencies**: none

[(back to index)](../index.md)

---

**Example**

We want to minimize a function ```func``` over a 2-dimensional vector ```x```. The variables are contrained to be between 10 and 0, and between 5 and -5. To find the minimum using the Nelder-Mead algorithm (see the [```simplex```](simplex.md) subroutine) we can write an auxiliry function that takes a 2-dimensional unbounded vector ```y```, maps ```y``` onto ```x```, and evaluate the funciton ```func```.

```fortran
program main

  use toolkit
  implicit none
  
  real(dp) :: x0(2),y0(2),f0
  real(dp) :: x1(2),y1(2),f1
  integer  :: numiter,exitcode

  ! variable x(1) should be between 10 and 0
  ! variable x(2) should be between 5 and -5
  
  ! transform the (constrained) initial point "x0" into the 
  ! (unconstrained) initial point "y"
  call normalize(y0(1),x0(1),10.0d0,0.0d0)
  call normalize(y0(2),x0(2),5.0d0,-50.0d0)

  ! use Nelder-Mead to minimize the function "unconstrained_func" 
  ! over "y", using y0 as initial point
  call simplex(unconstrained_func,y1,f1,numiter,exitcode,y0)
    
  ! transform variable "y1" (unconstrained) into "x1" (constrained)
  call denormalize(y1(1),x1(1),10.0d0,0.0d0)   ! x(1) is between 10 and 0
  call denormalize(y1(2),x1(2),5.0d0,-50.0d0)  ! x(2) is between 5 and -5
  
  ! these variable satisfy: unconstrained_func(y1) = func(x1)
  
  return
  contains
  
  ! the function (of "x") that we want to minimze
  function func(x) result(f)
    implicit none
    real(dp) :: x(:),f
    f =  ! some function with x(1) and (x(2)
    return
  end function func

  ! auxiliary function that takes an uncontrained vector "y"
  ! transforms it into the constrained vector "x" and get the 
  ! value of the function "func"
  function unconstrained_func(y) result(f)
    implicit none
    real(dp) :: y(:),f
    real(dp) :: x(2)
    
    ! get values of "x"
    call denormalize(y(1),x(1),10.0d0,0.0d0)
    call denormalize(y(2),x(2),5.0d0,-50.0d0)
    
    ! evaluate the function
    f = func(x)
    
    return
  end function unconstrained_func
  
end program main
```
