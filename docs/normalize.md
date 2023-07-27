<style TYPE="text/css">
code.has-jax {font: inherit; font-size: 100%; background: inherit; border: inherit;}
</style>
<script type="text/x-mathjax-config">
MathJax.Hub.Config({
    tex2jax: {
        inlineMath: [['$','$'], ['\\(','\\)']],
        skipTags: ['script', 'noscript', 'style', 'textarea', 'pre'] // removed 'code' entry
    }
});
MathJax.Hub.Queue(function() {
    var all = MathJax.Hub.getAllJax(), i;
    for(i = 0; i < all.length; i += 1) {
        all[i].SourceElement().parentNode.className += ' has-jax';
    }
});
</script>
<script type="text/javascript" src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.4/MathJax.js?config=TeX-AMS_HTML-full"></script>

### normalize & denormalize

The optimization subroutines included in this toolkit are design to handle unconstrained problems. However, the $\texttt{normalize}$ and $\texttt{denormalize}$ subroutines can transform a contrained optimization problem into an unconstrained one.

Imagine we want to maxmize a function $f(x):\mathbb{R}^n \to \mathbb{R}^m$, where $1\leq m\leq n$. Imagine we want to maximize the function $f$ contidional on $x\in\mathbb{X}^n$:

`$$\max_{x\in\mathbb{X}^n} \ f(x)$$`

None of the subroutines in the toolkit can solve this problem directly. However, we can define an intermediate function $g(z)$ such that $g(z):\mathbb{R}^n\to\mathbb{X}^n$. Using this function, we can rewrite out original problem as:

$$\max_{x\in(x0,x_1)} \ f(x) \ \equiv \ \max_{z} \ f(g(z)) \ \ \ \text{with} \ \ x = g(z)$$

The problem $\max_{z} \ f(g(z))$ is an unconstrained optimization problem that can be solve with the subrotuines included in this toolkit.

---

#### normalize

```fortran
subroutine normalize(y,x,xmax,xmin)
  implicit none
  real(kind=8) , intent(out) :: y     ! output: unbounded variable
  real(kind=8) , intent(in)  :: x     ! input: bounded variable
  real(kind=8) , intent(in)  :: xmax  ! input: upper-bound of x
  real(kind=8) , intent(in)  :: xmin  ! input: lower-bound of x
```

This subroutine, that corresponds to the inverse function $g^{-1}(z)$ explained above, takes a bounded varibale $\texttt{x}$, contrained to be between $\texttt{xmin}$ and $\texttt{xmax}$, and return an unbounded variable $\texttt{y}$ using the following transformation:

$$y = \log\left(\frac{\texttt{x} - \texttt{xmin}}{\texttt{xmax} - \texttt{x}}\right)$$

**Dependencies**: none

[(back to index)](../index.md)

--- 

#### denormalize

```fortran
subroutine denormalize(y,x,xmax,xmin)
  implicit none
  real(kind=8) , intent(in)  :: y     ! input: unbounded variable
  real(kind=8) , intent(in)  :: xmax  ! input: upper-bound of the bounded variable
  real(kind=8) , intent(in)  :: xmin  ! input: lower-bound of the bounded variable
  real(kind=8) , intent(out) :: x     ! output: bounded variable corresponding to y
```

This subroutine, that corresponds to the function $g(z)$ explained above, takes an unbounded varibale $\texttt{y}$ and applies the transformation:

$$ \texttt{x} = \texttt{xmin} + \left(\frac{\exp(y)}{1+\exp(y)} \right)\cdot(\texttt{xmax}-\texttt{xmin}) $$

to return a bounded variable $\texttt{x}$, contrained to be between $\texttt{xmin}$ and $\texttt{xmax}$.

**Dependencies**: none

[(back to index)](../index.md)

---

**Example**

We want to minimize a function $\texttt{func}$ over a 2-dimensional vector $\texttt{x}$. The variables are contrained to be between 10 and 0, and between 5 and -5. To find the minimum we want to use the Nelder-Mead algorithm (see the [```simplex```](simplex.md) subroutine) we can write an auxiliry function that takes a 2-dimensional unbounded vector $\texttt{y}$, maps $\texttt{y}$ onto $\texttt{x}$, and evaluate the funciton $\texttt{func}$.

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
