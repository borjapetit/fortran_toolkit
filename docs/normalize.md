
# Fortran toolkit

###### Borja Petit

---

## ```normalize``` & ```denormalize```

The optimization subroutines included in this toolkit are design to handle unconstrained problems. However, the $\texttt{normalize}$ and $\texttt{denormalize}$ subroutines can transform a constrained optimization problem into an unconstrained one.

Imagine we want to maximize a function $f(x):\mathbb{R}^n \to \mathbb{R}^m$, where $1\leq m\leq n$, conditional on $x\in\mathbb{X}^n$:

$$\max_{x\in\mathbb{X}^n} \ f(x)$$

To solve this problem with standard optimization algorithm, we can define an intermediate function $g(z)$ such that $g(z):\mathbb{R}^n\to\mathbb{X}^n$, and rewrite our original problem as:

$$\max_{x\in\mathbb{X}^n} \ f(x) \ \equiv \ \max_{z} \ f(g(z)) \ \ \ \text{with} \ \ x = g(z)$$

The problem $\max_{z} \ f(g(z))$ is an unconstrained optimization problem that can be solve with the subroutines included in this toolkit.

**How to use**: we have a function ```func``` that takes a variable (scalar or vector) ```x``` that should be between $a$ and $b$. To minimize ```func``` we need to:
- Use the ```normalize``` subroutine to define a new variable ```z```
- Define an auxiliary function ```func0``` that takes ```z```, applies the ```denormalize``` subroutine (to get ``` x``` back) and calls ```func(x)```.
- Use any of the optimization subroutines in this toolkit (```simplex```,```lmmin```) to minimize ```func0``` with respect to ```z```
- Apply the ```denormalize``` subroutine to the optimal ```z``` to get the corresponding value of ```x```.

Download a complete example, explained below, from this [link](https://borjapetit.github.io/fortran_toolkit/src/test_normalize.f90).

[(back to index)](../index.md)

---

#### normalize

```fortran
subroutine normalize(z,x,xmax,xmin)
  implicit none
  real(kind=8) , intent(out) :: z     ! output: unbounded variable
  real(kind=8) , intent(in)  :: x     ! input: bounded variable
  real(kind=8) , intent(in)  :: xmax  ! input: upper-bound of x
  real(kind=8) , intent(in)  :: xmin  ! input: lower-bound of x
```

This subroutine, that corresponds to the inverse function $g^{-1}(x)$ explained above, takes a bounded variable $\texttt{x}$, constrained to be between $\texttt{xmin}$ and $\texttt{xmax}$, and return an unbounded variable $\texttt{z}$ using the following transformation:

$$z = \log\left(\frac{\texttt{x} - \texttt{xmin}}{\texttt{xmax} - \texttt{x}}\right)$$

**Dependencies**: none

[(back to index)](../index.md)

--- 

#### denormalize

```fortran
subroutine denormalize(z,x,xmax,xmin)
  implicit none
  real(kind=8) , intent(in)  :: z     ! input: unbounded variable
  real(kind=8) , intent(in)  :: xmax  ! input: upper-bound of the bounded variable
  real(kind=8) , intent(in)  :: xmin  ! input: lower-bound of the bounded variable
  real(kind=8) , intent(out) :: x     ! output: bounded variable corresponding to z
```

This subroutine, that corresponds to the function $g(z)$ explained above, takes an unbounded variable $\texttt{z}$ and applies the transformation:

$$ \texttt{x} = \texttt{xmin} + \left(\frac{\exp(z)}{1+\exp(z)} \right)\cdot(\texttt{xmax}-\texttt{xmin}) $$

to return a bounded variable $\texttt{x}$, constrained to be between $\texttt{xmin}$ and $\texttt{xmax}$.

**Dependencies**: none

[(back to index)](../index.md)

---

**Example**

Imagine we want to find the minimum in the Booth function:

$$f(x_1,x_2) = ( x_2 + 2x_1 - 7 )^2 + ( 2x_2 + x_2 - 5 )^2, \ \ \ \text{s.t.} \ \ -10\leq x_1,x_2\leq 10$$

Since this is a single-valued function with 2 unknowns, we use the the Nelder-Mead algorithm (see the [```simplex```](simplex.md) subroutine).
To use this algorithm we first write an auxiliary function, ```unconstrained_booth```, that takes a 2-dimensional unbounded vector $\texttt{z}$, maps $\texttt{z}$ onto $\texttt{x}$, and evaluate the function ```booth```.

```fortran
program main

    use toolkit , only : dp,normalize,denormalize,simplex
    implicit none
    
    real(dp) :: x0(2),z0(2),f0
    real(dp) :: x1(2),z1(2),f1
    integer  :: numiter,exitcode

    ! x(1) and x(2) should be between 10 and -10  

    ! Initial guess:
    x0(1) = 5.0d0
    x0(2) = -3.0d0
  
    ! transform the (constrained) initial point "x0" into the 
    ! (unconstrained) initial point "z0"
    call normalize(z0(1),x0(1),10.0d0,-10.0d0)
    call normalize(z0(2),x0(2),10.0d0,-10.0d0)
  
    ! use Nelder-Mead to minimize the function "unconstrained_booth" 
    ! over "z", using "z0" as initial point
    call simplex(unconstrained_booth,z1,f1,numiter,exitcode,z0)
      
    ! transform "z1" (unconstrained) into "x1" (constrained)
    call denormalize(z1(1),x1(1),10.0d0,-10.0d0)
    call denormalize(z1(2),x1(2),10.0d0,-10.0d0)
    
    ! these variable satisfy: unconstrained_booth(z1) = booth(x1)
  
    write(*,'(a,3(f10.4))') ' Actual solution:        ',1.0d0,3.0d0, booth((/1.0d0,3.0d0/))
    write(*,'(a,3(f10.4))') ' Unconstrained solution: ',z1(:), unconstrained_booth(z1)
    write(*,'(a,3(f10.4))') ' Constrained solution:   ',x1(:), booth(x1)
    
    return
    contains
    
    ! the function (of "x") that we want to minimize
    function booth(x) result(f)
      implicit none
      real(dp) :: x(:),f
      f =  ( x(1) + dble(2.0)*x(2) - dble(7.0) )**2.0d0 + &
           ( dble(2.0)*x(1) + x(2) - dble(5.0) )**2.0d0
      return
    end function booth

    ! auxiliary function that takes an unconstrained vector "y"
    ! transforms it into the constrained vector "x" and get the 
    ! value of the function "func"
    function unconstrained_booth(z) result(f)
      implicit none
      real(dp) :: f
      real(dp) :: z(:)
      real(dp) :: x(2)
      
      ! get values of "x"
      call denormalize(z(1),x(1),10.0d0,-10.0d0)
      call denormalize(z(2),x(2),10.0d0,-10.0d0)
      
      ! evaluate the function
      f = booth(x)
      
      return
    end function unconstrained_booth
    
  end program main
```
The output of this program is:

```
  Actual solution:            1.0000    3.0000    0.0000
  Unconstrained solution:     0.2007    0.6190    0.0000
  Constrained solution:       1.0001    2.9998    0.0000
```