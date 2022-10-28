# A Toolkit for Fortran90
<a name="inicio"></a>

General propuse:

- [```grid```](#grid): generate a grid for a continuous variable.
- [```interpolation```](#interpolation): interpolate a value over a grid, returning position and distance.
- [```interpolate```](#interpolate): linearly interpolate a value over an n-dimensional grid, with n <= 6.
- [```timing```](#timing): returns the number of seconds since 00:00h of the 1st day of the month (_robust to parelalization_).
- [```multiplo```](#multiplo): returns 1 if an integer is a multiple of another user-provided integer.
- [```iseven```](#iseven): returns 1 if a user-provided integer is even.
- [```error```](#error): print error message and interrupt execution.

Statistics:

- [```varmean```](#varmean): returns the average of a variable, allowing for weigths
- [```varstd```](#varstd): returns the standard deviation of a variable, allowing for weigths
- [```correlation```](#correlation): returns the correlation of two variables, allowing for weigths
- [```percentile```](#percentile): returns the i-th percentile of a variables, allowing for weigths
- [```olsreg```](#olsreg): returns the ols coefficients of a n-var regression (with n<=8), allowing for weigths
- [```tauchen```](#tauchen): returns the transition matrix for a discretized ar(1) process
- [```normaldist```](#normaldist): returns the distribution for a nomral random variable
- [```randomnormal```](#randomnormal): returns a random draw for a nomal distribution
- [```cdfn```](#cdfn): retutns the cdf of a nomabl distribution.

Linear algebra:

- [```vect```](#vect): transform a matrix of NxM into a vector of NxM rows
- [```cumsum```](#cumsum): returns the vector with cummulative sum of a vector (as Matlab's cumsum function)
- [```diag```](#diag): returns the main diagonal of a matrix
- [```transmat```](#transmat): returns the transpose of a square matrix
- [```inverse```](#inverse): returns the invesrse of a sqaured matrix

Optimization

- [```simplex```](#simplex): Simplex algorithm
- [```lmmin```](#lmmin): Levenberg–Marquardt algorithm
- [```golden```](#golden): Golden search algorithm
- [```brent```](#brent): Brent method
- [```normalize```](#normalize): transform a bounded variable into an unbounded one [for optimizaton]
- [```denormalize```](#denormalize): transform a undonded variable into an counded one [for optimizaton]
- [```broyden```](#broyden): updates a Jacobian matrix using the Broyden's method

## General purpose

### grid
<a name="grid"></a>

```fortran
function grid(maxv,minv,n,s) result(v)
  implicit none
  real(kind=8)            :: maxv,minv,v(n)
  real(kind=8) , optional :: s
  integer                 :: n
```

_Dependencies_: none

This function creates a grid of ```n``` points between ```maxv``` and ```minv``` with a curvature of ```s```

- if ```s```$=1$: linear grid _(default)_
- if ```s```$>1$: more grids points around ```maxv```
- if ```s```$<1$: more grids points around ```min```

_Example_ 1: create a linear grid with 100 points between 0 and 10:

```fortran
vector = grid( 10.d0 , 0.0d0 , 100 )  
```

_Example_ 2: create a quadratic grid with 500 points between -1 and 1:

```fortran
vector = grid( 1.d0 , -1.0d0 , 400 , 2.0d0 )  
```

[(back to index)](#inicio)

### interpolation
<a name="interpolation"></a>

```fortran
subroutine interpolation(pos,wth,xnow,xgrid)
  implicit none
  real(kind=8), intent(in)  :: xnow,xgrid(:)
  real(kind=8), intent(out) :: wth
  integer     , intent(out) :: pos
```

_Dependencies_: none

This subroutine finds the closest point in a given grid and return its position and the relative distance.

_Example_: consider a vector ```vec = [1,2,3]```. We can use ```interpolation``` to find the closest point of ```xnow = 2.3```:

```fortran
call interpolation(pos,wth,2.3,vec)

! Results: pos = 3, wth = 0.3
! 2.3 = vec(pos)*wth + vec(pos-1)*(one-wth) = 3*0.3 + 2*0.7
```

This subroutine is mainly used by the function ```interpolate```.

[(back to index)](#inicio)


### interpolate
<a name="interpolate"></a>

```fortran
function interpolate(x1,x2,...,xn,y1,y2,...,yn,mat) result(xi)
  implicit none
  real(kind=8) :: xi
  real(kind=8) :: x1,x2,...,xn
  real(kind=8) :: y1(:),y2(:),...,yn(:)
  real(kind=8) :: mat(:,:,...,:)
```

_Dependencies_: ```interpolation```

This function returns the linearly interpolated value of an n-dimensional function. The variables ```x1```, ```x2```, ..., ```xn``` are the values of the variables to be interpolated over their coresponding grids ```y1```, ```y2```, ...., ```yn```, and ```mat``` is an n-dimensional array with the results.

- The array ```mat``` must have at most, dimension 6 (so ```n```$\leq6$).
- If ```xn```$<\min$(```yn```) the subroutine takes $\min($```yn```$)$ as the value of ```x```.

_Example 1_: We have a 2-dimensional array ```mat``` whose $(i,j)$-element is the value of some function $f$ evaluated at the $i$-element of the vector $x$, and the $j$-element of the vector $y$. Then, we can interpolate the value of $f$ at $(x_0,y_0)$ using the ```interpolate```:

```fortran
xi = interpolate(x_0,y_0,x,y,mat)
```

_Example 2_: We have a 3-dimensional array ```mat``` whose $(i,j,k)$-element is the value of some function $f$ evaluated at the $i$-element of the vector $x$, the $j$-element of the vector $y$, and the $k$-element of the vecor $z$. Then, we can interpolate the value of $f$ at $(x_0,y_0,z_0)$ using the ```interpolate```:

```fortran
xi = interpolate(x_0,y_0,z_0,x,y,z,mat)
```

[(back to index)](#inicio)


### timing
<a name="timing"></a>

```fortran
function timing(mode) result(time)
  implicit none
  integer , optional :: mode
  real(kind=8)       :: time
```

_Dependencies_: ```none```

This functions returns a timing number that is robust to parallel computing. In particular, it returns the number of seconds since 00:00h of the 1st day of the month. The variable ```mode``` controls how time is measured:

- If ```mode``` = 1, time is measured in seconds (default).
- If ```mode``` = 2, time is measures in minutes.
- If ```mode``` = 3, time is measured in hours


[(back to index)](#inicio)


### multiplo
<a name="multiplo"></a>

```fortran
elemental function multiplo(num,xx) result(mul)
  implicit none
  integer :: num,xx
  logical :: mul
```

_Dependencies_: ```none```

This function checks whether a number ```num0``` is a multiple of ```num1```. The result ```mul``` is a logical variable taking value ```.TRUE.``` if ```num0``` is a multiple of ```num1```, and ```.FALSE.``` otherwise.

_Example_: check whether 25 and 27 are multiples of 5:

```fortran
write(*,*) multiplo(5,25)  ! .TRUE.
write(*,*) multiplo(5,27)  ! .FALSE.
```

[(back to index)](#inicio)


### iseven
<a name="iseven"></a>

```fortran
function iseven(num) result(ise)
  implicit none
  integer :: num
  logical :: ise
```

_Dependencies_: none

This function checks whether a number ```num``` is even. The result ```ise``` is a logical variable taking value ```.TRUE.``` if ```num``` is even, and ```.FALSE.``` otherwise.

_Example_: check whether 3 and 6 are eve:

```fortran
write(*,*) iseven(3)  ! .FALSE.
write(*,*) iseven(6)  ! .TRUE.
```

[(back to index)](#inicio)

### error
<a name="error"></a>

```fortran
subroutine error(mess)
  implicit none
  character(len=*) , intent(in) :: mess
```

_Dependencies_: none

This subroutine prints an error message ```mess``` and interrupt the execution of the program until the user type an interger.

  
  

## Statistics

### varmean
<a name="varmean"></a>

```fortran
function varmean(var,w,mask) result(meanvar)
  implicit none
  real(kind=8)            :: var(:)
  real(kind=8) , optional :: w(:)
  real(kind=8)            :: meanvar
  logical      , optimal  :: mask
```

_Dependencies_: none

This function returns the mean of a variable ```var``` given some (optional) weigths ```w```. The user can also supply as ```mask``` to compute the conditional mean. If supplied, the vector ```w``` should have the same size as ```var```. If not supplied, the program assums uniform weigthing.

_Example_: compute the mean of a vector

```fortran
! without weigths
xvar = (/ 1.0, 4.0, 4.0, 9.0 /)
mean = varmean(xvar)  ! mean = 4.5

! without weigths and conditional on xvar>0
xvar = (/ 1.0, 4.0, 4.0, 9.0 /)
mean = varmean(xvar, mask = xvar.gt.2.0d0 )  ! mean = 5.66

! with weigths
xvar = (/ 1.0, 4.0, 4.0, 9.0 /)
wvar = (/ 2.0, 4.0, 5.0, 2.0 /)
mean = varmean(xvar,w = wvar)  ! mean = 4.3076
```

[(back to index)](#inicio)


### varstd
<a name="varstd"></a>

```fortran
function varstd(var,w,mask) result(stdvar)
  implicit none
  real(kind=8)            :: var(:)
  real(kind=8) , optional :: w(:)
  real(kind=8)            :: stdvar
  logical      , optimal  :: mask
```

_Dependencies_: ```varmean```

This function returns the stadard deviation of a variable ```var``` given some (optional) weigths ```w```. The user can also supply as ```mask``` to compute the conditional stadard deviation. If supplied, the vector ```w``` should have the same size as ```var```. If not supplied, the program assums uniform weigthing.

[(back to index)](#inicio)


### correlation
<a name="correlation"></a>

```fortran
function correlation(xvar1,xvar2,w,mask) result(corr)
  implicit none
  real(kind=8) , optional :: w(:)
  real(kind=8)            :: xvar1(:),xvar2(:)
  real(kind=8)            :: corr
  logical      , optimal  :: mask
```

_Dependencies_: ```varmean```, ```varstd```

This function returns the correlation coefficient between two variables ```xvar1``` and ```xvar2``` given some (optional) weigths ```w```. The user can also supply as ```mask``` to compute the conditional stadard deviation. If supplied, the vector ```w``` should have the same size as ```var```. If not supplied, the program assums uniform weigthing.

[(back to index)](#inicio)


### percentile
<a name="percentile"></a>

```fortran
function percentile(xvec,pct,w,mask) result(cutoff)
  implicit none
  real(kind=8)            :: xvec(:),pct
  real(kind=8) , optional :: w(:)
  real(kind=8)            :: cutoff
  logical      , optimal  :: mask
```

_Dependencies_: none

This function returns the percentile ```pct``` for a distribution ```xvec```, given some (optional) weigths ```w```. The user can also supply as ```mask``` to compute the conditional correlation. If supplied, the vector ```w``` should have the same size as ```var```. If not supplied, the program assums uniform weigthing.

_Example_: given a vector ```xvec``` with a sample of a variable ```x```, find the 60th percentile:

```fortran
pc60 = percentile(xvec,60.0d0)
```

[(back to index)](#inicio)


### olsreg
<a name="olsreg"></a>


```fortran
subroutine olsreg(coeffs,yvec,x1vec,x2vec,...,x8vec,w,mask,iprint)
  implicit none
  real(kind=8) , intent(out)           :: coeffs(:)
  real(kind=8) , intent(in)            :: yvec(:)
  real(kind=8) , intent(in)            :: x1vec(:)
  real(kind=8) , intent(in) , optional :: x2vec(:),x3vec(:),...,x8vec(:)
  real(kind=8) , intent(in) , optional :: w(:)
  logical      , intent(in) , optional :: mask(:)  
  integer      , intent(in) , optional :: iprint
```

_Dependencies_: ```varmean```, ```varvar```



[(back to index)](#inicio)


### tauchen
<a name="tauchen"></a>

```fortran
subroutine tauchen(xvec,rho,mu,sigma,n,pmat)
  implicit none
  integer     , intent(in)  :: n
  real(kind=8), intent(in)  :: rho,mu,sigma,xvec(n)
  real(kind=8), intent(out) :: pmat(n,n)
```

_Dependencies_: ```normaldist```

This function returns the transition matrix for a discretized AR(1) process of the form:
$$x' = \mu + \rho x + \sigma \epsilon, \hspace{0.2cm} \epsilon \sim N(0,1)$$
The vector with values of $x$, ```xvec```, if of dimension ```n``` and does not need to be equally spaced.

[(back to index)](#inicio)


### normaldist
<a name="normaldist"></a>

```fortran
subroutine normaldist(xvec,mu,sigma,n,dist)
  implicit none
  integer     , intent(in)  :: n
  real(kind=8), intent(in)  :: mu,sigma,xvec(n)
  real(kind=8), intent(out) :: dist(n)
```

_Dependencies_: ```cdfn```

This function returns the distribution of a normal random variable with mean ```mu``` and standadrd deviation ```sigma```.


[(back to index)](#inicio)


### randomnormal
<a name="randomnormal"></a>

```fortran
subroutine randomnormal(shock,mu,std)
  implicit none
  real(kind=8) , intent(in)  :: mu,std
  real(kind=8) , intent(out) :: shock ! or shock(:)
```

_Dependencies_: none

This function returns a random number draw from a distribution $N(\mu,\sigma)$. The output, ```shock```, can either be a scalar or a vector of dimension-$n$.

_Note_: ```randomnormal``` is an interface that calls ```randomnormal_scalar``` or ```randomnormal_vec``` depending on whether ```shock``` is a scalar or a vector.


[(back to index)](#inicio)


### cdfn
<a name="cdfn"></a>

```fortran
elemental function cdfn(x) result(f)
  implicit none
  real(kind=8), intent(in) :: x
  real(kind=8)             :: f
```

_Dependencies_: none

This function returns the cdf of a standard normal distribution, ```f``` $=\Phi(x)$. This subroutine is defined as ```elemental```, which implies that it can be call for both scalars and arrays.

_Example 1_: Scalar output

```fortran
print * , 'Result =', cdf(0.0)    ! Result = 0.500
print * , 'Result =', cdf(-1.0)   ! Result = 0.158
```

_Example 2_: Vector otput

```fortran
vec = (/ 0.0 , -1.0 /)

print * , 'Result =', cdf(vec)   ! Result = 0.500 0.158
```

[(back to index)](#inicio)

## Linear algebra

### vect
<a name="vect"></a>

```fortran
function vec(mat,byrow) result(vec)
  implicit none
  real(kind=8)       :: mat(:,:,...,:)
  real(kind=8)       :: vec(:)
  integer , optional :: byrow
```

_Dependencies_: none

This function returns a 1-dimensional array ```vec``` with all the elements of a user-supplied ```n```-dimensional array ```mat```, where ```n```$\leq5$. The optional argument ```byrow``` controls the order in which the elements of ```mat``` are stacked (see an example below).

_Example 1_:

```fortran
mat(:,1) = (/ 1 , 2 /)
mat(:,2) = (/ 3 , 4 /)

vec = vect(mat)

print * , 'vec =', vec   ! vec =  1.00  3.00  2.00  4.00

vec = vect(mat,byrow=1)

print * , 'vec =', vec   ! vec =  1.00  2.00  3.00  4.00

```


[(back to index)](#inicio)


### cumsum
<a name="cumsum"></a>

```fortran
function cumsum(vec0) result(vec1)
  implicit none
  real(kind=8) :: vec0(:),vec1(size(vec0))
```

_Dependencies_: none

Returns the cummulative sum of a vector ```vec0```.

_Example_:

```fortran
vec0 = (/ 1.0, 2.0, 1.0, 3.0 /) 

vec1 = cumsum(vec0) 

print * , 'vec1 = ', vec1   ! vec1 =  1.00  3.00  4.00  7.00
```

[(back to index)](#inicio)


### diag
<a name="diag"></a>

```fortran
function diag(mat) result(vec)
  implicit none
  real(kind=8) :: mat(:,:),vec(size(mat,dim=1))
```

_Dependencies_: none

This function returns the main diagonal of a matric ```mat```.

_Example_:

```fortran
mat(:,1) = (/ 1.0, 2.0, 1.0 /)
mat(:,2) = (/ 3.0, 3.0, 4.0 /)
mat(:,3) = (/ 5.0, 1.0, 3.0 /)

vec = diag(mat)

print * , 'vec = ', vec   ! vec =  1.00  3.00  3.00
```

[(back to index)](#inicio)


### transmat
<a name="transmat"></a>

```fortran
function transmat(mat) result(matt)
  implicit none
  real(kind=8) , intent(in)  :: mat(:,:)
  real(kind=8)               :: matt(size(mat,2),size(mat,1))
```

_Dependencies_: none

This function returns the transpose of a matrix ```mat```.



[(back to index)](#inicio)


### inverse
<a name="inverse"></a>

```fortran
function inverse(mat) result(imat)
  implicit none
  real(kind=8) :: mat(:,:)
  real(kind=8) :: imat(size(m,1),size(m,1))
```

_Dependencies_: none

This function returns the inverse of a squared matrix ```mat```.



[(back to index)](#inicio)


## Optimization

### simplex
<a name="simplex"></a>

```fortran
subroutine simplex(func,x,y,iy,ind,x0,itermax,tol,iprint)
  implicit none
  external                             :: func        ! user-supplied function to be minimize
  real(kind=8) , intent(out)           :: x(:)        ! values of "x" at minimum
  real(kind=8) , intent(out)           :: y           ! valuf of "func" at "x"
  integer      , intent(out)           :: iy          ! number of function evaluations
  integer      , intent(out)           :: ind         ! indicator of convergence
  real(kind=8) , intent(in)            :: x0(:)       ! initial guess
  real(kind=8) , intent(in) , optional :: tol         ! tolerance level
  integer      , intent(in) , optional :: itermax     ! max number of functione valuations
  integer      , intent(in) , optional :: iprint      ! indicator for printing behaviour
```

_Dependencies_: none

[(back to index)](#inicio)


### lmmin
<a name="lmmin"></a>

```fortran
subroutine lmmin_states_both(func,x,y,iy,ind,x0,itermax,damp,tol,toleach,shock,usebro,iprint)
  implicit none
  external                             :: func        ! user-supplied function to be minimize
  real(kind=8) , intent(out)           :: x(:)        ! values of "x" at minimum
  real(kind=8) , intent(out)           :: y(:)        ! valuf of "func" at "x"
  integer      , intent(out)           :: iy          ! number of function evaluations
  integer      , intent(out)           :: ind         ! indicator of convergence
  real(kind=8) , intent(in)            :: x0(:)       ! initial guess
  real(kind=8) , intent(in) , optional :: shock       ! shock to parameter values (as %)
  real(kind=8) , intent(in) , optional :: damp        ! damping factor
  real(kind=8) , intent(in) , optional :: tol         ! tolerance level
  real(kind=8) , intent(in) , optional :: toleach     ! tolerance level for each function
  integer      , intent(in) , optional :: itermax     ! max number of functione valuations
  integer      , intent(in) , optional :: iprint      ! indicator for printing behaviour
  integer      , intent(in) , optional :: usebro      ! indicator for the use of Broyden method to update Jacobian
```

_Dependencies_: ```broyden```, ```inverse```



[(back to index)](#inicio)


### golden
<a name="golden"></a>

```fortran
subroutine golden(func,x,y,xmax,xmin,itermax,tol)
  implicit none
  external                             :: func
  real(kind=8) , intent(out)           :: x
  real(kind=8) , intent(out)           :: y
  real(kind=8) , intent(in)            :: xmax
  real(kind=8) , intent(in)            :: xmin
  real(kind=8) , intent(in) , optional :: tol
  integer      , intent(in) , optional :: itermax
```

_Dependencies_: none




[(back to index)](#inicio)


### brent
<a name="brent"></a>

```fortran
subroutine brent(func,x,iy,ind,x0,x1,itermax,tol)
  implicit none
  external                             :: func
  real(kind=8) , intent(out)           :: x
  integer      , intent(out)           :: iy
  integer      , intent(out)           :: ind
  real(kind=8) , intent(in)            :: x0
  real(kind=8) , intent(in)            :: x1
  real(kind=8) , intent(in) , optional :: tol
  integer      , intent(in) , optional :: itermax
```

_Dependencies_: none




[(back to index)](#inicio)


### normalize
<a name="normalize"></a>

```fortran
subroutine normalize(y,x,xmax,xmin)
  implicit none
  real(kind=8) , intent(in)  :: xmax,xmin,x
  real(kind=8) , intent(out) :: y
```

_Dependencies_: none

This subroutine takes a bounded varibale  ```x```, contrained to be between ```xmin``` and ```xmax```, and applies the transformation 
$$ \texttt{y} = \log\left( \frac{ \texttt{x} - \texttt{xmin}}{\texttt{xmax} - \texttt{x}} \right) $$
to return an unbounded variable ```y```.

**Note**: This subroutine is useful to use stadard optimization algorithms for contrained problems. For isntance, one can use the Simplex method to minimize a function of ```x```, where ```x``` should be in the unit interval, by making use of ```normalize``` and ```denormalize```.

[(back to index)](#inicio)


### denormalize
<a name="denormalize"></a>

```fortran
subroutine denormalize(y,x,xmax,xmin)
  implicit none
  real(kind=8) , intent(in)  :: xmax,xmin,y
  real(kind=8) , intent(out) :: x
```

_Dependencies_: none

This subroutine takes an unbounded varibale ```y``` and applies the transformation 
$$ \texttt{x} = \texttt{xmin} + \left(\frac{\exp(y)}{1+\exp(y)} \right)(\texttt{xmax}-\texttt{xmin}) $$
to return a bounded variable ```x```, contrained to be between ```xmin``` and ```xmax```.

```fortran
call denormalize ( betau , beta , 1.0d0 , 0.0d0 )

! if betau = 3 --> beta = 0.952
! if betau = -2 --> beta = 0.119
```

**Note**: This subroutine is useful to use stadard optimization algorithms for contrained problems. For isntance, one can use the Simplex method to minimize a function of ```x```, where ```x``` should be in the unit interval, by making use of ```normalize``` and ```denormalize```.

[(back to index)](#inicio)


### broyden
<a name="broyden"></a>

```fortran
subroutine broyden(j1,j0,x1,x0,f1,f0)
  implicit none
  real(kind=8) , intent(in)  :: x1(:),f1(:)
  real(kind=8) , intent(in)  :: x0(:),f0(:)
  real(kind=8) , intent(in)  :: j0(:,:)
  real(kind=8) , intent(out) :: j1(:,:)
```

_Dependencies_: none

This subroutine applies the Boryden's method to update a Jacobian matrix.

Imagine we have an $m$-dimensional function $f$ in $n$ unknows. We evaluate two points $x_0$ and $x_1$, $f_1 = f(x_1)$ and $f_0 = f(x_0)$, and we compute the numerical jacobian of the function $f$ around $x=x_0$. This subroutine returns an opproximation to the jacobian matrix arounf the point $x=x_1$. The user must supply a pair of points ```x0``` and ```x1```, the value of the function evaluated at thos epoints ```f0``` and ```f1```, and the jacobian matrix ```j0``` evaluated at ```x0```.

You can learn more about this method in this [link](https://en.wikipedia.org/wiki/Broyden%27s_method).

[(back to index)](#inicio)


