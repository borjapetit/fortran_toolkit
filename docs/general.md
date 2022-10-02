# General purpose

  - [```grid```](#grid): generate a grid for a continuous varibale
  - [```interpolation```](#interpolation): interpolate a value over a grid, returning position and distance
  - [```interpolate```](#interpolate): linearly interpolate a value over an n-dimensional grid, with n <= 6
  - [```timing```](#timing): returns the number of seconds since 00:00h of the 1st day of the month [robust to parelalization]
  - [```multiplo```](#multiplo): returns 1 if an integer is a multiple of another user-provided integer
  - [```iseven```](#iseven): returns 1 if a user-provided integer is even
  
## grid <a name="grid"></a>

```fortran
function grid(maxv,minv,n,s) result(v)
  implicit none
  real(kind=8)            :: maxv,minv,v(n)
  real(kind=8) , optional :: s
  integer                 :: n

  ! Dependencies: none
```

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


## interpolation <a name="interpolation"></a>
```fortran
subroutine interpolation(pos,wth,xnow,xgrid)
  implicit none
  real(kind=8), intent(in)  :: xnow,xgrid(:)
  real(kind=8), intent(out) :: wth
  integer     , intent(out) :: pos

  ! Dependencies: none
```
This subroutine finds the closest point in a given grid and return its position and the relative distance.

_Example_: consider a vector ```vec = [1,2,3]```. We can use ```interpolation``` to find the closest point of ```xnow = 2.3```:
```fortran
call interpolation(pos,wth,2.3,vec)

! Results: pos = 3, wth = 0.3
! 2.3 = vec(pos)*wth + vec(pos-1)*(one-wth) = 3*0.3 + 2*0.7
```
This subroutine is mainly used by the function ```interpolate```.

[(back to index)](#inicio)


## interpolate <a name="interpolate"></a>
returns the average of a variable, allowing for weigths
```fortran
function interpolate(x1,x2,...,xn,y1,y2,...,yn,mat) result(xi)
  implicit none
  real(kind=8) :: xi
  real(kind=8) :: x1,x2,...,xn
  real(kind=8) :: y1(:),y2(:),...,yn(:)
  real(kind=8) :: mat(:,:,...,:)

  ! Dependencies: interpolation
```
This function returns the linearly interpolated value of an n-dimensional function. The variables ```x1```, ```x2```, ..., ```xn``` are the values of the variables to be interpolated over their coresponding grids ```y1```, ```y2```, ...., ```yn```, and ```mat``` is an n-dimensional array with the results. The array ```mat``` must have at most, dimension 6.

_Example_: consider 2-dimensional function $f(x,y)$. We have an array ```mat``` whose $(i,j)$-element is the value of $f$ evaluated at the $i$-element of a grid for $x$, and the $j$-element of a grid for $y$. Then, we can interpolate the value of $f$ at $(x_0,y_0)$ using the ```interpolate```:
```fortran
xi = interpolate(x_0,y_0,x_grid,y_grid,mat)
```
If $x_0<\min($```x_grid```$)$ the subroutine take $\min($```x_grid```$)$ as the value of x. The same happens with $y$. Then, if $x_0<\min($```x_grid```$)$ and $y_0<\min($```y_grid```$)$, the function would return ```xi = mat(1,1)```.

[(back to index)](#inicio)


## timing <a name="timing"></a>
```fortran
function timing(mode) result(time)
  implicit none
  integer , optional :: mode
  real(kind=8)       :: time

  ! Dependencies: none
```
This functions returns a timing number that is robust to parallel computing. In particular, it returns the number of seconds since 00:00h of the 1st day of the month. The variable ```mode``` controls how time is measured:
 - If ```mode``` = 1, time is measured in seconds (default).
 - If ```mode``` = 2, time is measures in minutes.
 - If ```mode``` = 3, time is measured in hours


[(back to index)](#inicio)


## multiplo <a name="multiplo"></a>

```fortran
elemental function multiplo(num,xx) result(mul)
  implicit none
  integer :: num,xx,mul

  ! Dependencies: none
```
This function checks whether a number ```num0``` is a multiple of ```num1```. The result ```mul``` takes value 1 if ```num0``` is a multiple of ```num1```, and 0 otherwise.

_Example_: check whether 25 and 27 are multiples of 5:
```fortran
print * , multiplo(5,25)  ! result 1
print * , multiplo(5,27)  ! result 0
```

[(back to index)](#inicio)


## iseven <a name="iseven"></a>
```fortran
function iseven(num) result(ise)
  implicit none
  integer :: num,ise

  ! Dependencies: multiplo
```
This function checks whether a number ```num``` is even. The result ```ise``` takes value 1 if ```num``` is even, and 0 otherwise.


_Example_: check whether 3 and 6 are eve:
```fortran
check = iseven(3)  ! check = 0
check = iseven(6)  ! check = 1
```

[(back to index)](#inicio)

