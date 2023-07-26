
### olsreg

```fortran
subroutine olsreg(coeffs,yvar,xvar1,xvar2,...,xvar8,w,mask,table)
  implicit none
  real(kind=8) , intent(out)           :: coeffs(:)
  real(kind=8) , intent(in)            :: yvar(:)
  real(kind=8) , intent(in)            :: xvar1(:)
  real(kind=8) , intent(in) , optional :: xvar2(:),xvar3(:),...,xvar8(:)
  real(kind=8) , intent(in) , optional :: w(:)     
  logical      , intent(in) , optional :: mask(:)
  integer      , intent(in) , optional :: table
```

This subroutine returns the OLS coefficients from a linear regression model:

$$\texttt{yvar} = \beta_0 + \beta_1 \cdot \texttt{xvar1} + \beta_2 \cdot\texttt{xvar2} + ... + \beta_8 \cdot\texttt{xvar8} + u$$

The subroutine allows up to 8 explanatory variables. The subroutine automatically includes a constant terms if ```size(coeffs)``` equals ```n```+1, where ```n``` is the number of explanatory variables. The user can supply a vector of weights, ```w``` of the same size as ```yvar```. If not supplied, the program assums uniform weigthing.

The user can also supply a ```mask``` to compute the impose a condition. The input ```mask``` is a logical array of the same size of ```yvar```. For example:

```fortran
call olsreg(coeffs,yvar,xvar1,xvar2,mask = xvar1.gt.0.0d0 .and. xvar2.lt.5.0d0)
```

This code computes the OLS coefficients from a linear regression of ```yvar``` on ```xvar1``` and ```xvar2``` conditional on ```xvar1``` being positive and ```xvar2``` being smaller than 5. If ```coeffs``` is of size 3, the model includes a contant term so that $\texttt{coeffs} = (\beta_0,\beta_1,\beta_2)$. If the size of ```coefs``` is 2, then the returned vector is $\texttt{coeffs} = (\beta_1,\beta_2)$.

Finally, the variable ```table``` controls the output of the subroutine. If ```table``` is 0 (or missing), the subroutine returns the coefficients in the vector ```coeffs```. If ```table``` is 1, the subroutine prints a regression table with the coefficients and some additional statistics (t-stats, R-squared, etc).

**Dependencies**: [```error```](error.md), [```varmean```](varmean.md),  [```varvar```](varvar.md)

[(back to index)](../index.md)

---

**Example**

Imagine we have four vectors (```x0```, ```x1```, ```x2```, ```x3```) each with 100 normal random numbers (all with mean zero). And we define a vector ```y``` such that

$$\texttt{y} = 0.10 + \texttt{x0} + 0.7\cdot\texttt{x1} - 0.5\cdot \texttt{x2} + 0.2\cdot\texttt{x3}$$

We want to estimate the following regression model 
$$\texttt{y} = \beta_0 + \beta_1 \cdot \texttt{x1} + \beta_2 \cdot\texttt{x2} +  \beta_3 \cdot\texttt{x3} + u$$

To do so we first define a vector ```coeffs``` of dimension 4 and the run:

```fortran
call olsreg(coeffs,y,x1,x2,x3)

! coeffs = (/ 0.1000 , 0.7265 , -0.4519 , 0.2535 /)
```
If we specify ```table = 1```, the command prints the following:
```fortran
call olsreg(coeffs,y,x1,x2,x3,table=1)

! Output:
!
!                               Number of variables =       4
!                             Number of observatios =     100
!                       Number of valid observatios =     100
!                                         R-squared =  0.6240
!                                Adjusted R-squared =  0.6122
!   
! -----------------------------------------------------------
!                beta     sd(b)      minb      maxb    t-stat
! -----------------------------------------------------------
! Constant     0.1000    0.0311    0.0390    0.1610    3.2141
! Var 1        0.8205    0.0839    0.6560    0.9850    9.7767
! Var 2       -0.4489    0.0664   -0.5791   -0.3187    6.7593
! Var 3        0.1603    0.0526    0.0571    0.2635    3.0456
! -----------------------------------------------------------
```

If we want to run the same regression model without a constant, we just define the vector ```coeffs``` to have dimension 3. 

```fortran
call olsreg(coeffs,y,x1,x2,x3)

! coeffs = (/ 0.7265 , -0.4519 , 0.2535 /)

call olsreg(coeffs,y,x1,x2,x3,table=1)

! Output:
!   
!                               Number of variables =       3
!                             Number of observatios =     100
!                       Number of valid observatios =     100
!                                         R-squared =  0.5997
!                                Adjusted R-squared =  0.5915
!   
! -----------------------------------------------------------
!                beta     sd(b)      minb      maxb    t-stat
! -----------------------------------------------------------
! Var 1        0.8205    0.0879    0.6483    0.9927    9.3379
! Var 2       -0.4489    0.0695   -0.5852   -0.3126    6.4560
! Var 3        0.1603    0.0551    0.0523    0.2683    2.9089
! -----------------------------------------------------------
```
_Note: the R-squared in a regression without a constant term is not well defined. This subroutine deals with this issue as Stata does. See the following [thread](https://www.statalist.org/forums/forum/general-stata-discussion/general/1324612-problem-with-sst-and-ssr-formula-in-a-regression-without-constant) for a discussion on this._
