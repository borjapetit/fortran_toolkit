# olsreg

```fortran
subroutine olsreg(coeffs,yvec,x1vec,x2vec,...,x8vec,w,mask,table)
  implicit none
  real(kind=8) , intent(out)           :: coeffs(:)
  real(kind=8) , intent(in)            :: yvec(:)
  real(kind=8) , intent(in)            :: x1vec(:)
  real(kind=8) , intent(in) , optional :: x2vec(:),x3vec(:),...,x8vec(:)
  real(kind=8) , intent(in) , optional :: w(:)     ! same length of "yvec"
  logical      , intent(in) , optional :: mask(:)  ! same length of "yvec"
  integer      , intent(in) , optional :: table
```

This subroutine returns the OLS coefficients from a linear regression model:
$$\texttt{yvec} = b_0 + b_1 \cdot\texttt{x1vec} + b_2\cdot \texttt{x2vec} + ... + b_8 \cdot\texttt{x8vec} + u$$
The subroutine allows up to 8 explanatory variables. The subroutine automatically includes a constant terms if ```size(coeffs)``` equals ```n```+1, where ```n``` be the number of explanatory variables.

The user can supply a vector of weigths, ```w``` of the same size as ```yvec```. If not supplied, the program assums uniform weigthing.

The user can also supply a ```mask``` to compute the impose a condition. The input ```mask``` is a logical array of the same size of ```yvec```. For example:

```fortran
call olsreg(yvec,x1vec,x2vex,mask = x1vecar.gt.0.0d0 .and. x2vex.lt.5.0d0)
```

This computes the OLS coefficients from a  linear regression of ```yvec``` on ```x1vec``` and ```x2vec``` conditional on ```x1vec``` being positive and ```x2vec``` being smaller than 5.

Finally, the variable ```table``` controls the output of the subroutine. If ```table``` is 0 (or missing), the subroutine returns the coefficients in the vector ```coeffs```. If ```table``` is 1, the subroutine prints a regression table with the coefficients and some additional statistics (t-stats, R-squared, etc).

**Dependencies**: [```varmean```](varmean.md),  [```varvar```](varvar.md)

[(back to index)](index.md)
