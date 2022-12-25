# simplex

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

**Note**: for constrained optimization problems, one can make use of the [```normalize and denormalize```](normalize.md) subroutines.

**Dependencies**: none

[(back to index)](../index.md)

---

**Example**
