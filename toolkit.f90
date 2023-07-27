
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! toolkit.f90, a toolkit for fortran90 programming
! Borja Petit, © 2022
!
! general purpose:
!   - grid: generate a grid for a continuous varibale
!   - interpolation: interpolate a value over a grid, returning position and distance
!   - interpolate: linearly interpolate a value over an n-dimensional grid, with n <= 6
!   - timing: returns the number of seconds since 00:00h of the 1st day of the month [robust to parelalization]
!   - multiplo: returns 1 if an integer is a multiple of another user-provided integer
!   - iseven: returns 1 if a user-provided integer is even
!
! statistics:
!   - varmean: returns the average of a variable, allowing for weigths
!   - varstd: returns the standard deviation of a variable, allowing for weigths
!   - correlation: returns the correlation of two variables, allowing for weigths
!   - percentile: returns the i-th percentile of a variables, allowing for weigths
!   - ols: returns the ols coefficients of a 1-var or 2-var regression, allowing for weigths
!   - fit2poli: fits a 2nd order polynomianl to a variable, allowing for weigths
!   - tachen: returns the transition matrix for a discretized ar(1) process
!   - randomnormal: returns a random draw for a nomal distribution
!   - cdfn: retutns the cdf of a nomabl distribution
!   - beta_noncentral_cdf: returns the cdf from a beta distribution
!
! linear algebra:
!   - vect: transform a matrix of dimension (n x m) into a vector of n·m rows
!   - cumsum: returns the vector with cummulative sum of a vector (as matlab's cumsum function)
!   - diag: returns the main diagonal of a matrix
!   - transmat: returns the transpose of a square matrix
!   - inverse: returns the invesrse of a sqaured matrix
!
! optimization (with and without states/fixed parameters)
!   - simplex: simplex algorithm
!   - lmmin: levenberg–marquardt algorithm
!   - golden: golden search algorithm
!   - brent: brent method
! 
! it also includes some other functions/subroutines used for optimization:
!   - broyden: updates a jacobian matrix using the broyden's method
!   - normalize: transform a bounded variable into an unbounded uno [for optimizaton]
!   - denormalize: transform a undonded variable into an counded uno [for optimizaton]
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module toolkit

  use, intrinsic :: IEEE_ARITHMETIC, only: IEEE_VALUE

  implicit none

  ! some general parameters
  integer  , parameter :: dp    = kind(1.0d00)
  real(dp) , parameter :: cero  = dble(0.000000000000000000)
  real(dp) , parameter :: medio = dble(0.500000000000000000)
  real(dp) , parameter :: uno   = dble(1.000000000000000000)
  real(dp) , parameter :: dos   = dble(2.000000000000000000)
  real(dp) , parameter :: cinco = dble(5.000000000000000000)
  real(dp) , parameter :: diez  = dble(10.00000000000000000)
  real(dp) , parameter :: cien  = dble(100.0000000000000000)
  real(dp) , parameter :: mil   = dble(1000.000000000000000)
  real(dp) , parameter :: tolvl = dble(0.000000000100000000)
  
  ! interface for the function "interpolate" in different dimessions.
  interface interpolate
    module procedure interpolate1d,interpolate2d,interpolate3d,interpolate4d,interpolate5d,interpolate6d
  end interface interpolate

  ! interface for the function "vect" (vectorize an array) in different dimessions.
  interface vect
    module procedure vectorize_in_2d,vectorize_in_3d,vectorize_in_4d,vectorize_in_5d,vectorize_in_6d,&
                     vectorize_dp_2d,vectorize_dp_3d,vectorize_dp_4d,vectorize_dp_5d,vectorize_dp_6d
  end interface vect

  ! interface for the function "randomnormal" to generate eithr an scalar or a vector of normal random numbers
  interface randomnormal
    module procedure randomnormal_scalar,randomnormal_vec
  end interface randomnormal

  contains









  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! general propuse functions and subroutines
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! ----------------------------------------------------------------------------
  ! this function creates a grid betwee "maxv" and "minv" with a curvature of "s"
  !   - if s>1: more grids points around "maxv"
  !   - if s<1: more grids points around "minv"

  function grid(maxv,minv,n,s) result(v)
    implicit none
    integer             :: i,n
    real(dp) , optional :: s
    real(dp)            :: maxv,minv,grid0(n),v(n),ss,xmin,xmax
    v(:) = cero ; ss = uno ; if (present(s)) ss = s
    xmin = min(maxv,minv)
    xmax = max(maxv,minv)
    if (ss.le.cero) call error(' errror in grid: spacing parameter is nonpositive')
    if (ss.gt.cero) then
      forall (i=1:n) grid0(i) = dble(i-1)/dble(n-1)
      forall (i=1:n) v(i) = (grid0(i)**ss)*(xmax-xmin) + minv
    end if
    return
  end function grid

  ! ----------------------------------------------------------------------------
  ! this subroutine finds the closest point in a given grid and return its
  ! position and the relative distance

  subroutine interpolation(pos,wth,xnow,xgrid)
    implicit none
    integer                :: j(1),n
    real(dp) , intent(in)  :: xnow,xgrid(:)
    real(dp) , intent(out) :: wth
    integer  , intent(out) :: pos
    real(dp)               :: xmin,xmax

    if (isnan(xnow)) then
      call error(' errror in interpolation: xnow is nan',1)
      return
    end if
    
    n = size(xgrid)
    if (n.eq.1) then
      call error(' errror in interpolation: xgrid is of dimension 1',1)
      return
    end if

    xmin = min(xgrid(1),xgrid(n))
    xmax = max(xgrid(1),xgrid(n))

    if (xnow.le.xmin) then
      if (xgrid(1).lt.xgrid(n)) then
        pos = 2 ; wth = cero
      else
        pos = n ; wth = uno
      end if
      return
    end if

    if (xnow.ge.xmax) then
      if (xgrid(n).gt.xgrid(1)) then
        pos = n ; wth = uno
      else
        pos = 2 ; wth = cero
      end if
      return
    end if 
    
    if (xnow.gt.xmin .and. xnow.lt.xmax) then
      j = minloc( xgrid - xnow , mask = xnow.le.xgrid )
      pos = j(1) ; wth = ( xnow - xgrid(pos-1) )/( xgrid(pos) - xgrid(pos-1) )
    end if
    if (pos.lt.1   ) call error(' errror in interpolation: pos < 1',1)
    if (pos.gt.n   ) call error(' errror in interpolation: pos > size',1)
    if (isnan(wth) ) call error(' errror in interpolation: wth is nan. is xgrid monotonic?',1)
    
    return
  end subroutine interpolation


  ! ----------------------------------------------------------------------------
  ! this set of functions return the linearly interpolated value of a
  ! n-dimensional vector with 6 >= n >= 1
  ! The user does not need oto call the specific subroutine, as the program will
  ! automatically call the appropiate uno depending on the inputs.
  ! To use these subroutine, the command is:
  !
  !     result = interpolate(val1,val2,...,valn,vector1,vector2,....vectorn,matrix,extra)
  !
  ! where "val1", "val2", ... are the values of the variables to be interpolated over
  ! their coresponding grids "vector1", "vector2, ...., and "matrix" is an
  ! n-dimensional array with the results.
  !
  ! The "extra" input is optional. If "extra" is set to 1, the subroutine will extrapolate 
  ! those  points "val" that are outside the grid "vector".
  ! This parameter is set to 0 by default
  
  function interpolate1d(x1,y1,m) result(xi)
    implicit none
    integer  :: pos
    real(dp) :: y1(:),m(:),x1,xi,wth
    if (size(y1).ne.size(m)) call error(' error in interpolate: 1st dimension incorrect')
    call interpolation(pos,wth,x1,y1)
    xi = m(pos)*wth + m(pos-1)*(uno-wth)
    return
  end function interpolate1d
  function interpolate2d(x1,x2,y1,y2,m) result(xi)
    implicit none
    integer  :: pos2
    real(dp) :: wth2
    real(dp) :: x1,x2,xi
    real(dp) :: y1(:),y2(:),m(:,:)
    if (size(m,2).ne.size(y2)) call error(' error in interpolate: 2nd dimension incorrect')
    call interpolation(pos2,wth2,x2,y2)
    xi = interpolate1d(x1,y1,m(:,pos2)  )*wth2 + &
         interpolate1d(x1,y1,m(:,pos2-1))*(uno-wth2)
    return
  end function interpolate2d
  function interpolate3d(x1,x2,x3,y1,y2,y3,m) result(xi)
    implicit none
    integer  :: pos3
    real(dp) :: y1(:),y2(:),y3(:),x1,x2,x3,xi
    real(dp) :: wth3,m(:,:,:)
    if (size(m,3).ne.size(y3)) call error(' error in interpolate: 3rd dimension incorrect')
    call interpolation(pos3,wth3,x3,y3)
    xi = interpolate2d(x1,x2,y1,y2,m(:,:,pos3)  )*wth3 + &
         interpolate2d(x1,x2,y1,y2,m(:,:,pos3-1))*(uno-wth3)
    return
  end function interpolate3d
  function interpolate4d(x1,x2,x3,x4,y1,y2,y3,y4,m) result(xi)
    implicit none
    integer  :: pos4
    real(dp) :: y1(:),y2(:),y3(:),y4(:),x1,x2,x3,x4,xi
    real(dp) :: wth4,m(:,:,:,:)
    if (size(m,4).ne.size(y4)) call error(' error in interpolate: 4th dimension incorrect')
    call interpolation(pos4,wth4,x4,y4)
    xi = interpolate3d(x1,x2,x3,y1,y2,y3,m(:,:,:,pos4)  )*wth4 + &
         interpolate3d(x1,x2,x3,y1,y2,y3,m(:,:,:,pos4-1))*(uno-wth4)
    return
  end function interpolate4d
  function interpolate5d(x1,x2,x3,x4,x5,y1,y2,y3,y4,y5,m) result(xi)
    implicit none
    integer  :: pos5
    real(dp) :: y1(:),y2(:),y3(:),y4(:),y5(:),x1,x2,x3,x4,x5,xi
    real(dp) :: wth5,m(:,:,:,:,:)
    if (size(m,5).ne.size(y5)) call error(' error in interpolate: 5th dimension incorrect')
    call interpolation(pos5,wth5,x5,y5)
    xi = interpolate4d(x1,x2,x3,x4,y1,y2,y3,y4,m(:,:,:,:,pos5))*wth5 + &
         interpolate4d(x1,x2,x3,x4,y1,y2,y3,y4,m(:,:,:,:,pos5-1))*(uno-wth5)
    return
  end function interpolate5d
  function interpolate6d(x1,x2,x3,x4,x5,x6,y1,y2,y3,y4,y5,y6,m) result(xi)
    implicit none
    integer  :: pos6
    real(dp) :: y1(:),y2(:),y3(:),y4(:),y5(:),y6(:),x1,x2,x3,x4,x5,x6,xi
    real(dp) :: wth6,m(:,:,:,:,:,:)
    if (size(m,5).ne.size(y5)) call error(' error in interpolate: 6th dimension incorrect')
    call interpolation(pos6,wth6,x6,y6)
    xi = interpolate5d(x1,x2,x3,x4,x5,y1,y2,y3,y4,y5,m(:,:,:,:,:,pos6))*wth6 + &
         interpolate5d(x1,x2,x3,x4,x5,y1,y2,y3,y4,y5,m(:,:,:,:,:,pos6-1))*(uno-wth6)
    return
  end function interpolate6d

  ! ----------------------------------------------------------------------------
  ! this functions returns a timing number that is robust to parallel computing
  ! it returns the number of seconds since 00:00h of the 1st day of the month.
  ! the variable "mode" controls how time is measured:
  !   - if mode = 1, time is measured in seconds (default).
  !   - if mode = 2, time is measures in minutes.
  !   - if mode = 3, time is measured in hours

  function timing(mode) result(time)
    implicit none
    integer , optional :: mode
    integer            :: v1(8),mod
    real(dp)           :: time0,time ; time0 = cero ; time = cero
    mod = 1 ; if (present(mode)) mod = mode
    call date_and_time(values=v1)  
    time0 = dble(v1(7)+60*(v1(6)+60*(v1(5)+24*(v1(3)+30*v1(2)))))
    if (mod.eq.3) then
      time = time0/dble(60*60)  ! measured in hours
    elseif (mod.eq.2) then
      time = time0/dble(60)     ! measured in minutes
    elseif (mod.eq.1) then
      time = time0              ! measured in seconds
    else
      call error("error in timing: invalid mode valua 8it should be between 1 and 3)")
    end if
    return
  end function timing

  ! ----------------------------------------------------------------------------
  ! this function checks whether a number "num" is a multiple of "xx".

  function multiplo(num,xx) result(mul)
    implicit none
    integer :: num,xx
    logical :: mul
    if (xx*floor(dble(num)/dble(xx)).eq.num) mul = .true.
    if (xx*floor(dble(num)/dble(xx)).ne.num) mul = .false.
    return
  end function multiplo

  ! ----------------------------------------------------------------------------
  ! this function checks whether a number "num" is an even number.

  function iseven(num) result(ise)
    implicit none
    integer :: num
    logical :: ise
    ise = multiplo(num,2)
    return
  end function iseven

  ! ----------------------------------------------------------------------------
  ! this subroutine prints an error message with the text "mess". If the optional 
  ! input "i" is present, the subroutine interrupts the execution until the user 
  ! types an integer number.

  subroutine error(mess,i)
    implicit none
    integer , optional , intent(in) :: i
    character(len=*)   , intent(in) :: mess
    character(len=10) :: typeany
    write(*,*) trim(adjustl(mess)) 
    if (present(i)) then
      write (*,'(a)',advance="no") ' Type any key to continue...'
      read (*,'(a)') typeany
    end if  
    return
  end subroutine error











  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! statistics
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function varmean(var,w,mask) result(meanvar)

    implicit none
    real(dp)               :: meanvar
    real(dp)               :: var(:)
    real(dp) , optional    :: w(:)
    logical  , optional    :: mask(:)
    logical  , allocatable :: mask1(:)
    real(dp) , allocatable :: weig(:)

    ! check mask
    allocate(mask1(size(var))) ; mask1 = .true.
    if (present(mask)) then
      if (size(var).eq.size(mask)) then
        mask1 = mask
      else
        call error('error in varmean: mask of incorrect size')
      end if
    end if

    ! check weigths
    allocate(weig(size(var))) ; weig(:) = uno
    if (present(w)) then
      if (size(var).ne.size(w)) call error('error in varmean!! var and w have different size')
      if (sum(w).lt.tolvl     ) call error('error in varmean!! w are cero')
      weig(:) = w(:)
    end if

    ! compute mean
    meanvar = sum(var(:)*weig(:),mask=mask1)/sum(weig,mask=mask1)

    deallocate(weig,mask1)

    return
  end function varmean

  ! ----------------------------------------------------------------------------
  ! this function returns the standard deviation of a variable "var"

  function varvar(var,w,mask) result(varian)

    implicit none
    real(dp)               :: varian
    real(dp)               :: var(:),mvar
    real(dp) , optional    :: w(:)
    logical  , optional    :: mask(:)
    logical  , allocatable :: mask1(:)
    real(dp) , allocatable :: weig(:)

    ! check mask
    allocate(mask1(size(var))) ; mask1 = .true.
    if (present(mask)) then
      if (size(var).eq.size(mask)) then
        mask1 = mask
      else
        call error('error in varvar: mask of incorrect size')
      end if
    end if

    ! check weigths
    allocate(weig(size(var))) ; weig(:) = uno
    if (present(w)) then
      if (size(var).ne.size(w)) call error('error in varvar!! var and w have different size')
      if (sum(w).lt.tolvl     ) call error('error in varvar!! w are cero')
      weig(:) = w(:)
    end if

    mvar   = varmean(var,weig,mask=mask1)
    varian = sum(weig*((var-mvar)**2.0d0),mask=mask1)/sum(weig,mask=mask1)

    deallocate(weig,mask1)

    return
  end function varvar

  ! ----------------------------------------------------------------------------
  ! this function returns the standard deviation of a variable "var"

  function varstd(var,w,mask) result(stdvar)

    implicit none
    real(dp)               :: stdvar
    real(dp)               :: var(:)
    real(dp) , optional    :: w(:)
    logical  , optional    :: mask(:)
    logical  , allocatable :: mask1(:)
    real(dp) , allocatable :: weig(:)

    ! check mask
    allocate(mask1(size(var))) ; mask1  = .true.
    if (present(mask)) then
      if (size(var).eq.size(mask)) then
        mask1 = mask
      else
        call error('error in varstd: mask of incorrect size')
      end if
    end if

    ! check weigths
    allocate(weig(size(var))) ; weig(:) = uno
    if (present(w)) then
      if (size(var).ne.size(w)) call error('error in varstd!! var and w have different size')
      if (sum(w).lt.tolvl     ) call error('error in varstd!! w are cero')
      weig(:) = w(:)
    end if

    ! compute standard deviation
    stdvar = sqrt(varvar(var,weig,mask=mask1))

    deallocate(weig,mask1)

    return
  end function varstd

  ! ----------------------------------------------------------------------------
  ! this function returns the correlation coefficient between two variables "xvar1" and "xvar2"

  function correlation(xvar1,xvar2,w,mask) result(corr)

    implicit none
    real(dp) , optional    :: w(:)
    real(dp)               :: xvar1(:),xvar2(:),corr
    real(dp)               :: aux1,aux2,aux3,aux4
    logical  , optional    :: mask(:)
    real(dp) , allocatable :: weig(:),weig1(:)
    real(dp) , allocatable :: vect1(:),vect2(:)

    ! initialize
    corr = cero
  
    ! check vector dimensions
    if (size(xvar1).ne.size(xvar2)) then
      call error('error in correaltion!! yvar and xvar of different sizes')
    end if

    ! check weigths
    allocate(weig(size(xvar1))) ; weig(:) = uno
    if (present(w)) then
      if (size(xvar1).ne.size(w)) then
        call error('error in correlation!! var and w have different size')
      end if
      if (sum(w).lt.tolvl) then
        call error('error in correlation!! w are cero')
      end if
      weig(:) = w(:)
    end if

    ! check mask
    if (present(mask)) then
      if (size(xvar1).eq.size(mask)) then
        allocate(vect1(count(mask))) ; vect1 = pack(xvar1,mask)
        allocate(vect2(count(mask))) ; vect2 = pack(xvar2,mask)
        allocate(weig1(count(mask))) ; vect1 = pack(weig ,mask)
      else
        call error('error in varmean: mask of incorrect size')
      end if
    else
      allocate(vect1(size(xvar1)),vect2(size(xvar1)),weig1(size(xvar1)))
      vect1 = xvar1
      vect2 = xvar2
      weig1 = weig
    end if

    ! compute means and stds
    aux1 = varmean(vect1,weig1)
    aux2 = varstd(vect1,weig1)
    aux3 = varmean(vect2,weig1)
    aux4 = varstd(vect2,weig1)

    ! check whether any of the variables is a constant
    if (aux2.lt.tolvl) then
      call error('error in correlation!! xvar1 is a constant')
    end if
    if (aux4.lt.tolvl) then
      call error('error in correlation!! xvar2 is a constant')
    end if

    ! compute correlation
    corr = (sum( weig1(:)*(vect1(:) - aux1)*(vect1(:) - aux3) )/sum(weig1))/(aux2*aux4)

    deallocate(weig,weig1,vect1,vect2)

    return
  end function correlation

  ! ----------------------------------------------------------------------------
  ! this function returns the percentile "pct" for a distribution "xvec"

  function percentile(xvec,pct,w,mask) result(cutoff)

    implicit none
    real(dp)               :: xvec(:),pct,cutoff
    real(dp)               :: aux1,aux2,aux3,aux4
    integer                :: iter
    real(dp) , optional    :: w(:)
    logical  , optional    :: mask(:)
    logical  , allocatable :: mask1(:)
    real(dp) , allocatable :: weig(:)
      
    ! check mask
    allocate(mask1(size(xvec))) ; mask1 = .true.
    if (present(mask)) then
      if (size(xvec).eq.size(mask)) then
        mask1 = mask
      else
        call error('error in percentile: mask of incorrect size')
      end if
    end if

    ! check weigths
    allocate(weig(size(xvec))) ; weig(:) = uno
    if (present(w)) then
      if (size(w).ne.size(w)) call error('error in percentile!! var and w have different size')
      if (sum(w).lt.tolvl   ) call error('error in percentile!! w are cero')
      weig(:) = w(:)
    end if

    ! check whether percentile is valid
    if (pct.gt.uno ) call error('error in percentile!! invalid percetile: larger than 100')
    if (pct.lt.cero) call error('error in percentile!! invalid percentile: negative value')

    ! find the percentile by bisection
    iter = 0 ; aux1 = maxval(xvec) ; aux2 = minval(xvec)
    do while ( abs(aux2-aux1).gt.tolvl .and. iter.lt.5000 ) ; iter = iter + 1
      aux3 = medio*(aux1+aux2)
      aux4 = sum(weig,mask = xvec.le.aux3 .and. mask1)/sum(weig,mask=mask1)
      if (aux4.le.pct) aux2 = aux3
      if (aux4.ge.pct) aux1 = aux3
    end do

    cutoff = medio*(aux1+aux2)

    deallocate(weig,mask1)

    return
  end function percentile

  ! ----------------------------------------------------------------------------
  ! these functions return the coefficients from a ols regresion.
  ! it allows for 1 or 2 explanatory variables, and returns the coefficients
  ! either as separate double precision variables or in a vector.
  ! if the vector coeffs has dimension equal to the number of variables, the subroutine 
  ! assumes that no conatnt is included in the regression. If the dimension of 
  ! coeffs is equal to 1 + the number of explanatory variables, the subroutine
  ! automathically includes a contant term.
  ! the user can specify the printing behaviour of the subroutine. if iprint is
  ! different from cero, the subroutine prints a regression table with additional
  ! statistics (t-stats, Rsquared, etc).

  subroutine olsreg(coeffs,yvec,x1vec,x2vec,x3vec,x4vec,x5vec,x6vec,x7vec,x8vec,w,mask,table)
    implicit none
    real(dp) , intent(out)           :: coeffs(:)
    real(dp) , intent(in)            :: yvec(:)
    real(dp) , intent(in)            :: x1vec(:)
    real(dp) , intent(in) , optional :: x2vec(:)
    real(dp) , intent(in) , optional :: x3vec(:)
    real(dp) , intent(in) , optional :: x4vec(:)
    real(dp) , intent(in) , optional :: x5vec(:)
    real(dp) , intent(in) , optional :: x6vec(:)
    real(dp) , intent(in) , optional :: x7vec(:)
    real(dp) , intent(in) , optional :: x8vec(:)
    real(dp) , intent(in) , optional :: w(:)
    logical  , intent(in) , optional :: mask(:)
    integer  , intent(in) , optional :: table
    logical  , allocatable           :: mask1(:)
    real(dp) , allocatable           :: yvar(:)
    real(dp) , allocatable           :: xvars(:,:)
    real(dp) , allocatable           :: zvar(:)
    real(dp) , allocatable           :: wvar(:)
    real(dp) , allocatable           :: xTx(:,:)
    real(dp) , allocatable           :: ixTx(:,:)
    real(dp) , allocatable           :: xTy(:)
    real(dp) , allocatable           :: wx(:,:)
    real(dp) , allocatable           :: evar(:)
    real(dp) , allocatable           :: sdbeta(:,:)
    real(dp) , allocatable           :: sdcoefs(:)
    real(dp) , allocatable           :: tstats(:)
    real(dp) , allocatable           :: inter(:,:)
    real(dp)                         :: rsq,arsq  
    integer                          :: i,n0,no,nc,wc,nx ; nx = 1

    ! number of observations
    n0 = size(yvec)
  
    ! mask and weights
    allocate(mask1(n0)) ; mask1 = .true.
    if (present(mask)) then
      if (size(mask).eq.n0) mask1 = mask
      if (size(mask).ne.n0) call error('error in olsreg!! mask of incorrect size',1)      
    end if
    allocate(wvar(n0)) ; wvar = uno
    if (present(w)) then
      if (size(w).eq.n0) wvar = w
      if (size(w).ne.n0) call error('error in olsreg!! weigths of incorrect size',1)
    end if

    ! check vector dimensions
    if (n0.ne.size(x1vec)) then ; nx = 1
      call error('error in olsreg!! yvec and x1vec different observations',1)
    end if
    if (present(x2vec)) then ; nx = 2
      if (n0.ne.size(x2vec)) call error('error in olsreg!! yvec and x2vec different observations',1)
    end if 
    if (present(x3vec)) then ; nx = 3
      if (n0.ne.size(x3vec)) call error('error in olsreg!! yvec and x3vec different observations',1)
    end if
    if (present(x4vec)) then ; nx = 4
      if (n0.ne.size(x4vec)) call error('error in olsreg!! yvec and x4vec different observations',1)
    end if
    if (present(x5vec)) then ; nx = 5
      if (n0.ne.size(x5vec)) call error('error in olsreg!! yvec and x5vec different observations',1)
    end if
    if (present(x6vec)) then ; nx = 6
      if (n0.ne.size(x6vec)) call error('error in olsreg!! yvec and x6vec different observations',1)
    end if
    if (present(x7vec)) then ; nx = 7
      if (n0.ne.size(x7vec)) call error('error in olsreg!! yvec and x7vec different observations',1)
    end if
    if (present(x8vec)) then ; nx = 8
      if (n0.ne.size(x8vec)) call error('error in olsreg!! yvec and x8vec different observations',1)
    end if 
    
    ! sizes
    no = count(mask1)  ! number of observations satisfying "mask"
    nc = size(coeffs)  ! number of coefficients
    wc = nc-nx         ! with constant = 1, 0 if not

    if (wc.ne.1 .and. wc.ne.0) call error('error in olsreg!! wc different from 1 and 0',1)

    ! allocate y-variable and weigths
    allocate(xvars(no,nc)) ; xvars = cero
    allocate(yvar(no))     ; yvar = pack(yvec,mask1) ! keep only valid observations
    allocate(zvar(no))     ; zvar = pack(wvar,mask1) ! keep only valid observations
    
    ! check weigths are positive    
    if (sum(zvar).lt.tolvl) call error('error in olsreg!! weigths are cero',1)
    
    ! fill constant (if exists)
    if (wc.eq.1) xvars(:,1) = zvar/zvar

    ! fill first explanatory variable
    xvars(:,1+wc) = pack(x1vec(:),mask1)

    ! fill remaining variables
    if (present(x2vec)) xvars(:,2+wc) = pack(x2vec(:),mask1)
    if (present(x3vec)) xvars(:,3+wc) = pack(x3vec(:),mask1)
    if (present(x4vec)) xvars(:,4+wc) = pack(x4vec(:),mask1)
    if (present(x5vec)) xvars(:,5+wc) = pack(x5vec(:),mask1)
    if (present(x6vec)) xvars(:,6+wc) = pack(x6vec(:),mask1)
    if (present(x7vec)) xvars(:,7+wc) = pack(x7vec(:),mask1)
    if (present(x8vec)) xvars(:,8+wc) = pack(x8vec(:),mask1)

    ! weigthed x-variables
    allocate(wx(no,nc))
    do i=1,no
      wx(i,:) = xvars(i,:)*zvar(i)
    end do

    ! allocate matrices for betas
    allocate(xTx(nc,nc),ixTx(nc,nc),xTy(nc))

    ! compute betas
    xTx    = matmul(transpose(xvars),wx)
    ixTx   = inverse(xTx)
    xTy    = matmul(transpose(xvars),yvar(:)*zvar(:))
    coeffs = matmul(ixTx,xTy)

    ! if iprint is different from cero, the subroutine prints a regression table
    if (present(table)) then
    if(table.eq.1) then

      allocate(evar(no),sdbeta(nc,nc),sdcoefs(nc),tstats(nc),inter(nc,2))

      evar    = yvar - matmul(xvars,coeffs)
      sdbeta  = ixTx*( sum(zvar*evar*evar) / dble(no-nc) )
      sdcoefs = sqrt(diag(sdbeta))
      tstats  = abs(coeffs/sdcoefs)

      inter(:,1) = coeffs - dble(1.9601068)*sdcoefs
      inter(:,2) = coeffs + dble(1.9601068)*sdcoefs

      if (wc.eq.1) then 
        rsq  = uno - varvar(evar,w=zvar)/varvar(yvar,w=zvar)
      else
        rsq  = uno - sum(zvar*(evar**2.0d0))/sum(zvar*(yvar**2.0d0))
      end if
      arsq = uno - ( (uno-rsq)*dble(no-1)/dble(no-nc) )

      write(*,99) '   '
      write(*,99) '   '
      write(*,'(a,i7  )') '                               Number of variables = ',nc
      write(*,'(a,i7  )') '                            Number of observations = ',n0
      write(*,'(a,i7  )') '                      Number of valid observations = ',no
      write(*,'(a,f7.4)') '                                         R-squared = ',rsq
      write(*,'(a,f7.4)') '                                Adjusted R-squared = ',arsq
      write(*,99) '   '
      write(*,99) ' -----------------------------------------------------------'
      write(*,99) '                beta     sd(b)      minb      maxb    t-stat'
      write(*,99) ' -----------------------------------------------------------'
      if (wc.eq.1) then
        write(*,99) ' Constant ' ,coeffs(wc),sdcoefs(wc),inter(wc,:),tstats(wc)
      end if
      do i=1,nx
        write(*,98) ' Var ' ,i, '    ',coeffs(i+wc),sdcoefs(i+wc),inter(i+wc,:),tstats(i+wc)
      end do
      write(*,99) ' -----------------------------------------------------------'
      write(*,99) '   '

      deallocate(evar,sdbeta,sdcoefs,tstats,inter)

    end if
    end if

    deallocate(xvars,xTx,ixTx,xTy,yvar,zvar,mask1,wvar)

    return
    99 format (a,f10.4,f10.4,f10.4,f10.4,f10.4)
    98 format (a,i1,a,f10.4,f10.4,f10.4,f10.4,f10.4)
  end subroutine olsreg

  ! ----------------------------------------------------------------------------
  ! this subroutine returns the transition matrix for a discretized ar(1) process
  ! with gaussian shocks

  subroutine tauchen(xvec,rho,mu,sigma,n,pmat)
    implicit none
    integer  , intent(in)  :: n
    real(dp) , intent(in)  :: rho,mu,sigma,xvec(n)
    real(dp) , intent(out) :: pmat(n,n)
    integer                :: i
    do i=1,n ; pmat(i,:) = cero
      call normaldist(xvec,mu+rho*xvec(i),sigma,n,pmat(i,:))
    end do
    return
  end subroutine tauchen

  ! ----------------------------------------------------------------------------
  ! this subroutine find the distribution of a variable "x" (discretized on a 
  ! "n"-dimensional grid "xvec") with  mean "mu" and stadard deviation "sigma"

  subroutine normaldist(xvec,mu,sigma,n,dist)
    implicit none
    integer  , intent(in)  :: n
    real(dp) , intent(in)  :: mu,sigma,xvec(n)
    real(dp) , intent(out) :: dist(n)
    real(dp)               :: xvec0(n),aux1
    integer                :: j
    xvec0 = (xvec-mu)/sigma
    dist(1) = cdfn( dble(0.5)*(xvec0(2)+xvec0(1)) )
    do j=2,n-1
      dist(j) = cdfn( dble(0.5)*(xvec0(j+1)+xvec0(j)) ) - cdfn( dble(0.5)*(xvec0(j)+xvec0(j-1)) )
    end do
    dist(n) = 1.0 - cdfn( dble(0.5)*(xvec0(n)+xvec0(n-1)) )
    aux1 = sum(dist)
    dist = dist(:)/aux1
    return
  end subroutine normaldist

  ! ----------------------------------------------------------------------------
  ! this function returns a random number from a normal distributino with mean
  ! "mu" and stadard deviation "std"

  subroutine randomnormal_scalar(shock,mu,std)
    implicit none
    real(dp) , intent(out) :: shock
    real(dp) , intent(in)  :: mu,std
    real(dp)               :: u,v,q
    real(dp) , parameter   :: s  = 0.449871 , t  = -0.386595 , a = 0.19600
    real(dp) , parameter   :: r1 = 0.275970 , r2 =  0.278460 , b = 0.25472
    do
      call random_number(u)
      call random_number(v)
      v = 1.7156 * (v - medio)
      q = (u-s)**2 + (abs(v)-t)*(a*(abs(v)-t) - b*(u-s))
      if (q < r1) exit
      if (q > r2) cycle
      if (v**2 < -dble(4.0)*log(u)*u**dble(2.00)) exit
    end do
    shock = std*v/u + mu
    return
  end subroutine randomnormal_scalar
  subroutine randomnormal_vec(shock,mu,std)
    implicit none
    integer                :: i
    real(dp) , intent(out) :: shock(:)
    real(dp) , intent(in)  :: mu,std
    real(dp)               :: q
    do i=1,size(shock)
      call randomnormal_scalar(shock(i),mu,std)
    end do
    q = sum(shock)/dble(size(shock))
    shock = shock - q + mu
    return
  end subroutine randomnormal_vec

  ! ----------------------------------------------------------------------------
  ! given a value "x" (an scalar or a vector), the function returns  f(x) = p(x<x),
  ! where f(·) is the cdf of a standard nomral distribution

  elemental function cdfn(x) result(f)
    implicit none
    real(dp) , intent(in)  :: x
    real(dp)               :: f,xabs,xsq
    real(dp) , parameter   :: a0  = 0.500000000000d0 , a1  = 0.398942280444d0
    real(dp) , parameter   :: a2  = 0.399903438504d0 , a3  = 5.758854804580d0
    real(dp) , parameter   :: a4  = 29.82135578080d0 , a5  = 2.624331216790d0
    real(dp) , parameter   :: a6  = 48.69599306920d0 , a7  = 5.928857244380d0
    real(dp) , parameter   :: b0  = 0.398942280385d0 , b1  = 3.8052d-8
    real(dp) , parameter   :: b2  = 1.000006153020d0 , b3  = 3.98064794d-4
    real(dp) , parameter   :: b4  = 1.986153813640d0 , b5  = 0.151679116635d0
    real(dp) , parameter   :: b6  = 5.293303249260d0 , b7  = 4.8385912808d0
    real(dp) , parameter   :: b8  = 15.15089724510d0 , b9  = 0.742380924027d0
    real(dp) , parameter   :: b10 = 30.78993303400d0 , b11 = 3.99019417011d0
    xabs = abs(x)
    xsq  = a0*(x**2)
    if (xabs <= 1.28d0)then
      f = a0-xabs*(a1-a2*xsq/(xsq+a3-a4/(xsq+a5+a6/(xsq+a7))))
    else if (xabs <= 12.7d0)then
      f = b0*exp(-xsq)/(xabs-b1+b2/(xabs+b3+b4/(xabs-b5+b6/(xabs+b7-b8/ &
          (xabs+b9+b10/(xabs+b11))))))
    else
      f = 0d0
    end if
    if (x > 0d0) f = 1d0-f
    return
  end function cdfn

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! linear algebra
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! ----------------------------------------------------------------------------
  ! these functions returns a vector with all the elements of a matrix.
  ! the function allows both double precission and integer matrices of
  ! 2-, 3- and 4-dimensions

  function vectorize_dp_2d(mat) result(vec)
    implicit none
    real(dp) :: mat(:,:)
    real(dp) :: vec(size(mat,1)*size(mat,2)) 
    vec = reshape(mat,(/size(vec)/))
    return
  end function vectorize_dp_2d
  function vectorize_dp_3d(mat) result(vec)
    implicit none
    real(dp) :: mat(:,:,:)
    real(dp) :: vec(size(mat,1)*size(mat,2)*size(mat,3)) 
    vec = reshape(mat,(/size(vec)/))
    return
  end function vectorize_dp_3d
  function vectorize_dp_4d(mat) result(vec)
    implicit none
    real(dp) :: mat(:,:,:,:)
    real(dp) :: vec(size(mat,1)*size(mat,2)*size(mat,3)*size(mat,4)) 
    vec = reshape(mat,(/size(vec)/))
    return
  end function vectorize_dp_4d
  function vectorize_dp_5d(mat) result(vec)
    implicit none
    real(dp) :: mat(:,:,:,:,:)
    real(dp) :: vec(size(mat,1)*size(mat,2)*size(mat,3)*size(mat,4)*size(mat,5))
    vec = reshape(mat,(/size(vec)/))
    return
  end function vectorize_dp_5d
  function vectorize_dp_6d(mat) result(vec)
    implicit none
    real(dp) :: mat(:,:,:,:,:,:)
    real(dp) :: vec(size(mat,1)*size(mat,2)*size(mat,3)*size(mat,4)*size(mat,5)*size(mat,6))
    vec = reshape(mat,(/size(vec)/))
    return
  end function vectorize_dp_6d
  function vectorize_in_2d(mat) result(vec)
    implicit none
    integer :: mat(:,:)
    integer :: vec(size(mat,1)*size(mat,2))
    vec = reshape(mat,(/size(vec)/))
    return
  end function vectorize_in_2d
  function vectorize_in_3d(mat) result(vec)
    implicit none
    integer :: mat(:,:,:)
    integer :: vec(size(mat,1)*size(mat,2)*size(mat,3)) 
    vec = reshape(mat,(/size(vec)/))
    return
  end function vectorize_in_3d
  function vectorize_in_4d(mat) result(vec)
    implicit none
    integer :: mat(:,:,:,:)
    integer :: vec(size(mat,1)*size(mat,2)*size(mat,3)*size(mat,4))
    vec = reshape(mat,(/size(vec)/))
    return
  end function vectorize_in_4d
  function vectorize_in_5d(mat) result(vec)
    implicit none
    integer :: mat(:,:,:,:,:)
    integer :: vec(size(mat,1)*size(mat,2)*size(mat,3)*size(mat,4)*size(mat,5))
    vec = reshape(mat,(/size(vec)/))
    return
  end function vectorize_in_5d
  function vectorize_in_6d(mat) result(vec)
    implicit none
    integer :: mat(:,:,:,:,:,:)
    integer :: vec(size(mat,1)*size(mat,2)*size(mat,3)*size(mat,4)*size(mat,5)*size(mat,6))
    vec = reshape(mat,(/size(vec)/))
    return
  end function vectorize_in_6d

  ! ----------------------------------------------------------------------------
  ! this function returns a vector "vec1" with the cummmulative sum of the elements
  ! of a vector "vec0"

  function cumsum(vec0) result(vec1)
    implicit none
    real(dp) :: vec0(:),vec1(size(vec0))
    integer  :: i
    vec1 = vec0(i)
    do i=2,size(vec0) 
      vec1(i) = sum(vec0(1:i))
    end do
    return
  end function cumsum

  ! ----------------------------------------------------------------------------
  ! this function returns an identity matriz of the same size as mat

  function eye(mat) result(mateye)
    implicit none
    real(dp) :: mat(:,:),mateye(size(mat,1),size(mat,2))
    integer  :: i
    if (size(mat,1).ne.size(mat,2)) call error(' error in diag: matrix not sqaure')
    mateye = cero ; forall (i=1:size(mat,1)) mateye(i,i) = uno
    return
  end function eye

  ! ----------------------------------------------------------------------------
  ! this function returns the main diagonal of a matric "mat"

  function diag(mat) result(vec)
    implicit none
    real(dp) :: mat(:,:),vec(size(mat,1))
    integer  :: i
    if (size(mat,1).ne.size(mat,2)) call error(' error in diag: matrix not sqaure')
    vec = cero ; forall (i=1:size(mat,1)) vec(i) = mat(i,i)
    return
  end function diag

  ! ----------------------------------------------------------------------------
  ! compute the inverse of a sqaure matrix "m"

  function inverse(m) result(im)
    implicit none
    integer  :: n,i,j,k
    real(dp) :: m(:,:),coeff
    real(dp) :: im(size(m,1),size(m,1)),b(size(m,1)),d(size(m,1)),x(size(m,1))
    real(dp) :: l(size(m,1),size(m,1)),u(size(m,1),size(m,1)),mb(size(m,1),size(m,1))
    if (size(m,1).ne.size(m,2)) call error(' error in inverse: matrix is not square')
    n = size(m,1) ; l = cero ; u = cero ; b = cero ; mb = m
    do k=1,n-1
      do i=k+1,n
        coeff  = m(i,k)/m(k,k)
        l(i,k) = coeff
        do j=k+1,n
          m(i,j) = m(i,j)-coeff*m(k,j)
        end do
      end do
    end do
    do i=1,n
      l(i,i) = 1.0
    end do
    do j=1,n
      do i=1,j
        u(i,j) = m(i,j)
      end do
    end do
    do k=1,n
      b(k) = 1.0
      d(1) = b(1)
      do i=2,n
        d(i) = b(i)
        do j=1,i-1
          d(i) = d(i) - l(i,j)*d(j)
        end do
      end do
      x(n)=d(n)/u(n,n)
      do i = n-1,1,-1
        x(i) = d(i)
        do j=n,i+1,-1
          x(i)=x(i)-u(i,j)*x(j)
        end do
        x(i) = x(i)/u(i,i)
        end do
       do i=1,n
         im(i,k) = x(i)
       end do
      b(k) = 0.0
    end do
    m = mb
    return
  end function inverse

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%











  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !
  ! optimization
  ! ^^^^^^^^^^^^
  !
  ! this code contains 4 optimization algorithms:
  !
  !  - golden search: single-value function with uno unknown
  !  - brent method: single-value function with uno unknown
  !  - simplex method: single-value function with "n" unknowns
  !  - levenberg–marquardt algorithm: system of "m" functions with "n" unknowns
  !
  ! the four algorithms can be used in problems with states (variables over which
  ! the algorithm does not optimize) which I found useful in many applications.
  ! for example, maximize a function f(x,y) over "x" for a given value of "y" (state)
  ! since I typically use this feature with state variables, the vector "states" is
  ! a vector of integers that indicates the position of each state variable on 
  ! its corresponding .
  !
  ! if your problem has states, please check the corresponding subroutine to make
  ! sure your function has the correct form.
  !
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



  ! ----------------------------------------------------------------------------
  ! golden search algorithm
  ! finds the maximum of a single-valued function with uno unknown
  ! the user should make sure the maximum of "func" is between "xmin" and "xmax"

  subroutine golden(func,x,y,numiter,xmax,xmin,itermax,tol,iprint)

    ! this subroutine finds the maximum of a single-valued function with uno unknown
    ! using the golden search algorithm
    !
    ! the inputs are:
    !
    !   - func: user-provided function to be minimize. this function should of the form:
    !
    !       function func(x) result(y)
    !         real(dp) :: x(:),y
    !       end function func
    !
    !   - xmin: lower bound of variable "x"
    !   - xmax: upper bound of variable "x"
    !
    ! there are also few optional inputs:
    !
    !   itermax:  max number of iterations (default = 500)
    !   tol:      tolerance level (default = 1.0d-8)
    !
    ! the outputs of the subroutine are:
    !
    !   x:        the value of x that maximizes "func"
    !   y:        the value of the function at the point x
    !   numiter:  number of function evaluations

    implicit none
    external               :: func
    real(dp) , intent(in)  :: xmax
    real(dp) , intent(in)  :: xmin
    real(dp) , intent(out) :: x
    real(dp) , intent(out) :: y
    integer  , intent(out) :: numiter
    
    real(dp) , intent(in) , optional :: tol
    integer  , intent(in) , optional :: itermax
    integer  , intent(in) , optional :: iprint
    real(dp) , parameter             :: alpha=0.61803399
    real(dp)                         :: x0,x1,x2,x3,f1,f2,toler
    integer                          :: maxiter,ipri

    interface
      function func(xx) result(ff)
      real(kind=8) :: xx,ff
      end function func
    end interface

    toler   = 1.0d-8 ; if (present(tol)    ) toler   = tol
    maxiter = 500    ; if (present(itermax)) maxiter = itermax
    ipri    = 0      ; if (present(iprint) ) ipri    = iprint

    if (ipri.ge.1) write(*,'(/,a,/)') ' starting golden search algorithm'

    x0 = min(xmax,xmin)
    x3 = max(xmax,xmin)

    if (abs(x0-x3).lt.toler) then
      if (ipri.ge.1) write(*,99) ' Not solved: xmax and xmin are too close'
      f1 = func(xmax) ; numiter = 1 ; x = xmax ; return
    end if

    x1 = alpha*x0 + (uno-alpha)*x3 ; f1 = func(x1) ; numiter = numiter + 1
    x2 = alpha*x3 + (uno-alpha)*x1 ; f2 = func(x2) ; numiter = numiter + 1

    do while (abs(x0-x3).gt.toler*(abs(x2)+abs(x1)) .and. abs(x0-x3).gt.toler .and. numiter.lt.maxiter)
      if (f2.gt.f1) then
        x0 = x1
        x1 = x2
        f1 = f2
        x2 = alpha*x1 + (uno-alpha)*x3 ; f2 = func(x2) ; numiter = numiter + 1
      else
        x3 = x2
        x2 = x1
        f2 = f1
        x1 = alpha*x2 + (uno-alpha)*x0 ; f1 = func(x1) ; numiter = numiter + 1 
      end if
      if (ipri.ge.1) then
        write(*,99) ' Iteration = ',numiter,'  | x1 = ',x1,'  x2 = ',x2,'  | f(x1) = ',f1,'  f(x2) = ',f2
      end if
    end do

    if (f1.gt.f2) then
      y = f1 ; x = x1
    else
      y = f2 ; x = x2
    end if

    if (ipri.ge.1) then
      write(*,99) ' '
      write(*,99) ' Solved: Iterations = ',numiter, '  | x = ',x,'  | f(x) = ',y
      write(*,99) ' '
    end if

    return
    99 format (a,i4,8(a,f10.4,a,f10.4))
  end subroutine golden

  ! ----------------------------------------------------------------------------
  ! brent's method
  ! root-finding algorithm combining the bisection and secant methods

  subroutine brent(func,x,numiter,exitcode,x0,x1,itermax,tol,iprint)

    ! this subroutine finds the root of a user-supplied single-valued function
    ! with one unknown using the brent's method
    !
    ! the inputs are:
    !
    !   - func: user-provided function to be minimize. this function should of the form:
    !
    !       function func(x) result(y)
    !         real(dp) :: x,y
    !       end function func
    !
    !   - x0: lower bound of variable "x" such that func(x0)<0
    !   - x1: upper bound of variable "x" such that func(x1)>0
    !
    ! there are also few optional inputs:
    !
    !   itermax:  max number of iterations (default = 500)
    !   tol:      tolerance level (default = 1.0d-8)
    !   iprint:   indicator for printing results (default = 0)
    !               0: print no results
    !               1: print main results
    !               2: print all the iterations
    !
    ! the outputs of the subroutine are:
    !
    !   x:         the value of x that maximizes "func"
    !   numiter:   number of function evaluations required
    !   exitcode:  exit indicator, taking value:
    !                0: the algorithm found a root
    !                1: either x0 or x1 is a root of func
    !                2: the root is not within the interval (x0,x1)
    !                3: the points x0 and x1 are too close
    !                4: the interval is too narrow
    !                9: maximum number of function evaluations reached

    implicit none
    external                         :: func
    real(dp) , intent(out)           :: x
    integer  , intent(out)           :: numiter
    integer  , intent(out)           :: exitcode
    real(dp) , intent(in)            :: x0,x1
    real(dp) , intent(in) , optional :: tol
    integer  , intent(in) , optional :: itermax
    integer  , intent(in) , optional :: iprint
    integer                          :: maxiter,ipri
    real(dp)                         :: toler,xa,ya,xb,yb,xc,yc,ys,xs

    interface
      function func(xx) result(resid)
        real(kind=8) :: xx,resid
      end function func
    end interface

    toler   = 1.0d-8 ; if (present(tol)    ) toler   = tol
    maxiter = 500    ; if (present(itermax)) maxiter = itermax
    ipri    = 0      ; if (present(iprint) ) ipri    = iprint

    if (ipri.ge.1) write(*,'(/,a,/)') ' starting brent algorithm'

    numiter = 0 ; ya = huge(uno) ; yb = ya ; yc = ya ; ys = ya

    ! if function is of equal sign at both initial points, return
    if (abs(x1-x0).lt.toler) then ; exitcode = 3
      if (ipri.ge.1) write(*,99) '   not solved: x1 and x0 too close'
      x = x1
      return
    end if

    xa = max(x0,x1) ; ya = func(xa) ; numiter = numiter + 1
    xb = min(x0,x1) ; yb = func(xb) ; numiter = numiter + 1

    ! if function is of equal sign at both initial points, existence is not guaranteed
    if (ya*yb.gt.cero .and. ipri.ge.1) then
      write(*,'(a,/)') '   root is not bracketed: existence is not guaranteed'
    end if
  
    ! if (ya*yb.gt.cero) then ; exitcode = 0 ; numiter = 2
    !   if (ipri.ge.1) write(*,99) '   not solved: invalid initial values'
    !   if (abs(ya).lt.abs(yb)) x = xa
    !   if (abs(ya).ge.abs(yb)) x = xb
    !   return
    ! end if

    ! if initial point is already a root, return
    if (abs(ya).lt.toler .or. abs(yb).lt.toler) then ; exitcode = 1 ; numiter = 0
      if (ipri.ge.1) write(*,99) '   not solved: initial point is already a root'
      if (abs(ya).lt.abs(yb)) x = xa
      if (abs(ya).ge.abs(yb)) x = xb
      return
    end if

    ! reorder initial points so that abs(ya) > abs(yb)
    if (abs(ya).lt.abs(yb)) then
      xc = xb ; xb = xa ; xa = min(x0,x1) ; yc = ya ; ya = yb
    end if
  
    xc = xa ; exitcode = 0

    do while (abs(xb-xa).gt.tolvl .and. numiter.lt.maxiter) 

      yc = func(xc) ; numiter = numiter + 1

      ! new guess
      if (abs(yc-ya).gt.toler .and. abs(yc-yb).gt.toler) then
        xs = xa*yb*yc/((ya-yb)*(ya-yc)) + xb*ya*yc/((yb-ya)*(yb-yc)) + xc*ya*yb/((yc-ya)*(yc-yb))
      else
        xs = xb - yb*(xb-xa)/(yb-ya)
      end if
      ys = func(xs) ; numiter = numiter + 1

      ! print iteration
      if (ipri.gt.1) write(*,98) numiter,xs,ys

      ! update points
      if (ya*ys.lt.cero) then
        xc = xb ; xb = xs ; yb = ys
      else
        xa = xs ; ya = ys
      end if

      ! re-order points
      if (abs(ya).lt.abs(yb)) then
        xc = xb ; xb = xa ; xa = xc ; yc = yb ; yb = ya ; ya = yc
      end if

      ! max iterations reached
      if (numiter.eq.maxiter) goto 9

      ! root found
      if (abs(ys).lt.toler) goto 10

      ! max iterations reached
      if (abs(xb-xa).le.tolvl) goto 11

    end do

    9 x = getbestpoint(xa,ya,xb,yb,xc,yc,xs,ys) ; exitcode = 9
    if (ipri.ge.1) write(*,97) 'interrupted: max iterations reached',numiter,x,func(x)
    return

    10 x = getbestpoint(xa,ya,xb,yb,xc,yc,xs,ys) ; exitcode = 0
    if (ipri.ge.1) write(*,97) 'solved: a root is found',numiter,x,func(x)
    return

    11 x = getbestpoint(xa,ya,xb,yb,xc,yc,xs,ys) ; exitcode = 4
    if (ipri.ge.1) write(*,97) 'interrupted: interval is too narrow',numiter,x,func(x)
    return

    99 format (a,f10.4,a,f10.4)
    98 format ('   Iteration = ',i4,' | xs = ',f10.4,' | f(xs) = ',f10.4)
    97 format ( /,'   ',a,/,/,'   iter = ',i4,/'   x    = ',f10.4,/,'   func = ',f10.4,/ )
    contains
    function getbestpoint(z1,f1,z2,f2,z3,f3,z4,f4) result (z5)
      implicit none
      real(dp) :: z1,f1,z2,f2,z3,f3,z4,f4,z5
      if (abs(f1).lt.min(abs(f2),abs(f3),abs(f4))) then
        z5 = z1 ; return
      elseif (abs(f2).lt.min(abs(f1),abs(f3),abs(f4))) then
        z5 = z2 ; return
      elseif (abs(f3).lt.min(abs(f1),abs(f2),abs(f4))) then
        z5 = z3 ; return
      elseif (abs(f4).lt.min(abs(f1),abs(f2),abs(f3))) then
        z5 = z4 ; return
      end if
    end function getbestpoint
  end subroutine brent

  ! ----------------------------------------------------------------------------
  ! simplex algorithm
  ! minimize a user-supplied single-valued function in "n" unknows

  subroutine simplex(func,x,y,numiter,exitcode,x0,itermax,tol,iprint)

    ! this subroutine minimizes a single-valued function with n unknowns using the
    ! simplex algorithm
    !
    ! the inputs are:
    !
    !   - func: user-provided function to be minimize. this function should of the form:
    !
    !       function func(x,states) result(y)
    !         real(dp) :: x(:),y
    !       end function func
    !
    !   - x0: the vector of dimension n with the initial guess
    !
    ! there are also few optional inputs:
    !
    !   itermax:  max number of iterations (default = 500)
    !   tol:      tolerance level (default = 1.0d-8)
    !   iprint:   indicator for printing results (default = 0)
    !               0: print no results
    !               1: print main results
    !               2: print all the iterations
    !
    ! the outputs of the subroutine are:
    !
    !   x:         the vector of x that solves the system
    !   y:         the value of the function at the point x
    !   numiter:   the number of total function evaluations
    !   exitcode:  exit indicator, taking value:
    !                0: system has converged, func(x) = 0
    !                1: simplex is too small
    !                9: max number of iterations reached

    implicit none

    external               :: func
    real(dp) , intent(out) :: x(:)
    real(dp) , intent(out) :: y
    integer  , intent(out) :: numiter
    integer  , intent(out) :: exitcode
    
    real(dp) , intent(in)            :: x0(:)
    real(dp) , intent(in) , optional :: tol
    integer  , intent(in) , optional :: itermax
    integer  , intent(in) , optional :: iprint

    integer              :: n,i,j,ilow,ihigh,ihigh2,maxiter,ipri
    real(dp)             :: xp(size(x,1),size(x,1)+1)
    real(dp)             :: yp(size(x,1)+1)
    real(dp)             :: xr(size(x,1))
    real(dp)             :: xm(size(x,1))
    real(dp)             :: xe(size(x,1))
    real(dp)             :: xc(size(x,1))
    real(dp)             :: y0,yr,ye,yc,cent,toler
    real(dp) , parameter :: alpha=1.0
    real(dp) , parameter :: beta=0.5
    real(dp) , parameter :: gamma=2.0

    interface
      function func(xvar) result(resid)
        real(kind=8) :: xvar(:),resid
      end function func
    end interface

    y0 = 0.0d-10 ; yp = 0.0d-10 ; ye = 0.0d-10 ; yc = 0.0d-10 ; yr = 0.0d-10
    xm = 0.0d-10 ; xp = 0.0d-10 ; xe = 0.0d-10 ; xc = 0.0d-10 ; xr = 0.0d-10

    toler   = 1.0d-8 ; if (present(tol)    ) toler   = tol      ! tolerance level; default: 1.0d-8
    maxiter = 500    ; if (present(itermax)) maxiter = itermax  ! max number of iterations; default: 500
    ipri    = 0      ; if (present(iprint) ) ipri    = iprint   ! printing behaviour; defualt: don't print

    if (ipri.ge.1) write(*,'(/,a,/)') ' starting simplex algorithm'

    ilow    = 1    
    n       = size(x,1)  ! number of variables
    y0      = func(x0)   ! value at initial point
    numiter = 1          ! interation counter
    yp(1)   = y0         ! initial value
    xp(:,1) = x0         ! initial simplex

    if (ipri.ge.1) write(*,93) y0

    if (y0.lt.toler) then
      yp(ilow) = y0 ; xp(:,ilow) = x0 ; goto 10
    end if

    if (ipri.ge.1) write(*,'(/,a,/)') ' evaluating initial simplex'
    
    xp = cero
    do j=1,n+1
      xp(:,j) = x0(:)
    end do

    do i=2,n+1
      if (x0(i-1).lt.cero) xp(i-1,i) = x0(i-1) - dble(0.20)
      if (x0(i-1).ge.cero) xp(i-1,i) = x0(i-1) + dble(0.20)
    end do

    do i=2,n+1
      yp(i) = func(xp(:,i))
      if (ipri.gt.1 .and. n.lt.5) write(*,97) i-1, yp(i) , xp(:,i)
      if (ipri.gt.1 .and. n.ge.5) write(*,97) i-1, yp(i)
      if (yp(i).lt.toler) then
        yp(ilow) = yp(i) ; xp(:,ilow) = xp(:,i) ; goto 10
      end if
    end do

    if (ipri.gt.1) write(*,'(/,a,/)') ' beginning iterations'

    1 continue

    ! sort oints
    ilow = 1
    do i=1,n+1
      if (yp(i).lt.yp(ilow)) ilow = i
    end do
    ihigh = 1
    do i=1,n+1
      if (yp(i).gt.yp(ihigh)) ihigh = i
    end do
    ihigh2 = ilow
    do i=1,n+1
      if (yp(i).gt.yp(ihigh2) .and. yp(i).lt.yp(ihigh)) ihigh2 = ihigh
    end do

    ! new iteration
    if (ipri.gt.1) write(*,97) numiter, yp(ilow)
    
    ! check if simplex is sufficiently large
    if (maxval(abs(xp(:,ihigh)-xp(:,ilow))).lt.toler) goto 11
    
    ! compute centroide
    do j=1,n
      cent = 0.0
      do i=1,n+1
        if (i.ne.ihigh) then
          cent = cent + xp(j,i)
        end if
      end do
      xm(j) = cent/dble(n)
    end do

    ! new point
    do j=1,n
      xr(j) = (uno+alpha)*xm(j) - alpha*xp(j,ihigh)
    end do
    yr = func(xr) ; numiter = numiter + 1

    ! check convergence
    if (yr.lt.toler) then
      yp(ilow) = yr ; xp(:,ilow) = xr(:) ; goto 10
    end if

    ! check number of function evaluations
    if (numiter.ge.maxiter) goto 19

    if (yr.lt.yp(ilow)) then

      ! new point
      do j=1,n
        xe(j) = (uno-gamma)*xm(j) - gamma*xp(j,ihigh)
      end do
      ye = func(xe) ; numiter = numiter + 1

      ! check convergence
      if (ye.lt.toler) then
        xp(:,ilow) = xe(:) ; yp(ilow) = ye ; goto 10
      end if

      if (ye.lt.yr) then
        xp(:,ihigh)  = xe(:) ; yp(ihigh)  = ye
      else
        xp(:,ihigh)  = xr(:) ; yp(ihigh)  = yr
      end if

      ! check number of function evaluations
      if (numiter.ge.maxiter) goto 19

    elseif (yr.gt.yp(ihigh2)) then

      ! if "yr" is worse, update worse point
      if (yr.lt.yp(ihigh)) then
        xp(:,ihigh) = xr(:) ; yp(ihigh) = yr
      end if

      ! new point
      do j=1,n
        xc(j) = (uno-beta)*xm(j) - beta*xp(j,ihigh)
      end do
      yc = func(xc) ; numiter = numiter + 1

      ! check convergence
      if (yc.lt.toler) then
        xp(:,ilow) = xc ; yp(ilow) = yc ; goto 10
      end if

      ! check number of function evaluations
      if (numiter.ge.maxiter) goto 19
      
      ! if "yc" is worse, update worse point
      if (yc.lt.yp(ihigh)) then

        xp(:,ihigh) = xc(:) ; yp(ihigh) = yc

      else

        do i=1,n+1

          if (i.ne.ilow) then

            ! new point
            do j=1,n
              xp(j,i) = medio*( xp(j,i) + xp(j,ilow) )
            end do
            yp(i) = func(xp(:,i)) ; numiter = numiter + 1

            ! check convergence
            if (yp(i).lt.toler) then
              xp(:,ilow) = xp(:,i) ; yp(ilow) = yp(i) ; goto 10
            end if

            ! check number of function evaluations
            if (numiter.ge.maxiter) goto 19

          end if
          
        end do

      end if

    else

      xp(:,ihigh) = xr(:)
      yp(ihigh)   = yr

    end if

    goto 1

    10 exitcode = 0 ; goto 99
    11 exitcode = 1 ; goto 99
    19 exitcode = 9 ; goto 99

    99 continue

    ilow = 1
    do i=1,n+1
      if (yp(i).lt.yp(ilow)) ilow = i
    end do
    
    if (ipri.ge.1) then
      if (exitcode.eq.0) write(*,'(/,a,/)') ' simplex finished: system has converged, func=0 '
      if (exitcode.eq.1) write(*,'(/,a,/)') ' simplex finished: system converged, simplex is very small'
      if (exitcode.eq.9) write(*,'(/,a,/)') ' simplex finished: max number of iterations'
      write (*,94) numiter
      write (*,93) yp(ilow)
      write (*,92) yp(ilow)/y0 - 1.0
      write (*, *) '  '
    end if

    x = xp(:,ilow)
    y = yp(ilow)

    return
    97 format ('  new point    ',i4,'  ',5(f10.4))
    94 format ('  iterations = ',i4)
    93 format ('  error      = ',f10.4)
    92 format ('  reduction  = ',f10.4)
  end subroutine simplex

  ! ----------------------------------------------------------------------------
  ! levenberg–marquardt algorithm
  ! minimize a user-supplied system of "m" equations in "n" unknows

  subroutine lmmin(func,x,y,numiter,exitcode,x0,itermax,damp,tol,toleach,shock,usebro,iprint)

    ! this subroutine minimizes the sum of squred errors of a system of m equations
    ! in n unknowns using the levenberg–marquardt algorithm. more about the algorithm
    ! can be found in https://en.wikipedia.org/wiki/Levenberg–Marquardt_algorithm.
    ! this method combines a newton-based algorithm and the gradient descent method.
    !
    ! the inputs are:
    !
    !   - func: user-provided function to be minimize. this function should of the form:
    !
    !       function func(x) result(y)
    !         real(dp) :: x(:)
    !         real(dp),allocatable :: y(:)
    !       end function func
    !
    !   - x0: the vector of dimension n with the initial guess
    !
    ! there are also few optional inputs:
    !
    !   itermax:  max number of iterations (default = 500)
    !   damp:     damping factor (default = 1)
    !   tol:      tolerance level (default = 1.0d-8)
    !   shock:    shock to compute numerical jacobian (default = 0.05)
    !   usebro:   whether to use broyden's method to update jacobian (default = 0)
    !             if broyden's method is used, the program will compute the
    !             numerical jacobian once, and then update it at every iteration
    !   iprint:   indicator for printing results (default = 0)
    !               0: print no results
    !               1: print main results
    !               2: print all the iterations
    !
    ! the outputs of the subroutine are:
    !
    !   x:         the vector of x that solves the system
    !   y:         the value of the function at the point x
    !   numiter:   the number of total function evaluations
    !   exitcode:  exit indicator, taking value:
    !                0: system has converged, func(x) = 0
    !                1: jacobian is closed to 0
    !                2: step in x is close to 0
    !                9: max number of iterations reached

    implicit none

    ! function to minimize
    external :: func
    
    ! outputs
    real(dp) , intent(out) :: x(:)
    real(dp) , intent(out) :: y(:)
    integer  , intent(out) :: numiter
    integer  , intent(out) :: exitcode

    ! inputs (including optionals)
    real(dp) , intent(in)            :: x0(:)
    real(dp) , intent(in) , optional :: shock
    real(dp)                         :: shck(size(x,1))
    real(dp) , intent(in) , optional :: damp
    real(dp)                         :: da
    real(dp) , intent(in) , optional :: tol
    real(dp)                         :: toler
    real(dp) , intent(in) , optional :: toleach
    integer  , intent(in) , optional :: itermax
    integer                          :: maxiter
    integer  , intent(in) , optional :: iprint
    integer                          :: ip
    integer  , intent(in) , optional :: usebro
    integer                          :: bro,br

    ! other variables
    integer  ::  i,k,n,m,qp
    real(dp) ::  j(size(y,1),size(x,1))
    real(dp) :: j0(size(y,1),size(x,1))
    real(dp) :: ij(size(x,1),size(x,1))
    real(dp) :: jj(size(x,1),size(x,1))
    real(dp) :: jt(size(x,1))
    real(dp) :: e0,y0(size(y,1))
    real(dp) :: eb,yb(size(y,1)),xb(size(x,1))
    real(dp) :: e1,y1(size(y,1)),x1(size(x,1))
    real(dp) :: ea,ya(size(y,1)),xa(size(x,1))
    real(dp) :: ej,yj(size(y,1),size(x,1)),xj(size(x,1),size(x,1))
    real(dp) :: dy(size(y,1)),dx(size(x,1)),ddx

    ! function
    interface
      function func(xvec) result(resid)
        real(kind=8)               :: xvec(:)
        real(kind=8) , allocatable :: resid(:)
      end function func
    end interface

    ! initialize optional inputs
    da        = 1.0d-0 ; if (present(damp)   ) da        = damp     ! damping factor
    toler     = 1.0d-8 ; if (present(tol)    ) toler     = tol      ! tolerance level
    shck(:)   = 5.0d-2 ; if (present(shock)  ) shck(:)   = shock    ! shock to parameters
    maxiter   = 500    ; if (present(itermax)) maxiter   = itermax  ! max number of iterations
    ip        = 0      ; if (present(iprint) ) ip        = iprint   ! indicator of printing behaviour
    bro       = 0      ; if (present(usebro) ) bro       = usebro   ! indicator for the use of Broyden Updating

    ! number of variables
    n = size(x,1)
    
    ! number of equations
    m = size(y,1)

    if (ip.ge.1) then
      write(*,'(/,a,/)') '  starting lmmin algorithm    '
      write(*,'(a,i2)' ) '  # equations   = ' , m
      write(*,'(a,i2)' ) '  # unknows     = ' , n
    end if

    ! initial point
    y0 = func(x0) ; e0 = sum(y0(:)*y0(:)) ; numiter = 1
    
    if (ip.ge.1) write(*,'(a,f10.4,/)') '  initial error = ' , e0

    ! best point so far = initial point
    yb = y0
    xb = x0
    eb = e0

    ! if converged, finish
    if (e0.lt.toler) goto 10

    ! if some function converged, return
    if (present(toleach)) then
      if (maxval(abs(y0)).lt.toleach) goto 13
    end if

    ! if reached max number of iterations, finish
    if (numiter.ge.maxiter) goto 19

    ! print results
    if ( ip.ge.1 ) write (*,91) eb
    if ( ip.eq.1 ) write (*, *) '  '

    qp = 0  ! this variable takes value 0 if no good point has been found with current jacobian, and 1 otherwise
    br = 0  ! this variable takes value 0 if new jacobian should be computed numerical (not use broyden)

    1 continue

    ! if converged, finish
    if (eb.lt.toler) goto 10

    ! if some function converged, return
    if (present(toleach)) then
      if (maxval(abs(yb)).lt.toleach) goto 13
    end if

    ! update numerical jacobian using broyden's method
    if ( bro.eq.1 .and. br.eq.1 .and. numiter.gt.1 ) then

      j0 = j
      dx(:) = x1(:)-xa(:)
      dy(:) = y1(:)-ya(:)
      do k = 1,m
        do i = 1,n
          j(k,i) = j0(k,i) + dx(i)*( dy(k) - sum(j0(k,:)*dx(:)) )/sum(dx(:)*dx(:))
        end do
      end do

    ! compute new numerical jacobian
    else ; if (ip.gt.1) write (*,*) '  '

      do i = 1,n

        ! shock the i-th parameter
        xj(:,i) = xb(:)
        if (x(i).le.cero) then
          xj(i,i) = xb(i) - shck(i)*abs(xb(i))
        elseif (x(i).ge.cero) then
          xj(i,i) = xb(i) + shck(i)*abs(xb(i))
        end if

        ! if shock is too small, increase it
        if (abs(xj(i,i)-xb(i)).lt.toler*dble(10.0)) then
          xj(i,i) = xb(i) + toler*dble(10.0)
        end if

        ! evaluate the function
        yj(:,i) = func(xj(:,i)) ; ej = sum(yj(:,i)*yj(:,i)) ; numiter = numiter + 1
        
        ! print result
        if (ip.gt.1) write (*,90) i,ej

        ! if converged, finish
        if (ej.lt.toler) then
          yb = yj(:,i) ; xb = xj(:,i) ; eb = ej ; goto 10
        end if

        ! if some function converged, return
        if (present(toleach)) then
          if (maxval(abs(yj(:,i))).lt.toleach) then
            yb = yj(:,i) ; xb = xj(:,i) ; eb = ej ; goto 13
          end if
        end if

        ! if reached max number of iterations, return best point so far
        if (numiter.ge.maxiter) then
          if (ej.lt.eb) then
            yb = yj(:,i) ; xb = xj(:,i) ; eb = ej
          end if
          goto 19
        end if

        ! if change in func is too large, compute new shocked parameter
        if (ej.gt.dble(50.0)*eb) then
          do while (ej.gt.dble(50.0)*eb .and. abs(xj(i,i)-xb(i)).gt.toler)
            xj(i,i) = (xb(i) + xj(i,i))/dble(2.0)
            yj(:,i) = func(xj(:,i)) ; ej = sum(yj(:,i)*yj(:,i)) ; numiter = numiter + 1
            if (ip.gt.1) write (*,90) i,ej
            if (numiter.ge.maxiter) then
              if (ej.lt.eb) then
                yb = yj(:,i) ; xb = xj(:,i) ; eb = ej
              end if
              goto 19
            end if
          end do
          shck(i) = abs( xb(i) - xj(i,i) ) / abs(xb(i))
        end if

        ! if change in func is too small, compute new shocked parameter
        if (abs(ej-eb).lt.tolvl) then ; ddx = xj(i,i) - xb(i)
          do while (abs(ej-eb).lt.tolvl)
            xj(i,i) = xj(i,i) + ddx
            yj(:,i) = func(xj(:,i)) ; ej = sum(yj(:,i)*yj(:,i)) ; numiter = numiter + 1
            if (ip.gt.1) write (*,90) i,ej
            if (numiter.ge.maxiter) then
              if (ej.lt.eb) then
                yb = yj(:,i) ; xb = xj(:,i) ; eb = ej
              end if
              goto 19
            end if
          end do
          shck(i) = abs( xb(i) - xj(i,i) ) / abs(xb(i))
        end if

        ! compute numerical jacobian
        do k = 1,m
          j(k,i) = ( yj(k,i) - yb(k) )/( xj(i,i) - xb(i) )
        end do

      end do
      qp = 0

      if (ip.gt.1) write (*,*) '  '
    end if

    ! for each iteration, starting point = best point
    ya = yb ; ea = eb ; xa = xb ; br = 1

    ! begin iteration
    2 continue

    ! compute new point
    jj = matmul(transpose(j),j) 
    jt = matmul(transpose(j),yb)
    do i=1,n
      jj(i,i) = (uno+da)*jj(i,i)
    end do
    ij = inverse(jj)
    x1 = xb - matmul(ij,jt)

    ! if jacobian is too close to cero, return
    if (maxval(abs(jt)).lt.toler) goto 11

    ! if the step in "x" is too small, return
    if (sum(abs(x1-xb)).le.toler) goto 12

    ! evaluate new point
    y1 = func(x1) ; e1 = sum(y1(:)*y1(:)) ; numiter = numiter + 1
  
    ! print iterations
    if (ip.ge.1) print 97 , numiter , e1 , eb , da

    ! check number of iterations
    if (numiter.ge.maxiter) then
      if (e1.lt.eb) then
        yb = y1 ; xb = x1 ; eb = e1   
      end if 
      goto 19
    end if

    ! if the point improves, reduce dumping factor and new jacobian
    if (e1.lt.eb) then
      yb = y1             ! if good point, update best point
      xb = x1             ! if good point, update best point
      eb = e1             ! if good point, update best point
      da = da/dble(2.0)   ! if good point, decrease the damping factor
      qp = 1              ! if good point, set parameter qp to 1
      goto 1              ! if good point, update jacobian (with broyden if first point, numerical if not)
    end if

    ! if the point is bad, keep increasing the duminping factor
    if (e1.ge.eb) then
      da = da*dble(5.00)      ! if bad point, increase the damping factor
      br = 0                  ! if bad point, don't update jacobian with broyden
      if (qp.eq.0) goto 2     ! if bad point, and the first after jacobian, get another point
      if (qp.eq.1) then       ! if bad point, but not the first after jacobian, compute new jacobian
        qp = 0
        goto 1
      end if
    end if

    10 exitcode = 0 ; goto 7
    11 exitcode = 1 ; goto 7
    12 exitcode = 2 ; goto 7
    13 exitcode = 3 ; goto 7
    19 exitcode = 9 ; goto 7

    7 if (ip.ge.1) then
      if (exitcode.eq.0) write(*,'(/,a,/)') ' lmmin finished: system has converged, func=0 '
      if (exitcode.eq.1) write(*,'(/,a,/)') ' lmmin finished: system has converged, jac close to 0'
      if (exitcode.eq.2) write(*,'(/,a,/)') ' lmmin finished: system has converged, step in x close to 0'
      if (exitcode.eq.3) write(*,'(/,a,/)') ' lmmin finished: some equations is 0'
      if (exitcode.eq.9) write(*,'(/,a,/)') ' lmmin finished: max number of iterations'
      write(*,*          ) ' iterations = ', numiter
      write(*,'(a,f10.4)') ' error      = ', eb
      write(*,'(a,f10.4)') ' reduction  = ', eb/e0 - 1.0
      write(*,*          ) '  '
    end if

    y = yb
    x = xb

    return
    97 format ('  new point    ',i4,'  ',5(f10.4))
    91 format ('  first point  ','    ','  ',5(f10.4))
    90 format ('   - jacobian, param =',i4,f10.4)
  end subroutine lmmin

  ! ----------------------------------------------------------------------------
  ! this subroutine applies the transofrmation 
  !   y = log( (x-xmin)/(xmax-xmin) / (1-(x-xmin)/(xmax-xmin)) )
  ! to return a unbounded variable "y" from a bounded variable "x"

  subroutine normalize(y,x,xmax,xmin)
    implicit none
    real(dp) , intent(in)  :: xmax,xmin,x
    real(dp) , intent(out) :: y
    if (xmax.lt.xmin) call error(' errror in normalize: xmax > xmin')
    if (x.gt.xmax   ) call error(' errror in normalize: x > xmax')
    if (x.lt.xmin   ) call error(' errror in normalize: x < xmin')
    y = log((x-xmin)/(xmax-xmin)/max(0.0001,uno-(x-xmin)/(xmax-xmin)))
    return
  end subroutine normalize

  ! ----------------------------------------------------------------------------
  ! this subroutine applies the transofrmation 
  !   x = xmin + ( exp(y)/(1+exp(y) )*(xmax-xmin) 
  ! to return a bounded variable "x" from an unbounded variable "y"

  subroutine denormalize(y,x,xmax,xmin)
    implicit none
    real(dp) , intent(in)  :: xmax,xmin,y
    real(dp) , intent(out) :: x
    if (xmax.lt.xmin) call error(' errror in denormalize: xmax > xmin')
    x = (exp(y)/(uno+exp(y)))*(xmax-xmin) + xmin
    if (isnan(x) .and. y.gt.cien) x = xmax
    return
  end subroutine denormalize

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module toolkit
