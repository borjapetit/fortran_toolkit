
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! toolkit.f90, a toolkit for fortran90 programming
! Borja Petit, © 2021
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
!   - normalize: transform a bounded variable into an unbounded one [for optimizaton]
!   - denormalize: transform a undonded variable into an counded one [for optimizaton]
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module toolkit

  implicit none

  integer  , parameter :: dp = kind(1.0d00)

  real(kind=8) , parameter :: zero  = dble(0.00000000000000)
  real(kind=8) , parameter :: one   = dble(1.00000000000000)
  real(kind=8) , parameter :: cien  = dble(100.000000000000)
  real(kind=8) , parameter :: mil   = dble(1000.00000000000)
  real(kind=8) , parameter :: tolvl = dble(0.00000000010000)
  real(kind=8) , parameter :: half  = dble(0.50000000000000)

  interface interpolate
    module procedure interpolate1d,interpolate2d,interpolate3d,interpolate4d,interpolate5d,interpolate6d
  end interface interpolate

  interface vect
    module procedure vectorize_int_2d,vectorize_int_3d,vectorize_int_4d,vectorize_int_5d,&
                     vectorize_dp_2d,vectorize_dp_3d,vectorize_dp_4d,vectorize_dp_5d
  end interface vect

  ! interface reg
  !   module procedure ols1var,ols2var,ols3var,ols4var
  ! end interface reg

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
    integer                 :: i,n
    real(kind=8) , optional :: s
    real(kind=8)            :: maxv,minv,grid0(n),v(n),ss,xmin,xmax
    v(:) = zero ; ss = one ; if (present(s)) ss = s
    if (ss.le.zero) call error(' errror in grid: spacing parameter is nonpositive')
    do i=1,n ; grid0(i)=dble(i-1)/dble(n-1) ; end do
    xmin = max(maxv,minv)
    xmax = min(maxv,minv)
    if (ss.gt.zero) then
      do i = 1,n ; v(i) = (grid0(i)**ss)*(xmax-xmin) + minv ; end do
    end if
    return
  end function grid

  ! ----------------------------------------------------------------------------
  ! this subroutine finds the closest point in a given grid and return its
  ! position and the relative distance

  subroutine interpolation(pos,wth,xnow,xgrid)
    implicit none
    integer                   :: j,n
    real(kind=8), intent(in)  :: xnow,xgrid(:)
    real(kind=8), intent(out) :: wth
    integer     , intent(out) :: pos
    if (isnan(xnow)) call error(' errror in interpolation: xnow is nan')
    n = size(xgrid)
    if (n.eq.1) then
      pos = 1
      wth = one
    else
      if (xnow.le.xgrid(1)) then
        pos = 2
        wth = zero
      else if (xnow.ge.xgrid(n)) then
        pos = n
        wth = one
      else if (xnow.gt.xgrid(1) .and. xnow.lt.xgrid(n)) then
        j = 2
        do while (xnow.gt.xgrid(j))
          j = j + 1
        end do
        pos = j
        wth = (xnow-xgrid(j-1))/(xgrid(j)-xgrid(j-1))
      end if
    end if
    if (pos.lt.1   ) call error(' errror in interpolation: pos < 1')
    if (pos.gt.n   ) call error(' errror in interpolation: pos > size')
    if (isnan(wth) ) call error(' errror in interpolation: wvar is nan')
    return
  end subroutine interpolation

  ! ----------------------------------------------------------------------------
  ! this set of functions return the linearly interpolated value of a
  ! n-dimensional vector with 6 >= n >= 1
  ! The user does not need oto call the specific subroutine, as the program will
  ! automatically call the appropiate one depending on the inputs.
  ! To use these subroutine, the command is:
  !
  !     result = interpolate(val1,val2,...,valn,vector1,vector2,....vectorn,matrix)
  !
  ! where "val1", "val2", ... are the values of the variables to be interpolated over
  ! their coresponding grids "vector1", "vector2, ...., and "matrix" is an
  ! n-dimensional array with the results.

  function interpolate1d(x1,y1,m) result(xi)
    implicit none
    integer      :: pos
    real(kind=8) :: y1(:),m(:),x1,xi,wth
    if (size(y1).ne.size(m)) call error(' error in interpolate: 1st dimension incorrect')
    call interpolation(pos,wth,x1,y1)
    xi = m(pos)*wth + m(pos-1)*(one-wth)
    return
  end function interpolate1d
  function interpolate2d(x1,x2,y1,y2,m) result(xi)
    implicit none
    integer      :: pos2
    real(kind=8) :: wth2
    real(kind=8) :: x1,x2,xi
    real(kind=8) :: y1(:),y2(:),m(:,:)
    if (size(m,2).ne.size(y2)) call error(' error in interpolate: 2nd dimension incorrect')
    call interpolation(pos2,wth2,x2,y2)
    xi = interpolate1d(x1,y1,m(:,pos2)  )*wth2 + &
         interpolate1d(x1,y1,m(:,pos2-1))*(one-wth2)
    return
  end function interpolate2d
  function interpolate3d(x1,x2,x3,y1,y2,y3,m) result(xi)
    implicit none
    integer      :: pos3
    real(kind=8) :: y1(:),y2(:),y3(:),x1,x2,x3,xi
    real(kind=8) :: wth3,m(:,:,:)
    if (size(m,3).ne.size(y3)) call error(' error in interpolate: 3rd dimension incorrect')
    call interpolation(pos3,wth3,x3,y3)
    xi = interpolate2d(x1,x2,y1,y2,m(:,:,pos3)  )*wth3 + &
         interpolate2d(x1,x2,y1,y2,m(:,:,pos3-1))*(one-wth3)
    return
  end function interpolate3d
  function interpolate4d(x1,x2,x3,x4,y1,y2,y3,y4,m) result(xi)
    implicit none
    integer      :: pos4
    real(kind=8) :: y1(:),y2(:),y3(:),y4(:),x1,x2,x3,x4,xi
    real(kind=8) :: wth4,m(:,:,:,:)
    if (size(m,4).ne.size(y4)) call error(' error in interpolate: 4th dimension incorrect')
    call interpolation(pos4,wth4,x4,y4)
    xi = interpolate3d(x1,x2,x3,y1,y2,y3,m(:,:,:,pos4)  )*wth4 + &
         interpolate3d(x1,x2,x3,y1,y2,y3,m(:,:,:,pos4-1))*(one-wth4)
    return
  end function interpolate4d
  function interpolate5d(x1,x2,x3,x4,x5,y1,y2,y3,y4,y5,m) result(xi)
    implicit none
    integer      :: pos5
    real(kind=8) :: y1(:),y2(:),y3(:),y4(:),y5(:),x1,x2,x3,x4,x5,xi
    real(kind=8) :: wth5,m(:,:,:,:,:)
    if (size(m,5).ne.size(y5)) call error(' error in interpolate: 5th dimension incorrect')
    call interpolation(pos5,wth5,x5,y5)
    xi = interpolate4d(x1,x2,x3,x4,y1,y2,y3,y4,m(:,:,:,:,pos5))*wth5 + &
         interpolate4d(x1,x2,x3,x4,y1,y2,y3,y4,m(:,:,:,:,pos5-1))*(one-wth5)
    return
  end function interpolate5d
  function interpolate6d(x1,x2,x3,x4,x5,x6,y1,y2,y3,y4,y5,y6,m) result(xi)
    implicit none
    integer      :: pos6
    real(kind=8) :: y1(:),y2(:),y3(:),y4(:),y5(:),y6(:),x1,x2,x3,x4,x5,x6,xi
    real(kind=8) :: wth6,m(:,:,:,:,:,:)
    if (size(m,5).ne.size(y5)) call error(' error in interpolate: 6th dimension incorrect')
    call interpolation(pos6,wth6,x6,y6)
    xi = interpolate5d(x1,x2,x3,x4,x5,y1,y2,y3,y4,y5,m(:,:,:,:,:,pos6))*wth6 + &
         interpolate5d(x1,x2,x3,x4,x5,y1,y2,y3,y4,y5,m(:,:,:,:,:,pos6-1))*(one-wth6)
    return
  end function interpolate6d

  ! ----------------------------------------------------------------------------
  ! this functions returns a timing number that is robust to parallel computing
  ! it returns the number of seconds since 00:00h of the 1st day of the month.
  ! the variable "mode" controls how time is measured:
  !   - if mode=1, time is measured in seconds (default).
  !   - if mode=2, time is measures in minutes.
  !   - if mode=3, time is measured in hours

  function timing(mode) result(time)

    implicit none
    integer , optional :: mode
    integer            :: v1(8)
    real(kind=8)       :: time,v2(8) ; time = zero ; v2 = zero
    
    call date_and_time(values=v1)  

    v2(2) = dble(v1(2)*30*24*60*60) ! months
    v2(3) = dble(v1(3)*   24*60*60) ! days
    v2(5) = dble(v1(5)*      60*60) ! hours
    v2(6) = dble(v1(6)*         60) ! minutes
    v2(7) = dble(v1(7)            ) ! seconds

    ! measured in hours
    if (present(mode) .and. mode.eq.3) then 
      time = sum(v2)/dble(60*60)
    ! measured in minutes
    else if (present(mode) .and. mode.eq.2) then 
      time = sum(v2)/dble(60)
    ! measured in seconds
    else 
      time = sum(v2)
    end if

    return
  end function timing

  ! ----------------------------------------------------------------------------
  ! this function checks whether a number "num" is a multiple of "xx".

  function multiplo(num,xx) result(mul)
    implicit none
    integer :: num,xx,mul
    if (xx*floor(dble(num)/dble(xx)).eq.num) mul = 1
    if (xx*floor(dble(num)/dble(xx)).ne.num) mul = 0
    return
  end function multiplo

  ! ----------------------------------------------------------------------------
  ! this function checks whether a number "num" is an even number.

  function iseven(num) result(ise)
    implicit none
    integer :: num,ise
    ise = multiplo(num,2)
    return
  end function iseven

  ! ----------------------------------------------------------------------------
  ! this subroutine prints an error message.

  subroutine error(mess)
    implicit none
    integer                       :: i
    character(len=*) , intent(in) :: mess
    write(*,*) trim(adjustl(mess)) ; read * , i
    return
  end subroutine













  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! statistics
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! ----------------------------------------------------------------------------
  ! this function returns the mean of a variable "var".

  function varmean(var,wvar) result(meanvar)

    implicit none
    real(kind=8) , optional :: wvar(:)
    real(kind=8)            :: var(:),meanvar,weig(size(var))
    
    weig(:) = one
    meanvar = zero

    if (present(wvar)) then
      if (size(var).ne.size(wvar)) call error('error in varmean!! var and wvar have different size')
      if (sum(wvar).lt.tolvl     ) call error('error in varmean!! wvar are zero')
      weig(:) = wvar(:)
    end if

    meanvar = sum(var(:)*weig(:))/sum(weig)

    return
  end function varmean

  ! ----------------------------------------------------------------------------
  ! this function returns the standard deviation of a variable "var"

  function varstd(var,wvar) result(stdvar)

    implicit none
    real(kind=8)            :: var(:)
    real(kind=8) , optional :: wvar(:)
    real(kind=8)            :: stdvar
    real(kind=8)            :: weig(size(var)),mvar

    weig(:) = one

    if (present(wvar)) then
      if (size(var).ne.size(wvar)) call error('error in varstd!! var and wvar have different size')
      if (sum(wvar).lt.tolvl     ) call error('error in varstd!! wvar are zero')
      weig(:) = wvar(:)
    end if

    mvar   = varmean(var,weig)
    stdvar = sqrt(sum(weig(:)*((var(:)-mvar)**dble(2.00)))/sum(weig))
    return

  end function varstd

  ! ----------------------------------------------------------------------------
  ! this function returns the correlation coefficient between two variables "xvar1" and "xvar2"

  function correlation(xvar1,xvar2,wvar) result(corr)

    implicit none
    real(kind=8) , optional :: wvar(:)
    real(kind=8)            :: xvar1(:),xvar2(:),corr,weig(size(xvar1))
    real(kind=8)            :: aux1,aux2,aux3,aux4,aux5

    weig(:) = one
    corr    = zero

    if (size(xvar1).ne.size(xvar2)) call error('error in correaltion!! yvar and xvar of different sizes')
    if (present(wvar)) then
      if (size(xvar1).ne.size(wvar)) call error('error in correlation!! var and wvar have different size')
      if (sum(wvar).lt.tolvl       ) call error('error in correlation!! wvar are zero')
      weig(:) = wvar(:)
    end if

    aux1 = varmean(xvar1,weig)
    aux2 = varstd(xvar1,weig)
    aux3 = varmean(xvar2,weig)
    aux4 = varstd(xvar2,weig)

    if (aux2.lt.tolvl) call error('error in correlation!! xvar1 is a constant')
    if (aux4.lt.tolvl) call error('error in correlation!! xvar2 is a constant')

    aux5 = sum(weig(:)*( (xvar1(:) - aux1)*(xvar2(:) - aux3) ))/sum(weig)
    corr = aux5/(aux2*aux4)

    return
  end function correlation

  ! ----------------------------------------------------------------------------
  ! this function returns the percentile "pct" for a distribution "xvec"

  function percentile(xvec,pct,wvar) result(cutoff)

    implicit none
    real(kind=8)            :: xvec(:),pct
    real(kind=8) , optional :: wvar(:)
    real(kind=8)            :: cutoff
    real(kind=8)            :: weig(size(xvec))
    real(kind=8)            :: aux1,aux2,aux3,aux4
    integer                 :: iter

    weig(:) = one
    cutoff  = zero

    if (present(wvar)) then
      if (size(xvec).ne.size(wvar)) call error('error in percentile!! var and wvar have different size')
      if (sum(wvar).lt.tolvl      ) call error('error in percentile!! wvar are zero')
      weig(:) = wvar(:)
    end if

    if (pct.gt.one ) call error('error in percentile!! invalid percetile: larger than 100')
    if (pct.lt.zero) call error('error in percentile!! invalid percentile: negative value')

    iter = 0
    aux1 = maxval(xvec)
    aux2 = minval(xvec)
    do while ( abs(aux2-aux1).gt.tolvl .and. iter.lt.5000 ) ; iter = iter + 1
      aux3 = half*(aux1+aux2)
      aux4 = sum(weig,mask = xvec.le.aux3)/sum(weig)
      if (aux4.le.pct) aux2 = aux3
      if (aux4.ge.pct) aux1 = aux3
    end do

    cutoff = half*(aux1+aux2)

    return
  end function percentile

  ! ----------------------------------------------------------------------------
  ! these functions return the coefficients from a ols regresion.
  ! it allows 1 and 2 explanatory variables, and returns the coefficients
  ! either as separate double precision variables or in a vector.
  ! The user can safely call the subroutine "ols(·)" without the need to specify 
  ! the subroutine. The program will automatically call the corresponding subroutine
  ! depending on the number of variables and on the format of output variables

  subroutine olsreg(coefs,yvec,x1vec,x2vec,x3vec,x4vec,x5vec,x6vec,x7vec,x8vec,wvec,iprint)
    implicit none
    real(kind=8) , intent(out)           :: coefs(:)
    real(kind=8) , intent(in)            :: yvec(:)
    real(kind=8) , intent(in)            :: x1vec(:)
    real(kind=8) , intent(in) , optional :: x2vec(:)
    real(kind=8) , intent(in) , optional :: x3vec(:)
    real(kind=8) , intent(in) , optional :: x4vec(:)
    real(kind=8) , intent(in) , optional :: x5vec(:)
    real(kind=8) , intent(in) , optional :: x6vec(:)
    real(kind=8) , intent(in) , optional :: x7vec(:)
    real(kind=8) , intent(in) , optional :: x8vec(:)
    real(kind=8) , intent(in) , optional :: wvec(:)
    integer      , intent(in) , optional :: iprint
    real(kind=8) , allocatable           :: xvars(:,:),xTx(:,:),ixTx(:,:),xTy(:)
    real(kind=8)                         :: wvar(size(yvec))
    real(kind=8)                         :: evar(size(yvec)),shat
    real(kind=8)                         :: sdbeta(size(coefs),size(coefs))
    real(kind=8)                         :: sdcoefs(size(coefs))
    real(kind=8)                         :: tstats(size(coefs))
    real(kind=8)                         :: inter(size(coefs),2)
    real(kind=8)                         :: pvals(size(coefs))
    integer                              :: j,i,nx ; nx = 1

    coefs = zero

    if (present(wvec)) then
      if (size(wvec).eq.size(yvec)) then
        wvar = wvec
      else
        call error('error in olsreg!! weigths of incorrect size')
      end if
    else
      wvar = one
    end if

    if (sum(wvar).lt.tolvl) then
      call error('error in olsreg!! weigths are zero')
    end if
    if (size(yvec).ne.size(x1vec)) then
      call error('error in olsreg!! yvec and x1vec different observations')
    end if
    if (present(x2vec)) then ; nx = 2 ; if (size(yvec).ne.size(x2vec)) then
      call error('error in olsreg!! yvec and x2vec different observations')
    end if ; end if
    if (present(x3vec)) then ; nx = 3 ; if (size(yvec).ne.size(x3vec)) then
      call error('error in olsreg!! yvec and x3vec different observations')
    end if ; end if
    if (present(x4vec)) then ; nx = 4 ; if (size(yvec).ne.size(x4vec)) then
      call error('error in olsreg!! yvec and x4vec different observations')
    end if ; end if
    if (present(x5vec)) then ; nx = 5 ; if (size(yvec).ne.size(x5vec)) then
      call error('error in olsreg!! yvec and x5vec different observations')
    end if ; end if
    if (present(x6vec)) then ; nx = 6 ; if (size(yvec).ne.size(x6vec)) then
      call error('error in olsreg!! yvec and x6vec different observations')
    end if ; end if
    if (present(x7vec)) then ; nx = 7 ; if (size(yvec).ne.size(x7vec)) then
      call error('error in olsreg!! yvec and x7vec different observations')
    end if ; end if
    if (present(x8vec)) then ; nx = 8 ; if (size(yvec).ne.size(x8vec)) then
      call error('error in olsreg!! yvec and x8vec different observations')
    end if ; end if
    
    if (size(coefs).eq.nx+1) then ; j = 1
      allocate(xvars(size(yvec),nx+1)) ; xvars(:,1) = wvar(:)
    elseif (size(coefs).eq.nx) then ; j = 0
      allocate(xvars(size(yvec),nx))
    else
      call error('error in olsreg!! coefs of incorrect size')
    end if

    xvars(:,1+j) = x1vec(:)*wvar(:)

    if (present(x2vec)) xvars(:,2+j) = x2vec(:)*wvar(:)
    if (present(x3vec)) xvars(:,3+j) = x3vec(:)*wvar(:)
    if (present(x4vec)) xvars(:,4+j) = x4vec(:)*wvar(:)
    if (present(x5vec)) xvars(:,5+j) = x5vec(:)*wvar(:)
    if (present(x6vec)) xvars(:,6+j) = x6vec(:)*wvar(:)
    if (present(x7vec)) xvars(:,7+j) = x7vec(:)*wvar(:)
    if (present(x8vec)) xvars(:,8+j) = x8vec(:)*wvar(:)

    allocate(xTx(size(xvars,2),size(xvars,2)))
    allocate(ixTx(size(xvars,2),size(xvars,2)))
    allocate(xTy(size(xvars,2)))

    xTx    = matmul(transmat(xvars),xvars)
    ixTx   = inverse(xTx)
    xTy    = matmul(transmat(xvars),yvec(:)*wvar(:))
    coefs  = matmul(ixTx,xTy)

    if (present(iprint) .and. iprint.eq.1) then

     evar    = yvec - matmul(xvars,coefs)
     sdbeta  = ixTx*( sum(evar(:)*evar(:)) / (dble(size(yvec)) - size(xvars,2)) )
     sdcoefs = sqrt(diag(sdbeta))
     tstats  = abs(coefs(:)/sdcoefs(:))
     inter(:,1) = coefs - dble(1.96)*sdcoefs(:)
     inter(:,2) = coefs + dble(1.96)*sdcoefs(:)

     write(*,99) '   '
     write(*,99) ' -----------------------------------------------------------'
     write(*,99) '                beta     sd(b)      minb      maxb    t-stat'
     write(*,99) ' -----------------------------------------------------------'
     if (j.eq.1) then
      write(*,99) ' Constant ' ,coefs(1),sdcoefs(1),inter(1,:),tstats(1)
     end if
     do i=1,nx
      write(*,98) ' Var ' ,i, '    ',coefs(i+j),sdcoefs(i+j),inter(i+j,:),tstats(i+j)
     end do
     write(*,99) ' -----------------------------------------------------------'
     write(*,99) '   '

    end if

    deallocate(xvars,xTx,ixTx,xTy)

    return
    99 format (a,f10.4,f10.4,f10.4,f10.4,f10.4)
    98 format (a,i1,a,f10.4,f10.4,f10.4,f10.4,f10.4)
  end subroutine olsreg

  ! ----------------------------------------------------------------------------
  ! this function returns the transition matrix for a discretized ar(1) process
  ! with gaussian shocks

  subroutine tauchen(xvec,rho,mu,sigma,n,pmat)
    implicit none
    integer     , intent(in)  :: n
    real(kind=8), intent(in)  :: rho,mu,sigma,xvec(n)
    real(kind=8), intent(out) :: pmat(n,n)
    integer                   :: i
    do i=1,n
      pmat(i,:) = zero
      call normaldist(xvec,mu+rho*xvec(i),sigma,n,pmat(i,:))
    end do
    return
  end subroutine tauchen

  ! ----------------------------------------------------------------------------
  ! this function returns the transition matrix for a discretized ar(1) process
  ! with gaussian shocks

  subroutine normaldist(xvec,mu,sigma,n,dist)
    implicit none
    integer     , intent(in)  :: n
    real(kind=8), intent(in)  :: mu,sigma,xvec(n)
    real(kind=8), intent(out) :: dist(n)
    real(kind=8)              :: xvec0(n),aux1
    integer                   :: j
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
    real(kind=8) , intent(out) :: shock
    real(kind=8) , intent(in)  :: mu,std
    real(kind=8)               :: u,v,q
    real(kind=8) , parameter   :: s  = 0.449871 , t  = -0.386595 , a = 0.19600
    real(kind=8) , parameter   :: r1 = 0.275970 , r2 =  0.278460 , b = 0.25472
    do
      call random_number(u)
      call random_number(v)
      v = 1.7156 * (v - half)
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
    integer                    :: i
    real(kind=8) , intent(out) :: shock(:)
    real(kind=8) , intent(in)  :: mu,std
    real(kind=8)               :: u,v,q
    real(kind=8) , parameter   :: s  = 0.449871 , t  = -0.386595 , a = 0.19600
    real(kind=8) , parameter   :: r1 = 0.275970 , r2 =  0.278460 , b = 0.25472
    do i=1,size(shock)
      do
        call random_number(u)
        call random_number(v)
        v = 1.7156 * (v - half)
        q = (u-s)**2 + (abs(v)-t)*(a*(abs(v)-t) - b*(u-s))
        if (q < r1) exit
        if (q > r2) cycle
        if (v**2 < -dble(4.0)*log(u)*u**dble(2.00)) exit
      end do
      shock(i) = std*v/u + mu
    end do
    return
  end subroutine randomnormal_vec

  ! ----------------------------------------------------------------------------
  ! given a value "x", the function returns the f(x) = p(x<x), where f(·) is the
  ! cdf of a standard nomral distribution

  elemental function cdfn(x) result(f)
    implicit none
    real(kind=8) , intent(in)  :: x
    real(kind=8)               :: f,xabs,xsq
    real(kind=8) , parameter   :: a0  = 0.500000000000d0 , a1  = 0.398942280444d0
    real(kind=8) , parameter   :: a2  = 0.399903438504d0 , a3  = 5.758854804580d0
    real(kind=8) , parameter   :: a4  = 29.82135578080d0 , a5  = 2.624331216790d0
    real(kind=8) , parameter   :: a6  = 48.69599306920d0 , a7  = 5.928857244380d0
    real(kind=8) , parameter   :: b0  = 0.398942280385d0 , b1  = 3.8052d-8
    real(kind=8) , parameter   :: b2  = 1.000006153020d0 , b3  = 3.98064794d-4
    real(kind=8) , parameter   :: b4  = 1.986153813640d0 , b5  = 0.151679116635d0
    real(kind=8) , parameter   :: b6  = 5.293303249260d0 , b7  = 4.8385912808d0
    real(kind=8) , parameter   :: b8  = 15.15089724510d0 , b9  = 0.742380924027d0
    real(kind=8) , parameter   :: b10 = 30.78993303400d0 , b11 = 3.99019417011d0
    xabs = abs(x)
    xsq  = a0*x**2
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
  ! these functions returns a vector with all the elements of a matrix. the
  ! parameter "byrow" controls how the matrix is vectorize:
  !
  !    byrow = 0 (default)           byrow = 1
  !
  !   | a b |       | a |       | a b |       | a |
  !   | c d |  -->  | c |       | c d |  -->  | b |
  !                 | b |                     | c |
  !                 | d |                     | d |
  !
  ! the function allows both double precission and integer matrices of
  ! 2-, 3- and 4-dimensions

  function vectorize_dp_2d(mat,byrow) result(vec)
    implicit none
    integer            :: i,j,k
    integer , optional :: byrow
    real(kind=8)       :: mat(:,:),vec(size(mat,dim=1)*size(mat,dim=2)) ; k = 0
    if (present(byrow) .and. byrow.eq.1) then
      do j=1,size(mat,dim=2) ; do i=1,size(mat,dim=1)
        k = k + 1 ; vec(k) = mat(i,j)
      end do ; end do
    else
      do i=1,size(mat,dim=1) ; do j=1,size(mat,dim=2)
        k = k + 1 ; vec(k) = mat(i,j)
      end do ; end do
    end if
    return
  end function vectorize_dp_2d
  function vectorize_dp_3d(mat,byrow) result(vec)
    implicit none
    integer            :: i,j,p,k
    integer , optional :: byrow
    real(kind=8)       :: mat(:,:,:),vec(size(mat,dim=1)*size(mat,dim=2)*size(mat,dim=3)) ; k = 0
    if (present(byrow) .and. byrow.eq.1) then
      do p=1,size(mat,dim=2) ; do j=1,size(mat,dim=2) ; do i=1,size(mat,dim=1)
        k = k + 1 ; vec(k) = mat(i,j,p)
      end do ; end do ; end do
    else
      do i=1,size(mat,dim=1) ; do j=1,size(mat,dim=2) ; do p=1,size(mat,dim=3)
        k = k + 1 ; vec(k) = mat(i,j,p)
      end do ; end do ; end do
    end if
    return
  end function vectorize_dp_3d
  function vectorize_dp_4d(mat,byrow) result(vec)
    implicit none
    integer            :: i,j,p,n,k
    integer , optional :: byrow
    real(kind=8)       :: mat(:,:,:,:),vec(size(mat,dim=1)*size(mat,dim=2)*size(mat,dim=3)*size(mat,dim=4)) ; k = 0
    if (present(byrow) .and. byrow.eq.1) then
      do n=1,size(mat,dim=4) ; do p=1,size(mat,dim=3) ; do j=1,size(mat,dim=2) ; do i=1,size(mat,dim=1)
        k = k + 1 ; vec(k) = mat(i,j,p,n)
      end do ; end do ; end do ; end do
    else
      do i=1,size(mat,dim=1) ; do j=1,size(mat,dim=2) ; do p=1,size(mat,dim=3) ; do n=1,size(mat,dim=4)
        k = k + 1 ; vec(k) = mat(i,j,p,n)
      end do ; end do ; end do ; end do
    end if
    return
  end function vectorize_dp_4d
  function vectorize_dp_5d(mat,byrow) result(vec)
    implicit none
    integer            :: i,j,p,n,k,w
    integer , optional :: byrow
    real(kind=8)       :: mat(:,:,:,:,:)
    real(kind=8)       :: vec(size(mat,dim=1)*size(mat,dim=2)*size(mat,dim=3)*size(mat,dim=4)*size(mat,dim=5)) ; k = 0
    if (present(byrow) .and. byrow.eq.1) then
      do w=1,size(mat,dim=4) ; do n=1,size(mat,dim=4) ; do p=1,size(mat,dim=3) ; do j=1,size(mat,dim=2) ; do i=1,size(mat,dim=1)
        k = k + 1 ; vec(k) = mat(i,j,p,n,w)
      end do ; end do ; end do ; end do ; end do
    else
      do i=1,size(mat,dim=1) ; do j=1,size(mat,dim=2) ; do p=1,size(mat,dim=3) ; do n=1,size(mat,dim=4) ; do w=1,size(mat,dim=5)
        k = k + 1 ; vec(k) = mat(i,j,p,n,w)
      end do ; end do ; end do ; end do ; end do
    end if
    return
  end function vectorize_dp_5d
  function vectorize_int_2d(mat,byrow) result(vec)
    implicit none
    integer            :: i,j,k
    integer , optional :: byrow
    integer            :: mat(:,:),vec(size(mat,dim=1)*size(mat,dim=2)) ; k = 0
    if (present(byrow) .and. byrow.eq.1) then
      do j=1,size(mat,dim=2) ; do i=1,size(mat,dim=1)
        k = k + 1 ; vec(k) = mat(i,j)
      end do ; end do
    else
      do i=1,size(mat,dim=1) ; do j=1,size(mat,dim=2)
        k = k + 1 ; vec(k) = mat(i,j)
      end do ; end do
    end if
    return
  end function vectorize_int_2d
  function vectorize_int_3d(mat,byrow) result(vec)
    implicit none
    integer            :: i,j,p,k
    integer , optional :: byrow
    integer            :: mat(:,:,:),vec(size(mat,dim=1)*size(mat,dim=2)*size(mat,dim=3)) ; k = 0
    if (present(byrow) .and. byrow.eq.1) then
      do p=1,size(mat,dim=3) ; do j=1,size(mat,dim=2) ; do i=1,size(mat,dim=1)
        k = k + 1 ; vec(k) = mat(i,j,p)
      end do ; end do ; end do
    else
      do i=1,size(mat,dim=1) ; do j=1,size(mat,dim=2) ; do p=1,size(mat,dim=3)
        k = k + 1 ; vec(k) = mat(i,j,p)
      end do ; end do ; end do
    end if
    return
  end function vectorize_int_3d
  function vectorize_int_4d(mat,byrow) result(vec)
    implicit none
    integer            :: i,j,p,n,k
    integer , optional :: byrow
    integer            :: mat(:,:,:,:),vec(size(mat,dim=1)*size(mat,dim=2)*size(mat,dim=3)*size(mat,dim=4)) ; k = 0
    if (present(byrow) .and. byrow.eq.1) then
      do n=1,size(mat,dim=4) ; do p=1,size(mat,dim=3) ; do j=1,size(mat,dim=2) ; do i=1,size(mat,dim=1)
        k = k + 1 ; vec(k) = mat(i,j,p,n)
      end do ; end do ; end do ; end do
    else
      do i=1,size(mat,dim=1) ; do j=1,size(mat,dim=2) ; do p=1,size(mat,dim=3) ; do n=1,size(mat,dim=4)
        k = k + 1 ; vec(k) = mat(i,j,p,n)
      end do ; end do ; end do ; end do
    end if
    return
  end function vectorize_int_4d
  function vectorize_int_5d(mat,byrow) result(vec)
    implicit none
    integer            :: i,j,p,n,k,w
    integer , optional :: byrow
    integer            :: mat(:,:,:,:,:)
    integer            :: vec(size(mat,dim=1)*size(mat,dim=2)*size(mat,dim=3)*size(mat,dim=4)*size(mat,dim=5)) ; k = 0
    if (present(byrow) .and. byrow.eq.1) then
      do w=1,size(mat,dim=4) ; do n=1,size(mat,dim=4) ; do p=1,size(mat,dim=3) ; do j=1,size(mat,dim=2) ; do i=1,size(mat,dim=1)
        k = k + 1 ; vec(k) = mat(i,j,p,n,w)
      end do ; end do ; end do ; end do ; end do
    else
      do i=1,size(mat,dim=1) ; do j=1,size(mat,dim=2) ; do p=1,size(mat,dim=3) ; do n=1,size(mat,dim=4) ; do w=1,size(mat,dim=5)
        k = k + 1 ; vec(k) = mat(i,j,p,n,w)
      end do ; end do ; end do ; end do ; end do
    end if
    return
  end function vectorize_int_5d

  ! ----------------------------------------------------------------------------
  ! this function returns a vector "vec1" with the cummmulative sum of the elements
  ! of a vector "vec0"

  function cumsum(vec0) result(vec1)
    implicit none
    real(kind=8) :: vec0(:),vec1(size(vec0))
    integer      :: i
    vec1(1) = vec0(1) ; do i=2,size(vec0) ; vec1(i) = vec1(i-1) + vec0(i) ; end do
    return
  end function cumsum

  ! ----------------------------------------------------------------------------
  ! this function returns the main diagonal of a matric "mat"

  function diag(mat) result(vec)
    implicit none
    real(kind=8) :: mat(:,:),vec(size(mat,dim=1))
    integer      :: i
    if (size(mat,1).ne.size(mat,2)) call error(' error in diag: matrix not sqaure')
    do i=1,size(mat,dim=1) ; vec(i) = mat(i,i) ; end do
    return
  end function diag

  ! ----------------------------------------------------------------------------
  ! returns the transpose of a matrix "mat"

  function transmat(mat) result(matt)
    implicit none
    integer                    :: i,j
    real(kind=8) , intent(in)  :: mat(:,:)
    real(kind=8)               :: matt(size(mat,2),size(mat,1))
    do i=1,size(mat,1) ; do j=1,size(mat,2)
      matt(j,i) = mat(i,j)
    end do ; end do
    return
  end function transmat

  ! ----------------------------------------------------------------------------
  ! compute the inverse of a sqaure matrix "m"

  function inverse(m) result(im)
    implicit none
    integer      :: n,i,j,k
    real(kind=8) :: m(:,:),coeff
    real(kind=8) :: im(size(m,1),size(m,1)),b(size(m,1)),d(size(m,1)),x(size(m,1))
    real(kind=8) :: l(size(m,1),size(m,1)),u(size(m,1),size(m,1)),mb(size(m,1),size(m,1))
    if (size(m,1).ne.size(m,2)) call error(' error in inverse: matrix is not square')
    n  = size(m,1)
    l  = 0.0
    u  = 0.0
    b  = 0.0
    mb = m
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
  !  - golden search: single-value function with one unknown
  !  - brent method: single-value function with one unknown
  !  - simplex method: single-value function with "n" unknowns
  !  - levenberg–marquardt algorithm: system of "m" functions with "n" unknowns
  !
  ! the four algorithms can be used in problems with states (variables over which
  ! the algorithm does not optimize) which I found useful in many applications.
  ! for example, maximize a function f(x,y) over "x" for a given value of "y" (state)
  ! since I typically use this feature with state variables, the vector "states" is
  ! a vector of integers that indicates the position of each state variable on 
  ! its corresponding grid.
  !
  ! if your problem has states, please check the corresponding subroutine to make
  ! sure your function has the correct form.
  !
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



  ! ----------------------------------------------------------------------------
  ! golden search algorithm
  ! finds the maximum of a single-valued function with one unknown
  ! the user should make sure the maximum of "func" is between "xmin" and "xmax"

  subroutine golden(func,x,y,xmax,xmin,itermax,tol)

    ! this subroutine finds the maximum of a single-valued function with one unknown
    ! using the golden search algorithm
    !
    ! the inputs are:
    !
    !   - func: user-provided function to be minimize. this function should of the form:
    !
    !       function func(x) result(y)
    !         real(kind=8) :: x(:),y
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
    !   x:    the value of x that maximizes "func"
    !   y:    the value of the function at the point x

    implicit none
    external                             :: func
    real(kind=8) , intent(in)            :: xmax,xmin
    real(kind=8) , intent(out)           :: x,y
    real(kind=8) , intent(in) , optional :: tol
    integer      , intent(in) , optional :: itermax
    real(kind=8) , parameter             :: alpha=0.61803399
    real(kind=8)                         :: x0,x1,x2,x3,f1,f2,tolgold
    integer                              :: it,maxiter

    interface
      function func(xx) result(ff)
      real(kind=8) :: xx,ff
      end function func
    end interface

    tolgold = 1.0d-8 ; if (present(tol)    ) tolgold = tol
    maxiter = 500    ; if (present(itermax)) maxiter = itermax

    x0 = xmin
    x3 = xmax
    x1 = alpha*x0 + (one-alpha)*x3
    x2 = alpha*x3 + (one-alpha)*x1 
    f1 = func(x1)
    f2 = func(x2)

    it = 0
    do while (abs(x0-x3).gt.tolgold*(abs(x2)+abs(x1)) .and. it.lt.maxiter)
      it = it + 1
      if (f2.gt.f1) then
        x0 = x1
        x1 = x2
        f1 = f2
        x2 = alpha*x1 + (one-alpha)*x3
        f2 = func(x2)
      else
        x3 = x2
        x2 = x1
        f2 = f1
        x1 = alpha*x2 + (one-alpha)*x0
        f1 = func(x1)
      end if
    end do

    if (f1.gt.f2) then
      y = f1
      x = x1
    else
      y = f2
      x = x2
    end if

    return
  end subroutine golden

  ! ----------------------------------------------------------------------------
  ! brent's method
  ! hybrid root-finding algorithm combining the bisection and the secant methods

  subroutine brent(func,x,iy,ind,x0,x1,itermax,tol)

    ! this subroutine finds the root of a user-supplied single-valued function
    ! with one unknown using the brent's method
    !
    ! the inputs are:
    !
    !   - func: user-provided function to be minimize. this function should of the form:
    !
    !       function func(x) result(y)
    !         real(kind=8) :: x(:),y
    !       end function func
    !
    !   - x0: lower bound of variable "x" such that func(x0)<0
    !   - x1: upper bound of variable "x" such that func(x1)>0
    !
    ! there are also few optional inputs:
    !
    !   itermax:  max number of iterations (default = 500)
    !   tol:      tolerance level (default = 1.0d-8)
    !
    ! the outputs of the subroutine are:
    !
    !   x:    the value of x that maximizes "func"
    !   iy:   number of function evaluations required
    !   ind:  exit indicator, taking value:
    !           0: solved, func(x) = 0
    !           1: not solved

    implicit none
    external                             :: func
    real(kind=8) , intent(out)           :: x
    integer      , intent(out)           :: iy
    integer      , intent(out)           :: ind
    real(kind=8) , intent(in)            :: x0,x1
    real(kind=8) , intent(in) , optional :: tol
    integer      , intent(in) , optional :: itermax
    integer                              :: maxiter
    real(kind=8)                         :: toler,xa,ya,xb,yb,xc,yc,ys,xs

    interface
      function func(xx) result(resid)
        real(kind=8) :: xx,resid
      end function func
    end interface

    toler   = 1.0d-8 ; if (present(tol)    ) toler   = tol
    maxiter = 500    ; if (present(itermax)) maxiter = itermax

    xa = x0 ; ya = func(xa)
    xb = x1 ; yb = func(xb)

    if (ya*yb.gt.zero) then
      ind = 1 ; x = xa ; return
    end if
    if (abs(ya).lt.abs(yb)) then
      xc = xb ; xb = xa ; xa = x0 ; yc = ya ; ya = yb
    end if
    if (abs(yb).lt.toler) then
      ind = 0 ; x = xb ; return
    end if

    xc = xa ; iy = 0 ; ind = 1

    do while (abs(xb-xa).gt.toler .and. iy.lt.maxiter)
      yc = func(xc)
      if (abs(yc-ya).gt.toler .and. abs(yc-yb).gt.toler) then
        xs = xa*yb*yc/((ya-yb)*(yc-yc)) + xb*ya*yc/((yb-ya)*(yb-yc)) + xc*ya*yb/((yc-ya)*(yc-yb))
      else
        xs = xb - yb*(xb-xa)/(yb-ya)
      end if
      ys = func(xs) ; iy = iy + 1
      if (abs(ys).lt.toler) then
        ind = 0 ; x = xs ; return
      end if
      if (ya*ys.lt.zero) then
        xc = xb ; xb = xs ; yb = ys
      else
        xa = xs ; ya = ys
      end if
      if (abs(ya).lt.abs(yb)) then
        xc = xb ; xb = xa ; xa = xc ; yc = yb ; yb = ya ; ya = yc
      end if
    end do

    x = xb

    return
  end subroutine brent

  ! ----------------------------------------------------------------------------
  ! simplex algorithm
  ! minimize a user-supplied single-valued function in "n" unknows

  subroutine simplex(func,x,y,iy,ind,x0,itermax,tol,iprint)

    ! this subroutine minimizes a single-valued function with n unknowns using the
    ! simplex algorithm
    !
    ! the inputs are:
    !
    !   - func: user-provided function to be minimize. this function should of the form:
    !
    !       function func(x,states) result(y)
    !         real(kind=8) :: x(:),y
    !       end function func
    !
    !   - x0:     the vector of dimension n with the initial guess
    !
    ! there are also few optional inputs:
    !
    !   itermax:  max number of iterations (default = 500)
    !   tol:      tolerance level (default = 1.0d-8)
    !   iprint:   indicator for printing results (default = 0)
    !             0: print no results
    !             1: print main results
    !             2: print all the iterations
    !
    ! the outputs of the subroutine are:
    !
    !   x:    the vector of x that solves the system
    !   y:    the value of the function at the point x
    !   iy:   the number of total function evaluations
    !   ind:  exit indicator, taking value:
    !           0: system has converged, func(x) = 0
    !           1: simplex is too small
    !           9: max number of iterations reached

    implicit none

    external func
    integer      , intent(out)           :: iy,ind
    real(kind=8) , intent(out)           :: x(:),y

    real(kind=8) , intent(in)            :: x0(:)
    real(kind=8) , intent(in) , optional :: tol
    integer      , intent(in) , optional :: itermax
    integer      , intent(in) , optional :: iprint

    integer                  :: n,i,j,ilow,ihigh,ihigh2,maxiter,ipri
    real(kind=8)             :: xp(size(x,1),size(x,1)+1),yp(size(x,1)+1),xr(size(x,1))
    real(kind=8)             :: xm(size(x,1)),xe(size(x,1)),xc(size(x,1))
    real(kind=8)             :: y0,yr,ye,yc,cent,toler
    real(kind=8) , parameter :: alpha=1.0,beta=0.5,gamma=2.0

    interface
      function func(xvar) result(resid)
        real(kind=8) :: xvar(:),resid
      end function func
    end interface

    y0 = 0.0d-10 ; yp = 0.0d-10 ; ye = 0.0d-10 ; yc = 0.0d-10 ; yr = 0.0d-10
    xm = 0.0d-10 ; xp = 0.0d-10 ; xe = 0.0d-10 ; xc = 0.0d-10 ; xr = 0.0d-10

    toler   = 1.0d-8 ; if (present(tol)    ) toler   = tol
    maxiter = 500    ; if (present(itermax)) maxiter = itermax
    ipri    = 0      ; if (present(iprint) ) ipri    = iprint

    if (ipri.ge.1) then
      write(*,*) '  '
      write(*,*) ' starting simplex algorithm'
      write(*,*) '  '
    end if

    ind = -9 ; ilow = 1 ; n = size(x,1) ; y0 = func(x0) ; iy = 1 ; yp(1) = y0 ; xp(:,1) = x0

    if (ipri.ge.1) print 93 , y0

    if (y0.lt.toler) then
      yp(ilow) = y0 ; xp(:,ilow) = x0 ; goto 10
    end if

    if (ipri.ge.1) then
      write(*,*) '  '
      write(*,*) ' evaluating initial simplex'
      write(*,*) '  '
    end if
    
    xp = zero
    do j=1,n+1
      xp(:,j) = x0(:)
    end do

    do i=2,n+1
      if (x0(i-1).lt.zero) xp(i-1,i) = x0(i-1) - dble(0.20)
      if (x0(i-1).ge.zero) xp(i-1,i) = x0(i-1) + dble(0.20)
    end do

    do i=2,n+1
      yp(i) = func(xp(:,i))
      if (ipri.gt.1 .and. n.lt.5) write(*,97) i-1, yp(i) , xp(:,i)
      if (ipri.gt.1 .and. n.ge.5) write(*,97) i-1, yp(i)
      if (yp(i).lt.toler) then
        yp(ilow) = yp(i) ; xp(:,ilow) = xp(:,i) ; goto 10
      end if
    end do

    if (ipri.gt.1) then
      write(*,*) '  '
      write(*,*) ' beginning iterations'
      write(*,*) '  '
    end if

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
    if (ipri.gt.1) write(*,97) iy, yp(ilow)
    
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
      xr(j) = (one+alpha)*xm(j) - alpha*xp(j,ihigh)
    end do
    yr = func(xr) ; iy = iy + 1

    ! check convergence
    if (yr.lt.toler) then
      yp(ilow) = yr ; xp(:,ilow) = xr(:) ; goto 10
    end if

    ! check number of function evaluations
    if (iy.ge.maxiter) goto 19

    if (yr.lt.yp(ilow)) then

      ! new point
      do j=1,n
        xe(j) = (one-gamma)*xm(j) - gamma*xp(j,ihigh)
      end do
      ye = func(xe) ; iy = iy + 1

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
      if (iy.ge.maxiter) goto 19

    elseif (yr.gt.yp(ihigh2)) then

      ! if "yr" is worse, update worse point
      if (yr.lt.yp(ihigh)) then
        xp(:,ihigh) = xr(:) ; yp(ihigh) = yr
      end if

      ! new point
      do j=1,n
        xc(j) = (one-beta)*xm(j) - beta*xp(j,ihigh)
      end do
      yc = func(xc) ; iy = iy + 1

      ! check convergence
      if (yc.lt.toler) then
        xp(:,ilow) = xc ; yp(ilow) = yc ; goto 10
      end if

      ! check number of function evaluations
      if (iy.ge.maxiter) goto 19
      
      ! if "yc" is worse, update worse point
      if (yc.lt.yp(ihigh)) then

        xp(:,ihigh) = xc(:) ; yp(ihigh) = yc

      else
        do i=1,n+1
          if (i.ne.ilow) then

            ! new point
            do j=1,n
              xp(j,i) = half*( xp(j,i) + xp(j,ilow) )
            end do
            yp(i) = func(xp(:,i)) ; iy = iy + 1

            ! check convergence
            if (yp(i).lt.toler) then
              xp(:,ilow) = xp(:,i) ; yp(ilow) = yp(i) ; goto 10
            end if

            ! check number of function evaluations
            if (iy.ge.maxiter) goto 19

          end if
        end do
      end if

    else

      xp(:,ihigh) = xr(:)
      yp(ihigh)   = yr

    end if

    goto 1

    10 ind = 0 ; goto 7
    11 ind = 1 ; goto 7
    19 ind = 9 ; goto 7

    7 continue

    ilow = 1
    do i=1,n+1
      if (yp(i).lt.yp(ilow)) ilow = i
    end do
    
    if (ipri.ge.1) then
      write(*,*) '  '
      if (ind.eq.0) write(*,*) ' simplex finished: system has converged, func=0 '
      if (ind.eq.1) write(*,*) ' simplex finished: system converged, simplex is very small'
      if (ind.eq.9) write(*,*) ' simplex finished: max number of iterations'
      write (*, *) '  '
      write (*,94) iy
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

  subroutine lmmin(func,x,y,iy,ind,x0,itermax,damp,tol,toleach,shock,usebro,iprint)

    ! this subroutine minimizes the sum of squred errors of a system of m equations
    ! in n unknowns using the levenberg–marquardt algorithm. more about the algorithm
    ! can be found in https://en.wikipedia.org/wiki/levenberg–marquardt_algorithm.
    ! this method combines a newton-based algorithm and the gradient descent method.
    !
    ! the inputs are:
    !
    !   - func: user-provided function to be minimize. this function should of the form:
    !
    !       function func(x) result(y)
    !         real(kind=8) :: x(:)
    !         real(kind=8),allocatable :: y(:)
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
    !             0: print no results
    !             1: print main results
    !             2: print all the iterations
    !
    ! the outputs of the subroutine are:
    !
    !   x:    the vector of x that solves the system
    !   y:    the value of the function at the point x
    !   iy:   the number of total function evaluations
    !   ind:  exit indicator, taking value:
    !           0: system has converged, func(x) = 0
    !           1: jacobian is closed to 0
    !           2: step in x is close to 0
    !           9: max number of iterations reached

    implicit none

    external                             :: func
    
    ! outputs
    real(kind=8) , intent(out)           :: x(:),y(:)
    integer      , intent(out)           :: iy,ind

    ! inputs (including optionals)
    real(kind=8) , intent(in)            :: x0(:)
    real(kind=8) , intent(in) , optional :: shock
    real(kind=8)                         :: shck(size(x,1))
    real(kind=8) , intent(in) , optional :: damp
    real(kind=8)                         :: da
    real(kind=8) , intent(in) , optional :: tol
    real(kind=8)                         :: toler
    real(kind=8) , intent(in) , optional :: toleach
    real(kind=8)                         :: tolereach
    integer      , intent(in) , optional :: itermax
    integer                              :: maxiter
    integer      , intent(in) , optional :: iprint
    integer                              :: ip
    integer      , intent(in) , optional :: usebro
    integer                              :: bro,br

    ! other variables
    integer      ::  i,k,n,m,qp
    real(kind=8) ::  j(size(y,1),size(x,1))
    real(kind=8) :: j0(size(y,1),size(x,1))
    real(kind=8) :: ij(size(x,1),size(x,1))
    real(kind=8) :: jj(size(x,1),size(x,1))
    real(kind=8) :: ja(size(x,1),size(x,1))
    real(kind=8) :: jt(size(x,1))
    real(kind=8) :: e0,y0(size(y,1))
    real(kind=8) :: eb,yb(size(y,1)),xb(size(x,1))
    real(kind=8) :: e1,y1(size(y,1)),x1(size(x,1))
    real(kind=8) :: ea,ya(size(y,1)),xa(size(x,1))
    real(kind=8) :: ej,yj(size(y,1),size(x,1)),xj(size(x,1),size(x,1))
    real(kind=8) :: dy(size(y,1)),dx(size(x,1)),ddx

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
    tolereach = 1.0d-8 ; if (present(toleach)) tolereach = toleach  ! tolerance level for each equation
    shck(:)   = 5.0d-2 ; if (present(shock)  ) shck(:)   = shock    ! shock to parameters
    maxiter   = 500    ; if (present(itermax)) maxiter   = itermax  ! max number of iterations
    ip        = 0      ; if (present(iprint) ) ip        = iprint   ! indicator of printing behaviour
    bro = 0 ; br = 0   ; if (present(usebro) ) bro       = usebro   ! indicator for the use of Broyden Updating

    if (ip.ge.1) then
      write(*,'(a)'   ) '                              '
      write(*,'(a)'   ) '  starting lmmin algorithm    '
      write(*,'(a)'   ) '                              '
      write(*,'(a,i2)') '  # equations = ' , size(y,1)
      write(*,'(a,i2)') '  # unknows   = ' , size(x,1)
      write(*,'(a)'   ) '                              '
    end if

    ! number of variables ad number of equations
    n = size(x,1) ; m = size(y,1)

    ! initial point
    y0 = func(x0) ; e0 = sum(y0(:)*y0(:)) ; iy = 1

    ! best point so far = initial point
    yb = y0 ; xb = x0 ; eb = e0

    ! if converged, finish
    if (e0.lt.toler .or. maxval(abs(y0)).lt.tolereach) goto 10

    ! if reached max number of iterations, finish
    if (iy.ge.maxiter) goto 19

    ! print results
    if ( ip.ge.1 ) write (*,91) eb
    if ( ip.eq.1 ) write (*, *) '  '

    qp = 0    ! this variable takes value 0 if no good point has een found with current jacobian, and 1 otherwise
    br = 0    ! this variable takes value 0 if new jacobian should be computed numerical (not use broyden)

    1 continue

    ! if converged, finish
    if (eb.lt.toler .or. maxval(abs(yb)).lt.tolereach) goto 10

    ! update numerical jacobian using broyden's method
    if ( bro.eq.1 .and. br.eq.1 .and. iy.gt.1 ) then

      j0 = j ; dx(:) = x1(:)-xa(:) ; dy(:) = y1(:)-ya(:)
      do k = 1,m ; do i = 1,n
        j(k,i) = j0(k,i) + dx(i)*( dy(k) - sum(j0(k,:)*dx(:)) )/sum(dx(:)*dx(:))
      end do ; end do

     ! compute new numerical jacobian
    else ; if (ip.gt.1) write (*,*) '  '

      do i = 1,n

        ! shock the i-th parameter
        xj(:,i) = xb(:)
        if (x(i).le.zero) then
          xj(i,i) = xb(i) - min(shck(i),shck(i)*max(0.30,eb))*abs(xb(i))
        elseif (x(i).ge.zero) then
          xj(i,i) = xb(i) + min(shck(i),shck(i)*max(0.30,eb))*abs(xb(i))
        end if
        if (abs(xj(i,i)-xb(i)).lt.toler*dble(10.0)) then
          xj(i,i) = xb(i) + toler*dble(10.0)
        end if

        ! evaluate the function
        yj(:,i) = func(xj(:,i))
        ej = sum(yj(:,i)*yj(:,i))
        iy = iy + 1

        ! print result
        if (ip.gt.1) write (*,90) i,ej

        ! if converged, finish
        if (ej.lt.toler .or. maxval(abs(yj(:,i))).lt.tolereach) then
          yb = yj(:,i) ; xb = xj(:,i) ; eb = ej ; goto 10
        end if

        ! if reached max number of iterations, finish
        if (iy.ge.maxiter) goto 19

        ! if change in func is too large, compute new shocked parameter
        if (ej.gt.dble(50.0)*eb) then
          do while (ej.gt.dble(50.0)*eb .and. abs(xj(i,i)-xb(i)).gt.toler)
            xj(i,i) = (xb(i) + xj(i,i))/dble(2.0)
            yj(:,i) = func(xj(:,i)) ; ej = sum(yj(:,i)*yj(:,i)) ; iy = iy + 1
            if (ip.gt.1) write (*,90) i,ej
            if (iy.ge.maxiter) goto 19
          end do
          shck(i) = abs( ( xb(i) - xj(i,i) )/abs(xb(i)) )
        end if

        ! if change in func is too small, compute new shocked parameter
        if (abs(ej-eb).lt.toler) then ; ddx = xj(i,i) - xb(i)
          do while (abs(ej-eb).lt.toler )
            xj(i,i) = xj(i,i) + ddx
            yj(:,i) = func(xj(:,i)) ; ej = sum(yj(:,i)*yj(:,i)) ; iy = iy + 1
            if (ip.gt.1) write (*,90) i,ej
            if (iy.ge.maxiter) goto 19
          end do
          shck(i) = abs( ( xb(i) - xj(i,i) )/abs(xb(i)) )
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
    do i=1,n
      do k=1,n
        jj(i,k) = sum(j(:,i)*j(:,k))
      end do
      jt(i) = sum(j(:,i)*yb(:))
    end do
    if (maxval(abs(jt)).lt.toler) goto 11
    ja = jj
    do i=1,n
      ja(i,i) = (one+da)*jj(i,i)
    end do
    ij = inverse(ja)
    do i=1,n
      x1(i) = xb(i) - sum(ij(i,:)*jt(:))
    end do
    if (sum(abs(x1-xb)).le.toler) goto 12

    ! evaluate new point
    y1 = func(x1) ; e1 = sum(y1(:)*y1(:)) ; iy = iy + 1

    ! print iterations
    if (ip.ge.1) print 97 , iy , e1 , eb , da

    ! check number of iterations
    if (iy.ge.maxiter) goto 19

    ! if the point improves, reduce dumping factor and new jacobian
    if (e1.lt.eb) then
      yb = y1 ; xb = x1 ; eb = e1   ! if good point, update best point
      da = da/dble(2.0)             ! if good point, decrease the damping factor
      qp = 1                        ! if good point, set parameter qp to 1
      goto 1                        ! if good point, update jacobian (with broyden if first point, numerical if not)
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

    10 ind = 0 ; goto 7
    11 ind = 1 ; goto 7
    12 ind = 2 ; goto 7
    19 ind = 9 ; goto 7

    7 if (ip.ge.1) then
      write(*,*) '  '
      if (ind.eq.0) write(*,*) ' lmmin finished: system has converged, func=0 '
      if (ind.eq.1) write(*,*) ' lmmin finished: system has converged, jac close to 0'
      if (ind.eq.2) write(*,*) ' lmmin finished: system has converged, step in x close to 0'
      if (ind.eq.9) write(*,*) ' lmmin finished: max number of iterations'
      write(*,*          ) '   '
      write(*,*          ) ' iterations = ', iy
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
    real(kind=8) , intent(in)  :: xmax,xmin,x
    real(kind=8) , intent(out) :: y
    if (xmax.lt.xmin) call error(' errror in normalize: xmax > xmin')
    if (x.gt.xmax   ) call error(' errror in normalize: x > xmax')
    if (x.lt.xmin   ) call error(' errror in normalize: x < xmin')
    y = log((x-xmin)/(xmax-xmin)/max(0.0001,one-(x-xmin)/(xmax-xmin)))
    return
  end subroutine normalize

  ! ----------------------------------------------------------------------------
  ! this subroutine applies the transofrmation 
  !   x = xmin + ( exp(y)/(1+exp(y) )*(xmax-xmin) 
  ! to return a bounded variable "x" from an unbounded variable "y"

  subroutine denormalize(y,x,xmax,xmin)
    implicit none
    real(kind=8) , intent(in)  :: xmax,xmin,y
    real(kind=8) , intent(out) :: x
    if (xmax.lt.xmin) call error(' errror in denormalize: xmax > xmin')
    x = (exp(y)/(one+exp(y)))*(xmax-xmin) + xmin
    if (isnan(x) .and. y.gt.cien) x = xmax
    return
  end subroutine denormalize

  ! ----------------------------------------------------------------------------
  ! this subroutines updates a jacobian matrix using the broyden's method.
  ! it is used by the levenberg–marquardt algorithm.

  subroutine broyden(j1,j0,x1,x0,f1,f0)
    implicit none
    real(kind=8) , intent(in)  :: x1(:),f1(:),x0(:),f0(:),j0(:,:)
    real(kind=8) , intent(out) :: j1(:,:)
    real(kind=8)               :: df(size(f0,1))
    real(kind=8)               :: dx(size(x1,1))
    integer                    :: i,k
    dx(:) = x1(:)-x0(:)
    df(:) = f1(:)-f0(:)
    do i=1,size(x1) ; do k=1,size(f1)
      j1(i,k) = j0(i,k) + dx(k)*( df(i)-sum(j0(i,:)*dx(:)) )/sum(dx(:)*dx(:))
    end do ; end do
    return
  end subroutine broyden

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module toolkit
