
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! toolkit.f90, a toolkit for fortran90 programming
! Borja Petit, CUNEF Universidad, © 2025
!
! The code is distributed under the MIT license.
! You are free to use and modify it as you wish acknowledging the source.
! If you modified version of this code, indicate the changes made along with your name and date.
!
! UPDATES:
! - October 2025 – Borja Petit:
!      * the functions "ces" and "crra" werre added
!      * some minor typos and inconsistencies were fixed. 
!      * modified the "error" call in different subroutines/functions to stop the execution when an error is found.
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! LIST OF FUNCTIONS AND SUBROUTINES INCLUDED IN THIS CODE:
! *******************************************************
!
! econ-related:
!   - crra: CRRA function
!   - ces: CES function
!
! general purpose:
!   - grid: generate a grid for a continuous varibale
!   - interpolation: interpolate a value over a grid, returning position and distance
!   - interpolate: linearly interpolate a value over an n-dimensional grid, with n <= 6
!   - timing: returns the number of seconds since 00:00h of the 1st day of the month [robust to parelalization]
!   - multiplo: returns 1 if an integer is a multiple of another user-provided integer
!   - iseven: returns 1 if a user-provided integer is even
!   - error: prints an error message
!   - num2text: returns a string with an integer or a real number
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

  integer  , parameter :: dp    = kind(1.0d00)
  real(dp) , parameter :: cero  = dble(0.00000000000000)
  real(dp) , parameter :: medio = dble(0.50000000000000)
  real(dp) , parameter :: uno   = dble(1.00000000000000)
  real(dp) , parameter :: dos   = dble(2.00000000000000)
  real(dp) , parameter :: cien  = dble(100.000000000000)
  real(dp) , parameter :: mil   = dble(1000.00000000000)
  real(dp) , parameter :: tolvl = dble(0.00000000010000)
  
  interface interpolate
    module procedure interpolate1d,&
                     interpolate2d,&
                     interpolate3d,&
                     interpolate4d,&
                     interpolate5d,&
                     interpolate6d
  end interface interpolate
  interface vect
    module procedure vectorize_int_2d,vectorize_dp_2d,&
                     vectorize_int_3d,vectorize_dp_3d,&
                     vectorize_int_4d,vectorize_dp_4d,&
                     vectorize_int_5d,vectorize_dp_5d
  end interface vect
  interface randomnormal
    module procedure randomnormal_scalar,randomnormal_vec
  end interface randomnormal
  interface error
    module procedure error_0,error_1
  end interface error
  interface num2text
    module procedure int2text,real2text
  end interface num2text
  interface shuffle_vect
    module procedure shuffle_vect_int,shuffle_vect_dp
  end interface shuffle_vect
  
  contains



  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! econ functions and subroutines
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   ! crra function with "x" being the argument and "b" being the rra parameter
  function crra(x,b) result(u)
    implicit none
    real(dp) :: x,b,u ; u = cero
    if (abs(b).gt.tolvl) then
      if (abs(x).gt.tolvl) u = (x**b)/b
      if (abs(x).le.tolvl) then
        if (b.gt.cero) u = cero
        if (b.lt.cero) u = -huge(x)
      end if
    else if (abs(b).le.tolvl) then
      if (x.gt.tolvl) u = log(x)
      if (x.le.tolvl) u = -huge(x)
    end if
    return
  end function crra

   ! CES function with "x1" and "x2" being the inputs, "b" being (inverse) elasticity of substitution and "a" being the share of input 1
  function ces(x1,x2,a,b) result(u)
    implicit none
    real(dp) :: x1,x2,a,b,u ; u = cero
    if (abs(b).gt.tolvl) then
      u = (a*(x1**b) + (uno-a)*(x2**b))**(uno/b)
    else if (abs(b).le.tolvl) then
      u = (x1**a)*(x2**(uno-a))
    end if
    return
  end function ces 


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
    if (ss.le.cero) call error(' errror in grid: spacing parameter is nonpositive',0)
    if (n.lt.2    ) call error(' errror in grid: required dimension of 1',0)
    xmin = min(maxv,minv)
    xmax = max(maxv,minv)
    if (ss.gt.cero) forall (i=1:n) grid0(i) = dble(i-1)/dble(n-1)
    if (ss.gt.cero) forall (i=1:n) v(i) = (grid0(i)**ss)*(xmax-xmin) + xmin
    return
  end function grid

  ! ----------------------------------------------------------------------------
  ! this subroutine finds the closest point in a given grid and return its
  ! position and the relative distance

  subroutine interpolation(pos,wth,xnow,xgrid,outof)
    implicit none
    integer                         :: j,n,ofs
    real(dp), intent(in)            :: xnow,xgrid(:)
    real(dp), intent(out)           :: wth
    integer , intent(out)           :: pos
    integer , intent(in) , optional :: outof

    pos = 0 ; wth = cero ; ofs = 0

    if (present(outof)) ofs = outof
    
    if (isnan(xnow)) then
      call error(' errror in interpolation: xnow is nan',0)
      return
    end if
    
    n = size(xgrid)
    if (n.eq.1) then
      pos = 1
      wth = uno
    else
      if (xnow.le.xgrid(1)) then
        pos = 2
        wth = cero
      else if (xnow.ge.xgrid(n)) then
        pos = n
        wth = uno
      else if (xnow.gt.xgrid(1) .and. xnow.lt.xgrid(n)) then
        j = 2
        do while (xnow.gt.xgrid(j))
          j = j + 1
        end do
        pos = j
        wth = (xnow-xgrid(j-1))/(xgrid(j)-xgrid(j-1))
      end if
    end if

    if (pos.lt.1   ) call error(' errror in interpolation: pos < 1',0)
    if (pos.gt.n   ) call error(' errror in interpolation: pos > size',0)
    if (isnan(wth) ) call error(' errror in interpolation: wvar is nan',0)
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
    real(dp) :: y1(:),m(:),x1,xi,wth
    if (size(y1).ne.size(m)) then
      call error(' error in interpolate: 1st dimension incorrect',0)
      return
    end if
    if (size(y1).eq.1      ) then
      call error(' error in interpolate: 1st dimension with single element',0)
      return
    end if
    call interpolation(pos,wth,x1,y1)
    xi = m(pos)*wth + m(pos-1)*(uno-wth)
    return
  end function interpolate1d
  function interpolate2d(x1,x2,y1,y2,m) result(xi)
    implicit none
    integer      :: pos2
    real(dp) :: wth2
    real(dp) :: x1,x2,xi
    real(dp) :: y1(:),y2(:),m(:,:)
    if (size(m,2).ne.size(y2)) then
      call error(' error in interpolate: 2nd dimension incorrect',0)
      return
    end if
    if (size(y2).eq.1        ) then
      call error(' error in interpolate: 2st dimension with single element',0)
      return
    end if
    call interpolation(pos2,wth2,x2,y2)
    xi = interpolate1d(x1,y1,m(:,pos2)  )*wth2 + &
         interpolate1d(x1,y1,m(:,pos2-1))*(uno-wth2)
    return
  end function interpolate2d
  function interpolate3d(x1,x2,x3,y1,y2,y3,m) result(xi)
    implicit none
    integer      :: pos3
    real(dp) :: y1(:),y2(:),y3(:),x1,x2,x3,xi
    real(dp) :: wth3,m(:,:,:)
    if (size(m,3).ne.size(y3)) then
      call error(' error in interpolate: 3rd dimension incorrect',0)
      return
    end if
    if (size(y3).eq.1        ) then
      call error(' error in interpolate: 3rd dimension with single element',0)
      return
    end if
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
    if (size(m,4).ne.size(y4)) then
      call error(' error in interpolate: 4th dimension incorrect',0)
      return
    end if
    if (size(y4).eq.1        ) then
      call error(' error in interpolate: 4th dimension with single element',0)
      return
    end if
    call interpolation(pos4,wth4,x4,y4)
    xi = interpolate3d(x1,x2,x3,y1,y2,y3,m(:,:,:,pos4)  )*wth4 + &
         interpolate3d(x1,x2,x3,y1,y2,y3,m(:,:,:,pos4-1))*(uno-wth4)
    return
  end function interpolate4d
  function interpolate5d(x1,x2,x3,x4,x5,y1,y2,y3,y4,y5,m) result(xi)
    implicit none
    integer      :: pos5
    real(dp) :: y1(:),y2(:),y3(:),y4(:),y5(:),x1,x2,x3,x4,x5,xi
    real(dp) :: wth5,m(:,:,:,:,:)
    if (size(m,5).ne.size(y5)) then
      call error(' error in interpolate: 5th dimension incorrect',0)
      return
    end if
    if (size(y5).eq.1        ) then
      call error(' error in interpolate: 5th dimension with single element',0)
      return
    end if
    call interpolation(pos5,wth5,x5,y5)
    xi = interpolate4d(x1,x2,x3,x4,y1,y2,y3,y4,m(:,:,:,:,pos5))*wth5 + &
         interpolate4d(x1,x2,x3,x4,y1,y2,y3,y4,m(:,:,:,:,pos5-1))*(uno-wth5)
    return
  end function interpolate5d
  function interpolate6d(x1,x2,x3,x4,x5,x6,y1,y2,y3,y4,y5,y6,m) result(xi)
    implicit none
    integer      :: pos6
    real(dp) :: y1(:),y2(:),y3(:),y4(:),y5(:),y6(:),x1,x2,x3,x4,x5,x6,xi
    real(dp) :: wth6,m(:,:,:,:,:,:)
    if (size(m,6).ne.size(y6)) then
      call error(' error in interpolate: 6th dimension incorrect',0)
      return
    end if
    if (size(y6).eq.1        ) then
      call error(' error in interpolate: 6th dimension with single element',0)
      return
    end if
    call interpolation(pos6,wth6,x6,y6)
    xi = interpolate5d(x1,x2,x3,x4,x5,y1,y2,y3,y4,y5,m(:,:,:,:,:,pos6))*wth6 + &
         interpolate5d(x1,x2,x3,x4,x5,y1,y2,y3,y4,y5,m(:,:,:,:,:,pos6-1))*(uno-wth6)
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
    real(dp)           :: mod,time,v2(8) ; time = cero ; v2 = cero
    
    mod = 1 ; if (present(mode)) mod = mode

    call date_and_time(values=v1)  

    v2(2) = dble(v1(2)*30*24*60*60) ! months
    v2(3) = dble(v1(3)*   24*60*60) ! days
    v2(5) = dble(v1(5)*      60*60) ! hours
    v2(6) = dble(v1(6)*         60) ! minutes
    v2(7) = dble(v1(7)            ) ! seconds

    ! measured in hours
    if (mod.eq.3) then 
      time = sum(v2)/dble(60*60)
    ! measured in minutes
    else if (mod.eq.2) then 
      time = sum(v2)/dble(60)
    ! measured in seconds
    else if (mod.eq.1) then 
      time = sum(v2)
    end if

    return
  end function timing

  ! ----------------------------------------------------------------------------
  ! this function checks whether a number "num" is a multiple of "xx".

  function multiplo(num,xx) result(mul)
    implicit none
    integer :: num,xx
    logical :: mul ; mul = .false.
    if (abs(xx).lt.tolvl) then
      call error(' error in multiplo: zero has no multiples',0)
      return
    end if
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
  ! this subroutine prints an error message.

  subroutine error_0(mess)
    implicit none
    character(len=*) , intent(in)    :: mess
    write(*,*) trim(adjustl(mess)) 
    return
  end subroutine error_0
  subroutine error_1(mess,i)
    implicit none
    character(len=*) , intent(in) :: mess
    integer          , intent(in) :: i
    integer                       :: j
    write(*,*) trim(adjustl(mess//' --> type any key to continue...')) ; read * , j
    return
  end subroutine error_1

  ! ----------------------------------------------------------------------------
  ! this returns a string with the interger supplied by the user

  function int2text(inter) result(txt)
    implicit none
    integer                        :: inter
    character(len=100)             :: txt0
    character(len=:) , allocatable :: txt
    write(txt0,'(i0)') inter 
    allocate(character(len=len(trim(adjustl(txt0)))) :: txt)
    txt = trim(adjustl(txt0))
    return
  end function int2text

  ! ----------------------------------------------------------------------------
  ! this returns a string with the real number supplied by the user (w/ 2 decimals)

  function real2text(realnum,decs) result(txt)
    implicit none
    real(dp)                       :: realnum
    integer                        :: decs0
    integer , optional             :: decs
    character(len=100)             :: fmt, txt0
    character(len=:) , allocatable :: txt
    decs0 = 2 ; if (present(decs)) decs0 = decs
    write(fmt, '(A,A,I0,A)') '(','F20.', decs0, ')'
    write(txt0, fmt=fmt) realnum
    allocate(character(len=len(trim(adjustl(txt0)))) :: txt)
    txt = trim(adjustl(txt0))
    return
  end function real2text

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  







  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! statistics
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! ----------------------------------------------------------------------------
  ! this function returns the avarage of a variable "var"

  function varmean(var,w,mask) result(meanvar)

    implicit none
    real(dp)               :: var(:),meanvar
    real(dp) , optional    :: w(:)
    logical  , optional    :: mask(:)
    logical  , allocatable :: mask1(:)
    real(dp) , allocatable :: weig(:)

    allocate(weig(size(var)))
    allocate(mask1(size(var)))

    weig(:) = uno
    mask1   = .true.
    meanvar = cero
    
    if (present(mask)) then
      if (size(var).eq.size(mask)) then
        mask1 = mask
      else
        call error('error in varmean: mask of incorrect size',0)
        return
      end if
    end if

    if (present(w)) then
      if (size(var).ne.size(w)) then
        call error('error in varmean!! var and w have different size',0)
        return
      else
        weig(:) = w(:)
      end if
      if (sum(weig).lt.tolvl) then
        call error('error in varmean!! w are cero',0)
        return
      end if
    end if

    meanvar = sum(var(:)*weig(:),mask=mask1)/sum(weig,mask=mask1)

    deallocate(weig)
    deallocate(mask1)

    return
  end function varmean

  ! ----------------------------------------------------------------------------
  ! this function returns the standard deviation of a variable "var"

  function varvar(var,w,mask) result(varian)

    implicit none
    real(dp)               :: var(:),varian,mvar
    real(dp) , optional    :: w(:)
    logical  , optional    :: mask(:)
    logical  , allocatable :: mask1(:)
    real(dp) , allocatable :: weig(:)

    allocate(weig(size(var)))
    allocate(mask1(size(var)))

    weig(:) = uno
    mask1   = .true.
    varian  = cero
    mvar    = cero

    if (present(mask)) then
      if (size(var).eq.size(mask)) then
        mask1 = mask
      else
        call error('error in varvar: mask of incorrect size',0)
        return
      end if
    end if

    if (present(w)) then
      if (size(var).ne.size(w)) then
        call error('error in varvar!! var and w have different size',0)
        return
      else
        weig(:) = w(:)
      end if
      if (sum(weig).lt.tolvl) then 
        call error('error in varvar!! w are cero',0)
        return
      end if
    end if

    mvar   = varmean(var,weig,mask=mask1.and..not.isnan(var))
    varian = sum(weig*((var-mvar)**dble(2.00)),mask=mask1)/sum(weig,mask=mask1)

    deallocate(weig)
    deallocate(mask1)

    return
  end function varvar

  ! ----------------------------------------------------------------------------
  ! this function returns the standard deviation of a variable "var"

  function varstd(var,w,mask) result(stdvar)

    implicit none
    real(dp)               :: var(:),stdvar
    real(dp) , optional    :: w(:)
    logical  , optional    :: mask(:)
    logical  , allocatable :: mask1(:)
    real(dp) , allocatable :: weig(:)

    allocate(weig(size(var)))
    allocate(mask1(size(var)))

    weig(:) = uno
    mask1   = .true.
    stdvar  = cero

    if (present(mask)) then
      if (size(var).eq.size(mask)) then
        mask1 = mask
      else
        call error('error in varstd: mask of incorrect size',0)
        return
      end if
    end if

    if (present(w)) then
      if (size(var).ne.size(w)) then
        call error('error in varstd!! var and w have different size',0)
        return
      else
        weig(:) = w(:)
      end if
      if (sum(weig).lt.tolvl) then
        call error('error in varstd!! w are cero',0)
        return
      end if
    end if

    stdvar = sqrt(varvar(var,weig,mask=mask1))

    deallocate(weig)
    deallocate(mask1)

    return
  end function varstd

  ! ----------------------------------------------------------------------------
  ! this function returns the correlation coefficient between two variables "xvar1" and "xvar2"

  function correlation(xvar1,xvar2,w,mask) result(corr)

    implicit none
    real(dp) , optional    :: w(:)
    real(dp)               :: xvar1(:),xvar2(:),corr
    real(dp)               :: aux1,aux2,aux3,aux4,aux5
    logical  , optional    :: mask(:)
    logical  , allocatable :: mask1(:)
    real(dp) , allocatable :: weig(:)

    allocate(weig(size(xvar1)))
    allocate(mask1(size(xvar1)))

    corr    = cero
    weig(:) = uno
    mask1   = .true.
  
    if (present(mask)) then
      if (size(xvar1).eq.size(mask)) then
        mask1 = mask
      else
        call error('error in varmean: mask of incorrect size',0)
        return
      end if
    end if

    if (size(xvar1).ne.size(xvar2)) then
      call error('error in correaltion!! yvar and xvar of different sizes',0)
      return
    end if
    if (present(w)) then
      if (size(xvar1).ne.size(w)) then
        call error('error in correlation!! var and w have different size',0)
        return
      end if
      if (sum(w).lt.tolvl) then
        call error('error in correlation!! w are cero',0)
        return
      end if
      weig(:) = w(:)
    end if

    aux1 = varmean(xvar1,weig,mask=mask1)
    aux2 = varstd(xvar1,weig,mask=mask1)
    aux3 = varmean(xvar2,weig,mask=mask1)
    aux4 = varstd(xvar2,weig,mask=mask1)

    if (aux2.lt.tolvl) then
      call error('error in correlation!! xvar1 is a constant',0)
      return
    end if
    if (aux4.lt.tolvl) then
      call error('error in correlation!! xvar2 is a constant',0)
      return
    end if

    aux5 = sum(weig(:)*( (xvar1(:) - aux1)*(xvar2(:) - aux3) ),mask=mask1)/sum(weig,mask=mask1)
    corr = aux5/(aux2*aux4)

    deallocate(weig)
    deallocate(mask1)

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
    real(dp)               :: perc

    allocate(weig(size(xvec)))
    allocate(mask1(size(xvec)))

    weig(:) = uno
    cutoff  = minval(xvec)
    mask1   = .true.
  
    if (present(mask)) then
      if (size(xvec).eq.size(mask)) then
        mask1 = mask
      else
        call error('error in percentile: mask of incorrect size',0)
        return
      end if
    end if

    if (present(w)) then
      if (size(w).ne.size(xvec)) then
        call error('error in percentile!! var and w have different size',0)
        return
      end if
      if (sum(w).lt.tolvl) then
        call error('error in percentile!! w are cero',0)
        return
      end if
      weig(:) = w(:)
    end if

    if (pct.lt.cero) then
      call error('error in percentile!! invalid percentile: negative value',0)
      return
    end if
    if (pct.gt.uno) then
      if (pct.lt.cien) then
        perc = pct/cien
      else
        call error('error in percentile!! invalid percetile: larger than 100',0)
        cutoff = maxval(xvec)
        return
      end if
    else
      perc = pct
    end if

    iter = 0
    aux1 = maxval(xvec)
    aux2 = minval(xvec)
    do while ( abs(aux2-aux1).gt.tolvl .and. iter.lt.5000 ) ; iter = iter + 1
      aux3 = medio*(aux1+aux2)
      aux4 = sum(weig,mask = xvec.le.aux3 .and. mask1)/sum(weig,mask=mask1)
      if (aux4.le.perc) aux2 = aux3
      if (aux4.ge.perc) aux1 = aux3
    end do

    cutoff = medio*(aux1+aux2)

    deallocate(weig)
    deallocate(mask1)

    return
  end function percentile

  ! ----------------------------------------------------------------------------
  ! these functions return the coefficients from a ols regresion.
  ! it allows 1 and 2 explanatory variables, and returns the coefficients
  ! either as separate double precision variables or in a vector.
  ! The user can safely call the subroutine "ols(·)" without the need to specify 
  ! the subroutine. The program will automatically call the corresponding subroutine
  ! depending on the number of variables and on the format of output variables

  subroutine olsreg(coeffs,yvec,x1vec,x2vec,x3vec,x4vec,x5vec,x6vec,x7vec,x8vec,w,mask,table,stds)
    implicit none
    real(dp) , intent(out)            :: coeffs(:)
    real(dp) , intent(out) , optional :: stds(:)
    real(dp) , intent(in)             :: yvec(:)
    real(dp) , intent(in) , optional  :: x1vec(:),x2vec(:),x3vec(:),x4vec(:),x5vec(:),x6vec(:),x7vec(:),x8vec(:)
    real(dp) , intent(in) , optional  :: w(:)
    logical  , intent(in) , optional  :: mask(:)
    logical  , intent(in) , optional  :: table
    logical  , allocatable            :: mask1(:)
    real(dp) , allocatable            :: yvar(:),xvars(:,:),zvar(:),wvar(:)
    real(dp) , allocatable            :: xTx(:,:),ixTx(:,:),xTy(:),wx(:,:)
    real(dp) , allocatable            :: evar(:),sdbeta(:,:),sdcoefs(:),tstats(:),inter(:,:)
    real(dp)                          :: rsq,arsq    
    integer                           :: i,n0,no,nc,wc,nx ; nx = 1

    n0     = size(yvec)
    coeffs = -huge(uno)

    ! mask – valid observations
    allocate(mask1(n0)) ; mask1 = .true.
    if (present(mask)) then
      if (size(mask).ne.n0) then
        call error('error in olsreg!! mask of incorrect size',0)   
        return
      else
        mask1 = mask
      end if
    end if

    ! number of valid observations
    no = count(mask1)  
    if (no.eq.0) then 
      call error('error in olsreg!! no (valid) observations',0)
      return
    end if

    ! weigths
    allocate(wvar(n0)) ; wvar = uno
    if (present(w)) then
      if (size(w).ne.n0) then
        call error('error in olsreg!! weigths of incorrect size',0)
        return
      else
        wvar = w
      end if
    end if

    ! check number of explanatory variabkes
    if (present(x1vec)) then ; nx = 1
      if (n0.ne.size(x1vec)) then
        call error('error in olsreg!! yvec and x1vec different observations',0)
        return
      end if
    end if
    if (present(x2vec)) then ; nx = 2
      if (n0.ne.size(x2vec)) then
        call error('error in olsreg!! yvec and x2vec different observations',0)
        return
      end if
    end if
    if (present(x3vec)) then ; nx = 3
      if (n0.ne.size(x3vec)) then
        call error('error in olsreg!! yvec and x3vec different observations',0)
        return
      end if
    end if
    if (present(x4vec)) then ; nx = 4
      if (n0.ne.size(x4vec)) then
        call error('error in olsreg!! yvec and x4vec different observations',0)
        return
      end if
    end if
    if (present(x5vec)) then ; nx = 5
      if (n0.ne.size(x5vec)) then
        call error('error in olsreg!! yvec and x5vec different observations',0)
        return
      end if
    end if
    if (present(x6vec)) then ; nx = 6
      if (n0.ne.size(x6vec)) then
        call error('error in olsreg!! yvec and x6vec different observations',0)
        return
      end if
    end if
    if (present(x7vec)) then ; nx = 7
      if (n0.ne.size(x7vec)) then
        call error('error in olsreg!! yvec and x7vec different observations',0)
        return
      end if
    end if
    if (present(x8vec)) then ; nx = 8
      if (n0.ne.size(x8vec)) then
        call error('error in olsreg!! yvec and x8vec different observations',0)
        return
      end if
    end if
    
    nc = size(coeffs)  ! number of coefficients
    wc = nc-nx         ! with constant = 1, 0 if not

    if (wc.ne.1 .and. wc.ne.0) then 
      print * , nc , nx , wc
      call error('error in olsreg!! '//num2text(nc)//' coefficients and '//num2text(nx)//' variables',0)
      return
    end if

    if (present(stds)) then
      if (size(stds).ne.size(coeffs)) call error('error in olsreg!! stds and coeffs of different size',0)
      return
    end if

    ! allocate y-variable and weigths
    allocate(xvars(no,nc)) ; xvars = cero
    allocate(yvar(no))     ; yvar = pack(yvec,mask1)
    allocate(zvar(no))     ; zvar = pack(wvar,mask1)
    
    ! check weigths are positive    
    if (sum(zvar).lt.tolvl) then
      call error('error in olsreg!! weigths are cero',0)
      return
    end if
    
    ! fill constant (if exists) and first variable
    if (wc.eq.1) xvars(:,1) = uno
    
    ! fill remaining variables
    if (present(x1vec)) xvars(:,1+wc) = pack(x1vec(:),mask1)
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
    xTx    = matmul(transmat(xvars),wx)
    ixTx   = inverse(xTx)
    xTy    = matmul(transmat(xvars),yvar(:)*zvar(:))
    coeffs = matmul(ixTx,xTy)

    allocate(evar(no),sdbeta(nc,nc),sdcoefs(nc),tstats(nc),inter(nc,2))

    evar    = yvar - matmul(xvars,coeffs)
    sdbeta  = ixTx*( sum(zvar*evar*evar) / dble(no-nc) )
    sdcoefs = sqrt(diag(sdbeta))
    tstats  = abs(coeffs/sdcoefs)

    if (present(stds)) then
      stds = sdcoefs
    end if

    if (present(table)) then 
      if (table .eqv. .true.) then

        inter(:,1) = coeffs - dble(1.9601068)*sdcoefs
        inter(:,2) = coeffs + dble(1.9601068)*sdcoefs

        rsq  = uno - varvar(evar,w=zvar)/varvar(yvar,w=zvar)
        arsq = uno - ( (uno-rsq)*dble(no-1)/dble(no-nc) )

        write(*,99) '   '
        write(*,99) '   '
        write(*,'(a,i7  )') '                               Number of variables = ',nc
        write(*,'(a,i7  )') '                             Number of observatios = ',n0
        write(*,'(a,i7  )') '                       Number of valid observatios = ',no
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

      end if
    end if

    deallocate(evar,sdbeta,sdcoefs,tstats,inter)
    deallocate(xvars,xTx,ixTx,xTy,yvar,zvar,mask1,wvar)

    return
    99 format (a,f10.4,f10.4,f10.4,f10.4,f10.4)
    98 format (a,i1,a,f10.4,f10.4,f10.4,f10.4,f10.4)
  end subroutine olsreg

  ! ----------------------------------------------------------------------------
  ! this function returns the transition matrix for a discretized ar(1) process
  ! with gaussian shocks

  subroutine tauchen(xvec,rho,mu,sigma,n,pmat)
    implicit none
    integer  , intent(in)  :: n
    real(dp) , intent(in)  :: rho,mu,sigma,xvec(n)
    real(dp) , intent(out) :: pmat(n,n)
    integer                :: i
    do i=1,n
      pmat(i,:) = cero
      call normaldist(xvec,mu+rho*xvec(i),sigma,n,pmat(i,:))
    end do
    return
  end subroutine tauchen

  ! ----------------------------------------------------------------------------
  ! this function returns the transition matrix for a discretized ar(1) process
  ! with gaussian shocks

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

  ! ----------------------------------------------------------------------------
  ! this subroutine takes an array "nums" of dimension "n" and fills it with
  ! integer numbers from 1 to n and them shuffled using the the Fisher–Yates
  ! shuffle algorithm

  subroutine fisheryates(nums,seed)
    implicit none
    integer , intent(inout) :: nums(:)
    integer , intent(in) , optional :: seed(:)
    real(dp) :: r
    integer :: i, j, temp , n ; n = size(nums)

    ! Initialize the array with values 1 to 100
    do i = 1,n ; nums(i) = i ; end do

    ! Seed the random number generator
    if (present(seed)) then
      call random_seed(put=seed)
    else
      call random_seed()
    end if

    ! Fisher–Yates shuffle
    do i = n, 2, -1
      call random_number(r)
      j = int(r * real(i)) + 1
      temp    = nums(i)
      nums(i) = nums(j)
      nums(j) = temp
    end do

    return
  end subroutine fisheryates

  ! ----------------------------------------------------------------------------
  ! this subroutine takes an input vector "input_vec" of size "n" and returns a
  ! vector "output_vec" with "m" shuffled values from "input_vec". The input vector
  ! is shuffled using the Fisher–Yates algorithm.
  ! if n > m: "output_vec" filled with "m" shuffled values of "input_vec"
  ! if n < m: "output_vec" filled with "n" shuffled values of "input_vec" until m. 

  subroutine shuffle_vect_int(input_vec, output_vec)
    implicit none
    integer , intent(in)  :: input_vec(:)
    integer , intent(out) :: output_vec(:)
    integer , allocatable :: indices(:)
    integer               :: i,n,m
    n = size(input_vec)
    m = size(output_vec)
    allocate(indices(n))
    call fisheryates(indices)
    do i = 1,m
      output_vec(i) = input_vec(indices(mod(i-1,n)+1))
    end do
    deallocate(indices)
    return
  end subroutine shuffle_vect_int
  subroutine shuffle_vect_dp(input_vec, output_vec)
    implicit none
    real(dp) , intent(in)  :: input_vec(:)
    real(dp) , intent(out) :: output_vec(:)
    integer  , allocatable :: indices(:)
    integer                :: i,n,m
    n = size(input_vec)
    m = size(output_vec)
    allocate(indices(n))
    call fisheryates(indices)
    do i = 1,size(output_vec)
      output_vec(i) = input_vec(indices(mod(i-1,n)+1))
    end do
    deallocate(indices)
    return
  end subroutine shuffle_vect_dp

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
  function vectorize_int_2d(mat) result(vec)
    implicit none
    integer :: mat(:,:)
    integer :: vec(size(mat,1)*size(mat,2))
    vec = reshape(mat,(/size(vec)/))
    return
  end function vectorize_int_2d
  function vectorize_int_3d(mat) result(vec)
    implicit none
    integer :: mat(:,:,:)
    integer :: vec(size(mat,1)*size(mat,2)*size(mat,3)) 
    vec = reshape(mat,(/size(vec)/))
    return
  end function vectorize_int_3d
  function vectorize_int_4d(mat) result(vec)
    implicit none
    integer :: mat(:,:,:,:)
    integer :: vec(size(mat,1)*size(mat,2)*size(mat,3)*size(mat,4))
    vec = reshape(mat,(/size(vec)/))
    return
  end function vectorize_int_4d
  function vectorize_int_5d(mat) result(vec)
    implicit none
    integer :: mat(:,:,:,:,:)
    integer :: vec(size(mat,1)*size(mat,2)*size(mat,3)*size(mat,4)*size(mat,5))
    vec = reshape(mat,(/size(vec)/))
    return
  end function vectorize_int_5d
  function vectorize_lo_2d(mat) result(vec)
    implicit none
    logical :: mat(:,:)
    logical :: vec(size(mat,1)*size(mat,2))
    vec = reshape(mat,(/size(vec)/))
    return
  end function vectorize_lo_2d
  function vectorize_lo_3d(mat) result(vec)
    implicit none
    logical :: mat(:,:,:)
    logical :: vec(size(mat,1)*size(mat,2)*size(mat,3)) 
    vec = reshape(mat,(/size(vec)/))
    return
  end function vectorize_lo_3d
  function vectorize_lo_4d(mat) result(vec)
    implicit none
    logical :: mat(:,:,:,:)
    logical :: vec(size(mat,1)*size(mat,2)*size(mat,3)*size(mat,4))
    vec = reshape(mat,(/size(vec)/))
    return
  end function vectorize_lo_4d
  function vectorize_lo_5d(mat) result(vec)
    implicit none
    logical :: mat(:,:,:,:,:)
    logical :: vec(size(mat,1)*size(mat,2)*size(mat,3)*size(mat,4)*size(mat,5))
    vec = reshape(mat,(/size(vec)/))
    return
  end function vectorize_lo_5d

  ! ----------------------------------------------------------------------------
  ! this function returns a vector "vec1" with the cummmulative sum of the elements
  ! of a vector "vec0"

  function cumsum(vec0) result(vec1)
    implicit none
    real(dp) :: vec0(:),vec1(size(vec0))
    integer  :: i ; i = 1
    vec1 = vec0(i) ; do i=2,size(vec0) ; vec1(i) = sum(vec0(1:i)) ; end do
    return
  end function cumsum

  ! ----------------------------------------------------------------------------
  ! this function returns the main diagonal of a matric "mat"

  function diag(mat) result(vec)
    implicit none
    real(dp) :: mat(:,:),vec(size(mat,1))
    integer  :: i
    if (size(mat,1).ne.size(mat,2)) then
      call error(' error in diag: matrix not square',0)
      return
    end if
    forall (i=1:size(mat,1)) vec(i) = mat(i,i) 
    return
  end function diag

  ! ----------------------------------------------------------------------------
  ! returns the transpose of a matrix "mat"

  function transmat(mat) result(matt)
    implicit none
    real(dp) :: mat(:,:),matt(size(mat,2),size(mat,1))
    integer  :: i
    forall (i=1:size(mat,1)) matt(:,i) = mat(i,:)
    return
  end function transmat

  ! ----------------------------------------------------------------------------
  ! compute the inverse of a sqaure matrix "m"

  function inverse(m) result(im)
    implicit none
    integer  :: n,i,j,k
    real(dp) :: m(:,:),coeff
    real(dp) :: b(size(m,1))
    real(dp) :: d(size(m,1))
    real(dp) :: x(size(m,1))
    real(dp) :: im(size(m,1),size(m,1))
    real(dp) :: l(size(m,1),size(m,1))
    real(dp) :: u(size(m,1),size(m,1))
    real(dp) :: mb(size(m,1),size(m,1))

    if (size(m,1).ne.size(m,2)) then
      call error(' error in inverse: matrix is not square',0)
      return
    end if
    n = size(m,1) ; l = cero ; u = cero ; b = cero ; mb = m
    do k=1,n-1 ; do i=k+1,n
      coeff  = m(i,k)/m(k,k)
      l(i,k) = coeff
      do j=k+1,n
        m(i,j) = m(i,j)-coeff*m(k,j)
      end do
    end do ; end do
    do i=1,n
      l(i,i) = 1.0
    end do
    do j=1,n ; do i=1,j
      u(i,j) = m(i,j)
    end do ; end do
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
  ! its corresponding .
  !
  ! if your problem has states, please check the corresponding subroutine to make
  ! sure your function has the correct form.
  !
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



  ! ----------------------------------------------------------------------------
  ! golden search algorithm
  
  subroutine golden(func,x,y,xmax,xmin,itermax,tol,iprint)

    ! this subroutine finds the maximum of a single-valued function with one unknown
    ! using the golden search algorithm
    !
    ! the inputs are:
    !
    !   - func: user-provided function to be minimize. this function should of the form:
    !
    !       function func(x) result(y)
    !         real(dp) :: x,y
    !       end function func
    !
    !   - xmin: lower bound of variable "x"
    !   - xmax: upper bound of variable "x"
    !
    ! there are also few optional inputs:
    !
    !   itermax:  max number of iterations (default = 500)
    !   tol:      tolerance level (default = 1.0d-8)
    !   iprint:   control the printing behaviour
    !       - iprint = 0: don't print anything (default)
    !       - iprint = 1: print warnings
    !       - iprint = 2: print warnings and every iteration
    !
    ! the outputs of the subroutine are:
    !
    !   x: the value of x that maximizes "func"
    !   y: the value of the function at the point x

    implicit none
    real(dp) , intent(in)            :: xmax,xmin
    real(dp) , intent(out)           :: x,y
    real(dp) , intent(in) , optional :: tol
    integer  , intent(in) , optional :: itermax
    integer  , intent(in) , optional :: iprint
    real(dp) , parameter             :: alpha=0.61803399
    real(dp)                         :: x0,x1,x2,x3,f1,f2,tolgold
    integer                          :: it,maxiter,ip

    interface
      function func(xx) result(ff)
      real(kind=8) :: xx,ff
      end function func
    end interface

    tolgold = 1.0d-8 ; if (present(tol)    ) tolgold = tol
    maxiter = 500    ; if (present(itermax)) maxiter = itermax
    ip      = 0      ; if (present(iprint) ) ip      = iprint

    if (ip.gt.0) write(*,*) ' starting golden search algorithm '
    if (ip.gt.0) write(*,*) '                                  '

    x0 = xmin
    x3 = xmax
    x1 = alpha*x0 + (uno-alpha)*x3
    x2 = alpha*x3 + (uno-alpha)*x1 
    f1 = func(x1)
    f2 = func(x2)

    it = 0
    do while (abs(x0-x3).gt.tolgold*(abs(x2)+abs(x1)) .and. it.lt.maxiter)
      it = it + 1
      if (f2.gt.f1) then
        x0 = x1
        x1 = x2
        f1 = f2
        x2 = alpha*x1 + (uno-alpha)*x3
        f2 = func(x2)
      else
        x3 = x2
        x2 = x1
        f2 = f1
        x1 = alpha*x2 + (uno-alpha)*x0
        f1 = func(x1)
      end if
      if (ip.gt.1) write(*,15) it , x1 , x2 , f1 , f2
    end do

    if (f1.gt.f2) then
      y = f1
      x = x1
    else
      y = f2
      x = x2
    end if

    if (ip.gt.0) write(*,*  ) '   '
    if (ip.gt.0) write(*,16) it , x , y

    return
    15 format( ' iteration = ',i4,' | x1 = ',f10.4,' | x2 = ',f10.4,' | f(x1) = ',f10.4,' | f(x2) = ',f10.4)
    16 format( ' solved: iterations = ',i4,' | x = ',f10.4,' | f(x) = ',f10.4)
  end subroutine golden

  ! ----------------------------------------------------------------------------
  ! brent's method

  subroutine brent(func,x,iy,ind,x0,x1,itermax,tol,iprint)

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
    !
    ! the outputs of the subroutine are:
    !
    !   x:    the value of x that maximizes "func"
    !   iy:   number of function evaluations required
    !   ind:  exit indicator, taking value:
    !           0: solved, func(x) = 0
    !           1: not solved

    implicit none
    real(dp) , intent(out)           :: x
    integer  , intent(out)           :: iy
    integer  , intent(out)           :: ind
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


    ! - ```exitcode``` = 0: the algorithm found a root
    ! - ```exitcode``` = 1: either ```x0``` or ```x1``` is a root of ```func```
    ! - ```exitcode``` = 2: the root is not within the interval (```x0```, ```x1```)
    ! - ```exitcode``` = 3: the points ```x0``` and ```x1``` are too close
    ! - ```exitcode``` = 9: maximum number of function evaluations reached

    toler   = 1.0d-8 ; if (present(tol)    ) toler   = tol
    maxiter = 500    ; if (present(itermax)) maxiter = itermax
    ipri    = 0      ; if (present(iprint) ) ipri    = iprint

    if (ipri.ge.0) then
      write(*,*) '  '
      write(*,*) ' starting brent algorithm'
      write(*,*) '  '
    end if

    if (abs(x1-x0).lt.tol) then
      if (ipri.ge.0) write(*,*) ' x0 and x1 are too close'
      ind = 3 ; return
    end if

    xa = x0 ; ya = func(xa)

    if (abs(ya).lt.toler) then
      if (ipri.ge.0) write(*,*) ' x0 is a root'
      ind = 1 ; x = xa ; return
    end if

    xb = x1 ; yb = func(xb)

    if (abs(yb).lt.toler) then
      if (ipri.ge.0) write(*,*) ' x1 is a root'
      ind = 1 ; x = xb ; return
    end if

    if (ya*yb.gt.cero) then
      if (ipri.ge.0) write(*,*) ' f(x0) and f(x1) have the same sign'
      ind = 2 ; x = xa ; return
    end if

    if (abs(ya).lt.abs(yb)) then
      xc = xb ; xb = xa ; xa = x0 ; yc = ya ; ya = yb
    end if

    xc = xa ; iy = 0 ; ind = 0

    do while (abs(xb-xa).gt.toler .and. iy.lt.maxiter)
      yc = func(xc)
      if (abs(yc-ya).gt.toler .and. abs(yc-yb).gt.toler) then
        xs = xa*yb*yc/((ya-yb)*(ya-yc)) + xb*ya*yc/((yb-ya)*(yb-yc)) + xc*ya*yb/((yc-ya)*(yc-yb))
      else
        xs = xb - yb*(xb-xa)/(yb-ya)
      end if
      ys = func(xs) ; iy = iy + 1
      if (abs(ys).lt.toler) then
        ind = 0 ; x = xs ; return
      end if
      if (iy.eq.maxiter) then
        ind = 9 ; x = xs ; return
      end if
      if (ya*ys.lt.cero) then
        xc = xb ; xb = xs ; yb = ys
      else
        xa = xs ; ya = ys
      end if
      if (abs(ya).lt.abs(yb)) then
        xc = xb ; xb = xa ; xa = xc ; yc = yb ; yb = ya ; ya = yc
      end if
      if (ipri.ge.1) then
        write(*,87) iy , xa , xb , xc , ya , yb , yc
      end if
    end do
    if (ipri.gt.0) then
      write(*,86) iy,xb,yb
    end if

    x = xb

    return
    87 format( ' iteration = ',i4,' | xa = ',f10.4,' | xb = ',f10.4,' | xc = ',f10.4,&
    ' | f(xa) = ',f10.4,' | f(xb) = ',f10.4,' | f(xc) = ',f10.4)
    86 format( ' solved: iterations = ',i4,' | x = ',f10.4,' | f(x) = ',f10.4)
  end subroutine brent

  ! ----------------------------------------------------------------------------
  ! simplex algorithm

  subroutine simplex(func,x,y,iy,ind,x0,itermax,tol,iprint)

    ! this subroutine minimizes a single-valued function with "n" unknowns using the
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

    integer  , intent(out)           :: iy,ind
    real(dp) , intent(out)           :: x(:),y

    real(dp) , intent(in)            :: x0(:)
    real(dp) , intent(in) , optional :: tol
    integer  , intent(in) , optional :: itermax
    integer  , intent(in) , optional :: iprint

    integer              :: n,i,j,ilow,ihigh,ihigh2,maxiter,ipri,ii(1)
    real(dp)             :: xp(size(x,1),size(x,1)+1),yp(size(x,1)+1),xr(size(x,1))
    real(dp)             :: xm(size(x,1)),xe(size(x,1)),xc(size(x,1))
    real(dp)             :: y0,yr,ye,yc,cent,toler
    real(dp) , parameter :: alpha=1.0,beta=0.5,gamma=2.0

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

    if (ipri.gt.1) then
      write(*,*) '  '
      write(*,*) ' evaluating initial simplex'
      write(*,*) '  '
    end if
    
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

    if (ipri.gt.1) then
      write(*,*) '  '
      write(*,*) ' beginning iterations'
      write(*,*) '  '
    end if

    1 continue

    ! sort points
    ii = minloc(yp) ; ilow = ii(1)
    ii = maxloc(yp) ; ihigh = ii(1)
    ihigh2 = ilow
    do i=1,n+1
      if (i.ne.ihigh .and. yp(i).ge.yp(ihigh2)) ihigh2 = i
    end do

    ! ilow = 1
    ! do i=1,n+1
    !   if (yp(i).lt.yp(ilow)) ilow = i
    ! end do
    ! ihigh = 1
    ! do i=1,n+1
    !   if (yp(i).gt.yp(ihigh)) ihigh = i
    ! end do
    ! ihigh2 = ilow
    ! do i=1,n+1
    !   if (yp(i).gt.yp(ihigh2) .and. yp(i).lt.yp(ihigh)) ihigh2 = ihigh
    ! end do

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
      xr(j) = (uno+alpha)*xm(j) - alpha*xp(j,ihigh)
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
        xe(j) = (uno-gamma)*xm(j) - gamma*xp(j,ihigh)
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
        xc(j) = (uno-beta)*xm(j) - beta*xp(j,ihigh)
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
              xp(j,i) = medio*( xp(j,i) + xp(j,ilow) )
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

  subroutine lmmin(func,x,y,iy,ind,x0,itermax,damp,tol,toleach,shock,usebro,iprint)

    ! this subroutine minimizes the sum of squred errors of a system of m equations
    ! in "n" unknowns using the levenberg–marquardt algorithm. more about the algorithm
    ! can be found in https://en.wikipedia.org/wiki/levenberg–marquardt_algorithm.
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
    !   - x0: a vector of dimension "n" with the initial guess
    !
    ! there are also some other optional inputs:
    !
    !   itermax:  max number of iterations (default = 500)
    !
    !   damp:     damping factor (default = 1)
    !
    !   tol:      tolerance level (default = 1.0d-8)
    !
    !   shock:    shock to compute numerical jacobian (default = 0.05)
    !
    !   usebro:   whether to use broyden's method to update jacobian (default = 0)
    !               if broyden's method is used, the program will compute the
    !               numerical jacobian once, and then update it at every iteration
    !
    !   iprint:   indicator for printing results (default = 0)
    !               0: print no results
    !               1: print main results
    !               2: print all the iterations
    !
    ! the outputs of the subroutine are:
    !
    !   x:    vector with the values of x that solves the system
    !
    !   y:    the value of the function at the point x
    !
    !   iy:   the number of total function evaluations
    !
    !   ind:  exit indicator, taking value:
    !           0: system has converged, func(x) = 0
    !           1: jacobian is closed to 0
    !           2: step in x is close to 0
    !           9: max number of iterations reached

    implicit none
    
    ! outputs
    real(dp) , intent(out)           :: x(:)
    real(dp) , intent(out)           :: y(:)
    integer  , intent(out)           :: iy
    integer  , intent(out)           :: ind

    ! inputs (including optionals)
    real(dp) , intent(in)            :: x0(:)
    real(dp) , intent(in) , optional :: shock
    real(dp)                         :: shck(size(x,1))
    real(dp) , intent(in) , optional :: damp
    real(dp)                         :: da
    real(dp) , intent(in) , optional :: tol
    real(dp)                         :: toler
    real(dp) , intent(in) , optional :: toleach
    real(dp)                         :: tolereach
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
    real(dp) :: ja(size(x,1),size(x,1))
    real(dp) :: jt(size(x,1))
    real(dp) :: e0,y0(size(y,1))
    real(dp) :: eb,yb(size(y,1)),xb(size(x,1))
    real(dp) :: e1,y1(size(y,1)),x1(size(x,1))
    real(dp) :: ea,ya(size(y,1)),xa(size(x,1))
    real(dp) :: ej,yj(size(y,1),size(x,1)),xj(size(x,1),size(x,1))
    real(dp) :: ddx

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

     ! compute new numerical jacobian      
    if (iy.eq.1 .or. bro.eq.0 .or. br.eq.0) then
      
      if (ip.gt.1) write (*,*) '  '
      do i = 1,n

        ! shock the i-th parameter
        xj(:,i) = xb(:)
        if (xb(i).lt.-tol) then
          xj(i,i) = xb(i) - min(shck(i),shck(i)*max(0.30,eb))*abs(xb(i))
        elseif (xb(i).gt.tol) then
          xj(i,i) = xb(i) + min(shck(i),shck(i)*max(0.30,eb))*abs(xb(i))
        elseif (abs(xb(i)).le.tol) then
          xj(i,i) = min(shck(i),shck(i)*max(0.30,eb))
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
          if (abs(xb(i)).gt.tol) shck(i) = abs( xb(i) - xj(i,i) ) / abs( xb(i) )
          if (abs(xb(i)).le.tol) shck(i) = abs( xj(i,i) )
        end if

        ! if change in func is too small, compute new shocked parameter
        if (abs(ej-eb).lt.tolvl) then ; ddx = xj(i,i) - xb(i)
          do while (abs(ej-eb).lt.tolvl)
            xj(i,i) = xj(i,i) + ddx
            yj(:,i) = func(xj(:,i)) ; ej = sum(yj(:,i)*yj(:,i)) ; iy = iy + 1
            if (ip.gt.1) write (*,90) i,ej
            if (iy.ge.maxiter) goto 19
          end do
          if (abs(xb(i)).gt.tol) shck(i) = abs( xb(i) - xj(i,i) ) / abs( xb(i) )
          if (abs(xb(i)).le.tol) shck(i) = abs( xj(i,i) )
        end if

        ! compute numerical jacobian
        do k = 1,m
          j(k,i) = ( yj(k,i) - yb(k) )/( xj(i,i) - xb(i) )
        end do

      end do
      qp = 0

      if (ip.gt.1) write (*,*) '  '
    
    ! if broyden's method is used, update jacobian with broyden
    else

      j0 = j ; call broyden(j,j0,x1,xa,y1,ya)

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
      ja(i,i) = (uno+da)*jj(i,i)
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
    real(dp) , intent(in)  :: xmax,xmin,x
    real(dp) , intent(out) :: y
    if (xmax.lt.xmin) then
      call error(' errror in normalize: xmax > xmin',0)
      return
    end if
    if (x.gt.xmax) then
      call error(' errror in normalize: x > xmax',0)
      return
    end if
    if (x.lt.xmin) then 
      call error(' errror in normalize: x < xmin',0)
      return
    end if
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
    if (xmax.lt.xmin) then
      call error(' errror in denormalize: xmax > xmin',0)
      return
    end if
    x = (exp(y)/(uno+exp(y)))*(xmax-xmin) + xmin
    if (isnan(x) .and. y.gt.cien) x = xmax
    return
  end subroutine denormalize

  ! ----------------------------------------------------------------------------
  ! this subroutines updates a jacobian matrix using the broyden's method.
  ! it is used by the levenberg–marquardt algorithm.

  subroutine broyden(j1,j0,x1,x0,f1,f0)
    implicit none
    real(dp) , intent(in)  :: x1(:),f1(:),x0(:),f0(:),j0(:,:)
    real(dp) , intent(out) :: j1(:,:)
    real(dp)               :: df(size(f0,1))
    real(dp)               :: dx(size(x1,1))
    real(dp)               :: den
    integer                :: i
    j1 = j0
    if (size(f1).ne.size(j0,1)) then
      call error(' error in broyden: size(y) /= size(j,1)',0)
      return
    end if
    if (size(x1).ne.size(j0,2)) then
      call error(' error in broyden: size(x) /= size(j,2)',0)
      return
    end if
    den = sum(dx(:)*dx(:))
    if (abs(den).lt.tolvl) then
      call error(' error in broyden: denominator close to 0',0)
      return
    end if
    dx(:) = x1(:)-x0(:)
    df(:) = f1(:)-f0(:)
    do i=1,size(df)
      j1(i,:) = j0(i,:) + dx(:)*( df(i)-sum(j0(i,:)*dx(:)) ) / den
    end do
    return
  end subroutine broyden

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module toolkit
