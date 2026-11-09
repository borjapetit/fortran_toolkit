
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! toolkit.f90, a toolkit for fortran90 programming
! Borja Petit, CUNEF Universidad, 2025
!
! The code is distributed under the MIT license.
! You are free to use and modify it as you wish acknowledging the source.
! If you modified version of this code, indicate the changes made along with your name and date.
!
! UPDATES:
! - October 2025 (Borja Petit):
!     * the functions "ces" and "crra" were added
!     * some minor typos and inconsistencies were fixed. 
!     * modified the "error" call in different subroutines/functions to stop the execution when an error is found.
! - April 2026 (Borja Petit)
!     * add a global variable "print_err" to deactivate error printing (by setting it to .false.)
!     * some minor errors
! - November 2026 (Borja Petit)
!     * simplex: corrected the expansion and contraction points; adaptive coefficients (Gao and Han, 2012);
!       convergence on the spread of func; optional inputs tolx, step and fzero.
!     * lmmin: bounded enlargement of the jacobian step; cholesky solve with scaling floor; gain-ratio damping
!       (Nielsen, 1999); stop after maxstall unproductive iterations (ind = 3); reuse of jacobian points;
!       optional inputs tolrel and maxstall.
!     * inverse: gauss-jordan with partial pivoting; the input matrix is not modified; error if singular.
!     * broyden: dimensions checked first; no update (instead of nan) when the step is zero.
!     * normalize/denormalize: symmetric bounds in normalize (no -infinity at xmin); overflow-safe
!       denormalize; corrected error messages.
!     * golden: initial interior point x2 placed at the golden ratio.
!     * percentile: pct = 100 accepted; convention for pct documented.
!     * randomnormal: the vector of draws is no longer re-centered to have mean exactly mu.
!     * timing: seconds since 1970-01-01 utc (with milliseconds); correct across months and years.
!     * grid, interpolation, interpolate: return after errors; nan and empty-grid handling;
!       unused argument "outof" removed from interpolation.
!     * error: the pause waits for enter (any input) instead of an integer.
!     * num2text: zero-padding of negative integers (-005).
!     * varvar: nan observations excluded from both the mean and the variance.
!     * correlation: corrected error messages.
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! MODULES INCLUDED IN THIS CODE:
! ******************************
!
!   - toolkit_pars:   parameters and constants
!   - toolkit_tools:  general purpose functions and subroutines
!   - toolkit_econ:   econ-related functions
!   - toolkit_linalg: linear algebra
!   - toolkit_stats:  statistics
!   - toolkit_optim:  optimization and root-finding
!   - toolkit:        loads all the modules above ("use toolkit" gives access to everything)
!
! the list of functions and subroutines of each module is at the top of the module.
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module toolkit_pars

  ! ----------------------------------------------------------------------------
  ! parameters and constants
  ! ^^^^^^^^^^^^^^^^^^^^^^^^
  !
  !   - dp: double precision kind
  !   - cero, medio, uno, dos, cien, mil: the numbers 0, 0.5, 1, 2, 100 and 1000 in double precision
  !   - tolvl: default tolerance level (1.0d-10)
  !   - print_err: if .false., error messages are not printed (default = .true.)
  ! ----------------------------------------------------------------------------
    
  implicit none

  integer  , parameter :: dp    = kind(1.0d00)
  real(dp) , parameter :: cero  = dble(0.00000000000000)
  real(dp) , parameter :: medio = dble(0.50000000000000)
  real(dp) , parameter :: uno   = dble(1.00000000000000)
  real(dp) , parameter :: dos   = dble(2.00000000000000)
  real(dp) , parameter :: cien  = dble(100.000000000000)
  real(dp) , parameter :: mil   = dble(1000.00000000000)
  real(dp) , parameter :: tolvl = dble(0.00000000010000)
  logical              :: print_err = .true.

end module toolkit_pars

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module toolkit_tools

  ! ----------------------------------------------------------------------------
  ! general propuse functions and subroutines
  ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  !
  !   - grid: generates a grid for a continuous variable
  !   - interpolation: interpolates a value over a grid, returning position and distance
  !   - interpolate: linearly interpolates a value over an n-dimensional grid, with n <= 6
  !   - vect: transforms a matrix of dimension (n x m) into a vector of n·m rows
  !   - timing: returns the number of seconds since 00:00h of the 1st day of the month [robust to parallelization]
  !   - multiplo: returns 1 if an integer is a multiple of another user-provided integer
  !   - iseven: returns 1 if a user-provided integer is even
  !   - error: prints an error message
  !   - num2text: returns a string with an integer or a real number
  !   - fisheryates: shuffles an integer vector using the Fisher–Yates algorithm
  !   - shuffle_vect: shuffles a vector of integers or reals using the Fisher–Yates algorithm
  ! ----------------------------------------------------------------------------

  use toolkit_pars  , only : dp,cero,medio,uno,cien,tolvl,print_err
  implicit none

  interface interpolate
    module procedure interpolate1d,interpolate2d,interpolate3d,&
                     interpolate4d,interpolate5d,interpolate6d
  end interface interpolate

  interface vect
    module procedure vectorize_int_2d,vectorize_dp_2d,vectorize_lo_2d,&
                     vectorize_int_3d,vectorize_dp_3d,vectorize_lo_3d,&
                     vectorize_int_4d,vectorize_dp_4d,vectorize_lo_4d,&
                     vectorize_int_5d,vectorize_dp_5d,vectorize_lo_5d
  end interface vect

  interface error
    module procedure error_0,error_1,error_2
  end interface error

  interface num2text
    module procedure int2text,real2text
  end interface num2text

  interface shuffle_vect
    module procedure shuffle_vect_int,shuffle_vect_dp
  end interface shuffle_vect

  contains

  ! ----------------------------------------------------------------------------
  ! this function creates a grid between "minv" and "maxv" with a curvature of "s"
  !   - if s>1: more grid points around "minv"
  !   - if s<1: more grid points around "maxv"
  function grid(maxv,minv,n,s) result(v)
    implicit none
    integer             :: i,n
    real(dp) , optional :: s
    real(dp)            :: maxv,minv,v(n),ss,xmin,xmax
    v(:) = cero ; ss = uno ; if (present(s)) ss = s
    if (ss.le.cero) then
      call error(' error in grid: curvature parameter s is nonpositive')
      return
    end if
    if (n.lt.2) then
      call error(' error in grid: the grid needs at least 2 points (n < 2)')
      return
    end if
    xmin = min(maxv,minv)
    xmax = max(maxv,minv)
    forall (i=1:n) v(i) = ((dble(i-1)/dble(n-1))**ss)*(xmax-xmin) + xmin
    return
  end function grid

  ! ----------------------------------------------------------------------------
  ! this subroutine linearly interpolates a variable !xnow" over a grid "xgrid".
  ! it finds the closest point in the grid and return its position and the relative distance
  subroutine interpolation(pos,wth,xnow,xgrid)
    implicit none
    integer               :: j,n
    real(dp), intent(in)  :: xnow,xgrid(:)
    real(dp), intent(out) :: wth
    integer , intent(out) :: pos

    pos = 0 ; wth = cero
    n   = size(xgrid)

    if (n.lt.1) then
      call error(' error in interpolation: xgrid is empty')
      return
    end if

    ! if xnow is nan, return a valid position and wth = nan, so that the
    ! interpolated value is nan (instead of accessing out-of-bounds elements)
    if (isnan(xnow)) then
      call error(' error in interpolation: xnow is nan')
      pos = min(2,n) ; wth = xnow
      return
    end if

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

    if (pos.lt.1   ) call error(' error in interpolation: pos < 1')
    if (pos.gt.n   ) call error(' error in interpolation: pos > size')
    if (isnan(wth) ) call error(' error in interpolation: wth is nan')
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
    real(dp) :: y1(:),m(:),x1,xi,wth ; xi = cero
    if (size(y1).ne.size(m)) then
      call error(' error in interpolate: 1st dimension incorrect')
      return
    end if
    if (size(y1).lt.2) then
      call error(' error in interpolate: 1st dimension with less than 2 elements')
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
    real(dp) :: y1(:),y2(:),m(:,:) ; xi = cero
    if (size(m,2).ne.size(y2)) then
      call error(' error in interpolate: 2nd dimension incorrect')
      return
    end if
    if (size(y2).lt.2) then
      call error(' error in interpolate: 2nd dimension with less than 2 elements')
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
    real(dp) :: wth3,m(:,:,:) ; xi = cero
    if (size(m,3).ne.size(y3)) then
      call error(' error in interpolate: 3rd dimension incorrect')
      return
    end if
    if (size(y3).lt.2) then
      call error(' error in interpolate: 3rd dimension with less than 2 elements')
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
    real(dp) :: wth4,m(:,:,:,:) ; xi = cero
    if (size(m,4).ne.size(y4)) then
      call error(' error in interpolate: 4th dimension incorrect')
      return
    end if
    if (size(y4).lt.2) then
      call error(' error in interpolate: 4th dimension with less than 2 elements')
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
    real(dp) :: wth5,m(:,:,:,:,:) ; xi = cero
    if (size(m,5).ne.size(y5)) then
      call error(' error in interpolate: 5th dimension incorrect')
      return
    end if
    if (size(y5).lt.2) then
      call error(' error in interpolate: 5th dimension with less than 2 elements')
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
    real(dp) :: wth6,m(:,:,:,:,:,:) ; xi = cero
    if (size(m,6).ne.size(y6)) then
      call error(' error in interpolate: 6th dimension incorrect')
      return
    end if
    if (size(y6).lt.2) then
      call error(' error in interpolate: 6th dimension with less than 2 elements')
      return
    end if
    call interpolation(pos6,wth6,x6,y6)
    xi = interpolate5d(x1,x2,x3,x4,x5,y1,y2,y3,y4,y5,m(:,:,:,:,:,pos6))*wth6 + &
         interpolate5d(x1,x2,x3,x4,x5,y1,y2,y3,y4,y5,m(:,:,:,:,:,pos6-1))*(uno-wth6)
    return
  end function interpolate6d

  ! ----------------------------------------------------------------------------
  ! this functions returns a timing number that is robust to parallel computing
  ! (wall-clock time). it returns the number of seconds (with milliseconds) since
  ! 00:00h UTC of january 1st, 1970, so that differences are correct across days,
  ! months, years and changes in the local time offset.
  ! the variable "mode" controls how time is measured:
  !   - if mode=1, time is measured in seconds (default).
  !   - if mode=2, time is measures in minutes.
  !   - if mode=3, time is measured in hours
  function timing(mode) result(time)
    implicit none
    integer , optional :: mode
    integer            :: v1(8),yy,era,yoe,doy,doe
    integer(8)         :: days
    real(dp)           :: time
    call date_and_time(values=v1)
    ! days since 1970-01-01 in the civil (gregorian) calendar
    yy = v1(1) ; if (v1(2).le.2) yy = yy - 1
    era  = yy/400
    yoe  = yy - era*400
    doy  = (153*(v1(2) + merge(-3,9,v1(2).gt.2)) + 2)/5 + v1(3) - 1
    doe  = yoe*365 + yoe/4 - yoe/100 + doy
    days = int(era,8)*146097_8 + int(doe,8) - 719468_8
    ! seconds since 1970-01-01 00:00h utc (v1(4) is the local offset in minutes)
    time = dble(days)*86400.0_dp + dble(v1(5)*3600 + v1(6)*60 + v1(7)) &
         + dble(v1(8))/1000.0_dp - dble(v1(4))*60.0_dp
    if (present(mode)) then
      if (mode.eq.3) time = time/dble(60*60)  ! measured in hours
      if (mode.eq.2) time = time/dble(60)     ! measured in minutes
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
      call error(' error in multiplo: zero has no multiples')
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
    if (print_err.eqv..false.) return
    write(*,*) trim(adjustl(mess)) 
    return
  end subroutine error_0
  subroutine error_1(mess,i)
    implicit none
    character(len=*) , intent(in) :: mess
    integer          , intent(in) :: i
    character(len=1)              :: ans
    if (print_err.eqv..false.) return
    if (i.eq.0) then
      call error_0(mess)
    else
      write(*,*) trim(adjustl(mess//' --> press enter to continue...')) ; read(*,'(a)') ans
    end if
    return
  end subroutine error_1
  subroutine error_2(mess,stpe)
    implicit none
    character(len=*) , intent(in) :: mess
    logical          , intent(in) :: stpe
    character(len=1)              :: ans
    if (print_err.eqv..false.) return
    if (stpe.eqv..false.) then
      call error_0(mess)
    else
      write(*,*) trim(adjustl(mess//' --> press enter to continue...')) ; read(*,'(a)') ans
    end if
    return
  end subroutine error_2

  ! ----------------------------------------------------------------------------
  ! this returns a string with the interger supplied by the user. if "width" is
  ! supplied, the number is padded on the left with "fill" (default '0'); when
  ! padding with '0', the minus sign of a negative number goes first (-005)

  function int2text(inter,width,fill) result(txt)
    implicit none
    integer                        :: inter
    integer , intent(in), optional :: width
    character(len=1) , optional    :: fill
    character(len=100)             :: txt0
    character(len=:) , allocatable :: txt
    character(len=1)               :: fillc
    integer                        :: nd
    write(txt0 , '(I0)') inter
    nd = len_trim(adjustl(txt0))
    if (present(width)) then
      if (width.gt.nd) then
        fillc = '0' ; if (present(fill)) fillc = fill
        allocate(character(len=width) :: txt)
        txt = repeat(fillc,width)
        if (inter.lt.0 .and. fillc.eq.'0') then
          txt(1:1) = '-'
          txt(width-nd+2:width) = txt0(2:nd)
        else
          txt(width-nd+1:width) = txt0(1:nd)
        end if
      else
        allocate(character(len=len(trim(adjustl(txt0)))) :: txt)
        txt = trim(adjustl(txt0))
      end if
      return
    else
      allocate(character(len=len(trim(adjustl(txt0)))) :: txt)
      txt = trim(adjustl(txt0))
      return
    end if
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

  ! ----------------------------------------------------------------------------
  ! this subroutine takes an array "nums" of dimension "n" and fills it with
  ! integer numbers from 1 to n and them shuffled using the the Fisher–Yates
  ! shuffle algorithm

  subroutine fisheryates(nums)
    implicit none
    integer , intent(inout) :: nums(:)
    real(dp)                :: r
    integer                 :: i,j,temp,n

    ! Initialize the array with values 1 to 100
    n = size(nums) ; do i = 1,n ; nums(i) = i ; end do

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

  subroutine shuffle_vect_int(input_vec,output_vec)
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

end module toolkit_tools

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module toolkit_econ

  ! ----------------------------------------------------------------------------
  ! econ functions and subroutines
  ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  !
  !   - crra: CRRA function
  !   - ces: CES function
  ! ----------------------------------------------------------------------------
  
  use toolkit_pars  , only : dp,cero,medio,uno,cien,tolvl
  use toolkit_tools , only : error
  implicit none

  interface ces
    module procedure ces_dos,ces_vec
  end interface ces

  contains

  ! ----------------------------------------------------------------------------
  ! crra function with "x" being the argument and "b" being the rra parameter
  ! if teh rra parameter is 0, the function returns log(x)
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

  ! ----------------------------------------------------------------------------
  ! CES function with inputs "x1" and "x2"; "b" is the (inverse) elasticity of subsitution; "a" is share of input "x1"
  function ces_dos(x1,x2,a,b) result(u)
    implicit none
    real(dp) :: x1,x2,a,b,u ; u = cero
    if (abs(b).gt.tolvl) then
      u = (a*(x1**b) + (uno-a)*(x2**b))**(uno/b)
    else if (abs(b).le.tolvl) then
      u = (x1**a)*(x2**(uno-a))
    end if
    return
  end function ces_dos
  function ces_vec(xx,aa,b) result(u)
    implicit none
    integer  :: i 
    real(dp) :: xx(:),aa(:),b,u
    u = cero
    if (sum(aa).lt.uno-tolvl .or. sum(aa).gt.uno+tolvl) then
      call error('error in ces_vec: shares do not sum to 1')
      return
    end if
    if (abs(b).gt.tolvl) then
      u = cero ; do i=1,size(xx)
        u = u + aa(i)*(xx(i)**b)
      end do
      u = u**(uno/b)
    else if (abs(b).le.tolvl) then
      u = uno; do i=1,size(xx)
        u = u*(xx(i)**aa(i))
      end do
    end if
    return
  end function ces_vec

end module toolkit_econ

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module toolkit_linalg

  ! ----------------------------------------------------------------------------
  ! linear algebra
  ! ^^^^^^^^^^^^^^
  !
  !   - cumsum: returns the vector with the cumulative sum of a vector (as matlab's cumsum function)
  !   - diag: returns the main diagonal of a matrix
  !   - eye: returns an identity matrix of size n x n
  !   - transmat: returns the transpose of a square matrix
  !   - inverse: returns the inverse of a square matrix
  ! ----------------------------------------------------------------------------
  
  use toolkit_pars  , only : dp,cero,medio,uno,cien,tolvl
  use toolkit_tools , only : error
  implicit none

  contains

  ! ----------------------------------------------------------------------------
  ! this function returns the main diagonal of a matric "mat"

  function diag(mat) result(vec)
    implicit none
    real(dp) :: mat(:,:),vec(size(mat,1))
    integer  :: i
    if (size(mat,1).ne.size(mat,2)) then
      call error(' error in diag: matrix not square')
      return
    end if
    forall (i=1:size(mat,1)) vec(i) = mat(i,i) 
    return
  end function diag

  ! ----------------------------------------------------------------------------
  ! this function returns a identity matrix of dimension nxn

  function eye(n) result(mat)
    implicit none
    integer  :: n,i
    real(dp) :: mat(n,n)
    mat = cero ; forall (i=1:n) mat(i,i) = uno
    return
  end function eye

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
  ! compute the inverse of a square matrix "m" by gauss-jordan elimination with
  ! partial pivoting. the input matrix is not modified. if the matrix is
  ! (numerically) singular, the function prints an error and returns zeros.

  function inverse(m) result(im)
    implicit none
    real(dp) , intent(in) :: m(:,:)
    real(dp)              :: im(size(m,1),size(m,1))
    real(dp)              :: a(size(m,1),size(m,1)),row(size(m,1)),piv,coeff,scl
    integer               :: n,i,k,p

    im = cero
    if (size(m,1).ne.size(m,2)) then
      call error(' error in inverse: matrix is not square')
      return
    end if
    n = size(m,1) ; a = m
    do i=1,n
      im(i,i) = uno
    end do
    scl = maxval(abs(m))
    do k=1,n
      ! pivot: largest element (in absolute value) of column k, rows k to n
      p   = maxloc(abs(a(k:n,k)),dim=1) + k - 1
      piv = a(p,k)
      if (abs(piv).le.epsilon(uno)*dble(n)*scl .or. scl.le.cero) then
        call error(' error in inverse: matrix is singular')
        im = cero
        return
      end if
      if (p.ne.k) then
        row = a(k,:)  ; a(k,:)  = a(p,:)  ; a(p,:)  = row
        row = im(k,:) ; im(k,:) = im(p,:) ; im(p,:) = row
      end if
      a(k,:)  = a(k,:)/piv
      im(k,:) = im(k,:)/piv
      do i=1,n
        if (i.ne.k) then
          coeff = a(i,k)
          if (coeff.ne.cero) then
            a(i,:)  = a(i,:)  - coeff*a(k,:)
            im(i,:) = im(i,:) - coeff*im(k,:)
          end if
        end if
      end do
    end do
    return
  end function inverse

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

end module toolkit_linalg

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module toolkit_stats

  ! ----------------------------------------------------------------------------
  ! statistics
  ! ^^^^^^^^^^
  !
  !   - varmean: returns the average of a variable, allowing for weights
  !   - varvar: returns the variance of a variable, allowing for weights
  !   - varstd: returns the standard deviation of a variable, allowing for weights
  !   - correlation: returns the correlation of two variables, allowing for weights
  !   - percentile: returns the i-th percentile of a variable, allowing for weights
  !   - olsreg: returns the ols coefficients of a regression with up to 8 regressors, allowing for weights
  !   - tauchen: returns the transition matrix for a discretized ar(1) process
  !   - normaldist: returns the distribution of a normal variable discretized over a grid
  !   - randomnormal: returns a random draw from a normal distribution
  !   - cdfn: returns the cdf of a standard normal distribution
  ! ----------------------------------------------------------------------------
  
  use toolkit_pars   , only : dp,cero,medio,uno,cien,tolvl
  use toolkit_tools  , only : error,num2text
  use toolkit_linalg , only : diag,inverse,transmat
  implicit none

  interface randomnormal
    module procedure randomnormal_scalar,randomnormal_vec
  end interface randomnormal

  contains

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
        call error('error in varmean: mask of incorrect size')
        return
      end if
    end if

    if (present(w)) then
      if (size(var).ne.size(w)) then
        call error('error in varmean!! var and w have different size')
        return
      else
        weig(:) = w(:)
      end if
      if (sum(weig).lt.tolvl) then
        call error('error in varmean!! w are cero')
        return
      end if
    end if

    meanvar = sum(var(:)*weig(:),mask=mask1)/sum(weig,mask=mask1)

    deallocate(weig)
    deallocate(mask1)

    return
  end function varmean

  ! ----------------------------------------------------------------------------
  ! this function returns the variance of a variable "var" (observations with
  ! var = nan are excluded)

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
        call error('error in varvar: mask of incorrect size')
        return
      end if
    end if

    if (present(w)) then
      if (size(var).ne.size(w)) then
        call error('error in varvar!! var and w have different size')
        return
      else
        weig(:) = w(:)
      end if
      if (sum(weig).lt.tolvl) then 
        call error('error in varvar!! w are cero')
        return
      end if
    end if

    ! the same observations (valid and not nan) for the mean and the variance
    mask1  = mask1 .and. .not.isnan(var)
    mvar   = varmean(var,weig,mask=mask1)
    varian = sum(weig*((var-mvar)**dos),mask=mask1)/sum(weig,mask=mask1)

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
        call error('error in varstd: mask of incorrect size')
        return
      end if
    end if

    if (present(w)) then
      if (size(var).ne.size(w)) then
        call error('error in varstd!! var and w have different size')
        return
      else
        weig(:) = w(:)
      end if
      if (sum(weig).lt.tolvl) then
        call error('error in varstd!! w are cero')
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
        call error('error in correlation: mask of incorrect size')
        return
      end if
    end if

    if (size(xvar1).ne.size(xvar2)) then
      call error('error in correlation!! xvar1 and xvar2 of different sizes')
      return
    end if
    if (present(w)) then
      if (size(xvar1).ne.size(w)) then
        call error('error in correlation!! var and w have different size')
        return
      end if
      if (sum(w).lt.tolvl) then
        call error('error in correlation!! w are cero')
        return
      end if
      weig(:) = w(:)
    end if

    aux1 = varmean(xvar1,weig,mask=mask1)
    aux2 = varstd(xvar1,weig,mask=mask1)
    aux3 = varmean(xvar2,weig,mask=mask1)
    aux4 = varstd(xvar2,weig,mask=mask1)

    if (aux2.lt.tolvl) then
      call error('error in correlation!! xvar1 is a constant')
      return
    end if
    if (aux4.lt.tolvl) then
      call error('error in correlation!! xvar2 is a constant')
      return
    end if

    aux5 = sum(weig(:)*( (xvar1(:) - aux1)*(xvar2(:) - aux3) ),mask=mask1)/sum(weig,mask=mask1)
    corr = aux5/(aux2*aux4)

    deallocate(weig)
    deallocate(mask1)

    return
  end function correlation

  ! ----------------------------------------------------------------------------
  ! this function returns the percentile "pct" for a distribution "xvec".
  ! "pct" can be given as a fraction (0 <= pct <= 1) or as a percentage
  ! (1 < pct <= 100). note that pct = 1 is read as a fraction (the maximum):
  ! for the 1st percentile use pct = 0.01.

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
        call error('error in percentile: mask of incorrect size')
        return
      end if
    end if

    if (present(w)) then
      if (size(w).ne.size(xvec)) then
        call error('error in percentile!! var and w have different size')
        return
      end if
      if (sum(w).lt.tolvl) then
        call error('error in percentile!! w are cero')
        return
      end if
      weig(:) = w(:)
    end if

    if (pct.lt.cero) then
      call error('error in percentile!! invalid percentile: negative value')
      return
    end if
    if (pct.gt.uno) then
      if (pct.le.cien) then
        perc = pct/cien
      else
        call error('error in percentile!! invalid percentile: larger than 100')
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
    integer                           :: i,n0,no,nc,wc,nx ; nx = 0

    n0     = size(yvec)
    coeffs = -huge(uno)

    ! mask – valid observations
    allocate(mask1(n0)) ; mask1 = .true.
    if (present(mask)) then
      if (size(mask).ne.n0) then
        call error('error in olsreg!! mask of incorrect size')   
        return
      else
        mask1 = mask
      end if
    end if

    ! number of valid observations
    no = count(mask1)  
    if (no.eq.0) then 
      call error('error in olsreg!! no (valid) observations')
      return
    end if

    ! weigths
    allocate(wvar(n0)) ; wvar = uno
    if (present(w)) then
      if (size(w).ne.n0) then
        call error('error in olsreg!! weigths of incorrect size')
        return
      else
        wvar = w
      end if
    end if

    ! check number of explanatory variabkes
    if (present(x1vec)) then ; nx = 1
      if (n0.ne.size(x1vec)) then
        call error('error in olsreg!! yvec and x1vec different observations')
        return
      end if
    end if
    if (present(x2vec)) then ; nx = 2
      if (n0.ne.size(x2vec)) then
        call error('error in olsreg!! yvec and x2vec different observations')
        return
      end if
    end if
    if (present(x3vec)) then ; nx = 3
      if (n0.ne.size(x3vec)) then
        call error('error in olsreg!! yvec and x3vec different observations')
        return
      end if
    end if
    if (present(x4vec)) then ; nx = 4
      if (n0.ne.size(x4vec)) then
        call error('error in olsreg!! yvec and x4vec different observations')
        return
      end if
    end if
    if (present(x5vec)) then ; nx = 5
      if (n0.ne.size(x5vec)) then
        call error('error in olsreg!! yvec and x5vec different observations')
        return
      end if
    end if
    if (present(x6vec)) then ; nx = 6
      if (n0.ne.size(x6vec)) then
        call error('error in olsreg!! yvec and x6vec different observations')
        return
      end if
    end if
    if (present(x7vec)) then ; nx = 7
      if (n0.ne.size(x7vec)) then
        call error('error in olsreg!! yvec and x7vec different observations')
        return
      end if
    end if
    if (present(x8vec)) then ; nx = 8
      if (n0.ne.size(x8vec)) then
        call error('error in olsreg!! yvec and x8vec different observations')
        return
      end if
    end if
    
    nc = size(coeffs)  ! number of coefficients
    wc = nc-nx         ! with constant = 1, 0 if not

    if (nc.ge.no) then 
      call error('error in olsreg!! '//num2text(nc)//' coefficients and '//num2text(no)//' observations')
      return
    end if

    if (wc.ne.1 .and. wc.ne.0) then 
      call error('error in olsreg!! '//num2text(nc)//' coefficients and '//num2text(nx)//' variables')
      return
    end if

    if (present(stds)) then
      if (size(stds).ne.size(coeffs)) then
        call error('error in olsreg!! stds and coeffs of different size')
        return
      end if
    end if

    ! allocate y-variable and weigths
    allocate(xvars(no,nc)) ; xvars = cero
    allocate(yvar(no))     ; yvar = pack(yvec,mask1)
    allocate(zvar(no))     ; zvar = pack(wvar,mask1)
    
    ! check weigths are positive    
    if (sum(zvar).lt.tolvl) then
      call error('error in olsreg!! weigths are cero')
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

        inter(:,1) = coeffs - dble(1.959964)*sdcoefs
        inter(:,2) = coeffs + dble(1.959964)*sdcoefs

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
    dist(1) = cdfn( medio*(xvec0(2)+xvec0(1)) )
    do j=2,n-1
      dist(j) = cdfn( medio*(xvec0(j+1)+xvec0(j)) ) - cdfn( medio*(xvec0(j)+xvec0(j-1)) )
    end do
    dist(n) = 1.0 - cdfn( medio*(xvec0(n)+xvec0(n-1)) )
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
      if (v**2 < -4.0d0*log(u)*u**dos) exit
    end do
    shock = std*v/u + mu
    return
  end subroutine randomnormal_scalar
  ! vector of independent draws (the sample is not re-centered, so its sample
  ! mean is not exactly "mu")
  subroutine randomnormal_vec(shock,mu,std)
    implicit none
    integer                :: i
    real(dp) , intent(out) :: shock(:)
    real(dp) , intent(in)  :: mu,std
    do i=1,size(shock)
      call randomnormal_scalar(shock(i),mu,std)
    end do
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

end module toolkit_stats

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module toolkit_optim

  ! ----------------------------------------------------------------------------
  ! optimization
  ! ^^^^^^^^^^^^
  !
  ! this code contains 5 optimization and root-finding algorithms:
  !
  !   - bisection: root of a single-value function with one unknown (bisection method)
  !   - golden: maximum of a single-value function with one unknown (golden search)
  !   - brent: root of a single-value function with one unknown (brent method)
  !   - simplex: minimum of a single-value function with "n" unknowns (nelder-mead)
  !   - lmmin: minimum of the sum of squares of a system of "m" functions with "n"
  !     unknowns (levenberg–marquardt algorithm)
  !
  ! it also includes some other functions/subroutines used for optimization:
  !
  !   - broyden: updates a jacobian matrix using the broyden's method
  !   - normalize: transforms a bounded variable into an unbounded one
  !   - denormalize: transforms an unbounded variable into a bounded one
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
  ! ----------------------------------------------------------------------------

  use toolkit_pars  , only : dp,cero,medio,uno,dos,cien,tolvl
  use toolkit_tools , only : error,num2text

  implicit none
  contains

  ! ----------------------------------------------------------------------------
  ! bisection's method

  subroutine bisection(func,x,iy,ind,x0,x1,itermax,tol,iprint)

    ! this subroutine finds the root of a user-supplied single-valued function
    ! with one unknown using the bisection method
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
    real(dp)                         :: toler,xa,ya,xb,yb,xc,yc

    interface
      function func(xx) result(resid)
        real(kind=8) :: xx,resid
      end function func
    end interface

    ! - ```exitcode``` = 0: the algorithm found a root
    ! - ```exitcode``` = 10: ```x0``` is a root of ```func```
    ! - ```exitcode``` = 11: ```x1``` is a root of ```func```
    ! - ```exitcode``` = 2: the root is not within the interval (```x0```, ```x1```)
    ! - ```exitcode``` = 3: the points ```x0``` and ```x1``` are too close
    ! - ```exitcode``` = 9: maximum number of function evaluations reached

    toler   = 1.0d-8 ; if (present(tol)    ) toler   = tol
    maxiter = 500    ; if (present(itermax)) maxiter = itermax
    ipri    = 0      ; if (present(iprint) ) ipri    = iprint

    if (ipri.gt.0) then
      write(*,*) '  '
      write(*,*) ' starting bisection algorithm'
      write(*,*) '  '
    end if

    if (abs(x1-x0).lt.toler) goto 91

    ! compute lower bound
    xa = x0 ; ya = func(xa)
    if (ipri.gt.0) write(*,*) '  lower bound: x0 = '//num2text(xa,2)//' | f(x0) = '//num2text(ya,2) 
    if (abs(ya).lt.toler) then
      x = xa ; goto 10 ; return
    end if
    xb = x1 ; yb = func(xb)
    if (ipri.gt.0) write(*,*) '  upper bound: x1 = '//num2text(xb,2)//' | f(x1) = '//num2text(yb,2) 
    if (abs(yb).lt.toler) then
      x = xb ; goto 11 ; return
    end if

    if (ya.gt.cero .and. yb.lt.cero) then
      xc = xa ; yc = ya ; xa = xb ; ya = yb ; xb = xc ; yb = yc
    else if (ya.gt.cero .and. yb.gt.cero) then
      goto 93
    elseif (yb.lt.cero .and. ya.lt.cero) then
      goto 94
    end if

    if (ipri.gt.1) then
      write(*,*) '                          '
      write(*,*) '  beginning iterations... '
      write(*,*) '                          '
    end if

    ! beginning iterations
    iy = 0 ; ind = 0

    do while (abs(xb-xa).gt.toler .and. iy.lt.maxiter)
      xc = medio*(xa+xb) ; yc = func(xc) ; iy = iy + 1 
      if (ipri.gt.1) then
        write(*,87) iy,xa,xb,xc,ya,yb,yc
      end if
      if (yc.gt.cero) then
        xb = xc ; yb = yc
      else
        xa = xc ; ya = yc
      end if
      if (abs(yc).lt.toler) then
        x = xc ; goto 20
      end if
      if (iy.eq.maxiter) then
        x = xc ; goto 99
      end if
    end do

    if (ipri.gt.0) then
      write(*,86) iy,xb,yb
    end if

    x = xb

    return

    91 ind = 3  ; if (ipri.gt.0) write(*,*) '  error: x0 and x1 are too close'            ; return
    92 ind = 2  ; if (ipri.gt.0) write(*,*) '  error: f(x0) and f(x1) have the same sign' ; return
    93 ind = 20 ; if (ipri.gt.0) write(*,*) '  error: f(x0) > 0 '                         ; return
    94 ind = 21 ; if (ipri.gt.0) write(*,*) '  error: f(x1) < 0 '                         ; return
    10 ind = 10 ; if (ipri.gt.0) write(*,*) '  solved: x0 is a root'                      ; return
    11 ind = 11 ; if (ipri.gt.0) write(*,*) '  solved: x1 is a root'                      ; return
    99 ind = 9  ; if (ipri.gt.0) write(*,*) '  max number of iterations reached'          ; return
    20 ind = 0  ; if (ipri.gt.0) write(*,*) '  solved: root found'                        ; return
         
    87 format( '  iteration = ',i4,' | xa = ',f10.4,' | xb = ',f10.4,' | xc = ',f10.4,&
    ' | f(xa) = ',f10.4,' | f(xb) = ',f10.4,' | f(xc) = ',f10.4)
    86 format( '  solved: iterations = ',i4,' | x = ',f10.4,' | f(x) = ',f10.4)
  end subroutine bisection

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
    real(dp) , parameter             :: alpha = (sqrt(5.d0) - 1.d0) / 2.d0
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
    x2 = alpha*x1 + (uno-alpha)*x3   ! = x0 + alpha*(x3-x0)
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
      if (ip.gt.1) write(*,15) it , x0 , x1 , x2 , x3 , f1 , f2
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
    15 format( '  iteration = ',i4,' | x0 = ',f10.4,' | x1 = ',f10.4,&
    ' | x2 = ',f10.4,' | x3 = ',f10.4,' | f(x1) = ',f10.4,' | f(x2) = ',f10.4)
    16 format( '  solved: iterations = ',i4,' | x = ',f10.4,' | f(x) = ',f10.4)
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
    ! - ```exitcode``` = 10: ```x0``` is a root of ```func```
    ! - ```exitcode``` = 11: ```x1``` is a root of ```func```
    ! - ```exitcode``` = 2: the root is not within the interval (```x0```, ```x1```)
    ! - ```exitcode``` = 3: the points ```x0``` and ```x1``` are too close
    ! - ```exitcode``` = 9: maximum number of function evaluations reached

    toler   = 1.0d-8 ; if (present(tol)    ) toler   = tol
    maxiter = 500    ; if (present(itermax)) maxiter = itermax
    ipri    = 0      ; if (present(iprint) ) ipri    = iprint

    if (ipri.gt.0) then
      write(*,*) '  '
      write(*,*) ' starting brent algorithm'
      write(*,*) '  '
    end if

    if (abs(x1-x0).lt.toler) goto 91

    ! compute lower bound
    xa = x0 ; ya = func(xa)
    if (ipri.gt.0) write(*,*) '  lower bound: x0 = '//num2text(xa,2)//' | f(x0) = '//num2text(ya,2) 
    if (abs(ya).lt.toler) then
      x = xa ; goto 10 ; return
    end if
    if (ya.gt.cero) goto 93

    ! compute upper bound
    xb = x1 ; yb = func(xb)
    if (ipri.gt.0) write(*,*) '  upper bound: x1 = '//num2text(xb,2)//' | f(x1) = '//num2text(yb,2) 
    if (abs(yb).lt.toler) then
      x = xb ; goto 11 ; return
    end if
    if (yb.lt.cero) goto 94

    ! check that root is within bounds
    if (ya*yb.gt.cero) goto 92

    ! reorder if bounds are reversed
    if (ya.gt.cero .and.yb.lt.cero) then ; 
      xc = xa ; yc = ya ; xa = xb ; ya = yb ; xb = xc ; yb = yc
    end if

    ! beginning iterations
    xc = medio*(xa+xb) ; iy = 0 ; ind = 0

    if (ipri.gt.1) then
      write(*,*) '                          '
      write(*,*) '  beginning iterations... '
      write(*,*) '                          '
    end if

    do while (abs(xb-xa).gt.toler .and. iy.lt.maxiter)

      yc = func(xc)

      if (abs(yc).lt.toler) then
        x = xc ; goto 20
      end if

      if (abs(yc-ya).gt.toler .and. abs(yc-yb).gt.toler) then
        xs = xa*yb*yc/((ya-yb)*(ya-yc)) + xb*ya*yc/((yb-ya)*(yb-yc)) + xc*ya*yb/((yc-ya)*(yc-yb))
      else
        xs = xb - yb*(xb-xa)/(yb-ya)
      end if

      ys = func(xs) ; iy = iy + 1

      if (ipri.gt.1) then
        write(*,*) '  iter = '//num2text(iy)//' | xs = '//num2text(xs,2)//' | f(xs) = '//num2text(ys,2) 
      end if

      if (ya*ys.lt.cero) then
        xc = xb ; xb = xs ; yb = ys
      else
        xa = xs ; ya = ys
      end if
      if (abs(ya).lt.abs(yb)) then
        xc = xb ; xb = xa ; xa = xc ; yc = yb ; yb = ya ; ya = yc
      end if
      if (abs(ys).lt.toler) then
        x = xs ; goto 20
      end if
      if (iy.eq.maxiter) then
        x = xs ; goto 99
      end if

    end do

    if (ipri.gt.0) then
      write(*,86) iy,xb,yb
    end if

    x = xb

    return

    91 ind = 3  ; if (ipri.gt.0) write(*,*) '  error: x0 and x1 are too close'            ; return
    92 ind = 2  ; if (ipri.gt.0) write(*,*) '  error: f(x0) and f(x1) have the same sign' ; return
    93 ind = 20 ; if (ipri.gt.0) write(*,*) '  error: f(x0) > 0 '                         ; return
    94 ind = 21 ; if (ipri.gt.0) write(*,*) '  error: f(x1) < 0 '                         ; return
    10 ind = 10 ; if (ipri.gt.0) write(*,*) '  solved: x0 is a root'                      ; return
    11 ind = 11 ; if (ipri.gt.0) write(*,*) '  solved: x1 is a root'                      ; return
    99 ind = 9  ; if (ipri.gt.0) write(*,*) '  max number of iterations reached'          ; return
    20 ind = 0  ; if (ipri.gt.0) write(*,*) '  solved: root found'                        ; return
         
    86 format( '  solved: iterations = ',i4,' | x = ',f10.4,' | f(x) = ',f10.4)
  end subroutine brent

  ! ----------------------------------------------------------------------------
  ! simplex algorithm

  subroutine simplex(func,x,y,iy,ind,x0,itermax,tol,iprint,tolx,step,fzero)

    ! this subroutine minimizes a single-valued function with "n" unknowns using the
    ! nelder-mead simplex algorithm, with the adaptive coefficients of gao and han
    ! (2012, computational optimization and applications 51, 259-277):
    !   reflection 1, expansion 1+2/n, contraction 0.75-1/(2n), shrink 1-1/n
    ! (with n = 1 the standard coefficients 1, 2, 0.5, 0.5 are used)
    !
    ! the inputs are:
    !
    !   - func: user-provided function to be minimize. this function should of the form:
    !
    !       function func(x) result(y)
    !         real(dp) :: x(:),y
    !       end function func
    !
    !   - x0: the vector of dimension n with the initial guess
    !
    ! there are also few optional inputs:
    !
    !   itermax:  max number of function evaluations (default = 500)
    !   tol:      tolerance level (default = 1.0d-8):
    !               - spread of the function values in the simplex
    !               - if fzero = .true., value of func that is considered zero
    !   tolx:     tolerance on the size of the simplex (default = sqrt(tol))
    !   step:     size of the initial simplex (default = 0.20)
    !   fzero:    .true. (default): stop as soon as func(x) < tol, as when func is a
    !             sum of squared residuals. use .false. for general minimization
    !             (e.g. if func can be negative)
    !   iprint:   indicator for printing results (default = 0)
    !             0: print no results
    !             1: print main results
    !             2: print all the iterations
    !
    ! the outputs of the subroutine are:
    !
    !   x:    the vector of x that minimizes func
    !   y:    the value of the function at the point x
    !   iy:   the number of total function evaluations
    !   ind:  exit indicator, taking value:
    !           0: func(x) < tol (only if fzero = .true.)
    !           1: converged: spread of func < tol and size of the simplex < tolx
    !              (or size of the simplex < tol)
    !           9: max number of function evaluations reached
    !
    ! november 2026 (Borja Petit): corrected the expansion point, previously
    ! (1-g)*xm - g*xhigh instead of xm + g*(xr-xm), and the contraction point,
    ! previously (1-b)*xm - b*xhigh instead of xm + b*(xhigh-xm); added the outside
    ! contraction, the adaptive coefficients, the convergence criterion on the
    ! spread of func, and the optional inputs tolx, step and fzero.

    implicit none

    integer  , intent(out)           :: iy,ind
    real(dp) , intent(out)           :: x(:),y

    real(dp) , intent(in)            :: x0(:)
    real(dp) , intent(in) , optional :: tol,tolx,step
    integer  , intent(in) , optional :: itermax
    integer  , intent(in) , optional :: iprint
    logical  , intent(in) , optional :: fzero

    integer  :: n,i,ib,iw,iw2,maxiter,ipri
    real(dp) :: xp(size(x0,1),size(x0,1)+1),yp(size(x0,1)+1)
    real(dp) :: xm(size(x0,1)),xr(size(x0,1)),xe(size(x0,1)),xc(size(x0,1))
    real(dp) :: y0,yr,ye,yc,toler,tolerx,stp,size0,ca,cg,cr,cs
    logical  :: fz,shrink

    interface
      function func(xvar) result(resid)
        import :: dp
        real(kind=dp) :: xvar(:),resid
      end function func
    end interface

    toler   = 1.0d-8      ; if (present(tol)    ) toler   = tol
    tolerx  = sqrt(toler) ; if (present(tolx)   ) tolerx  = tolx
    stp     = 0.20_dp     ; if (present(step)   ) stp     = step
    maxiter = 500         ; if (present(itermax)) maxiter = itermax
    ipri    = 0           ; if (present(iprint) ) ipri    = iprint
    fz      = .true.      ; if (present(fzero)  ) fz      = fzero

    n  = size(x0,1)
    ca = uno                                        ! reflection
    cg = uno + dos/max(dble(n),dos)                 ! expansion
    cr = 0.75_dp - uno/(dos*max(dble(n),dos))       ! contraction
    cs = uno - uno/max(dble(n),dos)                 ! shrink

    if (ipri.ge.1) then
      write(*,*) '  '
      write(*,*) ' starting simplex algorithm'
      write(*,*) '  '
    end if

    ind = 9 ; xp = cero ; yp = huge(uno)

    ! initial point
    y0 = func(x0) ; iy = 1 ; xp(:,1) = x0 ; yp(1) = y0
    if (ipri.ge.1) print 93 , y0
    if (fz .and. y0.lt.toler) then
      ind = 0 ; goto 7
    end if

    ! initial simplex
    if (ipri.gt.1) then
      write(*,*) '  '
      write(*,*) ' evaluating initial simplex'
      write(*,*) '  '
    end if
    do i=1,n
      xp(:,i+1) = x0
      if (x0(i).lt.cero) xp(i,i+1) = x0(i) - stp
      if (x0(i).ge.cero) xp(i,i+1) = x0(i) + stp
      yp(i+1) = func(xp(:,i+1)) ; iy = iy + 1
      if (ipri.gt.1 .and. n.lt.5) write(*,97) i, yp(i+1) , xp(:,i+1)
      if (ipri.gt.1 .and. n.ge.5) write(*,97) i, yp(i+1)
      if (fz .and. yp(i+1).lt.toler) then
        ind = 0 ; goto 7
      end if
    end do

    if (ipri.gt.1) then
      write(*,*) '  '
      write(*,*) ' beginning iterations'
      write(*,*) '  '
    end if

    do

      ! best (ib), worst (iw) and second worst (iw2) points
      ib = minloc(yp,dim=1)
      iw = maxloc(yp,dim=1)
      iw2 = ib
      do i=1,n+1
        if (i.ne.iw .and. yp(i).ge.yp(iw2)) iw2 = i
      end do

      if (ipri.gt.1) write(*,97) iy, yp(ib)

      ! convergence
      size0 = maxval(abs(xp-spread(xp(:,ib),2,n+1)))
      if (size0.lt.toler .or. (yp(iw)-yp(ib).lt.toler .and. size0.lt.tolerx)) then
        ind = 1 ; exit
      end if

      ! max number of function evaluations
      if (iy.ge.maxiter) then
        ind = 9 ; exit
      end if

      ! centroid of all points but the worst
      xm = (sum(xp,dim=2) - xp(:,iw))/dble(n)

      ! reflection
      xr = xm + ca*(xm-xp(:,iw)) ; yr = func(xr) ; iy = iy + 1
      shrink = .false.

      if (yr.lt.yp(ib)) then
        ! expansion
        xe = xm + cg*(xr-xm) ; ye = func(xe) ; iy = iy + 1
        if (ye.lt.yr) then
          xp(:,iw) = xe ; yp(iw) = ye
        else
          xp(:,iw) = xr ; yp(iw) = yr
        end if
      elseif (yr.lt.yp(iw2)) then
        ! accept reflection
        xp(:,iw) = xr ; yp(iw) = yr
      elseif (yr.lt.yp(iw)) then
        ! outside contraction
        xc = xm + cr*(xr-xm) ; yc = func(xc) ; iy = iy + 1
        if (yc.le.yr) then
          xp(:,iw) = xc ; yp(iw) = yc
        else
          shrink = .true.
        end if
      else
        ! inside contraction
        xc = xm + cr*(xp(:,iw)-xm) ; yc = func(xc) ; iy = iy + 1
        if (yc.lt.yp(iw)) then
          xp(:,iw) = xc ; yp(iw) = yc
        else
          shrink = .true.
        end if
      end if

      ! shrink towards the best point
      if (shrink) then
        do i=1,n+1
          if (i.ne.ib) then
            xp(:,i) = xp(:,ib) + cs*(xp(:,i)-xp(:,ib))
            yp(i)   = func(xp(:,i)) ; iy = iy + 1
          end if
        end do
      end if

      ! func is zero
      if (fz .and. minval(yp).lt.toler) then
        ind = 0 ; exit
      end if

    end do

    7 continue

    ib = minloc(yp,dim=1)
    x  = xp(:,ib)
    y  = yp(ib)

    if (ipri.ge.1) then
      write(*,*) '  '
      if (ind.eq.0) write(*,*) ' simplex finished: system has converged, func=0 '
      if (ind.eq.1) write(*,*) ' simplex finished: converged, simplex is very small'
      if (ind.eq.9) write(*,*) ' simplex finished: max number of iterations'
      write (*, *) '  '
      write (*,94) iy
      write (*,93) y
      write (*,92) y/y0 - 1.0
      write (*, *) '  '
    end if

    return
    97 format ('  new point    ',i4,'  ',5(f10.4))
    94 format ('  iterations = ',i4)
    93 format ('  error      = ',f10.4)
    92 format ('  reduction  = ',f10.4)

  end subroutine simplex

  ! ----------------------------------------------------------------------------
  ! levenberg–marquardt algorithm

  subroutine lmmin(func,x,y,iy,ind,x0,itermax,damp,tol,toleach,shock,usebro,iprint,tolrel,maxstall)

    ! this subroutine minimizes the sum of squared errors of a system of "m" equations
    ! in "n" unknowns using the levenberg-marquardt algorithm. given the current point
    ! xb, with residuals yb and jacobian J, the step dx solves
    !
    !     ( J'·J + lambda·D ) · dx = - J'·yb
    !
    ! where D is a diagonal scaling matrix (the largest diagonal of J'·J observed so far
    ! for each unknown, with a positive floor; more, 1978) and lambda is the damping
    ! factor. when lambda is large, the step is a short gradient-descent step; when it
    ! is small, it is a gauss-newton step. lambda is updated with the gain ratio
    ! rho = (actual reduction)/(predicted reduction) as in nielsen (1999):
    !   - if the step reduces the error: lambda = lambda·max(1/3, 1-(2·rho-1)^3)
    !   - otherwise: lambda = lambda·nu, nu = 2·nu
    ! the linear system is solved by cholesky decomposition.
    !
    ! the inputs are:
    !
    !   - func: user-provided function. this function should of the form:
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
    !   itermax:  max number of function evaluations (default = 500)
    !   damp:     initial damping factor lambda (default = 1)
    !   tol:      tolerance level for the sum of squared errors (default = 1.0d-8)
    !   toleach:  tolerance level for each equation, max|y| (default = 1.0d-8)
    !   shock:    relative shock to compute the numerical jacobian (default = 0.05)
    !   usebro:   whether to use broyden's method to update jacobian (default = 0)
    !               if broyden's method is used, the program computes the numerical
    !               jacobian once, and then updates it after every successful step;
    !               a new numerical jacobian is computed when a broyden-updated
    !               jacobian fails to produce a better point
    !   tolrel:   relative improvement in the error below which an iteration is
    !             considered unproductive (default = 1.0d-6)
    !   maxstall: number of consecutive unproductive iterations after which the
    !             algorithm stops (default = 10)
    !   iprint:   indicator for printing results (default = 0)
    !               0: print no results
    !               1: print main results
    !               2: print all the iterations
    !
    ! the outputs of the subroutine are:
    !
    !   x:    vector with the values of x with the lowest error found
    !   y:    the value of the function at the point x
    !   iy:   the number of total function evaluations
    !   ind:  exit indicator, taking value:
    !           0: system has converged, func(x) = 0
    !           1: gradient J'·y is close to 0 (or the linear system cannot be solved)
    !           2: step in x is close to 0
    !           3: no progress: maxstall consecutive unproductive iterations
    !           9: max number of function evaluations reached
    !
    ! november 2026 (Borja Petit): (1) the step used to compute a column of the
    ! jacobian is enlarged at most 5 times when the residuals do not change (it was
    ! enlarged without limit, until itermax, and the test was on the sum of squared
    ! errors instead of the residuals); if the residuals still do not change the
    ! column is set to 0; (2) scaling with a positive floor and cholesky solve instead
    ! of the inverse of (1+lambda)·diag(J'·J), which was singular with a zero column;
    ! (3) damping updated with the gain ratio; (4) stop after maxstall unproductive
    ! iterations (ind = 3) and relative test on the step; (5) points evaluated to
    ! compute the jacobian are used when they improve on the current point, and the
    ! best point evaluated is returned; (6) documentation and minor fixes.

    implicit none

    ! outputs
    real(dp) , intent(out)           :: x(:)
    real(dp) , intent(out)           :: y(:)
    integer  , intent(out)           :: iy
    integer  , intent(out)           :: ind

    ! inputs (including optionals)
    real(dp) , intent(in)            :: x0(:)
    real(dp) , intent(in) , optional :: damp,tol,toleach,shock,tolrel
    integer  , intent(in) , optional :: itermax,usebro,iprint,maxstall

    ! other variables
    integer                :: n,m,i,maxiter,ip,bro,mstall,nstall,itry
    real(dp)               :: lam,nu,toler,tolereach,tolr,shck0
    real(dp)               :: e0,eb,e1,ec,es,pred,rho,dfloor
    logical                :: newjac   ! compute a new numerical jacobian
    logical                :: brojac   ! jacobian has been updated with broyden
    logical                :: done     ! exit flag set inside the jacobian routine
    logical                :: ok
    real(dp) , allocatable :: j(:,:)          ! jacobian
    real(dp) , allocatable :: a(:,:)          ! J'·J
    real(dp) , allocatable :: g(:)            ! J'·y
    real(dp) , allocatable :: dg(:)           ! scaling: max diagonal of J'·J so far
    real(dp) , allocatable :: shck(:)         ! relative shock for each unknown
    real(dp) , allocatable :: xb(:),yb(:)     ! current point
    real(dp) , allocatable :: x1(:),y1(:)     ! trial point
    real(dp) , allocatable :: dx(:)           ! step
    real(dp) , allocatable :: xc(:),yc(:)     ! best point evaluated for the jacobian
    real(dp) , allocatable :: xs(:),ys(:)     ! best point evaluated overall
    real(dp) , allocatable :: xj(:),yj(:)     ! point evaluated for the jacobian
    real(dp)               :: ej              ! error at xj
    integer                :: ij              ! unknown being shocked

    interface
      function func(xvec) result(resid)
        real(kind=kind(1.0d00))               :: xvec(:)
        real(kind=kind(1.0d00)) , allocatable :: resid(:)
      end function func
    end interface

    ! initialize optional inputs
    lam       = 1.0_dp  ; if (present(damp)    ) lam       = damp
    toler     = 1.0d-8  ; if (present(tol)     ) toler     = tol
    tolereach = 1.0d-8  ; if (present(toleach) ) tolereach = toleach
    shck0     = 5.0d-2  ; if (present(shock)   ) shck0     = shock
    tolr      = 1.0d-6  ; if (present(tolrel)  ) tolr      = tolrel
    maxiter   = 500     ; if (present(itermax) ) maxiter   = itermax
    ip        = 0       ; if (present(iprint)  ) ip        = iprint
    bro       = 0       ; if (present(usebro)  ) bro       = usebro
    mstall    = 10      ; if (present(maxstall)) mstall    = maxstall

    n = size(x,1) ; m = size(y,1)

    if (size(x0,1).ne.n) then
      write(*,*) ' error in lmmin: size of x0 does not match size of x '
      stop
    end if

    allocate(j(m,n),a(n,n),g(n),dg(n),shck(n),xb(n),yb(m),x1(n),y1(m),dx(n))
    allocate(xc(n),yc(m),xs(n),ys(m),xj(n),yj(m))
    shck = shck0 ; dg = cero ; nu = dos ; nstall = 0 ; ind = 9

    if (ip.ge.1) then
      write(*,'(a)'   ) '                              '
      write(*,'(a)'   ) '  starting lmmin algorithm    '
      write(*,'(a)'   ) '                              '
      write(*,'(a,i4)') '  # equations = ' , m
      write(*,'(a,i4)') '  # unknowns  = ' , n
      write(*,'(a)'   ) '                              '
    end if

    ! initial point
    xb = x0 ; yb = func(xb) ; eb = sum(yb*yb) ; e0 = eb ; iy = 1
    xs = xb ; ys = yb ; es = eb ; ec = huge(uno)
    if (ip.ge.1) write (*,91) eb

    if (lm_converged(eb,yb)) then
      ind = 0 ; goto 7
    end if

    newjac = .true. ; brojac = .false.

    do

      if (iy.ge.maxiter) then
        ind = 9 ; exit
      end if

      ! numerical jacobian at xb
      if (newjac) then
        call lm_jacobian(done)
        if (done) exit
        newjac = .false. ; brojac = .false.
      end if

      ! gradient
      g = matmul(transpose(j),yb)
      if (maxval(abs(g)).lt.toler) then
        ind = 1 ; exit
      end if

      ! scaling (more, 1978) with a positive floor
      a = matmul(transpose(j),j)
      do i=1,n
        dg(i) = max(dg(i),a(i,i))
      end do
      dfloor = 1.0d-10*max(maxval(dg),uno)

      ! step: cholesky solve; if it fails, increase damping
      ok = .false.
      do itry=1,30
        call lm_cholsolve(a,max(dg,dfloor),lam,g,dx,ok)
        if (ok) exit
        lam = lam*10.0_dp
      end do
      if (.not.ok) then
        ind = 1 ; exit
      end if

      ! if step in x is too small, finish
      if (sum(abs(dx)).le.toler*(uno+sum(abs(xb)))) then
        ind = 2 ; exit
      end if

      ! predicted and actual reduction
      pred = eb - sum((yb+matmul(j,dx))**2)
      x1 = xb + dx
      y1 = func(x1) ; e1 = sum(y1*y1) ; iy = iy + 1
      call lm_keepbest(x1,y1,e1)
      rho = (eb-e1)/max(pred,tiny(uno))

      if (ip.ge.1) print 97 , iy , e1 , eb , lam

      if (e1.lt.eb) then

        ! good step: update damping, the jacobian, and the point
        lam = lam*max(uno/3.0_dp,uno-(dos*rho-uno)**3) ; nu = dos
        if ((eb-e1).lt.tolr*eb) then
          nstall = nstall + 1
        else
          nstall = 0
        end if
        if (bro.eq.1) then
          call lm_broyden(j,x1-xb,y1-yb) ; brojac = .true.
        else
          newjac = .true.
        end if
        xb = x1 ; yb = y1 ; eb = e1

      elseif (ec.lt.eb) then

        ! bad step, but a point evaluated for the jacobian is better: move there
        if ((eb-ec).lt.tolr*eb) then
          nstall = nstall + 1
        else
          nstall = 0
        end if
        if (bro.eq.1) then
          call lm_broyden(j,xc-xb,yc-yb) ; brojac = .true.
        else
          newjac = .true.
        end if
        xb = xc ; yb = yc ; eb = ec ; ec = huge(uno)

      else

        ! bad step: increase damping; if the jacobian was updated with broyden,
        ! compute a new numerical jacobian
        lam = lam*nu ; nu = dos*nu
        nstall = nstall + 1
        if (brojac) newjac = .true.

      end if

      if (lm_converged(eb,yb)) then
        ind = 0 ; exit
      end if
      if (nstall.ge.mstall .or. lam.gt.1.0d16) then
        ind = 3 ; exit
      end if

    end do

    7 continue

    ! return the best point evaluated
    if (es.lt.eb) then
      xb = xs ; yb = ys ; eb = es
    end if
    y = yb
    x = xb

    if (ip.ge.1) then
      write(*,*) '  '
      if (ind.eq.0) write(*,*) ' lmmin finished: system has converged, func=0 '
      if (ind.eq.1) write(*,*) ' lmmin finished: gradient close to 0'
      if (ind.eq.2) write(*,*) ' lmmin finished: step in x close to 0'
      if (ind.eq.3) write(*,*) ' lmmin finished: no progress'
      if (ind.eq.9) write(*,*) ' lmmin finished: max number of evaluations'
      write(*,*          ) '   '
      write(*,*          ) ' evaluations = ', iy
      write(*,'(a,es12.4)') ' error       = ', eb
      write(*,'(a,f10.4)' ) ' reduction   = ', eb/max(e0,tiny(uno)) - uno
      write(*,*          ) '  '
    end if

    deallocate(j,a,g,dg,shck,xb,yb,x1,y1,dx,xc,yc,xs,ys,xj,yj)

    return
    97 format ('  new point    ',i5,'  ',3(es12.4))
    91 format ('  first point  ','     ','  ',es12.4)

    contains

    ! convergence: sum of squared errors or each equation below tolerance
    logical function lm_converged(e,yv)
      real(dp) , intent(in) :: e,yv(:)
      lm_converged = (e.lt.toler .or. maxval(abs(yv)).lt.tolereach)
    end function lm_converged

    ! keep the best point evaluated overall
    subroutine lm_keepbest(xv,yv,e)
      real(dp) , intent(in) :: xv(:),yv(:),e
      if (e.lt.es) then
        xs = xv ; ys = yv ; es = e
      end if
    end subroutine lm_keepbest

    ! numerical jacobian at xb (forward differences). sets done = .true. (and ind)
    ! if the system converges or the max number of evaluations is reached.
    subroutine lm_jacobian(done)
      logical , intent(out) :: done
      real(dp) :: sh
      integer  :: k
      done = .false. ; ec = huge(uno)
      if (ip.gt.1) write (*,*) '  '
      do ij = 1,n

        ! shock the ij-th unknown
        xj = xb
        sh = min(shck(ij),shck(ij)*max(0.30_dp,eb))
        if (xb(ij).lt.-toler) then
          xj(ij) = xb(ij) - sh*abs(xb(ij))
        elseif (xb(ij).gt.toler) then
          xj(ij) = xb(ij) + sh*abs(xb(ij))
        else
          xj(ij) = sh
        end if
        if (abs(xj(ij)-xb(ij)).lt.10.0_dp*toler) xj(ij) = xb(ij) + 10.0_dp*toler
        call lm_evalj(done) ; if (done) return

        ! if change in the error is too large, reduce the shock (at most 10 times)
        if (ej.gt.50.0_dp*eb) then
          k = 0
          do while (ej.gt.50.0_dp*eb .and. abs(xj(ij)-xb(ij)).gt.toler .and. k.lt.10)
            xj(ij) = (xb(ij) + xj(ij))/dos ; k = k + 1
            call lm_evalj(done) ; if (done) return
          end do
          call lm_newshock( )
        end if

        ! if the residuals do not change, enlarge the shock (at most 5 times)
        if (maxval(abs(yj-yb)).lt.toler) then
          k = 0
          do while (maxval(abs(yj-yb)).lt.toler .and. k.lt.5)
            xj(ij) = xb(ij) + dos*(xj(ij)-xb(ij)) ; k = k + 1
            call lm_evalj(done) ; if (done) return
          end do
          if (maxval(abs(yj-yb)).ge.toler) call lm_newshock( )
        end if

        ! column of the jacobian (zero if the unknown does not affect the residuals)
        if (maxval(abs(yj-yb)).lt.toler) then
          j(:,ij) = cero
          if (ip.gt.1) write (*,'(a,i4,a)') '   - jacobian, param =',ij,': no effect on residuals'
        else
          j(:,ij) = (yj-yb)/(xj(ij)-xb(ij))
        end if

      end do
      if (ip.gt.1) write (*,*) '  '
    end subroutine lm_jacobian

    ! evaluate a point for the jacobian
    subroutine lm_evalj(done)
      logical , intent(out) :: done
      done = .false.
      yj = func(xj) ; ej = sum(yj*yj) ; iy = iy + 1
      call lm_keepbest(xj,yj,ej)
      if (ej.lt.ec) then
        xc = xj ; yc = yj ; ec = ej
      end if
      if (ip.gt.1) write (*,'(a,i4,es12.4)') '   - jacobian, param =',ij,ej
      if (lm_converged(ej,yj)) then
        xb = xj ; yb = yj ; eb = ej ; ind = 0 ; done = .true. ; return
      end if
      if (iy.ge.maxiter) then
        ind = 9 ; done = .true.
      end if
    end subroutine lm_evalj

    ! update the relative shock of the ij-th unknown after resizing it
    subroutine lm_newshock( )
      if (abs(xb(ij)).gt.toler) shck(ij) = abs(xb(ij)-xj(ij))/abs(xb(ij))
      if (abs(xb(ij)).le.toler) shck(ij) = abs(xj(ij))
    end subroutine lm_newshock

    ! solve (am + lm·diag(d))·dxv = -gv by cholesky decomposition
    subroutine lm_cholsolve(am,d,lm,gv,dxv,ok)
      real(dp) , intent(in)  :: am(:,:),d(:),lm,gv(:)
      real(dp) , intent(out) :: dxv(:)
      logical  , intent(out) :: ok
      real(dp) :: b(size(gv),size(gv)),l(size(gv),size(gv)),z(size(gv)),s
      integer  :: i,k,nn
      nn = size(gv) ; b = am ; l = cero ; ok = .false. ; dxv = cero
      do i=1,nn
        b(i,i) = am(i,i) + lm*d(i)
      end do
      do k=1,nn
        s = b(k,k) - sum(l(k,1:k-1)**2)
        if (s.le.cero .or. s.ne.s) return
        l(k,k) = sqrt(s)
        do i=k+1,nn
          l(i,k) = (b(i,k) - sum(l(i,1:k-1)*l(k,1:k-1)))/l(k,k)
        end do
      end do
      do i=1,nn
        z(i) = (-gv(i) - sum(l(i,1:i-1)*z(1:i-1)))/l(i,i)
      end do
      do i=nn,1,-1
        dxv(i) = (z(i) - sum(l(i+1:nn,i)*dxv(i+1:nn)))/l(i,i)
      end do
      ok = all(dxv.eq.dxv)
    end subroutine lm_cholsolve

    ! broyden's update of the jacobian given a step dxv and the change in residuals dyv
    subroutine lm_broyden(jac,dxv,dyv)
      real(dp) , intent(inout) :: jac(:,:)
      real(dp) , intent(in)    :: dxv(:),dyv(:)
      real(dp) :: den
      integer  :: k
      den = sum(dxv*dxv)
      if (den.le.tiny(uno)) return
      do k=1,size(dyv)
        jac(k,:) = jac(k,:) + dxv(:)*( dyv(k) - sum(jac(k,:)*dxv(:)) )/den
      end do
    end subroutine lm_broyden

  end subroutine lmmin

  ! ----------------------------------------------------------------------------
  ! this subroutine applies the transformation
  !   y = log( r/(1-r) ) , with r = (x-xmin)/(xmax-xmin)
  ! to return a unbounded variable "y" from a bounded variable "x".
  ! to avoid infinite values at the bounds, r is restricted (symmetrically) to
  ! [1.0d-4, 1-1.0d-4], so that |y| <= 9.21. inside this interval, normalize is
  ! the exact inverse of denormalize.

  subroutine normalize(y,x,xmax,xmin)
    implicit none
    real(dp) , intent(in)  :: xmax,xmin,x
    real(dp) , intent(out) :: y
    real(dp) , parameter   :: rcap = 1.0d-4
    real(dp)               :: r
    y = cero
    if (xmax.le.xmin) then
      call error(' error in normalize: xmax <= xmin')
      return
    end if
    if (x.gt.xmax) then
      call error(' error in normalize: x > xmax')
      return
    end if
    if (x.lt.xmin) then
      call error(' error in normalize: x < xmin')
      return
    end if
    r = min(max((x-xmin)/(xmax-xmin),rcap),uno-rcap)
    y = log(r/(uno-r))
    return
  end subroutine normalize

  ! ----------------------------------------------------------------------------
  ! this subroutine applies the transformation
  !   x = xmin + ( exp(y)/(1+exp(y) )*(xmax-xmin)
  ! to return a bounded variable "x" from an unbounded variable "y".
  ! the logistic function is evaluated in a way that does not overflow.

  subroutine denormalize(y,x,xmax,xmin)
    implicit none
    real(dp) , intent(in)  :: xmax,xmin,y
    real(dp) , intent(out) :: x
    x = xmin
    if (xmax.lt.xmin) then
      call error(' error in denormalize: xmax < xmin')
      return
    end if
    if (y.ge.cero) then
      x = xmin + (xmax-xmin)/(uno+exp(-y))
    else
      x = xmin + (xmax-xmin)*exp(y)/(uno+exp(y))
    end if
    return
  end subroutine denormalize

  ! ----------------------------------------------------------------------------
  ! this subroutines updates a jacobian matrix using the broyden's method.
  ! if the step x1-x0 is (close to) zero, the jacobian is not updated.

  subroutine broyden(jac,x1,x0,f1,f0)
    implicit none
    real(dp) , intent(in)    :: x1(:),f1(:),x0(:),f0(:)
    real(dp) , intent(inout) :: jac(:,:)
    real(dp)                 :: df(size(f0,1))
    real(dp)                 :: dx(size(x1,1))
    real(dp)                 :: j0(size(f0,1),size(x1,1))
    real(dp)                 :: den
    integer                  :: i
    ! check dimensions before any operation
    if (size(f1).ne.size(jac,1) .or. size(f0).ne.size(jac,1)) then
      call error(' error in broyden: size(f) /= size(j,1)')
      return
    end if
    if (size(x1).ne.size(jac,2) .or. size(x0).ne.size(jac,2)) then
      call error(' error in broyden: size(x) /= size(j,2)')
      return
    end if
    j0    = jac
    dx(:) = x1(:)-x0(:)
    df(:) = f1(:)-f0(:)
    den   = sum(dx(:)*dx(:))
    ! step too small (relative to x0, with an absolute floor): do not update
    if (maxval(abs(dx)).le.tolvl*max(uno,maxval(abs(x0))) .or. den.le.tiny(uno)) then
      call error(' error in broyden: step x1-x0 close to 0, jacobian not updated')
      return
    end if
    do i=1,size(df)
      jac(i,:) = j0(i,:) + dx(:)*( df(i)-sum(j0(i,:)*dx(:)) ) / den
    end do
    return
  end subroutine broyden

end module toolkit_optim


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module toolkit
  use toolkit_pars
  use toolkit_tools
  use toolkit_econ
  use toolkit_linalg
  use toolkit_stats
  use toolkit_optim
  implicit none
end module toolkit

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%