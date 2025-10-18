! this program tests different subroutines in the toolkit package.
! to compile the code, use the following command:
!   gfortran  toolkit.f90 example.f90 -o example

program toolkit_examples

  use toolkit
  !use, intrinsic :: ieee_arithmetic
  
  implicit none

  !real(dp) , parameter :: nadp = ieee_value(0.0_dp, ieee_quiet_nan)
  
  call system('clear')
 
  write(*,'(a)') ' '
  write(*,'(a)') ' ' 
  write(*,'(a)') ' '

  write(*,'(a)') ' *************************************************************'
  write(*,'(a)') ' TESTING THE LMMIN ALGORITHM WITH AND WITHOUT BROYDEN UPDATING'
  write(*,'(a)') ' *************************************************************'
  write(*,'(a)') ' '

  ! this subroutine tests the lmmin algorithm by minimizing a 2-dimension non-linear system.
  ! the system contains a dummy loop to make the computation slower, so that timing differences can be appreciated.

  call demo_lmmin()

  write(*,'(a)') ' '
  write(*,'(a)') ' '

  write(*,'(a)') ' *************************************************************'
  write(*,'(a)') ' TESTING THE SIMPLEX ALGORITHM USING THE MAYTAS FUNCTION '
  write(*,'(a)') ' *************************************************************'
  write(*,'(a)') ' '

  ! this subroutine tests the simples algorithm by minimizing athe maytas function.
  ! the function contains a dummy loop to make the computation slower, so that timing differences can be appreciated.

  call demo_simplex()

  write(*,'(a)') ' '
  write(*,'(a)') ' '

  contains

    subroutine demo_lmmin
      implicit none
      real(dp) , allocatable :: x(:),y(:)
      real(dp) :: x0(2)
      real(dp) :: residual_norm
      integer  :: iy, exitcode
      real(dp) :: t_start, t_end

      allocate(x(2),y(2))

      x0 = [1.0_dp, 4.0_dp]
      x  = cero
      y  = cero

      ! -----------------------------------------------------------------------------

      t_start = timing( )

      write(*,'(a)')         '  --- lmmin example w/o broyden updating -------------'
      write(*,'(a,2f12.6)')  '  Initial guess:            ', x0
      write(*,'(a,2f12.6)')  '  Solution:                 ', [6.0_dp,1.0_dp]

      call lmmin(nonlinear_system,x,y,iy,exitcode,x0)

      residual_norm = sqrt(sum(y**2))

      t_end = timing( )

      write(*,'(a,2f12.6)')  '  Solution found:           ', x
      write(*,'(a,2es12.4)') '  Residuals y:              ', y
      write(*,'(a,f12.6)')   '  Residual 2-norm:          ', residual_norm
      write(*,'(a,i6)')      '  Function evaluations (iy):', iy
      write(*,'(a,i6)')      '  Exit code (ind):          ', exitcode
      write(*,'(a,f12.6)')   '  Time elapsed (secs):      ', t_end - t_start
      write(*,'(a)')         '  ----------------------------------------------------'
      write(*,'(a)')         '                                                      '

      ! -----------------------------------------------------------------------------

      t_start = timing( )

      write(*,'(a)')         '  --- lmmin example w broyden updating ---------------'
      write(*,'(a,2f12.6)')  '  Initial guess:            ', x0
      write(*,'(a,2f12.6)')  '  Solution:                 ', [6.0_dp,1.0_dp]

      call lmmin(nonlinear_system,x,y,iy,exitcode,x0, iprint=0, usebro=1)

      residual_norm = sqrt(sum(y**2))

      t_end = timing( )

      write(*,'(a,2f12.6)')  '  Solution found:           ', x
      write(*,'(a,2es12.4)') '  Residuals y:              ', y
      write(*,'(a,f12.6)')   '  Residual 2-norm:          ', residual_norm
      write(*,'(a,i6)')      '  Function evaluations (iy):', iy
      write(*,'(a,i6)')      '  Exit code (ind):          ', exitcode
      write(*,'(a,f12.6)')   '  Time elapsed (secs):      ', t_end - t_start
      write(*,'(a)')         '  ----------------------------------------------------'

      deallocate(x,y)
      return
    end subroutine demo_lmmin

    function nonlinear_system(xvec) result(resid)
      real(dp) :: xvec(:)
      real(dp), allocatable :: resid(:)
      call dummy_loop()
      allocate(resid(2))
      resid(1) = xvec(1)**2 + xvec(2) - 37.0_dp
      resid(2) = xvec(1) - xvec(2)**2 - 5.0_dp
      return
    end function nonlinear_system

    subroutine demo_simplex( )
      implicit none
      real(dp) , allocatable :: x(:)
      real(dp) :: x0(2)
      real(dp) :: residual_norm
      integer  :: iy, exitcode
      real(dp) :: t_start, t_end

      allocate(x(2))

      x0 = [2.0_dp, 3.0_dp]
      x  = cero

      t_start = timing( )

      write(*,'(a)')         '  --- simplex example (Matyas function) --------------'
      write(*,'(a,2f12.6)')  '  Initial guess:            ', x0
      write(*,'(a,2f12.6)')  '  Solution:                 ', [cero,cero]

      call simplex(matyas,x,residual_norm,iy,exitcode,x0)

      t_end = timing( )

      write(*,'(a,2f12.6)')  '  Found solution:           ', x
      write(*,'(a,f12.6)')   '  Residual 2-norm:          ', residual_norm
      write(*,'(a,i6)')      '  Function evaluations (iy):', iy
      write(*,'(a,i6)')      '  Exit code (ind):          ', exitcode
      write(*,'(a,f12.6)')   '  Time elapsed (secs):      ', t_end - t_start
      write(*,'(a)')         '  ----------------------------------------------------'
      write(*,'(a)')         '                                                      '


      x0 = [1.0_dp, 4.0_dp]
      x  = cero

      t_start = timing( )

      write(*,'(a)')         '  --- simplex example (non liner system) -------------'
      write(*,'(a,2f12.6)')  '  Initial guess:            ', x0
      write(*,'(a,2f12.6)')  '  Solution:                 ', [6.0_dp,1.0_dp]

      call simplex(sum_nonlinear_system,x,residual_norm,iy,exitcode,x0)

      t_end = timing( )

      write(*,'(a,2f12.6)')  '  Found solution:           ', x
      write(*,'(a,f12.6)')   '  Residual 2-norm:          ', residual_norm
      write(*,'(a,i6)')      '  Function evaluations (iy):', iy
      write(*,'(a,i6)')      '  Exit code (ind):          ', exitcode
      write(*,'(a,f12.6)')   '  Time elapsed (secs):      ', t_end - t_start
      write(*,'(a)')         '  ----------------------------------------------------'

      deallocate(x)
      return
    end subroutine demo_simplex

    function matyas(x) result(y)
      implicit none
      real(dp) :: x(:),y
      call dummy_loop()
      y = dble(0.26)*( x(1)*x(1) + x(2)*x(2)) - dble(0.48)*x(1)*x(2)
      return
    end function matyas

    function sum_nonlinear_system(xvec) result(y)
      real(dp) :: xvec(:),ys(2),y
      call dummy_loop()
      ys(1) = xvec(1)**2 + xvec(2) - 37.0_dp
      ys(2) = xvec(1) - xvec(2)**2 - 5.0_dp
      y     = sum(ys**2)
      return
    end function sum_nonlinear_system

    subroutine dummy_loop()
      integer :: i
      do i=1,50000000
        ! do nothing
      end do
      return
    end subroutine dummy_loop

end program toolkit_examples
