program main

    use toolkit , only : dp,normalize,denormalize,simplex
    implicit none
    
    real(dp) :: x0(2),z0(2),f0
    real(dp) :: x1(2),z1(2),f1
    integer  :: numiter,exitcode
  
    ! Initial guess:
    x0(1) = 5.0d0
    x0(2) = -3.0d0

    ! x(1) and x(2) should be between 10 and -10  
    ! transform the (constrained) initial point "x0" into the 
    ! (unconstrained) initial point "z0"
    call normalize(z0(1),x0(1),10.0d0,-10.0d0)
    call normalize(z0(2),x0(2),10.0d0,-10.0d0)
  
    ! use Nelder-Mead to minimize the function "unconstrained_booth" 
    ! over "z", using "z0" as initial point
    call simplex(unconstrained_booth,z1,f1,numiter,exitcode,z0,iprint=2)
      
    ! transform "z1" (unconstrained) into "x1" (constrained)
    call denormalize(z1(1),x1(1),10.0d0,-10.0d0)
    call denormalize(z1(2),x1(2),10.0d0,-10.0d0)
    
    ! these variable satisfy: unconstrained_booth(z1) = booth(x1)
  
    write(*,'(a,3(f10.4))') ' Actual solution:        ',1.0d0,3.0d0, booth((/1.0d0,3.0d0/))
    write(*,'(a,3(f10.4))') ' Unconstrained solution: ',z1(:), unconstrained_booth(z1)
    write(*,'(a,3(f10.4))') ' Constrained solution:   ',x1(:), booth(x1)
    
    return
    contains
    
    ! the function (of "x") that we want to minimize
    function booth(x) result(f)
      implicit none
      real(dp) :: x(:),f
      f =  ( x(1) + 2.0d0*x(2) - 7.0d0 )**2.0d0 + &
           ( 2.0d0*x(1) + x(2) - 5.0d0 )**2.0d0
      return
    end function booth

    ! auxiliary function that takes an unconstrained vector "y"
    ! transforms it into the constrained vector "x" and get the 
    ! value of the function "func"
    function unconstrained_booth(z) result(f)
      implicit none
      real(dp) :: f
      real(dp) :: z(:)
      real(dp) :: x(2)
      
      ! get values of "x"
      call denormalize(z(1),x(1),10.0d0,-10.0d0)
      call denormalize(z(2),x(2),10.0d0,-10.0d0)
      
      ! evaluate the function
      f = booth(x)
      
      return
    end function unconstrained_booth
    
  end program main