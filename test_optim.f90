
! compilation command:
! gfortran toolkit.f90 test.f90 -o test
! gfortran -fcheck=all -o0 -fbacktrace toolkit.f90 test.f90 -o test

program test
use toolkit
implicit none
integer , parameter :: numsol = 1000
real(dp)  :: x0(2,numsol),xx
real(dp)  :: xs(2),x1(2)
real(dp)  :: ys(2),y1(2)
real(dp)  :: aveer,aveiter,aveconv
integer   :: iter1,j1,i

do i = 1,numsol
    call random_number(xx) ; x0(1,i) = dble(10.00)*(xx-dble(0.5))
    call random_number(xx) ; x0(2,i) = dble(10.00)*(xx-dble(0.5))
end do

write(*,99) '   '
write(*,99) '   '
write(*,99) ' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
write(*,99) '   '
write(*,99) ' testing minimization algorithm '
write(*,99) '   '
write(*,99) ' I take the functions from Wikipedia. Check:'
write(*,99) '   '
write(*,99) '    https://en.wikipedia.org/wiki/test_functions_for_optimization '
write(*,99) '   '
write(*,99) ' for more details'
write(*,99) '   '
write(*,99) ' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
write(*,99) '   '

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

write(*,99) '   '
write(*,99) ' ---------------------------------------------------------------------'
write(*,99) '   '
write(*,97) ' matyas (solution: x1 = 0.0, x2 = 0.0) '
write(*,99) '   '
write(*,97) '  solution : x1 =  0.0, x2 =  0.0' , matyas_1d( (/  0.0d0 , 0.0d0 /))
write(*,99) '   '

aveer   = zero
aveiter = zero
aveconv = zero
do i = 1,numsol
  call lmmin(matyas,x1,y1,iter1,j1,x0(:,i)*dble(2.0))
  aveer   = max(aveer,sum(y1(:)*y1(:)))
  aveiter = aveiter + dble(iter1)
  aveconv = aveconv + dble(min(0,j1))
end do
aveiter = aveiter/dble(numsol)
aveconv = - aveconv/dble(numsol) + dble(1.00)

write(*,98) ' matyas (jacobian):', aveer,aveiter

aveer   = zero
aveiter = zero
aveconv = zero
do i = 1,numsol
  call lmmin(matyas,x1,y1,iter1,j1,x0(:,i)*dble(2.0),usebro=1)
  aveer   = max(aveer,sum(y1(:)*y1(:)))
  aveiter = aveiter + dble(iter1)
  aveconv = aveconv + dble(min(0,j1))
end do
aveiter = aveiter/dble(numsol)
aveconv = - aveconv/dble(numsol) + dble(1.00)

write(*,98) ' matyas (broyden): ', aveer,aveiter

aveer   = zero
aveiter = zero
aveconv = zero
do i = 1,numsol
  call simplex(matyas_1d,x1,y1(1),iter1,j1,x0(:,i)*dble(2.0))
  aveer   = max(aveer,y1(1))
  aveiter = aveiter + dble(iter1)
  aveconv = aveconv + dble(min(0,j1))
end do
aveiter = aveiter/dble(numsol)
aveconv = - aveconv/dble(numsol) + dble(1.00)

write(*,98) ' matyas (simplex): ', aveer,aveiter

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

write(*,99) '   '
write(*,99) ' ---------------------------------------------------------------------'
write(*,99) '   '
write(*,97) ' booth'
write(*,99) '   '
write(*,97) '  solution : x1 =  1.0, x2 =  3.0' , booth_1d( (/  1.0d0 , 3.0d0 /))
write(*,99) '   '

aveer   = zero
aveiter = zero
aveconv = zero
do i = 1,numsol
  call lmmin(booth,x1,y1,iter1,j1,x0(:,i)*dble(2.0))
  aveer   = max(aveer,sum(y1(:)*y1(:)))
  aveiter = aveiter + dble(iter1)
  aveconv = aveconv + dble(min(0,j1))
end do
aveiter = aveiter/dble(numsol)
aveconv = - aveconv/dble(numsol) + dble(1.00)

write(*,98) ' booth (jacobian):', aveer,aveiter!,aveconv

aveer   = zero
aveiter = zero
aveconv = zero
do i = 1,numsol
  call lmmin(booth,x1,y1,iter1,j1,x0(:,i)*dble(2.0),usebro=1)
  aveer   = max(aveer,sum(y1(:)*y1(:)))
  aveiter = aveiter + dble(iter1)
  aveconv = aveconv + dble(min(0,j1))
end do
aveiter = aveiter/dble(numsol)
aveconv = - aveconv/dble(numsol) + dble(1.00)

write(*,98) ' booth (broyden): ', aveer,aveiter

aveer   = zero
aveiter = zero
aveconv = zero
do i = 1,numsol
  call simplex(booth_1d,x1,y1(1),iter1,j1,x0(:,i)*dble(2.0))
  aveer   = max(aveer,y1(1))
  aveiter = aveiter + dble(iter1)
  aveconv = aveconv + dble(min(0,j1))
end do
aveiter = aveiter/dble(numsol)
aveconv = - aveconv/dble(numsol) + dble(1.00)

write(*,98) ' booth (simplex): ', aveer,aveiter

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

write(*,99) '   '
write(*,99) ' ---------------------------------------------------------------------'
write(*,99) '   '
write(*,97) ' hummelblau'
write(*,99) '   '
write(*,97) '  solution 1: x1 =  3.000000, x2 =  2.000000' , hummelblau_1d( (/  3.0000d0 , 2.0000d0 /))
write(*,97) '  solution 2: x1 = -2.805118, x2 =  3.131312' , hummelblau_1d( (/ -2.805118d0,  3.131312d0 /))
write(*,97) '  solution 3: x1 = -3.779310, x2 = -3.283186' , hummelblau_1d( (/ -3.779310d0, -3.283186d0 /))
write(*,97) '  solution 3: x1 =  3.584428, x2 = -1.848126' , hummelblau_1d( (/  3.584428d0, -1.848126d0 /))
write(*,99) '   '

aveer   = zero
aveiter = zero
aveconv = zero
do i = 1,numsol
  call lmmin(hummelblau,x1,y1,iter1,j1,x0(:,i))
  aveer   = max(aveer,sum(y1(:)*y1(:)))
  aveiter = aveiter + dble(iter1)
  aveconv = aveconv + dble(min(0,j1))
end do
aveiter = aveiter/dble(numsol)
aveconv = - aveconv/dble(numsol) + dble(1.00)

write(*,98) ' lmmin (jacobian):', aveer,aveiter!,aveconv

aveer   = zero
aveiter = zero
aveconv = zero
do i = 1,numsol
  call lmmin(hummelblau,x1,y1,iter1,j1,x0(:,i),usebro=1)
  aveer   = max(aveer,sum(y1(:)*y1(:)))
  aveiter = aveiter + dble(iter1)
  aveconv = aveconv + dble(min(0,j1))
end do
aveiter = aveiter/dble(numsol)
aveconv = - aveconv/dble(numsol) + dble(1.00)

write(*,98) ' lmmin (broyden): ', aveer,aveiter!,aveconv

aveer   = zero
aveiter = zero
aveconv = zero
do i = 1,numsol
  call simplex(hummelblau_1d,x1,y1(1),iter1,j1,x0(:,i))
  aveer   = max(aveer,y1(1))
  aveiter = aveiter + dble(iter1)
  aveconv = aveconv + dble(min(0,j1))
end do
aveiter = aveiter/dble(numsol)
aveconv = - aveconv/dble(numsol) + dble(1.00)

write(*,98) ' lmmin (simplex): ', aveer,aveiter
write(*,99) '   '
write(*,99) ' ---------------------------------------------------------------------'
write(*,99) '   '
write(*,99) '   '
write(*,99) '   '
write(*,99) '   '

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

99 format (a,f10.4,f10.4,i10)
98 format (a,'  max error = ',f8.4,' ;  average iter = ',f8.4)
97 format (a,f10.4,f10.4,' error = ',f6.4)

contains

  function hummelblau(x) result(y)
    implicit none
    real(dp) :: x(:)
    real(dp) , allocatable :: y(:)
    allocate(y(2))
    y(1) = (x(1)*x(1) + x(2) - dble(11.00))**2.0d0
    y(2) = (x(1) + x(2)*x(2) - dble(7.00))**2.0d0
    return
  end function hummelblau
  function hummelblau_1d(x) result(y)
    implicit none
    real(dp) :: x(:),y,y0(2)
    y0 = hummelblau(x)
    y  = sum(y0)
    return
  end function hummelblau_1d

  function booth(x) result(y)
    implicit none
    real(dp) :: x(:)
    real(dp) , allocatable :: y(:)
    allocate(y(2))
    y(1) = x(1) + dble(2.0)*x(2) - dble(7.0)
    y(2) = dble(2.0)*x(1) + x(2) - dble(5.0)
    return
  end function booth
  function booth_1d(x) result(y)
    implicit none
    real(dp) :: x(:),y,y0(2)
    y0 = booth(x)
    y  = sum(y0)
    return
  end function booth_1d

  function matyas(x) result(y)
    implicit none
    real(dp) :: x(:)
    real(dp) , allocatable :: y(:)
    allocate(y(2))
    y(1) = dble(0.26)*( x(1)*x(1) + x(2)*x(2))
    y(2) = dble(0.48)*x(1)*x(2)
    return
  end function matyas
  function matyas_1d(x) result(y)
    implicit none
    real(dp) :: x(:),y
    y = dble(0.26)*( x(1)*x(1) + x(2)*x(2)) - dble(0.48)*x(1)*x(2)
    return
  end function matyas_1d

end program test
    