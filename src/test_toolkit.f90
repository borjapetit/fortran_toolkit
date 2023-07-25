
program test
use toolkit
implicit none
real(dp) :: mat0(2,2),vec(10)

mat0(:,1) = (/1,2/)
mat0(:,2) = (/3,4/)

vec = grid(10.0d0,1.0d0,10)

print * , diag(mat0)
print * , cumsum(vect(mat0))
print * , vec

print * , ' '

call test_interpolate()

call test_ols()

print * , ' '

return
contains

subroutine test_ols
  implicit none
  integer  , parameter :: nn = 100
  real(dp) :: y(nn),x1(nn),x2(nn),x3(nn),x0(nn),b(4)
  integer :: i
  
  call randomnormal(x0,dble(0.10),dble(0.30))
  call randomnormal(x1,dble(0.00),dble(0.40))
  call randomnormal(x2,dble(0.00),dble(0.50))
  call randomnormal(x3,dble(0.00),dble(0.60))
  
  y = x0 + dble(0.7)*x1 - dble(0.5)*x2 + dble(0.2)*x3

  open (unit=1, file="data.txt",action="write")
    do i=1,nn
      write(1,*) y(i),x0(i),x1(i),x2(i),x3(i)
    end do
  close (1)
  
  write(*,98) ' **********************************************************************'
  write(*,98) ' '
  write(*,98) ' TESTING STATISTICS '
  write(*,98) ' '
  write(*,98) ' Mean x0     =' , varmean(x0)
  write(*,98) ' Std of x1   =' , varstd(x1)
  write(*,98) ' Corr (y,x2) =' , correlation(y,x2)
  write(*,98) ' p50 of x3   =' , percentile(x3,0.50d0)
  write(*,98) ' '
  write(*,98) ' OLS REGRESSION '
  write(*,98) ' '
  write(*,98) ' True values:' , dble(0.10) , dble(0.7) , -dble(0.5) , dble(0.2)
  write(*,98) ' '
  write(*,97) ' Num. Obs:' , nn
  write(*,98) ' '
  write(*,98) ' With a constant ' ; call olsreg(b,y,x1,x2,x3,table=1)
  write(*,98) ' With no constant' ; call olsreg(b(1:3),y,x1,x2,x3,table=1)
  write(*,98) ' '
  98 format (a,4(f7.2))
  97 format (a,i7)
  return
end subroutine test_ols

subroutine test_interpolate
  implicit none
  real(dp) :: x1(10),x2(10),x3(10),x4(10),y(10,10,10,10)
  real(dp) :: a1,a2,a3,a4
  integer  :: i1,i2,i3,i4

  x1 = grid(1.0d0,0.0d0,10)
  x2 = grid(1.0d0,0.0d0,10)
  x3 = grid(1.0d0,0.0d0,10)
  x4 = grid(1.0d0,0.0d0,10)

  do i1=1,10
    do i2=1,10
      do i3=1,10
        do i4=1,10
          y(i1,i2,i3,i4)=func0((/x1(i1),x2(i2),x3(i3),x4(i4)/))
        enddo
      enddo
    enddo
  enddo

  call random_number(a1)
  call random_number(a2)
  call random_number(a3)
  call random_number(a4)
  
  write(*,98) ' **********************************************************************'
  write(*,98) ' '
  write(*,98) ' TESTING INTERPOLATE '
  write(*,98) ' '
  write(*,98) ' True value   =' , func0((/a1,a2,a3,a4/))
  write(*,98) ' Interpolated =' , interpolate(a1,a2,a3,a4,x1,x2,x3,x4,y)
  write(*,98) ' '

  98 format (a,4(f10.6))
  return
end subroutine test_interpolate

function func0(a) result(fun)
  implicit none
  real(dp) :: a(4),fun
  fun = (a(1)**(2)) - dble(2)*(a(2)**(2)) + dble(3)*(a(3)**(2)) - dble(4)*(a(4)**(2))
  return
end function func0

end program test
