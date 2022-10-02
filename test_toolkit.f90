
program test
use toolkit
implicit none

call test_ols()

return

contains

subroutine test_ols

    implicit none
    integer  , parameter :: nn = 2000
    real(dp) :: y(nn),x1(nn),x2(nn),x3(nn),x0(nn),b(4)
    
    call randomnormal(x0,dble(0.10),dble(0.30))
    call randomnormal(x1,dble(0.00),dble(0.40))
    call randomnormal(x2,dble(0.00),dble(0.50))
    call randomnormal(x3,dble(0.00),dble(0.60))
    
    y = x0 + dble(0.7)*x1 - dble(0.5)*x2 + dble(0.2)*x3
    
    write(*,98) ' '
    write(*,98) ' **********************************************************************'
    write(*,98) ' '
    write(*,98) ' TESTING OLS REGRESSION '
    write(*,98) ' '
    write(*,98) ' True values:' , dble(0.10) , dble(0.7) , -dble(0.5) , dble(0.2)
    write(*,98) ' '
    write(*,97) ' Num. Obs:' , nn
    write(*,98) ' '
    write(*,98) ' With a constant ' ; call olsreg(b,y,x1,x2,x3,iprint=1)
    write(*,98) ' With no constant' ; call olsreg(b(1:3),y,x1,x2,x3,iprint=1)
    write(*,98) ' '

    return
    98 format (a,4(f7.2))
    97 format (a,i7)
end subroutine test_ols

end program test