program example
  use toolkit , only : dp,brent
  implicit none
  real(dp) :: xvar
  integer  :: numiter
  integer  :: exitcode
    
  write(*,*) 'find the root: initial points 0 and 15, with f(0)<0 and f(15)>0'
  
  call brent(func,xvar,numiter,exitcode,0.0d0,15.0d0,iprint=2)

  write(*,*) 'find the root: initial points 5 and 10, with f(5)>0 and f(10)>0'

  call brent(func,xvar,numiter,exitcode,5.0d0,10.0d0,iprint=2)

  return

  contains

  function func(x) result(y)
    implicit none
    real(dp) :: x,y
    y = x**2.0 - 3.5d0*x - 5.0d0
    return
  end function func

end program example