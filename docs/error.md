### error

```fortran
subroutine error(mess,i)
  implicit none
  character(len=*) , intent(in) :: mess
  integer , optional , intent(in) :: i
```

This subroutine prints an error message ```mess```. If the user supplies any integer ```i``` the subroutine pauses the execution of the program until the user types any key.

**Dependencies**: none

[(back to index)](../index.md)

---

**Example**

```fortran
! error message with execution pause
call error(' error message',1)

! error message without execution pause
call error(' error message')
```


