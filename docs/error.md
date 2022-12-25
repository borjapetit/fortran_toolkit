# error

```fortran
subroutine error(mess)
  implicit none
  character(len=*) , intent(in) :: mess
  integer , optional , intent(inout) :: mess
```

This subroutine prints an error message ```mess```. If the user supplies an integer ```i``` the subroutine interrupts the execution of the program until the user type an interger.

**Dependencies**: none

[(back to index)](../index.md)
