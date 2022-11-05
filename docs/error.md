## error

```fortran
subroutine error(mess)
  implicit none
  character(len=*) , intent(in) :: mess
```

This subroutine prints an error message ```mess``` and interrupt the execution of the program until the user type an interger.

**Dependencies**: none

[(back to index)](inicio.md)