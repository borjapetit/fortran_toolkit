## error

```fortran
subroutine error(mess,i)
  implicit none
  character(len=*)   , intent(in) :: mess
  integer , optional , intent(in) :: i
```

This subroutine prints an error message $\texttt{mess}$. If the user supplies any integer $\texttt{i}$ the subroutine pauses the execution of the program until the user types any key.

**Dependencies**: none

[(back to index)](../index.md)

---

**Example**

Error message with execution pause:
```fortran
call error('error message',1)
```
The output would be:
```
error message
```

Error message without execution pause:
```fortran
call error('error message')
```
The output would be:
```
error message
Type any key to continue...
```









