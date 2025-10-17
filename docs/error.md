
<span style="text-align:right;display:block;">
<a href="https://borjapetit.github.io/fortran_toolkit/">Back to index</a>
</span>

## ```error```

```fortran
subroutine error(mess,i)
  implicit none
  character(len=*)   , intent(in) :: mess
  integer , optional , intent(in) :: i
```

This subroutine prints an error message `mess`. If the user supplies any integer `i` the subroutine pauses the execution of the program until the user types any key.

**Internal dependencies**: none

---

**Example**

Error message with execution pause:
```fortran
call error('error message',1)

! The output would be:
!   error message --> type any key to continue...
```

Error message without execution pause:
```fortran
call error('error message')

! The output would be:
!   error message
```









