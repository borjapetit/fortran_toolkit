
<span style="text-align:right;display:block;">
<a href="https://borjapetit.github.io/fortran_toolkit/">Back to index</a>
</span>

## ```error```

```fortran
subroutine error(mess,i)
  implicit none
  character(len=*)   , intent(in) :: mess
  integer , optional , intent(in) :: i     ! or: logical , optional , intent(in) :: i
```

This subroutine prints an error message `mess`. If the user supplies an integer `i` different from 0 (or the logical `.true.`), the subroutine pauses the execution of the program until the user presses Enter (any text typed before Enter is ignored).

Messages are not printed (and there is no pause) if the global variable `print_err` is set to `.false.`.

**Note**: `error` is an interface that calls `error_0` (only the message), `error_1` (message and integer) or `error_2` (message and logical). Since it writes to the standard output, a function that calls `error` should not be used inside a `print` or `write` statement (Fortran does not allow recursive I/O): assign its result to a variable first.

**Internal dependencies**: none

---

**Example**

Error message with execution pause:
```fortran
call error('error message',1)

! The output would be:
!   error message --> press enter to continue...
```

Error message without execution pause:
```fortran
call error('error message')

! The output would be:
!   error message
```









