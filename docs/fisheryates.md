
<span style="text-align:right;display:block;">
<a href="https://borjapetit.github.io/fortran_toolkit/">Back to index</a>
</span>

## ```fisheryates```

```fortran
subroutine fisheryates(nums)
  implicit none
  integer , intent(inout) :: nums(:)
```

This subroutine takes an array ```nums``` of dimension ```n```, fills it with the integers from 1 to ```n```, and then shuffles them using the Fisher–Yates shuffle algorithm. The initial values of ```nums``` are overwritten. The subroutine uses the intrinsic ```random_number```, so the user can fix the seed with the intrinsic ```random_seed``` before calling it.

**Internal dependencies**: none

---

**Example**:

```fortran
integer :: vect(4)

call fisheryates(vect)
print * , 'vect =', vect   ! vect =  3  1  4  2

call fisheryates(vect)
print * , 'vect =', vect   ! vect =  2  4  1  3
```


