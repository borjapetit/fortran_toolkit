
<span style="text-align:right;display:block;">
<a href="https://borjapetit.github.io/fortran_toolkit/">Back to index</a>
</span>

## ```fisheryates```

```fortran
subroutine fisheryates(vect,seed)
  implicit none
  integer , intent(inout)         :: vect(:)
  integer , intent(in) , optional :: seed(:)
```

This subroutine takes an array ```vect``` of dimension ```n``` and fills it with integers from 1 to ```n``` and them shuffled using the the Fisher–Yates shuffle algorithm. The user can fixed the seed by providing ```seed```.

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


