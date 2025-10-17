
<span style="text-align:right;display:block;">
<a href="https://borjapetit.github.io/fortran_toolkit/">Back to index</a>
</span>

## ```iseven```

```fortran
function iseven(num) result(ise)
  implicit none
  integer :: num
  logical :: ise
```

This function checks whether a number ```num``` is even. The result ```ise``` is a logical variable taking value ```.TRUE.``` if ```num``` is even, and ```.FALSE.``` otherwise.

**Internal dependencies**: [```multiplo```](multiplo.md)

---

**Example**

Check whether 3 and 6 are even:

```fortran
write(*,*) iseven(3)  ! .FALSE.
write(*,*) iseven(6)  ! .TRUE.
```
