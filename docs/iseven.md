## iseven

```fortran
function iseven(num) result(ise)
  implicit none
  integer :: num
  logical :: ise
```



This function checks whether a number ```num``` is even. The result ```ise``` is a logical variable taking value ```.TRUE.``` if ```num``` is even, and ```.FALSE.``` otherwise.

_Example_: check whether 3 and 6 are eve:

```fortran
write(*,*) iseven(3)  ! .FALSE.
write(*,*) iseven(6)  ! .TRUE.
```

**Dependencies**: 
 [```multiplo```](multiplo.md)

[(back to index)](index.md)