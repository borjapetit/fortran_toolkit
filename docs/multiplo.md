## multiplo

```fortran
elemental function multiplo(num,xx) result(mul)
  implicit none
  integer :: num,xx
  logical :: mul
```

_Dependencies_: ```none```

This function checks whether a number ```num0``` is a multiple of ```num1```. The result ```mul``` is a logical variable taking value ```.TRUE.``` if ```num0``` is a multiple of ```num1```, and ```.FALSE.``` otherwise.

_Example_: check whether 25 and 27 are multiples of 5:

```fortran
write(*,*) multiplo(5,25)  ! .TRUE.
write(*,*) multiplo(5,27)  ! .FALSE.
```

[(back to index)](inicio.md)