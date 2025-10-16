
## ```num2text```

```fortran
function num2text(realnum,dec) result(txt)
  implicit none
  real(kind=8)       :: realnum
  integer , optional :: decs
  character(len=:)   :: txt
```
```fortran
function num2text(integer) result(txt)
  implicit none
  integer          :: integer
  character(len=:) :: txt
```

This function returns a string ```txt``` with the user-provided number. The user can provide either a real number (double precision) or an integer. In the case of a real number, the user can also provide the number of decimals ```dec``` (default is 2).

**Dependencies**: none

**Note**: ```num2text``` is an interface that calls specific functions depending on the inputs provided . The specific functions are:
- ```real2text``` if the user provide a real number
- ```int2text``` if the user provide an integer

[(back to index)](../index.md)

---

**Example**:

```fortran
  ! print a real nummber 
  num = 1234.56789d0
  txt = num2text(num) ! default decimals = 2
  print * , 'num =', num, ' || The number is'//trim(txt)

  ! output:
  !     rnum =  1234.5678900000000 || The number is 1234.57
  
  ! print a real nummber setting the number of decimals
  num = 3.14159265359
  txt = num2text(rnum,5)
  print * , 'num =', num, ' || The number is'//trim(txt)

  ! output:
  !     rnum =  3.1415926535900 || The number is 3.14159

  ! print an integer
  int = 56789
  txt = num2text(int)
  print * , 'int =', int, ' || The number is'//trim(txt)
  
  ! output:
  !     rnum =  56789 || The number is 56789
```



