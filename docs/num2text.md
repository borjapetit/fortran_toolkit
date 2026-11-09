
<span style="text-align:right;display:block;">
<a href="https://borjapetit.github.io/fortran_toolkit/">Back to index</a>
</span>

## ```num2text```

```fortran
function num2text(realnum,decs) result(txt)
  implicit none
  real(kind=8)                   :: realnum
  integer , optional             :: decs
  character(len=:) , allocatable :: txt
```
```fortran
function num2text(inter,width,fill) result(txt)
  implicit none
  integer                        :: inter
  integer          , optional    :: width
  character(len=1) , optional    :: fill
  character(len=:) , allocatable :: txt
```

This function returns a string ```txt``` with the user-provided number. The user can provide either a real number (double precision) or an integer. In the case of a real number, the user can also provide the number of decimals ```decs``` (default is 2). In the case of an integer, the user can also provide a minimum width ```width```: if the number has fewer digits than ```width```, it is padded on the left with the character ```fill``` (default is ```'0'```). When padding a negative number with ```'0'```, the minus sign goes first (e.g., ```num2text(-5,4)``` returns ```-005```); with any other ```fill``` character the sign stays next to the digits (```num2text(-5,4,' ')``` returns ```'  -5'```).

The output string has no leading or trailing blanks.

**Internal dependencies**: none

**Note**: ```num2text``` is an interface that calls specific functions depending on the inputs provided . The specific functions are:
- ```real2text``` if the user provide a real number
- ```int2text``` if the user provide an integer

---

**Example**:

```fortran
  ! print a real nummber 
  num = 1234.56789d0
  txt = num2text(num) ! default decimals = 2
  print * , 'The number is '//txt

  ! output:
  !     The number is 1234.57
  
  ! print a real nummber setting the number of decimals
  num = 3.14159265359d0
  txt = num2text(num,5)
  print * , 'The number is '//txt

  ! output:
  !     The number is 3.14159

  ! print an integer
  inum = 56789
  txt = num2text(inum)
  print * , 'The number is '//txt
  
  ! output:
  !     The number is 56789

  ! print an integer with a width of 8 (padded with zeros)
  txt = num2text(inum,8)
  print * , 'The number is '//txt

  ! output:
  !     The number is 00056789
```



