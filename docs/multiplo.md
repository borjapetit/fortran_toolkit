
<span style="text-align:right;display:block;">
<a href="https://borjapetit.github.io/fortran_toolkit/">Back to index</a>
</span>

## ```multiplo```

```fortran
function multiplo(num,xx) result(mul)
  implicit none
  integer :: num,xx
  logical :: mul
```

This function checks whether a number ```num0``` is a multiple of ```num1```. The result ```mul``` is a logical variable taking value ```.TRUE.``` if ```num0``` is a multiple of ```num1```, and ```.FALSE.``` otherwise.

**Internal dependencies**: [```error```](error.md)

---

**Example**

We have an interative code and we want to print a message every 50 iterations

```fortran
do i = 1000
  ! ... some code ...
  if (multiplo(i,50)) print * , 'Iterations = ', i
end do
```
