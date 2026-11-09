
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

This function checks whether an integer ```num``` is a multiple of ```xx```. The result ```mul``` is a logical variable taking value ```.TRUE.``` if ```num``` is a multiple of ```xx```, and ```.FALSE.``` otherwise. If ```xx``` = 0, the function prints an error message and returns ```.FALSE.```.

**Internal dependencies**: [```error```](error.md)

---

**Example**

We have an iterative code and we want to print a message every 50 iterations

```fortran
do i = 1,1000
  ! ... some code ...
  if (multiplo(i,50)) print * , 'Iterations = ', i
end do
```
