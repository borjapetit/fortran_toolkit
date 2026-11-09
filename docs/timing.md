
<span style="text-align:right;display:block;">
<a href="https://borjapetit.github.io/fortran_toolkit/">Back to index</a>
</span>

## ```timing```

```fortran
function timing(mode) result(time)
  implicit none
  integer , optional :: mode
  real(kind=8)       :: time
```

This function returns a (wall-clock) timing measure that is robust to parallelization. In particular, it computes the number of seconds, with millisecond resolution, since 00:00h UTC of January 1st, 1970, so that differences between two calls are correct across days, months and years. The variable ```mode``` controls how ```time``` is measured:
- If ```mode``` = 1, ```time``` is measured in seconds (default).
- If ```mode``` = 2, ```time``` is measures in minutes.
- If ```mode``` = 3, ```time``` is measured in hours

**Internal dependencies**: none

---

**Example**

Measure execution time in seconds, minutes and hours:

```fortran
time0 = timing( )

! ... some code ...

time1 = timing( )

print * , time1 - time0 , 'secs'

time0 = timing(2)

! ... some code ...

time1 = timing(2)

print * , time1 - time0 , 'min'

time0 = timing(3)

! ... some code ...

time1 = timing(3)

print * , time1 - time0 , 'hours'
``` 