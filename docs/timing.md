
# Fortran toolkit

###### Borja Petit

---

## ```timing```

```fortran
function timing(mode) result(time)
  implicit none
  integer , optional :: mode
  real(kind=8)       :: time
```

This functions returns a timing measure that is robust to parallelization. In particular, it returns the number of seconds since 00:00h of the 1st day of the month. The variable ```mode``` controls how time is measured:
- If ```mode``` = 1, ```time``` is measured in seconds (default).
- If ```mode``` = 2, ```time``` is measures in minutes.
- If ```mode``` = 3, ```time``` is measured in hours

**Dependencies**: [```error```](error.md)

[(back to index)](../index.md)

---

**Example**

Measure execution time in minutes:

```fortran
time0 = timing(2)

! ... some code ...

time1 = timing(2)

print * , time1 - time2 , 'min'
``` 