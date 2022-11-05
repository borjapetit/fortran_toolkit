
## timing

```fortran
function timing(mode) result(time)
  implicit none
  integer , optional :: mode
  real(kind=8)       :: time
```

This functions returns a timing number that is robust to parallel computing. In particular, it returns the number of seconds since 00:00h of the 1st day of the month. The variable ```mode``` controls how time is measured:

- If ```mode``` = 1, time is measured in seconds (default).
- If ```mode``` = 2, time is measures in minutes.
- If ```mode``` = 3, time is measured in hours

**Dependencies**: none

[(back to index)](index.md)

