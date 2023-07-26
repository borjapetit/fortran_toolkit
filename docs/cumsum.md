# cumsum

```fortran
function cumsum(vec0) result(vec1)
  implicit none
  real(kind=8) :: vec0(:)
  real(kind=8) :: vec1(size(vec0))
```

Returns the cummulative sum of a vector ```vec0```.

**Dependencies**: none

[(back to index)](../index.md)

---

**Example**:

```fortran
vec0 = (/ 1.0, 2.0, 1.0, 3.0 /) 
vec1 = cumsum(vec0)

print * , 'vec1 = ', vec1   ! vec1 =  1.00  3.00  4.00  7.00
```


