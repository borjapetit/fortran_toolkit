
<span style="text-align:right;display:block;">
<a href="https://borjapetit.github.io/fortran_toolkit/">Back to index</a>
</span>

## ```vect```

```fortran
function vect(mat) result(vec)
  implicit none
  integer :: mat(:,:,...,:)
  integer :: vec(:)
```
```fortran
function vect(mat) result(vec)
  implicit none
  real(kind=8) :: mat(:,:,...,:)
  real(kind=8) :: vec(:)
```
```fortran
function vect(mat) result(vec)
  implicit none
  logical :: mat(:,:,...,:)
  logical :: vec(:)
```


This function returns a 1-dimensional array ```vec``` with all the elements of a user-supplied ```n```-dimensional array ```mat```, where $2\leq$ ```n``` $\leq5$. The input matrix can be double precision, integer or logical. The elements are taken in column-major order (as in ```reshape(mat,(/size(mat)/))```): first all the rows of the first column, then those of the second column, and so on.

**Internal dependencies**: none

**Note**: ```vect``` is an interface that calls specific functions depending on the type of ```mat``` (integer, real or logical) and depending on the dimensions of ```mat```. The specific functions are:
- ```mat``` is real $ \ \to \ $ ```vectorize_dp_2d```, ```vectorize_dp_3d```, ```vectorize_dp_4d```, ```vectorize_dp_5d```
- ```mat``` is integer $ \ \to \ $ ```vectorize_int_2d```, ```vectorize_int_3d```, ```vectorize_int_4d```, ```vectorize_int_5d```
- ```mat``` is logical $ \ \to \ $ ```vectorize_lo_2d```, ```vectorize_lo_3d```, ```vectorize_lo_4d```, ```vectorize_lo_5d```

---

**Example**:

```fortran
! matrix is integer

mat(:,1) = (/ 1 , 2 /)
mat(:,2) = (/ 3 , 4 /)

vec = vect(mat)

print * , 'vec =', vec   ! vec =  1  2  3  4
```

```fortran
! matrix is double precision

mat(:,1) = (/ 1.5d0 , 2.2d0 /)
mat(:,2) = (/ 2.4d0 , 0.7d0 /)

vec = vect(mat)

print * , 'vec =', vec   ! vec =  1.5000  2.2000  2.4000  0.7000
```


```fortran
! matrix is logical

mat(:,1) = (/ .true. , .false. /)
mat(:,2) = (/ .false. , .true. /)

vec = vect(mat)

print * , 'vec =', vec   ! vec =  T  F  F  T
```



