
## ```shuffle_vect```

```fortran
subroutine shuffle_vect(input_vec, output_vec)
  implicit none
  integer , intent(in)  :: input_vec(:)
  integer , intent(out) :: output_vec(:)
  ! or 
  real(kind=8) , intent(in)  :: input_vec(:)
  real(kind=8) , intent(out) :: output_vec(:)
```

This rubroutine takes an input vector ```input_vec``` and fills ```output_vec``` with the shuffled elements of ```input_vec``` using the Fisher-Yates algorithm.

**Dependencies**: ```fisheryates```

**Note**: ```shuffle_vect``` is an interface that calls specific functions depending on the type of ```input_vec``` and ```output_vec```(integer or double precision). The specific functions are:
- ```shuffle_vect_int``` if ```input_vec``` and ```output_vec``` are of type ```integer```
- ```shuffle_vect_dp``` if ```input_vec``` and ```output_vec``` are of type ```double precision```

[(back to index)](../index.md)

---

**Example**:

```fortran
integer :: vect0(4)
integer :: vect1(4)

vect0 = (/ 20 , 35 , 43 , 67 /)

call shuffle_vect(vect0,vect1)
print * , 'vect1 =', vect1   ! vect1 =  43  20  67  35

call shuffle_vect(vect)
print * , 'vect1 =', vect1   ! vect1 =  35  67  20  43
```