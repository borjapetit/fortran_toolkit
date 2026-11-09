
<span style="text-align:right;display:block;">
<a href="https://borjapetit.github.io/fortran_toolkit/">Back to index</a>
</span>

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

This subroutine takes an input vector ```input_vec``` and fills ```output_vec``` with the shuffled elements of ```input_vec``` using the Fisher-Yates algorithm. The two vectors do not need to have the same size:
- If ```size(output_vec)``` $\leq$ ```size(input_vec)```, ```output_vec``` contains the first ```size(output_vec)``` elements of a random permutation of ```input_vec``` (draws without replacement).
- If ```size(output_vec)``` > ```size(input_vec)```, the same random permutation of ```input_vec``` is repeated until ```output_vec``` is filled.

**Internal dependencies**: [```fisheryates```](fisheryates.md)

**Note**: ```shuffle_vect``` is an interface that calls specific functions depending on the type of ```input_vec``` and ```output_vec```(integer or double precision). The specific functions are:
- ```shuffle_vect_int``` if ```input_vec``` and ```output_vec``` are of type ```integer```
- ```shuffle_vect_dp``` if ```input_vec``` and ```output_vec``` are of type ```double precision```

---

**Example**:

```fortran
integer :: vect0(4)
integer :: vect1(4)

vect0 = (/ 20 , 35 , 43 , 67 /)

call shuffle_vect(vect0,vect1)
print * , 'vect1 =', vect1   ! vect1 =  43  20  67  35

call shuffle_vect(vect0,vect1)
print * , 'vect1 =', vect1   ! vect1 =  35  67  20  43
```