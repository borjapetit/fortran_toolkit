<style TYPE="text/css">
code.has-jax {font: inherit; font-size: 100%; background: inherit; border: inherit;}
</style>
<script type="text/x-mathjax-config">
MathJax.Hub.Config({
    tex2jax: {
        inlineMath: [['$','$'], ['\\(','\\)']],
        skipTags: ['script', 'noscript', 'style', 'textarea', 'pre'] // removed 'code' entry
    }
});
MathJax.Hub.Queue(function() {
    var all = MathJax.Hub.getAllJax(), i;
    for(i = 0; i < all.length; i += 1) {
        all[i].SourceElement().parentNode.className += ' has-jax';
    }
});
</script>
<script type="text/javascript" src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.4/MathJax.js?config=TeX-AMS_HTML-full"></script>

### vect

```fortran
function vec(mat) result(vec)
  implicit none
  real(kind=8)  :: mat(:,:,...,:)
  real(kind=8)  :: vec(:)
```

This function returns a 1-dimensional array ```vec``` with all the elements of a user-supplied ```n```-dimensional array ```mat```, where ```n```$\leq5$.

**Dependencies**: none

**Note**: ```vect``` is an interface that calls specific functions depending on the type of ```mat``` (integer or real) and depending on the dimensions of ```mat```. The specific functions are:
- ```mat``` is real $ \ \to \ $ ```vectorize_dp_2d```, ```vectorize_dp_3d```, ```...```, ```vectorize_in_6d```
- ```mat``` is integer $ \ \to \ $ ```vectorize_in_2d```, ```vectorize_in_3d```, ```...```, ```vectorize_in_6d```

[(back to index)](../index.md)

---

**Example**:

```fortran
mat(:,1) = (/ 1 , 2 /)
mat(:,2) = (/ 3 , 4 /)

vec = vect(mat)

print * , 'vec =', vec   ! vec =  1.00  3.00  2.00  4.00
```




