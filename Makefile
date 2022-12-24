
# Start of the makefile

model:
	gfortran -O3 src/toolkit.f90 src/test_toolkit.f90 -o src/test_gral
	gfortran -O3 src/toolkit.f90 src/test_optim.f90 -o src/test_optim

clean:
	rm -f *.mod
	rm -f test

# End of the makefile
