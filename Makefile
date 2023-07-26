
# Start of the makefile

model:
	gfortran -O3 src/toolkit.f90 src/test_toolkit.f90 -o src/test_gral
	gfortran -O3 src/toolkit.f90 src/test_optim.f90 -o src/test_optim

check: $(files)
	gfortran -fcheck=all -fbacktrace -Wall -g -O3 src/toolkit.f90 src/test_toolkit.f90 -o src/test_gral_check
	gfortran -fcheck=all -fbacktrace -Wall -g -O3 src/toolkit.f90 src/test_optim.f90 -o src/test_optim_check


clean:
	rm -f *.mod
	rm -f test

# End of the makefile
