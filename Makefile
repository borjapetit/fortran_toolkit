
# Start of the makefile

files = toolkit.f90 test_toolkit.f90

model: $(files)
	gfortran -O3 $(files) -o test

check: $(files)
	gfortran -fcheck=all -fbacktrace -Wall -g $(files) -o test

clean:
	rm -f *.mod
	rm -f test

# End of the makefile
