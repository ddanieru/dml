
all: gleg mc trap simp romberg

gleg: interfaces.f90 constants.f90 functions.f90 dml_int.f90 gleg.f90
	gfortran -o gleg interfaces.f90 constants.f90 functions.f90 dml_int.f90 gleg.f90
	./gleg

mc: interfaces.f90 constants.f90 functions.f90 dml_int.f90 mc.f90
	gfortran -o mc interfaces.f90 constants.f90 functions.f90 dml_int.f90 mc.f90
	./mc

trap:
	gfortran -o trap interfaces.f90 constants.f90 functions.f90 dml_int.f90 trapezoid.f90
	./trap

simp:
	gfortran -o simp interfaces.f90 constants.f90 functions.f90 dml_int.f90 simpson.f90
	./simp

romberg:
	gfortran -o romberg interfaces.f90 constants.f90 functions.f90 dml_int.f90 romberg.f90
	./romberg

clean:
	rm *mod gleg mc trap simp romberg #*out 
