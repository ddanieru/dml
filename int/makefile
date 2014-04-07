
all:
	gfortran -o gleg interfaces.f90 constants.f90 functions.f90 dml_int.f90 gleg.f90
	gfortran -o mc interfaces.f90 constants.f90 functions.f90 dml_int.f90 mc.f90
	gfortran -o trap interfaces.f90 constants.f90 functions.f90 dml_int.f90 trapezoid.f90
	gfortran -o simp interfaces.f90 constants.f90 functions.f90 dml_int.f90 simpson.f90
	gfortran -o romberg interfaces.f90 constants.f90 functions.f90 dml_int.f90 romberg.f90
	./gleg
	./mc
	./trap
	./simp
	./romberg

clean:
	rm *mod gleg mc trap simp romberg #*out 
