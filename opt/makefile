
all:
	gfortran -o bis interfaces.f90 functions.f90 dml_opt.f90 bisect.f90
	gfortran -o reg interfaces.f90 functions.f90 dml_opt.f90 regfalsi.f90
	gfortran -o qint interfaces.f90 functions.f90 dml_opt.f90 qinterp.f90
	gfortran -o nroots interfaces.f90 functions.f90 dml_opt.f90 nroots.f90
	gfortran -o nrmin interfaces.f90 functions.f90 dml_opt.f90 nrmin.f90
	gfortran -o nreqs interfaces.f90 functions.f90 dml_opt.f90 nreqs.f90
	./bis
	./reg
	./qint
	./nrmin
	./nreqs

clean:
	rm *mod bis reg qint nroots nrmin nreqs # *out
