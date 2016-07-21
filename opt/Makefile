
all: bis reg qint nroots nrmin nreqs

bis:
	gfortran -o bis interfaces.f90 functions.f90 dml_opt.f90 bisect.f90
	./bis
reg:
	gfortran -o reg interfaces.f90 functions.f90 dml_opt.f90 regfalsi.f90
	./reg
qint:
	gfortran -o qint interfaces.f90 functions.f90 dml_opt.f90 qinterp.f90
	./qint
nroots:
	gfortran -o nroots interfaces.f90 functions.f90 dml_opt.f90 nroots.f90
nrmin:
	gfortran -o nrmin interfaces.f90 functions.f90 dml_opt.f90 nrmin.f90
	./nrmin
nreqs:
	gfortran -o nreqs interfaces.f90 functions.f90 dml_opt.f90 nreqs.f90
	./nreqs

clean:
	rm *mod bis reg qint nroots nrmin nreqs # *out
