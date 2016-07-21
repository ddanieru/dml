
all: gram gaujor jacobi gaussel 

gram:
	gfortran -o gram dml_mat.f90 gram.f90

gaujor:
	gfortran -o gaujor dml_mat.f90 gaujor.f90
	
jacobi:
	gfortran -o jacobi dml_mat.f90 jacobi.f90

gaussel:
	gfortran -o gaussel dml_mat.f90 gaussel.f90

clean:
	rm *mod gram gaujor jacobi gaussel # *out 
