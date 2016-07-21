
all:
	cd int; make 
	cd mat; make 
	cd opt; make
	#cd ..	

clean:
	cd int; make clean
	cd mat; make clean
	cd opt; make clean
	#cd ..	
