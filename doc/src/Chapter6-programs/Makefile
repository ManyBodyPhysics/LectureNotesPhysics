	F95 = gfortran
	FC_OPT = -fno-range-check
	LINK_OPT = 
	LAPACK_ITYPE = 
	MKL_PATH = 
	MKLLIB = -L/usr/local/lib -llapack -lblas

# *************

OBJ = MT.o new_solvers.o waveinit.o transfermatrix.o nuclei.o 

nuclei: $(OBJ)
	$(F95) $(LINK_OPT) -o $@ $(OBJ) $(MKLLIB) 

# *************

nuclei.o: nuclei.f90 improve.f90 input.f90
	$(F95) $(FC_OPT) -c nuclei.f90 

transfermatrix.o: transfermatrix.f90 improve.f90 input.f90
	$(F95) $(FC_OPT) -c transfermatrix.f90 

waveinit.o: waveinit.f90 input.f90
	$(F95) $(FC_OPT) -c waveinit.f90

new_solvers.o: new_solvers.f90
	$(F95) $(FC_OPT) $(LAPACK_ITYPE) -c new_solvers.f90

MT.o: MT.f90
	$(F95) $(FC_OPT) -c MT.f90

# *************

clean :
	rm -v -f nuclei nuclei_omp *.o 


