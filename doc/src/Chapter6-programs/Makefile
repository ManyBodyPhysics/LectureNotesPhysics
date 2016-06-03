
# *************
# Standard compiler options (mkl version 10.1.0.015 and similar):

F95 = ifort
F77 = ifort 
MPI_FC = mpif95

#FC_OPT = -warn -openmp
#MPI_FC_OPT = -warn -openmp

FC_OPT = -openmp 
MPI_FC_OPT = -openmp
LINK_OPT = -openmp

LAPACK_ITYPE = -i8
#MKL_PATH = /opt/intel/mkl/10.1.0.015/lib/em64t
MKL_PATH = /usr/local/apps/intel/composer_xe_2013_sp1.2.144/mkl/lib/intel64

LIBSTA = -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lmkl_lapack
LIBDYN = -lm

LIBSTA_OMP = -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -lmkl_lapack
LIBDYN_OMP = -lpthread -lm

MKLLIB = -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lm 
MKLLIB_STA = -L$(MKL_PATH) -Bstatic $(LIBSTA) -Bdynamic $(LIBDYN)

MKLLIB_OMP = -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -lpthread -lm 
MKLLIB_OMP_STA = -L$(MKL_PATH) -Bstatic $(LIBSTA_OMP) -Bdynamic $(LIBDYN_OMP)

# *************
# Options for the NCSU cluster:

ifeq ($(shell echo $${HOSTNAME%5*}),login)

	F95 = ifort
	F77 = ifort 
	MPI_FC = mpif95 

	FC_OPT = -O1 -openmp 
	MPI_FC_OPT = -O1 -static
	LINK_OPT = -openmp

	LAPACK_ITYPE = -i4
	MKL_PATH = /usr/local/intel/mkl91023/lib_ilp64/em64t

	LIBSTA = 
	LIBDYN = -lm

	LIBSTA_OMP = 
	LIBDYN_OMP = -lpthread -lm

	MKLLIB = -Bdynamic -L$(MKL_PATH) -lmkl -lmkl_lapack -lm
	MKLLIB_STA = -L$(MKL_PATH) -Bstatic $(LIBSTA) -Bdynamic $(LIBDYN)

	MKLLIB_OMP = -L$(MKL_PATH) -lmkl -lmkl_lapack -lpthread -lm 
	MKLLIB_OMP_STA = -L$(MKL_PATH) -Bstatic $(LIBSTA_OMP) -Bdynamic $(LIBDYN_OMP)

endif

ifeq ($(shell echo $${HOSTNAME%0*}),login)

	F95 = /opt/intel/composer_xe_2013_sp1.2.144/bin/intel64/ifort
	F77 = /opt/intel/composer_xe_2013_sp1.2.144/bin/intel64/ifort 
	MPI_FC = /usr/local/apps/openmpi/1.8.4-intel2013/GiGE/bin/mpif77

	FC_OPT = -O1 -openmp 
	MPI_FC_OPT = -O1 -static
	LINK_OPT = -openmp

	LAPACK_ITYPE = -i4
	MKL_PATH = /usr/local/apps/intel/composer_xe_2013_sp1.2.144/mkl/lib/intel64

	LIBSTA = 
	LIBDYN = -lm

	LIBSTA_OMP = 
	LIBDYN_OMP = -lpthread -lm

	MKLLIB = -Bdynamic -L$(MKL_PATH) -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lm
	MKLLIB_STA = -L$(MKL_PATH) -Bstatic $(LIBSTA) -Bdynamic $(LIBDYN)

	MKLLIB_OMP = -L$(MKL_PATH)  -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lm 
	MKLLIB_OMP_STA = -L$(MKL_PATH) -Bstatic $(LIBSTA_OMP) -Bdynamic $(LIBDYN_OMP)

endif
# *************
# Options for the Aachen RWTH cluster:

HOST_NAME := $(shell echo $$HOSTNAME | awk -F '.' '{print $$2 "." $$3 "." $$4}')
ifeq ($(HOST_NAME),rz.RWTH-Aachen.DE)

    F95 = ifort
    F77 = ifort
    MPI_FC = mpif95

    FC_OPT = -openmp
    MPI_FC_OPT = -openmp
    LINK_OPT = -openmp

    LAPACK_ITYPE = -i8
    MKL_PATH = /opt/intel/mkl/10.1.0.015/lib/em64t

    LIBSTA = -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lmkl_lapack
    LIBDYN = -lm

    LIBSTA_OMP = -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -lmkl_lapack
    LIBDYN_OMP = -lpthread -lm

    MKLLIB = -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lm 
    MKLLIB_STA = -L$(MKL_PATH) -Bstatic $(LIBSTA) -Bdynamic $(LIBDYN)

    MKLLIB_OMP = -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -lpthread -lm 
    MKLLIB_OMP_STA = -L$(MKL_PATH) -Bstatic $(LIBSTA_OMP) -Bdynamic $(LIBDYN_OMP)

endif

# *************
# Options for JUGENE:

ifeq ($(shell echo $${HOSTNAME%%[0-9]*}),jugene)

	F95 = bgxlf95_r
	F77 = bgxlf_r
	MPI_FC = mpixlf95_r

	FC_OPT = -O3 -qstrict -qarch=450 -qtune=450 -qnosave
	MPI_FC_OPT = -O3 -qstrict -qarch=450 -qtune=450 -qnosave
	LINK_OPT = -static

	LAPACK_ITYPE = -qintsize=4

	MKLLIB = -L/bgsys/local/lapack/lib -llapack -L/bgsys/local/lib -lesslbg
	MKLLIB_STA = -L/bgsys/local/lapack/lib -llapack -L/bgsys/local/lib -lesslbg

	MKLLIB_OMP = -L/bgsys/local/lapack/lib -llapack -L/bgsys/local/lib -lesslbg
	MKLLIB_OMP_STA = -L/bgsys/local/lapack/lib -llapack -L/bgsys/local/lib -lesslbg

endif

# *************
# Options for JUQUEEN:

ifeq ($(shell echo $${HOSTNAME%%[0-9]*}),juqueen)

	F95 = bgxlf95_r
	F77 = bgxlf_r
	MPI_FC = mpixlf95_r

	FC_OPT = -O4 -qarch=qp -qtune=qp
	MPI_FC_OPT = -O4 -qarch=qp -qtune=qp
	LINK_OPT = -static

	LAPACK_ITYPE = -qintsize=4

	MKLLIB = -L/bgsys/local/lapack/3.3.0/lib -llapack -L/bgsys/local/lib -lesslbg
	MKLLIB_STA = -L/bgsys/local/lapack/3.3.0/lib -llapack -L/bgsys/local/lib -lesslbg

	MKLLIB_OMP = -L/bgsys/local/lapack/3.3.0/lib -llapack -L/bgsys/local/lib -lesslbg
	MKLLIB_OMP_STA = -L/bgsys/local/lapack/3.3.0/lib -llapack -L/bgsys/local/lib -lesslbg

endif

# *************

OBJ = MT.o new_solvers.o waveinit.o analysis.o transfermatrix.o nuclei.o 

nuclei: $(OBJ)
	$(MPI_FC) $(LINK_OPT) -o $@ $(OBJ) $(MKLLIB) 

nuclei_sta: $(OBJ)
	$(MPI_FC) $(LINK_OPT) -o $@ $(OBJ) $(MKLLIB_STA)

nuclei_omp: $(OBJ)
	$(MPI_FC) $(LINK_OPT) -o $@ $(OBJ) $(MKLLIB_OMP)

nuclei_omp_sta: $(OBJ)
	$(MPI_FC) $(LINK_OPT) -o $@ $(OBJ) $(MKLLIB_OMP_STA)

# *************
# Unused components(?): printdualwave.inc

nuclei.o: nuclei.f90 improve.f90 input.f90
	$(MPI_FC) $(MPI_FC_OPT) -c nuclei.f90 

transfermatrix.o: transfermatrix.f90 improve.f90 input.f90
	$(F95) $(FC_OPT) -c transfermatrix.f90 

analysis.o: analysis.f90
	$(MPI_FC) $(MPI_FC_OPT) -c analysis.f90

waveinit.o: waveinit.f90 input.f90
	$(F95) $(FC_OPT) -c waveinit.f90

new_solvers.o: new_solvers.f90
	$(F95) $(FC_OPT) $(LAPACK_ITYPE) -c new_solvers.f90

MT.o: MT.f90
	$(F95) $(FC_OPT) -c MT.f90

# *************

clean :
	rm -v -f nuclei nuclei_omp *.o 

nuclei1:  nuclei
	cp nuclei nuclei1

nuclei2:  nuclei
	cp nuclei nuclei2

nuclei3:  nuclei
	cp nuclei nuclei3

nuclei4:  nuclei
	cp nuclei nuclei4

nuclei5:  nuclei
	cp nuclei nuclei5

nuclei6:  nuclei
	cp nuclei nuclei6

nuclei7:  nuclei
	cp nuclei nuclei7

nuclei8:  nuclei
	cp nuclei nuclei8

nuclei9:  nuclei
	cp nuclei nuclei9

nuclei10:  nuclei
	cp nuclei nuclei10
