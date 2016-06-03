!***********************************************************

SUBROUTINE ave_err(datum,average,error,numprocs,myid)
  
  IMPLICIT integer(i-n)
  IMPLICIT double precision(a-h,o-y)
  IMPLICIT complex*16(z)
  INCLUDE "mpif.h"
  
  CALL MPI_REDUCE(datum,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  
  square = datum**2
  
  CALL MPI_REDUCE(square,sumsquares,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  
  IF (myid .eq. 0) THEN
     average = sum/numprocs
     avequad = sumsquares/numprocs
     error = dsqrt((avequad-average**2)/numprocs)
  END IF
  
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  
END SUBROUTINE ave_err

!***********************************************************
