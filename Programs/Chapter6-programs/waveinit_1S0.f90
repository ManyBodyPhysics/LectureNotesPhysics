!***********************************************************

SUBROUTINE waveinit(zwave,zdualwave,myid)

  IMPLICIT integer(i-n)
  IMPLICIT double precision(a-h,o-y)
  IMPLICIT complex*16(z)
  INCLUDE "input.f90"      
  DIMENSION zwave(0:L-1,0:L-1,0:L-1,0:1,0:1,0:n_f-1)
  DIMENSION zdualwave(0:L-1,0:L-1,0:L-1,0:1,0:1,0:n_f-1)
  DIMENSION px(0:n_f-1), py(0:n_f-1), pz(0:n_f-1)
  CHARACTER*3 trig(0:n_f-1)
  CHARACTER*2 particletype(0:n_f-1)
  
  IF (n_f .ne. 2) THEN
     WRITE(*,*)'error '
     WRITE(*,*)' n_f ',n_f
     STOP
  END IF
  
  DO npart = 0,n_f-1
     DO ni = 0,1; DO ns = 0,1
        DO nz = 0,L-1; DO ny = 0,L-1; DO nx = 0,L-1
           zwave(nx,ny,nz,ns,ni,npart) = 0.D0
        END DO; END DO; END DO
     END DO; END DO
  END DO
  
  DO npart = 0,1
     px(npart) = 0.D0
     py(npart) = 0.D0
     pz(npart) = 0.D0
     trig(npart) = 'one'
  END DO
  particletype(0) = 'n+'
  particletype(1) = 'n-'
  
  DO npart = 0, n_f-1         
     IF (particletype(npart) .eq. 'p+') THEN
        ns = 0
        ni = 0
     ELSE IF (particletype(npart) .eq. 'p-') THEN
        ns = 1
        ni = 0
     ELSE IF (particletype(npart) .eq. 'n+') THEN
        ns = 0
        ni = 1
     ELSE IF (particletype(npart) .eq. 'n-') THEN
        ns = 1
        ni = 1
     ELSE
        WRITE(*,*)'particletype error',npart,particletype(npart)
        STOP
     END IF
     IF (trig(npart) .eq. 'one') THEN                                      
        DO nz = 0,L-1; DO ny = 0,L-1; DO nx = 0,L-1                     
           zwave(nx,ny,nz,ns,ni,npart) = 1.D0
        END DO; END DO; END DO
     ELSE IF (trig(npart) .eq. 'cos') THEN         
        DO nz = 0,L-1; DO ny = 0,L-1; DO nx = 0,L-1                     
           zwave(nx,ny,nz,ns,ni,npart) = &
                dcos(nx*px(npart)*2*pi/L &
                + ny*py(npart)*2*pi/L &
                + nz*pz(npart)*2*pi/L)*dsqrt(2.D0) 
        END DO; END DO; END DO
     ELSE IF (trig(npart) .eq. 'sin') THEN
        DO nz = 0,L-1; DO ny = 0,L-1; DO nx = 0,L-1
           zwave(nx,ny,nz,ns,ni,npart) = & 
                dsin(nx*px(npart)*2*pi/L &
                + ny*py(npart)*2*pi/L &
                + nz*pz(npart)*2*pi/L)*dsqrt(2.D0) 
        END DO; END DO; END DO
     ELSE
        WRITE(*,*)'trig function error',npart
        STOP
     END IF
     
  END DO
  
  IF (myid .eq. 0) THEN
     DO npart = 0,n_f-1
        WRITE(*,*)''
        WRITE(*,*)npart,' ',particletype(npart),' ',trig(npart)
        WRITE(*,*)int(px(npart)),int(py(npart)),int(pz(npart))
     END DO
     WRITE(*,*)''
  END IF
  
  DO npart = 0,n_f-1
     DO ni = 0,1; DO ns = 0,1
        DO nz = 0,L-1; DO ny = 0,L-1; DO nx = 0,L-1
           zdualwave(nx,ny,nz,ns,ni,npart) &
                = dconjg(zwave(nx,ny,nz,ns,ni,npart))
        END DO; END DO; END DO
     END DO; END DO
  END DO
  
END SUBROUTINE waveinit

!***********************************************************
