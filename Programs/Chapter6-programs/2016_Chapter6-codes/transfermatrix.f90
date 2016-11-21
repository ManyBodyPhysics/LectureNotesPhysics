!***********************************************************

SUBROUTINE getzvecs(s,sI,zvecs,zvecsinit,nt2,nt1, &
     pion,ztau2x2,num)
  
  IMPLICIT INTEGER(i-n)
  IMPLICIT DOUBLE PRECISION(a-h,o-y)
  IMPLICIT COMPLEX*16(z)
  INCLUDE "input.f90"
  
  DIMENSION s(0:L-1,0:L-1,0:L-1,0:Lt-1)
  DIMENSION sI(0:L-1,0:L-1,0:L-1,0:Lt-1,1:3)
  DIMENSION zvecs(0:L-1,0:L-1,0:L-1,0:Lt,0:1,0:1,0:num-1)
  DIMENSION zvecsinit(0:L-1,0:L-1,0:L-1,0:1,0:1,0:num-1)
  DIMENSION pion(0:L-1,0:L-1,0:L-1,0:Lt-1,1:3)
  DIMENSION zpi2x2(0:1,0:1)
  DIMENSION zsI2x2(0:1,0:1)
  DIMENSION zsSI2x2(0:1,0:1,0:1,0:1)
  DIMENSION ztau2x2(0:1,0:1,0:3)
  
  !   nt2 > nt1
  
  INCLUDE "improve.f90"
  
  IF (mod(L,2) == 0) THEN
     pp = pi
  ELSE
     pp = (L-1)*pi/L
  END IF
  
  DO ni = 0,1; DO ns = 0,1
     DO nz = 0,L-1; DO ny = 0,L-1; DO nx = 0,L-1
        DO npart = 0,num-1
           zvecs(nx,ny,nz,nt1,ns,ni,npart) = zvecsinit(nx,ny,nz,ns,ni,npart)
        END DO
     END DO; END DO; END DO
  END DO; END DO
  
  DO nt = nt1+1, nt2     
     DO np = 0,num-1        
        DO nz = 0,L-1; DO ny = 0,L-1; DO nx = 0,L-1           
           DO ni = 0,1; DO ns = 0,1
              
              zvecs(nx,ny,nz,nt,ns,ni,np) = zvecs(nx,ny,nz,nt-1,ns,ni,np) &
                   * (1.D0-6.D0*w0_N*h+CDSQRT(-c0*atovera*(1.D0,0.D0))*s(nx,ny,nz,nt-1))
              
              zvecs(nx,ny,nz,nt,ns,ni,np) = zvecs(nx,ny,nz,nt,ns,ni,np) &
                   + w1_N*h*zvecs(MOD(nx+1,L),ny,nz,nt-1,ns,ni,np) &
                   + w1_N*h*zvecs(MOD(nx-1+L,L),ny,nz,nt-1,ns,ni,np) &
                   + w1_N*h*zvecs(nx,MOD(ny+1,L),nz,nt-1,ns,ni,np) &
                   + w1_N*h*zvecs(nx,MOD(ny-1+L,L),nz,nt-1,ns,ni,np) &
                   + w1_N*h*zvecs(nx,ny,MOD(nz+1,L),nt-1,ns,ni,np) &
                   + w1_N*h*zvecs(nx,ny,MOD(nz-1+L,L),nt-1,ns,ni,np)
              
              IF (improveN >= 1) THEN                 
                 zvecs(nx,ny,nz,nt,ns,ni,np) = zvecs(nx,ny,nz,nt,ns,ni,np) & 
                      - w2_N*h*zvecs(MOD(nx+2,L),ny,nz,nt-1,ns,ni,np) &
                      - w2_N*h*zvecs(MOD(nx-2+L,L),ny,nz,nt-1,ns,ni,np) &
                      - w2_N*h*zvecs(nx,MOD(ny+2,L),nz,nt-1,ns,ni,np) &
                      - w2_N*h*zvecs(nx,MOD(ny-2+L,L),nz,nt-1,ns,ni,np) &
                      - w2_N*h*zvecs(nx,ny,MOD(nz+2,L),nt-1,ns,ni,np) &
                      - w2_N*h*zvecs(nx,ny,MOD(nz-2+L,L),nt-1,ns,ni,np)
              END IF
              
              IF (improveN == 2) THEN
                 zvecs(nx,ny,nz,nt,ns,ni,np) = zvecs(nx,ny,nz,nt,ns,ni,np) &
                      + w3_N*h*zvecs(MOD(nx+3,L),ny,nz,nt-1,ns,ni,np) &
                      + w3_N*h*zvecs(MOD(nx-3+L,L),ny,nz,nt-1,ns,ni,np) &
                      + w3_N*h*zvecs(nx,MOD(ny+3,L),nz,nt-1,ns,ni,np) &
                      + w3_N*h*zvecs(nx,MOD(ny-3+L,L),nz,nt-1,ns,ni,np) &
                      + w3_N*h*zvecs(nx,ny,MOD(nz+3,L),nt-1,ns,ni,np) &
                      + w3_N*h*zvecs(nx,ny,MOD(nz-3+L,L),nt-1,ns,ni,np)
              END IF
              
           END DO; END DO
        END DO; END DO; END DO
     END DO
     
     DO nz = 0,L-1; DO ny = 0,L-1; DO nx = 0,L-1

        DO nii = 0,1; DO ni = 0,1
           DO nss = 0,1; DO ns = 0,1            
              zsSI2x2(ns,nss,ni,nii) = (0.D0,0.D0)
           END DO; END DO
        END DO; END DO
        
        zsI2x2(0,0) = sI(nx,ny,nz,nt-1,3)
        zsI2x2(1,1) = -sI(nx,ny,nz,nt-1,3)
        zsI2x2(0,1) = sI(nx,ny,nz,nt-1,1) - (0.D0,1.D0)*sI(nx,ny,nz,nt-1,2)
        zsI2x2(1,0) = sI(nx,ny,nz,nt-1,1) + (0.D0,1.D0)*sI(nx,ny,nz,nt-1,2)
        
        DO nii = 0,1; DO ni = 0,1
           DO ns = 0,1
              zsSI2x2(ns,ns,ni,nii) = zsSI2x2(ns,ns,ni,nii) &                    
                   + (0.D0,1.D0)*CDSQRT(cI*atovera*(1.D0,0.D0))*zsI2x2(ni,nii)
           END DO
        END DO; END DO
        
        DO iso = 1,3
           
           pi1 = 0.D0 &
                + o1/2.D0*pion(MOD(nx+1,L),ny,nz,nt-1,iso) &
                - o1/2.D0*pion(MOD(nx-1+L,L),ny,nz,nt-1,iso) &
                - o2/2.D0*pion(MOD(nx+2,L),ny,nz,nt-1,iso) &
                + o2/2.D0*pion(MOD(nx-2+L,L),ny,nz,nt-1,iso) &
                + o3/2.D0*pion(MOD(nx+3,L),ny,nz,nt-1,iso) &
                - o3/2.D0*pion(MOD(nx-3+L,L),ny,nz,nt-1,iso)
           
           pi2 = 0.D0 &
                + o1/2.D0*pion(nx,MOD(ny+1,L),nz,nt-1,iso) &
                - o1/2.D0*pion(nx,MOD(ny-1+L,L),nz,nt-1,iso) &
                - o2/2.D0*pion(nx,MOD(ny+2,L),nz,nt-1,iso) &
                + o2/2.D0*pion(nx,MOD(ny-2+L,L),nz,nt-1,iso) &
                + o3/2.D0*pion(nx,MOD(ny+3,L),nz,nt-1,iso) &
                - o3/2.D0*pion(nx,MOD(ny-3+L,L),nz,nt-1,iso)
           
           pi3 = 0.D0 &
                + o1/2.D0*pion(nx,ny,MOD(nz+1,L),nt-1,iso) &
                - o1/2.D0*pion(nx,ny,MOD(nz-1+L,L),nt-1,iso) &
                - o2/2.D0*pion(nx,ny,MOD(nz+2,L),nt-1,iso) &
                + o2/2.D0*pion(nx,ny,MOD(nz-2+L,L),nt-1,iso) &
                + o3/2.D0*pion(nx,ny,MOD(nz+3,L),nt-1,iso) &
                - o3/2.D0*pion(nx,ny,MOD(nz-3+L,L),nt-1,iso)
           
           zpi2x2(0,0) = pi3
           zpi2x2(1,1) = -pi3
           zpi2x2(0,1) = pi1 - (0.D0,1.D0)*pi2
           zpi2x2(1,0) = pi1 + (0.D0,1.D0)*pi2
           
           DO nii = 0,1; DO ni = 0,1
              DO nss = 0,1; DO ns = 0,1
                 zsSI2x2(ns,nss,ni,nii) = zsSI2x2(ns,nss,ni,nii) &
                      - gA*atovera/(2.D0*fpi*SQRT(qpi3))*zpi2x2(ns,nss) &
                      * ztau2x2(ni,nii,iso)                          
              END DO; END DO
           END DO; END DO
           
        END DO
        
        DO np = 0,num-1
           DO nii = 0,1; DO ni = 0,1
              DO nss = 0,1; DO ns = 0,1                                  
                 zvecs(nx,ny,nz,nt,ns,ni,np) = zvecs(nx,ny,nz,nt,ns,ni,np) &
                      + zsSI2x2(ns,nss,ni,nii)*zvecs(nx,ny,nz,nt-1,nss,nii,np)            
              END DO; END DO
           END DO; END DO      
        END DO
        
     END DO; END DO; END DO          
     
  END DO

END SUBROUTINE getzvecs

!***********************************************************

SUBROUTINE getzdualvecs(s,sI,zdualvecs,zdualvecsinit, &
           nt2,nt1,pion,ztau2x2,num)

  IMPLICIT INTEGER(i-n)
  IMPLICIT DOUBLE PRECISION(a-h,o-y)
  IMPLICIT COMPLEX*16(z)
  INCLUDE "input.f90"
  
  DIMENSION s(0:L-1,0:L-1,0:L-1,0:Lt-1)
  DIMENSION sI(0:L-1,0:L-1,0:L-1,0:Lt-1,1:3)
  DIMENSION zdualvecs(0:L-1,0:L-1,0:L-1,0:Lt,0:1,0:1,0:num-1)
  DIMENSION zdualvecsinit(0:L-1,0:L-1,0:L-1,0:1,0:1,0:num-1)
  DIMENSION pion(0:L-1,0:L-1,0:L-1,0:Lt-1,1:3)
  DIMENSION zpi2x2(0:1,0:1)
  DIMENSION zsI2x2(0:1,0:1)
  DIMENSION zsSI2x2(0:1,0:1,0:1,0:1)
  DIMENSION ztau2x2(0:1,0:1,0:9)
  
  !   nt2 > nt1
  
  INCLUDE "improve.f90"
  
  IF (mod(L,2) == 0) THEN
     pp = pi
  ELSE
     pp = (L-1)*pi/L
  END IF
  
  DO ni = 0,1; DO ns = 0,1
     DO nz = 0,L-1; DO ny = 0,L-1; DO nx = 0,L-1
        DO npart = 0,num-1
           zdualvecs(nx,ny,nz,nt2,ns,ni,npart) = zdualvecsinit(nx,ny,nz,ns,ni,npart)
        END DO
     END DO; END DO; END DO
  END DO; END DO
  
  DO nt = nt2,nt1+1,-1
     
     DO np = 0,num-1
         DO nz = 0,L-1; DO ny = 0,L-1; DO nx = 0,L-1
            DO ni = 0,1; DO ns = 0,1
               
               zdualvecs(nx,ny,nz,nt-1,ns,ni,np) = zdualvecs(nx,ny,nz,nt,ns,ni,np) &
                    * (1.D0-6.D0*w0_N*h+CDSQRT(-c0*atovera*(1.D0,0.D0))*s(nx,ny,nz,nt-1))
               
               zdualvecs(nx,ny,nz,nt-1,ns,ni,np) = zdualvecs(nx,ny,nz,nt-1,ns,ni,np) & 
                    + w1_N*h*zdualvecs(MOD(nx+1,L),ny,nz,nt,ns,ni,np) &
                    + w1_N*h*zdualvecs(MOD(nx-1+L,L),ny,nz,nt,ns,ni,np) &
                    + w1_N*h*zdualvecs(nx,MOD(ny+1,L),nz,nt,ns,ni,np) &
                    + w1_N*h*zdualvecs(nx,MOD(ny-1+L,L),nz,nt,ns,ni,np) &
                    + w1_N*h*zdualvecs(nx,ny,MOD(nz+1,L),nt,ns,ni,np) &
                    + w1_N*h*zdualvecs(nx,ny,MOD(nz-1+L,L),nt,ns,ni,np) 
               
               IF (improveN >= 1) THEN
                  zdualvecs(nx,ny,nz,nt-1,ns,ni,np) = zdualvecs(nx,ny,nz,nt-1,ns,ni,np) &
                       - w2_N*h*zdualvecs(MOD(nx+2,L),ny,nz,nt,ns,ni,np) &
                       - w2_N*h*zdualvecs(MOD(nx-2+L,L),ny,nz,nt,ns,ni,np) &
                       - w2_N*h*zdualvecs(nx,MOD(ny+2,L),nz,nt,ns,ni,np) &
                       - w2_N*h*zdualvecs(nx,MOD(ny-2+L,L),nz,nt,ns,ni,np) &
                       - w2_N*h*zdualvecs(nx,ny,MOD(nz+2,L),nt,ns,ni,np) &
                       - w2_N*h*zdualvecs(nx,ny,MOD(nz-2+L,L),nt,ns,ni,np) 
               END IF
               
               IF (improveN == 2) THEN
                  zdualvecs(nx,ny,nz,nt-1,ns,ni,np) = zdualvecs(nx,ny,nz,nt-1,ns,ni,np) &
                       + w3_N*h*zdualvecs(MOD(nx+3,L),ny,nz,nt,ns,ni,np) &
                       + w3_N*h*zdualvecs(MOD(nx-3+L,L),ny,nz,nt,ns,ni,np) &
                       + w3_N*h*zdualvecs(nx,MOD(ny+3,L),nz,nt,ns,ni,np) &
                       + w3_N*h*zdualvecs(nx,MOD(ny-3+L,L),nz,nt,ns,ni,np) &
                       + w3_N*h*zdualvecs(nx,ny,MOD(nz+3,L),nt,ns,ni,np) &
                       + w3_N*h*zdualvecs(nx,ny,MOD(nz-3+L,L),nt,ns,ni,np) 
               END IF
               
            END DO; END DO
         END DO; END DO; END DO
      END DO
   
      DO nz = 0,L-1; DO ny = 0,L-1; DO nx = 0,L-1
         
         DO nii = 0,1; DO ni = 0,1
            DO nss = 0,1; DO ns = 0,1            
               zsSI2x2(ns,nss,ni,nii) = (0.D0,0.D0)
            END DO; END DO
         END DO; END DO
         
         zsI2x2(0,0) = sI(nx,ny,nz,nt-1,3)
         zsI2x2(1,1) = -sI(nx,ny,nz,nt-1,3)
         zsI2x2(0,1) = sI(nx,ny,nz,nt-1,1) - (0.D0,1.D0)*sI(nx,ny,nz,nt-1,2)
         zsI2x2(1,0) = sI(nx,ny,nz,nt-1,1) + (0.D0,1.D0)*sI(nx,ny,nz,nt-1,2)
         
         DO nii = 0,1; DO ni = 0,1
            DO ns = 0,1
               zsSI2x2(ns,ns,ni,nii) = zsSI2x2(ns,ns,ni,nii) &                    
                    + (0.D0,1.D0)*CDSQRT(cI*atovera*(1.D0,0.D0))*zsI2x2(ni,nii)
            END DO
         END DO; END DO

         DO iso = 1,3
            
            pi1 = 0.D0 &
                 + o1/2.D0*pion(MOD(nx+1,L),ny,nz,nt-1,iso) &
                 - o1/2.D0*pion(MOD(nx-1+L,L),ny,nz,nt-1,iso) &
                 - o2/2.D0*pion(MOD(nx+2,L),ny,nz,nt-1,iso) &
                 + o2/2.D0*pion(MOD(nx-2+L,L),ny,nz,nt-1,iso) &
                 + o3/2.D0*pion(MOD(nx+3,L),ny,nz,nt-1,iso) &
                 - o3/2.D0*pion(MOD(nx-3+L,L),ny,nz,nt-1,iso)
            
            pi2 = 0.D0 &
                 + o1/2.D0*pion(nx,MOD(ny+1,L),nz,nt-1,iso) &
                 - o1/2.D0*pion(nx,MOD(ny-1+L,L),nz,nt-1,iso) &
                 - o2/2.D0*pion(nx,MOD(ny+2,L),nz,nt-1,iso) &
                 + o2/2.D0*pion(nx,MOD(ny-2+L,L),nz,nt-1,iso) &
                 + o3/2.D0*pion(nx,MOD(ny+3,L),nz,nt-1,iso) &
                 - o3/2.D0*pion(nx,MOD(ny-3+L,L),nz,nt-1,iso)
            
            pi3 = 0.D0 &
                 + o1/2.D0*pion(nx,ny,MOD(nz+1,L),nt-1,iso) &
                 - o1/2.D0*pion(nx,ny,MOD(nz-1+L,L),nt-1,iso) &
                 - o2/2.D0*pion(nx,ny,MOD(nz+2,L),nt-1,iso) &
                 + o2/2.D0*pion(nx,ny,MOD(nz-2+L,L),nt-1,iso) &
                 + o3/2.D0*pion(nx,ny,MOD(nz+3,L),nt-1,iso) &
                 - o3/2.D0*pion(nx,ny,MOD(nz-3+L,L),nt-1,iso)
            
            zpi2x2(0,0) = pi3
            zpi2x2(1,1) = -pi3
            zpi2x2(0,1) = pi1 - (0.D0,1.D0)*pi2
            zpi2x2(1,0) = pi1 + (0.D0,1.D0)*pi2
            
            DO nii = 0,1; DO ni = 0,1
               DO nss = 0,1; DO ns = 0,1
                  
                  zsSI2x2(ns,nss,ni,nii) = zsSI2x2(ns,nss,ni,nii) &
                       - gA*atovera/(2.D0*fpi*SQRT(qpi3))*zpi2x2(ns,nss) &
                       * ztau2x2(ni,nii,iso)                          
                  
               END DO; END DO
            END DO; END DO
            
         END DO

         DO np = 0,num-1
            DO nii = 0,1; DO ni = 0,1
               DO nss = 0,1; DO ns = 0,1                                  
                  zdualvecs(nx,ny,nz,nt-1,nss,nii,np) = zdualvecs(nx,ny,nz,nt-1,nss,nii,np) &
                       + zsSI2x2(ns,nss,ni,nii)*zdualvecs(nx,ny,nz,nt,ns,ni,np)            
               END DO; END DO
            END DO; END DO      
         END DO
         
      END DO; END DO; END DO          
      
   END DO
   
END SUBROUTINE getzdualvecs

!**********************************************************

SUBROUTINE getmidoverlap(zvecs,zdualvecs,zdeter,aldeterabs,nt1,nt2)

IMPLICIT integer(i-n)
IMPLICIT double precision(a-h,o-y)
IMPLICIT complex*16(z)
INCLUDE "input.f90"

DIMENSION zvecs(0:L-1,0:L-1,0:L-1,0:Lt,0:1,0:1,0:n_f-1)
DIMENSION zdualvecs(0:L-1,0:L-1,0:L-1,0:Lt,0:1,0:1,0:n_f-1)
DIMENSION zwave(0:L-1,0:L-1,0:L-1,0:1,0:1,0:n_f-1)
DIMENSION zcorrmatrix(0:n_f-1,0:n_f-1)
DIMENSION indx(0:n_f-1)

DO np1 = 0,n_f-1; DO np2 = 0,n_f-1
   zcorrmatrix(np2,np1) = 0.D0
   DO ni = 0,1; DO ns = 0,1
      DO nz = 0,L-1; DO ny = 0,L-1; DO nx = 0,L-1               
         zcorrmatrix(np2,np1) = &
              zcorrmatrix(np2,np1) + &
              zdualvecs(nx,ny,nz,nt2,ns,ni,np2) &
              *zvecs(nx,ny,nz,nt1,ns,ni,np1)/L**3
      END DO; END DO; END DO
   END DO; END DO
END DO; END DO

CALL zludcmp_ratio(zcorrmatrix,zdeter,aldeterabs,n_f)

END SUBROUTINE getmidoverlap

!**********************************************************

SUBROUTINE getinvcorr(zvecs,zdualvecs,zldeter,zcorrmatrix,zcorrinv,nt)

IMPLICIT integer(i-n)
IMPLICIT double precision(a-h,o-y)
IMPLICIT complex*16(z)
INCLUDE "input.f90"

DIMENSION zvecs(0:L-1,0:L-1,0:L-1,0:Lt,0:1,0:1,0:n_f-1)
DIMENSION zdualvecs(0:L-1,0:L-1,0:L-1,0:Lt,0:1,0:1,0:n_f-1)
DIMENSION zcorrmatrix(0:n_f-1,0:n_f-1)
DIMENSION zcorrinv(0:n_f-1,0:n_f-1)

DO np1 = 0,n_f-1; DO np2 = 0,n_f-1
   zcorrmatrix(np2,np1) = 0.D0
   DO ni = 0,1; DO ns = 0,1
      DO nz = 0,L-1; DO ny = 0,L-1; DO nx = 0,L-1               
         zcorrmatrix(np2,np1) = &
              zcorrmatrix(np2,np1) + &
              zdualvecs(nx,ny,nz,nt,ns,ni,np2)* &
              zvecs(nx,ny,nz,nt,ns,ni,np1)/L**3
      END DO; END DO; END DO
   END DO; END DO
END DO; END DO

CALL zdoinv_log(zcorrmatrix,zcorrinv,zldeter,n_f)      

END SUBROUTINE getinvcorr

! **********************************************************

SUBROUTINE dV(zvecs,zdualvecs,ztau2x2,zcorrinv,zdVall)

IMPLICIT integer(i-n)
IMPLICIT double precision(a-h,o-y)
IMPLICIT complex*16(z)
INCLUDE "input.f90"

DIMENSION zvecs(0:L-1,0:L-1,0:L-1,0:Lt,0:1,0:1,0:n_f-1)
DIMENSION zdualvecs(0:L-1,0:L-1,0:L-1,0:Lt,0:1,0:1,0:n_f-1)
DIMENSION ztau2x2(0:1,0:1,0:3)
DIMENSION zcorrinv(0:n_f-1,0:n_f-1)
DIMENSION zdVall(0:L-1,0:L-1,0:L-1,0:Lt-1,0:3,0:3)
DIMENSION zmat(0:1,0:1,0:1,0:1)

DO nt = 0,Lt-1
   
   DO nz = 0,L-1; DO ny = 0,L-1; DO nx = 0,L-1          
      
      DO nii = 0,1; DO ni = 0,1; DO nss = 0,1; DO ns = 0,1                              
         zmat(ns,nss,ni,nii) = 0.D0
         DO np2 = 0,n_f-1; DO np1 = 0,n_f-1
            zmat(ns,nss,ni,nii) = &
                 zmat(ns,nss,ni,nii) &  
                 + zdualvecs(nx,ny,nz,nt+1, &
                 ns,ni,np1) &
                 *zvecs(nx,ny,nz,nt,nss,nii,np2) &
                 *zcorrinv(np2,np1)
         END DO; END DO
      END DO; END DO; END DO; END DO
      
      DO iso = 0,3; DO nsvec = 0,3                      
         zdVall(nx,ny,nz,nt,nsvec,iso) = 0.D0
         DO nii = 0,1; DO ni = 0,1; DO nss = 0,1; DO ns = 0,1                              
            zdVall(nx,ny,nz,nt,nsvec,iso) = &
                 zdVall(nx,ny,nz,nt,nsvec,iso) & 
                 + zmat(ns,nss,ni,nii) &
                 *ztau2x2(ns,nss,nsvec) &
                 *ztau2x2(ni,nii,iso)                                     
         END DO; END DO; END DO; END DO
      END DO; END DO
      
   END DO; END DO; END DO
   
END DO

END SUBROUTINE dV

! **********************************************************

SUBROUTINE generalmidoverlap(zvecs,zdualvecs,zdeter,aldeterabs,nwhich1,nt1,nwhich2,nt2)
  
IMPLICIT INTEGER(i-n)
IMPLICIT DOUBLE PRECISION(a-h,o-y)
IMPLICIT COMPLEX*16(z)
INCLUDE "input.f90"

DIMENSION zvecs(0:L-1,0:L-1,0:L-1,0:Lt,0:1,0:1,0:n_f-1)
DIMENSION zdualvecs(0:L-1,0:L-1,0:L-1,0:Lt,0:1,0:1,0:n_f-1)
DIMENSION zcorrmatrix(0:n_f-1,0:n_f-1)
DIMENSION nwhich1(0:n_f-1)
DIMENSION nwhich2(0:n_f-1)
DIMENSION indx(0:n_f-1)
      
DO np1 = 0,n_f-1; DO np2 = 0,n_f-1
   
   zcorrmatrix(np2,np1) = 0.D0
   
   DO ni = 0,1; DO ns = 0,1
      DO nz = 0,L-1; DO ny = 0,L-1; DO nx = 0,L-1               
         zcorrmatrix(np2,np1) = zcorrmatrix(np2,np1) &
              + zdualvecs(nx,ny,nz,nt2,ns,ni,nwhich2(np2)) &
              * zvecs(nx,ny,nz,nt1,ns,ni,nwhich1(np1))/L**3
         
      END DO; END DO; END DO
   END DO; END DO
   
END DO; END DO

CALL zludcmp_ratio(zcorrmatrix,zdeter,aldeterabs,n_f)

END SUBROUTINE generalmidoverlap

! **********************************************************
