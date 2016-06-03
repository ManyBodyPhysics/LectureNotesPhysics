!***********************************************************

PROGRAM nuclei
  
  IMPLICIT integer(i-n)
  IMPLICIT double precision(a-h,o-y)
  IMPLICIT complex*16(z)
  
  INCLUDE 'input.f90'
  INCLUDE 'mpif.h'
      
  DIMENSION s(0:L-1,0:L-1,0:L-1,0:Lt-1)
  DIMENSION p_s(0:L-1,0:L-1,0:L-1,0:Lt-1)
  DIMENSION snew(0:L-1,0:L-1,0:L-1,0:Lt-1)
  DIMENSION p_snew(0:L-1,0:L-1,0:L-1,0:Lt-1)
  DIMENSION sHMC(0:L-1,0:L-1,0:L-1,0:Lt-1,0:nHMC)
  DIMENSION p_sHMC(0:L-1,0:L-1,0:L-1,0:Lt-1,0:nHMC)
  
  DIMENSION sI(0:L-1,0:L-1,0:L-1,0:Lt-1,1:3)
  DIMENSION p_sI(0:L-1,0:L-1,0:L-1,0:Lt-1,1:3)
  DIMENSION sInew(0:L-1,0:L-1,0:L-1,0:Lt-1,1:3)
  DIMENSION p_sInew(0:L-1,0:L-1,0:L-1,0:Lt-1,1:3)
  DIMENSION sIHMC(0:L-1,0:L-1,0:L-1,0:Lt-1,1:3,0:nHMC)
  DIMENSION p_sIHMC(0:L-1,0:L-1,0:L-1,0:Lt-1,1:3,0:nHMC)
  
  DIMENSION pion(0:L-1,0:L-1,0:L-1,0:Lt-1,1:3)
  DIMENSION pionnew(0:L-1,0:L-1,0:L-1,0:Lt-1,1:3)
  DIMENSION p_pion(0:L-1,0:L-1,0:L-1,0:Lt-1,1:3)
  DIMENSION p_pionnew(0:L-1,0:L-1,0:L-1,0:Lt-1,1:3)
  DIMENSION pionHMC(0:L-1,0:L-1,0:L-1,0:Lt-1,1:3,0:nHMC)
  DIMENSION p_pionHMC(0:L-1,0:L-1,0:L-1,0:Lt-1,1:3,0:nHMC)
  
  DIMENSION dVds(0:L-1,0:L-1,0:L-1,0:Lt-1)
  DIMENSION dVdsI(0:L-1,0:L-1,0:L-1,0:Lt-1,1:3)
  DIMENSION dVdpion(0:L-1,0:L-1,0:L-1,0:Lt-1,1:3)
  DIMENSION zdVall(0:L-1,0:L-1,0:L-1,0:Lt-1,0:3,0:3)
  DIMENSION zvecs(0:L-1,0:L-1,0:L-1,0:Lt,0:1,0:1,0:n_f-1)
  DIMENSION zdualvecs(0:L-1,0:L-1,0:L-1,0:Lt,0:1,0:1,0:n_f-1)
  DIMENSION zwave(0:L-1,0:L-1,0:L-1,0:1,0:1,0:n_f-1)
  DIMENSION zdualwave(0:L-1,0:L-1,0:L-1,0:1,0:1,0:n_f-1)
  DIMENSION zcorrmatrix(0:n_f-1,0:n_f-1)
  DIMENSION zcorrinv(0:n_f-1,0:n_f-1)
  DIMENSION zdcorrmatrix(0:n_f-1,0:n_f-1)
  DIMENSION ztau2x2(0:1,0:1,0:3)
  EXTERNAL zludcmp
  EXTERNAL zdoinv
  
!***********************************************************

  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)
  
  IF (myid .eq. 0) THEN
     CPUtime_0 = MPI_Wtime()
  END IF
  
  !     improve = 0:  standard lattice 
  !     improve = 1:  O(a**2)-improved kinetic action
  !     improve = 2:  O(a**4)-improved kinetic action
  !     +w0 is the coefficient at the center
  !     -w1 is the coefficient for the nearest neighbor hop
  !     +w2 is the coefficient for the next-nearest neighbor hop
  !     -w3 is the coefficient for the next-next-nearest neighbor hop
  
  INCLUDE "improve.f90"
  
  !     spin/isospin conventions:
  !
  !     spin-up ---- ns = 0
  !     spin-down -- ns = 1
  !
  !     proton ----- ni = 0
  !     neutron ---- ni = 1
  
  IF (myid .eq. 0) THEN
     
     WRITE(*,*) ' '
     WRITE(*,*) '**********'
     WRITE(*,'(1X,A,I0)') 'numprocs = ',numprocs
     WRITE(*,'(1X,A,I0)') 'n_f = ',n_f
     WRITE(*,'(1X,A,I0)') 'L = ',L
     WRITE(*,'(2(1X,A,F0.6,T25))') 'cutoff = ',cutoff,'temporalcutoff = ',temporalcutoff
     WRITE(*,'(1X,A,I0,T25)') 'Lt = ',Lt
     WRITE(*,'(1X,A,F0.6,1X,A)') '1/a = ',cutoff,'MeV'
     WRITE(*,'(1X,A,F0.6,1X,A)') '1/a_t = ',temporalcutoff,'MeV'
     WRITE(*,'(1X,A,F0.6)') 'atovera = ',atovera
     WRITE(*,'(1X,A,I0)') 'improveN = ',improveN
     WRITE(*,'(1X,A,I0)') 'improveP = ',improveP
     WRITE(*,'(1X,A,I0)') 'improveD = ',improveD
     
     WRITE(*,*) ' '
     WRITE(*,*) '**********'
     WRITE(*,'(1X,A,I0)') 'ntot = ',ntot
     WRITE(*,'(1X,A,I0)') 'myseed = ',myseed
     WRITE(*,'(1X,A,I0)') 'ntherm = ',ntherm
     WRITE(*,'(1X,A,I0)') 'nprintevery = ',nprintevery
     WRITE(*,'(1X,A,I0)') 'measureevery = ',measureevery
     WRITE(*,'(1X,A,I0)') 'ncheck = ',ncheck
     WRITE(*,'(1X,A,I0)') 'nHMC = ',nHMC
     WRITE(*,'(1X,A,E14.6)') 'eHMC =',eHMC
     WRITE(*,'(1X,A,E14.6)') 'epsilon =',epsilon
     WRITE(*,'(1X,A,E14.6)') 'epsprime =',epsprime
     WRITE(*,'(1X,A,F0.6)') 'startspread = ',startspread
     
     WRITE(*,*) ' '
     WRITE(*,*) '**********'
     WRITE(*,'(1X,A,E14.6,1X,A)') 'c1S0 =',c1S0_phys,'MeV^(-2)'
     WRITE(*,'(1X,A,E14.6,1X,A)') 'c3S1 =',c3S1_phys,'MeV^(-2)'
     WRITE(*,*) ' '
     WRITE(*,'(1X,A,E14.6,1X,A)') 'c0 =',c0_phys,'MeV^(-2)'
     WRITE(*,'(1X,A,E14.6,1X,A)') 'cI =',cI_phys,'MeV^(-2)'
     
     WRITE(*,*) ' '
     WRITE(*,*) '**********'
     WRITE(*,'(1X,A,F0.6)') 'gA = ',gA
     WRITE(*,'(1X,A,E14.6)') 'h =',h
     WRITE(*,'(1X,A,F0.6,1X,A)') 'amnu = ',amnu_phys,'MeV'
     WRITE(*,'(1X,A,F0.6,1X,A)') 'ampi3 = ',ampi3_phys,'MeV'
     WRITE(*,'(1X,A,F0.6,1X,A)') 'fpi = ',fpi_phys,'MeV'
     WRITE(*,'(1X,A,F0.6)') 'qpi3 = ',qpi3
     WRITE(*,'(1X,A,I0)') 'nHMC = ',nHMC
     WRITE(*,'(1X,A,E14.6)') 'eHMC =',eHMC
     WRITE(*,*)
  END IF

!************************************************************
  
  CALL sgrnd(myseed+10*myid)
  
  accept = 0.D0
  zdeterphasebin = 0.D0
  
  DO ndir = 0,3
     DO nii = 0,1; DO ni = 0,1
        ztau2x2(ni,nii,ndir) = 0.D0
     END DO; END DO
  END DO
  
  !     ndir = 0,1,2,3 -> 1,x,y,z
  
  ztau2x2(0,0,0) = 1.D0
  ztau2x2(1,1,0) = 1.D0
  ztau2x2(0,1,1) = 1.D0
  ztau2x2(1,0,1) = 1.D0
  ztau2x2(0,1,2) = (0.D0,-1.D0)
  ztau2x2(1,0,2) = (0.D0,1.D0)
  ztau2x2(0,0,3) = 1.D0
  ztau2x2(1,1,3) = -1.D0
  
  DO npart = 0,n_f-1         
     DO ni = 0,1; DO ns = 0,1
        DO nz = 0,L-1; DO ny = 0,L-1; DO nx = 0,L-1
           zvecs(nx,ny,nz,0,ns,ni,npart) = 0.D0
           zdualvecs(nx,ny,nz,Lt,ns,ni,npart) = 0.D0
        END DO; END DO; END DO
     END DO; END DO
  END DO
  
  CALL waveinit(zwave,zdualwave,myid)
  
!************************************************************
  
  DO nt = 0,Lt-1
     DO nz = 0,L-1; DO ny = 0,L-1; DO nx = 0,L-1
        CALL gaussrnd(gr)
        s(nx,ny,nz,nt) = gr*dsqrt(startspread) & 
             + startspread
        DO iso = 1,3
           CALL gaussrnd(gr)
           sI(nx,ny,nz,nt,iso) = gr*dsqrt(startspread)
        END DO
        DO iso = 1,3
           CALL gaussrnd(gr)
           pion(nx,ny,nz,nt,iso) = gr*dsqrt(startspread)
        END DO
     END DO; END DO; END DO
  END DO
  
  DO ntrial = 1,ntot
     
     nconfig = ntrial-ntherm
     mconfig = nconfig/measureevery
     
     DO nt = 0,Lt-1
        DO nz = 0,L-1; DO ny = 0,L-1; DO nx = 0,L-1
           CALL gaussrnd(gr)
           p_s(nx,ny,nz,nt) = gr
           DO iso = 1,3
              CALL gaussrnd(gr)
              p_sI(nx,ny,nz,nt,iso) = gr
           END DO
           DO iso = 1,3
              CALL gaussrnd(gr)
              p_pion(nx,ny,nz,nt,iso) = gr
           END DO
        END DO; END DO; END DO
     END DO
     
     bose = 0.D0
     DO nt = 0,Lt-1
        DO nz = 0,L-1; DO ny = 0,L-1; DO nx = 0,L-1      
           bose = bose &
                + s(nx,ny,nz,nt)**2.D0/2.D0 &
                + p_s(nx,ny,nz,nt)**2.D0/2.D0
           DO iso = 1,3
              bose = bose &
                   + sI(nx,ny,nz,nt,iso)**2.D0/2.D0 &
                   + p_sI(nx,ny,nz,nt,iso)**2.D0/2.D0 &
                   + pion(nx,ny,nz,nt,iso)**2.D0/2.D0 &
                   + atovera/qpi3*pion(nx,ny,nz,nt,iso)*( &
                   - w1_P*pion(MOD(nx+1,L),ny,nz,nt,iso) &
                   - w1_P*pion(nx,MOD(ny+1,L),nz,nt,iso) &
                   - w1_P*pion(nx,ny,MOD(nz+1,L),nt,iso) &
                   + w2_P*pion(MOD(nx+2,L),ny,nz,nt,iso) &
                   + w2_P*pion(nx,MOD(ny+2,L),nz,nt,iso) &
                   + w2_P*pion(nx,ny,MOD(nz+2,L),nt,iso) &
                   - w3_P*pion(MOD(nx+3,L),ny,nz,nt,iso) &
                   - w3_P*pion(nx,MOD(ny+3,L),nz,nt,iso) &
                   - w3_P*pion(nx,ny,MOD(nz+3,L),nt,iso)) &
                   + p_pion(nx,ny,nz,nt,iso)**2.D0/2.D0
           END DO
        END DO; END DO; END DO
     END DO
     
     DO nt = 0,Lt-1
        DO nz = 0,L-1; DO ny = 0,L-1; DO nx = 0,L-1      
           sHMC(nx,ny,nz,nt,0) = s(nx,ny,nz,nt)
           DO iso = 1,3
              sIHMC(nx,ny,nz,nt,iso,0) = & 
                   sI(nx,ny,nz,nt,iso)
              pionHMC(nx,ny,nz,nt,iso,0) = &
                   pion(nx,ny,nz,nt,iso)
           END DO
        END DO; END DO; END DO
     END DO
     
     !     initial half step for p_sHMC, p_sIHMC, p_pionHMC
     
     CALL getzvecs(s,sI,zvecs,zwave,Lt,0, &
          pion,ztau2x2,n_f)            
     CALL getzdualvecs(s,sI,zdualvecs,zdualwave, &
          Lt,0,pion,ztau2x2,n_f)         
     CALL getinvcorr(zvecs,zdualvecs,zldeter, &
          zcorrmatrix,zcorrinv,Lt)
     aldeterabs = DBLE(zldeter)
     zdeterphase = CDEXP((0.D0,1.D0)*DIMAG(zldeter))
     act = bose - aldeterabs
     
     CALL dV(zvecs,zdualvecs,ztau2x2,zcorrinv,zdVall)
     
     DO nt = 0,Lt-1        
        DO nz = 0,L-1; DO ny = 0,L-1; DO nx = 0,L-1      
           
           DO npart1 = 0,n_f-1; DO npart2 = 0,n_f-1
              zdcorrmatrix(npart2,npart1) = 0.D0
              DO ni = 0,1; DO ns = 0,1
                 zdcorrmatrix(npart2,npart1) = &
                      zdcorrmatrix(npart2,npart1) + &
                      zdualvecs(nx,ny,nz,nt+1,ns,ni,npart2) &
                      *zvecs(nx,ny,nz,nt,ns,ni,npart1) &
                      *CDSQRT(-c0*atovera*(1.D0,0.D0))/L**3
              END DO; END DO
           END DO; END DO
           
           dVds(nx,ny,nz,nt) = s(nx,ny,nz,nt) 
           DO npart1 = 0,n_f-1; DO npart2 = 0,n_f-1
              dVds(nx,ny,nz,nt) = dVds(nx,ny,nz,nt) &
                   - DBLE(zdcorrmatrix(npart2,npart1) &
                   *zcorrinv(npart1,npart2))
           END DO; END DO
                 
           p_sHMC(nx,ny,nz,nt,0) = &
                p_s(nx,ny,nz,nt) - 0.5D0*eHMC*dVds(nx,ny,nz,nt)
           
           DO iso = 1,3
              dVdsI(nx,ny,nz,nt,iso) = &
                   sI(nx,ny,nz,nt,iso) &
                   - DBLE((0.D0,1.D0)*CDSQRT(cI*atovera*(1.D0,0.D0)) &
                   *zdVall(nx,ny,nz,nt,0,iso)/L**3)
           END DO
           
           DO iso = 1,3
              
              dVdpion(nx,ny,nz,nt,iso) = &
                   pion(nx,ny,nz,nt,iso) &
                   + atovera/qpi3*( &
                   - w1_P*pion(MOD(nx+1,L),ny,nz,nt,iso) &
                   - w1_P*pion(MOD(nx-1+L,L),ny,nz,nt,iso) &
                   - w1_P*pion(nx,MOD(ny+1,L),nz,nt,iso) &
                   - w1_P*pion(nx,MOD(ny-1+L,L),nz,nt,iso) &
                   - w1_P*pion(nx,ny,MOD(nz+1,L),nt,iso) &
                   - w1_P*pion(nx,ny,MOD(nz-1+L,L),nt,iso) &
                   + w2_P*pion(MOD(nx+2,L),ny,nz,nt,iso) &
                   + w2_P*pion(MOD(nx-2+L,L),ny,nz,nt,iso) &
                   + w2_P*pion(nx,MOD(ny+2,L),nz,nt,iso) &
                   + w2_P*pion(nx,MOD(ny-2+L,L),nz,nt,iso) &
                   + w2_P*pion(nx,ny,MOD(nz+2,L),nt,iso) &
                   + w2_P*pion(nx,ny,MOD(nz-2+L,L),nt,iso) &
                   - w3_P*pion(MOD(nx+3,L),ny,nz,nt,iso) &
                   - w3_P*pion(MOD(nx-3+L,L),ny,nz,nt,iso) &
                   - w3_P*pion(nx,MOD(ny+3,L),nz,nt,iso) &
                   - w3_P*pion(nx,MOD(ny-3+L,L),nz,nt,iso) &
                   - w3_P*pion(nx,ny,MOD(nz+3,L),nt,iso) &
                   - w3_P*pion(nx,ny,MOD(nz-3+L,L),nt,iso))
              
              dVdpion(nx,ny,nz,nt,iso) = &
                   dVdpion(nx,ny,nz,nt,iso) &
                   + gA*atovera &
                   /(2.D0*fpi*dsqrt(qpi3)*L**3) &
                   *DBLE(0.D0 &
                   + o1/2.D0*zdVall(MOD(nx-1+L,L),ny,nz,nt,1,iso) &
                   + o1/2.D0*zdVall(nx,MOD(ny-1+L,L),nz,nt,2,iso) &
                   + o1/2.D0*zdVall(nx,ny,MOD(nz-1+L,L),nt,3,iso) &
                   - o1/2.D0*zdVall(MOD(nx+1,L),ny,nz,nt,1,iso) &
                   - o1/2.D0*zdVall(nx,MOD(ny+1,L),nz,nt,2,iso) &
                   - o1/2.D0*zdVall(nx,ny,MOD(nz+1,L),nt,3,iso) &
                   - o2/2.D0*zdVall(MOD(nx-2+L,L),ny,nz,nt,1,iso) &
                   - o2/2.D0*zdVall(nx,MOD(ny-2+L,L),nz,nt,2,iso) &
                   - o2/2.D0*zdVall(nx,ny,MOD(nz-2+L,L),nt,3,iso) &
                   + o2/2.D0*zdVall(MOD(nx+2,L),ny,nz,nt,1,iso) &
                   + o2/2.D0*zdVall(nx,MOD(ny+2,L),nz,nt,2,iso) &
                   + o2/2.D0*zdVall(nx,ny,MOD(nz+2,L),nt,3,iso) &
                   + o3/2.D0*zdVall(MOD(nx-3+L,L),ny,nz,nt,1,iso) &
                   + o3/2.D0*zdVall(nx,MOD(ny-3+L,L),nz,nt,2,iso) &
                   + o3/2.D0*zdVall(nx,ny,MOD(nz-3+L,L),nt,3,iso) &
                   - o3/2.D0*zdVall(MOD(nx+3,L),ny,nz,nt,1,iso) &
                   - o3/2.D0*zdVall(nx,MOD(ny+3,L),nz,nt,2,iso) &
                   - o3/2.D0*zdVall(nx,ny,MOD(nz+3,L),nt,3,iso))
              
           END DO
           
           DO iso = 1,3
              p_sIHMC(nx,ny,nz,nt,iso,0) = &
                   p_sI(nx,ny,nz,nt,iso) &
                   - 0.5D0*eHMC*dVdsI(nx,ny,nz,nt,iso)
              p_pionHMC(nx,ny,nz,nt,iso,0) = &
                   p_pion(nx,ny,nz,nt,iso) &
                   - 0.5D0*eHMC*dVdpion(nx,ny,nz,nt,iso)
           END DO
           
        END DO; END DO; END DO
     END DO
     
!*************************************************************
     
     ! full steps for sHMC, sIHMC, pionHMC and then 
     ! p_sHMC, p_sIHMC, p_pionHMC in leapfrog succession
     ! last step for p_sHMC, p_sIHMC, p_pionHMC should be a half step
     
     DO nstep = 0,nHMC-1
        
        IF (nstep .eq. nHMC-1) THEN
           step = 0.5D0
        ELSE
           step = 1.D0
        END IF
        
        DO nt = 0,Lt-1
           DO nz = 0,L-1; DO ny = 0,L-1; DO nx = 0,L-1      
              sHMC(nx,ny,nz,nt,nstep+1) = sHMC(nx,ny,nz,nt,nstep) &
                   + eHMC*p_sHMC(nx,ny,nz,nt,nstep)
           END DO; END DO; END DO
        END DO
        DO nt = 0,Lt-1
           DO nz = 0,L-1; DO ny = 0,L-1; DO nx = 0,L-1                  
              DO iso = 1,3
                 sIHMC(nx,ny,nz,nt,iso,nstep+1) = &
                      sIHMC(nx,ny,nz,nt,iso,nstep) &
                      + eHMC*p_sIHMC(nx,ny,nz,nt,iso,nstep)
                 pionHMC(nx,ny,nz,nt,iso,nstep+1) = &
                      pionHMC(nx,ny,nz,nt,iso,nstep) &
                      + eHMC*p_pionHMC(nx,ny,nz,nt,iso,nstep)
              END DO
           END DO; END DO; END DO
        END DO
        
        CALL getzvecs(sHMC(0,0,0,0,nstep+1), &
             sIHMC(0,0,0,Ltouter,1,nstep+1), &
             zvecs,zwave,Lt,0, &
             pionHMC(0,0,0,Ltouter,1,nstep+1), &
             ztau2x2,n_f)
        
        CALL getzdualvecs(sHMC(0,0,0,0,nstep+1), &
             sIHMC(0,0,0,Ltouter,1,nstep+1), &
             zdualvecs,zdualwave,Lt,0, &
             pionHMC(0,0,0,Ltouter,1,nstep+1), &
             ztau2x2,n_f)
        
        CALL getinvcorr(zvecs,zdualvecs,zldeter_HMC, &
             zcorrmatrix,zcorrinv,Lt-1)
        
        CALL dV(zvecs,zdualvecs,ztau2x2, &
             zcorrinv,zdVall)

        DO nt = 0,Lt-1                    
           DO nz = 0,L-1; DO ny = 0,L-1; DO nx = 0,L-1      
              
              DO npart1 = 0,n_f-1; DO npart2 = 0,n_f-1
                 zdcorrmatrix(npart2,npart1) = 0.D0
                 DO ni = 0,1; DO ns = 0,1
                    zdcorrmatrix(npart2,npart1) = &
                         zdcorrmatrix(npart2,npart1) &
                         + zdualvecs(nx,ny,nz,nt+1,ns,ni,npart2) &
                         *zvecs(nx,ny,nz,nt,ns,ni,npart1) &
                         *CDSQRT(-c0*atovera*(1.D0,0.D0))/L**3
                 END DO; END DO; END DO
              END DO
                    
              dVds(nx,ny,nz,nt) = sHMC(nx,ny,nz,nt,nstep+1) 
              DO npart1 = 0,n_f-1; DO npart2 = 0,n_f-1
                 dVds(nx,ny,nz,nt) = dVds(nx,ny,nz,nt) & 
                      - DBLE(zdcorrmatrix(npart2,npart1) &
                      *zcorrinv(npart1,npart2))
              END DO; END DO
              
              p_sHMC(nx,ny,nz,nt,nstep+1) = &
                   p_sHMC(nx,ny,nz,nt,nstep) &
                   - step*eHMC*dVds(nx,ny,nz,nt)
              
              DO iso = 1,3                           
                 dVdsI(nx,ny,nz,nt,iso) = &
                      sIHMC(nx,ny,nz,nt,iso,nstep+1) &
                      - DBLE((0.D0,1.D0)*CDSQRT(cI*atovera*(1.D0,0.D0)) &
                      *zdVall(nx,ny,nz,nt,0,iso)/L**3)
              END DO
              
              DO iso = 1,3
                       
                 dVdpion(nx,ny,nz,nt,iso) = &
                      pionHMC(nx,ny,nz,nt,iso,nstep+1) &
                      + atovera/qpi3*( &
                      - w1_P*pionHMC(MOD(nx+1,L),ny,nz,nt,iso,nstep+1) &
                      - w1_P*pionHMC(MOD(nx-1+L,L),ny,nz,nt,iso,nstep+1) &
                      - w1_P*pionHMC(nx,MOD(ny+1,L),nz,nt,iso,nstep+1) &
                      - w1_P*pionHMC(nx,MOD(ny-1+L,L),nz,nt,iso,nstep+1) &
                      - w1_P*pionHMC(nx,ny,MOD(nz+1,L),nt,iso,nstep+1) &
                      - w1_P*pionHMC(nx,ny,MOD(nz-1+L,L),nt,iso,nstep+1) &
                      + w2_P*pionHMC(MOD(nx+2,L),ny,nz,nt,iso,nstep+1) &
                      + w2_P*pionHMC(MOD(nx-2+L,L),ny,nz,nt,iso,nstep+1) &
                      + w2_P*pionHMC(nx,MOD(ny+2,L),nz,nt,iso,nstep+1) &
                      + w2_P*pionHMC(nx,MOD(ny-2+L,L),nz,nt,iso,nstep+1) &
                      + w2_P*pionHMC(nx,ny,MOD(nz+2,L),nt,iso,nstep+1) &
                      + w2_P*pionHMC(nx,ny,MOD(nz-2+L,L),nt,iso,nstep+1) &
                      - w3_P*pionHMC(MOD(nx+3,L),ny,nz,nt,iso,nstep+1) &
                      - w3_P*pionHMC(MOD(nx-3+L,L),ny,nz,nt,iso,nstep+1) &
                      - w3_P*pionHMC(nx,MOD(ny+3,L),nz,nt,iso,nstep+1) &
                      - w3_P*pionHMC(nx,MOD(ny-3+L,L),nz,nt,iso,nstep+1) &
                      - w3_P*pionHMC(nx,ny,MOD(nz+3,L),nt,iso,nstep+1) &
                      - w3_P*pionHMC(nx,ny,MOD(nz-3+L,L),nt,iso,nstep+1))
                 
                 dVdpion(nx,ny,nz,nt,iso) = &
                      dVdpion(nx,ny,nz,nt,iso) &
                      + gA*atovera/(2.D0*fpi*dsqrt(qpi3)*L**3) &
                      *DBLE( &
                      + o1/2.D0*zdVall(MOD(nx-1+L,L),ny,nz,nt,1,iso) &
                      + o1/2.D0*zdVall(nx,MOD(ny-1+L,L),nz,nt,2,iso) &
                      + o1/2.D0*zdVall(nx,ny,MOD(nz-1+L,L),nt,3,iso) &
                      - o1/2.D0*zdVall(MOD(nx+1,L),ny,nz,nt,1,iso) &
                      - o1/2.D0*zdVall(nx,MOD(ny+1,L),nz,nt,2,iso) &
                      - o1/2.D0*zdVall(nx,ny,MOD(nz+1,L),nt,3,iso) &
                      - o2/2.D0*zdVall(MOD(nx-2+L,L),ny,nz,nt,1,iso) &
                      - o2/2.D0*zdVall(nx,MOD(ny-2+L,L),nz,nt,2,iso) &
                      - o2/2.D0*zdVall(nx,ny,MOD(nz-2+L,L),nt,3,iso) &
                      + o2/2.D0*zdVall(MOD(nx+2,L),ny,nz,nt,1,iso) &
                      + o2/2.D0*zdVall(nx,MOD(ny+2,L),nz,nt,2,iso) &
                      + o2/2.D0*zdVall(nx,ny,MOD(nz+2,L),nt,3,iso) &
                      + o3/2.D0*zdVall(MOD(nx-3+L,L),ny,nz,nt,1,iso) &
                      + o3/2.D0*zdVall(nx,MOD(ny-3+L,L),nz,nt,2,iso) &
                      + o3/2.D0*zdVall(nx,ny,MOD(nz-3+L,L),nt,3,iso) &
                      - o3/2.D0*zdVall(MOD(nx+3,L),ny,nz,nt,1,iso) &
                      - o3/2.D0*zdVall(nx,MOD(ny+3,L),nz,nt,2,iso) &
                      - o3/2.D0*zdVall(nx,ny,MOD(nz+3,L),nt,3,iso))
                 
              END DO
              
              DO iso = 1,3
                 p_sIHMC(nx,ny,nz,nt,iso,nstep+1) = &
                      p_sIHMC(nx,ny,nz,nt,iso,nstep) &
                      - step*eHMC*dVdsI(nx,ny,nz,nt,iso)
                 p_pionHMC(nx,ny,nz,nt,iso,nstep+1) = &
                      p_pionHMC(nx,ny,nz,nt,iso,nstep) &
                      - step*eHMC*dVdpion(nx,ny,nz,nt,iso)
              END DO

           END DO
        END DO; END DO; END DO

     END DO
         
!***********************************************************
         
         DO nt = 0,Lt-1
            DO nz = 0,L-1; DO ny = 0,L-1; DO nx = 0,L-1                           
               snew(nx,ny,nz,nt) = sHMC(nx,ny,nz,nt,nHMC)
               p_snew(nx,ny,nz,nt) = p_sHMC(nx,ny,nz,nt,nHMC)
               DO iso = 1,3
                  sInew(nx,ny,nz,nt,iso) = sIHMC(nx,ny,nz,nt,iso,nHMC)
                  p_sInew(nx,ny,nz,nt,iso) = p_sIHMC(nx,ny,nz,nt,iso,nHMC)
                  pionnew(nx,ny,nz,nt,iso) = pionHMC(nx,ny,nz,nt,iso,nHMC)
                  p_pionnew(nx,ny,nz,nt,iso) = p_pionHMC(nx,ny,nz,nt,iso,nHMC)
               END DO
            END DO; END DO; END DO
         END DO
         
         bosenew = 0.D0
         DO nt = 0,Lt-1
            DO nz = 0,L-1; DO ny = 0,L-1; DO nx = 0,L-1      
               bosenew = bosenew &
                    + snew(nx,ny,nz,nt)**2.D0/2.D0 &
                    + p_snew(nx,ny,nz,nt)**2.D0/2.D0
            END DO; END DO; END DO
         END DO
         
         DO nt = 0,Lt-1
            DO nz = 0,L-1; DO ny = 0,L-1; DO nx = 0,L-1      
               DO iso = 1,3
                  bosenew = bosenew &
                       + sInew(nx,ny,nz,nt,iso)**2.D0/2.D0 &
                       + p_sInew(nx,ny,nz,nt,iso)**2.D0/2.D0 &
                       + pionnew(nx,ny,nz,nt,iso)**2.D0/2.D0 &
                       + atovera/qpi3*pionnew(nx,ny,nz,nt,iso)*( &
                       - w1_P*pionnew(MOD(nx+1,L),ny,nz,nt,iso) &
                       - w1_P*pionnew(nx,MOD(ny+1,L),nz,nt,iso) &
                       - w1_P*pionnew(nx,ny,MOD(nz+1,L),nt,iso) &
                       + w2_P*pionnew(MOD(nx+2,L),ny,nz,nt,iso) &
                       + w2_P*pionnew(nx,MOD(ny+2,L),nz,nt,iso) &
                       + w2_P*pionnew(nx,ny,MOD(nz+2,L),nt,iso) &
                       - w3_P*pionnew(MOD(nx+3,L),ny,nz,nt,iso) &
                       - w3_P*pionnew(nx,MOD(ny+3,L),nz,nt,iso) &
                       - w3_P*pionnew(nx,ny,MOD(nz+3,L),nt,iso)) &
                             + p_pionnew(nx,ny,nz,nt,iso)**2.D0/2.D0
               END DO
            END DO; END DO; END DO
         END DO
         
         CALL getzvecs(snew,sInew,zvecs, &
              zwave,Lt,0,pionnew,ztau2x2,n_f)
         CALL getzdualvecs(snew,sInew, &
              zdualvecs,zdualwave,Lt,0,pionnew, &
              ztau2x2,n_f)
         
         CALL getinvcorr(zvecs,zdualvecs,zldeternew, &
              zcorrmatrix,zcorrinv,Lt-1)
         aldeternewabs = DBLE(zldeternew)
         zdeternewphase = CDEXP((0.D0,1.D0)*DIMAG(zldeternew))
         actnew = bosenew - aldeternewabs
         
         IF (ncheck .eq. 1) THEN
            WRITE(*,*)'act = ',act,' actnew = ',actnew
         END IF
         
         IF (ntrial .eq. 1 .or. grnd() .lt. DEXP(-actnew+act)) THEN
            
            accept = accept + 1.
            
            DO nt = 0,Lt-1
               DO nz = 0,L-1; DO ny = 0,L-1; DO nx = 0,L-1      
                  s(nx,ny,nz,nt) = snew(nx,ny,nz,nt)
               END DO; END DO; END DO
            END DO
            DO nt = 0,Lt-1
               DO nz = 0,L-1; DO ny = 0,L-1; DO nx = 0,L-1      
                  DO iso = 1,3
                     sI(nx,ny,nz,nt,iso) = sInew(nx,ny,nz,nt,iso)
                     pion(nx,ny,nz,nt,iso) = pionnew(nx,ny,nz,nt,iso)
                  END DO
               END DO; END DO; END DO
            END DO
            aldeterabs = aldeternewabs
            zdeterphase = zdeternewphase
            
         END IF
         
         IF (nconfig .ge. 1 .and. MOD(nconfig,measureevery) .eq. 0) THEN
            
            zdeterphasebin = zdeterphasebin + zdeterphase
            CALL getzvecs(s,sI,zvecs,zwave,Lt,0,pion,ztau2x2,n_f)
            CALL getzdualvecs(s,sI,zdualvecs,zdualwave,Lt,0,pion,ztau2x2,n_f)

            deterabs = DEXP(aldeterabs)

            zoverlap = 0.D0
            DO nt = 0,Lt-1
               CALL getmidoverlap(zvecs,zdualvecs,ztemp,aldeterabs,nt,nt+1)
               zoverlap = zoverlap + ztemp/Lt
            END DO 
           
            zamp = zoverlap
            zampbin = zampbin + zamp
            ampbin = DBLE(zampbin)
            
         END IF
         
         IF (MOD(nconfig,nprintevery) .eq. 0 .or. nreadonly .eq. 1) THEN
            
            IF (myid .eq. 0) THEN
               WRITE(*,*)
               IF (nreadonly .eq. 0) THEN
                  WRITE(*,'(A,I0)')'** nconfig ',nconfig
                  WRITE(*,'(A,I0,2X,A,F9.6)')'** mconfig ',mconfig, &
                       'acceptance',accept/ntrial
               ELSE
                  WRITE(*,'(A,I0,2X,A,I0)')'** nconfig ',nconfig, &
                      '** mconfig ',mconfig
               END IF
            END IF
            
            IF (mconfig .ge. 1) THEN
               
               realpartbin = DBLE(zdeterphasebin)
               aimagpartbin = DIMAG(zdeterphasebin)

               IF (myid .eq. 0) WRITE(*,*)' '
               IF (myid .eq. 0) WRITE(*,*)'**********'

               CALL ave_err(realpartbin,average,error,numprocs,myid)
               IF (myid .eq. 0) THEN
                  realpart_ave = average/mconfig
                  realpart_err = error/mconfig
                  WRITE(*,'(1X,A30,E14.6,E14.6)') &
                      'real part of phase', &
                      realpart_ave,realpart_err
               END IF

               CALL ave_err(aimagpartbin,average,error,numprocs,myid) 
               IF (myid .eq. 0) THEN
                  aimagpart_ave = average/mconfig
                  aimagpart_err = error/mconfig
                  WRITE(*,'(1X,A30,E14.6,E14.6)') &
                      'imaginary part of phase', &
                      aimagpart_ave,aimagpart_err
               END IF
                             
               CALL ave_err(ampbin,average,error,numprocs,myid)
               IF (myid .eq. 0) THEN
                  raw_ave = average/mconfig
                  raw_err = error/mconfig
                  energy_ave = cutoff*dlog(raw_ave/realpart_ave) &
                       /atovera
                  rel_err = dsqrt((raw_err/raw_ave)**2.D0 &
                       + (realpart_err/realpart_ave)**2.D0)
                  energy_err = cutoff*dlog(raw_ave/realpart_ave &
                       *(1.D0 + rel_err))/atovera &
                       - energy_ave                  
                  WRITE(*,*)'**********'
                  WRITE(*,'(1X,A,I0,1X)') 'Lt = ',Lt
                  WRITE(*,'(1X,A30,E14.6,E14.6)') &
                       'energy (MeV)',energy_ave,energy_err
                  WRITE(*,'(1X,A30,E14.6,E14.6)') &
                       'raw amplitude',raw_ave,raw_err
                  WRITE(*,*)''
               END IF
               
            ELSE
               
               IF (myid .eq. 0) THEN               
                  WRITE(*,'(3X,A,I0,2X,A,I0)') 'nconfig = ',nconfig, &
                      'mconfig = ',mconfig
               END IF

            END IF

         END IF 

         IF (myid .eq. 0) THEN
            CPUtime_1 = MPI_Wtime()
            CPUtime = CPUtime_1 - CPUtime_0
            IF (MOD(ntrial,100) .eq. 0) THEN
               WRITE(*,*) ' '
               WRITE(*,'(3X,A,I0,2X,A,F0.2)') 'ntrial ',ntrial,'CPUtime ',CPUtime
            END IF
         END IF
        
      END DO

!************************************************************
      
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      CALL MPI_FINALIZE(ierr)

      IF (myid .eq. 0) THEN
         WRITE(*,*)'done'
      END IF

      STOP
      END PROGRAM nuclei

!************************************************************
