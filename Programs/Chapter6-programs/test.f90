!***********************************************************

PROGRAM test
  
  IMPLICIT integer(i-n)
  IMPLICIT double precision(a-h,o-y)
  IMPLICIT complex*16(z)
  
  INCLUDE 'input.f90'
      
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
  
!************************************************************

  call srand(123)
  write(*,*)rand()
  write(*,*)rand()
  write(*,*)rand()

  x1 = 1.D0
  w = 1.3D0
  write(*,*)w
  v = dlog(1.3D0)/w;
  write(*,*)v
  v = dlog(w)/w;
  write(*,*)v
  rr = x1*v
  
  stop
end PROGRAM test

!************************************************************
