!***********************************************************

  if (improveN .eq. 0) then
     w0_N = 1.D0
     w1_N = 1.D0
     w2_N = 0.D0
     w3_N = 0.D0
  elseif (improveN .eq. 1) then
     w0_N = 5.D0/4.D0
     w1_N = 4.D0/3.D0
     w2_N = 1.D0/12.D0
     w3_N = 0.D0
  elseif (improveN .eq. 2) then
     w0_N = 49.D0/36.D0
     w1_N = 3.D0/2.D0
     w2_N = 3.D0/20.D0
     w3_N = 1.D0/90.D0
  else			 
     if (myid .eq. 0) then
        write(*,*)'error in improveN'
     endif
     stop
  endif
  
  if (improveP .eq. 0) then
     w0_P = 1.D0
     w1_P = 1.D0
     w2_P = 0.D0
     w3_P = 0.D0
  elseif (improveP .eq. 1) then
     w0_P = 5.D0/4.D0
     w1_P = 4.D0/3.D0
     w2_P = 1.D0/12.D0
     w3_P = 0.D0
  elseif (improveP .eq. 2) then
     w0_P = 49.D0/36.D0
     w1_P = 3.D0/2.D0
     w2_P = 3.D0/20.D0
     w3_P = 1.D0/90.D0
  else			 
     if (myid .eq. 0) then
        write(*,*)'error in improveP'
     endif
     stop
  endif
  
  if (improveD .eq. 0) then
     o1 = 1.D0
     o2 = 0.D0
     o3 = 0.D0
  elseif (improveD .eq. 1) then
     o1 = 4.D0/3.D0
     o2 = 1.D0/6.D0
     o3 = 0.D0
  elseif (improveD .eq. 2) then
     o1 = 3.D0/2.D0
     o2 = 3.D0/10.D0
     o3 = 1.D0/30.D0
  else			 
     if (myid .eq. 0) then
        write(*,*)'error in improveD'
     endif
     stop
  endif
  
  !     parameter for renormalizing the pion field
  
  qpi3 = atovera*(ampi3**2+6.D0*w0_P)
  
!***********************************************************
