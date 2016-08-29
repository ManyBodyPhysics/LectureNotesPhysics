  ! use D for double precision exponents
  
  parameter(n_f = 2)
  parameter(L = 4)
  parameter(Lt = 6)  
  parameter(cutoff = 100.D0, temporalcutoff = 150.D0)
  parameter(improveN = 2)
  parameter(improveP = 0)
  parameter(improveD = 0)
  parameter(c1S0_phys = -5.615D-5)
  parameter(c3S1_phys = -6.543D-5)
  
  ! run parameters
  
  parameter(ntot = 400000)
  parameter(myseed = 123981)
  parameter(ntherm = 1000) 
  parameter(nprintevery = 500)
  parameter(measureevery = 2)
  parameter(ncheck = 0)
  parameter(nHMC = 10)
  parameter(eHMC = 0.1D0)
  parameter(startspread = 1.0D0)
  
  ! constants 
  
  parameter(amnu_phys = 938.92D0)
  parameter(ampi3_phys = 134.98D0)
  parameter(fpi_phys = 92.2D0)
  parameter(gA = 1.29D0)
  parameter(pi = 3.14159265358979324D0)
  
  ! derived parameters
  
  parameter(atovera = cutoff/temporalcutoff)
  parameter(amnu = amnu_phys/cutoff)
  parameter(ampi3 = ampi3_phys/cutoff)
  parameter(fpi = fpi_phys/cutoff)
  parameter(c0_phys = (3.D0*c1S0_phys+1.D0*c3S1_phys)/4.D0)
  parameter(cI_phys = (1.D0*c1S0_phys-1.D0*c3S1_phys)/4.D0)
  parameter(c0 = c0_phys*cutoff**2)
  parameter(cI = cI_phys*cutoff**2)
  parameter(h = atovera/(2.D0*amnu))
  
