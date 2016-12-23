!
!
!
!   "Green's function and coupled cluster calculations of spectral function
!  for a pairing model"
!
!
!   This code is based on the lectures of the TALENT school, course #2, held at
!  GANIL on July 6th-24th, 2015 and it has been later improved to perform the
!  benchmark calcuations presented in chapters 8 and 11 of the Lecture Notes
!  in Physics
!
!
!   This version uses the SCGF and CC methods to solve the toy pairing model
!  with 4 levels and 4 spin-1/2 particles that was introduced during the first
!  week of the Talent lectures (see
!   http://nucleartalent.github.io/Course2ManyBodyMethods/doc/web/course.html )
!  and in chapters 8 and 11 of the above book. The present version on the code
!  is publishes under the GNU general public licence and it is available at:
!  https://github.com/ManyBodyPhysics/LectureNotesPhysics/blob/master/doc/src/Chapter11-programs/Pair_Model
!
!   If you are using this software (or any part of it) for publication, we ask
!  you that its origin is formally acknowledged by citing the following chapter:
!
!   C. Barbieri and A. Carbone, "Self-consistent Green's function approaches",
!  chapter 11 of Lecture Notes in Physics "An advanced course in computational
!  nuclear physics: Bridging the scales from quarks to neutron stars",
!  Edited by Hjorth-Jensen, M. P. Lombardo, U. van Kolck, Editors
!  ISBN:
!  http://arxiv.org/abs/1611.03923
!
!
!
!
!   Calculations can be done at different levels of SCGF theory by using the ADC(2)
!  approach at second order, by adding all order resummations of the two-particle
!  and two-hole ladders [the so-called 2p1h-TDA, which is an extension of ADC(2)],
!  in the full ADC(3) approximation and finally by adding CCD corrections to
!  obtain the ADC(3)-D scheme.
!
!   It also iterates the energy independent self-energy (i.e. the 'correlated HF'
!  diagram with the dressed propagator to achieve self-consistency at the sc0
!  level.
!
!   The coupled cluster equations can be solved by either including on only pp
!  and hh ladders or in the full CCD approximation scheme.
!
!   Note that for this particular model does not allow for ph (ring) TDA-like
!  resummation in any of the methods.
!
!
!   The ground state is then calculated in FCI to obtain the exact fundamental
!  energy and correlation energy, for comparison with the Koltun Sr result. It
!  is up to you, as an exercise, to do FCI for the N=3 and N=5 particle systems
!  and to extract the exact spectral function to compare with the results from
!  the ADC(n) method.
!
!   The coupled cluster equations are also solved to resum the pp and hh ladders.
!  These diagrams are a subset of the full CCD method that could be compared to the
!  correlations added in the 2p1h-ADC(2).
!
!   The particular pairing Hamiltonian used in this exercise allows us to use
!  the following two simplifications:
!
!   - The spin-up and spin-down parts of the Dyson equation completely decouple
!    and are equal because of symmetry. Hence, we need to calculate only one
!    of the two.
!
!   - The particle-hole contribution from the pairing interactions is zero in
!    the ring summation. So, we don't need to worry about rings.
!
!
!
!
!  Compile as:    gfortran Pairing_GF_CC.f90  -llapack -o Pairing.exe
!
!
!   You will need to link to the LAPACK library to compile this program.
!
!
!  (c) C. Barbieri and A. Carbone,   Surrey/Darmstadt,   October 2016.
!
!


module Bases
  !
  !  Global array to store the tables for the single-particle basis (model
  ! space) and the (quantum numbers of the) intermediate state configurations
  !

  implicit none


  ! Max dimension to be allocated for the mod. sp.
  integer, parameter :: MaxVectsDim = 16 ! Well redundant...

  integer ::   N_1b_mdsp  !  Actual sp basis size (all levels including spin)
  integer ::   isp_lv(MaxVectsDim),  isp_sp(MaxVectsDim),  isp_oc(MaxVectsDim)
  real*8  ::   xsp_en(MaxVectsDim) ! sp energies for the unperturbed Hamiltoniann

  real*8  ::   esp0_en(MaxVectsDim) ! unperturbed sp energies
  real*8  ::   Utld_en(MaxVectsDim) ! unperturbed sp energies
  real*8  ::   vhf_en(MaxVectsDim)  ! mtx els of HF potential
  real*8  ::   vhf_3rd(MaxVectsDim) ! mtx els 3rd order tadpole corrections


  integer :: N_1b_bas  !  Actual mod. sp. size (only # of levels)
  integer :: i_1b_bas(MaxVectsDim),i_1b_occ(MaxVectsDim)

  integer :: N_2p1h_bas, N_2h1p_bas  !  ADC(3) ISCs
  integer ::  innk_n1(MaxVectsDim), innk_n2(MaxVectsDim), innk_k3(MaxVectsDim)
  integer ::  ikkn_k1(MaxVectsDim), ikkn_k2(MaxVectsDim), ikkn_n3(MaxVectsDim)
  real*8  ::  xnnk_en(MaxVectsDim), xkkn_en(MaxVectsDim)

  integer :: N_2p_bas, N_2h_bas  !  partial 2p and 2h configurations
  integer ::  i_2p_n1(MaxVectsDim), i_2h_k1(MaxVectsDim)
  integer ::  i_2p_n2(MaxVectsDim), i_2h_k2(MaxVectsDim)
  real*8  ::  x_2p_en(MaxVectsDim), x_2h_en(MaxVectsDim)

  !
  !  We will use the following to store the poles and overlap
  ! amplitudes of the propagator
  integer, parameter :: MaxDySols = 20 ! Well redundant...
  integer :: N_bk_poles,  N_fw_poles
  real*8 ::  E_bk_poles(MaxVectsDim),            E_fw_poles(MaxVectsDim)
  real*8 :: SF_bk_poles(MaxVectsDim),           SF_fw_poles(MaxVectsDim)
  real*8 ::  Y_bk_poles(MaxVectsDim,MaxDySols),  X_fw_poles(MaxVectsDim,MaxDySols)
  real*8 :: E_Fermi

  !
  !  We will used the following global flags to select the level
  ! of approximation to be used.
  logical :: ADC2TDA, ADC3, ADC3CCD, TsqLadd, FullCCD

contains


subroutine Build_bases(xdlt)
  !
  !  Builds the single particle basis for the multi-level system. Since each
  ! energy level has spin 1/2, the total number of states will be twice the
  ! number of levels.
  !
  !  This can be controlled with the variable N_1b_mdsp which is the total number
  ! of basis states and
  !
  !  xdlt is the level spacing (the \xi of chapter 11)
  !
  implicit none
  real*8, intent(in) :: xdlt

  integer :: i, i1, i2, i3, k
  integer :: n1,n2

  N_1b_bas = 0
  N_2p_bas = 0
  N_2h_bas = 0
  N_2p1h_bas = 0
  N_2h1p_bas = 0

  write(6,*)
  write(6,*) ' Building the single particle basis...  (''p'' is the level number)'

  N_1b_mdsp=8  !  2 (spin orientations) * 4 (energy levels)

  write(6,*) '          #           p          m_s         occ     e(p) '

  i1 = 0
  do i = 1 , N_1b_mdsp

    ! Sets the occupation. Note that the number of particle must be even.
    isp_oc(i) = 0
    if (i <= 4) isp_oc(i) = 1 ! Here, 4 is the # of particles

    isp_lv(i) = (i-1)/2           ! level number
    isp_sp(i) = -1 + 2*MOD(i,2)   ! spin orientation
    xsp_en(i) = isp_lv(i) * xdlt  ! level energy


    write(6,*) i, isp_lv(i), isp_sp(i), isp_oc(i),  xsp_en(i)

    !
    ! Counts the energy levels in the basis
    if (0 < isp_sp(i)) then
     i1=i1+1
     i_1b_bas(i1) = i
     i_1b_occ(i1) = isp_oc(i)
    end if

  end do

  N_1b_bas = i1  ! eventually it must be   N_1b_mdsp = 2 * N_1b_bas


  !
  ! Builds tables for all the 2p1h and 2h1p intermediate state configurations (ISCs)
  !
  i1 = 0
  i2 = 0
  do n1 =  1 , N_1b_mdsp
  do n2 = n1+1 , N_1b_mdsp
   do k  = 1 , N_1b_mdsp

    if (isp_lv(n1) /=  isp_lv(n2)) cycle
    if (isp_sp(n1) /= -isp_sp(n2)) cycle
    if (isp_sp(k)  > 0) cycle

    if ( (0 == isp_oc(n1)) .AND. (0 == isp_oc(n2)) .AND. (1 == isp_oc(k)) ) then
      i1=i1+1
      innk_n1(i1) = n1
      innk_n2(i1) = n2
      innk_k3(i1) = k
      xnnk_en(i1) = xsp_en(n1) + xsp_en(n2) - xsp_en(k)
    end if

    if ( (1 == isp_oc(n1)) .AND. (1 == isp_oc(n2)) .AND. (0 == isp_oc(k)) ) then
      i2=i2+1
      ikkn_k1(i2) = n1
      ikkn_k2(i2) = n2
      ikkn_n3(i2) = k
      xkkn_en(i2) = xsp_en(n1) + xsp_en(n2) - xsp_en(k)
    end if

   end do
  end do
  end do

  N_2p1h_bas = i1
  N_2h1p_bas = i2



  !
  !  Build the (Pauli ordered) 2qp and 2qh ISCs. These are used to calculate
  ! the corrections for the ADC(3) coupling and the CCD(ladder) resummations
  !if (ADC3) then

    i1 = 0
    i2 = 0
    do n1 =  1 , N_1b_mdsp
     do n2 = n1+1 , N_1b_mdsp

      if (isp_lv(n1) /=  isp_lv(n2)) cycle
      if (isp_sp(n1) /= -isp_sp(n2)) cycle

      if ( (0 == isp_oc(n1)) .AND. (0 == isp_oc(n2)) ) then
        i1=i1+1
        i_2p_n1(i1) = n1
        i_2p_n2(i1) = n2
        x_2p_en(i1) = xsp_en(n1) + xsp_en(n2)
      end if

      if ( (1 == isp_oc(n1)) .AND. (1 == isp_oc(n2)) ) then
        i2=i2+1
        i_2h_k1(i2) = n1
        i_2h_k2(i2) = n2
        x_2h_en(i2) = xsp_en(n1) + xsp_en(n2)
      end if

     end do
    end do
    N_2p_bas = i1
    N_2h_bas = i2

  !end if

end subroutine Build_bases


subroutine Refresh_ISCs_energies(verbose)
  !
  !  Recalculate the unperturbed energies of ISCs, to be used in SCGF iterations
  !
  implicit none
  logical, intent(in) :: verbose

  integer :: i, i1, i2, i3, k
  integer :: n1,n2


  do i1 = 1 , N_2p1h_bas
    n1 = innk_n1(i1)
    n2 = innk_n2(i1)
    k  = innk_k3(i1)
    xnnk_en(i1) = esp0_en(n1) + esp0_en(n2) - esp0_en(k)
  end do

  do i2 = 1 , N_2h1p_bas
    n1 = ikkn_k1(i2)
    n2 = ikkn_k2(i2)
    k  = ikkn_n3(i2)
    xkkn_en(i2) = esp0_en(n1) + esp0_en(n2) - esp0_en(k)
  end do


  if (verbose) then
    !
    !  'verbose' controls the amount of output messages
    !
    write(6,*)
    write(6,*) '  I am recalculating the unperturbed energies doe the ISCs:'
    write(6,*) '  ---------------------------------------------------------'
    write(6,*)
    write(6,*) ' Single particle states in the spin-up channel...'
    write(6,*)
    do i = 1 , N_1b_bas
      i1 = i_1b_bas(i)
      write(6,*) i,i_1b_bas(i),i_1b_occ(i),isp_lv(i1),isp_sp(i1),isp_oc(i1),esp0_en(i1)
    end do
    !
    write(6,*)
    write(6,*)
    write(6,*) ' 2p1h basis...'
    write(6,*)
    do i = 1 , N_2p1h_bas
      i1 = innk_n1(i)
      i2 = innk_n2(i)
      i3 = innk_k3(i)
      write(6,101) i,i1,i2,i3,isp_lv(i1),isp_sp(i1),isp_lv(i2),isp_sp(i2),isp_lv(i3),isp_sp(i3),xnnk_en(i)
    end do
101 format(i4,2x,3I3,' - (',i2,',',i2') (',i2,',',i2') (',i2,',',i2')',f10.4)
    !
    write(6,*)
    write(6,*)
    write(6,*) ' 2h1p basis...'
    write(6,*)
    do i = 1 , N_2h1p_bas
      i1 = ikkn_k1(i)
      i2 = ikkn_k2(i)
      i3 = ikkn_n3(i)
      write(6,101) i,i1,i2,i3,isp_lv(i1),isp_sp(i1),isp_lv(i2),isp_sp(i2),isp_lv(i3),isp_sp(i3),xkkn_en(i)
    end do
    write(6,*)
    !
  end if ! end if (verbose)...

  do i1 = 1 , N_2p_bas
    n1 = i_2p_n1(i1)
    n2 = i_2p_n2(i1)
    x_2p_en(i1) = esp0_en(n1) + esp0_en(n2)
  end do

  do i2 = 1 , N_2h_bas
    n1 = i_2h_k1(i2)
    n2 = i_2h_k2(i2)
    x_2h_en(i2) = esp0_en(n1) + esp0_en(n2)
  end do

  if (verbose) then
    !
    !  'verbose' controls the amount of output messages
    !
    write(6,*)
    write(6,*)
    write(6,*) ' 2p basis...'
    write(6,*)
    do i = 1 , N_2p_bas
      i1 = i_2p_n1(i)
      i2 = i_2p_n2(i)
      write(6,106) i,i1,i2,isp_lv(i1),isp_sp(i1),isp_lv(i2),isp_sp(i2),x_2p_en(i)
    end do
106 format(i4,2x,2I3,' - (',i2,',',i2') (',i2,',',i2')',f10.4)
    !
    write(6,*)
    write(6,*)
    write(6,*) ' 2h basis...'
    write(6,*)
    do i = 1 , N_2h_bas
      i1 = i_2h_k1(i)
      i2 = i_2h_k2(i)
      write(6,106) i,i1,i2,isp_lv(i1),isp_sp(i1),isp_lv(i2),isp_sp(i2),x_2h_en(i)
    end do
  end if

end subroutine Refresh_ISCs_energies


end module Bases


module CoupledCluster
  !
  !  Here we store the T2 CCD amplitudes and the relative 2p and 2h configurations
  ! so that they can be used in the calculation of ADC(3)-D vertices.
  !
  !  In principle exactly the same 2-particle and the 2-hole configurations
  ! enter *both* the CCD equations and the formulas for the ADC(3) vertices. However,
  ! the configurations built in the module Bases are optimized to take into account
  ! the symmetries of the pairing Hamiltonian.
  !  For CCD, we take a more straightforward approach and and release these constraints.
  ! This simplifies coding the CCD calculation (of course) but at the price that we will
  ! explicitly sum over many zero contributions and that we need to construct the most
  ! general 2p and 2h bases as well.  For a a problem of very small dimension as this
  ! 4-level pairing the extra computing time remains negligible.
  !
  !  This module also stores the CCD amplitudes to be used for calculating ADC(3)-D
  ! later on.
  !

  use Bases
  implicit none

  ! T2 are the CCD amplitudes
  real*8,  allocatable :: index_xsng(:,:) , T2(:,:)
  integer, allocatable :: index_ijab(:,:)

  ! These are the most general 2p and 2h configuration, to be constructed in an
  ! analogous way to the ones in Bases::Build_bases
  integer              :: N_CC2p_bas, N_CC2h_bas
  integer, allocatable ::  i_CC2p_n1(:), i_CC2h_k1(:)
  integer, allocatable ::  i_CC2p_n2(:), i_CC2h_k2(:)
  real*8,  allocatable ::  xe_CC2p(:), xe_CC2h(:)


contains

subroutine Build_CCD_bases

  implicit none

  integer :: i1, i2, n1, n2, n2bas

  !  mod sp indices are ordered because of Pauli, these tables are to recover the correct
  ! indices in the list of configuration give the two separate sp indices. If they are 
  ! inverted w.r.t. the ones in the list, index_xsng=-1 to account for antisymmetry.
  if ( allocated(index_ijab) ) deallocate(index_ijab)
  if ( allocated(index_xsng) ) deallocate(index_xsng)

  allocate ( index_ijab(N_1b_mdsp,N_1b_mdsp) )
  allocate ( index_xsng(N_1b_mdsp,N_1b_mdsp) )

  n2bas = ((N_1b_mdsp+1) * N_1b_mdsp) / 2  ! max possible numebr of 2p  of 2h

  if ( allocated( i_CC2p_n1 ) ) deallocate( i_CC2p_n1 ) ;    if ( allocated( i_CC2h_k1 ) ) deallocate( i_CC2h_k1 )
  if ( allocated( i_CC2p_n2 ) ) deallocate( i_CC2p_n2 ) ;    if ( allocated( i_CC2h_k2 ) ) deallocate( i_CC2h_k2 )
  if ( allocated(  xe_CC2p  ) ) deallocate(  xe_CC2p  ) ;    if ( allocated(  xe_CC2h  ) ) deallocate(  xe_CC2h  )

  allocate ( i_CC2p_n1(n2bas) ) ;                            allocate ( i_CC2h_k1(n2bas) )
  allocate ( i_CC2p_n2(n2bas) ) ;                            allocate ( i_CC2h_k2(n2bas) )
  allocate (   xe_CC2p(n2bas) ) ;                            allocate (   xe_CC2h(n2bas) )


  ! Now let's generate the 2p and 2h bases without symmetry constraints:
  i1 = 0
  i2 = 0
  index_ijab(:,:) = -100 ! < 0  means the configuration is not in the basis (e.g. n1=n2 is forbidden by Pauli)
  index_xsng(:,:) = 0.d0
  do n1 =  1 , N_1b_mdsp
    do n2 = n1+1 , N_1b_mdsp

    if ( (0 == isp_oc(n1)) .AND. (0 == isp_oc(n2)) ) then
      i1=i1+1
      i_CC2p_n1(i1) = n1
      i_CC2p_n2(i1) = n2
      xe_CC2p(i1) = esp0_en(n1) + esp0_en(n2)
      index_ijab(n1,n2) = i1;    index_xsng(n1,n2) =  1.d0
      index_ijab(n2,n1) = i1;    index_xsng(n2,n1) = -1.d0
    end if

    if ( (1 == isp_oc(n1)) .AND. (1 == isp_oc(n2)) ) then
      i2=i2+1
      i_CC2h_k1(i2) = n1
      i_CC2h_k2(i2) = n2
      xe_CC2h(i2) = esp0_en(n1) + esp0_en(n2)
      index_ijab(n1,n2) = i2;    index_xsng(n1,n2) =  1.d0
      index_ijab(n2,n1) = i2;    index_xsng(n2,n1) = -1.d0
    end if

    end do
  end do
  N_CC2p_bas = i1
  N_CC2h_bas = i2


  if ( allocated( T2 ) ) deallocate(  T2  );    allocate ( T2(N_CC2p_bas,N_CC2h_bas) )

  T2(:,:) = 0.d0


end subroutine Build_CCD_bases

end module CoupledCluster


program pairing
  !
  ! This is the main program
  !

  use Bases
  implicit none


  integer :: n1

  real*8  :: xdlt, xg, x1, x2, x3, x4

  integer :: n_sc_itrs,  i_sc_type

  real*8, external :: Exact_energy, Dyson, SolveCCD

  logical :: OptnAccptd


  xdlt = 1.d0
  xg = -1.0d0

  write(6,*)
  write(6,*) ' Value of g ? '
  read(5,*) xg


  ! pre initialise running flags
  ADC2TDA  = .false.  !  to summ the 2p1h/2h1p TDA ladders
  ADC3     = .false.  !  to add the ADC(3) coupling corrections
  ADC3CCD  = .false.  !  to add the ADC(3) coupling corrections in the CCD approx
  TsqLadd  = .false.  !  pure ladders in the CC equations
  FullCCD  = .false.  !  quatratic ladders in the CC equations (these are those generating a ladder-RPA)

  OptnAccptd = .false.
  do while(OptnAccptd)
    !
    OptnAccptd = .true.  ! Because we belive in the presumption of innocence
    !
    write(6,*)
    write(6,*) ' Type  1 for ADC(2),  2 for 2p1h-TDA,  3 for ADC(3)  and  4 for ADC(3)-D ? '
    read(5,*) n1
    if (1 == n1) then
      ADC2TDA = .false.
      ADC3    = .false.
      ADC3CCD = .false.
    else if (2 == n1) then
      ADC2TDA = .true.
      ADC3    = .false.
      ADC3CCD = .false.
    else if (3 == n1) then
      ADC2TDA = .true.
      ADC3    = .true.
     !ADC3CCD = .false.
    else if (4 == n1) then
      ADC2TDA = .true.
     !ADC3    = .false.
      ADC3CCD = .true.
    else
      ! input is not recognized
      OptnAccptd = .false.
    end if
  end do

  if (ADC3CCD) then ADC3 = .false. ! For safety: it must be EITHER one OR the other (see 'Build_ADCmtx()' below)

  OptnAccptd = .false.
  do while(OptnAccptd)
    !
    OptnAccptd = .true.  ! Because we belive in the presumption of innocence
    !
    write(6,*)
    write(6,*) ' Type 1 for CC(ladder)(2),  2 adding the TVT(ladders) or 3 for Full CCD ? '
    read(5,*) n1
    if (1 == n1) then
      TsqLadd = .false.
      FullCCD = .false.
    else if (2 == n1) then
      TsqLadd = .true.
      FullCCD = .false.
    else if (3 == n1) then
      TsqLadd = .true.
      FullCCD = .true.
    else
      ! input is not recognized
      OptnAccptd = .false.
    end if
  end do

  write(6,*)
  write(6,*) ' Type  0  for using unperturbed s.p. energies from H_0, anything else to use HF as the reference state: '
  read(5,*)  i_sc_type


  ! Build bases and ISCs for SCGF
  call Build_bases(xdlt)

  !
  ! If one want to plot the results to a file:
  !
  !open(UNIT=8,FILE='Pair_sol_CCM.dat')
  !open(UNIT=9,FILE='Pair_sol_Dyson.dat')
  !open(UNIT=10,FILE='Pair_sol_Exact.dat')

  !do n1 = 0, 52
  ! xg = -1.3d0 + n1* 0.05d0

  !n_sc_itrs = 0   ! Just do one dyagonalization with the unperturbed HF diagram
  n_sc_itrs = 10  ! Iterate the sc0 scheme a given number of times

  !
  !  Remember if ADC3CCD is set we MUST call CCD before doing Dyson because the
  ! coupling vertices in ADC(3)-D need the T2 amplitudes
  !
  x4 = SolveCCD(xg,x2,i_sc_type) ! x2 returns MBPT2

  x1 = Dyson(xg,n_sc_itrs, i_sc_type)

  x3 = Exact_energy(xdlt,xg)


  !
  ! If one want to plot the results to a file:
  !
  !write( 8,*) xg,      x4,     x2
  !write( 9,*) xg,   x1+xg-2.0, x2
  !write(10,*) xg,   x3+xg-2.0, x2
  !
  !end do
  !
  !close(UNIT=8)
  !close(UNIT=9)
  !close(UNIT=10)

end program pairing


real*8 function Dyson(xg,n_sc_itrs, i_sc_type)

  use Bases
  implicit none
  real*8,  intent(in) :: xg  !  xdlt is alrady set in module Bases by calling 'Build_bases(xdlt)'
  integer, intent(in) :: n_sc_itrs, i_sc_type

  ! This is the stuff needed by the LAPACK eigenvalue pakage:
  ! For 'DSYEV':
  integer, parameter :: LDA = 40
  integer, parameter :: LWORK  = 2+(6+2*LDA)*(LDA)  ! For 'DSYEVD' only...
  integer, parameter :: LIWORK = 3+5*LDA            ! For 'DSYEVD' only...
  real*8  :: A(LDA,LDA), B(LDA,LDA)
  real*8  :: W(LDA), WORK(LWORK)
  integer :: IWORK(LIWORK)

  integer :: INFO, LWopt, LIWopt

  integer :: i, i1, i2, i3, i4, j, k, n1,   ih1, ih2

  real*8  :: x1, x2, xA, xSF, xNorm, xKolt, xA_run, xKolt_run, x3,x4 ,x5

  integer :: Ntot, itr, max_itrs


  real*8, external :: Vpair


  write(6,*)
  write(6,*) ' --- --- ---'
  write(6,*)
  write(6,*) '  Start calculations with the GF method...'
  write(6,*)

  !
  ! This is the HF potential, which is diagonal in s.p. basis
  !
  vhf_en(:) = 0.d0
  do j = 1 , N_1b_mdsp
   x1 = 0.d0
   do k = 1 , N_1b_mdsp
    if (1 == isp_oc(k)) x1 = x1 + Vpair(j,k,j,k,xg)
   end do
   vhf_en(j) = x1
  end do


  Ntot = N_1b_bas + N_2p1h_bas + N_2h1p_bas

  esp0_en(:) = xsp_en(:)
  Utld_en(:) = 0.d0
  if (0 == i_sc_type) then
    !
    !  Use  the energy levels for the unprturbed H_0:
    esp0_en(:) = xsp_en(:)
    Utld_en(:) = vhf_en(:)

  else if (1 <= i_sc_type) then
    !
    !  Use the Hf basis to define H_0 (here Vhf is diagonal!)
    esp0_en(:) = xsp_en(:) + vhf_en
    Utld_en(:) = 0.d0

  end if

  call Refresh_ISCs_energies(.true.)
  call Build_ADCmtx(B,LDA,xg)


  x4 = 0.d0; x5 = 0.d0
  vhf_3rd(:) = 0.d0
  do k = 1 , N_1b_mdsp

    i1 = (k-1)/2 + 1;
    i4 = i_1b_bas(i1)
    if ( (k /= i4 ).and.(k /= 1+i4 ) ) then; write(6,*) ' Some trouble occurred ',i1, j,i_1b_bas(i1)  ; stop; end if
    if (1 == isp_oc(k)) then
      !
      ! Must connect with the 2p1h states
      !
      x1 = 0.d0
      do j = 1 , N_2p1h_bas
        x3 = Vpair(i4,innk_k3(j),innk_n1(j),innk_n2(j),xg)
        x2 = x3 / (-esp0_en(k) + xnnk_en(j) )
        x1 = x1 + x2*x2
        x4 = x4 - x2*x3
      end do


    else
      !
      ! Must connect with the 2h1p states
      !
      x1 = 0.d0
      do j = 1 , N_2h1p_bas
        x3 = Vpair(ikkn_k1(j),ikkn_k2(j),i4,ikkn_n3(j),xg)
        x2 = x3 / (-esp0_en(k) + xkkn_en(j) )
        x1 = x1 - x2*x2
        x5 = x5 + x3*x2
      end do

    end if
    vhf_3rd(k) = x1
  end do

  write(6,*)
  write(6,*) 'MBPT2 = ', x4/2.d0, x5/2.d0
  write(6,*)

  !
  !  Seek for an estimation of the Fermi energy,  to be used later
  ! for separating qps and qhs
  x1 = -1.d20
  x2 = +1.d20
  do  i = 1 , N_1b_bas
    k = i_1b_bas(i)
    if (1 == isp_oc(k)) x1 = MAX(x1,xsp_en(k)+vhf_en(k))
    if (0 == isp_oc(k)) x2 = MIN(x2,xsp_en(k)+vhf_en(k))
  end do
  E_Fermi = (x1+x2)/2.d0
  write(6,'(A,f10.4)') ' The Fermi energy will be taken to be ',E_Fermi




  max_itrs = MAX(0, n_sc_itrs)
  do itr = 0 , max_itrs

    ! load the part fothe Dyson matrix that will not change with iterations
    A = B

    !
    ! Add the HF(cHF) self energy:
    !
    do i = 1 , N_1b_bas
     i1 = i_1b_bas(i)
     do j = 1 , N_1b_bas
      i2 = i_1b_bas(j)

      if (0 == itr) then
        !
        !  At the first iterations, just use the 1st
        ! order HF diagram
        x1 = 0.d0
        do k = 1 , N_1b_mdsp
          if (1 == isp_oc(k)) x1 = x1 + Vpair(i1,k,i2,k,xg)
        end do
      else
        !
        !  Now we have a propagator and thus we can do
        ! the fully correlated 1
        x1 = 0.d0
        do i3 = 1 , N_1b_bas
         do i4 = 1 , N_1b_bas
          do k  = 1 , N_bk_poles
            ! 
            !  We are using a trick here: in principle we should be calculating also the
            ! spin down part of the spectral function (which we are skipping because it is
            ! exactly the same as the the spin up part--in this model). In fact, both spin
            ! up and down spectral functions can contribute to the static self-energy. For
            ! this particular pairing model, it is only the spin down that contributes to the
            ! spin up HF diagram (because of the interaction we have chosen)!!!
            !  Here we use the calculated spin-up part in which is stored in Y_bk_poles but
            ! then we need to use the corresponding spin-down s.p. states in for the
            ! interaction matrix elements (hence the i_1b_bas(i...)+1 indices).
            x1 = x1 + Vpair(i1,i_1b_bas(i3)+1,i2,i_1b_bas(i4)+1,xg) * Y_bk_poles(k,i3)*Y_bk_poles(k,i4)
          end do
         end do
        end do
      end if

      A(i,j) = A(i,j) + x1
     end do

     !if ((0 == itr) .and. (ADC3 .or. ADC3CCD)) A(i,i) = A(i,i) + vhf_3rd(i1)
    end do


    !
    ! Show the initial and final Dyson matrix
    if ((max_itrs == itr) .OR. (0 == itr)) then
      write(6,*)
      write(6,*)
      if (0 == itr) then
        write(6,*) ' The initial Dyson matrix is:'
      else
        write(6,*) ' The Dyson matrix is:'
      end if
        write(6,*)
        do i = 1 , Ntot
          write(6,110) i, A(i,1:Ntot)
          if (0 == MOD(i,4)) write(6,*)
        end do
    end if
110 format(i4,3(3X,4f7.3))




    ! LAPACK library (w/ DSYEVD):
    INFO = 0;
    LWopt=0
    LIWopt=0
    call dsyevd('Vectors','Upper',Ntot,A,LDA,W,WORK, LWORK,IWORK,LIWORK,INFO);
    if (0 /= INFO) write(6,*) 'Dyson (DSYEV), Wrong value of IFAIL: INFO= ',INFO
    if (0 == INFO) then
      LWopt = MAX(int(WORK(1) + 0.01) , LWopt)
      LIWopt = MAX(IWORK(1) , LIWopt)
    else
      write(6,*) ' ''DEEGV'' gave INFO = ', INFO
      if (INFO < 0) write(6,*) ' The ',-INFO,'-th aggument in line 408 was illegal:'
      write(6,*) ' Program has been aborted since IERR != 0.'
      write(6,*) ' Aborted at the line 410, Dyson.f!'
      write(6,*) '   --> stop.'
      stop
    end if

    if ((max_itrs == itr) .OR. (0 == itr)) then
      write(6,*)
      write(6,*) '              e^+/-      \sum A        SF       \sum KRS          KSR      Corr. En                Check norm'
    end if

    xA = 0.d0
    xKolt = 0.d0
    xA_run = 0.d0
    xKolt_run = 0.d0
    N_bk_poles = 0
    N_fw_poles = 0
    do i = 1 , Ntot

      xNorm = 0.d0; xSF = 0.d0
      x2 = 0.d0
      do j = 1 , Ntot
        x1 = A(j,i)
        xNorm = xNorm + x1*x1
        if (j <= N_1b_bas) then
          xSF = xSF + x1*x1
!          write(6,*) j,xsp_en(i_1b_bas(j))
          x2  = x2  + x1*x1*(xsp_en(i_1b_bas(j)) + W(i))
        end if
      end do
      x2 = x2 / xNorm / 2.0
      xSF   = xSF / xNorm
      xNorm = DSQRT(xNorm)

      if (W(i) < E_Fermi) then
        ! quasiholes
         N_bk_poles = N_bk_poles + 1
         E_bk_poles(N_bk_poles) = W(i)
        SF_bk_poles(N_bk_poles) = xSF
         Y_bk_poles(N_bk_poles,1:N_1b_bas) = A(1:N_1b_bas,i) / xNorm
        xKolt = xKolt + x2
        xA    = xA  + xSF
      else
        ! quasiparticles
         N_fw_poles = N_fw_poles + 1
         E_fw_poles(N_fw_poles) = W(i)
        SF_fw_poles(N_fw_poles) = xSF
         X_fw_poles(N_fw_poles,1:N_1b_bas) = A(1:N_1b_bas,i) / xNorm
      end if

      xKolt_run = xKolt_run + x2
      xA_run    = xA_run  + xSF


      if ((max_itrs == itr) .OR. (0 == itr)) then
         write(6,120) i, W(i),  xA_run, xSF, xKolt_run, xKolt_run*2.d0, xKolt_run*2.d0 - 2.d0 + xg, x2, xNorm
         if (0 == MOD(i,6)) write(6,*)
      end if
120   format(i4,3(4X,4f12.5))

    end do

    !
    !  We have calculated only the spectral function for particles with
    ! spin up  but the spin spin down would be exactly the same and give
    ! another equal contribution. The factor 2.d0 is to account for this.
    !
    write(6,121) ' Particle number and Koltun SR :  ', 2.d0*xA, 2.d0*xKolt,   &
                        '    (corr. en. = ',2*(xKolt-1.0)+xg,')'
121 format(A,2f12.5,A,f12.5,A)

    Dyson = 2.d0*xKolt
  end do ! itr loop

end function Dyson


subroutine Build_ADCmtx(ADCmtx,LDA,xg)
  use Bases
  use CoupledCluster
  implicit none
  integer, intent(in)  :: LDA
  real*8,  intent(out) :: ADCmtx(LDA,LDA)
  real*8,  intent(in)  :: xg

  integer :: i, i4, j, k, n1, Ntot,  ih1, ih2,  i_ab

  real*8  :: x1, x_ab

  real*8, external :: Vpair


  ADCmtx(:,:) = 0.d0
  Ntot = N_1b_bas + N_2p1h_bas + N_2h1p_bas
  do i = 1 , N_1b_bas

    i4 = i_1b_bas(i)

    ADCmtx(i,i) = ADCmtx(i,i) + xsp_en(i4)

    do j = 1 , N_2p1h_bas
      x1 = Vpair(i4,innk_k3(j),innk_n1(j),innk_n2(j),xg)
      ADCmtx(i,N_1b_bas+j) = x1
      ADCmtx(N_1b_bas+j,i) = x1
    end do

    do j = 1 , N_2h1p_bas
      x1 = Vpair(ikkn_k1(j),ikkn_k2(j),i4,ikkn_n3(j),xg)
      ADCmtx(i,N_1b_bas+N_2p1h_bas+j) = x1
      ADCmtx(N_1b_bas+N_2p1h_bas+j,i) = x1
    end do

  end do

  do j = 1 , N_2p1h_bas
    n1 = N_1b_bas+j
    ADCmtx( n1 , n1 ) = xnnk_en(j) + Utld_en(innk_n1(j)) + Utld_en(innk_n2(j)) - Utld_en(innk_k3(j))
  end do
  do j = 1 , N_2h1p_bas
    n1 = N_1b_bas+N_2p1h_bas+j
    ADCmtx( n1 , n1 ) = xkkn_en(j) + Utld_en(ikkn_k1(j)) + Utld_en(ikkn_k2(j)) - Utld_en(ikkn_n3(j))
  end do


  if (ADC2TDA) then
    !
    !  Add the pp and hh ladders that define the extended-ADC(2) approximation
    ! and are also included in the ADC(3) method.
    !
    !  A comment: diagonalizing these sub matrices separately, would let you
    ! calculate the very same pp or hh ladders that were done during the
    ! coupled cluster CCD exercises.
    !
    n1 = N_1b_bas
    do i = 1 , N_2p1h_bas
     do j = 1 , N_2p1h_bas
      if (innk_k3(i) /= innk_k3(j)) cycle
      x1 = Vpair(innk_n1(i),innk_n2(i),innk_n1(j),innk_n2(j),xg)
      ADCmtx( n1 + i , n1 + j ) = ADCmtx( n1 + i , n1 + j ) + x1
     end do
    end do
    !
    n1 = N_1b_bas + N_2p1h_bas
    do i = 1 , N_2h1p_bas
     do j = 1 , N_2h1p_bas
      if (ikkn_n3(i) /= ikkn_n3(j)) cycle
      x1 = Vpair(ikkn_k1(i),ikkn_k2(i),ikkn_k1(j),ikkn_k2(j),xg)
      ADCmtx( n1 + i , n1 + j ) = ADCmtx( n1 + i , n1 + j ) - x1
     end do
    end do
    !
  end if

  if (ADC3) then
    !
    !  Add the ADC(3) corrections for the \alpha-2p1h and \alpha-2h1p
    ! coupling terms.
    !  ONLY the pp/hh ladder corrections are included here because the
    ! rings do not contribute to this model...(?)
    !
    !
    do i = 1 , N_1b_bas
     i4 = i_1b_bas(i)
     do j = 1 ,  N_2p1h_bas

      x1 = 0.d0
      do k = 1 , N_2h_bas
        ih1 = i_2h_k1(k)
        ih2 = i_2h_k2(k)
        x1 = x1 + Vpair(i4,innk_k3(j),ih1,ih2,xg) * Vpair(ih1,ih2,innk_n1(j),innk_n2(j),xg) /   &
            (x_2h_en(k) - esp0_en(innk_n1(j)) - esp0_en(innk_n2(j)) )
      end do

      ADCmtx(i,N_1b_bas+j) = ADCmtx(i,N_1b_bas+j) + x1
      ADCmtx(N_1b_bas+j,i) = ADCmtx(N_1b_bas+j,i) + x1

     end do
    end do
    !
    do i = 1 , N_1b_bas
     i4 = i_1b_bas(i)
     do j = 1 , N_2h1p_bas

      x1 = 0.d0
      do k = 1 , N_2p_bas
        ih1 = i_2p_n1(k)
        ih2 = i_2p_n2(k)
        x1 = x1 + Vpair(ikkn_k1(j),ikkn_k2(j),ih1,ih2,xg) * Vpair(ih1,ih2,i4,ikkn_n3(j),xg) /   &
            (esp0_en(ikkn_k1(j)) + esp0_en(ikkn_k2(j)) - x_2p_en(k))
      end do

      ADCmtx(i,N_1b_bas+N_2p1h_bas+j) = ADCmtx(i,N_1b_bas+N_2p1h_bas+j) + x1
      ADCmtx(N_1b_bas+N_2p1h_bas+j,i) = ADCmtx(N_1b_bas+N_2p1h_bas+j,i) + x1
     end do
    end do
    !
  else if (ADC3CCD) then
    !
    !  Add the ADC(3)-D corrections for the \alpha-2p1h and \alpha-2h1p
    ! coupling terms.
    !  ONLY the pp/hh ladder corrections are included here because the
    ! rings do not contribute to this model...(?)
    !
    !
    do i = 1 , N_1b_bas
     i4 = i_1b_bas(i)
     do j = 1 ,  N_2p1h_bas

       !  Get the index in the extended 2p/2h bases (used by the CCD program) and also
       ! read the symmetry factor to account for the -1 if the fermion n1 and n2 are inverted
       i_ab = index_ijab(innk_n1(j),innk_n2(j));   x_ab = index_xsng(innk_n1(j),innk_n2(j))
       if (i_ab < 0) cycle  ! if 1_ab < 0 that means that this configuration does not exists.

       x1 = 0.d0
       do k = 1 , N_CC2h_bas
         ih1 = i_CC2h_k1(k)
         ih2 = i_CC2h_k2(k)
         x1 = x1 + Vpair(i4,innk_k3(j),ih1,ih2,xg) * T2(i_ab,k) * x_ab
       end do

       ADCmtx(i,N_1b_bas+j) = ADCmtx(i,N_1b_bas+j) + x1
       ADCmtx(N_1b_bas+j,i) = ADCmtx(N_1b_bas+j,i) + x1

     end do
    end do
    !
    do i = 1 , N_1b_bas
     i4 = i_1b_bas(i)
     do j = 1 , N_2h1p_bas

       i_ab = index_ijab(ikkn_k1(j),ikkn_k2(j));   x_ab = index_xsng(ikkn_k1(j),ikkn_k2(j))
       if (i_ab < 0) cycle

      x1 = 0.d0
      do k = 1 , N_CC2p_bas
        ih1 = i_CC2p_n1(k)
        ih2 = i_CC2p_n2(k)
        x1 = x1 + x_ab * T2(k,i_ab) * Vpair(ih1,ih2,i4,ikkn_n3(j),xg)
      end do

      ADCmtx(i,N_1b_bas+N_2p1h_bas+j) = ADCmtx(i,N_1b_bas+N_2p1h_bas+j) + x1
      ADCmtx(N_1b_bas+N_2p1h_bas+j,i) = ADCmtx(N_1b_bas+N_2p1h_bas+j,i) + x1
     end do
    end do
    !
   end if

end subroutine Build_ADCmtx


real*8 function Exact_energy(xdlt,xg)
  !
  ! Calculates all the exact eigenstates by full CI.
  !
  implicit none
  real*8, intent(in) :: xdlt, xg

  integer :: Ntot

  ! The stuff needed by the LAPACK eigenvalue pakage:
  ! For 'DSYEV':
  integer, parameter :: LDA = 6
  integer, parameter :: LWORK  = 2+(6+2*LDA)*(LDA)  ! For 'DSYEVD' only...
  integer, parameter :: LIWORK = 3+5*LDA              ! For 'DSYEVD' only...
  real*8  :: A(LDA,LDA)
  real*8  :: W(LDA), WORK(LWORK)
  integer :: IWORK(LIWORK)

  integer :: INFO, LWopt, LIWopt, i


  write(6,*)
  write(6,*) ' --- --- ---'
  write(6,*)
  write(6,*) '  Solving the exact ground state energy by diagonalising the '
  write(6,*) ' following hamiltonian:'
  write(6,*)


  A(:,:) = -xg / 2.d0
  Ntot = 6
  do i = 1 , Ntot
    A(i, Ntot+1-i) = 0.d0
    A(i,i) = 2.d0 * i - xg
    if (i > Ntot/2) A(i,i) = xdlt * 2.d0 * (i-1) - xg
  end do

  do i = 1 , Ntot
  write(6,130) i, A(i,1:Ntot)
  if (0 == MOD(i,3)) write(6,*)
  end do
130 format(i4,3(3X,3f7.3))



  ! LAPACK library (w/ DSYEVD):
  INFO = 0;
  LWopt=0
  LIWopt=0
  Ntot = 6
  call dsyevd('Vectors','Upper',Ntot,A,LDA,W,WORK, LWORK,IWORK,LIWORK,INFO);
  if (0 /= INFO) write(6,*) 'Dyson (DSYEV), Wrong value of IFAIL: INFO= ',INFO
  if (0 == INFO) then
    LWopt = MAX(int(WORK(1) + 0.01) , LWopt)
    LIWopt = MAX(IWORK(1) , LIWopt)
  else
    write(6,*) ' ''DEEGV'' gave INFO = ', INFO
    if (INFO < 0) write(6,*) ' The ',-INFO,'-th aggument in line 408 was illegal:'
    write(6,*) ' Program has been aborted since IERR != 0.'
    write(6,*) ' Aborted at the line 410, Dyson.f!'
    write(6,*) '   --> stop.'
    stop
  end if

  Exact_energy = W(1)


  write(6,*)
  write(6,*)
  write(6,142) ' The exact correlation energy is: ', W(1) - 2.d0 + xg
  write(6,*)
  write(6,*)
  write(6,*) 'Eigenvalues:'
  write(6,140) W(1:Ntot)
  write(6,*)
  write(6,*) 'Eivenvectors:'
  do i = 1 , Ntot
    write(6,141) i, A(i,1:Ntot)
  end do
  write(6,*)

140 format(4X,5(3X,3f10.4))
141 format(i4,5(3X,3f10.4))
142 format(A,f12.5)

  return
end function Exact_energy



real*8 function SolveCCD(xg,dEmbpt2,i_sc_type)
  use Bases
  use CoupledCluster

  implicit none
  real*8, intent(in)  :: xg  !  xdlt is alrady set in module Bases by calling 'Build_bases(xdlt)'
  real*8, intent(out) :: dEmbpt2
  integer, intent(in) :: i_sc_type

  real*8, allocatable :: V2(:,:), E0(:,:), Vpp(:,:), Vhh(:,:), TT(:,:)
  real*8, allocatable :: TVX1(:,:), TVX2(:,:), u1b_no(:)

  integer ::  r, q, itr, ia, ib, ii, ij
  real*8  ::  xv, xe, xe_new
  real*8  ::  x1, x2, x3, x4, x5
  integer ::  r2,q2, ik, il, ic, id
  integer :: i1, i2, n1, n2

  integer ::  i_ac, i_ad, i_bc, i_bd,  i_ik, i_il, i_jk, i_jl
  real*8  ::  x_ac, x_ad, x_bc, x_bd,  x_ik, x_il, x_jk, x_jl


  real*8, external :: Vpair


  write(6,*)
  write(6,*) ' --- --- ---'
  write(6,*)
  write(6,*) '  Solving the CC equation for the correlation. This will also yield'
  write(6,*) ' the MBPT2 result as its 0-th iteration...'
  write(6,*)
       if (FullCCD) then; write(6,*) ' The full CCD equations will be solved.'
  else if (TsqLadd) then; write(6,*) ' The pp-hh RPA term will be added.'
  else;                   write(6,*) ' Ony the pp/hh ladders will be included.'
  end if
  write(6,*)


  allocate (u1b_no(MaxVectsDim))

  !
  ! Calculate the HF potential (here Vhf is diagonal!)
  vhf_en(:) = 0.d0
  do r = 1 , N_1b_mdsp
   x1 = 0.d0
   do q = 1 , N_1b_mdsp
    if (1 == isp_oc(q)) x1 = x1 + Vpair(r,q,r,q,xg)
   end do
   vhf_en(r) = x1
  end do

  !
  !  Normal ordered 1b mtx elment is:
  esp0_en(:) = xsp_en(:)
  u1b_no(:)  = vhf_en(:)
  if (1 <= i_sc_type) then
    !
    !  In this case, use the HF basis to define H_0
    esp0_en(:) = xsp_en(:) + vhf_en(:)
    u1b_no(:) = 0.d0
  end if


  call Build_CCD_bases

  ! T2(N_CC2p_bas,N_CC2h_bas)  has already been allocated by 'Build_CCD_bases'
  allocate (  V2(N_CC2h_bas,N_CC2p_bas) )
  allocate (  E0(N_CC2p_bas,N_CC2h_bas) )

  allocate ( TVX1(N_CC2p_bas,N_CC2p_bas) )
  allocate ( TVX2(N_CC2h_bas,N_CC2h_bas) )

  allocate ( Vpp(N_CC2p_bas,N_CC2p_bas) )
  allocate ( Vhh(N_CC2h_bas,N_CC2h_bas) )

  dEmbpt2 = 0.d0
  do r = 1 , N_CC2p_bas
    ia = i_CC2p_n1(r)
    ib = i_CC2p_n2(r)
    do q = 1 , N_CC2h_bas
      ii = i_CC2h_k1(q)
      ij = i_CC2h_k2(q)

      xv = Vpair(ia,ib,ii,ij,xg)
      xe = xe_CC2h(q) - xe_CC2p(r)
      V2(q,r) = xv
      E0(r,q) = xe
      T2(r,q) = xv / xe

      dEmbpt2 = dEmbpt2 + xv*xv/xe
    end do
  end do

  do r = 1 , N_CC2p_bas
   ia = i_CC2p_n1(r)
   ib = i_CC2p_n2(r)
   do q = 1 , N_CC2p_bas
    ii = i_CC2p_n1(q)
    ij = i_CC2p_n2(q)
    Vpp(r,q) = Vpair(ia,ib,ii,ij,xg)
   end do
  end do
  do r = 1 , N_CC2h_bas
   ia = i_CC2h_k1(r)
   ib = i_CC2h_k2(r)
   do q = 1 , N_CC2h_bas
    ii = i_CC2h_k1(q)
    ij = i_CC2h_k2(q)
    Vhh(r,q) = Vpair(ia,ib,ii,ij,xg)
   end do
  end do

  write(6,*)
  write(6,*) ' E_corr from MBPT2 : ',dEmbpt2
  write(6,*)

  write(6,*) '                     dE_old           dE_new'

  xe_new = dEmbpt2
  itr = 0
  do while (ABS(xe - xe_new) > 1.d-7)

    xe = xe_new
    itr = itr + 1

    TT = TRANSPOSE( V2 ) + MATMUL(T2, Vhh ) + MATMUL(Vpp, T2)

    if (TsqLadd .or. FullCCD) TT = TT + MATMUL( MATMUL (T2 , V2) , T2 )

    !
    !  The TVT ring term is zero for this model because the pairing interaciton
    ! does not connect p-h states
    !

    if (FullCCD) then
      !
      ! The remaining thems needed to complete the CCD approximation
      !
      TVX1 = MATMUL (T2 , V2)
      TVX2 = MATMUL (V2 , T2)
      !
      do r = 1 , N_CC2p_bas
        ia = i_CC2p_n1(r)
        ib = i_CC2p_n2(r)
        do q = 1 , N_CC2h_bas
          ii = i_CC2h_k1(q)
          ij = i_CC2h_k2(q)
          x1 = 0.d0
          x2 = 0.d0
          x3 = 0.d0
          do r2 = 1 , N_CC2p_bas
            ic = i_CC2p_n1(r2)
            id = i_CC2p_n2(r2)
            i_ac = index_ijab(ia,ic);   x_ac = index_xsng(ia,ic)
            i_ad = index_ijab(ia,id);   x_ad = index_xsng(ia,id)
            i_bc = index_ijab(ib,ic);   x_bc = index_xsng(ib,ic)
            i_bd = index_ijab(ib,id);   x_bd = index_xsng(ib,id)

            if ((i_ac > 0) .and. (i_bd > 0)) x1 = x1 - x_ac * x_bd * TVX1(i_ac,r2) * T2(i_bd,q)
            if ((i_ad > 0) .and. (i_bc > 0)) x1 = x1 + x_ad * x_bc * TVX1(i_ad,r2) * T2(i_bc,q)
            if ((i_bc > 0) .and. (i_ad > 0)) x1 = x1 + x_bc * x_ad * TVX1(i_bc,r2) * T2(i_ad,q)
            if ((i_bd > 0) .and. (i_ac > 0)) x1 = x1 - x_bd * x_ac * TVX1(i_bd,r2) * T2(i_ac,q)
            do q2 = 1 , N_CC2h_bas
              ik = i_CC2h_k1(q2)
              il = i_CC2h_k2(q2)
              i_ik = index_ijab(ii,ik);   x_ik = index_xsng(ii,ik)
              i_il = index_ijab(ii,il);   x_il = index_xsng(ii,il)
              i_jk = index_ijab(ij,ik);   x_jk = index_xsng(ij,ik)
              i_jl = index_ijab(ij,il);   x_jl = index_xsng(ij,il)

              x4 = 0.d0; x5 = 0.d0
              if ((i_ac > 0) .and. (i_bd > 0)) then
                if ((i_ik > 0) .and. (i_jl > 0)) x4 = x4 + x_ac * x_bd * x_ik * x_jl * T2(i_ac,i_ik) * T2(i_bd,i_jl)
                if ((i_il > 0) .and. (i_jk > 0)) x4 = x4 - x_ac * x_bd * x_il * x_jk * T2(i_ac,i_il) * T2(i_bd,i_jk)
                if ((i_jk > 0) .and. (i_il > 0)) x4 = x4 - x_ac * x_bd * x_jk * x_il * T2(i_ac,i_jk) * T2(i_bd,i_il)
                if ((i_jl > 0) .and. (i_ik > 0)) x4 = x4 + x_ac * x_bd * x_jl * x_ik * T2(i_ac,i_jl) * T2(i_bd,i_ik)
              end if
              if ((i_bc > 0) .and. (i_ad > 0)) then
                if ((i_ik > 0) .and. (i_jl > 0)) x5 = x5 - x_bc * x_ad * x_ik * x_jl * T2(i_bc,i_ik) * T2(i_ad,i_jl)
                if ((i_il > 0) .and. (i_jk > 0)) x5 = x5 + x_bc * x_ad * x_il * x_jk * T2(i_bc,i_il) * T2(i_ad,i_jk)
                if ((i_jk > 0) .and. (i_il > 0)) x5 = x5 + x_bc * x_ad * x_jk * x_il * T2(i_bc,i_jk) * T2(i_ad,i_il)
                if ((i_jl > 0) .and. (i_ik > 0)) x5 = x5 - x_bc * x_ad * x_jl * x_ik * T2(i_bc,i_jl) * T2(i_ad,i_ik)
              end if
              x3 = x3 + V2(q2,r2) * (x4 + x5)
            end do
          end do
          do q2 = 1 , N_CC2h_bas
            ik = i_CC2h_k1(q2)
            il = i_CC2h_k2(q2)
            i_ik = index_ijab(ii,ik);   x_ik = index_xsng(ii,ik)
            i_il = index_ijab(ii,il);   x_il = index_xsng(ii,il)
            i_jk = index_ijab(ij,ik);   x_jk = index_xsng(ij,ik)
            i_jl = index_ijab(ij,il);   x_jl = index_xsng(ij,il)

            if ((i_ik > 0) .and. (i_jl > 0)) x2 = x2 - x_ik * x_jl * T2(r,i_jl) * TVX2(q2,i_ik)
            if ((i_il > 0) .and. (i_jk > 0)) x2 = x2 + x_il * x_jk * T2(r,i_jk) * TVX2(q2,i_il)
            if ((i_jk > 0) .and. (i_il > 0)) x2 = x2 + x_jk * x_il * T2(r,i_il) * TVX2(q2,i_jk)
            if ((i_jl > 0) .and. (i_ik > 0)) x2 = x2 - x_jl * x_ik * T2(r,i_ik) * TVX2(q2,i_jl)

          end do

          x4 = (u1b_no(ia) + u1b_no(ib) - u1b_no(ii) - u1b_no(ij) ) * T2(r,q)

          TT(r,q) = TT(r,q) + (x1 + x2 + x3 + x4)

        end do
      end do
    end if


    xe_new = 0.d0
    do r = 1 , N_CC2p_bas
     do q = 1 , N_CC2h_bas
      TT(r,q) = TT(r,q) / E0(r,q)
      T2(r,q) = 0.5*T2(r,q) + 0.5*TT(r,q)
      xe_new = xe_new + V2(q,r)*T2(r,q)
     end do
    end do

    write(6,'(A,I5,2f17.7)') 'Itr n, ',itr, xe, xe_new

  end do

  write(6,*)
  write(6,*) ' The converged CC corr. energy is : ',xe_new
  write(6,*)

  SolveCCD = xe_new

  deallocate (u1b_no)

  !!!deallocate (  T2 )   this is needed by ADC(3)-D
  deallocate (  V2 )
  deallocate (  E0 )
  deallocate ( Vpp )
  deallocate ( Vhh )
  deallocate ( TVX1 )
  deallocate ( TVX2 )
  !!!deallocate ( index_ijab )   this is needed by ADC(3)-D
  !!!deallocate ( index_xsng )   this is needed by ADC(3)-D

  return
end function SolveCCD


real*8 function Vpair(a,b,c,d,xg)
  !
  ! Matrix elements of the pairing interaction
  !
  use Bases
  implicit none
  integer, intent(in) :: a,b,c,d
  real*8,  intent(in) :: xg
  Vpair = 0.d0
  if (isp_lv(a) /=  isp_lv(b)) return
  if (isp_sp(a) /= -isp_sp(b)) return
  if (isp_lv(c) /=  isp_lv(d)) return
  if (isp_sp(c) /= -isp_sp(d)) return
  Vpair = -xg/2.d0
  if (isp_sp(a) /=  isp_sp(c)) Vpair = -Vpair
  return
end function Vpair
