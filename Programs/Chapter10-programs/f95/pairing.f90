!             program pairing.f90
!
!             author:   morten hjorth-jensen
!             e-mail:   hjensen@msu.edu
!             language: f90/f95  
!             last upgrade : July 2016
!
!

MODULE constants
  IMPLICIT NONE
  PUBLIC
  INTEGER,  PARAMETER :: dp = KIND(1.0D0)
  INTEGER, PARAMETER :: dpc = KIND((1.0D0,1.0D0))
  REAL(dp), PUBLIC :: k_lambda, start_point
  REAL(DP), PUBLIC, ALLOCATABLE :: onebodyenergy(:)
  INTEGER , PUBLIC :: n_total
END MODULE constants

!Module to read in name/value pairs from a file, with each line of the form line 'name = value'

MODULE inifile
  IMPLICIT NONE
  PUBLIC
  INTEGER, PARAMETER :: ini_max_name_len = 128
  INTEGER, PARAMETER :: ini_max_string_len = 1024
  LOGICAL :: ini_fail_on_not_found = .FALSE.
  LOGICAL :: ini_echo_read = .FALSE.
  TYPE tnamevalue
     !no known way to make character string pointers..
     CHARACTER(ini_max_name_len)  :: name
     CHARACTER(ini_max_string_len):: value
  END TYPE tnamevalue

  TYPE tnamevalue_pointer
     TYPE(tnamevalue), POINTER :: p
  END TYPE tnamevalue_pointer

  TYPE tnamevaluelist
     INTEGER count
     INTEGER delta
     INTEGER capacity
     TYPE(tnamevalue_pointer), DIMENSION(:), POINTER :: items
  END TYPE tnamevaluelist

  TYPE tinifile
     LOGICAL slashcomments
     TYPE (tnamevaluelist) :: l, readvalues
  END TYPE tinifile

  TYPE(tinifile) :: defini

CONTAINS

  SUBROUTINE tnamevaluelist_init(l)
    TYPE (tnamevaluelist) :: l

    l%count = 0
    l%capacity = 0
    l%delta = 128
    NULLIFY(l%items)

  END SUBROUTINE tnamevaluelist_init

  SUBROUTINE tnamevaluelist_clear(l)
    TYPE (tnamevaluelist) :: l
    INTEGER i, status

    DO i=l%count,1,-1
       DEALLOCATE (l%items(i)%p, stat = status)
    END DO
    DEALLOCATE (l%items, stat = status)
    CALL tnamevaluelist_init(l)

  END SUBROUTINE tnamevaluelist_clear

  SUBROUTINE tnamevaluelist_valueof(l, aname, avalue)
    TYPE (tnamevaluelist) :: l
    CHARACTER(len=*), INTENT(in) :: aname
    CHARACTER(len=*) :: avalue
    INTEGER i

    DO i=1, l%count
       IF (l%items(i)%p%name == aname) THEN
          avalue = l%items(i)%p%value 
          RETURN
       END IF
    END DO
    avalue = ''

  END SUBROUTINE tnamevaluelist_valueof

  SUBROUTINE tnamevaluelist_add(l, aname, avalue)
    TYPE (tnamevaluelist) :: l
    CHARACTER(len=*), INTENT(in) :: aname, avalue

    IF (l%count == l%capacity) CALL tnamevaluelist_setcapacity(l, l%capacity + l%delta)
    l%count = l%count + 1
    ALLOCATE(l%items(l%count)%p)
    l%items(l%count)%p%name = aname
    l%items(l%count)%p%value = avalue

  END SUBROUTINE tnamevaluelist_add

  SUBROUTINE tnamevaluelist_setcapacity(l, c)
    TYPE (tnamevaluelist) :: l
    INTEGER c
    TYPE(tnamevalue_pointer), DIMENSION(:), POINTER :: tmpitems

    IF (l%count > 0) THEN
       IF (c < l%count) STOP 'tnamevaluelist_setcapacity: smaller than count'
       ALLOCATE(tmpitems(l%count))
       tmpitems = l%items(1:l%count)
       DEALLOCATE(l%items)
       ALLOCATE(l%items(c))
       l%items(1:l%count) = tmpitems
       DEALLOCATE(tmpitems)
    ELSE
       ALLOCATE(l%items(c))
    END IF
    l%capacity = c

  END SUBROUTINE tnamevaluelist_setcapacity

  SUBROUTINE tnamevaluelist_delete(l, i)
    TYPE (tnamevaluelist) :: l
    INTEGER, INTENT(in) :: i

    DEALLOCATE(l%items(i)%p)
    IF (l%count > 1) l%items(i:l%count-1) = l%items(i+1:l%count)
    l%count = l%count -1

  END SUBROUTINE tnamevaluelist_delete

  SUBROUTINE ini_namevalue_add(ini,ainline)
    TYPE(tinifile) :: ini
    CHARACTER (len=*), INTENT(in) :: ainline
    INTEGER eqpos, slashpos, lastpos
    CHARACTER (len=LEN(ainline)) :: aname, s, inline

    inline=TRIM(ADJUSTL(ainline))
    eqpos = SCAN(inline,'=')
    IF (eqpos/=0 .AND. inline(1:1)/='#' .AND. inline(1:7) /= 'comment' ) THEN
       aname = TRIM(inline(1:eqpos-1))
       s = ADJUSTL(inline(eqpos+1:)) 
       IF (ini%slashcomments) THEN
          slashpos=SCAN(s,'/')
          IF (slashpos /= 0) THEN
             s  = s(1:slashpos-1)
          END IF
       END IF
       lastpos=LEN_TRIM(s)
       IF (lastpos>1) THEN
          IF (s(1:1)=='''' .AND. s(lastpos:lastpos)=='''') THEN
             s = s(2:lastpos-1)
          END IF
       END IF
       CALL tnamevaluelist_add(ini%l, aname, s)
    END IF

  END SUBROUTINE ini_namevalue_add


  SUBROUTINE ini_open(filename, unit_id,  error, slash_comments)
    CHARACTER (len=*), INTENT(in) :: filename
    INTEGER, INTENT(in) :: unit_id
    LOGICAL, OPTIONAL, INTENT(out) :: error
    LOGICAL, OPTIONAL, INTENT(in) :: slash_comments
    LOGICAL aerror

    CALL tnamevaluelist_init(defini%l)
    CALL tnamevaluelist_init(defini%readvalues)
    IF (PRESENT(slash_comments)) THEN
       CALL ini_open_file(defini,filename,unit_id,aerror,slash_comments)
    ELSE
       CALL ini_open_file(defini,filename,unit_id,aerror)
    END IF
    IF (PRESENT(error)) THEN
       error = aerror
    ELSE
       IF (aerror) THEN
          WRITE (*,*) 'ini_open: error opening file ' // TRIM(filename)
          STOP
       END IF
    END IF

  END SUBROUTINE ini_open


  SUBROUTINE ini_open_file(ini, filename, unit_id,  error, slash_comments)
    TYPE(tinifile) :: ini
    CHARACTER (len=*), INTENT(in) :: filename
    INTEGER, INTENT(in) :: unit_id
    LOGICAL, INTENT(out) :: error
    LOGICAL, OPTIONAL, INTENT(in) :: slash_comments
    CHARACTER (len=120) :: inline

    CALL tnamevaluelist_init(ini%l)
    CALL tnamevaluelist_init(ini%readvalues)

    IF (PRESENT(slash_comments)) THEN
       ini%slashcomments = slash_comments
    ELSE
       ini%slashcomments = .FALSE.
    END IF
    OPEN(unit=unit_id,file=filename,form='formatted',status='old', err=500)
    DO 
       READ (unit_id,'(a)',END=400) inline
       IF (inline == 'end') EXIT;
       IF (inline /= '') CALL ini_namevalue_add(ini,inline) 
    END DO
400 CLOSE(unit_id)
    error=.FALSE.
    RETURN
500 error=.TRUE.

  END SUBROUTINE ini_open_file

  SUBROUTINE ini_open_fromlines(ini, lines, numlines, slash_comments)
    TYPE(tinifile) :: ini
    INTEGER, INTENT(in) :: numlines
    CHARACTER (len=*), DIMENSION(numlines), INTENT(in) :: lines
    LOGICAL, INTENT(in) :: slash_comments
    INTEGER i

    CALL tnamevaluelist_init(ini%l)
    ini%slashcomments = slash_comments
    DO i=1,numlines
       CALL ini_namevalue_add(ini,lines(i))
    END DO

  END  SUBROUTINE ini_open_fromlines

  SUBROUTINE ini_close

    CALL ini_close_file(defini)

  END SUBROUTINE ini_close


  SUBROUTINE ini_close_file(ini)
    TYPE(tinifile) :: ini

    CALL tnamevaluelist_clear(ini%l)
    CALL tnamevaluelist_clear(ini%readvalues)

  END  SUBROUTINE ini_close_file


  FUNCTION ini_read_string(key, notfoundfail) RESULT(avalue)
    CHARACTER (len=*), INTENT(in) :: key
    LOGICAL, OPTIONAL, INTENT(in) :: notfoundfail
    CHARACTER(len=ini_max_string_len) :: avalue

    IF (PRESENT(notfoundfail)) THEN
       avalue = ini_read_string_file(defini, key, notfoundfail)
    ELSE
       avalue = ini_read_string_file(defini, key)
    END IF

  END FUNCTION ini_read_string


  FUNCTION ini_read_string_file(ini, key, notfoundfail) RESULT(avalue)
    TYPE(tinifile) :: ini
    CHARACTER (len=*), INTENT(in) :: key
    LOGICAL, OPTIONAL, INTENT(in) :: notfoundfail
    CHARACTER(len=ini_max_string_len) :: avalue

    CALL tnamevaluelist_valueof(ini%l, key, avalue)
    IF (avalue/='') THEN
       CALL  tnamevaluelist_add(ini%readvalues, key, avalue)
       IF (ini_echo_read) WRITE (*,*) TRIM(key)//' = ',TRIM(avalue)
       RETURN
    END IF
    IF (ini_fail_on_not_found) THEN
       WRITE(*,*) 'key not found : '//key
       STOP
    END IF
    IF (PRESENT(notfoundfail)) THEN
       IF (notfoundfail) THEN
          WRITE(*,*) 'key not found : '//key
          STOP
       END IF
    END IF

  END FUNCTION ini_read_string_file


  FUNCTION ini_read_int(key, default)
    INTEGER, OPTIONAL, INTENT(in) :: default
    CHARACTER (len=*), INTENT(in) :: key
    INTEGER ini_read_int

    IF (PRESENT(default)) THEN
       ini_read_int = ini_read_int_file(defini, key, default)
    ELSE
       ini_read_int = ini_read_int_file(defini, key)
    END IF

  END FUNCTION ini_read_int


  FUNCTION ini_read_int_file(ini, key, default)
    TYPE(tinifile) :: ini
    INTEGER ini_read_int_file
    INTEGER, OPTIONAL, INTENT(in) :: default
    CHARACTER  (len=*), INTENT(in) :: key
    CHARACTER(len=ini_max_string_len) :: s

    s = ini_read_string_file(ini, key,.NOT. PRESENT(default))
    IF (s == '') THEN
       IF (.NOT. PRESENT(default)) THEN
          WRITE(*,*) 'no value for key: '//key
          STOP
       END IF
       ini_read_int_file = default
       WRITE (s,*) default
       CALL  tnamevaluelist_add(ini%readvalues, key, s)
    ELSE
       IF (VERIFY(TRIM(s),'-+0123456789') /= 0) GOTO 10
       READ (s,*, err = 10) ini_read_int_file
    END IF
    RETURN
10  WRITE (*,*) 'error reading integer for key: '//key
    STOP

  END FUNCTION ini_read_int_file

  FUNCTION ini_read_double(key, default)
    DOUBLE PRECISION, OPTIONAL, INTENT(in) :: default
    CHARACTER (len=*), INTENT(in) :: key
    DOUBLE PRECISION ini_read_double

    IF (PRESENT(default)) THEN
       ini_read_double = ini_read_double_file(defini, key, default)
    ELSE
       ini_read_double = ini_read_double_file(defini, key)
    END IF

  END FUNCTION ini_read_double


  FUNCTION ini_read_double_file(ini,key, default)
    TYPE(tinifile) :: ini
    DOUBLE PRECISION ini_read_double_file 
    DOUBLE PRECISION, OPTIONAL, INTENT(in) :: default
    CHARACTER (len=*), INTENT(in) :: key
    CHARACTER(len=ini_max_string_len) :: s

    s = ini_read_string_file(ini,key,.NOT. PRESENT(default))
    IF (s == '') THEN
       IF (.NOT. PRESENT(default)) THEN
          WRITE(*,*) 'no value for key: '//key
          STOP
       END IF
       ini_read_double_file = default
       WRITE (s,*) default
       CALL  tnamevaluelist_add(ini%readvalues, key, s)
    ELSE
       READ (s,*, err=10) ini_read_double_file
    END IF
    RETURN
10  WRITE (*,*) 'error reading double for key: '//key
    STOP

  END FUNCTION ini_read_double_file


  FUNCTION ini_read_real(key, default)
    REAL, OPTIONAL, INTENT(in) :: default
    CHARACTER (len=*), INTENT(in) :: key
    REAL ini_read_real

    IF (PRESENT(default)) THEN
       ini_read_real = ini_read_real_file(defini, key, default)
    ELSE
       ini_read_real = ini_read_real_file(defini, key)
    END IF

  END FUNCTION ini_read_real

  FUNCTION ini_read_real_file(ini,key, default)
    TYPE(tinifile) :: ini
    REAL ini_read_real_file 
    REAL, OPTIONAL, INTENT(in) :: default
    CHARACTER (len=*), INTENT(in) :: key
    CHARACTER(len=ini_max_string_len) :: s

    s = ini_read_string_file(ini,key,.NOT. PRESENT(default))
    IF (s == '') THEN
       IF (.NOT. PRESENT(default)) THEN
          WRITE(*,*) 'no value for key: '//key
          STOP
       END IF
       ini_read_real_file = default
       WRITE (s,*) default
       CALL  tnamevaluelist_add(ini%readvalues, key, s)
    ELSE
       READ (s,*, err=10) ini_read_real_file
    END IF
    RETURN
10  WRITE (*,*) 'error reading double for key: '//key
    STOP

  END FUNCTION ini_read_real_file


  FUNCTION ini_read_logical(key, default)
    LOGICAL, OPTIONAL, INTENT(in) :: default
    CHARACTER (len=*), INTENT(in) :: key
    LOGICAL ini_read_logical

    IF (PRESENT(default)) THEN
       ini_read_logical = ini_read_logical_file(defini, key, default)
    ELSE
       ini_read_logical = ini_read_logical_file(defini, key)
    END IF

  END FUNCTION ini_read_logical

  FUNCTION ini_read_logical_file(ini, key, default)
    TYPE(tinifile) :: ini
    LOGICAL ini_read_logical_file
    LOGICAL, OPTIONAL, INTENT(in) :: default
    CHARACTER  (len=*), INTENT(in) :: key

    CHARACTER(len=ini_max_string_len) :: s

    s = ini_read_string_file(ini,key,.NOT. PRESENT(default))
    IF (s == '') THEN
       IF (.NOT. PRESENT(default)) THEN
          WRITE(*,*) 'no value for key: '//key
          STOP
       END IF
       ini_read_logical_file = default
       WRITE (s,*) default
       CALL  tnamevaluelist_add(ini%readvalues, key, s)
    ELSE
       IF (VERIFY(TRIM(s),'10tf') /= 0) GOTO 10  
       READ (s,*, err = 10) ini_read_logical_file
    END IF

    RETURN

10  WRITE (*,*) 'error reading logical for key: '//key
    STOP
  END FUNCTION ini_read_logical_file


  SUBROUTINE ini_savereadvalues(afile,unit_id)
    CHARACTER(len=*)  :: afile
    INTEGER, INTENT(in) :: unit_id

    CALL ini_savereadvalues_file(defini, afile, unit_id)

  END SUBROUTINE ini_savereadvalues


  SUBROUTINE ini_savereadvalues_file(ini, afile, unit_id)
    TYPE(tinifile) :: ini
    CHARACTER(len=*), INTENT(in) :: afile
    INTEGER, INTENT(in) :: unit_id
    INTEGER i

    OPEN(unit=unit_id,file=afile,form='formatted',status='replace', err=500)

    DO i=1, ini%readvalues%count
       WRITE (unit_id,'(a)') TRIM(ini%readvalues%items(i)%p%name) // ' = ' &
            //TRIM(ini%readvalues%items(i)%p%value)

    END DO

    CLOSE(unit_id)
    RETURN

500 WRITE(*,*) 'ini_savereadvalues_file: error creating '//TRIM(afile)

  END SUBROUTINE ini_savereadvalues_file

END MODULE inifile



PROGRAM  pairing
  USE constants
  USE inifile
  IMPLICIT NONE
  INTEGER :: i, j, pq_confs
  REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: voper, toper, hoper
  REAL(dp) :: delta, g
  CHARACTER (LEN=120) :: infilename, outputfile
  LOGICAL :: fail

  infilename   = 'pairing.ini'
  CALL ini_open(infilename, 5, fail, .FALSE.)
  IF (fail) STOP 'Error opening parameter file, probably wrong file, use pairing.ini as name'
  !     open the output file
  outputfile = ini_read_string('output_run')
  OPEN(unit=6,file=outputfile)
  pq_confs = Ini_Read_Int('dimension_matrix')
  delta = Ini_Read_Double('value_delta');   g = Ini_Read_Double('value_g')
  k_lambda = Ini_Read_Double('value_lambda')
  start_point = Ini_Read_Double('value_startpoint')
  ALLOCATE( voper(pq_confs,pq_confs), toper(pq_confs,pq_confs)) 

  WRITE(6,'(7h delta=,f7.3,4h  g=,f7.3)')delta, g
  !     setup the interaction part of the matrix
  voper = 0.0d0;   toper = 0.0d0
  DO i=1,pq_confs
     voper(i,i) = -g
     DO  j=i+1, pq_confs
         voper(i,j) = -g*0.5d0
     ENDDO 
  ENDDO
  ! Setting up by hand zero elements
  voper(1,6) =0.0d0; voper(2,5) = 0.0d0;  voper(3,4) = 0.0d0; voper(2,5) = 0.0d0;   
  !     setup lower triangular part of h
  DO i=1,pq_confs-1
     DO  j=i+1, pq_confs
        voper(j,i)=voper(i,j)
     ENDDO
  ENDDO
  !     then the diagonal 3body unperturbed part
  !     two of the 2p2h excitations are degenerate
  DO i=1,3
     toper(i,i) = delta*i*2.d0
     toper(i+3,i+3)=delta*(i+2)*2.d0
  ENDDO
  hoper = toper+voper
  WRITE(6,*) 'Hamiltonian matrix'
  DO i = 1, pq_confs
     WRITE(6,'(6(4X,E12.6))') (hoper(j,i),j = 1, pq_confs) 
  ENDDO
  CALL flow_equations(voper,toper, pq_confs)
  DEALLOCATE(hoper,toper,voper)

END PROGRAM pairing

!
!           set up the interaction using
!           a flow renormalization group approach. 
!
SUBROUTINE flow_equations(vint, unperturbed, pq_confs)
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(in) ::  pq_confs
  INTEGER :: i, j
  REAL(kind=8), INTENT(in) :: vint(pq_confs,pq_confs), unperturbed(pq_confs,pq_confs)
  REAL(kind=8), DIMENSION(pq_confs,pq_confs) :: heff
  WRITE(6,*) 'RENORMALIZATION FLOW METHOD STARTS HERE'
  heff = 0.d0
  CALL vsrg(pq_confs,vint,heff,unperturbed)
  WRITE(6,*) 'Hamiltonian matrix  with the renormalization flow method'
  DO i = 1, pq_confs
     WRITE(6,'(6(4X,E12.6))') (heff(j,i)+unperturbed(j,i),j = 1, pq_confs) 
  ENDDO

END SUBROUTINE flow_equations
!
!  
!
SUBROUTINE vsrg(nconfs,vzz,heff, unperturbed)
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(in) :: nconfs
  REAL(kind=8), DIMENSION(nconfs,nconfs), INTENT(in) :: vzz, unperturbed
  REAL(kind=8), DIMENSION(nconfs,nconfs), INTENT(inout) :: heff
  INTEGER(kind=4) :: i, j, n_ode, ij, iflag, i1, j1
  INTEGER(kind=4), ALLOCATABLE, DIMENSION(:) :: iwork
  REAL(dp), ALLOCATABLE, DIMENSION(:) :: rg_vec, work
  REAL(dp) :: relerr, abserr, lambda, END
  INTERFACE
     SUBROUTINE derivative(lambda, v, dv)
       USE constants
       INTEGER ::  ij, i, j, i1, j2, kk, k
       REAL(dp) :: sum, v(:), dv(:), k1, k2, p
       REAL(dp), ALLOCATABLE :: vij(:,:)
     END SUBROUTINE derivative
  END INTERFACE

  relerr = 1.0e-8_dp; abserr =  1.0e-8_dp; lambda = k_lambda; END = start_point
  ! dimension of vectors and matrices
  n_total = nconfs; n_ode = n_total*(n_total+1)/2
  ! total dim used by derivatives function
  ALLOCATE(rg_vec(n_ode)); ALLOCATE( work(100+21*n_ode)); ALLOCATE (iwork(5))
  ALLOCATE ( onebodyenergy(n_total))
  rg_vec = 0.0_dp; work = 0.0_dp;  iwork = 0
  !  transform the matrix -v- to a one-dim vector of dim ntot*(ntot+1)/2
  !  which is also the number of ordinary differential equations
  ij = 0  
  DO i = 1, n_total
     onebodyenergy(i) = unperturbed(i,i)
     DO j = i, n_total
        ij = ij + 1
        rg_vec(ij) = vzz(j,i)
     ENDDO
  ENDDO
  iflag = 1
  !   oscillator energy of max space, nlmas typically around 200-300 (2n+l)
  CALL ode(derivative,n_ode,rg_vec,END,lambda,relerr,abserr,iflag,work,iwork)
  WRITE(6,*) 'iflag=', iflag
  IF ( iflag /= 2) WRITE(6,*) 'error in ode, iflag not equal 2'
  !  now transform back and get final effective interaction
  ij = 0  
  DO i = 1, n_total
     DO j = i, n_total
        ij = ij + 1
        heff(j,i) = rg_vec(ij)
        heff(i,j) = rg_vec(ij)
     ENDDO
  ENDDO
  DEALLOCATE(rg_vec); DEALLOCATE(work); DEALLOCATE(iwork)

  DEALLOCATE(onebodyenergy)

END SUBROUTINE vsrg

SUBROUTINE derivative(lambda, v, dv)
  USE constants
  IMPLICIT NONE
  INTEGER ::  ij, i, j, i1, j2, kk, k, i2
  REAL(dp) :: sum, lambda, v(:), dv(:), k1, k2, p
  REAL(dp), ALLOCATABLE :: vij(:,:)

  ALLOCATE(vij(n_total,n_total)) 
  ij = 0
  DO j = 1, n_total
     DO i = j, n_total
        ij = ij + 1
        vij(i,j) = v(ij)
        vij(j,i) = v(ij)
     ENDDO
  ENDDO
  ij = 0
  DO i = 1, n_total
     k1 = onebodyenergy(i) 
     DO j = i, n_total
        k2 = onebodyenergy(j) 
        ij = ij + 1
        sum = 0_dp
        DO k = 1, n_total
           p = onebodyenergy(k) 
           sum = sum +  (k1+k2-2.0_dp*p)*vij(j,k)*vij(k,i) 
        ENDDO
        dv(ij) =  sum - (k2-k1)*(k2-k1)*vij(j,i)
     ENDDO
  ENDDO
  dv = -2.0_dp*dv/(lambda**3)
  DEALLOCATE(vij)

END SUBROUTINE derivative




!
!  This function is taken from www.netlib.org
!
SUBROUTINE ode(f,neqn,y,t,tout,relerr,abserr,iflag,work,iwork)
  IMPLICIT REAL*8(a-h,o-z)
  !
  !   double precision subroutine ode integrates a system of neqn
  !   first order ordinary differential equations of the form:
  !             dy(i)/dt = f(t,y(1),y(2),...,y(neqn))
  !             y(i) given at  t .
  !   the subroutine integrates from  t  to  tout .  on return the
  !   parameters in the call list are set for continuing the integration.
  !   the user has only to define a new value  tout  and call  ode  again.
  !
  !   the differential equations are actually solved by a suite of codes
  !   de ,  step , and  intrp .  ode  allocates virtual storage in the
  !   arrays  work  and  iwork  and calls  de .  de  is a supervisor which
  !   directs the solution.  it calls on the routines  step  and  intrp
  !   to advance the integration and to interpolate at output points.
  !   step  uses a modified divided difference form of the adams pece
  !   formulas and local extrapolation.  it adjusts the order and step
  !   size to control the local error per unit step in a generalized
  !   sense.  normally each call to  step  advances the solution one step
  !   in the direction of  tout .  for reasons of efficiency  de
  !   integrates beyond  tout  internally, though never beyond
  !   t+10*(tout-t), and calls  intrp  to interpolate the solution at
  !   tout .  an option is provided to stop the integration at  tout  but
  !   it should be used only if it is impossible to continue the
  !   integration beyond  tout .
  !
  !   this code is completely explained and documented in the text,
  !   computer solution of ordinary differential equations:  the initial
  !   value problem  by l. f. shampine and m. k. gordon.
  !
  !   the parameters represent:
  !      f -- double precision subroutine f(t,y,yp) to evaluate
  !                derivatives yp(i)=dy(i)/dt
  !      neqn -- number of equations to be integrated (integer*4)
  !      y(*) -- solution vector at t                 (real*8)
  !      t -- independent variable                    (real*8)
  !      tout -- point at which solution is desired   (real*8)
  !      relerr,abserr -- relative and absolute error tolerances for local
  !           error test (real*8).  at each step the code requires
  !             dabs(local error) .le. dabs(y)*relerr + abserr
  !           for each component of the local error and solution vectors
  !      iflag -- indicates status of integration     (integer*4)
  !      work(*)  (real*8)  -- arrays to hold information internal to
  !      iwork(*) (integer*4)    which is necessary for subsequent calls
  !
  !   first call to ode --
  !
  !   the user must provide storage in his calling program for the arrays
  !   in the call list,
  !      y(neqn), work(100+21*neqn), iwork(5),
  !   declare  f  in an external statement, supply the double precision
  !   subroutine f(t,y,yp)  to evaluate
  !      dy(i)/dt = yp(i) = f(t,y(1),y(2),...,y(neqn))
  !   and initialize the parameters:
  !      neqn -- number of equations to be integrated
  !      y(*) -- vector of initial conditions
  !      t -- starting point of integration
  !      tout -- point at which solution is desired
  !      relerr,abserr -- relative and absolute local error tolerances
  !      iflag -- +1,-1.  indicator to initialize the code.  normal input
  !           is +1.  the user should set iflag=-1 only if it is
  !           impossible to continue the integration beyond  tout .
  !   all parameters except  f ,  neqn  and  tout  may be altered by the
  !   code on output so must be variables in the calling program.
  !
  !   output from  ode  --
  !
  !      neqn -- unchanged
  !      y(*) -- solution at  t
  !      t -- last point reached in integration.  normal return has
  !           t = tout .
  !      tout -- unchanged
  !      relerr,abserr -- normal return has tolerances unchanged.  iflag=3
  !           signals tolerances increased
  !      iflag = 2 -- normal return.  integration reached  tout
  !            = 3 -- integration did not reach  tout  because error
  !                   tolerances too small.  relerr ,  abserr  increased
  !                   appropriately for continuing
  !            = 4 -- integration did not reach  tout  because more than
  !                   500 steps needed
  !            = 5 -- integration did not reach  tout  because equations
  !                   appear to be stiff
  !            = 6 -- invalid input parameters (fatal error)
  !           the value of  iflag  is returned negative when the input
  !           value is negative and the integration does not reach  tout ,
  !           i.e., -3, -4, -5.
  !      work(*),iwork(*) -- information generally of no interest to the
  !           user but necessary for subsequent calls.
  !
  !   subsequent calls to  ode --
  !
  !   subroutine  ode  returns with all information needed to continue
  !   the integration.  if the integration reached  tout , the user need
  !   only define a new  tout  and call again.  if the integration did not
  !   reach  tout  and the user wants to continue, he just calls again.
  !   the output value of  iflag  is the appropriate input value for
  !   subsequent calls.  the only situation in which it should be altered
  !   is to stop the integration internally at the new  tout , i.e.,
  !   change output  iflag=2  to input  iflag=-2 .  error tolerances may
  !   be changed by the user before continuing.  all other parameters must
  !   remain unchanged.
  !
  !***********************************************************************
  !*  subroutines  de  and  step  contain machine dependent constants. *
  !*  be sure they are set before using  ode .                          *
  !***********************************************************************
  !
  LOGICAL start,phase1,nornd
  DIMENSION y(neqn),work(1),iwork(5)
  DATA ialpha,ibeta,isig,iv,iw,ig,iphase,ipsi,ix,ih,ihold,istart,  &
       itold,idelsn/1,13,25,38,50,62,75,76,88,89,90,91,92,93/

  INTERFACE
     SUBROUTINE f(lam, v, dv)
       USE constants
       INTEGER ::  ij, i, j, ii, jj, kk, k
       REAL(DP) :: ans, lam, v(:), dv(:), k1, k2, p
       REAL(DP), ALLOCATABLE :: vij(:,:)
     END SUBROUTINE f
  END INTERFACE

  iyy = 100
  iwt = iyy + neqn
  ip = iwt + neqn
  iyp = ip + neqn
  iypout = iyp + neqn
  iphi = iypout + neqn
  IF(iabs(iflag) /= 1) THEN
     start = work(istart) > 0.0d0
     phase1 = work(iphase) > 0.0d0
     nornd = iwork(2) /= -1
  ELSE
     CALL de(f,neqn,y,t,tout,relerr,abserr,iflag,work(iyy), & 
          work(iwt),work(ip),work(iyp),work(iypout),work(iphi),  &
          work(ialpha),work(ibeta),work(isig),work(iv),work(iw),work(ig), &
          phase1,work(ipsi),work(ix),work(ih),work(ihold),start,  &
          work(itold),work(idelsn),iwork(1),nornd,iwork(3),iwork(4),  &
          iwork(5))
  ENDIF
  work(istart) = -1.0d0
  IF(start) work(istart) = 1.0d0
  work(iphase) = -1.0d0
  IF(phase1) work(iphase) = 1.0d0
  iwork(2) = -1
  IF(nornd) iwork(2) = 1

END SUBROUTINE  ode


SUBROUTINE de(f,neqn,y,t,tout,relerr,abserr,iflag,  &
     yy,wt,p,yp,ypout,phi,alpha,beta,sig,v,w,g,phase1,psi,x,h,hold, &
     start,told,delsgn,ns,nornd,k,kold,isnold)
  IMPLICIT REAL*8(a-h,o-z)
  !
  !   ode  merely allocates storage for  de  to relieve the user of the
  !   inconvenience of a long call list.  consequently  de  is used as
  !   described in the comments for  ode .
  !
  !   this code is completely explained and documented in the text,
  !   computer solution of ordinary differential equations:  the initial
  !   value problem  by l. f. shampine and m. k. gordon.
  !
  LOGICAL stiff,crash,start,phase1,nornd
  DIMENSION y(neqn),yy(neqn),wt(neqn),phi(neqn,16),p(neqn),yp(neqn),  &
       ypout(neqn),psi(12),alpha(12),beta(12),sig(13),v(12),w(12),g(13)

  INTERFACE
     SUBROUTINE f(lam, v, dv)
       USE constants
       INTEGER ::  ij, i, j, ii, jj, kk, k
       REAL(DP) :: ans, lam, v(:), dv(:), k1, k2, p
       REAL(DP), ALLOCATABLE :: vij(:,:)
     END SUBROUTINE f
  END INTERFACE
  !
  !***********************************************************************
  !*  the only machine dependent constant is based on the machine unit   *
  !*  roundoff error  u  which is the smallest positive number such that *
  !*  1.0+u .gt. 1.0 .  u  must be calculated and  fouru=4.0*u  inserted *
  !*  in the following data statement before using  de .  the routine    *
  !*  machin  calculates  u .  fouru  and  twou=2.0*u  must also be      *
  !*  inserted in subroutine  step  before calling  de .                 *
  !     data fouru/.888d-15/                                              ***
  !***********************************************************************
  !
  !   the constant  maxnum  is the maximum number of steps allowed in one
  !   call to  de .  the user may change this limit by altering the
  !   following statement
  DATA maxnum/5000000/


  fouru = 4.0 * d1mach(4)                                          ! ***
  IF(neqn <  1) go to 10
  IF(t ==  tout) go to 10
  IF(relerr  <  0.0d0  .OR.  abserr .LT. 0.0d0) go to 10
  eps = dmax1(relerr,abserr)
  IF(eps  <=  0.0d0) go to 10
  IF(iflag ==  0) go to 10
  isn = isign(1,iflag)
  iflag = iabs(iflag)
  IF(iflag ==  1) go to 20
  IF(t /=  told) go to 10
  IF(iflag .GE. 2  .AND.  iflag .LE. 5) go to 20
10 iflag = 6
  RETURN
  !
  !   on each call set interval of integration and counter for number of
  !   steps.  adjust input error tolerances to define weight vector for
  !   subroutine  step
  !
20 del = tout - t
  absdel = dabs(del)
  tend = t + 10.0d0*del
  IF(isn .LT. 0) tend = tout
  nostep = 0
  kle4 = 0
  stiff = .FALSE.
  releps = relerr/eps
  abseps = abserr/eps
  IF(iflag .EQ. 1) go to 30
  IF(isnold .LT. 0) go to 30
  IF(delsgn*del .GT. 0.0d0) go to 50
  !
  !   on start and restart also set work variables x and yy(*), store the
  !   direction of integration and initialize the step size
  !
30 start = .TRUE.
  x = t
  DO l = 1,neqn
     yy(l) = y(l)
  ENDDO
  delsgn = dsign(1.0d0,del)
  h = dsign(dmax1(dabs(tout-x),fouru*dabs(x)),tout-x)
  !
  !   if already past output point, interpolate and return
  !
50 IF(dabs(x-t) .LT. absdel) go to 60
  CALL intrp(x,yy,tout,y,ypout,neqn,kold,phi,psi)
  iflag = 2
  t = tout
  told = t
  isnold = isn
  RETURN
  !
  !   if cannot go past output point and sufficiently close,
  !   extrapolate and return
  !
60 IF(isn .GT. 0  .OR.  dabs(tout-x) .GE. fouru*dabs(x)) go to 80
  h = tout - x
  CALL f(x,yy,yp)
  DO l = 1,neqn
     y(l) = yy(l) + h*yp(l)
  ENDDO
  iflag = 2
  t = tout
  told = t
  isnold = isn
  RETURN
  !
  !   test for too many steps
  !
80 IF(nostep .LT. maxnum) go to 100
  iflag = isn*4
  IF(stiff) iflag = isn*5
  DO l = 1,neqn
     y(l) = yy(l)
  ENDDO
  t = x
  told = t
  isnold = 1
  RETURN
  !
  !   limit step size, set weight vector and take a step
  !
100 h = dsign(dmin1(dabs(h),dabs(tend-x)),h)
  DO l = 1,neqn
     wt(l) = releps*dabs(yy(l)) + abseps
  ENDDO
  CALL step(x,yy,f,neqn,h,eps,wt,start,  &
       hold,k,kold,crash,phi,p,yp,psi,         &
       alpha,beta,sig,v,w,g,phase1,ns,nornd)
  !
  !   test for tolerances too small
  !
  IF(.NOT.crash) go to 130
  iflag = isn*3
  relerr = eps*releps
  abserr = eps*abseps
  DO l = 1,neqn
     y(l) = yy(l)
  ENDDO
  t = x
  told = t
  isnold = 1
  RETURN
  !
  !   augment counter on number of steps and test for stiffness
  !
130 nostep = nostep + 1
  kle4 = kle4 + 1
  IF(kold .GT. 4) kle4 = 0
  IF(kle4 .GE. 50) stiff = .TRUE.
  go to 50

END SUBROUTINE  de



  SUBROUTINE step(x,y,f,neqn,h,eps,wt,start,  &
       hold,k,kold,crash,phi,p,yp,psi,  &
       alpha,beta,sig,v,w,g,phase1,ns,nornd)
    IMPLICIT REAL*8(a-h,o-z)
    !
    !   double precision subroutine  step
    !   integrates a system of first order ordinary
    !   differential equations one step, normally from x to x+h, using a
    !   modified divided difference form of the adams pece formulas.  local
    !   extrapolation is used to improve absolute stability and accuracy.
    !   the code adjusts its order and step size to control the local error
    !   per unit step in a generalized sense.  special devices are included
    !   to control roundoff error and to detect when the user is requesting
    !   too much accuracy.
    !
    !   this code is completely explained and documented in the text,
    !   computer solution of ordinary differential equations:  the initial
    !   value problem  by l. f. shampine and m. k. gordon.
    !
    !
    !   the parameters represent:
    !      x -- independent variable             (real*8)
    !      y(*) -- solution vector at x          (real*8)
    !      yp(*) -- derivative of solution vector at  x  after successful
    !           step                             (real*8)
    !      neqn -- number of equations to be integrated (integer*4)
    !      h -- appropriate step size for next step.  normally determined by
    !           code                             (real*8)
    !      eps -- local error tolerance.  must be variable  (real*8)
    !      wt(*) -- vector of weights for error criterion   (real*8)
    !      start -- logical variable set .true. for first step,  .false.
    !           otherwise                        (logical*4)
    !      hold -- step size used for last successful step  (real*8)
    !      k -- appropriate order for next step (determined by code)
    !      kold -- order used for last successful step
    !      crash -- logical variable set .true. when no step can be taken,
    !           .false. otherwise.
    !   the arrays  phi, psi  are required for the interpolation subroutine
    !   intrp.  the array p is internal to the code.  all are real*8
    !
    !   input to  step
    !
    !      first call --
    !
    !   the user must provide storage in his driver program for all arrays
    !   in the call list, namely
    !
    !     dimension y(neqn),wt(neqn),phi(neqn,16),p(neqn),yp(neqn),psi(12)
    !
    !   the user must also declare  start  and  crash  logical variables
    !   and  f  an external subroutine, supply the subroutine  f(x,y,yp)
    !   to evaluate
    !      dy(i)/dx = yp(i) = f(x,y(1),y(2),...,y(neqn))
    !   and initialize only the following parameters:
    !      x -- initial value of the independent variable
    !      y(*) -- vector of initial values of dependent variables
    !      neqn -- number of equations to be integrated
    !      h -- nominal step size indicating direction of integration
    !           and maximum size of step.  must be variable
    !      eps -- local error tolerance per step.  must be variable
    !      wt(*) -- vector of non-zero weights for error criterion
    !      start -- .true.
    !
    !   step  requires the l2 norm of the vector with components
    !   local error(l)/wt(l)  be less than  eps  for a successful step.  the
    !   array  wt  allows the user to specify an error test appropriate
    !   for his problem.  for example,
    !      wt(l) = 1.0  specifies absolute error,
    !            = dabs(y(l))  error relative to the most recent value of
    !                 the l-th component of the solution,
    !            = dabs(yp(l))  error relative to the most recent value of
    !                 the l-th component of the derivative,
    !            = dmax1(wt(l),dabs(y(l)))  error relative to the largest
    !                 magnitude of l-th component obtained so far,
    !            = dabs(y(l))*relerr/eps + abserr/eps  specifies a mixed
    !                 relative-absolute test where  relerr  is relative
    !                 error,  abserr  is absolute error and  eps =
    !                 dmax1(relerr,abserr) .
    !
    !      subsequent calls --
    !
    !   subroutine  step  is designed so that all information needed to
    !   continue the integration, including the step size  h  and the order
    !   k , is returned with each step.  with the exception of the step
    !   size, the error tolerance, and the weights, none of the parameters
    !   should be altered.  the array  wt  must be updated after each step
    !   to maintain relative error tests like those above.  normally the
    !   integration is continued just beyond the desired endpoint and the
    !   solution interpolated there with subroutine  intrp .  if it is
    !   impossible to integrate beyond the endpoint, the step size may be
    !   reduced to hit the endpoint since the code will not take a step
    !   larger than the  h  input.  changing the direction of integration,
    !   i.e., the sign of  h , requires the user set  start = .true. before
    !   calling  step  again.  this is the only situation in which  start
    !   should be altered.
    !
    !   output from  step
    !
    !      successful step --
    !
    !   the subroutine returns after each successful step with  start  and
    !   crash  set .false. .  x  represents the independent variable
    !   advanced one step of length  hold  from its value on input and  y
    !   the solution vector at the new value of  x .  all other parameters
    !   represent information corresponding to the new  x  needed to
    !   continue the integration.
    !
    !      unsuccessful step --
    !
    !   when the error tolerance is too small for the machine precision,
    !   the subroutine returns without taking a step and  crash = .true. .
    !   an appropriate step size and error tolerance for continuing are
    !   estimated and all other information is restored as upon input
    !   before returning.  to continue with the larger tolerance, the user
    !   just calls the code again.  a restart is neither required nor
    !   desirable.
    !
    LOGICAL start,crash,phase1,nornd
    DIMENSION y(neqn),wt(neqn),phi(neqn,16),p(neqn),yp(neqn),psi(12)
    DIMENSION alpha(12),beta(12),sig(13),w(12),v(12),g(13),  &
         gstr(13),two(13)

    INTERFACE
       SUBROUTINE f(lam, v, dv)
         USE constants
         INTEGER ::  ij, i, j, ii, jj, kk, k
         REAL(DP) :: ans, lam, v(:), dv(:), k1, k2, p
         REAL(DP), ALLOCATABLE :: vij(:,:)
       END SUBROUTINE f
    END INTERFACE
    !***********************************************************************
    !*  the only machine dependent constants are based on the machine unit *
    !*  roundoff error  u  which is the smallest positive number such that *
    !*  1.0+u .gt. 1.0  .  the user must calculate  u  and insert          *
    !*  twou=2.0*u  and  fouru=4.0*u  in the data statement before calling *
    !*  the code.  the routine  machin  calculates  u .                    *
    !     data twou,fouru/.444d-15,.888d-15/                                ***
    !***********************************************************************
    DATA two/2.0d0,4.0d0,8.0d0,16.0d0,32.0d0,64.0d0,128.0d0,256.0d0,  &
         512.0d0,1024.0d0,2048.0d0,4096.0d0,8192.0d0/
    DATA gstr/0.500d0,0.0833d0,0.0417d0,0.0264d0,0.0188d0,0.0143d0, &
         0.0114d0,0.00936d0,0.00789d0,0.00679d0,0.00592d0,0.00524d0,  &
         0.00468d0/

    twou = 2.0 * d1mach(4)                                           ! ***
    fouru = 2.0 * twou                                               ! ***
    !
    !   if step size is too small, determine an acceptable one
    !
    crash = .TRUE.
      IF(dabs(h) .GE. fouru*dabs(x)) go to 5
      h = dsign(fouru*dabs(x),h)
      RETURN
    5 p5eps = 0.5d0*eps
!
!   if error tolerance is too small, increase it to an acceptable value
!
      round = 0.0d0
      DO l = 1,neqn
      round = round + (y(l)/wt(l))**2
      ENDDO
      round = twou*dsqrt(round)
      IF(p5eps .GE. round) go to 15
      eps = 2.0*round*(1.0d0 + fouru)
      RETURN
   15 crash = .FALSE.
      g(1)=1.0d0
      g(2)=0.5d0
      sig(1)=1.0d0
      IF(.NOT.start) go to 99
!
!   initialize.  compute appropriate step size for first step
!
      CALL f(x,y,yp)
      sum = 0.0d0
      DO l = 1,neqn
        phi(l,1) = yp(l)
        phi(l,2) = 0.0d0
        sum = sum + (yp(l)/wt(l))**2
      ENDDO
      sum = dsqrt(sum)
      absh = dabs(h)
      IF(eps <  16.0d0*sum*h*h) absh = 0.25d0*dsqrt(eps/sum)
      h = dsign(dmax1(absh,fouru*dabs(x)),h)
      hold = 0.0d0
      k = 1
      kold = 0
      start = .FALSE.
      phase1 = .TRUE.
      nornd = .TRUE.
      IF(p5eps .GT. 100.0d0*round) go to 99
      nornd = .FALSE.
      DO l = 1,neqn
       phi(l,15) = 0.0d0
      ENDDO
   99 ifail = 0
!       ***     end block 0     ***
!
!       ***     begin block 1     ***
!   compute coefficients of formulas for this step.  avoid computing
!   those quantities not changed when step size is not changed.
!                   ***
!
  100 kp1 = k+1
      kp2 = k+2
      km1 = k-1
      km2 = k-2
!
!   ns is the number of steps taken with size h, including the current
!   one.  when k.lt.ns, no coefficients change
!
      IF(h .NE. hold) ns = 0
      IF(ns.LE.kold)   ns=ns+1
      nsp1 = ns+1
      IF (k .LT. ns) go to 199
!
!   compute those components of alpha(*),beta(*),psi(*),sig(*) which
!   are changed
!
      beta(ns) = 1.0d0
      realns = ns
      alpha(ns) = 1.0d0/realns
      temp1 = h*realns
      sig(nsp1) = 1.0d0
      IF(k <  nsp1) go to 110
      DO i = nsp1,k
        im1 = i-1
        temp2 = psi(im1)
        psi(im1) = temp1
        beta(i) = beta(im1)*psi(im1)/temp2
        temp1 = temp2 + h
        alpha(i) = h/temp1
        reali = i
        sig(i+1) = reali*alpha(i)*sig(i)
      ENDDO
  110 psi(k) = temp1
!
!   compute coefficients g(*)
!
!   initialize v(*) and set w(*).  g(2) is set in data statement
!
      IF(ns >  1) go to 120
      DO iq = 1,k
        temp3 = iq*(iq+1)
        v(iq) = 1.0d0/temp3
        w(iq) = v(iq)
      ENDDO
      go to 140
!
!   if order was raised, update diagonal part of v(*)
!
  120 IF(k <=  kold) go to 130
      temp4 = k*kp1
      v(k) = 1.0d0/temp4
      nsm2 = ns-2
      IF(nsm2 .LT. 1) go to 130
      DO j = 1,nsm2
        i = k-j
       v(i) = v(i) - alpha(j+1)*v(i+1)
      ENDDO
!
!   update v(*) and set w(*)
!
  130 limit1 = kp1 - ns
      temp5 = alpha(ns)
      DO iq = 1,limit1
        v(iq) = v(iq) - temp5*v(iq+1)
        w(iq) = v(iq)
      ENDDO
      g(nsp1) = w(1)
!
!   compute the g(*) in the work vector w(*)
!
  140 nsp2 = ns + 2
      IF(kp1 .LT. nsp2) go to 199
      DO 150 i = nsp2,kp1
        limit2 = kp2 - i
        temp6 = alpha(i-1)
        DO 145 iq = 1,limit2
  145     w(iq) = w(iq) - temp6*w(iq+1)
  150   g(i) = w(1)
  199   CONTINUE
!       ***     end block 1     ***
!
!       ***     begin block 2     ***
!   predict a solution p(*), evaluate derivatives using predicted
!   solution, estimate local error at order k and errors at orders k,
!   k-1, k-2 as if constant step size were used.
!                   ***
!
!   change phi to phi star
!
      IF(k .LT. nsp1) go to 215
      DO 210 i = nsp1,k
        temp1 = beta(i)
        DO 205 l = 1,neqn
  205     phi(l,i) = temp1*phi(l,i)
  210   CONTINUE
!
!   predict solution and differences
!
  215 DO 220 l = 1,neqn
        phi(l,kp2) = phi(l,kp1)
        phi(l,kp1) = 0.0d0
  220   p(l) = 0.0d0
      DO 230 j = 1,k
        i = kp1 - j
        ip1 = i+1
        temp2 = g(i)
        DO 225 l = 1,neqn
          p(l) = p(l) + temp2*phi(l,i)
  225     phi(l,i) = phi(l,i) + phi(l,ip1)
  230   CONTINUE
      IF(nornd) go to 240
      DO 235 l = 1,neqn
        tau = h*p(l) - phi(l,15)
        p(l) = y(l) + tau
  235   phi(l,16) = (p(l) - y(l)) - tau
      go to 250
  240 DO 245 l = 1,neqn
  245   p(l) = y(l) + h*p(l)
  250 xold = x
      x = x + h
      absh = dabs(h)
      CALL f(x,p,yp)
!
!   estimate errors at orders k,k-1,k-2
!
      erkm2 = 0.0d0
      erkm1 = 0.0d0
      erk = 0.0d0
      DO 265 l = 1,neqn
        temp3 = 1.0d0/wt(l)
        temp4 = yp(l) - phi(l,1)
        IF(km2)265,260,255
  255   erkm2 = erkm2 + ((phi(l,km1)+temp4)*temp3)**2
  260   erkm1 = erkm1 + ((phi(l,k)+temp4)*temp3)**2
  265   erk = erk + (temp4*temp3)**2
      IF(km2)280,275,270
  270 erkm2 = absh*sig(km1)*gstr(km2)*dsqrt(erkm2)
  275 erkm1 = absh*sig(k)*gstr(km1)*dsqrt(erkm1)
  280 temp5 = absh*dsqrt(erk)
      err = temp5*(g(k)-g(kp1))
      erk = temp5*sig(kp1)*gstr(k)
      knew = k
!
!   test if order should be lowered
!
      IF(km2)299,290,285
  285 IF(dmax1(erkm1,erkm2) .LE. erk) knew = km1
      go to 299
  290 IF(erkm1 .LE. 0.5d0*erk) knew = km1
!
!   test if step successful
!
  299 IF(err .LE. eps) go to 400
!       ***     end block 2     ***
!
!       ***     begin block 3     ***
!   the step is unsuccessful.  restore  x, phi(*,*), psi(*) .
!   if third consecutive failure, set order to one.  if step fails more
!   than three times, consider an optimal step size.  double error
!   tolerance and return if estimated step size is too small for machine
!   precision.
!                   ***

!   restore x, phi(*,*) and psi(*)

      phase1 = .FALSE.
      x = xold
      DO 310 i = 1,k
        temp1 = 1.0d0/beta(i)
        ip1 = i+1
        DO 305 l = 1,neqn
  305     phi(l,i) = temp1*(phi(l,i) - phi(l,ip1))
  310   CONTINUE
      IF(k .LT. 2) go to 320
      DO 315 i = 2,k
  315   psi(i-1) = psi(i) - h
!
!   on third failure, set order to one.  thereafter, use optimal step
!   size
!
  320 ifail = ifail + 1
      temp2 = 0.5d0
      IF(ifail - 3) 335,330,325
  325 IF(p5eps .LT. 0.25d0*erk) temp2 = dsqrt(p5eps/erk)
  330 knew = 1
  335 h = temp2*h
      k = knew
      IF(dabs(h) .GE. fouru*dabs(x)) go to 340
      crash = .TRUE.
      h = dsign(fouru*dabs(x),h)
      eps = eps + eps
      RETURN
  340 go to 100
!       ***     end block 3     ***
!
!       ***     begin block 4     ***
!   the step is successful.  correct the predicted solution, evaluate
!   the derivatives using the corrected solution and update the
!   differences.  determine best order and step size for next step.
!                   ***
  400 kold = k
      hold = h
!
!   correct and evaluate
!
      temp1 = h*g(kp1)
      IF(nornd) go to 410
      DO 405 l = 1,neqn
        rho = temp1*(yp(l) - phi(l,1)) - phi(l,16)
        y(l) = p(l) + rho
  405   phi(l,15) = (y(l) - p(l)) - rho
      go to 420
  410 DO 415 l = 1,neqn
  415   y(l) = p(l) + temp1*(yp(l) - phi(l,1))
  420 CALL f(x,y,yp)
!
!   update differences for next step
!
      DO 425 l = 1,neqn
        phi(l,kp1) = yp(l) - phi(l,1)
  425   phi(l,kp2) = phi(l,kp1) - phi(l,kp2)
      DO 435 i = 1,k
        DO 430 l = 1,neqn
  430     phi(l,i) = phi(l,i) + phi(l,kp1)
  435   CONTINUE
!
!   estimate error at order k+1 unless:
!     in first phase when always raise order,
!     already decided to lower order,
!     step size not constant so estimate unreliable
!
      erkp1 = 0.0d0
      IF(knew .EQ. km1  .OR.  k .EQ. 12) phase1 = .FALSE.
      IF(phase1) go to 450
      IF(knew .EQ. km1) go to 455
      IF(kp1 .GT. ns) go to 460
      DO 440 l = 1,neqn
  440   erkp1 = erkp1 + (phi(l,kp2)/wt(l))**2
      erkp1 = absh*gstr(kp1)*dsqrt(erkp1)
!
!   using estimated error at order k+1, determine appropriate order
!   for next step
!
      IF(k .GT. 1) go to 445
      IF(erkp1 .GE. 0.5d0*erk) go to 460
      go to 450
  445 IF(erkm1 .LE. dmin1(erk,erkp1)) go to 455
      IF(erkp1 .GE. erk  .OR.  k .EQ. 12) go to 460
!
!   here erkp1 .lt. erk .lt. dmax1(erkm1,erkm2) else order would have
!   been lowered in block 2.  thus order is to be raised
!
!   raise order
!
  450 k = kp1
      erk = erkp1
      go to 460
!
!   lower order
!
  455 k = km1
      erk = erkm1
!
!   with new order determine appropriate step size for next step
!
  460 hnew = h + h
      IF(phase1) go to 465
      IF(p5eps .GE. erk*two(k+1)) go to 465
      hnew = h
      IF(p5eps .GE. erk) go to 465
      temp2 = k+1
      r = (p5eps/erk)**(1.0d0/temp2)
      hnew = absh*dmax1(0.5d0,dmin1(0.9d0,r))
      hnew = dsign(dmax1(hnew,fouru*dabs(x)),h)
  465 h = hnew
      RETURN
!       ***     end block 4     ***

  END SUBROUTINE  step


  SUBROUTINE intrp(x,y,xout,yout,ypout,neqn,kold,phi,psi)
    IMPLICIT REAL*8(a-h,o-z)
    !
    !   the methods in subroutine  step  approximate the solution near  x
    !   by a polynomial.  subroutine  intrp  approximates the solution at
    !   xout  by evaluating the polynomial there.  information defining this
    !   polynomial is passed from  step  so  intrp  cannot be used alone.
    !
    !   this code is completely explained and documented in the text,
    !   computer solution of ordinary differential equations:  the initial
    !   value problem  by l. f. shampine and m. k. gordon.
    !
    !   input to intrp --
    !
    !   all floating point variables are double precision
    !   the user provides storage in the calling program for the arrays in
    !   the call list
    DIMENSION y(neqn),yout(neqn),ypout(neqn),phi(neqn,16),psi(12)
    !   and defines
    !      xout -- point at which solution is desired.
    !   the remaining parameters are defined in  step  and passed to  intrp
    !   from that subroutine
    !
    !   output from  intrp --
    !
    !      yout(*) -- solution at  xout
    !      ypout(*) -- derivative of solution at  xout
    !   the remaining parameters are returned unaltered from their input
    !   values.  integration with  step  may be continued.
    !
    DIMENSION g(13),w(13),rho(13)
    DATA g(1)/1.0d0/,rho(1)/1.0d0/
    !
    hi = xout - x
    ki = kold + 1
    kip1 = ki + 1
    !
    !   initialize w(*) for computing g(*)
    !
    DO  i = 1,ki
       temp1 = i
       w(i) = 1.0d0/temp1
    ENDDO
    term = 0.0d0
    !
    !   compute g(*)
    !
    DO j = 2,ki
       jm1 = j - 1
       psijm1 = psi(jm1)
       gamma = (hi + term)/psijm1
       eta = hi/psijm1
       limit1 = kip1 - j
       DO  i = 1,limit1
          w(i) = gamma*w(i) - eta*w(i+1)
       ENDDO
       g(j) = w(1)
       rho(j) = gamma*rho(jm1)
       term = psijm1
    ENDDO
    !
    !   interpolate
    !
    DO  l = 1,neqn
       ypout(l) = 0.0d0
       yout(l) = 0.0d0
    ENDDO
    DO  j = 1,ki
       i = kip1 - j
       temp2 = g(i)
       temp3 = rho(i)
       DO l = 1,neqn
          yout(l) = yout(l) + temp2*phi(l,i)
          ypout(l) = ypout(l) + temp3*phi(l,i)
       ENDDO
    ENDDO
    DO l = 1,neqn
       yout(l) = y(l) + hi*yout(l)
    ENDDO

  END SUBROUTINE  intrp


  DOUBLE PRECISION FUNCTION d1mach (i)
    IMPLICIT NONE
    INTEGER :: i
    DOUBLE PRECISION :: b, x
    !***begin prologue  d1mach
    !***purpose  return floating point machine dependent constants.
    !***library   slatec
    !***category  r1
    !***type      single precision (d1mach-s, d1mach-d)
    !***keywords  machine constants
    !***author  fox, p. a., (bell labs)
    !           hall, a. d., (bell labs)
    !           schryer, n. l., (bell labs)
    !***description
    !
    !   d1mach can be used to obtain machine-dependent parameters for the
    !   local machine environment.  it is a function subprogram with one
    !   (input) argument, and can be referenced as follows:
    !
    !        a = d1mach(i)
    !
    !   where i=1,...,5.  the (output) value of a above is determined by
    !   the (input) value of i.  the results for various values of i are
    !   discussed below.
    !
    !   d1mach(1) = b**(emin-1), the smallest positive magnitude.
    !   d1mach(2) = b**emax*(1 - b**(-t)), the largest magnitude.
    !   d1mach(3) = b**(-t), the smallest relative spacing.
    !   d1mach(4) = b**(1-t), the largest relative spacing.
    !   d1mach(5) = log10(b)
    !
    !   assume single precision numbers are represented in the t-digit,
    !   base-b form
    !
    !              sign (b**e)*( (x(1)/b) + ... + (x(t)/b**t) )
    !
    !   where 0 .le. x(i) .lt. b for i=1,...,t, 0 .lt. x(1), and
    !   emin .le. e .le. emax.
    !
    !   the values of b, t, emin and emax are provided in i1mach as
    !   follows:
    !   i1mach(10) = b, the base.
    !   i1mach(11) = t, the number of base-b digits.
    !   i1mach(12) = emin, the smallest exponent e.
    !   i1mach(13) = emax, the largest exponent e.
    !
    !
    !***references  p. a. fox, a. d. hall and n. l. schryer, framework for
    !                 a portable library, acm transactions on mathematical
    !                 software 4, 2 (june 1978), pp. 177-188.
    !***routines called  xermsg
    !***revision history  (yymmdd)
    !   790101  date written
    !   960329  modified for fortran 90 (be after suggestions by ehg)      
    !***end prologue  d1mach
    !      
    x = 1.0d0
    b = RADIX(x)
    SELECT CASE (i)
    CASE (1)
       d1mach = b**(MINEXPONENT(x)-1) ! the smallest positive magnitude.
    CASE (2)
       d1mach = HUGE(x)               ! the largest magnitude.
    CASE (3)
       d1mach = b**(-DIGITS(x))       ! the smallest relative spacing.
    CASE (4)
       d1mach = b**(1-DIGITS(x))      ! the largest relative spacing.
    CASE (5)
       d1mach = LOG10(b)
    CASE default
       WRITE (*, fmt = 9000)
9000   FORMAT ('1error    1 in d1mach - i out of bounds')
       STOP
    END SELECT

  END FUNCTION  d1mach







