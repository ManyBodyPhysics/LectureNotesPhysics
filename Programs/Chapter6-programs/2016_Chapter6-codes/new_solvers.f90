!************************************************************

SUBROUTINE zludcmp(A,det,n)
IMPLICIT NONE

INTEGER*4, INTENT(IN) :: n
COMPLEX*16, DIMENSION(0:n-1,0:n-1), INTENT(IN) :: A
COMPLEX*16, INTENT(OUT) :: det

! ******
! For input variables with zero-based indexing:
! Calculation of det(A) using LU decomposition. Note that in general,
! numerical overflow may result unless LOG(det(A)) is formed instead.

COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: A_aux
INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv
INTEGER :: m,i,info

ALLOCATE(A_aux(1:n,1:n))
ALLOCATE(ipiv(1:n))

A_aux(1:n,1:n) = A(0:n-1,0:n-1)

m = n
CALL ZGETRF(m,m,A_aux,m,ipiv,info)

IF (info /= 0) THEN
	
	WRITE (*,*) 'ZGETRF failed in zludcmp, code: ',info
	STOP

END IF

det = 1.d0

DO i = 1,n
                
	IF (ipiv(i) /= i) det = -det * A_aux(i,i)
	IF (ipiv(i) == i) det = det * A_aux(i,i)

END DO

DEALLOCATE(A_aux,ipiv)

END SUBROUTINE zludcmp

!************************************************************

SUBROUTINE zludcmp_ratio(A,detratio,ldetref,n)
IMPLICIT NONE

INTEGER*4, INTENT(IN) :: n
COMPLEX*16, DIMENSION(0:n-1,0:n-1), INTENT(IN) :: A
DOUBLE PRECISION, INTENT(IN) :: ldetref
COMPLEX*16 :: ldet
COMPLEX*16, INTENT(OUT) :: detratio

! ******
! For input variables with zero-based indexing:
! Calculation of det(A)/detref = det(A)/exp(ldetref) using LU decomposition,
! where detref is positive

COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: A_aux
INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv
INTEGER :: m,i,info

ALLOCATE(A_aux(1:n,1:n))
ALLOCATE(ipiv(1:n))

A_aux(1:n,1:n) = A(0:n-1,0:n-1)

m = n
CALL ZGETRF(m,m,A_aux,m,ipiv,info)

IF (info /= 0) THEN
	
	WRITE (*,*) 'ZGETRF failed in zludcmp, code: ',info
	STOP

END IF

ldet = 0.d0

DO i = 1,n

	IF (ipiv(i) /= i) ldet = ldet + LOG(-A_aux(i,i))
	IF (ipiv(i) == i) ldet = ldet + LOG(A_aux(i,i))

END DO

detratio = exp(ldet-ldetref)

!det = 1.d0
!
!DO i = 1,n
!                
!	IF (ipiv(i) /= i) det = -det * A_aux(i,i)
!	IF (ipiv(i) == i) det = det * A_aux(i,i)
!
!END DO

DEALLOCATE(A_aux,ipiv)

END SUBROUTINE zludcmp_ratio

!************************************************************

SUBROUTINE zludcmp_log(A,ldet,n)
IMPLICIT NONE

INTEGER*4, INTENT(IN) :: n
COMPLEX*16, DIMENSION(0:n-1,0:n-1), INTENT(IN) :: A
COMPLEX*16, INTENT(OUT) :: ldet

! ******
! For input variables with zero-based indexing:
! Calculation of log(det(A)) using LU decomposition.

COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: A_aux
INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv
INTEGER :: m,i,info

ALLOCATE(A_aux(1:n,1:n))
ALLOCATE(ipiv(1:n))

A_aux(1:n,1:n) = A(0:n-1,0:n-1)

m = n
CALL ZGETRF(m,m,A_aux,m,ipiv,info)

IF (info /= 0) THEN
	
	WRITE (*,*) 'ZGETRF failed in zludcmp, code: ',info
	STOP

END IF

ldet = 0.d0

DO i = 1,n

	IF (ipiv(i) /= i) ldet = ldet + LOG(-A_aux(i,i))
	IF (ipiv(i) == i) ldet = ldet + LOG(A_aux(i,i))

END DO

!det = 1.d0
!
!DO i = 1,n
!                
!	IF (ipiv(i) /= i) det = -det * A_aux(i,i)
!	IF (ipiv(i) == i) det = det * A_aux(i,i)
!
!END DO

DEALLOCATE(A_aux,ipiv)

END SUBROUTINE zludcmp_log

!***********************************************************

SUBROUTINE zdoinv(A,x,det,n)
IMPLICIT NONE

INTEGER*4, INTENT(IN) :: n
COMPLEX*16, DIMENSION(0:n-1,0:n-1), INTENT(IN) :: A
COMPLEX*16, DIMENSION(0:n-1,0:n-1), INTENT(OUT) :: x
COMPLEX*16, INTENT(OUT) :: det

! ******
! For input variables with zero-based indexing:
! Calculation of A^-1 from A*x = I including det(A).

COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: A_aux,B_aux
INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv
INTEGER :: m,nRhs,i,info

ALLOCATE(A_aux(1:n,1:n))
ALLOCATE(ipiv(1:n))

A_aux(1:n,1:n) = A(0:n-1,0:n-1)

m = n
CALL ZGETRF(m,m,A_aux,m,ipiv,info)

IF (info /= 0) THEN
	
	WRITE (*,*) 'ZGETRF failed in zdoinv, code: ',info
	STOP

END IF

det = 1.d0

DO i = 1,n
                
	IF (ipiv(i) /= i) det = -det * A_aux(i,i)
	IF (ipiv(i) == i) det = det * A_aux(i,i)

END DO

nRhs = n
ALLOCATE(B_aux(1:n,1:nRhs))
B_aux = 0.d0

DO i = 1,nRhs

	B_aux(i,i) = 1.d0                

END DO

CALL ZGETRS('N',m,nRhs,A_aux,m,ipiv,B_aux,m,info)

IF (info /= 0) THEN
	
	WRITE (*,*) 'ZGETRS failed in zdoinv, code: ',info
	STOP

END IF

x(0:n-1,0:n-1) = B_aux(1:n,1:n)

DEALLOCATE(A_aux,B_aux,ipiv)

END SUBROUTINE zdoinv

!***********************************************************

SUBROUTINE zdoinv_log(A,x,ldet,n)
IMPLICIT NONE

INTEGER*4, INTENT(IN) :: n
COMPLEX*16, DIMENSION(0:n-1,0:n-1), INTENT(IN) :: A
COMPLEX*16, DIMENSION(0:n-1,0:n-1), INTENT(OUT) :: x
COMPLEX*16, INTENT(OUT) :: ldet

! ******
! For input variables with zero-based indexing:
! Calculation of A^-1 from A*x = I including det(A).

COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: A_aux,B_aux
INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv
INTEGER :: m,nRhs,i,info

ALLOCATE(A_aux(1:n,1:n))
ALLOCATE(ipiv(1:n))

A_aux(1:n,1:n) = A(0:n-1,0:n-1)

m = n
CALL ZGETRF(m,m,A_aux,m,ipiv,info)

IF (info /= 0) THEN
	
	WRITE (*,*) 'ZGETRF failed in zdoinv, code: ',info
	STOP

END IF

ldet = 0.d0

DO i = 1,n

	IF (ipiv(i) /= i) ldet = ldet + LOG(-A_aux(i,i))
	IF (ipiv(i) == i) ldet = ldet + LOG(A_aux(i,i))

END DO

! det = 1.d0
!
! DO i = 1,n
!                
!	IF (ipiv(i) /= i) det = -det * A_aux(i,i)
!	IF (ipiv(i) == i) det = det * A_aux(i,i)
!
! END DO

nRhs = n
ALLOCATE(B_aux(1:n,1:nRhs))
B_aux = 0.d0

DO i = 1,nRhs

	B_aux(i,i) = 1.d0                

END DO

CALL ZGETRS('N',m,nRhs,A_aux,m,ipiv,B_aux,m,info)

IF (info /= 0) THEN
	
	WRITE (*,*) 'ZGETRS failed in zdoinv, code: ',info
	STOP

END IF

x(0:n-1,0:n-1) = B_aux(1:n,1:n)

DEALLOCATE(A_aux,B_aux,ipiv)

END SUBROUTINE zdoinv_log

!*************************************************************
