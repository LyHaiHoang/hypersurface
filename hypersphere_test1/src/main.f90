!> @brief
!! This programm show how to call the equipartition library

PROGRAM main

USE equipartition
    
IMPLICIT NONE

! Variables
INTEGER                                 :: nv, nv_fix      ! Number of vectors and of fixed vectors
INTEGER                                 :: D, nat          ! Dimension, number of atoms
REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: at_pos          ! Coordinates of each atom
REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: X               ! The different vectors
REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: sym_matrices    ! The different symmetry matrix
REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: perm_matrices   ! The different permutation matrix
REAL(DP), DIMENSION(:),     ALLOCATABLE :: line            ! The line read ineach file
INTEGER,  DIMENSION(:),     ALLOCATABLE :: atm_types       ! The different permutation matrix
REAL(DP), DIMENSION(:),     ALLOCATABLE :: nnangl          ! The final angles between nearest neighboors
REAL(DP)                                :: dmin            ! min of nnangl
REAL(DP)                                :: dmax            ! max of nnangl
REAL(DP)                                :: dmean           ! mean value of nnangl
REAL(DP)                                :: dstddev         ! standard deviation of nnangl
REAL(DP)                                :: Etot            ! Energy of the positions
INTEGER                                 :: Ndim            ! Total number of dimensions = nat*D 
INTEGER                                 :: nb_mat          ! number of symmetry matrix (or permutation matrix)
INTEGER                                 :: init_algo       ! Technique for the initialization
INTEGER                                 :: opt_algo        ! Technique for the optimization
INTEGER                                 :: i, j, k, l      ! counter
CHARACTER(len=250)                      :: filenameperm    ! file containing the permutation matrix
CHARACTER(len=250)                      :: filenamesym     ! file containing the symmetry matrix
CHARACTER(len=250)                      :: filenamelistvec ! file containing the fixed vectors       
CHARACTER(len=250)                      :: filenameatmpos  ! file containing the atomic positions 
CHARACTER(len=250)                      :: filenameoutput  ! file containing the vectors (output)
LOGICAL                                 :: flag            ! 
!
! ... Read the parameters
OPEN(1, FILE='input_variables.dat', STATUS='old')
READ(1, *) nat
READ(1, *) D
READ(1, *) nv
READ(1, *) init_algo
READ(1, *) opt_algo
READ(1, *) filenameoutput
READ(1, *) filenameatmpos
READ(1, *) filenamelistvec
READ(1, *) filenamesym
READ(1, *) filenameperm
CLOSE(1)
Ndim = nat*D
!
! ... Read the symmetry matrix
INQUIRE(file = TRIM(filenamesym), EXIST =  flag)
IF (flag) THEN
   OPEN(2, FILE=TRIM(filenamesym), STATUS='old')
   READ(2, *) nb_mat
   
   REWIND(2)
ELSE
   nb_mat=0
ENDIF
  
ALLOCATE(sym_matrices(D,D,nb_mat))
ALLOCATE(line(D**2))
DO k=1,nb_mat
    READ(2, *, IOSTAT=i) nb_mat
    IF (i/=0) EXIT
    READ(2, *) line
    DO l=1,D
        DO j=1,D
            sym_matrices(l,j,k) = line((l-1)*D+j)
        ENDDO
    ENDDO
ENDDO
DEALLOCATE(line)
CLOSE(2)
!
! ... Read the permutation matrix 
INQUIRE(file=TRIM(filenameperm), EXIST= flag)
IF (flag) THEN
    OPEN(3, FILE=TRIM(filenameperm), STATUS='old')
    READ(3, *) nb_mat
    
    REWIND(3)
ELSE
    nb_mat=0
ENDIF

ALLOCATE(perm_matrices(nat,nat,nb_mat))
ALLOCATE(line(nat**2))
DO k=1,nb_mat
    READ(3, *, IOSTAT=i) nb_mat 
    IF (i/=0) EXIT
    READ(3, *) line
    DO l=1,nat
        DO j=1,nat
            perm_matrices(l,j,k) = line((l-1)*nat+j)
        ENDDO
    ENDDO
ENDDO
DEALLOCATE(line)
CLOSE(3)
!
! ... Ddeterminate nv_fix
INQUIRE(file=TRIM(filenamelistvec), EXIST= flag)
IF (flag) THEN
    ALLOCATE(line(Ndim))
    OPEN(4, FILE=TRIM(filenamelistvec), STATUS='old')
    nv_fix = 0
    DO WHILE (.TRUE.)
        READ(4, *, IOSTAT=i) nat
        IF (i/=0) EXIT
        READ(4, *) line
        IF (SIZE(line)==0) EXIT
        nv_fix = nv_fix + 1
    ENDDO
    DEALLOCATE(line)
    REWIND(4)
ELSE
    nv_fix=0
ENDIF
!
! ... Read the vectors to fix
ALLOCATE(X(nat,D,(nb_mat+1)*(nv+nv_fix)))
ALLOCATE(nnangl((nb_mat+1)*(nv+nv_fix)))
ALLOCATE(line(Ndim))
X      = 99.0
nnangl = 99.0
DO k=1,nv_fix
    READ(4, *, IOSTAT=i) nat
    IF (i/=0) EXIT
    READ(4, *) line
     DO l=1,nat
         DO j=1,D
             X(l,j,k) = line((l-1)*D+j)
         ENDDO
     ENDDO
ENDDO
DEALLOCATE(line)
!
! ... To normalize the vectors
DO k=1, nv_fix
     X(:,:,k) = X(:,:,k)/NORM2(X(:,:,k))
 ENDDO
CLOSE(4)
!
! ... Read the atom positions
ALLOCATE(at_pos(nat,D))
ALLOCATE(atm_types(nat))
OPEN(5,FILE=TRIM(filenameatmpos), STATUS='old')
READ(5, *) 
READ(5, *) 
DO k = 1, nat
    READ(5, *) atm_types(k), (at_pos(k,j), j=1,D)
ENDDO
CLOSE(5)
!
! ... Call of equipartition_ND program
CALL equipartition_ND(nv,D,nat,nb_mat, nv_fix,sym_matrices,perm_matrices, init_algo, &
                      opt_algo, X, nnangl, dmean, dstddev, dmin, dmax, Etot )
!
! ... Write angles statistics                  
OPEN(10, FILE="angles.dat", STATUS="unknown", POSITION="rewind")
WRITE(10,*) 0, dmean 
WRITE(10,*) (nv+nv_fix)*(nb_mat+1), dmean
WRITE(10,*)  
WRITE(10,*) 0, dmean+dstddev
WRITE(10,*) (nv+nv_fix)*(nb_mat+1), dmean+dstddev
WRITE(10,*)  
WRITE(10,*) 0, dmean-dstddev
WRITE(10,*) (nv+nv_fix)*(nb_mat+1), dmean-dstddev
WRITE(10,*)  

DO i=1,(nv+nv_fix)*(nb_mat+1)
   WRITE(10,*) nnangl(i) 
ENDDO
CLOSE(10)

!
! ... Write vectors in Ovito's format
OPEN(10, FILE=TRIM(filenameoutput), STATUS="unknown", POSITION="rewind")
!
! ... OR All in one line
WRITE(10,'(I6)') (nv+nv_fix)*(nb_mat+1)
WRITE(10,*)
DO i=1,(nv+nv_fix)*(nb_mat+1)
   WRITE(10,*) 1, (X(j,:,i),j=1,nat)
ENDDO
! ... Or all separated
!  DO i=1,(nv+nv_fix)*(nb_mat+1)
!      WRITE(10,'(I3)') nat
!      WRITE(10,*) 
!      DO j=1,nat
!          WRITE(10,*) atm_types(j), at_pos(j,:) ,X(j,:,i)
!      ENDDO
!  ENDDO
CLOSE(10)

DEALLOCATE(at_pos)
DEALLOCATE(sym_matrices)
DEALLOCATE(perm_matrices)
DEALLOCATE(X)

END PROGRAM main
