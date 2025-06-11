MODULE equipartition						! Create the module to calculate the distribution
  USE omp_lib							! Use this library to calculate parallel processing
  !
  PRIVATE
  !
  PUBLIC :: DP, PI, r, PHI, Riesz,                     &
            equipartition_ND, sort,                    &
            random_init_vec, golden_spiral,            &
            mat_curv,  val_graph,                      &
            force_elec, symetries, add_constrain,      & 
            min_force, min_dotproduct, add_constrain2D,&
            cart_to_hs, hs_to_cart
  INTEGER,  PARAMETER :: DP    = selected_real_kind(14,200)  !< @brief double precision
  REAL(DP), PARAMETER :: PI    = 3.14159265358979323846_DP   !< @brief pi number 
  REAL(DP), PARAMETER :: PHI   = 1.61803398874989484820_DP   !< @brief golden ratio=(1+sqrt(5))/2 
  REAL(DP), PARAMETER :: r     = 1.0_DP                      !< @briet hyperadius of the sphere
  REAL(DP), PARAMETER :: Riesz = 1.0_DP                      !< @briet Riesz parameter for the calculation of E/F

CONTAINS

!----------------------------------------------------------------------
SUBROUTINE equipartition_ND(nv, D, nat, nb_mat, nv_fix, sym_mat, perm_mat, init_algo, &
                            opt_algo, X, nnangl, dmean, dstddev, dmin, dmax, Etot )
!----------------------------------------------------------------------
!> @brief
!! This program equidistribute vectors on a hypersphere and returns
!! a sorted list of vectors and their average distance between first neighbors
  IMPLICIT NONE
  !
  ! ... Arguments
  INTEGER,                    INTENT(IN)    :: D, nat         ! Dimension, number of atoms
  INTEGER,                    INTENT(IN)    :: nv, nv_fix     ! Number of vectors and fixed vectors
  INTEGER,                    INTENT(IN)    :: nb_mat         ! Number of permutation matrices
  REAL(DP), DIMENSION(:,:,:), INTENT(IN)    :: sym_mat        ! The different symmetry matrix size=[D,D,nb_mat]
  REAL(DP), DIMENSION(:,:,:), INTENT(IN)    :: perm_mat       ! The different permutation matrix size=[nat,nat,nb_mat]
  INTEGER,                    INTENT(IN)    :: init_algo      ! The initialization method
  INTEGER,                    INTENT(IN)    :: opt_algo       ! The optization method
  REAL(DP), DIMENSION(:,:,:), INTENT(INOUT) :: X              ! The different vectors size=[nat,D, (nv+nv_fix)*(nb_mat+1)] 
  REAL(DP), DIMENSION(:),     INTENT(INOUT) :: nnangl         ! List of the angles with nearest neighboor size=(nv+nv_fix)*(nb_mat+1)
  REAL(DP),                   INTENT(INOUT) :: dmean          ! The average distance between first neighbors
  REAL(DP),                   INTENT(INOUT) :: dstddev        ! The standard deviation of distance between first neighbors
  REAL(DP),                   INTENT(INOUT) :: dmin           ! The minimum distance between first neighbors
  REAL(DP),                   INTENT(INOUT) :: dmax           ! The maximal distance between first neighbors
  REAL(DP),                   INTENT(INOUT) :: Etot           ! The total energy
  !
  ! ... Local Variables 
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE   :: mat            ! Distance matrix
  INTEGER                                   :: count_rate     ! Counters for knowing wall time
  INTEGER                                   :: count_start    ! Counters for knowing wall time
  INTEGER                                   :: count_end      ! Counters for knowing wall time
  !
  ! ... For mesuring wall time
  CALL SYSTEM_CLOCK(count_rate=count_rate)
  CALL SYSTEM_CLOCK(count_start)
  !
  ! ... Initialize distance matrix
  ALLOCATE (mat((nv+nv_fix)*(nb_mat+1),(nv+nv_fix)*(nb_mat+1))) 
  !
  ! ... Starting positions using optimal method or random
  IF (init_algo==1) THEN
     CALL golden_spiral(nat, D, nv, nv_fix, X)
  ELSE   
     CALL random_init_vec(0, D, nat, nv, X, nv_fix, init_algo)
  ENDIF  
  !
  ! ... Optimize the starting positions by minimizing forces
  IF (opt_algo==1 .OR. opt_algo==3) &
      CALL min_force(nv,D,nat,nb_mat,X,mat,nv_fix,sym_mat,perm_mat, Etot,Riesz)
  !
  ! ... Optimize the starting positions by minimizing dotproduct
  IF (opt_algo==2 .OR. opt_algo==3) &
      CALL min_dotproduct(nv,D,nat,nb_mat,X,nv_fix,sym_mat,perm_mat)
  !
  ! ... Calculate distances and statistiques for the output print 
  CALL mat_curv(nat, D, nv, nb_mat, X, mat, nv_fix)
  CALL val_graph(nv, mat, nv_fix, dmean, dstddev, dmin, dmax, nnangl)
  !
  ! ... Finalize
  DEALLOCATE (mat)
  !
  ! ... Printing wall time
  CALL SYSTEM_CLOCK(count_end)
  PRINT *, 'Wall time is: ', REAL(count_end-count_start,8)/REAL(count_rate,8), 'seconds.'
  ! 
END SUBROUTINE equipartition_ND

!----------------------------------------------------------------------
SUBROUTINE golden_spiral(nat, D, nv, nv_fix, X) 
!----------------------------------------------------------------------
!> @brief
!! Try to optimally sample the vectors on the hs using the golden spiral 
  IMPLICIT NONE
  !
  ! ... Arguments
  INTEGER,                    INTENT(IN)    :: D, nat                  ! Dimension, number of atoms
  INTEGER,                    INTENT(IN)    :: nv, nv_fix              ! Number of vectors and fixed vectors 
  REAL(DP), DIMENSION(:,:,:), INTENT(INOUT) :: X                       ! The different vectors
  ! ... Local Variables 
  INTEGER                                   :: i, k, j                 ! Counters 
  INTEGER                                   :: ndim                    ! Total dimension=nat*D
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE   :: hcx_                    ! One colone vectors in the cartesian hypercube
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE   :: x_, theta               ! One colone vectors in cartesian and hs coords
  REAL(DP), DIMENSION(:),     ALLOCATABLE   :: r3, r5                  ! One colone vectors in cartesian and hs coords
  REAL(DP)                                  :: shift                   ! Value of the shift for the 3D sphere
  REAL(DP)                                  :: inv,l,v,a,Delta,q,p     ! Values used for solving the inverse function
  REAL(DP)                                  :: r38, r36, r34, r32, r30 ! Coefs for the non-linear fit of F_3^-1
  REAL(DP)                                  :: r58, r56, r54, r52, r50 ! Coefs fot the non-linear fit of F_5^-1
  REAL(DP),                  DIMENSION(100) :: GR= &                   ! List of the first golden ratios
 [1.0, 1.618033988749894848204587, 1.324717957244746025960909, 1.220744084605759475361685, &
  1.167303978261418684256045, 1.134724138401519492605446, 1.112775684278705470629704, 1.096981557798559817908278, &
  1.085070245491450828336895, 1.075766066086837158059599, 1.068297188920841276369429, &
  1.062169167864255148458944, 1.057050575221228384881686, 1.052710920147558222713506, &
  1.048984934757034374979192, 1.045751024156341777656803, 1.042917732301786579722142, &
  1.040414947781847334028689, 1.038188019436449963433441, 1.036193717130683432452140, &
  1.034397396133807063097024, 1.032770966441042909329492, 1.031291412479247519263551, &
  1.029939696670523363272615, 1.028699935528251113143435, 1.027558772393021350509075, &
  1.026504894146267274183850, 1.025528654764581375071588, 1.024621779136584304596140, &
  1.023777127861464928136924, 1.022988508866290595641538, 1.022250525317837113978621, &
  1.021558451924347143715484, 1.020908133630808862993757, 1.020295902116472488880329, &
  1.019718506548559334613726, 1.019173055831054180109913, 1.018656970182196412191996, &
  1.018167940328672157219673, 1.017703892954414601856849, 1.017262961313372503366779, &
  1.016843460127677118547567, 1.016443864059417072092280, 1.016062789176217485797117, &
  1.015698976935896138235701, 1.015351280299596889581619, 1.015018651650519247869991, &
  1.014700132250153509917328, 1.014394843008478786591146, 1.014101976380970708487466, &
  1.013820789235123744334286, 1.013550596553796932718770, 1.013290765863049093579816, &
  1.013040712289039336241387, 1.012799894162667653651290, 1.012567809102428307734323, &
  1.012343990515855419328763, 1.012128004468286632919072, 1.011919446874725273464815, &
  1.011717940976562483615863, 1.011523135070006907496794, 1.011334700457406945242343, &
  1.011152329596359904470946, 1.010975734424682999403399, 1.010804644842055514450301, &
  1.010638807331498058952426, 1.010477983705890803429202, 1.010321949966495659158857, &
  1.010170495261977587669925, 1.010023420937751279574821, 1.009880539666639717116116, &
  1.009741674652844635726009, 1.009606658902115970842961, 1.009475334551785301766556, &
  1.009347552255011804755798, 1.009223170614190868388147, 1.009102055659006055687836, &
  1.008984080365073737140878, 1.008869124209544386411662, 1.008757072760392045587508, &
  1.008647817296449712394330, 1.008541254455538463814505, 1.008437285908296413641377, &
  1.008335818055543927734378, 1.008236761747227212475401, 1.008140032021166342336675, &
  1.008045547859998525292864, 1.007953231964855125229010, 1.007863010544443598506126, &
  1.007774813118324749296975, 1.007688572333283050600664, 1.007604223791784521126158, &
  1.007521705891603942722325, 1.007440959675782063113415, 1.007361928692144750242677, &
  1.007284558861680639830231, 1.007208798355132348466486, 1.007134597477209432893627, &
  1.007061908557879514544852, 1.006990685850237848089913]   
  !
  ! ... Initialization
  shift = 0.5                                                       ! shift of the pole for the 3D case
  r38   = 620.86318                                                 ! coef of the 8 order polynome for the non linear fit
  r36   = -297.50546                                                ! used to find F_3^-1 
  r34   = 47.76329
  r32   = -2.69733
  r30   = 0.399142
  r58   = 6151.5217                                                 ! coef of the 8 order polynome for the non linear fit
  r56   = -3212.62911                                               ! used to fin F_5^-1  
  r54   = 561.29499
  r52   = -35.58149
  r50   = 0.619622
  ndim  = nat*D                                                     ! Give the number of total dimension
  ALLOCATE(r3(nv))                                                  ! F_3^-1 values of the 8 order polynome fot non linear fit 
  ALLOCATE(r5(nv))                                                  ! F_5^-1 values of the 8 order polynome fot non linear fit 
  ALLOCATE(hcx_(ndim-1,nv))                                         ! Size hcx_ according to total of dimension and vectors
  ALLOCATE(x_(ndim,nv))                                             ! Size x_ according to total of dimension and vectors
  ALLOCATE(theta(ndim-1,nv))                                        ! Size theta according to total of dimension and vectors, in hs
  hcx_(:,:)=0.5                                                    
  x_(:,:)=99.0                                                       
  theta(:,:)=99.0                                                    
  !
  ! ... Equidistribution on the ndim hypercube => try to find better
  ! DO i=1, ndim-1
  !    DO k=1, nv 
  !       hcx_(i,k)   = MOD(k/GR(ndim)**i,1.0_DP)                   ! Use equidistribution theorem with the powers of n-th golden ratio
  !    ENDDO
  ! ENDDO
  DO i=1, ndim-2
     DO k=1, nv 
        hcx_(i,k)   = MOD(k/GR(ndim-1)**i,1.0_DP)                   ! Use equidistribution theorem with the powers of n-th golden ratio
     ENDDO
  ENDDO
  ! ... Angle of a spherical cap with area= 1/N the area of the hypersphere:
  !shift=ASIN(((ndim-1.0_DP)*GAMMA(0.5*(ndim-1.0_DP))*GAMMA(0.5)/(2.0*GAMMA(0.5*ndim)*nv))**(1/(ndim-1.0_DP)))
  shift=1.0_DP/2.0                                                  ! Correct approximated shift
  DO k=1, nv 
     hcx_(ndim-1,k) = (k*1.0_DP+shift)/(nv*1.0_DP+2*shift)        ! last coordinates is simply k/n  
  ENDDO
  !
  ! ... For debug, used to print the theta functions with hcx all linear from 0 to 1
  !DO k=1, nv 
  !   hcx_(:,k) = REAL(k-1)/REAL(nv-1)
  !ENDDO
  !
  ! ... Calculate angles in hs coordinates. Dim 1 to 6 have an exact value of the reciproque F
  DO i=1, ndim-1
     SELECT CASE (i)
        ! ... F_1(t)=t/ 2pi -> exact inverse exists
        CASE (1)
           theta(1,:) = 2*PI*hcx_(1,:)
        !   
        ! ... F_2(t)=0.5 - 0.5cos(t) -> exact inverse exists  
        CASE (2)
           theta(2,:) = ACOS(1-2*hcx_(2,:))
        !   
        ! ... F_3(t)=t/pi + sin(2t)/ 2pi -> Inverse expressed as a non linear combination of f_2^-1 and f_4^-1
        CASE (3)
           r3(:)= r38*(hcx_(3,:)-0.5)**8+r36*(hcx_(3,:)-0.5)**6+r34*(hcx_(3,:)-0.5)**4+r32*(hcx_(3,:)-0.5)**2+r30
           theta(3,:) = r3(:)*( ACOS(1-2*hcx_(3,:)) )+&
                    (1-r3(:))*( ACOS(2*COS( (ACOS(2*hcx_(3,:)-1))/3.0 +4*PI/3.0  ) ))
        !        
        ! ... F_4(t)=0.5 - 9cos(t)/16 +cos(3t)/16 -> exact inverse exists
        CASE (4)
           theta(4,:) = ACOS(2*COS( (ACOS(2*hcx_(4,:)-1))/3.0 +4*PI/3.0  ) )
        !   
        ! ... F_5(t)=0.5 - 9cos(t)/16 +cos(3t)/16 -> as for case 3 but of f_4^-1 and as f_6^-1 but for 5
        CASE (5)
           l=GAMMA((i*1.0_DP+1)/2)/(SQRT(PI)*GAMMA(i*1.0_DP/2))     
           r5(:)= r58*(hcx_(5,:)-0.5)**8+r56*(hcx_(5,:)-0.5)**6+r54*(hcx_(5,:)-0.5)**4+r52*(hcx_(5,:)-0.5)**2+r50
           theta(5,:) = r5(:)*( ACOS(2*COS( (ACOS(2*hcx_(5,:)-1))/3.0 +4*PI/3.0  ) ) )+ &
                    (1-r5(:))*( 0.5*PI+ ATAN(0.5/l* ATANH(2*hcx_(5,:)-1.0)) )
        !        
        ! ... F_6(t)=0.5 - 9cos(t)/16 +cos(3t)/16 -> as general but with simplified formulae
        CASE (6)
           l=GAMMA((i*1.0_DP+1)/2)/(SQRT(PI)*GAMMA(i*1.0_DP/2))     
           theta(6,:) = 0.5*PI+ ATAN(0.5/l* ATANH(2*hcx_(i,:)-1.0))
        !
        ! ... Approxiamtion of the reciprocal function
        CASE DEFAULT                                                
           l=GAMMA((i*1.0_DP+1)/2)/(SQRT(PI)*GAMMA(i*1.0_DP/2))     ! l=F_i'(pi/2)
           v=1-i*1.0_DP+8*l**2
           a = 2*SQRT( (192 + 8*PI**2*v) / (3*PI**6) ) *  &
              COS((ACOS(  ( 1152/(192*PI**2+8*PI**4*v)) * &
              SQRT(3*PI**6/(192+8*PI**2*v))) + 4*PI )/3) + 4/PI**2
           p          = 1.0/a -PI**2/4.0
           DO k=1, nv 
              inv = 0.5*PI+ ATAN(0.5*(1-0.25*a*PI**2)/l* ATANH(2*hcx_(i,k)-1.0))
              ! This is the inverse function of the third order polynome:
              ! P(x)=(1-a*pi^2/4)(x-pi/2)+a(x-pi/2)^3 determined with the Cardan Formulae
              q          = (0.5*PI-inv)/a
              Delta      = -4*p**3 -27.0*q**2
              theta(i,k) = 0.5*PI+ (0.5*(-q+SQRT(-Delta/27.0)))**(1.0/3.0) &
                                 - (0.5*( q+SQRT(-Delta/27.0)))**(1.0/3.0)
           ENDDO
        !   
     END SELECT 
  ENDDO
  !
  ! ... Avoid NaN in dim 3,4,5 due to ACOS
  IF (ndim>3) THEN
     WHERE (hcx_(3:min(5,ndim-1),:) ==1)
         theta(3:min(5,ndim-1),:)=PI
     ENDWHERE
  ENDIF
  !
  ! ... Switch the hyperspherical coords to cartesian in ndim
  CALL hs_to_cart(nv, ndim, theta, x_)
  !
  ! ... Create a D*nat table instead of ndim column vector
  DO k=nv_fix+1, nv_fix+nv  ! the nv_fix first vectors are already into X(:,:,1...nv_fix)
      DO i=1, nat
          DO j=1, D
              X(i,j,k)=x_((i-1)*D+j,k)
          ENDDO
      ENDDO
  ENDDO
  ! 
  ! ... PRINT for debug
  !DO k=1, nv                                                      ! Prints 1 line vectors
  !  PRINT*, "v=",k, "HC    = ", hcx_(:,k)
  !  PRINT*, "v=",k, "Theta = ", theta(:,k)
  !  PRINT*, "v=",k, "Final = ", x_(:,k)
  !ENDDO
  !DO i=1, ndim-1                                                      ! Prints one componant per line
  !   DO k=1, nv 
  !     PRINT *,theta(i,k)
  !    ! PRINT *,hcx_(i,k)
  !   ENDDO  
  !   PRINT *,""
  !ENDDO 
  !v=999.0                                                           
  !DO i=1, nv                                                      ! Prints distance between nearest neighbors 
  !   DO k=1, nv
  !      l = NORM2(hcx_(:,i)-hcx_(:,k))                              ! Calculation distance between vectors
  !      IF (l<v .AND. l>0) v=l
  !   ENDDO
  !   PRINT*,v!*((REAL(nv))**(1.0/REAL(ndim-1)))                   ! normalization by average dist in the hypercube
  !ENDDO
  !
  ! ... Final dealocation
  DEALLOCATE(x_)
  DEALLOCATE(theta)
  !
  ! ... Normalize each vector so that it stays on the hypersphere (for numerical only)
  !DO k=nv_fix+1,nv_fix+nv
  !    X(:,:,k) = X(:,:,k)/NORM2(X(:,:,k))
  !ENDDO
  !PRINT *, 'End of initialization'
  ! 
END SUBROUTINE golden_spiral 

!----------------------------------------------------------------------
SUBROUTINE random_init_vec(zseed, D, nat, nv, X, nv_fix, init_algo)
!----------------------------------------------------------------------
!> @brief
!! Create a tabular of normalized random vectors
  IMPLICIT NONE
  !
  ! ... Arguments
  INTEGER,                    INTENT(IN)    :: D, nat      ! Dimension, number of atoms
  INTEGER,                    INTENT(IN)    :: nv, nv_fix  ! Number of vectors and fixed vectors 
  REAL(DP), DIMENSION(:,:,:), INTENT(INOUT) :: X           ! The different vectors
  INTEGER,                    INTENT(IN)    :: zseed       ! The random seed
  INTEGER,                    INTENT(IN)    :: init_algo   ! The method used (2= uniform, 3=normal)
  ! ... Local Variables 
  INTEGER                                   :: i, k        ! Counters 
  REAL(DP)                                  :: zrand       ! Random number used create the first seed
  INTEGER                                   :: state_size  ! Number of seeds to initialize random numbers
  INTEGER                                   :: zseed_tmp   ! Initial seed used to fill the array of seeds
  INTEGER,                    ALLOCATABLE   :: state(:)    ! Array of seeds to initialize random numbers 
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE   :: u1,u2       ! Usefull for the Gaussian distribution
  !
  ! ... Initialize random seed array with the seed or not
  zseed_tmp = zseed
  IF( zseed_tmp == 0) THEN            ! Value is processor dependant and different for each run
      CALL RANDOM_SEED()
      CALL RANDOM_NUMBER(zrand)  
      zseed_tmp = INT(zrand *10e8)
  ENDIF
  CALL RANDOM_SEED(size=state_size)   ! In fortran, the seed is an array of size state_size 
  ALLOCATE(state(state_size))
  DO i=1, state_size
    state(i)=zseed_tmp**(i+5)         ! Put some entropy in the state
  ENDDO
  CALL RANDOM_SEED(put=state)         ! Initialize random state
  !
  ! ... Fill the vectors with a uniform distribution between -1 and 1 (BAD choice)
  IF (init_algo==2) THEN                   
     CALL RANDOM_NUMBER(X(:,:,nv_fix+1:nv_fix+nv)) 
     X(:,:,nv_fix+1:nv_fix+nv) = 2.0*X(:,:,nv_fix+1:nv_fix+nv) - 1
  ENDIF
  !
  ! ... Fill the vectors with a normal distribution between -1 and 1 (Box-Muller transformation) BETTER
  IF (init_algo==3) THEN                   
     ALLOCATE(u1(nat,D,nv))
     ALLOCATE(u2(nat,D,nv))
     CALL RANDOM_NUMBER(u1(:,:,:))
     CALL RANDOM_NUMBER(u2(:,:,:))
     X(:,:,nv_fix+1:nv_fix+nv) = SQRT(-2.0*LOG(u1(:,:,:))) * COS(2.0*3.14*u2(:,:,:))
     DEALLOCATE(u1)
     DEALLOCATE(u2)
  ENDIF
  !
  ! ... Normalize each vector so that it stays on the hypersphere
  DO k=nv_fix+1,nv+nv_fix
      X(:,:,k) = X(:,:,k)/NORM2(X(:,:,k))
  ENDDO
  !
END SUBROUTINE random_init_vec

!----------------------------------------------------------------------
SUBROUTINE min_force(nv,D,nat,nb_mat,X,mat,nv_fix,sym_matrices,perm_matrices, Etot, s)
!----------------------------------------------------------------------
!> @brief
!! Modifies the vectors positions on the hypersphere such that their interaction energy is minimum
  IMPLICIT NONE
  !
  ! ... Arguments
  INTEGER,                    INTENT(IN)    :: D, nat                    ! Dimension, number of atoms
  INTEGER,                    INTENT(IN)    :: nv, nv_fix                ! Number of vectors and fixed vectors
  INTEGER,                    INTENT(IN)    :: nb_mat                    ! Number of permutation matrix
  REAL(DP), DIMENSION(:,:,:), INTENT(INOUT) :: X                         ! The different vectors
  REAL(DP), DIMENSION(:,:) ,  INTENT(INOUT) :: mat                       ! Distance matrix
  REAL(DP), DIMENSION(:,:,:), INTENT(IN)    :: sym_matrices              ! The symmetry matrices
  REAL(DP), DIMENSION(:,:,:), INTENT(IN)    :: perm_matrices             ! The permutation matrices
  REAL(DP),                   INTENT(OUT)   :: Etot                      ! The final total energy
  REAL(DP),                   INTENT(IN)    :: s                         ! Riesz parameter for calculating E/F
  !
  ! ... Local Variables 
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE   :: Xold                      ! The different vectors
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE   :: F, F_old                  ! Force vector (nat,D,nv)
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE   :: F_perp, Fp_old, F_rad     ! Componants of the forces (nat,D,nv)
  REAL(DP)                                  :: eps                       ! Energy threshold to stop min
  REAL(DP)                                  :: alpha                     ! Force factor for displacement
  REAL(DP)                                  :: Dx_max                    ! Maximum allowed displacement
  REAL(DP)                                  :: Etot_old                  ! To save previous total energy
  INTEGER                                   :: p, try, nsteps            ! Counters
  INTEGER                                   :: nvtot                     ! Total number of vectors
  LOGICAL                                   :: do_exit                   ! Exit condition
  !
  ! ... Allocate tabulars
  nvtot = (nv+nv_fix)*(nb_mat+1)
  ALLOCATE (Xold(nat,D,nvtot))
  ALLOCATE (F(nat,D,nvtot))
  ALLOCATE (F_old(nat,D,nvtot))
  ALLOCATE (Fp_old(nat,D,nvtot))
  ALLOCATE (F_rad(nat,D,nvtot))
  ALLOCATE (F_perp(nat,D,nvtot))
  !
  ! ... Initialization
  eps       = 0.0000000000001_DP*(nv**2)
  alpha     = 10.0_DP*D*nat/(nvtot*1.0_DP)
  Dx_max    = ((2.0_DP*PI**(D*nat/2.0_DP))/(GAMMA(D*nat/2.0_DP)*nvtot*1.0_DP))**(1.0_DP/(D*nat*1.0_DP-1.0)) ! (Surface/nvec)^(1/(dim-1))
  try       = 0
  nsteps    = 0
  do_exit   = .FALSE.
  F_perp    = 0.0_DP
  Fp_old    = 0.0_DP
  F_old     = 0.0_DP
  F         = 0.0_DP
  Etot_old  = 99999999999999.0_DP
  Xold      = X
  !
  ! ... Run the minimization
  DO WHILE ( .NOT. do_exit )
      !
      ! ... Displace vectors according to forces, constrain, symetries and norm
      X = X + F_perp*alpha
      CALL add_constrain(X, Xold, nat, D, nv, nv_fix, nb_mat)
      CALL symetries(nv, nat, X, sym_matrices, perm_matrices, nv_fix)
      DO p=1,(nv+nv_fix)*(nb_mat+1)
          X(:,:,p) = X(:,:,p)/NORM2(X(:,:,p))
      ENDDO
      !
      ! ... Calculate Forces and their radial and tengeantial componants
      CALL mat_curv(nat, D, nv, nb_mat, X, mat, nv_fix)
      CALL force_elec(nv, nb_mat, X, F, F_perp, Etot, mat, nv_fix, D, nat, s)
      !
      ! ... For debug, print if needed
      PRINT*,nsteps,"E=",Etot, 'a=',alpha,'Fpmax=',MAXVAL(F_perp),'DE=',Etot_old-Etot,'DX=',MAXVAL(X-Xold),'try=',try 
      ! 
      ! ... Check if exit
      !do_exit = ((( (Etot_old-Etot)<eps .OR. MAXVAL(X-Xold)<(PI/4000.0)) &
      !            .AND. (try>3000) )                 &
      !            .OR. MAXVAL(F_perp)<0.00000000001) &
      do_exit = (nsteps>500) 
      nsteps  = nsteps+1
      !
      ! ... Modify alpha according to the new energy
      IF ( (Etot-Etot_old) .GT. 0.0_DP )  THEN
         ! ... Move is not accepted 
         X        = Xold
         F_perp   = Fp_old
         F        = F_old
         Etot     = Etot_old
         alpha    = alpha*0.7_DP
         try      = try+1
         !IF (nsteps <100) try =0 ! time to find the correct alpha
      ELSE
         ! ... Move is accepted 
         Etot_old = Etot
         Fp_old   = F_perp
         F_old    = F
         Xold     = X
         alpha    = alpha*1.1_DP
         IF ( MAXVAL(F_perp)*alpha > Dx_max ) alpha = Dx_max/MAXVAL(F_perp)
      ENDIF
      !
  ENDDO
  !
  ! ... Finalize
  DEALLOCATE(Xold)
  DEALLOCATE(F)
  DEALLOCATE(Fp_old)
  !
END SUBROUTINE min_force

!----------------------------------------------------------------------
SUBROUTINE add_constrain(X, Xold, nat, D, nv, nv_fix, nb_mat)
!----------------------------------------------------------------------
!> @brief
!! Add a constrain on the vectors positions during the minimization
  IMPLICIT NONE
  !
  ! ... Arguments
  REAL(DP), DIMENSION(:,:,:), INTENT(INOUT) :: X, Xold    ! The different vectors and their previous position
  INTEGER,                    INTENT(IN)    :: nat, D     ! Number of atoms, Dimension 
  INTEGER,                    INTENT(IN)    :: nv, nv_fix ! Number of vectors and fixed vectors
  INTEGER,                    INTENT(IN)    :: nb_mat     ! Number of permutation matrices
  ! ... Local Variables 
  INTEGER                                   :: i          ! Counters
  ! 
  ! ... Reset the position of the fixed atoms
  DO i=1, nv_fix
    X(:,:,i) = Xold(:,:,i)
  ENDDO
  i=nat+D+nv+nv_fix+nb_mat
  ! ... Here the user can write explicitelly what he wants
  !i=nb_mat+nv+nv_fix+nat+D ! dummy staff
  ! X(:,2,:)=0
  !
END SUBROUTINE add_constrain 

!----------------------------------------------------------------------
SUBROUTINE force_elec(nv, nb_mat, X, F, F_perp, Etot, mat, nv_fix, D, nat, s)
!----------------------------------------------------------------------
!> @brief
!! Calculate the repulsive force between vectors using an electromagnetic
!! interaction. The Force field is in 1/r^(s+1)
  !
  ! ... Arguments
  INTEGER,                    INTENT(IN)  :: nv, nv_fix ! Number of vectors and of fixed vectors
  INTEGER,                    INTENT(IN)  :: nb_mat     ! Number of permutation matrix
  REAL(DP), DIMENSION(:,:,:), INTENT(IN)  :: X          ! The different vectors
  REAL(DP), DIMENSION(:,:,:), INTENT(OUT) :: F          ! Force vector (nat,D,nv)
  REAL(DP), DIMENSION(:,:,:), INTENT(OUT) :: F_perp     ! Force tangeantielle (nat,D,nv)
  REAL(DP),                   INTENT(OUT) :: Etot       ! Total energy 
  REAL(DP), DIMENSION(:,:) ,  INTENT(IN)  :: mat        ! Distance matrix  
  INTEGER,                    INTENT(IN)  :: D, nat     ! Number dimensions and atoms
  REAL(DP),                   INTENT(IN)  :: s          ! Riesz parameter for E/F 
  ! ... Local Variables 
  INTEGER                                 :: ivec, jvec ! Counters over the vectors
  INTEGER                                 :: iat, iD    ! Counters over atoms and dimensions
  REAL(DP)                                :: dotp       ! Dot product
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: F_rad      ! Radial force vector (nat,D,nv)
  !
  ! ... Initialize
  ALLOCATE (F_rad(nat,D,(nv+nv_fix)*(nb_mat+1)))
  !
  ! ... Calculate Total force and total energy
  F = 0
  Etot = 0
  !$OMP PARALLEL DO PRIVATE(ivec) REDUCTION(+:Etot) 
  DO ivec=1, (nv+nv_fix)*(nb_mat+1)
     DO jvec=1, (nv+nv_fix)*(nb_mat+1)
        IF (ivec/=jvec) THEN
           F(:,:,ivec) = F(:,:,ivec)+2*s*(X(:,:,ivec)-X(:,:,jvec))/(mat(ivec,jvec)**(s+2))
           Etot = Etot + 1.0_DP/(mat(ivec,jvec)**s)   
        ENDIF   
     ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  !
  ! ... Calculate tangeantial and radial componants
  !$OMP PARALLEL DO PRIVATE(ivec, iD, iat, dotp)
  DO ivec=1,(nv+nv_fix)*(nb_mat+1)
     dotp=0.0
     DO iD=1, D
        DO iat=1, nat
           dotp=dotp+F(iat,iD,ivec)*X(iat,iD,ivec)
        ENDDO
     ENDDO
     F_rad(:,:,ivec)=dotp*X(:,:,ivec)
  ENDDO
  !$OMP END PARALLEL DO
  F_perp=F-F_rad
  DEALLOCATE (F_rad)
  !
END SUBROUTINE force_elec

!----------------------------------------------------------------------
SUBROUTINE mat_curv(nat, D, nv,nb_mat,X,mat,nv_fix)
!----------------------------------------------------------------------
!> @brief
!! Create the matrice of curvilign distance between atoms
!! For this the vectors must be in single column
  IMPLICIT NONE
  !
  ! ... Arguments
  INTEGER,                    INTENT(IN)  :: nv, nv_fix ! Number of vectors, number of fixed vectors
  INTEGER,                    INTENT(IN)  :: nat, D     ! Number of vectors, number of fixed vectors
  INTEGER,                    INTENT(IN)  :: nb_mat     ! Number of symetry matrices
  REAL(DP), DIMENSION(:,:,:), INTENT(IN)  :: X          ! the different vectors
  REAL(DP), DIMENSION(:,:) ,  INTENT(OUT) :: mat        ! distance matrix
  ! ... Local Variables 
  INTEGER                                 :: ivec, jvec ! Counters over vectors
  INTEGER                                 :: iat, id    ! Counters over atoms and dim
  REAL(DP)                                :: prod       ! Intermediate product
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: x_         ! Ndim cartesian coordinates of the vectors
  !
  ALLOCATE(x_(nat*D,(nv+nv_fix)*(nb_mat+1)))
  !
  ! ... Create a column vector x_ instead of a D*nat table
  !$OMP PARALLEL DO 
  DO ivec=1, (nv+nv_fix)*(nb_mat+1)
     DO iat=1, nat
        DO id=1, D
            x_((iat-1)*D+id,ivec) = X(iat,id,ivec)
        ENDDO
     ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  !
  ! ... Fill the distance matrix
  !$OMP PARALLEL DO PRIVATE(jvec, ivec, prod) SHARED(x_, mat)
  DO jvec=1,(nv+nv_fix)*(nb_mat+1)
     DO ivec=jvec+1,(nv+nv_fix)*(nb_mat+1)
        prod=DOT_PRODUCT(x_(:,ivec),x_(:,jvec))
        !prod=NORM2(x_(:,ivec)-x_(:,jvec))       ! Euclidian distance (not relevant on an HS) 
        IF (prod>1) prod=1
        IF (prod<-1) prod=-1
        mat(ivec,jvec)=ACOS(prod)
        mat(jvec,ivec) = mat(ivec,jvec)          ! the matrix is symetric 
     ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  !
  DEALLOCATE(x_)
  !
END SUBROUTINE mat_curv

!----------------------------------------------------------------------
SUBROUTINE cart_to_hs(nvec, ndim, x_, theta)
!----------------------------------------------------------------------
!> @brief
!! Convert a vector from cartesian to hyperspherical coordinates
!! Here the input vector is a column vector instead of a table
  IMPLICIT NONE
  !
  ! ... Arguments
  INTEGER,                   INTENT(IN)    :: nvec       ! Number of vectors
  INTEGER,                   INTENT(IN)    :: ndim       ! Number of Dimension
  REAL(DP), DIMENSION(:,:),  INTENT(IN)    :: x_         ! the different vectors
  REAL(DP), DIMENSION(:,:),  INTENT(INOUT) :: theta      ! the different vectors
  ! ... Local Variables 
  INTEGER                                  :: ivec, id   ! Counters
  REAL(DP)                                 :: prod       ! Intermediate product
  !
  ! ... x_ is in cartesian coordinates, switch it to hyperspherical coordinates  
  DO ivec=1, nvec                                        
      IF (ndim>2) THEN
          ! ... Calculation of azimuthal angles (between 0 and PI)
          theta(1,ivec) = ACOS(x_(1,ivec)/r)                     ! Calculation of 1st theta
          prod = 1
          IF (ndim>3) THEN                                       ! After 3 dimensions, a recursive definition required
              DO id=2, ndim-2                                   
                  prod = prod*SIN(theta(id-1,ivec))       
                  theta(id,ivec) = ACOS(x_(id,ivec)/(r*prod))    ! Calculation of next dependent theta (recursive)
              ENDDO
          ENDIF
      ENDIF
      ! ... Calculation of polar angle (between 0 and 2PI)
      theta(ndim-1,ivec) = ATAN2(x_(ndim,ivec),x_(ndim-1,ivec))
      IF (theta(ndim-1,ivec) < 0) THEN
          theta(ndim-1,ivec) = 2*PI + theta(ndim-1,ivec)         ! So the value is between 0 and 2*PI
      ENDIF
  ENDDO
  !
END SUBROUTINE cart_to_hs

!----------------------------------------------------------------------
SUBROUTINE hs_to_cart(nvec, ndim, theta, x_)
!----------------------------------------------------------------------
!> @brief
!! Convert the vector theta from hyperspherical coordinates to x_ that is in cartesian coordinates
!! Here the input vector is a column vector instead of a table
  IMPLICIT NONE
  !
  ! ... Arguments
  INTEGER,                   INTENT(IN)    :: nvec       ! Number of vectors
  INTEGER,                   INTENT(IN)    :: ndim       ! Number of Dimension
  REAL(DP), DIMENSION(:,:),  INTENT(IN)    :: theta      ! the different vectors in hs coord
  REAL(DP), DIMENSION(:,:),  INTENT(INOUT) :: x_         ! the different vectors in cart coord
  ! ... Local Variables 
  INTEGER                                  :: ivec, id   ! Counters
  REAL(DP)                                 :: prod       ! Intermediate sinus product
  !
  ! ... theta is in hyperspherical coordinates, switch to cartesian 
  DO ivec=1, nvec
     x_(1,ivec) = r*COS(theta(ndim-1,ivec))              ! Calculation of 1st x_ 
     prod = 1
     DO id=2, ndim-1
         prod = prod*SIN(theta(ndim-id+1,ivec))
         x_(id,ivec) = r*prod*COS(theta(ndim-id,ivec))   ! Calculation of intermediate x_ 
     ENDDO
     x_(ndim,ivec) = r*prod*SIN(theta(1,ivec))           ! Last x_ is a product of all sin(old_theta) 
  ENDDO
  !
END SUBROUTINE hs_to_cart

!----------------------------------------------------------------------
SUBROUTINE symetries(nv,nat,X,sym_matrices,perm_matrices,nv_fix)
!----------------------------------------------------------------------
!> @brief
!! Apply the symetry and permutation matrices on the list of vectors
  IMPLICIT NONE
  !
  ! ... Arguments
  INTEGER,                    INTENT(IN)        :: nv, nat, nv_fix        ! number of vectors, number of atoms, number of fixed vectors
  REAL(DP), DIMENSION(:,:,:), INTENT(IN)        :: sym_matrices           ! list of the rotation matrix
  REAL(DP), DIMENSION(:,:,:), INTENT(IN)        :: perm_matrices          ! list of the permutation matrix
  REAL(DP), DIMENSION(:,:,:), INTENT(INOUT)     :: X                      ! the different vectors
  ! ... Local Variables 
  INTEGER                                       :: ivec, jmat, kat        ! counter above vectors, sym matrices and atoms
  !
  !$OMP PARALLEL DO PRIVATE(ivec, jmat, kat) SHARED(X, sym_matrices, perm_matrices)
  DO ivec=1,nv+nv_fix                                                     ! For all vectors,...
      DO jmat=1, SIZE(sym_matrices, 3)                                    ! from all dimension ...
          DO kat=1,nat                                                    ! & for all atoms
              ! ... Apply the rotation matrices
              X(kat,:,jmat*(nv+nv_fix)+ivec) = MATMUL(sym_matrices(:,:,jmat), X(kat,:,ivec))
          ENDDO
          ! ... Apply the permutation matrices
          X(:,:,jmat*(nv+nv_fix)+ivec) = MATMUL(perm_matrices(:,:,jmat), X(:,:,jmat*(nv+nv_fix)+ivec))
      ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  !
END SUBROUTINE symetries

!----------------------------------------------------------------------
SUBROUTINE sort(nv,D,nat,mat,X,nv_fix, nb_mat)
!----------------------------------------------------------------------
!> @brief
!! Sort the Vectors such that the last vector is as far as possible from the previous ones
  !
  ! ... Arguments
  INTEGER,                    INTENT(IN)    :: D, nat           ! Dimension, number of atoms
  INTEGER,                    INTENT(IN)    :: nv, nv_fix       ! Number of vectors, number of fixed vectors
  INTEGER,                    INTENT(IN)    :: nb_mat           ! Number of symetric matrix
  REAL(DP), DIMENSION(:,:) ,  INTENT(INOUT) :: mat              ! Distance matrix 
  REAL(DP), DIMENSION(:,:,:), INTENT(INOUT) :: X                ! The different vectors
  ! ... Local Variables 
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE   :: C                ! Temporary vector used for the exchange
  REAL(DP), DIMENSION(:),     ALLOCATABLE   :: dist             ! Temporary distance used for the exchange   
  REAL(DP), DIMENSION(:),     ALLOCATABLE   :: E                ! Energy 
  INTEGER                                   :: ivec, jvec, kvec ! Counters above vectors
  INTEGER                                   :: lmat             ! Counters above symetric matrix
  INTEGER                                   :: posmax           ! Position of  
  !
  ! ... Initialization
  ALLOCATE(E(nv+nv_fix))
  ALLOCATE(C(nat,D))
  ALLOCATE(dist((nv+nv_fix)*(nb_mat+1)))
  !
  ! ... Sort the symetric vectors with respect to the 1rst vector
  DO ivec=2, (nv_fix+nv)
     !
     ! ... Find the symetric vector that is the closest to the 1rst vector
     posmax=0
     DO lmat=1, nb_mat
        IF (mat(1,ivec+(nv+nv_fix)*posmax) > mat(1,ivec+(nv+nv_fix)*lmat)) posmax=lmat
     ENDDO
     !
     ! ... Exchange its position and the line and column of the distances matrice
     C(:,:)                         = X(:,:,ivec)
     X(:,:,ivec)                    = X(:,:,ivec+(nv+nv_fix)*posmax)
     X(:,:,ivec+(nv+nv_fix)*posmax) = C(:,:)
     dist(:)                        = mat(:,ivec)
     mat(:,ivec)                    = mat(:,ivec+(nv+nv_fix)*posmax)
     mat(:,ivec+(nv+nv_fix)*posmax) = dist(:)
     dist(:)                        = mat(ivec,:)
     mat(ivec,:)                    = mat(ivec+(nv+nv_fix)*posmax,:)
     mat(ivec+(nv+nv_fix)*posmax,:) = dist(:)
  ENDDO
  !
  ! ... Sort all the vectors such that they are as far as possible from the previous ones
  DO ivec=2, (nv_fix+nv-1)
     !
     ! ... For each vector, calculate the sum of the distances from the others
     E(:) = 0
     DO jvec=ivec, (nv+nv_fix)
         DO kvec=1, ivec-1
            DO lmat=0, nb_mat
               E(jvec) = E(jvec) + 1.0/mat(jvec,kvec+(nv+nv_fix)*lmat)
            ENDDO   
         ENDDO
     ENDDO
     !
     ! ... Find the vector that is the farthest from all the previous ones 
     posmax = ivec
     DO jvec=ivec+1, nv+nv_fix
         IF (E(jvec) < E(posmax))  posmax = jvec
     ENDDO
     !
     ! ... Exchange this vector and its symetric ones with the one at the ieme place.
     IF (.NOT.(posmax==ivec) ) THEN
        DO lmat=0, nb_mat
           C(:,:)                         = X(:,:,ivec+(nv+nv_fix)*lmat)
           X(:,:,ivec+(nv+nv_fix)*lmat)   = X(:,:,posmax+(nv+nv_fix)*lmat)
           X(:,:,posmax+(nv+nv_fix)*lmat) = C(:,:)
           dist(:)                        = mat(:,ivec+(nv+nv_fix)*lmat)
           mat(:,ivec+(nv+nv_fix)*lmat)   = mat(:,posmax+(nv+nv_fix)*lmat)
           mat(:,posmax+(nv+nv_fix)*lmat) = dist(:)
           dist(:)                        = mat(ivec+(nv+nv_fix)*lmat,:)
           mat(ivec+(nv+nv_fix)*lmat,:)   = mat(posmax+(nv+nv_fix)*lmat,:)
           mat(posmax+(nv+nv_fix)*lmat,:) = dist(:)
        ENDDO
     ENDIF
     !
  ENDDO
  !
  ! ...Finalization
  DEALLOCATE(C)
  DEALLOCATE(E)
  DEALLOCATE(dist)
  !
END SUBROUTINE sort

!----------------------------------------------------------------------
SUBROUTINE val_graph(nv, mat, nv_fix, dmean, dstddev, dmin, dmax, nnangl)
!----------------------------------------------------------------------
!> @brief
!! This routine returns the average value of distance between closest vectors
  IMPLICIT NONE
  !
  ! ... Arguments
  INTEGER,                    INTENT(IN)  :: nv, nv_fix   ! Number of vectors, number of fixed vectors
  REAL(DP), DIMENSION(:,:) ,  INTENT(IN)  :: mat          ! Distance matrix 
  REAL(DP),                   INTENT(OUT) :: dmean        ! Average distance between first neighbors
  REAL(DP),                   INTENT(OUT) :: dmax         ! maximum distance between first neighbors
  REAL(DP),                   INTENT(OUT) :: dmin         ! minimum distance between first neighbors
  REAL(DP),                   INTENT(OUT) :: dstddev      ! standard deviation from mean
  REAL(DP), DIMENSION(:),     INTENT(OUT) :: nnangl       ! list of distances with nearest neighboor
  ! ... Local Variables 
  INTEGER                                 :: ivec, jvec   ! Counter                   
  REAL(DP)                                :: tmp_min      ! Minimum distance between two vectors  
  !
  ! ... Initialization
  dmean   = 0.0_DP
  dmax    = 0.0_DP
  dmin    = 9.0_DP
  dstddev = 0.0_DP
  !
  ! ... Parse the distance matrix
  DO jvec=1, nv+nv_fix
     tmp_min = 9999.0_DP 
     DO ivec=1, nv+nv_fix
        ! ... Extract the minimum distance between vectors
        IF (mat(ivec,jvec) < tmp_min .AND. ivec/=jvec)  tmp_min = mat(ivec,jvec)
     ENDDO
     nnangl(jvec)=tmp_min
     dmean = dmean + tmp_min
     IF (tmp_min>dmax) dmax=tmp_min
     IF (tmp_min<dmin) dmin=tmp_min
  ENDDO
  !
  ! ... Take the average
  dmean = dmean/REAL(nv+nv_fix)
  !
  ! ... Calculate the standard deviation from mean value
  DO ivec=1, nv+nv_fix
      dstddev=dstddev+(nnangl(ivec)-dmean)**2
  ENDDO
  dstddev=(dstddev/REAL(nv+nv_fix))**0.5
  !
  ! ... All in degree
  dmean   = dmean   *180.0_DP/PI
  dmin    = dmin    *180.0_DP/PI
  dmax    = dmax    *180.0_DP/PI
  dstddev = dstddev *180.0_DP/PI
  nnangl  = nnangl  *180.0_DP/PI
  !
END SUBROUTINE val_graph

!----------------------------------------------------------------------
SUBROUTINE min_dotproduct(nv,D,nat,nb_mat,X, nv_fix,sym_mat,perm_mat)
!----------------------------------------------------------------------
!> @brief
!! Modifies the vectors on the hypersphere such as to minimize the dot product between them
  IMPLICIT NONE
  !
  ! ... Arguments
  INTEGER,                    INTENT(IN)    :: D, nat         ! Dimension, number of atoms
  INTEGER,                    INTENT(IN)    :: nv, nv_fix     ! Number of vectors and fixed vectors
  INTEGER,                    INTENT(IN)    :: nb_mat         ! Number of permutation matrix
  REAL(DP), DIMENSION(:,:,:), INTENT(INOUT) :: X              ! The different vectors
  REAL(DP), DIMENSION(:,:,:), INTENT(IN)    :: sym_mat        ! The symmetry matrices
  REAL(DP), DIMENSION(:,:,:), INTENT(IN)    :: perm_mat       ! The permutation matrices
  ! ... Local Variables 
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE   :: x_, x_old       ! The different in single column
  REAL(DP), DIMENSION(:),     ALLOCATABLE   :: x_tmp          ! Temporary vector
  REAL(DP)                                  :: PS, PS_max     ! Scalar product and its max value
  REAL(DP)                                  :: PS_max_new     ! New max scalar product after displacement
  REAL(DP)                                  :: exit_thr       ! Angle threshold for exit
  REAL(DP)                                  :: alpha          ! Force factor for displacement
  INTEGER                                   :: ndim           ! total dimension
  INTEGER                                   :: try, nsteps    ! Counters
  INTEGER                                   :: nvec           ! Total number of vectors
  INTEGER                                   :: i,j,k          ! Counters
  INTEGER                                   :: i_max, imaxnew ! Counters
  INTEGER                                   :: j_max, jmaxnew ! Counters
  LOGICAL                                   :: do_exit        ! Exit condition
  !
  ! ... Allocate tabulars
  ndim = nat*D
  nvec = (nv+nv_fix)*(nb_mat+1)
  ALLOCATE( x_(ndim,nvec)     )
  ALLOCATE( x_old(ndim,nvec) )
  ALLOCATE( x_tmp(ndim)      )
  !
  ! ... Create a column vector x_ instead of a D*nat table
  DO k=1, nvec
      DO i=1, nat
          DO j=1, D
              x_((i-1)*D+j,k) = X(i,j,k)
          ENDDO
      ENDDO
  ENDDO
  !
  ! ... Initialization
  x_old(:,:) = x_(:,:)
  alpha      = 0.01_DP/(nv**2)
  PS_max     = -2.0
  exit_thr   = -3.0
  try        = 0
  nsteps     = 0
  do_exit    = .FALSE.
  !
  ! ... Run the minimization
  DO WHILE ( .NOT. do_exit )
     !
     ! ... Extract the 2 vectors with maximal dot product
     PS_max = -2.0
     i_max   = -1
     j_max   = -1
     DO i=1, nvec
        DO j=i+1, nvec
           PS= DOT_PRODUCT( x_(:,i), x_(:,j) )
           !PRINT*, 'test all',i,j, PS
           IF (PS> PS_max) THEN
               PS_max = PS
               i_max  = i
               j_max  = j
           ENDIF    
        ENDDO
     ENDDO    
 !    PRINT*, 'Results', i_max, j_max, PS_max
     !PRINT*, i_max, x_(:, i_max)
     !PRINT*, j_max, x_(:, j_max)
     !
     ! ... Check if exit
     do_exit = ( (PS_max < exit_thr) .OR. (try>2000) .OR. (nsteps>6000) )    
     nsteps  = nsteps+1
     !
     ! ... Modify these vectors and normalize
     x_tmp(:)   = x_(:,i_max)
     x_(:,i_max) = x_(:,i_max) -alpha*x_(:,j_max)
     x_(:,j_max) = x_(:,j_max) -alpha*x_tmp(:)
     x_(:,i_max) = x_(:,i_max)/NORM2(x_(:,i_max))
     x_(:,j_max) = x_(:,j_max)/NORM2(x_(:,j_max))
     !PRINT*, 'After',  DOT_PRODUCT( x_(:,i_max), x_(:,j_max) )
     !PRINT*, i_max, x_(:, i_max)
     !PRINT*, j_max, x_(:, j_max)
     !
     CALL add_constrain2D(x_, x_old, ndim, nvec, nv_fix)
     PRINT*, 8
     PRINT*, ""
     DO k=1, nvec
        PRINT*,1, (x_(j,k), j=1,D)
     ENDDO

     ! ... Extract the new maximal dot product
     PS_max_new = -2.0
     DO k=1, nvec
        IF ( (k .NE. i_max) .AND. (k .NE. j_max) ) THEN
           PS = DOT_PRODUCT( x_(:,i_max), x_(:,k) )
           IF (PS> PS_max_new) THEN
              PS_max_new = PS
              imaxnew    = k
              jmaxnew    = i_max
           ENDIF   
           PS= DOT_PRODUCT( x_(:,i_max), x_(:,k) )
           IF (PS> PS_max_new) THEN
              PS_max_new = PS
              imaxnew    = k
              jmaxnew    = j_max
           ENDIF   
        ENDIF    
     ENDDO
     !PRINT*, 'PS_max_new=', PS_max_new, 'between', imaxnew, jmaxnew
     !
     ! ... Modify alpha according to the new energy
     IF ( DOT_PRODUCT( x_(:,i_max), x_(:,j_max) ) > PS_max )  THEN
     !IF ( PS_max_new > PS_max )  THEN
        ! ... Move is not accepted, reduce alpha 
        x_       = x_old
        alpha    = alpha*0.7_DP
        try      = try+1
        IF (nsteps <10) try =0 ! time to find the correct alpha
     ELSE
        ! ... Move is accepted, increase alpha 
        x_old    = x_
        alpha    = alpha*1.1_DP
        IF ( alpha> 0.05 ) alpha = 0.05
     ENDIF
     PRINT*, 8
     PRINT *, 'alpha=', alpha, 'try=', try,  (PS_max_new > PS_max )
     DO k=1, nvec
        PRINT*,1, (x_(j,k), j=1,D)
     ENDDO
     !
  ENDDO
  !
  ! ... Recreate a D*nat table instead of ndim column vector
  DO k=1, nvec
      DO i=1, nat
          DO j=1, D
              X(i,j,k)=x_((i-1)*D+j,k)
          ENDDO
      ENDDO
  ENDDO
  !
  DEALLOCATE(x_)
  DEALLOCATE(x_old)
  DEALLOCATE(x_tmp)
  !
END SUBROUTINE min_dotproduct

!----------------------------------------------------------------------
SUBROUTINE add_constrain2D(X, Xold, D, nvtot, nv_fix)
!----------------------------------------------------------------------
!> @brief
!! Add a constrain on the vectors positions during the minimization
  IMPLICIT NONE
  !
  ! ... Arguments
  REAL(DP), DIMENSION(:,:), INTENT(INOUT)   :: X, Xold    ! The different vectors and their previous position
  INTEGER,                    INTENT(IN)    :: D          ! Dimension 
  INTEGER,                    INTENT(IN)    :: nvtot, nv_fix ! Number of vectors and fixed vectors
  ! ... Local Variables 
  INTEGER                                   :: i          ! Counters
  ! 
  ! ... Reset the position of the fixed atoms
  DO i=1, nv_fix
    X(:,i) = Xold(:,i)
  ENDDO
  i=nvtot+D
  !
  ! ... Here the user can write explicitelly what he wants
  !i=nb_mat+nv+nv_fix+nat+D ! dummy staff
  ! X(:,2,:)=0
  !
END SUBROUTINE add_constrain2D

END MODULE equipartition
