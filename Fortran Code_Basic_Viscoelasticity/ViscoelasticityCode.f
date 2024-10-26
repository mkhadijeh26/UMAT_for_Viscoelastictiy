!!! This one is working 100%

		    SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
      1 RPL,DDSDDT,DRPLDE,DRPLDT,
      2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
      3 NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
      4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
 !.. 
               INCLUDE 'ABA_PARAM.INC'
 !..
       CHARACTER*8 CMNAME
       DIMENSION STRESS(NTENS),STATEV(NSTATEV),
      1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
      2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
      3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
       
       ! Dynamic arrays based on NPT (number of Prony terms)
       REAL*8, ALLOCATABLE :: G_i(:), TAU_i(:), M_i(:)
       REAL*8, ALLOCATABLE :: SM_OLD(:,:), SM(:,:), SM_DOT(:,:)
       REAL*8, ALLOCATABLE :: DEV_STRAIN(:), DEV_DSTRAIN(:)
       DIMENSION G(6,6)
       REAL*8 VOL_STRAIN, VOL_DSTRAIN, HYDRO_STRESS
       REAL*8 G_INF
       
 !.. Material Properties
       E0 = PROPS(1)    ! Instantaneous Young's Modulus
       NU = PROPS(2)    ! Poisson's Ratio
       NPT = INT(PROPS(3))  ! Number of Prony series terms
       
       ! Allocate arrays based on NPT
       ALLOCATE(G_i(NPT))
       ALLOCATE(TAU_i(NPT))
       ALLOCATE(M_i(NPT))
       ALLOCATE(SM_OLD(NPT,6))
       ALLOCATE(SM(NPT,6))
       ALLOCATE(SM_DOT(NPT,6))
       ALLOCATE(DEV_STRAIN(6))
       ALLOCATE(DEV_DSTRAIN(6))
       
       ! Read g_i and tau_i values
       DO I = 1, NPT
           G_i(I) = PROPS(3+2*I-1)    ! g_i values
           TAU_i(I) = PROPS(3+2*I)    ! tau_i values
       END DO
       
 !.. Calculate G_infinity (long-term shear modulus ratio)
       G_INF = 1.0D0
       DO I = 1, NPT
           G_INF = G_INF - G_i(I)
       END DO
       
 !.. Calculate instantaneous shear and bulk moduli
       MU0 = E0/(2.0D0*(1.0D0 + NU))  ! Instantaneous shear modulus
       K0 = E0/(3.0D0*(1.0D0 - 2.0D0*NU))  ! Instantaneous bulk modulus
       
 !.. Read state variables
       DO N = 1, NPT
           DO I = 1, 6
               SM_OLD(N,I) = STATEV((N-1)*6 + I)
           END DO
       END DO
       
 !.. Calculate M_i terms
       DO N = 1, NPT
           M_i(N) = (TAU_i(N)*G_i(N)*MU0 - 
      1    TAU_i(N)*G_i(N)*MU0*EXP(-DTIME/TAU_i(N)))/(MU0*DTIME)
       END DO
       
 !.. Decompose strain into volumetric and deviatoric parts
       ! Calculate volumetric strain
       VOL_STRAIN = (STRAN(1) + STRAN(2) + STRAN(3))/3.D0
       VOL_DSTRAIN = (DSTRAN(1) + DSTRAN(2) + DSTRAN(3))/3.D0
       
       ! Calculate deviatoric strain
       DO I = 1, 3
           DEV_STRAIN(I) = STRAN(I) - VOL_STRAIN
           DEV_DSTRAIN(I) = DSTRAN(I) - VOL_DSTRAIN
       END DO
       DO I = 4, 6
           DEV_STRAIN(I) = STRAN(I)
           DEV_DSTRAIN(I) = DSTRAN(I)
       END DO
       
 !.. Calculate hydrostatic stress (elastic only)
       HYDRO_STRESS = 3.D0*K0*(VOL_STRAIN + VOL_DSTRAIN)
       
 !.. Calculate total M term for deviatoric part
       TOTAL_M = 1.D0
       DO N = 1, NPT
           TOTAL_M = TOTAL_M + M_i(N)
       END DO
       
 !.. Calculate stresses
       ! Normal stresses (11, 22, 33)
       DO I = 1, 3
           ! Hydrostatic contribution (elastic)
           STRESS(I) = HYDRO_STRESS
           
           ! Long-term elastic deviatoric contribution
           STRESS(I) = STRESS(I) + 2.D0*MU0*G_INF*DEV_STRAIN(I)
           
           ! Add viscoelastic history contribution
           DO N = 1, NPT
               STRESS(I) = STRESS(I) + SM_OLD(N,I)
           END DO
           
           ! Add viscoelastic increment contribution
           STRESS(I) = STRESS(I) + 
      1      2.D0*MU0*TOTAL_M*DEV_DSTRAIN(I)
       END DO
       
       ! Shear stresses (12, 23, 13)
       DO I = 4, 6
           ! Long-term elastic shear contribution
           STRESS(I) = 2.D0*MU0*G_INF*DEV_STRAIN(I)
           
           ! Add viscoelastic history contribution
           DO N = 1, NPT
               STRESS(I) = STRESS(I) + SM_OLD(N,I)
           END DO
           
           ! Add viscoelastic increment contribution
           STRESS(I) = STRESS(I) + 
      1      2.D0*MU0*TOTAL_M*DEV_DSTRAIN(I)
       END DO
       
 !.. Update internal state variables
       DO N = 1, NPT
           ! Normal components
           DO I = 1, 3
               SM(N,I) = SM_OLD(N,I) + 
      1          2.D0*MU0*M_i(N)*DEV_DSTRAIN(I)
           END DO
           
           ! Shear components
           DO I = 4, 6
               SM(N,I) = SM_OLD(N,I) + 
      1          2.D0*MU0*M_i(N)*DEV_DSTRAIN(I)
           END DO
       END DO
       
 !.. Calculate final state variables
       DO N = 1, NPT
           DO I = 1, 6
               SM_DOT(N,I) = EXP(-DTIME/TAU_i(N))*SM(N,I)
               STATEV((N-1)*6 + I) = SM_DOT(N,I)
           END DO
       END DO
       
 !.. Create Jacobian matrix
       ! First, set everything to zero
       DO K1 = 1, NTENS
           DO K2 = 1, NTENS
               DDSDDE(K2,K1) = 0.D0
           END DO
       END DO
       
       ! Add bulk (volumetric) contribution
       DO I = 1, 3
           DO J = 1, 3
               DDSDDE(I,J) = K0
           END DO
       END DO
       
       ! Add shear (deviatoric) contribution
       DO I = 1, 3
           DDSDDE(I,I) = DDSDDE(I,I) + 2.D0*MU0*TOTAL_M
       END DO
       
       ! Add shear terms
       DO I = 4, 6
           DDSDDE(I,I) = 2.D0*MU0*TOTAL_M
       END DO
       
 !.. Deallocate arrays
       DEALLOCATE(G_i, TAU_i, M_i, SM_OLD, SM, SM_DOT)
       DEALLOCATE(DEV_STRAIN, DEV_DSTRAIN)
       
       RETURN
       END