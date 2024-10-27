!---------------------------------------------------------------------------
! UMAT Implementation for Viscoelastic Material Model with Prony Series
!---------------------------------------------------------------------------
! This subroutine implements a viscoelastic material model using Prony series
! for ABAQUS/Standard. The model accounts for both volumetric (elastic) and 
! deviatoric (viscoelastic) deformation.
!
! Input Properties (PROPS array):
! PROPS(1) - E0: Instantaneous Young's Modulus
! PROPS(2) - NU: Poisson's Ratio
! PROPS(3) - NPT: Number of Prony series terms
! PROPS(4:end) - Alternating g_i and tau_i values for each Prony term
!
! State Variables (STATEV array):
! Stores the stress-like internal variables (SM) for each Prony term
! Size = NPT * 6 (6 stress components per Prony term)
!---------------------------------------------------------------------------

      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
        1 RPL,DDSDDT,DRPLDE,DRPLDT,
        2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
        3 NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
        4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
   
         INCLUDE 'ABA_PARAM.INC'
         
   !---------------------------------------------------------------------------
   ! Variable Declarations
   !---------------------------------------------------------------------------
         CHARACTER*8 CMNAME
         DIMENSION STRESS(NTENS),STATEV(NSTATEV),
        1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
        2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
        3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
         
         ! Dynamic arrays for Prony series terms
         REAL*8, ALLOCATABLE :: G_i(:), TAU_i(:), M_i(:)
         REAL*8, ALLOCATABLE :: SM_OLD(:,:), SM(:,:), SM_DOT(:,:)
         REAL*8, ALLOCATABLE :: DEV_STRAIN(:), DEV_DSTRAIN(:)
         DIMENSION G(6,6)
         
         ! Additional variables
         REAL*8 VOL_STRAIN, VOL_DSTRAIN, HYDRO_STRESS
         REAL*8 G_INF, MU0, K0, TOTAL_M
         
   !---------------------------------------------------------------------------
   ! 1. Initialize Material Properties
   !---------------------------------------------------------------------------
         ! Read basic material properties
         E0 = PROPS(1)    ! Instantaneous Young's Modulus
         NU = PROPS(2)    ! Poisson's Ratio
         NPT = INT(PROPS(3))  ! Number of Prony series terms
         
         ! Allocate dynamic arrays
         ALLOCATE(G_i(NPT), TAU_i(NPT), M_i(NPT))
         ALLOCATE(SM_OLD(NPT,6), SM(NPT,6), SM_DOT(NPT,6))
         ALLOCATE(DEV_STRAIN(6), DEV_DSTRAIN(6))
         
         ! Read Prony series parameters
         DO I = 1, NPT
             G_i(I) = PROPS(3+2*I-1)    ! Relative moduli
             TAU_i(I) = PROPS(3+2*I)    ! Relaxation times
         END DO
         
   !---------------------------------------------------------------------------
   ! 2. Calculate Material Constants
   !---------------------------------------------------------------------------
         ! Calculate long-term shear modulus ratio
         G_INF = 1.0D0
         DO I = 1, NPT
             G_INF = G_INF - G_i(I)
         END DO
         
         ! Calculate elastic constants
         MU0 = E0/(2.0D0*(1.0D0 + NU))  ! Instantaneous shear modulus
         K0 = E0/(3.0D0*(1.0D0 - 2.0D0*NU))  ! Bulk modulus
         
   !---------------------------------------------------------------------------
   ! 3. Process State Variables and Calculate M_i Terms
   !---------------------------------------------------------------------------
         ! Read previous state
         DO N = 1, NPT
             DO I = 1, 6
                 SM_OLD(N,I) = STATEV((N-1)*6 + I)
             END DO
         END DO
         
         ! Calculate M_i terms for viscoelastic contribution
         DO N = 1, NPT
             M_i(N) = (TAU_i(N)*G_i(N)*MU0 - 
        1       TAU_i(N)*G_i(N)*MU0*EXP(-DTIME/TAU_i(N)))/(MU0*DTIME)
         END DO
         
   !---------------------------------------------------------------------------
   ! 4. Strain Decomposition
   !---------------------------------------------------------------------------
         ! Calculate volumetric components
         VOL_STRAIN = (STRAN(1) + STRAN(2) + STRAN(3))/3.D0
         VOL_DSTRAIN = (DSTRAN(1) + DSTRAN(2) + DSTRAN(3))/3.D0
         
         ! Calculate deviatoric components
         DO I = 1, 3
             DEV_STRAIN(I) = STRAN(I) - VOL_STRAIN
             DEV_DSTRAIN(I) = DSTRAN(I) - VOL_DSTRAIN
         END DO
         DO I = 4, 6
             DEV_STRAIN(I) = STRAN(I)
             DEV_DSTRAIN(I) = DSTRAN(I)
         END DO
         
   !---------------------------------------------------------------------------
   ! 5. Stress Calculation
   !---------------------------------------------------------------------------
         ! Calculate hydrostatic stress (elastic)
         HYDRO_STRESS = 3.D0*K0*(VOL_STRAIN + VOL_DSTRAIN)
         
         ! Calculate total M term for deviatoric part
         TOTAL_M = 1.D0
         DO N = 1, NPT
             TOTAL_M = TOTAL_M + M_i(N)
         END DO
         
         ! Calculate normal stresses (11, 22, 33)
         DO I = 1, 3
             STRESS(I) = HYDRO_STRESS +  ! Hydrostatic part
        1                2.D0*MU0*G_INF*DEV_STRAIN(I)  ! Long-term elastic
             
             ! Add viscoelastic contributions
             DO N = 1, NPT
                 STRESS(I) = STRESS(I) + SM_OLD(N,I)
             END DO
             STRESS(I) = STRESS(I) + 2.D0*MU0*TOTAL_M*DEV_DSTRAIN(I)
         END DO
         
         ! Calculate shear stresses (12, 23, 13)
         DO I = 4, 6
             STRESS(I) = 2.D0*MU0*G_INF*DEV_STRAIN(I)  ! Long-term elastic
             
             ! Add viscoelastic contributions
             DO N = 1, NPT
                 STRESS(I) = STRESS(I) + SM_OLD(N,I)
             END DO
             STRESS(I) = STRESS(I) + 2.D0*MU0*TOTAL_M*DEV_DSTRAIN(I)
         END DO
         
   !---------------------------------------------------------------------------
   ! 6. Update Internal Variables
   !---------------------------------------------------------------------------
         DO N = 1, NPT
             ! Update all components
             DO I = 1, 6
                 SM(N,I) = SM_OLD(N,I) + 2.D0*MU0*M_i(N)*DEV_DSTRAIN(I)
                 SM_DOT(N,I) = EXP(-DTIME/TAU_i(N))*SM(N,I)
                 STATEV((N-1)*6 + I) = SM_DOT(N,I)
             END DO
         END DO
         
   !---------------------------------------------------------------------------
   ! 7. Calculate Material Jacobian
   !---------------------------------------------------------------------------
         ! Initialize Jacobian matrix
         DO K1 = 1, NTENS
             DO K2 = 1, NTENS
                 DDSDDE(K2,K1) = 0.D0
             END DO
         END DO
         
         ! Add volumetric contribution
         DO I = 1, 3
             DO J = 1, 3
                 DDSDDE(I,J) = K0
             END DO
         END DO
         
         ! Add deviatoric contribution
         DO I = 1, 3
             DDSDDE(I,I) = DDSDDE(I,I) + 2.D0*MU0*TOTAL_M
         END DO
         
         ! Add shear terms
         DO I = 4, 6
             DDSDDE(I,I) = 2.D0*MU0*TOTAL_M
         END DO
         
   !---------------------------------------------------------------------------
   ! 8. Cleanup
   !---------------------------------------------------------------------------
         DEALLOCATE(G_i, TAU_i, M_i, SM_OLD, SM, SM_DOT)
         DEALLOCATE(DEV_STRAIN, DEV_DSTRAIN)
         
         RETURN
         END