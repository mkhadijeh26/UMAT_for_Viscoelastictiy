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
   ! 1. Variable Declarations
   !---------------------------------------------------------------------------
         CHARACTER*8 CMNAME
         DIMENSION STRESS(NTENS),STATEV(NSTATEV),
        1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
        2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
        3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
         
         ! Arrays for Prony series
         REAL*8, ALLOCATABLE :: G_i(:), TAU_i(:), M_i(:)
         REAL*8, ALLOCATABLE :: SM_OLD(:,:), SM(:,:), SM_DOT(:,:)
         
         ! Arrays for strain components
         REAL*8, ALLOCATABLE :: DEV_STRAIN(:), DEV_DSTRAIN(:)
         
         ! Arrays for stress components
         REAL*8, ALLOCATABLE :: STRESS_HYDROSTATIC(:)
         REAL*8, ALLOCATABLE :: STRESS_LONGTERM(:)
         REAL*8, ALLOCATABLE :: STRESS_HISTORY(:)
         REAL*8, ALLOCATABLE :: STRESS_INSTANTANEOUS(:)
         
         ! Material constants
         REAL*8 E0, NU, K0, MU0, G_INF, TOTAL_M
         REAL*8 VOL_STRAIN, VOL_DSTRAIN
         INTEGER NUM_PRONY_TERMS
         
   !---------------------------------------------------------------------------
   ! 2. Initialize Arrays and Material Properties
   !---------------------------------------------------------------------------
         ! Read basic material properties
         E0 = PROPS(1)    ! Instantaneous Young's Modulus
         NU = PROPS(2)    ! Poisson's Ratio
         NUM_PRONY_TERMS = INT(PROPS(3))  ! Number of Prony series terms
         
         ! Allocate all arrays
         ALLOCATE(G_i(NUM_PRONY_TERMS), TAU_i(NUM_PRONY_TERMS))
         ALLOCATE(M_i(NUM_PRONY_TERMS))
         ALLOCATE(SM_OLD(NUM_PRONY_TERMS,6))
         ALLOCATE(SM(NUM_PRONY_TERMS,6))
         ALLOCATE(SM_DOT(NUM_PRONY_TERMS,6))
         ALLOCATE(DEV_STRAIN(6), DEV_DSTRAIN(6))
         ALLOCATE(STRESS_HYDROSTATIC(3))
         ALLOCATE(STRESS_LONGTERM(NTENS))
         ALLOCATE(STRESS_HISTORY(NTENS))
         ALLOCATE(STRESS_INSTANTANEOUS(NTENS))
         
         ! Read Prony series parameters
         DO I = 1, NUM_PRONY_TERMS
             G_i(I) = PROPS(3+2*I-1)    ! Relative moduli
             TAU_i(I) = PROPS(3+2*I)    ! Relaxation times
         END DO
         
   !---------------------------------------------------------------------------
   ! 3. Calculate Material Constants
   !---------------------------------------------------------------------------
         ! Calculate elastic constants
         MU0 = E0/(2.0D0*(1.0D0 + NU))  ! Instantaneous shear modulus
         K0 = E0/(3.0D0*(1.0D0 - 2.0D0*NU))  ! Bulk modulus
         
         ! Calculate long-term shear modulus ratio
         G_INF = 1.0D0
         DO I = 1, NUM_PRONY_TERMS
             G_INF = G_INF - G_i(I)
         END DO
         
         ! Calculate M_i terms for viscoelastic contribution
         DO N = 1, NUM_PRONY_TERMS
             M_i(N) = (TAU_i(N)*G_i(N)*MU0 - 
        1        TAU_i(N)*G_i(N)*MU0*EXP(-DTIME/TAU_i(N)))/(MU0*DTIME)
         END DO
         
         ! Calculate total M for instantaneous response
         TOTAL_M = 1.D0
         DO N = 1, NUM_PRONY_TERMS
             TOTAL_M = TOTAL_M + M_i(N)
         END DO
         
   !---------------------------------------------------------------------------
   ! 4. Process State Variables
   !---------------------------------------------------------------------------
         ! Read previous state variables
         DO N = 1, NUM_PRONY_TERMS
             DO I = 1, 6
                 SM_OLD(N,I) = STATEV((N-1)*6 + I)
             END DO
         END DO
         
   !---------------------------------------------------------------------------
   ! 5. Strain Decomposition
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
   ! 6. Stress Calculation
   !---------------------------------------------------------------------------
   !------------------------
   ! 6.1 Hydrostatic Stress
   !------------------------
         ! Pure elastic volumetric response (only for normal components 1,2,3)
         HYDRO_PRESSURE = 3.D0 * K0 * (VOL_STRAIN + VOL_DSTRAIN)
         DO I = 1, 3
             STRESS_HYDROSTATIC(I) = HYDRO_PRESSURE
         END DO
         
   !------------------------
   ! 6.2 Long-term Elastic Stress
   !------------------------
         ! Calculate for all components (both normal and shear)
         DO I = 1, NTENS
             STRESS_LONGTERM(I) = 2.D0 * MU0 * G_INF * DEV_STRAIN(I)
         END DO
         
   !------------------------
   ! 6.3 Viscoelastic History Stress
   !------------------------
         ! Initialize history stress
         DO I = 1, NTENS
             STRESS_HISTORY(I) = 0.D0
         END DO
         
         ! Sum up contributions from all Prony terms
         DO N = 1, NUM_PRONY_TERMS
             DO I = 1, NTENS
                 STRESS_HISTORY(I) = STRESS_HISTORY(I) + SM_OLD(N,I)
             END DO
         END DO
         
   !------------------------
   ! 6.4 Instantaneous Stress Response
   !------------------------
         ! Calculate for all components
         DO I = 1, NTENS
             STRESS_INSTANTANEOUS(I) = 2.D0 * MU0 * TOTAL_M * DEV_DSTRAIN(I)
         END DO
         
   !------------------------
   ! 6.5 Combine All Components
   !------------------------
         ! Normal components (1,2,3): Include all terms
         DO I = 1, 3
             STRESS(I) = STRESS_HYDROSTATIC(I) +    ! Hydrostatic (volumetric)
        1                STRESS_LONGTERM(I) +        ! Long-term elastic
        2                STRESS_HISTORY(I) +         ! Viscoelastic history
        3                STRESS_INSTANTANEOUS(I)     ! Instantaneous response
         END DO
         
         ! Shear components (4,5,6): No hydrostatic term
         DO I = 4, 6
             STRESS(I) = STRESS_LONGTERM(I) +        ! Long-term elastic
        1                STRESS_HISTORY(I) +         ! Viscoelastic history
        2                STRESS_INSTANTANEOUS(I)     ! Instantaneous response
         END DO
         
   !---------------------------------------------------------------------------
   ! 7. Update Internal Variables
   !---------------------------------------------------------------------------
         DO N = 1, NUM_PRONY_TERMS
             DO I = 1, 6
                 ! Calculate new internal stress variable
                 SM(N,I) = SM_OLD(N,I) + 2.D0*MU0*M_i(N)*DEV_DSTRAIN(I)
                 
                 ! Apply time-dependent relaxation
                 SM_DOT(N,I) = EXP(-DTIME/TAU_i(N))*SM(N,I)
                 
                 ! Store for next increment
                 STATEV((N-1)*6 + I) = SM_DOT(N,I)
             END DO
         END DO
         
   !---------------------------------------------------------------------------
   ! 8. Calculate Material Jacobian
   !---------------------------------------------------------------------------
         ! Initialize Jacobian matrix
         DO I = 1, NTENS
             DO J = 1, NTENS
                 DDSDDE(I,J) = 0.D0
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
   ! 9. Cleanup
   !---------------------------------------------------------------------------
         DEALLOCATE(G_i, TAU_i, M_i)
         DEALLOCATE(SM_OLD, SM, SM_DOT)
         DEALLOCATE(DEV_STRAIN, DEV_DSTRAIN)
         DEALLOCATE(STRESS_HYDROSTATIC)
         DEALLOCATE(STRESS_LONGTERM)
         DEALLOCATE(STRESS_HISTORY)
         DEALLOCATE(STRESS_INSTANTANEOUS)
         
         RETURN
         END