SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
    1 RPL,DDSDDT,DRPLDE,DRPLDT,
    2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
    3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
    4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

!-----------------------------------------------------------------------
! Description:
! This UMAT implements a linear viscoelastic material model with aging effects.
! The aging is modeled through a carbonyl area-based aging index that 
! modifies the material stiffness.
!
! Input Parameters in PROPS array:
! PROPS(1) - E0: Initial Young's modulus
! PROPS(2) - NU: Poisson's ratio
! PROPS(3) - N: Number of Prony series terms
! PROPS(4:3+2N) - Pairs of (G_i, TAU_i) for each Prony term
! PROPS(4+2N) - k: Aging scaling factor
! PROPS(5+2N) - CA0: Initial carbonyl area
! PROPS(6+2N) - CAINF: Long-term carbonyl area
! PROPS(7+2N) - RT: Reaction rate
!-----------------------------------------------------------------------

     INCLUDE 'ABA_PARAM.INC'
     
     CHARACTER*80 CMNAME
     DIMENSION STRESS(NTENS),STATEV(NSTATV),
    1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
    2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
    3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
     
     ! Arrays for Prony series and internal variables
     REAL*8, ALLOCATABLE :: G_i(:), TAU_i(:), H_i(:)
     REAL*8, ALLOCATABLE :: SM_OLD(:,:), SM_NEW(:,:)
     REAL*8, ALLOCATABLE :: DEV_STRAIN(:), DEV_DSTRAIN(:)
     REAL*8, ALLOCATABLE :: STRESS_VOL(:), STRESS_DEV(:)
     
     ! Material constants and computational variables
     REAL*8 E0, E_AGED, NU, K0, G0, G_INF
     REAL*8 VOL_STRAIN, VOL_DSTRAIN
     REAL*8 ALPHA, BETA, EXP_TERM
     REAL*8 TOTAL_G_EFF, DELTA_IJ
     INTEGER NUM_PRONY_TERMS
     
     ! Aging parameters
     REAL*8 K_AGING        ! Scaling factor k
     REAL*8 CA_0          ! Initial carbonyl area
     REAL*8 CA_INF        ! Long-term carbonyl area
     REAL*8 R_T           ! Reaction rate
     REAL*8 CURRENT_TIME  ! Current simulation time
     REAL*8 CA_T          ! Current carbonyl area
     REAL*8 AGING_INDEX   ! Current aging index A(t)
     
!---------------------------------------------------------------------------
! Initialize Arrays and Material Properties
!---------------------------------------------------------------------------
     ! Read basic material properties
     E0 = PROPS(1)    ! Initial Young's Modulus
     NU = PROPS(2)    ! Poisson's Ratio
     NUM_PRONY_TERMS = INT(PROPS(3))  ! Number of Prony terms
     
     ! Read aging parameters
     K_AGING = PROPS(3 + 2*NUM_PRONY_TERMS + 1)  ! Scaling factor k
     CA_0 = PROPS(3 + 2*NUM_PRONY_TERMS + 2)     ! Initial CA
     CA_INF = PROPS(3 + 2*NUM_PRONY_TERMS + 3)   ! Long-term CA
     R_T = PROPS(3 + 2*NUM_PRONY_TERMS + 4)      ! Reaction rate
     
     ! Calculate current time and aging parameters
     CURRENT_TIME = TIME(1)
     
     ! Calculate current carbonyl area CA(t)
     CA_T = CA_0 + CA_INF*(1.0D0 - EXP(-R_T*CURRENT_TIME))
     
     ! Calculate aging index A(t)
     AGING_INDEX = 1.0D0 + K_AGING*CA_T
     
     ! Calculate aged Young's modulus
     E_AGED = E0*AGING_INDEX
     
     ! Calculate elastic constants with aged modulus
     G0 = E_AGED/(2.0D0*(1.0D0 + NU))    ! Initial shear modulus
     K0 = E_AGED/(3.0D0*(1.0D0 - 2.0D0*NU))  ! Bulk modulus
     
     ! Allocate arrays based on problem size
     ALLOCATE(G_i(NUM_PRONY_TERMS))
     ALLOCATE(TAU_i(NUM_PRONY_TERMS))
     ALLOCATE(H_i(NUM_PRONY_TERMS))
     ALLOCATE(SM_OLD(NUM_PRONY_TERMS,NTENS))
     ALLOCATE(SM_NEW(NUM_PRONY_TERMS,NTENS))
     ALLOCATE(DEV_STRAIN(NTENS))
     ALLOCATE(DEV_DSTRAIN(NTENS))
     ALLOCATE(STRESS_VOL(3))
     ALLOCATE(STRESS_DEV(NTENS))
     
     ! Initialize and process Prony series parameters
     G_INF = 1.0D0    ! Initialize long-term modulus ratio
     TOTAL_G_EFF = 0.0D0
     
     DO I = 1, NUM_PRONY_TERMS
         G_i(I) = PROPS(3+2*I-1)    ! Relative moduli
         TAU_i(I) = PROPS(3+2*I)    ! Relaxation times
         G_INF = G_INF - G_i(I)     ! Update long-term modulus
         
         ! Calculate integration parameters
         ALPHA = DTIME/TAU_i(I)
         BETA = 1.0D0/(1.0D0 + ALPHA)
         
         ! Calculate effective modulus contribution
         H_i(I) = G_i(I)*BETA
         TOTAL_G_EFF = TOTAL_G_EFF + H_i(I)
     END DO
     
!---------------------------------------------------------------------------
! Process State Variables and Calculate Strains
!---------------------------------------------------------------------------
     ! Retrieve internal variables from previous increment
     DO N = 1, NUM_PRONY_TERMS
         DO I = 1, NTENS
             SM_OLD(N,I) = STATEV((N-1)*NTENS + I)
         END DO
     END DO
     
     ! Calculate volumetric strain components
     VOL_STRAIN = (STRAN(1) + STRAN(2) + STRAN(3))/3.0D0
     VOL_DSTRAIN = (DSTRAN(1) + DSTRAN(2) + DSTRAN(3))/3.0D0
     
     ! Calculate deviatoric strain components
     DO I = 1, 3
         DEV_STRAIN(I) = STRAN(I) - VOL_STRAIN
         DEV_DSTRAIN(I) = DSTRAN(I) - VOL_DSTRAIN
     END DO
     DO I = 4, NTENS
         DEV_STRAIN(I) = STRAN(I)/2.0D0
         DEV_DSTRAIN(I) = DSTRAN(I)/2.0D0
     END DO
     
!---------------------------------------------------------------------------
! Calculate Stresses
!---------------------------------------------------------------------------
     ! Calculate volumetric stress (purely elastic)
     DO I = 1, 3
         STRESS_VOL(I) = 3.0D0*K0*(VOL_STRAIN + VOL_DSTRAIN)
     END DO
     
     ! Initialize deviatoric stress with long-term response
     DO I = 1, NTENS
         STRESS_DEV(I) = 2.0D0*G0*G_INF*(DEV_STRAIN(I) + DEV_DSTRAIN(I))
     END DO
     
     ! Add viscoelastic contributions from each Prony term
     DO N = 1, NUM_PRONY_TERMS
         ALPHA = DTIME/TAU_i(N)
         BETA = 1.0D0/(1.0D0 + ALPHA)
         EXP_TERM = EXP(-ALPHA)
         
         DO I = 1, NTENS
             ! Update internal variables using exponential integration
             SM_NEW(N,I) = EXP_TERM*SM_OLD(N,I) + 
    1                      2.0D0*G0*G_i(N)*BETA*DEV_DSTRAIN(I)
             
             ! Add contribution to total deviatoric stress
             STRESS_DEV(I) = STRESS_DEV(I) + SM_NEW(N,I)
         END DO
     END DO
     
     ! Combine volumetric and deviatoric stresses
     DO I = 1, 3
         STRESS(I) = STRESS_VOL(I) + STRESS_DEV(I)
     END DO
     DO I = 4, NTENS
         STRESS(I) = STRESS_DEV(I)
     END DO
     
!---------------------------------------------------------------------------
! Calculate Material Jacobian (Tangent Stiffness Matrix)
!---------------------------------------------------------------------------
     ! Initialize Jacobian
     DO I = 1, NTENS
         DO J = 1, NTENS
             DDSDDE(I,J) = 0.0D0
         END DO
     END DO
     
     ! Calculate effective shear modulus
     TOTAL_G_EFF = G0*(G_INF + TOTAL_G_EFF)
     
     ! Build material Jacobian
     DO I = 1, 3
         DO J = 1, 3
             ! Set Kronecker delta
             IF (I .EQ. J) THEN
                 DELTA_IJ = 1.0D0
             ELSE
                 DELTA_IJ = 0.0D0
             END IF
             
             ! Combined volumetric and deviatoric response
             DDSDDE(I,J) = K0 - 2.0D0*TOTAL_G_EFF/3.0D0 + 
    1                      2.0D0*TOTAL_G_EFF*DELTA_IJ
         END DO
     END DO
     
     ! Set shear components
     DO I = 4, NTENS
         DDSDDE(I,I) = TOTAL_G_EFF
     END DO
     
!---------------------------------------------------------------------------
! Update State Variables
!---------------------------------------------------------------------------
     ! Update Prony series internal variables
     DO N = 1, NUM_PRONY_TERMS
         DO I = 1, NTENS
             STATEV((N-1)*NTENS + I) = SM_NEW(N,I)
         END DO
     END DO
     
     ! Store aging-related variables in state variables
     STATEV(NSTATV-3) = CA_T          ! Current carbonyl area
     STATEV(NSTATV-2) = AGING_INDEX   ! Current aging index
     STATEV(NSTATV-1) = E_AGED        ! Current aged Young's modulus
     STATEV(NSTATV) = CURRENT_TIME    ! Current time
     
!---------------------------------------------------------------------------
! Cleanup
!---------------------------------------------------------------------------
     DEALLOCATE(G_i, TAU_i, H_i)
     DEALLOCATE(SM_OLD, SM_NEW)
     DEALLOCATE(DEV_STRAIN, DEV_DSTRAIN)
     DEALLOCATE(STRESS_VOL, STRESS_DEV)
     
     RETURN
     END