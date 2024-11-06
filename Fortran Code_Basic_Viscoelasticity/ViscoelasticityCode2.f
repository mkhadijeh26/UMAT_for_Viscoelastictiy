!---------------------------------------------------------------------------
! UMAT Implementation for Viscoelastic Material Model with Prony Series
!---------------------------------------------------------------------------
! Input Properties (PROPS array):
! PROPS(1) - E0: Instantaneous Young's Modulus
! PROPS(2) - NU: Poisson's Ratio  
! PROPS(3) - NPT: Number of Prony series terms
! PROPS(4:end) - Alternating g_i and tau_i values for each Prony term
!
! State Variables (STATEV array):
! Stores internal variables for stress history of each Prony term
! Size = NPT * 6 (6 stress components per Prony term)
!---------------------------------------------------------------------------

      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
        1 RPL,DDSDDT,DRPLDE,DRPLDT,
        2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
        3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
        4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
   
         INCLUDE 'ABA_PARAM.INC'
         
   !---------------------------------------------------------------------------
   ! 1. Variable Declarations
   !---------------------------------------------------------------------------
         CHARACTER*80 CMNAME
         DIMENSION STRESS(NTENS),STATEV(NSTATV),
        1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
        2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
        3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
         
         ! Arrays for Prony series and internal variables
         REAL*8, ALLOCATABLE :: G_i(:), TAU_i(:)
         REAL*8, ALLOCATABLE :: EVN(:,:), EN(:)
         REAL*8, ALLOCATABLE :: DEN(:,:), DE(:,:)
         
         ! Arrays for strain and stress components
         REAL*8, ALLOCATABLE :: DEV_STRAIN(:), DEV_DSTRAIN(:)
         
         ! Material constants
         REAL*8 E0, NU, K0, G0, G_INF, GT, KT
         REAL*8 VOL_STRAIN, VOL_DSTRAIN, HYDRO_PRESSURE
         INTEGER NUM_PRONY_TERMS
         
   !---------------------------------------------------------------------------
   ! 2. Initialize Arrays and Material Properties
   !---------------------------------------------------------------------------
         ! Read basic material properties
         E0 = PROPS(1)    ! Instantaneous Young's Modulus
         NU = PROPS(2)    ! Poisson's Ratio
         NUM_PRONY_TERMS = INT(PROPS(3))  ! Number of Prony series terms
         
         ! Calculate elastic constants
         G0 = E0/(2.0D0*(1.0D0 + NU))  ! Instantaneous shear modulus
         K0 = E0/(3.0D0*(1.0D0 - 2.0D0*NU))  ! Bulk modulus
         
         ! Allocate arrays
         ALLOCATE(G_i(NUM_PRONY_TERMS), TAU_i(NUM_PRONY_TERMS))
         ALLOCATE(EVN(NUM_PRONY_TERMS,NTENS))
         ALLOCATE(EN(NTENS))
         ALLOCATE(DEN(NUM_PRONY_TERMS,NTENS))
         ALLOCATE(DE(NUM_PRONY_TERMS,NTENS))
         ALLOCATE(DEV_STRAIN(NTENS), DEV_DSTRAIN(NTENS))
         
         ! Read Prony series parameters
         DO I = 1, NUM_PRONY_TERMS
             G_i(I) = PROPS(3+2*I-1)    ! Relative moduli
             TAU_i(I) = PROPS(3+2*I)    ! Relaxation times
         END DO
         
   !---------------------------------------------------------------------------
   ! 3. Process State Variables
   !---------------------------------------------------------------------------
         ! Read previous state variables
         DO I = 1, NTENS
             EN(I) = STATEV(I)
         END DO
         
         DO N = 1, NUM_PRONY_TERMS
             DO I = 1, NTENS
                 EVN(N,I) = STATEV(NTENS*N + I)
             END DO
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
         DO I = 4, NTENS
             DEV_STRAIN(I) = STRAN(I)/2.D0
             DEV_DSTRAIN(I) = DSTRAN(I)/2.D0
         END DO
         
   !---------------------------------------------------------------------------
   ! 5. Calculate Internal Variables
   !---------------------------------------------------------------------------
         ! Calculate DEN terms for each Prony component
         DO N = 1, NUM_PRONY_TERMS
             DO I = 1, NTENS
                 DEN(N,I) = DEV_DSTRAIN(I)/DTIME*(DTIME-TAU_i(N)*
        1                   (1-EXP(-DTIME/TAU_i(N))))
             END DO
         END DO
         
         ! Calculate viscous strain increments
         DO N = 1, NUM_PRONY_TERMS
             DO I = 1, NTENS
                 DE(N,I) = (1-EXP(-DTIME/TAU_i(N)))*EN(I)
        1                  -(1-EXP(-DTIME/TAU_i(N)))*EVN(N,I)
        2                  + DEN(N,I)
             END DO
         END DO
         
   !---------------------------------------------------------------------------
   ! 6. Stress Update
   !---------------------------------------------------------------------------
         ! Calculate deviatoric stress increment
         DO I = 1, NTENS
             DGSTRES = 2.D0*G0*DEV_DSTRAIN(I)
             DO N = 1, NUM_PRONY_TERMS
                 DGSTRES = DGSTRES - 2.D0*G0*G_i(N)*DE(N,I)
             END DO
             STRESS(I) = STRESS(I) + DGSTRES
         END DO
         
         ! Add hydrostatic pressure (volumetric contribution)
         HYDRO_PRESSURE = 3.D0*K0*VOL_DSTRAIN
         DO I = 1, NDI
             STRESS(I) = STRESS(I) + HYDRO_PRESSURE
         END DO
         
   !---------------------------------------------------------------------------
   ! 7. Update Material Jacobian
   !---------------------------------------------------------------------------
         ! Calculate effective moduli
         GT = G0
         DO N = 1, NUM_PRONY_TERMS
             GT = GT - G0*G_i(N)*(1-TAU_i(N)/DTIME*
        1         (1-EXP(-DTIME/TAU_i(N))))
         END DO
         
         KT = K0
         
         ! Build Jacobian matrix
         DO I = 1, NTENS
             DO J = 1, NTENS
                 DDSDDE(I,J) = 0.D0
             END DO
         END DO
         
         DO I = 1, NDI
             DO J = 1, NDI
                 DDSDDE(I,J) = KT - 2.D0*GT/3.D0
             END DO
             DDSDDE(I,I) = DDSDDE(I,I) + 2.D0*GT
         END DO
         
         DO I = NDI+1, NTENS
             DDSDDE(I,I) = GT
         END DO
         
   !---------------------------------------------------------------------------
   ! 8. Update State Variables
   !---------------------------------------------------------------------------
         ! Update EN
         DO I = 1, NTENS
             STATEV(I) = EN(I) + DEV_DSTRAIN(I)
         END DO
         
         ! Update EVN for each Prony term
         DO N = 1, NUM_PRONY_TERMS
             DO I = 1, NTENS
                 STATEV(NTENS*N + I) = EVN(N,I) + DE(N,I)
             END DO
         END DO
         
   !---------------------------------------------------------------------------
   ! 9. Cleanup
   !---------------------------------------------------------------------------
         DEALLOCATE(G_i, TAU_i)
         DEALLOCATE(EVN, EN)
         DEALLOCATE(DEN, DE)
         DEALLOCATE(DEV_STRAIN, DEV_DSTRAIN)
         
         RETURN
         END