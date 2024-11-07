      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

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
      REAL*8 E0, NU, K0, G0, G_INF
      REAL*8 VOL_STRAIN, VOL_DSTRAIN
      REAL*8 ALPHA, BETA, EXP_TERM
      REAL*8 TOTAL_G_EFF, DELTA_IJ
      INTEGER NUM_PRONY_TERMS
      
!---------------------------------------------------------------------------
! Initialize Arrays and Material Properties
!---------------------------------------------------------------------------
      ! Read material properties
      E0 = PROPS(1)    ! Initial Young's Modulus
      NU = PROPS(2)    ! Poisson's Ratio
      NUM_PRONY_TERMS = INT(PROPS(3))  ! Number of Prony terms
      
      ! Calculate elastic constants
      G0 = E0/(2.0D0*(1.0D0 + NU))    ! Initial shear modulus
      K0 = E0/(3.0D0*(1.0D0 - 2.0D0*NU))  ! Bulk modulus
      
      ! Allocate arrays
      ALLOCATE(G_i(NUM_PRONY_TERMS))
      ALLOCATE(TAU_i(NUM_PRONY_TERMS))
      ALLOCATE(H_i(NUM_PRONY_TERMS))
      ALLOCATE(SM_OLD(NUM_PRONY_TERMS,NTENS))
      ALLOCATE(SM_NEW(NUM_PRONY_TERMS,NTENS))
      ALLOCATE(DEV_STRAIN(NTENS))
      ALLOCATE(DEV_DSTRAIN(NTENS))
      ALLOCATE(STRESS_VOL(3))
      ALLOCATE(STRESS_DEV(NTENS))
      
      ! Initialize Prony series parameters
      G_INF = 1.0D0
      TOTAL_G_EFF = 0.0D0
      
      ! Read and process Prony series parameters
      DO I = 1, NUM_PRONY_TERMS
          G_i(I) = PROPS(3+2*I-1)    ! Relative moduli
          TAU_i(I) = PROPS(3+2*I)    ! Relaxation times
          G_INF = G_INF - G_i(I)
          
          ! Calculate integration parameters
          ALPHA = DTIME/TAU_i(I)
          BETA = 1.0D0/(1.0D0 + ALPHA)
          EXP_TERM = EXP(-ALPHA)
          
          ! Calculate effective modulus contribution
          H_i(I) = G_i(I)*BETA
          TOTAL_G_EFF = TOTAL_G_EFF + H_i(I)
      END DO
      
!---------------------------------------------------------------------------
! Process State Variables and Calculate Strains
!---------------------------------------------------------------------------
      ! Retrieve previous internal variables
      DO N = 1, NUM_PRONY_TERMS
          DO I = 1, NTENS
              SM_OLD(N,I) = STATEV((N-1)*NTENS + I)
          END DO
      END DO
      
      ! Calculate strain components
      VOL_STRAIN = (STRAN(1) + STRAN(2) + STRAN(3))/3.0D0
      VOL_DSTRAIN = (DSTRAN(1) + DSTRAN(2) + DSTRAN(3))/3.0D0
      
      ! Calculate deviatoric strains
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
      ! Volumetric stress (purely elastic)
      DO I = 1, 3
          STRESS_VOL(I) = 3.0D0*K0*(VOL_STRAIN + VOL_DSTRAIN)
      END DO
      
      ! Initialize deviatoric stress with long-term response
      DO I = 1, NTENS
          STRESS_DEV(I) = 2.0D0*G0*G_INF*(DEV_STRAIN(I) + DEV_DSTRAIN(I))
      END DO
      
      ! Add viscoelastic contributions
      DO N = 1, NUM_PRONY_TERMS
          ALPHA = DTIME/TAU_i(N)
          BETA = 1.0D0/(1.0D0 + ALPHA)
          EXP_TERM = EXP(-ALPHA)
          
          DO I = 1, NTENS
              ! Update internal variables with improved integration
              SM_NEW(N,I) = EXP_TERM*SM_OLD(N,I) + 
     1                      2.0D0*G0*G_i(N)*BETA*DEV_DSTRAIN(I)
              
              ! Add contribution to total deviatoric stress
              STRESS_DEV(I) = STRESS_DEV(I) + SM_NEW(N,I)
          END DO
      END DO
      
      ! Combine stresses
      DO I = 1, 3
          STRESS(I) = STRESS_VOL(I) + STRESS_DEV(I)
      END DO
      DO I = 4, NTENS
          STRESS(I) = STRESS_DEV(I)
      END DO
      
!---------------------------------------------------------------------------
! Calculate Material Jacobian
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
              ! Kronecker delta
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
      
      ! Shear components
      DO I = 4, NTENS
          DDSDDE(I,I) = TOTAL_G_EFF
      END DO
      
!---------------------------------------------------------------------------
! Update State Variables
!---------------------------------------------------------------------------
      DO N = 1, NUM_PRONY_TERMS
          DO I = 1, NTENS
              STATEV((N-1)*NTENS + I) = SM_NEW(N,I)
          END DO
      END DO
      
!---------------------------------------------------------------------------
! Cleanup
!---------------------------------------------------------------------------
      DEALLOCATE(G_i, TAU_i, H_i)
      DEALLOCATE(SM_OLD, SM_NEW)
      DEALLOCATE(DEV_STRAIN, DEV_DSTRAIN)
      DEALLOCATE(STRESS_VOL, STRESS_DEV)
      
      RETURN
      END