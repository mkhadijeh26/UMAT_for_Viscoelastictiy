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
        REAL*8, ALLOCATABLE :: DEV_STRAIN(:), DEV_DSTRAIN(:)
        REAL*8, ALLOCATABLE :: STRESS_HYDROSTATIC(:)
        REAL*8, ALLOCATABLE :: STRESS_LONGTERM(:)
        REAL*8, ALLOCATABLE :: STRESS_HISTORY(:)
        REAL*8, ALLOCATABLE :: STRESS_INSTANTANEOUS(:)
        
        ! Material constants
        REAL*8 E0, NU, K0, MU0, G_INF, TOTAL_M
        REAL*8 VOL_STRAIN, VOL_DSTRAIN
        REAL*8 CA0, CA_INF, R_T, K_AGING, CA_t, AGING_INDEX, E0_AGED
        INTEGER NUM_PRONY_TERMS
        
  !---------------------------------------------------------------------------
  ! 2. Calculate Aged Young's Modulus
  !---------------------------------------------------------------------------
        ! Read basic material properties
        E0 = PROPS(1)    ! Initial Young's Modulus
        NU = PROPS(2)    ! Poisson's Ratio
        NUM_PRONY_TERMS = INT(PROPS(3))  ! Number of Prony series terms
        
        ! Read aging parameters
        CA0 = PROPS(NPROPS-3)     ! Initial carbonyl area
        CA_INF = PROPS(NPROPS-2)  ! Long-term carbonyl area
        R_T = PROPS(NPROPS-1)     ! Reaction rate
        K_AGING = PROPS(NPROPS)   ! Aging scaling factor
        
        ! Calculate current carbonyl area and aging index
        CA_t = CA0 + CA_INF * (1.0D0 - EXP(-R_T * TIME(1)))
        AGING_INDEX = 1.0D0 + K_AGING * CA_t
        
        ! Calculate aged Young's modulus - THIS IS THE KEY LINE
        E0_AGED = E0 * AGING_INDEX
        
  !---------------------------------------------------------------------------
  ! 3. Initialize Arrays
  !---------------------------------------------------------------------------
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
        
  !---------------------------------------------------------------------------
  ! 4. Calculate Material Constants using E0_AGED
  !---------------------------------------------------------------------------
        ! Calculate elastic constants
        MU0 = E0_AGED/(2.0D0*(1.0D0 + NU))  ! Using E0_AGED
        K0 = E0_AGED/(3.0D0*(1.0D0 - 2.0D0*NU))  ! Using E0_AGED
        
        ! Read and modify Prony series parameters
        DO I = 1, NUM_PRONY_TERMS
            G_i(I) = PROPS(3+2*I-1)    ! Relative moduli
            TAU_i(I) = PROPS(3+2*I) * (1.0D0/AGING_INDEX)  ! Aged relaxation times
        END DO
        
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
        
        TOTAL_M = 1.D0
        DO N = 1, NUM_PRONY_TERMS
            TOTAL_M = TOTAL_M + M_i(N)
        END DO
        
  !---------------------------------------------------------------------------
  ! 5. Process State Variables
  !---------------------------------------------------------------------------
        DO N = 1, NUM_PRONY_TERMS
            DO I = 1, 6
                SM_OLD(N,I) = STATEV((N-1)*6 + I)
            END DO
        END DO
        
  !---------------------------------------------------------------------------
  ! 6. Calculate Strains
  !---------------------------------------------------------------------------
        VOL_STRAIN = (STRAN(1) + STRAN(2) + STRAN(3))/3.D0
        VOL_DSTRAIN = (DSTRAN(1) + DSTRAN(2) + DSTRAN(3))/3.D0
        
        DO I = 1, 3
            DEV_STRAIN(I) = STRAN(I) - VOL_STRAIN
            DEV_DSTRAIN(I) = DSTRAN(I) - VOL_DSTRAIN
        END DO
        DO I = 4, 6
            DEV_STRAIN(I) = STRAN(I)
            DEV_DSTRAIN(I) = DSTRAN(I)
        END DO
        
  !---------------------------------------------------------------------------
  ! 7. Calculate Stresses
  !---------------------------------------------------------------------------
        ! Hydrostatic stress
        HYDRO_PRESSURE = 3.D0 * K0 * (VOL_STRAIN + VOL_DSTRAIN)
        DO I = 1, 3
            STRESS_HYDROSTATIC(I) = HYDRO_PRESSURE
        END DO
        
        ! Long-term elastic stress
        DO I = 1, NTENS
            STRESS_LONGTERM(I) = 2.D0 * MU0 * G_INF * DEV_STRAIN(I)
        END DO
        
        ! History stress
        DO I = 1, NTENS
            STRESS_HISTORY(I) = 0.D0
            DO N = 1, NUM_PRONY_TERMS
                STRESS_HISTORY(I) = STRESS_HISTORY(I) + SM_OLD(N,I)
            END DO
        END DO
        
        ! Instantaneous stress
        DO I = 1, NTENS
            STRESS_INSTANTANEOUS(I) = 2.D0 * MU0 * TOTAL_M * DEV_DSTRAIN(I)
        END DO
        
        ! Combine all stress components
        DO I = 1, 3
            STRESS(I) = STRESS_HYDROSTATIC(I) + 
    1                   STRESS_LONGTERM(I) +
    2                   STRESS_HISTORY(I) +
    3                   STRESS_INSTANTANEOUS(I)
        END DO
        
        DO I = 4, 6
            STRESS(I) = STRESS_LONGTERM(I) +
    1                   STRESS_HISTORY(I) +
    2                   STRESS_INSTANTANEOUS(I)
        END DO
        
  !---------------------------------------------------------------------------
  ! 8. Update Internal Variables
  !---------------------------------------------------------------------------
        DO N = 1, NUM_PRONY_TERMS
            DO I = 1, 6
                SM(N,I) = SM_OLD(N,I) + 2.D0*MU0*M_i(N)*DEV_DSTRAIN(I)
                SM_DOT(N,I) = EXP(-DTIME/TAU_i(N))*SM(N,I)
                STATEV((N-1)*6 + I) = SM_DOT(N,I)
            END DO
        END DO
        
  !---------------------------------------------------------------------------
  ! 9. Calculate Material Jacobian
  !---------------------------------------------------------------------------
        DO I = 1, NTENS
            DO J = 1, NTENS
                DDSDDE(I,J) = 0.D0
            END DO
        END DO
        
        DO I = 1, 3
            DO J = 1, 3
                DDSDDE(I,J) = K0
            END DO
        END DO
        
        DO I = 1, 3
            DDSDDE(I,I) = DDSDDE(I,I) + 2.D0*MU0*TOTAL_M
        END DO
        
        DO I = 4, 6
            DDSDDE(I,I) = 2.D0*MU0*TOTAL_M
        END DO
        
  !---------------------------------------------------------------------------
  ! 10. Store Aging Variables
  !---------------------------------------------------------------------------
        STATEV(NSTATEV-3) = CA_t
        STATEV(NSTATEV-2) = AGING_INDEX
        STATEV(NSTATEV-1) = E0_AGED
        
  !---------------------------------------------------------------------------
  ! 11. Cleanup
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