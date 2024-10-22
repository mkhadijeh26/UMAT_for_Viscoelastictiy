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

     ! Local variables for viscoelastic behavior
     REAL*8 E0, E0_AGED, NU, G0, K
     REAL*8 GI(10), TAUI(10), TAUI_AGED(10), VISCOSITY(10)
     INTEGER NTERMS, I, J
     REAL*8 DVOLSTR, DDEVSTR(6), VOLSTR, DEVSTR(6)
     REAL*8 QOLD(60), QNEW(60)
     REAL*8 DSTRESS(6), STRESS_VE(6)

     ! Variables for aging
     REAL*8 K_AGING           ! Aging scaling factor
     REAL*8 CA_0             ! Initial carbonyl area
     REAL*8 CA_INF           ! Long-term carbonyl area
     REAL*8 R_REACTION       ! Reaction rate
     REAL*8 CA_CURRENT       ! Current carbonyl area
     REAL*8 AGING_INDEX      ! Current aging index

     ! Material properties (viscoelastic)
     E0 = PROPS(1)
     NU = PROPS(2)
     NTERMS = NINT(PROPS(3))
     
     ! Store initial Prony series parameters
     DO I = 1, NTERMS
       GI(I) = PROPS(3+2*I-1)
       TAUI(I) = PROPS(3+2*I)
       ! Calculate and store initial viscosity for each term
       VISCOSITY(I) = TAUI(I) * E0
     END DO

     ! Aging parameters
     K_AGING = PROPS(3+2*NTERMS+1)
     CA_0 = PROPS(3+2*NTERMS+2)
     CA_INF = PROPS(3+2*NTERMS+3)
     R_REACTION = PROPS(3+2*NTERMS+4)

     ! Calculate current carbonyl area (CA)
     CA_CURRENT = CA_0 + (CA_INF - CA_0) * 
    1 (1.0D0 - EXP(-R_REACTION * TIME(1)))

     ! Calculate aging index A(t)
     AGING_INDEX = 1.0D0 + K_AGING * CA_CURRENT

     ! Calculate aged elastic modulus
     E0_AGED = E0 * AGING_INDEX

     ! Calculate new relaxation times based on aged E0
     DO I = 1, NTERMS
       ! τᵢ_aged = η/E0_aged
       TAUI_AGED(I) = VISCOSITY(I) / E0_AGED
     END DO

     ! Calculate aged initial moduli
     G0 = E0_AGED / (2.0D0 * (1.0D0 + NU))
     K = E0_AGED / (3.0D0 * (1.0D0 - 2.0D0 * NU))

     ! Calculate volumetric and deviatoric strain increments
     DVOLSTR = DSTRAN(1) + DSTRAN(2) + DSTRAN(3)
     DO I = 1, NDI
       DDEVSTR(I) = DSTRAN(I) - DVOLSTR/3.0D0
     END DO
     DO I = NDI+1, NTENS
       DDEVSTR(I) = DSTRAN(I)
     END DO

     ! Update internal variables using aged relaxation times
     DO I = 1, NTERMS
       DO J = 1, NTENS
         QOLD(6*(I-1)+J) = STATEV(6*(I-1)+J)
         QNEW(6*(I-1)+J) = QOLD(6*(I-1)+J) * 
    1      EXP(-DTIME/TAUI_AGED(I)) +
    2      2.0D0 * GI(I) * DDEVSTR(J) * 
    3      (1.0D0 - EXP(-DTIME/TAUI_AGED(I)))
         STATEV(6*(I-1)+J) = QNEW(6*(I-1)+J)
       END DO
     END DO

     ! Store aging state variables
     STATEV(6*NTERMS+1) = CA_CURRENT
     STATEV(6*NTERMS+2) = AGING_INDEX
     STATEV(6*NTERMS+3) = E0_AGED

     ! Calculate stress increment using aged properties
     DO I = 1, NDI
       DSTRESS(I) = K * DVOLSTR + 2.0D0 * G0 * DDEVSTR(I)
     END DO
     DO I = NDI+1, NTENS
       DSTRESS(I) = 2.0D0 * G0 * DDEVSTR(I)
     END DO

     ! Update total stress
     DO I = 1, NTENS
       STRESS(I) = STRESS(I) + DSTRESS(I)
     END DO

     ! Calculate viscoelastic stress
     DO I = 1, NTENS
       STRESS_VE(I) = 0.0D0
       DO J = 1, NTERMS
         STRESS_VE(I) = STRESS_VE(I) + QNEW(6*(J-1)+I)
       END DO
       STRESS(I) = STRESS(I) - STRESS_VE(I)
     END DO

     ! Calculate material Jacobian with aged properties
     DO I = 1, NTENS
       DO J = 1, NTENS
         DDSDDE(I,J) = 0.0D0
       END DO
     END DO

     DO I = 1, NDI
       DDSDDE(I,I) = K + 4.0D0*G0/3.0D0
       DO J = 1, NDI
         IF (I .NE. J) THEN
           DDSDDE(I,J) = K - 2.0D0*G0/3.0D0
         END IF
       END DO
     END DO

     DO I = NDI+1, NTENS
       DDSDDE(I,I) = G0
     END DO

     RETURN
     END