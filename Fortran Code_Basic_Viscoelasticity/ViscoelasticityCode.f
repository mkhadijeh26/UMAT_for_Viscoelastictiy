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

      ! Local variables
      REAL*8 E0, NU, G0, K     !! E0 is the Young's modulus, NU is the Poisson's ratio, G0 is the shear modulus, and K is the bulk modulus
      REAL*8 GI(10), TAUI(10)  !! GI and TAUI are arrays to store the Prony series coefficients
      INTEGER NTERMS, I, J     !! NTERMS is the number of terms in the Prony series
      REAL*8 DVOLSTR, DDEVSTR(6), VOLSTR, DEVSTR(6)  !! DVOLSTR is the volumetric strain increment, DDEVSTR is the deviatoric strain increment
      REAL*8 QOLD(60), QNEW(60)  !! QOLD and QNEW are arrays to store the internal variables
      REAL*8 DSTRESS(6), STRESS_VE(6) !! DSTRESS is the stress increment, STRESS_VE is the viscoelastic stress

      ! Material properties
      E0 = PROPS(1)
      NU = PROPS(2)
      NTERMS = NINT(PROPS(3))  
      ! NINT is a Fortran intrinsic function that stands for "Nearest INTeger." It's used to round a real number to the nearest whole number (integer).
      DO I = 1, NTERMS   !g_1 , Tau_1 , g_2, Tau_2 , g_3, Tau_3 .... etc
        GI(I) = PROPS(3+2*I-1)
        TAUI(I) = PROPS(3+2*I)
      END DO

      ! Calculate initial shear and bulk moduli
      G0 = E0 / (2.0D0 * (1.0D0 + NU))
      K = E0 / (3.0D0 * (1.0D0 - 2.0D0 * NU))

      !! Volumetric Strain (DVOLSTR) represents the change in volume of a material element. 
      !! -- calculated by summing the normal components of the strain tensor
      !! Deviatoric Strain (DDEVSTR) represents the change in shape of a material element.
      !! It is calculated by subtracting the volumetric strain 
      !!  -- (divided equally across all directions) from the total strain in each direction.
      ! Calculate volumetric and deviatoric strain increments
      DVOLSTR = DSTRAN(1) + DSTRAN(2) + DSTRAN(3) !! Volumetric strain increment (trace of strain tensor) (sum of the normal strain components)
      DO I = 1, NDI
        DDEVSTR(I) = DSTRAN(I) - DVOLSTR/3.0D0    !! Deviatoric strain increment (deviatoric part of strain tensor) (difference between normal strain component and volumetric strain)
      END DO
      DO I = NDI+1, NTENS
        DDEVSTR(I) = DSTRAN(I)                    !! Deviatoric strain increment (deviatoric part of strain tensor)
      END DO

      ! Update internal variables (QOLD and QNEW)
      !! Q represents the portion of stress that relaxes over time.
      DO I = 1, NTERMS    !! 
        DO J = 1, NTENS   !! NTENS = 6 (3 normal strain components and 3 shear strain components)
          QOLD(6*(I-1)+J) = STATEV(6*(I-1)+J) !! QOLD is the value of the internal variable at the previous time step
          QNEW(6*(I-1)+J) = QOLD(6*(I-1)+J) * EXP(-DTIME/TAUI(I)) +  !! QNEW is the updated value of the internal variable
     1      2.0D0 * GI(I) * DDEVSTR(J) * (1.0D0 - EXP(-DTIME/TAUI(I)))
          STATEV(6*(I-1)+J) = QNEW(6*(I-1)+J) !! Update the internal variable in the state variable array
        END DO
      END DO

      ! Calculate stress increment
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
        STRESS_VE(I) = 0.0D0 !!The viscoelastic stress represents the part of the stress that decays over time due to the time-dependent nature of the material (as in the Maxwell model). Subtracting it from the total stress gives the final stress in the material at the current time step, which reflects both the elastic and time-dependent behaviors.
        DO J = 1, NTERMS
          STRESS_VE(I) = STRESS_VE(I) + QNEW(6*(J-1)+I) !! The viscoelastic stress is calculated by summing the contributions from each Prony series term
        END DO
        STRESS(I) = STRESS(I) - STRESS_VE(I)
      END DO

      ! Calculate material Jacobian
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
