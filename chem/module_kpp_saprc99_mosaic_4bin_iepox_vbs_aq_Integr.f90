
































MODULE saprc99_mosaic_4bin_iepox_vbs_aq_Integrator

 USE saprc99_mosaic_4bin_iepox_vbs_aq_Parameters
 USE saprc99_mosaic_4bin_iepox_vbs_aq_Precision
 USE saprc99_mosaic_4bin_iepox_vbs_aq_JacobianSP

  IMPLICIT NONE
 


  INTEGER, PARAMETER :: ifun=1, ijac=2, istp=3, iacc=4, &
    irej=5, idec=6, isol=7, isng=8, itexit=1, ihexit=2
    

  
  CHARACTER(LEN=50), PARAMETER, DIMENSION(-8:1) :: IERR_NAMES = (/ &
    'Matrix is repeatedly singular                     ', & 
    'Step size too small                               ', & 
    'No of steps exceeds maximum bound                 ', & 
    'Improper tolerance values                         ', & 
    'FacMin/FacMax/FacRej must be positive             ', & 
    'Hmin/Hmax/Hstart must be positive                 ', & 
    'Selected Rosenbrock method not implemented        ', & 
    'Improper value for maximal no of steps            ', & 
    '                                                  ', & 
    'Success                                           ' /) 

CONTAINS

SUBROUTINE  saprc99_mosaic_4bin_iepox_vbs_aq_INTEGRATE( TIN, TOUT, &
  FIX, VAR,  RCONST, ATOL, RTOL, IRR_WRK,  &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U  )

   USE saprc99_mosaic_4bin_iepox_vbs_aq_Parameters

   IMPLICIT NONE
   REAL(kind=dp), INTENT(INOUT), DIMENSION(NFIX) :: FIX
   REAL(kind=dp), INTENT(INOUT), DIMENSION(NVAR) :: VAR
   REAL(kind=dp), INTENT(INOUT) :: IRR_WRK(NREACT)
   REAL(kind=dp), INTENT(IN), DIMENSION(NSPEC) :: ATOL, RTOL
   REAL(kind=dp), INTENT(IN), DIMENSION(NREACT) :: RCONST
   REAL(kind=dp), INTENT(IN) :: TIN  
   REAL(kind=dp), INTENT(IN) :: TOUT 
   
   INTEGER,  INTENT(IN),  OPTIONAL :: ICNTRL_U(20)
   REAL(kind=dp), INTENT(IN),  OPTIONAL :: RCNTRL_U(20)
   INTEGER,  INTENT(OUT), OPTIONAL :: ISTATUS_U(20)
   REAL(kind=dp), INTENT(OUT), OPTIONAL :: RSTATUS_U(20)
   INTEGER,  INTENT(OUT), OPTIONAL :: IERR_U

   REAL(kind=dp)  :: STEPMIN


   INTEGER :: N_stp, N_acc, N_rej, N_sng
   SAVE N_stp, N_acc, N_rej, N_sng
   INTEGER :: i, IERR
   REAL(kind=dp) :: RCNTRL(20), RSTATUS(20)
   INTEGER :: ICNTRL(20), ISTATUS(20)


   ICNTRL(:)  = 0
   RCNTRL(:)  = 0.0_dp
   ISTATUS(:) = 0
   RSTATUS(:) = 0.0_dp

   
   
   IF (PRESENT(ICNTRL_U)) THEN
     WHERE(ICNTRL_U(:) > 0) ICNTRL(:) = ICNTRL_U(:)
   END IF
   IF (PRESENT(RCNTRL_U)) THEN
     WHERE(RCNTRL_U(:) > 0) RCNTRL(:) = RCNTRL_U(:)
   END IF

   CALL saprc99_mosaic_4bin_iepox_vbs_aq_Rosenbrock(VAR, FIX, RCONST, TIN,TOUT,   &
         ATOL,RTOL,               &
         RCNTRL,ICNTRL,RSTATUS,ISTATUS,IRR_WRK,IERR)

   STEPMIN = RCNTRL(ihexit)
   
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U))    IERR_U       = IERR

END SUBROUTINE  saprc99_mosaic_4bin_iepox_vbs_aq_INTEGRATE


SUBROUTINE  saprc99_mosaic_4bin_iepox_vbs_aq_Rosenbrock(Y, FIX, RCONST, Tstart,Tend, &
           AbsTol,RelTol,            &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IRR_WRK,IERR)







































































































  USE saprc99_mosaic_4bin_iepox_vbs_aq_Parameters

  IMPLICIT NONE


   REAL(kind=dp), INTENT(INOUT) :: Y(NVAR)
   REAL(kind=dp), INTENT(INOUT) :: IRR_WRK(NREACT)
   REAL(kind=dp), INTENT(IN), DIMENSION(NFIX) :: FIX
   REAL(kind=dp), INTENT(IN), DIMENSION(NREACT) :: RCONST
   REAL(kind=dp), INTENT(IN)   :: Tstart,Tend
   REAL(kind=dp), INTENT(IN)   :: AbsTol(NVAR),RelTol(NVAR)
   INTEGER, INTENT(IN)    :: ICNTRL(20)
   REAL(kind=dp), INTENT(IN)   :: RCNTRL(20)
   INTEGER, INTENT(INOUT) :: ISTATUS(20)
   REAL(kind=dp), INTENT(INOUT) :: RSTATUS(20)
   INTEGER, INTENT(OUT)   :: IERR

   INTEGER, PARAMETER :: Smax = 6
   INTEGER  :: Method, ros_S
   REAL(kind=dp), DIMENSION(Smax) :: ros_M, ros_E, ros_Alpha, ros_Gamma
   REAL(kind=dp), DIMENSION(Smax*(Smax-1)/2) :: ros_A, ros_C
   REAL(kind=dp) :: ros_ELO
   LOGICAL, DIMENSION(Smax) :: ros_NewF
   CHARACTER(LEN=12) :: ros_Name


  INTEGER :: Nfun,Njac,Nstp,Nacc,Nrej,Ndec,Nsol,Nsng



   REAL(kind=dp) :: Roundoff, FacMin, FacMax, FacRej, FacSafe
   REAL(kind=dp) :: Hmin, Hmax, Hstart, Hexit
   REAL(kind=dp) :: Texit
   INTEGER :: i, UplimTol, Max_no_steps
   LOGICAL :: Autonomous, VectorTol

   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp, ONE  = 1.0_dp
   REAL(kind=dp), PARAMETER :: DeltaMin = 1.0E-5_dp


   Nfun = ISTATUS(ifun)
   Njac = ISTATUS(ijac)
   Nstp = ISTATUS(istp)
   Nacc = ISTATUS(iacc)
   Nrej = ISTATUS(irej)
   Ndec = ISTATUS(idec)
   Nsol = ISTATUS(isol)
   Nsng = ISTATUS(isng)


   Autonomous = .NOT.(ICNTRL(1) == 0)







   IF (ICNTRL(2) == 0) THEN
      VectorTol = .TRUE.
         UplimTol  = NVAR
   ELSE
      VectorTol = .FALSE.
         UplimTol  = 1
   END IF


   IF (ICNTRL(3) == 0) THEN
      Method = 4
   ELSEIF ( (ICNTRL(3) >= 1).AND.(ICNTRL(3) <= 5) ) THEN
      Method = ICNTRL(3)
   ELSE
      PRINT * , 'User-selected Rosenbrock method: ICNTRL(3)=', Method
      CALL saprc99_mosaic_4bin_iepox_vbs_aq_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
      RETURN
   END IF


   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 100000
   ELSEIF (ICNTRL(4) > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL saprc99_mosaic_4bin_iepox_vbs_aq_ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF


   Roundoff = saprc99_mosaic_4bin_iepox_vbs_aq_WLAMCH('E')


   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      Hmin = RCNTRL(1)
   ELSE
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL saprc99_mosaic_4bin_iepox_vbs_aq_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL saprc99_mosaic_4bin_iepox_vbs_aq_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL saprc99_mosaic_4bin_iepox_vbs_aq_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2_dp
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL saprc99_mosaic_4bin_iepox_vbs_aq_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0_dp
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL saprc99_mosaic_4bin_iepox_vbs_aq_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1_dp
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL saprc99_mosaic_4bin_iepox_vbs_aq_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9_dp
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL saprc99_mosaic_4bin_iepox_vbs_aq_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.0_dp*Roundoff) &
         .OR. (RelTol(i) >= 1.0_dp) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL saprc99_mosaic_4bin_iepox_vbs_aq_ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO



   SELECT CASE (Method)
     CASE (1)
       CALL saprc99_mosaic_4bin_iepox_vbs_aq_Ros2(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (2)
       CALL saprc99_mosaic_4bin_iepox_vbs_aq_Ros3(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (3)
       CALL saprc99_mosaic_4bin_iepox_vbs_aq_Ros4(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (4)
       CALL saprc99_mosaic_4bin_iepox_vbs_aq_Rodas3(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (5)
       CALL saprc99_mosaic_4bin_iepox_vbs_aq_Rodas4(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(4)=', Method
       CALL saprc99_mosaic_4bin_iepox_vbs_aq_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       RETURN
   END SELECT


   CALL saprc99_mosaic_4bin_iepox_vbs_aq_ros_Integrator(Y,Tstart,Tend,Texit,      &
        AbsTol, RelTol,                          &

        ros_S, ros_M, ros_E, ros_A, ros_C,       &
        ros_Alpha, ros_Gamma, ros_ELO, ros_NewF, &

        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart, Hexit,     &
        FacMin, FacMax, FacRej, FacSafe,         &

        IRR_WRK,IERR,                            &

         Nfun,Njac,Nstp,Nacc,Nrej,Ndec,Nsol,Nsng,&

         RCONST, FIX &
)



   ISTATUS(ifun) = Nfun
   ISTATUS(ijac) = Njac
   ISTATUS(istp) = Nstp
   ISTATUS(iacc) = Nacc
   ISTATUS(irej) = Nrej
   ISTATUS(idec) = Ndec
   ISTATUS(isol) = Nsol
   ISTATUS(isng) = Nsng

   RSTATUS(itexit) = Texit
   RSTATUS(ihexit) = Hexit


CONTAINS 



 SUBROUTINE  saprc99_mosaic_4bin_iepox_vbs_aq_ros_ErrorMsg(Code,T,H,IERR)



   USE saprc99_mosaic_4bin_iepox_vbs_aq_Precision

   REAL(kind=dp), INTENT(IN) :: T, H
   INTEGER, INTENT(IN)  :: Code
   INTEGER, INTENT(OUT) :: IERR

   IERR = Code
   PRINT * , &
     'Forced exit from Rosenbrock due to the following error:'
   IF ((Code>=-8).AND.(Code<=-1)) THEN
     PRINT *, IERR_NAMES(Code)
   ELSE
     PRINT *, 'Unknown Error code: ', Code
   ENDIF

   PRINT *, "T=", T, "and H=", H

 END SUBROUTINE  saprc99_mosaic_4bin_iepox_vbs_aq_ros_ErrorMsg


 SUBROUTINE  saprc99_mosaic_4bin_iepox_vbs_aq_ros_Integrator (Y, Tstart, Tend, T,     &
        AbsTol, RelTol,                          &

        ros_S, ros_M, ros_E, ros_A, ros_C,       &
        ros_Alpha, ros_Gamma, ros_ELO, ros_NewF, &

        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart, Hexit,     &
        FacMin, FacMax, FacRej, FacSafe,         &

        IRR_WRK,IERR,                            &

        Nfun,Njac,Nstp,Nacc,Nrej,Ndec,Nsol,Nsng, &

        RCONST, FIX )






  IMPLICIT NONE


   REAL(kind=dp), INTENT(INOUT) :: Y(NVAR)

   REAL(kind=dp), INTENT(INOUT) :: IRR_WRK(NREACT)

   REAL(kind=dp), INTENT(IN) :: Tstart,Tend

   REAL(kind=dp), INTENT(OUT) ::  T

   REAL(kind=dp), INTENT(IN) ::  AbsTol(NVAR), RelTol(NVAR)

   INTEGER, INTENT(IN) ::  ros_S
   REAL(kind=dp), INTENT(IN) :: ros_M(ros_S), ros_E(ros_S),  &
       ros_Alpha(ros_S), ros_A(ros_S*(ros_S-1)/2), &
       ros_Gamma(ros_S), ros_C(ros_S*(ros_S-1)/2), ros_ELO
   LOGICAL, INTENT(IN) :: ros_NewF(ros_S)

   LOGICAL, INTENT(IN) :: Autonomous, VectorTol
   REAL(kind=dp), INTENT(IN) :: Hstart, Hmin, Hmax
   INTEGER, INTENT(IN) :: Max_no_steps
   REAL(kind=dp), INTENT(IN) :: Roundoff, FacMin, FacMax, FacRej, FacSafe

   REAL(kind=dp), INTENT(OUT) :: Hexit

   INTEGER, INTENT(OUT) :: IERR

   REAL(kind=dp), INTENT(IN), DIMENSION(NFIX) :: FIX

   REAL(kind=dp), INTENT(IN), DIMENSION(NREACT) :: RCONST


  INTEGER, INTENT(INOUT)  :: Nfun,Njac,Nstp,Nacc,Nrej,Ndec,Nsol,Nsng


   REAL(kind=dp) :: Ynew(NVAR), Fcn0(NVAR), Fcn(NVAR)
   REAL(kind=dp) :: K(NVAR*ros_S), dFdT(NVAR)
   REAL(kind=dp) :: Jac0(LU_NONZERO), Ghimj(LU_NONZERO)
   REAL(kind=dp) :: H, Hnew, HC, HG, Fac, Tau
   REAL(kind=dp) :: Err, Yerr(NVAR)
   INTEGER :: Pivot(NVAR), Direction, ioffset, j, istage
   LOGICAL :: RejectLastH, RejectMoreH, Singular

   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp, ONE  = 1.0_dp
   REAL(kind=dp), PARAMETER :: DeltaMin = 1.0E-5_dp







   T = Tstart
   Hexit = 0.0_dp
   H = MIN(Hstart,Hmax)
   IF (ABS(H) <= 10.0_dp*Roundoff) H = DeltaMin

   IF (Tend  >=  Tstart) THEN
     Direction = +1
   ELSE
     Direction = -1
   END IF

   RejectLastH=.FALSE.
   RejectMoreH=.FALSE.



TimeLoop: DO WHILE ( (Direction > 0).AND.((T-Tend)+Roundoff <= ZERO) &
       .OR. (Direction < 0).AND.((Tend-T)+Roundoff <= ZERO) )

   IF ( Nstp > Max_no_steps ) THEN  
      CALL saprc99_mosaic_4bin_iepox_vbs_aq_ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1_dp*H) == T).OR.(H <= Roundoff) ) THEN  
      CALL saprc99_mosaic_4bin_iepox_vbs_aq_ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF


   Hexit = H
   H = MIN(H,ABS(Tend-T))


   CALL saprc99_mosaic_4bin_iepox_vbs_aq_FunTemplate(T,Y,Fcn0, RCONST, FIX, Nfun)
   IF( T == Tstart ) THEN
     CALL saprc99_mosaic_4bin_iepox_vbs_aq_IRRFun( Y, FIX, RCONST, IRR_WRK )
   ENDIF


   IF (.NOT.Autonomous) THEN
      CALL saprc99_mosaic_4bin_iepox_vbs_aq_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )
   END IF


   CALL saprc99_mosaic_4bin_iepox_vbs_aq_JacTemplate(T,Y,Jac0, FIX, Njac, RCONST)


UntilAccepted: DO

   CALL saprc99_mosaic_4bin_iepox_vbs_aq_ros_PrepareMatrix(H,Direction,ros_Gamma(1), &
          Jac0,Ghimj,Pivot,Singular, Ndec,  Nsng )
   IF (Singular) THEN 
       CALL saprc99_mosaic_4bin_iepox_vbs_aq_ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF


Stage: DO istage = 1, ros_S

      
       ioffset = NVAR*(istage-1)

      
       IF ( istage == 1 ) THEN
         CALL saprc99_mosaic_4bin_iepox_vbs_aq_WCOPY(NVAR,Fcn0,1,Fcn,1)
      
       ELSEIF ( ros_NewF(istage) ) THEN
         CALL saprc99_mosaic_4bin_iepox_vbs_aq_WCOPY(NVAR,Y,1,Ynew,1)
         DO j = 1, istage-1
           CALL saprc99_mosaic_4bin_iepox_vbs_aq_WAXPY(NVAR,ros_A((istage-1)*(istage-2)/2+j), &
            K(NVAR*(j-1)+1),1,Ynew,1)
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL saprc99_mosaic_4bin_iepox_vbs_aq_FunTemplate(Tau,Ynew,Fcn, RCONST, FIX, Nfun)
       END IF 
       CALL saprc99_mosaic_4bin_iepox_vbs_aq_WCOPY(NVAR,Fcn,1,K(ioffset+1),1)
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         CALL saprc99_mosaic_4bin_iepox_vbs_aq_WAXPY(NVAR,HC,K(NVAR*(j-1)+1),1,K(ioffset+1),1)
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         CALL saprc99_mosaic_4bin_iepox_vbs_aq_WAXPY(NVAR,HG,dFdT,1,K(ioffset+1),1)
       END IF
       CALL saprc99_mosaic_4bin_iepox_vbs_aq_ros_Solve(Ghimj, Pivot, K(ioffset+1), Nsol)

   END DO Stage



   CALL saprc99_mosaic_4bin_iepox_vbs_aq_WCOPY(NVAR,Y,1,Ynew,1)
   DO j=1,ros_S
         CALL saprc99_mosaic_4bin_iepox_vbs_aq_WAXPY(NVAR,ros_M(j),K(NVAR*(j-1)+1),1,Ynew,1)
   END DO


   CALL saprc99_mosaic_4bin_iepox_vbs_aq_WSCAL(NVAR,ZERO,Yerr,1)
   DO j=1,ros_S
        CALL saprc99_mosaic_4bin_iepox_vbs_aq_WAXPY(NVAR,ros_E(j),K(NVAR*(j-1)+1),1,Yerr,1)
   END DO
   Err = saprc99_mosaic_4bin_iepox_vbs_aq_ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )


   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac


   Nstp = Nstp+1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  
      Nacc = Nacc+1
      CALL saprc99_mosaic_4bin_iepox_vbs_aq_WCOPY(NVAR,Ynew,1,Y,1)
      T = T + Direction*H
      Hnew = MAX(Hmin,MIN(Hnew,Hmax))
      IF (RejectLastH) THEN  
         Hnew = MIN(Hnew,H)
      END IF
      RejectLastH = .FALSE.
      RejectMoreH = .FALSE.
      H = Hnew
      EXIT UntilAccepted 
   ELSE           
      IF (RejectMoreH) THEN
         Hnew = H*FacRej
      END IF
      RejectMoreH = RejectLastH
      RejectLastH = .TRUE.
      H = Hnew
      IF (Nacc >= 1) THEN
         Nrej = Nrej+1
      END IF
   END IF 

   END DO UntilAccepted

   END DO TimeLoop


   IERR = 1  

  END SUBROUTINE  saprc99_mosaic_4bin_iepox_vbs_aq_ros_Integrator



  REAL(kind=dp) FUNCTION  saprc99_mosaic_4bin_iepox_vbs_aq_ros_ErrorNorm ( Y, Ynew, Yerr, &
               AbsTol, RelTol, VectorTol )



   IMPLICIT NONE


   REAL(kind=dp), INTENT(IN) :: Y(NVAR), Ynew(NVAR), &
          Yerr(NVAR), AbsTol(NVAR), RelTol(NVAR)
   LOGICAL, INTENT(IN) ::  VectorTol

   REAL(kind=dp) :: Err, Scale, Ymax
   INTEGER  :: i
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp

   Err = ZERO
   DO i=1,NVAR
     Ymax = MAX(ABS(Y(i)),ABS(Ynew(i)))
     IF (VectorTol) THEN
       Scale = AbsTol(i)+RelTol(i)*Ymax
     ELSE
       Scale = AbsTol(1)+RelTol(1)*Ymax
     END IF
     Err = Err+(Yerr(i)/Scale)**2
   END DO
   Err  = SQRT(Err/NVAR)

    saprc99_mosaic_4bin_iepox_vbs_aq_ros_ErrorNorm = Err

  END FUNCTION  saprc99_mosaic_4bin_iepox_vbs_aq_ros_ErrorNorm



  SUBROUTINE saprc99_mosaic_4bin_iepox_vbs_aq_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )



   IMPLICIT NONE


   REAL(kind=dp), INTENT(IN) :: T, Roundoff, Y(NVAR), Fcn0(NVAR)
   REAL(kind=dp), INTENT(IN) :: RCONST(NREACT), FIX(NFIX)

   REAL(kind=dp), INTENT(OUT) :: dFdT(NVAR)

   INTEGER, INTENT(INOUT) ::Nfun

   REAL(kind=dp) :: Delta
   REAL(kind=dp), PARAMETER :: ONE = 1.0_dp, DeltaMin = 1.0E-6_dp

   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL saprc99_mosaic_4bin_iepox_vbs_aq_FunTemplate(T+Delta,Y,dFdT, RCONST, FIX, Nfun)
   CALL saprc99_mosaic_4bin_iepox_vbs_aq_WAXPY(NVAR,(-ONE),Fcn0,1,dFdT,1)
   CALL saprc99_mosaic_4bin_iepox_vbs_aq_WSCAL(NVAR,(ONE/Delta),dFdT,1)

  END SUBROUTINE  saprc99_mosaic_4bin_iepox_vbs_aq_ros_FunTimeDeriv



  SUBROUTINE  saprc99_mosaic_4bin_iepox_vbs_aq_ros_PrepareMatrix ( H, Direction, gam, &
             Jac0, Ghimj, Pivot, Singular, Ndec,  Nsng  )








   IMPLICIT NONE


   REAL(kind=dp), INTENT(IN) ::  Jac0(LU_NONZERO)
   REAL(kind=dp), INTENT(IN) ::  gam
   INTEGER, INTENT(IN) ::  Direction

   REAL(kind=dp), INTENT(OUT) :: Ghimj(LU_NONZERO)
   LOGICAL, INTENT(OUT) ::  Singular
   INTEGER, INTENT(OUT) ::  Pivot(NVAR)

   REAL(kind=dp), INTENT(INOUT) :: H   
   INTEGER, INTENT(INOUT) ::  Ndec, Nsng

   INTEGER  :: i, ising, Nconsecutive
   REAL(kind=dp) :: ghinv
   REAL(kind=dp), PARAMETER :: ONE  = 1.0_dp, HALF = 0.5_dp

   Nconsecutive = 0
   Singular = .TRUE.

   DO WHILE (Singular)


     CALL saprc99_mosaic_4bin_iepox_vbs_aq_WCOPY(LU_NONZERO,Jac0,1,Ghimj,1)
     CALL saprc99_mosaic_4bin_iepox_vbs_aq_WSCAL(LU_NONZERO,(-ONE),Ghimj,1)
     ghinv = ONE/(Direction*H*gam)
     DO i=1,NVAR
       Ghimj(LU_DIAG(i)) = Ghimj(LU_DIAG(i))+ghinv
     END DO

     CALL saprc99_mosaic_4bin_iepox_vbs_aq_ros_Decomp( Ghimj, Pivot, ising, Ndec )
     IF (ising == 0) THEN

        Singular = .FALSE.
     ELSE 

        Nsng = Nsng+1
        Nconsecutive = Nconsecutive+1
        Singular = .TRUE.
        PRINT*,'Warning: LU Decomposition returned ising = ',ising
        IF (Nconsecutive <= 5) THEN 
           H = H*HALF
        ELSE  
           RETURN
        END IF  
      END IF    

   END DO 

  END SUBROUTINE  saprc99_mosaic_4bin_iepox_vbs_aq_ros_PrepareMatrix



  SUBROUTINE  saprc99_mosaic_4bin_iepox_vbs_aq_ros_Decomp( A, Pivot, ising, Ndec )



   IMPLICIT NONE

   REAL(kind=dp), INTENT(INOUT) :: A(LU_NONZERO)

   INTEGER, INTENT(OUT) :: Pivot(NVAR), ising
   INTEGER, INTENT(INOUT) :: Ndec 

   CALL saprc99_mosaic_4bin_iepox_vbs_aq_KppDecomp ( A, ising )
   Pivot(1) = 1
   Ndec = Ndec + 1

  END SUBROUTINE  saprc99_mosaic_4bin_iepox_vbs_aq_ros_Decomp



  SUBROUTINE  saprc99_mosaic_4bin_iepox_vbs_aq_ros_Solve( A, Pivot, b, Nsol )



   IMPLICIT NONE

   REAL(kind=dp), INTENT(IN) :: A(LU_NONZERO)
   INTEGER, INTENT(IN) :: Pivot(NVAR)

   INTEGER, INTENT(INOUT) :: nsol 

   REAL(kind=dp), INTENT(INOUT) :: b(NVAR)


   CALL saprc99_mosaic_4bin_iepox_vbs_aq_KppSolve( A, b )

   Nsol = Nsol+1

  END SUBROUTINE  saprc99_mosaic_4bin_iepox_vbs_aq_ros_Solve




  SUBROUTINE  saprc99_mosaic_4bin_iepox_vbs_aq_Ros2 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
            ros_Gamma,ros_NewF,ros_ELO,ros_Name)




  IMPLICIT NONE

   INTEGER, PARAMETER :: S=2
   INTEGER, INTENT(OUT) ::  ros_S
   REAL(kind=dp), DIMENSION(S), INTENT(OUT) :: ros_M,ros_E,ros_Alpha,ros_Gamma
   REAL(kind=dp), DIMENSION(S*(S-1)/2), INTENT(OUT) :: ros_A, ros_C
   REAL(kind=dp), INTENT(OUT) :: ros_ELO
   LOGICAL, DIMENSION(S), INTENT(OUT) :: ros_NewF
   CHARACTER(LEN=12), INTENT(OUT) :: ros_Name

    REAL(kind=dp) :: g

    g = 1.0_dp + 1.0_dp/SQRT(2.0_dp)


    ros_Name = 'ROS-2'

    ros_S = S








    ros_A(1) = (1.0_dp)/g
    ros_C(1) = (-2.0_dp)/g


    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.

    ros_M(1)= (3.0_dp)/(2.0_dp*g)
    ros_M(2)= (1.0_dp)/(2.0_dp*g)

    ros_E(1) = 1.0_dp/(2.0_dp*g)
    ros_E(2) = 1.0_dp/(2.0_dp*g)


    ros_ELO = 2.0_dp

    ros_Alpha(1) = 0.0_dp
    ros_Alpha(2) = 1.0_dp

    ros_Gamma(1) = g
    ros_Gamma(2) =-g

 END SUBROUTINE  saprc99_mosaic_4bin_iepox_vbs_aq_Ros2



  SUBROUTINE  saprc99_mosaic_4bin_iepox_vbs_aq_Ros3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
           ros_Gamma,ros_NewF,ros_ELO,ros_Name)




  IMPLICIT NONE

   INTEGER, PARAMETER :: S=3
   INTEGER, INTENT(OUT) ::  ros_S
   REAL(kind=dp), DIMENSION(S), INTENT(OUT) :: ros_M,ros_E,ros_Alpha,ros_Gamma
   REAL(kind=dp), DIMENSION(S*(S-1)/2), INTENT(OUT) :: ros_A, ros_C
   REAL(kind=dp), INTENT(OUT) :: ros_ELO
   LOGICAL, DIMENSION(S), INTENT(OUT) :: ros_NewF
   CHARACTER(LEN=12), INTENT(OUT) :: ros_Name


   ros_Name = 'ROS-3'

   ros_S = S








   ros_A(1)= 1.0_dp
   ros_A(2)= 1.0_dp
   ros_A(3)= 0.0_dp

   ros_C(1) = -0.10156171083877702091975600115545E+01_dp
   ros_C(2) =  0.40759956452537699824805835358067E+01_dp
   ros_C(3) =  0.92076794298330791242156818474003E+01_dp


   ros_NewF(1) = .TRUE.
   ros_NewF(2) = .TRUE.
   ros_NewF(3) = .FALSE.

   ros_M(1) =  0.1E+01_dp
   ros_M(2) =  0.61697947043828245592553615689730E+01_dp
   ros_M(3) = -0.42772256543218573326238373806514E+00_dp

   ros_E(1) =  0.5E+00_dp
   ros_E(2) = -0.29079558716805469821718236208017E+01_dp
   ros_E(3) =  0.22354069897811569627360909276199E+00_dp


   ros_ELO = 3.0_dp

   ros_Alpha(1)= 0.0E+00_dp
   ros_Alpha(2)= 0.43586652150845899941601945119356E+00_dp
   ros_Alpha(3)= 0.43586652150845899941601945119356E+00_dp

   ros_Gamma(1)= 0.43586652150845899941601945119356E+00_dp
   ros_Gamma(2)= 0.24291996454816804366592249683314E+00_dp
   ros_Gamma(3)= 0.21851380027664058511513169485832E+01_dp

  END SUBROUTINE  saprc99_mosaic_4bin_iepox_vbs_aq_Ros3





  SUBROUTINE  saprc99_mosaic_4bin_iepox_vbs_aq_Ros4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
           ros_Gamma,ros_NewF,ros_ELO,ros_Name)










  IMPLICIT NONE

   INTEGER, PARAMETER :: S=4
   INTEGER, INTENT(OUT) ::  ros_S
   REAL(kind=dp), DIMENSION(4), INTENT(OUT) :: ros_M,ros_E,ros_Alpha,ros_Gamma
   REAL(kind=dp), DIMENSION(6), INTENT(OUT) :: ros_A, ros_C
   REAL(kind=dp), INTENT(OUT) :: ros_ELO
   LOGICAL, DIMENSION(4), INTENT(OUT) :: ros_NewF
   CHARACTER(LEN=12), INTENT(OUT) :: ros_Name

   REAL(kind=dp) :: g



   ros_Name = 'ROS-4'

   ros_S = S








   ros_A(1) = 0.2000000000000000E+01_dp
   ros_A(2) = 0.1867943637803922E+01_dp
   ros_A(3) = 0.2344449711399156E+00_dp
   ros_A(4) = ros_A(2)
   ros_A(5) = ros_A(3)
   ros_A(6) = 0.0_dp

   ros_C(1) =-0.7137615036412310E+01_dp
   ros_C(2) = 0.2580708087951457E+01_dp
   ros_C(3) = 0.6515950076447975E+00_dp
   ros_C(4) =-0.2137148994382534E+01_dp
   ros_C(5) =-0.3214669691237626E+00_dp
   ros_C(6) =-0.6949742501781779E+00_dp


   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .TRUE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .FALSE.

   ros_M(1) = 0.2255570073418735E+01_dp
   ros_M(2) = 0.2870493262186792E+00_dp
   ros_M(3) = 0.4353179431840180E+00_dp
   ros_M(4) = 0.1093502252409163E+01_dp

   ros_E(1) =-0.2815431932141155E+00_dp
   ros_E(2) =-0.7276199124938920E-01_dp
   ros_E(3) =-0.1082196201495311E+00_dp
   ros_E(4) =-0.1093502252409163E+01_dp


   ros_ELO  = 4.0_dp

   ros_Alpha(1) = 0.0_dp
   ros_Alpha(2) = 0.1145640000000000E+01_dp
   ros_Alpha(3) = 0.6552168638155900E+00_dp
   ros_Alpha(4) = ros_Alpha(3)

   ros_Gamma(1) = 0.5728200000000000E+00_dp
   ros_Gamma(2) =-0.1769193891319233E+01_dp
   ros_Gamma(3) = 0.7592633437920482E+00_dp
   ros_Gamma(4) =-0.1049021087100450E+00_dp

  END SUBROUTINE  saprc99_mosaic_4bin_iepox_vbs_aq_Ros4


  SUBROUTINE  saprc99_mosaic_4bin_iepox_vbs_aq_Rodas3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
            ros_Gamma,ros_NewF,ros_ELO,ros_Name)




  IMPLICIT NONE

   INTEGER, PARAMETER :: S=4
   INTEGER, INTENT(OUT) ::  ros_S
   REAL(kind=dp), DIMENSION(S), INTENT(OUT) :: ros_M,ros_E,ros_Alpha,ros_Gamma
   REAL(kind=dp), DIMENSION(S*(S-1)/2), INTENT(OUT) :: ros_A, ros_C
   REAL(kind=dp), INTENT(OUT) :: ros_ELO
   LOGICAL, DIMENSION(S), INTENT(OUT) :: ros_NewF
   CHARACTER(LEN=12), INTENT(OUT) :: ros_Name

   REAL(kind=dp) :: g


   ros_Name = 'RODAS-3'

   ros_S = S








   ros_A(1) = 0.0E+00_dp
   ros_A(2) = 2.0E+00_dp
   ros_A(3) = 0.0E+00_dp
   ros_A(4) = 2.0E+00_dp
   ros_A(5) = 0.0E+00_dp
   ros_A(6) = 1.0E+00_dp

   ros_C(1) = 4.0E+00_dp
   ros_C(2) = 1.0E+00_dp
   ros_C(3) =-1.0E+00_dp
   ros_C(4) = 1.0E+00_dp
   ros_C(5) =-1.0E+00_dp
   ros_C(6) =-(8.0E+00_dp/3.0E+00_dp)



   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .FALSE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .TRUE.

   ros_M(1) = 2.0E+00_dp
   ros_M(2) = 0.0E+00_dp
   ros_M(3) = 1.0E+00_dp
   ros_M(4) = 1.0E+00_dp

   ros_E(1) = 0.0E+00_dp
   ros_E(2) = 0.0E+00_dp
   ros_E(3) = 0.0E+00_dp
   ros_E(4) = 1.0E+00_dp


   ros_ELO  = 3.0E+00_dp

   ros_Alpha(1) = 0.0E+00_dp
   ros_Alpha(2) = 0.0E+00_dp
   ros_Alpha(3) = 1.0E+00_dp
   ros_Alpha(4) = 1.0E+00_dp

   ros_Gamma(1) = 0.5E+00_dp
   ros_Gamma(2) = 1.5E+00_dp
   ros_Gamma(3) = 0.0E+00_dp
   ros_Gamma(4) = 0.0E+00_dp

  END SUBROUTINE  saprc99_mosaic_4bin_iepox_vbs_aq_Rodas3


  SUBROUTINE  saprc99_mosaic_4bin_iepox_vbs_aq_Rodas4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
             ros_Gamma,ros_NewF,ros_ELO,ros_Name)









  IMPLICIT NONE

   INTEGER, PARAMETER :: S=6
   INTEGER, INTENT(OUT) ::  ros_S
   REAL(kind=dp), DIMENSION(S), INTENT(OUT) :: ros_M,ros_E,ros_Alpha,ros_Gamma
   REAL(kind=dp), DIMENSION(S*(S-1)/2), INTENT(OUT) :: ros_A, ros_C
   REAL(kind=dp), INTENT(OUT) :: ros_ELO
   LOGICAL, DIMENSION(S), INTENT(OUT) :: ros_NewF
   CHARACTER(LEN=12), INTENT(OUT) :: ros_Name

    REAL(kind=dp) :: g


    ros_Name = 'RODAS-4'

    ros_S = 6


    ros_Alpha(1) = 0.000_dp
    ros_Alpha(2) = 0.386_dp
    ros_Alpha(3) = 0.210_dp
    ros_Alpha(4) = 0.630_dp
    ros_Alpha(5) = 1.000_dp
    ros_Alpha(6) = 1.000_dp


    ros_Gamma(1) = 0.2500000000000000E+00_dp
    ros_Gamma(2) =-0.1043000000000000E+00_dp
    ros_Gamma(3) = 0.1035000000000000E+00_dp
    ros_Gamma(4) =-0.3620000000000023E-01_dp
    ros_Gamma(5) = 0.0_dp
    ros_Gamma(6) = 0.0_dp







    ros_A(1) = 0.1544000000000000E+01_dp
    ros_A(2) = 0.9466785280815826E+00_dp
    ros_A(3) = 0.2557011698983284E+00_dp
    ros_A(4) = 0.3314825187068521E+01_dp
    ros_A(5) = 0.2896124015972201E+01_dp
    ros_A(6) = 0.9986419139977817E+00_dp
    ros_A(7) = 0.1221224509226641E+01_dp
    ros_A(8) = 0.6019134481288629E+01_dp
    ros_A(9) = 0.1253708332932087E+02_dp
    ros_A(10) =-0.6878860361058950E+00_dp
    ros_A(11) = ros_A(7)
    ros_A(12) = ros_A(8)
    ros_A(13) = ros_A(9)
    ros_A(14) = ros_A(10)
    ros_A(15) = 1.0E+00_dp

    ros_C(1) =-0.5668800000000000E+01_dp
    ros_C(2) =-0.2430093356833875E+01_dp
    ros_C(3) =-0.2063599157091915E+00_dp
    ros_C(4) =-0.1073529058151375E+00_dp
    ros_C(5) =-0.9594562251023355E+01_dp
    ros_C(6) =-0.2047028614809616E+02_dp
    ros_C(7) = 0.7496443313967647E+01_dp
    ros_C(8) =-0.1024680431464352E+02_dp
    ros_C(9) =-0.3399990352819905E+02_dp
    ros_C(10) = 0.1170890893206160E+02_dp
    ros_C(11) = 0.8083246795921522E+01_dp
    ros_C(12) =-0.7981132988064893E+01_dp
    ros_C(13) =-0.3152159432874371E+02_dp
    ros_C(14) = 0.1631930543123136E+02_dp
    ros_C(15) =-0.6058818238834054E+01_dp


    ros_M(1) = ros_A(7)
    ros_M(2) = ros_A(8)
    ros_M(3) = ros_A(9)
    ros_M(4) = ros_A(10)
    ros_M(5) = 1.0E+00_dp
    ros_M(6) = 1.0E+00_dp


    ros_E(1) = 0.0E+00_dp
    ros_E(2) = 0.0E+00_dp
    ros_E(3) = 0.0E+00_dp
    ros_E(4) = 0.0E+00_dp
    ros_E(5) = 0.0E+00_dp
    ros_E(6) = 1.0E+00_dp



    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
    ros_NewF(3) = .TRUE.
    ros_NewF(4) = .TRUE.
    ros_NewF(5) = .TRUE.
    ros_NewF(6) = .TRUE.



    ros_ELO = 4.0_dp

  END SUBROUTINE  saprc99_mosaic_4bin_iepox_vbs_aq_Rodas4




END SUBROUTINE  saprc99_mosaic_4bin_iepox_vbs_aq_Rosenbrock




SUBROUTINE  saprc99_mosaic_4bin_iepox_vbs_aq_FunTemplate( T, Y, Ydot, RCONST, FIX, Nfun )




   USE saprc99_mosaic_4bin_iepox_vbs_aq_Parameters




   REAL(kind=dp) :: T, Y(NVAR)
   REAL(kind=dp) :: RCONST(NREACT)
   REAL(kind=dp) :: FIX(NFIX)

   REAL(kind=dp) :: Ydot(NVAR)
   INTEGER :: Nfun









   CALL saprc99_mosaic_4bin_iepox_vbs_aq_Fun( Y, FIX, RCONST, Ydot )


   Nfun = Nfun+1

END SUBROUTINE  saprc99_mosaic_4bin_iepox_vbs_aq_FunTemplate



SUBROUTINE  saprc99_mosaic_4bin_iepox_vbs_aq_JacTemplate( T, Y, Jcb, FIX, Njac, RCONST )




 USE saprc99_mosaic_4bin_iepox_vbs_aq_Parameters
 
 USE saprc99_mosaic_4bin_iepox_vbs_aq_Jacobian



    REAL(kind=dp) :: T, Y(NVAR)
    REAL(kind=dp) :: FIX(NFIX)
    REAL(kind=dp) :: RCONST(NREACT)

    INTEGER :: Njac


    REAL(kind=dp) :: Jcb(LU_NONZERO)

    REAL(kind=dp) :: Told





    CALL saprc99_mosaic_4bin_iepox_vbs_aq_Jac_SP( Y, FIX, RCONST, Jcb )


    Njac = Njac+1

END SUBROUTINE  saprc99_mosaic_4bin_iepox_vbs_aq_JacTemplate

















SUBROUTINE saprc99_mosaic_4bin_iepox_vbs_aq_Fun ( V, F, RCT, Vdot )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: Vdot(NVAR)




  REAL(kind=dp) :: A(NREACT)


  A(1) = RCT(1)*V(109)
  A(2) = RCT(2)*V(102)*F(2)
  A(3) = RCT(3)*V(102)*V(106)
  A(4) = RCT(4)*V(102)*V(115)*F(2)
  A(5) = RCT(5)*V(102)*V(109)
  A(6) = RCT(6)*V(102)*V(109)
  A(7) = RCT(7)*V(106)*V(115)
  A(8) = RCT(8)*V(106)*V(109)
  A(9) = RCT(9)*V(115)*V(118)
  A(10) = RCT(10)*V(115)*V(115)*F(2)
  A(11) = RCT(11)*V(109)*V(118)
  A(12) = RCT(12)*V(58)
  A(13) = RCT(13)*V(58)*F(1)
  A(14) = RCT(14)*V(109)*V(118)
  A(15) = RCT(15)*V(118)
  A(16) = RCT(16)*V(118)
  A(17) = RCT(17)*V(106)
  A(18) = RCT(18)*V(106)
  A(19) = RCT(19)*V(31)*F(1)
  A(20) = RCT(20)*V(31)*F(2)
  A(21) = RCT(21)*V(114)*V(115)
  A(22) = RCT(22)*V(60)
  A(23) = RCT(23)*V(60)
  A(24) = RCT(24)*V(60)*V(114)
  A(25) = RCT(25)*V(109)*V(114)
  A(26) = RCT(26)*V(114)*V(118)
  A(27) = RCT(27)*V(84)*V(114)
  A(28) = RCT(28)*V(84)
  A(29) = RCT(29)*V(83)*V(114)
  A(30) = RCT(30)*V(106)*V(114)
  A(31) = RCT(31)*V(112)*V(115)
  A(32) = RCT(32)*V(109)*V(112)
  A(33) = RCT(33)*V(70)
  A(34) = RCT(34)*V(70)
  A(35) = RCT(35)*V(70)*V(114)
  A(36) = RCT(36)*V(106)*V(112)
  A(37) = RCT(37)*V(112)*V(112)
  A(38) = RCT(38)*V(112)*V(112)*F(1)
  A(39) = RCT(39)*V(112)*V(118)
  A(40) = RCT(40)*V(118)*V(118)
  A(41) = RCT(41)*V(48)
  A(42) = RCT(42)*V(48)*V(114)
  A(43) = RCT(43)*V(112)*V(114)
  A(44) = RCT(44)*V(59)*V(114)
  A(45) = RCT(45)*V(114)*F(2)
  A(46) = RCT(46)*V(110)*V(115)
  A(47) = RCT(47)*V(110)*V(112)
  A(48) = RCT(48)*V(110)*V(118)
  A(49) = RCT(49)*V(110)*V(110)
  A(50) = RCT(50)*V(110)*V(110)
  A(51) = RCT(51)*V(107)*V(115)
  A(52) = RCT(52)*V(107)*V(112)
  A(53) = RCT(53)*V(107)*V(118)
  A(54) = RCT(54)*V(107)*V(110)
  A(55) = RCT(55)*V(107)*V(107)
  A(56) = RCT(56)*V(92)*V(115)
  A(57) = RCT(57)*V(92)*V(112)
  A(58) = RCT(58)*V(92)*V(118)
  A(59) = RCT(59)*V(92)*V(110)
  A(60) = RCT(60)*V(92)*V(107)
  A(61) = RCT(61)*V(92)*V(92)
  A(62) = RCT(62)*V(113)*V(115)
  A(63) = RCT(63)*V(112)*V(113)
  A(64) = RCT(64)*V(110)*V(113)
  A(65) = RCT(65)*V(113)*V(118)
  A(66) = RCT(66)*V(107)*V(113)
  A(67) = RCT(67)*V(92)*V(113)
  A(68) = RCT(68)*V(113)*V(113)
  A(69) = RCT(69)*V(109)*V(111)
  A(70) = RCT(70)*V(41)
  A(71) = RCT(71)*V(111)*V(115)
  A(72) = RCT(72)*V(111)*V(112)
  A(73) = RCT(73)*V(111)*V(118)
  A(74) = RCT(74)*V(110)*V(111)
  A(75) = RCT(75)*V(107)*V(111)
  A(76) = RCT(76)*V(92)*V(111)
  A(77) = RCT(77)*V(111)*V(113)
  A(78) = RCT(78)*V(111)*V(111)
  A(79) = RCT(79)*V(108)*V(109)
  A(80) = RCT(80)*V(42)
  A(81) = RCT(81)*V(108)*V(115)
  A(82) = RCT(82)*V(108)*V(112)
  A(83) = RCT(83)*V(108)*V(118)
  A(84) = RCT(84)*V(108)*V(110)
  A(85) = RCT(85)*V(107)*V(108)
  A(86) = RCT(86)*V(92)*V(108)
  A(87) = RCT(87)*V(108)*V(113)
  A(88) = RCT(88)*V(108)*V(111)
  A(89) = RCT(89)*V(108)*V(108)
  A(90) = RCT(90)*V(109)*V(117)
  A(91) = RCT(91)*V(43)
  A(92) = RCT(92)*V(115)*V(117)
  A(93) = RCT(93)*V(112)*V(117)
  A(94) = RCT(94)*V(117)*V(118)
  A(95) = RCT(95)*V(110)*V(117)
  A(96) = RCT(96)*V(107)*V(117)
  A(97) = RCT(97)*V(92)*V(117)
  A(98) = RCT(98)*V(113)*V(117)
  A(99) = RCT(99)*V(111)*V(117)
  A(100) = RCT(100)*V(108)*V(117)
  A(101) = RCT(101)*V(117)*V(117)
  A(102) = RCT(102)*V(109)*V(116)
  A(103) = RCT(103)*V(44)
  A(104) = RCT(104)*V(115)*V(116)
  A(105) = RCT(105)*V(112)*V(116)
  A(106) = RCT(106)*V(116)*V(118)
  A(107) = RCT(107)*V(110)*V(116)
  A(108) = RCT(108)*V(107)*V(116)
  A(109) = RCT(109)*V(92)*V(116)
  A(110) = RCT(110)*V(113)*V(116)
  A(111) = RCT(111)*V(111)*V(116)
  A(112) = RCT(112)*V(108)*V(116)
  A(113) = RCT(113)*V(116)*V(117)
  A(114) = RCT(114)*V(116)*V(116)
  A(115) = RCT(115)*V(64)*V(109)
  A(116) = RCT(116)*V(64)
  A(117) = RCT(117)*V(89)*V(109)
  A(118) = RCT(118)*V(89)*V(112)
  A(119) = RCT(119)*V(89)
  A(120) = RCT(120)*V(67)*V(109)
  A(121) = RCT(121)*V(67)*V(112)
  A(122) = RCT(122)*V(67)
  A(123) = RCT(123)*V(100)
  A(124) = RCT(124)*V(100)
  A(125) = RCT(125)*V(100)*V(114)
  A(126) = RCT(126)*V(100)*V(112)
  A(127) = RCT(127)*V(66)
  A(128) = RCT(128)*V(66)*V(115)
  A(129) = RCT(129)*V(100)*V(118)
  A(130) = RCT(130)*V(97)*V(114)
  A(131) = RCT(131)*V(97)
  A(132) = RCT(132)*V(97)*V(118)
  A(133) = RCT(133)*V(103)*V(114)
  A(134) = RCT(134)*V(103)
  A(135) = RCT(135)*V(103)*V(118)
  A(136) = RCT(136)*V(86)*V(114)
  A(137) = RCT(137)*V(86)
  A(138) = RCT(138)*V(104)*V(114)
  A(139) = RCT(139)*V(104)
  A(140) = RCT(140)*V(69)*V(114)
  A(141) = RCT(141)*V(55)*V(114)
  A(142) = RCT(142)*V(65)*V(114)
  A(143) = RCT(143)*V(65)
  A(144) = RCT(144)*V(79)*V(114)
  A(145) = RCT(145)*V(79)
  A(146) = RCT(146)*V(88)
  A(147) = RCT(147)*V(88)
  A(148) = RCT(148)*V(88)*V(114)
  A(149) = RCT(149)*V(88)*V(118)
  A(150) = RCT(150)*V(82)
  A(151) = RCT(151)*V(82)*V(114)
  A(152) = RCT(152)*V(82)*V(118)
  A(153) = RCT(153)*V(57)
  A(154) = RCT(154)*V(81)*V(114)
  A(155) = RCT(155)*V(81)*V(118)
  A(156) = RCT(156)*V(74)*V(114)
  A(157) = RCT(157)*V(74)*V(118)
  A(158) = RCT(158)*V(78)*V(118)
  A(159) = RCT(159)*V(80)*V(114)
  A(160) = RCT(160)*V(80)
  A(161) = RCT(161)*V(80)*V(118)
  A(162) = RCT(162)*V(95)*V(114)
  A(163) = RCT(163)*V(95)*V(106)
  A(164) = RCT(164)*V(95)*V(118)
  A(165) = RCT(165)*V(95)*V(102)
  A(166) = RCT(166)*V(95)
  A(167) = RCT(167)*V(99)*V(114)
  A(168) = RCT(168)*V(99)*V(106)
  A(169) = RCT(169)*V(99)*V(102)
  A(170) = RCT(170)*V(99)
  A(171) = RCT(171)*V(96)*V(114)
  A(172) = RCT(172)*V(96)*V(106)
  A(173) = RCT(173)*V(96)*V(118)
  A(174) = RCT(174)*V(96)
  A(175) = RCT(175)*V(105)*V(114)
  A(176) = RCT(176)*V(105)
  A(177) = RCT(177)*V(101)*V(114)
  A(178) = RCT(178)*V(101)
  A(179) = RCT(179)*V(77)*V(114)
  A(180) = RCT(180)*V(77)*V(106)
  A(181) = RCT(181)*V(72)*V(114)
  A(182) = RCT(182)*V(72)
  A(183) = RCT(183)*V(73)*V(114)
  A(184) = RCT(184)*V(73)
  A(185) = RCT(185)*V(32)*V(114)
  A(186) = RCT(186)*V(85)*V(114)
  A(187) = RCT(187)*V(85)*V(106)
  A(188) = RCT(188)*V(85)*V(118)
  A(189) = RCT(189)*V(85)*V(102)
  A(190) = RCT(190)*V(94)*V(114)
  A(191) = RCT(191)*V(34)*V(112)
  A(192) = RCT(192)*V(34)*V(115)
  A(193) = RCT(193)*V(34)*V(34)
  A(194) = RCT(194)*V(33)*V(114)
  A(195) = RCT(195)*V(19)*V(114)
  A(196) = RCT(196)*V(20)*V(114)
  A(197) = RCT(197)*V(18)*V(112)
  A(198) = RCT(198)*V(18)*V(115)
  A(199) = RCT(199)*V(18)
  A(200) = RCT(200)*V(94)*V(106)
  A(201) = RCT(201)*V(94)*V(118)
  A(202) = RCT(202)*V(94)*V(102)
  A(203) = RCT(203)*V(91)*V(114)
  A(204) = RCT(204)*V(91)*V(106)
  A(205) = RCT(205)*V(91)*V(118)
  A(206) = RCT(206)*V(91)*V(102)
  A(207) = RCT(207)*V(93)*V(114)
  A(208) = RCT(208)*V(93)*V(106)
  A(209) = RCT(209)*V(93)*V(118)
  A(210) = RCT(210)*V(93)*V(102)
  A(211) = RCT(211)*V(35)*V(114)
  A(212) = RCT(212)*V(56)*V(114)
  A(213) = RCT(213)*V(76)*V(114)
  A(214) = RCT(214)*V(63)*V(114)
  A(215) = RCT(215)*V(75)*V(114)
  A(216) = RCT(216)*V(62)*V(114)
  A(217) = RCT(217)*V(71)*V(114)
  A(218) = RCT(218)*V(68)*V(114)
  A(219) = RCT(219)*V(90)*V(114)
  A(220) = RCT(220)*V(90)*V(106)
  A(221) = RCT(221)*V(90)*V(118)
  A(222) = RCT(222)*V(90)*V(102)
  A(223) = RCT(223)*V(98)*V(114)
  A(224) = RCT(224)*V(98)*V(106)
  A(225) = RCT(225)*V(98)*V(118)
  A(226) = RCT(226)*V(98)*V(102)
  A(227) = RCT(227)*V(76)*V(106)
  A(228) = RCT(228)*V(87)*V(114)
  A(229) = RCT(229)*V(87)*V(106)
  A(230) = RCT(230)*V(87)*V(118)
  A(231) = RCT(231)*V(87)*V(102)
  A(232) = RCT(232)*V(61)*V(114)
  A(233) = RCT(233)*V(61)*V(114)
  A(234) = RCT(234)*V(61)*V(118)
  A(235) = RCT(235)*V(94)*V(114)
  A(236) = RCT(236)*V(94)*V(114)
  A(237) = RCT(237)*V(94)*V(114)
  A(238) = RCT(238)*V(94)*V(114)
  A(239) = RCT(239)*V(91)*V(114)
  A(240) = RCT(240)*V(91)*V(114)
  A(241) = RCT(241)*V(91)*V(114)
  A(242) = RCT(242)*V(93)*V(114)
  A(243) = RCT(243)*V(93)*V(114)
  A(244) = RCT(244)*V(94)*V(106)
  A(245) = RCT(245)*V(91)*V(106)
  A(246) = RCT(246)*V(91)*V(106)
  A(247) = RCT(247)*V(91)*V(106)
  A(248) = RCT(248)*V(91)*V(106)
  A(249) = RCT(249)*V(93)*V(106)
  A(250) = RCT(250)*V(93)*V(106)
  A(251) = RCT(251)*V(93)*V(106)
  A(252) = RCT(252)*V(94)*V(118)
  A(253) = RCT(253)*V(94)*V(118)
  A(254) = RCT(254)*V(94)*V(118)
  A(255) = RCT(255)*V(91)*V(118)
  A(256) = RCT(256)*V(91)*V(118)
  A(257) = RCT(257)*V(91)*V(118)
  A(258) = RCT(258)*V(93)*V(118)
  A(259) = RCT(259)*V(93)*V(118)
  A(260) = RCT(260)*V(93)*V(118)
  A(261) = RCT(261)*V(36)*V(114)
  A(262) = RCT(262)*V(40)*V(114)
  A(264) = RCT(264)*V(39)*V(114)
  A(265) = RCT(265)*V(39)*V(114)
  A(266) = RCT(266)*V(38)*V(114)
  A(267) = RCT(267)*V(38)*V(114)
  A(268) = RCT(268)*V(37)*V(114)
  A(269) = RCT(269)*V(37)*V(114)
  A(270) = RCT(270)*V(46)*V(114)
  A(272) = RCT(272)*V(47)*V(114)
  A(273) = RCT(273)*V(47)*V(114)
  A(274) = RCT(274)*V(45)*V(114)
  A(275) = RCT(275)*V(45)*V(114)
  A(276) = RCT(276)*V(54)*V(114)
  A(278) = RCT(278)*V(53)*V(114)
  A(279) = RCT(279)*V(53)*V(114)
  A(280) = RCT(280)*V(52)*V(114)
  A(281) = RCT(281)*V(52)*V(114)
  A(282) = RCT(282)*V(51)*V(114)
  A(284) = RCT(284)*V(50)*V(114)
  A(285) = RCT(285)*V(50)*V(114)
  A(286) = RCT(286)*V(49)*V(114)
  A(287) = RCT(287)*V(49)*V(114)
  A(288) = RCT(288)*V(26)*V(114)
  A(289) = RCT(289)*V(30)*V(114)


  Vdot(1) = A(44)
  Vdot(2) = A(128)+0.333*A(163)+0.351*A(168)+0.1*A(172)+0.37*A(187)+0.204*A(200)+0.103*A(204)+0.103*A(208)+0.297*A(213)&
              &+0.185*A(220)+0.073*A(224)+0.185*A(229)
  Vdot(3) = 0.25*A(72)+A(74)+A(75)+A(77)+0.05*A(220)+0.129*A(224)+0.17*A(229)
  Vdot(4) = 0.25*A(82)+A(84)+A(85)+A(87)+0.25*A(93)+A(95)+A(96)+A(98)+0.25*A(105)+A(107)+A(108)+2*A(110)+0.372*A(172)&
              &+0.15*A(200)+0.189*A(204)+0.189*A(208)+0.119*A(220)+0.247*A(224)
  Vdot(5) = 0.5*A(227)+0.135*A(229)
  Vdot(6) = 0.75*A(72)
  Vdot(7) = 0.75*A(82)+0.75*A(93)+0.75*A(105)
  Vdot(8) = 2*A(120)+A(230)
  Vdot(9) = 6*A(120)+7*A(160)+0.048*A(228)+0.07*A(229)+2.693*A(230)+0.55*A(231)
  Vdot(10) = A(46)+A(48)+A(51)+A(53)+A(56)+A(58)+A(62)+A(65)
  Vdot(11) = A(47)+A(49)+A(50)+A(52)+A(54)+A(55)+A(57)+A(59)+A(60)+A(61)+A(63)+A(64)+A(66)+A(67)+A(68)
  Vdot(12) = A(249)+A(258)+A(286)
  Vdot(13) = A(197)
  Vdot(14) = A(199)
  Vdot(15) = A(235)+A(274)
  Vdot(16) = A(245)+A(255)+A(280)
  Vdot(17) = A(193)+A(195)+A(196)+A(198)
  Vdot(18) = 0.12*A(194)-A(197)-A(198)-A(199)
  Vdot(19) = 0.75*A(194)-A(195)
  Vdot(20) = -A(196)
  Vdot(21) = 0.044*A(261)+A(268)
  Vdot(22) = 0.005757*A(288)
  Vdot(23) = 0.074142*A(288)
  Vdot(24) = 0.132523*A(288)
  Vdot(25) = 0.478008*A(288)
  Vdot(26) = -A(288)
  Vdot(27) = 0.081486*A(289)
  Vdot(28) = 0.5*A(289)
  Vdot(29) = 0.312248*A(289)
  Vdot(30) = -A(289)
  Vdot(31) = A(18)-A(19)-A(20)
  Vdot(32) = -A(185)
  Vdot(33) = 0.88*A(191)-A(194)
  Vdot(34) = A(190)-A(191)-A(192)-2*A(193)+0.13*A(194)
  Vdot(35) = -A(211)
  Vdot(36) = -A(261)
  Vdot(37) = 0.071*A(261)+A(266)-A(268)-A(269)
  Vdot(38) = 0.41*A(261)+A(264)-A(266)-A(267)
  Vdot(39) = 0.3*A(261)+A(262)-A(264)-A(265)
  Vdot(40) = -A(262)+A(265)+A(267)+A(269)
  Vdot(41) = A(69)-A(70)
  Vdot(42) = A(79)-A(80)
  Vdot(43) = A(90)-A(91)
  Vdot(44) = A(102)-A(103)
  Vdot(45) = A(236)+A(252)+A(272)-A(274)-A(275)
  Vdot(46) = A(238)+A(254)-A(270)+A(273)+A(275)
  Vdot(47) = A(237)+A(244)+A(253)+A(270)-A(272)-A(273)
  Vdot(48) = A(37)+A(38)-A(41)-A(42)
  Vdot(49) = A(242)+A(250)+A(259)+A(284)-A(286)-A(287)
  Vdot(50) = A(243)+A(251)+A(282)-A(284)-A(285)
  Vdot(51) = A(260)-A(282)+A(285)+A(287)
  Vdot(52) = A(239)+A(246)+A(256)+A(278)-A(280)-A(281)
  Vdot(53) = A(240)+A(247)+A(276)-A(278)-A(279)
  Vdot(54) = A(241)+A(248)+A(257)-A(276)+A(279)+A(281)
  Vdot(55) = -A(141)
  Vdot(56) = -A(212)
  Vdot(57) = -A(153)+0.031*A(204)+0.031*A(208)+0.087*A(218)
  Vdot(58) = A(11)-A(12)-A(13)
  Vdot(59) = -A(44)+A(232)+0.5*A(233)+A(234)
  Vdot(60) = A(21)-A(22)-A(23)-A(24)
  Vdot(61) = -A(232)-A(233)-A(234)
  Vdot(62) = -A(216)
  Vdot(63) = -A(214)
  Vdot(64) = -A(115)-A(116)+0.236*A(214)
  Vdot(65) = A(47)-A(142)-A(143)
  Vdot(66) = A(126)-A(127)-A(128)
  Vdot(67) = -A(120)-A(121)-A(122)+A(158)
  Vdot(68) = -A(218)
  Vdot(69) = A(49)+0.25*A(54)+0.25*A(64)-A(140)
  Vdot(70) = A(32)-A(33)-A(34)-A(35)
  Vdot(71) = -A(217)
  Vdot(72) = -A(181)-A(182)+0.108*A(217)+0.099*A(218)
  Vdot(73) = -A(183)-A(184)+0.051*A(217)+0.093*A(218)
  Vdot(74) = -A(156)-A(157)+0.207*A(217)+0.187*A(218)
  Vdot(75) = -A(215)
  Vdot(76) = -A(213)-A(227)
  Vdot(77) = -A(179)-A(180)+0.491*A(217)+0.561*A(218)
  Vdot(78) = A(117)+A(121)+A(122)-A(158)
  Vdot(79) = A(52)+A(63)-A(144)-A(145)
  Vdot(80) = -A(159)-A(160)-A(161)+0.059*A(217)+0.05*A(218)+0.061*A(223)+0.042*A(224)+0.015*A(225)
  Vdot(81) = A(118)+A(119)-A(154)-A(155)+0.017*A(217)
  Vdot(82) = -A(150)-A(151)-A(152)+0.23*A(156)+0.084*A(162)+0.9*A(163)+0.3*A(167)+0.95*A(168)+0.174*A(171)+0.742*A(172)&
               &+0.008*A(173)+0.5*A(182)+0.5*A(184)+0.119*A(217)+0.287*A(218)
  Vdot(83) = -A(29)+A(123)+A(124)+A(125)+A(129)+A(131)+0.034*A(133)+A(134)+2*A(146)+A(147)+1.26*A(148)+1.26*A(149)&
               &+A(150)+A(151)+A(152)+0.416*A(162)+0.45*A(163)+0.5*A(164)+0.67*A(166)+0.475*A(168)+0.7*A(170)+0.336*A(171)&
               &+0.498*A(172)+0.572*A(173)+1.233*A(174)+A(179)+1.5*A(180)+A(182)+A(184)+0.5*A(187)+0.491*A(189)+0.275*A(200)&
               &+0.157*A(204)+0.157*A(208)+0.393*A(213)+0.002*A(215)+0.345*A(220)+0.265*A(224)+0.012*A(226)+1.5*A(227)+0.51&
               &*A(229)
  Vdot(84) = 2*A(13)+A(25)-A(27)-A(28)+0.2*A(39)+A(129)+A(132)+A(135)+A(149)+A(152)+A(155)+A(157)+A(158)+A(161)+0.5&
               &*A(164)+0.15*A(173)+A(234)
  Vdot(85) = -A(186)-A(187)-A(188)-A(189)
  Vdot(86) = A(116)-A(136)-A(137)+0.006*A(177)+0.02*A(178)+0.13*A(204)+0.13*A(208)+0.704*A(212)+0.024*A(214)+0.452&
               &*A(215)+0.072*A(216)+0.005*A(219)+0.001*A(220)+0.024*A(221)+0.127*A(223)+0.045*A(224)+0.102*A(225)
  Vdot(87) = -A(228)-A(229)-A(230)-A(231)
  Vdot(88) = -A(146)-A(147)-A(148)-A(149)+0.23*A(154)+0.15*A(171)+0.023*A(172)+A(180)+0.5*A(182)+0.5*A(184)+0.009*A(189)&
               &+0.001*A(204)+0.001*A(208)+0.607*A(213)+0.118*A(217)+0.097*A(218)
  Vdot(89) = A(92)+A(94)+A(99)+A(100)+2*A(101)+A(113)-A(117)-A(118)-A(119)+0.24*A(154)+A(155)+0.24*A(156)+A(157)
  Vdot(90) = -A(219)-A(220)-A(221)-A(222)
  Vdot(91) = -A(203)-A(204)-A(205)-A(206)
  Vdot(92) = -A(56)-A(57)-A(58)-A(59)-A(60)-A(67)-A(76)-A(86)+A(92)+A(94)-A(97)+A(99)+A(100)+2*A(101)-A(109)+A(113)&
               &+A(136)+0.616*A(138)+0.675*A(167)+0.515*A(176)+0.596*A(177)+0.152*A(178)+A(181)+A(182)+A(183)+A(184)+0.079&
               &*A(190)+0.126*A(200)+0.187*A(201)+0.24*A(202)+0.5*A(203)+0.729*A(204)+0.75*A(205)+0.5*A(207)+0.729*A(208)&
               &+0.75*A(209)+0.559*A(214)+0.936*A(215)+0.948*A(216)+0.205*A(219)+0.488*A(221)+0.001*A(223)+0.137*A(224)&
               &+0.711*A(225)
  Vdot(93) = -A(207)-A(208)-A(209)-A(210)
  Vdot(94) = -A(190)-A(200)-A(201)-A(202)
  Vdot(95) = -A(162)-A(163)-A(164)-A(165)-A(166)+0.23*A(190)+0.39*A(200)+0.025*A(223)+0.026*A(224)+0.012*A(226)
  Vdot(96) = -A(171)-A(172)-A(173)-A(174)+0.357*A(190)+0.936*A(201)+0.025*A(223)
  Vdot(97) = A(81)+A(83)+A(88)+2*A(89)+A(100)+A(112)-A(130)-A(131)-A(132)+0.034*A(133)+A(134)+0.482*A(138)+A(139)+0.96&
               &*A(141)+0.129*A(171)+0.047*A(172)+0.467*A(174)+0.084*A(175)+0.246*A(176)+0.439*A(177)+0.431*A(178)+0.195&
               &*A(186)+0.25*A(189)+A(211)+0.445*A(214)+0.455*A(215)+0.099*A(216)+0.294*A(219)+0.154*A(220)+0.009*A(221)&
               &+0.732*A(223)+0.456*A(224)+0.507*A(225)+0.984*A(228)+0.5*A(229)
  Vdot(98) = -A(223)-A(224)-A(225)-A(226)
  Vdot(99) = -A(167)-A(168)-A(169)-A(170)+0.32*A(190)+0.16*A(200)+0.019*A(224)+0.048*A(225)
  Vdot(100) = A(46)+A(48)+A(49)+2*A(50)+0.75*A(54)+0.75*A(64)+A(74)+A(84)+A(95)+A(104)+A(106)+A(107)+A(111)+A(112)&
                &+A(113)+2*A(114)-A(123)-A(124)-A(125)-A(126)+A(127)-A(129)+A(136)+0.115*A(138)+A(140)+0.081*A(141)+0.35&
                &*A(142)+A(143)+A(147)+0.084*A(162)+0.2*A(163)+0.67*A(166)+0.3*A(167)+0.1*A(168)+0.055*A(171)+0.125*A(172)&
                &+0.227*A(173)+0.3*A(174)+0.213*A(175)+0.506*A(176)+0.01*A(177)+0.134*A(178)+1.61*A(186)+A(187)+0.191*A(189)&
                &+0.624*A(190)+0.592*A(200)+0.24*A(202)+0.276*A(203)+0.235*A(204)+0.276*A(207)+0.235*A(208)+0.096*A(213)&
                &+0.026*A(214)+0.024*A(215)+0.026*A(216)+0.732*A(219)+0.5*A(220)+0.244*A(223)+0.269*A(224)+0.079*A(225)&
                &+0.984*A(228)+0.5*A(229)
  Vdot(101) = A(62)+A(115)+0.572*A(173)-0.69*A(177)-A(178)+0.276*A(205)+0.276*A(209)+0.511*A(221)+0.321*A(225)
  Vdot(102) = A(1)-A(2)-A(3)-A(4)-A(5)-A(6)+A(16)+A(17)+A(20)-A(165)-A(169)-A(189)-A(202)-A(206)-A(210)-A(222)-A(226)&
                &-A(231)
  Vdot(103) = -A(133)-A(134)-A(135)+0.37*A(138)+A(144)+A(145)+A(165)+0.675*A(167)+0.45*A(169)+0.013*A(171)+0.218*A(173)&
                &+0.558*A(175)+0.71*A(176)+0.213*A(177)+0.147*A(178)+A(179)+A(181)+A(183)+A(188)+0.474*A(203)+0.205*A(204)&
                &+0.474*A(205)+0.147*A(206)+0.474*A(207)+0.205*A(208)+0.474*A(209)+0.147*A(210)+0.261*A(212)+0.122*A(214)&
                &+0.244*A(215)+0.204*A(216)+0.497*A(219)+0.363*A(220)+0.037*A(221)+0.45*A(222)+0.511*A(223)+0.305*A(224)&
                &+0.151*A(225)+0.069*A(226)+0.45*A(231)
  Vdot(104) = 0.5*A(64)+A(65)+0.5*A(66)+A(68)-A(138)-A(139)+0.416*A(162)+0.55*A(169)+0.15*A(171)+0.21*A(172)+0.233&
                &*A(174)+0.115*A(175)+0.177*A(177)+0.243*A(178)+0.332*A(214)+0.11*A(215)+0.089*A(216)+0.437*A(222)+0.072&
                &*A(223)+0.026*A(224)+0.001*A(225)+0.659*A(226)+0.55*A(231)
  Vdot(105) = 0.5*A(64)+0.5*A(66)+A(68)+A(77)+A(87)+A(98)+0.7*A(170)+0.332*A(171)-0.671*A(175)-A(176)+0.048*A(177)+0.435&
                &*A(178)+0.1*A(200)+0.75*A(202)+0.276*A(203)+0.276*A(204)+0.853*A(206)+0.276*A(207)+0.276*A(208)+0.853&
                &*A(210)+0.125*A(215)+0.417*A(216)+0.055*A(217)+0.119*A(219)+0.215*A(220)+0.113*A(222)+0.043*A(224)+0.259&
                &*A(226)
  Vdot(106) = A(2)-A(3)-A(7)-A(8)-A(17)-A(18)-A(30)-A(36)+0.25*A(72)+0.25*A(82)+0.25*A(93)+0.25*A(105)-A(163)-A(168)&
                &-A(172)-A(180)-A(187)-A(200)-A(204)-A(208)-A(220)-A(224)-A(227)-A(229)
  Vdot(107) = -A(51)-A(52)-A(53)-A(54)-2*A(55)-A(66)-A(75)+A(81)+A(83)-A(85)+A(88)+2*A(89)-A(96)+A(100)-A(108)+A(112)&
                &+0.034*A(133)+A(134)+0.37*A(138)+A(139)+0.05*A(141)+0.34*A(144)+0.76*A(154)+0.76*A(156)+0.5*A(162)+0.1&
                &*A(163)+0.5*A(164)+0.33*A(166)+0.3*A(167)+0.05*A(168)+0.67*A(171)+0.048*A(172)+0.799*A(173)+0.473*A(175)&
                &+0.96*A(176)+0.376*A(177)+0.564*A(178)+A(179)+A(182)+A(184)+A(186)+A(188)+0.2*A(189)+0.907*A(190)+0.066&
                &*A(200)+0.749*A(201)+0.75*A(203)+0.031*A(204)+0.276*A(205)+0.75*A(207)+0.031*A(208)+0.276*A(209)+A(211)&
                &+0.965*A(212)+0.1*A(213)+0.695*A(214)+0.835*A(215)+0.653*A(216)+0.765*A(217)+0.804*A(218)+0.91*A(219)+0.022&
                &*A(220)+0.824*A(221)+0.918*A(223)+0.033*A(224)+0.442*A(225)+0.012*A(226)+0.984*A(228)+0.949*A(230)
  Vdot(108) = -A(79)+A(80)-A(81)-A(82)-A(83)-A(84)-A(85)-A(87)-A(88)-2*A(89)-A(100)-A(112)+0.965*A(133)+A(135)+0.096&
                &*A(138)+0.37*A(148)+0.37*A(149)+0.1*A(163)+0.05*A(168)+0.048*A(172)+0.3*A(174)+0.049*A(175)+0.333*A(176)&
                &+0.201*A(204)+0.201*A(208)+0.006*A(224)
  Vdot(109) = -A(1)+A(4)-A(5)-A(6)+A(7)-A(8)+2*A(9)+2*A(10)-A(11)+A(12)+A(16)+A(23)+A(24)-A(25)+A(26)+A(28)+A(31)-A(32)&
                &+A(33)+0.61*A(34)+A(35)+0.8*A(39)+2*A(40)+A(46)+A(48)+A(51)+A(53)+A(56)+A(58)+A(65)-A(69)+A(70)+A(71)+A(73)&
                &-A(79)+A(80)+A(81)+A(83)-A(90)+A(91)+A(92)+A(94)-A(102)+A(103)+A(104)+A(106)-A(115)-A(117)-A(120)+A(128)&
                &+0.338*A(177)+A(178)+0.187*A(201)+0.474*A(205)+0.474*A(209)+0.391*A(225)
  Vdot(110) = -A(46)-A(47)-A(48)-2*A(49)-2*A(50)-A(54)-A(64)+A(71)+A(73)-A(74)+2*A(78)-A(84)+A(88)-A(95)+A(99)-A(107)&
                &+A(111)+A(116)+A(131)+A(137)+0.65*A(142)+0.3*A(170)+A(185)+0.3*A(189)+0.25*A(202)+0.011*A(215)+0.076*A(220)&
                &+0.197*A(224)+0.03*A(225)+0.26*A(229)
  Vdot(111) = -A(69)+A(70)-A(71)-A(72)-A(73)-A(74)-A(75)-A(77)-2*A(78)-A(88)-A(99)+A(104)+A(106)+A(112)+A(113)+2*A(114)&
                &+A(130)+A(132)+A(136)+A(137)+0.492*A(138)+A(139)+A(150)+A(151)+A(152)+2*A(153)+0.67*A(166)+0.675*A(167)&
                &+0.467*A(174)+0.029*A(175)+0.667*A(176)+A(181)+0.5*A(182)+A(183)+0.5*A(184)+0.123*A(204)+0.123*A(208)+0.011&
                &*A(215)+0.137*A(224)
  Vdot(112) = A(23)+A(26)+A(29)+A(30)-A(31)-A(32)+A(33)+0.61*A(34)-A(36)-2*A(37)-2*A(38)-A(39)+A(42)-A(43)+A(44)+A(45)&
                &+A(46)-A(47)+A(48)+2*A(50)+A(51)-A(52)+A(53)+A(54)+A(55)-A(63)+A(64)+A(65)+A(66)+A(68)-A(72)-A(82)-A(93)&
                &-A(105)-A(118)-A(121)+2*A(123)+A(125)-A(126)+A(127)+A(128)+A(129)+A(131)+A(134)+A(140)+0.95*A(141)+A(143)&
                &+A(145)+2*A(146)+0.63*A(148)+0.63*A(149)+A(150)+0.008*A(163)+0.34*A(166)+0.064*A(168)+0.4*A(172)+1.233&
                &*A(174)+0.379*A(175)+0.113*A(177)+0.341*A(178)+1.5*A(180)+0.5*A(182)+0.5*A(184)+0.12*A(187)+0.5*A(189)&
                &+0.907*A(190)+0.033*A(204)+0.033*A(208)+0.297*A(213)+0.224*A(217)+0.187*A(218)+0.056*A(220)+0.003*A(224)&
                &+0.013*A(226)+1.5*A(227)+0.06*A(229)+0.5*A(233)
  Vdot(113) = -A(62)-A(63)-A(64)-A(65)-A(66)-2*A(68)-A(77)-A(87)-A(98)-A(110)+0.001*A(133)+0.042*A(138)+0.025*A(167)&
                &+0.041*A(171)+0.051*A(173)+0.07*A(175)+0.04*A(176)+0.173*A(177)+0.095*A(178)+0.093*A(190)+0.008*A(200)&
                &+0.064*A(201)+0.01*A(202)+0.25*A(203)+0.18*A(204)+0.25*A(205)+0.25*A(207)+0.18*A(208)+0.25*A(209)+0.035&
                &*A(212)+0.07*A(214)+0.143*A(215)+0.347*A(216)+0.011*A(217)+0.009*A(218)+0.09*A(219)+0.001*A(220)+0.176&
                &*A(221)+0.082*A(223)+0.002*A(224)+0.136*A(225)+0.001*A(226)+0.016*A(228)+0.051*A(230)
  Vdot(114) = 2*A(19)-A(21)+A(22)-A(24)-A(25)-A(26)-A(27)+A(28)-A(29)-A(30)+A(31)+0.39*A(34)-A(35)+A(36)+0.8*A(39)+2&
                &*A(41)-A(42)-A(43)-A(44)-A(45)-A(125)-A(130)-A(133)-A(136)-A(138)-A(140)-A(141)-0.65*A(142)+A(143)-0.34&
                &*A(144)+A(145)-A(148)-A(151)-A(154)-A(156)-A(159)-A(162)+0.208*A(163)+0.33*A(166)-A(167)+0.164*A(168)&
                &-A(171)+0.285*A(172)-A(175)-A(177)-A(179)+0.5*A(180)-A(181)-A(183)-A(185)-A(186)+0.12*A(187)-A(190)+0.266&
                &*A(200)-A(203)+0.567*A(204)-A(207)+0.567*A(208)-A(211)-A(212)-0.397*A(213)-A(214)-A(215)-A(216)-A(217)&
                &-A(218)-A(219)+0.155*A(220)-A(223)+0.378*A(224)+0.5*A(227)-A(228)+0.32*A(229)-A(232)-A(233)
  Vdot(115) = A(1)-A(4)+A(5)-A(7)-A(9)-2*A(10)+A(14)+A(15)-A(21)+A(22)-A(31)-A(46)-A(51)-A(56)-A(62)-A(71)-A(81)-A(92)&
                &-A(104)-A(128)
  Vdot(116) = -A(102)+A(103)-A(104)-A(105)-A(106)-A(107)-A(108)-A(110)-A(111)-A(112)-A(113)-2*A(114)+0.5*A(162)+0.5&
                &*A(164)+0.33*A(166)+0.3*A(170)+0.289*A(171)+0.15*A(173)+0.192*A(200)+0.24*A(202)
  Vdot(117) = -A(90)+A(91)-A(92)-A(93)-A(94)-A(95)-A(96)-A(98)-A(99)-A(100)-2*A(101)-A(113)+A(159)+A(161)
  Vdot(118) = A(6)+A(8)-A(9)-A(11)+A(12)-A(14)-A(15)-A(16)-A(26)+A(27)+0.39*A(34)-A(39)-2*A(40)-A(48)-A(53)-A(58)-A(65)&
                &-A(73)-A(83)-A(94)-A(106)-A(129)-A(132)-A(135)-A(149)-A(152)-A(155)-A(157)-A(158)-A(161)-A(164)-A(173)&
                &-A(188)-A(201)-A(205)-A(209)-A(221)-A(225)-A(230)-A(234)
      
END SUBROUTINE saprc99_mosaic_4bin_iepox_vbs_aq_Fun
















SUBROUTINE saprc99_mosaic_4bin_iepox_vbs_aq_IRRFun ( V, F, RCT, IRR )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: IRR(NREACT)



  IRR(1) = RCT(1)*V(109)
  IRR(2) = RCT(2)*V(102)*F(2)
  IRR(3) = RCT(3)*V(102)*V(106)
  IRR(4) = RCT(4)*V(102)*V(115)*F(2)
  IRR(5) = RCT(5)*V(102)*V(109)
  IRR(6) = RCT(6)*V(102)*V(109)
  IRR(7) = RCT(7)*V(106)*V(115)
  IRR(8) = RCT(8)*V(106)*V(109)
  IRR(9) = RCT(9)*V(115)*V(118)
  IRR(10) = RCT(10)*V(115)*V(115)*F(2)
  IRR(11) = RCT(11)*V(109)*V(118)
  IRR(12) = RCT(12)*V(58)
  IRR(13) = RCT(13)*V(58)*F(1)
  IRR(14) = RCT(14)*V(109)*V(118)
  IRR(15) = RCT(15)*V(118)
  IRR(16) = RCT(16)*V(118)
  IRR(17) = RCT(17)*V(106)
  IRR(18) = RCT(18)*V(106)
  IRR(19) = RCT(19)*V(31)*F(1)
  IRR(20) = RCT(20)*V(31)*F(2)
  IRR(21) = RCT(21)*V(114)*V(115)
  IRR(22) = RCT(22)*V(60)
  IRR(23) = RCT(23)*V(60)
  IRR(24) = RCT(24)*V(60)*V(114)
  IRR(25) = RCT(25)*V(109)*V(114)
  IRR(26) = RCT(26)*V(114)*V(118)
  IRR(27) = RCT(27)*V(84)*V(114)
  IRR(28) = RCT(28)*V(84)
  IRR(29) = RCT(29)*V(83)*V(114)
  IRR(30) = RCT(30)*V(106)*V(114)
  IRR(31) = RCT(31)*V(112)*V(115)
  IRR(32) = RCT(32)*V(109)*V(112)
  IRR(33) = RCT(33)*V(70)
  IRR(34) = RCT(34)*V(70)
  IRR(35) = RCT(35)*V(70)*V(114)
  IRR(36) = RCT(36)*V(106)*V(112)
  IRR(37) = RCT(37)*V(112)*V(112)
  IRR(38) = RCT(38)*V(112)*V(112)*F(1)
  IRR(39) = RCT(39)*V(112)*V(118)
  IRR(40) = RCT(40)*V(118)*V(118)
  IRR(41) = RCT(41)*V(48)
  IRR(42) = RCT(42)*V(48)*V(114)
  IRR(43) = RCT(43)*V(112)*V(114)
  IRR(44) = RCT(44)*V(59)*V(114)
  IRR(45) = RCT(45)*V(114)*F(2)
  IRR(46) = RCT(46)*V(110)*V(115)
  IRR(47) = RCT(47)*V(110)*V(112)
  IRR(48) = RCT(48)*V(110)*V(118)
  IRR(49) = RCT(49)*V(110)*V(110)
  IRR(50) = RCT(50)*V(110)*V(110)
  IRR(51) = RCT(51)*V(107)*V(115)
  IRR(52) = RCT(52)*V(107)*V(112)
  IRR(53) = RCT(53)*V(107)*V(118)
  IRR(54) = RCT(54)*V(107)*V(110)
  IRR(55) = RCT(55)*V(107)*V(107)
  IRR(56) = RCT(56)*V(92)*V(115)
  IRR(57) = RCT(57)*V(92)*V(112)
  IRR(58) = RCT(58)*V(92)*V(118)
  IRR(59) = RCT(59)*V(92)*V(110)
  IRR(60) = RCT(60)*V(92)*V(107)
  IRR(61) = RCT(61)*V(92)*V(92)
  IRR(62) = RCT(62)*V(113)*V(115)
  IRR(63) = RCT(63)*V(112)*V(113)
  IRR(64) = RCT(64)*V(110)*V(113)
  IRR(65) = RCT(65)*V(113)*V(118)
  IRR(66) = RCT(66)*V(107)*V(113)
  IRR(67) = RCT(67)*V(92)*V(113)
  IRR(68) = RCT(68)*V(113)*V(113)
  IRR(69) = RCT(69)*V(109)*V(111)
  IRR(70) = RCT(70)*V(41)
  IRR(71) = RCT(71)*V(111)*V(115)
  IRR(72) = RCT(72)*V(111)*V(112)
  IRR(73) = RCT(73)*V(111)*V(118)
  IRR(74) = RCT(74)*V(110)*V(111)
  IRR(75) = RCT(75)*V(107)*V(111)
  IRR(76) = RCT(76)*V(92)*V(111)
  IRR(77) = RCT(77)*V(111)*V(113)
  IRR(78) = RCT(78)*V(111)*V(111)
  IRR(79) = RCT(79)*V(108)*V(109)
  IRR(80) = RCT(80)*V(42)
  IRR(81) = RCT(81)*V(108)*V(115)
  IRR(82) = RCT(82)*V(108)*V(112)
  IRR(83) = RCT(83)*V(108)*V(118)
  IRR(84) = RCT(84)*V(108)*V(110)
  IRR(85) = RCT(85)*V(107)*V(108)
  IRR(86) = RCT(86)*V(92)*V(108)
  IRR(87) = RCT(87)*V(108)*V(113)
  IRR(88) = RCT(88)*V(108)*V(111)
  IRR(89) = RCT(89)*V(108)*V(108)
  IRR(90) = RCT(90)*V(109)*V(117)
  IRR(91) = RCT(91)*V(43)
  IRR(92) = RCT(92)*V(115)*V(117)
  IRR(93) = RCT(93)*V(112)*V(117)
  IRR(94) = RCT(94)*V(117)*V(118)
  IRR(95) = RCT(95)*V(110)*V(117)
  IRR(96) = RCT(96)*V(107)*V(117)
  IRR(97) = RCT(97)*V(92)*V(117)
  IRR(98) = RCT(98)*V(113)*V(117)
  IRR(99) = RCT(99)*V(111)*V(117)
  IRR(100) = RCT(100)*V(108)*V(117)
  IRR(101) = RCT(101)*V(117)*V(117)
  IRR(102) = RCT(102)*V(109)*V(116)
  IRR(103) = RCT(103)*V(44)
  IRR(104) = RCT(104)*V(115)*V(116)
  IRR(105) = RCT(105)*V(112)*V(116)
  IRR(106) = RCT(106)*V(116)*V(118)
  IRR(107) = RCT(107)*V(110)*V(116)
  IRR(108) = RCT(108)*V(107)*V(116)
  IRR(109) = RCT(109)*V(92)*V(116)
  IRR(110) = RCT(110)*V(113)*V(116)
  IRR(111) = RCT(111)*V(111)*V(116)
  IRR(112) = RCT(112)*V(108)*V(116)
  IRR(113) = RCT(113)*V(116)*V(117)
  IRR(114) = RCT(114)*V(116)*V(116)
  IRR(115) = RCT(115)*V(64)*V(109)
  IRR(116) = RCT(116)*V(64)
  IRR(117) = RCT(117)*V(89)*V(109)
  IRR(118) = RCT(118)*V(89)*V(112)
  IRR(119) = RCT(119)*V(89)
  IRR(120) = RCT(120)*V(67)*V(109)
  IRR(121) = RCT(121)*V(67)*V(112)
  IRR(122) = RCT(122)*V(67)
  IRR(123) = RCT(123)*V(100)
  IRR(124) = RCT(124)*V(100)
  IRR(125) = RCT(125)*V(100)*V(114)
  IRR(126) = RCT(126)*V(100)*V(112)
  IRR(127) = RCT(127)*V(66)
  IRR(128) = RCT(128)*V(66)*V(115)
  IRR(129) = RCT(129)*V(100)*V(118)
  IRR(130) = RCT(130)*V(97)*V(114)
  IRR(131) = RCT(131)*V(97)
  IRR(132) = RCT(132)*V(97)*V(118)
  IRR(133) = RCT(133)*V(103)*V(114)
  IRR(134) = RCT(134)*V(103)
  IRR(135) = RCT(135)*V(103)*V(118)
  IRR(136) = RCT(136)*V(86)*V(114)
  IRR(137) = RCT(137)*V(86)
  IRR(138) = RCT(138)*V(104)*V(114)
  IRR(139) = RCT(139)*V(104)
  IRR(140) = RCT(140)*V(69)*V(114)
  IRR(141) = RCT(141)*V(55)*V(114)
  IRR(142) = RCT(142)*V(65)*V(114)
  IRR(143) = RCT(143)*V(65)
  IRR(144) = RCT(144)*V(79)*V(114)
  IRR(145) = RCT(145)*V(79)
  IRR(146) = RCT(146)*V(88)
  IRR(147) = RCT(147)*V(88)
  IRR(148) = RCT(148)*V(88)*V(114)
  IRR(149) = RCT(149)*V(88)*V(118)
  IRR(150) = RCT(150)*V(82)
  IRR(151) = RCT(151)*V(82)*V(114)
  IRR(152) = RCT(152)*V(82)*V(118)
  IRR(153) = RCT(153)*V(57)
  IRR(154) = RCT(154)*V(81)*V(114)
  IRR(155) = RCT(155)*V(81)*V(118)
  IRR(156) = RCT(156)*V(74)*V(114)
  IRR(157) = RCT(157)*V(74)*V(118)
  IRR(158) = RCT(158)*V(78)*V(118)
  IRR(159) = RCT(159)*V(80)*V(114)
  IRR(160) = RCT(160)*V(80)
  IRR(161) = RCT(161)*V(80)*V(118)
  IRR(162) = RCT(162)*V(95)*V(114)
  IRR(163) = RCT(163)*V(95)*V(106)
  IRR(164) = RCT(164)*V(95)*V(118)
  IRR(165) = RCT(165)*V(95)*V(102)
  IRR(166) = RCT(166)*V(95)
  IRR(167) = RCT(167)*V(99)*V(114)
  IRR(168) = RCT(168)*V(99)*V(106)
  IRR(169) = RCT(169)*V(99)*V(102)
  IRR(170) = RCT(170)*V(99)
  IRR(171) = RCT(171)*V(96)*V(114)
  IRR(172) = RCT(172)*V(96)*V(106)
  IRR(173) = RCT(173)*V(96)*V(118)
  IRR(174) = RCT(174)*V(96)
  IRR(175) = RCT(175)*V(105)*V(114)
  IRR(176) = RCT(176)*V(105)
  IRR(177) = RCT(177)*V(101)*V(114)
  IRR(178) = RCT(178)*V(101)
  IRR(179) = RCT(179)*V(77)*V(114)
  IRR(180) = RCT(180)*V(77)*V(106)
  IRR(181) = RCT(181)*V(72)*V(114)
  IRR(182) = RCT(182)*V(72)
  IRR(183) = RCT(183)*V(73)*V(114)
  IRR(184) = RCT(184)*V(73)
  IRR(185) = RCT(185)*V(32)*V(114)
  IRR(186) = RCT(186)*V(85)*V(114)
  IRR(187) = RCT(187)*V(85)*V(106)
  IRR(188) = RCT(188)*V(85)*V(118)
  IRR(189) = RCT(189)*V(85)*V(102)
  IRR(190) = RCT(190)*V(94)*V(114)
  IRR(191) = RCT(191)*V(34)*V(112)
  IRR(192) = RCT(192)*V(34)*V(115)
  IRR(193) = RCT(193)*V(34)*V(34)
  IRR(194) = RCT(194)*V(33)*V(114)
  IRR(195) = RCT(195)*V(19)*V(114)
  IRR(196) = RCT(196)*V(20)*V(114)
  IRR(197) = RCT(197)*V(18)*V(112)
  IRR(198) = RCT(198)*V(18)*V(115)
  IRR(199) = RCT(199)*V(18)
  IRR(200) = RCT(200)*V(94)*V(106)
  IRR(201) = RCT(201)*V(94)*V(118)
  IRR(202) = RCT(202)*V(94)*V(102)
  IRR(203) = RCT(203)*V(91)*V(114)
  IRR(204) = RCT(204)*V(91)*V(106)
  IRR(205) = RCT(205)*V(91)*V(118)
  IRR(206) = RCT(206)*V(91)*V(102)
  IRR(207) = RCT(207)*V(93)*V(114)
  IRR(208) = RCT(208)*V(93)*V(106)
  IRR(209) = RCT(209)*V(93)*V(118)
  IRR(210) = RCT(210)*V(93)*V(102)
  IRR(211) = RCT(211)*V(35)*V(114)
  IRR(212) = RCT(212)*V(56)*V(114)
  IRR(213) = RCT(213)*V(76)*V(114)
  IRR(214) = RCT(214)*V(63)*V(114)
  IRR(215) = RCT(215)*V(75)*V(114)
  IRR(216) = RCT(216)*V(62)*V(114)
  IRR(217) = RCT(217)*V(71)*V(114)
  IRR(218) = RCT(218)*V(68)*V(114)
  IRR(219) = RCT(219)*V(90)*V(114)
  IRR(220) = RCT(220)*V(90)*V(106)
  IRR(221) = RCT(221)*V(90)*V(118)
  IRR(222) = RCT(222)*V(90)*V(102)
  IRR(223) = RCT(223)*V(98)*V(114)
  IRR(224) = RCT(224)*V(98)*V(106)
  IRR(225) = RCT(225)*V(98)*V(118)
  IRR(226) = RCT(226)*V(98)*V(102)
  IRR(227) = RCT(227)*V(76)*V(106)
  IRR(228) = RCT(228)*V(87)*V(114)
  IRR(229) = RCT(229)*V(87)*V(106)
  IRR(230) = RCT(230)*V(87)*V(118)
  IRR(231) = RCT(231)*V(87)*V(102)
  IRR(232) = RCT(232)*V(61)*V(114)
  IRR(233) = RCT(233)*V(61)*V(114)
  IRR(234) = RCT(234)*V(61)*V(118)
  IRR(235) = RCT(235)*V(94)*V(114)
  IRR(236) = RCT(236)*V(94)*V(114)
  IRR(237) = RCT(237)*V(94)*V(114)
  IRR(238) = RCT(238)*V(94)*V(114)
  IRR(239) = RCT(239)*V(91)*V(114)
  IRR(240) = RCT(240)*V(91)*V(114)
  IRR(241) = RCT(241)*V(91)*V(114)
  IRR(242) = RCT(242)*V(93)*V(114)
  IRR(243) = RCT(243)*V(93)*V(114)
  IRR(244) = RCT(244)*V(94)*V(106)
  IRR(245) = RCT(245)*V(91)*V(106)
  IRR(246) = RCT(246)*V(91)*V(106)
  IRR(247) = RCT(247)*V(91)*V(106)
  IRR(248) = RCT(248)*V(91)*V(106)
  IRR(249) = RCT(249)*V(93)*V(106)
  IRR(250) = RCT(250)*V(93)*V(106)
  IRR(251) = RCT(251)*V(93)*V(106)
  IRR(252) = RCT(252)*V(94)*V(118)
  IRR(253) = RCT(253)*V(94)*V(118)
  IRR(254) = RCT(254)*V(94)*V(118)
  IRR(255) = RCT(255)*V(91)*V(118)
  IRR(256) = RCT(256)*V(91)*V(118)
  IRR(257) = RCT(257)*V(91)*V(118)
  IRR(258) = RCT(258)*V(93)*V(118)
  IRR(259) = RCT(259)*V(93)*V(118)
  IRR(260) = RCT(260)*V(93)*V(118)
  IRR(261) = RCT(261)*V(36)*V(114)
  IRR(262) = RCT(262)*V(40)*V(114)
  IRR(264) = RCT(264)*V(39)*V(114)
  IRR(265) = RCT(265)*V(39)*V(114)
  IRR(266) = RCT(266)*V(38)*V(114)
  IRR(267) = RCT(267)*V(38)*V(114)
  IRR(268) = RCT(268)*V(37)*V(114)
  IRR(269) = RCT(269)*V(37)*V(114)
  IRR(270) = RCT(270)*V(46)*V(114)
  IRR(272) = RCT(272)*V(47)*V(114)
  IRR(273) = RCT(273)*V(47)*V(114)
  IRR(274) = RCT(274)*V(45)*V(114)
  IRR(275) = RCT(275)*V(45)*V(114)
  IRR(276) = RCT(276)*V(54)*V(114)
  IRR(278) = RCT(278)*V(53)*V(114)
  IRR(279) = RCT(279)*V(53)*V(114)
  IRR(280) = RCT(280)*V(52)*V(114)
  IRR(281) = RCT(281)*V(52)*V(114)
  IRR(282) = RCT(282)*V(51)*V(114)
  IRR(284) = RCT(284)*V(50)*V(114)
  IRR(285) = RCT(285)*V(50)*V(114)
  IRR(286) = RCT(286)*V(49)*V(114)
  IRR(287) = RCT(287)*V(49)*V(114)
  IRR(288) = RCT(288)*V(26)*V(114)
  IRR(289) = RCT(289)*V(30)*V(114)
      
END SUBROUTINE saprc99_mosaic_4bin_iepox_vbs_aq_IRRFun
















SUBROUTINE saprc99_mosaic_4bin_iepox_vbs_aq_Jac_SP ( V, F, RCT, JVS )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: JVS(LU_NONZERO)




  REAL(kind=dp) :: B(526)


  B(1) = RCT(1)

  B(2) = RCT(2)*F(2)

  B(4) = RCT(3)*V(106)

  B(5) = RCT(3)*V(102)

  B(6) = RCT(4)*V(115)*F(2)

  B(7) = RCT(4)*V(102)*F(2)

  B(9) = RCT(5)*V(109)

  B(10) = RCT(5)*V(102)

  B(11) = RCT(6)*V(109)

  B(12) = RCT(6)*V(102)

  B(13) = RCT(7)*V(115)

  B(14) = RCT(7)*V(106)

  B(15) = RCT(8)*V(109)

  B(16) = RCT(8)*V(106)

  B(17) = RCT(9)*V(118)

  B(18) = RCT(9)*V(115)

  B(19) = RCT(10)*2*V(115)*F(2)

  B(21) = RCT(11)*V(118)

  B(22) = RCT(11)*V(109)

  B(23) = RCT(12)

  B(24) = RCT(13)*F(1)

  B(26) = RCT(14)*V(118)

  B(27) = RCT(14)*V(109)

  B(28) = RCT(15)

  B(29) = RCT(16)

  B(30) = RCT(17)

  B(31) = RCT(18)

  B(32) = RCT(19)*F(1)

  B(34) = RCT(20)*F(2)

  B(36) = RCT(21)*V(115)

  B(37) = RCT(21)*V(114)

  B(38) = RCT(22)

  B(39) = RCT(23)

  B(40) = RCT(24)*V(114)

  B(41) = RCT(24)*V(60)

  B(42) = RCT(25)*V(114)

  B(43) = RCT(25)*V(109)

  B(44) = RCT(26)*V(118)

  B(45) = RCT(26)*V(114)

  B(46) = RCT(27)*V(114)

  B(47) = RCT(27)*V(84)

  B(48) = RCT(28)

  B(49) = RCT(29)*V(114)

  B(50) = RCT(29)*V(83)

  B(51) = RCT(30)*V(114)

  B(52) = RCT(30)*V(106)

  B(53) = RCT(31)*V(115)

  B(54) = RCT(31)*V(112)

  B(55) = RCT(32)*V(112)

  B(56) = RCT(32)*V(109)

  B(57) = RCT(33)

  B(58) = RCT(34)

  B(59) = RCT(35)*V(114)

  B(60) = RCT(35)*V(70)

  B(61) = RCT(36)*V(112)

  B(62) = RCT(36)*V(106)

  B(63) = RCT(37)*2*V(112)

  B(64) = RCT(38)*2*V(112)*F(1)

  B(66) = RCT(39)*V(118)

  B(67) = RCT(39)*V(112)

  B(68) = RCT(40)*2*V(118)

  B(69) = RCT(41)

  B(70) = RCT(42)*V(114)

  B(71) = RCT(42)*V(48)

  B(72) = RCT(43)*V(114)

  B(73) = RCT(43)*V(112)

  B(74) = RCT(44)*V(114)

  B(75) = RCT(44)*V(59)

  B(76) = RCT(45)*F(2)

  B(78) = RCT(46)*V(115)

  B(79) = RCT(46)*V(110)

  B(80) = RCT(47)*V(112)

  B(81) = RCT(47)*V(110)

  B(82) = RCT(48)*V(118)

  B(83) = RCT(48)*V(110)

  B(84) = RCT(49)*2*V(110)

  B(85) = RCT(50)*2*V(110)

  B(86) = RCT(51)*V(115)

  B(87) = RCT(51)*V(107)

  B(88) = RCT(52)*V(112)

  B(89) = RCT(52)*V(107)

  B(90) = RCT(53)*V(118)

  B(91) = RCT(53)*V(107)

  B(92) = RCT(54)*V(110)

  B(93) = RCT(54)*V(107)

  B(94) = RCT(55)*2*V(107)

  B(95) = RCT(56)*V(115)

  B(96) = RCT(56)*V(92)

  B(97) = RCT(57)*V(112)

  B(98) = RCT(57)*V(92)

  B(99) = RCT(58)*V(118)

  B(100) = RCT(58)*V(92)

  B(101) = RCT(59)*V(110)

  B(102) = RCT(59)*V(92)

  B(103) = RCT(60)*V(107)

  B(104) = RCT(60)*V(92)

  B(105) = RCT(61)*2*V(92)

  B(106) = RCT(62)*V(115)

  B(107) = RCT(62)*V(113)

  B(108) = RCT(63)*V(113)

  B(109) = RCT(63)*V(112)

  B(110) = RCT(64)*V(113)

  B(111) = RCT(64)*V(110)

  B(112) = RCT(65)*V(118)

  B(113) = RCT(65)*V(113)

  B(114) = RCT(66)*V(113)

  B(115) = RCT(66)*V(107)

  B(116) = RCT(67)*V(113)

  B(117) = RCT(67)*V(92)

  B(118) = RCT(68)*2*V(113)

  B(119) = RCT(69)*V(111)

  B(120) = RCT(69)*V(109)

  B(121) = RCT(70)

  B(122) = RCT(71)*V(115)

  B(123) = RCT(71)*V(111)

  B(124) = RCT(72)*V(112)

  B(125) = RCT(72)*V(111)

  B(126) = RCT(73)*V(118)

  B(127) = RCT(73)*V(111)

  B(128) = RCT(74)*V(111)

  B(129) = RCT(74)*V(110)

  B(130) = RCT(75)*V(111)

  B(131) = RCT(75)*V(107)

  B(132) = RCT(76)*V(111)

  B(133) = RCT(76)*V(92)

  B(134) = RCT(77)*V(113)

  B(135) = RCT(77)*V(111)

  B(136) = RCT(78)*2*V(111)

  B(137) = RCT(79)*V(109)

  B(138) = RCT(79)*V(108)

  B(139) = RCT(80)

  B(140) = RCT(81)*V(115)

  B(141) = RCT(81)*V(108)

  B(142) = RCT(82)*V(112)

  B(143) = RCT(82)*V(108)

  B(144) = RCT(83)*V(118)

  B(145) = RCT(83)*V(108)

  B(146) = RCT(84)*V(110)

  B(147) = RCT(84)*V(108)

  B(148) = RCT(85)*V(108)

  B(149) = RCT(85)*V(107)

  B(150) = RCT(86)*V(108)

  B(151) = RCT(86)*V(92)

  B(152) = RCT(87)*V(113)

  B(153) = RCT(87)*V(108)

  B(154) = RCT(88)*V(111)

  B(155) = RCT(88)*V(108)

  B(156) = RCT(89)*2*V(108)

  B(157) = RCT(90)*V(117)

  B(158) = RCT(90)*V(109)

  B(159) = RCT(91)

  B(160) = RCT(92)*V(117)

  B(161) = RCT(92)*V(115)

  B(162) = RCT(93)*V(117)

  B(163) = RCT(93)*V(112)

  B(164) = RCT(94)*V(118)

  B(165) = RCT(94)*V(117)

  B(166) = RCT(95)*V(117)

  B(167) = RCT(95)*V(110)

  B(168) = RCT(96)*V(117)

  B(169) = RCT(96)*V(107)

  B(170) = RCT(97)*V(117)

  B(171) = RCT(97)*V(92)

  B(172) = RCT(98)*V(117)

  B(173) = RCT(98)*V(113)

  B(174) = RCT(99)*V(117)

  B(175) = RCT(99)*V(111)

  B(176) = RCT(100)*V(117)

  B(177) = RCT(100)*V(108)

  B(178) = RCT(101)*2*V(117)

  B(179) = RCT(102)*V(116)

  B(180) = RCT(102)*V(109)

  B(181) = RCT(103)

  B(182) = RCT(104)*V(116)

  B(183) = RCT(104)*V(115)

  B(184) = RCT(105)*V(116)

  B(185) = RCT(105)*V(112)

  B(186) = RCT(106)*V(118)

  B(187) = RCT(106)*V(116)

  B(188) = RCT(107)*V(116)

  B(189) = RCT(107)*V(110)

  B(190) = RCT(108)*V(116)

  B(191) = RCT(108)*V(107)

  B(192) = RCT(109)*V(116)

  B(193) = RCT(109)*V(92)

  B(194) = RCT(110)*V(116)

  B(195) = RCT(110)*V(113)

  B(196) = RCT(111)*V(116)

  B(197) = RCT(111)*V(111)

  B(198) = RCT(112)*V(116)

  B(199) = RCT(112)*V(108)

  B(200) = RCT(113)*V(117)

  B(201) = RCT(113)*V(116)

  B(202) = RCT(114)*2*V(116)

  B(203) = RCT(115)*V(109)

  B(204) = RCT(115)*V(64)

  B(205) = RCT(116)

  B(206) = RCT(117)*V(109)

  B(207) = RCT(117)*V(89)

  B(208) = RCT(118)*V(112)

  B(209) = RCT(118)*V(89)

  B(210) = RCT(119)

  B(211) = RCT(120)*V(109)

  B(212) = RCT(120)*V(67)

  B(213) = RCT(121)*V(112)

  B(214) = RCT(121)*V(67)

  B(215) = RCT(122)

  B(216) = RCT(123)

  B(217) = RCT(124)

  B(218) = RCT(125)*V(114)

  B(219) = RCT(125)*V(100)

  B(220) = RCT(126)*V(112)

  B(221) = RCT(126)*V(100)

  B(222) = RCT(127)

  B(223) = RCT(128)*V(115)

  B(224) = RCT(128)*V(66)

  B(225) = RCT(129)*V(118)

  B(226) = RCT(129)*V(100)

  B(227) = RCT(130)*V(114)

  B(228) = RCT(130)*V(97)

  B(229) = RCT(131)

  B(230) = RCT(132)*V(118)

  B(231) = RCT(132)*V(97)

  B(232) = RCT(133)*V(114)

  B(233) = RCT(133)*V(103)

  B(234) = RCT(134)

  B(235) = RCT(135)*V(118)

  B(236) = RCT(135)*V(103)

  B(237) = RCT(136)*V(114)

  B(238) = RCT(136)*V(86)

  B(239) = RCT(137)

  B(240) = RCT(138)*V(114)

  B(241) = RCT(138)*V(104)

  B(242) = RCT(139)

  B(243) = RCT(140)*V(114)

  B(244) = RCT(140)*V(69)

  B(245) = RCT(141)*V(114)

  B(246) = RCT(141)*V(55)

  B(247) = RCT(142)*V(114)

  B(248) = RCT(142)*V(65)

  B(249) = RCT(143)

  B(250) = RCT(144)*V(114)

  B(251) = RCT(144)*V(79)

  B(252) = RCT(145)

  B(253) = RCT(146)

  B(254) = RCT(147)

  B(255) = RCT(148)*V(114)

  B(256) = RCT(148)*V(88)

  B(257) = RCT(149)*V(118)

  B(258) = RCT(149)*V(88)

  B(259) = RCT(150)

  B(260) = RCT(151)*V(114)

  B(261) = RCT(151)*V(82)

  B(262) = RCT(152)*V(118)

  B(263) = RCT(152)*V(82)

  B(264) = RCT(153)

  B(265) = RCT(154)*V(114)

  B(266) = RCT(154)*V(81)

  B(267) = RCT(155)*V(118)

  B(268) = RCT(155)*V(81)

  B(269) = RCT(156)*V(114)

  B(270) = RCT(156)*V(74)

  B(271) = RCT(157)*V(118)

  B(272) = RCT(157)*V(74)

  B(273) = RCT(158)*V(118)

  B(274) = RCT(158)*V(78)

  B(275) = RCT(159)*V(114)

  B(276) = RCT(159)*V(80)

  B(277) = RCT(160)

  B(278) = RCT(161)*V(118)

  B(279) = RCT(161)*V(80)

  B(280) = RCT(162)*V(114)

  B(281) = RCT(162)*V(95)

  B(282) = RCT(163)*V(106)

  B(283) = RCT(163)*V(95)

  B(284) = RCT(164)*V(118)

  B(285) = RCT(164)*V(95)

  B(286) = RCT(165)*V(102)

  B(287) = RCT(165)*V(95)

  B(288) = RCT(166)

  B(289) = RCT(167)*V(114)

  B(290) = RCT(167)*V(99)

  B(291) = RCT(168)*V(106)

  B(292) = RCT(168)*V(99)

  B(293) = RCT(169)*V(102)

  B(294) = RCT(169)*V(99)

  B(295) = RCT(170)

  B(296) = RCT(171)*V(114)

  B(297) = RCT(171)*V(96)

  B(298) = RCT(172)*V(106)

  B(299) = RCT(172)*V(96)

  B(300) = RCT(173)*V(118)

  B(301) = RCT(173)*V(96)

  B(302) = RCT(174)

  B(303) = RCT(175)*V(114)

  B(304) = RCT(175)*V(105)

  B(305) = RCT(176)

  B(306) = RCT(177)*V(114)

  B(307) = RCT(177)*V(101)

  B(308) = RCT(178)

  B(309) = RCT(179)*V(114)

  B(310) = RCT(179)*V(77)

  B(311) = RCT(180)*V(106)

  B(312) = RCT(180)*V(77)

  B(313) = RCT(181)*V(114)

  B(314) = RCT(181)*V(72)

  B(315) = RCT(182)

  B(316) = RCT(183)*V(114)

  B(317) = RCT(183)*V(73)

  B(318) = RCT(184)

  B(319) = RCT(185)*V(114)

  B(320) = RCT(185)*V(32)

  B(321) = RCT(186)*V(114)

  B(322) = RCT(186)*V(85)

  B(323) = RCT(187)*V(106)

  B(324) = RCT(187)*V(85)

  B(325) = RCT(188)*V(118)

  B(326) = RCT(188)*V(85)

  B(327) = RCT(189)*V(102)

  B(328) = RCT(189)*V(85)

  B(329) = RCT(190)*V(114)

  B(330) = RCT(190)*V(94)

  B(331) = RCT(191)*V(112)

  B(332) = RCT(191)*V(34)

  B(333) = RCT(192)*V(115)

  B(334) = RCT(192)*V(34)

  B(335) = RCT(193)*2*V(34)

  B(336) = RCT(194)*V(114)

  B(337) = RCT(194)*V(33)

  B(338) = RCT(195)*V(114)

  B(339) = RCT(195)*V(19)

  B(340) = RCT(196)*V(114)

  B(341) = RCT(196)*V(20)

  B(342) = RCT(197)*V(112)

  B(343) = RCT(197)*V(18)

  B(344) = RCT(198)*V(115)

  B(345) = RCT(198)*V(18)

  B(346) = RCT(199)

  B(347) = RCT(200)*V(106)

  B(348) = RCT(200)*V(94)

  B(349) = RCT(201)*V(118)

  B(350) = RCT(201)*V(94)

  B(351) = RCT(202)*V(102)

  B(352) = RCT(202)*V(94)

  B(353) = RCT(203)*V(114)

  B(354) = RCT(203)*V(91)

  B(355) = RCT(204)*V(106)

  B(356) = RCT(204)*V(91)

  B(357) = RCT(205)*V(118)

  B(358) = RCT(205)*V(91)

  B(359) = RCT(206)*V(102)

  B(360) = RCT(206)*V(91)

  B(361) = RCT(207)*V(114)

  B(362) = RCT(207)*V(93)

  B(363) = RCT(208)*V(106)

  B(364) = RCT(208)*V(93)

  B(365) = RCT(209)*V(118)

  B(366) = RCT(209)*V(93)

  B(367) = RCT(210)*V(102)

  B(368) = RCT(210)*V(93)

  B(369) = RCT(211)*V(114)

  B(370) = RCT(211)*V(35)

  B(371) = RCT(212)*V(114)

  B(372) = RCT(212)*V(56)

  B(373) = RCT(213)*V(114)

  B(374) = RCT(213)*V(76)

  B(375) = RCT(214)*V(114)

  B(376) = RCT(214)*V(63)

  B(377) = RCT(215)*V(114)

  B(378) = RCT(215)*V(75)

  B(379) = RCT(216)*V(114)

  B(380) = RCT(216)*V(62)

  B(381) = RCT(217)*V(114)

  B(382) = RCT(217)*V(71)

  B(383) = RCT(218)*V(114)

  B(384) = RCT(218)*V(68)

  B(385) = RCT(219)*V(114)

  B(386) = RCT(219)*V(90)

  B(387) = RCT(220)*V(106)

  B(388) = RCT(220)*V(90)

  B(389) = RCT(221)*V(118)

  B(390) = RCT(221)*V(90)

  B(391) = RCT(222)*V(102)

  B(392) = RCT(222)*V(90)

  B(393) = RCT(223)*V(114)

  B(394) = RCT(223)*V(98)

  B(395) = RCT(224)*V(106)

  B(396) = RCT(224)*V(98)

  B(397) = RCT(225)*V(118)

  B(398) = RCT(225)*V(98)

  B(399) = RCT(226)*V(102)

  B(400) = RCT(226)*V(98)

  B(401) = RCT(227)*V(106)

  B(402) = RCT(227)*V(76)

  B(403) = RCT(228)*V(114)

  B(404) = RCT(228)*V(87)

  B(405) = RCT(229)*V(106)

  B(406) = RCT(229)*V(87)

  B(407) = RCT(230)*V(118)

  B(408) = RCT(230)*V(87)

  B(409) = RCT(231)*V(102)

  B(410) = RCT(231)*V(87)

  B(411) = RCT(232)*V(114)

  B(412) = RCT(232)*V(61)

  B(413) = RCT(233)*V(114)

  B(414) = RCT(233)*V(61)

  B(415) = RCT(234)*V(118)

  B(416) = RCT(234)*V(61)

  B(417) = RCT(235)*V(114)

  B(418) = RCT(235)*V(94)

  B(419) = RCT(236)*V(114)

  B(420) = RCT(236)*V(94)

  B(421) = RCT(237)*V(114)

  B(422) = RCT(237)*V(94)

  B(423) = RCT(238)*V(114)

  B(424) = RCT(238)*V(94)

  B(425) = RCT(239)*V(114)

  B(426) = RCT(239)*V(91)

  B(427) = RCT(240)*V(114)

  B(428) = RCT(240)*V(91)

  B(429) = RCT(241)*V(114)

  B(430) = RCT(241)*V(91)

  B(431) = RCT(242)*V(114)

  B(432) = RCT(242)*V(93)

  B(433) = RCT(243)*V(114)

  B(434) = RCT(243)*V(93)

  B(435) = RCT(244)*V(106)

  B(436) = RCT(244)*V(94)

  B(437) = RCT(245)*V(106)

  B(438) = RCT(245)*V(91)

  B(439) = RCT(246)*V(106)

  B(440) = RCT(246)*V(91)

  B(441) = RCT(247)*V(106)

  B(442) = RCT(247)*V(91)

  B(443) = RCT(248)*V(106)

  B(444) = RCT(248)*V(91)

  B(445) = RCT(249)*V(106)

  B(446) = RCT(249)*V(93)

  B(447) = RCT(250)*V(106)

  B(448) = RCT(250)*V(93)

  B(449) = RCT(251)*V(106)

  B(450) = RCT(251)*V(93)

  B(451) = RCT(252)*V(118)

  B(452) = RCT(252)*V(94)

  B(453) = RCT(253)*V(118)

  B(454) = RCT(253)*V(94)

  B(455) = RCT(254)*V(118)

  B(456) = RCT(254)*V(94)

  B(457) = RCT(255)*V(118)

  B(458) = RCT(255)*V(91)

  B(459) = RCT(256)*V(118)

  B(460) = RCT(256)*V(91)

  B(461) = RCT(257)*V(118)

  B(462) = RCT(257)*V(91)

  B(463) = RCT(258)*V(118)

  B(464) = RCT(258)*V(93)

  B(465) = RCT(259)*V(118)

  B(466) = RCT(259)*V(93)

  B(467) = RCT(260)*V(118)

  B(468) = RCT(260)*V(93)

  B(469) = RCT(261)*V(114)

  B(470) = RCT(261)*V(36)

  B(471) = RCT(262)*V(114)

  B(472) = RCT(262)*V(40)

  B(473) = RCT(263)*V(114)

  B(474) = RCT(263)*V(40)

  B(475) = RCT(264)*V(114)

  B(476) = RCT(264)*V(39)

  B(477) = RCT(265)*V(114)

  B(478) = RCT(265)*V(39)

  B(479) = RCT(266)*V(114)

  B(480) = RCT(266)*V(38)

  B(481) = RCT(267)*V(114)

  B(482) = RCT(267)*V(38)

  B(483) = RCT(268)*V(114)

  B(484) = RCT(268)*V(37)

  B(485) = RCT(269)*V(114)

  B(486) = RCT(269)*V(37)

  B(487) = RCT(270)*V(114)

  B(488) = RCT(270)*V(46)

  B(489) = RCT(271)*V(114)

  B(490) = RCT(271)*V(46)

  B(491) = RCT(272)*V(114)

  B(492) = RCT(272)*V(47)

  B(493) = RCT(273)*V(114)

  B(494) = RCT(273)*V(47)

  B(495) = RCT(274)*V(114)

  B(496) = RCT(274)*V(45)

  B(497) = RCT(275)*V(114)

  B(498) = RCT(275)*V(45)

  B(499) = RCT(276)*V(114)

  B(500) = RCT(276)*V(54)

  B(501) = RCT(277)*V(114)

  B(502) = RCT(277)*V(54)

  B(503) = RCT(278)*V(114)

  B(504) = RCT(278)*V(53)

  B(505) = RCT(279)*V(114)

  B(506) = RCT(279)*V(53)

  B(507) = RCT(280)*V(114)

  B(508) = RCT(280)*V(52)

  B(509) = RCT(281)*V(114)

  B(510) = RCT(281)*V(52)

  B(511) = RCT(282)*V(114)

  B(512) = RCT(282)*V(51)

  B(513) = RCT(283)*V(114)

  B(514) = RCT(283)*V(51)

  B(515) = RCT(284)*V(114)

  B(516) = RCT(284)*V(50)

  B(517) = RCT(285)*V(114)

  B(518) = RCT(285)*V(50)

  B(519) = RCT(286)*V(114)

  B(520) = RCT(286)*V(49)

  B(521) = RCT(287)*V(114)

  B(522) = RCT(287)*V(49)

  B(523) = RCT(288)*V(114)

  B(524) = RCT(288)*V(26)

  B(525) = RCT(289)*V(114)

  B(526) = RCT(289)*V(30)



  JVS(1) = 0

  JVS(2) = B(74)

  JVS(3) = B(75)

  JVS(4) = 0

  JVS(5) = B(223)

  JVS(6) = 0.297*B(373)

  JVS(7) = 0.37*B(323)

  JVS(8) = 0.185*B(405)

  JVS(9) = 0.185*B(387)

  JVS(10) = 0.103*B(355)

  JVS(11) = 0.103*B(363)

  JVS(12) = 0.204*B(347)

  JVS(13) = 0.333*B(282)

  JVS(14) = 0.1*B(298)

  JVS(15) = 0.073*B(395)

  JVS(16) = 0.351*B(291)

  JVS(17) = 0.333*B(283)+0.351*B(292)+0.1*B(299)+0.37*B(324)+0.204*B(348)+0.103*B(356)+0.103*B(364)+0.185*B(388)+0.073&
              &*B(396)+0.185*B(406)

  JVS(18) = 0.297*B(374)

  JVS(19) = B(224)

  JVS(20) = 0

  JVS(21) = 0.17*B(405)

  JVS(22) = 0.05*B(387)

  JVS(23) = 0.129*B(395)

  JVS(24) = 0.05*B(388)+0.129*B(396)+0.17*B(406)

  JVS(25) = B(130)

  JVS(26) = B(128)

  JVS(27) = 0.25*B(124)+B(129)+B(131)+B(134)

  JVS(28) = 0.25*B(125)

  JVS(29) = B(135)

  JVS(30) = 0

  JVS(31) = 0.119*B(387)

  JVS(32) = 0.189*B(355)

  JVS(33) = 0.189*B(363)

  JVS(34) = 0.15*B(347)

  JVS(35) = 0.372*B(298)

  JVS(36) = 0.247*B(395)

  JVS(37) = 0.372*B(299)+0.15*B(348)+0.189*B(356)+0.189*B(364)+0.119*B(388)+0.247*B(396)

  JVS(38) = B(148)+B(168)+B(190)

  JVS(39) = 0.25*B(142)+B(146)+B(149)+B(152)

  JVS(40) = B(147)+B(166)+B(188)

  JVS(41) = 0.25*B(143)+0.25*B(162)+0.25*B(184)

  JVS(42) = B(153)+B(172)+2*B(194)

  JVS(43) = 0.25*B(185)+B(189)+B(191)+2*B(195)

  JVS(44) = 0.25*B(163)+B(167)+B(169)+B(173)

  JVS(45) = 0

  JVS(46) = 0.5*B(401)

  JVS(47) = 0.135*B(405)

  JVS(48) = 0.5*B(402)+0.135*B(406)

  JVS(49) = 0

  JVS(50) = 0.75*B(124)

  JVS(51) = 0.75*B(125)

  JVS(52) = 0

  JVS(53) = 0.75*B(142)

  JVS(54) = 0.75*B(143)+0.75*B(162)+0.75*B(184)

  JVS(55) = 0.75*B(185)

  JVS(56) = 0.75*B(163)

  JVS(57) = 0

  JVS(58) = 2*B(211)

  JVS(59) = B(407)

  JVS(60) = 2*B(212)

  JVS(61) = B(408)

  JVS(62) = 0

  JVS(63) = 6*B(211)

  JVS(64) = 7*B(277)

  JVS(65) = 0.048*B(403)+0.07*B(405)+2.693*B(407)+0.55*B(409)

  JVS(66) = 0.55*B(410)

  JVS(67) = 0.07*B(406)

  JVS(68) = 6*B(212)

  JVS(69) = 0.048*B(404)

  JVS(70) = 2.693*B(408)

  JVS(71) = 0

  JVS(72) = B(95)+B(99)

  JVS(73) = B(86)+B(90)

  JVS(74) = B(78)+B(82)

  JVS(75) = B(106)+B(112)

  JVS(76) = B(79)+B(87)+B(96)+B(107)

  JVS(77) = B(83)+B(91)+B(100)+B(113)

  JVS(78) = 0

  JVS(79) = B(97)+B(101)+B(103)+B(105)+B(116)

  JVS(80) = B(88)+B(92)+B(94)+B(104)+B(114)

  JVS(81) = B(80)+B(84)+B(85)+B(93)+B(102)+B(110)

  JVS(82) = B(81)+B(89)+B(98)+B(108)

  JVS(83) = B(109)+B(111)+B(115)+B(117)+B(118)

  JVS(84) = 0

  JVS(85) = B(519)

  JVS(86) = B(445)+B(463)

  JVS(87) = B(446)

  JVS(88) = B(520)

  JVS(89) = B(464)

  JVS(90) = 0

  JVS(91) = B(342)

  JVS(92) = B(343)

  JVS(93) = 0

  JVS(94) = B(346)

  JVS(95) = 0

  JVS(96) = B(495)

  JVS(97) = B(417)

  JVS(98) = B(418)+B(496)

  JVS(99) = 0

  JVS(100) = B(507)

  JVS(101) = B(437)+B(457)

  JVS(102) = B(438)

  JVS(103) = B(508)

  JVS(104) = B(458)

  JVS(105) = 0

  JVS(106) = B(344)

  JVS(107) = B(338)

  JVS(108) = B(340)

  JVS(109) = B(335)

  JVS(110) = B(339)+B(341)

  JVS(111) = B(345)

  JVS(112) = -B(342)-B(344)-B(346)

  JVS(113) = 0.12*B(336)

  JVS(114) = -B(343)

  JVS(115) = 0.12*B(337)

  JVS(116) = -B(345)

  JVS(117) = -B(338)

  JVS(118) = 0.75*B(336)

  JVS(119) = 0.75*B(337)-B(339)

  JVS(120) = -B(340)

  JVS(121) = -B(341)

  JVS(122) = 0

  JVS(123) = 0.044*B(469)

  JVS(124) = B(483)

  JVS(125) = 0.044*B(470)+B(484)

  JVS(126) = 0

  JVS(127) = 0.005757*B(523)

  JVS(128) = 0.005757*B(524)

  JVS(129) = 0

  JVS(130) = 0.074142*B(523)

  JVS(131) = 0.074142*B(524)

  JVS(132) = 0

  JVS(133) = 0.132523*B(523)

  JVS(134) = 0.132523*B(524)

  JVS(135) = 0

  JVS(136) = 0.478008*B(523)

  JVS(137) = 0.478008*B(524)

  JVS(138) = -B(523)

  JVS(139) = -B(524)

  JVS(140) = 0

  JVS(141) = 0.081486*B(525)

  JVS(142) = 0.081486*B(526)

  JVS(143) = 0

  JVS(144) = 0.5*B(525)

  JVS(145) = 0.5*B(526)

  JVS(146) = 0

  JVS(147) = 0.312248*B(525)

  JVS(148) = 0.312248*B(526)

  JVS(149) = -B(525)

  JVS(150) = -B(526)

  JVS(151) = -B(32)-B(34)

  JVS(152) = B(31)

  JVS(153) = -B(319)

  JVS(154) = -B(320)

  JVS(155) = -B(336)

  JVS(156) = 0.88*B(331)

  JVS(157) = 0.88*B(332)

  JVS(158) = -B(337)

  JVS(159) = 0.13*B(336)

  JVS(160) = -B(331)-B(333)-2*B(335)

  JVS(161) = B(329)

  JVS(162) = -B(332)

  JVS(163) = B(330)+0.13*B(337)

  JVS(164) = -B(334)

  JVS(165) = -B(369)

  JVS(166) = -B(370)

  JVS(167) = -B(469)

  JVS(168) = -B(470)

  JVS(169) = 0.071*B(469)

  JVS(170) = -B(483)-B(485)

  JVS(171) = B(479)

  JVS(172) = 0.071*B(470)+B(480)-B(484)-B(486)

  JVS(173) = 0.41*B(469)

  JVS(174) = -B(479)-B(481)

  JVS(175) = B(475)

  JVS(176) = 0.41*B(470)+B(476)-B(480)-B(482)

  JVS(177) = 0.3*B(469)

  JVS(178) = -B(475)-B(477)

  JVS(179) = B(471)

  JVS(180) = 0.3*B(470)+B(472)-B(476)-B(478)

  JVS(181) = B(485)

  JVS(182) = B(481)

  JVS(183) = B(477)

  JVS(184) = -B(471)

  JVS(185) = -B(472)+B(478)+B(482)+B(486)

  JVS(186) = -B(121)

  JVS(187) = B(119)

  JVS(188) = B(120)

  JVS(189) = -B(139)

  JVS(190) = B(137)

  JVS(191) = B(138)

  JVS(192) = -B(159)

  JVS(193) = B(157)

  JVS(194) = B(158)

  JVS(195) = -B(181)

  JVS(196) = B(179)

  JVS(197) = B(180)

  JVS(198) = -B(495)-B(497)

  JVS(199) = B(491)

  JVS(200) = B(419)+B(451)

  JVS(201) = B(420)+B(492)-B(496)-B(498)

  JVS(202) = B(452)

  JVS(203) = B(497)

  JVS(204) = -B(487)

  JVS(205) = B(493)

  JVS(206) = B(423)+B(455)

  JVS(207) = B(424)-B(488)+B(494)+B(498)

  JVS(208) = B(456)

  JVS(209) = B(487)

  JVS(210) = -B(491)-B(493)

  JVS(211) = B(421)+B(435)+B(453)

  JVS(212) = B(436)

  JVS(213) = B(422)+B(488)-B(492)-B(494)

  JVS(214) = B(454)

  JVS(215) = -B(69)-B(70)

  JVS(216) = B(63)+B(64)

  JVS(217) = -B(71)

  JVS(218) = -B(519)-B(521)

  JVS(219) = B(515)

  JVS(220) = B(431)+B(447)+B(465)

  JVS(221) = B(448)

  JVS(222) = B(432)+B(516)-B(520)-B(522)

  JVS(223) = B(466)

  JVS(224) = -B(515)-B(517)

  JVS(225) = B(511)

  JVS(226) = B(433)+B(449)

  JVS(227) = B(450)

  JVS(228) = B(434)+B(512)-B(516)-B(518)

  JVS(229) = B(521)

  JVS(230) = B(517)

  JVS(231) = -B(511)

  JVS(232) = B(467)

  JVS(233) = 0

  JVS(234) = -B(512)+B(518)+B(522)

  JVS(235) = B(468)

  JVS(236) = -B(507)-B(509)

  JVS(237) = B(503)

  JVS(238) = B(425)+B(439)+B(459)

  JVS(239) = B(440)

  JVS(240) = B(426)+B(504)-B(508)-B(510)

  JVS(241) = B(460)

  JVS(242) = -B(503)-B(505)

  JVS(243) = B(499)

  JVS(244) = B(427)+B(441)

  JVS(245) = B(442)

  JVS(246) = B(428)+B(500)-B(504)-B(506)

  JVS(247) = B(509)

  JVS(248) = B(505)

  JVS(249) = -B(499)

  JVS(250) = B(429)+B(443)+B(461)

  JVS(251) = B(444)

  JVS(252) = B(430)-B(500)+B(506)+B(510)

  JVS(253) = B(462)

  JVS(254) = -B(245)

  JVS(255) = -B(246)

  JVS(256) = -B(371)

  JVS(257) = -B(372)

  JVS(258) = -B(264)

  JVS(259) = 0.087*B(383)

  JVS(260) = 0.031*B(355)

  JVS(261) = 0.031*B(363)

  JVS(262) = 0.031*B(356)+0.031*B(364)

  JVS(263) = 0.087*B(384)

  JVS(264) = -B(23)-B(24)

  JVS(265) = B(21)

  JVS(266) = B(22)

  JVS(267) = -B(74)

  JVS(268) = B(411)+0.5*B(413)+B(415)

  JVS(269) = -B(75)+B(412)+0.5*B(414)

  JVS(270) = B(416)

  JVS(271) = -B(38)-B(39)-B(40)

  JVS(272) = B(36)-B(41)

  JVS(273) = B(37)

  JVS(274) = -B(411)-B(413)-B(415)

  JVS(275) = -B(412)-B(414)

  JVS(276) = -B(416)

  JVS(277) = -B(379)

  JVS(278) = -B(380)

  JVS(279) = -B(375)

  JVS(280) = -B(376)

  JVS(281) = 0.236*B(375)

  JVS(282) = -B(203)-B(205)

  JVS(283) = -B(204)

  JVS(284) = 0.236*B(376)

  JVS(285) = -B(247)-B(249)

  JVS(286) = B(80)

  JVS(287) = B(81)

  JVS(288) = -B(248)

  JVS(289) = -B(222)-B(223)

  JVS(290) = B(220)

  JVS(291) = B(221)

  JVS(292) = -B(224)

  JVS(293) = -B(211)-B(213)-B(215)

  JVS(294) = B(273)

  JVS(295) = -B(212)

  JVS(296) = -B(214)

  JVS(297) = B(274)

  JVS(298) = -B(383)

  JVS(299) = -B(384)

  JVS(300) = -B(243)

  JVS(301) = 0.25*B(92)

  JVS(302) = B(84)+0.25*B(93)+0.25*B(110)

  JVS(303) = 0.25*B(111)

  JVS(304) = -B(244)

  JVS(305) = -B(57)-B(58)-B(59)

  JVS(306) = B(55)

  JVS(307) = B(56)

  JVS(308) = -B(60)

  JVS(309) = -B(381)

  JVS(310) = -B(382)

  JVS(311) = 0.099*B(383)

  JVS(312) = 0.108*B(381)

  JVS(313) = -B(313)-B(315)

  JVS(314) = -B(314)+0.108*B(382)+0.099*B(384)

  JVS(315) = 0.093*B(383)

  JVS(316) = 0.051*B(381)

  JVS(317) = -B(316)-B(318)

  JVS(318) = -B(317)+0.051*B(382)+0.093*B(384)

  JVS(319) = 0.187*B(383)

  JVS(320) = 0.207*B(381)

  JVS(321) = -B(269)-B(271)

  JVS(322) = -B(270)+0.207*B(382)+0.187*B(384)

  JVS(323) = -B(272)

  JVS(324) = -B(377)

  JVS(325) = -B(378)

  JVS(326) = -B(373)-B(401)

  JVS(327) = -B(402)

  JVS(328) = -B(374)

  JVS(329) = 0.561*B(383)

  JVS(330) = 0.491*B(381)

  JVS(331) = -B(309)-B(311)

  JVS(332) = -B(312)

  JVS(333) = -B(310)+0.491*B(382)+0.561*B(384)

  JVS(334) = B(213)+B(215)

  JVS(335) = -B(273)

  JVS(336) = B(206)

  JVS(337) = B(207)

  JVS(338) = B(214)

  JVS(339) = -B(274)

  JVS(340) = -B(250)-B(252)

  JVS(341) = B(88)

  JVS(342) = B(89)+B(108)

  JVS(343) = B(109)

  JVS(344) = -B(251)

  JVS(345) = 0.05*B(383)

  JVS(346) = 0.059*B(381)

  JVS(347) = -B(275)-B(277)-B(278)

  JVS(348) = 0.061*B(393)+0.042*B(395)+0.015*B(397)

  JVS(349) = 0.042*B(396)

  JVS(350) = -B(276)+0.059*B(382)+0.05*B(384)+0.061*B(394)

  JVS(351) = -B(279)+0.015*B(398)

  JVS(352) = 0.017*B(381)

  JVS(353) = -B(265)-B(267)

  JVS(354) = B(208)+B(210)

  JVS(355) = B(209)

  JVS(356) = -B(266)+0.017*B(382)

  JVS(357) = -B(268)

  JVS(358) = 0.287*B(383)

  JVS(359) = 0.119*B(381)

  JVS(360) = 0.5*B(315)

  JVS(361) = 0.5*B(318)

  JVS(362) = 0.23*B(269)

  JVS(363) = -B(259)-B(260)-B(262)

  JVS(364) = 0.084*B(280)+0.9*B(282)

  JVS(365) = 0.174*B(296)+0.742*B(298)+0.008*B(300)

  JVS(366) = 0.3*B(289)+0.95*B(291)

  JVS(367) = 0.9*B(283)+0.95*B(292)+0.742*B(299)

  JVS(368) = -B(261)+0.23*B(270)+0.084*B(281)+0.3*B(290)+0.174*B(297)+0.119*B(382)+0.287*B(384)

  JVS(369) = -B(263)+0.008*B(301)

  JVS(370) = B(315)

  JVS(371) = B(318)

  JVS(372) = 0.002*B(377)

  JVS(373) = 0.393*B(373)+1.5*B(401)

  JVS(374) = B(309)+1.5*B(311)

  JVS(375) = B(259)+B(260)+B(262)

  JVS(376) = -B(49)

  JVS(377) = 0.5*B(323)+0.491*B(327)

  JVS(378) = 0.51*B(405)

  JVS(379) = 2*B(253)+B(254)+1.26*B(255)+1.26*B(257)

  JVS(380) = 0.345*B(387)

  JVS(381) = 0.157*B(355)

  JVS(382) = 0.157*B(363)

  JVS(383) = 0.275*B(347)

  JVS(384) = 0.416*B(280)+0.45*B(282)+0.5*B(284)+0.67*B(288)

  JVS(385) = 0.336*B(296)+0.498*B(298)+0.572*B(300)+1.233*B(302)

  JVS(386) = B(229)

  JVS(387) = 0.265*B(395)+0.012*B(399)

  JVS(388) = 0.475*B(291)+0.7*B(295)

  JVS(389) = B(216)+B(217)+B(218)+B(225)

  JVS(390) = 0.491*B(328)+0.012*B(400)

  JVS(391) = 0.034*B(232)+B(234)

  JVS(392) = 0.45*B(283)+0.475*B(292)+0.498*B(299)+1.5*B(312)+0.5*B(324)+0.275*B(348)+0.157*B(356)+0.157*B(364)+0.345&
               &*B(388)+0.265*B(396)+1.5*B(402)+0.51*B(406)

  JVS(393) = -B(50)+B(219)+0.034*B(233)+1.26*B(256)+B(261)+0.416*B(281)+0.336*B(297)+B(310)+0.393*B(374)+0.002*B(378)

  JVS(394) = B(226)+1.26*B(258)+B(263)+0.5*B(285)+0.572*B(301)

  JVS(395) = 2*B(24)

  JVS(396) = B(415)

  JVS(397) = B(271)

  JVS(398) = B(273)

  JVS(399) = B(278)

  JVS(400) = B(267)

  JVS(401) = B(262)

  JVS(402) = -B(46)-B(48)

  JVS(403) = B(257)

  JVS(404) = 0

  JVS(405) = 0.5*B(284)

  JVS(406) = 0.15*B(300)

  JVS(407) = B(230)

  JVS(408) = 0

  JVS(409) = 0

  JVS(410) = B(225)

  JVS(411) = B(235)

  JVS(412) = 0

  JVS(413) = B(42)

  JVS(414) = 0.2*B(66)

  JVS(415) = B(43)-B(47)

  JVS(416) = 0.2*B(67)+B(226)+B(231)+B(236)+B(258)+B(263)+B(268)+B(272)+B(274)+B(279)+0.5*B(285)+0.15*B(301)+B(416)

  JVS(417) = -B(321)-B(323)-B(325)-B(327)

  JVS(418) = -B(328)

  JVS(419) = -B(324)

  JVS(420) = -B(322)

  JVS(421) = -B(326)

  JVS(422) = 0.704*B(371)

  JVS(423) = 0.072*B(379)

  JVS(424) = 0.024*B(375)

  JVS(425) = B(205)

  JVS(426) = 0.452*B(377)

  JVS(427) = -B(237)-B(239)

  JVS(428) = 0.005*B(385)+0.001*B(387)+0.024*B(389)

  JVS(429) = 0.13*B(355)

  JVS(430) = 0.13*B(363)

  JVS(431) = 0.127*B(393)+0.045*B(395)+0.102*B(397)

  JVS(432) = 0.006*B(306)+0.02*B(308)

  JVS(433) = 0.13*B(356)+0.13*B(364)+0.001*B(388)+0.045*B(396)

  JVS(434) = 0

  JVS(435) = -B(238)+0.006*B(307)+0.704*B(372)+0.024*B(376)+0.452*B(378)+0.072*B(380)+0.005*B(386)+0.127*B(394)

  JVS(436) = 0.024*B(390)+0.102*B(398)

  JVS(437) = -B(403)-B(405)-B(407)-B(409)

  JVS(438) = -B(410)

  JVS(439) = -B(406)

  JVS(440) = -B(404)

  JVS(441) = -B(408)

  JVS(442) = 0.097*B(383)

  JVS(443) = 0.118*B(381)

  JVS(444) = 0.5*B(315)

  JVS(445) = 0.5*B(318)

  JVS(446) = 0.607*B(373)

  JVS(447) = B(311)

  JVS(448) = 0.23*B(265)

  JVS(449) = 0.009*B(327)

  JVS(450) = -B(253)-B(254)-B(255)-B(257)

  JVS(451) = 0

  JVS(452) = 0.001*B(355)

  JVS(453) = 0.001*B(363)

  JVS(454) = 0.15*B(296)+0.023*B(298)

  JVS(455) = 0.009*B(328)

  JVS(456) = 0.023*B(299)+B(312)+0.001*B(356)+0.001*B(364)

  JVS(457) = 0

  JVS(458) = -B(256)+0.23*B(266)+0.15*B(297)+0.607*B(374)+0.118*B(382)+0.097*B(384)

  JVS(459) = -B(258)

  JVS(460) = 0.24*B(269)+B(271)

  JVS(461) = 0.24*B(265)+B(267)

  JVS(462) = -B(206)-B(208)-B(210)

  JVS(463) = B(176)

  JVS(464) = -B(207)

  JVS(465) = B(174)

  JVS(466) = -B(209)

  JVS(467) = 0.24*B(266)+0.24*B(270)

  JVS(468) = B(160)

  JVS(469) = B(200)

  JVS(470) = B(161)+B(164)+B(175)+B(177)+2*B(178)+B(201)

  JVS(471) = B(165)+B(268)+B(272)

  JVS(472) = -B(385)-B(387)-B(389)-B(391)

  JVS(473) = -B(392)

  JVS(474) = -B(388)

  JVS(475) = -B(386)

  JVS(476) = -B(390)

  JVS(477) = -B(353)-B(355)-B(357)-B(359)

  JVS(478) = -B(360)

  JVS(479) = -B(356)

  JVS(480) = -B(354)

  JVS(481) = -B(358)

  JVS(482) = 0.948*B(379)

  JVS(483) = 0.559*B(375)

  JVS(484) = B(313)+B(315)

  JVS(485) = B(316)+B(318)

  JVS(486) = 0.936*B(377)

  JVS(487) = B(237)

  JVS(488) = 0.205*B(385)+0.488*B(389)

  JVS(489) = 0.5*B(353)+0.729*B(355)+0.75*B(357)

  JVS(490) = -B(95)-B(97)-B(99)-B(101)-B(103)-B(116)-B(132)-B(150)-B(170)-B(192)

  JVS(491) = 0.5*B(361)+0.729*B(363)+0.75*B(365)

  JVS(492) = 0.079*B(329)+0.126*B(347)+0.187*B(349)+0.24*B(351)

  JVS(493) = 0.001*B(393)+0.137*B(395)+0.711*B(397)

  JVS(494) = 0.675*B(289)

  JVS(495) = 0.596*B(306)+0.152*B(308)

  JVS(496) = 0.24*B(352)

  JVS(497) = 0.616*B(240)

  JVS(498) = 0.515*B(305)

  JVS(499) = 0.126*B(348)+0.729*B(356)+0.729*B(364)+0.137*B(396)

  JVS(500) = -B(104)

  JVS(501) = -B(151)+B(176)

  JVS(502) = 0

  JVS(503) = -B(102)

  JVS(504) = -B(133)+B(174)

  JVS(505) = -B(98)

  JVS(506) = -B(117)

  JVS(507) = B(238)+0.616*B(241)+0.675*B(290)+0.596*B(307)+B(314)+B(317)+0.079*B(330)+0.5*B(354)+0.5*B(362)+0.559*B(376)&
               &+0.936*B(378)+0.948*B(380)+0.205*B(386)+0.001*B(394)

  JVS(508) = -B(96)+B(160)

  JVS(509) = -B(193)+B(200)

  JVS(510) = B(161)+B(164)-B(171)+B(175)+B(177)+2*B(178)+B(201)

  JVS(511) = -B(100)+B(165)+0.187*B(350)+0.75*B(358)+0.75*B(366)+0.488*B(390)+0.711*B(398)

  JVS(512) = -B(361)-B(363)-B(365)-B(367)

  JVS(513) = -B(368)

  JVS(514) = -B(364)

  JVS(515) = -B(362)

  JVS(516) = -B(366)

  JVS(517) = -B(329)-B(347)-B(349)-B(351)

  JVS(518) = -B(352)

  JVS(519) = -B(348)

  JVS(520) = -B(330)

  JVS(521) = -B(350)

  JVS(522) = 0.23*B(329)+0.39*B(347)

  JVS(523) = -B(280)-B(282)-B(284)-B(286)-B(288)

  JVS(524) = 0.025*B(393)+0.026*B(395)+0.012*B(399)

  JVS(525) = -B(287)+0.012*B(400)

  JVS(526) = -B(283)+0.39*B(348)+0.026*B(396)

  JVS(527) = -B(281)+0.23*B(330)+0.025*B(394)

  JVS(528) = -B(285)

  JVS(529) = 0.357*B(329)+0.936*B(349)

  JVS(530) = -B(296)-B(298)-B(300)-B(302)

  JVS(531) = 0.025*B(393)

  JVS(532) = 0

  JVS(533) = -B(299)

  JVS(534) = -B(297)+0.357*B(330)+0.025*B(394)

  JVS(535) = -B(301)+0.936*B(350)

  JVS(536) = B(369)

  JVS(537) = 0.96*B(245)

  JVS(538) = 0.099*B(379)

  JVS(539) = 0.445*B(375)

  JVS(540) = 0.455*B(377)

  JVS(541) = 0.195*B(321)+0.25*B(327)

  JVS(542) = 0.984*B(403)+0.5*B(405)

  JVS(543) = 0.294*B(385)+0.154*B(387)+0.009*B(389)

  JVS(544) = 0.129*B(296)+0.047*B(298)+0.467*B(302)

  JVS(545) = -B(227)-B(229)-B(230)

  JVS(546) = 0.732*B(393)+0.456*B(395)+0.507*B(397)

  JVS(547) = 0.439*B(306)+0.431*B(308)

  JVS(548) = 0.25*B(328)

  JVS(549) = 0.034*B(232)+B(234)

  JVS(550) = 0.482*B(240)+B(242)

  JVS(551) = 0.084*B(303)+0.246*B(305)

  JVS(552) = 0.047*B(299)+0.154*B(388)+0.456*B(396)+0.5*B(406)

  JVS(553) = B(140)+B(144)+B(154)+2*B(156)+B(176)+B(198)

  JVS(554) = B(155)

  JVS(555) = -B(228)+0.034*B(233)+0.482*B(241)+0.96*B(246)+0.129*B(297)+0.084*B(304)+0.439*B(307)+0.195*B(322)+B(370)&
               &+0.445*B(376)+0.455*B(378)+0.099*B(380)+0.294*B(386)+0.732*B(394)+0.984*B(404)

  JVS(556) = B(141)

  JVS(557) = B(199)

  JVS(558) = B(177)

  JVS(559) = B(145)-B(231)+0.009*B(390)+0.507*B(398)

  JVS(560) = -B(393)-B(395)-B(397)-B(399)

  JVS(561) = -B(400)

  JVS(562) = -B(396)

  JVS(563) = -B(394)

  JVS(564) = -B(398)

  JVS(565) = 0.32*B(329)+0.16*B(347)

  JVS(566) = 0.019*B(395)+0.048*B(397)

  JVS(567) = -B(289)-B(291)-B(293)-B(295)

  JVS(568) = -B(294)

  JVS(569) = -B(292)+0.16*B(348)+0.019*B(396)

  JVS(570) = -B(290)+0.32*B(330)

  JVS(571) = 0.048*B(398)

  JVS(572) = 0.081*B(245)

  JVS(573) = 0.026*B(379)

  JVS(574) = 0.026*B(375)

  JVS(575) = 0.35*B(247)+B(249)

  JVS(576) = B(222)

  JVS(577) = B(243)

  JVS(578) = 0.024*B(377)

  JVS(579) = 0.096*B(373)

  JVS(580) = 1.61*B(321)+B(323)+0.191*B(327)

  JVS(581) = B(237)

  JVS(582) = 0.984*B(403)+0.5*B(405)

  JVS(583) = B(254)

  JVS(584) = 0

  JVS(585) = 0.732*B(385)+0.5*B(387)

  JVS(586) = 0.276*B(353)+0.235*B(355)

  JVS(587) = 0.276*B(361)+0.235*B(363)

  JVS(588) = 0.624*B(329)+0.592*B(347)+0.24*B(351)

  JVS(589) = 0.084*B(280)+0.2*B(282)+0.67*B(288)

  JVS(590) = 0.055*B(296)+0.125*B(298)+0.227*B(300)+0.3*B(302)

  JVS(591) = 0.244*B(393)+0.269*B(395)+0.079*B(397)

  JVS(592) = 0.3*B(289)+0.1*B(291)

  JVS(593) = -B(216)-B(217)-B(218)-B(220)-B(225)

  JVS(594) = 0.01*B(306)+0.134*B(308)

  JVS(595) = 0.191*B(328)+0.24*B(352)

  JVS(596) = 0.115*B(240)

  JVS(597) = 0.213*B(303)+0.506*B(305)

  JVS(598) = 0.2*B(283)+0.1*B(292)+0.125*B(299)+B(324)+0.592*B(348)+0.235*B(356)+0.235*B(364)+0.5*B(388)+0.269*B(396)&
               &+0.5*B(406)

  JVS(599) = 0.75*B(92)

  JVS(600) = B(146)+B(198)

  JVS(601) = 0

  JVS(602) = B(78)+B(82)+B(84)+2*B(85)+0.75*B(93)+0.75*B(110)+B(128)+B(147)+B(166)+B(188)

  JVS(603) = B(129)+B(196)

  JVS(604) = -B(221)

  JVS(605) = 0.75*B(111)

  JVS(606) = -B(219)+B(238)+0.115*B(241)+B(244)+0.081*B(246)+0.35*B(248)+0.084*B(281)+0.3*B(290)+0.055*B(297)+0.213&
               &*B(304)+0.01*B(307)+1.61*B(322)+0.624*B(330)+0.276*B(354)+0.276*B(362)+0.096*B(374)+0.026*B(376)+0.024&
               &*B(378)+0.026*B(380)+0.732*B(386)+0.244*B(394)+0.984*B(404)

  JVS(607) = B(79)+B(182)

  JVS(608) = B(183)+B(186)+B(189)+B(197)+B(199)+B(200)+2*B(202)

  JVS(609) = B(167)+B(201)

  JVS(610) = B(83)+B(187)-B(226)+0.227*B(301)+0.079*B(398)

  JVS(611) = B(203)

  JVS(612) = 0.511*B(389)

  JVS(613) = 0.276*B(357)

  JVS(614) = 0.276*B(365)

  JVS(615) = 0.572*B(300)

  JVS(616) = 0.321*B(397)

  JVS(617) = -0.69*B(306)-B(308)

  JVS(618) = 0

  JVS(619) = 0

  JVS(620) = B(204)

  JVS(621) = B(106)

  JVS(622) = -0.69*B(307)

  JVS(623) = B(107)

  JVS(624) = 0.572*B(301)+0.276*B(358)+0.276*B(366)+0.511*B(390)+0.321*B(398)

  JVS(625) = B(34)

  JVS(626) = -B(327)

  JVS(627) = -B(409)

  JVS(628) = -B(391)

  JVS(629) = -B(359)

  JVS(630) = -B(367)

  JVS(631) = -B(351)

  JVS(632) = -B(286)

  JVS(633) = -B(399)

  JVS(634) = -B(293)

  JVS(635) = -B(2)-B(4)-B(6)-B(9)-B(11)-B(287)-B(294)-B(328)-B(352)-B(360)-B(368)-B(392)-B(400)-B(410)

  JVS(636) = -B(5)+B(30)

  JVS(637) = B(1)-B(10)-B(12)

  JVS(638) = 0

  JVS(639) = -B(7)

  JVS(640) = B(29)

  JVS(641) = 0.261*B(371)

  JVS(642) = 0.204*B(379)

  JVS(643) = 0.122*B(375)

  JVS(644) = B(313)

  JVS(645) = B(316)

  JVS(646) = 0.244*B(377)

  JVS(647) = B(309)

  JVS(648) = B(250)+B(252)

  JVS(649) = B(325)

  JVS(650) = 0.45*B(409)

  JVS(651) = 0.497*B(385)+0.363*B(387)+0.037*B(389)+0.45*B(391)

  JVS(652) = 0.474*B(353)+0.205*B(355)+0.474*B(357)+0.147*B(359)

  JVS(653) = 0.474*B(361)+0.205*B(363)+0.474*B(365)+0.147*B(367)

  JVS(654) = B(286)

  JVS(655) = 0.013*B(296)+0.218*B(300)

  JVS(656) = 0.511*B(393)+0.305*B(395)+0.151*B(397)+0.069*B(399)

  JVS(657) = 0.675*B(289)+0.45*B(293)

  JVS(658) = 0.213*B(306)+0.147*B(308)

  JVS(659) = B(287)+0.45*B(294)+0.147*B(360)+0.147*B(368)+0.45*B(392)+0.069*B(400)+0.45*B(410)

  JVS(660) = -B(232)-B(234)-B(235)

  JVS(661) = 0.37*B(240)

  JVS(662) = 0.558*B(303)+0.71*B(305)

  JVS(663) = 0.205*B(356)+0.205*B(364)+0.363*B(388)+0.305*B(396)

  JVS(664) = 0

  JVS(665) = 0

  JVS(666) = 0

  JVS(667) = 0

  JVS(668) = -B(233)+0.37*B(241)+B(251)+0.675*B(290)+0.013*B(297)+0.558*B(304)+0.213*B(307)+B(310)+B(314)+B(317)+0.474&
               &*B(354)+0.474*B(362)+0.261*B(372)+0.122*B(376)+0.244*B(378)+0.204*B(380)+0.497*B(386)+0.511*B(394)

  JVS(669) = 0

  JVS(670) = -B(236)+0.218*B(301)+B(326)+0.474*B(358)+0.474*B(366)+0.037*B(390)+0.151*B(398)

  JVS(671) = 0.089*B(379)

  JVS(672) = 0.332*B(375)

  JVS(673) = 0.11*B(377)

  JVS(674) = 0.55*B(409)

  JVS(675) = 0.437*B(391)

  JVS(676) = 0.416*B(280)

  JVS(677) = 0.15*B(296)+0.21*B(298)+0.233*B(302)

  JVS(678) = 0.072*B(393)+0.026*B(395)+0.001*B(397)+0.659*B(399)

  JVS(679) = 0.55*B(293)

  JVS(680) = 0.177*B(306)+0.243*B(308)

  JVS(681) = 0.55*B(294)+0.437*B(392)+0.659*B(400)+0.55*B(410)

  JVS(682) = -B(240)-B(242)

  JVS(683) = 0.115*B(303)

  JVS(684) = 0.21*B(299)+0.026*B(396)

  JVS(685) = 0.5*B(114)

  JVS(686) = 0

  JVS(687) = 0.5*B(110)

  JVS(688) = 0.5*B(111)+B(112)+0.5*B(115)+B(118)

  JVS(689) = -B(241)+0.416*B(281)+0.15*B(297)+0.115*B(304)+0.177*B(307)+0.332*B(376)+0.11*B(378)+0.089*B(380)+0.072&
               &*B(394)

  JVS(690) = 0

  JVS(691) = B(113)+0.001*B(398)

  JVS(692) = 0.417*B(379)

  JVS(693) = 0.055*B(381)

  JVS(694) = 0.125*B(377)

  JVS(695) = 0.119*B(385)+0.215*B(387)+0.113*B(391)

  JVS(696) = 0.276*B(353)+0.276*B(355)+0.853*B(359)

  JVS(697) = 0.276*B(361)+0.276*B(363)+0.853*B(367)

  JVS(698) = 0.1*B(347)+0.75*B(351)

  JVS(699) = 0.332*B(296)

  JVS(700) = 0.043*B(395)+0.259*B(399)

  JVS(701) = 0.7*B(295)

  JVS(702) = 0.048*B(306)+0.435*B(308)

  JVS(703) = 0.75*B(352)+0.853*B(360)+0.853*B(368)+0.113*B(392)+0.259*B(400)

  JVS(704) = -0.671*B(303)-B(305)

  JVS(705) = 0.1*B(348)+0.276*B(356)+0.276*B(364)+0.215*B(388)+0.043*B(396)

  JVS(706) = 0.5*B(114)

  JVS(707) = B(152)

  JVS(708) = 0

  JVS(709) = 0.5*B(110)

  JVS(710) = B(134)

  JVS(711) = 0.5*B(111)+0.5*B(115)+B(118)+B(135)+B(153)+B(172)

  JVS(712) = 0.332*B(297)-0.671*B(304)+0.048*B(307)+0.276*B(354)+0.276*B(362)+0.125*B(378)+0.417*B(380)+0.055*B(382)&
               &+0.119*B(386)

  JVS(713) = 0

  JVS(714) = B(173)

  JVS(715) = 0

  JVS(716) = -B(401)

  JVS(717) = -B(311)

  JVS(718) = -B(323)

  JVS(719) = -B(405)

  JVS(720) = -B(387)

  JVS(721) = -B(355)

  JVS(722) = -B(363)

  JVS(723) = -B(347)

  JVS(724) = -B(282)

  JVS(725) = -B(298)

  JVS(726) = -B(395)

  JVS(727) = -B(291)

  JVS(728) = B(2)-B(4)

  JVS(729) = -B(5)-B(13)-B(15)-B(30)-B(31)-B(51)-B(61)-B(283)-B(292)-B(299)-B(312)-B(324)-B(348)-B(356)-B(364)-B(388)&
               &-B(396)-B(402)-B(406)

  JVS(730) = 0.25*B(142)

  JVS(731) = -B(16)

  JVS(732) = 0.25*B(124)

  JVS(733) = -B(62)+0.25*B(125)+0.25*B(143)+0.25*B(162)+0.25*B(184)

  JVS(734) = -B(52)

  JVS(735) = -B(14)

  JVS(736) = 0.25*B(185)

  JVS(737) = 0.25*B(163)

  JVS(738) = 0

  JVS(739) = B(369)

  JVS(740) = 0.05*B(245)

  JVS(741) = 0.965*B(371)

  JVS(742) = 0.653*B(379)

  JVS(743) = 0.695*B(375)

  JVS(744) = 0.804*B(383)

  JVS(745) = 0.765*B(381)

  JVS(746) = B(315)

  JVS(747) = B(318)

  JVS(748) = 0.76*B(269)

  JVS(749) = 0.835*B(377)

  JVS(750) = 0.1*B(373)

  JVS(751) = B(309)

  JVS(752) = 0.34*B(250)

  JVS(753) = 0.76*B(265)

  JVS(754) = B(321)+B(325)+0.2*B(327)

  JVS(755) = 0.984*B(403)+0.949*B(407)

  JVS(756) = 0

  JVS(757) = 0.91*B(385)+0.022*B(387)+0.824*B(389)

  JVS(758) = 0.75*B(353)+0.031*B(355)+0.276*B(357)

  JVS(759) = 0.75*B(361)+0.031*B(363)+0.276*B(365)

  JVS(760) = 0.907*B(329)+0.066*B(347)+0.749*B(349)

  JVS(761) = 0.5*B(280)+0.1*B(282)+0.5*B(284)+0.33*B(288)

  JVS(762) = 0.67*B(296)+0.048*B(298)+0.799*B(300)

  JVS(763) = 0.918*B(393)+0.033*B(395)+0.442*B(397)+0.012*B(399)

  JVS(764) = 0.3*B(289)+0.05*B(291)

  JVS(765) = 0.376*B(306)+0.564*B(308)

  JVS(766) = 0.2*B(328)+0.012*B(400)

  JVS(767) = 0.034*B(232)+B(234)

  JVS(768) = 0.37*B(240)+B(242)

  JVS(769) = 0.473*B(303)+0.96*B(305)

  JVS(770) = 0.1*B(283)+0.05*B(292)+0.048*B(299)+0.066*B(348)+0.031*B(356)+0.031*B(364)+0.022*B(388)+0.033*B(396)

  JVS(771) = -B(86)-B(88)-B(90)-B(92)-2*B(94)-B(114)-B(130)-B(148)-B(168)-B(190)

  JVS(772) = B(140)+B(144)-B(149)+B(154)+2*B(156)+B(176)+B(198)

  JVS(773) = 0

  JVS(774) = -B(93)

  JVS(775) = -B(131)+B(155)

  JVS(776) = -B(89)

  JVS(777) = -B(115)

  JVS(778) = 0.034*B(233)+0.37*B(241)+0.05*B(246)+0.34*B(251)+0.76*B(266)+0.76*B(270)+0.5*B(281)+0.3*B(290)+0.67*B(297)&
               &+0.473*B(304)+0.376*B(307)+B(310)+B(322)+0.907*B(330)+0.75*B(354)+0.75*B(362)+B(370)+0.965*B(372)+0.1*B(374)&
               &+0.695*B(376)+0.835*B(378)+0.653*B(380)+0.765*B(382)+0.804*B(384)+0.91*B(386)+0.918*B(394)+0.984*B(404)

  JVS(779) = -B(87)+B(141)

  JVS(780) = -B(191)+B(199)

  JVS(781) = -B(169)+B(177)

  JVS(782) = -B(91)+B(145)+0.5*B(285)+0.799*B(301)+B(326)+0.749*B(350)+0.276*B(358)+0.276*B(366)+0.824*B(390)+0.442&
               &*B(398)+0.949*B(408)

  JVS(783) = B(139)

  JVS(784) = 0.37*B(255)+0.37*B(257)

  JVS(785) = 0

  JVS(786) = 0.201*B(355)

  JVS(787) = 0.201*B(363)

  JVS(788) = 0.1*B(282)

  JVS(789) = 0.048*B(298)+0.3*B(302)

  JVS(790) = 0.006*B(395)

  JVS(791) = 0.05*B(291)

  JVS(792) = 0

  JVS(793) = 0.965*B(232)+B(235)

  JVS(794) = 0.096*B(240)

  JVS(795) = 0.049*B(303)+0.333*B(305)

  JVS(796) = 0.1*B(283)+0.05*B(292)+0.048*B(299)+0.201*B(356)+0.201*B(364)+0.006*B(396)

  JVS(797) = -B(148)

  JVS(798) = -B(137)-B(140)-B(142)-B(144)-B(146)-B(149)-B(152)-B(154)-2*B(156)-B(176)-B(198)

  JVS(799) = -B(138)

  JVS(800) = -B(147)

  JVS(801) = -B(155)

  JVS(802) = -B(143)

  JVS(803) = -B(153)

  JVS(804) = 0.965*B(233)+0.096*B(241)+0.37*B(256)+0.049*B(304)

  JVS(805) = -B(141)

  JVS(806) = -B(199)

  JVS(807) = -B(177)

  JVS(808) = -B(145)+B(236)+0.37*B(258)

  JVS(809) = B(121)

  JVS(810) = B(139)

  JVS(811) = B(159)

  JVS(812) = B(181)

  JVS(813) = B(23)

  JVS(814) = B(39)+B(40)

  JVS(815) = -B(203)

  JVS(816) = B(223)

  JVS(817) = -B(211)

  JVS(818) = B(57)+0.61*B(58)+B(59)

  JVS(819) = 0

  JVS(820) = B(48)

  JVS(821) = 0

  JVS(822) = -B(206)

  JVS(823) = 0.474*B(357)

  JVS(824) = B(95)+B(99)

  JVS(825) = 0.474*B(365)

  JVS(826) = 0.187*B(349)

  JVS(827) = 0

  JVS(828) = 0

  JVS(829) = 0

  JVS(830) = 0.391*B(397)

  JVS(831) = 0

  JVS(832) = 0

  JVS(833) = 0.338*B(306)+B(308)

  JVS(834) = B(6)-B(9)-B(11)

  JVS(835) = 0

  JVS(836) = 0

  JVS(837) = 0

  JVS(838) = B(13)-B(15)

  JVS(839) = B(86)+B(90)

  JVS(840) = -B(137)+B(140)+B(144)

  JVS(841) = -B(1)-B(10)-B(12)-B(16)-B(21)-B(42)-B(55)-B(119)-B(138)-B(157)-B(179)-B(204)-B(207)-B(212)

  JVS(842) = B(78)+B(82)

  JVS(843) = -B(120)+B(122)+B(126)

  JVS(844) = B(53)-B(56)+0.8*B(66)

  JVS(845) = B(112)

  JVS(846) = B(41)-B(43)+B(44)+B(60)+0.338*B(307)

  JVS(847) = B(7)+B(14)+2*B(17)+2*B(19)+B(54)+B(79)+B(87)+B(96)+B(123)+B(141)+B(160)+B(182)+B(224)

  JVS(848) = -B(180)+B(183)+B(186)

  JVS(849) = -B(158)+B(161)+B(164)

  JVS(850) = 2*B(18)-B(22)+B(29)+B(45)+0.8*B(67)+2*B(68)+B(83)+B(91)+B(100)+B(113)+B(127)+B(145)+B(165)+B(187)+0.187&
               &*B(350)+0.474*B(358)+0.474*B(366)+0.391*B(398)

  JVS(851) = B(319)

  JVS(852) = B(205)

  JVS(853) = 0.65*B(247)

  JVS(854) = 0.011*B(377)

  JVS(855) = 0.3*B(327)

  JVS(856) = B(239)

  JVS(857) = 0.26*B(405)

  JVS(858) = 0.076*B(387)

  JVS(859) = 0

  JVS(860) = 0

  JVS(861) = 0.25*B(351)

  JVS(862) = B(229)

  JVS(863) = 0.197*B(395)+0.03*B(397)

  JVS(864) = 0.3*B(295)

  JVS(865) = 0

  JVS(866) = 0.3*B(328)+0.25*B(352)

  JVS(867) = 0

  JVS(868) = 0

  JVS(869) = 0

  JVS(870) = 0.076*B(388)+0.197*B(396)+0.26*B(406)

  JVS(871) = -B(92)

  JVS(872) = -B(146)+B(154)

  JVS(873) = 0

  JVS(874) = -B(78)-B(80)-B(82)-2*B(84)-2*B(85)-B(93)-B(110)-B(128)-B(147)-B(166)-B(188)

  JVS(875) = B(122)+B(126)-B(129)+2*B(136)+B(155)+B(174)+B(196)

  JVS(876) = -B(81)

  JVS(877) = -B(111)

  JVS(878) = 0.65*B(248)+B(320)+0.011*B(378)

  JVS(879) = -B(79)+B(123)

  JVS(880) = -B(189)+B(197)

  JVS(881) = -B(167)+B(175)

  JVS(882) = -B(83)+B(127)+0.03*B(398)

  JVS(883) = B(121)

  JVS(884) = 2*B(264)

  JVS(885) = 0

  JVS(886) = B(313)+0.5*B(315)

  JVS(887) = B(316)+0.5*B(318)

  JVS(888) = 0.011*B(377)

  JVS(889) = B(259)+B(260)+B(262)

  JVS(890) = B(237)+B(239)

  JVS(891) = 0

  JVS(892) = 0.123*B(355)

  JVS(893) = 0.123*B(363)

  JVS(894) = 0.67*B(288)

  JVS(895) = 0.467*B(302)

  JVS(896) = B(227)+B(230)

  JVS(897) = 0.137*B(395)

  JVS(898) = 0.675*B(289)

  JVS(899) = 0

  JVS(900) = 0

  JVS(901) = 0

  JVS(902) = 0.492*B(240)+B(242)

  JVS(903) = 0.029*B(303)+0.667*B(305)

  JVS(904) = 0.123*B(356)+0.123*B(364)+0.137*B(396)

  JVS(905) = -B(130)

  JVS(906) = -B(154)+B(198)

  JVS(907) = -B(119)

  JVS(908) = -B(128)

  JVS(909) = -B(120)-B(122)-B(124)-B(126)-B(129)-B(131)-B(134)-2*B(136)-B(155)-B(174)

  JVS(910) = -B(125)

  JVS(911) = -B(135)

  JVS(912) = B(228)+B(238)+0.492*B(241)+B(261)+0.675*B(290)+0.029*B(304)+B(314)+B(317)+0.011*B(378)

  JVS(913) = -B(123)+B(182)

  JVS(914) = B(183)+B(186)+B(199)+B(200)+2*B(202)

  JVS(915) = -B(175)+B(201)

  JVS(916) = -B(127)+B(187)+B(231)+B(263)

  JVS(917) = B(70)

  JVS(918) = 0.95*B(245)

  JVS(919) = B(74)

  JVS(920) = B(39)

  JVS(921) = 0.5*B(413)

  JVS(922) = B(249)

  JVS(923) = B(222)+B(223)

  JVS(924) = -B(213)

  JVS(925) = 0.187*B(383)

  JVS(926) = B(243)

  JVS(927) = B(57)+0.61*B(58)

  JVS(928) = 0.224*B(381)

  JVS(929) = 0.5*B(315)

  JVS(930) = 0.5*B(318)

  JVS(931) = 0.297*B(373)+1.5*B(401)

  JVS(932) = 1.5*B(311)

  JVS(933) = 0

  JVS(934) = B(252)

  JVS(935) = B(259)

  JVS(936) = B(49)

  JVS(937) = 0.12*B(323)+0.5*B(327)

  JVS(938) = 0.06*B(405)

  JVS(939) = 2*B(253)+0.63*B(255)+0.63*B(257)

  JVS(940) = -B(208)

  JVS(941) = 0.056*B(387)

  JVS(942) = 0.033*B(355)

  JVS(943) = 0.033*B(363)

  JVS(944) = 0.907*B(329)

  JVS(945) = 0.008*B(282)+0.34*B(288)

  JVS(946) = 0.4*B(298)+1.233*B(302)

  JVS(947) = B(229)

  JVS(948) = 0.003*B(395)+0.013*B(399)

  JVS(949) = 0.064*B(291)

  JVS(950) = 2*B(216)+B(218)-B(220)+B(225)

  JVS(951) = 0.113*B(306)+0.341*B(308)

  JVS(952) = 0.5*B(328)+0.013*B(400)

  JVS(953) = B(234)

  JVS(954) = 0

  JVS(955) = 0.379*B(303)

  JVS(956) = B(51)-B(61)+0.008*B(283)+0.064*B(292)+0.4*B(299)+1.5*B(312)+0.12*B(324)+0.033*B(356)+0.033*B(364)+0.056&
               &*B(388)+0.003*B(396)+1.5*B(402)+0.06*B(406)

  JVS(957) = B(86)-B(88)+B(90)+B(92)+B(94)+B(114)

  JVS(958) = -B(142)

  JVS(959) = -B(55)

  JVS(960) = B(78)-B(80)+B(82)+2*B(85)+B(93)+B(110)

  JVS(961) = -B(124)

  JVS(962) = -B(53)-B(56)-B(62)-2*B(63)-2*B(64)-B(66)-B(72)-B(81)-B(89)-B(108)-B(125)-B(143)-B(162)-B(184)-B(209)-B(214)&
               &-B(221)

  JVS(963) = -B(109)+B(111)+B(112)+B(115)+B(118)

  JVS(964) = B(44)+B(50)+B(52)+B(71)-B(73)+B(75)+B(76)+B(219)+B(244)+0.95*B(246)+0.63*B(256)+0.379*B(304)+0.113*B(307)&
               &+0.907*B(330)+0.297*B(374)+0.224*B(382)+0.187*B(384)+0.5*B(414)

  JVS(965) = -B(54)+B(79)+B(87)+B(224)

  JVS(966) = -B(185)

  JVS(967) = -B(163)

  JVS(968) = B(45)-B(67)+B(83)+B(91)+B(113)+B(226)+0.63*B(258)

  JVS(969) = 0.035*B(371)

  JVS(970) = 0.347*B(379)

  JVS(971) = 0.07*B(375)

  JVS(972) = 0.009*B(383)

  JVS(973) = 0.011*B(381)

  JVS(974) = 0.143*B(377)

  JVS(975) = 0.016*B(403)+0.051*B(407)

  JVS(976) = 0.09*B(385)+0.001*B(387)+0.176*B(389)

  JVS(977) = 0.25*B(353)+0.18*B(355)+0.25*B(357)

  JVS(978) = 0.25*B(361)+0.18*B(363)+0.25*B(365)

  JVS(979) = 0.093*B(329)+0.008*B(347)+0.064*B(349)+0.01*B(351)

  JVS(980) = 0.041*B(296)+0.051*B(300)

  JVS(981) = 0.082*B(393)+0.002*B(395)+0.136*B(397)+0.001*B(399)

  JVS(982) = 0.025*B(289)

  JVS(983) = 0.173*B(306)+0.095*B(308)

  JVS(984) = 0.01*B(352)+0.001*B(400)

  JVS(985) = 0.001*B(232)

  JVS(986) = 0.042*B(240)

  JVS(987) = 0.07*B(303)+0.04*B(305)

  JVS(988) = 0.008*B(348)+0.18*B(356)+0.18*B(364)+0.001*B(388)+0.002*B(396)

  JVS(989) = -B(114)

  JVS(990) = -B(152)

  JVS(991) = 0

  JVS(992) = -B(110)

  JVS(993) = -B(134)

  JVS(994) = -B(108)

  JVS(995) = -B(106)-B(109)-B(111)-B(112)-B(115)-2*B(118)-B(135)-B(153)-B(172)-B(194)

  JVS(996) = 0.001*B(233)+0.042*B(241)+0.025*B(290)+0.041*B(297)+0.07*B(304)+0.173*B(307)+0.093*B(330)+0.25*B(354)+0.25&
               &*B(362)+0.035*B(372)+0.07*B(376)+0.143*B(378)+0.347*B(380)+0.011*B(382)+0.009*B(384)+0.09*B(386)+0.082&
               &*B(394)+0.016*B(404)

  JVS(997) = -B(107)

  JVS(998) = -B(195)

  JVS(999) = -B(173)

  JVS(1000) = -B(113)+0.051*B(301)+0.064*B(350)+0.25*B(358)+0.25*B(366)+0.176*B(390)+0.136*B(398)+0.051*B(408)

  JVS(1001) = 2*B(32)

  JVS(1002) = -B(319)

  JVS(1003) = -B(369)

  JVS(1004) = 2*B(69)-B(70)

  JVS(1005) = -B(245)

  JVS(1006) = -B(371)

  JVS(1007) = -B(74)

  JVS(1008) = B(38)-B(40)

  JVS(1009) = -B(411)-B(413)

  JVS(1010) = -B(379)

  JVS(1011) = -B(375)

  JVS(1012) = -0.65*B(247)+B(249)

  JVS(1013) = -B(383)

  JVS(1014) = -B(243)

  JVS(1015) = 0.39*B(58)-B(59)

  JVS(1016) = -B(381)

  JVS(1017) = -B(313)

  JVS(1018) = -B(316)

  JVS(1019) = -B(269)

  JVS(1020) = -B(377)

  JVS(1021) = -0.397*B(373)+0.5*B(401)

  JVS(1022) = -B(309)+0.5*B(311)

  JVS(1023) = -0.34*B(250)+B(252)

  JVS(1024) = -B(275)

  JVS(1025) = -B(265)

  JVS(1026) = -B(260)

  JVS(1027) = -B(49)

  JVS(1028) = -B(46)+B(48)

  JVS(1029) = -B(321)+0.12*B(323)

  JVS(1030) = -B(237)

  JVS(1031) = -B(403)+0.32*B(405)

  JVS(1032) = -B(255)

  JVS(1033) = 0

  JVS(1034) = -B(385)+0.155*B(387)

  JVS(1035) = -B(353)+0.567*B(355)

  JVS(1036) = -B(361)+0.567*B(363)

  JVS(1037) = -B(329)+0.266*B(347)

  JVS(1038) = -B(280)+0.208*B(282)+0.33*B(288)

  JVS(1039) = -B(296)+0.285*B(298)

  JVS(1040) = -B(227)

  JVS(1041) = -B(393)+0.378*B(395)

  JVS(1042) = -B(289)+0.164*B(291)

  JVS(1043) = -B(218)

  JVS(1044) = -B(306)

  JVS(1045) = 0

  JVS(1046) = -B(232)

  JVS(1047) = -B(240)

  JVS(1048) = -B(303)

  JVS(1049) = -B(51)+B(61)+0.208*B(283)+0.164*B(292)+0.285*B(299)+0.5*B(312)+0.12*B(324)+0.266*B(348)+0.567*B(356)+0.567&
                &*B(364)+0.155*B(388)+0.378*B(396)+0.5*B(402)+0.32*B(406)

  JVS(1050) = 0

  JVS(1051) = 0

  JVS(1052) = -B(42)

  JVS(1053) = 0

  JVS(1054) = 0

  JVS(1055) = B(53)+B(62)+0.8*B(66)-B(72)

  JVS(1056) = 0

  JVS(1057) = -B(36)-B(41)-B(43)-B(44)-B(47)-B(50)-B(52)-B(60)-B(71)-B(73)-B(75)-B(76)-B(219)-B(228)-B(233)-B(238)&
                &-B(241)-B(244)-B(246)-0.65*B(248)-0.34*B(251)-B(256)-B(261)-B(266)-B(270)-B(276)-B(281)-B(290)-B(297)&
                &-B(304)-B(307)-B(310)-B(314)-B(317)-B(320)-B(322)-B(330)-B(354)-B(362)-B(370)-B(372)-0.397*B(374)-B(376)&
                &-B(378)-B(380)-B(382)-B(384)-B(386)-B(394)-B(404)-B(412)-B(414)

  JVS(1058) = -B(37)+B(54)

  JVS(1059) = 0

  JVS(1060) = 0

  JVS(1061) = -B(45)+0.8*B(67)

  JVS(1062) = B(38)

  JVS(1063) = -B(223)

  JVS(1064) = -B(95)

  JVS(1065) = 0

  JVS(1066) = 0

  JVS(1067) = 0

  JVS(1068) = 0

  JVS(1069) = 0

  JVS(1070) = 0

  JVS(1071) = -B(6)+B(9)

  JVS(1072) = 0

  JVS(1073) = 0

  JVS(1074) = -B(13)

  JVS(1075) = -B(86)

  JVS(1076) = -B(140)

  JVS(1077) = B(1)+B(10)+B(26)

  JVS(1078) = -B(78)

  JVS(1079) = -B(122)

  JVS(1080) = -B(53)

  JVS(1081) = -B(106)

  JVS(1082) = -B(36)

  JVS(1083) = -B(7)-B(14)-B(17)-2*B(19)-B(37)-B(54)-B(79)-B(87)-B(96)-B(107)-B(123)-B(141)-B(160)-B(182)-B(224)

  JVS(1084) = -B(183)

  JVS(1085) = -B(161)

  JVS(1086) = -B(18)+B(27)+B(28)

  JVS(1087) = B(181)

  JVS(1088) = 0.192*B(347)+0.24*B(351)

  JVS(1089) = 0.5*B(280)+0.5*B(284)+0.33*B(288)

  JVS(1090) = 0.289*B(296)+0.15*B(300)

  JVS(1091) = 0

  JVS(1092) = 0.3*B(295)

  JVS(1093) = 0.24*B(352)

  JVS(1094) = 0.192*B(348)

  JVS(1095) = -B(190)

  JVS(1096) = -B(198)

  JVS(1097) = -B(179)

  JVS(1098) = -B(188)

  JVS(1099) = -B(196)

  JVS(1100) = -B(184)

  JVS(1101) = -B(194)

  JVS(1102) = 0.5*B(281)+0.289*B(297)

  JVS(1103) = -B(182)

  JVS(1104) = -B(180)-B(183)-B(185)-B(186)-B(189)-B(191)-B(195)-B(197)-B(199)-B(200)-2*B(202)

  JVS(1105) = -B(201)

  JVS(1106) = -B(187)+0.5*B(285)+0.15*B(301)

  JVS(1107) = B(159)

  JVS(1108) = B(275)+B(278)

  JVS(1109) = 0

  JVS(1110) = 0

  JVS(1111) = 0

  JVS(1112) = -B(168)

  JVS(1113) = -B(176)

  JVS(1114) = -B(157)

  JVS(1115) = -B(166)

  JVS(1116) = -B(174)

  JVS(1117) = -B(162)

  JVS(1118) = -B(172)

  JVS(1119) = B(276)

  JVS(1120) = -B(160)

  JVS(1121) = -B(200)

  JVS(1122) = -B(158)-B(161)-B(163)-B(164)-B(167)-B(169)-B(173)-B(175)-B(177)-2*B(178)-B(201)

  JVS(1123) = -B(165)+B(279)

  JVS(1124) = B(23)

  JVS(1125) = -B(415)

  JVS(1126) = 0.39*B(58)

  JVS(1127) = -B(271)

  JVS(1128) = -B(273)

  JVS(1129) = -B(278)

  JVS(1130) = -B(267)

  JVS(1131) = -B(262)

  JVS(1132) = B(46)

  JVS(1133) = -B(325)

  JVS(1134) = -B(407)

  JVS(1135) = -B(257)

  JVS(1136) = 0

  JVS(1137) = -B(389)

  JVS(1138) = -B(357)

  JVS(1139) = -B(99)

  JVS(1140) = -B(365)

  JVS(1141) = -B(349)

  JVS(1142) = -B(284)

  JVS(1143) = -B(300)

  JVS(1144) = -B(230)

  JVS(1145) = -B(397)

  JVS(1146) = 0

  JVS(1147) = -B(225)

  JVS(1148) = 0

  JVS(1149) = B(11)

  JVS(1150) = -B(235)

  JVS(1151) = 0

  JVS(1152) = 0

  JVS(1153) = B(15)

  JVS(1154) = -B(90)

  JVS(1155) = -B(144)

  JVS(1156) = B(12)+B(16)-B(21)-B(26)

  JVS(1157) = -B(82)

  JVS(1158) = -B(126)

  JVS(1159) = -B(66)

  JVS(1160) = -B(112)

  JVS(1161) = -B(44)+B(47)

  JVS(1162) = -B(17)

  JVS(1163) = -B(186)

  JVS(1164) = -B(164)

  JVS(1165) = -B(18)-B(22)-B(27)-B(28)-B(29)-B(45)-B(67)-2*B(68)-B(83)-B(91)-B(100)-B(113)-B(127)-B(145)-B(165)-B(187)&
                &-B(226)-B(231)-B(236)-B(258)-B(263)-B(268)-B(272)-B(274)-B(279)-B(285)-B(301)-B(326)-B(350)-B(358)-B(366)&
                &-B(390)-B(398)-B(408)-B(416)
      
END SUBROUTINE saprc99_mosaic_4bin_iepox_vbs_aq_Jac_SP














SUBROUTINE saprc99_mosaic_4bin_iepox_vbs_aq_KppDecomp( JVS, IER )







      INTEGER  :: IER
      REAL(kind=dp) :: JVS(1165), W(118), a
      INTEGER  :: k, kk, j, jj

      a = 0. 
      IER = 0
      DO k=1,NVAR

        
        IF ( ABS(JVS(LU_DIAG(k))) < TINY(a) ) THEN
            IER = k
            RETURN
        END IF
        DO kk = LU_CROW(k), LU_CROW(k+1)-1
              W( LU_ICOL(kk) ) = JVS(kk)
        END DO
        DO kk = LU_CROW(k), LU_DIAG(k)-1
            j = LU_ICOL(kk)
            a = -W(j) / JVS( LU_DIAG(j) )
            W(j) = -a
            DO jj = LU_DIAG(j)+1, LU_CROW(j+1)-1
               W( LU_ICOL(jj) ) = W( LU_ICOL(jj) ) + a*JVS(jj)
            END DO
         END DO
         DO kk = LU_CROW(k), LU_CROW(k+1)-1
            JVS(kk) = W( LU_ICOL(kk) )
         END DO
      END DO
      
END SUBROUTINE saprc99_mosaic_4bin_iepox_vbs_aq_KppDecomp



SUBROUTINE saprc99_mosaic_4bin_iepox_vbs_aq_KppDecompCmplx( JVS, IER )







      INTEGER  :: IER
      DOUBLE COMPLEX :: JVS(1165), W(118), a
      INTEGER  :: k, kk, j, jj

      IER = 0
      DO k=1,NVAR
        IF ( JVS( LU_DIAG(k) ) .EQ. 0. ) THEN
            IER = k
            RETURN
        END IF
        DO kk = LU_CROW(k), LU_CROW(k+1)-1
              W( LU_ICOL(kk) ) = JVS(kk)
        END DO
        DO kk = LU_CROW(k), LU_DIAG(k)-1
            j = LU_ICOL(kk)
            a = -W(j) / JVS( LU_DIAG(j) )
            W(j) = -a
            DO jj = LU_DIAG(j)+1, LU_CROW(j+1)-1
               W( LU_ICOL(jj) ) = W( LU_ICOL(jj) ) + a*JVS(jj)
            END DO
         END DO
         DO kk = LU_CROW(k), LU_CROW(k+1)-1
            JVS(kk) = W( LU_ICOL(kk) )
         END DO
      END DO
      
END SUBROUTINE saprc99_mosaic_4bin_iepox_vbs_aq_KppDecompCmplx


SUBROUTINE saprc99_mosaic_4bin_iepox_vbs_aq_KppSolveIndirect( JVS, X )







      INTEGER i, j
      REAL(kind=dp) JVS(1165), X(118), sum

      DO i=1,NVAR
         DO j = LU_CROW(i), LU_DIAG(i)-1 
             X(i) = X(i) - JVS(j)*X(LU_ICOL(j));
         END DO  
      END DO

      DO i=NVAR,1,-1
        sum = X(i);
        DO j = LU_DIAG(i)+1, LU_CROW(i+1)-1
          sum = sum - JVS(j)*X(LU_ICOL(j));
        END DO
        X(i) = sum/JVS(LU_DIAG(i));
      END DO
      
END SUBROUTINE saprc99_mosaic_4bin_iepox_vbs_aq_KppSolveIndirect


SUBROUTINE saprc99_mosaic_4bin_iepox_vbs_aq_KppSolveCmplx( JVS, X )







      INTEGER i, j
      DOUBLE COMPLEX JVS(1165), X(118), sum

      DO i=1,NVAR
         DO j = LU_CROW(i), LU_DIAG(i)-1 
             X(i) = X(i) - JVS(j)*X(LU_ICOL(j));
         END DO  
      END DO

      DO i=NVAR,1,-1
        sum = X(i);
        DO j = LU_DIAG(i)+1, LU_CROW(i+1)-1
          sum = sum - JVS(j)*X(LU_ICOL(j));
        END DO
        X(i) = sum/JVS(LU_DIAG(i));
      END DO
      
END SUBROUTINE saprc99_mosaic_4bin_iepox_vbs_aq_KppSolveCmplx













SUBROUTINE saprc99_mosaic_4bin_iepox_vbs_aq_KppSolve ( JVS, X )


  REAL(kind=dp) :: JVS(LU_NONZERO)

  REAL(kind=dp) :: X(NVAR)

  X(34) = X(34)-JVS(159)*X(33)
  X(37) = X(37)-JVS(169)*X(36)
  X(38) = X(38)-JVS(173)*X(36)
  X(39) = X(39)-JVS(177)*X(36)
  X(40) = X(40)-JVS(181)*X(37)-JVS(182)*X(38)-JVS(183)*X(39)
  X(46) = X(46)-JVS(203)*X(45)
  X(47) = X(47)-JVS(209)*X(46)
  X(51) = X(51)-JVS(229)*X(49)-JVS(230)*X(50)
  X(54) = X(54)-JVS(247)*X(52)-JVS(248)*X(53)
  X(64) = X(64)-JVS(281)*X(63)
  X(72) = X(72)-JVS(311)*X(68)-JVS(312)*X(71)
  X(73) = X(73)-JVS(315)*X(68)-JVS(316)*X(71)
  X(74) = X(74)-JVS(319)*X(68)-JVS(320)*X(71)
  X(77) = X(77)-JVS(329)*X(68)-JVS(330)*X(71)
  X(78) = X(78)-JVS(334)*X(67)
  X(80) = X(80)-JVS(345)*X(68)-JVS(346)*X(71)
  X(81) = X(81)-JVS(352)*X(71)
  X(82) = X(82)-JVS(358)*X(68)-JVS(359)*X(71)-JVS(360)*X(72)-JVS(361)*X(73)-JVS(362)*X(74)
  X(83) = X(83)-JVS(370)*X(72)-JVS(371)*X(73)-JVS(372)*X(75)-JVS(373)*X(76)-JVS(374)*X(77)-JVS(375)*X(82)
  X(84) = X(84)-JVS(395)*X(58)-JVS(396)*X(61)-JVS(397)*X(74)-JVS(398)*X(78)-JVS(399)*X(80)-JVS(400)*X(81)-JVS(401)*X(82)
  X(86) = X(86)-JVS(422)*X(56)-JVS(423)*X(62)-JVS(424)*X(63)-JVS(425)*X(64)-JVS(426)*X(75)
  X(88) = X(88)-JVS(442)*X(68)-JVS(443)*X(71)-JVS(444)*X(72)-JVS(445)*X(73)-JVS(446)*X(76)-JVS(447)*X(77)-JVS(448)*X(81)&
            &-JVS(449)*X(85)
  X(89) = X(89)-JVS(460)*X(74)-JVS(461)*X(81)
  X(92) = X(92)-JVS(482)*X(62)-JVS(483)*X(63)-JVS(484)*X(72)-JVS(485)*X(73)-JVS(486)*X(75)-JVS(487)*X(86)-JVS(488)*X(90)&
            &-JVS(489)*X(91)
  X(95) = X(95)-JVS(522)*X(94)
  X(96) = X(96)-JVS(529)*X(94)
  X(97) = X(97)-JVS(536)*X(35)-JVS(537)*X(55)-JVS(538)*X(62)-JVS(539)*X(63)-JVS(540)*X(75)-JVS(541)*X(85)-JVS(542)*X(87)&
            &-JVS(543)*X(90)-JVS(544)*X(96)
  X(99) = X(99)-JVS(565)*X(94)-JVS(566)*X(98)
  X(100) = X(100)-JVS(572)*X(55)-JVS(573)*X(62)-JVS(574)*X(63)-JVS(575)*X(65)-JVS(576)*X(66)-JVS(577)*X(69)-JVS(578)&
             &*X(75)-JVS(579)*X(76)-JVS(580)*X(85)-JVS(581)*X(86)-JVS(582)*X(87)-JVS(583)*X(88)-JVS(584)*X(89)-JVS(585)&
             &*X(90)-JVS(586)*X(91)-JVS(587)*X(93)-JVS(588)*X(94)-JVS(589)*X(95)-JVS(590)*X(96)-JVS(591)*X(98)-JVS(592)&
             &*X(99)
  X(101) = X(101)-JVS(611)*X(64)-JVS(612)*X(90)-JVS(613)*X(91)-JVS(614)*X(93)-JVS(615)*X(96)-JVS(616)*X(98)
  X(102) = X(102)-JVS(625)*X(31)-JVS(626)*X(85)-JVS(627)*X(87)-JVS(628)*X(90)-JVS(629)*X(91)-JVS(630)*X(93)-JVS(631)&
             &*X(94)-JVS(632)*X(95)-JVS(633)*X(98)-JVS(634)*X(99)
  X(103) = X(103)-JVS(641)*X(56)-JVS(642)*X(62)-JVS(643)*X(63)-JVS(644)*X(72)-JVS(645)*X(73)-JVS(646)*X(75)-JVS(647)&
             &*X(77)-JVS(648)*X(79)-JVS(649)*X(85)-JVS(650)*X(87)-JVS(651)*X(90)-JVS(652)*X(91)-JVS(653)*X(93)-JVS(654)&
             &*X(95)-JVS(655)*X(96)-JVS(656)*X(98)-JVS(657)*X(99)-JVS(658)*X(101)-JVS(659)*X(102)
  X(104) = X(104)-JVS(671)*X(62)-JVS(672)*X(63)-JVS(673)*X(75)-JVS(674)*X(87)-JVS(675)*X(90)-JVS(676)*X(95)-JVS(677)&
             &*X(96)-JVS(678)*X(98)-JVS(679)*X(99)-JVS(680)*X(101)-JVS(681)*X(102)
  X(105) = X(105)-JVS(692)*X(62)-JVS(693)*X(71)-JVS(694)*X(75)-JVS(695)*X(90)-JVS(696)*X(91)-JVS(697)*X(93)-JVS(698)&
             &*X(94)-JVS(699)*X(96)-JVS(700)*X(98)-JVS(701)*X(99)-JVS(702)*X(101)-JVS(703)*X(102)
  X(106) = X(106)-JVS(716)*X(76)-JVS(717)*X(77)-JVS(718)*X(85)-JVS(719)*X(87)-JVS(720)*X(90)-JVS(721)*X(91)-JVS(722)&
             &*X(93)-JVS(723)*X(94)-JVS(724)*X(95)-JVS(725)*X(96)-JVS(726)*X(98)-JVS(727)*X(99)-JVS(728)*X(102)
  X(107) = X(107)-JVS(739)*X(35)-JVS(740)*X(55)-JVS(741)*X(56)-JVS(742)*X(62)-JVS(743)*X(63)-JVS(744)*X(68)-JVS(745)&
             &*X(71)-JVS(746)*X(72)-JVS(747)*X(73)-JVS(748)*X(74)-JVS(749)*X(75)-JVS(750)*X(76)-JVS(751)*X(77)-JVS(752)&
             &*X(79)-JVS(753)*X(81)-JVS(754)*X(85)-JVS(755)*X(87)-JVS(756)*X(89)-JVS(757)*X(90)-JVS(758)*X(91)-JVS(759)&
             &*X(93)-JVS(760)*X(94)-JVS(761)*X(95)-JVS(762)*X(96)-JVS(763)*X(98)-JVS(764)*X(99)-JVS(765)*X(101)-JVS(766)&
             &*X(102)-JVS(767)*X(103)-JVS(768)*X(104)-JVS(769)*X(105)-JVS(770)*X(106)
  X(108) = X(108)-JVS(783)*X(42)-JVS(784)*X(88)-JVS(785)*X(89)-JVS(786)*X(91)-JVS(787)*X(93)-JVS(788)*X(95)-JVS(789)&
             &*X(96)-JVS(790)*X(98)-JVS(791)*X(99)-JVS(792)*X(102)-JVS(793)*X(103)-JVS(794)*X(104)-JVS(795)*X(105)-JVS(796)&
             &*X(106)-JVS(797)*X(107)
  X(109) = X(109)-JVS(809)*X(41)-JVS(810)*X(42)-JVS(811)*X(43)-JVS(812)*X(44)-JVS(813)*X(58)-JVS(814)*X(60)-JVS(815)&
             &*X(64)-JVS(816)*X(66)-JVS(817)*X(67)-JVS(818)*X(70)-JVS(819)*X(78)-JVS(820)*X(84)-JVS(821)*X(88)-JVS(822)&
             &*X(89)-JVS(823)*X(91)-JVS(824)*X(92)-JVS(825)*X(93)-JVS(826)*X(94)-JVS(827)*X(95)-JVS(828)*X(96)-JVS(829)&
             &*X(97)-JVS(830)*X(98)-JVS(831)*X(99)-JVS(832)*X(100)-JVS(833)*X(101)-JVS(834)*X(102)-JVS(835)*X(103)-JVS(836)&
             &*X(104)-JVS(837)*X(105)-JVS(838)*X(106)-JVS(839)*X(107)-JVS(840)*X(108)
  X(110) = X(110)-JVS(851)*X(32)-JVS(852)*X(64)-JVS(853)*X(65)-JVS(854)*X(75)-JVS(855)*X(85)-JVS(856)*X(86)-JVS(857)&
             &*X(87)-JVS(858)*X(90)-JVS(859)*X(91)-JVS(860)*X(93)-JVS(861)*X(94)-JVS(862)*X(97)-JVS(863)*X(98)-JVS(864)&
             &*X(99)-JVS(865)*X(101)-JVS(866)*X(102)-JVS(867)*X(103)-JVS(868)*X(104)-JVS(869)*X(105)-JVS(870)*X(106)&
             &-JVS(871)*X(107)-JVS(872)*X(108)-JVS(873)*X(109)
  X(111) = X(111)-JVS(883)*X(41)-JVS(884)*X(57)-JVS(885)*X(68)-JVS(886)*X(72)-JVS(887)*X(73)-JVS(888)*X(75)-JVS(889)&
             &*X(82)-JVS(890)*X(86)-JVS(891)*X(90)-JVS(892)*X(91)-JVS(893)*X(93)-JVS(894)*X(95)-JVS(895)*X(96)-JVS(896)&
             &*X(97)-JVS(897)*X(98)-JVS(898)*X(99)-JVS(899)*X(101)-JVS(900)*X(102)-JVS(901)*X(103)-JVS(902)*X(104)-JVS(903)&
             &*X(105)-JVS(904)*X(106)-JVS(905)*X(107)-JVS(906)*X(108)-JVS(907)*X(109)-JVS(908)*X(110)
  X(112) = X(112)-JVS(917)*X(48)-JVS(918)*X(55)-JVS(919)*X(59)-JVS(920)*X(60)-JVS(921)*X(61)-JVS(922)*X(65)-JVS(923)&
             &*X(66)-JVS(924)*X(67)-JVS(925)*X(68)-JVS(926)*X(69)-JVS(927)*X(70)-JVS(928)*X(71)-JVS(929)*X(72)-JVS(930)&
             &*X(73)-JVS(931)*X(76)-JVS(932)*X(77)-JVS(933)*X(78)-JVS(934)*X(79)-JVS(935)*X(82)-JVS(936)*X(83)-JVS(937)&
             &*X(85)-JVS(938)*X(87)-JVS(939)*X(88)-JVS(940)*X(89)-JVS(941)*X(90)-JVS(942)*X(91)-JVS(943)*X(93)-JVS(944)&
             &*X(94)-JVS(945)*X(95)-JVS(946)*X(96)-JVS(947)*X(97)-JVS(948)*X(98)-JVS(949)*X(99)-JVS(950)*X(100)-JVS(951)&
             &*X(101)-JVS(952)*X(102)-JVS(953)*X(103)-JVS(954)*X(104)-JVS(955)*X(105)-JVS(956)*X(106)-JVS(957)*X(107)&
             &-JVS(958)*X(108)-JVS(959)*X(109)-JVS(960)*X(110)-JVS(961)*X(111)
  X(113) = X(113)-JVS(969)*X(56)-JVS(970)*X(62)-JVS(971)*X(63)-JVS(972)*X(68)-JVS(973)*X(71)-JVS(974)*X(75)-JVS(975)&
             &*X(87)-JVS(976)*X(90)-JVS(977)*X(91)-JVS(978)*X(93)-JVS(979)*X(94)-JVS(980)*X(96)-JVS(981)*X(98)-JVS(982)&
             &*X(99)-JVS(983)*X(101)-JVS(984)*X(102)-JVS(985)*X(103)-JVS(986)*X(104)-JVS(987)*X(105)-JVS(988)*X(106)&
             &-JVS(989)*X(107)-JVS(990)*X(108)-JVS(991)*X(109)-JVS(992)*X(110)-JVS(993)*X(111)-JVS(994)*X(112)
  X(114) = X(114)-JVS(1001)*X(31)-JVS(1002)*X(32)-JVS(1003)*X(35)-JVS(1004)*X(48)-JVS(1005)*X(55)-JVS(1006)*X(56)&
             &-JVS(1007)*X(59)-JVS(1008)*X(60)-JVS(1009)*X(61)-JVS(1010)*X(62)-JVS(1011)*X(63)-JVS(1012)*X(65)-JVS(1013)&
             &*X(68)-JVS(1014)*X(69)-JVS(1015)*X(70)-JVS(1016)*X(71)-JVS(1017)*X(72)-JVS(1018)*X(73)-JVS(1019)*X(74)&
             &-JVS(1020)*X(75)-JVS(1021)*X(76)-JVS(1022)*X(77)-JVS(1023)*X(79)-JVS(1024)*X(80)-JVS(1025)*X(81)-JVS(1026)&
             &*X(82)-JVS(1027)*X(83)-JVS(1028)*X(84)-JVS(1029)*X(85)-JVS(1030)*X(86)-JVS(1031)*X(87)-JVS(1032)*X(88)&
             &-JVS(1033)*X(89)-JVS(1034)*X(90)-JVS(1035)*X(91)-JVS(1036)*X(93)-JVS(1037)*X(94)-JVS(1038)*X(95)-JVS(1039)&
             &*X(96)-JVS(1040)*X(97)-JVS(1041)*X(98)-JVS(1042)*X(99)-JVS(1043)*X(100)-JVS(1044)*X(101)-JVS(1045)*X(102)&
             &-JVS(1046)*X(103)-JVS(1047)*X(104)-JVS(1048)*X(105)-JVS(1049)*X(106)-JVS(1050)*X(107)-JVS(1051)*X(108)&
             &-JVS(1052)*X(109)-JVS(1053)*X(110)-JVS(1054)*X(111)-JVS(1055)*X(112)-JVS(1056)*X(113)
  X(115) = X(115)-JVS(1062)*X(60)-JVS(1063)*X(66)-JVS(1064)*X(92)-JVS(1065)*X(93)-JVS(1066)*X(94)-JVS(1067)*X(98)&
             &-JVS(1068)*X(99)-JVS(1069)*X(100)-JVS(1070)*X(101)-JVS(1071)*X(102)-JVS(1072)*X(104)-JVS(1073)*X(105)&
             &-JVS(1074)*X(106)-JVS(1075)*X(107)-JVS(1076)*X(108)-JVS(1077)*X(109)-JVS(1078)*X(110)-JVS(1079)*X(111)&
             &-JVS(1080)*X(112)-JVS(1081)*X(113)-JVS(1082)*X(114)
  X(116) = X(116)-JVS(1087)*X(44)-JVS(1088)*X(94)-JVS(1089)*X(95)-JVS(1090)*X(96)-JVS(1091)*X(98)-JVS(1092)*X(99)&
             &-JVS(1093)*X(102)-JVS(1094)*X(106)-JVS(1095)*X(107)-JVS(1096)*X(108)-JVS(1097)*X(109)-JVS(1098)*X(110)&
             &-JVS(1099)*X(111)-JVS(1100)*X(112)-JVS(1101)*X(113)-JVS(1102)*X(114)-JVS(1103)*X(115)
  X(117) = X(117)-JVS(1107)*X(43)-JVS(1108)*X(80)-JVS(1109)*X(98)-JVS(1110)*X(102)-JVS(1111)*X(106)-JVS(1112)*X(107)&
             &-JVS(1113)*X(108)-JVS(1114)*X(109)-JVS(1115)*X(110)-JVS(1116)*X(111)-JVS(1117)*X(112)-JVS(1118)*X(113)&
             &-JVS(1119)*X(114)-JVS(1120)*X(115)-JVS(1121)*X(116)
  X(118) = X(118)-JVS(1124)*X(58)-JVS(1125)*X(61)-JVS(1126)*X(70)-JVS(1127)*X(74)-JVS(1128)*X(78)-JVS(1129)*X(80)&
             &-JVS(1130)*X(81)-JVS(1131)*X(82)-JVS(1132)*X(84)-JVS(1133)*X(85)-JVS(1134)*X(87)-JVS(1135)*X(88)-JVS(1136)&
             &*X(89)-JVS(1137)*X(90)-JVS(1138)*X(91)-JVS(1139)*X(92)-JVS(1140)*X(93)-JVS(1141)*X(94)-JVS(1142)*X(95)&
             &-JVS(1143)*X(96)-JVS(1144)*X(97)-JVS(1145)*X(98)-JVS(1146)*X(99)-JVS(1147)*X(100)-JVS(1148)*X(101)-JVS(1149)&
             &*X(102)-JVS(1150)*X(103)-JVS(1151)*X(104)-JVS(1152)*X(105)-JVS(1153)*X(106)-JVS(1154)*X(107)-JVS(1155)*X(108)&
             &-JVS(1156)*X(109)-JVS(1157)*X(110)-JVS(1158)*X(111)-JVS(1159)*X(112)-JVS(1160)*X(113)-JVS(1161)*X(114)&
             &-JVS(1162)*X(115)-JVS(1163)*X(116)-JVS(1164)*X(117)
  X(118) = X(118)/JVS(1165)
  X(117) = (X(117)-JVS(1123)*X(118))/(JVS(1122))
  X(116) = (X(116)-JVS(1105)*X(117)-JVS(1106)*X(118))/(JVS(1104))
  X(115) = (X(115)-JVS(1084)*X(116)-JVS(1085)*X(117)-JVS(1086)*X(118))/(JVS(1083))
  X(114) = (X(114)-JVS(1058)*X(115)-JVS(1059)*X(116)-JVS(1060)*X(117)-JVS(1061)*X(118))/(JVS(1057))
  X(113) = (X(113)-JVS(996)*X(114)-JVS(997)*X(115)-JVS(998)*X(116)-JVS(999)*X(117)-JVS(1000)*X(118))/(JVS(995))
  X(112) = (X(112)-JVS(963)*X(113)-JVS(964)*X(114)-JVS(965)*X(115)-JVS(966)*X(116)-JVS(967)*X(117)-JVS(968)*X(118))&
             &/(JVS(962))
  X(111) = (X(111)-JVS(910)*X(112)-JVS(911)*X(113)-JVS(912)*X(114)-JVS(913)*X(115)-JVS(914)*X(116)-JVS(915)*X(117)&
             &-JVS(916)*X(118))/(JVS(909))
  X(110) = (X(110)-JVS(875)*X(111)-JVS(876)*X(112)-JVS(877)*X(113)-JVS(878)*X(114)-JVS(879)*X(115)-JVS(880)*X(116)&
             &-JVS(881)*X(117)-JVS(882)*X(118))/(JVS(874))
  X(109) = (X(109)-JVS(842)*X(110)-JVS(843)*X(111)-JVS(844)*X(112)-JVS(845)*X(113)-JVS(846)*X(114)-JVS(847)*X(115)&
             &-JVS(848)*X(116)-JVS(849)*X(117)-JVS(850)*X(118))/(JVS(841))
  X(108) = (X(108)-JVS(799)*X(109)-JVS(800)*X(110)-JVS(801)*X(111)-JVS(802)*X(112)-JVS(803)*X(113)-JVS(804)*X(114)&
             &-JVS(805)*X(115)-JVS(806)*X(116)-JVS(807)*X(117)-JVS(808)*X(118))/(JVS(798))
  X(107) = (X(107)-JVS(772)*X(108)-JVS(773)*X(109)-JVS(774)*X(110)-JVS(775)*X(111)-JVS(776)*X(112)-JVS(777)*X(113)&
             &-JVS(778)*X(114)-JVS(779)*X(115)-JVS(780)*X(116)-JVS(781)*X(117)-JVS(782)*X(118))/(JVS(771))
  X(106) = (X(106)-JVS(730)*X(108)-JVS(731)*X(109)-JVS(732)*X(111)-JVS(733)*X(112)-JVS(734)*X(114)-JVS(735)*X(115)&
             &-JVS(736)*X(116)-JVS(737)*X(117)-JVS(738)*X(118))/(JVS(729))
  X(105) = (X(105)-JVS(705)*X(106)-JVS(706)*X(107)-JVS(707)*X(108)-JVS(708)*X(109)-JVS(709)*X(110)-JVS(710)*X(111)&
             &-JVS(711)*X(113)-JVS(712)*X(114)-JVS(713)*X(115)-JVS(714)*X(117)-JVS(715)*X(118))/(JVS(704))
  X(104) = (X(104)-JVS(683)*X(105)-JVS(684)*X(106)-JVS(685)*X(107)-JVS(686)*X(109)-JVS(687)*X(110)-JVS(688)*X(113)&
             &-JVS(689)*X(114)-JVS(690)*X(115)-JVS(691)*X(118))/(JVS(682))
  X(103) = (X(103)-JVS(661)*X(104)-JVS(662)*X(105)-JVS(663)*X(106)-JVS(664)*X(107)-JVS(665)*X(109)-JVS(666)*X(112)&
             &-JVS(667)*X(113)-JVS(668)*X(114)-JVS(669)*X(115)-JVS(670)*X(118))/(JVS(660))
  X(102) = (X(102)-JVS(636)*X(106)-JVS(637)*X(109)-JVS(638)*X(114)-JVS(639)*X(115)-JVS(640)*X(118))/(JVS(635))
  X(101) = (X(101)-JVS(618)*X(102)-JVS(619)*X(106)-JVS(620)*X(109)-JVS(621)*X(113)-JVS(622)*X(114)-JVS(623)*X(115)&
             &-JVS(624)*X(118))/(JVS(617))
  X(100) = (X(100)-JVS(594)*X(101)-JVS(595)*X(102)-JVS(596)*X(104)-JVS(597)*X(105)-JVS(598)*X(106)-JVS(599)*X(107)&
             &-JVS(600)*X(108)-JVS(601)*X(109)-JVS(602)*X(110)-JVS(603)*X(111)-JVS(604)*X(112)-JVS(605)*X(113)-JVS(606)&
             &*X(114)-JVS(607)*X(115)-JVS(608)*X(116)-JVS(609)*X(117)-JVS(610)*X(118))/(JVS(593))
  X(99) = (X(99)-JVS(568)*X(102)-JVS(569)*X(106)-JVS(570)*X(114)-JVS(571)*X(118))/(JVS(567))
  X(98) = (X(98)-JVS(561)*X(102)-JVS(562)*X(106)-JVS(563)*X(114)-JVS(564)*X(118))/(JVS(560))
  X(97) = (X(97)-JVS(546)*X(98)-JVS(547)*X(101)-JVS(548)*X(102)-JVS(549)*X(103)-JVS(550)*X(104)-JVS(551)*X(105)-JVS(552)&
            &*X(106)-JVS(553)*X(108)-JVS(554)*X(111)-JVS(555)*X(114)-JVS(556)*X(115)-JVS(557)*X(116)-JVS(558)*X(117)&
            &-JVS(559)*X(118))/(JVS(545))
  X(96) = (X(96)-JVS(531)*X(98)-JVS(532)*X(102)-JVS(533)*X(106)-JVS(534)*X(114)-JVS(535)*X(118))/(JVS(530))
  X(95) = (X(95)-JVS(524)*X(98)-JVS(525)*X(102)-JVS(526)*X(106)-JVS(527)*X(114)-JVS(528)*X(118))/(JVS(523))
  X(94) = (X(94)-JVS(518)*X(102)-JVS(519)*X(106)-JVS(520)*X(114)-JVS(521)*X(118))/(JVS(517))
  X(93) = (X(93)-JVS(513)*X(102)-JVS(514)*X(106)-JVS(515)*X(114)-JVS(516)*X(118))/(JVS(512))
  X(92) = (X(92)-JVS(491)*X(93)-JVS(492)*X(94)-JVS(493)*X(98)-JVS(494)*X(99)-JVS(495)*X(101)-JVS(496)*X(102)-JVS(497)&
            &*X(104)-JVS(498)*X(105)-JVS(499)*X(106)-JVS(500)*X(107)-JVS(501)*X(108)-JVS(502)*X(109)-JVS(503)*X(110)&
            &-JVS(504)*X(111)-JVS(505)*X(112)-JVS(506)*X(113)-JVS(507)*X(114)-JVS(508)*X(115)-JVS(509)*X(116)-JVS(510)&
            &*X(117)-JVS(511)*X(118))/(JVS(490))
  X(91) = (X(91)-JVS(478)*X(102)-JVS(479)*X(106)-JVS(480)*X(114)-JVS(481)*X(118))/(JVS(477))
  X(90) = (X(90)-JVS(473)*X(102)-JVS(474)*X(106)-JVS(475)*X(114)-JVS(476)*X(118))/(JVS(472))
  X(89) = (X(89)-JVS(463)*X(108)-JVS(464)*X(109)-JVS(465)*X(111)-JVS(466)*X(112)-JVS(467)*X(114)-JVS(468)*X(115)&
            &-JVS(469)*X(116)-JVS(470)*X(117)-JVS(471)*X(118))/(JVS(462))
  X(88) = (X(88)-JVS(451)*X(89)-JVS(452)*X(91)-JVS(453)*X(93)-JVS(454)*X(96)-JVS(455)*X(102)-JVS(456)*X(106)-JVS(457)&
            &*X(112)-JVS(458)*X(114)-JVS(459)*X(118))/(JVS(450))
  X(87) = (X(87)-JVS(438)*X(102)-JVS(439)*X(106)-JVS(440)*X(114)-JVS(441)*X(118))/(JVS(437))
  X(86) = (X(86)-JVS(428)*X(90)-JVS(429)*X(91)-JVS(430)*X(93)-JVS(431)*X(98)-JVS(432)*X(101)-JVS(433)*X(106)-JVS(434)&
            &*X(109)-JVS(435)*X(114)-JVS(436)*X(118))/(JVS(427))
  X(85) = (X(85)-JVS(418)*X(102)-JVS(419)*X(106)-JVS(420)*X(114)-JVS(421)*X(118))/(JVS(417))
  X(84) = (X(84)-JVS(403)*X(88)-JVS(404)*X(89)-JVS(405)*X(95)-JVS(406)*X(96)-JVS(407)*X(97)-JVS(408)*X(98)-JVS(409)&
            &*X(99)-JVS(410)*X(100)-JVS(411)*X(103)-JVS(412)*X(106)-JVS(413)*X(109)-JVS(414)*X(112)-JVS(415)*X(114)-JVS(416)&
            &*X(118))/(JVS(402))
  X(83) = (X(83)-JVS(377)*X(85)-JVS(378)*X(87)-JVS(379)*X(88)-JVS(380)*X(90)-JVS(381)*X(91)-JVS(382)*X(93)-JVS(383)&
            &*X(94)-JVS(384)*X(95)-JVS(385)*X(96)-JVS(386)*X(97)-JVS(387)*X(98)-JVS(388)*X(99)-JVS(389)*X(100)-JVS(390)&
            &*X(102)-JVS(391)*X(103)-JVS(392)*X(106)-JVS(393)*X(114)-JVS(394)*X(118))/(JVS(376))
  X(82) = (X(82)-JVS(364)*X(95)-JVS(365)*X(96)-JVS(366)*X(99)-JVS(367)*X(106)-JVS(368)*X(114)-JVS(369)*X(118))&
            &/(JVS(363))
  X(81) = (X(81)-JVS(354)*X(89)-JVS(355)*X(112)-JVS(356)*X(114)-JVS(357)*X(118))/(JVS(353))
  X(80) = (X(80)-JVS(348)*X(98)-JVS(349)*X(106)-JVS(350)*X(114)-JVS(351)*X(118))/(JVS(347))
  X(79) = (X(79)-JVS(341)*X(107)-JVS(342)*X(112)-JVS(343)*X(113)-JVS(344)*X(114))/(JVS(340))
  X(78) = (X(78)-JVS(336)*X(89)-JVS(337)*X(109)-JVS(338)*X(112)-JVS(339)*X(118))/(JVS(335))
  X(77) = (X(77)-JVS(332)*X(106)-JVS(333)*X(114))/(JVS(331))
  X(76) = (X(76)-JVS(327)*X(106)-JVS(328)*X(114))/(JVS(326))
  X(75) = (X(75)-JVS(325)*X(114))/(JVS(324))
  X(74) = (X(74)-JVS(322)*X(114)-JVS(323)*X(118))/(JVS(321))
  X(73) = (X(73)-JVS(318)*X(114))/(JVS(317))
  X(72) = (X(72)-JVS(314)*X(114))/(JVS(313))
  X(71) = (X(71)-JVS(310)*X(114))/(JVS(309))
  X(70) = (X(70)-JVS(306)*X(109)-JVS(307)*X(112)-JVS(308)*X(114))/(JVS(305))
  X(69) = (X(69)-JVS(301)*X(107)-JVS(302)*X(110)-JVS(303)*X(113)-JVS(304)*X(114))/(JVS(300))
  X(68) = (X(68)-JVS(299)*X(114))/(JVS(298))
  X(67) = (X(67)-JVS(294)*X(78)-JVS(295)*X(109)-JVS(296)*X(112)-JVS(297)*X(118))/(JVS(293))
  X(66) = (X(66)-JVS(290)*X(100)-JVS(291)*X(112)-JVS(292)*X(115))/(JVS(289))
  X(65) = (X(65)-JVS(286)*X(110)-JVS(287)*X(112)-JVS(288)*X(114))/(JVS(285))
  X(64) = (X(64)-JVS(283)*X(109)-JVS(284)*X(114))/(JVS(282))
  X(63) = (X(63)-JVS(280)*X(114))/(JVS(279))
  X(62) = (X(62)-JVS(278)*X(114))/(JVS(277))
  X(61) = (X(61)-JVS(275)*X(114)-JVS(276)*X(118))/(JVS(274))
  X(60) = (X(60)-JVS(272)*X(114)-JVS(273)*X(115))/(JVS(271))
  X(59) = (X(59)-JVS(268)*X(61)-JVS(269)*X(114)-JVS(270)*X(118))/(JVS(267))
  X(58) = (X(58)-JVS(265)*X(109)-JVS(266)*X(118))/(JVS(264))
  X(57) = (X(57)-JVS(259)*X(68)-JVS(260)*X(91)-JVS(261)*X(93)-JVS(262)*X(106)-JVS(263)*X(114))/(JVS(258))
  X(56) = (X(56)-JVS(257)*X(114))/(JVS(256))
  X(55) = (X(55)-JVS(255)*X(114))/(JVS(254))
  X(54) = (X(54)-JVS(250)*X(91)-JVS(251)*X(106)-JVS(252)*X(114)-JVS(253)*X(118))/(JVS(249))
  X(53) = (X(53)-JVS(243)*X(54)-JVS(244)*X(91)-JVS(245)*X(106)-JVS(246)*X(114))/(JVS(242))
  X(52) = (X(52)-JVS(237)*X(53)-JVS(238)*X(91)-JVS(239)*X(106)-JVS(240)*X(114)-JVS(241)*X(118))/(JVS(236))
  X(51) = (X(51)-JVS(232)*X(93)-JVS(233)*X(106)-JVS(234)*X(114)-JVS(235)*X(118))/(JVS(231))
  X(50) = (X(50)-JVS(225)*X(51)-JVS(226)*X(93)-JVS(227)*X(106)-JVS(228)*X(114))/(JVS(224))
  X(49) = (X(49)-JVS(219)*X(50)-JVS(220)*X(93)-JVS(221)*X(106)-JVS(222)*X(114)-JVS(223)*X(118))/(JVS(218))
  X(48) = (X(48)-JVS(216)*X(112)-JVS(217)*X(114))/(JVS(215))
  X(47) = (X(47)-JVS(211)*X(94)-JVS(212)*X(106)-JVS(213)*X(114)-JVS(214)*X(118))/(JVS(210))
  X(46) = (X(46)-JVS(205)*X(47)-JVS(206)*X(94)-JVS(207)*X(114)-JVS(208)*X(118))/(JVS(204))
  X(45) = (X(45)-JVS(199)*X(47)-JVS(200)*X(94)-JVS(201)*X(114)-JVS(202)*X(118))/(JVS(198))
  X(44) = (X(44)-JVS(196)*X(109)-JVS(197)*X(116))/(JVS(195))
  X(43) = (X(43)-JVS(193)*X(109)-JVS(194)*X(117))/(JVS(192))
  X(42) = (X(42)-JVS(190)*X(108)-JVS(191)*X(109))/(JVS(189))
  X(41) = (X(41)-JVS(187)*X(109)-JVS(188)*X(111))/(JVS(186))
  X(40) = (X(40)-JVS(185)*X(114))/(JVS(184))
  X(39) = (X(39)-JVS(179)*X(40)-JVS(180)*X(114))/(JVS(178))
  X(38) = (X(38)-JVS(175)*X(39)-JVS(176)*X(114))/(JVS(174))
  X(37) = (X(37)-JVS(171)*X(38)-JVS(172)*X(114))/(JVS(170))
  X(36) = (X(36)-JVS(168)*X(114))/(JVS(167))
  X(35) = (X(35)-JVS(166)*X(114))/(JVS(165))
  X(34) = (X(34)-JVS(161)*X(94)-JVS(162)*X(112)-JVS(163)*X(114)-JVS(164)*X(115))/(JVS(160))
  X(33) = (X(33)-JVS(156)*X(34)-JVS(157)*X(112)-JVS(158)*X(114))/(JVS(155))
  X(32) = (X(32)-JVS(154)*X(114))/(JVS(153))
  X(31) = (X(31)-JVS(152)*X(106))/(JVS(151))
  X(30) = (X(30)-JVS(150)*X(114))/(JVS(149))
  X(29) = (X(29)-JVS(147)*X(30)-JVS(148)*X(114))/(JVS(146))
  X(28) = (X(28)-JVS(144)*X(30)-JVS(145)*X(114))/(JVS(143))
  X(27) = (X(27)-JVS(141)*X(30)-JVS(142)*X(114))/(JVS(140))
  X(26) = (X(26)-JVS(139)*X(114))/(JVS(138))
  X(25) = (X(25)-JVS(136)*X(26)-JVS(137)*X(114))/(JVS(135))
  X(24) = (X(24)-JVS(133)*X(26)-JVS(134)*X(114))/(JVS(132))
  X(23) = (X(23)-JVS(130)*X(26)-JVS(131)*X(114))/(JVS(129))
  X(22) = (X(22)-JVS(127)*X(26)-JVS(128)*X(114))/(JVS(126))
  X(21) = (X(21)-JVS(123)*X(36)-JVS(124)*X(37)-JVS(125)*X(114))/(JVS(122))
  X(20) = (X(20)-JVS(121)*X(114))/(JVS(120))
  X(19) = (X(19)-JVS(118)*X(33)-JVS(119)*X(114))/(JVS(117))
  X(18) = (X(18)-JVS(113)*X(33)-JVS(114)*X(112)-JVS(115)*X(114)-JVS(116)*X(115))/(JVS(112))
  X(17) = (X(17)-JVS(106)*X(18)-JVS(107)*X(19)-JVS(108)*X(20)-JVS(109)*X(34)-JVS(110)*X(114)-JVS(111)*X(115))/(JVS(105))
  X(16) = (X(16)-JVS(100)*X(52)-JVS(101)*X(91)-JVS(102)*X(106)-JVS(103)*X(114)-JVS(104)*X(118))/(JVS(99))
  X(15) = (X(15)-JVS(96)*X(45)-JVS(97)*X(94)-JVS(98)*X(114))/(JVS(95))
  X(14) = (X(14)-JVS(94)*X(18))/(JVS(93))
  X(13) = (X(13)-JVS(91)*X(18)-JVS(92)*X(112))/(JVS(90))
  X(12) = (X(12)-JVS(85)*X(49)-JVS(86)*X(93)-JVS(87)*X(106)-JVS(88)*X(114)-JVS(89)*X(118))/(JVS(84))
  X(11) = (X(11)-JVS(79)*X(92)-JVS(80)*X(107)-JVS(81)*X(110)-JVS(82)*X(112)-JVS(83)*X(113))/(JVS(78))
  X(10) = (X(10)-JVS(72)*X(92)-JVS(73)*X(107)-JVS(74)*X(110)-JVS(75)*X(113)-JVS(76)*X(115)-JVS(77)*X(118))/(JVS(71))
  X(9) = (X(9)-JVS(63)*X(67)-JVS(64)*X(80)-JVS(65)*X(87)-JVS(66)*X(102)-JVS(67)*X(106)-JVS(68)*X(109)-JVS(69)*X(114)&
           &-JVS(70)*X(118))/(JVS(62))
  X(8) = (X(8)-JVS(58)*X(67)-JVS(59)*X(87)-JVS(60)*X(109)-JVS(61)*X(118))/(JVS(57))
  X(7) = (X(7)-JVS(53)*X(108)-JVS(54)*X(112)-JVS(55)*X(116)-JVS(56)*X(117))/(JVS(52))
  X(6) = (X(6)-JVS(50)*X(111)-JVS(51)*X(112))/(JVS(49))
  X(5) = (X(5)-JVS(46)*X(76)-JVS(47)*X(87)-JVS(48)*X(106))/(JVS(45))
  X(4) = (X(4)-JVS(31)*X(90)-JVS(32)*X(91)-JVS(33)*X(93)-JVS(34)*X(94)-JVS(35)*X(96)-JVS(36)*X(98)-JVS(37)*X(106)&
           &-JVS(38)*X(107)-JVS(39)*X(108)-JVS(40)*X(110)-JVS(41)*X(112)-JVS(42)*X(113)-JVS(43)*X(116)-JVS(44)*X(117))&
           &/(JVS(30))
  X(3) = (X(3)-JVS(21)*X(87)-JVS(22)*X(90)-JVS(23)*X(98)-JVS(24)*X(106)-JVS(25)*X(107)-JVS(26)*X(110)-JVS(27)*X(111)&
           &-JVS(28)*X(112)-JVS(29)*X(113))/(JVS(20))
  X(2) = (X(2)-JVS(5)*X(66)-JVS(6)*X(76)-JVS(7)*X(85)-JVS(8)*X(87)-JVS(9)*X(90)-JVS(10)*X(91)-JVS(11)*X(93)-JVS(12)&
           &*X(94)-JVS(13)*X(95)-JVS(14)*X(96)-JVS(15)*X(98)-JVS(16)*X(99)-JVS(17)*X(106)-JVS(18)*X(114)-JVS(19)*X(115))&
           &/(JVS(4))
  X(1) = (X(1)-JVS(2)*X(59)-JVS(3)*X(114))/(JVS(1))
      
END SUBROUTINE saprc99_mosaic_4bin_iepox_vbs_aq_KppSolve
























      SUBROUTINE saprc99_mosaic_4bin_iepox_vbs_aq_WCOPY(N,X,incX,Y,incY)








      
      INTEGER i,incX,incY,M,MP1,N
      REAL(kind=dp) X(N),Y(N)

      IF (N.LE.0) RETURN

      M = MOD(N,8)
      IF( M .NE. 0 ) THEN
        DO i = 1,M
          Y(i) = X(i)
        END DO
        IF( N .LT. 8 ) RETURN
      END IF    
      MP1 = M+1
      DO i = MP1,N,8
        Y(i) = X(i)
        Y(i + 1) = X(i + 1)
        Y(i + 2) = X(i + 2)
        Y(i + 3) = X(i + 3)
        Y(i + 4) = X(i + 4)
        Y(i + 5) = X(i + 5)
        Y(i + 6) = X(i + 6)
        Y(i + 7) = X(i + 7)
      END DO

      END SUBROUTINE saprc99_mosaic_4bin_iepox_vbs_aq_WCOPY



      SUBROUTINE saprc99_mosaic_4bin_iepox_vbs_aq_WAXPY(N,Alpha,X,incX,Y,incY)









      INTEGER i,incX,incY,M,MP1,N
      REAL(kind=dp) X(N),Y(N),Alpha
      REAL(kind=dp) ZERO
      PARAMETER( ZERO = 0.0_dp )

      IF (Alpha .EQ. ZERO) RETURN
      IF (N .LE. 0) RETURN

      M = MOD(N,4)
      IF( M .NE. 0 ) THEN
        DO i = 1,M
          Y(i) = Y(i) + Alpha*X(i)
        END DO
        IF( N .LT. 4 ) RETURN
      END IF
      MP1 = M + 1
      DO i = MP1,N,4
        Y(i) = Y(i) + Alpha*X(i)
        Y(i + 1) = Y(i + 1) + Alpha*X(i + 1)
        Y(i + 2) = Y(i + 2) + Alpha*X(i + 2)
        Y(i + 3) = Y(i + 3) + Alpha*X(i + 3)
      END DO
      
      END SUBROUTINE saprc99_mosaic_4bin_iepox_vbs_aq_WAXPY




      SUBROUTINE saprc99_mosaic_4bin_iepox_vbs_aq_WSCAL(N,Alpha,X,incX)









      INTEGER i,incX,M,MP1,N
      REAL(kind=dp) X(N),Alpha
      REAL(kind=dp) ZERO, ONE
      PARAMETER( ZERO = 0.0_dp ) 
      PARAMETER( ONE  = 1.0_dp )

      IF (Alpha .EQ. ONE) RETURN
      IF (N .LE. 0) RETURN

      M = MOD(N,5)
      IF( M .NE. 0 ) THEN
        IF (Alpha .EQ. (-ONE)) THEN
          DO i = 1,M
            X(i) = -X(i)
          END DO
        ELSEIF (Alpha .EQ. ZERO) THEN
          DO i = 1,M
            X(i) = ZERO
          END DO
        ELSE
          DO i = 1,M
            X(i) = Alpha*X(i)
          END DO
        END IF
        IF( N .LT. 5 ) RETURN
      END IF
      MP1 = M + 1
      IF (Alpha .EQ. (-ONE)) THEN
        DO i = MP1,N,5
          X(i)     = -X(i)
          X(i + 1) = -X(i + 1)
          X(i + 2) = -X(i + 2)
          X(i + 3) = -X(i + 3)
          X(i + 4) = -X(i + 4)
        END DO
      ELSEIF (Alpha .EQ. ZERO) THEN
        DO i = MP1,N,5
          X(i)     = ZERO
          X(i + 1) = ZERO
          X(i + 2) = ZERO
          X(i + 3) = ZERO
          X(i + 4) = ZERO
        END DO
      ELSE
        DO i = MP1,N,5
          X(i)     = Alpha*X(i)
          X(i + 1) = Alpha*X(i + 1)
          X(i + 2) = Alpha*X(i + 2)
          X(i + 3) = Alpha*X(i + 3)
          X(i + 4) = Alpha*X(i + 4)
        END DO
      END IF

      END SUBROUTINE saprc99_mosaic_4bin_iepox_vbs_aq_WSCAL


      REAL(kind=dp) FUNCTION saprc99_mosaic_4bin_iepox_vbs_aq_WLAMCH( C )








      CHARACTER C
      INTEGER   i
      REAL(kind=dp)  ONE, HALF, Eps, Sum
      PARAMETER (ONE  = 1.0_dp)
      PARAMETER (HALF = 0.5_dp)
      LOGICAL   First
      SAVE     First, Eps
      DATA     First /.TRUE./
      
      IF (First) THEN
        First = .FALSE.
        Eps = HALF**(16)
        DO i = 17, 80
          Eps = Eps*HALF
          CALL saprc99_mosaic_4bin_iepox_vbs_aq_WLAMCH_ADD(ONE,Eps,Sum)
          IF (Sum.LE.ONE) GOTO 10
        END DO
        PRINT*,'ERROR IN WLAMCH. EPS < ',Eps
        RETURN
10      Eps = Eps*2
        i = i-1      
      END IF

      saprc99_mosaic_4bin_iepox_vbs_aq_WLAMCH = Eps

      END FUNCTION saprc99_mosaic_4bin_iepox_vbs_aq_WLAMCH
     
      SUBROUTINE saprc99_mosaic_4bin_iepox_vbs_aq_WLAMCH_ADD( A, B, Sum )

      
      REAL(kind=dp) A, B, Sum
      Sum = A + B

      END SUBROUTINE saprc99_mosaic_4bin_iepox_vbs_aq_WLAMCH_ADD




      SUBROUTINE saprc99_mosaic_4bin_iepox_vbs_aq_SET2ZERO(N,Y)




      
      INTEGER ::  i,M,MP1,N
      REAL(kind=dp) ::  Y(N)
      REAL(kind=dp), PARAMETER :: ZERO = 0.0d0

      IF (N.LE.0) RETURN

      M = MOD(N,8)
      IF( M .NE. 0 ) THEN
        DO i = 1,M
          Y(i) = ZERO
        END DO
        IF( N .LT. 8 ) RETURN
      END IF    
      MP1 = M+1
      DO i = MP1,N,8
        Y(i)     = ZERO
        Y(i + 1) = ZERO
        Y(i + 2) = ZERO
        Y(i + 3) = ZERO
        Y(i + 4) = ZERO
        Y(i + 5) = ZERO
        Y(i + 6) = ZERO
        Y(i + 7) = ZERO
      END DO

      END SUBROUTINE saprc99_mosaic_4bin_iepox_vbs_aq_SET2ZERO



      REAL(kind=dp) FUNCTION saprc99_mosaic_4bin_iepox_vbs_aq_WDOT (N, DX, incX, DY, incY) 









      IMPLICIT NONE
      INTEGER :: N, incX, incY
      REAL(kind=dp) :: DX(N), DY(N) 

      INTEGER :: i, IX, IY, M, MP1, NS
                                 
      saprc99_mosaic_4bin_iepox_vbs_aq_WDOT = 0.0D0 
      IF (N .LE. 0) RETURN 
      IF (incX .EQ. incY) IF (incX-1) 5,20,60 



    5 IX = 1 
      IY = 1 
      IF (incX .LT. 0) IX = (-N+1)*incX + 1 
      IF (incY .LT. 0) IY = (-N+1)*incY + 1 
      DO i = 1,N 
        saprc99_mosaic_4bin_iepox_vbs_aq_WDOT = saprc99_mosaic_4bin_iepox_vbs_aq_WDOT + DX(IX)*DY(IY) 
        IX = IX + incX 
        IY = IY + incY 
      END DO 
      RETURN 





   20 M = MOD(N,5) 
      IF (M .EQ. 0) GO TO 40 
      DO i = 1,M 
         saprc99_mosaic_4bin_iepox_vbs_aq_WDOT = saprc99_mosaic_4bin_iepox_vbs_aq_WDOT + DX(i)*DY(i) 
      END DO 
      IF (N .LT. 5) RETURN 
   40 MP1 = M + 1 
      DO i = MP1,N,5 
          saprc99_mosaic_4bin_iepox_vbs_aq_WDOT = saprc99_mosaic_4bin_iepox_vbs_aq_WDOT + DX(i)*DY(i) + DX(i+1)*DY(i+1) +&
                   DX(i+2)*DY(i+2) +  &
                   DX(i+3)*DY(i+3) + DX(i+4)*DY(i+4)                   
      END DO 
      RETURN 



   60 NS = N*incX 
      DO i = 1,NS,incX 
        saprc99_mosaic_4bin_iepox_vbs_aq_WDOT = saprc99_mosaic_4bin_iepox_vbs_aq_WDOT + DX(i)*DY(i) 
      END DO 

      END FUNCTION saprc99_mosaic_4bin_iepox_vbs_aq_WDOT                                          




END MODULE saprc99_mosaic_4bin_iepox_vbs_aq_Integrator
