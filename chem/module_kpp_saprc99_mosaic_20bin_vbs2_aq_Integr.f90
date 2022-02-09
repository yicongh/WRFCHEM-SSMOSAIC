
































MODULE saprc99_mosaic_20bin_vbs2_aq_Integrator

 USE saprc99_mosaic_20bin_vbs2_aq_Parameters
 USE saprc99_mosaic_20bin_vbs2_aq_Precision
 USE saprc99_mosaic_20bin_vbs2_aq_JacobianSP

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

SUBROUTINE  saprc99_mosaic_20bin_vbs2_aq_INTEGRATE( TIN, TOUT, &
  FIX, VAR,  RCONST, ATOL, RTOL, IRR_WRK,  &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U  )

   USE saprc99_mosaic_20bin_vbs2_aq_Parameters

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

   CALL saprc99_mosaic_20bin_vbs2_aq_Rosenbrock(VAR, FIX, RCONST, TIN,TOUT,   &
         ATOL,RTOL,               &
         RCNTRL,ICNTRL,RSTATUS,ISTATUS,IRR_WRK,IERR)

   STEPMIN = RCNTRL(ihexit)
   
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U))    IERR_U       = IERR

END SUBROUTINE  saprc99_mosaic_20bin_vbs2_aq_INTEGRATE


SUBROUTINE  saprc99_mosaic_20bin_vbs2_aq_Rosenbrock(Y, FIX, RCONST, Tstart,Tend, &
           AbsTol,RelTol,            &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IRR_WRK,IERR)







































































































  USE saprc99_mosaic_20bin_vbs2_aq_Parameters

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
      CALL saprc99_mosaic_20bin_vbs2_aq_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
      RETURN
   END IF


   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 100000
   ELSEIF (ICNTRL(4) > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL saprc99_mosaic_20bin_vbs2_aq_ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF


   Roundoff = saprc99_mosaic_20bin_vbs2_aq_WLAMCH('E')


   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      Hmin = RCNTRL(1)
   ELSE
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL saprc99_mosaic_20bin_vbs2_aq_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL saprc99_mosaic_20bin_vbs2_aq_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL saprc99_mosaic_20bin_vbs2_aq_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2_dp
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL saprc99_mosaic_20bin_vbs2_aq_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0_dp
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL saprc99_mosaic_20bin_vbs2_aq_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1_dp
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL saprc99_mosaic_20bin_vbs2_aq_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9_dp
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL saprc99_mosaic_20bin_vbs2_aq_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.0_dp*Roundoff) &
         .OR. (RelTol(i) >= 1.0_dp) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL saprc99_mosaic_20bin_vbs2_aq_ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO



   SELECT CASE (Method)
     CASE (1)
       CALL saprc99_mosaic_20bin_vbs2_aq_Ros2(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (2)
       CALL saprc99_mosaic_20bin_vbs2_aq_Ros3(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (3)
       CALL saprc99_mosaic_20bin_vbs2_aq_Ros4(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (4)
       CALL saprc99_mosaic_20bin_vbs2_aq_Rodas3(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (5)
       CALL saprc99_mosaic_20bin_vbs2_aq_Rodas4(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(4)=', Method
       CALL saprc99_mosaic_20bin_vbs2_aq_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       RETURN
   END SELECT


   CALL saprc99_mosaic_20bin_vbs2_aq_ros_Integrator(Y,Tstart,Tend,Texit,      &
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



 SUBROUTINE  saprc99_mosaic_20bin_vbs2_aq_ros_ErrorMsg(Code,T,H,IERR)



   USE saprc99_mosaic_20bin_vbs2_aq_Precision

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

 END SUBROUTINE  saprc99_mosaic_20bin_vbs2_aq_ros_ErrorMsg


 SUBROUTINE  saprc99_mosaic_20bin_vbs2_aq_ros_Integrator (Y, Tstart, Tend, T,     &
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
      CALL saprc99_mosaic_20bin_vbs2_aq_ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1_dp*H) == T).OR.(H <= Roundoff) ) THEN  
      CALL saprc99_mosaic_20bin_vbs2_aq_ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF


   Hexit = H
   H = MIN(H,ABS(Tend-T))


   CALL saprc99_mosaic_20bin_vbs2_aq_FunTemplate(T,Y,Fcn0, RCONST, FIX, Nfun)
   IF( T == Tstart ) THEN
     CALL saprc99_mosaic_20bin_vbs2_aq_IRRFun( Y, FIX, RCONST, IRR_WRK )
   ENDIF


   IF (.NOT.Autonomous) THEN
      CALL saprc99_mosaic_20bin_vbs2_aq_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )
   END IF


   CALL saprc99_mosaic_20bin_vbs2_aq_JacTemplate(T,Y,Jac0, FIX, Njac, RCONST)


UntilAccepted: DO

   CALL saprc99_mosaic_20bin_vbs2_aq_ros_PrepareMatrix(H,Direction,ros_Gamma(1), &
          Jac0,Ghimj,Pivot,Singular, Ndec,  Nsng )
   IF (Singular) THEN 
       CALL saprc99_mosaic_20bin_vbs2_aq_ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF


Stage: DO istage = 1, ros_S

      
       ioffset = NVAR*(istage-1)

      
       IF ( istage == 1 ) THEN
         CALL saprc99_mosaic_20bin_vbs2_aq_WCOPY(NVAR,Fcn0,1,Fcn,1)
      
       ELSEIF ( ros_NewF(istage) ) THEN
         CALL saprc99_mosaic_20bin_vbs2_aq_WCOPY(NVAR,Y,1,Ynew,1)
         DO j = 1, istage-1
           CALL saprc99_mosaic_20bin_vbs2_aq_WAXPY(NVAR,ros_A((istage-1)*(istage-2)/2+j), &
            K(NVAR*(j-1)+1),1,Ynew,1)
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL saprc99_mosaic_20bin_vbs2_aq_FunTemplate(Tau,Ynew,Fcn, RCONST, FIX, Nfun)
       END IF 
       CALL saprc99_mosaic_20bin_vbs2_aq_WCOPY(NVAR,Fcn,1,K(ioffset+1),1)
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         CALL saprc99_mosaic_20bin_vbs2_aq_WAXPY(NVAR,HC,K(NVAR*(j-1)+1),1,K(ioffset+1),1)
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         CALL saprc99_mosaic_20bin_vbs2_aq_WAXPY(NVAR,HG,dFdT,1,K(ioffset+1),1)
       END IF
       CALL saprc99_mosaic_20bin_vbs2_aq_ros_Solve(Ghimj, Pivot, K(ioffset+1), Nsol)

   END DO Stage



   CALL saprc99_mosaic_20bin_vbs2_aq_WCOPY(NVAR,Y,1,Ynew,1)
   DO j=1,ros_S
         CALL saprc99_mosaic_20bin_vbs2_aq_WAXPY(NVAR,ros_M(j),K(NVAR*(j-1)+1),1,Ynew,1)
   END DO


   CALL saprc99_mosaic_20bin_vbs2_aq_WSCAL(NVAR,ZERO,Yerr,1)
   DO j=1,ros_S
        CALL saprc99_mosaic_20bin_vbs2_aq_WAXPY(NVAR,ros_E(j),K(NVAR*(j-1)+1),1,Yerr,1)
   END DO
   Err = saprc99_mosaic_20bin_vbs2_aq_ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )


   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac


   Nstp = Nstp+1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  
      Nacc = Nacc+1
      CALL saprc99_mosaic_20bin_vbs2_aq_WCOPY(NVAR,Ynew,1,Y,1)
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

  END SUBROUTINE  saprc99_mosaic_20bin_vbs2_aq_ros_Integrator



  REAL(kind=dp) FUNCTION  saprc99_mosaic_20bin_vbs2_aq_ros_ErrorNorm ( Y, Ynew, Yerr, &
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

    saprc99_mosaic_20bin_vbs2_aq_ros_ErrorNorm = Err

  END FUNCTION  saprc99_mosaic_20bin_vbs2_aq_ros_ErrorNorm



  SUBROUTINE saprc99_mosaic_20bin_vbs2_aq_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )



   IMPLICIT NONE


   REAL(kind=dp), INTENT(IN) :: T, Roundoff, Y(NVAR), Fcn0(NVAR)
   REAL(kind=dp), INTENT(IN) :: RCONST(NREACT), FIX(NFIX)

   REAL(kind=dp), INTENT(OUT) :: dFdT(NVAR)

   INTEGER, INTENT(INOUT) ::Nfun

   REAL(kind=dp) :: Delta
   REAL(kind=dp), PARAMETER :: ONE = 1.0_dp, DeltaMin = 1.0E-6_dp

   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL saprc99_mosaic_20bin_vbs2_aq_FunTemplate(T+Delta,Y,dFdT, RCONST, FIX, Nfun)
   CALL saprc99_mosaic_20bin_vbs2_aq_WAXPY(NVAR,(-ONE),Fcn0,1,dFdT,1)
   CALL saprc99_mosaic_20bin_vbs2_aq_WSCAL(NVAR,(ONE/Delta),dFdT,1)

  END SUBROUTINE  saprc99_mosaic_20bin_vbs2_aq_ros_FunTimeDeriv



  SUBROUTINE  saprc99_mosaic_20bin_vbs2_aq_ros_PrepareMatrix ( H, Direction, gam, &
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


     CALL saprc99_mosaic_20bin_vbs2_aq_WCOPY(LU_NONZERO,Jac0,1,Ghimj,1)
     CALL saprc99_mosaic_20bin_vbs2_aq_WSCAL(LU_NONZERO,(-ONE),Ghimj,1)
     ghinv = ONE/(Direction*H*gam)
     DO i=1,NVAR
       Ghimj(LU_DIAG(i)) = Ghimj(LU_DIAG(i))+ghinv
     END DO

     CALL saprc99_mosaic_20bin_vbs2_aq_ros_Decomp( Ghimj, Pivot, ising, Ndec )
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

  END SUBROUTINE  saprc99_mosaic_20bin_vbs2_aq_ros_PrepareMatrix



  SUBROUTINE  saprc99_mosaic_20bin_vbs2_aq_ros_Decomp( A, Pivot, ising, Ndec )



   IMPLICIT NONE

   REAL(kind=dp), INTENT(INOUT) :: A(LU_NONZERO)

   INTEGER, INTENT(OUT) :: Pivot(NVAR), ising
   INTEGER, INTENT(INOUT) :: Ndec 

   CALL saprc99_mosaic_20bin_vbs2_aq_KppDecomp ( A, ising )
   Pivot(1) = 1
   Ndec = Ndec + 1

  END SUBROUTINE  saprc99_mosaic_20bin_vbs2_aq_ros_Decomp



  SUBROUTINE  saprc99_mosaic_20bin_vbs2_aq_ros_Solve( A, Pivot, b, Nsol )



   IMPLICIT NONE

   REAL(kind=dp), INTENT(IN) :: A(LU_NONZERO)
   INTEGER, INTENT(IN) :: Pivot(NVAR)

   INTEGER, INTENT(INOUT) :: nsol 

   REAL(kind=dp), INTENT(INOUT) :: b(NVAR)


   CALL saprc99_mosaic_20bin_vbs2_aq_KppSolve( A, b )

   Nsol = Nsol+1

  END SUBROUTINE  saprc99_mosaic_20bin_vbs2_aq_ros_Solve




  SUBROUTINE  saprc99_mosaic_20bin_vbs2_aq_Ros2 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

 END SUBROUTINE  saprc99_mosaic_20bin_vbs2_aq_Ros2



  SUBROUTINE  saprc99_mosaic_20bin_vbs2_aq_Ros3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  saprc99_mosaic_20bin_vbs2_aq_Ros3





  SUBROUTINE  saprc99_mosaic_20bin_vbs2_aq_Ros4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  saprc99_mosaic_20bin_vbs2_aq_Ros4


  SUBROUTINE  saprc99_mosaic_20bin_vbs2_aq_Rodas3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  saprc99_mosaic_20bin_vbs2_aq_Rodas3


  SUBROUTINE  saprc99_mosaic_20bin_vbs2_aq_Rodas4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  saprc99_mosaic_20bin_vbs2_aq_Rodas4




END SUBROUTINE  saprc99_mosaic_20bin_vbs2_aq_Rosenbrock




SUBROUTINE  saprc99_mosaic_20bin_vbs2_aq_FunTemplate( T, Y, Ydot, RCONST, FIX, Nfun )




   USE saprc99_mosaic_20bin_vbs2_aq_Parameters




   REAL(kind=dp) :: T, Y(NVAR)
   REAL(kind=dp) :: RCONST(NREACT)
   REAL(kind=dp) :: FIX(NFIX)

   REAL(kind=dp) :: Ydot(NVAR)
   INTEGER :: Nfun









   CALL saprc99_mosaic_20bin_vbs2_aq_Fun( Y, FIX, RCONST, Ydot )


   Nfun = Nfun+1

END SUBROUTINE  saprc99_mosaic_20bin_vbs2_aq_FunTemplate



SUBROUTINE  saprc99_mosaic_20bin_vbs2_aq_JacTemplate( T, Y, Jcb, FIX, Njac, RCONST )




 USE saprc99_mosaic_20bin_vbs2_aq_Parameters
 
 USE saprc99_mosaic_20bin_vbs2_aq_Jacobian



    REAL(kind=dp) :: T, Y(NVAR)
    REAL(kind=dp) :: FIX(NFIX)
    REAL(kind=dp) :: RCONST(NREACT)

    INTEGER :: Njac


    REAL(kind=dp) :: Jcb(LU_NONZERO)

    REAL(kind=dp) :: Told





    CALL saprc99_mosaic_20bin_vbs2_aq_Jac_SP( Y, FIX, RCONST, Jcb )


    Njac = Njac+1

END SUBROUTINE  saprc99_mosaic_20bin_vbs2_aq_JacTemplate

















SUBROUTINE saprc99_mosaic_20bin_vbs2_aq_Fun ( V, F, RCT, Vdot )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: Vdot(NVAR)




  REAL(kind=dp) :: A(NREACT)


  A(1) = RCT(1)*V(93)
  A(2) = RCT(2)*V(85)*F(2)
  A(3) = RCT(3)*V(85)*V(89)
  A(4) = RCT(4)*V(85)*V(96)*F(2)
  A(5) = RCT(5)*V(85)*V(93)
  A(6) = RCT(6)*V(85)*V(93)
  A(7) = RCT(7)*V(89)*V(96)
  A(8) = RCT(8)*V(89)*V(93)
  A(9) = RCT(9)*V(94)*V(96)
  A(10) = RCT(10)*V(96)*V(96)*F(2)
  A(11) = RCT(11)*V(93)*V(94)
  A(12) = RCT(12)*V(41)
  A(13) = RCT(13)*V(41)*F(1)
  A(14) = RCT(14)*V(93)*V(94)
  A(15) = RCT(15)*V(94)
  A(16) = RCT(16)*V(94)
  A(17) = RCT(17)*V(89)
  A(18) = RCT(18)*V(89)
  A(19) = RCT(19)*V(26)*F(1)
  A(20) = RCT(20)*V(26)*F(2)
  A(21) = RCT(21)*V(95)*V(96)
  A(22) = RCT(22)*V(44)
  A(23) = RCT(23)*V(44)
  A(24) = RCT(24)*V(44)*V(95)
  A(25) = RCT(25)*V(93)*V(95)
  A(26) = RCT(26)*V(94)*V(95)
  A(27) = RCT(27)*V(67)*V(95)
  A(28) = RCT(28)*V(67)
  A(29) = RCT(29)*V(66)*V(95)
  A(30) = RCT(30)*V(89)*V(95)
  A(31) = RCT(31)*V(96)*V(97)
  A(32) = RCT(32)*V(93)*V(97)
  A(33) = RCT(33)*V(51)
  A(34) = RCT(34)*V(51)
  A(35) = RCT(35)*V(51)*V(95)
  A(36) = RCT(36)*V(89)*V(97)
  A(37) = RCT(37)*V(97)*V(97)
  A(38) = RCT(38)*V(97)*V(97)*F(1)
  A(39) = RCT(39)*V(94)*V(97)
  A(40) = RCT(40)*V(94)*V(94)
  A(41) = RCT(41)*V(37)
  A(42) = RCT(42)*V(37)*V(95)
  A(43) = RCT(43)*V(95)*V(97)
  A(44) = RCT(44)*V(42)*V(95)
  A(45) = RCT(45)*V(95)*F(2)
  A(46) = RCT(46)*V(96)*V(98)
  A(47) = RCT(47)*V(97)*V(98)
  A(48) = RCT(48)*V(94)*V(98)
  A(49) = RCT(49)*V(98)*V(98)
  A(50) = RCT(50)*V(98)*V(98)
  A(51) = RCT(51)*V(96)*V(99)
  A(52) = RCT(52)*V(97)*V(99)
  A(53) = RCT(53)*V(94)*V(99)
  A(54) = RCT(54)*V(98)*V(99)
  A(55) = RCT(55)*V(99)*V(99)
  A(56) = RCT(56)*V(74)*V(96)
  A(57) = RCT(57)*V(74)*V(97)
  A(58) = RCT(58)*V(74)*V(94)
  A(59) = RCT(59)*V(74)*V(98)
  A(60) = RCT(60)*V(74)*V(99)
  A(61) = RCT(61)*V(74)*V(74)
  A(62) = RCT(62)*V(96)*V(101)
  A(63) = RCT(63)*V(97)*V(101)
  A(64) = RCT(64)*V(98)*V(101)
  A(65) = RCT(65)*V(94)*V(101)
  A(66) = RCT(66)*V(99)*V(101)
  A(67) = RCT(67)*V(74)*V(101)
  A(68) = RCT(68)*V(101)*V(101)
  A(69) = RCT(69)*V(90)*V(93)
  A(70) = RCT(70)*V(33)
  A(71) = RCT(71)*V(90)*V(96)
  A(72) = RCT(72)*V(90)*V(97)
  A(73) = RCT(73)*V(90)*V(94)
  A(74) = RCT(74)*V(90)*V(98)
  A(75) = RCT(75)*V(90)*V(99)
  A(76) = RCT(76)*V(74)*V(90)
  A(77) = RCT(77)*V(90)*V(101)
  A(78) = RCT(78)*V(90)*V(90)
  A(79) = RCT(79)*V(92)*V(93)
  A(80) = RCT(80)*V(34)
  A(81) = RCT(81)*V(92)*V(96)
  A(82) = RCT(82)*V(92)*V(97)
  A(83) = RCT(83)*V(92)*V(94)
  A(84) = RCT(84)*V(92)*V(98)
  A(85) = RCT(85)*V(92)*V(99)
  A(86) = RCT(86)*V(74)*V(92)
  A(87) = RCT(87)*V(92)*V(101)
  A(88) = RCT(88)*V(90)*V(92)
  A(89) = RCT(89)*V(92)*V(92)
  A(90) = RCT(90)*V(91)*V(93)
  A(91) = RCT(91)*V(35)
  A(92) = RCT(92)*V(91)*V(96)
  A(93) = RCT(93)*V(91)*V(97)
  A(94) = RCT(94)*V(91)*V(94)
  A(95) = RCT(95)*V(91)*V(98)
  A(96) = RCT(96)*V(91)*V(99)
  A(97) = RCT(97)*V(74)*V(91)
  A(98) = RCT(98)*V(91)*V(101)
  A(99) = RCT(99)*V(90)*V(91)
  A(100) = RCT(100)*V(91)*V(92)
  A(101) = RCT(101)*V(91)*V(91)
  A(102) = RCT(102)*V(93)*V(100)
  A(103) = RCT(103)*V(36)
  A(104) = RCT(104)*V(96)*V(100)
  A(105) = RCT(105)*V(97)*V(100)
  A(106) = RCT(106)*V(94)*V(100)
  A(107) = RCT(107)*V(98)*V(100)
  A(108) = RCT(108)*V(99)*V(100)
  A(109) = RCT(109)*V(74)*V(100)
  A(110) = RCT(110)*V(100)*V(101)
  A(111) = RCT(111)*V(90)*V(100)
  A(112) = RCT(112)*V(92)*V(100)
  A(113) = RCT(113)*V(91)*V(100)
  A(114) = RCT(114)*V(100)*V(100)
  A(115) = RCT(115)*V(46)*V(93)
  A(116) = RCT(116)*V(46)
  A(117) = RCT(117)*V(71)*V(93)
  A(118) = RCT(118)*V(71)*V(97)
  A(119) = RCT(119)*V(71)
  A(120) = RCT(120)*V(50)*V(93)
  A(121) = RCT(121)*V(50)*V(97)
  A(122) = RCT(122)*V(50)
  A(123) = RCT(123)*V(83)
  A(124) = RCT(124)*V(83)
  A(125) = RCT(125)*V(83)*V(95)
  A(126) = RCT(126)*V(83)*V(97)
  A(127) = RCT(127)*V(49)
  A(128) = RCT(128)*V(49)*V(96)
  A(129) = RCT(129)*V(83)*V(94)
  A(130) = RCT(130)*V(82)*V(95)
  A(131) = RCT(131)*V(82)
  A(132) = RCT(132)*V(82)*V(94)
  A(133) = RCT(133)*V(86)*V(95)
  A(134) = RCT(134)*V(86)
  A(135) = RCT(135)*V(86)*V(94)
  A(136) = RCT(136)*V(69)*V(95)
  A(137) = RCT(137)*V(69)
  A(138) = RCT(138)*V(87)*V(95)
  A(139) = RCT(139)*V(87)
  A(140) = RCT(140)*V(53)*V(95)
  A(141) = RCT(141)*V(40)*V(95)
  A(142) = RCT(142)*V(48)*V(95)
  A(143) = RCT(143)*V(48)
  A(144) = RCT(144)*V(62)*V(95)
  A(145) = RCT(145)*V(62)
  A(146) = RCT(146)*V(78)
  A(147) = RCT(147)*V(78)
  A(148) = RCT(148)*V(78)*V(95)
  A(149) = RCT(149)*V(78)*V(94)
  A(150) = RCT(150)*V(65)
  A(151) = RCT(151)*V(65)*V(95)
  A(152) = RCT(152)*V(65)*V(94)
  A(153) = RCT(153)*V(39)
  A(154) = RCT(154)*V(64)*V(95)
  A(155) = RCT(155)*V(64)*V(94)
  A(156) = RCT(156)*V(58)*V(95)
  A(157) = RCT(157)*V(58)*V(94)
  A(158) = RCT(158)*V(61)*V(94)
  A(159) = RCT(159)*V(63)*V(95)
  A(160) = RCT(160)*V(63)
  A(161) = RCT(161)*V(63)*V(94)
  A(162) = RCT(162)*V(75)*V(95)
  A(163) = RCT(163)*V(75)*V(89)
  A(164) = RCT(164)*V(75)*V(94)
  A(165) = RCT(165)*V(75)*V(85)
  A(166) = RCT(166)*V(75)
  A(167) = RCT(167)*V(81)*V(95)
  A(168) = RCT(168)*V(81)*V(89)
  A(169) = RCT(169)*V(81)*V(85)
  A(170) = RCT(170)*V(81)
  A(171) = RCT(171)*V(79)*V(95)
  A(172) = RCT(172)*V(79)*V(89)
  A(173) = RCT(173)*V(79)*V(94)
  A(174) = RCT(174)*V(79)
  A(175) = RCT(175)*V(88)*V(95)
  A(176) = RCT(176)*V(88)
  A(177) = RCT(177)*V(84)*V(95)
  A(178) = RCT(178)*V(84)
  A(179) = RCT(179)*V(59)*V(95)
  A(180) = RCT(180)*V(59)*V(89)
  A(181) = RCT(181)*V(56)*V(95)
  A(182) = RCT(182)*V(56)
  A(183) = RCT(183)*V(57)*V(95)
  A(184) = RCT(184)*V(57)
  A(185) = RCT(185)*V(31)*V(95)
  A(186) = RCT(186)*V(68)*V(95)
  A(187) = RCT(187)*V(68)*V(89)
  A(188) = RCT(188)*V(68)*V(94)
  A(189) = RCT(189)*V(68)*V(85)
  A(190) = RCT(190)*V(73)*V(95)
  A(191) = RCT(191)*V(73)*V(89)
  A(192) = RCT(192)*V(73)*V(94)
  A(193) = RCT(193)*V(73)*V(85)
  A(194) = RCT(194)*V(77)*V(95)
  A(195) = RCT(195)*V(77)*V(89)
  A(196) = RCT(196)*V(77)*V(94)
  A(197) = RCT(197)*V(77)*V(85)
  A(198) = RCT(198)*V(76)*V(95)
  A(199) = RCT(199)*V(76)*V(89)
  A(200) = RCT(200)*V(76)*V(94)
  A(201) = RCT(201)*V(76)*V(85)
  A(202) = RCT(202)*V(32)*V(95)
  A(203) = RCT(203)*V(38)*V(95)
  A(204) = RCT(204)*V(60)*V(95)
  A(205) = RCT(205)*V(45)*V(95)
  A(206) = RCT(206)*V(54)*V(95)
  A(207) = RCT(207)*V(47)*V(95)
  A(208) = RCT(208)*V(55)*V(95)
  A(209) = RCT(209)*V(52)*V(95)
  A(210) = RCT(210)*V(72)*V(95)
  A(211) = RCT(211)*V(72)*V(89)
  A(212) = RCT(212)*V(72)*V(94)
  A(213) = RCT(213)*V(72)*V(85)
  A(214) = RCT(214)*V(80)*V(95)
  A(215) = RCT(215)*V(80)*V(89)
  A(216) = RCT(216)*V(80)*V(94)
  A(217) = RCT(217)*V(80)*V(85)
  A(218) = RCT(218)*V(60)*V(89)
  A(219) = RCT(219)*V(70)*V(95)
  A(220) = RCT(220)*V(70)*V(89)
  A(221) = RCT(221)*V(70)*V(94)
  A(222) = RCT(222)*V(70)*V(85)
  A(223) = RCT(223)*V(42)
  A(224) = RCT(224)*V(97)
  A(225) = RCT(225)*V(42)
  A(226) = RCT(226)*V(1)
  A(227) = RCT(227)*V(67)
  A(228) = RCT(228)*V(37)
  A(229) = RCT(229)*V(2)
  A(230) = RCT(230)*V(54)*V(95)
  A(231) = RCT(231)*V(47)*V(95)
  A(232) = RCT(232)*V(72)*V(95)
  A(233) = RCT(233)*V(80)*V(95)
  A(234) = RCT(234)*V(55)*V(95)
  A(235) = RCT(235)*V(52)*V(95)
  A(236) = RCT(236)*V(73)*V(95)
  A(237) = RCT(237)*V(73)*V(89)
  A(238) = RCT(238)*V(73)*V(94)
  A(239) = RCT(239)*V(77)*V(95)
  A(240) = RCT(240)*V(77)*V(89)
  A(241) = RCT(241)*V(77)*V(94)
  A(242) = RCT(242)*V(76)*V(95)
  A(243) = RCT(243)*V(76)*V(89)
  A(244) = RCT(244)*V(76)*V(94)
  A(245) = RCT(245)*V(3)*V(95)
  A(246) = RCT(246)*V(4)*V(95)
  A(247) = RCT(247)*V(28)*V(95)
  A(248) = RCT(248)*V(22)*V(95)
  A(249) = RCT(249)*V(5)*V(95)
  A(250) = RCT(250)*V(6)*V(95)
  A(251) = RCT(251)*V(30)*V(95)
  A(252) = RCT(252)*V(24)*V(95)
  A(253) = RCT(253)*V(27)*V(95)
  A(254) = RCT(254)*V(23)*V(95)
  A(255) = RCT(255)*V(29)*V(95)
  A(256) = RCT(256)*V(25)*V(95)
  A(257) = RCT(257)*V(43)*V(95)
  A(258) = RCT(258)*V(43)*V(95)
  A(259) = RCT(259)*V(43)*V(94)


  Vdot(1) = A(44)+A(223)-A(226)
  Vdot(2) = 0.5*A(218)+0.135*A(220)-A(229)
  Vdot(3) = 0
  Vdot(4) = 0
  Vdot(5) = 0
  Vdot(6) = 0
  Vdot(7) = A(128)+0.333*A(163)+0.351*A(168)+0.1*A(172)+0.37*A(187)+0.204*A(191)+0.103*A(195)+0.103*A(199)+0.297*A(204)&
              &+0.185*A(211)+0.073*A(215)+0.185*A(220)
  Vdot(8) = 0.25*A(72)+A(74)+A(75)+A(77)+0.05*A(211)+0.129*A(215)+0.17*A(220)
  Vdot(9) = 0.25*A(82)+A(84)+A(85)+A(87)+0.25*A(93)+A(95)+A(96)+A(98)+0.25*A(105)+A(107)+A(108)+2*A(110)+0.372*A(172)&
              &+0.15*A(191)+0.189*A(195)+0.189*A(199)+0.119*A(211)+0.247*A(215)
  Vdot(10) = 0.75*A(72)
  Vdot(11) = 0.75*A(82)+0.75*A(93)+0.75*A(105)
  Vdot(12) = 2*A(120)+A(221)
  Vdot(13) = 6*A(120)+7*A(160)+0.048*A(219)+0.07*A(220)+2.693*A(221)+0.55*A(222)
  Vdot(14) = A(46)+A(48)+A(51)+A(53)+A(56)+A(58)+A(62)+A(65)
  Vdot(15) = A(47)+A(49)+A(50)+A(52)+A(54)+A(55)+A(57)+A(59)+A(60)+A(61)+A(63)+A(64)+A(66)+A(67)+A(68)
  Vdot(16) = A(230)+A(231)+A(232)+A(233)+A(234)+A(235)
  Vdot(17) = A(239)
  Vdot(18) = A(240)
  Vdot(19) = A(236)+A(237)+0.0327*A(238)+0.0545*A(241)+A(242)+A(243)+0.0816*A(244)
  Vdot(20) = A(21)+A(24)+A(25)+A(26)+A(27)+A(29)+A(30)+A(43)+A(44)+A(45)+A(125)+A(130)+A(133)+A(136)+A(138)+A(140)&
               &+A(141)+A(142)+A(144)+A(148)+A(151)+A(154)+A(156)+A(159)+A(162)+A(167)+A(171)+A(175)+A(177)+A(179)+A(181)&
               &+A(183)+A(185)+A(186)+A(190)+A(194)+A(198)+A(202)+A(203)+A(204)+A(205)+A(206)+A(207)+A(208)+A(209)+A(210)&
               &+A(214)+A(219)
  Vdot(21) = A(245)+A(246)+A(247)+A(248)+A(249)+A(250)+A(251)+A(252)+A(253)+A(254)+A(255)+A(256)
  Vdot(22) = 0.5*A(253)+A(254)
  Vdot(23) = -A(254)
  Vdot(24) = 0.5*A(255)+A(256)
  Vdot(25) = -A(256)
  Vdot(26) = A(18)-A(19)-A(20)
  Vdot(27) = -A(253)
  Vdot(28) = A(253)
  Vdot(29) = -A(255)
  Vdot(30) = A(255)
  Vdot(31) = -A(185)
  Vdot(32) = -A(202)
  Vdot(33) = A(69)-A(70)
  Vdot(34) = A(79)-A(80)
  Vdot(35) = A(90)-A(91)
  Vdot(36) = A(102)-A(103)
  Vdot(37) = A(37)+A(38)-A(41)-A(42)-A(228)
  Vdot(38) = -A(203)
  Vdot(39) = -A(153)+0.031*A(195)+0.031*A(199)+0.087*A(209)
  Vdot(40) = -A(141)
  Vdot(41) = A(11)-A(12)-A(13)
  Vdot(42) = -A(44)-A(223)-A(225)+A(257)+0.5*A(258)+A(259)
  Vdot(43) = -A(257)-A(258)-A(259)
  Vdot(44) = A(21)-A(22)-A(23)-A(24)
  Vdot(45) = -A(205)
  Vdot(46) = -A(115)-A(116)+0.236*A(205)
  Vdot(47) = -A(207)
  Vdot(48) = A(47)-A(142)-A(143)
  Vdot(49) = A(126)-A(127)-A(128)
  Vdot(50) = -A(120)-A(121)-A(122)+A(158)
  Vdot(51) = A(32)-A(33)-A(34)-A(35)
  Vdot(52) = -A(209)
  Vdot(53) = A(49)+0.25*A(54)+0.25*A(64)-A(140)
  Vdot(54) = -A(206)
  Vdot(55) = -A(208)
  Vdot(56) = -A(181)-A(182)+0.108*A(208)+0.099*A(209)
  Vdot(57) = -A(183)-A(184)+0.051*A(208)+0.093*A(209)
  Vdot(58) = -A(156)-A(157)+0.207*A(208)+0.187*A(209)
  Vdot(59) = -A(179)-A(180)+0.491*A(208)+0.561*A(209)
  Vdot(60) = -A(204)-A(218)
  Vdot(61) = A(117)+A(121)+A(122)-A(158)
  Vdot(62) = A(52)+A(63)-A(144)-A(145)
  Vdot(63) = -A(159)-A(160)-A(161)+0.059*A(208)+0.05*A(209)+0.061*A(214)+0.042*A(215)+0.015*A(216)
  Vdot(64) = A(118)+A(119)-A(154)-A(155)+0.017*A(208)
  Vdot(65) = -A(150)-A(151)-A(152)+0.23*A(156)+0.084*A(162)+0.9*A(163)+0.3*A(167)+0.95*A(168)+0.174*A(171)+0.742*A(172)&
               &+0.008*A(173)+0.5*A(182)+0.5*A(184)+0.119*A(208)+0.287*A(209)
  Vdot(66) = -A(29)+A(123)+A(124)+A(125)+A(129)+A(131)+0.034*A(133)+A(134)+2*A(146)+A(147)+1.26*A(148)+1.26*A(149)&
               &+A(150)+A(151)+A(152)+0.416*A(162)+0.45*A(163)+0.5*A(164)+0.67*A(166)+0.475*A(168)+0.7*A(170)+0.336*A(171)&
               &+0.498*A(172)+0.572*A(173)+1.233*A(174)+A(179)+1.5*A(180)+A(182)+A(184)+0.5*A(187)+0.491*A(189)+0.275*A(191)&
               &+0.157*A(195)+0.157*A(199)+0.393*A(204)+0.002*A(206)+0.345*A(211)+0.265*A(215)+0.012*A(217)+1.5*A(218)+0.51&
               &*A(220)
  Vdot(67) = 2*A(13)+A(25)-A(27)-A(28)+0.2*A(39)+A(129)+A(132)+A(135)+A(149)+A(152)+A(155)+A(157)+A(158)+A(161)+0.5&
               &*A(164)+0.15*A(173)-A(227)+A(259)
  Vdot(68) = -A(186)-A(187)-A(188)-A(189)
  Vdot(69) = A(116)-A(136)-A(137)+0.006*A(177)+0.02*A(178)+0.13*A(195)+0.13*A(199)+0.704*A(203)+0.024*A(205)+0.452&
               &*A(206)+0.072*A(207)+0.005*A(210)+0.001*A(211)+0.024*A(212)+0.127*A(214)+0.045*A(215)+0.102*A(216)
  Vdot(70) = -A(219)-A(220)-A(221)-A(222)
  Vdot(71) = A(92)+A(94)+A(99)+A(100)+2*A(101)+A(113)-A(117)-A(118)-A(119)+0.24*A(154)+A(155)+0.24*A(156)+A(157)
  Vdot(72) = -A(210)-A(211)-A(212)-A(213)
  Vdot(73) = -A(190)-A(191)-A(192)-A(193)
  Vdot(74) = -A(56)-A(57)-A(58)-A(59)-A(60)-A(67)-A(76)-A(86)+A(92)+A(94)-A(97)+A(99)+A(100)+2*A(101)-A(109)+A(113)&
               &+A(136)+0.616*A(138)+0.675*A(167)+0.515*A(176)+0.596*A(177)+0.152*A(178)+A(181)+A(182)+A(183)+A(184)+0.079&
               &*A(190)+0.126*A(191)+0.187*A(192)+0.24*A(193)+0.5*A(194)+0.729*A(195)+0.75*A(196)+0.5*A(198)+0.729*A(199)&
               &+0.75*A(200)+0.559*A(205)+0.936*A(206)+0.948*A(207)+0.205*A(210)+0.488*A(212)+0.001*A(214)+0.137*A(215)&
               &+0.711*A(216)
  Vdot(75) = -A(162)-A(163)-A(164)-A(165)-A(166)+0.23*A(190)+0.39*A(191)+0.025*A(214)+0.026*A(215)+0.012*A(217)
  Vdot(76) = -A(198)-A(199)-A(200)-A(201)
  Vdot(77) = -A(194)-A(195)-A(196)-A(197)
  Vdot(78) = -A(146)-A(147)-A(148)-A(149)+0.23*A(154)+0.15*A(171)+0.023*A(172)+A(180)+0.5*A(182)+0.5*A(184)+0.009*A(189)&
               &+0.001*A(195)+0.001*A(199)+0.607*A(204)+0.118*A(208)+0.097*A(209)
  Vdot(79) = -A(171)-A(172)-A(173)-A(174)+0.357*A(190)+0.936*A(192)+0.025*A(214)
  Vdot(80) = -A(214)-A(215)-A(216)-A(217)
  Vdot(81) = -A(167)-A(168)-A(169)-A(170)+0.32*A(190)+0.16*A(191)+0.019*A(215)+0.048*A(216)
  Vdot(82) = A(81)+A(83)+A(88)+2*A(89)+A(100)+A(112)-A(130)-A(131)-A(132)+0.034*A(133)+A(134)+0.482*A(138)+A(139)+0.96&
               &*A(141)+0.129*A(171)+0.047*A(172)+0.467*A(174)+0.084*A(175)+0.246*A(176)+0.439*A(177)+0.431*A(178)+0.195&
               &*A(186)+0.25*A(189)+A(202)+0.445*A(205)+0.455*A(206)+0.099*A(207)+0.294*A(210)+0.154*A(211)+0.009*A(212)&
               &+0.732*A(214)+0.456*A(215)+0.507*A(216)+0.984*A(219)+0.5*A(220)
  Vdot(83) = A(46)+A(48)+A(49)+2*A(50)+0.75*A(54)+0.75*A(64)+A(74)+A(84)+A(95)+A(104)+A(106)+A(107)+A(111)+A(112)+A(113)&
               &+2*A(114)-A(123)-A(124)-A(125)-A(126)+A(127)-A(129)+A(136)+0.115*A(138)+A(140)+0.081*A(141)+0.35*A(142)&
               &+A(143)+A(147)+0.084*A(162)+0.2*A(163)+0.67*A(166)+0.3*A(167)+0.1*A(168)+0.055*A(171)+0.125*A(172)+0.227&
               &*A(173)+0.3*A(174)+0.213*A(175)+0.506*A(176)+0.01*A(177)+0.134*A(178)+1.61*A(186)+A(187)+0.191*A(189)+0.624&
               &*A(190)+0.592*A(191)+0.24*A(193)+0.276*A(194)+0.235*A(195)+0.276*A(198)+0.235*A(199)+0.096*A(204)+0.026&
               &*A(205)+0.024*A(206)+0.026*A(207)+0.732*A(210)+0.5*A(211)+0.244*A(214)+0.269*A(215)+0.079*A(216)+0.984&
               &*A(219)+0.5*A(220)
  Vdot(84) = A(62)+A(115)+0.572*A(173)-0.69*A(177)-A(178)+0.276*A(196)+0.276*A(200)+0.511*A(212)+0.321*A(216)
  Vdot(85) = A(1)-A(2)-A(3)-A(4)-A(5)-A(6)+A(16)+A(17)+A(20)-A(165)-A(169)-A(189)-A(193)-A(197)-A(201)-A(213)-A(217)&
               &-A(222)
  Vdot(86) = -A(133)-A(134)-A(135)+0.37*A(138)+A(144)+A(145)+A(165)+0.675*A(167)+0.45*A(169)+0.013*A(171)+0.218*A(173)&
               &+0.558*A(175)+0.71*A(176)+0.213*A(177)+0.147*A(178)+A(179)+A(181)+A(183)+A(188)+0.474*A(194)+0.205*A(195)&
               &+0.474*A(196)+0.147*A(197)+0.474*A(198)+0.205*A(199)+0.474*A(200)+0.147*A(201)+0.261*A(203)+0.122*A(205)&
               &+0.244*A(206)+0.204*A(207)+0.497*A(210)+0.363*A(211)+0.037*A(212)+0.45*A(213)+0.511*A(214)+0.305*A(215)&
               &+0.151*A(216)+0.069*A(217)+0.45*A(222)
  Vdot(87) = 0.5*A(64)+A(65)+0.5*A(66)+A(68)-A(138)-A(139)+0.416*A(162)+0.55*A(169)+0.15*A(171)+0.21*A(172)+0.233*A(174)&
               &+0.115*A(175)+0.177*A(177)+0.243*A(178)+0.332*A(205)+0.11*A(206)+0.089*A(207)+0.437*A(213)+0.072*A(214)&
               &+0.026*A(215)+0.001*A(216)+0.659*A(217)+0.55*A(222)
  Vdot(88) = 0.5*A(64)+0.5*A(66)+A(68)+A(77)+A(87)+A(98)+0.7*A(170)+0.332*A(171)-0.671*A(175)-A(176)+0.048*A(177)+0.435&
               &*A(178)+0.1*A(191)+0.75*A(193)+0.276*A(194)+0.276*A(195)+0.853*A(197)+0.276*A(198)+0.276*A(199)+0.853*A(201)&
               &+0.125*A(206)+0.417*A(207)+0.055*A(208)+0.119*A(210)+0.215*A(211)+0.113*A(213)+0.043*A(215)+0.259*A(217)
  Vdot(89) = A(2)-A(3)-A(7)-A(8)-A(17)-A(18)-A(30)-A(36)+0.25*A(72)+0.25*A(82)+0.25*A(93)+0.25*A(105)-A(163)-A(168)&
               &-A(172)-A(180)-A(187)-A(191)-A(195)-A(199)-A(211)-A(215)-A(218)-A(220)
  Vdot(90) = -A(69)+A(70)-A(71)-A(72)-A(73)-A(74)-A(75)-A(77)-2*A(78)-A(88)-A(99)+A(104)+A(106)+A(112)+A(113)+2*A(114)&
               &+A(130)+A(132)+A(136)+A(137)+0.492*A(138)+A(139)+A(150)+A(151)+A(152)+2*A(153)+0.67*A(166)+0.675*A(167)&
               &+0.467*A(174)+0.029*A(175)+0.667*A(176)+A(181)+0.5*A(182)+A(183)+0.5*A(184)+0.123*A(195)+0.123*A(199)+0.011&
               &*A(206)+0.137*A(215)
  Vdot(91) = -A(90)+A(91)-A(92)-A(93)-A(94)-A(95)-A(96)-A(98)-A(99)-A(100)-2*A(101)-A(113)+A(159)+A(161)
  Vdot(92) = -A(79)+A(80)-A(81)-A(82)-A(83)-A(84)-A(85)-A(87)-A(88)-2*A(89)-A(100)-A(112)+0.965*A(133)+A(135)+0.096&
               &*A(138)+0.37*A(148)+0.37*A(149)+0.1*A(163)+0.05*A(168)+0.048*A(172)+0.3*A(174)+0.049*A(175)+0.333*A(176)&
               &+0.201*A(195)+0.201*A(199)+0.006*A(215)
  Vdot(93) = -A(1)+A(4)-A(5)-A(6)+A(7)-A(8)+2*A(9)+2*A(10)-A(11)+A(12)+A(16)+A(23)+A(24)-A(25)+A(26)+A(28)+A(31)-A(32)&
               &+A(33)+0.61*A(34)+A(35)+0.8*A(39)+2*A(40)+A(46)+A(48)+A(51)+A(53)+A(56)+A(58)+A(65)-A(69)+A(70)+A(71)+A(73)&
               &-A(79)+A(80)+A(81)+A(83)-A(90)+A(91)+A(92)+A(94)-A(102)+A(103)+A(104)+A(106)-A(115)-A(117)-A(120)+A(128)&
               &+0.338*A(177)+A(178)+0.187*A(192)+0.474*A(196)+0.474*A(200)+0.391*A(216)
  Vdot(94) = A(6)+A(8)-A(9)-A(11)+A(12)-A(14)-A(15)-A(16)-A(26)+A(27)+0.39*A(34)-A(39)-2*A(40)-A(48)-A(53)-A(58)-A(65)&
               &-A(73)-A(83)-A(94)-A(106)-A(129)-A(132)-A(135)-A(149)-A(152)-A(155)-A(157)-A(158)-A(161)-A(164)-A(173)&
               &-A(188)-A(192)-A(196)-A(200)-A(212)-A(216)-A(221)-A(259)
  Vdot(95) = 2*A(19)-A(21)+A(22)-A(24)-A(25)-A(26)-A(27)+A(28)-A(29)-A(30)+A(31)+0.39*A(34)-A(35)+A(36)+0.8*A(39)+2&
               &*A(41)-A(42)-A(43)-A(44)-A(45)-A(125)-A(130)-A(133)-A(136)-A(138)-A(140)-A(141)-0.65*A(142)+A(143)-0.34&
               &*A(144)+A(145)-A(148)-A(151)-A(154)-A(156)-A(159)-A(162)+0.208*A(163)+0.33*A(166)-A(167)+0.164*A(168)-A(171)&
               &+0.285*A(172)-A(175)-A(177)-A(179)+0.5*A(180)-A(181)-A(183)-A(185)-A(186)+0.12*A(187)-A(190)+0.266*A(191)&
               &-A(194)+0.567*A(195)-A(198)+0.567*A(199)-A(202)-A(203)-0.397*A(204)-A(205)-A(206)-A(207)-A(208)-A(209)&
               &-A(210)+0.155*A(211)-A(214)+0.378*A(215)+0.5*A(218)-A(219)+0.32*A(220)-A(245)-A(247)-A(249)-A(251)-A(253)&
               &-A(255)-A(257)-A(258)
  Vdot(96) = A(1)-A(4)+A(5)-A(7)-A(9)-2*A(10)+A(14)+A(15)-A(21)+A(22)-A(31)-A(46)-A(51)-A(56)-A(62)-A(71)-A(81)-A(92)&
               &-A(104)-A(128)
  Vdot(97) = A(23)+A(26)+A(29)+A(30)-A(31)-A(32)+A(33)+0.61*A(34)-A(36)-2*A(37)-2*A(38)-A(39)+A(42)-A(43)+A(44)+A(45)&
               &+A(46)-A(47)+A(48)+2*A(50)+A(51)-A(52)+A(53)+A(54)+A(55)-A(63)+A(64)+A(65)+A(66)+A(68)-A(72)-A(82)-A(93)&
               &-A(105)-A(118)-A(121)+2*A(123)+A(125)-A(126)+A(127)+A(128)+A(129)+A(131)+A(134)+A(140)+0.95*A(141)+A(143)&
               &+A(145)+2*A(146)+0.63*A(148)+0.63*A(149)+A(150)+0.008*A(163)+0.34*A(166)+0.064*A(168)+0.4*A(172)+1.233&
               &*A(174)+0.379*A(175)+0.113*A(177)+0.341*A(178)+1.5*A(180)+0.5*A(182)+0.5*A(184)+0.12*A(187)+0.5*A(189)+0.033&
               &*A(195)+0.033*A(199)+0.297*A(204)+0.224*A(208)+0.187*A(209)+0.056*A(211)+0.003*A(215)+0.013*A(217)+1.5&
               &*A(218)+0.06*A(220)-A(224)+0.5*A(258)
  Vdot(98) = -A(46)-A(47)-A(48)-2*A(49)-2*A(50)-A(54)-A(64)+A(71)+A(73)-A(74)+2*A(78)-A(84)+A(88)-A(95)+A(99)-A(107)&
               &+A(111)+A(116)+A(131)+A(137)+0.65*A(142)+0.3*A(170)+A(185)+0.3*A(189)+0.25*A(193)+0.011*A(206)+0.076*A(211)&
               &+0.197*A(215)+0.03*A(216)+0.26*A(220)
  Vdot(99) = -A(51)-A(52)-A(53)-A(54)-2*A(55)-A(66)-A(75)+A(81)+A(83)-A(85)+A(88)+2*A(89)-A(96)+A(100)-A(108)+A(112)&
               &+0.034*A(133)+A(134)+0.37*A(138)+A(139)+0.05*A(141)+0.34*A(144)+0.76*A(154)+0.76*A(156)+0.5*A(162)+0.1&
               &*A(163)+0.5*A(164)+0.33*A(166)+0.3*A(167)+0.05*A(168)+0.67*A(171)+0.048*A(172)+0.799*A(173)+0.473*A(175)&
               &+0.96*A(176)+0.376*A(177)+0.564*A(178)+A(179)+A(182)+A(184)+A(186)+A(188)+0.2*A(189)+0.907*A(190)+0.066&
               &*A(191)+0.749*A(192)+0.75*A(194)+0.031*A(195)+0.276*A(196)+0.75*A(198)+0.031*A(199)+0.276*A(200)+A(202)&
               &+0.965*A(203)+0.1*A(204)+0.695*A(205)+0.835*A(206)+0.653*A(207)+0.765*A(208)+0.804*A(209)+0.91*A(210)+0.022&
               &*A(211)+0.824*A(212)+0.918*A(214)+0.033*A(215)+0.442*A(216)+0.012*A(217)+0.984*A(219)+0.949*A(221)
  Vdot(100) = -A(102)+A(103)-A(104)-A(105)-A(106)-A(107)-A(108)-A(110)-A(111)-A(112)-A(113)-2*A(114)+0.5*A(162)+0.5&
                &*A(164)+0.33*A(166)+0.3*A(170)+0.289*A(171)+0.15*A(173)+0.192*A(191)+0.24*A(193)
  Vdot(101) = -A(62)-A(63)-A(64)-A(65)-A(66)-2*A(68)-A(77)-A(87)-A(98)-A(110)+0.001*A(133)+0.042*A(138)+0.025*A(167)&
                &+0.041*A(171)+0.051*A(173)+0.07*A(175)+0.04*A(176)+0.173*A(177)+0.095*A(178)+0.093*A(190)+0.008*A(191)&
                &+0.064*A(192)+0.01*A(193)+0.25*A(194)+0.18*A(195)+0.25*A(196)+0.25*A(198)+0.18*A(199)+0.25*A(200)+0.035&
                &*A(203)+0.07*A(205)+0.143*A(206)+0.347*A(207)+0.011*A(208)+0.009*A(209)+0.09*A(210)+0.001*A(211)+0.176&
                &*A(212)+0.082*A(214)+0.002*A(215)+0.136*A(216)+0.001*A(217)+0.016*A(219)+0.051*A(221)
      
END SUBROUTINE saprc99_mosaic_20bin_vbs2_aq_Fun
















SUBROUTINE saprc99_mosaic_20bin_vbs2_aq_IRRFun ( V, F, RCT, IRR )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: IRR(NREACT)



  IRR(1) = RCT(1)*V(93)
  IRR(2) = RCT(2)*V(85)*F(2)
  IRR(3) = RCT(3)*V(85)*V(89)
  IRR(4) = RCT(4)*V(85)*V(96)*F(2)
  IRR(5) = RCT(5)*V(85)*V(93)
  IRR(6) = RCT(6)*V(85)*V(93)
  IRR(7) = RCT(7)*V(89)*V(96)
  IRR(8) = RCT(8)*V(89)*V(93)
  IRR(9) = RCT(9)*V(94)*V(96)
  IRR(10) = RCT(10)*V(96)*V(96)*F(2)
  IRR(11) = RCT(11)*V(93)*V(94)
  IRR(12) = RCT(12)*V(41)
  IRR(13) = RCT(13)*V(41)*F(1)
  IRR(14) = RCT(14)*V(93)*V(94)
  IRR(15) = RCT(15)*V(94)
  IRR(16) = RCT(16)*V(94)
  IRR(17) = RCT(17)*V(89)
  IRR(18) = RCT(18)*V(89)
  IRR(19) = RCT(19)*V(26)*F(1)
  IRR(20) = RCT(20)*V(26)*F(2)
  IRR(21) = RCT(21)*V(95)*V(96)
  IRR(22) = RCT(22)*V(44)
  IRR(23) = RCT(23)*V(44)
  IRR(24) = RCT(24)*V(44)*V(95)
  IRR(25) = RCT(25)*V(93)*V(95)
  IRR(26) = RCT(26)*V(94)*V(95)
  IRR(27) = RCT(27)*V(67)*V(95)
  IRR(28) = RCT(28)*V(67)
  IRR(29) = RCT(29)*V(66)*V(95)
  IRR(30) = RCT(30)*V(89)*V(95)
  IRR(31) = RCT(31)*V(96)*V(97)
  IRR(32) = RCT(32)*V(93)*V(97)
  IRR(33) = RCT(33)*V(51)
  IRR(34) = RCT(34)*V(51)
  IRR(35) = RCT(35)*V(51)*V(95)
  IRR(36) = RCT(36)*V(89)*V(97)
  IRR(37) = RCT(37)*V(97)*V(97)
  IRR(38) = RCT(38)*V(97)*V(97)*F(1)
  IRR(39) = RCT(39)*V(94)*V(97)
  IRR(40) = RCT(40)*V(94)*V(94)
  IRR(41) = RCT(41)*V(37)
  IRR(42) = RCT(42)*V(37)*V(95)
  IRR(43) = RCT(43)*V(95)*V(97)
  IRR(44) = RCT(44)*V(42)*V(95)
  IRR(45) = RCT(45)*V(95)*F(2)
  IRR(46) = RCT(46)*V(96)*V(98)
  IRR(47) = RCT(47)*V(97)*V(98)
  IRR(48) = RCT(48)*V(94)*V(98)
  IRR(49) = RCT(49)*V(98)*V(98)
  IRR(50) = RCT(50)*V(98)*V(98)
  IRR(51) = RCT(51)*V(96)*V(99)
  IRR(52) = RCT(52)*V(97)*V(99)
  IRR(53) = RCT(53)*V(94)*V(99)
  IRR(54) = RCT(54)*V(98)*V(99)
  IRR(55) = RCT(55)*V(99)*V(99)
  IRR(56) = RCT(56)*V(74)*V(96)
  IRR(57) = RCT(57)*V(74)*V(97)
  IRR(58) = RCT(58)*V(74)*V(94)
  IRR(59) = RCT(59)*V(74)*V(98)
  IRR(60) = RCT(60)*V(74)*V(99)
  IRR(61) = RCT(61)*V(74)*V(74)
  IRR(62) = RCT(62)*V(96)*V(101)
  IRR(63) = RCT(63)*V(97)*V(101)
  IRR(64) = RCT(64)*V(98)*V(101)
  IRR(65) = RCT(65)*V(94)*V(101)
  IRR(66) = RCT(66)*V(99)*V(101)
  IRR(67) = RCT(67)*V(74)*V(101)
  IRR(68) = RCT(68)*V(101)*V(101)
  IRR(69) = RCT(69)*V(90)*V(93)
  IRR(70) = RCT(70)*V(33)
  IRR(71) = RCT(71)*V(90)*V(96)
  IRR(72) = RCT(72)*V(90)*V(97)
  IRR(73) = RCT(73)*V(90)*V(94)
  IRR(74) = RCT(74)*V(90)*V(98)
  IRR(75) = RCT(75)*V(90)*V(99)
  IRR(76) = RCT(76)*V(74)*V(90)
  IRR(77) = RCT(77)*V(90)*V(101)
  IRR(78) = RCT(78)*V(90)*V(90)
  IRR(79) = RCT(79)*V(92)*V(93)
  IRR(80) = RCT(80)*V(34)
  IRR(81) = RCT(81)*V(92)*V(96)
  IRR(82) = RCT(82)*V(92)*V(97)
  IRR(83) = RCT(83)*V(92)*V(94)
  IRR(84) = RCT(84)*V(92)*V(98)
  IRR(85) = RCT(85)*V(92)*V(99)
  IRR(86) = RCT(86)*V(74)*V(92)
  IRR(87) = RCT(87)*V(92)*V(101)
  IRR(88) = RCT(88)*V(90)*V(92)
  IRR(89) = RCT(89)*V(92)*V(92)
  IRR(90) = RCT(90)*V(91)*V(93)
  IRR(91) = RCT(91)*V(35)
  IRR(92) = RCT(92)*V(91)*V(96)
  IRR(93) = RCT(93)*V(91)*V(97)
  IRR(94) = RCT(94)*V(91)*V(94)
  IRR(95) = RCT(95)*V(91)*V(98)
  IRR(96) = RCT(96)*V(91)*V(99)
  IRR(97) = RCT(97)*V(74)*V(91)
  IRR(98) = RCT(98)*V(91)*V(101)
  IRR(99) = RCT(99)*V(90)*V(91)
  IRR(100) = RCT(100)*V(91)*V(92)
  IRR(101) = RCT(101)*V(91)*V(91)
  IRR(102) = RCT(102)*V(93)*V(100)
  IRR(103) = RCT(103)*V(36)
  IRR(104) = RCT(104)*V(96)*V(100)
  IRR(105) = RCT(105)*V(97)*V(100)
  IRR(106) = RCT(106)*V(94)*V(100)
  IRR(107) = RCT(107)*V(98)*V(100)
  IRR(108) = RCT(108)*V(99)*V(100)
  IRR(109) = RCT(109)*V(74)*V(100)
  IRR(110) = RCT(110)*V(100)*V(101)
  IRR(111) = RCT(111)*V(90)*V(100)
  IRR(112) = RCT(112)*V(92)*V(100)
  IRR(113) = RCT(113)*V(91)*V(100)
  IRR(114) = RCT(114)*V(100)*V(100)
  IRR(115) = RCT(115)*V(46)*V(93)
  IRR(116) = RCT(116)*V(46)
  IRR(117) = RCT(117)*V(71)*V(93)
  IRR(118) = RCT(118)*V(71)*V(97)
  IRR(119) = RCT(119)*V(71)
  IRR(120) = RCT(120)*V(50)*V(93)
  IRR(121) = RCT(121)*V(50)*V(97)
  IRR(122) = RCT(122)*V(50)
  IRR(123) = RCT(123)*V(83)
  IRR(124) = RCT(124)*V(83)
  IRR(125) = RCT(125)*V(83)*V(95)
  IRR(126) = RCT(126)*V(83)*V(97)
  IRR(127) = RCT(127)*V(49)
  IRR(128) = RCT(128)*V(49)*V(96)
  IRR(129) = RCT(129)*V(83)*V(94)
  IRR(130) = RCT(130)*V(82)*V(95)
  IRR(131) = RCT(131)*V(82)
  IRR(132) = RCT(132)*V(82)*V(94)
  IRR(133) = RCT(133)*V(86)*V(95)
  IRR(134) = RCT(134)*V(86)
  IRR(135) = RCT(135)*V(86)*V(94)
  IRR(136) = RCT(136)*V(69)*V(95)
  IRR(137) = RCT(137)*V(69)
  IRR(138) = RCT(138)*V(87)*V(95)
  IRR(139) = RCT(139)*V(87)
  IRR(140) = RCT(140)*V(53)*V(95)
  IRR(141) = RCT(141)*V(40)*V(95)
  IRR(142) = RCT(142)*V(48)*V(95)
  IRR(143) = RCT(143)*V(48)
  IRR(144) = RCT(144)*V(62)*V(95)
  IRR(145) = RCT(145)*V(62)
  IRR(146) = RCT(146)*V(78)
  IRR(147) = RCT(147)*V(78)
  IRR(148) = RCT(148)*V(78)*V(95)
  IRR(149) = RCT(149)*V(78)*V(94)
  IRR(150) = RCT(150)*V(65)
  IRR(151) = RCT(151)*V(65)*V(95)
  IRR(152) = RCT(152)*V(65)*V(94)
  IRR(153) = RCT(153)*V(39)
  IRR(154) = RCT(154)*V(64)*V(95)
  IRR(155) = RCT(155)*V(64)*V(94)
  IRR(156) = RCT(156)*V(58)*V(95)
  IRR(157) = RCT(157)*V(58)*V(94)
  IRR(158) = RCT(158)*V(61)*V(94)
  IRR(159) = RCT(159)*V(63)*V(95)
  IRR(160) = RCT(160)*V(63)
  IRR(161) = RCT(161)*V(63)*V(94)
  IRR(162) = RCT(162)*V(75)*V(95)
  IRR(163) = RCT(163)*V(75)*V(89)
  IRR(164) = RCT(164)*V(75)*V(94)
  IRR(165) = RCT(165)*V(75)*V(85)
  IRR(166) = RCT(166)*V(75)
  IRR(167) = RCT(167)*V(81)*V(95)
  IRR(168) = RCT(168)*V(81)*V(89)
  IRR(169) = RCT(169)*V(81)*V(85)
  IRR(170) = RCT(170)*V(81)
  IRR(171) = RCT(171)*V(79)*V(95)
  IRR(172) = RCT(172)*V(79)*V(89)
  IRR(173) = RCT(173)*V(79)*V(94)
  IRR(174) = RCT(174)*V(79)
  IRR(175) = RCT(175)*V(88)*V(95)
  IRR(176) = RCT(176)*V(88)
  IRR(177) = RCT(177)*V(84)*V(95)
  IRR(178) = RCT(178)*V(84)
  IRR(179) = RCT(179)*V(59)*V(95)
  IRR(180) = RCT(180)*V(59)*V(89)
  IRR(181) = RCT(181)*V(56)*V(95)
  IRR(182) = RCT(182)*V(56)
  IRR(183) = RCT(183)*V(57)*V(95)
  IRR(184) = RCT(184)*V(57)
  IRR(185) = RCT(185)*V(31)*V(95)
  IRR(186) = RCT(186)*V(68)*V(95)
  IRR(187) = RCT(187)*V(68)*V(89)
  IRR(188) = RCT(188)*V(68)*V(94)
  IRR(189) = RCT(189)*V(68)*V(85)
  IRR(190) = RCT(190)*V(73)*V(95)
  IRR(191) = RCT(191)*V(73)*V(89)
  IRR(192) = RCT(192)*V(73)*V(94)
  IRR(193) = RCT(193)*V(73)*V(85)
  IRR(194) = RCT(194)*V(77)*V(95)
  IRR(195) = RCT(195)*V(77)*V(89)
  IRR(196) = RCT(196)*V(77)*V(94)
  IRR(197) = RCT(197)*V(77)*V(85)
  IRR(198) = RCT(198)*V(76)*V(95)
  IRR(199) = RCT(199)*V(76)*V(89)
  IRR(200) = RCT(200)*V(76)*V(94)
  IRR(201) = RCT(201)*V(76)*V(85)
  IRR(202) = RCT(202)*V(32)*V(95)
  IRR(203) = RCT(203)*V(38)*V(95)
  IRR(204) = RCT(204)*V(60)*V(95)
  IRR(205) = RCT(205)*V(45)*V(95)
  IRR(206) = RCT(206)*V(54)*V(95)
  IRR(207) = RCT(207)*V(47)*V(95)
  IRR(208) = RCT(208)*V(55)*V(95)
  IRR(209) = RCT(209)*V(52)*V(95)
  IRR(210) = RCT(210)*V(72)*V(95)
  IRR(211) = RCT(211)*V(72)*V(89)
  IRR(212) = RCT(212)*V(72)*V(94)
  IRR(213) = RCT(213)*V(72)*V(85)
  IRR(214) = RCT(214)*V(80)*V(95)
  IRR(215) = RCT(215)*V(80)*V(89)
  IRR(216) = RCT(216)*V(80)*V(94)
  IRR(217) = RCT(217)*V(80)*V(85)
  IRR(218) = RCT(218)*V(60)*V(89)
  IRR(219) = RCT(219)*V(70)*V(95)
  IRR(220) = RCT(220)*V(70)*V(89)
  IRR(221) = RCT(221)*V(70)*V(94)
  IRR(222) = RCT(222)*V(70)*V(85)
  IRR(223) = RCT(223)*V(42)
  IRR(224) = RCT(224)*V(97)
  IRR(225) = RCT(225)*V(42)
  IRR(226) = RCT(226)*V(1)
  IRR(227) = RCT(227)*V(67)
  IRR(228) = RCT(228)*V(37)
  IRR(229) = RCT(229)*V(2)
  IRR(230) = RCT(230)*V(54)*V(95)
  IRR(231) = RCT(231)*V(47)*V(95)
  IRR(232) = RCT(232)*V(72)*V(95)
  IRR(233) = RCT(233)*V(80)*V(95)
  IRR(234) = RCT(234)*V(55)*V(95)
  IRR(235) = RCT(235)*V(52)*V(95)
  IRR(236) = RCT(236)*V(73)*V(95)
  IRR(237) = RCT(237)*V(73)*V(89)
  IRR(238) = RCT(238)*V(73)*V(94)
  IRR(239) = RCT(239)*V(77)*V(95)
  IRR(240) = RCT(240)*V(77)*V(89)
  IRR(241) = RCT(241)*V(77)*V(94)
  IRR(242) = RCT(242)*V(76)*V(95)
  IRR(243) = RCT(243)*V(76)*V(89)
  IRR(244) = RCT(244)*V(76)*V(94)
  IRR(245) = RCT(245)*V(3)*V(95)
  IRR(246) = RCT(246)*V(4)*V(95)
  IRR(247) = RCT(247)*V(28)*V(95)
  IRR(248) = RCT(248)*V(22)*V(95)
  IRR(249) = RCT(249)*V(5)*V(95)
  IRR(250) = RCT(250)*V(6)*V(95)
  IRR(251) = RCT(251)*V(30)*V(95)
  IRR(252) = RCT(252)*V(24)*V(95)
  IRR(253) = RCT(253)*V(27)*V(95)
  IRR(254) = RCT(254)*V(23)*V(95)
  IRR(255) = RCT(255)*V(29)*V(95)
  IRR(256) = RCT(256)*V(25)*V(95)
  IRR(257) = RCT(257)*V(43)*V(95)
  IRR(258) = RCT(258)*V(43)*V(95)
  IRR(259) = RCT(259)*V(43)*V(94)
      
END SUBROUTINE saprc99_mosaic_20bin_vbs2_aq_IRRFun
















SUBROUTINE saprc99_mosaic_20bin_vbs2_aq_Jac_SP ( V, F, RCT, JVS )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: JVS(LU_NONZERO)




  REAL(kind=dp) :: B(461)


  B(1) = RCT(1)

  B(2) = RCT(2)*F(2)

  B(4) = RCT(3)*V(89)

  B(5) = RCT(3)*V(85)

  B(6) = RCT(4)*V(96)*F(2)

  B(7) = RCT(4)*V(85)*F(2)

  B(9) = RCT(5)*V(93)

  B(10) = RCT(5)*V(85)

  B(11) = RCT(6)*V(93)

  B(12) = RCT(6)*V(85)

  B(13) = RCT(7)*V(96)

  B(14) = RCT(7)*V(89)

  B(15) = RCT(8)*V(93)

  B(16) = RCT(8)*V(89)

  B(17) = RCT(9)*V(96)

  B(18) = RCT(9)*V(94)

  B(19) = RCT(10)*2*V(96)*F(2)

  B(21) = RCT(11)*V(94)

  B(22) = RCT(11)*V(93)

  B(23) = RCT(12)

  B(24) = RCT(13)*F(1)

  B(26) = RCT(14)*V(94)

  B(27) = RCT(14)*V(93)

  B(28) = RCT(15)

  B(29) = RCT(16)

  B(30) = RCT(17)

  B(31) = RCT(18)

  B(32) = RCT(19)*F(1)

  B(34) = RCT(20)*F(2)

  B(36) = RCT(21)*V(96)

  B(37) = RCT(21)*V(95)

  B(38) = RCT(22)

  B(39) = RCT(23)

  B(40) = RCT(24)*V(95)

  B(41) = RCT(24)*V(44)

  B(42) = RCT(25)*V(95)

  B(43) = RCT(25)*V(93)

  B(44) = RCT(26)*V(95)

  B(45) = RCT(26)*V(94)

  B(46) = RCT(27)*V(95)

  B(47) = RCT(27)*V(67)

  B(48) = RCT(28)

  B(49) = RCT(29)*V(95)

  B(50) = RCT(29)*V(66)

  B(51) = RCT(30)*V(95)

  B(52) = RCT(30)*V(89)

  B(53) = RCT(31)*V(97)

  B(54) = RCT(31)*V(96)

  B(55) = RCT(32)*V(97)

  B(56) = RCT(32)*V(93)

  B(57) = RCT(33)

  B(58) = RCT(34)

  B(59) = RCT(35)*V(95)

  B(60) = RCT(35)*V(51)

  B(61) = RCT(36)*V(97)

  B(62) = RCT(36)*V(89)

  B(63) = RCT(37)*2*V(97)

  B(64) = RCT(38)*2*V(97)*F(1)

  B(66) = RCT(39)*V(97)

  B(67) = RCT(39)*V(94)

  B(68) = RCT(40)*2*V(94)

  B(69) = RCT(41)

  B(70) = RCT(42)*V(95)

  B(71) = RCT(42)*V(37)

  B(72) = RCT(43)*V(97)

  B(73) = RCT(43)*V(95)

  B(74) = RCT(44)*V(95)

  B(75) = RCT(44)*V(42)

  B(76) = RCT(45)*F(2)

  B(78) = RCT(46)*V(98)

  B(79) = RCT(46)*V(96)

  B(80) = RCT(47)*V(98)

  B(81) = RCT(47)*V(97)

  B(82) = RCT(48)*V(98)

  B(83) = RCT(48)*V(94)

  B(84) = RCT(49)*2*V(98)

  B(85) = RCT(50)*2*V(98)

  B(86) = RCT(51)*V(99)

  B(87) = RCT(51)*V(96)

  B(88) = RCT(52)*V(99)

  B(89) = RCT(52)*V(97)

  B(90) = RCT(53)*V(99)

  B(91) = RCT(53)*V(94)

  B(92) = RCT(54)*V(99)

  B(93) = RCT(54)*V(98)

  B(94) = RCT(55)*2*V(99)

  B(95) = RCT(56)*V(96)

  B(96) = RCT(56)*V(74)

  B(97) = RCT(57)*V(97)

  B(98) = RCT(57)*V(74)

  B(99) = RCT(58)*V(94)

  B(100) = RCT(58)*V(74)

  B(101) = RCT(59)*V(98)

  B(102) = RCT(59)*V(74)

  B(103) = RCT(60)*V(99)

  B(104) = RCT(60)*V(74)

  B(105) = RCT(61)*2*V(74)

  B(106) = RCT(62)*V(101)

  B(107) = RCT(62)*V(96)

  B(108) = RCT(63)*V(101)

  B(109) = RCT(63)*V(97)

  B(110) = RCT(64)*V(101)

  B(111) = RCT(64)*V(98)

  B(112) = RCT(65)*V(101)

  B(113) = RCT(65)*V(94)

  B(114) = RCT(66)*V(101)

  B(115) = RCT(66)*V(99)

  B(116) = RCT(67)*V(101)

  B(117) = RCT(67)*V(74)

  B(118) = RCT(68)*2*V(101)

  B(119) = RCT(69)*V(93)

  B(120) = RCT(69)*V(90)

  B(121) = RCT(70)

  B(122) = RCT(71)*V(96)

  B(123) = RCT(71)*V(90)

  B(124) = RCT(72)*V(97)

  B(125) = RCT(72)*V(90)

  B(126) = RCT(73)*V(94)

  B(127) = RCT(73)*V(90)

  B(128) = RCT(74)*V(98)

  B(129) = RCT(74)*V(90)

  B(130) = RCT(75)*V(99)

  B(131) = RCT(75)*V(90)

  B(132) = RCT(76)*V(90)

  B(133) = RCT(76)*V(74)

  B(134) = RCT(77)*V(101)

  B(135) = RCT(77)*V(90)

  B(136) = RCT(78)*2*V(90)

  B(137) = RCT(79)*V(93)

  B(138) = RCT(79)*V(92)

  B(139) = RCT(80)

  B(140) = RCT(81)*V(96)

  B(141) = RCT(81)*V(92)

  B(142) = RCT(82)*V(97)

  B(143) = RCT(82)*V(92)

  B(144) = RCT(83)*V(94)

  B(145) = RCT(83)*V(92)

  B(146) = RCT(84)*V(98)

  B(147) = RCT(84)*V(92)

  B(148) = RCT(85)*V(99)

  B(149) = RCT(85)*V(92)

  B(150) = RCT(86)*V(92)

  B(151) = RCT(86)*V(74)

  B(152) = RCT(87)*V(101)

  B(153) = RCT(87)*V(92)

  B(154) = RCT(88)*V(92)

  B(155) = RCT(88)*V(90)

  B(156) = RCT(89)*2*V(92)

  B(157) = RCT(90)*V(93)

  B(158) = RCT(90)*V(91)

  B(159) = RCT(91)

  B(160) = RCT(92)*V(96)

  B(161) = RCT(92)*V(91)

  B(162) = RCT(93)*V(97)

  B(163) = RCT(93)*V(91)

  B(164) = RCT(94)*V(94)

  B(165) = RCT(94)*V(91)

  B(166) = RCT(95)*V(98)

  B(167) = RCT(95)*V(91)

  B(168) = RCT(96)*V(99)

  B(169) = RCT(96)*V(91)

  B(170) = RCT(97)*V(91)

  B(171) = RCT(97)*V(74)

  B(172) = RCT(98)*V(101)

  B(173) = RCT(98)*V(91)

  B(174) = RCT(99)*V(91)

  B(175) = RCT(99)*V(90)

  B(176) = RCT(100)*V(92)

  B(177) = RCT(100)*V(91)

  B(178) = RCT(101)*2*V(91)

  B(179) = RCT(102)*V(100)

  B(180) = RCT(102)*V(93)

  B(181) = RCT(103)

  B(182) = RCT(104)*V(100)

  B(183) = RCT(104)*V(96)

  B(184) = RCT(105)*V(100)

  B(185) = RCT(105)*V(97)

  B(186) = RCT(106)*V(100)

  B(187) = RCT(106)*V(94)

  B(188) = RCT(107)*V(100)

  B(189) = RCT(107)*V(98)

  B(190) = RCT(108)*V(100)

  B(191) = RCT(108)*V(99)

  B(192) = RCT(109)*V(100)

  B(193) = RCT(109)*V(74)

  B(194) = RCT(110)*V(101)

  B(195) = RCT(110)*V(100)

  B(196) = RCT(111)*V(100)

  B(197) = RCT(111)*V(90)

  B(198) = RCT(112)*V(100)

  B(199) = RCT(112)*V(92)

  B(200) = RCT(113)*V(100)

  B(201) = RCT(113)*V(91)

  B(202) = RCT(114)*2*V(100)

  B(203) = RCT(115)*V(93)

  B(204) = RCT(115)*V(46)

  B(205) = RCT(116)

  B(206) = RCT(117)*V(93)

  B(207) = RCT(117)*V(71)

  B(208) = RCT(118)*V(97)

  B(209) = RCT(118)*V(71)

  B(210) = RCT(119)

  B(211) = RCT(120)*V(93)

  B(212) = RCT(120)*V(50)

  B(213) = RCT(121)*V(97)

  B(214) = RCT(121)*V(50)

  B(215) = RCT(122)

  B(216) = RCT(123)

  B(217) = RCT(124)

  B(218) = RCT(125)*V(95)

  B(219) = RCT(125)*V(83)

  B(220) = RCT(126)*V(97)

  B(221) = RCT(126)*V(83)

  B(222) = RCT(127)

  B(223) = RCT(128)*V(96)

  B(224) = RCT(128)*V(49)

  B(225) = RCT(129)*V(94)

  B(226) = RCT(129)*V(83)

  B(227) = RCT(130)*V(95)

  B(228) = RCT(130)*V(82)

  B(229) = RCT(131)

  B(230) = RCT(132)*V(94)

  B(231) = RCT(132)*V(82)

  B(232) = RCT(133)*V(95)

  B(233) = RCT(133)*V(86)

  B(234) = RCT(134)

  B(235) = RCT(135)*V(94)

  B(236) = RCT(135)*V(86)

  B(237) = RCT(136)*V(95)

  B(238) = RCT(136)*V(69)

  B(239) = RCT(137)

  B(240) = RCT(138)*V(95)

  B(241) = RCT(138)*V(87)

  B(242) = RCT(139)

  B(243) = RCT(140)*V(95)

  B(244) = RCT(140)*V(53)

  B(245) = RCT(141)*V(95)

  B(246) = RCT(141)*V(40)

  B(247) = RCT(142)*V(95)

  B(248) = RCT(142)*V(48)

  B(249) = RCT(143)

  B(250) = RCT(144)*V(95)

  B(251) = RCT(144)*V(62)

  B(252) = RCT(145)

  B(253) = RCT(146)

  B(254) = RCT(147)

  B(255) = RCT(148)*V(95)

  B(256) = RCT(148)*V(78)

  B(257) = RCT(149)*V(94)

  B(258) = RCT(149)*V(78)

  B(259) = RCT(150)

  B(260) = RCT(151)*V(95)

  B(261) = RCT(151)*V(65)

  B(262) = RCT(152)*V(94)

  B(263) = RCT(152)*V(65)

  B(264) = RCT(153)

  B(265) = RCT(154)*V(95)

  B(266) = RCT(154)*V(64)

  B(267) = RCT(155)*V(94)

  B(268) = RCT(155)*V(64)

  B(269) = RCT(156)*V(95)

  B(270) = RCT(156)*V(58)

  B(271) = RCT(157)*V(94)

  B(272) = RCT(157)*V(58)

  B(273) = RCT(158)*V(94)

  B(274) = RCT(158)*V(61)

  B(275) = RCT(159)*V(95)

  B(276) = RCT(159)*V(63)

  B(277) = RCT(160)

  B(278) = RCT(161)*V(94)

  B(279) = RCT(161)*V(63)

  B(280) = RCT(162)*V(95)

  B(281) = RCT(162)*V(75)

  B(282) = RCT(163)*V(89)

  B(283) = RCT(163)*V(75)

  B(284) = RCT(164)*V(94)

  B(285) = RCT(164)*V(75)

  B(286) = RCT(165)*V(85)

  B(287) = RCT(165)*V(75)

  B(288) = RCT(166)

  B(289) = RCT(167)*V(95)

  B(290) = RCT(167)*V(81)

  B(291) = RCT(168)*V(89)

  B(292) = RCT(168)*V(81)

  B(293) = RCT(169)*V(85)

  B(294) = RCT(169)*V(81)

  B(295) = RCT(170)

  B(296) = RCT(171)*V(95)

  B(297) = RCT(171)*V(79)

  B(298) = RCT(172)*V(89)

  B(299) = RCT(172)*V(79)

  B(300) = RCT(173)*V(94)

  B(301) = RCT(173)*V(79)

  B(302) = RCT(174)

  B(303) = RCT(175)*V(95)

  B(304) = RCT(175)*V(88)

  B(305) = RCT(176)

  B(306) = RCT(177)*V(95)

  B(307) = RCT(177)*V(84)

  B(308) = RCT(178)

  B(309) = RCT(179)*V(95)

  B(310) = RCT(179)*V(59)

  B(311) = RCT(180)*V(89)

  B(312) = RCT(180)*V(59)

  B(313) = RCT(181)*V(95)

  B(314) = RCT(181)*V(56)

  B(315) = RCT(182)

  B(316) = RCT(183)*V(95)

  B(317) = RCT(183)*V(57)

  B(318) = RCT(184)

  B(319) = RCT(185)*V(95)

  B(320) = RCT(185)*V(31)

  B(321) = RCT(186)*V(95)

  B(322) = RCT(186)*V(68)

  B(323) = RCT(187)*V(89)

  B(324) = RCT(187)*V(68)

  B(325) = RCT(188)*V(94)

  B(326) = RCT(188)*V(68)

  B(327) = RCT(189)*V(85)

  B(328) = RCT(189)*V(68)

  B(329) = RCT(190)*V(95)

  B(330) = RCT(190)*V(73)

  B(331) = RCT(191)*V(89)

  B(332) = RCT(191)*V(73)

  B(333) = RCT(192)*V(94)

  B(334) = RCT(192)*V(73)

  B(335) = RCT(193)*V(85)

  B(336) = RCT(193)*V(73)

  B(337) = RCT(194)*V(95)

  B(338) = RCT(194)*V(77)

  B(339) = RCT(195)*V(89)

  B(340) = RCT(195)*V(77)

  B(341) = RCT(196)*V(94)

  B(342) = RCT(196)*V(77)

  B(343) = RCT(197)*V(85)

  B(344) = RCT(197)*V(77)

  B(345) = RCT(198)*V(95)

  B(346) = RCT(198)*V(76)

  B(347) = RCT(199)*V(89)

  B(348) = RCT(199)*V(76)

  B(349) = RCT(200)*V(94)

  B(350) = RCT(200)*V(76)

  B(351) = RCT(201)*V(85)

  B(352) = RCT(201)*V(76)

  B(353) = RCT(202)*V(95)

  B(354) = RCT(202)*V(32)

  B(355) = RCT(203)*V(95)

  B(356) = RCT(203)*V(38)

  B(357) = RCT(204)*V(95)

  B(358) = RCT(204)*V(60)

  B(359) = RCT(205)*V(95)

  B(360) = RCT(205)*V(45)

  B(361) = RCT(206)*V(95)

  B(362) = RCT(206)*V(54)

  B(363) = RCT(207)*V(95)

  B(364) = RCT(207)*V(47)

  B(365) = RCT(208)*V(95)

  B(366) = RCT(208)*V(55)

  B(367) = RCT(209)*V(95)

  B(368) = RCT(209)*V(52)

  B(369) = RCT(210)*V(95)

  B(370) = RCT(210)*V(72)

  B(371) = RCT(211)*V(89)

  B(372) = RCT(211)*V(72)

  B(373) = RCT(212)*V(94)

  B(374) = RCT(212)*V(72)

  B(375) = RCT(213)*V(85)

  B(376) = RCT(213)*V(72)

  B(377) = RCT(214)*V(95)

  B(378) = RCT(214)*V(80)

  B(379) = RCT(215)*V(89)

  B(380) = RCT(215)*V(80)

  B(381) = RCT(216)*V(94)

  B(382) = RCT(216)*V(80)

  B(383) = RCT(217)*V(85)

  B(384) = RCT(217)*V(80)

  B(385) = RCT(218)*V(89)

  B(386) = RCT(218)*V(60)

  B(387) = RCT(219)*V(95)

  B(388) = RCT(219)*V(70)

  B(389) = RCT(220)*V(89)

  B(390) = RCT(220)*V(70)

  B(391) = RCT(221)*V(94)

  B(392) = RCT(221)*V(70)

  B(393) = RCT(222)*V(85)

  B(394) = RCT(222)*V(70)

  B(395) = RCT(223)

  B(396) = RCT(224)

  B(397) = RCT(225)

  B(398) = RCT(226)

  B(399) = RCT(227)

  B(400) = RCT(228)

  B(401) = RCT(229)

  B(402) = RCT(230)*V(95)

  B(403) = RCT(230)*V(54)

  B(404) = RCT(231)*V(95)

  B(405) = RCT(231)*V(47)

  B(406) = RCT(232)*V(95)

  B(407) = RCT(232)*V(72)

  B(408) = RCT(233)*V(95)

  B(409) = RCT(233)*V(80)

  B(410) = RCT(234)*V(95)

  B(411) = RCT(234)*V(55)

  B(412) = RCT(235)*V(95)

  B(413) = RCT(235)*V(52)

  B(414) = RCT(236)*V(95)

  B(415) = RCT(236)*V(73)

  B(416) = RCT(237)*V(89)

  B(417) = RCT(237)*V(73)

  B(418) = RCT(238)*V(94)

  B(419) = RCT(238)*V(73)

  B(420) = RCT(239)*V(95)

  B(421) = RCT(239)*V(77)

  B(422) = RCT(240)*V(89)

  B(423) = RCT(240)*V(77)

  B(424) = RCT(241)*V(94)

  B(425) = RCT(241)*V(77)

  B(426) = RCT(242)*V(95)

  B(427) = RCT(242)*V(76)

  B(428) = RCT(243)*V(89)

  B(429) = RCT(243)*V(76)

  B(430) = RCT(244)*V(94)

  B(431) = RCT(244)*V(76)

  B(432) = RCT(245)*V(95)

  B(433) = RCT(245)*V(3)

  B(434) = RCT(246)*V(95)

  B(435) = RCT(246)*V(4)

  B(436) = RCT(247)*V(95)

  B(437) = RCT(247)*V(28)

  B(438) = RCT(248)*V(95)

  B(439) = RCT(248)*V(22)

  B(440) = RCT(249)*V(95)

  B(441) = RCT(249)*V(5)

  B(442) = RCT(250)*V(95)

  B(443) = RCT(250)*V(6)

  B(444) = RCT(251)*V(95)

  B(445) = RCT(251)*V(30)

  B(446) = RCT(252)*V(95)

  B(447) = RCT(252)*V(24)

  B(448) = RCT(253)*V(95)

  B(449) = RCT(253)*V(27)

  B(450) = RCT(254)*V(95)

  B(451) = RCT(254)*V(23)

  B(452) = RCT(255)*V(95)

  B(453) = RCT(255)*V(29)

  B(454) = RCT(256)*V(95)

  B(455) = RCT(256)*V(25)

  B(456) = RCT(257)*V(95)

  B(457) = RCT(257)*V(43)

  B(458) = RCT(258)*V(95)

  B(459) = RCT(258)*V(43)

  B(460) = RCT(259)*V(94)

  B(461) = RCT(259)*V(43)



  JVS(1) = -B(398)

  JVS(2) = B(74)+B(395)

  JVS(3) = B(75)

  JVS(4) = -B(401)

  JVS(5) = 0.5*B(385)

  JVS(6) = 0.135*B(389)

  JVS(7) = 0.5*B(386)+0.135*B(390)

  JVS(8) = 0

  JVS(9) = 0

  JVS(10) = 0

  JVS(11) = 0

  JVS(12) = 0

  JVS(13) = B(223)

  JVS(14) = 0.297*B(357)

  JVS(15) = 0.37*B(323)

  JVS(16) = 0.185*B(389)

  JVS(17) = 0.185*B(371)

  JVS(18) = 0.204*B(331)

  JVS(19) = 0.333*B(282)

  JVS(20) = 0.103*B(347)

  JVS(21) = 0.103*B(339)

  JVS(22) = 0.1*B(298)

  JVS(23) = 0.073*B(379)

  JVS(24) = 0.351*B(291)

  JVS(25) = 0.333*B(283)+0.351*B(292)+0.1*B(299)+0.37*B(324)+0.204*B(332)+0.103*B(340)+0.103*B(348)+0.185*B(372)+0.073&
              &*B(380)+0.185*B(390)

  JVS(26) = 0.297*B(358)

  JVS(27) = B(224)

  JVS(28) = 0

  JVS(29) = 0.17*B(389)

  JVS(30) = 0.05*B(371)

  JVS(31) = 0.129*B(379)

  JVS(32) = 0.05*B(372)+0.129*B(380)+0.17*B(390)

  JVS(33) = 0.25*B(124)+B(128)+B(130)+B(134)

  JVS(34) = 0.25*B(125)

  JVS(35) = B(129)

  JVS(36) = B(131)

  JVS(37) = B(135)

  JVS(38) = 0

  JVS(39) = 0.119*B(371)

  JVS(40) = 0.15*B(331)

  JVS(41) = 0.189*B(347)

  JVS(42) = 0.189*B(339)

  JVS(43) = 0.372*B(298)

  JVS(44) = 0.247*B(379)

  JVS(45) = 0.372*B(299)+0.15*B(332)+0.189*B(340)+0.189*B(348)+0.119*B(372)+0.247*B(380)

  JVS(46) = 0.25*B(162)+B(166)+B(168)+B(172)

  JVS(47) = 0.25*B(142)+B(146)+B(148)+B(152)

  JVS(48) = 0.25*B(143)+0.25*B(163)+0.25*B(184)

  JVS(49) = B(147)+B(167)+B(188)

  JVS(50) = B(149)+B(169)+B(190)

  JVS(51) = 0.25*B(185)+B(189)+B(191)+2*B(194)

  JVS(52) = B(153)+B(173)+2*B(195)

  JVS(53) = 0

  JVS(54) = 0.75*B(124)

  JVS(55) = 0.75*B(125)

  JVS(56) = 0

  JVS(57) = 0.75*B(162)

  JVS(58) = 0.75*B(142)

  JVS(59) = 0.75*B(143)+0.75*B(163)+0.75*B(184)

  JVS(60) = 0.75*B(185)

  JVS(61) = 0

  JVS(62) = 2*B(211)

  JVS(63) = B(391)

  JVS(64) = 2*B(212)

  JVS(65) = B(392)

  JVS(66) = 0

  JVS(67) = 6*B(211)

  JVS(68) = 7*B(277)

  JVS(69) = 0.048*B(387)+0.07*B(389)+2.693*B(391)+0.55*B(393)

  JVS(70) = 0.55*B(394)

  JVS(71) = 0.07*B(390)

  JVS(72) = 6*B(212)

  JVS(73) = 2.693*B(392)

  JVS(74) = 0.048*B(388)

  JVS(75) = 0

  JVS(76) = B(95)+B(99)

  JVS(77) = B(82)+B(90)+B(100)+B(112)

  JVS(78) = B(78)+B(86)+B(96)+B(106)

  JVS(79) = B(79)+B(83)

  JVS(80) = B(87)+B(91)

  JVS(81) = B(107)+B(113)

  JVS(82) = 0

  JVS(83) = B(97)+B(101)+B(103)+B(105)+B(116)

  JVS(84) = B(80)+B(88)+B(98)+B(108)

  JVS(85) = B(81)+B(84)+B(85)+B(92)+B(102)+B(110)

  JVS(86) = B(89)+B(93)+B(94)+B(104)+B(114)

  JVS(87) = B(109)+B(111)+B(115)+B(117)+B(118)

  JVS(88) = 0

  JVS(89) = B(404)

  JVS(90) = B(412)

  JVS(91) = B(402)

  JVS(92) = B(410)

  JVS(93) = B(406)

  JVS(94) = B(408)

  JVS(95) = B(403)+B(405)+B(407)+B(409)+B(411)+B(413)

  JVS(96) = 0

  JVS(97) = B(420)

  JVS(98) = B(421)

  JVS(99) = 0

  JVS(100) = B(422)

  JVS(101) = B(423)

  JVS(102) = 0

  JVS(103) = B(414)+B(416)+0.0327*B(418)

  JVS(104) = B(426)+B(428)+0.0816*B(430)

  JVS(105) = 0.0545*B(424)

  JVS(106) = B(417)+B(429)

  JVS(107) = 0.0327*B(419)+0.0545*B(425)+0.0816*B(431)

  JVS(108) = B(415)+B(427)

  JVS(109) = 0

  JVS(110) = B(319)

  JVS(111) = B(353)

  JVS(112) = B(355)

  JVS(113) = B(245)

  JVS(114) = B(74)

  JVS(115) = B(40)

  JVS(116) = B(359)

  JVS(117) = B(363)

  JVS(118) = B(247)

  JVS(119) = B(367)

  JVS(120) = B(243)

  JVS(121) = B(361)

  JVS(122) = B(365)

  JVS(123) = B(313)

  JVS(124) = B(316)

  JVS(125) = B(269)

  JVS(126) = B(309)

  JVS(127) = B(357)

  JVS(128) = B(250)

  JVS(129) = B(275)

  JVS(130) = B(265)

  JVS(131) = B(260)

  JVS(132) = B(49)

  JVS(133) = B(46)

  JVS(134) = B(321)

  JVS(135) = B(237)

  JVS(136) = B(387)

  JVS(137) = B(369)

  JVS(138) = B(329)

  JVS(139) = B(280)

  JVS(140) = B(345)

  JVS(141) = B(337)

  JVS(142) = B(255)

  JVS(143) = B(296)

  JVS(144) = B(377)

  JVS(145) = B(289)

  JVS(146) = B(227)

  JVS(147) = B(218)

  JVS(148) = B(306)

  JVS(149) = B(232)

  JVS(150) = B(240)

  JVS(151) = B(303)

  JVS(152) = B(51)

  JVS(153) = B(42)

  JVS(154) = B(44)

  JVS(155) = B(36)+B(41)+B(43)+B(45)+B(47)+B(50)+B(52)+B(72)+B(75)+B(76)+B(219)+B(228)+B(233)+B(238)+B(241)+B(244)&
               &+B(246)+B(248)+B(251)+B(256)+B(261)+B(266)+B(270)+B(276)+B(281)+B(290)+B(297)+B(304)+B(307)+B(310)+B(314)&
               &+B(317)+B(320)+B(322)+B(330)+B(338)+B(346)+B(354)+B(356)+B(358)+B(360)+B(362)+B(364)+B(366)+B(368)+B(370)&
               &+B(378)+B(388)

  JVS(156) = B(37)

  JVS(157) = B(73)

  JVS(158) = B(432)

  JVS(159) = B(434)

  JVS(160) = B(440)

  JVS(161) = B(442)

  JVS(162) = 0

  JVS(163) = B(438)

  JVS(164) = B(450)

  JVS(165) = B(446)

  JVS(166) = B(454)

  JVS(167) = B(448)

  JVS(168) = B(436)

  JVS(169) = B(452)

  JVS(170) = B(444)

  JVS(171) = B(433)+B(435)+B(437)+B(439)+B(441)+B(443)+B(445)+B(447)+B(449)+B(451)+B(453)+B(455)

  JVS(172) = 0

  JVS(173) = B(450)

  JVS(174) = 0.5*B(448)

  JVS(175) = 0.5*B(449)+B(451)

  JVS(176) = -B(450)

  JVS(177) = -B(451)

  JVS(178) = 0

  JVS(179) = B(454)

  JVS(180) = 0.5*B(452)

  JVS(181) = 0.5*B(453)+B(455)

  JVS(182) = -B(454)

  JVS(183) = -B(455)

  JVS(184) = -B(32)-B(34)

  JVS(185) = B(31)

  JVS(186) = -B(448)

  JVS(187) = -B(449)

  JVS(188) = B(448)

  JVS(189) = 0

  JVS(190) = B(449)

  JVS(191) = -B(452)

  JVS(192) = -B(453)

  JVS(193) = B(452)

  JVS(194) = 0

  JVS(195) = B(453)

  JVS(196) = -B(319)

  JVS(197) = -B(320)

  JVS(198) = -B(353)

  JVS(199) = -B(354)

  JVS(200) = -B(121)

  JVS(201) = B(119)

  JVS(202) = B(120)

  JVS(203) = -B(139)

  JVS(204) = B(137)

  JVS(205) = B(138)

  JVS(206) = -B(159)

  JVS(207) = B(157)

  JVS(208) = B(158)

  JVS(209) = -B(181)

  JVS(210) = B(179)

  JVS(211) = B(180)

  JVS(212) = -B(69)-B(70)-B(400)

  JVS(213) = -B(71)

  JVS(214) = B(63)+B(64)

  JVS(215) = -B(355)

  JVS(216) = -B(356)

  JVS(217) = -B(264)

  JVS(218) = 0.087*B(367)

  JVS(219) = 0.031*B(347)

  JVS(220) = 0.031*B(339)

  JVS(221) = 0.031*B(340)+0.031*B(348)

  JVS(222) = 0.087*B(368)

  JVS(223) = -B(245)

  JVS(224) = -B(246)

  JVS(225) = -B(23)-B(24)

  JVS(226) = B(21)

  JVS(227) = B(22)

  JVS(228) = -B(74)-B(395)-B(397)

  JVS(229) = B(456)+0.5*B(458)+B(460)

  JVS(230) = B(461)

  JVS(231) = -B(75)+B(457)+0.5*B(459)

  JVS(232) = -B(456)-B(458)-B(460)

  JVS(233) = -B(461)

  JVS(234) = -B(457)-B(459)

  JVS(235) = -B(38)-B(39)-B(40)

  JVS(236) = B(36)-B(41)

  JVS(237) = B(37)

  JVS(238) = -B(359)

  JVS(239) = -B(360)

  JVS(240) = 0.236*B(359)

  JVS(241) = -B(203)-B(205)

  JVS(242) = -B(204)

  JVS(243) = 0.236*B(360)

  JVS(244) = -B(363)

  JVS(245) = -B(364)

  JVS(246) = -B(247)-B(249)

  JVS(247) = -B(248)

  JVS(248) = B(80)

  JVS(249) = B(81)

  JVS(250) = -B(222)-B(223)

  JVS(251) = B(220)

  JVS(252) = -B(224)

  JVS(253) = B(221)

  JVS(254) = -B(211)-B(213)-B(215)

  JVS(255) = B(273)

  JVS(256) = -B(212)

  JVS(257) = B(274)

  JVS(258) = -B(214)

  JVS(259) = -B(57)-B(58)-B(59)

  JVS(260) = B(55)

  JVS(261) = -B(60)

  JVS(262) = B(56)

  JVS(263) = -B(367)

  JVS(264) = -B(368)

  JVS(265) = -B(243)

  JVS(266) = -B(244)

  JVS(267) = B(84)+0.25*B(92)+0.25*B(110)

  JVS(268) = 0.25*B(93)

  JVS(269) = 0.25*B(111)

  JVS(270) = -B(361)

  JVS(271) = -B(362)

  JVS(272) = -B(365)

  JVS(273) = -B(366)

  JVS(274) = 0.099*B(367)

  JVS(275) = 0.108*B(365)

  JVS(276) = -B(313)-B(315)

  JVS(277) = -B(314)+0.108*B(366)+0.099*B(368)

  JVS(278) = 0.093*B(367)

  JVS(279) = 0.051*B(365)

  JVS(280) = -B(316)-B(318)

  JVS(281) = -B(317)+0.051*B(366)+0.093*B(368)

  JVS(282) = 0.187*B(367)

  JVS(283) = 0.207*B(365)

  JVS(284) = -B(269)-B(271)

  JVS(285) = -B(272)

  JVS(286) = -B(270)+0.207*B(366)+0.187*B(368)

  JVS(287) = 0.561*B(367)

  JVS(288) = 0.491*B(365)

  JVS(289) = -B(309)-B(311)

  JVS(290) = -B(312)

  JVS(291) = -B(310)+0.491*B(366)+0.561*B(368)

  JVS(292) = -B(357)-B(385)

  JVS(293) = -B(386)

  JVS(294) = -B(358)

  JVS(295) = B(213)+B(215)

  JVS(296) = -B(273)

  JVS(297) = B(206)

  JVS(298) = B(207)

  JVS(299) = -B(274)

  JVS(300) = B(214)

  JVS(301) = -B(250)-B(252)

  JVS(302) = -B(251)

  JVS(303) = B(88)+B(108)

  JVS(304) = B(89)

  JVS(305) = B(109)

  JVS(306) = 0.05*B(367)

  JVS(307) = 0.059*B(365)

  JVS(308) = -B(275)-B(277)-B(278)

  JVS(309) = 0.061*B(377)+0.042*B(379)+0.015*B(381)

  JVS(310) = 0.042*B(380)

  JVS(311) = -B(279)+0.015*B(382)

  JVS(312) = -B(276)+0.059*B(366)+0.05*B(368)+0.061*B(378)

  JVS(313) = 0.017*B(365)

  JVS(314) = -B(265)-B(267)

  JVS(315) = B(208)+B(210)

  JVS(316) = -B(268)

  JVS(317) = -B(266)+0.017*B(366)

  JVS(318) = B(209)

  JVS(319) = 0.287*B(367)

  JVS(320) = 0.119*B(365)

  JVS(321) = 0.5*B(315)

  JVS(322) = 0.5*B(318)

  JVS(323) = 0.23*B(269)

  JVS(324) = -B(259)-B(260)-B(262)

  JVS(325) = 0.084*B(280)+0.9*B(282)

  JVS(326) = 0.174*B(296)+0.742*B(298)+0.008*B(300)

  JVS(327) = 0.3*B(289)+0.95*B(291)

  JVS(328) = 0.9*B(283)+0.95*B(292)+0.742*B(299)

  JVS(329) = -B(263)+0.008*B(301)

  JVS(330) = -B(261)+0.23*B(270)+0.084*B(281)+0.3*B(290)+0.174*B(297)+0.119*B(366)+0.287*B(368)

  JVS(331) = 0.002*B(361)

  JVS(332) = B(315)

  JVS(333) = B(318)

  JVS(334) = B(309)+1.5*B(311)

  JVS(335) = 0.393*B(357)+1.5*B(385)

  JVS(336) = B(259)+B(260)+B(262)

  JVS(337) = -B(49)

  JVS(338) = 0.5*B(323)+0.491*B(327)

  JVS(339) = 0.51*B(389)

  JVS(340) = 0.345*B(371)

  JVS(341) = 0.275*B(331)

  JVS(342) = 0.416*B(280)+0.45*B(282)+0.5*B(284)+0.67*B(288)

  JVS(343) = 0.157*B(347)

  JVS(344) = 0.157*B(339)

  JVS(345) = 2*B(253)+B(254)+1.26*B(255)+1.26*B(257)

  JVS(346) = 0.336*B(296)+0.498*B(298)+0.572*B(300)+1.233*B(302)

  JVS(347) = 0.265*B(379)+0.012*B(383)

  JVS(348) = 0.475*B(291)+0.7*B(295)

  JVS(349) = B(229)

  JVS(350) = B(216)+B(217)+B(218)+B(225)

  JVS(351) = 0.491*B(328)+0.012*B(384)

  JVS(352) = 0.034*B(232)+B(234)

  JVS(353) = 0.45*B(283)+0.475*B(292)+0.498*B(299)+1.5*B(312)+0.5*B(324)+0.275*B(332)+0.157*B(340)+0.157*B(348)+0.345&
               &*B(372)+0.265*B(380)+1.5*B(386)+0.51*B(390)

  JVS(354) = B(226)+1.26*B(258)+B(263)+0.5*B(285)+0.572*B(301)

  JVS(355) = -B(50)+B(219)+0.034*B(233)+1.26*B(256)+B(261)+0.416*B(281)+0.336*B(297)+B(310)+0.393*B(358)+0.002*B(362)

  JVS(356) = 2*B(24)

  JVS(357) = B(460)

  JVS(358) = B(271)

  JVS(359) = B(273)

  JVS(360) = B(278)

  JVS(361) = B(267)

  JVS(362) = B(262)

  JVS(363) = -B(46)-B(48)-B(399)

  JVS(364) = 0

  JVS(365) = 0.5*B(284)

  JVS(366) = B(257)

  JVS(367) = 0.15*B(300)

  JVS(368) = 0

  JVS(369) = 0

  JVS(370) = B(230)

  JVS(371) = B(225)

  JVS(372) = B(235)

  JVS(373) = 0

  JVS(374) = B(42)

  JVS(375) = 0.2*B(66)+B(226)+B(231)+B(236)+B(258)+B(263)+B(268)+B(272)+B(274)+B(279)+0.5*B(285)+0.15*B(301)+B(461)

  JVS(376) = B(43)-B(47)

  JVS(377) = 0.2*B(67)

  JVS(378) = -B(321)-B(323)-B(325)-B(327)

  JVS(379) = -B(328)

  JVS(380) = -B(324)

  JVS(381) = -B(326)

  JVS(382) = -B(322)

  JVS(383) = 0.704*B(355)

  JVS(384) = 0.024*B(359)

  JVS(385) = B(205)

  JVS(386) = 0.072*B(363)

  JVS(387) = 0.452*B(361)

  JVS(388) = -B(237)-B(239)

  JVS(389) = 0.005*B(369)+0.001*B(371)+0.024*B(373)

  JVS(390) = 0.13*B(347)

  JVS(391) = 0.13*B(339)

  JVS(392) = 0.127*B(377)+0.045*B(379)+0.102*B(381)

  JVS(393) = 0.006*B(306)+0.02*B(308)

  JVS(394) = 0.13*B(340)+0.13*B(348)+0.001*B(372)+0.045*B(380)

  JVS(395) = 0

  JVS(396) = 0.024*B(374)+0.102*B(382)

  JVS(397) = -B(238)+0.006*B(307)+0.704*B(356)+0.024*B(360)+0.452*B(362)+0.072*B(364)+0.005*B(370)+0.127*B(378)

  JVS(398) = -B(387)-B(389)-B(391)-B(393)

  JVS(399) = -B(394)

  JVS(400) = -B(390)

  JVS(401) = -B(392)

  JVS(402) = -B(388)

  JVS(403) = 0.24*B(269)+B(271)

  JVS(404) = 0.24*B(265)+B(267)

  JVS(405) = -B(206)-B(208)-B(210)

  JVS(406) = B(174)

  JVS(407) = B(160)+B(164)+B(175)+B(176)+2*B(178)+B(200)

  JVS(408) = B(177)

  JVS(409) = -B(207)

  JVS(410) = B(165)+B(268)+B(272)

  JVS(411) = 0.24*B(266)+0.24*B(270)

  JVS(412) = B(161)

  JVS(413) = -B(209)

  JVS(414) = B(201)

  JVS(415) = -B(369)-B(371)-B(373)-B(375)

  JVS(416) = -B(376)

  JVS(417) = -B(372)

  JVS(418) = -B(374)

  JVS(419) = -B(370)

  JVS(420) = -B(329)-B(331)-B(333)-B(335)

  JVS(421) = -B(336)

  JVS(422) = -B(332)

  JVS(423) = -B(334)

  JVS(424) = -B(330)

  JVS(425) = 0.559*B(359)

  JVS(426) = 0.948*B(363)

  JVS(427) = 0.936*B(361)

  JVS(428) = B(313)+B(315)

  JVS(429) = B(316)+B(318)

  JVS(430) = B(237)

  JVS(431) = 0.205*B(369)+0.488*B(373)

  JVS(432) = 0.079*B(329)+0.126*B(331)+0.187*B(333)+0.24*B(335)

  JVS(433) = -B(95)-B(97)-B(99)-B(101)-B(103)-B(116)-B(132)-B(150)-B(170)-B(192)

  JVS(434) = 0.5*B(345)+0.729*B(347)+0.75*B(349)

  JVS(435) = 0.5*B(337)+0.729*B(339)+0.75*B(341)

  JVS(436) = 0.001*B(377)+0.137*B(379)+0.711*B(381)

  JVS(437) = 0.675*B(289)

  JVS(438) = 0.596*B(306)+0.152*B(308)

  JVS(439) = 0.24*B(336)

  JVS(440) = 0.616*B(240)

  JVS(441) = 0.515*B(305)

  JVS(442) = 0.126*B(332)+0.729*B(340)+0.729*B(348)+0.137*B(380)

  JVS(443) = -B(133)+B(174)

  JVS(444) = B(160)+B(164)-B(171)+B(175)+B(176)+2*B(178)+B(200)

  JVS(445) = -B(151)+B(177)

  JVS(446) = 0

  JVS(447) = -B(100)+B(165)+0.187*B(334)+0.75*B(342)+0.75*B(350)+0.488*B(374)+0.711*B(382)

  JVS(448) = B(238)+0.616*B(241)+0.675*B(290)+0.596*B(307)+B(314)+B(317)+0.079*B(330)+0.5*B(338)+0.5*B(346)+0.559*B(360)&
               &+0.936*B(362)+0.948*B(364)+0.205*B(370)+0.001*B(378)

  JVS(449) = -B(96)+B(161)

  JVS(450) = -B(98)

  JVS(451) = -B(102)

  JVS(452) = -B(104)

  JVS(453) = -B(193)+B(201)

  JVS(454) = -B(117)

  JVS(455) = 0.23*B(329)+0.39*B(331)

  JVS(456) = -B(280)-B(282)-B(284)-B(286)-B(288)

  JVS(457) = 0.025*B(377)+0.026*B(379)+0.012*B(383)

  JVS(458) = -B(287)+0.012*B(384)

  JVS(459) = -B(283)+0.39*B(332)+0.026*B(380)

  JVS(460) = -B(285)

  JVS(461) = -B(281)+0.23*B(330)+0.025*B(378)

  JVS(462) = -B(345)-B(347)-B(349)-B(351)

  JVS(463) = -B(352)

  JVS(464) = -B(348)

  JVS(465) = -B(350)

  JVS(466) = -B(346)

  JVS(467) = -B(337)-B(339)-B(341)-B(343)

  JVS(468) = -B(344)

  JVS(469) = -B(340)

  JVS(470) = -B(342)

  JVS(471) = -B(338)

  JVS(472) = 0.097*B(367)

  JVS(473) = 0.118*B(365)

  JVS(474) = 0.5*B(315)

  JVS(475) = 0.5*B(318)

  JVS(476) = B(311)

  JVS(477) = 0.607*B(357)

  JVS(478) = 0.23*B(265)

  JVS(479) = 0.009*B(327)

  JVS(480) = 0

  JVS(481) = 0.001*B(347)

  JVS(482) = 0.001*B(339)

  JVS(483) = -B(253)-B(254)-B(255)-B(257)

  JVS(484) = 0.15*B(296)+0.023*B(298)

  JVS(485) = 0.009*B(328)

  JVS(486) = 0.023*B(299)+B(312)+0.001*B(340)+0.001*B(348)

  JVS(487) = 0

  JVS(488) = 0

  JVS(489) = 0

  JVS(490) = 0

  JVS(491) = -B(258)

  JVS(492) = -B(256)+0.23*B(266)+0.15*B(297)+0.607*B(358)+0.118*B(366)+0.097*B(368)

  JVS(493) = 0

  JVS(494) = 0

  JVS(495) = 0

  JVS(496) = 0.357*B(329)+0.936*B(333)

  JVS(497) = -B(296)-B(298)-B(300)-B(302)

  JVS(498) = 0.025*B(377)

  JVS(499) = 0

  JVS(500) = -B(299)

  JVS(501) = -B(301)+0.936*B(334)

  JVS(502) = -B(297)+0.357*B(330)+0.025*B(378)

  JVS(503) = -B(377)-B(379)-B(381)-B(383)

  JVS(504) = -B(384)

  JVS(505) = -B(380)

  JVS(506) = -B(382)

  JVS(507) = -B(378)

  JVS(508) = 0.32*B(329)+0.16*B(331)

  JVS(509) = 0.019*B(379)+0.048*B(381)

  JVS(510) = -B(289)-B(291)-B(293)-B(295)

  JVS(511) = -B(294)

  JVS(512) = -B(292)+0.16*B(332)+0.019*B(380)

  JVS(513) = 0.048*B(382)

  JVS(514) = -B(290)+0.32*B(330)

  JVS(515) = B(353)

  JVS(516) = 0.96*B(245)

  JVS(517) = 0.445*B(359)

  JVS(518) = 0.099*B(363)

  JVS(519) = 0.455*B(361)

  JVS(520) = 0.195*B(321)+0.25*B(327)

  JVS(521) = 0.984*B(387)+0.5*B(389)

  JVS(522) = 0.294*B(369)+0.154*B(371)+0.009*B(373)

  JVS(523) = 0.129*B(296)+0.047*B(298)+0.467*B(302)

  JVS(524) = 0.732*B(377)+0.456*B(379)+0.507*B(381)

  JVS(525) = -B(227)-B(229)-B(230)

  JVS(526) = 0.439*B(306)+0.431*B(308)

  JVS(527) = 0.25*B(328)

  JVS(528) = 0.034*B(232)+B(234)

  JVS(529) = 0.482*B(240)+B(242)

  JVS(530) = 0.084*B(303)+0.246*B(305)

  JVS(531) = 0.047*B(299)+0.154*B(372)+0.456*B(380)+0.5*B(390)

  JVS(532) = B(154)

  JVS(533) = B(176)

  JVS(534) = B(140)+B(144)+B(155)+2*B(156)+B(177)+B(198)

  JVS(535) = B(145)-B(231)+0.009*B(374)+0.507*B(382)

  JVS(536) = -B(228)+0.034*B(233)+0.482*B(241)+0.96*B(246)+0.129*B(297)+0.084*B(304)+0.439*B(307)+0.195*B(322)+B(354)&
               &+0.445*B(360)+0.455*B(362)+0.099*B(364)+0.294*B(370)+0.732*B(378)+0.984*B(388)

  JVS(537) = B(141)

  JVS(538) = B(199)

  JVS(539) = 0.081*B(245)

  JVS(540) = 0.026*B(359)

  JVS(541) = 0.026*B(363)

  JVS(542) = 0.35*B(247)+B(249)

  JVS(543) = B(222)

  JVS(544) = B(243)

  JVS(545) = 0.024*B(361)

  JVS(546) = 0.096*B(357)

  JVS(547) = 1.61*B(321)+B(323)+0.191*B(327)

  JVS(548) = B(237)

  JVS(549) = 0.984*B(387)+0.5*B(389)

  JVS(550) = 0.732*B(369)+0.5*B(371)

  JVS(551) = 0.624*B(329)+0.592*B(331)+0.24*B(335)

  JVS(552) = 0.084*B(280)+0.2*B(282)+0.67*B(288)

  JVS(553) = 0.276*B(345)+0.235*B(347)

  JVS(554) = 0.276*B(337)+0.235*B(339)

  JVS(555) = B(254)

  JVS(556) = 0.055*B(296)+0.125*B(298)+0.227*B(300)+0.3*B(302)

  JVS(557) = 0.244*B(377)+0.269*B(379)+0.079*B(381)

  JVS(558) = 0.3*B(289)+0.1*B(291)

  JVS(559) = -B(216)-B(217)-B(218)-B(220)-B(225)

  JVS(560) = 0.01*B(306)+0.134*B(308)

  JVS(561) = 0.191*B(328)+0.24*B(336)

  JVS(562) = 0.115*B(240)

  JVS(563) = 0.213*B(303)+0.506*B(305)

  JVS(564) = 0.2*B(283)+0.1*B(292)+0.125*B(299)+B(324)+0.592*B(332)+0.235*B(340)+0.235*B(348)+0.5*B(372)+0.269*B(380)&
               &+0.5*B(390)

  JVS(565) = B(128)+B(196)

  JVS(566) = B(166)+B(200)

  JVS(567) = B(146)+B(198)

  JVS(568) = 0

  JVS(569) = B(82)+B(186)-B(226)+0.227*B(301)+0.079*B(382)

  JVS(570) = -B(219)+B(238)+0.115*B(241)+B(244)+0.081*B(246)+0.35*B(248)+0.084*B(281)+0.3*B(290)+0.055*B(297)+0.213&
               &*B(304)+0.01*B(307)+1.61*B(322)+0.624*B(330)+0.276*B(338)+0.276*B(346)+0.096*B(358)+0.026*B(360)+0.024&
               &*B(362)+0.026*B(364)+0.732*B(370)+0.244*B(378)+0.984*B(388)

  JVS(571) = B(78)+B(182)

  JVS(572) = -B(221)

  JVS(573) = B(79)+B(83)+B(84)+2*B(85)+0.75*B(92)+0.75*B(110)+B(129)+B(147)+B(167)+B(188)

  JVS(574) = 0.75*B(93)

  JVS(575) = B(183)+B(187)+B(189)+B(197)+B(199)+B(201)+2*B(202)

  JVS(576) = 0.75*B(111)

  JVS(577) = B(203)

  JVS(578) = 0.511*B(373)

  JVS(579) = 0.276*B(349)

  JVS(580) = 0.276*B(341)

  JVS(581) = 0.572*B(300)

  JVS(582) = 0.321*B(381)

  JVS(583) = -0.69*B(306)-B(308)

  JVS(584) = 0

  JVS(585) = 0

  JVS(586) = B(204)

  JVS(587) = 0.572*B(301)+0.276*B(342)+0.276*B(350)+0.511*B(374)+0.321*B(382)

  JVS(588) = -0.69*B(307)

  JVS(589) = B(106)

  JVS(590) = B(107)

  JVS(591) = B(34)

  JVS(592) = -B(327)

  JVS(593) = -B(393)

  JVS(594) = -B(375)

  JVS(595) = -B(335)

  JVS(596) = -B(286)

  JVS(597) = -B(351)

  JVS(598) = -B(343)

  JVS(599) = -B(383)

  JVS(600) = -B(293)

  JVS(601) = -B(2)-B(4)-B(6)-B(9)-B(11)-B(287)-B(294)-B(328)-B(336)-B(344)-B(352)-B(376)-B(384)-B(394)

  JVS(602) = -B(5)+B(30)

  JVS(603) = B(1)-B(10)-B(12)

  JVS(604) = B(29)

  JVS(605) = 0

  JVS(606) = -B(7)

  JVS(607) = 0.261*B(355)

  JVS(608) = 0.122*B(359)

  JVS(609) = 0.204*B(363)

  JVS(610) = 0.244*B(361)

  JVS(611) = B(313)

  JVS(612) = B(316)

  JVS(613) = B(309)

  JVS(614) = B(250)+B(252)

  JVS(615) = B(325)

  JVS(616) = 0.45*B(393)

  JVS(617) = 0.497*B(369)+0.363*B(371)+0.037*B(373)+0.45*B(375)

  JVS(618) = B(286)

  JVS(619) = 0.474*B(345)+0.205*B(347)+0.474*B(349)+0.147*B(351)

  JVS(620) = 0.474*B(337)+0.205*B(339)+0.474*B(341)+0.147*B(343)

  JVS(621) = 0.013*B(296)+0.218*B(300)

  JVS(622) = 0.511*B(377)+0.305*B(379)+0.151*B(381)+0.069*B(383)

  JVS(623) = 0.675*B(289)+0.45*B(293)

  JVS(624) = 0.213*B(306)+0.147*B(308)

  JVS(625) = B(287)+0.45*B(294)+0.147*B(344)+0.147*B(352)+0.45*B(376)+0.069*B(384)+0.45*B(394)

  JVS(626) = -B(232)-B(234)-B(235)

  JVS(627) = 0.37*B(240)

  JVS(628) = 0.558*B(303)+0.71*B(305)

  JVS(629) = 0.205*B(340)+0.205*B(348)+0.363*B(372)+0.305*B(380)

  JVS(630) = 0

  JVS(631) = -B(236)+0.218*B(301)+B(326)+0.474*B(342)+0.474*B(350)+0.037*B(374)+0.151*B(382)

  JVS(632) = -B(233)+0.37*B(241)+B(251)+0.675*B(290)+0.013*B(297)+0.558*B(304)+0.213*B(307)+B(310)+B(314)+B(317)+0.474&
               &*B(338)+0.474*B(346)+0.261*B(356)+0.122*B(360)+0.244*B(362)+0.204*B(364)+0.497*B(370)+0.511*B(378)

  JVS(633) = 0

  JVS(634) = 0

  JVS(635) = 0

  JVS(636) = 0

  JVS(637) = 0.332*B(359)

  JVS(638) = 0.089*B(363)

  JVS(639) = 0.11*B(361)

  JVS(640) = 0.55*B(393)

  JVS(641) = 0.437*B(375)

  JVS(642) = 0.416*B(280)

  JVS(643) = 0.15*B(296)+0.21*B(298)+0.233*B(302)

  JVS(644) = 0.072*B(377)+0.026*B(379)+0.001*B(381)+0.659*B(383)

  JVS(645) = 0.55*B(293)

  JVS(646) = 0.177*B(306)+0.243*B(308)

  JVS(647) = 0.55*B(294)+0.437*B(376)+0.659*B(384)+0.55*B(394)

  JVS(648) = -B(240)-B(242)

  JVS(649) = 0.115*B(303)

  JVS(650) = 0.21*B(299)+0.026*B(380)

  JVS(651) = 0

  JVS(652) = B(112)+0.001*B(382)

  JVS(653) = -B(241)+0.416*B(281)+0.15*B(297)+0.115*B(304)+0.177*B(307)+0.332*B(360)+0.11*B(362)+0.089*B(364)+0.072&
               &*B(378)

  JVS(654) = 0

  JVS(655) = 0.5*B(110)

  JVS(656) = 0.5*B(114)

  JVS(657) = 0.5*B(111)+B(113)+0.5*B(115)+B(118)

  JVS(658) = 0.417*B(363)

  JVS(659) = 0.125*B(361)

  JVS(660) = 0.055*B(365)

  JVS(661) = 0.119*B(369)+0.215*B(371)+0.113*B(375)

  JVS(662) = 0.1*B(331)+0.75*B(335)

  JVS(663) = 0.276*B(345)+0.276*B(347)+0.853*B(351)

  JVS(664) = 0.276*B(337)+0.276*B(339)+0.853*B(343)

  JVS(665) = 0.332*B(296)

  JVS(666) = 0.043*B(379)+0.259*B(383)

  JVS(667) = 0.7*B(295)

  JVS(668) = 0.048*B(306)+0.435*B(308)

  JVS(669) = 0.75*B(336)+0.853*B(344)+0.853*B(352)+0.113*B(376)+0.259*B(384)

  JVS(670) = -0.671*B(303)-B(305)

  JVS(671) = 0.1*B(332)+0.276*B(340)+0.276*B(348)+0.215*B(372)+0.043*B(380)

  JVS(672) = B(134)

  JVS(673) = B(172)

  JVS(674) = B(152)

  JVS(675) = 0

  JVS(676) = 0

  JVS(677) = 0.332*B(297)-0.671*B(304)+0.048*B(307)+0.276*B(338)+0.276*B(346)+0.125*B(362)+0.417*B(364)+0.055*B(366)&
               &+0.119*B(370)

  JVS(678) = 0

  JVS(679) = 0.5*B(110)

  JVS(680) = 0.5*B(114)

  JVS(681) = 0.5*B(111)+0.5*B(115)+B(118)+B(135)+B(153)+B(173)

  JVS(682) = -B(311)

  JVS(683) = -B(385)

  JVS(684) = -B(323)

  JVS(685) = -B(389)

  JVS(686) = -B(371)

  JVS(687) = -B(331)

  JVS(688) = -B(282)

  JVS(689) = -B(347)

  JVS(690) = -B(339)

  JVS(691) = -B(298)

  JVS(692) = -B(379)

  JVS(693) = -B(291)

  JVS(694) = B(2)-B(4)

  JVS(695) = -B(5)-B(13)-B(15)-B(30)-B(31)-B(51)-B(61)-B(283)-B(292)-B(299)-B(312)-B(324)-B(332)-B(340)-B(348)-B(372)&
               &-B(380)-B(386)-B(390)

  JVS(696) = 0.25*B(124)

  JVS(697) = 0.25*B(162)

  JVS(698) = 0.25*B(142)

  JVS(699) = -B(16)

  JVS(700) = 0

  JVS(701) = -B(52)

  JVS(702) = -B(14)

  JVS(703) = -B(62)+0.25*B(125)+0.25*B(143)+0.25*B(163)+0.25*B(184)

  JVS(704) = 0.25*B(185)

  JVS(705) = B(121)

  JVS(706) = 2*B(264)

  JVS(707) = 0

  JVS(708) = 0.011*B(361)

  JVS(709) = B(313)+0.5*B(315)

  JVS(710) = B(316)+0.5*B(318)

  JVS(711) = B(259)+B(260)+B(262)

  JVS(712) = B(237)+B(239)

  JVS(713) = 0

  JVS(714) = 0.67*B(288)

  JVS(715) = 0.123*B(347)

  JVS(716) = 0.123*B(339)

  JVS(717) = 0.467*B(302)

  JVS(718) = 0.137*B(379)

  JVS(719) = 0.675*B(289)

  JVS(720) = B(227)+B(230)

  JVS(721) = 0

  JVS(722) = 0

  JVS(723) = 0

  JVS(724) = 0.492*B(240)+B(242)

  JVS(725) = 0.029*B(303)+0.667*B(305)

  JVS(726) = 0.123*B(340)+0.123*B(348)+0.137*B(380)

  JVS(727) = -B(119)-B(122)-B(124)-B(126)-B(128)-B(130)-B(134)-2*B(136)-B(154)-B(174)

  JVS(728) = -B(175)+B(200)

  JVS(729) = -B(155)+B(198)

  JVS(730) = -B(120)

  JVS(731) = -B(127)+B(186)+B(231)+B(263)

  JVS(732) = B(228)+B(238)+0.492*B(241)+B(261)+0.675*B(290)+0.029*B(304)+B(314)+B(317)+0.011*B(362)

  JVS(733) = -B(123)+B(182)

  JVS(734) = -B(125)

  JVS(735) = -B(129)

  JVS(736) = -B(131)

  JVS(737) = B(183)+B(187)+B(199)+B(201)+2*B(202)

  JVS(738) = -B(135)

  JVS(739) = B(159)

  JVS(740) = B(275)+B(278)

  JVS(741) = 0

  JVS(742) = 0

  JVS(743) = 0

  JVS(744) = -B(174)

  JVS(745) = -B(157)-B(160)-B(162)-B(164)-B(166)-B(168)-B(172)-B(175)-B(176)-2*B(178)-B(200)

  JVS(746) = -B(177)

  JVS(747) = -B(158)

  JVS(748) = -B(165)+B(279)

  JVS(749) = B(276)

  JVS(750) = -B(161)

  JVS(751) = -B(163)

  JVS(752) = -B(167)

  JVS(753) = -B(169)

  JVS(754) = -B(201)

  JVS(755) = -B(173)

  JVS(756) = B(139)

  JVS(757) = 0.1*B(282)

  JVS(758) = 0.201*B(347)

  JVS(759) = 0.201*B(339)

  JVS(760) = 0.37*B(255)+0.37*B(257)

  JVS(761) = 0.048*B(298)+0.3*B(302)

  JVS(762) = 0.006*B(379)

  JVS(763) = 0.05*B(291)

  JVS(764) = 0

  JVS(765) = 0.965*B(232)+B(235)

  JVS(766) = 0.096*B(240)

  JVS(767) = 0.049*B(303)+0.333*B(305)

  JVS(768) = 0.1*B(283)+0.05*B(292)+0.048*B(299)+0.201*B(340)+0.201*B(348)+0.006*B(380)

  JVS(769) = -B(154)

  JVS(770) = -B(176)

  JVS(771) = -B(137)-B(140)-B(142)-B(144)-B(146)-B(148)-B(152)-B(155)-2*B(156)-B(177)-B(198)

  JVS(772) = -B(138)

  JVS(773) = -B(145)+B(236)+0.37*B(258)

  JVS(774) = 0.965*B(233)+0.096*B(241)+0.37*B(256)+0.049*B(304)

  JVS(775) = -B(141)

  JVS(776) = -B(143)

  JVS(777) = -B(147)

  JVS(778) = -B(149)

  JVS(779) = -B(199)

  JVS(780) = -B(153)

  JVS(781) = B(121)

  JVS(782) = B(139)

  JVS(783) = B(159)

  JVS(784) = B(181)

  JVS(785) = B(23)

  JVS(786) = B(39)+B(40)

  JVS(787) = -B(203)

  JVS(788) = B(223)

  JVS(789) = -B(211)

  JVS(790) = B(57)+0.61*B(58)+B(59)

  JVS(791) = 0

  JVS(792) = B(48)

  JVS(793) = -B(206)

  JVS(794) = 0.187*B(333)

  JVS(795) = B(95)+B(99)

  JVS(796) = 0

  JVS(797) = 0.474*B(349)

  JVS(798) = 0.474*B(341)

  JVS(799) = 0

  JVS(800) = 0

  JVS(801) = 0.391*B(381)

  JVS(802) = 0

  JVS(803) = 0

  JVS(804) = 0

  JVS(805) = 0.338*B(306)+B(308)

  JVS(806) = B(6)-B(9)-B(11)

  JVS(807) = 0

  JVS(808) = 0

  JVS(809) = 0

  JVS(810) = B(13)-B(15)

  JVS(811) = -B(119)+B(122)+B(126)

  JVS(812) = -B(157)+B(160)+B(164)

  JVS(813) = -B(137)+B(140)+B(144)

  JVS(814) = -B(1)-B(10)-B(12)-B(16)-B(21)-B(42)-B(55)-B(120)-B(138)-B(158)-B(179)-B(204)-B(207)-B(212)

  JVS(815) = 2*B(17)-B(22)+B(29)+B(44)+0.8*B(66)+2*B(68)+B(82)+B(90)+B(100)+B(112)+B(127)+B(145)+B(165)+B(186)+0.187&
               &*B(334)+0.474*B(342)+0.474*B(350)+0.391*B(382)

  JVS(816) = B(41)-B(43)+B(45)+B(60)+0.338*B(307)

  JVS(817) = B(7)+B(14)+2*B(18)+2*B(19)+B(53)+B(78)+B(86)+B(96)+B(123)+B(141)+B(161)+B(182)+B(224)

  JVS(818) = B(54)-B(56)+0.8*B(67)

  JVS(819) = B(79)+B(83)

  JVS(820) = B(87)+B(91)

  JVS(821) = -B(180)+B(183)+B(187)

  JVS(822) = B(113)

  JVS(823) = B(23)

  JVS(824) = -B(460)

  JVS(825) = 0.39*B(58)

  JVS(826) = -B(271)

  JVS(827) = -B(273)

  JVS(828) = -B(278)

  JVS(829) = -B(267)

  JVS(830) = -B(262)

  JVS(831) = B(46)

  JVS(832) = -B(325)

  JVS(833) = -B(391)

  JVS(834) = 0

  JVS(835) = -B(373)

  JVS(836) = -B(333)

  JVS(837) = -B(99)

  JVS(838) = -B(284)

  JVS(839) = -B(349)

  JVS(840) = -B(341)

  JVS(841) = -B(257)

  JVS(842) = -B(300)

  JVS(843) = -B(381)

  JVS(844) = 0

  JVS(845) = -B(230)

  JVS(846) = -B(225)

  JVS(847) = 0

  JVS(848) = B(11)

  JVS(849) = -B(235)

  JVS(850) = 0

  JVS(851) = 0

  JVS(852) = B(15)

  JVS(853) = -B(126)

  JVS(854) = -B(164)

  JVS(855) = -B(144)

  JVS(856) = B(12)+B(16)-B(21)-B(26)

  JVS(857) = -B(17)-B(22)-B(27)-B(28)-B(29)-B(44)-B(66)-2*B(68)-B(82)-B(90)-B(100)-B(112)-B(127)-B(145)-B(165)-B(186)&
               &-B(226)-B(231)-B(236)-B(258)-B(263)-B(268)-B(272)-B(274)-B(279)-B(285)-B(301)-B(326)-B(334)-B(342)-B(350)&
               &-B(374)-B(382)-B(392)-B(461)

  JVS(858) = -B(45)+B(47)

  JVS(859) = -B(18)

  JVS(860) = -B(67)

  JVS(861) = -B(83)

  JVS(862) = -B(91)

  JVS(863) = -B(187)

  JVS(864) = -B(113)

  JVS(865) = -B(432)

  JVS(866) = -B(440)

  JVS(867) = 2*B(32)

  JVS(868) = -B(448)

  JVS(869) = -B(436)

  JVS(870) = -B(452)

  JVS(871) = -B(444)

  JVS(872) = -B(319)

  JVS(873) = -B(353)

  JVS(874) = 2*B(69)-B(70)

  JVS(875) = -B(355)

  JVS(876) = -B(245)

  JVS(877) = -B(74)

  JVS(878) = -B(456)-B(458)

  JVS(879) = B(38)-B(40)

  JVS(880) = -B(359)

  JVS(881) = -B(363)

  JVS(882) = -0.65*B(247)+B(249)

  JVS(883) = 0.39*B(58)-B(59)

  JVS(884) = -B(367)

  JVS(885) = -B(243)

  JVS(886) = -B(361)

  JVS(887) = -B(365)

  JVS(888) = -B(313)

  JVS(889) = -B(316)

  JVS(890) = -B(269)

  JVS(891) = -B(309)+0.5*B(311)

  JVS(892) = -0.397*B(357)+0.5*B(385)

  JVS(893) = -0.34*B(250)+B(252)

  JVS(894) = -B(275)

  JVS(895) = -B(265)

  JVS(896) = -B(260)

  JVS(897) = -B(49)

  JVS(898) = -B(46)+B(48)

  JVS(899) = -B(321)+0.12*B(323)

  JVS(900) = -B(237)

  JVS(901) = -B(387)+0.32*B(389)

  JVS(902) = 0

  JVS(903) = -B(369)+0.155*B(371)

  JVS(904) = -B(329)+0.266*B(331)

  JVS(905) = -B(280)+0.208*B(282)+0.33*B(288)

  JVS(906) = -B(345)+0.567*B(347)

  JVS(907) = -B(337)+0.567*B(339)

  JVS(908) = -B(255)

  JVS(909) = -B(296)+0.285*B(298)

  JVS(910) = -B(377)+0.378*B(379)

  JVS(911) = -B(289)+0.164*B(291)

  JVS(912) = -B(227)

  JVS(913) = -B(218)

  JVS(914) = -B(306)

  JVS(915) = 0

  JVS(916) = -B(232)

  JVS(917) = -B(240)

  JVS(918) = -B(303)

  JVS(919) = -B(51)+B(61)+0.208*B(283)+0.164*B(292)+0.285*B(299)+0.5*B(312)+0.12*B(324)+0.266*B(332)+0.567*B(340)+0.567&
               &*B(348)+0.155*B(372)+0.378*B(380)+0.5*B(386)+0.32*B(390)

  JVS(920) = 0

  JVS(921) = 0

  JVS(922) = 0

  JVS(923) = -B(42)

  JVS(924) = -B(44)+0.8*B(66)

  JVS(925) = -B(36)-B(41)-B(43)-B(45)-B(47)-B(50)-B(52)-B(60)-B(71)-B(72)-B(75)-B(76)-B(219)-B(228)-B(233)-B(238)-B(241)&
               &-B(244)-B(246)-0.65*B(248)-0.34*B(251)-B(256)-B(261)-B(266)-B(270)-B(276)-B(281)-B(290)-B(297)-B(304)-B(307)&
               &-B(310)-B(314)-B(317)-B(320)-B(322)-B(330)-B(338)-B(346)-B(354)-B(356)-0.397*B(358)-B(360)-B(362)-B(364)&
               &-B(366)-B(368)-B(370)-B(378)-B(388)-B(433)-B(437)-B(441)-B(445)-B(449)-B(453)-B(457)-B(459)

  JVS(926) = -B(37)+B(53)

  JVS(927) = B(54)+B(62)+0.8*B(67)-B(73)

  JVS(928) = 0

  JVS(929) = 0

  JVS(930) = 0

  JVS(931) = 0

  JVS(932) = B(38)

  JVS(933) = -B(223)

  JVS(934) = -B(95)

  JVS(935) = 0

  JVS(936) = 0

  JVS(937) = 0

  JVS(938) = 0

  JVS(939) = 0

  JVS(940) = 0

  JVS(941) = -B(6)+B(9)

  JVS(942) = 0

  JVS(943) = 0

  JVS(944) = -B(13)

  JVS(945) = -B(122)

  JVS(946) = -B(160)

  JVS(947) = -B(140)

  JVS(948) = B(1)+B(10)+B(26)

  JVS(949) = -B(17)+B(27)+B(28)

  JVS(950) = -B(36)

  JVS(951) = -B(7)-B(14)-B(18)-2*B(19)-B(37)-B(53)-B(78)-B(86)-B(96)-B(106)-B(123)-B(141)-B(161)-B(182)-B(224)

  JVS(952) = -B(54)

  JVS(953) = -B(79)

  JVS(954) = -B(87)

  JVS(955) = -B(183)

  JVS(956) = -B(107)

  JVS(957) = B(70)

  JVS(958) = 0.95*B(245)

  JVS(959) = B(74)

  JVS(960) = 0.5*B(458)

  JVS(961) = B(39)

  JVS(962) = B(249)

  JVS(963) = B(222)+B(223)

  JVS(964) = -B(213)

  JVS(965) = B(57)+0.61*B(58)

  JVS(966) = 0.187*B(367)

  JVS(967) = B(243)

  JVS(968) = 0.224*B(365)

  JVS(969) = 0.5*B(315)

  JVS(970) = 0.5*B(318)

  JVS(971) = 1.5*B(311)

  JVS(972) = 0.297*B(357)+1.5*B(385)

  JVS(973) = 0

  JVS(974) = B(252)

  JVS(975) = B(259)

  JVS(976) = B(49)

  JVS(977) = 0.12*B(323)+0.5*B(327)

  JVS(978) = 0.06*B(389)

  JVS(979) = -B(208)

  JVS(980) = 0.056*B(371)

  JVS(981) = 0

  JVS(982) = 0.008*B(282)+0.34*B(288)

  JVS(983) = 0.033*B(347)

  JVS(984) = 0.033*B(339)

  JVS(985) = 2*B(253)+0.63*B(255)+0.63*B(257)

  JVS(986) = 0.4*B(298)+1.233*B(302)

  JVS(987) = 0.003*B(379)+0.013*B(383)

  JVS(988) = 0.064*B(291)

  JVS(989) = B(229)

  JVS(990) = 2*B(216)+B(218)-B(220)+B(225)

  JVS(991) = 0.113*B(306)+0.341*B(308)

  JVS(992) = 0.5*B(328)+0.013*B(384)

  JVS(993) = B(234)

  JVS(994) = 0

  JVS(995) = 0.379*B(303)

  JVS(996) = B(51)-B(61)+0.008*B(283)+0.064*B(292)+0.4*B(299)+1.5*B(312)+0.12*B(324)+0.033*B(340)+0.033*B(348)+0.056&
               &*B(372)+0.003*B(380)+1.5*B(386)+0.06*B(390)

  JVS(997) = -B(124)

  JVS(998) = -B(162)

  JVS(999) = -B(142)

  JVS(1000) = -B(55)

  JVS(1001) = B(44)-B(66)+B(82)+B(90)+B(112)+B(226)+0.63*B(258)

  JVS(1002) = B(45)+B(50)+B(52)+B(71)-B(72)+B(75)+B(76)+B(219)+B(244)+0.95*B(246)+0.63*B(256)+0.379*B(304)+0.113*B(307)&
                &+0.297*B(358)+0.224*B(366)+0.187*B(368)+0.5*B(459)

  JVS(1003) = -B(53)+B(78)+B(86)+B(224)

  JVS(1004) = -B(54)-B(56)-B(62)-2*B(63)-2*B(64)-B(67)-B(73)-B(80)-B(88)-B(108)-B(125)-B(143)-B(163)-B(184)-B(209)&
                &-B(214)-B(221)-B(396)

  JVS(1005) = B(79)-B(81)+B(83)+2*B(85)+B(92)+B(110)

  JVS(1006) = B(87)-B(89)+B(91)+B(93)+B(94)+B(114)

  JVS(1007) = -B(185)

  JVS(1008) = -B(109)+B(111)+B(113)+B(115)+B(118)

  JVS(1009) = B(319)

  JVS(1010) = B(205)

  JVS(1011) = 0.65*B(247)

  JVS(1012) = 0.011*B(361)

  JVS(1013) = 0.3*B(327)

  JVS(1014) = B(239)

  JVS(1015) = 0.26*B(389)

  JVS(1016) = 0.076*B(371)

  JVS(1017) = 0.25*B(335)

  JVS(1018) = 0

  JVS(1019) = 0

  JVS(1020) = 0.197*B(379)+0.03*B(381)

  JVS(1021) = 0.3*B(295)

  JVS(1022) = B(229)

  JVS(1023) = 0

  JVS(1024) = 0.3*B(328)+0.25*B(336)

  JVS(1025) = 0

  JVS(1026) = 0

  JVS(1027) = 0

  JVS(1028) = 0.076*B(372)+0.197*B(380)+0.26*B(390)

  JVS(1029) = B(122)+B(126)-B(128)+2*B(136)+B(154)+B(174)+B(196)

  JVS(1030) = -B(166)+B(175)

  JVS(1031) = -B(146)+B(155)

  JVS(1032) = 0

  JVS(1033) = -B(82)+B(127)+0.03*B(382)

  JVS(1034) = 0.65*B(248)+B(320)+0.011*B(362)

  JVS(1035) = -B(78)+B(123)

  JVS(1036) = -B(80)

  JVS(1037) = -B(79)-B(81)-B(83)-2*B(84)-2*B(85)-B(92)-B(110)-B(129)-B(147)-B(167)-B(188)

  JVS(1038) = -B(93)

  JVS(1039) = -B(189)+B(197)

  JVS(1040) = -B(111)

  JVS(1041) = B(353)

  JVS(1042) = 0.965*B(355)

  JVS(1043) = 0.05*B(245)

  JVS(1044) = 0.695*B(359)

  JVS(1045) = 0.653*B(363)

  JVS(1046) = 0.804*B(367)

  JVS(1047) = 0.835*B(361)

  JVS(1048) = 0.765*B(365)

  JVS(1049) = B(315)

  JVS(1050) = B(318)

  JVS(1051) = 0.76*B(269)

  JVS(1052) = B(309)

  JVS(1053) = 0.1*B(357)

  JVS(1054) = 0.34*B(250)

  JVS(1055) = 0.76*B(265)

  JVS(1056) = B(321)+B(325)+0.2*B(327)

  JVS(1057) = 0.984*B(387)+0.949*B(391)

  JVS(1058) = 0

  JVS(1059) = 0.91*B(369)+0.022*B(371)+0.824*B(373)

  JVS(1060) = 0.907*B(329)+0.066*B(331)+0.749*B(333)

  JVS(1061) = 0.5*B(280)+0.1*B(282)+0.5*B(284)+0.33*B(288)

  JVS(1062) = 0.75*B(345)+0.031*B(347)+0.276*B(349)

  JVS(1063) = 0.75*B(337)+0.031*B(339)+0.276*B(341)

  JVS(1064) = 0.67*B(296)+0.048*B(298)+0.799*B(300)

  JVS(1065) = 0.918*B(377)+0.033*B(379)+0.442*B(381)+0.012*B(383)

  JVS(1066) = 0.3*B(289)+0.05*B(291)

  JVS(1067) = 0.376*B(306)+0.564*B(308)

  JVS(1068) = 0.2*B(328)+0.012*B(384)

  JVS(1069) = 0.034*B(232)+B(234)

  JVS(1070) = 0.37*B(240)+B(242)

  JVS(1071) = 0.473*B(303)+0.96*B(305)

  JVS(1072) = 0.1*B(283)+0.05*B(292)+0.048*B(299)+0.066*B(332)+0.031*B(340)+0.031*B(348)+0.022*B(372)+0.033*B(380)

  JVS(1073) = -B(130)+B(154)

  JVS(1074) = -B(168)+B(176)

  JVS(1075) = B(140)+B(144)-B(148)+B(155)+2*B(156)+B(177)+B(198)

  JVS(1076) = 0

  JVS(1077) = -B(90)+B(145)+0.5*B(285)+0.799*B(301)+B(326)+0.749*B(334)+0.276*B(342)+0.276*B(350)+0.824*B(374)+0.442&
                &*B(382)+0.949*B(392)

  JVS(1078) = 0.034*B(233)+0.37*B(241)+0.05*B(246)+0.34*B(251)+0.76*B(266)+0.76*B(270)+0.5*B(281)+0.3*B(290)+0.67*B(297)&
                &+0.473*B(304)+0.376*B(307)+B(310)+B(322)+0.907*B(330)+0.75*B(338)+0.75*B(346)+B(354)+0.965*B(356)+0.1&
                &*B(358)+0.695*B(360)+0.835*B(362)+0.653*B(364)+0.765*B(366)+0.804*B(368)+0.91*B(370)+0.918*B(378)+0.984&
                &*B(388)

  JVS(1079) = -B(86)+B(141)

  JVS(1080) = -B(88)

  JVS(1081) = -B(92)

  JVS(1082) = -B(87)-B(89)-B(91)-B(93)-2*B(94)-B(114)-B(131)-B(149)-B(169)-B(190)

  JVS(1083) = -B(191)+B(199)

  JVS(1084) = -B(115)

  JVS(1085) = B(181)

  JVS(1086) = 0.192*B(331)+0.24*B(335)

  JVS(1087) = 0.5*B(280)+0.5*B(284)+0.33*B(288)

  JVS(1088) = 0.289*B(296)+0.15*B(300)

  JVS(1089) = 0

  JVS(1090) = 0.3*B(295)

  JVS(1091) = 0.24*B(336)

  JVS(1092) = 0.192*B(332)

  JVS(1093) = -B(196)

  JVS(1094) = -B(200)

  JVS(1095) = -B(198)

  JVS(1096) = -B(179)

  JVS(1097) = -B(186)+0.5*B(285)+0.15*B(301)

  JVS(1098) = 0.5*B(281)+0.289*B(297)

  JVS(1099) = -B(182)

  JVS(1100) = -B(184)

  JVS(1101) = -B(188)

  JVS(1102) = -B(190)

  JVS(1103) = -B(180)-B(183)-B(185)-B(187)-B(189)-B(191)-B(194)-B(197)-B(199)-B(201)-2*B(202)

  JVS(1104) = -B(195)

  JVS(1105) = 0.035*B(355)

  JVS(1106) = 0.07*B(359)

  JVS(1107) = 0.347*B(363)

  JVS(1108) = 0.009*B(367)

  JVS(1109) = 0.143*B(361)

  JVS(1110) = 0.011*B(365)

  JVS(1111) = 0.016*B(387)+0.051*B(391)

  JVS(1112) = 0.09*B(369)+0.001*B(371)+0.176*B(373)

  JVS(1113) = 0.093*B(329)+0.008*B(331)+0.064*B(333)+0.01*B(335)

  JVS(1114) = 0.25*B(345)+0.18*B(347)+0.25*B(349)

  JVS(1115) = 0.25*B(337)+0.18*B(339)+0.25*B(341)

  JVS(1116) = 0.041*B(296)+0.051*B(300)

  JVS(1117) = 0.082*B(377)+0.002*B(379)+0.136*B(381)+0.001*B(383)

  JVS(1118) = 0.025*B(289)

  JVS(1119) = 0.173*B(306)+0.095*B(308)

  JVS(1120) = 0.01*B(336)+0.001*B(384)

  JVS(1121) = 0.001*B(232)

  JVS(1122) = 0.042*B(240)

  JVS(1123) = 0.07*B(303)+0.04*B(305)

  JVS(1124) = 0.008*B(332)+0.18*B(340)+0.18*B(348)+0.001*B(372)+0.002*B(380)

  JVS(1125) = -B(134)

  JVS(1126) = -B(172)

  JVS(1127) = -B(152)

  JVS(1128) = 0

  JVS(1129) = -B(112)+0.051*B(301)+0.064*B(334)+0.25*B(342)+0.25*B(350)+0.176*B(374)+0.136*B(382)+0.051*B(392)

  JVS(1130) = 0.001*B(233)+0.042*B(241)+0.025*B(290)+0.041*B(297)+0.07*B(304)+0.173*B(307)+0.093*B(330)+0.25*B(338)+0.25&
                &*B(346)+0.035*B(356)+0.07*B(360)+0.143*B(362)+0.347*B(364)+0.011*B(366)+0.009*B(368)+0.09*B(370)+0.082&
                &*B(378)+0.016*B(388)

  JVS(1131) = -B(106)

  JVS(1132) = -B(108)

  JVS(1133) = -B(110)

  JVS(1134) = -B(114)

  JVS(1135) = -B(194)

  JVS(1136) = -B(107)-B(109)-B(111)-B(113)-B(115)-2*B(118)-B(135)-B(153)-B(173)-B(195)
      
END SUBROUTINE saprc99_mosaic_20bin_vbs2_aq_Jac_SP














SUBROUTINE saprc99_mosaic_20bin_vbs2_aq_KppDecomp( JVS, IER )







      INTEGER  :: IER
      REAL(kind=dp) :: JVS(1136), W(101), a
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
      
END SUBROUTINE saprc99_mosaic_20bin_vbs2_aq_KppDecomp



SUBROUTINE saprc99_mosaic_20bin_vbs2_aq_KppDecompCmplx( JVS, IER )







      INTEGER  :: IER
      DOUBLE COMPLEX :: JVS(1136), W(101), a
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
      
END SUBROUTINE saprc99_mosaic_20bin_vbs2_aq_KppDecompCmplx


SUBROUTINE saprc99_mosaic_20bin_vbs2_aq_KppSolveIndirect( JVS, X )







      INTEGER i, j
      REAL(kind=dp) JVS(1136), X(101), sum

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
      
END SUBROUTINE saprc99_mosaic_20bin_vbs2_aq_KppSolveIndirect


SUBROUTINE saprc99_mosaic_20bin_vbs2_aq_KppSolveCmplx( JVS, X )







      INTEGER i, j
      DOUBLE COMPLEX JVS(1136), X(101), sum

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
      
END SUBROUTINE saprc99_mosaic_20bin_vbs2_aq_KppSolveCmplx













SUBROUTINE saprc99_mosaic_20bin_vbs2_aq_KppSolve ( JVS, X )


  REAL(kind=dp) :: JVS(LU_NONZERO)

  REAL(kind=dp) :: X(NVAR)

  X(21) = X(21)-JVS(158)*X(3)-JVS(159)*X(4)-JVS(160)*X(5)-JVS(161)*X(6)
  X(28) = X(28)-JVS(188)*X(27)
  X(30) = X(30)-JVS(193)*X(29)
  X(46) = X(46)-JVS(240)*X(45)
  X(56) = X(56)-JVS(274)*X(52)-JVS(275)*X(55)
  X(57) = X(57)-JVS(278)*X(52)-JVS(279)*X(55)
  X(58) = X(58)-JVS(282)*X(52)-JVS(283)*X(55)
  X(59) = X(59)-JVS(287)*X(52)-JVS(288)*X(55)
  X(61) = X(61)-JVS(295)*X(50)
  X(63) = X(63)-JVS(306)*X(52)-JVS(307)*X(55)
  X(64) = X(64)-JVS(313)*X(55)
  X(65) = X(65)-JVS(319)*X(52)-JVS(320)*X(55)-JVS(321)*X(56)-JVS(322)*X(57)-JVS(323)*X(58)
  X(66) = X(66)-JVS(331)*X(54)-JVS(332)*X(56)-JVS(333)*X(57)-JVS(334)*X(59)-JVS(335)*X(60)-JVS(336)*X(65)
  X(67) = X(67)-JVS(356)*X(41)-JVS(357)*X(43)-JVS(358)*X(58)-JVS(359)*X(61)-JVS(360)*X(63)-JVS(361)*X(64)-JVS(362)*X(65)
  X(69) = X(69)-JVS(383)*X(38)-JVS(384)*X(45)-JVS(385)*X(46)-JVS(386)*X(47)-JVS(387)*X(54)
  X(71) = X(71)-JVS(403)*X(58)-JVS(404)*X(64)
  X(74) = X(74)-JVS(425)*X(45)-JVS(426)*X(47)-JVS(427)*X(54)-JVS(428)*X(56)-JVS(429)*X(57)-JVS(430)*X(69)-JVS(431)*X(72)&
            &-JVS(432)*X(73)
  X(75) = X(75)-JVS(455)*X(73)
  X(78) = X(78)-JVS(472)*X(52)-JVS(473)*X(55)-JVS(474)*X(56)-JVS(475)*X(57)-JVS(476)*X(59)-JVS(477)*X(60)-JVS(478)*X(64)&
            &-JVS(479)*X(68)-JVS(480)*X(71)-JVS(481)*X(76)-JVS(482)*X(77)
  X(79) = X(79)-JVS(496)*X(73)
  X(81) = X(81)-JVS(508)*X(73)-JVS(509)*X(80)
  X(82) = X(82)-JVS(515)*X(32)-JVS(516)*X(40)-JVS(517)*X(45)-JVS(518)*X(47)-JVS(519)*X(54)-JVS(520)*X(68)-JVS(521)*X(70)&
            &-JVS(522)*X(72)-JVS(523)*X(79)-JVS(524)*X(80)
  X(83) = X(83)-JVS(539)*X(40)-JVS(540)*X(45)-JVS(541)*X(47)-JVS(542)*X(48)-JVS(543)*X(49)-JVS(544)*X(53)-JVS(545)*X(54)&
            &-JVS(546)*X(60)-JVS(547)*X(68)-JVS(548)*X(69)-JVS(549)*X(70)-JVS(550)*X(72)-JVS(551)*X(73)-JVS(552)*X(75)&
            &-JVS(553)*X(76)-JVS(554)*X(77)-JVS(555)*X(78)-JVS(556)*X(79)-JVS(557)*X(80)-JVS(558)*X(81)
  X(84) = X(84)-JVS(577)*X(46)-JVS(578)*X(72)-JVS(579)*X(76)-JVS(580)*X(77)-JVS(581)*X(79)-JVS(582)*X(80)
  X(85) = X(85)-JVS(591)*X(26)-JVS(592)*X(68)-JVS(593)*X(70)-JVS(594)*X(72)-JVS(595)*X(73)-JVS(596)*X(75)-JVS(597)*X(76)&
            &-JVS(598)*X(77)-JVS(599)*X(80)-JVS(600)*X(81)
  X(86) = X(86)-JVS(607)*X(38)-JVS(608)*X(45)-JVS(609)*X(47)-JVS(610)*X(54)-JVS(611)*X(56)-JVS(612)*X(57)-JVS(613)*X(59)&
            &-JVS(614)*X(62)-JVS(615)*X(68)-JVS(616)*X(70)-JVS(617)*X(72)-JVS(618)*X(75)-JVS(619)*X(76)-JVS(620)*X(77)&
            &-JVS(621)*X(79)-JVS(622)*X(80)-JVS(623)*X(81)-JVS(624)*X(84)-JVS(625)*X(85)
  X(87) = X(87)-JVS(637)*X(45)-JVS(638)*X(47)-JVS(639)*X(54)-JVS(640)*X(70)-JVS(641)*X(72)-JVS(642)*X(75)-JVS(643)*X(79)&
            &-JVS(644)*X(80)-JVS(645)*X(81)-JVS(646)*X(84)-JVS(647)*X(85)
  X(88) = X(88)-JVS(658)*X(47)-JVS(659)*X(54)-JVS(660)*X(55)-JVS(661)*X(72)-JVS(662)*X(73)-JVS(663)*X(76)-JVS(664)*X(77)&
            &-JVS(665)*X(79)-JVS(666)*X(80)-JVS(667)*X(81)-JVS(668)*X(84)-JVS(669)*X(85)
  X(89) = X(89)-JVS(682)*X(59)-JVS(683)*X(60)-JVS(684)*X(68)-JVS(685)*X(70)-JVS(686)*X(72)-JVS(687)*X(73)-JVS(688)*X(75)&
            &-JVS(689)*X(76)-JVS(690)*X(77)-JVS(691)*X(79)-JVS(692)*X(80)-JVS(693)*X(81)-JVS(694)*X(85)
  X(90) = X(90)-JVS(705)*X(33)-JVS(706)*X(39)-JVS(707)*X(52)-JVS(708)*X(54)-JVS(709)*X(56)-JVS(710)*X(57)-JVS(711)*X(65)&
            &-JVS(712)*X(69)-JVS(713)*X(72)-JVS(714)*X(75)-JVS(715)*X(76)-JVS(716)*X(77)-JVS(717)*X(79)-JVS(718)*X(80)&
            &-JVS(719)*X(81)-JVS(720)*X(82)-JVS(721)*X(84)-JVS(722)*X(85)-JVS(723)*X(86)-JVS(724)*X(87)-JVS(725)*X(88)&
            &-JVS(726)*X(89)
  X(91) = X(91)-JVS(739)*X(35)-JVS(740)*X(63)-JVS(741)*X(80)-JVS(742)*X(85)-JVS(743)*X(89)-JVS(744)*X(90)
  X(92) = X(92)-JVS(756)*X(34)-JVS(757)*X(75)-JVS(758)*X(76)-JVS(759)*X(77)-JVS(760)*X(78)-JVS(761)*X(79)-JVS(762)*X(80)&
            &-JVS(763)*X(81)-JVS(764)*X(85)-JVS(765)*X(86)-JVS(766)*X(87)-JVS(767)*X(88)-JVS(768)*X(89)-JVS(769)*X(90)&
            &-JVS(770)*X(91)
  X(93) = X(93)-JVS(781)*X(33)-JVS(782)*X(34)-JVS(783)*X(35)-JVS(784)*X(36)-JVS(785)*X(41)-JVS(786)*X(44)-JVS(787)*X(46)&
            &-JVS(788)*X(49)-JVS(789)*X(50)-JVS(790)*X(51)-JVS(791)*X(61)-JVS(792)*X(67)-JVS(793)*X(71)-JVS(794)*X(73)&
            &-JVS(795)*X(74)-JVS(796)*X(75)-JVS(797)*X(76)-JVS(798)*X(77)-JVS(799)*X(78)-JVS(800)*X(79)-JVS(801)*X(80)&
            &-JVS(802)*X(81)-JVS(803)*X(82)-JVS(804)*X(83)-JVS(805)*X(84)-JVS(806)*X(85)-JVS(807)*X(86)-JVS(808)*X(87)&
            &-JVS(809)*X(88)-JVS(810)*X(89)-JVS(811)*X(90)-JVS(812)*X(91)-JVS(813)*X(92)
  X(94) = X(94)-JVS(823)*X(41)-JVS(824)*X(43)-JVS(825)*X(51)-JVS(826)*X(58)-JVS(827)*X(61)-JVS(828)*X(63)-JVS(829)*X(64)&
            &-JVS(830)*X(65)-JVS(831)*X(67)-JVS(832)*X(68)-JVS(833)*X(70)-JVS(834)*X(71)-JVS(835)*X(72)-JVS(836)*X(73)&
            &-JVS(837)*X(74)-JVS(838)*X(75)-JVS(839)*X(76)-JVS(840)*X(77)-JVS(841)*X(78)-JVS(842)*X(79)-JVS(843)*X(80)&
            &-JVS(844)*X(81)-JVS(845)*X(82)-JVS(846)*X(83)-JVS(847)*X(84)-JVS(848)*X(85)-JVS(849)*X(86)-JVS(850)*X(87)&
            &-JVS(851)*X(88)-JVS(852)*X(89)-JVS(853)*X(90)-JVS(854)*X(91)-JVS(855)*X(92)-JVS(856)*X(93)
  X(95) = X(95)-JVS(865)*X(3)-JVS(866)*X(5)-JVS(867)*X(26)-JVS(868)*X(27)-JVS(869)*X(28)-JVS(870)*X(29)-JVS(871)*X(30)&
            &-JVS(872)*X(31)-JVS(873)*X(32)-JVS(874)*X(37)-JVS(875)*X(38)-JVS(876)*X(40)-JVS(877)*X(42)-JVS(878)*X(43)&
            &-JVS(879)*X(44)-JVS(880)*X(45)-JVS(881)*X(47)-JVS(882)*X(48)-JVS(883)*X(51)-JVS(884)*X(52)-JVS(885)*X(53)&
            &-JVS(886)*X(54)-JVS(887)*X(55)-JVS(888)*X(56)-JVS(889)*X(57)-JVS(890)*X(58)-JVS(891)*X(59)-JVS(892)*X(60)&
            &-JVS(893)*X(62)-JVS(894)*X(63)-JVS(895)*X(64)-JVS(896)*X(65)-JVS(897)*X(66)-JVS(898)*X(67)-JVS(899)*X(68)&
            &-JVS(900)*X(69)-JVS(901)*X(70)-JVS(902)*X(71)-JVS(903)*X(72)-JVS(904)*X(73)-JVS(905)*X(75)-JVS(906)*X(76)&
            &-JVS(907)*X(77)-JVS(908)*X(78)-JVS(909)*X(79)-JVS(910)*X(80)-JVS(911)*X(81)-JVS(912)*X(82)-JVS(913)*X(83)&
            &-JVS(914)*X(84)-JVS(915)*X(85)-JVS(916)*X(86)-JVS(917)*X(87)-JVS(918)*X(88)-JVS(919)*X(89)-JVS(920)*X(90)&
            &-JVS(921)*X(91)-JVS(922)*X(92)-JVS(923)*X(93)-JVS(924)*X(94)
  X(96) = X(96)-JVS(932)*X(44)-JVS(933)*X(49)-JVS(934)*X(74)-JVS(935)*X(76)-JVS(936)*X(77)-JVS(937)*X(80)-JVS(938)*X(81)&
            &-JVS(939)*X(83)-JVS(940)*X(84)-JVS(941)*X(85)-JVS(942)*X(87)-JVS(943)*X(88)-JVS(944)*X(89)-JVS(945)*X(90)&
            &-JVS(946)*X(91)-JVS(947)*X(92)-JVS(948)*X(93)-JVS(949)*X(94)-JVS(950)*X(95)
  X(97) = X(97)-JVS(957)*X(37)-JVS(958)*X(40)-JVS(959)*X(42)-JVS(960)*X(43)-JVS(961)*X(44)-JVS(962)*X(48)-JVS(963)*X(49)&
            &-JVS(964)*X(50)-JVS(965)*X(51)-JVS(966)*X(52)-JVS(967)*X(53)-JVS(968)*X(55)-JVS(969)*X(56)-JVS(970)*X(57)&
            &-JVS(971)*X(59)-JVS(972)*X(60)-JVS(973)*X(61)-JVS(974)*X(62)-JVS(975)*X(65)-JVS(976)*X(66)-JVS(977)*X(68)&
            &-JVS(978)*X(70)-JVS(979)*X(71)-JVS(980)*X(72)-JVS(981)*X(73)-JVS(982)*X(75)-JVS(983)*X(76)-JVS(984)*X(77)&
            &-JVS(985)*X(78)-JVS(986)*X(79)-JVS(987)*X(80)-JVS(988)*X(81)-JVS(989)*X(82)-JVS(990)*X(83)-JVS(991)*X(84)&
            &-JVS(992)*X(85)-JVS(993)*X(86)-JVS(994)*X(87)-JVS(995)*X(88)-JVS(996)*X(89)-JVS(997)*X(90)-JVS(998)*X(91)&
            &-JVS(999)*X(92)-JVS(1000)*X(93)-JVS(1001)*X(94)-JVS(1002)*X(95)-JVS(1003)*X(96)
  X(98) = X(98)-JVS(1009)*X(31)-JVS(1010)*X(46)-JVS(1011)*X(48)-JVS(1012)*X(54)-JVS(1013)*X(68)-JVS(1014)*X(69)&
            &-JVS(1015)*X(70)-JVS(1016)*X(72)-JVS(1017)*X(73)-JVS(1018)*X(76)-JVS(1019)*X(77)-JVS(1020)*X(80)-JVS(1021)&
            &*X(81)-JVS(1022)*X(82)-JVS(1023)*X(84)-JVS(1024)*X(85)-JVS(1025)*X(86)-JVS(1026)*X(87)-JVS(1027)*X(88)&
            &-JVS(1028)*X(89)-JVS(1029)*X(90)-JVS(1030)*X(91)-JVS(1031)*X(92)-JVS(1032)*X(93)-JVS(1033)*X(94)-JVS(1034)&
            &*X(95)-JVS(1035)*X(96)-JVS(1036)*X(97)
  X(99) = X(99)-JVS(1041)*X(32)-JVS(1042)*X(38)-JVS(1043)*X(40)-JVS(1044)*X(45)-JVS(1045)*X(47)-JVS(1046)*X(52)&
            &-JVS(1047)*X(54)-JVS(1048)*X(55)-JVS(1049)*X(56)-JVS(1050)*X(57)-JVS(1051)*X(58)-JVS(1052)*X(59)-JVS(1053)&
            &*X(60)-JVS(1054)*X(62)-JVS(1055)*X(64)-JVS(1056)*X(68)-JVS(1057)*X(70)-JVS(1058)*X(71)-JVS(1059)*X(72)&
            &-JVS(1060)*X(73)-JVS(1061)*X(75)-JVS(1062)*X(76)-JVS(1063)*X(77)-JVS(1064)*X(79)-JVS(1065)*X(80)-JVS(1066)&
            &*X(81)-JVS(1067)*X(84)-JVS(1068)*X(85)-JVS(1069)*X(86)-JVS(1070)*X(87)-JVS(1071)*X(88)-JVS(1072)*X(89)&
            &-JVS(1073)*X(90)-JVS(1074)*X(91)-JVS(1075)*X(92)-JVS(1076)*X(93)-JVS(1077)*X(94)-JVS(1078)*X(95)-JVS(1079)&
            &*X(96)-JVS(1080)*X(97)-JVS(1081)*X(98)
  X(100) = X(100)-JVS(1085)*X(36)-JVS(1086)*X(73)-JVS(1087)*X(75)-JVS(1088)*X(79)-JVS(1089)*X(80)-JVS(1090)*X(81)&
             &-JVS(1091)*X(85)-JVS(1092)*X(89)-JVS(1093)*X(90)-JVS(1094)*X(91)-JVS(1095)*X(92)-JVS(1096)*X(93)-JVS(1097)&
             &*X(94)-JVS(1098)*X(95)-JVS(1099)*X(96)-JVS(1100)*X(97)-JVS(1101)*X(98)-JVS(1102)*X(99)
  X(101) = X(101)-JVS(1105)*X(38)-JVS(1106)*X(45)-JVS(1107)*X(47)-JVS(1108)*X(52)-JVS(1109)*X(54)-JVS(1110)*X(55)&
             &-JVS(1111)*X(70)-JVS(1112)*X(72)-JVS(1113)*X(73)-JVS(1114)*X(76)-JVS(1115)*X(77)-JVS(1116)*X(79)-JVS(1117)&
             &*X(80)-JVS(1118)*X(81)-JVS(1119)*X(84)-JVS(1120)*X(85)-JVS(1121)*X(86)-JVS(1122)*X(87)-JVS(1123)*X(88)&
             &-JVS(1124)*X(89)-JVS(1125)*X(90)-JVS(1126)*X(91)-JVS(1127)*X(92)-JVS(1128)*X(93)-JVS(1129)*X(94)-JVS(1130)&
             &*X(95)-JVS(1131)*X(96)-JVS(1132)*X(97)-JVS(1133)*X(98)-JVS(1134)*X(99)-JVS(1135)*X(100)
  X(101) = X(101)/JVS(1136)
  X(100) = (X(100)-JVS(1104)*X(101))/(JVS(1103))
  X(99) = (X(99)-JVS(1083)*X(100)-JVS(1084)*X(101))/(JVS(1082))
  X(98) = (X(98)-JVS(1038)*X(99)-JVS(1039)*X(100)-JVS(1040)*X(101))/(JVS(1037))
  X(97) = (X(97)-JVS(1005)*X(98)-JVS(1006)*X(99)-JVS(1007)*X(100)-JVS(1008)*X(101))/(JVS(1004))
  X(96) = (X(96)-JVS(952)*X(97)-JVS(953)*X(98)-JVS(954)*X(99)-JVS(955)*X(100)-JVS(956)*X(101))/(JVS(951))
  X(95) = (X(95)-JVS(926)*X(96)-JVS(927)*X(97)-JVS(928)*X(98)-JVS(929)*X(99)-JVS(930)*X(100)-JVS(931)*X(101))/(JVS(925))
  X(94) = (X(94)-JVS(858)*X(95)-JVS(859)*X(96)-JVS(860)*X(97)-JVS(861)*X(98)-JVS(862)*X(99)-JVS(863)*X(100)-JVS(864)&
            &*X(101))/(JVS(857))
  X(93) = (X(93)-JVS(815)*X(94)-JVS(816)*X(95)-JVS(817)*X(96)-JVS(818)*X(97)-JVS(819)*X(98)-JVS(820)*X(99)-JVS(821)&
            &*X(100)-JVS(822)*X(101))/(JVS(814))
  X(92) = (X(92)-JVS(772)*X(93)-JVS(773)*X(94)-JVS(774)*X(95)-JVS(775)*X(96)-JVS(776)*X(97)-JVS(777)*X(98)-JVS(778)&
            &*X(99)-JVS(779)*X(100)-JVS(780)*X(101))/(JVS(771))
  X(91) = (X(91)-JVS(746)*X(92)-JVS(747)*X(93)-JVS(748)*X(94)-JVS(749)*X(95)-JVS(750)*X(96)-JVS(751)*X(97)-JVS(752)&
            &*X(98)-JVS(753)*X(99)-JVS(754)*X(100)-JVS(755)*X(101))/(JVS(745))
  X(90) = (X(90)-JVS(728)*X(91)-JVS(729)*X(92)-JVS(730)*X(93)-JVS(731)*X(94)-JVS(732)*X(95)-JVS(733)*X(96)-JVS(734)&
            &*X(97)-JVS(735)*X(98)-JVS(736)*X(99)-JVS(737)*X(100)-JVS(738)*X(101))/(JVS(727))
  X(89) = (X(89)-JVS(696)*X(90)-JVS(697)*X(91)-JVS(698)*X(92)-JVS(699)*X(93)-JVS(700)*X(94)-JVS(701)*X(95)-JVS(702)&
            &*X(96)-JVS(703)*X(97)-JVS(704)*X(100))/(JVS(695))
  X(88) = (X(88)-JVS(671)*X(89)-JVS(672)*X(90)-JVS(673)*X(91)-JVS(674)*X(92)-JVS(675)*X(93)-JVS(676)*X(94)-JVS(677)&
            &*X(95)-JVS(678)*X(96)-JVS(679)*X(98)-JVS(680)*X(99)-JVS(681)*X(101))/(JVS(670))
  X(87) = (X(87)-JVS(649)*X(88)-JVS(650)*X(89)-JVS(651)*X(93)-JVS(652)*X(94)-JVS(653)*X(95)-JVS(654)*X(96)-JVS(655)&
            &*X(98)-JVS(656)*X(99)-JVS(657)*X(101))/(JVS(648))
  X(86) = (X(86)-JVS(627)*X(87)-JVS(628)*X(88)-JVS(629)*X(89)-JVS(630)*X(93)-JVS(631)*X(94)-JVS(632)*X(95)-JVS(633)&
            &*X(96)-JVS(634)*X(97)-JVS(635)*X(99)-JVS(636)*X(101))/(JVS(626))
  X(85) = (X(85)-JVS(602)*X(89)-JVS(603)*X(93)-JVS(604)*X(94)-JVS(605)*X(95)-JVS(606)*X(96))/(JVS(601))
  X(84) = (X(84)-JVS(584)*X(85)-JVS(585)*X(89)-JVS(586)*X(93)-JVS(587)*X(94)-JVS(588)*X(95)-JVS(589)*X(96)-JVS(590)&
            &*X(101))/(JVS(583))
  X(83) = (X(83)-JVS(560)*X(84)-JVS(561)*X(85)-JVS(562)*X(87)-JVS(563)*X(88)-JVS(564)*X(89)-JVS(565)*X(90)-JVS(566)&
            &*X(91)-JVS(567)*X(92)-JVS(568)*X(93)-JVS(569)*X(94)-JVS(570)*X(95)-JVS(571)*X(96)-JVS(572)*X(97)-JVS(573)*X(98)&
            &-JVS(574)*X(99)-JVS(575)*X(100)-JVS(576)*X(101))/(JVS(559))
  X(82) = (X(82)-JVS(526)*X(84)-JVS(527)*X(85)-JVS(528)*X(86)-JVS(529)*X(87)-JVS(530)*X(88)-JVS(531)*X(89)-JVS(532)&
            &*X(90)-JVS(533)*X(91)-JVS(534)*X(92)-JVS(535)*X(94)-JVS(536)*X(95)-JVS(537)*X(96)-JVS(538)*X(100))/(JVS(525))
  X(81) = (X(81)-JVS(511)*X(85)-JVS(512)*X(89)-JVS(513)*X(94)-JVS(514)*X(95))/(JVS(510))
  X(80) = (X(80)-JVS(504)*X(85)-JVS(505)*X(89)-JVS(506)*X(94)-JVS(507)*X(95))/(JVS(503))
  X(79) = (X(79)-JVS(498)*X(80)-JVS(499)*X(85)-JVS(500)*X(89)-JVS(501)*X(94)-JVS(502)*X(95))/(JVS(497))
  X(78) = (X(78)-JVS(484)*X(79)-JVS(485)*X(85)-JVS(486)*X(89)-JVS(487)*X(90)-JVS(488)*X(91)-JVS(489)*X(92)-JVS(490)&
            &*X(93)-JVS(491)*X(94)-JVS(492)*X(95)-JVS(493)*X(96)-JVS(494)*X(97)-JVS(495)*X(100))/(JVS(483))
  X(77) = (X(77)-JVS(468)*X(85)-JVS(469)*X(89)-JVS(470)*X(94)-JVS(471)*X(95))/(JVS(467))
  X(76) = (X(76)-JVS(463)*X(85)-JVS(464)*X(89)-JVS(465)*X(94)-JVS(466)*X(95))/(JVS(462))
  X(75) = (X(75)-JVS(457)*X(80)-JVS(458)*X(85)-JVS(459)*X(89)-JVS(460)*X(94)-JVS(461)*X(95))/(JVS(456))
  X(74) = (X(74)-JVS(434)*X(76)-JVS(435)*X(77)-JVS(436)*X(80)-JVS(437)*X(81)-JVS(438)*X(84)-JVS(439)*X(85)-JVS(440)&
            &*X(87)-JVS(441)*X(88)-JVS(442)*X(89)-JVS(443)*X(90)-JVS(444)*X(91)-JVS(445)*X(92)-JVS(446)*X(93)-JVS(447)*X(94)&
            &-JVS(448)*X(95)-JVS(449)*X(96)-JVS(450)*X(97)-JVS(451)*X(98)-JVS(452)*X(99)-JVS(453)*X(100)-JVS(454)*X(101))&
            &/(JVS(433))
  X(73) = (X(73)-JVS(421)*X(85)-JVS(422)*X(89)-JVS(423)*X(94)-JVS(424)*X(95))/(JVS(420))
  X(72) = (X(72)-JVS(416)*X(85)-JVS(417)*X(89)-JVS(418)*X(94)-JVS(419)*X(95))/(JVS(415))
  X(71) = (X(71)-JVS(406)*X(90)-JVS(407)*X(91)-JVS(408)*X(92)-JVS(409)*X(93)-JVS(410)*X(94)-JVS(411)*X(95)-JVS(412)&
            &*X(96)-JVS(413)*X(97)-JVS(414)*X(100))/(JVS(405))
  X(70) = (X(70)-JVS(399)*X(85)-JVS(400)*X(89)-JVS(401)*X(94)-JVS(402)*X(95))/(JVS(398))
  X(69) = (X(69)-JVS(389)*X(72)-JVS(390)*X(76)-JVS(391)*X(77)-JVS(392)*X(80)-JVS(393)*X(84)-JVS(394)*X(89)-JVS(395)&
            &*X(93)-JVS(396)*X(94)-JVS(397)*X(95))/(JVS(388))
  X(68) = (X(68)-JVS(379)*X(85)-JVS(380)*X(89)-JVS(381)*X(94)-JVS(382)*X(95))/(JVS(378))
  X(67) = (X(67)-JVS(364)*X(71)-JVS(365)*X(75)-JVS(366)*X(78)-JVS(367)*X(79)-JVS(368)*X(80)-JVS(369)*X(81)-JVS(370)&
            &*X(82)-JVS(371)*X(83)-JVS(372)*X(86)-JVS(373)*X(89)-JVS(374)*X(93)-JVS(375)*X(94)-JVS(376)*X(95)-JVS(377)&
            &*X(97))/(JVS(363))
  X(66) = (X(66)-JVS(338)*X(68)-JVS(339)*X(70)-JVS(340)*X(72)-JVS(341)*X(73)-JVS(342)*X(75)-JVS(343)*X(76)-JVS(344)&
            &*X(77)-JVS(345)*X(78)-JVS(346)*X(79)-JVS(347)*X(80)-JVS(348)*X(81)-JVS(349)*X(82)-JVS(350)*X(83)-JVS(351)*X(85)&
            &-JVS(352)*X(86)-JVS(353)*X(89)-JVS(354)*X(94)-JVS(355)*X(95))/(JVS(337))
  X(65) = (X(65)-JVS(325)*X(75)-JVS(326)*X(79)-JVS(327)*X(81)-JVS(328)*X(89)-JVS(329)*X(94)-JVS(330)*X(95))/(JVS(324))
  X(64) = (X(64)-JVS(315)*X(71)-JVS(316)*X(94)-JVS(317)*X(95)-JVS(318)*X(97))/(JVS(314))
  X(63) = (X(63)-JVS(309)*X(80)-JVS(310)*X(89)-JVS(311)*X(94)-JVS(312)*X(95))/(JVS(308))
  X(62) = (X(62)-JVS(302)*X(95)-JVS(303)*X(97)-JVS(304)*X(99)-JVS(305)*X(101))/(JVS(301))
  X(61) = (X(61)-JVS(297)*X(71)-JVS(298)*X(93)-JVS(299)*X(94)-JVS(300)*X(97))/(JVS(296))
  X(60) = (X(60)-JVS(293)*X(89)-JVS(294)*X(95))/(JVS(292))
  X(59) = (X(59)-JVS(290)*X(89)-JVS(291)*X(95))/(JVS(289))
  X(58) = (X(58)-JVS(285)*X(94)-JVS(286)*X(95))/(JVS(284))
  X(57) = (X(57)-JVS(281)*X(95))/(JVS(280))
  X(56) = (X(56)-JVS(277)*X(95))/(JVS(276))
  X(55) = (X(55)-JVS(273)*X(95))/(JVS(272))
  X(54) = (X(54)-JVS(271)*X(95))/(JVS(270))
  X(53) = (X(53)-JVS(266)*X(95)-JVS(267)*X(98)-JVS(268)*X(99)-JVS(269)*X(101))/(JVS(265))
  X(52) = (X(52)-JVS(264)*X(95))/(JVS(263))
  X(51) = (X(51)-JVS(260)*X(93)-JVS(261)*X(95)-JVS(262)*X(97))/(JVS(259))
  X(50) = (X(50)-JVS(255)*X(61)-JVS(256)*X(93)-JVS(257)*X(94)-JVS(258)*X(97))/(JVS(254))
  X(49) = (X(49)-JVS(251)*X(83)-JVS(252)*X(96)-JVS(253)*X(97))/(JVS(250))
  X(48) = (X(48)-JVS(247)*X(95)-JVS(248)*X(97)-JVS(249)*X(98))/(JVS(246))
  X(47) = (X(47)-JVS(245)*X(95))/(JVS(244))
  X(46) = (X(46)-JVS(242)*X(93)-JVS(243)*X(95))/(JVS(241))
  X(45) = (X(45)-JVS(239)*X(95))/(JVS(238))
  X(44) = (X(44)-JVS(236)*X(95)-JVS(237)*X(96))/(JVS(235))
  X(43) = (X(43)-JVS(233)*X(94)-JVS(234)*X(95))/(JVS(232))
  X(42) = (X(42)-JVS(229)*X(43)-JVS(230)*X(94)-JVS(231)*X(95))/(JVS(228))
  X(41) = (X(41)-JVS(226)*X(93)-JVS(227)*X(94))/(JVS(225))
  X(40) = (X(40)-JVS(224)*X(95))/(JVS(223))
  X(39) = (X(39)-JVS(218)*X(52)-JVS(219)*X(76)-JVS(220)*X(77)-JVS(221)*X(89)-JVS(222)*X(95))/(JVS(217))
  X(38) = (X(38)-JVS(216)*X(95))/(JVS(215))
  X(37) = (X(37)-JVS(213)*X(95)-JVS(214)*X(97))/(JVS(212))
  X(36) = (X(36)-JVS(210)*X(93)-JVS(211)*X(100))/(JVS(209))
  X(35) = (X(35)-JVS(207)*X(91)-JVS(208)*X(93))/(JVS(206))
  X(34) = (X(34)-JVS(204)*X(92)-JVS(205)*X(93))/(JVS(203))
  X(33) = (X(33)-JVS(201)*X(90)-JVS(202)*X(93))/(JVS(200))
  X(32) = (X(32)-JVS(199)*X(95))/(JVS(198))
  X(31) = (X(31)-JVS(197)*X(95))/(JVS(196))
  X(30) = (X(30)-JVS(195)*X(95))/(JVS(194))
  X(29) = (X(29)-JVS(192)*X(95))/(JVS(191))
  X(28) = (X(28)-JVS(190)*X(95))/(JVS(189))
  X(27) = (X(27)-JVS(187)*X(95))/(JVS(186))
  X(26) = (X(26)-JVS(185)*X(89))/(JVS(184))
  X(25) = (X(25)-JVS(183)*X(95))/(JVS(182))
  X(24) = (X(24)-JVS(179)*X(25)-JVS(180)*X(29)-JVS(181)*X(95))/(JVS(178))
  X(23) = (X(23)-JVS(177)*X(95))/(JVS(176))
  X(22) = (X(22)-JVS(173)*X(23)-JVS(174)*X(27)-JVS(175)*X(95))/(JVS(172))
  X(21) = (X(21)-JVS(163)*X(22)-JVS(164)*X(23)-JVS(165)*X(24)-JVS(166)*X(25)-JVS(167)*X(27)-JVS(168)*X(28)-JVS(169)&
            &*X(29)-JVS(170)*X(30)-JVS(171)*X(95))/(JVS(162))
  X(20) = (X(20)-JVS(110)*X(31)-JVS(111)*X(32)-JVS(112)*X(38)-JVS(113)*X(40)-JVS(114)*X(42)-JVS(115)*X(44)-JVS(116)&
            &*X(45)-JVS(117)*X(47)-JVS(118)*X(48)-JVS(119)*X(52)-JVS(120)*X(53)-JVS(121)*X(54)-JVS(122)*X(55)-JVS(123)*X(56)&
            &-JVS(124)*X(57)-JVS(125)*X(58)-JVS(126)*X(59)-JVS(127)*X(60)-JVS(128)*X(62)-JVS(129)*X(63)-JVS(130)*X(64)&
            &-JVS(131)*X(65)-JVS(132)*X(66)-JVS(133)*X(67)-JVS(134)*X(68)-JVS(135)*X(69)-JVS(136)*X(70)-JVS(137)*X(72)&
            &-JVS(138)*X(73)-JVS(139)*X(75)-JVS(140)*X(76)-JVS(141)*X(77)-JVS(142)*X(78)-JVS(143)*X(79)-JVS(144)*X(80)&
            &-JVS(145)*X(81)-JVS(146)*X(82)-JVS(147)*X(83)-JVS(148)*X(84)-JVS(149)*X(86)-JVS(150)*X(87)-JVS(151)*X(88)&
            &-JVS(152)*X(89)-JVS(153)*X(93)-JVS(154)*X(94)-JVS(155)*X(95)-JVS(156)*X(96)-JVS(157)*X(97))/(JVS(109))
  X(19) = (X(19)-JVS(103)*X(73)-JVS(104)*X(76)-JVS(105)*X(77)-JVS(106)*X(89)-JVS(107)*X(94)-JVS(108)*X(95))/(JVS(102))
  X(18) = (X(18)-JVS(100)*X(77)-JVS(101)*X(89))/(JVS(99))
  X(17) = (X(17)-JVS(97)*X(77)-JVS(98)*X(95))/(JVS(96))
  X(16) = (X(16)-JVS(89)*X(47)-JVS(90)*X(52)-JVS(91)*X(54)-JVS(92)*X(55)-JVS(93)*X(72)-JVS(94)*X(80)-JVS(95)*X(95))&
            &/(JVS(88))
  X(15) = (X(15)-JVS(83)*X(74)-JVS(84)*X(97)-JVS(85)*X(98)-JVS(86)*X(99)-JVS(87)*X(101))/(JVS(82))
  X(14) = (X(14)-JVS(76)*X(74)-JVS(77)*X(94)-JVS(78)*X(96)-JVS(79)*X(98)-JVS(80)*X(99)-JVS(81)*X(101))/(JVS(75))
  X(13) = (X(13)-JVS(67)*X(50)-JVS(68)*X(63)-JVS(69)*X(70)-JVS(70)*X(85)-JVS(71)*X(89)-JVS(72)*X(93)-JVS(73)*X(94)&
            &-JVS(74)*X(95))/(JVS(66))
  X(12) = (X(12)-JVS(62)*X(50)-JVS(63)*X(70)-JVS(64)*X(93)-JVS(65)*X(94))/(JVS(61))
  X(11) = (X(11)-JVS(57)*X(91)-JVS(58)*X(92)-JVS(59)*X(97)-JVS(60)*X(100))/(JVS(56))
  X(10) = (X(10)-JVS(54)*X(90)-JVS(55)*X(97))/(JVS(53))
  X(9) = (X(9)-JVS(39)*X(72)-JVS(40)*X(73)-JVS(41)*X(76)-JVS(42)*X(77)-JVS(43)*X(79)-JVS(44)*X(80)-JVS(45)*X(89)-JVS(46)&
           &*X(91)-JVS(47)*X(92)-JVS(48)*X(97)-JVS(49)*X(98)-JVS(50)*X(99)-JVS(51)*X(100)-JVS(52)*X(101))/(JVS(38))
  X(8) = (X(8)-JVS(29)*X(70)-JVS(30)*X(72)-JVS(31)*X(80)-JVS(32)*X(89)-JVS(33)*X(90)-JVS(34)*X(97)-JVS(35)*X(98)-JVS(36)&
           &*X(99)-JVS(37)*X(101))/(JVS(28))
  X(7) = (X(7)-JVS(13)*X(49)-JVS(14)*X(60)-JVS(15)*X(68)-JVS(16)*X(70)-JVS(17)*X(72)-JVS(18)*X(73)-JVS(19)*X(75)-JVS(20)&
           &*X(76)-JVS(21)*X(77)-JVS(22)*X(79)-JVS(23)*X(80)-JVS(24)*X(81)-JVS(25)*X(89)-JVS(26)*X(95)-JVS(27)*X(96))&
           &/(JVS(12))
  X(6) = X(6)/JVS(11)
  X(5) = X(5)/JVS(10)
  X(4) = X(4)/JVS(9)
  X(3) = X(3)/JVS(8)
  X(2) = (X(2)-JVS(5)*X(60)-JVS(6)*X(70)-JVS(7)*X(89))/(JVS(4))
  X(1) = (X(1)-JVS(2)*X(42)-JVS(3)*X(95))/(JVS(1))
      
END SUBROUTINE saprc99_mosaic_20bin_vbs2_aq_KppSolve
























      SUBROUTINE saprc99_mosaic_20bin_vbs2_aq_WCOPY(N,X,incX,Y,incY)








      
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

      END SUBROUTINE saprc99_mosaic_20bin_vbs2_aq_WCOPY



      SUBROUTINE saprc99_mosaic_20bin_vbs2_aq_WAXPY(N,Alpha,X,incX,Y,incY)









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
      
      END SUBROUTINE saprc99_mosaic_20bin_vbs2_aq_WAXPY




      SUBROUTINE saprc99_mosaic_20bin_vbs2_aq_WSCAL(N,Alpha,X,incX)









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

      END SUBROUTINE saprc99_mosaic_20bin_vbs2_aq_WSCAL


      REAL(kind=dp) FUNCTION saprc99_mosaic_20bin_vbs2_aq_WLAMCH( C )








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
          CALL saprc99_mosaic_20bin_vbs2_aq_WLAMCH_ADD(ONE,Eps,Sum)
          IF (Sum.LE.ONE) GOTO 10
        END DO
        PRINT*,'ERROR IN WLAMCH. EPS < ',Eps
        RETURN
10      Eps = Eps*2
        i = i-1      
      END IF

      saprc99_mosaic_20bin_vbs2_aq_WLAMCH = Eps

      END FUNCTION saprc99_mosaic_20bin_vbs2_aq_WLAMCH
     
      SUBROUTINE saprc99_mosaic_20bin_vbs2_aq_WLAMCH_ADD( A, B, Sum )

      
      REAL(kind=dp) A, B, Sum
      Sum = A + B

      END SUBROUTINE saprc99_mosaic_20bin_vbs2_aq_WLAMCH_ADD




      SUBROUTINE saprc99_mosaic_20bin_vbs2_aq_SET2ZERO(N,Y)




      
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

      END SUBROUTINE saprc99_mosaic_20bin_vbs2_aq_SET2ZERO



      REAL(kind=dp) FUNCTION saprc99_mosaic_20bin_vbs2_aq_WDOT (N, DX, incX, DY, incY) 









      IMPLICIT NONE
      INTEGER :: N, incX, incY
      REAL(kind=dp) :: DX(N), DY(N) 

      INTEGER :: i, IX, IY, M, MP1, NS
                                 
      saprc99_mosaic_20bin_vbs2_aq_WDOT = 0.0D0 
      IF (N .LE. 0) RETURN 
      IF (incX .EQ. incY) IF (incX-1) 5,20,60 



    5 IX = 1 
      IY = 1 
      IF (incX .LT. 0) IX = (-N+1)*incX + 1 
      IF (incY .LT. 0) IY = (-N+1)*incY + 1 
      DO i = 1,N 
        saprc99_mosaic_20bin_vbs2_aq_WDOT = saprc99_mosaic_20bin_vbs2_aq_WDOT + DX(IX)*DY(IY) 
        IX = IX + incX 
        IY = IY + incY 
      END DO 
      RETURN 





   20 M = MOD(N,5) 
      IF (M .EQ. 0) GO TO 40 
      DO i = 1,M 
         saprc99_mosaic_20bin_vbs2_aq_WDOT = saprc99_mosaic_20bin_vbs2_aq_WDOT + DX(i)*DY(i) 
      END DO 
      IF (N .LT. 5) RETURN 
   40 MP1 = M + 1 
      DO i = MP1,N,5 
          saprc99_mosaic_20bin_vbs2_aq_WDOT = saprc99_mosaic_20bin_vbs2_aq_WDOT + DX(i)*DY(i) + DX(i+1)*DY(i+1) +&
                   DX(i+2)*DY(i+2) +  &
                   DX(i+3)*DY(i+3) + DX(i+4)*DY(i+4)                   
      END DO 
      RETURN 



   60 NS = N*incX 
      DO i = 1,NS,incX 
        saprc99_mosaic_20bin_vbs2_aq_WDOT = saprc99_mosaic_20bin_vbs2_aq_WDOT + DX(i)*DY(i) 
      END DO 

      END FUNCTION saprc99_mosaic_20bin_vbs2_aq_WDOT                                          




END MODULE saprc99_mosaic_20bin_vbs2_aq_Integrator
