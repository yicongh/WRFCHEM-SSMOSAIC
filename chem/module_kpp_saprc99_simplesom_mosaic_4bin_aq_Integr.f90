
































MODULE saprc99_simplesom_mosaic_4bin_aq_Integrator

 USE saprc99_simplesom_mosaic_4bin_aq_Parameters
 USE saprc99_simplesom_mosaic_4bin_aq_Precision
 USE saprc99_simplesom_mosaic_4bin_aq_JacobianSP

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

SUBROUTINE  saprc99_simplesom_mosaic_4bin_aq_INTEGRATE( TIN, TOUT, &
  FIX, VAR,  RCONST, ATOL, RTOL, IRR_WRK,  &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U  )

   USE saprc99_simplesom_mosaic_4bin_aq_Parameters

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

   CALL saprc99_simplesom_mosaic_4bin_aq_Rosenbrock(VAR, FIX, RCONST, TIN,TOUT,   &
         ATOL,RTOL,               &
         RCNTRL,ICNTRL,RSTATUS,ISTATUS,IRR_WRK,IERR)

   STEPMIN = RCNTRL(ihexit)
   
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U))    IERR_U       = IERR

END SUBROUTINE  saprc99_simplesom_mosaic_4bin_aq_INTEGRATE


SUBROUTINE  saprc99_simplesom_mosaic_4bin_aq_Rosenbrock(Y, FIX, RCONST, Tstart,Tend, &
           AbsTol,RelTol,            &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IRR_WRK,IERR)







































































































  USE saprc99_simplesom_mosaic_4bin_aq_Parameters

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
      CALL saprc99_simplesom_mosaic_4bin_aq_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
      RETURN
   END IF


   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 100000
   ELSEIF (ICNTRL(4) > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL saprc99_simplesom_mosaic_4bin_aq_ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF


   Roundoff = saprc99_simplesom_mosaic_4bin_aq_WLAMCH('E')


   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      Hmin = RCNTRL(1)
   ELSE
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL saprc99_simplesom_mosaic_4bin_aq_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL saprc99_simplesom_mosaic_4bin_aq_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL saprc99_simplesom_mosaic_4bin_aq_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2_dp
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL saprc99_simplesom_mosaic_4bin_aq_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0_dp
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL saprc99_simplesom_mosaic_4bin_aq_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1_dp
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL saprc99_simplesom_mosaic_4bin_aq_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9_dp
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL saprc99_simplesom_mosaic_4bin_aq_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.0_dp*Roundoff) &
         .OR. (RelTol(i) >= 1.0_dp) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL saprc99_simplesom_mosaic_4bin_aq_ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO



   SELECT CASE (Method)
     CASE (1)
       CALL saprc99_simplesom_mosaic_4bin_aq_Ros2(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (2)
       CALL saprc99_simplesom_mosaic_4bin_aq_Ros3(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (3)
       CALL saprc99_simplesom_mosaic_4bin_aq_Ros4(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (4)
       CALL saprc99_simplesom_mosaic_4bin_aq_Rodas3(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (5)
       CALL saprc99_simplesom_mosaic_4bin_aq_Rodas4(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(4)=', Method
       CALL saprc99_simplesom_mosaic_4bin_aq_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       RETURN
   END SELECT


   CALL saprc99_simplesom_mosaic_4bin_aq_ros_Integrator(Y,Tstart,Tend,Texit,      &
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



 SUBROUTINE  saprc99_simplesom_mosaic_4bin_aq_ros_ErrorMsg(Code,T,H,IERR)



   USE saprc99_simplesom_mosaic_4bin_aq_Precision

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

 END SUBROUTINE  saprc99_simplesom_mosaic_4bin_aq_ros_ErrorMsg


 SUBROUTINE  saprc99_simplesom_mosaic_4bin_aq_ros_Integrator (Y, Tstart, Tend, T,     &
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
      CALL saprc99_simplesom_mosaic_4bin_aq_ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1_dp*H) == T).OR.(H <= Roundoff) ) THEN  
      CALL saprc99_simplesom_mosaic_4bin_aq_ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF


   Hexit = H
   H = MIN(H,ABS(Tend-T))


   CALL saprc99_simplesom_mosaic_4bin_aq_FunTemplate(T,Y,Fcn0, RCONST, FIX, Nfun)
   IF( T == Tstart ) THEN
     CALL saprc99_simplesom_mosaic_4bin_aq_IRRFun( Y, FIX, RCONST, IRR_WRK )
   ENDIF


   IF (.NOT.Autonomous) THEN
      CALL saprc99_simplesom_mosaic_4bin_aq_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )
   END IF


   CALL saprc99_simplesom_mosaic_4bin_aq_JacTemplate(T,Y,Jac0, FIX, Njac, RCONST)


UntilAccepted: DO

   CALL saprc99_simplesom_mosaic_4bin_aq_ros_PrepareMatrix(H,Direction,ros_Gamma(1), &
          Jac0,Ghimj,Pivot,Singular, Ndec,  Nsng )
   IF (Singular) THEN 
       CALL saprc99_simplesom_mosaic_4bin_aq_ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF


Stage: DO istage = 1, ros_S

      
       ioffset = NVAR*(istage-1)

      
       IF ( istage == 1 ) THEN
         CALL saprc99_simplesom_mosaic_4bin_aq_WCOPY(NVAR,Fcn0,1,Fcn,1)
      
       ELSEIF ( ros_NewF(istage) ) THEN
         CALL saprc99_simplesom_mosaic_4bin_aq_WCOPY(NVAR,Y,1,Ynew,1)
         DO j = 1, istage-1
           CALL saprc99_simplesom_mosaic_4bin_aq_WAXPY(NVAR,ros_A((istage-1)*(istage-2)/2+j), &
            K(NVAR*(j-1)+1),1,Ynew,1)
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL saprc99_simplesom_mosaic_4bin_aq_FunTemplate(Tau,Ynew,Fcn, RCONST, FIX, Nfun)
       END IF 
       CALL saprc99_simplesom_mosaic_4bin_aq_WCOPY(NVAR,Fcn,1,K(ioffset+1),1)
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         CALL saprc99_simplesom_mosaic_4bin_aq_WAXPY(NVAR,HC,K(NVAR*(j-1)+1),1,K(ioffset+1),1)
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         CALL saprc99_simplesom_mosaic_4bin_aq_WAXPY(NVAR,HG,dFdT,1,K(ioffset+1),1)
       END IF
       CALL saprc99_simplesom_mosaic_4bin_aq_ros_Solve(Ghimj, Pivot, K(ioffset+1), Nsol)

   END DO Stage



   CALL saprc99_simplesom_mosaic_4bin_aq_WCOPY(NVAR,Y,1,Ynew,1)
   DO j=1,ros_S
         CALL saprc99_simplesom_mosaic_4bin_aq_WAXPY(NVAR,ros_M(j),K(NVAR*(j-1)+1),1,Ynew,1)
   END DO


   CALL saprc99_simplesom_mosaic_4bin_aq_WSCAL(NVAR,ZERO,Yerr,1)
   DO j=1,ros_S
        CALL saprc99_simplesom_mosaic_4bin_aq_WAXPY(NVAR,ros_E(j),K(NVAR*(j-1)+1),1,Yerr,1)
   END DO
   Err = saprc99_simplesom_mosaic_4bin_aq_ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )


   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac


   Nstp = Nstp+1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  
      Nacc = Nacc+1
      CALL saprc99_simplesom_mosaic_4bin_aq_WCOPY(NVAR,Ynew,1,Y,1)
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

  END SUBROUTINE  saprc99_simplesom_mosaic_4bin_aq_ros_Integrator



  REAL(kind=dp) FUNCTION  saprc99_simplesom_mosaic_4bin_aq_ros_ErrorNorm ( Y, Ynew, Yerr, &
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

    saprc99_simplesom_mosaic_4bin_aq_ros_ErrorNorm = Err

  END FUNCTION  saprc99_simplesom_mosaic_4bin_aq_ros_ErrorNorm



  SUBROUTINE saprc99_simplesom_mosaic_4bin_aq_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )



   IMPLICIT NONE


   REAL(kind=dp), INTENT(IN) :: T, Roundoff, Y(NVAR), Fcn0(NVAR)
   REAL(kind=dp), INTENT(IN) :: RCONST(NREACT), FIX(NFIX)

   REAL(kind=dp), INTENT(OUT) :: dFdT(NVAR)

   INTEGER, INTENT(INOUT) ::Nfun

   REAL(kind=dp) :: Delta
   REAL(kind=dp), PARAMETER :: ONE = 1.0_dp, DeltaMin = 1.0E-6_dp

   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL saprc99_simplesom_mosaic_4bin_aq_FunTemplate(T+Delta,Y,dFdT, RCONST, FIX, Nfun)
   CALL saprc99_simplesom_mosaic_4bin_aq_WAXPY(NVAR,(-ONE),Fcn0,1,dFdT,1)
   CALL saprc99_simplesom_mosaic_4bin_aq_WSCAL(NVAR,(ONE/Delta),dFdT,1)

  END SUBROUTINE  saprc99_simplesom_mosaic_4bin_aq_ros_FunTimeDeriv



  SUBROUTINE  saprc99_simplesom_mosaic_4bin_aq_ros_PrepareMatrix ( H, Direction, gam, &
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


     CALL saprc99_simplesom_mosaic_4bin_aq_WCOPY(LU_NONZERO,Jac0,1,Ghimj,1)
     CALL saprc99_simplesom_mosaic_4bin_aq_WSCAL(LU_NONZERO,(-ONE),Ghimj,1)
     ghinv = ONE/(Direction*H*gam)
     DO i=1,NVAR
       Ghimj(LU_DIAG(i)) = Ghimj(LU_DIAG(i))+ghinv
     END DO

     CALL saprc99_simplesom_mosaic_4bin_aq_ros_Decomp( Ghimj, Pivot, ising, Ndec )
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

  END SUBROUTINE  saprc99_simplesom_mosaic_4bin_aq_ros_PrepareMatrix



  SUBROUTINE  saprc99_simplesom_mosaic_4bin_aq_ros_Decomp( A, Pivot, ising, Ndec )



   IMPLICIT NONE

   REAL(kind=dp), INTENT(INOUT) :: A(LU_NONZERO)

   INTEGER, INTENT(OUT) :: Pivot(NVAR), ising
   INTEGER, INTENT(INOUT) :: Ndec 

   CALL saprc99_simplesom_mosaic_4bin_aq_KppDecomp ( A, ising )
   Pivot(1) = 1
   Ndec = Ndec + 1

  END SUBROUTINE  saprc99_simplesom_mosaic_4bin_aq_ros_Decomp



  SUBROUTINE  saprc99_simplesom_mosaic_4bin_aq_ros_Solve( A, Pivot, b, Nsol )



   IMPLICIT NONE

   REAL(kind=dp), INTENT(IN) :: A(LU_NONZERO)
   INTEGER, INTENT(IN) :: Pivot(NVAR)

   INTEGER, INTENT(INOUT) :: nsol 

   REAL(kind=dp), INTENT(INOUT) :: b(NVAR)


   CALL saprc99_simplesom_mosaic_4bin_aq_KppSolve( A, b )

   Nsol = Nsol+1

  END SUBROUTINE  saprc99_simplesom_mosaic_4bin_aq_ros_Solve




  SUBROUTINE  saprc99_simplesom_mosaic_4bin_aq_Ros2 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

 END SUBROUTINE  saprc99_simplesom_mosaic_4bin_aq_Ros2



  SUBROUTINE  saprc99_simplesom_mosaic_4bin_aq_Ros3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  saprc99_simplesom_mosaic_4bin_aq_Ros3





  SUBROUTINE  saprc99_simplesom_mosaic_4bin_aq_Ros4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  saprc99_simplesom_mosaic_4bin_aq_Ros4


  SUBROUTINE  saprc99_simplesom_mosaic_4bin_aq_Rodas3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  saprc99_simplesom_mosaic_4bin_aq_Rodas3


  SUBROUTINE  saprc99_simplesom_mosaic_4bin_aq_Rodas4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  saprc99_simplesom_mosaic_4bin_aq_Rodas4




END SUBROUTINE  saprc99_simplesom_mosaic_4bin_aq_Rosenbrock




SUBROUTINE  saprc99_simplesom_mosaic_4bin_aq_FunTemplate( T, Y, Ydot, RCONST, FIX, Nfun )




   USE saprc99_simplesom_mosaic_4bin_aq_Parameters




   REAL(kind=dp) :: T, Y(NVAR)
   REAL(kind=dp) :: RCONST(NREACT)
   REAL(kind=dp) :: FIX(NFIX)

   REAL(kind=dp) :: Ydot(NVAR)
   INTEGER :: Nfun









   CALL saprc99_simplesom_mosaic_4bin_aq_Fun( Y, FIX, RCONST, Ydot )


   Nfun = Nfun+1

END SUBROUTINE  saprc99_simplesom_mosaic_4bin_aq_FunTemplate



SUBROUTINE  saprc99_simplesom_mosaic_4bin_aq_JacTemplate( T, Y, Jcb, FIX, Njac, RCONST )




 USE saprc99_simplesom_mosaic_4bin_aq_Parameters
 
 USE saprc99_simplesom_mosaic_4bin_aq_Jacobian



    REAL(kind=dp) :: T, Y(NVAR)
    REAL(kind=dp) :: FIX(NFIX)
    REAL(kind=dp) :: RCONST(NREACT)

    INTEGER :: Njac


    REAL(kind=dp) :: Jcb(LU_NONZERO)

    REAL(kind=dp) :: Told





    CALL saprc99_simplesom_mosaic_4bin_aq_Jac_SP( Y, FIX, RCONST, Jcb )


    Njac = Njac+1

END SUBROUTINE  saprc99_simplesom_mosaic_4bin_aq_JacTemplate

















SUBROUTINE saprc99_simplesom_mosaic_4bin_aq_Fun ( V, F, RCT, Vdot )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: Vdot(NVAR)




  REAL(kind=dp) :: A(NREACT)


  A(1) = RCT(1)*V(94)
  A(2) = RCT(2)*V(88)*F(2)
  A(3) = RCT(3)*V(88)*V(92)
  A(4) = RCT(4)*V(88)*V(93)*F(2)
  A(5) = RCT(5)*V(88)*V(94)
  A(6) = RCT(6)*V(88)*V(94)
  A(7) = RCT(7)*V(92)*V(93)
  A(8) = RCT(8)*V(92)*V(94)
  A(9) = RCT(9)*V(93)*V(95)
  A(10) = RCT(10)*V(93)*V(93)*F(2)
  A(11) = RCT(11)*V(94)*V(95)
  A(12) = RCT(12)*V(31)
  A(13) = RCT(13)*V(31)*F(1)
  A(14) = RCT(14)*V(94)*V(95)
  A(15) = RCT(15)*V(95)
  A(16) = RCT(16)*V(95)
  A(17) = RCT(17)*V(92)
  A(18) = RCT(18)*V(92)
  A(19) = RCT(19)*V(18)*F(1)
  A(20) = RCT(20)*V(18)*F(2)
  A(21) = RCT(21)*V(93)*V(102)
  A(22) = RCT(22)*V(33)
  A(23) = RCT(23)*V(33)
  A(24) = RCT(24)*V(33)*V(102)
  A(25) = RCT(25)*V(94)*V(102)
  A(26) = RCT(26)*V(95)*V(102)
  A(27) = RCT(27)*V(70)*V(102)
  A(28) = RCT(28)*V(70)
  A(29) = RCT(29)*V(69)*V(102)
  A(30) = RCT(30)*V(92)*V(102)
  A(31) = RCT(31)*V(93)*V(103)
  A(32) = RCT(32)*V(94)*V(103)
  A(33) = RCT(33)*V(43)
  A(34) = RCT(34)*V(43)
  A(35) = RCT(35)*V(43)*V(102)
  A(36) = RCT(36)*V(92)*V(103)
  A(37) = RCT(37)*V(103)*V(103)
  A(38) = RCT(38)*V(103)*V(103)*F(1)
  A(39) = RCT(39)*V(95)*V(103)
  A(40) = RCT(40)*V(95)*V(95)
  A(41) = RCT(41)*V(27)
  A(42) = RCT(42)*V(27)*V(102)
  A(43) = RCT(43)*V(102)*V(103)
  A(44) = RCT(44)*V(32)*V(102)
  A(45) = RCT(45)*V(102)*F(2)
  A(46) = RCT(46)*V(93)*V(99)
  A(47) = RCT(47)*V(99)*V(103)
  A(48) = RCT(48)*V(95)*V(99)
  A(49) = RCT(49)*V(99)*V(99)
  A(50) = RCT(50)*V(99)*V(99)
  A(51) = RCT(51)*V(93)*V(104)
  A(52) = RCT(52)*V(103)*V(104)
  A(53) = RCT(53)*V(95)*V(104)
  A(54) = RCT(54)*V(99)*V(104)
  A(55) = RCT(55)*V(104)*V(104)
  A(56) = RCT(56)*V(78)*V(93)
  A(57) = RCT(57)*V(78)*V(103)
  A(58) = RCT(58)*V(78)*V(95)
  A(59) = RCT(59)*V(78)*V(99)
  A(60) = RCT(60)*V(78)*V(104)
  A(61) = RCT(61)*V(78)*V(78)
  A(62) = RCT(62)*V(93)*V(101)
  A(63) = RCT(63)*V(101)*V(103)
  A(64) = RCT(64)*V(99)*V(101)
  A(65) = RCT(65)*V(95)*V(101)
  A(66) = RCT(66)*V(101)*V(104)
  A(67) = RCT(67)*V(78)*V(101)
  A(68) = RCT(68)*V(101)*V(101)
  A(69) = RCT(69)*V(94)*V(100)
  A(70) = RCT(70)*V(23)
  A(71) = RCT(71)*V(93)*V(100)
  A(72) = RCT(72)*V(100)*V(103)
  A(73) = RCT(73)*V(95)*V(100)
  A(74) = RCT(74)*V(99)*V(100)
  A(75) = RCT(75)*V(100)*V(104)
  A(76) = RCT(76)*V(78)*V(100)
  A(77) = RCT(77)*V(100)*V(101)
  A(78) = RCT(78)*V(100)*V(100)
  A(79) = RCT(79)*V(94)*V(97)
  A(80) = RCT(80)*V(24)
  A(81) = RCT(81)*V(93)*V(97)
  A(82) = RCT(82)*V(97)*V(103)
  A(83) = RCT(83)*V(95)*V(97)
  A(84) = RCT(84)*V(97)*V(99)
  A(85) = RCT(85)*V(97)*V(104)
  A(86) = RCT(86)*V(78)*V(97)
  A(87) = RCT(87)*V(97)*V(101)
  A(88) = RCT(88)*V(97)*V(100)
  A(89) = RCT(89)*V(97)*V(97)
  A(90) = RCT(90)*V(94)*V(96)
  A(91) = RCT(91)*V(25)
  A(92) = RCT(92)*V(93)*V(96)
  A(93) = RCT(93)*V(96)*V(103)
  A(94) = RCT(94)*V(95)*V(96)
  A(95) = RCT(95)*V(96)*V(99)
  A(96) = RCT(96)*V(96)*V(104)
  A(97) = RCT(97)*V(78)*V(96)
  A(98) = RCT(98)*V(96)*V(101)
  A(99) = RCT(99)*V(96)*V(100)
  A(100) = RCT(100)*V(96)*V(97)
  A(101) = RCT(101)*V(96)*V(96)
  A(102) = RCT(102)*V(94)*V(98)
  A(103) = RCT(103)*V(26)
  A(104) = RCT(104)*V(93)*V(98)
  A(105) = RCT(105)*V(98)*V(103)
  A(106) = RCT(106)*V(95)*V(98)
  A(107) = RCT(107)*V(98)*V(99)
  A(108) = RCT(108)*V(98)*V(104)
  A(109) = RCT(109)*V(78)*V(98)
  A(110) = RCT(110)*V(98)*V(101)
  A(111) = RCT(111)*V(98)*V(100)
  A(112) = RCT(112)*V(97)*V(98)
  A(113) = RCT(113)*V(96)*V(98)
  A(114) = RCT(114)*V(98)*V(98)
  A(115) = RCT(115)*V(36)*V(94)
  A(116) = RCT(116)*V(36)
  A(117) = RCT(117)*V(76)*V(94)
  A(118) = RCT(118)*V(76)*V(103)
  A(119) = RCT(119)*V(76)
  A(120) = RCT(120)*V(42)*V(94)
  A(121) = RCT(121)*V(42)*V(103)
  A(122) = RCT(122)*V(42)
  A(123) = RCT(123)*V(86)
  A(124) = RCT(124)*V(86)
  A(125) = RCT(125)*V(86)*V(102)
  A(126) = RCT(126)*V(86)*V(103)
  A(127) = RCT(127)*V(41)
  A(128) = RCT(128)*V(41)*V(93)
  A(129) = RCT(129)*V(86)*V(95)
  A(130) = RCT(130)*V(85)*V(102)
  A(131) = RCT(131)*V(85)
  A(132) = RCT(132)*V(85)*V(95)
  A(133) = RCT(133)*V(89)*V(102)
  A(134) = RCT(134)*V(89)
  A(135) = RCT(135)*V(89)*V(95)
  A(136) = RCT(136)*V(72)*V(102)
  A(137) = RCT(137)*V(72)
  A(138) = RCT(138)*V(90)*V(102)
  A(139) = RCT(139)*V(90)
  A(140) = RCT(140)*V(39)*V(102)
  A(141) = RCT(141)*V(30)*V(102)
  A(142) = RCT(142)*V(40)*V(102)
  A(143) = RCT(143)*V(40)
  A(144) = RCT(144)*V(53)*V(102)
  A(145) = RCT(145)*V(53)
  A(146) = RCT(146)*V(75)
  A(147) = RCT(147)*V(75)
  A(148) = RCT(148)*V(75)*V(102)
  A(149) = RCT(149)*V(75)*V(95)
  A(150) = RCT(150)*V(68)
  A(151) = RCT(151)*V(68)*V(102)
  A(152) = RCT(152)*V(68)*V(95)
  A(153) = RCT(153)*V(28)
  A(154) = RCT(154)*V(55)*V(102)
  A(155) = RCT(155)*V(55)*V(95)
  A(156) = RCT(156)*V(48)*V(102)
  A(157) = RCT(157)*V(48)*V(95)
  A(158) = RCT(158)*V(51)*V(95)
  A(159) = RCT(159)*V(52)*V(102)
  A(160) = RCT(160)*V(52)
  A(161) = RCT(161)*V(52)*V(95)
  A(162) = RCT(162)*V(81)*V(102)
  A(163) = RCT(163)*V(81)*V(92)
  A(164) = RCT(164)*V(81)*V(95)
  A(165) = RCT(165)*V(81)*V(88)
  A(166) = RCT(166)*V(81)
  A(167) = RCT(167)*V(84)*V(102)
  A(168) = RCT(168)*V(84)*V(92)
  A(169) = RCT(169)*V(84)*V(88)
  A(170) = RCT(170)*V(84)
  A(171) = RCT(171)*V(82)*V(102)
  A(172) = RCT(172)*V(82)*V(92)
  A(173) = RCT(173)*V(82)*V(95)
  A(174) = RCT(174)*V(82)
  A(175) = RCT(175)*V(91)*V(102)
  A(176) = RCT(176)*V(91)
  A(177) = RCT(177)*V(87)*V(102)
  A(178) = RCT(178)*V(87)
  A(179) = RCT(179)*V(49)*V(102)
  A(180) = RCT(180)*V(49)*V(92)
  A(181) = RCT(181)*V(46)*V(102)
  A(182) = RCT(182)*V(46)
  A(183) = RCT(183)*V(47)*V(102)
  A(184) = RCT(184)*V(47)
  A(185) = RCT(185)*V(19)*V(102)
  A(186) = RCT(186)*V(71)*V(102)
  A(187) = RCT(187)*V(71)*V(92)
  A(188) = RCT(188)*V(71)*V(95)
  A(189) = RCT(189)*V(71)*V(88)
  A(190) = RCT(190)*V(80)*V(102)
  A(191) = RCT(191)*V(21)*V(103)
  A(192) = RCT(192)*V(21)*V(93)
  A(193) = RCT(193)*V(21)*V(21)
  A(194) = RCT(194)*V(20)*V(102)
  A(195) = RCT(195)*V(16)*V(102)
  A(196) = RCT(196)*V(17)*V(102)
  A(197) = RCT(197)*V(15)*V(103)
  A(198) = RCT(198)*V(15)*V(93)
  A(199) = RCT(199)*V(15)
  A(200) = RCT(200)*V(80)*V(92)
  A(201) = RCT(201)*V(80)*V(95)
  A(202) = RCT(202)*V(80)*V(88)
  A(203) = RCT(203)*V(74)*V(92)
  A(204) = RCT(204)*V(74)*V(95)
  A(205) = RCT(205)*V(74)*V(88)
  A(206) = RCT(206)*V(79)*V(102)
  A(207) = RCT(207)*V(79)*V(92)
  A(208) = RCT(208)*V(79)*V(95)
  A(209) = RCT(209)*V(79)*V(88)
  A(210) = RCT(210)*V(22)*V(102)
  A(211) = RCT(211)*V(29)*V(102)
  A(212) = RCT(212)*V(50)*V(102)
  A(213) = RCT(213)*V(35)*V(102)
  A(214) = RCT(214)*V(44)*V(102)
  A(215) = RCT(215)*V(37)*V(102)
  A(216) = RCT(216)*V(45)*V(102)
  A(217) = RCT(217)*V(38)*V(102)
  A(218) = RCT(218)*V(77)*V(102)
  A(219) = RCT(219)*V(77)*V(92)
  A(220) = RCT(220)*V(77)*V(95)
  A(221) = RCT(221)*V(77)*V(88)
  A(222) = RCT(222)*V(83)*V(102)
  A(223) = RCT(223)*V(83)*V(92)
  A(224) = RCT(224)*V(83)*V(95)
  A(225) = RCT(225)*V(83)*V(88)
  A(226) = RCT(226)*V(50)*V(92)
  A(227) = RCT(227)*V(73)*V(102)
  A(228) = RCT(228)*V(73)*V(92)
  A(229) = RCT(229)*V(73)*V(95)
  A(230) = RCT(230)*V(73)*V(88)
  A(231) = RCT(231)*V(34)*V(102)
  A(232) = RCT(232)*V(34)*V(102)
  A(233) = RCT(233)*V(34)*V(95)
  A(234) = RCT(234)*V(64)*V(102)
  A(235) = RCT(235)*V(67)*V(102)
  A(236) = RCT(236)*V(66)*V(102)
  A(237) = RCT(237)*V(65)*V(102)
  A(238) = RCT(238)*V(63)*V(102)
  A(239) = RCT(239)*V(62)*V(102)
  A(240) = RCT(240)*V(61)*V(102)
  A(241) = RCT(241)*V(60)*V(102)
  A(242) = RCT(242)*V(59)*V(102)
  A(243) = RCT(243)*V(58)*V(102)
  A(244) = RCT(244)*V(56)*V(102)
  A(245) = RCT(245)*V(54)*V(102)
  A(246) = RCT(246)*V(57)*V(102)
  A(247) = RCT(247)*V(74)*V(102)


  Vdot(1) = A(44)
  Vdot(2) = A(128)+0.333*A(163)+0.351*A(168)+0.1*A(172)+0.37*A(187)+0.204*A(200)+0.103*A(203)+0.103*A(207)+0.297*A(212)&
              &+0.185*A(219)+0.073*A(223)+0.185*A(228)
  Vdot(3) = 0.25*A(72)+A(74)+A(75)+A(77)+0.05*A(219)+0.129*A(223)+0.17*A(228)
  Vdot(4) = 0.25*A(82)+A(84)+A(85)+A(87)+0.25*A(93)+A(95)+A(96)+A(98)+0.25*A(105)+A(107)+A(108)+2*A(110)+0.372*A(172)&
              &+0.15*A(200)+0.189*A(203)+0.189*A(207)+0.119*A(219)+0.247*A(223)
  Vdot(5) = 0.5*A(226)+0.135*A(228)
  Vdot(6) = 0.75*A(72)
  Vdot(7) = 0.75*A(82)+0.75*A(93)+0.75*A(105)
  Vdot(8) = 2*A(120)+A(229)
  Vdot(9) = 6*A(120)+7*A(160)+0.048*A(227)+0.07*A(228)+2.693*A(229)+0.55*A(230)
  Vdot(10) = A(46)+A(48)+A(51)+A(53)+A(56)+A(58)+A(62)+A(65)
  Vdot(11) = A(47)+A(49)+A(50)+A(52)+A(54)+A(55)+A(57)+A(59)+A(60)+A(61)+A(63)+A(64)+A(66)+A(67)+A(68)
  Vdot(12) = A(197)
  Vdot(13) = A(199)
  Vdot(14) = A(193)+A(195)+A(196)+A(198)
  Vdot(15) = 0.12*A(194)-A(197)-A(198)-A(199)
  Vdot(16) = 0.75*A(194)-A(195)
  Vdot(17) = -A(196)
  Vdot(18) = A(18)-A(19)-A(20)
  Vdot(19) = -A(185)
  Vdot(20) = 0.88*A(191)-A(194)
  Vdot(21) = A(190)-A(191)-A(192)-2*A(193)+0.13*A(194)
  Vdot(22) = -A(210)
  Vdot(23) = A(69)-A(70)
  Vdot(24) = A(79)-A(80)
  Vdot(25) = A(90)-A(91)
  Vdot(26) = A(102)-A(103)
  Vdot(27) = A(37)+A(38)-A(41)-A(42)
  Vdot(28) = -A(153)+0.031*A(203)+0.031*A(207)+0.087*A(217)
  Vdot(29) = -A(211)
  Vdot(30) = -A(141)
  Vdot(31) = A(11)-A(12)-A(13)
  Vdot(32) = -A(44)+A(231)+0.5*A(232)+A(233)
  Vdot(33) = A(21)-A(22)-A(23)-A(24)
  Vdot(34) = -A(231)-A(232)-A(233)
  Vdot(35) = -A(213)
  Vdot(36) = -A(115)-A(116)+0.236*A(213)
  Vdot(37) = -A(215)
  Vdot(38) = -A(217)
  Vdot(39) = A(49)+0.25*A(54)+0.25*A(64)-A(140)
  Vdot(40) = A(47)-A(142)-A(143)
  Vdot(41) = A(126)-A(127)-A(128)
  Vdot(42) = -A(120)-A(121)-A(122)+A(158)
  Vdot(43) = A(32)-A(33)-A(34)-A(35)
  Vdot(44) = -A(214)
  Vdot(45) = -A(216)
  Vdot(46) = -A(181)-A(182)+0.108*A(216)+0.099*A(217)
  Vdot(47) = -A(183)-A(184)+0.051*A(216)+0.093*A(217)
  Vdot(48) = -A(156)-A(157)+0.207*A(216)+0.187*A(217)
  Vdot(49) = -A(179)-A(180)+0.491*A(216)+0.561*A(217)
  Vdot(50) = -A(212)-A(226)
  Vdot(51) = A(117)+A(121)+A(122)-A(158)
  Vdot(52) = -A(159)-A(160)-A(161)+0.059*A(216)+0.05*A(217)+0.061*A(222)+0.042*A(223)+0.015*A(224)
  Vdot(53) = A(52)+A(63)-A(144)-A(145)
  Vdot(54) = 0.006*A(239)+0.022*A(240)+0.072*A(241)+0.098*A(242)+0.021*A(243)-A(245)+0.005*A(246)
  Vdot(55) = A(118)+A(119)-A(154)-A(155)+0.017*A(216)
  Vdot(56) = 0.008*A(238)+0.027*A(239)+0.088*A(240)+0.119*A(241)+0.026*A(242)-A(244)+0.004*A(245)+0.005*A(246)
  Vdot(57) = 0.012*A(234)+0.01*A(235)+0.008*A(236)+0.007*A(237)+0.005*A(238)+0.005*A(239)+0.009*A(240)+0.025*A(241)+0.08&
               &*A(242)+0.145*A(243)+0.136*A(244)+0.112*A(245)-0.908*A(246)+0.019*A(247)
  Vdot(58) = 0.009*A(237)+0.033*A(238)+0.107*A(239)+0.146*A(240)+0.031*A(241)-A(243)+0.004*A(244)+0.004*A(245)+0.006&
               &*A(247)
  Vdot(59) = 0.011*A(236)+0.04*A(237)+0.13*A(238)+0.178*A(239)+0.038*A(240)-A(242)+0.004*A(243)+0.004*A(244)+0.028&
               &*A(247)
  Vdot(60) = 0.014*A(235)+0.049*A(236)+0.159*A(237)+0.217*A(238)+0.047*A(239)-A(241)+0.004*A(242)+0.004*A(243)+0.114&
               &*A(247)
  Vdot(61) = 0.017*A(234)+0.06*A(235)+0.194*A(236)+0.265*A(237)+0.057*A(238)-A(240)+0.004*A(241)+0.004*A(242)+0.234&
               &*A(247)
  Vdot(62) = 0.073*A(234)+0.237*A(235)+0.323*A(236)+0.07*A(237)-A(239)+0.003*A(240)+0.004*A(241)+0.169*A(247)
  Vdot(63) = 0.289*A(234)+0.394*A(235)+0.085*A(236)-A(238)+0.003*A(239)+0.003*A(240)+0.034*A(247)
  Vdot(64) = -A(234)+0.001*A(235)+0.002*A(236)+0.001*A(247)
  Vdot(65) = 0.481*A(234)+0.104*A(235)-A(237)+0.003*A(238)+0.003*A(239)
  Vdot(66) = 0.127*A(234)-A(236)+0.002*A(237)+0.003*A(238)+0.001*A(247)
  Vdot(67) = -A(235)+0.002*A(236)+0.002*A(237)+0.002*A(247)
  Vdot(68) = -A(150)-A(151)-A(152)+0.23*A(156)+0.084*A(162)+0.9*A(163)+0.3*A(167)+0.95*A(168)+0.174*A(171)+0.742*A(172)&
               &+0.008*A(173)+0.5*A(182)+0.5*A(184)+0.119*A(216)+0.287*A(217)
  Vdot(69) = -A(29)+A(123)+A(124)+A(125)+A(129)+A(131)+0.034*A(133)+A(134)+2*A(146)+A(147)+1.26*A(148)+1.26*A(149)&
               &+A(150)+A(151)+A(152)+0.416*A(162)+0.45*A(163)+0.5*A(164)+0.67*A(166)+0.475*A(168)+0.7*A(170)+0.336*A(171)&
               &+0.498*A(172)+0.572*A(173)+1.233*A(174)+A(179)+1.5*A(180)+A(182)+A(184)+0.5*A(187)+0.491*A(189)+0.275*A(200)&
               &+0.157*A(203)+0.157*A(207)+0.393*A(212)+0.002*A(214)+0.345*A(219)+0.265*A(223)+0.012*A(225)+1.5*A(226)+0.51&
               &*A(228)
  Vdot(70) = 2*A(13)+A(25)-A(27)-A(28)+0.2*A(39)+A(129)+A(132)+A(135)+A(149)+A(152)+A(155)+A(157)+A(158)+A(161)+0.5&
               &*A(164)+0.15*A(173)+A(233)
  Vdot(71) = -A(186)-A(187)-A(188)-A(189)
  Vdot(72) = A(116)-A(136)-A(137)+0.006*A(177)+0.02*A(178)+0.13*A(203)+0.13*A(207)+0.704*A(211)+0.024*A(213)+0.452&
               &*A(214)+0.072*A(215)+0.005*A(218)+0.001*A(219)+0.024*A(220)+0.127*A(222)+0.045*A(223)+0.102*A(224)
  Vdot(73) = -A(227)-A(228)-A(229)-A(230)
  Vdot(74) = -A(203)-A(204)-A(205)
  Vdot(75) = -A(146)-A(147)-A(148)-A(149)+0.23*A(154)+0.15*A(171)+0.023*A(172)+A(180)+0.5*A(182)+0.5*A(184)+0.009*A(189)&
               &+0.001*A(203)+0.001*A(207)+0.607*A(212)+0.118*A(216)+0.097*A(217)
  Vdot(76) = A(92)+A(94)+A(99)+A(100)+2*A(101)+A(113)-A(117)-A(118)-A(119)+0.24*A(154)+A(155)+0.24*A(156)+A(157)
  Vdot(77) = -A(218)-A(219)-A(220)-A(221)
  Vdot(78) = -A(56)-A(57)-A(58)-A(59)-A(60)-A(67)-A(76)-A(86)+A(92)+A(94)-A(97)+A(99)+A(100)+2*A(101)-A(109)+A(113)&
               &+A(136)+0.616*A(138)+0.675*A(167)+0.515*A(176)+0.596*A(177)+0.152*A(178)+A(181)+A(182)+A(183)+A(184)+0.079&
               &*A(190)+0.126*A(200)+0.187*A(201)+0.24*A(202)+0.729*A(203)+0.75*A(204)+0.5*A(206)+0.729*A(207)+0.75*A(208)&
               &+0.559*A(213)+0.936*A(214)+0.948*A(215)+0.205*A(218)+0.488*A(220)+0.001*A(222)+0.137*A(223)+0.711*A(224)
  Vdot(79) = -A(206)-A(207)-A(208)-A(209)
  Vdot(80) = -A(190)-A(200)-A(201)-A(202)
  Vdot(81) = -A(162)-A(163)-A(164)-A(165)-A(166)+0.23*A(190)+0.39*A(200)+0.025*A(222)+0.026*A(223)+0.012*A(225)
  Vdot(82) = -A(171)-A(172)-A(173)-A(174)+0.357*A(190)+0.936*A(201)+0.025*A(222)
  Vdot(83) = -A(222)-A(223)-A(224)-A(225)
  Vdot(84) = -A(167)-A(168)-A(169)-A(170)+0.32*A(190)+0.16*A(200)+0.019*A(223)+0.048*A(224)
  Vdot(85) = A(81)+A(83)+A(88)+2*A(89)+A(100)+A(112)-A(130)-A(131)-A(132)+0.034*A(133)+A(134)+0.482*A(138)+A(139)+0.96&
               &*A(141)+0.129*A(171)+0.047*A(172)+0.467*A(174)+0.084*A(175)+0.246*A(176)+0.439*A(177)+0.431*A(178)+0.195&
               &*A(186)+0.25*A(189)+A(210)+0.445*A(213)+0.455*A(214)+0.099*A(215)+0.294*A(218)+0.154*A(219)+0.009*A(220)&
               &+0.732*A(222)+0.456*A(223)+0.507*A(224)+0.984*A(227)+0.5*A(228)
  Vdot(86) = A(46)+A(48)+A(49)+2*A(50)+0.75*A(54)+0.75*A(64)+A(74)+A(84)+A(95)+A(104)+A(106)+A(107)+A(111)+A(112)+A(113)&
               &+2*A(114)-A(123)-A(124)-A(125)-A(126)+A(127)-A(129)+A(136)+0.115*A(138)+A(140)+0.081*A(141)+0.35*A(142)&
               &+A(143)+A(147)+0.084*A(162)+0.2*A(163)+0.67*A(166)+0.3*A(167)+0.1*A(168)+0.055*A(171)+0.125*A(172)+0.227&
               &*A(173)+0.3*A(174)+0.213*A(175)+0.506*A(176)+0.01*A(177)+0.134*A(178)+1.61*A(186)+A(187)+0.191*A(189)+0.624&
               &*A(190)+0.592*A(200)+0.24*A(202)+0.235*A(203)+0.276*A(206)+0.235*A(207)+0.096*A(212)+0.026*A(213)+0.024&
               &*A(214)+0.026*A(215)+0.732*A(218)+0.5*A(219)+0.244*A(222)+0.269*A(223)+0.079*A(224)+0.984*A(227)+0.5*A(228)
  Vdot(87) = A(62)+A(115)+0.572*A(173)-0.69*A(177)-A(178)+0.276*A(204)+0.276*A(208)+0.511*A(220)+0.321*A(224)
  Vdot(88) = A(1)-A(2)-A(3)-A(4)-A(5)-A(6)+A(16)+A(17)+A(20)-A(165)-A(169)-A(189)-A(202)-A(205)-A(209)-A(221)-A(225)&
               &-A(230)
  Vdot(89) = -A(133)-A(134)-A(135)+0.37*A(138)+A(144)+A(145)+A(165)+0.675*A(167)+0.45*A(169)+0.013*A(171)+0.218*A(173)&
               &+0.558*A(175)+0.71*A(176)+0.213*A(177)+0.147*A(178)+A(179)+A(181)+A(183)+A(188)+0.205*A(203)+0.474*A(204)&
               &+0.147*A(205)+0.474*A(206)+0.205*A(207)+0.474*A(208)+0.147*A(209)+0.261*A(211)+0.122*A(213)+0.244*A(214)&
               &+0.204*A(215)+0.497*A(218)+0.363*A(219)+0.037*A(220)+0.45*A(221)+0.511*A(222)+0.305*A(223)+0.151*A(224)&
               &+0.069*A(225)+0.45*A(230)
  Vdot(90) = 0.5*A(64)+A(65)+0.5*A(66)+A(68)-A(138)-A(139)+0.416*A(162)+0.55*A(169)+0.15*A(171)+0.21*A(172)+0.233*A(174)&
               &+0.115*A(175)+0.177*A(177)+0.243*A(178)+0.332*A(213)+0.11*A(214)+0.089*A(215)+0.437*A(221)+0.072*A(222)&
               &+0.026*A(223)+0.001*A(224)+0.659*A(225)+0.55*A(230)
  Vdot(91) = 0.5*A(64)+0.5*A(66)+A(68)+A(77)+A(87)+A(98)+0.7*A(170)+0.332*A(171)-0.671*A(175)-A(176)+0.048*A(177)+0.435&
               &*A(178)+0.1*A(200)+0.75*A(202)+0.276*A(203)+0.853*A(205)+0.276*A(206)+0.276*A(207)+0.853*A(209)+0.125*A(214)&
               &+0.417*A(215)+0.055*A(216)+0.119*A(218)+0.215*A(219)+0.113*A(221)+0.043*A(223)+0.259*A(225)
  Vdot(92) = A(2)-A(3)-A(7)-A(8)-A(17)-A(18)-A(30)-A(36)+0.25*A(72)+0.25*A(82)+0.25*A(93)+0.25*A(105)-A(163)-A(168)&
               &-A(172)-A(180)-A(187)-A(200)-A(203)-A(207)-A(219)-A(223)-A(226)-A(228)
  Vdot(93) = A(1)-A(4)+A(5)-A(7)-A(9)-2*A(10)+A(14)+A(15)-A(21)+A(22)-A(31)-A(46)-A(51)-A(56)-A(62)-A(71)-A(81)-A(92)&
               &-A(104)-A(128)
  Vdot(94) = -A(1)+A(4)-A(5)-A(6)+A(7)-A(8)+2*A(9)+2*A(10)-A(11)+A(12)+A(16)+A(23)+A(24)-A(25)+A(26)+A(28)+A(31)-A(32)&
               &+A(33)+0.61*A(34)+A(35)+0.8*A(39)+2*A(40)+A(46)+A(48)+A(51)+A(53)+A(56)+A(58)+A(65)-A(69)+A(70)+A(71)+A(73)&
               &-A(79)+A(80)+A(81)+A(83)-A(90)+A(91)+A(92)+A(94)-A(102)+A(103)+A(104)+A(106)-A(115)-A(117)-A(120)+A(128)&
               &+0.338*A(177)+A(178)+0.187*A(201)+0.474*A(204)+0.474*A(208)+0.391*A(224)
  Vdot(95) = A(6)+A(8)-A(9)-A(11)+A(12)-A(14)-A(15)-A(16)-A(26)+A(27)+0.39*A(34)-A(39)-2*A(40)-A(48)-A(53)-A(58)-A(65)&
               &-A(73)-A(83)-A(94)-A(106)-A(129)-A(132)-A(135)-A(149)-A(152)-A(155)-A(157)-A(158)-A(161)-A(164)-A(173)&
               &-A(188)-A(201)-A(204)-A(208)-A(220)-A(224)-A(229)-A(233)
  Vdot(96) = -A(90)+A(91)-A(92)-A(93)-A(94)-A(95)-A(96)-A(98)-A(99)-A(100)-2*A(101)-A(113)+A(159)+A(161)
  Vdot(97) = -A(79)+A(80)-A(81)-A(82)-A(83)-A(84)-A(85)-A(87)-A(88)-2*A(89)-A(100)-A(112)+0.965*A(133)+A(135)+0.096&
               &*A(138)+0.37*A(148)+0.37*A(149)+0.1*A(163)+0.05*A(168)+0.048*A(172)+0.3*A(174)+0.049*A(175)+0.333*A(176)&
               &+0.201*A(203)+0.201*A(207)+0.006*A(223)
  Vdot(98) = -A(102)+A(103)-A(104)-A(105)-A(106)-A(107)-A(108)-A(110)-A(111)-A(112)-A(113)-2*A(114)+0.5*A(162)+0.5&
               &*A(164)+0.33*A(166)+0.3*A(170)+0.289*A(171)+0.15*A(173)+0.192*A(200)+0.24*A(202)
  Vdot(99) = -A(46)-A(47)-A(48)-2*A(49)-2*A(50)-A(54)-A(64)+A(71)+A(73)-A(74)+2*A(78)-A(84)+A(88)-A(95)+A(99)-A(107)&
               &+A(111)+A(116)+A(131)+A(137)+0.65*A(142)+0.3*A(170)+A(185)+0.3*A(189)+0.25*A(202)+0.011*A(214)+0.076*A(219)&
               &+0.197*A(223)+0.03*A(224)+0.26*A(228)
  Vdot(100) = -A(69)+A(70)-A(71)-A(72)-A(73)-A(74)-A(75)-A(77)-2*A(78)-A(88)-A(99)+A(104)+A(106)+A(112)+A(113)+2*A(114)&
                &+A(130)+A(132)+A(136)+A(137)+0.492*A(138)+A(139)+A(150)+A(151)+A(152)+2*A(153)+0.67*A(166)+0.675*A(167)&
                &+0.467*A(174)+0.029*A(175)+0.667*A(176)+A(181)+0.5*A(182)+A(183)+0.5*A(184)+0.123*A(203)+0.123*A(207)+0.011&
                &*A(214)+0.137*A(223)
  Vdot(101) = -A(62)-A(63)-A(64)-A(65)-A(66)-2*A(68)-A(77)-A(87)-A(98)-A(110)+0.001*A(133)+0.042*A(138)+0.025*A(167)&
                &+0.041*A(171)+0.051*A(173)+0.07*A(175)+0.04*A(176)+0.173*A(177)+0.095*A(178)+0.093*A(190)+0.008*A(200)&
                &+0.064*A(201)+0.01*A(202)+0.18*A(203)+0.25*A(204)+0.25*A(206)+0.18*A(207)+0.25*A(208)+0.035*A(211)+0.07&
                &*A(213)+0.143*A(214)+0.347*A(215)+0.011*A(216)+0.009*A(217)+0.09*A(218)+0.001*A(219)+0.176*A(220)+0.082&
                &*A(222)+0.002*A(223)+0.136*A(224)+0.001*A(225)+0.016*A(227)+0.051*A(229)
  Vdot(102) = 2*A(19)-A(21)+A(22)-A(24)-A(25)-A(26)-A(27)+A(28)-A(29)-A(30)+A(31)+0.39*A(34)-A(35)+A(36)+0.8*A(39)+2&
                &*A(41)-A(42)-A(43)-A(44)-A(45)-A(125)-A(130)-A(133)-A(136)-A(138)-A(140)-A(141)-0.65*A(142)+A(143)-0.34&
                &*A(144)+A(145)-A(148)-A(151)-A(154)-A(156)-A(159)-A(162)+0.208*A(163)+0.33*A(166)-A(167)+0.164*A(168)&
                &-A(171)+0.285*A(172)-A(175)-A(177)-A(179)+0.5*A(180)-A(181)-A(183)-A(185)-A(186)+0.12*A(187)-A(190)+0.266&
                &*A(200)+0.567*A(203)-A(206)+0.567*A(207)-A(210)-A(211)-0.397*A(212)-A(213)-A(214)-A(215)-A(216)-A(217)&
                &-A(218)+0.155*A(219)-A(222)+0.378*A(223)+0.5*A(226)-A(227)+0.32*A(228)-A(231)-A(232)
  Vdot(103) = A(23)+A(26)+A(29)+A(30)-A(31)-A(32)+A(33)+0.61*A(34)-A(36)-2*A(37)-2*A(38)-A(39)+A(42)-A(43)+A(44)+A(45)&
                &+A(46)-A(47)+A(48)+2*A(50)+A(51)-A(52)+A(53)+A(54)+A(55)-A(63)+A(64)+A(65)+A(66)+A(68)-A(72)-A(82)-A(93)&
                &-A(105)-A(118)-A(121)+2*A(123)+A(125)-A(126)+A(127)+A(128)+A(129)+A(131)+A(134)+A(140)+0.95*A(141)+A(143)&
                &+A(145)+2*A(146)+0.63*A(148)+0.63*A(149)+A(150)+0.008*A(163)+0.34*A(166)+0.064*A(168)+0.4*A(172)+1.233&
                &*A(174)+0.379*A(175)+0.113*A(177)+0.341*A(178)+1.5*A(180)+0.5*A(182)+0.5*A(184)+0.12*A(187)+0.5*A(189)&
                &+0.907*A(190)+0.033*A(203)+0.033*A(207)+0.297*A(212)+0.224*A(216)+0.187*A(217)+0.056*A(219)+0.003*A(223)&
                &+0.013*A(225)+1.5*A(226)+0.06*A(228)+0.5*A(232)
  Vdot(104) = -A(51)-A(52)-A(53)-A(54)-2*A(55)-A(66)-A(75)+A(81)+A(83)-A(85)+A(88)+2*A(89)-A(96)+A(100)-A(108)+A(112)&
                &+0.034*A(133)+A(134)+0.37*A(138)+A(139)+0.05*A(141)+0.34*A(144)+0.76*A(154)+0.76*A(156)+0.5*A(162)+0.1&
                &*A(163)+0.5*A(164)+0.33*A(166)+0.3*A(167)+0.05*A(168)+0.67*A(171)+0.048*A(172)+0.799*A(173)+0.473*A(175)&
                &+0.96*A(176)+0.376*A(177)+0.564*A(178)+A(179)+A(182)+A(184)+A(186)+A(188)+0.2*A(189)+0.907*A(190)+0.066&
                &*A(200)+0.749*A(201)+0.031*A(203)+0.276*A(204)+0.75*A(206)+0.031*A(207)+0.276*A(208)+A(210)+0.965*A(211)&
                &+0.1*A(212)+0.695*A(213)+0.835*A(214)+0.653*A(215)+0.765*A(216)+0.804*A(217)+0.91*A(218)+0.022*A(219)+0.824&
                &*A(220)+0.918*A(222)+0.033*A(223)+0.442*A(224)+0.012*A(225)+0.984*A(227)+0.949*A(229)
      
END SUBROUTINE saprc99_simplesom_mosaic_4bin_aq_Fun
















SUBROUTINE saprc99_simplesom_mosaic_4bin_aq_IRRFun ( V, F, RCT, IRR )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: IRR(NREACT)



  IRR(1) = RCT(1)*V(94)
  IRR(2) = RCT(2)*V(88)*F(2)
  IRR(3) = RCT(3)*V(88)*V(92)
  IRR(4) = RCT(4)*V(88)*V(93)*F(2)
  IRR(5) = RCT(5)*V(88)*V(94)
  IRR(6) = RCT(6)*V(88)*V(94)
  IRR(7) = RCT(7)*V(92)*V(93)
  IRR(8) = RCT(8)*V(92)*V(94)
  IRR(9) = RCT(9)*V(93)*V(95)
  IRR(10) = RCT(10)*V(93)*V(93)*F(2)
  IRR(11) = RCT(11)*V(94)*V(95)
  IRR(12) = RCT(12)*V(31)
  IRR(13) = RCT(13)*V(31)*F(1)
  IRR(14) = RCT(14)*V(94)*V(95)
  IRR(15) = RCT(15)*V(95)
  IRR(16) = RCT(16)*V(95)
  IRR(17) = RCT(17)*V(92)
  IRR(18) = RCT(18)*V(92)
  IRR(19) = RCT(19)*V(18)*F(1)
  IRR(20) = RCT(20)*V(18)*F(2)
  IRR(21) = RCT(21)*V(93)*V(102)
  IRR(22) = RCT(22)*V(33)
  IRR(23) = RCT(23)*V(33)
  IRR(24) = RCT(24)*V(33)*V(102)
  IRR(25) = RCT(25)*V(94)*V(102)
  IRR(26) = RCT(26)*V(95)*V(102)
  IRR(27) = RCT(27)*V(70)*V(102)
  IRR(28) = RCT(28)*V(70)
  IRR(29) = RCT(29)*V(69)*V(102)
  IRR(30) = RCT(30)*V(92)*V(102)
  IRR(31) = RCT(31)*V(93)*V(103)
  IRR(32) = RCT(32)*V(94)*V(103)
  IRR(33) = RCT(33)*V(43)
  IRR(34) = RCT(34)*V(43)
  IRR(35) = RCT(35)*V(43)*V(102)
  IRR(36) = RCT(36)*V(92)*V(103)
  IRR(37) = RCT(37)*V(103)*V(103)
  IRR(38) = RCT(38)*V(103)*V(103)*F(1)
  IRR(39) = RCT(39)*V(95)*V(103)
  IRR(40) = RCT(40)*V(95)*V(95)
  IRR(41) = RCT(41)*V(27)
  IRR(42) = RCT(42)*V(27)*V(102)
  IRR(43) = RCT(43)*V(102)*V(103)
  IRR(44) = RCT(44)*V(32)*V(102)
  IRR(45) = RCT(45)*V(102)*F(2)
  IRR(46) = RCT(46)*V(93)*V(99)
  IRR(47) = RCT(47)*V(99)*V(103)
  IRR(48) = RCT(48)*V(95)*V(99)
  IRR(49) = RCT(49)*V(99)*V(99)
  IRR(50) = RCT(50)*V(99)*V(99)
  IRR(51) = RCT(51)*V(93)*V(104)
  IRR(52) = RCT(52)*V(103)*V(104)
  IRR(53) = RCT(53)*V(95)*V(104)
  IRR(54) = RCT(54)*V(99)*V(104)
  IRR(55) = RCT(55)*V(104)*V(104)
  IRR(56) = RCT(56)*V(78)*V(93)
  IRR(57) = RCT(57)*V(78)*V(103)
  IRR(58) = RCT(58)*V(78)*V(95)
  IRR(59) = RCT(59)*V(78)*V(99)
  IRR(60) = RCT(60)*V(78)*V(104)
  IRR(61) = RCT(61)*V(78)*V(78)
  IRR(62) = RCT(62)*V(93)*V(101)
  IRR(63) = RCT(63)*V(101)*V(103)
  IRR(64) = RCT(64)*V(99)*V(101)
  IRR(65) = RCT(65)*V(95)*V(101)
  IRR(66) = RCT(66)*V(101)*V(104)
  IRR(67) = RCT(67)*V(78)*V(101)
  IRR(68) = RCT(68)*V(101)*V(101)
  IRR(69) = RCT(69)*V(94)*V(100)
  IRR(70) = RCT(70)*V(23)
  IRR(71) = RCT(71)*V(93)*V(100)
  IRR(72) = RCT(72)*V(100)*V(103)
  IRR(73) = RCT(73)*V(95)*V(100)
  IRR(74) = RCT(74)*V(99)*V(100)
  IRR(75) = RCT(75)*V(100)*V(104)
  IRR(76) = RCT(76)*V(78)*V(100)
  IRR(77) = RCT(77)*V(100)*V(101)
  IRR(78) = RCT(78)*V(100)*V(100)
  IRR(79) = RCT(79)*V(94)*V(97)
  IRR(80) = RCT(80)*V(24)
  IRR(81) = RCT(81)*V(93)*V(97)
  IRR(82) = RCT(82)*V(97)*V(103)
  IRR(83) = RCT(83)*V(95)*V(97)
  IRR(84) = RCT(84)*V(97)*V(99)
  IRR(85) = RCT(85)*V(97)*V(104)
  IRR(86) = RCT(86)*V(78)*V(97)
  IRR(87) = RCT(87)*V(97)*V(101)
  IRR(88) = RCT(88)*V(97)*V(100)
  IRR(89) = RCT(89)*V(97)*V(97)
  IRR(90) = RCT(90)*V(94)*V(96)
  IRR(91) = RCT(91)*V(25)
  IRR(92) = RCT(92)*V(93)*V(96)
  IRR(93) = RCT(93)*V(96)*V(103)
  IRR(94) = RCT(94)*V(95)*V(96)
  IRR(95) = RCT(95)*V(96)*V(99)
  IRR(96) = RCT(96)*V(96)*V(104)
  IRR(97) = RCT(97)*V(78)*V(96)
  IRR(98) = RCT(98)*V(96)*V(101)
  IRR(99) = RCT(99)*V(96)*V(100)
  IRR(100) = RCT(100)*V(96)*V(97)
  IRR(101) = RCT(101)*V(96)*V(96)
  IRR(102) = RCT(102)*V(94)*V(98)
  IRR(103) = RCT(103)*V(26)
  IRR(104) = RCT(104)*V(93)*V(98)
  IRR(105) = RCT(105)*V(98)*V(103)
  IRR(106) = RCT(106)*V(95)*V(98)
  IRR(107) = RCT(107)*V(98)*V(99)
  IRR(108) = RCT(108)*V(98)*V(104)
  IRR(109) = RCT(109)*V(78)*V(98)
  IRR(110) = RCT(110)*V(98)*V(101)
  IRR(111) = RCT(111)*V(98)*V(100)
  IRR(112) = RCT(112)*V(97)*V(98)
  IRR(113) = RCT(113)*V(96)*V(98)
  IRR(114) = RCT(114)*V(98)*V(98)
  IRR(115) = RCT(115)*V(36)*V(94)
  IRR(116) = RCT(116)*V(36)
  IRR(117) = RCT(117)*V(76)*V(94)
  IRR(118) = RCT(118)*V(76)*V(103)
  IRR(119) = RCT(119)*V(76)
  IRR(120) = RCT(120)*V(42)*V(94)
  IRR(121) = RCT(121)*V(42)*V(103)
  IRR(122) = RCT(122)*V(42)
  IRR(123) = RCT(123)*V(86)
  IRR(124) = RCT(124)*V(86)
  IRR(125) = RCT(125)*V(86)*V(102)
  IRR(126) = RCT(126)*V(86)*V(103)
  IRR(127) = RCT(127)*V(41)
  IRR(128) = RCT(128)*V(41)*V(93)
  IRR(129) = RCT(129)*V(86)*V(95)
  IRR(130) = RCT(130)*V(85)*V(102)
  IRR(131) = RCT(131)*V(85)
  IRR(132) = RCT(132)*V(85)*V(95)
  IRR(133) = RCT(133)*V(89)*V(102)
  IRR(134) = RCT(134)*V(89)
  IRR(135) = RCT(135)*V(89)*V(95)
  IRR(136) = RCT(136)*V(72)*V(102)
  IRR(137) = RCT(137)*V(72)
  IRR(138) = RCT(138)*V(90)*V(102)
  IRR(139) = RCT(139)*V(90)
  IRR(140) = RCT(140)*V(39)*V(102)
  IRR(141) = RCT(141)*V(30)*V(102)
  IRR(142) = RCT(142)*V(40)*V(102)
  IRR(143) = RCT(143)*V(40)
  IRR(144) = RCT(144)*V(53)*V(102)
  IRR(145) = RCT(145)*V(53)
  IRR(146) = RCT(146)*V(75)
  IRR(147) = RCT(147)*V(75)
  IRR(148) = RCT(148)*V(75)*V(102)
  IRR(149) = RCT(149)*V(75)*V(95)
  IRR(150) = RCT(150)*V(68)
  IRR(151) = RCT(151)*V(68)*V(102)
  IRR(152) = RCT(152)*V(68)*V(95)
  IRR(153) = RCT(153)*V(28)
  IRR(154) = RCT(154)*V(55)*V(102)
  IRR(155) = RCT(155)*V(55)*V(95)
  IRR(156) = RCT(156)*V(48)*V(102)
  IRR(157) = RCT(157)*V(48)*V(95)
  IRR(158) = RCT(158)*V(51)*V(95)
  IRR(159) = RCT(159)*V(52)*V(102)
  IRR(160) = RCT(160)*V(52)
  IRR(161) = RCT(161)*V(52)*V(95)
  IRR(162) = RCT(162)*V(81)*V(102)
  IRR(163) = RCT(163)*V(81)*V(92)
  IRR(164) = RCT(164)*V(81)*V(95)
  IRR(165) = RCT(165)*V(81)*V(88)
  IRR(166) = RCT(166)*V(81)
  IRR(167) = RCT(167)*V(84)*V(102)
  IRR(168) = RCT(168)*V(84)*V(92)
  IRR(169) = RCT(169)*V(84)*V(88)
  IRR(170) = RCT(170)*V(84)
  IRR(171) = RCT(171)*V(82)*V(102)
  IRR(172) = RCT(172)*V(82)*V(92)
  IRR(173) = RCT(173)*V(82)*V(95)
  IRR(174) = RCT(174)*V(82)
  IRR(175) = RCT(175)*V(91)*V(102)
  IRR(176) = RCT(176)*V(91)
  IRR(177) = RCT(177)*V(87)*V(102)
  IRR(178) = RCT(178)*V(87)
  IRR(179) = RCT(179)*V(49)*V(102)
  IRR(180) = RCT(180)*V(49)*V(92)
  IRR(181) = RCT(181)*V(46)*V(102)
  IRR(182) = RCT(182)*V(46)
  IRR(183) = RCT(183)*V(47)*V(102)
  IRR(184) = RCT(184)*V(47)
  IRR(185) = RCT(185)*V(19)*V(102)
  IRR(186) = RCT(186)*V(71)*V(102)
  IRR(187) = RCT(187)*V(71)*V(92)
  IRR(188) = RCT(188)*V(71)*V(95)
  IRR(189) = RCT(189)*V(71)*V(88)
  IRR(190) = RCT(190)*V(80)*V(102)
  IRR(191) = RCT(191)*V(21)*V(103)
  IRR(192) = RCT(192)*V(21)*V(93)
  IRR(193) = RCT(193)*V(21)*V(21)
  IRR(194) = RCT(194)*V(20)*V(102)
  IRR(195) = RCT(195)*V(16)*V(102)
  IRR(196) = RCT(196)*V(17)*V(102)
  IRR(197) = RCT(197)*V(15)*V(103)
  IRR(198) = RCT(198)*V(15)*V(93)
  IRR(199) = RCT(199)*V(15)
  IRR(200) = RCT(200)*V(80)*V(92)
  IRR(201) = RCT(201)*V(80)*V(95)
  IRR(202) = RCT(202)*V(80)*V(88)
  IRR(203) = RCT(203)*V(74)*V(92)
  IRR(204) = RCT(204)*V(74)*V(95)
  IRR(205) = RCT(205)*V(74)*V(88)
  IRR(206) = RCT(206)*V(79)*V(102)
  IRR(207) = RCT(207)*V(79)*V(92)
  IRR(208) = RCT(208)*V(79)*V(95)
  IRR(209) = RCT(209)*V(79)*V(88)
  IRR(210) = RCT(210)*V(22)*V(102)
  IRR(211) = RCT(211)*V(29)*V(102)
  IRR(212) = RCT(212)*V(50)*V(102)
  IRR(213) = RCT(213)*V(35)*V(102)
  IRR(214) = RCT(214)*V(44)*V(102)
  IRR(215) = RCT(215)*V(37)*V(102)
  IRR(216) = RCT(216)*V(45)*V(102)
  IRR(217) = RCT(217)*V(38)*V(102)
  IRR(218) = RCT(218)*V(77)*V(102)
  IRR(219) = RCT(219)*V(77)*V(92)
  IRR(220) = RCT(220)*V(77)*V(95)
  IRR(221) = RCT(221)*V(77)*V(88)
  IRR(222) = RCT(222)*V(83)*V(102)
  IRR(223) = RCT(223)*V(83)*V(92)
  IRR(224) = RCT(224)*V(83)*V(95)
  IRR(225) = RCT(225)*V(83)*V(88)
  IRR(226) = RCT(226)*V(50)*V(92)
  IRR(227) = RCT(227)*V(73)*V(102)
  IRR(228) = RCT(228)*V(73)*V(92)
  IRR(229) = RCT(229)*V(73)*V(95)
  IRR(230) = RCT(230)*V(73)*V(88)
  IRR(231) = RCT(231)*V(34)*V(102)
  IRR(232) = RCT(232)*V(34)*V(102)
  IRR(233) = RCT(233)*V(34)*V(95)
  IRR(234) = RCT(234)*V(64)*V(102)
  IRR(235) = RCT(235)*V(67)*V(102)
  IRR(236) = RCT(236)*V(66)*V(102)
  IRR(237) = RCT(237)*V(65)*V(102)
  IRR(238) = RCT(238)*V(63)*V(102)
  IRR(239) = RCT(239)*V(62)*V(102)
  IRR(240) = RCT(240)*V(61)*V(102)
  IRR(241) = RCT(241)*V(60)*V(102)
  IRR(242) = RCT(242)*V(59)*V(102)
  IRR(243) = RCT(243)*V(58)*V(102)
  IRR(244) = RCT(244)*V(56)*V(102)
  IRR(245) = RCT(245)*V(54)*V(102)
  IRR(246) = RCT(246)*V(57)*V(102)
  IRR(247) = RCT(247)*V(74)*V(102)
      
END SUBROUTINE saprc99_simplesom_mosaic_4bin_aq_IRRFun
















SUBROUTINE saprc99_simplesom_mosaic_4bin_aq_Jac_SP ( V, F, RCT, JVS )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: JVS(LU_NONZERO)




  REAL(kind=dp) :: B(442)


  B(1) = RCT(1)

  B(2) = RCT(2)*F(2)

  B(4) = RCT(3)*V(92)

  B(5) = RCT(3)*V(88)

  B(6) = RCT(4)*V(93)*F(2)

  B(7) = RCT(4)*V(88)*F(2)

  B(9) = RCT(5)*V(94)

  B(10) = RCT(5)*V(88)

  B(11) = RCT(6)*V(94)

  B(12) = RCT(6)*V(88)

  B(13) = RCT(7)*V(93)

  B(14) = RCT(7)*V(92)

  B(15) = RCT(8)*V(94)

  B(16) = RCT(8)*V(92)

  B(17) = RCT(9)*V(95)

  B(18) = RCT(9)*V(93)

  B(19) = RCT(10)*2*V(93)*F(2)

  B(21) = RCT(11)*V(95)

  B(22) = RCT(11)*V(94)

  B(23) = RCT(12)

  B(24) = RCT(13)*F(1)

  B(26) = RCT(14)*V(95)

  B(27) = RCT(14)*V(94)

  B(28) = RCT(15)

  B(29) = RCT(16)

  B(30) = RCT(17)

  B(31) = RCT(18)

  B(32) = RCT(19)*F(1)

  B(34) = RCT(20)*F(2)

  B(36) = RCT(21)*V(102)

  B(37) = RCT(21)*V(93)

  B(38) = RCT(22)

  B(39) = RCT(23)

  B(40) = RCT(24)*V(102)

  B(41) = RCT(24)*V(33)

  B(42) = RCT(25)*V(102)

  B(43) = RCT(25)*V(94)

  B(44) = RCT(26)*V(102)

  B(45) = RCT(26)*V(95)

  B(46) = RCT(27)*V(102)

  B(47) = RCT(27)*V(70)

  B(48) = RCT(28)

  B(49) = RCT(29)*V(102)

  B(50) = RCT(29)*V(69)

  B(51) = RCT(30)*V(102)

  B(52) = RCT(30)*V(92)

  B(53) = RCT(31)*V(103)

  B(54) = RCT(31)*V(93)

  B(55) = RCT(32)*V(103)

  B(56) = RCT(32)*V(94)

  B(57) = RCT(33)

  B(58) = RCT(34)

  B(59) = RCT(35)*V(102)

  B(60) = RCT(35)*V(43)

  B(61) = RCT(36)*V(103)

  B(62) = RCT(36)*V(92)

  B(63) = RCT(37)*2*V(103)

  B(64) = RCT(38)*2*V(103)*F(1)

  B(66) = RCT(39)*V(103)

  B(67) = RCT(39)*V(95)

  B(68) = RCT(40)*2*V(95)

  B(69) = RCT(41)

  B(70) = RCT(42)*V(102)

  B(71) = RCT(42)*V(27)

  B(72) = RCT(43)*V(103)

  B(73) = RCT(43)*V(102)

  B(74) = RCT(44)*V(102)

  B(75) = RCT(44)*V(32)

  B(76) = RCT(45)*F(2)

  B(78) = RCT(46)*V(99)

  B(79) = RCT(46)*V(93)

  B(80) = RCT(47)*V(103)

  B(81) = RCT(47)*V(99)

  B(82) = RCT(48)*V(99)

  B(83) = RCT(48)*V(95)

  B(84) = RCT(49)*2*V(99)

  B(85) = RCT(50)*2*V(99)

  B(86) = RCT(51)*V(104)

  B(87) = RCT(51)*V(93)

  B(88) = RCT(52)*V(104)

  B(89) = RCT(52)*V(103)

  B(90) = RCT(53)*V(104)

  B(91) = RCT(53)*V(95)

  B(92) = RCT(54)*V(104)

  B(93) = RCT(54)*V(99)

  B(94) = RCT(55)*2*V(104)

  B(95) = RCT(56)*V(93)

  B(96) = RCT(56)*V(78)

  B(97) = RCT(57)*V(103)

  B(98) = RCT(57)*V(78)

  B(99) = RCT(58)*V(95)

  B(100) = RCT(58)*V(78)

  B(101) = RCT(59)*V(99)

  B(102) = RCT(59)*V(78)

  B(103) = RCT(60)*V(104)

  B(104) = RCT(60)*V(78)

  B(105) = RCT(61)*2*V(78)

  B(106) = RCT(62)*V(101)

  B(107) = RCT(62)*V(93)

  B(108) = RCT(63)*V(103)

  B(109) = RCT(63)*V(101)

  B(110) = RCT(64)*V(101)

  B(111) = RCT(64)*V(99)

  B(112) = RCT(65)*V(101)

  B(113) = RCT(65)*V(95)

  B(114) = RCT(66)*V(104)

  B(115) = RCT(66)*V(101)

  B(116) = RCT(67)*V(101)

  B(117) = RCT(67)*V(78)

  B(118) = RCT(68)*2*V(101)

  B(119) = RCT(69)*V(100)

  B(120) = RCT(69)*V(94)

  B(121) = RCT(70)

  B(122) = RCT(71)*V(100)

  B(123) = RCT(71)*V(93)

  B(124) = RCT(72)*V(103)

  B(125) = RCT(72)*V(100)

  B(126) = RCT(73)*V(100)

  B(127) = RCT(73)*V(95)

  B(128) = RCT(74)*V(100)

  B(129) = RCT(74)*V(99)

  B(130) = RCT(75)*V(104)

  B(131) = RCT(75)*V(100)

  B(132) = RCT(76)*V(100)

  B(133) = RCT(76)*V(78)

  B(134) = RCT(77)*V(101)

  B(135) = RCT(77)*V(100)

  B(136) = RCT(78)*2*V(100)

  B(137) = RCT(79)*V(97)

  B(138) = RCT(79)*V(94)

  B(139) = RCT(80)

  B(140) = RCT(81)*V(97)

  B(141) = RCT(81)*V(93)

  B(142) = RCT(82)*V(103)

  B(143) = RCT(82)*V(97)

  B(144) = RCT(83)*V(97)

  B(145) = RCT(83)*V(95)

  B(146) = RCT(84)*V(99)

  B(147) = RCT(84)*V(97)

  B(148) = RCT(85)*V(104)

  B(149) = RCT(85)*V(97)

  B(150) = RCT(86)*V(97)

  B(151) = RCT(86)*V(78)

  B(152) = RCT(87)*V(101)

  B(153) = RCT(87)*V(97)

  B(154) = RCT(88)*V(100)

  B(155) = RCT(88)*V(97)

  B(156) = RCT(89)*2*V(97)

  B(157) = RCT(90)*V(96)

  B(158) = RCT(90)*V(94)

  B(159) = RCT(91)

  B(160) = RCT(92)*V(96)

  B(161) = RCT(92)*V(93)

  B(162) = RCT(93)*V(103)

  B(163) = RCT(93)*V(96)

  B(164) = RCT(94)*V(96)

  B(165) = RCT(94)*V(95)

  B(166) = RCT(95)*V(99)

  B(167) = RCT(95)*V(96)

  B(168) = RCT(96)*V(104)

  B(169) = RCT(96)*V(96)

  B(170) = RCT(97)*V(96)

  B(171) = RCT(97)*V(78)

  B(172) = RCT(98)*V(101)

  B(173) = RCT(98)*V(96)

  B(174) = RCT(99)*V(100)

  B(175) = RCT(99)*V(96)

  B(176) = RCT(100)*V(97)

  B(177) = RCT(100)*V(96)

  B(178) = RCT(101)*2*V(96)

  B(179) = RCT(102)*V(98)

  B(180) = RCT(102)*V(94)

  B(181) = RCT(103)

  B(182) = RCT(104)*V(98)

  B(183) = RCT(104)*V(93)

  B(184) = RCT(105)*V(103)

  B(185) = RCT(105)*V(98)

  B(186) = RCT(106)*V(98)

  B(187) = RCT(106)*V(95)

  B(188) = RCT(107)*V(99)

  B(189) = RCT(107)*V(98)

  B(190) = RCT(108)*V(104)

  B(191) = RCT(108)*V(98)

  B(192) = RCT(109)*V(98)

  B(193) = RCT(109)*V(78)

  B(194) = RCT(110)*V(101)

  B(195) = RCT(110)*V(98)

  B(196) = RCT(111)*V(100)

  B(197) = RCT(111)*V(98)

  B(198) = RCT(112)*V(98)

  B(199) = RCT(112)*V(97)

  B(200) = RCT(113)*V(98)

  B(201) = RCT(113)*V(96)

  B(202) = RCT(114)*2*V(98)

  B(203) = RCT(115)*V(94)

  B(204) = RCT(115)*V(36)

  B(205) = RCT(116)

  B(206) = RCT(117)*V(94)

  B(207) = RCT(117)*V(76)

  B(208) = RCT(118)*V(103)

  B(209) = RCT(118)*V(76)

  B(210) = RCT(119)

  B(211) = RCT(120)*V(94)

  B(212) = RCT(120)*V(42)

  B(213) = RCT(121)*V(103)

  B(214) = RCT(121)*V(42)

  B(215) = RCT(122)

  B(216) = RCT(123)

  B(217) = RCT(124)

  B(218) = RCT(125)*V(102)

  B(219) = RCT(125)*V(86)

  B(220) = RCT(126)*V(103)

  B(221) = RCT(126)*V(86)

  B(222) = RCT(127)

  B(223) = RCT(128)*V(93)

  B(224) = RCT(128)*V(41)

  B(225) = RCT(129)*V(95)

  B(226) = RCT(129)*V(86)

  B(227) = RCT(130)*V(102)

  B(228) = RCT(130)*V(85)

  B(229) = RCT(131)

  B(230) = RCT(132)*V(95)

  B(231) = RCT(132)*V(85)

  B(232) = RCT(133)*V(102)

  B(233) = RCT(133)*V(89)

  B(234) = RCT(134)

  B(235) = RCT(135)*V(95)

  B(236) = RCT(135)*V(89)

  B(237) = RCT(136)*V(102)

  B(238) = RCT(136)*V(72)

  B(239) = RCT(137)

  B(240) = RCT(138)*V(102)

  B(241) = RCT(138)*V(90)

  B(242) = RCT(139)

  B(243) = RCT(140)*V(102)

  B(244) = RCT(140)*V(39)

  B(245) = RCT(141)*V(102)

  B(246) = RCT(141)*V(30)

  B(247) = RCT(142)*V(102)

  B(248) = RCT(142)*V(40)

  B(249) = RCT(143)

  B(250) = RCT(144)*V(102)

  B(251) = RCT(144)*V(53)

  B(252) = RCT(145)

  B(253) = RCT(146)

  B(254) = RCT(147)

  B(255) = RCT(148)*V(102)

  B(256) = RCT(148)*V(75)

  B(257) = RCT(149)*V(95)

  B(258) = RCT(149)*V(75)

  B(259) = RCT(150)

  B(260) = RCT(151)*V(102)

  B(261) = RCT(151)*V(68)

  B(262) = RCT(152)*V(95)

  B(263) = RCT(152)*V(68)

  B(264) = RCT(153)

  B(265) = RCT(154)*V(102)

  B(266) = RCT(154)*V(55)

  B(267) = RCT(155)*V(95)

  B(268) = RCT(155)*V(55)

  B(269) = RCT(156)*V(102)

  B(270) = RCT(156)*V(48)

  B(271) = RCT(157)*V(95)

  B(272) = RCT(157)*V(48)

  B(273) = RCT(158)*V(95)

  B(274) = RCT(158)*V(51)

  B(275) = RCT(159)*V(102)

  B(276) = RCT(159)*V(52)

  B(277) = RCT(160)

  B(278) = RCT(161)*V(95)

  B(279) = RCT(161)*V(52)

  B(280) = RCT(162)*V(102)

  B(281) = RCT(162)*V(81)

  B(282) = RCT(163)*V(92)

  B(283) = RCT(163)*V(81)

  B(284) = RCT(164)*V(95)

  B(285) = RCT(164)*V(81)

  B(286) = RCT(165)*V(88)

  B(287) = RCT(165)*V(81)

  B(288) = RCT(166)

  B(289) = RCT(167)*V(102)

  B(290) = RCT(167)*V(84)

  B(291) = RCT(168)*V(92)

  B(292) = RCT(168)*V(84)

  B(293) = RCT(169)*V(88)

  B(294) = RCT(169)*V(84)

  B(295) = RCT(170)

  B(296) = RCT(171)*V(102)

  B(297) = RCT(171)*V(82)

  B(298) = RCT(172)*V(92)

  B(299) = RCT(172)*V(82)

  B(300) = RCT(173)*V(95)

  B(301) = RCT(173)*V(82)

  B(302) = RCT(174)

  B(303) = RCT(175)*V(102)

  B(304) = RCT(175)*V(91)

  B(305) = RCT(176)

  B(306) = RCT(177)*V(102)

  B(307) = RCT(177)*V(87)

  B(308) = RCT(178)

  B(309) = RCT(179)*V(102)

  B(310) = RCT(179)*V(49)

  B(311) = RCT(180)*V(92)

  B(312) = RCT(180)*V(49)

  B(313) = RCT(181)*V(102)

  B(314) = RCT(181)*V(46)

  B(315) = RCT(182)

  B(316) = RCT(183)*V(102)

  B(317) = RCT(183)*V(47)

  B(318) = RCT(184)

  B(319) = RCT(185)*V(102)

  B(320) = RCT(185)*V(19)

  B(321) = RCT(186)*V(102)

  B(322) = RCT(186)*V(71)

  B(323) = RCT(187)*V(92)

  B(324) = RCT(187)*V(71)

  B(325) = RCT(188)*V(95)

  B(326) = RCT(188)*V(71)

  B(327) = RCT(189)*V(88)

  B(328) = RCT(189)*V(71)

  B(329) = RCT(190)*V(102)

  B(330) = RCT(190)*V(80)

  B(331) = RCT(191)*V(103)

  B(332) = RCT(191)*V(21)

  B(333) = RCT(192)*V(93)

  B(334) = RCT(192)*V(21)

  B(335) = RCT(193)*2*V(21)

  B(336) = RCT(194)*V(102)

  B(337) = RCT(194)*V(20)

  B(338) = RCT(195)*V(102)

  B(339) = RCT(195)*V(16)

  B(340) = RCT(196)*V(102)

  B(341) = RCT(196)*V(17)

  B(342) = RCT(197)*V(103)

  B(343) = RCT(197)*V(15)

  B(344) = RCT(198)*V(93)

  B(345) = RCT(198)*V(15)

  B(346) = RCT(199)

  B(347) = RCT(200)*V(92)

  B(348) = RCT(200)*V(80)

  B(349) = RCT(201)*V(95)

  B(350) = RCT(201)*V(80)

  B(351) = RCT(202)*V(88)

  B(352) = RCT(202)*V(80)

  B(353) = RCT(203)*V(92)

  B(354) = RCT(203)*V(74)

  B(355) = RCT(204)*V(95)

  B(356) = RCT(204)*V(74)

  B(357) = RCT(205)*V(88)

  B(358) = RCT(205)*V(74)

  B(359) = RCT(206)*V(102)

  B(360) = RCT(206)*V(79)

  B(361) = RCT(207)*V(92)

  B(362) = RCT(207)*V(79)

  B(363) = RCT(208)*V(95)

  B(364) = RCT(208)*V(79)

  B(365) = RCT(209)*V(88)

  B(366) = RCT(209)*V(79)

  B(367) = RCT(210)*V(102)

  B(368) = RCT(210)*V(22)

  B(369) = RCT(211)*V(102)

  B(370) = RCT(211)*V(29)

  B(371) = RCT(212)*V(102)

  B(372) = RCT(212)*V(50)

  B(373) = RCT(213)*V(102)

  B(374) = RCT(213)*V(35)

  B(375) = RCT(214)*V(102)

  B(376) = RCT(214)*V(44)

  B(377) = RCT(215)*V(102)

  B(378) = RCT(215)*V(37)

  B(379) = RCT(216)*V(102)

  B(380) = RCT(216)*V(45)

  B(381) = RCT(217)*V(102)

  B(382) = RCT(217)*V(38)

  B(383) = RCT(218)*V(102)

  B(384) = RCT(218)*V(77)

  B(385) = RCT(219)*V(92)

  B(386) = RCT(219)*V(77)

  B(387) = RCT(220)*V(95)

  B(388) = RCT(220)*V(77)

  B(389) = RCT(221)*V(88)

  B(390) = RCT(221)*V(77)

  B(391) = RCT(222)*V(102)

  B(392) = RCT(222)*V(83)

  B(393) = RCT(223)*V(92)

  B(394) = RCT(223)*V(83)

  B(395) = RCT(224)*V(95)

  B(396) = RCT(224)*V(83)

  B(397) = RCT(225)*V(88)

  B(398) = RCT(225)*V(83)

  B(399) = RCT(226)*V(92)

  B(400) = RCT(226)*V(50)

  B(401) = RCT(227)*V(102)

  B(402) = RCT(227)*V(73)

  B(403) = RCT(228)*V(92)

  B(404) = RCT(228)*V(73)

  B(405) = RCT(229)*V(95)

  B(406) = RCT(229)*V(73)

  B(407) = RCT(230)*V(88)

  B(408) = RCT(230)*V(73)

  B(409) = RCT(231)*V(102)

  B(410) = RCT(231)*V(34)

  B(411) = RCT(232)*V(102)

  B(412) = RCT(232)*V(34)

  B(413) = RCT(233)*V(95)

  B(414) = RCT(233)*V(34)

  B(415) = RCT(234)*V(102)

  B(416) = RCT(234)*V(64)

  B(417) = RCT(235)*V(102)

  B(418) = RCT(235)*V(67)

  B(419) = RCT(236)*V(102)

  B(420) = RCT(236)*V(66)

  B(421) = RCT(237)*V(102)

  B(422) = RCT(237)*V(65)

  B(423) = RCT(238)*V(102)

  B(424) = RCT(238)*V(63)

  B(425) = RCT(239)*V(102)

  B(426) = RCT(239)*V(62)

  B(427) = RCT(240)*V(102)

  B(428) = RCT(240)*V(61)

  B(429) = RCT(241)*V(102)

  B(430) = RCT(241)*V(60)

  B(431) = RCT(242)*V(102)

  B(432) = RCT(242)*V(59)

  B(433) = RCT(243)*V(102)

  B(434) = RCT(243)*V(58)

  B(435) = RCT(244)*V(102)

  B(436) = RCT(244)*V(56)

  B(437) = RCT(245)*V(102)

  B(438) = RCT(245)*V(54)

  B(439) = RCT(246)*V(102)

  B(440) = RCT(246)*V(57)

  B(441) = RCT(247)*V(102)

  B(442) = RCT(247)*V(74)



  JVS(1) = 0

  JVS(2) = B(74)

  JVS(3) = B(75)

  JVS(4) = 0

  JVS(5) = B(223)

  JVS(6) = 0.297*B(371)

  JVS(7) = 0.37*B(323)

  JVS(8) = 0.185*B(403)

  JVS(9) = 0.103*B(353)

  JVS(10) = 0.185*B(385)

  JVS(11) = 0.103*B(361)

  JVS(12) = 0.204*B(347)

  JVS(13) = 0.333*B(282)

  JVS(14) = 0.1*B(298)

  JVS(15) = 0.073*B(393)

  JVS(16) = 0.351*B(291)

  JVS(17) = 0.333*B(283)+0.351*B(292)+0.1*B(299)+0.37*B(324)+0.204*B(348)+0.103*B(354)+0.103*B(362)+0.185*B(386)+0.073&
              &*B(394)+0.185*B(404)

  JVS(18) = B(224)

  JVS(19) = 0.297*B(372)

  JVS(20) = 0

  JVS(21) = 0.17*B(403)

  JVS(22) = 0.05*B(385)

  JVS(23) = 0.129*B(393)

  JVS(24) = 0.05*B(386)+0.129*B(394)+0.17*B(404)

  JVS(25) = B(128)

  JVS(26) = 0.25*B(124)+B(129)+B(130)+B(134)

  JVS(27) = B(135)

  JVS(28) = 0.25*B(125)

  JVS(29) = B(131)

  JVS(30) = 0

  JVS(31) = 0.189*B(353)

  JVS(32) = 0.119*B(385)

  JVS(33) = 0.189*B(361)

  JVS(34) = 0.15*B(347)

  JVS(35) = 0.372*B(298)

  JVS(36) = 0.247*B(393)

  JVS(37) = 0.372*B(299)+0.15*B(348)+0.189*B(354)+0.189*B(362)+0.119*B(386)+0.247*B(394)

  JVS(38) = 0.25*B(162)+B(166)+B(168)+B(172)

  JVS(39) = 0.25*B(142)+B(146)+B(148)+B(152)

  JVS(40) = 0.25*B(184)+B(188)+B(190)+2*B(194)

  JVS(41) = B(147)+B(167)+B(189)

  JVS(42) = B(153)+B(173)+2*B(195)

  JVS(43) = 0.25*B(143)+0.25*B(163)+0.25*B(185)

  JVS(44) = B(149)+B(169)+B(191)

  JVS(45) = 0

  JVS(46) = 0.5*B(399)

  JVS(47) = 0.135*B(403)

  JVS(48) = 0.5*B(400)+0.135*B(404)

  JVS(49) = 0

  JVS(50) = 0.75*B(124)

  JVS(51) = 0.75*B(125)

  JVS(52) = 0

  JVS(53) = 0.75*B(162)

  JVS(54) = 0.75*B(142)

  JVS(55) = 0.75*B(184)

  JVS(56) = 0.75*B(143)+0.75*B(163)+0.75*B(185)

  JVS(57) = 0

  JVS(58) = 2*B(211)

  JVS(59) = B(405)

  JVS(60) = 2*B(212)

  JVS(61) = B(406)

  JVS(62) = 0

  JVS(63) = 6*B(211)

  JVS(64) = 7*B(277)

  JVS(65) = 0.048*B(401)+0.07*B(403)+2.693*B(405)+0.55*B(407)

  JVS(66) = 0.55*B(408)

  JVS(67) = 0.07*B(404)

  JVS(68) = 6*B(212)

  JVS(69) = 2.693*B(406)

  JVS(70) = 0.048*B(402)

  JVS(71) = 0

  JVS(72) = B(95)+B(99)

  JVS(73) = B(78)+B(86)+B(96)+B(106)

  JVS(74) = B(82)+B(90)+B(100)+B(112)

  JVS(75) = B(79)+B(83)

  JVS(76) = B(107)+B(113)

  JVS(77) = B(87)+B(91)

  JVS(78) = 0

  JVS(79) = B(97)+B(101)+B(103)+B(105)+B(116)

  JVS(80) = B(80)+B(84)+B(85)+B(92)+B(102)+B(110)

  JVS(81) = B(108)+B(111)+B(114)+B(117)+B(118)

  JVS(82) = B(81)+B(88)+B(98)+B(109)

  JVS(83) = B(89)+B(93)+B(94)+B(104)+B(115)

  JVS(84) = 0

  JVS(85) = B(342)

  JVS(86) = B(343)

  JVS(87) = 0

  JVS(88) = B(346)

  JVS(89) = 0

  JVS(90) = B(344)

  JVS(91) = B(338)

  JVS(92) = B(340)

  JVS(93) = B(335)

  JVS(94) = B(345)

  JVS(95) = B(339)+B(341)

  JVS(96) = -B(342)-B(344)-B(346)

  JVS(97) = 0.12*B(336)

  JVS(98) = -B(345)

  JVS(99) = 0.12*B(337)

  JVS(100) = -B(343)

  JVS(101) = -B(338)

  JVS(102) = 0.75*B(336)

  JVS(103) = 0.75*B(337)-B(339)

  JVS(104) = -B(340)

  JVS(105) = -B(341)

  JVS(106) = -B(32)-B(34)

  JVS(107) = B(31)

  JVS(108) = -B(319)

  JVS(109) = -B(320)

  JVS(110) = -B(336)

  JVS(111) = 0.88*B(331)

  JVS(112) = -B(337)

  JVS(113) = 0.88*B(332)

  JVS(114) = 0.13*B(336)

  JVS(115) = -B(331)-B(333)-2*B(335)

  JVS(116) = B(329)

  JVS(117) = -B(334)

  JVS(118) = B(330)+0.13*B(337)

  JVS(119) = -B(332)

  JVS(120) = -B(367)

  JVS(121) = -B(368)

  JVS(122) = -B(121)

  JVS(123) = B(119)

  JVS(124) = B(120)

  JVS(125) = -B(139)

  JVS(126) = B(137)

  JVS(127) = B(138)

  JVS(128) = -B(159)

  JVS(129) = B(157)

  JVS(130) = B(158)

  JVS(131) = -B(181)

  JVS(132) = B(179)

  JVS(133) = B(180)

  JVS(134) = -B(69)-B(70)

  JVS(135) = -B(71)

  JVS(136) = B(63)+B(64)

  JVS(137) = -B(264)

  JVS(138) = 0.087*B(381)

  JVS(139) = 0.031*B(353)

  JVS(140) = 0.031*B(361)

  JVS(141) = 0.031*B(354)+0.031*B(362)

  JVS(142) = 0.087*B(382)

  JVS(143) = -B(369)

  JVS(144) = -B(370)

  JVS(145) = -B(245)

  JVS(146) = -B(246)

  JVS(147) = -B(23)-B(24)

  JVS(148) = B(21)

  JVS(149) = B(22)

  JVS(150) = -B(74)

  JVS(151) = B(409)+0.5*B(411)+B(413)

  JVS(152) = B(414)

  JVS(153) = -B(75)+B(410)+0.5*B(412)

  JVS(154) = -B(38)-B(39)-B(40)

  JVS(155) = B(36)

  JVS(156) = B(37)-B(41)

  JVS(157) = -B(409)-B(411)-B(413)

  JVS(158) = -B(414)

  JVS(159) = -B(410)-B(412)

  JVS(160) = -B(373)

  JVS(161) = -B(374)

  JVS(162) = 0.236*B(373)

  JVS(163) = -B(203)-B(205)

  JVS(164) = -B(204)

  JVS(165) = 0.236*B(374)

  JVS(166) = -B(377)

  JVS(167) = -B(378)

  JVS(168) = -B(381)

  JVS(169) = -B(382)

  JVS(170) = -B(243)

  JVS(171) = B(84)+0.25*B(92)+0.25*B(110)

  JVS(172) = 0.25*B(111)

  JVS(173) = -B(244)

  JVS(174) = 0.25*B(93)

  JVS(175) = -B(247)-B(249)

  JVS(176) = B(80)

  JVS(177) = -B(248)

  JVS(178) = B(81)

  JVS(179) = -B(222)-B(223)

  JVS(180) = B(220)

  JVS(181) = -B(224)

  JVS(182) = B(221)

  JVS(183) = -B(211)-B(213)-B(215)

  JVS(184) = B(273)

  JVS(185) = -B(212)

  JVS(186) = B(274)

  JVS(187) = -B(214)

  JVS(188) = -B(57)-B(58)-B(59)

  JVS(189) = B(55)

  JVS(190) = -B(60)

  JVS(191) = B(56)

  JVS(192) = -B(375)

  JVS(193) = -B(376)

  JVS(194) = -B(379)

  JVS(195) = -B(380)

  JVS(196) = 0.099*B(381)

  JVS(197) = 0.108*B(379)

  JVS(198) = -B(313)-B(315)

  JVS(199) = -B(314)+0.108*B(380)+0.099*B(382)

  JVS(200) = 0.093*B(381)

  JVS(201) = 0.051*B(379)

  JVS(202) = -B(316)-B(318)

  JVS(203) = -B(317)+0.051*B(380)+0.093*B(382)

  JVS(204) = 0.187*B(381)

  JVS(205) = 0.207*B(379)

  JVS(206) = -B(269)-B(271)

  JVS(207) = -B(272)

  JVS(208) = -B(270)+0.207*B(380)+0.187*B(382)

  JVS(209) = 0.561*B(381)

  JVS(210) = 0.491*B(379)

  JVS(211) = -B(309)-B(311)

  JVS(212) = -B(312)

  JVS(213) = -B(310)+0.491*B(380)+0.561*B(382)

  JVS(214) = -B(371)-B(399)

  JVS(215) = -B(400)

  JVS(216) = -B(372)

  JVS(217) = B(213)+B(215)

  JVS(218) = -B(273)

  JVS(219) = B(206)

  JVS(220) = B(207)

  JVS(221) = -B(274)

  JVS(222) = B(214)

  JVS(223) = 0.05*B(381)

  JVS(224) = 0.059*B(379)

  JVS(225) = -B(275)-B(277)-B(278)

  JVS(226) = 0.061*B(391)+0.042*B(393)+0.015*B(395)

  JVS(227) = 0.042*B(394)

  JVS(228) = -B(279)+0.015*B(396)

  JVS(229) = -B(276)+0.059*B(380)+0.05*B(382)+0.061*B(392)

  JVS(230) = -B(250)-B(252)

  JVS(231) = B(108)

  JVS(232) = -B(251)

  JVS(233) = B(88)+B(109)

  JVS(234) = B(89)

  JVS(235) = -B(437)

  JVS(236) = 0.005*B(439)

  JVS(237) = 0.021*B(433)

  JVS(238) = 0.098*B(431)

  JVS(239) = 0.072*B(429)

  JVS(240) = 0.022*B(427)

  JVS(241) = 0.006*B(425)

  JVS(242) = 0.006*B(426)+0.022*B(428)+0.072*B(430)+0.098*B(432)+0.021*B(434)-B(438)+0.005*B(440)

  JVS(243) = 0.017*B(379)

  JVS(244) = -B(265)-B(267)

  JVS(245) = B(208)+B(210)

  JVS(246) = -B(268)

  JVS(247) = -B(266)+0.017*B(380)

  JVS(248) = B(209)

  JVS(249) = 0.004*B(437)

  JVS(250) = -B(435)

  JVS(251) = 0.005*B(439)

  JVS(252) = 0

  JVS(253) = 0.026*B(431)

  JVS(254) = 0.119*B(429)

  JVS(255) = 0.088*B(427)

  JVS(256) = 0.027*B(425)

  JVS(257) = 0.008*B(423)

  JVS(258) = 0.008*B(424)+0.027*B(426)+0.088*B(428)+0.119*B(430)+0.026*B(432)-B(436)+0.004*B(438)+0.005*B(440)

  JVS(259) = 0.112*B(437)

  JVS(260) = 0.136*B(435)

  JVS(261) = -0.908*B(439)

  JVS(262) = 0.145*B(433)

  JVS(263) = 0.08*B(431)

  JVS(264) = 0.025*B(429)

  JVS(265) = 0.009*B(427)

  JVS(266) = 0.005*B(425)

  JVS(267) = 0.005*B(423)

  JVS(268) = 0.012*B(415)

  JVS(269) = 0.007*B(421)

  JVS(270) = 0.008*B(419)

  JVS(271) = 0.01*B(417)

  JVS(272) = 0.019*B(441)

  JVS(273) = 0.012*B(416)+0.01*B(418)+0.008*B(420)+0.007*B(422)+0.005*B(424)+0.005*B(426)+0.009*B(428)+0.025*B(430)+0.08&
               &*B(432)+0.145*B(434)+0.136*B(436)+0.112*B(438)-0.908*B(440)+0.019*B(442)

  JVS(274) = 0.004*B(437)

  JVS(275) = 0.004*B(435)

  JVS(276) = 0

  JVS(277) = -B(433)

  JVS(278) = 0

  JVS(279) = 0.031*B(429)

  JVS(280) = 0.146*B(427)

  JVS(281) = 0.107*B(425)

  JVS(282) = 0.033*B(423)

  JVS(283) = 0

  JVS(284) = 0.009*B(421)

  JVS(285) = 0

  JVS(286) = 0

  JVS(287) = 0.006*B(441)

  JVS(288) = 0.009*B(422)+0.033*B(424)+0.107*B(426)+0.146*B(428)+0.031*B(430)-B(434)+0.004*B(436)+0.004*B(438)+0.006&
               &*B(442)

  JVS(289) = 0.004*B(435)

  JVS(290) = 0

  JVS(291) = 0.004*B(433)

  JVS(292) = -B(431)

  JVS(293) = 0

  JVS(294) = 0.038*B(427)

  JVS(295) = 0.178*B(425)

  JVS(296) = 0.13*B(423)

  JVS(297) = 0

  JVS(298) = 0.04*B(421)

  JVS(299) = 0.011*B(419)

  JVS(300) = 0

  JVS(301) = 0.028*B(441)

  JVS(302) = 0.011*B(420)+0.04*B(422)+0.13*B(424)+0.178*B(426)+0.038*B(428)-B(432)+0.004*B(434)+0.004*B(436)+0.028&
               &*B(442)

  JVS(303) = 0.004*B(433)

  JVS(304) = 0.004*B(431)

  JVS(305) = -B(429)

  JVS(306) = 0

  JVS(307) = 0.047*B(425)

  JVS(308) = 0.217*B(423)

  JVS(309) = 0

  JVS(310) = 0.159*B(421)

  JVS(311) = 0.049*B(419)

  JVS(312) = 0.014*B(417)

  JVS(313) = 0.114*B(441)

  JVS(314) = 0.014*B(418)+0.049*B(420)+0.159*B(422)+0.217*B(424)+0.047*B(426)-B(430)+0.004*B(432)+0.004*B(434)+0.114&
               &*B(442)

  JVS(315) = 0.004*B(431)

  JVS(316) = 0.004*B(429)

  JVS(317) = -B(427)

  JVS(318) = 0

  JVS(319) = 0.057*B(423)

  JVS(320) = 0.017*B(415)

  JVS(321) = 0.265*B(421)

  JVS(322) = 0.194*B(419)

  JVS(323) = 0.06*B(417)

  JVS(324) = 0.234*B(441)

  JVS(325) = 0.017*B(416)+0.06*B(418)+0.194*B(420)+0.265*B(422)+0.057*B(424)-B(428)+0.004*B(430)+0.004*B(432)+0.234&
               &*B(442)

  JVS(326) = 0.004*B(429)

  JVS(327) = 0.003*B(427)

  JVS(328) = -B(425)

  JVS(329) = 0

  JVS(330) = 0.073*B(415)

  JVS(331) = 0.07*B(421)

  JVS(332) = 0.323*B(419)

  JVS(333) = 0.237*B(417)

  JVS(334) = 0.169*B(441)

  JVS(335) = 0.073*B(416)+0.237*B(418)+0.323*B(420)+0.07*B(422)-B(426)+0.003*B(428)+0.004*B(430)+0.169*B(442)

  JVS(336) = 0.003*B(427)

  JVS(337) = 0.003*B(425)

  JVS(338) = -B(423)

  JVS(339) = 0.289*B(415)

  JVS(340) = 0

  JVS(341) = 0.085*B(419)

  JVS(342) = 0.394*B(417)

  JVS(343) = 0.034*B(441)

  JVS(344) = 0.289*B(416)+0.394*B(418)+0.085*B(420)-B(424)+0.003*B(426)+0.003*B(428)+0.034*B(442)

  JVS(345) = -B(415)

  JVS(346) = 0.002*B(419)

  JVS(347) = 0.001*B(417)

  JVS(348) = 0.001*B(441)

  JVS(349) = -B(416)+0.001*B(418)+0.002*B(420)+0.001*B(442)

  JVS(350) = 0.003*B(425)

  JVS(351) = 0.003*B(423)

  JVS(352) = 0.481*B(415)

  JVS(353) = -B(421)

  JVS(354) = 0

  JVS(355) = 0.104*B(417)

  JVS(356) = 0

  JVS(357) = 0.481*B(416)+0.104*B(418)-B(422)+0.003*B(424)+0.003*B(426)

  JVS(358) = 0.003*B(423)

  JVS(359) = 0.127*B(415)

  JVS(360) = 0.002*B(421)

  JVS(361) = -B(419)

  JVS(362) = 0

  JVS(363) = 0.001*B(441)

  JVS(364) = 0.127*B(416)-B(420)+0.002*B(422)+0.003*B(424)+0.001*B(442)

  JVS(365) = 0.002*B(421)

  JVS(366) = 0.002*B(419)

  JVS(367) = -B(417)

  JVS(368) = 0.002*B(441)

  JVS(369) = -B(418)+0.002*B(420)+0.002*B(422)+0.002*B(442)

  JVS(370) = 0.287*B(381)

  JVS(371) = 0.119*B(379)

  JVS(372) = 0.5*B(315)

  JVS(373) = 0.5*B(318)

  JVS(374) = 0.23*B(269)

  JVS(375) = -B(259)-B(260)-B(262)

  JVS(376) = 0.084*B(280)+0.9*B(282)

  JVS(377) = 0.174*B(296)+0.742*B(298)+0.008*B(300)

  JVS(378) = 0.3*B(289)+0.95*B(291)

  JVS(379) = 0.9*B(283)+0.95*B(292)+0.742*B(299)

  JVS(380) = -B(263)+0.008*B(301)

  JVS(381) = -B(261)+0.23*B(270)+0.084*B(281)+0.3*B(290)+0.174*B(297)+0.119*B(380)+0.287*B(382)

  JVS(382) = 0.002*B(375)

  JVS(383) = B(315)

  JVS(384) = B(318)

  JVS(385) = B(309)+1.5*B(311)

  JVS(386) = 0.393*B(371)+1.5*B(399)

  JVS(387) = B(259)+B(260)+B(262)

  JVS(388) = -B(49)

  JVS(389) = 0.5*B(323)+0.491*B(327)

  JVS(390) = 0.51*B(403)

  JVS(391) = 0.157*B(353)

  JVS(392) = 2*B(253)+B(254)+1.26*B(255)+1.26*B(257)

  JVS(393) = 0.345*B(385)

  JVS(394) = 0.157*B(361)

  JVS(395) = 0.275*B(347)

  JVS(396) = 0.416*B(280)+0.45*B(282)+0.5*B(284)+0.67*B(288)

  JVS(397) = 0.336*B(296)+0.498*B(298)+0.572*B(300)+1.233*B(302)

  JVS(398) = 0.265*B(393)+0.012*B(397)

  JVS(399) = 0.475*B(291)+0.7*B(295)

  JVS(400) = B(229)

  JVS(401) = B(216)+B(217)+B(218)+B(225)

  JVS(402) = 0.491*B(328)+0.012*B(398)

  JVS(403) = 0.034*B(232)+B(234)

  JVS(404) = 0.45*B(283)+0.475*B(292)+0.498*B(299)+1.5*B(312)+0.5*B(324)+0.275*B(348)+0.157*B(354)+0.157*B(362)+0.345&
               &*B(386)+0.265*B(394)+1.5*B(400)+0.51*B(404)

  JVS(405) = B(226)+1.26*B(258)+B(263)+0.5*B(285)+0.572*B(301)

  JVS(406) = -B(50)+B(219)+0.034*B(233)+1.26*B(256)+B(261)+0.416*B(281)+0.336*B(297)+B(310)+0.393*B(372)+0.002*B(376)

  JVS(407) = 2*B(24)

  JVS(408) = B(413)

  JVS(409) = B(271)

  JVS(410) = B(273)

  JVS(411) = B(278)

  JVS(412) = B(267)

  JVS(413) = B(262)

  JVS(414) = -B(46)-B(48)

  JVS(415) = B(257)

  JVS(416) = 0

  JVS(417) = 0.5*B(284)

  JVS(418) = 0.15*B(300)

  JVS(419) = 0

  JVS(420) = 0

  JVS(421) = B(230)

  JVS(422) = B(225)

  JVS(423) = B(235)

  JVS(424) = 0

  JVS(425) = B(42)

  JVS(426) = 0.2*B(66)+B(226)+B(231)+B(236)+B(258)+B(263)+B(268)+B(272)+B(274)+B(279)+0.5*B(285)+0.15*B(301)+B(414)

  JVS(427) = B(43)-B(47)

  JVS(428) = 0.2*B(67)

  JVS(429) = -B(321)-B(323)-B(325)-B(327)

  JVS(430) = -B(328)

  JVS(431) = -B(324)

  JVS(432) = -B(326)

  JVS(433) = -B(322)

  JVS(434) = 0.704*B(369)

  JVS(435) = 0.024*B(373)

  JVS(436) = B(205)

  JVS(437) = 0.072*B(377)

  JVS(438) = 0.452*B(375)

  JVS(439) = -B(237)-B(239)

  JVS(440) = 0.13*B(353)

  JVS(441) = 0.005*B(383)+0.001*B(385)+0.024*B(387)

  JVS(442) = 0.13*B(361)

  JVS(443) = 0.127*B(391)+0.045*B(393)+0.102*B(395)

  JVS(444) = 0.006*B(306)+0.02*B(308)

  JVS(445) = 0.13*B(354)+0.13*B(362)+0.001*B(386)+0.045*B(394)

  JVS(446) = 0

  JVS(447) = 0.024*B(388)+0.102*B(396)

  JVS(448) = -B(238)+0.006*B(307)+0.704*B(370)+0.024*B(374)+0.452*B(376)+0.072*B(378)+0.005*B(384)+0.127*B(392)

  JVS(449) = -B(401)-B(403)-B(405)-B(407)

  JVS(450) = -B(408)

  JVS(451) = -B(404)

  JVS(452) = -B(406)

  JVS(453) = -B(402)

  JVS(454) = -B(353)-B(355)-B(357)

  JVS(455) = -B(358)

  JVS(456) = -B(354)

  JVS(457) = -B(356)

  JVS(458) = 0.097*B(381)

  JVS(459) = 0.118*B(379)

  JVS(460) = 0.5*B(315)

  JVS(461) = 0.5*B(318)

  JVS(462) = B(311)

  JVS(463) = 0.607*B(371)

  JVS(464) = 0.23*B(265)

  JVS(465) = 0.009*B(327)

  JVS(466) = 0.001*B(353)

  JVS(467) = -B(253)-B(254)-B(255)-B(257)

  JVS(468) = 0

  JVS(469) = 0.001*B(361)

  JVS(470) = 0.15*B(296)+0.023*B(298)

  JVS(471) = 0.009*B(328)

  JVS(472) = 0.023*B(299)+B(312)+0.001*B(354)+0.001*B(362)

  JVS(473) = -B(258)

  JVS(474) = -B(256)+0.23*B(266)+0.15*B(297)+0.607*B(372)+0.118*B(380)+0.097*B(382)

  JVS(475) = 0

  JVS(476) = 0.24*B(269)+B(271)

  JVS(477) = 0.24*B(265)+B(267)

  JVS(478) = -B(206)-B(208)-B(210)

  JVS(479) = B(160)

  JVS(480) = -B(207)

  JVS(481) = B(164)+B(268)+B(272)

  JVS(482) = B(161)+B(165)+B(174)+B(176)+2*B(178)+B(200)

  JVS(483) = B(177)

  JVS(484) = B(201)

  JVS(485) = B(175)

  JVS(486) = 0.24*B(266)+0.24*B(270)

  JVS(487) = -B(209)

  JVS(488) = -B(383)-B(385)-B(387)-B(389)

  JVS(489) = -B(390)

  JVS(490) = -B(386)

  JVS(491) = -B(388)

  JVS(492) = -B(384)

  JVS(493) = 0.559*B(373)

  JVS(494) = 0.948*B(377)

  JVS(495) = 0.936*B(375)

  JVS(496) = B(313)+B(315)

  JVS(497) = B(316)+B(318)

  JVS(498) = B(237)

  JVS(499) = 0.729*B(353)+0.75*B(355)

  JVS(500) = 0.205*B(383)+0.488*B(387)

  JVS(501) = -B(95)-B(97)-B(99)-B(101)-B(103)-B(116)-B(132)-B(150)-B(170)-B(192)

  JVS(502) = 0.5*B(359)+0.729*B(361)+0.75*B(363)

  JVS(503) = 0.079*B(329)+0.126*B(347)+0.187*B(349)+0.24*B(351)

  JVS(504) = 0.001*B(391)+0.137*B(393)+0.711*B(395)

  JVS(505) = 0.675*B(289)

  JVS(506) = 0.596*B(306)+0.152*B(308)

  JVS(507) = 0.24*B(352)

  JVS(508) = 0.616*B(240)

  JVS(509) = 0.515*B(305)

  JVS(510) = 0.126*B(348)+0.729*B(354)+0.729*B(362)+0.137*B(394)

  JVS(511) = -B(96)+B(160)

  JVS(512) = 0

  JVS(513) = -B(100)+B(164)+0.187*B(350)+0.75*B(356)+0.75*B(364)+0.488*B(388)+0.711*B(396)

  JVS(514) = B(161)+B(165)-B(171)+B(174)+B(176)+2*B(178)+B(200)

  JVS(515) = -B(151)+B(177)

  JVS(516) = -B(193)+B(201)

  JVS(517) = -B(102)

  JVS(518) = -B(133)+B(175)

  JVS(519) = -B(117)

  JVS(520) = B(238)+0.616*B(241)+0.675*B(290)+0.596*B(307)+B(314)+B(317)+0.079*B(330)+0.5*B(360)+0.559*B(374)+0.936&
               &*B(376)+0.948*B(378)+0.205*B(384)+0.001*B(392)

  JVS(521) = -B(98)

  JVS(522) = -B(104)

  JVS(523) = -B(359)-B(361)-B(363)-B(365)

  JVS(524) = -B(366)

  JVS(525) = -B(362)

  JVS(526) = -B(364)

  JVS(527) = -B(360)

  JVS(528) = -B(329)-B(347)-B(349)-B(351)

  JVS(529) = -B(352)

  JVS(530) = -B(348)

  JVS(531) = -B(350)

  JVS(532) = -B(330)

  JVS(533) = 0.23*B(329)+0.39*B(347)

  JVS(534) = -B(280)-B(282)-B(284)-B(286)-B(288)

  JVS(535) = 0.025*B(391)+0.026*B(393)+0.012*B(397)

  JVS(536) = -B(287)+0.012*B(398)

  JVS(537) = -B(283)+0.39*B(348)+0.026*B(394)

  JVS(538) = -B(285)

  JVS(539) = -B(281)+0.23*B(330)+0.025*B(392)

  JVS(540) = 0.357*B(329)+0.936*B(349)

  JVS(541) = -B(296)-B(298)-B(300)-B(302)

  JVS(542) = 0.025*B(391)

  JVS(543) = 0

  JVS(544) = -B(299)

  JVS(545) = -B(301)+0.936*B(350)

  JVS(546) = -B(297)+0.357*B(330)+0.025*B(392)

  JVS(547) = -B(391)-B(393)-B(395)-B(397)

  JVS(548) = -B(398)

  JVS(549) = -B(394)

  JVS(550) = -B(396)

  JVS(551) = -B(392)

  JVS(552) = 0.32*B(329)+0.16*B(347)

  JVS(553) = 0.019*B(393)+0.048*B(395)

  JVS(554) = -B(289)-B(291)-B(293)-B(295)

  JVS(555) = -B(294)

  JVS(556) = -B(292)+0.16*B(348)+0.019*B(394)

  JVS(557) = 0.048*B(396)

  JVS(558) = -B(290)+0.32*B(330)

  JVS(559) = B(367)

  JVS(560) = 0.96*B(245)

  JVS(561) = 0.445*B(373)

  JVS(562) = 0.099*B(377)

  JVS(563) = 0.455*B(375)

  JVS(564) = 0.195*B(321)+0.25*B(327)

  JVS(565) = 0.984*B(401)+0.5*B(403)

  JVS(566) = 0.294*B(383)+0.154*B(385)+0.009*B(387)

  JVS(567) = 0.129*B(296)+0.047*B(298)+0.467*B(302)

  JVS(568) = 0.732*B(391)+0.456*B(393)+0.507*B(395)

  JVS(569) = -B(227)-B(229)-B(230)

  JVS(570) = 0.439*B(306)+0.431*B(308)

  JVS(571) = 0.25*B(328)

  JVS(572) = 0.034*B(232)+B(234)

  JVS(573) = 0.482*B(240)+B(242)

  JVS(574) = 0.084*B(303)+0.246*B(305)

  JVS(575) = 0.047*B(299)+0.154*B(386)+0.456*B(394)+0.5*B(404)

  JVS(576) = B(140)

  JVS(577) = B(144)-B(231)+0.009*B(388)+0.507*B(396)

  JVS(578) = B(176)

  JVS(579) = B(141)+B(145)+B(154)+2*B(156)+B(177)+B(198)

  JVS(580) = B(199)

  JVS(581) = B(155)

  JVS(582) = -B(228)+0.034*B(233)+0.482*B(241)+0.96*B(246)+0.129*B(297)+0.084*B(304)+0.439*B(307)+0.195*B(322)+B(368)&
               &+0.445*B(374)+0.455*B(376)+0.099*B(378)+0.294*B(384)+0.732*B(392)+0.984*B(402)

  JVS(583) = 0.081*B(245)

  JVS(584) = 0.026*B(373)

  JVS(585) = 0.026*B(377)

  JVS(586) = B(243)

  JVS(587) = 0.35*B(247)+B(249)

  JVS(588) = B(222)

  JVS(589) = 0.024*B(375)

  JVS(590) = 0.096*B(371)

  JVS(591) = 1.61*B(321)+B(323)+0.191*B(327)

  JVS(592) = B(237)

  JVS(593) = 0.984*B(401)+0.5*B(403)

  JVS(594) = 0.235*B(353)

  JVS(595) = B(254)

  JVS(596) = 0

  JVS(597) = 0.732*B(383)+0.5*B(385)

  JVS(598) = 0.276*B(359)+0.235*B(361)

  JVS(599) = 0.624*B(329)+0.592*B(347)+0.24*B(351)

  JVS(600) = 0.084*B(280)+0.2*B(282)+0.67*B(288)

  JVS(601) = 0.055*B(296)+0.125*B(298)+0.227*B(300)+0.3*B(302)

  JVS(602) = 0.244*B(391)+0.269*B(393)+0.079*B(395)

  JVS(603) = 0.3*B(289)+0.1*B(291)

  JVS(604) = -B(216)-B(217)-B(218)-B(220)-B(225)

  JVS(605) = 0.01*B(306)+0.134*B(308)

  JVS(606) = 0.191*B(328)+0.24*B(352)

  JVS(607) = 0.115*B(240)

  JVS(608) = 0.213*B(303)+0.506*B(305)

  JVS(609) = 0.2*B(283)+0.1*B(292)+0.125*B(299)+B(324)+0.592*B(348)+0.235*B(354)+0.235*B(362)+0.5*B(386)+0.269*B(394)&
               &+0.5*B(404)

  JVS(610) = B(78)+B(182)

  JVS(611) = 0

  JVS(612) = B(82)+B(186)-B(226)+0.227*B(301)+0.079*B(396)

  JVS(613) = B(166)+B(200)

  JVS(614) = B(146)+B(198)

  JVS(615) = B(183)+B(187)+B(188)+B(196)+B(199)+B(201)+2*B(202)

  JVS(616) = B(79)+B(83)+B(84)+2*B(85)+0.75*B(92)+0.75*B(110)+B(128)+B(147)+B(167)+B(189)

  JVS(617) = B(129)+B(197)

  JVS(618) = 0.75*B(111)

  JVS(619) = -B(219)+B(238)+0.115*B(241)+B(244)+0.081*B(246)+0.35*B(248)+0.084*B(281)+0.3*B(290)+0.055*B(297)+0.213&
               &*B(304)+0.01*B(307)+1.61*B(322)+0.624*B(330)+0.276*B(360)+0.096*B(372)+0.026*B(374)+0.024*B(376)+0.026&
               &*B(378)+0.732*B(384)+0.244*B(392)+0.984*B(402)

  JVS(620) = -B(221)

  JVS(621) = 0.75*B(93)

  JVS(622) = B(203)

  JVS(623) = 0.276*B(355)

  JVS(624) = 0.511*B(387)

  JVS(625) = 0.276*B(363)

  JVS(626) = 0.572*B(300)

  JVS(627) = 0.321*B(395)

  JVS(628) = -0.69*B(306)-B(308)

  JVS(629) = 0

  JVS(630) = 0

  JVS(631) = B(106)

  JVS(632) = B(204)

  JVS(633) = 0.572*B(301)+0.276*B(356)+0.276*B(364)+0.511*B(388)+0.321*B(396)

  JVS(634) = B(107)

  JVS(635) = -0.69*B(307)

  JVS(636) = B(34)

  JVS(637) = -B(327)

  JVS(638) = -B(407)

  JVS(639) = -B(357)

  JVS(640) = -B(389)

  JVS(641) = -B(365)

  JVS(642) = -B(351)

  JVS(643) = -B(286)

  JVS(644) = -B(397)

  JVS(645) = -B(293)

  JVS(646) = -B(2)-B(4)-B(6)-B(9)-B(11)-B(287)-B(294)-B(328)-B(352)-B(358)-B(366)-B(390)-B(398)-B(408)

  JVS(647) = -B(5)+B(30)

  JVS(648) = -B(7)

  JVS(649) = B(1)-B(10)-B(12)

  JVS(650) = B(29)

  JVS(651) = 0

  JVS(652) = 0.261*B(369)

  JVS(653) = 0.122*B(373)

  JVS(654) = 0.204*B(377)

  JVS(655) = 0.244*B(375)

  JVS(656) = B(313)

  JVS(657) = B(316)

  JVS(658) = B(309)

  JVS(659) = B(250)+B(252)

  JVS(660) = B(325)

  JVS(661) = 0.45*B(407)

  JVS(662) = 0.205*B(353)+0.474*B(355)+0.147*B(357)

  JVS(663) = 0.497*B(383)+0.363*B(385)+0.037*B(387)+0.45*B(389)

  JVS(664) = 0.474*B(359)+0.205*B(361)+0.474*B(363)+0.147*B(365)

  JVS(665) = B(286)

  JVS(666) = 0.013*B(296)+0.218*B(300)

  JVS(667) = 0.511*B(391)+0.305*B(393)+0.151*B(395)+0.069*B(397)

  JVS(668) = 0.675*B(289)+0.45*B(293)

  JVS(669) = 0.213*B(306)+0.147*B(308)

  JVS(670) = B(287)+0.45*B(294)+0.147*B(358)+0.147*B(366)+0.45*B(390)+0.069*B(398)+0.45*B(408)

  JVS(671) = -B(232)-B(234)-B(235)

  JVS(672) = 0.37*B(240)

  JVS(673) = 0.558*B(303)+0.71*B(305)

  JVS(674) = 0.205*B(354)+0.205*B(362)+0.363*B(386)+0.305*B(394)

  JVS(675) = 0

  JVS(676) = 0

  JVS(677) = -B(236)+0.218*B(301)+B(326)+0.474*B(356)+0.474*B(364)+0.037*B(388)+0.151*B(396)

  JVS(678) = 0

  JVS(679) = -B(233)+0.37*B(241)+B(251)+0.675*B(290)+0.013*B(297)+0.558*B(304)+0.213*B(307)+B(310)+B(314)+B(317)+0.474&
               &*B(360)+0.261*B(370)+0.122*B(374)+0.244*B(376)+0.204*B(378)+0.497*B(384)+0.511*B(392)

  JVS(680) = 0

  JVS(681) = 0

  JVS(682) = 0.332*B(373)

  JVS(683) = 0.089*B(377)

  JVS(684) = 0.11*B(375)

  JVS(685) = 0.55*B(407)

  JVS(686) = 0.437*B(389)

  JVS(687) = 0.416*B(280)

  JVS(688) = 0.15*B(296)+0.21*B(298)+0.233*B(302)

  JVS(689) = 0.072*B(391)+0.026*B(393)+0.001*B(395)+0.659*B(397)

  JVS(690) = 0.55*B(293)

  JVS(691) = 0.177*B(306)+0.243*B(308)

  JVS(692) = 0.55*B(294)+0.437*B(390)+0.659*B(398)+0.55*B(408)

  JVS(693) = -B(240)-B(242)

  JVS(694) = 0.115*B(303)

  JVS(695) = 0.21*B(299)+0.026*B(394)

  JVS(696) = 0

  JVS(697) = 0

  JVS(698) = B(112)+0.001*B(396)

  JVS(699) = 0.5*B(110)

  JVS(700) = 0.5*B(111)+B(113)+0.5*B(114)+B(118)

  JVS(701) = -B(241)+0.416*B(281)+0.15*B(297)+0.115*B(304)+0.177*B(307)+0.332*B(374)+0.11*B(376)+0.089*B(378)+0.072&
               &*B(392)

  JVS(702) = 0.5*B(115)

  JVS(703) = 0.417*B(377)

  JVS(704) = 0.125*B(375)

  JVS(705) = 0.055*B(379)

  JVS(706) = 0.276*B(353)+0.853*B(357)

  JVS(707) = 0.119*B(383)+0.215*B(385)+0.113*B(389)

  JVS(708) = 0.276*B(359)+0.276*B(361)+0.853*B(365)

  JVS(709) = 0.1*B(347)+0.75*B(351)

  JVS(710) = 0.332*B(296)

  JVS(711) = 0.043*B(393)+0.259*B(397)

  JVS(712) = 0.7*B(295)

  JVS(713) = 0.048*B(306)+0.435*B(308)

  JVS(714) = 0.75*B(352)+0.853*B(358)+0.853*B(366)+0.113*B(390)+0.259*B(398)

  JVS(715) = -0.671*B(303)-B(305)

  JVS(716) = 0.1*B(348)+0.276*B(354)+0.276*B(362)+0.215*B(386)+0.043*B(394)

  JVS(717) = 0

  JVS(718) = 0

  JVS(719) = 0

  JVS(720) = B(172)

  JVS(721) = B(152)

  JVS(722) = 0.5*B(110)

  JVS(723) = B(134)

  JVS(724) = 0.5*B(111)+0.5*B(114)+B(118)+B(135)+B(153)+B(173)

  JVS(725) = 0.332*B(297)-0.671*B(304)+0.048*B(307)+0.276*B(360)+0.125*B(376)+0.417*B(378)+0.055*B(380)+0.119*B(384)

  JVS(726) = 0.5*B(115)

  JVS(727) = -B(311)

  JVS(728) = -B(399)

  JVS(729) = -B(323)

  JVS(730) = -B(403)

  JVS(731) = -B(353)

  JVS(732) = -B(385)

  JVS(733) = -B(361)

  JVS(734) = -B(347)

  JVS(735) = -B(282)

  JVS(736) = -B(298)

  JVS(737) = -B(393)

  JVS(738) = -B(291)

  JVS(739) = B(2)-B(4)

  JVS(740) = -B(5)-B(13)-B(15)-B(30)-B(31)-B(51)-B(61)-B(283)-B(292)-B(299)-B(312)-B(324)-B(348)-B(354)-B(362)-B(386)&
               &-B(394)-B(400)-B(404)

  JVS(741) = -B(14)

  JVS(742) = -B(16)

  JVS(743) = 0

  JVS(744) = 0.25*B(162)

  JVS(745) = 0.25*B(142)

  JVS(746) = 0.25*B(184)

  JVS(747) = 0.25*B(124)

  JVS(748) = -B(52)

  JVS(749) = -B(62)+0.25*B(125)+0.25*B(143)+0.25*B(163)+0.25*B(185)

  JVS(750) = B(38)

  JVS(751) = -B(223)

  JVS(752) = -B(95)

  JVS(753) = 0

  JVS(754) = 0

  JVS(755) = 0

  JVS(756) = 0

  JVS(757) = 0

  JVS(758) = 0

  JVS(759) = -B(6)+B(9)

  JVS(760) = 0

  JVS(761) = 0

  JVS(762) = -B(13)

  JVS(763) = -B(7)-B(14)-B(17)-2*B(19)-B(36)-B(53)-B(78)-B(86)-B(96)-B(106)-B(122)-B(140)-B(160)-B(182)-B(224)

  JVS(764) = B(1)+B(10)+B(26)

  JVS(765) = -B(18)+B(27)+B(28)

  JVS(766) = -B(161)

  JVS(767) = -B(141)

  JVS(768) = -B(183)

  JVS(769) = -B(79)

  JVS(770) = -B(123)

  JVS(771) = -B(107)

  JVS(772) = -B(37)

  JVS(773) = -B(54)

  JVS(774) = -B(87)

  JVS(775) = B(121)

  JVS(776) = B(139)

  JVS(777) = B(159)

  JVS(778) = B(181)

  JVS(779) = B(23)

  JVS(780) = B(39)+B(40)

  JVS(781) = -B(203)

  JVS(782) = B(223)

  JVS(783) = -B(211)

  JVS(784) = B(57)+0.61*B(58)+B(59)

  JVS(785) = 0

  JVS(786) = B(48)

  JVS(787) = 0.474*B(355)

  JVS(788) = 0

  JVS(789) = -B(206)

  JVS(790) = B(95)+B(99)

  JVS(791) = 0.474*B(363)

  JVS(792) = 0.187*B(349)

  JVS(793) = 0

  JVS(794) = 0

  JVS(795) = 0.391*B(395)

  JVS(796) = 0

  JVS(797) = 0

  JVS(798) = 0

  JVS(799) = 0.338*B(306)+B(308)

  JVS(800) = B(6)-B(9)-B(11)

  JVS(801) = 0

  JVS(802) = 0

  JVS(803) = 0

  JVS(804) = B(13)-B(15)

  JVS(805) = B(7)+B(14)+2*B(17)+2*B(19)+B(53)+B(78)+B(86)+B(96)+B(122)+B(140)+B(160)+B(182)+B(224)

  JVS(806) = -B(1)-B(10)-B(12)-B(16)-B(21)-B(42)-B(55)-B(119)-B(137)-B(157)-B(179)-B(204)-B(207)-B(212)

  JVS(807) = 2*B(18)-B(22)+B(29)+B(44)+0.8*B(66)+2*B(68)+B(82)+B(90)+B(100)+B(112)+B(126)+B(144)+B(164)+B(186)+0.187&
               &*B(350)+0.474*B(356)+0.474*B(364)+0.391*B(396)

  JVS(808) = -B(158)+B(161)+B(165)

  JVS(809) = -B(138)+B(141)+B(145)

  JVS(810) = -B(180)+B(183)+B(187)

  JVS(811) = B(79)+B(83)

  JVS(812) = -B(120)+B(123)+B(127)

  JVS(813) = B(113)

  JVS(814) = B(41)-B(43)+B(45)+B(60)+0.338*B(307)

  JVS(815) = B(54)-B(56)+0.8*B(67)

  JVS(816) = B(87)+B(91)

  JVS(817) = B(23)

  JVS(818) = -B(413)

  JVS(819) = 0.39*B(58)

  JVS(820) = -B(271)

  JVS(821) = -B(273)

  JVS(822) = -B(278)

  JVS(823) = -B(267)

  JVS(824) = -B(262)

  JVS(825) = B(46)

  JVS(826) = -B(325)

  JVS(827) = -B(405)

  JVS(828) = -B(355)

  JVS(829) = -B(257)

  JVS(830) = 0

  JVS(831) = -B(387)

  JVS(832) = -B(99)

  JVS(833) = -B(363)

  JVS(834) = -B(349)

  JVS(835) = -B(284)

  JVS(836) = -B(300)

  JVS(837) = -B(395)

  JVS(838) = 0

  JVS(839) = -B(230)

  JVS(840) = -B(225)

  JVS(841) = 0

  JVS(842) = B(11)

  JVS(843) = -B(235)

  JVS(844) = 0

  JVS(845) = 0

  JVS(846) = B(15)

  JVS(847) = -B(17)

  JVS(848) = B(12)+B(16)-B(21)-B(26)

  JVS(849) = -B(18)-B(22)-B(27)-B(28)-B(29)-B(44)-B(66)-2*B(68)-B(82)-B(90)-B(100)-B(112)-B(126)-B(144)-B(164)-B(186)&
               &-B(226)-B(231)-B(236)-B(258)-B(263)-B(268)-B(272)-B(274)-B(279)-B(285)-B(301)-B(326)-B(350)-B(356)-B(364)&
               &-B(388)-B(396)-B(406)-B(414)

  JVS(850) = -B(165)

  JVS(851) = -B(145)

  JVS(852) = -B(187)

  JVS(853) = -B(83)

  JVS(854) = -B(127)

  JVS(855) = -B(113)

  JVS(856) = -B(45)+B(47)

  JVS(857) = -B(67)

  JVS(858) = -B(91)

  JVS(859) = B(159)

  JVS(860) = B(275)+B(278)

  JVS(861) = 0

  JVS(862) = 0

  JVS(863) = 0

  JVS(864) = -B(160)

  JVS(865) = -B(157)

  JVS(866) = -B(164)+B(279)

  JVS(867) = -B(158)-B(161)-B(162)-B(165)-B(166)-B(168)-B(172)-B(174)-B(176)-2*B(178)-B(200)

  JVS(868) = -B(177)

  JVS(869) = -B(201)

  JVS(870) = -B(167)

  JVS(871) = -B(175)

  JVS(872) = -B(173)

  JVS(873) = B(276)

  JVS(874) = -B(163)

  JVS(875) = -B(169)

  JVS(876) = B(139)

  JVS(877) = 0.201*B(353)

  JVS(878) = 0.37*B(255)+0.37*B(257)

  JVS(879) = 0

  JVS(880) = 0.201*B(361)

  JVS(881) = 0.1*B(282)

  JVS(882) = 0.048*B(298)+0.3*B(302)

  JVS(883) = 0.006*B(393)

  JVS(884) = 0.05*B(291)

  JVS(885) = 0

  JVS(886) = 0.965*B(232)+B(235)

  JVS(887) = 0.096*B(240)

  JVS(888) = 0.049*B(303)+0.333*B(305)

  JVS(889) = 0.1*B(283)+0.05*B(292)+0.048*B(299)+0.201*B(354)+0.201*B(362)+0.006*B(394)

  JVS(890) = -B(140)

  JVS(891) = -B(137)

  JVS(892) = -B(144)+B(236)+0.37*B(258)

  JVS(893) = -B(176)

  JVS(894) = -B(138)-B(141)-B(142)-B(145)-B(146)-B(148)-B(152)-B(154)-2*B(156)-B(177)-B(198)

  JVS(895) = -B(199)

  JVS(896) = -B(147)

  JVS(897) = -B(155)

  JVS(898) = -B(153)

  JVS(899) = 0.965*B(233)+0.096*B(241)+0.37*B(256)+0.049*B(304)

  JVS(900) = -B(143)

  JVS(901) = -B(149)

  JVS(902) = B(181)

  JVS(903) = 0.192*B(347)+0.24*B(351)

  JVS(904) = 0.5*B(280)+0.5*B(284)+0.33*B(288)

  JVS(905) = 0.289*B(296)+0.15*B(300)

  JVS(906) = 0

  JVS(907) = 0.3*B(295)

  JVS(908) = 0.24*B(352)

  JVS(909) = 0.192*B(348)

  JVS(910) = -B(182)

  JVS(911) = -B(179)

  JVS(912) = -B(186)+0.5*B(285)+0.15*B(301)

  JVS(913) = -B(200)

  JVS(914) = -B(198)

  JVS(915) = -B(180)-B(183)-B(184)-B(187)-B(188)-B(190)-B(194)-B(196)-B(199)-B(201)-2*B(202)

  JVS(916) = -B(189)

  JVS(917) = -B(197)

  JVS(918) = -B(195)

  JVS(919) = 0.5*B(281)+0.289*B(297)

  JVS(920) = -B(185)

  JVS(921) = -B(191)

  JVS(922) = B(319)

  JVS(923) = B(205)

  JVS(924) = 0.65*B(247)

  JVS(925) = 0.011*B(375)

  JVS(926) = 0.3*B(327)

  JVS(927) = B(239)

  JVS(928) = 0.26*B(403)

  JVS(929) = 0

  JVS(930) = 0.076*B(385)

  JVS(931) = 0

  JVS(932) = 0.25*B(351)

  JVS(933) = 0.197*B(393)+0.03*B(395)

  JVS(934) = 0.3*B(295)

  JVS(935) = B(229)

  JVS(936) = 0

  JVS(937) = 0.3*B(328)+0.25*B(352)

  JVS(938) = 0

  JVS(939) = 0

  JVS(940) = 0

  JVS(941) = 0.076*B(386)+0.197*B(394)+0.26*B(404)

  JVS(942) = -B(78)+B(122)

  JVS(943) = 0

  JVS(944) = -B(82)+B(126)+0.03*B(396)

  JVS(945) = -B(166)+B(174)

  JVS(946) = -B(146)+B(154)

  JVS(947) = -B(188)+B(196)

  JVS(948) = -B(79)-B(80)-B(83)-2*B(84)-2*B(85)-B(92)-B(110)-B(128)-B(147)-B(167)-B(189)

  JVS(949) = B(123)+B(127)-B(129)+2*B(136)+B(155)+B(175)+B(197)

  JVS(950) = -B(111)

  JVS(951) = 0.65*B(248)+B(320)+0.011*B(376)

  JVS(952) = -B(81)

  JVS(953) = -B(93)

  JVS(954) = B(121)

  JVS(955) = 2*B(264)

  JVS(956) = 0

  JVS(957) = 0.011*B(375)

  JVS(958) = B(313)+0.5*B(315)

  JVS(959) = B(316)+0.5*B(318)

  JVS(960) = B(259)+B(260)+B(262)

  JVS(961) = B(237)+B(239)

  JVS(962) = 0.123*B(353)

  JVS(963) = 0

  JVS(964) = 0.123*B(361)

  JVS(965) = 0.67*B(288)

  JVS(966) = 0.467*B(302)

  JVS(967) = 0.137*B(393)

  JVS(968) = 0.675*B(289)

  JVS(969) = B(227)+B(230)

  JVS(970) = 0

  JVS(971) = 0

  JVS(972) = 0

  JVS(973) = 0.492*B(240)+B(242)

  JVS(974) = 0.029*B(303)+0.667*B(305)

  JVS(975) = 0.123*B(354)+0.123*B(362)+0.137*B(394)

  JVS(976) = -B(122)+B(182)

  JVS(977) = -B(119)

  JVS(978) = -B(126)+B(186)+B(231)+B(263)

  JVS(979) = -B(174)+B(200)

  JVS(980) = -B(154)+B(198)

  JVS(981) = B(183)+B(187)+B(199)+B(201)+2*B(202)

  JVS(982) = -B(128)

  JVS(983) = -B(120)-B(123)-B(124)-B(127)-B(129)-B(130)-B(134)-2*B(136)-B(155)-B(175)

  JVS(984) = -B(135)

  JVS(985) = B(228)+B(238)+0.492*B(241)+B(261)+0.675*B(290)+0.029*B(304)+B(314)+B(317)+0.011*B(376)

  JVS(986) = -B(125)

  JVS(987) = -B(131)

  JVS(988) = 0.035*B(369)

  JVS(989) = 0.07*B(373)

  JVS(990) = 0.347*B(377)

  JVS(991) = 0.009*B(381)

  JVS(992) = 0.143*B(375)

  JVS(993) = 0.011*B(379)

  JVS(994) = 0.016*B(401)+0.051*B(405)

  JVS(995) = 0.18*B(353)+0.25*B(355)

  JVS(996) = 0.09*B(383)+0.001*B(385)+0.176*B(387)

  JVS(997) = 0.25*B(359)+0.18*B(361)+0.25*B(363)

  JVS(998) = 0.093*B(329)+0.008*B(347)+0.064*B(349)+0.01*B(351)

  JVS(999) = 0.041*B(296)+0.051*B(300)

  JVS(1000) = 0.082*B(391)+0.002*B(393)+0.136*B(395)+0.001*B(397)

  JVS(1001) = 0.025*B(289)

  JVS(1002) = 0.173*B(306)+0.095*B(308)

  JVS(1003) = 0.01*B(352)+0.001*B(398)

  JVS(1004) = 0.001*B(232)

  JVS(1005) = 0.042*B(240)

  JVS(1006) = 0.07*B(303)+0.04*B(305)

  JVS(1007) = 0.008*B(348)+0.18*B(354)+0.18*B(362)+0.001*B(386)+0.002*B(394)

  JVS(1008) = -B(106)

  JVS(1009) = 0

  JVS(1010) = -B(112)+0.051*B(301)+0.064*B(350)+0.25*B(356)+0.25*B(364)+0.176*B(388)+0.136*B(396)+0.051*B(406)

  JVS(1011) = -B(172)

  JVS(1012) = -B(152)

  JVS(1013) = -B(194)

  JVS(1014) = -B(110)

  JVS(1015) = -B(134)

  JVS(1016) = -B(107)-B(108)-B(111)-B(113)-B(114)-2*B(118)-B(135)-B(153)-B(173)-B(195)

  JVS(1017) = 0.001*B(233)+0.042*B(241)+0.025*B(290)+0.041*B(297)+0.07*B(304)+0.173*B(307)+0.093*B(330)+0.25*B(360)&
                &+0.035*B(370)+0.07*B(374)+0.143*B(376)+0.347*B(378)+0.011*B(380)+0.009*B(382)+0.09*B(384)+0.082*B(392)&
                &+0.016*B(402)

  JVS(1018) = -B(109)

  JVS(1019) = -B(115)

  JVS(1020) = 2*B(32)

  JVS(1021) = -B(319)

  JVS(1022) = -B(367)

  JVS(1023) = 2*B(69)-B(70)

  JVS(1024) = -B(369)

  JVS(1025) = -B(245)

  JVS(1026) = -B(74)

  JVS(1027) = B(38)-B(40)

  JVS(1028) = -B(409)-B(411)

  JVS(1029) = -B(373)

  JVS(1030) = -B(377)

  JVS(1031) = -B(381)

  JVS(1032) = -B(243)

  JVS(1033) = -0.65*B(247)+B(249)

  JVS(1034) = 0.39*B(58)-B(59)

  JVS(1035) = -B(375)

  JVS(1036) = -B(379)

  JVS(1037) = -B(313)

  JVS(1038) = -B(316)

  JVS(1039) = -B(269)

  JVS(1040) = -B(309)+0.5*B(311)

  JVS(1041) = -0.397*B(371)+0.5*B(399)

  JVS(1042) = -B(275)

  JVS(1043) = -0.34*B(250)+B(252)

  JVS(1044) = -B(265)

  JVS(1045) = -B(260)

  JVS(1046) = -B(49)

  JVS(1047) = -B(46)+B(48)

  JVS(1048) = -B(321)+0.12*B(323)

  JVS(1049) = -B(237)

  JVS(1050) = -B(401)+0.32*B(403)

  JVS(1051) = 0.567*B(353)

  JVS(1052) = -B(255)

  JVS(1053) = 0

  JVS(1054) = -B(383)+0.155*B(385)

  JVS(1055) = -B(359)+0.567*B(361)

  JVS(1056) = -B(329)+0.266*B(347)

  JVS(1057) = -B(280)+0.208*B(282)+0.33*B(288)

  JVS(1058) = -B(296)+0.285*B(298)

  JVS(1059) = -B(391)+0.378*B(393)

  JVS(1060) = -B(289)+0.164*B(291)

  JVS(1061) = -B(227)

  JVS(1062) = -B(218)

  JVS(1063) = -B(306)

  JVS(1064) = 0

  JVS(1065) = -B(232)

  JVS(1066) = -B(240)

  JVS(1067) = -B(303)

  JVS(1068) = -B(51)+B(61)+0.208*B(283)+0.164*B(292)+0.285*B(299)+0.5*B(312)+0.12*B(324)+0.266*B(348)+0.567*B(354)+0.567&
                &*B(362)+0.155*B(386)+0.378*B(394)+0.5*B(400)+0.32*B(404)

  JVS(1069) = -B(36)+B(53)

  JVS(1070) = -B(42)

  JVS(1071) = -B(44)+0.8*B(66)

  JVS(1072) = 0

  JVS(1073) = 0

  JVS(1074) = 0

  JVS(1075) = 0

  JVS(1076) = 0

  JVS(1077) = 0

  JVS(1078) = -B(37)-B(41)-B(43)-B(45)-B(47)-B(50)-B(52)-B(60)-B(71)-B(72)-B(75)-B(76)-B(219)-B(228)-B(233)-B(238)&
                &-B(241)-B(244)-B(246)-0.65*B(248)-0.34*B(251)-B(256)-B(261)-B(266)-B(270)-B(276)-B(281)-B(290)-B(297)&
                &-B(304)-B(307)-B(310)-B(314)-B(317)-B(320)-B(322)-B(330)-B(360)-B(368)-B(370)-0.397*B(372)-B(374)-B(376)&
                &-B(378)-B(380)-B(382)-B(384)-B(392)-B(402)-B(410)-B(412)

  JVS(1079) = B(54)+B(62)+0.8*B(67)-B(73)

  JVS(1080) = 0

  JVS(1081) = B(70)

  JVS(1082) = 0.95*B(245)

  JVS(1083) = B(74)

  JVS(1084) = B(39)

  JVS(1085) = 0.5*B(411)

  JVS(1086) = 0.187*B(381)

  JVS(1087) = B(243)

  JVS(1088) = B(249)

  JVS(1089) = B(222)+B(223)

  JVS(1090) = -B(213)

  JVS(1091) = B(57)+0.61*B(58)

  JVS(1092) = 0.224*B(379)

  JVS(1093) = 0.5*B(315)

  JVS(1094) = 0.5*B(318)

  JVS(1095) = 1.5*B(311)

  JVS(1096) = 0.297*B(371)+1.5*B(399)

  JVS(1097) = 0

  JVS(1098) = B(252)

  JVS(1099) = B(259)

  JVS(1100) = B(49)

  JVS(1101) = 0.12*B(323)+0.5*B(327)

  JVS(1102) = 0.06*B(403)

  JVS(1103) = 0.033*B(353)

  JVS(1104) = 2*B(253)+0.63*B(255)+0.63*B(257)

  JVS(1105) = -B(208)

  JVS(1106) = 0.056*B(385)

  JVS(1107) = 0.033*B(361)

  JVS(1108) = 0.907*B(329)

  JVS(1109) = 0.008*B(282)+0.34*B(288)

  JVS(1110) = 0.4*B(298)+1.233*B(302)

  JVS(1111) = 0.003*B(393)+0.013*B(397)

  JVS(1112) = 0.064*B(291)

  JVS(1113) = B(229)

  JVS(1114) = 2*B(216)+B(218)-B(220)+B(225)

  JVS(1115) = 0.113*B(306)+0.341*B(308)

  JVS(1116) = 0.5*B(328)+0.013*B(398)

  JVS(1117) = B(234)

  JVS(1118) = 0

  JVS(1119) = 0.379*B(303)

  JVS(1120) = B(51)-B(61)+0.008*B(283)+0.064*B(292)+0.4*B(299)+1.5*B(312)+0.12*B(324)+0.033*B(354)+0.033*B(362)+0.056&
                &*B(386)+0.003*B(394)+1.5*B(400)+0.06*B(404)

  JVS(1121) = -B(53)+B(78)+B(86)+B(224)

  JVS(1122) = -B(55)

  JVS(1123) = B(44)-B(66)+B(82)+B(90)+B(112)+B(226)+0.63*B(258)

  JVS(1124) = -B(162)

  JVS(1125) = -B(142)

  JVS(1126) = -B(184)

  JVS(1127) = B(79)-B(80)+B(83)+2*B(85)+B(92)+B(110)

  JVS(1128) = -B(124)

  JVS(1129) = -B(108)+B(111)+B(113)+B(114)+B(118)

  JVS(1130) = B(45)+B(50)+B(52)+B(71)-B(72)+B(75)+B(76)+B(219)+B(244)+0.95*B(246)+0.63*B(256)+0.379*B(304)+0.113*B(307)&
                &+0.907*B(330)+0.297*B(372)+0.224*B(380)+0.187*B(382)+0.5*B(412)

  JVS(1131) = -B(54)-B(56)-B(62)-2*B(63)-2*B(64)-B(67)-B(73)-B(81)-B(88)-B(109)-B(125)-B(143)-B(163)-B(185)-B(209)&
                &-B(214)-B(221)

  JVS(1132) = B(87)-B(89)+B(91)+B(93)+B(94)+B(115)

  JVS(1133) = B(367)

  JVS(1134) = 0.965*B(369)

  JVS(1135) = 0.05*B(245)

  JVS(1136) = 0.695*B(373)

  JVS(1137) = 0.653*B(377)

  JVS(1138) = 0.804*B(381)

  JVS(1139) = 0.835*B(375)

  JVS(1140) = 0.765*B(379)

  JVS(1141) = B(315)

  JVS(1142) = B(318)

  JVS(1143) = 0.76*B(269)

  JVS(1144) = B(309)

  JVS(1145) = 0.1*B(371)

  JVS(1146) = 0.34*B(250)

  JVS(1147) = 0.76*B(265)

  JVS(1148) = B(321)+B(325)+0.2*B(327)

  JVS(1149) = 0.984*B(401)+0.949*B(405)

  JVS(1150) = 0.031*B(353)+0.276*B(355)

  JVS(1151) = 0

  JVS(1152) = 0.91*B(383)+0.022*B(385)+0.824*B(387)

  JVS(1153) = 0.75*B(359)+0.031*B(361)+0.276*B(363)

  JVS(1154) = 0.907*B(329)+0.066*B(347)+0.749*B(349)

  JVS(1155) = 0.5*B(280)+0.1*B(282)+0.5*B(284)+0.33*B(288)

  JVS(1156) = 0.67*B(296)+0.048*B(298)+0.799*B(300)

  JVS(1157) = 0.918*B(391)+0.033*B(393)+0.442*B(395)+0.012*B(397)

  JVS(1158) = 0.3*B(289)+0.05*B(291)

  JVS(1159) = 0.376*B(306)+0.564*B(308)

  JVS(1160) = 0.2*B(328)+0.012*B(398)

  JVS(1161) = 0.034*B(232)+B(234)

  JVS(1162) = 0.37*B(240)+B(242)

  JVS(1163) = 0.473*B(303)+0.96*B(305)

  JVS(1164) = 0.1*B(283)+0.05*B(292)+0.048*B(299)+0.066*B(348)+0.031*B(354)+0.031*B(362)+0.022*B(386)+0.033*B(394)

  JVS(1165) = -B(86)+B(140)

  JVS(1166) = 0

  JVS(1167) = -B(90)+B(144)+0.5*B(285)+0.799*B(301)+B(326)+0.749*B(350)+0.276*B(356)+0.276*B(364)+0.824*B(388)+0.442&
                &*B(396)+0.949*B(406)

  JVS(1168) = -B(168)+B(176)

  JVS(1169) = B(141)+B(145)-B(148)+B(154)+2*B(156)+B(177)+B(198)

  JVS(1170) = -B(190)+B(199)

  JVS(1171) = -B(92)

  JVS(1172) = -B(130)+B(155)

  JVS(1173) = -B(114)

  JVS(1174) = 0.034*B(233)+0.37*B(241)+0.05*B(246)+0.34*B(251)+0.76*B(266)+0.76*B(270)+0.5*B(281)+0.3*B(290)+0.67*B(297)&
                &+0.473*B(304)+0.376*B(307)+B(310)+B(322)+0.907*B(330)+0.75*B(360)+B(368)+0.965*B(370)+0.1*B(372)+0.695&
                &*B(374)+0.835*B(376)+0.653*B(378)+0.765*B(380)+0.804*B(382)+0.91*B(384)+0.918*B(392)+0.984*B(402)

  JVS(1175) = -B(88)

  JVS(1176) = -B(87)-B(89)-B(91)-B(93)-2*B(94)-B(115)-B(131)-B(149)-B(169)-B(191)
      
END SUBROUTINE saprc99_simplesom_mosaic_4bin_aq_Jac_SP














SUBROUTINE saprc99_simplesom_mosaic_4bin_aq_KppDecomp( JVS, IER )







      INTEGER  :: IER
      REAL(kind=dp) :: JVS(1176), W(104), a
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
      
END SUBROUTINE saprc99_simplesom_mosaic_4bin_aq_KppDecomp



SUBROUTINE saprc99_simplesom_mosaic_4bin_aq_KppDecompCmplx( JVS, IER )







      INTEGER  :: IER
      DOUBLE COMPLEX :: JVS(1176), W(104), a
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
      
END SUBROUTINE saprc99_simplesom_mosaic_4bin_aq_KppDecompCmplx


SUBROUTINE saprc99_simplesom_mosaic_4bin_aq_KppSolveIndirect( JVS, X )







      INTEGER i, j
      REAL(kind=dp) JVS(1176), X(104), sum

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
      
END SUBROUTINE saprc99_simplesom_mosaic_4bin_aq_KppSolveIndirect


SUBROUTINE saprc99_simplesom_mosaic_4bin_aq_KppSolveCmplx( JVS, X )







      INTEGER i, j
      DOUBLE COMPLEX JVS(1176), X(104), sum

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
      
END SUBROUTINE saprc99_simplesom_mosaic_4bin_aq_KppSolveCmplx













SUBROUTINE saprc99_simplesom_mosaic_4bin_aq_KppSolve ( JVS, X )


  REAL(kind=dp) :: JVS(LU_NONZERO)

  REAL(kind=dp) :: X(NVAR)

  X(21) = X(21)-JVS(114)*X(20)
  X(36) = X(36)-JVS(162)*X(35)
  X(46) = X(46)-JVS(196)*X(38)-JVS(197)*X(45)
  X(47) = X(47)-JVS(200)*X(38)-JVS(201)*X(45)
  X(48) = X(48)-JVS(204)*X(38)-JVS(205)*X(45)
  X(49) = X(49)-JVS(209)*X(38)-JVS(210)*X(45)
  X(51) = X(51)-JVS(217)*X(42)
  X(52) = X(52)-JVS(223)*X(38)-JVS(224)*X(45)
  X(55) = X(55)-JVS(243)*X(45)
  X(56) = X(56)-JVS(249)*X(54)
  X(57) = X(57)-JVS(259)*X(54)-JVS(260)*X(56)
  X(58) = X(58)-JVS(274)*X(54)-JVS(275)*X(56)-JVS(276)*X(57)
  X(59) = X(59)-JVS(289)*X(56)-JVS(290)*X(57)-JVS(291)*X(58)
  X(60) = X(60)-JVS(303)*X(58)-JVS(304)*X(59)
  X(61) = X(61)-JVS(315)*X(59)-JVS(316)*X(60)
  X(62) = X(62)-JVS(326)*X(60)-JVS(327)*X(61)
  X(63) = X(63)-JVS(336)*X(61)-JVS(337)*X(62)
  X(65) = X(65)-JVS(350)*X(62)-JVS(351)*X(63)-JVS(352)*X(64)
  X(66) = X(66)-JVS(358)*X(63)-JVS(359)*X(64)-JVS(360)*X(65)
  X(67) = X(67)-JVS(365)*X(65)-JVS(366)*X(66)
  X(68) = X(68)-JVS(370)*X(38)-JVS(371)*X(45)-JVS(372)*X(46)-JVS(373)*X(47)-JVS(374)*X(48)
  X(69) = X(69)-JVS(382)*X(44)-JVS(383)*X(46)-JVS(384)*X(47)-JVS(385)*X(49)-JVS(386)*X(50)-JVS(387)*X(68)
  X(70) = X(70)-JVS(407)*X(31)-JVS(408)*X(34)-JVS(409)*X(48)-JVS(410)*X(51)-JVS(411)*X(52)-JVS(412)*X(55)-JVS(413)*X(68)
  X(72) = X(72)-JVS(434)*X(29)-JVS(435)*X(35)-JVS(436)*X(36)-JVS(437)*X(37)-JVS(438)*X(44)
  X(75) = X(75)-JVS(458)*X(38)-JVS(459)*X(45)-JVS(460)*X(46)-JVS(461)*X(47)-JVS(462)*X(49)-JVS(463)*X(50)-JVS(464)*X(55)&
            &-JVS(465)*X(71)-JVS(466)*X(74)
  X(76) = X(76)-JVS(476)*X(48)-JVS(477)*X(55)
  X(78) = X(78)-JVS(493)*X(35)-JVS(494)*X(37)-JVS(495)*X(44)-JVS(496)*X(46)-JVS(497)*X(47)-JVS(498)*X(72)-JVS(499)*X(74)&
            &-JVS(500)*X(77)
  X(81) = X(81)-JVS(533)*X(80)
  X(82) = X(82)-JVS(540)*X(80)
  X(84) = X(84)-JVS(552)*X(80)-JVS(553)*X(83)
  X(85) = X(85)-JVS(559)*X(22)-JVS(560)*X(30)-JVS(561)*X(35)-JVS(562)*X(37)-JVS(563)*X(44)-JVS(564)*X(71)-JVS(565)*X(73)&
            &-JVS(566)*X(77)-JVS(567)*X(82)-JVS(568)*X(83)
  X(86) = X(86)-JVS(583)*X(30)-JVS(584)*X(35)-JVS(585)*X(37)-JVS(586)*X(39)-JVS(587)*X(40)-JVS(588)*X(41)-JVS(589)*X(44)&
            &-JVS(590)*X(50)-JVS(591)*X(71)-JVS(592)*X(72)-JVS(593)*X(73)-JVS(594)*X(74)-JVS(595)*X(75)-JVS(596)*X(76)&
            &-JVS(597)*X(77)-JVS(598)*X(79)-JVS(599)*X(80)-JVS(600)*X(81)-JVS(601)*X(82)-JVS(602)*X(83)-JVS(603)*X(84)
  X(87) = X(87)-JVS(622)*X(36)-JVS(623)*X(74)-JVS(624)*X(77)-JVS(625)*X(79)-JVS(626)*X(82)-JVS(627)*X(83)
  X(88) = X(88)-JVS(636)*X(18)-JVS(637)*X(71)-JVS(638)*X(73)-JVS(639)*X(74)-JVS(640)*X(77)-JVS(641)*X(79)-JVS(642)*X(80)&
            &-JVS(643)*X(81)-JVS(644)*X(83)-JVS(645)*X(84)
  X(89) = X(89)-JVS(652)*X(29)-JVS(653)*X(35)-JVS(654)*X(37)-JVS(655)*X(44)-JVS(656)*X(46)-JVS(657)*X(47)-JVS(658)*X(49)&
            &-JVS(659)*X(53)-JVS(660)*X(71)-JVS(661)*X(73)-JVS(662)*X(74)-JVS(663)*X(77)-JVS(664)*X(79)-JVS(665)*X(81)&
            &-JVS(666)*X(82)-JVS(667)*X(83)-JVS(668)*X(84)-JVS(669)*X(87)-JVS(670)*X(88)
  X(90) = X(90)-JVS(682)*X(35)-JVS(683)*X(37)-JVS(684)*X(44)-JVS(685)*X(73)-JVS(686)*X(77)-JVS(687)*X(81)-JVS(688)*X(82)&
            &-JVS(689)*X(83)-JVS(690)*X(84)-JVS(691)*X(87)-JVS(692)*X(88)
  X(91) = X(91)-JVS(703)*X(37)-JVS(704)*X(44)-JVS(705)*X(45)-JVS(706)*X(74)-JVS(707)*X(77)-JVS(708)*X(79)-JVS(709)*X(80)&
            &-JVS(710)*X(82)-JVS(711)*X(83)-JVS(712)*X(84)-JVS(713)*X(87)-JVS(714)*X(88)
  X(92) = X(92)-JVS(727)*X(49)-JVS(728)*X(50)-JVS(729)*X(71)-JVS(730)*X(73)-JVS(731)*X(74)-JVS(732)*X(77)-JVS(733)*X(79)&
            &-JVS(734)*X(80)-JVS(735)*X(81)-JVS(736)*X(82)-JVS(737)*X(83)-JVS(738)*X(84)-JVS(739)*X(88)
  X(93) = X(93)-JVS(750)*X(33)-JVS(751)*X(41)-JVS(752)*X(78)-JVS(753)*X(79)-JVS(754)*X(80)-JVS(755)*X(83)-JVS(756)*X(84)&
            &-JVS(757)*X(86)-JVS(758)*X(87)-JVS(759)*X(88)-JVS(760)*X(90)-JVS(761)*X(91)-JVS(762)*X(92)
  X(94) = X(94)-JVS(775)*X(23)-JVS(776)*X(24)-JVS(777)*X(25)-JVS(778)*X(26)-JVS(779)*X(31)-JVS(780)*X(33)-JVS(781)*X(36)&
            &-JVS(782)*X(41)-JVS(783)*X(42)-JVS(784)*X(43)-JVS(785)*X(51)-JVS(786)*X(70)-JVS(787)*X(74)-JVS(788)*X(75)&
            &-JVS(789)*X(76)-JVS(790)*X(78)-JVS(791)*X(79)-JVS(792)*X(80)-JVS(793)*X(81)-JVS(794)*X(82)-JVS(795)*X(83)&
            &-JVS(796)*X(84)-JVS(797)*X(85)-JVS(798)*X(86)-JVS(799)*X(87)-JVS(800)*X(88)-JVS(801)*X(89)-JVS(802)*X(90)&
            &-JVS(803)*X(91)-JVS(804)*X(92)-JVS(805)*X(93)
  X(95) = X(95)-JVS(817)*X(31)-JVS(818)*X(34)-JVS(819)*X(43)-JVS(820)*X(48)-JVS(821)*X(51)-JVS(822)*X(52)-JVS(823)*X(55)&
            &-JVS(824)*X(68)-JVS(825)*X(70)-JVS(826)*X(71)-JVS(827)*X(73)-JVS(828)*X(74)-JVS(829)*X(75)-JVS(830)*X(76)&
            &-JVS(831)*X(77)-JVS(832)*X(78)-JVS(833)*X(79)-JVS(834)*X(80)-JVS(835)*X(81)-JVS(836)*X(82)-JVS(837)*X(83)&
            &-JVS(838)*X(84)-JVS(839)*X(85)-JVS(840)*X(86)-JVS(841)*X(87)-JVS(842)*X(88)-JVS(843)*X(89)-JVS(844)*X(90)&
            &-JVS(845)*X(91)-JVS(846)*X(92)-JVS(847)*X(93)-JVS(848)*X(94)
  X(96) = X(96)-JVS(859)*X(25)-JVS(860)*X(52)-JVS(861)*X(83)-JVS(862)*X(88)-JVS(863)*X(92)-JVS(864)*X(93)-JVS(865)*X(94)&
            &-JVS(866)*X(95)
  X(97) = X(97)-JVS(876)*X(24)-JVS(877)*X(74)-JVS(878)*X(75)-JVS(879)*X(76)-JVS(880)*X(79)-JVS(881)*X(81)-JVS(882)*X(82)&
            &-JVS(883)*X(83)-JVS(884)*X(84)-JVS(885)*X(88)-JVS(886)*X(89)-JVS(887)*X(90)-JVS(888)*X(91)-JVS(889)*X(92)&
            &-JVS(890)*X(93)-JVS(891)*X(94)-JVS(892)*X(95)-JVS(893)*X(96)
  X(98) = X(98)-JVS(902)*X(26)-JVS(903)*X(80)-JVS(904)*X(81)-JVS(905)*X(82)-JVS(906)*X(83)-JVS(907)*X(84)-JVS(908)*X(88)&
            &-JVS(909)*X(92)-JVS(910)*X(93)-JVS(911)*X(94)-JVS(912)*X(95)-JVS(913)*X(96)-JVS(914)*X(97)
  X(99) = X(99)-JVS(922)*X(19)-JVS(923)*X(36)-JVS(924)*X(40)-JVS(925)*X(44)-JVS(926)*X(71)-JVS(927)*X(72)-JVS(928)*X(73)&
            &-JVS(929)*X(74)-JVS(930)*X(77)-JVS(931)*X(79)-JVS(932)*X(80)-JVS(933)*X(83)-JVS(934)*X(84)-JVS(935)*X(85)&
            &-JVS(936)*X(87)-JVS(937)*X(88)-JVS(938)*X(89)-JVS(939)*X(90)-JVS(940)*X(91)-JVS(941)*X(92)-JVS(942)*X(93)&
            &-JVS(943)*X(94)-JVS(944)*X(95)-JVS(945)*X(96)-JVS(946)*X(97)-JVS(947)*X(98)
  X(100) = X(100)-JVS(954)*X(23)-JVS(955)*X(28)-JVS(956)*X(38)-JVS(957)*X(44)-JVS(958)*X(46)-JVS(959)*X(47)-JVS(960)&
             &*X(68)-JVS(961)*X(72)-JVS(962)*X(74)-JVS(963)*X(77)-JVS(964)*X(79)-JVS(965)*X(81)-JVS(966)*X(82)-JVS(967)&
             &*X(83)-JVS(968)*X(84)-JVS(969)*X(85)-JVS(970)*X(87)-JVS(971)*X(88)-JVS(972)*X(89)-JVS(973)*X(90)-JVS(974)&
             &*X(91)-JVS(975)*X(92)-JVS(976)*X(93)-JVS(977)*X(94)-JVS(978)*X(95)-JVS(979)*X(96)-JVS(980)*X(97)-JVS(981)&
             &*X(98)-JVS(982)*X(99)
  X(101) = X(101)-JVS(988)*X(29)-JVS(989)*X(35)-JVS(990)*X(37)-JVS(991)*X(38)-JVS(992)*X(44)-JVS(993)*X(45)-JVS(994)&
             &*X(73)-JVS(995)*X(74)-JVS(996)*X(77)-JVS(997)*X(79)-JVS(998)*X(80)-JVS(999)*X(82)-JVS(1000)*X(83)-JVS(1001)&
             &*X(84)-JVS(1002)*X(87)-JVS(1003)*X(88)-JVS(1004)*X(89)-JVS(1005)*X(90)-JVS(1006)*X(91)-JVS(1007)*X(92)&
             &-JVS(1008)*X(93)-JVS(1009)*X(94)-JVS(1010)*X(95)-JVS(1011)*X(96)-JVS(1012)*X(97)-JVS(1013)*X(98)-JVS(1014)&
             &*X(99)-JVS(1015)*X(100)
  X(102) = X(102)-JVS(1020)*X(18)-JVS(1021)*X(19)-JVS(1022)*X(22)-JVS(1023)*X(27)-JVS(1024)*X(29)-JVS(1025)*X(30)&
             &-JVS(1026)*X(32)-JVS(1027)*X(33)-JVS(1028)*X(34)-JVS(1029)*X(35)-JVS(1030)*X(37)-JVS(1031)*X(38)-JVS(1032)&
             &*X(39)-JVS(1033)*X(40)-JVS(1034)*X(43)-JVS(1035)*X(44)-JVS(1036)*X(45)-JVS(1037)*X(46)-JVS(1038)*X(47)&
             &-JVS(1039)*X(48)-JVS(1040)*X(49)-JVS(1041)*X(50)-JVS(1042)*X(52)-JVS(1043)*X(53)-JVS(1044)*X(55)-JVS(1045)&
             &*X(68)-JVS(1046)*X(69)-JVS(1047)*X(70)-JVS(1048)*X(71)-JVS(1049)*X(72)-JVS(1050)*X(73)-JVS(1051)*X(74)&
             &-JVS(1052)*X(75)-JVS(1053)*X(76)-JVS(1054)*X(77)-JVS(1055)*X(79)-JVS(1056)*X(80)-JVS(1057)*X(81)-JVS(1058)&
             &*X(82)-JVS(1059)*X(83)-JVS(1060)*X(84)-JVS(1061)*X(85)-JVS(1062)*X(86)-JVS(1063)*X(87)-JVS(1064)*X(88)&
             &-JVS(1065)*X(89)-JVS(1066)*X(90)-JVS(1067)*X(91)-JVS(1068)*X(92)-JVS(1069)*X(93)-JVS(1070)*X(94)-JVS(1071)&
             &*X(95)-JVS(1072)*X(96)-JVS(1073)*X(97)-JVS(1074)*X(98)-JVS(1075)*X(99)-JVS(1076)*X(100)-JVS(1077)*X(101)
  X(103) = X(103)-JVS(1081)*X(27)-JVS(1082)*X(30)-JVS(1083)*X(32)-JVS(1084)*X(33)-JVS(1085)*X(34)-JVS(1086)*X(38)&
             &-JVS(1087)*X(39)-JVS(1088)*X(40)-JVS(1089)*X(41)-JVS(1090)*X(42)-JVS(1091)*X(43)-JVS(1092)*X(45)-JVS(1093)&
             &*X(46)-JVS(1094)*X(47)-JVS(1095)*X(49)-JVS(1096)*X(50)-JVS(1097)*X(51)-JVS(1098)*X(53)-JVS(1099)*X(68)&
             &-JVS(1100)*X(69)-JVS(1101)*X(71)-JVS(1102)*X(73)-JVS(1103)*X(74)-JVS(1104)*X(75)-JVS(1105)*X(76)-JVS(1106)&
             &*X(77)-JVS(1107)*X(79)-JVS(1108)*X(80)-JVS(1109)*X(81)-JVS(1110)*X(82)-JVS(1111)*X(83)-JVS(1112)*X(84)&
             &-JVS(1113)*X(85)-JVS(1114)*X(86)-JVS(1115)*X(87)-JVS(1116)*X(88)-JVS(1117)*X(89)-JVS(1118)*X(90)-JVS(1119)&
             &*X(91)-JVS(1120)*X(92)-JVS(1121)*X(93)-JVS(1122)*X(94)-JVS(1123)*X(95)-JVS(1124)*X(96)-JVS(1125)*X(97)&
             &-JVS(1126)*X(98)-JVS(1127)*X(99)-JVS(1128)*X(100)-JVS(1129)*X(101)-JVS(1130)*X(102)
  X(104) = X(104)-JVS(1133)*X(22)-JVS(1134)*X(29)-JVS(1135)*X(30)-JVS(1136)*X(35)-JVS(1137)*X(37)-JVS(1138)*X(38)&
             &-JVS(1139)*X(44)-JVS(1140)*X(45)-JVS(1141)*X(46)-JVS(1142)*X(47)-JVS(1143)*X(48)-JVS(1144)*X(49)-JVS(1145)&
             &*X(50)-JVS(1146)*X(53)-JVS(1147)*X(55)-JVS(1148)*X(71)-JVS(1149)*X(73)-JVS(1150)*X(74)-JVS(1151)*X(76)&
             &-JVS(1152)*X(77)-JVS(1153)*X(79)-JVS(1154)*X(80)-JVS(1155)*X(81)-JVS(1156)*X(82)-JVS(1157)*X(83)-JVS(1158)&
             &*X(84)-JVS(1159)*X(87)-JVS(1160)*X(88)-JVS(1161)*X(89)-JVS(1162)*X(90)-JVS(1163)*X(91)-JVS(1164)*X(92)&
             &-JVS(1165)*X(93)-JVS(1166)*X(94)-JVS(1167)*X(95)-JVS(1168)*X(96)-JVS(1169)*X(97)-JVS(1170)*X(98)-JVS(1171)&
             &*X(99)-JVS(1172)*X(100)-JVS(1173)*X(101)-JVS(1174)*X(102)-JVS(1175)*X(103)
  X(104) = X(104)/JVS(1176)
  X(103) = (X(103)-JVS(1132)*X(104))/(JVS(1131))
  X(102) = (X(102)-JVS(1079)*X(103)-JVS(1080)*X(104))/(JVS(1078))
  X(101) = (X(101)-JVS(1017)*X(102)-JVS(1018)*X(103)-JVS(1019)*X(104))/(JVS(1016))
  X(100) = (X(100)-JVS(984)*X(101)-JVS(985)*X(102)-JVS(986)*X(103)-JVS(987)*X(104))/(JVS(983))
  X(99) = (X(99)-JVS(949)*X(100)-JVS(950)*X(101)-JVS(951)*X(102)-JVS(952)*X(103)-JVS(953)*X(104))/(JVS(948))
  X(98) = (X(98)-JVS(916)*X(99)-JVS(917)*X(100)-JVS(918)*X(101)-JVS(919)*X(102)-JVS(920)*X(103)-JVS(921)*X(104))&
            &/(JVS(915))
  X(97) = (X(97)-JVS(895)*X(98)-JVS(896)*X(99)-JVS(897)*X(100)-JVS(898)*X(101)-JVS(899)*X(102)-JVS(900)*X(103)-JVS(901)&
            &*X(104))/(JVS(894))
  X(96) = (X(96)-JVS(868)*X(97)-JVS(869)*X(98)-JVS(870)*X(99)-JVS(871)*X(100)-JVS(872)*X(101)-JVS(873)*X(102)-JVS(874)&
            &*X(103)-JVS(875)*X(104))/(JVS(867))
  X(95) = (X(95)-JVS(850)*X(96)-JVS(851)*X(97)-JVS(852)*X(98)-JVS(853)*X(99)-JVS(854)*X(100)-JVS(855)*X(101)-JVS(856)&
            &*X(102)-JVS(857)*X(103)-JVS(858)*X(104))/(JVS(849))
  X(94) = (X(94)-JVS(807)*X(95)-JVS(808)*X(96)-JVS(809)*X(97)-JVS(810)*X(98)-JVS(811)*X(99)-JVS(812)*X(100)-JVS(813)&
            &*X(101)-JVS(814)*X(102)-JVS(815)*X(103)-JVS(816)*X(104))/(JVS(806))
  X(93) = (X(93)-JVS(764)*X(94)-JVS(765)*X(95)-JVS(766)*X(96)-JVS(767)*X(97)-JVS(768)*X(98)-JVS(769)*X(99)-JVS(770)&
            &*X(100)-JVS(771)*X(101)-JVS(772)*X(102)-JVS(773)*X(103)-JVS(774)*X(104))/(JVS(763))
  X(92) = (X(92)-JVS(741)*X(93)-JVS(742)*X(94)-JVS(743)*X(95)-JVS(744)*X(96)-JVS(745)*X(97)-JVS(746)*X(98)-JVS(747)&
            &*X(100)-JVS(748)*X(102)-JVS(749)*X(103))/(JVS(740))
  X(91) = (X(91)-JVS(716)*X(92)-JVS(717)*X(93)-JVS(718)*X(94)-JVS(719)*X(95)-JVS(720)*X(96)-JVS(721)*X(97)-JVS(722)&
            &*X(99)-JVS(723)*X(100)-JVS(724)*X(101)-JVS(725)*X(102)-JVS(726)*X(104))/(JVS(715))
  X(90) = (X(90)-JVS(694)*X(91)-JVS(695)*X(92)-JVS(696)*X(93)-JVS(697)*X(94)-JVS(698)*X(95)-JVS(699)*X(99)-JVS(700)&
            &*X(101)-JVS(701)*X(102)-JVS(702)*X(104))/(JVS(693))
  X(89) = (X(89)-JVS(672)*X(90)-JVS(673)*X(91)-JVS(674)*X(92)-JVS(675)*X(93)-JVS(676)*X(94)-JVS(677)*X(95)-JVS(678)&
            &*X(101)-JVS(679)*X(102)-JVS(680)*X(103)-JVS(681)*X(104))/(JVS(671))
  X(88) = (X(88)-JVS(647)*X(92)-JVS(648)*X(93)-JVS(649)*X(94)-JVS(650)*X(95)-JVS(651)*X(102))/(JVS(646))
  X(87) = (X(87)-JVS(629)*X(88)-JVS(630)*X(92)-JVS(631)*X(93)-JVS(632)*X(94)-JVS(633)*X(95)-JVS(634)*X(101)-JVS(635)&
            &*X(102))/(JVS(628))
  X(86) = (X(86)-JVS(605)*X(87)-JVS(606)*X(88)-JVS(607)*X(90)-JVS(608)*X(91)-JVS(609)*X(92)-JVS(610)*X(93)-JVS(611)&
            &*X(94)-JVS(612)*X(95)-JVS(613)*X(96)-JVS(614)*X(97)-JVS(615)*X(98)-JVS(616)*X(99)-JVS(617)*X(100)-JVS(618)&
            &*X(101)-JVS(619)*X(102)-JVS(620)*X(103)-JVS(621)*X(104))/(JVS(604))
  X(85) = (X(85)-JVS(570)*X(87)-JVS(571)*X(88)-JVS(572)*X(89)-JVS(573)*X(90)-JVS(574)*X(91)-JVS(575)*X(92)-JVS(576)&
            &*X(93)-JVS(577)*X(95)-JVS(578)*X(96)-JVS(579)*X(97)-JVS(580)*X(98)-JVS(581)*X(100)-JVS(582)*X(102))/(JVS(569))
  X(84) = (X(84)-JVS(555)*X(88)-JVS(556)*X(92)-JVS(557)*X(95)-JVS(558)*X(102))/(JVS(554))
  X(83) = (X(83)-JVS(548)*X(88)-JVS(549)*X(92)-JVS(550)*X(95)-JVS(551)*X(102))/(JVS(547))
  X(82) = (X(82)-JVS(542)*X(83)-JVS(543)*X(88)-JVS(544)*X(92)-JVS(545)*X(95)-JVS(546)*X(102))/(JVS(541))
  X(81) = (X(81)-JVS(535)*X(83)-JVS(536)*X(88)-JVS(537)*X(92)-JVS(538)*X(95)-JVS(539)*X(102))/(JVS(534))
  X(80) = (X(80)-JVS(529)*X(88)-JVS(530)*X(92)-JVS(531)*X(95)-JVS(532)*X(102))/(JVS(528))
  X(79) = (X(79)-JVS(524)*X(88)-JVS(525)*X(92)-JVS(526)*X(95)-JVS(527)*X(102))/(JVS(523))
  X(78) = (X(78)-JVS(502)*X(79)-JVS(503)*X(80)-JVS(504)*X(83)-JVS(505)*X(84)-JVS(506)*X(87)-JVS(507)*X(88)-JVS(508)&
            &*X(90)-JVS(509)*X(91)-JVS(510)*X(92)-JVS(511)*X(93)-JVS(512)*X(94)-JVS(513)*X(95)-JVS(514)*X(96)-JVS(515)*X(97)&
            &-JVS(516)*X(98)-JVS(517)*X(99)-JVS(518)*X(100)-JVS(519)*X(101)-JVS(520)*X(102)-JVS(521)*X(103)-JVS(522)*X(104))&
            &/(JVS(501))
  X(77) = (X(77)-JVS(489)*X(88)-JVS(490)*X(92)-JVS(491)*X(95)-JVS(492)*X(102))/(JVS(488))
  X(76) = (X(76)-JVS(479)*X(93)-JVS(480)*X(94)-JVS(481)*X(95)-JVS(482)*X(96)-JVS(483)*X(97)-JVS(484)*X(98)-JVS(485)&
            &*X(100)-JVS(486)*X(102)-JVS(487)*X(103))/(JVS(478))
  X(75) = (X(75)-JVS(468)*X(76)-JVS(469)*X(79)-JVS(470)*X(82)-JVS(471)*X(88)-JVS(472)*X(92)-JVS(473)*X(95)-JVS(474)&
            &*X(102)-JVS(475)*X(103))/(JVS(467))
  X(74) = (X(74)-JVS(455)*X(88)-JVS(456)*X(92)-JVS(457)*X(95))/(JVS(454))
  X(73) = (X(73)-JVS(450)*X(88)-JVS(451)*X(92)-JVS(452)*X(95)-JVS(453)*X(102))/(JVS(449))
  X(72) = (X(72)-JVS(440)*X(74)-JVS(441)*X(77)-JVS(442)*X(79)-JVS(443)*X(83)-JVS(444)*X(87)-JVS(445)*X(92)-JVS(446)&
            &*X(94)-JVS(447)*X(95)-JVS(448)*X(102))/(JVS(439))
  X(71) = (X(71)-JVS(430)*X(88)-JVS(431)*X(92)-JVS(432)*X(95)-JVS(433)*X(102))/(JVS(429))
  X(70) = (X(70)-JVS(415)*X(75)-JVS(416)*X(76)-JVS(417)*X(81)-JVS(418)*X(82)-JVS(419)*X(83)-JVS(420)*X(84)-JVS(421)&
            &*X(85)-JVS(422)*X(86)-JVS(423)*X(89)-JVS(424)*X(92)-JVS(425)*X(94)-JVS(426)*X(95)-JVS(427)*X(102)-JVS(428)&
            &*X(103))/(JVS(414))
  X(69) = (X(69)-JVS(389)*X(71)-JVS(390)*X(73)-JVS(391)*X(74)-JVS(392)*X(75)-JVS(393)*X(77)-JVS(394)*X(79)-JVS(395)&
            &*X(80)-JVS(396)*X(81)-JVS(397)*X(82)-JVS(398)*X(83)-JVS(399)*X(84)-JVS(400)*X(85)-JVS(401)*X(86)-JVS(402)*X(88)&
            &-JVS(403)*X(89)-JVS(404)*X(92)-JVS(405)*X(95)-JVS(406)*X(102))/(JVS(388))
  X(68) = (X(68)-JVS(376)*X(81)-JVS(377)*X(82)-JVS(378)*X(84)-JVS(379)*X(92)-JVS(380)*X(95)-JVS(381)*X(102))/(JVS(375))
  X(67) = (X(67)-JVS(368)*X(74)-JVS(369)*X(102))/(JVS(367))
  X(66) = (X(66)-JVS(362)*X(67)-JVS(363)*X(74)-JVS(364)*X(102))/(JVS(361))
  X(65) = (X(65)-JVS(354)*X(66)-JVS(355)*X(67)-JVS(356)*X(74)-JVS(357)*X(102))/(JVS(353))
  X(64) = (X(64)-JVS(346)*X(66)-JVS(347)*X(67)-JVS(348)*X(74)-JVS(349)*X(102))/(JVS(345))
  X(63) = (X(63)-JVS(339)*X(64)-JVS(340)*X(65)-JVS(341)*X(66)-JVS(342)*X(67)-JVS(343)*X(74)-JVS(344)*X(102))/(JVS(338))
  X(62) = (X(62)-JVS(329)*X(63)-JVS(330)*X(64)-JVS(331)*X(65)-JVS(332)*X(66)-JVS(333)*X(67)-JVS(334)*X(74)-JVS(335)&
            &*X(102))/(JVS(328))
  X(61) = (X(61)-JVS(318)*X(62)-JVS(319)*X(63)-JVS(320)*X(64)-JVS(321)*X(65)-JVS(322)*X(66)-JVS(323)*X(67)-JVS(324)&
            &*X(74)-JVS(325)*X(102))/(JVS(317))
  X(60) = (X(60)-JVS(306)*X(61)-JVS(307)*X(62)-JVS(308)*X(63)-JVS(309)*X(64)-JVS(310)*X(65)-JVS(311)*X(66)-JVS(312)&
            &*X(67)-JVS(313)*X(74)-JVS(314)*X(102))/(JVS(305))
  X(59) = (X(59)-JVS(293)*X(60)-JVS(294)*X(61)-JVS(295)*X(62)-JVS(296)*X(63)-JVS(297)*X(64)-JVS(298)*X(65)-JVS(299)&
            &*X(66)-JVS(300)*X(67)-JVS(301)*X(74)-JVS(302)*X(102))/(JVS(292))
  X(58) = (X(58)-JVS(278)*X(59)-JVS(279)*X(60)-JVS(280)*X(61)-JVS(281)*X(62)-JVS(282)*X(63)-JVS(283)*X(64)-JVS(284)&
            &*X(65)-JVS(285)*X(66)-JVS(286)*X(67)-JVS(287)*X(74)-JVS(288)*X(102))/(JVS(277))
  X(57) = (X(57)-JVS(262)*X(58)-JVS(263)*X(59)-JVS(264)*X(60)-JVS(265)*X(61)-JVS(266)*X(62)-JVS(267)*X(63)-JVS(268)&
            &*X(64)-JVS(269)*X(65)-JVS(270)*X(66)-JVS(271)*X(67)-JVS(272)*X(74)-JVS(273)*X(102))/(JVS(261))
  X(56) = (X(56)-JVS(251)*X(57)-JVS(252)*X(58)-JVS(253)*X(59)-JVS(254)*X(60)-JVS(255)*X(61)-JVS(256)*X(62)-JVS(257)&
            &*X(63)-JVS(258)*X(102))/(JVS(250))
  X(55) = (X(55)-JVS(245)*X(76)-JVS(246)*X(95)-JVS(247)*X(102)-JVS(248)*X(103))/(JVS(244))
  X(54) = (X(54)-JVS(236)*X(57)-JVS(237)*X(58)-JVS(238)*X(59)-JVS(239)*X(60)-JVS(240)*X(61)-JVS(241)*X(62)-JVS(242)&
            &*X(102))/(JVS(235))
  X(53) = (X(53)-JVS(231)*X(101)-JVS(232)*X(102)-JVS(233)*X(103)-JVS(234)*X(104))/(JVS(230))
  X(52) = (X(52)-JVS(226)*X(83)-JVS(227)*X(92)-JVS(228)*X(95)-JVS(229)*X(102))/(JVS(225))
  X(51) = (X(51)-JVS(219)*X(76)-JVS(220)*X(94)-JVS(221)*X(95)-JVS(222)*X(103))/(JVS(218))
  X(50) = (X(50)-JVS(215)*X(92)-JVS(216)*X(102))/(JVS(214))
  X(49) = (X(49)-JVS(212)*X(92)-JVS(213)*X(102))/(JVS(211))
  X(48) = (X(48)-JVS(207)*X(95)-JVS(208)*X(102))/(JVS(206))
  X(47) = (X(47)-JVS(203)*X(102))/(JVS(202))
  X(46) = (X(46)-JVS(199)*X(102))/(JVS(198))
  X(45) = (X(45)-JVS(195)*X(102))/(JVS(194))
  X(44) = (X(44)-JVS(193)*X(102))/(JVS(192))
  X(43) = (X(43)-JVS(189)*X(94)-JVS(190)*X(102)-JVS(191)*X(103))/(JVS(188))
  X(42) = (X(42)-JVS(184)*X(51)-JVS(185)*X(94)-JVS(186)*X(95)-JVS(187)*X(103))/(JVS(183))
  X(41) = (X(41)-JVS(180)*X(86)-JVS(181)*X(93)-JVS(182)*X(103))/(JVS(179))
  X(40) = (X(40)-JVS(176)*X(99)-JVS(177)*X(102)-JVS(178)*X(103))/(JVS(175))
  X(39) = (X(39)-JVS(171)*X(99)-JVS(172)*X(101)-JVS(173)*X(102)-JVS(174)*X(104))/(JVS(170))
  X(38) = (X(38)-JVS(169)*X(102))/(JVS(168))
  X(37) = (X(37)-JVS(167)*X(102))/(JVS(166))
  X(36) = (X(36)-JVS(164)*X(94)-JVS(165)*X(102))/(JVS(163))
  X(35) = (X(35)-JVS(161)*X(102))/(JVS(160))
  X(34) = (X(34)-JVS(158)*X(95)-JVS(159)*X(102))/(JVS(157))
  X(33) = (X(33)-JVS(155)*X(93)-JVS(156)*X(102))/(JVS(154))
  X(32) = (X(32)-JVS(151)*X(34)-JVS(152)*X(95)-JVS(153)*X(102))/(JVS(150))
  X(31) = (X(31)-JVS(148)*X(94)-JVS(149)*X(95))/(JVS(147))
  X(30) = (X(30)-JVS(146)*X(102))/(JVS(145))
  X(29) = (X(29)-JVS(144)*X(102))/(JVS(143))
  X(28) = (X(28)-JVS(138)*X(38)-JVS(139)*X(74)-JVS(140)*X(79)-JVS(141)*X(92)-JVS(142)*X(102))/(JVS(137))
  X(27) = (X(27)-JVS(135)*X(102)-JVS(136)*X(103))/(JVS(134))
  X(26) = (X(26)-JVS(132)*X(94)-JVS(133)*X(98))/(JVS(131))
  X(25) = (X(25)-JVS(129)*X(94)-JVS(130)*X(96))/(JVS(128))
  X(24) = (X(24)-JVS(126)*X(94)-JVS(127)*X(97))/(JVS(125))
  X(23) = (X(23)-JVS(123)*X(94)-JVS(124)*X(100))/(JVS(122))
  X(22) = (X(22)-JVS(121)*X(102))/(JVS(120))
  X(21) = (X(21)-JVS(116)*X(80)-JVS(117)*X(93)-JVS(118)*X(102)-JVS(119)*X(103))/(JVS(115))
  X(20) = (X(20)-JVS(111)*X(21)-JVS(112)*X(102)-JVS(113)*X(103))/(JVS(110))
  X(19) = (X(19)-JVS(109)*X(102))/(JVS(108))
  X(18) = (X(18)-JVS(107)*X(92))/(JVS(106))
  X(17) = (X(17)-JVS(105)*X(102))/(JVS(104))
  X(16) = (X(16)-JVS(102)*X(20)-JVS(103)*X(102))/(JVS(101))
  X(15) = (X(15)-JVS(97)*X(20)-JVS(98)*X(93)-JVS(99)*X(102)-JVS(100)*X(103))/(JVS(96))
  X(14) = (X(14)-JVS(90)*X(15)-JVS(91)*X(16)-JVS(92)*X(17)-JVS(93)*X(21)-JVS(94)*X(93)-JVS(95)*X(102))/(JVS(89))
  X(13) = (X(13)-JVS(88)*X(15))/(JVS(87))
  X(12) = (X(12)-JVS(85)*X(15)-JVS(86)*X(103))/(JVS(84))
  X(11) = (X(11)-JVS(79)*X(78)-JVS(80)*X(99)-JVS(81)*X(101)-JVS(82)*X(103)-JVS(83)*X(104))/(JVS(78))
  X(10) = (X(10)-JVS(72)*X(78)-JVS(73)*X(93)-JVS(74)*X(95)-JVS(75)*X(99)-JVS(76)*X(101)-JVS(77)*X(104))/(JVS(71))
  X(9) = (X(9)-JVS(63)*X(42)-JVS(64)*X(52)-JVS(65)*X(73)-JVS(66)*X(88)-JVS(67)*X(92)-JVS(68)*X(94)-JVS(69)*X(95)-JVS(70)&
           &*X(102))/(JVS(62))
  X(8) = (X(8)-JVS(58)*X(42)-JVS(59)*X(73)-JVS(60)*X(94)-JVS(61)*X(95))/(JVS(57))
  X(7) = (X(7)-JVS(53)*X(96)-JVS(54)*X(97)-JVS(55)*X(98)-JVS(56)*X(103))/(JVS(52))
  X(6) = (X(6)-JVS(50)*X(100)-JVS(51)*X(103))/(JVS(49))
  X(5) = (X(5)-JVS(46)*X(50)-JVS(47)*X(73)-JVS(48)*X(92))/(JVS(45))
  X(4) = (X(4)-JVS(31)*X(74)-JVS(32)*X(77)-JVS(33)*X(79)-JVS(34)*X(80)-JVS(35)*X(82)-JVS(36)*X(83)-JVS(37)*X(92)-JVS(38)&
           &*X(96)-JVS(39)*X(97)-JVS(40)*X(98)-JVS(41)*X(99)-JVS(42)*X(101)-JVS(43)*X(103)-JVS(44)*X(104))/(JVS(30))
  X(3) = (X(3)-JVS(21)*X(73)-JVS(22)*X(77)-JVS(23)*X(83)-JVS(24)*X(92)-JVS(25)*X(99)-JVS(26)*X(100)-JVS(27)*X(101)&
           &-JVS(28)*X(103)-JVS(29)*X(104))/(JVS(20))
  X(2) = (X(2)-JVS(5)*X(41)-JVS(6)*X(50)-JVS(7)*X(71)-JVS(8)*X(73)-JVS(9)*X(74)-JVS(10)*X(77)-JVS(11)*X(79)-JVS(12)&
           &*X(80)-JVS(13)*X(81)-JVS(14)*X(82)-JVS(15)*X(83)-JVS(16)*X(84)-JVS(17)*X(92)-JVS(18)*X(93)-JVS(19)*X(102))&
           &/(JVS(4))
  X(1) = (X(1)-JVS(2)*X(32)-JVS(3)*X(102))/(JVS(1))
      
END SUBROUTINE saprc99_simplesom_mosaic_4bin_aq_KppSolve
























      SUBROUTINE saprc99_simplesom_mosaic_4bin_aq_WCOPY(N,X,incX,Y,incY)








      
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

      END SUBROUTINE saprc99_simplesom_mosaic_4bin_aq_WCOPY



      SUBROUTINE saprc99_simplesom_mosaic_4bin_aq_WAXPY(N,Alpha,X,incX,Y,incY)









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
      
      END SUBROUTINE saprc99_simplesom_mosaic_4bin_aq_WAXPY




      SUBROUTINE saprc99_simplesom_mosaic_4bin_aq_WSCAL(N,Alpha,X,incX)









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

      END SUBROUTINE saprc99_simplesom_mosaic_4bin_aq_WSCAL


      REAL(kind=dp) FUNCTION saprc99_simplesom_mosaic_4bin_aq_WLAMCH( C )








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
          CALL saprc99_simplesom_mosaic_4bin_aq_WLAMCH_ADD(ONE,Eps,Sum)
          IF (Sum.LE.ONE) GOTO 10
        END DO
        PRINT*,'ERROR IN WLAMCH. EPS < ',Eps
        RETURN
10      Eps = Eps*2
        i = i-1      
      END IF

      saprc99_simplesom_mosaic_4bin_aq_WLAMCH = Eps

      END FUNCTION saprc99_simplesom_mosaic_4bin_aq_WLAMCH
     
      SUBROUTINE saprc99_simplesom_mosaic_4bin_aq_WLAMCH_ADD( A, B, Sum )

      
      REAL(kind=dp) A, B, Sum
      Sum = A + B

      END SUBROUTINE saprc99_simplesom_mosaic_4bin_aq_WLAMCH_ADD




      SUBROUTINE saprc99_simplesom_mosaic_4bin_aq_SET2ZERO(N,Y)




      
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

      END SUBROUTINE saprc99_simplesom_mosaic_4bin_aq_SET2ZERO



      REAL(kind=dp) FUNCTION saprc99_simplesom_mosaic_4bin_aq_WDOT (N, DX, incX, DY, incY) 









      IMPLICIT NONE
      INTEGER :: N, incX, incY
      REAL(kind=dp) :: DX(N), DY(N) 

      INTEGER :: i, IX, IY, M, MP1, NS
                                 
      saprc99_simplesom_mosaic_4bin_aq_WDOT = 0.0D0 
      IF (N .LE. 0) RETURN 
      IF (incX .EQ. incY) IF (incX-1) 5,20,60 



    5 IX = 1 
      IY = 1 
      IF (incX .LT. 0) IX = (-N+1)*incX + 1 
      IF (incY .LT. 0) IY = (-N+1)*incY + 1 
      DO i = 1,N 
        saprc99_simplesom_mosaic_4bin_aq_WDOT = saprc99_simplesom_mosaic_4bin_aq_WDOT + DX(IX)*DY(IY) 
        IX = IX + incX 
        IY = IY + incY 
      END DO 
      RETURN 





   20 M = MOD(N,5) 
      IF (M .EQ. 0) GO TO 40 
      DO i = 1,M 
         saprc99_simplesom_mosaic_4bin_aq_WDOT = saprc99_simplesom_mosaic_4bin_aq_WDOT + DX(i)*DY(i) 
      END DO 
      IF (N .LT. 5) RETURN 
   40 MP1 = M + 1 
      DO i = MP1,N,5 
          saprc99_simplesom_mosaic_4bin_aq_WDOT = saprc99_simplesom_mosaic_4bin_aq_WDOT + DX(i)*DY(i) + DX(i+1)*DY(i+1) +&
                   DX(i+2)*DY(i+2) +  &
                   DX(i+3)*DY(i+3) + DX(i+4)*DY(i+4)                   
      END DO 
      RETURN 



   60 NS = N*incX 
      DO i = 1,NS,incX 
        saprc99_simplesom_mosaic_4bin_aq_WDOT = saprc99_simplesom_mosaic_4bin_aq_WDOT + DX(i)*DY(i) 
      END DO 

      END FUNCTION saprc99_simplesom_mosaic_4bin_aq_WDOT                                          




END MODULE saprc99_simplesom_mosaic_4bin_aq_Integrator
