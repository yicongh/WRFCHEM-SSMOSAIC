
































MODULE saprc99_mosaic_8bin_vbs2_Integrator

 USE saprc99_mosaic_8bin_vbs2_Parameters
 USE saprc99_mosaic_8bin_vbs2_Precision
 USE saprc99_mosaic_8bin_vbs2_JacobianSP

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

SUBROUTINE  saprc99_mosaic_8bin_vbs2_INTEGRATE( TIN, TOUT, &
  FIX, VAR,  RCONST, ATOL, RTOL, IRR_WRK,  &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U  )

   USE saprc99_mosaic_8bin_vbs2_Parameters

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

   CALL saprc99_mosaic_8bin_vbs2_Rosenbrock(VAR, FIX, RCONST, TIN,TOUT,   &
         ATOL,RTOL,               &
         RCNTRL,ICNTRL,RSTATUS,ISTATUS,IRR_WRK,IERR)

   STEPMIN = RCNTRL(ihexit)
   
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U))    IERR_U       = IERR

END SUBROUTINE  saprc99_mosaic_8bin_vbs2_INTEGRATE


SUBROUTINE  saprc99_mosaic_8bin_vbs2_Rosenbrock(Y, FIX, RCONST, Tstart,Tend, &
           AbsTol,RelTol,            &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IRR_WRK,IERR)







































































































  USE saprc99_mosaic_8bin_vbs2_Parameters

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
      CALL saprc99_mosaic_8bin_vbs2_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
      RETURN
   END IF


   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 100000
   ELSEIF (ICNTRL(4) > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL saprc99_mosaic_8bin_vbs2_ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF


   Roundoff = saprc99_mosaic_8bin_vbs2_WLAMCH('E')


   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      Hmin = RCNTRL(1)
   ELSE
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL saprc99_mosaic_8bin_vbs2_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL saprc99_mosaic_8bin_vbs2_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL saprc99_mosaic_8bin_vbs2_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2_dp
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL saprc99_mosaic_8bin_vbs2_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0_dp
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL saprc99_mosaic_8bin_vbs2_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1_dp
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL saprc99_mosaic_8bin_vbs2_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9_dp
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL saprc99_mosaic_8bin_vbs2_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.0_dp*Roundoff) &
         .OR. (RelTol(i) >= 1.0_dp) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL saprc99_mosaic_8bin_vbs2_ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO



   SELECT CASE (Method)
     CASE (1)
       CALL saprc99_mosaic_8bin_vbs2_Ros2(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (2)
       CALL saprc99_mosaic_8bin_vbs2_Ros3(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (3)
       CALL saprc99_mosaic_8bin_vbs2_Ros4(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (4)
       CALL saprc99_mosaic_8bin_vbs2_Rodas3(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (5)
       CALL saprc99_mosaic_8bin_vbs2_Rodas4(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(4)=', Method
       CALL saprc99_mosaic_8bin_vbs2_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       RETURN
   END SELECT


   CALL saprc99_mosaic_8bin_vbs2_ros_Integrator(Y,Tstart,Tend,Texit,      &
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



 SUBROUTINE  saprc99_mosaic_8bin_vbs2_ros_ErrorMsg(Code,T,H,IERR)



   USE saprc99_mosaic_8bin_vbs2_Precision

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

 END SUBROUTINE  saprc99_mosaic_8bin_vbs2_ros_ErrorMsg


 SUBROUTINE  saprc99_mosaic_8bin_vbs2_ros_Integrator (Y, Tstart, Tend, T,     &
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
      CALL saprc99_mosaic_8bin_vbs2_ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1_dp*H) == T).OR.(H <= Roundoff) ) THEN  
      CALL saprc99_mosaic_8bin_vbs2_ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF


   Hexit = H
   H = MIN(H,ABS(Tend-T))


   CALL saprc99_mosaic_8bin_vbs2_FunTemplate(T,Y,Fcn0, RCONST, FIX, Nfun)
   IF( T == Tstart ) THEN
     CALL saprc99_mosaic_8bin_vbs2_IRRFun( Y, FIX, RCONST, IRR_WRK )
   ENDIF


   IF (.NOT.Autonomous) THEN
      CALL saprc99_mosaic_8bin_vbs2_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )
   END IF


   CALL saprc99_mosaic_8bin_vbs2_JacTemplate(T,Y,Jac0, FIX, Njac, RCONST)


UntilAccepted: DO

   CALL saprc99_mosaic_8bin_vbs2_ros_PrepareMatrix(H,Direction,ros_Gamma(1), &
          Jac0,Ghimj,Pivot,Singular, Ndec,  Nsng )
   IF (Singular) THEN 
       CALL saprc99_mosaic_8bin_vbs2_ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF


Stage: DO istage = 1, ros_S

      
       ioffset = NVAR*(istage-1)

      
       IF ( istage == 1 ) THEN
         CALL saprc99_mosaic_8bin_vbs2_WCOPY(NVAR,Fcn0,1,Fcn,1)
      
       ELSEIF ( ros_NewF(istage) ) THEN
         CALL saprc99_mosaic_8bin_vbs2_WCOPY(NVAR,Y,1,Ynew,1)
         DO j = 1, istage-1
           CALL saprc99_mosaic_8bin_vbs2_WAXPY(NVAR,ros_A((istage-1)*(istage-2)/2+j), &
            K(NVAR*(j-1)+1),1,Ynew,1)
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL saprc99_mosaic_8bin_vbs2_FunTemplate(Tau,Ynew,Fcn, RCONST, FIX, Nfun)
       END IF 
       CALL saprc99_mosaic_8bin_vbs2_WCOPY(NVAR,Fcn,1,K(ioffset+1),1)
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         CALL saprc99_mosaic_8bin_vbs2_WAXPY(NVAR,HC,K(NVAR*(j-1)+1),1,K(ioffset+1),1)
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         CALL saprc99_mosaic_8bin_vbs2_WAXPY(NVAR,HG,dFdT,1,K(ioffset+1),1)
       END IF
       CALL saprc99_mosaic_8bin_vbs2_ros_Solve(Ghimj, Pivot, K(ioffset+1), Nsol)

   END DO Stage



   CALL saprc99_mosaic_8bin_vbs2_WCOPY(NVAR,Y,1,Ynew,1)
   DO j=1,ros_S
         CALL saprc99_mosaic_8bin_vbs2_WAXPY(NVAR,ros_M(j),K(NVAR*(j-1)+1),1,Ynew,1)
   END DO


   CALL saprc99_mosaic_8bin_vbs2_WSCAL(NVAR,ZERO,Yerr,1)
   DO j=1,ros_S
        CALL saprc99_mosaic_8bin_vbs2_WAXPY(NVAR,ros_E(j),K(NVAR*(j-1)+1),1,Yerr,1)
   END DO
   Err = saprc99_mosaic_8bin_vbs2_ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )


   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac


   Nstp = Nstp+1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  
      Nacc = Nacc+1
      CALL saprc99_mosaic_8bin_vbs2_WCOPY(NVAR,Ynew,1,Y,1)
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

  END SUBROUTINE  saprc99_mosaic_8bin_vbs2_ros_Integrator



  REAL(kind=dp) FUNCTION  saprc99_mosaic_8bin_vbs2_ros_ErrorNorm ( Y, Ynew, Yerr, &
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

    saprc99_mosaic_8bin_vbs2_ros_ErrorNorm = Err

  END FUNCTION  saprc99_mosaic_8bin_vbs2_ros_ErrorNorm



  SUBROUTINE saprc99_mosaic_8bin_vbs2_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )



   IMPLICIT NONE


   REAL(kind=dp), INTENT(IN) :: T, Roundoff, Y(NVAR), Fcn0(NVAR)
   REAL(kind=dp), INTENT(IN) :: RCONST(NREACT), FIX(NFIX)

   REAL(kind=dp), INTENT(OUT) :: dFdT(NVAR)

   INTEGER, INTENT(INOUT) ::Nfun

   REAL(kind=dp) :: Delta
   REAL(kind=dp), PARAMETER :: ONE = 1.0_dp, DeltaMin = 1.0E-6_dp

   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL saprc99_mosaic_8bin_vbs2_FunTemplate(T+Delta,Y,dFdT, RCONST, FIX, Nfun)
   CALL saprc99_mosaic_8bin_vbs2_WAXPY(NVAR,(-ONE),Fcn0,1,dFdT,1)
   CALL saprc99_mosaic_8bin_vbs2_WSCAL(NVAR,(ONE/Delta),dFdT,1)

  END SUBROUTINE  saprc99_mosaic_8bin_vbs2_ros_FunTimeDeriv



  SUBROUTINE  saprc99_mosaic_8bin_vbs2_ros_PrepareMatrix ( H, Direction, gam, &
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


     CALL saprc99_mosaic_8bin_vbs2_WCOPY(LU_NONZERO,Jac0,1,Ghimj,1)
     CALL saprc99_mosaic_8bin_vbs2_WSCAL(LU_NONZERO,(-ONE),Ghimj,1)
     ghinv = ONE/(Direction*H*gam)
     DO i=1,NVAR
       Ghimj(LU_DIAG(i)) = Ghimj(LU_DIAG(i))+ghinv
     END DO

     CALL saprc99_mosaic_8bin_vbs2_ros_Decomp( Ghimj, Pivot, ising, Ndec )
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

  END SUBROUTINE  saprc99_mosaic_8bin_vbs2_ros_PrepareMatrix



  SUBROUTINE  saprc99_mosaic_8bin_vbs2_ros_Decomp( A, Pivot, ising, Ndec )



   IMPLICIT NONE

   REAL(kind=dp), INTENT(INOUT) :: A(LU_NONZERO)

   INTEGER, INTENT(OUT) :: Pivot(NVAR), ising
   INTEGER, INTENT(INOUT) :: Ndec 

   CALL saprc99_mosaic_8bin_vbs2_KppDecomp ( A, ising )
   Pivot(1) = 1
   Ndec = Ndec + 1

  END SUBROUTINE  saprc99_mosaic_8bin_vbs2_ros_Decomp



  SUBROUTINE  saprc99_mosaic_8bin_vbs2_ros_Solve( A, Pivot, b, Nsol )



   IMPLICIT NONE

   REAL(kind=dp), INTENT(IN) :: A(LU_NONZERO)
   INTEGER, INTENT(IN) :: Pivot(NVAR)

   INTEGER, INTENT(INOUT) :: nsol 

   REAL(kind=dp), INTENT(INOUT) :: b(NVAR)


   CALL saprc99_mosaic_8bin_vbs2_KppSolve( A, b )

   Nsol = Nsol+1

  END SUBROUTINE  saprc99_mosaic_8bin_vbs2_ros_Solve




  SUBROUTINE  saprc99_mosaic_8bin_vbs2_Ros2 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

 END SUBROUTINE  saprc99_mosaic_8bin_vbs2_Ros2



  SUBROUTINE  saprc99_mosaic_8bin_vbs2_Ros3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  saprc99_mosaic_8bin_vbs2_Ros3





  SUBROUTINE  saprc99_mosaic_8bin_vbs2_Ros4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  saprc99_mosaic_8bin_vbs2_Ros4


  SUBROUTINE  saprc99_mosaic_8bin_vbs2_Rodas3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  saprc99_mosaic_8bin_vbs2_Rodas3


  SUBROUTINE  saprc99_mosaic_8bin_vbs2_Rodas4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  saprc99_mosaic_8bin_vbs2_Rodas4




END SUBROUTINE  saprc99_mosaic_8bin_vbs2_Rosenbrock




SUBROUTINE  saprc99_mosaic_8bin_vbs2_FunTemplate( T, Y, Ydot, RCONST, FIX, Nfun )




   USE saprc99_mosaic_8bin_vbs2_Parameters




   REAL(kind=dp) :: T, Y(NVAR)
   REAL(kind=dp) :: RCONST(NREACT)
   REAL(kind=dp) :: FIX(NFIX)

   REAL(kind=dp) :: Ydot(NVAR)
   INTEGER :: Nfun









   CALL saprc99_mosaic_8bin_vbs2_Fun( Y, FIX, RCONST, Ydot )


   Nfun = Nfun+1

END SUBROUTINE  saprc99_mosaic_8bin_vbs2_FunTemplate



SUBROUTINE  saprc99_mosaic_8bin_vbs2_JacTemplate( T, Y, Jcb, FIX, Njac, RCONST )




 USE saprc99_mosaic_8bin_vbs2_Parameters
 
 USE saprc99_mosaic_8bin_vbs2_Jacobian



    REAL(kind=dp) :: T, Y(NVAR)
    REAL(kind=dp) :: FIX(NFIX)
    REAL(kind=dp) :: RCONST(NREACT)

    INTEGER :: Njac


    REAL(kind=dp) :: Jcb(LU_NONZERO)

    REAL(kind=dp) :: Told





    CALL saprc99_mosaic_8bin_vbs2_Jac_SP( Y, FIX, RCONST, Jcb )


    Njac = Njac+1

END SUBROUTINE  saprc99_mosaic_8bin_vbs2_JacTemplate

















SUBROUTINE saprc99_mosaic_8bin_vbs2_Fun ( V, F, RCT, Vdot )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: Vdot(NVAR)




  REAL(kind=dp) :: A(NREACT)


  A(1) = RCT(1)*V(92)
  A(2) = RCT(2)*V(83)*F(2)
  A(3) = RCT(3)*V(83)*V(87)
  A(4) = RCT(4)*V(83)*V(97)*F(2)
  A(5) = RCT(5)*V(83)*V(92)
  A(6) = RCT(6)*V(83)*V(92)
  A(7) = RCT(7)*V(87)*V(97)
  A(8) = RCT(8)*V(87)*V(92)
  A(9) = RCT(9)*V(97)*V(99)
  A(10) = RCT(10)*V(97)*V(97)*F(2)
  A(11) = RCT(11)*V(92)*V(99)
  A(12) = RCT(12)*V(39)
  A(13) = RCT(13)*V(39)*F(1)
  A(14) = RCT(14)*V(92)*V(99)
  A(15) = RCT(15)*V(99)
  A(16) = RCT(16)*V(99)
  A(17) = RCT(17)*V(87)
  A(18) = RCT(18)*V(87)
  A(19) = RCT(19)*V(24)*F(1)
  A(20) = RCT(20)*V(24)*F(2)
  A(21) = RCT(21)*V(97)*V(98)
  A(22) = RCT(22)*V(42)
  A(23) = RCT(23)*V(42)
  A(24) = RCT(24)*V(42)*V(98)
  A(25) = RCT(25)*V(92)*V(98)
  A(26) = RCT(26)*V(98)*V(99)
  A(27) = RCT(27)*V(65)*V(98)
  A(28) = RCT(28)*V(65)
  A(29) = RCT(29)*V(64)*V(98)
  A(30) = RCT(30)*V(87)*V(98)
  A(31) = RCT(31)*V(91)*V(97)
  A(32) = RCT(32)*V(91)*V(92)
  A(33) = RCT(33)*V(50)
  A(34) = RCT(34)*V(50)
  A(35) = RCT(35)*V(50)*V(98)
  A(36) = RCT(36)*V(87)*V(91)
  A(37) = RCT(37)*V(91)*V(91)
  A(38) = RCT(38)*V(91)*V(91)*F(1)
  A(39) = RCT(39)*V(91)*V(99)
  A(40) = RCT(40)*V(99)*V(99)
  A(41) = RCT(41)*V(35)
  A(42) = RCT(42)*V(35)*V(98)
  A(43) = RCT(43)*V(91)*V(98)
  A(44) = RCT(44)*V(40)*V(98)
  A(45) = RCT(45)*V(98)*F(2)
  A(46) = RCT(46)*V(93)*V(97)
  A(47) = RCT(47)*V(91)*V(93)
  A(48) = RCT(48)*V(93)*V(99)
  A(49) = RCT(49)*V(93)*V(93)
  A(50) = RCT(50)*V(93)*V(93)
  A(51) = RCT(51)*V(94)*V(97)
  A(52) = RCT(52)*V(91)*V(94)
  A(53) = RCT(53)*V(94)*V(99)
  A(54) = RCT(54)*V(93)*V(94)
  A(55) = RCT(55)*V(94)*V(94)
  A(56) = RCT(56)*V(72)*V(97)
  A(57) = RCT(57)*V(72)*V(91)
  A(58) = RCT(58)*V(72)*V(99)
  A(59) = RCT(59)*V(72)*V(93)
  A(60) = RCT(60)*V(72)*V(94)
  A(61) = RCT(61)*V(72)*V(72)
  A(62) = RCT(62)*V(88)*V(97)
  A(63) = RCT(63)*V(88)*V(91)
  A(64) = RCT(64)*V(88)*V(93)
  A(65) = RCT(65)*V(88)*V(99)
  A(66) = RCT(66)*V(88)*V(94)
  A(67) = RCT(67)*V(72)*V(88)
  A(68) = RCT(68)*V(88)*V(88)
  A(69) = RCT(69)*V(92)*V(96)
  A(70) = RCT(70)*V(31)
  A(71) = RCT(71)*V(96)*V(97)
  A(72) = RCT(72)*V(91)*V(96)
  A(73) = RCT(73)*V(96)*V(99)
  A(74) = RCT(74)*V(93)*V(96)
  A(75) = RCT(75)*V(94)*V(96)
  A(76) = RCT(76)*V(72)*V(96)
  A(77) = RCT(77)*V(88)*V(96)
  A(78) = RCT(78)*V(96)*V(96)
  A(79) = RCT(79)*V(89)*V(92)
  A(80) = RCT(80)*V(32)
  A(81) = RCT(81)*V(89)*V(97)
  A(82) = RCT(82)*V(89)*V(91)
  A(83) = RCT(83)*V(89)*V(99)
  A(84) = RCT(84)*V(89)*V(93)
  A(85) = RCT(85)*V(89)*V(94)
  A(86) = RCT(86)*V(72)*V(89)
  A(87) = RCT(87)*V(88)*V(89)
  A(88) = RCT(88)*V(89)*V(96)
  A(89) = RCT(89)*V(89)*V(89)
  A(90) = RCT(90)*V(92)*V(95)
  A(91) = RCT(91)*V(33)
  A(92) = RCT(92)*V(95)*V(97)
  A(93) = RCT(93)*V(91)*V(95)
  A(94) = RCT(94)*V(95)*V(99)
  A(95) = RCT(95)*V(93)*V(95)
  A(96) = RCT(96)*V(94)*V(95)
  A(97) = RCT(97)*V(72)*V(95)
  A(98) = RCT(98)*V(88)*V(95)
  A(99) = RCT(99)*V(95)*V(96)
  A(100) = RCT(100)*V(89)*V(95)
  A(101) = RCT(101)*V(95)*V(95)
  A(102) = RCT(102)*V(90)*V(92)
  A(103) = RCT(103)*V(34)
  A(104) = RCT(104)*V(90)*V(97)
  A(105) = RCT(105)*V(90)*V(91)
  A(106) = RCT(106)*V(90)*V(99)
  A(107) = RCT(107)*V(90)*V(93)
  A(108) = RCT(108)*V(90)*V(94)
  A(109) = RCT(109)*V(72)*V(90)
  A(110) = RCT(110)*V(88)*V(90)
  A(111) = RCT(111)*V(90)*V(96)
  A(112) = RCT(112)*V(89)*V(90)
  A(113) = RCT(113)*V(90)*V(95)
  A(114) = RCT(114)*V(90)*V(90)
  A(115) = RCT(115)*V(45)*V(92)
  A(116) = RCT(116)*V(45)
  A(117) = RCT(117)*V(69)*V(92)
  A(118) = RCT(118)*V(69)*V(91)
  A(119) = RCT(119)*V(69)
  A(120) = RCT(120)*V(48)*V(92)
  A(121) = RCT(121)*V(48)*V(91)
  A(122) = RCT(122)*V(48)
  A(123) = RCT(123)*V(81)
  A(124) = RCT(124)*V(81)
  A(125) = RCT(125)*V(81)*V(98)
  A(126) = RCT(126)*V(81)*V(91)
  A(127) = RCT(127)*V(47)
  A(128) = RCT(128)*V(47)*V(97)
  A(129) = RCT(129)*V(81)*V(99)
  A(130) = RCT(130)*V(78)*V(98)
  A(131) = RCT(131)*V(78)
  A(132) = RCT(132)*V(78)*V(99)
  A(133) = RCT(133)*V(84)*V(98)
  A(134) = RCT(134)*V(84)
  A(135) = RCT(135)*V(84)*V(99)
  A(136) = RCT(136)*V(67)*V(98)
  A(137) = RCT(137)*V(67)
  A(138) = RCT(138)*V(85)*V(98)
  A(139) = RCT(139)*V(85)
  A(140) = RCT(140)*V(51)*V(98)
  A(141) = RCT(141)*V(38)*V(98)
  A(142) = RCT(142)*V(46)*V(98)
  A(143) = RCT(143)*V(46)
  A(144) = RCT(144)*V(59)*V(98)
  A(145) = RCT(145)*V(59)
  A(146) = RCT(146)*V(76)
  A(147) = RCT(147)*V(76)
  A(148) = RCT(148)*V(76)*V(98)
  A(149) = RCT(149)*V(76)*V(99)
  A(150) = RCT(150)*V(63)
  A(151) = RCT(151)*V(63)*V(98)
  A(152) = RCT(152)*V(63)*V(99)
  A(153) = RCT(153)*V(37)
  A(154) = RCT(154)*V(62)*V(98)
  A(155) = RCT(155)*V(62)*V(99)
  A(156) = RCT(156)*V(56)*V(98)
  A(157) = RCT(157)*V(56)*V(99)
  A(158) = RCT(158)*V(60)*V(99)
  A(159) = RCT(159)*V(61)*V(98)
  A(160) = RCT(160)*V(61)
  A(161) = RCT(161)*V(61)*V(99)
  A(162) = RCT(162)*V(73)*V(98)
  A(163) = RCT(163)*V(73)*V(87)
  A(164) = RCT(164)*V(73)*V(99)
  A(165) = RCT(165)*V(73)*V(83)
  A(166) = RCT(166)*V(73)
  A(167) = RCT(167)*V(80)*V(98)
  A(168) = RCT(168)*V(80)*V(87)
  A(169) = RCT(169)*V(80)*V(83)
  A(170) = RCT(170)*V(80)
  A(171) = RCT(171)*V(77)*V(98)
  A(172) = RCT(172)*V(77)*V(87)
  A(173) = RCT(173)*V(77)*V(99)
  A(174) = RCT(174)*V(77)
  A(175) = RCT(175)*V(86)*V(98)
  A(176) = RCT(176)*V(86)
  A(177) = RCT(177)*V(82)*V(98)
  A(178) = RCT(178)*V(82)
  A(179) = RCT(179)*V(57)*V(98)
  A(180) = RCT(180)*V(57)*V(87)
  A(181) = RCT(181)*V(54)*V(98)
  A(182) = RCT(182)*V(54)
  A(183) = RCT(183)*V(55)*V(98)
  A(184) = RCT(184)*V(55)
  A(185) = RCT(185)*V(29)*V(98)
  A(186) = RCT(186)*V(66)*V(98)
  A(187) = RCT(187)*V(66)*V(87)
  A(188) = RCT(188)*V(66)*V(99)
  A(189) = RCT(189)*V(66)*V(83)
  A(190) = RCT(190)*V(70)*V(98)
  A(191) = RCT(191)*V(70)*V(87)
  A(192) = RCT(192)*V(70)*V(99)
  A(193) = RCT(193)*V(70)*V(83)
  A(194) = RCT(194)*V(75)*V(98)
  A(195) = RCT(195)*V(75)*V(87)
  A(196) = RCT(196)*V(75)*V(99)
  A(197) = RCT(197)*V(75)*V(83)
  A(198) = RCT(198)*V(74)*V(98)
  A(199) = RCT(199)*V(74)*V(87)
  A(200) = RCT(200)*V(74)*V(99)
  A(201) = RCT(201)*V(74)*V(83)
  A(202) = RCT(202)*V(30)*V(98)
  A(203) = RCT(203)*V(36)*V(98)
  A(204) = RCT(204)*V(58)*V(98)
  A(205) = RCT(205)*V(44)*V(98)
  A(206) = RCT(206)*V(52)*V(98)
  A(207) = RCT(207)*V(43)*V(98)
  A(208) = RCT(208)*V(53)*V(98)
  A(209) = RCT(209)*V(49)*V(98)
  A(210) = RCT(210)*V(71)*V(98)
  A(211) = RCT(211)*V(71)*V(87)
  A(212) = RCT(212)*V(71)*V(99)
  A(213) = RCT(213)*V(71)*V(83)
  A(214) = RCT(214)*V(79)*V(98)
  A(215) = RCT(215)*V(79)*V(87)
  A(216) = RCT(216)*V(79)*V(99)
  A(217) = RCT(217)*V(79)*V(83)
  A(218) = RCT(218)*V(58)*V(87)
  A(219) = RCT(219)*V(68)*V(98)
  A(220) = RCT(220)*V(68)*V(87)
  A(221) = RCT(221)*V(68)*V(99)
  A(222) = RCT(222)*V(68)*V(83)
  A(223) = RCT(223)*V(40)
  A(224) = RCT(224)*V(91)
  A(225) = RCT(225)*V(40)
  A(226) = RCT(226)*V(1)
  A(227) = RCT(227)*V(65)
  A(228) = RCT(228)*V(35)
  A(229) = RCT(229)*V(2)
  A(230) = RCT(230)*V(52)*V(98)
  A(231) = RCT(231)*V(43)*V(98)
  A(232) = RCT(232)*V(71)*V(98)
  A(233) = RCT(233)*V(79)*V(98)
  A(234) = RCT(234)*V(53)*V(98)
  A(235) = RCT(235)*V(49)*V(98)
  A(236) = RCT(236)*V(70)*V(98)
  A(237) = RCT(237)*V(70)*V(87)
  A(238) = RCT(238)*V(70)*V(99)
  A(239) = RCT(239)*V(75)*V(98)
  A(240) = RCT(240)*V(75)*V(87)
  A(241) = RCT(241)*V(75)*V(99)
  A(242) = RCT(242)*V(74)*V(98)
  A(243) = RCT(243)*V(74)*V(87)
  A(244) = RCT(244)*V(74)*V(99)
  A(245) = RCT(245)*V(3)*V(98)
  A(246) = RCT(246)*V(4)*V(98)
  A(247) = RCT(247)*V(26)*V(98)
  A(248) = RCT(248)*V(20)*V(98)
  A(249) = RCT(249)*V(5)*V(98)
  A(250) = RCT(250)*V(6)*V(98)
  A(251) = RCT(251)*V(28)*V(98)
  A(252) = RCT(252)*V(22)*V(98)
  A(253) = RCT(253)*V(25)*V(98)
  A(254) = RCT(254)*V(21)*V(98)
  A(255) = RCT(255)*V(27)*V(98)
  A(256) = RCT(256)*V(23)*V(98)
  A(257) = RCT(257)*V(41)*V(98)
  A(258) = RCT(258)*V(41)*V(98)
  A(259) = RCT(259)*V(41)*V(99)


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
  Vdot(17) = A(236)+A(237)+0.0327*A(238)+A(239)+A(240)+0.0545*A(241)+A(242)+A(243)+0.0816*A(244)
  Vdot(18) = A(21)+A(24)+A(25)+A(26)+A(27)+A(29)+A(30)+A(43)+A(44)+A(45)+A(125)+A(130)+A(133)+A(136)+A(138)+A(140)&
               &+A(141)+A(142)+A(144)+A(148)+A(151)+A(154)+A(156)+A(159)+A(162)+A(167)+A(171)+A(175)+A(177)+A(179)+A(181)&
               &+A(183)+A(185)+A(186)+A(190)+A(194)+A(198)+A(202)+A(203)+A(204)+A(205)+A(206)+A(207)+A(208)+A(209)+A(210)&
               &+A(214)+A(219)
  Vdot(19) = A(245)+A(246)+A(247)+A(248)+A(249)+A(250)+A(251)+A(252)+A(253)+A(254)+A(255)+A(256)
  Vdot(20) = 0.5*A(253)+A(254)
  Vdot(21) = -A(254)
  Vdot(22) = 0.5*A(255)+A(256)
  Vdot(23) = -A(256)
  Vdot(24) = A(18)-A(19)-A(20)
  Vdot(25) = -A(253)
  Vdot(26) = A(253)
  Vdot(27) = -A(255)
  Vdot(28) = A(255)
  Vdot(29) = -A(185)
  Vdot(30) = -A(202)
  Vdot(31) = A(69)-A(70)
  Vdot(32) = A(79)-A(80)
  Vdot(33) = A(90)-A(91)
  Vdot(34) = A(102)-A(103)
  Vdot(35) = A(37)+A(38)-A(41)-A(42)-A(228)
  Vdot(36) = -A(203)
  Vdot(37) = -A(153)+0.031*A(195)+0.031*A(199)+0.087*A(209)
  Vdot(38) = -A(141)
  Vdot(39) = A(11)-A(12)-A(13)
  Vdot(40) = -A(44)-A(223)-A(225)+A(257)+0.5*A(258)+A(259)
  Vdot(41) = -A(257)-A(258)-A(259)
  Vdot(42) = A(21)-A(22)-A(23)-A(24)
  Vdot(43) = -A(207)
  Vdot(44) = -A(205)
  Vdot(45) = -A(115)-A(116)+0.236*A(205)
  Vdot(46) = A(47)-A(142)-A(143)
  Vdot(47) = A(126)-A(127)-A(128)
  Vdot(48) = -A(120)-A(121)-A(122)+A(158)
  Vdot(49) = -A(209)
  Vdot(50) = A(32)-A(33)-A(34)-A(35)
  Vdot(51) = A(49)+0.25*A(54)+0.25*A(64)-A(140)
  Vdot(52) = -A(206)
  Vdot(53) = -A(208)
  Vdot(54) = -A(181)-A(182)+0.108*A(208)+0.099*A(209)
  Vdot(55) = -A(183)-A(184)+0.051*A(208)+0.093*A(209)
  Vdot(56) = -A(156)-A(157)+0.207*A(208)+0.187*A(209)
  Vdot(57) = -A(179)-A(180)+0.491*A(208)+0.561*A(209)
  Vdot(58) = -A(204)-A(218)
  Vdot(59) = A(52)+A(63)-A(144)-A(145)
  Vdot(60) = A(117)+A(121)+A(122)-A(158)
  Vdot(61) = -A(159)-A(160)-A(161)+0.059*A(208)+0.05*A(209)+0.061*A(214)+0.042*A(215)+0.015*A(216)
  Vdot(62) = A(118)+A(119)-A(154)-A(155)+0.017*A(208)
  Vdot(63) = -A(150)-A(151)-A(152)+0.23*A(156)+0.084*A(162)+0.9*A(163)+0.3*A(167)+0.95*A(168)+0.174*A(171)+0.742*A(172)&
               &+0.008*A(173)+0.5*A(182)+0.5*A(184)+0.119*A(208)+0.287*A(209)
  Vdot(64) = -A(29)+A(123)+A(124)+A(125)+A(129)+A(131)+0.034*A(133)+A(134)+2*A(146)+A(147)+1.26*A(148)+1.26*A(149)&
               &+A(150)+A(151)+A(152)+0.416*A(162)+0.45*A(163)+0.5*A(164)+0.67*A(166)+0.475*A(168)+0.7*A(170)+0.336*A(171)&
               &+0.498*A(172)+0.572*A(173)+1.233*A(174)+A(179)+1.5*A(180)+A(182)+A(184)+0.5*A(187)+0.491*A(189)+0.275*A(191)&
               &+0.157*A(195)+0.157*A(199)+0.393*A(204)+0.002*A(206)+0.345*A(211)+0.265*A(215)+0.012*A(217)+1.5*A(218)+0.51&
               &*A(220)
  Vdot(65) = 2*A(13)+A(25)-A(27)-A(28)+0.2*A(39)+A(129)+A(132)+A(135)+A(149)+A(152)+A(155)+A(157)+A(158)+A(161)+0.5&
               &*A(164)+0.15*A(173)-A(227)+A(259)
  Vdot(66) = -A(186)-A(187)-A(188)-A(189)
  Vdot(67) = A(116)-A(136)-A(137)+0.006*A(177)+0.02*A(178)+0.13*A(195)+0.13*A(199)+0.704*A(203)+0.024*A(205)+0.452&
               &*A(206)+0.072*A(207)+0.005*A(210)+0.001*A(211)+0.024*A(212)+0.127*A(214)+0.045*A(215)+0.102*A(216)
  Vdot(68) = -A(219)-A(220)-A(221)-A(222)
  Vdot(69) = A(92)+A(94)+A(99)+A(100)+2*A(101)+A(113)-A(117)-A(118)-A(119)+0.24*A(154)+A(155)+0.24*A(156)+A(157)
  Vdot(70) = -A(190)-A(191)-A(192)-A(193)
  Vdot(71) = -A(210)-A(211)-A(212)-A(213)
  Vdot(72) = -A(56)-A(57)-A(58)-A(59)-A(60)-A(67)-A(76)-A(86)+A(92)+A(94)-A(97)+A(99)+A(100)+2*A(101)-A(109)+A(113)&
               &+A(136)+0.616*A(138)+0.675*A(167)+0.515*A(176)+0.596*A(177)+0.152*A(178)+A(181)+A(182)+A(183)+A(184)+0.079&
               &*A(190)+0.126*A(191)+0.187*A(192)+0.24*A(193)+0.5*A(194)+0.729*A(195)+0.75*A(196)+0.5*A(198)+0.729*A(199)&
               &+0.75*A(200)+0.559*A(205)+0.936*A(206)+0.948*A(207)+0.205*A(210)+0.488*A(212)+0.001*A(214)+0.137*A(215)&
               &+0.711*A(216)
  Vdot(73) = -A(162)-A(163)-A(164)-A(165)-A(166)+0.23*A(190)+0.39*A(191)+0.025*A(214)+0.026*A(215)+0.012*A(217)
  Vdot(74) = -A(198)-A(199)-A(200)-A(201)
  Vdot(75) = -A(194)-A(195)-A(196)-A(197)
  Vdot(76) = -A(146)-A(147)-A(148)-A(149)+0.23*A(154)+0.15*A(171)+0.023*A(172)+A(180)+0.5*A(182)+0.5*A(184)+0.009*A(189)&
               &+0.001*A(195)+0.001*A(199)+0.607*A(204)+0.118*A(208)+0.097*A(209)
  Vdot(77) = -A(171)-A(172)-A(173)-A(174)+0.357*A(190)+0.936*A(192)+0.025*A(214)
  Vdot(78) = A(81)+A(83)+A(88)+2*A(89)+A(100)+A(112)-A(130)-A(131)-A(132)+0.034*A(133)+A(134)+0.482*A(138)+A(139)+0.96&
               &*A(141)+0.129*A(171)+0.047*A(172)+0.467*A(174)+0.084*A(175)+0.246*A(176)+0.439*A(177)+0.431*A(178)+0.195&
               &*A(186)+0.25*A(189)+A(202)+0.445*A(205)+0.455*A(206)+0.099*A(207)+0.294*A(210)+0.154*A(211)+0.009*A(212)&
               &+0.732*A(214)+0.456*A(215)+0.507*A(216)+0.984*A(219)+0.5*A(220)
  Vdot(79) = -A(214)-A(215)-A(216)-A(217)
  Vdot(80) = -A(167)-A(168)-A(169)-A(170)+0.32*A(190)+0.16*A(191)+0.019*A(215)+0.048*A(216)
  Vdot(81) = A(46)+A(48)+A(49)+2*A(50)+0.75*A(54)+0.75*A(64)+A(74)+A(84)+A(95)+A(104)+A(106)+A(107)+A(111)+A(112)+A(113)&
               &+2*A(114)-A(123)-A(124)-A(125)-A(126)+A(127)-A(129)+A(136)+0.115*A(138)+A(140)+0.081*A(141)+0.35*A(142)&
               &+A(143)+A(147)+0.084*A(162)+0.2*A(163)+0.67*A(166)+0.3*A(167)+0.1*A(168)+0.055*A(171)+0.125*A(172)+0.227&
               &*A(173)+0.3*A(174)+0.213*A(175)+0.506*A(176)+0.01*A(177)+0.134*A(178)+1.61*A(186)+A(187)+0.191*A(189)+0.624&
               &*A(190)+0.592*A(191)+0.24*A(193)+0.276*A(194)+0.235*A(195)+0.276*A(198)+0.235*A(199)+0.096*A(204)+0.026&
               &*A(205)+0.024*A(206)+0.026*A(207)+0.732*A(210)+0.5*A(211)+0.244*A(214)+0.269*A(215)+0.079*A(216)+0.984&
               &*A(219)+0.5*A(220)
  Vdot(82) = A(62)+A(115)+0.572*A(173)-0.69*A(177)-A(178)+0.276*A(196)+0.276*A(200)+0.511*A(212)+0.321*A(216)
  Vdot(83) = A(1)-A(2)-A(3)-A(4)-A(5)-A(6)+A(16)+A(17)+A(20)-A(165)-A(169)-A(189)-A(193)-A(197)-A(201)-A(213)-A(217)&
               &-A(222)
  Vdot(84) = -A(133)-A(134)-A(135)+0.37*A(138)+A(144)+A(145)+A(165)+0.675*A(167)+0.45*A(169)+0.013*A(171)+0.218*A(173)&
               &+0.558*A(175)+0.71*A(176)+0.213*A(177)+0.147*A(178)+A(179)+A(181)+A(183)+A(188)+0.474*A(194)+0.205*A(195)&
               &+0.474*A(196)+0.147*A(197)+0.474*A(198)+0.205*A(199)+0.474*A(200)+0.147*A(201)+0.261*A(203)+0.122*A(205)&
               &+0.244*A(206)+0.204*A(207)+0.497*A(210)+0.363*A(211)+0.037*A(212)+0.45*A(213)+0.511*A(214)+0.305*A(215)&
               &+0.151*A(216)+0.069*A(217)+0.45*A(222)
  Vdot(85) = 0.5*A(64)+A(65)+0.5*A(66)+A(68)-A(138)-A(139)+0.416*A(162)+0.55*A(169)+0.15*A(171)+0.21*A(172)+0.233*A(174)&
               &+0.115*A(175)+0.177*A(177)+0.243*A(178)+0.332*A(205)+0.11*A(206)+0.089*A(207)+0.437*A(213)+0.072*A(214)&
               &+0.026*A(215)+0.001*A(216)+0.659*A(217)+0.55*A(222)
  Vdot(86) = 0.5*A(64)+0.5*A(66)+A(68)+A(77)+A(87)+A(98)+0.7*A(170)+0.332*A(171)-0.671*A(175)-A(176)+0.048*A(177)+0.435&
               &*A(178)+0.1*A(191)+0.75*A(193)+0.276*A(194)+0.276*A(195)+0.853*A(197)+0.276*A(198)+0.276*A(199)+0.853*A(201)&
               &+0.125*A(206)+0.417*A(207)+0.055*A(208)+0.119*A(210)+0.215*A(211)+0.113*A(213)+0.043*A(215)+0.259*A(217)
  Vdot(87) = A(2)-A(3)-A(7)-A(8)-A(17)-A(18)-A(30)-A(36)+0.25*A(72)+0.25*A(82)+0.25*A(93)+0.25*A(105)-A(163)-A(168)&
               &-A(172)-A(180)-A(187)-A(191)-A(195)-A(199)-A(211)-A(215)-A(218)-A(220)
  Vdot(88) = -A(62)-A(63)-A(64)-A(65)-A(66)-2*A(68)-A(77)-A(87)-A(98)-A(110)+0.001*A(133)+0.042*A(138)+0.025*A(167)&
               &+0.041*A(171)+0.051*A(173)+0.07*A(175)+0.04*A(176)+0.173*A(177)+0.095*A(178)+0.093*A(190)+0.008*A(191)+0.064&
               &*A(192)+0.01*A(193)+0.25*A(194)+0.18*A(195)+0.25*A(196)+0.25*A(198)+0.18*A(199)+0.25*A(200)+0.035*A(203)&
               &+0.07*A(205)+0.143*A(206)+0.347*A(207)+0.011*A(208)+0.009*A(209)+0.09*A(210)+0.001*A(211)+0.176*A(212)+0.082&
               &*A(214)+0.002*A(215)+0.136*A(216)+0.001*A(217)+0.016*A(219)+0.051*A(221)
  Vdot(89) = -A(79)+A(80)-A(81)-A(82)-A(83)-A(84)-A(85)-A(87)-A(88)-2*A(89)-A(100)-A(112)+0.965*A(133)+A(135)+0.096&
               &*A(138)+0.37*A(148)+0.37*A(149)+0.1*A(163)+0.05*A(168)+0.048*A(172)+0.3*A(174)+0.049*A(175)+0.333*A(176)&
               &+0.201*A(195)+0.201*A(199)+0.006*A(215)
  Vdot(90) = -A(102)+A(103)-A(104)-A(105)-A(106)-A(107)-A(108)-A(110)-A(111)-A(112)-A(113)-2*A(114)+0.5*A(162)+0.5&
               &*A(164)+0.33*A(166)+0.3*A(170)+0.289*A(171)+0.15*A(173)+0.192*A(191)+0.24*A(193)
  Vdot(91) = A(23)+A(26)+A(29)+A(30)-A(31)-A(32)+A(33)+0.61*A(34)-A(36)-2*A(37)-2*A(38)-A(39)+A(42)-A(43)+A(44)+A(45)&
               &+A(46)-A(47)+A(48)+2*A(50)+A(51)-A(52)+A(53)+A(54)+A(55)-A(63)+A(64)+A(65)+A(66)+A(68)-A(72)-A(82)-A(93)&
               &-A(105)-A(118)-A(121)+2*A(123)+A(125)-A(126)+A(127)+A(128)+A(129)+A(131)+A(134)+A(140)+0.95*A(141)+A(143)&
               &+A(145)+2*A(146)+0.63*A(148)+0.63*A(149)+A(150)+0.008*A(163)+0.34*A(166)+0.064*A(168)+0.4*A(172)+1.233&
               &*A(174)+0.379*A(175)+0.113*A(177)+0.341*A(178)+1.5*A(180)+0.5*A(182)+0.5*A(184)+0.12*A(187)+0.5*A(189)+0.033&
               &*A(195)+0.033*A(199)+0.297*A(204)+0.224*A(208)+0.187*A(209)+0.056*A(211)+0.003*A(215)+0.013*A(217)+1.5&
               &*A(218)+0.06*A(220)-A(224)+0.5*A(258)
  Vdot(92) = -A(1)+A(4)-A(5)-A(6)+A(7)-A(8)+2*A(9)+2*A(10)-A(11)+A(12)+A(16)+A(23)+A(24)-A(25)+A(26)+A(28)+A(31)-A(32)&
               &+A(33)+0.61*A(34)+A(35)+0.8*A(39)+2*A(40)+A(46)+A(48)+A(51)+A(53)+A(56)+A(58)+A(65)-A(69)+A(70)+A(71)+A(73)&
               &-A(79)+A(80)+A(81)+A(83)-A(90)+A(91)+A(92)+A(94)-A(102)+A(103)+A(104)+A(106)-A(115)-A(117)-A(120)+A(128)&
               &+0.338*A(177)+A(178)+0.187*A(192)+0.474*A(196)+0.474*A(200)+0.391*A(216)
  Vdot(93) = -A(46)-A(47)-A(48)-2*A(49)-2*A(50)-A(54)-A(64)+A(71)+A(73)-A(74)+2*A(78)-A(84)+A(88)-A(95)+A(99)-A(107)&
               &+A(111)+A(116)+A(131)+A(137)+0.65*A(142)+0.3*A(170)+A(185)+0.3*A(189)+0.25*A(193)+0.011*A(206)+0.076*A(211)&
               &+0.197*A(215)+0.03*A(216)+0.26*A(220)
  Vdot(94) = -A(51)-A(52)-A(53)-A(54)-2*A(55)-A(66)-A(75)+A(81)+A(83)-A(85)+A(88)+2*A(89)-A(96)+A(100)-A(108)+A(112)&
               &+0.034*A(133)+A(134)+0.37*A(138)+A(139)+0.05*A(141)+0.34*A(144)+0.76*A(154)+0.76*A(156)+0.5*A(162)+0.1&
               &*A(163)+0.5*A(164)+0.33*A(166)+0.3*A(167)+0.05*A(168)+0.67*A(171)+0.048*A(172)+0.799*A(173)+0.473*A(175)&
               &+0.96*A(176)+0.376*A(177)+0.564*A(178)+A(179)+A(182)+A(184)+A(186)+A(188)+0.2*A(189)+0.907*A(190)+0.066&
               &*A(191)+0.749*A(192)+0.75*A(194)+0.031*A(195)+0.276*A(196)+0.75*A(198)+0.031*A(199)+0.276*A(200)+A(202)&
               &+0.965*A(203)+0.1*A(204)+0.695*A(205)+0.835*A(206)+0.653*A(207)+0.765*A(208)+0.804*A(209)+0.91*A(210)+0.022&
               &*A(211)+0.824*A(212)+0.918*A(214)+0.033*A(215)+0.442*A(216)+0.012*A(217)+0.984*A(219)+0.949*A(221)
  Vdot(95) = -A(90)+A(91)-A(92)-A(93)-A(94)-A(95)-A(96)-A(98)-A(99)-A(100)-2*A(101)-A(113)+A(159)+A(161)
  Vdot(96) = -A(69)+A(70)-A(71)-A(72)-A(73)-A(74)-A(75)-A(77)-2*A(78)-A(88)-A(99)+A(104)+A(106)+A(112)+A(113)+2*A(114)&
               &+A(130)+A(132)+A(136)+A(137)+0.492*A(138)+A(139)+A(150)+A(151)+A(152)+2*A(153)+0.67*A(166)+0.675*A(167)&
               &+0.467*A(174)+0.029*A(175)+0.667*A(176)+A(181)+0.5*A(182)+A(183)+0.5*A(184)+0.123*A(195)+0.123*A(199)+0.011&
               &*A(206)+0.137*A(215)
  Vdot(97) = A(1)-A(4)+A(5)-A(7)-A(9)-2*A(10)+A(14)+A(15)-A(21)+A(22)-A(31)-A(46)-A(51)-A(56)-A(62)-A(71)-A(81)-A(92)&
               &-A(104)-A(128)
  Vdot(98) = 2*A(19)-A(21)+A(22)-A(24)-A(25)-A(26)-A(27)+A(28)-A(29)-A(30)+A(31)+0.39*A(34)-A(35)+A(36)+0.8*A(39)+2&
               &*A(41)-A(42)-A(43)-A(44)-A(45)-A(125)-A(130)-A(133)-A(136)-A(138)-A(140)-A(141)-0.65*A(142)+A(143)-0.34&
               &*A(144)+A(145)-A(148)-A(151)-A(154)-A(156)-A(159)-A(162)+0.208*A(163)+0.33*A(166)-A(167)+0.164*A(168)-A(171)&
               &+0.285*A(172)-A(175)-A(177)-A(179)+0.5*A(180)-A(181)-A(183)-A(185)-A(186)+0.12*A(187)-A(190)+0.266*A(191)&
               &-A(194)+0.567*A(195)-A(198)+0.567*A(199)-A(202)-A(203)-0.397*A(204)-A(205)-A(206)-A(207)-A(208)-A(209)&
               &-A(210)+0.155*A(211)-A(214)+0.378*A(215)+0.5*A(218)-A(219)+0.32*A(220)-A(245)-A(247)-A(249)-A(251)-A(253)&
               &-A(255)-A(257)-A(258)
  Vdot(99) = A(6)+A(8)-A(9)-A(11)+A(12)-A(14)-A(15)-A(16)-A(26)+A(27)+0.39*A(34)-A(39)-2*A(40)-A(48)-A(53)-A(58)-A(65)&
               &-A(73)-A(83)-A(94)-A(106)-A(129)-A(132)-A(135)-A(149)-A(152)-A(155)-A(157)-A(158)-A(161)-A(164)-A(173)&
               &-A(188)-A(192)-A(196)-A(200)-A(212)-A(216)-A(221)-A(259)
      
END SUBROUTINE saprc99_mosaic_8bin_vbs2_Fun
















SUBROUTINE saprc99_mosaic_8bin_vbs2_IRRFun ( V, F, RCT, IRR )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: IRR(NREACT)



  IRR(1) = RCT(1)*V(92)
  IRR(2) = RCT(2)*V(83)*F(2)
  IRR(3) = RCT(3)*V(83)*V(87)
  IRR(4) = RCT(4)*V(83)*V(97)*F(2)
  IRR(5) = RCT(5)*V(83)*V(92)
  IRR(6) = RCT(6)*V(83)*V(92)
  IRR(7) = RCT(7)*V(87)*V(97)
  IRR(8) = RCT(8)*V(87)*V(92)
  IRR(9) = RCT(9)*V(97)*V(99)
  IRR(10) = RCT(10)*V(97)*V(97)*F(2)
  IRR(11) = RCT(11)*V(92)*V(99)
  IRR(12) = RCT(12)*V(39)
  IRR(13) = RCT(13)*V(39)*F(1)
  IRR(14) = RCT(14)*V(92)*V(99)
  IRR(15) = RCT(15)*V(99)
  IRR(16) = RCT(16)*V(99)
  IRR(17) = RCT(17)*V(87)
  IRR(18) = RCT(18)*V(87)
  IRR(19) = RCT(19)*V(24)*F(1)
  IRR(20) = RCT(20)*V(24)*F(2)
  IRR(21) = RCT(21)*V(97)*V(98)
  IRR(22) = RCT(22)*V(42)
  IRR(23) = RCT(23)*V(42)
  IRR(24) = RCT(24)*V(42)*V(98)
  IRR(25) = RCT(25)*V(92)*V(98)
  IRR(26) = RCT(26)*V(98)*V(99)
  IRR(27) = RCT(27)*V(65)*V(98)
  IRR(28) = RCT(28)*V(65)
  IRR(29) = RCT(29)*V(64)*V(98)
  IRR(30) = RCT(30)*V(87)*V(98)
  IRR(31) = RCT(31)*V(91)*V(97)
  IRR(32) = RCT(32)*V(91)*V(92)
  IRR(33) = RCT(33)*V(50)
  IRR(34) = RCT(34)*V(50)
  IRR(35) = RCT(35)*V(50)*V(98)
  IRR(36) = RCT(36)*V(87)*V(91)
  IRR(37) = RCT(37)*V(91)*V(91)
  IRR(38) = RCT(38)*V(91)*V(91)*F(1)
  IRR(39) = RCT(39)*V(91)*V(99)
  IRR(40) = RCT(40)*V(99)*V(99)
  IRR(41) = RCT(41)*V(35)
  IRR(42) = RCT(42)*V(35)*V(98)
  IRR(43) = RCT(43)*V(91)*V(98)
  IRR(44) = RCT(44)*V(40)*V(98)
  IRR(45) = RCT(45)*V(98)*F(2)
  IRR(46) = RCT(46)*V(93)*V(97)
  IRR(47) = RCT(47)*V(91)*V(93)
  IRR(48) = RCT(48)*V(93)*V(99)
  IRR(49) = RCT(49)*V(93)*V(93)
  IRR(50) = RCT(50)*V(93)*V(93)
  IRR(51) = RCT(51)*V(94)*V(97)
  IRR(52) = RCT(52)*V(91)*V(94)
  IRR(53) = RCT(53)*V(94)*V(99)
  IRR(54) = RCT(54)*V(93)*V(94)
  IRR(55) = RCT(55)*V(94)*V(94)
  IRR(56) = RCT(56)*V(72)*V(97)
  IRR(57) = RCT(57)*V(72)*V(91)
  IRR(58) = RCT(58)*V(72)*V(99)
  IRR(59) = RCT(59)*V(72)*V(93)
  IRR(60) = RCT(60)*V(72)*V(94)
  IRR(61) = RCT(61)*V(72)*V(72)
  IRR(62) = RCT(62)*V(88)*V(97)
  IRR(63) = RCT(63)*V(88)*V(91)
  IRR(64) = RCT(64)*V(88)*V(93)
  IRR(65) = RCT(65)*V(88)*V(99)
  IRR(66) = RCT(66)*V(88)*V(94)
  IRR(67) = RCT(67)*V(72)*V(88)
  IRR(68) = RCT(68)*V(88)*V(88)
  IRR(69) = RCT(69)*V(92)*V(96)
  IRR(70) = RCT(70)*V(31)
  IRR(71) = RCT(71)*V(96)*V(97)
  IRR(72) = RCT(72)*V(91)*V(96)
  IRR(73) = RCT(73)*V(96)*V(99)
  IRR(74) = RCT(74)*V(93)*V(96)
  IRR(75) = RCT(75)*V(94)*V(96)
  IRR(76) = RCT(76)*V(72)*V(96)
  IRR(77) = RCT(77)*V(88)*V(96)
  IRR(78) = RCT(78)*V(96)*V(96)
  IRR(79) = RCT(79)*V(89)*V(92)
  IRR(80) = RCT(80)*V(32)
  IRR(81) = RCT(81)*V(89)*V(97)
  IRR(82) = RCT(82)*V(89)*V(91)
  IRR(83) = RCT(83)*V(89)*V(99)
  IRR(84) = RCT(84)*V(89)*V(93)
  IRR(85) = RCT(85)*V(89)*V(94)
  IRR(86) = RCT(86)*V(72)*V(89)
  IRR(87) = RCT(87)*V(88)*V(89)
  IRR(88) = RCT(88)*V(89)*V(96)
  IRR(89) = RCT(89)*V(89)*V(89)
  IRR(90) = RCT(90)*V(92)*V(95)
  IRR(91) = RCT(91)*V(33)
  IRR(92) = RCT(92)*V(95)*V(97)
  IRR(93) = RCT(93)*V(91)*V(95)
  IRR(94) = RCT(94)*V(95)*V(99)
  IRR(95) = RCT(95)*V(93)*V(95)
  IRR(96) = RCT(96)*V(94)*V(95)
  IRR(97) = RCT(97)*V(72)*V(95)
  IRR(98) = RCT(98)*V(88)*V(95)
  IRR(99) = RCT(99)*V(95)*V(96)
  IRR(100) = RCT(100)*V(89)*V(95)
  IRR(101) = RCT(101)*V(95)*V(95)
  IRR(102) = RCT(102)*V(90)*V(92)
  IRR(103) = RCT(103)*V(34)
  IRR(104) = RCT(104)*V(90)*V(97)
  IRR(105) = RCT(105)*V(90)*V(91)
  IRR(106) = RCT(106)*V(90)*V(99)
  IRR(107) = RCT(107)*V(90)*V(93)
  IRR(108) = RCT(108)*V(90)*V(94)
  IRR(109) = RCT(109)*V(72)*V(90)
  IRR(110) = RCT(110)*V(88)*V(90)
  IRR(111) = RCT(111)*V(90)*V(96)
  IRR(112) = RCT(112)*V(89)*V(90)
  IRR(113) = RCT(113)*V(90)*V(95)
  IRR(114) = RCT(114)*V(90)*V(90)
  IRR(115) = RCT(115)*V(45)*V(92)
  IRR(116) = RCT(116)*V(45)
  IRR(117) = RCT(117)*V(69)*V(92)
  IRR(118) = RCT(118)*V(69)*V(91)
  IRR(119) = RCT(119)*V(69)
  IRR(120) = RCT(120)*V(48)*V(92)
  IRR(121) = RCT(121)*V(48)*V(91)
  IRR(122) = RCT(122)*V(48)
  IRR(123) = RCT(123)*V(81)
  IRR(124) = RCT(124)*V(81)
  IRR(125) = RCT(125)*V(81)*V(98)
  IRR(126) = RCT(126)*V(81)*V(91)
  IRR(127) = RCT(127)*V(47)
  IRR(128) = RCT(128)*V(47)*V(97)
  IRR(129) = RCT(129)*V(81)*V(99)
  IRR(130) = RCT(130)*V(78)*V(98)
  IRR(131) = RCT(131)*V(78)
  IRR(132) = RCT(132)*V(78)*V(99)
  IRR(133) = RCT(133)*V(84)*V(98)
  IRR(134) = RCT(134)*V(84)
  IRR(135) = RCT(135)*V(84)*V(99)
  IRR(136) = RCT(136)*V(67)*V(98)
  IRR(137) = RCT(137)*V(67)
  IRR(138) = RCT(138)*V(85)*V(98)
  IRR(139) = RCT(139)*V(85)
  IRR(140) = RCT(140)*V(51)*V(98)
  IRR(141) = RCT(141)*V(38)*V(98)
  IRR(142) = RCT(142)*V(46)*V(98)
  IRR(143) = RCT(143)*V(46)
  IRR(144) = RCT(144)*V(59)*V(98)
  IRR(145) = RCT(145)*V(59)
  IRR(146) = RCT(146)*V(76)
  IRR(147) = RCT(147)*V(76)
  IRR(148) = RCT(148)*V(76)*V(98)
  IRR(149) = RCT(149)*V(76)*V(99)
  IRR(150) = RCT(150)*V(63)
  IRR(151) = RCT(151)*V(63)*V(98)
  IRR(152) = RCT(152)*V(63)*V(99)
  IRR(153) = RCT(153)*V(37)
  IRR(154) = RCT(154)*V(62)*V(98)
  IRR(155) = RCT(155)*V(62)*V(99)
  IRR(156) = RCT(156)*V(56)*V(98)
  IRR(157) = RCT(157)*V(56)*V(99)
  IRR(158) = RCT(158)*V(60)*V(99)
  IRR(159) = RCT(159)*V(61)*V(98)
  IRR(160) = RCT(160)*V(61)
  IRR(161) = RCT(161)*V(61)*V(99)
  IRR(162) = RCT(162)*V(73)*V(98)
  IRR(163) = RCT(163)*V(73)*V(87)
  IRR(164) = RCT(164)*V(73)*V(99)
  IRR(165) = RCT(165)*V(73)*V(83)
  IRR(166) = RCT(166)*V(73)
  IRR(167) = RCT(167)*V(80)*V(98)
  IRR(168) = RCT(168)*V(80)*V(87)
  IRR(169) = RCT(169)*V(80)*V(83)
  IRR(170) = RCT(170)*V(80)
  IRR(171) = RCT(171)*V(77)*V(98)
  IRR(172) = RCT(172)*V(77)*V(87)
  IRR(173) = RCT(173)*V(77)*V(99)
  IRR(174) = RCT(174)*V(77)
  IRR(175) = RCT(175)*V(86)*V(98)
  IRR(176) = RCT(176)*V(86)
  IRR(177) = RCT(177)*V(82)*V(98)
  IRR(178) = RCT(178)*V(82)
  IRR(179) = RCT(179)*V(57)*V(98)
  IRR(180) = RCT(180)*V(57)*V(87)
  IRR(181) = RCT(181)*V(54)*V(98)
  IRR(182) = RCT(182)*V(54)
  IRR(183) = RCT(183)*V(55)*V(98)
  IRR(184) = RCT(184)*V(55)
  IRR(185) = RCT(185)*V(29)*V(98)
  IRR(186) = RCT(186)*V(66)*V(98)
  IRR(187) = RCT(187)*V(66)*V(87)
  IRR(188) = RCT(188)*V(66)*V(99)
  IRR(189) = RCT(189)*V(66)*V(83)
  IRR(190) = RCT(190)*V(70)*V(98)
  IRR(191) = RCT(191)*V(70)*V(87)
  IRR(192) = RCT(192)*V(70)*V(99)
  IRR(193) = RCT(193)*V(70)*V(83)
  IRR(194) = RCT(194)*V(75)*V(98)
  IRR(195) = RCT(195)*V(75)*V(87)
  IRR(196) = RCT(196)*V(75)*V(99)
  IRR(197) = RCT(197)*V(75)*V(83)
  IRR(198) = RCT(198)*V(74)*V(98)
  IRR(199) = RCT(199)*V(74)*V(87)
  IRR(200) = RCT(200)*V(74)*V(99)
  IRR(201) = RCT(201)*V(74)*V(83)
  IRR(202) = RCT(202)*V(30)*V(98)
  IRR(203) = RCT(203)*V(36)*V(98)
  IRR(204) = RCT(204)*V(58)*V(98)
  IRR(205) = RCT(205)*V(44)*V(98)
  IRR(206) = RCT(206)*V(52)*V(98)
  IRR(207) = RCT(207)*V(43)*V(98)
  IRR(208) = RCT(208)*V(53)*V(98)
  IRR(209) = RCT(209)*V(49)*V(98)
  IRR(210) = RCT(210)*V(71)*V(98)
  IRR(211) = RCT(211)*V(71)*V(87)
  IRR(212) = RCT(212)*V(71)*V(99)
  IRR(213) = RCT(213)*V(71)*V(83)
  IRR(214) = RCT(214)*V(79)*V(98)
  IRR(215) = RCT(215)*V(79)*V(87)
  IRR(216) = RCT(216)*V(79)*V(99)
  IRR(217) = RCT(217)*V(79)*V(83)
  IRR(218) = RCT(218)*V(58)*V(87)
  IRR(219) = RCT(219)*V(68)*V(98)
  IRR(220) = RCT(220)*V(68)*V(87)
  IRR(221) = RCT(221)*V(68)*V(99)
  IRR(222) = RCT(222)*V(68)*V(83)
  IRR(223) = RCT(223)*V(40)
  IRR(224) = RCT(224)*V(91)
  IRR(225) = RCT(225)*V(40)
  IRR(226) = RCT(226)*V(1)
  IRR(227) = RCT(227)*V(65)
  IRR(228) = RCT(228)*V(35)
  IRR(229) = RCT(229)*V(2)
  IRR(230) = RCT(230)*V(52)*V(98)
  IRR(231) = RCT(231)*V(43)*V(98)
  IRR(232) = RCT(232)*V(71)*V(98)
  IRR(233) = RCT(233)*V(79)*V(98)
  IRR(234) = RCT(234)*V(53)*V(98)
  IRR(235) = RCT(235)*V(49)*V(98)
  IRR(236) = RCT(236)*V(70)*V(98)
  IRR(237) = RCT(237)*V(70)*V(87)
  IRR(238) = RCT(238)*V(70)*V(99)
  IRR(239) = RCT(239)*V(75)*V(98)
  IRR(240) = RCT(240)*V(75)*V(87)
  IRR(241) = RCT(241)*V(75)*V(99)
  IRR(242) = RCT(242)*V(74)*V(98)
  IRR(243) = RCT(243)*V(74)*V(87)
  IRR(244) = RCT(244)*V(74)*V(99)
  IRR(245) = RCT(245)*V(3)*V(98)
  IRR(246) = RCT(246)*V(4)*V(98)
  IRR(247) = RCT(247)*V(26)*V(98)
  IRR(248) = RCT(248)*V(20)*V(98)
  IRR(249) = RCT(249)*V(5)*V(98)
  IRR(250) = RCT(250)*V(6)*V(98)
  IRR(251) = RCT(251)*V(28)*V(98)
  IRR(252) = RCT(252)*V(22)*V(98)
  IRR(253) = RCT(253)*V(25)*V(98)
  IRR(254) = RCT(254)*V(21)*V(98)
  IRR(255) = RCT(255)*V(27)*V(98)
  IRR(256) = RCT(256)*V(23)*V(98)
  IRR(257) = RCT(257)*V(41)*V(98)
  IRR(258) = RCT(258)*V(41)*V(98)
  IRR(259) = RCT(259)*V(41)*V(99)
      
END SUBROUTINE saprc99_mosaic_8bin_vbs2_IRRFun
















SUBROUTINE saprc99_mosaic_8bin_vbs2_Jac_SP ( V, F, RCT, JVS )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: JVS(LU_NONZERO)




  REAL(kind=dp) :: B(461)


  B(1) = RCT(1)

  B(2) = RCT(2)*F(2)

  B(4) = RCT(3)*V(87)

  B(5) = RCT(3)*V(83)

  B(6) = RCT(4)*V(97)*F(2)

  B(7) = RCT(4)*V(83)*F(2)

  B(9) = RCT(5)*V(92)

  B(10) = RCT(5)*V(83)

  B(11) = RCT(6)*V(92)

  B(12) = RCT(6)*V(83)

  B(13) = RCT(7)*V(97)

  B(14) = RCT(7)*V(87)

  B(15) = RCT(8)*V(92)

  B(16) = RCT(8)*V(87)

  B(17) = RCT(9)*V(99)

  B(18) = RCT(9)*V(97)

  B(19) = RCT(10)*2*V(97)*F(2)

  B(21) = RCT(11)*V(99)

  B(22) = RCT(11)*V(92)

  B(23) = RCT(12)

  B(24) = RCT(13)*F(1)

  B(26) = RCT(14)*V(99)

  B(27) = RCT(14)*V(92)

  B(28) = RCT(15)

  B(29) = RCT(16)

  B(30) = RCT(17)

  B(31) = RCT(18)

  B(32) = RCT(19)*F(1)

  B(34) = RCT(20)*F(2)

  B(36) = RCT(21)*V(98)

  B(37) = RCT(21)*V(97)

  B(38) = RCT(22)

  B(39) = RCT(23)

  B(40) = RCT(24)*V(98)

  B(41) = RCT(24)*V(42)

  B(42) = RCT(25)*V(98)

  B(43) = RCT(25)*V(92)

  B(44) = RCT(26)*V(99)

  B(45) = RCT(26)*V(98)

  B(46) = RCT(27)*V(98)

  B(47) = RCT(27)*V(65)

  B(48) = RCT(28)

  B(49) = RCT(29)*V(98)

  B(50) = RCT(29)*V(64)

  B(51) = RCT(30)*V(98)

  B(52) = RCT(30)*V(87)

  B(53) = RCT(31)*V(97)

  B(54) = RCT(31)*V(91)

  B(55) = RCT(32)*V(92)

  B(56) = RCT(32)*V(91)

  B(57) = RCT(33)

  B(58) = RCT(34)

  B(59) = RCT(35)*V(98)

  B(60) = RCT(35)*V(50)

  B(61) = RCT(36)*V(91)

  B(62) = RCT(36)*V(87)

  B(63) = RCT(37)*2*V(91)

  B(64) = RCT(38)*2*V(91)*F(1)

  B(66) = RCT(39)*V(99)

  B(67) = RCT(39)*V(91)

  B(68) = RCT(40)*2*V(99)

  B(69) = RCT(41)

  B(70) = RCT(42)*V(98)

  B(71) = RCT(42)*V(35)

  B(72) = RCT(43)*V(98)

  B(73) = RCT(43)*V(91)

  B(74) = RCT(44)*V(98)

  B(75) = RCT(44)*V(40)

  B(76) = RCT(45)*F(2)

  B(78) = RCT(46)*V(97)

  B(79) = RCT(46)*V(93)

  B(80) = RCT(47)*V(93)

  B(81) = RCT(47)*V(91)

  B(82) = RCT(48)*V(99)

  B(83) = RCT(48)*V(93)

  B(84) = RCT(49)*2*V(93)

  B(85) = RCT(50)*2*V(93)

  B(86) = RCT(51)*V(97)

  B(87) = RCT(51)*V(94)

  B(88) = RCT(52)*V(94)

  B(89) = RCT(52)*V(91)

  B(90) = RCT(53)*V(99)

  B(91) = RCT(53)*V(94)

  B(92) = RCT(54)*V(94)

  B(93) = RCT(54)*V(93)

  B(94) = RCT(55)*2*V(94)

  B(95) = RCT(56)*V(97)

  B(96) = RCT(56)*V(72)

  B(97) = RCT(57)*V(91)

  B(98) = RCT(57)*V(72)

  B(99) = RCT(58)*V(99)

  B(100) = RCT(58)*V(72)

  B(101) = RCT(59)*V(93)

  B(102) = RCT(59)*V(72)

  B(103) = RCT(60)*V(94)

  B(104) = RCT(60)*V(72)

  B(105) = RCT(61)*2*V(72)

  B(106) = RCT(62)*V(97)

  B(107) = RCT(62)*V(88)

  B(108) = RCT(63)*V(91)

  B(109) = RCT(63)*V(88)

  B(110) = RCT(64)*V(93)

  B(111) = RCT(64)*V(88)

  B(112) = RCT(65)*V(99)

  B(113) = RCT(65)*V(88)

  B(114) = RCT(66)*V(94)

  B(115) = RCT(66)*V(88)

  B(116) = RCT(67)*V(88)

  B(117) = RCT(67)*V(72)

  B(118) = RCT(68)*2*V(88)

  B(119) = RCT(69)*V(96)

  B(120) = RCT(69)*V(92)

  B(121) = RCT(70)

  B(122) = RCT(71)*V(97)

  B(123) = RCT(71)*V(96)

  B(124) = RCT(72)*V(96)

  B(125) = RCT(72)*V(91)

  B(126) = RCT(73)*V(99)

  B(127) = RCT(73)*V(96)

  B(128) = RCT(74)*V(96)

  B(129) = RCT(74)*V(93)

  B(130) = RCT(75)*V(96)

  B(131) = RCT(75)*V(94)

  B(132) = RCT(76)*V(96)

  B(133) = RCT(76)*V(72)

  B(134) = RCT(77)*V(96)

  B(135) = RCT(77)*V(88)

  B(136) = RCT(78)*2*V(96)

  B(137) = RCT(79)*V(92)

  B(138) = RCT(79)*V(89)

  B(139) = RCT(80)

  B(140) = RCT(81)*V(97)

  B(141) = RCT(81)*V(89)

  B(142) = RCT(82)*V(91)

  B(143) = RCT(82)*V(89)

  B(144) = RCT(83)*V(99)

  B(145) = RCT(83)*V(89)

  B(146) = RCT(84)*V(93)

  B(147) = RCT(84)*V(89)

  B(148) = RCT(85)*V(94)

  B(149) = RCT(85)*V(89)

  B(150) = RCT(86)*V(89)

  B(151) = RCT(86)*V(72)

  B(152) = RCT(87)*V(89)

  B(153) = RCT(87)*V(88)

  B(154) = RCT(88)*V(96)

  B(155) = RCT(88)*V(89)

  B(156) = RCT(89)*2*V(89)

  B(157) = RCT(90)*V(95)

  B(158) = RCT(90)*V(92)

  B(159) = RCT(91)

  B(160) = RCT(92)*V(97)

  B(161) = RCT(92)*V(95)

  B(162) = RCT(93)*V(95)

  B(163) = RCT(93)*V(91)

  B(164) = RCT(94)*V(99)

  B(165) = RCT(94)*V(95)

  B(166) = RCT(95)*V(95)

  B(167) = RCT(95)*V(93)

  B(168) = RCT(96)*V(95)

  B(169) = RCT(96)*V(94)

  B(170) = RCT(97)*V(95)

  B(171) = RCT(97)*V(72)

  B(172) = RCT(98)*V(95)

  B(173) = RCT(98)*V(88)

  B(174) = RCT(99)*V(96)

  B(175) = RCT(99)*V(95)

  B(176) = RCT(100)*V(95)

  B(177) = RCT(100)*V(89)

  B(178) = RCT(101)*2*V(95)

  B(179) = RCT(102)*V(92)

  B(180) = RCT(102)*V(90)

  B(181) = RCT(103)

  B(182) = RCT(104)*V(97)

  B(183) = RCT(104)*V(90)

  B(184) = RCT(105)*V(91)

  B(185) = RCT(105)*V(90)

  B(186) = RCT(106)*V(99)

  B(187) = RCT(106)*V(90)

  B(188) = RCT(107)*V(93)

  B(189) = RCT(107)*V(90)

  B(190) = RCT(108)*V(94)

  B(191) = RCT(108)*V(90)

  B(192) = RCT(109)*V(90)

  B(193) = RCT(109)*V(72)

  B(194) = RCT(110)*V(90)

  B(195) = RCT(110)*V(88)

  B(196) = RCT(111)*V(96)

  B(197) = RCT(111)*V(90)

  B(198) = RCT(112)*V(90)

  B(199) = RCT(112)*V(89)

  B(200) = RCT(113)*V(95)

  B(201) = RCT(113)*V(90)

  B(202) = RCT(114)*2*V(90)

  B(203) = RCT(115)*V(92)

  B(204) = RCT(115)*V(45)

  B(205) = RCT(116)

  B(206) = RCT(117)*V(92)

  B(207) = RCT(117)*V(69)

  B(208) = RCT(118)*V(91)

  B(209) = RCT(118)*V(69)

  B(210) = RCT(119)

  B(211) = RCT(120)*V(92)

  B(212) = RCT(120)*V(48)

  B(213) = RCT(121)*V(91)

  B(214) = RCT(121)*V(48)

  B(215) = RCT(122)

  B(216) = RCT(123)

  B(217) = RCT(124)

  B(218) = RCT(125)*V(98)

  B(219) = RCT(125)*V(81)

  B(220) = RCT(126)*V(91)

  B(221) = RCT(126)*V(81)

  B(222) = RCT(127)

  B(223) = RCT(128)*V(97)

  B(224) = RCT(128)*V(47)

  B(225) = RCT(129)*V(99)

  B(226) = RCT(129)*V(81)

  B(227) = RCT(130)*V(98)

  B(228) = RCT(130)*V(78)

  B(229) = RCT(131)

  B(230) = RCT(132)*V(99)

  B(231) = RCT(132)*V(78)

  B(232) = RCT(133)*V(98)

  B(233) = RCT(133)*V(84)

  B(234) = RCT(134)

  B(235) = RCT(135)*V(99)

  B(236) = RCT(135)*V(84)

  B(237) = RCT(136)*V(98)

  B(238) = RCT(136)*V(67)

  B(239) = RCT(137)

  B(240) = RCT(138)*V(98)

  B(241) = RCT(138)*V(85)

  B(242) = RCT(139)

  B(243) = RCT(140)*V(98)

  B(244) = RCT(140)*V(51)

  B(245) = RCT(141)*V(98)

  B(246) = RCT(141)*V(38)

  B(247) = RCT(142)*V(98)

  B(248) = RCT(142)*V(46)

  B(249) = RCT(143)

  B(250) = RCT(144)*V(98)

  B(251) = RCT(144)*V(59)

  B(252) = RCT(145)

  B(253) = RCT(146)

  B(254) = RCT(147)

  B(255) = RCT(148)*V(98)

  B(256) = RCT(148)*V(76)

  B(257) = RCT(149)*V(99)

  B(258) = RCT(149)*V(76)

  B(259) = RCT(150)

  B(260) = RCT(151)*V(98)

  B(261) = RCT(151)*V(63)

  B(262) = RCT(152)*V(99)

  B(263) = RCT(152)*V(63)

  B(264) = RCT(153)

  B(265) = RCT(154)*V(98)

  B(266) = RCT(154)*V(62)

  B(267) = RCT(155)*V(99)

  B(268) = RCT(155)*V(62)

  B(269) = RCT(156)*V(98)

  B(270) = RCT(156)*V(56)

  B(271) = RCT(157)*V(99)

  B(272) = RCT(157)*V(56)

  B(273) = RCT(158)*V(99)

  B(274) = RCT(158)*V(60)

  B(275) = RCT(159)*V(98)

  B(276) = RCT(159)*V(61)

  B(277) = RCT(160)

  B(278) = RCT(161)*V(99)

  B(279) = RCT(161)*V(61)

  B(280) = RCT(162)*V(98)

  B(281) = RCT(162)*V(73)

  B(282) = RCT(163)*V(87)

  B(283) = RCT(163)*V(73)

  B(284) = RCT(164)*V(99)

  B(285) = RCT(164)*V(73)

  B(286) = RCT(165)*V(83)

  B(287) = RCT(165)*V(73)

  B(288) = RCT(166)

  B(289) = RCT(167)*V(98)

  B(290) = RCT(167)*V(80)

  B(291) = RCT(168)*V(87)

  B(292) = RCT(168)*V(80)

  B(293) = RCT(169)*V(83)

  B(294) = RCT(169)*V(80)

  B(295) = RCT(170)

  B(296) = RCT(171)*V(98)

  B(297) = RCT(171)*V(77)

  B(298) = RCT(172)*V(87)

  B(299) = RCT(172)*V(77)

  B(300) = RCT(173)*V(99)

  B(301) = RCT(173)*V(77)

  B(302) = RCT(174)

  B(303) = RCT(175)*V(98)

  B(304) = RCT(175)*V(86)

  B(305) = RCT(176)

  B(306) = RCT(177)*V(98)

  B(307) = RCT(177)*V(82)

  B(308) = RCT(178)

  B(309) = RCT(179)*V(98)

  B(310) = RCT(179)*V(57)

  B(311) = RCT(180)*V(87)

  B(312) = RCT(180)*V(57)

  B(313) = RCT(181)*V(98)

  B(314) = RCT(181)*V(54)

  B(315) = RCT(182)

  B(316) = RCT(183)*V(98)

  B(317) = RCT(183)*V(55)

  B(318) = RCT(184)

  B(319) = RCT(185)*V(98)

  B(320) = RCT(185)*V(29)

  B(321) = RCT(186)*V(98)

  B(322) = RCT(186)*V(66)

  B(323) = RCT(187)*V(87)

  B(324) = RCT(187)*V(66)

  B(325) = RCT(188)*V(99)

  B(326) = RCT(188)*V(66)

  B(327) = RCT(189)*V(83)

  B(328) = RCT(189)*V(66)

  B(329) = RCT(190)*V(98)

  B(330) = RCT(190)*V(70)

  B(331) = RCT(191)*V(87)

  B(332) = RCT(191)*V(70)

  B(333) = RCT(192)*V(99)

  B(334) = RCT(192)*V(70)

  B(335) = RCT(193)*V(83)

  B(336) = RCT(193)*V(70)

  B(337) = RCT(194)*V(98)

  B(338) = RCT(194)*V(75)

  B(339) = RCT(195)*V(87)

  B(340) = RCT(195)*V(75)

  B(341) = RCT(196)*V(99)

  B(342) = RCT(196)*V(75)

  B(343) = RCT(197)*V(83)

  B(344) = RCT(197)*V(75)

  B(345) = RCT(198)*V(98)

  B(346) = RCT(198)*V(74)

  B(347) = RCT(199)*V(87)

  B(348) = RCT(199)*V(74)

  B(349) = RCT(200)*V(99)

  B(350) = RCT(200)*V(74)

  B(351) = RCT(201)*V(83)

  B(352) = RCT(201)*V(74)

  B(353) = RCT(202)*V(98)

  B(354) = RCT(202)*V(30)

  B(355) = RCT(203)*V(98)

  B(356) = RCT(203)*V(36)

  B(357) = RCT(204)*V(98)

  B(358) = RCT(204)*V(58)

  B(359) = RCT(205)*V(98)

  B(360) = RCT(205)*V(44)

  B(361) = RCT(206)*V(98)

  B(362) = RCT(206)*V(52)

  B(363) = RCT(207)*V(98)

  B(364) = RCT(207)*V(43)

  B(365) = RCT(208)*V(98)

  B(366) = RCT(208)*V(53)

  B(367) = RCT(209)*V(98)

  B(368) = RCT(209)*V(49)

  B(369) = RCT(210)*V(98)

  B(370) = RCT(210)*V(71)

  B(371) = RCT(211)*V(87)

  B(372) = RCT(211)*V(71)

  B(373) = RCT(212)*V(99)

  B(374) = RCT(212)*V(71)

  B(375) = RCT(213)*V(83)

  B(376) = RCT(213)*V(71)

  B(377) = RCT(214)*V(98)

  B(378) = RCT(214)*V(79)

  B(379) = RCT(215)*V(87)

  B(380) = RCT(215)*V(79)

  B(381) = RCT(216)*V(99)

  B(382) = RCT(216)*V(79)

  B(383) = RCT(217)*V(83)

  B(384) = RCT(217)*V(79)

  B(385) = RCT(218)*V(87)

  B(386) = RCT(218)*V(58)

  B(387) = RCT(219)*V(98)

  B(388) = RCT(219)*V(68)

  B(389) = RCT(220)*V(87)

  B(390) = RCT(220)*V(68)

  B(391) = RCT(221)*V(99)

  B(392) = RCT(221)*V(68)

  B(393) = RCT(222)*V(83)

  B(394) = RCT(222)*V(68)

  B(395) = RCT(223)

  B(396) = RCT(224)

  B(397) = RCT(225)

  B(398) = RCT(226)

  B(399) = RCT(227)

  B(400) = RCT(228)

  B(401) = RCT(229)

  B(402) = RCT(230)*V(98)

  B(403) = RCT(230)*V(52)

  B(404) = RCT(231)*V(98)

  B(405) = RCT(231)*V(43)

  B(406) = RCT(232)*V(98)

  B(407) = RCT(232)*V(71)

  B(408) = RCT(233)*V(98)

  B(409) = RCT(233)*V(79)

  B(410) = RCT(234)*V(98)

  B(411) = RCT(234)*V(53)

  B(412) = RCT(235)*V(98)

  B(413) = RCT(235)*V(49)

  B(414) = RCT(236)*V(98)

  B(415) = RCT(236)*V(70)

  B(416) = RCT(237)*V(87)

  B(417) = RCT(237)*V(70)

  B(418) = RCT(238)*V(99)

  B(419) = RCT(238)*V(70)

  B(420) = RCT(239)*V(98)

  B(421) = RCT(239)*V(75)

  B(422) = RCT(240)*V(87)

  B(423) = RCT(240)*V(75)

  B(424) = RCT(241)*V(99)

  B(425) = RCT(241)*V(75)

  B(426) = RCT(242)*V(98)

  B(427) = RCT(242)*V(74)

  B(428) = RCT(243)*V(87)

  B(429) = RCT(243)*V(74)

  B(430) = RCT(244)*V(99)

  B(431) = RCT(244)*V(74)

  B(432) = RCT(245)*V(98)

  B(433) = RCT(245)*V(3)

  B(434) = RCT(246)*V(98)

  B(435) = RCT(246)*V(4)

  B(436) = RCT(247)*V(98)

  B(437) = RCT(247)*V(26)

  B(438) = RCT(248)*V(98)

  B(439) = RCT(248)*V(20)

  B(440) = RCT(249)*V(98)

  B(441) = RCT(249)*V(5)

  B(442) = RCT(250)*V(98)

  B(443) = RCT(250)*V(6)

  B(444) = RCT(251)*V(98)

  B(445) = RCT(251)*V(28)

  B(446) = RCT(252)*V(98)

  B(447) = RCT(252)*V(22)

  B(448) = RCT(253)*V(98)

  B(449) = RCT(253)*V(25)

  B(450) = RCT(254)*V(98)

  B(451) = RCT(254)*V(21)

  B(452) = RCT(255)*V(98)

  B(453) = RCT(255)*V(27)

  B(454) = RCT(256)*V(98)

  B(455) = RCT(256)*V(23)

  B(456) = RCT(257)*V(98)

  B(457) = RCT(257)*V(41)

  B(458) = RCT(258)*V(98)

  B(459) = RCT(258)*V(41)

  B(460) = RCT(259)*V(99)

  B(461) = RCT(259)*V(41)



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

  JVS(17) = 0.204*B(331)

  JVS(18) = 0.185*B(371)

  JVS(19) = 0.333*B(282)

  JVS(20) = 0.103*B(347)

  JVS(21) = 0.103*B(339)

  JVS(22) = 0.1*B(298)

  JVS(23) = 0.073*B(379)

  JVS(24) = 0.351*B(291)

  JVS(25) = 0.333*B(283)+0.351*B(292)+0.1*B(299)+0.37*B(324)+0.204*B(332)+0.103*B(340)+0.103*B(348)+0.185*B(372)+0.073&
              &*B(380)+0.185*B(390)

  JVS(26) = B(224)

  JVS(27) = 0.297*B(358)

  JVS(28) = 0

  JVS(29) = 0.17*B(389)

  JVS(30) = 0.05*B(371)

  JVS(31) = 0.129*B(379)

  JVS(32) = 0.05*B(372)+0.129*B(380)+0.17*B(390)

  JVS(33) = B(134)

  JVS(34) = 0.25*B(124)

  JVS(35) = B(128)

  JVS(36) = B(130)

  JVS(37) = 0.25*B(125)+B(129)+B(131)+B(135)

  JVS(38) = 0

  JVS(39) = 0.15*B(331)

  JVS(40) = 0.119*B(371)

  JVS(41) = 0.189*B(347)

  JVS(42) = 0.189*B(339)

  JVS(43) = 0.372*B(298)

  JVS(44) = 0.247*B(379)

  JVS(45) = 0.372*B(299)+0.15*B(332)+0.189*B(340)+0.189*B(348)+0.119*B(372)+0.247*B(380)

  JVS(46) = B(152)+B(172)+2*B(194)

  JVS(47) = 0.25*B(142)+B(146)+B(148)+B(153)

  JVS(48) = 0.25*B(184)+B(188)+B(190)+2*B(195)

  JVS(49) = 0.25*B(143)+0.25*B(162)+0.25*B(185)

  JVS(50) = B(147)+B(166)+B(189)

  JVS(51) = B(149)+B(168)+B(191)

  JVS(52) = 0.25*B(163)+B(167)+B(169)+B(173)

  JVS(53) = 0

  JVS(54) = 0.75*B(124)

  JVS(55) = 0.75*B(125)

  JVS(56) = 0

  JVS(57) = 0.75*B(142)

  JVS(58) = 0.75*B(184)

  JVS(59) = 0.75*B(143)+0.75*B(162)+0.75*B(185)

  JVS(60) = 0.75*B(163)

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

  JVS(73) = 0.048*B(388)

  JVS(74) = 2.693*B(392)

  JVS(75) = 0

  JVS(76) = B(95)+B(99)

  JVS(77) = B(106)+B(112)

  JVS(78) = B(78)+B(82)

  JVS(79) = B(86)+B(90)

  JVS(80) = B(79)+B(87)+B(96)+B(107)

  JVS(81) = B(83)+B(91)+B(100)+B(113)

  JVS(82) = 0

  JVS(83) = B(97)+B(101)+B(103)+B(105)+B(116)

  JVS(84) = B(108)+B(110)+B(114)+B(117)+B(118)

  JVS(85) = B(80)+B(88)+B(98)+B(109)

  JVS(86) = B(81)+B(84)+B(85)+B(92)+B(102)+B(111)

  JVS(87) = B(89)+B(93)+B(94)+B(104)+B(115)

  JVS(88) = 0

  JVS(89) = B(404)

  JVS(90) = B(412)

  JVS(91) = B(402)

  JVS(92) = B(410)

  JVS(93) = B(406)

  JVS(94) = B(408)

  JVS(95) = B(403)+B(405)+B(407)+B(409)+B(411)+B(413)

  JVS(96) = 0

  JVS(97) = B(414)+B(416)+0.0327*B(418)

  JVS(98) = B(426)+B(428)+0.0816*B(430)

  JVS(99) = B(420)+B(422)+0.0545*B(424)

  JVS(100) = B(417)+B(423)+B(429)

  JVS(101) = B(415)+B(421)+B(427)

  JVS(102) = 0.0327*B(419)+0.0545*B(425)+0.0816*B(431)

  JVS(103) = 0

  JVS(104) = B(319)

  JVS(105) = B(353)

  JVS(106) = B(355)

  JVS(107) = B(245)

  JVS(108) = B(74)

  JVS(109) = B(40)

  JVS(110) = B(363)

  JVS(111) = B(359)

  JVS(112) = B(247)

  JVS(113) = B(367)

  JVS(114) = B(243)

  JVS(115) = B(361)

  JVS(116) = B(365)

  JVS(117) = B(313)

  JVS(118) = B(316)

  JVS(119) = B(269)

  JVS(120) = B(309)

  JVS(121) = B(357)

  JVS(122) = B(250)

  JVS(123) = B(275)

  JVS(124) = B(265)

  JVS(125) = B(260)

  JVS(126) = B(49)

  JVS(127) = B(46)

  JVS(128) = B(321)

  JVS(129) = B(237)

  JVS(130) = B(387)

  JVS(131) = B(329)

  JVS(132) = B(369)

  JVS(133) = B(280)

  JVS(134) = B(345)

  JVS(135) = B(337)

  JVS(136) = B(255)

  JVS(137) = B(296)

  JVS(138) = B(227)

  JVS(139) = B(377)

  JVS(140) = B(289)

  JVS(141) = B(218)

  JVS(142) = B(306)

  JVS(143) = B(232)

  JVS(144) = B(240)

  JVS(145) = B(303)

  JVS(146) = B(51)

  JVS(147) = B(72)

  JVS(148) = B(42)

  JVS(149) = B(36)

  JVS(150) = B(37)+B(41)+B(43)+B(44)+B(47)+B(50)+B(52)+B(73)+B(75)+B(76)+B(219)+B(228)+B(233)+B(238)+B(241)+B(244)&
               &+B(246)+B(248)+B(251)+B(256)+B(261)+B(266)+B(270)+B(276)+B(281)+B(290)+B(297)+B(304)+B(307)+B(310)+B(314)&
               &+B(317)+B(320)+B(322)+B(330)+B(338)+B(346)+B(354)+B(356)+B(358)+B(360)+B(362)+B(364)+B(366)+B(368)+B(370)&
               &+B(378)+B(388)

  JVS(151) = B(45)

  JVS(152) = B(432)

  JVS(153) = B(434)

  JVS(154) = B(440)

  JVS(155) = B(442)

  JVS(156) = 0

  JVS(157) = B(438)

  JVS(158) = B(450)

  JVS(159) = B(446)

  JVS(160) = B(454)

  JVS(161) = B(448)

  JVS(162) = B(436)

  JVS(163) = B(452)

  JVS(164) = B(444)

  JVS(165) = B(433)+B(435)+B(437)+B(439)+B(441)+B(443)+B(445)+B(447)+B(449)+B(451)+B(453)+B(455)

  JVS(166) = 0

  JVS(167) = B(450)

  JVS(168) = 0.5*B(448)

  JVS(169) = 0.5*B(449)+B(451)

  JVS(170) = -B(450)

  JVS(171) = -B(451)

  JVS(172) = 0

  JVS(173) = B(454)

  JVS(174) = 0.5*B(452)

  JVS(175) = 0.5*B(453)+B(455)

  JVS(176) = -B(454)

  JVS(177) = -B(455)

  JVS(178) = -B(32)-B(34)

  JVS(179) = B(31)

  JVS(180) = -B(448)

  JVS(181) = -B(449)

  JVS(182) = B(448)

  JVS(183) = 0

  JVS(184) = B(449)

  JVS(185) = -B(452)

  JVS(186) = -B(453)

  JVS(187) = B(452)

  JVS(188) = 0

  JVS(189) = B(453)

  JVS(190) = -B(319)

  JVS(191) = -B(320)

  JVS(192) = -B(353)

  JVS(193) = -B(354)

  JVS(194) = -B(121)

  JVS(195) = B(119)

  JVS(196) = B(120)

  JVS(197) = -B(139)

  JVS(198) = B(137)

  JVS(199) = B(138)

  JVS(200) = -B(159)

  JVS(201) = B(157)

  JVS(202) = B(158)

  JVS(203) = -B(181)

  JVS(204) = B(179)

  JVS(205) = B(180)

  JVS(206) = -B(69)-B(70)-B(400)

  JVS(207) = B(63)+B(64)

  JVS(208) = -B(71)

  JVS(209) = -B(355)

  JVS(210) = -B(356)

  JVS(211) = -B(264)

  JVS(212) = 0.087*B(367)

  JVS(213) = 0.031*B(347)

  JVS(214) = 0.031*B(339)

  JVS(215) = 0.031*B(340)+0.031*B(348)

  JVS(216) = 0.087*B(368)

  JVS(217) = -B(245)

  JVS(218) = -B(246)

  JVS(219) = -B(23)-B(24)

  JVS(220) = B(21)

  JVS(221) = B(22)

  JVS(222) = -B(74)-B(395)-B(397)

  JVS(223) = B(456)+0.5*B(458)+B(460)

  JVS(224) = -B(75)+B(457)+0.5*B(459)

  JVS(225) = B(461)

  JVS(226) = -B(456)-B(458)-B(460)

  JVS(227) = -B(457)-B(459)

  JVS(228) = -B(461)

  JVS(229) = -B(38)-B(39)-B(40)

  JVS(230) = B(36)

  JVS(231) = B(37)-B(41)

  JVS(232) = -B(363)

  JVS(233) = -B(364)

  JVS(234) = -B(359)

  JVS(235) = -B(360)

  JVS(236) = 0.236*B(359)

  JVS(237) = -B(203)-B(205)

  JVS(238) = -B(204)

  JVS(239) = 0.236*B(360)

  JVS(240) = -B(247)-B(249)

  JVS(241) = B(80)

  JVS(242) = B(81)

  JVS(243) = -B(248)

  JVS(244) = -B(222)-B(223)

  JVS(245) = B(220)

  JVS(246) = B(221)

  JVS(247) = -B(224)

  JVS(248) = -B(211)-B(213)-B(215)

  JVS(249) = B(273)

  JVS(250) = -B(214)

  JVS(251) = -B(212)

  JVS(252) = B(274)

  JVS(253) = -B(367)

  JVS(254) = -B(368)

  JVS(255) = -B(57)-B(58)-B(59)

  JVS(256) = B(55)

  JVS(257) = B(56)

  JVS(258) = -B(60)

  JVS(259) = -B(243)

  JVS(260) = 0.25*B(110)

  JVS(261) = B(84)+0.25*B(92)+0.25*B(111)

  JVS(262) = 0.25*B(93)

  JVS(263) = -B(244)

  JVS(264) = -B(361)

  JVS(265) = -B(362)

  JVS(266) = -B(365)

  JVS(267) = -B(366)

  JVS(268) = 0.099*B(367)

  JVS(269) = 0.108*B(365)

  JVS(270) = -B(313)-B(315)

  JVS(271) = -B(314)+0.108*B(366)+0.099*B(368)

  JVS(272) = 0.093*B(367)

  JVS(273) = 0.051*B(365)

  JVS(274) = -B(316)-B(318)

  JVS(275) = -B(317)+0.051*B(366)+0.093*B(368)

  JVS(276) = 0.187*B(367)

  JVS(277) = 0.207*B(365)

  JVS(278) = -B(269)-B(271)

  JVS(279) = -B(270)+0.207*B(366)+0.187*B(368)

  JVS(280) = -B(272)

  JVS(281) = 0.561*B(367)

  JVS(282) = 0.491*B(365)

  JVS(283) = -B(309)-B(311)

  JVS(284) = -B(312)

  JVS(285) = -B(310)+0.491*B(366)+0.561*B(368)

  JVS(286) = -B(357)-B(385)

  JVS(287) = -B(386)

  JVS(288) = -B(358)

  JVS(289) = -B(250)-B(252)

  JVS(290) = B(108)

  JVS(291) = B(88)+B(109)

  JVS(292) = B(89)

  JVS(293) = -B(251)

  JVS(294) = B(213)+B(215)

  JVS(295) = -B(273)

  JVS(296) = B(206)

  JVS(297) = B(214)

  JVS(298) = B(207)

  JVS(299) = -B(274)

  JVS(300) = 0.05*B(367)

  JVS(301) = 0.059*B(365)

  JVS(302) = -B(275)-B(277)-B(278)

  JVS(303) = 0.061*B(377)+0.042*B(379)+0.015*B(381)

  JVS(304) = 0.042*B(380)

  JVS(305) = -B(276)+0.059*B(366)+0.05*B(368)+0.061*B(378)

  JVS(306) = -B(279)+0.015*B(382)

  JVS(307) = 0.017*B(365)

  JVS(308) = -B(265)-B(267)

  JVS(309) = B(208)+B(210)

  JVS(310) = B(209)

  JVS(311) = -B(266)+0.017*B(366)

  JVS(312) = -B(268)

  JVS(313) = 0.287*B(367)

  JVS(314) = 0.119*B(365)

  JVS(315) = 0.5*B(315)

  JVS(316) = 0.5*B(318)

  JVS(317) = 0.23*B(269)

  JVS(318) = -B(259)-B(260)-B(262)

  JVS(319) = 0.084*B(280)+0.9*B(282)

  JVS(320) = 0.174*B(296)+0.742*B(298)+0.008*B(300)

  JVS(321) = 0.3*B(289)+0.95*B(291)

  JVS(322) = 0.9*B(283)+0.95*B(292)+0.742*B(299)

  JVS(323) = -B(261)+0.23*B(270)+0.084*B(281)+0.3*B(290)+0.174*B(297)+0.119*B(366)+0.287*B(368)

  JVS(324) = -B(263)+0.008*B(301)

  JVS(325) = 0.002*B(361)

  JVS(326) = B(315)

  JVS(327) = B(318)

  JVS(328) = B(309)+1.5*B(311)

  JVS(329) = 0.393*B(357)+1.5*B(385)

  JVS(330) = B(259)+B(260)+B(262)

  JVS(331) = -B(49)

  JVS(332) = 0.5*B(323)+0.491*B(327)

  JVS(333) = 0.51*B(389)

  JVS(334) = 0.275*B(331)

  JVS(335) = 0.345*B(371)

  JVS(336) = 0.416*B(280)+0.45*B(282)+0.5*B(284)+0.67*B(288)

  JVS(337) = 0.157*B(347)

  JVS(338) = 0.157*B(339)

  JVS(339) = 2*B(253)+B(254)+1.26*B(255)+1.26*B(257)

  JVS(340) = 0.336*B(296)+0.498*B(298)+0.572*B(300)+1.233*B(302)

  JVS(341) = B(229)

  JVS(342) = 0.265*B(379)+0.012*B(383)

  JVS(343) = 0.475*B(291)+0.7*B(295)

  JVS(344) = B(216)+B(217)+B(218)+B(225)

  JVS(345) = 0.491*B(328)+0.012*B(384)

  JVS(346) = 0.034*B(232)+B(234)

  JVS(347) = 0.45*B(283)+0.475*B(292)+0.498*B(299)+1.5*B(312)+0.5*B(324)+0.275*B(332)+0.157*B(340)+0.157*B(348)+0.345&
               &*B(372)+0.265*B(380)+1.5*B(386)+0.51*B(390)

  JVS(348) = -B(50)+B(219)+0.034*B(233)+1.26*B(256)+B(261)+0.416*B(281)+0.336*B(297)+B(310)+0.393*B(358)+0.002*B(362)

  JVS(349) = B(226)+1.26*B(258)+B(263)+0.5*B(285)+0.572*B(301)

  JVS(350) = 2*B(24)

  JVS(351) = B(460)

  JVS(352) = B(271)

  JVS(353) = B(273)

  JVS(354) = B(278)

  JVS(355) = B(267)

  JVS(356) = B(262)

  JVS(357) = -B(46)-B(48)-B(399)

  JVS(358) = 0

  JVS(359) = 0.5*B(284)

  JVS(360) = B(257)

  JVS(361) = 0.15*B(300)

  JVS(362) = B(230)

  JVS(363) = 0

  JVS(364) = 0

  JVS(365) = B(225)

  JVS(366) = B(235)

  JVS(367) = 0

  JVS(368) = 0.2*B(66)

  JVS(369) = B(42)

  JVS(370) = B(43)-B(47)

  JVS(371) = 0.2*B(67)+B(226)+B(231)+B(236)+B(258)+B(263)+B(268)+B(272)+B(274)+B(279)+0.5*B(285)+0.15*B(301)+B(461)

  JVS(372) = -B(321)-B(323)-B(325)-B(327)

  JVS(373) = -B(328)

  JVS(374) = -B(324)

  JVS(375) = -B(322)

  JVS(376) = -B(326)

  JVS(377) = 0.704*B(355)

  JVS(378) = 0.072*B(363)

  JVS(379) = 0.024*B(359)

  JVS(380) = B(205)

  JVS(381) = 0.452*B(361)

  JVS(382) = -B(237)-B(239)

  JVS(383) = 0.005*B(369)+0.001*B(371)+0.024*B(373)

  JVS(384) = 0.13*B(347)

  JVS(385) = 0.13*B(339)

  JVS(386) = 0.127*B(377)+0.045*B(379)+0.102*B(381)

  JVS(387) = 0.006*B(306)+0.02*B(308)

  JVS(388) = 0.13*B(340)+0.13*B(348)+0.001*B(372)+0.045*B(380)

  JVS(389) = 0

  JVS(390) = -B(238)+0.006*B(307)+0.704*B(356)+0.024*B(360)+0.452*B(362)+0.072*B(364)+0.005*B(370)+0.127*B(378)

  JVS(391) = 0.024*B(374)+0.102*B(382)

  JVS(392) = -B(387)-B(389)-B(391)-B(393)

  JVS(393) = -B(394)

  JVS(394) = -B(390)

  JVS(395) = -B(388)

  JVS(396) = -B(392)

  JVS(397) = 0.24*B(269)+B(271)

  JVS(398) = 0.24*B(265)+B(267)

  JVS(399) = -B(206)-B(208)-B(210)

  JVS(400) = B(176)

  JVS(401) = B(200)

  JVS(402) = -B(209)

  JVS(403) = -B(207)

  JVS(404) = B(160)+B(164)+B(174)+B(177)+2*B(178)+B(201)

  JVS(405) = B(175)

  JVS(406) = B(161)

  JVS(407) = 0.24*B(266)+0.24*B(270)

  JVS(408) = B(165)+B(268)+B(272)

  JVS(409) = -B(329)-B(331)-B(333)-B(335)

  JVS(410) = -B(336)

  JVS(411) = -B(332)

  JVS(412) = -B(330)

  JVS(413) = -B(334)

  JVS(414) = -B(369)-B(371)-B(373)-B(375)

  JVS(415) = -B(376)

  JVS(416) = -B(372)

  JVS(417) = -B(370)

  JVS(418) = -B(374)

  JVS(419) = 0.948*B(363)

  JVS(420) = 0.559*B(359)

  JVS(421) = 0.936*B(361)

  JVS(422) = B(313)+B(315)

  JVS(423) = B(316)+B(318)

  JVS(424) = B(237)

  JVS(425) = 0.079*B(329)+0.126*B(331)+0.187*B(333)+0.24*B(335)

  JVS(426) = 0.205*B(369)+0.488*B(373)

  JVS(427) = -B(95)-B(97)-B(99)-B(101)-B(103)-B(116)-B(132)-B(150)-B(170)-B(192)

  JVS(428) = 0.5*B(345)+0.729*B(347)+0.75*B(349)

  JVS(429) = 0.5*B(337)+0.729*B(339)+0.75*B(341)

  JVS(430) = 0.001*B(377)+0.137*B(379)+0.711*B(381)

  JVS(431) = 0.675*B(289)

  JVS(432) = 0.596*B(306)+0.152*B(308)

  JVS(433) = 0.24*B(336)

  JVS(434) = 0.616*B(240)

  JVS(435) = 0.515*B(305)

  JVS(436) = 0.126*B(332)+0.729*B(340)+0.729*B(348)+0.137*B(380)

  JVS(437) = -B(117)

  JVS(438) = -B(151)+B(176)

  JVS(439) = -B(193)+B(200)

  JVS(440) = -B(98)

  JVS(441) = 0

  JVS(442) = -B(102)

  JVS(443) = -B(104)

  JVS(444) = B(160)+B(164)-B(171)+B(174)+B(177)+2*B(178)+B(201)

  JVS(445) = -B(133)+B(175)

  JVS(446) = -B(96)+B(161)

  JVS(447) = B(238)+0.616*B(241)+0.675*B(290)+0.596*B(307)+B(314)+B(317)+0.079*B(330)+0.5*B(338)+0.5*B(346)+0.559*B(360)&
               &+0.936*B(362)+0.948*B(364)+0.205*B(370)+0.001*B(378)

  JVS(448) = -B(100)+B(165)+0.187*B(334)+0.75*B(342)+0.75*B(350)+0.488*B(374)+0.711*B(382)

  JVS(449) = 0.23*B(329)+0.39*B(331)

  JVS(450) = -B(280)-B(282)-B(284)-B(286)-B(288)

  JVS(451) = 0.025*B(377)+0.026*B(379)+0.012*B(383)

  JVS(452) = -B(287)+0.012*B(384)

  JVS(453) = -B(283)+0.39*B(332)+0.026*B(380)

  JVS(454) = -B(281)+0.23*B(330)+0.025*B(378)

  JVS(455) = -B(285)

  JVS(456) = -B(345)-B(347)-B(349)-B(351)

  JVS(457) = -B(352)

  JVS(458) = -B(348)

  JVS(459) = -B(346)

  JVS(460) = -B(350)

  JVS(461) = -B(337)-B(339)-B(341)-B(343)

  JVS(462) = -B(344)

  JVS(463) = -B(340)

  JVS(464) = -B(338)

  JVS(465) = -B(342)

  JVS(466) = 0.097*B(367)

  JVS(467) = 0.118*B(365)

  JVS(468) = 0.5*B(315)

  JVS(469) = 0.5*B(318)

  JVS(470) = B(311)

  JVS(471) = 0.607*B(357)

  JVS(472) = 0.23*B(265)

  JVS(473) = 0.009*B(327)

  JVS(474) = 0

  JVS(475) = 0.001*B(347)

  JVS(476) = 0.001*B(339)

  JVS(477) = -B(253)-B(254)-B(255)-B(257)

  JVS(478) = 0.15*B(296)+0.023*B(298)

  JVS(479) = 0.009*B(328)

  JVS(480) = 0.023*B(299)+B(312)+0.001*B(340)+0.001*B(348)

  JVS(481) = 0

  JVS(482) = 0

  JVS(483) = 0

  JVS(484) = 0

  JVS(485) = 0

  JVS(486) = 0

  JVS(487) = 0

  JVS(488) = -B(256)+0.23*B(266)+0.15*B(297)+0.607*B(358)+0.118*B(366)+0.097*B(368)

  JVS(489) = -B(258)

  JVS(490) = 0.357*B(329)+0.936*B(333)

  JVS(491) = -B(296)-B(298)-B(300)-B(302)

  JVS(492) = 0.025*B(377)

  JVS(493) = 0

  JVS(494) = -B(299)

  JVS(495) = -B(297)+0.357*B(330)+0.025*B(378)

  JVS(496) = -B(301)+0.936*B(334)

  JVS(497) = B(353)

  JVS(498) = 0.96*B(245)

  JVS(499) = 0.099*B(363)

  JVS(500) = 0.445*B(359)

  JVS(501) = 0.455*B(361)

  JVS(502) = 0.195*B(321)+0.25*B(327)

  JVS(503) = 0.984*B(387)+0.5*B(389)

  JVS(504) = 0.294*B(369)+0.154*B(371)+0.009*B(373)

  JVS(505) = 0.129*B(296)+0.047*B(298)+0.467*B(302)

  JVS(506) = -B(227)-B(229)-B(230)

  JVS(507) = 0.732*B(377)+0.456*B(379)+0.507*B(381)

  JVS(508) = 0.439*B(306)+0.431*B(308)

  JVS(509) = 0.25*B(328)

  JVS(510) = 0.034*B(232)+B(234)

  JVS(511) = 0.482*B(240)+B(242)

  JVS(512) = 0.084*B(303)+0.246*B(305)

  JVS(513) = 0.047*B(299)+0.154*B(372)+0.456*B(380)+0.5*B(390)

  JVS(514) = B(140)+B(144)+B(154)+2*B(156)+B(176)+B(198)

  JVS(515) = B(199)

  JVS(516) = B(177)

  JVS(517) = B(155)

  JVS(518) = B(141)

  JVS(519) = -B(228)+0.034*B(233)+0.482*B(241)+0.96*B(246)+0.129*B(297)+0.084*B(304)+0.439*B(307)+0.195*B(322)+B(354)&
               &+0.445*B(360)+0.455*B(362)+0.099*B(364)+0.294*B(370)+0.732*B(378)+0.984*B(388)

  JVS(520) = B(145)-B(231)+0.009*B(374)+0.507*B(382)

  JVS(521) = -B(377)-B(379)-B(381)-B(383)

  JVS(522) = -B(384)

  JVS(523) = -B(380)

  JVS(524) = -B(378)

  JVS(525) = -B(382)

  JVS(526) = 0.32*B(329)+0.16*B(331)

  JVS(527) = 0.019*B(379)+0.048*B(381)

  JVS(528) = -B(289)-B(291)-B(293)-B(295)

  JVS(529) = -B(294)

  JVS(530) = -B(292)+0.16*B(332)+0.019*B(380)

  JVS(531) = -B(290)+0.32*B(330)

  JVS(532) = 0.048*B(382)

  JVS(533) = 0.081*B(245)

  JVS(534) = 0.026*B(363)

  JVS(535) = 0.026*B(359)

  JVS(536) = 0.35*B(247)+B(249)

  JVS(537) = B(222)

  JVS(538) = B(243)

  JVS(539) = 0.024*B(361)

  JVS(540) = 0.096*B(357)

  JVS(541) = 1.61*B(321)+B(323)+0.191*B(327)

  JVS(542) = B(237)

  JVS(543) = 0.984*B(387)+0.5*B(389)

  JVS(544) = 0.624*B(329)+0.592*B(331)+0.24*B(335)

  JVS(545) = 0.732*B(369)+0.5*B(371)

  JVS(546) = 0.084*B(280)+0.2*B(282)+0.67*B(288)

  JVS(547) = 0.276*B(345)+0.235*B(347)

  JVS(548) = 0.276*B(337)+0.235*B(339)

  JVS(549) = B(254)

  JVS(550) = 0.055*B(296)+0.125*B(298)+0.227*B(300)+0.3*B(302)

  JVS(551) = 0.244*B(377)+0.269*B(379)+0.079*B(381)

  JVS(552) = 0.3*B(289)+0.1*B(291)

  JVS(553) = -B(216)-B(217)-B(218)-B(220)-B(225)

  JVS(554) = 0.01*B(306)+0.134*B(308)

  JVS(555) = 0.191*B(328)+0.24*B(336)

  JVS(556) = 0.115*B(240)

  JVS(557) = 0.213*B(303)+0.506*B(305)

  JVS(558) = 0.2*B(283)+0.1*B(292)+0.125*B(299)+B(324)+0.592*B(332)+0.235*B(340)+0.235*B(348)+0.5*B(372)+0.269*B(380)&
               &+0.5*B(390)

  JVS(559) = 0.75*B(110)

  JVS(560) = B(146)+B(198)

  JVS(561) = B(182)+B(186)+B(188)+B(196)+B(199)+B(200)+2*B(202)

  JVS(562) = -B(221)

  JVS(563) = 0

  JVS(564) = B(78)+B(82)+B(84)+2*B(85)+0.75*B(92)+0.75*B(111)+B(128)+B(147)+B(166)+B(189)

  JVS(565) = 0.75*B(93)

  JVS(566) = B(167)+B(201)

  JVS(567) = B(129)+B(197)

  JVS(568) = B(79)+B(183)

  JVS(569) = -B(219)+B(238)+0.115*B(241)+B(244)+0.081*B(246)+0.35*B(248)+0.084*B(281)+0.3*B(290)+0.055*B(297)+0.213&
               &*B(304)+0.01*B(307)+1.61*B(322)+0.624*B(330)+0.276*B(338)+0.276*B(346)+0.096*B(358)+0.026*B(360)+0.024&
               &*B(362)+0.026*B(364)+0.732*B(370)+0.244*B(378)+0.984*B(388)

  JVS(570) = B(83)+B(187)-B(226)+0.227*B(301)+0.079*B(382)

  JVS(571) = B(203)

  JVS(572) = 0.511*B(373)

  JVS(573) = 0.276*B(349)

  JVS(574) = 0.276*B(341)

  JVS(575) = 0.572*B(300)

  JVS(576) = 0.321*B(381)

  JVS(577) = -0.69*B(306)-B(308)

  JVS(578) = 0

  JVS(579) = 0

  JVS(580) = B(106)

  JVS(581) = B(204)

  JVS(582) = B(107)

  JVS(583) = -0.69*B(307)

  JVS(584) = 0.572*B(301)+0.276*B(342)+0.276*B(350)+0.511*B(374)+0.321*B(382)

  JVS(585) = B(34)

  JVS(586) = -B(327)

  JVS(587) = -B(393)

  JVS(588) = -B(335)

  JVS(589) = -B(375)

  JVS(590) = -B(286)

  JVS(591) = -B(351)

  JVS(592) = -B(343)

  JVS(593) = -B(383)

  JVS(594) = -B(293)

  JVS(595) = -B(2)-B(4)-B(6)-B(9)-B(11)-B(287)-B(294)-B(328)-B(336)-B(344)-B(352)-B(376)-B(384)-B(394)

  JVS(596) = -B(5)+B(30)

  JVS(597) = B(1)-B(10)-B(12)

  JVS(598) = -B(7)

  JVS(599) = 0

  JVS(600) = B(29)

  JVS(601) = 0.261*B(355)

  JVS(602) = 0.204*B(363)

  JVS(603) = 0.122*B(359)

  JVS(604) = 0.244*B(361)

  JVS(605) = B(313)

  JVS(606) = B(316)

  JVS(607) = B(309)

  JVS(608) = B(250)+B(252)

  JVS(609) = B(325)

  JVS(610) = 0.45*B(393)

  JVS(611) = 0.497*B(369)+0.363*B(371)+0.037*B(373)+0.45*B(375)

  JVS(612) = B(286)

  JVS(613) = 0.474*B(345)+0.205*B(347)+0.474*B(349)+0.147*B(351)

  JVS(614) = 0.474*B(337)+0.205*B(339)+0.474*B(341)+0.147*B(343)

  JVS(615) = 0.013*B(296)+0.218*B(300)

  JVS(616) = 0.511*B(377)+0.305*B(379)+0.151*B(381)+0.069*B(383)

  JVS(617) = 0.675*B(289)+0.45*B(293)

  JVS(618) = 0.213*B(306)+0.147*B(308)

  JVS(619) = B(287)+0.45*B(294)+0.147*B(344)+0.147*B(352)+0.45*B(376)+0.069*B(384)+0.45*B(394)

  JVS(620) = -B(232)-B(234)-B(235)

  JVS(621) = 0.37*B(240)

  JVS(622) = 0.558*B(303)+0.71*B(305)

  JVS(623) = 0.205*B(340)+0.205*B(348)+0.363*B(372)+0.305*B(380)

  JVS(624) = 0

  JVS(625) = 0

  JVS(626) = 0

  JVS(627) = 0

  JVS(628) = 0

  JVS(629) = -B(233)+0.37*B(241)+B(251)+0.675*B(290)+0.013*B(297)+0.558*B(304)+0.213*B(307)+B(310)+B(314)+B(317)+0.474&
               &*B(338)+0.474*B(346)+0.261*B(356)+0.122*B(360)+0.244*B(362)+0.204*B(364)+0.497*B(370)+0.511*B(378)

  JVS(630) = -B(236)+0.218*B(301)+B(326)+0.474*B(342)+0.474*B(350)+0.037*B(374)+0.151*B(382)

  JVS(631) = 0.089*B(363)

  JVS(632) = 0.332*B(359)

  JVS(633) = 0.11*B(361)

  JVS(634) = 0.55*B(393)

  JVS(635) = 0.437*B(375)

  JVS(636) = 0.416*B(280)

  JVS(637) = 0.15*B(296)+0.21*B(298)+0.233*B(302)

  JVS(638) = 0.072*B(377)+0.026*B(379)+0.001*B(381)+0.659*B(383)

  JVS(639) = 0.55*B(293)

  JVS(640) = 0.177*B(306)+0.243*B(308)

  JVS(641) = 0.55*B(294)+0.437*B(376)+0.659*B(384)+0.55*B(394)

  JVS(642) = -B(240)-B(242)

  JVS(643) = 0.115*B(303)

  JVS(644) = 0.21*B(299)+0.026*B(380)

  JVS(645) = 0.5*B(110)+B(112)+0.5*B(114)+B(118)

  JVS(646) = 0

  JVS(647) = 0.5*B(111)

  JVS(648) = 0.5*B(115)

  JVS(649) = 0

  JVS(650) = -B(241)+0.416*B(281)+0.15*B(297)+0.115*B(304)+0.177*B(307)+0.332*B(360)+0.11*B(362)+0.089*B(364)+0.072&
               &*B(378)

  JVS(651) = B(113)+0.001*B(382)

  JVS(652) = 0.417*B(363)

  JVS(653) = 0.125*B(361)

  JVS(654) = 0.055*B(365)

  JVS(655) = 0.1*B(331)+0.75*B(335)

  JVS(656) = 0.119*B(369)+0.215*B(371)+0.113*B(375)

  JVS(657) = 0.276*B(345)+0.276*B(347)+0.853*B(351)

  JVS(658) = 0.276*B(337)+0.276*B(339)+0.853*B(343)

  JVS(659) = 0.332*B(296)

  JVS(660) = 0.043*B(379)+0.259*B(383)

  JVS(661) = 0.7*B(295)

  JVS(662) = 0.048*B(306)+0.435*B(308)

  JVS(663) = 0.75*B(336)+0.853*B(344)+0.853*B(352)+0.113*B(376)+0.259*B(384)

  JVS(664) = -0.671*B(303)-B(305)

  JVS(665) = 0.1*B(332)+0.276*B(340)+0.276*B(348)+0.215*B(372)+0.043*B(380)

  JVS(666) = 0.5*B(110)+0.5*B(114)+B(118)+B(134)+B(152)+B(172)

  JVS(667) = B(153)

  JVS(668) = 0

  JVS(669) = 0.5*B(111)

  JVS(670) = 0.5*B(115)

  JVS(671) = B(173)

  JVS(672) = B(135)

  JVS(673) = 0

  JVS(674) = 0.332*B(297)-0.671*B(304)+0.048*B(307)+0.276*B(338)+0.276*B(346)+0.125*B(362)+0.417*B(364)+0.055*B(366)&
               &+0.119*B(370)

  JVS(675) = 0

  JVS(676) = -B(311)

  JVS(677) = -B(385)

  JVS(678) = -B(323)

  JVS(679) = -B(389)

  JVS(680) = -B(331)

  JVS(681) = -B(371)

  JVS(682) = -B(282)

  JVS(683) = -B(347)

  JVS(684) = -B(339)

  JVS(685) = -B(298)

  JVS(686) = -B(379)

  JVS(687) = -B(291)

  JVS(688) = B(2)-B(4)

  JVS(689) = -B(5)-B(13)-B(15)-B(30)-B(31)-B(51)-B(61)-B(283)-B(292)-B(299)-B(312)-B(324)-B(332)-B(340)-B(348)-B(372)&
               &-B(380)-B(386)-B(390)

  JVS(690) = 0.25*B(142)

  JVS(691) = 0.25*B(184)

  JVS(692) = -B(62)+0.25*B(124)+0.25*B(143)+0.25*B(162)+0.25*B(185)

  JVS(693) = -B(16)

  JVS(694) = 0.25*B(163)

  JVS(695) = 0.25*B(125)

  JVS(696) = -B(14)

  JVS(697) = -B(52)

  JVS(698) = 0

  JVS(699) = 0.035*B(355)

  JVS(700) = 0.347*B(363)

  JVS(701) = 0.07*B(359)

  JVS(702) = 0.009*B(367)

  JVS(703) = 0.143*B(361)

  JVS(704) = 0.011*B(365)

  JVS(705) = 0.016*B(387)+0.051*B(391)

  JVS(706) = 0.093*B(329)+0.008*B(331)+0.064*B(333)+0.01*B(335)

  JVS(707) = 0.09*B(369)+0.001*B(371)+0.176*B(373)

  JVS(708) = 0.25*B(345)+0.18*B(347)+0.25*B(349)

  JVS(709) = 0.25*B(337)+0.18*B(339)+0.25*B(341)

  JVS(710) = 0.041*B(296)+0.051*B(300)

  JVS(711) = 0.082*B(377)+0.002*B(379)+0.136*B(381)+0.001*B(383)

  JVS(712) = 0.025*B(289)

  JVS(713) = 0.173*B(306)+0.095*B(308)

  JVS(714) = 0.01*B(336)+0.001*B(384)

  JVS(715) = 0.001*B(232)

  JVS(716) = 0.042*B(240)

  JVS(717) = 0.07*B(303)+0.04*B(305)

  JVS(718) = 0.008*B(332)+0.18*B(340)+0.18*B(348)+0.001*B(372)+0.002*B(380)

  JVS(719) = -B(106)-B(108)-B(110)-B(112)-B(114)-2*B(118)-B(134)-B(152)-B(172)-B(194)

  JVS(720) = -B(153)

  JVS(721) = -B(195)

  JVS(722) = -B(109)

  JVS(723) = 0

  JVS(724) = -B(111)

  JVS(725) = -B(115)

  JVS(726) = -B(173)

  JVS(727) = -B(135)

  JVS(728) = -B(107)

  JVS(729) = 0.001*B(233)+0.042*B(241)+0.025*B(290)+0.041*B(297)+0.07*B(304)+0.173*B(307)+0.093*B(330)+0.25*B(338)+0.25&
               &*B(346)+0.035*B(356)+0.07*B(360)+0.143*B(362)+0.347*B(364)+0.011*B(366)+0.009*B(368)+0.09*B(370)+0.082&
               &*B(378)+0.016*B(388)

  JVS(730) = -B(113)+0.051*B(301)+0.064*B(334)+0.25*B(342)+0.25*B(350)+0.176*B(374)+0.136*B(382)+0.051*B(392)

  JVS(731) = B(139)

  JVS(732) = 0.1*B(282)

  JVS(733) = 0.201*B(347)

  JVS(734) = 0.201*B(339)

  JVS(735) = 0.37*B(255)+0.37*B(257)

  JVS(736) = 0.048*B(298)+0.3*B(302)

  JVS(737) = 0.006*B(379)

  JVS(738) = 0.05*B(291)

  JVS(739) = 0

  JVS(740) = 0.965*B(232)+B(235)

  JVS(741) = 0.096*B(240)

  JVS(742) = 0.049*B(303)+0.333*B(305)

  JVS(743) = 0.1*B(283)+0.05*B(292)+0.048*B(299)+0.201*B(340)+0.201*B(348)+0.006*B(380)

  JVS(744) = -B(152)

  JVS(745) = -B(137)-B(140)-B(142)-B(144)-B(146)-B(148)-B(153)-B(154)-2*B(156)-B(176)-B(198)

  JVS(746) = -B(199)

  JVS(747) = -B(143)

  JVS(748) = -B(138)

  JVS(749) = -B(147)

  JVS(750) = -B(149)

  JVS(751) = -B(177)

  JVS(752) = -B(155)

  JVS(753) = -B(141)

  JVS(754) = 0.965*B(233)+0.096*B(241)+0.37*B(256)+0.049*B(304)

  JVS(755) = -B(145)+B(236)+0.37*B(258)

  JVS(756) = B(181)

  JVS(757) = 0.192*B(331)+0.24*B(335)

  JVS(758) = 0.5*B(280)+0.5*B(284)+0.33*B(288)

  JVS(759) = 0.289*B(296)+0.15*B(300)

  JVS(760) = 0

  JVS(761) = 0.3*B(295)

  JVS(762) = 0.24*B(336)

  JVS(763) = 0.192*B(332)

  JVS(764) = -B(194)

  JVS(765) = -B(198)

  JVS(766) = -B(179)-B(182)-B(184)-B(186)-B(188)-B(190)-B(195)-B(196)-B(199)-B(200)-2*B(202)

  JVS(767) = -B(185)

  JVS(768) = -B(180)

  JVS(769) = -B(189)

  JVS(770) = -B(191)

  JVS(771) = -B(201)

  JVS(772) = -B(197)

  JVS(773) = -B(183)

  JVS(774) = 0.5*B(281)+0.289*B(297)

  JVS(775) = -B(187)+0.5*B(285)+0.15*B(301)

  JVS(776) = B(70)

  JVS(777) = 0.95*B(245)

  JVS(778) = B(74)

  JVS(779) = 0.5*B(458)

  JVS(780) = B(39)

  JVS(781) = B(249)

  JVS(782) = B(222)+B(223)

  JVS(783) = -B(213)

  JVS(784) = 0.187*B(367)

  JVS(785) = B(57)+0.61*B(58)

  JVS(786) = B(243)

  JVS(787) = 0.224*B(365)

  JVS(788) = 0.5*B(315)

  JVS(789) = 0.5*B(318)

  JVS(790) = 1.5*B(311)

  JVS(791) = 0.297*B(357)+1.5*B(385)

  JVS(792) = B(252)

  JVS(793) = 0

  JVS(794) = B(259)

  JVS(795) = B(49)

  JVS(796) = 0.12*B(323)+0.5*B(327)

  JVS(797) = 0.06*B(389)

  JVS(798) = -B(208)

  JVS(799) = 0

  JVS(800) = 0.056*B(371)

  JVS(801) = 0.008*B(282)+0.34*B(288)

  JVS(802) = 0.033*B(347)

  JVS(803) = 0.033*B(339)

  JVS(804) = 2*B(253)+0.63*B(255)+0.63*B(257)

  JVS(805) = 0.4*B(298)+1.233*B(302)

  JVS(806) = B(229)

  JVS(807) = 0.003*B(379)+0.013*B(383)

  JVS(808) = 0.064*B(291)

  JVS(809) = 2*B(216)+B(218)-B(220)+B(225)

  JVS(810) = 0.113*B(306)+0.341*B(308)

  JVS(811) = 0.5*B(328)+0.013*B(384)

  JVS(812) = B(234)

  JVS(813) = 0

  JVS(814) = 0.379*B(303)

  JVS(815) = B(51)-B(61)+0.008*B(283)+0.064*B(292)+0.4*B(299)+1.5*B(312)+0.12*B(324)+0.033*B(340)+0.033*B(348)+0.056&
               &*B(372)+0.003*B(380)+1.5*B(386)+0.06*B(390)

  JVS(816) = -B(108)+B(110)+B(112)+B(114)+B(118)

  JVS(817) = -B(142)

  JVS(818) = -B(184)

  JVS(819) = -B(53)-B(55)-B(62)-2*B(63)-2*B(64)-B(66)-B(72)-B(80)-B(88)-B(109)-B(124)-B(143)-B(162)-B(185)-B(209)-B(214)&
               &-B(221)-B(396)

  JVS(820) = -B(56)

  JVS(821) = B(78)-B(81)+B(82)+2*B(85)+B(92)+B(111)

  JVS(822) = B(86)-B(89)+B(90)+B(93)+B(94)+B(115)

  JVS(823) = -B(163)

  JVS(824) = -B(125)

  JVS(825) = -B(54)+B(79)+B(87)+B(224)

  JVS(826) = B(44)+B(50)+B(52)+B(71)-B(73)+B(75)+B(76)+B(219)+B(244)+0.95*B(246)+0.63*B(256)+0.379*B(304)+0.113*B(307)&
               &+0.297*B(358)+0.224*B(366)+0.187*B(368)+0.5*B(459)

  JVS(827) = B(45)-B(67)+B(83)+B(91)+B(113)+B(226)+0.63*B(258)

  JVS(828) = B(121)

  JVS(829) = B(139)

  JVS(830) = B(159)

  JVS(831) = B(181)

  JVS(832) = B(23)

  JVS(833) = B(39)+B(40)

  JVS(834) = -B(203)

  JVS(835) = B(223)

  JVS(836) = -B(211)

  JVS(837) = B(57)+0.61*B(58)+B(59)

  JVS(838) = 0

  JVS(839) = B(48)

  JVS(840) = -B(206)

  JVS(841) = 0.187*B(333)

  JVS(842) = B(95)+B(99)

  JVS(843) = 0

  JVS(844) = 0.474*B(349)

  JVS(845) = 0.474*B(341)

  JVS(846) = 0

  JVS(847) = 0

  JVS(848) = 0

  JVS(849) = 0.391*B(381)

  JVS(850) = 0

  JVS(851) = 0

  JVS(852) = 0.338*B(306)+B(308)

  JVS(853) = B(6)-B(9)-B(11)

  JVS(854) = 0

  JVS(855) = 0

  JVS(856) = 0

  JVS(857) = B(13)-B(15)

  JVS(858) = B(112)

  JVS(859) = -B(137)+B(140)+B(144)

  JVS(860) = -B(179)+B(182)+B(186)

  JVS(861) = B(53)-B(55)+0.8*B(66)

  JVS(862) = -B(1)-B(10)-B(12)-B(16)-B(21)-B(42)-B(56)-B(119)-B(138)-B(157)-B(180)-B(204)-B(207)-B(212)

  JVS(863) = B(78)+B(82)

  JVS(864) = B(86)+B(90)

  JVS(865) = -B(158)+B(160)+B(164)

  JVS(866) = -B(120)+B(122)+B(126)

  JVS(867) = B(7)+B(14)+2*B(17)+2*B(19)+B(54)+B(79)+B(87)+B(96)+B(123)+B(141)+B(161)+B(183)+B(224)

  JVS(868) = B(41)-B(43)+B(44)+B(60)+0.338*B(307)

  JVS(869) = 2*B(18)-B(22)+B(29)+B(45)+0.8*B(67)+2*B(68)+B(83)+B(91)+B(100)+B(113)+B(127)+B(145)+B(165)+B(187)+0.187&
               &*B(334)+0.474*B(342)+0.474*B(350)+0.391*B(382)

  JVS(870) = B(319)

  JVS(871) = B(205)

  JVS(872) = 0.65*B(247)

  JVS(873) = 0.011*B(361)

  JVS(874) = 0.3*B(327)

  JVS(875) = B(239)

  JVS(876) = 0.26*B(389)

  JVS(877) = 0.25*B(335)

  JVS(878) = 0.076*B(371)

  JVS(879) = 0

  JVS(880) = 0

  JVS(881) = B(229)

  JVS(882) = 0.197*B(379)+0.03*B(381)

  JVS(883) = 0.3*B(295)

  JVS(884) = 0

  JVS(885) = 0.3*B(328)+0.25*B(336)

  JVS(886) = 0

  JVS(887) = 0

  JVS(888) = 0

  JVS(889) = 0.076*B(372)+0.197*B(380)+0.26*B(390)

  JVS(890) = -B(110)

  JVS(891) = -B(146)+B(154)

  JVS(892) = -B(188)+B(196)

  JVS(893) = -B(80)

  JVS(894) = 0

  JVS(895) = -B(78)-B(81)-B(82)-2*B(84)-2*B(85)-B(92)-B(111)-B(128)-B(147)-B(166)-B(189)

  JVS(896) = -B(93)

  JVS(897) = -B(167)+B(174)

  JVS(898) = B(122)+B(126)-B(129)+2*B(136)+B(155)+B(175)+B(197)

  JVS(899) = -B(79)+B(123)

  JVS(900) = 0.65*B(248)+B(320)+0.011*B(362)

  JVS(901) = -B(83)+B(127)+0.03*B(382)

  JVS(902) = B(353)

  JVS(903) = 0.965*B(355)

  JVS(904) = 0.05*B(245)

  JVS(905) = 0.653*B(363)

  JVS(906) = 0.695*B(359)

  JVS(907) = 0.804*B(367)

  JVS(908) = 0.835*B(361)

  JVS(909) = 0.765*B(365)

  JVS(910) = B(315)

  JVS(911) = B(318)

  JVS(912) = 0.76*B(269)

  JVS(913) = B(309)

  JVS(914) = 0.1*B(357)

  JVS(915) = 0.34*B(250)

  JVS(916) = 0.76*B(265)

  JVS(917) = B(321)+B(325)+0.2*B(327)

  JVS(918) = 0.984*B(387)+0.949*B(391)

  JVS(919) = 0

  JVS(920) = 0.907*B(329)+0.066*B(331)+0.749*B(333)

  JVS(921) = 0.91*B(369)+0.022*B(371)+0.824*B(373)

  JVS(922) = 0.5*B(280)+0.1*B(282)+0.5*B(284)+0.33*B(288)

  JVS(923) = 0.75*B(345)+0.031*B(347)+0.276*B(349)

  JVS(924) = 0.75*B(337)+0.031*B(339)+0.276*B(341)

  JVS(925) = 0.67*B(296)+0.048*B(298)+0.799*B(300)

  JVS(926) = 0.918*B(377)+0.033*B(379)+0.442*B(381)+0.012*B(383)

  JVS(927) = 0.3*B(289)+0.05*B(291)

  JVS(928) = 0.376*B(306)+0.564*B(308)

  JVS(929) = 0.2*B(328)+0.012*B(384)

  JVS(930) = 0.034*B(232)+B(234)

  JVS(931) = 0.37*B(240)+B(242)

  JVS(932) = 0.473*B(303)+0.96*B(305)

  JVS(933) = 0.1*B(283)+0.05*B(292)+0.048*B(299)+0.066*B(332)+0.031*B(340)+0.031*B(348)+0.022*B(372)+0.033*B(380)

  JVS(934) = -B(114)

  JVS(935) = B(140)+B(144)-B(148)+B(154)+2*B(156)+B(176)+B(198)

  JVS(936) = -B(190)+B(199)

  JVS(937) = -B(88)

  JVS(938) = 0

  JVS(939) = -B(92)

  JVS(940) = -B(86)-B(89)-B(90)-B(93)-2*B(94)-B(115)-B(130)-B(149)-B(168)-B(191)

  JVS(941) = -B(169)+B(177)

  JVS(942) = -B(131)+B(155)

  JVS(943) = -B(87)+B(141)

  JVS(944) = 0.034*B(233)+0.37*B(241)+0.05*B(246)+0.34*B(251)+0.76*B(266)+0.76*B(270)+0.5*B(281)+0.3*B(290)+0.67*B(297)&
               &+0.473*B(304)+0.376*B(307)+B(310)+B(322)+0.907*B(330)+0.75*B(338)+0.75*B(346)+B(354)+0.965*B(356)+0.1*B(358)&
               &+0.695*B(360)+0.835*B(362)+0.653*B(364)+0.765*B(366)+0.804*B(368)+0.91*B(370)+0.918*B(378)+0.984*B(388)

  JVS(945) = -B(91)+B(145)+0.5*B(285)+0.799*B(301)+B(326)+0.749*B(334)+0.276*B(342)+0.276*B(350)+0.824*B(374)+0.442&
               &*B(382)+0.949*B(392)

  JVS(946) = B(159)

  JVS(947) = B(275)+B(278)

  JVS(948) = 0

  JVS(949) = 0

  JVS(950) = 0

  JVS(951) = -B(172)

  JVS(952) = -B(176)

  JVS(953) = -B(200)

  JVS(954) = -B(162)

  JVS(955) = -B(157)

  JVS(956) = -B(166)

  JVS(957) = -B(168)

  JVS(958) = -B(158)-B(160)-B(163)-B(164)-B(167)-B(169)-B(173)-B(174)-B(177)-2*B(178)-B(201)

  JVS(959) = -B(175)

  JVS(960) = -B(161)

  JVS(961) = B(276)

  JVS(962) = -B(165)+B(279)

  JVS(963) = B(121)

  JVS(964) = 2*B(264)

  JVS(965) = 0

  JVS(966) = 0.011*B(361)

  JVS(967) = B(313)+0.5*B(315)

  JVS(968) = B(316)+0.5*B(318)

  JVS(969) = B(259)+B(260)+B(262)

  JVS(970) = B(237)+B(239)

  JVS(971) = 0

  JVS(972) = 0.67*B(288)

  JVS(973) = 0.123*B(347)

  JVS(974) = 0.123*B(339)

  JVS(975) = 0.467*B(302)

  JVS(976) = B(227)+B(230)

  JVS(977) = 0.137*B(379)

  JVS(978) = 0.675*B(289)

  JVS(979) = 0

  JVS(980) = 0

  JVS(981) = 0

  JVS(982) = 0.492*B(240)+B(242)

  JVS(983) = 0.029*B(303)+0.667*B(305)

  JVS(984) = 0.123*B(340)+0.123*B(348)+0.137*B(380)

  JVS(985) = -B(134)

  JVS(986) = -B(154)+B(198)

  JVS(987) = B(182)+B(186)+B(199)+B(200)+2*B(202)

  JVS(988) = -B(124)

  JVS(989) = -B(119)

  JVS(990) = -B(128)

  JVS(991) = -B(130)

  JVS(992) = -B(174)+B(201)

  JVS(993) = -B(120)-B(122)-B(125)-B(126)-B(129)-B(131)-B(135)-2*B(136)-B(155)-B(175)

  JVS(994) = -B(123)+B(183)

  JVS(995) = B(228)+B(238)+0.492*B(241)+B(261)+0.675*B(290)+0.029*B(304)+B(314)+B(317)+0.011*B(362)

  JVS(996) = -B(127)+B(187)+B(231)+B(263)

  JVS(997) = B(38)

  JVS(998) = -B(223)

  JVS(999) = -B(95)

  JVS(1000) = 0

  JVS(1001) = 0

  JVS(1002) = 0

  JVS(1003) = 0

  JVS(1004) = 0

  JVS(1005) = 0

  JVS(1006) = -B(6)+B(9)

  JVS(1007) = 0

  JVS(1008) = 0

  JVS(1009) = -B(13)

  JVS(1010) = -B(106)

  JVS(1011) = -B(140)

  JVS(1012) = -B(182)

  JVS(1013) = -B(53)

  JVS(1014) = B(1)+B(10)+B(26)

  JVS(1015) = -B(78)

  JVS(1016) = -B(86)

  JVS(1017) = -B(160)

  JVS(1018) = -B(122)

  JVS(1019) = -B(7)-B(14)-B(17)-2*B(19)-B(36)-B(54)-B(79)-B(87)-B(96)-B(107)-B(123)-B(141)-B(161)-B(183)-B(224)

  JVS(1020) = -B(37)

  JVS(1021) = -B(18)+B(27)+B(28)

  JVS(1022) = -B(432)

  JVS(1023) = -B(440)

  JVS(1024) = 2*B(32)

  JVS(1025) = -B(448)

  JVS(1026) = -B(436)

  JVS(1027) = -B(452)

  JVS(1028) = -B(444)

  JVS(1029) = -B(319)

  JVS(1030) = -B(353)

  JVS(1031) = 2*B(69)-B(70)

  JVS(1032) = -B(355)

  JVS(1033) = -B(245)

  JVS(1034) = -B(74)

  JVS(1035) = -B(456)-B(458)

  JVS(1036) = B(38)-B(40)

  JVS(1037) = -B(363)

  JVS(1038) = -B(359)

  JVS(1039) = -0.65*B(247)+B(249)

  JVS(1040) = -B(367)

  JVS(1041) = 0.39*B(58)-B(59)

  JVS(1042) = -B(243)

  JVS(1043) = -B(361)

  JVS(1044) = -B(365)

  JVS(1045) = -B(313)

  JVS(1046) = -B(316)

  JVS(1047) = -B(269)

  JVS(1048) = -B(309)+0.5*B(311)

  JVS(1049) = -0.397*B(357)+0.5*B(385)

  JVS(1050) = -0.34*B(250)+B(252)

  JVS(1051) = -B(275)

  JVS(1052) = -B(265)

  JVS(1053) = -B(260)

  JVS(1054) = -B(49)

  JVS(1055) = -B(46)+B(48)

  JVS(1056) = -B(321)+0.12*B(323)

  JVS(1057) = -B(237)

  JVS(1058) = -B(387)+0.32*B(389)

  JVS(1059) = 0

  JVS(1060) = -B(329)+0.266*B(331)

  JVS(1061) = -B(369)+0.155*B(371)

  JVS(1062) = -B(280)+0.208*B(282)+0.33*B(288)

  JVS(1063) = -B(345)+0.567*B(347)

  JVS(1064) = -B(337)+0.567*B(339)

  JVS(1065) = -B(255)

  JVS(1066) = -B(296)+0.285*B(298)

  JVS(1067) = -B(227)

  JVS(1068) = -B(377)+0.378*B(379)

  JVS(1069) = -B(289)+0.164*B(291)

  JVS(1070) = -B(218)

  JVS(1071) = -B(306)

  JVS(1072) = 0

  JVS(1073) = -B(232)

  JVS(1074) = -B(240)

  JVS(1075) = -B(303)

  JVS(1076) = -B(51)+B(61)+0.208*B(283)+0.164*B(292)+0.285*B(299)+0.5*B(312)+0.12*B(324)+0.266*B(332)+0.567*B(340)+0.567&
                &*B(348)+0.155*B(372)+0.378*B(380)+0.5*B(386)+0.32*B(390)

  JVS(1077) = 0

  JVS(1078) = 0

  JVS(1079) = 0

  JVS(1080) = B(53)+B(62)+0.8*B(66)-B(72)

  JVS(1081) = -B(42)

  JVS(1082) = 0

  JVS(1083) = 0

  JVS(1084) = 0

  JVS(1085) = 0

  JVS(1086) = -B(36)+B(54)

  JVS(1087) = -B(37)-B(41)-B(43)-B(44)-B(47)-B(50)-B(52)-B(60)-B(71)-B(73)-B(75)-B(76)-B(219)-B(228)-B(233)-B(238)&
                &-B(241)-B(244)-B(246)-0.65*B(248)-0.34*B(251)-B(256)-B(261)-B(266)-B(270)-B(276)-B(281)-B(290)-B(297)&
                &-B(304)-B(307)-B(310)-B(314)-B(317)-B(320)-B(322)-B(330)-B(338)-B(346)-B(354)-B(356)-0.397*B(358)-B(360)&
                &-B(362)-B(364)-B(366)-B(368)-B(370)-B(378)-B(388)-B(433)-B(437)-B(441)-B(445)-B(449)-B(453)-B(457)-B(459)

  JVS(1088) = -B(45)+0.8*B(67)

  JVS(1089) = B(23)

  JVS(1090) = -B(460)

  JVS(1091) = 0.39*B(58)

  JVS(1092) = -B(271)

  JVS(1093) = -B(273)

  JVS(1094) = -B(278)

  JVS(1095) = -B(267)

  JVS(1096) = -B(262)

  JVS(1097) = B(46)

  JVS(1098) = -B(325)

  JVS(1099) = -B(391)

  JVS(1100) = 0

  JVS(1101) = -B(333)

  JVS(1102) = -B(373)

  JVS(1103) = -B(99)

  JVS(1104) = -B(284)

  JVS(1105) = -B(349)

  JVS(1106) = -B(341)

  JVS(1107) = -B(257)

  JVS(1108) = -B(300)

  JVS(1109) = -B(230)

  JVS(1110) = -B(381)

  JVS(1111) = 0

  JVS(1112) = -B(225)

  JVS(1113) = 0

  JVS(1114) = B(11)

  JVS(1115) = -B(235)

  JVS(1116) = 0

  JVS(1117) = 0

  JVS(1118) = B(15)

  JVS(1119) = -B(112)

  JVS(1120) = -B(144)

  JVS(1121) = -B(186)

  JVS(1122) = -B(66)

  JVS(1123) = B(12)+B(16)-B(21)-B(26)

  JVS(1124) = -B(82)

  JVS(1125) = -B(90)

  JVS(1126) = -B(164)

  JVS(1127) = -B(126)

  JVS(1128) = -B(17)

  JVS(1129) = -B(44)+B(47)

  JVS(1130) = -B(18)-B(22)-B(27)-B(28)-B(29)-B(45)-B(67)-2*B(68)-B(83)-B(91)-B(100)-B(113)-B(127)-B(145)-B(165)-B(187)&
                &-B(226)-B(231)-B(236)-B(258)-B(263)-B(268)-B(272)-B(274)-B(279)-B(285)-B(301)-B(326)-B(334)-B(342)-B(350)&
                &-B(374)-B(382)-B(392)-B(461)
      
END SUBROUTINE saprc99_mosaic_8bin_vbs2_Jac_SP














SUBROUTINE saprc99_mosaic_8bin_vbs2_KppDecomp( JVS, IER )







      INTEGER  :: IER
      REAL(kind=dp) :: JVS(1130), W(99), a
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
      
END SUBROUTINE saprc99_mosaic_8bin_vbs2_KppDecomp



SUBROUTINE saprc99_mosaic_8bin_vbs2_KppDecompCmplx( JVS, IER )







      INTEGER  :: IER
      DOUBLE COMPLEX :: JVS(1130), W(99), a
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
      
END SUBROUTINE saprc99_mosaic_8bin_vbs2_KppDecompCmplx


SUBROUTINE saprc99_mosaic_8bin_vbs2_KppSolveIndirect( JVS, X )







      INTEGER i, j
      REAL(kind=dp) JVS(1130), X(99), sum

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
      
END SUBROUTINE saprc99_mosaic_8bin_vbs2_KppSolveIndirect


SUBROUTINE saprc99_mosaic_8bin_vbs2_KppSolveCmplx( JVS, X )







      INTEGER i, j
      DOUBLE COMPLEX JVS(1130), X(99), sum

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
      
END SUBROUTINE saprc99_mosaic_8bin_vbs2_KppSolveCmplx













SUBROUTINE saprc99_mosaic_8bin_vbs2_KppSolve ( JVS, X )


  REAL(kind=dp) :: JVS(LU_NONZERO)

  REAL(kind=dp) :: X(NVAR)

  X(19) = X(19)-JVS(152)*X(3)-JVS(153)*X(4)-JVS(154)*X(5)-JVS(155)*X(6)
  X(26) = X(26)-JVS(182)*X(25)
  X(28) = X(28)-JVS(187)*X(27)
  X(45) = X(45)-JVS(236)*X(44)
  X(54) = X(54)-JVS(268)*X(49)-JVS(269)*X(53)
  X(55) = X(55)-JVS(272)*X(49)-JVS(273)*X(53)
  X(56) = X(56)-JVS(276)*X(49)-JVS(277)*X(53)
  X(57) = X(57)-JVS(281)*X(49)-JVS(282)*X(53)
  X(60) = X(60)-JVS(294)*X(48)
  X(61) = X(61)-JVS(300)*X(49)-JVS(301)*X(53)
  X(62) = X(62)-JVS(307)*X(53)
  X(63) = X(63)-JVS(313)*X(49)-JVS(314)*X(53)-JVS(315)*X(54)-JVS(316)*X(55)-JVS(317)*X(56)
  X(64) = X(64)-JVS(325)*X(52)-JVS(326)*X(54)-JVS(327)*X(55)-JVS(328)*X(57)-JVS(329)*X(58)-JVS(330)*X(63)
  X(65) = X(65)-JVS(350)*X(39)-JVS(351)*X(41)-JVS(352)*X(56)-JVS(353)*X(60)-JVS(354)*X(61)-JVS(355)*X(62)-JVS(356)*X(63)
  X(67) = X(67)-JVS(377)*X(36)-JVS(378)*X(43)-JVS(379)*X(44)-JVS(380)*X(45)-JVS(381)*X(52)
  X(69) = X(69)-JVS(397)*X(56)-JVS(398)*X(62)
  X(72) = X(72)-JVS(419)*X(43)-JVS(420)*X(44)-JVS(421)*X(52)-JVS(422)*X(54)-JVS(423)*X(55)-JVS(424)*X(67)-JVS(425)*X(70)&
            &-JVS(426)*X(71)
  X(73) = X(73)-JVS(449)*X(70)
  X(76) = X(76)-JVS(466)*X(49)-JVS(467)*X(53)-JVS(468)*X(54)-JVS(469)*X(55)-JVS(470)*X(57)-JVS(471)*X(58)-JVS(472)*X(62)&
            &-JVS(473)*X(66)-JVS(474)*X(69)-JVS(475)*X(74)-JVS(476)*X(75)
  X(77) = X(77)-JVS(490)*X(70)
  X(78) = X(78)-JVS(497)*X(30)-JVS(498)*X(38)-JVS(499)*X(43)-JVS(500)*X(44)-JVS(501)*X(52)-JVS(502)*X(66)-JVS(503)*X(68)&
            &-JVS(504)*X(71)-JVS(505)*X(77)
  X(80) = X(80)-JVS(526)*X(70)-JVS(527)*X(79)
  X(81) = X(81)-JVS(533)*X(38)-JVS(534)*X(43)-JVS(535)*X(44)-JVS(536)*X(46)-JVS(537)*X(47)-JVS(538)*X(51)-JVS(539)*X(52)&
            &-JVS(540)*X(58)-JVS(541)*X(66)-JVS(542)*X(67)-JVS(543)*X(68)-JVS(544)*X(70)-JVS(545)*X(71)-JVS(546)*X(73)&
            &-JVS(547)*X(74)-JVS(548)*X(75)-JVS(549)*X(76)-JVS(550)*X(77)-JVS(551)*X(79)-JVS(552)*X(80)
  X(82) = X(82)-JVS(571)*X(45)-JVS(572)*X(71)-JVS(573)*X(74)-JVS(574)*X(75)-JVS(575)*X(77)-JVS(576)*X(79)
  X(83) = X(83)-JVS(585)*X(24)-JVS(586)*X(66)-JVS(587)*X(68)-JVS(588)*X(70)-JVS(589)*X(71)-JVS(590)*X(73)-JVS(591)*X(74)&
            &-JVS(592)*X(75)-JVS(593)*X(79)-JVS(594)*X(80)
  X(84) = X(84)-JVS(601)*X(36)-JVS(602)*X(43)-JVS(603)*X(44)-JVS(604)*X(52)-JVS(605)*X(54)-JVS(606)*X(55)-JVS(607)*X(57)&
            &-JVS(608)*X(59)-JVS(609)*X(66)-JVS(610)*X(68)-JVS(611)*X(71)-JVS(612)*X(73)-JVS(613)*X(74)-JVS(614)*X(75)&
            &-JVS(615)*X(77)-JVS(616)*X(79)-JVS(617)*X(80)-JVS(618)*X(82)-JVS(619)*X(83)
  X(85) = X(85)-JVS(631)*X(43)-JVS(632)*X(44)-JVS(633)*X(52)-JVS(634)*X(68)-JVS(635)*X(71)-JVS(636)*X(73)-JVS(637)*X(77)&
            &-JVS(638)*X(79)-JVS(639)*X(80)-JVS(640)*X(82)-JVS(641)*X(83)
  X(86) = X(86)-JVS(652)*X(43)-JVS(653)*X(52)-JVS(654)*X(53)-JVS(655)*X(70)-JVS(656)*X(71)-JVS(657)*X(74)-JVS(658)*X(75)&
            &-JVS(659)*X(77)-JVS(660)*X(79)-JVS(661)*X(80)-JVS(662)*X(82)-JVS(663)*X(83)
  X(87) = X(87)-JVS(676)*X(57)-JVS(677)*X(58)-JVS(678)*X(66)-JVS(679)*X(68)-JVS(680)*X(70)-JVS(681)*X(71)-JVS(682)*X(73)&
            &-JVS(683)*X(74)-JVS(684)*X(75)-JVS(685)*X(77)-JVS(686)*X(79)-JVS(687)*X(80)-JVS(688)*X(83)
  X(88) = X(88)-JVS(699)*X(36)-JVS(700)*X(43)-JVS(701)*X(44)-JVS(702)*X(49)-JVS(703)*X(52)-JVS(704)*X(53)-JVS(705)*X(68)&
            &-JVS(706)*X(70)-JVS(707)*X(71)-JVS(708)*X(74)-JVS(709)*X(75)-JVS(710)*X(77)-JVS(711)*X(79)-JVS(712)*X(80)&
            &-JVS(713)*X(82)-JVS(714)*X(83)-JVS(715)*X(84)-JVS(716)*X(85)-JVS(717)*X(86)-JVS(718)*X(87)
  X(89) = X(89)-JVS(731)*X(32)-JVS(732)*X(73)-JVS(733)*X(74)-JVS(734)*X(75)-JVS(735)*X(76)-JVS(736)*X(77)-JVS(737)*X(79)&
            &-JVS(738)*X(80)-JVS(739)*X(83)-JVS(740)*X(84)-JVS(741)*X(85)-JVS(742)*X(86)-JVS(743)*X(87)-JVS(744)*X(88)
  X(90) = X(90)-JVS(756)*X(34)-JVS(757)*X(70)-JVS(758)*X(73)-JVS(759)*X(77)-JVS(760)*X(79)-JVS(761)*X(80)-JVS(762)*X(83)&
            &-JVS(763)*X(87)-JVS(764)*X(88)-JVS(765)*X(89)
  X(91) = X(91)-JVS(776)*X(35)-JVS(777)*X(38)-JVS(778)*X(40)-JVS(779)*X(41)-JVS(780)*X(42)-JVS(781)*X(46)-JVS(782)*X(47)&
            &-JVS(783)*X(48)-JVS(784)*X(49)-JVS(785)*X(50)-JVS(786)*X(51)-JVS(787)*X(53)-JVS(788)*X(54)-JVS(789)*X(55)&
            &-JVS(790)*X(57)-JVS(791)*X(58)-JVS(792)*X(59)-JVS(793)*X(60)-JVS(794)*X(63)-JVS(795)*X(64)-JVS(796)*X(66)&
            &-JVS(797)*X(68)-JVS(798)*X(69)-JVS(799)*X(70)-JVS(800)*X(71)-JVS(801)*X(73)-JVS(802)*X(74)-JVS(803)*X(75)&
            &-JVS(804)*X(76)-JVS(805)*X(77)-JVS(806)*X(78)-JVS(807)*X(79)-JVS(808)*X(80)-JVS(809)*X(81)-JVS(810)*X(82)&
            &-JVS(811)*X(83)-JVS(812)*X(84)-JVS(813)*X(85)-JVS(814)*X(86)-JVS(815)*X(87)-JVS(816)*X(88)-JVS(817)*X(89)&
            &-JVS(818)*X(90)
  X(92) = X(92)-JVS(828)*X(31)-JVS(829)*X(32)-JVS(830)*X(33)-JVS(831)*X(34)-JVS(832)*X(39)-JVS(833)*X(42)-JVS(834)*X(45)&
            &-JVS(835)*X(47)-JVS(836)*X(48)-JVS(837)*X(50)-JVS(838)*X(60)-JVS(839)*X(65)-JVS(840)*X(69)-JVS(841)*X(70)&
            &-JVS(842)*X(72)-JVS(843)*X(73)-JVS(844)*X(74)-JVS(845)*X(75)-JVS(846)*X(76)-JVS(847)*X(77)-JVS(848)*X(78)&
            &-JVS(849)*X(79)-JVS(850)*X(80)-JVS(851)*X(81)-JVS(852)*X(82)-JVS(853)*X(83)-JVS(854)*X(84)-JVS(855)*X(85)&
            &-JVS(856)*X(86)-JVS(857)*X(87)-JVS(858)*X(88)-JVS(859)*X(89)-JVS(860)*X(90)-JVS(861)*X(91)
  X(93) = X(93)-JVS(870)*X(29)-JVS(871)*X(45)-JVS(872)*X(46)-JVS(873)*X(52)-JVS(874)*X(66)-JVS(875)*X(67)-JVS(876)*X(68)&
            &-JVS(877)*X(70)-JVS(878)*X(71)-JVS(879)*X(74)-JVS(880)*X(75)-JVS(881)*X(78)-JVS(882)*X(79)-JVS(883)*X(80)&
            &-JVS(884)*X(82)-JVS(885)*X(83)-JVS(886)*X(84)-JVS(887)*X(85)-JVS(888)*X(86)-JVS(889)*X(87)-JVS(890)*X(88)&
            &-JVS(891)*X(89)-JVS(892)*X(90)-JVS(893)*X(91)-JVS(894)*X(92)
  X(94) = X(94)-JVS(902)*X(30)-JVS(903)*X(36)-JVS(904)*X(38)-JVS(905)*X(43)-JVS(906)*X(44)-JVS(907)*X(49)-JVS(908)*X(52)&
            &-JVS(909)*X(53)-JVS(910)*X(54)-JVS(911)*X(55)-JVS(912)*X(56)-JVS(913)*X(57)-JVS(914)*X(58)-JVS(915)*X(59)&
            &-JVS(916)*X(62)-JVS(917)*X(66)-JVS(918)*X(68)-JVS(919)*X(69)-JVS(920)*X(70)-JVS(921)*X(71)-JVS(922)*X(73)&
            &-JVS(923)*X(74)-JVS(924)*X(75)-JVS(925)*X(77)-JVS(926)*X(79)-JVS(927)*X(80)-JVS(928)*X(82)-JVS(929)*X(83)&
            &-JVS(930)*X(84)-JVS(931)*X(85)-JVS(932)*X(86)-JVS(933)*X(87)-JVS(934)*X(88)-JVS(935)*X(89)-JVS(936)*X(90)&
            &-JVS(937)*X(91)-JVS(938)*X(92)-JVS(939)*X(93)
  X(95) = X(95)-JVS(946)*X(33)-JVS(947)*X(61)-JVS(948)*X(79)-JVS(949)*X(83)-JVS(950)*X(87)-JVS(951)*X(88)-JVS(952)*X(89)&
            &-JVS(953)*X(90)-JVS(954)*X(91)-JVS(955)*X(92)-JVS(956)*X(93)-JVS(957)*X(94)
  X(96) = X(96)-JVS(963)*X(31)-JVS(964)*X(37)-JVS(965)*X(49)-JVS(966)*X(52)-JVS(967)*X(54)-JVS(968)*X(55)-JVS(969)*X(63)&
            &-JVS(970)*X(67)-JVS(971)*X(71)-JVS(972)*X(73)-JVS(973)*X(74)-JVS(974)*X(75)-JVS(975)*X(77)-JVS(976)*X(78)&
            &-JVS(977)*X(79)-JVS(978)*X(80)-JVS(979)*X(82)-JVS(980)*X(83)-JVS(981)*X(84)-JVS(982)*X(85)-JVS(983)*X(86)&
            &-JVS(984)*X(87)-JVS(985)*X(88)-JVS(986)*X(89)-JVS(987)*X(90)-JVS(988)*X(91)-JVS(989)*X(92)-JVS(990)*X(93)&
            &-JVS(991)*X(94)-JVS(992)*X(95)
  X(97) = X(97)-JVS(997)*X(42)-JVS(998)*X(47)-JVS(999)*X(72)-JVS(1000)*X(74)-JVS(1001)*X(75)-JVS(1002)*X(79)-JVS(1003)&
            &*X(80)-JVS(1004)*X(81)-JVS(1005)*X(82)-JVS(1006)*X(83)-JVS(1007)*X(85)-JVS(1008)*X(86)-JVS(1009)*X(87)&
            &-JVS(1010)*X(88)-JVS(1011)*X(89)-JVS(1012)*X(90)-JVS(1013)*X(91)-JVS(1014)*X(92)-JVS(1015)*X(93)-JVS(1016)&
            &*X(94)-JVS(1017)*X(95)-JVS(1018)*X(96)
  X(98) = X(98)-JVS(1022)*X(3)-JVS(1023)*X(5)-JVS(1024)*X(24)-JVS(1025)*X(25)-JVS(1026)*X(26)-JVS(1027)*X(27)-JVS(1028)&
            &*X(28)-JVS(1029)*X(29)-JVS(1030)*X(30)-JVS(1031)*X(35)-JVS(1032)*X(36)-JVS(1033)*X(38)-JVS(1034)*X(40)&
            &-JVS(1035)*X(41)-JVS(1036)*X(42)-JVS(1037)*X(43)-JVS(1038)*X(44)-JVS(1039)*X(46)-JVS(1040)*X(49)-JVS(1041)&
            &*X(50)-JVS(1042)*X(51)-JVS(1043)*X(52)-JVS(1044)*X(53)-JVS(1045)*X(54)-JVS(1046)*X(55)-JVS(1047)*X(56)&
            &-JVS(1048)*X(57)-JVS(1049)*X(58)-JVS(1050)*X(59)-JVS(1051)*X(61)-JVS(1052)*X(62)-JVS(1053)*X(63)-JVS(1054)&
            &*X(64)-JVS(1055)*X(65)-JVS(1056)*X(66)-JVS(1057)*X(67)-JVS(1058)*X(68)-JVS(1059)*X(69)-JVS(1060)*X(70)&
            &-JVS(1061)*X(71)-JVS(1062)*X(73)-JVS(1063)*X(74)-JVS(1064)*X(75)-JVS(1065)*X(76)-JVS(1066)*X(77)-JVS(1067)&
            &*X(78)-JVS(1068)*X(79)-JVS(1069)*X(80)-JVS(1070)*X(81)-JVS(1071)*X(82)-JVS(1072)*X(83)-JVS(1073)*X(84)&
            &-JVS(1074)*X(85)-JVS(1075)*X(86)-JVS(1076)*X(87)-JVS(1077)*X(88)-JVS(1078)*X(89)-JVS(1079)*X(90)-JVS(1080)&
            &*X(91)-JVS(1081)*X(92)-JVS(1082)*X(93)-JVS(1083)*X(94)-JVS(1084)*X(95)-JVS(1085)*X(96)-JVS(1086)*X(97)
  X(99) = X(99)-JVS(1089)*X(39)-JVS(1090)*X(41)-JVS(1091)*X(50)-JVS(1092)*X(56)-JVS(1093)*X(60)-JVS(1094)*X(61)&
            &-JVS(1095)*X(62)-JVS(1096)*X(63)-JVS(1097)*X(65)-JVS(1098)*X(66)-JVS(1099)*X(68)-JVS(1100)*X(69)-JVS(1101)&
            &*X(70)-JVS(1102)*X(71)-JVS(1103)*X(72)-JVS(1104)*X(73)-JVS(1105)*X(74)-JVS(1106)*X(75)-JVS(1107)*X(76)&
            &-JVS(1108)*X(77)-JVS(1109)*X(78)-JVS(1110)*X(79)-JVS(1111)*X(80)-JVS(1112)*X(81)-JVS(1113)*X(82)-JVS(1114)&
            &*X(83)-JVS(1115)*X(84)-JVS(1116)*X(85)-JVS(1117)*X(86)-JVS(1118)*X(87)-JVS(1119)*X(88)-JVS(1120)*X(89)&
            &-JVS(1121)*X(90)-JVS(1122)*X(91)-JVS(1123)*X(92)-JVS(1124)*X(93)-JVS(1125)*X(94)-JVS(1126)*X(95)-JVS(1127)&
            &*X(96)-JVS(1128)*X(97)-JVS(1129)*X(98)
  X(99) = X(99)/JVS(1130)
  X(98) = (X(98)-JVS(1088)*X(99))/(JVS(1087))
  X(97) = (X(97)-JVS(1020)*X(98)-JVS(1021)*X(99))/(JVS(1019))
  X(96) = (X(96)-JVS(994)*X(97)-JVS(995)*X(98)-JVS(996)*X(99))/(JVS(993))
  X(95) = (X(95)-JVS(959)*X(96)-JVS(960)*X(97)-JVS(961)*X(98)-JVS(962)*X(99))/(JVS(958))
  X(94) = (X(94)-JVS(941)*X(95)-JVS(942)*X(96)-JVS(943)*X(97)-JVS(944)*X(98)-JVS(945)*X(99))/(JVS(940))
  X(93) = (X(93)-JVS(896)*X(94)-JVS(897)*X(95)-JVS(898)*X(96)-JVS(899)*X(97)-JVS(900)*X(98)-JVS(901)*X(99))/(JVS(895))
  X(92) = (X(92)-JVS(863)*X(93)-JVS(864)*X(94)-JVS(865)*X(95)-JVS(866)*X(96)-JVS(867)*X(97)-JVS(868)*X(98)-JVS(869)&
            &*X(99))/(JVS(862))
  X(91) = (X(91)-JVS(820)*X(92)-JVS(821)*X(93)-JVS(822)*X(94)-JVS(823)*X(95)-JVS(824)*X(96)-JVS(825)*X(97)-JVS(826)&
            &*X(98)-JVS(827)*X(99))/(JVS(819))
  X(90) = (X(90)-JVS(767)*X(91)-JVS(768)*X(92)-JVS(769)*X(93)-JVS(770)*X(94)-JVS(771)*X(95)-JVS(772)*X(96)-JVS(773)&
            &*X(97)-JVS(774)*X(98)-JVS(775)*X(99))/(JVS(766))
  X(89) = (X(89)-JVS(746)*X(90)-JVS(747)*X(91)-JVS(748)*X(92)-JVS(749)*X(93)-JVS(750)*X(94)-JVS(751)*X(95)-JVS(752)&
            &*X(96)-JVS(753)*X(97)-JVS(754)*X(98)-JVS(755)*X(99))/(JVS(745))
  X(88) = (X(88)-JVS(720)*X(89)-JVS(721)*X(90)-JVS(722)*X(91)-JVS(723)*X(92)-JVS(724)*X(93)-JVS(725)*X(94)-JVS(726)&
            &*X(95)-JVS(727)*X(96)-JVS(728)*X(97)-JVS(729)*X(98)-JVS(730)*X(99))/(JVS(719))
  X(87) = (X(87)-JVS(690)*X(89)-JVS(691)*X(90)-JVS(692)*X(91)-JVS(693)*X(92)-JVS(694)*X(95)-JVS(695)*X(96)-JVS(696)&
            &*X(97)-JVS(697)*X(98)-JVS(698)*X(99))/(JVS(689))
  X(86) = (X(86)-JVS(665)*X(87)-JVS(666)*X(88)-JVS(667)*X(89)-JVS(668)*X(92)-JVS(669)*X(93)-JVS(670)*X(94)-JVS(671)&
            &*X(95)-JVS(672)*X(96)-JVS(673)*X(97)-JVS(674)*X(98)-JVS(675)*X(99))/(JVS(664))
  X(85) = (X(85)-JVS(643)*X(86)-JVS(644)*X(87)-JVS(645)*X(88)-JVS(646)*X(92)-JVS(647)*X(93)-JVS(648)*X(94)-JVS(649)&
            &*X(97)-JVS(650)*X(98)-JVS(651)*X(99))/(JVS(642))
  X(84) = (X(84)-JVS(621)*X(85)-JVS(622)*X(86)-JVS(623)*X(87)-JVS(624)*X(88)-JVS(625)*X(91)-JVS(626)*X(92)-JVS(627)&
            &*X(94)-JVS(628)*X(97)-JVS(629)*X(98)-JVS(630)*X(99))/(JVS(620))
  X(83) = (X(83)-JVS(596)*X(87)-JVS(597)*X(92)-JVS(598)*X(97)-JVS(599)*X(98)-JVS(600)*X(99))/(JVS(595))
  X(82) = (X(82)-JVS(578)*X(83)-JVS(579)*X(87)-JVS(580)*X(88)-JVS(581)*X(92)-JVS(582)*X(97)-JVS(583)*X(98)-JVS(584)&
            &*X(99))/(JVS(577))
  X(81) = (X(81)-JVS(554)*X(82)-JVS(555)*X(83)-JVS(556)*X(85)-JVS(557)*X(86)-JVS(558)*X(87)-JVS(559)*X(88)-JVS(560)&
            &*X(89)-JVS(561)*X(90)-JVS(562)*X(91)-JVS(563)*X(92)-JVS(564)*X(93)-JVS(565)*X(94)-JVS(566)*X(95)-JVS(567)*X(96)&
            &-JVS(568)*X(97)-JVS(569)*X(98)-JVS(570)*X(99))/(JVS(553))
  X(80) = (X(80)-JVS(529)*X(83)-JVS(530)*X(87)-JVS(531)*X(98)-JVS(532)*X(99))/(JVS(528))
  X(79) = (X(79)-JVS(522)*X(83)-JVS(523)*X(87)-JVS(524)*X(98)-JVS(525)*X(99))/(JVS(521))
  X(78) = (X(78)-JVS(507)*X(79)-JVS(508)*X(82)-JVS(509)*X(83)-JVS(510)*X(84)-JVS(511)*X(85)-JVS(512)*X(86)-JVS(513)&
            &*X(87)-JVS(514)*X(89)-JVS(515)*X(90)-JVS(516)*X(95)-JVS(517)*X(96)-JVS(518)*X(97)-JVS(519)*X(98)-JVS(520)&
            &*X(99))/(JVS(506))
  X(77) = (X(77)-JVS(492)*X(79)-JVS(493)*X(83)-JVS(494)*X(87)-JVS(495)*X(98)-JVS(496)*X(99))/(JVS(491))
  X(76) = (X(76)-JVS(478)*X(77)-JVS(479)*X(83)-JVS(480)*X(87)-JVS(481)*X(89)-JVS(482)*X(90)-JVS(483)*X(91)-JVS(484)&
            &*X(92)-JVS(485)*X(95)-JVS(486)*X(96)-JVS(487)*X(97)-JVS(488)*X(98)-JVS(489)*X(99))/(JVS(477))
  X(75) = (X(75)-JVS(462)*X(83)-JVS(463)*X(87)-JVS(464)*X(98)-JVS(465)*X(99))/(JVS(461))
  X(74) = (X(74)-JVS(457)*X(83)-JVS(458)*X(87)-JVS(459)*X(98)-JVS(460)*X(99))/(JVS(456))
  X(73) = (X(73)-JVS(451)*X(79)-JVS(452)*X(83)-JVS(453)*X(87)-JVS(454)*X(98)-JVS(455)*X(99))/(JVS(450))
  X(72) = (X(72)-JVS(428)*X(74)-JVS(429)*X(75)-JVS(430)*X(79)-JVS(431)*X(80)-JVS(432)*X(82)-JVS(433)*X(83)-JVS(434)&
            &*X(85)-JVS(435)*X(86)-JVS(436)*X(87)-JVS(437)*X(88)-JVS(438)*X(89)-JVS(439)*X(90)-JVS(440)*X(91)-JVS(441)*X(92)&
            &-JVS(442)*X(93)-JVS(443)*X(94)-JVS(444)*X(95)-JVS(445)*X(96)-JVS(446)*X(97)-JVS(447)*X(98)-JVS(448)*X(99))&
            &/(JVS(427))
  X(71) = (X(71)-JVS(415)*X(83)-JVS(416)*X(87)-JVS(417)*X(98)-JVS(418)*X(99))/(JVS(414))
  X(70) = (X(70)-JVS(410)*X(83)-JVS(411)*X(87)-JVS(412)*X(98)-JVS(413)*X(99))/(JVS(409))
  X(69) = (X(69)-JVS(400)*X(89)-JVS(401)*X(90)-JVS(402)*X(91)-JVS(403)*X(92)-JVS(404)*X(95)-JVS(405)*X(96)-JVS(406)&
            &*X(97)-JVS(407)*X(98)-JVS(408)*X(99))/(JVS(399))
  X(68) = (X(68)-JVS(393)*X(83)-JVS(394)*X(87)-JVS(395)*X(98)-JVS(396)*X(99))/(JVS(392))
  X(67) = (X(67)-JVS(383)*X(71)-JVS(384)*X(74)-JVS(385)*X(75)-JVS(386)*X(79)-JVS(387)*X(82)-JVS(388)*X(87)-JVS(389)&
            &*X(92)-JVS(390)*X(98)-JVS(391)*X(99))/(JVS(382))
  X(66) = (X(66)-JVS(373)*X(83)-JVS(374)*X(87)-JVS(375)*X(98)-JVS(376)*X(99))/(JVS(372))
  X(65) = (X(65)-JVS(358)*X(69)-JVS(359)*X(73)-JVS(360)*X(76)-JVS(361)*X(77)-JVS(362)*X(78)-JVS(363)*X(79)-JVS(364)&
            &*X(80)-JVS(365)*X(81)-JVS(366)*X(84)-JVS(367)*X(87)-JVS(368)*X(91)-JVS(369)*X(92)-JVS(370)*X(98)-JVS(371)&
            &*X(99))/(JVS(357))
  X(64) = (X(64)-JVS(332)*X(66)-JVS(333)*X(68)-JVS(334)*X(70)-JVS(335)*X(71)-JVS(336)*X(73)-JVS(337)*X(74)-JVS(338)&
            &*X(75)-JVS(339)*X(76)-JVS(340)*X(77)-JVS(341)*X(78)-JVS(342)*X(79)-JVS(343)*X(80)-JVS(344)*X(81)-JVS(345)*X(83)&
            &-JVS(346)*X(84)-JVS(347)*X(87)-JVS(348)*X(98)-JVS(349)*X(99))/(JVS(331))
  X(63) = (X(63)-JVS(319)*X(73)-JVS(320)*X(77)-JVS(321)*X(80)-JVS(322)*X(87)-JVS(323)*X(98)-JVS(324)*X(99))/(JVS(318))
  X(62) = (X(62)-JVS(309)*X(69)-JVS(310)*X(91)-JVS(311)*X(98)-JVS(312)*X(99))/(JVS(308))
  X(61) = (X(61)-JVS(303)*X(79)-JVS(304)*X(87)-JVS(305)*X(98)-JVS(306)*X(99))/(JVS(302))
  X(60) = (X(60)-JVS(296)*X(69)-JVS(297)*X(91)-JVS(298)*X(92)-JVS(299)*X(99))/(JVS(295))
  X(59) = (X(59)-JVS(290)*X(88)-JVS(291)*X(91)-JVS(292)*X(94)-JVS(293)*X(98))/(JVS(289))
  X(58) = (X(58)-JVS(287)*X(87)-JVS(288)*X(98))/(JVS(286))
  X(57) = (X(57)-JVS(284)*X(87)-JVS(285)*X(98))/(JVS(283))
  X(56) = (X(56)-JVS(279)*X(98)-JVS(280)*X(99))/(JVS(278))
  X(55) = (X(55)-JVS(275)*X(98))/(JVS(274))
  X(54) = (X(54)-JVS(271)*X(98))/(JVS(270))
  X(53) = (X(53)-JVS(267)*X(98))/(JVS(266))
  X(52) = (X(52)-JVS(265)*X(98))/(JVS(264))
  X(51) = (X(51)-JVS(260)*X(88)-JVS(261)*X(93)-JVS(262)*X(94)-JVS(263)*X(98))/(JVS(259))
  X(50) = (X(50)-JVS(256)*X(91)-JVS(257)*X(92)-JVS(258)*X(98))/(JVS(255))
  X(49) = (X(49)-JVS(254)*X(98))/(JVS(253))
  X(48) = (X(48)-JVS(249)*X(60)-JVS(250)*X(91)-JVS(251)*X(92)-JVS(252)*X(99))/(JVS(248))
  X(47) = (X(47)-JVS(245)*X(81)-JVS(246)*X(91)-JVS(247)*X(97))/(JVS(244))
  X(46) = (X(46)-JVS(241)*X(91)-JVS(242)*X(93)-JVS(243)*X(98))/(JVS(240))
  X(45) = (X(45)-JVS(238)*X(92)-JVS(239)*X(98))/(JVS(237))
  X(44) = (X(44)-JVS(235)*X(98))/(JVS(234))
  X(43) = (X(43)-JVS(233)*X(98))/(JVS(232))
  X(42) = (X(42)-JVS(230)*X(97)-JVS(231)*X(98))/(JVS(229))
  X(41) = (X(41)-JVS(227)*X(98)-JVS(228)*X(99))/(JVS(226))
  X(40) = (X(40)-JVS(223)*X(41)-JVS(224)*X(98)-JVS(225)*X(99))/(JVS(222))
  X(39) = (X(39)-JVS(220)*X(92)-JVS(221)*X(99))/(JVS(219))
  X(38) = (X(38)-JVS(218)*X(98))/(JVS(217))
  X(37) = (X(37)-JVS(212)*X(49)-JVS(213)*X(74)-JVS(214)*X(75)-JVS(215)*X(87)-JVS(216)*X(98))/(JVS(211))
  X(36) = (X(36)-JVS(210)*X(98))/(JVS(209))
  X(35) = (X(35)-JVS(207)*X(91)-JVS(208)*X(98))/(JVS(206))
  X(34) = (X(34)-JVS(204)*X(90)-JVS(205)*X(92))/(JVS(203))
  X(33) = (X(33)-JVS(201)*X(92)-JVS(202)*X(95))/(JVS(200))
  X(32) = (X(32)-JVS(198)*X(89)-JVS(199)*X(92))/(JVS(197))
  X(31) = (X(31)-JVS(195)*X(92)-JVS(196)*X(96))/(JVS(194))
  X(30) = (X(30)-JVS(193)*X(98))/(JVS(192))
  X(29) = (X(29)-JVS(191)*X(98))/(JVS(190))
  X(28) = (X(28)-JVS(189)*X(98))/(JVS(188))
  X(27) = (X(27)-JVS(186)*X(98))/(JVS(185))
  X(26) = (X(26)-JVS(184)*X(98))/(JVS(183))
  X(25) = (X(25)-JVS(181)*X(98))/(JVS(180))
  X(24) = (X(24)-JVS(179)*X(87))/(JVS(178))
  X(23) = (X(23)-JVS(177)*X(98))/(JVS(176))
  X(22) = (X(22)-JVS(173)*X(23)-JVS(174)*X(27)-JVS(175)*X(98))/(JVS(172))
  X(21) = (X(21)-JVS(171)*X(98))/(JVS(170))
  X(20) = (X(20)-JVS(167)*X(21)-JVS(168)*X(25)-JVS(169)*X(98))/(JVS(166))
  X(19) = (X(19)-JVS(157)*X(20)-JVS(158)*X(21)-JVS(159)*X(22)-JVS(160)*X(23)-JVS(161)*X(25)-JVS(162)*X(26)-JVS(163)&
            &*X(27)-JVS(164)*X(28)-JVS(165)*X(98))/(JVS(156))
  X(18) = (X(18)-JVS(104)*X(29)-JVS(105)*X(30)-JVS(106)*X(36)-JVS(107)*X(38)-JVS(108)*X(40)-JVS(109)*X(42)-JVS(110)&
            &*X(43)-JVS(111)*X(44)-JVS(112)*X(46)-JVS(113)*X(49)-JVS(114)*X(51)-JVS(115)*X(52)-JVS(116)*X(53)-JVS(117)*X(54)&
            &-JVS(118)*X(55)-JVS(119)*X(56)-JVS(120)*X(57)-JVS(121)*X(58)-JVS(122)*X(59)-JVS(123)*X(61)-JVS(124)*X(62)&
            &-JVS(125)*X(63)-JVS(126)*X(64)-JVS(127)*X(65)-JVS(128)*X(66)-JVS(129)*X(67)-JVS(130)*X(68)-JVS(131)*X(70)&
            &-JVS(132)*X(71)-JVS(133)*X(73)-JVS(134)*X(74)-JVS(135)*X(75)-JVS(136)*X(76)-JVS(137)*X(77)-JVS(138)*X(78)&
            &-JVS(139)*X(79)-JVS(140)*X(80)-JVS(141)*X(81)-JVS(142)*X(82)-JVS(143)*X(84)-JVS(144)*X(85)-JVS(145)*X(86)&
            &-JVS(146)*X(87)-JVS(147)*X(91)-JVS(148)*X(92)-JVS(149)*X(97)-JVS(150)*X(98)-JVS(151)*X(99))/(JVS(103))
  X(17) = (X(17)-JVS(97)*X(70)-JVS(98)*X(74)-JVS(99)*X(75)-JVS(100)*X(87)-JVS(101)*X(98)-JVS(102)*X(99))/(JVS(96))
  X(16) = (X(16)-JVS(89)*X(43)-JVS(90)*X(49)-JVS(91)*X(52)-JVS(92)*X(53)-JVS(93)*X(71)-JVS(94)*X(79)-JVS(95)*X(98))&
            &/(JVS(88))
  X(15) = (X(15)-JVS(83)*X(72)-JVS(84)*X(88)-JVS(85)*X(91)-JVS(86)*X(93)-JVS(87)*X(94))/(JVS(82))
  X(14) = (X(14)-JVS(76)*X(72)-JVS(77)*X(88)-JVS(78)*X(93)-JVS(79)*X(94)-JVS(80)*X(97)-JVS(81)*X(99))/(JVS(75))
  X(13) = (X(13)-JVS(67)*X(48)-JVS(68)*X(61)-JVS(69)*X(68)-JVS(70)*X(83)-JVS(71)*X(87)-JVS(72)*X(92)-JVS(73)*X(98)&
            &-JVS(74)*X(99))/(JVS(66))
  X(12) = (X(12)-JVS(62)*X(48)-JVS(63)*X(68)-JVS(64)*X(92)-JVS(65)*X(99))/(JVS(61))
  X(11) = (X(11)-JVS(57)*X(89)-JVS(58)*X(90)-JVS(59)*X(91)-JVS(60)*X(95))/(JVS(56))
  X(10) = (X(10)-JVS(54)*X(91)-JVS(55)*X(96))/(JVS(53))
  X(9) = (X(9)-JVS(39)*X(70)-JVS(40)*X(71)-JVS(41)*X(74)-JVS(42)*X(75)-JVS(43)*X(77)-JVS(44)*X(79)-JVS(45)*X(87)-JVS(46)&
           &*X(88)-JVS(47)*X(89)-JVS(48)*X(90)-JVS(49)*X(91)-JVS(50)*X(93)-JVS(51)*X(94)-JVS(52)*X(95))/(JVS(38))
  X(8) = (X(8)-JVS(29)*X(68)-JVS(30)*X(71)-JVS(31)*X(79)-JVS(32)*X(87)-JVS(33)*X(88)-JVS(34)*X(91)-JVS(35)*X(93)-JVS(36)&
           &*X(94)-JVS(37)*X(96))/(JVS(28))
  X(7) = (X(7)-JVS(13)*X(47)-JVS(14)*X(58)-JVS(15)*X(66)-JVS(16)*X(68)-JVS(17)*X(70)-JVS(18)*X(71)-JVS(19)*X(73)-JVS(20)&
           &*X(74)-JVS(21)*X(75)-JVS(22)*X(77)-JVS(23)*X(79)-JVS(24)*X(80)-JVS(25)*X(87)-JVS(26)*X(97)-JVS(27)*X(98))&
           &/(JVS(12))
  X(6) = X(6)/JVS(11)
  X(5) = X(5)/JVS(10)
  X(4) = X(4)/JVS(9)
  X(3) = X(3)/JVS(8)
  X(2) = (X(2)-JVS(5)*X(58)-JVS(6)*X(68)-JVS(7)*X(87))/(JVS(4))
  X(1) = (X(1)-JVS(2)*X(40)-JVS(3)*X(98))/(JVS(1))
      
END SUBROUTINE saprc99_mosaic_8bin_vbs2_KppSolve
























      SUBROUTINE saprc99_mosaic_8bin_vbs2_WCOPY(N,X,incX,Y,incY)








      
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

      END SUBROUTINE saprc99_mosaic_8bin_vbs2_WCOPY



      SUBROUTINE saprc99_mosaic_8bin_vbs2_WAXPY(N,Alpha,X,incX,Y,incY)









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
      
      END SUBROUTINE saprc99_mosaic_8bin_vbs2_WAXPY




      SUBROUTINE saprc99_mosaic_8bin_vbs2_WSCAL(N,Alpha,X,incX)









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

      END SUBROUTINE saprc99_mosaic_8bin_vbs2_WSCAL


      REAL(kind=dp) FUNCTION saprc99_mosaic_8bin_vbs2_WLAMCH( C )








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
          CALL saprc99_mosaic_8bin_vbs2_WLAMCH_ADD(ONE,Eps,Sum)
          IF (Sum.LE.ONE) GOTO 10
        END DO
        PRINT*,'ERROR IN WLAMCH. EPS < ',Eps
        RETURN
10      Eps = Eps*2
        i = i-1      
      END IF

      saprc99_mosaic_8bin_vbs2_WLAMCH = Eps

      END FUNCTION saprc99_mosaic_8bin_vbs2_WLAMCH
     
      SUBROUTINE saprc99_mosaic_8bin_vbs2_WLAMCH_ADD( A, B, Sum )

      
      REAL(kind=dp) A, B, Sum
      Sum = A + B

      END SUBROUTINE saprc99_mosaic_8bin_vbs2_WLAMCH_ADD




      SUBROUTINE saprc99_mosaic_8bin_vbs2_SET2ZERO(N,Y)




      
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

      END SUBROUTINE saprc99_mosaic_8bin_vbs2_SET2ZERO



      REAL(kind=dp) FUNCTION saprc99_mosaic_8bin_vbs2_WDOT (N, DX, incX, DY, incY) 









      IMPLICIT NONE
      INTEGER :: N, incX, incY
      REAL(kind=dp) :: DX(N), DY(N) 

      INTEGER :: i, IX, IY, M, MP1, NS
                                 
      saprc99_mosaic_8bin_vbs2_WDOT = 0.0D0 
      IF (N .LE. 0) RETURN 
      IF (incX .EQ. incY) IF (incX-1) 5,20,60 



    5 IX = 1 
      IY = 1 
      IF (incX .LT. 0) IX = (-N+1)*incX + 1 
      IF (incY .LT. 0) IY = (-N+1)*incY + 1 
      DO i = 1,N 
        saprc99_mosaic_8bin_vbs2_WDOT = saprc99_mosaic_8bin_vbs2_WDOT + DX(IX)*DY(IY) 
        IX = IX + incX 
        IY = IY + incY 
      END DO 
      RETURN 





   20 M = MOD(N,5) 
      IF (M .EQ. 0) GO TO 40 
      DO i = 1,M 
         saprc99_mosaic_8bin_vbs2_WDOT = saprc99_mosaic_8bin_vbs2_WDOT + DX(i)*DY(i) 
      END DO 
      IF (N .LT. 5) RETURN 
   40 MP1 = M + 1 
      DO i = MP1,N,5 
          saprc99_mosaic_8bin_vbs2_WDOT = saprc99_mosaic_8bin_vbs2_WDOT + DX(i)*DY(i) + DX(i+1)*DY(i+1) +&
                   DX(i+2)*DY(i+2) +  &
                   DX(i+3)*DY(i+3) + DX(i+4)*DY(i+4)                   
      END DO 
      RETURN 



   60 NS = N*incX 
      DO i = 1,NS,incX 
        saprc99_mosaic_8bin_vbs2_WDOT = saprc99_mosaic_8bin_vbs2_WDOT + DX(i)*DY(i) 
      END DO 

      END FUNCTION saprc99_mosaic_8bin_vbs2_WDOT                                          




END MODULE saprc99_mosaic_8bin_vbs2_Integrator
