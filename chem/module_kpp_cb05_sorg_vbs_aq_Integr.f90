
































MODULE cb05_sorg_vbs_aq_Integrator

 USE cb05_sorg_vbs_aq_Parameters
 USE cb05_sorg_vbs_aq_Precision
 USE cb05_sorg_vbs_aq_JacobianSP

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

SUBROUTINE  cb05_sorg_vbs_aq_INTEGRATE( TIN, TOUT, &
  FIX, VAR,  RCONST, ATOL, RTOL, IRR_WRK,  &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U  )

   USE cb05_sorg_vbs_aq_Parameters

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

   CALL cb05_sorg_vbs_aq_Rosenbrock(VAR, FIX, RCONST, TIN,TOUT,   &
         ATOL,RTOL,               &
         RCNTRL,ICNTRL,RSTATUS,ISTATUS,IRR_WRK,IERR)

   STEPMIN = RCNTRL(ihexit)
   
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U))    IERR_U       = IERR

END SUBROUTINE  cb05_sorg_vbs_aq_INTEGRATE


SUBROUTINE  cb05_sorg_vbs_aq_Rosenbrock(Y, FIX, RCONST, Tstart,Tend, &
           AbsTol,RelTol,            &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IRR_WRK,IERR)







































































































  USE cb05_sorg_vbs_aq_Parameters

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
      CALL cb05_sorg_vbs_aq_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
      RETURN
   END IF


   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 100000
   ELSEIF (ICNTRL(4) > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL cb05_sorg_vbs_aq_ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF


   Roundoff = cb05_sorg_vbs_aq_WLAMCH('E')


   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      Hmin = RCNTRL(1)
   ELSE
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL cb05_sorg_vbs_aq_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL cb05_sorg_vbs_aq_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL cb05_sorg_vbs_aq_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2_dp
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL cb05_sorg_vbs_aq_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0_dp
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL cb05_sorg_vbs_aq_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1_dp
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL cb05_sorg_vbs_aq_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9_dp
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL cb05_sorg_vbs_aq_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.0_dp*Roundoff) &
         .OR. (RelTol(i) >= 1.0_dp) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL cb05_sorg_vbs_aq_ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO



   SELECT CASE (Method)
     CASE (1)
       CALL cb05_sorg_vbs_aq_Ros2(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (2)
       CALL cb05_sorg_vbs_aq_Ros3(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (3)
       CALL cb05_sorg_vbs_aq_Ros4(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (4)
       CALL cb05_sorg_vbs_aq_Rodas3(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (5)
       CALL cb05_sorg_vbs_aq_Rodas4(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(4)=', Method
       CALL cb05_sorg_vbs_aq_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       RETURN
   END SELECT


   CALL cb05_sorg_vbs_aq_ros_Integrator(Y,Tstart,Tend,Texit,      &
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



 SUBROUTINE  cb05_sorg_vbs_aq_ros_ErrorMsg(Code,T,H,IERR)



   USE cb05_sorg_vbs_aq_Precision

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

 END SUBROUTINE  cb05_sorg_vbs_aq_ros_ErrorMsg


 SUBROUTINE  cb05_sorg_vbs_aq_ros_Integrator (Y, Tstart, Tend, T,     &
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
      CALL cb05_sorg_vbs_aq_ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1_dp*H) == T).OR.(H <= Roundoff) ) THEN  
      CALL cb05_sorg_vbs_aq_ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF


   Hexit = H
   H = MIN(H,ABS(Tend-T))


   CALL cb05_sorg_vbs_aq_FunTemplate(T,Y,Fcn0, RCONST, FIX, Nfun)
   IF( T == Tstart ) THEN
     CALL cb05_sorg_vbs_aq_IRRFun( Y, FIX, RCONST, IRR_WRK )
   ENDIF


   IF (.NOT.Autonomous) THEN
      CALL cb05_sorg_vbs_aq_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )
   END IF


   CALL cb05_sorg_vbs_aq_JacTemplate(T,Y,Jac0, FIX, Njac, RCONST)


UntilAccepted: DO

   CALL cb05_sorg_vbs_aq_ros_PrepareMatrix(H,Direction,ros_Gamma(1), &
          Jac0,Ghimj,Pivot,Singular, Ndec,  Nsng )
   IF (Singular) THEN 
       CALL cb05_sorg_vbs_aq_ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF


Stage: DO istage = 1, ros_S

      
       ioffset = NVAR*(istage-1)

      
       IF ( istage == 1 ) THEN
         CALL cb05_sorg_vbs_aq_WCOPY(NVAR,Fcn0,1,Fcn,1)
      
       ELSEIF ( ros_NewF(istage) ) THEN
         CALL cb05_sorg_vbs_aq_WCOPY(NVAR,Y,1,Ynew,1)
         DO j = 1, istage-1
           CALL cb05_sorg_vbs_aq_WAXPY(NVAR,ros_A((istage-1)*(istage-2)/2+j), &
            K(NVAR*(j-1)+1),1,Ynew,1)
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL cb05_sorg_vbs_aq_FunTemplate(Tau,Ynew,Fcn, RCONST, FIX, Nfun)
       END IF 
       CALL cb05_sorg_vbs_aq_WCOPY(NVAR,Fcn,1,K(ioffset+1),1)
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         CALL cb05_sorg_vbs_aq_WAXPY(NVAR,HC,K(NVAR*(j-1)+1),1,K(ioffset+1),1)
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         CALL cb05_sorg_vbs_aq_WAXPY(NVAR,HG,dFdT,1,K(ioffset+1),1)
       END IF
       CALL cb05_sorg_vbs_aq_ros_Solve(Ghimj, Pivot, K(ioffset+1), Nsol)

   END DO Stage



   CALL cb05_sorg_vbs_aq_WCOPY(NVAR,Y,1,Ynew,1)
   DO j=1,ros_S
         CALL cb05_sorg_vbs_aq_WAXPY(NVAR,ros_M(j),K(NVAR*(j-1)+1),1,Ynew,1)
   END DO


   CALL cb05_sorg_vbs_aq_WSCAL(NVAR,ZERO,Yerr,1)
   DO j=1,ros_S
        CALL cb05_sorg_vbs_aq_WAXPY(NVAR,ros_E(j),K(NVAR*(j-1)+1),1,Yerr,1)
   END DO
   Err = cb05_sorg_vbs_aq_ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )


   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac


   Nstp = Nstp+1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  
      Nacc = Nacc+1
      CALL cb05_sorg_vbs_aq_WCOPY(NVAR,Ynew,1,Y,1)
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

  END SUBROUTINE  cb05_sorg_vbs_aq_ros_Integrator



  REAL(kind=dp) FUNCTION  cb05_sorg_vbs_aq_ros_ErrorNorm ( Y, Ynew, Yerr, &
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

    cb05_sorg_vbs_aq_ros_ErrorNorm = Err

  END FUNCTION  cb05_sorg_vbs_aq_ros_ErrorNorm



  SUBROUTINE cb05_sorg_vbs_aq_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )



   IMPLICIT NONE


   REAL(kind=dp), INTENT(IN) :: T, Roundoff, Y(NVAR), Fcn0(NVAR)
   REAL(kind=dp), INTENT(IN) :: RCONST(NREACT), FIX(NFIX)

   REAL(kind=dp), INTENT(OUT) :: dFdT(NVAR)

   INTEGER, INTENT(INOUT) ::Nfun

   REAL(kind=dp) :: Delta
   REAL(kind=dp), PARAMETER :: ONE = 1.0_dp, DeltaMin = 1.0E-6_dp

   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL cb05_sorg_vbs_aq_FunTemplate(T+Delta,Y,dFdT, RCONST, FIX, Nfun)
   CALL cb05_sorg_vbs_aq_WAXPY(NVAR,(-ONE),Fcn0,1,dFdT,1)
   CALL cb05_sorg_vbs_aq_WSCAL(NVAR,(ONE/Delta),dFdT,1)

  END SUBROUTINE  cb05_sorg_vbs_aq_ros_FunTimeDeriv



  SUBROUTINE  cb05_sorg_vbs_aq_ros_PrepareMatrix ( H, Direction, gam, &
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


     CALL cb05_sorg_vbs_aq_WCOPY(LU_NONZERO,Jac0,1,Ghimj,1)
     CALL cb05_sorg_vbs_aq_WSCAL(LU_NONZERO,(-ONE),Ghimj,1)
     ghinv = ONE/(Direction*H*gam)
     DO i=1,NVAR
       Ghimj(LU_DIAG(i)) = Ghimj(LU_DIAG(i))+ghinv
     END DO

     CALL cb05_sorg_vbs_aq_ros_Decomp( Ghimj, Pivot, ising, Ndec )
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

  END SUBROUTINE  cb05_sorg_vbs_aq_ros_PrepareMatrix



  SUBROUTINE  cb05_sorg_vbs_aq_ros_Decomp( A, Pivot, ising, Ndec )



   IMPLICIT NONE

   REAL(kind=dp), INTENT(INOUT) :: A(LU_NONZERO)

   INTEGER, INTENT(OUT) :: Pivot(NVAR), ising
   INTEGER, INTENT(INOUT) :: Ndec 

   CALL cb05_sorg_vbs_aq_KppDecomp ( A, ising )
   Pivot(1) = 1
   Ndec = Ndec + 1

  END SUBROUTINE  cb05_sorg_vbs_aq_ros_Decomp



  SUBROUTINE  cb05_sorg_vbs_aq_ros_Solve( A, Pivot, b, Nsol )



   IMPLICIT NONE

   REAL(kind=dp), INTENT(IN) :: A(LU_NONZERO)
   INTEGER, INTENT(IN) :: Pivot(NVAR)

   INTEGER, INTENT(INOUT) :: nsol 

   REAL(kind=dp), INTENT(INOUT) :: b(NVAR)


   CALL cb05_sorg_vbs_aq_KppSolve( A, b )

   Nsol = Nsol+1

  END SUBROUTINE  cb05_sorg_vbs_aq_ros_Solve




  SUBROUTINE  cb05_sorg_vbs_aq_Ros2 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

 END SUBROUTINE  cb05_sorg_vbs_aq_Ros2



  SUBROUTINE  cb05_sorg_vbs_aq_Ros3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  cb05_sorg_vbs_aq_Ros3





  SUBROUTINE  cb05_sorg_vbs_aq_Ros4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  cb05_sorg_vbs_aq_Ros4


  SUBROUTINE  cb05_sorg_vbs_aq_Rodas3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  cb05_sorg_vbs_aq_Rodas3


  SUBROUTINE  cb05_sorg_vbs_aq_Rodas4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  cb05_sorg_vbs_aq_Rodas4




END SUBROUTINE  cb05_sorg_vbs_aq_Rosenbrock




SUBROUTINE  cb05_sorg_vbs_aq_FunTemplate( T, Y, Ydot, RCONST, FIX, Nfun )




   USE cb05_sorg_vbs_aq_Parameters




   REAL(kind=dp) :: T, Y(NVAR)
   REAL(kind=dp) :: RCONST(NREACT)
   REAL(kind=dp) :: FIX(NFIX)

   REAL(kind=dp) :: Ydot(NVAR)
   INTEGER :: Nfun









   CALL cb05_sorg_vbs_aq_Fun( Y, FIX, RCONST, Ydot )


   Nfun = Nfun+1

END SUBROUTINE  cb05_sorg_vbs_aq_FunTemplate



SUBROUTINE  cb05_sorg_vbs_aq_JacTemplate( T, Y, Jcb, FIX, Njac, RCONST )




 USE cb05_sorg_vbs_aq_Parameters
 
 USE cb05_sorg_vbs_aq_Jacobian



    REAL(kind=dp) :: T, Y(NVAR)
    REAL(kind=dp) :: FIX(NFIX)
    REAL(kind=dp) :: RCONST(NREACT)

    INTEGER :: Njac


    REAL(kind=dp) :: Jcb(LU_NONZERO)

    REAL(kind=dp) :: Told





    CALL cb05_sorg_vbs_aq_Jac_SP( Y, FIX, RCONST, Jcb )


    Njac = Njac+1

END SUBROUTINE  cb05_sorg_vbs_aq_JacTemplate

















SUBROUTINE cb05_sorg_vbs_aq_Fun ( V, F, RCT, Vdot )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: Vdot(NVAR)




  REAL(kind=dp) :: A(NREACT)


  A(1) = RCT(1)*V(105)
  A(2) = RCT(2)*V(101)
  A(3) = RCT(3)*V(94)*V(99)
  A(4) = RCT(4)*V(101)*V(105)
  A(5) = RCT(5)*V(101)*V(105)
  A(6) = RCT(6)*V(99)*V(101)
  A(7) = RCT(7)*V(94)*V(105)
  A(8) = RCT(8)*V(94)
  A(9) = RCT(9)*V(94)
  A(10) = RCT(10)*V(53)*F(2)
  A(11) = RCT(11)*V(53)*F(1)
  A(12) = RCT(12)*V(94)*V(98)
  A(13) = RCT(13)*V(94)*V(102)
  A(14) = RCT(14)*V(95)
  A(15) = RCT(15)*V(95)
  A(16) = RCT(16)*V(95)*V(99)
  A(17) = RCT(17)*V(95)*V(105)
  A(18) = RCT(18)*V(95)*V(105)
  A(19) = RCT(19)*V(54)*F(1)
  A(20) = RCT(20)*V(54)*F(1)*F(1)
  A(21) = RCT(21)*V(54)
  A(22) = RCT(22)*V(99)*V(99)
  A(23) = RCT(23)*V(99)*V(105)*F(1)
  A(24) = RCT(24)*V(98)*V(99)
  A(25) = RCT(25)*V(58)
  A(26) = RCT(26)*V(58)*V(98)
  A(27) = RCT(27)*V(58)*V(58)
  A(28) = RCT(28)*V(98)*V(105)
  A(29) = RCT(29)*V(77)*V(98)
  A(30) = RCT(30)*V(99)*V(102)
  A(31) = RCT(31)*V(102)*V(105)
  A(32) = RCT(32)*V(64)
  A(33) = RCT(33)*V(64)*V(98)
  A(34) = RCT(34)*V(102)*V(102)
  A(35) = RCT(35)*V(102)*V(102)*F(1)
  A(36) = RCT(36)*V(68)
  A(37) = RCT(37)*V(68)*V(98)
  A(38) = RCT(38)*V(52)*V(53)
  A(39) = RCT(39)*V(52)*V(98)
  A(40) = RCT(40)*V(98)*V(101)
  A(41) = RCT(41)*V(98)*V(98)
  A(42) = RCT(42)*V(98)*V(98)
  A(43) = RCT(43)*V(98)*V(102)
  A(44) = RCT(44)*V(101)*V(102)
  A(45) = RCT(45)*V(68)*V(101)
  A(46) = RCT(46)*V(95)*V(101)
  A(47) = RCT(47)*V(95)*V(98)
  A(48) = RCT(48)*V(95)*V(102)
  A(49) = RCT(49)*V(94)*V(95)
  A(50) = RCT(50)*V(95)*V(95)
  A(51) = RCT(51)*V(64)
  A(52) = RCT(52)*V(77)
  A(53) = RCT(53)*V(54)
  A(54) = RCT(54)*V(96)*V(99)
  A(55) = RCT(55)*V(86)*V(99)
  A(56) = RCT(56)*V(96)*V(102)
  A(57) = RCT(57)*V(86)*V(102)
  A(58) = RCT(58)*V(96)*V(96)
  A(59) = RCT(59)*V(86)*V(86)
  A(60) = RCT(60)*V(86)*V(96)
  A(61) = RCT(61)*V(93)*V(98)
  A(62) = RCT(62)*V(93)
  A(63) = RCT(63)*V(72)*V(98)
  A(64) = RCT(64)*V(72)
  A(65) = RCT(65)*V(78)*V(98)
  A(66) = RCT(66)*V(56)*V(98)
  A(67) = RCT(67)*V(99)*V(104)
  A(68) = RCT(68)*V(102)*V(104)
  A(69) = RCT(69)*V(104)*V(104)
  A(70) = RCT(70)*V(76)*V(98)
  A(71) = RCT(71)*V(76)
  A(72) = RCT(72)*V(66)*V(98)
  A(73) = RCT(73)*V(92)*V(98)
  A(74) = RCT(74)*V(92)
  A(75) = RCT(75)*V(92)
  A(76) = RCT(76)*V(92)*V(101)
  A(77) = RCT(77)*V(92)*V(95)
  A(78) = RCT(78)*V(92)*V(102)
  A(79) = RCT(79)*V(69)
  A(80) = RCT(80)*V(69)*V(99)
  A(81) = RCT(81)*V(69)*V(102)
  A(82) = RCT(82)*V(59)*V(98)
  A(83) = RCT(83)*V(91)*V(101)
  A(84) = RCT(84)*V(91)*V(98)
  A(85) = RCT(85)*V(91)*V(95)
  A(86) = RCT(86)*V(91)
  A(87) = RCT(87)*V(97)*V(99)
  A(88) = RCT(88)*V(97)*V(105)
  A(89) = RCT(89)*V(49)
  A(90) = RCT(90)*V(49)
  A(91) = RCT(91)*V(97)*V(102)
  A(92) = RCT(92)*V(97)*V(104)
  A(93) = RCT(93)*V(96)*V(97)
  A(94) = RCT(94)*V(97)*V(97)
  A(95) = RCT(95)*V(61)*V(98)
  A(96) = RCT(96)*V(61)
  A(97) = RCT(97)*V(62)*V(98)
  A(98) = RCT(98)*V(90)*V(101)
  A(99) = RCT(99)*V(90)*V(98)
  A(100) = RCT(100)*V(90)*V(95)
  A(101) = RCT(101)*V(90)
  A(102) = RCT(102)*V(99)*V(103)
  A(103) = RCT(103)*V(103)*V(105)
  A(104) = RCT(104)*V(63)
  A(105) = RCT(105)*V(63)
  A(106) = RCT(106)*V(63)*V(98)
  A(107) = RCT(107)*V(102)*V(103)
  A(108) = RCT(108)*V(103)*V(104)
  A(109) = RCT(109)*V(96)*V(103)
  A(110) = RCT(110)*V(103)*V(103)
  A(111) = RCT(111)*V(97)*V(103)
  A(112) = RCT(112)*V(87)*V(98)
  A(113) = RCT(113)*V(80)
  A(114) = 1600*V(80)
  A(115) = RCT(115)*V(80)*V(105)
  A(116) = RCT(116)*V(85)*V(101)
  A(117) = RCT(117)*V(85)*V(98)
  A(118) = RCT(118)*V(85)*V(94)
  A(119) = RCT(119)*V(85)*V(95)
  A(120) = RCT(120)*V(81)*V(101)
  A(121) = RCT(121)*V(81)*V(98)
  A(122) = RCT(122)*V(81)*V(94)
  A(123) = RCT(123)*V(81)*V(95)
  A(124) = RCT(124)*V(84)*V(101)
  A(125) = RCT(125)*V(84)*V(98)
  A(126) = RCT(126)*V(84)*V(94)
  A(127) = RCT(127)*V(84)*V(95)
  A(128) = RCT(128)*V(51)*V(98)
  A(129) = RCT(129)*V(60)*V(99)
  A(130) = 4.2*V(60)
  A(131) = RCT(131)*V(82)*V(98)
  A(132) = RCT(132)*V(82)*V(95)
  A(133) = RCT(133)*V(71)*V(105)
  A(134) = RCT(134)*V(71)*V(102)
  A(135) = RCT(135)*V(79)
  A(136) = RCT(136)*V(79)*V(98)
  A(137) = RCT(137)*V(79)*V(94)
  A(138) = RCT(138)*V(55)*V(98)
  A(139) = RCT(139)*V(74)*V(98)
  A(140) = RCT(140)*V(74)
  A(141) = RCT(141)*V(88)*V(101)
  A(142) = RCT(142)*V(88)*V(98)
  A(143) = RCT(143)*V(88)*V(94)
  A(144) = RCT(144)*V(88)*V(95)
  A(145) = RCT(145)*V(89)*V(98)
  A(146) = RCT(146)*V(89)*V(94)
  A(147) = RCT(147)*V(89)*V(95)
  A(148) = RCT(148)*V(89)
  A(149) = RCT(149)*V(83)*V(101)
  A(150) = RCT(150)*V(83)*V(98)
  A(151) = RCT(151)*V(83)*V(94)
  A(152) = RCT(152)*V(83)*V(95)
  A(153) = RCT(153)*V(48)*V(98)
  A(154) = RCT(154)*V(67)*V(98)
  A(155) = RCT(155)*V(65)*V(98)
  A(156) = RCT(156)*V(88)*V(105)
  A(157) = RCT(157)*V(47)
  A(158) = RCT(158)*V(50)
  A(159) = RCT(159)*V(94)*V(100)
  A(160) = RCT(160)*V(73)*V(73)
  A(161) = RCT(161)*V(73)*V(99)
  A(162) = RCT(162)*V(73)*V(102)
  A(163) = RCT(163)*V(75)*V(98)
  A(164) = RCT(164)*V(75)
  A(165) = RCT(165)*V(56)*V(100)
  A(166) = RCT(166)*V(87)*V(100)
  A(167) = RCT(167)*V(65)*V(100)
  A(168) = RCT(168)*V(81)*V(100)
  A(169) = RCT(169)*V(85)*V(100)
  A(170) = RCT(170)*V(84)*V(100)
  A(171) = RCT(171)*V(88)*V(100)
  A(172) = RCT(172)*V(92)*V(100)
  A(173) = RCT(173)*V(91)*V(100)
  A(174) = RCT(174)*V(90)*V(100)
  A(175) = RCT(175)*V(66)*V(100)
  A(176) = RCT(176)*V(67)*V(100)
  A(177) = RCT(177)*V(70)*V(98)
  A(178) = RCT(178)*V(57)*V(94)
  A(179) = RCT(179)*V(57)*V(98)
  A(180) = RCT(180)*V(57)*V(68)
  A(181) = RCT(181)*V(12)*V(98)
  A(182) = RCT(182)*V(15)*V(98)
  A(183) = RCT(183)*V(18)*V(98)
  A(184) = RCT(184)*V(23)*V(98)
  A(185) = RCT(185)*V(23)*V(94)
  A(186) = RCT(186)*V(29)*V(98)
  A(187) = RCT(187)*V(29)*V(94)
  A(188) = RCT(188)*V(29)*V(95)
  A(189) = RCT(189)*V(32)*V(98)
  A(190) = RCT(190)*V(34)*V(98)
  A(191) = RCT(191)*V(37)*V(98)
  A(192) = RCT(192)*V(42)*V(98)
  A(193) = RCT(193)*V(41)*V(98)
  A(194) = RCT(194)*V(40)*V(98)
  A(195) = RCT(195)*V(46)*V(98)
  A(196) = RCT(196)*V(45)*V(98)
  A(197) = RCT(197)*V(44)*V(98)
  A(198) = RCT(198)*V(48)


  Vdot(1) = 0.071*A(128)
  Vdot(2) = 0.138*A(128)
  Vdot(3) = A(131)+A(132)
  Vdot(4) = 0.038*A(138)
  Vdot(5) = 0.167*A(138)
  Vdot(6) = 0.232*A(141)+0.232*A(142)+0.232*A(143)+0.232*A(144)
  Vdot(7) = 0.0228*A(141)+0.0288*A(142)+0.0288*A(143)+0.0288*A(144)
  Vdot(8) = A(153)+A(198)
  Vdot(9) = A(153)
  Vdot(10) = A(149)+A(150)+A(151)+A(152)
  Vdot(11) = A(181)
  Vdot(12) = -A(181)
  Vdot(13) = 0.239*A(182)
  Vdot(14) = 0.363*A(182)
  Vdot(15) = -A(182)
  Vdot(16) = 0.045*A(183)
  Vdot(17) = 0.149*A(183)
  Vdot(18) = -A(183)
  Vdot(19) = 0.038*A(184)
  Vdot(20) = 0.326*A(184)
  Vdot(21) = 0.125*A(185)
  Vdot(22) = 0.102*A(185)
  Vdot(23) = -A(184)-A(185)
  Vdot(24) = 0.13*A(186)
  Vdot(25) = 0.0406*A(186)
  Vdot(26) = 0.026*A(187)
  Vdot(27) = 0.485*A(187)
  Vdot(28) = A(188)
  Vdot(29) = -A(186)-A(187)-A(188)
  Vdot(30) = 0.091*A(189)
  Vdot(31) = 0.367*A(189)
  Vdot(32) = -A(189)
  Vdot(33) = 1.173*A(190)
  Vdot(34) = -A(190)
  Vdot(35) = 0.156*A(191)
  Vdot(36) = 0.777*A(191)
  Vdot(37) = -A(191)
  Vdot(38) = A(178)+A(179)+A(180)
  Vdot(39) = 1.075*A(194)
  Vdot(40) = 1.075*A(193)-A(194)
  Vdot(41) = 1.075*A(192)-A(193)
  Vdot(42) = -A(192)
  Vdot(43) = 1.075*A(197)
  Vdot(44) = 1.075*A(196)-A(197)
  Vdot(45) = 1.075*A(195)-A(196)
  Vdot(46) = -A(195)
  Vdot(47) = -A(157)+0.3*A(160)
  Vdot(48) = -A(153)-A(198)
  Vdot(49) = A(88)-A(89)-A(90)
  Vdot(50) = -A(158)+A(162)
  Vdot(51) = -A(128)
  Vdot(52) = -A(38)-A(39)
  Vdot(53) = A(9)-A(10)-A(11)-A(38)
  Vdot(54) = A(18)-A(19)-A(20)-A(21)-A(53)
  Vdot(55) = -A(138)
  Vdot(56) = -A(66)-A(165)
  Vdot(57) = -A(178)-A(179)-A(180)
  Vdot(58) = 2*A(23)+A(24)-A(25)-A(26)-2*A(27)
  Vdot(59) = A(80)-A(82)+0.37*A(122)
  Vdot(60) = 0.56*A(128)-A(129)-A(130)+0.3*A(138)
  Vdot(61) = 0.8*A(91)-A(95)-A(96)+0.8*A(107)
  Vdot(62) = 0.2*A(91)+0.1*A(92)+0.1*A(93)-A(97)+0.2*A(107)+0.1*A(108)+0.1*A(109)
  Vdot(63) = A(103)-A(104)-A(105)-A(106)
  Vdot(64) = A(31)-A(32)-A(33)-A(51)
  Vdot(65) = -A(155)-A(167)
  Vdot(66) = 0.63*A(69)-A(72)-A(175)
  Vdot(67) = -A(154)-A(176)
  Vdot(68) = A(34)+A(35)-A(36)-A(37)+A(42)-A(45)-A(180)
  Vdot(69) = A(78)-A(79)-A(80)-A(81)
  Vdot(70) = A(165)+A(166)+A(167)+0.3*A(170)+0.15*A(171)+A(172)+A(173)+A(174)+A(175)+A(176)-A(177)
  Vdot(71) = 0.4*A(131)+A(132)-A(133)-A(134)
  Vdot(72) = A(56)+A(57)-A(63)-A(64)
  Vdot(73) = A(159)-2*A(160)-A(161)-A(162)
  Vdot(74) = 0.2*A(137)+0.8*A(138)-A(139)-A(140)+0.168*A(145)+0.85*A(146)
  Vdot(75) = -A(163)-A(164)+A(168)+A(169)+0.7*A(170)+0.85*A(171)
  Vdot(76) = A(68)-A(70)-A(71)+A(81)
  Vdot(77) = 2*A(19)+2*A(20)+A(28)-A(29)+A(48)-A(52)+A(61)+A(77)+A(85)+A(100)+A(132)+0.15*A(147)
  Vdot(78) = -A(65)+A(73)+A(74)+A(75)+A(76)+A(77)+A(86)+A(101)+0.2*A(116)+0.33*A(118)+A(120)+0.63*A(122)+0.1*A(124)+0.25&
               &*A(126)+A(135)+2*A(136)+0.69*A(137)+A(140)+0.066*A(143)+0.334*A(145)+0.225*A(146)+0.643*A(147)+0.333*A(148)&
               &+0.001*A(151)+A(163)+A(164)+A(172)
  Vdot(79) = 0.9*A(129)+0.3*A(131)-A(135)-A(136)-A(137)
  Vdot(80) = 0.76*A(112)-0.98*A(113)-A(114)-A(115)+0.76*A(166)
  Vdot(81) = -A(120)-A(121)-A(122)-A(123)-A(168)
  Vdot(82) = 0.36*A(128)+A(130)-A(131)-A(132)+A(134)+0.2*A(138)
  Vdot(83) = -A(149)-A(150)-A(151)-A(152)
  Vdot(84) = -A(124)-A(125)-A(126)-A(127)-A(170)
  Vdot(85) = -A(116)-A(117)-A(118)-A(119)-A(169)+0.3*A(170)
  Vdot(86) = -A(55)-A(57)-2*A(59)-A(60)+0.13*A(112)+0.04*A(113)+0.01*A(116)+0.09*A(119)+0.088*A(142)+0.25*A(150)+0.18&
               &*A(151)+0.25*A(152)+0.009*A(155)+0.13*A(166)+0.009*A(167)
  Vdot(87) = -0.66*A(61)-0.66*A(62)-1.11*A(112)-2.1*A(113)+0.2*A(116)-0.7*A(117)-A(118)-A(119)+0.1*A(124)+1.1*A(138)&
               &+0.35*A(143)+1.565*A(145)+0.36*A(146)+1.282*A(147)+0.832*A(148)+5.12*A(149)+1.66*A(150)+7*A(151)+2.4*A(156)&
               &-1.11*A(166)-A(169)+0.3*A(170)
  Vdot(88) = -A(141)-A(142)-A(143)-A(144)-A(156)-A(171)
  Vdot(89) = 0.75*A(141)+0.912*A(142)+0.65*A(143)+0.2*A(144)-A(145)-A(146)-A(147)-A(148)+0.2*A(156)+A(171)
  Vdot(90) = 0.33*A(61)+0.33*A(62)+0.5*A(63)+0.5*A(64)-A(98)-A(99)-A(100)-A(101)+0.05*A(112)+0.5*A(113)+0.3*A(116)+0.62&
               &*A(117)+0.32*A(118)+0.56*A(119)+0.22*A(121)+0.66*A(124)+0.7*A(125)+0.35*A(126)+0.64*A(127)+0.03*A(137)+0.15&
               &*A(143)+0.8*A(144)+0.12*A(145)+0.357*A(147)+0.47*A(150)+0.21*A(151)+0.47*A(152)+0.05*A(154)+0.8*A(156)+0.05&
               &*A(166)+0.67*A(169)+0.55*A(170)-A(174)
  Vdot(91) = 0.33*A(61)+0.33*A(62)+0.5*A(63)+0.5*A(64)-A(83)-A(84)-A(85)-A(86)+A(102)+A(106)+0.9*A(108)+0.9*A(109)+2&
               &*A(110)+A(111)+0.06*A(112)+0.6*A(113)+0.2*A(116)+0.33*A(117)+0.18*A(118)+0.35*A(119)+1.24*A(124)+1.3*A(125)&
               &+0.65*A(126)+1.18*A(127)+0.252*A(145)+0.02*A(146)+0.067*A(148)+0.9*A(154)+0.991*A(155)+0.06*A(166)+0.991&
               &*A(167)+0.33*A(169)+0.45*A(170)-A(173)+A(176)
  Vdot(92) = 0.33*A(61)+0.33*A(62)+A(67)+1.37*A(69)+A(71)+A(72)-A(73)-A(74)-A(75)-A(76)-A(77)-A(78)+A(79)+A(92)+0.1&
               &*A(108)+0.2*A(116)+0.8*A(117)+0.74*A(118)+A(119)+A(120)+1.56*A(121)+A(122)+2*A(123)+0.25*A(126)+A(136)+0.7&
               &*A(137)+0.5*A(141)+0.629*A(142)+0.6*A(143)+0.167*A(145)+0.15*A(146)+0.282*A(147)+0.9*A(148)+0.28*A(150)+0.24&
               &*A(151)+0.1*A(154)+A(168)-A(172)+A(175)
  Vdot(93) = A(55)-A(61)-A(62)+A(115)+0.1*A(129)+A(133)+0.8*A(144)+0.85*A(147)+0.53*A(152)+0.8*A(156)
  Vdot(94) = A(2)-A(3)-A(7)-A(8)-A(9)-A(12)-A(13)-A(49)+0.2*A(91)+0.2*A(107)-A(118)-A(122)-A(126)-A(137)-A(143)-A(146)&
               &-A(151)-A(159)-A(178)
  Vdot(95) = A(5)+A(7)-A(14)-A(15)-A(16)-A(17)-A(18)+A(21)+A(29)-A(46)-A(47)-A(48)-A(49)-2*A(50)+0.39*A(51)+A(53)-A(77)&
               &-A(85)-A(100)-A(119)-A(123)-A(127)-A(132)-A(144)-A(147)-A(152)
  Vdot(96) = -A(54)-A(56)-2*A(58)-A(60)+A(63)+0.3*A(70)-A(93)+A(102)+0.9*A(108)-A(109)+2*A(110)+A(111)+0.87*A(112)+0.96&
               &*A(113)+0.2*A(116)+0.8*A(117)+0.22*A(118)+0.91*A(119)+0.7*A(120)+A(121)+A(123)+0.1*A(124)+A(125)+0.08*A(128)&
               &+0.6*A(131)+A(136)+0.03*A(137)+0.5*A(138)+A(139)+0.25*A(141)+0.991*A(142)+0.2*A(143)+A(144)+0.713*A(145)&
               &+0.064*A(146)+0.075*A(147)+0.7*A(148)+1.25*A(150)+0.76*A(151)+1.03*A(152)+0.1*A(154)+0.991*A(155)+A(156)&
               &+0.87*A(166)+0.991*A(167)+2*A(168)+2*A(169)+1.7*A(170)+A(171)
  Vdot(97) = A(83)+A(84)+A(85)-A(87)-A(88)+A(89)+A(90)-A(91)-A(92)-A(93)-2*A(94)+A(95)-A(111)+A(135)+A(136)+0.62*A(137)&
               &+A(139)+A(140)+0.21*A(145)+0.114*A(146)+0.967*A(148)+A(173)
  Vdot(98) = 2*A(11)-A(12)+A(13)-A(24)+A(25)-A(26)-A(28)-A(29)+A(30)-A(33)+2*A(36)-A(37)+A(38)-A(39)-A(40)-2*A(41)-2&
               &*A(42)-A(43)+A(44)+A(45)-A(47)+0.39*A(51)+A(52)-A(61)-A(63)+A(64)-A(65)-A(66)-A(70)+A(71)-A(72)-A(73)+A(76)&
               &-A(82)+A(83)-A(84)-A(95)+A(96)-A(97)+A(98)-A(99)-A(106)-A(112)+0.1*A(116)-A(117)+0.1*A(118)+0.3*A(120)&
               &-A(121)+0.13*A(122)-A(125)+0.5*A(126)-A(128)-A(131)-A(136)+0.08*A(137)-A(138)-A(139)-A(142)+0.266*A(143)&
               &-A(145)+0.268*A(146)+0.15*A(149)-A(150)+0.57*A(151)-A(153)-A(154)-A(155)+A(158)-A(163)-A(177)-A(179)
  Vdot(99) = A(1)-A(3)+A(4)-A(6)+A(15)-A(16)+A(17)-2*A(22)-A(23)-A(24)+A(25)+A(27)-A(30)-A(54)-A(55)-A(67)-A(80)-A(87)&
               &-A(102)-A(129)+0.2*A(156)-A(161)
  Vdot(100) = 2*A(157)+A(158)-A(159)+1.4*A(160)+A(161)+A(163)+A(164)-A(165)-A(166)-A(167)-A(168)-A(169)-A(170)-A(171)&
                &-A(172)-A(173)-A(174)-A(175)-A(176)+A(177)
  Vdot(101) = A(1)-A(2)-A(4)-A(5)-A(6)+A(8)+A(10)+A(14)-A(40)+A(41)-A(44)-A(45)-A(46)-A(76)-A(83)-A(98)-A(116)-A(120)&
                &-A(124)+0.5*A(126)-A(141)-A(149)
  Vdot(102) = A(12)-A(13)-A(30)-A(31)+A(32)-2*A(34)-2*A(35)+A(37)+A(38)+A(39)+A(40)-A(43)-A(44)+A(45)+A(47)-A(48)+0.61&
                &*A(51)-A(56)-A(57)+A(61)+A(62)+A(64)+A(65)+A(67)-A(68)+0.74*A(69)+0.3*A(70)+A(71)+A(72)+A(73)+2*A(74)+A(76)&
                &+A(77)-A(78)+A(79)+A(80)-A(81)+A(82)+A(86)-A(91)+0.9*A(92)+A(101)+A(102)-A(107)+A(108)+2*A(110)+A(111)+0.11&
                &*A(112)+0.94*A(113)+A(114)+0.3*A(116)+0.95*A(117)+0.44*A(118)+1.7*A(120)+A(121)+0.13*A(122)+0.1*A(124)&
                &+A(125)+0.5*A(126)+A(127)+0.44*A(128)+0.9*A(129)+A(130)+0.6*A(131)-A(134)+A(135)+2*A(136)+0.76*A(137)+0.7&
                &*A(138)+A(140)+0.25*A(141)+0.912*A(142)+0.066*A(143)+0.8*A(144)+0.503*A(145)+0.154*A(146)+0.925*A(147)&
                &+1.033*A(148)+0.75*A(150)+0.07*A(151)+0.28*A(152)+A(153)+A(154)+A(155)+0.8*A(156)-A(162)+A(164)+0.11*A(166)&
                &+A(167)+A(168)+A(169)+A(170)+A(171)+A(172)+A(175)+A(176)
  Vdot(103) = A(98)+A(99)+A(100)-A(102)-A(103)+A(104)+A(105)-A(107)-A(108)-A(109)-2*A(110)-A(111)+0.25*A(141)+0.2*A(143)&
                &+0.25*A(145)+0.075*A(147)+0.39*A(151)+A(174)
  Vdot(104) = A(66)-A(67)-A(68)-2*A(69)+0.7*A(70)+A(86)+A(87)-0.1*A(92)+0.9*A(93)+2*A(94)+A(96)+A(97)+A(101)-A(108)&
                &+A(111)+A(165)
  Vdot(105) = -A(1)+A(3)-A(4)-A(5)+A(6)-A(7)+A(14)+2*A(16)-A(18)+A(21)+2*A(22)-A(23)+A(26)+A(27)-A(28)+A(30)-A(31)+A(32)&
                &+A(33)+A(46)+A(47)+A(49)+2*A(50)+0.61*A(51)+A(52)+A(53)+A(54)+A(62)+A(67)+A(80)+A(87)-A(88)+A(89)+A(90)&
                &+A(102)-A(103)+A(104)+A(105)+A(106)-A(115)+A(119)+A(123)+A(127)+0.9*A(129)-A(133)+0.2*A(144)+0.47*A(152)&
                &-A(156)+A(161)
      
END SUBROUTINE cb05_sorg_vbs_aq_Fun
















SUBROUTINE cb05_sorg_vbs_aq_IRRFun ( V, F, RCT, IRR )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: IRR(NREACT)



  IRR(1) = RCT(1)*V(105)
  IRR(2) = RCT(2)*V(101)
  IRR(3) = RCT(3)*V(94)*V(99)
  IRR(4) = RCT(4)*V(101)*V(105)
  IRR(5) = RCT(5)*V(101)*V(105)
  IRR(6) = RCT(6)*V(99)*V(101)
  IRR(7) = RCT(7)*V(94)*V(105)
  IRR(8) = RCT(8)*V(94)
  IRR(9) = RCT(9)*V(94)
  IRR(10) = RCT(10)*V(53)*F(2)
  IRR(11) = RCT(11)*V(53)*F(1)
  IRR(12) = RCT(12)*V(94)*V(98)
  IRR(13) = RCT(13)*V(94)*V(102)
  IRR(14) = RCT(14)*V(95)
  IRR(15) = RCT(15)*V(95)
  IRR(16) = RCT(16)*V(95)*V(99)
  IRR(17) = RCT(17)*V(95)*V(105)
  IRR(18) = RCT(18)*V(95)*V(105)
  IRR(19) = RCT(19)*V(54)*F(1)
  IRR(20) = RCT(20)*V(54)*F(1)*F(1)
  IRR(21) = RCT(21)*V(54)
  IRR(22) = RCT(22)*V(99)*V(99)
  IRR(23) = RCT(23)*V(99)*V(105)*F(1)
  IRR(24) = RCT(24)*V(98)*V(99)
  IRR(25) = RCT(25)*V(58)
  IRR(26) = RCT(26)*V(58)*V(98)
  IRR(27) = RCT(27)*V(58)*V(58)
  IRR(28) = RCT(28)*V(98)*V(105)
  IRR(29) = RCT(29)*V(77)*V(98)
  IRR(30) = RCT(30)*V(99)*V(102)
  IRR(31) = RCT(31)*V(102)*V(105)
  IRR(32) = RCT(32)*V(64)
  IRR(33) = RCT(33)*V(64)*V(98)
  IRR(34) = RCT(34)*V(102)*V(102)
  IRR(35) = RCT(35)*V(102)*V(102)*F(1)
  IRR(36) = RCT(36)*V(68)
  IRR(37) = RCT(37)*V(68)*V(98)
  IRR(38) = RCT(38)*V(52)*V(53)
  IRR(39) = RCT(39)*V(52)*V(98)
  IRR(40) = RCT(40)*V(98)*V(101)
  IRR(41) = RCT(41)*V(98)*V(98)
  IRR(42) = RCT(42)*V(98)*V(98)
  IRR(43) = RCT(43)*V(98)*V(102)
  IRR(44) = RCT(44)*V(101)*V(102)
  IRR(45) = RCT(45)*V(68)*V(101)
  IRR(46) = RCT(46)*V(95)*V(101)
  IRR(47) = RCT(47)*V(95)*V(98)
  IRR(48) = RCT(48)*V(95)*V(102)
  IRR(49) = RCT(49)*V(94)*V(95)
  IRR(50) = RCT(50)*V(95)*V(95)
  IRR(51) = RCT(51)*V(64)
  IRR(52) = RCT(52)*V(77)
  IRR(53) = RCT(53)*V(54)
  IRR(54) = RCT(54)*V(96)*V(99)
  IRR(55) = RCT(55)*V(86)*V(99)
  IRR(56) = RCT(56)*V(96)*V(102)
  IRR(57) = RCT(57)*V(86)*V(102)
  IRR(58) = RCT(58)*V(96)*V(96)
  IRR(59) = RCT(59)*V(86)*V(86)
  IRR(60) = RCT(60)*V(86)*V(96)
  IRR(61) = RCT(61)*V(93)*V(98)
  IRR(62) = RCT(62)*V(93)
  IRR(63) = RCT(63)*V(72)*V(98)
  IRR(64) = RCT(64)*V(72)
  IRR(65) = RCT(65)*V(78)*V(98)
  IRR(66) = RCT(66)*V(56)*V(98)
  IRR(67) = RCT(67)*V(99)*V(104)
  IRR(68) = RCT(68)*V(102)*V(104)
  IRR(69) = RCT(69)*V(104)*V(104)
  IRR(70) = RCT(70)*V(76)*V(98)
  IRR(71) = RCT(71)*V(76)
  IRR(72) = RCT(72)*V(66)*V(98)
  IRR(73) = RCT(73)*V(92)*V(98)
  IRR(74) = RCT(74)*V(92)
  IRR(75) = RCT(75)*V(92)
  IRR(76) = RCT(76)*V(92)*V(101)
  IRR(77) = RCT(77)*V(92)*V(95)
  IRR(78) = RCT(78)*V(92)*V(102)
  IRR(79) = RCT(79)*V(69)
  IRR(80) = RCT(80)*V(69)*V(99)
  IRR(81) = RCT(81)*V(69)*V(102)
  IRR(82) = RCT(82)*V(59)*V(98)
  IRR(83) = RCT(83)*V(91)*V(101)
  IRR(84) = RCT(84)*V(91)*V(98)
  IRR(85) = RCT(85)*V(91)*V(95)
  IRR(86) = RCT(86)*V(91)
  IRR(87) = RCT(87)*V(97)*V(99)
  IRR(88) = RCT(88)*V(97)*V(105)
  IRR(89) = RCT(89)*V(49)
  IRR(90) = RCT(90)*V(49)
  IRR(91) = RCT(91)*V(97)*V(102)
  IRR(92) = RCT(92)*V(97)*V(104)
  IRR(93) = RCT(93)*V(96)*V(97)
  IRR(94) = RCT(94)*V(97)*V(97)
  IRR(95) = RCT(95)*V(61)*V(98)
  IRR(96) = RCT(96)*V(61)
  IRR(97) = RCT(97)*V(62)*V(98)
  IRR(98) = RCT(98)*V(90)*V(101)
  IRR(99) = RCT(99)*V(90)*V(98)
  IRR(100) = RCT(100)*V(90)*V(95)
  IRR(101) = RCT(101)*V(90)
  IRR(102) = RCT(102)*V(99)*V(103)
  IRR(103) = RCT(103)*V(103)*V(105)
  IRR(104) = RCT(104)*V(63)
  IRR(105) = RCT(105)*V(63)
  IRR(106) = RCT(106)*V(63)*V(98)
  IRR(107) = RCT(107)*V(102)*V(103)
  IRR(108) = RCT(108)*V(103)*V(104)
  IRR(109) = RCT(109)*V(96)*V(103)
  IRR(110) = RCT(110)*V(103)*V(103)
  IRR(111) = RCT(111)*V(97)*V(103)
  IRR(112) = RCT(112)*V(87)*V(98)
  IRR(113) = RCT(113)*V(80)
  IRR(114) = 1600*V(80)
  IRR(115) = RCT(115)*V(80)*V(105)
  IRR(116) = RCT(116)*V(85)*V(101)
  IRR(117) = RCT(117)*V(85)*V(98)
  IRR(118) = RCT(118)*V(85)*V(94)
  IRR(119) = RCT(119)*V(85)*V(95)
  IRR(120) = RCT(120)*V(81)*V(101)
  IRR(121) = RCT(121)*V(81)*V(98)
  IRR(122) = RCT(122)*V(81)*V(94)
  IRR(123) = RCT(123)*V(81)*V(95)
  IRR(124) = RCT(124)*V(84)*V(101)
  IRR(125) = RCT(125)*V(84)*V(98)
  IRR(126) = RCT(126)*V(84)*V(94)
  IRR(127) = RCT(127)*V(84)*V(95)
  IRR(128) = RCT(128)*V(51)*V(98)
  IRR(129) = RCT(129)*V(60)*V(99)
  IRR(130) = 4.2*V(60)
  IRR(131) = RCT(131)*V(82)*V(98)
  IRR(132) = RCT(132)*V(82)*V(95)
  IRR(133) = RCT(133)*V(71)*V(105)
  IRR(134) = RCT(134)*V(71)*V(102)
  IRR(135) = RCT(135)*V(79)
  IRR(136) = RCT(136)*V(79)*V(98)
  IRR(137) = RCT(137)*V(79)*V(94)
  IRR(138) = RCT(138)*V(55)*V(98)
  IRR(139) = RCT(139)*V(74)*V(98)
  IRR(140) = RCT(140)*V(74)
  IRR(141) = RCT(141)*V(88)*V(101)
  IRR(142) = RCT(142)*V(88)*V(98)
  IRR(143) = RCT(143)*V(88)*V(94)
  IRR(144) = RCT(144)*V(88)*V(95)
  IRR(145) = RCT(145)*V(89)*V(98)
  IRR(146) = RCT(146)*V(89)*V(94)
  IRR(147) = RCT(147)*V(89)*V(95)
  IRR(148) = RCT(148)*V(89)
  IRR(149) = RCT(149)*V(83)*V(101)
  IRR(150) = RCT(150)*V(83)*V(98)
  IRR(151) = RCT(151)*V(83)*V(94)
  IRR(152) = RCT(152)*V(83)*V(95)
  IRR(153) = RCT(153)*V(48)*V(98)
  IRR(154) = RCT(154)*V(67)*V(98)
  IRR(155) = RCT(155)*V(65)*V(98)
  IRR(156) = RCT(156)*V(88)*V(105)
  IRR(157) = RCT(157)*V(47)
  IRR(158) = RCT(158)*V(50)
  IRR(159) = RCT(159)*V(94)*V(100)
  IRR(160) = RCT(160)*V(73)*V(73)
  IRR(161) = RCT(161)*V(73)*V(99)
  IRR(162) = RCT(162)*V(73)*V(102)
  IRR(163) = RCT(163)*V(75)*V(98)
  IRR(164) = RCT(164)*V(75)
  IRR(165) = RCT(165)*V(56)*V(100)
  IRR(166) = RCT(166)*V(87)*V(100)
  IRR(167) = RCT(167)*V(65)*V(100)
  IRR(168) = RCT(168)*V(81)*V(100)
  IRR(169) = RCT(169)*V(85)*V(100)
  IRR(170) = RCT(170)*V(84)*V(100)
  IRR(171) = RCT(171)*V(88)*V(100)
  IRR(172) = RCT(172)*V(92)*V(100)
  IRR(173) = RCT(173)*V(91)*V(100)
  IRR(174) = RCT(174)*V(90)*V(100)
  IRR(175) = RCT(175)*V(66)*V(100)
  IRR(176) = RCT(176)*V(67)*V(100)
  IRR(177) = RCT(177)*V(70)*V(98)
  IRR(178) = RCT(178)*V(57)*V(94)
  IRR(179) = RCT(179)*V(57)*V(98)
  IRR(180) = RCT(180)*V(57)*V(68)
  IRR(181) = RCT(181)*V(12)*V(98)
  IRR(182) = RCT(182)*V(15)*V(98)
  IRR(183) = RCT(183)*V(18)*V(98)
  IRR(184) = RCT(184)*V(23)*V(98)
  IRR(185) = RCT(185)*V(23)*V(94)
  IRR(186) = RCT(186)*V(29)*V(98)
  IRR(187) = RCT(187)*V(29)*V(94)
  IRR(188) = RCT(188)*V(29)*V(95)
  IRR(189) = RCT(189)*V(32)*V(98)
  IRR(190) = RCT(190)*V(34)*V(98)
  IRR(191) = RCT(191)*V(37)*V(98)
  IRR(192) = RCT(192)*V(42)*V(98)
  IRR(193) = RCT(193)*V(41)*V(98)
  IRR(194) = RCT(194)*V(40)*V(98)
  IRR(195) = RCT(195)*V(46)*V(98)
  IRR(196) = RCT(196)*V(45)*V(98)
  IRR(197) = RCT(197)*V(44)*V(98)
  IRR(198) = RCT(198)*V(48)
      
END SUBROUTINE cb05_sorg_vbs_aq_IRRFun
















SUBROUTINE cb05_sorg_vbs_aq_Jac_SP ( V, F, RCT, JVS )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: JVS(LU_NONZERO)




  REAL(kind=dp) :: B(349)


  B(1) = RCT(1)

  B(2) = RCT(2)

  B(3) = RCT(3)*V(99)

  B(4) = RCT(3)*V(94)

  B(5) = RCT(4)*V(105)

  B(6) = RCT(4)*V(101)

  B(7) = RCT(5)*V(105)

  B(8) = RCT(5)*V(101)

  B(9) = RCT(6)*V(101)

  B(10) = RCT(6)*V(99)

  B(11) = RCT(7)*V(105)

  B(12) = RCT(7)*V(94)

  B(13) = RCT(8)

  B(14) = RCT(9)

  B(15) = RCT(10)*F(2)

  B(17) = RCT(11)*F(1)

  B(19) = RCT(12)*V(98)

  B(20) = RCT(12)*V(94)

  B(21) = RCT(13)*V(102)

  B(22) = RCT(13)*V(94)

  B(23) = RCT(14)

  B(24) = RCT(15)

  B(25) = RCT(16)*V(99)

  B(26) = RCT(16)*V(95)

  B(27) = RCT(17)*V(105)

  B(28) = RCT(17)*V(95)

  B(29) = RCT(18)*V(105)

  B(30) = RCT(18)*V(95)

  B(31) = RCT(19)*F(1)

  B(33) = RCT(20)*F(1)*F(1)

  B(35) = RCT(21)

  B(36) = RCT(22)*2*V(99)

  B(37) = RCT(23)*V(105)*F(1)

  B(38) = RCT(23)*V(99)*F(1)

  B(40) = RCT(24)*V(99)

  B(41) = RCT(24)*V(98)

  B(42) = RCT(25)

  B(43) = RCT(26)*V(98)

  B(44) = RCT(26)*V(58)

  B(45) = RCT(27)*2*V(58)

  B(46) = RCT(28)*V(105)

  B(47) = RCT(28)*V(98)

  B(48) = RCT(29)*V(98)

  B(49) = RCT(29)*V(77)

  B(50) = RCT(30)*V(102)

  B(51) = RCT(30)*V(99)

  B(52) = RCT(31)*V(105)

  B(53) = RCT(31)*V(102)

  B(54) = RCT(32)

  B(55) = RCT(33)*V(98)

  B(56) = RCT(33)*V(64)

  B(57) = RCT(34)*2*V(102)

  B(58) = RCT(35)*2*V(102)*F(1)

  B(60) = RCT(36)

  B(61) = RCT(37)*V(98)

  B(62) = RCT(37)*V(68)

  B(63) = RCT(38)*V(53)

  B(64) = RCT(38)*V(52)

  B(65) = RCT(39)*V(98)

  B(66) = RCT(39)*V(52)

  B(67) = RCT(40)*V(101)

  B(68) = RCT(40)*V(98)

  B(69) = RCT(41)*2*V(98)

  B(70) = RCT(42)*2*V(98)

  B(71) = RCT(43)*V(102)

  B(72) = RCT(43)*V(98)

  B(73) = RCT(44)*V(102)

  B(74) = RCT(44)*V(101)

  B(75) = RCT(45)*V(101)

  B(76) = RCT(45)*V(68)

  B(77) = RCT(46)*V(101)

  B(78) = RCT(46)*V(95)

  B(79) = RCT(47)*V(98)

  B(80) = RCT(47)*V(95)

  B(81) = RCT(48)*V(102)

  B(82) = RCT(48)*V(95)

  B(83) = RCT(49)*V(95)

  B(84) = RCT(49)*V(94)

  B(85) = RCT(50)*2*V(95)

  B(86) = RCT(51)

  B(87) = RCT(52)

  B(88) = RCT(53)

  B(89) = RCT(54)*V(99)

  B(90) = RCT(54)*V(96)

  B(91) = RCT(55)*V(99)

  B(92) = RCT(55)*V(86)

  B(93) = RCT(56)*V(102)

  B(94) = RCT(56)*V(96)

  B(95) = RCT(57)*V(102)

  B(96) = RCT(57)*V(86)

  B(97) = RCT(58)*2*V(96)

  B(98) = RCT(59)*2*V(86)

  B(99) = RCT(60)*V(96)

  B(100) = RCT(60)*V(86)

  B(101) = RCT(61)*V(98)

  B(102) = RCT(61)*V(93)

  B(103) = RCT(62)

  B(104) = RCT(63)*V(98)

  B(105) = RCT(63)*V(72)

  B(106) = RCT(64)

  B(107) = RCT(65)*V(98)

  B(108) = RCT(65)*V(78)

  B(109) = RCT(66)*V(98)

  B(110) = RCT(66)*V(56)

  B(111) = RCT(67)*V(104)

  B(112) = RCT(67)*V(99)

  B(113) = RCT(68)*V(104)

  B(114) = RCT(68)*V(102)

  B(115) = RCT(69)*2*V(104)

  B(116) = RCT(70)*V(98)

  B(117) = RCT(70)*V(76)

  B(118) = RCT(71)

  B(119) = RCT(72)*V(98)

  B(120) = RCT(72)*V(66)

  B(121) = RCT(73)*V(98)

  B(122) = RCT(73)*V(92)

  B(123) = RCT(74)

  B(124) = RCT(75)

  B(125) = RCT(76)*V(101)

  B(126) = RCT(76)*V(92)

  B(127) = RCT(77)*V(95)

  B(128) = RCT(77)*V(92)

  B(129) = RCT(78)*V(102)

  B(130) = RCT(78)*V(92)

  B(131) = RCT(79)

  B(132) = RCT(80)*V(99)

  B(133) = RCT(80)*V(69)

  B(134) = RCT(81)*V(102)

  B(135) = RCT(81)*V(69)

  B(136) = RCT(82)*V(98)

  B(137) = RCT(82)*V(59)

  B(138) = RCT(83)*V(101)

  B(139) = RCT(83)*V(91)

  B(140) = RCT(84)*V(98)

  B(141) = RCT(84)*V(91)

  B(142) = RCT(85)*V(95)

  B(143) = RCT(85)*V(91)

  B(144) = RCT(86)

  B(145) = RCT(87)*V(99)

  B(146) = RCT(87)*V(97)

  B(147) = RCT(88)*V(105)

  B(148) = RCT(88)*V(97)

  B(149) = RCT(89)

  B(150) = RCT(90)

  B(151) = RCT(91)*V(102)

  B(152) = RCT(91)*V(97)

  B(153) = RCT(92)*V(104)

  B(154) = RCT(92)*V(97)

  B(155) = RCT(93)*V(97)

  B(156) = RCT(93)*V(96)

  B(157) = RCT(94)*2*V(97)

  B(158) = RCT(95)*V(98)

  B(159) = RCT(95)*V(61)

  B(160) = RCT(96)

  B(161) = RCT(97)*V(98)

  B(162) = RCT(97)*V(62)

  B(163) = RCT(98)*V(101)

  B(164) = RCT(98)*V(90)

  B(165) = RCT(99)*V(98)

  B(166) = RCT(99)*V(90)

  B(167) = RCT(100)*V(95)

  B(168) = RCT(100)*V(90)

  B(169) = RCT(101)

  B(170) = RCT(102)*V(103)

  B(171) = RCT(102)*V(99)

  B(172) = RCT(103)*V(105)

  B(173) = RCT(103)*V(103)

  B(174) = RCT(104)

  B(175) = RCT(105)

  B(176) = RCT(106)*V(98)

  B(177) = RCT(106)*V(63)

  B(178) = RCT(107)*V(103)

  B(179) = RCT(107)*V(102)

  B(180) = RCT(108)*V(104)

  B(181) = RCT(108)*V(103)

  B(182) = RCT(109)*V(103)

  B(183) = RCT(109)*V(96)

  B(184) = RCT(110)*2*V(103)

  B(185) = RCT(111)*V(103)

  B(186) = RCT(111)*V(97)

  B(187) = RCT(112)*V(98)

  B(188) = RCT(112)*V(87)

  B(189) = RCT(113)

  B(190) = 1600

  B(191) = RCT(115)*V(105)

  B(192) = RCT(115)*V(80)

  B(193) = RCT(116)*V(101)

  B(194) = RCT(116)*V(85)

  B(195) = RCT(117)*V(98)

  B(196) = RCT(117)*V(85)

  B(197) = RCT(118)*V(94)

  B(198) = RCT(118)*V(85)

  B(199) = RCT(119)*V(95)

  B(200) = RCT(119)*V(85)

  B(201) = RCT(120)*V(101)

  B(202) = RCT(120)*V(81)

  B(203) = RCT(121)*V(98)

  B(204) = RCT(121)*V(81)

  B(205) = RCT(122)*V(94)

  B(206) = RCT(122)*V(81)

  B(207) = RCT(123)*V(95)

  B(208) = RCT(123)*V(81)

  B(209) = RCT(124)*V(101)

  B(210) = RCT(124)*V(84)

  B(211) = RCT(125)*V(98)

  B(212) = RCT(125)*V(84)

  B(213) = RCT(126)*V(94)

  B(214) = RCT(126)*V(84)

  B(215) = RCT(127)*V(95)

  B(216) = RCT(127)*V(84)

  B(217) = RCT(128)*V(98)

  B(218) = RCT(128)*V(51)

  B(219) = RCT(129)*V(99)

  B(220) = RCT(129)*V(60)

  B(221) = 4.2

  B(222) = RCT(131)*V(98)

  B(223) = RCT(131)*V(82)

  B(224) = RCT(132)*V(95)

  B(225) = RCT(132)*V(82)

  B(226) = RCT(133)*V(105)

  B(227) = RCT(133)*V(71)

  B(228) = RCT(134)*V(102)

  B(229) = RCT(134)*V(71)

  B(230) = RCT(135)

  B(231) = RCT(136)*V(98)

  B(232) = RCT(136)*V(79)

  B(233) = RCT(137)*V(94)

  B(234) = RCT(137)*V(79)

  B(235) = RCT(138)*V(98)

  B(236) = RCT(138)*V(55)

  B(237) = RCT(139)*V(98)

  B(238) = RCT(139)*V(74)

  B(239) = RCT(140)

  B(240) = RCT(141)*V(101)

  B(241) = RCT(141)*V(88)

  B(242) = RCT(142)*V(98)

  B(243) = RCT(142)*V(88)

  B(244) = RCT(143)*V(94)

  B(245) = RCT(143)*V(88)

  B(246) = RCT(144)*V(95)

  B(247) = RCT(144)*V(88)

  B(248) = RCT(145)*V(98)

  B(249) = RCT(145)*V(89)

  B(250) = RCT(146)*V(94)

  B(251) = RCT(146)*V(89)

  B(252) = RCT(147)*V(95)

  B(253) = RCT(147)*V(89)

  B(254) = RCT(148)

  B(255) = RCT(149)*V(101)

  B(256) = RCT(149)*V(83)

  B(257) = RCT(150)*V(98)

  B(258) = RCT(150)*V(83)

  B(259) = RCT(151)*V(94)

  B(260) = RCT(151)*V(83)

  B(261) = RCT(152)*V(95)

  B(262) = RCT(152)*V(83)

  B(263) = RCT(153)*V(98)

  B(264) = RCT(153)*V(48)

  B(265) = RCT(154)*V(98)

  B(266) = RCT(154)*V(67)

  B(267) = RCT(155)*V(98)

  B(268) = RCT(155)*V(65)

  B(269) = RCT(156)*V(105)

  B(270) = RCT(156)*V(88)

  B(271) = RCT(157)

  B(272) = RCT(158)

  B(273) = RCT(159)*V(100)

  B(274) = RCT(159)*V(94)

  B(275) = RCT(160)*2*V(73)

  B(276) = RCT(161)*V(99)

  B(277) = RCT(161)*V(73)

  B(278) = RCT(162)*V(102)

  B(279) = RCT(162)*V(73)

  B(280) = RCT(163)*V(98)

  B(281) = RCT(163)*V(75)

  B(282) = RCT(164)

  B(283) = RCT(165)*V(100)

  B(284) = RCT(165)*V(56)

  B(285) = RCT(166)*V(100)

  B(286) = RCT(166)*V(87)

  B(287) = RCT(167)*V(100)

  B(288) = RCT(167)*V(65)

  B(289) = RCT(168)*V(100)

  B(290) = RCT(168)*V(81)

  B(291) = RCT(169)*V(100)

  B(292) = RCT(169)*V(85)

  B(293) = RCT(170)*V(100)

  B(294) = RCT(170)*V(84)

  B(295) = RCT(171)*V(100)

  B(296) = RCT(171)*V(88)

  B(297) = RCT(172)*V(100)

  B(298) = RCT(172)*V(92)

  B(299) = RCT(173)*V(100)

  B(300) = RCT(173)*V(91)

  B(301) = RCT(174)*V(100)

  B(302) = RCT(174)*V(90)

  B(303) = RCT(175)*V(100)

  B(304) = RCT(175)*V(66)

  B(305) = RCT(176)*V(100)

  B(306) = RCT(176)*V(67)

  B(307) = RCT(177)*V(98)

  B(308) = RCT(177)*V(70)

  B(309) = RCT(178)*V(94)

  B(310) = RCT(178)*V(57)

  B(311) = RCT(179)*V(98)

  B(312) = RCT(179)*V(57)

  B(313) = RCT(180)*V(68)

  B(314) = RCT(180)*V(57)

  B(315) = RCT(181)*V(98)

  B(316) = RCT(181)*V(12)

  B(317) = RCT(182)*V(98)

  B(318) = RCT(182)*V(15)

  B(319) = RCT(183)*V(98)

  B(320) = RCT(183)*V(18)

  B(321) = RCT(184)*V(98)

  B(322) = RCT(184)*V(23)

  B(323) = RCT(185)*V(94)

  B(324) = RCT(185)*V(23)

  B(325) = RCT(186)*V(98)

  B(326) = RCT(186)*V(29)

  B(327) = RCT(187)*V(94)

  B(328) = RCT(187)*V(29)

  B(329) = RCT(188)*V(95)

  B(330) = RCT(188)*V(29)

  B(331) = RCT(189)*V(98)

  B(332) = RCT(189)*V(32)

  B(333) = RCT(190)*V(98)

  B(334) = RCT(190)*V(34)

  B(335) = RCT(191)*V(98)

  B(336) = RCT(191)*V(37)

  B(337) = RCT(192)*V(98)

  B(338) = RCT(192)*V(42)

  B(339) = RCT(193)*V(98)

  B(340) = RCT(193)*V(41)

  B(341) = RCT(194)*V(98)

  B(342) = RCT(194)*V(40)

  B(343) = RCT(195)*V(98)

  B(344) = RCT(195)*V(46)

  B(345) = RCT(196)*V(98)

  B(346) = RCT(196)*V(45)

  B(347) = RCT(197)*V(98)

  B(348) = RCT(197)*V(44)

  B(349) = RCT(198)



  JVS(1) = 0

  JVS(2) = 0.071*B(217)

  JVS(3) = 0.071*B(218)

  JVS(4) = 0

  JVS(5) = 0.138*B(217)

  JVS(6) = 0.138*B(218)

  JVS(7) = 0

  JVS(8) = B(222)+B(224)

  JVS(9) = B(225)

  JVS(10) = B(223)

  JVS(11) = 0

  JVS(12) = 0.038*B(235)

  JVS(13) = 0.038*B(236)

  JVS(14) = 0

  JVS(15) = 0.167*B(235)

  JVS(16) = 0.167*B(236)

  JVS(17) = 0

  JVS(18) = 0.232*B(240)+0.232*B(242)+0.232*B(244)+0.232*B(246)

  JVS(19) = 0.232*B(245)

  JVS(20) = 0.232*B(247)

  JVS(21) = 0.232*B(243)

  JVS(22) = 0.232*B(241)

  JVS(23) = 0

  JVS(24) = 0.0228*B(240)+0.0288*B(242)+0.0288*B(244)+0.0288*B(246)

  JVS(25) = 0.0288*B(245)

  JVS(26) = 0.0288*B(247)

  JVS(27) = 0.0288*B(243)

  JVS(28) = 0.0228*B(241)

  JVS(29) = 0

  JVS(30) = B(263)+B(349)

  JVS(31) = B(264)

  JVS(32) = 0

  JVS(33) = B(263)

  JVS(34) = B(264)

  JVS(35) = 0

  JVS(36) = B(255)+B(257)+B(259)+B(261)

  JVS(37) = B(260)

  JVS(38) = B(262)

  JVS(39) = B(258)

  JVS(40) = B(256)

  JVS(41) = 0

  JVS(42) = B(315)

  JVS(43) = B(316)

  JVS(44) = -B(315)

  JVS(45) = -B(316)

  JVS(46) = 0

  JVS(47) = 0.239*B(317)

  JVS(48) = 0.239*B(318)

  JVS(49) = 0

  JVS(50) = 0.363*B(317)

  JVS(51) = 0.363*B(318)

  JVS(52) = -B(317)

  JVS(53) = -B(318)

  JVS(54) = 0

  JVS(55) = 0.045*B(319)

  JVS(56) = 0.045*B(320)

  JVS(57) = 0

  JVS(58) = 0.149*B(319)

  JVS(59) = 0.149*B(320)

  JVS(60) = -B(319)

  JVS(61) = -B(320)

  JVS(62) = 0

  JVS(63) = 0.038*B(321)

  JVS(64) = 0.038*B(322)

  JVS(65) = 0

  JVS(66) = 0.326*B(321)

  JVS(67) = 0.326*B(322)

  JVS(68) = 0

  JVS(69) = 0.125*B(323)

  JVS(70) = 0.125*B(324)

  JVS(71) = 0

  JVS(72) = 0.102*B(323)

  JVS(73) = 0.102*B(324)

  JVS(74) = -B(321)-B(323)

  JVS(75) = -B(324)

  JVS(76) = -B(322)

  JVS(77) = 0

  JVS(78) = 0.13*B(325)

  JVS(79) = 0.13*B(326)

  JVS(80) = 0

  JVS(81) = 0.0406*B(325)

  JVS(82) = 0.0406*B(326)

  JVS(83) = 0

  JVS(84) = 0.026*B(327)

  JVS(85) = 0.026*B(328)

  JVS(86) = 0

  JVS(87) = 0.485*B(327)

  JVS(88) = 0.485*B(328)

  JVS(89) = 0

  JVS(90) = B(329)

  JVS(91) = B(330)

  JVS(92) = -B(325)-B(327)-B(329)

  JVS(93) = -B(328)

  JVS(94) = -B(330)

  JVS(95) = -B(326)

  JVS(96) = 0

  JVS(97) = 0.091*B(331)

  JVS(98) = 0.091*B(332)

  JVS(99) = 0

  JVS(100) = 0.367*B(331)

  JVS(101) = 0.367*B(332)

  JVS(102) = -B(331)

  JVS(103) = -B(332)

  JVS(104) = 0

  JVS(105) = 1.173*B(333)

  JVS(106) = 1.173*B(334)

  JVS(107) = -B(333)

  JVS(108) = -B(334)

  JVS(109) = 0

  JVS(110) = 0.156*B(335)

  JVS(111) = 0.156*B(336)

  JVS(112) = 0

  JVS(113) = 0.777*B(335)

  JVS(114) = 0.777*B(336)

  JVS(115) = -B(335)

  JVS(116) = -B(336)

  JVS(117) = 0

  JVS(118) = B(309)+B(311)+B(313)

  JVS(119) = B(314)

  JVS(120) = B(310)

  JVS(121) = B(312)

  JVS(122) = 0

  JVS(123) = 1.075*B(341)

  JVS(124) = 1.075*B(342)

  JVS(125) = -B(341)

  JVS(126) = 1.075*B(339)

  JVS(127) = 1.075*B(340)-B(342)

  JVS(128) = -B(339)

  JVS(129) = 1.075*B(337)

  JVS(130) = 1.075*B(338)-B(340)

  JVS(131) = -B(337)

  JVS(132) = -B(338)

  JVS(133) = 0

  JVS(134) = 1.075*B(347)

  JVS(135) = 1.075*B(348)

  JVS(136) = -B(347)

  JVS(137) = 1.075*B(345)

  JVS(138) = 1.075*B(346)-B(348)

  JVS(139) = -B(345)

  JVS(140) = 1.075*B(343)

  JVS(141) = 1.075*B(344)-B(346)

  JVS(142) = -B(343)

  JVS(143) = -B(344)

  JVS(144) = -B(271)

  JVS(145) = 0.3*B(275)

  JVS(146) = -B(263)-B(349)

  JVS(147) = -B(264)

  JVS(148) = -B(149)-B(150)

  JVS(149) = B(147)

  JVS(150) = B(148)

  JVS(151) = -B(272)

  JVS(152) = B(278)

  JVS(153) = B(279)

  JVS(154) = -B(217)

  JVS(155) = -B(218)

  JVS(156) = -B(63)-B(65)

  JVS(157) = -B(64)

  JVS(158) = -B(66)

  JVS(159) = -B(63)

  JVS(160) = -B(15)-B(17)-B(64)

  JVS(161) = B(14)

  JVS(162) = 0

  JVS(163) = -B(31)-B(33)-B(35)-B(88)

  JVS(164) = B(29)

  JVS(165) = B(30)

  JVS(166) = -B(235)

  JVS(167) = -B(236)

  JVS(168) = -B(109)-B(283)

  JVS(169) = -B(110)

  JVS(170) = -B(284)

  JVS(171) = -B(309)-B(311)-B(313)

  JVS(172) = -B(314)

  JVS(173) = -B(310)

  JVS(174) = -B(312)

  JVS(175) = -B(42)-B(43)-2*B(45)

  JVS(176) = B(40)-B(44)

  JVS(177) = 2*B(37)+B(41)

  JVS(178) = 2*B(38)

  JVS(179) = -B(136)

  JVS(180) = B(132)

  JVS(181) = 0.37*B(205)

  JVS(182) = 0.37*B(206)

  JVS(183) = -B(137)

  JVS(184) = B(133)

  JVS(185) = 0.56*B(217)

  JVS(186) = 0.3*B(235)

  JVS(187) = -B(219)-B(221)

  JVS(188) = 0.56*B(218)+0.3*B(236)

  JVS(189) = -B(220)

  JVS(190) = -B(158)-B(160)

  JVS(191) = 0.8*B(151)

  JVS(192) = -B(159)

  JVS(193) = 0.8*B(152)+0.8*B(178)

  JVS(194) = 0.8*B(179)

  JVS(195) = -B(161)

  JVS(196) = 0.1*B(155)+0.1*B(182)

  JVS(197) = 0.2*B(151)+0.1*B(153)+0.1*B(156)

  JVS(198) = -B(162)

  JVS(199) = 0.2*B(152)+0.2*B(178)

  JVS(200) = 0.2*B(179)+0.1*B(180)+0.1*B(183)

  JVS(201) = 0.1*B(154)+0.1*B(181)

  JVS(202) = -B(174)-B(175)-B(176)

  JVS(203) = -B(177)

  JVS(204) = B(172)

  JVS(205) = B(173)

  JVS(206) = -B(54)-B(55)-B(86)

  JVS(207) = -B(56)

  JVS(208) = B(52)

  JVS(209) = B(53)

  JVS(210) = -B(267)-B(287)

  JVS(211) = -B(268)

  JVS(212) = -B(288)

  JVS(213) = -B(119)-B(303)

  JVS(214) = -B(120)

  JVS(215) = -B(304)

  JVS(216) = 0.63*B(115)

  JVS(217) = -B(265)-B(305)

  JVS(218) = -B(266)

  JVS(219) = -B(306)

  JVS(220) = -B(313)

  JVS(221) = -B(60)-B(61)-B(75)-B(314)

  JVS(222) = 0

  JVS(223) = -B(62)+B(70)

  JVS(224) = -B(76)

  JVS(225) = B(57)+B(58)

  JVS(226) = -B(131)-B(132)-B(134)

  JVS(227) = B(129)

  JVS(228) = -B(133)

  JVS(229) = B(130)-B(135)

  JVS(230) = B(283)

  JVS(231) = B(287)

  JVS(232) = B(303)

  JVS(233) = B(305)

  JVS(234) = -B(307)

  JVS(235) = 0.3*B(293)

  JVS(236) = B(285)

  JVS(237) = 0.15*B(295)

  JVS(238) = B(301)

  JVS(239) = B(299)

  JVS(240) = B(297)

  JVS(241) = -B(308)

  JVS(242) = B(284)+B(286)+B(288)+0.3*B(294)+0.15*B(296)+B(298)+B(300)+B(302)+B(304)+B(306)

  JVS(243) = 0

  JVS(244) = -B(226)-B(228)

  JVS(245) = 0.4*B(222)+B(224)

  JVS(246) = B(225)

  JVS(247) = 0.4*B(223)

  JVS(248) = -B(229)

  JVS(249) = -B(227)

  JVS(250) = -B(104)-B(106)

  JVS(251) = B(95)

  JVS(252) = B(93)

  JVS(253) = -B(105)

  JVS(254) = B(94)+B(96)

  JVS(255) = -2*B(275)-B(276)-B(278)

  JVS(256) = B(273)

  JVS(257) = -B(277)

  JVS(258) = B(274)

  JVS(259) = -B(279)

  JVS(260) = 0.8*B(235)

  JVS(261) = -B(237)-B(239)

  JVS(262) = 0.2*B(233)

  JVS(263) = 0.168*B(248)+0.85*B(250)

  JVS(264) = 0.2*B(234)+0.85*B(251)

  JVS(265) = 0.8*B(236)-B(238)+0.168*B(249)

  JVS(266) = -B(280)-B(282)

  JVS(267) = B(289)

  JVS(268) = 0.7*B(293)

  JVS(269) = B(291)

  JVS(270) = 0.85*B(295)

  JVS(271) = -B(281)

  JVS(272) = B(290)+B(292)+0.7*B(294)+0.85*B(296)

  JVS(273) = B(134)

  JVS(274) = -B(116)-B(118)

  JVS(275) = 0

  JVS(276) = -B(117)

  JVS(277) = 0

  JVS(278) = B(113)+B(135)

  JVS(279) = B(114)

  JVS(280) = 2*B(31)+2*B(33)

  JVS(281) = -B(48)-B(87)

  JVS(282) = B(224)

  JVS(283) = 0.15*B(252)

  JVS(284) = B(167)

  JVS(285) = B(142)

  JVS(286) = B(127)

  JVS(287) = B(101)

  JVS(288) = B(81)+B(128)+B(143)+B(168)+B(225)+0.15*B(253)

  JVS(289) = B(46)-B(49)+B(102)

  JVS(290) = B(82)

  JVS(291) = B(47)

  JVS(292) = B(239)

  JVS(293) = B(280)+B(282)

  JVS(294) = -B(107)

  JVS(295) = B(230)+2*B(231)+0.69*B(233)

  JVS(296) = B(201)+0.63*B(205)

  JVS(297) = 0.001*B(259)

  JVS(298) = 0.1*B(209)+0.25*B(213)

  JVS(299) = 0.2*B(193)+0.33*B(197)

  JVS(300) = 0.066*B(244)

  JVS(301) = 0.334*B(248)+0.225*B(250)+0.643*B(252)+0.333*B(254)

  JVS(302) = B(169)

  JVS(303) = B(144)

  JVS(304) = B(121)+B(123)+B(124)+B(125)+B(127)+B(297)

  JVS(305) = 0.33*B(198)+0.63*B(206)+0.25*B(214)+0.69*B(234)+0.066*B(245)+0.225*B(251)+0.001*B(260)

  JVS(306) = B(128)+0.643*B(253)

  JVS(307) = -B(108)+B(122)+2*B(232)+0.334*B(249)+B(281)

  JVS(308) = B(298)

  JVS(309) = B(126)+0.2*B(194)+B(202)+0.1*B(210)

  JVS(310) = 0.9*B(219)

  JVS(311) = -B(230)-B(231)-B(233)

  JVS(312) = 0.3*B(222)

  JVS(313) = -B(234)

  JVS(314) = 0.3*B(223)-B(232)

  JVS(315) = 0.9*B(220)

  JVS(316) = -0.98*B(189)-B(190)-B(191)

  JVS(317) = 0.76*B(187)+0.76*B(285)

  JVS(318) = 0.76*B(188)

  JVS(319) = 0.76*B(286)

  JVS(320) = -B(192)

  JVS(321) = -B(201)-B(203)-B(205)-B(207)-B(289)

  JVS(322) = -B(206)

  JVS(323) = -B(208)

  JVS(324) = -B(204)

  JVS(325) = -B(290)

  JVS(326) = -B(202)

  JVS(327) = 0.36*B(217)

  JVS(328) = 0.2*B(235)

  JVS(329) = B(221)

  JVS(330) = B(228)

  JVS(331) = -B(222)-B(224)

  JVS(332) = -B(225)

  JVS(333) = 0.36*B(218)-B(223)+0.2*B(236)

  JVS(334) = 0

  JVS(335) = B(229)

  JVS(336) = 0

  JVS(337) = -B(255)-B(257)-B(259)-B(261)

  JVS(338) = -B(260)

  JVS(339) = -B(262)

  JVS(340) = -B(258)

  JVS(341) = -B(256)

  JVS(342) = -B(209)-B(211)-B(213)-B(215)-B(293)

  JVS(343) = -B(214)

  JVS(344) = -B(216)

  JVS(345) = -B(212)

  JVS(346) = -B(294)

  JVS(347) = -B(210)

  JVS(348) = 0.3*B(293)

  JVS(349) = -B(193)-B(195)-B(197)-B(199)-B(291)

  JVS(350) = -B(198)

  JVS(351) = -B(200)

  JVS(352) = -B(196)

  JVS(353) = -B(292)+0.3*B(294)

  JVS(354) = -B(194)

  JVS(355) = 0.009*B(267)+0.009*B(287)

  JVS(356) = 0.04*B(189)

  JVS(357) = 0.25*B(257)+0.18*B(259)+0.25*B(261)

  JVS(358) = 0.01*B(193)+0.09*B(199)

  JVS(359) = -B(91)-B(95)-2*B(98)-B(99)

  JVS(360) = 0.13*B(187)+0.13*B(285)

  JVS(361) = 0.088*B(242)

  JVS(362) = 0.18*B(260)

  JVS(363) = 0.09*B(200)+0.25*B(262)

  JVS(364) = -B(100)

  JVS(365) = 0.13*B(188)+0.088*B(243)+0.25*B(258)+0.009*B(268)

  JVS(366) = -B(92)

  JVS(367) = 0.13*B(286)+0.009*B(288)

  JVS(368) = 0.01*B(194)

  JVS(369) = -B(96)

  JVS(370) = 0

  JVS(371) = 1.1*B(235)

  JVS(372) = -2.1*B(189)

  JVS(373) = 5.12*B(255)+1.66*B(257)+7*B(259)

  JVS(374) = 0.1*B(209)+0.3*B(293)

  JVS(375) = 0.2*B(193)-0.7*B(195)-B(197)-B(199)-B(291)

  JVS(376) = -1.11*B(187)-1.11*B(285)

  JVS(377) = 0.35*B(244)+2.4*B(269)

  JVS(378) = 1.565*B(248)+0.36*B(250)+1.282*B(252)+0.832*B(254)

  JVS(379) = -0.66*B(101)-0.66*B(103)

  JVS(380) = -B(198)+0.35*B(245)+0.36*B(251)+7*B(260)

  JVS(381) = -B(200)+1.282*B(253)

  JVS(382) = -0.66*B(102)-1.11*B(188)-0.7*B(196)+1.1*B(236)+1.565*B(249)+1.66*B(258)

  JVS(383) = -1.11*B(286)-B(292)+0.3*B(294)

  JVS(384) = 0.2*B(194)+0.1*B(210)+5.12*B(256)

  JVS(385) = 2.4*B(270)

  JVS(386) = -B(240)-B(242)-B(244)-B(246)-B(269)-B(295)

  JVS(387) = -B(245)

  JVS(388) = -B(247)

  JVS(389) = -B(243)

  JVS(390) = -B(296)

  JVS(391) = -B(241)

  JVS(392) = -B(270)

  JVS(393) = 0.75*B(240)+0.912*B(242)+0.65*B(244)+0.2*B(246)+0.2*B(269)+B(295)

  JVS(394) = -B(248)-B(250)-B(252)-B(254)

  JVS(395) = 0.65*B(245)-B(251)

  JVS(396) = 0.2*B(247)-B(253)

  JVS(397) = 0.912*B(243)-B(249)

  JVS(398) = B(296)

  JVS(399) = 0.75*B(241)

  JVS(400) = 0.2*B(270)

  JVS(401) = 0.05*B(265)

  JVS(402) = 0.5*B(104)+0.5*B(106)

  JVS(403) = 0.03*B(233)

  JVS(404) = 0.5*B(189)

  JVS(405) = 0.22*B(203)

  JVS(406) = 0

  JVS(407) = 0.47*B(257)+0.21*B(259)+0.47*B(261)

  JVS(408) = 0.66*B(209)+0.7*B(211)+0.35*B(213)+0.64*B(215)+0.55*B(293)

  JVS(409) = 0.3*B(193)+0.62*B(195)+0.32*B(197)+0.56*B(199)+0.67*B(291)

  JVS(410) = 0

  JVS(411) = 0.05*B(187)+0.05*B(285)

  JVS(412) = 0.15*B(244)+0.8*B(246)+0.8*B(269)

  JVS(413) = 0.12*B(248)+0.357*B(252)

  JVS(414) = -B(163)-B(165)-B(167)-B(169)-B(301)

  JVS(415) = 0.33*B(101)+0.33*B(103)

  JVS(416) = 0.32*B(198)+0.35*B(214)+0.03*B(234)+0.15*B(245)+0.21*B(260)

  JVS(417) = -B(168)+0.56*B(200)+0.64*B(216)+0.8*B(247)+0.357*B(253)+0.47*B(262)

  JVS(418) = 0

  JVS(419) = 0.33*B(102)+0.5*B(105)-B(166)+0.05*B(188)+0.62*B(196)+0.22*B(204)+0.7*B(212)+0.12*B(249)+0.47*B(258)+0.05&
               &*B(266)

  JVS(420) = 0

  JVS(421) = 0.05*B(286)+0.67*B(292)+0.55*B(294)-B(302)

  JVS(422) = -B(164)+0.3*B(194)+0.66*B(210)

  JVS(423) = 0

  JVS(424) = 0.8*B(270)

  JVS(425) = B(176)

  JVS(426) = 0.991*B(267)+0.991*B(287)

  JVS(427) = 0.9*B(265)+B(305)

  JVS(428) = 0.5*B(104)+0.5*B(106)

  JVS(429) = 0.6*B(189)

  JVS(430) = 1.24*B(209)+1.3*B(211)+0.65*B(213)+1.18*B(215)+0.45*B(293)

  JVS(431) = 0.2*B(193)+0.33*B(195)+0.18*B(197)+0.35*B(199)+0.33*B(291)

  JVS(432) = 0

  JVS(433) = 0.06*B(187)+0.06*B(285)

  JVS(434) = 0

  JVS(435) = 0.252*B(248)+0.02*B(250)+0.067*B(254)

  JVS(436) = -B(138)-B(140)-B(142)-B(144)-B(299)

  JVS(437) = 0.33*B(101)+0.33*B(103)

  JVS(438) = 0.18*B(198)+0.65*B(214)+0.02*B(251)

  JVS(439) = -B(143)+0.35*B(200)+1.18*B(216)

  JVS(440) = 0.9*B(182)

  JVS(441) = B(185)

  JVS(442) = 0.33*B(102)+0.5*B(105)-B(141)+B(177)+0.06*B(188)+0.33*B(196)+1.3*B(212)+0.252*B(249)+0.9*B(266)+0.991&
               &*B(268)

  JVS(443) = B(170)

  JVS(444) = 0.06*B(286)+0.991*B(288)+0.33*B(292)+0.45*B(294)-B(300)+B(306)

  JVS(445) = -B(139)+0.2*B(194)+1.24*B(210)

  JVS(446) = 0

  JVS(447) = B(171)+0.9*B(180)+0.9*B(183)+2*B(184)+B(186)

  JVS(448) = 0.9*B(181)

  JVS(449) = 0

  JVS(450) = B(119)+B(303)

  JVS(451) = 0.1*B(265)

  JVS(452) = B(131)

  JVS(453) = B(118)

  JVS(454) = B(231)+0.7*B(233)

  JVS(455) = B(201)+1.56*B(203)+B(205)+2*B(207)+B(289)

  JVS(456) = 0

  JVS(457) = 0.28*B(257)+0.24*B(259)

  JVS(458) = 0.25*B(213)

  JVS(459) = 0.2*B(193)+0.8*B(195)+0.74*B(197)+B(199)

  JVS(460) = 0.5*B(240)+0.629*B(242)+0.6*B(244)

  JVS(461) = 0.167*B(248)+0.15*B(250)+0.282*B(252)+0.9*B(254)

  JVS(462) = -B(121)-B(123)-B(124)-B(125)-B(127)-B(129)-B(297)

  JVS(463) = 0.33*B(101)+0.33*B(103)

  JVS(464) = 0.74*B(198)+B(206)+0.25*B(214)+0.7*B(234)+0.6*B(245)+0.15*B(251)+0.24*B(260)

  JVS(465) = -B(128)+B(200)+2*B(208)+0.282*B(253)

  JVS(466) = B(153)

  JVS(467) = 0.33*B(102)+B(120)-B(122)+0.8*B(196)+1.56*B(204)+B(232)+0.629*B(243)+0.167*B(249)+0.28*B(258)+0.1*B(266)

  JVS(468) = B(111)

  JVS(469) = B(290)-B(298)+B(304)

  JVS(470) = -B(126)+0.2*B(194)+B(202)+0.5*B(241)

  JVS(471) = -B(130)

  JVS(472) = 0.1*B(180)

  JVS(473) = B(112)+1.37*B(115)+B(154)+0.1*B(181)

  JVS(474) = 0

  JVS(475) = 0.1*B(219)

  JVS(476) = B(226)

  JVS(477) = B(191)

  JVS(478) = 0

  JVS(479) = 0.53*B(261)

  JVS(480) = B(91)

  JVS(481) = 0

  JVS(482) = 0.8*B(246)+0.8*B(269)

  JVS(483) = 0.85*B(252)

  JVS(484) = -B(101)-B(103)

  JVS(485) = 0

  JVS(486) = 0.8*B(247)+0.85*B(253)+0.53*B(262)

  JVS(487) = 0

  JVS(488) = -B(102)

  JVS(489) = B(92)+0.1*B(220)

  JVS(490) = 0

  JVS(491) = 0

  JVS(492) = 0

  JVS(493) = B(192)+B(227)+0.8*B(270)

  JVS(494) = -B(309)

  JVS(495) = 0

  JVS(496) = -B(233)

  JVS(497) = -B(205)

  JVS(498) = 0

  JVS(499) = -B(259)

  JVS(500) = -B(213)

  JVS(501) = -B(197)

  JVS(502) = -B(244)

  JVS(503) = -B(250)

  JVS(504) = -B(3)-B(11)-B(13)-B(14)-B(19)-B(21)-B(83)-B(198)-B(206)-B(214)-B(234)-B(245)-B(251)-B(260)-B(273)-B(310)

  JVS(505) = -B(84)

  JVS(506) = 0.2*B(151)

  JVS(507) = -B(20)

  JVS(508) = -B(4)

  JVS(509) = -B(274)

  JVS(510) = B(2)

  JVS(511) = -B(22)+0.2*B(152)+0.2*B(178)

  JVS(512) = 0.2*B(179)

  JVS(513) = -B(12)

  JVS(514) = B(35)+B(88)

  JVS(515) = 0.39*B(86)

  JVS(516) = B(48)

  JVS(517) = -B(207)

  JVS(518) = -B(224)

  JVS(519) = -B(261)

  JVS(520) = -B(215)

  JVS(521) = -B(199)

  JVS(522) = -B(246)

  JVS(523) = -B(252)

  JVS(524) = -B(167)

  JVS(525) = -B(142)

  JVS(526) = -B(127)

  JVS(527) = 0

  JVS(528) = B(11)-B(83)

  JVS(529) = -B(23)-B(24)-B(25)-B(27)-B(29)-B(77)-B(79)-B(81)-B(84)-2*B(85)-B(128)-B(143)-B(168)-B(200)-B(208)-B(216)&
               &-B(225)-B(247)-B(253)-B(262)

  JVS(530) = 0

  JVS(531) = 0

  JVS(532) = B(49)-B(80)

  JVS(533) = -B(26)

  JVS(534) = 0

  JVS(535) = B(7)-B(78)

  JVS(536) = -B(82)

  JVS(537) = 0

  JVS(538) = 0

  JVS(539) = B(8)+B(12)-B(28)-B(30)

  JVS(540) = 0.08*B(217)

  JVS(541) = 0.5*B(235)

  JVS(542) = 0.991*B(267)+0.991*B(287)

  JVS(543) = 0.1*B(265)

  JVS(544) = B(104)

  JVS(545) = B(237)

  JVS(546) = 0.3*B(116)

  JVS(547) = B(231)+0.03*B(233)

  JVS(548) = 0.96*B(189)

  JVS(549) = 0.7*B(201)+B(203)+B(207)+2*B(289)

  JVS(550) = 0.6*B(222)

  JVS(551) = 1.25*B(257)+0.76*B(259)+1.03*B(261)

  JVS(552) = 0.1*B(209)+B(211)+1.7*B(293)

  JVS(553) = 0.2*B(193)+0.8*B(195)+0.22*B(197)+0.91*B(199)+2*B(291)

  JVS(554) = -B(99)

  JVS(555) = 0.87*B(187)+0.87*B(285)

  JVS(556) = 0.25*B(240)+0.991*B(242)+0.2*B(244)+B(246)+B(269)+B(295)

  JVS(557) = 0.713*B(248)+0.064*B(250)+0.075*B(252)+0.7*B(254)

  JVS(558) = 0

  JVS(559) = 0

  JVS(560) = 0.22*B(198)+0.03*B(234)+0.2*B(245)+0.064*B(251)+0.76*B(260)

  JVS(561) = 0.91*B(200)+B(208)+B(247)+0.075*B(253)+1.03*B(262)

  JVS(562) = -B(89)-B(93)-2*B(97)-B(100)-B(155)-B(182)

  JVS(563) = -B(156)+B(185)

  JVS(564) = B(105)+0.3*B(117)+0.87*B(188)+0.8*B(196)+B(204)+B(212)+0.08*B(218)+0.6*B(223)+B(232)+0.5*B(236)+B(238)&
               &+0.991*B(243)+0.713*B(249)+1.25*B(258)+0.1*B(266)+0.991*B(268)

  JVS(565) = -B(90)+B(170)

  JVS(566) = 0.87*B(286)+0.991*B(288)+2*B(290)+2*B(292)+1.7*B(294)+B(296)

  JVS(567) = 0.2*B(194)+0.7*B(202)+0.1*B(210)+0.25*B(241)

  JVS(568) = -B(94)

  JVS(569) = B(171)+0.9*B(180)-B(183)+2*B(184)+B(186)

  JVS(570) = 0.9*B(181)

  JVS(571) = B(270)

  JVS(572) = B(149)+B(150)

  JVS(573) = B(158)

  JVS(574) = B(237)+B(239)

  JVS(575) = B(230)+B(231)+0.62*B(233)

  JVS(576) = 0

  JVS(577) = 0.21*B(248)+0.114*B(250)+0.967*B(254)

  JVS(578) = B(138)+B(140)+B(142)+B(299)

  JVS(579) = 0

  JVS(580) = 0.62*B(234)+0.114*B(251)

  JVS(581) = B(143)

  JVS(582) = -B(155)

  JVS(583) = -B(145)-B(147)-B(151)-B(153)-B(156)-2*B(157)-B(185)

  JVS(584) = B(141)+B(159)+B(232)+B(238)+0.21*B(249)

  JVS(585) = -B(146)

  JVS(586) = B(300)

  JVS(587) = B(139)

  JVS(588) = -B(152)

  JVS(589) = -B(186)

  JVS(590) = -B(154)

  JVS(591) = -B(148)

  JVS(592) = -B(263)

  JVS(593) = B(272)

  JVS(594) = -B(217)

  JVS(595) = B(63)-B(65)

  JVS(596) = 2*B(17)+B(64)

  JVS(597) = -B(235)

  JVS(598) = -B(109)

  JVS(599) = -B(311)

  JVS(600) = B(42)-B(43)

  JVS(601) = -B(136)

  JVS(602) = -B(158)+B(160)

  JVS(603) = -B(161)

  JVS(604) = -B(176)

  JVS(605) = -B(55)+0.39*B(86)

  JVS(606) = -B(267)

  JVS(607) = -B(119)

  JVS(608) = -B(265)

  JVS(609) = 2*B(60)-B(61)+B(75)

  JVS(610) = 0

  JVS(611) = -B(307)

  JVS(612) = -B(104)+B(106)

  JVS(613) = 0

  JVS(614) = -B(237)

  JVS(615) = -B(280)

  JVS(616) = -B(116)+B(118)

  JVS(617) = -B(48)+B(87)

  JVS(618) = -B(107)

  JVS(619) = -B(231)+0.08*B(233)

  JVS(620) = 0.3*B(201)-B(203)+0.13*B(205)

  JVS(621) = -B(222)

  JVS(622) = 0.15*B(255)-B(257)+0.57*B(259)

  JVS(623) = -B(211)+0.5*B(213)

  JVS(624) = 0.1*B(193)-B(195)+0.1*B(197)

  JVS(625) = 0

  JVS(626) = -B(187)

  JVS(627) = -B(242)+0.266*B(244)

  JVS(628) = -B(248)+0.268*B(250)

  JVS(629) = B(163)-B(165)

  JVS(630) = B(138)-B(140)

  JVS(631) = -B(121)+B(125)

  JVS(632) = -B(101)

  JVS(633) = -B(19)+B(21)+0.1*B(198)+0.13*B(206)+0.5*B(214)+0.08*B(234)+0.266*B(245)+0.268*B(251)+0.57*B(260)

  JVS(634) = -B(79)

  JVS(635) = 0

  JVS(636) = 0

  JVS(637) = -B(20)-B(40)-B(44)-B(46)-B(49)-B(56)-B(62)-B(66)-B(67)-2*B(69)-2*B(70)-B(71)-B(80)-B(102)-B(105)-B(108)&
               &-B(110)-B(117)-B(120)-B(122)-B(137)-B(141)-B(159)-B(162)-B(166)-B(177)-B(188)-B(196)-B(204)-B(212)-B(218)&
               &-B(223)-B(232)-B(236)-B(238)-B(243)-B(249)-B(258)-B(264)-B(266)-B(268)-B(281)-B(308)-B(312)

  JVS(638) = -B(41)+B(50)

  JVS(639) = 0

  JVS(640) = -B(68)+B(73)+B(76)+B(126)+B(139)+B(164)+0.1*B(194)+0.3*B(202)+0.15*B(256)

  JVS(641) = B(22)+B(51)-B(72)+B(74)

  JVS(642) = 0

  JVS(643) = 0

  JVS(644) = -B(47)

  JVS(645) = B(42)+B(45)

  JVS(646) = -B(219)

  JVS(647) = -B(132)

  JVS(648) = -B(276)

  JVS(649) = -B(91)

  JVS(650) = 0

  JVS(651) = 0.2*B(269)

  JVS(652) = 0

  JVS(653) = 0

  JVS(654) = 0

  JVS(655) = -B(3)

  JVS(656) = B(24)-B(25)+B(27)

  JVS(657) = -B(89)

  JVS(658) = -B(145)

  JVS(659) = -B(40)

  JVS(660) = -B(4)-B(9)-B(26)-2*B(36)-B(37)-B(41)-B(50)-B(90)-B(92)-B(111)-B(133)-B(146)-B(170)-B(220)-B(277)

  JVS(661) = 0

  JVS(662) = B(5)-B(10)

  JVS(663) = -B(51)

  JVS(664) = -B(171)

  JVS(665) = -B(112)

  JVS(666) = B(1)+B(6)+B(28)-B(38)+0.2*B(270)

  JVS(667) = 2*B(271)

  JVS(668) = B(272)

  JVS(669) = -B(283)

  JVS(670) = -B(287)

  JVS(671) = -B(303)

  JVS(672) = -B(305)

  JVS(673) = B(307)

  JVS(674) = 1.4*B(275)+B(276)

  JVS(675) = B(280)+B(282)

  JVS(676) = -B(289)

  JVS(677) = -B(293)

  JVS(678) = -B(291)

  JVS(679) = -B(285)

  JVS(680) = -B(295)

  JVS(681) = 0

  JVS(682) = -B(301)

  JVS(683) = -B(299)

  JVS(684) = -B(297)

  JVS(685) = 0

  JVS(686) = -B(273)

  JVS(687) = 0

  JVS(688) = 0

  JVS(689) = 0

  JVS(690) = B(281)+B(308)

  JVS(691) = B(277)

  JVS(692) = -B(274)-B(284)-B(286)-B(288)-B(290)-B(292)-B(294)-B(296)-B(298)-B(300)-B(302)-B(304)-B(306)

  JVS(693) = 0

  JVS(694) = 0

  JVS(695) = 0

  JVS(696) = 0

  JVS(697) = 0

  JVS(698) = B(15)

  JVS(699) = -B(75)

  JVS(700) = -B(201)

  JVS(701) = -B(255)

  JVS(702) = -B(209)+0.5*B(213)

  JVS(703) = -B(193)

  JVS(704) = -B(240)

  JVS(705) = -B(163)

  JVS(706) = -B(138)

  JVS(707) = -B(125)

  JVS(708) = 0

  JVS(709) = B(13)+0.5*B(214)

  JVS(710) = B(23)-B(77)

  JVS(711) = 0

  JVS(712) = 0

  JVS(713) = -B(67)+B(69)

  JVS(714) = -B(9)

  JVS(715) = 0

  JVS(716) = -B(2)-B(5)-B(7)-B(10)-B(68)-B(73)-B(76)-B(78)-B(126)-B(139)-B(164)-B(194)-B(202)-B(210)-B(241)-B(256)

  JVS(717) = -B(74)

  JVS(718) = 0

  JVS(719) = 0

  JVS(720) = B(1)-B(6)-B(8)

  JVS(721) = B(263)

  JVS(722) = 0.44*B(217)

  JVS(723) = B(63)+B(65)

  JVS(724) = B(64)

  JVS(725) = 0.7*B(235)

  JVS(726) = B(136)

  JVS(727) = 0.9*B(219)+B(221)

  JVS(728) = B(54)+0.61*B(86)

  JVS(729) = B(267)+B(287)

  JVS(730) = B(119)+B(303)

  JVS(731) = B(265)+B(305)

  JVS(732) = B(61)+B(75)

  JVS(733) = B(131)+B(132)-B(134)

  JVS(734) = -B(228)

  JVS(735) = B(106)

  JVS(736) = -B(278)

  JVS(737) = B(239)

  JVS(738) = B(282)

  JVS(739) = 0.3*B(116)+B(118)

  JVS(740) = B(107)

  JVS(741) = B(230)+2*B(231)+0.76*B(233)

  JVS(742) = 0.94*B(189)+B(190)

  JVS(743) = 1.7*B(201)+B(203)+0.13*B(205)+B(289)

  JVS(744) = 0.6*B(222)

  JVS(745) = 0.75*B(257)+0.07*B(259)+0.28*B(261)

  JVS(746) = 0.1*B(209)+B(211)+0.5*B(213)+B(215)+B(293)

  JVS(747) = 0.3*B(193)+0.95*B(195)+0.44*B(197)+B(291)

  JVS(748) = -B(95)

  JVS(749) = 0.11*B(187)+0.11*B(285)

  JVS(750) = 0.25*B(240)+0.912*B(242)+0.066*B(244)+0.8*B(246)+0.8*B(269)+B(295)

  JVS(751) = 0.503*B(248)+0.154*B(250)+0.925*B(252)+1.033*B(254)

  JVS(752) = B(169)

  JVS(753) = B(144)

  JVS(754) = B(121)+2*B(123)+B(125)+B(127)-B(129)+B(297)

  JVS(755) = B(101)+B(103)

  JVS(756) = B(19)-B(21)+0.44*B(198)+0.13*B(206)+0.5*B(214)+0.76*B(234)+0.066*B(245)+0.154*B(251)+0.07*B(260)

  JVS(757) = B(79)-B(81)+B(128)+B(216)+0.8*B(247)+0.925*B(253)+0.28*B(262)

  JVS(758) = -B(93)

  JVS(759) = -B(151)+0.9*B(153)+B(185)

  JVS(760) = B(20)+B(62)+B(66)+B(67)-B(71)+B(80)+B(102)+B(108)+0.3*B(117)+B(120)+B(122)+B(137)+0.11*B(188)+0.95*B(196)&
               &+B(204)+B(212)+0.44*B(218)+0.6*B(223)+2*B(232)+0.7*B(236)+0.912*B(243)+0.503*B(249)+0.75*B(258)+B(264)&
               &+B(266)+B(268)

  JVS(761) = -B(50)+B(111)+B(133)+B(170)+0.9*B(220)

  JVS(762) = 0.11*B(286)+B(288)+B(290)+B(292)+B(294)+B(296)+B(298)+B(304)+B(306)

  JVS(763) = B(68)-B(73)+B(76)+B(126)+0.3*B(194)+1.7*B(202)+0.1*B(210)+0.25*B(241)

  JVS(764) = -B(22)-B(51)-B(52)-2*B(57)-2*B(58)-B(72)-B(74)-B(82)-B(94)-B(96)-B(113)-B(130)-B(135)-B(152)-B(178)-B(229)&
               &-B(279)

  JVS(765) = B(171)-B(179)+B(180)+2*B(184)+B(186)

  JVS(766) = B(112)-B(114)+0.74*B(115)+0.9*B(154)+B(181)

  JVS(767) = -B(53)+0.8*B(270)

  JVS(768) = B(174)+B(175)

  JVS(769) = 0.39*B(259)

  JVS(770) = 0.25*B(240)+0.2*B(244)

  JVS(771) = 0.25*B(248)+0.075*B(252)

  JVS(772) = B(163)+B(165)+B(167)+B(301)

  JVS(773) = 0

  JVS(774) = 0.2*B(245)+0.39*B(260)

  JVS(775) = B(168)+0.075*B(253)

  JVS(776) = -B(182)

  JVS(777) = -B(185)

  JVS(778) = B(166)+0.25*B(249)

  JVS(779) = -B(170)

  JVS(780) = B(302)

  JVS(781) = B(164)+0.25*B(241)

  JVS(782) = -B(178)

  JVS(783) = -B(171)-B(172)-B(179)-B(180)-B(183)-2*B(184)-B(186)

  JVS(784) = -B(181)

  JVS(785) = -B(173)

  JVS(786) = B(109)+B(283)

  JVS(787) = B(160)

  JVS(788) = B(161)

  JVS(789) = 0.7*B(116)

  JVS(790) = B(169)

  JVS(791) = B(144)

  JVS(792) = 0

  JVS(793) = 0

  JVS(794) = 0

  JVS(795) = 0

  JVS(796) = 0.9*B(155)

  JVS(797) = B(145)-0.1*B(153)+0.9*B(156)+2*B(157)+B(185)

  JVS(798) = B(110)+0.7*B(117)+B(162)

  JVS(799) = -B(111)+B(146)

  JVS(800) = B(284)

  JVS(801) = 0

  JVS(802) = -B(113)

  JVS(803) = -B(180)+B(186)

  JVS(804) = -B(112)-B(114)-2*B(115)-0.1*B(154)-B(181)

  JVS(805) = 0

  JVS(806) = B(149)+B(150)

  JVS(807) = B(35)+B(88)

  JVS(808) = B(43)+B(45)

  JVS(809) = 0.9*B(219)

  JVS(810) = B(174)+B(175)+B(176)

  JVS(811) = B(54)+B(55)+0.61*B(86)

  JVS(812) = B(132)

  JVS(813) = -B(226)

  JVS(814) = B(276)

  JVS(815) = B(87)

  JVS(816) = -B(191)

  JVS(817) = B(207)

  JVS(818) = 0

  JVS(819) = 0.47*B(261)

  JVS(820) = B(215)

  JVS(821) = B(199)

  JVS(822) = 0

  JVS(823) = 0.2*B(246)-B(269)

  JVS(824) = 0

  JVS(825) = 0

  JVS(826) = 0

  JVS(827) = 0

  JVS(828) = B(103)

  JVS(829) = B(3)-B(11)+B(83)

  JVS(830) = B(23)+2*B(25)-B(29)+B(77)+B(79)+B(84)+2*B(85)+B(200)+B(208)+B(216)+0.2*B(247)+0.47*B(262)

  JVS(831) = B(89)

  JVS(832) = B(145)-B(147)

  JVS(833) = B(44)-B(46)+B(56)+B(80)+B(177)

  JVS(834) = B(4)+B(9)+2*B(26)+2*B(36)-B(37)+B(50)+B(90)+B(111)+B(133)+B(146)+B(170)+0.9*B(220)+B(277)

  JVS(835) = 0

  JVS(836) = -B(5)-B(7)+B(10)+B(78)

  JVS(837) = B(51)-B(52)

  JVS(838) = B(171)-B(172)

  JVS(839) = B(112)

  JVS(840) = -B(1)-B(6)-B(8)-B(12)-B(30)-B(38)-B(47)-B(53)-B(148)-B(173)-B(192)-B(227)-B(270)
      
END SUBROUTINE cb05_sorg_vbs_aq_Jac_SP














SUBROUTINE cb05_sorg_vbs_aq_KppDecomp( JVS, IER )







      INTEGER  :: IER
      REAL(kind=dp) :: JVS(840), W(105), a
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
      
END SUBROUTINE cb05_sorg_vbs_aq_KppDecomp



SUBROUTINE cb05_sorg_vbs_aq_KppDecompCmplx( JVS, IER )







      INTEGER  :: IER
      DOUBLE COMPLEX :: JVS(840), W(105), a
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
      
END SUBROUTINE cb05_sorg_vbs_aq_KppDecompCmplx


SUBROUTINE cb05_sorg_vbs_aq_KppSolveIndirect( JVS, X )







      INTEGER i, j
      REAL(kind=dp) JVS(840), X(105), sum

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
      
END SUBROUTINE cb05_sorg_vbs_aq_KppSolveIndirect


SUBROUTINE cb05_sorg_vbs_aq_KppSolveCmplx( JVS, X )







      INTEGER i, j
      DOUBLE COMPLEX JVS(840), X(105), sum

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
      
END SUBROUTINE cb05_sorg_vbs_aq_KppSolveCmplx













SUBROUTINE cb05_sorg_vbs_aq_KppSolve ( JVS, X )


  REAL(kind=dp) :: JVS(LU_NONZERO)

  REAL(kind=dp) :: X(NVAR)

  X(53) = X(53)-JVS(159)*X(52)
  X(60) = X(60)-JVS(185)*X(51)-JVS(186)*X(55)
  X(68) = X(68)-JVS(220)*X(57)
  X(70) = X(70)-JVS(230)*X(56)-JVS(231)*X(65)-JVS(232)*X(66)-JVS(233)*X(67)
  X(74) = X(74)-JVS(260)*X(55)
  X(76) = X(76)-JVS(273)*X(69)
  X(77) = X(77)-JVS(280)*X(54)
  X(78) = X(78)-JVS(292)*X(74)-JVS(293)*X(75)
  X(79) = X(79)-JVS(310)*X(60)
  X(82) = X(82)-JVS(327)*X(51)-JVS(328)*X(55)-JVS(329)*X(60)-JVS(330)*X(71)
  X(85) = X(85)-JVS(348)*X(84)
  X(86) = X(86)-JVS(355)*X(65)-JVS(356)*X(80)-JVS(357)*X(83)-JVS(358)*X(85)
  X(87) = X(87)-JVS(371)*X(55)-JVS(372)*X(80)-JVS(373)*X(83)-JVS(374)*X(84)-JVS(375)*X(85)
  X(89) = X(89)-JVS(393)*X(88)
  X(90) = X(90)-JVS(401)*X(67)-JVS(402)*X(72)-JVS(403)*X(79)-JVS(404)*X(80)-JVS(405)*X(81)-JVS(406)*X(82)-JVS(407)*X(83)&
            &-JVS(408)*X(84)-JVS(409)*X(85)-JVS(410)*X(86)-JVS(411)*X(87)-JVS(412)*X(88)-JVS(413)*X(89)
  X(91) = X(91)-JVS(425)*X(63)-JVS(426)*X(65)-JVS(427)*X(67)-JVS(428)*X(72)-JVS(429)*X(80)-JVS(430)*X(84)-JVS(431)*X(85)&
            &-JVS(432)*X(86)-JVS(433)*X(87)-JVS(434)*X(88)-JVS(435)*X(89)
  X(92) = X(92)-JVS(450)*X(66)-JVS(451)*X(67)-JVS(452)*X(69)-JVS(453)*X(76)-JVS(454)*X(79)-JVS(455)*X(81)-JVS(456)*X(82)&
            &-JVS(457)*X(83)-JVS(458)*X(84)-JVS(459)*X(85)-JVS(460)*X(88)-JVS(461)*X(89)
  X(93) = X(93)-JVS(475)*X(60)-JVS(476)*X(71)-JVS(477)*X(80)-JVS(478)*X(82)-JVS(479)*X(83)-JVS(480)*X(86)-JVS(481)*X(87)&
            &-JVS(482)*X(88)-JVS(483)*X(89)
  X(94) = X(94)-JVS(494)*X(57)-JVS(495)*X(68)-JVS(496)*X(79)-JVS(497)*X(81)-JVS(498)*X(82)-JVS(499)*X(83)-JVS(500)*X(84)&
            &-JVS(501)*X(85)-JVS(502)*X(88)-JVS(503)*X(89)
  X(95) = X(95)-JVS(514)*X(54)-JVS(515)*X(64)-JVS(516)*X(77)-JVS(517)*X(81)-JVS(518)*X(82)-JVS(519)*X(83)-JVS(520)*X(84)&
            &-JVS(521)*X(85)-JVS(522)*X(88)-JVS(523)*X(89)-JVS(524)*X(90)-JVS(525)*X(91)-JVS(526)*X(92)-JVS(527)*X(93)&
            &-JVS(528)*X(94)
  X(96) = X(96)-JVS(540)*X(51)-JVS(541)*X(55)-JVS(542)*X(65)-JVS(543)*X(67)-JVS(544)*X(72)-JVS(545)*X(74)-JVS(546)*X(76)&
            &-JVS(547)*X(79)-JVS(548)*X(80)-JVS(549)*X(81)-JVS(550)*X(82)-JVS(551)*X(83)-JVS(552)*X(84)-JVS(553)*X(85)&
            &-JVS(554)*X(86)-JVS(555)*X(87)-JVS(556)*X(88)-JVS(557)*X(89)-JVS(558)*X(92)-JVS(559)*X(93)-JVS(560)*X(94)&
            &-JVS(561)*X(95)
  X(97) = X(97)-JVS(572)*X(49)-JVS(573)*X(61)-JVS(574)*X(74)-JVS(575)*X(79)-JVS(576)*X(82)-JVS(577)*X(89)-JVS(578)*X(91)&
            &-JVS(579)*X(93)-JVS(580)*X(94)-JVS(581)*X(95)-JVS(582)*X(96)
  X(98) = X(98)-JVS(592)*X(48)-JVS(593)*X(50)-JVS(594)*X(51)-JVS(595)*X(52)-JVS(596)*X(53)-JVS(597)*X(55)-JVS(598)*X(56)&
            &-JVS(599)*X(57)-JVS(600)*X(58)-JVS(601)*X(59)-JVS(602)*X(61)-JVS(603)*X(62)-JVS(604)*X(63)-JVS(605)*X(64)&
            &-JVS(606)*X(65)-JVS(607)*X(66)-JVS(608)*X(67)-JVS(609)*X(68)-JVS(610)*X(69)-JVS(611)*X(70)-JVS(612)*X(72)&
            &-JVS(613)*X(73)-JVS(614)*X(74)-JVS(615)*X(75)-JVS(616)*X(76)-JVS(617)*X(77)-JVS(618)*X(78)-JVS(619)*X(79)&
            &-JVS(620)*X(81)-JVS(621)*X(82)-JVS(622)*X(83)-JVS(623)*X(84)-JVS(624)*X(85)-JVS(625)*X(86)-JVS(626)*X(87)&
            &-JVS(627)*X(88)-JVS(628)*X(89)-JVS(629)*X(90)-JVS(630)*X(91)-JVS(631)*X(92)-JVS(632)*X(93)-JVS(633)*X(94)&
            &-JVS(634)*X(95)-JVS(635)*X(96)-JVS(636)*X(97)
  X(99) = X(99)-JVS(645)*X(58)-JVS(646)*X(60)-JVS(647)*X(69)-JVS(648)*X(73)-JVS(649)*X(86)-JVS(650)*X(87)-JVS(651)*X(88)&
            &-JVS(652)*X(89)-JVS(653)*X(92)-JVS(654)*X(93)-JVS(655)*X(94)-JVS(656)*X(95)-JVS(657)*X(96)-JVS(658)*X(97)&
            &-JVS(659)*X(98)
  X(100) = X(100)-JVS(667)*X(47)-JVS(668)*X(50)-JVS(669)*X(56)-JVS(670)*X(65)-JVS(671)*X(66)-JVS(672)*X(67)-JVS(673)&
             &*X(70)-JVS(674)*X(73)-JVS(675)*X(75)-JVS(676)*X(81)-JVS(677)*X(84)-JVS(678)*X(85)-JVS(679)*X(87)-JVS(680)&
             &*X(88)-JVS(681)*X(89)-JVS(682)*X(90)-JVS(683)*X(91)-JVS(684)*X(92)-JVS(685)*X(93)-JVS(686)*X(94)-JVS(687)&
             &*X(95)-JVS(688)*X(96)-JVS(689)*X(97)-JVS(690)*X(98)-JVS(691)*X(99)
  X(101) = X(101)-JVS(698)*X(53)-JVS(699)*X(68)-JVS(700)*X(81)-JVS(701)*X(83)-JVS(702)*X(84)-JVS(703)*X(85)-JVS(704)&
             &*X(88)-JVS(705)*X(90)-JVS(706)*X(91)-JVS(707)*X(92)-JVS(708)*X(93)-JVS(709)*X(94)-JVS(710)*X(95)-JVS(711)&
             &*X(96)-JVS(712)*X(97)-JVS(713)*X(98)-JVS(714)*X(99)-JVS(715)*X(100)
  X(102) = X(102)-JVS(721)*X(48)-JVS(722)*X(51)-JVS(723)*X(52)-JVS(724)*X(53)-JVS(725)*X(55)-JVS(726)*X(59)-JVS(727)&
             &*X(60)-JVS(728)*X(64)-JVS(729)*X(65)-JVS(730)*X(66)-JVS(731)*X(67)-JVS(732)*X(68)-JVS(733)*X(69)-JVS(734)&
             &*X(71)-JVS(735)*X(72)-JVS(736)*X(73)-JVS(737)*X(74)-JVS(738)*X(75)-JVS(739)*X(76)-JVS(740)*X(78)-JVS(741)&
             &*X(79)-JVS(742)*X(80)-JVS(743)*X(81)-JVS(744)*X(82)-JVS(745)*X(83)-JVS(746)*X(84)-JVS(747)*X(85)-JVS(748)&
             &*X(86)-JVS(749)*X(87)-JVS(750)*X(88)-JVS(751)*X(89)-JVS(752)*X(90)-JVS(753)*X(91)-JVS(754)*X(92)-JVS(755)&
             &*X(93)-JVS(756)*X(94)-JVS(757)*X(95)-JVS(758)*X(96)-JVS(759)*X(97)-JVS(760)*X(98)-JVS(761)*X(99)-JVS(762)&
             &*X(100)-JVS(763)*X(101)
  X(103) = X(103)-JVS(768)*X(63)-JVS(769)*X(83)-JVS(770)*X(88)-JVS(771)*X(89)-JVS(772)*X(90)-JVS(773)*X(93)-JVS(774)&
             &*X(94)-JVS(775)*X(95)-JVS(776)*X(96)-JVS(777)*X(97)-JVS(778)*X(98)-JVS(779)*X(99)-JVS(780)*X(100)-JVS(781)&
             &*X(101)-JVS(782)*X(102)
  X(104) = X(104)-JVS(786)*X(56)-JVS(787)*X(61)-JVS(788)*X(62)-JVS(789)*X(76)-JVS(790)*X(90)-JVS(791)*X(91)-JVS(792)&
             &*X(92)-JVS(793)*X(93)-JVS(794)*X(94)-JVS(795)*X(95)-JVS(796)*X(96)-JVS(797)*X(97)-JVS(798)*X(98)-JVS(799)&
             &*X(99)-JVS(800)*X(100)-JVS(801)*X(101)-JVS(802)*X(102)-JVS(803)*X(103)
  X(105) = X(105)-JVS(806)*X(49)-JVS(807)*X(54)-JVS(808)*X(58)-JVS(809)*X(60)-JVS(810)*X(63)-JVS(811)*X(64)-JVS(812)&
             &*X(69)-JVS(813)*X(71)-JVS(814)*X(73)-JVS(815)*X(77)-JVS(816)*X(80)-JVS(817)*X(81)-JVS(818)*X(82)-JVS(819)&
             &*X(83)-JVS(820)*X(84)-JVS(821)*X(85)-JVS(822)*X(87)-JVS(823)*X(88)-JVS(824)*X(89)-JVS(825)*X(90)-JVS(826)&
             &*X(91)-JVS(827)*X(92)-JVS(828)*X(93)-JVS(829)*X(94)-JVS(830)*X(95)-JVS(831)*X(96)-JVS(832)*X(97)-JVS(833)&
             &*X(98)-JVS(834)*X(99)-JVS(835)*X(100)-JVS(836)*X(101)-JVS(837)*X(102)-JVS(838)*X(103)-JVS(839)*X(104)
  X(105) = X(105)/JVS(840)
  X(104) = (X(104)-JVS(805)*X(105))/(JVS(804))
  X(103) = (X(103)-JVS(784)*X(104)-JVS(785)*X(105))/(JVS(783))
  X(102) = (X(102)-JVS(765)*X(103)-JVS(766)*X(104)-JVS(767)*X(105))/(JVS(764))
  X(101) = (X(101)-JVS(717)*X(102)-JVS(718)*X(103)-JVS(719)*X(104)-JVS(720)*X(105))/(JVS(716))
  X(100) = (X(100)-JVS(693)*X(101)-JVS(694)*X(102)-JVS(695)*X(103)-JVS(696)*X(104)-JVS(697)*X(105))/(JVS(692))
  X(99) = (X(99)-JVS(661)*X(100)-JVS(662)*X(101)-JVS(663)*X(102)-JVS(664)*X(103)-JVS(665)*X(104)-JVS(666)*X(105))&
            &/(JVS(660))
  X(98) = (X(98)-JVS(638)*X(99)-JVS(639)*X(100)-JVS(640)*X(101)-JVS(641)*X(102)-JVS(642)*X(103)-JVS(643)*X(104)-JVS(644)&
            &*X(105))/(JVS(637))
  X(97) = (X(97)-JVS(584)*X(98)-JVS(585)*X(99)-JVS(586)*X(100)-JVS(587)*X(101)-JVS(588)*X(102)-JVS(589)*X(103)-JVS(590)&
            &*X(104)-JVS(591)*X(105))/(JVS(583))
  X(96) = (X(96)-JVS(563)*X(97)-JVS(564)*X(98)-JVS(565)*X(99)-JVS(566)*X(100)-JVS(567)*X(101)-JVS(568)*X(102)-JVS(569)&
            &*X(103)-JVS(570)*X(104)-JVS(571)*X(105))/(JVS(562))
  X(95) = (X(95)-JVS(530)*X(96)-JVS(531)*X(97)-JVS(532)*X(98)-JVS(533)*X(99)-JVS(534)*X(100)-JVS(535)*X(101)-JVS(536)&
            &*X(102)-JVS(537)*X(103)-JVS(538)*X(104)-JVS(539)*X(105))/(JVS(529))
  X(94) = (X(94)-JVS(505)*X(95)-JVS(506)*X(97)-JVS(507)*X(98)-JVS(508)*X(99)-JVS(509)*X(100)-JVS(510)*X(101)-JVS(511)&
            &*X(102)-JVS(512)*X(103)-JVS(513)*X(105))/(JVS(504))
  X(93) = (X(93)-JVS(485)*X(94)-JVS(486)*X(95)-JVS(487)*X(96)-JVS(488)*X(98)-JVS(489)*X(99)-JVS(490)*X(100)-JVS(491)&
            &*X(101)-JVS(492)*X(102)-JVS(493)*X(105))/(JVS(484))
  X(92) = (X(92)-JVS(463)*X(93)-JVS(464)*X(94)-JVS(465)*X(95)-JVS(466)*X(97)-JVS(467)*X(98)-JVS(468)*X(99)-JVS(469)&
            &*X(100)-JVS(470)*X(101)-JVS(471)*X(102)-JVS(472)*X(103)-JVS(473)*X(104)-JVS(474)*X(105))/(JVS(462))
  X(91) = (X(91)-JVS(437)*X(93)-JVS(438)*X(94)-JVS(439)*X(95)-JVS(440)*X(96)-JVS(441)*X(97)-JVS(442)*X(98)-JVS(443)&
            &*X(99)-JVS(444)*X(100)-JVS(445)*X(101)-JVS(446)*X(102)-JVS(447)*X(103)-JVS(448)*X(104)-JVS(449)*X(105))&
            &/(JVS(436))
  X(90) = (X(90)-JVS(415)*X(93)-JVS(416)*X(94)-JVS(417)*X(95)-JVS(418)*X(96)-JVS(419)*X(98)-JVS(420)*X(99)-JVS(421)&
            &*X(100)-JVS(422)*X(101)-JVS(423)*X(102)-JVS(424)*X(105))/(JVS(414))
  X(89) = (X(89)-JVS(395)*X(94)-JVS(396)*X(95)-JVS(397)*X(98)-JVS(398)*X(100)-JVS(399)*X(101)-JVS(400)*X(105))&
            &/(JVS(394))
  X(88) = (X(88)-JVS(387)*X(94)-JVS(388)*X(95)-JVS(389)*X(98)-JVS(390)*X(100)-JVS(391)*X(101)-JVS(392)*X(105))&
            &/(JVS(386))
  X(87) = (X(87)-JVS(377)*X(88)-JVS(378)*X(89)-JVS(379)*X(93)-JVS(380)*X(94)-JVS(381)*X(95)-JVS(382)*X(98)-JVS(383)&
            &*X(100)-JVS(384)*X(101)-JVS(385)*X(105))/(JVS(376))
  X(86) = (X(86)-JVS(360)*X(87)-JVS(361)*X(88)-JVS(362)*X(94)-JVS(363)*X(95)-JVS(364)*X(96)-JVS(365)*X(98)-JVS(366)&
            &*X(99)-JVS(367)*X(100)-JVS(368)*X(101)-JVS(369)*X(102)-JVS(370)*X(105))/(JVS(359))
  X(85) = (X(85)-JVS(350)*X(94)-JVS(351)*X(95)-JVS(352)*X(98)-JVS(353)*X(100)-JVS(354)*X(101))/(JVS(349))
  X(84) = (X(84)-JVS(343)*X(94)-JVS(344)*X(95)-JVS(345)*X(98)-JVS(346)*X(100)-JVS(347)*X(101))/(JVS(342))
  X(83) = (X(83)-JVS(338)*X(94)-JVS(339)*X(95)-JVS(340)*X(98)-JVS(341)*X(101))/(JVS(337))
  X(82) = (X(82)-JVS(332)*X(95)-JVS(333)*X(98)-JVS(334)*X(99)-JVS(335)*X(102)-JVS(336)*X(105))/(JVS(331))
  X(81) = (X(81)-JVS(322)*X(94)-JVS(323)*X(95)-JVS(324)*X(98)-JVS(325)*X(100)-JVS(326)*X(101))/(JVS(321))
  X(80) = (X(80)-JVS(317)*X(87)-JVS(318)*X(98)-JVS(319)*X(100)-JVS(320)*X(105))/(JVS(316))
  X(79) = (X(79)-JVS(312)*X(82)-JVS(313)*X(94)-JVS(314)*X(98)-JVS(315)*X(99))/(JVS(311))
  X(78) = (X(78)-JVS(295)*X(79)-JVS(296)*X(81)-JVS(297)*X(83)-JVS(298)*X(84)-JVS(299)*X(85)-JVS(300)*X(88)-JVS(301)&
            &*X(89)-JVS(302)*X(90)-JVS(303)*X(91)-JVS(304)*X(92)-JVS(305)*X(94)-JVS(306)*X(95)-JVS(307)*X(98)-JVS(308)&
            &*X(100)-JVS(309)*X(101))/(JVS(294))
  X(77) = (X(77)-JVS(282)*X(82)-JVS(283)*X(89)-JVS(284)*X(90)-JVS(285)*X(91)-JVS(286)*X(92)-JVS(287)*X(93)-JVS(288)&
            &*X(95)-JVS(289)*X(98)-JVS(290)*X(102)-JVS(291)*X(105))/(JVS(281))
  X(76) = (X(76)-JVS(275)*X(92)-JVS(276)*X(98)-JVS(277)*X(99)-JVS(278)*X(102)-JVS(279)*X(104))/(JVS(274))
  X(75) = (X(75)-JVS(267)*X(81)-JVS(268)*X(84)-JVS(269)*X(85)-JVS(270)*X(88)-JVS(271)*X(98)-JVS(272)*X(100))/(JVS(266))
  X(74) = (X(74)-JVS(262)*X(79)-JVS(263)*X(89)-JVS(264)*X(94)-JVS(265)*X(98))/(JVS(261))
  X(73) = (X(73)-JVS(256)*X(94)-JVS(257)*X(99)-JVS(258)*X(100)-JVS(259)*X(102))/(JVS(255))
  X(72) = (X(72)-JVS(251)*X(86)-JVS(252)*X(96)-JVS(253)*X(98)-JVS(254)*X(102))/(JVS(250))
  X(71) = (X(71)-JVS(245)*X(82)-JVS(246)*X(95)-JVS(247)*X(98)-JVS(248)*X(102)-JVS(249)*X(105))/(JVS(244))
  X(70) = (X(70)-JVS(235)*X(84)-JVS(236)*X(87)-JVS(237)*X(88)-JVS(238)*X(90)-JVS(239)*X(91)-JVS(240)*X(92)-JVS(241)&
            &*X(98)-JVS(242)*X(100)-JVS(243)*X(104))/(JVS(234))
  X(69) = (X(69)-JVS(227)*X(92)-JVS(228)*X(99)-JVS(229)*X(102))/(JVS(226))
  X(68) = (X(68)-JVS(222)*X(94)-JVS(223)*X(98)-JVS(224)*X(101)-JVS(225)*X(102))/(JVS(221))
  X(67) = (X(67)-JVS(218)*X(98)-JVS(219)*X(100))/(JVS(217))
  X(66) = (X(66)-JVS(214)*X(98)-JVS(215)*X(100)-JVS(216)*X(104))/(JVS(213))
  X(65) = (X(65)-JVS(211)*X(98)-JVS(212)*X(100))/(JVS(210))
  X(64) = (X(64)-JVS(207)*X(98)-JVS(208)*X(102)-JVS(209)*X(105))/(JVS(206))
  X(63) = (X(63)-JVS(203)*X(98)-JVS(204)*X(103)-JVS(205)*X(105))/(JVS(202))
  X(62) = (X(62)-JVS(196)*X(96)-JVS(197)*X(97)-JVS(198)*X(98)-JVS(199)*X(102)-JVS(200)*X(103)-JVS(201)*X(104))&
            &/(JVS(195))
  X(61) = (X(61)-JVS(191)*X(97)-JVS(192)*X(98)-JVS(193)*X(102)-JVS(194)*X(103))/(JVS(190))
  X(60) = (X(60)-JVS(188)*X(98)-JVS(189)*X(99))/(JVS(187))
  X(59) = (X(59)-JVS(180)*X(69)-JVS(181)*X(81)-JVS(182)*X(94)-JVS(183)*X(98)-JVS(184)*X(99))/(JVS(179))
  X(58) = (X(58)-JVS(176)*X(98)-JVS(177)*X(99)-JVS(178)*X(105))/(JVS(175))
  X(57) = (X(57)-JVS(172)*X(68)-JVS(173)*X(94)-JVS(174)*X(98))/(JVS(171))
  X(56) = (X(56)-JVS(169)*X(98)-JVS(170)*X(100))/(JVS(168))
  X(55) = (X(55)-JVS(167)*X(98))/(JVS(166))
  X(54) = (X(54)-JVS(164)*X(95)-JVS(165)*X(105))/(JVS(163))
  X(53) = (X(53)-JVS(161)*X(94)-JVS(162)*X(98))/(JVS(160))
  X(52) = (X(52)-JVS(157)*X(53)-JVS(158)*X(98))/(JVS(156))
  X(51) = (X(51)-JVS(155)*X(98))/(JVS(154))
  X(50) = (X(50)-JVS(152)*X(73)-JVS(153)*X(102))/(JVS(151))
  X(49) = (X(49)-JVS(149)*X(97)-JVS(150)*X(105))/(JVS(148))
  X(48) = (X(48)-JVS(147)*X(98))/(JVS(146))
  X(47) = (X(47)-JVS(145)*X(73))/(JVS(144))
  X(46) = (X(46)-JVS(143)*X(98))/(JVS(142))
  X(45) = (X(45)-JVS(140)*X(46)-JVS(141)*X(98))/(JVS(139))
  X(44) = (X(44)-JVS(137)*X(45)-JVS(138)*X(98))/(JVS(136))
  X(43) = (X(43)-JVS(134)*X(44)-JVS(135)*X(98))/(JVS(133))
  X(42) = (X(42)-JVS(132)*X(98))/(JVS(131))
  X(41) = (X(41)-JVS(129)*X(42)-JVS(130)*X(98))/(JVS(128))
  X(40) = (X(40)-JVS(126)*X(41)-JVS(127)*X(98))/(JVS(125))
  X(39) = (X(39)-JVS(123)*X(40)-JVS(124)*X(98))/(JVS(122))
  X(38) = (X(38)-JVS(118)*X(57)-JVS(119)*X(68)-JVS(120)*X(94)-JVS(121)*X(98))/(JVS(117))
  X(37) = (X(37)-JVS(116)*X(98))/(JVS(115))
  X(36) = (X(36)-JVS(113)*X(37)-JVS(114)*X(98))/(JVS(112))
  X(35) = (X(35)-JVS(110)*X(37)-JVS(111)*X(98))/(JVS(109))
  X(34) = (X(34)-JVS(108)*X(98))/(JVS(107))
  X(33) = (X(33)-JVS(105)*X(34)-JVS(106)*X(98))/(JVS(104))
  X(32) = (X(32)-JVS(103)*X(98))/(JVS(102))
  X(31) = (X(31)-JVS(100)*X(32)-JVS(101)*X(98))/(JVS(99))
  X(30) = (X(30)-JVS(97)*X(32)-JVS(98)*X(98))/(JVS(96))
  X(29) = (X(29)-JVS(93)*X(94)-JVS(94)*X(95)-JVS(95)*X(98))/(JVS(92))
  X(28) = (X(28)-JVS(90)*X(29)-JVS(91)*X(95))/(JVS(89))
  X(27) = (X(27)-JVS(87)*X(29)-JVS(88)*X(94))/(JVS(86))
  X(26) = (X(26)-JVS(84)*X(29)-JVS(85)*X(94))/(JVS(83))
  X(25) = (X(25)-JVS(81)*X(29)-JVS(82)*X(98))/(JVS(80))
  X(24) = (X(24)-JVS(78)*X(29)-JVS(79)*X(98))/(JVS(77))
  X(23) = (X(23)-JVS(75)*X(94)-JVS(76)*X(98))/(JVS(74))
  X(22) = (X(22)-JVS(72)*X(23)-JVS(73)*X(94))/(JVS(71))
  X(21) = (X(21)-JVS(69)*X(23)-JVS(70)*X(94))/(JVS(68))
  X(20) = (X(20)-JVS(66)*X(23)-JVS(67)*X(98))/(JVS(65))
  X(19) = (X(19)-JVS(63)*X(23)-JVS(64)*X(98))/(JVS(62))
  X(18) = (X(18)-JVS(61)*X(98))/(JVS(60))
  X(17) = (X(17)-JVS(58)*X(18)-JVS(59)*X(98))/(JVS(57))
  X(16) = (X(16)-JVS(55)*X(18)-JVS(56)*X(98))/(JVS(54))
  X(15) = (X(15)-JVS(53)*X(98))/(JVS(52))
  X(14) = (X(14)-JVS(50)*X(15)-JVS(51)*X(98))/(JVS(49))
  X(13) = (X(13)-JVS(47)*X(15)-JVS(48)*X(98))/(JVS(46))
  X(12) = (X(12)-JVS(45)*X(98))/(JVS(44))
  X(11) = (X(11)-JVS(42)*X(12)-JVS(43)*X(98))/(JVS(41))
  X(10) = (X(10)-JVS(36)*X(83)-JVS(37)*X(94)-JVS(38)*X(95)-JVS(39)*X(98)-JVS(40)*X(101))/(JVS(35))
  X(9) = (X(9)-JVS(33)*X(48)-JVS(34)*X(98))/(JVS(32))
  X(8) = (X(8)-JVS(30)*X(48)-JVS(31)*X(98))/(JVS(29))
  X(7) = (X(7)-JVS(24)*X(88)-JVS(25)*X(94)-JVS(26)*X(95)-JVS(27)*X(98)-JVS(28)*X(101))/(JVS(23))
  X(6) = (X(6)-JVS(18)*X(88)-JVS(19)*X(94)-JVS(20)*X(95)-JVS(21)*X(98)-JVS(22)*X(101))/(JVS(17))
  X(5) = (X(5)-JVS(15)*X(55)-JVS(16)*X(98))/(JVS(14))
  X(4) = (X(4)-JVS(12)*X(55)-JVS(13)*X(98))/(JVS(11))
  X(3) = (X(3)-JVS(8)*X(82)-JVS(9)*X(95)-JVS(10)*X(98))/(JVS(7))
  X(2) = (X(2)-JVS(5)*X(51)-JVS(6)*X(98))/(JVS(4))
  X(1) = (X(1)-JVS(2)*X(51)-JVS(3)*X(98))/(JVS(1))
      
END SUBROUTINE cb05_sorg_vbs_aq_KppSolve
























      SUBROUTINE cb05_sorg_vbs_aq_WCOPY(N,X,incX,Y,incY)








      
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

      END SUBROUTINE cb05_sorg_vbs_aq_WCOPY



      SUBROUTINE cb05_sorg_vbs_aq_WAXPY(N,Alpha,X,incX,Y,incY)









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
      
      END SUBROUTINE cb05_sorg_vbs_aq_WAXPY




      SUBROUTINE cb05_sorg_vbs_aq_WSCAL(N,Alpha,X,incX)









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

      END SUBROUTINE cb05_sorg_vbs_aq_WSCAL


      REAL(kind=dp) FUNCTION cb05_sorg_vbs_aq_WLAMCH( C )








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
          CALL cb05_sorg_vbs_aq_WLAMCH_ADD(ONE,Eps,Sum)
          IF (Sum.LE.ONE) GOTO 10
        END DO
        PRINT*,'ERROR IN WLAMCH. EPS < ',Eps
        RETURN
10      Eps = Eps*2
        i = i-1      
      END IF

      cb05_sorg_vbs_aq_WLAMCH = Eps

      END FUNCTION cb05_sorg_vbs_aq_WLAMCH
     
      SUBROUTINE cb05_sorg_vbs_aq_WLAMCH_ADD( A, B, Sum )

      
      REAL(kind=dp) A, B, Sum
      Sum = A + B

      END SUBROUTINE cb05_sorg_vbs_aq_WLAMCH_ADD




      SUBROUTINE cb05_sorg_vbs_aq_SET2ZERO(N,Y)




      
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

      END SUBROUTINE cb05_sorg_vbs_aq_SET2ZERO



      REAL(kind=dp) FUNCTION cb05_sorg_vbs_aq_WDOT (N, DX, incX, DY, incY) 









      IMPLICIT NONE
      INTEGER :: N, incX, incY
      REAL(kind=dp) :: DX(N), DY(N) 

      INTEGER :: i, IX, IY, M, MP1, NS
                                 
      cb05_sorg_vbs_aq_WDOT = 0.0D0 
      IF (N .LE. 0) RETURN 
      IF (incX .EQ. incY) IF (incX-1) 5,20,60 



    5 IX = 1 
      IY = 1 
      IF (incX .LT. 0) IX = (-N+1)*incX + 1 
      IF (incY .LT. 0) IY = (-N+1)*incY + 1 
      DO i = 1,N 
        cb05_sorg_vbs_aq_WDOT = cb05_sorg_vbs_aq_WDOT + DX(IX)*DY(IY) 
        IX = IX + incX 
        IY = IY + incY 
      END DO 
      RETURN 





   20 M = MOD(N,5) 
      IF (M .EQ. 0) GO TO 40 
      DO i = 1,M 
         cb05_sorg_vbs_aq_WDOT = cb05_sorg_vbs_aq_WDOT + DX(i)*DY(i) 
      END DO 
      IF (N .LT. 5) RETURN 
   40 MP1 = M + 1 
      DO i = MP1,N,5 
          cb05_sorg_vbs_aq_WDOT = cb05_sorg_vbs_aq_WDOT + DX(i)*DY(i) + DX(i+1)*DY(i+1) +&
                   DX(i+2)*DY(i+2) +  &
                   DX(i+3)*DY(i+3) + DX(i+4)*DY(i+4)                   
      END DO 
      RETURN 



   60 NS = N*incX 
      DO i = 1,NS,incX 
        cb05_sorg_vbs_aq_WDOT = cb05_sorg_vbs_aq_WDOT + DX(i)*DY(i) 
      END DO 

      END FUNCTION cb05_sorg_vbs_aq_WDOT                                          




END MODULE cb05_sorg_vbs_aq_Integrator
