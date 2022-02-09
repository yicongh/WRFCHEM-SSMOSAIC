
































MODULE crimech_Integrator

 USE crimech_Parameters
 USE crimech_Precision
 USE crimech_JacobianSP

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

SUBROUTINE  crimech_INTEGRATE( TIN, TOUT, &
  FIX, VAR,  RCONST, ATOL, RTOL, IRR_WRK,  &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U  )

   USE crimech_Parameters

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

   CALL crimech_Rosenbrock(VAR, FIX, RCONST, TIN,TOUT,   &
         ATOL,RTOL,               &
         RCNTRL,ICNTRL,RSTATUS,ISTATUS,IRR_WRK,IERR)

   STEPMIN = RCNTRL(ihexit)
   
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U))    IERR_U       = IERR

END SUBROUTINE  crimech_INTEGRATE


SUBROUTINE  crimech_Rosenbrock(Y, FIX, RCONST, Tstart,Tend, &
           AbsTol,RelTol,            &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IRR_WRK,IERR)







































































































  USE crimech_Parameters

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
      CALL crimech_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
      RETURN
   END IF


   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 100000
   ELSEIF (ICNTRL(4) > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL crimech_ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF


   Roundoff = crimech_WLAMCH('E')


   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      Hmin = RCNTRL(1)
   ELSE
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL crimech_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL crimech_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL crimech_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2_dp
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL crimech_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0_dp
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL crimech_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1_dp
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL crimech_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9_dp
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL crimech_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.0_dp*Roundoff) &
         .OR. (RelTol(i) >= 1.0_dp) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL crimech_ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO



   SELECT CASE (Method)
     CASE (1)
       CALL crimech_Ros2(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (2)
       CALL crimech_Ros3(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (3)
       CALL crimech_Ros4(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (4)
       CALL crimech_Rodas3(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (5)
       CALL crimech_Rodas4(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(4)=', Method
       CALL crimech_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       RETURN
   END SELECT


   CALL crimech_ros_Integrator(Y,Tstart,Tend,Texit,      &
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



 SUBROUTINE  crimech_ros_ErrorMsg(Code,T,H,IERR)



   USE crimech_Precision

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

 END SUBROUTINE  crimech_ros_ErrorMsg


 SUBROUTINE  crimech_ros_Integrator (Y, Tstart, Tend, T,     &
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
      CALL crimech_ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1_dp*H) == T).OR.(H <= Roundoff) ) THEN  
      CALL crimech_ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF


   Hexit = H
   H = MIN(H,ABS(Tend-T))


   CALL crimech_FunTemplate(T,Y,Fcn0, RCONST, FIX, Nfun)
   IF( T == Tstart ) THEN
     CALL crimech_IRRFun( Y, FIX, RCONST, IRR_WRK )
   ENDIF


   IF (.NOT.Autonomous) THEN
      CALL crimech_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )
   END IF


   CALL crimech_JacTemplate(T,Y,Jac0, FIX, Njac, RCONST)


UntilAccepted: DO

   CALL crimech_ros_PrepareMatrix(H,Direction,ros_Gamma(1), &
          Jac0,Ghimj,Pivot,Singular, Ndec,  Nsng )
   IF (Singular) THEN 
       CALL crimech_ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF


Stage: DO istage = 1, ros_S

      
       ioffset = NVAR*(istage-1)

      
       IF ( istage == 1 ) THEN
         CALL crimech_WCOPY(NVAR,Fcn0,1,Fcn,1)
      
       ELSEIF ( ros_NewF(istage) ) THEN
         CALL crimech_WCOPY(NVAR,Y,1,Ynew,1)
         DO j = 1, istage-1
           CALL crimech_WAXPY(NVAR,ros_A((istage-1)*(istage-2)/2+j), &
            K(NVAR*(j-1)+1),1,Ynew,1)
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL crimech_FunTemplate(Tau,Ynew,Fcn, RCONST, FIX, Nfun)
       END IF 
       CALL crimech_WCOPY(NVAR,Fcn,1,K(ioffset+1),1)
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         CALL crimech_WAXPY(NVAR,HC,K(NVAR*(j-1)+1),1,K(ioffset+1),1)
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         CALL crimech_WAXPY(NVAR,HG,dFdT,1,K(ioffset+1),1)
       END IF
       CALL crimech_ros_Solve(Ghimj, Pivot, K(ioffset+1), Nsol)

   END DO Stage



   CALL crimech_WCOPY(NVAR,Y,1,Ynew,1)
   DO j=1,ros_S
         CALL crimech_WAXPY(NVAR,ros_M(j),K(NVAR*(j-1)+1),1,Ynew,1)
   END DO


   CALL crimech_WSCAL(NVAR,ZERO,Yerr,1)
   DO j=1,ros_S
        CALL crimech_WAXPY(NVAR,ros_E(j),K(NVAR*(j-1)+1),1,Yerr,1)
   END DO
   Err = crimech_ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )


   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac


   Nstp = Nstp+1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  
      Nacc = Nacc+1
      CALL crimech_WCOPY(NVAR,Ynew,1,Y,1)
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

  END SUBROUTINE  crimech_ros_Integrator



  REAL(kind=dp) FUNCTION  crimech_ros_ErrorNorm ( Y, Ynew, Yerr, &
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

    crimech_ros_ErrorNorm = Err

  END FUNCTION  crimech_ros_ErrorNorm



  SUBROUTINE crimech_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )



   IMPLICIT NONE


   REAL(kind=dp), INTENT(IN) :: T, Roundoff, Y(NVAR), Fcn0(NVAR)
   REAL(kind=dp), INTENT(IN) :: RCONST(NREACT), FIX(NFIX)

   REAL(kind=dp), INTENT(OUT) :: dFdT(NVAR)

   INTEGER, INTENT(INOUT) ::Nfun

   REAL(kind=dp) :: Delta
   REAL(kind=dp), PARAMETER :: ONE = 1.0_dp, DeltaMin = 1.0E-6_dp

   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL crimech_FunTemplate(T+Delta,Y,dFdT, RCONST, FIX, Nfun)
   CALL crimech_WAXPY(NVAR,(-ONE),Fcn0,1,dFdT,1)
   CALL crimech_WSCAL(NVAR,(ONE/Delta),dFdT,1)

  END SUBROUTINE  crimech_ros_FunTimeDeriv



  SUBROUTINE  crimech_ros_PrepareMatrix ( H, Direction, gam, &
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


     CALL crimech_WCOPY(LU_NONZERO,Jac0,1,Ghimj,1)
     CALL crimech_WSCAL(LU_NONZERO,(-ONE),Ghimj,1)
     ghinv = ONE/(Direction*H*gam)
     DO i=1,NVAR
       Ghimj(LU_DIAG(i)) = Ghimj(LU_DIAG(i))+ghinv
     END DO

     CALL crimech_ros_Decomp( Ghimj, Pivot, ising, Ndec )
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

  END SUBROUTINE  crimech_ros_PrepareMatrix



  SUBROUTINE  crimech_ros_Decomp( A, Pivot, ising, Ndec )



   IMPLICIT NONE

   REAL(kind=dp), INTENT(INOUT) :: A(LU_NONZERO)

   INTEGER, INTENT(OUT) :: Pivot(NVAR), ising
   INTEGER, INTENT(INOUT) :: Ndec 

   CALL crimech_KppDecomp ( A, ising )
   Pivot(1) = 1
   Ndec = Ndec + 1

  END SUBROUTINE  crimech_ros_Decomp



  SUBROUTINE  crimech_ros_Solve( A, Pivot, b, Nsol )



   IMPLICIT NONE

   REAL(kind=dp), INTENT(IN) :: A(LU_NONZERO)
   INTEGER, INTENT(IN) :: Pivot(NVAR)

   INTEGER, INTENT(INOUT) :: nsol 

   REAL(kind=dp), INTENT(INOUT) :: b(NVAR)


   CALL crimech_KppSolve( A, b )

   Nsol = Nsol+1

  END SUBROUTINE  crimech_ros_Solve




  SUBROUTINE  crimech_Ros2 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

 END SUBROUTINE  crimech_Ros2



  SUBROUTINE  crimech_Ros3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  crimech_Ros3





  SUBROUTINE  crimech_Ros4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  crimech_Ros4


  SUBROUTINE  crimech_Rodas3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  crimech_Rodas3


  SUBROUTINE  crimech_Rodas4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  crimech_Rodas4




END SUBROUTINE  crimech_Rosenbrock




SUBROUTINE  crimech_FunTemplate( T, Y, Ydot, RCONST, FIX, Nfun )




   USE crimech_Parameters




   REAL(kind=dp) :: T, Y(NVAR)
   REAL(kind=dp) :: RCONST(NREACT)
   REAL(kind=dp) :: FIX(NFIX)

   REAL(kind=dp) :: Ydot(NVAR)
   INTEGER :: Nfun









   CALL crimech_Fun( Y, FIX, RCONST, Ydot )


   Nfun = Nfun+1

END SUBROUTINE  crimech_FunTemplate



SUBROUTINE  crimech_JacTemplate( T, Y, Jcb, FIX, Njac, RCONST )




 USE crimech_Parameters
 
 USE crimech_Jacobian



    REAL(kind=dp) :: T, Y(NVAR)
    REAL(kind=dp) :: FIX(NFIX)
    REAL(kind=dp) :: RCONST(NREACT)

    INTEGER :: Njac


    REAL(kind=dp) :: Jcb(LU_NONZERO)

    REAL(kind=dp) :: Told





    CALL crimech_Jac_SP( Y, FIX, RCONST, Jcb )


    Njac = Njac+1

END SUBROUTINE  crimech_JacTemplate

















SUBROUTINE crimech_Fun ( V, F, RCT, Vdot )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: Vdot(NVAR)




  REAL(kind=dp) :: A(NREACT)


  A(1) = RCT(1)*V(222)
  A(2) = RCT(2)*V(222)
  A(3) = RCT(3)*V(106)
  A(4) = RCT(4)*V(229)
  A(5) = RCT(5)*V(230)
  A(6) = RCT(6)*V(230)
  A(7) = RCT(7)*V(40)
  A(8) = RCT(8)*V(185)
  A(9) = RCT(9)*V(206)
  A(10) = RCT(10)*V(206)
  A(11) = RCT(11)*V(231)
  A(12) = RCT(12)*V(226)
  A(13) = RCT(13)*V(187)
  A(14) = RCT(14)*V(31)
  A(15) = RCT(15)*V(139)
  A(16) = RCT(16)*V(131)
  A(17) = RCT(17)*V(149)
  A(18) = RCT(18)*V(179)
  A(19) = RCT(19)*V(140)
  A(20) = RCT(20)*V(169)
  A(21) = RCT(21)*V(180)
  A(22) = RCT(22)*V(196)
  A(23) = RCT(23)*V(201)
  A(24) = RCT(24)*V(150)
  A(25) = RCT(25)*V(193)
  A(26) = RCT(26)*V(148)
  A(27) = RCT(27)*V(44)
  A(28) = RCT(28)*V(37)
  A(29) = RCT(29)*V(194)
  A(30) = RCT(30)*V(151)
  A(31) = RCT(31)*V(161)
  A(32) = RCT(32)*V(154)
  A(33) = RCT(33)*V(182)
  A(34) = RCT(34)*V(167)
  A(35) = RCT(35)*V(165)
  A(36) = RCT(36)*V(158)
  A(37) = RCT(37)*V(90)
  A(38) = RCT(38)*V(92)
  A(39) = RCT(39)*V(98)
  A(40) = RCT(40)*V(85)
  A(41) = RCT(41)*V(134)
  A(42) = RCT(42)*V(134)
  A(43) = RCT(43)*V(79)
  A(44) = RCT(44)*V(103)
  A(45) = RCT(45)*V(118)
  A(46) = RCT(46)*V(120)
  A(47) = RCT(47)*V(138)
  A(48) = RCT(48)*V(97)
  A(49) = RCT(49)*V(101)
  A(50) = RCT(50)*V(73)
  A(51) = RCT(51)*V(56)
  A(52) = RCT(52)*V(75)
  A(53) = RCT(53)*V(136)
  A(54) = RCT(54)*V(136)
  A(55) = RCT(55)*V(80)
  A(56) = RCT(56)*V(57)
  A(57) = RCT(57)*V(53)
  A(58) = RCT(58)*V(62)
  A(59) = RCT(59)*V(84)
  A(60) = RCT(60)*V(58)
  A(61) = RCT(61)*V(59)
  A(62) = RCT(62)*V(115)
  A(63) = RCT(63)*V(81)
  A(64) = RCT(64)*V(113)
  A(65) = RCT(65)*V(113)
  A(66) = RCT(66)*V(114)
  A(67) = RCT(67)*V(86)
  A(68) = RCT(68)*V(64)
  A(69) = RCT(69)*V(88)
  A(70) = RCT(70)*V(105)
  A(71) = RCT(71)*V(117)
  A(72) = RCT(72)*V(82)
  A(73) = RCT(73)*V(137)
  A(74) = RCT(74)*V(129)
  A(75) = RCT(75)*V(60)
  A(76) = RCT(76)*V(83)
  A(77) = RCT(77)*V(61)
  A(78) = RCT(78)*V(89)
  A(79) = RCT(79)*V(91)
  A(80) = RCT(80)*V(130)
  A(81) = RCT(81)*V(67)
  A(82) = RCT(82)*V(69)
  A(83) = RCT(83)*V(68)
  A(84) = RCT(84)*V(70)
  A(85) = RCT(85)*V(71)
  A(86) = RCT(86)*V(93)
  A(87) = RCT(87)*V(123)
  A(88) = RCT(88)*V(94)
  A(89) = RCT(89)*V(124)
  A(90) = RCT(90)*V(72)
  A(91) = RCT(91)*V(100)
  A(92) = RCT(92)*V(125)
  A(93) = RCT(93)*V(154)
  A(94) = RCT(94)*V(182)
  A(95) = RCT(95)*V(167)
  A(96) = RCT(96)*V(173)*F(2)
  A(97) = RCT(97)*V(34)
  A(98) = RCT(98)*V(173)*V(222)
  A(99) = RCT(99)*V(173)*V(227)
  A(100) = RCT(100)*V(173)*V(229)
  A(101) = RCT(101)*V(173)*V(229)
  A(102) = RCT(102)*V(7)*F(2)
  A(103) = RCT(103)*V(222)*V(227)
  A(104) = RCT(104)*V(222)*V(229)
  A(105) = RCT(105)*V(227)*V(227)*F(2)
  A(106) = RCT(106)*V(227)*V(230)
  A(107) = RCT(107)*V(229)*V(230)
  A(108) = RCT(108)*V(229)*V(230)
  A(109) = RCT(109)*V(32)
  A(110) = RCT(110)*V(7)*F(1)
  A(111) = RCT(111)*V(222)*V(228)
  A(112) = RCT(112)*V(228)*F(2)
  A(113) = RCT(113)*V(189)*V(228)
  A(114) = RCT(114)*V(106)*V(228)
  A(115) = RCT(115)*V(222)*V(232)
  A(116) = RCT(116)*V(228)*V(232)
  A(117) = RCT(117)*V(232)*V(232)
  A(118) = RCT(118)*V(232)*V(232)*F(1)
  A(119) = RCT(119)*V(227)*V(228)
  A(120) = RCT(120)*V(228)*V(229)
  A(121) = RCT(121)*V(228)*V(230)
  A(122) = RCT(122)*V(227)*V(232)
  A(123) = RCT(123)*V(229)*V(232)
  A(124) = RCT(124)*V(74)
  A(125) = RCT(125)*V(74)*V(228)
  A(126) = RCT(126)*V(230)*V(232)
  A(127) = RCT(127)*V(40)*V(228)
  A(128) = RCT(128)*V(185)*V(228)
  A(129) = RCT(129)*V(127)*V(173)
  A(130) = RCT(130)*V(127)*V(228)
  A(131) = RCT(131)*V(9)*F(2)
  A(132) = RCT(132)*V(4)*F(1)*F(1)
  A(133) = RCT(133)*V(19)*V(228)
  A(134) = RCT(134)*V(8)*V(228)
  A(135) = RCT(135)*V(25)*V(228)
  A(136) = RCT(136)*V(25)*V(228)
  A(137) = RCT(137)*V(10)*V(228)
  A(138) = RCT(138)*V(153)*V(228)
  A(139) = RCT(139)*V(143)*V(228)
  A(140) = RCT(140)*V(162)*V(228)
  A(141) = RCT(141)*V(153)*V(230)
  A(142) = RCT(142)*V(143)*V(230)
  A(143) = RCT(143)*V(162)*V(230)
  A(144) = RCT(144)*V(153)*V(222)
  A(145) = RCT(145)*V(153)*V(222)
  A(146) = RCT(146)*V(143)*V(222)
  A(147) = RCT(147)*V(143)*V(222)
  A(148) = RCT(148)*V(162)*V(222)
  A(149) = RCT(149)*V(162)*V(222)
  A(150) = RCT(150)*V(175)*V(228)
  A(151) = RCT(151)*V(175)*V(230)
  A(152) = RCT(152)*V(175)*V(222)
  A(153) = RCT(153)*V(175)*V(222)
  A(154) = RCT(154)*V(164)*V(228)
  A(155) = RCT(155)*V(164)*V(230)
  A(156) = RCT(156)*V(164)*V(222)
  A(157) = RCT(157)*V(164)*V(222)
  A(158) = RCT(158)*V(164)*V(222)
  A(159) = RCT(159)*V(159)*V(228)
  A(160) = RCT(160)*V(159)*V(230)
  A(161) = RCT(161)*V(159)*V(222)
  A(162) = RCT(162)*V(159)*V(222)
  A(163) = RCT(163)*V(159)*V(222)
  A(164) = RCT(164)*V(159)*V(222)
  A(165) = RCT(165)*V(35)*V(228)
  A(166) = RCT(166)*V(35)*V(228)
  A(167) = RCT(167)*V(29)*V(228)
  A(168) = RCT(168)*V(29)*V(228)
  A(169) = RCT(169)*V(30)*V(228)
  A(170) = RCT(170)*V(30)*V(228)
  A(171) = RCT(171)*V(26)*V(228)
  A(172) = RCT(172)*V(26)*V(228)
  A(173) = RCT(173)*V(206)*V(228)
  A(174) = RCT(174)*V(228)*V(231)
  A(175) = RCT(175)*V(226)*V(228)
  A(176) = RCT(176)*V(206)*V(230)
  A(177) = RCT(177)*V(230)*V(231)
  A(178) = RCT(178)*V(226)*V(230)
  A(179) = RCT(179)*V(187)*V(228)
  A(180) = RCT(180)*V(31)*V(228)
  A(181) = RCT(181)*V(36)*V(228)
  A(182) = RCT(182)*V(42)*V(228)
  A(183) = RCT(183)*V(42)*V(228)
  A(184) = RCT(184)*V(50)*V(228)
  A(185) = RCT(185)*V(50)*V(228)
  A(186) = RCT(186)*V(43)*V(228)
  A(187) = RCT(187)*V(43)*V(228)
  A(188) = RCT(188)*V(45)*V(228)
  A(189) = RCT(189)*V(51)*V(228)
  A(190) = RCT(190)*V(211)*V(227)
  A(191) = RCT(191)*V(209)*V(227)
  A(192) = RCT(192)*V(215)*V(227)
  A(193) = RCT(193)*V(156)*V(227)
  A(194) = RCT(194)*V(225)*V(227)
  A(195) = RCT(195)*V(225)*V(227)
  A(196) = RCT(196)*V(223)*V(227)
  A(197) = RCT(197)*V(205)*V(227)
  A(198) = RCT(198)*V(176)*V(227)
  A(199) = RCT(199)*V(177)*V(227)
  A(200) = RCT(200)*V(163)*V(227)
  A(201) = RCT(201)*V(181)*V(227)
  A(202) = RCT(202)*V(181)*V(227)
  A(203) = RCT(203)*V(192)*V(227)
  A(204) = RCT(204)*V(202)*V(227)
  A(205) = RCT(205)*V(199)*V(227)
  A(206) = RCT(206)*V(199)*V(227)
  A(207) = RCT(207)*V(217)*V(227)
  A(208) = RCT(208)*V(213)*V(227)
  A(209) = RCT(209)*V(207)*V(227)
  A(210) = RCT(210)*V(170)*V(227)
  A(211) = RCT(211)*V(224)*V(227)
  A(212) = RCT(212)*V(195)*V(227)
  A(213) = RCT(213)*V(219)*V(227)
  A(214) = RCT(214)*V(203)*V(227)
  A(215) = RCT(215)*V(191)*V(227)
  A(216) = RCT(216)*V(220)*V(227)
  A(217) = RCT(217)*V(208)*V(227)
  A(218) = RCT(218)*V(216)*V(227)
  A(219) = RCT(219)*V(218)*V(227)
  A(220) = RCT(220)*V(178)*V(227)
  A(221) = RCT(221)*V(178)*V(227)
  A(222) = RCT(222)*V(198)*V(227)
  A(223) = RCT(223)*V(198)*V(227)
  A(224) = RCT(224)*V(197)*V(227)
  A(225) = RCT(225)*V(197)*V(227)
  A(226) = RCT(226)*V(197)*V(227)
  A(227) = RCT(227)*V(147)*V(227)
  A(228) = RCT(228)*V(166)*V(227)
  A(229) = RCT(229)*V(145)*V(227)
  A(230) = RCT(230)*V(144)*V(227)
  A(231) = RCT(231)*V(188)*V(227)
  A(232) = RCT(232)*V(186)*V(227)
  A(233) = RCT(233)*V(186)*V(227)
  A(234) = RCT(234)*V(146)*V(227)
  A(235) = RCT(235)*V(221)*V(227)
  A(236) = RCT(236)*V(212)*V(227)
  A(237) = RCT(237)*V(157)*V(227)
  A(238) = RCT(238)*V(204)*V(227)
  A(239) = RCT(239)*V(190)*V(227)
  A(240) = RCT(240)*V(183)*V(227)
  A(241) = RCT(241)*V(184)*V(227)
  A(242) = RCT(242)*V(184)*V(227)
  A(243) = RCT(243)*V(168)*V(227)
  A(244) = RCT(244)*V(214)*V(227)
  A(245) = RCT(245)*V(214)*V(227)
  A(246) = RCT(246)*V(210)*V(227)
  A(247) = RCT(247)*V(211)*V(227)
  A(248) = RCT(248)*V(209)*V(227)
  A(249) = RCT(249)*V(215)*V(227)
  A(250) = RCT(250)*V(156)*V(227)
  A(251) = RCT(251)*V(225)*V(227)
  A(252) = RCT(252)*V(223)*V(227)
  A(253) = RCT(253)*V(205)*V(227)
  A(254) = RCT(254)*V(199)*V(227)
  A(255) = RCT(255)*V(217)*V(227)
  A(256) = RCT(256)*V(213)*V(227)
  A(257) = RCT(257)*V(207)*V(227)
  A(258) = RCT(258)*V(170)*V(227)
  A(259) = RCT(259)*V(224)*V(227)
  A(260) = RCT(260)*V(195)*V(227)
  A(261) = RCT(261)*V(178)*V(227)
  A(262) = RCT(262)*V(163)*V(227)
  A(263) = RCT(263)*V(181)*V(227)
  A(264) = RCT(264)*V(192)*V(227)
  A(265) = RCT(265)*V(202)*V(227)
  A(266) = RCT(266)*V(186)*V(227)
  A(267) = RCT(267)*V(212)*V(227)
  A(268) = RCT(268)*V(184)*V(227)
  A(269) = RCT(269)*V(214)*V(227)
  A(270) = RCT(270)*V(210)*V(227)
  A(271) = RCT(271)*V(211)*V(230)
  A(272) = RCT(272)*V(209)*V(230)
  A(273) = RCT(273)*V(215)*V(230)
  A(274) = RCT(274)*V(156)*V(230)
  A(275) = RCT(275)*V(225)*V(230)
  A(276) = RCT(276)*V(225)*V(230)
  A(277) = RCT(277)*V(223)*V(230)
  A(278) = RCT(278)*V(205)*V(230)
  A(279) = RCT(279)*V(176)*V(230)
  A(280) = RCT(280)*V(177)*V(230)
  A(281) = RCT(281)*V(163)*V(230)
  A(282) = RCT(282)*V(181)*V(230)
  A(283) = RCT(283)*V(181)*V(230)
  A(284) = RCT(284)*V(192)*V(230)
  A(285) = RCT(285)*V(202)*V(230)
  A(286) = RCT(286)*V(199)*V(230)
  A(287) = RCT(287)*V(199)*V(230)
  A(288) = RCT(288)*V(217)*V(230)
  A(289) = RCT(289)*V(213)*V(230)
  A(290) = RCT(290)*V(207)*V(230)
  A(291) = RCT(291)*V(170)*V(230)
  A(292) = RCT(292)*V(224)*V(230)
  A(293) = RCT(293)*V(195)*V(230)
  A(294) = RCT(294)*V(219)*V(230)
  A(295) = RCT(295)*V(203)*V(230)
  A(296) = RCT(296)*V(191)*V(230)
  A(297) = RCT(297)*V(220)*V(230)
  A(298) = RCT(298)*V(208)*V(230)
  A(299) = RCT(299)*V(216)*V(230)
  A(300) = RCT(300)*V(218)*V(230)
  A(301) = RCT(301)*V(178)*V(230)
  A(302) = RCT(302)*V(178)*V(230)
  A(303) = RCT(303)*V(198)*V(230)
  A(304) = RCT(304)*V(198)*V(230)
  A(305) = RCT(305)*V(197)*V(230)
  A(306) = RCT(306)*V(197)*V(230)
  A(307) = RCT(307)*V(197)*V(230)
  A(308) = RCT(308)*V(147)*V(230)
  A(309) = RCT(309)*V(166)*V(230)
  A(310) = RCT(310)*V(145)*V(230)
  A(311) = RCT(311)*V(144)*V(230)
  A(312) = RCT(312)*V(188)*V(230)
  A(313) = RCT(313)*V(186)*V(230)
  A(314) = RCT(314)*V(146)*V(230)
  A(315) = RCT(315)*V(221)*V(230)
  A(316) = RCT(316)*V(212)*V(230)
  A(317) = RCT(317)*V(157)*V(230)
  A(318) = RCT(318)*V(204)*V(230)
  A(319) = RCT(319)*V(190)*V(230)
  A(320) = RCT(320)*V(183)*V(230)
  A(321) = RCT(321)*V(184)*V(230)
  A(322) = RCT(322)*V(214)*V(230)
  A(323) = RCT(323)*V(210)*V(230)
  A(324) = RCT(324)*V(168)*V(230)
  A(325) = RCT(325)*V(211)*V(232)
  A(326) = RCT(326)*V(209)*V(232)
  A(327) = RCT(327)*V(215)*V(232)
  A(328) = RCT(328)*V(156)*V(232)
  A(329) = RCT(329)*V(225)*V(232)
  A(330) = RCT(330)*V(223)*V(232)
  A(331) = RCT(331)*V(205)*V(232)
  A(332) = RCT(332)*V(176)*V(232)
  A(333) = RCT(333)*V(177)*V(232)
  A(334) = RCT(334)*V(163)*V(232)
  A(335) = RCT(335)*V(181)*V(232)
  A(336) = RCT(336)*V(192)*V(232)
  A(337) = RCT(337)*V(202)*V(232)
  A(338) = RCT(338)*V(199)*V(232)
  A(339) = RCT(339)*V(217)*V(232)
  A(340) = RCT(340)*V(213)*V(232)
  A(341) = RCT(341)*V(207)*V(232)
  A(342) = RCT(342)*V(170)*V(232)
  A(343) = RCT(343)*V(224)*V(232)
  A(344) = RCT(344)*V(195)*V(232)
  A(345) = RCT(345)*V(219)*V(232)
  A(346) = RCT(346)*V(203)*V(232)
  A(347) = RCT(347)*V(191)*V(232)
  A(348) = RCT(348)*V(220)*V(232)
  A(349) = RCT(349)*V(208)*V(232)
  A(350) = RCT(350)*V(216)*V(232)
  A(351) = RCT(351)*V(218)*V(232)
  A(352) = RCT(352)*V(178)*V(232)
  A(353) = RCT(353)*V(198)*V(232)
  A(354) = RCT(354)*V(197)*V(232)
  A(355) = RCT(355)*V(147)*V(232)
  A(356) = RCT(356)*V(166)*V(232)
  A(357) = RCT(357)*V(145)*V(232)
  A(358) = RCT(358)*V(144)*V(232)
  A(359) = RCT(359)*V(188)*V(232)
  A(360) = RCT(360)*V(186)*V(232)
  A(361) = RCT(361)*V(146)*V(232)
  A(362) = RCT(362)*V(221)*V(232)
  A(363) = RCT(363)*V(212)*V(232)
  A(364) = RCT(364)*V(157)*V(232)
  A(365) = RCT(365)*V(204)*V(232)
  A(366) = RCT(366)*V(190)*V(232)
  A(367) = RCT(367)*V(183)*V(232)
  A(368) = RCT(368)*V(184)*V(232)
  A(369) = RCT(369)*V(214)*V(232)
  A(370) = RCT(370)*V(210)*V(232)
  A(371) = RCT(371)*V(168)*V(232)
  A(372) = RCT(372)*V(211)
  A(373) = RCT(373)*V(211)
  A(374) = RCT(374)*V(211)
  A(375) = RCT(375)*V(209)
  A(376) = RCT(376)*V(209)
  A(377) = RCT(377)*V(209)
  A(378) = RCT(378)*V(215)
  A(379) = RCT(379)*V(215)
  A(380) = RCT(380)*V(215)
  A(381) = RCT(381)*V(156)
  A(382) = RCT(382)*V(156)
  A(383) = RCT(383)*V(156)
  A(384) = RCT(384)*V(225)
  A(385) = RCT(385)*V(225)
  A(386) = RCT(386)*V(176)
  A(387) = RCT(387)*V(177)
  A(388) = RCT(388)*V(163)
  A(389) = RCT(389)*V(181)
  A(390) = RCT(390)*V(181)
  A(391) = RCT(391)*V(192)
  A(392) = RCT(392)*V(202)
  A(393) = RCT(393)*V(223)
  A(394) = RCT(394)*V(205)
  A(395) = RCT(395)*V(199)
  A(396) = RCT(396)*V(199)
  A(397) = RCT(397)*V(217)
  A(398) = RCT(398)*V(213)
  A(399) = RCT(399)*V(207)
  A(400) = RCT(400)*V(170)
  A(401) = RCT(401)*V(224)
  A(402) = RCT(402)*V(195)
  A(403) = RCT(403)*V(219)
  A(404) = RCT(404)*V(203)
  A(405) = RCT(405)*V(191)
  A(406) = RCT(406)*V(220)
  A(407) = RCT(407)*V(208)
  A(408) = RCT(408)*V(216)
  A(409) = RCT(409)*V(218)
  A(410) = RCT(410)*V(178)
  A(411) = RCT(411)*V(178)
  A(412) = RCT(412)*V(198)
  A(413) = RCT(413)*V(198)
  A(414) = RCT(414)*V(197)
  A(415) = RCT(415)*V(197)
  A(416) = RCT(416)*V(197)
  A(417) = RCT(417)*V(147)
  A(418) = RCT(418)*V(166)
  A(419) = RCT(419)*V(145)
  A(420) = RCT(420)*V(144)
  A(421) = RCT(421)*V(188)
  A(422) = RCT(422)*V(186)
  A(423) = RCT(423)*V(146)
  A(424) = RCT(424)*V(221)
  A(425) = RCT(425)*V(212)
  A(426) = RCT(426)*V(157)
  A(427) = RCT(427)*V(204)
  A(428) = RCT(428)*V(190)
  A(429) = RCT(429)*V(183)
  A(430) = RCT(430)*V(184)
  A(431) = RCT(431)*V(214)
  A(432) = RCT(432)*V(210)
  A(433) = RCT(433)*V(168)
  A(434) = RCT(434)*V(139)*V(228)
  A(435) = RCT(435)*V(131)*V(228)
  A(436) = RCT(436)*V(149)*V(228)
  A(437) = RCT(437)*V(179)*V(228)
  A(438) = RCT(438)*V(140)*V(228)
  A(439) = RCT(439)*V(169)*V(228)
  A(440) = RCT(440)*V(180)*V(228)
  A(441) = RCT(441)*V(201)*V(228)
  A(442) = RCT(442)*V(201)*V(230)
  A(443) = RCT(443)*V(201)*V(222)
  A(444) = RCT(444)*V(201)*V(222)
  A(445) = RCT(445)*V(196)*V(228)
  A(446) = RCT(446)*V(196)*V(230)
  A(447) = RCT(447)*V(150)*V(228)
  A(448) = RCT(448)*V(193)*V(228)
  A(449) = RCT(449)*V(148)*V(228)
  A(450) = RCT(450)*V(44)*V(228)
  A(451) = RCT(451)*V(37)*V(228)
  A(452) = RCT(452)*V(38)*V(228)
  A(453) = RCT(453)*V(194)*V(228)
  A(454) = RCT(454)*V(194)*V(230)
  A(455) = RCT(455)*V(194)*V(222)
  A(456) = RCT(456)*V(194)*V(222)
  A(457) = RCT(457)*V(151)*V(228)
  A(458) = RCT(458)*V(161)*V(228)
  A(459) = RCT(459)*V(154)*V(228)
  A(460) = RCT(460)*V(182)*V(228)
  A(461) = RCT(461)*V(167)*V(228)
  A(462) = RCT(462)*V(165)*V(228)
  A(463) = RCT(463)*V(33)*V(228)
  A(464) = RCT(464)*V(158)*V(228)
  A(465) = RCT(465)*V(165)*V(230)
  A(466) = RCT(466)*V(158)*V(230)
  A(467) = RCT(467)*V(39)*V(228)
  A(468) = RCT(468)*V(132)*V(228)
  A(469) = RCT(469)*V(107)*V(228)
  A(470) = RCT(470)*V(90)*V(228)
  A(471) = RCT(471)*V(92)*V(228)
  A(472) = RCT(472)*V(98)*V(228)
  A(473) = RCT(473)*V(85)*V(228)
  A(474) = RCT(474)*V(134)*V(228)
  A(475) = RCT(475)*V(79)*V(228)
  A(476) = RCT(476)*V(103)*V(228)
  A(477) = RCT(477)*V(65)*V(228)
  A(478) = RCT(478)*V(54)*V(228)
  A(479) = RCT(479)*V(55)*V(228)
  A(480) = RCT(480)*V(102)*V(228)
  A(481) = RCT(481)*V(78)*V(228)
  A(482) = RCT(482)*V(63)*V(228)
  A(483) = RCT(483)*V(118)*V(228)
  A(484) = RCT(484)*V(120)*V(228)
  A(485) = RCT(485)*V(138)*V(228)
  A(486) = RCT(486)*V(66)*V(228)
  A(487) = RCT(487)*V(95)*V(228)
  A(488) = RCT(488)*V(96)*V(228)
  A(489) = RCT(489)*V(97)*V(228)
  A(490) = RCT(490)*V(99)*V(228)
  A(491) = RCT(491)*V(46)*V(228)
  A(492) = RCT(492)*V(46)*V(230)
  A(493) = RCT(493)*V(47)*V(229)
  A(494) = RCT(494)*V(119)*V(228)
  A(495) = RCT(495)*V(119)*V(230)
  A(496) = RCT(496)*V(48)*V(228)
  A(497) = RCT(497)*V(48)*V(230)
  A(498) = RCT(498)*V(49)*V(229)
  A(499) = RCT(499)*V(121)*V(228)
  A(500) = RCT(500)*V(121)*V(230)
  A(501) = RCT(501)*V(101)*V(228)
  A(502) = RCT(502)*V(101)*V(228)
  A(503) = RCT(503)*V(73)*V(228)
  A(504) = RCT(504)*V(56)*V(228)
  A(505) = RCT(505)*V(75)*V(228)
  A(506) = RCT(506)*V(136)*V(228)
  A(507) = RCT(507)*V(80)*V(228)
  A(508) = RCT(508)*V(57)*V(228)
  A(509) = RCT(509)*V(53)*V(228)
  A(510) = RCT(510)*V(62)*V(228)
  A(511) = RCT(511)*V(84)*V(228)
  A(512) = RCT(512)*V(58)*V(228)
  A(513) = RCT(513)*V(59)*V(228)
  A(514) = RCT(514)*V(115)*V(228)
  A(515) = RCT(515)*V(81)*V(228)
  A(516) = RCT(516)*V(113)*V(228)
  A(517) = RCT(517)*V(114)*V(228)
  A(518) = RCT(518)*V(86)*V(228)
  A(519) = RCT(519)*V(64)*V(228)
  A(520) = RCT(520)*V(88)*V(228)
  A(521) = RCT(521)*V(105)*V(228)
  A(522) = RCT(522)*V(117)*V(228)
  A(523) = RCT(523)*V(82)*V(228)
  A(524) = RCT(524)*V(137)*V(228)
  A(525) = RCT(525)*V(129)*V(228)
  A(526) = RCT(526)*V(60)*V(228)
  A(527) = RCT(527)*V(83)*V(228)
  A(528) = RCT(528)*V(61)*V(228)
  A(529) = RCT(529)*V(89)*V(228)
  A(530) = RCT(530)*V(91)*V(228)
  A(531) = RCT(531)*V(130)*V(228)
  A(532) = RCT(532)*V(67)*V(228)
  A(533) = RCT(533)*V(68)*V(228)
  A(534) = RCT(534)*V(69)*V(228)
  A(535) = RCT(535)*V(70)*V(228)
  A(536) = RCT(536)*V(71)*V(228)
  A(537) = RCT(537)*V(93)*V(228)
  A(538) = RCT(538)*V(123)*V(228)
  A(539) = RCT(539)*V(94)*V(228)
  A(540) = RCT(540)*V(124)*V(228)
  A(541) = RCT(541)*V(72)*V(228)
  A(542) = RCT(542)*V(100)*V(228)
  A(543) = RCT(543)*V(125)*V(228)
  A(544) = RCT(544)*V(219)*V(229)
  A(545) = RCT(545)*V(126)
  A(546) = RCT(546)*V(203)*V(229)
  A(547) = RCT(547)*V(111)
  A(548) = RCT(548)*V(191)*V(229)
  A(549) = RCT(549)*V(112)
  A(550) = RCT(550)*V(126)*V(228)
  A(551) = RCT(551)*V(111)*V(228)
  A(552) = RCT(552)*V(112)*V(228)
  A(553) = RCT(553)*V(198)*V(229)
  A(554) = RCT(554)*V(87)
  A(555) = RCT(555)*V(197)*V(229)
  A(556) = RCT(556)*V(116)
  A(557) = RCT(557)*V(116)*V(228)
  A(558) = RCT(558)*V(87)*V(228)
  A(559) = RCT(559)*V(221)*V(229)
  A(560) = RCT(560)*V(122)
  A(561) = RCT(561)*V(122)*V(228)
  A(562) = RCT(562)*V(76)*V(228)
  A(563) = RCT(563)*V(154)*V(228)
  A(564) = RCT(564)*V(182)*V(228)
  A(565) = RCT(565)*V(167)*V(228)
  A(566) = RCT(566)*V(11)*V(228)
  A(567) = RCT(567)*V(12)*V(228)
  A(568) = RCT(568)*V(13)*V(228)
  A(569) = RCT(569)*V(14)*V(228)
  A(570) = RCT(570)*V(18)*V(228)
  A(571) = RCT(571)*V(17)*V(228)
  A(572) = RCT(572)*V(15)*V(228)
  A(573) = RCT(573)*V(16)*V(228)
  A(574) = RCT(574)*V(211)*V(229)
  A(575) = RCT(575)*V(77)
  A(576) = RCT(576)*V(135)
  A(577) = RCT(577)*V(141)
  A(578) = RCT(578)*V(133)
  A(579) = RCT(579)*V(128)
  A(580) = RCT(580)*V(104)
  A(581) = RCT(581)*V(152)
  A(582) = RCT(582)*V(152)
  A(583) = RCT(583)*V(204)*V(227)
  A(584) = RCT(584)*V(160)*V(228)
  A(585) = RCT(585)*V(77)*V(228)
  A(586) = RCT(586)*V(135)*V(228)
  A(587) = RCT(587)*V(135)*V(230)
  A(588) = RCT(588)*V(27)*V(228)
  A(589) = RCT(589)*V(27)*V(228)
  A(590) = RCT(590)*V(28)*V(228)
  A(591) = RCT(591)*V(28)*V(228)
  A(592) = RCT(592)*V(20)*V(228)
  A(593) = RCT(593)*V(21)*V(228)
  A(594) = RCT(594)*V(22)*V(228)
  A(595) = RCT(595)*V(23)*V(228)
  A(596) = RCT(596)*V(155)*V(227)
  A(597) = RCT(597)*V(200)*V(227)
  A(598) = RCT(598)*V(155)*V(227)
  A(599) = RCT(599)*V(200)*V(227)
  A(600) = RCT(600)*V(155)*V(230)
  A(601) = RCT(601)*V(200)*V(230)
  A(602) = RCT(602)*V(155)*V(232)
  A(603) = RCT(603)*V(200)*V(232)
  A(604) = RCT(604)*V(155)
  A(605) = RCT(605)*V(200)
  A(606) = RCT(606)*V(141)*V(228)
  A(607) = RCT(607)*V(133)*V(228)
  A(608) = RCT(608)*V(24)*V(228)
  A(609) = RCT(609)*V(174)*V(227)
  A(610) = RCT(610)*V(174)*V(227)
  A(611) = RCT(611)*V(174)*V(230)
  A(612) = RCT(612)*V(174)*V(232)
  A(613) = RCT(613)*V(174)
  A(614) = RCT(614)*V(128)*V(228)
  A(615) = RCT(615)*V(104)*V(228)
  A(616) = RCT(616)*V(152)*V(228)
  A(617) = RCT(617)*V(152)*V(228)
  A(620) = RCT(620)*V(108)*V(228)
  A(621) = RCT(621)*V(108)*V(228)
  A(622) = RCT(622)*V(108)*V(230)
  A(623) = RCT(623)*V(110)*V(227)
  A(624) = RCT(624)*V(110)*V(110)
  A(625) = RCT(625)*V(142)*V(222)
  A(626) = RCT(626)*V(142)*V(229)
  A(627) = RCT(627)*V(172)*V(229)
  A(628) = RCT(628)*V(172)*V(222)
  A(629) = RCT(629)*V(171)
  A(630) = RCT(630)*V(171)*V(229)
  A(631) = RCT(631)*V(171)*V(222)
  A(632) = RCT(632)*V(52)*V(232)
  A(633) = RCT(633)*V(52)
  A(634) = RCT(634)*V(41)*V(228)
  A(635) = RCT(635)*V(109)*V(228)
  A(636) = RCT(636)*V(109)*V(230)


  Vdot(1) = 0
  Vdot(2) = 0
  Vdot(3) = A(132)+A(633)
  Vdot(4) = A(129)+A(131)-A(132)
  Vdot(5) = A(632)+0.05*A(635)
  Vdot(6) = 0.05*A(634)
  Vdot(7) = A(1)-A(102)-A(110)
  Vdot(8) = -A(134)
  Vdot(9) = A(130)-A(131)
  Vdot(10) = -A(137)
  Vdot(11) = -A(566)
  Vdot(12) = -A(567)
  Vdot(13) = -A(568)
  Vdot(14) = -A(569)
  Vdot(15) = -A(572)
  Vdot(16) = -A(573)
  Vdot(17) = -A(571)
  Vdot(18) = -A(570)
  Vdot(19) = -A(133)
  Vdot(20) = -A(592)
  Vdot(21) = -A(593)
  Vdot(22) = -A(594)
  Vdot(23) = -A(595)
  Vdot(24) = -A(608)
  Vdot(25) = -A(135)-A(136)
  Vdot(26) = -A(171)-A(172)
  Vdot(27) = -A(588)-A(589)
  Vdot(28) = -A(590)-A(591)
  Vdot(29) = -A(167)-A(168)
  Vdot(30) = -A(169)-A(170)
  Vdot(31) = -A(14)-A(180)
  Vdot(32) = A(108)-A(109)
  Vdot(33) = -A(463)+A(487)
  Vdot(34) = -A(97)+A(574)
  Vdot(35) = -A(165)-A(166)
  Vdot(36) = -A(181)+A(374)
  Vdot(37) = -A(28)-A(451)+A(515)
  Vdot(38) = -A(452)+A(490)+A(542)
  Vdot(39) = A(158)-A(467)
  Vdot(40) = -A(7)+A(119)-A(127)
  Vdot(41) = A(621)-A(634)
  Vdot(42) = -A(182)-A(183)+A(377)
  Vdot(43) = -A(186)-A(187)+A(383)
  Vdot(44) = -A(27)-A(450)+A(514)
  Vdot(45) = A(145)+A(153)+A(165)-A(188)
  Vdot(46) = A(168)-A(491)-A(492)
  Vdot(47) = A(491)+A(492)-A(493)
  Vdot(48) = A(170)-A(496)-A(497)
  Vdot(49) = A(496)+A(497)-A(498)
  Vdot(50) = -A(184)-A(185)+A(380)
  Vdot(51) = A(147)+A(149)-A(189)
  Vdot(52) = A(630)+A(631)-A(632)-A(633)
  Vdot(53) = -A(57)+A(345)-A(509)
  Vdot(54) = A(255)-A(478)
  Vdot(55) = A(256)-A(479)
  Vdot(56) = -A(51)+A(327)-A(504)
  Vdot(57) = -A(56)+A(331)-A(508)
  Vdot(58) = -A(60)+A(348)-A(512)
  Vdot(59) = -A(61)+A(349)-A(513)
  Vdot(60) = -A(75)+A(355)-A(526)
  Vdot(61) = -A(77)+A(357)-A(528)
  Vdot(62) = -A(58)+A(346)-A(510)
  Vdot(63) = A(261)-A(482)
  Vdot(64) = -A(68)+A(358)-A(519)
  Vdot(65) = A(254)-A(477)
  Vdot(66) = A(266)-A(486)
  Vdot(67) = -A(81)+A(360)-A(532)
  Vdot(68) = -A(83)+A(362)-A(533)
  Vdot(69) = -A(82)+A(361)-A(534)
  Vdot(70) = -A(84)+A(363)-A(535)
  Vdot(71) = -A(85)+A(364)-A(536)
  Vdot(72) = -A(90)+A(369)-A(541)
  Vdot(73) = -A(50)+A(326)-A(503)
  Vdot(74) = A(123)-A(124)-A(125)
  Vdot(75) = -A(52)+A(328)-A(505)
  Vdot(76) = A(93)+A(94)+A(95)-A(562)+A(563)+A(564)+A(565)+A(581)+A(617)
  Vdot(77) = -A(575)+A(584)-A(585)
  Vdot(78) = A(258)+A(260)-A(481)
  Vdot(79) = -A(43)+A(252)-A(475)
  Vdot(80) = -A(55)+A(330)+A(333)-A(507)
  Vdot(81) = -A(63)+A(351)-A(515)
  Vdot(82) = -A(72)+A(340)-A(523)
  Vdot(83) = -A(76)+A(356)-A(527)
  Vdot(84) = -A(59)+A(347)-A(511)
  Vdot(85) = -A(40)+A(250)-A(473)
  Vdot(86) = -A(67)+A(354)-A(518)
  Vdot(87) = A(553)-A(554)-A(558)
  Vdot(88) = -A(69)+A(359)-A(520)
  Vdot(89) = -A(78)+A(334)-A(529)
  Vdot(90) = -A(37)+A(247)-A(470)
  Vdot(91) = -A(79)+A(335)-A(530)
  Vdot(92) = -A(38)+A(248)-A(471)
  Vdot(93) = -A(86)+A(365)-A(537)
  Vdot(94) = -A(88)+A(367)-A(539)
  Vdot(95) = A(267)-A(487)
  Vdot(96) = A(268)-A(488)
  Vdot(97) = -A(48)+A(269)-A(489)
  Vdot(98) = -A(39)+A(249)-A(472)
  Vdot(99) = A(270)-A(490)
  Vdot(100) = -A(91)+A(370)-A(542)
  Vdot(101) = -A(49)+A(325)-A(501)-A(502)
  Vdot(102) = A(257)+A(259)-A(480)
  Vdot(103) = -A(44)+A(253)-A(476)
  Vdot(104) = -A(580)+A(612)-A(615)
  Vdot(105) = -A(70)+A(338)-A(521)
  Vdot(106) = -A(3)-A(114)+A(117)+A(118)+A(157)+A(162)+A(444)+A(456)
  Vdot(107) = A(48)+A(90)+A(163)+A(244)+A(322)+A(431)-A(469)+A(489)+A(541)
  Vdot(108) = -A(620)-A(621)-A(622)
  Vdot(109) = 0.95*A(634)-A(635)-A(636)
  Vdot(110) = A(620)+A(622)-A(623)-2*A(624)
  Vdot(111) = A(546)-A(547)-A(551)
  Vdot(112) = A(548)-A(549)-A(552)
  Vdot(113) = -A(64)-A(65)+A(352)-A(516)
  Vdot(114) = -A(66)+A(353)-A(517)
  Vdot(115) = -A(62)+A(350)-A(514)
  Vdot(116) = A(555)-A(556)-A(557)
  Vdot(117) = -A(71)+A(339)-A(522)
  Vdot(118) = -A(45)+A(262)-A(483)
  Vdot(119) = A(493)-A(494)-A(495)
  Vdot(120) = -A(46)+A(263)-A(484)
  Vdot(121) = A(498)-A(499)-A(500)
  Vdot(122) = A(559)-A(560)-A(561)
  Vdot(123) = -A(87)+A(366)-A(538)
  Vdot(124) = -A(89)+A(368)-A(540)
  Vdot(125) = -A(92)+A(371)-A(543)
  Vdot(126) = A(544)-A(545)-A(550)
  Vdot(127) = -A(129)-A(130)+0.18*A(627)+A(629)
  Vdot(128) = -A(579)+A(610)-A(614)
  Vdot(129) = -A(74)+A(342)+A(344)-A(525)
  Vdot(130) = -A(80)+A(336)+A(337)-A(531)
  Vdot(131) = -A(16)-A(435)+A(476)+A(508)
  Vdot(132) = A(89)+A(92)+A(162)+A(164)+A(241)+A(243)+A(321)+A(324)+A(430)+A(433)-A(468)+A(488)
  Vdot(133) = -A(578)+A(602)+A(603)-A(607)
  Vdot(134) = -A(41)-A(42)+A(251)-A(474)
  Vdot(135) = -A(576)+A(585)-A(586)-A(587)
  Vdot(136) = -A(53)-A(54)+A(329)+A(332)-A(506)
  Vdot(137) = -A(73)+A(341)+A(343)-A(524)
  Vdot(138) = -A(47)+A(264)+A(265)-A(485)
  Vdot(139) = -A(15)-A(434)+A(475)+A(507)
  Vdot(140) = -A(19)-A(438)+A(479)+A(523)
  Vdot(141) = -A(577)+A(598)+A(599)-A(606)
  Vdot(142) = A(623)+2*A(624)-A(625)-A(626)
  Vdot(143) = -A(139)-A(142)-A(146)-A(147)
  Vdot(144) = A(151)-A(230)-A(311)-A(358)-A(420)
  Vdot(145) = A(143)-A(229)-A(310)-A(357)-A(419)
  Vdot(146) = A(155)-A(234)-A(314)-A(361)-A(423)
  Vdot(147) = A(141)-A(227)-A(308)-A(355)-A(417)
  Vdot(148) = -A(26)+A(204)+A(285)+A(438)-A(449)+A(513)+A(597)+A(601)+A(605)
  Vdot(149) = -A(17)+A(42)+A(54)+A(195)+A(276)+A(385)-A(436)+A(474)+A(506)
  Vdot(150) = -A(24)+A(45)+A(46)+A(78)+A(79)+A(166)+A(200)+A(201)+A(203)+A(281)+A(282)+A(284)+A(388)+A(389)+A(391)&
                &+A(392)-A(447)+A(483)+A(484)+A(529)+A(530)
  Vdot(151) = -A(30)+A(68)+A(230)+A(311)+A(420)-A(457)+A(519)
  Vdot(152) = A(579)+A(580)-A(581)-A(582)+A(609)+A(611)+A(613)+A(614)+A(615)-A(616)-A(617)
  Vdot(153) = -A(138)-A(141)-A(144)-A(145)
  Vdot(154) = -A(32)+A(45)+A(78)-A(93)+A(200)+A(202)+A(204)+A(281)+A(283)+A(285)+A(388)+A(390)-A(459)+A(483)+A(529)&
                &-A(563)
  Vdot(155) = A(588)+A(590)+A(592)+A(593)+A(594)+A(595)-A(596)-A(598)-A(600)-A(602)-A(604)
  Vdot(156) = A(135)-A(193)-A(250)-A(274)-A(328)-A(381)-A(382)-A(383)
  Vdot(157) = A(84)+A(236)-A(237)+A(316)-A(317)-A(364)+A(425)-A(426)+A(536)
  Vdot(158) = -A(36)+A(87)+A(239)+A(319)+A(428)-A(464)-A(466)
  Vdot(159) = -A(159)-A(160)-A(161)-A(162)-A(163)-A(164)
  Vdot(160) = A(583)-A(584)
  Vdot(161) = A(30)-A(31)+A(69)+A(231)+A(312)+A(421)-A(458)+A(520)
  Vdot(162) = -A(140)-A(143)-A(148)-A(149)
  Vdot(163) = A(167)-A(200)-A(262)-A(281)-A(334)-A(388)
  Vdot(164) = -A(154)-A(155)-A(156)-A(157)-A(158)
  Vdot(165) = -A(35)+A(81)+A(82)+A(157)+A(232)+A(234)+A(313)+A(314)+A(422)+A(423)-A(462)-A(465)+A(486)+A(532)+A(534)
  Vdot(166) = A(142)-A(228)-A(309)-A(356)-A(418)
  Vdot(167) = -A(34)-A(95)+A(203)+A(284)+A(391)+A(392)-A(461)-A(565)+A(577)+A(578)+A(596)+A(600)+A(604)+A(606)+A(607)
  Vdot(168) = A(160)-A(243)-A(324)-A(371)-A(433)+A(543)
  Vdot(169) = -A(20)+A(211)+A(292)+A(401)-A(439)+A(480)+A(494)+A(495)+A(524)
  Vdot(170) = A(44)-A(210)-A(258)-A(291)-A(342)-A(400)
  Vdot(171) = 0.82*A(627)+A(628)-A(629)-A(630)-A(631)+0.95*A(635)+A(636)
  Vdot(172) = A(625)+A(626)-A(627)-A(628)
  Vdot(173) = A(2)+A(4)+A(6)-A(96)-A(98)-A(99)-A(100)-A(101)+A(102)-A(129)
  Vdot(174) = A(608)-A(609)-A(610)-A(611)-A(612)-A(613)
  Vdot(175) = -A(150)-A(151)-A(152)-A(153)
  Vdot(176) = -A(198)+A(245)-A(279)-A(332)-A(386)
  Vdot(177) = -A(199)+A(219)-A(280)+A(300)-A(333)-A(387)+A(409)
  Vdot(178) = A(150)-A(220)-A(221)-A(261)-A(301)-A(302)-A(352)-A(410)-A(411)
  Vdot(179) = -A(18)+A(223)+A(226)+A(304)+A(307)+A(413)+A(416)-A(437)+A(478)+A(522)+A(557)
  Vdot(180) = -A(21)+A(212)+A(293)+A(402)-A(440)+A(481)+A(499)+A(500)+A(525)+A(561)
  Vdot(181) = A(169)-A(201)-A(202)-A(263)-A(282)-A(283)-A(335)-A(389)-A(390)
  Vdot(182) = -A(33)+A(46)+A(47)+A(79)+A(80)-A(94)+A(201)+A(282)+A(389)-A(460)+A(484)+A(485)+A(530)+A(531)-A(564)+A(597)&
                &+A(601)+A(605)
  Vdot(183) = -A(240)-A(320)-A(367)-A(429)+A(464)+A(466)+A(539)+A(576)+A(586)+A(587)
  Vdot(184) = A(159)-A(241)-A(242)-A(268)-A(321)-A(368)-A(430)+A(540)
  Vdot(185) = -A(8)+A(120)-A(128)+A(176)+A(177)+A(178)+A(442)+A(446)+A(454)+A(465)+A(466)+A(492)+A(495)+A(497)+A(500)&
                &+A(587)+A(622)+A(636)
  Vdot(186) = A(154)-A(232)-A(233)-A(266)-A(313)-A(360)-A(422)
  Vdot(187) = -A(13)+A(40)+A(52)+A(86)+A(91)+A(156)-A(179)+A(186)+A(193)+A(233)+A(238)+A(242)+A(245)+A(246)+A(274)&
                &+A(318)+A(323)+A(381)+A(382)+A(427)+A(432)+A(473)+A(487)+A(490)+A(505)+A(542)+A(561)+A(584)
  Vdot(188) = -A(231)-A(312)-A(359)-A(421)+A(457)
  Vdot(189) = A(9)+A(10)+A(11)+A(12)+A(22)+2*A(24)+A(25)+A(29)+2*A(30)+A(36)+A(69)+A(88)-A(113)+A(144)+A(146)+A(148)&
                &+A(152)+A(164)+A(165)+A(173)+A(176)+A(223)+A(231)+A(240)+A(304)+A(312)+A(320)+A(421)+A(429)+A(443)+2*A(447)&
                &+A(448)+A(455)+A(520)+A(550)+A(551)+A(552)+A(557)+A(576)+A(586)+A(587)
  Vdot(190) = A(86)+A(238)-A(239)+A(318)-A(319)-A(366)+A(427)-A(428)+A(538)
  Vdot(191) = -A(215)-A(296)-A(347)-A(405)+A(445)+A(446)+A(511)-A(548)+A(549)+A(575)
  Vdot(192) = A(171)-A(203)-A(264)-A(284)-A(336)-A(391)
  Vdot(193) = -A(25)+A(47)+A(66)+A(80)+A(202)+A(225)+A(283)+A(306)+A(390)+A(415)+A(437)+A(444)-A(448)+A(456)+A(458)&
                &+A(485)+A(512)+A(531)+A(577)+A(578)+A(579)+A(580)+A(596)+A(600)+A(604)+A(606)+A(607)+A(609)+A(611)+A(613)&
                &+A(614)+A(615)
  Vdot(194) = -A(29)+A(64)+A(220)+A(301)+A(410)-A(453)-A(454)-A(455)-A(456)+A(482)+A(516)
  Vdot(195) = A(56)+A(156)+A(197)-A(212)-A(260)+A(278)-A(293)-A(344)+A(394)-A(402)
  Vdot(196) = -A(22)+A(29)+A(66)+A(67)+A(206)+A(222)+A(224)+A(287)+A(303)+A(305)+A(396)+A(412)+A(413)+A(414)-A(445)&
                &-A(446)+A(455)+A(456)+A(477)+A(521)
  Vdot(197) = -A(224)-A(225)-A(226)-A(305)-A(306)-A(307)-A(354)-A(414)-A(415)-A(416)+A(441)+A(442)+A(518)-A(555)+A(556)
  Vdot(198) = -A(222)-A(223)-A(303)-A(304)-A(353)-A(412)-A(413)+A(453)+A(454)+A(517)-A(553)+A(554)
  Vdot(199) = A(138)+A(183)-A(205)-A(206)-A(254)-A(286)-A(287)-A(338)-A(395)-A(396)+A(562)+A(570)+A(571)+A(572)+A(573)
  Vdot(200) = A(589)+A(591)-A(597)-A(599)-A(601)-A(603)-A(605)
  Vdot(201) = -A(23)+A(65)+A(152)+A(153)+A(221)+A(302)+A(411)-A(441)-A(442)-A(443)-A(444)+A(558)
  Vdot(202) = A(172)-A(204)-A(265)-A(285)-A(337)-A(392)
  Vdot(203) = A(175)+A(178)-A(214)+A(218)-A(295)+A(299)-A(346)-A(404)+A(408)+A(510)-A(546)+A(547)
  Vdot(204) = A(85)+A(237)-A(238)+A(317)-A(318)-A(365)+A(426)-A(427)+A(537)-A(583)
  Vdot(205) = -A(197)+A(233)+A(242)-A(253)-A(278)-A(331)-A(394)
  Vdot(206) = -A(9)-A(10)+A(18)+A(22)+A(23)+A(31)+A(37)+A(49)+A(59)+A(65)+2*A(70)+A(71)+2*A(75)+A(76)+A(87)+A(89)+A(92)&
                &+A(144)+A(145)+A(146)+A(147)+A(162)+A(163)-A(173)-A(176)+A(181)+A(190)+2*A(205)+A(207)+A(215)+A(216)+A(221)&
                &+A(225)+A(226)+2*A(227)+A(228)+A(239)+A(241)+A(243)+A(245)+A(271)+2*A(286)+A(288)+A(296)+A(297)+A(302)&
                &+A(306)+A(307)+2*A(308)+A(309)+A(319)+A(321)+A(324)+A(372)+A(373)+2*A(395)+A(397)+A(405)+A(406)+A(411)&
                &+A(415)+A(416)+2*A(417)+A(418)+A(428)+A(430)+A(433)+A(443)+A(444)+A(470)+A(488)+A(502)+2*A(526)+A(527)&
                &+A(550)+A(552)+A(623)+2*A(624)
  Vdot(207) = A(43)+A(199)-A(209)-A(257)+A(280)-A(290)-A(341)+A(387)-A(399)+A(451)
  Vdot(208) = A(180)-A(217)-A(298)-A(349)-A(407)+A(436)
  Vdot(209) = A(12)+A(14)+A(17)+A(32)+A(41)+A(53)+A(58)+A(60)+A(95)+A(134)-A(191)+A(194)+A(214)-A(248)-A(272)+A(275)&
                &+A(295)-A(326)-A(375)-A(376)-A(377)+A(384)+A(404)+A(459)+A(565)+A(569)
  Vdot(210) = -A(246)-A(270)-A(323)-A(370)-A(432)+A(469)
  Vdot(211) = A(11)+A(13)+A(57)+A(94)+A(97)+A(133)+A(146)+A(148)+A(189)-A(190)+A(213)-A(247)-A(271)+A(294)-A(325)-A(372)&
                &-A(373)-A(374)+A(403)+A(501)+A(564)+A(566)+A(567)+A(568)-A(574)+0.18*A(627)+A(629)+A(633)+0.95*A(634)
  Vdot(212) = A(83)+A(235)-A(236)-A(267)+A(315)-A(316)-A(363)+A(424)-A(425)+A(467)+A(535)
  Vdot(213) = A(140)+A(198)-A(208)-A(256)+A(279)-A(289)-A(340)+A(386)-A(398)+A(450)+A(452)
  Vdot(214) = A(161)-A(244)-A(245)-A(269)-A(322)-A(369)-A(431)+A(468)
  Vdot(215) = A(15)+A(16)+A(33)+A(61)+A(136)-A(192)-A(249)-A(273)-A(327)-A(378)-A(379)-A(380)+A(460)+A(581)+A(617)
  Vdot(216) = -A(218)-A(299)-A(350)-A(408)+A(434)
  Vdot(217) = A(139)+A(185)+A(187)-A(207)-A(255)-A(288)-A(339)-A(397)+A(449)+A(575)
  Vdot(218) = -A(219)-A(300)-A(351)-A(409)+A(435)
  Vdot(219) = A(13)+A(14)+A(15)+A(17)+A(18)+A(19)+A(23)+A(25)+2*A(26)+A(27)+A(29)+A(31)+2*A(36)+A(67)+A(174)+A(177)&
                &-A(213)+A(216)+A(217)+A(222)+A(224)-A(294)+A(297)+A(298)+A(303)+A(305)-A(345)-A(403)+A(406)+A(407)+A(412)&
                &+A(414)+A(443)+A(448)+A(455)+A(509)-A(544)+A(545)
  Vdot(220) = A(16)+A(20)+A(21)+A(27)+2*A(28)+A(88)+A(179)-A(216)+A(240)-A(297)+A(320)-A(348)-A(406)+A(429)
  Vdot(221) = A(35)-A(235)-A(315)-A(362)-A(424)+A(462)+A(465)+A(533)-A(559)+A(560)
  Vdot(222) = -A(1)-A(2)+A(96)-A(98)-A(103)-A(104)-A(111)-A(115)-A(144)-A(145)-A(146)-A(147)-A(148)-A(149)-A(152)-A(153)&
                &-A(156)-A(157)-A(158)-A(161)-A(162)-A(163)-A(164)-A(443)-A(444)-A(455)-A(456)-A(625)-A(628)-A(631)
  Vdot(223) = A(63)-A(196)-A(252)-A(277)-A(330)-A(393)+A(440)+A(582)+A(616)
  Vdot(224) = A(55)+A(196)-A(211)-A(259)+A(277)-A(292)-A(343)+A(393)-A(401)+A(463)
  Vdot(225) = A(34)+A(62)+A(91)+A(137)-A(194)-A(195)+A(246)-A(251)-A(275)-A(276)+A(323)-A(329)-A(384)-A(385)+A(432)&
                &+A(439)+A(461)
  Vdot(226) = -A(12)+A(21)+A(39)+A(51)+A(73)+2*A(74)-A(175)-A(178)+A(184)+A(192)+A(209)+2*A(210)+A(273)+A(290)+2*A(291)&
                &+A(378)+A(379)+A(399)+2*A(400)+A(472)+A(504)
  Vdot(227) = A(4)+A(5)+A(7)-A(99)+A(100)-A(103)-2*A(105)-A(106)+A(107)-A(119)-A(122)-A(190)-A(191)-A(192)-A(193)-A(194)&
                &-A(195)-A(196)-A(197)-A(198)-A(199)-A(200)-A(201)-A(202)-A(203)-A(204)-A(205)-A(206)-A(207)-A(208)-A(209)&
                &-A(210)-A(211)-A(212)-A(213)-A(214)-A(215)-A(216)-A(217)-A(218)-A(219)-A(220)-A(221)-A(222)-A(223)-A(224)&
                &-A(225)-A(226)-A(227)-A(228)-A(229)-A(230)-A(231)-A(232)-A(233)-A(234)-A(235)-A(236)-A(237)-A(238)-A(239)&
                &-A(240)-A(241)-A(242)-A(243)-A(244)-A(245)-A(246)-A(247)-A(248)-A(249)-A(250)-A(251)-A(252)-A(253)-A(254)&
                &-A(255)-A(256)-A(257)-A(258)-A(259)-A(260)-A(261)-A(262)-A(263)-A(264)-A(265)-A(266)-A(267)-A(268)-A(269)&
                &-A(270)-A(583)-A(596)-A(597)-A(598)-A(599)-A(609)-A(610)-A(623)+A(626)+A(627)+A(630)
  Vdot(228) = 2*A(3)+A(7)+A(8)+A(49)+A(50)+A(51)+A(52)+A(53)+A(54)+A(55)+A(56)+A(57)+A(58)+A(59)+A(60)+A(61)+A(62)+A(63)&
                &+A(64)+A(65)+A(66)+A(67)+A(68)+A(69)+A(70)+A(71)+A(72)+A(73)+A(74)+A(75)+A(76)+A(77)+A(78)+A(79)+A(80)&
                &+A(81)+A(82)+A(83)+A(84)+A(85)+A(86)+A(87)+A(88)+A(89)+A(90)+A(91)+A(92)+2*A(110)-A(111)-A(112)-A(113)&
                &-A(114)+A(115)-A(116)-A(119)-A(120)-A(121)+A(122)-A(125)+A(126)-A(127)-A(128)-A(130)-A(133)-A(134)-A(135)&
                &-A(136)-A(137)-A(138)-A(139)-A(140)+A(144)+A(146)+A(148)-A(150)+A(152)-A(154)+A(156)-A(159)+A(161)-A(165)&
                &-A(167)-A(168)-A(169)-A(170)-A(171)-A(172)-A(173)-A(174)-A(175)-A(179)-A(180)-A(181)-A(182)-A(183)-A(184)&
                &-A(185)-A(186)-A(187)-A(188)-A(189)-A(434)-A(435)-A(436)-A(437)-A(438)-A(439)-A(440)-A(441)+A(443)-A(445)&
                &-A(447)-A(448)-A(449)-A(450)-A(451)-A(452)-A(453)+A(455)-A(457)-A(458)-A(459)-A(460)-A(461)-A(462)-A(463)&
                &-A(464)-A(467)-A(468)-A(469)-A(470)-A(471)-A(472)-A(473)-A(474)-A(475)-A(476)-A(477)-A(478)-A(479)-A(480)&
                &-A(481)-A(482)-A(483)-A(484)-A(485)-A(486)-A(487)-A(488)-A(489)-A(490)-A(491)-A(494)-A(496)-A(499)-A(501)&
                &-A(509)-A(510)-A(511)-A(517)-A(518)-A(533)-A(535)-A(536)-A(537)-A(538)-A(539)-A(540)-A(543)-A(550)-A(551)&
                &-A(552)-A(557)-A(558)-A(561)-A(562)-A(563)-A(564)-A(565)-A(566)-A(567)-A(568)-A(569)-A(570)-A(571)-A(572)&
                &-A(573)+A(578)+A(580)-A(584)-A(585)-A(586)-A(588)-A(589)-A(590)-A(591)-A(592)-A(593)-A(594)-A(595)-A(606)&
                &-A(608)-A(614)-A(616)-A(617)-A(620)-A(621)-A(634)-A(635)
  Vdot(229) = -A(4)+A(6)+A(8)+A(31)+A(37)+A(38)+A(39)+A(40)+A(41)+A(42)+A(43)+A(44)+A(45)+A(46)+A(47)+A(48)+A(75)+A(76)&
                &+A(77)+A(82)+A(92)+A(97)+A(99)-A(100)-A(101)+A(103)-A(104)+2*A(105)+2*A(106)-A(108)+A(109)-A(120)+A(121)&
                &+A(122)-A(123)+A(124)+A(125)+A(126)+A(127)+A(190)+A(191)+A(192)+A(193)+A(194)+A(195)+A(196)+A(197)+A(198)&
                &+A(199)+A(200)+A(201)+A(202)+A(203)+A(204)+A(205)+A(206)+A(207)+A(208)+A(209)+A(210)+A(211)+A(212)+A(213)&
                &+A(214)+A(215)+A(216)+A(217)+A(218)+A(219)+A(220)+A(221)+A(222)+A(223)+A(224)+A(225)+A(226)+2*A(227)+2&
                &*A(228)+2*A(229)+A(230)+A(231)+A(232)+A(233)+2*A(234)+A(235)+A(236)+A(237)+A(238)+A(239)+A(240)+A(241)&
                &+A(242)+2*A(243)+A(244)+A(245)+A(246)+A(271)+A(272)+A(273)+A(274)+A(275)+A(276)+A(277)+A(278)+A(279)+A(280)&
                &+A(281)+A(282)+A(283)+A(284)+A(285)+A(286)+A(287)+A(288)+A(289)+A(290)+A(291)+A(292)+A(293)+A(294)+A(295)&
                &+A(296)+A(297)+A(298)+A(299)+A(300)+A(301)+A(302)+A(303)+A(304)+A(305)+A(306)+A(307)+2*A(308)+2*A(309)+2&
                &*A(310)+A(311)+A(312)+A(313)+2*A(314)+A(315)+A(316)+A(317)+A(318)+A(319)+A(320)+A(321)+A(322)+A(323)+2&
                &*A(324)+A(417)+A(418)+A(419)+A(423)+A(433)+A(458)+A(470)+A(471)+A(472)+A(473)+A(474)+A(475)+A(476)+A(477)&
                &+A(478)+A(479)+A(480)+A(481)+A(482)+A(483)+A(484)+A(485)+A(486)+A(487)+A(488)+A(489)+A(490)-A(493)+A(494)&
                &+A(495)-A(498)+A(499)+A(500)+A(526)+A(527)+A(528)+A(534)-A(544)+A(545)-A(546)+A(547)-A(548)+A(549)+A(550)&
                &+A(551)+A(552)-A(553)+A(554)-A(555)+A(556)+A(557)+A(558)-A(559)+A(560)+A(561)-A(574)+A(577)+A(579)+A(584)&
                &+A(596)+A(597)+A(600)+A(601)+A(606)+A(609)+A(611)+A(614)+A(623)-A(626)-A(627)-A(630)
  Vdot(230) = -A(5)-A(6)+A(101)+A(104)-A(106)-A(107)-A(108)+A(109)-A(121)-A(126)+A(128)-A(141)-A(142)-A(143)-A(151)&
                &-A(155)-A(160)-A(176)-A(177)-A(178)-A(271)-A(272)-A(273)-A(274)-A(275)-A(276)-A(277)-A(278)-A(279)-A(280)&
                &-A(281)-A(282)-A(283)-A(284)-A(285)-A(286)-A(287)-A(288)-A(289)-A(290)-A(291)-A(292)-A(293)-A(294)-A(295)&
                &-A(296)-A(297)-A(298)-A(299)-A(300)-A(301)-A(302)-A(303)-A(304)-A(305)-A(306)-A(307)-A(308)-A(309)-A(310)&
                &-A(311)-A(312)-A(313)-A(314)-A(315)-A(316)-A(317)-A(318)-A(319)-A(320)-A(321)-A(322)-A(323)-A(324)-A(442)&
                &-A(446)-A(454)-A(465)-A(466)-A(492)-A(495)-A(497)-A(500)-A(587)-A(600)-A(601)-A(611)-A(622)-A(636)
  Vdot(231) = -A(11)+A(19)+A(20)+A(38)+A(41)+A(50)+A(53)+A(71)+2*A(72)+A(73)+A(76)+2*A(77)+A(148)+A(149)-A(174)-A(177)&
                &+A(182)+A(191)+A(194)+A(207)+2*A(208)+A(209)+A(217)+A(218)+A(228)+2*A(229)+A(272)+A(275)+A(288)+2*A(289)&
                &+A(290)+A(298)+A(299)+A(309)+2*A(310)+A(375)+A(376)+A(384)+A(397)+2*A(398)+A(399)+A(407)+A(408)+A(418)+2&
                &*A(419)+A(471)+A(503)+A(527)+2*A(528)+A(551)
  Vdot(232) = 2*A(9)+A(11)+A(12)+A(18)+A(19)+A(20)+A(21)+2*A(22)+A(23)+2*A(24)+A(25)+A(29)+2*A(30)+A(32)+A(33)+A(34)&
                &+A(35)+A(37)+A(38)+A(39)+A(40)+A(42)+A(45)+A(46)+A(47)+A(48)+A(49)+A(50)+A(51)+A(52)+A(54)+A(59)+A(64)&
                &+A(65)+A(66)+A(68)+A(69)+A(70)+A(71)+A(72)+A(73)+A(74)+A(78)+A(79)+A(80)+A(81)+A(87)+A(89)+A(90)+2*A(93)&
                &+A(94)+A(95)+A(111)+A(112)+A(113)+A(114)-A(115)-A(116)-2*A(117)-2*A(118)+A(121)-A(122)-A(123)+A(124)-A(126)&
                &+A(131)+A(144)+A(152)+A(165)+A(168)+A(170)+A(173)+A(176)+A(181)+A(182)+A(184)+A(186)+A(188)+A(190)+A(191)&
                &+A(192)+A(193)+A(195)+A(200)+A(201)+A(202)+A(203)+A(204)+A(205)+A(206)+A(207)+A(208)+A(209)+A(210)+A(211)&
                &+A(212)+A(215)+A(220)+A(221)+A(223)+A(225)+A(226)+A(230)+A(231)+A(232)+A(239)+A(241)+A(244)+A(271)+A(272)&
                &+A(273)+A(274)+A(276)+A(281)+A(282)+A(283)+A(284)+A(285)+A(286)+A(287)+A(288)+A(289)+A(290)+A(291)+A(292)&
                &+A(293)+A(296)+A(301)+A(302)+A(304)+A(306)+A(307)+A(311)+A(312)+A(313)+A(319)+A(321)+A(322)-A(325)-A(326)&
                &-A(327)-A(328)-A(329)-A(330)-A(331)-A(332)-A(333)-A(334)-A(335)-A(336)-A(337)-A(338)-A(339)-A(340)-A(341)&
                &-A(342)-A(343)-A(344)-A(345)-A(346)-A(347)-A(348)-A(349)-A(350)-A(351)-A(352)-A(353)-A(354)-A(355)-A(356)&
                &-A(357)-A(358)-A(359)-A(360)-A(361)-A(362)-A(363)-A(364)-A(365)-A(366)-A(367)-A(368)-A(369)-A(370)-A(371)&
                &+A(372)+A(375)+A(378)+A(381)+A(385)+A(388)+A(389)+A(390)+A(391)+A(392)+A(395)+A(396)+A(397)+A(398)+A(399)&
                &+A(400)+A(401)+A(402)+A(405)+A(410)+A(411)+A(413)+A(415)+A(416)+A(420)+A(421)+A(422)+A(428)+A(430)+A(431)&
                &+A(437)+A(438)+A(447)+A(563)+A(576)+A(577)+A(578)+A(579)+A(580)+A(581)+A(582)+A(585)+A(596)+A(597)+A(600)&
                &+A(601)-A(602)-A(603)+A(604)+A(605)+A(609)+A(611)-A(612)+A(613)+A(621)-A(632)+0.05*A(635)
      
END SUBROUTINE crimech_Fun
















SUBROUTINE crimech_IRRFun ( V, F, RCT, IRR )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: IRR(NREACT)



  IRR(1) = RCT(1)*V(222)
  IRR(2) = RCT(2)*V(222)
  IRR(3) = RCT(3)*V(106)
  IRR(4) = RCT(4)*V(229)
  IRR(5) = RCT(5)*V(230)
  IRR(6) = RCT(6)*V(230)
  IRR(7) = RCT(7)*V(40)
  IRR(8) = RCT(8)*V(185)
  IRR(9) = RCT(9)*V(206)
  IRR(10) = RCT(10)*V(206)
  IRR(11) = RCT(11)*V(231)
  IRR(12) = RCT(12)*V(226)
  IRR(13) = RCT(13)*V(187)
  IRR(14) = RCT(14)*V(31)
  IRR(15) = RCT(15)*V(139)
  IRR(16) = RCT(16)*V(131)
  IRR(17) = RCT(17)*V(149)
  IRR(18) = RCT(18)*V(179)
  IRR(19) = RCT(19)*V(140)
  IRR(20) = RCT(20)*V(169)
  IRR(21) = RCT(21)*V(180)
  IRR(22) = RCT(22)*V(196)
  IRR(23) = RCT(23)*V(201)
  IRR(24) = RCT(24)*V(150)
  IRR(25) = RCT(25)*V(193)
  IRR(26) = RCT(26)*V(148)
  IRR(27) = RCT(27)*V(44)
  IRR(28) = RCT(28)*V(37)
  IRR(29) = RCT(29)*V(194)
  IRR(30) = RCT(30)*V(151)
  IRR(31) = RCT(31)*V(161)
  IRR(32) = RCT(32)*V(154)
  IRR(33) = RCT(33)*V(182)
  IRR(34) = RCT(34)*V(167)
  IRR(35) = RCT(35)*V(165)
  IRR(36) = RCT(36)*V(158)
  IRR(37) = RCT(37)*V(90)
  IRR(38) = RCT(38)*V(92)
  IRR(39) = RCT(39)*V(98)
  IRR(40) = RCT(40)*V(85)
  IRR(41) = RCT(41)*V(134)
  IRR(42) = RCT(42)*V(134)
  IRR(43) = RCT(43)*V(79)
  IRR(44) = RCT(44)*V(103)
  IRR(45) = RCT(45)*V(118)
  IRR(46) = RCT(46)*V(120)
  IRR(47) = RCT(47)*V(138)
  IRR(48) = RCT(48)*V(97)
  IRR(49) = RCT(49)*V(101)
  IRR(50) = RCT(50)*V(73)
  IRR(51) = RCT(51)*V(56)
  IRR(52) = RCT(52)*V(75)
  IRR(53) = RCT(53)*V(136)
  IRR(54) = RCT(54)*V(136)
  IRR(55) = RCT(55)*V(80)
  IRR(56) = RCT(56)*V(57)
  IRR(57) = RCT(57)*V(53)
  IRR(58) = RCT(58)*V(62)
  IRR(59) = RCT(59)*V(84)
  IRR(60) = RCT(60)*V(58)
  IRR(61) = RCT(61)*V(59)
  IRR(62) = RCT(62)*V(115)
  IRR(63) = RCT(63)*V(81)
  IRR(64) = RCT(64)*V(113)
  IRR(65) = RCT(65)*V(113)
  IRR(66) = RCT(66)*V(114)
  IRR(67) = RCT(67)*V(86)
  IRR(68) = RCT(68)*V(64)
  IRR(69) = RCT(69)*V(88)
  IRR(70) = RCT(70)*V(105)
  IRR(71) = RCT(71)*V(117)
  IRR(72) = RCT(72)*V(82)
  IRR(73) = RCT(73)*V(137)
  IRR(74) = RCT(74)*V(129)
  IRR(75) = RCT(75)*V(60)
  IRR(76) = RCT(76)*V(83)
  IRR(77) = RCT(77)*V(61)
  IRR(78) = RCT(78)*V(89)
  IRR(79) = RCT(79)*V(91)
  IRR(80) = RCT(80)*V(130)
  IRR(81) = RCT(81)*V(67)
  IRR(82) = RCT(82)*V(69)
  IRR(83) = RCT(83)*V(68)
  IRR(84) = RCT(84)*V(70)
  IRR(85) = RCT(85)*V(71)
  IRR(86) = RCT(86)*V(93)
  IRR(87) = RCT(87)*V(123)
  IRR(88) = RCT(88)*V(94)
  IRR(89) = RCT(89)*V(124)
  IRR(90) = RCT(90)*V(72)
  IRR(91) = RCT(91)*V(100)
  IRR(92) = RCT(92)*V(125)
  IRR(93) = RCT(93)*V(154)
  IRR(94) = RCT(94)*V(182)
  IRR(95) = RCT(95)*V(167)
  IRR(96) = RCT(96)*V(173)*F(2)
  IRR(97) = RCT(97)*V(34)
  IRR(98) = RCT(98)*V(173)*V(222)
  IRR(99) = RCT(99)*V(173)*V(227)
  IRR(100) = RCT(100)*V(173)*V(229)
  IRR(101) = RCT(101)*V(173)*V(229)
  IRR(102) = RCT(102)*V(7)*F(2)
  IRR(103) = RCT(103)*V(222)*V(227)
  IRR(104) = RCT(104)*V(222)*V(229)
  IRR(105) = RCT(105)*V(227)*V(227)*F(2)
  IRR(106) = RCT(106)*V(227)*V(230)
  IRR(107) = RCT(107)*V(229)*V(230)
  IRR(108) = RCT(108)*V(229)*V(230)
  IRR(109) = RCT(109)*V(32)
  IRR(110) = RCT(110)*V(7)*F(1)
  IRR(111) = RCT(111)*V(222)*V(228)
  IRR(112) = RCT(112)*V(228)*F(2)
  IRR(113) = RCT(113)*V(189)*V(228)
  IRR(114) = RCT(114)*V(106)*V(228)
  IRR(115) = RCT(115)*V(222)*V(232)
  IRR(116) = RCT(116)*V(228)*V(232)
  IRR(117) = RCT(117)*V(232)*V(232)
  IRR(118) = RCT(118)*V(232)*V(232)*F(1)
  IRR(119) = RCT(119)*V(227)*V(228)
  IRR(120) = RCT(120)*V(228)*V(229)
  IRR(121) = RCT(121)*V(228)*V(230)
  IRR(122) = RCT(122)*V(227)*V(232)
  IRR(123) = RCT(123)*V(229)*V(232)
  IRR(124) = RCT(124)*V(74)
  IRR(125) = RCT(125)*V(74)*V(228)
  IRR(126) = RCT(126)*V(230)*V(232)
  IRR(127) = RCT(127)*V(40)*V(228)
  IRR(128) = RCT(128)*V(185)*V(228)
  IRR(129) = RCT(129)*V(127)*V(173)
  IRR(130) = RCT(130)*V(127)*V(228)
  IRR(131) = RCT(131)*V(9)*F(2)
  IRR(132) = RCT(132)*V(4)*F(1)*F(1)
  IRR(133) = RCT(133)*V(19)*V(228)
  IRR(134) = RCT(134)*V(8)*V(228)
  IRR(135) = RCT(135)*V(25)*V(228)
  IRR(136) = RCT(136)*V(25)*V(228)
  IRR(137) = RCT(137)*V(10)*V(228)
  IRR(138) = RCT(138)*V(153)*V(228)
  IRR(139) = RCT(139)*V(143)*V(228)
  IRR(140) = RCT(140)*V(162)*V(228)
  IRR(141) = RCT(141)*V(153)*V(230)
  IRR(142) = RCT(142)*V(143)*V(230)
  IRR(143) = RCT(143)*V(162)*V(230)
  IRR(144) = RCT(144)*V(153)*V(222)
  IRR(145) = RCT(145)*V(153)*V(222)
  IRR(146) = RCT(146)*V(143)*V(222)
  IRR(147) = RCT(147)*V(143)*V(222)
  IRR(148) = RCT(148)*V(162)*V(222)
  IRR(149) = RCT(149)*V(162)*V(222)
  IRR(150) = RCT(150)*V(175)*V(228)
  IRR(151) = RCT(151)*V(175)*V(230)
  IRR(152) = RCT(152)*V(175)*V(222)
  IRR(153) = RCT(153)*V(175)*V(222)
  IRR(154) = RCT(154)*V(164)*V(228)
  IRR(155) = RCT(155)*V(164)*V(230)
  IRR(156) = RCT(156)*V(164)*V(222)
  IRR(157) = RCT(157)*V(164)*V(222)
  IRR(158) = RCT(158)*V(164)*V(222)
  IRR(159) = RCT(159)*V(159)*V(228)
  IRR(160) = RCT(160)*V(159)*V(230)
  IRR(161) = RCT(161)*V(159)*V(222)
  IRR(162) = RCT(162)*V(159)*V(222)
  IRR(163) = RCT(163)*V(159)*V(222)
  IRR(164) = RCT(164)*V(159)*V(222)
  IRR(165) = RCT(165)*V(35)*V(228)
  IRR(166) = RCT(166)*V(35)*V(228)
  IRR(167) = RCT(167)*V(29)*V(228)
  IRR(168) = RCT(168)*V(29)*V(228)
  IRR(169) = RCT(169)*V(30)*V(228)
  IRR(170) = RCT(170)*V(30)*V(228)
  IRR(171) = RCT(171)*V(26)*V(228)
  IRR(172) = RCT(172)*V(26)*V(228)
  IRR(173) = RCT(173)*V(206)*V(228)
  IRR(174) = RCT(174)*V(228)*V(231)
  IRR(175) = RCT(175)*V(226)*V(228)
  IRR(176) = RCT(176)*V(206)*V(230)
  IRR(177) = RCT(177)*V(230)*V(231)
  IRR(178) = RCT(178)*V(226)*V(230)
  IRR(179) = RCT(179)*V(187)*V(228)
  IRR(180) = RCT(180)*V(31)*V(228)
  IRR(181) = RCT(181)*V(36)*V(228)
  IRR(182) = RCT(182)*V(42)*V(228)
  IRR(183) = RCT(183)*V(42)*V(228)
  IRR(184) = RCT(184)*V(50)*V(228)
  IRR(185) = RCT(185)*V(50)*V(228)
  IRR(186) = RCT(186)*V(43)*V(228)
  IRR(187) = RCT(187)*V(43)*V(228)
  IRR(188) = RCT(188)*V(45)*V(228)
  IRR(189) = RCT(189)*V(51)*V(228)
  IRR(190) = RCT(190)*V(211)*V(227)
  IRR(191) = RCT(191)*V(209)*V(227)
  IRR(192) = RCT(192)*V(215)*V(227)
  IRR(193) = RCT(193)*V(156)*V(227)
  IRR(194) = RCT(194)*V(225)*V(227)
  IRR(195) = RCT(195)*V(225)*V(227)
  IRR(196) = RCT(196)*V(223)*V(227)
  IRR(197) = RCT(197)*V(205)*V(227)
  IRR(198) = RCT(198)*V(176)*V(227)
  IRR(199) = RCT(199)*V(177)*V(227)
  IRR(200) = RCT(200)*V(163)*V(227)
  IRR(201) = RCT(201)*V(181)*V(227)
  IRR(202) = RCT(202)*V(181)*V(227)
  IRR(203) = RCT(203)*V(192)*V(227)
  IRR(204) = RCT(204)*V(202)*V(227)
  IRR(205) = RCT(205)*V(199)*V(227)
  IRR(206) = RCT(206)*V(199)*V(227)
  IRR(207) = RCT(207)*V(217)*V(227)
  IRR(208) = RCT(208)*V(213)*V(227)
  IRR(209) = RCT(209)*V(207)*V(227)
  IRR(210) = RCT(210)*V(170)*V(227)
  IRR(211) = RCT(211)*V(224)*V(227)
  IRR(212) = RCT(212)*V(195)*V(227)
  IRR(213) = RCT(213)*V(219)*V(227)
  IRR(214) = RCT(214)*V(203)*V(227)
  IRR(215) = RCT(215)*V(191)*V(227)
  IRR(216) = RCT(216)*V(220)*V(227)
  IRR(217) = RCT(217)*V(208)*V(227)
  IRR(218) = RCT(218)*V(216)*V(227)
  IRR(219) = RCT(219)*V(218)*V(227)
  IRR(220) = RCT(220)*V(178)*V(227)
  IRR(221) = RCT(221)*V(178)*V(227)
  IRR(222) = RCT(222)*V(198)*V(227)
  IRR(223) = RCT(223)*V(198)*V(227)
  IRR(224) = RCT(224)*V(197)*V(227)
  IRR(225) = RCT(225)*V(197)*V(227)
  IRR(226) = RCT(226)*V(197)*V(227)
  IRR(227) = RCT(227)*V(147)*V(227)
  IRR(228) = RCT(228)*V(166)*V(227)
  IRR(229) = RCT(229)*V(145)*V(227)
  IRR(230) = RCT(230)*V(144)*V(227)
  IRR(231) = RCT(231)*V(188)*V(227)
  IRR(232) = RCT(232)*V(186)*V(227)
  IRR(233) = RCT(233)*V(186)*V(227)
  IRR(234) = RCT(234)*V(146)*V(227)
  IRR(235) = RCT(235)*V(221)*V(227)
  IRR(236) = RCT(236)*V(212)*V(227)
  IRR(237) = RCT(237)*V(157)*V(227)
  IRR(238) = RCT(238)*V(204)*V(227)
  IRR(239) = RCT(239)*V(190)*V(227)
  IRR(240) = RCT(240)*V(183)*V(227)
  IRR(241) = RCT(241)*V(184)*V(227)
  IRR(242) = RCT(242)*V(184)*V(227)
  IRR(243) = RCT(243)*V(168)*V(227)
  IRR(244) = RCT(244)*V(214)*V(227)
  IRR(245) = RCT(245)*V(214)*V(227)
  IRR(246) = RCT(246)*V(210)*V(227)
  IRR(247) = RCT(247)*V(211)*V(227)
  IRR(248) = RCT(248)*V(209)*V(227)
  IRR(249) = RCT(249)*V(215)*V(227)
  IRR(250) = RCT(250)*V(156)*V(227)
  IRR(251) = RCT(251)*V(225)*V(227)
  IRR(252) = RCT(252)*V(223)*V(227)
  IRR(253) = RCT(253)*V(205)*V(227)
  IRR(254) = RCT(254)*V(199)*V(227)
  IRR(255) = RCT(255)*V(217)*V(227)
  IRR(256) = RCT(256)*V(213)*V(227)
  IRR(257) = RCT(257)*V(207)*V(227)
  IRR(258) = RCT(258)*V(170)*V(227)
  IRR(259) = RCT(259)*V(224)*V(227)
  IRR(260) = RCT(260)*V(195)*V(227)
  IRR(261) = RCT(261)*V(178)*V(227)
  IRR(262) = RCT(262)*V(163)*V(227)
  IRR(263) = RCT(263)*V(181)*V(227)
  IRR(264) = RCT(264)*V(192)*V(227)
  IRR(265) = RCT(265)*V(202)*V(227)
  IRR(266) = RCT(266)*V(186)*V(227)
  IRR(267) = RCT(267)*V(212)*V(227)
  IRR(268) = RCT(268)*V(184)*V(227)
  IRR(269) = RCT(269)*V(214)*V(227)
  IRR(270) = RCT(270)*V(210)*V(227)
  IRR(271) = RCT(271)*V(211)*V(230)
  IRR(272) = RCT(272)*V(209)*V(230)
  IRR(273) = RCT(273)*V(215)*V(230)
  IRR(274) = RCT(274)*V(156)*V(230)
  IRR(275) = RCT(275)*V(225)*V(230)
  IRR(276) = RCT(276)*V(225)*V(230)
  IRR(277) = RCT(277)*V(223)*V(230)
  IRR(278) = RCT(278)*V(205)*V(230)
  IRR(279) = RCT(279)*V(176)*V(230)
  IRR(280) = RCT(280)*V(177)*V(230)
  IRR(281) = RCT(281)*V(163)*V(230)
  IRR(282) = RCT(282)*V(181)*V(230)
  IRR(283) = RCT(283)*V(181)*V(230)
  IRR(284) = RCT(284)*V(192)*V(230)
  IRR(285) = RCT(285)*V(202)*V(230)
  IRR(286) = RCT(286)*V(199)*V(230)
  IRR(287) = RCT(287)*V(199)*V(230)
  IRR(288) = RCT(288)*V(217)*V(230)
  IRR(289) = RCT(289)*V(213)*V(230)
  IRR(290) = RCT(290)*V(207)*V(230)
  IRR(291) = RCT(291)*V(170)*V(230)
  IRR(292) = RCT(292)*V(224)*V(230)
  IRR(293) = RCT(293)*V(195)*V(230)
  IRR(294) = RCT(294)*V(219)*V(230)
  IRR(295) = RCT(295)*V(203)*V(230)
  IRR(296) = RCT(296)*V(191)*V(230)
  IRR(297) = RCT(297)*V(220)*V(230)
  IRR(298) = RCT(298)*V(208)*V(230)
  IRR(299) = RCT(299)*V(216)*V(230)
  IRR(300) = RCT(300)*V(218)*V(230)
  IRR(301) = RCT(301)*V(178)*V(230)
  IRR(302) = RCT(302)*V(178)*V(230)
  IRR(303) = RCT(303)*V(198)*V(230)
  IRR(304) = RCT(304)*V(198)*V(230)
  IRR(305) = RCT(305)*V(197)*V(230)
  IRR(306) = RCT(306)*V(197)*V(230)
  IRR(307) = RCT(307)*V(197)*V(230)
  IRR(308) = RCT(308)*V(147)*V(230)
  IRR(309) = RCT(309)*V(166)*V(230)
  IRR(310) = RCT(310)*V(145)*V(230)
  IRR(311) = RCT(311)*V(144)*V(230)
  IRR(312) = RCT(312)*V(188)*V(230)
  IRR(313) = RCT(313)*V(186)*V(230)
  IRR(314) = RCT(314)*V(146)*V(230)
  IRR(315) = RCT(315)*V(221)*V(230)
  IRR(316) = RCT(316)*V(212)*V(230)
  IRR(317) = RCT(317)*V(157)*V(230)
  IRR(318) = RCT(318)*V(204)*V(230)
  IRR(319) = RCT(319)*V(190)*V(230)
  IRR(320) = RCT(320)*V(183)*V(230)
  IRR(321) = RCT(321)*V(184)*V(230)
  IRR(322) = RCT(322)*V(214)*V(230)
  IRR(323) = RCT(323)*V(210)*V(230)
  IRR(324) = RCT(324)*V(168)*V(230)
  IRR(325) = RCT(325)*V(211)*V(232)
  IRR(326) = RCT(326)*V(209)*V(232)
  IRR(327) = RCT(327)*V(215)*V(232)
  IRR(328) = RCT(328)*V(156)*V(232)
  IRR(329) = RCT(329)*V(225)*V(232)
  IRR(330) = RCT(330)*V(223)*V(232)
  IRR(331) = RCT(331)*V(205)*V(232)
  IRR(332) = RCT(332)*V(176)*V(232)
  IRR(333) = RCT(333)*V(177)*V(232)
  IRR(334) = RCT(334)*V(163)*V(232)
  IRR(335) = RCT(335)*V(181)*V(232)
  IRR(336) = RCT(336)*V(192)*V(232)
  IRR(337) = RCT(337)*V(202)*V(232)
  IRR(338) = RCT(338)*V(199)*V(232)
  IRR(339) = RCT(339)*V(217)*V(232)
  IRR(340) = RCT(340)*V(213)*V(232)
  IRR(341) = RCT(341)*V(207)*V(232)
  IRR(342) = RCT(342)*V(170)*V(232)
  IRR(343) = RCT(343)*V(224)*V(232)
  IRR(344) = RCT(344)*V(195)*V(232)
  IRR(345) = RCT(345)*V(219)*V(232)
  IRR(346) = RCT(346)*V(203)*V(232)
  IRR(347) = RCT(347)*V(191)*V(232)
  IRR(348) = RCT(348)*V(220)*V(232)
  IRR(349) = RCT(349)*V(208)*V(232)
  IRR(350) = RCT(350)*V(216)*V(232)
  IRR(351) = RCT(351)*V(218)*V(232)
  IRR(352) = RCT(352)*V(178)*V(232)
  IRR(353) = RCT(353)*V(198)*V(232)
  IRR(354) = RCT(354)*V(197)*V(232)
  IRR(355) = RCT(355)*V(147)*V(232)
  IRR(356) = RCT(356)*V(166)*V(232)
  IRR(357) = RCT(357)*V(145)*V(232)
  IRR(358) = RCT(358)*V(144)*V(232)
  IRR(359) = RCT(359)*V(188)*V(232)
  IRR(360) = RCT(360)*V(186)*V(232)
  IRR(361) = RCT(361)*V(146)*V(232)
  IRR(362) = RCT(362)*V(221)*V(232)
  IRR(363) = RCT(363)*V(212)*V(232)
  IRR(364) = RCT(364)*V(157)*V(232)
  IRR(365) = RCT(365)*V(204)*V(232)
  IRR(366) = RCT(366)*V(190)*V(232)
  IRR(367) = RCT(367)*V(183)*V(232)
  IRR(368) = RCT(368)*V(184)*V(232)
  IRR(369) = RCT(369)*V(214)*V(232)
  IRR(370) = RCT(370)*V(210)*V(232)
  IRR(371) = RCT(371)*V(168)*V(232)
  IRR(372) = RCT(372)*V(211)
  IRR(373) = RCT(373)*V(211)
  IRR(374) = RCT(374)*V(211)
  IRR(375) = RCT(375)*V(209)
  IRR(376) = RCT(376)*V(209)
  IRR(377) = RCT(377)*V(209)
  IRR(378) = RCT(378)*V(215)
  IRR(379) = RCT(379)*V(215)
  IRR(380) = RCT(380)*V(215)
  IRR(381) = RCT(381)*V(156)
  IRR(382) = RCT(382)*V(156)
  IRR(383) = RCT(383)*V(156)
  IRR(384) = RCT(384)*V(225)
  IRR(385) = RCT(385)*V(225)
  IRR(386) = RCT(386)*V(176)
  IRR(387) = RCT(387)*V(177)
  IRR(388) = RCT(388)*V(163)
  IRR(389) = RCT(389)*V(181)
  IRR(390) = RCT(390)*V(181)
  IRR(391) = RCT(391)*V(192)
  IRR(392) = RCT(392)*V(202)
  IRR(393) = RCT(393)*V(223)
  IRR(394) = RCT(394)*V(205)
  IRR(395) = RCT(395)*V(199)
  IRR(396) = RCT(396)*V(199)
  IRR(397) = RCT(397)*V(217)
  IRR(398) = RCT(398)*V(213)
  IRR(399) = RCT(399)*V(207)
  IRR(400) = RCT(400)*V(170)
  IRR(401) = RCT(401)*V(224)
  IRR(402) = RCT(402)*V(195)
  IRR(403) = RCT(403)*V(219)
  IRR(404) = RCT(404)*V(203)
  IRR(405) = RCT(405)*V(191)
  IRR(406) = RCT(406)*V(220)
  IRR(407) = RCT(407)*V(208)
  IRR(408) = RCT(408)*V(216)
  IRR(409) = RCT(409)*V(218)
  IRR(410) = RCT(410)*V(178)
  IRR(411) = RCT(411)*V(178)
  IRR(412) = RCT(412)*V(198)
  IRR(413) = RCT(413)*V(198)
  IRR(414) = RCT(414)*V(197)
  IRR(415) = RCT(415)*V(197)
  IRR(416) = RCT(416)*V(197)
  IRR(417) = RCT(417)*V(147)
  IRR(418) = RCT(418)*V(166)
  IRR(419) = RCT(419)*V(145)
  IRR(420) = RCT(420)*V(144)
  IRR(421) = RCT(421)*V(188)
  IRR(422) = RCT(422)*V(186)
  IRR(423) = RCT(423)*V(146)
  IRR(424) = RCT(424)*V(221)
  IRR(425) = RCT(425)*V(212)
  IRR(426) = RCT(426)*V(157)
  IRR(427) = RCT(427)*V(204)
  IRR(428) = RCT(428)*V(190)
  IRR(429) = RCT(429)*V(183)
  IRR(430) = RCT(430)*V(184)
  IRR(431) = RCT(431)*V(214)
  IRR(432) = RCT(432)*V(210)
  IRR(433) = RCT(433)*V(168)
  IRR(434) = RCT(434)*V(139)*V(228)
  IRR(435) = RCT(435)*V(131)*V(228)
  IRR(436) = RCT(436)*V(149)*V(228)
  IRR(437) = RCT(437)*V(179)*V(228)
  IRR(438) = RCT(438)*V(140)*V(228)
  IRR(439) = RCT(439)*V(169)*V(228)
  IRR(440) = RCT(440)*V(180)*V(228)
  IRR(441) = RCT(441)*V(201)*V(228)
  IRR(442) = RCT(442)*V(201)*V(230)
  IRR(443) = RCT(443)*V(201)*V(222)
  IRR(444) = RCT(444)*V(201)*V(222)
  IRR(445) = RCT(445)*V(196)*V(228)
  IRR(446) = RCT(446)*V(196)*V(230)
  IRR(447) = RCT(447)*V(150)*V(228)
  IRR(448) = RCT(448)*V(193)*V(228)
  IRR(449) = RCT(449)*V(148)*V(228)
  IRR(450) = RCT(450)*V(44)*V(228)
  IRR(451) = RCT(451)*V(37)*V(228)
  IRR(452) = RCT(452)*V(38)*V(228)
  IRR(453) = RCT(453)*V(194)*V(228)
  IRR(454) = RCT(454)*V(194)*V(230)
  IRR(455) = RCT(455)*V(194)*V(222)
  IRR(456) = RCT(456)*V(194)*V(222)
  IRR(457) = RCT(457)*V(151)*V(228)
  IRR(458) = RCT(458)*V(161)*V(228)
  IRR(459) = RCT(459)*V(154)*V(228)
  IRR(460) = RCT(460)*V(182)*V(228)
  IRR(461) = RCT(461)*V(167)*V(228)
  IRR(462) = RCT(462)*V(165)*V(228)
  IRR(463) = RCT(463)*V(33)*V(228)
  IRR(464) = RCT(464)*V(158)*V(228)
  IRR(465) = RCT(465)*V(165)*V(230)
  IRR(466) = RCT(466)*V(158)*V(230)
  IRR(467) = RCT(467)*V(39)*V(228)
  IRR(468) = RCT(468)*V(132)*V(228)
  IRR(469) = RCT(469)*V(107)*V(228)
  IRR(470) = RCT(470)*V(90)*V(228)
  IRR(471) = RCT(471)*V(92)*V(228)
  IRR(472) = RCT(472)*V(98)*V(228)
  IRR(473) = RCT(473)*V(85)*V(228)
  IRR(474) = RCT(474)*V(134)*V(228)
  IRR(475) = RCT(475)*V(79)*V(228)
  IRR(476) = RCT(476)*V(103)*V(228)
  IRR(477) = RCT(477)*V(65)*V(228)
  IRR(478) = RCT(478)*V(54)*V(228)
  IRR(479) = RCT(479)*V(55)*V(228)
  IRR(480) = RCT(480)*V(102)*V(228)
  IRR(481) = RCT(481)*V(78)*V(228)
  IRR(482) = RCT(482)*V(63)*V(228)
  IRR(483) = RCT(483)*V(118)*V(228)
  IRR(484) = RCT(484)*V(120)*V(228)
  IRR(485) = RCT(485)*V(138)*V(228)
  IRR(486) = RCT(486)*V(66)*V(228)
  IRR(487) = RCT(487)*V(95)*V(228)
  IRR(488) = RCT(488)*V(96)*V(228)
  IRR(489) = RCT(489)*V(97)*V(228)
  IRR(490) = RCT(490)*V(99)*V(228)
  IRR(491) = RCT(491)*V(46)*V(228)
  IRR(492) = RCT(492)*V(46)*V(230)
  IRR(493) = RCT(493)*V(47)*V(229)
  IRR(494) = RCT(494)*V(119)*V(228)
  IRR(495) = RCT(495)*V(119)*V(230)
  IRR(496) = RCT(496)*V(48)*V(228)
  IRR(497) = RCT(497)*V(48)*V(230)
  IRR(498) = RCT(498)*V(49)*V(229)
  IRR(499) = RCT(499)*V(121)*V(228)
  IRR(500) = RCT(500)*V(121)*V(230)
  IRR(501) = RCT(501)*V(101)*V(228)
  IRR(502) = RCT(502)*V(101)*V(228)
  IRR(503) = RCT(503)*V(73)*V(228)
  IRR(504) = RCT(504)*V(56)*V(228)
  IRR(505) = RCT(505)*V(75)*V(228)
  IRR(506) = RCT(506)*V(136)*V(228)
  IRR(507) = RCT(507)*V(80)*V(228)
  IRR(508) = RCT(508)*V(57)*V(228)
  IRR(509) = RCT(509)*V(53)*V(228)
  IRR(510) = RCT(510)*V(62)*V(228)
  IRR(511) = RCT(511)*V(84)*V(228)
  IRR(512) = RCT(512)*V(58)*V(228)
  IRR(513) = RCT(513)*V(59)*V(228)
  IRR(514) = RCT(514)*V(115)*V(228)
  IRR(515) = RCT(515)*V(81)*V(228)
  IRR(516) = RCT(516)*V(113)*V(228)
  IRR(517) = RCT(517)*V(114)*V(228)
  IRR(518) = RCT(518)*V(86)*V(228)
  IRR(519) = RCT(519)*V(64)*V(228)
  IRR(520) = RCT(520)*V(88)*V(228)
  IRR(521) = RCT(521)*V(105)*V(228)
  IRR(522) = RCT(522)*V(117)*V(228)
  IRR(523) = RCT(523)*V(82)*V(228)
  IRR(524) = RCT(524)*V(137)*V(228)
  IRR(525) = RCT(525)*V(129)*V(228)
  IRR(526) = RCT(526)*V(60)*V(228)
  IRR(527) = RCT(527)*V(83)*V(228)
  IRR(528) = RCT(528)*V(61)*V(228)
  IRR(529) = RCT(529)*V(89)*V(228)
  IRR(530) = RCT(530)*V(91)*V(228)
  IRR(531) = RCT(531)*V(130)*V(228)
  IRR(532) = RCT(532)*V(67)*V(228)
  IRR(533) = RCT(533)*V(68)*V(228)
  IRR(534) = RCT(534)*V(69)*V(228)
  IRR(535) = RCT(535)*V(70)*V(228)
  IRR(536) = RCT(536)*V(71)*V(228)
  IRR(537) = RCT(537)*V(93)*V(228)
  IRR(538) = RCT(538)*V(123)*V(228)
  IRR(539) = RCT(539)*V(94)*V(228)
  IRR(540) = RCT(540)*V(124)*V(228)
  IRR(541) = RCT(541)*V(72)*V(228)
  IRR(542) = RCT(542)*V(100)*V(228)
  IRR(543) = RCT(543)*V(125)*V(228)
  IRR(544) = RCT(544)*V(219)*V(229)
  IRR(545) = RCT(545)*V(126)
  IRR(546) = RCT(546)*V(203)*V(229)
  IRR(547) = RCT(547)*V(111)
  IRR(548) = RCT(548)*V(191)*V(229)
  IRR(549) = RCT(549)*V(112)
  IRR(550) = RCT(550)*V(126)*V(228)
  IRR(551) = RCT(551)*V(111)*V(228)
  IRR(552) = RCT(552)*V(112)*V(228)
  IRR(553) = RCT(553)*V(198)*V(229)
  IRR(554) = RCT(554)*V(87)
  IRR(555) = RCT(555)*V(197)*V(229)
  IRR(556) = RCT(556)*V(116)
  IRR(557) = RCT(557)*V(116)*V(228)
  IRR(558) = RCT(558)*V(87)*V(228)
  IRR(559) = RCT(559)*V(221)*V(229)
  IRR(560) = RCT(560)*V(122)
  IRR(561) = RCT(561)*V(122)*V(228)
  IRR(562) = RCT(562)*V(76)*V(228)
  IRR(563) = RCT(563)*V(154)*V(228)
  IRR(564) = RCT(564)*V(182)*V(228)
  IRR(565) = RCT(565)*V(167)*V(228)
  IRR(566) = RCT(566)*V(11)*V(228)
  IRR(567) = RCT(567)*V(12)*V(228)
  IRR(568) = RCT(568)*V(13)*V(228)
  IRR(569) = RCT(569)*V(14)*V(228)
  IRR(570) = RCT(570)*V(18)*V(228)
  IRR(571) = RCT(571)*V(17)*V(228)
  IRR(572) = RCT(572)*V(15)*V(228)
  IRR(573) = RCT(573)*V(16)*V(228)
  IRR(574) = RCT(574)*V(211)*V(229)
  IRR(575) = RCT(575)*V(77)
  IRR(576) = RCT(576)*V(135)
  IRR(577) = RCT(577)*V(141)
  IRR(578) = RCT(578)*V(133)
  IRR(579) = RCT(579)*V(128)
  IRR(580) = RCT(580)*V(104)
  IRR(581) = RCT(581)*V(152)
  IRR(582) = RCT(582)*V(152)
  IRR(583) = RCT(583)*V(204)*V(227)
  IRR(584) = RCT(584)*V(160)*V(228)
  IRR(585) = RCT(585)*V(77)*V(228)
  IRR(586) = RCT(586)*V(135)*V(228)
  IRR(587) = RCT(587)*V(135)*V(230)
  IRR(588) = RCT(588)*V(27)*V(228)
  IRR(589) = RCT(589)*V(27)*V(228)
  IRR(590) = RCT(590)*V(28)*V(228)
  IRR(591) = RCT(591)*V(28)*V(228)
  IRR(592) = RCT(592)*V(20)*V(228)
  IRR(593) = RCT(593)*V(21)*V(228)
  IRR(594) = RCT(594)*V(22)*V(228)
  IRR(595) = RCT(595)*V(23)*V(228)
  IRR(596) = RCT(596)*V(155)*V(227)
  IRR(597) = RCT(597)*V(200)*V(227)
  IRR(598) = RCT(598)*V(155)*V(227)
  IRR(599) = RCT(599)*V(200)*V(227)
  IRR(600) = RCT(600)*V(155)*V(230)
  IRR(601) = RCT(601)*V(200)*V(230)
  IRR(602) = RCT(602)*V(155)*V(232)
  IRR(603) = RCT(603)*V(200)*V(232)
  IRR(604) = RCT(604)*V(155)
  IRR(605) = RCT(605)*V(200)
  IRR(606) = RCT(606)*V(141)*V(228)
  IRR(607) = RCT(607)*V(133)*V(228)
  IRR(608) = RCT(608)*V(24)*V(228)
  IRR(609) = RCT(609)*V(174)*V(227)
  IRR(610) = RCT(610)*V(174)*V(227)
  IRR(611) = RCT(611)*V(174)*V(230)
  IRR(612) = RCT(612)*V(174)*V(232)
  IRR(613) = RCT(613)*V(174)
  IRR(614) = RCT(614)*V(128)*V(228)
  IRR(615) = RCT(615)*V(104)*V(228)
  IRR(616) = RCT(616)*V(152)*V(228)
  IRR(617) = RCT(617)*V(152)*V(228)
  IRR(620) = RCT(620)*V(108)*V(228)
  IRR(621) = RCT(621)*V(108)*V(228)
  IRR(622) = RCT(622)*V(108)*V(230)
  IRR(623) = RCT(623)*V(110)*V(227)
  IRR(624) = RCT(624)*V(110)*V(110)
  IRR(625) = RCT(625)*V(142)*V(222)
  IRR(626) = RCT(626)*V(142)*V(229)
  IRR(627) = RCT(627)*V(172)*V(229)
  IRR(628) = RCT(628)*V(172)*V(222)
  IRR(629) = RCT(629)*V(171)
  IRR(630) = RCT(630)*V(171)*V(229)
  IRR(631) = RCT(631)*V(171)*V(222)
  IRR(632) = RCT(632)*V(52)*V(232)
  IRR(633) = RCT(633)*V(52)
  IRR(634) = RCT(634)*V(41)*V(228)
  IRR(635) = RCT(635)*V(109)*V(228)
  IRR(636) = RCT(636)*V(109)*V(230)
      
END SUBROUTINE crimech_IRRFun
















SUBROUTINE crimech_Jac_SP ( V, F, RCT, JVS )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: JVS(LU_NONZERO)




  REAL(kind=dp) :: B(1089)


  B(1) = RCT(1)

  B(2) = RCT(2)

  B(3) = RCT(3)

  B(4) = RCT(4)

  B(5) = RCT(5)

  B(6) = RCT(6)

  B(7) = RCT(7)

  B(8) = RCT(8)

  B(9) = RCT(9)

  B(10) = RCT(10)

  B(11) = RCT(11)

  B(12) = RCT(12)

  B(13) = RCT(13)

  B(14) = RCT(14)

  B(15) = RCT(15)

  B(16) = RCT(16)

  B(17) = RCT(17)

  B(18) = RCT(18)

  B(19) = RCT(19)

  B(20) = RCT(20)

  B(21) = RCT(21)

  B(22) = RCT(22)

  B(23) = RCT(23)

  B(24) = RCT(24)

  B(25) = RCT(25)

  B(26) = RCT(26)

  B(27) = RCT(27)

  B(28) = RCT(28)

  B(29) = RCT(29)

  B(30) = RCT(30)

  B(31) = RCT(31)

  B(32) = RCT(32)

  B(33) = RCT(33)

  B(34) = RCT(34)

  B(35) = RCT(35)

  B(36) = RCT(36)

  B(37) = RCT(37)

  B(38) = RCT(38)

  B(39) = RCT(39)

  B(40) = RCT(40)

  B(41) = RCT(41)

  B(42) = RCT(42)

  B(43) = RCT(43)

  B(44) = RCT(44)

  B(45) = RCT(45)

  B(46) = RCT(46)

  B(47) = RCT(47)

  B(48) = RCT(48)

  B(49) = RCT(49)

  B(50) = RCT(50)

  B(51) = RCT(51)

  B(52) = RCT(52)

  B(53) = RCT(53)

  B(54) = RCT(54)

  B(55) = RCT(55)

  B(56) = RCT(56)

  B(57) = RCT(57)

  B(58) = RCT(58)

  B(59) = RCT(59)

  B(60) = RCT(60)

  B(61) = RCT(61)

  B(62) = RCT(62)

  B(63) = RCT(63)

  B(64) = RCT(64)

  B(65) = RCT(65)

  B(66) = RCT(66)

  B(67) = RCT(67)

  B(68) = RCT(68)

  B(69) = RCT(69)

  B(70) = RCT(70)

  B(71) = RCT(71)

  B(72) = RCT(72)

  B(73) = RCT(73)

  B(74) = RCT(74)

  B(75) = RCT(75)

  B(76) = RCT(76)

  B(77) = RCT(77)

  B(78) = RCT(78)

  B(79) = RCT(79)

  B(80) = RCT(80)

  B(81) = RCT(81)

  B(82) = RCT(82)

  B(83) = RCT(83)

  B(84) = RCT(84)

  B(85) = RCT(85)

  B(86) = RCT(86)

  B(87) = RCT(87)

  B(88) = RCT(88)

  B(89) = RCT(89)

  B(90) = RCT(90)

  B(91) = RCT(91)

  B(92) = RCT(92)

  B(93) = RCT(93)

  B(94) = RCT(94)

  B(95) = RCT(95)

  B(96) = RCT(96)*F(2)

  B(98) = RCT(97)

  B(99) = RCT(98)*V(222)

  B(100) = RCT(98)*V(173)

  B(101) = RCT(99)*V(227)

  B(102) = RCT(99)*V(173)

  B(103) = RCT(100)*V(229)

  B(104) = RCT(100)*V(173)

  B(105) = RCT(101)*V(229)

  B(106) = RCT(101)*V(173)

  B(107) = RCT(102)*F(2)

  B(109) = RCT(103)*V(227)

  B(110) = RCT(103)*V(222)

  B(111) = RCT(104)*V(229)

  B(112) = RCT(104)*V(222)

  B(113) = RCT(105)*2*V(227)*F(2)

  B(115) = RCT(106)*V(230)

  B(116) = RCT(106)*V(227)

  B(117) = RCT(107)*V(230)

  B(118) = RCT(107)*V(229)

  B(119) = RCT(108)*V(230)

  B(120) = RCT(108)*V(229)

  B(121) = RCT(109)

  B(122) = RCT(110)*F(1)

  B(124) = RCT(111)*V(228)

  B(125) = RCT(111)*V(222)

  B(126) = RCT(112)*F(2)

  B(128) = RCT(113)*V(228)

  B(129) = RCT(113)*V(189)

  B(130) = RCT(114)*V(228)

  B(131) = RCT(114)*V(106)

  B(132) = RCT(115)*V(232)

  B(133) = RCT(115)*V(222)

  B(134) = RCT(116)*V(232)

  B(135) = RCT(116)*V(228)

  B(136) = RCT(117)*2*V(232)

  B(137) = RCT(118)*2*V(232)*F(1)

  B(139) = RCT(119)*V(228)

  B(140) = RCT(119)*V(227)

  B(141) = RCT(120)*V(229)

  B(142) = RCT(120)*V(228)

  B(143) = RCT(121)*V(230)

  B(144) = RCT(121)*V(228)

  B(145) = RCT(122)*V(232)

  B(146) = RCT(122)*V(227)

  B(147) = RCT(123)*V(232)

  B(148) = RCT(123)*V(229)

  B(149) = RCT(124)

  B(150) = RCT(125)*V(228)

  B(151) = RCT(125)*V(74)

  B(152) = RCT(126)*V(232)

  B(153) = RCT(126)*V(230)

  B(154) = RCT(127)*V(228)

  B(155) = RCT(127)*V(40)

  B(156) = RCT(128)*V(228)

  B(157) = RCT(128)*V(185)

  B(158) = RCT(129)*V(173)

  B(159) = RCT(129)*V(127)

  B(160) = RCT(130)*V(228)

  B(161) = RCT(130)*V(127)

  B(162) = RCT(131)*F(2)

  B(164) = RCT(132)*F(1)*F(1)

  B(166) = RCT(133)*V(228)

  B(167) = RCT(133)*V(19)

  B(168) = RCT(134)*V(228)

  B(169) = RCT(134)*V(8)

  B(170) = RCT(135)*V(228)

  B(171) = RCT(135)*V(25)

  B(172) = RCT(136)*V(228)

  B(173) = RCT(136)*V(25)

  B(174) = RCT(137)*V(228)

  B(175) = RCT(137)*V(10)

  B(176) = RCT(138)*V(228)

  B(177) = RCT(138)*V(153)

  B(178) = RCT(139)*V(228)

  B(179) = RCT(139)*V(143)

  B(180) = RCT(140)*V(228)

  B(181) = RCT(140)*V(162)

  B(182) = RCT(141)*V(230)

  B(183) = RCT(141)*V(153)

  B(184) = RCT(142)*V(230)

  B(185) = RCT(142)*V(143)

  B(186) = RCT(143)*V(230)

  B(187) = RCT(143)*V(162)

  B(188) = RCT(144)*V(222)

  B(189) = RCT(144)*V(153)

  B(190) = RCT(145)*V(222)

  B(191) = RCT(145)*V(153)

  B(192) = RCT(146)*V(222)

  B(193) = RCT(146)*V(143)

  B(194) = RCT(147)*V(222)

  B(195) = RCT(147)*V(143)

  B(196) = RCT(148)*V(222)

  B(197) = RCT(148)*V(162)

  B(198) = RCT(149)*V(222)

  B(199) = RCT(149)*V(162)

  B(200) = RCT(150)*V(228)

  B(201) = RCT(150)*V(175)

  B(202) = RCT(151)*V(230)

  B(203) = RCT(151)*V(175)

  B(204) = RCT(152)*V(222)

  B(205) = RCT(152)*V(175)

  B(206) = RCT(153)*V(222)

  B(207) = RCT(153)*V(175)

  B(208) = RCT(154)*V(228)

  B(209) = RCT(154)*V(164)

  B(210) = RCT(155)*V(230)

  B(211) = RCT(155)*V(164)

  B(212) = RCT(156)*V(222)

  B(213) = RCT(156)*V(164)

  B(214) = RCT(157)*V(222)

  B(215) = RCT(157)*V(164)

  B(216) = RCT(158)*V(222)

  B(217) = RCT(158)*V(164)

  B(218) = RCT(159)*V(228)

  B(219) = RCT(159)*V(159)

  B(220) = RCT(160)*V(230)

  B(221) = RCT(160)*V(159)

  B(222) = RCT(161)*V(222)

  B(223) = RCT(161)*V(159)

  B(224) = RCT(162)*V(222)

  B(225) = RCT(162)*V(159)

  B(226) = RCT(163)*V(222)

  B(227) = RCT(163)*V(159)

  B(228) = RCT(164)*V(222)

  B(229) = RCT(164)*V(159)

  B(230) = RCT(165)*V(228)

  B(231) = RCT(165)*V(35)

  B(232) = RCT(166)*V(228)

  B(233) = RCT(166)*V(35)

  B(234) = RCT(167)*V(228)

  B(235) = RCT(167)*V(29)

  B(236) = RCT(168)*V(228)

  B(237) = RCT(168)*V(29)

  B(238) = RCT(169)*V(228)

  B(239) = RCT(169)*V(30)

  B(240) = RCT(170)*V(228)

  B(241) = RCT(170)*V(30)

  B(242) = RCT(171)*V(228)

  B(243) = RCT(171)*V(26)

  B(244) = RCT(172)*V(228)

  B(245) = RCT(172)*V(26)

  B(246) = RCT(173)*V(228)

  B(247) = RCT(173)*V(206)

  B(248) = RCT(174)*V(231)

  B(249) = RCT(174)*V(228)

  B(250) = RCT(175)*V(228)

  B(251) = RCT(175)*V(226)

  B(252) = RCT(176)*V(230)

  B(253) = RCT(176)*V(206)

  B(254) = RCT(177)*V(231)

  B(255) = RCT(177)*V(230)

  B(256) = RCT(178)*V(230)

  B(257) = RCT(178)*V(226)

  B(258) = RCT(179)*V(228)

  B(259) = RCT(179)*V(187)

  B(260) = RCT(180)*V(228)

  B(261) = RCT(180)*V(31)

  B(262) = RCT(181)*V(228)

  B(263) = RCT(181)*V(36)

  B(264) = RCT(182)*V(228)

  B(265) = RCT(182)*V(42)

  B(266) = RCT(183)*V(228)

  B(267) = RCT(183)*V(42)

  B(268) = RCT(184)*V(228)

  B(269) = RCT(184)*V(50)

  B(270) = RCT(185)*V(228)

  B(271) = RCT(185)*V(50)

  B(272) = RCT(186)*V(228)

  B(273) = RCT(186)*V(43)

  B(274) = RCT(187)*V(228)

  B(275) = RCT(187)*V(43)

  B(276) = RCT(188)*V(228)

  B(277) = RCT(188)*V(45)

  B(278) = RCT(189)*V(228)

  B(279) = RCT(189)*V(51)

  B(280) = RCT(190)*V(227)

  B(281) = RCT(190)*V(211)

  B(282) = RCT(191)*V(227)

  B(283) = RCT(191)*V(209)

  B(284) = RCT(192)*V(227)

  B(285) = RCT(192)*V(215)

  B(286) = RCT(193)*V(227)

  B(287) = RCT(193)*V(156)

  B(288) = RCT(194)*V(227)

  B(289) = RCT(194)*V(225)

  B(290) = RCT(195)*V(227)

  B(291) = RCT(195)*V(225)

  B(292) = RCT(196)*V(227)

  B(293) = RCT(196)*V(223)

  B(294) = RCT(197)*V(227)

  B(295) = RCT(197)*V(205)

  B(296) = RCT(198)*V(227)

  B(297) = RCT(198)*V(176)

  B(298) = RCT(199)*V(227)

  B(299) = RCT(199)*V(177)

  B(300) = RCT(200)*V(227)

  B(301) = RCT(200)*V(163)

  B(302) = RCT(201)*V(227)

  B(303) = RCT(201)*V(181)

  B(304) = RCT(202)*V(227)

  B(305) = RCT(202)*V(181)

  B(306) = RCT(203)*V(227)

  B(307) = RCT(203)*V(192)

  B(308) = RCT(204)*V(227)

  B(309) = RCT(204)*V(202)

  B(310) = RCT(205)*V(227)

  B(311) = RCT(205)*V(199)

  B(312) = RCT(206)*V(227)

  B(313) = RCT(206)*V(199)

  B(314) = RCT(207)*V(227)

  B(315) = RCT(207)*V(217)

  B(316) = RCT(208)*V(227)

  B(317) = RCT(208)*V(213)

  B(318) = RCT(209)*V(227)

  B(319) = RCT(209)*V(207)

  B(320) = RCT(210)*V(227)

  B(321) = RCT(210)*V(170)

  B(322) = RCT(211)*V(227)

  B(323) = RCT(211)*V(224)

  B(324) = RCT(212)*V(227)

  B(325) = RCT(212)*V(195)

  B(326) = RCT(213)*V(227)

  B(327) = RCT(213)*V(219)

  B(328) = RCT(214)*V(227)

  B(329) = RCT(214)*V(203)

  B(330) = RCT(215)*V(227)

  B(331) = RCT(215)*V(191)

  B(332) = RCT(216)*V(227)

  B(333) = RCT(216)*V(220)

  B(334) = RCT(217)*V(227)

  B(335) = RCT(217)*V(208)

  B(336) = RCT(218)*V(227)

  B(337) = RCT(218)*V(216)

  B(338) = RCT(219)*V(227)

  B(339) = RCT(219)*V(218)

  B(340) = RCT(220)*V(227)

  B(341) = RCT(220)*V(178)

  B(342) = RCT(221)*V(227)

  B(343) = RCT(221)*V(178)

  B(344) = RCT(222)*V(227)

  B(345) = RCT(222)*V(198)

  B(346) = RCT(223)*V(227)

  B(347) = RCT(223)*V(198)

  B(348) = RCT(224)*V(227)

  B(349) = RCT(224)*V(197)

  B(350) = RCT(225)*V(227)

  B(351) = RCT(225)*V(197)

  B(352) = RCT(226)*V(227)

  B(353) = RCT(226)*V(197)

  B(354) = RCT(227)*V(227)

  B(355) = RCT(227)*V(147)

  B(356) = RCT(228)*V(227)

  B(357) = RCT(228)*V(166)

  B(358) = RCT(229)*V(227)

  B(359) = RCT(229)*V(145)

  B(360) = RCT(230)*V(227)

  B(361) = RCT(230)*V(144)

  B(362) = RCT(231)*V(227)

  B(363) = RCT(231)*V(188)

  B(364) = RCT(232)*V(227)

  B(365) = RCT(232)*V(186)

  B(366) = RCT(233)*V(227)

  B(367) = RCT(233)*V(186)

  B(368) = RCT(234)*V(227)

  B(369) = RCT(234)*V(146)

  B(370) = RCT(235)*V(227)

  B(371) = RCT(235)*V(221)

  B(372) = RCT(236)*V(227)

  B(373) = RCT(236)*V(212)

  B(374) = RCT(237)*V(227)

  B(375) = RCT(237)*V(157)

  B(376) = RCT(238)*V(227)

  B(377) = RCT(238)*V(204)

  B(378) = RCT(239)*V(227)

  B(379) = RCT(239)*V(190)

  B(380) = RCT(240)*V(227)

  B(381) = RCT(240)*V(183)

  B(382) = RCT(241)*V(227)

  B(383) = RCT(241)*V(184)

  B(384) = RCT(242)*V(227)

  B(385) = RCT(242)*V(184)

  B(386) = RCT(243)*V(227)

  B(387) = RCT(243)*V(168)

  B(388) = RCT(244)*V(227)

  B(389) = RCT(244)*V(214)

  B(390) = RCT(245)*V(227)

  B(391) = RCT(245)*V(214)

  B(392) = RCT(246)*V(227)

  B(393) = RCT(246)*V(210)

  B(394) = RCT(247)*V(227)

  B(395) = RCT(247)*V(211)

  B(396) = RCT(248)*V(227)

  B(397) = RCT(248)*V(209)

  B(398) = RCT(249)*V(227)

  B(399) = RCT(249)*V(215)

  B(400) = RCT(250)*V(227)

  B(401) = RCT(250)*V(156)

  B(402) = RCT(251)*V(227)

  B(403) = RCT(251)*V(225)

  B(404) = RCT(252)*V(227)

  B(405) = RCT(252)*V(223)

  B(406) = RCT(253)*V(227)

  B(407) = RCT(253)*V(205)

  B(408) = RCT(254)*V(227)

  B(409) = RCT(254)*V(199)

  B(410) = RCT(255)*V(227)

  B(411) = RCT(255)*V(217)

  B(412) = RCT(256)*V(227)

  B(413) = RCT(256)*V(213)

  B(414) = RCT(257)*V(227)

  B(415) = RCT(257)*V(207)

  B(416) = RCT(258)*V(227)

  B(417) = RCT(258)*V(170)

  B(418) = RCT(259)*V(227)

  B(419) = RCT(259)*V(224)

  B(420) = RCT(260)*V(227)

  B(421) = RCT(260)*V(195)

  B(422) = RCT(261)*V(227)

  B(423) = RCT(261)*V(178)

  B(424) = RCT(262)*V(227)

  B(425) = RCT(262)*V(163)

  B(426) = RCT(263)*V(227)

  B(427) = RCT(263)*V(181)

  B(428) = RCT(264)*V(227)

  B(429) = RCT(264)*V(192)

  B(430) = RCT(265)*V(227)

  B(431) = RCT(265)*V(202)

  B(432) = RCT(266)*V(227)

  B(433) = RCT(266)*V(186)

  B(434) = RCT(267)*V(227)

  B(435) = RCT(267)*V(212)

  B(436) = RCT(268)*V(227)

  B(437) = RCT(268)*V(184)

  B(438) = RCT(269)*V(227)

  B(439) = RCT(269)*V(214)

  B(440) = RCT(270)*V(227)

  B(441) = RCT(270)*V(210)

  B(442) = RCT(271)*V(230)

  B(443) = RCT(271)*V(211)

  B(444) = RCT(272)*V(230)

  B(445) = RCT(272)*V(209)

  B(446) = RCT(273)*V(230)

  B(447) = RCT(273)*V(215)

  B(448) = RCT(274)*V(230)

  B(449) = RCT(274)*V(156)

  B(450) = RCT(275)*V(230)

  B(451) = RCT(275)*V(225)

  B(452) = RCT(276)*V(230)

  B(453) = RCT(276)*V(225)

  B(454) = RCT(277)*V(230)

  B(455) = RCT(277)*V(223)

  B(456) = RCT(278)*V(230)

  B(457) = RCT(278)*V(205)

  B(458) = RCT(279)*V(230)

  B(459) = RCT(279)*V(176)

  B(460) = RCT(280)*V(230)

  B(461) = RCT(280)*V(177)

  B(462) = RCT(281)*V(230)

  B(463) = RCT(281)*V(163)

  B(464) = RCT(282)*V(230)

  B(465) = RCT(282)*V(181)

  B(466) = RCT(283)*V(230)

  B(467) = RCT(283)*V(181)

  B(468) = RCT(284)*V(230)

  B(469) = RCT(284)*V(192)

  B(470) = RCT(285)*V(230)

  B(471) = RCT(285)*V(202)

  B(472) = RCT(286)*V(230)

  B(473) = RCT(286)*V(199)

  B(474) = RCT(287)*V(230)

  B(475) = RCT(287)*V(199)

  B(476) = RCT(288)*V(230)

  B(477) = RCT(288)*V(217)

  B(478) = RCT(289)*V(230)

  B(479) = RCT(289)*V(213)

  B(480) = RCT(290)*V(230)

  B(481) = RCT(290)*V(207)

  B(482) = RCT(291)*V(230)

  B(483) = RCT(291)*V(170)

  B(484) = RCT(292)*V(230)

  B(485) = RCT(292)*V(224)

  B(486) = RCT(293)*V(230)

  B(487) = RCT(293)*V(195)

  B(488) = RCT(294)*V(230)

  B(489) = RCT(294)*V(219)

  B(490) = RCT(295)*V(230)

  B(491) = RCT(295)*V(203)

  B(492) = RCT(296)*V(230)

  B(493) = RCT(296)*V(191)

  B(494) = RCT(297)*V(230)

  B(495) = RCT(297)*V(220)

  B(496) = RCT(298)*V(230)

  B(497) = RCT(298)*V(208)

  B(498) = RCT(299)*V(230)

  B(499) = RCT(299)*V(216)

  B(500) = RCT(300)*V(230)

  B(501) = RCT(300)*V(218)

  B(502) = RCT(301)*V(230)

  B(503) = RCT(301)*V(178)

  B(504) = RCT(302)*V(230)

  B(505) = RCT(302)*V(178)

  B(506) = RCT(303)*V(230)

  B(507) = RCT(303)*V(198)

  B(508) = RCT(304)*V(230)

  B(509) = RCT(304)*V(198)

  B(510) = RCT(305)*V(230)

  B(511) = RCT(305)*V(197)

  B(512) = RCT(306)*V(230)

  B(513) = RCT(306)*V(197)

  B(514) = RCT(307)*V(230)

  B(515) = RCT(307)*V(197)

  B(516) = RCT(308)*V(230)

  B(517) = RCT(308)*V(147)

  B(518) = RCT(309)*V(230)

  B(519) = RCT(309)*V(166)

  B(520) = RCT(310)*V(230)

  B(521) = RCT(310)*V(145)

  B(522) = RCT(311)*V(230)

  B(523) = RCT(311)*V(144)

  B(524) = RCT(312)*V(230)

  B(525) = RCT(312)*V(188)

  B(526) = RCT(313)*V(230)

  B(527) = RCT(313)*V(186)

  B(528) = RCT(314)*V(230)

  B(529) = RCT(314)*V(146)

  B(530) = RCT(315)*V(230)

  B(531) = RCT(315)*V(221)

  B(532) = RCT(316)*V(230)

  B(533) = RCT(316)*V(212)

  B(534) = RCT(317)*V(230)

  B(535) = RCT(317)*V(157)

  B(536) = RCT(318)*V(230)

  B(537) = RCT(318)*V(204)

  B(538) = RCT(319)*V(230)

  B(539) = RCT(319)*V(190)

  B(540) = RCT(320)*V(230)

  B(541) = RCT(320)*V(183)

  B(542) = RCT(321)*V(230)

  B(543) = RCT(321)*V(184)

  B(544) = RCT(322)*V(230)

  B(545) = RCT(322)*V(214)

  B(546) = RCT(323)*V(230)

  B(547) = RCT(323)*V(210)

  B(548) = RCT(324)*V(230)

  B(549) = RCT(324)*V(168)

  B(550) = RCT(325)*V(232)

  B(551) = RCT(325)*V(211)

  B(552) = RCT(326)*V(232)

  B(553) = RCT(326)*V(209)

  B(554) = RCT(327)*V(232)

  B(555) = RCT(327)*V(215)

  B(556) = RCT(328)*V(232)

  B(557) = RCT(328)*V(156)

  B(558) = RCT(329)*V(232)

  B(559) = RCT(329)*V(225)

  B(560) = RCT(330)*V(232)

  B(561) = RCT(330)*V(223)

  B(562) = RCT(331)*V(232)

  B(563) = RCT(331)*V(205)

  B(564) = RCT(332)*V(232)

  B(565) = RCT(332)*V(176)

  B(566) = RCT(333)*V(232)

  B(567) = RCT(333)*V(177)

  B(568) = RCT(334)*V(232)

  B(569) = RCT(334)*V(163)

  B(570) = RCT(335)*V(232)

  B(571) = RCT(335)*V(181)

  B(572) = RCT(336)*V(232)

  B(573) = RCT(336)*V(192)

  B(574) = RCT(337)*V(232)

  B(575) = RCT(337)*V(202)

  B(576) = RCT(338)*V(232)

  B(577) = RCT(338)*V(199)

  B(578) = RCT(339)*V(232)

  B(579) = RCT(339)*V(217)

  B(580) = RCT(340)*V(232)

  B(581) = RCT(340)*V(213)

  B(582) = RCT(341)*V(232)

  B(583) = RCT(341)*V(207)

  B(584) = RCT(342)*V(232)

  B(585) = RCT(342)*V(170)

  B(586) = RCT(343)*V(232)

  B(587) = RCT(343)*V(224)

  B(588) = RCT(344)*V(232)

  B(589) = RCT(344)*V(195)

  B(590) = RCT(345)*V(232)

  B(591) = RCT(345)*V(219)

  B(592) = RCT(346)*V(232)

  B(593) = RCT(346)*V(203)

  B(594) = RCT(347)*V(232)

  B(595) = RCT(347)*V(191)

  B(596) = RCT(348)*V(232)

  B(597) = RCT(348)*V(220)

  B(598) = RCT(349)*V(232)

  B(599) = RCT(349)*V(208)

  B(600) = RCT(350)*V(232)

  B(601) = RCT(350)*V(216)

  B(602) = RCT(351)*V(232)

  B(603) = RCT(351)*V(218)

  B(604) = RCT(352)*V(232)

  B(605) = RCT(352)*V(178)

  B(606) = RCT(353)*V(232)

  B(607) = RCT(353)*V(198)

  B(608) = RCT(354)*V(232)

  B(609) = RCT(354)*V(197)

  B(610) = RCT(355)*V(232)

  B(611) = RCT(355)*V(147)

  B(612) = RCT(356)*V(232)

  B(613) = RCT(356)*V(166)

  B(614) = RCT(357)*V(232)

  B(615) = RCT(357)*V(145)

  B(616) = RCT(358)*V(232)

  B(617) = RCT(358)*V(144)

  B(618) = RCT(359)*V(232)

  B(619) = RCT(359)*V(188)

  B(620) = RCT(360)*V(232)

  B(621) = RCT(360)*V(186)

  B(622) = RCT(361)*V(232)

  B(623) = RCT(361)*V(146)

  B(624) = RCT(362)*V(232)

  B(625) = RCT(362)*V(221)

  B(626) = RCT(363)*V(232)

  B(627) = RCT(363)*V(212)

  B(628) = RCT(364)*V(232)

  B(629) = RCT(364)*V(157)

  B(630) = RCT(365)*V(232)

  B(631) = RCT(365)*V(204)

  B(632) = RCT(366)*V(232)

  B(633) = RCT(366)*V(190)

  B(634) = RCT(367)*V(232)

  B(635) = RCT(367)*V(183)

  B(636) = RCT(368)*V(232)

  B(637) = RCT(368)*V(184)

  B(638) = RCT(369)*V(232)

  B(639) = RCT(369)*V(214)

  B(640) = RCT(370)*V(232)

  B(641) = RCT(370)*V(210)

  B(642) = RCT(371)*V(232)

  B(643) = RCT(371)*V(168)

  B(644) = RCT(372)

  B(645) = RCT(373)

  B(646) = RCT(374)

  B(647) = RCT(375)

  B(648) = RCT(376)

  B(649) = RCT(377)

  B(650) = RCT(378)

  B(651) = RCT(379)

  B(652) = RCT(380)

  B(653) = RCT(381)

  B(654) = RCT(382)

  B(655) = RCT(383)

  B(656) = RCT(384)

  B(657) = RCT(385)

  B(658) = RCT(386)

  B(659) = RCT(387)

  B(660) = RCT(388)

  B(661) = RCT(389)

  B(662) = RCT(390)

  B(663) = RCT(391)

  B(664) = RCT(392)

  B(665) = RCT(393)

  B(666) = RCT(394)

  B(667) = RCT(395)

  B(668) = RCT(396)

  B(669) = RCT(397)

  B(670) = RCT(398)

  B(671) = RCT(399)

  B(672) = RCT(400)

  B(673) = RCT(401)

  B(674) = RCT(402)

  B(675) = RCT(403)

  B(676) = RCT(404)

  B(677) = RCT(405)

  B(678) = RCT(406)

  B(679) = RCT(407)

  B(680) = RCT(408)

  B(681) = RCT(409)

  B(682) = RCT(410)

  B(683) = RCT(411)

  B(684) = RCT(412)

  B(685) = RCT(413)

  B(686) = RCT(414)

  B(687) = RCT(415)

  B(688) = RCT(416)

  B(689) = RCT(417)

  B(690) = RCT(418)

  B(691) = RCT(419)

  B(692) = RCT(420)

  B(693) = RCT(421)

  B(694) = RCT(422)

  B(695) = RCT(423)

  B(696) = RCT(424)

  B(697) = RCT(425)

  B(698) = RCT(426)

  B(699) = RCT(427)

  B(700) = RCT(428)

  B(701) = RCT(429)

  B(702) = RCT(430)

  B(703) = RCT(431)

  B(704) = RCT(432)

  B(705) = RCT(433)

  B(706) = RCT(434)*V(228)

  B(707) = RCT(434)*V(139)

  B(708) = RCT(435)*V(228)

  B(709) = RCT(435)*V(131)

  B(710) = RCT(436)*V(228)

  B(711) = RCT(436)*V(149)

  B(712) = RCT(437)*V(228)

  B(713) = RCT(437)*V(179)

  B(714) = RCT(438)*V(228)

  B(715) = RCT(438)*V(140)

  B(716) = RCT(439)*V(228)

  B(717) = RCT(439)*V(169)

  B(718) = RCT(440)*V(228)

  B(719) = RCT(440)*V(180)

  B(720) = RCT(441)*V(228)

  B(721) = RCT(441)*V(201)

  B(722) = RCT(442)*V(230)

  B(723) = RCT(442)*V(201)

  B(724) = RCT(443)*V(222)

  B(725) = RCT(443)*V(201)

  B(726) = RCT(444)*V(222)

  B(727) = RCT(444)*V(201)

  B(728) = RCT(445)*V(228)

  B(729) = RCT(445)*V(196)

  B(730) = RCT(446)*V(230)

  B(731) = RCT(446)*V(196)

  B(732) = RCT(447)*V(228)

  B(733) = RCT(447)*V(150)

  B(734) = RCT(448)*V(228)

  B(735) = RCT(448)*V(193)

  B(736) = RCT(449)*V(228)

  B(737) = RCT(449)*V(148)

  B(738) = RCT(450)*V(228)

  B(739) = RCT(450)*V(44)

  B(740) = RCT(451)*V(228)

  B(741) = RCT(451)*V(37)

  B(742) = RCT(452)*V(228)

  B(743) = RCT(452)*V(38)

  B(744) = RCT(453)*V(228)

  B(745) = RCT(453)*V(194)

  B(746) = RCT(454)*V(230)

  B(747) = RCT(454)*V(194)

  B(748) = RCT(455)*V(222)

  B(749) = RCT(455)*V(194)

  B(750) = RCT(456)*V(222)

  B(751) = RCT(456)*V(194)

  B(752) = RCT(457)*V(228)

  B(753) = RCT(457)*V(151)

  B(754) = RCT(458)*V(228)

  B(755) = RCT(458)*V(161)

  B(756) = RCT(459)*V(228)

  B(757) = RCT(459)*V(154)

  B(758) = RCT(460)*V(228)

  B(759) = RCT(460)*V(182)

  B(760) = RCT(461)*V(228)

  B(761) = RCT(461)*V(167)

  B(762) = RCT(462)*V(228)

  B(763) = RCT(462)*V(165)

  B(764) = RCT(463)*V(228)

  B(765) = RCT(463)*V(33)

  B(766) = RCT(464)*V(228)

  B(767) = RCT(464)*V(158)

  B(768) = RCT(465)*V(230)

  B(769) = RCT(465)*V(165)

  B(770) = RCT(466)*V(230)

  B(771) = RCT(466)*V(158)

  B(772) = RCT(467)*V(228)

  B(773) = RCT(467)*V(39)

  B(774) = RCT(468)*V(228)

  B(775) = RCT(468)*V(132)

  B(776) = RCT(469)*V(228)

  B(777) = RCT(469)*V(107)

  B(778) = RCT(470)*V(228)

  B(779) = RCT(470)*V(90)

  B(780) = RCT(471)*V(228)

  B(781) = RCT(471)*V(92)

  B(782) = RCT(472)*V(228)

  B(783) = RCT(472)*V(98)

  B(784) = RCT(473)*V(228)

  B(785) = RCT(473)*V(85)

  B(786) = RCT(474)*V(228)

  B(787) = RCT(474)*V(134)

  B(788) = RCT(475)*V(228)

  B(789) = RCT(475)*V(79)

  B(790) = RCT(476)*V(228)

  B(791) = RCT(476)*V(103)

  B(792) = RCT(477)*V(228)

  B(793) = RCT(477)*V(65)

  B(794) = RCT(478)*V(228)

  B(795) = RCT(478)*V(54)

  B(796) = RCT(479)*V(228)

  B(797) = RCT(479)*V(55)

  B(798) = RCT(480)*V(228)

  B(799) = RCT(480)*V(102)

  B(800) = RCT(481)*V(228)

  B(801) = RCT(481)*V(78)

  B(802) = RCT(482)*V(228)

  B(803) = RCT(482)*V(63)

  B(804) = RCT(483)*V(228)

  B(805) = RCT(483)*V(118)

  B(806) = RCT(484)*V(228)

  B(807) = RCT(484)*V(120)

  B(808) = RCT(485)*V(228)

  B(809) = RCT(485)*V(138)

  B(810) = RCT(486)*V(228)

  B(811) = RCT(486)*V(66)

  B(812) = RCT(487)*V(228)

  B(813) = RCT(487)*V(95)

  B(814) = RCT(488)*V(228)

  B(815) = RCT(488)*V(96)

  B(816) = RCT(489)*V(228)

  B(817) = RCT(489)*V(97)

  B(818) = RCT(490)*V(228)

  B(819) = RCT(490)*V(99)

  B(820) = RCT(491)*V(228)

  B(821) = RCT(491)*V(46)

  B(822) = RCT(492)*V(230)

  B(823) = RCT(492)*V(46)

  B(824) = RCT(493)*V(229)

  B(825) = RCT(493)*V(47)

  B(826) = RCT(494)*V(228)

  B(827) = RCT(494)*V(119)

  B(828) = RCT(495)*V(230)

  B(829) = RCT(495)*V(119)

  B(830) = RCT(496)*V(228)

  B(831) = RCT(496)*V(48)

  B(832) = RCT(497)*V(230)

  B(833) = RCT(497)*V(48)

  B(834) = RCT(498)*V(229)

  B(835) = RCT(498)*V(49)

  B(836) = RCT(499)*V(228)

  B(837) = RCT(499)*V(121)

  B(838) = RCT(500)*V(230)

  B(839) = RCT(500)*V(121)

  B(840) = RCT(501)*V(228)

  B(841) = RCT(501)*V(101)

  B(842) = RCT(502)*V(228)

  B(843) = RCT(502)*V(101)

  B(844) = RCT(503)*V(228)

  B(845) = RCT(503)*V(73)

  B(846) = RCT(504)*V(228)

  B(847) = RCT(504)*V(56)

  B(848) = RCT(505)*V(228)

  B(849) = RCT(505)*V(75)

  B(850) = RCT(506)*V(228)

  B(851) = RCT(506)*V(136)

  B(852) = RCT(507)*V(228)

  B(853) = RCT(507)*V(80)

  B(854) = RCT(508)*V(228)

  B(855) = RCT(508)*V(57)

  B(856) = RCT(509)*V(228)

  B(857) = RCT(509)*V(53)

  B(858) = RCT(510)*V(228)

  B(859) = RCT(510)*V(62)

  B(860) = RCT(511)*V(228)

  B(861) = RCT(511)*V(84)

  B(862) = RCT(512)*V(228)

  B(863) = RCT(512)*V(58)

  B(864) = RCT(513)*V(228)

  B(865) = RCT(513)*V(59)

  B(866) = RCT(514)*V(228)

  B(867) = RCT(514)*V(115)

  B(868) = RCT(515)*V(228)

  B(869) = RCT(515)*V(81)

  B(870) = RCT(516)*V(228)

  B(871) = RCT(516)*V(113)

  B(872) = RCT(517)*V(228)

  B(873) = RCT(517)*V(114)

  B(874) = RCT(518)*V(228)

  B(875) = RCT(518)*V(86)

  B(876) = RCT(519)*V(228)

  B(877) = RCT(519)*V(64)

  B(878) = RCT(520)*V(228)

  B(879) = RCT(520)*V(88)

  B(880) = RCT(521)*V(228)

  B(881) = RCT(521)*V(105)

  B(882) = RCT(522)*V(228)

  B(883) = RCT(522)*V(117)

  B(884) = RCT(523)*V(228)

  B(885) = RCT(523)*V(82)

  B(886) = RCT(524)*V(228)

  B(887) = RCT(524)*V(137)

  B(888) = RCT(525)*V(228)

  B(889) = RCT(525)*V(129)

  B(890) = RCT(526)*V(228)

  B(891) = RCT(526)*V(60)

  B(892) = RCT(527)*V(228)

  B(893) = RCT(527)*V(83)

  B(894) = RCT(528)*V(228)

  B(895) = RCT(528)*V(61)

  B(896) = RCT(529)*V(228)

  B(897) = RCT(529)*V(89)

  B(898) = RCT(530)*V(228)

  B(899) = RCT(530)*V(91)

  B(900) = RCT(531)*V(228)

  B(901) = RCT(531)*V(130)

  B(902) = RCT(532)*V(228)

  B(903) = RCT(532)*V(67)

  B(904) = RCT(533)*V(228)

  B(905) = RCT(533)*V(68)

  B(906) = RCT(534)*V(228)

  B(907) = RCT(534)*V(69)

  B(908) = RCT(535)*V(228)

  B(909) = RCT(535)*V(70)

  B(910) = RCT(536)*V(228)

  B(911) = RCT(536)*V(71)

  B(912) = RCT(537)*V(228)

  B(913) = RCT(537)*V(93)

  B(914) = RCT(538)*V(228)

  B(915) = RCT(538)*V(123)

  B(916) = RCT(539)*V(228)

  B(917) = RCT(539)*V(94)

  B(918) = RCT(540)*V(228)

  B(919) = RCT(540)*V(124)

  B(920) = RCT(541)*V(228)

  B(921) = RCT(541)*V(72)

  B(922) = RCT(542)*V(228)

  B(923) = RCT(542)*V(100)

  B(924) = RCT(543)*V(228)

  B(925) = RCT(543)*V(125)

  B(926) = RCT(544)*V(229)

  B(927) = RCT(544)*V(219)

  B(928) = RCT(545)

  B(929) = RCT(546)*V(229)

  B(930) = RCT(546)*V(203)

  B(931) = RCT(547)

  B(932) = RCT(548)*V(229)

  B(933) = RCT(548)*V(191)

  B(934) = RCT(549)

  B(935) = RCT(550)*V(228)

  B(936) = RCT(550)*V(126)

  B(937) = RCT(551)*V(228)

  B(938) = RCT(551)*V(111)

  B(939) = RCT(552)*V(228)

  B(940) = RCT(552)*V(112)

  B(941) = RCT(553)*V(229)

  B(942) = RCT(553)*V(198)

  B(943) = RCT(554)

  B(944) = RCT(555)*V(229)

  B(945) = RCT(555)*V(197)

  B(946) = RCT(556)

  B(947) = RCT(557)*V(228)

  B(948) = RCT(557)*V(116)

  B(949) = RCT(558)*V(228)

  B(950) = RCT(558)*V(87)

  B(951) = RCT(559)*V(229)

  B(952) = RCT(559)*V(221)

  B(953) = RCT(560)

  B(954) = RCT(561)*V(228)

  B(955) = RCT(561)*V(122)

  B(956) = RCT(562)*V(228)

  B(957) = RCT(562)*V(76)

  B(958) = RCT(563)*V(228)

  B(959) = RCT(563)*V(154)

  B(960) = RCT(564)*V(228)

  B(961) = RCT(564)*V(182)

  B(962) = RCT(565)*V(228)

  B(963) = RCT(565)*V(167)

  B(964) = RCT(566)*V(228)

  B(965) = RCT(566)*V(11)

  B(966) = RCT(567)*V(228)

  B(967) = RCT(567)*V(12)

  B(968) = RCT(568)*V(228)

  B(969) = RCT(568)*V(13)

  B(970) = RCT(569)*V(228)

  B(971) = RCT(569)*V(14)

  B(972) = RCT(570)*V(228)

  B(973) = RCT(570)*V(18)

  B(974) = RCT(571)*V(228)

  B(975) = RCT(571)*V(17)

  B(976) = RCT(572)*V(228)

  B(977) = RCT(572)*V(15)

  B(978) = RCT(573)*V(228)

  B(979) = RCT(573)*V(16)

  B(980) = RCT(574)*V(229)

  B(981) = RCT(574)*V(211)

  B(982) = RCT(575)

  B(983) = RCT(576)

  B(984) = RCT(577)

  B(985) = RCT(578)

  B(986) = RCT(579)

  B(987) = RCT(580)

  B(988) = RCT(581)

  B(989) = RCT(582)

  B(990) = RCT(583)*V(227)

  B(991) = RCT(583)*V(204)

  B(992) = RCT(584)*V(228)

  B(993) = RCT(584)*V(160)

  B(994) = RCT(585)*V(228)

  B(995) = RCT(585)*V(77)

  B(996) = RCT(586)*V(228)

  B(997) = RCT(586)*V(135)

  B(998) = RCT(587)*V(230)

  B(999) = RCT(587)*V(135)

  B(1000) = RCT(588)*V(228)

  B(1001) = RCT(588)*V(27)

  B(1002) = RCT(589)*V(228)

  B(1003) = RCT(589)*V(27)

  B(1004) = RCT(590)*V(228)

  B(1005) = RCT(590)*V(28)

  B(1006) = RCT(591)*V(228)

  B(1007) = RCT(591)*V(28)

  B(1008) = RCT(592)*V(228)

  B(1009) = RCT(592)*V(20)

  B(1010) = RCT(593)*V(228)

  B(1011) = RCT(593)*V(21)

  B(1012) = RCT(594)*V(228)

  B(1013) = RCT(594)*V(22)

  B(1014) = RCT(595)*V(228)

  B(1015) = RCT(595)*V(23)

  B(1016) = RCT(596)*V(227)

  B(1017) = RCT(596)*V(155)

  B(1018) = RCT(597)*V(227)

  B(1019) = RCT(597)*V(200)

  B(1020) = RCT(598)*V(227)

  B(1021) = RCT(598)*V(155)

  B(1022) = RCT(599)*V(227)

  B(1023) = RCT(599)*V(200)

  B(1024) = RCT(600)*V(230)

  B(1025) = RCT(600)*V(155)

  B(1026) = RCT(601)*V(230)

  B(1027) = RCT(601)*V(200)

  B(1028) = RCT(602)*V(232)

  B(1029) = RCT(602)*V(155)

  B(1030) = RCT(603)*V(232)

  B(1031) = RCT(603)*V(200)

  B(1032) = RCT(604)

  B(1033) = RCT(605)

  B(1034) = RCT(606)*V(228)

  B(1035) = RCT(606)*V(141)

  B(1036) = RCT(607)*V(228)

  B(1037) = RCT(607)*V(133)

  B(1038) = RCT(608)*V(228)

  B(1039) = RCT(608)*V(24)

  B(1040) = RCT(609)*V(227)

  B(1041) = RCT(609)*V(174)

  B(1042) = RCT(610)*V(227)

  B(1043) = RCT(610)*V(174)

  B(1044) = RCT(611)*V(230)

  B(1045) = RCT(611)*V(174)

  B(1046) = RCT(612)*V(232)

  B(1047) = RCT(612)*V(174)

  B(1048) = RCT(613)

  B(1049) = RCT(614)*V(228)

  B(1050) = RCT(614)*V(128)

  B(1051) = RCT(615)*V(228)

  B(1052) = RCT(615)*V(104)

  B(1053) = RCT(616)*V(228)

  B(1054) = RCT(616)*V(152)

  B(1055) = RCT(617)*V(228)

  B(1056) = RCT(617)*V(152)

  B(1057) = RCT(618)

  B(1058) = RCT(619)

  B(1059) = RCT(620)*V(228)

  B(1060) = RCT(620)*V(108)

  B(1061) = RCT(621)*V(228)

  B(1062) = RCT(621)*V(108)

  B(1063) = RCT(622)*V(230)

  B(1064) = RCT(622)*V(108)

  B(1065) = RCT(623)*V(227)

  B(1066) = RCT(623)*V(110)

  B(1067) = RCT(624)*2*V(110)

  B(1068) = RCT(625)*V(222)

  B(1069) = RCT(625)*V(142)

  B(1070) = RCT(626)*V(229)

  B(1071) = RCT(626)*V(142)

  B(1072) = RCT(627)*V(229)

  B(1073) = RCT(627)*V(172)

  B(1074) = RCT(628)*V(222)

  B(1075) = RCT(628)*V(172)

  B(1076) = RCT(629)

  B(1077) = RCT(630)*V(229)

  B(1078) = RCT(630)*V(171)

  B(1079) = RCT(631)*V(222)

  B(1080) = RCT(631)*V(171)

  B(1081) = RCT(632)*V(232)

  B(1082) = RCT(632)*V(52)

  B(1083) = RCT(633)

  B(1084) = RCT(634)*V(228)

  B(1085) = RCT(634)*V(41)

  B(1086) = RCT(635)*V(228)

  B(1087) = RCT(635)*V(109)

  B(1088) = RCT(636)*V(230)

  B(1089) = RCT(636)*V(109)



  JVS(1) = 0

  JVS(2) = 0

  JVS(3) = 0

  JVS(4) = B(164)

  JVS(5) = B(1083)

  JVS(6) = -B(164)

  JVS(7) = B(162)

  JVS(8) = B(158)

  JVS(9) = B(159)

  JVS(10) = 0

  JVS(11) = B(1081)

  JVS(12) = 0.05*B(1086)

  JVS(13) = 0.05*B(1087)

  JVS(14) = B(1082)

  JVS(15) = 0

  JVS(16) = 0.05*B(1084)

  JVS(17) = 0.05*B(1085)

  JVS(18) = -B(107)-B(122)

  JVS(19) = B(1)

  JVS(20) = -B(168)

  JVS(21) = -B(169)

  JVS(22) = -B(162)

  JVS(23) = B(160)

  JVS(24) = B(161)

  JVS(25) = -B(174)

  JVS(26) = -B(175)

  JVS(27) = -B(964)

  JVS(28) = -B(965)

  JVS(29) = -B(966)

  JVS(30) = -B(967)

  JVS(31) = -B(968)

  JVS(32) = -B(969)

  JVS(33) = -B(970)

  JVS(34) = -B(971)

  JVS(35) = -B(976)

  JVS(36) = -B(977)

  JVS(37) = -B(978)

  JVS(38) = -B(979)

  JVS(39) = -B(974)

  JVS(40) = -B(975)

  JVS(41) = -B(972)

  JVS(42) = -B(973)

  JVS(43) = -B(166)

  JVS(44) = -B(167)

  JVS(45) = -B(1008)

  JVS(46) = -B(1009)

  JVS(47) = -B(1010)

  JVS(48) = -B(1011)

  JVS(49) = -B(1012)

  JVS(50) = -B(1013)

  JVS(51) = -B(1014)

  JVS(52) = -B(1015)

  JVS(53) = -B(1038)

  JVS(54) = -B(1039)

  JVS(55) = -B(170)-B(172)

  JVS(56) = -B(171)-B(173)

  JVS(57) = -B(242)-B(244)

  JVS(58) = -B(243)-B(245)

  JVS(59) = -B(1000)-B(1002)

  JVS(60) = -B(1001)-B(1003)

  JVS(61) = -B(1004)-B(1006)

  JVS(62) = -B(1005)-B(1007)

  JVS(63) = -B(234)-B(236)

  JVS(64) = -B(235)-B(237)

  JVS(65) = -B(238)-B(240)

  JVS(66) = -B(239)-B(241)

  JVS(67) = -B(14)-B(260)

  JVS(68) = -B(261)

  JVS(69) = -B(121)

  JVS(70) = B(119)

  JVS(71) = B(120)

  JVS(72) = -B(764)

  JVS(73) = B(812)

  JVS(74) = -B(765)+B(813)

  JVS(75) = -B(98)

  JVS(76) = B(980)

  JVS(77) = B(981)

  JVS(78) = -B(230)-B(232)

  JVS(79) = -B(231)-B(233)

  JVS(80) = -B(262)

  JVS(81) = B(646)

  JVS(82) = -B(263)

  JVS(83) = -B(28)-B(740)

  JVS(84) = B(868)

  JVS(85) = -B(741)+B(869)

  JVS(86) = -B(742)

  JVS(87) = B(818)

  JVS(88) = B(922)

  JVS(89) = -B(743)+B(819)+B(923)

  JVS(90) = -B(772)

  JVS(91) = B(216)

  JVS(92) = B(217)

  JVS(93) = -B(773)

  JVS(94) = -B(7)-B(154)

  JVS(95) = B(139)

  JVS(96) = B(140)-B(155)

  JVS(97) = -B(1084)

  JVS(98) = B(1061)

  JVS(99) = B(1062)-B(1085)

  JVS(100) = -B(264)-B(266)

  JVS(101) = B(649)

  JVS(102) = -B(265)-B(267)

  JVS(103) = -B(272)-B(274)

  JVS(104) = B(655)

  JVS(105) = -B(273)-B(275)

  JVS(106) = -B(27)-B(738)

  JVS(107) = B(866)

  JVS(108) = -B(739)+B(867)

  JVS(109) = B(230)

  JVS(110) = -B(276)

  JVS(111) = B(190)

  JVS(112) = B(206)

  JVS(113) = B(191)+B(207)

  JVS(114) = B(231)-B(277)

  JVS(115) = B(236)

  JVS(116) = -B(820)-B(822)

  JVS(117) = B(237)-B(821)

  JVS(118) = -B(823)

  JVS(119) = B(820)+B(822)

  JVS(120) = -B(824)

  JVS(121) = B(821)

  JVS(122) = -B(825)

  JVS(123) = B(823)

  JVS(124) = B(240)

  JVS(125) = -B(830)-B(832)

  JVS(126) = B(241)-B(831)

  JVS(127) = -B(833)

  JVS(128) = B(830)+B(832)

  JVS(129) = -B(834)

  JVS(130) = B(831)

  JVS(131) = -B(835)

  JVS(132) = B(833)

  JVS(133) = -B(268)-B(270)

  JVS(134) = B(652)

  JVS(135) = -B(269)-B(271)

  JVS(136) = -B(278)

  JVS(137) = B(194)

  JVS(138) = B(198)

  JVS(139) = B(195)+B(199)

  JVS(140) = -B(279)

  JVS(141) = -B(1081)-B(1083)

  JVS(142) = B(1077)+B(1079)

  JVS(143) = B(1080)

  JVS(144) = B(1078)

  JVS(145) = -B(1082)

  JVS(146) = -B(57)-B(856)

  JVS(147) = B(590)

  JVS(148) = -B(857)

  JVS(149) = B(591)

  JVS(150) = -B(794)

  JVS(151) = B(410)

  JVS(152) = B(411)

  JVS(153) = -B(795)

  JVS(154) = -B(796)

  JVS(155) = B(412)

  JVS(156) = B(413)

  JVS(157) = -B(797)

  JVS(158) = -B(51)-B(846)

  JVS(159) = B(554)

  JVS(160) = -B(847)

  JVS(161) = B(555)

  JVS(162) = -B(56)-B(854)

  JVS(163) = B(562)

  JVS(164) = -B(855)

  JVS(165) = B(563)

  JVS(166) = -B(60)-B(862)

  JVS(167) = B(596)

  JVS(168) = -B(863)

  JVS(169) = B(597)

  JVS(170) = -B(61)-B(864)

  JVS(171) = B(598)

  JVS(172) = -B(865)

  JVS(173) = B(599)

  JVS(174) = -B(75)-B(890)

  JVS(175) = B(610)

  JVS(176) = -B(891)

  JVS(177) = B(611)

  JVS(178) = -B(77)-B(894)

  JVS(179) = B(614)

  JVS(180) = -B(895)

  JVS(181) = B(615)

  JVS(182) = -B(58)-B(858)

  JVS(183) = B(592)

  JVS(184) = -B(859)

  JVS(185) = B(593)

  JVS(186) = -B(802)

  JVS(187) = B(422)

  JVS(188) = B(423)

  JVS(189) = -B(803)

  JVS(190) = -B(68)-B(876)

  JVS(191) = B(616)

  JVS(192) = -B(877)

  JVS(193) = B(617)

  JVS(194) = -B(792)

  JVS(195) = B(408)

  JVS(196) = B(409)

  JVS(197) = -B(793)

  JVS(198) = -B(810)

  JVS(199) = B(432)

  JVS(200) = B(433)

  JVS(201) = -B(811)

  JVS(202) = -B(81)-B(902)

  JVS(203) = B(620)

  JVS(204) = -B(903)

  JVS(205) = B(621)

  JVS(206) = -B(83)-B(904)

  JVS(207) = B(624)

  JVS(208) = -B(905)

  JVS(209) = B(625)

  JVS(210) = -B(82)-B(906)

  JVS(211) = B(622)

  JVS(212) = -B(907)

  JVS(213) = B(623)

  JVS(214) = -B(84)-B(908)

  JVS(215) = B(626)

  JVS(216) = -B(909)

  JVS(217) = B(627)

  JVS(218) = -B(85)-B(910)

  JVS(219) = B(628)

  JVS(220) = -B(911)

  JVS(221) = B(629)

  JVS(222) = -B(90)-B(920)

  JVS(223) = B(638)

  JVS(224) = -B(921)

  JVS(225) = B(639)

  JVS(226) = -B(50)-B(844)

  JVS(227) = B(552)

  JVS(228) = -B(845)

  JVS(229) = B(553)

  JVS(230) = -B(149)-B(150)

  JVS(231) = -B(151)

  JVS(232) = B(147)

  JVS(233) = B(148)

  JVS(234) = -B(52)-B(848)

  JVS(235) = B(556)

  JVS(236) = -B(849)

  JVS(237) = B(557)

  JVS(238) = -B(956)

  JVS(239) = B(988)+B(1055)

  JVS(240) = B(93)+B(958)

  JVS(241) = B(95)+B(962)

  JVS(242) = B(94)+B(960)

  JVS(243) = -B(957)+B(959)+B(961)+B(963)+B(1056)

  JVS(244) = -B(982)-B(994)

  JVS(245) = B(992)

  JVS(246) = B(993)-B(995)

  JVS(247) = -B(800)

  JVS(248) = B(416)

  JVS(249) = B(420)

  JVS(250) = B(417)+B(421)

  JVS(251) = -B(801)

  JVS(252) = -B(43)-B(788)

  JVS(253) = B(404)

  JVS(254) = B(405)

  JVS(255) = -B(789)

  JVS(256) = -B(55)-B(852)

  JVS(257) = B(566)

  JVS(258) = B(560)

  JVS(259) = -B(853)

  JVS(260) = B(561)+B(567)

  JVS(261) = -B(63)-B(868)

  JVS(262) = B(602)

  JVS(263) = -B(869)

  JVS(264) = B(603)

  JVS(265) = -B(72)-B(884)

  JVS(266) = B(580)

  JVS(267) = -B(885)

  JVS(268) = B(581)

  JVS(269) = -B(76)-B(892)

  JVS(270) = B(612)

  JVS(271) = -B(893)

  JVS(272) = B(613)

  JVS(273) = -B(59)-B(860)

  JVS(274) = B(594)

  JVS(275) = -B(861)

  JVS(276) = B(595)

  JVS(277) = -B(40)-B(784)

  JVS(278) = B(400)

  JVS(279) = B(401)

  JVS(280) = -B(785)

  JVS(281) = -B(67)-B(874)

  JVS(282) = B(608)

  JVS(283) = -B(875)

  JVS(284) = B(609)

  JVS(285) = -B(943)-B(949)

  JVS(286) = B(941)

  JVS(287) = -B(950)

  JVS(288) = B(942)

  JVS(289) = -B(69)-B(878)

  JVS(290) = B(618)

  JVS(291) = -B(879)

  JVS(292) = B(619)

  JVS(293) = -B(78)-B(896)

  JVS(294) = B(568)

  JVS(295) = -B(897)

  JVS(296) = B(569)

  JVS(297) = -B(37)-B(778)

  JVS(298) = B(394)

  JVS(299) = B(395)

  JVS(300) = -B(779)

  JVS(301) = -B(79)-B(898)

  JVS(302) = B(570)

  JVS(303) = -B(899)

  JVS(304) = B(571)

  JVS(305) = -B(38)-B(780)

  JVS(306) = B(396)

  JVS(307) = B(397)

  JVS(308) = -B(781)

  JVS(309) = -B(86)-B(912)

  JVS(310) = B(630)

  JVS(311) = -B(913)

  JVS(312) = B(631)

  JVS(313) = -B(88)-B(916)

  JVS(314) = B(634)

  JVS(315) = -B(917)

  JVS(316) = B(635)

  JVS(317) = -B(812)

  JVS(318) = B(434)

  JVS(319) = B(435)

  JVS(320) = -B(813)

  JVS(321) = -B(814)

  JVS(322) = B(436)

  JVS(323) = B(437)

  JVS(324) = -B(815)

  JVS(325) = -B(48)-B(816)

  JVS(326) = B(438)

  JVS(327) = B(439)

  JVS(328) = -B(817)

  JVS(329) = -B(39)-B(782)

  JVS(330) = B(398)

  JVS(331) = B(399)

  JVS(332) = -B(783)

  JVS(333) = -B(818)

  JVS(334) = B(440)

  JVS(335) = B(441)

  JVS(336) = -B(819)

  JVS(337) = -B(91)-B(922)

  JVS(338) = B(640)

  JVS(339) = -B(923)

  JVS(340) = B(641)

  JVS(341) = -B(49)-B(840)-B(842)

  JVS(342) = B(550)

  JVS(343) = -B(841)-B(843)

  JVS(344) = B(551)

  JVS(345) = -B(798)

  JVS(346) = B(414)

  JVS(347) = B(418)

  JVS(348) = B(415)+B(419)

  JVS(349) = -B(799)

  JVS(350) = -B(44)-B(790)

  JVS(351) = B(406)

  JVS(352) = B(407)

  JVS(353) = -B(791)

  JVS(354) = -B(987)-B(1051)

  JVS(355) = B(1046)

  JVS(356) = -B(1052)

  JVS(357) = B(1047)

  JVS(358) = -B(70)-B(880)

  JVS(359) = B(576)

  JVS(360) = -B(881)

  JVS(361) = B(577)

  JVS(362) = -B(3)-B(130)

  JVS(363) = B(224)

  JVS(364) = B(214)

  JVS(365) = B(750)

  JVS(366) = B(726)

  JVS(367) = B(215)+B(225)+B(727)+B(751)

  JVS(368) = -B(131)

  JVS(369) = B(136)+B(137)

  JVS(370) = B(90)+B(920)

  JVS(371) = B(48)+B(816)

  JVS(372) = -B(776)

  JVS(373) = B(226)

  JVS(374) = B(388)+B(544)+B(703)

  JVS(375) = B(227)

  JVS(376) = B(389)

  JVS(377) = -B(777)+B(817)+B(921)

  JVS(378) = B(545)

  JVS(379) = 0

  JVS(380) = -B(1059)-B(1061)-B(1063)

  JVS(381) = -B(1060)-B(1062)

  JVS(382) = -B(1064)

  JVS(383) = 0.95*B(1084)

  JVS(384) = 0

  JVS(385) = -B(1086)-B(1088)

  JVS(386) = 0.95*B(1085)-B(1087)

  JVS(387) = -B(1089)

  JVS(388) = B(1059)+B(1063)

  JVS(389) = -B(1065)-2*B(1067)

  JVS(390) = -B(1066)

  JVS(391) = B(1060)

  JVS(392) = B(1064)

  JVS(393) = -B(931)-B(937)

  JVS(394) = B(929)

  JVS(395) = -B(938)

  JVS(396) = B(930)

  JVS(397) = -B(934)-B(939)

  JVS(398) = B(932)

  JVS(399) = -B(940)

  JVS(400) = B(933)

  JVS(401) = -B(64)-B(65)-B(870)

  JVS(402) = B(604)

  JVS(403) = -B(871)

  JVS(404) = B(605)

  JVS(405) = -B(66)-B(872)

  JVS(406) = B(606)

  JVS(407) = -B(873)

  JVS(408) = B(607)

  JVS(409) = -B(62)-B(866)

  JVS(410) = B(600)

  JVS(411) = -B(867)

  JVS(412) = B(601)

  JVS(413) = -B(946)-B(947)

  JVS(414) = B(944)

  JVS(415) = -B(948)

  JVS(416) = B(945)

  JVS(417) = -B(71)-B(882)

  JVS(418) = B(578)

  JVS(419) = -B(883)

  JVS(420) = B(579)

  JVS(421) = -B(45)-B(804)

  JVS(422) = B(424)

  JVS(423) = B(425)

  JVS(424) = -B(805)

  JVS(425) = B(824)

  JVS(426) = -B(826)-B(828)

  JVS(427) = -B(827)

  JVS(428) = B(825)

  JVS(429) = -B(829)

  JVS(430) = -B(46)-B(806)

  JVS(431) = B(426)

  JVS(432) = B(427)

  JVS(433) = -B(807)

  JVS(434) = B(834)

  JVS(435) = -B(836)-B(838)

  JVS(436) = -B(837)

  JVS(437) = B(835)

  JVS(438) = -B(839)

  JVS(439) = -B(953)-B(954)

  JVS(440) = B(951)

  JVS(441) = -B(955)

  JVS(442) = B(952)

  JVS(443) = -B(87)-B(914)

  JVS(444) = B(632)

  JVS(445) = -B(915)

  JVS(446) = B(633)

  JVS(447) = -B(89)-B(918)

  JVS(448) = B(636)

  JVS(449) = -B(919)

  JVS(450) = B(637)

  JVS(451) = -B(92)-B(924)

  JVS(452) = B(642)

  JVS(453) = -B(925)

  JVS(454) = B(643)

  JVS(455) = -B(928)-B(935)

  JVS(456) = B(926)

  JVS(457) = -B(936)

  JVS(458) = B(927)

  JVS(459) = -B(158)-B(160)

  JVS(460) = B(1076)

  JVS(461) = 0.18*B(1072)

  JVS(462) = -B(159)

  JVS(463) = -B(161)

  JVS(464) = 0.18*B(1073)

  JVS(465) = -B(986)-B(1049)

  JVS(466) = B(1042)

  JVS(467) = B(1043)

  JVS(468) = -B(1050)

  JVS(469) = -B(74)-B(888)

  JVS(470) = B(584)

  JVS(471) = B(588)

  JVS(472) = -B(889)

  JVS(473) = B(585)+B(589)

  JVS(474) = -B(80)-B(900)

  JVS(475) = B(572)

  JVS(476) = B(574)

  JVS(477) = -B(901)

  JVS(478) = B(573)+B(575)

  JVS(479) = B(854)

  JVS(480) = B(790)

  JVS(481) = -B(16)-B(708)

  JVS(482) = 0

  JVS(483) = 0

  JVS(484) = -B(709)+B(791)+B(855)

  JVS(485) = 0

  JVS(486) = B(814)

  JVS(487) = B(89)

  JVS(488) = B(92)

  JVS(489) = -B(774)

  JVS(490) = B(224)+B(228)

  JVS(491) = B(386)+B(548)+B(705)

  JVS(492) = B(382)+B(542)+B(702)

  JVS(493) = B(225)+B(229)

  JVS(494) = B(383)+B(387)

  JVS(495) = -B(775)+B(815)

  JVS(496) = B(543)+B(549)

  JVS(497) = 0

  JVS(498) = -B(985)-B(1036)

  JVS(499) = B(1028)

  JVS(500) = B(1030)

  JVS(501) = -B(1037)

  JVS(502) = B(1029)+B(1031)

  JVS(503) = -B(41)-B(42)-B(786)

  JVS(504) = B(402)

  JVS(505) = B(403)

  JVS(506) = -B(787)

  JVS(507) = B(994)

  JVS(508) = -B(983)-B(996)-B(998)

  JVS(509) = 0

  JVS(510) = B(995)-B(997)

  JVS(511) = -B(999)

  JVS(512) = -B(53)-B(54)-B(850)

  JVS(513) = B(564)

  JVS(514) = B(558)

  JVS(515) = -B(851)

  JVS(516) = B(559)+B(565)

  JVS(517) = -B(73)-B(886)

  JVS(518) = B(582)

  JVS(519) = B(586)

  JVS(520) = -B(887)

  JVS(521) = B(583)+B(587)

  JVS(522) = -B(47)-B(808)

  JVS(523) = B(428)

  JVS(524) = B(430)

  JVS(525) = B(429)+B(431)

  JVS(526) = -B(809)

  JVS(527) = B(788)

  JVS(528) = B(852)

  JVS(529) = -B(15)-B(706)

  JVS(530) = 0

  JVS(531) = 0

  JVS(532) = 0

  JVS(533) = -B(707)+B(789)+B(853)

  JVS(534) = 0

  JVS(535) = B(796)

  JVS(536) = B(884)

  JVS(537) = -B(19)-B(714)

  JVS(538) = 0

  JVS(539) = 0

  JVS(540) = -B(715)+B(797)+B(885)

  JVS(541) = 0

  JVS(542) = -B(984)-B(1034)

  JVS(543) = B(1020)

  JVS(544) = B(1022)

  JVS(545) = B(1021)+B(1023)

  JVS(546) = -B(1035)

  JVS(547) = B(1065)+2*B(1067)

  JVS(548) = -B(1068)-B(1070)

  JVS(549) = -B(1069)

  JVS(550) = B(1066)

  JVS(551) = 0

  JVS(552) = -B(1071)

  JVS(553) = 0

  JVS(554) = -B(178)-B(184)-B(192)-B(194)

  JVS(555) = -B(193)-B(195)

  JVS(556) = -B(179)

  JVS(557) = -B(185)

  JVS(558) = -B(360)-B(522)-B(616)-B(692)

  JVS(559) = B(202)

  JVS(560) = -B(361)

  JVS(561) = B(203)-B(523)

  JVS(562) = -B(617)

  JVS(563) = -B(358)-B(520)-B(614)-B(691)

  JVS(564) = B(186)

  JVS(565) = -B(359)

  JVS(566) = B(187)-B(521)

  JVS(567) = -B(615)

  JVS(568) = -B(368)-B(528)-B(622)-B(695)

  JVS(569) = B(210)

  JVS(570) = -B(369)

  JVS(571) = B(211)-B(529)

  JVS(572) = -B(623)

  JVS(573) = -B(354)-B(516)-B(610)-B(689)

  JVS(574) = B(182)

  JVS(575) = -B(355)

  JVS(576) = B(183)-B(517)

  JVS(577) = -B(611)

  JVS(578) = B(864)

  JVS(579) = B(714)

  JVS(580) = -B(26)-B(736)

  JVS(581) = B(1018)+B(1026)+B(1033)

  JVS(582) = B(308)+B(470)

  JVS(583) = 0

  JVS(584) = 0

  JVS(585) = B(309)+B(1019)

  JVS(586) = B(715)-B(737)+B(865)

  JVS(587) = B(471)+B(1027)

  JVS(588) = 0

  JVS(589) = B(42)+B(786)

  JVS(590) = B(54)+B(850)

  JVS(591) = -B(17)-B(710)

  JVS(592) = 0

  JVS(593) = B(290)+B(452)+B(657)

  JVS(594) = B(291)

  JVS(595) = -B(711)+B(787)+B(851)

  JVS(596) = B(453)

  JVS(597) = 0

  JVS(598) = B(232)

  JVS(599) = B(78)+B(896)

  JVS(600) = B(79)+B(898)

  JVS(601) = B(45)+B(804)

  JVS(602) = B(46)+B(806)

  JVS(603) = -B(24)-B(732)

  JVS(604) = B(300)+B(462)+B(660)

  JVS(605) = B(302)+B(464)+B(661)

  JVS(606) = B(306)+B(468)+B(663)

  JVS(607) = B(664)

  JVS(608) = B(301)+B(303)+B(307)

  JVS(609) = B(233)-B(733)+B(805)+B(807)+B(897)+B(899)

  JVS(610) = B(463)+B(465)+B(469)

  JVS(611) = 0

  JVS(612) = B(68)+B(876)

  JVS(613) = B(360)+B(522)+B(692)

  JVS(614) = -B(30)-B(752)

  JVS(615) = 0

  JVS(616) = B(361)

  JVS(617) = -B(753)+B(877)

  JVS(618) = B(523)

  JVS(619) = 0

  JVS(620) = B(987)+B(1051)

  JVS(621) = B(986)+B(1049)

  JVS(622) = -B(988)-B(989)-B(1053)-B(1055)

  JVS(623) = B(1040)+B(1044)+B(1048)

  JVS(624) = B(1041)

  JVS(625) = B(1050)+B(1052)-B(1054)-B(1056)

  JVS(626) = B(1045)

  JVS(627) = 0

  JVS(628) = -B(176)-B(182)-B(188)-B(190)

  JVS(629) = -B(189)-B(191)

  JVS(630) = -B(177)

  JVS(631) = -B(183)

  JVS(632) = B(78)+B(896)

  JVS(633) = B(45)+B(804)

  JVS(634) = -B(32)-B(93)-B(756)-B(958)

  JVS(635) = B(300)+B(462)+B(660)

  JVS(636) = B(304)+B(466)+B(662)

  JVS(637) = B(308)+B(470)

  JVS(638) = B(301)+B(305)+B(309)

  JVS(639) = -B(757)+B(805)+B(897)-B(959)

  JVS(640) = B(463)+B(467)+B(471)

  JVS(641) = 0

  JVS(642) = B(1008)

  JVS(643) = B(1010)

  JVS(644) = B(1012)

  JVS(645) = B(1014)

  JVS(646) = B(1000)

  JVS(647) = B(1004)

  JVS(648) = -B(1016)-B(1020)-B(1024)-B(1028)-B(1032)

  JVS(649) = -B(1017)-B(1021)

  JVS(650) = B(1001)+B(1005)+B(1009)+B(1011)+B(1013)+B(1015)

  JVS(651) = -B(1025)

  JVS(652) = -B(1029)

  JVS(653) = B(170)

  JVS(654) = -B(286)-B(400)-B(448)-B(556)-B(653)-B(654)-B(655)

  JVS(655) = -B(287)-B(401)

  JVS(656) = B(171)

  JVS(657) = -B(449)

  JVS(658) = -B(557)

  JVS(659) = B(84)

  JVS(660) = B(910)

  JVS(661) = -B(374)-B(534)-B(628)-B(698)

  JVS(662) = B(372)+B(532)+B(697)

  JVS(663) = B(373)-B(375)

  JVS(664) = B(911)

  JVS(665) = B(533)-B(535)

  JVS(666) = -B(629)

  JVS(667) = B(87)

  JVS(668) = -B(36)-B(766)-B(770)

  JVS(669) = B(378)+B(538)+B(700)

  JVS(670) = B(379)

  JVS(671) = -B(767)

  JVS(672) = B(539)-B(771)

  JVS(673) = 0

  JVS(674) = -B(218)-B(220)-B(222)-B(224)-B(226)-B(228)

  JVS(675) = -B(223)-B(225)-B(227)-B(229)

  JVS(676) = -B(219)

  JVS(677) = -B(221)

  JVS(678) = -B(992)

  JVS(679) = B(990)

  JVS(680) = B(991)

  JVS(681) = -B(993)

  JVS(682) = B(69)+B(878)

  JVS(683) = B(30)

  JVS(684) = -B(31)-B(754)

  JVS(685) = 0

  JVS(686) = B(362)+B(524)+B(693)

  JVS(687) = B(363)

  JVS(688) = -B(755)+B(879)

  JVS(689) = B(525)

  JVS(690) = 0

  JVS(691) = -B(180)-B(186)-B(196)-B(198)

  JVS(692) = -B(197)-B(199)

  JVS(693) = -B(181)

  JVS(694) = -B(187)

  JVS(695) = B(234)

  JVS(696) = -B(300)-B(424)-B(462)-B(568)-B(660)

  JVS(697) = -B(301)-B(425)

  JVS(698) = B(235)

  JVS(699) = -B(463)

  JVS(700) = -B(569)

  JVS(701) = -B(208)-B(210)-B(212)-B(214)-B(216)

  JVS(702) = -B(213)-B(215)-B(217)

  JVS(703) = -B(209)

  JVS(704) = -B(211)

  JVS(705) = B(810)

  JVS(706) = B(81)+B(902)

  JVS(707) = B(82)+B(906)

  JVS(708) = B(368)+B(528)+B(695)

  JVS(709) = B(214)

  JVS(710) = -B(35)-B(762)-B(768)

  JVS(711) = B(364)+B(526)+B(694)

  JVS(712) = B(215)

  JVS(713) = B(365)+B(369)

  JVS(714) = -B(763)+B(811)+B(903)+B(907)

  JVS(715) = B(527)+B(529)-B(769)

  JVS(716) = 0

  JVS(717) = B(184)

  JVS(718) = -B(356)-B(518)-B(612)-B(690)

  JVS(719) = 0

  JVS(720) = -B(357)

  JVS(721) = 0

  JVS(722) = B(185)-B(519)

  JVS(723) = -B(613)

  JVS(724) = B(985)+B(1036)

  JVS(725) = B(984)+B(1034)

  JVS(726) = B(1016)+B(1024)+B(1032)

  JVS(727) = -B(34)-B(95)-B(760)-B(962)

  JVS(728) = B(306)+B(468)+B(663)

  JVS(729) = 0

  JVS(730) = B(664)

  JVS(731) = B(307)+B(1017)

  JVS(732) = -B(761)-B(963)+B(1035)+B(1037)

  JVS(733) = B(469)+B(1025)

  JVS(734) = 0

  JVS(735) = B(924)

  JVS(736) = B(220)

  JVS(737) = -B(386)-B(548)-B(642)-B(705)

  JVS(738) = 0

  JVS(739) = -B(387)

  JVS(740) = B(925)

  JVS(741) = B(221)-B(549)

  JVS(742) = -B(643)

  JVS(743) = B(798)

  JVS(744) = B(826)+B(828)

  JVS(745) = B(886)

  JVS(746) = -B(20)-B(716)

  JVS(747) = 0

  JVS(748) = B(322)+B(484)+B(673)

  JVS(749) = B(323)

  JVS(750) = -B(717)+B(799)+B(827)+B(887)

  JVS(751) = 0

  JVS(752) = B(485)+B(829)

  JVS(753) = 0

  JVS(754) = B(44)

  JVS(755) = -B(320)-B(416)-B(482)-B(584)-B(672)

  JVS(756) = 0

  JVS(757) = -B(321)-B(417)

  JVS(758) = 0

  JVS(759) = -B(483)

  JVS(760) = -B(585)

  JVS(761) = 0.95*B(1086)+B(1088)

  JVS(762) = -B(1076)-B(1077)-B(1079)

  JVS(763) = 0.82*B(1072)+B(1074)

  JVS(764) = B(1075)-B(1080)

  JVS(765) = 0.95*B(1087)

  JVS(766) = 0.82*B(1073)-B(1078)

  JVS(767) = B(1089)

  JVS(768) = B(1068)+B(1070)

  JVS(769) = -B(1072)-B(1074)

  JVS(770) = B(1069)-B(1075)

  JVS(771) = 0

  JVS(772) = 0

  JVS(773) = B(1071)-B(1073)

  JVS(774) = 0

  JVS(775) = B(107)

  JVS(776) = -B(158)

  JVS(777) = 0

  JVS(778) = 0

  JVS(779) = -B(96)-B(99)-B(101)-B(103)-B(105)-B(159)

  JVS(780) = B(2)-B(100)

  JVS(781) = -B(102)

  JVS(782) = 0

  JVS(783) = B(4)-B(104)-B(106)

  JVS(784) = B(6)

  JVS(785) = B(1038)

  JVS(786) = -B(1040)-B(1042)-B(1044)-B(1046)-B(1048)

  JVS(787) = -B(1041)-B(1043)

  JVS(788) = B(1039)

  JVS(789) = -B(1045)

  JVS(790) = -B(1047)

  JVS(791) = -B(200)-B(202)-B(204)-B(206)

  JVS(792) = -B(205)-B(207)

  JVS(793) = -B(201)

  JVS(794) = -B(203)

  JVS(795) = -B(296)-B(458)-B(564)-B(658)

  JVS(796) = B(390)

  JVS(797) = -B(297)+B(391)

  JVS(798) = -B(459)

  JVS(799) = -B(565)

  JVS(800) = -B(298)-B(460)-B(566)-B(659)

  JVS(801) = B(338)+B(500)+B(681)

  JVS(802) = -B(299)+B(339)

  JVS(803) = -B(461)+B(501)

  JVS(804) = -B(567)

  JVS(805) = B(200)

  JVS(806) = -B(340)-B(342)-B(422)-B(502)-B(504)-B(604)-B(682)-B(683)

  JVS(807) = 0

  JVS(808) = -B(341)-B(343)-B(423)

  JVS(809) = B(201)

  JVS(810) = -B(503)-B(505)

  JVS(811) = -B(605)

  JVS(812) = B(794)

  JVS(813) = B(947)

  JVS(814) = B(882)

  JVS(815) = -B(18)-B(712)

  JVS(816) = B(352)+B(514)+B(688)

  JVS(817) = B(346)+B(508)+B(685)

  JVS(818) = 0

  JVS(819) = B(347)+B(353)

  JVS(820) = -B(713)+B(795)+B(883)+B(948)

  JVS(821) = 0

  JVS(822) = B(509)+B(515)

  JVS(823) = 0

  JVS(824) = B(800)

  JVS(825) = B(836)+B(838)

  JVS(826) = B(954)

  JVS(827) = B(888)

  JVS(828) = 0

  JVS(829) = -B(21)-B(718)

  JVS(830) = B(324)+B(486)+B(674)

  JVS(831) = 0

  JVS(832) = 0

  JVS(833) = B(325)

  JVS(834) = -B(719)+B(801)+B(837)+B(889)+B(955)

  JVS(835) = 0

  JVS(836) = B(487)+B(839)

  JVS(837) = 0

  JVS(838) = B(238)

  JVS(839) = -B(302)-B(304)-B(426)-B(464)-B(466)-B(570)-B(661)-B(662)

  JVS(840) = -B(303)-B(305)-B(427)

  JVS(841) = B(239)

  JVS(842) = -B(465)-B(467)

  JVS(843) = -B(571)

  JVS(844) = B(79)+B(898)

  JVS(845) = B(46)+B(806)

  JVS(846) = B(80)+B(900)

  JVS(847) = B(47)+B(808)

  JVS(848) = B(302)+B(464)+B(661)

  JVS(849) = -B(33)-B(94)-B(758)-B(960)

  JVS(850) = 0

  JVS(851) = B(1018)+B(1026)+B(1033)

  JVS(852) = 0

  JVS(853) = B(303)+B(1019)

  JVS(854) = -B(759)+B(807)+B(809)+B(899)+B(901)-B(961)

  JVS(855) = B(465)+B(1027)

  JVS(856) = 0

  JVS(857) = B(916)

  JVS(858) = B(983)+B(996)+B(998)

  JVS(859) = B(766)+B(770)

  JVS(860) = 0

  JVS(861) = -B(380)-B(540)-B(634)-B(701)

  JVS(862) = 0

  JVS(863) = 0

  JVS(864) = -B(381)

  JVS(865) = B(767)+B(917)+B(997)

  JVS(866) = -B(541)+B(771)+B(999)

  JVS(867) = -B(635)

  JVS(868) = B(918)

  JVS(869) = B(218)

  JVS(870) = -B(382)-B(384)-B(436)-B(542)-B(636)-B(702)

  JVS(871) = 0

  JVS(872) = -B(383)-B(385)-B(437)

  JVS(873) = B(219)+B(919)

  JVS(874) = -B(543)

  JVS(875) = -B(637)

  JVS(876) = B(822)

  JVS(877) = B(832)

  JVS(878) = B(1063)

  JVS(879) = B(1088)

  JVS(880) = B(828)

  JVS(881) = B(838)

  JVS(882) = B(998)

  JVS(883) = B(770)

  JVS(884) = 0

  JVS(885) = B(768)

  JVS(886) = -B(8)-B(156)

  JVS(887) = 0

  JVS(888) = 0

  JVS(889) = B(746)

  JVS(890) = B(730)

  JVS(891) = B(722)

  JVS(892) = 0

  JVS(893) = B(252)

  JVS(894) = 0

  JVS(895) = B(256)

  JVS(896) = 0

  JVS(897) = B(141)-B(157)

  JVS(898) = B(142)

  JVS(899) = B(253)+B(254)+B(257)+B(723)+B(731)+B(747)+B(769)+B(771)+B(823)+B(829)+B(833)+B(839)+B(999)+B(1064)+B(1089)

  JVS(900) = B(255)

  JVS(901) = 0

  JVS(902) = B(208)

  JVS(903) = -B(364)-B(366)-B(432)-B(526)-B(620)-B(694)

  JVS(904) = 0

  JVS(905) = -B(365)-B(367)-B(433)

  JVS(906) = B(209)

  JVS(907) = -B(527)

  JVS(908) = -B(621)

  JVS(909) = B(272)

  JVS(910) = B(52)+B(848)

  JVS(911) = B(40)+B(784)

  JVS(912) = B(86)

  JVS(913) = B(812)

  JVS(914) = B(818)

  JVS(915) = B(91)+B(922)

  JVS(916) = B(954)

  JVS(917) = B(286)+B(448)+B(653)+B(654)

  JVS(918) = B(992)

  JVS(919) = B(212)

  JVS(920) = B(384)

  JVS(921) = B(366)

  JVS(922) = -B(13)-B(258)

  JVS(923) = B(376)+B(536)+B(699)

  JVS(924) = B(392)+B(546)+B(704)

  JVS(925) = 0

  JVS(926) = B(390)

  JVS(927) = 0

  JVS(928) = B(213)

  JVS(929) = B(287)+B(367)+B(377)+B(385)+B(391)+B(393)

  JVS(930) = -B(259)+B(273)+B(785)+B(813)+B(819)+B(849)+B(923)+B(955)+B(993)

  JVS(931) = 0

  JVS(932) = B(449)+B(537)+B(547)

  JVS(933) = 0

  JVS(934) = B(752)

  JVS(935) = 0

  JVS(936) = -B(362)-B(524)-B(618)-B(693)

  JVS(937) = 0

  JVS(938) = -B(363)

  JVS(939) = B(753)

  JVS(940) = -B(525)

  JVS(941) = -B(619)

  JVS(942) = B(230)

  JVS(943) = B(69)+B(878)

  JVS(944) = B(88)

  JVS(945) = B(937)

  JVS(946) = B(939)

  JVS(947) = B(947)

  JVS(948) = B(935)

  JVS(949) = B(983)+B(996)+B(998)

  JVS(950) = B(192)

  JVS(951) = 2*B(24)+2*B(732)

  JVS(952) = 2*B(30)

  JVS(953) = B(188)

  JVS(954) = B(36)

  JVS(955) = B(228)

  JVS(956) = 0

  JVS(957) = B(196)

  JVS(958) = 0

  JVS(959) = B(204)

  JVS(960) = 0

  JVS(961) = B(380)+B(540)+B(701)

  JVS(962) = B(362)+B(524)+B(693)

  JVS(963) = -B(128)

  JVS(964) = 0

  JVS(965) = 0

  JVS(966) = 0

  JVS(967) = B(25)+B(734)

  JVS(968) = B(29)+B(748)

  JVS(969) = B(22)

  JVS(970) = 0

  JVS(971) = B(346)+B(508)

  JVS(972) = B(724)

  JVS(973) = 0

  JVS(974) = 0

  JVS(975) = 0

  JVS(976) = B(9)+B(10)+B(246)+B(252)

  JVS(977) = 0

  JVS(978) = B(189)+B(193)+B(197)+B(205)+B(229)+B(725)+B(749)

  JVS(979) = B(12)

  JVS(980) = B(347)+B(363)+B(381)

  JVS(981) = -B(129)+B(231)+B(247)+2*B(733)+B(735)+B(879)+B(936)+B(938)+B(940)+B(948)+B(997)

  JVS(982) = 0

  JVS(983) = B(253)+B(509)+B(525)+B(541)+B(999)

  JVS(984) = B(11)

  JVS(985) = 0

  JVS(986) = B(86)

  JVS(987) = B(914)

  JVS(988) = -B(378)-B(538)-B(632)-B(700)

  JVS(989) = B(376)+B(536)+B(699)

  JVS(990) = B(377)-B(379)

  JVS(991) = B(915)

  JVS(992) = B(537)-B(539)

  JVS(993) = -B(633)

  JVS(994) = B(982)

  JVS(995) = B(860)

  JVS(996) = B(934)

  JVS(997) = 0

  JVS(998) = -B(330)-B(492)-B(594)-B(677)-B(932)

  JVS(999) = B(728)+B(730)

  JVS(1000) = 0

  JVS(1001) = -B(331)

  JVS(1002) = B(729)+B(861)

  JVS(1003) = -B(933)

  JVS(1004) = -B(493)+B(731)

  JVS(1005) = -B(595)

  JVS(1006) = B(242)

  JVS(1007) = -B(306)-B(428)-B(468)-B(572)-B(663)

  JVS(1008) = -B(307)-B(429)

  JVS(1009) = B(243)

  JVS(1010) = -B(469)

  JVS(1011) = -B(573)

  JVS(1012) = B(862)

  JVS(1013) = B(987)+B(1051)

  JVS(1014) = B(66)

  JVS(1015) = B(986)+B(1049)

  JVS(1016) = B(80)+B(900)

  JVS(1017) = B(985)+B(1036)

  JVS(1018) = B(47)+B(808)

  JVS(1019) = B(984)+B(1034)

  JVS(1020) = B(1016)+B(1024)+B(1032)

  JVS(1021) = B(754)

  JVS(1022) = B(1040)+B(1044)+B(1048)

  JVS(1023) = 0

  JVS(1024) = B(712)

  JVS(1025) = B(304)+B(466)+B(662)

  JVS(1026) = 0

  JVS(1027) = 0

  JVS(1028) = -B(25)-B(734)

  JVS(1029) = B(750)

  JVS(1030) = B(350)+B(512)+B(687)

  JVS(1031) = 0

  JVS(1032) = 0

  JVS(1033) = B(726)

  JVS(1034) = 0

  JVS(1035) = 0

  JVS(1036) = 0

  JVS(1037) = B(727)+B(751)

  JVS(1038) = B(305)+B(351)+B(1017)+B(1041)

  JVS(1039) = B(713)-B(735)+B(755)+B(809)+B(863)+B(901)+B(1035)+B(1037)+B(1050)+B(1052)

  JVS(1040) = 0

  JVS(1041) = B(467)+B(513)+B(1025)+B(1045)

  JVS(1042) = 0

  JVS(1043) = B(802)

  JVS(1044) = B(64)+B(870)

  JVS(1045) = B(340)+B(502)+B(682)

  JVS(1046) = -B(29)-B(744)-B(746)-B(748)-B(750)

  JVS(1047) = -B(749)-B(751)

  JVS(1048) = B(341)

  JVS(1049) = -B(745)+B(803)+B(871)

  JVS(1050) = B(503)-B(747)

  JVS(1051) = 0

  JVS(1052) = B(56)

  JVS(1053) = B(212)

  JVS(1054) = -B(324)-B(420)-B(486)-B(588)-B(674)

  JVS(1055) = B(294)+B(456)+B(666)

  JVS(1056) = B(213)

  JVS(1057) = B(295)-B(325)-B(421)

  JVS(1058) = 0

  JVS(1059) = B(457)-B(487)

  JVS(1060) = -B(589)

  JVS(1061) = B(792)

  JVS(1062) = B(67)

  JVS(1063) = B(880)

  JVS(1064) = B(66)

  JVS(1065) = B(29)+B(748)+B(750)

  JVS(1066) = -B(22)-B(728)-B(730)

  JVS(1067) = B(348)+B(510)+B(686)

  JVS(1068) = B(344)+B(506)+B(684)+B(685)

  JVS(1069) = B(312)+B(474)+B(668)

  JVS(1070) = B(749)+B(751)

  JVS(1071) = B(313)+B(345)+B(349)

  JVS(1072) = -B(729)+B(793)+B(881)

  JVS(1073) = B(475)+B(507)+B(511)-B(731)

  JVS(1074) = 0

  JVS(1075) = B(874)

  JVS(1076) = B(946)

  JVS(1077) = -B(348)-B(350)-B(352)-B(510)-B(512)-B(514)-B(608)-B(686)-B(687)-B(688)-B(944)

  JVS(1078) = B(720)+B(722)

  JVS(1079) = -B(349)-B(351)-B(353)

  JVS(1080) = B(721)+B(875)

  JVS(1081) = -B(945)

  JVS(1082) = -B(511)-B(513)-B(515)+B(723)

  JVS(1083) = -B(609)

  JVS(1084) = B(943)

  JVS(1085) = B(872)

  JVS(1086) = B(744)+B(746)

  JVS(1087) = -B(344)-B(346)-B(506)-B(508)-B(606)-B(684)-B(685)-B(941)

  JVS(1088) = 0

  JVS(1089) = -B(345)-B(347)

  JVS(1090) = B(745)+B(873)

  JVS(1091) = -B(942)

  JVS(1092) = -B(507)-B(509)+B(747)

  JVS(1093) = -B(607)

  JVS(1094) = B(976)

  JVS(1095) = B(978)

  JVS(1096) = B(974)

  JVS(1097) = B(972)

  JVS(1098) = B(266)

  JVS(1099) = B(956)

  JVS(1100) = 0

  JVS(1101) = B(176)

  JVS(1102) = 0

  JVS(1103) = 0

  JVS(1104) = 0

  JVS(1105) = 0

  JVS(1106) = 0

  JVS(1107) = 0

  JVS(1108) = 0

  JVS(1109) = -B(310)-B(312)-B(408)-B(472)-B(474)-B(576)-B(667)-B(668)

  JVS(1110) = 0

  JVS(1111) = 0

  JVS(1112) = 0

  JVS(1113) = 0

  JVS(1114) = -B(311)-B(313)-B(409)

  JVS(1115) = B(177)+B(267)+B(957)+B(973)+B(975)+B(977)+B(979)

  JVS(1116) = -B(473)-B(475)

  JVS(1117) = -B(577)

  JVS(1118) = B(1002)

  JVS(1119) = B(1006)

  JVS(1120) = -B(1018)-B(1022)-B(1026)-B(1030)-B(1033)

  JVS(1121) = -B(1019)-B(1023)

  JVS(1122) = B(1003)+B(1007)

  JVS(1123) = -B(1027)

  JVS(1124) = -B(1031)

  JVS(1125) = B(949)

  JVS(1126) = B(65)

  JVS(1127) = B(204)+B(206)

  JVS(1128) = B(342)+B(504)+B(683)

  JVS(1129) = 0

  JVS(1130) = -B(23)-B(720)-B(722)-B(724)-B(726)

  JVS(1131) = B(205)+B(207)-B(725)-B(727)

  JVS(1132) = B(343)

  JVS(1133) = -B(721)+B(950)

  JVS(1134) = 0

  JVS(1135) = B(505)-B(723)

  JVS(1136) = 0

  JVS(1137) = B(244)

  JVS(1138) = -B(308)-B(430)-B(470)-B(574)-B(664)

  JVS(1139) = -B(309)-B(431)

  JVS(1140) = B(245)

  JVS(1141) = -B(471)

  JVS(1142) = -B(575)

  JVS(1143) = B(858)

  JVS(1144) = B(931)

  JVS(1145) = -B(328)-B(490)-B(592)-B(676)-B(929)

  JVS(1146) = B(336)+B(498)+B(680)

  JVS(1147) = B(250)+B(256)

  JVS(1148) = -B(329)+B(337)

  JVS(1149) = B(251)+B(859)

  JVS(1150) = -B(930)

  JVS(1151) = B(257)-B(491)+B(499)

  JVS(1152) = -B(593)

  JVS(1153) = B(85)

  JVS(1154) = B(912)

  JVS(1155) = B(374)+B(534)+B(698)

  JVS(1156) = -B(376)-B(536)-B(630)-B(699)-B(990)

  JVS(1157) = 0

  JVS(1158) = B(375)-B(377)-B(991)

  JVS(1159) = B(913)

  JVS(1160) = B(535)-B(537)

  JVS(1161) = -B(631)

  JVS(1162) = B(384)

  JVS(1163) = B(366)

  JVS(1164) = -B(294)-B(406)-B(456)-B(562)-B(666)

  JVS(1165) = 0

  JVS(1166) = -B(295)+B(367)+B(385)-B(407)

  JVS(1167) = 0

  JVS(1168) = -B(457)

  JVS(1169) = -B(563)

  JVS(1170) = B(262)

  JVS(1171) = 2*B(75)+2*B(890)

  JVS(1172) = B(76)+B(892)

  JVS(1173) = B(59)

  JVS(1174) = B(37)+B(778)

  JVS(1175) = B(814)

  JVS(1176) = B(49)+B(842)

  JVS(1177) = 2*B(70)

  JVS(1178) = B(1065)+2*B(1067)

  JVS(1179) = B(939)

  JVS(1180) = B(65)

  JVS(1181) = B(71)

  JVS(1182) = B(87)

  JVS(1183) = B(89)

  JVS(1184) = B(92)

  JVS(1185) = B(935)

  JVS(1186) = B(192)+B(194)

  JVS(1187) = 2*B(354)+2*B(516)+2*B(689)

  JVS(1188) = B(188)+B(190)

  JVS(1189) = B(224)+B(226)

  JVS(1190) = B(31)

  JVS(1191) = B(356)+B(518)+B(690)

  JVS(1192) = B(386)+B(548)+B(705)

  JVS(1193) = 0

  JVS(1194) = B(342)+B(504)+B(683)

  JVS(1195) = B(18)

  JVS(1196) = B(382)+B(542)+B(702)

  JVS(1197) = 0

  JVS(1198) = B(378)+B(538)+B(700)

  JVS(1199) = B(330)+B(492)+B(677)

  JVS(1200) = B(22)

  JVS(1201) = B(350)+B(352)+B(512)+B(514)+B(687)+B(688)

  JVS(1202) = 0

  JVS(1203) = 2*B(310)+2*B(472)+2*B(667)

  JVS(1204) = 0

  JVS(1205) = B(23)+B(724)+B(726)

  JVS(1206) = 0

  JVS(1207) = 0

  JVS(1208) = -B(9)-B(10)-B(246)-B(252)

  JVS(1209) = 0

  JVS(1210) = B(280)+B(442)+B(644)+B(645)

  JVS(1211) = 0

  JVS(1212) = B(390)

  JVS(1213) = B(314)+B(476)+B(669)

  JVS(1214) = 0

  JVS(1215) = B(332)+B(494)+B(678)

  JVS(1216) = B(189)+B(191)+B(193)+B(195)+B(225)+B(227)+B(725)+B(727)

  JVS(1217) = B(281)+2*B(311)+B(315)+B(331)+B(333)+B(343)+B(351)+B(353)+2*B(355)+B(357)+B(379)+B(383)+B(387)+B(391)&
                &+B(1066)

  JVS(1218) = -B(247)+B(263)+B(779)+B(815)+B(843)+2*B(891)+B(893)+B(936)+B(940)

  JVS(1219) = 0

  JVS(1220) = -B(253)+B(443)+2*B(473)+B(477)+B(493)+B(495)+B(505)+B(513)+B(515)+2*B(517)+B(519)+B(539)+B(543)+B(549)

  JVS(1221) = 0

  JVS(1222) = B(740)

  JVS(1223) = B(43)

  JVS(1224) = 0

  JVS(1225) = B(298)+B(460)+B(659)

  JVS(1226) = -B(318)-B(414)-B(480)-B(582)-B(671)

  JVS(1227) = 0

  JVS(1228) = 0

  JVS(1229) = B(299)-B(319)-B(415)

  JVS(1230) = B(741)

  JVS(1231) = B(461)-B(481)

  JVS(1232) = -B(583)

  JVS(1233) = B(260)

  JVS(1234) = B(710)

  JVS(1235) = 0

  JVS(1236) = -B(334)-B(496)-B(598)-B(679)

  JVS(1237) = 0

  JVS(1238) = 0

  JVS(1239) = -B(335)

  JVS(1240) = B(261)+B(711)

  JVS(1241) = -B(497)

  JVS(1242) = -B(599)

  JVS(1243) = B(168)

  JVS(1244) = B(970)

  JVS(1245) = B(14)

  JVS(1246) = B(60)

  JVS(1247) = B(58)

  JVS(1248) = B(41)

  JVS(1249) = B(53)

  JVS(1250) = B(17)

  JVS(1251) = B(32)+B(756)

  JVS(1252) = 0

  JVS(1253) = B(95)+B(962)

  JVS(1254) = 0

  JVS(1255) = 0

  JVS(1256) = 0

  JVS(1257) = 0

  JVS(1258) = 0

  JVS(1259) = B(328)+B(490)+B(676)

  JVS(1260) = -B(282)-B(396)-B(444)-B(552)-B(647)-B(648)-B(649)

  JVS(1261) = 0

  JVS(1262) = 0

  JVS(1263) = 0

  JVS(1264) = B(288)+B(450)+B(656)

  JVS(1265) = B(12)

  JVS(1266) = -B(283)+B(289)+B(329)-B(397)

  JVS(1267) = B(169)+B(757)+B(963)+B(971)

  JVS(1268) = 0

  JVS(1269) = -B(445)+B(451)+B(491)

  JVS(1270) = -B(553)

  JVS(1271) = B(776)

  JVS(1272) = 0

  JVS(1273) = -B(392)-B(440)-B(546)-B(640)-B(704)

  JVS(1274) = 0

  JVS(1275) = 0

  JVS(1276) = -B(393)-B(441)

  JVS(1277) = B(777)

  JVS(1278) = -B(547)

  JVS(1279) = -B(641)

  JVS(1280) = B(964)

  JVS(1281) = B(966)

  JVS(1282) = B(968)

  JVS(1283) = B(166)

  JVS(1284) = B(98)

  JVS(1285) = 0.95*B(1084)

  JVS(1286) = B(278)

  JVS(1287) = B(1083)

  JVS(1288) = B(57)

  JVS(1289) = B(840)

  JVS(1290) = 0

  JVS(1291) = B(192)

  JVS(1292) = B(196)

  JVS(1293) = B(1076)

  JVS(1294) = 0.18*B(1072)

  JVS(1295) = B(94)+B(960)

  JVS(1296) = B(13)

  JVS(1297) = 0

  JVS(1298) = 0

  JVS(1299) = 0

  JVS(1300) = 0

  JVS(1301) = 0

  JVS(1302) = -B(280)-B(394)-B(442)-B(550)-B(644)-B(645)-B(646)-B(980)

  JVS(1303) = 0

  JVS(1304) = 0

  JVS(1305) = B(326)+B(488)+B(675)

  JVS(1306) = 0

  JVS(1307) = B(193)+B(197)

  JVS(1308) = -B(281)+B(327)-B(395)

  JVS(1309) = B(167)+B(279)+B(841)+B(961)+B(965)+B(967)+B(969)+0.95*B(1085)

  JVS(1310) = -B(981)+0.18*B(1073)

  JVS(1311) = -B(443)+B(489)

  JVS(1312) = B(11)

  JVS(1313) = -B(551)

  JVS(1314) = B(772)

  JVS(1315) = B(83)

  JVS(1316) = B(908)

  JVS(1317) = 0

  JVS(1318) = -B(372)-B(434)-B(532)-B(626)-B(697)

  JVS(1319) = B(370)+B(530)+B(696)

  JVS(1320) = 0

  JVS(1321) = B(371)-B(373)-B(435)

  JVS(1322) = B(773)+B(909)

  JVS(1323) = B(531)-B(533)

  JVS(1324) = -B(627)

  JVS(1325) = B(742)

  JVS(1326) = B(738)

  JVS(1327) = 0

  JVS(1328) = 0

  JVS(1329) = 0

  JVS(1330) = B(180)

  JVS(1331) = B(296)+B(458)+B(658)

  JVS(1332) = 0

  JVS(1333) = -B(316)-B(412)-B(478)-B(580)-B(670)

  JVS(1334) = 0

  JVS(1335) = 0

  JVS(1336) = 0

  JVS(1337) = B(297)-B(317)-B(413)

  JVS(1338) = B(181)+B(739)+B(743)

  JVS(1339) = B(459)-B(479)

  JVS(1340) = -B(581)

  JVS(1341) = B(774)

  JVS(1342) = B(222)

  JVS(1343) = 0

  JVS(1344) = 0

  JVS(1345) = -B(388)-B(390)-B(438)-B(544)-B(638)-B(703)

  JVS(1346) = B(223)

  JVS(1347) = -B(389)-B(391)-B(439)

  JVS(1348) = B(775)

  JVS(1349) = -B(545)

  JVS(1350) = -B(639)

  JVS(1351) = B(172)

  JVS(1352) = B(61)

  JVS(1353) = B(16)

  JVS(1354) = B(15)

  JVS(1355) = B(988)+B(1055)

  JVS(1356) = 0

  JVS(1357) = 0

  JVS(1358) = B(33)+B(758)

  JVS(1359) = 0

  JVS(1360) = 0

  JVS(1361) = 0

  JVS(1362) = 0

  JVS(1363) = 0

  JVS(1364) = 0

  JVS(1365) = -B(284)-B(398)-B(446)-B(554)-B(650)-B(651)-B(652)

  JVS(1366) = 0

  JVS(1367) = 0

  JVS(1368) = 0

  JVS(1369) = 0

  JVS(1370) = -B(285)-B(399)

  JVS(1371) = B(173)+B(759)+B(1056)

  JVS(1372) = -B(447)

  JVS(1373) = -B(555)

  JVS(1374) = B(706)

  JVS(1375) = 0

  JVS(1376) = -B(336)-B(498)-B(600)-B(680)

  JVS(1377) = 0

  JVS(1378) = 0

  JVS(1379) = -B(337)

  JVS(1380) = B(707)

  JVS(1381) = -B(499)

  JVS(1382) = -B(601)

  JVS(1383) = B(274)

  JVS(1384) = B(270)

  JVS(1385) = B(982)

  JVS(1386) = B(178)

  JVS(1387) = B(736)

  JVS(1388) = 0

  JVS(1389) = 0

  JVS(1390) = 0

  JVS(1391) = 0

  JVS(1392) = 0

  JVS(1393) = 0

  JVS(1394) = 0

  JVS(1395) = 0

  JVS(1396) = 0

  JVS(1397) = 0

  JVS(1398) = 0

  JVS(1399) = -B(314)-B(410)-B(476)-B(578)-B(669)

  JVS(1400) = 0

  JVS(1401) = 0

  JVS(1402) = 0

  JVS(1403) = 0

  JVS(1404) = 0

  JVS(1405) = -B(315)-B(411)

  JVS(1406) = B(179)+B(271)+B(275)+B(737)

  JVS(1407) = -B(477)

  JVS(1408) = -B(579)

  JVS(1409) = B(708)

  JVS(1410) = 0

  JVS(1411) = -B(338)-B(500)-B(602)-B(681)

  JVS(1412) = 0

  JVS(1413) = -B(339)

  JVS(1414) = B(709)

  JVS(1415) = -B(501)

  JVS(1416) = -B(603)

  JVS(1417) = B(14)

  JVS(1418) = B(27)

  JVS(1419) = B(856)

  JVS(1420) = B(67)

  JVS(1421) = 0

  JVS(1422) = B(928)

  JVS(1423) = B(15)

  JVS(1424) = B(19)

  JVS(1425) = 2*B(26)

  JVS(1426) = B(17)

  JVS(1427) = 2*B(36)

  JVS(1428) = B(31)

  JVS(1429) = 0

  JVS(1430) = 0

  JVS(1431) = 0

  JVS(1432) = B(18)

  JVS(1433) = B(13)

  JVS(1434) = 0

  JVS(1435) = 0

  JVS(1436) = B(25)+B(734)

  JVS(1437) = B(29)+B(748)

  JVS(1438) = B(348)+B(510)+B(686)

  JVS(1439) = B(344)+B(506)+B(684)

  JVS(1440) = 0

  JVS(1441) = B(23)+B(724)

  JVS(1442) = 0

  JVS(1443) = 0

  JVS(1444) = B(334)+B(496)+B(679)

  JVS(1445) = 0

  JVS(1446) = 0

  JVS(1447) = 0

  JVS(1448) = 0

  JVS(1449) = 0

  JVS(1450) = 0

  JVS(1451) = 0

  JVS(1452) = -B(326)-B(488)-B(590)-B(675)-B(926)

  JVS(1453) = B(332)+B(494)+B(678)

  JVS(1454) = 0

  JVS(1455) = B(725)+B(749)

  JVS(1456) = 0

  JVS(1457) = 0

  JVS(1458) = -B(327)+B(333)+B(335)+B(345)+B(349)

  JVS(1459) = B(248)+B(735)+B(857)

  JVS(1460) = -B(927)

  JVS(1461) = B(254)-B(489)+B(495)+B(497)+B(507)+B(511)

  JVS(1462) = B(249)+B(255)

  JVS(1463) = -B(591)

  JVS(1464) = 2*B(28)

  JVS(1465) = B(27)

  JVS(1466) = 0

  JVS(1467) = B(88)

  JVS(1468) = 0

  JVS(1469) = B(16)

  JVS(1470) = B(20)

  JVS(1471) = B(21)

  JVS(1472) = B(380)+B(540)+B(701)

  JVS(1473) = B(258)

  JVS(1474) = 0

  JVS(1475) = 0

  JVS(1476) = 0

  JVS(1477) = 0

  JVS(1478) = 0

  JVS(1479) = 0

  JVS(1480) = 0

  JVS(1481) = 0

  JVS(1482) = 0

  JVS(1483) = 0

  JVS(1484) = -B(332)-B(494)-B(596)-B(678)

  JVS(1485) = 0

  JVS(1486) = 0

  JVS(1487) = 0

  JVS(1488) = 0

  JVS(1489) = -B(333)+B(381)

  JVS(1490) = B(259)

  JVS(1491) = 0

  JVS(1492) = -B(495)+B(541)

  JVS(1493) = -B(597)

  JVS(1494) = B(904)

  JVS(1495) = B(953)

  JVS(1496) = B(35)+B(762)+B(768)

  JVS(1497) = 0

  JVS(1498) = -B(370)-B(530)-B(624)-B(696)-B(951)

  JVS(1499) = 0

  JVS(1500) = -B(371)

  JVS(1501) = B(763)+B(905)

  JVS(1502) = -B(952)

  JVS(1503) = -B(531)+B(769)

  JVS(1504) = -B(625)

  JVS(1505) = -B(1068)

  JVS(1506) = -B(192)-B(194)

  JVS(1507) = -B(188)-B(190)

  JVS(1508) = -B(222)-B(224)-B(226)-B(228)

  JVS(1509) = -B(196)-B(198)

  JVS(1510) = -B(212)-B(214)-B(216)

  JVS(1511) = -B(1079)

  JVS(1512) = -B(1074)

  JVS(1513) = B(96)-B(99)

  JVS(1514) = -B(204)-B(206)

  JVS(1515) = -B(748)-B(750)

  JVS(1516) = -B(724)-B(726)

  JVS(1517) = -B(1)-B(2)-B(100)-B(109)-B(111)-B(124)-B(132)-B(189)-B(191)-B(193)-B(195)-B(197)-B(199)-B(205)-B(207)&
                &-B(213)-B(215)-B(217)-B(223)-B(225)-B(227)-B(229)-B(725)-B(727)-B(749)-B(751)-B(1069)-B(1075)-B(1080)

  JVS(1518) = -B(110)

  JVS(1519) = -B(125)

  JVS(1520) = -B(112)

  JVS(1521) = 0

  JVS(1522) = -B(133)

  JVS(1523) = B(63)

  JVS(1524) = B(989)+B(1053)

  JVS(1525) = 0

  JVS(1526) = B(718)

  JVS(1527) = 0

  JVS(1528) = 0

  JVS(1529) = 0

  JVS(1530) = 0

  JVS(1531) = 0

  JVS(1532) = -B(292)-B(404)-B(454)-B(560)-B(665)

  JVS(1533) = -B(293)-B(405)

  JVS(1534) = B(719)+B(1054)

  JVS(1535) = 0

  JVS(1536) = -B(455)

  JVS(1537) = -B(561)

  JVS(1538) = B(764)

  JVS(1539) = B(55)

  JVS(1540) = 0

  JVS(1541) = 0

  JVS(1542) = 0

  JVS(1543) = 0

  JVS(1544) = 0

  JVS(1545) = 0

  JVS(1546) = B(292)+B(454)+B(665)

  JVS(1547) = -B(322)-B(418)-B(484)-B(586)-B(673)

  JVS(1548) = B(293)-B(323)-B(419)

  JVS(1549) = B(765)

  JVS(1550) = 0

  JVS(1551) = B(455)-B(485)

  JVS(1552) = -B(587)

  JVS(1553) = B(174)

  JVS(1554) = B(91)

  JVS(1555) = B(62)

  JVS(1556) = B(34)+B(760)

  JVS(1557) = B(716)

  JVS(1558) = 0

  JVS(1559) = 0

  JVS(1560) = 0

  JVS(1561) = 0

  JVS(1562) = B(392)+B(546)+B(704)

  JVS(1563) = 0

  JVS(1564) = 0

  JVS(1565) = 0

  JVS(1566) = 0

  JVS(1567) = 0

  JVS(1568) = 0

  JVS(1569) = -B(288)-B(290)-B(402)-B(450)-B(452)-B(558)-B(656)-B(657)

  JVS(1570) = -B(289)-B(291)+B(393)-B(403)

  JVS(1571) = B(175)+B(717)+B(761)

  JVS(1572) = 0

  JVS(1573) = -B(451)-B(453)+B(547)

  JVS(1574) = -B(559)

  JVS(1575) = B(268)

  JVS(1576) = B(51)+B(846)

  JVS(1577) = B(39)+B(782)

  JVS(1578) = 2*B(74)

  JVS(1579) = B(73)

  JVS(1580) = 2*B(320)+2*B(482)+2*B(672)

  JVS(1581) = B(21)

  JVS(1582) = 0

  JVS(1583) = 0

  JVS(1584) = B(318)+B(480)+B(671)

  JVS(1585) = B(284)+B(446)+B(650)+B(651)

  JVS(1586) = 0

  JVS(1587) = 0

  JVS(1588) = 0

  JVS(1589) = 0

  JVS(1590) = 0

  JVS(1591) = 0

  JVS(1592) = -B(12)-B(250)-B(256)

  JVS(1593) = B(285)+B(319)+2*B(321)

  JVS(1594) = -B(251)+B(269)+B(783)+B(847)

  JVS(1595) = 0

  JVS(1596) = -B(257)+B(447)+B(481)+2*B(483)

  JVS(1597) = 0

  JVS(1598) = B(7)

  JVS(1599) = -B(1065)

  JVS(1600) = B(1070)

  JVS(1601) = -B(360)

  JVS(1602) = -B(358)

  JVS(1603) = -B(368)

  JVS(1604) = -B(354)

  JVS(1605) = 0

  JVS(1606) = -B(1016)-B(1020)

  JVS(1607) = -B(286)-B(400)

  JVS(1608) = -B(374)

  JVS(1609) = 0

  JVS(1610) = -B(300)-B(424)

  JVS(1611) = 0

  JVS(1612) = -B(356)

  JVS(1613) = -B(386)

  JVS(1614) = -B(320)-B(416)

  JVS(1615) = B(1077)

  JVS(1616) = B(1072)

  JVS(1617) = -B(101)+B(103)

  JVS(1618) = -B(1040)-B(1042)

  JVS(1619) = 0

  JVS(1620) = -B(296)

  JVS(1621) = -B(298)

  JVS(1622) = -B(340)-B(342)-B(422)

  JVS(1623) = -B(302)-B(304)-B(426)

  JVS(1624) = -B(380)

  JVS(1625) = -B(382)-B(384)-B(436)

  JVS(1626) = -B(364)-B(366)-B(432)

  JVS(1627) = -B(362)

  JVS(1628) = -B(378)

  JVS(1629) = -B(330)

  JVS(1630) = -B(306)-B(428)

  JVS(1631) = -B(324)-B(420)

  JVS(1632) = 0

  JVS(1633) = -B(348)-B(350)-B(352)

  JVS(1634) = -B(344)-B(346)

  JVS(1635) = -B(310)-B(312)-B(408)

  JVS(1636) = -B(1018)-B(1022)

  JVS(1637) = 0

  JVS(1638) = -B(308)-B(430)

  JVS(1639) = -B(328)

  JVS(1640) = -B(376)-B(990)

  JVS(1641) = -B(294)-B(406)

  JVS(1642) = -B(318)-B(414)

  JVS(1643) = -B(334)

  JVS(1644) = -B(282)-B(396)

  JVS(1645) = -B(392)-B(440)

  JVS(1646) = -B(280)-B(394)

  JVS(1647) = -B(372)-B(434)

  JVS(1648) = -B(316)-B(412)

  JVS(1649) = -B(388)-B(390)-B(438)

  JVS(1650) = -B(284)-B(398)

  JVS(1651) = -B(336)

  JVS(1652) = -B(314)-B(410)

  JVS(1653) = -B(338)

  JVS(1654) = -B(326)

  JVS(1655) = -B(332)

  JVS(1656) = -B(370)

  JVS(1657) = -B(109)

  JVS(1658) = -B(292)-B(404)

  JVS(1659) = -B(322)-B(418)

  JVS(1660) = -B(288)-B(290)-B(402)

  JVS(1661) = 0

  JVS(1662) = -B(102)-B(110)-2*B(113)-B(115)-B(139)-B(145)-B(281)-B(283)-B(285)-B(287)-B(289)-B(291)-B(293)-B(295)&
                &-B(297)-B(299)-B(301)-B(303)-B(305)-B(307)-B(309)-B(311)-B(313)-B(315)-B(317)-B(319)-B(321)-B(323)-B(325)&
                &-B(327)-B(329)-B(331)-B(333)-B(335)-B(337)-B(339)-B(341)-B(343)-B(345)-B(347)-B(349)-B(351)-B(353)-B(355)&
                &-B(357)-B(359)-B(361)-B(363)-B(365)-B(367)-B(369)-B(371)-B(373)-B(375)-B(377)-B(379)-B(381)-B(383)-B(385)&
                &-B(387)-B(389)-B(391)-B(393)-B(395)-B(397)-B(399)-B(401)-B(403)-B(405)-B(407)-B(409)-B(411)-B(413)-B(415)&
                &-B(417)-B(419)-B(421)-B(423)-B(425)-B(427)-B(429)-B(431)-B(433)-B(435)-B(437)-B(439)-B(441)-B(991)-B(1017)&
                &-B(1019)-B(1021)-B(1023)-B(1041)-B(1043)-B(1066)

  JVS(1663) = -B(140)

  JVS(1664) = B(4)+B(104)+B(117)+B(1071)+B(1073)+B(1078)

  JVS(1665) = B(5)-B(116)+B(118)

  JVS(1666) = 0

  JVS(1667) = -B(146)

  JVS(1668) = 2*B(122)

  JVS(1669) = -B(168)

  JVS(1670) = -B(174)

  JVS(1671) = -B(964)

  JVS(1672) = -B(966)

  JVS(1673) = -B(968)

  JVS(1674) = -B(970)

  JVS(1675) = -B(976)

  JVS(1676) = -B(978)

  JVS(1677) = -B(974)

  JVS(1678) = -B(972)

  JVS(1679) = -B(166)

  JVS(1680) = -B(1008)

  JVS(1681) = -B(1010)

  JVS(1682) = -B(1012)

  JVS(1683) = -B(1014)

  JVS(1684) = -B(1038)

  JVS(1685) = -B(170)-B(172)

  JVS(1686) = -B(242)-B(244)

  JVS(1687) = -B(1000)-B(1002)

  JVS(1688) = -B(1004)-B(1006)

  JVS(1689) = -B(234)-B(236)

  JVS(1690) = -B(238)-B(240)

  JVS(1691) = -B(260)

  JVS(1692) = -B(764)

  JVS(1693) = -B(230)

  JVS(1694) = -B(262)

  JVS(1695) = -B(740)

  JVS(1696) = -B(742)

  JVS(1697) = -B(772)

  JVS(1698) = B(7)-B(154)

  JVS(1699) = -B(1084)

  JVS(1700) = -B(264)-B(266)

  JVS(1701) = -B(272)-B(274)

  JVS(1702) = -B(738)

  JVS(1703) = -B(276)

  JVS(1704) = -B(820)

  JVS(1705) = -B(830)

  JVS(1706) = -B(268)-B(270)

  JVS(1707) = -B(278)

  JVS(1708) = B(57)-B(856)

  JVS(1709) = -B(794)

  JVS(1710) = -B(796)

  JVS(1711) = B(51)

  JVS(1712) = B(56)

  JVS(1713) = B(60)

  JVS(1714) = B(61)

  JVS(1715) = B(75)

  JVS(1716) = B(77)

  JVS(1717) = B(58)-B(858)

  JVS(1718) = -B(802)

  JVS(1719) = B(68)

  JVS(1720) = -B(792)

  JVS(1721) = -B(810)

  JVS(1722) = B(81)

  JVS(1723) = B(83)-B(904)

  JVS(1724) = B(82)

  JVS(1725) = B(84)-B(908)

  JVS(1726) = B(85)-B(910)

  JVS(1727) = B(90)

  JVS(1728) = B(50)

  JVS(1729) = -B(150)

  JVS(1730) = B(52)

  JVS(1731) = -B(956)

  JVS(1732) = -B(994)

  JVS(1733) = -B(800)

  JVS(1734) = -B(788)

  JVS(1735) = B(55)

  JVS(1736) = B(63)

  JVS(1737) = B(72)

  JVS(1738) = B(76)

  JVS(1739) = B(59)-B(860)

  JVS(1740) = -B(784)

  JVS(1741) = B(67)-B(874)

  JVS(1742) = -B(949)

  JVS(1743) = B(69)

  JVS(1744) = B(78)

  JVS(1745) = -B(778)

  JVS(1746) = B(79)

  JVS(1747) = -B(780)

  JVS(1748) = B(86)-B(912)

  JVS(1749) = B(88)-B(916)

  JVS(1750) = -B(812)

  JVS(1751) = -B(814)

  JVS(1752) = -B(816)

  JVS(1753) = -B(782)

  JVS(1754) = -B(818)

  JVS(1755) = B(91)

  JVS(1756) = B(49)-B(840)

  JVS(1757) = -B(798)

  JVS(1758) = -B(790)

  JVS(1759) = B(987)

  JVS(1760) = B(70)

  JVS(1761) = 2*B(3)-B(130)

  JVS(1762) = -B(776)

  JVS(1763) = -B(1059)-B(1061)

  JVS(1764) = -B(1086)

  JVS(1765) = -B(937)

  JVS(1766) = -B(939)

  JVS(1767) = B(64)+B(65)

  JVS(1768) = B(66)-B(872)

  JVS(1769) = B(62)

  JVS(1770) = -B(947)

  JVS(1771) = B(71)

  JVS(1772) = -B(804)

  JVS(1773) = -B(826)

  JVS(1774) = -B(806)

  JVS(1775) = -B(836)

  JVS(1776) = -B(954)

  JVS(1777) = B(87)-B(914)

  JVS(1778) = B(89)-B(918)

  JVS(1779) = B(92)-B(924)

  JVS(1780) = -B(935)

  JVS(1781) = -B(160)

  JVS(1782) = -B(1049)

  JVS(1783) = B(74)

  JVS(1784) = B(80)

  JVS(1785) = -B(708)

  JVS(1786) = -B(774)

  JVS(1787) = B(985)

  JVS(1788) = -B(786)

  JVS(1789) = -B(996)

  JVS(1790) = B(53)+B(54)

  JVS(1791) = B(73)

  JVS(1792) = -B(808)

  JVS(1793) = -B(706)

  JVS(1794) = -B(714)

  JVS(1795) = -B(1034)

  JVS(1796) = -B(178)+B(192)

  JVS(1797) = 0

  JVS(1798) = 0

  JVS(1799) = 0

  JVS(1800) = 0

  JVS(1801) = -B(736)

  JVS(1802) = -B(710)

  JVS(1803) = -B(732)

  JVS(1804) = -B(752)

  JVS(1805) = -B(1053)-B(1055)

  JVS(1806) = -B(176)+B(188)

  JVS(1807) = -B(756)-B(958)

  JVS(1808) = 0

  JVS(1809) = 0

  JVS(1810) = 0

  JVS(1811) = -B(766)

  JVS(1812) = -B(218)+B(222)

  JVS(1813) = -B(992)

  JVS(1814) = -B(754)

  JVS(1815) = -B(180)+B(196)

  JVS(1816) = 0

  JVS(1817) = -B(208)+B(212)

  JVS(1818) = -B(762)

  JVS(1819) = 0

  JVS(1820) = -B(760)-B(962)

  JVS(1821) = 0

  JVS(1822) = -B(716)

  JVS(1823) = 0

  JVS(1824) = 0

  JVS(1825) = 0

  JVS(1826) = 0

  JVS(1827) = 0

  JVS(1828) = -B(200)+B(204)

  JVS(1829) = 0

  JVS(1830) = 0

  JVS(1831) = 0

  JVS(1832) = -B(712)

  JVS(1833) = -B(718)

  JVS(1834) = 0

  JVS(1835) = -B(758)-B(960)

  JVS(1836) = 0

  JVS(1837) = 0

  JVS(1838) = B(8)-B(156)

  JVS(1839) = 0

  JVS(1840) = -B(258)

  JVS(1841) = 0

  JVS(1842) = -B(128)

  JVS(1843) = 0

  JVS(1844) = 0

  JVS(1845) = 0

  JVS(1846) = -B(734)

  JVS(1847) = -B(744)+B(748)

  JVS(1848) = 0

  JVS(1849) = -B(728)

  JVS(1850) = 0

  JVS(1851) = 0

  JVS(1852) = 0

  JVS(1853) = 0

  JVS(1854) = -B(720)+B(724)

  JVS(1855) = 0

  JVS(1856) = 0

  JVS(1857) = 0

  JVS(1858) = 0

  JVS(1859) = -B(246)

  JVS(1860) = 0

  JVS(1861) = 0

  JVS(1862) = 0

  JVS(1863) = 0

  JVS(1864) = 0

  JVS(1865) = 0

  JVS(1866) = 0

  JVS(1867) = 0

  JVS(1868) = 0

  JVS(1869) = 0

  JVS(1870) = 0

  JVS(1871) = 0

  JVS(1872) = 0

  JVS(1873) = 0

  JVS(1874) = 0

  JVS(1875) = -B(124)+B(132)+B(189)+B(193)+B(197)+B(205)+B(213)+B(223)+B(725)+B(749)

  JVS(1876) = 0

  JVS(1877) = 0

  JVS(1878) = 0

  JVS(1879) = -B(250)

  JVS(1880) = -B(139)+B(145)

  JVS(1881) = -B(125)-B(126)-B(129)-B(131)-B(134)-B(140)-B(141)-B(143)-B(151)-B(155)-B(157)-B(161)-B(167)-B(169)-B(171)&
                &-B(173)-B(175)-B(177)-B(179)-B(181)-B(201)-B(209)-B(219)-B(231)-B(235)-B(237)-B(239)-B(241)-B(243)-B(245)&
                &-B(247)-B(248)-B(251)-B(259)-B(261)-B(263)-B(265)-B(267)-B(269)-B(271)-B(273)-B(275)-B(277)-B(279)-B(707)&
                &-B(709)-B(711)-B(713)-B(715)-B(717)-B(719)-B(721)-B(729)-B(733)-B(735)-B(737)-B(739)-B(741)-B(743)-B(745)&
                &-B(753)-B(755)-B(757)-B(759)-B(761)-B(763)-B(765)-B(767)-B(773)-B(775)-B(777)-B(779)-B(781)-B(783)-B(785)&
                &-B(787)-B(789)-B(791)-B(793)-B(795)-B(797)-B(799)-B(801)-B(803)-B(805)-B(807)-B(809)-B(811)-B(813)-B(815)&
                &-B(817)-B(819)-B(821)-B(827)-B(831)-B(837)-B(841)-B(857)-B(859)-B(861)-B(873)-B(875)-B(905)-B(909)-B(911)&
                &-B(913)-B(915)-B(917)-B(919)-B(925)-B(936)-B(938)-B(940)-B(948)-B(950)-B(955)-B(957)-B(959)-B(961)-B(963)&
                &-B(965)-B(967)-B(969)-B(971)-B(973)-B(975)-B(977)-B(979)-B(993)-B(995)-B(997)-B(1001)-B(1003)-B(1005)&
                &-B(1007)-B(1009)-B(1011)-B(1013)-B(1015)-B(1035)-B(1039)-B(1050)-B(1054)-B(1056)-B(1060)-B(1062)-B(1085)&
                &-B(1087)

  JVS(1882) = -B(142)

  JVS(1883) = -B(144)+B(152)

  JVS(1884) = -B(249)

  JVS(1885) = B(133)-B(135)+B(146)+B(153)

  JVS(1886) = B(121)

  JVS(1887) = B(98)

  JVS(1888) = B(154)

  JVS(1889) = -B(824)

  JVS(1890) = -B(834)

  JVS(1891) = B(794)

  JVS(1892) = B(796)

  JVS(1893) = B(75)+B(890)

  JVS(1894) = B(77)+B(894)

  JVS(1895) = B(802)

  JVS(1896) = B(792)

  JVS(1897) = B(810)

  JVS(1898) = B(82)+B(906)

  JVS(1899) = B(149)+B(150)

  JVS(1900) = B(800)

  JVS(1901) = B(43)+B(788)

  JVS(1902) = B(76)+B(892)

  JVS(1903) = B(40)+B(784)

  JVS(1904) = B(943)+B(949)

  JVS(1905) = B(37)+B(778)

  JVS(1906) = B(38)+B(780)

  JVS(1907) = B(812)

  JVS(1908) = B(814)

  JVS(1909) = B(48)+B(816)

  JVS(1910) = B(39)+B(782)

  JVS(1911) = B(818)

  JVS(1912) = B(798)

  JVS(1913) = B(44)+B(790)

  JVS(1914) = B(1065)

  JVS(1915) = B(931)+B(937)

  JVS(1916) = B(934)+B(939)

  JVS(1917) = B(946)+B(947)

  JVS(1918) = B(45)+B(804)

  JVS(1919) = B(826)+B(828)

  JVS(1920) = B(46)+B(806)

  JVS(1921) = B(836)+B(838)

  JVS(1922) = B(953)+B(954)

  JVS(1923) = B(92)

  JVS(1924) = B(928)+B(935)

  JVS(1925) = B(986)+B(1049)

  JVS(1926) = B(41)+B(42)+B(786)

  JVS(1927) = B(47)+B(808)

  JVS(1928) = B(984)+B(1034)

  JVS(1929) = -B(1070)

  JVS(1930) = B(360)+B(522)

  JVS(1931) = 2*B(358)+2*B(520)+B(691)

  JVS(1932) = 2*B(368)+2*B(528)+B(695)

  JVS(1933) = 2*B(354)+2*B(516)+B(689)

  JVS(1934) = 0

  JVS(1935) = B(1016)+B(1024)

  JVS(1936) = B(286)+B(448)

  JVS(1937) = B(374)+B(534)

  JVS(1938) = B(992)

  JVS(1939) = B(31)+B(754)

  JVS(1940) = 0

  JVS(1941) = B(300)+B(462)

  JVS(1942) = 0

  JVS(1943) = 2*B(356)+2*B(518)+B(690)

  JVS(1944) = 2*B(386)+2*B(548)+B(705)

  JVS(1945) = B(320)+B(482)

  JVS(1946) = -B(1077)

  JVS(1947) = -B(1072)

  JVS(1948) = B(101)-B(103)-B(105)

  JVS(1949) = B(1040)+B(1044)

  JVS(1950) = 0

  JVS(1951) = B(296)+B(458)

  JVS(1952) = B(298)+B(460)

  JVS(1953) = B(340)+B(342)+B(502)+B(504)

  JVS(1954) = B(302)+B(304)+B(464)+B(466)

  JVS(1955) = B(380)+B(540)

  JVS(1956) = B(382)+B(384)+B(542)

  JVS(1957) = B(8)

  JVS(1958) = B(364)+B(366)+B(526)

  JVS(1959) = B(362)+B(524)

  JVS(1960) = B(378)+B(538)

  JVS(1961) = B(330)+B(492)-B(932)

  JVS(1962) = B(306)+B(468)

  JVS(1963) = 0

  JVS(1964) = B(324)+B(486)

  JVS(1965) = 0

  JVS(1966) = B(348)+B(350)+B(352)+B(510)+B(512)+B(514)-B(944)

  JVS(1967) = B(344)+B(346)+B(506)+B(508)-B(941)

  JVS(1968) = B(310)+B(312)+B(472)+B(474)

  JVS(1969) = B(1018)+B(1026)

  JVS(1970) = 0

  JVS(1971) = B(308)+B(470)

  JVS(1972) = B(328)+B(490)-B(929)

  JVS(1973) = B(376)+B(536)

  JVS(1974) = B(294)+B(456)

  JVS(1975) = 0

  JVS(1976) = B(318)+B(480)

  JVS(1977) = B(334)+B(496)

  JVS(1978) = B(282)+B(444)

  JVS(1979) = B(392)+B(546)

  JVS(1980) = B(280)+B(442)-B(980)

  JVS(1981) = B(372)+B(532)

  JVS(1982) = B(316)+B(478)

  JVS(1983) = B(388)+B(390)+B(544)

  JVS(1984) = B(284)+B(446)

  JVS(1985) = B(336)+B(498)

  JVS(1986) = B(314)+B(476)

  JVS(1987) = B(338)+B(500)

  JVS(1988) = B(326)+B(488)-B(926)

  JVS(1989) = B(332)+B(494)

  JVS(1990) = B(370)+B(530)-B(951)

  JVS(1991) = B(109)-B(111)

  JVS(1992) = B(292)+B(454)

  JVS(1993) = B(322)+B(484)

  JVS(1994) = B(288)+B(290)+B(450)+B(452)

  JVS(1995) = 0

  JVS(1996) = B(102)+B(110)+2*B(113)+2*B(115)+B(145)+B(281)+B(283)+B(285)+B(287)+B(289)+B(291)+B(293)+B(295)+B(297)&
                &+B(299)+B(301)+B(303)+B(305)+B(307)+B(309)+B(311)+B(313)+B(315)+B(317)+B(319)+B(321)+B(323)+B(325)+B(327)&
                &+B(329)+B(331)+B(333)+B(335)+B(337)+B(339)+B(341)+B(343)+B(345)+B(347)+B(349)+B(351)+B(353)+2*B(355)+2&
                &*B(357)+2*B(359)+B(361)+B(363)+B(365)+B(367)+2*B(369)+B(371)+B(373)+B(375)+B(377)+B(379)+B(381)+B(383)&
                &+B(385)+2*B(387)+B(389)+B(391)+B(393)+B(1017)+B(1019)+B(1041)+B(1066)

  JVS(1997) = -B(141)+B(143)+B(151)+B(155)+B(755)+B(779)+B(781)+B(783)+B(785)+B(787)+B(789)+B(791)+B(793)+B(795)+B(797)&
                &+B(799)+B(801)+B(803)+B(805)+B(807)+B(809)+B(811)+B(813)+B(815)+B(817)+B(819)+B(827)+B(837)+B(891)+B(893)&
                &+B(895)+B(907)+B(936)+B(938)+B(940)+B(948)+B(950)+B(955)+B(993)+B(1035)+B(1050)

  JVS(1998) = -B(4)-B(104)-B(106)-B(112)-B(119)-B(142)-B(147)-B(825)-B(835)-B(927)-B(930)-B(933)-B(942)-B(945)-B(952)&
                &-B(981)-B(1071)-B(1073)-B(1078)

  JVS(1999) = B(6)+2*B(116)-B(120)+B(144)+B(152)+B(443)+B(445)+B(447)+B(449)+B(451)+B(453)+B(455)+B(457)+B(459)+B(461)&
                &+B(463)+B(465)+B(467)+B(469)+B(471)+B(473)+B(475)+B(477)+B(479)+B(481)+B(483)+B(485)+B(487)+B(489)+B(491)&
                &+B(493)+B(495)+B(497)+B(499)+B(501)+B(503)+B(505)+B(507)+B(509)+B(511)+B(513)+B(515)+2*B(517)+2*B(519)+2&
                &*B(521)+B(523)+B(525)+B(527)+2*B(529)+B(531)+B(533)+B(535)+B(537)+B(539)+B(541)+B(543)+B(545)+B(547)+2&
                &*B(549)+B(829)+B(839)+B(1025)+B(1027)+B(1045)

  JVS(2000) = 0

  JVS(2001) = B(146)-B(148)+B(153)

  JVS(2002) = B(121)

  JVS(2003) = -B(822)

  JVS(2004) = -B(832)

  JVS(2005) = -B(1063)

  JVS(2006) = -B(1088)

  JVS(2007) = -B(828)

  JVS(2008) = -B(838)

  JVS(2009) = -B(998)

  JVS(2010) = -B(184)

  JVS(2011) = -B(522)

  JVS(2012) = -B(520)

  JVS(2013) = -B(528)

  JVS(2014) = -B(516)

  JVS(2015) = -B(182)

  JVS(2016) = -B(1024)

  JVS(2017) = -B(448)

  JVS(2018) = -B(534)

  JVS(2019) = -B(770)

  JVS(2020) = -B(220)

  JVS(2021) = 0

  JVS(2022) = -B(186)

  JVS(2023) = -B(462)

  JVS(2024) = -B(210)

  JVS(2025) = -B(768)

  JVS(2026) = -B(518)

  JVS(2027) = -B(548)

  JVS(2028) = -B(482)

  JVS(2029) = B(105)

  JVS(2030) = -B(1044)

  JVS(2031) = -B(202)

  JVS(2032) = -B(458)

  JVS(2033) = -B(460)

  JVS(2034) = -B(502)-B(504)

  JVS(2035) = -B(464)-B(466)

  JVS(2036) = -B(540)

  JVS(2037) = -B(542)

  JVS(2038) = B(156)

  JVS(2039) = -B(526)

  JVS(2040) = -B(524)

  JVS(2041) = -B(538)

  JVS(2042) = -B(492)

  JVS(2043) = -B(468)

  JVS(2044) = -B(746)

  JVS(2045) = -B(486)

  JVS(2046) = -B(730)

  JVS(2047) = -B(510)-B(512)-B(514)

  JVS(2048) = -B(506)-B(508)

  JVS(2049) = -B(472)-B(474)

  JVS(2050) = -B(1026)

  JVS(2051) = -B(722)

  JVS(2052) = -B(470)

  JVS(2053) = -B(490)

  JVS(2054) = -B(536)

  JVS(2055) = -B(456)

  JVS(2056) = -B(252)

  JVS(2057) = -B(480)

  JVS(2058) = -B(496)

  JVS(2059) = -B(444)

  JVS(2060) = -B(546)

  JVS(2061) = -B(442)

  JVS(2062) = -B(532)

  JVS(2063) = -B(478)

  JVS(2064) = -B(544)

  JVS(2065) = -B(446)

  JVS(2066) = -B(498)

  JVS(2067) = -B(476)

  JVS(2068) = -B(500)

  JVS(2069) = -B(488)

  JVS(2070) = -B(494)

  JVS(2071) = -B(530)

  JVS(2072) = B(111)

  JVS(2073) = -B(454)

  JVS(2074) = -B(484)

  JVS(2075) = -B(450)-B(452)

  JVS(2076) = -B(256)

  JVS(2077) = -B(115)

  JVS(2078) = -B(143)+B(157)

  JVS(2079) = B(106)+B(112)-B(117)-B(119)

  JVS(2080) = -B(5)-B(6)-B(116)-B(118)-B(120)-B(144)-B(152)-B(183)-B(185)-B(187)-B(203)-B(211)-B(221)-B(253)-B(254)&
                &-B(257)-B(443)-B(445)-B(447)-B(449)-B(451)-B(453)-B(455)-B(457)-B(459)-B(461)-B(463)-B(465)-B(467)-B(469)&
                &-B(471)-B(473)-B(475)-B(477)-B(479)-B(481)-B(483)-B(485)-B(487)-B(489)-B(491)-B(493)-B(495)-B(497)-B(499)&
                &-B(501)-B(503)-B(505)-B(507)-B(509)-B(511)-B(513)-B(515)-B(517)-B(519)-B(521)-B(523)-B(525)-B(527)-B(529)&
                &-B(531)-B(533)-B(535)-B(537)-B(539)-B(541)-B(543)-B(545)-B(547)-B(549)-B(723)-B(731)-B(747)-B(769)-B(771)&
                &-B(823)-B(829)-B(833)-B(839)-B(999)-B(1025)-B(1027)-B(1045)-B(1064)-B(1089)

  JVS(2081) = -B(255)

  JVS(2082) = -B(153)

  JVS(2083) = B(264)

  JVS(2084) = 2*B(77)+2*B(894)

  JVS(2085) = B(50)+B(844)

  JVS(2086) = 2*B(72)

  JVS(2087) = B(76)+B(892)

  JVS(2088) = B(38)+B(780)

  JVS(2089) = B(937)

  JVS(2090) = B(71)

  JVS(2091) = B(41)

  JVS(2092) = B(53)

  JVS(2093) = B(73)

  JVS(2094) = B(19)

  JVS(2095) = 2*B(358)+2*B(520)+2*B(691)

  JVS(2096) = B(196)+B(198)

  JVS(2097) = B(356)+B(518)+B(690)

  JVS(2098) = B(20)

  JVS(2099) = 0

  JVS(2100) = 0

  JVS(2101) = B(318)+B(480)+B(671)

  JVS(2102) = B(334)+B(496)+B(679)

  JVS(2103) = B(282)+B(444)+B(647)+B(648)

  JVS(2104) = 2*B(316)+2*B(478)+2*B(670)

  JVS(2105) = 0

  JVS(2106) = B(336)+B(498)+B(680)

  JVS(2107) = B(314)+B(476)+B(669)

  JVS(2108) = 0

  JVS(2109) = 0

  JVS(2110) = 0

  JVS(2111) = B(197)+B(199)

  JVS(2112) = 0

  JVS(2113) = 0

  JVS(2114) = B(288)+B(450)+B(656)

  JVS(2115) = 0

  JVS(2116) = B(283)+B(289)+B(315)+2*B(317)+B(319)+B(335)+B(337)+B(357)+2*B(359)

  JVS(2117) = -B(248)+B(265)+B(781)+B(845)+B(893)+2*B(895)+B(938)

  JVS(2118) = 0

  JVS(2119) = -B(254)+B(445)+B(451)+B(477)+2*B(479)+B(481)+B(497)+B(499)+B(519)+2*B(521)

  JVS(2120) = -B(11)-B(249)-B(255)

  JVS(2121) = 0

  JVS(2122) = B(162)

  JVS(2123) = B(236)

  JVS(2124) = B(240)

  JVS(2125) = B(230)

  JVS(2126) = B(262)

  JVS(2127) = B(264)

  JVS(2128) = B(272)

  JVS(2129) = B(276)

  JVS(2130) = B(268)

  JVS(2131) = -B(1081)

  JVS(2132) = B(51)

  JVS(2133) = B(68)

  JVS(2134) = B(81)

  JVS(2135) = B(90)

  JVS(2136) = B(50)

  JVS(2137) = B(149)

  JVS(2138) = B(52)

  JVS(2139) = B(994)

  JVS(2140) = B(72)

  JVS(2141) = B(59)

  JVS(2142) = B(40)

  JVS(2143) = B(69)

  JVS(2144) = B(78)

  JVS(2145) = B(37)

  JVS(2146) = B(79)

  JVS(2147) = B(38)

  JVS(2148) = B(48)

  JVS(2149) = B(39)

  JVS(2150) = B(49)

  JVS(2151) = B(987)

  JVS(2152) = B(70)

  JVS(2153) = B(130)

  JVS(2154) = B(1061)

  JVS(2155) = 0.05*B(1086)

  JVS(2156) = B(64)+B(65)

  JVS(2157) = B(66)

  JVS(2158) = B(71)

  JVS(2159) = B(45)

  JVS(2160) = B(46)

  JVS(2161) = B(87)

  JVS(2162) = B(89)

  JVS(2163) = 0

  JVS(2164) = B(986)

  JVS(2165) = B(74)

  JVS(2166) = B(80)

  JVS(2167) = B(985)

  JVS(2168) = B(42)

  JVS(2169) = B(983)

  JVS(2170) = B(54)

  JVS(2171) = B(73)

  JVS(2172) = B(47)

  JVS(2173) = B(19)+B(714)

  JVS(2174) = B(984)

  JVS(2175) = B(360)+B(522)-B(616)+B(692)

  JVS(2176) = -B(614)

  JVS(2177) = -B(622)

  JVS(2178) = -B(610)

  JVS(2179) = 2*B(24)+B(732)

  JVS(2180) = 2*B(30)

  JVS(2181) = B(988)+B(989)

  JVS(2182) = B(188)

  JVS(2183) = B(32)+2*B(93)+B(958)

  JVS(2184) = B(1016)+B(1024)-B(1028)+B(1032)

  JVS(2185) = B(286)+B(448)-B(556)+B(653)

  JVS(2186) = -B(628)

  JVS(2187) = 0

  JVS(2188) = 0

  JVS(2189) = 0

  JVS(2190) = B(300)+B(462)-B(568)+B(660)

  JVS(2191) = 0

  JVS(2192) = B(35)

  JVS(2193) = -B(612)

  JVS(2194) = B(34)+B(95)

  JVS(2195) = -B(642)

  JVS(2196) = B(20)

  JVS(2197) = B(320)+B(482)-B(584)+B(672)

  JVS(2198) = 0

  JVS(2199) = 0

  JVS(2200) = 0

  JVS(2201) = B(1040)+B(1044)-B(1046)+B(1048)

  JVS(2202) = B(204)

  JVS(2203) = -B(564)

  JVS(2204) = -B(566)

  JVS(2205) = B(340)+B(342)+B(502)+B(504)-B(604)+B(682)+B(683)

  JVS(2206) = B(18)+B(712)

  JVS(2207) = B(21)

  JVS(2208) = B(302)+B(304)+B(464)+B(466)-B(570)+B(661)+B(662)

  JVS(2209) = B(33)+B(94)

  JVS(2210) = -B(634)

  JVS(2211) = B(382)+B(542)-B(636)+B(702)

  JVS(2212) = B(364)+B(526)-B(620)+B(694)

  JVS(2213) = B(362)+B(524)-B(618)+B(693)

  JVS(2214) = B(128)

  JVS(2215) = B(378)+B(538)-B(632)+B(700)

  JVS(2216) = B(330)+B(492)-B(594)+B(677)

  JVS(2217) = B(306)+B(468)-B(572)+B(663)

  JVS(2218) = B(25)

  JVS(2219) = B(29)

  JVS(2220) = B(324)+B(486)-B(588)+B(674)

  JVS(2221) = 2*B(22)

  JVS(2222) = B(350)+B(352)+B(512)+B(514)-B(608)+B(687)+B(688)

  JVS(2223) = B(346)+B(508)-B(606)+B(685)

  JVS(2224) = B(310)+B(312)+B(472)+B(474)-B(576)+B(667)+B(668)

  JVS(2225) = B(1018)+B(1026)-B(1030)+B(1033)

  JVS(2226) = B(23)

  JVS(2227) = B(308)+B(470)-B(574)+B(664)

  JVS(2228) = -B(592)

  JVS(2229) = -B(630)

  JVS(2230) = -B(562)

  JVS(2231) = 2*B(9)+B(246)+B(252)

  JVS(2232) = B(318)+B(480)-B(582)+B(671)

  JVS(2233) = -B(598)

  JVS(2234) = B(282)+B(444)-B(552)+B(647)

  JVS(2235) = -B(640)

  JVS(2236) = B(280)+B(442)-B(550)+B(644)

  JVS(2237) = -B(626)

  JVS(2238) = B(316)+B(478)-B(580)+B(670)

  JVS(2239) = B(388)+B(544)-B(638)+B(703)

  JVS(2240) = B(284)+B(446)-B(554)+B(650)

  JVS(2241) = -B(600)

  JVS(2242) = B(314)+B(476)-B(578)+B(669)

  JVS(2243) = -B(602)

  JVS(2244) = -B(590)

  JVS(2245) = -B(596)

  JVS(2246) = -B(624)

  JVS(2247) = B(124)-B(132)+B(189)+B(205)

  JVS(2248) = -B(560)

  JVS(2249) = B(322)+B(484)-B(586)+B(673)

  JVS(2250) = B(290)+B(452)-B(558)+B(657)

  JVS(2251) = B(12)

  JVS(2252) = -B(145)+B(281)+B(283)+B(285)+B(287)+B(291)+B(301)+B(303)+B(305)+B(307)+B(309)+B(311)+B(313)+B(315)+B(317)&
                &+B(319)+B(321)+B(323)+B(325)+B(331)+B(341)+B(343)+B(347)+B(351)+B(353)+B(361)+B(363)+B(365)+B(379)+B(383)&
                &+B(389)+B(1017)+B(1019)+B(1041)

  JVS(2253) = B(125)+B(126)+B(129)+B(131)-B(134)+B(143)+B(231)+B(237)+B(241)+B(247)+B(263)+B(265)+B(269)+B(273)+B(277)&
                &+B(713)+B(715)+B(733)+B(959)+B(995)+B(1062)+0.05*B(1087)

  JVS(2254) = -B(147)

  JVS(2255) = B(144)-B(152)+B(253)+B(443)+B(445)+B(447)+B(449)+B(453)+B(463)+B(465)+B(467)+B(469)+B(471)+B(473)+B(475)&
                &+B(477)+B(479)+B(481)+B(483)+B(485)+B(487)+B(493)+B(503)+B(505)+B(509)+B(513)+B(515)+B(523)+B(525)+B(527)&
                &+B(539)+B(543)+B(545)+B(1025)+B(1027)+B(1045)

  JVS(2256) = B(11)

  JVS(2257) = -B(133)-B(135)-2*B(136)-2*B(137)-B(146)-B(148)-B(153)-B(551)-B(553)-B(555)-B(557)-B(559)-B(561)-B(563)&
                &-B(565)-B(567)-B(569)-B(571)-B(573)-B(575)-B(577)-B(579)-B(581)-B(583)-B(585)-B(587)-B(589)-B(591)-B(593)&
                &-B(595)-B(597)-B(599)-B(601)-B(603)-B(605)-B(607)-B(609)-B(611)-B(613)-B(615)-B(617)-B(619)-B(621)-B(623)&
                &-B(625)-B(627)-B(629)-B(631)-B(633)-B(635)-B(637)-B(639)-B(641)-B(643)-B(1029)-B(1031)-B(1047)-B(1082)
      
END SUBROUTINE crimech_Jac_SP














SUBROUTINE crimech_KppDecomp( JVS, IER )







      INTEGER  :: IER
      REAL(kind=dp) :: JVS(2257), W(232), a
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
      
END SUBROUTINE crimech_KppDecomp



SUBROUTINE crimech_KppDecompCmplx( JVS, IER )







      INTEGER  :: IER
      DOUBLE COMPLEX :: JVS(2257), W(232), a
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
      
END SUBROUTINE crimech_KppDecompCmplx


SUBROUTINE crimech_KppSolveIndirect( JVS, X )







      INTEGER i, j
      REAL(kind=dp) JVS(2257), X(232), sum

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
      
END SUBROUTINE crimech_KppSolveIndirect


SUBROUTINE crimech_KppSolveCmplx( JVS, X )







      INTEGER i, j
      DOUBLE COMPLEX JVS(2257), X(232), sum

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
      
END SUBROUTINE crimech_KppSolveCmplx













SUBROUTINE crimech_KppSolve ( JVS, X )


  REAL(kind=dp) :: JVS(LU_NONZERO)

  REAL(kind=dp) :: X(NVAR)

  X(45) = X(45)-JVS(109)*X(35)
  X(46) = X(46)-JVS(115)*X(29)
  X(47) = X(47)-JVS(119)*X(46)
  X(48) = X(48)-JVS(124)*X(30)
  X(49) = X(49)-JVS(128)*X(48)
  X(107) = X(107)-JVS(370)*X(72)-JVS(371)*X(97)
  X(109) = X(109)-JVS(383)*X(41)-JVS(384)*X(108)
  X(110) = X(110)-JVS(388)*X(108)
  X(119) = X(119)-JVS(425)*X(47)
  X(121) = X(121)-JVS(434)*X(49)
  X(131) = X(131)-JVS(479)*X(57)-JVS(480)*X(103)
  X(132) = X(132)-JVS(486)*X(96)-JVS(487)*X(124)-JVS(488)*X(125)
  X(135) = X(135)-JVS(507)*X(77)
  X(139) = X(139)-JVS(527)*X(79)-JVS(528)*X(80)
  X(140) = X(140)-JVS(535)*X(55)-JVS(536)*X(82)
  X(142) = X(142)-JVS(547)*X(110)
  X(148) = X(148)-JVS(578)*X(59)-JVS(579)*X(140)
  X(149) = X(149)-JVS(589)*X(134)-JVS(590)*X(136)
  X(150) = X(150)-JVS(598)*X(35)-JVS(599)*X(89)-JVS(600)*X(91)-JVS(601)*X(118)-JVS(602)*X(120)
  X(151) = X(151)-JVS(612)*X(64)-JVS(613)*X(144)
  X(152) = X(152)-JVS(620)*X(104)-JVS(621)*X(128)
  X(154) = X(154)-JVS(632)*X(89)-JVS(633)*X(118)
  X(155) = X(155)-JVS(642)*X(20)-JVS(643)*X(21)-JVS(644)*X(22)-JVS(645)*X(23)-JVS(646)*X(27)-JVS(647)*X(28)
  X(156) = X(156)-JVS(653)*X(25)
  X(157) = X(157)-JVS(659)*X(70)-JVS(660)*X(71)
  X(158) = X(158)-JVS(667)*X(123)
  X(161) = X(161)-JVS(682)*X(88)-JVS(683)*X(151)
  X(163) = X(163)-JVS(695)*X(29)
  X(165) = X(165)-JVS(705)*X(66)-JVS(706)*X(67)-JVS(707)*X(69)-JVS(708)*X(146)-JVS(709)*X(164)
  X(166) = X(166)-JVS(717)*X(143)
  X(167) = X(167)-JVS(724)*X(133)-JVS(725)*X(141)-JVS(726)*X(155)
  X(168) = X(168)-JVS(735)*X(125)-JVS(736)*X(159)
  X(169) = X(169)-JVS(743)*X(102)-JVS(744)*X(119)-JVS(745)*X(137)
  X(170) = X(170)-JVS(754)*X(103)
  X(171) = X(171)-JVS(761)*X(109)
  X(172) = X(172)-JVS(768)*X(142)
  X(173) = X(173)-JVS(775)*X(7)-JVS(776)*X(127)-JVS(777)*X(171)-JVS(778)*X(172)
  X(174) = X(174)-JVS(785)*X(24)
  X(178) = X(178)-JVS(805)*X(175)
  X(179) = X(179)-JVS(812)*X(54)-JVS(813)*X(116)-JVS(814)*X(117)
  X(180) = X(180)-JVS(824)*X(78)-JVS(825)*X(121)-JVS(826)*X(122)-JVS(827)*X(129)-JVS(828)*X(170)
  X(181) = X(181)-JVS(838)*X(30)
  X(182) = X(182)-JVS(844)*X(91)-JVS(845)*X(120)-JVS(846)*X(130)-JVS(847)*X(138)-JVS(848)*X(181)
  X(183) = X(183)-JVS(857)*X(94)-JVS(858)*X(135)-JVS(859)*X(158)-JVS(860)*X(160)
  X(184) = X(184)-JVS(868)*X(124)-JVS(869)*X(159)
  X(185) = X(185)-JVS(876)*X(46)-JVS(877)*X(48)-JVS(878)*X(108)-JVS(879)*X(109)-JVS(880)*X(119)-JVS(881)*X(121)-JVS(882)&
             &*X(135)-JVS(883)*X(158)-JVS(884)*X(160)-JVS(885)*X(165)
  X(186) = X(186)-JVS(902)*X(164)
  X(187) = X(187)-JVS(909)*X(43)-JVS(910)*X(75)-JVS(911)*X(85)-JVS(912)*X(93)-JVS(913)*X(95)-JVS(914)*X(99)-JVS(915)&
             &*X(100)-JVS(916)*X(122)-JVS(917)*X(156)-JVS(918)*X(160)-JVS(919)*X(164)-JVS(920)*X(184)-JVS(921)*X(186)
  X(188) = X(188)-JVS(934)*X(151)-JVS(935)*X(175)
  X(189) = X(189)-JVS(942)*X(35)-JVS(943)*X(88)-JVS(944)*X(94)-JVS(945)*X(111)-JVS(946)*X(112)-JVS(947)*X(116)-JVS(948)&
             &*X(126)-JVS(949)*X(135)-JVS(950)*X(143)-JVS(951)*X(150)-JVS(952)*X(151)-JVS(953)*X(153)-JVS(954)*X(158)&
             &-JVS(955)*X(159)-JVS(956)*X(160)-JVS(957)*X(162)-JVS(958)*X(163)-JVS(959)*X(175)-JVS(960)*X(181)-JVS(961)&
             &*X(183)-JVS(962)*X(188)
  X(190) = X(190)-JVS(986)*X(93)-JVS(987)*X(123)
  X(191) = X(191)-JVS(994)*X(77)-JVS(995)*X(84)-JVS(996)*X(112)-JVS(997)*X(160)
  X(192) = X(192)-JVS(1006)*X(26)
  X(193) = X(193)-JVS(1012)*X(58)-JVS(1013)*X(104)-JVS(1014)*X(114)-JVS(1015)*X(128)-JVS(1016)*X(130)-JVS(1017)*X(133)&
             &-JVS(1018)*X(138)-JVS(1019)*X(141)-JVS(1020)*X(155)-JVS(1021)*X(161)-JVS(1022)*X(174)-JVS(1023)*X(175)&
             &-JVS(1024)*X(179)-JVS(1025)*X(181)-JVS(1026)*X(188)-JVS(1027)*X(192)
  X(194) = X(194)-JVS(1043)*X(63)-JVS(1044)*X(113)-JVS(1045)*X(178)
  X(195) = X(195)-JVS(1052)*X(57)-JVS(1053)*X(164)
  X(196) = X(196)-JVS(1061)*X(65)-JVS(1062)*X(86)-JVS(1063)*X(105)-JVS(1064)*X(114)-JVS(1065)*X(194)
  X(197) = X(197)-JVS(1075)*X(86)-JVS(1076)*X(116)
  X(198) = X(198)-JVS(1084)*X(87)-JVS(1085)*X(114)-JVS(1086)*X(194)
  X(199) = X(199)-JVS(1094)*X(15)-JVS(1095)*X(16)-JVS(1096)*X(17)-JVS(1097)*X(18)-JVS(1098)*X(42)-JVS(1099)*X(76)&
             &-JVS(1100)*X(152)-JVS(1101)*X(153)-JVS(1102)*X(154)-JVS(1103)*X(163)-JVS(1104)*X(167)-JVS(1105)*X(174)&
             &-JVS(1106)*X(181)-JVS(1107)*X(182)-JVS(1108)*X(192)
  X(200) = X(200)-JVS(1118)*X(27)-JVS(1119)*X(28)
  X(201) = X(201)-JVS(1125)*X(87)-JVS(1126)*X(113)-JVS(1127)*X(175)-JVS(1128)*X(178)-JVS(1129)*X(198)
  X(202) = X(202)-JVS(1137)*X(26)
  X(203) = X(203)-JVS(1143)*X(62)-JVS(1144)*X(111)
  X(204) = X(204)-JVS(1153)*X(71)-JVS(1154)*X(93)-JVS(1155)*X(157)
  X(205) = X(205)-JVS(1162)*X(184)-JVS(1163)*X(186)
  X(206) = X(206)-JVS(1170)*X(36)-JVS(1171)*X(60)-JVS(1172)*X(83)-JVS(1173)*X(84)-JVS(1174)*X(90)-JVS(1175)*X(96)&
             &-JVS(1176)*X(101)-JVS(1177)*X(105)-JVS(1178)*X(110)-JVS(1179)*X(112)-JVS(1180)*X(113)-JVS(1181)*X(117)&
             &-JVS(1182)*X(123)-JVS(1183)*X(124)-JVS(1184)*X(125)-JVS(1185)*X(126)-JVS(1186)*X(143)-JVS(1187)*X(147)&
             &-JVS(1188)*X(153)-JVS(1189)*X(159)-JVS(1190)*X(161)-JVS(1191)*X(166)-JVS(1192)*X(168)-JVS(1193)*X(175)&
             &-JVS(1194)*X(178)-JVS(1195)*X(179)-JVS(1196)*X(184)-JVS(1197)*X(188)-JVS(1198)*X(190)-JVS(1199)*X(191)&
             &-JVS(1200)*X(196)-JVS(1201)*X(197)-JVS(1202)*X(198)-JVS(1203)*X(199)-JVS(1204)*X(200)-JVS(1205)*X(201)&
             &-JVS(1206)*X(202)-JVS(1207)*X(204)
  X(207) = X(207)-JVS(1222)*X(37)-JVS(1223)*X(79)-JVS(1224)*X(81)-JVS(1225)*X(177)
  X(208) = X(208)-JVS(1233)*X(31)-JVS(1234)*X(149)-JVS(1235)*X(176)
  X(209) = X(209)-JVS(1243)*X(8)-JVS(1244)*X(14)-JVS(1245)*X(31)-JVS(1246)*X(58)-JVS(1247)*X(62)-JVS(1248)*X(134)&
             &-JVS(1249)*X(136)-JVS(1250)*X(149)-JVS(1251)*X(154)-JVS(1252)*X(163)-JVS(1253)*X(167)-JVS(1254)*X(176)&
             &-JVS(1255)*X(181)-JVS(1256)*X(192)-JVS(1257)*X(200)-JVS(1258)*X(202)-JVS(1259)*X(203)
  X(210) = X(210)-JVS(1271)*X(107)-JVS(1272)*X(159)
  X(211) = X(211)-JVS(1280)*X(11)-JVS(1281)*X(12)-JVS(1282)*X(13)-JVS(1283)*X(19)-JVS(1284)*X(34)-JVS(1285)*X(41)&
             &-JVS(1286)*X(51)-JVS(1287)*X(52)-JVS(1288)*X(53)-JVS(1289)*X(101)-JVS(1290)*X(108)-JVS(1291)*X(143)-JVS(1292)&
             &*X(162)-JVS(1293)*X(171)-JVS(1294)*X(172)-JVS(1295)*X(182)-JVS(1296)*X(187)-JVS(1297)*X(192)-JVS(1298)*X(200)&
             &-JVS(1299)*X(202)-JVS(1300)*X(204)-JVS(1301)*X(210)
  X(212) = X(212)-JVS(1314)*X(39)-JVS(1315)*X(68)-JVS(1316)*X(70)-JVS(1317)*X(164)
  X(213) = X(213)-JVS(1325)*X(38)-JVS(1326)*X(44)-JVS(1327)*X(99)-JVS(1328)*X(100)-JVS(1329)*X(115)-JVS(1330)*X(162)&
             &-JVS(1331)*X(176)-JVS(1332)*X(210)
  X(214) = X(214)-JVS(1341)*X(132)-JVS(1342)*X(159)-JVS(1343)*X(168)-JVS(1344)*X(184)
  X(215) = X(215)-JVS(1351)*X(25)-JVS(1352)*X(59)-JVS(1353)*X(131)-JVS(1354)*X(139)-JVS(1355)*X(152)-JVS(1356)*X(174)&
             &-JVS(1357)*X(177)-JVS(1358)*X(182)-JVS(1359)*X(192)-JVS(1360)*X(200)-JVS(1361)*X(202)-JVS(1362)*X(205)&
             &-JVS(1363)*X(208)-JVS(1364)*X(214)
  X(216) = X(216)-JVS(1374)*X(139)-JVS(1375)*X(177)
  X(217) = X(217)-JVS(1383)*X(43)-JVS(1384)*X(50)-JVS(1385)*X(77)-JVS(1386)*X(143)-JVS(1387)*X(148)-JVS(1388)*X(156)&
             &-JVS(1389)*X(160)-JVS(1390)*X(200)-JVS(1391)*X(202)-JVS(1392)*X(204)-JVS(1393)*X(208)-JVS(1394)*X(212)&
             &-JVS(1395)*X(213)-JVS(1396)*X(214)-JVS(1397)*X(215)-JVS(1398)*X(216)
  X(218) = X(218)-JVS(1409)*X(131)-JVS(1410)*X(205)
  X(219) = X(219)-JVS(1417)*X(31)-JVS(1418)*X(44)-JVS(1419)*X(53)-JVS(1420)*X(86)-JVS(1421)*X(115)-JVS(1422)*X(126)&
             &-JVS(1423)*X(139)-JVS(1424)*X(140)-JVS(1425)*X(148)-JVS(1426)*X(149)-JVS(1427)*X(158)-JVS(1428)*X(161)&
             &-JVS(1429)*X(175)-JVS(1430)*X(176)-JVS(1431)*X(177)-JVS(1432)*X(179)-JVS(1433)*X(187)-JVS(1434)*X(188)&
             &-JVS(1435)*X(190)-JVS(1436)*X(193)-JVS(1437)*X(194)-JVS(1438)*X(197)-JVS(1439)*X(198)-JVS(1440)*X(200)&
             &-JVS(1441)*X(201)-JVS(1442)*X(202)-JVS(1443)*X(204)-JVS(1444)*X(208)-JVS(1445)*X(210)-JVS(1446)*X(212)&
             &-JVS(1447)*X(213)-JVS(1448)*X(214)-JVS(1449)*X(216)-JVS(1450)*X(217)-JVS(1451)*X(218)
  X(220) = X(220)-JVS(1464)*X(37)-JVS(1465)*X(44)-JVS(1466)*X(81)-JVS(1467)*X(94)-JVS(1468)*X(115)-JVS(1469)*X(131)&
             &-JVS(1470)*X(169)-JVS(1471)*X(180)-JVS(1472)*X(183)-JVS(1473)*X(187)-JVS(1474)*X(190)-JVS(1475)*X(195)&
             &-JVS(1476)*X(204)-JVS(1477)*X(205)-JVS(1478)*X(207)-JVS(1479)*X(210)-JVS(1480)*X(212)-JVS(1481)*X(214)&
             &-JVS(1482)*X(216)-JVS(1483)*X(218)
  X(221) = X(221)-JVS(1494)*X(68)-JVS(1495)*X(122)-JVS(1496)*X(165)-JVS(1497)*X(186)
  X(222) = X(222)-JVS(1505)*X(142)-JVS(1506)*X(143)-JVS(1507)*X(153)-JVS(1508)*X(159)-JVS(1509)*X(162)-JVS(1510)*X(164)&
             &-JVS(1511)*X(171)-JVS(1512)*X(172)-JVS(1513)*X(173)-JVS(1514)*X(175)-JVS(1515)*X(194)-JVS(1516)*X(201)
  X(223) = X(223)-JVS(1523)*X(81)-JVS(1524)*X(152)-JVS(1525)*X(174)-JVS(1526)*X(180)-JVS(1527)*X(195)-JVS(1528)*X(205)&
             &-JVS(1529)*X(218)-JVS(1530)*X(221)-JVS(1531)*X(222)
  X(224) = X(224)-JVS(1538)*X(33)-JVS(1539)*X(80)-JVS(1540)*X(95)-JVS(1541)*X(177)-JVS(1542)*X(212)-JVS(1543)*X(218)&
             &-JVS(1544)*X(221)-JVS(1545)*X(222)-JVS(1546)*X(223)
  X(225) = X(225)-JVS(1553)*X(10)-JVS(1554)*X(100)-JVS(1555)*X(115)-JVS(1556)*X(167)-JVS(1557)*X(169)-JVS(1558)*X(192)&
             &-JVS(1559)*X(200)-JVS(1560)*X(202)-JVS(1561)*X(207)-JVS(1562)*X(210)-JVS(1563)*X(214)-JVS(1564)*X(216)&
             &-JVS(1565)*X(218)-JVS(1566)*X(222)-JVS(1567)*X(223)-JVS(1568)*X(224)
  X(226) = X(226)-JVS(1575)*X(50)-JVS(1576)*X(56)-JVS(1577)*X(98)-JVS(1578)*X(129)-JVS(1579)*X(137)-JVS(1580)*X(170)&
             &-JVS(1581)*X(180)-JVS(1582)*X(195)-JVS(1583)*X(205)-JVS(1584)*X(207)-JVS(1585)*X(215)-JVS(1586)*X(218)&
             &-JVS(1587)*X(221)-JVS(1588)*X(222)-JVS(1589)*X(223)-JVS(1590)*X(224)-JVS(1591)*X(225)
  X(227) = X(227)-JVS(1598)*X(40)-JVS(1599)*X(110)-JVS(1600)*X(142)-JVS(1601)*X(144)-JVS(1602)*X(145)-JVS(1603)*X(146)&
             &-JVS(1604)*X(147)-JVS(1605)*X(153)-JVS(1606)*X(155)-JVS(1607)*X(156)-JVS(1608)*X(157)-JVS(1609)*X(162)&
             &-JVS(1610)*X(163)-JVS(1611)*X(164)-JVS(1612)*X(166)-JVS(1613)*X(168)-JVS(1614)*X(170)-JVS(1615)*X(171)&
             &-JVS(1616)*X(172)-JVS(1617)*X(173)-JVS(1618)*X(174)-JVS(1619)*X(175)-JVS(1620)*X(176)-JVS(1621)*X(177)&
             &-JVS(1622)*X(178)-JVS(1623)*X(181)-JVS(1624)*X(183)-JVS(1625)*X(184)-JVS(1626)*X(186)-JVS(1627)*X(188)&
             &-JVS(1628)*X(190)-JVS(1629)*X(191)-JVS(1630)*X(192)-JVS(1631)*X(195)-JVS(1632)*X(196)-JVS(1633)*X(197)&
             &-JVS(1634)*X(198)-JVS(1635)*X(199)-JVS(1636)*X(200)-JVS(1637)*X(201)-JVS(1638)*X(202)-JVS(1639)*X(203)&
             &-JVS(1640)*X(204)-JVS(1641)*X(205)-JVS(1642)*X(207)-JVS(1643)*X(208)-JVS(1644)*X(209)-JVS(1645)*X(210)&
             &-JVS(1646)*X(211)-JVS(1647)*X(212)-JVS(1648)*X(213)-JVS(1649)*X(214)-JVS(1650)*X(215)-JVS(1651)*X(216)&
             &-JVS(1652)*X(217)-JVS(1653)*X(218)-JVS(1654)*X(219)-JVS(1655)*X(220)-JVS(1656)*X(221)-JVS(1657)*X(222)&
             &-JVS(1658)*X(223)-JVS(1659)*X(224)-JVS(1660)*X(225)-JVS(1661)*X(226)
  X(228) = X(228)-JVS(1668)*X(7)-JVS(1669)*X(8)-JVS(1670)*X(10)-JVS(1671)*X(11)-JVS(1672)*X(12)-JVS(1673)*X(13)&
             &-JVS(1674)*X(14)-JVS(1675)*X(15)-JVS(1676)*X(16)-JVS(1677)*X(17)-JVS(1678)*X(18)-JVS(1679)*X(19)-JVS(1680)&
             &*X(20)-JVS(1681)*X(21)-JVS(1682)*X(22)-JVS(1683)*X(23)-JVS(1684)*X(24)-JVS(1685)*X(25)-JVS(1686)*X(26)&
             &-JVS(1687)*X(27)-JVS(1688)*X(28)-JVS(1689)*X(29)-JVS(1690)*X(30)-JVS(1691)*X(31)-JVS(1692)*X(33)-JVS(1693)&
             &*X(35)-JVS(1694)*X(36)-JVS(1695)*X(37)-JVS(1696)*X(38)-JVS(1697)*X(39)-JVS(1698)*X(40)-JVS(1699)*X(41)&
             &-JVS(1700)*X(42)-JVS(1701)*X(43)-JVS(1702)*X(44)-JVS(1703)*X(45)-JVS(1704)*X(46)-JVS(1705)*X(48)-JVS(1706)&
             &*X(50)-JVS(1707)*X(51)-JVS(1708)*X(53)-JVS(1709)*X(54)-JVS(1710)*X(55)-JVS(1711)*X(56)-JVS(1712)*X(57)&
             &-JVS(1713)*X(58)-JVS(1714)*X(59)-JVS(1715)*X(60)-JVS(1716)*X(61)-JVS(1717)*X(62)-JVS(1718)*X(63)-JVS(1719)&
             &*X(64)-JVS(1720)*X(65)-JVS(1721)*X(66)-JVS(1722)*X(67)-JVS(1723)*X(68)-JVS(1724)*X(69)-JVS(1725)*X(70)&
             &-JVS(1726)*X(71)-JVS(1727)*X(72)-JVS(1728)*X(73)-JVS(1729)*X(74)-JVS(1730)*X(75)-JVS(1731)*X(76)-JVS(1732)&
             &*X(77)-JVS(1733)*X(78)-JVS(1734)*X(79)-JVS(1735)*X(80)-JVS(1736)*X(81)-JVS(1737)*X(82)-JVS(1738)*X(83)&
             &-JVS(1739)*X(84)-JVS(1740)*X(85)-JVS(1741)*X(86)-JVS(1742)*X(87)-JVS(1743)*X(88)-JVS(1744)*X(89)-JVS(1745)&
             &*X(90)-JVS(1746)*X(91)-JVS(1747)*X(92)-JVS(1748)*X(93)-JVS(1749)*X(94)-JVS(1750)*X(95)-JVS(1751)*X(96)&
             &-JVS(1752)*X(97)-JVS(1753)*X(98)-JVS(1754)*X(99)-JVS(1755)*X(100)-JVS(1756)*X(101)-JVS(1757)*X(102)-JVS(1758)&
             &*X(103)-JVS(1759)*X(104)-JVS(1760)*X(105)-JVS(1761)*X(106)-JVS(1762)*X(107)-JVS(1763)*X(108)-JVS(1764)*X(109)&
             &-JVS(1765)*X(111)-JVS(1766)*X(112)-JVS(1767)*X(113)-JVS(1768)*X(114)-JVS(1769)*X(115)-JVS(1770)*X(116)&
             &-JVS(1771)*X(117)-JVS(1772)*X(118)-JVS(1773)*X(119)-JVS(1774)*X(120)-JVS(1775)*X(121)-JVS(1776)*X(122)&
             &-JVS(1777)*X(123)-JVS(1778)*X(124)-JVS(1779)*X(125)-JVS(1780)*X(126)-JVS(1781)*X(127)-JVS(1782)*X(128)&
             &-JVS(1783)*X(129)-JVS(1784)*X(130)-JVS(1785)*X(131)-JVS(1786)*X(132)-JVS(1787)*X(133)-JVS(1788)*X(134)&
             &-JVS(1789)*X(135)-JVS(1790)*X(136)-JVS(1791)*X(137)-JVS(1792)*X(138)-JVS(1793)*X(139)-JVS(1794)*X(140)&
             &-JVS(1795)*X(141)-JVS(1796)*X(143)-JVS(1797)*X(144)-JVS(1798)*X(145)-JVS(1799)*X(146)-JVS(1800)*X(147)&
             &-JVS(1801)*X(148)-JVS(1802)*X(149)-JVS(1803)*X(150)-JVS(1804)*X(151)-JVS(1805)*X(152)-JVS(1806)*X(153)&
             &-JVS(1807)*X(154)-JVS(1808)*X(155)-JVS(1809)*X(156)-JVS(1810)*X(157)-JVS(1811)*X(158)-JVS(1812)*X(159)&
             &-JVS(1813)*X(160)-JVS(1814)*X(161)-JVS(1815)*X(162)-JVS(1816)*X(163)-JVS(1817)*X(164)-JVS(1818)*X(165)&
             &-JVS(1819)*X(166)-JVS(1820)*X(167)-JVS(1821)*X(168)-JVS(1822)*X(169)-JVS(1823)*X(170)-JVS(1824)*X(171)&
             &-JVS(1825)*X(172)-JVS(1826)*X(173)-JVS(1827)*X(174)-JVS(1828)*X(175)-JVS(1829)*X(176)-JVS(1830)*X(177)&
             &-JVS(1831)*X(178)-JVS(1832)*X(179)-JVS(1833)*X(180)-JVS(1834)*X(181)-JVS(1835)*X(182)-JVS(1836)*X(183)&
             &-JVS(1837)*X(184)-JVS(1838)*X(185)-JVS(1839)*X(186)-JVS(1840)*X(187)-JVS(1841)*X(188)-JVS(1842)*X(189)&
             &-JVS(1843)*X(190)-JVS(1844)*X(191)-JVS(1845)*X(192)-JVS(1846)*X(193)-JVS(1847)*X(194)-JVS(1848)*X(195)&
             &-JVS(1849)*X(196)-JVS(1850)*X(197)-JVS(1851)*X(198)-JVS(1852)*X(199)-JVS(1853)*X(200)-JVS(1854)*X(201)&
             &-JVS(1855)*X(202)-JVS(1856)*X(203)-JVS(1857)*X(204)-JVS(1858)*X(205)-JVS(1859)*X(206)-JVS(1860)*X(207)&
             &-JVS(1861)*X(208)-JVS(1862)*X(209)-JVS(1863)*X(210)-JVS(1864)*X(211)-JVS(1865)*X(212)-JVS(1866)*X(213)&
             &-JVS(1867)*X(214)-JVS(1868)*X(215)-JVS(1869)*X(216)-JVS(1870)*X(217)-JVS(1871)*X(218)-JVS(1872)*X(219)&
             &-JVS(1873)*X(220)-JVS(1874)*X(221)-JVS(1875)*X(222)-JVS(1876)*X(223)-JVS(1877)*X(224)-JVS(1878)*X(225)&
             &-JVS(1879)*X(226)-JVS(1880)*X(227)
  X(229) = X(229)-JVS(1886)*X(32)-JVS(1887)*X(34)-JVS(1888)*X(40)-JVS(1889)*X(47)-JVS(1890)*X(49)-JVS(1891)*X(54)&
             &-JVS(1892)*X(55)-JVS(1893)*X(60)-JVS(1894)*X(61)-JVS(1895)*X(63)-JVS(1896)*X(65)-JVS(1897)*X(66)-JVS(1898)&
             &*X(69)-JVS(1899)*X(74)-JVS(1900)*X(78)-JVS(1901)*X(79)-JVS(1902)*X(83)-JVS(1903)*X(85)-JVS(1904)*X(87)&
             &-JVS(1905)*X(90)-JVS(1906)*X(92)-JVS(1907)*X(95)-JVS(1908)*X(96)-JVS(1909)*X(97)-JVS(1910)*X(98)-JVS(1911)&
             &*X(99)-JVS(1912)*X(102)-JVS(1913)*X(103)-JVS(1914)*X(110)-JVS(1915)*X(111)-JVS(1916)*X(112)-JVS(1917)*X(116)&
             &-JVS(1918)*X(118)-JVS(1919)*X(119)-JVS(1920)*X(120)-JVS(1921)*X(121)-JVS(1922)*X(122)-JVS(1923)*X(125)&
             &-JVS(1924)*X(126)-JVS(1925)*X(128)-JVS(1926)*X(134)-JVS(1927)*X(138)-JVS(1928)*X(141)-JVS(1929)*X(142)&
             &-JVS(1930)*X(144)-JVS(1931)*X(145)-JVS(1932)*X(146)-JVS(1933)*X(147)-JVS(1934)*X(153)-JVS(1935)*X(155)&
             &-JVS(1936)*X(156)-JVS(1937)*X(157)-JVS(1938)*X(160)-JVS(1939)*X(161)-JVS(1940)*X(162)-JVS(1941)*X(163)&
             &-JVS(1942)*X(164)-JVS(1943)*X(166)-JVS(1944)*X(168)-JVS(1945)*X(170)-JVS(1946)*X(171)-JVS(1947)*X(172)&
             &-JVS(1948)*X(173)-JVS(1949)*X(174)-JVS(1950)*X(175)-JVS(1951)*X(176)-JVS(1952)*X(177)-JVS(1953)*X(178)&
             &-JVS(1954)*X(181)-JVS(1955)*X(183)-JVS(1956)*X(184)-JVS(1957)*X(185)-JVS(1958)*X(186)-JVS(1959)*X(188)&
             &-JVS(1960)*X(190)-JVS(1961)*X(191)-JVS(1962)*X(192)-JVS(1963)*X(194)-JVS(1964)*X(195)-JVS(1965)*X(196)&
             &-JVS(1966)*X(197)-JVS(1967)*X(198)-JVS(1968)*X(199)-JVS(1969)*X(200)-JVS(1970)*X(201)-JVS(1971)*X(202)&
             &-JVS(1972)*X(203)-JVS(1973)*X(204)-JVS(1974)*X(205)-JVS(1975)*X(206)-JVS(1976)*X(207)-JVS(1977)*X(208)&
             &-JVS(1978)*X(209)-JVS(1979)*X(210)-JVS(1980)*X(211)-JVS(1981)*X(212)-JVS(1982)*X(213)-JVS(1983)*X(214)&
             &-JVS(1984)*X(215)-JVS(1985)*X(216)-JVS(1986)*X(217)-JVS(1987)*X(218)-JVS(1988)*X(219)-JVS(1989)*X(220)&
             &-JVS(1990)*X(221)-JVS(1991)*X(222)-JVS(1992)*X(223)-JVS(1993)*X(224)-JVS(1994)*X(225)-JVS(1995)*X(226)&
             &-JVS(1996)*X(227)-JVS(1997)*X(228)
  X(230) = X(230)-JVS(2002)*X(32)-JVS(2003)*X(46)-JVS(2004)*X(48)-JVS(2005)*X(108)-JVS(2006)*X(109)-JVS(2007)*X(119)&
             &-JVS(2008)*X(121)-JVS(2009)*X(135)-JVS(2010)*X(143)-JVS(2011)*X(144)-JVS(2012)*X(145)-JVS(2013)*X(146)&
             &-JVS(2014)*X(147)-JVS(2015)*X(153)-JVS(2016)*X(155)-JVS(2017)*X(156)-JVS(2018)*X(157)-JVS(2019)*X(158)&
             &-JVS(2020)*X(159)-JVS(2021)*X(160)-JVS(2022)*X(162)-JVS(2023)*X(163)-JVS(2024)*X(164)-JVS(2025)*X(165)&
             &-JVS(2026)*X(166)-JVS(2027)*X(168)-JVS(2028)*X(170)-JVS(2029)*X(173)-JVS(2030)*X(174)-JVS(2031)*X(175)&
             &-JVS(2032)*X(176)-JVS(2033)*X(177)-JVS(2034)*X(178)-JVS(2035)*X(181)-JVS(2036)*X(183)-JVS(2037)*X(184)&
             &-JVS(2038)*X(185)-JVS(2039)*X(186)-JVS(2040)*X(188)-JVS(2041)*X(190)-JVS(2042)*X(191)-JVS(2043)*X(192)&
             &-JVS(2044)*X(194)-JVS(2045)*X(195)-JVS(2046)*X(196)-JVS(2047)*X(197)-JVS(2048)*X(198)-JVS(2049)*X(199)&
             &-JVS(2050)*X(200)-JVS(2051)*X(201)-JVS(2052)*X(202)-JVS(2053)*X(203)-JVS(2054)*X(204)-JVS(2055)*X(205)&
             &-JVS(2056)*X(206)-JVS(2057)*X(207)-JVS(2058)*X(208)-JVS(2059)*X(209)-JVS(2060)*X(210)-JVS(2061)*X(211)&
             &-JVS(2062)*X(212)-JVS(2063)*X(213)-JVS(2064)*X(214)-JVS(2065)*X(215)-JVS(2066)*X(216)-JVS(2067)*X(217)&
             &-JVS(2068)*X(218)-JVS(2069)*X(219)-JVS(2070)*X(220)-JVS(2071)*X(221)-JVS(2072)*X(222)-JVS(2073)*X(223)&
             &-JVS(2074)*X(224)-JVS(2075)*X(225)-JVS(2076)*X(226)-JVS(2077)*X(227)-JVS(2078)*X(228)-JVS(2079)*X(229)
  X(231) = X(231)-JVS(2083)*X(42)-JVS(2084)*X(61)-JVS(2085)*X(73)-JVS(2086)*X(82)-JVS(2087)*X(83)-JVS(2088)*X(92)&
             &-JVS(2089)*X(111)-JVS(2090)*X(117)-JVS(2091)*X(134)-JVS(2092)*X(136)-JVS(2093)*X(137)-JVS(2094)*X(140)&
             &-JVS(2095)*X(145)-JVS(2096)*X(162)-JVS(2097)*X(166)-JVS(2098)*X(169)-JVS(2099)*X(176)-JVS(2100)*X(203)&
             &-JVS(2101)*X(207)-JVS(2102)*X(208)-JVS(2103)*X(209)-JVS(2104)*X(213)-JVS(2105)*X(214)-JVS(2106)*X(216)&
             &-JVS(2107)*X(217)-JVS(2108)*X(218)-JVS(2109)*X(220)-JVS(2110)*X(221)-JVS(2111)*X(222)-JVS(2112)*X(223)&
             &-JVS(2113)*X(224)-JVS(2114)*X(225)-JVS(2115)*X(226)-JVS(2116)*X(227)-JVS(2117)*X(228)-JVS(2118)*X(229)&
             &-JVS(2119)*X(230)
  X(232) = X(232)-JVS(2122)*X(9)-JVS(2123)*X(29)-JVS(2124)*X(30)-JVS(2125)*X(35)-JVS(2126)*X(36)-JVS(2127)*X(42)&
             &-JVS(2128)*X(43)-JVS(2129)*X(45)-JVS(2130)*X(50)-JVS(2131)*X(52)-JVS(2132)*X(56)-JVS(2133)*X(64)-JVS(2134)&
             &*X(67)-JVS(2135)*X(72)-JVS(2136)*X(73)-JVS(2137)*X(74)-JVS(2138)*X(75)-JVS(2139)*X(77)-JVS(2140)*X(82)&
             &-JVS(2141)*X(84)-JVS(2142)*X(85)-JVS(2143)*X(88)-JVS(2144)*X(89)-JVS(2145)*X(90)-JVS(2146)*X(91)-JVS(2147)&
             &*X(92)-JVS(2148)*X(97)-JVS(2149)*X(98)-JVS(2150)*X(101)-JVS(2151)*X(104)-JVS(2152)*X(105)-JVS(2153)*X(106)&
             &-JVS(2154)*X(108)-JVS(2155)*X(109)-JVS(2156)*X(113)-JVS(2157)*X(114)-JVS(2158)*X(117)-JVS(2159)*X(118)&
             &-JVS(2160)*X(120)-JVS(2161)*X(123)-JVS(2162)*X(124)-JVS(2163)*X(127)-JVS(2164)*X(128)-JVS(2165)*X(129)&
             &-JVS(2166)*X(130)-JVS(2167)*X(133)-JVS(2168)*X(134)-JVS(2169)*X(135)-JVS(2170)*X(136)-JVS(2171)*X(137)&
             &-JVS(2172)*X(138)-JVS(2173)*X(140)-JVS(2174)*X(141)-JVS(2175)*X(144)-JVS(2176)*X(145)-JVS(2177)*X(146)&
             &-JVS(2178)*X(147)-JVS(2179)*X(150)-JVS(2180)*X(151)-JVS(2181)*X(152)-JVS(2182)*X(153)-JVS(2183)*X(154)&
             &-JVS(2184)*X(155)-JVS(2185)*X(156)-JVS(2186)*X(157)-JVS(2187)*X(159)-JVS(2188)*X(160)-JVS(2189)*X(162)&
             &-JVS(2190)*X(163)-JVS(2191)*X(164)-JVS(2192)*X(165)-JVS(2193)*X(166)-JVS(2194)*X(167)-JVS(2195)*X(168)&
             &-JVS(2196)*X(169)-JVS(2197)*X(170)-JVS(2198)*X(171)-JVS(2199)*X(172)-JVS(2200)*X(173)-JVS(2201)*X(174)&
             &-JVS(2202)*X(175)-JVS(2203)*X(176)-JVS(2204)*X(177)-JVS(2205)*X(178)-JVS(2206)*X(179)-JVS(2207)*X(180)&
             &-JVS(2208)*X(181)-JVS(2209)*X(182)-JVS(2210)*X(183)-JVS(2211)*X(184)-JVS(2212)*X(186)-JVS(2213)*X(188)&
             &-JVS(2214)*X(189)-JVS(2215)*X(190)-JVS(2216)*X(191)-JVS(2217)*X(192)-JVS(2218)*X(193)-JVS(2219)*X(194)&
             &-JVS(2220)*X(195)-JVS(2221)*X(196)-JVS(2222)*X(197)-JVS(2223)*X(198)-JVS(2224)*X(199)-JVS(2225)*X(200)&
             &-JVS(2226)*X(201)-JVS(2227)*X(202)-JVS(2228)*X(203)-JVS(2229)*X(204)-JVS(2230)*X(205)-JVS(2231)*X(206)&
             &-JVS(2232)*X(207)-JVS(2233)*X(208)-JVS(2234)*X(209)-JVS(2235)*X(210)-JVS(2236)*X(211)-JVS(2237)*X(212)&
             &-JVS(2238)*X(213)-JVS(2239)*X(214)-JVS(2240)*X(215)-JVS(2241)*X(216)-JVS(2242)*X(217)-JVS(2243)*X(218)&
             &-JVS(2244)*X(219)-JVS(2245)*X(220)-JVS(2246)*X(221)-JVS(2247)*X(222)-JVS(2248)*X(223)-JVS(2249)*X(224)&
             &-JVS(2250)*X(225)-JVS(2251)*X(226)-JVS(2252)*X(227)-JVS(2253)*X(228)-JVS(2254)*X(229)-JVS(2255)*X(230)&
             &-JVS(2256)*X(231)
  X(232) = X(232)/JVS(2257)
  X(231) = (X(231)-JVS(2121)*X(232))/(JVS(2120))
  X(230) = (X(230)-JVS(2081)*X(231)-JVS(2082)*X(232))/(JVS(2080))
  X(229) = (X(229)-JVS(1999)*X(230)-JVS(2000)*X(231)-JVS(2001)*X(232))/(JVS(1998))
  X(228) = (X(228)-JVS(1882)*X(229)-JVS(1883)*X(230)-JVS(1884)*X(231)-JVS(1885)*X(232))/(JVS(1881))
  X(227) = (X(227)-JVS(1663)*X(228)-JVS(1664)*X(229)-JVS(1665)*X(230)-JVS(1666)*X(231)-JVS(1667)*X(232))/(JVS(1662))
  X(226) = (X(226)-JVS(1593)*X(227)-JVS(1594)*X(228)-JVS(1595)*X(229)-JVS(1596)*X(230)-JVS(1597)*X(232))/(JVS(1592))
  X(225) = (X(225)-JVS(1570)*X(227)-JVS(1571)*X(228)-JVS(1572)*X(229)-JVS(1573)*X(230)-JVS(1574)*X(232))/(JVS(1569))
  X(224) = (X(224)-JVS(1548)*X(227)-JVS(1549)*X(228)-JVS(1550)*X(229)-JVS(1551)*X(230)-JVS(1552)*X(232))/(JVS(1547))
  X(223) = (X(223)-JVS(1533)*X(227)-JVS(1534)*X(228)-JVS(1535)*X(229)-JVS(1536)*X(230)-JVS(1537)*X(232))/(JVS(1532))
  X(222) = (X(222)-JVS(1518)*X(227)-JVS(1519)*X(228)-JVS(1520)*X(229)-JVS(1521)*X(230)-JVS(1522)*X(232))/(JVS(1517))
  X(221) = (X(221)-JVS(1499)*X(222)-JVS(1500)*X(227)-JVS(1501)*X(228)-JVS(1502)*X(229)-JVS(1503)*X(230)-JVS(1504)&
             &*X(232))/(JVS(1498))
  X(220) = (X(220)-JVS(1485)*X(221)-JVS(1486)*X(222)-JVS(1487)*X(223)-JVS(1488)*X(224)-JVS(1489)*X(227)-JVS(1490)*X(228)&
             &-JVS(1491)*X(229)-JVS(1492)*X(230)-JVS(1493)*X(232))/(JVS(1484))
  X(219) = (X(219)-JVS(1453)*X(220)-JVS(1454)*X(221)-JVS(1455)*X(222)-JVS(1456)*X(223)-JVS(1457)*X(225)-JVS(1458)*X(227)&
             &-JVS(1459)*X(228)-JVS(1460)*X(229)-JVS(1461)*X(230)-JVS(1462)*X(231)-JVS(1463)*X(232))/(JVS(1452))
  X(218) = (X(218)-JVS(1412)*X(222)-JVS(1413)*X(227)-JVS(1414)*X(228)-JVS(1415)*X(230)-JVS(1416)*X(232))/(JVS(1411))
  X(217) = (X(217)-JVS(1400)*X(218)-JVS(1401)*X(221)-JVS(1402)*X(222)-JVS(1403)*X(223)-JVS(1404)*X(225)-JVS(1405)*X(227)&
             &-JVS(1406)*X(228)-JVS(1407)*X(230)-JVS(1408)*X(232))/(JVS(1399))
  X(216) = (X(216)-JVS(1377)*X(218)-JVS(1378)*X(223)-JVS(1379)*X(227)-JVS(1380)*X(228)-JVS(1381)*X(230)-JVS(1382)&
             &*X(232))/(JVS(1376))
  X(215) = (X(215)-JVS(1366)*X(218)-JVS(1367)*X(222)-JVS(1368)*X(223)-JVS(1369)*X(225)-JVS(1370)*X(227)-JVS(1371)*X(228)&
             &-JVS(1372)*X(230)-JVS(1373)*X(232))/(JVS(1365))
  X(214) = (X(214)-JVS(1346)*X(222)-JVS(1347)*X(227)-JVS(1348)*X(228)-JVS(1349)*X(230)-JVS(1350)*X(232))/(JVS(1345))
  X(213) = (X(213)-JVS(1334)*X(214)-JVS(1335)*X(216)-JVS(1336)*X(222)-JVS(1337)*X(227)-JVS(1338)*X(228)-JVS(1339)*X(230)&
             &-JVS(1340)*X(232))/(JVS(1333))
  X(212) = (X(212)-JVS(1319)*X(221)-JVS(1320)*X(222)-JVS(1321)*X(227)-JVS(1322)*X(228)-JVS(1323)*X(230)-JVS(1324)&
             &*X(232))/(JVS(1318))
  X(211) = (X(211)-JVS(1303)*X(212)-JVS(1304)*X(214)-JVS(1305)*X(219)-JVS(1306)*X(221)-JVS(1307)*X(222)-JVS(1308)*X(227)&
             &-JVS(1309)*X(228)-JVS(1310)*X(229)-JVS(1311)*X(230)-JVS(1312)*X(231)-JVS(1313)*X(232))/(JVS(1302))
  X(210) = (X(210)-JVS(1274)*X(214)-JVS(1275)*X(222)-JVS(1276)*X(227)-JVS(1277)*X(228)-JVS(1278)*X(230)-JVS(1279)&
             &*X(232))/(JVS(1273))
  X(209) = (X(209)-JVS(1261)*X(214)-JVS(1262)*X(216)-JVS(1263)*X(220)-JVS(1264)*X(225)-JVS(1265)*X(226)-JVS(1266)*X(227)&
             &-JVS(1267)*X(228)-JVS(1268)*X(229)-JVS(1269)*X(230)-JVS(1270)*X(232))/(JVS(1260))
  X(208) = (X(208)-JVS(1237)*X(214)-JVS(1238)*X(225)-JVS(1239)*X(227)-JVS(1240)*X(228)-JVS(1241)*X(230)-JVS(1242)&
             &*X(232))/(JVS(1236))
  X(207) = (X(207)-JVS(1227)*X(218)-JVS(1228)*X(223)-JVS(1229)*X(227)-JVS(1230)*X(228)-JVS(1231)*X(230)-JVS(1232)&
             &*X(232))/(JVS(1226))
  X(206) = (X(206)-JVS(1209)*X(209)-JVS(1210)*X(211)-JVS(1211)*X(212)-JVS(1212)*X(214)-JVS(1213)*X(217)-JVS(1214)*X(219)&
             &-JVS(1215)*X(220)-JVS(1216)*X(222)-JVS(1217)*X(227)-JVS(1218)*X(228)-JVS(1219)*X(229)-JVS(1220)*X(230)&
             &-JVS(1221)*X(232))/(JVS(1208))
  X(205) = (X(205)-JVS(1165)*X(222)-JVS(1166)*X(227)-JVS(1167)*X(228)-JVS(1168)*X(230)-JVS(1169)*X(232))/(JVS(1164))
  X(204) = (X(204)-JVS(1157)*X(212)-JVS(1158)*X(227)-JVS(1159)*X(228)-JVS(1160)*X(230)-JVS(1161)*X(232))/(JVS(1156))
  X(203) = (X(203)-JVS(1146)*X(216)-JVS(1147)*X(226)-JVS(1148)*X(227)-JVS(1149)*X(228)-JVS(1150)*X(229)-JVS(1151)*X(230)&
             &-JVS(1152)*X(232))/(JVS(1145))
  X(202) = (X(202)-JVS(1139)*X(227)-JVS(1140)*X(228)-JVS(1141)*X(230)-JVS(1142)*X(232))/(JVS(1138))
  X(201) = (X(201)-JVS(1131)*X(222)-JVS(1132)*X(227)-JVS(1133)*X(228)-JVS(1134)*X(229)-JVS(1135)*X(230)-JVS(1136)&
             &*X(232))/(JVS(1130))
  X(200) = (X(200)-JVS(1121)*X(227)-JVS(1122)*X(228)-JVS(1123)*X(230)-JVS(1124)*X(232))/(JVS(1120))
  X(199) = (X(199)-JVS(1110)*X(200)-JVS(1111)*X(202)-JVS(1112)*X(209)-JVS(1113)*X(222)-JVS(1114)*X(227)-JVS(1115)*X(228)&
             &-JVS(1116)*X(230)-JVS(1117)*X(232))/(JVS(1109))
  X(198) = (X(198)-JVS(1088)*X(222)-JVS(1089)*X(227)-JVS(1090)*X(228)-JVS(1091)*X(229)-JVS(1092)*X(230)-JVS(1093)&
             &*X(232))/(JVS(1087))
  X(197) = (X(197)-JVS(1078)*X(201)-JVS(1079)*X(227)-JVS(1080)*X(228)-JVS(1081)*X(229)-JVS(1082)*X(230)-JVS(1083)&
             &*X(232))/(JVS(1077))
  X(196) = (X(196)-JVS(1067)*X(197)-JVS(1068)*X(198)-JVS(1069)*X(199)-JVS(1070)*X(222)-JVS(1071)*X(227)-JVS(1072)*X(228)&
             &-JVS(1073)*X(230)-JVS(1074)*X(232))/(JVS(1066))
  X(195) = (X(195)-JVS(1055)*X(205)-JVS(1056)*X(222)-JVS(1057)*X(227)-JVS(1058)*X(228)-JVS(1059)*X(230)-JVS(1060)&
             &*X(232))/(JVS(1054))
  X(194) = (X(194)-JVS(1047)*X(222)-JVS(1048)*X(227)-JVS(1049)*X(228)-JVS(1050)*X(230)-JVS(1051)*X(232))/(JVS(1046))
  X(193) = (X(193)-JVS(1029)*X(194)-JVS(1030)*X(197)-JVS(1031)*X(198)-JVS(1032)*X(200)-JVS(1033)*X(201)-JVS(1034)*X(202)&
             &-JVS(1035)*X(217)-JVS(1036)*X(220)-JVS(1037)*X(222)-JVS(1038)*X(227)-JVS(1039)*X(228)-JVS(1040)*X(229)&
             &-JVS(1041)*X(230)-JVS(1042)*X(232))/(JVS(1028))
  X(192) = (X(192)-JVS(1008)*X(227)-JVS(1009)*X(228)-JVS(1010)*X(230)-JVS(1011)*X(232))/(JVS(1007))
  X(191) = (X(191)-JVS(999)*X(196)-JVS(1000)*X(204)-JVS(1001)*X(227)-JVS(1002)*X(228)-JVS(1003)*X(229)-JVS(1004)*X(230)&
             &-JVS(1005)*X(232))/(JVS(998))
  X(190) = (X(190)-JVS(989)*X(204)-JVS(990)*X(227)-JVS(991)*X(228)-JVS(992)*X(230)-JVS(993)*X(232))/(JVS(988))
  X(189) = (X(189)-JVS(964)*X(190)-JVS(965)*X(191)-JVS(966)*X(192)-JVS(967)*X(193)-JVS(968)*X(194)-JVS(969)*X(196)&
             &-JVS(970)*X(197)-JVS(971)*X(198)-JVS(972)*X(201)-JVS(973)*X(202)-JVS(974)*X(203)-JVS(975)*X(204)-JVS(976)&
             &*X(206)-JVS(977)*X(219)-JVS(978)*X(222)-JVS(979)*X(226)-JVS(980)*X(227)-JVS(981)*X(228)-JVS(982)*X(229)&
             &-JVS(983)*X(230)-JVS(984)*X(231)-JVS(985)*X(232))/(JVS(963))
  X(188) = (X(188)-JVS(937)*X(222)-JVS(938)*X(227)-JVS(939)*X(228)-JVS(940)*X(230)-JVS(941)*X(232))/(JVS(936))
  X(187) = (X(187)-JVS(923)*X(204)-JVS(924)*X(210)-JVS(925)*X(212)-JVS(926)*X(214)-JVS(927)*X(221)-JVS(928)*X(222)&
             &-JVS(929)*X(227)-JVS(930)*X(228)-JVS(931)*X(229)-JVS(932)*X(230)-JVS(933)*X(232))/(JVS(922))
  X(186) = (X(186)-JVS(904)*X(222)-JVS(905)*X(227)-JVS(906)*X(228)-JVS(907)*X(230)-JVS(908)*X(232))/(JVS(903))
  X(185) = (X(185)-JVS(887)*X(186)-JVS(888)*X(190)-JVS(889)*X(194)-JVS(890)*X(196)-JVS(891)*X(201)-JVS(892)*X(204)&
             &-JVS(893)*X(206)-JVS(894)*X(222)-JVS(895)*X(226)-JVS(896)*X(227)-JVS(897)*X(228)-JVS(898)*X(229)-JVS(899)&
             &*X(230)-JVS(900)*X(231)-JVS(901)*X(232))/(JVS(886))
  X(184) = (X(184)-JVS(871)*X(222)-JVS(872)*X(227)-JVS(873)*X(228)-JVS(874)*X(230)-JVS(875)*X(232))/(JVS(870))
  X(183) = (X(183)-JVS(862)*X(190)-JVS(863)*X(204)-JVS(864)*X(227)-JVS(865)*X(228)-JVS(866)*X(230)-JVS(867)*X(232))&
             &/(JVS(861))
  X(182) = (X(182)-JVS(850)*X(192)-JVS(851)*X(200)-JVS(852)*X(202)-JVS(853)*X(227)-JVS(854)*X(228)-JVS(855)*X(230)&
             &-JVS(856)*X(232))/(JVS(849))
  X(181) = (X(181)-JVS(840)*X(227)-JVS(841)*X(228)-JVS(842)*X(230)-JVS(843)*X(232))/(JVS(839))
  X(180) = (X(180)-JVS(830)*X(195)-JVS(831)*X(205)-JVS(832)*X(221)-JVS(833)*X(227)-JVS(834)*X(228)-JVS(835)*X(229)&
             &-JVS(836)*X(230)-JVS(837)*X(232))/(JVS(829))
  X(179) = (X(179)-JVS(816)*X(197)-JVS(817)*X(198)-JVS(818)*X(217)-JVS(819)*X(227)-JVS(820)*X(228)-JVS(821)*X(229)&
             &-JVS(822)*X(230)-JVS(823)*X(232))/(JVS(815))
  X(178) = (X(178)-JVS(807)*X(222)-JVS(808)*X(227)-JVS(809)*X(228)-JVS(810)*X(230)-JVS(811)*X(232))/(JVS(806))
  X(177) = (X(177)-JVS(801)*X(218)-JVS(802)*X(227)-JVS(803)*X(230)-JVS(804)*X(232))/(JVS(800))
  X(176) = (X(176)-JVS(796)*X(214)-JVS(797)*X(227)-JVS(798)*X(230)-JVS(799)*X(232))/(JVS(795))
  X(175) = (X(175)-JVS(792)*X(222)-JVS(793)*X(228)-JVS(794)*X(230))/(JVS(791))
  X(174) = (X(174)-JVS(787)*X(227)-JVS(788)*X(228)-JVS(789)*X(230)-JVS(790)*X(232))/(JVS(786))
  X(173) = (X(173)-JVS(780)*X(222)-JVS(781)*X(227)-JVS(782)*X(228)-JVS(783)*X(229)-JVS(784)*X(230))/(JVS(779))
  X(172) = (X(172)-JVS(770)*X(222)-JVS(771)*X(227)-JVS(772)*X(228)-JVS(773)*X(229)-JVS(774)*X(230))/(JVS(769))
  X(171) = (X(171)-JVS(763)*X(172)-JVS(764)*X(222)-JVS(765)*X(228)-JVS(766)*X(229)-JVS(767)*X(230))/(JVS(762))
  X(170) = (X(170)-JVS(756)*X(205)-JVS(757)*X(227)-JVS(758)*X(228)-JVS(759)*X(230)-JVS(760)*X(232))/(JVS(755))
  X(169) = (X(169)-JVS(747)*X(207)-JVS(748)*X(224)-JVS(749)*X(227)-JVS(750)*X(228)-JVS(751)*X(229)-JVS(752)*X(230)&
             &-JVS(753)*X(232))/(JVS(746))
  X(168) = (X(168)-JVS(738)*X(222)-JVS(739)*X(227)-JVS(740)*X(228)-JVS(741)*X(230)-JVS(742)*X(232))/(JVS(737))
  X(167) = (X(167)-JVS(728)*X(192)-JVS(729)*X(200)-JVS(730)*X(202)-JVS(731)*X(227)-JVS(732)*X(228)-JVS(733)*X(230)&
             &-JVS(734)*X(232))/(JVS(727))
  X(166) = (X(166)-JVS(719)*X(222)-JVS(720)*X(227)-JVS(721)*X(228)-JVS(722)*X(230)-JVS(723)*X(232))/(JVS(718))
  X(165) = (X(165)-JVS(711)*X(186)-JVS(712)*X(222)-JVS(713)*X(227)-JVS(714)*X(228)-JVS(715)*X(230)-JVS(716)*X(232))&
             &/(JVS(710))
  X(164) = (X(164)-JVS(702)*X(222)-JVS(703)*X(228)-JVS(704)*X(230))/(JVS(701))
  X(163) = (X(163)-JVS(697)*X(227)-JVS(698)*X(228)-JVS(699)*X(230)-JVS(700)*X(232))/(JVS(696))
  X(162) = (X(162)-JVS(692)*X(222)-JVS(693)*X(228)-JVS(694)*X(230))/(JVS(691))
  X(161) = (X(161)-JVS(685)*X(175)-JVS(686)*X(188)-JVS(687)*X(227)-JVS(688)*X(228)-JVS(689)*X(230)-JVS(690)*X(232))&
             &/(JVS(684))
  X(160) = (X(160)-JVS(679)*X(204)-JVS(680)*X(227)-JVS(681)*X(228))/(JVS(678))
  X(159) = (X(159)-JVS(675)*X(222)-JVS(676)*X(228)-JVS(677)*X(230))/(JVS(674))
  X(158) = (X(158)-JVS(669)*X(190)-JVS(670)*X(227)-JVS(671)*X(228)-JVS(672)*X(230)-JVS(673)*X(232))/(JVS(668))
  X(157) = (X(157)-JVS(662)*X(212)-JVS(663)*X(227)-JVS(664)*X(228)-JVS(665)*X(230)-JVS(666)*X(232))/(JVS(661))
  X(156) = (X(156)-JVS(655)*X(227)-JVS(656)*X(228)-JVS(657)*X(230)-JVS(658)*X(232))/(JVS(654))
  X(155) = (X(155)-JVS(649)*X(227)-JVS(650)*X(228)-JVS(651)*X(230)-JVS(652)*X(232))/(JVS(648))
  X(154) = (X(154)-JVS(635)*X(163)-JVS(636)*X(181)-JVS(637)*X(202)-JVS(638)*X(227)-JVS(639)*X(228)-JVS(640)*X(230)&
             &-JVS(641)*X(232))/(JVS(634))
  X(153) = (X(153)-JVS(629)*X(222)-JVS(630)*X(228)-JVS(631)*X(230))/(JVS(628))
  X(152) = (X(152)-JVS(623)*X(174)-JVS(624)*X(227)-JVS(625)*X(228)-JVS(626)*X(230)-JVS(627)*X(232))/(JVS(622))
  X(151) = (X(151)-JVS(615)*X(175)-JVS(616)*X(227)-JVS(617)*X(228)-JVS(618)*X(230)-JVS(619)*X(232))/(JVS(614))
  X(150) = (X(150)-JVS(604)*X(163)-JVS(605)*X(181)-JVS(606)*X(192)-JVS(607)*X(202)-JVS(608)*X(227)-JVS(609)*X(228)&
             &-JVS(610)*X(230)-JVS(611)*X(232))/(JVS(603))
  X(149) = (X(149)-JVS(592)*X(176)-JVS(593)*X(225)-JVS(594)*X(227)-JVS(595)*X(228)-JVS(596)*X(230)-JVS(597)*X(232))&
             &/(JVS(591))
  X(148) = (X(148)-JVS(581)*X(200)-JVS(582)*X(202)-JVS(583)*X(208)-JVS(584)*X(213)-JVS(585)*X(227)-JVS(586)*X(228)&
             &-JVS(587)*X(230)-JVS(588)*X(232))/(JVS(580))
  X(147) = (X(147)-JVS(574)*X(153)-JVS(575)*X(227)-JVS(576)*X(230)-JVS(577)*X(232))/(JVS(573))
  X(146) = (X(146)-JVS(569)*X(164)-JVS(570)*X(227)-JVS(571)*X(230)-JVS(572)*X(232))/(JVS(568))
  X(145) = (X(145)-JVS(564)*X(162)-JVS(565)*X(227)-JVS(566)*X(230)-JVS(567)*X(232))/(JVS(563))
  X(144) = (X(144)-JVS(559)*X(175)-JVS(560)*X(227)-JVS(561)*X(230)-JVS(562)*X(232))/(JVS(558))
  X(143) = (X(143)-JVS(555)*X(222)-JVS(556)*X(228)-JVS(557)*X(230))/(JVS(554))
  X(142) = (X(142)-JVS(549)*X(222)-JVS(550)*X(227)-JVS(551)*X(228)-JVS(552)*X(229)-JVS(553)*X(230))/(JVS(548))
  X(141) = (X(141)-JVS(543)*X(155)-JVS(544)*X(200)-JVS(545)*X(227)-JVS(546)*X(228))/(JVS(542))
  X(140) = (X(140)-JVS(538)*X(213)-JVS(539)*X(227)-JVS(540)*X(228)-JVS(541)*X(232))/(JVS(537))
  X(139) = (X(139)-JVS(530)*X(177)-JVS(531)*X(223)-JVS(532)*X(227)-JVS(533)*X(228)-JVS(534)*X(232))/(JVS(529))
  X(138) = (X(138)-JVS(523)*X(192)-JVS(524)*X(202)-JVS(525)*X(227)-JVS(526)*X(228))/(JVS(522))
  X(137) = (X(137)-JVS(518)*X(207)-JVS(519)*X(224)-JVS(520)*X(228)-JVS(521)*X(232))/(JVS(517))
  X(136) = (X(136)-JVS(513)*X(176)-JVS(514)*X(225)-JVS(515)*X(228)-JVS(516)*X(232))/(JVS(512))
  X(135) = (X(135)-JVS(509)*X(160)-JVS(510)*X(228)-JVS(511)*X(230))/(JVS(508))
  X(134) = (X(134)-JVS(504)*X(225)-JVS(505)*X(227)-JVS(506)*X(228))/(JVS(503))
  X(133) = (X(133)-JVS(499)*X(155)-JVS(500)*X(200)-JVS(501)*X(228)-JVS(502)*X(232))/(JVS(498))
  X(132) = (X(132)-JVS(490)*X(159)-JVS(491)*X(168)-JVS(492)*X(184)-JVS(493)*X(222)-JVS(494)*X(227)-JVS(495)*X(228)&
             &-JVS(496)*X(230)-JVS(497)*X(232))/(JVS(489))
  X(131) = (X(131)-JVS(482)*X(205)-JVS(483)*X(227)-JVS(484)*X(228)-JVS(485)*X(232))/(JVS(481))
  X(130) = (X(130)-JVS(475)*X(192)-JVS(476)*X(202)-JVS(477)*X(228)-JVS(478)*X(232))/(JVS(474))
  X(129) = (X(129)-JVS(470)*X(170)-JVS(471)*X(195)-JVS(472)*X(228)-JVS(473)*X(232))/(JVS(469))
  X(128) = (X(128)-JVS(466)*X(174)-JVS(467)*X(227)-JVS(468)*X(228))/(JVS(465))
  X(127) = (X(127)-JVS(460)*X(171)-JVS(461)*X(172)-JVS(462)*X(173)-JVS(463)*X(228)-JVS(464)*X(229))/(JVS(459))
  X(126) = (X(126)-JVS(456)*X(219)-JVS(457)*X(228)-JVS(458)*X(229))/(JVS(455))
  X(125) = (X(125)-JVS(452)*X(168)-JVS(453)*X(228)-JVS(454)*X(232))/(JVS(451))
  X(124) = (X(124)-JVS(448)*X(184)-JVS(449)*X(228)-JVS(450)*X(232))/(JVS(447))
  X(123) = (X(123)-JVS(444)*X(190)-JVS(445)*X(228)-JVS(446)*X(232))/(JVS(443))
  X(122) = (X(122)-JVS(440)*X(221)-JVS(441)*X(228)-JVS(442)*X(229))/(JVS(439))
  X(121) = (X(121)-JVS(436)*X(228)-JVS(437)*X(229)-JVS(438)*X(230))/(JVS(435))
  X(120) = (X(120)-JVS(431)*X(181)-JVS(432)*X(227)-JVS(433)*X(228))/(JVS(430))
  X(119) = (X(119)-JVS(427)*X(228)-JVS(428)*X(229)-JVS(429)*X(230))/(JVS(426))
  X(118) = (X(118)-JVS(422)*X(163)-JVS(423)*X(227)-JVS(424)*X(228))/(JVS(421))
  X(117) = (X(117)-JVS(418)*X(217)-JVS(419)*X(228)-JVS(420)*X(232))/(JVS(417))
  X(116) = (X(116)-JVS(414)*X(197)-JVS(415)*X(228)-JVS(416)*X(229))/(JVS(413))
  X(115) = (X(115)-JVS(410)*X(216)-JVS(411)*X(228)-JVS(412)*X(232))/(JVS(409))
  X(114) = (X(114)-JVS(406)*X(198)-JVS(407)*X(228)-JVS(408)*X(232))/(JVS(405))
  X(113) = (X(113)-JVS(402)*X(178)-JVS(403)*X(228)-JVS(404)*X(232))/(JVS(401))
  X(112) = (X(112)-JVS(398)*X(191)-JVS(399)*X(228)-JVS(400)*X(229))/(JVS(397))
  X(111) = (X(111)-JVS(394)*X(203)-JVS(395)*X(228)-JVS(396)*X(229))/(JVS(393))
  X(110) = (X(110)-JVS(390)*X(227)-JVS(391)*X(228)-JVS(392)*X(230))/(JVS(389))
  X(109) = (X(109)-JVS(386)*X(228)-JVS(387)*X(230))/(JVS(385))
  X(108) = (X(108)-JVS(381)*X(228)-JVS(382)*X(230))/(JVS(380))
  X(107) = (X(107)-JVS(373)*X(159)-JVS(374)*X(214)-JVS(375)*X(222)-JVS(376)*X(227)-JVS(377)*X(228)-JVS(378)*X(230)&
             &-JVS(379)*X(232))/(JVS(372))
  X(106) = (X(106)-JVS(363)*X(159)-JVS(364)*X(164)-JVS(365)*X(194)-JVS(366)*X(201)-JVS(367)*X(222)-JVS(368)*X(228)&
             &-JVS(369)*X(232))/(JVS(362))
  X(105) = (X(105)-JVS(359)*X(199)-JVS(360)*X(228)-JVS(361)*X(232))/(JVS(358))
  X(104) = (X(104)-JVS(355)*X(174)-JVS(356)*X(228)-JVS(357)*X(232))/(JVS(354))
  X(103) = (X(103)-JVS(351)*X(205)-JVS(352)*X(227)-JVS(353)*X(228))/(JVS(350))
  X(102) = (X(102)-JVS(346)*X(207)-JVS(347)*X(224)-JVS(348)*X(227)-JVS(349)*X(228))/(JVS(345))
  X(101) = (X(101)-JVS(342)*X(211)-JVS(343)*X(228)-JVS(344)*X(232))/(JVS(341))
  X(100) = (X(100)-JVS(338)*X(210)-JVS(339)*X(228)-JVS(340)*X(232))/(JVS(337))
  X(99) = (X(99)-JVS(334)*X(210)-JVS(335)*X(227)-JVS(336)*X(228))/(JVS(333))
  X(98) = (X(98)-JVS(330)*X(215)-JVS(331)*X(227)-JVS(332)*X(228))/(JVS(329))
  X(97) = (X(97)-JVS(326)*X(214)-JVS(327)*X(227)-JVS(328)*X(228))/(JVS(325))
  X(96) = (X(96)-JVS(322)*X(184)-JVS(323)*X(227)-JVS(324)*X(228))/(JVS(321))
  X(95) = (X(95)-JVS(318)*X(212)-JVS(319)*X(227)-JVS(320)*X(228))/(JVS(317))
  X(94) = (X(94)-JVS(314)*X(183)-JVS(315)*X(228)-JVS(316)*X(232))/(JVS(313))
  X(93) = (X(93)-JVS(310)*X(204)-JVS(311)*X(228)-JVS(312)*X(232))/(JVS(309))
  X(92) = (X(92)-JVS(306)*X(209)-JVS(307)*X(227)-JVS(308)*X(228))/(JVS(305))
  X(91) = (X(91)-JVS(302)*X(181)-JVS(303)*X(228)-JVS(304)*X(232))/(JVS(301))
  X(90) = (X(90)-JVS(298)*X(211)-JVS(299)*X(227)-JVS(300)*X(228))/(JVS(297))
  X(89) = (X(89)-JVS(294)*X(163)-JVS(295)*X(228)-JVS(296)*X(232))/(JVS(293))
  X(88) = (X(88)-JVS(290)*X(188)-JVS(291)*X(228)-JVS(292)*X(232))/(JVS(289))
  X(87) = (X(87)-JVS(286)*X(198)-JVS(287)*X(228)-JVS(288)*X(229))/(JVS(285))
  X(86) = (X(86)-JVS(282)*X(197)-JVS(283)*X(228)-JVS(284)*X(232))/(JVS(281))
  X(85) = (X(85)-JVS(278)*X(156)-JVS(279)*X(227)-JVS(280)*X(228))/(JVS(277))
  X(84) = (X(84)-JVS(274)*X(191)-JVS(275)*X(228)-JVS(276)*X(232))/(JVS(273))
  X(83) = (X(83)-JVS(270)*X(166)-JVS(271)*X(228)-JVS(272)*X(232))/(JVS(269))
  X(82) = (X(82)-JVS(266)*X(213)-JVS(267)*X(228)-JVS(268)*X(232))/(JVS(265))
  X(81) = (X(81)-JVS(262)*X(218)-JVS(263)*X(228)-JVS(264)*X(232))/(JVS(261))
  X(80) = (X(80)-JVS(257)*X(177)-JVS(258)*X(223)-JVS(259)*X(228)-JVS(260)*X(232))/(JVS(256))
  X(79) = (X(79)-JVS(253)*X(223)-JVS(254)*X(227)-JVS(255)*X(228))/(JVS(252))
  X(78) = (X(78)-JVS(248)*X(170)-JVS(249)*X(195)-JVS(250)*X(227)-JVS(251)*X(228))/(JVS(247))
  X(77) = (X(77)-JVS(245)*X(160)-JVS(246)*X(228))/(JVS(244))
  X(76) = (X(76)-JVS(239)*X(152)-JVS(240)*X(154)-JVS(241)*X(167)-JVS(242)*X(182)-JVS(243)*X(228))/(JVS(238))
  X(75) = (X(75)-JVS(235)*X(156)-JVS(236)*X(228)-JVS(237)*X(232))/(JVS(234))
  X(74) = (X(74)-JVS(231)*X(228)-JVS(232)*X(229)-JVS(233)*X(232))/(JVS(230))
  X(73) = (X(73)-JVS(227)*X(209)-JVS(228)*X(228)-JVS(229)*X(232))/(JVS(226))
  X(72) = (X(72)-JVS(223)*X(214)-JVS(224)*X(228)-JVS(225)*X(232))/(JVS(222))
  X(71) = (X(71)-JVS(219)*X(157)-JVS(220)*X(228)-JVS(221)*X(232))/(JVS(218))
  X(70) = (X(70)-JVS(215)*X(212)-JVS(216)*X(228)-JVS(217)*X(232))/(JVS(214))
  X(69) = (X(69)-JVS(211)*X(146)-JVS(212)*X(228)-JVS(213)*X(232))/(JVS(210))
  X(68) = (X(68)-JVS(207)*X(221)-JVS(208)*X(228)-JVS(209)*X(232))/(JVS(206))
  X(67) = (X(67)-JVS(203)*X(186)-JVS(204)*X(228)-JVS(205)*X(232))/(JVS(202))
  X(66) = (X(66)-JVS(199)*X(186)-JVS(200)*X(227)-JVS(201)*X(228))/(JVS(198))
  X(65) = (X(65)-JVS(195)*X(199)-JVS(196)*X(227)-JVS(197)*X(228))/(JVS(194))
  X(64) = (X(64)-JVS(191)*X(144)-JVS(192)*X(228)-JVS(193)*X(232))/(JVS(190))
  X(63) = (X(63)-JVS(187)*X(178)-JVS(188)*X(227)-JVS(189)*X(228))/(JVS(186))
  X(62) = (X(62)-JVS(183)*X(203)-JVS(184)*X(228)-JVS(185)*X(232))/(JVS(182))
  X(61) = (X(61)-JVS(179)*X(145)-JVS(180)*X(228)-JVS(181)*X(232))/(JVS(178))
  X(60) = (X(60)-JVS(175)*X(147)-JVS(176)*X(228)-JVS(177)*X(232))/(JVS(174))
  X(59) = (X(59)-JVS(171)*X(208)-JVS(172)*X(228)-JVS(173)*X(232))/(JVS(170))
  X(58) = (X(58)-JVS(167)*X(220)-JVS(168)*X(228)-JVS(169)*X(232))/(JVS(166))
  X(57) = (X(57)-JVS(163)*X(205)-JVS(164)*X(228)-JVS(165)*X(232))/(JVS(162))
  X(56) = (X(56)-JVS(159)*X(215)-JVS(160)*X(228)-JVS(161)*X(232))/(JVS(158))
  X(55) = (X(55)-JVS(155)*X(213)-JVS(156)*X(227)-JVS(157)*X(228))/(JVS(154))
  X(54) = (X(54)-JVS(151)*X(217)-JVS(152)*X(227)-JVS(153)*X(228))/(JVS(150))
  X(53) = (X(53)-JVS(147)*X(219)-JVS(148)*X(228)-JVS(149)*X(232))/(JVS(146))
  X(52) = (X(52)-JVS(142)*X(171)-JVS(143)*X(222)-JVS(144)*X(229)-JVS(145)*X(232))/(JVS(141))
  X(51) = (X(51)-JVS(137)*X(143)-JVS(138)*X(162)-JVS(139)*X(222)-JVS(140)*X(228))/(JVS(136))
  X(50) = (X(50)-JVS(134)*X(215)-JVS(135)*X(228))/(JVS(133))
  X(49) = (X(49)-JVS(130)*X(228)-JVS(131)*X(229)-JVS(132)*X(230))/(JVS(129))
  X(48) = (X(48)-JVS(126)*X(228)-JVS(127)*X(230))/(JVS(125))
  X(47) = (X(47)-JVS(121)*X(228)-JVS(122)*X(229)-JVS(123)*X(230))/(JVS(120))
  X(46) = (X(46)-JVS(117)*X(228)-JVS(118)*X(230))/(JVS(116))
  X(45) = (X(45)-JVS(111)*X(153)-JVS(112)*X(175)-JVS(113)*X(222)-JVS(114)*X(228))/(JVS(110))
  X(44) = (X(44)-JVS(107)*X(115)-JVS(108)*X(228))/(JVS(106))
  X(43) = (X(43)-JVS(104)*X(156)-JVS(105)*X(228))/(JVS(103))
  X(42) = (X(42)-JVS(101)*X(209)-JVS(102)*X(228))/(JVS(100))
  X(41) = (X(41)-JVS(98)*X(108)-JVS(99)*X(228))/(JVS(97))
  X(40) = (X(40)-JVS(95)*X(227)-JVS(96)*X(228))/(JVS(94))
  X(39) = (X(39)-JVS(91)*X(164)-JVS(92)*X(222)-JVS(93)*X(228))/(JVS(90))
  X(38) = (X(38)-JVS(87)*X(99)-JVS(88)*X(100)-JVS(89)*X(228))/(JVS(86))
  X(37) = (X(37)-JVS(84)*X(81)-JVS(85)*X(228))/(JVS(83))
  X(36) = (X(36)-JVS(81)*X(211)-JVS(82)*X(228))/(JVS(80))
  X(35) = (X(35)-JVS(79)*X(228))/(JVS(78))
  X(34) = (X(34)-JVS(76)*X(211)-JVS(77)*X(229))/(JVS(75))
  X(33) = (X(33)-JVS(73)*X(95)-JVS(74)*X(228))/(JVS(72))
  X(32) = (X(32)-JVS(70)*X(229)-JVS(71)*X(230))/(JVS(69))
  X(31) = (X(31)-JVS(68)*X(228))/(JVS(67))
  X(30) = (X(30)-JVS(66)*X(228))/(JVS(65))
  X(29) = (X(29)-JVS(64)*X(228))/(JVS(63))
  X(28) = (X(28)-JVS(62)*X(228))/(JVS(61))
  X(27) = (X(27)-JVS(60)*X(228))/(JVS(59))
  X(26) = (X(26)-JVS(58)*X(228))/(JVS(57))
  X(25) = (X(25)-JVS(56)*X(228))/(JVS(55))
  X(24) = (X(24)-JVS(54)*X(228))/(JVS(53))
  X(23) = (X(23)-JVS(52)*X(228))/(JVS(51))
  X(22) = (X(22)-JVS(50)*X(228))/(JVS(49))
  X(21) = (X(21)-JVS(48)*X(228))/(JVS(47))
  X(20) = (X(20)-JVS(46)*X(228))/(JVS(45))
  X(19) = (X(19)-JVS(44)*X(228))/(JVS(43))
  X(18) = (X(18)-JVS(42)*X(228))/(JVS(41))
  X(17) = (X(17)-JVS(40)*X(228))/(JVS(39))
  X(16) = (X(16)-JVS(38)*X(228))/(JVS(37))
  X(15) = (X(15)-JVS(36)*X(228))/(JVS(35))
  X(14) = (X(14)-JVS(34)*X(228))/(JVS(33))
  X(13) = (X(13)-JVS(32)*X(228))/(JVS(31))
  X(12) = (X(12)-JVS(30)*X(228))/(JVS(29))
  X(11) = (X(11)-JVS(28)*X(228))/(JVS(27))
  X(10) = (X(10)-JVS(26)*X(228))/(JVS(25))
  X(9) = (X(9)-JVS(23)*X(127)-JVS(24)*X(228))/(JVS(22))
  X(8) = (X(8)-JVS(21)*X(228))/(JVS(20))
  X(7) = (X(7)-JVS(19)*X(222))/(JVS(18))
  X(6) = (X(6)-JVS(16)*X(41)-JVS(17)*X(228))/(JVS(15))
  X(5) = (X(5)-JVS(11)*X(52)-JVS(12)*X(109)-JVS(13)*X(228)-JVS(14)*X(232))/(JVS(10))
  X(4) = (X(4)-JVS(7)*X(9)-JVS(8)*X(127)-JVS(9)*X(173))/(JVS(6))
  X(3) = (X(3)-JVS(4)*X(4)-JVS(5)*X(52))/(JVS(3))
  X(2) = X(2)/JVS(2)
  X(1) = X(1)/JVS(1)
      
END SUBROUTINE crimech_KppSolve
























      SUBROUTINE crimech_WCOPY(N,X,incX,Y,incY)








      
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

      END SUBROUTINE crimech_WCOPY



      SUBROUTINE crimech_WAXPY(N,Alpha,X,incX,Y,incY)









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
      
      END SUBROUTINE crimech_WAXPY




      SUBROUTINE crimech_WSCAL(N,Alpha,X,incX)









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

      END SUBROUTINE crimech_WSCAL


      REAL(kind=dp) FUNCTION crimech_WLAMCH( C )








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
          CALL crimech_WLAMCH_ADD(ONE,Eps,Sum)
          IF (Sum.LE.ONE) GOTO 10
        END DO
        PRINT*,'ERROR IN WLAMCH. EPS < ',Eps
        RETURN
10      Eps = Eps*2
        i = i-1      
      END IF

      crimech_WLAMCH = Eps

      END FUNCTION crimech_WLAMCH
     
      SUBROUTINE crimech_WLAMCH_ADD( A, B, Sum )

      
      REAL(kind=dp) A, B, Sum
      Sum = A + B

      END SUBROUTINE crimech_WLAMCH_ADD




      SUBROUTINE crimech_SET2ZERO(N,Y)




      
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

      END SUBROUTINE crimech_SET2ZERO



      REAL(kind=dp) FUNCTION crimech_WDOT (N, DX, incX, DY, incY) 









      IMPLICIT NONE
      INTEGER :: N, incX, incY
      REAL(kind=dp) :: DX(N), DY(N) 

      INTEGER :: i, IX, IY, M, MP1, NS
                                 
      crimech_WDOT = 0.0D0 
      IF (N .LE. 0) RETURN 
      IF (incX .EQ. incY) IF (incX-1) 5,20,60 



    5 IX = 1 
      IY = 1 
      IF (incX .LT. 0) IX = (-N+1)*incX + 1 
      IF (incY .LT. 0) IY = (-N+1)*incY + 1 
      DO i = 1,N 
        crimech_WDOT = crimech_WDOT + DX(IX)*DY(IY) 
        IX = IX + incX 
        IY = IY + incY 
      END DO 
      RETURN 





   20 M = MOD(N,5) 
      IF (M .EQ. 0) GO TO 40 
      DO i = 1,M 
         crimech_WDOT = crimech_WDOT + DX(i)*DY(i) 
      END DO 
      IF (N .LT. 5) RETURN 
   40 MP1 = M + 1 
      DO i = MP1,N,5 
          crimech_WDOT = crimech_WDOT + DX(i)*DY(i) + DX(i+1)*DY(i+1) +&
                   DX(i+2)*DY(i+2) +  &
                   DX(i+3)*DY(i+3) + DX(i+4)*DY(i+4)                   
      END DO 
      RETURN 



   60 NS = N*incX 
      DO i = 1,NS,incX 
        crimech_WDOT = crimech_WDOT + DX(i)*DY(i) 
      END DO 

      END FUNCTION crimech_WDOT                                          




END MODULE crimech_Integrator
