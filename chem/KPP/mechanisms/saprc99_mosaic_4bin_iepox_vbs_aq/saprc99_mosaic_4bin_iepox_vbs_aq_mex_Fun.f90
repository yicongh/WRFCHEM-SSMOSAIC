
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE mexFunction( nlhs, plhs, nrhs, prhs )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                 Matlab Gateway for the Derivative Function Fun
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 USE saprc99_mosaic_4bin_iepox_vbs_aq_Model

      INTEGER nlhs, nrhs
      INTEGER plhs(*), prhs(*)
      INTEGER mxGetPr, mxCreateFull, mxGetM, mxgetN
      INTEGER VPtr, FPtr, RPtr, VdotPtr
      REAL(kind=dp) V(118), F(2), RCT(289)
      REAL(kind=dp) Vdot(118)

! Check for the right number of input arguments
      IF ( nrhs .ne. 3 ) THEN
         CALL mexErrMsgTxt('Fun requires 3 input vectors: &
     &V(118), F(2), RCT(289)')
      END IF 
! Check for the right number of output arguments
      IF ( nlhs .ne. 1 ) THEN
         CALL mexErrMsgTxt('Fun requires 1 output vector: &
     &Vdot(118)')
      END IF 

      plhs(1) = mxCreateDoubleMatrix(118,1,0)

      VPtr = mxGetPr(prhs(1))
      CALL mxCopyPtrToReal8(VPtr,V,118)
      
      FPtr = mxGetPr(prhs(2))
      CALL mxCopyPtrToReal8(FPtr,F,2)
      
      RPtr = mxGetPr(prhs(3))
      CALL mxCopyPtrToReal8(RPtr,RCT,289)

      VdotPtr = mxGetPr(plhs(1))

      CALL Fun( V, F, RCT, Vdot )

      CALL mxCopyReal8ToPtr(Vdot, VdotPtr, 118)

 END SUBROUTINE mexFunction
