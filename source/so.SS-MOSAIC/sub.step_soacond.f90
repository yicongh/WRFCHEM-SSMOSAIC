!     =======================================================================
!               THIS SUBROUINTE STEPS FORWARD FOR SOA CONDENSATION
!     =======================================================================

      MODULE mod_STEP_SOACOND
      
!     MODULE INTERFACE:
!     =======================================================================      
!     POINTERS FOR CHEM IN WRF:
      USE mod_SSMOSAIC_PARAMS, ONLY: CSAT,PI
      
      CONTAINS
      
!     =======================================================================
!               THIS SUBROUINTE STEPS FORWARD FOR SOA CONDENSATION
!     =======================================================================
      
      SUBROUTINE STEP_SOACOND(moGAS, moAERO, nCOMP, nBINS, DIAM, NUM, mmdt)
      
!     DECLARATIONS:
!     =======================================================================
      IMPLICIT NONE
      
!     NUMBER OF GAS SPECIES AND
!     AEROSOL BINS:
      INTEGER,INTENT(IN) :: nCOMP
      INTEGER,INTENT(IN) :: nBINS

!     GAS- AND AERSOL-PHASE
!     ARRAYS:
      REAL,INTENT(INOUT) :: moGAS(nCOMP)
      REAL,INTENT(INOUT) :: moAERO(nCOMP,nBINS)

!     AEROSOl SIZE BIN DIAMETERS [nm]:
      REAL,DIMENSION(nBINS),INTENT(INOUT) :: DIAM

!     AEROSOL NUMBER CONC. [cm-3]:
      REAL,DIMENSION(nBINS),INTENT(IN) :: NUM

!     INTEGRATION TIMESTEP:
      REAL,INTENT(IN) :: mmdt

!     COUNTERS:
      INTEGER :: i,j,k

!     SATURATION RATIO:
      REAL,DIMENSION(nBINS) :: SRATIO
      
!     GAS- AND PARTICLE-PHASE 
!     MASS-TRANSFER COEFF. [m s-1]:
      REAL :: KTg
      REAL :: KTp

!     OVERALL MASS-TRANSFER COEFF. 
!     ACCOUNTING FOR AREA [s-1]:
      REAL,DIMENSION(nBINS) :: KOTAg
      
!     GAS- AND PARTICLE-PHASE
!     DIFFUSIVITIES [m2 s-1]:
      REAL :: DIFFg = 4e-6
      REAL :: DIFFb = 1e-10

!     FOR INTEGRATION:      
      REAL :: TOT, SUM1, SUM2
      REAL :: GAS_NEXT
      REAL :: AERO_NEXT(nBINS)      
      
!     STEP FORWARD:
!     =======================================================================
      LOOP1: DO j = 1,nCOMP
         
!     SATURATION RATIO:
      DO k = 1,nBINS
         SRATIO(k) = CSAT(j)/MAX(1e-3,2.0 + SUM(moAERO(:,k)))  !!! ASSUME BACKGROUND OA MASS ??
      END DO

!     MASS-TRANSFER COEFFICIENT:
      DO k = 1,nBINS
         
         KTg = DIFFg/(0.5*DIAM(k)*1e-9)
         KTp = DIFFb/(0.5*DIAM(k)*1e-9)/5.0
         
         KOTAg(k) = 1.0/(1.0/KTg + 1.0/KTp*SRATIO(k))* &
                    NUM(k)*1e6*PI*((DIAM(k)*1e-9)**2.0)
      END DO
      
!     TOTAL SPECIES MASS:
      TOT = moGAS(j) + SUM(moAERO(j,1:nBINS))

!     CALCULATE SUM TERMS:
      SUM1 = 0.0
      SUM2 = 0.0

      DO k = 1,nBINS
         SUM1 = &
         SUM1 + moAERO(j,k)/(1.0 + mmdt*KOTAg(k)*SRATIO(k))

         SUM2 = &
         SUM2 + mmdt*KOTAg(k)/(1.0 + mmdt*KOTAg(k)*SRATIO(k))
      END DO      

!     GAS AT NEXT TIMESTEP:
      GAS_NEXT = (TOT - SUM1)/(1 + SUM2)

!     AEROSOL AT NEXT TIMESTEP:
      DO k = 1,nBINS
         AERO_NEXT(k) = (moAERO(j,k) + mmdt*KOTAg(k)*GAS_NEXT)/ &
                        (1.0 + mmdt*KOTAg(k)*SRATIO(k))
      END DO

!     UPDATE GAS AND AERO ARRAYS:
      moGAS(j) = GAS_NEXT
      
      DO k = 1,nBINS
         moAERO(j,k) = AERO_NEXT(k)
      END DO
      
      END DO LOOP1
      
      RETURN
      END SUBROUTINE STEP_SOACOND
      
      END MODULE mod_STEP_SOACOND
