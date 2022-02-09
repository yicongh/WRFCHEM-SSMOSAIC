!     =======================================================================
!           THIS SUBROUINTE STEPS FORWARD FOR HETEROGENEOUS OXIDATION
!     =======================================================================

      MODULE mod_STEP_HETCHEM
      
!     MODULE INTERFACE:
!     =======================================================================
      USE mod_SSMOSAIC_PARAMS, ONLY: MW_SS, NAVO, PI
      
      CONTAINS
      
!     =======================================================================
!           THIS SUBROUINTE STEPS FORWARD FOR HETEROGENEOUS OXIDATION
!     =======================================================================
      
      SUBROUTINE STEP_HETCHEM(moAERO, nCOMP, nBINS, DIAM, NUM, Gamma_OH, OH , mmdt)
      
!     DECLARATIONS:
!     =======================================================================
      IMPLICIT NONE
      
!     NUMBER OF GAS SPECIES AND
!     AEROSOL BINS:
      INTEGER,INTENT(IN) :: nCOMP
      INTEGER,INTENT(IN) :: nBINS

!     AERSOL-PHASE ARRAY:
      REAL,INTENT(INOUT) :: moAERO(nCOMP,nBINS)
      
!     FOR INTEGRATION:
      REAL :: moAERO_PREV(nCOMP,nBINS)

!     FOR INTEGRATION:
      REAL :: binMASS(nBINS)

!     AEROSOL NUMBER CONC. [cm-3]:
      REAL,INTENT(IN) :: NUM(nBINS)

!     AEROSOL DIAMETER [nm]:
      REAL,INTENT(IN) :: DIAM(nBINS)

!     OH UPTAKE COEFFICIENT:
      REAL,INTENT(IN) :: Gamma_OH

!     OH CONCENTRATION [cm-3]:
      REAL,INTENT(IN) :: OH

!     INTEGRATION TIMESTEP:
      REAL,INTENT(IN) :: mmdt

!     OH GAS-PHASE DIFFUSIVITY [m2 s-1],
!     MOLECULAR WEIGHT [g mol-1], AND
!     GAS CONSTANT:
      REAL,PARAMETER :: DG_OH = 3e-5
      REAL,PARAMETER :: MW_OH = 17.
      REAL,PARAMETER :: RGAS  = 8.314

!     FUNCTIONALIZATION AND
!     FRAGMENTATION PROBABILITIES:
      REAL,PARAMETER :: pFRAG = 0.9
      REAL,PARAMETER :: pFUNC = 1 - pFRAG

!     MEAN FREE PATH AND
!     MEAN SPEED:
      REAL :: MFP_OH
      REAL :: RMS_OH

!     KNUDSEN NUMBER,
!     FUCHS NUMBER, AND
!     OH REACTION COEFFICIENT:
      REAL :: KNUD_OH(nBINS)
      REAL :: FUCH_OH(nBINS)
      REAL :: kRXN_OH(nCOMP,nBINS)

!     TEMPERATURE [K]:
      REAL,PARAMETER :: TEMP = 298.

!     COUNTERS:
      INTEGER :: i,j,k
      
!     CALCULATE THE HET. CHEM. RATE COEFFICIENT:
!     =======================================================================
!     MEAN FREE PATH OF OH:
      MFP_OH = 3.*DG_OH/SQRT(8.*RGAS*Temp/PI/(MW_OH/1000.))

!     ROOT-MEAN-SQUARE SPEED OF OH:
      RMS_OH = SQRT(RGAS*Temp/(MW_OH/1000.)/2./PI)

!     KNUDSEN NUMBER:
      KNUD_OH = 2.*MFP_OH/(DIAM*1e-9)

!     FUCHS-SUTUGIN FACTOR:
      FUCH_OH = (1. + KNUD_OH)/(1. + KNUD_OH + 0.75/KNUD_OH)

!     HET. CHEM. RATE COEFFICIENT [# cm-3 s-1] to [ug m-3]:
      DO j = 1,nCOMP
      DO k = 1,nBINS
         kRXN_OH(j,k) = &
         Gamma_OH*(NUM(k)*PI*((DIAM(k)*1e-9)**2.)*FUCH_OH(k))*(0.25*OH*1e6*RMS_OH)* &
         MW_SS*(1e6/NAVO)*1e6
      END DO
      END DO

!     INTEGRATION:
!     =======================================================================
      moAERO_PREV(:,:) = 0.0

      DO j = 1,nCOMP
      DO k = 1,nBINS
         moAERO_PREV(j,k) = moAERO(j,k)
      END DO
      END DO

      binMASS = MAX(SUM(moAERO,1),0.5)

!     LOOP OVER SPECIES AND BINS:
      DO j = 1,nCOMP
      DO k = 1,nBINS

      IF (j.EQ.nCOMP) THEN
         moAERO(j,k) = (moAERO_PREV(j,k) + &
                        moAERO_PREV(j-1,k)*kRXN_OH(j,k)*mmdt/binMASS(k)*pFRAG)/ &
                       (1. + kRXN_OH(j,k)*mmdt/binMASS(k))

      ELSE IF (j.EQ.1) THEN
         moAERO(j,k) = (moAERO_PREV(j,k) + &
                        moAERO_PREV(j,k)*kRXN_OH(j,k)*mmdt/binMASS(k)*pFUNC +   &
                        moAERO_PREV(j+1,k)*kRXN_OH(j,k)*mmdt/binMASS(k)*pFUNC)/ &
                       (1. + kRXN_OH(j,k)*mmdt/binMASS(k))
      ELSE
         moAERO(j,k) = (moAERO_PREV(j,k) + &
                        moAERO_PREV(j+1,k)*kRXN_OH(j,k)*mmdt/binMASS(k)*pFUNC + &
                        moAERO_PREV(j-1,k)*kRXN_OH(j,k)*mmdt/binMASS(k)*pFRAG)/ &
                       (1. + kRXN_OH(j,k)*mmdt/binMASS(k))
      END IF

      END DO
      END DO

      RETURN
      END SUBROUTINE STEP_HETCHEM
      
      END MODULE mod_STEP_HETCHEM
