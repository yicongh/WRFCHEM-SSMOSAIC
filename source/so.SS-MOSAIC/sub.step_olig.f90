!     =======================================================================
!               THIS SUBROUINTE STEPS FORWARD FOR OLIGOMERIZATION
!     =======================================================================

      MODULE mod_STEP_OLIG
      
!     MODULE INTERFACE:
!     =======================================================================
      USE mod_SSMOSAIC_PARAMS, ONLY: MW_SS, NAVO, PI
      
      CONTAINS
      
!     =======================================================================
!               THIS SUBROUINTE STEPS FORWARD FOR OLIGOMERIZATION
!     =======================================================================
      
      SUBROUTINE STEP_OLIG(moAERO, moOLIG, nCOMP, nBINS, DIAM, NUM, kf, kr, mmdt)
      
!     DECLARATIONS:
!     =======================================================================
      IMPLICIT NONE
      
!     NUMBER OF GAS SPECIES AND
!     AEROSOL BINS:
      INTEGER,INTENT(IN) :: nCOMP
      INTEGER,INTENT(IN) :: nBINS

!     AERSOL-PHASE MONOMER AND
!     OLIGOMER ARRAYS:
      REAL,INTENT(INOUT) :: moAERO(nCOMP,nBINS)
      REAL,INTENT(INOUT) :: moOLIG(nCOMP,nBINS)

!     AERSOL-PHASE MONOMER AND
!     OLIGOMER ARRAYS [cm-3]:     
      REAL :: jkAERO(nCOMP,nBINS)
      REAL :: jkOLIG(nCOMP,nBINS)
      REAL :: jkTOTO(nCOMP,nBINS)      

      REAL :: jkAERO_NXT(nCOMP,nBINS)
      REAL :: jkOLIG_NXT(nCOMP,nBINS)

      REAL :: kAERO(nBINS)

!     AEROSOL NUMBER CONC. [cm-3]:
      REAL,INTENT(IN) :: NUM(nBINS)

!     AEROSOL DIAMETER [nm]:
      REAL,INTENT(IN) :: DIAM(nBINS)

!     OLIGOMERIZATION FORMATION AND
!     DISSOCIATION RATES:
      REAL,INTENT(IN) :: kf
      REAL,INTENT(IN) :: kr

!     EFFECTIVE FORMATION 
!     RATE [cm3 s-1]:
      REAL :: kf_EFF(nBINS)

!     INTEGRATION TIMESTEP:
      REAL,INTENT(IN) :: mmdt

!     COUNTERS:
      INTEGER :: i,j,k
            
!     CONVERT UNITS AND MAP-IN:
!     =======================================================================
      jkAERO(:,:) = 0.0
      jkOLIG(:,:) = 0.0

      DO j = 1,nCOMP
      DO k = 1,nBINS
         jkAERO(j,k) = moAERO(j,k)*1e-6*1e-6/MW_SS*NAVO
         jkOLIG(j,k) = moOLIG(j,k)*1e-6*1e-6/MW_SS*NAVO
      END DO
      END DO
      
!     EFFECTIVE OLIGOMERIZATION RATE:
!     =======================================================================
      DO k = 1,nBINS
         kf_EFF(k) = kf/(NUM(k)*(1e-21*PI/6.0*(DIAM(k)**3.0)))
      END DO

!     STEP FORWARD:
!     =======================================================================
      kAERO = SUM(jkAERO,1)

      jkAERO_NXT(:,:) = 0.0
      jkOLIG_NXT(:,:) = 0.0
      
      jkTOTO = jkAERO + jkOLIG

      DO j = 1,nCOMP
      DO k = 1,nBINS
        jkAERO_NXT(j,k) = (jkAERO(j,k) + mmdt*0.5*kr*jkTOTO(j,k))/ &
                          (1.0 + mmdt*kf_EFF(k)*kAERO(k) + mmdt*0.5*kr)

        jkOLIG_NXT(j,k) = (jkTOTO(j,k)*(1.0 + mmdt*kf_EFF(k)*kAERO(k)) - jkAERO(j,k))/ &
                          (1.0 + mmdt*kf_EFF(k)*kAERO(k) + mmdt*0.5*kr)
      END DO
      END DO

!     CONVERT UNITS AND MAP-OUT:
!     =======================================================================
      DO j = 1,nCOMP
      DO k = 1,nBINS
         moAERO(j,k) = jkAERO_NXT(j,k)/(1e-6*1e-6/MW_SS*NAVO) 
         moOLIG(j,k) = jkOLIG_NXT(j,k)/(1e-6*1e-6/MW_SS*NAVO)
      END DO
      END DO      
      
      RETURN
      END SUBROUTINE STEP_OLIG
      
      END MODULE mod_STEP_OLIG
