!     =======================================================================
!                 THIS IS THE INTERFACE FOR THE SS-MOSAIC MODEL
!     =======================================================================

      MODULE mod_SSMOSAIC

!     MODULE INTERFACE:            
!     =======================================================================
      USE mod_SSMOSAIC_PARAMS, &
      ONLY: nCOMP, nBINS, num_chem, p_HO
      
      USE mod_STEP_SOACOND, ONLY: STEP_SOACOND
      USE mod_STEP_OLIG,    ONLY: STEP_OLIG
      USE mod_STEP_HETCHEM, ONLY: STEP_HETCHEM

      USE mod_MAP_INOUT, ONLY: MAP_INOUT
      
      CONTAINS
      
!     =======================================================================
!                          THIS IS THE DRIVER FOR MOSAIC
!     =======================================================================

      SUBROUTINE drive_SSMOSAIC(CHEM, iD1,iD2, jD1,jD2, kD1,kD2, &
                                      iM1,iM2, jM1,jM2, kM1,kM2, & 
                                      iT1,iT2, jT1,jT2, kT1,kT2, mdt)
      
      IMPLICIT NONE

      INTEGER,INTENT(IN) :: iD1, iD2, jD1, jD2, kD1, kD2
      INTEGER,INTENT(IN) :: iM1, iM2, jM1, jM2, kM1, kM2
      INTEGER,INTENT(IN) :: iT1, iT2, jT1, jT2, kT1, kT2

!     WRF ARRAY FOR TRACER SPECIES 
!     FOR GAS [ppm] AND FOR AEROSOL SPECIES [ug per kg-air]:
      REAL,INTENT(INOUT), &
      DIMENSION( iM1:iM2, kM1:kM2, jM1:jM2, 1:num_chem ) :: CHEM
      
!     MODEL TIMESTEP,
!     INTERNAL CLOCK, AND 
!     INTEGRATION TIMESTEP [s]:
      REAL,INTENT(IN) :: mdt
      REAL :: clock
      REAL :: mmdt
      
!     MOSAIC GAS [ug m-3], 
!     AEROSOL [ug m-3] AND OLIGOMER [ug m-3] ARRAYS:
      REAL :: moGAS(nCOMP)
      REAL :: moAERO(nCOMP,nBINS)
      REAL :: moOLIG(nCOMP,nBINS)

!     AEROSOl SIZE BIN DIAMETERS [nm]:
      REAL,DIMENSION(nBINS) :: &
      DIAM = (/10.,100.,1000.,2000./)

!     AEROSOL NUMBER CONC. [cm-3]:
      REAL,DIMENSION(nBINS) :: &
      NUM = (/15000.,15000.,15000.,15000./)
      
!     OLIGOMERIZATION FORWARD [cm3 s-1] AND
!     REVERSE RATES [s-1]:
      REAL :: kf = 0.0 !1e-24
      REAL :: kr = 0.0 !1.6e-2

!     HET. OH UPTAKE COEFFICIENT:
      REAL :: Gamma_OH = 0.0 !1.0

!     OH CONCENTRATION:
      REAL :: OH

!     COUNTERS:
      INTEGER :: i,j,k,ii,jj,kk
      
!     LOOP OVER TILES:
!     =======================================================================
      LOOP1: DO i = iT1,iT2
      LOOP2: DO j = jT1,jT2
      LOOP3: DO k = kT1,kT2

!     FLUSH ARRAYS:
!     =======================================================================
      moGAS(:) = 0.0
      
      moAERO(:,:) = 0.0
      moOLIG(:,:) = 0.0

!     MAPPING IN:
!     =======================================================================
      CALL MAP_INOUT(i,k,j,0,CHEM, &
                     iM1,iM2, jM1,jM2, kM1,kM2, num_chem, &
                     moGAS, moAERO, moOLIG, nCOMP, nBINS)

!     LOOP OVER TIMESTEPS FOR PARTITIONING:
!     =======================================================================
!     ZERO CLOCK:
      clock = 0.0

      LOOP6: DO WHILE (clock < mdt)

!     CALCULATE INTEGRATION TIMESTEP (ASTEM):
      mmdt = mdt/1.0; clock = clock + mmdt
      
!     STEP FOR SOA CONDENSATION:
!     =======================================================================
      CALL STEP_SOACOND(moGAS, moAERO, nCOMP, nBINS, DIAM, NUM, mmdt)
      
!     STEP FOR OLIGOMERIZATION:
!     =======================================================================
      CALL STEP_OLIG(moAERO, moOLIG, nCOMP, nBINS, DIAM, NUM, kf, kr, mmdt)

!     STEP FOR HETEROGENEOUS OXIDATION:
!     =======================================================================
      OH = CHEM(i,k,j,p_HO)*1e-6*(101325./8.314/298.)*6.022e23*1e-6

      CALL STEP_HETCHEM(moAERO, nCOMP, nBINS, DIAM, NUM, Gamma_OH, OH , mmdt)

!     MAPPING OUT:
!     =======================================================================
      CALL MAP_INOUT(i,k,j,1,CHEM, &
                     iM1,iM2, jM1,jM2, kM1,kM2, num_chem, &
                     moGAS, moAERO, moOLIG, nCOMP, nBINS)
      
      END DO LOOP6

      END DO LOOP3
      END DO LOOP2
      END DO LOOP1      

      RETURN
      END SUBROUTINE
      
      END MODULE mod_SSMOSAIC
