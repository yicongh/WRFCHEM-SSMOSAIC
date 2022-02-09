!     =======================================================================
!            THIS MODULE CONTAINS THE GLOBAL PARAMETERS FOR SS-MOSAIC
!     =======================================================================

      MODULE mod_SSMOSAIC_PARAMS
      
!     MODULE INTERFACE:
!     =======================================================================
      USE MODULE_STATE_DESCRIPTION, ONLY: num_moist, num_chem
      
!     ACCESS THE POINTERS FOR CHEM IN WRF:
      USE MODULE_CONFIGURE, ONLY: & 
      p_TERP_CN3_g, p_TERP_CN2_g, p_TERP_CN1_g, p_TERP_C00_g, p_TERP_C01_g, &
      p_TERP_C02_g, p_TERP_C03_g, p_TERP_C04_g, p_TERP_C05_g, p_TERP_C06_g, &
      p_TERP_C07_g, p_TERP_C08_g, p_TERP_C09_g, &
      p_TERP_CN3_a01, p_TERP_CN3_a02, p_TERP_CN3_a03, p_TERP_CN3_a04, &
      p_TERP_CN2_a01, p_TERP_CN2_a02, p_TERP_CN2_a03, p_TERP_CN2_a04, &
      p_TERP_CN1_a01, p_TERP_CN1_a02, p_TERP_CN1_a03, p_TERP_CN1_a04, &
      p_TERP_C00_a01, p_TERP_C00_a02, p_TERP_C00_a03, p_TERP_C00_a04, &
      p_TERP_C01_a01, p_TERP_C01_a02, p_TERP_C01_a03, p_TERP_C01_a04, &
      p_TERP_C02_a01, p_TERP_C02_a02, p_TERP_C02_a03, p_TERP_C02_a04, &
      p_TERP_C03_a01, p_TERP_C03_a02, p_TERP_C03_a03, p_TERP_C03_a04, &
      p_TERP_C04_a01, p_TERP_C04_a02, p_TERP_C04_a03, p_TERP_C04_a04, &
      p_TERP_C05_a01, p_TERP_C05_a02, p_TERP_C05_a03, p_TERP_C05_a04, &
      p_TERP_C06_a01, p_TERP_C06_a02, p_TERP_C06_a03, p_TERP_C06_a04, &
      p_TERP_C07_a01, p_TERP_C07_a02, p_TERP_C07_a03, p_TERP_C07_a04, &
      p_TERP_C08_a01, p_TERP_C08_a02, p_TERP_C08_a03, p_TERP_C08_a04, &
      p_TERP_C09_a01, p_TERP_C09_a02, p_TERP_C09_a03, p_TERP_C09_a04, &      
      p_TERP_CN3_d01, p_TERP_CN3_d02, p_TERP_CN3_d03, p_TERP_CN3_d04, &
      p_TERP_CN2_d01, p_TERP_CN2_d02, p_TERP_CN2_d03, p_TERP_CN2_d04, &
      p_TERP_CN1_d01, p_TERP_CN1_d02, p_TERP_CN1_d03, p_TERP_CN1_d04, &
      p_TERP_C00_d01, p_TERP_C00_d02, p_TERP_C00_d03, p_TERP_C00_d04, &
      p_TERP_C01_d01, p_TERP_C01_d02, p_TERP_C01_d03, p_TERP_C01_d04, &
      p_TERP_C02_d01, p_TERP_C02_d02, p_TERP_C02_d03, p_TERP_C02_d04, &
      p_TERP_C03_d01, p_TERP_C03_d02, p_TERP_C03_d03, p_TERP_C03_d04, &
      p_TERP_C04_d01, p_TERP_C04_d02, p_TERP_C04_d03, p_TERP_C04_d04, &
      p_TERP_C05_d01, p_TERP_C05_d02, p_TERP_C05_d03, p_TERP_C05_d04, &
      p_TERP_C06_d01, p_TERP_C06_d02, p_TERP_C06_d03, p_TERP_C06_d04, &
      p_TERP_C07_d01, p_TERP_C07_d02, p_TERP_C07_d03, p_TERP_C07_d04, &
      p_TERP_C08_d01, p_TERP_C08_d02, p_TERP_C08_d03, p_TERP_C08_d04, &
      p_TERP_C09_d01, p_TERP_C09_d02, p_TERP_C09_d03, p_TERP_C09_d04, &
      p_HO, p_NUM_a01, p_NUM_a02, p_NUM_a03, p_NUM_a04
      
!     DECLARATIONS:
!     =======================================================================
!     MOLECULAR WEIGHT OF 
!     SIMPLE-SOM SPECIES [g mol-1]:
      REAL,PARAMETER :: MW_SS = 250.      

!     NUMBER OF ACTIVE SPECIES:
      INTEGER,PARAMETER :: nCOMP = 13

!     NUMBER OF AEROSOL SIZE BINS:
      INTEGER,PARAMETER :: nBINS = 4
      
!     SATURATION CONC. [ug m-3]:
      REAL,PARAMETER,DIMENSION(nCOMP) :: &
      CSAT = (/1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9/)
      
!     PI VALUE:
      REAL,PARAMETER :: PI = 3.1415926

!     AVOGADRO'S NUMBER:
      REAL,PARAMETER :: NAVO = 6.022e23
      
      END MODULE mod_SSMOSAIC_PARAMS
