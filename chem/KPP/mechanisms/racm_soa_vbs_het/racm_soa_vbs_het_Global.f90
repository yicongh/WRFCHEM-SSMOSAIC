! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Global Data Module File
! 
! Generated by KPP-2.1 symbolic chemistry Kinetics PreProcessor
!       (http://www.cs.vt.edu/~asandu/Software/KPP)
! KPP is distributed under GPL, the general public licence
!       (http://www.gnu.org/copyleft/gpl.html)
! (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa
! (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech
!     With important contributions from:
!        M. Damian, Villanova University, USA
!        R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany
! 
! File                 : racm_soa_vbs_het_Global.f90
! Time                 : Tue Jan 11 14:33:08 2022
! Working directory    : /jathar-scratch/yicongh/WRF_GoAmazon_gasaerochem.DEV/chem/KPP/mechanisms/racm_soa_vbs_het
! Equation file        : racm_soa_vbs_het.kpp
! Output root filename : racm_soa_vbs_het
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE racm_soa_vbs_het_Global

  USE racm_soa_vbs_het_Parameters, ONLY: dp, NSPEC, NVAR, NFIX, NREACT
  PUBLIC
  SAVE


! Declaration of global variables

! C - Concentration of all species
  REAL(kind=dp) :: C(NSPEC)
! VAR - Concentrations of variable species (global)
  REAL(kind=dp) :: VAR(NVAR)
! FIX - Concentrations of fixed species (global)
  REAL(kind=dp) :: FIX(NFIX)
! VAR, FIX are chunks of array C
      EQUIVALENCE( C(1),VAR(1) )
      EQUIVALENCE( C(98),FIX(1) )
! RCONST - Rate constants (global)
  REAL(kind=dp) :: RCONST(NREACT)
! TIME - Current integration time
  REAL(kind=dp) :: TIME
! SUN - Sunlight intensity between [0,1]
  REAL(kind=dp) :: SUN
! TEMP - Temperature
  REAL(kind=dp) :: TEMP
! RTOLS - (scalar) Relative tolerance
  REAL(kind=dp) :: RTOLS
! TSTART - Integration start time
  REAL(kind=dp) :: TSTART
! TEND - Integration end time
  REAL(kind=dp) :: TEND
! DT - Integration step
  REAL(kind=dp) :: DT
! ATOL - Absolute tolerance
  REAL(kind=dp) :: ATOL(NSPEC)
! RTOL - Relative tolerance
  REAL(kind=dp) :: RTOL(NSPEC)
! STEPMIN - Lower bound for integration step
  REAL(kind=dp) :: STEPMIN
! STEPMAX - Upper bound for integration step
  REAL(kind=dp) :: STEPMAX
! CFACTOR - Conversion factor for concentration units
  REAL(kind=dp) :: CFACTOR

! INLINED global variable declarations

! INLINED global variable declarations


END MODULE racm_soa_vbs_het_Global

