! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Parameter Module File
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
! File                 : racm_mim_Parameters.f90
! Time                 : Tue Jan 11 14:33:02 2022
! Working directory    : /jathar-scratch/yicongh/WRF_GoAmazon_gasaerochem.DEV/chem/KPP/mechanisms/racm_mim
! Equation file        : racm_mim.kpp
! Output root filename : racm_mim
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE racm_mim_Parameters

  USE racm_mim_Precision
  PUBLIC
  SAVE


! NSPEC - Number of chemical species
  INTEGER, PARAMETER :: NSPEC = 82 
! NVAR - Number of Variable species
  INTEGER, PARAMETER :: NVAR = 80 
! NVARACT - Number of Active species
  INTEGER, PARAMETER :: NVARACT = 76 
! NFIX - Number of Fixed species
  INTEGER, PARAMETER :: NFIX = 2 
! NREACT - Number of reactions
  INTEGER, PARAMETER :: NREACT = 245 
! NVARST - Starting of variables in conc. vect.
  INTEGER, PARAMETER :: NVARST = 1 
! NFIXST - Starting of fixed in conc. vect.
  INTEGER, PARAMETER :: NFIXST = 81 
! NONZERO - Number of nonzero entries in Jacobian
  INTEGER, PARAMETER :: NONZERO = 950 
! LU_NONZERO - Number of nonzero entries in LU factoriz. of Jacobian
  INTEGER, PARAMETER :: LU_NONZERO = 1087 
! CNVAR - (NVAR+1) Number of elements in compressed row format
  INTEGER, PARAMETER :: CNVAR = 81 
! NLOOKAT - Number of species to look at
  INTEGER, PARAMETER :: NLOOKAT = 0 
! NMONITOR - Number of species to monitor
  INTEGER, PARAMETER :: NMONITOR = 0 
! NMASS - Number of atoms to check mass balance
  INTEGER, PARAMETER :: NMASS = 1 
! PI - Value of pi
  REAL(kind=dp), PARAMETER :: PI = 3.14159265358979 

! Index declaration for variable species in C and VAR
!   VAR(ind_spc) = C(ind_spc)

  INTEGER, PARAMETER :: ind_SULF = 1 
  INTEGER, PARAMETER :: ind_CO2 = 2 
  INTEGER, PARAMETER :: ind_ORA1 = 3 
  INTEGER, PARAMETER :: ind_ORA2 = 4 
  INTEGER, PARAMETER :: ind_SO2 = 5 
  INTEGER, PARAMETER :: ind_O1D = 6 
  INTEGER, PARAMETER :: ind_ISHP = 7 
  INTEGER, PARAMETER :: ind_HC5 = 8 
  INTEGER, PARAMETER :: ind_TOL = 9 
  INTEGER, PARAMETER :: ind_XYL = 10 
  INTEGER, PARAMETER :: ind_N2O5 = 11 
  INTEGER, PARAMETER :: ind_HC8 = 12 
  INTEGER, PARAMETER :: ind_MAHP = 13 
  INTEGER, PARAMETER :: ind_HC3 = 14 
  INTEGER, PARAMETER :: ind_CH4 = 15 
  INTEGER, PARAMETER :: ind_UDD = 16 
  INTEGER, PARAMETER :: ind_ETH = 17 
  INTEGER, PARAMETER :: ind_NALD = 18 
  INTEGER, PARAMETER :: ind_HNO4 = 19 
  INTEGER, PARAMETER :: ind_OP1 = 20 
  INTEGER, PARAMETER :: ind_MPAN = 21 
  INTEGER, PARAMETER :: ind_HACE = 22 
  INTEGER, PARAMETER :: ind_HONO = 23 
  INTEGER, PARAMETER :: ind_H2O2 = 24 
  INTEGER, PARAMETER :: ind_HKET = 25 
  INTEGER, PARAMETER :: ind_O3P = 26 
  INTEGER, PARAMETER :: ind_PHO = 27 
  INTEGER, PARAMETER :: ind_ADDT = 28 
  INTEGER, PARAMETER :: ind_ADDX = 29 
  INTEGER, PARAMETER :: ind_ETE = 30 
  INTEGER, PARAMETER :: ind_ADDC = 31 
  INTEGER, PARAMETER :: ind_HNO3 = 32 
  INTEGER, PARAMETER :: ind_PAA = 33 
  INTEGER, PARAMETER :: ind_ISON = 34 
  INTEGER, PARAMETER :: ind_PAN = 35 
  INTEGER, PARAMETER :: ind_API = 36 
  INTEGER, PARAMETER :: ind_CO = 37 
  INTEGER, PARAMETER :: ind_LIM = 38 
  INTEGER, PARAMETER :: ind_ISO = 39 
  INTEGER, PARAMETER :: ind_CSL = 40 
  INTEGER, PARAMETER :: ind_DIEN = 41 
  INTEGER, PARAMETER :: ind_MACP = 42 
  INTEGER, PARAMETER :: ind_TPAN = 43 
  INTEGER, PARAMETER :: ind_GLY = 44 
  INTEGER, PARAMETER :: ind_ETEP = 45 
  INTEGER, PARAMETER :: ind_OLTP = 46 
  INTEGER, PARAMETER :: ind_OLIP = 47 
  INTEGER, PARAMETER :: ind_MACR = 48 
  INTEGER, PARAMETER :: ind_KET = 49 
  INTEGER, PARAMETER :: ind_ISOP = 50 
  INTEGER, PARAMETER :: ind_MGLY = 51 
  INTEGER, PARAMETER :: ind_CSLP = 52 
  INTEGER, PARAMETER :: ind_HC8P = 53 
  INTEGER, PARAMETER :: ind_LIMP = 54 
  INTEGER, PARAMETER :: ind_HC5P = 55 
  INTEGER, PARAMETER :: ind_HCHO = 56 
  INTEGER, PARAMETER :: ind_TOLP = 57 
  INTEGER, PARAMETER :: ind_XYLP = 58 
  INTEGER, PARAMETER :: ind_APIP = 59 
  INTEGER, PARAMETER :: ind_ONIT = 60 
  INTEGER, PARAMETER :: ind_DCB = 61 
  INTEGER, PARAMETER :: ind_XO2 = 62 
  INTEGER, PARAMETER :: ind_OLT = 63 
  INTEGER, PARAMETER :: ind_ALD = 64 
  INTEGER, PARAMETER :: ind_OLI = 65 
  INTEGER, PARAMETER :: ind_OLNN = 66 
  INTEGER, PARAMETER :: ind_OLND = 67 
  INTEGER, PARAMETER :: ind_ETHP = 68 
  INTEGER, PARAMETER :: ind_KETP = 69 
  INTEGER, PARAMETER :: ind_OP2 = 70 
  INTEGER, PARAMETER :: ind_HC3P = 71 
  INTEGER, PARAMETER :: ind_ACO3 = 72 
  INTEGER, PARAMETER :: ind_HO = 73 
  INTEGER, PARAMETER :: ind_NO3 = 74 
  INTEGER, PARAMETER :: ind_HO2 = 75 
  INTEGER, PARAMETER :: ind_MO2 = 76 
  INTEGER, PARAMETER :: ind_O3 = 77 
  INTEGER, PARAMETER :: ind_TCO3 = 78 
  INTEGER, PARAMETER :: ind_NO = 79 
  INTEGER, PARAMETER :: ind_NO2 = 80 

! Index declaration for fixed species in C
!   C(ind_spc)

  INTEGER, PARAMETER :: ind_H2O = 81 
  INTEGER, PARAMETER :: ind_M = 82 

! Index declaration for fixed species in FIX
!    FIX(indf_spc) = C(ind_spc) = C(NVAR+indf_spc)

  INTEGER, PARAMETER :: indf_H2O = 1 
  INTEGER, PARAMETER :: indf_M = 2 

END MODULE racm_mim_Parameters

