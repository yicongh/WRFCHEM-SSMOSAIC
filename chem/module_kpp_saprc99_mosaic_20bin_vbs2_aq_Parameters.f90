























MODULE saprc99_mosaic_20bin_vbs2_aq_Parameters

  USE saprc99_mosaic_20bin_vbs2_aq_Precision
  PUBLIC
  SAVE



  INTEGER, PARAMETER :: NSPEC = 103 

  INTEGER, PARAMETER :: NVAR = 101 

  INTEGER, PARAMETER :: NVARACT = 86 

  INTEGER, PARAMETER :: NFIX = 2 

  INTEGER, PARAMETER :: NREACT = 259 

  INTEGER, PARAMETER :: NVARST = 1 

  INTEGER, PARAMETER :: NFIXST = 102 

  INTEGER, PARAMETER :: NONZERO = 1051 

  INTEGER, PARAMETER :: LU_NONZERO = 1136 

  INTEGER, PARAMETER :: CNVAR = 102 

  INTEGER, PARAMETER :: NLOOKAT = 0 

  INTEGER, PARAMETER :: NMONITOR = 0 

  INTEGER, PARAMETER :: NMASS = 1 

  REAL(kind=dp), PARAMETER :: PI = 3.14159265358979 




  INTEGER, PARAMETER :: ind_H2SO4 = 1 
  INTEGER, PARAMETER :: ind_CO2 = 2 
  INTEGER, PARAMETER :: ind_PCG1_B_C = 3 
  INTEGER, PARAMETER :: ind_PCG1_B_O = 4 
  INTEGER, PARAMETER :: ind_PCG1_F_C = 5 
  INTEGER, PARAMETER :: ind_PCG1_F_O = 6 
  INTEGER, PARAMETER :: ind_HCOOH = 7 
  INTEGER, PARAMETER :: ind_CCO_OH = 8 
  INTEGER, PARAMETER :: ind_RCO_OH = 9 
  INTEGER, PARAMETER :: ind_CCO_OOH = 10 
  INTEGER, PARAMETER :: ind_RCO_OOH = 11 
  INTEGER, PARAMETER :: ind_XN = 12 
  INTEGER, PARAMETER :: ind_XC = 13 
  INTEGER, PARAMETER :: ind_NUME = 14 
  INTEGER, PARAMETER :: ind_DEN = 15 
  INTEGER, PARAMETER :: ind_ANT1_c = 16 
  INTEGER, PARAMETER :: ind_BIOG1_c = 17 
  INTEGER, PARAMETER :: ind_BIOG2_c = 18 
  INTEGER, PARAMETER :: ind_BIOG1_o = 19 
  INTEGER, PARAMETER :: ind_PSD1 = 20 
  INTEGER, PARAMETER :: ind_PSD2 = 21 
  INTEGER, PARAMETER :: ind_OPCG1_B_O = 22 
  INTEGER, PARAMETER :: ind_PCG2_B_O = 23 
  INTEGER, PARAMETER :: ind_OPCG1_F_O = 24 
  INTEGER, PARAMETER :: ind_PCG2_F_O = 25 
  INTEGER, PARAMETER :: ind_O1D = 26 
  INTEGER, PARAMETER :: ind_PCG2_B_C = 27 
  INTEGER, PARAMETER :: ind_OPCG1_B_C = 28 
  INTEGER, PARAMETER :: ind_PCG2_F_C = 29 
  INTEGER, PARAMETER :: ind_OPCG1_F_C = 30 
  INTEGER, PARAMETER :: ind_CH4 = 31 
  INTEGER, PARAMETER :: ind_C2H6 = 32 
  INTEGER, PARAMETER :: ind_PAN = 33 
  INTEGER, PARAMETER :: ind_PAN2 = 34 
  INTEGER, PARAMETER :: ind_PBZN = 35 
  INTEGER, PARAMETER :: ind_MA_PAN = 36 
  INTEGER, PARAMETER :: ind_H2O2 = 37 
  INTEGER, PARAMETER :: ind_C3H8 = 38 
  INTEGER, PARAMETER :: ind_BACL = 39 
  INTEGER, PARAMETER :: ind_ETOH = 40 
  INTEGER, PARAMETER :: ind_N2O5 = 41 
  INTEGER, PARAMETER :: ind_SO2 = 42 
  INTEGER, PARAMETER :: ind_DMS = 43 
  INTEGER, PARAMETER :: ind_HONO = 44 
  INTEGER, PARAMETER :: ind_ALK3 = 45 
  INTEGER, PARAMETER :: ind_TBU_O = 46 
  INTEGER, PARAMETER :: ind_ALK5 = 47 
  INTEGER, PARAMETER :: ind_COOH = 48 
  INTEGER, PARAMETER :: ind_HOCOO = 49 
  INTEGER, PARAMETER :: ind_BZNO2_O = 50 
  INTEGER, PARAMETER :: ind_HNO4 = 51 
  INTEGER, PARAMETER :: ind_ARO2 = 52 
  INTEGER, PARAMETER :: ind_MEOH = 53 
  INTEGER, PARAMETER :: ind_ALK4 = 54 
  INTEGER, PARAMETER :: ind_ARO1 = 55 
  INTEGER, PARAMETER :: ind_DCB2 = 56 
  INTEGER, PARAMETER :: ind_DCB3 = 57 
  INTEGER, PARAMETER :: ind_CRES = 58 
  INTEGER, PARAMETER :: ind_DCB1 = 59 
  INTEGER, PARAMETER :: ind_C2H2 = 60 
  INTEGER, PARAMETER :: ind_NPHE = 61 
  INTEGER, PARAMETER :: ind_ROOH = 62 
  INTEGER, PARAMETER :: ind_BALD = 63 
  INTEGER, PARAMETER :: ind_PHEN = 64 
  INTEGER, PARAMETER :: ind_MGLY = 65 
  INTEGER, PARAMETER :: ind_CO = 66 
  INTEGER, PARAMETER :: ind_HNO3 = 67 
  INTEGER, PARAMETER :: ind_ETHENE = 68 
  INTEGER, PARAMETER :: ind_ACET = 69 
  INTEGER, PARAMETER :: ind_C3H6 = 70 
  INTEGER, PARAMETER :: ind_BZ_O = 71 
  INTEGER, PARAMETER :: ind_OLE1 = 72 
  INTEGER, PARAMETER :: ind_ISOPRENE = 73 
  INTEGER, PARAMETER :: ind_R2O2 = 74 
  INTEGER, PARAMETER :: ind_METHACRO = 75 
  INTEGER, PARAMETER :: ind_SESQ = 76 
  INTEGER, PARAMETER :: ind_TERP = 77 
  INTEGER, PARAMETER :: ind_GLY = 78 
  INTEGER, PARAMETER :: ind_ISOPROD = 79 
  INTEGER, PARAMETER :: ind_OLE2 = 80 
  INTEGER, PARAMETER :: ind_MVK = 81 
  INTEGER, PARAMETER :: ind_CCHO = 82 
  INTEGER, PARAMETER :: ind_HCHO = 83 
  INTEGER, PARAMETER :: ind_RNO3 = 84 
  INTEGER, PARAMETER :: ind_O3P = 85 
  INTEGER, PARAMETER :: ind_RCHO = 86 
  INTEGER, PARAMETER :: ind_MEK = 87 
  INTEGER, PARAMETER :: ind_PROD2 = 88 
  INTEGER, PARAMETER :: ind_O3 = 89 
  INTEGER, PARAMETER :: ind_CCO_O2 = 90 
  INTEGER, PARAMETER :: ind_BZCO_O2 = 91 
  INTEGER, PARAMETER :: ind_RCO_O2 = 92 
  INTEGER, PARAMETER :: ind_NO2 = 93 
  INTEGER, PARAMETER :: ind_NO3 = 94 
  INTEGER, PARAMETER :: ind_OH = 95 
  INTEGER, PARAMETER :: ind_NO = 96 
  INTEGER, PARAMETER :: ind_HO2 = 97 
  INTEGER, PARAMETER :: ind_C_O2 = 98 
  INTEGER, PARAMETER :: ind_RO2_R = 99 
  INTEGER, PARAMETER :: ind_MA_RCO3 = 100 
  INTEGER, PARAMETER :: ind_RO2_N = 101 




  INTEGER, PARAMETER :: ind_H2O = 102 
  INTEGER, PARAMETER :: ind_M = 103 




  INTEGER, PARAMETER :: indf_H2O = 1 
  INTEGER, PARAMETER :: indf_M = 2 

END MODULE saprc99_mosaic_20bin_vbs2_aq_Parameters

