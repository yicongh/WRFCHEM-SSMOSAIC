























MODULE saprc99_simplesom_mosaic_4bin_aq_Parameters

  USE saprc99_simplesom_mosaic_4bin_aq_Precision
  PUBLIC
  SAVE



  INTEGER, PARAMETER :: NSPEC = 106 

  INTEGER, PARAMETER :: NVAR = 104 

  INTEGER, PARAMETER :: NVARACT = 90 

  INTEGER, PARAMETER :: NFIX = 2 

  INTEGER, PARAMETER :: NREACT = 247 

  INTEGER, PARAMETER :: NVARST = 1 

  INTEGER, PARAMETER :: NFIXST = 105 

  INTEGER, PARAMETER :: NONZERO = 1078 

  INTEGER, PARAMETER :: LU_NONZERO = 1176 

  INTEGER, PARAMETER :: CNVAR = 105 

  INTEGER, PARAMETER :: NLOOKAT = 0 

  INTEGER, PARAMETER :: NMONITOR = 0 

  INTEGER, PARAMETER :: NMASS = 1 

  REAL(kind=dp), PARAMETER :: PI = 3.14159265358979 




  INTEGER, PARAMETER :: ind_H2SO4 = 1 
  INTEGER, PARAMETER :: ind_HCOOH = 2 
  INTEGER, PARAMETER :: ind_CCO_OH = 3 
  INTEGER, PARAMETER :: ind_RCO_OH = 4 
  INTEGER, PARAMETER :: ind_CO2 = 5 
  INTEGER, PARAMETER :: ind_CCO_OOH = 6 
  INTEGER, PARAMETER :: ind_RCO_OOH = 7 
  INTEGER, PARAMETER :: ind_XN = 8 
  INTEGER, PARAMETER :: ind_XC = 9 
  INTEGER, PARAMETER :: ind_NUME = 10 
  INTEGER, PARAMETER :: ind_DEN = 11 
  INTEGER, PARAMETER :: ind_ANT2_o = 12 
  INTEGER, PARAMETER :: ind_ANT3_o = 13 
  INTEGER, PARAMETER :: ind_PSD1 = 14 
  INTEGER, PARAMETER :: ind_ISOP2O2 = 15 
  INTEGER, PARAMETER :: ind_IEPOX = 16 
  INTEGER, PARAMETER :: ind_TETROL = 17 
  INTEGER, PARAMETER :: ind_O1D = 18 
  INTEGER, PARAMETER :: ind_CH4 = 19 
  INTEGER, PARAMETER :: ind_ISOPOOH = 20 
  INTEGER, PARAMETER :: ind_ISOPOO = 21 
  INTEGER, PARAMETER :: ind_C2H6 = 22 
  INTEGER, PARAMETER :: ind_PAN = 23 
  INTEGER, PARAMETER :: ind_PAN2 = 24 
  INTEGER, PARAMETER :: ind_PBZN = 25 
  INTEGER, PARAMETER :: ind_MA_PAN = 26 
  INTEGER, PARAMETER :: ind_H2O2 = 27 
  INTEGER, PARAMETER :: ind_BACL = 28 
  INTEGER, PARAMETER :: ind_C3H8 = 29 
  INTEGER, PARAMETER :: ind_ETOH = 30 
  INTEGER, PARAMETER :: ind_N2O5 = 31 
  INTEGER, PARAMETER :: ind_SO2 = 32 
  INTEGER, PARAMETER :: ind_HONO = 33 
  INTEGER, PARAMETER :: ind_DMS = 34 
  INTEGER, PARAMETER :: ind_ALK3 = 35 
  INTEGER, PARAMETER :: ind_TBU_O = 36 
  INTEGER, PARAMETER :: ind_ALK5 = 37 
  INTEGER, PARAMETER :: ind_ARO2 = 38 
  INTEGER, PARAMETER :: ind_MEOH = 39 
  INTEGER, PARAMETER :: ind_COOH = 40 
  INTEGER, PARAMETER :: ind_HOCOO = 41 
  INTEGER, PARAMETER :: ind_BZNO2_O = 42 
  INTEGER, PARAMETER :: ind_HNO4 = 43 
  INTEGER, PARAMETER :: ind_ALK4 = 44 
  INTEGER, PARAMETER :: ind_ARO1 = 45 
  INTEGER, PARAMETER :: ind_DCB2 = 46 
  INTEGER, PARAMETER :: ind_DCB3 = 47 
  INTEGER, PARAMETER :: ind_CRES = 48 
  INTEGER, PARAMETER :: ind_DCB1 = 49 
  INTEGER, PARAMETER :: ind_C2H2 = 50 
  INTEGER, PARAMETER :: ind_NPHE = 51 
  INTEGER, PARAMETER :: ind_BALD = 52 
  INTEGER, PARAMETER :: ind_ROOH = 53 
  INTEGER, PARAMETER :: ind_TERP_CN2_g = 54 
  INTEGER, PARAMETER :: ind_PHEN = 55 
  INTEGER, PARAMETER :: ind_TERP_CN1_g = 56 
  INTEGER, PARAMETER :: ind_TERP_CN3_g = 57 
  INTEGER, PARAMETER :: ind_TERP_C00_g = 58 
  INTEGER, PARAMETER :: ind_TERP_C01_g = 59 
  INTEGER, PARAMETER :: ind_TERP_C02_g = 60 
  INTEGER, PARAMETER :: ind_TERP_C03_g = 61 
  INTEGER, PARAMETER :: ind_TERP_C04_g = 62 
  INTEGER, PARAMETER :: ind_TERP_C05_g = 63 
  INTEGER, PARAMETER :: ind_TERP_C09_g = 64 
  INTEGER, PARAMETER :: ind_TERP_C06_g = 65 
  INTEGER, PARAMETER :: ind_TERP_C07_g = 66 
  INTEGER, PARAMETER :: ind_TERP_C08_g = 67 
  INTEGER, PARAMETER :: ind_MGLY = 68 
  INTEGER, PARAMETER :: ind_CO = 69 
  INTEGER, PARAMETER :: ind_HNO3 = 70 
  INTEGER, PARAMETER :: ind_ETHENE = 71 
  INTEGER, PARAMETER :: ind_ACET = 72 
  INTEGER, PARAMETER :: ind_C3H6 = 73 
  INTEGER, PARAMETER :: ind_TERP = 74 
  INTEGER, PARAMETER :: ind_GLY = 75 
  INTEGER, PARAMETER :: ind_BZ_O = 76 
  INTEGER, PARAMETER :: ind_OLE1 = 77 
  INTEGER, PARAMETER :: ind_R2O2 = 78 
  INTEGER, PARAMETER :: ind_SESQ = 79 
  INTEGER, PARAMETER :: ind_ISOPRENE = 80 
  INTEGER, PARAMETER :: ind_METHACRO = 81 
  INTEGER, PARAMETER :: ind_ISOPROD = 82 
  INTEGER, PARAMETER :: ind_OLE2 = 83 
  INTEGER, PARAMETER :: ind_MVK = 84 
  INTEGER, PARAMETER :: ind_CCHO = 85 
  INTEGER, PARAMETER :: ind_HCHO = 86 
  INTEGER, PARAMETER :: ind_RNO3 = 87 
  INTEGER, PARAMETER :: ind_O3P = 88 
  INTEGER, PARAMETER :: ind_RCHO = 89 
  INTEGER, PARAMETER :: ind_MEK = 90 
  INTEGER, PARAMETER :: ind_PROD2 = 91 
  INTEGER, PARAMETER :: ind_O3 = 92 
  INTEGER, PARAMETER :: ind_NO = 93 
  INTEGER, PARAMETER :: ind_NO2 = 94 
  INTEGER, PARAMETER :: ind_NO3 = 95 
  INTEGER, PARAMETER :: ind_BZCO_O2 = 96 
  INTEGER, PARAMETER :: ind_RCO_O2 = 97 
  INTEGER, PARAMETER :: ind_MA_RCO3 = 98 
  INTEGER, PARAMETER :: ind_C_O2 = 99 
  INTEGER, PARAMETER :: ind_CCO_O2 = 100 
  INTEGER, PARAMETER :: ind_RO2_N = 101 
  INTEGER, PARAMETER :: ind_OH = 102 
  INTEGER, PARAMETER :: ind_HO2 = 103 
  INTEGER, PARAMETER :: ind_RO2_R = 104 




  INTEGER, PARAMETER :: ind_H2O = 105 
  INTEGER, PARAMETER :: ind_M = 106 




  INTEGER, PARAMETER :: indf_H2O = 1 
  INTEGER, PARAMETER :: indf_M = 2 

END MODULE saprc99_simplesom_mosaic_4bin_aq_Parameters
