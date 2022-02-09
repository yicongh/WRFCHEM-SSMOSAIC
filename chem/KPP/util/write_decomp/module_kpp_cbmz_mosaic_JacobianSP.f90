! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Sparse Jacobian Data Structures File
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
! File                 : cbmz_mosaic_JacobianSP.f90
! Time                 : Tue Jan 11 14:32:44 2022
! Working directory    : /jathar-scratch/yicongh/WRF_GoAmazon_gasaerochem.DEV/chem/KPP/mechanisms/cbmz_mosaic
! Equation file        : cbmz_mosaic.kpp
! Output root filename : cbmz_mosaic
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE cbmz_mosaic_JacobianSP

  PUBLIC
  SAVE


! Sparse Jacobian Data


  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_0 = (/ &
       1,  1,  1,  2,  3,  4,  4,  4,  4,  5,  5,  5, &
       5,  5,  5,  5,  5,  6,  6,  6,  6,  6,  6,  7, &
       7,  7,  7,  7,  7,  8,  8,  8,  9,  9,  9,  9, &
       9,  9, 10, 10, 10, 10, 11, 11, 11, 11, 12, 12, &
      12, 12, 13, 13, 13, 13, 14, 14, 15, 15, 16, 16, &
      17, 17, 17, 18, 18, 18, 19, 19, 20, 20, 20, 21, &
      21, 22, 22, 22, 22, 22, 23, 23, 23, 23, 23, 24, &
      24, 24, 24, 25, 25, 25, 25, 26, 26, 26, 26, 27, &
      27, 27, 27, 27, 28, 28, 28, 28, 28, 28, 29, 29, &
      29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 30, 30, &
      30, 30, 31, 31, 31, 31, 31, 32, 32, 32, 32, 33, &
      33, 33, 34, 34, 34, 34, 34, 34, 35, 35, 35, 35, &
      35, 36, 36, 36, 36, 36, 36, 37, 37, 37, 37, 37, &
      37, 37, 37, 37, 37, 37, 38, 38, 38, 38, 38, 38, &
      38, 38, 38, 38, 38, 38, 38, 38, 38, 39, 39, 39, &
      39, 39, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, &
      40, 40, 40, 40, 40, 40, 40, 40, 41, 41, 41, 41, &
      41, 41, 41, 42, 42, 42, 42, 42, 43, 43, 43, 43, &
      43, 44, 44, 44, 44, 45, 45, 45, 45, 46, 46, 46, &
      46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, &
      46, 46, 46, 46, 47, 47, 47, 47, 47, 47, 47, 47, &
      47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 48, &
      48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, &
      49, 49, 49, 49, 50, 50, 50, 50, 50, 50, 50, 50, &
      50, 50, 50, 50, 50, 50, 51, 51, 51, 51, 51, 51, &
      51, 51, 51, 51, 51, 51, 51, 52, 52, 52, 52, 52, &
      52, 52, 52, 52, 53, 53, 53, 53, 53, 53, 53, 53, &
      53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, &
      53, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, &
      54, 54, 54, 54, 54, 54, 54, 55, 55, 55, 55, 55 /)
  INTEGER, PARAMETER, DIMENSION(266) :: LU_IROW_1 = (/ &
      55, 55, 55, 55, 56, 56, 56, 56, 56, 56, 56, 56, &
      56, 56, 56, 56, 57, 57, 57, 57, 57, 57, 57, 57, &
      57, 57, 57, 57, 58, 58, 58, 58, 58, 58, 58, 58, &
      58, 58, 58, 58, 58, 58, 58, 58, 59, 59, 59, 59, &
      59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, &
      60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, &
      60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, &
      60, 60, 60, 60, 61, 61, 61, 61, 61, 61, 61, 61, &
      61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, &
      61, 61, 61, 61, 61, 62, 62, 62, 62, 62, 62, 62, &
      62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, &
      62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, &
      62, 62, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, &
      63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, &
      63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, &
      63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 64, &
      64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, &
      64, 64, 64, 65, 65, 65, 65, 65, 65, 65, 65, 65, &
      65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, &
      65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, &
      65, 65, 65, 65, 65, 65, 66, 66, 66, 66, 66, 66, &
      66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, &
      66, 66 /)
  INTEGER, PARAMETER, DIMENSION(626) :: LU_IROW = (/&
    LU_IROW_0, LU_IROW_1 /)

  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_0 = (/ &
       1, 14, 63,  2,  3,  4, 33, 44, 64,  5, 44, 45, &
      49, 55, 64, 65, 66,  6, 19, 21, 35, 60, 63,  7, &
      19, 21, 35, 60, 63,  8, 40, 63,  9, 44, 49, 60, &
      63, 64, 10, 24, 63, 64, 11, 24, 63, 64, 12, 25, &
      63, 64, 13, 25, 63, 64, 14, 63, 15, 64, 16, 63, &
      17, 63, 65, 18, 62, 66, 19, 63, 20, 60, 62, 21, &
      63, 22, 44, 49, 63, 64, 23, 35, 60, 62, 63, 24, &
      60, 63, 64, 25, 60, 63, 64, 26, 62, 63, 65, 19, &
      21, 27, 61, 63, 28, 44, 49, 51, 63, 64, 29, 40, &
      44, 49, 52, 56, 58, 59, 60, 61, 63, 64, 30, 51, &
      63, 65, 31, 61, 62, 63, 65, 32, 54, 63, 65, 33, &
      63, 64, 34, 44, 49, 54, 63, 64, 19, 21, 35, 60, &
      63, 15, 36, 60, 61, 62, 64, 20, 35, 37, 46, 50, &
      53, 55, 60, 62, 63, 65, 33, 38, 41, 43, 44, 45, &
      46, 49, 50, 53, 55, 60, 61, 63, 64, 39, 45, 60, &
      61, 65, 21, 29, 39, 40, 42, 44, 45, 49, 52, 55, &
      56, 58, 59, 60, 61, 63, 64, 65, 27, 35, 41, 60, &
      61, 63, 64, 42, 45, 61, 63, 65, 43, 55, 61, 63, &
      65, 44, 60, 63, 64, 45, 60, 63, 64, 32, 33, 34, &
      41, 42, 43, 44, 45, 46, 49, 52, 54, 55, 57, 60, &
      61, 63, 64, 65, 19, 21, 33, 35, 41, 44, 45, 47, &
      49, 50, 55, 56, 58, 59, 60, 61, 63, 64, 65, 43, &
      48, 49, 55, 56, 57, 58, 59, 60, 61, 63, 64, 65, &
      49, 60, 63, 64, 21, 41, 43, 44, 49, 50, 55, 56, &
      57, 60, 61, 63, 64, 65, 28, 30, 44, 49, 51, 56, &
      58, 59, 60, 61, 63, 64, 65, 44, 49, 52, 59, 60, &
      61, 63, 64, 65, 16, 30, 33, 39, 41, 43, 44, 45, &
      49, 51, 52, 53, 55, 56, 58, 59, 60, 61, 63, 64, &
      65, 22, 32, 44, 48, 49, 53, 54, 55, 56, 57, 58, &
      59, 60, 61, 63, 64, 65, 66, 39, 42, 45, 55, 60 /)
  INTEGER, PARAMETER, DIMENSION(266) :: LU_ICOL_1 = (/ &
      61, 63, 64, 65, 42, 43, 45, 55, 56, 57, 58, 60, &
      61, 63, 64, 65, 48, 49, 55, 56, 57, 58, 59, 60, &
      61, 63, 64, 65, 40, 42, 44, 45, 49, 52, 55, 56, &
      57, 58, 59, 60, 61, 63, 64, 65, 23, 27, 35, 39, &
      42, 45, 52, 55, 58, 59, 60, 61, 62, 63, 64, 65, &
      20, 24, 25, 35, 36, 37, 44, 45, 46, 47, 49, 50, &
      51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, &
      63, 64, 65, 66, 27, 31, 36, 39, 42, 43, 45, 47, &
      49, 50, 51, 52, 54, 55, 56, 57, 58, 59, 60, 61, &
      62, 63, 64, 65, 66, 18, 20, 23, 26, 27, 31, 35, &
      36, 37, 39, 42, 43, 45, 46, 47, 49, 50, 51, 52, &
      53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, &
      65, 66, 14, 15, 16, 17, 19, 21, 22, 24, 25, 26, &
      28, 30, 31, 32, 33, 34, 35, 37, 38, 40, 41, 42, &
      43, 44, 45, 46, 48, 49, 50, 51, 52, 53, 54, 55, &
      56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 24, &
      25, 33, 36, 41, 44, 45, 49, 55, 60, 61, 62, 63, &
      64, 65, 66, 14, 16, 17, 19, 21, 26, 27, 30, 32, &
      33, 34, 35, 38, 39, 41, 42, 43, 44, 45, 46, 47, &
      49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, &
      61, 62, 63, 64, 65, 66, 18, 41, 44, 45, 48, 49, &
      50, 53, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, &
      65, 66 /)
  INTEGER, PARAMETER, DIMENSION(626) :: LU_ICOL = (/&
    LU_ICOL_0, LU_ICOL_1 /)

  INTEGER, PARAMETER, DIMENSION(67) :: LU_CROW = (/ &
       1,  4,  5,  6, 10, 18, 24, 30, 33, 39, 43, 47, &
      51, 55, 57, 59, 61, 64, 67, 69, 72, 74, 79, 84, &
      88, 92, 96,101,107,119,123,128,132,135,141,146, &
     152,163,178,183,201,208,213,218,222,226,245,264, &
     277,281,295,308,317,338,356,365,377,389,405,421, &
     449,474,507,552,568,607,627 /)

  INTEGER, PARAMETER, DIMENSION(67) :: LU_DIAG = (/ &
       1,  4,  5,  6, 10, 18, 24, 30, 33, 39, 43, 47, &
      51, 55, 57, 59, 61, 64, 67, 69, 72, 74, 79, 84, &
      88, 92, 98,101,107,119,123,128,132,135,143,147, &
     154,164,178,186,203,208,213,218,222,234,252,265, &
     277,286,299,310,328,344,359,369,381,398,414,442, &
     468,502,548,565,605,626,627 /)


END MODULE cbmz_mosaic_JacobianSP

