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
! File                 : saprc99_mosaic_4bin_vbs2_JacobianSP.f90
! Time                 : Tue Jan 11 14:33:24 2022
! Working directory    : /jathar-scratch/yicongh/WRF_GoAmazon_gasaerochem.DEV/chem/KPP/mechanisms/saprc99_mosaic_4bin_vbs2
! Equation file        : saprc99_mosaic_4bin_vbs2.kpp
! Output root filename : saprc99_mosaic_4bin_vbs2
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE saprc99_mosaic_4bin_vbs2_JacobianSP

  PUBLIC
  SAVE


! Sparse Jacobian Data


  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_0 = (/ &
       1,  1,  1,  2,  2,  2,  2,  3,  4,  5,  6,  7, &
       7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7, &
       7,  7,  7,  8,  8,  8,  8,  8,  8,  8,  8,  8, &
       8,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9, &
       9,  9,  9,  9, 10, 10, 10, 11, 11, 11, 11, 11, &
      12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 13, &
      13, 13, 14, 14, 14, 14, 14, 14, 14, 15, 15, 15, &
      15, 15, 15, 16, 16, 16, 16, 16, 16, 17, 17, 17, &
      17, 18, 18, 18, 19, 19, 19, 19, 20, 20, 20, 20, &
      20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, &
      20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, &
      20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, &
      20, 20, 20, 20, 20, 20, 20, 20, 20, 21, 21, 21, &
      21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 22, &
      22, 22, 22, 23, 23, 24, 24, 24, 24, 25, 25, 26, &
      26, 27, 27, 28, 28, 28, 29, 29, 30, 30, 30, 31, &
      31, 32, 32, 33, 33, 34, 34, 34, 35, 35, 35, 36, &
      36, 36, 37, 37, 37, 38, 38, 38, 39, 39, 40, 40, &
      40, 40, 40, 40, 41, 41, 42, 42, 42, 43, 43, 43, &
      44, 44, 45, 45, 45, 45, 46, 46, 47, 47, 47, 47, &
      48, 48, 48, 48, 49, 49, 49, 49, 49, 50, 50, 51, &
      51, 51, 51, 52, 52, 52, 52, 52, 53, 53, 54, 54, &
      55, 55, 55, 55, 56, 56, 56, 56, 57, 57, 57, 57, &
      57, 58, 58, 58, 58, 58, 59, 59, 59, 60, 60, 60, &
      60, 60, 61, 61, 61, 61, 61, 61, 62, 62, 62, 62, &
      62, 62, 62, 63, 63, 63, 63, 63, 63, 64, 64, 64, &
      64, 64, 64, 64, 64, 64, 64, 64, 64, 65, 65, 65, &
      65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, &
      65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 66, 66, &
      66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_1 = (/ &
      66, 66, 66, 66, 66, 66, 66, 67, 67, 67, 67, 67, &
      68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, &
      68, 68, 68, 69, 69, 69, 69, 69, 70, 70, 70, 70, &
      70, 70, 70, 70, 70, 70, 70, 70, 71, 71, 71, 71, &
      71, 72, 72, 72, 72, 72, 73, 73, 73, 73, 73, 73, &
      73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, &
      73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, &
      74, 74, 74, 74, 74, 74, 74, 75, 75, 75, 75, 75, &
      76, 76, 76, 76, 76, 77, 77, 77, 77, 77, 77, 77, &
      77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77, &
      77, 77, 77, 77, 77, 78, 78, 78, 78, 78, 78, 78, &
      79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, &
      79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, &
      80, 80, 80, 80, 80, 81, 81, 81, 81, 81, 81, 81, &
      82, 82, 82, 82, 82, 82, 82, 82, 82, 82, 82, 82, &
      82, 82, 82, 82, 82, 82, 82, 82, 82, 82, 82, 82, &
      82, 82, 82, 82, 82, 82, 82, 82, 82, 82, 82, 82, &
      82, 82, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, &
      83, 83, 83, 83, 84, 84, 84, 84, 84, 84, 84, 84, &
      84, 84, 84, 84, 84, 84, 84, 84, 85, 85, 85, 85, &
      85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, &
      85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, &
      85, 85, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, &
      86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 87, &
      87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, &
      87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 88, &
      88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, &
      88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 89, 89, &
      89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, &
      89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_2 = (/ &
      89, 89, 89, 89, 89, 89, 89, 89, 90, 90, 90, 90, &
      90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, &
      90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, &
      90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, &
      90, 90, 90, 90, 91, 91, 91, 91, 91, 91, 91, 91, &
      91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, &
      91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, &
      92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, &
      92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, &
      92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, &
      92, 92, 92, 92, 92, 92, 93, 93, 93, 93, 93, 93, &
      93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, &
      93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, &
      93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 94, &
      94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, &
      94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, &
      94, 94, 94, 94, 94, 94, 94, 95, 95, 95, 95, 95, &
      95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, &
      96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, &
      96, 96, 96, 96, 96, 96, 96, 96, 97, 97, 97, 97, &
      97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, &
      97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, &
      97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, &
      97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 98, &
      98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, &
      98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, &
      98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, &
      98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, &
      98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, &
      98, 98, 98, 98, 98, 99, 99, 99, 99, 99, 99, 99 /)
  INTEGER, PARAMETER, DIMENSION(43) :: LU_IROW_3 = (/ &
      99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, &
      99, 99, 99, 99, 99, 99,100,100,100,100,100,100, &
     100,100,100,100,100,100,100,100,100,100,100,100, &
     100,100,100,100,100,100,100 /)
  INTEGER, PARAMETER, DIMENSION(1123) :: LU_IROW = (/&
    LU_IROW_0, LU_IROW_1, LU_IROW_2, LU_IROW_3 /)

  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_0 = (/ &
       1, 32, 98,  2, 59, 69, 88,  3,  4,  5,  6,  7, &
      48, 59, 67, 69, 71, 72, 74, 75, 76, 78, 80, 81, &
      88, 98,100,  8, 69, 72, 80, 88, 89, 90, 91, 94, &
      97,  9, 71, 72, 75, 76, 78, 80, 88, 90, 91, 94, &
      95, 96, 97, 99, 10, 89, 97, 11, 95, 96, 97, 99, &
      12, 49, 69, 92, 93, 13, 49, 62, 69, 84, 88, 92, &
      93, 98, 14, 73, 90, 91, 93, 94,100, 15, 73, 90, &
      91, 94, 97, 16, 46, 53, 72, 80, 98, 17, 50, 54, &
      98, 18, 71, 98, 19, 75, 76, 98, 20, 31, 32, 33, &
      39, 41, 43, 44, 46, 47, 50, 52, 53, 54, 55, 56, &
      57, 58, 59, 60, 62, 63, 64, 65, 66, 67, 68, 69, &
      71, 72, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, &
      85, 86, 87, 88, 92, 93, 97, 98,100,  3,  4,  5, &
       6, 21, 22, 23, 24, 25, 27, 28, 29, 30, 98, 22, &
      23, 27, 98, 23, 98, 24, 25, 29, 98, 25, 98, 26, &
      88, 27, 98, 27, 28, 98, 29, 98, 29, 30, 98, 31, &
      98, 32, 98, 33, 98, 34, 89, 92, 35, 92, 99, 36, &
      92, 95, 37, 92, 96, 38, 97, 98, 39, 98, 40, 50, &
      75, 76, 88, 98, 41, 98, 42, 92, 93, 43, 98,100, &
      44, 98, 44, 45, 92, 98, 46, 98, 47, 94, 97, 98, &
      48, 82, 97,100, 49, 61, 92, 93, 97, 50, 98, 51, &
      92, 97, 98, 52, 90, 91, 94, 98, 53, 98, 54, 98, &
      50, 54, 55, 98, 50, 54, 56, 98, 50, 54, 57, 93, &
      98, 50, 54, 58, 88, 98, 59, 88, 98, 60, 90, 91, &
      97, 98, 49, 61, 70, 92, 93, 97, 50, 54, 62, 80, &
      88, 93, 98, 54, 63, 70, 93, 97, 98, 50, 54, 55, &
      56, 57, 64, 74, 78, 81, 88, 93, 98, 53, 55, 56, &
      58, 59, 64, 65, 67, 69, 71, 72, 74, 75, 76, 77, &
      78, 79, 80, 81, 82, 84, 85, 88, 93, 98, 42, 57, &
      61, 62, 63, 64, 66, 70, 74, 77, 78, 79, 80, 81 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_1 = (/ &
      82, 85, 88, 92, 93, 97, 98, 67, 84, 88, 93, 98, &
      39, 44, 45, 46, 53, 68, 72, 75, 76, 80, 83, 88, &
      92, 93, 98, 69, 84, 88, 93, 98, 57, 63, 70, 89, &
      92, 93, 95, 96, 97, 98, 99,100, 71, 84, 88, 93, &
      98, 72, 84, 88, 93, 98, 44, 46, 53, 55, 56, 68, &
      71, 72, 73, 75, 76, 80, 81, 83, 84, 86, 87, 88, &
      89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99,100, &
      71, 74, 80, 84, 88, 93, 98, 75, 84, 88, 93, 98, &
      76, 84, 88, 93, 98, 50, 54, 55, 56, 58, 59, 63, &
      67, 70, 75, 76, 77, 78, 84, 88, 89, 92, 93, 95, &
      96, 97, 98, 99,100, 71, 78, 80, 84, 88, 93, 98, &
      33, 41, 44, 46, 53, 67, 69, 72, 78, 79, 80, 83, &
      84, 85, 86, 87, 88, 89, 93, 95, 96, 98, 99,100, &
      80, 84, 88, 93, 98, 71, 80, 81, 84, 88, 93, 98, &
      41, 44, 46, 47, 48, 52, 53, 59, 67, 68, 69, 71, &
      72, 74, 75, 76, 77, 78, 80, 81, 82, 83, 84, 86, &
      87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, &
      99,100, 45, 72, 75, 76, 78, 80, 83, 84, 88, 91, &
      92, 93, 98,100, 26, 67, 69, 71, 72, 74, 75, 76, &
      80, 81, 84, 88, 92, 93, 98,100, 39, 44, 46, 53, &
      55, 56, 58, 60, 67, 69, 72, 74, 75, 76, 78, 80, &
      81, 83, 84, 85, 86, 87, 88, 90, 91, 92, 93, 97, &
      98,100, 44, 46, 53, 69, 72, 74, 78, 80, 81, 83, &
      84, 86, 87, 88, 90, 91, 92, 93, 94, 98,100, 46, &
      53, 54, 71, 72, 75, 76, 78, 80, 81, 83, 84, 87, &
      88, 89, 90, 91, 92, 93, 94, 95, 98, 99,100, 58, &
      59, 67, 69, 71, 72, 74, 75, 76, 78, 80, 81, 84, &
      88, 89, 92, 93, 95, 96, 97, 98, 99,100, 34, 40, &
      50, 53, 55, 56, 64, 68, 72, 74, 75, 76, 78, 79, &
      80, 81, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_2 = (/ &
      93, 94, 95, 96, 97, 98, 99,100, 33, 39, 41, 44, &
      46, 50, 53, 54, 55, 56, 57, 58, 59, 60, 63, 67, &
      69, 70, 71, 72, 74, 75, 76, 78, 80, 81, 83, 84, &
      85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, &
      97, 98, 99,100, 39, 44, 46, 50, 53, 54, 69, 71, &
      72, 75, 76, 78, 80, 81, 83, 84, 85, 86, 87, 88, &
      89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99,100, &
      34, 35, 36, 37, 42, 43, 45, 48, 49, 51, 61, 66, &
      70, 71, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, &
      83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, &
      95, 96, 97, 98, 99,100, 42, 51, 57, 61, 62, 63, &
      64, 66, 67, 69, 70, 71, 72, 73, 74, 75, 76, 77, &
      78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, &
      90, 91, 92, 93, 94, 95, 96, 97, 98, 99,100, 31, &
      45, 47, 53, 67, 68, 69, 71, 72, 75, 76, 79, 80, &
      81, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, &
      94, 95, 96, 97, 98, 99,100, 36, 62, 80, 84, 88, &
      89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99,100, &
      37, 71, 74, 78, 80, 81, 84, 88, 89, 90, 91, 92, &
      93, 94, 95, 96, 97, 98, 99,100, 32, 38, 41, 43, &
      47, 48, 49, 50, 51, 52, 54, 55, 56, 58, 59, 60, &
      61, 64, 65, 67, 69, 70, 71, 72, 74, 75, 76, 77, &
      78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, &
      90, 91, 92, 93, 94, 95, 96, 97, 98, 99,100,  3, &
       5, 26, 27, 28, 29, 30, 31, 32, 33, 38, 39, 41, &
      43, 44, 46, 47, 50, 51, 52, 53, 54, 55, 56, 57, &
      58, 59, 60, 62, 63, 64, 65, 66, 67, 68, 69, 70, &
      71, 72, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, &
      84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, &
      96, 97, 98, 99,100, 35, 74, 75, 76, 77, 78, 80 /)
  INTEGER, PARAMETER, DIMENSION(43) :: LU_ICOL_3 = (/ &
      81, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, &
      95, 96, 97, 98, 99,100, 43, 48, 73, 75, 76, 80, &
      81, 82, 83, 84, 86, 87, 88, 89, 90, 91, 92, 93, &
      94, 95, 96, 97, 98, 99,100 /)
  INTEGER, PARAMETER, DIMENSION(1123) :: LU_ICOL = (/&
    LU_ICOL_0, LU_ICOL_1, LU_ICOL_2, LU_ICOL_3 /)

  INTEGER, PARAMETER, DIMENSION(101) :: LU_CROW = (/ &
       1,  4,  8,  9, 10, 11, 12, 28, 38, 53, 56, 61, &
      66, 75, 82, 88, 94, 98,101,105,154,168,172,174, &
     178,180,182,184,187,189,192,194,196,198,201,204, &
     207,210,213,215,221,223,226,229,231,235,237,241, &
     245,250,252,256,261,263,265,269,273,278,283,286, &
     291,297,304,310,322,347,368,373,388,393,405,410, &
     415,445,452,457,462,486,493,517,522,529,567,581, &
     597,627,648,672,695,729,773,805,847,888,920,937, &
     957,1008,1074,1099,1124 /)

  INTEGER, PARAMETER, DIMENSION(101) :: LU_DIAG = (/ &
       1,  4,  8,  9, 10, 11, 12, 28, 38, 53, 56, 61, &
      66, 75, 82, 88, 94, 98,101,105,158,168,172,174, &
     178,180,182,185,187,190,192,194,196,198,201,204, &
     207,210,213,215,221,223,226,229,232,235,237,241, &
     245,250,252,256,261,263,267,271,275,280,283,286, &
     292,299,305,315,328,353,368,378,388,395,405,410, &
     423,446,452,457,473,487,502,517,524,549,573,591, &
     616,638,660,685,717,762,795,838,880,913,931,952, &
     1004,1071,1097,1123,1124 /)


END MODULE saprc99_mosaic_4bin_vbs2_JacobianSP

