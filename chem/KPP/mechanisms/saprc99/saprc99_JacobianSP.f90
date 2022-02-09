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
! File                 : saprc99_JacobianSP.f90
! Time                 : Tue Jan 11 14:33:29 2022
! Working directory    : /jathar-scratch/yicongh/WRF_GoAmazon_gasaerochem.DEV/chem/KPP/mechanisms/saprc99
! Equation file        : saprc99.kpp
! Output root filename : saprc99
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE saprc99_JacobianSP

  PUBLIC
  SAVE


! Sparse Jacobian Data


  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_0 = (/ &
       1,  1,  1,  2,  2,  2,  2,  2,  2,  2,  2,  2, &
       2,  2,  2,  2,  2,  2,  2,  3,  3,  3,  3,  3, &
       3,  3,  3,  3,  3,  4,  4,  4,  4,  4,  4,  4, &
       4,  4,  4,  4,  4,  4,  4,  4,  5,  5,  5,  6, &
       6,  6,  6,  6,  7,  7,  7,  8,  8,  8,  8,  8, &
       9,  9,  9,  9, 10, 10, 10, 11, 11, 11, 11, 11, &
      12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 13, &
      13, 13, 14, 14, 15, 15, 16, 16, 17, 17, 17, 17, &
      18, 18, 18, 18, 18, 18, 19, 19, 20, 20, 20, 21, &
      21, 21, 22, 22, 22, 23, 23, 23, 24, 24, 24, 25, &
      25, 26, 26, 27, 27, 27, 27, 27, 27, 28, 28, 28, &
      29, 29, 29, 30, 30, 31, 31, 31, 31, 32, 32, 33, &
      33, 33, 33, 33, 34, 34, 35, 35, 35, 35, 36, 36, &
      36, 36, 37, 37, 37, 37, 37, 38, 38, 38, 38, 39, &
      39, 40, 40, 41, 41, 41, 41, 42, 42, 42, 42, 43, &
      43, 43, 43, 43, 44, 44, 44, 45, 45, 45, 45, 45, &
      46, 46, 46, 46, 46, 46, 47, 47, 47, 47, 47, 47, &
      47, 48, 48, 48, 48, 48, 49, 49, 49, 49, 49, 49, &
      50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, &
      51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, &
      51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, &
      51, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, &
      52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 53, 53, &
      53, 53, 53, 54, 54, 54, 54, 54, 54, 54, 54, 54, &
      54, 54, 54, 54, 54, 54, 55, 55, 55, 55, 55, 56, &
      56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, &
      56, 56, 56, 56, 56, 57, 57, 57, 57, 57, 57, 57, &
      57, 57, 57, 57, 57, 58, 58, 58, 58, 58, 59, 59, &
      59, 59, 59, 60, 60, 60, 60, 60, 60, 60, 60, 60, &
      60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_1 = (/ &
      60, 60, 60, 60, 60, 60, 60, 60, 60, 61, 61, 61, &
      61, 61, 62, 62, 62, 62, 62, 62, 62, 63, 63, 63, &
      63, 63, 64, 64, 64, 64, 64, 64, 64, 65, 65, 65, &
      65, 65, 66, 66, 66, 66, 66, 66, 66, 67, 67, 67, &
      67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, &
      67, 67, 67, 67, 67, 67, 67, 67, 67, 68, 68, 68, &
      68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, &
      68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, &
      68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, &
      69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, &
      69, 69, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, &
      70, 70, 70, 70, 70, 70, 71, 71, 71, 71, 71, 71, &
      71, 71, 71, 71, 71, 71, 71, 71, 71, 71, 71, 71, &
      71, 71, 71, 71, 71, 71, 71, 71, 71, 71, 71, 71, &
      72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, &
      72, 72, 72, 72, 72, 72, 72, 72, 72, 73, 73, 73, &
      73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, &
      73, 73, 73, 73, 73, 73, 73, 73, 73, 74, 74, 74, &
      74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, &
      74, 74, 74, 74, 74, 74, 74, 74, 75, 75, 75, 75, &
      75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, &
      75, 75, 75, 75, 76, 76, 76, 76, 76, 76, 76, 76, &
      76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, &
      76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, &
      76, 76, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77, &
      77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77, &
      77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77, &
      77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 78, 78, &
      78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, &
      78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 79 /)
  INTEGER, PARAMETER, DIMENSION(301) :: LU_IROW_2 = (/ &
      79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, &
      79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, &
      79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, &
      79, 79, 79, 79, 79, 79, 80, 80, 80, 80, 80, 80, &
      80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, &
      80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, &
      80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 81, &
      81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, &
      81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, &
      81, 81, 81, 81, 81, 81, 81, 82, 82, 82, 82, 82, &
      82, 82, 82, 82, 82, 82, 82, 82, 82, 82, 82, 82, &
      82, 82, 82, 82, 82, 82, 82, 82, 82, 82, 82, 82, &
      82, 82, 82, 82, 82, 82, 82, 82, 82, 82, 82, 82, &
      82, 82, 82, 82, 82, 82, 82, 82, 82, 82, 82, 82, &
      82, 82, 82, 82, 82, 82, 82, 83, 83, 83, 83, 83, &
      83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, &
      83, 83, 83, 83, 83, 83, 83, 83, 83, 84, 84, 84, &
      84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, &
      84, 84, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, &
      85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, &
      85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, &
      85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, &
      85, 85, 85, 85, 85, 86, 86, 86, 86, 86, 86, 86, &
      86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, &
      86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, &
      86 /)
  INTEGER, PARAMETER, DIMENSION(1021) :: LU_IROW = (/&
    LU_IROW_0, LU_IROW_1, LU_IROW_2 /)

  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_0 = (/ &
       1, 16, 82,  2, 36, 44, 53, 55, 58, 59, 61, 62, &
      63, 64, 65, 66, 74, 78, 82,  3, 55, 63, 65, 74, &
      76, 77, 81, 85, 86,  4, 58, 59, 61, 63, 64, 65, &
      74, 75, 77, 81, 83, 84, 85, 86,  5,  8, 85,  6, &
       7,  8, 18, 82,  7, 17, 82,  8, 17, 78, 82, 85, &
       9, 44, 55, 74, 10, 76, 85, 11, 75, 83, 84, 85, &
      12, 37, 55, 79, 80, 13, 37, 47, 55, 70, 74, 79, &
      80, 82, 14, 74, 15, 82, 16, 82, 17, 18, 82, 85, &
      17, 18, 58, 78, 82, 85, 19, 82, 20, 76, 79, 21, &
      79, 83, 22, 79, 84, 23, 75, 79, 24, 82, 85, 25, &
      82, 26, 82, 27, 34, 59, 61, 74, 82, 28, 79, 80, &
      29, 78, 82, 30, 82, 30, 31, 79, 82, 32, 82, 33, &
      77, 81, 82, 86, 34, 82, 35, 81, 82, 85, 36, 68, &
      78, 85, 37, 46, 79, 80, 85, 38, 79, 82, 85, 39, &
      82, 40, 82, 34, 40, 41, 82, 34, 40, 42, 82, 34, &
      40, 43, 80, 82, 44, 74, 82, 34, 40, 45, 74, 82, &
      37, 46, 57, 79, 80, 85, 34, 40, 47, 65, 74, 80, &
      82, 48, 77, 82, 85, 86, 40, 49, 57, 80, 82, 85, &
      34, 40, 41, 42, 43, 50, 62, 64, 66, 74, 80, 82, &
      39, 41, 42, 44, 45, 50, 51, 53, 55, 56, 58, 59, &
      61, 62, 63, 64, 65, 66, 67, 68, 70, 71, 74, 80, &
      82, 28, 43, 46, 47, 49, 50, 52, 56, 57, 62, 64, &
      65, 66, 67, 68, 71, 74, 79, 80, 82, 85, 53, 70, &
      74, 80, 82, 26, 30, 31, 32, 39, 54, 59, 61, 63, &
      65, 69, 74, 79, 80, 82, 55, 70, 74, 80, 82, 34, &
      40, 41, 42, 44, 45, 49, 53, 56, 57, 59, 61, 64, &
      70, 74, 80, 82, 85, 43, 49, 57, 75, 76, 78, 79, &
      80, 82, 83, 84, 85, 58, 70, 74, 80, 82, 59, 70, &
      74, 80, 82, 30, 32, 39, 41, 42, 54, 58, 59, 60, &
      61, 63, 65, 66, 69, 70, 72, 73, 74, 75, 76, 77 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_1 = (/ &
      78, 79, 80, 81, 82, 83, 84, 85, 86, 61, 70, 74, &
      80, 82, 58, 62, 65, 70, 74, 80, 82, 63, 70, 74, &
      80, 82, 58, 64, 65, 70, 74, 80, 82, 65, 70, 74, &
      80, 82, 58, 65, 66, 70, 74, 80, 82, 19, 25, 30, &
      32, 39, 53, 55, 63, 64, 65, 67, 69, 70, 71, 72, &
      73, 74, 75, 76, 78, 80, 82, 83, 84, 25, 30, 32, &
      33, 35, 36, 39, 44, 53, 54, 55, 56, 57, 58, 59, &
      61, 62, 63, 64, 65, 66, 68, 69, 70, 72, 73, 74, &
      75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, &
      31, 59, 61, 63, 64, 65, 69, 70, 74, 78, 79, 80, &
      82, 86, 14, 53, 55, 58, 59, 61, 62, 63, 65, 66, &
      70, 74, 78, 79, 80, 82, 26, 30, 32, 39, 41, 42, &
      45, 48, 53, 55, 59, 61, 62, 63, 64, 65, 66, 69, &
      70, 71, 72, 73, 74, 77, 78, 79, 80, 82, 85, 86, &
      30, 32, 39, 55, 62, 63, 64, 65, 66, 69, 70, 72, &
      73, 74, 77, 78, 79, 80, 81, 82, 86, 32, 39, 40, &
      58, 59, 61, 63, 64, 65, 66, 69, 70, 73, 74, 76, &
      77, 78, 79, 80, 81, 82, 83, 84, 86, 44, 45, 53, &
      55, 58, 59, 61, 62, 63, 64, 65, 66, 70, 74, 75, &
      76, 78, 79, 80, 82, 83, 84, 85, 23, 58, 62, 64, &
      65, 66, 70, 74, 75, 76, 77, 78, 79, 80, 81, 82, &
      83, 84, 85, 86, 20, 27, 34, 39, 41, 42, 50, 54, &
      59, 61, 62, 63, 64, 65, 66, 67, 69, 70, 71, 72, &
      73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, &
      85, 86, 19, 25, 26, 30, 32, 34, 39, 40, 41, 42, &
      43, 44, 45, 48, 49, 53, 55, 57, 58, 59, 61, 62, &
      63, 64, 65, 66, 69, 70, 71, 72, 73, 74, 75, 76, &
      77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 29, 36, &
      60, 61, 63, 65, 66, 68, 69, 70, 72, 73, 74, 75, &
      76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 20 /)
  INTEGER, PARAMETER, DIMENSION(301) :: LU_ICOL_2 = (/ &
      21, 22, 23, 28, 29, 31, 36, 37, 38, 46, 52, 56, &
      57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, &
      69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, &
      81, 82, 83, 84, 85, 86, 28, 38, 43, 46, 47, 49, &
      50, 52, 53, 55, 56, 57, 58, 59, 60, 61, 62, 63, &
      64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, &
      76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 15, &
      31, 35, 39, 53, 54, 55, 58, 59, 61, 63, 65, 66, &
      67, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, &
      80, 81, 82, 83, 84, 85, 86, 14, 15, 16, 19, 24, &
      25, 26, 29, 30, 32, 33, 34, 35, 38, 39, 40, 41, &
      42, 43, 44, 45, 47, 48, 49, 50, 51, 52, 53, 54, &
      55, 56, 57, 58, 59, 61, 62, 63, 64, 65, 66, 67, &
      68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, &
      80, 81, 82, 83, 84, 85, 86, 21, 56, 57, 59, 61, &
      62, 64, 65, 66, 70, 71, 72, 73, 74, 75, 76, 77, &
      78, 79, 80, 81, 82, 83, 84, 85, 86, 22, 47, 65, &
      70, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, &
      85, 86, 16, 24, 25, 29, 33, 34, 35, 36, 37, 38, &
      40, 41, 42, 44, 45, 46, 48, 50, 51, 53, 55, 56, &
      57, 58, 59, 61, 62, 63, 64, 65, 66, 67, 68, 69, &
      70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, &
      82, 83, 84, 85, 86, 26, 30, 32, 34, 39, 40, 55, &
      58, 59, 61, 63, 64, 65, 66, 69, 70, 71, 72, 73, &
      74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, &
      86 /)
  INTEGER, PARAMETER, DIMENSION(1021) :: LU_ICOL = (/&
    LU_ICOL_0, LU_ICOL_1, LU_ICOL_2 /)

  INTEGER, PARAMETER, DIMENSION(87) :: LU_CROW = (/ &
       1,  4, 20, 30, 45, 48, 53, 56, 61, 65, 68, 73, &
      78, 87, 89, 91, 93, 97,103,105,108,111,114,117, &
     120,122,124,130,133,136,138,142,144,149,151,155, &
     159,164,168,170,172,176,180,185,188,193,199,206, &
     211,217,229,254,275,280,295,300,318,330,335,340, &
     370,375,382,387,394,399,406,430,469,483,499,529, &
     550,574,597,617,651,695,720,763,804,836,896,922, &
     939,990,1022 /)

  INTEGER, PARAMETER, DIMENSION(87) :: LU_DIAG = (/ &
       1,  4, 20, 30, 45, 48, 53, 56, 61, 65, 68, 73, &
      78, 87, 89, 91, 93, 98,103,105,108,111,114,117, &
     120,122,124,130,133,136,139,142,144,149,151,155, &
     159,164,168,170,174,178,182,185,190,194,201,206, &
     212,222,235,260,275,285,295,308,320,330,335,348, &
     370,376,382,388,394,401,416,451,475,493,518,540, &
     562,587,605,640,685,711,755,797,830,891,918,936, &
     988,1021,1022 /)


END MODULE saprc99_JacobianSP

