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
! File                 : cbm4_JacobianSP.f90
! Time                 : Tue Jan 11 14:32:41 2022
! Working directory    : /jathar-scratch/yicongh/WRF_GoAmazon_gasaerochem.DEV/chem/KPP/mechanisms/cbm4
! Equation file        : cbm4.kpp
! Output root filename : cbm4
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE cbm4_JacobianSP

  PUBLIC
  SAVE


! Sparse Jacobian Data


  INTEGER, PARAMETER, DIMENSION(300) :: LU_IROW = (/ &
       1,  1,  2,  2,  2,  3,  3,  3,  4,  4,  4,  4, &
       4,  5,  5,  6,  6,  6,  7,  7,  8,  8,  8,  8, &
       8,  8,  8,  8,  8,  9,  9,  9,  9, 10, 10, 10, &
      10, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, &
      12, 12, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14, &
      14, 15, 15, 15, 15, 15, 15, 16, 16, 16, 16, 16, &
      16, 16, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17, &
      18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, &
      18, 18, 18, 18, 18, 18, 18, 18, 19, 19, 19, 19, &
      19, 19, 19, 20, 20, 20, 20, 20, 20, 20, 20, 20, &
      20, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, &
      21, 21, 22, 22, 22, 22, 22, 23, 23, 23, 23, 23, &
      23, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, &
      24, 24, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, &
      25, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, &
      26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 27, &
      27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, &
      27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, &
      28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, &
      28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, &
      29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, &
      29, 29, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, &
      30, 30, 30, 30, 30, 31, 31, 31, 31, 31, 31, 31, &
      31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 32, &
      32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32 /)

  INTEGER, PARAMETER, DIMENSION(300) :: LU_ICOL = (/ &
       1, 25,  2, 27, 28,  3, 26, 32,  4, 14, 26, 27, &
      30,  5, 27,  6, 26, 30,  7, 27,  8, 13, 20, 22, &
      23, 27, 29, 30, 31,  9, 26, 27, 31, 10, 26, 27, &
      28,  5,  7, 11, 27, 31,  6, 12, 14, 21, 24, 26, &
      27, 30, 13, 20, 26, 27,  5,  7, 11, 14, 27, 30, &
      31,  7, 15, 19, 22, 25, 27, 15, 16, 17, 19, 21, &
      22, 23, 24, 25, 27, 29, 30, 17, 22, 25, 27, 29, &
       5,  7, 13, 14, 15, 17, 18, 19, 20, 22, 23, 24, &
      25, 26, 27, 28, 29, 30, 31, 32, 11, 14, 19, 25, &
      27, 30, 31,  7, 13, 20, 22, 23, 25, 26, 27, 29, &
      30, 17, 19, 21, 22, 23, 24, 25, 27, 28, 29, 30, &
      31, 32, 22, 25, 27, 29, 30, 22, 23, 25, 27, 29, &
      30, 13, 17, 19, 20, 22, 23, 24, 25, 26, 27, 29, &
      30, 31, 17, 19, 22, 23, 25, 26, 27, 28, 29, 30, &
      31,  3,  4,  6,  9, 10, 11, 13, 14, 18, 19, 20, &
      22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,  1, &
       2,  5,  7,  9, 10, 12, 14, 15, 16, 17, 19, 20, &
      21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, &
       2,  5,  7, 10, 11, 13, 14, 15, 16, 17, 19, 20, &
      21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, &
       1, 17, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, &
      31, 32,  6, 12, 14, 21, 22, 23, 24, 25, 26, 27, &
      28, 29, 30, 31, 32,  8,  9, 11, 13, 18, 19, 20, &
      22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,  3, &
      15, 19, 22, 24, 25, 26, 27, 28, 29, 30, 31, 32 /)

  INTEGER, PARAMETER, DIMENSION(33) :: LU_CROW = (/ &
       1,  3,  6,  9, 14, 16, 19, 21, 30, 34, 38, 43, &
      51, 55, 62, 68, 80, 85,105,112,122,135,140,146, &
     159,170,192,217,241,255,270,288,301 /)

  INTEGER, PARAMETER, DIMENSION(33) :: LU_DIAG = (/ &
       1,  3,  6,  9, 14, 16, 19, 21, 30, 34, 40, 44, &
      51, 58, 63, 69, 80, 91,107,114,124,135,141,152, &
     163,185,211,236,251,267,286,300,301 /)


END MODULE cbm4_JacobianSP

