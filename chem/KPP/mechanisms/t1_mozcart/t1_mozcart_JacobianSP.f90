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
! File                 : t1_mozcart_JacobianSP.f90
! Time                 : Tue Jan 11 14:33:32 2022
! Working directory    : /jathar-scratch/yicongh/WRF_GoAmazon_gasaerochem.DEV/chem/KPP/mechanisms/t1_mozcart
! Equation file        : t1_mozcart.kpp
! Output root filename : t1_mozcart
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE t1_mozcart_JacobianSP

  PUBLIC
  SAVE


! Sparse Jacobian Data


  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_0 = (/ &
       1,  1,  1,  2,  2,  3,  3,  4,  4,  5,  5,  6, &
       6,  7,  7,  7,  7,  8,  8,  8,  9,  9,  9, 10, &
      10, 10, 11, 11, 12, 12, 13, 13, 13, 14, 14, 14, &
      15, 15, 15, 16, 16, 16, 17, 17, 18, 18, 18, 19, &
      19, 19, 19, 20, 20, 20, 20, 21, 21, 22, 22, 22, &
      23, 23, 23, 24, 24, 24, 25, 25, 25, 25, 25, 26, &
      26, 26, 26, 26, 27, 27, 27, 27, 28, 28, 28, 28, &
      29, 29, 29, 29, 30, 30, 30, 30, 30, 31, 31, 31, &
      31, 31, 31, 32, 32, 32, 32, 33, 33, 33, 33, 34, &
      34, 34, 34, 35, 35, 35, 35, 35, 35, 35, 36, 36, &
      36, 36, 37, 37, 37, 37, 38, 38, 38, 38, 39, 39, &
      39, 39, 40, 40, 40, 40, 41, 41, 41, 41, 42, 42, &
      42, 42, 43, 43, 43, 43, 43, 43, 43, 43, 44, 44, &
      44, 45, 45, 45, 45, 45, 45, 45, 45, 46, 46, 46, &
      47, 47, 47, 47, 48, 48, 48, 48, 49, 49, 49, 49, &
      50, 50, 50, 50, 51, 51, 51, 51, 52, 52, 52, 52, &
      52, 52, 52, 53, 53, 53, 53, 54, 54, 54, 54, 54, &
      54, 54, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, &
      55, 56, 56, 56, 56, 56, 57, 57, 57, 57, 58, 58, &
      58, 58, 58, 58, 58, 58, 58, 58, 58, 59, 59, 59, &
      59, 60, 60, 60, 60, 61, 61, 61, 61, 61, 61, 62, &
      62, 62, 62, 62, 62, 63, 63, 63, 63, 63, 63, 64, &
      64, 64, 64, 64, 65, 65, 65, 65, 65, 65, 65, 65, &
      65, 66, 66, 66, 66, 66, 66, 67, 67, 67, 67, 68, &
      68, 68, 68, 69, 69, 69, 69, 70, 70, 70, 70, 70, &
      71, 71, 71, 71, 71, 72, 72, 72, 72, 73, 73, 73, &
      73, 73, 73, 73, 73, 73, 73, 74, 74, 74, 74, 75, &
      75, 75, 75, 75, 75, 76, 76, 76, 76, 76, 76, 77, &
      77, 77, 77, 77, 77, 77, 77, 77, 78, 78, 78, 78, &
      79, 79, 79, 79, 80, 80, 80, 80, 81, 81, 81, 81 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_1 = (/ &
      81, 81, 81, 81, 81, 82, 82, 82, 82, 82, 82, 83, &
      83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, &
      83, 83, 83, 83, 83, 83, 84, 84, 84, 84, 84, 84, &
      84, 85, 85, 85, 85, 86, 86, 86, 86, 86, 86, 86, &
      86, 86, 86, 86, 87, 87, 87, 87, 87, 87, 88, 88, &
      88, 88, 88, 88, 88, 89, 89, 89, 89, 90, 90, 90, &
      90, 90, 90, 91, 91, 91, 91, 91, 91, 91, 91, 91, &
      91, 92, 92, 92, 92, 92, 92, 92, 93, 93, 93, 93, &
      93, 93, 93, 93, 93, 93, 94, 94, 94, 94, 94, 94, &
      94, 94, 94, 94, 94, 94, 94, 94, 95, 95, 95, 95, &
      95, 95, 95, 95, 95, 95, 95, 96, 96, 96, 96, 96, &
      96, 96, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, &
      97, 98, 98, 98, 98, 99, 99, 99, 99, 99, 99, 99, &
      99, 99, 99, 99, 99, 99,100,100,100,100,101,101, &
     101,101,102,102,102,102,102,102,102,102,102,102, &
     102,102,102,102,102,102,102,102,102,102,102,102, &
     102,102,102,102,102,102,102,102,102,102,103,103, &
     103,103,103,103,103,103,103,103,103,103,104,104, &
     104,104,105,105,105,105,106,106,106,106,107,107, &
     107,107,107,107,107,107,107,107,107,107,107,108, &
     108,108,108,108,108,108,108,108,108,108,108,108, &
     108,108,108,108,108,108,108,108,108,108,108,108, &
     108,108,108,108,108,108,109,109,109,109,109,109, &
     109,109,109,109,109,109,109,109,109,109,109,109, &
     109,109,109,109,109,109,109,109,109,109,109,109, &
     109,109,109,109,109,109,109,109,109,109,109,109, &
     109,109,109,109,109,109,109,109,109,110,110,110, &
     110,110,110,110,110,110,110,110,110,110,110,110, &
     110,110,110,110,110,110,110,110,111,111,111,111, &
     111,111,111,111,111,112,112,112,112,112,112,112 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_2 = (/ &
     112,112,112,112,112,112,112,112,113,113,113,113, &
     113,113,113,113,113,113,113,113,113,113,113,113, &
     113,113,113,113,113,113,113,113,113,113,114,114, &
     114,114,114,114,115,115,115,115,115,115,115,115, &
     115,115,115,115,115,115,115,115,115,115,115,115, &
     115,115,115,116,116,116,116,116,116,116,116,117, &
     117,117,117,117,117,117,118,118,118,118,118,118, &
     118,118,118,118,118,118,118,118,119,119,119,119, &
     119,119,119,119,119,119,119,119,119,120,120,120, &
     120,120,120,120,120,120,120,120,120,120,120,120, &
     120,120,121,121,121,121,121,121,121,121,121,121, &
     121,121,121,121,121,121,121,121,121,121,122,122, &
     122,122,122,122,122,122,122,122,122,122,122,122, &
     123,123,123,123,123,123,123,123,123,123,123,123, &
     123,123,123,123,123,123,123,123,123,123,123,124, &
     124,124,124,124,124,124,124,124,124,124,124,124, &
     125,125,125,125,125,125,125,125,125,125,126,126, &
     126,126,126,126,126,126,126,126,126,126,126,126, &
     126,126,126,126,126,126,126,126,126,126,126,126, &
     126,126,126,126,126,126,126,126,126,126,126,126, &
     127,127,127,127,127,127,127,127,127,127,127,127, &
     127,127,127,128,128,128,128,128,128,128,128,128, &
     129,129,129,129,129,129,129,129,129,130,130,130, &
     130,130,130,130,130,130,131,131,131,131,131,131, &
     131,131,131,131,131,131,131,132,132,132,132,132, &
     132,132,132,132,132,132,132,132,132,132,132,132, &
     132,132,132,132,132,132,132,132,132,132,132,132, &
     132,132,132,132,132,132,132,132,132,132,132,132, &
     132,132,132,132,132,132,132,132,132,132,132,132, &
     132,132,132,132,132,132,132,132,132,132,132,132 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_3 = (/ &
     133,133,133,133,133,133,133,133,133,133,133,133, &
     133,134,134,134,134,134,134,134,134,134,134,134, &
     134,134,134,134,134,134,134,134,134,134,134,134, &
     134,134,134,134,134,134,134,134,134,134,134,134, &
     134,134,134,134,134,134,134,135,135,135,135,135, &
     135,135,135,135,135,135,135,135,135,135,135,135, &
     135,135,135,135,135,135,135,135,136,136,136,136, &
     136,136,136,136,136,136,136,136,136,136,136,136, &
     136,136,136,136,136,136,136,136,136,136,136,136, &
     136,136,136,136,136,136,136,136,136,136,136,136, &
     136,136,136,136,136,136,136,136,136,136,136,137, &
     137,137,137,137,137,137,137,137,137,137,137,137, &
     137,137,137,137,137,137,137,137,137,137,137,137, &
     137,137,137,137,137,137,137,137,137,137,137,137, &
     137,137,137,137,137,137,137,137,137,137,137,137, &
     137,137,137,137,137,137,137,137,137,137,137,137, &
     137,137,137,137,137,138,138,138,138,138,138,138, &
     138,138,138,138,138,138,138,138,138,138,138,138, &
     138,138,138,138,138,138,138,138,138,138,138,138, &
     138,138,138,138,138,138,138,138,138,138,138,138, &
     138,138,138,138,138,138,138,138,138,138,138,138, &
     138,138,138,138,138,138,138,138,138,138,138,138, &
     138,138,138,138,138,138,138,138,138,138,138,138, &
     138,138,138,138,138,138,138,138,138,138,138,138, &
     138,138,138,138,138,138,138,138,138,138,138,138, &
     138,138,138,138,138,138,138,138,138,138,138,138, &
     138,138,138,138,138,138,138,138,138,139,139,139, &
     139,139,139,139,139,139,139,139,139,139,139,140, &
     140,140,140,140,140,140,140,140,140,140,140,140, &
     140,140,140,140,140,140,140,140,140,140,140,140 /)
  INTEGER, PARAMETER, DIMENSION(182) :: LU_IROW_4 = (/ &
     140,140,140,140,140,140,140,140,140,140,140,140, &
     140,140,140,140,140,141,141,141,141,141,141,141, &
     141,141,141,141,141,141,141,141,141,141,141,141, &
     141,141,141,141,141,141,141,141,141,141,141,141, &
     141,141,141,141,141,141,141,141,141,141,141,141, &
     141,141,141,141,141,141,141,141,141,141,141,141, &
     141,141,141,141,141,141,141,141,141,141,141,141, &
     141,141,141,141,141,141,141,141,141,141,141,141, &
     141,141,141,141,141,141,141,141,141,141,141,141, &
     141,141,141,141,141,141,141,141,141,141,141,141, &
     141,141,141,141,141,141,141,141,141,141,141,142, &
     142,142,142,142,142,142,142,142,142,142,142,142, &
     142,142,142,142,142,142,142,142,142,142,142,142, &
     142,142,142,142,142,142,142,142,142,142,142,142, &
     142,142,142,142,142,142,142,142,142,142,142,142, &
     142,142 /)
  INTEGER, PARAMETER, DIMENSION(1622) :: LU_IROW = (/&
    LU_IROW_0, LU_IROW_1, LU_IROW_2, LU_IROW_3, LU_IROW_4 /)

  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_0 = (/ &
       1,  7,138,  2,138,  3,138,  4,138,  5,138,  6, &
      93,  7, 23,134,138,  8, 75,141,  9, 88,138, 10, &
      65,137, 11,138, 12,138, 12, 13,138, 12, 14,138, &
      15, 93,138, 16,134,137, 17,138, 17, 18,138, 17, &
      19, 21,138, 20,127,138,141, 21,138, 21, 22,138, &
      23,134,138, 24,130,138, 25, 40, 52,138,142, 26, &
     116,124,138,142, 27,122,138,141, 28, 77,138,141, &
      29, 95,138,140,  8, 30, 75,141,142, 31, 79, 85, &
      87, 90,142, 32,138,140,141, 33,137,138,141, 34, &
      92,138,141, 35,130,134,136,138,140,142, 36,125, &
     138,141, 37, 81,138,141, 38, 54,138,141, 39,121, &
     138,141, 40, 52,138,141, 41, 85, 90,142, 42, 95, &
     138,141, 43,114,135,136,138,139,140,141, 44,134, &
     138, 45,100,101,104,105,106,135,138, 46,135,138, &
      47,129,138,141, 48, 61,138,141, 49,129,138,140, &
      50,118,138,141, 51, 62,138,141, 17, 21, 40, 52, &
     138,141,142, 53,132,141,142, 14, 18, 38, 54,138, &
     141,142, 13, 48, 55, 61, 79, 85, 87, 90,138,141, &
     142, 56,136,138,139,141, 57,136,137,138, 14, 18, &
      22, 28, 58, 77,135,137,138,141,142, 59, 82,138, &
     141, 60,119,138,141, 12, 48, 61,138,141,142, 22, &
      51, 62,138,141,142, 63, 78, 80,107,138,142, 64, &
      97,138,141,142, 10, 25, 40, 52, 65,137,138,141, &
     142, 66, 71, 93,122,132,138, 67,111,138,141, 68, &
     137,138,139, 69,130,138,142, 44, 70,134,138,142, &
      71, 93,114,135,138, 72,128,138,142, 24, 73, 79, &
      85, 87, 90,117,130,138,142, 74,116,138,141, 46, &
      75,135,138,141,142, 76,118,119,138,140,142, 28, &
      58, 65, 77,135,137,138,141,142, 78,107,138,141, &
      79, 87,138,141, 80,107,138,142, 37, 63, 78, 80 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_1 = (/ &
      81,107,138,141,142, 59, 82,114,138,141,142, 11, &
      46, 53, 83, 89, 98,100,101,104,105,106,114,131, &
     132,133,135,138,141,142, 84,100,101,104,105,106, &
     135, 85, 90,138,141, 55, 61, 79, 85, 86, 87, 90, &
     137,138,141,142, 17, 79, 87,138,141,142, 47, 88, &
     128,129,130,138,141, 89,134,135,138, 21, 85, 90, &
     138,141,142, 31, 79, 85, 87, 90, 91,137,138,141, &
     142,  4, 34, 92,138,140,141,142,  6, 15, 66, 71, &
      93,114,122,132,135,138, 94, 95,111,116,118,119, &
     121,125,127,128,129,130,138,140,  3, 42, 63, 78, &
      80, 95,107,138,140,141,142, 89, 96,134,135,138, &
     141,142, 67, 89, 96, 97,111,134,135,138,140,141, &
     142, 98,134,135,138, 47, 49, 69, 99,114,117,129, &
     130,134,138,140,141,142,100,134,135,138,101,134, &
     135,138, 11, 38, 48, 51, 54, 61, 62, 79, 84, 85, &
      86, 87, 90,100,101,102,103,104,105,106,110,117, &
     127,130,134,135,136,137,138,140,141,142, 73, 79, &
      85, 87, 90,103,117,130,137,138,141,142,104,134, &
     135,138,105,134,135,138,106,134,135,138, 45, 78, &
     100,101,104,105,106,107,134,135,138,141,142, 34, &
      44, 60, 64, 67, 70, 74, 78, 80, 89, 92, 96, 97, &
     100,101,104,105,106,107,108,111,112,116,119,124, &
     134,135,138,140,141,142, 11, 13, 19, 21, 41, 46, &
      73, 74, 79, 84, 85, 86, 87, 89, 90, 91, 97, 98, &
     100,101,102,103,104,105,106,109,110,111,112,114, &
     115,116,117,120,122,124,125,126,127,130,131,132, &
     133,134,135,136,137,138,140,141,142, 30, 49, 67, &
      69, 72, 74, 75,110,111,116,124,125,127,128,129, &
     130,134,135,136,138,140,141,142, 67, 89,111,134, &
     135,138,140,141,142, 74,100,101,104,105,106,112 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_2 = (/ &
     116,120,134,135,138,140,141,142, 16, 23, 26, 49, &
      50, 69, 72, 76,113,115,116,117,118,119,124,126, &
     128,129,130,132,134,137,138,140,141,142, 98,114, &
     133,134,135,138, 29, 37, 42, 44, 59, 70, 78, 80, &
      81, 82, 92, 95,107,114,115,124,133,134,135,138, &
     140,141,142, 74,116,120,134,138,140,141,142,117, &
     129,134,136,138,140,142, 50,100,101,104,105,106, &
     118,120,134,135,138,140,141,142, 60,100,101,104, &
     105,106,119,134,135,138,140,141,142, 50, 60, 76, &
     100,101,104,105,106,118,119,120,134,135,138,140, &
     141,142, 39,100,101,104,105,106,108,111,112,116, &
     119,120,121,124,134,135,138,140,141,142, 27, 66, &
      71, 93,114,122,132,133,134,135,137,138,141,142, &
      59, 68, 69, 72, 82,114,121,123,124,125,127,128, &
     130,133,134,135,136,137,138,139,140,141,142, 69, &
      70, 72, 96,124,125,128,130,134,135,138,141,142, &
      36,125,131,133,134,136,138,140,141,142, 41, 51, &
      62, 79, 84, 85, 87, 90, 91, 99,100,101,103,104, &
     105,106,114,117,121,123,124,125,126,127,128,129, &
     130,131,133,134,135,136,137,138,139,140,141,142, &
       9, 20, 24, 35, 88,127,128,129,130,134,136,138, &
     140,141,142, 98,128,134,135,136,138,140,141,142, &
      98,129,134,135,136,138,140,141,142, 98,130,134, &
     135,136,138,140,141,142, 88, 98,128,129,130,131, &
     134,135,136,138,140,141,142, 30, 32, 39, 44, 46, &
      53, 56, 57, 59, 60, 67, 68, 70, 71, 72, 74, 75, &
      78, 80, 82, 88, 89, 92, 93, 94, 95, 96, 98, 99, &
     100,101,104,105,106,107,110,111,112,114,116,117, &
     118,119,120,121,122,123,124,125,127,128,129,130, &
     131,132,133,134,135,136,137,138,139,140,141,142 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_3 = (/ &
      88, 98,128,129,130,133,134,135,136,138,140,141, &
     142, 16, 23, 33, 44, 57, 68, 89, 96, 98,100,101, &
     104,105,106,113,114,115,116,117,118,119,120,122, &
     124,125,126,127,128,129,130,131,132,133,134,135, &
     136,137,138,139,140,141,142, 46, 58, 77, 89, 93, &
      98,100,101,104,105,106,114,122,131,132,133,134, &
     135,136,137,138,139,140,141,142, 19, 21, 37, 39, &
      41, 56, 57, 63, 78, 80, 81, 84, 85, 90, 98, 99, &
     100,101,104,105,106,107,108,111,112,114,115,116, &
     117,119,120,121,123,124,125,126,127,128,129,130, &
     131,133,134,135,136,137,138,139,140,141,142, 10, &
      16, 26, 33, 44, 47, 50, 52, 53, 54, 57, 58, 61, &
      62, 64, 65, 68, 70, 72, 75, 76, 77, 80, 81, 82, &
      86, 87, 90, 91, 92, 95, 96, 97, 99,103,107,111, &
     113,114,115,116,117,118,119,120,121,122,124,125, &
     126,127,128,129,130,131,132,133,134,135,136,137, &
     138,139,140,141,142,  2,  3,  4,  5,  7,  8,  9, &
      11, 12, 14, 15, 17, 18, 20, 21, 22, 23, 24, 25, &
      27, 28, 29, 32, 33, 34, 35, 36, 37, 38, 39, 40, &
      42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 54, &
      56, 57, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, &
      69, 71, 72, 74, 75, 76, 77, 78, 79, 80, 81, 82, &
      83, 85, 87, 88, 89, 90, 91, 92, 93, 94, 95, 97, &
      98, 99,100,101,102,103,104,105,106,107,108,109, &
     110,111,112,113,114,115,116,117,118,119,120,121, &
     122,123,124,125,126,127,128,129,130,131,132,133, &
     134,135,136,137,138,139,140,141,142, 36, 68,125, &
     131,133,134,135,136,137,138,139,140,141,142, 32, &
      43, 56, 57, 71, 91, 92, 93, 95, 98,103,107,108, &
     111,112,114,115,116,117,118,119,120,121,122,124 /)
  INTEGER, PARAMETER, DIMENSION(182) :: LU_ICOL_4 = (/ &
     125,127,128,129,130,131,132,133,134,135,136,137, &
     138,139,140,141,142,  5, 11, 12, 13, 14, 15, 17, &
      18, 19, 21, 22, 23, 24, 27, 29, 30, 31, 32, 33, &
      34, 36, 38, 40, 41, 42, 46, 47, 48, 49, 51, 52, &
      53, 54, 55, 59, 60, 61, 62, 64, 65, 66, 67, 68, &
      69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, &
      81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, &
      93, 94, 95, 96, 97, 98,100,101,102,103,104,105, &
     106,107,109,110,111,112,114,115,116,117,118,119, &
     120,121,122,123,124,125,126,127,128,129,130,131, &
     132,133,134,135,136,137,138,139,140,141,142,  6, &
      52, 53, 54, 61, 62, 64, 65, 70, 75, 77, 81, 82, &
      86, 87, 90, 91, 92, 93, 95, 96, 97,103,107,111, &
     114,116,117,118,119,120,121,122,124,125,127,128, &
     129,130,131,132,133,134,135,136,137,138,139,140, &
     141,142 /)
  INTEGER, PARAMETER, DIMENSION(1622) :: LU_ICOL = (/&
    LU_ICOL_0, LU_ICOL_1, LU_ICOL_2, LU_ICOL_3, LU_ICOL_4 /)

  INTEGER, PARAMETER, DIMENSION(143) :: LU_CROW = (/ &
       1,  4,  6,  8, 10, 12, 14, 18, 21, 24, 27, 29, &
      31, 34, 37, 40, 43, 45, 48, 52, 56, 58, 61, 64, &
      67, 72, 77, 81, 85, 89, 94,100,104,108,112,119, &
     123,127,131,135,139,143,147,155,158,166,169,173, &
     177,181,185,189,196,200,207,218,223,227,238,242, &
     246,252,258,264,269,278,284,288,292,296,301,306, &
     310,320,324,330,336,345,349,353,357,366,372,391, &
     398,402,413,419,426,430,436,446,453,463,477,488, &
     495,506,510,523,527,531,563,575,579,583,587,600, &
     631,682,705,714,729,755,761,784,792,799,813,826, &
     843,863,877,900,913,923,961,976,985,994,1003,1016, &
     1081,1094,1136,1161,1212,1278,1402,1416,1458,1572,1623 /)

  INTEGER, PARAMETER, DIMENSION(143) :: LU_DIAG = (/ &
       1,  4,  6,  8, 10, 12, 14, 18, 21, 24, 27, 29, &
      32, 35, 37, 40, 43, 46, 49, 52, 56, 59, 61, 64, &
      67, 72, 77, 81, 85, 90, 94,100,104,108,112,119, &
     123,127,131,135,139,143,147,155,158,166,169,173, &
     177,181,185,192,196,203,209,218,223,231,238,242, &
     248,254,258,264,273,278,284,288,292,297,301,306, &
     311,320,325,330,339,345,349,353,361,367,375,391, &
     398,406,415,420,426,432,441,448,457,463,482,489, &
     498,506,513,523,527,546,568,575,579,583,594,619, &
     656,689,707,720,737,756,775,785,792,805,819,836, &
     855,868,884,904,914,945,966,977,986,995,1008,1070, &
     1086,1127,1153,1205,1272,1397,1412,1455,1570,1622,1623 /)


END MODULE t1_mozcart_JacobianSP

