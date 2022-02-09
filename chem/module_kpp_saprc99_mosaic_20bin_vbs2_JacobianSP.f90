























MODULE saprc99_mosaic_20bin_vbs2_JacobianSP

  PUBLIC
  SAVE





  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_0 = (/ &
       1,  1,  1,  2,  2,  2,  2,  3,  4,  5,  6,  7, &
       7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7, &
       7,  7,  7,  8,  8,  8,  8,  8,  8,  8,  8,  8, &
       8,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9, &
       9,  9,  9,  9, 10, 10, 10, 11, 11, 11, 11, 11, &
      12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 13, &
      13, 13, 14, 14, 14, 14, 14, 14, 14, 15, 15, 15, &
      15, 15, 15, 16, 16, 16, 16, 16, 16, 16, 16, 17, &
      17, 17, 18, 18, 18, 19, 19, 19, 19, 19, 19, 19, &
      20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, &
      20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, &
      20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, &
      20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, &
      20, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, &
      21, 21, 21, 22, 22, 22, 22, 23, 23, 24, 24, 24, &
      24, 25, 25, 26, 26, 27, 27, 28, 28, 28, 29, 29, &
      30, 30, 30, 31, 31, 32, 32, 33, 33, 33, 34, 34, &
      34, 35, 35, 35, 36, 36, 36, 37, 37, 37, 38, 38, &
      39, 39, 39, 39, 39, 39, 40, 40, 41, 41, 41, 42, &
      42, 42, 42, 43, 43, 43, 44, 44, 44, 45, 45, 46, &
      46, 46, 46, 47, 47, 48, 48, 48, 48, 49, 49, 49, &
      49, 50, 50, 50, 50, 50, 51, 51, 51, 51, 52, 52, &
      53, 53, 53, 53, 53, 54, 54, 55, 55, 56, 56, 56, &
      56, 57, 57, 57, 57, 58, 58, 58, 58, 58, 59, 59, &
      59, 59, 59, 60, 60, 60, 61, 61, 61, 61, 61, 61, &
      62, 62, 62, 62, 62, 63, 63, 63, 63, 63, 63, 63, &
      64, 64, 64, 64, 64, 64, 65, 65, 65, 65, 65, 65, &
      65, 65, 65, 65, 65, 65, 66, 66, 66, 66, 66, 66, &
      66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, &
      66, 66, 66, 66, 66, 66, 66, 67, 67, 67, 67, 67 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_1 = (/ &
      67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, &
      67, 67, 67, 67, 67, 68, 68, 68, 68, 68, 69, 69, &
      69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, &
      69, 70, 70, 70, 70, 70, 71, 71, 71, 71, 71, 71, &
      71, 71, 71, 71, 71, 71, 72, 72, 72, 72, 72, 73, &
      73, 73, 73, 73, 74, 74, 74, 74, 74, 74, 74, 74, &
      74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, &
      74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 75, 75, &
      75, 75, 75, 75, 75, 76, 76, 76, 76, 76, 77, 77, &
      77, 77, 77, 78, 78, 78, 78, 78, 78, 78, 78, 78, &
      78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, &
      78, 78, 78, 79, 79, 79, 79, 79, 79, 79, 80, 80, &
      80, 80, 80, 81, 81, 81, 81, 81, 81, 81, 82, 82, &
      82, 82, 82, 82, 82, 82, 82, 82, 82, 82, 82, 82, &
      82, 82, 82, 82, 82, 82, 82, 82, 82, 82, 83, 83, &
      83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, &
      83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, &
      83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, &
      84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, &
      84, 84, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, &
      85, 85, 85, 85, 85, 85, 86, 86, 86, 86, 86, 86, &
      86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, &
      86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, &
      87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, &
      87, 87, 87, 87, 87, 87, 87, 87, 87, 88, 88, 88, &
      88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, &
      88, 88, 88, 88, 88, 88, 88, 88, 88, 89, 89, 89, &
      89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, &
      89, 89, 89, 89, 89, 89, 89, 89, 90, 90, 90, 90, &
      90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_2 = (/ &
      90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, &
      90, 90, 90, 90, 90, 90, 91, 91, 91, 91, 91, 91, &
      91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 92, &
      92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, &
      92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, &
      93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, &
      93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, &
      93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, &
      93, 93, 93, 93, 93, 93, 94, 94, 94, 94, 94, 94, &
      94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, &
      94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, &
      94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, &
      95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, &
      95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, &
      95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, &
      95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, &
      95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, &
      95, 95, 95, 95, 95, 95, 95, 96, 96, 96, 96, 96, &
      96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, &
      96, 96, 96, 96, 96, 96, 96, 96, 97, 97, 97, 97, &
      97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, &
      97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, &
      97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, &
      97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, &
      98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, &
      98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, &
      98, 98, 98, 98, 98, 98, 98, 98, 99, 99, 99, 99, &
      99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, &
      99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, &
      99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99 /)
  INTEGER, PARAMETER, DIMENSION(56) :: LU_IROW_3 = (/ &
      99, 99, 99, 99,100,100,100,100,100,100,100,100, &
     100,100,100,100,100,100,100,100,100,100,100,100, &
     101,101,101,101,101,101,101,101,101,101,101,101, &
     101,101,101,101,101,101,101,101,101,101,101,101, &
     101,101,101,101,101,101,101,101 /)
  INTEGER, PARAMETER, DIMENSION(1136) :: LU_IROW = (/&
    LU_IROW_0, LU_IROW_1, LU_IROW_2, LU_IROW_3 /)

  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_0 = (/ &
       1, 42, 95,  2, 60, 70, 89,  3,  4,  5,  6,  7, &
      49, 60, 68, 70, 72, 73, 75, 76, 77, 79, 80, 81, &
      89, 95, 96,  8, 70, 72, 80, 89, 90, 97, 98, 99, &
     101,  9, 72, 73, 76, 77, 79, 80, 89, 91, 92, 97, &
      98, 99,100,101, 10, 90, 97, 11, 91, 92, 97,100, &
      12, 50, 70, 93, 94, 13, 50, 63, 70, 85, 89, 93, &
      94, 95, 14, 74, 94, 96, 98, 99,101, 15, 74, 97, &
      98, 99,101, 16, 47, 52, 54, 55, 72, 80, 95, 17, &
      77, 95, 18, 77, 89, 19, 73, 76, 77, 89, 94, 95, &
      20, 31, 32, 38, 40, 42, 44, 45, 47, 48, 52, 53, &
      54, 55, 56, 57, 58, 59, 60, 62, 63, 64, 65, 66, &
      67, 68, 69, 70, 72, 73, 75, 76, 77, 78, 79, 80, &
      81, 82, 83, 84, 86, 87, 88, 89, 93, 94, 95, 96, &
      97,  3,  4,  5,  6, 21, 22, 23, 24, 25, 27, 28, &
      29, 30, 95, 22, 23, 27, 95, 23, 95, 24, 25, 29, &
      95, 25, 95, 26, 89, 27, 95, 27, 28, 95, 29, 95, &
      29, 30, 95, 31, 95, 32, 95, 33, 90, 93, 34, 92, &
      93, 35, 91, 93, 36, 93,100, 37, 95, 97, 38, 95, &
      39, 52, 76, 77, 89, 95, 40, 95, 41, 93, 94, 42, &
      43, 94, 95, 43, 94, 95, 44, 95, 96, 45, 95, 45, &
      46, 93, 95, 47, 95, 48, 95, 97, 98, 49, 83, 96, &
      97, 50, 61, 93, 94, 97, 51, 93, 95, 97, 52, 95, &
      53, 95, 98, 99,101, 54, 95, 55, 95, 52, 55, 56, &
      95, 52, 55, 57, 95, 52, 55, 58, 94, 95, 52, 55, &
      59, 89, 95, 60, 89, 95, 50, 61, 71, 93, 94, 97, &
      62, 95, 97, 99,101, 52, 55, 63, 80, 89, 94, 95, &
      55, 64, 71, 94, 95, 97, 52, 55, 56, 57, 58, 65, &
      75, 79, 81, 89, 94, 95, 54, 56, 57, 59, 60, 65, &
      66, 68, 70, 72, 73, 75, 76, 77, 78, 79, 80, 81, &
      82, 83, 85, 86, 89, 94, 95, 41, 43, 58, 61, 63 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_1 = (/ &
      64, 65, 67, 71, 75, 78, 79, 80, 81, 82, 83, 86, &
      89, 93, 94, 95, 97, 68, 85, 89, 94, 95, 38, 45, &
      46, 47, 54, 69, 72, 76, 77, 80, 84, 89, 93, 94, &
      95, 70, 85, 89, 94, 95, 58, 64, 71, 90, 91, 92, &
      93, 94, 95, 96, 97,100, 72, 85, 89, 94, 95, 73, &
      85, 89, 94, 95, 45, 47, 54, 56, 57, 69, 72, 73, &
      74, 76, 77, 80, 81, 84, 85, 87, 88, 89, 90, 91, &
      92, 93, 94, 95, 96, 97, 98, 99,100,101, 73, 75, &
      80, 85, 89, 94, 95, 76, 85, 89, 94, 95, 77, 85, &
      89, 94, 95, 52, 55, 56, 57, 59, 60, 64, 68, 71, &
      76, 77, 78, 79, 85, 89, 90, 91, 92, 93, 94, 95, &
      96, 97,100, 73, 79, 80, 85, 89, 94, 95, 80, 85, &
      89, 94, 95, 73, 80, 81, 85, 89, 94, 95, 32, 40, &
      45, 47, 54, 68, 70, 72, 79, 80, 82, 84, 85, 86, &
      87, 88, 89, 90, 91, 92, 94, 95, 96,100, 40, 45, &
      47, 48, 49, 53, 54, 60, 68, 69, 70, 72, 73, 75, &
      76, 77, 78, 79, 80, 81, 83, 84, 85, 87, 88, 89, &
      90, 91, 92, 93, 94, 95, 96, 97, 98, 99,100,101, &
      46, 72, 76, 77, 79, 80, 84, 85, 89, 93, 94, 95, &
      96,101, 26, 68, 70, 72, 73, 75, 76, 77, 80, 81, &
      85, 89, 93, 94, 95, 96, 38, 45, 47, 54, 56, 57, &
      59, 62, 68, 70, 72, 75, 76, 77, 79, 80, 81, 84, &
      85, 86, 87, 88, 89, 93, 94, 95, 96, 97, 99,101, &
      45, 47, 54, 70, 72, 75, 79, 80, 81, 84, 85, 87, &
      88, 89, 93, 94, 95, 96, 98, 99,101, 47, 54, 55, &
      72, 73, 76, 77, 79, 80, 81, 84, 85, 88, 89, 90, &
      91, 92, 93, 94, 95, 96, 98, 99,101, 59, 60, 68, &
      70, 72, 73, 75, 76, 77, 79, 80, 81, 85, 89, 90, &
      91, 92, 93, 94, 95, 96, 97,100, 33, 39, 52, 54, &
      56, 57, 65, 69, 72, 75, 76, 77, 79, 80, 81, 82 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_2 = (/ &
      84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, &
      96, 97, 98, 99,100,101, 35, 63, 80, 85, 89, 90, &
      91, 92, 93, 94, 95, 96, 97, 98, 99,100,101, 34, &
      75, 76, 77, 78, 79, 80, 81, 85, 86, 87, 88, 89, &
      90, 91, 92, 93, 94, 95, 96, 97, 98, 99,100,101, &
      33, 34, 35, 36, 41, 44, 46, 49, 50, 51, 61, 67, &
      71, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, &
      84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, &
      96, 97, 98, 99,100,101, 41, 43, 51, 58, 61, 63, &
      64, 65, 67, 68, 70, 71, 72, 73, 74, 75, 76, 77, &
      78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, &
      90, 91, 92, 93, 94, 95, 96, 97, 98, 99,100,101, &
       3,  5, 26, 27, 28, 29, 30, 31, 32, 37, 38, 40, &
      42, 43, 44, 45, 47, 48, 51, 52, 53, 54, 55, 56, &
      57, 58, 59, 60, 62, 63, 64, 65, 66, 67, 68, 69, &
      70, 71, 72, 73, 75, 76, 77, 78, 79, 80, 81, 82, &
      83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, &
      95, 96, 97, 98, 99,100,101, 44, 49, 74, 76, 77, &
      80, 81, 83, 84, 85, 87, 88, 89, 90, 91, 92, 93, &
      94, 95, 96, 97, 98, 99,100,101, 37, 40, 42, 43, &
      44, 48, 49, 50, 51, 52, 53, 55, 56, 57, 59, 60, &
      61, 62, 65, 66, 68, 70, 71, 72, 73, 75, 76, 77, &
      78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, &
      90, 91, 92, 93, 94, 95, 96, 97, 98, 99,100,101, &
      31, 46, 48, 54, 68, 69, 70, 72, 73, 76, 77, 80, &
      81, 82, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, &
      94, 95, 96, 97, 98, 99,100,101, 32, 38, 40, 45, &
      47, 52, 54, 55, 56, 57, 58, 59, 60, 62, 64, 68, &
      70, 71, 72, 73, 75, 76, 77, 79, 80, 81, 84, 85, &
      86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97 /)
  INTEGER, PARAMETER, DIMENSION(56) :: LU_ICOL_3 = (/ &
      98, 99,100,101, 36, 73, 75, 79, 80, 81, 85, 89, &
      90, 91, 92, 93, 94, 95, 96, 97, 98, 99,100,101, &
      38, 45, 47, 52, 54, 55, 70, 72, 73, 76, 77, 79, &
      80, 81, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, &
      94, 95, 96, 97, 98, 99,100,101 /)
  INTEGER, PARAMETER, DIMENSION(1136) :: LU_ICOL = (/&
    LU_ICOL_0, LU_ICOL_1, LU_ICOL_2, LU_ICOL_3 /)

  INTEGER, PARAMETER, DIMENSION(102) :: LU_CROW = (/ &
       1,  4,  8,  9, 10, 11, 12, 28, 38, 53, 56, 61, &
      66, 75, 82, 88, 96, 99,102,109,158,172,176,178, &
     182,184,186,188,191,193,196,198,200,203,206,209, &
     212,215,217,223,225,228,232,235,238,240,244,246, &
     250,254,259,263,265,270,272,274,278,282,287,292, &
     295,301,306,313,319,331,356,378,383,398,403,415, &
     420,425,455,462,467,472,496,503,508,515,539,577, &
     591,607,637,658,682,705,739,756,781,823,865,932, &
     957,1009,1041,1085,1105,1137 /)

  INTEGER, PARAMETER, DIMENSION(102) :: LU_DIAG = (/ &
       1,  4,  8,  9, 10, 11, 12, 28, 38, 53, 56, 61, &
      66, 75, 82, 88, 96, 99,102,109,162,172,176,178, &
     182,184,186,189,191,194,196,198,200,203,206,209, &
     212,215,217,223,225,228,232,235,238,241,244,246, &
     250,254,259,263,265,270,272,276,280,284,289,292, &
     296,301,308,314,324,337,363,378,388,398,405,415, &
     420,433,456,462,467,483,497,503,510,525,559,583, &
     601,626,648,670,695,727,745,771,814,857,925,951, &
     1004,1037,1082,1103,1136,1137 /)


END MODULE saprc99_mosaic_20bin_vbs2_JacobianSP

