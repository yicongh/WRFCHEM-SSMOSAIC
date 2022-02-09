























MODULE saprc99_mosaic_8bin_vbs2_aq_JacobianSP

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
      17, 17, 17, 17, 17, 17, 18, 18, 18, 18, 18, 18, &
      18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, &
      18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, &
      18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, &
      18, 18, 18, 18, 18, 18, 18, 19, 19, 19, 19, 19, &
      19, 19, 19, 19, 19, 19, 19, 19, 19, 20, 20, 20, &
      20, 21, 21, 22, 22, 22, 22, 23, 23, 24, 24, 25, &
      25, 26, 26, 26, 27, 27, 28, 28, 28, 29, 29, 30, &
      30, 31, 31, 31, 32, 32, 32, 33, 33, 33, 34, 34, &
      34, 35, 35, 35, 36, 36, 37, 37, 37, 37, 37, 37, &
      38, 38, 39, 39, 39, 40, 40, 40, 40, 41, 41, 41, &
      42, 42, 42, 43, 43, 44, 44, 45, 45, 45, 45, 46, &
      46, 46, 46, 47, 47, 47, 47, 48, 48, 48, 48, 48, &
      49, 49, 50, 50, 50, 50, 51, 51, 51, 51, 51, 52, &
      52, 53, 53, 54, 54, 54, 54, 55, 55, 55, 55, 56, &
      56, 56, 56, 56, 57, 57, 57, 57, 57, 58, 58, 58, &
      59, 59, 59, 59, 59, 60, 60, 60, 60, 60, 60, 61, &
      61, 61, 61, 61, 61, 61, 62, 62, 62, 62, 62, 62, &
      63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, &
      64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, &
      64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, &
      64, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_1 = (/ &
      65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 66, &
      66, 66, 66, 66, 67, 67, 67, 67, 67, 67, 67, 67, &
      67, 67, 67, 67, 67, 67, 67, 68, 68, 68, 68, 68, &
      69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, &
      70, 70, 70, 70, 70, 71, 71, 71, 71, 71, 72, 72, &
      72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, &
      72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, &
      72, 72, 72, 72, 73, 73, 73, 73, 73, 73, 73, 74, &
      74, 74, 74, 74, 75, 75, 75, 75, 75, 76, 76, 76, &
      76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, &
      76, 76, 76, 76, 76, 76, 76, 76, 76, 77, 77, 77, &
      77, 77, 77, 77, 78, 78, 78, 78, 78, 78, 78, 78, &
      78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, &
      78, 78, 78, 78, 79, 79, 79, 79, 79, 80, 80, 80, &
      80, 80, 80, 80, 81, 81, 81, 81, 81, 81, 81, 81, &
      81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, &
      81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, &
      81, 81, 81, 81, 81, 81, 82, 82, 82, 82, 82, 82, &
      82, 82, 82, 82, 82, 82, 82, 82, 83, 83, 83, 83, &
      83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, &
      84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, &
      84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, &
      84, 84, 84, 84, 84, 84, 85, 85, 85, 85, 85, 85, &
      85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, &
      85, 85, 85, 86, 86, 86, 86, 86, 86, 86, 86, 86, &
      86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, &
      86, 86, 86, 87, 87, 87, 87, 87, 87, 87, 87, 87, &
      87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, &
      87, 87, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, &
      88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_2 = (/ &
      88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 89, 89, &
      89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, &
      89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 90, &
      90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, &
      90, 90, 90, 90, 90, 90, 90, 91, 91, 91, 91, 91, &
      91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, &
      91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, &
      91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, &
      91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 92, &
      92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, &
      92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, &
      92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, &
      92, 92, 92, 92, 92, 93, 93, 93, 93, 93, 93, 93, &
      93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, &
      93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, &
      93, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, &
      94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, &
      94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, &
      94, 94, 94, 94, 94, 94, 94, 94, 94, 95, 95, 95, &
      95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, &
      95, 95, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, &
      96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, &
      96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, &
      97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, &
      97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, &
      97, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, &
      98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, &
      98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, &
      98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, &
      98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98 /)
  INTEGER, PARAMETER, DIMENSION(50) :: LU_IROW_3 = (/ &
      98, 98, 98, 98, 98, 98, 98, 98, 99, 99, 99, 99, &
      99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, &
      99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, &
      99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, &
      99, 99 /)
  INTEGER, PARAMETER, DIMENSION(1130) :: LU_IROW = (/&
    LU_IROW_0, LU_IROW_1, LU_IROW_2, LU_IROW_3 /)

  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_0 = (/ &
       1, 40, 98,  2, 58, 68, 87,  3,  4,  5,  6,  7, &
      47, 58, 66, 68, 70, 71, 73, 74, 75, 77, 79, 80, &
      87, 97, 98,  8, 68, 71, 79, 87, 88, 91, 93, 94, &
      96,  9, 70, 71, 74, 75, 77, 79, 87, 88, 89, 90, &
      91, 93, 94, 95, 10, 91, 96, 11, 89, 90, 91, 95, &
      12, 48, 68, 92, 99, 13, 48, 61, 68, 83, 87, 92, &
      98, 99, 14, 72, 88, 93, 94, 97, 99, 15, 72, 88, &
      91, 93, 94, 16, 43, 49, 52, 53, 71, 79, 98, 17, &
      70, 74, 75, 87, 98, 99, 18, 29, 30, 36, 38, 40, &
      42, 43, 44, 46, 49, 51, 52, 53, 54, 55, 56, 57, &
      58, 59, 61, 62, 63, 64, 65, 66, 67, 68, 70, 71, &
      73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 84, 85, &
      86, 87, 91, 92, 97, 98, 99,  3,  4,  5,  6, 19, &
      20, 21, 22, 23, 25, 26, 27, 28, 98, 20, 21, 25, &
      98, 21, 98, 22, 23, 27, 98, 23, 98, 24, 87, 25, &
      98, 25, 26, 98, 27, 98, 27, 28, 98, 29, 98, 30, &
      98, 31, 92, 96, 32, 89, 92, 33, 92, 95, 34, 90, &
      92, 35, 91, 98, 36, 98, 37, 49, 74, 75, 87, 98, &
      38, 98, 39, 92, 99, 40, 41, 98, 99, 41, 98, 99, &
      42, 97, 98, 43, 98, 44, 98, 44, 45, 92, 98, 46, &
      91, 93, 98, 47, 81, 91, 97, 48, 60, 91, 92, 99, &
      49, 98, 50, 91, 92, 98, 51, 88, 93, 94, 98, 52, &
      98, 53, 98, 49, 53, 54, 98, 49, 53, 55, 98, 49, &
      53, 56, 98, 99, 49, 53, 57, 87, 98, 58, 87, 98, &
      59, 88, 91, 94, 98, 48, 60, 69, 91, 92, 99, 49, &
      53, 61, 79, 87, 98, 99, 53, 62, 69, 91, 98, 99, &
      49, 53, 54, 55, 56, 63, 73, 77, 80, 87, 98, 99, &
      52, 54, 55, 57, 58, 63, 64, 66, 68, 70, 71, 73, &
      74, 75, 76, 77, 78, 79, 80, 81, 83, 84, 87, 98, &
      99, 39, 41, 56, 60, 61, 62, 63, 65, 69, 73, 76 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_1 = (/ &
      77, 78, 79, 80, 81, 84, 87, 91, 92, 98, 99, 66, &
      83, 87, 98, 99, 36, 43, 44, 45, 52, 67, 71, 74, &
      75, 79, 82, 87, 92, 98, 99, 68, 83, 87, 98, 99, &
      56, 62, 69, 89, 90, 91, 92, 95, 96, 97, 98, 99, &
      70, 83, 87, 98, 99, 71, 83, 87, 98, 99, 43, 44, &
      52, 54, 55, 67, 70, 71, 72, 74, 75, 79, 80, 82, &
      83, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, &
      96, 97, 98, 99, 70, 73, 79, 83, 87, 98, 99, 74, &
      83, 87, 98, 99, 75, 83, 87, 98, 99, 49, 53, 54, &
      55, 57, 58, 62, 66, 69, 74, 75, 76, 77, 83, 87, &
      89, 90, 91, 92, 95, 96, 97, 98, 99, 70, 77, 79, &
      83, 87, 98, 99, 30, 38, 43, 44, 52, 66, 68, 71, &
      77, 78, 79, 82, 83, 84, 85, 86, 87, 89, 90, 95, &
      96, 97, 98, 99, 79, 83, 87, 98, 99, 70, 79, 80, &
      83, 87, 98, 99, 38, 43, 44, 46, 47, 51, 52, 58, &
      66, 67, 68, 70, 71, 73, 74, 75, 76, 77, 79, 80, &
      81, 82, 83, 85, 86, 87, 88, 89, 90, 91, 92, 93, &
      94, 95, 96, 97, 98, 99, 45, 71, 74, 75, 77, 79, &
      82, 83, 87, 88, 92, 97, 98, 99, 24, 66, 68, 70, &
      71, 73, 74, 75, 79, 80, 83, 87, 92, 97, 98, 99, &
      36, 43, 44, 52, 54, 55, 57, 59, 66, 68, 71, 73, &
      74, 75, 77, 79, 80, 82, 83, 84, 85, 86, 87, 88, &
      91, 92, 94, 97, 98, 99, 43, 44, 52, 68, 71, 73, &
      77, 79, 80, 82, 83, 85, 86, 87, 88, 92, 93, 94, &
      97, 98, 99, 43, 52, 53, 70, 71, 74, 75, 77, 79, &
      80, 82, 83, 86, 87, 88, 89, 92, 93, 94, 95, 96, &
      97, 98, 99, 57, 58, 66, 68, 70, 71, 73, 74, 75, &
      77, 79, 80, 83, 87, 89, 90, 91, 92, 95, 96, 97, &
      98, 99, 36, 43, 44, 49, 52, 53, 68, 70, 71, 74, &
      75, 77, 79, 80, 82, 83, 84, 85, 86, 87, 88, 89 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_2 = (/ &
      90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 32, 73, &
      74, 75, 76, 77, 79, 80, 83, 84, 85, 86, 87, 88, &
      89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 34, &
      70, 73, 77, 79, 80, 83, 87, 88, 89, 90, 91, 92, &
      93, 94, 95, 96, 97, 98, 99, 35, 38, 40, 41, 42, &
      46, 47, 48, 49, 50, 51, 53, 54, 55, 57, 58, 59, &
      60, 63, 64, 66, 68, 69, 70, 71, 73, 74, 75, 76, &
      77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, &
      89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 31, &
      32, 33, 34, 39, 42, 45, 47, 48, 50, 60, 65, 69, &
      70, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, &
      83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, &
      95, 96, 97, 98, 99, 29, 45, 46, 52, 66, 67, 68, &
      70, 71, 74, 75, 78, 79, 80, 82, 83, 84, 85, 86, &
      87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, &
      99, 30, 36, 38, 43, 44, 49, 52, 53, 54, 55, 56, &
      57, 58, 59, 62, 66, 68, 69, 70, 71, 73, 74, 75, &
      77, 79, 80, 82, 83, 84, 85, 86, 87, 88, 89, 90, &
      91, 92, 93, 94, 95, 96, 97, 98, 99, 33, 61, 79, &
      83, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, &
      98, 99, 31, 37, 49, 52, 54, 55, 63, 67, 71, 73, &
      74, 75, 77, 78, 79, 80, 82, 83, 84, 85, 86, 87, &
      88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, &
      42, 47, 72, 74, 75, 79, 80, 81, 82, 83, 85, 86, &
      87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, &
      99,  3,  5, 24, 25, 26, 27, 28, 29, 30, 35, 36, &
      38, 40, 41, 42, 43, 44, 46, 49, 50, 51, 52, 53, &
      54, 55, 56, 57, 58, 59, 61, 62, 63, 64, 65, 66, &
      67, 68, 69, 70, 71, 73, 74, 75, 76, 77, 78, 79, &
      80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91 /)
  INTEGER, PARAMETER, DIMENSION(50) :: LU_ICOL_3 = (/ &
      92, 93, 94, 95, 96, 97, 98, 99, 39, 41, 50, 56, &
      60, 61, 62, 63, 65, 66, 68, 69, 70, 71, 72, 73, &
      74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, &
      86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, &
      98, 99 /)
  INTEGER, PARAMETER, DIMENSION(1130) :: LU_ICOL = (/&
    LU_ICOL_0, LU_ICOL_1, LU_ICOL_2, LU_ICOL_3 /)

  INTEGER, PARAMETER, DIMENSION(100) :: LU_CROW = (/ &
       1,  4,  8,  9, 10, 11, 12, 28, 38, 53, 56, 61, &
      66, 75, 82, 88, 96,103,152,166,170,172,176,178, &
     180,182,185,187,190,192,194,197,200,203,206,209, &
     211,217,219,222,226,229,232,234,236,240,244,248, &
     253,255,259,264,266,268,272,276,281,286,289,294, &
     300,307,313,325,350,372,377,392,397,409,414,419, &
     449,456,461,466,490,497,521,526,533,571,585,601, &
     631,652,676,699,731,756,776,828,870,902,946,963, &
     997,1022,1089,1131 /)

  INTEGER, PARAMETER, DIMENSION(100) :: LU_DIAG = (/ &
       1,  4,  8,  9, 10, 11, 12, 28, 38, 53, 56, 61, &
      66, 75, 82, 88, 96,103,156,166,170,172,176,178, &
     180,183,185,188,190,192,194,197,200,203,206,209, &
     211,217,219,222,226,229,232,234,237,240,244,248, &
     253,255,259,264,266,270,274,278,283,286,289,295, &
     302,308,318,331,357,372,382,392,399,409,414,427, &
     450,456,461,477,491,506,521,528,553,577,595,620, &
     642,664,689,719,745,766,819,862,895,940,958,993, &
     1019,1087,1130,1131 /)


END MODULE saprc99_mosaic_8bin_vbs2_aq_JacobianSP

