























MODULE saprc99_simplesom_mosaic_4bin_aq_JacobianSP

  PUBLIC
  SAVE





  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_0 = (/ &
       1,  1,  1,  2,  2,  2,  2,  2,  2,  2,  2,  2, &
       2,  2,  2,  2,  2,  2,  2,  3,  3,  3,  3,  3, &
       3,  3,  3,  3,  3,  4,  4,  4,  4,  4,  4,  4, &
       4,  4,  4,  4,  4,  4,  4,  4,  5,  5,  5,  5, &
       6,  6,  6,  7,  7,  7,  7,  7,  8,  8,  8,  8, &
       8,  9,  9,  9,  9,  9,  9,  9,  9,  9, 10, 10, &
      10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 12, &
      12, 12, 13, 13, 14, 14, 14, 14, 14, 14, 14, 15, &
      15, 15, 15, 15, 16, 16, 16, 17, 17, 18, 18, 19, &
      19, 20, 20, 20, 20, 21, 21, 21, 21, 21, 21, 22, &
      22, 23, 23, 23, 24, 24, 24, 25, 25, 25, 26, 26, &
      26, 27, 27, 27, 28, 28, 28, 28, 28, 28, 29, 29, &
      30, 30, 31, 31, 31, 32, 32, 32, 32, 33, 33, 33, &
      34, 34, 34, 35, 35, 36, 36, 36, 36, 37, 37, 38, &
      38, 39, 39, 39, 39, 39, 40, 40, 40, 40, 41, 41, &
      41, 41, 42, 42, 42, 42, 42, 43, 43, 43, 43, 44, &
      44, 45, 45, 46, 46, 46, 46, 47, 47, 47, 47, 48, &
      48, 48, 48, 48, 49, 49, 49, 49, 49, 50, 50, 50, &
      51, 51, 51, 51, 51, 51, 52, 52, 52, 52, 52, 52, &
      52, 53, 53, 53, 53, 53, 54, 54, 54, 54, 54, 54, &
      54, 54, 55, 55, 55, 55, 55, 55, 56, 56, 56, 56, &
      56, 56, 56, 56, 56, 56, 57, 57, 57, 57, 57, 57, &
      57, 57, 57, 57, 57, 57, 57, 57, 57, 58, 58, 58, &
      58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, &
      59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, &
      59, 59, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, &
      60, 60, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, &
      61, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 63, &
      63, 63, 63, 63, 63, 63, 63, 63, 64, 64, 64, 64, &
      64, 65, 65, 65, 65, 65, 65, 65, 65, 66, 66, 66 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_1 = (/ &
      66, 66, 66, 66, 67, 67, 67, 67, 67, 68, 68, 68, &
      68, 68, 68, 68, 68, 68, 68, 68, 68, 69, 69, 69, &
      69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, &
      69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 70, 70, &
      70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, &
      70, 70, 70, 70, 70, 70, 70, 70, 71, 71, 71, 71, &
      71, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, &
      72, 72, 72, 72, 73, 73, 73, 73, 73, 74, 74, 74, &
      74, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, &
      75, 75, 75, 75, 75, 75, 75, 76, 76, 76, 76, 76, &
      76, 76, 76, 76, 76, 76, 76, 77, 77, 77, 77, 77, &
      78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, &
      78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, &
      78, 78, 78, 78, 78, 78, 79, 79, 79, 79, 79, 80, &
      80, 80, 80, 80, 81, 81, 81, 81, 81, 81, 81, 82, &
      82, 82, 82, 82, 82, 82, 83, 83, 83, 83, 83, 84, &
      84, 84, 84, 84, 84, 84, 85, 85, 85, 85, 85, 85, &
      85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, &
      85, 85, 85, 85, 85, 85, 86, 86, 86, 86, 86, 86, &
      86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, &
      86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, &
      86, 86, 86, 86, 86, 86, 86, 86, 86, 87, 87, 87, &
      87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 88, &
      88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, &
      88, 88, 88, 89, 89, 89, 89, 89, 89, 89, 89, 89, &
      89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, &
      89, 89, 89, 89, 89, 89, 89, 89, 89, 90, 90, 90, &
      90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, &
      90, 90, 90, 90, 90, 90, 91, 91, 91, 91, 91, 91, &
      91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_2 = (/ &
      91, 91, 91, 91, 91, 91, 92, 92, 92, 92, 92, 92, &
      92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, &
      92, 92, 92, 92, 92, 93, 93, 93, 93, 93, 93, 93, &
      93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, &
      93, 93, 93, 93, 93, 93, 94, 94, 94, 94, 94, 94, &
      94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, &
      94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, &
      94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, &
      95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, &
      95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, &
      95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, &
      95, 95, 95, 95, 95, 95, 96, 96, 96, 96, 96, 96, &
      96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 97, &
      97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, &
      97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, &
      97, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, &
      98, 98, 98, 98, 98, 98, 98, 98, 98, 99, 99, 99, &
      99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, &
      99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, &
      99, 99, 99, 99, 99,100,100,100,100,100,100,100, &
     100,100,100,100,100,100,100,100,100,100,100,100, &
     100,100,100,100,100,100,100,100,100,100,100,100, &
     100,100,100,101,101,101,101,101,101,101,101,101, &
     101,101,101,101,101,101,101,101,101,101,101,101, &
     101,101,101,101,101,101,101,101,101,101,101,102, &
     102,102,102,102,102,102,102,102,102,102,102,102, &
     102,102,102,102,102,102,102,102,102,102,102,102, &
     102,102,102,102,102,102,102,102,102,102,102,102, &
     102,102,102,102,102,102,102,102,102,102,102,102, &
     102,102,102,102,102,102,102,102,102,102,102,102 /)
  INTEGER, PARAMETER, DIMENSION(96) :: LU_IROW_3 = (/ &
     103,103,103,103,103,103,103,103,103,103,103,103, &
     103,103,103,103,103,103,103,103,103,103,103,103, &
     103,103,103,103,103,103,103,103,103,103,103,103, &
     103,103,103,103,103,103,103,103,103,103,103,103, &
     103,103,103,103,104,104,104,104,104,104,104,104, &
     104,104,104,104,104,104,104,104,104,104,104,104, &
     104,104,104,104,104,104,104,104,104,104,104,104, &
     104,104,104,104,104,104,104,104,104,104,104,104 /)
  INTEGER, PARAMETER, DIMENSION(1176) :: LU_IROW = (/&
    LU_IROW_0, LU_IROW_1, LU_IROW_2, LU_IROW_3 /)

  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_0 = (/ &
       1, 32,102,  2, 41, 50, 71, 73, 74, 77, 79, 80, &
      81, 82, 83, 84, 92, 93,102,  3, 73, 77, 83, 92, &
      99,100,101,103,104,  4, 74, 77, 79, 80, 82, 83, &
      92, 96, 97, 98, 99,101,103,104,  5, 50, 73, 92, &
       6,100,103,  7, 96, 97, 98,103,  8, 42, 73, 94, &
      95,  9, 42, 52, 73, 88, 92, 94, 95,102, 10, 78, &
      93, 95, 99,101,104, 11, 78, 99,101,103,104, 12, &
      15,103, 13, 15, 14, 15, 16, 17, 21, 93,102, 15, &
      20, 93,102,103, 16, 20,102, 17,102, 18, 92, 19, &
     102, 20, 21,102,103, 20, 21, 80, 93,102,103, 22, &
     102, 23, 94,100, 24, 94, 97, 25, 94, 96, 26, 94, &
      98, 27,102,103, 28, 38, 74, 79, 92,102, 29,102, &
      30,102, 31, 94, 95, 32, 34, 95,102, 33, 93,102, &
      34, 95,102, 35,102, 35, 36, 94,102, 37,102, 38, &
     102, 39, 99,101,102,104, 40, 99,102,103, 41, 86, &
      93,103, 42, 51, 94, 95,103, 43, 94,102,103, 44, &
     102, 45,102, 38, 45, 46,102, 38, 45, 47,102, 38, &
      45, 48, 95,102, 38, 45, 49, 92,102, 50, 92,102, &
      42, 51, 76, 94, 95,103, 38, 45, 52, 83, 92, 95, &
     102, 53,101,102,103,104, 54, 57, 58, 59, 60, 61, &
      62,102, 45, 55, 76, 95,102,103, 54, 56, 57, 58, &
      59, 60, 61, 62, 63,102, 54, 56, 57, 58, 59, 60, &
      61, 62, 63, 64, 65, 66, 67, 74,102, 54, 56, 57, &
      58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 74,102, &
      56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, &
      74,102, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, &
      74,102, 59, 60, 61, 62, 63, 64, 65, 66, 67, 74, &
     102, 60, 61, 62, 63, 64, 65, 66, 67, 74,102, 61, &
      62, 63, 64, 65, 66, 67, 74,102, 64, 66, 67, 74, &
     102, 62, 63, 64, 65, 66, 67, 74,102, 63, 64, 65 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_1 = (/ &
      66, 67, 74,102, 65, 66, 67, 74,102, 38, 45, 46, &
      47, 48, 68, 81, 82, 84, 92, 95,102, 44, 46, 47, &
      49, 50, 68, 69, 71, 73, 74, 75, 77, 79, 80, 81, &
      82, 83, 84, 85, 86, 88, 89, 92, 95,102, 31, 34, &
      48, 51, 52, 55, 68, 70, 75, 76, 81, 82, 83, 84, &
      85, 86, 89, 92, 94, 95,102,103, 71, 88, 92, 95, &
     102, 29, 35, 36, 37, 44, 72, 74, 77, 79, 83, 87, &
      92, 94, 95,102, 73, 88, 92, 95,102, 74, 88, 92, &
      95, 38, 45, 46, 47, 49, 50, 55, 71, 74, 75, 76, &
      79, 82, 88, 92, 95,102,103, 48, 55, 76, 93, 94, &
      95, 96, 97, 98,100,102,103, 77, 88, 92, 95,102, &
      35, 37, 44, 46, 47, 72, 74, 77, 78, 79, 80, 83, &
      84, 87, 88, 90, 91, 92, 93, 94, 95, 96, 97, 98, &
      99,100,101,102,103,104, 79, 88, 92, 95,102, 80, &
      88, 92, 95,102, 80, 81, 83, 88, 92, 95,102, 80, &
      82, 83, 88, 92, 95,102, 83, 88, 92, 95,102, 80, &
      83, 84, 88, 92, 95,102, 22, 30, 35, 37, 44, 71, &
      73, 77, 82, 83, 85, 87, 88, 89, 90, 91, 92, 93, &
      95, 96, 97, 98,100,102, 30, 35, 37, 39, 40, 41, &
      44, 50, 71, 72, 73, 74, 75, 76, 77, 79, 80, 81, &
      82, 83, 84, 86, 87, 88, 90, 91, 92, 93, 94, 95, &
      96, 97, 98, 99,100,101,102,103,104, 36, 74, 77, &
      79, 82, 83, 87, 88, 92, 93, 94, 95,101,102, 18, &
      71, 73, 74, 77, 79, 80, 81, 83, 84, 88, 92, 93, &
      94, 95,102, 29, 35, 37, 44, 46, 47, 49, 53, 71, &
      73, 74, 77, 79, 81, 82, 83, 84, 87, 88, 89, 90, &
      91, 92, 93, 94, 95,101,102,103,104, 35, 37, 44, &
      73, 77, 81, 82, 83, 84, 87, 88, 90, 91, 92, 93, &
      94, 95, 99,101,102,104, 37, 44, 45, 74, 77, 79, &
      80, 82, 83, 84, 87, 88, 91, 92, 93, 94, 95, 96 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_2 = (/ &
      97, 99,100,101,102,104, 49, 50, 71, 73, 74, 77, &
      79, 80, 81, 82, 83, 84, 88, 92, 93, 94, 95, 96, &
      97, 98,100,102,103, 33, 41, 78, 79, 80, 83, 84, &
      86, 87, 88, 90, 91, 92, 93, 94, 95, 96, 97, 98, &
      99,100,101,102,103,104, 23, 24, 25, 26, 31, 33, &
      36, 41, 42, 43, 51, 70, 74, 75, 76, 78, 79, 80, &
      81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, &
      93, 94, 95, 96, 97, 98, 99,100,101,102,103,104, &
      31, 34, 43, 48, 51, 52, 55, 68, 70, 71, 73, 74, &
      75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, &
      87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, &
      99,100,101,102,103,104, 25, 52, 83, 88, 92, 93, &
      94, 95, 96, 97, 98, 99,100,101,102,103,104, 24, &
      74, 75, 76, 79, 81, 82, 83, 84, 88, 89, 90, 91, &
      92, 93, 94, 95, 96, 97, 98, 99,100,101,102,103, &
     104, 26, 80, 81, 82, 83, 84, 88, 92, 93, 94, 95, &
      96, 97, 98, 99,100,101,102,103,104, 19, 36, 40, &
      44, 71, 72, 73, 74, 77, 79, 80, 83, 84, 85, 87, &
      88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, &
     100,101,102,103,104, 23, 28, 38, 44, 46, 47, 68, &
      72, 74, 77, 79, 81, 82, 83, 84, 85, 87, 88, 89, &
      90, 91, 92, 93, 94, 95, 96, 97, 98, 99,100,101, &
     102,103,104, 29, 35, 37, 38, 44, 45, 73, 74, 77, &
      79, 80, 82, 83, 84, 87, 88, 89, 90, 91, 92, 93, &
      94, 95, 96, 97, 98, 99,100,101,102,103,104, 18, &
      19, 22, 27, 29, 30, 32, 33, 34, 35, 37, 38, 39, &
      40, 43, 44, 45, 46, 47, 48, 49, 50, 52, 53, 55, &
      68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 79, 80, &
      81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, &
      93, 94, 95, 96, 97, 98, 99,100,101,102,103,104 /)
  INTEGER, PARAMETER, DIMENSION(96) :: LU_ICOL_3 = (/ &
      27, 30, 32, 33, 34, 38, 39, 40, 41, 42, 43, 45, &
      46, 47, 49, 50, 51, 53, 68, 69, 71, 73, 74, 75, &
      76, 77, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, &
      89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99,100, &
     101,102,103,104, 22, 29, 30, 35, 37, 38, 44, 45, &
      46, 47, 48, 49, 50, 53, 55, 71, 73, 74, 76, 77, &
      79, 80, 81, 82, 83, 84, 87, 88, 89, 90, 91, 92, &
      93, 94, 95, 96, 97, 98, 99,100,101,102,103,104 /)
  INTEGER, PARAMETER, DIMENSION(1176) :: LU_ICOL = (/&
    LU_ICOL_0, LU_ICOL_1, LU_ICOL_2, LU_ICOL_3 /)

  INTEGER, PARAMETER, DIMENSION(105) :: LU_CROW = (/ &
       1,  4, 20, 30, 45, 49, 52, 57, 62, 71, 78, 84, &
      87, 89, 96,101,104,106,108,110,114,120,122,125, &
     128,131,134,137,143,145,147,150,154,157,160,162, &
     166,168,170,175,179,183,188,192,194,196,200,204, &
     209,214,217,223,230,235,243,249,259,274,289,303, &
     315,326,336,345,350,358,365,370,382,407,429,434, &
     449,454,458,476,488,493,523,528,533,540,547,552, &
     559,583,622,636,652,682,703,727,750,775,817,859, &
     876,902,922,954,988,1020,1081,1133,1177 /)

  INTEGER, PARAMETER, DIMENSION(105) :: LU_DIAG = (/ &
       1,  4, 20, 30, 45, 49, 52, 57, 62, 71, 78, 84, &
      87, 89, 96,101,104,106,108,110,115,120,122,125, &
     128,131,134,137,143,145,147,150,154,157,160,163, &
     166,168,170,175,179,183,188,192,194,198,202,206, &
     211,214,218,225,230,235,244,250,261,277,292,305, &
     317,328,338,345,353,361,367,375,388,414,429,439, &
     449,454,467,478,488,501,523,528,534,541,547,554, &
     569,604,628,646,671,693,715,740,763,806,849,867, &
     894,915,948,983,1016,1078,1131,1176,1177 /)


END MODULE saprc99_simplesom_mosaic_4bin_aq_JacobianSP

