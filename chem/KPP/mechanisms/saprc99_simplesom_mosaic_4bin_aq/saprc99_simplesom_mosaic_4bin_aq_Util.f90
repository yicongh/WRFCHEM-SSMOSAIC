! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Auxiliary Routines File
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
! File                 : saprc99_simplesom_mosaic_4bin_aq_Util.f90
! Time                 : Tue Feb  8 14:06:57 2022
! Working directory    : /jathar-scratch/yicongh/WRF_GoAmazon_gasaerochem.DEV/chem/KPP/mechanisms/saprc99_simplesom_mosaic_4bin_aq
! Equation file        : saprc99_simplesom_mosaic_4bin_aq.kpp
! Output root filename : saprc99_simplesom_mosaic_4bin_aq
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE saprc99_simplesom_mosaic_4bin_aq_Util

  USE saprc99_simplesom_mosaic_4bin_aq_Parameters
  IMPLICIT NONE

CONTAINS



! User INLINED Utility Functions

! End INLINED Utility Functions

! Utility Functions from KPP_HOME/util/util
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! UTIL - Utility functions
!   Arguments :
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! ****************************************************************
!                            
! InitSaveData - Opens the data file for writing
!   Parameters :                                                  
!
! ****************************************************************

      SUBROUTINE InitSaveData ()

      USE saprc99_simplesom_mosaic_4bin_aq_Parameters

      open(10, file='saprc99_simplesom_mosaic_4bin_aq.dat')

      END SUBROUTINE InitSaveData

! End of InitSaveData function
! ****************************************************************

! ****************************************************************
!                            
! SaveData - Write LOOKAT species in the data file 
!   Parameters :                                                  
!
! ****************************************************************

      SUBROUTINE SaveData ()

      USE saprc99_simplesom_mosaic_4bin_aq_Global
      USE saprc99_simplesom_mosaic_4bin_aq_Monitor

      INTEGER i

      WRITE(10,999) (TIME-TSTART)/3600.D0,  &
                   (C(LOOKAT(i))/CFACTOR, i=1,NLOOKAT)
999   FORMAT(E24.16,100(1X,E24.16))

      END SUBROUTINE SaveData

! End of SaveData function
! ****************************************************************

! ****************************************************************
!                            
! CloseSaveData - Close the data file 
!   Parameters :                                                  
!
! ****************************************************************

      SUBROUTINE CloseSaveData ()

      USE saprc99_simplesom_mosaic_4bin_aq_Parameters

      CLOSE(10)

      END SUBROUTINE CloseSaveData

! End of CloseSaveData function
! ****************************************************************

! ****************************************************************
!                            
! GenerateMatlab - Generates MATLAB file to load the data file 
!   Parameters : 
!                It will have a character string to prefix each 
!                species name with.                                                 
!
! ****************************************************************

      SUBROUTINE GenerateMatlab ( PREFIX )

      USE saprc99_simplesom_mosaic_4bin_aq_Parameters
      USE saprc99_simplesom_mosaic_4bin_aq_Global
      USE saprc99_simplesom_mosaic_4bin_aq_Monitor

      
      CHARACTER(LEN=8) PREFIX 
      INTEGER i

      open(20, file='saprc99_simplesom_mosaic_4bin_aq.m')
      write(20,*) 'load saprc99_simplesom_mosaic_4bin_aq.dat;'
      write(20,990) PREFIX
990   FORMAT(A1,'c = saprc99_simplesom_mosaic_4bin_aq;')
      write(20,*) 'clear saprc99_simplesom_mosaic_4bin_aq;'
      write(20,991) PREFIX, PREFIX
991   FORMAT(A1,'t=',A1,'c(:,1);')
      write(20,992) PREFIX
992   FORMAT(A1,'c(:,1)=[];')

      do i=1,NLOOKAT
        write(20,993) PREFIX, SPC_NAMES(LOOKAT(i)), PREFIX, i
993     FORMAT(A1,A6,' = ',A1,'c(:,',I2,');')
      end do
      
      CLOSE(20)

      END SUBROUTINE GenerateMatlab

! End of GenerateMatlab function
! ****************************************************************


! End Utility Functions from KPP_HOME/util/util
! End of UTIL function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Shuffle_user2kpp - function to copy concentrations from USER to KPP
!   Arguments :
!      V_USER    - Concentration of variable species in USER's order
!      V         - Concentrations of variable species (local)
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE Shuffle_user2kpp ( V_USER, V )

! V_USER - Concentration of variable species in USER's order
  REAL(kind=dp) :: V_USER(NVAR)
! V - Concentrations of variable species (local)
  REAL(kind=dp) :: V(NVAR)

  V(92) = V_USER(1)
  V(27) = V_USER(2)
  V(93) = V_USER(3)
  V(94) = V_USER(4)
  V(95) = V_USER(5)
  V(31) = V_USER(6)
  V(33) = V_USER(7)
  V(70) = V_USER(8)
  V(43) = V_USER(9)
  V(32) = V_USER(10)
  V(1) = V_USER(11)
  V(69) = V_USER(12)
  V(86) = V_USER(13)
  V(85) = V_USER(14)
  V(34) = V_USER(15)
  V(89) = V_USER(16)
  V(72) = V_USER(17)
  V(90) = V_USER(18)
  V(2) = V_USER(19)
  V(39) = V_USER(20)
  V(30) = V_USER(21)
  V(3) = V_USER(22)
  V(4) = V_USER(23)
  V(75) = V_USER(24)
  V(68) = V_USER(25)
  V(28) = V_USER(26)
  V(48) = V_USER(27)
  V(52) = V_USER(28)
  V(82) = V_USER(29)
  V(81) = V_USER(30)
  V(84) = V_USER(31)
  V(91) = V_USER(32)
  V(49) = V_USER(33)
  V(46) = V_USER(34)
  V(47) = V_USER(35)
  V(71) = V_USER(36)
  V(80) = V_USER(37)
  V(21) = V_USER(38)
  V(15) = V_USER(39)
  V(20) = V_USER(40)
  V(16) = V_USER(41)
  V(17) = V_USER(42)
  V(22) = V_USER(43)
  V(29) = V_USER(44)
  V(50) = V_USER(45)
  V(73) = V_USER(46)
  V(35) = V_USER(47)
  V(44) = V_USER(48)
  V(37) = V_USER(49)
  V(45) = V_USER(50)
  V(38) = V_USER(51)
  V(77) = V_USER(52)
  V(83) = V_USER(53)
  V(74) = V_USER(54)
  V(79) = V_USER(55)
  V(87) = V_USER(56)
  V(51) = V_USER(57)
  V(55) = V_USER(58)
  V(23) = V_USER(59)
  V(24) = V_USER(60)
  V(25) = V_USER(61)
  V(26) = V_USER(62)
  V(5) = V_USER(63)
  V(6) = V_USER(64)
  V(97) = V_USER(65)
  V(7) = V_USER(66)
  V(8) = V_USER(67)
  V(9) = V_USER(68)
  V(88) = V_USER(69)
  V(18) = V_USER(70)
  V(102) = V_USER(71)
  V(103) = V_USER(72)
  V(99) = V_USER(73)
  V(40) = V_USER(74)
  V(53) = V_USER(75)
  V(104) = V_USER(76)
  V(78) = V_USER(77)
  V(101) = V_USER(78)
  V(41) = V_USER(79)
  V(100) = V_USER(80)
  V(96) = V_USER(81)
  V(42) = V_USER(82)
  V(76) = V_USER(83)
  V(98) = V_USER(84)
  V(36) = V_USER(85)
  V(10) = V_USER(86)
  V(11) = V_USER(87)
  V(12) = V_USER(88)
  V(13) = V_USER(89)
  V(14) = V_USER(90)
  V(57) = V_USER(91)
  V(54) = V_USER(92)
  V(56) = V_USER(93)
  V(58) = V_USER(94)
  V(59) = V_USER(95)
  V(60) = V_USER(96)
  V(61) = V_USER(97)
  V(62) = V_USER(98)
  V(63) = V_USER(99)
  V(65) = V_USER(100)
  V(66) = V_USER(101)
  V(67) = V_USER(102)
  V(64) = V_USER(103)
  V(19) = V_USER(104)
      
END SUBROUTINE Shuffle_user2kpp

! End of Shuffle_user2kpp function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Shuffle_kpp2user - function to restore concentrations from KPP to USER
!   Arguments :
!      V         - Concentrations of variable species (local)
!      V_USER    - Concentration of variable species in USER's order
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE Shuffle_kpp2user ( V, V_USER )

! V - Concentrations of variable species (local)
  REAL(kind=dp) :: V(NVAR)
! V_USER - Concentration of variable species in USER's order
  REAL(kind=dp) :: V_USER(NVAR)

  V_USER(1) = V(92)
  V_USER(2) = V(27)
  V_USER(3) = V(93)
  V_USER(4) = V(94)
  V_USER(5) = V(95)
  V_USER(6) = V(31)
  V_USER(7) = V(33)
  V_USER(8) = V(70)
  V_USER(9) = V(43)
  V_USER(10) = V(32)
  V_USER(11) = V(1)
  V_USER(12) = V(69)
  V_USER(13) = V(86)
  V_USER(14) = V(85)
  V_USER(15) = V(34)
  V_USER(16) = V(89)
  V_USER(17) = V(72)
  V_USER(18) = V(90)
  V_USER(19) = V(2)
  V_USER(20) = V(39)
  V_USER(21) = V(30)
  V_USER(22) = V(3)
  V_USER(23) = V(4)
  V_USER(24) = V(75)
  V_USER(25) = V(68)
  V_USER(26) = V(28)
  V_USER(27) = V(48)
  V_USER(28) = V(52)
  V_USER(29) = V(82)
  V_USER(30) = V(81)
  V_USER(31) = V(84)
  V_USER(32) = V(91)
  V_USER(33) = V(49)
  V_USER(34) = V(46)
  V_USER(35) = V(47)
  V_USER(36) = V(71)
  V_USER(37) = V(80)
  V_USER(38) = V(21)
  V_USER(39) = V(15)
  V_USER(40) = V(20)
  V_USER(41) = V(16)
  V_USER(42) = V(17)
  V_USER(43) = V(22)
  V_USER(44) = V(29)
  V_USER(45) = V(50)
  V_USER(46) = V(73)
  V_USER(47) = V(35)
  V_USER(48) = V(44)
  V_USER(49) = V(37)
  V_USER(50) = V(45)
  V_USER(51) = V(38)
  V_USER(52) = V(77)
  V_USER(53) = V(83)
  V_USER(54) = V(74)
  V_USER(55) = V(79)
  V_USER(56) = V(87)
  V_USER(57) = V(51)
  V_USER(58) = V(55)
  V_USER(59) = V(23)
  V_USER(60) = V(24)
  V_USER(61) = V(25)
  V_USER(62) = V(26)
  V_USER(63) = V(5)
  V_USER(64) = V(6)
  V_USER(65) = V(97)
  V_USER(66) = V(7)
  V_USER(67) = V(8)
  V_USER(68) = V(9)
  V_USER(69) = V(88)
  V_USER(70) = V(18)
  V_USER(71) = V(102)
  V_USER(72) = V(103)
  V_USER(73) = V(99)
  V_USER(74) = V(40)
  V_USER(75) = V(53)
  V_USER(76) = V(104)
  V_USER(77) = V(78)
  V_USER(78) = V(101)
  V_USER(79) = V(41)
  V_USER(80) = V(100)
  V_USER(81) = V(96)
  V_USER(82) = V(42)
  V_USER(83) = V(76)
  V_USER(84) = V(98)
  V_USER(85) = V(36)
  V_USER(86) = V(10)
  V_USER(87) = V(11)
  V_USER(88) = V(12)
  V_USER(89) = V(13)
  V_USER(90) = V(14)
  V_USER(91) = V(57)
  V_USER(92) = V(54)
  V_USER(93) = V(56)
  V_USER(94) = V(58)
  V_USER(95) = V(59)
  V_USER(96) = V(60)
  V_USER(97) = V(61)
  V_USER(98) = V(62)
  V_USER(99) = V(63)
  V_USER(100) = V(65)
  V_USER(101) = V(66)
  V_USER(102) = V(67)
  V_USER(103) = V(64)
  V_USER(104) = V(19)
      
END SUBROUTINE Shuffle_kpp2user

! End of Shuffle_kpp2user function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! GetMass - compute total mass of selected atoms
!   Arguments :
!      CL        - Concentration of all species (local)
!      Mass      - value of mass balance
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE GetMass ( CL, Mass )

! CL - Concentration of all species (local)
  REAL(kind=dp) :: CL(NSPEC)
! Mass - value of mass balance
  REAL(kind=dp) :: Mass(1)

      
END SUBROUTINE GetMass

! End of GetMass function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



END MODULE saprc99_simplesom_mosaic_4bin_aq_Util

