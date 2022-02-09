! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! The Reaction Rates File
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
! File                 : cbm4_Rates.f90
! Time                 : Tue Jan 11 14:32:41 2022
! Working directory    : /jathar-scratch/yicongh/WRF_GoAmazon_gasaerochem.DEV/chem/KPP/mechanisms/cbm4
! Equation file        : cbm4.kpp
! Output root filename : cbm4
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE cbm4_Rates

  USE cbm4_Parameters
  USE cbm4_Global
  IMPLICIT NONE

CONTAINS



! Begin Rate Law Functions from KPP_HOME/util/UserRateLaws

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  User-defined Rate Law functions
!  Note: the default argument type for rate laws, as read from the equations file, is single precision
!        but all the internal calculations are performed in double precision
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~>  Arrhenius
   REAL(kind=dp) FUNCTION ARR( A0,B0,C0 )
      REAL A0,B0,C0      
      ARR =  DBLE(A0) * EXP(-DBLE(B0)/TEMP) * (TEMP/300.0_dp)**DBLE(C0)
   END FUNCTION ARR        

!~~~> Simplified Arrhenius, with two arguments
!~~~> Note: The argument B0 has a changed sign when compared to ARR
   REAL(kind=dp) FUNCTION ARR2( A0,B0 )
      REAL A0,B0           
      ARR2 =  DBLE(A0) * EXP( DBLE(B0)/TEMP )              
   END FUNCTION ARR2          

   REAL(kind=dp) FUNCTION EP2(A0,C0,A2,C2,A3,C3)
      REAL A0,C0,A2,C2,A3,C3
      REAL(dp) K0,K2,K3            
      K0 = DBLE(A0) * EXP(-DBLE(C0)/TEMP)
      K2 = DBLE(A2) * EXP(-DBLE(C2)/TEMP)
      K3 = DBLE(A3) * EXP(-DBLE(C3)/TEMP)
      K3 = K3*CFACTOR*1.0E6_dp
      EP2 = K0 + K3/(1.0_dp+K3/K2 )
   END FUNCTION EP2

   REAL(kind=dp) FUNCTION EP3(A1,C1,A2,C2) 
      REAL A1, C1, A2, C2
      REAL(dp) K1, K2      
      K1 = DBLE(A1) * EXP(-DBLE(C1)/TEMP)
      K2 = DBLE(A2) * EXP(-DBLE(C2)/TEMP)
      EP3 = K1 + K2*(1.0E6_dp*CFACTOR)
   END FUNCTION EP3 

   REAL(kind=dp) FUNCTION FALL ( A0,B0,C0,A1,B1,C1,CF)
      REAL A0,B0,C0,A1,B1,C1,CF
      REAL(dp) K0, K1     
      K0 = DBLE(A0) * EXP(-DBLE(B0)/TEMP)* (TEMP/300.0_dp)**DBLE(C0)
      K1 = DBLE(A1) * EXP(-DBLE(B1)/TEMP)* (TEMP/300.0_dp)**DBLE(C1)
      K0 = K0*CFACTOR*1.0E6_dp
      K1 = K0/K1
      FALL = (K0/(1.0_dp+K1))*   &
           DBLE(CF)**(1.0_dp/(1.0_dp+(LOG10(K1))**2))
   END FUNCTION FALL

  !---------------------------------------------------------------------------

  ELEMENTAL REAL(dp) FUNCTION k_3rd(temp,cair,k0_300K,n,kinf_300K,m,fc)

    INTRINSIC LOG10

    REAL(dp), INTENT(IN) :: temp      ! temperature [K]
    REAL(dp), INTENT(IN) :: cair      ! air concentration [molecules/cm3]
    REAL,     INTENT(IN) :: k0_300K   ! low pressure limit at 300 K
    REAL,     INTENT(IN) :: n         ! exponent for low pressure limit
    REAL,     INTENT(IN) :: kinf_300K ! high pressure limit at 300 K
    REAL,     INTENT(IN) :: m         ! exponent for high pressure limit
    REAL,     INTENT(IN) :: fc        ! broadening factor (usually fc=0.6)
    REAL                 :: zt_help, k0_T, kinf_T, k_ratio

    zt_help = 300._dp/temp
    k0_T    = k0_300K   * zt_help**(n) * cair ! k_0   at current T
    kinf_T  = kinf_300K * zt_help**(m)        ! k_inf at current T
    k_ratio = k0_T/kinf_T
    k_3rd   = k0_T/(1._dp+k_ratio)*fc**(1._dp/(1._dp+LOG10(k_ratio)**2))

  END FUNCTION k_3rd

  !---------------------------------------------------------------------------

  ELEMENTAL REAL(dp) FUNCTION k_arr (k_298,tdep,temp)
    ! Arrhenius function

    REAL,     INTENT(IN) :: k_298 ! k at T = 298.15K
    REAL,     INTENT(IN) :: tdep  ! temperature dependence
    REAL(dp), INTENT(IN) :: temp  ! temperature

    INTRINSIC EXP

    k_arr = k_298 * EXP(tdep*(1._dp/temp-3.3540E-3_dp)) ! 1/298.15=3.3540e-3

  END FUNCTION k_arr

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  End of User-defined Rate Law functions
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! End Rate Law Functions from KPP_HOME/util/UserRateLaws


! Begin INLINED Rate Law Functions


! End INLINED Rate Law Functions

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Update_SUN - update SUN light using TIME
!   Arguments :
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE Update_SUN()
      !USE cbm4_Parameters
      !USE cbm4_Global

    IMPLICIT NONE

    REAL(kind=dp) SunRise, SunSet
    REAL(kind=dp) Thour, Tlocal, Ttmp 
   
    SunRise = 4.5_dp 
    SunSet  = 19.5_dp 
    Thour = TIME/3600.0_dp 
    Tlocal = Thour - (INT(Thour)/24)*24

    IF ((Tlocal>=SunRise).AND.(Tlocal<=SunSet)) THEN
       Ttmp = (2.0*Tlocal-SunRise-SunSet)/(SunSet-SunRise)
       IF (Ttmp.GT.0) THEN
          Ttmp =  Ttmp*Ttmp
       ELSE
          Ttmp = -Ttmp*Ttmp
       END IF
       SUN = ( 1.0_dp + COS(PI*Ttmp) )/2.0_dp 
    ELSE
       SUN = 0.0_dp 
    END IF

 END SUBROUTINE Update_SUN

! End of Update_SUN function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Update_RCONST - function to update rate constants
!   Arguments :
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE Update_RCONST ( )




! Begin INLINED RCONST


! End INLINED RCONST

  RCONST(1) = (j(pj_no2))
  RCONST(2) = (j(pj_o33p))
  RCONST(3) = (j(pj_o31d))
  RCONST(4) = (j(pj_no3o)+j(pj_no3o2))
  RCONST(5) = (j(pj_hno2))
  RCONST(6) = (j(pj_h2o2))
  RCONST(7) = (j(pj_ch2or))
  RCONST(8) = (j(pj_ch2om))
  RCONST(9) = (4.6D-4*j(pj_no2))
  RCONST(10) = (9.04_dp*j(pj_ch2or))
  RCONST(11) = (9.64_dp*j(pj_ch2or))
  RCONST(12) = (ARR2(1.4D+3,1175.0_dp,TEMP))
  RCONST(13) = (ARR2(1.8D-12,-1370.0_dp,TEMP))
  RCONST(14) = (9.3D-12)
  RCONST(15) = (ARR2(1.6D-13,687.0_dp,TEMP))
  RCONST(16) = (ARR2(2.2D-13,602.0_dp,TEMP))
  RCONST(17) = (ARR2(1.2D-13,-2450.0_dp,TEMP))
  RCONST(18) = (ARR2(1.9D+8,390.0_dp,TEMP))
  RCONST(19) = (2.2D-10)
  RCONST(20) = (ARR2(1.6D-12,-940.0_dp,TEMP))
  RCONST(21) = (ARR2(1.4D-14,-580.0_dp,TEMP))
  RCONST(22) = (ARR2(1.3D-11,250.0_dp,TEMP))
  RCONST(23) = (ARR2(2.5D-14,-1230.0_dp,TEMP))
  RCONST(24) = (ARR2(5.3D-13,256.0_dp,TEMP))
  RCONST(25) = (1.3D-21)
  RCONST(26) = (ARR2(3.5D+14,-10897.0_dp,TEMP))
  RCONST(27) = (ARR2(1.8D-20,530.0_dp,TEMP))
  RCONST(28) = (4.4D-40)
  RCONST(29) = (ARR2(4.5D-13,806.0_dp,TEMP))
  RCONST(30) = (6.6D-12)
  RCONST(31) = (1.0D-20)
  RCONST(32) = (ARR2(1.0D-12,713.0_dp,TEMP))
  RCONST(33) = (ARR2(5.1D-15,1000.0_dp,TEMP))
  RCONST(34) = (ARR2(3.7D-12,240.0_dp,TEMP))
  RCONST(35) = (ARR2(1.2D-13,749.0_dp,TEMP))
  RCONST(36) = (ARR2(4.8D+13,-10121.0_dp,TEMP))
  RCONST(37) = (ARR2(1.3D-12,380.0_dp,TEMP))
  RCONST(38) = (ARR2(5.9D-14,1150.0_dp,TEMP))
  RCONST(39) = (ARR2(2.2D-38,5800.0_dp,TEMP))
  RCONST(40) = (ARR2(3.1D-12,-187.0_dp,TEMP))
  RCONST(41) = (2.2D-13)
  RCONST(42) = (1.0D-11)
  RCONST(43) = (ARR2(3.0D-11,-1550.0_dp,TEMP))
  RCONST(44) = (6.3D-16)
  RCONST(45) = (ARR2(1.2D-11,-986.0_dp,TEMP))
  RCONST(46) = (ARR2(7.0D-12,250.0_dp,TEMP))
  RCONST(47) = (2.5D-15)
  RCONST(48) = (ARR2(5.4D-12,250.0_dp,TEMP))
  RCONST(49) = (ARR2(8.0D-20,5500.0_dp,TEMP))
  RCONST(50) = (ARR2(9.4D+16,-14000.0_dp,TEMP))
  RCONST(51) = (2.0D-12)
  RCONST(52) = (6.5D-12)
  RCONST(53) = (ARR2(1.1D+2,-1710.0_dp,TEMP))
  RCONST(54) = (8.1D-13)
  RCONST(55) = (ARR2(1.0D+15,-8000.0_dp,TEMP))
  RCONST(56) = (1.6D+03)
  RCONST(57) = (1.5D-11)
  RCONST(58) = (ARR2(1.2D-11,-324.0_dp,TEMP))
  RCONST(59) = (ARR2(5.2D-12,504.0_dp,TEMP))
  RCONST(60) = (ARR2(1.4D-14,-2105.0_dp,TEMP))
  RCONST(61) = (7.7D-15)
  RCONST(62) = (ARR2(1.0D-11,-792.0_dp,TEMP))
  RCONST(63) = (ARR2(2.0D-12,411.0_dp,TEMP))
  RCONST(64) = (ARR2(1.3D-14,-2633.0_dp,TEMP))
  RCONST(65) = (ARR2(2.1D-12,322.0_dp,TEMP))
  RCONST(66) = (8.1D-12)
! RCONST(67) = constant rate coefficient
  RCONST(68) = (4.1D-11)
  RCONST(69) = (2.2D-11)
  RCONST(70) = (1.4D-11)
  RCONST(71) = (ARR2(1.7D-11,116.0_dp,TEMP))
  RCONST(72) = (3.0D-11)
  RCONST(73) = (ARR2(5.4D-17,-500.0_dp,TEMP))
  RCONST(74) = (1.70D-11)
  RCONST(75) = (1.80D-11)
  RCONST(76) = (9.6D-11)
  RCONST(77) = (1.2D-17)
  RCONST(78) = (3.2D-13)
  RCONST(79) = (8.1D-12)
  RCONST(80) = (ARR2(1.7D-14,1300.0_dp,TEMP))
  RCONST(81) = (6.8D-13)
      
END SUBROUTINE Update_RCONST

! End of Update_RCONST function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Update_PHOTO - function to update photolytical rate constants
!   Arguments :
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE Update_PHOTO ( )


   USE cbm4_Global

  RCONST(1) = (j(pj_no2))
  RCONST(2) = (j(pj_o33p))
  RCONST(3) = (j(pj_o31d))
  RCONST(4) = (j(pj_no3o)+j(pj_no3o2))
  RCONST(5) = (j(pj_hno2))
  RCONST(6) = (j(pj_h2o2))
  RCONST(7) = (j(pj_ch2or))
  RCONST(8) = (j(pj_ch2om))
  RCONST(9) = (4.6D-4*j(pj_no2))
  RCONST(10) = (9.04_dp*j(pj_ch2or))
  RCONST(11) = (9.64_dp*j(pj_ch2or))
      
END SUBROUTINE Update_PHOTO

! End of Update_PHOTO function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



END MODULE cbm4_Rates

