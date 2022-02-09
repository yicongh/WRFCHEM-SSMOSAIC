! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! The ODE Jacobian of Chemical Model File
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
! File                 : cbmz_mosaic_Jacobian.f90
! Time                 : Tue Jan 11 14:32:44 2022
! Working directory    : /jathar-scratch/yicongh/WRF_GoAmazon_gasaerochem.DEV/chem/KPP/mechanisms/cbmz_mosaic
! Equation file        : cbmz_mosaic.kpp
! Output root filename : cbmz_mosaic
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE cbmz_mosaic_Jacobian

  USE cbmz_mosaic_Parameters
  USE cbmz_mosaic_JacobianSP

  IMPLICIT NONE

CONTAINS


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Jac_SP_Vec - function for sparse multiplication: sparse Jacobian times vector
!   Arguments :
!      JVS       - sparse Jacobian of variables
!      UV        - User vector for variables
!      JUV       - Jacobian times user vector
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE Jac_SP_Vec ( JVS, UV, JUV )

! JVS - sparse Jacobian of variables
  REAL(kind=dp) :: JVS(LU_NONZERO)
! UV - User vector for variables
  REAL(kind=dp) :: UV(NVAR)
! JUV - Jacobian times user vector
  REAL(kind=dp) :: JUV(NVAR)

  JUV(1) = JVS(1)*UV(1)+JVS(2)*UV(14)+JVS(3)*UV(63)
  JUV(2) = JVS(4)*UV(2)
  JUV(3) = JVS(5)*UV(3)
  JUV(4) = JVS(6)*UV(4)+JVS(7)*UV(33)+JVS(8)*UV(44)+JVS(9)*UV(64)
  JUV(5) = JVS(10)*UV(5)+JVS(11)*UV(44)+JVS(12)*UV(45)+JVS(13)*UV(49)+JVS(14)*UV(55)+JVS(15)*UV(64)+JVS(16)*UV(65)&
             &+JVS(17)*UV(66)
  JUV(6) = JVS(18)*UV(6)+JVS(19)*UV(19)+JVS(20)*UV(21)+JVS(21)*UV(35)+JVS(22)*UV(60)+JVS(23)*UV(63)
  JUV(7) = JVS(24)*UV(7)+JVS(25)*UV(19)+JVS(26)*UV(21)+JVS(27)*UV(35)+JVS(28)*UV(60)+JVS(29)*UV(63)
  JUV(8) = JVS(30)*UV(8)+JVS(31)*UV(40)+JVS(32)*UV(63)
  JUV(9) = JVS(33)*UV(9)+JVS(34)*UV(44)+JVS(35)*UV(49)+JVS(36)*UV(60)+JVS(37)*UV(63)+JVS(38)*UV(64)
  JUV(10) = JVS(39)*UV(10)+JVS(40)*UV(24)+JVS(41)*UV(63)+JVS(42)*UV(64)
  JUV(11) = JVS(43)*UV(11)+JVS(44)*UV(24)+JVS(45)*UV(63)+JVS(46)*UV(64)
  JUV(12) = JVS(47)*UV(12)+JVS(48)*UV(25)+JVS(49)*UV(63)+JVS(50)*UV(64)
  JUV(13) = JVS(51)*UV(13)+JVS(52)*UV(25)+JVS(53)*UV(63)+JVS(54)*UV(64)
  JUV(14) = JVS(55)*UV(14)+JVS(56)*UV(63)
  JUV(15) = JVS(57)*UV(15)+JVS(58)*UV(64)
  JUV(16) = JVS(59)*UV(16)+JVS(60)*UV(63)
  JUV(17) = JVS(61)*UV(17)+JVS(62)*UV(63)+JVS(63)*UV(65)
  JUV(18) = JVS(64)*UV(18)+JVS(65)*UV(62)+JVS(66)*UV(66)
  JUV(19) = JVS(67)*UV(19)+JVS(68)*UV(63)
  JUV(20) = JVS(69)*UV(20)+JVS(70)*UV(60)+JVS(71)*UV(62)
  JUV(21) = JVS(72)*UV(21)+JVS(73)*UV(63)
  JUV(22) = JVS(74)*UV(22)+JVS(75)*UV(44)+JVS(76)*UV(49)+JVS(77)*UV(63)+JVS(78)*UV(64)
  JUV(23) = JVS(79)*UV(23)+JVS(80)*UV(35)+JVS(81)*UV(60)+JVS(82)*UV(62)+JVS(83)*UV(63)
  JUV(24) = JVS(84)*UV(24)+JVS(85)*UV(60)+JVS(86)*UV(63)+JVS(87)*UV(64)
  JUV(25) = JVS(88)*UV(25)+JVS(89)*UV(60)+JVS(90)*UV(63)+JVS(91)*UV(64)
  JUV(26) = JVS(92)*UV(26)+JVS(93)*UV(62)+JVS(94)*UV(63)+JVS(95)*UV(65)
  JUV(27) = JVS(96)*UV(19)+JVS(97)*UV(21)+JVS(98)*UV(27)+JVS(99)*UV(61)+JVS(100)*UV(63)
  JUV(28) = JVS(101)*UV(28)+JVS(102)*UV(44)+JVS(103)*UV(49)+JVS(104)*UV(51)+JVS(105)*UV(63)+JVS(106)*UV(64)
  JUV(29) = JVS(107)*UV(29)+JVS(108)*UV(40)+JVS(109)*UV(44)+JVS(110)*UV(49)+JVS(111)*UV(52)+JVS(112)*UV(56)+JVS(113)&
              &*UV(58)+JVS(114)*UV(59)+JVS(115)*UV(60)+JVS(116)*UV(61)+JVS(117)*UV(63)+JVS(118)*UV(64)
  JUV(30) = JVS(119)*UV(30)+JVS(120)*UV(51)+JVS(121)*UV(63)+JVS(122)*UV(65)
  JUV(31) = JVS(123)*UV(31)+JVS(124)*UV(61)+JVS(125)*UV(62)+JVS(126)*UV(63)+JVS(127)*UV(65)
  JUV(32) = JVS(128)*UV(32)+JVS(129)*UV(54)+JVS(130)*UV(63)+JVS(131)*UV(65)
  JUV(33) = JVS(132)*UV(33)+JVS(133)*UV(63)+JVS(134)*UV(64)
  JUV(34) = JVS(135)*UV(34)+JVS(136)*UV(44)+JVS(137)*UV(49)+JVS(138)*UV(54)+JVS(139)*UV(63)+JVS(140)*UV(64)
  JUV(35) = JVS(141)*UV(19)+JVS(142)*UV(21)+JVS(143)*UV(35)+JVS(144)*UV(60)+JVS(145)*UV(63)
  JUV(36) = JVS(146)*UV(15)+JVS(147)*UV(36)+JVS(148)*UV(60)+JVS(149)*UV(61)+JVS(150)*UV(62)+JVS(151)*UV(64)
  JUV(37) = JVS(152)*UV(20)+JVS(153)*UV(35)+JVS(154)*UV(37)+JVS(155)*UV(46)+JVS(156)*UV(50)+JVS(157)*UV(53)+JVS(158)&
              &*UV(55)+JVS(159)*UV(60)+JVS(160)*UV(62)+JVS(161)*UV(63)+JVS(162)*UV(65)
  JUV(38) = JVS(163)*UV(33)+JVS(164)*UV(38)+JVS(165)*UV(41)+JVS(166)*UV(43)+JVS(167)*UV(44)+JVS(168)*UV(45)+JVS(169)&
              &*UV(46)+JVS(170)*UV(49)+JVS(171)*UV(50)+JVS(172)*UV(53)+JVS(173)*UV(55)+JVS(174)*UV(60)+JVS(175)*UV(61)&
              &+JVS(176)*UV(63)+JVS(177)*UV(64)
  JUV(39) = JVS(178)*UV(39)+JVS(179)*UV(45)+JVS(180)*UV(60)+JVS(181)*UV(61)+JVS(182)*UV(65)
  JUV(40) = JVS(183)*UV(21)+JVS(184)*UV(29)+JVS(185)*UV(39)+JVS(186)*UV(40)+JVS(187)*UV(42)+JVS(192)*UV(55)+JVS(196)&
              &*UV(60)+JVS(197)*UV(61)+JVS(198)*UV(63)+JVS(200)*UV(65)
  JUV(41) = JVS(201)*UV(27)+JVS(202)*UV(35)+JVS(203)*UV(41)+JVS(205)*UV(61)+JVS(206)*UV(63)+JVS(207)*UV(64)
  JUV(42) = JVS(208)*UV(42)+JVS(209)*UV(45)+JVS(210)*UV(61)+JVS(211)*UV(63)+JVS(212)*UV(65)
  JUV(43) = JVS(213)*UV(43)+JVS(214)*UV(55)+JVS(215)*UV(61)+JVS(216)*UV(63)+JVS(217)*UV(65)
  JUV(44) = JVS(218)*UV(44)+JVS(219)*UV(60)+JVS(220)*UV(63)+JVS(221)*UV(64)
  JUV(45) = JVS(222)*UV(45)+JVS(223)*UV(60)+JVS(224)*UV(63)+JVS(225)*UV(64)
  JUV(46) = JVS(226)*UV(32)+JVS(227)*UV(33)+JVS(228)*UV(34)+JVS(229)*UV(41)+JVS(230)*UV(42)+JVS(231)*UV(43)+JVS(232)&
              &*UV(44)+JVS(233)*UV(45)+JVS(234)*UV(46)+JVS(236)*UV(52)+JVS(237)*UV(54)+JVS(238)*UV(55)+JVS(239)*UV(57)&
              &+JVS(240)*UV(60)+JVS(241)*UV(61)+JVS(242)*UV(63)+JVS(243)*UV(64)
  JUV(47) = JVS(245)*UV(19)+JVS(246)*UV(21)+JVS(247)*UV(33)+JVS(248)*UV(35)+JVS(249)*UV(41)+JVS(250)*UV(44)+JVS(251)&
              &*UV(45)+JVS(252)*UV(47)+JVS(253)*UV(49)+JVS(254)*UV(50)+JVS(255)*UV(55)+JVS(256)*UV(56)+JVS(257)*UV(58)&
              &+JVS(258)*UV(59)+JVS(259)*UV(60)+JVS(260)*UV(61)+JVS(261)*UV(63)+JVS(262)*UV(64)+JVS(263)*UV(65)
  JUV(48) = JVS(264)*UV(43)+JVS(265)*UV(48)+JVS(266)*UV(49)+JVS(267)*UV(55)+JVS(268)*UV(56)+JVS(269)*UV(57)+JVS(270)&
              &*UV(58)+JVS(271)*UV(59)+JVS(272)*UV(60)+JVS(273)*UV(61)+JVS(274)*UV(63)+JVS(275)*UV(64)
  JUV(49) = JVS(277)*UV(49)+JVS(278)*UV(60)+JVS(279)*UV(63)+JVS(280)*UV(64)
  JUV(50) = JVS(281)*UV(21)+JVS(282)*UV(41)+JVS(283)*UV(43)+JVS(284)*UV(44)+JVS(285)*UV(49)+JVS(286)*UV(50)+JVS(287)&
              &*UV(55)+JVS(288)*UV(56)+JVS(289)*UV(57)+JVS(290)*UV(60)+JVS(291)*UV(61)+JVS(292)*UV(63)+JVS(293)*UV(64)
  JUV(51) = JVS(295)*UV(28)+JVS(296)*UV(30)+JVS(297)*UV(44)+JVS(298)*UV(49)+JVS(299)*UV(51)+JVS(300)*UV(56)+JVS(301)&
              &*UV(58)+JVS(302)*UV(59)+JVS(303)*UV(60)+JVS(304)*UV(61)+JVS(305)*UV(63)+JVS(306)*UV(64)+JVS(307)*UV(65)
  JUV(52) = JVS(308)*UV(44)+JVS(309)*UV(49)+JVS(310)*UV(52)+JVS(311)*UV(59)+JVS(312)*UV(60)+JVS(313)*UV(61)+JVS(314)&
              &*UV(63)+JVS(316)*UV(65)
  JUV(53) = JVS(317)*UV(16)+JVS(318)*UV(30)+JVS(319)*UV(33)+JVS(320)*UV(39)+JVS(321)*UV(41)+JVS(322)*UV(43)+JVS(323)&
              &*UV(44)+JVS(324)*UV(45)+JVS(325)*UV(49)+JVS(326)*UV(51)+JVS(327)*UV(52)+JVS(328)*UV(53)+JVS(329)*UV(55)&
              &+JVS(330)*UV(56)+JVS(331)*UV(58)+JVS(332)*UV(59)+JVS(333)*UV(60)+JVS(334)*UV(61)+JVS(335)*UV(63)+JVS(336)&
              &*UV(64)
  JUV(54) = JVS(338)*UV(22)+JVS(339)*UV(32)+JVS(340)*UV(44)+JVS(341)*UV(48)+JVS(342)*UV(49)+JVS(343)*UV(53)+JVS(344)&
              &*UV(54)+JVS(345)*UV(55)+JVS(350)*UV(60)+JVS(351)*UV(61)+JVS(352)*UV(63)+JVS(353)*UV(64)+JVS(354)*UV(65)&
              &+JVS(355)*UV(66)
  JUV(55) = JVS(356)*UV(39)+JVS(357)*UV(42)+JVS(358)*UV(45)+JVS(359)*UV(55)+JVS(360)*UV(60)+JVS(361)*UV(61)+JVS(362)&
              &*UV(63)+JVS(363)*UV(64)
  JUV(56) = JVS(365)*UV(42)+JVS(366)*UV(43)+JVS(369)*UV(56)+JVS(370)*UV(57)+JVS(371)*UV(58)+JVS(374)*UV(63)+JVS(376)&
              &*UV(65)
  JUV(57) = JVS(377)*UV(48)+JVS(378)*UV(49)+JVS(381)*UV(57)+JVS(384)*UV(60)+JVS(385)*UV(61)+JVS(386)*UV(63)+JVS(387)&
              &*UV(64)+JVS(388)*UV(65)
  JUV(58) = JVS(389)*UV(40)+JVS(391)*UV(44)+JVS(393)*UV(49)+JVS(396)*UV(56)+JVS(398)*UV(58)+JVS(400)*UV(60)+JVS(401)&
              &*UV(61)+JVS(402)*UV(63)+JVS(403)*UV(64)+JVS(404)*UV(65)
  JUV(59) = JVS(405)*UV(23)+JVS(406)*UV(27)+JVS(408)*UV(39)+JVS(409)*UV(42)+JVS(411)*UV(52)+JVS(412)*UV(55)+JVS(413)&
              &*UV(58)+JVS(414)*UV(59)+JVS(415)*UV(60)+JVS(416)*UV(61)+JVS(417)*UV(62)+JVS(418)*UV(63)+JVS(420)*UV(65)
  JUV(60) = JVS(421)*UV(20)+JVS(422)*UV(24)+JVS(423)*UV(25)+JVS(424)*UV(35)+JVS(425)*UV(36)+JVS(426)*UV(37)+JVS(427)&
              &*UV(44)+JVS(428)*UV(45)+JVS(429)*UV(46)+JVS(430)*UV(47)+JVS(431)*UV(49)+JVS(432)*UV(50)+JVS(433)*UV(51)&
              &+JVS(434)*UV(52)+JVS(435)*UV(53)+JVS(436)*UV(54)+JVS(437)*UV(55)+JVS(439)*UV(57)+JVS(440)*UV(58)+JVS(442)&
              &*UV(60)+JVS(443)*UV(61)+JVS(444)*UV(62)+JVS(445)*UV(63)+JVS(446)*UV(64)+JVS(447)*UV(65)+JVS(448)*UV(66)
  JUV(61) = JVS(449)*UV(27)+JVS(450)*UV(31)+JVS(451)*UV(36)+JVS(452)*UV(39)+JVS(453)*UV(42)+JVS(454)*UV(43)+JVS(456)&
              &*UV(47)+JVS(459)*UV(51)+JVS(460)*UV(52)+JVS(461)*UV(54)+JVS(464)*UV(57)+JVS(465)*UV(58)+JVS(467)*UV(60)&
              &+JVS(468)*UV(61)+JVS(469)*UV(62)+JVS(470)*UV(63)+JVS(471)*UV(64)+JVS(472)*UV(65)+JVS(473)*UV(66)
  JUV(62) = JVS(474)*UV(18)+JVS(475)*UV(20)+JVS(476)*UV(23)+JVS(477)*UV(26)+JVS(478)*UV(27)+JVS(479)*UV(31)+JVS(481)&
              &*UV(36)+JVS(482)*UV(37)+JVS(483)*UV(39)+JVS(484)*UV(42)+JVS(485)*UV(43)+JVS(488)*UV(47)+JVS(491)*UV(51)&
              &+JVS(492)*UV(52)+JVS(494)*UV(54)+JVS(497)*UV(57)+JVS(498)*UV(58)+JVS(499)*UV(59)+JVS(500)*UV(60)+JVS(501)&
              &*UV(61)+JVS(502)*UV(62)+JVS(503)*UV(63)+JVS(504)*UV(64)+JVS(505)*UV(65)+JVS(506)*UV(66)
  JUV(63) = JVS(507)*UV(14)+JVS(508)*UV(15)+JVS(509)*UV(16)+JVS(510)*UV(17)+JVS(511)*UV(19)+JVS(512)*UV(21)+JVS(513)&
              &*UV(22)+JVS(514)*UV(24)+JVS(515)*UV(25)+JVS(516)*UV(26)+JVS(517)*UV(28)+JVS(518)*UV(30)+JVS(519)*UV(31)&
              &+JVS(520)*UV(32)+JVS(521)*UV(33)+JVS(522)*UV(34)+JVS(523)*UV(35)+JVS(524)*UV(37)+JVS(525)*UV(38)+JVS(526)&
              &*UV(40)+JVS(527)*UV(41)+JVS(530)*UV(44)+JVS(531)*UV(45)+JVS(532)*UV(46)+JVS(533)*UV(48)+JVS(534)*UV(49)&
              &+JVS(535)*UV(50)+JVS(538)*UV(53)+JVS(540)*UV(55)+JVS(541)*UV(56)+JVS(544)*UV(59)+JVS(545)*UV(60)+JVS(546)&
              &*UV(61)+JVS(547)*UV(62)+JVS(548)*UV(63)+JVS(549)*UV(64)+JVS(550)*UV(65)
  JUV(64) = JVS(552)*UV(24)+JVS(553)*UV(25)+JVS(554)*UV(33)+JVS(555)*UV(36)+JVS(556)*UV(41)+JVS(557)*UV(44)+JVS(558)&
              &*UV(45)+JVS(559)*UV(49)+JVS(560)*UV(55)+JVS(562)*UV(61)+JVS(563)*UV(62)+JVS(564)*UV(63)+JVS(565)*UV(64)&
              &+JVS(566)*UV(65)+JVS(567)*UV(66)
  JUV(65) = JVS(568)*UV(14)+JVS(569)*UV(16)+JVS(570)*UV(17)+JVS(571)*UV(19)+JVS(572)*UV(21)+JVS(573)*UV(26)+JVS(574)&
              &*UV(27)+JVS(575)*UV(30)+JVS(576)*UV(32)+JVS(577)*UV(33)+JVS(578)*UV(34)+JVS(579)*UV(35)+JVS(580)*UV(38)&
              &+JVS(581)*UV(39)+JVS(582)*UV(41)+JVS(583)*UV(42)+JVS(584)*UV(43)+JVS(585)*UV(44)+JVS(586)*UV(45)+JVS(587)&
              &*UV(46)+JVS(588)*UV(47)+JVS(589)*UV(49)+JVS(590)*UV(50)+JVS(591)*UV(51)+JVS(592)*UV(52)+JVS(593)*UV(53)&
              &+JVS(594)*UV(54)+JVS(595)*UV(55)+JVS(596)*UV(56)+JVS(597)*UV(57)+JVS(598)*UV(58)+JVS(599)*UV(59)+JVS(600)&
              &*UV(60)+JVS(601)*UV(61)+JVS(602)*UV(62)+JVS(603)*UV(63)+JVS(604)*UV(64)+JVS(605)*UV(65)+JVS(606)*UV(66)
  JUV(66) = JVS(607)*UV(18)+JVS(608)*UV(41)+JVS(609)*UV(44)+JVS(610)*UV(45)+JVS(611)*UV(48)+JVS(612)*UV(49)+JVS(613)&
              &*UV(50)+JVS(614)*UV(53)+JVS(615)*UV(55)+JVS(617)*UV(57)+JVS(620)*UV(60)+JVS(621)*UV(61)+JVS(622)*UV(62)&
              &+JVS(623)*UV(63)+JVS(624)*UV(64)+JVS(625)*UV(65)+JVS(626)*UV(66)
      
END SUBROUTINE Jac_SP_Vec

! End of Jac_SP_Vec function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! JacTR_SP_Vec - sparse multiplication: sparse Jacobian transposed times vector
!   Arguments :
!      JVS       - sparse Jacobian of variables
!      UV        - User vector for variables
!      JTUV      - Jacobian transposed times user vector
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE JacTR_SP_Vec ( JVS, UV, JTUV )

! JVS - sparse Jacobian of variables
  REAL(kind=dp) :: JVS(LU_NONZERO)
! UV - User vector for variables
  REAL(kind=dp) :: UV(NVAR)
! JTUV - Jacobian transposed times user vector
  REAL(kind=dp) :: JTUV(NVAR)

  JTUV(1) = JVS(1)*UV(1)
  JTUV(2) = JVS(4)*UV(2)
  JTUV(3) = JVS(5)*UV(3)
  JTUV(4) = JVS(6)*UV(4)
  JTUV(5) = JVS(10)*UV(5)
  JTUV(6) = JVS(18)*UV(6)
  JTUV(7) = JVS(24)*UV(7)
  JTUV(8) = JVS(30)*UV(8)
  JTUV(9) = JVS(33)*UV(9)
  JTUV(10) = JVS(39)*UV(10)
  JTUV(11) = JVS(43)*UV(11)
  JTUV(12) = JVS(47)*UV(12)
  JTUV(13) = JVS(51)*UV(13)
  JTUV(14) = JVS(2)*UV(1)+JVS(55)*UV(14)+JVS(507)*UV(63)+JVS(568)*UV(65)
  JTUV(15) = JVS(57)*UV(15)+JVS(146)*UV(36)+JVS(508)*UV(63)
  JTUV(16) = JVS(59)*UV(16)+JVS(317)*UV(53)+JVS(509)*UV(63)+JVS(569)*UV(65)
  JTUV(17) = JVS(61)*UV(17)+JVS(510)*UV(63)+JVS(570)*UV(65)
  JTUV(18) = JVS(64)*UV(18)+JVS(474)*UV(62)+JVS(607)*UV(66)
  JTUV(19) = JVS(19)*UV(6)+JVS(25)*UV(7)+JVS(67)*UV(19)+JVS(96)*UV(27)+JVS(141)*UV(35)+JVS(245)*UV(47)+JVS(511)*UV(63)&
               &+JVS(571)*UV(65)
  JTUV(20) = JVS(69)*UV(20)+JVS(152)*UV(37)+JVS(421)*UV(60)+JVS(475)*UV(62)
  JTUV(21) = JVS(20)*UV(6)+JVS(26)*UV(7)+JVS(72)*UV(21)+JVS(97)*UV(27)+JVS(142)*UV(35)+JVS(183)*UV(40)+JVS(246)*UV(47)&
               &+JVS(281)*UV(50)+JVS(512)*UV(63)+JVS(572)*UV(65)
  JTUV(22) = JVS(74)*UV(22)+JVS(338)*UV(54)+JVS(513)*UV(63)
  JTUV(23) = JVS(79)*UV(23)+JVS(405)*UV(59)+JVS(476)*UV(62)
  JTUV(24) = JVS(40)*UV(10)+JVS(44)*UV(11)+JVS(84)*UV(24)+JVS(422)*UV(60)+JVS(514)*UV(63)+JVS(552)*UV(64)
  JTUV(25) = JVS(48)*UV(12)+JVS(52)*UV(13)+JVS(88)*UV(25)+JVS(423)*UV(60)+JVS(515)*UV(63)+JVS(553)*UV(64)
  JTUV(26) = JVS(92)*UV(26)+JVS(477)*UV(62)+JVS(516)*UV(63)+JVS(573)*UV(65)
  JTUV(27) = JVS(98)*UV(27)+JVS(201)*UV(41)+JVS(406)*UV(59)+JVS(449)*UV(61)+JVS(478)*UV(62)+JVS(574)*UV(65)
  JTUV(28) = JVS(101)*UV(28)+JVS(295)*UV(51)+JVS(517)*UV(63)
  JTUV(29) = JVS(107)*UV(29)+JVS(184)*UV(40)
  JTUV(30) = JVS(119)*UV(30)+JVS(296)*UV(51)+JVS(318)*UV(53)+JVS(518)*UV(63)+JVS(575)*UV(65)
  JTUV(31) = JVS(123)*UV(31)+JVS(450)*UV(61)+JVS(479)*UV(62)+JVS(519)*UV(63)
  JTUV(32) = JVS(128)*UV(32)+JVS(226)*UV(46)+JVS(339)*UV(54)+JVS(520)*UV(63)+JVS(576)*UV(65)
  JTUV(33) = JVS(7)*UV(4)+JVS(132)*UV(33)+JVS(163)*UV(38)+JVS(227)*UV(46)+JVS(247)*UV(47)+JVS(319)*UV(53)+JVS(521)&
               &*UV(63)+JVS(554)*UV(64)+JVS(577)*UV(65)
  JTUV(34) = JVS(135)*UV(34)+JVS(228)*UV(46)+JVS(522)*UV(63)+JVS(578)*UV(65)
  JTUV(35) = JVS(21)*UV(6)+JVS(27)*UV(7)+JVS(80)*UV(23)+JVS(143)*UV(35)+JVS(153)*UV(37)+JVS(202)*UV(41)+JVS(248)*UV(47)&
               &+JVS(424)*UV(60)+JVS(523)*UV(63)+JVS(579)*UV(65)
  JTUV(36) = JVS(147)*UV(36)+JVS(425)*UV(60)+JVS(451)*UV(61)+JVS(481)*UV(62)+JVS(555)*UV(64)
  JTUV(37) = JVS(154)*UV(37)+JVS(426)*UV(60)+JVS(482)*UV(62)+JVS(524)*UV(63)
  JTUV(38) = JVS(164)*UV(38)+JVS(525)*UV(63)+JVS(580)*UV(65)
  JTUV(39) = JVS(178)*UV(39)+JVS(185)*UV(40)+JVS(320)*UV(53)+JVS(356)*UV(55)+JVS(408)*UV(59)+JVS(452)*UV(61)+JVS(483)&
               &*UV(62)+JVS(581)*UV(65)
  JTUV(40) = JVS(31)*UV(8)+JVS(108)*UV(29)+JVS(186)*UV(40)+JVS(389)*UV(58)+JVS(526)*UV(63)
  JTUV(41) = JVS(165)*UV(38)+JVS(203)*UV(41)+JVS(229)*UV(46)+JVS(249)*UV(47)+JVS(282)*UV(50)+JVS(321)*UV(53)+JVS(527)&
               &*UV(63)+JVS(556)*UV(64)+JVS(582)*UV(65)+JVS(608)*UV(66)
  JTUV(42) = JVS(187)*UV(40)+JVS(208)*UV(42)+JVS(230)*UV(46)+JVS(357)*UV(55)+JVS(365)*UV(56)+JVS(409)*UV(59)+JVS(453)&
               &*UV(61)+JVS(484)*UV(62)+JVS(583)*UV(65)
  JTUV(43) = JVS(166)*UV(38)+JVS(213)*UV(43)+JVS(231)*UV(46)+JVS(264)*UV(48)+JVS(283)*UV(50)+JVS(322)*UV(53)+JVS(366)&
               &*UV(56)+JVS(454)*UV(61)+JVS(485)*UV(62)+JVS(584)*UV(65)
  JTUV(44) = JVS(8)*UV(4)+JVS(11)*UV(5)+JVS(34)*UV(9)+JVS(75)*UV(22)+JVS(102)*UV(28)+JVS(109)*UV(29)+JVS(136)*UV(34)&
               &+JVS(167)*UV(38)+JVS(218)*UV(44)+JVS(232)*UV(46)+JVS(250)*UV(47)+JVS(284)*UV(50)+JVS(297)*UV(51)+JVS(308)&
               &*UV(52)+JVS(323)*UV(53)+JVS(340)*UV(54)+JVS(391)*UV(58)+JVS(427)*UV(60)+JVS(530)*UV(63)+JVS(557)*UV(64)&
               &+JVS(585)*UV(65)+JVS(609)*UV(66)
  JTUV(45) = JVS(12)*UV(5)+JVS(168)*UV(38)+JVS(179)*UV(39)+JVS(209)*UV(42)+JVS(222)*UV(45)+JVS(233)*UV(46)+JVS(251)&
               &*UV(47)+JVS(324)*UV(53)+JVS(358)*UV(55)+JVS(428)*UV(60)+JVS(531)*UV(63)+JVS(558)*UV(64)+JVS(586)*UV(65)&
               &+JVS(610)*UV(66)
  JTUV(46) = JVS(155)*UV(37)+JVS(169)*UV(38)+JVS(234)*UV(46)+JVS(429)*UV(60)+JVS(532)*UV(63)+JVS(587)*UV(65)
  JTUV(47) = JVS(252)*UV(47)+JVS(430)*UV(60)+JVS(456)*UV(61)+JVS(488)*UV(62)+JVS(588)*UV(65)
  JTUV(48) = JVS(265)*UV(48)+JVS(341)*UV(54)+JVS(377)*UV(57)+JVS(533)*UV(63)+JVS(611)*UV(66)
  JTUV(49) = JVS(13)*UV(5)+JVS(35)*UV(9)+JVS(76)*UV(22)+JVS(103)*UV(28)+JVS(110)*UV(29)+JVS(137)*UV(34)+JVS(170)*UV(38)&
               &+JVS(253)*UV(47)+JVS(266)*UV(48)+JVS(277)*UV(49)+JVS(285)*UV(50)+JVS(298)*UV(51)+JVS(309)*UV(52)+JVS(325)&
               &*UV(53)+JVS(342)*UV(54)+JVS(378)*UV(57)+JVS(393)*UV(58)+JVS(431)*UV(60)+JVS(534)*UV(63)+JVS(559)*UV(64)&
               &+JVS(589)*UV(65)+JVS(612)*UV(66)
  JTUV(50) = JVS(156)*UV(37)+JVS(171)*UV(38)+JVS(254)*UV(47)+JVS(286)*UV(50)+JVS(432)*UV(60)+JVS(535)*UV(63)+JVS(590)&
               &*UV(65)+JVS(613)*UV(66)
  JTUV(51) = JVS(104)*UV(28)+JVS(120)*UV(30)+JVS(299)*UV(51)+JVS(326)*UV(53)+JVS(433)*UV(60)+JVS(459)*UV(61)+JVS(491)&
               &*UV(62)+JVS(591)*UV(65)
  JTUV(52) = JVS(111)*UV(29)+JVS(236)*UV(46)+JVS(310)*UV(52)+JVS(327)*UV(53)+JVS(411)*UV(59)+JVS(434)*UV(60)+JVS(460)&
               &*UV(61)+JVS(492)*UV(62)+JVS(592)*UV(65)
  JTUV(53) = JVS(157)*UV(37)+JVS(172)*UV(38)+JVS(328)*UV(53)+JVS(343)*UV(54)+JVS(435)*UV(60)+JVS(538)*UV(63)+JVS(593)&
               &*UV(65)+JVS(614)*UV(66)
  JTUV(54) = JVS(129)*UV(32)+JVS(138)*UV(34)+JVS(237)*UV(46)+JVS(344)*UV(54)+JVS(436)*UV(60)+JVS(461)*UV(61)+JVS(494)&
               &*UV(62)+JVS(594)*UV(65)
  JTUV(55) = JVS(14)*UV(5)+JVS(158)*UV(37)+JVS(173)*UV(38)+JVS(192)*UV(40)+JVS(214)*UV(43)+JVS(238)*UV(46)+JVS(255)&
               &*UV(47)+JVS(267)*UV(48)+JVS(287)*UV(50)+JVS(329)*UV(53)+JVS(345)*UV(54)+JVS(359)*UV(55)+JVS(412)*UV(59)&
               &+JVS(437)*UV(60)+JVS(540)*UV(63)+JVS(560)*UV(64)+JVS(595)*UV(65)+JVS(615)*UV(66)
  JTUV(56) = JVS(112)*UV(29)+JVS(256)*UV(47)+JVS(268)*UV(48)+JVS(288)*UV(50)+JVS(300)*UV(51)+JVS(330)*UV(53)+JVS(369)&
               &*UV(56)+JVS(396)*UV(58)+JVS(541)*UV(63)+JVS(596)*UV(65)
  JTUV(57) = JVS(239)*UV(46)+JVS(269)*UV(48)+JVS(289)*UV(50)+JVS(370)*UV(56)+JVS(381)*UV(57)+JVS(439)*UV(60)+JVS(464)&
               &*UV(61)+JVS(497)*UV(62)+JVS(597)*UV(65)+JVS(617)*UV(66)
  JTUV(58) = JVS(113)*UV(29)+JVS(257)*UV(47)+JVS(270)*UV(48)+JVS(301)*UV(51)+JVS(331)*UV(53)+JVS(371)*UV(56)+JVS(398)&
               &*UV(58)+JVS(413)*UV(59)+JVS(440)*UV(60)+JVS(465)*UV(61)+JVS(498)*UV(62)+JVS(598)*UV(65)
  JTUV(59) = JVS(114)*UV(29)+JVS(258)*UV(47)+JVS(271)*UV(48)+JVS(302)*UV(51)+JVS(311)*UV(52)+JVS(332)*UV(53)+JVS(414)&
               &*UV(59)+JVS(499)*UV(62)+JVS(544)*UV(63)+JVS(599)*UV(65)
  JTUV(60) = JVS(22)*UV(6)+JVS(28)*UV(7)+JVS(36)*UV(9)+JVS(70)*UV(20)+JVS(81)*UV(23)+JVS(85)*UV(24)+JVS(89)*UV(25)&
               &+JVS(115)*UV(29)+JVS(144)*UV(35)+JVS(148)*UV(36)+JVS(159)*UV(37)+JVS(174)*UV(38)+JVS(180)*UV(39)+JVS(196)&
               &*UV(40)+JVS(219)*UV(44)+JVS(223)*UV(45)+JVS(240)*UV(46)+JVS(259)*UV(47)+JVS(272)*UV(48)+JVS(278)*UV(49)&
               &+JVS(290)*UV(50)+JVS(303)*UV(51)+JVS(312)*UV(52)+JVS(333)*UV(53)+JVS(350)*UV(54)+JVS(360)*UV(55)+JVS(384)&
               &*UV(57)+JVS(400)*UV(58)+JVS(415)*UV(59)+JVS(442)*UV(60)+JVS(467)*UV(61)+JVS(500)*UV(62)+JVS(545)*UV(63)&
               &+JVS(600)*UV(65)+JVS(620)*UV(66)
  JTUV(61) = JVS(99)*UV(27)+JVS(116)*UV(29)+JVS(124)*UV(31)+JVS(149)*UV(36)+JVS(175)*UV(38)+JVS(181)*UV(39)+JVS(197)&
               &*UV(40)+JVS(205)*UV(41)+JVS(210)*UV(42)+JVS(215)*UV(43)+JVS(241)*UV(46)+JVS(260)*UV(47)+JVS(273)*UV(48)&
               &+JVS(291)*UV(50)+JVS(304)*UV(51)+JVS(313)*UV(52)+JVS(334)*UV(53)+JVS(351)*UV(54)+JVS(361)*UV(55)+JVS(385)&
               &*UV(57)+JVS(401)*UV(58)+JVS(416)*UV(59)+JVS(443)*UV(60)+JVS(468)*UV(61)+JVS(501)*UV(62)+JVS(546)*UV(63)&
               &+JVS(562)*UV(64)+JVS(601)*UV(65)+JVS(621)*UV(66)
  JTUV(62) = JVS(65)*UV(18)+JVS(71)*UV(20)+JVS(82)*UV(23)+JVS(93)*UV(26)+JVS(125)*UV(31)+JVS(150)*UV(36)+JVS(160)*UV(37)&
               &+JVS(417)*UV(59)+JVS(444)*UV(60)+JVS(469)*UV(61)+JVS(502)*UV(62)+JVS(547)*UV(63)+JVS(563)*UV(64)+JVS(602)&
               &*UV(65)+JVS(622)*UV(66)
  JTUV(63) = JVS(3)*UV(1)+JVS(23)*UV(6)+JVS(29)*UV(7)+JVS(32)*UV(8)+JVS(37)*UV(9)+JVS(41)*UV(10)+JVS(45)*UV(11)+JVS(49)&
               &*UV(12)+JVS(53)*UV(13)+JVS(56)*UV(14)+JVS(60)*UV(16)+JVS(62)*UV(17)+JVS(68)*UV(19)+JVS(73)*UV(21)+JVS(77)&
               &*UV(22)+JVS(83)*UV(23)+JVS(86)*UV(24)+JVS(90)*UV(25)+JVS(94)*UV(26)+JVS(100)*UV(27)+JVS(105)*UV(28)+JVS(117)&
               &*UV(29)+JVS(121)*UV(30)+JVS(126)*UV(31)+JVS(130)*UV(32)+JVS(133)*UV(33)+JVS(139)*UV(34)+JVS(145)*UV(35)&
               &+JVS(161)*UV(37)+JVS(176)*UV(38)+JVS(198)*UV(40)+JVS(206)*UV(41)+JVS(211)*UV(42)+JVS(216)*UV(43)+JVS(220)&
               &*UV(44)+JVS(224)*UV(45)+JVS(242)*UV(46)+JVS(261)*UV(47)+JVS(274)*UV(48)+JVS(279)*UV(49)+JVS(292)*UV(50)&
               &+JVS(305)*UV(51)+JVS(314)*UV(52)+JVS(335)*UV(53)+JVS(352)*UV(54)+JVS(362)*UV(55)+JVS(374)*UV(56)+JVS(386)&
               &*UV(57)+JVS(402)*UV(58)+JVS(418)*UV(59)+JVS(445)*UV(60)+JVS(470)*UV(61)+JVS(503)*UV(62)+JVS(548)*UV(63)&
               &+JVS(564)*UV(64)+JVS(603)*UV(65)+JVS(623)*UV(66)
  JTUV(64) = JVS(9)*UV(4)+JVS(15)*UV(5)+JVS(38)*UV(9)+JVS(42)*UV(10)+JVS(46)*UV(11)+JVS(50)*UV(12)+JVS(54)*UV(13)&
               &+JVS(58)*UV(15)+JVS(78)*UV(22)+JVS(87)*UV(24)+JVS(91)*UV(25)+JVS(106)*UV(28)+JVS(118)*UV(29)+JVS(134)*UV(33)&
               &+JVS(140)*UV(34)+JVS(151)*UV(36)+JVS(177)*UV(38)+JVS(207)*UV(41)+JVS(221)*UV(44)+JVS(225)*UV(45)+JVS(243)&
               &*UV(46)+JVS(262)*UV(47)+JVS(275)*UV(48)+JVS(280)*UV(49)+JVS(293)*UV(50)+JVS(306)*UV(51)+JVS(336)*UV(53)&
               &+JVS(353)*UV(54)+JVS(363)*UV(55)+JVS(387)*UV(57)+JVS(403)*UV(58)+JVS(446)*UV(60)+JVS(471)*UV(61)+JVS(504)&
               &*UV(62)+JVS(549)*UV(63)+JVS(565)*UV(64)+JVS(604)*UV(65)+JVS(624)*UV(66)
  JTUV(65) = JVS(16)*UV(5)+JVS(63)*UV(17)+JVS(95)*UV(26)+JVS(122)*UV(30)+JVS(127)*UV(31)+JVS(131)*UV(32)+JVS(162)*UV(37)&
               &+JVS(182)*UV(39)+JVS(200)*UV(40)+JVS(212)*UV(42)+JVS(217)*UV(43)+JVS(263)*UV(47)+JVS(307)*UV(51)+JVS(316)&
               &*UV(52)+JVS(354)*UV(54)+JVS(376)*UV(56)+JVS(388)*UV(57)+JVS(404)*UV(58)+JVS(420)*UV(59)+JVS(447)*UV(60)&
               &+JVS(472)*UV(61)+JVS(505)*UV(62)+JVS(550)*UV(63)+JVS(566)*UV(64)+JVS(605)*UV(65)+JVS(625)*UV(66)
  JTUV(66) = JVS(17)*UV(5)+JVS(66)*UV(18)+JVS(355)*UV(54)+JVS(448)*UV(60)+JVS(473)*UV(61)+JVS(506)*UV(62)+JVS(567)&
               &*UV(64)+JVS(606)*UV(65)+JVS(626)*UV(66)
      
END SUBROUTINE JacTR_SP_Vec

! End of JacTR_SP_Vec function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



END MODULE cbmz_mosaic_Jacobian

