! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Linear Algebra Data and Routines File
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
! File                 : nmhc9_LinearAlgebra.f90
! Time                 : Tue Jan 11 14:32:58 2022
! Working directory    : /jathar-scratch/yicongh/WRF_GoAmazon_gasaerochem.DEV/chem/KPP/mechanisms/nmhc9
! Equation file        : nmhc9.kpp
! Output root filename : nmhc9
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE nmhc9_LinearAlgebra

  USE nmhc9_Parameters
  USE nmhc9_JacobianSP

  IMPLICIT NONE

CONTAINS


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! KppSolveTR - sparse, transposed back substitution
!   Arguments :
!      JVS       - sparse Jacobian of variables
!      X         - Vector for variables
!      XX        - Vector for output variables
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE KppSolveTR ( JVS, X, XX )

! JVS - sparse Jacobian of variables
  REAL(kind=dp) :: JVS(LU_NONZERO)
! X - Vector for variables
  REAL(kind=dp) :: X(NVAR)
! XX - Vector for output variables
  REAL(kind=dp) :: XX(NVAR)

  XX(1) = X(1)/JVS(1)
  XX(2) = X(2)/JVS(3)
  XX(3) = X(3)/JVS(5)
  XX(4) = X(4)/JVS(8)
  XX(5) = X(5)/JVS(10)
  XX(6) = X(6)/JVS(13)
  XX(7) = X(7)/JVS(16)
  XX(8) = X(8)/JVS(19)
  XX(9) = X(9)/JVS(24)
  XX(10) = X(10)/JVS(28)
  XX(11) = X(11)/JVS(34)
  XX(12) = X(12)/JVS(38)
  XX(13) = X(13)/JVS(42)
  XX(14) = X(14)/JVS(46)
  XX(15) = X(15)/JVS(50)
  XX(16) = X(16)/JVS(57)
  XX(17) = X(17)/JVS(61)
  XX(18) = X(18)/JVS(65)
  XX(19) = (X(19)-JVS(51)*XX(15))/(JVS(69))
  XX(20) = (X(20)-JVS(6)*XX(3))/(JVS(72))
  XX(21) = X(21)/JVS(76)
  XX(22) = X(22)/JVS(79)
  XX(23) = X(23)/JVS(83)
  XX(24) = (X(24)-JVS(77)*XX(21))/(JVS(88))
  XX(25) = X(25)/JVS(93)
  XX(26) = X(26)/JVS(99)
  XX(27) = X(27)/JVS(103)
  XX(28) = X(28)/JVS(109)
  XX(29) = X(29)/JVS(115)
  XX(30) = X(30)/JVS(120)
  XX(31) = X(31)/JVS(130)
  XX(32) = X(32)/JVS(145)
  XX(33) = (X(33)-JVS(17)*XX(7)-JVS(131)*XX(31))/(JVS(156))
  XX(34) = (X(34)-JVS(62)*XX(17)-JVS(66)*XX(18)-JVS(110)*XX(28)-JVS(146)*XX(32))/(JVS(164))
  XX(35) = (X(35)-JVS(132)*XX(31))/(JVS(171))
  XX(36) = (X(36)-JVS(84)*XX(23)-JVS(147)*XX(32)-JVS(172)*XX(35))/(JVS(186))
  XX(37) = (X(37)-JVS(52)*XX(15)-JVS(89)*XX(24)-JVS(133)*XX(31)-JVS(148)*XX(32)-JVS(173)*XX(35))/(JVS(191))
  XX(38) = (X(38)-JVS(20)*XX(8)-JVS(53)*XX(15)-JVS(134)*XX(31)-JVS(157)*XX(33))/(JVS(195))
  XX(39) = (X(39)-JVS(174)*XX(35))/(JVS(205))
  XX(40) = (X(40)-JVS(100)*XX(26)-JVS(206)*XX(39))/(JVS(219))
  XX(41) = (X(41)-JVS(54)*XX(15)-JVS(135)*XX(31)-JVS(175)*XX(35))/(JVS(228))
  XX(42) = (X(42)-JVS(73)*XX(20))/(JVS(238))
  XX(43) = (X(43)-JVS(104)*XX(27)-JVS(121)*XX(30)-JVS(149)*XX(32)-JVS(239)*XX(42))/(JVS(247))
  XX(44) = (X(44)-JVS(122)*XX(30)-JVS(240)*XX(42))/(JVS(256))
  XX(45) = (X(45)-JVS(94)*XX(25)-JVS(136)*XX(31))/(JVS(283))
  XX(46) = (X(46)-JVS(95)*XX(25)-JVS(137)*XX(31))/(JVS(304))
  XX(47) = (X(47)-JVS(29)*XX(10)-JVS(58)*XX(16)-JVS(150)*XX(32)-JVS(284)*XX(45)-JVS(305)*XX(46))/(JVS(321))
  XX(48) = (X(48)-JVS(35)*XX(11)-JVS(151)*XX(32)-JVS(158)*XX(33)-JVS(176)*XX(35)-JVS(207)*XX(39)-JVS(229)*XX(41)&
             &-JVS(285)*XX(45))/(JVS(330))
  XX(49) = (X(49)-JVS(47)*XX(14)-JVS(116)*XX(29)-JVS(138)*XX(31)-JVS(152)*XX(32)-JVS(177)*XX(35)-JVS(208)*XX(39)&
             &-JVS(286)*XX(45))/(JVS(342))
  XX(50) = (X(50)-JVS(14)*XX(6)-JVS(30)*XX(10)-JVS(39)*XX(12)-JVS(111)*XX(28)-JVS(123)*XX(30)-JVS(139)*XX(31)-JVS(153)&
             &*XX(32)-JVS(165)*XX(34)-JVS(178)*XX(35)-JVS(187)*XX(36)-JVS(209)*XX(39)-JVS(230)*XX(41)-JVS(241)*XX(42)&
             &-JVS(248)*XX(43)-JVS(257)*XX(44)-JVS(287)*XX(45)-JVS(306)*XX(46)-JVS(322)*XX(47)-JVS(331)*XX(48)-JVS(343)&
             &*XX(49))/(JVS(366))
  XX(51) = (X(51)-JVS(25)*XX(9)-JVS(31)*XX(10)-JVS(43)*XX(13)-JVS(288)*XX(45)-JVS(307)*XX(46)-JVS(323)*XX(47)-JVS(367)&
             &*XX(50))/(JVS(398))
  XX(52) = (X(52)-JVS(11)*XX(5)-JVS(96)*XX(25)-JVS(140)*XX(31)-JVS(159)*XX(33)-JVS(192)*XX(37)-JVS(196)*XX(38)-JVS(210)&
             &*XX(39)-JVS(220)*XX(40)-JVS(231)*XX(41)-JVS(258)*XX(44)-JVS(289)*XX(45)-JVS(308)*XX(46)-JVS(324)*XX(47)&
             &-JVS(332)*XX(48)-JVS(344)*XX(49)-JVS(368)*XX(50)-JVS(399)*XX(51))/(JVS(418))
  XX(53) = (X(53)-JVS(21)*XX(8)-JVS(55)*XX(15)-JVS(70)*XX(19)-JVS(78)*XX(21)-JVS(90)*XX(24)-JVS(141)*XX(31)-JVS(154)&
             &*XX(32)-JVS(179)*XX(35)-JVS(193)*XX(37)-JVS(197)*XX(38)-JVS(211)*XX(39)-JVS(221)*XX(40)-JVS(232)*XX(41)&
             &-JVS(259)*XX(44)-JVS(290)*XX(45)-JVS(309)*XX(46)-JVS(325)*XX(47)-JVS(333)*XX(48)-JVS(345)*XX(49)-JVS(369)&
             &*XX(50)-JVS(400)*XX(51)-JVS(419)*XX(52))/(JVS(434))
  XX(54) = (X(54)-JVS(12)*XX(5)-JVS(15)*XX(6)-JVS(44)*XX(13)-JVS(48)*XX(14)-JVS(80)*XX(22)-JVS(97)*XX(25)-JVS(212)&
             &*XX(39)-JVS(291)*XX(45)-JVS(346)*XX(49)-JVS(370)*XX(50)-JVS(401)*XX(51)-JVS(420)*XX(52)-JVS(435)*XX(53))&
             &/(JVS(464))
  XX(55) = (X(55)-JVS(2)*XX(1)-JVS(4)*XX(2)-JVS(7)*XX(3)-JVS(9)*XX(4)-JVS(18)*XX(7)-JVS(22)*XX(8)-JVS(26)*XX(9)-JVS(32)&
             &*XX(10)-JVS(36)*XX(11)-JVS(40)*XX(12)-JVS(45)*XX(13)-JVS(49)*XX(14)-JVS(56)*XX(15)-JVS(59)*XX(16)-JVS(63)&
             &*XX(17)-JVS(67)*XX(18)-JVS(71)*XX(19)-JVS(74)*XX(20)-JVS(81)*XX(22)-JVS(85)*XX(23)-JVS(91)*XX(24)-JVS(98)&
             &*XX(25)-JVS(101)*XX(26)-JVS(105)*XX(27)-JVS(112)*XX(28)-JVS(117)*XX(29)-JVS(124)*XX(30)-JVS(142)*XX(31)&
             &-JVS(155)*XX(32)-JVS(160)*XX(33)-JVS(166)*XX(34)-JVS(180)*XX(35)-JVS(188)*XX(36)-JVS(194)*XX(37)-JVS(198)&
             &*XX(38)-JVS(213)*XX(39)-JVS(222)*XX(40)-JVS(233)*XX(41)-JVS(242)*XX(42)-JVS(249)*XX(43)-JVS(260)*XX(44)&
             &-JVS(292)*XX(45)-JVS(310)*XX(46)-JVS(326)*XX(47)-JVS(334)*XX(48)-JVS(347)*XX(49)-JVS(371)*XX(50)-JVS(402)&
             &*XX(51)-JVS(421)*XX(52)-JVS(436)*XX(53)-JVS(465)*XX(54))/(JVS(519))
  XX(56) = (X(56)-JVS(68)*XX(18)-JVS(113)*XX(28)-JVS(125)*XX(30)-JVS(143)*XX(31)-JVS(161)*XX(33)-JVS(167)*XX(34)&
             &-JVS(181)*XX(35)-JVS(189)*XX(36)-JVS(214)*XX(39)-JVS(223)*XX(40)-JVS(234)*XX(41)-JVS(243)*XX(42)-JVS(250)&
             &*XX(43)-JVS(261)*XX(44)-JVS(293)*XX(45)-JVS(311)*XX(46)-JVS(327)*XX(47)-JVS(335)*XX(48)-JVS(348)*XX(49)&
             &-JVS(372)*XX(50)-JVS(403)*XX(51)-JVS(422)*XX(52)-JVS(437)*XX(53)-JVS(466)*XX(54)-JVS(520)*XX(55))/(JVS(537))
  XX(57) = (X(57)-JVS(23)*XX(8)-JVS(27)*XX(9)-JVS(33)*XX(10)-JVS(37)*XX(11)-JVS(41)*XX(12)-JVS(60)*XX(16)-JVS(64)*XX(17)&
             &-JVS(75)*XX(20)-JVS(82)*XX(22)-JVS(86)*XX(23)-JVS(102)*XX(26)-JVS(106)*XX(27)-JVS(114)*XX(28)-JVS(118)*XX(29)&
             &-JVS(126)*XX(30)-JVS(144)*XX(31)-JVS(168)*XX(34)-JVS(182)*XX(35)-JVS(190)*XX(36)-JVS(215)*XX(39)-JVS(224)&
             &*XX(40)-JVS(235)*XX(41)-JVS(244)*XX(42)-JVS(251)*XX(43)-JVS(262)*XX(44)-JVS(294)*XX(45)-JVS(312)*XX(46)&
             &-JVS(328)*XX(47)-JVS(336)*XX(48)-JVS(349)*XX(49)-JVS(373)*XX(50)-JVS(404)*XX(51)-JVS(423)*XX(52)-JVS(438)&
             &*XX(53)-JVS(467)*XX(54)-JVS(521)*XX(55)-JVS(538)*XX(56))/(JVS(580))
  XX(57) = XX(57)
  XX(56) = XX(56)-JVS(579)*XX(57)
  XX(55) = XX(55)-JVS(536)*XX(56)-JVS(578)*XX(57)
  XX(54) = XX(54)-JVS(518)*XX(55)-JVS(535)*XX(56)-JVS(577)*XX(57)
  XX(53) = XX(53)-JVS(463)*XX(54)-JVS(517)*XX(55)-JVS(534)*XX(56)-JVS(576)*XX(57)
  XX(52) = XX(52)-JVS(433)*XX(53)-JVS(462)*XX(54)-JVS(516)*XX(55)-JVS(533)*XX(56)-JVS(575)*XX(57)
  XX(51) = XX(51)-JVS(417)*XX(52)-JVS(432)*XX(53)-JVS(461)*XX(54)-JVS(515)*XX(55)-JVS(532)*XX(56)-JVS(574)*XX(57)
  XX(50) = XX(50)-JVS(397)*XX(51)-JVS(416)*XX(52)-JVS(431)*XX(53)-JVS(460)*XX(54)-JVS(514)*XX(55)-JVS(531)*XX(56)&
             &-JVS(573)*XX(57)
  XX(49) = XX(49)-JVS(365)*XX(50)-JVS(396)*XX(51)-JVS(415)*XX(52)-JVS(459)*XX(54)-JVS(513)*XX(55)-JVS(530)*XX(56)&
             &-JVS(572)*XX(57)
  XX(48) = XX(48)-JVS(341)*XX(49)-JVS(364)*XX(50)-JVS(395)*XX(51)-JVS(414)*XX(52)-JVS(430)*XX(53)-JVS(458)*XX(54)&
             &-JVS(512)*XX(55)-JVS(529)*XX(56)-JVS(571)*XX(57)
  XX(47) = XX(47)-JVS(363)*XX(50)-JVS(394)*XX(51)-JVS(413)*XX(52)-JVS(457)*XX(54)-JVS(511)*XX(55)-JVS(528)*XX(56)&
             &-JVS(570)*XX(57)
  XX(46) = XX(46)-JVS(362)*XX(50)-JVS(393)*XX(51)-JVS(412)*XX(52)-JVS(456)*XX(54)-JVS(510)*XX(55)-JVS(569)*XX(57)
  XX(45) = XX(45)-JVS(411)*XX(52)-JVS(455)*XX(54)-JVS(509)*XX(55)-JVS(568)*XX(57)
  XX(44) = XX(44)-JVS(303)*XX(46)-JVS(320)*XX(47)-JVS(392)*XX(51)-JVS(454)*XX(54)-JVS(508)*XX(55)-JVS(527)*XX(56)&
             &-JVS(567)*XX(57)
  XX(43) = XX(43)-JVS(255)*XX(44)-JVS(282)*XX(45)-JVS(302)*XX(46)-JVS(319)*XX(47)-JVS(361)*XX(50)-JVS(391)*XX(51)&
             &-JVS(453)*XX(54)-JVS(507)*XX(55)-JVS(526)*XX(56)-JVS(566)*XX(57)
  XX(42) = XX(42)-JVS(254)*XX(44)-JVS(301)*XX(46)-JVS(390)*XX(51)-JVS(452)*XX(54)-JVS(506)*XX(55)-JVS(525)*XX(56)&
             &-JVS(565)*XX(57)
  XX(41) = XX(41)-JVS(281)*XX(45)-JVS(340)*XX(49)-JVS(389)*XX(51)-JVS(429)*XX(53)-JVS(505)*XX(55)-JVS(564)*XX(57)
  XX(40) = XX(40)-JVS(253)*XX(44)-JVS(280)*XX(45)-JVS(300)*XX(46)-JVS(388)*XX(51)-JVS(451)*XX(54)-JVS(504)*XX(55)&
             &-JVS(524)*XX(56)-JVS(563)*XX(57)
  XX(39) = XX(39)-JVS(279)*XX(45)-JVS(387)*XX(51)-JVS(503)*XX(55)-JVS(562)*XX(57)
  XX(38) = XX(38)-JVS(204)*XX(39)-JVS(227)*XX(41)-JVS(278)*XX(45)-JVS(329)*XX(48)-JVS(339)*XX(49)-JVS(360)*XX(50)&
             &-JVS(386)*XX(51)-JVS(410)*XX(52)-JVS(428)*XX(53)-JVS(450)*XX(54)-JVS(502)*XX(55)-JVS(561)*XX(57)
  XX(37) = XX(37)-JVS(218)*XX(40)-JVS(252)*XX(44)-JVS(277)*XX(45)-JVS(299)*XX(46)-JVS(318)*XX(47)-JVS(359)*XX(50)&
             &-JVS(385)*XX(51)-JVS(409)*XX(52)-JVS(427)*XX(53)-JVS(501)*XX(55)-JVS(560)*XX(57)
  XX(36) = XX(36)-JVS(203)*XX(39)-JVS(276)*XX(45)-JVS(358)*XX(50)-JVS(384)*XX(51)-JVS(449)*XX(54)-JVS(500)*XX(55)&
             &-JVS(523)*XX(56)-JVS(559)*XX(57)
  XX(35) = XX(35)-JVS(383)*XX(51)-JVS(499)*XX(55)-JVS(558)*XX(57)
  XX(34) = XX(34)-JVS(185)*XX(36)-JVS(275)*XX(45)-JVS(357)*XX(50)-JVS(382)*XX(51)-JVS(448)*XX(54)-JVS(498)*XX(55)&
             &-JVS(522)*XX(56)-JVS(557)*XX(57)
  XX(33) = XX(33)-JVS(202)*XX(39)-JVS(226)*XX(41)-JVS(274)*XX(45)-JVS(447)*XX(54)-JVS(497)*XX(55)-JVS(556)*XX(57)
  XX(32) = XX(32)-JVS(273)*XX(45)-JVS(496)*XX(55)-JVS(555)*XX(57)
  XX(31) = XX(31)-JVS(495)*XX(55)-JVS(554)*XX(57)
  XX(30) = XX(30)-JVS(237)*XX(42)-JVS(317)*XX(47)-JVS(381)*XX(51)-JVS(494)*XX(55)
  XX(29) = XX(29)-JVS(129)*XX(31)-JVS(170)*XX(35)-JVS(201)*XX(39)-JVS(272)*XX(45)-JVS(338)*XX(49)-JVS(380)*XX(51)&
             &-JVS(493)*XX(55)-JVS(553)*XX(57)
  XX(28) = XX(28)-JVS(184)*XX(36)-JVS(356)*XX(50)-JVS(379)*XX(51)-JVS(492)*XX(55)
  XX(27) = XX(27)-JVS(119)*XX(30)-JVS(246)*XX(43)-JVS(298)*XX(46)-JVS(316)*XX(47)-JVS(491)*XX(55)-JVS(552)*XX(57)
  XX(26) = XX(26)-JVS(200)*XX(39)-JVS(217)*XX(40)-JVS(271)*XX(45)-JVS(297)*XX(46)-JVS(490)*XX(55)-JVS(551)*XX(57)
  XX(25) = XX(25)-JVS(408)*XX(52)-JVS(446)*XX(54)-JVS(489)*XX(55)
  XX(24) = XX(24)-JVS(270)*XX(45)-JVS(355)*XX(50)-JVS(426)*XX(53)-JVS(488)*XX(55)-JVS(550)*XX(57)
  XX(23) = XX(23)-JVS(169)*XX(35)-JVS(183)*XX(36)-JVS(378)*XX(51)-JVS(487)*XX(55)-JVS(549)*XX(57)
  XX(22) = XX(22)-JVS(407)*XX(52)-JVS(445)*XX(54)-JVS(486)*XX(55)-JVS(548)*XX(57)
  XX(21) = XX(21)-JVS(87)*XX(24)-JVS(269)*XX(45)-JVS(354)*XX(50)-JVS(425)*XX(53)-JVS(485)*XX(55)-JVS(547)*XX(57)
  XX(20) = XX(20)-JVS(236)*XX(42)-JVS(296)*XX(46)-JVS(377)*XX(51)-JVS(484)*XX(55)
  XX(19) = XX(19)-JVS(128)*XX(31)-JVS(216)*XX(40)-JVS(268)*XX(45)-JVS(424)*XX(53)-JVS(483)*XX(55)-JVS(546)*XX(57)
  XX(18) = XX(18)-JVS(108)*XX(28)-JVS(444)*XX(54)-JVS(482)*XX(55)-JVS(545)*XX(57)
  XX(17) = XX(17)-JVS(107)*XX(28)-JVS(163)*XX(34)-JVS(481)*XX(55)-JVS(544)*XX(57)
  XX(16) = XX(16)-JVS(295)*XX(46)-JVS(315)*XX(47)-JVS(480)*XX(55)-JVS(543)*XX(57)
  XX(15) = XX(15)-JVS(479)*XX(55)-JVS(542)*XX(57)
  XX(14) = XX(14)-JVS(199)*XX(39)-JVS(337)*XX(49)-JVS(443)*XX(54)-JVS(478)*XX(55)
  XX(13) = XX(13)-JVS(267)*XX(45)-JVS(376)*XX(51)-JVS(442)*XX(54)-JVS(477)*XX(55)
  XX(12) = XX(12)-JVS(266)*XX(45)-JVS(353)*XX(50)-JVS(476)*XX(55)-JVS(541)*XX(57)
  XX(11) = XX(11)-JVS(225)*XX(41)-JVS(265)*XX(45)-JVS(475)*XX(55)-JVS(540)*XX(57)
  XX(10) = XX(10)-JVS(352)*XX(50)-JVS(474)*XX(55)
  XX(9) = XX(9)-JVS(351)*XX(50)-JVS(375)*XX(51)-JVS(473)*XX(55)
  XX(8) = XX(8)-JVS(472)*XX(55)-JVS(539)*XX(57)
  XX(7) = XX(7)-JVS(127)*XX(31)-JVS(264)*XX(45)-JVS(441)*XX(54)-JVS(471)*XX(55)
  XX(6) = XX(6)-JVS(263)*XX(45)-JVS(350)*XX(50)-JVS(406)*XX(52)-JVS(440)*XX(54)
  XX(5) = XX(5)-JVS(92)*XX(25)-JVS(405)*XX(52)-JVS(439)*XX(54)
  XX(4) = XX(4)-JVS(162)*XX(34)-JVS(314)*XX(47)-JVS(470)*XX(55)
  XX(3) = XX(3)-JVS(374)*XX(51)
  XX(2) = XX(2)-JVS(245)*XX(43)-JVS(469)*XX(55)
  XX(1) = XX(1)-JVS(313)*XX(47)-JVS(468)*XX(55)
      
END SUBROUTINE KppSolveTR

! End of KppSolveTR function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



END MODULE nmhc9_LinearAlgebra

