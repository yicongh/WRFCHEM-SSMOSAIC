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
! File                 : mozart_LinearAlgebra.f90
! Time                 : Tue Jan 11 14:32:55 2022
! Working directory    : /jathar-scratch/yicongh/WRF_GoAmazon_gasaerochem.DEV/chem/KPP/mechanisms/mozart
! Equation file        : mozart.kpp
! Output root filename : mozart
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE mozart_LinearAlgebra

  USE mozart_Parameters
  USE mozart_JacobianSP

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
  XX(2) = X(2)/JVS(4)
  XX(3) = X(3)/JVS(6)
  XX(4) = X(4)/JVS(8)
  XX(5) = X(5)/JVS(10)
  XX(6) = X(6)/JVS(12)
  XX(7) = X(7)/JVS(14)
  XX(8) = (X(8)-JVS(2)*XX(1))/(JVS(16))
  XX(9) = X(9)/JVS(20)
  XX(10) = X(10)/JVS(23)
  XX(11) = X(11)/JVS(26)
  XX(12) = X(12)/JVS(28)
  XX(13) = X(13)/JVS(31)
  XX(14) = X(14)/JVS(35)
  XX(15) = X(15)/JVS(38)
  XX(16) = (X(16)-JVS(17)*XX(8))/(JVS(42))
  XX(17) = X(17)/JVS(45)
  XX(18) = X(18)/JVS(49)
  XX(19) = X(19)/JVS(53)
  XX(20) = X(20)/JVS(58)
  XX(21) = X(21)/JVS(62)
  XX(22) = X(22)/JVS(66)
  XX(23) = X(23)/JVS(70)
  XX(24) = X(24)/JVS(75)
  XX(25) = X(25)/JVS(78)
  XX(26) = X(26)/JVS(82)
  XX(27) = X(27)/JVS(86)
  XX(28) = X(28)/JVS(89)
  XX(29) = X(29)/JVS(97)
  XX(30) = (X(30)-JVS(29)*XX(12))/(JVS(104))
  XX(31) = X(31)/JVS(108)
  XX(32) = X(32)/JVS(112)
  XX(33) = X(33)/JVS(117)
  XX(34) = X(34)/JVS(123)
  XX(35) = X(35)/JVS(132)
  XX(36) = X(36)/JVS(136)
  XX(37) = X(37)/JVS(141)
  XX(38) = X(38)/JVS(148)
  XX(39) = X(39)/JVS(154)
  XX(40) = (X(40)-JVS(54)*XX(19))/(JVS(161))
  XX(41) = X(41)/JVS(166)
  XX(42) = X(42)/JVS(174)
  XX(43) = X(43)/JVS(178)
  XX(44) = (X(44)-JVS(109)*XX(31)-JVS(142)*XX(37))/(JVS(184))
  XX(45) = (X(45)-JVS(143)*XX(37))/(JVS(191))
  XX(46) = X(46)/JVS(196)
  XX(47) = (X(47)-JVS(7)*XX(3)-JVS(55)*XX(19)-JVS(162)*XX(40))/(JVS(203))
  XX(48) = (X(48)-JVS(118)*XX(33))/(JVS(208))
  XX(49) = (X(49)-JVS(79)*XX(25))/(JVS(215))
  XX(50) = (X(50)-JVS(133)*XX(35))/(JVS(221))
  XX(51) = X(51)/JVS(230)
  XX(52) = X(52)/JVS(239)
  XX(53) = (X(53)-JVS(144)*XX(37))/(JVS(248))
  XX(54) = (X(54)-JVS(98)*XX(29)-JVS(119)*XX(33)-JVS(209)*XX(48)-JVS(216)*XX(49)-JVS(231)*XX(51))/(JVS(259))
  XX(55) = (X(55)-JVS(46)*XX(17)-JVS(63)*XX(21)-JVS(167)*XX(41))/(JVS(268))
  XX(56) = X(56)/JVS(273)
  XX(57) = (X(57)-JVS(67)*XX(22)-JVS(232)*XX(51))/(JVS(280))
  XX(58) = X(58)/JVS(290)
  XX(59) = (X(59)-JVS(124)*XX(34)-JVS(274)*XX(56)-JVS(291)*XX(58))/(JVS(309))
  XX(60) = (X(60)-JVS(90)*XX(28)-JVS(292)*XX(58))/(JVS(315))
  XX(61) = (X(61)-JVS(99)*XX(29)-JVS(125)*XX(34)-JVS(163)*XX(40)-JVS(204)*XX(47)-JVS(222)*XX(50)-JVS(240)*XX(52)&
             &-JVS(293)*XX(58))/(JVS(324))
  XX(62) = (X(62)-JVS(155)*XX(39)-JVS(294)*XX(58))/(JVS(341))
  XX(63) = (X(63)-JVS(71)*XX(23)-JVS(168)*XX(41))/(JVS(353))
  XX(64) = (X(64)-JVS(149)*XX(38)-JVS(233)*XX(51)-JVS(354)*XX(63))/(JVS(361))
  XX(65) = X(65)/JVS(373)
  XX(66) = (X(66)-JVS(295)*XX(58))/(JVS(392))
  XX(67) = (X(67)-JVS(156)*XX(39)-JVS(296)*XX(58))/(JVS(411))
  XX(68) = (X(68)-JVS(39)*XX(15)-JVS(169)*XX(41)-JVS(249)*XX(53)-JVS(297)*XX(58)-JVS(374)*XX(65)-JVS(412)*XX(67))&
             &/(JVS(429))
  XX(69) = (X(69)-JVS(298)*XX(58)-JVS(325)*XX(61)-JVS(342)*XX(62)-JVS(375)*XX(65)-JVS(413)*XX(67))/(JVS(445))
  XX(70) = (X(70)-JVS(59)*XX(20)-JVS(170)*XX(41)-JVS(250)*XX(53)-JVS(299)*XX(58)-JVS(316)*XX(60)-JVS(376)*XX(65)&
             &-JVS(414)*XX(67)-JVS(430)*XX(68))/(JVS(457))
  XX(71) = (X(71)-JVS(91)*XX(28)-JVS(171)*XX(41)-JVS(175)*XX(42)-JVS(317)*XX(60)-JVS(393)*XX(66)-JVS(415)*XX(67)&
             &-JVS(431)*XX(68)-JVS(446)*XX(69)-JVS(458)*XX(70))/(JVS(468))
  XX(72) = (X(72)-JVS(56)*XX(19)-JVS(157)*XX(39)-JVS(205)*XX(47)-JVS(241)*XX(52)-JVS(300)*XX(58))/(JVS(508))
  XX(73) = (X(73)-JVS(32)*XX(13)-JVS(36)*XX(14)-JVS(83)*XX(26)-JVS(137)*XX(36)-JVS(158)*XX(39)-JVS(179)*XX(43)-JVS(192)&
             &*XX(45)-JVS(242)*XX(52)-JVS(301)*XX(58)-JVS(377)*XX(65)-JVS(416)*XX(67)-JVS(509)*XX(72))/(JVS(548))
  XX(74) = (X(74)-JVS(113)*XX(32)-JVS(126)*XX(34)-JVS(180)*XX(43)-JVS(378)*XX(65)-JVS(417)*XX(67)-JVS(510)*XX(72)&
             &-JVS(549)*XX(73))/(JVS(564))
  XX(75) = (X(75)-JVS(3)*XX(1)-JVS(5)*XX(2)-JVS(9)*XX(4)-JVS(11)*XX(5)-JVS(13)*XX(6)-JVS(15)*XX(7)-JVS(18)*XX(8)-JVS(21)&
             &*XX(9)-JVS(24)*XX(10)-JVS(27)*XX(11)-JVS(37)*XX(14)-JVS(40)*XX(15)-JVS(43)*XX(16)-JVS(47)*XX(17)-JVS(50)&
             &*XX(18)-JVS(57)*XX(19)-JVS(60)*XX(20)-JVS(64)*XX(21)-JVS(68)*XX(22)-JVS(72)*XX(23)-JVS(76)*XX(24)-JVS(80)&
             &*XX(25)-JVS(84)*XX(26)-JVS(87)*XX(27)-JVS(92)*XX(28)-JVS(100)*XX(29)-JVS(105)*XX(30)-JVS(110)*XX(31)-JVS(114)&
             &*XX(32)-JVS(120)*XX(33)-JVS(127)*XX(34)-JVS(134)*XX(35)-JVS(138)*XX(36)-JVS(145)*XX(37)-JVS(150)*XX(38)&
             &-JVS(159)*XX(39)-JVS(164)*XX(40)-JVS(172)*XX(41)-JVS(176)*XX(42)-JVS(181)*XX(43)-JVS(185)*XX(44)-JVS(193)&
             &*XX(45)-JVS(197)*XX(46)-JVS(206)*XX(47)-JVS(210)*XX(48)-JVS(217)*XX(49)-JVS(223)*XX(50)-JVS(234)*XX(51)&
             &-JVS(243)*XX(52)-JVS(251)*XX(53)-JVS(260)*XX(54)-JVS(269)*XX(55)-JVS(281)*XX(57)-JVS(302)*XX(58)-JVS(310)&
             &*XX(59)-JVS(318)*XX(60)-JVS(326)*XX(61)-JVS(343)*XX(62)-JVS(355)*XX(63)-JVS(362)*XX(64)-JVS(379)*XX(65)&
             &-JVS(394)*XX(66)-JVS(418)*XX(67)-JVS(432)*XX(68)-JVS(447)*XX(69)-JVS(459)*XX(70)-JVS(469)*XX(71)-JVS(511)&
             &*XX(72)-JVS(550)*XX(73)-JVS(565)*XX(74))/(JVS(639))
  XX(76) = (X(76)-JVS(30)*XX(12)-JVS(77)*XX(24)-JVS(93)*XX(28)-JVS(101)*XX(29)-JVS(106)*XX(30)-JVS(121)*XX(33)-JVS(146)&
             &*XX(37)-JVS(186)*XX(44)-JVS(194)*XX(45)-JVS(218)*XX(49)-JVS(224)*XX(50)-JVS(235)*XX(51)-JVS(252)*XX(53)&
             &-JVS(261)*XX(54)-JVS(270)*XX(55)-JVS(275)*XX(56)-JVS(282)*XX(57)-JVS(303)*XX(58)-JVS(319)*XX(60)-JVS(344)&
             &*XX(62)-JVS(356)*XX(63)-JVS(363)*XX(64)-JVS(380)*XX(65)-JVS(395)*XX(66)-JVS(419)*XX(67)-JVS(433)*XX(68)&
             &-JVS(448)*XX(69)-JVS(460)*XX(70)-JVS(470)*XX(71)-JVS(512)*XX(72)-JVS(551)*XX(73)-JVS(566)*XX(74)-JVS(640)&
             &*XX(75))/(JVS(670))
  XX(77) = (X(77)-JVS(48)*XX(17)-JVS(51)*XX(18)-JVS(94)*XX(28)-JVS(128)*XX(34)-JVS(173)*XX(41)-JVS(236)*XX(51)-JVS(253)&
             &*XX(53)-JVS(271)*XX(55)-JVS(283)*XX(57)-JVS(304)*XX(58)-JVS(345)*XX(62)-JVS(357)*XX(63)-JVS(381)*XX(65)&
             &-JVS(396)*XX(66)-JVS(420)*XX(67)-JVS(434)*XX(68)-JVS(449)*XX(69)-JVS(461)*XX(70)-JVS(471)*XX(71)-JVS(513)&
             &*XX(72)-JVS(552)*XX(73)-JVS(567)*XX(74)-JVS(641)*XX(75)-JVS(671)*XX(76))/(JVS(700))
  XX(78) = (X(78)-JVS(95)*XX(28)-JVS(115)*XX(32)-JVS(129)*XX(34)-JVS(139)*XX(36)-JVS(254)*XX(53)-JVS(305)*XX(58)&
             &-JVS(382)*XX(65)-JVS(397)*XX(66)-JVS(421)*XX(67)-JVS(435)*XX(68)-JVS(450)*XX(69)-JVS(462)*XX(70)-JVS(472)&
             &*XX(71)-JVS(514)*XX(72)-JVS(553)*XX(73)-JVS(568)*XX(74)-JVS(642)*XX(75)-JVS(672)*XX(76)-JVS(701)*XX(77))&
             &/(JVS(731))
  XX(79) = (X(79)-JVS(22)*XX(9)-JVS(41)*XX(15)-JVS(52)*XX(18)-JVS(61)*XX(20)-JVS(65)*XX(21)-JVS(69)*XX(22)-JVS(73)&
             &*XX(23)-JVS(81)*XX(25)-JVS(85)*XX(26)-JVS(111)*XX(31)-JVS(116)*XX(32)-JVS(130)*XX(34)-JVS(135)*XX(35)-JVS(147)&
             &*XX(37)-JVS(151)*XX(38)-JVS(177)*XX(42)-JVS(187)*XX(44)-JVS(195)*XX(45)-JVS(211)*XX(48)-JVS(219)*XX(49)&
             &-JVS(225)*XX(50)-JVS(237)*XX(51)-JVS(244)*XX(52)-JVS(262)*XX(54)-JVS(272)*XX(55)-JVS(276)*XX(56)-JVS(284)&
             &*XX(57)-JVS(306)*XX(58)-JVS(320)*XX(60)-JVS(346)*XX(62)-JVS(358)*XX(63)-JVS(364)*XX(64)-JVS(383)*XX(65)&
             &-JVS(398)*XX(66)-JVS(422)*XX(67)-JVS(436)*XX(68)-JVS(451)*XX(69)-JVS(463)*XX(70)-JVS(473)*XX(71)-JVS(515)&
             &*XX(72)-JVS(554)*XX(73)-JVS(569)*XX(74)-JVS(643)*XX(75)-JVS(673)*XX(76)-JVS(702)*XX(77)-JVS(732)*XX(78))&
             &/(JVS(791))
  XX(80) = (X(80)-JVS(19)*XX(8)-JVS(33)*XX(13)-JVS(44)*XX(16)-JVS(96)*XX(28)-JVS(102)*XX(29)-JVS(160)*XX(39)-JVS(198)&
             &*XX(46)-JVS(255)*XX(53)-JVS(277)*XX(56)-JVS(307)*XX(58)-JVS(311)*XX(59)-JVS(321)*XX(60)-JVS(327)*XX(61)&
             &-JVS(347)*XX(62)-JVS(365)*XX(64)-JVS(384)*XX(65)-JVS(399)*XX(66)-JVS(423)*XX(67)-JVS(437)*XX(68)-JVS(452)&
             &*XX(69)-JVS(464)*XX(70)-JVS(474)*XX(71)-JVS(516)*XX(72)-JVS(555)*XX(73)-JVS(570)*XX(74)-JVS(644)*XX(75)&
             &-JVS(674)*XX(76)-JVS(703)*XX(77)-JVS(733)*XX(78)-JVS(792)*XX(79))/(JVS(819))
  XX(81) = (X(81)-JVS(88)*XX(27)-JVS(107)*XX(30)-JVS(131)*XX(34)-JVS(165)*XX(40)-JVS(199)*XX(46)-JVS(207)*XX(47)&
             &-JVS(245)*XX(52)-JVS(256)*XX(53)-JVS(308)*XX(58)-JVS(312)*XX(59)-JVS(322)*XX(60)-JVS(328)*XX(61)-JVS(348)&
             &*XX(62)-JVS(366)*XX(64)-JVS(385)*XX(65)-JVS(400)*XX(66)-JVS(424)*XX(67)-JVS(438)*XX(68)-JVS(453)*XX(69)&
             &-JVS(465)*XX(70)-JVS(475)*XX(71)-JVS(517)*XX(72)-JVS(556)*XX(73)-JVS(571)*XX(74)-JVS(645)*XX(75)-JVS(675)&
             &*XX(76)-JVS(704)*XX(77)-JVS(734)*XX(78)-JVS(793)*XX(79)-JVS(820)*XX(80))/(JVS(838))
  XX(81) = XX(81)
  XX(80) = XX(80)-JVS(837)*XX(81)
  XX(79) = XX(79)-JVS(818)*XX(80)-JVS(836)*XX(81)
  XX(78) = XX(78)-JVS(790)*XX(79)-JVS(817)*XX(80)-JVS(835)*XX(81)
  XX(77) = XX(77)-JVS(730)*XX(78)-JVS(789)*XX(79)-JVS(816)*XX(80)-JVS(834)*XX(81)
  XX(76) = XX(76)-JVS(699)*XX(77)-JVS(729)*XX(78)-JVS(788)*XX(79)-JVS(815)*XX(80)-JVS(833)*XX(81)
  XX(75) = XX(75)-JVS(669)*XX(76)-JVS(698)*XX(77)-JVS(728)*XX(78)-JVS(787)*XX(79)-JVS(814)*XX(80)-JVS(832)*XX(81)
  XX(74) = XX(74)-JVS(638)*XX(75)-JVS(668)*XX(76)-JVS(697)*XX(77)-JVS(727)*XX(78)-JVS(786)*XX(79)-JVS(813)*XX(80)&
             &-JVS(831)*XX(81)
  XX(73) = XX(73)-JVS(563)*XX(74)-JVS(637)*XX(75)-JVS(667)*XX(76)-JVS(696)*XX(77)-JVS(726)*XX(78)-JVS(785)*XX(79)&
             &-JVS(812)*XX(80)-JVS(830)*XX(81)
  XX(72) = XX(72)-JVS(547)*XX(73)-JVS(636)*XX(75)-JVS(666)*XX(76)-JVS(695)*XX(77)-JVS(784)*XX(79)-JVS(811)*XX(80)&
             &-JVS(829)*XX(81)
  XX(71) = XX(71)-JVS(507)*XX(72)-JVS(546)*XX(73)-JVS(562)*XX(74)-JVS(635)*XX(75)-JVS(665)*XX(76)-JVS(694)*XX(77)&
             &-JVS(725)*XX(78)-JVS(783)*XX(79)-JVS(810)*XX(80)-JVS(828)*XX(81)
  XX(70) = XX(70)-JVS(506)*XX(72)-JVS(545)*XX(73)-JVS(561)*XX(74)-JVS(634)*XX(75)-JVS(664)*XX(76)-JVS(693)*XX(77)&
             &-JVS(724)*XX(78)-JVS(782)*XX(79)-JVS(809)*XX(80)
  XX(69) = XX(69)-JVS(456)*XX(70)-JVS(505)*XX(72)-JVS(544)*XX(73)-JVS(633)*XX(75)-JVS(663)*XX(76)-JVS(692)*XX(77)&
             &-JVS(723)*XX(78)-JVS(781)*XX(79)-JVS(808)*XX(80)-JVS(827)*XX(81)
  XX(68) = XX(68)-JVS(504)*XX(72)-JVS(543)*XX(73)-JVS(632)*XX(75)-JVS(662)*XX(76)-JVS(691)*XX(77)-JVS(722)*XX(78)&
             &-JVS(780)*XX(79)-JVS(807)*XX(80)
  XX(67) = XX(67)-JVS(542)*XX(73)-JVS(631)*XX(75)-JVS(721)*XX(78)-JVS(779)*XX(79)-JVS(806)*XX(80)
  XX(66) = XX(66)-JVS(410)*XX(67)-JVS(455)*XX(70)-JVS(503)*XX(72)-JVS(560)*XX(74)-JVS(630)*XX(75)-JVS(720)*XX(78)&
             &-JVS(778)*XX(79)-JVS(826)*XX(81)
  XX(65) = XX(65)-JVS(409)*XX(67)-JVS(502)*XX(72)-JVS(629)*XX(75)-JVS(719)*XX(78)-JVS(777)*XX(79)
  XX(64) = XX(64)-JVS(372)*XX(65)-JVS(391)*XX(66)-JVS(408)*XX(67)-JVS(444)*XX(69)-JVS(501)*XX(72)-JVS(541)*XX(73)&
             &-JVS(628)*XX(75)-JVS(661)*XX(76)-JVS(690)*XX(77)-JVS(718)*XX(78)-JVS(776)*XX(79)
  XX(63) = XX(63)-JVS(371)*XX(65)-JVS(407)*XX(67)-JVS(500)*XX(72)-JVS(540)*XX(73)-JVS(627)*XX(75)-JVS(660)*XX(76)&
             &-JVS(689)*XX(77)-JVS(717)*XX(78)-JVS(775)*XX(79)
  XX(62) = XX(62)-JVS(539)*XX(73)-JVS(626)*XX(75)-JVS(688)*XX(77)-JVS(716)*XX(78)-JVS(774)*XX(79)-JVS(805)*XX(80)
  XX(61) = XX(61)-JVS(340)*XX(62)-JVS(370)*XX(65)-JVS(406)*XX(67)-JVS(499)*XX(72)-JVS(538)*XX(73)-JVS(625)*XX(75)&
             &-JVS(659)*XX(76)-JVS(687)*XX(77)-JVS(773)*XX(79)-JVS(804)*XX(80)-JVS(825)*XX(81)
  XX(60) = XX(60)-JVS(428)*XX(68)-JVS(498)*XX(72)-JVS(537)*XX(73)-JVS(624)*XX(75)-JVS(772)*XX(79)-JVS(803)*XX(80)
  XX(59) = XX(59)-JVS(314)*XX(60)-JVS(323)*XX(61)-JVS(390)*XX(66)-JVS(443)*XX(69)-JVS(467)*XX(71)-JVS(497)*XX(72)&
             &-JVS(536)*XX(73)-JVS(559)*XX(74)-JVS(623)*XX(75)-JVS(658)*XX(76)-JVS(686)*XX(77)-JVS(771)*XX(79)-JVS(802)&
             &*XX(80)-JVS(824)*XX(81)
  XX(58) = XX(58)-JVS(622)*XX(75)-JVS(770)*XX(79)
  XX(57) = XX(57)-JVS(339)*XX(62)-JVS(352)*XX(63)-JVS(496)*XX(72)-JVS(535)*XX(73)-JVS(621)*XX(75)-JVS(657)*XX(76)&
             &-JVS(685)*XX(77)-JVS(715)*XX(78)-JVS(769)*XX(79)
  XX(56) = XX(56)-JVS(313)*XX(60)-JVS(389)*XX(66)-JVS(442)*XX(69)-JVS(495)*XX(72)-JVS(534)*XX(73)-JVS(620)*XX(75)&
             &-JVS(656)*XX(76)-JVS(768)*XX(79)-JVS(801)*XX(80)
  XX(55) = XX(55)-JVS(338)*XX(62)-JVS(494)*XX(72)-JVS(533)*XX(73)-JVS(619)*XX(75)-JVS(655)*XX(76)-JVS(684)*XX(77)&
             &-JVS(767)*XX(79)
  XX(54) = XX(54)-JVS(267)*XX(55)-JVS(337)*XX(62)-JVS(351)*XX(63)-JVS(405)*XX(67)-JVS(493)*XX(72)-JVS(532)*XX(73)&
             &-JVS(618)*XX(75)-JVS(654)*XX(76)-JVS(683)*XX(77)-JVS(714)*XX(78)-JVS(766)*XX(79)
  XX(53) = XX(53)-JVS(289)*XX(58)-JVS(492)*XX(72)-JVS(617)*XX(75)-JVS(765)*XX(79)
  XX(52) = XX(52)-JVS(531)*XX(73)-JVS(616)*XX(75)-JVS(653)*XX(76)-JVS(764)*XX(79)-JVS(823)*XX(81)
  XX(51) = XX(51)-JVS(350)*XX(63)-JVS(615)*XX(75)-JVS(682)*XX(77)-JVS(713)*XX(78)
  XX(50) = XX(50)-JVS(336)*XX(62)-JVS(369)*XX(65)-JVS(491)*XX(72)-JVS(530)*XX(73)-JVS(614)*XX(75)-JVS(652)*XX(76)&
             &-JVS(763)*XX(79)
  XX(49) = XX(49)-JVS(335)*XX(62)-JVS(529)*XX(73)-JVS(613)*XX(75)-JVS(651)*XX(76)-JVS(712)*XX(78)-JVS(762)*XX(79)
  XX(48) = XX(48)-JVS(214)*XX(49)-JVS(229)*XX(51)-JVS(258)*XX(54)-JVS(266)*XX(55)-JVS(334)*XX(62)-JVS(490)*XX(72)&
             &-JVS(612)*XX(75)-JVS(711)*XX(78)-JVS(761)*XX(79)
  XX(47) = XX(47)-JVS(238)*XX(52)-JVS(489)*XX(72)-JVS(611)*XX(75)-JVS(650)*XX(76)-JVS(681)*XX(77)-JVS(760)*XX(79)
  XX(46) = XX(46)-JVS(360)*XX(64)-JVS(388)*XX(66)-JVS(441)*XX(69)-JVS(528)*XX(73)-JVS(610)*XX(75)-JVS(759)*XX(79)&
             &-JVS(800)*XX(80)-JVS(822)*XX(81)
  XX(45) = XX(45)-JVS(288)*XX(58)-JVS(404)*XX(67)-JVS(609)*XX(75)-JVS(710)*XX(78)-JVS(758)*XX(79)
  XX(44) = XX(44)-JVS(190)*XX(45)-JVS(287)*XX(58)-JVS(403)*XX(67)-JVS(527)*XX(73)-JVS(608)*XX(75)-JVS(649)*XX(76)&
             &-JVS(757)*XX(79)
  XX(43) = XX(43)-JVS(368)*XX(65)-JVS(488)*XX(72)-JVS(526)*XX(73)-JVS(558)*XX(74)-JVS(607)*XX(75)-JVS(756)*XX(79)&
             &-JVS(799)*XX(80)
  XX(42) = XX(42)-JVS(387)*XX(66)-JVS(427)*XX(68)-JVS(440)*XX(69)-JVS(466)*XX(71)-JVS(487)*XX(72)-JVS(606)*XX(75)&
             &-JVS(755)*XX(79)
  XX(41) = XX(41)-JVS(486)*XX(72)-JVS(605)*XX(75)-JVS(754)*XX(79)
  XX(40) = XX(40)-JVS(202)*XX(47)-JVS(485)*XX(72)-JVS(604)*XX(75)-JVS(680)*XX(77)-JVS(753)*XX(79)
  XX(39) = XX(39)-JVS(525)*XX(73)-JVS(603)*XX(75)-JVS(798)*XX(80)
  XX(38) = XX(38)-JVS(228)*XX(51)-JVS(359)*XX(64)-JVS(386)*XX(66)-JVS(439)*XX(69)-JVS(602)*XX(75)-JVS(752)*XX(79)
  XX(37) = XX(37)-JVS(286)*XX(58)-JVS(601)*XX(75)-JVS(751)*XX(79)
  XX(36) = XX(36)-JVS(484)*XX(72)-JVS(524)*XX(73)-JVS(600)*XX(75)-JVS(679)*XX(77)-JVS(709)*XX(78)-JVS(797)*XX(80)
  XX(35) = XX(35)-JVS(220)*XX(50)-JVS(333)*XX(62)-JVS(367)*XX(65)-JVS(483)*XX(72)-JVS(599)*XX(75)-JVS(750)*XX(79)
  XX(34) = XX(34)-JVS(598)*XX(75)-JVS(678)*XX(77)
  XX(33) = XX(33)-JVS(213)*XX(49)-JVS(265)*XX(55)-JVS(597)*XX(75)-JVS(708)*XX(78)
  XX(32) = XX(32)-JVS(482)*XX(72)-JVS(596)*XX(75)-JVS(677)*XX(77)-JVS(707)*XX(78)
  XX(31) = XX(31)-JVS(140)*XX(37)-JVS(183)*XX(44)-JVS(189)*XX(45)-JVS(402)*XX(67)-JVS(595)*XX(75)
  XX(30) = XX(30)-JVS(247)*XX(53)-JVS(481)*XX(72)-JVS(523)*XX(73)-JVS(648)*XX(76)-JVS(749)*XX(79)
  XX(29) = XX(29)-JVS(401)*XX(67)-JVS(522)*XX(73)-JVS(594)*XX(75)
  XX(28) = XX(28)-JVS(426)*XX(68)-JVS(593)*XX(75)
  XX(27) = XX(27)-JVS(103)*XX(30)-JVS(122)*XX(34)-JVS(285)*XX(58)-JVS(480)*XX(72)-JVS(592)*XX(75)-JVS(748)*XX(79)&
             &-JVS(821)*XX(81)
  XX(26) = XX(26)-JVS(521)*XX(73)-JVS(591)*XX(75)-JVS(747)*XX(79)-JVS(796)*XX(80)
  XX(25) = XX(25)-JVS(212)*XX(49)-JVS(332)*XX(62)-JVS(590)*XX(75)-JVS(706)*XX(78)
  XX(24) = XX(24)-JVS(227)*XX(51)-JVS(331)*XX(62)-JVS(479)*XX(72)-JVS(520)*XX(73)-JVS(647)*XX(76)-JVS(746)*XX(79)
  XX(23) = XX(23)-JVS(349)*XX(63)-JVS(478)*XX(72)-JVS(589)*XX(75)-JVS(705)*XX(78)
  XX(22) = XX(22)-JVS(226)*XX(51)-JVS(279)*XX(57)-JVS(588)*XX(75)-JVS(745)*XX(79)
  XX(21) = XX(21)-JVS(264)*XX(55)-JVS(330)*XX(62)-JVS(587)*XX(75)-JVS(744)*XX(79)
  XX(20) = XX(20)-JVS(454)*XX(70)-JVS(557)*XX(74)-JVS(586)*XX(75)-JVS(743)*XX(79)
  XX(19) = XX(19)-JVS(201)*XX(47)-JVS(585)*XX(75)-JVS(742)*XX(79)
  XX(18) = XX(18)-JVS(477)*XX(72)-JVS(584)*XX(75)-JVS(676)*XX(77)-JVS(741)*XX(79)
  XX(17) = XX(17)-JVS(329)*XX(62)-JVS(583)*XX(75)-JVS(740)*XX(79)
  XX(16) = XX(16)-JVS(153)*XX(39)-JVS(582)*XX(75)-JVS(739)*XX(79)-JVS(795)*XX(80)
  XX(15) = XX(15)-JVS(425)*XX(68)-JVS(581)*XX(75)
  XX(14) = XX(14)-JVS(188)*XX(45)-JVS(519)*XX(73)-JVS(738)*XX(79)
  XX(13) = XX(13)-JVS(152)*XX(39)-JVS(518)*XX(73)-JVS(794)*XX(80)
  XX(12) = XX(12)-JVS(246)*XX(53)-JVS(476)*XX(72)-JVS(737)*XX(79)
  XX(11) = XX(11)-JVS(34)*XX(14)-JVS(580)*XX(75)
  XX(10) = XX(10)-JVS(25)*XX(11)-JVS(182)*XX(44)-JVS(579)*XX(75)-JVS(736)*XX(79)
  XX(9) = XX(9)-JVS(578)*XX(75)-JVS(735)*XX(79)
  XX(8) = XX(8)-JVS(577)*XX(75)
  XX(7) = XX(7)-JVS(257)*XX(54)-JVS(576)*XX(75)
  XX(6) = XX(6)-JVS(74)*XX(24)-JVS(575)*XX(75)
  XX(5) = XX(5)-JVS(278)*XX(57)-JVS(574)*XX(75)
  XX(4) = XX(4)-JVS(263)*XX(55)-JVS(573)*XX(75)
  XX(3) = XX(3)-JVS(200)*XX(47)-JVS(646)*XX(76)
  XX(2) = XX(2)-JVS(572)*XX(75)
  XX(1) = XX(1)
      
END SUBROUTINE KppSolveTR

! End of KppSolveTR function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



END MODULE mozart_LinearAlgebra

