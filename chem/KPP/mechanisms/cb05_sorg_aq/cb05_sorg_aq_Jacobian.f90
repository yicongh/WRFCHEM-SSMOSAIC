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
! File                 : cb05_sorg_aq_Jacobian.f90
! Time                 : Tue Jan 11 14:32:38 2022
! Working directory    : /jathar-scratch/yicongh/WRF_GoAmazon_gasaerochem.DEV/chem/KPP/mechanisms/cb05_sorg_aq
! Equation file        : cb05_sorg_aq.kpp
! Output root filename : cb05_sorg_aq
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE cb05_sorg_aq_Jacobian

  USE cb05_sorg_aq_Parameters
  USE cb05_sorg_aq_JacobianSP

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

  JUV(1) = JVS(1)*UV(1)+JVS(2)*UV(43)+JVS(3)*UV(96)
  JUV(2) = JVS(4)*UV(2)+JVS(5)*UV(43)+JVS(6)*UV(96)
  JUV(3) = JVS(7)*UV(3)+JVS(8)*UV(73)+JVS(9)*UV(96)+JVS(10)*UV(97)
  JUV(4) = JVS(11)*UV(4)+JVS(12)*UV(47)+JVS(13)*UV(96)
  JUV(5) = JVS(14)*UV(5)+JVS(15)*UV(47)+JVS(16)*UV(96)
  JUV(6) = JVS(17)*UV(6)+JVS(18)*UV(80)+JVS(19)*UV(86)+JVS(20)*UV(93)+JVS(21)*UV(96)+JVS(22)*UV(97)
  JUV(7) = JVS(23)*UV(7)+JVS(24)*UV(80)+JVS(25)*UV(86)+JVS(26)*UV(93)+JVS(27)*UV(96)+JVS(28)*UV(97)
  JUV(8) = JVS(29)*UV(8)+JVS(30)*UV(40)+JVS(31)*UV(96)
  JUV(9) = JVS(32)*UV(9)+JVS(33)*UV(40)+JVS(34)*UV(96)
  JUV(10) = JVS(35)*UV(10)+JVS(36)*UV(75)+JVS(37)*UV(86)+JVS(38)*UV(93)+JVS(39)*UV(96)+JVS(40)*UV(97)
  JUV(11) = JVS(41)*UV(11)+JVS(42)*UV(12)+JVS(43)*UV(96)
  JUV(12) = JVS(44)*UV(12)+JVS(45)*UV(96)
  JUV(13) = JVS(46)*UV(13)+JVS(47)*UV(15)+JVS(48)*UV(96)
  JUV(14) = JVS(49)*UV(14)+JVS(50)*UV(15)+JVS(51)*UV(96)
  JUV(15) = JVS(52)*UV(15)+JVS(53)*UV(96)
  JUV(16) = JVS(54)*UV(16)+JVS(55)*UV(18)+JVS(56)*UV(96)
  JUV(17) = JVS(57)*UV(17)+JVS(58)*UV(18)+JVS(59)*UV(96)
  JUV(18) = JVS(60)*UV(18)+JVS(61)*UV(96)
  JUV(19) = JVS(62)*UV(19)+JVS(63)*UV(23)+JVS(64)*UV(96)
  JUV(20) = JVS(65)*UV(20)+JVS(66)*UV(23)+JVS(67)*UV(96)
  JUV(21) = JVS(68)*UV(21)+JVS(69)*UV(23)+JVS(70)*UV(86)
  JUV(22) = JVS(71)*UV(22)+JVS(72)*UV(23)+JVS(73)*UV(86)
  JUV(23) = JVS(74)*UV(23)+JVS(75)*UV(86)+JVS(76)*UV(96)
  JUV(24) = JVS(77)*UV(24)+JVS(78)*UV(29)+JVS(79)*UV(96)
  JUV(25) = JVS(80)*UV(25)+JVS(81)*UV(29)+JVS(82)*UV(96)
  JUV(26) = JVS(83)*UV(26)+JVS(84)*UV(29)+JVS(85)*UV(86)
  JUV(27) = JVS(86)*UV(27)+JVS(87)*UV(29)+JVS(88)*UV(86)
  JUV(28) = JVS(89)*UV(28)+JVS(90)*UV(29)+JVS(91)*UV(97)
  JUV(29) = JVS(92)*UV(29)+JVS(93)*UV(86)+JVS(94)*UV(96)+JVS(95)*UV(97)
  JUV(30) = JVS(96)*UV(30)+JVS(97)*UV(32)+JVS(98)*UV(96)
  JUV(31) = JVS(99)*UV(31)+JVS(100)*UV(32)+JVS(101)*UV(96)
  JUV(32) = JVS(102)*UV(32)+JVS(103)*UV(96)
  JUV(33) = JVS(104)*UV(33)+JVS(105)*UV(34)+JVS(106)*UV(96)
  JUV(34) = JVS(107)*UV(34)+JVS(108)*UV(96)
  JUV(35) = JVS(109)*UV(35)+JVS(110)*UV(37)+JVS(111)*UV(96)
  JUV(36) = JVS(112)*UV(36)+JVS(113)*UV(37)+JVS(114)*UV(96)
  JUV(37) = JVS(115)*UV(37)+JVS(116)*UV(96)
  JUV(38) = JVS(117)*UV(38)+JVS(118)*UV(49)+JVS(119)*UV(60)+JVS(120)*UV(86)+JVS(121)*UV(96)
  JUV(39) = JVS(122)*UV(39)+JVS(123)*UV(64)
  JUV(40) = JVS(124)*UV(40)+JVS(125)*UV(96)
  JUV(41) = JVS(126)*UV(41)+JVS(127)*UV(87)+JVS(128)*UV(91)
  JUV(42) = JVS(129)*UV(42)+JVS(130)*UV(64)+JVS(131)*UV(90)
  JUV(43) = JVS(132)*UV(43)+JVS(133)*UV(96)
  JUV(44) = JVS(134)*UV(44)+JVS(135)*UV(45)+JVS(136)*UV(96)
  JUV(45) = JVS(137)*UV(44)+JVS(138)*UV(45)+JVS(139)*UV(86)
  JUV(46) = JVS(141)*UV(46)+JVS(142)*UV(87)+JVS(143)*UV(97)
  JUV(47) = JVS(144)*UV(47)+JVS(145)*UV(96)
  JUV(48) = JVS(146)*UV(48)+JVS(147)*UV(88)+JVS(148)*UV(96)
  JUV(49) = JVS(149)*UV(49)+JVS(150)*UV(60)+JVS(151)*UV(86)+JVS(152)*UV(96)
  JUV(50) = JVS(153)*UV(50)+JVS(154)*UV(87)+JVS(155)*UV(95)+JVS(156)*UV(96)
  JUV(51) = JVS(157)*UV(51)+JVS(158)*UV(62)+JVS(159)*UV(74)+JVS(160)*UV(86)+JVS(161)*UV(95)+JVS(162)*UV(96)
  JUV(52) = JVS(163)*UV(43)+JVS(164)*UV(47)+JVS(165)*UV(52)+JVS(166)*UV(95)+JVS(167)*UV(96)
  JUV(53) = JVS(168)*UV(53)+JVS(169)*UV(87)+JVS(170)*UV(90)+JVS(171)*UV(96)
  JUV(54) = JVS(172)*UV(54)+JVS(173)*UV(90)+JVS(174)*UV(91)+JVS(175)*UV(94)+JVS(176)*UV(96)
  JUV(55) = JVS(177)*UV(55)+JVS(178)*UV(89)+JVS(179)*UV(90)+JVS(180)*UV(91)+JVS(181)*UV(92)+JVS(182)*UV(94)+JVS(183)&
              &*UV(96)
  JUV(56) = JVS(184)*UV(56)+JVS(185)*UV(87)+JVS(186)*UV(94)+JVS(187)*UV(96)
  JUV(57) = JVS(188)*UV(57)+JVS(189)*UV(88)+JVS(190)*UV(96)
  JUV(58) = JVS(191)*UV(58)+JVS(192)*UV(88)+JVS(193)*UV(92)+JVS(194)*UV(96)
  JUV(59) = JVS(195)*UV(59)+JVS(196)*UV(88)+JVS(197)*UV(96)
  JUV(60) = JVS(198)*UV(49)+JVS(199)*UV(60)+JVS(201)*UV(90)+JVS(202)*UV(93)+JVS(203)*UV(96)
  JUV(61) = JVS(204)*UV(48)+JVS(205)*UV(57)+JVS(206)*UV(58)+JVS(207)*UV(59)+JVS(208)*UV(61)+JVS(209)*UV(76)+JVS(210)&
              &*UV(79)+JVS(211)*UV(80)+JVS(212)*UV(82)+JVS(213)*UV(83)+JVS(214)*UV(85)+JVS(215)*UV(88)+JVS(217)*UV(96)
  JUV(62) = JVS(218)*UV(62)+JVS(219)*UV(85)+JVS(220)*UV(90)+JVS(221)*UV(95)
  JUV(63) = JVS(222)*UV(63)+JVS(223)*UV(73)+JVS(224)*UV(87)+JVS(225)*UV(90)+JVS(226)*UV(96)+JVS(227)*UV(97)
  JUV(64) = JVS(228)*UV(64)+JVS(229)*UV(86)+JVS(230)*UV(88)+JVS(231)*UV(90)+JVS(232)*UV(95)
  JUV(65) = JVS(233)*UV(65)+JVS(234)*UV(78)+JVS(235)*UV(89)+JVS(236)*UV(90)+JVS(237)*UV(96)
  JUV(66) = JVS(238)*UV(47)+JVS(239)*UV(66)+JVS(240)*UV(71)+JVS(241)*UV(81)+JVS(242)*UV(86)+JVS(243)*UV(96)
  JUV(67) = JVS(244)*UV(67)+JVS(245)*UV(74)+JVS(246)*UV(76)+JVS(247)*UV(77)+JVS(248)*UV(80)+JVS(249)*UV(88)+JVS(250)&
              &*UV(96)
  JUV(68) = JVS(251)*UV(62)+JVS(252)*UV(68)+JVS(254)*UV(90)+JVS(255)*UV(92)+JVS(257)*UV(96)
  JUV(69) = JVS(258)*UV(46)+JVS(259)*UV(69)+JVS(260)*UV(73)+JVS(261)*UV(81)+JVS(262)*UV(82)+JVS(263)*UV(83)+JVS(264)&
              &*UV(84)+JVS(265)*UV(85)+JVS(266)*UV(87)+JVS(267)*UV(90)+JVS(268)*UV(96)+JVS(269)*UV(97)
  JUV(70) = JVS(270)*UV(66)+JVS(271)*UV(67)+JVS(272)*UV(70)+JVS(273)*UV(71)+JVS(274)*UV(74)+JVS(275)*UV(75)+JVS(276)&
              &*UV(76)+JVS(277)*UV(77)+JVS(278)*UV(80)+JVS(279)*UV(81)+JVS(280)*UV(82)+JVS(281)*UV(83)+JVS(282)*UV(85)&
              &+JVS(283)*UV(86)+JVS(284)*UV(88)+JVS(285)*UV(93)+JVS(286)*UV(96)+JVS(287)*UV(97)
  JUV(71) = JVS(288)*UV(52)+JVS(289)*UV(71)+JVS(290)*UV(73)+JVS(291)*UV(86)+JVS(292)*UV(95)+JVS(293)*UV(96)
  JUV(72) = JVS(294)*UV(72)+JVS(295)*UV(79)+JVS(296)*UV(87)+JVS(297)*UV(88)+JVS(298)*UV(96)
  JUV(73) = JVS(299)*UV(43)+JVS(300)*UV(47)+JVS(301)*UV(52)+JVS(302)*UV(63)+JVS(303)*UV(73)+JVS(305)*UV(90)+JVS(307)&
              &*UV(96)+JVS(308)*UV(97)
  JUV(74) = JVS(309)*UV(74)+JVS(310)*UV(86)+JVS(311)*UV(88)+JVS(312)*UV(93)+JVS(313)*UV(96)+JVS(314)*UV(97)
  JUV(75) = JVS(315)*UV(75)+JVS(316)*UV(86)+JVS(317)*UV(93)+JVS(318)*UV(96)+JVS(319)*UV(97)
  JUV(76) = JVS(320)*UV(76)+JVS(321)*UV(86)+JVS(322)*UV(88)+JVS(323)*UV(93)+JVS(324)*UV(96)+JVS(325)*UV(97)
  JUV(77) = JVS(326)*UV(76)+JVS(327)*UV(77)+JVS(328)*UV(86)+JVS(329)*UV(88)+JVS(330)*UV(93)+JVS(331)*UV(96)+JVS(332)&
              &*UV(97)
  JUV(78) = JVS(333)*UV(57)+JVS(334)*UV(72)+JVS(335)*UV(75)+JVS(336)*UV(77)+JVS(337)*UV(78)+JVS(338)*UV(79)+JVS(339)&
              &*UV(80)+JVS(340)*UV(86)+JVS(342)*UV(88)+JVS(343)*UV(89)+JVS(344)*UV(90)+JVS(345)*UV(93)+JVS(346)*UV(95)&
              &+JVS(347)*UV(96)+JVS(348)*UV(97)
  JUV(79) = JVS(349)*UV(47)+JVS(350)*UV(72)+JVS(351)*UV(75)+JVS(352)*UV(76)+JVS(353)*UV(77)+JVS(354)*UV(79)+JVS(355)&
              &*UV(80)+JVS(356)*UV(81)+JVS(357)*UV(84)+JVS(358)*UV(86)+JVS(359)*UV(87)+JVS(360)*UV(88)+JVS(361)*UV(93)&
              &+JVS(362)*UV(96)+JVS(363)*UV(97)
  JUV(80) = JVS(364)*UV(80)+JVS(365)*UV(86)+JVS(366)*UV(87)+JVS(367)*UV(88)+JVS(368)*UV(93)+JVS(369)*UV(96)+JVS(370)&
              &*UV(97)
  JUV(81) = JVS(371)*UV(80)+JVS(372)*UV(81)+JVS(373)*UV(86)+JVS(374)*UV(87)+JVS(375)*UV(88)+JVS(376)*UV(93)+JVS(377)&
              &*UV(96)+JVS(378)*UV(97)
  JUV(82) = JVS(379)*UV(59)+JVS(380)*UV(65)+JVS(381)*UV(71)+JVS(382)*UV(72)+JVS(384)*UV(74)+JVS(385)*UV(75)+JVS(386)&
              &*UV(76)+JVS(387)*UV(77)+JVS(389)*UV(79)+JVS(390)*UV(80)+JVS(391)*UV(81)+JVS(392)*UV(82)+JVS(393)*UV(84)&
              &+JVS(394)*UV(86)+JVS(395)*UV(87)+JVS(396)*UV(88)+JVS(399)*UV(93)+JVS(401)*UV(96)+JVS(402)*UV(97)
  JUV(83) = JVS(403)*UV(56)+JVS(404)*UV(57)+JVS(405)*UV(59)+JVS(406)*UV(65)+JVS(407)*UV(72)+JVS(408)*UV(76)+JVS(409)&
              &*UV(77)+JVS(411)*UV(79)+JVS(413)*UV(81)+JVS(414)*UV(83)+JVS(415)*UV(84)+JVS(416)*UV(86)+JVS(418)*UV(88)&
              &+JVS(419)*UV(89)+JVS(421)*UV(91)+JVS(422)*UV(92)+JVS(423)*UV(93)+JVS(424)*UV(94)+JVS(425)*UV(95)+JVS(426)&
              &*UV(96)+JVS(427)*UV(97)
  JUV(84) = JVS(428)*UV(52)+JVS(429)*UV(63)+JVS(430)*UV(72)+JVS(432)*UV(75)+JVS(433)*UV(78)+JVS(435)*UV(80)+JVS(436)&
              &*UV(81)+JVS(437)*UV(84)+JVS(439)*UV(87)+JVS(444)*UV(95)+JVS(445)*UV(96)+JVS(446)*UV(97)
  JUV(85) = JVS(447)*UV(58)+JVS(448)*UV(59)+JVS(449)*UV(62)+JVS(450)*UV(68)+JVS(451)*UV(71)+JVS(453)*UV(74)+JVS(454)&
              &*UV(75)+JVS(455)*UV(76)+JVS(456)*UV(77)+JVS(457)*UV(80)+JVS(458)*UV(81)+JVS(459)*UV(84)+JVS(460)*UV(85)&
              &+JVS(461)*UV(86)+JVS(463)*UV(88)+JVS(465)*UV(90)+JVS(466)*UV(91)+JVS(467)*UV(92)+JVS(468)*UV(93)+JVS(469)&
              &*UV(94)+JVS(470)*UV(95)+JVS(471)*UV(96)+JVS(472)*UV(97)
  JUV(86) = JVS(473)*UV(49)+JVS(475)*UV(71)+JVS(477)*UV(74)+JVS(478)*UV(75)+JVS(479)*UV(76)+JVS(480)*UV(77)+JVS(481)&
              &*UV(80)+JVS(482)*UV(81)+JVS(483)*UV(86)+JVS(484)*UV(87)+JVS(485)*UV(88)+JVS(486)*UV(90)+JVS(487)*UV(91)&
              &+JVS(488)*UV(93)+JVS(489)*UV(94)+JVS(490)*UV(95)+JVS(491)*UV(96)+JVS(492)*UV(97)
  JUV(87) = JVS(493)*UV(41)+JVS(494)*UV(46)+JVS(495)*UV(50)+JVS(496)*UV(52)+JVS(497)*UV(53)+JVS(498)*UV(56)+JVS(499)&
              &*UV(62)+JVS(500)*UV(63)+JVS(501)*UV(64)+JVS(502)*UV(69)+JVS(503)*UV(72)+JVS(505)*UV(74)+JVS(506)*UV(75)&
              &+JVS(507)*UV(76)+JVS(508)*UV(77)+JVS(510)*UV(80)+JVS(514)*UV(84)+JVS(516)*UV(86)+JVS(517)*UV(87)+JVS(519)&
              &*UV(89)+JVS(520)*UV(90)+JVS(521)*UV(91)+JVS(522)*UV(92)+JVS(523)*UV(93)+JVS(524)*UV(94)+JVS(525)*UV(95)&
              &+JVS(526)*UV(96)+JVS(527)*UV(97)
  JUV(88) = JVS(528)*UV(39)+JVS(529)*UV(42)+JVS(530)*UV(48)+JVS(531)*UV(57)+JVS(532)*UV(58)+JVS(533)*UV(59)+JVS(534)&
              &*UV(61)+JVS(535)*UV(64)+JVS(536)*UV(67)+JVS(537)*UV(74)+JVS(538)*UV(76)+JVS(539)*UV(77)+JVS(540)*UV(79)&
              &+JVS(541)*UV(80)+JVS(543)*UV(82)+JVS(544)*UV(83)+JVS(546)*UV(85)+JVS(547)*UV(86)+JVS(549)*UV(88)+JVS(556)&
              &*UV(95)+JVS(557)*UV(96)
  JUV(89) = JVS(559)*UV(43)+JVS(560)*UV(47)+JVS(561)*UV(57)+JVS(562)*UV(59)+JVS(563)*UV(65)+JVS(564)*UV(66)+JVS(565)&
              &*UV(68)+JVS(566)*UV(71)+JVS(567)*UV(72)+JVS(568)*UV(73)+JVS(569)*UV(74)+JVS(570)*UV(75)+JVS(571)*UV(76)&
              &+JVS(572)*UV(77)+JVS(573)*UV(78)+JVS(574)*UV(79)+JVS(575)*UV(80)+JVS(576)*UV(81)+JVS(579)*UV(86)+JVS(580)&
              &*UV(87)+JVS(581)*UV(88)+JVS(582)*UV(89)+JVS(583)*UV(90)+JVS(584)*UV(91)+JVS(585)*UV(92)+JVS(586)*UV(93)&
              &+JVS(587)*UV(94)+JVS(588)*UV(95)+JVS(589)*UV(96)+JVS(590)*UV(97)
  JUV(90) = JVS(591)*UV(40)+JVS(592)*UV(43)+JVS(593)*UV(44)+JVS(594)*UV(45)+JVS(595)*UV(47)+JVS(596)*UV(51)+JVS(597)&
              &*UV(52)+JVS(598)*UV(53)+JVS(599)*UV(57)+JVS(600)*UV(58)+JVS(601)*UV(59)+JVS(602)*UV(60)+JVS(603)*UV(62)&
              &+JVS(604)*UV(63)+JVS(605)*UV(64)+JVS(606)*UV(65)+JVS(607)*UV(66)+JVS(608)*UV(67)+JVS(609)*UV(68)+JVS(610)&
              &*UV(70)+JVS(611)*UV(71)+JVS(612)*UV(72)+JVS(613)*UV(73)+JVS(614)*UV(74)+JVS(615)*UV(75)+JVS(616)*UV(76)&
              &+JVS(617)*UV(77)+JVS(618)*UV(78)+JVS(619)*UV(79)+JVS(620)*UV(80)+JVS(621)*UV(81)+JVS(622)*UV(82)+JVS(623)&
              &*UV(83)+JVS(624)*UV(84)+JVS(625)*UV(85)+JVS(626)*UV(86)+JVS(627)*UV(87)+JVS(628)*UV(88)+JVS(629)*UV(89)&
              &+JVS(630)*UV(90)+JVS(631)*UV(91)+JVS(632)*UV(92)+JVS(633)*UV(93)+JVS(634)*UV(94)+JVS(635)*UV(95)+JVS(636)&
              &*UV(96)+JVS(637)*UV(97)
  JUV(91) = JVS(638)*UV(41)+JVS(639)*UV(54)+JVS(640)*UV(66)+JVS(641)*UV(71)+JVS(643)*UV(81)+JVS(644)*UV(83)+JVS(646)&
              &*UV(86)+JVS(647)*UV(87)+JVS(648)*UV(88)+JVS(649)*UV(89)+JVS(650)*UV(90)+JVS(651)*UV(91)+JVS(652)*UV(92)&
              &+JVS(653)*UV(93)+JVS(654)*UV(94)+JVS(655)*UV(95)+JVS(656)*UV(96)+JVS(657)*UV(97)
  JUV(92) = JVS(658)*UV(48)+JVS(659)*UV(54)+JVS(660)*UV(55)+JVS(661)*UV(68)+JVS(662)*UV(82)+JVS(663)*UV(83)+JVS(668)&
              &*UV(88)+JVS(669)*UV(89)+JVS(670)*UV(90)+JVS(671)*UV(91)+JVS(672)*UV(92)+JVS(674)*UV(94)+JVS(675)*UV(95)&
              &+JVS(676)*UV(96)
  JUV(93) = JVS(678)*UV(45)+JVS(679)*UV(60)+JVS(680)*UV(74)+JVS(681)*UV(75)+JVS(682)*UV(76)+JVS(683)*UV(77)+JVS(684)&
              &*UV(80)+JVS(685)*UV(82)+JVS(686)*UV(83)+JVS(688)*UV(85)+JVS(689)*UV(86)+JVS(690)*UV(87)+JVS(693)*UV(90)&
              &+JVS(696)*UV(93)+JVS(698)*UV(95)+JVS(699)*UV(96)+JVS(700)*UV(97)
  JUV(94) = JVS(701)*UV(56)+JVS(702)*UV(75)+JVS(703)*UV(80)+JVS(704)*UV(81)+JVS(705)*UV(82)+JVS(707)*UV(86)+JVS(708)&
              &*UV(87)+JVS(709)*UV(88)+JVS(710)*UV(89)+JVS(711)*UV(90)+JVS(712)*UV(91)+JVS(713)*UV(92)+JVS(714)*UV(93)&
              &+JVS(715)*UV(94)+JVS(716)*UV(95)+JVS(717)*UV(96)+JVS(718)*UV(97)
  JUV(95) = JVS(719)*UV(50)+JVS(720)*UV(52)+JVS(721)*UV(62)+JVS(722)*UV(64)+JVS(723)*UV(78)+JVS(725)*UV(80)+JVS(729)&
              &*UV(86)+JVS(730)*UV(87)+JVS(732)*UV(89)+JVS(733)*UV(90)+JVS(734)*UV(91)+JVS(735)*UV(92)+JVS(736)*UV(93)&
              &+JVS(737)*UV(94)+JVS(738)*UV(95)+JVS(739)*UV(96)+JVS(740)*UV(97)
  JUV(96) = JVS(741)*UV(40)+JVS(742)*UV(42)+JVS(743)*UV(43)+JVS(744)*UV(44)+JVS(745)*UV(45)+JVS(746)*UV(47)+JVS(747)&
              &*UV(48)+JVS(748)*UV(49)+JVS(749)*UV(50)+JVS(750)*UV(51)+JVS(751)*UV(53)+JVS(752)*UV(54)+JVS(753)*UV(55)&
              &+JVS(754)*UV(56)+JVS(755)*UV(57)+JVS(756)*UV(58)+JVS(757)*UV(59)+JVS(758)*UV(60)+JVS(759)*UV(61)+JVS(762)&
              &*UV(65)+JVS(763)*UV(66)+JVS(764)*UV(67)+JVS(765)*UV(68)+JVS(766)*UV(69)+JVS(767)*UV(70)+JVS(768)*UV(71)&
              &+JVS(769)*UV(73)+JVS(770)*UV(74)+JVS(771)*UV(75)+JVS(772)*UV(76)+JVS(773)*UV(77)+JVS(775)*UV(79)+JVS(776)&
              &*UV(80)+JVS(777)*UV(81)+JVS(778)*UV(82)+JVS(779)*UV(83)+JVS(780)*UV(84)+JVS(781)*UV(85)+JVS(782)*UV(86)&
              &+JVS(783)*UV(87)+JVS(786)*UV(90)+JVS(789)*UV(93)+JVS(791)*UV(95)+JVS(792)*UV(96)+JVS(793)*UV(97)
  JUV(97) = JVS(794)*UV(46)+JVS(795)*UV(53)+JVS(796)*UV(69)+JVS(797)*UV(73)+JVS(798)*UV(74)+JVS(799)*UV(75)+JVS(800)&
              &*UV(76)+JVS(801)*UV(77)+JVS(802)*UV(80)+JVS(803)*UV(81)+JVS(804)*UV(82)+JVS(805)*UV(83)+JVS(807)*UV(85)&
              &+JVS(808)*UV(86)+JVS(809)*UV(87)+JVS(812)*UV(90)+JVS(815)*UV(93)+JVS(817)*UV(95)+JVS(818)*UV(96)+JVS(819)&
              &*UV(97)
      
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
  JTUV(3) = JVS(7)*UV(3)
  JTUV(4) = JVS(11)*UV(4)
  JTUV(5) = JVS(14)*UV(5)
  JTUV(6) = JVS(17)*UV(6)
  JTUV(7) = JVS(23)*UV(7)
  JTUV(8) = JVS(29)*UV(8)
  JTUV(9) = JVS(32)*UV(9)
  JTUV(10) = JVS(35)*UV(10)
  JTUV(11) = JVS(41)*UV(11)
  JTUV(12) = JVS(42)*UV(11)+JVS(44)*UV(12)
  JTUV(13) = JVS(46)*UV(13)
  JTUV(14) = JVS(49)*UV(14)
  JTUV(15) = JVS(47)*UV(13)+JVS(50)*UV(14)+JVS(52)*UV(15)
  JTUV(16) = JVS(54)*UV(16)
  JTUV(17) = JVS(57)*UV(17)
  JTUV(18) = JVS(55)*UV(16)+JVS(58)*UV(17)+JVS(60)*UV(18)
  JTUV(19) = JVS(62)*UV(19)
  JTUV(20) = JVS(65)*UV(20)
  JTUV(21) = JVS(68)*UV(21)
  JTUV(22) = JVS(71)*UV(22)
  JTUV(23) = JVS(63)*UV(19)+JVS(66)*UV(20)+JVS(69)*UV(21)+JVS(72)*UV(22)+JVS(74)*UV(23)
  JTUV(24) = JVS(77)*UV(24)
  JTUV(25) = JVS(80)*UV(25)
  JTUV(26) = JVS(83)*UV(26)
  JTUV(27) = JVS(86)*UV(27)
  JTUV(28) = JVS(89)*UV(28)
  JTUV(29) = JVS(78)*UV(24)+JVS(81)*UV(25)+JVS(84)*UV(26)+JVS(87)*UV(27)+JVS(90)*UV(28)+JVS(92)*UV(29)
  JTUV(30) = JVS(96)*UV(30)
  JTUV(31) = JVS(99)*UV(31)
  JTUV(32) = JVS(97)*UV(30)+JVS(100)*UV(31)+JVS(102)*UV(32)
  JTUV(33) = JVS(104)*UV(33)
  JTUV(34) = JVS(105)*UV(33)+JVS(107)*UV(34)
  JTUV(35) = JVS(109)*UV(35)
  JTUV(36) = JVS(112)*UV(36)
  JTUV(37) = JVS(110)*UV(35)+JVS(113)*UV(36)+JVS(115)*UV(37)
  JTUV(38) = JVS(117)*UV(38)
  JTUV(39) = JVS(122)*UV(39)+JVS(528)*UV(88)
  JTUV(40) = JVS(30)*UV(8)+JVS(33)*UV(9)+JVS(124)*UV(40)+JVS(591)*UV(90)+JVS(741)*UV(96)
  JTUV(41) = JVS(126)*UV(41)+JVS(493)*UV(87)+JVS(638)*UV(91)
  JTUV(42) = JVS(129)*UV(42)+JVS(529)*UV(88)+JVS(742)*UV(96)
  JTUV(43) = JVS(2)*UV(1)+JVS(5)*UV(2)+JVS(132)*UV(43)+JVS(163)*UV(52)+JVS(299)*UV(73)+JVS(559)*UV(89)+JVS(592)*UV(90)&
               &+JVS(743)*UV(96)
  JTUV(44) = JVS(134)*UV(44)+JVS(137)*UV(45)+JVS(593)*UV(90)+JVS(744)*UV(96)
  JTUV(45) = JVS(135)*UV(44)+JVS(138)*UV(45)+JVS(594)*UV(90)+JVS(678)*UV(93)+JVS(745)*UV(96)
  JTUV(46) = JVS(141)*UV(46)+JVS(258)*UV(69)+JVS(494)*UV(87)+JVS(794)*UV(97)
  JTUV(47) = JVS(12)*UV(4)+JVS(15)*UV(5)+JVS(144)*UV(47)+JVS(164)*UV(52)+JVS(238)*UV(66)+JVS(300)*UV(73)+JVS(349)*UV(79)&
               &+JVS(560)*UV(89)+JVS(595)*UV(90)+JVS(746)*UV(96)
  JTUV(48) = JVS(146)*UV(48)+JVS(204)*UV(61)+JVS(530)*UV(88)+JVS(658)*UV(92)+JVS(747)*UV(96)
  JTUV(49) = JVS(118)*UV(38)+JVS(149)*UV(49)+JVS(198)*UV(60)+JVS(473)*UV(86)+JVS(748)*UV(96)
  JTUV(50) = JVS(153)*UV(50)+JVS(495)*UV(87)+JVS(719)*UV(95)+JVS(749)*UV(96)
  JTUV(51) = JVS(157)*UV(51)+JVS(596)*UV(90)+JVS(750)*UV(96)
  JTUV(52) = JVS(165)*UV(52)+JVS(288)*UV(71)+JVS(301)*UV(73)+JVS(428)*UV(84)+JVS(496)*UV(87)+JVS(597)*UV(90)+JVS(720)&
               &*UV(95)
  JTUV(53) = JVS(168)*UV(53)+JVS(497)*UV(87)+JVS(598)*UV(90)+JVS(751)*UV(96)+JVS(795)*UV(97)
  JTUV(54) = JVS(172)*UV(54)+JVS(639)*UV(91)+JVS(659)*UV(92)+JVS(752)*UV(96)
  JTUV(55) = JVS(177)*UV(55)+JVS(660)*UV(92)+JVS(753)*UV(96)
  JTUV(56) = JVS(184)*UV(56)+JVS(403)*UV(83)+JVS(498)*UV(87)+JVS(701)*UV(94)+JVS(754)*UV(96)
  JTUV(57) = JVS(188)*UV(57)+JVS(205)*UV(61)+JVS(333)*UV(78)+JVS(404)*UV(83)+JVS(531)*UV(88)+JVS(561)*UV(89)+JVS(599)&
               &*UV(90)+JVS(755)*UV(96)
  JTUV(58) = JVS(191)*UV(58)+JVS(206)*UV(61)+JVS(447)*UV(85)+JVS(532)*UV(88)+JVS(600)*UV(90)+JVS(756)*UV(96)
  JTUV(59) = JVS(195)*UV(59)+JVS(207)*UV(61)+JVS(379)*UV(82)+JVS(405)*UV(83)+JVS(448)*UV(85)+JVS(533)*UV(88)+JVS(562)&
               &*UV(89)+JVS(601)*UV(90)+JVS(757)*UV(96)
  JTUV(60) = JVS(119)*UV(38)+JVS(150)*UV(49)+JVS(199)*UV(60)+JVS(602)*UV(90)+JVS(679)*UV(93)+JVS(758)*UV(96)
  JTUV(61) = JVS(208)*UV(61)+JVS(534)*UV(88)+JVS(759)*UV(96)
  JTUV(62) = JVS(158)*UV(51)+JVS(218)*UV(62)+JVS(251)*UV(68)+JVS(449)*UV(85)+JVS(499)*UV(87)+JVS(603)*UV(90)+JVS(721)&
               &*UV(95)
  JTUV(63) = JVS(222)*UV(63)+JVS(302)*UV(73)+JVS(429)*UV(84)+JVS(500)*UV(87)+JVS(604)*UV(90)
  JTUV(64) = JVS(123)*UV(39)+JVS(130)*UV(42)+JVS(228)*UV(64)+JVS(501)*UV(87)+JVS(535)*UV(88)+JVS(605)*UV(90)+JVS(722)&
               &*UV(95)
  JTUV(65) = JVS(233)*UV(65)+JVS(380)*UV(82)+JVS(406)*UV(83)+JVS(563)*UV(89)+JVS(606)*UV(90)+JVS(762)*UV(96)
  JTUV(66) = JVS(239)*UV(66)+JVS(270)*UV(70)+JVS(564)*UV(89)+JVS(607)*UV(90)+JVS(640)*UV(91)+JVS(763)*UV(96)
  JTUV(67) = JVS(244)*UV(67)+JVS(271)*UV(70)+JVS(536)*UV(88)+JVS(608)*UV(90)+JVS(764)*UV(96)
  JTUV(68) = JVS(252)*UV(68)+JVS(450)*UV(85)+JVS(565)*UV(89)+JVS(609)*UV(90)+JVS(661)*UV(92)+JVS(765)*UV(96)
  JTUV(69) = JVS(259)*UV(69)+JVS(502)*UV(87)+JVS(766)*UV(96)+JVS(796)*UV(97)
  JTUV(70) = JVS(272)*UV(70)+JVS(610)*UV(90)+JVS(767)*UV(96)
  JTUV(71) = JVS(240)*UV(66)+JVS(273)*UV(70)+JVS(289)*UV(71)+JVS(381)*UV(82)+JVS(451)*UV(85)+JVS(475)*UV(86)+JVS(566)&
               &*UV(89)+JVS(611)*UV(90)+JVS(641)*UV(91)+JVS(768)*UV(96)
  JTUV(72) = JVS(294)*UV(72)+JVS(334)*UV(78)+JVS(350)*UV(79)+JVS(382)*UV(82)+JVS(407)*UV(83)+JVS(430)*UV(84)+JVS(503)&
               &*UV(87)+JVS(567)*UV(89)+JVS(612)*UV(90)
  JTUV(73) = JVS(8)*UV(3)+JVS(223)*UV(63)+JVS(260)*UV(69)+JVS(290)*UV(71)+JVS(303)*UV(73)+JVS(568)*UV(89)+JVS(613)&
               &*UV(90)+JVS(769)*UV(96)+JVS(797)*UV(97)
  JTUV(74) = JVS(159)*UV(51)+JVS(245)*UV(67)+JVS(274)*UV(70)+JVS(309)*UV(74)+JVS(384)*UV(82)+JVS(453)*UV(85)+JVS(477)&
               &*UV(86)+JVS(505)*UV(87)+JVS(537)*UV(88)+JVS(569)*UV(89)+JVS(614)*UV(90)+JVS(680)*UV(93)+JVS(770)*UV(96)&
               &+JVS(798)*UV(97)
  JTUV(75) = JVS(36)*UV(10)+JVS(275)*UV(70)+JVS(315)*UV(75)+JVS(335)*UV(78)+JVS(351)*UV(79)+JVS(385)*UV(82)+JVS(432)&
               &*UV(84)+JVS(454)*UV(85)+JVS(478)*UV(86)+JVS(506)*UV(87)+JVS(570)*UV(89)+JVS(615)*UV(90)+JVS(681)*UV(93)&
               &+JVS(702)*UV(94)+JVS(771)*UV(96)+JVS(799)*UV(97)
  JTUV(76) = JVS(209)*UV(61)+JVS(246)*UV(67)+JVS(276)*UV(70)+JVS(320)*UV(76)+JVS(326)*UV(77)+JVS(352)*UV(79)+JVS(386)&
               &*UV(82)+JVS(408)*UV(83)+JVS(455)*UV(85)+JVS(479)*UV(86)+JVS(507)*UV(87)+JVS(538)*UV(88)+JVS(571)*UV(89)&
               &+JVS(616)*UV(90)+JVS(682)*UV(93)+JVS(772)*UV(96)+JVS(800)*UV(97)
  JTUV(77) = JVS(247)*UV(67)+JVS(277)*UV(70)+JVS(327)*UV(77)+JVS(336)*UV(78)+JVS(353)*UV(79)+JVS(387)*UV(82)+JVS(409)&
               &*UV(83)+JVS(456)*UV(85)+JVS(480)*UV(86)+JVS(508)*UV(87)+JVS(539)*UV(88)+JVS(572)*UV(89)+JVS(617)*UV(90)&
               &+JVS(683)*UV(93)+JVS(773)*UV(96)+JVS(801)*UV(97)
  JTUV(78) = JVS(234)*UV(65)+JVS(337)*UV(78)+JVS(433)*UV(84)+JVS(573)*UV(89)+JVS(618)*UV(90)+JVS(723)*UV(95)
  JTUV(79) = JVS(210)*UV(61)+JVS(295)*UV(72)+JVS(338)*UV(78)+JVS(354)*UV(79)+JVS(389)*UV(82)+JVS(411)*UV(83)+JVS(540)&
               &*UV(88)+JVS(574)*UV(89)+JVS(619)*UV(90)+JVS(775)*UV(96)
  JTUV(80) = JVS(18)*UV(6)+JVS(24)*UV(7)+JVS(211)*UV(61)+JVS(248)*UV(67)+JVS(278)*UV(70)+JVS(339)*UV(78)+JVS(355)*UV(79)&
               &+JVS(364)*UV(80)+JVS(371)*UV(81)+JVS(390)*UV(82)+JVS(435)*UV(84)+JVS(457)*UV(85)+JVS(481)*UV(86)+JVS(510)&
               &*UV(87)+JVS(541)*UV(88)+JVS(575)*UV(89)+JVS(620)*UV(90)+JVS(684)*UV(93)+JVS(703)*UV(94)+JVS(725)*UV(95)&
               &+JVS(776)*UV(96)+JVS(802)*UV(97)
  JTUV(81) = JVS(241)*UV(66)+JVS(261)*UV(69)+JVS(279)*UV(70)+JVS(356)*UV(79)+JVS(372)*UV(81)+JVS(391)*UV(82)+JVS(413)&
               &*UV(83)+JVS(436)*UV(84)+JVS(458)*UV(85)+JVS(482)*UV(86)+JVS(576)*UV(89)+JVS(621)*UV(90)+JVS(643)*UV(91)&
               &+JVS(704)*UV(94)+JVS(777)*UV(96)+JVS(803)*UV(97)
  JTUV(82) = JVS(212)*UV(61)+JVS(262)*UV(69)+JVS(280)*UV(70)+JVS(392)*UV(82)+JVS(543)*UV(88)+JVS(622)*UV(90)+JVS(662)&
               &*UV(92)+JVS(685)*UV(93)+JVS(705)*UV(94)+JVS(778)*UV(96)+JVS(804)*UV(97)
  JTUV(83) = JVS(213)*UV(61)+JVS(263)*UV(69)+JVS(281)*UV(70)+JVS(414)*UV(83)+JVS(544)*UV(88)+JVS(623)*UV(90)+JVS(644)&
               &*UV(91)+JVS(663)*UV(92)+JVS(686)*UV(93)+JVS(779)*UV(96)+JVS(805)*UV(97)
  JTUV(84) = JVS(264)*UV(69)+JVS(357)*UV(79)+JVS(393)*UV(82)+JVS(415)*UV(83)+JVS(437)*UV(84)+JVS(459)*UV(85)+JVS(514)&
               &*UV(87)+JVS(624)*UV(90)+JVS(780)*UV(96)
  JTUV(85) = JVS(214)*UV(61)+JVS(219)*UV(62)+JVS(265)*UV(69)+JVS(282)*UV(70)+JVS(460)*UV(85)+JVS(546)*UV(88)+JVS(625)&
               &*UV(90)+JVS(688)*UV(93)+JVS(781)*UV(96)+JVS(807)*UV(97)
  JTUV(86) = JVS(19)*UV(6)+JVS(25)*UV(7)+JVS(37)*UV(10)+JVS(70)*UV(21)+JVS(73)*UV(22)+JVS(75)*UV(23)+JVS(85)*UV(26)&
               &+JVS(88)*UV(27)+JVS(93)*UV(29)+JVS(120)*UV(38)+JVS(139)*UV(45)+JVS(151)*UV(49)+JVS(160)*UV(51)+JVS(229)&
               &*UV(64)+JVS(242)*UV(66)+JVS(283)*UV(70)+JVS(291)*UV(71)+JVS(310)*UV(74)+JVS(316)*UV(75)+JVS(321)*UV(76)&
               &+JVS(328)*UV(77)+JVS(340)*UV(78)+JVS(358)*UV(79)+JVS(365)*UV(80)+JVS(373)*UV(81)+JVS(394)*UV(82)+JVS(416)&
               &*UV(83)+JVS(461)*UV(85)+JVS(483)*UV(86)+JVS(516)*UV(87)+JVS(547)*UV(88)+JVS(579)*UV(89)+JVS(626)*UV(90)&
               &+JVS(646)*UV(91)+JVS(689)*UV(93)+JVS(707)*UV(94)+JVS(729)*UV(95)+JVS(782)*UV(96)+JVS(808)*UV(97)
  JTUV(87) = JVS(127)*UV(41)+JVS(142)*UV(46)+JVS(154)*UV(50)+JVS(169)*UV(53)+JVS(185)*UV(56)+JVS(224)*UV(63)+JVS(266)&
               &*UV(69)+JVS(296)*UV(72)+JVS(359)*UV(79)+JVS(366)*UV(80)+JVS(374)*UV(81)+JVS(395)*UV(82)+JVS(439)*UV(84)&
               &+JVS(484)*UV(86)+JVS(517)*UV(87)+JVS(580)*UV(89)+JVS(627)*UV(90)+JVS(647)*UV(91)+JVS(690)*UV(93)+JVS(708)&
               &*UV(94)+JVS(730)*UV(95)+JVS(783)*UV(96)+JVS(809)*UV(97)
  JTUV(88) = JVS(147)*UV(48)+JVS(189)*UV(57)+JVS(192)*UV(58)+JVS(196)*UV(59)+JVS(215)*UV(61)+JVS(230)*UV(64)+JVS(249)&
               &*UV(67)+JVS(284)*UV(70)+JVS(297)*UV(72)+JVS(311)*UV(74)+JVS(322)*UV(76)+JVS(329)*UV(77)+JVS(342)*UV(78)&
               &+JVS(360)*UV(79)+JVS(367)*UV(80)+JVS(375)*UV(81)+JVS(396)*UV(82)+JVS(418)*UV(83)+JVS(463)*UV(85)+JVS(485)&
               &*UV(86)+JVS(549)*UV(88)+JVS(581)*UV(89)+JVS(628)*UV(90)+JVS(648)*UV(91)+JVS(668)*UV(92)+JVS(709)*UV(94)
  JTUV(89) = JVS(178)*UV(55)+JVS(235)*UV(65)+JVS(343)*UV(78)+JVS(419)*UV(83)+JVS(519)*UV(87)+JVS(582)*UV(89)+JVS(629)&
               &*UV(90)+JVS(649)*UV(91)+JVS(669)*UV(92)+JVS(710)*UV(94)+JVS(732)*UV(95)
  JTUV(90) = JVS(131)*UV(42)+JVS(170)*UV(53)+JVS(173)*UV(54)+JVS(179)*UV(55)+JVS(201)*UV(60)+JVS(220)*UV(62)+JVS(225)&
               &*UV(63)+JVS(231)*UV(64)+JVS(236)*UV(65)+JVS(254)*UV(68)+JVS(267)*UV(69)+JVS(305)*UV(73)+JVS(344)*UV(78)&
               &+JVS(465)*UV(85)+JVS(486)*UV(86)+JVS(520)*UV(87)+JVS(583)*UV(89)+JVS(630)*UV(90)+JVS(650)*UV(91)+JVS(670)&
               &*UV(92)+JVS(693)*UV(93)+JVS(711)*UV(94)+JVS(733)*UV(95)+JVS(786)*UV(96)+JVS(812)*UV(97)
  JTUV(91) = JVS(128)*UV(41)+JVS(174)*UV(54)+JVS(180)*UV(55)+JVS(421)*UV(83)+JVS(466)*UV(85)+JVS(487)*UV(86)+JVS(521)&
               &*UV(87)+JVS(584)*UV(89)+JVS(631)*UV(90)+JVS(651)*UV(91)+JVS(671)*UV(92)+JVS(712)*UV(94)+JVS(734)*UV(95)
  JTUV(92) = JVS(181)*UV(55)+JVS(193)*UV(58)+JVS(255)*UV(68)+JVS(422)*UV(83)+JVS(467)*UV(85)+JVS(522)*UV(87)+JVS(585)&
               &*UV(89)+JVS(632)*UV(90)+JVS(652)*UV(91)+JVS(672)*UV(92)+JVS(713)*UV(94)+JVS(735)*UV(95)
  JTUV(93) = JVS(20)*UV(6)+JVS(26)*UV(7)+JVS(38)*UV(10)+JVS(202)*UV(60)+JVS(285)*UV(70)+JVS(312)*UV(74)+JVS(317)*UV(75)&
               &+JVS(323)*UV(76)+JVS(330)*UV(77)+JVS(345)*UV(78)+JVS(361)*UV(79)+JVS(368)*UV(80)+JVS(376)*UV(81)+JVS(399)&
               &*UV(82)+JVS(423)*UV(83)+JVS(468)*UV(85)+JVS(488)*UV(86)+JVS(523)*UV(87)+JVS(586)*UV(89)+JVS(633)*UV(90)&
               &+JVS(653)*UV(91)+JVS(696)*UV(93)+JVS(714)*UV(94)+JVS(736)*UV(95)+JVS(789)*UV(96)+JVS(815)*UV(97)
  JTUV(94) = JVS(175)*UV(54)+JVS(182)*UV(55)+JVS(186)*UV(56)+JVS(424)*UV(83)+JVS(469)*UV(85)+JVS(489)*UV(86)+JVS(524)&
               &*UV(87)+JVS(587)*UV(89)+JVS(634)*UV(90)+JVS(654)*UV(91)+JVS(674)*UV(92)+JVS(715)*UV(94)+JVS(737)*UV(95)
  JTUV(95) = JVS(155)*UV(50)+JVS(161)*UV(51)+JVS(166)*UV(52)+JVS(221)*UV(62)+JVS(232)*UV(64)+JVS(292)*UV(71)+JVS(346)&
               &*UV(78)+JVS(425)*UV(83)+JVS(444)*UV(84)+JVS(470)*UV(85)+JVS(490)*UV(86)+JVS(525)*UV(87)+JVS(556)*UV(88)&
               &+JVS(588)*UV(89)+JVS(635)*UV(90)+JVS(655)*UV(91)+JVS(675)*UV(92)+JVS(698)*UV(93)+JVS(716)*UV(94)+JVS(738)&
               &*UV(95)+JVS(791)*UV(96)+JVS(817)*UV(97)
  JTUV(96) = JVS(3)*UV(1)+JVS(6)*UV(2)+JVS(9)*UV(3)+JVS(13)*UV(4)+JVS(16)*UV(5)+JVS(21)*UV(6)+JVS(27)*UV(7)+JVS(31)&
               &*UV(8)+JVS(34)*UV(9)+JVS(39)*UV(10)+JVS(43)*UV(11)+JVS(45)*UV(12)+JVS(48)*UV(13)+JVS(51)*UV(14)+JVS(53)&
               &*UV(15)+JVS(56)*UV(16)+JVS(59)*UV(17)+JVS(61)*UV(18)+JVS(64)*UV(19)+JVS(67)*UV(20)+JVS(76)*UV(23)+JVS(79)&
               &*UV(24)+JVS(82)*UV(25)+JVS(94)*UV(29)+JVS(98)*UV(30)+JVS(101)*UV(31)+JVS(103)*UV(32)+JVS(106)*UV(33)&
               &+JVS(108)*UV(34)+JVS(111)*UV(35)+JVS(114)*UV(36)+JVS(116)*UV(37)+JVS(121)*UV(38)+JVS(125)*UV(40)+JVS(133)&
               &*UV(43)+JVS(136)*UV(44)+JVS(145)*UV(47)+JVS(148)*UV(48)+JVS(152)*UV(49)+JVS(156)*UV(50)+JVS(162)*UV(51)&
               &+JVS(167)*UV(52)+JVS(171)*UV(53)+JVS(176)*UV(54)+JVS(183)*UV(55)+JVS(187)*UV(56)+JVS(190)*UV(57)+JVS(194)&
               &*UV(58)+JVS(197)*UV(59)+JVS(203)*UV(60)+JVS(217)*UV(61)+JVS(226)*UV(63)+JVS(237)*UV(65)+JVS(243)*UV(66)&
               &+JVS(250)*UV(67)+JVS(257)*UV(68)+JVS(268)*UV(69)+JVS(286)*UV(70)+JVS(293)*UV(71)+JVS(298)*UV(72)+JVS(307)&
               &*UV(73)+JVS(313)*UV(74)+JVS(318)*UV(75)+JVS(324)*UV(76)+JVS(331)*UV(77)+JVS(347)*UV(78)+JVS(362)*UV(79)&
               &+JVS(369)*UV(80)+JVS(377)*UV(81)+JVS(401)*UV(82)+JVS(426)*UV(83)+JVS(445)*UV(84)+JVS(471)*UV(85)+JVS(491)&
               &*UV(86)+JVS(526)*UV(87)+JVS(557)*UV(88)+JVS(589)*UV(89)+JVS(636)*UV(90)+JVS(656)*UV(91)+JVS(676)*UV(92)&
               &+JVS(699)*UV(93)+JVS(717)*UV(94)+JVS(739)*UV(95)+JVS(792)*UV(96)+JVS(818)*UV(97)
  JTUV(97) = JVS(10)*UV(3)+JVS(22)*UV(6)+JVS(28)*UV(7)+JVS(40)*UV(10)+JVS(91)*UV(28)+JVS(95)*UV(29)+JVS(143)*UV(46)&
               &+JVS(227)*UV(63)+JVS(269)*UV(69)+JVS(287)*UV(70)+JVS(308)*UV(73)+JVS(314)*UV(74)+JVS(319)*UV(75)+JVS(325)&
               &*UV(76)+JVS(332)*UV(77)+JVS(348)*UV(78)+JVS(363)*UV(79)+JVS(370)*UV(80)+JVS(378)*UV(81)+JVS(402)*UV(82)&
               &+JVS(427)*UV(83)+JVS(446)*UV(84)+JVS(472)*UV(85)+JVS(492)*UV(86)+JVS(527)*UV(87)+JVS(590)*UV(89)+JVS(637)&
               &*UV(90)+JVS(657)*UV(91)+JVS(700)*UV(93)+JVS(718)*UV(94)+JVS(740)*UV(95)+JVS(793)*UV(96)+JVS(819)*UV(97)
      
END SUBROUTINE JacTR_SP_Vec

! End of JacTR_SP_Vec function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



END MODULE cb05_sorg_aq_Jacobian

