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
! File                 : radm2_Jacobian.f90
! Time                 : Tue Jan 11 14:33:14 2022
! Working directory    : /jathar-scratch/yicongh/WRF_GoAmazon_gasaerochem.DEV/chem/KPP/mechanisms/radm2
! Equation file        : radm2.kpp
! Output root filename : radm2
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE radm2_Jacobian

  USE radm2_Parameters
  USE radm2_JacobianSP

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

  JUV(1) = JVS(1)*UV(1)+JVS(2)*UV(5)+JVS(3)*UV(53)
  JUV(2) = JVS(4)*UV(2)+JVS(5)*UV(24)+JVS(6)*UV(27)+JVS(7)*UV(28)+JVS(8)*UV(29)+JVS(9)*UV(46)
  JUV(3) = JVS(10)*UV(3)+JVS(11)*UV(27)+JVS(12)*UV(28)+JVS(13)*UV(29)+JVS(14)*UV(37)+JVS(15)*UV(38)+JVS(16)*UV(40)&
             &+JVS(17)*UV(41)+JVS(18)*UV(44)+JVS(19)*UV(45)+JVS(20)*UV(46)+JVS(21)*UV(49)+JVS(22)*UV(50)+JVS(23)*UV(52)&
             &+JVS(24)*UV(55)+JVS(25)*UV(56)+JVS(26)*UV(57)
  JUV(4) = JVS(27)*UV(4)+JVS(28)*UV(26)+JVS(29)*UV(53)
  JUV(5) = JVS(30)*UV(5)+JVS(31)*UV(53)
  JUV(6) = JVS(32)*UV(6)+JVS(33)*UV(53)
  JUV(7) = JVS(34)*UV(7)+JVS(35)*UV(46)
  JUV(8) = JVS(36)*UV(8)+JVS(37)*UV(53)
  JUV(9) = JVS(38)*UV(9)+JVS(39)*UV(53)
  JUV(10) = JVS(40)*UV(10)+JVS(41)*UV(53)
  JUV(11) = JVS(42)*UV(11)+JVS(43)*UV(53)
  JUV(12) = JVS(44)*UV(12)+JVS(45)*UV(45)+JVS(46)*UV(59)
  JUV(13) = JVS(47)*UV(13)+JVS(48)*UV(53)+JVS(49)*UV(58)
  JUV(14) = JVS(50)*UV(14)+JVS(51)*UV(51)+JVS(52)*UV(53)
  JUV(15) = JVS(53)*UV(15)+JVS(54)*UV(54)+JVS(55)*UV(59)
  JUV(16) = JVS(56)*UV(16)+JVS(57)*UV(53)
  JUV(17) = JVS(58)*UV(17)+JVS(59)*UV(28)+JVS(60)*UV(29)+JVS(61)*UV(46)+JVS(62)*UV(53)
  JUV(18) = JVS(63)*UV(18)+JVS(64)*UV(51)+JVS(65)*UV(53)+JVS(66)*UV(55)
  JUV(19) = JVS(67)*UV(7)+JVS(68)*UV(19)+JVS(69)*UV(46)+JVS(70)*UV(54)+JVS(71)*UV(59)
  JUV(20) = JVS(72)*UV(20)+JVS(73)*UV(51)+JVS(74)*UV(53)+JVS(75)*UV(59)
  JUV(21) = JVS(76)*UV(21)+JVS(77)*UV(51)+JVS(78)*UV(52)+JVS(79)*UV(53)
  JUV(22) = JVS(80)*UV(10)+JVS(81)*UV(11)+JVS(82)*UV(22)+JVS(83)*UV(53)+JVS(84)*UV(54)
  JUV(23) = JVS(85)*UV(23)+JVS(86)*UV(53)+JVS(87)*UV(55)+JVS(88)*UV(59)
  JUV(24) = JVS(89)*UV(24)+JVS(90)*UV(46)+JVS(91)*UV(53)+JVS(92)*UV(54)
  JUV(25) = JVS(93)*UV(15)+JVS(94)*UV(22)+JVS(95)*UV(25)+JVS(96)*UV(30)+JVS(97)*UV(31)+JVS(98)*UV(34)+JVS(99)*UV(43)&
              &+JVS(100)*UV(48)+JVS(101)*UV(51)+JVS(102)*UV(53)+JVS(103)*UV(54)+JVS(104)*UV(59)
  JUV(26) = JVS(105)*UV(24)+JVS(106)*UV(26)+JVS(107)*UV(27)+JVS(108)*UV(28)+JVS(109)*UV(29)+JVS(110)*UV(31)+JVS(111)&
              &*UV(34)+JVS(112)*UV(43)+JVS(113)*UV(45)+JVS(114)*UV(46)+JVS(115)*UV(48)+JVS(116)*UV(52)+JVS(117)*UV(53)&
              &+JVS(118)*UV(54)+JVS(119)*UV(55)+JVS(120)*UV(58)
  JUV(27) = JVS(121)*UV(27)+JVS(122)*UV(46)+JVS(123)*UV(53)+JVS(124)*UV(54)
  JUV(28) = JVS(125)*UV(28)+JVS(126)*UV(46)+JVS(127)*UV(53)+JVS(128)*UV(54)
  JUV(29) = JVS(129)*UV(29)+JVS(130)*UV(46)+JVS(131)*UV(53)+JVS(132)*UV(54)
  JUV(30) = JVS(133)*UV(30)+JVS(134)*UV(35)+JVS(135)*UV(36)+JVS(136)*UV(52)+JVS(137)*UV(53)+JVS(138)*UV(54)+JVS(139)&
              &*UV(55)+JVS(140)*UV(58)
  JUV(31) = JVS(141)*UV(31)+JVS(142)*UV(35)+JVS(143)*UV(45)+JVS(144)*UV(52)+JVS(145)*UV(53)+JVS(146)*UV(54)+JVS(147)&
              &*UV(55)+JVS(148)*UV(58)
  JUV(32) = JVS(149)*UV(22)+JVS(150)*UV(32)+JVS(151)*UV(51)+JVS(152)*UV(52)+JVS(154)*UV(54)+JVS(155)*UV(55)+JVS(156)&
              &*UV(59)
  JUV(33) = JVS(157)*UV(16)+JVS(158)*UV(29)+JVS(159)*UV(33)+JVS(160)*UV(41)+JVS(161)*UV(44)+JVS(162)*UV(46)+JVS(163)&
              &*UV(47)+JVS(164)*UV(49)+JVS(165)*UV(52)+JVS(166)*UV(53)+JVS(168)*UV(55)+JVS(169)*UV(56)+JVS(170)*UV(58)
  JUV(34) = JVS(171)*UV(34)+JVS(172)*UV(35)+JVS(173)*UV(36)+JVS(174)*UV(45)+JVS(175)*UV(50)+JVS(176)*UV(52)+JVS(177)&
              &*UV(53)+JVS(178)*UV(54)+JVS(179)*UV(55)+JVS(180)*UV(58)
  JUV(35) = JVS(181)*UV(10)+JVS(182)*UV(35)+JVS(183)*UV(51)+JVS(184)*UV(52)+JVS(185)*UV(53)+JVS(186)*UV(55)+JVS(187)&
              &*UV(58)
  JUV(36) = JVS(188)*UV(11)+JVS(189)*UV(36)+JVS(190)*UV(51)+JVS(191)*UV(52)+JVS(192)*UV(53)+JVS(193)*UV(55)+JVS(194)&
              &*UV(58)
  JUV(37) = JVS(195)*UV(27)+JVS(196)*UV(28)+JVS(197)*UV(37)+JVS(199)*UV(51)+JVS(200)*UV(52)+JVS(201)*UV(53)+JVS(203)&
              &*UV(55)+JVS(204)*UV(58)
  JUV(38) = JVS(205)*UV(24)+JVS(206)*UV(27)+JVS(207)*UV(28)+JVS(208)*UV(29)+JVS(209)*UV(38)+JVS(211)*UV(51)+JVS(212)&
              &*UV(52)+JVS(214)*UV(54)+JVS(215)*UV(55)+JVS(216)*UV(58)
  JUV(39) = JVS(217)*UV(8)+JVS(218)*UV(9)+JVS(219)*UV(22)+JVS(220)*UV(23)+JVS(221)*UV(39)+JVS(222)*UV(45)+JVS(223)&
              &*UV(51)+JVS(224)*UV(52)+JVS(225)*UV(53)+JVS(227)*UV(55)+JVS(228)*UV(58)
  JUV(40) = JVS(230)*UV(24)+JVS(231)*UV(40)+JVS(233)*UV(51)+JVS(234)*UV(52)+JVS(235)*UV(53)+JVS(237)*UV(55)+JVS(238)&
              &*UV(58)
  JUV(41) = JVS(239)*UV(8)+JVS(240)*UV(41)+JVS(241)*UV(51)+JVS(242)*UV(52)+JVS(243)*UV(53)+JVS(244)*UV(55)+JVS(245)&
              &*UV(58)
  JUV(42) = JVS(246)*UV(32)+JVS(247)*UV(35)+JVS(248)*UV(36)+JVS(249)*UV(37)+JVS(250)*UV(39)+JVS(251)*UV(40)+JVS(252)&
              &*UV(41)+JVS(253)*UV(42)+JVS(254)*UV(44)+JVS(255)*UV(45)+JVS(257)*UV(49)+JVS(258)*UV(50)+JVS(259)*UV(51)&
              &+JVS(261)*UV(53)+JVS(264)*UV(56)+JVS(265)*UV(57)
  JUV(43) = JVS(268)*UV(16)+JVS(269)*UV(21)+JVS(270)*UV(23)+JVS(271)*UV(24)+JVS(272)*UV(27)+JVS(273)*UV(28)+JVS(274)&
              &*UV(29)+JVS(275)*UV(31)+JVS(276)*UV(32)+JVS(277)*UV(35)+JVS(278)*UV(36)+JVS(279)*UV(37)+JVS(280)*UV(38)&
              &+JVS(281)*UV(39)+JVS(282)*UV(40)+JVS(283)*UV(41)+JVS(284)*UV(43)+JVS(285)*UV(44)+JVS(286)*UV(45)+JVS(287)&
              &*UV(46)+JVS(288)*UV(49)+JVS(289)*UV(50)+JVS(291)*UV(52)+JVS(292)*UV(53)+JVS(293)*UV(54)+JVS(294)*UV(55)&
              &+JVS(295)*UV(56)+JVS(296)*UV(57)+JVS(297)*UV(58)
  JUV(44) = JVS(299)*UV(9)+JVS(300)*UV(44)+JVS(301)*UV(51)+JVS(302)*UV(52)+JVS(303)*UV(53)+JVS(304)*UV(55)+JVS(305)&
              &*UV(58)
  JUV(45) = JVS(306)*UV(12)+JVS(307)*UV(22)+JVS(308)*UV(30)+JVS(311)*UV(45)+JVS(312)*UV(51)+JVS(313)*UV(52)+JVS(314)&
              &*UV(53)+JVS(315)*UV(54)+JVS(316)*UV(55)+JVS(317)*UV(58)+JVS(318)*UV(59)
  JUV(46) = JVS(319)*UV(19)+JVS(320)*UV(24)+JVS(321)*UV(27)+JVS(322)*UV(28)+JVS(323)*UV(29)+JVS(324)*UV(46)+JVS(325)&
              &*UV(51)+JVS(326)*UV(53)+JVS(328)*UV(58)+JVS(329)*UV(59)
  JUV(47) = JVS(330)*UV(32)+JVS(331)*UV(38)+JVS(332)*UV(41)+JVS(333)*UV(44)+JVS(335)*UV(47)+JVS(336)*UV(51)+JVS(338)&
              &*UV(53)+JVS(341)*UV(56)+JVS(342)*UV(58)+JVS(343)*UV(59)
  JUV(48) = JVS(344)*UV(16)+JVS(345)*UV(27)+JVS(346)*UV(28)+JVS(347)*UV(29)+JVS(348)*UV(37)+JVS(349)*UV(38)+JVS(350)&
              &*UV(40)+JVS(351)*UV(41)+JVS(352)*UV(42)+JVS(353)*UV(44)+JVS(355)*UV(46)+JVS(356)*UV(47)+JVS(357)*UV(48)&
              &+JVS(358)*UV(49)+JVS(361)*UV(52)+JVS(362)*UV(53)+JVS(363)*UV(54)+JVS(364)*UV(55)+JVS(365)*UV(56)+JVS(366)&
              &*UV(57)+JVS(367)*UV(58)
  JUV(49) = JVS(369)*UV(29)+JVS(371)*UV(49)+JVS(372)*UV(51)+JVS(373)*UV(52)+JVS(374)*UV(53)+JVS(376)*UV(55)+JVS(377)&
              &*UV(58)
  JUV(50) = JVS(379)*UV(33)+JVS(385)*UV(50)+JVS(386)*UV(51)+JVS(387)*UV(52)+JVS(388)*UV(53)+JVS(390)*UV(55)+JVS(392)&
              &*UV(58)
  JUV(51) = JVS(394)*UV(5)+JVS(395)*UV(10)+JVS(396)*UV(11)+JVS(397)*UV(14)+JVS(398)*UV(16)+JVS(399)*UV(20)+JVS(400)&
              &*UV(21)+JVS(401)*UV(22)+JVS(402)*UV(24)+JVS(403)*UV(26)+JVS(404)*UV(27)+JVS(405)*UV(28)+JVS(406)*UV(29)&
              &+JVS(407)*UV(30)+JVS(408)*UV(31)+JVS(409)*UV(32)+JVS(410)*UV(34)+JVS(411)*UV(35)+JVS(412)*UV(36)+JVS(413)&
              &*UV(37)+JVS(414)*UV(38)+JVS(415)*UV(39)+JVS(416)*UV(40)+JVS(417)*UV(41)+JVS(418)*UV(42)+JVS(419)*UV(43)&
              &+JVS(420)*UV(44)+JVS(421)*UV(45)+JVS(422)*UV(46)+JVS(423)*UV(47)+JVS(424)*UV(48)+JVS(425)*UV(49)+JVS(426)&
              &*UV(50)+JVS(427)*UV(51)+JVS(428)*UV(52)+JVS(429)*UV(53)+JVS(430)*UV(54)+JVS(431)*UV(55)+JVS(432)*UV(56)&
              &+JVS(433)*UV(57)+JVS(434)*UV(58)+JVS(435)*UV(59)
  JUV(52) = JVS(436)*UV(17)+JVS(437)*UV(18)+JVS(438)*UV(21)+JVS(439)*UV(27)+JVS(440)*UV(28)+JVS(441)*UV(29)+JVS(442)&
              &*UV(32)+JVS(443)*UV(35)+JVS(444)*UV(36)+JVS(445)*UV(37)+JVS(446)*UV(38)+JVS(447)*UV(39)+JVS(448)*UV(40)&
              &+JVS(449)*UV(41)+JVS(450)*UV(44)+JVS(451)*UV(45)+JVS(452)*UV(46)+JVS(453)*UV(48)+JVS(454)*UV(49)+JVS(455)&
              &*UV(50)+JVS(456)*UV(51)+JVS(457)*UV(52)+JVS(458)*UV(53)+JVS(460)*UV(55)+JVS(461)*UV(56)+JVS(462)*UV(57)&
              &+JVS(463)*UV(58)
  JUV(53) = JVS(465)*UV(5)+JVS(466)*UV(6)+JVS(467)*UV(7)+JVS(468)*UV(8)+JVS(469)*UV(9)+JVS(470)*UV(10)+JVS(471)*UV(11)&
              &+JVS(472)*UV(13)+JVS(473)*UV(14)+JVS(474)*UV(16)+JVS(475)*UV(17)+JVS(476)*UV(18)+JVS(477)*UV(20)+JVS(478)&
              &*UV(21)+JVS(479)*UV(22)+JVS(480)*UV(23)+JVS(481)*UV(24)+JVS(482)*UV(25)+JVS(483)*UV(26)+JVS(484)*UV(27)&
              &+JVS(485)*UV(28)+JVS(486)*UV(29)+JVS(487)*UV(30)+JVS(488)*UV(31)+JVS(489)*UV(33)+JVS(490)*UV(34)+JVS(494)&
              &*UV(42)+JVS(495)*UV(43)+JVS(498)*UV(46)+JVS(499)*UV(47)+JVS(500)*UV(48)+JVS(503)*UV(51)+JVS(505)*UV(53)&
              &+JVS(510)*UV(58)+JVS(511)*UV(59)
  JUV(54) = JVS(512)*UV(15)+JVS(513)*UV(20)+JVS(514)*UV(22)+JVS(515)*UV(23)+JVS(516)*UV(24)+JVS(517)*UV(25)+JVS(518)&
              &*UV(27)+JVS(519)*UV(28)+JVS(520)*UV(29)+JVS(521)*UV(30)+JVS(522)*UV(31)+JVS(523)*UV(34)+JVS(526)*UV(43)&
              &+JVS(529)*UV(46)+JVS(530)*UV(48)+JVS(533)*UV(51)+JVS(535)*UV(53)+JVS(536)*UV(54)+JVS(540)*UV(58)+JVS(541)&
              &*UV(59)
  JUV(55) = JVS(542)*UV(18)+JVS(543)*UV(23)+JVS(544)*UV(32)+JVS(545)*UV(33)+JVS(546)*UV(34)+JVS(547)*UV(35)+JVS(548)&
              &*UV(36)+JVS(549)*UV(37)+JVS(550)*UV(38)+JVS(551)*UV(39)+JVS(552)*UV(40)+JVS(553)*UV(41)+JVS(554)*UV(44)&
              &+JVS(555)*UV(45)+JVS(558)*UV(48)+JVS(559)*UV(49)+JVS(560)*UV(50)+JVS(561)*UV(51)+JVS(562)*UV(52)+JVS(563)&
              &*UV(53)+JVS(564)*UV(54)+JVS(565)*UV(55)+JVS(566)*UV(56)+JVS(567)*UV(57)+JVS(568)*UV(58)+JVS(569)*UV(59)
  JUV(56) = JVS(570)*UV(16)+JVS(571)*UV(42)+JVS(575)*UV(47)+JVS(578)*UV(51)+JVS(579)*UV(52)+JVS(580)*UV(53)+JVS(582)&
              &*UV(55)+JVS(583)*UV(56)+JVS(585)*UV(58)
  JUV(57) = JVS(587)*UV(6)+JVS(588)*UV(33)+JVS(594)*UV(51)+JVS(595)*UV(52)+JVS(596)*UV(53)+JVS(598)*UV(55)+JVS(600)&
              &*UV(57)+JVS(601)*UV(58)
  JUV(58) = JVS(603)*UV(13)+JVS(604)*UV(19)+JVS(605)*UV(35)+JVS(606)*UV(36)+JVS(607)*UV(37)+JVS(608)*UV(38)+JVS(609)&
              &*UV(39)+JVS(610)*UV(40)+JVS(611)*UV(41)+JVS(612)*UV(44)+JVS(613)*UV(45)+JVS(614)*UV(46)+JVS(615)*UV(49)&
              &+JVS(616)*UV(50)+JVS(617)*UV(51)+JVS(618)*UV(52)+JVS(619)*UV(53)+JVS(620)*UV(54)+JVS(621)*UV(55)+JVS(622)&
              &*UV(56)+JVS(623)*UV(57)+JVS(624)*UV(58)+JVS(625)*UV(59)
  JUV(59) = JVS(626)*UV(12)+JVS(627)*UV(15)+JVS(628)*UV(19)+JVS(629)*UV(20)+JVS(630)*UV(23)+JVS(631)*UV(25)+JVS(634)&
              &*UV(32)+JVS(636)*UV(35)+JVS(637)*UV(36)+JVS(638)*UV(37)+JVS(639)*UV(38)+JVS(640)*UV(39)+JVS(641)*UV(40)&
              &+JVS(642)*UV(41)+JVS(644)*UV(44)+JVS(645)*UV(45)+JVS(646)*UV(46)+JVS(647)*UV(47)+JVS(649)*UV(49)+JVS(650)&
              &*UV(50)+JVS(651)*UV(51)+JVS(652)*UV(52)+JVS(653)*UV(53)+JVS(654)*UV(54)+JVS(655)*UV(55)+JVS(656)*UV(56)&
              &+JVS(657)*UV(57)+JVS(658)*UV(58)+JVS(659)*UV(59)
      
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
  JTUV(3) = JVS(10)*UV(3)
  JTUV(4) = JVS(27)*UV(4)
  JTUV(5) = JVS(2)*UV(1)+JVS(30)*UV(5)+JVS(394)*UV(51)+JVS(465)*UV(53)
  JTUV(6) = JVS(32)*UV(6)+JVS(466)*UV(53)+JVS(587)*UV(57)
  JTUV(7) = JVS(34)*UV(7)+JVS(67)*UV(19)+JVS(467)*UV(53)
  JTUV(8) = JVS(36)*UV(8)+JVS(217)*UV(39)+JVS(239)*UV(41)+JVS(468)*UV(53)
  JTUV(9) = JVS(38)*UV(9)+JVS(218)*UV(39)+JVS(299)*UV(44)+JVS(469)*UV(53)
  JTUV(10) = JVS(40)*UV(10)+JVS(80)*UV(22)+JVS(181)*UV(35)+JVS(395)*UV(51)+JVS(470)*UV(53)
  JTUV(11) = JVS(42)*UV(11)+JVS(81)*UV(22)+JVS(188)*UV(36)+JVS(396)*UV(51)+JVS(471)*UV(53)
  JTUV(12) = JVS(44)*UV(12)+JVS(306)*UV(45)+JVS(626)*UV(59)
  JTUV(13) = JVS(47)*UV(13)+JVS(472)*UV(53)+JVS(603)*UV(58)
  JTUV(14) = JVS(50)*UV(14)+JVS(397)*UV(51)+JVS(473)*UV(53)
  JTUV(15) = JVS(53)*UV(15)+JVS(93)*UV(25)+JVS(512)*UV(54)+JVS(627)*UV(59)
  JTUV(16) = JVS(56)*UV(16)+JVS(157)*UV(33)+JVS(268)*UV(43)+JVS(344)*UV(48)+JVS(398)*UV(51)+JVS(474)*UV(53)+JVS(570)&
               &*UV(56)
  JTUV(17) = JVS(58)*UV(17)+JVS(436)*UV(52)+JVS(475)*UV(53)
  JTUV(18) = JVS(63)*UV(18)+JVS(437)*UV(52)+JVS(476)*UV(53)+JVS(542)*UV(55)
  JTUV(19) = JVS(68)*UV(19)+JVS(319)*UV(46)+JVS(604)*UV(58)+JVS(628)*UV(59)
  JTUV(20) = JVS(72)*UV(20)+JVS(399)*UV(51)+JVS(477)*UV(53)+JVS(513)*UV(54)+JVS(629)*UV(59)
  JTUV(21) = JVS(76)*UV(21)+JVS(269)*UV(43)+JVS(400)*UV(51)+JVS(438)*UV(52)+JVS(478)*UV(53)
  JTUV(22) = JVS(82)*UV(22)+JVS(94)*UV(25)+JVS(149)*UV(32)+JVS(219)*UV(39)+JVS(307)*UV(45)+JVS(401)*UV(51)+JVS(479)&
               &*UV(53)+JVS(514)*UV(54)
  JTUV(23) = JVS(85)*UV(23)+JVS(220)*UV(39)+JVS(270)*UV(43)+JVS(480)*UV(53)+JVS(515)*UV(54)+JVS(543)*UV(55)+JVS(630)&
               &*UV(59)
  JTUV(24) = JVS(5)*UV(2)+JVS(89)*UV(24)+JVS(105)*UV(26)+JVS(205)*UV(38)+JVS(230)*UV(40)+JVS(271)*UV(43)+JVS(320)*UV(46)&
               &+JVS(402)*UV(51)+JVS(481)*UV(53)+JVS(516)*UV(54)
  JTUV(25) = JVS(95)*UV(25)+JVS(482)*UV(53)+JVS(517)*UV(54)+JVS(631)*UV(59)
  JTUV(26) = JVS(28)*UV(4)+JVS(106)*UV(26)+JVS(403)*UV(51)+JVS(483)*UV(53)
  JTUV(27) = JVS(6)*UV(2)+JVS(11)*UV(3)+JVS(107)*UV(26)+JVS(121)*UV(27)+JVS(195)*UV(37)+JVS(206)*UV(38)+JVS(272)*UV(43)&
               &+JVS(321)*UV(46)+JVS(345)*UV(48)+JVS(404)*UV(51)+JVS(439)*UV(52)+JVS(484)*UV(53)+JVS(518)*UV(54)
  JTUV(28) = JVS(7)*UV(2)+JVS(12)*UV(3)+JVS(59)*UV(17)+JVS(108)*UV(26)+JVS(125)*UV(28)+JVS(196)*UV(37)+JVS(207)*UV(38)&
               &+JVS(273)*UV(43)+JVS(322)*UV(46)+JVS(346)*UV(48)+JVS(405)*UV(51)+JVS(440)*UV(52)+JVS(485)*UV(53)+JVS(519)&
               &*UV(54)
  JTUV(29) = JVS(8)*UV(2)+JVS(13)*UV(3)+JVS(60)*UV(17)+JVS(109)*UV(26)+JVS(129)*UV(29)+JVS(158)*UV(33)+JVS(208)*UV(38)&
               &+JVS(274)*UV(43)+JVS(323)*UV(46)+JVS(347)*UV(48)+JVS(369)*UV(49)+JVS(406)*UV(51)+JVS(441)*UV(52)+JVS(486)&
               &*UV(53)+JVS(520)*UV(54)
  JTUV(30) = JVS(96)*UV(25)+JVS(133)*UV(30)+JVS(308)*UV(45)+JVS(407)*UV(51)+JVS(487)*UV(53)+JVS(521)*UV(54)
  JTUV(31) = JVS(97)*UV(25)+JVS(110)*UV(26)+JVS(141)*UV(31)+JVS(275)*UV(43)+JVS(408)*UV(51)+JVS(488)*UV(53)+JVS(522)&
               &*UV(54)
  JTUV(32) = JVS(150)*UV(32)+JVS(246)*UV(42)+JVS(276)*UV(43)+JVS(330)*UV(47)+JVS(409)*UV(51)+JVS(442)*UV(52)+JVS(544)&
               &*UV(55)+JVS(634)*UV(59)
  JTUV(33) = JVS(159)*UV(33)+JVS(379)*UV(50)+JVS(489)*UV(53)+JVS(545)*UV(55)+JVS(588)*UV(57)
  JTUV(34) = JVS(98)*UV(25)+JVS(111)*UV(26)+JVS(171)*UV(34)+JVS(410)*UV(51)+JVS(490)*UV(53)+JVS(523)*UV(54)+JVS(546)&
               &*UV(55)
  JTUV(35) = JVS(134)*UV(30)+JVS(142)*UV(31)+JVS(172)*UV(34)+JVS(182)*UV(35)+JVS(247)*UV(42)+JVS(277)*UV(43)+JVS(411)&
               &*UV(51)+JVS(443)*UV(52)+JVS(547)*UV(55)+JVS(605)*UV(58)+JVS(636)*UV(59)
  JTUV(36) = JVS(135)*UV(30)+JVS(173)*UV(34)+JVS(189)*UV(36)+JVS(248)*UV(42)+JVS(278)*UV(43)+JVS(412)*UV(51)+JVS(444)&
               &*UV(52)+JVS(548)*UV(55)+JVS(606)*UV(58)+JVS(637)*UV(59)
  JTUV(37) = JVS(14)*UV(3)+JVS(197)*UV(37)+JVS(249)*UV(42)+JVS(279)*UV(43)+JVS(348)*UV(48)+JVS(413)*UV(51)+JVS(445)&
               &*UV(52)+JVS(549)*UV(55)+JVS(607)*UV(58)+JVS(638)*UV(59)
  JTUV(38) = JVS(15)*UV(3)+JVS(209)*UV(38)+JVS(280)*UV(43)+JVS(331)*UV(47)+JVS(349)*UV(48)+JVS(414)*UV(51)+JVS(446)&
               &*UV(52)+JVS(550)*UV(55)+JVS(608)*UV(58)+JVS(639)*UV(59)
  JTUV(39) = JVS(221)*UV(39)+JVS(250)*UV(42)+JVS(281)*UV(43)+JVS(415)*UV(51)+JVS(447)*UV(52)+JVS(551)*UV(55)+JVS(609)&
               &*UV(58)+JVS(640)*UV(59)
  JTUV(40) = JVS(16)*UV(3)+JVS(231)*UV(40)+JVS(251)*UV(42)+JVS(282)*UV(43)+JVS(350)*UV(48)+JVS(416)*UV(51)+JVS(448)&
               &*UV(52)+JVS(552)*UV(55)+JVS(610)*UV(58)+JVS(641)*UV(59)
  JTUV(41) = JVS(17)*UV(3)+JVS(160)*UV(33)+JVS(240)*UV(41)+JVS(252)*UV(42)+JVS(283)*UV(43)+JVS(332)*UV(47)+JVS(351)&
               &*UV(48)+JVS(417)*UV(51)+JVS(449)*UV(52)+JVS(553)*UV(55)+JVS(611)*UV(58)+JVS(642)*UV(59)
  JTUV(42) = JVS(253)*UV(42)+JVS(352)*UV(48)+JVS(418)*UV(51)+JVS(494)*UV(53)+JVS(571)*UV(56)
  JTUV(43) = JVS(99)*UV(25)+JVS(112)*UV(26)+JVS(284)*UV(43)+JVS(419)*UV(51)+JVS(495)*UV(53)+JVS(526)*UV(54)
  JTUV(44) = JVS(18)*UV(3)+JVS(161)*UV(33)+JVS(254)*UV(42)+JVS(285)*UV(43)+JVS(300)*UV(44)+JVS(333)*UV(47)+JVS(353)&
               &*UV(48)+JVS(420)*UV(51)+JVS(450)*UV(52)+JVS(554)*UV(55)+JVS(612)*UV(58)+JVS(644)*UV(59)
  JTUV(45) = JVS(19)*UV(3)+JVS(45)*UV(12)+JVS(113)*UV(26)+JVS(143)*UV(31)+JVS(174)*UV(34)+JVS(222)*UV(39)+JVS(255)&
               &*UV(42)+JVS(286)*UV(43)+JVS(311)*UV(45)+JVS(421)*UV(51)+JVS(451)*UV(52)+JVS(555)*UV(55)+JVS(613)*UV(58)&
               &+JVS(645)*UV(59)
  JTUV(46) = JVS(9)*UV(2)+JVS(20)*UV(3)+JVS(35)*UV(7)+JVS(61)*UV(17)+JVS(69)*UV(19)+JVS(90)*UV(24)+JVS(114)*UV(26)&
               &+JVS(122)*UV(27)+JVS(126)*UV(28)+JVS(130)*UV(29)+JVS(162)*UV(33)+JVS(287)*UV(43)+JVS(324)*UV(46)+JVS(355)&
               &*UV(48)+JVS(422)*UV(51)+JVS(452)*UV(52)+JVS(498)*UV(53)+JVS(529)*UV(54)+JVS(614)*UV(58)+JVS(646)*UV(59)
  JTUV(47) = JVS(163)*UV(33)+JVS(335)*UV(47)+JVS(356)*UV(48)+JVS(423)*UV(51)+JVS(499)*UV(53)+JVS(575)*UV(56)+JVS(647)&
               &*UV(59)
  JTUV(48) = JVS(100)*UV(25)+JVS(115)*UV(26)+JVS(357)*UV(48)+JVS(424)*UV(51)+JVS(453)*UV(52)+JVS(500)*UV(53)+JVS(530)&
               &*UV(54)+JVS(558)*UV(55)
  JTUV(49) = JVS(21)*UV(3)+JVS(164)*UV(33)+JVS(257)*UV(42)+JVS(288)*UV(43)+JVS(358)*UV(48)+JVS(371)*UV(49)+JVS(425)&
               &*UV(51)+JVS(454)*UV(52)+JVS(559)*UV(55)+JVS(615)*UV(58)+JVS(649)*UV(59)
  JTUV(50) = JVS(22)*UV(3)+JVS(175)*UV(34)+JVS(258)*UV(42)+JVS(289)*UV(43)+JVS(385)*UV(50)+JVS(426)*UV(51)+JVS(455)&
               &*UV(52)+JVS(560)*UV(55)+JVS(616)*UV(58)+JVS(650)*UV(59)
  JTUV(51) = JVS(51)*UV(14)+JVS(64)*UV(18)+JVS(73)*UV(20)+JVS(77)*UV(21)+JVS(101)*UV(25)+JVS(151)*UV(32)+JVS(183)*UV(35)&
               &+JVS(190)*UV(36)+JVS(199)*UV(37)+JVS(211)*UV(38)+JVS(223)*UV(39)+JVS(233)*UV(40)+JVS(241)*UV(41)+JVS(259)&
               &*UV(42)+JVS(301)*UV(44)+JVS(312)*UV(45)+JVS(325)*UV(46)+JVS(336)*UV(47)+JVS(372)*UV(49)+JVS(386)*UV(50)&
               &+JVS(427)*UV(51)+JVS(456)*UV(52)+JVS(503)*UV(53)+JVS(533)*UV(54)+JVS(561)*UV(55)+JVS(578)*UV(56)+JVS(594)&
               &*UV(57)+JVS(617)*UV(58)+JVS(651)*UV(59)
  JTUV(52) = JVS(23)*UV(3)+JVS(78)*UV(21)+JVS(116)*UV(26)+JVS(136)*UV(30)+JVS(144)*UV(31)+JVS(152)*UV(32)+JVS(165)&
               &*UV(33)+JVS(176)*UV(34)+JVS(184)*UV(35)+JVS(191)*UV(36)+JVS(200)*UV(37)+JVS(212)*UV(38)+JVS(224)*UV(39)&
               &+JVS(234)*UV(40)+JVS(242)*UV(41)+JVS(291)*UV(43)+JVS(302)*UV(44)+JVS(313)*UV(45)+JVS(361)*UV(48)+JVS(373)&
               &*UV(49)+JVS(387)*UV(50)+JVS(428)*UV(51)+JVS(457)*UV(52)+JVS(562)*UV(55)+JVS(579)*UV(56)+JVS(595)*UV(57)&
               &+JVS(618)*UV(58)+JVS(652)*UV(59)
  JTUV(53) = JVS(3)*UV(1)+JVS(29)*UV(4)+JVS(31)*UV(5)+JVS(33)*UV(6)+JVS(37)*UV(8)+JVS(39)*UV(9)+JVS(41)*UV(10)+JVS(43)&
               &*UV(11)+JVS(48)*UV(13)+JVS(52)*UV(14)+JVS(57)*UV(16)+JVS(62)*UV(17)+JVS(65)*UV(18)+JVS(74)*UV(20)+JVS(79)&
               &*UV(21)+JVS(83)*UV(22)+JVS(86)*UV(23)+JVS(91)*UV(24)+JVS(102)*UV(25)+JVS(117)*UV(26)+JVS(123)*UV(27)&
               &+JVS(127)*UV(28)+JVS(131)*UV(29)+JVS(137)*UV(30)+JVS(145)*UV(31)+JVS(166)*UV(33)+JVS(177)*UV(34)+JVS(185)&
               &*UV(35)+JVS(192)*UV(36)+JVS(201)*UV(37)+JVS(225)*UV(39)+JVS(235)*UV(40)+JVS(243)*UV(41)+JVS(261)*UV(42)&
               &+JVS(292)*UV(43)+JVS(303)*UV(44)+JVS(314)*UV(45)+JVS(326)*UV(46)+JVS(338)*UV(47)+JVS(362)*UV(48)+JVS(374)&
               &*UV(49)+JVS(388)*UV(50)+JVS(429)*UV(51)+JVS(458)*UV(52)+JVS(505)*UV(53)+JVS(535)*UV(54)+JVS(563)*UV(55)&
               &+JVS(580)*UV(56)+JVS(596)*UV(57)+JVS(619)*UV(58)+JVS(653)*UV(59)
  JTUV(54) = JVS(54)*UV(15)+JVS(70)*UV(19)+JVS(84)*UV(22)+JVS(92)*UV(24)+JVS(103)*UV(25)+JVS(118)*UV(26)+JVS(124)*UV(27)&
               &+JVS(128)*UV(28)+JVS(132)*UV(29)+JVS(138)*UV(30)+JVS(146)*UV(31)+JVS(154)*UV(32)+JVS(178)*UV(34)+JVS(214)&
               &*UV(38)+JVS(293)*UV(43)+JVS(315)*UV(45)+JVS(363)*UV(48)+JVS(430)*UV(51)+JVS(536)*UV(54)+JVS(564)*UV(55)&
               &+JVS(620)*UV(58)+JVS(654)*UV(59)
  JTUV(55) = JVS(24)*UV(3)+JVS(66)*UV(18)+JVS(87)*UV(23)+JVS(119)*UV(26)+JVS(139)*UV(30)+JVS(147)*UV(31)+JVS(155)*UV(32)&
               &+JVS(168)*UV(33)+JVS(179)*UV(34)+JVS(186)*UV(35)+JVS(193)*UV(36)+JVS(203)*UV(37)+JVS(215)*UV(38)+JVS(227)&
               &*UV(39)+JVS(237)*UV(40)+JVS(244)*UV(41)+JVS(294)*UV(43)+JVS(304)*UV(44)+JVS(316)*UV(45)+JVS(364)*UV(48)&
               &+JVS(376)*UV(49)+JVS(390)*UV(50)+JVS(431)*UV(51)+JVS(460)*UV(52)+JVS(565)*UV(55)+JVS(582)*UV(56)+JVS(598)&
               &*UV(57)+JVS(621)*UV(58)+JVS(655)*UV(59)
  JTUV(56) = JVS(25)*UV(3)+JVS(169)*UV(33)+JVS(264)*UV(42)+JVS(295)*UV(43)+JVS(341)*UV(47)+JVS(365)*UV(48)+JVS(432)&
               &*UV(51)+JVS(461)*UV(52)+JVS(566)*UV(55)+JVS(583)*UV(56)+JVS(622)*UV(58)+JVS(656)*UV(59)
  JTUV(57) = JVS(26)*UV(3)+JVS(265)*UV(42)+JVS(296)*UV(43)+JVS(366)*UV(48)+JVS(433)*UV(51)+JVS(462)*UV(52)+JVS(567)&
               &*UV(55)+JVS(600)*UV(57)+JVS(623)*UV(58)+JVS(657)*UV(59)
  JTUV(58) = JVS(49)*UV(13)+JVS(120)*UV(26)+JVS(140)*UV(30)+JVS(148)*UV(31)+JVS(170)*UV(33)+JVS(180)*UV(34)+JVS(187)&
               &*UV(35)+JVS(194)*UV(36)+JVS(204)*UV(37)+JVS(216)*UV(38)+JVS(228)*UV(39)+JVS(238)*UV(40)+JVS(245)*UV(41)&
               &+JVS(297)*UV(43)+JVS(305)*UV(44)+JVS(317)*UV(45)+JVS(328)*UV(46)+JVS(342)*UV(47)+JVS(367)*UV(48)+JVS(377)&
               &*UV(49)+JVS(392)*UV(50)+JVS(434)*UV(51)+JVS(463)*UV(52)+JVS(510)*UV(53)+JVS(540)*UV(54)+JVS(568)*UV(55)&
               &+JVS(585)*UV(56)+JVS(601)*UV(57)+JVS(624)*UV(58)+JVS(658)*UV(59)
  JTUV(59) = JVS(46)*UV(12)+JVS(55)*UV(15)+JVS(71)*UV(19)+JVS(75)*UV(20)+JVS(88)*UV(23)+JVS(104)*UV(25)+JVS(156)*UV(32)&
               &+JVS(318)*UV(45)+JVS(329)*UV(46)+JVS(343)*UV(47)+JVS(435)*UV(51)+JVS(511)*UV(53)+JVS(541)*UV(54)+JVS(569)&
               &*UV(55)+JVS(625)*UV(58)+JVS(659)*UV(59)
      
END SUBROUTINE JacTR_SP_Vec

! End of JacTR_SP_Vec function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



END MODULE radm2_Jacobian

