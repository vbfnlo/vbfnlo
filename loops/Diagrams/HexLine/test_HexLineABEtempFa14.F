c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineAbe14(p1mup2,p1mup3,p1mup4,p1m
     -   up5,p2mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,
     -   p4mup2,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mu
     -   p2,p6mup3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,
     -   mup3mup5,mup4mup5,Fa14)
       IMPLICIT NONE
      Complex*16   p1mup2, p1mup3, p1mup4, p1mup5, p2mup2, p2mup3, 
     -          p2mup4, p2mup5, p3mup2, p3mup3, p3mup4, p3mup5, 
     -          p4mup2, p4mup3, p4mup4, p4mup5, p5mup2, p5mup3, 
     -          p5mup4, p5mup5, p6mup2, p6mup3, p6mup4, p6mup5
       Complex*16   mup2mup3, mup2mup4, mup2mup5, mup3mup4, mup3mup5, 
     -          mup4mup5
        common/test_FhlFunctionsAbe/F
       COMMON/Invariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12,s23,s34,s45
     -   ,s56,s16,s123,s234,s345
       Real*8   p1sq 
       Real*8   p2sq 
       Real*8   p3sq 
       Real*8   p4sq 
       Real*8   p5sq 
       Real*8   p6sq 
       Real*8   s12, s16 
       Real*8   s23 
       Real*8   s34 
       Real*8   s45 
       Real*8   s56 
       Real*8   s123,s234,s345
      Complex*16 Fa(445),F(5930)
      Real*8 P(11674) 
        COMMON/test_PFunctionsAbe/P
       Complex*16 Fa14(144:154)
       COMMON/test_FaFunctionsAbe/Fa
       Fa14(144) = p2mup4*F(459)+p3mup4*F(460)+p4mup4*F(461)+p1mup4*F(
     -   462)
       Fa14(145) = 4*mup3mup4*F(464)+p2mup4*p3mup3*F(471)-p3mup3*p3mup
     -   4*F(472)-p3mup3*p4mup4*F(473)+p1mup4*p4mup3*F(474)+p2mup3*(p
     -   2mup4*F(468)+8*(p3mup4*F(469)+p4mup4*F(470)+p1mup4*F(476)))+
     -   8*(p2mup4*p4mup3*F(475)+p3mup4*p4mup3*F(477))+p4mup3*p4mup4*
     -   F(478)+p1mup3*(p2mup4*F(465)-8*(p3mup4*F(466)+p4mup4*F(467))
     -   +p1mup4*F(479))+p1mup4*p3mup3*F(480)
       Fa14(146) = 4*mup3mup4*F(481)+p1mup4*p2mup3*F(483)+p2mup4*p3mup
     -   3*F(488)+8*(p1mup4*p3mup3*F(484)-p2mup3*p2mup4*F(485)-p2mup3
     -   *(p3mup4*F(486)+p4mup4*F(487))-p3mup3*p3mup4*F(489))+p3mup3*
     -   p4mup4*F(490)+p4mup3*p4mup4*F(493)+p2mup4*p4mup3*F(494)+p3mu
     -   p4*p4mup3*F(495)+p1mup3*(p3mup4*F(482)+p2mup4*F(491)+p4mup4*
     -   F(492)+p1mup4*F(496))+p1mup4*p4mup3*F(497)
       Fa14(147) = mup2mup4*F(498)+p3mup2*p3mup4*F(501)+p1mup2*p3mup4*
     -   F(503)+p2mup2*p3mup4*F(504)+p3mup4*p4mup2*F(509)+p2mup2*p4mu
     -   p4*F(510)+p1mup4*(p3mup2*F(499)+p2mup2*F(506)+p4mup2*F(508)-
     -   8*p1mup2*F(512))+p2mup4*(p3mup2*F(500)+p1mup2*F(505)+p4mup2*
     -   F(511)-8*p2mup2*F(513))+8*(p3mup2*p4mup4*F(502)+p1mup2*p4mup
     -   4*F(507)-p4mup2*p4mup4*F(514))
       Fa14(148) = p3mup4*F(515)+p1mup4*F(516)+p2mup4*F(517)+p4mup4*F(
     -   518)
       Fa14(149) = mup2mup4*(p4mup3*F(519)+p1mup3*F(525)-p2mup3*F(526)
     -   -p3mup3*F(527))+mup3mup4*(p4mup2*F(522)+p2mup2*F(524)+p1mup2
     -   *F(528)+p3mup2*F(529))+mup2mup3*(p3mup4*F(521)+8*(p2mup4*F(5
     -   20)+p4mup4*F(523)+p1mup4*F(530)))-16*(p1mup4*p2mup2*p2mup3*F
     -   (547)+p1mup4*p2mup2*p3mup3*F(548)+p2mup2*p2mup3*p2mup4*F(549
     -   )+p2mup2*p2mup4*p3mup3*F(550)+p2mup2*p2mup3*p3mup4*F(550)+p2
     -   mup2*p2mup3*p4mup4*F(551)+p2mup2*p3mup3*p3mup4*F(552)+p2mup2
     -   *p3mup3*p4mup4*F(553)-p1mup4*p2mup3*p4mup2*F(557)+p1mup2*(p1
     -   mup3*(p1mup4*F(531)+p2mup4*F(532)+p3mup4*F(533)+p4mup4*F(534
     -   ))+p1mup4*p2mup3*F(535)+p1mup4*p3mup3*F(536)+p2mup3*p2mup4*F
     -   (537)+p2mup3*p3mup4*F(538)+p2mup3*p4mup4*F(539)+p2mup4*p3mup
     -   3*F(540)+p3mup3*p3mup4*F(541)+p3mup3*p4mup4*F(542)-p1mup4*p4
     -   mup3*F(554)-p2mup4*p4mup3*F(555)-p3mup4*p4mup3*F(558)-p4mup3
     -   *p4mup4*F(561))-p2mup2*p4mup3*p4mup4*F(562)-p1mup4*p2mup2*p4
     -   mup3*F(563)-p2mup2*p2mup4*p4mup3*F(564)-p2mup2*p3mup4*p4mup3
     -   *F(565)-p3mup2*p3mup4*p4mup3*F(566)-p2mup3*p4mup2*p4mup4*F(5
     -   67)-p3mup3*p4mup2*p4mup4*F(568)+p1mup4*p4mup2*p4mup3*F(570)-
     -   p2mup3*p2mup4*p4mup2*F(571)-p2mup3*p3mup4*p4mup2*F(572)-p2mu
     -   p4*p3mup2*p4mup3*F(573)-p1mup4*p3mup2*p4mup3*F(574)-p1mup4*p
     -   3mup2*p3mup3*F(576)-p2mup4*p3mup2*p3mup3*F(577)-p2mup3*p3mup
     -   2*p3mup4*F(577)-p3mup2*p3mup3*p3mup4*F(578)-p3mup2*p3mup3*p4
     -   mup4*F(579)-p3mup3*p3mup4*p4mup2*F(580)-p1mup4*p2mup3*p3mup2
     -   *F(582)-p2mup3*p2mup4*p3mup2*F(583)+p1mup4*p3mup3*p4mup2*F(5
     -   85)+p1mup3*(p1mup4*p2mup2*F(543)+p2mup2*p2mup4*F(544)+p2mup2
     -   *p3mup4*F(545)+p2mup2*p4mup4*F(546)+p2mup4*p4mup2*F(556)+p3m
     -   up4*p4mup2*F(559)+p1mup4*p4mup2*F(560)-p4mup2*p4mup4*F(569)-
     -   p3mup2*p3mup4*F(575)-p2mup4*p3mup2*F(581)-p1mup4*p3mup2*F(58
     -   4)-p3mup2*p4mup4*F(586))-p2mup4*p3mup3*p4mup2*F(587)+p2mup4*
     -   p4mup2*p4mup3*F(588)-p3mup4*p4mup2*p4mup3*F(589)-p2mup3*p3mu
     -   p2*p4mup4*F(590)+p4mup2*p4mup3*p4mup4*F(591)-p3mup2*p4mup3*p
     -   4mup4*F(592))
       Fa14(150) = mup2mup4*(p2mup3*F(593)+p3mup3*F(594)+p1mup3*F(596)
     -   +p4mup3*F(599))+mup3mup4*(p3mup2*F(595)+p1mup2*F(597)+p2mup2
     -   *F(600)+p4mup2*F(602))+mup2mup3*(p1mup4*F(598)+p2mup4*F(601)
     -   +p4mup4*F(603)+p3mup4*F(604))-16*(p1mup4*p2mup2*p2mup3*F(621
     -   )+p1mup4*p2mup2*p3mup3*F(622)+p2mup2*p2mup3*p2mup4*F(623)+p2
     -   mup2*p2mup3*p3mup4*F(624)+p2mup2*p2mup3*p4mup4*F(625)+p2mup2
     -   *p2mup4*p3mup3*F(626)+p2mup2*p3mup3*p3mup4*F(627)+p2mup2*p3m
     -   up3*p4mup4*F(628)+p1mup4*p2mup3*p3mup2*F(630)-p1mup4*p3mup2*
     -   p3mup3*F(634)-p2mup3*p3mup2*p3mup4*F(635)-p2mup4*p3mup2*p3mu
     -   p3*F(636)-p3mup2*p3mup3*p3mup4*F(637)-p3mup2*p3mup3*p4mup4*F
     -   (638)+p3mup2*p3mup4*p4mup3*F(639)-p3mup3*p3mup4*p4mup2*F(640
     -   )-p2mup3*p2mup4*p3mup2*F(641)-p2mup3*p3mup4*p4mup2*F(642)-p1
     -   mup4*p3mup2*p4mup3*F(644)-p2mup3*p4mup2*p4mup4*F(646)-p2mup2
     -   *p4mup3*p4mup4*F(646)+p1mup2*(p1mup3*(p1mup4*F(605)+p2mup4*F
     -   (606)+p3mup4*F(607)+p4mup4*F(608))+p1mup4*p2mup3*F(609)+p1mu
     -   p4*p3mup3*F(610)+p2mup3*p2mup4*F(611)+p2mup3*p3mup4*F(612)+p
     -   2mup3*p4mup4*F(613)+p2mup4*p3mup3*F(614)+p3mup3*p3mup4*F(615
     -   )+p3mup3*p4mup4*F(616)+p3mup4*p4mup3*F(633)-p4mup3*p4mup4*F(
     -   645)-p1mup4*p4mup3*F(647)-p2mup4*p4mup3*F(648))-p1mup4*p2mup
     -   2*p4mup3*F(649)-p2mup3*p2mup4*p4mup2*F(650)-p2mup2*p2mup4*p4
     -   mup3*F(650)-p2mup2*p3mup4*p4mup3*F(651)-p3mup3*p4mup2*p4mup4
     -   *F(652)-p1mup4*p2mup3*p4mup2*F(655)-p2mup4*p3mup2*p4mup3*F(6
     -   56)+p1mup3*(p1mup4*p2mup2*F(617)+p2mup2*p2mup4*F(618)+p2mup2
     -   *p3mup4*F(619)+p2mup2*p4mup4*F(620)-p2mup4*p3mup2*F(629)-p3m
     -   up2*p3mup4*F(631)-p1mup4*p3mup2*F(632)-p3mup4*p4mup2*F(643)-
     -   p4mup2*p4mup4*F(653)-p2mup4*p4mup2*F(654)+p3mup2*p4mup4*F(65
     -   7)-p1mup4*p4mup2*F(658))-p1mup4*p3mup3*p4mup2*F(659)+p2mup3*
     -   p3mup2*p4mup4*F(660)-p1mup4*p4mup2*p4mup3*F(661)+p2mup4*p4mu
     -   p2*p4mup3*F(662)+p3mup4*p4mup2*p4mup3*F(663)+p4mup2*p4mup3*p
     -   4mup4*F(664)-p2mup4*p3mup3*p4mup2*F(665)+p3mup2*p4mup3*p4mup
     -   4*F(666))
       Fa14(151) = -(mup3mup4*(p1mup2*F(668)+p2mup2*F(669)+p3mup2*F(67
     -   3)+p4mup2*F(675)))+mup2mup3*(-(p2mup4*F(672))+8*(p1mup4*F(67
     -   0)-p3mup4*F(674)+p4mup4*F(676)))-8*mup2mup4*(p1mup3*F(667)-p
     -   2mup3*F(671)+p4mup3*F(701)+p3mup3*F(733))+16*(-(p1mup4*p2mup
     -   2*p2mup3*F(680))+p1mup4*p2mup2*p3mup3*F(694)+p2mup2*p2mup3*p
     -   2mup4*F(695)+p2mup2*p2mup3*p3mup4*F(696)+p2mup2*p2mup3*p4mup
     -   4*F(697)+p2mup2*p2mup4*p3mup3*F(698)+p2mup2*p3mup3*p3mup4*F(
     -   699)+p2mup2*p3mup3*p4mup4*F(700)+p2mup2*p4mup3*p4mup4*F(703)
     -   +p1mup2*(p1mup3*(p2mup4*F(677)+p1mup4*F(681)+p3mup4*F(682)+p
     -   4mup4*F(683))+p2mup3*p2mup4*F(685)+p2mup3*p3mup4*F(686)+p2mu
     -   p3*p4mup4*F(687)+p2mup4*p3mup3*F(688)+p3mup3*p3mup4*F(689)+p
     -   3mup3*p4mup4*F(690)+p4mup3*p4mup4*F(702)+p1mup4*(-(p2mup3*F(
     -   678))+p3mup3*F(684)+p4mup3*F(704))+p2mup4*p4mup3*F(705)+p3mu
     -   p4*p4mup3*F(706))+p1mup4*p2mup2*p4mup3*F(707)+p2mup2*p2mup4*
     -   p4mup3*F(708)-p1mup4*p3mup2*p3mup3*F(710)-p2mup3*p3mup2*p3mu
     -   p4*F(711)-p2mup4*p3mup2*p3mup3*F(712)-p3mup2*p3mup3*p3mup4*F
     -   (713)-p3mup2*p3mup3*p4mup4*F(714)+p3mup3*p3mup4*p4mup2*F(715
     -   )+p2mup3*p4mup2*p4mup4*F(716)+p3mup3*p4mup2*p4mup4*F(717)-p1
     -   mup4*p2mup3*p3mup2*F(721)+p1mup4*p2mup3*p4mup2*F(722)-p2mup3
     -   *p2mup4*p3mup2*F(723)+p2mup3*p2mup4*p4mup2*F(724)-p2mup3*p3m
     -   up2*p4mup4*F(725)+p2mup3*p3mup4*p4mup2*F(726)+p2mup4*p3mup3*
     -   p4mup2*F(727)+p1mup4*p3mup3*p4mup2*F(731)+p1mup3*(p2mup2*(p2
     -   mup4*F(679)+p1mup4*F(691)+p3mup4*F(692)+p4mup4*F(693))-p3mup
     -   2*(p3mup4*F(709)+p2mup4*F(719)+p1mup4*F(728)+p4mup4*F(729))+
     -   p4mup2*(p4mup4*F(718)+p2mup4*F(720)+p3mup4*F(730)+p1mup4*F(7
     -   32)))+p3mup2*p3mup4*p4mup3*F(734)+p3mup2*p4mup3*p4mup4*F(735
     -   )+p2mup4*p3mup2*p4mup3*F(736)+p1mup4*p3mup2*p4mup3*F(737)-p1
     -   mup4*p4mup2*p4mup3*F(738)-p2mup2*p3mup4*p4mup3*F(739)-p2mup4
     -   *p4mup2*p4mup3*F(740)+p3mup4*p4mup2*p4mup3*F(741)+p4mup2*p4m
     -   up3*p4mup4*F(742))
       Fa14(152) = mup3mup4*F(743)+8*(-(p1mup4*p2mup3*F(745))+p1mup3*(
     -   p2mup4*F(744)+p1mup4*F(746)+p3mup4*F(747)+p4mup4*F(748))+p1m
     -   up4*p3mup3*F(749)+p2mup3*p2mup4*F(750)+p2mup3*p3mup4*F(751)+
     -   p2mup3*p4mup4*F(752)+p2mup4*p3mup3*F(753)+p3mup3*p3mup4*F(75
     -   4)+p3mup3*p4mup4*F(755)+p4mup3*p4mup4*F(756)+p1mup4*p4mup3*F
     -   (757)+p2mup4*p4mup3*F(758)-p3mup4*p4mup3*F(759))
       Fa14(153) = p2mup3*F(761)-4*(p1mup3*F(760)+p3mup3*F(762)+p4mup3
     -   *F(763))
       Fa14(154) = mup2mup3*F(764)+8*(p1mup3*p2mup2*F(766)+p2mup2*p2mu
     -   p3*F(768)+p2mup2*p3mup3*F(770)+p1mup2*(p1mup3*F(765)+p2mup3*
     -   F(767)+p3mup3*F(769)+p4mup3*F(771))+p2mup2*p4mup3*F(772)+p3m
     -   up2*p4mup3*F(773)-p1mup3*p4mup2*F(774)-p2mup3*p4mup2*F(775)-
     -   p3mup3*p4mup2*F(776)-p4mup2*p4mup3*F(777)+p1mup3*p3mup2*F(77
     -   8)+p2mup3*p3mup2*F(779)-p3mup2*p3mup3*F(780))
       Return
       End
