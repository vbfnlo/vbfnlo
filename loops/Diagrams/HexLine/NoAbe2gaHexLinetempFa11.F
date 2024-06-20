c************** Calling the Fa functions*************************
       subroutine NoAbe2gaHexLineFa11(p1mup2,p1mup3,p1mup5,p1mup6,p2mu
     -   p2,p2mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p
     -   4mup3,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup
     -   3,p6mup5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6
     -   ,mup5mup6,Fa11)
       IMPLICIT NONE
      Complex*16   p1mup2, p1mup3, p1mup5, p1mup6, p2mup2, p2mup3, 
     -          p2mup5, p2mup6, p3mup2, p3mup3, p3mup5, p3mup6, 
     -          p4mup2, p4mup3, p4mup5, p4mup6, p5mup2, p5mup3, 
     -          p5mup5, p5mup6, p6mup2, p6mup3, p6mup5, p6mup6
       Complex*16   mup2mup3, mup2mup5, mup2mup6, mup3mup5, mup3mup6, 
     -          mup5mup6
        common/NoAbe2gaHexLineFhlFunctions/F
       COMMON/NoAbe2gaHexLineInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,
     -   s12,s23,s34,s45,s56,s16,s123,s234,s345
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
      Complex*16 Fa(872),F(9706)
      Real*8 P(9687) 
        COMMON/NoAbe2gaHexLinePFunctions/P
       Complex*16 Fa11(241:264)
       COMMON/NoAbe2gaHexLineFaFunctions/Fa
       Fa11(241) = (p1mup2+p2mup2)*(p4mup6*F(678)+p6mup6*F(679)+(p1mup
     -   6+p2mup6+p3mup6)*F(680)+p5mup6*F(682))
       Fa11(242) = p1mup5+p2mup5+p3mup5+p4mup5+p5mup5+p6mup5
       Fa11(243) = 2*p6mup2*(2*((p1mup5+p2mup5+p3mup5)*(p1mup6+p2mup6+
     -   p3mup6)*F(632)+p2mup5*p4mup6*F(637)+p3mup5*p4mup6*F(637)-(p1
     -   mup6+p2mup6+p3mup6)*p4mup5*F(638)+p4mup5*p5mup6*F(639)+p2mup
     -   5*p5mup6*F(641)+p3mup5*p5mup6*F(641)-p4mup5*p4mup6*F(642)-p4
     -   mup6*p6mup5*F(643)-p1mup6*p6mup5*F(644)-p2mup6*p6mup5*F(644)
     -   -p3mup6*p6mup5*F(644)+p5mup6*p6mup5*F(645)-p2mup5*p6mup6*F(6
     -   46)-p3mup5*p6mup6*F(646)+p1mup5*(p4mup6*F(637)+p5mup6*F(641)
     -   -p6mup6*F(646))+p5mup5*(-(p5mup6*F(631))+(p1mup6+p2mup6+p3mu
     -   p6)*F(633)-p4mup6*F(636)+p6mup6*F(648))+p4mup5*p6mup6*F(649)
     -   +2*p6mup5*p6mup6*F(650))-mup5mup6*F(685))
       Fa11(244) = p5mup2*(4*((p1mup5+p2mup5+p3mup5)*(p1mup6+p2mup6+p3
     -   mup6)*F(632)+p2mup5*p4mup6*F(637)+p3mup5*p4mup6*F(637)-(p1mu
     -   p6+p2mup6+p3mup6)*p4mup5*F(638)+p4mup5*p5mup6*F(639)+p2mup5*
     -   p5mup6*F(641)+p3mup5*p5mup6*F(641)-p4mup5*p4mup6*F(642)-p4mu
     -   p6*p6mup5*F(643)-p1mup6*p6mup5*F(644)-p2mup6*p6mup5*F(644)-p
     -   3mup6*p6mup5*F(644)+p5mup6*p6mup5*F(645)-p2mup5*p6mup6*F(646
     -   )-p3mup5*p6mup6*F(646)+p1mup5*(p4mup6*F(637)+p5mup6*F(641)-p
     -   6mup6*F(646))+p5mup5*(-(p5mup6*F(631))+(p1mup6+p2mup6+p3mup6
     -   )*F(633)-p4mup6*F(636)+p6mup6*F(648))+p4mup5*p6mup6*F(649)+2
     -   *p6mup5*p6mup6*F(650))-mup5mup6*F(684))
       Fa11(245) = 2*(p1mup5*p5mup2*F(653)+p2mup5*p5mup2*F(653)+p3mup5
     -   *p5mup2*F(653)+p5mup2*p5mup5*F(655)+p5mup5*p6mup2*F(656)+p6m
     -   up2*p6mup5*F(658))+p5mup2*p6mup5*F(705)+p1mup5*p6mup2*F(706)
     -   +p2mup5*p6mup2*F(706)+p3mup5*p6mup2*F(706)-p4mup5*p5mup2*F(7
     -   07)+p4mup5*p6mup2*F(708)
       Fa11(246) = -((p1mup2+p2mup2)*(2*p5mup5*F(655)-(p1mup5+p2mup5+p
     -   3mup5)*F(686)-p6mup5*F(688)-p4mup5*F(690)))
       Fa11(247) = -((p1mup2+p2mup2)*(-((p1mup5+p2mup5+p3mup5)*F(687))
     -   +2*(p5mup5*F(656)-p6mup5*F(689))-p4mup5*F(691)))
       Fa11(248) = -2*(p1mup2+p2mup2)*(-2*((p1mup5+p2mup5+p3mup5)*(p1m
     -   up6+p2mup6+p3mup6)*F(632)+p2mup5*p4mup6*F(637)+p3mup5*p4mup6
     -   *F(637)-(p1mup6+p2mup6+p3mup6)*p4mup5*F(638)+p4mup5*p5mup6*F
     -   (639)+p2mup5*p5mup6*F(641)+p3mup5*p5mup6*F(641)-p4mup5*p4mup
     -   6*F(642)-p4mup6*p6mup5*F(643)-p1mup6*p6mup5*F(644)-p2mup6*p6
     -   mup5*F(644)-p3mup6*p6mup5*F(644)+p5mup6*p6mup5*F(645)-p2mup5
     -   *p6mup6*F(646)-p3mup5*p6mup6*F(646)+p1mup5*(p4mup6*F(637)+p5
     -   mup6*F(641)-p6mup6*F(646))+p5mup5*(-(p5mup6*F(631))+(p1mup6+
     -   p2mup6+p3mup6)*F(633)-p4mup6*F(636)+p6mup6*F(648))+p4mup5*p6
     -   mup6*F(649)+2*p6mup5*p6mup6*F(650))+mup5mup6*F(685))
       Fa11(249) = 2*(p2mup5*p4mup6*F(637)+p3mup5*p4mup6*F(637)+p4mup5
     -   *p5mup6*F(639)-p2mup5*p6mup6*F(646)-p3mup5*p6mup6*F(646)+p5m
     -   up5*(-(p5mup6*F(631))+(p1mup6+p2mup6+p3mup6)*F(633)-p4mup6*F
     -   (636)+p6mup6*F(648))+p4mup5*p6mup6*F(649))+4*p6mup5*p6mup6*F
     -   (650)-mup5mup6*F(709)-p1mup6*p4mup5*F(710)-p2mup6*p4mup5*F(7
     -   10)-p3mup6*p4mup5*F(710)-p4mup6*p6mup5*F(711)+p2mup5*p5mup6*
     -   F(712)+p3mup5*p5mup6*F(712)-p1mup6*p6mup5*F(713)-p2mup6*p6mu
     -   p5*F(713)-p3mup6*p6mup5*F(713)-p5mup6*p6mup5*F(714)+p1mup6*p
     -   2mup5*F(715)+p2mup5*p2mup6*F(715)+p1mup6*p3mup5*F(715)+p2mup
     -   6*p3mup5*F(715)+p2mup5*p3mup6*F(715)+p3mup5*p3mup6*F(715)+p1
     -   mup5*(2*(p4mup6*F(637)-p6mup6*F(646))+p5mup6*F(712)+(p1mup6+
     -   p2mup6+p3mup6)*F(715))-p4mup5*p4mup6*F(716)
       Fa11(250) = (p1mup5+p2mup5+p3mup5)*F(717)+p5mup5*F(718)+p6mup5*
     -   F(719)+p4mup5*F(720)
       Fa11(251) = (p1mup5+p2mup5+p3mup5)*F(721)+p5mup5*F(722)+p6mup5*
     -   F(723)+p4mup5*F(724)
       Fa11(252) = (p1mup6+p2mup6+p3mup6)*F(725)+p6mup6*F(726)+p4mup6*
     -   F(727)+p5mup6*F(728)
       Fa11(253) = (p1mup6+p2mup6+p3mup6)*F(729)+p6mup6*F(730)+p4mup6*
     -   F(731)+p5mup6*F(732)
       Fa11(254) = mup3mup6*(p1mup5*F(733)+p2mup5*F(735)+p5mup5*F(737)
     -   +p6mup5*F(739)+(p3mup5+p4mup5)*F(741))+mup3mup5*(p1mup6*F(73
     -   4)+p2mup6*F(736)+p5mup6*F(738)+2*p6mup6*F(740)+(p3mup6+p4mup
     -   6)*F(742))
       Fa11(255) = mup3mup6*(p2mup2*F(743)+4*(p5mup2*F(744)+p1mup2*F(7
     -   45)+p6mup2*F(746)+(p3mup2+p4mup2)*F(747)))
       Fa11(256) = mup5mup6*(p2mup2*F(748)+p5mup2*F(749)+2*(p3mup2+p4m
     -   up2)*F(750))
       Fa11(257) = p2mup2*p2mup6*(4*p5mup5*F(753)-p6mup5*F(754)+p2mup5
     -   *F(756))+p5mup2*p5mup5*p5mup6*F(765)+p1mup5*(4*p5mup2*p5mup6
     -   *F(761)+p2mup2*(p1mup6*F(755)-4*(p2mup6*F(752)-p5mup6*F(764)
     -   ))+(p3mup2+p4mup2)*(p1mup6*F(757)+4*(p2mup6*F(766)+p5mup6*F(
     -   767)+p6mup6*F(768))))+4*(p1mup6*p2mup5*p3mup2*F(758)+p1mup6*
     -   p2mup5*p4mup2*F(758)+p1mup6*p3mup2*p5mup5*F(759)+p1mup6*p4mu
     -   p2*p5mup5*F(759)+p1mup6*p3mup2*p6mup5*F(760)+p1mup6*p4mup2*p
     -   6mup5*F(760)+p2mup5*p5mup2*p5mup6*F(762)-p5mup2*p5mup6*p6mup
     -   5*F(763)+p2mup5*p2mup6*p3mup2*F(769)+p2mup5*p2mup6*p4mup2*F(
     -   769)+p2mup5*p3mup2*p5mup6*F(770)+p2mup5*p4mup2*p5mup6*F(770)
     -   +p2mup5*p3mup2*p6mup6*F(771)+p2mup5*p4mup2*p6mup6*F(771)+p2m
     -   up6*p3mup2*p5mup5*F(772)+p2mup6*p4mup2*p5mup5*F(772)-p2mup6*
     -   p3mup2*p6mup5*F(773)-p2mup6*p4mup2*p6mup5*F(773)+p3mup2*p5mu
     -   p5*p5mup6*F(774)+p4mup2*p5mup5*p5mup6*F(774)+p3mup2*p5mup5*p
     -   6mup6*F(775)+p4mup2*p5mup5*p6mup6*F(775)+p3mup2*p5mup6*p6mup
     -   5*F(776)+p4mup2*p5mup6*p6mup5*F(776)-(p3mup2+p4mup2)*p6mup5*
     -   p6mup6*F(777))
       Fa11(258) = mup5mup6*p6mup2*F(781)+mup2mup6*(p1mup5*F(778)+p2mu
     -   p5*F(779)+p5mup5*F(780)+p6mup5*F(782))+p1mup6*p2mup2*p5mup5*
     -   F(785)+p2mup2*p5mup6*p6mup5*F(789)+p1mup6*p2mup5*p5mup2*F(79
     -   4)+p2mup5*p5mup2*p6mup6*F(797)+p5mup2*p5mup5*p6mup6*F(798)+p
     -   2mup5*p2mup6*p5mup2*F(802)+p1mup5*(-4*p2mup2*p6mup6*F(783)+p
     -   5mup2*(-4*p6mup6*F(793)+p2mup6*F(801)+p1mup6*F(803)))+4*(p1m
     -   up6*p2mup2*p2mup5*F(784)-p1mup6*p2mup2*p6mup5*F(786)+p2mup2*
     -   p2mup5*p5mup6*F(787)-p2mup2*p5mup5*p6mup6*F(788)-p2mup2*p6mu
     -   p5*p6mup6*F(790)+p2mup6*p5mup2*p5mup5*F(791)-p2mup6*p5mup2*p
     -   6mup5*F(792)-p1mup6*p5mup2*p5mup5*F(795)+p1mup6*p5mup2*p6mup
     -   5*F(796)+p2mup2*p2mup5*p6mup6*F(799)+p2mup2*p5mup5*p5mup6*F(
     -   800)-p5mup2*p6mup5*p6mup6*F(804))
       Fa11(259) = mup2mup6*p4mup5*F(820)+4*p2mup2*p4mup5*(p5mup6*F(82
     -   1)+p6mup6*F(823)+4*p4mup6*F(826))+p1mup2*(p1mup6*(p5mup5*F(8
     -   06)-4*(p6mup5*F(807)-p2mup5*F(813)))+4*(p5mup5*p5mup6*F(809)
     -   -p5mup5*p6mup6*F(810)-p5mup6*p6mup5*F(811)-p6mup5*p6mup6*F(8
     -   12)+p2mup5*p6mup6*F(815)-p2mup6*p5mup5*F(816)-p2mup6*p6mup5*
     -   F(817)+p1mup5*(p1mup6*F(805)-p6mup6*F(808)+p2mup6*F(814)+p5m
     -   up6*F(818))+p2mup5*p5mup6*F(819)+p4mup5*p6mup6*F(822)-p2mup5
     -   *p2mup6*F(824)+4*p4mup5*p4mup6*F(825)+p4mup5*p5mup6*F(827)))
       Fa11(260) = mup2mup6*p3mup5*F(820)+mup2mup5*p6mup6*F(829)+p1mup
     -   6*p2mup5*p6mup2*F(830)+p2mup5*p5mup6*p6mup2*F(831)+p2mup5*p6
     -   mup2*p6mup6*F(832)+p1mup6*p3mup5*p5mup2*F(834)+p1mup6*p4mup5
     -   *p5mup2*F(834)+4*(p1mup2*p3mup5*p6mup6*F(822)+p1mup2*p3mup5*
     -   p5mup6*F(827)+mup2mup5*p1mup6*F(828)-mup2mup5*p2mup6*F(833)+
     -   p3mup5*p5mup2*p6mup6*F(835)+p4mup5*p5mup2*p6mup6*F(835)+p2mu
     -   p2*(p3mup5*p5mup6*F(821)+p3mup5*p6mup6*F(823)+p1mup6*p3mup5*
     -   F(836)+p1mup6*p4mup5*F(836)+p2mup6*(p3mup5+p4mup5)*F(837))+p
     -   1mup2*p1mup6*p3mup5*F(838)+p1mup2*p1mup6*p4mup5*F(838)-mup2m
     -   up5*p5mup6*F(839)+p2mup6*p3mup5*p5mup2*F(840)+p2mup6*p4mup5*
     -   p5mup2*F(840)+p3mup5*p5mup2*p5mup6*F(841)+p4mup5*p5mup2*p5mu
     -   p6*F(841)+p1mup2*p2mup6*p3mup5*F(842)+p1mup2*p2mup6*p4mup5*F
     -   (842)+p3mup2*p3mup5*p5mup6*F(843)+p3mup5*p4mup2*p5mup6*F(843
     -   )+p3mup2*p4mup5*p5mup6*F(843)+p4mup2*p4mup5*p5mup6*F(843)+p2
     -   mup6*p3mup2*p3mup5*F(844)+p2mup6*p3mup5*p4mup2*F(844)+p2mup6
     -   *p3mup2*p4mup5*F(844)+p2mup6*p4mup2*p4mup5*F(844)+p1mup6*p3m
     -   up2*p3mup5*F(845)+p1mup6*p3mup5*p4mup2*F(845)+p1mup6*p3mup2*
     -   p4mup5*F(845)+p1mup6*p4mup2*p4mup5*F(845)+(p3mup2+p4mup2)*(p
     -   3mup5+p4mup5)*p6mup6*F(846))
       Fa11(261) = p4mup6*(p2mup2*(16*p3mup5*F(826)+4*p5mup5*F(848)+p1
     -   mup5*F(850)+p2mup5*F(851))+4*p1mup2*(4*p3mup5*F(825)+p5mup5*
     -   F(847)+p1mup5*F(849)+p2mup5*F(852)))
       Fa11(262) = p6mup2*(p5mup5*(p5mup6*F(853)+p2mup6*F(854))+p2mup5
     -   *p2mup6*F(855)+4*(p3mup5+p4mup5)*(p2mup6*F(856)-p5mup6*F(857
     -   )+p6mup6*F(858)-p1mup6*F(859)))
       Fa11(263) = p2mup2*p3mup6*(16*p4mup5*F(826)+4*p5mup5*F(848)+p1m
     -   up5*F(850)+p2mup5*F(851))+p1mup6*p5mup5*p6mup2*F(860)+p1mup5
     -   *p2mup6*p6mup2*F(862)+p1mup5*p5mup6*p6mup2*F(863)-p1mup5*p6m
     -   up2*p6mup6*F(864)+p1mup5*p1mup6*p6mup2*F(865)-4*(p5mup5*p6mu
     -   p2*p6mup6*F(861)+(p3mup6+p4mup6)*(mup2mup5*F(866)-p2mup5*p5m
     -   up2*F(867))-p3mup6*p5mup2*p5mup5*F(868)-p4mup6*p5mup2*p5mup5
     -   *F(868)-p2mup2*(p3mup6+p4mup6)*p6mup5*F(869)-p1mup2*(p3mup6*
     -   (4*p4mup5*F(825)+p5mup5*F(847)+p1mup5*F(849)+p2mup5*F(852))+
     -   (p3mup6+p4mup6)*p6mup5*F(870))-p1mup5*p3mup6*p5mup2*F(871)-p
     -   1mup5*p4mup6*p5mup2*F(871)+p3mup6*p5mup2*p6mup5*F(872)+p4mup
     -   6*p5mup2*p6mup5*F(872)-p3mup2*p3mup6*p5mup5*F(873)-p3mup6*p4
     -   mup2*p5mup5*F(873)-p3mup2*p4mup6*p5mup5*F(873)-p4mup2*p4mup6
     -   *p5mup5*F(873)-p2mup5*p3mup2*p3mup6*F(874)-p2mup5*p3mup6*p4m
     -   up2*F(874)-p2mup5*p3mup2*p4mup6*F(874)-p2mup5*p4mup2*p4mup6*
     -   F(874)-p1mup5*p3mup2*p3mup6*F(875)-p1mup5*p3mup6*p4mup2*F(87
     -   5)-p1mup5*p3mup2*p4mup6*F(875)-p1mup5*p4mup2*p4mup6*F(875)+(
     -   p3mup2+p4mup2)*(p3mup6+p4mup6)*p6mup5*F(876))
       Fa11(264) = 4*(p2mup6*p6mup2*p6mup5*F(877)+p1mup5*(p3mup6+p4mup
     -   6)*p6mup2*F(878)+p3mup6*p5mup5*p6mup2*F(879)+p4mup6*p5mup5*p
     -   6mup2*F(879)+p3mup6*p6mup2*p6mup5*F(880)+p4mup6*p6mup2*p6mup
     -   5*F(880)+p5mup6*p6mup2*p6mup5*F(881)+p1mup6*p6mup2*p6mup5*F(
     -   882))+p2mup5*p3mup6*p6mup2*F(883)+p2mup5*p4mup6*p6mup2*F(883
     -   )-16*(p3mup5*p3mup6*p5mup2*F(884)-p4mup2*(p3mup6*p4mup5+(p3m
     -   up5+p4mup5)*p4mup6)*F(885))
       Return
       End
