c************** Calling the Fa functions*************************
       subroutine HexBoxCross1Fa14(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,
     &   p3mup4,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa
     &   14)
       IMPLICIT NONE
      Complex*16   p1mup3, p1mup4, p2mup3, 
     -          p2mup4, p3mup3, p3mup4, 
     -          p4mup3, p4mup4, p5mup3, 
     -          p5mup4, p6mup3, p6mup4
       Complex*16   mup3mup4

        common/HexBoxCross1FhlFunctions/F
       COMMON/HexBoxCross1Invariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12
     &   ,s23,s34,s45,s56,s16,s123,s234,s345
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
      Complex*16 Fa(272),F(1461)
      Real*8 P(1825) 
        COMMON/HexBoxCross1PFunctions/P
       Complex*16 Fa14(170:182)
       COMMON/HexBoxCross1FaFunctions/Fa
       Fa14(170) = mup3mup4*F(686)-2*(-(p1mup4*p6mup3*F(690))-p1mup4*p
     &   2mup3*F(692)+p5mup4*p6mup3*F(695)-p2mup3*p4mup4*F(696)-p2mup
     &   3*p6mup4*F(697)-p2mup3*p5mup4*F(698)+p6mup3*p6mup4*F(699)-p2
     &   mup3*p2mup4*F(701)-p2mup4*p6mup3*F(703)+p5mup3*(p5mup4*F(637
     &   )-p1mup4*F(689)+p6mup4*F(693)-p2mup4*F(702)+p4mup4*F(705))+p
     &   4mup4*p6mup3*F(706)+p4mup3*(p5mup4*F(636)+p6mup4*F(694)+p1mu
     &   p4*F(700)+p4mup4*F(704)-p2mup4*F(708))+p1mup3*(p5mup4*F(687)
     &   +p6mup4*F(688)+p2mup4*F(691)-p1mup4*F(707)-p4mup4*F(709)))
       Fa14(171) = mup3mup4*F(710)+2*(p2mup4*p4mup3*F(717)+p2mup4*p5mu
     &   p3*F(718)+p2mup4*p6mup3*F(719)-p2mup3*(p6mup4*F(720)+p5mup4*
     &   F(721))+p4mup3*p4mup4*F(722)+p4mup3*p5mup4*F(723)+p4mup4*p5m
     &   up3*F(724)-p5mup4*p6mup3*F(725)+p4mup4*p6mup3*F(726)+p2mup3*
     &   p4mup4*F(727)-p6mup3*p6mup4*F(728)+p1mup4*(p2mup3*F(712)-p4m
     &   up3*F(713)+p6mup3*F(716)+p5mup3*F(730))+p5mup3*p6mup4*F(731)
     &   -p4mup3*p6mup4*F(732)-p1mup3*(p2mup4*F(711)+p6mup4*F(714)-p4
     &   mup4*F(715)+p5mup4*F(729)-p1mup4*F(733))-p2mup3*p2mup4*F(734
     &   ))
       Fa14(172) = 8*(p4mup4*F(736)+p1mup4*F(737)+p5mup4*F(738)+p6mup4
     &   *F(739))+p2mup4*F(740)
       Fa14(173) = mup3mup4*F(741)+2*(-(p5mup4*p6mup3*F(745))-p4mup4*p
     &   6mup3*F(746)-p6mup3*p6mup4*F(747)-p2mup3*p4mup4*F(750)+p1mup
     &   4*p2mup3*F(751)-p5mup3*(p5mup4*F(742)+p6mup4*F(743)+p4mup4*F
     &   (744)-p1mup4*F(754))-p1mup3*p4mup4*F(755)+p4mup3*(p4mup4*F(7
     &   22)+p1mup4*F(752)+p5mup4*F(753)+p6mup4*F(756))+p1mup4*p6mup3
     &   *F(757)+p2mup3*p6mup4*F(758)+p2mup3*p5mup4*F(759)+p2mup4*(p4
     &   mup3*F(521)-p1mup3*F(525)+p5mup3*F(748)+p6mup3*F(749)-p2mup3
     &   *F(760)))
       Fa14(174) = -8*(p4mup4*F(762)-p1mup4*F(763)-p5mup4*F(764)-p6mup
     &   4*F(765))+p2mup4*F(766)
       Fa14(175) = p4mup3*F(767)+p1mup3*F(768)-p5mup3*F(769)-p2mup3*F(
     &   770)-p6mup3*F(771)
       Fa14(176) = mup3mup4*F(772)+2*(p1mup4*(p2mup3*F(774)-p4mup3*F(7
     &   75))+p2mup4*p4mup3*F(782)+p2mup4*p5mup3*F(783)+p2mup4*p6mup3
     &   *F(784)-p2mup3*(p4mup4*F(779)+p6mup4*F(780)+p5mup4*F(781)-p2
     &   mup4*F(785))+p4mup3*p4mup4*F(786)+p4mup3*p5mup4*F(787)+p1mup
     &   3*(-(p2mup4*F(773))+p5mup4*F(776)+p6mup4*F(777)+p4mup4*F(778
     &   )-p1mup4*F(788))+p4mup4*p5mup3*F(789)+p4mup3*p6mup4*F(790)+p
     &   4mup4*p6mup3*F(791)+p5mup3*p5mup4*F(792)+p5mup3*p6mup4*F(793
     &   )-p1mup4*p5mup3*F(794)-p5mup4*p6mup3*F(795)-p6mup3*p6mup4*F(
     &   796)+p1mup4*p6mup3*F(797))
       Fa14(177) = mup3mup4*F(798)-2*(-(p1mup3*(p5mup4*F(800)+p6mup4*F
     &   (801)+p4mup4*F(802)-p2mup4*F(806)))+p2mup3*(p4mup4*F(803)+p6
     &   mup4*F(804)+p5mup4*F(805)-p2mup4*F(807))-p4mup4*(p5mup3*F(81
     &   2)+p6mup3*F(813))+p1mup4*(p4mup3*F(799)+p1mup3*F(811)-p2mup3
     &   *F(814))-p4mup3*(p4mup4*F(808)+p5mup4*F(809)+p6mup4*F(810)-p
     &   2mup4*F(815)))
       Fa14(178) = -(p2mup4*F(818))+p4mup4*F(819)+p5mup4*F(820)+p6mup4
     &   *F(821)+p1mup4*F(822)
       Fa14(179) = p2mup3*F(824)+p1mup3*F(825)+p4mup3*F(826)+p5mup3*F(
     &   827)+p6mup3*F(828)
       Fa14(180) = 4*mup3mup4*F(829)+8*(p2mup3*p6mup4*F(834)-p2mup3*p5
     &   mup4*F(835)+p4mup3*p4mup4*F(837)+p4mup3*p5mup4*F(838)+p5mup3
     &   *p5mup4*F(839)+p5mup3*p6mup4*F(840)+p4mup4*p5mup3*F(841)-p1m
     &   up3*p4mup4*F(842)+p4mup3*p6mup4*F(843)-p1mup3*p5mup4*F(844)+
     &   p5mup4*p6mup3*F(848)+p4mup4*p6mup3*F(849)+p6mup3*p6mup4*F(85
     &   0)-p1mup3*p6mup4*F(851)+p2mup4*(p4mup3*F(830)-p1mup3*F(831)+
     &   p5mup3*F(832)-p6mup3*F(846)-p2mup3*F(852))+p2mup3*p4mup4*F(8
     &   53)+p1mup4*(p2mup3*F(833)+p4mup3*F(836)+p5mup3*F(845)+p6mup3
     &   *F(847)-p1mup3*F(854)))
       Fa14(181) = mup3mup4*F(855)-8*(p4mup4*p6mup3*F(791)-p6mup3*p6mu
     &   p4*F(796)+p1mup4*p2mup3*F(857)+p1mup3*p5mup4*F(859)-p1mup4*p
     &   5mup3*F(860)+p1mup3*p4mup4*F(861)-p2mup3*p6mup4*F(863)-p2mup
     &   3*p5mup4*F(864)+p4mup3*(p4mup4*F(786)+p6mup4*F(790)-p1mup4*F
     &   (858)+p5mup4*F(865))+p1mup3*p6mup4*F(866)+p5mup3*p5mup4*F(86
     &   7)+p5mup3*p6mup4*F(868)+p4mup4*p5mup3*F(869)+p5mup4*p6mup3*F
     &   (870)-p1mup4*p6mup3*F(871)+p1mup3*p1mup4*F(872)+p2mup4*(p4mu
     &   p3*F(782)+p6mup3*F(784)-p1mup3*F(856)+p5mup3*F(862)-p2mup3*F
     &   (873))+p2mup3*p4mup4*F(874))
       Fa14(182) = p4mup4*F(879)+p5mup4*F(880)+p6mup4*F(881)-p1mup4*F(
     &   882)-p2mup4*F(883)
       Return
       End
