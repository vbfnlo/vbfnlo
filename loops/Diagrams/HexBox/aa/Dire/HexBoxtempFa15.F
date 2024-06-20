c************** Calling the Fa functions*************************
       subroutine HexBoxFa15(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4
     &   ,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa15)
       IMPLICIT NONE
      Complex*16   p1mup3, p1mup4, p2mup3, 
     -          p2mup4, p3mup3, p3mup4, 
     -          p4mup3, p4mup4, p5mup3, 
     -          p5mup4, p6mup3, p6mup4
       Complex*16   mup3mup4

        common/HexBoxFhlFunctions/F
       COMMON/HexBoxInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12,s23,s
     &   34,s45,s56,s16,s123,s234,s345
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
      Complex*16 Fa(259),F(1314)
      Real*8 P(1362) 
        COMMON/HexBoxPFunctions/P
       Complex*16 Fa15(169:180)
       COMMON/HexBoxFaFunctions/Fa
       Fa15(169) = p6mup3*F(724)+p4mup3*F(725)+p5mup3*F(726)+p2mup3*F(
     &   727)+p1mup3*F(728)
       Fa15(170) = p1mup3*F(729)-4*(p2mup3*F(730)-p6mup3*F(731)-p4mup3
     &   *F(732))+p5mup3*F(733)
       Fa15(171) = p1mup3*F(734)+p6mup3*F(735)+p4mup3*F(736)+p5mup3*F(
     &   737)+p2mup3*F(738)
       Fa15(172) = 4*mup3mup4*F(739)+8*(-(p4mup4*p5mup3*F(746))-p5mup3
     &   *p5mup4*F(747)-p5mup3*p6mup4*F(748)-p1mup3*p4mup4*F(749)-p1m
     &   up3*p5mup4*F(750)-p1mup3*p6mup4*F(751)+p4mup4*p6mup3*F(752)+
     &   p5mup4*p6mup3*F(753)+p6mup3*p6mup4*F(754)-p1mup3*p2mup4*F(75
     &   5)+p2mup4*p5mup3*F(757)+p2mup4*p6mup3*F(758)+p1mup3*p1mup4*F
     &   (759)-p4mup3*(p4mup4*F(743)+p5mup4*F(744)+p6mup4*F(745)-p2mu
     &   p4*F(756)+p1mup4*F(760))-p1mup4*p5mup3*F(761)+p1mup4*p6mup3*
     &   F(762)+p2mup3*(p4mup4*F(740)+p5mup4*F(741)+p6mup4*F(742)-p1m
     &   up4*F(755)-p2mup4*F(763)))
       Fa15(173) = 4*mup3mup4*F(764)-8*(p4mup3*p5mup4*F(657)+p4mup4*(p
     &   4mup3*F(656)+p5mup3*F(659))+p5mup3*p5mup4*F(661)+p5mup3*p6mu
     &   p4*F(662)+p5mup4*p6mup3*F(663)+p6mup3*p6mup4*F(664)+p2mup4*(
     &   p4mup3*F(617)+p5mup3*F(619)+p6mup3*F(620)-p2mup3*F(628)-p1mu
     &   p3*F(765))-p2mup3*(p5mup4*F(654)+p6mup4*F(655)+p4mup4*F(668)
     &   -p1mup4*F(766))-p1mup3*p4mup4*F(767)-p1mup4*p4mup3*F(768)-p1
     &   mup3*p5mup4*F(769)-p1mup3*p6mup4*F(770)-p1mup4*p5mup3*F(771)
     &   +p1mup4*p6mup3*F(772)-p1mup3*p1mup4*F(773)-p4mup3*p6mup4*F(7
     &   74)-p4mup4*p6mup3*F(775))
       Fa15(174) = p6mup4*F(778)+p4mup4*F(779)+p5mup4*F(780)-p1mup4*F(
     &   781)-p2mup4*F(782)
       Fa15(175) = p1mup4*F(783)+p2mup4*F(784)+p4mup4*F(785)+p5mup4*F(
     &   786)-p6mup4*F(1313)
       Fa15(176) = p6mup4*F(787)+p4mup4*F(788)+p5mup4*F(789)-p1mup4*F(
     &   790)-p2mup4*F(791)
       Fa15(177) = p2mup3*F(782)+p1mup3*F(792)+p6mup3*F(793)+p4mup3*F(
     &   794)+p5mup3*F(795)
       Fa15(178) = -(p1mup3*F(783))+p2mup3*F(796)+p4mup3*F(797)+p5mup3
     &   *F(798)+p6mup3*F(1313)
       Fa15(179) = p2mup3*F(795)-p1mup3*F(799)+p6mup3*F(800)+p4mup3*F(
     &   801)+p5mup3*F(802)
       Fa15(180) = 2*(mup3mup4*F(824)+2*(-(p2mup4*p6mup3*F(793))-p2mup
     &   4*p4mup3*F(794)-p2mup4*p5mup3*F(795)-p5mup4*p6mup3*F(800)-p4
     &   mup3*p5mup4*F(801)-p5mup3*p5mup4*F(802)+p1mup3*(-(p2mup4*F(7
     &   92))+p5mup4*F(799)+p4mup4*F(825)+p6mup4*F(826)-p1mup4*F(827)
     &   )-p4mup3*p4mup4*F(830)-p4mup3*p6mup4*F(831)-p4mup4*p5mup3*F(
     &   832)-p4mup4*p6mup3*F(833)-p5mup3*p6mup4*F(834)-p6mup3*p6mup4
     &   *F(835)+p2mup3*(p5mup4*F(791)+p4mup4*F(828)+p6mup4*F(829)-p2
     &   mup4*F(836)-p1mup4*F(837))-p1mup4*p4mup3*F(838)-p1mup4*p5mup
     &   3*F(839)-p1mup4*p6mup3*F(840)))
       Return
       End
