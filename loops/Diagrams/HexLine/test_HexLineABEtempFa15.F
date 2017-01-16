c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineAbe15(p1mup2,p1mup3,p1mup4,p1m
     -   up5,p2mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,
     -   p4mup2,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mu
     -   p2,p6mup3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,
     -   mup3mup5,mup4mup5,Fa15)
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
       Complex*16 Fa15(155:165)
       COMMON/test_FaFunctionsAbe/Fa
       Fa15(155) = p1mup4*F(781)+p2mup4*F(782)+4*(p3mup4*F(783)+p4mup4
     -   *F(784))
       Fa15(156) = mup2mup4*F(785)-8*(p1mup2*(p1mup4*F(786)+p2mup4*F(7
     -   87)+p3mup4*F(788)+p4mup4*F(789))+p1mup4*p2mup2*F(790)+p2mup2
     -   *p2mup4*F(791)+p2mup2*p3mup4*F(792)+p2mup2*p4mup4*F(793)+p1m
     -   up4*p3mup2*F(794)+p2mup4*p3mup2*F(795)+p3mup2*p3mup4*F(796)+
     -   p3mup2*p4mup4*F(797)-p3mup4*p4mup2*F(798)-p4mup2*p4mup4*F(79
     -   9)+p1mup4*p4mup2*F(800)+p2mup4*p4mup2*F(801))
       Fa15(157) = (p1mup2+p2mup2)*F(802)+s12*(p4mup2*F(803)+p3mup2*F(
     -   804))
       Fa15(158) = p4mup3*F(805)+p3mup3*F(806)-4*((p1mup3+p2mup3)*F(80
     -   7)-p6mup3*F(808))
       Fa15(159) = p1mup2*(p4mup3*F(809)+8*((p1mup3+p2mup3)*F(810)-p6m
     -   up3*F(811)-p3mup3*F(812)))+p2mup2*(p4mup3*F(809)+8*((p1mup3+
     -   p2mup3)*F(810)-p6mup3*F(811)-p3mup3*F(812)))+8*p4mup2*s12*(p
     -   6mup3*F(813)+3*((p1mup3+p2mup3)*F(814)-p3mup3*F(815)+p4mup3*
     -   F(816)))
       Fa15(160) = p4mup2*((p1mup3+p2mup3)*F(818)+p3mup3*F(819)+p6mup3
     -   *F(820)+p4mup3*F(823))+p3mup2*(p4mup3*F(817)+8*(-((p1mup3+p2
     -   mup3)*F(821))+p6mup3*F(822)+p3mup3*F(824)))
       Fa15(161) = 4*((p1mup3+p2mup3)*F(810)-p6mup3*F(811)-p3mup3*F(81
     -   2))+p4mup3*F(825)
       Fa15(162) = (p1mup2+p2mup2)*(p1mup3*p4mup4*F(826)+p2mup3*p4mup4
     -   *F(826)+p3mup3*p4mup4*F(827)+p1mup3*p3mup4*F(828)+p2mup3*p3m
     -   up4*F(828)+p1mup4*p3mup3*F(829)+p2mup4*p3mup3*F(829)+p3mup3*
     -   p3mup4*F(830)+(p1mup3+p2mup3)*(p1mup4+p2mup4)*F(831))
       Fa15(163) = p6mup4*(p1mup2*(p4mup3*F(832)-16*((p1mup3+p2mup3)*F
     -   (833)-p3mup3*F(838)))+p2mup2*(p4mup3*F(832)+p6mup3*F(834)-16
     -   *((p1mup3+p2mup3)*F(833)-p3mup3*F(838)))+16*s12*(p2mup3*p4mu
     -   p2*F(835)+p3mup3*p4mup2*F(837)+p4mup2*p4mup3*F(839)-p2mup3*p
     -   3mup2*F(841)-p3mup2*(p4mup3*F(836)+p6mup3*F(840)+p1mup3*F(84
     -   1))+p3mup2*p3mup3*F(842)))
       Fa15(164) = p2mup2*(-8*mup3mup4*F(843)+p1mup4*p4mup3*F(844)+p2m
     -   up4*p4mup3*F(844)+p3mup4*p4mup3*F(845)+p3mup4*p6mup3*F(846)+
     -   p4mup4*p6mup3*F(847)+p4mup3*p4mup4*F(848)+(p1mup4+p2mup4)*p6
     -   mup3*F(849))+p1mup2*(-8*mup3mup4*F(843)+p4mup3*((p1mup4+p2mu
     -   p4)*F(844)+p3mup4*F(845)+p4mup4*F(848))+p6mup3*(p6mup4*F(834
     -   )+p3mup4*F(846)+p4mup4*F(847)+(p1mup4+p2mup4)*F(849)))+16*s1
     -   2*(p2mup3*p3mup2*p3mup4*F(850)+p1mup3*(p4mup2*p6mup4*F(835)+
     -   p3mup2*p3mup4*F(850))-p4mup2*p6mup3*p6mup4*F(851)-p3mup2*p3m
     -   up3*((p1mup4+p2mup4)*F(852)+p3mup4*F(853)))
       Fa15(165) = mup3mup4*(-(p4mup2*F(854))+p3mup2*F(855))-2*(p3mup2
     -   *(p2mup4*p6mup3*F(857)+p4mup4*p6mup3*F(858)+p1mup3*p2mup4*F(
     -   859)+p2mup3*p2mup4*F(859)+p1mup3*p4mup4*F(861)+p2mup3*p4mup4
     -   *F(861)-p2mup4*p4mup3*F(863)+p1mup4*(p6mup3*F(857)+(p1mup3+p
     -   2mup3)*F(859)-p4mup3*F(863))+p3mup4*(p6mup3*F(856)-p4mup3*F(
     -   865))-p4mup3*p4mup4*F(868)-p3mup3*p4mup4*F(873))+p4mup2*(p2m
     -   up3*p2mup4*F(860)-p2mup3*p3mup4*F(862)+p2mup4*p3mup3*F(864)+
     -   p3mup3*p3mup4*F(866)+p2mup4*p4mup3*F(867)+p2mup3*p4mup4*F(86
     -   7)+p1mup3*((p1mup4+p2mup4)*F(860)-p3mup4*F(862)+p4mup4*F(867
     -   ))+p3mup3*p4mup4*F(869)+p3mup4*p4mup3*F(870)+p2mup4*p6mup3*F
     -   (871)+p1mup4*(p2mup3*F(860)+p3mup3*F(864)+p4mup3*F(867)+p6mu
     -   p3*F(871))+p3mup4*p6mup3*F(872)+p4mup3*p4mup4*F(874)+p4mup4*
     -   p6mup3*F(875)))
       Return
       End
