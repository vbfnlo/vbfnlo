c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineNoAbe1g18(p1mup2,p1mup3,p1mup4
     -   ,p1mup6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3m
     -   up6,p4mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,
     -   p6mup2,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3m
     -   up4,mup3mup6,mup4mup6,Fa18)
       IMPLICIT NONE
      Complex*16   p1mup2, p1mup3, p1mup4, p1mup6, p2mup2, p2mup3, 
     -          p2mup4, p2mup6, p3mup2, p3mup3, p3mup4, p3mup6, 
     -          p4mup2, p4mup3, p4mup4, p4mup6, p5mup2, p5mup3, 
     -          p5mup4, p5mup6, p6mup2, p6mup3, p6mup4, p6mup6
       Complex*16   mup2mup3, mup2mup4, mup2mup6, mup3mup4, mup3mup6, 
     -          mup4mup6
        common/test_FhlFunctionsNoAbe1g/F
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
      Complex*16 Fa(587),F(6395)
      Real*8 P(7129) 
        COMMON/test_PFunctionsNoAbe1g/P
       Complex*16 Fa18(239:252)
       COMMON/test_FaFunctionsNoAbe1g/Fa
       Fa18(239) = p2mup3*p4mup2*F(813)+p1mup3*p4mup2*F(815)+4*(p3mup2
     -   *p5mup3*s45*F(714)+p5mup2*p5mup3*s45*F(716)+p1mup2*p1mup3*s4
     -   5*F(721)+p1mup2*p2mup3*s45*F(722)+p2mup3*p5mup2*s45*F(727)+p
     -   1mup2*p5mup3*s45*F(732)+p1mup3*p5mup2*s45*F(734)-p3mup2*p4mu
     -   p3*F(809)+p4mup2*p5mup3*F(810)-p4mup3*p5mup2*F(811)+p2mup2*(
     -   p2mup3*s45*F(713)+p1mup3*s45*F(723)+p5mup3*s45*F(725)-p4mup3
     -   *F(812))-p1mup2*p4mup3*F(814)+p4mup2*p4mup3*F(816)+p2mup3*p3
     -   mup2*s45*F(817)+p3mup3*(p4mup2*F(710)+s45*(p5mup2*F(715)+p2m
     -   up2*F(720)+p1mup2*F(724)+p3mup2*F(818)))+p1mup3*p3mup2*s45*F
     -   (819)+mup2mup3*s45*F(820))
       Fa18(240) = mup4mup6*(8*(p3mup2*(p4mup3+p5mup3)*F(822)+p3mup3*(
     -   p4mup2+p5mup2)*F(823))+p1mup2*p3mup3*F(827)+p1mup3*p3mup2*F(
     -   828)+p2mup3*((p4mup2+p5mup2)*F(825)+p1mup2*F(826)+p3mup2*F(8
     -   29))+p3mup2*p3mup3*F(830)+(p4mup2+p5mup2)*(8*(p4mup3+p5mup3)
     -   *F(831)+p1mup3*F(833))+mup2mup3*F(834)+p1mup2*p1mup3*F(835)+
     -   p2mup2*(p2mup3*F(821)+p3mup3*F(824)+8*(p4mup3+p5mup3)*F(832)
     -   +p1mup3*F(836))+p1mup2*(p4mup3+p5mup3)*F(837))
       Fa18(241) = p3mup2*(p4mup3+p5mup3)*F(839)+p3mup3*(p4mup2+p5mup2
     -   )*F(840)+p1mup2*p3mup3*F(844)+p1mup3*p3mup2*F(845)+p2mup3*((
     -   p4mup2+p5mup2)*F(842)+p1mup2*F(843)+p3mup2*F(846))+p3mup2*p3
     -   mup3*F(847)+(p4mup2+p5mup2)*((p4mup3+p5mup3)*F(848)+p1mup3*F
     -   (850))+mup2mup3*F(851)+p1mup2*p1mup3*F(852)+p2mup2*(p2mup3*F
     -   (838)+p3mup3*F(841)+(p4mup3+p5mup3)*F(849)+p1mup3*F(853))+p1
     -   mup2*(p4mup3+p5mup3)*F(854)
       Fa18(242) = p5mup4*(p3mup2*(p4mup3+p5mup3)*F(856)+p1mup3*p3mup2
     -   *F(862)+p2mup3*((p4mup2+p5mup2)*F(855)+p2mup2*F(859)+p1mup2*
     -   F(860)+p3mup2*F(863))+p3mup3*((p4mup2+p5mup2)*F(857)+p2mup2*
     -   F(858)+p1mup2*F(861)+p3mup2*F(864))+p4mup2*p4mup3*F(865)+p4m
     -   up3*p5mup2*F(865)+p4mup2*p5mup3*F(865)+p5mup2*p5mup3*F(865)+
     -   p2mup2*p4mup3*F(866)+p2mup2*p5mup3*F(866)+p1mup3*p4mup2*F(86
     -   7)+p1mup3*p5mup2*F(867)+mup2mup3*F(868)+p1mup2*p1mup3*F(869)
     -   +p1mup3*p2mup2*F(870)+p1mup2*(p4mup3+p5mup3)*F(871))
       Fa18(243) = mup3mup6*p5mup4*F(872)+mup4mup6*p4mup3*F(873)+mup4m
     -   up6*p2mup3*s45*F(879)+mup4mup6*p1mup3*s45*F(880)+mup4mup6*p3
     -   mup3*s45*F(881)-mup4mup6*p5mup3*s45*F(889)+2*(-(p1mup3*p5mup
     -   4*((p4mup6+p5mup6)*F(874)+p1mup6*F(875)+p3mup6*F(877)))+p1mu
     -   p6*p2mup3*p5mup4*F(878)+p2mup3*p4mup6*p5mup4*F(882)+p2mup3*p
     -   5mup4*p5mup6*F(882)+p3mup3*p4mup6*p5mup4*F(883)+p3mup3*p5mup
     -   4*p5mup6*F(883)+p3mup3*p3mup6*p5mup4*F(884)+p2mup3*p3mup6*p5
     -   mup4*F(885)+p2mup6*p5mup4*((p4mup3+p5mup3)*F(306)-p1mup3*F(8
     -   76)+p3mup3*F(886)+p2mup3*F(887))+p1mup6*p3mup3*p5mup4*F(888)
     -   +p4mup3*p4mup6*p5mup4*F(890)+p4mup6*p5mup3*p5mup4*F(890)+p4m
     -   up3*p5mup4*p5mup6*F(890)+p5mup3*p5mup4*p5mup6*F(890)+p3mup6*
     -   p4mup3*p5mup4*F(891)+p3mup6*p5mup3*p5mup4*F(891)+p1mup6*(p4m
     -   up3+p5mup3)*p5mup4*F(892))
       Fa18(244) = 2*p5mup3*F(889)+(p4mup3*F(893))/s45+p2mup3*F(894)+p
     -   3mup3*F(895)+p1mup3*F(896)
       Fa18(245) = mup4mup6*(p2mup3*F(897)+p3mup3*F(898)+(p4mup3+p5mup
     -   3)*F(899)+p1mup3*F(900))
       Fa18(246) = p2mup3*F(901)+(p4mup3+p5mup3)*F(903)-2*(p3mup3*F(90
     -   2)-p1mup3*F(904))
       Fa18(247) = p5mup4*(p2mup3*F(897)+p3mup3*F(898)+(p4mup3+p5mup3)
     -   *F(899)+p1mup3*F(900))
       Fa18(248) = mup4mup6*(p1mup2*F(320)+p2mup2*F(321)+(p4mup2+p5mup
     -   2)*F(323)+p3mup2*F(905))
       Fa18(249) = p3mup2*F(906)+p2mup2*F(907)+(p4mup2+p5mup2)*F(908)+
     -   p1mup2*F(909)
       Fa18(250) = mup4mup6*((p4mup2+p5mup2)*F(910)+4*(p1mup2*F(911)+p
     -   2mup2*F(912)-p3mup2*F(913)))
       Fa18(251) = p2mup2*F(914)+(p4mup2+p5mup2)*F(915)+p3mup2*F(916)+
     -   p1mup2*F(917)
       Fa18(252) = p5mup4*(p1mup2*F(320)+p2mup2*F(321)+(p4mup2+p5mup2)
     -   *F(323)+p3mup2*F(905))
       Return
       End
