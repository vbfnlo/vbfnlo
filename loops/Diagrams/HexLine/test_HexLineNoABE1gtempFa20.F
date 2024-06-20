c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineNoAbe1g20(p1mup2,p1mup3,p1mup4
     -   ,p1mup6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3m
     -   up6,p4mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,
     -   p6mup2,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3m
     -   up4,mup3mup6,mup4mup6,Fa20)
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
       Complex*16 Fa20(267:280)
       COMMON/test_FaFunctionsNoAbe1g/Fa
       Fa20(267) = 4*p5mup6*(mup2mup3*F(740)+2*(p3mup3*p5mup2*F(744)+p
     -   3mup2*p5mup3*F(745)+p4mup3*p5mup2*F(747)+p5mup2*p5mup3*F(747
     -   )-p1mup3*p3mup2*F(769)-p2mup3*p3mup2*F(787)-p1mup3*p5mup2*F(
     -   788)+p2mup2*p5mup3*F(793)-p2mup3*p5mup2*F(794)+p1mup2*p5mup3
     -   *F(796)-p3mup2*p3mup3*F(806)))
       Fa20(268) = mup4mup6*(p3mup3*F(341)+(p4mup3+p5mup3)*F(342)+p1mu
     -   p3*F(343)+p2mup3*F(946))
       Fa20(269) = (p4mup3+p5mup3)*F(957)+p2mup3*F(958)+2*(p1mup3*F(34
     -   3)+p3mup3*F(959))
       Fa20(270) = mup2mup6*p4mup3*F(960)-mup3mup6*(-(p4mup2*F(961))+4
     -   *s45*(p1mup2*F(327)+p2mup2*F(328)+p5mup2*F(332)+p3mup2*F(980
     -   )))+4*(mup2mup3*p5mup6*s45*F(330)+mup2mup6*p1mup3*s45*F(346)
     -   +mup2mup6*p3mup3*s45*F(348)+mup2mup3*p1mup6*s45*F(349)+mup2m
     -   up6*p5mup3*s45*F(350)+mup2mup3*p3mup6*s45*F(351)+mup2mup6*p2
     -   mup3*s45*F(973)+mup2mup3*p4mup6*F(984)+p2mup6*(mup2mup3*s45*
     -   F(331)-2*(p4mup3*(p3mup2*F(887)-p1mup2*F(922)+p2mup2*F(928)+
     -   p5mup2*F(934))-p4mup2*(p5mup3*F(306)-p1mup3*F(876)+p3mup3*F(
     -   886)+p2mup3*F(887)+p4mup3*F(985)))))-8*(p3mup6*p4mup3*p5mup2
     -   *F(369)+p1mup3*p4mup2*p5mup6*F(874)+p1mup3*p1mup6*p4mup2*F(8
     -   75)+p1mup3*p3mup6*p4mup2*F(877)-p1mup6*p2mup3*p4mup2*F(878)-
     -   p2mup3*p4mup2*p5mup6*F(882)-p3mup3*p4mup2*p5mup6*F(883)-p3mu
     -   p3*p3mup6*p4mup2*F(884)-p2mup3*p3mup6*p4mup2*F(885)+p3mup2*p
     -   3mup6*p4mup3*F(885)-p1mup6*p3mup3*p4mup2*F(888)-p4mup2*p5mup
     -   3*p5mup6*F(890)-p3mup6*p4mup2*p5mup3*F(891)-p1mup6*p4mup2*p5
     -   mup3*F(892)-p1mup2*p4mup3*p5mup6*F(920)-p1mup2*p1mup6*p4mup3
     -   *F(921)-p1mup2*p3mup6*p4mup3*F(923)+p2mup2*p4mup3*p5mup6*F(9
     -   26)+p1mup6*p2mup2*p4mup3*F(927)+p2mup2*p3mup6*p4mup3*F(929)+
     -   p1mup6*p3mup2*p4mup3*F(931)+p3mup2*p4mup3*p5mup6*F(932)+p1mu
     -   p6*p4mup3*p5mup2*F(935)+p4mup3*p5mup2*p5mup6*F(936)-p1mup2*p
     -   4mup3*p4mup6*F(962)+p1mup3*p4mup2*p4mup6*F(963)+p2mup2*p4mup
     -   3*p4mup6*F(964)-p2mup3*p4mup2*p4mup6*F(974)+p3mup2*p4mup3*p4
     -   mup6*F(975)-p3mup3*p4mup2*p4mup6*F(976)-p4mup2*p4mup6*p5mup3
     -   *F(982)+p4mup3*p4mup6*p5mup2*F(983)-p3mup6*p4mup2*p4mup3*F(9
     -   86)-p4mup2*p4mup3*p4mup6*F(988)+p1mup6*p4mup2*p4mup3*F(989)-
     -   p4mup2*p4mup3*p5mup6*F(990))-32*p2mup3*p3mup2*p4mup6*F(991)+
     -   16*(p2mup2*p3mup3*p4mup6*F(965)+p2mup2*p4mup6*p5mup3*F(966)+
     -   p2mup2*p2mup3*p4mup6*F(967)+p1mup2*p1mup3*p4mup6*F(968)+p1mu
     -   p2*p2mup3*p4mup6*F(969)+p1mup3*p2mup2*p4mup6*F(970)+p1mup2*p
     -   3mup3*p4mup6*F(971)+p1mup2*p4mup6*p5mup3*F(972)-p2mup3*p4mup
     -   6*p5mup2*F(977)+p3mup2*p4mup6*p5mup3*F(978)-p3mup3*p4mup6*p5
     -   mup2*F(979)-p1mup3*p4mup6*p5mup2*F(981)-p4mup6*p5mup2*p5mup3
     -   *F(987)-p3mup2*p3mup3*p4mup6*F(992)-p1mup3*p3mup2*p4mup6*F(9
     -   93))
       Fa20(271) = -(mup3mup6*F(994))-4*(p1mup3*((p4mup6+p5mup6)*F(874
     -   )+p1mup6*F(875)+p3mup6*F(877))-p2mup6*((p4mup3+p5mup3)*F(306
     -   )-p1mup3*F(876)+p3mup3*F(886)+p2mup3*F(887))-p3mup3*((p4mup6
     -   +p5mup6)*F(883)+p3mup6*F(884)+p1mup6*F(888))-(p4mup3+p5mup3)
     -   *((p4mup6+p5mup6)*F(890)+p3mup6*F(891)+p1mup6*F(892))-p2mup3
     -   *((p4mup6+p5mup6)*F(882)+p3mup6*F(885)-p1mup6*F(995)))
       Fa20(272) = -2*p1mup6*F(284)+p3mup6*F(996)+p2mup6*F(997)+(p4mup
     -   6+p5mup6)*F(998)
       Fa20(273) = -2*(p1mup6*F(336)-p3mup6*F(338))+p2mup6*F(999)+(p4m
     -   up6+p5mup6)*F(1000)
       Fa20(274) = -4*(p1mup2*(p4mup6+p5mup6)*F(920)+p1mup2*p1mup6*F(9
     -   21)-p2mup2*((p4mup6+p5mup6)*F(926)+p1mup6*F(927))-p2mup6*(p3
     -   mup2*F(887)-p1mup2*F(922)+p2mup2*F(928))-p3mup6*((p4mup2+p5m
     -   up2)*F(369)+p3mup2*F(885)-p1mup2*F(923)+p2mup2*F(929))-p1mup
     -   6*p3mup2*F(931)-p3mup2*p4mup6*F(932)-p3mup2*p5mup6*F(932))+m
     -   up2mup6*F(1001)+p4mup2*p4mup6*F(1002)+p4mup6*p5mup2*F(1002)+
     -   p4mup2*p5mup6*F(1002)+p5mup2*p5mup6*F(1002)+p2mup6*(p4mup2+p
     -   5mup2)*F(1003)+p1mup6*(p4mup2+p5mup2)*F(1004)
       Fa20(275) = p5mup4*(p2mup2*(p3mup3*F(965)+(p4mup3+p5mup3)*F(966
     -   )+p2mup3*F(967)+p1mup3*F(970))+p1mup2*(p1mup3*F(968)+p2mup3*
     -   F(969)+p3mup3*F(971)+(p4mup3+p5mup3)*F(972))-p2mup3*p4mup2*F
     -   (977)-p2mup3*p5mup2*F(977)+p3mup2*p4mup3*F(978)+p3mup2*p5mup
     -   3*F(978)-p3mup3*p4mup2*F(979)-p3mup3*p5mup2*F(979)-p1mup3*p4
     -   mup2*F(981)-p1mup3*p5mup2*F(981)-p4mup2*p4mup3*F(987)-p4mup3
     -   *p5mup2*F(987)-p4mup2*p5mup3*F(987)-p5mup2*p5mup3*F(987)-2*p
     -   2mup3*p3mup2*F(991)-p3mup2*p3mup3*F(992)-p1mup3*p3mup2*F(993
     -   )+mup2mup3*F(1005))
       Fa20(276) = mup4mup6*(p2mup2*(p3mup3*F(965)+(p4mup3+p5mup3)*F(9
     -   66)+p2mup3*F(967)+p1mup3*F(970))+p1mup2*(p1mup3*F(968)+p2mup
     -   3*F(969)+p3mup3*F(971)+(p4mup3+p5mup3)*F(972))-p2mup3*p4mup2
     -   *F(977)-p2mup3*p5mup2*F(977)+p3mup2*p4mup3*F(978)+p3mup2*p5m
     -   up3*F(978)-p3mup3*p4mup2*F(979)-p3mup3*p5mup2*F(979)-p1mup3*
     -   p4mup2*F(981)-p1mup3*p5mup2*F(981)-p4mup2*p4mup3*F(987)-p4mu
     -   p3*p5mup2*F(987)-p4mup2*p5mup3*F(987)-p5mup2*p5mup3*F(987)-2
     -   *p2mup3*p3mup2*F(991)-p3mup2*p3mup3*F(992)-p1mup3*p3mup2*F(9
     -   93)+mup2mup3*F(1005))
       Fa20(277) = p2mup2*(p3mup3*F(965)+(p4mup3+p5mup3)*F(966)+p2mup3
     -   *F(967)+p1mup3*F(970))+p1mup2*(p1mup3*F(968)+p2mup3*F(969)+p
     -   3mup3*F(971)+(p4mup3+p5mup3)*F(972))-p2mup3*p4mup2*F(977)-p2
     -   mup3*p5mup2*F(977)+p3mup2*p4mup3*F(978)+p3mup2*p5mup3*F(978)
     -   -p3mup3*p4mup2*F(979)-p3mup3*p5mup2*F(979)-p1mup3*p4mup2*F(9
     -   81)-p1mup3*p5mup2*F(981)-p4mup2*p4mup3*F(987)-p4mup3*p5mup2*
     -   F(987)-p4mup2*p5mup3*F(987)-p5mup2*p5mup3*F(987)-2*p2mup3*p3
     -   mup2*F(991)-p3mup2*p3mup3*F(992)+mup2mup3*F(1005)+p1mup3*p3m
     -   up2*F(1006)
       Fa20(278) = -2*mup4mup6*F(1007)+(p1mup4+p2mup4)*p4mup6*F(1015)+
     -   p1mup4*p3mup6*F(1016)+p2mup4*p3mup6*F(1016)-p5mup4*p5mup6*F(
     -   1017)-(p1mup6+p2mup6)*p5mup4*F(1018)-p3mup6*p5mup4*F(1019)-p
     -   4mup6*p5mup4*F(1020)+p1mup4*p5mup6*F(1021)+p2mup4*p5mup6*F(1
     -   021)+p4mup4*(p5mup6*F(1009)+4*((p1mup6+p2mup6)*F(1013)+p3mup
     -   6*F(1014)+p4mup6*F(1022)))+4*(p3mup4*(p5mup6*F(1008)-(p1mup6
     -   +p2mup6)*F(1010)+p3mup6*F(1011)+p4mup6*F(1012))+(p1mup4+p2mu
     -   p4)*(p1mup6+p2mup6)*F(1023))
       Fa20(279) = mup3mup6*F(1024)+p1mup3*p4mup6*F(1027)+p2mup3*p4mup
     -   6*F(1027)+(p1mup6+p2mup6)*p4mup3*F(1028)+p4mup3*p5mup6*F(102
     -   9)+p3mup3*(p4mup6*F(1025)+p3mup6*F(1030)-4*(p1mup6+p2mup6)*F
     -   (1032)+p5mup6*F(1033))+4*(p4mup3*p4mup6*F(1034)+p3mup6*(p4mu
     -   p3*F(1026)+(p1mup3+p2mup3)*F(1031)+p5mup3*F(1035))+p4mup6*p5
     -   mup3*F(1036)+p1mup6*p5mup3*F(1037)+p2mup6*p5mup3*F(1037)+p1m
     -   up3*p1mup6*F(1038)+p1mup6*p2mup3*F(1038)+p1mup3*p2mup6*F(103
     -   8)+p2mup3*p2mup6*F(1038)+p1mup3*p5mup6*F(1039)+p2mup3*p5mup6
     -   *F(1039)-p5mup3*p5mup6*F(1040))
       Fa20(280) = 2*(p4mup6*F(1041)+p3mup6*F(1042)+p5mup6*F(1043))+(p
     -   1mup6+p2mup6)*F(1044)
       Return
       End
