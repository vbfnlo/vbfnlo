c************** Calling the Fa functions*************************
       subroutine HexBox2Fa15(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup
     &   6,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa15)
       IMPLICIT NONE
      Complex*16   p1mup3, p1mup6, p2mup3, p2mup6, p3mup3, p3mup6, 
     -          p4mup3, p4mup6, p5mup3, p5mup6, p6mup3, p6mup6
       Complex*16   mup3mup6
        common/HexBox2FhlFunctions/F
       COMMON/HexBox2Invariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12,s23,
     &   s34,s45,s56,s16,s123,s234,s345
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
      Complex*16 Fa(385),F(1646)
      Real*8 P(1185) 
        COMMON/HexBox2PFunctions/P
       Complex*16 Fa15(267:285)
       COMMON/HexBox2FaFunctions/Fa
       Fa15(267) = (p5mup6+p6mup6)*(4*((p5mup3+p6mup3)*F(650)+p2mup3*F
     &   (667)+p1mup3*F(726))+p4mup3*F(871))
       Fa15(268) = -(p1mup3*F(867))-p2mup3*F(868)+p4mup3*F(869)+(p5mup
     &   3+p6mup3)*F(1641)
       Fa15(269) = p4mup6*((p5mup3+p6mup3)*F(650)+p2mup3*F(667)+p1mup3
     &   *F(726)-p4mup3*F(728))-p2mup6*(p1mup3*F(725)-p4mup3*F(727)+(
     &   p5mup3+p6mup3)*F(730)+p2mup3*F(731))+mup3mup6*F(803)
       Fa15(270) = p1mup3*F(725)-p4mup3*F(727)+(p5mup3+p6mup3)*F(730)+
     &   p2mup3*F(731)
       Fa15(271) = (p5mup3+p6mup3)*F(650)+p2mup3*F(667)+p1mup3*F(726)-
     &   p4mup3*F(728)
       Fa15(272) = p2mup6*F(878)+p1mup6*F(879)+p6mup6*F(880)+p4mup6*F(
     &   881)+p5mup6*F(882)
       Fa15(273) = p1mup6*(p2mup3*F(883)+p1mup3*F(886))+p1mup3*(-(p2mu
     &   p6*F(883))+p5mup6*F(884)+p6mup6*F(885)+p4mup6*F(887))+p2mup3
     &   *(p2mup6*F(888)-p5mup6*F(889)-p6mup6*F(890)-p4mup6*F(891))
       Fa15(274) = p1mup3*(p2mup6*F(895)+p4mup6*F(896)+p5mup6*F(897)+p
     &   6mup6*F(898)+p1mup6*F(899))-p2mup6*p4mup3*F(900)-p2mup6*p5mu
     &   p3*F(901)-p2mup6*p6mup3*F(902)-p4mup6*p5mup3*F(903)-p4mup3*p
     &   5mup6*F(903)-p4mup6*p6mup3*F(904)-p4mup3*p6mup6*F(904)-p4mup
     &   3*p4mup6*F(905)+p5mup3*p5mup6*F(906)+p5mup6*p6mup3*F(907)+p5
     &   mup3*p6mup6*F(907)+p6mup3*p6mup6*F(908)+p2mup3*(p2mup6*F(893
     &   )-p1mup6*F(894)+p5mup6*F(909)+p6mup6*F(910)+p4mup6*F(911))+p
     &   1mup6*p5mup3*F(912)+p1mup6*p6mup3*F(913)+p1mup6*p4mup3*F(914
     &   )+mup3mup6*F(915)
       Fa15(275) = p2mup3*F(919)+p1mup3*F(920)+p6mup3*F(921)+p4mup3*F(
     &   922)+p5mup3*F(923)
       Fa15(276) = p4mup6*F(924)-p5mup6*F(925)+p6mup6*F(926)-p2mup6*F(
     &   927)+p1mup6*F(928)
       Fa15(277) = p6mup3*p6mup6*F(910)+p4mup6*p6mup3*F(911)+p4mup3*p6
     &   mup6*F(911)-p1mup6*p2mup3*F(929)-p1mup3*p1mup6*F(930)-p2mup3
     &   *p2mup6*F(931)-p1mup3*p2mup6*F(932)+p2mup3*p6mup6*F(934)+p2m
     &   up3*p4mup6*F(935)+p1mup3*p6mup6*F(937)+p1mup3*p4mup6*F(938)-
     &   mup3mup6*F(939)+p1mup6*p6mup3*F(941)+p1mup6*p4mup3*F(942)+p2
     &   mup6*p4mup3*F(943)+p2mup6*p6mup3*F(945)+p5mup3*(p6mup6*F(909
     &   )+p1mup6*F(940)+p2mup6*F(944)+p4mup6*F(946))+p5mup6*(p6mup3*
     &   F(909)+p2mup3*F(933)+p1mup3*F(936)+p4mup3*F(946)+p5mup3*F(94
     &   7))+p4mup3*p4mup6*F(948)
       Fa15(278) = p4mup6*F(949)+p5mup6*F(950)+p6mup6*F(951)+p2mup6*F(
     &   952)+p1mup6*F(953)
       Fa15(279) = p2mup6*F(954)+p4mup6*F(955)+p5mup6*F(956)+p6mup6*F(
     &   957)+p1mup6*F(958)
       Fa15(280) = p2mup6*F(959)+p1mup6*F(960)+p4mup6*F(961)+p5mup6*F(
     &   962)+p6mup6*F(963)
       Fa15(281) = p2mup6*F(964)+p4mup6*F(965)+p5mup6*F(966)+p6mup6*F(
     &   967)+p1mup6*F(968)
       Fa15(282) = -(p1mup6*F(885))-p2mup6*F(969)+p4mup6*F(970)+p5mup6
     &   *F(971)+p6mup6*F(972)
       Fa15(283) = p2mup6*F(973)-p1mup6*F(974)+p4mup6*F(975)-p5mup6*F(
     &   976)+p6mup6*F(977)
       Fa15(284) = -(p2mup6*F(978))+p4mup6*F(979)+p5mup6*F(980)+p6mup6
     &   *F(981)+p1mup6*F(982)
       Fa15(285) = p2mup6*F(983)-p1mup6*F(984)
       Return
       End
