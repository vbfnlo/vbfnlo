c************** Calling the Fa functions*************************
       subroutine HexBoxCross1Fa15(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,
     &   p3mup4,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa
     &   15)
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
       Complex*16 Fa15(183:195)
       COMMON/HexBoxCross1FaFunctions/Fa
       Fa15(183) = p2mup3*F(888)+4*(p4mup3*F(885)+p5mup3*F(886)+p6mup3
     &   *F(887)+p1mup3*F(889))
       Fa15(184) = p2mup3*F(893)+4*(p1mup3*F(890)+p4mup3*F(891)-p5mup3
     &   *F(892)-p6mup3*F(894))
       Fa15(185) = p1mup4*F(552)+p2mup4*F(898)+p4mup4*F(899)+p5mup4*F(
     &   900)-p6mup4*F(1460)
       Fa15(186) = p1mup3*F(901)+4*(p2mup3*F(902)-p4mup3*F(903)-p5mup3
     &   *F(904)+p6mup3*F(905))
       Fa15(187) = mup3mup4*F(910)+2*(-(p1mup3*(p5mup4*F(912)+p6mup4*F
     &   (913)+p1mup4*F(915)+p2mup4*F(917)))+p4mup4*(p6mup3*F(784)-p1
     &   mup3*F(914)+p5mup3*F(916)+p2mup3*F(918))-p5mup3*p5mup4*F(923
     &   )-p5mup3*p6mup4*F(924)+p1mup4*p5mup3*F(926)-p2mup3*(-(p6mup4
     &   *F(919))-p5mup4*F(920)+p1mup4*F(925)+p2mup4*F(927))+p4mup3*(
     &   p4mup4*F(782)+p1mup4*F(911)-p5mup4*F(921)-p6mup4*F(922)+p2mu
     &   p4*F(928))+p2mup4*p5mup3*F(929)-p5mup4*p6mup3*F(930)-p6mup3*
     &   p6mup4*F(931)+p1mup4*p6mup3*F(932)+p2mup4*p6mup3*F(933))
       Fa15(188) = p5mup4*F(742)+p6mup4*F(743)-p1mup4*F(754)+p4mup4*F(
     &   935)-p2mup4*F(936)
       Fa15(189) = p5mup4*F(723)-p1mup4*F(937)+p6mup4*F(938)-p2mup4*F(
     &   939)+p4mup4*F(940)
       Fa15(190) = p4mup4*F(941)+p5mup4*F(942)+p6mup4*F(943)-p1mup4*F(
     &   944)-p2mup4*F(945)
       Fa15(191) = p2mup3*F(947)-p4mup3*F(948)+p6mup3*F(949)+p1mup3*F(
     &   950)+p5mup3*F(951)
       Fa15(192) = 2*(p5mup3*F(812)+p6mup3*F(813))+p1mup3*F(953)-p2mup
     &   3*F(954)+p4mup3*F(955)
       Fa15(193) = p4mup4*F(959)+4*p1mup4*F(960)+p5mup4*F(961)+p6mup4*
     &   F(962)+p2mup4*F(963)
       Fa15(194) = p4mup4*F(964)-4*p1mup4*F(965)+p5mup4*F(966)+p6mup4*
     &   F(967)+p2mup4*F(968)
       Fa15(195) = p2mup4*F(972)-p4mup4*F(973)+p6mup4*F(974)+p5mup4*F(
     &   975)+p1mup4*F(976)
       Return
       End
