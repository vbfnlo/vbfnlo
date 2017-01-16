c************** Calling the Fa functions*************************
       subroutine HexBoxCross1Fa16(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,
     &   p3mup4,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa
     &   16)
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
       Complex*16 Fa16(196:208)
       COMMON/HexBoxCross1FaFunctions/Fa
       Fa16(196) = -(p1mup4*F(977))+p2mup4*F(978)+2*(p5mup4*F(725)+p6m
     &   up4*F(979))+p4mup4*F(980)
       Fa16(197) = mup3mup4*F(981)-2*(p4mup3*p5mup4*F(986)+p5mup3*p5mu
     &   p4*F(987)+p5mup3*p6mup4*F(988)+p4mup3*p6mup4*F(989)+p5mup4*p
     &   6mup3*F(990)-p2mup3*p6mup4*F(994)-p2mup3*p5mup4*F(995)+p6mup
     &   3*p6mup4*F(996)+p1mup4*(p5mup3*F(984)+p6mup3*F(985)-p2mup3*F
     &   (993)-p4mup3*F(997))+p2mup3*p2mup4*F(998)+p2mup4*p5mup3*F(99
     &   9)+p2mup4*p6mup3*F(1000)+p1mup3*(-(p5mup4*F(982))-p6mup4*F(9
     &   83)+p2mup4*F(991)+p1mup4*F(1003))+p4mup4*(p5mup3*F(649)-p2mu
     &   p3*F(992)+p4mup3*F(1001)+p6mup3*F(1002)-p1mup3*F(1004))-p2mu
     &   p4*p4mup3*F(1005))
       Fa16(198) = 2*(p1mup3*F(525)+p2mup3*F(530)+p4mup3*F(760))-p5mup
     &   3*F(1006)-p6mup3*F(1007)
       Fa16(199) = -(p1mup3*F(1008))+p4mup3*F(1009)+p5mup3*F(1010)+p2m
     &   up3*F(1011)+p6mup3*F(1012)
       Fa16(200) = -(p1mup3*F(1013))+p4mup3*F(1014)+p5mup3*F(1015)+p6m
     &   up3*F(1016)-p2mup3*F(1017)
       Fa16(201) = -(p5mup4*F(1018))-p6mup4*F(1019)+p4mup4*F(1020)+p1m
     &   up4*F(1021)+p2mup4*F(1022)
       Fa16(202) = -(p1mup4*F(1023))+p4mup4*F(1024)+p5mup4*F(1025)+p6m
     &   up4*F(1026)-p2mup4*F(1027)
       Fa16(203) = p1mup4*F(752)+p5mup4*F(1028)+p6mup4*F(1029)+p2mup4*
     &   F(1030)+p4mup4*F(1031)
       Fa16(204) = -(p1mup4*F(1032))+p2mup4*F(1033)+p5mup4*F(1034)+p6m
     &   up4*F(1035)+p4mup4*F(1036)
       Fa16(205) = -(p1mup4*F(1039))+p4mup4*F(1040)+p5mup4*F(1041)+p6m
     &   up4*F(1042)-p2mup4*F(1043)
       Fa16(206) = -(p2mup4*F(598))+p4mup4*F(1047)+p5mup4*F(1048)+p6mu
     &   p4*F(1049)-p1mup4*F(1050)
       Fa16(207) = p5mup4*F(514)+p6mup4*F(515)+p1mup4*F(524)+p2mup4*F(
     &   528)+p4mup4*F(1051)
       Fa16(208) = p4mup4*F(939)-p1mup4*F(1052)+p5mup4*F(1053)+p6mup4*
     &   F(1054)-p2mup4*F(1055)
       Return
       End
