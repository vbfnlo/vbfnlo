c************** Calling the Fa functions*************************
       subroutine HexBoxCross1Fa17(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,
     &   p3mup4,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa
     &   17)
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
       Complex*16 Fa17(209:221)
       COMMON/HexBoxCross1FaFunctions/Fa
       Fa17(209) = p4mup4*F(1056)+p5mup4*F(1057)+p6mup4*F(1058)-p1mup4
     &   *F(1059)-p2mup4*F(1060)
       Fa17(210) = p1mup3*F(755)+p5mup3*F(1061)+p6mup3*F(1062)+p2mup3*
     &   F(1063)-p4mup3*F(1064)
       Fa17(211) = p1mup3*F(1065)-p2mup3*F(1066)+p4mup3*F(1067)+p5mup3
     &   *F(1068)+p6mup3*F(1069)
       Fa17(212) = 4*mup3mup4*F(1071)-8*(p5mup4*p6mup3*F(651)+p5mup3*(
     &   p5mup4*F(650)+p4mup4*F(652)+p6mup4*F(653)-p1mup4*F(654)-p2mu
     &   p4*F(655))-p2mup3*p4mup4*F(992)+p2mup3*p2mup4*F(998)+p4mup3*
     &   (p5mup4*F(649)+p4mup4*F(1001)+p6mup4*F(1002)-p1mup4*F(1004)-
     &   p2mup4*F(1005))-p2mup3*p6mup4*F(1074)-p2mup3*p5mup4*F(1075)+
     &   p1mup4*p2mup3*F(1077)-p1mup3*(p5mup4*F(635)+p6mup4*F(1072)+p
     &   4mup4*F(1073)-p1mup4*F(1076)-p2mup4*F(1078))+p4mup4*p6mup3*F
     &   (1079)+p6mup3*p6mup4*F(1080)-p1mup4*p6mup3*F(1081)-p2mup4*p6
     &   mup3*F(1082))
       Fa17(213) = p2mup3*F(824)+p1mup3*F(825)+p5mup3*F(827)+p6mup3*F(
     &   828)-p4mup3*F(1062)
       Fa17(214) = p5mup3*F(724)+p1mup3*F(1085)+p4mup3*F(1086)+p6mup3*
     &   F(1087)+p2mup3*F(1088)
       Fa17(215) = p4mup4*F(1089)+p5mup4*F(1090)+p6mup4*F(1091)+p1mup4
     &   *F(1092)+p2mup4*F(1093)
       Fa17(216) = -(p2mup4*F(1094))+p4mup4*F(1095)+p5mup4*F(1096)+p6m
     &   up4*F(1097)-p1mup4*F(1098)
       Fa17(217) = p4mup4*F(1099)+p5mup4*F(1100)+p6mup4*F(1101)-p2mup4
     &   *F(1102)-p1mup4*F(1103)
       Fa17(218) = 4*mup3mup4*F(1105)-8*(p5mup3*p5mup4*F(1111)+p5mup3*
     &   p6mup4*F(1112)+p4mup4*p5mup3*F(1113)-p2mup4*p5mup3*F(1116)-p
     &   2mup3*(p4mup4*F(1106)+p6mup4*F(1107)+p5mup4*F(1108)+p1mup4*F
     &   (1117)+p2mup4*F(1118))+p4mup3*(p4mup4*F(1109)+p5mup4*F(1110)
     &   +p6mup4*F(1114)-p2mup4*F(1115)-p1mup4*F(1119))-p1mup4*p5mup3
     &   *F(1120)-p1mup3*p5mup4*F(1121)-p1mup3*p6mup4*F(1122)-p1mup3*
     &   p4mup4*F(1123)+p5mup4*p6mup3*F(1124)+p4mup4*p6mup3*F(1125)+p
     &   6mup3*p6mup4*F(1126)-p1mup3*p2mup4*F(1127)+p2mup4*p6mup3*F(1
     &   128)+p1mup3*p1mup4*F(1129)-p1mup4*p6mup3*F(1130))
       Fa17(219) = p4mup4*F(1133)+p1mup4*F(1134)+p5mup4*F(1135)+p6mup4
     &   *F(1136)+p2mup4*F(1137)
       Fa17(220) = -2*(p4mup4*F(1138)-p1mup4*F(1139)-p5mup4*F(1140)-p6
     &   mup4*F(1141))+p2mup4*F(1142)
       Fa17(221) = 2*(p4mup4*F(1143)+p1mup4*F(1144)-p5mup4*F(1145)+p6m
     &   up4*F(1146))+p2mup4*F(1147)
       Return
       End
