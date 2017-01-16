c************** Calling the Fa functions*************************
       subroutine HexBoxFa18(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4
     &   ,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa18)
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
       Complex*16 Fa18(205:216)
       COMMON/HexBoxFaFunctions/Fa
       Fa18(205) = -(p2mup3*F(687))+p4mup3*F(1070)+p5mup3*F(1071)-p6mu
     &   p3*F(1072)+p1mup3*F(1073)
       Fa18(206) = p6mup4*F(1075)+p4mup4*F(1076)+p5mup4*F(1077)+p1mup4
     &   *F(1078)+p2mup4*F(1079)
       Fa18(207) = p6mup4*F(1080)+p4mup4*F(1081)+p5mup4*F(1082)-p2mup4
     &   *F(1083)-p1mup4*F(1084)
       Fa18(208) = p1mup4*F(1086)-p6mup4*F(1087)-p5mup4*F(1088)+p2mup4
     &   *F(1089)-p4mup4*F(1090)
       Fa18(209) = -(p2mup4*F(1090))-p1mup4*F(1091)+p6mup4*F(1092)+p5m
     &   up4*F(1093)+p4mup4*F(1094)
       Fa18(210) = p4mup4*F(672)+p5mup4*F(673)+p6mup4*F(674)-p1mup4*F(
     &   678)-p2mup4*F(680)
       Fa18(211) = p6mup4*F(1095)+p4mup4*F(1096)+p5mup4*F(1097)-p1mup4
     &   *F(1098)-p2mup4*F(1099)
       Fa18(212) = p6mup4*F(1101)+p4mup4*F(1102)+p5mup4*F(1103)+p1mup4
     &   *F(1104)+p2mup4*F(1105)
       Fa18(213) = p1mup4*F(1106)+p2mup4*F(1107)-p6mup4*F(1108)-p5mup4
     &   *F(1109)+p4mup4*F(1110)
       Fa18(214) = p6mup4*F(1111)+p4mup4*F(1112)+p5mup4*F(1113)-p2mup4
     &   *F(1114)-p1mup4*F(1115)
       Fa18(215) = p6mup4*F(1117)+p4mup4*F(1118)+p5mup4*F(1119)-p2mup4
     &   *F(1120)-p1mup4*F(1121)
       Fa18(216) = p2mup3*F(1122)+p1mup3*F(1123)-p6mup3*F(1124)-p4mup3
     &   *F(1125)-p5mup3*F(1126)
       Return
       End
