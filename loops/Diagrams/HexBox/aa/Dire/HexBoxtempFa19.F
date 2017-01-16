c************** Calling the Fa functions*************************
       subroutine HexBoxFa19(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4
     &   ,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa19)
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
       Complex*16 Fa19(217:228)
       COMMON/HexBoxFaFunctions/Fa
       Fa19(217) = p1mup3*F(1128)-p2mup3*F(1129)+p4mup3*F(1130)+p6mup3
     &   *F(1131)-p5mup3*F(1132)
       Fa19(218) = p1mup3*F(1133)+p6mup3*F(1134)+p4mup3*F(1135)+p5mup3
     &   *F(1136)+p2mup3*F(1137)
       Fa19(219) = -(p1mup3*F(1139))-p2mup3*F(1140)+p6mup3*F(1141)+p4m
     &   up3*F(1142)+p5mup3*F(1143)
       Fa19(220) = p6mup4*F(1145)+p4mup4*F(1146)+p5mup4*F(1147)-p2mup4
     &   *F(1148)-p1mup4*F(1149)
       Fa19(221) = p6mup4*F(1152)+p4mup4*F(1153)+p5mup4*F(1154)+p1mup4
     &   *F(1155)+p2mup4*F(1156)
       Fa19(222) = -(p1mup4*F(1157))+p2mup4*F(1158)+p6mup4*F(1159)+p5m
     &   up4*F(1160)+p4mup4*F(1161)
       Fa19(223) = p1mup3*F(1162)-(p2mup3+p6mup3)*F(1163)-p4mup3*F(116
     &   4)-p5mup3*F(1165)
       Fa19(224) = p1mup4*F(1167)-p6mup4*F(1168)-p5mup4*F(1169)+p2mup4
     &   *F(1170)-p4mup4*F(1171)
       Fa19(225) = p1mup4*F(1172)+p2mup4*F(1173)+p6mup4*F(1174)+p4mup4
     &   *F(1175)+p5mup4*F(1176)
       Fa19(226) = p1mup4*F(1177)+p2mup4*F(1178)+p6mup4*F(1179)+p4mup4
     &   *F(1180)+p5mup4*F(1181)
       Fa19(227) = p5mup3*F(1137)-p1mup3*F(1182)-p2mup3*F(1183)+p6mup3
     &   *F(1184)+p4mup3*F(1185)
       Fa19(228) = p1mup3*F(1187)-p2mup3*F(1188)+p6mup3*F(1189)+p4mup3
     &   *F(1190)+p5mup3*F(1191)
       Return
       End
