c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineNoAbe1g21(p1mup2,p1mup3,p1mup4
     -   ,p1mup6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3m
     -   up6,p4mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,
     -   p6mup2,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3m
     -   up4,mup3mup6,mup4mup6,Fa21)
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
       Complex*16 Fa21(281:294)
       COMMON/test_FaFunctionsNoAbe1g/Fa
       Fa21(281) = p4mup2*(p4mup6*F(1046)+4*((p1mup6+p2mup6)*F(1047)+p
     -   3mup6*F(1048))+p5mup6*F(1049))+p3mup2*(-(p4mup6*F(1045))+4*(
     -   p5mup6*F(1050)-(p1mup6+p2mup6)*F(1051))+p3mup6*F(1052))
       Fa21(282) = 2*mup4mup6*F(1054)-4*(p1mup4*p5mup6*F(1057)+p2mup4*
     -   p5mup6*F(1057)+p3mup4*(p5mup6*F(1055)+(p1mup6+p2mup6)*F(1058
     -   )+p3mup6*F(1059)+p4mup6*F(1060))+p4mup4*(p5mup6*F(1056)+(p1m
     -   up6+p2mup6)*F(1061)+p3mup6*F(1062)+p4mup6*F(1063))+p1mup4*p1
     -   mup6*F(1064)+p1mup6*p2mup4*F(1064)+p1mup4*p2mup6*F(1064)+p2m
     -   up4*p2mup6*F(1064)+p1mup4*p3mup6*F(1065)+p2mup4*p3mup6*F(106
     -   5)+p1mup4*p4mup6*F(1066)+p2mup4*p4mup6*F(1066)+p5mup4*p5mup6
     -   *F(1067)+p1mup6*p5mup4*F(1068)+p2mup6*p5mup4*F(1068)+p3mup6*
     -   p5mup4*F(1069)+p4mup6*p5mup4*F(1070))
       Fa21(283) = mup3mup6*F(1071)+4*p4mup3*p5mup6*F(1075)+8*(p4mup3*
     -   (p4mup6*F(1072)+p3mup6*F(1073)+(p1mup6+p2mup6)*F(1074))-p4mu
     -   p6*((p1mup3+p2mup3)*F(1076)-p3mup3*F(1077)+p5mup3*F(1078)))
       Fa21(284) = -12*p3mup6*p4mup2*s12*F(1073)+(p1mup2+p2mup2)*(p4mu
     -   p6*F(1079)+4*p5mup6*F(1080))+(p1mup2+p2mup2)*p3mup6*F(1082)-
     -   4*(3*(p1mup6*p4mup2*s12*F(1074)+p2mup6*p4mup2*s12*F(1074))+p
     -   1mup6*(p1mup2+p2mup2)*F(1081)+(p1mup2+p2mup2)*p2mup6*F(1081)
     -   +p4mup2*s12*(3*p4mup6*F(1083)+p5mup6*F(1084)))
       Fa21(285) = p1mup6*p3mup4*F(1074)+p2mup6*p3mup4*F(1074)+p3mup4*
     -   p4mup6*F(1083)+mup4mup6*F(1086)+p3mup4*p5mup6*F(1087)+p3mup4
     -   *p3mup6*F(1088)-2*(p1mup4*p3mup6*F(1089)+p2mup4*p3mup6*F(108
     -   9)-p3mup6*p5mup4*F(1090)-p3mup6*p4mup4*F(1091))
       Fa21(286) = p5mup3*p5mup6*F(1067)+p1mup6*p5mup3*F(1068)+p2mup6*
     -   p5mup3*F(1068)-mup3mup6*F(1092)+p3mup3*p3mup6*F(1093)+p3mup3
     -   *p4mup6*F(1094)+p1mup3*(p5mup6*F(1057)+(p1mup6+p2mup6)*F(106
     -   4)+p4mup6*F(1095)+p3mup6*F(1096))+p2mup3*(p5mup6*F(1057)+(p1
     -   mup6+p2mup6)*F(1064)+p4mup6*F(1095)+p3mup6*F(1096))+p1mup6*p
     -   3mup3*F(1097)+p2mup6*p3mup3*F(1097)+p3mup3*p5mup6*F(1098)+p4
     -   mup3*p4mup6*F(1099)+p3mup6*p4mup3*F(1100)+p1mup6*p4mup3*F(11
     -   01)+p2mup6*p4mup3*F(1101)+p4mup3*p5mup6*F(1102)+p3mup6*p5mup
     -   3*F(1103)+p4mup6*p5mup3*F(1104)
       Fa21(287) = 4*p3mup2*s12*(3*(p3mup6*F(1073)+(p1mup6+p2mup6)*F(1
     -   074)+p4mup6*F(1083))+p5mup6*F(1084))+(p1mup2+p2mup2)*(p4mup6
     -   *F(1105)+p5mup6*F(1106)+4*(p1mup6+p2mup6)*F(1107)+p3mup6*F(1
     -   108))
       Fa21(288) = mup4mup6*(p4mup3*F(1110)+4*(p5mup3*F(1111)+(p1mup3+
     -   p2mup3)*F(1112))+p3mup3*F(1113))+8*mup3mup6*p5mup4*F(1114)
       Fa21(289) = mup3mup6*(p3mup4*F(1115)+p4mup4*F(1116))
       Fa21(290) = p5mup6*(mup3mup4*F(1123)+2*(-(p3mup4*p4mup3*F(1122)
     -   )+p4mup3*p4mup4*F(1125)+p1mup4*p4mup3*F(1126)+p2mup4*p4mup3*
     -   F(1126)+p1mup4*p5mup3*F(1127)+p2mup4*p5mup3*F(1127)+p1mup3*(
     -   p3mup4*F(1118)+p4mup4*F(1120)+(p1mup4+p2mup4)*F(1124)+p5mup4
     -   *F(1127))+p2mup3*(p3mup4*F(1118)+p4mup4*F(1120)+(p1mup4+p2mu
     -   p4)*F(1124)+p5mup4*F(1127))+p3mup3*(p3mup4*F(1117)+(p1mup4+p
     -   2mup4)*F(1119)+p4mup4*F(1121)+p5mup4*F(1128))+p3mup4*p5mup3*
     -   F(1129)+p4mup4*p5mup3*F(1130)+p4mup3*p5mup4*F(1131)+p5mup3*p
     -   5mup4*F(1132)))
       Fa21(291) = mup3mup6*(p1mup4+p2mup4)*F(1133)+p5mup4*(p1mup3*p4m
     -   up6*F(1134)+p2mup3*p4mup6*F(1134)+p4mup3*p4mup6*F(1135)+p3mu
     -   p3*p4mup6*F(1136)+p1mup6*p4mup3*F(1137)+p2mup6*p4mup3*F(1137
     -   )+p3mup6*p4mup3*F(1138)+p4mup6*p5mup3*F(1139)+p2mup6*p5mup3*
     -   F(1140)+p3mup6*p5mup3*F(1141))
       Fa21(292) = -4*(mup3mup4*p3mup6*F(1142)+mup3mup4*p4mup6*F(1143)
     -   +mup3mup4*p1mup6*F(1144)+mup3mup4*p2mup6*F(1144))+8*(p4mup4*
     -   p4mup6*p5mup3*F(1135)+p2mup6*p4mup4*p5mup3*F(1137)+p1mup6*p5
     -   mup3*p5mup4*F(1140)+p1mup3*p1mup6*p3mup4*F(1146)+p2mup3*p2mu
     -   p6*p3mup4*F(1146)+p2mup3*p3mup4*p3mup6*F(1148)-p2mup6*p3mup3
     -   *p4mup4*F(1149)+p2mup3*p3mup4*p4mup6*F(1151)+p2mup6*p3mup4*p
     -   4mup3*F(1152)-p3mup3*p4mup4*p4mup6*F(1153)+p3mup4*p3mup6*p4m
     -   up3*F(1154)+p2mup3*p3mup6*p4mup4*F(1155)+p2mup6*p3mup3*p3mup
     -   4*F(1159)+p3mup3*p3mup4*p4mup6*F(1160)+p3mup3*p3mup4*p3mup6*
     -   F(1162)+p1mup3*p1mup6*p4mup4*F(1163)+p2mup3*p2mup6*p4mup4*F(
     -   1163)+p2mup3*p4mup4*p4mup6*F(1164)+p2mup6*p3mup4*p5mup3*F(11
     -   66)+p2mup4*(p4mup6*p5mup3*F(1134)+p1mup3*p1mup6*F(1145)+p2mu
     -   p3*p2mup6*F(1145)+p2mup3*p3mup6*F(1147)+p2mup3*p4mup6*F(1150
     -   )+p3mup3*p4mup6*F(1156)-p3mup6*p4mup3*F(1157)+p2mup6*p3mup3*
     -   F(1158)+p3mup3*p3mup6*F(1161)+p2mup6*p4mup3*F(1163)+p2mup6*p
     -   5mup3*F(1165)+p3mup6*p5mup3*F(1167))+p1mup4*(p1mup3*p1mup6*F
     -   (1145)+p2mup3*p2mup6*F(1145)+p2mup3*p3mup6*F(1147)+p4mup6*(p
     -   5mup3*F(1134)+p2mup3*F(1150)+p3mup3*F(1156))-p3mup6*p4mup3*F
     -   (1157)+p2mup6*p3mup3*F(1158)+p3mup3*p3mup6*F(1161)+p2mup6*p4
     -   mup3*F(1163)+p2mup6*p5mup3*F(1165)+p3mup6*p5mup3*F(1167))+p3
     -   mup4*p3mup6*p5mup3*F(1168)+p3mup4*p4mup6*p5mup3*F(1169)+p3mu
     -   p6*p4mup4*p5mup3*F(1170)+p2mup6*p4mup3*p4mup4*F(1171)+p3mup6
     -   *p4mup3*p4mup4*F(1172))
       Fa21(293) = p1mup3*(p2mup6*p3mup4*F(1146)+p3mup4*p3mup6*F(1148)
     -   +p1mup4*(p2mup6*F(1145)+p3mup6*F(1147)+p4mup6*F(1150))+p2mup
     -   4*(p2mup6*F(1145)+p3mup6*F(1147)+p4mup6*F(1150))+p3mup4*p4mu
     -   p6*F(1151)+p3mup6*p4mup4*F(1155)+p2mup6*p4mup4*F(1163)+p4mup
     -   4*p4mup6*F(1164))+p1mup6*(p2mup3*p2mup4*F(1145)+p2mup3*p3mup
     -   4*F(1146)+p3mup4*p4mup3*F(1152)+p2mup4*p3mup3*F(1158)+p3mup3
     -   *p3mup4*F(1159)+p2mup4*p4mup3*F(1163)+p2mup4*p5mup3*F(1165)+
     -   p1mup4*(p2mup3*F(1145)+p3mup3*F(1158)+p4mup3*F(1163)+p5mup3*
     -   F(1165))+p3mup4*p5mup3*F(1166)+p4mup4*(p5mup3*F(1137)-p3mup3
     -   *F(1149)+p2mup3*F(1163)+p4mup3*F(1171)))-p3mup3*p3mup6*p4mup
     -   4*F(1173)+p4mup3*p4mup6*((p1mup4+p2mup4)*F(1164)+p3mup4*F(11
     -   74)+p4mup4*F(1175))
       Fa21(294) = p5mup4*(p1mup6*p2mup3*F(1165)+p2mup3*p2mup6*F(1165)
     -   +p1mup3*(p1mup6+p2mup6)*F(1165)+p1mup3*p3mup6*F(1167)+p2mup3
     -   *p3mup6*F(1167)+p3mup3*p3mup6*F(1176)+(p1mup6+p2mup6)*p3mup3
     -   *F(1177))
       Return
       End
