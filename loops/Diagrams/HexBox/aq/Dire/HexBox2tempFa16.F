c************** Calling the Fa functions*************************
       subroutine HexBox2Fa16(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup
     &   6,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa16)
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
       Complex*16 Fa16(286:304)
       COMMON/HexBox2FaFunctions/Fa
       Fa16(286) = p2mup3*F(988)+p1mup3*F(989)-p6mup3*F(990)-p4mup3*F(
     &   991)-p5mup3*F(992)
       Fa16(287) = p2mup3*F(993)+p1mup3*F(994)
       Fa16(288) = -(p1mup3*F(997))+p2mup3*F(998)+p6mup3*F(999)+p4mup3
     &   *F(1000)-p5mup3*F(1001)
       Fa16(289) = p5mup3*p5mup6*F(1022)-8*(-(p2mup3*p2mup6*F(1004))+p
     &   1mup3*(p1mup6*F(1003)-p2mup6*F(1005)+p5mup6*F(1006)+p6mup6*F
     &   (1007)+p4mup6*F(1008))+p2mup6*(p4mup3*F(1009)+p5mup3*F(1010)
     &   +p6mup3*F(1011))-mup3mup6*F(1015)+p2mup3*(p5mup6*F(1012)+p6m
     &   up6*F(1013)+p4mup6*F(1014)-p1mup6*F(1016))+p1mup6*p5mup3*F(1
     &   017)-p1mup6*p6mup3*F(1018)-p1mup6*p4mup3*F(1019)-p4mup6*p5mu
     &   p3*F(1020)-p4mup3*p5mup6*F(1020)+p4mup6*p6mup3*F(1021)+p4mup
     &   3*p6mup6*F(1021)-p5mup6*p6mup3*F(1023)-p5mup3*p6mup6*F(1023)
     &   +p6mup3*p6mup6*F(1024)+p4mup3*p4mup6*F(1025))
       Fa16(290) = p2mup6*F(1026)+p5mup6*F(1029)+4*(p1mup6*F(1027)+p4m
     &   up6*F(1028)+p6mup6*F(1030))
       Fa16(291) = p2mup6*F(1031)+p4mup6*F(1032)+p5mup6*F(1033)+p6mup6
     &   *F(1034)-p1mup6*F(1035)
       Fa16(292) = p6mup3*F(915)+p1mup3*F(1036)+p2mup3*F(1037)+p4mup3*
     &   F(1038)+p5mup3*F(1039)
       Fa16(293) = -4*mup3mup6*F(1040)+8*(-(p2mup6*p4mup3*F(1042))-p2m
     &   up6*p5mup3*F(1043)-p2mup6*p6mup3*F(1044)+p4mup6*p5mup3*F(104
     &   7)+p4mup3*p5mup6*F(1047)+p4mup6*p6mup3*F(1048)+p4mup3*p6mup6
     &   *F(1048)+p5mup3*p5mup6*F(1049)+p5mup6*p6mup3*F(1050)+p5mup3*
     &   p6mup6*F(1050)+p6mup3*p6mup6*F(1051)+p1mup6*(p2mup3*F(1041)-
     &   p1mup3*F(1052)+p5mup3*F(1053)+p6mup3*F(1054)+p4mup3*F(1055))
     &   +p1mup3*(p2mup6*F(1041)-p5mup6*F(1045)-p6mup6*F(1046)-p4mup6
     &   *F(1056))+p4mup3*p4mup6*F(1057)+p2mup3*p2mup6*F(1058)-p2mup3
     &   *p5mup6*F(1059)-p2mup3*p6mup6*F(1060)-p2mup3*p4mup6*F(1061))
       Fa16(294) = 4*mup3mup6*F(1062)-8*(-(p5mup6*p6mup3*F(907))-p6mup
     &   3*p6mup6*F(908)-p2mup3*p5mup6*F(909)-p2mup3*p6mup6*F(910)-p1
     &   mup6*p2mup3*F(941)-p2mup3*p2mup6*F(945)-p1mup3*p2mup6*F(1063
     &   )+p4mup6*(p5mup3*F(903)+p6mup3*F(904)+p4mup3*F(905)-p2mup3*F
     &   (911)-p1mup3*F(1064))-p1mup3*p5mup6*F(1065)-p1mup3*p6mup6*F(
     &   1066)-p1mup3*p1mup6*F(1067)+p1mup6*p6mup3*F(1069)+p4mup3*(p5
     &   mup6*F(903)+p6mup6*F(904)+p1mup6*F(1070)+p2mup6*F(1071))-p5m
     &   up3*(p5mup6*F(906)+p6mup6*F(907)-p1mup6*F(1068)-p2mup6*F(107
     &   2))+p2mup6*p6mup3*F(1073))
       Fa16(295) = p2mup6*F(1074)+p4mup6*F(1075)-p5mup6*F(1076)-p6mup6
     &   *F(1077)+p1mup6*F(1078)
       Fa16(296) = p6mup3*F(1081)+p4mup3*F(1082)+p1mup3*F(1083)+p5mup3
     &   *F(1084)+p2mup3*F(1085)
       Fa16(297) = p2mup3*(p2mup6*F(1086)+p5mup6*F(1087)+p6mup6*F(1088
     &   )+p4mup6*F(1089)+p1mup6*F(1090))+p2mup6*p4mup3*F(1091)+p2mup
     &   6*p5mup3*F(1092)+p2mup6*p6mup3*F(1093)-p4mup6*p5mup3*F(1094)
     &   -p4mup3*p5mup6*F(1094)-p4mup6*p6mup3*F(1095)-p4mup3*p6mup6*F
     &   (1095)-p5mup3*p5mup6*F(1096)-p5mup6*p6mup3*F(1097)-p5mup3*p6
     &   mup6*F(1097)-p6mup3*p6mup6*F(1098)-p1mup6*p5mup3*F(1099)-p1m
     &   up6*p6mup3*F(1100)-p1mup6*p4mup3*F(1101)-p4mup3*p4mup6*F(110
     &   2)-p1mup3*p2mup6*F(1103)+p1mup3*p5mup6*F(1104)+p1mup3*p6mup6
     &   *F(1105)+p1mup3*p1mup6*F(1106)+p1mup3*p4mup6*F(1107)+mup3mup
     &   6*F(1108)
       Fa16(298) = p1mup3*(-(p2mup6*F(1109))+p5mup6*F(1110)+p6mup6*F(1
     &   111)+p4mup6*F(1112))+p4mup6*p5mup3*F(1116)+p4mup3*p5mup6*F(1
     &   116)+p4mup6*p6mup3*F(1117)+p4mup3*p6mup6*F(1117)+p4mup3*p4mu
     &   p6*F(1118)+p5mup3*p5mup6*F(1119)+p5mup6*p6mup3*F(1120)+p5mup
     &   3*p6mup6*F(1120)-p6mup3*p6mup6*F(1121)+p1mup6*(p4mup3*F(1095
     &   )+p5mup3*F(1097)+p6mup3*F(1098)+p1mup3*F(1122)-p2mup3*F(1123
     &   ))+p2mup6*(p4mup3*F(1113)+p5mup3*F(1114)+p6mup3*F(1115)-p2mu
     &   p3*F(1124))+p2mup3*p5mup6*F(1125)+p2mup3*p6mup6*F(1126)+p2mu
     &   p3*p4mup6*F(1127)-mup3mup6*F(1128)
       Fa16(299) = -(p2mup3*F(989))-p1mup3*F(1132)+p6mup3*F(1133)+p4mu
     &   p3*F(1134)+p5mup3*F(1135)
       Fa16(300) = p2mup3*F(990)+p1mup3*F(1136)+p4mup3*F(1137)+p6mup3*
     &   F(1138)-p5mup3*F(1139)
       Fa16(301) = p1mup3*F(1141)+p6mup3*F(1143)-4*(p4mup3*F(1140)+p5m
     &   up3*F(1142)+p2mup3*F(1144))
       Fa16(302) = p1mup3*F(1145)+p4mup3*F(1148)-4*(p2mup3*F(1146)+p6m
     &   up3*F(1147)-p5mup3*F(1149))
       Fa16(303) = p6mup3*F(1150)+p4mup3*F(1151)+p1mup3*F(1152)+p5mup3
     &   *F(1153)+p2mup3*F(1154)
       Fa16(304) = p2mup3*F(991)+p1mup3*F(1159)+p6mup3*F(1160)+p4mup3*
     &   F(1161)+p5mup3*F(1162)
       Return
       End
