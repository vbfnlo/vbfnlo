c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineNoAbe1g22(p1mup2,p1mup3,p1mup4
     -   ,p1mup6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3m
     -   up6,p4mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,
     -   p6mup2,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3m
     -   up4,mup3mup6,mup4mup6,Fa22)
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
       Complex*16 Fa22(295:308)
       COMMON/test_FaFunctionsNoAbe1g/Fa
       Fa22(295) = p1mup2*(mup4mup6*F(1178)+8*((p1mup6+p2mup6)*(p4mup4
     -   *F(1013)+(p1mup4+p2mup4)*F(1023))+p3mup4*(p4mup6*F(1012)+p5m
     -   up6*F(1180)+(p1mup6+p2mup6)*F(1182)+p3mup6*F(1183)))+p4mup4*
     -   (p5mup6*F(1181)+8*(p4mup6*F(1022)+p3mup6*F(1184)))+p1mup4*p4
     -   mup6*F(1185)+p2mup4*p4mup6*F(1185)+p1mup4*p3mup6*F(1187)+p2m
     -   up4*p3mup6*F(1187)+p5mup4*p5mup6*F(1188)+(p1mup6+p2mup6)*p5m
     -   up4*F(1189)-p3mup6*p5mup4*F(1190)-p4mup6*p5mup4*F(1191)+(p1m
     -   up4+p2mup4)*p5mup6*F(1192))+p2mup2*(mup4mup6*F(1178)+8*((p1m
     -   up6+p2mup6)*(p4mup4*F(1013)+(p1mup4+p2mup4)*F(1023))+p3mup4*
     -   (p4mup6*F(1012)+p5mup6*F(1180)+(p1mup6+p2mup6)*F(1182)+p3mup
     -   6*F(1183)))+p4mup4*(p5mup6*F(1181)+8*(p4mup6*F(1022)+p3mup6*
     -   F(1184)))+p1mup4*p4mup6*F(1185)+p2mup4*p4mup6*F(1185)+p1mup4
     -   *p3mup6*F(1187)+p2mup4*p3mup6*F(1187)+p5mup4*p5mup6*F(1188)+
     -   (p1mup6+p2mup6)*p5mup4*F(1189)-p3mup6*p5mup4*F(1190)-p4mup6*
     -   p5mup4*F(1191)+(p1mup4+p2mup4)*p5mup6*F(1192))+s12*(mup4mup6
     -   *p3mup2*F(1186)+8*(p3mup2*(p1mup4*p5mup6*F(1057)+p2mup4*p5mu
     -   p6*F(1057)+p3mup4*(p5mup6*F(1055)+(p1mup6+p2mup6)*F(1058)+p3
     -   mup6*F(1059)+p4mup6*F(1060))+p4mup4*(p5mup6*F(1056)+(p1mup6+
     -   p2mup6)*F(1061)+p3mup6*F(1062)+p4mup6*F(1063))+p1mup4*p1mup6
     -   *F(1064)+p1mup6*p2mup4*F(1064)+p1mup4*p2mup6*F(1064)+p2mup4*
     -   p2mup6*F(1064)+p1mup4*p3mup6*F(1065)+p2mup4*p3mup6*F(1065)+p
     -   1mup4*p4mup6*F(1066)+p2mup4*p4mup6*F(1066)+p5mup4*p5mup6*F(1
     -   067)+p1mup6*p5mup4*F(1068)+p2mup6*p5mup4*F(1068)+p3mup6*p5mu
     -   p4*F(1069)+p4mup6*p5mup4*F(1070))-p4mup2*(p1mup6*p3mup4*F(10
     -   74)+p2mup6*p3mup4*F(1074)+p3mup4*p4mup6*F(1083)+p3mup4*p5mup
     -   6*F(1087)+p3mup4*p3mup6*F(1088)-mup4mup6*F(1179)-2*(p1mup4*p
     -   3mup6*F(1089)+p2mup4*p3mup6*F(1089)+p3mup6*p5mup4*F(1193)+p3
     -   mup6*p4mup4*F(1194)))))
       Fa22(296) = -(p1mup2*(-4*mup3mup6*F(1195)+p4mup3*(8*(p4mup6*F(1
     -   034)-(p1mup6+p2mup6)*F(1199))+p5mup6*F(1200))-8*(p1mup6*p3mu
     -   p3*F(1032)+p2mup6*p3mup3*F(1032)-p1mup6*(p1mup3+p2mup3)*F(10
     -   38)-(p1mup3+p2mup3)*p2mup6*F(1038)-(p1mup3+p2mup3)*p5mup6*F(
     -   1039)+p5mup3*p5mup6*F(1040)-(p1mup3+p2mup3)*p4mup6*F(1198)-p
     -   3mup3*(p4mup6*F(1196)+p5mup6*F(1204))-p3mup6*(p4mup3*F(1026)
     -   +(p1mup3+p2mup3)*F(1202)+p3mup3*F(1203)-p5mup3*F(1205))+p4mu
     -   p6*p5mup3*F(1206)+(p1mup6+p2mup6)*p5mup3*F(1207))))-p2mup2*(
     -   -4*mup3mup6*F(1195)+p4mup3*(8*(p4mup6*F(1034)-(p1mup6+p2mup6
     -   )*F(1199))+p5mup6*F(1200))-8*(p1mup6*p3mup3*F(1032)+p2mup6*p
     -   3mup3*F(1032)-p1mup6*(p1mup3+p2mup3)*F(1038)-(p1mup3+p2mup3)
     -   *p2mup6*F(1038)-(p1mup3+p2mup3)*p5mup6*F(1039)+p5mup3*p5mup6
     -   *F(1040)-(p1mup3+p2mup3)*p4mup6*F(1198)-p3mup3*(p4mup6*F(119
     -   6)+p5mup6*F(1204))-p3mup6*(p4mup3*F(1026)+(p1mup3+p2mup3)*F(
     -   1202)+p3mup3*F(1203)-p5mup3*F(1205))+p4mup6*p5mup3*F(1206)+(
     -   p1mup6+p2mup6)*p5mup3*F(1207)))+4*s12*(2*(p4mup2*p5mup3*p5mu
     -   p6*F(1067)+p1mup6*p4mup2*p5mup3*F(1068)+p2mup6*p4mup2*p5mup3
     -   *F(1068)-p3mup2*p4mup3*p5mup6*F(1075)+p3mup3*p3mup6*p4mup2*F
     -   (1093)+p3mup3*p4mup2*p4mup6*F(1094)+p1mup3*(2*p3mup2*p4mup6*
     -   F(1076)+p4mup2*(p5mup6*F(1057)+(p1mup6+p2mup6)*F(1064)+p4mup
     -   6*F(1095)+p3mup6*F(1096)))+p2mup3*(2*p3mup2*p4mup6*F(1076)+p
     -   4mup2*(p5mup6*F(1057)+(p1mup6+p2mup6)*F(1064)+p4mup6*F(1095)
     -   +p3mup6*F(1096)))+p1mup6*p3mup3*p4mup2*F(1097)+p2mup6*p3mup3
     -   *p4mup2*F(1097)+p3mup3*p4mup2*p5mup6*F(1098)+p4mup2*p4mup3*p
     -   4mup6*F(1099)+p3mup6*p4mup2*p4mup3*F(1100)+p1mup6*p4mup2*p4m
     -   up3*F(1101)+p2mup6*p4mup2*p4mup3*F(1101)+p4mup2*p4mup3*p5mup
     -   6*F(1102)+p3mup6*p4mup2*p5mup3*F(1103)+p4mup2*p4mup6*p5mup3*
     -   F(1104)+mup3mup6*p4mup2*F(1197))+mup3mup6*p3mup2*F(1201)-4*(
     -   p3mup2*p4mup3*(p4mup6*F(1072)+p3mup6*F(1073))+p1mup6*p3mup2*
     -   p4mup3*F(1074)+p2mup6*p3mup2*p4mup3*F(1074)-p3mup2*p4mup6*p5
     -   mup3*F(1078)-p3mup2*p3mup3*p4mup6*F(1208)))
       Fa22(297) = p4mup4*F(1209)+p3mup4*F(1210)+p5mup4*F(1211)+(p1mup
     -   4+p2mup4)*F(1212)
       Fa22(298) = p4mup3*F(1213)+p5mup3*F(1214)+(p1mup3+p2mup3)*F(121
     -   5)+p3mup3*F(1216)
       Fa22(299) = (p1mup2+p2mup2)*F(1217)+s12*(p4mup2*F(1218)+p3mup2*
     -   F(1219))
       Fa22(300) = mup4mup6*(p3mup3*p4mup2*F(1221)+p4mup2*p4mup3*F(122
     -   2)+8*(p3mup2*p4mup3*F(1220)-p3mup2*p5mup3*F(1223))+p4mup2*p5
     -   mup3*F(1224))
       Fa22(301) = mup4mup6*(p2mup3*p4mup2*F(1225)-8*(p2mup3*p3mup2*F(
     -   1226)+p3mup2*p3mup3*F(1227)))
       Fa22(302) = mup4mup6*p1mup3*(p4mup2*F(1225)-8*p3mup2*F(1226))
       Fa22(303) = mup3mup6*(8*p3mup2*(-2*(p3mup4*F(1228)-p4mup4*F(122
     -   9))+p5mup4*F(1231))+p4mup2*(p3mup4*F(1230)+p5mup4*F(1232)+p4
     -   mup4*F(1233)))
       Fa22(304) = mup3mup6*(p1mup4+p2mup4)*(p4mup2*F(1234)+p3mup2*F(1
     -   235))
       Fa22(305) = mup3mup4*(p4mup2*(p4mup6*F(1236)+p3mup6*F(1237)+p5m
     -   up6*F(1238)+(p1mup6+p2mup6)*F(1239))+8*p3mup2*(p3mup6*F(1240
     -   )+p4mup6*F(1241)+p5mup6*F(1242)-(p1mup6+p2mup6)*F(1243)))
       Fa22(306) = p3mup2*(p2mup3*p2mup4*p2mup6*F(1244)+p2mup3*p2mup6*
     -   p3mup4*F(1245)+p2mup3*p2mup6*p4mup4*F(1246)+p1mup3*p1mup6*((
     -   p1mup4+p2mup4)*F(1244)+p3mup4*F(1245)+p4mup4*F(1246))+p2mup3
     -   *p2mup4*p3mup6*F(1247)+p2mup3*p3mup4*p3mup6*F(1248)+p2mup3*p
     -   3mup6*p4mup4*F(1249)+p2mup3*p2mup4*p4mup6*F(1250)-p2mup4*p2m
     -   up6*p4mup3*F(1251)+p2mup3*p3mup4*p4mup6*F(1252)+p2mup3*p4mup
     -   4*p4mup6*F(1253)+p2mup3*p2mup4*p5mup6*F(1254)+p2mup3*p3mup4*
     -   p5mup6*F(1255)+p2mup3*p4mup4*p5mup6*F(1256)-p2mup4*p3mup6*p4
     -   mup3*F(1257)-p2mup4*p4mup3*p5mup6*F(1258)-p2mup4*p4mup3*p4mu
     -   p6*F(1259)+p2mup4*p2mup6*p3mup3*F(1260)+p2mup6*p3mup3*p3mup4
     -   *F(1261)+p2mup6*p3mup4*p4mup3*F(1262)+p2mup4*p3mup3*p4mup6*F
     -   (1263)+p3mup3*p3mup4*p4mup6*F(1264)-p2mup4*p3mup3*p5mup6*F(1
     -   265)-p3mup3*p3mup4*p5mup6*F(1266)-p3mup3*p4mup4*p5mup6*F(126
     -   7)+p2mup4*p3mup3*p3mup6*F(1268)+p3mup3*p3mup4*p3mup6*F(1269)
     -   +p3mup4*p3mup6*p4mup3*F(1270)+p3mup4*p4mup3*p5mup6*F(1271)+p
     -   3mup4*p4mup3*p4mup6*F(1272)-p2mup6*p4mup4*p5mup3*F(1273)-p3m
     -   up6*p4mup4*p5mup3*F(1274)-p4mup4*p4mup6*p5mup3*F(1275)-p2mup
     -   4*p4mup6*p5mup3*F(1276)-p3mup4*p4mup6*p5mup3*F(1277)+p4mup3*
     -   p4mup4*p5mup6*F(1278)+p2mup6*p3mup3*p4mup4*F(1279)-p2mup6*p4
     -   mup3*p4mup4*F(1280)+p3mup6*p4mup3*p4mup4*F(1281)+p4mup3*p4mu
     -   p4*p4mup6*F(1282)+p2mup4*p5mup3*p5mup6*F(1283)+p3mup4*p5mup3
     -   *p5mup6*F(1284)-p2mup4*p2mup6*p5mup3*F(1285)-p2mup6*p3mup4*p
     -   5mup3*F(1286)+p3mup3*p4mup4*p4mup6*F(1287)+p3mup3*p3mup6*p4m
     -   up4*F(1288)+p2mup4*p3mup6*p5mup3*F(1289)+p1mup4*(-(p2mup6*p4
     -   mup3*F(1251))+p2mup3*(p2mup6*F(1244)+p3mup6*F(1247)+p4mup6*F
     -   (1250)+p5mup6*F(1254))-p3mup6*p4mup3*F(1257)-p4mup3*p5mup6*F
     -   (1258)-p4mup3*p4mup6*F(1259)+p2mup6*p3mup3*F(1260)+p3mup3*p4
     -   mup6*F(1263)-p3mup3*p5mup6*F(1265)+p3mup3*p3mup6*F(1268)-p4m
     -   up6*p5mup3*F(1276)+p5mup3*p5mup6*F(1283)-p2mup6*p5mup3*F(128
     -   5)+p3mup6*p5mup3*F(1289))+p3mup4*p3mup6*p5mup3*F(1290))
       Fa22(307) = p3mup2*(p1mup3*(p2mup6*p3mup4*F(1245)+p2mup6*p4mup4
     -   *F(1246)+p3mup4*p3mup6*F(1248)+p3mup6*p4mup4*F(1249)+p3mup4*
     -   p4mup6*F(1252)+p4mup4*p4mup6*F(1253)+p2mup4*(p2mup6*F(1244)+
     -   p3mup6*F(1247)+p4mup6*F(1250)+p5mup6*F(1254))+p3mup4*p5mup6*
     -   F(1255)+p4mup4*p5mup6*F(1256))+p1mup4*(p1mup3*(p2mup6*F(1244
     -   )+p3mup6*F(1247)+p4mup6*F(1250)+p5mup6*F(1254))+p1mup6*(p2mu
     -   p3*F(1244)-p4mup3*F(1251)+p3mup3*F(1260)-p5mup3*F(1285)))+p1
     -   mup6*(p2mup3*(p2mup4*F(1244)+p3mup4*F(1245)+p4mup4*F(1246))-
     -   p2mup4*p4mup3*F(1251)+p2mup4*p3mup3*F(1260)+p3mup3*p3mup4*F(
     -   1261)+p3mup4*p4mup3*F(1262)-p4mup4*p5mup3*F(1273)+p3mup3*p4m
     -   up4*F(1279)-p4mup3*p4mup4*F(1280)-p2mup4*p5mup3*F(1285)-p3mu
     -   p4*p5mup3*F(1286))-p4mup4*p5mup3*p5mup6*F(1291))
       Fa22(308) = p1mup3*p2mup6*p4mup2*((p1mup4+p2mup4)*F(1292)-p3mup
     -   4*F(1293))
       Return
       End
