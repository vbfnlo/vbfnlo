c************** Calling the Fa functions*************************
       subroutine HexBoxCross2Fa16(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,
     &   p3mup6,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa
     &   16)
       IMPLICIT NONE
      Complex*16   p1mup3, p1mup6, p2mup3, p2mup6, p3mup3, p3mup6, 
     -          p4mup3, p4mup6, p5mup3, p5mup6, p6mup3, p6mup6
       Complex*16   mup3mup6
        common/HexBoxCross2FhlFunctions/F
       COMMON/HexBoxCross2Invariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12
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
      Complex*16 Fa(392),F(1853)
      Real*8 P(1355) 
        COMMON/HexBoxCross2PFunctions/P
       Complex*16 Fa16(286:304)
       COMMON/HexBoxCross2FaFunctions/Fa
       Fa16(286) = p1mup6*F(1143)+p6mup6*F(1144)+p2mup6*F(1145)+p4mup6
     &   *F(1146)+p5mup6*F(1147)
       Fa16(287) = p1mup6*F(1149)+p6mup6*F(1150)+p2mup6*F(1151)+p4mup6
     &   *F(1152)+p5mup6*F(1153)
       Fa16(288) = -(p4mup3*p6mup6*F(1159))+p2mup6*p4mup3*F(1162)+p2mu
     &   p3*(p6mup6*F(1157)+p5mup6*F(1163)+p4mup6*F(1164)+p1mup6*F(11
     &   65)+p2mup6*F(1166))+p1mup6*p4mup3*F(1168)+p4mup3*p4mup6*F(11
     &   69)+p4mup3*p5mup6*F(1170)-p1mup3*p6mup6*F(1171)+p5mup3*p6mup
     &   6*F(1172)-p6mup3*(p6mup6*F(1158)+p5mup6*F(1160)-p2mup6*F(116
     &   1)-p4mup6*F(1167)+p1mup6*F(1173))+p1mup3*p2mup6*F(1174)-p2mu
     &   p6*p5mup3*F(1175)-p1mup3*p4mup6*F(1176)+p4mup6*p5mup3*F(1177
     &   )-p1mup3*p1mup6*F(1178)+p1mup6*p5mup3*F(1179)-p1mup3*p5mup6*
     &   F(1180)+p5mup3*p5mup6*F(1181)+mup3mup6*F(1182)
       Fa16(289) = p2mup3*p2mup6*F(1053)-p2mup3*p5mup6*F(1058)+p1mup3*
     &   p6mup6*F(1183)-p2mup3*p6mup6*F(1185)+p4mup3*p6mup6*F(1187)-p
     &   1mup3*(p2mup6*F(1184)+p4mup6*F(1186)+p1mup6*F(1188))-p1mup3*
     &   p5mup6*F(1190)-mup3mup6*F(1192)+p1mup6*p2mup3*F(1193)+p1mup6
     &   *p4mup3*F(1196)+p2mup6*p4mup3*F(1199)-p2mup3*p4mup6*F(1200)+
     &   p5mup3*(-(p5mup6*F(984))-p6mup6*F(1189)+p1mup6*F(1194)+p2mup
     &   6*F(1197)+p4mup6*F(1201))+p6mup3*(p6mup6*F(980)-p5mup6*F(119
     &   1)+p1mup6*F(1195)+p2mup6*F(1198)+p4mup6*F(1202))+p4mup3*p4mu
     &   p6*F(1203)+p4mup3*p5mup6*F(1204)
       Fa16(290) = 4*(p1mup3*F(1207)+p4mup3*F(1208)+p5mup3*F(1209)+p6m
     &   up3*F(1210))+p2mup3*F(1211)
       Fa16(291) = p2mup6*F(1214)+4*(p1mup6*F(1212)-p6mup6*F(1213)-p4m
     &   up6*F(1215))+p5mup6*F(1216)
       Fa16(292) = p6mup6*F(1218)+p4mup6*F(1220)+4*(p1mup6*F(1217)-p2m
     &   up6*F(1219)+p5mup6*F(1221))
       Fa16(293) = p1mup6*F(1222)+p6mup6*F(1223)+p2mup6*F(1224)+p4mup6
     &   *F(1225)+p5mup6*F(1226)
       Fa16(294) = p2mup3*F(1234)-p1mup3*F(1235)+p4mup3*F(1236)+p5mup3
     &   *F(1237)-p6mup3*F(1238)
       Fa16(295) = p1mup3*F(1239)+4*(p2mup3*F(1240)-p4mup3*F(1241)-p5m
     &   up3*F(1242)+p6mup3*F(1243))
       Fa16(296) = p5mup6*F(1257)+p6mup6*F(1258)-p2mup6*F(1259)-p1mup6
     &   *F(1260)+p4mup6*F(1261)
       Fa16(297) = p2mup3*F(1262)-p4mup3*F(1263)-p6mup3*F(1264)+p5mup3
     &   *F(1265)+p1mup3*F(1266)
       Fa16(298) = -(p6mup3*p6mup6*F(1268))-p4mup3*p6mup6*F(1269)-p5mu
     &   p6*p6mup3*F(1270)+p2mup6*p6mup3*F(1271)+p2mup6*p4mup3*F(1272
     &   )+p2mup3*(p6mup6*F(1267)+p5mup6*F(1273)+p4mup6*F(1274)+p1mup
     &   6*F(1275)+p2mup6*F(1276))+p4mup6*p6mup3*F(1277)+p1mup6*p4mup
     &   3*F(1278)+p4mup3*p4mup6*F(1279)+p4mup3*p5mup6*F(1280)+mup3mu
     &   p6*F(1282)-p1mup6*p6mup3*F(1283)+p1mup3*p6mup6*F(1286)+p5mup
     &   3*(p5mup6*F(1194)+p6mup6*F(1281)-p2mup6*F(1284)+p4mup6*F(128
     &   5)+p1mup6*F(1287))+p1mup3*p2mup6*F(1288)-p1mup3*p4mup6*F(128
     &   9)-p1mup3*p1mup6*F(1290)-p1mup3*p5mup6*F(1291)
       Fa16(299) = p6mup6*F(1298)-p1mup6*F(1299)-p2mup6*F(1300)+p4mup6
     &   *F(1301)+p5mup6*F(1302)
       Fa16(300) = p5mup6*F(1296)+p4mup6*F(1303)+p6mup6*F(1304)
       Fa16(301) = p1mup6*F(1296)+p4mup6*F(1305)-p6mup6*F(1306)
       Fa16(302) = p1mup6*F(1304)+p5mup6*F(1306)+p4mup6*F(1307)
       Fa16(303) = -(p2mup3*F(1308))+p4mup3*F(1309)+p6mup3*F(1310)+p1m
     &   up3*F(1311)+p5mup3*F(1312)
       Fa16(304) = p1mup3*F(1313)+p2mup3*F(1314)-p4mup3*F(1315)-p6mup3
     &   *F(1316)+p5mup3*F(1317)
       Return
       End
