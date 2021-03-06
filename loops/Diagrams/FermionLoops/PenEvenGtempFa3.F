c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenEvenG3(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup
     -   1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3
     -   mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3
     -   ,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,
     -   mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa3)
       IMPLICIT NONE
      Complex*16   p1mup1, p1mup2, p1mup3, p1mup4, p1mup5, p2mup1, 
     -          p2mup2, p2mup3, p2mup4, p2mup5, p3mup1, p3mup2, 
     -          p3mup3, p3mup4, p3mup5, p4mup1, p4mup2, p4mup3, 
     -          p4mup4, p4mup5, p5mup1, p5mup2, p5mup3, p5mup4, p5mup5
       Complex*16   mup1mup2, mup1mup3, mup1mup4, mup1mup5, mup2mup3, 
     -          mup2mup4, mup2mup5, mup3mup4, mup3mup5, mup4mup5
        common/FPenEvenG/F
      Complex*16 Fa(150),F(4300)
      Real*8 P(2563) 
        COMMON/PPenEvenG/P
       Complex*16 Fa3(21:30)
       COMMON/FaPenEvenG/Fa
       Fa3(21) = -(mup4mup5*p1mup1*p2mup2*p3mup3*F(1167))+mup3mup5*p1m
     -   up1*p2mup2*p3mup4*F(1167)+mup3mup4*p1mup1*p2mup2*p3mup5*F(11
     -   67)-mup4mup5*p1mup1*p2mup3*p3mup2*F(1168)+mup3mup5*p1mup1*p2
     -   mup4*p3mup2*F(1168)-mup2mup5*p1mup1*p2mup4*p3mup3*F(1169)+mu
     -   p2mup4*p1mup1*p2mup5*p3mup3*F(1169)+mup2mup5*p1mup1*p2mup3*p
     -   3mup4*F(1169)-mup2mup3*p1mup1*p2mup5*p3mup4*F(1169)+mup2mup4
     -   *p1mup1*p2mup3*p3mup5*F(1169)-mup2mup3*p1mup1*p2mup4*p3mup5*
     -   F(1169)+mup3mup4*p1mup1*p2mup5*p3mup2*F(1170)-mup4mup5*p1mup
     -   2*p2mup1*p3mup3*F(1172)+mup3mup5*p1mup2*p2mup1*p3mup4*F(1172
     -   )+mup3mup4*p1mup2*p2mup1*p3mup5*F(1172)-mup4mup5*p1mup2*p2mu
     -   p3*p3mup1*F(1173)+mup3mup5*p1mup2*p2mup4*p3mup1*F(1173)+mup2
     -   mup5*p1mup4*p2mup1*p3mup3*F(1174)-mup2mup4*p1mup5*p2mup1*p3m
     -   up3*F(1174)-mup1mup5*p1mup2*p2mup4*p3mup3*F(1174)-mup1mup2*p
     -   1mup4*p2mup5*p3mup3*F(1174)-mup2mup5*p1mup3*p2mup1*p3mup4*F(
     -   1174)+mup2mup3*p1mup5*p2mup1*p3mup4*F(1174)+mup1mup5*p1mup2*
     -   p2mup3*p3mup4*F(1174)+mup1mup2*p1mup3*p2mup5*p3mup4*F(1174)-
     -   mup2mup4*p1mup3*p2mup1*p3mup5*F(1174)+mup2mup3*p1mup4*p2mup1
     -   *p3mup5*F(1174)-mup1mup2*p1mup4*p2mup3*p3mup5*F(1174)+mup1mu
     -   p2*p1mup3*p2mup4*p3mup5*F(1174)+mup3mup4*p1mup2*p2mup5*p3mup
     -   1*F(1175)+mup3mup5*p1mup4*p2mup1*p3mup2*F(1177)-mup4mup5*p1m
     -   up3*p2mup2*p3mup1*F(1178)+mup3mup5*p1mup4*p2mup2*p3mup1*F(11
     -   78)-mup1mup5*p1mup4*p2mup2*p3mup3*F(1179)+mup1mup5*p1mup3*p2
     -   mup2*p3mup4*F(1179)+mup2mup5*p1mup4*p2mup3*p3mup1*F(1180)-mu
     -   p2mup5*p1mup3*p2mup4*p3mup1*F(1180)+mup2mup4*p1mup3*p2mup5*p
     -   3mup1*F(1180)-mup2mup3*p1mup4*p2mup5*p3mup1*F(1180)+mup3mup4
     -   *p1mup5*p2mup1*p3mup2*F(1183)-mup1mup5*p1mup4*p2mup3*p3mup2*
     -   F(1186)+mup1mup5*p1mup3*p2mup4*p3mup2*F(1186)-mup1mup3*(2*(p
     -   1mup2*p2mup4*p2mup5*F(1163)+p1mup5*p2mup2*p2mup4*F(1165))+p1
     -   mup2*(p2mup5*p3mup4+p2mup4*p3mup5)*F(1174)+p1mup4*(2*p2mup2*
     -   p2mup5*F(932)+p2mup2*p3mup5*F(1179)-p2mup5*p3mup2*F(1182)+4*
     -   p1mup5*p4mup2*F(1192)))+mup1mup4*(2*p1mup5*p2mup2*p2mup3*F(1
     -   165)+p1mup2*((p2mup5*p3mup3+p2mup3*p3mup5)*F(1174)+2*(p2mup3
     -   *p2mup5*F(1163)+p1mup5*p4mup3*F(1201)))+p1mup3*(p2mup2*p3mup
     -   5*F(1179)+p2mup5*p3mup2*F(1181)+4*p1mup5*p4mup2*F(1192)+2*(p
     -   2mup2*p2mup5*F(932)+p1mup2*p4mup5*F(1201))))+2*(mup2mup4*p1m
     -   up3*p2mup1*p2mup5*F(1049)-mup2mup3*p1mup4*p2mup1*p2mup5*F(10
     -   49)+p1mup1*p2mup2*(-(mup4mup5*p2mup3)+mup3mup5*p2mup4)*F(116
     -   0)+mup2mup4*p1mup1*p2mup3*p2mup5*F(1161)-mup2mup3*p1mup1*p2m
     -   up4*p2mup5*F(1161)-mup4mup5*p1mup2*p2mup1*p2mup3*F(1162)+mup
     -   3mup5*p1mup2*p2mup1*p2mup4*F(1162)-mup1mup2*p1mup4*p2mup3*p2
     -   mup5*F(1163)+mup1mup2*p1mup3*p2mup4*p2mup5*F(1163)+mup2mup4*
     -   p1mup5*p2mup1*p2mup3*F(1164)-mup1mup5*p1mup4*p2mup2*p2mup3*F
     -   (1165)+mup1mup5*p1mup3*p2mup2*p2mup4*F(1165)+mup2mup3*p1mup5
     -   *p2mup1*p2mup4*F(1166)+mup3mup4*p1mup1*p2mup2*p2mup5*F(1188)
     -   +mup3mup4*p1mup2*p2mup1*p2mup5*F(1189)-mup4mup5*p1mup3*p2mup
     -   1*p2mup2*F(1190)+mup3mup5*p1mup4*p2mup1*p2mup2*F(1190)+mup3m
     -   up4*p1mup5*p2mup1*p2mup2*F(1193)-mup4mup5*p1mup2*p1mup3*p4mu
     -   p1*F(1200)+mup3mup5*p1mup2*p1mup4*p4mup1*F(1200)-mup1mup5*p1
     -   mup2*p1mup4*p4mup3*F(1201)+mup1mup5*p1mup2*p1mup3*p4mup4*F(1
     -   201)-mup1mup3*p1mup2*(p1mup5*p4mup4+p1mup4*p4mup5)*F(1201)-m
     -   up3mup4*p1mup2*p1mup5*p4mup1*F(1202)+mup2mup5*p1mup4*p2mup1*
     -   p2mup3*F(1203)-mup2mup5*p1mup3*p2mup1*p2mup4*F(1203))+mup4mu
     -   p5*p1mup3*p2mup1*p3mup2*F(1204)
       Fa3(22) = mup1mup3*p1mup5*p2mup2*p3mup4*F(1207)+mup2mup4*p1mup5
     -   *p2mup3*p3mup1*F(1208)+mup1mup2*p1mup5*p2mup3*p3mup4*F(1210)
     -   +mup2mup3*p1mup5*p2mup4*p3mup1*F(1211)+mup1mup3*p1mup5*p2mup
     -   4*p3mup2*F(1212)+mup1mup2*p1mup5*p2mup4*p3mup3*F(1213)+mup4m
     -   up5*p1mup1*p2mup2*p4mup3*F(1214)+mup3mup5*p1mup1*p2mup2*p4mu
     -   p4*F(1215)+mup4mup5*p1mup1*p2mup3*p4mup2*F(1217)+mup2mup5*p1
     -   mup1*p2mup3*p4mup4*F(1218)+mup2mup4*p1mup1*p2mup3*p4mup5*F(1
     -   219)+mup3mup5*p1mup1*p2mup4*p4mup2*F(1220)+mup2mup5*p1mup1*p
     -   2mup4*p4mup3*F(1221)+mup2mup3*p1mup1*p2mup4*p4mup5*F(1222)+m
     -   up2mup4*p1mup1*p2mup5*p4mup3*F(1224)+mup2mup3*p1mup1*p2mup5*
     -   p4mup4*F(1225)+mup4mup5*p1mup2*p2mup1*p4mup3*F(1226)+mup3mup
     -   5*p1mup2*p2mup1*p4mup4*F(1227)+mup4mup5*p1mup2*p2mup3*p4mup1
     -   *F(1229)+mup1mup5*p1mup2*p2mup3*p4mup4*F(1230)+mup3mup5*p1mu
     -   p2*p2mup4*p4mup1*F(1232)+mup1mup5*p1mup2*p2mup4*p4mup3*F(123
     -   3)+mup1mup3*p1mup2*p2mup4*p4mup5*F(1234)+mup1mup3*p1mup2*p2m
     -   up5*p4mup4*F(1237)+mup4mup5*p1mup3*p2mup1*p4mup2*F(1238)+mup
     -   2mup5*p1mup3*p2mup1*p4mup4*F(1239)+mup2mup4*p1mup3*p2mup1*p4
     -   mup5*F(1240)+mup4mup5*p1mup3*p2mup2*p4mup1*F(1241)+mup1mup5*
     -   p1mup3*p2mup2*p4mup4*F(1242)+mup2mup5*p1mup3*p2mup4*p4mup1*F
     -   (1244)+mup1mup5*p1mup3*p2mup4*p4mup2*F(1245)+mup1mup2*p1mup3
     -   *p2mup4*p4mup5*F(1246)+mup2mup4*p1mup3*p2mup5*p4mup1*F(1247)
     -   +mup1mup2*p1mup3*p2mup5*p4mup4*F(1249)+mup3mup5*p1mup4*p2mup
     -   1*p4mup2*F(1250)+mup2mup5*p1mup4*p2mup1*p4mup3*F(1251)+mup2m
     -   up3*p1mup4*p2mup1*p4mup5*F(1252)+mup3mup5*p1mup4*p2mup2*p4mu
     -   p1*F(1253)+mup1mup5*p1mup4*p2mup2*p4mup3*F(1254)+mup1mup3*p1
     -   mup4*p2mup2*p4mup5*F(1255)+mup2mup5*p1mup4*p2mup3*p4mup1*F(1
     -   256)+mup1mup5*p1mup4*p2mup3*p4mup2*F(1257)+mup1mup2*p1mup4*p
     -   2mup3*p4mup5*F(1258)+mup2mup3*p1mup4*p2mup5*p4mup1*F(1259)+m
     -   up1mup3*p1mup4*p2mup5*p4mup2*F(1260)+mup1mup2*p1mup4*p2mup5*
     -   p4mup3*F(1261)+mup2mup4*p1mup5*p2mup1*p4mup3*F(1263)+mup2mup
     -   3*p1mup5*p2mup1*p4mup4*F(1264)+mup1mup3*p1mup5*p2mup2*p4mup4
     -   *F(1267)+mup2mup4*p1mup5*p2mup3*p4mup1*F(1268)+mup1mup2*p1mu
     -   p5*p2mup3*p4mup4*F(1270)+mup2mup3*p1mup5*p2mup4*p4mup1*F(127
     -   1)+mup1mup3*p1mup5*p2mup4*p4mup2*F(1272)+mup1mup2*p1mup5*p2m
     -   up4*p4mup3*F(1273)+mup2mup5*p1mup1*p3mup3*p3mup4*F(1275)+mup
     -   2mup3*p1mup1*p3mup4*p3mup5*F(1276)+mup4mup5*p1mup1*p3mup2*p4
     -   mup3*F(1286)+mup1mup2*p1mup4*p3mup3*p3mup5*F(1292)+mup3mup4*
     -   (p1mup1*p2mup2*p4mup5*F(1216)+p1mup1*p2mup5*p4mup2*F(1223)+p
     -   1mup2*p2mup1*p4mup5*F(1228)+p1mup2*p2mup5*p4mup1*F(1235)+2*(
     -   p1mup1*p3mup2*p3mup5*F(1294)+p1mup2*p3mup1*p3mup5*F(1296))+p
     -   1mup5*(p2mup2*p3mup1*F(1205)+p2mup1*p4mup2*F(1262)+p2mup2*p4
     -   mup1*F(1265)+2*p3mup1*p3mup2*F(1299)))+mup1mup4*(p1mup2*p2mu
     -   p3*p4mup5*F(1231)+p1mup2*p2mup5*p4mup3*F(1236)+p1mup3*p2mup2
     -   *p4mup5*F(1243)+p1mup3*p2mup5*p4mup2*F(1248)+2*p1mup3*p3mup2
     -   *p3mup5*F(1298)+p1mup5*(p2mup2*p3mup3*F(1206)+p2mup3*p3mup2*
     -   F(1209)+p2mup2*p4mup3*F(1266)+p2mup3*p4mup2*F(1269)+2*p3mup2
     -   *p3mup3*F(1300))+p1mup2*p3mup3*p3mup5*F(1302))+2*(mup3mup5*p
     -   1mup1*p3mup2*p3mup4*F(1274)+mup3mup5*p1mup2*p3mup1*p3mup4*F(
     -   1277)+mup1mup5*p1mup2*p3mup3*p3mup4*F(1278)+mup1mup3*p1mup2*
     -   p3mup4*p3mup5*F(1279)+mup2mup5*p1mup3*p3mup1*p3mup4*F(1280)+
     -   mup1mup5*p1mup3*p3mup2*p3mup4*F(1281)+mup1mup2*p1mup3*p3mup4
     -   *p3mup5*F(1282)+mup2mup3*p1mup5*p3mup1*p3mup4*F(1283)+mup1mu
     -   p3*p1mup5*p3mup2*p3mup4*F(1284)+mup1mup2*p1mup5*p3mup3*p3mup
     -   4*F(1285)+mup3mup5*p1mup4*p3mup1*p3mup2*F(1287)-mup2mup5*p1m
     -   up4*p3mup1*p3mup3*F(1288)+mup2mup3*p1mup4*p3mup1*p3mup5*F(12
     -   89)+mup1mup5*p1mup4*p3mup2*p3mup3*F(1290)+mup1mup3*p1mup4*p3
     -   mup2*p3mup5*F(1291)+mup4mup5*p1mup1*p3mup2*p3mup3*F(1293)+mu
     -   p4mup5*p1mup2*p3mup1*p3mup3*F(1295)+mup4mup5*p1mup3*p3mup1*p
     -   3mup2*F(1297)-mup2mup4*p1mup1*p3mup3*p3mup5*F(1301)-mup2mup4
     -   *p1mup3*p3mup1*p3mup5*F(1303)-mup2mup4*p1mup5*p3mup1*p3mup3*
     -   F(1304))
       Fa3(23) = -(mup1mup2*p1mup5*(p2mup4*p3mup3+p2mup3*p3mup4)*F(117
     -   4))-mup1mup3*p1mup5*p2mup2*p3mup4*F(1179)+mup1mup3*p1mup5*p2
     -   mup4*p3mup2*F(1182)+mup2mup3*p1mup5*p2mup4*p3mup1*F(1185)+mu
     -   p2mup4*p1mup5*p2mup3*p3mup1*F(1187)+mup3mup4*p1mup5*p2mup2*p
     -   3mup1*F(1305)+4*(mup4mup5*p1mup2*p3mup1*p3mup3*F(1307)+mup3m
     -   up5*p1mup2*p3mup1*p3mup4*F(1307)-mup3mup4*p1mup2*p3mup1*p3mu
     -   p5*F(1307))+mup3mup5*p1mup4*p2mup2*p4mup1*F(1312)-mup3mup4*p
     -   1mup5*p2mup2*p4mup1*F(1312)+mup3mup5*p1mup1*p2mup4*p4mup2*F(
     -   1312)-mup3mup4*p1mup1*p2mup5*p4mup2*F(1312)+mup4mup5*p1mup1*
     -   p2mup2*p4mup3*F(1312)+mup3mup5*p1mup1*p2mup2*p4mup4*F(1312)-
     -   mup3mup4*p1mup1*p2mup2*p4mup5*F(1312)+mup4mup5*p1mup1*p2mup3
     -   *p4mup2*F(1313)+mup2mup5*p1mup4*p2mup3*p4mup1*F(1314)-mup2mu
     -   p4*p1mup5*p2mup3*p4mup1*F(1314)+mup2mup3*p1mup5*p2mup4*p4mup
     -   1*F(1314)-mup2mup3*p1mup4*p2mup5*p4mup1*F(1314)+mup2mup5*p1m
     -   up1*p2mup4*p4mup3*F(1314)-mup2mup4*p1mup1*p2mup5*p4mup3*F(13
     -   14)+mup2mup5*p1mup1*p2mup3*p4mup4*F(1314)-mup2mup3*p1mup1*p2
     -   mup5*p4mup4*F(1314)-mup2mup4*p1mup1*p2mup3*p4mup5*F(1314)+mu
     -   p2mup3*p1mup1*p2mup4*p4mup5*F(1314)+mup4mup5*p1mup2*p2mup1*p
     -   4mup3*F(1315)+mup3mup5*p1mup2*p2mup1*p4mup4*F(1315)-mup3mup4
     -   *p1mup2*p2mup1*p4mup5*F(1315)+mup4mup5*p1mup2*p2mup3*p4mup1*
     -   F(1316)+mup3mup5*p1mup2*p2mup4*p4mup1*F(1317)-mup3mup4*p1mup
     -   2*p2mup5*p4mup1*F(1317)+mup4mup5*p1mup3*p2mup1*p4mup2*F(1318
     -   )+mup4mup5*p1mup3*p2mup2*p4mup1*F(1319)+mup1mup5*p1mup3*p2mu
     -   p4*p4mup2*F(1320)-mup1mup3*p1mup5*p2mup4*p4mup2*F(1320)+mup1
     -   mup3*p1mup4*p2mup5*p4mup2*F(1320)+mup1mup5*p1mup4*p2mup2*p4m
     -   up3*F(1320)+mup1mup5*p1mup3*p2mup2*p4mup4*F(1320)-mup1mup3*p
     -   1mup5*p2mup2*p4mup4*F(1320)+mup1mup3*p1mup4*p2mup2*p4mup5*F(
     -   1320)+mup2mup5*p1mup3*p2mup4*p4mup1*F(1321)-mup2mup4*p1mup3*
     -   p2mup5*p4mup1*F(1321)-mup3mup5*p1mup4*p2mup1*p4mup2*F(1322)+
     -   mup3mup4*p1mup5*p2mup1*p4mup2*F(1322)+mup1mup5*p1mup4*p2mup3
     -   *p4mup2*F(1323)+2*(mup4mup5*p1mup1*p3mup2*p3mup3*F(1306)+mup
     -   3mup5*p1mup1*p3mup2*p3mup4*F(1306)-mup3mup4*p1mup1*p3mup2*p3
     -   mup5*F(1306)+mup4mup5*p1mup3*p3mup1*p3mup2*F(1308)+mup1mup5*
     -   p1mup4*p3mup2*p3mup3*F(1309)+mup1mup5*p1mup3*p3mup2*p3mup4*F
     -   (1309)-mup1mup3*p1mup5*p3mup2*p3mup4*F(1309)+mup1mup3*p1mup4
     -   *p3mup2*p3mup5*F(1309)+mup3mup5*p1mup4*p3mup1*p3mup2*F(1310)
     -   +mup3mup4*p1mup5*p3mup1*p3mup2*F(1311)+mup4mup5*p1mup1*p3mup
     -   2*p4mup3*F(1324))-mup1mup2*(p1mup5*p2mup4*p4mup3-p1mup4*p2mu
     -   p5*p4mup3+p1mup5*p2mup3*p4mup4-p1mup3*p2mup5*p4mup4-p1mup4*p
     -   2mup3*p4mup5+p1mup3*p2mup4*p4mup5)*F(4285)+(-(mup2mup5*p2mup
     -   1*(p1mup4*p4mup3+p1mup3*p4mup4))+mup2mup4*p2mup1*(p1mup5*p4m
     -   up3+p1mup3*p4mup5)+mup2mup3*p2mup1*(p1mup5*p4mup4-p1mup4*p4m
     -   up5)+p1mup2*(mup1mup5*p2mup4*p4mup3+mup1mup5*p2mup3*p4mup4-m
     -   up1mup3*p2mup5*p4mup4+mup1mup3*p2mup4*p4mup5))*F(4285)-mup1m
     -   up4*(2*p1mup3*p3mup2*p3mup5*F(1309)+p1mup3*(p2mup5*p4mup2+p2
     -   mup2*p4mup5)*F(1320)+p1mup5*(p2mup2*p3mup3*F(1179)+p2mup3*p3
     -   mup2*F(1186)+2*p3mup2*p3mup3*F(1309)+p2mup2*p4mup3*F(1320)+p
     -   2mup3*p4mup2*F(1323))+p1mup2*(p2mup5*p4mup3+p2mup3*p4mup5)*F
     -   (4285))
       Fa3(24) = -(mup2mup3*p1mup5*p2mup4*p3mup1*F(1195))-mup1mup4*p1m
     -   up5*p2mup3*p3mup2*F(1199)+mup1mup3*p1mup5*p2mup4*p3mup2*F(11
     -   99)-mup1mup4*p1mup5*p2mup2*p3mup3*F(1199)+mup1mup3*p1mup5*p2
     -   mup2*p3mup4*F(1199)+mup4mup5*p1mup3*p2mup2*p4mup1*F(1312)-mu
     -   p3mup5*p1mup4*p2mup2*p4mup1*F(1312)+mup4mup5*p1mup1*p2mup3*p
     -   4mup2*F(1312)-mup3mup5*p1mup1*p2mup4*p4mup2*F(1312)+mup4mup5
     -   *p1mup1*p2mup2*p4mup3*F(1312)-mup3mup5*p1mup1*p2mup2*p4mup4*
     -   F(1312)-mup2mup5*p1mup4*p2mup3*p4mup1*F(1314)+mup2mup5*p1mup
     -   3*p2mup4*p4mup1*F(1314)-mup2mup3*p1mup5*p2mup4*p4mup1*F(1314
     -   )+mup2mup3*p1mup4*p2mup5*p4mup1*F(1314)+mup2mup5*p1mup1*p2mu
     -   p4*p4mup3*F(1314)-mup2mup5*p1mup1*p2mup3*p4mup4*F(1314)+mup2
     -   mup3*p1mup1*p2mup5*p4mup4*F(1314)-mup2mup3*p1mup1*p2mup4*p4m
     -   up5*F(1314)+mup4mup5*p1mup2*p2mup1*p4mup3*F(1315)-mup3mup5*p
     -   1mup2*p2mup1*p4mup4*F(1315)+mup4mup5*p1mup2*p2mup3*p4mup1*F(
     -   1317)-mup3mup5*p1mup2*p2mup4*p4mup1*F(1317)+mup1mup5*p1mup4*
     -   p2mup3*p4mup2*F(1320)-mup1mup4*p1mup5*p2mup3*p4mup2*F(1320)-
     -   mup1mup5*p1mup3*p2mup4*p4mup2*F(1320)+mup1mup3*p1mup5*p2mup4
     -   *p4mup2*F(1320)+mup1mup4*p1mup3*p2mup5*p4mup2*F(1320)-mup1mu
     -   p3*p1mup4*p2mup5*p4mup2*F(1320)+mup1mup5*p1mup4*p2mup2*p4mup
     -   3*F(1320)-mup1mup4*p1mup5*p2mup2*p4mup3*F(1320)-mup1mup5*p1m
     -   up3*p2mup2*p4mup4*F(1320)+mup1mup3*p1mup5*p2mup2*p4mup4*F(13
     -   20)+mup1mup4*p1mup3*p2mup2*p4mup5*F(1320)-mup1mup3*p1mup4*p2
     -   mup2*p4mup5*F(1320)-mup4mup5*p1mup3*p2mup1*p4mup2*F(1322)+mu
     -   p3mup5*p1mup4*p2mup1*p4mup2*F(1322)+4*(mup4mup5*p1mup2*p3mup
     -   1*p3mup3*F(1326)-mup3mup5*p1mup2*p3mup1*p3mup4*F(1326))+mup3
     -   mup4*(p1mup1*(p2mup5*p4mup2+p2mup2*p4mup5)*F(1312)+2*p1mup1*
     -   p3mup2*p3mup5*F(1325)+p1mup5*(p2mup2*p3mup1*F(1194)+p2mup2*p
     -   4mup1*F(1312)-p2mup1*p4mup2*F(1322)+2*p3mup1*p3mup2*F(1325))
     -   +p1mup2*(p2mup1*p4mup5*F(1315)+p2mup5*p4mup1*F(1317)+4*p3mup
     -   1*p3mup5*F(1326)))-mup1mup2*p1mup5*p2mup4*p3mup3*F(4284)+mup
     -   1mup2*p1mup5*p2mup3*p3mup4*F(4284)-mup2mup5*p1mup4*p2mup1*p4
     -   mup3*F(4285)+mup1mup5*p1mup2*p2mup4*p4mup3*F(4285)-mup1mup2*
     -   p1mup5*p2mup4*p4mup3*F(4285)-mup1mup4*p1mup2*p2mup5*p4mup3*F
     -   (4285)+mup1mup2*p1mup4*p2mup5*p4mup3*F(4285)+mup2mup5*p1mup3
     -   *p2mup1*p4mup4*F(4285)-mup2mup3*p1mup5*p2mup1*p4mup4*F(4285)
     -   -mup1mup5*p1mup2*p2mup3*p4mup4*F(4285)+mup1mup2*p1mup5*p2mup
     -   3*p4mup4*F(4285)+mup1mup3*p1mup2*p2mup5*p4mup4*F(4285)-mup1m
     -   up2*p1mup3*p2mup5*p4mup4*F(4285)+mup2mup3*p1mup4*p2mup1*p4mu
     -   p5*F(4285)+mup1mup4*p1mup2*p2mup3*p4mup5*F(4285)-mup1mup2*p1
     -   mup4*p2mup3*p4mup5*F(4285)-mup1mup3*p1mup2*p2mup4*p4mup5*F(4
     -   285)+mup1mup2*p1mup3*p2mup4*p4mup5*F(4285)+mup2mup4*(-((p1mu
     -   p3*p2mup5*p4mup1+p1mup1*p2mup5*p4mup3-p1mup1*p2mup3*p4mup5)*
     -   F(1314))-p1mup3*p2mup1*p4mup5*F(4285)+p1mup5*(p2mup3*p3mup1*
     -   F(1195)+p2mup3*p4mup1*F(1314)+p2mup1*p4mup3*F(4285)))+2*(mup
     -   4mup5*p1mup1*p3mup2*p4mup3*F(1324)+mup4mup5*p1mup3*p3mup1*p3
     -   mup2*F(1325)-mup3mup5*p1mup4*p3mup1*p3mup2*F(1325)+mup4mup5*
     -   p1mup1*p3mup2*p3mup3*F(1325)-mup3mup5*p1mup1*p3mup2*p3mup4*F
     -   (1325)+p3mup2*(mup1mup5*p1mup4*p3mup3-mup1mup4*p1mup5*p3mup3
     -   -mup1mup5*p1mup3*p3mup4+mup1mup3*p1mup5*p3mup4+mup1mup4*p1mu
     -   p3*p3mup5-mup1mup3*p1mup4*p3mup5)*F(4287))
       Fa3(25) = mup1mup2*p1mup5*(p2mup4*p3mup3-p2mup3*p3mup4)*F(1174)
     -   -mup1mup3*p1mup5*p2mup2*p3mup4*F(1179)+mup3mup4*p1mup5*p2mup
     -   2*p3mup1*F(1184)+mup2mup4*p1mup5*p2mup3*p3mup1*F(1185)-mup1m
     -   up3*p1mup5*p2mup4*p3mup2*F(1186)+mup2mup3*p1mup5*p2mup4*p3mu
     -   p1*F(1187)+4*(-(mup4mup5*p1mup2*p3mup1*p3mup3*F(1307))+mup3m
     -   up5*p1mup2*p3mup1*p3mup4*F(1307)+mup3mup4*p1mup2*p3mup1*p3mu
     -   p5*F(1307))-mup4mup5*p1mup1*p2mup2*p4mup3*F(1327)+mup3mup5*p
     -   1mup1*p2mup2*p4mup4*F(1327)+mup3mup4*p1mup1*p2mup2*p4mup5*F(
     -   1327)-mup4mup5*p1mup1*p2mup3*p4mup2*F(1328)+mup3mup5*p1mup1*
     -   p2mup4*p4mup2*F(1328)-mup2mup5*p1mup1*p2mup4*p4mup3*F(1329)+
     -   mup2mup4*p1mup1*p2mup5*p4mup3*F(1329)+mup2mup5*p1mup1*p2mup3
     -   *p4mup4*F(1329)-mup2mup3*p1mup1*p2mup5*p4mup4*F(1329)+mup2mu
     -   p4*p1mup1*p2mup3*p4mup5*F(1329)-mup2mup3*p1mup1*p2mup4*p4mup
     -   5*F(1329)+mup3mup4*p1mup1*p2mup5*p4mup2*F(1330)-mup4mup5*p1m
     -   up2*p2mup1*p4mup3*F(1331)+mup3mup5*p1mup2*p2mup1*p4mup4*F(13
     -   31)+mup3mup4*p1mup2*p2mup1*p4mup5*F(1331)-mup4mup5*p1mup2*p2
     -   mup3*p4mup1*F(1332)+mup3mup5*p1mup2*p2mup4*p4mup1*F(1332)+mu
     -   p2mup5*p1mup4*p2mup1*p4mup3*F(1333)-mup2mup4*p1mup5*p2mup1*p
     -   4mup3*F(1333)-mup1mup5*p1mup2*p2mup4*p4mup3*F(1333)-mup2mup5
     -   *p1mup3*p2mup1*p4mup4*F(1333)+mup2mup3*p1mup5*p2mup1*p4mup4*
     -   F(1333)+mup1mup5*p1mup2*p2mup3*p4mup4*F(1333)-mup1mup3*p1mup
     -   2*p2mup5*p4mup4*F(1333)-mup2mup4*p1mup3*p2mup1*p4mup5*F(1333
     -   )+mup2mup3*p1mup4*p2mup1*p4mup5*F(1333)-mup1mup3*p1mup2*p2mu
     -   p4*p4mup5*F(1333)+mup1mup2*(p1mup5*p2mup4*p4mup3-p1mup4*p2mu
     -   p5*p4mup3-p1mup5*p2mup3*p4mup4+p1mup3*p2mup5*p4mup4-p1mup4*p
     -   2mup3*p4mup5+p1mup3*p2mup4*p4mup5)*F(1333)+mup3mup4*p1mup2*p
     -   2mup5*p4mup1*F(1334)+mup4mup5*p1mup3*p2mup1*p4mup2*F(1335)-m
     -   up4mup5*p1mup3*p2mup2*p4mup1*F(1336)+mup3mup5*p1mup4*p2mup2*
     -   p4mup1*F(1336)-mup1mup5*p1mup4*p2mup2*p4mup3*F(1337)+mup1mup
     -   5*p1mup3*p2mup2*p4mup4*F(1337)-mup1mup3*p1mup5*p2mup2*p4mup4
     -   *F(1337)-mup1mup3*p1mup4*p2mup2*p4mup5*F(1337)+mup2mup5*p1mu
     -   p4*p2mup3*p4mup1*F(1338)-mup2mup5*p1mup3*p2mup4*p4mup1*F(133
     -   8)+mup2mup4*p1mup3*p2mup5*p4mup1*F(1338)-mup2mup3*p1mup4*p2m
     -   up5*p4mup1*F(1338)-mup1mup5*p1mup4*p2mup3*p4mup2*F(1339)+mup
     -   1mup5*p1mup3*p2mup4*p4mup2*F(1339)-mup1mup3*p1mup5*p2mup4*p4
     -   mup2*F(1339)+mup1mup4*(p1mup2*(p2mup5*p4mup3+p2mup3*p4mup5)*
     -   F(1333)+p1mup5*(p2mup2*p3mup3*F(1179)+p2mup3*p3mup2*F(1186)+
     -   2*p3mup2*p3mup3*F(1309)+p2mup2*p4mup3*F(1337)+p2mup3*p4mup2*
     -   F(1339))+p1mup3*(2*p3mup2*p3mup5*F(1309)+p2mup2*p4mup5*F(133
     -   7)+p2mup5*p4mup2*F(1340)))+mup3mup5*p1mup4*p2mup1*p4mup2*F(1
     -   341)+mup1mup3*p1mup4*p2mup5*p4mup2*F(1342)+mup3mup4*p1mup5*p
     -   2mup1*p4mup2*F(1343)+mup3mup4*p1mup5*p2mup2*p4mup1*F(1344)+m
     -   up2mup4*p1mup5*p2mup3*p4mup1*F(1345)+mup2mup3*p1mup5*p2mup4*
     -   p4mup1*F(1346)+2*(-(mup4mup5*p1mup1*p3mup2*p3mup3*F(1306))+m
     -   up3mup5*p1mup1*p3mup2*p3mup4*F(1306)+mup3mup4*p1mup1*p3mup2*
     -   p3mup5*F(1306)-mup4mup5*p1mup3*p3mup1*p3mup2*F(1308)+mup3mup
     -   5*p1mup4*p3mup1*p3mup2*F(1308)-mup1mup5*p1mup4*p3mup2*p3mup3
     -   *F(1309)+mup1mup5*p1mup3*p3mup2*p3mup4*F(1309)-mup1mup3*p1mu
     -   p5*p3mup2*p3mup4*F(1309)-mup1mup3*p1mup4*p3mup2*p3mup5*F(130
     -   9)+mup3mup4*p1mup5*p3mup1*p3mup2*F(1310)-mup4mup5*p1mup1*p3m
     -   up2*p4mup3*F(1347))
       Fa3(26) = mup4mup5*p1mup1*p3mup3*p4mup2*F(1350)+mup2mup5*p1mup1
     -   *p3mup3*p4mup4*F(1351)+mup2mup4*p1mup1*p3mup3*p4mup5*F(1352)
     -   +mup2mup5*p1mup1*p3mup4*p4mup3*F(1354)+mup2mup3*p1mup1*p3mup
     -   4*p4mup5*F(1355)+mup2mup4*p1mup1*p3mup5*p4mup3*F(1357)+mup2m
     -   up3*p1mup1*p3mup5*p4mup4*F(1358)+mup4mup5*p1mup2*p3mup1*p4mu
     -   p3*F(1359)+mup4mup5*p1mup2*p3mup3*p4mup1*F(1362)+mup1mup5*p1
     -   mup2*p3mup3*p4mup4*F(1363)+mup1mup4*p1mup2*p3mup3*p4mup5*F(1
     -   364)+mup1mup5*p1mup2*p3mup4*p4mup3*F(1366)+mup1mup3*p1mup2*p
     -   3mup4*p4mup5*F(1367)+mup1mup4*p1mup2*p3mup5*p4mup3*F(1369)+m
     -   up1mup3*p1mup2*p3mup5*p4mup4*F(1370)+mup4mup5*p1mup3*p3mup1*
     -   p4mup2*F(1371)+mup2mup5*p1mup3*p3mup1*p4mup4*F(1372)+mup2mup
     -   4*p1mup3*p3mup1*p4mup5*F(1373)+mup4mup5*p1mup3*p3mup2*p4mup1
     -   *F(1374)+mup1mup5*p1mup3*p3mup2*p4mup4*F(1375)+mup1mup4*p1mu
     -   p3*p3mup2*p4mup5*F(1376)+mup2mup5*p1mup3*p3mup4*p4mup1*F(137
     -   7)+mup1mup5*p1mup3*p3mup4*p4mup2*F(1378)+mup1mup2*p1mup3*p3m
     -   up4*p4mup5*F(1379)+mup2mup4*p1mup3*p3mup5*p4mup1*F(1380)+mup
     -   1mup4*p1mup3*p3mup5*p4mup2*F(1381)+mup1mup2*p1mup3*p3mup5*p4
     -   mup4*F(1382)+mup2mup5*p1mup4*p3mup1*p4mup3*F(1384)+mup2mup3*
     -   p1mup4*p3mup1*p4mup5*F(1385)+mup1mup5*p1mup4*p3mup2*p4mup3*F
     -   (1387)+mup1mup3*p1mup4*p3mup2*p4mup5*F(1388)+mup2mup5*p1mup4
     -   *p3mup3*p4mup1*F(1389)+mup1mup5*p1mup4*p3mup3*p4mup2*F(1390)
     -   +mup1mup2*p1mup4*p3mup3*p4mup5*F(1391)+mup2mup3*p1mup4*p3mup
     -   5*p4mup1*F(1392)+mup1mup3*p1mup4*p3mup5*p4mup2*F(1393)+mup1m
     -   up2*p1mup4*p3mup5*p4mup3*F(1394)+mup2mup4*p1mup5*p3mup1*p4mu
     -   p3*F(1396)+mup2mup3*p1mup5*p3mup1*p4mup4*F(1397)+mup1mup4*p1
     -   mup5*p3mup2*p4mup3*F(1399)+mup1mup3*p1mup5*p3mup2*p4mup4*F(1
     -   400)+mup2mup4*p1mup5*p3mup3*p4mup1*F(1401)+mup1mup4*p1mup5*p
     -   3mup3*p4mup2*F(1402)+mup1mup2*p1mup5*p3mup3*p4mup4*F(1403)+m
     -   up2mup3*p1mup5*p3mup4*p4mup1*F(1404)+mup1mup3*p1mup5*p3mup4*
     -   p4mup2*F(1405)+mup1mup2*p1mup5*p3mup4*p4mup3*F(1406)+mup3mup
     -   5*(p1mup2*p3mup1*p4mup4*F(1360)+p1mup2*p3mup4*p4mup1*F(1365)
     -   +p1mup4*p3mup1*p4mup2*F(1383)+p1mup4*p3mup2*p4mup1*F(1386)+p
     -   1mup1*(p3mup2*p4mup4*F(1348)+p3mup4*p4mup2*F(1353)-2*p4mup2*
     -   p4mup4*F(1414))-2*(p2mup1*p2mup2*p3mup4*F(1409)+p1mup2*p4mup
     -   1*p4mup4*F(1417)+p1mup4*p4mup1*p4mup2*F(1422)))+mup3mup4*(p1
     -   mup1*p3mup2*p4mup5*F(1349)+p1mup1*p3mup5*p4mup2*F(1356)+p1mu
     -   p2*p3mup1*p4mup5*F(1361)+p1mup2*p3mup5*p4mup1*F(1368)+p1mup5
     -   *p3mup1*p4mup2*F(1395)+p1mup5*p3mup2*p4mup1*F(1398)+2*(p2mup
     -   1*p2mup2*p3mup5*F(1410)-p1mup1*p4mup2*p4mup5*F(1415)-p1mup2*
     -   p4mup1*p4mup5*F(1418)+p1mup5*p4mup1*p4mup2*F(1425)))+mup1mup
     -   3*p1mup2*p4mup4*p4mup5*F(1433)-2*(mup2mup5*p2mup1*p2mup3*p3m
     -   up4*F(1407)-mup4mup5*p2mup1*p2mup2*p3mup3*F(1408)-mup4mup5*p
     -   2mup1*p2mup3*p3mup2*F(1411)-mup2mup4*p2mup1*p2mup3*p3mup5*F(
     -   1412)-mup4mup5*p1mup1*p4mup2*p4mup3*F(1413)-mup4mup5*p1mup2*
     -   p4mup1*p4mup3*F(1416)-mup4mup5*p1mup3*p4mup1*p4mup2*F(1419)+
     -   mup1mup5*p1mup3*p4mup2*p4mup4*F(1420)+mup1mup4*p1mup3*p4mup2
     -   *p4mup5*F(1421)-mup1mup5*p1mup4*p4mup2*p4mup3*F(1423)-mup1mu
     -   p3*p1mup4*p4mup2*p4mup5*F(1424)+mup1mup4*p1mup5*p4mup2*p4mup
     -   3*F(1426)-mup1mup3*p1mup5*p4mup2*p4mup4*F(1427)+mup2mup5*p1m
     -   up1*p4mup3*p4mup4*F(1428)-mup2mup4*p1mup1*p4mup3*p4mup5*F(14
     -   29)-mup2mup3*p1mup1*p4mup4*p4mup5*F(1430)-mup1mup5*p1mup2*p4
     -   mup3*p4mup4*F(1431)-mup1mup4*p1mup2*p4mup3*p4mup5*F(1432)-mu
     -   p2mup5*p1mup3*p4mup1*p4mup4*F(1434)-mup1mup2*p1mup3*p4mup4*p
     -   4mup5*F(1435)+mup2mup5*p1mup4*p4mup1*p4mup3*F(1436)-mup2mup3
     -   *p1mup4*p4mup1*p4mup5*F(1437)-mup1mup2*p1mup4*p4mup3*p4mup5*
     -   F(1438)-mup2mup4*p1mup5*p4mup1*p4mup3*F(1439)-mup1mup2*p1mup
     -   5*p4mup3*p4mup4*F(1440)-mup2mup4*p1mup3*p4mup1*p4mup5*F(1441
     -   )-mup2mup3*p1mup5*p4mup1*p4mup4*F(1442))+16*(p1mup1*p1mup2*p
     -   1mup3*p1mup4*p4mup5*F(1443)+p1mup1*p1mup2*p1mup3*p1mup5*p4mu
     -   p4*F(1444)+p1mup2*p1mup3*p1mup4*p1mup5*p4mup1*F(1445))+32*p1
     -   mup1*p1mup4*p1mup5*(p1mup3*p4mup2+p1mup2*p4mup3)*F(1446)
       Fa3(27) = mup4mup5*p2mup1*p2mup2*p3mup3*F(1447)+mup4mup5*p2mup1
     -   *p2mup3*p3mup2*F(1448)+mup2mup5*p2mup1*p2mup3*p3mup4*F(1449)
     -   -mup2mup4*p2mup1*p2mup3*p3mup5*F(1449)+mup4mup5*p1mup1*p3mup
     -   3*p4mup2*F(1450)+mup4mup5*p1mup3*p3mup1*p4mup2*F(1452)+mup1m
     -   up5*p1mup4*p3mup3*p4mup2*F(1453)-mup1mup4*p1mup5*p3mup3*p4mu
     -   p2*F(1453)+mup1mup5*p1mup3*p3mup4*p4mup2*F(1453)-mup1mup3*p1
     -   mup5*p3mup4*p4mup2*F(1453)-mup1mup4*p1mup3*p3mup5*p4mup2*F(1
     -   453)+mup1mup3*p1mup4*p3mup5*p4mup2*F(1453)+mup4mup5*p1mup3*p
     -   3mup2*p4mup1*F(1457)+mup4mup5*p1mup1*p4mup2*p4mup3*F(1458)+2
     -   *(mup4mup5*p1mup2*p3mup3*p4mup1*F(1451)+mup4mup5*p1mup2*p3mu
     -   p1*p4mup3*F(1456)+mup4mup5*p1mup2*p4mup1*p4mup3*F(1459))+mup
     -   3mup5*(p2mup1*p2mup2*p3mup4*F(1447)+p1mup4*(p3mup2*p4mup1*F(
     -   1324)+p3mup1*p4mup2*F(1454)+p4mup1*p4mup2*F(1458))+p1mup1*(p
     -   3mup2*p4mup4*F(1324)+p3mup4*p4mup2*F(1450)+p4mup2*p4mup4*F(1
     -   458))+2*p1mup2*(p3mup4*p4mup1*F(1451)+p3mup1*p4mup4*F(1456)+
     -   p4mup1*p4mup4*F(1459)))-mup3mup4*(p2mup1*p2mup2*p3mup5*F(144
     -   7)+p1mup5*(p3mup2*p4mup1*F(1324)-p3mup1*p4mup2*F(1455)+p4mup
     -   1*p4mup2*F(1458))+p1mup1*(p3mup2*p4mup5*F(1324)+p3mup5*p4mup
     -   2*F(1450)+p4mup2*p4mup5*F(1458))+2*p1mup2*(p3mup5*p4mup1*F(1
     -   451)+p3mup1*p4mup5*F(1456)+p4mup1*p4mup5*F(1459)))+mup4mup5*
     -   p1mup3*p4mup1*p4mup2*F(1460)+mup1mup5*p1mup4*p3mup2*p4mup3*F
     -   (4286)-mup1mup4*p1mup5*p3mup2*p4mup3*F(4286)+mup1mup5*p1mup3
     -   *p3mup2*p4mup4*F(4286)-mup1mup3*p1mup5*p3mup2*p4mup4*F(4286)
     -   -mup1mup4*p1mup3*p3mup2*p4mup5*F(4286)+mup1mup3*p1mup4*p3mup
     -   2*p4mup5*F(4286)+p4mup2*(mup1mup5*p1mup4*p4mup3-mup1mup4*p1m
     -   up5*p4mup3+mup1mup5*p1mup3*p4mup4-mup1mup3*p1mup5*p4mup4-mup
     -   1mup4*p1mup3*p4mup5+mup1mup3*p1mup4*p4mup5)*F(4288)
       Fa3(28) = -(mup3mup4*p1mup5*p3mup2*p4mup1*F(1324))-mup3mup4*p1m
     -   up5*p3mup1*p4mup2*F(1324)-mup3mup4*p1mup1*p3mup5*p4mup2*F(13
     -   24)-mup3mup4*p1mup1*p3mup2*p4mup5*F(1324)-mup3mup4*p1mup5*p4
     -   mup1*p4mup2*F(1458)-mup3mup4*p1mup1*p4mup2*p4mup5*F(1458)-2*
     -   (mup3mup4*p1mup2*p3mup5*p4mup1*F(1456)+mup3mup4*p1mup2*p3mup
     -   1*p4mup5*F(1456)+mup3mup4*p1mup2*p4mup1*p4mup5*F(1459))-mup3
     -   mup4*p2mup1*p2mup2*p3mup5*F(1461)-mup4mup5*(p1mup3*(p3mup2*p
     -   4mup1+p3mup1*p4mup2)*F(1324)+p1mup3*p4mup1*p4mup2*F(1458)+p1
     -   mup1*p4mup2*(p3mup3*F(1324)+p4mup3*F(1458))+2*p1mup2*((p3mup
     -   3*p4mup1+p3mup1*p4mup3)*F(1456)+p4mup1*p4mup3*F(1459))+p2mup
     -   1*(p2mup3*p3mup2+p2mup2*p3mup3)*F(1461))+mup3mup5*(p1mup4*(p
     -   3mup2*p4mup1+p3mup1*p4mup2)*F(1324)+p1mup1*(p3mup4*p4mup2+p3
     -   mup2*p4mup4)*F(1324)+p1mup4*p4mup1*p4mup2*F(1458)+p1mup1*p4m
     -   up2*p4mup4*F(1458)+2*p1mup2*((p3mup4*p4mup1+p3mup1*p4mup4)*F
     -   (1456)+p4mup1*p4mup4*F(1459))+p2mup1*p2mup2*p3mup4*F(1461))-
     -   mup1mup5*p1mup4*p3mup3*p4mup2*F(4286)+mup1mup4*p1mup5*p3mup3
     -   *p4mup2*F(4286)+mup1mup5*p1mup3*p3mup4*p4mup2*F(4286)-mup1mu
     -   p3*p1mup5*p3mup4*p4mup2*F(4286)-mup1mup4*p1mup3*p3mup5*p4mup
     -   2*F(4286)+mup1mup3*p1mup4*p3mup5*p4mup2*F(4286)-mup1mup5*p1m
     -   up4*p3mup2*p4mup3*F(4286)+mup1mup4*p1mup5*p3mup2*p4mup3*F(42
     -   86)+mup1mup5*p1mup3*p3mup2*p4mup4*F(4286)-mup1mup3*p1mup5*p3
     -   mup2*p4mup4*F(4286)-mup1mup4*p1mup3*p3mup2*p4mup5*F(4286)+mu
     -   p1mup3*p1mup4*p3mup2*p4mup5*F(4286)-mup1mup5*p1mup4*p4mup2*p
     -   4mup3*F(4288)+mup1mup4*p1mup5*p4mup2*p4mup3*F(4288)+mup1mup5
     -   *p1mup3*p4mup2*p4mup4*F(4288)-mup1mup3*p1mup5*p4mup2*p4mup4*
     -   F(4288)-mup1mup4*p1mup3*p4mup2*p4mup5*F(4288)+mup1mup3*p1mup
     -   4*p4mup2*p4mup5*F(4288)+p2mup1*p2mup3*(mup2mup5*p3mup4-mup2m
     -   up4*p3mup5)*F(4289)
       Fa3(29) = -(mup4mup5*p2mup1*p2mup2*p3mup3*F(1447))-mup4mup5*p2m
     -   up1*p2mup3*p3mup2*F(1448)+mup2mup5*p2mup1*p2mup3*p3mup4*F(14
     -   49)+mup2mup4*p2mup1*p2mup3*p3mup5*F(1449)-mup4mup5*p1mup1*p3
     -   mup3*p4mup2*F(1450)-mup4mup5*p1mup3*p3mup1*p4mup2*F(1452)-mu
     -   p1mup5*p1mup4*p3mup3*p4mup2*F(1453)+mup1mup4*p1mup5*p3mup3*p
     -   4mup2*F(1453)+mup1mup5*p1mup3*p3mup4*p4mup2*F(1453)-mup1mup3
     -   *p1mup5*p3mup4*p4mup2*F(1453)+mup1mup4*p1mup3*p3mup5*p4mup2*
     -   F(1453)-mup1mup3*p1mup4*p3mup5*p4mup2*F(1453)-mup4mup5*p1mup
     -   3*p3mup2*p4mup1*F(1463)-mup1mup5*p1mup4*p3mup2*p4mup3*F(1464
     -   )+mup1mup4*p1mup5*p3mup2*p4mup3*F(1464)+mup1mup5*p1mup3*p3mu
     -   p2*p4mup4*F(1464)-mup1mup3*p1mup5*p3mup2*p4mup4*F(1464)+mup1
     -   mup4*p1mup3*p3mup2*p4mup5*F(1464)-mup1mup3*p1mup4*p3mup2*p4m
     -   up5*F(1464)-mup4mup5*p1mup1*p4mup2*p4mup3*F(1466)-2*(mup4mup
     -   5*p1mup2*p3mup3*p4mup1*F(1451)+mup4mup5*p1mup2*p3mup1*p4mup3
     -   *F(1462)+mup4mup5*p1mup2*p4mup1*p4mup3*F(1467))-mup4mup5*p1m
     -   up3*p4mup1*p4mup2*F(1468)+mup3mup5*(p2mup1*p2mup2*p3mup4*F(1
     -   447)+p1mup1*(p3mup2*p4mup4*F(1347)+p3mup4*p4mup2*F(1450)+p4m
     -   up2*p4mup4*F(1466))+2*p1mup2*(p3mup4*p4mup1*F(1451)+p3mup1*p
     -   4mup4*F(1462)+p4mup1*p4mup4*F(1467))+p1mup4*(p3mup1*p4mup2*F
     -   (1452)+p3mup2*p4mup1*F(1463)+p4mup1*p4mup2*F(1468)))-p4mup2*
     -   (mup1mup5*p1mup4*p4mup3-mup1mup4*p1mup5*p4mup3-mup1mup5*p1mu
     -   p3*p4mup4+mup1mup3*p1mup5*p4mup4-mup1mup4*p1mup3*p4mup5+mup1
     -   mup3*p1mup4*p4mup5)*F(1469)+mup3mup4*(p2mup1*p2mup2*p3mup5*F
     -   (1447)+p1mup1*(p3mup2*p4mup5*F(1347)+p3mup5*p4mup2*F(1450)+p
     -   4mup2*p4mup5*F(1466))+2*p1mup2*(p3mup5*p4mup1*F(1451)+p3mup1
     -   *p4mup5*F(1462)+p4mup1*p4mup5*F(1467))+p1mup5*(p3mup1*p4mup2
     -   *F(1454)+p3mup2*p4mup1*F(1465)+p4mup1*p4mup2*F(1470)))
       Fa3(30) = mup2mup4*p2mup3*p2mup5*p4mup1*F(1494)+mup2mup3*p2mup1
     -   *p3mup4*p3mup5*F(1499)+mup4mup5*p2mup1*p3mup2*p4mup3*F(1509)
     -   +mup3mup4*p2mup1*p3mup2*p4mup5*F(1511)+mup4mup5*p2mup1*p3mup
     -   3*p4mup2*F(1512)+mup2mup4*p2mup1*p3mup3*p4mup5*F(1514)+mup2m
     -   up3*p2mup1*p3mup4*p4mup5*F(1517)+mup3mup4*p2mup1*p3mup5*p4mu
     -   p2*F(1518)+mup2mup4*p2mup1*p3mup5*p4mup3*F(1519)+mup2mup3*p2
     -   mup1*p3mup5*p4mup4*F(1520)+mup4mup5*p2mup2*p3mup1*p4mup3*F(1
     -   521)+mup3mup4*p2mup2*p3mup1*p4mup5*F(1523)+mup4mup5*p2mup2*p
     -   3mup3*p4mup1*F(1524)+mup2mup5*(p2mup1*(p3mup3*p3mup4*F(1498)
     -   +p3mup3*p4mup4*F(1513)+p3mup4*p4mup3*F(1516)+2*(p2mup4*p3mup
     -   3*F(1472)+p2mup3*p4mup4*F(1486)+p2mup4*p4mup3*F(1532)))+2*(p
     -   2mup3*(p2mup4*p3mup1*F(1478)+p2mup4*p4mup1*F(1491)-p3mup1*p3
     -   mup4*F(1503))-p2mup4*p3mup1*p3mup3*F(1538)))+mup1mup2*p2mup4
     -   *p3mup3*p3mup5*F(1542)+mup2mup4*p2mup3*p2mup5*p3mup1*F(1544)
     -   +mup3mup5*(p2mup1*(p3mup2*p4mup4*F(1510)+p3mup4*p4mup2*F(151
     -   5)-2*(p3mup2*p3mup4*F(1497)+p2mup2*p4mup4*F(1529)+p2mup4*(p3
     -   mup2*F(1471)+p4mup2*F(1531))))+p3mup1*(p2mup2*(2*(p2mup4*F(1
     -   475)+p3mup4*F(1500))+p4mup4*F(1522))+2*p2mup4*p3mup2*F(1537)
     -   )+2*p2mup2*p2mup4*p4mup1*F(1551))+mup3mup4*p2mup1*p3mup2*p3m
     -   up5*F(1556)+mup3mup4*p2mup2*p2mup5*p3mup1*F(1557)+mup3mup4*p
     -   2mup2*p2mup5*p4mup1*F(1558)+mup1mup3*p2mup2*p2mup5*p3mup4*F(
     -   1561)+mup1mup4*p2mup3*p2mup5*p3mup2*F(1565)+mup1mup3*p2mup4*
     -   p2mup5*p4mup2*F(1569)+2*(mup2mup3*p2mup1*p2mup4*p3mup5*F(147
     -   3)+mup1mup5*p2mup2*p2mup3*p3mup4*F(1474)+mup1mup5*p2mup2*p2m
     -   up4*p3mup3*F(1476)-mup1mup3*p2mup2*p2mup4*p3mup5*F(1477)+mup
     -   1mup5*p2mup3*p2mup4*p3mup2*F(1479)+mup1mup2*p2mup3*p2mup4*p3
     -   mup5*F(1480)+mup1mup2*p2mup3*p2mup5*p3mup4*F(1481)+mup2mup3*
     -   p2mup4*p2mup5*p3mup1*F(1482)+mup1mup3*p2mup4*p2mup5*p3mup2*F
     -   (1483)+mup1mup2*p2mup4*p2mup5*p3mup3*F(1484)+mup4mup5*p2mup1
     -   *p2mup3*p4mup2*F(1485)+mup2mup4*p2mup1*p2mup3*p4mup5*F(1487)
     -   +mup4mup5*p2mup2*p2mup3*p4mup1*F(1488)+mup1mup5*p2mup2*p2mup
     -   3*p4mup4*F(1489)-mup1mup4*p2mup2*p2mup3*p4mup5*F(1490)+mup1m
     -   up5*p2mup3*p2mup4*p4mup2*F(1492)+mup1mup2*p2mup3*p2mup4*p4mu
     -   p5*F(1493)+mup1mup4*p2mup3*p2mup5*p4mup2*F(1495)+mup1mup2*p2
     -   mup3*p2mup5*p4mup4*F(1496)+mup1mup5*p2mup2*p3mup3*p3mup4*F(1
     -   501)-mup1mup3*p2mup2*p3mup4*p3mup5*F(1502)+mup1mup5*p2mup3*p
     -   3mup2*p3mup4*F(1504)+mup1mup2*p2mup3*p3mup4*p3mup5*F(1505)+m
     -   up2mup3*p2mup5*p3mup1*p3mup4*F(1506)+mup1mup3*p2mup5*p3mup2*
     -   p3mup4*F(1507)+mup1mup2*p2mup5*p3mup3*p3mup4*F(1508)+mup3mup
     -   4*p2mup1*p2mup5*p3mup2*F(1525)+mup2mup4*p2mup1*p2mup5*p3mup3
     -   *F(1526)+mup2mup3*p2mup1*p2mup5*p3mup4*F(1527)+mup4mup5*p2mu
     -   p1*p2mup2*p4mup3*F(1528)-mup3mup4*p2mup1*p2mup2*p4mup5*F(153
     -   0)+mup2mup3*p2mup1*p2mup4*p4mup5*F(1533)+mup3mup4*p2mup1*p2m
     -   up5*p4mup2*F(1534)-mup2mup4*p2mup1*p2mup5*p4mup3*F(1535)+mup
     -   2mup3*p2mup1*p2mup5*p4mup4*F(1536)+mup2mup3*p2mup4*p3mup1*p3
     -   mup5*F(1539)+mup1mup5*p2mup4*p3mup2*p3mup3*F(1540)+mup1mup3*
     -   p2mup4*p3mup2*p3mup5*F(1541)+mup4mup5*p2mup2*p2mup3*p3mup1*F
     -   (1543)+mup4mup5*p2mup2*p3mup1*p3mup3*F(1545)+mup3mup4*p2mup2
     -   *p3mup1*p3mup5*F(1546)+mup4mup5*p2mup3*p3mup1*p3mup2*F(1547)
     -   +mup2mup4*p2mup3*p3mup1*p3mup5*F(1548)+mup3mup4*p2mup5*p3mup
     -   1*p3mup2*F(1549)+mup2mup4*p2mup5*p3mup1*p3mup3*F(1550)+mup2m
     -   up3*p2mup4*p2mup5*p4mup1*F(1552)-mup1mup5*p2mup2*p2mup4*p4mu
     -   p3*F(1553)+mup1mup2*p2mup4*p2mup5*p4mup3*F(1554)+mup4mup5*p2
     -   mup1*p3mup2*p3mup3*F(1555)-mup1mup4*p2mup2*p2mup3*p3mup5*F(1
     -   559)-mup1mup4*p2mup2*p2mup5*p3mup3*F(1560)-mup1mup3*p2mup2*p
     -   2mup4*p4mup5*F(1562)+mup1mup4*p2mup2*p2mup5*p4mup3*F(1563)-m
     -   up1mup3*p2mup2*p2mup5*p4mup4*F(1564)+mup2mup4*p2mup1*p3mup3*
     -   p3mup5*F(1566)-mup1mup4*p2mup2*p3mup3*p3mup5*F(1567)-mup1mup
     -   4*p2mup3*p3mup2*p3mup5*F(1568)-mup1mup4*p2mup5*p3mup2*p3mup3
     -   *F(1570))
       Return
       End
