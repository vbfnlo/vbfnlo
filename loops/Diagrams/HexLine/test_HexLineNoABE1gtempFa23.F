c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineNoAbe1g23(p1mup2,p1mup3,p1mup4
     -   ,p1mup6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3m
     -   up6,p4mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,
     -   p6mup2,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3m
     -   up4,mup3mup6,mup4mup6,Fa23)
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
       Complex*16 Fa23(309:322)
       COMMON/test_FaFunctionsNoAbe1g/Fa
       Fa23(309) = p4mup2*(p2mup4*p2mup6*p3mup3*F(1294)+p1mup3*p2mup4*
     -   p3mup6*F(1295)+p2mup4*p3mup3*p3mup6*F(1296)+p1mup3*p2mup4*p5
     -   mup6*F(1297)+p2mup4*p3mup3*p5mup6*F(1298)+p1mup3*p2mup4*p4mu
     -   p6*F(1299)+p2mup4*p3mup3*p4mup6*F(1300)+p1mup4*(p1mup3*(p3mu
     -   p6*F(1295)+p5mup6*F(1297)+p4mup6*F(1299))+p2mup3*(p2mup6*F(1
     -   292)+p3mup6*F(1295)+p5mup6*F(1297)+p4mup6*F(1299))+p3mup3*(p
     -   2mup6*F(1294)+p3mup6*F(1296)+p5mup6*F(1298)+p4mup6*F(1300)))
     -   -p2mup6*p3mup3*p3mup4*F(1301)+p1mup3*p3mup4*p3mup6*F(1302)+p
     -   3mup3*p3mup4*p3mup6*F(1303)+p1mup3*p3mup4*p5mup6*F(1304)+p3m
     -   up3*p3mup4*p5mup6*F(1305)+p1mup3*p3mup4*p4mup6*F(1306)+p3mup
     -   3*p3mup4*p4mup6*F(1307)+p1mup3*p2mup6*p4mup4*F(1308)+p1mup3*
     -   p3mup6*p4mup4*F(1309)-p1mup3*p4mup4*p5mup6*F(1310)-p3mup3*p4
     -   mup4*p5mup6*F(1311)+p1mup3*p4mup4*p4mup6*F(1312)+p2mup3*(-(p
     -   2mup6*p3mup4*F(1293))+p2mup4*(p2mup6*F(1292)+p3mup6*F(1295)+
     -   p5mup6*F(1297)+p4mup6*F(1299))+p3mup4*p3mup6*F(1302)+p3mup4*
     -   p5mup6*F(1304)+p3mup4*p4mup6*F(1306)+p2mup6*p4mup4*F(1308)+p
     -   3mup6*p4mup4*F(1309)-p4mup4*p5mup6*F(1310)+p4mup4*p4mup6*F(1
     -   312))+p2mup6*p3mup3*p4mup4*F(1313)+p3mup3*p3mup6*p4mup4*F(13
     -   14)+p3mup3*p4mup4*p4mup6*F(1315))
       Fa23(310) = p1mup6*p2mup3*p2mup4*p4mup2*F(1292)-p1mup6*p2mup3*p
     -   3mup4*p4mup2*F(1293)+p1mup6*p2mup4*p3mup3*p4mup2*F(1294)+p1m
     -   up4*p1mup6*p4mup2*(p2mup3*F(1292)+p3mup3*F(1294))-p1mup6*p3m
     -   up3*p3mup4*p4mup2*F(1301)+p1mup6*p2mup3*p4mup2*p4mup4*F(1308
     -   )+p1mup3*p1mup6*p4mup2*((p1mup4+p2mup4)*F(1292)-p3mup4*F(129
     -   3)+p4mup4*F(1308))+p1mup6*p3mup3*p4mup2*p4mup4*F(1313)+p3mup
     -   4*p4mup2*p4mup3*p4mup6*F(1316)+p3mup4*p4mup2*p4mup3*p5mup6*F
     -   (1317)+p3mup2*p4mup3*p4mup6*p5mup4*F(1318)-p3mup3*p4mup2*p4m
     -   up6*p5mup4*F(1319)+p3mup4*p4mup2*p4mup6*p5mup3*F(1320)+p3mup
     -   2*p4mup3*p5mup4*p5mup6*F(1321)-p3mup3*p4mup2*p5mup4*p5mup6*F
     -   (1322)+p3mup4*p4mup2*p5mup3*p5mup6*F(1323)+p3mup6*p4mup2*p4m
     -   up4*p5mup3*F(1324)+p4mup2*p4mup3*p4mup4*p5mup6*F(1325)+p4mup
     -   2*p4mup4*p4mup6*p5mup3*F(1326)-p4mup2*p4mup4*p5mup3*p5mup6*F
     -   (1327)+p3mup6*p4mup2*p4mup3*p4mup4*F(1328)+p4mup2*p4mup3*p4m
     -   up4*p4mup6*F(1329)+p3mup6*p4mup2*p4mup3*p5mup4*F(1330)+p4mup
     -   2*p4mup3*p4mup6*p5mup4*F(1331)-p4mup2*p4mup3*p5mup4*p5mup6*F
     -   (1332)+p3mup2*p4mup6*p5mup3*p5mup4*F(1333)-p3mup6*p4mup2*p5m
     -   up3*p5mup4*F(1334)-p3mup2*p5mup3*p5mup4*p5mup6*F(1335)+p4mup
     -   2*p4mup6*p5mup3*p5mup4*F(1336)-p4mup2*p5mup3*p5mup4*p5mup6*F
     -   (1337)
       Fa23(311) = p1mup4*p4mup2*(p1mup6*p4mup3*F(1308)+p2mup6*p4mup3*
     -   F(1308)+p3mup6*p4mup3*F(1338)+p1mup6*p5mup3*F(1341)+p2mup6*p
     -   5mup3*F(1341)+p3mup6*p5mup3*F(1344))+p1mup3*p5mup4*(p2mup6*p
     -   3mup2*F(1340)+p3mup2*p3mup6*F(1342)+p3mup2*p4mup6*F(1343)+p2
     -   mup6*p4mup2*F(1345)+p3mup6*p4mup2*F(1346))+p1mup6*(p3mup4*p4
     -   mup2*p4mup3*F(1339)+p2mup4*p4mup2*(p4mup3*F(1308)+p5mup3*F(1
     -   341))+(p1mup3+p2mup3)*p5mup4*(p3mup2*F(1340)+p4mup2*F(1345))
     -   +p3mup2*p3mup3*p5mup4*F(1347))
       Fa23(312) = p4mup2*p4mup3*(p2mup4*p2mup6*F(1308)+(p1mup4+p2mup4
     -   )*p4mup6*F(1312))+p3mup2*p5mup4*(p3mup3*(p3mup6*F(1290)+p2mu
     -   p6*F(1347))+p1mup3*p5mup6*F(1353)+p2mup3*(p2mup6*F(1340)+p3m
     -   up6*F(1342)+p4mup6*F(1343)+p5mup6*F(1353))+(p1mup6+p2mup6)*(
     -   p4mup3*F(1350)-p5mup3*F(1361)))+p4mup2*(p2mup6*p3mup4*p4mup3
     -   *F(1339)+p2mup3*p2mup6*p5mup4*F(1345)+p2mup3*p3mup6*p5mup4*F
     -   (1346)+p1mup4*p4mup3*p5mup6*F(1348)+p3mup4*p3mup6*p4mup3*F(1
     -   349)-p1mup6*p3mup3*p5mup4*F(1351)-p2mup6*p3mup3*p5mup4*F(135
     -   1)+p1mup6*p3mup4*p5mup3*F(1352)+p2mup6*p3mup4*p5mup3*F(1352)
     -   +p1mup4*p4mup6*p5mup3*F(1354)+p1mup4*p5mup3*p5mup6*F(1355)+p
     -   2mup4*(p3mup6*p4mup3*F(1338)+p2mup6*p5mup3*F(1341)+p3mup6*p5
     -   mup3*F(1344)+p4mup3*p5mup6*F(1348)+p4mup6*p5mup3*F(1354)+p5m
     -   up3*p5mup6*F(1355))+p1mup3*p4mup6*p5mup4*F(1356)+p2mup3*p4mu
     -   p6*p5mup4*F(1356)+p1mup6*p4mup4*p5mup3*F(1357)+p2mup6*p4mup4
     -   *p5mup3*F(1357)-p1mup3*p5mup4*p5mup6*F(1358)-p2mup3*p5mup4*p
     -   5mup6*F(1358)+p1mup6*p4mup3*p4mup4*F(1359)+p2mup6*p4mup3*p4m
     -   up4*F(1359)+p1mup6*p4mup3*p5mup4*F(1360)+p2mup6*p4mup3*p5mup
     -   4*F(1360)-(p1mup6+p2mup6)*p5mup3*p5mup4*F(1362))
       Fa23(313) = p3mup6*p4mup2*(-(p3mup3*p5mup4*F(1364))+p3mup4*p5mu
     -   p3*F(1365))+p3mup2*p5mup4*(p3mup3*p5mup6*F(1284)+p3mup6*p4mu
     -   p3*F(1363)+p3mup3*p4mup6*F(1366)+p3mup6*p5mup3*F(1367))
       Fa23(314) = p4mup2*(-(p4mup3*F(1368))-p3mup3*F(1370)+p5mup3*F(1
     -   371)+4*(p1mup3+p2mup3)*F(1372))+p3mup2*(-4*p4mup3*F(1369)-p5
     -   mup3*F(1373)+(p1mup3+p2mup3)*F(1374)+p3mup3*F(1375))
       Fa23(315) = mup4mup6*(p2mup2*p2mup3*F(1226)+2*p4mup2*s12*(p5mup
     -   3*F(1078)+p3mup3*F(1208)+p4mup3*F(1376)))
       Fa23(316) = mup4mup6*(p2mup2*(-(p4mup3*F(1220))+p5mup3*F(1223))
     -   +p1mup3*(p1mup2+p2mup2)*F(1226)+p2mup3*(2*p4mup2*s12*F(1076)
     -   +p1mup2*F(1226))+(p1mup2+p2mup2)*p3mup3*F(1227))
       Fa23(317) = mup4mup6*(2*p1mup3*p4mup2*s12*F(1076)-p1mup2*p4mup3
     -   *F(1220)+p1mup2*p5mup3*F(1223))
       Fa23(318) = mup3mup6*((p1mup2+p2mup2)*(2*(p3mup4*F(1228)-p4mup4
     -   *F(1229))+p5mup4*F(1377))+2*p4mup2*s12*(p5mup4*F(1090)+p4mup
     -   4*F(1091)-p3mup4*F(1378)))
       Fa23(319) = -(mup3mup6*(p1mup4+p2mup4)*(16*p4mup2*s12*F(1089)-(
     -   p1mup2+p2mup2)*F(1379)))
       Fa23(320) = mup3mup4*(p3mup6*p4mup2*s12*F(1073)+p1mup6*p4mup2*s
     -   12*F(1074)+p2mup6*p4mup2*s12*F(1074)+p4mup2*p4mup6*s12*F(108
     -   3)+p4mup2*p5mup6*s12*F(1087)-(p1mup2+p2mup2)*p3mup6*F(1240)-
     -   p1mup2*p5mup6*F(1242)-p2mup2*p5mup6*F(1242)-p1mup2*p4mup6*F(
     -   1380)-p2mup2*p4mup6*F(1380)+(p1mup2+p2mup2)*(p1mup6+p2mup6)*
     -   F(1381))
       Fa23(321) = p2mup2*(p1mup3*p1mup6*((p1mup4+p2mup4)*F(1244)+p3mu
     -   p4*F(1245)+p4mup4*F(1246))+p2mup3*p2mup4*p3mup6*F(1247)+p2mu
     -   p3*p3mup4*p3mup6*F(1248)+p2mup3*p3mup6*p4mup4*F(1249)+p2mup3
     -   *p2mup4*p4mup6*F(1250)-p2mup4*p2mup6*p4mup3*F(1251)+p2mup3*p
     -   3mup4*p4mup6*F(1252)+p2mup3*p4mup4*p4mup6*F(1253)+p2mup3*p2m
     -   up4*p5mup6*F(1254)+p2mup3*p3mup4*p5mup6*F(1255)+p2mup3*p4mup
     -   4*p5mup6*F(1256)+p2mup4*p2mup6*p3mup3*F(1260)+p2mup6*p3mup3*
     -   p3mup4*F(1261)+p2mup6*p3mup4*p4mup3*F(1262)+p2mup4*p3mup3*p4
     -   mup6*F(1263)-p2mup6*p4mup4*p5mup3*F(1273)+p2mup6*p3mup3*p4mu
     -   p4*F(1279)-p2mup6*p4mup3*p4mup4*F(1280)-p2mup4*p2mup6*p5mup3
     -   *F(1285)+p1mup4*(p2mup3*p3mup6*F(1247)+p2mup3*p4mup6*F(1250)
     -   -p2mup6*p4mup3*F(1251)+p2mup3*p5mup6*F(1254)-p3mup6*p4mup3*F
     -   (1257)+p2mup6*p3mup3*F(1260)+p3mup3*p4mup6*F(1263)-p2mup6*p5
     -   mup3*F(1285))+p2mup6*p3mup4*p5mup3*F(1382))+p1mup2*(p1mup3*p
     -   1mup6*((p1mup4+p2mup4)*F(1244)+p3mup4*F(1245)+p4mup4*F(1246)
     -   )+p2mup3*p2mup4*p3mup6*F(1247)+p2mup3*p3mup4*p3mup6*F(1248)+
     -   p2mup3*p3mup6*p4mup4*F(1249)+p2mup3*p2mup4*p4mup6*F(1250)-p2
     -   mup4*p2mup6*p4mup3*F(1251)+p2mup3*p3mup4*p4mup6*F(1252)+p2mu
     -   p3*p4mup4*p4mup6*F(1253)+p2mup3*p2mup4*p5mup6*F(1254)+p2mup3
     -   *p3mup4*p5mup6*F(1255)+p2mup3*p4mup4*p5mup6*F(1256)-p2mup4*p
     -   3mup6*p4mup3*F(1257)+p2mup4*p2mup6*p3mup3*F(1260)+p2mup6*p3m
     -   up3*p3mup4*F(1261)+p2mup6*p3mup4*p4mup3*F(1262)+p2mup4*p3mup
     -   3*p4mup6*F(1263)+p3mup3*p3mup4*p4mup6*F(1264)+p3mup4*p3mup6*
     -   p4mup3*F(1270)-p2mup6*p4mup4*p5mup3*F(1273)+p2mup6*p3mup3*p4
     -   mup4*F(1279)-p2mup6*p4mup3*p4mup4*F(1280)+p3mup6*p4mup3*p4mu
     -   p4*F(1281)-p2mup4*p2mup6*p5mup3*F(1285)+p1mup4*(p2mup3*p3mup
     -   6*F(1247)+p2mup3*p4mup6*F(1250)-p2mup6*p4mup3*F(1251)+p2mup3
     -   *p5mup6*F(1254)-p3mup6*p4mup3*F(1257)+p2mup6*p3mup3*F(1260)+
     -   p3mup3*p4mup6*F(1263)-p2mup6*p5mup3*F(1285))+p3mup3*p4mup4*p
     -   4mup6*F(1287)+p2mup6*p3mup4*p5mup3*F(1382))
       Fa23(322) = p1mup2*(p1mup3*p2mup4*p2mup6*F(1244)+p2mup3*p2mup4*
     -   p2mup6*F(1244)+p1mup3*p2mup6*p3mup4*F(1245)+p2mup3*p2mup6*p3
     -   mup4*F(1245)+p1mup3*p2mup6*p4mup4*F(1246)+p2mup3*p2mup6*p4mu
     -   p4*F(1246)+p1mup3*p2mup4*p3mup6*F(1247)+p1mup3*p3mup4*p3mup6
     -   *F(1248)+p1mup3*p3mup6*p4mup4*F(1249)-p2mup4*p4mup3*p5mup6*F
     -   (1258)-p2mup4*p4mup3*p4mup6*F(1259)-p2mup4*p3mup3*p5mup6*F(1
     -   265)-p3mup3*p3mup4*p5mup6*F(1266)-p3mup3*p4mup4*p5mup6*F(126
     -   7)+p2mup4*p3mup3*p3mup6*F(1268)+p3mup3*p3mup4*p3mup6*F(1269)
     -   +p3mup4*p4mup3*p5mup6*F(1271)+p3mup4*p4mup3*p4mup6*F(1272)-p
     -   3mup6*p4mup4*p5mup3*F(1274)-p4mup4*p4mup6*p5mup3*F(1275)-p2m
     -   up4*p4mup6*p5mup3*F(1276)+p4mup3*p4mup4*p5mup6*F(1278)+p1mup
     -   6*(p2mup3*p2mup4*F(1244)+p2mup3*p3mup4*F(1245)+p2mup3*p4mup4
     -   *F(1246)+p2mup4*p3mup3*F(1260)+p3mup3*p3mup4*F(1261)+p3mup3*
     -   p4mup4*F(1279))+p4mup3*p4mup4*p4mup6*F(1282)+p2mup4*p5mup3*p
     -   5mup6*F(1283)+p3mup4*p5mup3*p5mup6*F(1284)+p3mup3*p3mup6*p4m
     -   up4*F(1288)+p2mup4*p3mup6*p5mup3*F(1289)+p1mup4*(p1mup6*p2mu
     -   p3*F(1244)+p1mup3*p2mup6*F(1244)+p2mup3*p2mup6*F(1244)+p1mup
     -   3*p3mup6*F(1247)-p4mup3*p5mup6*F(1258)-p4mup3*p4mup6*F(1259)
     -   +p1mup6*p3mup3*F(1260)-p3mup3*p5mup6*F(1265)+p3mup3*p3mup6*F
     -   (1268)-p4mup6*p5mup3*F(1276)+p5mup3*p5mup6*F(1283)+p3mup6*p5
     -   mup3*F(1289))+p3mup4*p3mup6*p5mup3*F(1290)-p4mup4*p5mup3*p5m
     -   up6*F(1291)+p3mup4*p4mup6*p5mup3*F(1383))+p2mup2*(p1mup3*p2m
     -   up4*p2mup6*F(1244)+p2mup3*p2mup4*p2mup6*F(1244)+p1mup3*p2mup
     -   6*p3mup4*F(1245)+p2mup3*p2mup6*p3mup4*F(1245)+p1mup3*p2mup6*
     -   p4mup4*F(1246)+p2mup3*p2mup6*p4mup4*F(1246)+p1mup3*p2mup4*p3
     -   mup6*F(1247)+p1mup3*p3mup4*p3mup6*F(1248)+p1mup3*p3mup6*p4mu
     -   p4*F(1249)-p2mup4*p3mup6*p4mup3*F(1257)-p2mup4*p4mup3*p5mup6
     -   *F(1258)-p2mup4*p4mup3*p4mup6*F(1259)+p1mup6*(p2mup3*p2mup4*
     -   F(1244)+p2mup3*p3mup4*F(1245)+p2mup3*p4mup4*F(1246)+p2mup4*p
     -   3mup3*F(1260)+p3mup3*p3mup4*F(1261))+p3mup3*p3mup4*p4mup6*F(
     -   1264)-p2mup4*p3mup3*p5mup6*F(1265)-p3mup3*p3mup4*p5mup6*F(12
     -   66)-p3mup3*p4mup4*p5mup6*F(1267)+p2mup4*p3mup3*p3mup6*F(1268
     -   )+p3mup3*p3mup4*p3mup6*F(1269)+p3mup4*p3mup6*p4mup3*F(1270)+
     -   p3mup4*p4mup3*p5mup6*F(1271)+p3mup4*p4mup3*p4mup6*F(1272)-p3
     -   mup6*p4mup4*p5mup3*F(1274)-p4mup4*p4mup6*p5mup3*F(1275)-p2mu
     -   p4*p4mup6*p5mup3*F(1276)+p4mup3*p4mup4*p5mup6*F(1278)+p3mup6
     -   *p4mup3*p4mup4*F(1281)+p4mup3*p4mup4*p4mup6*F(1282)+p2mup4*p
     -   5mup3*p5mup6*F(1283)+p3mup4*p5mup3*p5mup6*F(1284)+p3mup3*p4m
     -   up4*p4mup6*F(1287)+p3mup3*p3mup6*p4mup4*F(1288)+p2mup4*p3mup
     -   6*p5mup3*F(1289)+p1mup4*(p1mup6*p2mup3*F(1244)+p1mup3*p2mup6
     -   *F(1244)+p2mup3*p2mup6*F(1244)+p1mup3*p3mup6*F(1247)-p4mup3*
     -   p5mup6*F(1258)-p4mup3*p4mup6*F(1259)+p1mup6*p3mup3*F(1260)-p
     -   3mup3*p5mup6*F(1265)+p3mup3*p3mup6*F(1268)-p4mup6*p5mup3*F(1
     -   276)+p5mup3*p5mup6*F(1283)+p3mup6*p5mup3*F(1289))+p3mup4*p3m
     -   up6*p5mup3*F(1290)-p4mup4*p5mup3*p5mup6*F(1291)+p3mup4*p4mup
     -   6*p5mup3*F(1383))
       Return
       End
