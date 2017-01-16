c************** Calling the Fa functions*************************
       subroutine HexBox2Fa18(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup
     &   6,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa18)
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
       Complex*16 Fa18(324:342)
       COMMON/HexBox2FaFunctions/Fa
       Fa18(324) = -(p2mup6*F(1288))+p4mup6*F(1289)+p5mup6*F(1290)+p6m
     &   up6*F(1291)+p1mup6*F(1292)
       Fa18(325) = -(p2mup6*F(1293))+p4mup6*F(1294)+p5mup6*F(1295)+p6m
     &   up6*F(1296)+p1mup6*F(1297)
       Fa18(326) = -(p2mup6*F(1298))+p4mup6*F(1299)+p5mup6*F(1300)+p6m
     &   up6*F(1301)-p1mup6*F(1302)
       Fa18(327) = p2mup6*F(1303)+p4mup6*F(1304)+p5mup6*F(1305)+p6mup6
     &   *F(1306)-p1mup6*F(1307)
       Fa18(328) = p2mup6*F(1308)+p4mup6*F(1309)+p5mup6*F(1310)+p6mup6
     &   *F(1311)-p1mup6*F(1312)
       Fa18(329) = -(p1mup6*F(1319))+p4mup6*F(1320)+p5mup6*F(1321)+p6m
     &   up6*F(1322)-p2mup6*F(1323)
       Fa18(330) = p6mup3*F(1325)+4*(p2mup3*F(1324)+p4mup3*F(1326)+p1m
     &   up3*F(1327)+p5mup3*F(1328))
       Fa18(331) = p1mup3*(p2mup6*F(1329)+p5mup6*F(1330)+p6mup6*F(1331
     &   )+p1mup6*F(1332)+p4mup6*F(1333))-p1mup6*p2mup3*F(1334)+p1mup
     &   6*p5mup3*F(1335)+p1mup6*p6mup3*F(1336)+p1mup6*p4mup3*F(1337)
     &   -p2mup3*p2mup6*F(1338)+p2mup6*p4mup3*F(1339)+p2mup6*p5mup3*F
     &   (1340)+p2mup6*p6mup3*F(1341)+p2mup3*p5mup6*F(1342)+p2mup3*p6
     &   mup6*F(1343)+p2mup3*p4mup6*F(1344)-p4mup6*p5mup3*F(1345)-p4m
     &   up3*p5mup6*F(1345)-p4mup6*p6mup3*F(1346)-p4mup3*p6mup6*F(134
     &   6)-p5mup3*p5mup6*F(1347)-p5mup6*p6mup3*F(1348)-p5mup3*p6mup6
     &   *F(1348)-p6mup3*p6mup6*F(1349)-p4mup3*p4mup6*F(1350)-mup3mup
     &   6*F(1351)
       Fa18(332) = p4mup6*p5mup3*F(1071)+p4mup3*p5mup6*F(1071)+p5mup3*
     &   p5mup6*F(1072)+p5mup6*p6mup3*F(1073)+p5mup3*p6mup6*F(1073)+p
     &   4mup6*p6mup3*F(1210)+p4mup3*p6mup6*F(1210)+p4mup3*p4mup6*F(1
     &   211)-p1mup3*p5mup6*F(1278)-p1mup3*p6mup6*F(1279)-p1mup3*p4mu
     &   p6*F(1280)-p2mup3*p5mup6*F(1352)-p2mup3*p6mup6*F(1353)-p2mup
     &   3*p4mup6*F(1354)+mup3mup6*F(1355)-p1mup6*(p5mup3*F(941)+p6mu
     &   p3*F(1192)+p4mup3*F(1193)-p1mup3*F(1286)-p2mup3*F(1356))+p6m
     &   up3*p6mup6*F(1359)-p2mup6*(p5mup3*F(945)-p2mup3*F(1357)-p1mu
     &   p3*F(1358)+p4mup3*F(1360)+p6mup3*F(1361))
       Fa18(333) = p6mup3*F(1230)-p2mup3*F(1365)+p1mup3*F(1366)+p4mup3
     &   *F(1367)+p5mup3*F(1368)
       Fa18(334) = p5mup6*F(939)-p1mup6*F(1371)-p2mup6*F(1372)+p4mup6*
     &   F(1373)+p6mup6*F(1374)
       Fa18(335) = p1mup6*F(1375)+p2mup6*F(1376)+p4mup6*F(1377)+p5mup6
     &   *F(1378)+p6mup6*F(1379)
       Fa18(336) = p5mup3*p5mup6*F(940)+p5mup6*p6mup3*F(941)+p5mup3*p6
     &   mup6*F(941)+p4mup6*p5mup3*F(942)+p4mup3*p5mup6*F(942)-p1mup6
     &   *p4mup3*F(1103)-p1mup6*p6mup3*F(1187)-p2mup6*p6mup3*F(1191)+
     &   p6mup3*p6mup6*F(1192)+p4mup6*p6mup3*F(1193)+p4mup3*p6mup6*F(
     &   1193)+p4mup3*p4mup6*F(1252)-p1mup3*(p5mup6*F(930)+p4mup6*F(1
     &   243)+p6mup6*F(1286)-p2mup6*F(1380)+p1mup6*F(1381))-p2mup3*(p
     &   5mup6*F(932)+p4mup6*F(1248)+p6mup6*F(1358)-p1mup6*F(1380)-p2
     &   mup6*F(1382))+mup3mup6*F(1383)-p1mup6*p5mup3*F(1384)-p2mup6*
     &   p4mup3*F(1385)-p2mup6*p5mup3*F(1386)
       Fa18(337) = p4mup6*F(1258)+p5mup6*F(1387)+p6mup6*F(1388)+p2mup6
     &   *F(1389)+p1mup6*F(1390)
       Fa18(338) = p2mup6*F(1393)+2*(p4mup6*F(1263)+p5mup6*F(1391)+p6m
     &   up6*F(1392)+p1mup6*F(1394))
       Fa18(339) = p4mup6*F(1395)+p5mup6*F(1396)+p6mup6*F(1397)+p2mup6
     &   *F(1398)+p1mup6*F(1399)
       Fa18(340) = p2mup6*F(1400)+p4mup6*F(1401)+p5mup6*F(1402)-p6mup6
     &   *F(1403)+p1mup6*F(1404)
       Fa18(341) = p2mup6*F(1405)-p4mup6*F(1406)+p5mup6*F(1407)-p6mup6
     &   *F(1408)+p1mup6*F(1409)
       Fa18(342) = -(p5mup6*F(958))+p2mup6*F(1411)+p1mup6*F(1412)-p4mu
     &   p6*F(1413)+p6mup6*F(1414)
       Return
       End
