c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineAbe24(p1mup2,p1mup3,p1mup4,p1m
     -   up5,p2mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,
     -   p4mup2,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mu
     -   p2,p6mup3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,
     -   mup3mup5,mup4mup5,Fa24)
       IMPLICIT NONE
      Complex*16   p1mup2, p1mup3, p1mup4, p1mup5, p2mup2, p2mup3, 
     -          p2mup4, p2mup5, p3mup2, p3mup3, p3mup4, p3mup5, 
     -          p4mup2, p4mup3, p4mup4, p4mup5, p5mup2, p5mup3, 
     -          p5mup4, p5mup5, p6mup2, p6mup3, p6mup4, p6mup5
       Complex*16   mup2mup3, mup2mup4, mup2mup5, mup3mup4, mup3mup5, 
     -          mup4mup5
        common/test_FhlFunctionsAbe/F
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
      Complex*16 Fa(445),F(5930)
      Real*8 P(11674) 
        COMMON/test_PFunctionsAbe/P
       Complex*16 Fa24(254:264)
       COMMON/test_FaFunctionsAbe/Fa
       Fa24(254) = (p1mup2+p2mup2)*(mup3mup4*F(1328)-2*(-(p3mup3*p6mup
     -   4*F(1330))-p4mup3*p6mup4*F(1333)+p4mup4*p6mup3*F(1334)+p1mup
     -   4*p6mup3*F(1335)+p2mup4*p6mup3*F(1335)+p3mup4*p6mup3*F(1336)
     -   -p6mup3*p6mup4*F(1337)-p3mup3*p3mup4*F(1338)+p1mup3*(-(p6mup
     -   4*F(1329))+(p1mup4+p2mup4)*F(1331)+p3mup4*F(1332)+p4mup4*F(1
     -   339))+p2mup3*(-(p6mup4*F(1329))+(p1mup4+p2mup4)*F(1331)+p3mu
     -   p4*F(1332)+p4mup4*F(1339))-p1mup4*p4mup3*F(1340)-p2mup4*p4mu
     -   p3*F(1340)-p3mup4*p4mup3*F(1341)-p1mup4*p3mup3*F(1342)-p2mup
     -   4*p3mup3*F(1342)+p3mup3*p4mup4*F(1343)-p4mup3*p4mup4*F(1344)
     -   ))
       Fa24(255) = 2*(p6mup3*(p4mup4*F(1334)+(p1mup4+p2mup4)*F(1335)+p
     -   3mup4*F(1336))-p6mup3*p6mup4*F(1337)+p1mup3*(-(p6mup4*F(1329
     -   ))+(p1mup4+p2mup4)*F(1331)+p3mup4*F(1332)+p4mup4*F(1339))+p2
     -   mup3*(-(p6mup4*F(1329))+(p1mup4+p2mup4)*F(1331)+p3mup4*F(133
     -   2)+p4mup4*F(1339))-p3mup3*(p6mup4*F(1330)+p3mup4*F(1338)+(p1
     -   mup4+p2mup4)*F(1342)-p4mup4*F(1343))-p4mup3*(p6mup4*F(1333)+
     -   (p1mup4+p2mup4)*F(1340)+p3mup4*F(1341)+p4mup4*F(1344)))+mup3
     -   mup4*F(1345)
       Fa24(256) = p6mup3*F(813)+3*((p1mup3+p2mup3)*F(814)-p3mup3*F(81
     -   5)+p4mup3*F(816))
       Fa24(257) = p6mup4*F(912)+3*((p1mup4+p2mup4)*F(913)+p3mup4*F(91
     -   4)+p4mup4*F(915))
       Fa24(258) = p1mup2*F(1347)+p2mup2*F(1348)+4*(p4mup2*F(1346)-p6m
     -   up2*F(1349))+p3mup2*F(1350)
       Fa24(259) = p1mup3*F(1351)+p2mup3*F(1352)+p4mup3*F(1353)+4*p6mu
     -   p3*F(1354)+p3mup3*F(1355)
       Fa24(260) = -4*mup2mup3*F(1356)+p1mup2*(p1mup3*F(1372)+8*(p4mup
     -   3*F(1357)+p2mup3*F(1364)+p6mup3*F(1375)+p3mup3*F(1378)))+8*(
     -   p1mup3*p4mup2*F(1359)+p2mup3*p4mup2*F(1360)-p4mup2*p6mup3*F(
     -   1361)+p4mup3*p6mup2*F(1362)-p3mup2*p4mup3*F(1365)+p4mup2*p4m
     -   up3*F(1367)+p1mup3*p3mup2*F(1369)+p2mup3*p3mup2*F(1370)+p3mu
     -   p2*p6mup3*F(1371)+p2mup2*(p4mup3*F(1358)+p1mup3*F(1363)+p6mu
     -   p3*F(1366)-p3mup3*F(1368)+p2mup3*F(1373))-p2mup3*p6mup2*F(13
     -   74)+p1mup3*p6mup2*F(1376)-p6mup2*p6mup3*F(1377)-p3mup2*p3mup
     -   3*F(1379)+p3mup3*p4mup2*F(1380)+p3mup3*p6mup2*F(1381))
       Fa24(261) = 4*mup2mup3*F(1382)+p3mup3*p4mup2*F(1385)+p1mup2*(p1
     -   mup3*F(1399)-8*(p4mup3*F(1387)-p3mup3*F(1392)+p6mup3*F(1395)
     -   -p2mup3*F(1400)))+8*(p1mup3*p2mup2*F(1386)+p1mup3*p4mup2*F(1
     -   388)-p2mup2*p6mup3*F(1389)-p4mup2*p6mup3*F(1390)-p2mup3*(p4m
     -   up2*F(1384)+p2mup2*F(1397))-p4mup3*(p6mup2*F(1391)+p4mup2*F(
     -   1394)-p2mup2*F(1398))-p2mup2*p3mup3*F(1401)+p6mup2*(p1mup3*F
     -   (1396)+p3mup3*F(1404)+p2mup3*F(1405)-p6mup3*F(1406))+p3mup2*
     -   (p4mup3*F(1383)-p1mup3*F(1393)+p2mup3*F(1402)+p6mup3*F(1403)
     -   +p3mup3*F(1407)))
       Fa24(262) = p4mup3*F(1408)-p1mup3*F(1409)+p6mup3*F(1410)+p3mup3
     -   *F(1411)-p2mup3*F(1412)
       Fa24(263) = -16*p2mup2*(p2mup4*p3mup3*F(1418)+p3mup3*p3mup4*F(1
     -   419)+p3mup3*p4mup4*F(1420)+p1mup3*p3mup4*F(1423)+p1mup3*p4mu
     -   p4*F(1424)+p2mup3*(p4mup4*F(1417)+p1mup4*F(1425)+p2mup4*F(14
     -   26)+p3mup4*F(1427))-p1mup4*p3mup3*F(1428)-p1mup3*p1mup4*F(14
     -   32)+p1mup3*p2mup4*F(1433))+p1mup2*(-16*(p3mup3*(p2mup4*F(141
     -   4)+p3mup4*F(1415))+p3mup3*p4mup4*F(1416)-p1mup3*p4mup4*F(142
     -   9))+p1mup4*p3mup3*F(1430)+p2mup3*(p2mup4*F(1421)+p3mup4*F(14
     -   22)-16*(p4mup4*F(1413)-p1mup4*F(1431)))+p1mup3*p2mup4*F(1434
     -   )+p1mup3*p3mup4*F(1435)+p1mup3*p1mup4*F(1436))
       Fa24(264) = p6mup4*(p1mup2*p1mup3*F(1452)+p2mup3*(p1mup2*F(1441
     -   )-16*(p4mup2*F(1437)-p3mup2*F(1445)+p6mup2*F(1446)-p2mup2*F(
     -   1453)))+16*(p4mup2*p4mup3*F(1440)+p1mup3*p2mup2*F(1443)+p3mu
     -   p3*(-(p4mup2*F(1439))+p1mup2*F(1442)+p2mup2*F(1444)-p6mup2*F
     -   (1449))+p1mup2*p4mup3*F(1450)+p2mup2*p4mup3*F(1451)+p3mup2*(
     -   p4mup3*F(1438)+p1mup3*F(1447)+p6mup3*F(1448)+p3mup3*F(1454))
     -   -p2mup2*p6mup3*F(1455)))
       Return
       End
