c************** Calling the Fa functions*************************
       subroutine HexBoxCross1Fa20(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,
     &   p3mup4,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa
     &   20)
       IMPLICIT NONE
      Complex*16   p1mup3, p1mup4, p2mup3, 
     -          p2mup4, p3mup3, p3mup4, 
     -          p4mup3, p4mup4, p5mup3, 
     -          p5mup4, p6mup3, p6mup4
       Complex*16   mup3mup4

        common/HexBoxCross1FhlFunctions/F
       COMMON/HexBoxCross1Invariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12
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
      Complex*16 Fa(272),F(1461)
      Real*8 P(1825) 
        COMMON/HexBoxCross1PFunctions/P
       Complex*16 Fa20(248:272)
       COMMON/HexBoxCross1FaFunctions/Fa
       Fa20(248) = -(p1mup3*F(1374))-p2mup3*F(1375)+p4mup3*F(1376)+p5m
     &   up3*F(1377)+p6mup3*F(1378)
       Fa20(249) = p1mup3*F(525)+p2mup3*F(530)+p4mup3*F(760)+p5mup3*F(
     &   1359)+p6mup3*F(1360)
       Fa20(250) = -(p1mup3*F(755))+p6mup3*F(826)+p4mup3*F(1064)+p5mup
     &   3*F(1339)-p2mup3*F(1361)
       Fa20(251) = -(p1mup3*F(1379))+p4mup3*F(1380)+p5mup3*F(1381)+p6m
     &   up3*F(1382)-p2mup3*F(1383)
       Fa20(252) = 4*(p5mup3*F(728)+p6mup3*F(1069)+p2mup3*F(1364)+p1mu
     &   p3*F(1384))+p4mup3*F(1385)
       Fa20(253) = p5mup3*F(1387)+4*(p4mup3*F(728)+p6mup3*F(813)+p1mup
     &   3*F(1386)+p2mup3*F(1388))
       Fa20(254) = p2mup3*F(600)-p1mup3*F(1169)+p5mup3*F(1175)+p4mup3*
     &   F(1177)+p6mup3*F(1179)
       Fa20(255) = 4*(p2mup3*F(1393)-p4mup3*F(1394)-p5mup3*F(1395)+p1m
     &   up3*F(1396))+p6mup3*F(1397)
       Fa20(256) = p4mup4*F(1398)+p5mup4*F(1399)+p6mup4*F(1400)-p1mup4
     &   *F(1401)-p2mup4*F(1402)
       Fa20(257) = p1mup4*F(1403)+p2mup4*F(1404)-p4mup4*F(1405)+p5mup4
     &   *F(1406)+p6mup4*F(1407)
       Fa20(258) = p1mup4*F(520)+p5mup4*F(522)+p6mup4*F(526)+p4mup4*F(
     &   1408)-p2mup4*F(1409)
       Fa20(259) = p4mup4*F(1410)+p5mup4*F(1411)+p6mup4*F(1412)-p1mup4
     &   *F(1413)-p2mup4*F(1414)
       Fa20(260) = p5mup4*F(517)+p6mup4*F(519)+p1mup4*F(527)+p4mup4*F(
     &   1415)-p2mup4*F(1416)
       Fa20(261) = -(p1mup4*F(1417))+p5mup4*F(1418)+p6mup4*F(1419)-p2m
     &   up4*F(1420)+p4mup4*F(1421)
       Fa20(262) = p1mup4*F(752)+p5mup4*F(753)+p2mup4*F(1408)-p4mup4*F
     &   (1424)-p6mup4*F(1425)
       Fa20(263) = p2mup4*F(1426)+p4mup4*F(1427)+p5mup4*F(1428)+p6mup4
     &   *F(1429)+p1mup4*F(1430)
       Fa20(264) = p5mup4*F(745)+p6mup4*F(747)-p1mup4*F(757)+p4mup4*F(
     &   1431)-p2mup4*F(1432)
       Fa20(265) = p4mup4*F(1433)+p5mup4*F(1434)+p6mup4*F(1435)-p2mup4
     &   *F(1436)-p1mup4*F(1437)
       Fa20(266) = p6mup4*F(1035)+p4mup4*F(1036)-p2mup4*F(1421)-p1mup4
     &   *F(1438)+p5mup4*F(1439)
       Fa20(267) = p5mup4*F(725)+p6mup4*F(979)-p1mup4*F(1440)-p2mup4*F
     &   (1441)+p4mup4*F(1442)
       Fa20(268) = -(p2mup4*F(898))+p4mup4*F(1444)+p5mup4*F(1445)+p6mu
     &   p4*F(1446)-p1mup4*F(1447)
       Fa20(269) = -(p1mup4*F(552))+p4mup4*F(554)+p5mup4*F(555)-p2mup4
     &   *F(898)+p6mup4*F(1460)
       Fa20(270) = -(p1mup3*F(552))+p4mup3*F(554)+p5mup3*F(555)-p2mup3
     &   *F(898)+p6mup3*F(1460)
       Fa20(271) = 2*(p5mup3*F(728)+p6mup3*F(1069)+p2mup3*F(1364)+p1mu
     &   p3*F(1384))+p4mup3*F(1448)
       Fa20(272) = 2*(p4mup3*F(728)+p6mup3*F(813)+p1mup3*F(1386)+p2mup
     &   3*F(1388))+p5mup3*F(1450)
       Return
       End
