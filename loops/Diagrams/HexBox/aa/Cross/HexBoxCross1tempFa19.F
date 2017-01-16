c************** Calling the Fa functions*************************
       subroutine HexBoxCross1Fa19(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,
     &   p3mup4,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa
     &   19)
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
       Complex*16 Fa19(235:247)
       COMMON/HexBoxCross1FaFunctions/Fa
       Fa19(235) = p4mup4*F(1322)+p5mup4*F(1323)+p6mup4*F(1324)+p1mup4
     &   *F(1325)+p2mup4*F(1326)
       Fa19(236) = p4mup4*F(1327)+p5mup4*F(1328)+p6mup4*F(1329)-p2mup4
     &   *F(1330)-p1mup4*F(1331)
       Fa19(237) = -(p5mup4*F(721))+p1mup4*F(1332)+p2mup4*F(1333)-p6mu
     &   p4*F(1334)+p4mup4*F(1335)
       Fa19(238) = p6mup4*F(804)+p5mup4*F(805)-p2mup4*F(807)-p1mup4*F(
     &   814)+p4mup4*F(1336)
       Fa19(239) = p6mup3*F(827)-p2mup3*F(1337)-p1mup3*F(1338)+p4mup3*
     &   F(1339)+p5mup3*F(1340)
       Fa19(240) = p5mup3*F(718)-p1mup3*F(1341)+p4mup3*F(1342)-p2mup3*
     &   F(1343)+p6mup3*F(1344)
       Fa19(241) = p1mup4*F(751)+p6mup4*F(758)+p5mup4*F(759)+p2mup4*F(
     &   1346)+p4mup4*F(1347)
       Fa19(242) = p1mup4*F(1348)+p2mup4*F(1349)+p4mup4*F(1350)+p5mup4
     &   *F(1351)+p6mup4*F(1352)
       Fa19(243) = p4mup4*F(1353)+p5mup4*F(1354)+p6mup4*F(1355)-p2mup4
     &   *F(1356)-p1mup4*F(1357)
       Fa19(244) = p1mup3*F(525)+p4mup3*F(760)-p2mup3*F(1358)+p5mup3*F
     &   (1359)+p6mup3*F(1360)
       Fa19(245) = p5mup3*F(1088)-p1mup3*F(1361)-p2mup3*F(1362)+p4mup3
     &   *F(1363)+p6mup3*F(1364)
       Fa19(246) = p6mup4*F(1372)-4*(p2mup4*F(1369)-p4mup4*F(1370)-p5m
     &   up4*F(1371)+p1mup4*F(1373))
       Fa19(247) = p2mup3*F(597)+p1mup3*F(1171)+p4mup3*F(1181)+p5mup3*
     &   F(1182)+p6mup3*F(1183)
       Return
       End
