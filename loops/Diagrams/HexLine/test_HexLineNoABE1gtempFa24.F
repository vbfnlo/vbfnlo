c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineNoAbe1g24(p1mup2,p1mup3,p1mup4
     -   ,p1mup6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3m
     -   up6,p4mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,
     -   p6mup2,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3m
     -   up4,mup3mup6,mup4mup6,Fa24)
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
       Complex*16 Fa24(323:336)
       COMMON/test_FaFunctionsNoAbe1g/Fa
       Fa24(323) = p1mup2*(p1mup3*(p1mup4*p4mup6*F(1250)+p2mup4*p4mup6
     -   *F(1250)+p3mup4*p4mup6*F(1252)+p4mup4*p4mup6*F(1253)+p1mup4*
     -   p5mup6*F(1254)+p2mup4*p5mup6*F(1254)+p3mup4*p5mup6*F(1255)+p
     -   4mup4*p5mup6*F(1256))-p1mup6*(p1mup4*p4mup3*F(1251)+p2mup4*p
     -   4mup3*F(1251)-p3mup4*p4mup3*F(1262)+p4mup4*p5mup3*F(1273)+p4
     -   mup3*p4mup4*F(1280)+p1mup4*p5mup3*F(1285)+p2mup4*p5mup3*F(12
     -   85)-p3mup4*p5mup3*F(1382)))+p2mup2*(p1mup3*(p1mup4*p4mup6*F(
     -   1250)+p2mup4*p4mup6*F(1250)+p3mup4*p4mup6*F(1252)+p4mup4*p4m
     -   up6*F(1253)+p1mup4*p5mup6*F(1254)+p2mup4*p5mup6*F(1254)+p3mu
     -   p4*p5mup6*F(1255)+p4mup4*p5mup6*F(1256))-p1mup6*(p1mup4*p4mu
     -   p3*F(1251)+p2mup4*p4mup3*F(1251)-p3mup4*p4mup3*F(1262)+p4mup
     -   4*p5mup3*F(1273)-p3mup3*p4mup4*F(1279)+p4mup3*p4mup4*F(1280)
     -   +p1mup4*p5mup3*F(1285)+p2mup4*p5mup3*F(1285)-p3mup4*p5mup3*F
     -   (1382)))
       Fa24(324) = p5mup4*(16*(p2mup2*(p1mup6*(p2mup3*F(1340)+p3mup3*F
     -   (1347)+p4mup3*F(1350))+p1mup3*((p1mup6+p2mup6)*F(1340)+p3mup
     -   6*F(1342)+p4mup6*F(1343)+p5mup6*F(1353)))+p1mup2*((p1mup3+p2
     -   mup3)*(p1mup6+p2mup6)*F(1340)+(p1mup6+p2mup6)*p4mup3*F(1350)
     -   +(p1mup3+p2mup3)*(p3mup6*F(1342)+p4mup6*F(1343)+p5mup6*F(135
     -   3))+p3mup6*p4mup3*F(1363)+p3mup3*(p3mup6*F(1290)+(p1mup6+p2m
     -   up6)*F(1347)+p4mup6*F(1366))))+(p1mup6*p2mup2+p1mup2*(p1mup6
     -   +p2mup6))*p5mup3*F(1384))
       Fa24(325) = p5mup4*(16*(p1mup2*(p3mup3*p5mup6*F(1284)+p4mup3*p4
     -   mup6*F(1318)+p4mup3*p5mup6*F(1321)+p4mup6*p5mup3*F(1333)-p5m
     -   up3*p5mup6*F(1335)+p3mup6*p5mup3*F(1367))+p2mup2*(p4mup6*p5m
     -   up3*F(1333)-p5mup3*p5mup6*F(1335)+p2mup3*p2mup6*F(1340)+p2mu
     -   p3*p3mup6*F(1342)+p2mup3*p4mup6*F(1343)+p2mup3*p5mup6*F(1353
     -   )+p4mup3*(p4mup6*F(1318)+p5mup6*F(1321)+p2mup6*F(1350)+p3mup
     -   6*F(1363))+p3mup3*(p5mup6*F(1284)+p3mup6*F(1290)+p2mup6*F(13
     -   47)+p4mup6*F(1366))+p3mup6*p5mup3*F(1367)))+p2mup2*p2mup6*p5
     -   mup3*F(1384))
       Fa24(326) = p1mup2*(p4mup3*F(1385)+(p1mup3+p2mup3)*F(1386)+p3mu
     -   p3*F(1387)+p5mup3*F(1388))+p2mup2*(p4mup3*F(1385)+(p1mup3+p2
     -   mup3)*F(1386)+p3mup3*F(1387)+p5mup3*F(1388))-4*p4mup2*s12*(p
     -   4mup3*F(1389)+p3mup3*F(1390)+p5mup3*F(1391)+(p1mup3+p2mup3)*
     -   F(1392))
       Fa24(327) = (mup4mup6*(p2mup2*p2mup3*F(1393)+16*p3mup2*s12*(p4m
     -   up3*F(1394)+p5mup3*F(1395))))/s12
       Fa24(328) = mup4mup6*(-16*(p2mup3*p3mup2*s12*F(1076)-p3mup2*p3m
     -   up3*s12*F(1077))+p1mup3*(p1mup2+p2mup2)*F(1393)+p1mup2*p2mup
     -   3*F(1393)+p1mup2*p3mup3*F(1396)+p2mup2*p3mup3*F(1396)+8*(p2m
     -   up2*p4mup3*F(1397)+p2mup2*p5mup3*F(1398)))
       Fa24(329) = mup4mup6*(-2*p1mup3*p3mup2*s12*F(1076)+p1mup2*p4mup
     -   3*F(1397)+p1mup2*p5mup3*F(1398))
       Fa24(330) = mup3mup6*(16*p3mup2*s12*(p5mup4*F(1193)+p4mup4*F(11
     -   94)+p3mup4*F(1378))+(p1mup2+p2mup2)*(p3mup4*F(1399)+p5mup4*F
     -   (1400)+p4mup4*F(1401)))
       Fa24(331) = mup3mup6*(p1mup4+p2mup4)*(16*p3mup2*s12*F(1089)+(p1
     -   mup2+p2mup2)*F(1402))
       Fa24(332) = mup3mup4*(8*p3mup2*s12*(p3mup6*F(1073)+(p1mup6+p2mu
     -   p6)*F(1074)+p4mup6*F(1083)+p5mup6*F(1087))+(p1mup2+p2mup2)*(
     -   p3mup6*F(1403)+p4mup6*F(1404)+p5mup6*F(1405)-(p1mup6+p2mup6)
     -   *F(1406)))
       Fa24(333) = p2mup2*(p2mup4*p2mup6*p3mup3*F(1294)+p2mup3*p2mup4*
     -   p3mup6*F(1295)+p2mup3*p2mup4*p5mup6*F(1297)+p2mup3*p2mup4*p4
     -   mup6*F(1299)+p2mup4*p3mup3*p4mup6*F(1300)-p2mup6*p3mup3*p3mu
     -   p4*F(1301)+p2mup3*p3mup4*p3mup6*F(1302)+p2mup3*p3mup4*p5mup6
     -   *F(1304)+p2mup3*p3mup4*p4mup6*F(1306)+p2mup4*p2mup6*p4mup3*F
     -   (1308)+p1mup3*p1mup6*((p1mup4+p2mup4)*F(1292)-p3mup4*F(1293)
     -   +p4mup4*F(1308))+p2mup3*p3mup6*p4mup4*F(1309)-p2mup3*p4mup4*
     -   p5mup6*F(1310)+p2mup3*p4mup4*p4mup6*F(1312)+p2mup6*p3mup3*p4
     -   mup4*F(1313)+p2mup6*p3mup4*p4mup3*F(1339)+p2mup4*p2mup6*p5mu
     -   p3*F(1341)+p1mup4*(p2mup6*p3mup3*F(1294)+p2mup3*p3mup6*F(129
     -   5)+p2mup3*p5mup6*F(1297)+p2mup3*p4mup6*F(1299)+p3mup3*p4mup6
     -   *F(1300)+p2mup6*p4mup3*F(1308)+p3mup6*p4mup3*F(1338)+p2mup6*
     -   p5mup3*F(1341))+p2mup6*p3mup4*p5mup3*F(1352)+p2mup6*p4mup4*p
     -   5mup3*F(1357)+p2mup6*p4mup3*p4mup4*F(1359))+p1mup2*(p2mup4*p
     -   2mup6*p3mup3*F(1294)+p2mup3*p2mup4*p3mup6*F(1295)+p2mup3*p2m
     -   up4*p5mup6*F(1297)+p2mup3*p2mup4*p4mup6*F(1299)+p2mup4*p3mup
     -   3*p4mup6*F(1300)-p2mup6*p3mup3*p3mup4*F(1301)+p2mup3*p3mup4*
     -   p3mup6*F(1302)+p2mup3*p3mup4*p5mup6*F(1304)+p2mup3*p3mup4*p4
     -   mup6*F(1306)+p3mup3*p3mup4*p4mup6*F(1307)+p2mup4*p2mup6*p4mu
     -   p3*F(1308)+p1mup3*p1mup6*((p1mup4+p2mup4)*F(1292)-p3mup4*F(1
     -   293)+p4mup4*F(1308))+p2mup3*p3mup6*p4mup4*F(1309)-p2mup3*p4m
     -   up4*p5mup6*F(1310)+p2mup3*p4mup4*p4mup6*F(1312)+p2mup6*p3mup
     -   3*p4mup4*F(1313)+p3mup3*p4mup4*p4mup6*F(1315)+p3mup6*p4mup3*
     -   p4mup4*F(1328)+p2mup4*p3mup6*p4mup3*F(1338)+p2mup6*p3mup4*p4
     -   mup3*F(1339)+p2mup4*p2mup6*p5mup3*F(1341)+p1mup4*(p2mup6*p3m
     -   up3*F(1294)+p2mup3*p3mup6*F(1295)+p2mup3*p5mup6*F(1297)+p2mu
     -   p3*p4mup6*F(1299)+p3mup3*p4mup6*F(1300)+p2mup6*p4mup3*F(1308
     -   )+p3mup6*p4mup3*F(1338)+p2mup6*p5mup3*F(1341))+p3mup4*p3mup6
     -   *p4mup3*F(1349)+p2mup6*p3mup4*p5mup3*F(1352)+p2mup6*p4mup4*p
     -   5mup3*F(1357)+p2mup6*p4mup3*p4mup4*F(1359))
       Fa24(334) = p1mup2*(p1mup3*p2mup4*p2mup6*F(1292)+p2mup3*p2mup4*
     -   p2mup6*F(1292)-p1mup3*p2mup6*p3mup4*F(1293)-p2mup3*p2mup6*p3
     -   mup4*F(1293)+p1mup3*p2mup4*p3mup6*F(1295)+p2mup4*p3mup3*p3mu
     -   p6*F(1296)+p2mup4*p3mup3*p5mup6*F(1298)+p1mup3*p3mup4*p3mup6
     -   *F(1302)+p3mup3*p3mup4*p3mup6*F(1303)+p3mup3*p3mup4*p5mup6*F
     -   (1305)+p1mup3*p2mup6*p4mup4*F(1308)+p2mup3*p2mup6*p4mup4*F(1
     -   308)+p1mup3*p3mup6*p4mup4*F(1309)-p3mup3*p4mup4*p5mup6*F(131
     -   1)+p2mup4*p4mup3*p4mup6*F(1312)+p1mup6*(p2mup3*p2mup4*F(1292
     -   )-p2mup3*p3mup4*F(1293)+p2mup4*p3mup3*F(1294)-p3mup3*p3mup4*
     -   F(1301)+p2mup3*p4mup4*F(1308)+p3mup3*p4mup4*F(1313))+p3mup3*
     -   p3mup6*p4mup4*F(1314)+p3mup4*p4mup3*p4mup6*F(1316)+p3mup4*p4
     -   mup3*p5mup6*F(1317)+p3mup4*p4mup6*p5mup3*F(1320)+p3mup4*p5mu
     -   p3*p5mup6*F(1323)+p3mup6*p4mup4*p5mup3*F(1324)+p4mup3*p4mup4
     -   *p5mup6*F(1325)+p4mup4*p4mup6*p5mup3*F(1326)-p4mup4*p5mup3*p
     -   5mup6*F(1327)+p4mup3*p4mup4*p4mup6*F(1329)+p2mup4*p3mup6*p5m
     -   up3*F(1344)+p2mup4*p4mup3*p5mup6*F(1348)+p2mup4*p4mup6*p5mup
     -   3*F(1354)+p2mup4*p5mup3*p5mup6*F(1355)+p1mup4*(p1mup6*p2mup3
     -   *F(1292)+p1mup3*p2mup6*F(1292)+p2mup3*p2mup6*F(1292)+p1mup6*
     -   p3mup3*F(1294)+p1mup3*p3mup6*F(1295)+p3mup3*p3mup6*F(1296)+p
     -   3mup3*p5mup6*F(1298)+p4mup3*p4mup6*F(1312)+p3mup6*p5mup3*F(1
     -   344)+p4mup3*p5mup6*F(1348)+p4mup6*p5mup3*F(1354)+p5mup3*p5mu
     -   p6*F(1355))+p3mup4*p3mup6*p5mup3*F(1365))+p2mup2*(p1mup3*p2m
     -   up4*p2mup6*F(1292)+p2mup3*p2mup4*p2mup6*F(1292)-p1mup3*p2mup
     -   6*p3mup4*F(1293)-p2mup3*p2mup6*p3mup4*F(1293)+p1mup3*p2mup4*
     -   p3mup6*F(1295)+p2mup4*p3mup3*p3mup6*F(1296)+p2mup4*p3mup3*p5
     -   mup6*F(1298)+p1mup3*p3mup4*p3mup6*F(1302)+p3mup3*p3mup4*p3mu
     -   p6*F(1303)+p3mup3*p3mup4*p5mup6*F(1305)+p3mup3*p3mup4*p4mup6
     -   *F(1307)+p1mup3*p2mup6*p4mup4*F(1308)+p2mup3*p2mup6*p4mup4*F
     -   (1308)+p1mup6*(p2mup3*p2mup4*F(1292)-p2mup3*p3mup4*F(1293)+p
     -   2mup4*p3mup3*F(1294)-p3mup3*p3mup4*F(1301)+p2mup3*p4mup4*F(1
     -   308))+p1mup3*p3mup6*p4mup4*F(1309)-p3mup3*p4mup4*p5mup6*F(13
     -   11)+p2mup4*p4mup3*p4mup6*F(1312)+p3mup3*p3mup6*p4mup4*F(1314
     -   )+p3mup3*p4mup4*p4mup6*F(1315)+p3mup4*p4mup3*p4mup6*F(1316)+
     -   p3mup4*p4mup3*p5mup6*F(1317)+p3mup4*p4mup6*p5mup3*F(1320)+p3
     -   mup4*p5mup3*p5mup6*F(1323)+p3mup6*p4mup4*p5mup3*F(1324)+p4mu
     -   p3*p4mup4*p5mup6*F(1325)+p4mup4*p4mup6*p5mup3*F(1326)-p4mup4
     -   *p5mup3*p5mup6*F(1327)+p3mup6*p4mup3*p4mup4*F(1328)+p4mup3*p
     -   4mup4*p4mup6*F(1329)+p2mup4*p3mup6*p4mup3*F(1338)+p2mup4*p3m
     -   up6*p5mup3*F(1344)+p2mup4*p4mup3*p5mup6*F(1348)+p3mup4*p3mup
     -   6*p4mup3*F(1349)+p2mup4*p4mup6*p5mup3*F(1354)+p2mup4*p5mup3*
     -   p5mup6*F(1355)+p1mup4*(p1mup6*p2mup3*F(1292)+p1mup3*p2mup6*F
     -   (1292)+p2mup3*p2mup6*F(1292)+p1mup6*p3mup3*F(1294)+p1mup3*p3
     -   mup6*F(1295)+p3mup3*p3mup6*F(1296)+p3mup3*p5mup6*F(1298)+p4m
     -   up3*p4mup6*F(1312)+p3mup6*p5mup3*F(1344)+p4mup3*p5mup6*F(134
     -   8)+p4mup6*p5mup3*F(1354)+p5mup3*p5mup6*F(1355))+p3mup4*p3mup
     -   6*p5mup3*F(1365))
       Fa24(335) = p1mup2*(p1mup3*(p1mup4*p5mup6*F(1297)+p2mup4*p5mup6
     -   *F(1297)+p1mup4*p4mup6*F(1299)+p2mup4*p4mup6*F(1299)+p3mup4*
     -   p5mup6*F(1304)+p3mup4*p4mup6*F(1306)-p4mup4*p5mup6*F(1310)+p
     -   4mup4*p4mup6*F(1312))+p1mup6*(p1mup4*p4mup3*F(1308)+p2mup4*p
     -   4mup3*F(1308)+p3mup4*p4mup3*F(1339)+p1mup4*p5mup3*F(1341)+p2
     -   mup4*p5mup3*F(1341)+p3mup4*p5mup3*F(1352)+p4mup4*p5mup3*F(13
     -   57)+p4mup3*p4mup4*F(1359)))+p2mup2*(p1mup3*(p1mup4*p5mup6*F(
     -   1297)+p2mup4*p5mup6*F(1297)+p1mup4*p4mup6*F(1299)+p2mup4*p4m
     -   up6*F(1299)+p3mup4*p5mup6*F(1304)+p3mup4*p4mup6*F(1306)-p4mu
     -   p4*p5mup6*F(1310)+p4mup4*p4mup6*F(1312))+p1mup6*(p1mup4*p4mu
     -   p3*F(1308)+p2mup4*p4mup3*F(1308)+p3mup3*p4mup4*F(1313)+p3mup
     -   4*p4mup3*F(1339)+p1mup4*p5mup3*F(1341)+p2mup4*p5mup3*F(1341)
     -   +p3mup4*p5mup3*F(1352)+p4mup4*p5mup3*F(1357)+p4mup3*p4mup4*F
     -   (1359)))
       Fa24(336) = p5mup4*(16*p1mup2*(p3mup6*p4mup3*F(1330)+(p1mup3+p2
     -   mup3)*p3mup6*F(1346)+(p1mup3+p2mup3)*((p1mup6+p2mup6)*F(1345
     -   )+p4mup6*F(1356)-p5mup6*F(1358))+(p1mup6+p2mup6)*p4mup3*F(13
     -   60)-p3mup3*(p4mup6*F(1319)+(p1mup6+p2mup6)*F(1351)+p3mup6*F(
     -   1364)))+p1mup2*(p1mup6+p2mup6)*p5mup3*F(1407)+p2mup2*(16*(p1
     -   mup3*((p1mup6+p2mup6)*F(1345)+p3mup6*F(1346)+p4mup6*F(1356)-
     -   p5mup6*F(1358))+p1mup6*(p2mup3*F(1345)-p3mup3*F(1351)+p4mup3
     -   *F(1360)))+p1mup6*p5mup3*F(1407)))
       Return
       End
