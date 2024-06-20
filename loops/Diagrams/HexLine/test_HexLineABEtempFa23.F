c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineAbe23(p1mup2,p1mup3,p1mup4,p1m
     -   up5,p2mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,
     -   p4mup2,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mu
     -   p2,p6mup3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,
     -   mup3mup5,mup4mup5,Fa23)
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
       Complex*16 Fa23(243:253)
       COMMON/test_FaFunctionsAbe/Fa
       Fa23(243) = 8*(-(p1mup2*(p3mup4*F(887)+p4mup4*F(891)))-p2mup2*(
     -   p3mup4*F(887)+p4mup4*F(891))+p3mup2*s12*(p6mup4*F(912)+3*(p1
     -   mup4+p2mup4)*F(913))+3*p3mup2*s12*(p3mup4*F(914)+p4mup4*F(91
     -   5))+(p1mup2+p2mup2)*(p1mup4+p2mup4)*F(1307))+(p1mup2+p2mup2)
     -   *p6mup4*F(1308)
       Fa23(244) = -4*(p3mup4*F(887)+p4mup4*F(891)-(p1mup4+p2mup4)*F(1
     -   307))+p6mup4*F(1309)
       Fa23(245) = 4*p3mup5*F(1302)+(p1mup5+p2mup5)*F(1310)-p6mup5*F(1
     -   311)+p4mup5*F(1312)
       Fa23(246) = p4mup2*(mup3mup5*F(1313)+2*((p1mup3+p2mup3)*(-(p6mu
     -   p5*F(1314))+(p1mup5+p2mup5)*F(1316))+p1mup3*p4mup5*F(1319)+p
     -   2mup3*p4mup5*F(1319)+p4mup5*p6mup3*F(1320)+p4mup3*(p4mup5*F(
     -   1278)-p6mup5*F(1318)+(p1mup5+p2mup5)*F(1321))+p3mup3*(p4mup5
     -   *F(1277)-p6mup5*F(1315)+(p1mup5+p2mup5)*F(1322))-p6mup3*p6mu
     -   p5*F(1324)-(p1mup5+p2mup5)*p6mup3*F(1325)+p3mup5*(p4mup3*F(1
     -   277)+(p1mup3+p2mup3)*F(1317)-p3mup3*F(1323)-p6mup3*F(1326)))
     -   )
       Fa23(247) = p3mup2*(2*(p4mup3*p4mup5*F(1278)+(p1mup3+p2mup3)*(-
     -   (p6mup5*F(1314))+(p1mup5+p2mup5)*F(1316))-p4mup3*p6mup5*F(13
     -   18)+p1mup3*p4mup5*F(1319)+p2mup3*p4mup5*F(1319)+p4mup5*p6mup
     -   3*F(1320)+p1mup5*p4mup3*F(1321)+p2mup5*p4mup3*F(1321)+p3mup3
     -   *(p4mup5*F(1277)-p6mup5*F(1315)+(p1mup5+p2mup5)*F(1322))-p6m
     -   up3*p6mup5*F(1324)-(p1mup5+p2mup5)*p6mup3*F(1325)+p3mup5*(p4
     -   mup3*F(1277)+(p1mup3+p2mup3)*F(1317)-p3mup3*F(1323)-p6mup3*F
     -   (1326)))-mup3mup5*F(1327))
       Fa23(248) = (p1mup2+p2mup2)*(mup3mup5*F(1313)+2*((p1mup3+p2mup3
     -   )*(-(p6mup5*F(1314))+(p1mup5+p2mup5)*F(1316))+p1mup3*p4mup5*
     -   F(1319)+p2mup3*p4mup5*F(1319)+p4mup5*p6mup3*F(1320)+p4mup3*(
     -   p4mup5*F(1278)-p6mup5*F(1318)+(p1mup5+p2mup5)*F(1321))+p3mup
     -   3*(p4mup5*F(1277)-p6mup5*F(1315)+(p1mup5+p2mup5)*F(1322))-p6
     -   mup3*p6mup5*F(1324)-(p1mup5+p2mup5)*p6mup3*F(1325)+p3mup5*(p
     -   4mup3*F(1277)+(p1mup3+p2mup3)*F(1317)-p3mup3*F(1323)-p6mup3*
     -   F(1326))))
       Fa23(249) = p6mup5*F(990)+3*((p1mup5+p2mup5)*F(950)+p3mup5*F(95
     -   4)+p4mup5*F(991))
       Fa23(250) = 8*(p4mup3*p4mup5*F(1278)+(p1mup3+p2mup3)*(-(p6mup5*
     -   F(1314))+(p1mup5+p2mup5)*F(1316))-p4mup3*p6mup5*F(1318)+p1mu
     -   p3*p4mup5*F(1319)+p2mup3*p4mup5*F(1319)+p4mup5*p6mup3*F(1320
     -   )+p1mup5*p4mup3*F(1321)+p2mup5*p4mup3*F(1321)+p3mup3*(p4mup5
     -   *F(1277)-p6mup5*F(1315)+(p1mup5+p2mup5)*F(1322))-p6mup3*p6mu
     -   p5*F(1324)-(p1mup5+p2mup5)*p6mup3*F(1325)+p3mup5*(p4mup3*F(1
     -   277)+(p1mup3+p2mup3)*F(1317)-p3mup3*F(1323)-p6mup3*F(1326)))
     -   -4*mup3mup5*F(1327)
       Fa23(251) = mup4mup5*((p1mup3+p2mup3)*F(814)-p3mup3*F(815)+p4mu
     -   p3*F(816)-p6mup3*F(968))-mup3mup5*((p1mup4+p2mup4)*F(913)+p3
     -   mup4*F(914)+p4mup4*F(915)+p6mup4*F(1139))+mup3mup4*((p1mup5+
     -   p2mup5)*F(950)+p3mup5*F(954)+p4mup5*F(991)+p6mup5*F(1146))
       Fa23(252) = p4mup2*(mup3mup4*F(1328)-2*(-(p3mup3*p6mup4*F(1330)
     -   )-p4mup3*p6mup4*F(1333)+p4mup4*p6mup3*F(1334)+p1mup4*p6mup3*
     -   F(1335)+p2mup4*p6mup3*F(1335)+p3mup4*p6mup3*F(1336)-p6mup3*p
     -   6mup4*F(1337)-p3mup3*p3mup4*F(1338)+p1mup3*(-(p6mup4*F(1329)
     -   )+(p1mup4+p2mup4)*F(1331)+p3mup4*F(1332)+p4mup4*F(1339))+p2m
     -   up3*(-(p6mup4*F(1329))+(p1mup4+p2mup4)*F(1331)+p3mup4*F(1332
     -   )+p4mup4*F(1339))-p1mup4*p4mup3*F(1340)-p2mup4*p4mup3*F(1340
     -   )-p3mup4*p4mup3*F(1341)-p1mup4*p3mup3*F(1342)-p2mup4*p3mup3*
     -   F(1342)+p3mup3*p4mup4*F(1343)-p4mup3*p4mup4*F(1344)))
       Fa23(253) = p3mup2*(2*(p4mup3*p6mup4*F(1333)-p6mup3*(p4mup4*F(1
     -   334)+(p1mup4+p2mup4)*F(1335)+p3mup4*F(1336))+p6mup3*p6mup4*F
     -   (1337)+p1mup3*(p6mup4*F(1329)-(p1mup4+p2mup4)*F(1331)-p3mup4
     -   *F(1332)-p4mup4*F(1339))+p2mup3*(p6mup4*F(1329)-(p1mup4+p2mu
     -   p4)*F(1331)-p3mup4*F(1332)-p4mup4*F(1339))+p1mup4*p4mup3*F(1
     -   340)+p2mup4*p4mup3*F(1340)+p3mup4*p4mup3*F(1341)+p3mup3*(p6m
     -   up4*F(1330)+p3mup4*F(1338)+(p1mup4+p2mup4)*F(1342)-p4mup4*F(
     -   1343))+p4mup3*p4mup4*F(1344))-mup3mup4*F(1345))
       Return
       End
