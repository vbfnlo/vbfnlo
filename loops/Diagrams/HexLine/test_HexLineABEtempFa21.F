c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineAbe21(p1mup2,p1mup3,p1mup4,p1m
     -   up5,p2mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,
     -   p4mup2,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mu
     -   p2,p6mup3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,
     -   mup3mup5,mup4mup5,Fa21)
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
       Complex*16 Fa21(221:231)
       COMMON/test_FaFunctionsAbe/Fa
       Fa21(221) = 8*(p3mup4*p6mup5*F(1234)+p4mup4*p6mup5*F(1235)+p1mu
     -   p5*p3mup4*F(1236)+p2mup5*p3mup4*F(1236)+p1mup5*p4mup4*F(1237
     -   )+p2mup5*p4mup4*F(1237)+p1mup4*(p6mup5*F(945)+(p1mup5+p2mup5
     -   )*F(948)-p4mup5*F(1238)+p3mup5*F(1240))+p2mup4*(p6mup5*F(945
     -   )+(p1mup5+p2mup5)*F(948)-p4mup5*F(1238)+p3mup5*F(1240))+p3mu
     -   p4*p3mup5*F(1241)+p3mup5*p4mup4*F(1242)-p3mup4*p4mup5*F(1243
     -   )+p4mup4*p4mup5*F(1245)+p6mup4*(-(p6mup5*F(960))+(p1mup5+p2m
     -   up5)*F(966)-p4mup5*F(1247)+p3mup5*F(1248)))+4*mup4mup5*F(128
     -   3)
       Fa21(222) = mup4mup5*(p3mup3*F(1001)+p4mup3*F(1284)-(p1mup3+p2m
     -   up3)*F(1285))
       Fa21(223) = mup3mup4*p6mup5*F(1009)+mup3mup5*p6mup4*F(1132)+mup
     -   4mup5*p6mup3*F(1286)
       Fa21(224) = mup3mup4*(p3mup5*F(1004)+(p1mup5+p2mup5)*F(1006)-p4
     -   mup5*F(1008))
       Fa21(225) = p1mup3*(p1mup5+p2mup5)*(p4mup4*F(1012)+(p1mup4+p2mu
     -   p4)*F(1081)+p3mup4*F(1082))+p2mup3*(p1mup4*p2mup5*F(1081)+p1
     -   mup5*(p4mup4*F(1012)+(p1mup4+p2mup4)*F(1081)+p3mup4*F(1082))
     -   )+p1mup4*(p1mup5+p2mup5)*p3mup3*F(1087)
       Fa21(226) = -(p1mup4*p4mup3*p4mup5*F(1017))-p2mup4*p4mup3*p4mup
     -   5*F(1017)+p1mup5*p4mup3*p4mup4*F(1029)+p2mup5*p4mup3*p4mup4*
     -   F(1029)-p3mup4*p4mup3*p4mup5*F(1032)-p3mup3*p4mup4*p4mup5*F(
     -   1043)+p4mup3*p4mup4*p4mup5*F(1044)+p3mup5*p4mup3*p4mup4*F(10
     -   45)-p3mup3*p4mup4*p6mup5*F(1053)-p3mup4*p3mup5*p4mup3*F(1061
     -   )+p1mup5*p3mup3*p4mup4*F(1066)+p2mup5*p3mup3*p4mup4*F(1066)+
     -   p3mup3*p3mup4*p4mup5*F(1067)+p1mup4*p3mup3*p4mup5*F(1068)+p2
     -   mup4*p3mup3*p4mup5*F(1068)-p3mup3*p3mup5*p4mup4*F(1076)+p1mu
     -   p3*(-(p4mup4*p6mup5*F(1015))+p4mup4*p4mup5*F(1016)+p1mup4*p4
     -   mup5*F(1020)+p2mup4*p4mup5*F(1020)+p3mup4*p4mup5*F(1021)-p1m
     -   up4*p6mup5*F(1083)-p2mup4*p6mup5*F(1083)+p3mup5*(p4mup4*F(10
     -   14)+(p1mup4+p2mup4)*F(1082)+p3mup4*F(1084))-p3mup4*p6mup5*F(
     -   1085))+p2mup3*(-(p4mup4*p6mup5*F(1015))+p4mup4*p4mup5*F(1016
     -   )+p1mup4*p4mup5*F(1020)+p2mup4*p4mup5*F(1020)+p3mup4*p4mup5*
     -   F(1021)+p2mup5*(p4mup4*F(1012)+p2mup4*F(1081)+p3mup4*F(1082)
     -   )-p1mup4*p6mup5*F(1083)-p2mup4*p6mup5*F(1083)+p3mup5*(p4mup4
     -   *F(1014)+(p1mup4+p2mup4)*F(1082)+p3mup4*F(1084))-p3mup4*p6mu
     -   p5*F(1085))-p1mup4*p3mup3*p6mup5*F(1086)-p2mup4*p3mup3*p6mup
     -   5*F(1086)+p1mup5*p2mup4*p3mup3*F(1087)+p2mup4*p2mup5*p3mup3*
     -   F(1087)+p1mup5*p3mup3*p3mup4*F(1088)+p2mup5*p3mup3*p3mup4*F(
     -   1088)+p1mup4*p3mup3*p3mup5*F(1088)+p2mup4*p3mup3*p3mup5*F(10
     -   88)-p3mup3*p3mup4*p6mup5*F(1095)+p3mup3*p3mup4*p3mup5*F(1100
     -   )-p3mup5*(p3mup4*p6mup3+p3mup3*p6mup4)*F(1103)
       Fa21(227) = 16*(p2mup4*p2mup5*p4mup3*F(1023)-p4mup3*p4mup4*p6mu
     -   p5*F(1033)+p4mup3*p4mup5*p6mup4*F(1034)-p4mup4*p4mup5*p6mup3
     -   *F(1035)+p1mup3*p4mup5*p6mup4*F(1040)+p2mup3*p4mup5*p6mup4*F
     -   (1040)+p2mup4*p3mup5*p4mup3*F(1049)-p2mup5*p3mup4*p4mup3*F(1
     -   051)-p2mup4*p4mup5*p6mup3*F(1056)-p2mup5*p4mup3*p6mup4*F(105
     -   8)-p3mup4*p4mup5*p6mup3*F(1062)+p2mup4*p4mup3*p6mup5*F(1069)
     -   +p2mup5*p4mup4*p6mup3*F(1070)+p3mup4*p4mup3*p6mup5*F(1071)-p
     -   3mup5*p4mup3*p6mup4*F(1072)+p3mup5*p4mup4*p6mup3*F(1073)+p4m
     -   up4*p6mup3*p6mup5*F(1074)+p3mup3*p4mup5*p6mup4*F(1077)+p4mup
     -   3*p6mup4*p6mup5*F(1078)-p4mup5*p6mup3*p6mup4*F(1079)-p1mup3*
     -   p2mup5*p6mup4*F(1089)-p2mup3*p2mup5*p6mup4*F(1089)-p1mup3*p3
     -   mup5*p6mup4*F(1090)-p2mup3*p3mup5*p6mup4*F(1090)+p1mup3*p6mu
     -   p4*p6mup5*F(1091)+p2mup3*p6mup4*p6mup5*F(1091)+p2mup4*p6mup3
     -   *p6mup5*F(1092)-p2mup4*p2mup5*p6mup3*F(1093)+p2mup4*p3mup5*p
     -   6mup3*F(1094)+p2mup5*p3mup3*p6mup4*F(1094)+p1mup4*(p1mup5*p4
     -   mup3*F(1023)+p2mup5*p4mup3*F(1023)+p3mup5*p4mup3*F(1049)-p4m
     -   up5*p6mup3*F(1056)+p4mup3*p6mup5*F(1069)+p6mup3*p6mup5*F(109
     -   2)-p1mup5*p6mup3*F(1093)-p2mup5*p6mup3*F(1093)+p3mup5*p6mup3
     -   *F(1094))+p3mup4*p6mup3*p6mup5*F(1101)+p3mup3*p6mup4*p6mup5*
     -   F(1101)-p2mup5*p3mup4*p6mup3*F(1102)+p2mup5*p6mup3*p6mup4*F(
     -   1107)+p1mup5*(p4mup4*p6mup3*F(1070)+p2mup4*(p4mup3*F(1023)-p
     -   6mup3*F(1093))-p3mup4*(p4mup3*F(1051)+p6mup3*F(1102))-p6mup4
     -   *(p4mup3*F(1058)+(p1mup3+p2mup3)*F(1089)-p3mup3*F(1094)-p6mu
     -   p3*F(1107)))+p6mup3*p6mup4*p6mup5*F(1108)+p3mup5*p6mup3*p6mu
     -   p4*F(1136))+8*(mup3mup5*p4mup4*F(1129)+mup3mup5*p3mup4*F(113
     -   7))+mup3mup5*(p1mup4+p2mup4)*F(1287)
       Fa21(228) = -8*(p3mup4*((p1mup5+p2mup5)*F(950)+p4mup5*F(991)+p6
     -   mup5*F(1146))-p3mup5*((p1mup4+p2mup4)*F(913)+p4mup4*F(915)+p
     -   6mup4*F(1139)-p3mup4*F(1246)))+4*mup4mup5*F(1288)
       Fa21(229) = p3mup2*(2*(-((p1mup5+p2mup5)*p3mup4*F(1268))-p4mup4
     -   *p6mup5*F(1270)+p1mup4*(p3mup5*F(1267)+p4mup5*F(1269)+(p1mup
     -   5+p2mup5)*F(1271)-p6mup5*F(1272))+p2mup4*(p3mup5*F(1267)+p4m
     -   up5*F(1269)+(p1mup5+p2mup5)*F(1271)-p6mup5*F(1272))+p4mup5*p
     -   6mup4*F(1274)-(p1mup5+p2mup5)*p6mup4*F(1275)+p3mup4*(-(p6mup
     -   5*F(1273))+p3mup5*F(1277)+p4mup5*F(1278))+p4mup4*(p3mup5*F(1
     -   279)+p4mup5*F(1280)+(p1mup5+p2mup5)*F(1281))+p6mup4*(p6mup5*
     -   F(1276)-p3mup5*F(1282)))+mup4mup5*F(1289))
       Fa21(230) = mup4mup5*(-8*p4mup3*F(999)+(p1mup3+p2mup3)*F(1290)+
     -   p3mup3*F(1291))
       Fa21(231) = 8*(mup4mup5*p6mup3*F(1000)+mup3mup4*p6mup5*F(1010))
     -   +mup3mup5*p6mup4*F(1292)
       Return
       End
