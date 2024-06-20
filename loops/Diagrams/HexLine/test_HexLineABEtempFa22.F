c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineAbe22(p1mup2,p1mup3,p1mup4,p1m
     -   up5,p2mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,
     -   p4mup2,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mu
     -   p2,p6mup3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,
     -   mup3mup5,mup4mup5,Fa22)
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
       Complex*16 Fa22(232:242)
       COMMON/test_FaFunctionsAbe/Fa
       Fa22(232) = mup3mup4*((p1mup5+p2mup5)*F(1293)+p3mup5*F(1294)+p4
     -   mup5*F(1295))
       Fa22(233) = -(p1mup5*p2mup3*p3mup4*F(1013))+p1mup4*p1mup5*p3mup
     -   3*F(1024)+p1mup4*p2mup5*p3mup3*F(1024)-p1mup5*p2mup3*p4mup4*
     -   F(1099)+p2mup3*(p1mup5*p2mup4+p1mup4*(p1mup5+p2mup5))*F(1112
     -   )+p1mup3*(p1mup5+p2mup5)*(-(p3mup4*F(1013))-p4mup4*F(1099)+(
     -   p1mup4+p2mup4)*F(1112))
       Fa22(234) = p1mup3*p1mup4*p3mup5*F(1011)+p1mup4*p2mup3*p3mup5*F
     -   (1011)+p1mup3*p2mup4*p3mup5*F(1011)+p2mup3*p2mup4*p3mup5*F(1
     -   011)-p2mup3*p2mup5*p3mup4*F(1013)+p1mup4*p3mup3*p4mup5*F(101
     -   8)+p2mup4*p3mup3*p4mup5*F(1018)-p3mup3*p4mup4*p6mup5*F(1019)
     -   -p1mup3*p3mup4*p6mup5*F(1022)-p2mup3*p3mup4*p6mup5*F(1022)+p
     -   1mup5*p2mup4*p3mup3*F(1024)+p2mup4*p2mup5*p3mup3*F(1024)-p3m
     -   up4*p4mup3*p4mup5*F(1025)-p1mup3*p3mup4*p4mup5*F(1026)-p2mup
     -   3*p3mup4*p4mup5*F(1026)-p1mup3*p3mup5*p4mup4*F(1027)-p2mup3*
     -   p3mup5*p4mup4*F(1027)-p1mup5*p3mup3*p4mup4*F(1030)-p2mup5*p3
     -   mup3*p4mup4*F(1030)-p3mup3*p3mup4*p4mup5*F(1032)-p3mup4*p3mu
     -   p5*p4mup3*F(1037)-p1mup3*p3mup4*p3mup5*F(1041)-p2mup3*p3mup4
     -   *p3mup5*F(1041)+p3mup3*p4mup4*p4mup5*F(1044)+p3mup3*p3mup5*p
     -   4mup4*F(1045)-p3mup5*p4mup3*p4mup4*F(1046)+p1mup4*p3mup3*p3m
     -   up5*F(1050)+p2mup4*p3mup3*p3mup5*F(1050)-p1mup5*p3mup3*p3mup
     -   4*F(1052)-p2mup5*p3mup3*p3mup4*F(1052)-p1mup4*p3mup3*p6mup5*
     -   F(1054)-p2mup4*p3mup3*p6mup5*F(1054)-p3mup3*p3mup4*p6mup5*F(
     -   1055)-p3mup3*p3mup4*p3mup5*F(1061)-p3mup4*p3mup5*p6mup3*F(10
     -   64)-p3mup3*p3mup5*p6mup4*F(1075)+p1mup3*p1mup4*p4mup5*F(1096
     -   )+p1mup4*p2mup3*p4mup5*F(1096)+p1mup3*p2mup4*p4mup5*F(1096)+
     -   p2mup3*p2mup4*p4mup5*F(1096)-p1mup4*p4mup3*p4mup5*F(1097)-p2
     -   mup4*p4mup3*p4mup5*F(1097)-p1mup3*p4mup4*p4mup5*F(1097)-p2mu
     -   p3*p4mup4*p4mup5*F(1097)-p2mup3*p2mup5*p4mup4*F(1099)-p4mup3
     -   *p4mup4*p4mup5*F(1104)-p1mup5*p4mup3*p4mup4*F(1105)-p2mup5*p
     -   4mup3*p4mup4*F(1105)-p1mup3*p4mup4*p6mup5*F(1111)-p2mup3*p4m
     -   up4*p6mup5*F(1111)+p2mup3*p2mup4*p2mup5*F(1112)-(p1mup3+p2mu
     -   p3)*(p1mup4+p2mup4)*p6mup5*F(1113)
       Fa22(235) = 16*(p1mup5*p3mup4*p4mup3*F(1031)+p2mup5*p3mup4*p4mu
     -   p3*F(1031)-p3mup4*p4mup3*p6mup5*F(1038)+p3mup4*p4mup5*p6mup3
     -   *F(1039)+p6mup4*(-(p3mup3*p4mup5*F(1036))+(p1mup3+p2mup3)*p3
     -   mup5*F(1042))+p3mup5*p4mup3*p6mup4*F(1047)+p3mup5*p4mup4*p6m
     -   up3*F(1048)+p1mup5*p3mup3*p6mup4*F(1059)+p2mup5*p3mup3*p6mup
     -   4*F(1059)+p1mup5*p3mup4*p6mup3*F(1060)+p2mup5*p3mup4*p6mup3*
     -   F(1060)-p3mup3*p6mup4*p6mup5*F(1063)-p3mup4*p6mup3*p6mup5*F(
     -   1065)+p3mup5*p6mup3*p6mup4*F(1080)-p4mup3*p4mup4*p6mup5*F(10
     -   98)+p4mup4*p4mup5*p6mup3*F(1106)+p4mup3*p4mup5*p6mup4*F(1109
     -   )-p4mup4*p6mup3*p6mup5*F(1116)+p1mup3*p1mup5*p6mup4*F(1117)+
     -   p1mup5*p2mup3*p6mup4*F(1117)+p1mup3*p2mup5*p6mup4*F(1117)+p2
     -   mup3*p2mup5*p6mup4*F(1117)-p1mup3*p6mup4*p6mup5*F(1118)-p2mu
     -   p3*p6mup4*p6mup5*F(1118)+p1mup5*p4mup4*p6mup3*F(1120)+p2mup5
     -   *p4mup4*p6mup3*F(1120)+p1mup3*p4mup5*p6mup4*F(1120)+p2mup3*p
     -   4mup5*p6mup4*F(1120)+p1mup5*p4mup3*p6mup4*F(1121)+p2mup5*p4m
     -   up3*p6mup4*F(1121)-p4mup3*p6mup4*p6mup5*F(1122)+p1mup4*(p3mu
     -   p5*p4mup3*F(1028)+p3mup5*p6mup3*F(1057)+p1mup5*p4mup3*F(1099
     -   )+p2mup5*p4mup3*F(1099)-p4mup3*p6mup5*F(1114)+p4mup5*p6mup3*
     -   F(1115)+p1mup5*p6mup3*F(1119)+p2mup5*p6mup3*F(1119)-p6mup3*p
     -   6mup5*F(1124))+p2mup4*(p3mup5*p4mup3*F(1028)+p3mup5*p6mup3*F
     -   (1057)+p1mup5*p4mup3*F(1099)+p2mup5*p4mup3*F(1099)-p4mup3*p6
     -   mup5*F(1114)+p4mup5*p6mup3*F(1115)+p1mup5*p6mup3*F(1119)+p2m
     -   up5*p6mup3*F(1119)-p6mup3*p6mup5*F(1124))+p1mup5*p6mup3*p6mu
     -   p4*F(1125)+p2mup5*p6mup3*p6mup4*F(1125)+p6mup3*p6mup4*p6mup5
     -   *F(1126)+p4mup5*p6mup3*p6mup4*F(1147))+mup3mup5*p1mup4*F(129
     -   6)+mup3mup5*p2mup4*F(1296)+mup3mup5*p3mup4*F(1297)+mup3mup5*
     -   p4mup4*F(1298)
       Fa22(236) = (p1mup2+p2mup2)*(mup4mup5*F(1266)-2*(-((p1mup5+p2mu
     -   p5)*p3mup4*F(1268))-p4mup4*p6mup5*F(1270)+p1mup4*(p3mup5*F(1
     -   267)+p4mup5*F(1269)+(p1mup5+p2mup5)*F(1271)-p6mup5*F(1272))+
     -   p2mup4*(p3mup5*F(1267)+p4mup5*F(1269)+(p1mup5+p2mup5)*F(1271
     -   )-p6mup5*F(1272))-p3mup4*p6mup5*F(1273)+p4mup5*p6mup4*F(1274
     -   )-(p1mup5+p2mup5)*p6mup4*F(1275)+p6mup4*p6mup5*F(1276)+p3mup
     -   4*p3mup5*F(1277)+p3mup4*p4mup5*F(1278)+p3mup5*p4mup4*F(1279)
     -   +p4mup4*p4mup5*F(1280)+p1mup5*p4mup4*F(1281)+p2mup5*p4mup4*F
     -   (1281)-p3mup5*p6mup4*F(1282)))
       Fa22(237) = 8*(-((p1mup5+p2mup5)*p3mup4*F(1268))-p4mup4*p6mup5*
     -   F(1270)+p1mup4*(p3mup5*F(1267)+p4mup5*F(1269)+(p1mup5+p2mup5
     -   )*F(1271)-p6mup5*F(1272))+p2mup4*(p3mup5*F(1267)+p4mup5*F(12
     -   69)+(p1mup5+p2mup5)*F(1271)-p6mup5*F(1272))+p4mup5*p6mup4*F(
     -   1274)-(p1mup5+p2mup5)*p6mup4*F(1275)+p3mup4*(-(p6mup5*F(1273
     -   ))+p3mup5*F(1277)+p4mup5*F(1278))+p4mup4*(p3mup5*F(1279)+p4m
     -   up5*F(1280)+(p1mup5+p2mup5)*F(1281))+p6mup4*(p6mup5*F(1276)-
     -   p3mup5*F(1282)))+4*mup4mup5*F(1289)
       Fa22(238) = -2*(p4mup3*(p6mup5*F(946)+2*((p1mup5+p2mup5)*F(950)
     -   +p3mup5*F(954)))-p4mup5*((p1mup3+p2mup3)*F(814)-p3mup3*F(815
     -   )-p4mup3*F(958)-p6mup3*F(968)))-mup3mup5*F(1299)
       Fa22(239) = 8*(p3mup3*p6mup5*F(947)+p1mup5*p3mup3*F(951)+p2mup5
     -   *p3mup3*F(951)+p3mup3*p3mup5*F(955)+p4mup3*p6mup5*F(956)+p1m
     -   up3*(p6mup5*F(945)+(p1mup5+p2mup5)*F(948)-p3mup5*F(949)+p4mu
     -   p5*F(957))+p2mup3*(p6mup5*F(945)+(p1mup5+p2mup5)*F(948)-p3mu
     -   p5*F(949)+p4mup5*F(957))+p3mup3*p4mup5*F(959)-p6mup3*p6mup5*
     -   F(960)-p4mup5*p6mup3*F(962)-p1mup5*p4mup3*F(963)-p2mup5*p4mu
     -   p3*F(963)-p4mup3*p4mup5*F(964)-p3mup5*p4mup3*F(965)+p1mup5*p
     -   6mup3*F(966)+p2mup5*p6mup3*F(966)+p3mup5*p6mup3*F(967))+4*mu
     -   p3mup5*F(1300)
       Fa22(240) = -24*(p1mup5*p3mup2*s12*F(950)+p2mup5*p3mup2*s12*F(9
     -   50)+p3mup2*p3mup5*s12*F(954)+p3mup2*p4mup5*s12*F(991))+p1mup
     -   5*(p1mup2+p2mup2)*F(1301)+(p1mup2+p2mup2)*p2mup5*F(1301)+8*(
     -   -(p3mup2*p6mup5*s12*F(990))+p1mup2*p3mup5*F(1302)+p2mup2*p3m
     -   up5*F(1302))-p1mup2*p6mup5*F(1303)-p2mup2*p6mup5*F(1303)+(p1
     -   mup2+p2mup2)*p4mup5*F(1304)
       Fa22(241) = 8*(p2mup3*p6mup4*F(841)-p3mup3*p6mup4*F(842)-p2mup3
     -   *p3mup4*F(850)+p1mup4*p3mup3*F(852)+p2mup4*p3mup3*F(852)+p3m
     -   up3*p3mup4*F(853)-p6mup3*(-(p6mup4*F(840))+p3mup4*F(856)+(p1
     -   mup4+p2mup4)*F(857)+p4mup4*F(858))-p1mup4*p2mup3*F(859)-p2mu
     -   p3*p2mup4*F(859)-p2mup3*p4mup4*F(861)-p1mup3*(-(p6mup4*F(841
     -   ))+p3mup4*F(850)+(p1mup4+p2mup4)*F(859)+p4mup4*F(861))+p4mup
     -   3*(p6mup4*F(836)+(p1mup4+p2mup4)*F(863)+p3mup4*F(865)+p4mup4
     -   *F(868))+p3mup3*p4mup4*F(873))+4*mup3mup4*F(1305)
       Fa22(242) = 2*(p3mup3*p6mup4*F(837)+p4mup3*p6mup4*F(839)-p6mup3
     -   *p6mup4*F(851)+p1mup4*p3mup3*F(864)+p2mup4*p3mup3*F(864)+p3m
     -   up3*p3mup4*F(866)+p1mup4*p4mup3*F(867)+p2mup4*p4mup3*F(867)+
     -   p1mup3*(p6mup4*F(835)+(p1mup4+p2mup4)*F(860)-p3mup4*F(862)+p
     -   4mup4*F(867))+p2mup3*(p6mup4*F(835)+(p1mup4+p2mup4)*F(860)-p
     -   3mup4*F(862)+p4mup4*F(867))+p3mup3*p4mup4*F(869)+p3mup4*p4mu
     -   p3*F(870)+p1mup4*p6mup3*F(871)+p2mup4*p6mup3*F(871)+p3mup4*p
     -   6mup3*F(872)+p4mup3*p4mup4*F(874)+p4mup4*p6mup3*F(875))-mup3
     -   mup4*F(1306)
       Return
       End
