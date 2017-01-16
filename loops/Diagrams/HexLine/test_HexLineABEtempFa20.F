c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineAbe20(p1mup2,p1mup3,p1mup4,p1m
     -   up5,p2mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,
     -   p4mup2,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mu
     -   p2,p6mup3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,
     -   mup3mup5,mup4mup5,Fa20)
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
       Complex*16 Fa20(210:220)
       COMMON/test_FaFunctionsAbe/Fa
       Fa20(210) = p1mup2*(-(p1mup4*p2mup3*p2mup5*F(1112))-p2mup3*p2mu
     -   p4*p2mup5*F(1112)+p1mup4*p2mup3*p6mup5*F(1113)+p2mup3*p2mup4
     -   *p6mup5*F(1113)+p2mup3*p2mup5*p6mup4*F(1117)-p2mup3*p6mup4*p
     -   6mup5*F(1118)+p1mup3*(p1mup4*(-(p1mup5*F(1112))+p6mup5*F(111
     -   3))+p2mup4*(-(p2mup5*F(1112))+p6mup5*F(1113))+(p1mup5+p2mup5
     -   )*p6mup4*F(1117)-p6mup4*p6mup5*F(1118))+p1mup4*p2mup5*p6mup3
     -   *F(1119)+p2mup4*p2mup5*p6mup3*F(1119)-p1mup4*p6mup3*p6mup5*F
     -   (1124)-p2mup4*p6mup3*p6mup5*F(1124)+p2mup5*p6mup3*p6mup4*F(1
     -   125)+p1mup5*(p2mup3*(-(p2mup4*F(1112))+p6mup4*F(1117))+(p1mu
     -   p4+p2mup4)*p6mup3*F(1119)+p6mup3*p6mup4*F(1125))+p6mup3*p6mu
     -   p4*p6mup5*F(1126))+p2mup2*(p2mup4*p2mup5*p4mup3*F(1099)-p1mu
     -   p3*p1mup4*p2mup5*F(1112)-p1mup3*p2mup4*p2mup5*F(1112)+p1mup3
     -   *p1mup4*p6mup5*F(1113)+p1mup3*p2mup4*p6mup5*F(1113)-p2mup4*p
     -   4mup3*p6mup5*F(1114)+p2mup4*p4mup5*p6mup3*F(1115)-p4mup4*p6m
     -   up3*p6mup5*F(1116)+p1mup3*p2mup5*p6mup4*F(1117)-p1mup3*p6mup
     -   4*p6mup5*F(1118)+p1mup4*p2mup5*p6mup3*F(1119)+p2mup5*p4mup4*
     -   p6mup3*F(1120)+p2mup3*(-(p2mup4*p4mup5*F(1096))+p1mup5*p4mup
     -   4*F(1099)+p2mup5*p4mup4*F(1099)+p4mup4*p6mup5*F(1111)-p1mup4
     -   *p1mup5*F(1112)-p1mup4*p2mup5*F(1112)+p1mup4*p6mup5*F(1113)+
     -   p4mup5*p6mup4*F(1120))+p2mup5*p4mup3*p6mup4*F(1121)+p1mup5*(
     -   p2mup4*p4mup3*F(1099)-p1mup3*p2mup4*F(1112)+p1mup3*p6mup4*F(
     -   1117)+p1mup4*p6mup3*F(1119)+p4mup4*p6mup3*F(1120)+p4mup3*p6m
     -   up4*F(1121))-p4mup3*p6mup4*p6mup5*F(1122)-p1mup4*p6mup3*p6mu
     -   p5*F(1124)+p4mup5*p6mup3*p6mup4*F(1147))
       Fa20(211) = 32*p2mup2*(p2mup3*(-(p2mup4*p2mup5*F(1112))+p2mup4*
     -   p6mup5*F(1113)+p2mup5*p6mup4*F(1117)-p6mup4*p6mup5*F(1118))+
     -   p1mup5*(-(p2mup3*p2mup4*F(1112))+p2mup3*p6mup4*F(1117)+p2mup
     -   4*p6mup3*F(1119)+p6mup3*p6mup4*F(1125))+p6mup3*(p2mup4*p2mup
     -   5*F(1119)-p2mup4*p6mup5*F(1124)+p2mup5*p6mup4*F(1125)+p6mup4
     -   *p6mup5*F(1126)))+mup3mup5*(p1mup4*(p1mup2+p2mup2)*F(1148)+(
     -   p1mup2+p2mup2)*p2mup4*F(1148)+p1mup2*p3mup4*F(1149)+p2mup2*p
     -   3mup4*F(1149)+p1mup2*p4mup4*F(1150)+p2mup2*p4mup4*F(1150)+p1
     -   mup2*p6mup4*F(1151)+p2mup2*p6mup4*F(1151)-16*(p1mup4*p3mup2*
     -   s12*F(913)+p2mup4*p3mup2*s12*F(913)+p3mup2*p3mup4*s12*F(914)
     -   +p3mup2*p4mup4*s12*F(915)-p3mup2*p6mup4*s12*F(1152)))
       Fa20(212) = mup4mup5*(p4mup3*F(1153)+(p1mup3+p2mup3)*F(1154)+p3
     -   mup3*F(1155))
       Fa20(213) = mup4mup5*p6mup3*F(1156)+mup3mup4*p6mup5*F(1157)+mup
     -   3mup5*p6mup4*F(1158)
       Fa20(214) = mup3mup4*((p1mup5+p2mup5)*F(1159)+p3mup5*F(1160)+p4
     -   mup5*F(1161))
       Fa20(215) = p1mup5*p2mup3*(p2mup4*F(1162)+p3mup4*F(1163)-p4mup4
     -   *F(1164))+p1mup3*(p1mup5+p2mup5)*((p1mup4+p2mup4)*F(1162)+p3
     -   mup4*F(1163)-p4mup4*F(1164))+p1mup4*(p1mup5+p2mup5)*(p2mup3*
     -   F(1162)+p3mup3*F(1165))
       Fa20(216) = p1mup5*p2mup4*p3mup3*F(1165)+p2mup4*p2mup5*p3mup3*F
     -   (1165)-p1mup4*p3mup3*p6mup5*F(1169)-p2mup4*p3mup3*p6mup5*F(1
     -   169)-p3mup3*p3mup4*p6mup5*F(1170)-p3mup3*p4mup4*p6mup5*F(117
     -   1)-p1mup5*p3mup3*p4mup4*F(1172)-p2mup5*p3mup3*p4mup4*F(1172)
     -   +p1mup4*p3mup3*p4mup5*F(1176)+p2mup4*p3mup3*p4mup5*F(1176)-p
     -   3mup3*p3mup5*p4mup4*F(1177)+p3mup3*p3mup4*p4mup5*F(1179)-p3m
     -   up3*p4mup4*p4mup5*F(1180)+p3mup4*p3mup5*p4mup3*F(1181)+p3mup
     -   4*p4mup3*p4mup5*F(1182)-p1mup4*p4mup3*p4mup5*F(1183)-p2mup4*
     -   p4mup3*p4mup5*F(1183)+p1mup3*(p1mup4*p3mup5*F(1163)+p2mup4*p
     -   3mup5*F(1163)-p1mup4*p6mup5*F(1166)-p2mup4*p6mup5*F(1166)-p3
     -   mup4*p6mup5*F(1167)-p4mup4*p6mup5*F(1168)+p3mup4*p3mup5*F(11
     -   73)+p3mup5*p4mup4*F(1174)+p1mup4*p4mup5*F(1175)+p2mup4*p4mup
     -   5*F(1175)+p3mup4*p4mup5*F(1178)-p4mup4*p4mup5*F(1183))+p2mup
     -   3*(p2mup5*p3mup4*F(1163)+p1mup4*p3mup5*F(1163)-p2mup5*p4mup4
     -   *F(1164)-p1mup4*p6mup5*F(1166)-p3mup4*p6mup5*F(1167)-p4mup4*
     -   p6mup5*F(1168)+p3mup4*p3mup5*F(1173)+p3mup5*p4mup4*F(1174)+p
     -   1mup4*p4mup5*F(1175)+p2mup4*(p2mup5*F(1162)+p3mup5*F(1163)-p
     -   6mup5*F(1166)+p4mup5*F(1175))+p3mup4*p4mup5*F(1178)-p4mup4*p
     -   4mup5*F(1183))-p3mup5*p4mup3*p4mup4*F(1184)+p1mup5*p3mup3*p3
     -   mup4*F(1185)+p2mup5*p3mup3*p3mup4*F(1185)+p1mup4*p3mup3*p3mu
     -   p5*F(1185)+p2mup4*p3mup3*p3mup5*F(1185)+p3mup3*p3mup4*p3mup5
     -   *F(1186)+p3mup4*p3mup5*p6mup3*F(1187)+p1mup5*p4mup3*p4mup4*F
     -   (1188)+p2mup5*p4mup3*p4mup4*F(1188)+p4mup3*p4mup4*p4mup5*F(1
     -   189)+p3mup3*p3mup5*p6mup4*F(1190)
       Fa20(217) = mup3mup5*(p3mup4*F(1191)-p4mup4*F(1192))+mup3mup5*(
     -   p1mup4+p2mup4)*F(1195)+2*(p1mup3*p6mup4*p6mup5*F(1199)+p2mup
     -   3*p6mup4*p6mup5*F(1199)+p3mup3*p6mup4*p6mup5*F(1200)+p3mup4*
     -   p6mup3*p6mup5*F(1201)+p4mup4*p6mup3*p6mup5*F(1202)+p3mup3*p4
     -   mup5*p6mup4*F(1203)+p1mup3*p1mup5*p6mup4*F(1205)+p1mup5*p2mu
     -   p3*p6mup4*F(1205)+p1mup3*p2mup5*p6mup4*F(1205)+p2mup3*p2mup5
     -   *p6mup4*F(1205)+p1mup5*p3mup4*p6mup3*F(1206)+p2mup5*p3mup4*p
     -   6mup3*F(1206)+p1mup3*p3mup5*p6mup4*F(1206)+p2mup3*p3mup5*p6m
     -   up4*F(1206)+p1mup3*p4mup5*p6mup4*F(1207)+p2mup3*p4mup5*p6mup
     -   4*F(1207)+p1mup4*p6mup3*(p6mup5*F(1199)+(p1mup5+p2mup5)*F(12
     -   05)+p3mup5*F(1206)+p4mup5*F(1207))+p2mup4*p6mup3*(p6mup5*F(1
     -   199)+(p1mup5+p2mup5)*F(1205)+p3mup5*F(1206)+p4mup5*F(1207))+
     -   p3mup4*p4mup5*p6mup3*F(1208)-p6mup3*p6mup4*p6mup5*F(1209)-p3
     -   mup5*p4mup4*p6mup3*F(1210)+p4mup4*p4mup5*p6mup3*F(1211)+p1mu
     -   p5*p3mup3*p6mup4*F(1212)+p2mup5*p3mup3*p6mup4*F(1212)+p1mup5
     -   *p4mup4*p6mup3*F(1213)+p2mup5*p4mup4*p6mup3*F(1213)+p4mup3*(
     -   -((p1mup4+p2mup4)*((p1mup5+p2mup5)*F(1164)+p3mup5*F(1193)))+
     -   (p1mup5+p2mup5)*p3mup4*F(1194)-p6mup5*(p4mup4*F(1196)+(p1mup
     -   4+p2mup4)*F(1197)+p3mup4*F(1198))+p6mup4*(p6mup5*F(1118)-p3m
     -   up5*F(1204)+p4mup5*F(1211)+(p1mup5+p2mup5)*F(1213)))+p1mup5*
     -   p6mup3*p6mup4*F(1214)+p2mup5*p6mup3*p6mup4*F(1214)+p3mup5*p6
     -   mup3*p6mup4*F(1215)+p4mup5*p6mup3*p6mup4*F(1216))
       Fa20(218) = p1mup2*p1mup5*p4mup4*F(1221)+p1mup5*p2mup2*p4mup4*F
     -   (1221)+p1mup2*p2mup5*p4mup4*F(1221)+p2mup2*p2mup5*p4mup4*F(1
     -   221)-p1mup2*p3mup4*p6mup5*F(1222)-p2mup2*p3mup4*p6mup5*F(122
     -   2)+p1mup2*p4mup4*p6mup5*F(1224)+p2mup2*p4mup4*p6mup5*F(1224)
     -   -8*(mup4mup5*(p1mup2+p2mup2)*F(1217)-mup4mup5*p4mup2*s12*F(1
     -   218)+mup4mup5*p3mup2*s12*F(1225))+p1mup2*p1mup5*p3mup4*F(122
     -   7)+p1mup5*p2mup2*p3mup4*F(1227)+p1mup2*p2mup5*p3mup4*F(1227)
     -   +p2mup2*p2mup5*p3mup4*F(1227)+p1mup2*p3mup4*p4mup5*F(1228)+p
     -   2mup2*p3mup4*p4mup5*F(1228)+p1mup2*p1mup5*p6mup4*F(1229)+p1m
     -   up5*p2mup2*p6mup4*F(1229)+p1mup2*p2mup5*p6mup4*F(1229)+p2mup
     -   2*p2mup5*p6mup4*F(1229)+p1mup2*p3mup5*p6mup4*F(1230)+p2mup2*
     -   p3mup5*p6mup4*F(1230)+p1mup2*p6mup4*p6mup5*F(1231)+p2mup2*p6
     -   mup4*p6mup5*F(1231)+p1mup2*p4mup5*p6mup4*F(1232)+p2mup2*p4mu
     -   p5*p6mup4*F(1232)+p1mup4*((p1mup2+p2mup2)*(p4mup5*F(1220)+p6
     -   mup5*F(1233)-16*(p1mup5+p2mup5)*F(1239))-16*(p3mup2*s12*(p6m
     -   up5*F(945)+(p1mup5+p2mup5)*F(948)-p4mup5*F(1238))+p3mup5*(p4
     -   mup2*s12*F(913)-(p1mup2+p2mup2)*F(1226)+p3mup2*s12*F(1240)))
     -   )+p2mup4*((p1mup2+p2mup2)*(p4mup5*F(1220)+p6mup5*F(1233)-16*
     -   (p1mup5+p2mup5)*F(1239))-16*(p3mup2*s12*(p6mup5*F(945)+(p1mu
     -   p5+p2mup5)*F(948)-p4mup5*F(1238))+p3mup5*(p4mup2*s12*F(913)-
     -   (p1mup2+p2mup2)*F(1226)+p3mup2*s12*F(1240))))+16*(-(p3mup5*p
     -   4mup2*p4mup4*s12*F(915))+p1mup5*p3mup4*p4mup2*s12*F(950)+p2m
     -   up5*p3mup4*p4mup2*s12*F(950)+p3mup2*p6mup4*p6mup5*s12*F(960)
     -   -p1mup5*p3mup2*p6mup4*s12*F(966)-p2mup5*p3mup2*p6mup4*s12*F(
     -   966)+p3mup4*p4mup2*p4mup5*s12*F(991)+p3mup4*p4mup2*p6mup5*s1
     -   2*F(1146)+p3mup5*p4mup2*p6mup4*s12*F(1152)+p1mup2*p3mup4*p3m
     -   up5*F(1219)+p2mup2*p3mup4*p3mup5*F(1219)+p1mup2*p3mup5*p4mup
     -   4*F(1223)+p2mup2*p3mup5*p4mup4*F(1223)-p3mup2*p3mup4*p6mup5*
     -   s12*F(1234)-p3mup2*p4mup4*p6mup5*s12*F(1235)-p1mup5*p3mup2*p
     -   3mup4*s12*F(1236)-p2mup5*p3mup2*p3mup4*s12*F(1236)-p1mup5*p3
     -   mup2*p4mup4*s12*F(1237)-p2mup5*p3mup2*p4mup4*s12*F(1237)-p3m
     -   up2*p3mup4*p3mup5*s12*F(1241)-p3mup2*p3mup5*p4mup4*s12*F(124
     -   2)+p3mup2*p3mup4*p4mup5*s12*F(1243)-p1mup2*p4mup4*p4mup5*F(1
     -   244)-p2mup2*p4mup4*p4mup5*F(1244)-p3mup2*p4mup4*p4mup5*s12*F
     -   (1245)+p3mup4*p3mup5*p4mup2*s12*F(1246)+p3mup2*p4mup5*p6mup4
     -   *s12*F(1247)-p3mup2*p3mup5*p6mup4*s12*F(1248))
       Fa20(219) = 4*mup4mup5*F(1249)+(p1mup4+p2mup4)*p4mup5*F(1251)+(
     -   p1mup5+p2mup5)*p4mup4*F(1252)+p4mup4*p6mup5*F(1255)+p3mup4*(
     -   8*p3mup5*F(1250)+p6mup5*F(1253)+(p1mup5+p2mup5)*F(1257)+p4mu
     -   p5*F(1258))+(p1mup5+p2mup5)*p6mup4*F(1259)+p3mup5*(-8*(p4mup
     -   4*F(1254)-(p1mup4+p2mup4)*F(1256))+p6mup4*F(1260))+p6mup4*p6
     -   mup5*F(1261)+p4mup5*p6mup4*F(1262)+p1mup4*p6mup5*F(1263)+p2m
     -   up4*p6mup5*F(1263)+8*(p1mup4*p1mup5*F(1264)+p1mup5*p2mup4*F(
     -   1264)+p1mup4*p2mup5*F(1264)+p2mup4*p2mup5*F(1264)+p4mup4*p4m
     -   up5*F(1265))
       Fa20(220) = p4mup2*(mup4mup5*F(1266)-2*(-((p1mup5+p2mup5)*p3mup
     -   4*F(1268))-p4mup4*p6mup5*F(1270)+p1mup4*(p3mup5*F(1267)+p4mu
     -   p5*F(1269)+(p1mup5+p2mup5)*F(1271)-p6mup5*F(1272))+p2mup4*(p
     -   3mup5*F(1267)+p4mup5*F(1269)+(p1mup5+p2mup5)*F(1271)-p6mup5*
     -   F(1272))-p3mup4*p6mup5*F(1273)+p4mup5*p6mup4*F(1274)-(p1mup5
     -   +p2mup5)*p6mup4*F(1275)+p6mup4*p6mup5*F(1276)+p3mup4*p3mup5*
     -   F(1277)+p3mup4*p4mup5*F(1278)+p3mup5*p4mup4*F(1279)+p4mup4*p
     -   4mup5*F(1280)+p1mup5*p4mup4*F(1281)+p2mup5*p4mup4*F(1281)-p3
     -   mup5*p6mup4*F(1282)))
       Return
       End
