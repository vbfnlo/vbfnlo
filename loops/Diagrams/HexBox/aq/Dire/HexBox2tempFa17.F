c************** Calling the Fa functions*************************
       subroutine HexBox2Fa17(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup
     &   6,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa17)
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
       Complex*16 Fa17(305:323)
       COMMON/HexBox2FaFunctions/Fa
       Fa17(305) = -(p2mup3*F(1166))+p4mup3*F(1167)-p1mup3*F(1168)+p5m
     &   up3*F(1169)+p6mup3*F(1170)
       Fa17(306) = p2mup3*F(990)+p4mup3*F(1137)+p6mup3*F(1138)-p5mup3*
     &   F(1139)+p1mup3*F(1171)
       Fa17(307) = p2mup3*F(989)+p1mup3*F(1132)-p6mup3*F(1133)-p4mup3*
     &   F(1134)-p5mup3*F(1135)
       Fa17(308) = p4mup3*F(1108)+p1mup3*F(1173)+p2mup3*F(1174)+p6mup3
     &   *F(1175)+p5mup3*F(1176)
       Fa17(309) = 4*mup3mup6*F(1177)+8*(-(p1mup3*p4mup6*F(1105))-p2mu
     &   p6*(p4mup3*F(1178)+p6mup3*F(1180))+p4mup6*p6mup3*F(1181)+p4m
     &   up3*p6mup6*F(1181)+p6mup3*p6mup6*F(1182)+p5mup3*(p5mup6*F(10
     &   68)+p6mup6*F(1069)+p4mup6*F(1070)-p2mup6*F(1179)+p1mup6*F(11
     &   83))+p1mup6*p6mup3*F(1184)+p1mup6*p4mup3*F(1185)+p4mup3*p4mu
     &   p6*F(1186)+p1mup3*p2mup6*F(1187)+p5mup6*(p6mup3*F(1069)+p4mu
     &   p3*F(1070)-p1mup3*F(1188))-p1mup3*p6mup6*F(1189)-p1mup3*p1mu
     &   p6*F(1190)-p2mup3*(p5mup6*F(941)-p1mup6*F(1187)-p2mup6*F(119
     &   1)+p6mup6*F(1192)+p4mup6*F(1193)))
       Fa17(310) = -(p2mup6*F(1195))+p4mup6*F(1196)+p5mup6*F(1197)+p6m
     &   up6*F(1198)+p1mup6*F(1199)
       Fa17(311) = p4mup6*p5mup3*F(1116)+p4mup3*p5mup6*F(1116)+p4mup3*
     &   p4mup6*F(1118)+p5mup3*p5mup6*F(1119)+p5mup6*p6mup3*F(1120)+p
     &   5mup3*p6mup6*F(1120)-p6mup3*p6mup6*F(1121)-mup3mup6*F(1128)-
     &   p1mup3*p2mup6*F(1200)+p1mup3*p5mup6*F(1201)+p1mup3*p6mup6*F(
     &   1202)+p1mup3*p4mup6*F(1203)+p2mup6*p4mup3*F(1205)+p2mup6*p5m
     &   up3*F(1206)+p2mup6*p6mup3*F(1207)-(p4mup6*p6mup3+p4mup3*p6mu
     &   p6)*F(1208)+p1mup6*(p4mup3*F(1095)+p5mup3*F(1097)+p6mup3*F(1
     &   098)+p1mup3*F(1209))+p2mup3*(p5mup6*F(1071)-p1mup6*F(1193)-p
     &   2mup6*F(1204)+p6mup6*F(1210)+p4mup6*F(1211))
       Fa17(312) = p2mup3*F(1212)-p4mup3*F(1213)-p5mup3*F(1214)+p6mup3
     &   *F(1215)-p1mup3*F(1216)
       Fa17(313) = p6mup3*F(1217)+p4mup3*F(1218)+p1mup3*F(1219)+p5mup3
     &   *F(1220)+p2mup3*F(1221)
       Fa17(314) = p5mup3*F(1134)-p2mup3*F(1223)-p1mup3*F(1224)+p6mup3
     &   *F(1225)+p4mup3*F(1226)
       Fa17(315) = 4*(p4mup3*F(1227)+p5mup3*F(1228)-p6mup3*F(1229)-p2m
     &   up3*F(1230))+p1mup3*F(1231)
       Fa17(316) = -(p1mup6*F(1233))-p2mup6*F(1234)+p4mup6*F(1235)+p5m
     &   up6*F(1236)+p6mup6*F(1237)
       Fa17(317) = p1mup6*F(1238)+p5mup6*F(1239)+p4mup6*F(1240)+p6mup6
     &   *F(1241)+p2mup6*F(1242)
       Fa17(318) = p1mup6*p6mup3*F(1105)+p1mup6*p4mup3*F(1107)+p2mup6*
     &   p6mup3*F(1193)-p1mup6*p2mup3*F(1243)-p1mup3*p2mup6*F(1243)+p
     &   1mup3*p6mup6*F(1245)+p1mup3*p1mup6*F(1246)+p1mup3*p4mup6*F(1
     &   247)-p2mup3*p2mup6*F(1248)+p2mup3*p6mup6*F(1249)+p2mup3*p4mu
     &   p6*F(1250)-mup3mup6*F(1251)+p2mup6*p4mup3*F(1252)+p5mup3*(p6
     &   mup6*F(914)+p2mup6*F(942)+p1mup6*F(1104)+p4mup6*F(1253))+p4m
     &   up6*p6mup3*F(1254)+p4mup3*p6mup6*F(1254)+p5mup6*(p6mup3*F(91
     &   4)+p2mup3*F(938)+p1mup3*F(1244)+p4mup3*F(1253)+p5mup3*F(1255
     &   ))+p6mup3*p6mup6*F(1256)+p4mup3*p4mup6*F(1257)
       Fa17(319) = p2mup6*F(1258)+p4mup6*F(1259)+p5mup6*F(1260)+p6mup6
     &   *F(1261)+p1mup6*F(1262)
       Fa17(320) = p4mup6*F(1264)+2*(p2mup6*F(1263)+p5mup6*F(1265)+p6m
     &   up6*F(1266)+p1mup6*F(1267))
       Fa17(321) = p2mup6*F(1268)+2*(-(p1mup6*F(1269))+p4mup6*F(1270)+
     &   p5mup6*F(1271)+p6mup6*F(1644))
       Fa17(322) = p2mup6*F(1272)+p4mup6*F(1273)-2*p5mup6*F(1274)+p6mu
     &   p6*F(1275)+p1mup6*F(1276)
       Fa17(323) = p5mup6*p6mup3*F(1069)+p4mup3*p5mup6*F(1070)-p1mup6*
     &   p4mup3*F(1105)-p2mup6*p4mup3*F(1178)-p2mup6*p6mup3*F(1180)+p
     &   4mup6*p6mup3*F(1181)+p4mup3*p6mup6*F(1181)+p6mup3*p6mup6*F(1
     &   182)+p4mup3*p4mup6*F(1186)+p5mup3*(p5mup6*F(1068)+p6mup6*F(1
     &   069)+p4mup6*F(1070)-p2mup6*F(1179)-p1mup6*F(1188))-p1mup6*p6
     &   mup3*F(1189)-p2mup3*p2mup6*F(1277)-p2mup3*p5mup6*F(1278)-p2m
     &   up3*p6mup6*F(1279)-p2mup3*p4mup6*F(1280)-p1mup3*p2mup6*F(128
     &   1)-p1mup3*p4mup6*F(1282)-p1mup3*p5mup6*F(1283)+mup3mup6*F(12
     &   84)-p1mup3*p6mup6*F(1285)+p1mup6*p2mup3*F(1286)-p1mup3*p1mup
     &   6*F(1287)
       Return
       End
