c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenOddGSim27(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4
     -   ,p3mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5m
     -   up3,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mu
     -   p3,mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa27)
       IMPLICIT NONE
      Complex*16   p1mup1, p1mup2, p1mup3, p1mup4, p1mup5, p2mup1, 
     -          p2mup2, p2mup3, p2mup4, p2mup5, p3mup1, p3mup2, 
     -          p3mup3, p3mup4, p3mup5, p4mup1, p4mup2, p4mup3, 
     -          p4mup4, p4mup5, p5mup1, p5mup2, p5mup3, p5mup4, p5mup5
       Complex*16   mup1mup2, mup1mup3, mup1mup4, mup1mup5, mup2mup3, 
     -          mup2mup4, mup2mup5, mup3mup4, mup3mup5, mup4mup5
        common/FPenOddGSim/F
      Complex*16 EpsF(126)
       common/EPSFPenOddGSim/EPSF
       COMMON/InvariantsPenOddGSim/p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s3
     -   4,s45,s234,m1,m1sq,m3,m3sq,m4,m4sq
       Real*8   p1sq 
        Real*8   p2sq 
       Real*8   p3sq 
       Real*8   p4sq 
       Real*8   p5sq 
       Real*8   s12 
       Real*8   s23 
       Real*8   s34 
       Real*8   s45 
       Real*8   s234,m1,m1sq,m3,m3sq,m4,m4sq
      Complex*16 Fa(3196),F(4476)
      Real*8 P(1383) 
        COMMON/PPenOddGSim/P
       Complex*16 Fa27(2575:2673)
       COMMON/FaPenOddGSim/Fa
       Fa27(2575) = -2*mup1mup2*p1mup3*F(1251)+Fa(685)+Fa(769)+4*Fa(96
     -   8)
       Fa27(2576) = Fa(687)+Fa(770)+4*Fa(969)+2*mup1mup2*p1mup3*F(258)
     -   *P(7)
       Fa27(2577) = 2*(mup1mup2*p1mup3*F(628)+mup1mup3*p4mup2*F(1749))
     -   +Fa(692)+Fa(975)
       Fa27(2578) = Fa(759)+Fa(867)-4*(mup1mup2*p4mup3*p4mup5*F(277)+F
     -   a(929)-Fa(985)+Fa(1011))
       Fa27(2579) = Fa(768)+Fa(913)+Fa(967)+4*Fa(1004)
       Fa27(2580) = Fa(773)+4*(Fa(918)+Fa(971)+Fa(994)+Fa(1017))
       Fa27(2581) = Fa(774)-4*(Fa(921)-Fa(973)-Fa(1005)+Fa(1018))
       Fa27(2582) = Fa(775)+Fa(974)-4*(Fa(900)-Fa(1006)-Fa(1020))
       Fa27(2583) = Fa(1033)+Fa(1134)
       Fa27(2584) = Fa(1034)+Fa(1135)
       Fa27(2585) = Fa(1035)+Fa(1136)
       Fa27(2586) = Fa(1036)+Fa(1138)
       Fa27(2587) = Fa(1039)+Fa(1147)
       Fa27(2588) = Fa(1040)-Fa(1148)
       Fa27(2589) = Fa(1045)-Fa(1149)
       Fa27(2590) = mup3mup4*F(2553)+Fa(1141)
       Fa27(2591) = mup3mup4*F(2554)+4*Fa(1142)
       Fa27(2592) = mup3mup4*F(2555)+Fa(1143)
       Fa27(2593) = mup3mup4*F(2556)+Fa(1144)
       Fa27(2594) = mup3mup4*F(2557)-4*Fa(1145)
       Fa27(2595) = mup3mup4*F(2558)+Fa(1146)
       Fa27(2596) = Fa(1184)-Fa(1218)
       Fa27(2597) = Fa(1185)-Fa(1219)
       Fa27(2598) = Fa(1186)-Fa(1220)
       Fa27(2599) = mup1mup4*F(2563)-4*Fa(1188)
       Fa27(2600) = Fa(1189)-Fa(1216)
       Fa27(2601) = Fa(1190)-Fa(1217)
       Fa27(2602) = mup1mup3*F(2571)-4*(p3mup1*p4mup3*F(2845)+Fa(1193)
     -   )
       Fa27(2603) = mup1mup3*F(2572)+4*(p2mup1*p4mup3*F(2845)+Fa(1194)
     -   )
       Fa27(2604) = Fa(1195)-Fa(1221)
       Fa27(2605) = mup1mup4*F(2574)-4*Fa(1196)
       Fa27(2606) = mup1mup4*F(2573)+4*Fa(1197)
       Fa27(2607) = -(mup1mup4*p4mup3*F(2845))+Fa(1198)
       Fa27(2608) = mup1mup4*F(2580)+4*Fa(1199)
       Fa27(2609) = mup1mup4*F(2581)+4*Fa(1200)
       Fa27(2610) = mup1mup4*F(2582)-4*Fa(1201)
       Fa27(2611) = -(mup1mup4*p4mup3*F(2846))+Fa(1202)
       Fa27(2612) = -(mup1mup4*p4mup3*F(2848))+Fa(1203)
       Fa27(2613) = mup1mup3*F(2524)+4*(p4mup3*(p2mup1*F(2846)+p3mup1*
     -   F(2848))+Fa(1204))
       Fa27(2614) = mup1mup5*F(2568)-4*Fa(1209)
       Fa27(2615) = mup1mup5*F(2570)+4*Fa(1210)
       Fa27(2616) = mup1mup5*F(2569)-4*Fa(1211)
       Fa27(2617) = mup1mup3*F(2570)+4*Fa(1212)
       Fa27(2618) = mup1mup3*F(2568)+4*Fa(1213)
       Fa27(2619) = mup1mup3*F(2569)-4*Fa(1214)
       Fa27(2620) = Fa(1052)+Fa(1150)
       Fa27(2621) = Fa(1053)+Fa(1151)
       Fa27(2622) = Fa(1054)+2*Fa(1152)
       Fa27(2623) = Fa(1055)+2*Fa(1153)
       Fa27(2624) = Fa(1065)+Fa(1154)
       Fa27(2625) = Fa(1066)-Fa(1155)
       Fa27(2626) = Fa(1067)+2*Fa(1156)
       Fa27(2627) = -Fa(1077)+Fa(1157)
       Fa27(2628) = -Fa(1053)+Fa(1158)
       Fa27(2629) = Fa(1078)+2*Fa(1159)
       Fa27(2630) = Fa(1079)+2*Fa(1160)
       Fa27(2631) = Fa(1052)+Fa(1161)
       Fa27(2632) = Fa(1053)-Fa(1158)
       Fa27(2633) = Fa(1087)+2*Fa(1162)
       Fa27(2634) = Fa(1097)+Fa(1163)
       Fa27(2635) = Fa(1098)+Fa(1155)
       Fa27(2636) = Fa(1099)+2*Fa(1164)
       Fa27(2637) = Fa(1100)+2*Fa(1165)
       Fa27(2638) = -Fa(1077)+Fa(1166)
       Fa27(2639) = -Fa(1053)+Fa(1167)
       Fa27(2640) = Fa(1121)+Fa(1171)
       Fa27(2641) = Fa(1122)+Fa(1172)
       Fa27(2642) = Fa(1123)+Fa(1173)
       Fa27(2643) = Fa(1124)+Fa(1175)
       Fa27(2644) = mup3mup4*F(2701)+Fa(1178)
       Fa27(2645) = mup3mup4*F(2702)+4*Fa(1179)
       Fa27(2646) = mup3mup4*F(2703)+Fa(1180)
       Fa27(2647) = mup3mup4*F(2704)+Fa(1181)
       Fa27(2648) = mup3mup4*F(2705)+Fa(1182)
       Fa27(2649) = mup3mup4*F(2706)+4*Fa(1183)
       Fa27(2650) = Fa(1184)+Fa(1222)
       Fa27(2651) = mup1mup4*F(2570)+4*Fa(1188)
       Fa27(2652) = mup1mup3*F(2573)+4*(p3mup1*p4mup3*F(2845)+Fa(1193)
     -   )
       Fa27(2653) = mup1mup3*F(2574)-4*(p2mup1*p4mup3*F(2845)+Fa(1194)
     -   )
       Fa27(2654) = mup1mup4*F(2572)+4*Fa(1196)
       Fa27(2655) = mup1mup4*F(2571)-4*Fa(1197)
       Fa27(2656) = mup1mup4*F(2569)-4*Fa(1199)
       Fa27(2657) = mup1mup4*F(2568)-4*Fa(1200)
       Fa27(2658) = mup1mup4*F(2524)+4*Fa(1201)
       Fa27(2659) = mup1mup3*F(2582)-4*(p4mup3*(p2mup1*F(2846)+p3mup1*
     -   F(2848))-Fa(1208))
       Fa27(2660) = mup1mup5*F(2581)+4*Fa(1209)
       Fa27(2661) = mup1mup5*F(2563)-4*Fa(1210)
       Fa27(2662) = mup1mup5*F(2580)+4*Fa(1211)
       Fa27(2663) = mup1mup3*F(2563)-4*Fa(1212)
       Fa27(2664) = mup1mup3*F(2581)+4*Fa(1215)
       Fa27(2665) = mup1mup3*F(2580)+4*Fa(1214)
       Fa27(2666) = mup2mup4*F(2993)+Fa(1428)
       Fa27(2667) = mup2mup4*F(2994)+Fa(1429)
       Fa27(2668) = mup2mup4*F(2995)+Fa(1430)
       Fa27(2669) = mup2mup4*F(2996)+Fa(1431)
       Fa27(2670) = mup2mup4*F(2997)+Fa(1432)
       Fa27(2671) = mup2mup4*F(2998)-2*Fa(1433)
       Fa27(2672) = 2*(mup1mup5*p4mup2*F(3481)-mup2mup5*p4mup1*F(3529)
     -   )+Fa(1466)+4*Fa(1605)
       Fa27(2673) = 4*mup1mup5*p4mup2*F(3490)+2*mup2mup5*p3mup1*F(3529
     -   )+Fa(1467)+Fa(1606)
       Return
       End
