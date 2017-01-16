c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenOddGSim11(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4
     -   ,p3mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5m
     -   up3,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mu
     -   p3,mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa11)
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
       Complex*16 Fa11(991:1089)
       COMMON/FaPenOddGSim/Fa
       Fa11(991) = mup2mup5*p2mup1*(p1mup4*F(346)-p4mup4*F(362)+p3mup4
     -   *F(2188)+p2mup4*F(2318))
       Fa11(992) = (mup1mup5*mup2mup3-mup1mup3*mup2mup5+mup1mup2*mup3m
     -   up5)*(p1mup4*F(346)-p4mup4*F(362)+p3mup4*F(2188)+p2mup4*F(23
     -   18))
       Fa11(993) = mup1mup4*(p1mup3*p3mup2*F(1660)+8*(-(p3mup2*p3mup3*
     -   F(1306))+p3mup3*p4mup2*F(1307)+p3mup2*p4mup3*F(1307)-p4mup2*
     -   p4mup3*F(1308)+p2mup3*(p3mup2*F(1299)+p4mup2*F(2417)))-4*p1m
     -   up2*(p1mup3*F(1289)+p2mup3*F(1290)-p3mup3*F(1300)-p4mup3*F(2
     -   418))+p1mup3*p4mup2*F(2419))
       Fa11(994) = mup2mup3*p4mup1*(p1mup4*F(346)-p4mup4*F(362)+p3mup4
     -   *F(2188)+p2mup4*F(2318))-mup1mup4*(-(p1mup3*p4mup2*F(258))+p
     -   3mup3*p4mup2*F(282)+p1mup2*p4mup3*F(1262)-p2mup3*p4mup2*F(23
     -   38)+p3mup2*p4mup3*F(2368)+p4mup2*p4mup3*F(2416))
       Fa11(995) = 4*mup2mup3*p3mup1*(p1mup4*F(346)-p4mup4*F(362)+p3mu
     -   p4*F(2188)+p2mup4*F(2318))+mup1mup4*(4*(p1mup3*p3mup2*F(1263
     -   )+p2mup3*p3mup2*F(1265)-p1mup2*(p3mup3*F(1261)+2*(p1mup3*F(6
     -   33)+p2mup3*F(635)+p4mup3*F(1276))))+p3mup2*p3mup3*F(1646)+p3
     -   mup2*p4mup3*F(2420)+p3mup3*p4mup2*F(2422)-8*(p1mup3*p4mup2*F
     -   (1276)+p2mup3*p4mup2*F(2421)+p4mup2*p4mup3*F(2423)))
       Fa11(996) = -4*mup2mup3*p2mup1*(p1mup4*F(346)-p4mup4*F(362)+p3m
     -   up4*F(2188)+p2mup4*F(2318))+mup1mup4*(-4*p1mup2*p2mup3*F(149
     -   9)+p2mup3*p3mup2*F(1778)+p2mup3*p4mup2*F(2425)+8*(p1mup2*p3m
     -   up3*F(279)+p1mup2*p4mup3*F(1503)+p1mup3*(p3mup2*F(279)-p1mup
     -   2*F(882)+p4mup2*F(1503))-p3mup2*p3mup3*F(1507)+p3mup3*p4mup2
     -   *F(1508)+p3mup2*p4mup3*F(1508)-p4mup2*p4mup3*F(2426)))
       Fa11(997) = p4mup1*(-4*(p1mup4*p3mup5*F(346)-p3mup5*p4mup4*F(36
     -   2))+mup4mup5*F(2319))
       Fa11(998) = p4mup1*(4*p1mup4*(p3mup5*F(272)+p2mup5*F(641))+8*(p
     -   3mup5*p4mup4*F(308)+p2mup5*p4mup4*F(691))+mup4mup5*F(1830))
       Fa11(999) = p4mup1*(4*(p1mup4*p2mup5*F(346)-p2mup5*p4mup4*F(362
     -   ))+mup4mup5*F(2187))
       Fa11(1000) = p3mup1*(-(mup3mup5*(p1mup4*F(346)-p4mup4*F(362)+p3
     -   mup4*F(2188)+p2mup4*F(2318)))+mup4mup5*(p1mup3*F(258)-p3mup3
     -   *F(282)+p2mup3*F(2338)-p4mup3*F(2345)))
       Fa11(1001) = 4*mup2mup5*p3mup1*(-(p1mup3*F(258))+p3mup3*F(282)-
     -   p2mup3*F(2338)+p4mup3*F(2345))+mup1mup5*(4*(p1mup3*p3mup2*F(
     -   1263)+p2mup3*p3mup2*F(1265))-8*(p1mup3*p4mup2*F(1276)+p1mup2
     -   *(p2mup3*F(635)+p4mup3*F(1276)))+p3mup2*p4mup3*F(2420)+p3mup
     -   3*p4mup2*F(2422))
       Fa11(1002) = mup1mup5*(8*(p2mup3*p3mup2*F(1299)+p3mup3*p4mup2*F
     -   (1307)+p3mup2*p4mup3*F(1307))+p1mup3*p3mup2*F(1660)-4*(p1mup
     -   2*p2mup3*F(1290)-p1mup2*p4mup3*F(2418))+p1mup3*p4mup2*F(2419
     -   ))
       Fa11(1003) = mup1mup5*(-4*p1mup2*p2mup3*F(1499)+8*(p1mup3*p3mup
     -   2*F(279)+p1mup3*p4mup2*F(1503)+p1mup2*p4mup3*F(1503)+(p3mup3
     -   *p4mup2+p3mup2*p4mup3)*F(1508))+p2mup3*p3mup2*F(1778))+4*mup
     -   2mup5*p2mup1*(p1mup3*F(258)-p3mup3*F(282)+p2mup3*F(2338)-p4m
     -   up3*F(2345))
       Fa11(1004) = mup1mup4*(p1mup3*p3mup2*F(1301)-2*(-(p3mup2*p3mup3
     -   *F(1306))+p3mup3*p4mup2*F(1307)+p3mup2*p4mup3*F(1307)-p4mup2
     -   *p4mup3*F(1308)+p2mup3*(p3mup2*F(1299)+p4mup2*F(2417)))+p1mu
     -   p2*(p1mup3*F(1289)+p2mup3*F(1290)-p3mup3*F(1300)-p4mup3*F(24
     -   18))+p1mup3*p4mup2*F(2424))
       Fa11(1005) = -(mup2mup3*p3mup1*(p1mup4*F(346)-p4mup4*F(362)+p3m
     -   up4*F(2188)+p2mup4*F(2318)))+mup1mup4*(-(p1mup3*p3mup2*F(126
     -   3))-p2mup3*p3mup2*F(1265)+p3mup2*p3mup3*F(1274)+p1mup2*(p3mu
     -   p3*F(1261)+2*(p1mup3*F(633)+p2mup3*F(635)+p4mup3*F(1276)))+2
     -   *(p1mup3*p4mup2*F(1276)+p2mup3*p4mup2*F(2421)+p4mup2*p4mup3*
     -   F(2423))+p3mup2*p4mup3*F(2428)+p3mup3*p4mup2*F(2429))
       Fa11(1006) = mup2mup3*p2mup1*(p1mup4*F(346)-p4mup4*F(362)+p3mup
     -   4*F(2188)+p2mup4*F(2318))+mup1mup4*(p1mup2*p2mup3*F(1499)+p2
     -   mup3*p3mup2*F(1505)-2*(p1mup2*p3mup3*F(279)+p1mup2*p4mup3*F(
     -   1503)+p1mup3*(p3mup2*F(279)-p1mup2*F(882)+p4mup2*F(1503))-p3
     -   mup2*p3mup3*F(1507)+p3mup3*p4mup2*F(1508)+p3mup2*p4mup3*F(15
     -   08)-p4mup2*p4mup3*F(2426))+p2mup3*p4mup2*F(2427))
       Fa11(1007) = p3mup1*(p1mup3*p3mup5*F(293)+p1mup5*p3mup3*F(310)+
     -   2*p3mup3*p3mup5*F(312)+p1mup3*p1mup5*F(828))
       Fa11(1008) = p4mup1*(p1mup4*p4mup5*F(298)+p1mup5*p4mup4*F(313)+
     -   2*p4mup4*p4mup5*F(315)+p1mup4*p1mup5*F(684))
       Fa11(1009) = mup4mup5*p4mup1*(p1mup2*F(1262)+p3mup2*F(2368)+p4m
     -   up2*F(2369))
       Fa11(1010) = mup4mup5*p3mup1*(p1mup2*F(1262)+p3mup2*F(2368)+p4m
     -   up2*F(2369))
       Fa11(1011) = mup3mup5*p4mup1*(p1mup2*F(1262)+p3mup2*F(2368)+p4m
     -   up2*F(2369))
       Fa11(1012) = p1mup2*p2mup1*(p1mup5*F(678)+p2mup5*F(680)+p3mup5*
     -   F(1292)+p4mup5*F(1399))
       Fa11(1013) = mup3mup5*p3mup1*(p1mup2*F(1262)+p3mup2*F(2368)+p4m
     -   up2*F(2369))
       Fa11(1014) = mup3mup5*p2mup1*(p1mup2*F(1262)+p3mup2*F(2368)+p4m
     -   up2*F(2369))
       Fa11(1015) = mup4mup5*p2mup1*(p1mup2*F(1262)+p3mup2*F(2368)+p4m
     -   up2*F(2369))
       Fa11(1016) = mup4mup5*(mup1mup2*(p1mup3*F(258)-p3mup3*F(282)+p2
     -   mup3*F(2338)-p4mup3*F(2345))-mup1mup3*(p1mup2*F(1262)+p3mup2
     -   *F(2368)+p4mup2*F(2369)))+mup1mup5*(mup2mup4*(p1mup3*F(258)-
     -   p3mup3*F(282)+p2mup3*F(2338)-p4mup3*F(2345))-mup3mup4*(p1mup
     -   2*F(1262)+p3mup2*F(2368)+p4mup2*F(2369)))+mup1mup4*(mup2mup5
     -   *(-(p1mup3*F(258))+p3mup3*F(282)-p2mup3*F(2338)+p4mup3*F(234
     -   5))+mup3mup5*(p1mup2*F(1262)+p3mup2*F(2368)+p4mup2*F(2369)))
       Fa11(1017) = p4mup1*(mup2mup4*(-(p1mup3*F(258))+p3mup3*F(282)-p
     -   2mup3*F(2338)+p4mup3*F(2345))+mup3mup4*(p1mup2*F(1262)+p3mup
     -   2*F(2368)+p4mup2*F(2369)))
       Fa11(1018) = p3mup1*(mup2mup4*(-(p1mup3*F(258))+p3mup3*F(282)-p
     -   2mup3*F(2338)+p4mup3*F(2345))+mup3mup4*(p1mup2*F(1262)+p3mup
     -   2*F(2368)+p4mup2*F(2369)))
       Fa11(1019) = p2mup1*(mup2mup4*(p1mup3*F(258)-p3mup3*F(282)+p2mu
     -   p3*F(2338)-p4mup3*F(2345))-mup3mup4*(p1mup2*F(1262)+p3mup2*F
     -   (2368)+p4mup2*F(2369)))
       Fa11(1020) = p2mup1*(mup2mup4*(-(p1mup3*F(258))+p3mup3*F(282)-p
     -   2mup3*F(2338)+p4mup3*F(2345))+mup3mup4*(p1mup2*F(1262)+p3mup
     -   2*F(2368)+p4mup2*F(2369)))
       Fa11(1021) = mup4mup5*F(2430)-4*(2*p1mup5*(p2mup4+p3mup4)*F(243
     -   2)+p1mup5*p4mup4*F(2434)+(p2mup5+p3mup5)*p4mup4*F(2436)+p1mu
     -   p4*(p4mup5*F(2433)+2*(p1mup5*F(2431)+(p2mup5+p3mup5)*F(2437)
     -   ))+(p2mup4+p3mup4)*(p4mup5*F(2435)+2*(p2mup5+p3mup5)*F(2438)
     -   )+p4mup4*p4mup5*F(2439))
       Fa11(1022) = mup4mup5*F(2440)+4*(p2mup4*p4mup5*F(2442)+p1mup5*p
     -   4mup4*F(2443)+p2mup5*p4mup4*F(2443)+p3mup4*p4mup5*F(2445)+p3
     -   mup5*p4mup4*F(2446)+p1mup4*(p4mup5*F(2442)+2*((p1mup5+p2mup5
     -   )*F(2441)+p3mup5*F(2447)))+2*(p1mup5*p2mup4*F(2441)+p2mup4*p
     -   2mup5*F(2441)+p1mup5*p3mup4*F(2444)+p2mup5*p3mup4*F(2444)+p2
     -   mup4*p3mup5*F(2447)+p3mup4*p3mup5*F(2448))+p4mup4*p4mup5*F(2
     -   449))
       Fa11(1023) = mup4mup5*F(2450)+4*(p1mup5*p3mup4*F(2452)+p2mup5*p
     -   3mup4*F(2454)+p1mup4*(p2mup5*F(2436)+p3mup5*F(2451)+p1mup5*F
     -   (2456)+p4mup5*F(2457))+p3mup4*p4mup5*F(2458)+p3mup4*p3mup5*F
     -   (2459)+p1mup5*p4mup4*F(2461)+p3mup5*p4mup4*F(2463)+2*(p4mup4
     -   *p4mup5*F(2462)+p2mup5*p4mup4*F(2464))+p2mup4*(p3mup5*F(2453
     -   )+p1mup5*F(2455)+2*(p4mup5*F(2460)+p2mup5*F(2465))))
       Fa11(1024) = p2mup5*F(2466)+p3mup5*F(2467)+p4mup5*F(2468)+2*p1m
     -   up5*F(2469)
       Fa11(1025) = p1mup5*F(2470)+p3mup5*F(2471)+p4mup5*F(2472)+2*p2m
     -   up5*F(2473)
       Fa11(1026) = p1mup5*F(2474)+p2mup5*F(2475)+p4mup5*F(2476)+2*p3m
     -   up5*F(2477)
       Fa11(1027) = p1mup5*F(2478)+p2mup5*F(2479)+p3mup5*F(2480)+2*p4m
     -   up5*F(2481)
       Fa11(1028) = p1mup4*F(2483)+p2mup4*F(2484)+p3mup4*F(2485)+p4mup
     -   4*F(2486)
       Fa11(1029) = mup4mup5*F(2487)-4*(p1mup5*p4mup4*F(2492)-p2mup4*p
     -   4mup5*F(2493)-p2mup5*p4mup4*F(2494)+p3mup4*p4mup5*F(2498)+p3
     -   mup5*p4mup4*F(2499)+p1mup4*(p4mup5*F(2491)+2*(p1mup5*F(2488)
     -   -p2mup5*F(2495)+p3mup5*F(2500)))+2*(p1mup5*p2mup4*F(2489)+p1
     -   mup5*p3mup4*F(2490)-p2mup4*p2mup5*F(2496)-p2mup5*p3mup4*F(24
     -   97)+p2mup4*p3mup5*F(2501)+p3mup4*p3mup5*F(2502))+p4mup4*p4mu
     -   p5*F(2503))
       Fa11(1030) = mup4mup5*F(2504)-4*(p1mup5*(p2mup4*F(2506)+p3mup4*
     -   F(2507))-p2mup5*p3mup4*F(2509)-p1mup4*(p3mup5*F(2499)+p2mup5
     -   *F(2505)+p1mup5*F(2510)+p4mup5*F(2511))+p2mup4*(p3mup5*F(250
     -   8)+p4mup5*F(2512)+p2mup5*F(2513))-p1mup5*p4mup4*F(2515)+p2mu
     -   p5*p4mup4*F(2517)-2*(p3mup4*p4mup5*F(2514)+p4mup4*p4mup5*F(2
     -   516)+p3mup5*p4mup4*F(2518)+p3mup4*p3mup5*F(2519)))
       Fa11(1031) = p1mup4*F(2520)+p3mup4*F(2521)+p2mup4*F(2522)+p4mup
     -   4*F(2523)
       Fa11(1032) = mup4mup5*F(2524)+4*(p1mup4*p2mup5*F(2525)+p1mup4*p
     -   3mup5*F(2526)-p1mup5*(p2mup4*F(2527)+p3mup4*F(2528))-p2mup4*
     -   p3mup5*F(2529)+p2mup5*p3mup4*F(2530)-p2mup4*p4mup5*F(2531)-p
     -   3mup4*p4mup5*F(2532)-p2mup4*p2mup5*F(2533)-p3mup4*p3mup5*F(2
     -   534)-p2mup5*p4mup4*F(2535)-p3mup5*p4mup4*F(2536))
       Fa11(1033) = p2mup3*F(2537)-2*p1mup3*F(2538)
       Fa11(1034) = p1mup3*F(2539)+2*p2mup3*F(2540)
       Fa11(1035) = p1mup3*F(2541)+p2mup3*F(2542)
       Fa11(1036) = p1mup3*F(2543)+p2mup3*F(2544)
       Fa11(1037) = p2mup4*F(2545)+p3mup4*F(2546)-2*p1mup4*F(2547)+p4m
     -   up4*F(2548)
       Fa11(1038) = p1mup4*F(2549)+p2mup4*F(2550)+p3mup4*F(2551)+2*p4m
     -   up4*F(2552)
       Fa11(1039) = mup4mup5*(p1mup3+p2mup3)*F(2559)-mup3mup5*(p1mup4+
     -   p2mup4)*F(2560)+mup3mup4*(p1mup5+p2mup5)*F(2561)-mup4mup5*p3
     -   mup3*F(2562)
       Fa11(1040) = -(mup4mup5*p1mup3*F(2564))+mup3mup5*(p1mup4*F(2525
     -   )+p2mup4*F(2565))-mup3mup4*(p1mup5*F(2527)+p2mup5*F(2566))+m
     -   up4mup5*(p2mup3+p3mup3)*F(2567)
       Fa11(1041) = p4mup1*F(2524)+p2mup1*F(2568)+p3mup1*F(2569)
       Fa11(1042) = p3mup1*F(2570)+p4mup1*F(2571)
       Fa11(1043) = p2mup1*F(2563)+p4mup1*F(2572)
       Fa11(1044) = p2mup1*F(2573)+p3mup1*F(2574)
       Fa11(1045) = mup3mup5*p1mup4*F(2526)-mup3mup4*p1mup5*F(2528)-mu
     -   p4mup5*p1mup3*F(2575)+mup4mup5*p2mup3*F(2576)-mup3mup5*p2mup
     -   4*F(2577)+mup3mup4*p2mup5*F(2578)+mup4mup5*p3mup3*F(2579)
       Fa11(1046) = 2*p1mup5*F(2528)-p2mup5*F(2585)+p3mup5*F(2586)+p4m
     -   up5*F(2587)
       Fa11(1047) = p1mup5*F(2588)+p3mup5*F(2589)+2*p2mup5*F(2590)+p4m
     -   up5*F(2591)
       Fa11(1048) = p1mup5*F(2592)-p2mup5*F(2593)+2*(p4mup5*F(2534)+p3
     -   mup5*F(2594))
       Fa11(1049) = p1mup5*F(2587)-p2mup5*F(2595)+2*(p4mup5*F(2596)+p3
     -   mup5*F(2597))
       Fa11(1050) = p2mup4*F(2593)+p1mup4*F(2599)+2*(p3mup4*F(2600)+p4
     -   mup4*F(2601))
       Fa11(1051) = p1mup4*F(2588)+p3mup4*F(2589)-2*p2mup4*F(2603)-p4m
     -   up4*F(2604)
       Fa11(1052) = p2mup3*F(2606)+2*p1mup3*F(2607)
       Fa11(1053) = p1mup3*F(2608)-2*p2mup3*F(2609)
       Fa11(1054) = p1mup3*F(2592)-p2mup3*F(2593)
       Fa11(1055) = p1mup3*F(2587)-p2mup3*F(2595)
       Fa11(1056) = p2mup4*F(2585)+p3mup4*F(2610)+2*p1mup4*F(2611)+p4m
     -   up4*F(2612)
       Fa11(1057) = p2mup4*F(2595)+p1mup4*F(2612)+2*(p4mup4*F(2613)+p3
     -   mup4*F(2614))
       Fa11(1058) = -(p2mup1*F(2583))+p3mup1*F(2602)+p4mup1*F(2605)
       Fa11(1059) = p3mup1*F(2584)+p4mup1*F(4462)
       Fa11(1060) = 2*p1mup5*F(2528)+p4mup5*F(2587)+p2mup5*F(2617)+p3m
     -   up5*F(2618)
       Fa11(1061) = p4mup5*F(2595)-p1mup5*F(2617)+2*p2mup5*F(2619)+p3m
     -   up5*F(2620)
       Fa11(1062) = p1mup5*F(2618)-p2mup5*F(2620)+2*(p4mup5*F(2597)+p3
     -   mup5*F(2621))
       Fa11(1063) = p2mup4*F(2620)+p1mup4*F(2623)+2*(p4mup4*F(2614)+p3
     -   mup4*F(2624))
       Fa11(1064) = p4mup4*F(2595)+2*p2mup4*F(2619)+p3mup4*F(2620)+p1m
     -   up4*F(2626)
       Fa11(1065) = 2*p1mup3*F(2528)+p2mup3*F(2617)
       Fa11(1066) = p1mup3*F(2617)-2*p2mup3*F(2619)
       Fa11(1067) = p1mup3*F(2618)-p2mup3*F(2620)
       Fa11(1068) = 2*p1mup4*F(2611)+p4mup4*F(2612)+p3mup4*F(2623)+p2m
     -   up4*F(2626)
       Fa11(1069) = p4mup1*F(2605)+p3mup1*F(2625)-p2mup1*F(4460)
       Fa11(1070) = p3mup1*F(4461)+p4mup1*F(4462)
       Fa11(1071) = p2mup5*F(2585)+p3mup5*F(2610)+2*p1mup5*F(2611)+p4m
     -   up5*F(2629)
       Fa11(1072) = p1mup5*F(2588)+p3mup5*F(2589)-2*p2mup5*F(2603)+p4m
     -   up5*F(2630)
       Fa11(1073) = p2mup5*F(2593)+p1mup5*F(2599)+2*(p3mup5*F(2600)+p4
     -   mup5*F(2601))
       Fa11(1074) = p1mup5*F(2631)+p2mup5*F(2632)+2*(p3mup5*F(2633)+p4
     -   mup5*F(2634))
       Fa11(1075) = p1mup4*F(2592)-p2mup4*F(2593)+2*(p4mup4*F(2534)+p3
     -   mup4*F(2594))
       Fa11(1076) = -(p1mup4*F(2608))+2*p2mup4*F(2609)+p3mup4*F(2636)+
     -   p4mup4*F(2637)
       Fa11(1077) = p2mup3*F(2606)-2*p1mup3*F(2639)
       Fa11(1078) = p2mup3*F(2593)+p1mup3*F(2599)
       Fa11(1079) = p1mup3*F(2631)+p2mup3*F(2632)
       Fa11(1080) = p2mup4*F(2606)+2*p1mup4*F(2607)+p3mup4*F(2640)+p4m
     -   up4*F(2641)
       Fa11(1081) = -(p2mup4*F(2632))+p1mup4*F(2642)+2*(p3mup4*F(2643)
     -   +p4mup4*F(2644))
       Fa11(1082) = p2mup1*F(2583)+p3mup1*F(2615)+p4mup1*F(2638)
       Fa11(1083) = p3mup1*F(2584)+p4mup1*F(2628)
       Fa11(1084) = 2*p1mup5*F(2528)-p2mup5*F(2585)+p3mup5*F(2586)+p4m
     -   up5*F(2646)
       Fa11(1085) = p1mup5*F(2588)+p3mup5*F(2589)+2*p2mup5*F(2590)+p4m
     -   up5*F(2630)
       Fa11(1086) = -(p2mup5*F(2632))+p1mup5*F(2642)+2*(p3mup5*F(2643)
     -   +p4mup5*F(2644))
       Fa11(1087) = -(p2mup3*F(2632))+p1mup3*F(2642)
       Fa11(1088) = -(p2mup4*F(2606))+2*p1mup4*F(2639)+p3mup4*F(2648)+
     -   p4mup4*F(2649)
       Fa11(1089) = p1mup4*F(2631)+p2mup4*F(2632)+2*(p3mup4*F(2633)+p4
     -   mup4*F(2634))
       Return
       End
