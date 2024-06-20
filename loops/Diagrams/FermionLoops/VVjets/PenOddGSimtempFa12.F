c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenOddGSim12(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4
     -   ,p3mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5m
     -   up3,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mu
     -   p3,mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa12)
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
       Complex*16 Fa12(1090:1188)
       COMMON/FaPenOddGSim/Fa
       Fa12(1090) = -(p2mup1*F(2583))+p3mup1*F(2602)+p4mup1*F(2645)
       Fa12(1091) = 2*p1mup5*F(2611)+p4mup5*F(2612)+p3mup5*F(2623)+p2m
     -   up5*F(2626)
       Fa12(1092) = p4mup5*F(2595)+2*p2mup5*F(2619)+p3mup5*F(2620)+p1m
     -   up5*F(2626)
       Fa12(1093) = p2mup5*F(2620)+p1mup5*F(2623)+2*(p4mup5*F(2614)+p3
     -   mup5*F(2624))
       Fa12(1094) = p2mup5*F(2595)+p1mup5*F(2612)+2*(p4mup5*F(2613)+p3
     -   mup5*F(2614))
       Fa12(1095) = p1mup4*F(2618)-p2mup4*F(2620)+2*(p4mup4*F(2597)+p3
     -   mup4*F(2621))
       Fa12(1096) = p4mup4*F(2595)-p1mup4*F(2617)+2*p2mup4*F(2619)+p3m
     -   up4*F(2620)
       Fa12(1097) = 2*p1mup3*F(2611)+p2mup3*F(2626)
       Fa12(1098) = 2*p2mup3*F(2619)+p1mup3*F(2626)
       Fa12(1099) = p2mup3*F(2620)+p1mup3*F(2623)
       Fa12(1100) = p2mup3*F(2595)+p1mup3*F(2612)
       Fa12(1101) = 2*p1mup4*F(2528)+p4mup4*F(2587)+p2mup4*F(2617)+p3m
     -   up4*F(2618)
       Fa12(1102) = p1mup4*F(2587)-p2mup4*F(2595)+2*(p4mup4*F(2596)+p3
     -   mup4*F(2597))
       Fa12(1103) = p4mup1*F(2616)+p3mup1*F(2627)+p2mup1*F(4460)
       Fa12(1104) = p2mup5*F(2585)+p3mup5*F(2610)+2*p1mup5*F(2611)+p4m
     -   up5*F(2612)
       Fa12(1105) = p1mup5*F(2588)+p3mup5*F(2589)-2*p2mup5*F(2603)-p4m
     -   up5*F(2604)
       Fa12(1106) = p1mup4*F(2588)+p3mup4*F(2589)+2*p2mup4*F(2590)+p4m
     -   up4*F(2591)
       Fa12(1107) = 2*p1mup4*F(2528)-p2mup4*F(2585)+p3mup4*F(2586)+p4m
     -   up4*F(2587)
       Fa12(1108) = p2mup1*F(2583)+p3mup1*F(2615)+p4mup1*F(2616)
       Fa12(1109) = 4*(2*p1mup5*(p2mup4+p3mup4)*F(2432)+p1mup5*p4mup4*
     -   F(2434)+(p2mup5+p3mup5)*p4mup4*F(2436)+p1mup4*(p4mup5*F(2433
     -   )+2*(p1mup5*F(2431)+(p2mup5+p3mup5)*F(2437)))+(p2mup4+p3mup4
     -   )*(p4mup5*F(2435)+2*(p2mup5+p3mup5)*F(2438))+p4mup4*p4mup5*F
     -   (2439))+mup4mup5*F(2652)
       Fa12(1110) = -4*(p2mup4*p4mup5*F(2442)+p2mup5*p4mup4*F(2443)+p1
     -   mup5*(p4mup4*F(2443)+2*(p2mup4*F(2441)+p3mup4*F(2444)))+p3mu
     -   p4*p4mup5*F(2445)+p3mup5*p4mup4*F(2446)+p1mup4*(p4mup5*F(244
     -   2)+2*((p1mup5+p2mup5)*F(2441)+p3mup5*F(2447)))+2*(p2mup4*p2m
     -   up5*F(2441)+p2mup5*p3mup4*F(2444)+p2mup4*p3mup5*F(2447)+p3mu
     -   p4*p3mup5*F(2448))+p4mup4*p4mup5*F(2449))+mup4mup5*F(2653)
       Fa12(1111) = mup4mup5*F(2654)-4*(p2mup5*p3mup4*F(2454)+p3mup4*p
     -   4mup5*F(2458)+p3mup4*p3mup5*F(2459)+p3mup5*p4mup4*F(2463)+2*
     -   (p2mup4*p4mup5*F(2460)+p4mup4*p4mup5*F(2462)+p2mup5*p4mup4*F
     -   (2464)+p2mup4*p2mup5*F(2465))+p1mup4*(p2mup5*F(2436)+p1mup5*
     -   F(2456)+p4mup5*F(2457)-p3mup5*F(2655))-p2mup4*p3mup5*F(2656)
     -   +p1mup5*(p3mup4*F(2452)+p4mup4*F(2461)-p2mup4*F(2657)))
       Fa12(1112) = p2mup5*F(2658)+p3mup5*F(2659)+p4mup5*F(2660)+2*p1m
     -   up5*F(2661)
       Fa12(1113) = p1mup5*F(2662)+p3mup5*F(2663)+p4mup5*F(2664)+2*p2m
     -   up5*F(2665)
       Fa12(1114) = p1mup5*F(2666)+p2mup5*F(2667)+p4mup5*F(2668)+2*p3m
     -   up5*F(2669)
       Fa12(1115) = p1mup5*F(2670)+p2mup5*F(2671)+p3mup5*F(2672)+2*p4m
     -   up5*F(2673)
       Fa12(1116) = p1mup4*F(2675)+p2mup4*F(2676)+p3mup4*F(2677)+p4mup
     -   4*F(2678)
       Fa12(1117) = 4*(p1mup5*(2*(p2mup4*F(2489)+p3mup4*F(2490))+p4mup
     -   4*F(2492))-p2mup4*p4mup5*F(2493)-p2mup5*p4mup4*F(2494)+p3mup
     -   4*p4mup5*F(2498)+p3mup5*p4mup4*F(2499)+p1mup4*(p4mup5*F(2491
     -   )+2*(p1mup5*F(2488)-p2mup5*F(2495)+p3mup5*F(2500)))-2*(p2mup
     -   4*p2mup5*F(2496)+p2mup5*p3mup4*F(2497)-p2mup4*p3mup5*F(2501)
     -   -p3mup4*p3mup5*F(2502))+p4mup4*p4mup5*F(2503))+mup4mup5*F(26
     -   79)
       Fa12(1118) = 4*(p1mup5*(p2mup4*F(2506)+p3mup4*F(2507))-p2mup5*p
     -   3mup4*F(2509)-p1mup4*(p3mup5*F(2499)+p2mup5*F(2505)+p1mup5*F
     -   (2510)+p4mup5*F(2511))+p2mup4*(p3mup5*F(2508)+p4mup5*F(2512)
     -   +p2mup5*F(2513))-p1mup5*p4mup4*F(2515)+p2mup5*p4mup4*F(2517)
     -   -2*(p3mup4*p4mup5*F(2514)+p4mup4*p4mup5*F(2516)+p3mup5*p4mup
     -   4*F(2518)+p3mup4*p3mup5*F(2519)))+mup4mup5*F(2680)
       Fa12(1119) = p1mup4*F(2681)+p3mup4*F(2682)+p2mup4*F(2683)+p4mup
     -   4*F(2684)
       Fa12(1120) = -4*(p1mup4*(p2mup5*F(2525)+p3mup5*F(2526))-p1mup5*
     -   p2mup4*F(2527)-p1mup5*p3mup4*F(2528)-p2mup4*p3mup5*F(2529)+p
     -   2mup5*p3mup4*F(2530)-p2mup4*p4mup5*F(2531)-p3mup4*p4mup5*F(2
     -   532)-p2mup4*p2mup5*F(2533)-p3mup4*p3mup5*F(2534)-p2mup5*p4mu
     -   p4*F(2535)-p3mup5*p4mup4*F(2536))+mup4mup5*F(2582)
       Fa12(1121) = p2mup3*F(2685)+2*p1mup3*F(2686)
       Fa12(1122) = p1mup3*F(2687)+2*p2mup3*F(2688)
       Fa12(1123) = p1mup3*F(2689)+p2mup3*F(2690)
       Fa12(1124) = p1mup3*F(2691)+p2mup3*F(2692)
       Fa12(1125) = p2mup4*F(2693)+p3mup4*F(2694)+2*p1mup4*F(2695)+p4m
     -   up4*F(2696)
       Fa12(1126) = p1mup4*F(2697)+p2mup4*F(2698)+p3mup4*F(2699)+2*p4m
     -   up4*F(2700)
       Fa12(1127) = p3mup1*F(2580)+p2mup1*F(2581)+p4mup1*F(2582)
       Fa12(1128) = p3mup1*F(2563)+p4mup1*F(2573)
       Fa12(1129) = p2mup1*F(2570)+p4mup1*F(2574)
       Fa12(1130) = p2mup1*F(2571)+p3mup1*F(2572)
       Fa12(1131) = mup3mup5*F(2581)+4*(p1mup3*p4mup5*F(2564)-(p2mup3+
     -   p3mup3)*p4mup5*F(2567)+p4mup3*(p1mup5*F(2527)+(p2mup5+p3mup5
     -   )*F(2566)+p4mup5*F(2707)))
       Fa12(1132) = mup3mup5*F(2563)-4*((p1mup5+p2mup5)*p4mup3*F(2561)
     -   +p3mup5*p4mup3*F(2708)+p4mup5*((p1mup3+p2mup3)*F(2559)-p3mup
     -   3*F(2562)+p4mup3*F(2709)))
       Fa12(1133) = mup3mup5*F(2710)+p1mup3*(p2mup5*F(2713)+4*(p1mup5*
     -   F(2711)+p3mup5*F(2715)-p4mup5*F(2720)))+p2mup3*(p1mup5*F(271
     -   4)+4*p3mup5*F(2717)+8*(p2mup5*F(2712)+p4mup5*F(2721)))+4*(p1
     -   mup5*p3mup3*F(2716)+p2mup5*p3mup3*F(2718)+p3mup3*p3mup5*F(27
     -   19)+p3mup3*p4mup5*F(2722)-p1mup5*p4mup3*F(2724)-2*(p4mup3*p4
     -   mup5*F(2723)-p2mup5*p4mup3*F(2725))+p3mup5*p4mup3*F(2726))
       Fa12(1134) = p3mup3*F(2727)+p4mup3*F(2728)
       Fa12(1135) = p3mup3*F(2729)+p4mup3*F(2730)
       Fa12(1136) = p4mup3*F(2731)+2*p3mup3*F(2732)
       Fa12(1137) = mup3mup5*F(2580)+4*(p1mup3*p4mup5*F(2575)-p4mup5*(
     -   p2mup3*F(2576)+p3mup3*F(2579))+p4mup3*(p1mup5*F(2528)-p2mup5
     -   *F(2578)+p3mup5*F(2733)+p4mup5*F(2734)))
       Fa12(1138) = p3mup3*F(2735)+2*p4mup3*F(2736)
       Fa12(1139) = mup3mup5*F(2737)+p2mup3*(p2mup5*F(2742)-4*(p4mup5*
     -   F(2720)+p1mup5*F(2739)-p3mup5*F(2743)))-4*(p2mup5*p4mup3*F(2
     -   724)+p1mup5*p3mup3*F(2741)-p2mup5*p3mup3*F(2744)+p3mup3*p3mu
     -   p5*F(2745)+p1mup3*(p2mup5*F(2738)+p3mup5*F(2740)+2*(p1mup5*F
     -   (2488)+p4mup5*F(2746)))+p3mup3*p4mup5*F(2747)+2*(p4mup3*p4mu
     -   p5*F(2748)+p1mup5*p4mup3*F(2749))+p3mup5*p4mup3*F(2750))
       Fa12(1140) = mup3mup5*F(2751)-4*(p4mup3*(2*(p4mup5*F(2514)+p3mu
     -   p5*F(2519))-p2mup5*F(2726)+p1mup5*F(2750))+p1mup3*(p4mup5*F(
     -   2747)+p1mup5*F(2752)+p2mup5*F(2753)+p3mup5*F(2754))-p2mup3*(
     -   p4mup5*F(2722)+p2mup5*F(2755)+p1mup5*F(2756)+p3mup5*F(2757))
     -   +p3mup3*(2*(p4mup5*F(2519)+p3mup5*F(2758))+p1mup5*F(2759)-p2
     -   mup5*F(2760)))
       Fa12(1141) = -8*(p2mup3*p2mup4*F(2766)+p2mup4*p3mup3*F(2766)+p2
     -   mup3*p3mup4*F(2766)+p3mup3*p3mup4*F(2766)+p1mup4*(p2mup3+p3m
     -   up3)*F(2767))+p2mup3*p4mup4*F(2768)+p3mup3*p4mup4*F(2768)-4*
     -   (-(p1mup3*(2*(p1mup4*F(2761)+(p2mup4+p3mup4)*F(2762))+p4mup4
     -   *F(2763)))+p1mup4*p4mup3*F(2764)+p2mup4*p4mup3*F(2765)+p3mup
     -   4*p4mup3*F(2765)-p4mup3*p4mup4*F(2769))
       Fa12(1142) = p1mup3*(2*(p1mup4*F(2770)+p2mup4*F(2771)+p3mup4*F(
     -   2772))+p4mup4*F(2773))-p1mup4*p4mup3*F(2774)-p2mup3*p4mup4*F
     -   (2778)+p2mup4*p4mup3*F(2779)-p3mup4*p4mup3*F(2780)-2*(p2mup3
     -   *p2mup4*F(2775)+p1mup4*p2mup3*F(2776)+p2mup3*p3mup4*F(2777)+
     -   p3mup3*p3mup4*F(2781)+p1mup4*p3mup3*F(2782)+p2mup4*p3mup3*F(
     -   2783))-p3mup3*p4mup4*F(2784)+p4mup3*p4mup4*F(2785)
       Fa12(1143) = p2mup3*p4mup4*F(2796)-4*(p2mup3*p2mup4*F(2791)+p2m
     -   up3*p3mup4*F(2792)+p2mup4*p3mup3*F(2793)-p3mup3*p3mup4*F(279
     -   4)-p1mup4*(p2mup3*F(2788)+p3mup3*F(2790)-2*p4mup3*F(2795))-p
     -   1mup3*(p2mup4*F(2787)+p3mup4*F(2789)+2*(p1mup4*F(2786)-p4mup
     -   4*F(2795)))-p2mup4*p4mup3*F(2797)+p3mup3*p4mup4*F(2798)+p3mu
     -   p4*p4mup3*F(2799))-8*p4mup3*p4mup4*F(2800)
       Fa12(1144) = -8*(p2mup3*p2mup4*F(2801)+p2mup3*p3mup4*F(2802)-p3
     -   mup3*p3mup4*F(2806)-p2mup4*p3mup3*F(2807))+p1mup4*(4*p4mup3*
     -   F(2804)-8*(p2mup3*F(2801)-p3mup3*F(2807)))-4*(p2mup3*p4mup4*
     -   F(2803)+p1mup3*(2*((p1mup4+p2mup4)*F(2801)+p3mup4*F(2802))+p
     -   4mup4*F(2803))-p2mup4*p4mup3*F(2804)-p3mup4*p4mup3*F(2805)-p
     -   3mup3*p4mup4*F(2808))+p4mup3*p4mup4*F(2809)
       Fa12(1145) = p2mup3*p3mup4*F(2816)+p2mup4*p3mup3*F(2817)+p3mup3
     -   *p3mup4*F(2818)+p1mup3*(p1mup4*F(2810)+p2mup4*F(2812)+p3mup4
     -   *F(2814)+p4mup4*F(2819))-p1mup4*(p2mup3*F(2813)+p3mup3*F(281
     -   5)+p4mup3*F(2820))+p3mup4*p4mup3*F(2822)-p3mup3*p4mup4*F(282
     -   4)+2*(p2mup3*p2mup4*F(2811)+p2mup4*p4mup3*F(2821)-p2mup3*p4m
     -   up4*F(2823)+p4mup3*p4mup4*F(2825))
       Fa12(1146) = p2mup4*p3mup3*F(2834)+4*(p1mup4*p3mup3*F(2833)+p1m
     -   up3*(p1mup4*F(2826)+p2mup4*F(2827)+p3mup4*F(2828)-p4mup4*F(2
     -   835))-p1mup4*p4mup3*F(2836)-p2mup4*p4mup3*F(2837)-p2mup3*(p2
     -   mup4*F(2829)-p1mup4*F(2830)+p3mup4*F(2831)-p4mup4*F(2839)))+
     -   8*(p3mup3*p3mup4*F(2832)-p3mup4*p4mup3*F(2838)-p4mup3*p4mup4
     -   *F(2840)+p3mup3*p4mup4*F(2841))
       Fa12(1147) = mup3mup4*p3mup5*F(2708)-mup3mup5*p3mup4*F(2842)+mu
     -   p3mup4*p4mup5*F(2843)+mup3mup5*p4mup4*F(2844)+mup4mup5*p4mup
     -   3*F(2845)
       Fa12(1148) = mup3mup4*p4mup5*F(2531)+mup3mup5*p4mup4*F(2535)-mu
     -   p3mup5*p3mup4*F(2565)+mup3mup4*p3mup5*F(2566)+mup4mup5*p4mup
     -   3*F(2846)
       Fa12(1149) = mup3mup4*p4mup5*F(2532)+mup3mup5*p4mup4*F(2536)+mu
     -   p3mup4*p3mup5*F(2733)-mup3mup5*p3mup4*F(2847)+mup4mup5*p4mup
     -   3*F(2848)
       Fa12(1150) = p3mup3*F(2640)+p4mup3*F(2849)
       Fa12(1151) = -(p3mup3*F(2636))+p4mup3*F(2850)
       Fa12(1152) = p4mup3*F(2534)+p3mup3*F(2594)
       Fa12(1153) = p4mup3*F(2596)+p3mup3*F(2597)
       Fa12(1154) = p4mup3*F(2587)+p3mup3*F(2618)
       Fa12(1155) = p4mup3*F(2595)+p3mup3*F(2620)
       Fa12(1156) = p4mup3*F(2597)+p3mup3*F(2621)
       Fa12(1157) = p3mup3*F(2648)+p4mup3*F(2649)
       Fa12(1158) = p3mup3*F(2636)+p4mup3*F(2637)
       Fa12(1159) = p3mup3*F(2600)+p4mup3*F(2601)
       Fa12(1160) = p3mup3*F(2633)+p4mup3*F(2634)
       Fa12(1161) = p3mup3*F(2640)+p4mup3*F(2641)
       Fa12(1162) = p3mup3*F(2643)+p4mup3*F(2644)
       Fa12(1163) = p4mup3*F(2612)+p3mup3*F(2623)
       Fa12(1164) = p4mup3*F(2614)+p3mup3*F(2624)
       Fa12(1165) = p4mup3*F(2613)+p3mup3*F(2614)
       Fa12(1166) = p3mup3*F(2648)+p4mup3*F(2851)
       Fa12(1167) = p3mup3*F(2636)+p4mup3*F(2852)
       Fa12(1168) = mup3mup5*F(2568)-4*(p1mup3*p4mup5*F(2564)-(p2mup3+
     -   p3mup3)*p4mup5*F(2567)+p4mup3*(p1mup5*F(2527)+(p2mup5+p3mup5
     -   )*F(2566)+p4mup5*F(2707)))
       Fa12(1169) = mup3mup5*F(2570)+4*(p4mup3*((p1mup5+p2mup5)*F(2561
     -   )+p3mup5*F(2708))+p4mup5*((p1mup3+p2mup3)*F(2559)-p3mup3*F(2
     -   562)+p4mup3*F(2709)))
       Fa12(1170) = mup3mup5*F(2853)-4*(p1mup5*p3mup3*F(2716)+p2mup5*p
     -   3mup3*F(2718)+p3mup3*p3mup5*F(2719)+p3mup3*p4mup5*F(2722)-p1
     -   mup5*p4mup3*F(2724)-2*(p4mup3*p4mup5*F(2723)-p2mup5*p4mup3*F
     -   (2725))+p3mup5*p4mup3*F(2726)-p1mup3*(-(p3mup5*F(2715))+p4mu
     -   p5*F(2720)+p1mup5*F(2854)+p2mup5*F(2855))+p2mup3*(p3mup5*F(2
     -   717)+2*(p2mup5*F(2712)+p4mup5*F(2721))-p1mup5*F(2856)))
       Fa12(1171) = p3mup3*F(2857)+p4mup3*F(2858)
       Fa12(1172) = p3mup3*F(2859)+p4mup3*F(2860)
       Fa12(1173) = p4mup3*F(2861)+2*p3mup3*F(2862)
       Fa12(1174) = mup3mup5*F(2569)+4*(p4mup5*(-(p1mup3*F(2575))+p2mu
     -   p3*F(2576)+p3mup3*F(2579))-p4mup3*(p1mup5*F(2528)-p2mup5*F(2
     -   578)+p3mup5*F(2733)+p4mup5*F(2734)))
       Fa12(1175) = p3mup3*F(2863)-2*p4mup3*F(2864)
       Fa12(1176) = mup3mup5*F(2865)+4*(p2mup5*p4mup3*F(2724)+p1mup5*p
     -   3mup3*F(2741)-p2mup5*p3mup3*F(2744)+p3mup3*p3mup5*F(2745)+p1
     -   mup3*(p2mup5*F(2738)+p3mup5*F(2740)+2*(p1mup5*F(2488)+p4mup5
     -   *F(2746)))+p3mup3*p4mup5*F(2747)+2*(p4mup3*p4mup5*F(2748)+p1
     -   mup5*p4mup3*F(2749))+p3mup5*p4mup3*F(2750)+p2mup3*(p4mup5*F(
     -   2720)+p1mup5*F(2739)-p3mup5*F(2743)+p2mup5*F(2866)))
       Fa12(1177) = mup3mup5*F(2867)+4*(p4mup3*(2*(p4mup5*F(2514)+p3mu
     -   p5*F(2519))-p2mup5*F(2726)+p1mup5*F(2750))+p1mup3*(p4mup5*F(
     -   2747)+p1mup5*F(2752)+p2mup5*F(2753)+p3mup5*F(2754))-p2mup3*(
     -   p4mup5*F(2722)+p2mup5*F(2755)+p3mup5*F(2757)-p1mup5*F(2868))
     -   )+p3mup3*(8*(p4mup5*F(2519)+p3mup5*F(2758))+4*p1mup5*F(2759)
     -   +p2mup5*F(2869))
       Fa12(1178) = 8*p1mup4*(p2mup3+p3mup3)*F(2767)-4*(p1mup3*(2*(p1m
     -   up4*F(2761)+(p2mup4+p3mup4)*F(2762))+p4mup4*F(2763))-p1mup4*
     -   p4mup3*F(2764)-(p2mup4+p3mup4)*(p4mup3*F(2765)+2*(p2mup3+p3m
     -   up3)*F(2766))-(p2mup3+p3mup3)*p4mup4*F(2870))+p4mup3*p4mup4*
     -   F(2871)
       Fa12(1179) = -(p1mup3*(2*(p1mup4*F(2770)+p2mup4*F(2771)+p3mup4*
     -   F(2772))+p4mup4*F(2773)))+p2mup3*p4mup4*F(2778)-p2mup4*p4mup
     -   3*F(2779)+p3mup4*p4mup3*F(2780)+p1mup4*(p4mup3*F(2774)+2*(p2
     -   mup3*F(2776)+p3mup3*F(2782)))+2*(p2mup3*p2mup4*F(2775)+p2mup
     -   3*p3mup4*F(2777)+p3mup3*p3mup4*F(2781)+p2mup4*p3mup3*F(2783)
     -   )+p3mup3*p4mup4*F(2784)-p4mup3*p4mup4*F(2785)
       Fa12(1180) = -4*(-(p2mup3*p2mup4*F(2791))-p2mup3*p3mup4*F(2792)
     -   -p2mup4*p3mup3*F(2793)+p3mup3*p3mup4*F(2794)+p1mup4*(p2mup3*
     -   F(2788)+p3mup3*F(2790)-2*p4mup3*F(2795))+p1mup3*(p2mup4*F(27
     -   87)+p3mup4*F(2789)+2*(p1mup4*F(2786)-p4mup4*F(2795)))-p3mup3
     -   *p4mup4*F(2798)-p3mup4*p4mup3*F(2799)-2*p4mup3*p4mup4*F(2800
     -   )-p2mup3*p4mup4*F(2872))+p2mup4*p4mup3*F(2873)
       Fa12(1181) = 8*(p2mup3*p2mup4*F(2801)+p2mup3*p3mup4*F(2802)-p3m
     -   up3*p3mup4*F(2806)-p2mup4*p3mup3*F(2807))+p1mup4*(-4*p4mup3*
     -   F(2804)+8*(p2mup3*F(2801)-p3mup3*F(2807)))+p3mup3*p4mup4*F(2
     -   874)+4*(p2mup3*p4mup4*F(2803)+p1mup3*(2*((p1mup4+p2mup4)*F(2
     -   801)+p3mup4*F(2802))+p4mup4*F(2803))-p2mup4*p4mup3*F(2804)-p
     -   3mup4*p4mup3*F(2805)+p4mup3*p4mup4*F(2875))
       Fa12(1182) = 8*(p2mup4*p4mup3*F(2821)+p4mup3*p4mup4*F(2825))+p1
     -   mup4*p4mup3*F(2878)+4*(p2mup4*p3mup3*F(2817)+p3mup3*p3mup4*F
     -   (2818)+p1mup3*(p1mup4*F(2810)+p2mup4*F(2812)+p3mup4*F(2814)+
     -   p4mup4*F(2819))+p3mup4*p4mup3*F(2822)+p2mup3*(p3mup4*F(2816)
     -   +2*(p2mup4*F(2811)-p4mup4*F(2823))+p1mup4*F(2876))+p1mup4*p3
     -   mup3*F(2877)+p3mup3*p4mup4*F(2879))
       Fa12(1183) = -(p1mup4*p3mup3*F(2833))-p1mup3*(p1mup4*F(2826)+p2
     -   mup4*F(2827)+p3mup4*F(2828)-p4mup4*F(2835))+p1mup4*p4mup3*F(
     -   2836)+p2mup4*p4mup3*F(2837)-2*(p3mup3*p3mup4*F(2832)-p3mup4*
     -   p4mup3*F(2838)-p4mup3*p4mup4*F(2840)+p3mup3*p4mup4*F(2841))+
     -   p2mup3*(p2mup4*F(2829)+p3mup4*F(2831)-p4mup4*F(2839)+p1mup4*
     -   F(2880))+p2mup4*p3mup3*F(2881)
       Fa12(1184) = mup1mup3*(p1mup5*F(2527)+p4mup5*F(2531)+(p2mup5+p3
     -   mup5)*F(2566))
       Fa12(1185) = mup1mup3*((p1mup5+p2mup5)*F(2561)+p3mup5*F(2708)+p
     -   4mup5*F(2843))
       Fa12(1186) = mup1mup3*(p1mup5*F(2528)+p4mup5*F(2532)-p2mup5*F(2
     -   578)+p3mup5*F(2733))
       Fa12(1187) = mup1mup5*F(2571)-4*p3mup1*((p1mup5+p2mup5)*F(2561)
     -   +p3mup5*F(2708)+p4mup5*F(2843))
       Fa12(1188) = p4mup1*((p1mup4+p2mup4)*F(2560)+p3mup4*F(2842)-p4m
     -   up4*F(2844))
       Return
       End
