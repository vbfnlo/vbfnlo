c************** Calling the Fa functions*************************
       subroutine NoAbe2gaHexLineFa22(p1mup2,p1mup3,p1mup5,p1mup6,p2mu
     -   p2,p2mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p
     -   4mup3,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup
     -   3,p6mup5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6
     -   ,mup5mup6,Fa22)
       IMPLICIT NONE
      Complex*16   p1mup2, p1mup3, p1mup5, p1mup6, p2mup2, p2mup3, 
     -          p2mup5, p2mup6, p3mup2, p3mup3, p3mup5, p3mup6, 
     -          p4mup2, p4mup3, p4mup5, p4mup6, p5mup2, p5mup3, 
     -          p5mup5, p5mup6, p6mup2, p6mup3, p6mup5, p6mup6
       Complex*16   mup2mup3, mup2mup5, mup2mup6, mup3mup5, mup3mup6, 
     -          mup5mup6
        common/NoAbe2gaHexLineFhlFunctions/F
       COMMON/NoAbe2gaHexLineInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,
     -   s12,s23,s34,s45,s56,s16,s123,s234,s345
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
      Complex*16 Fa(872),F(9706)
      Real*8 P(9687) 
        COMMON/NoAbe2gaHexLinePFunctions/P
       Complex*16 Fa22(505:528)
       COMMON/NoAbe2gaHexLineFaFunctions/Fa
       Fa22(505) = p5mup2*(8*(p1mup5*p4mup3*F(2320)+p2mup5*p4mup3*F(23
     -   20)+p1mup5*p5mup3*F(2341)+p2mup5*p5mup3*F(2341)+p1mup5*p3mup
     -   3*F(2346)+p2mup5*p3mup3*F(2346)+p3mup5*(p5mup3*F(2276)+p3mup
     -   3*F(2279)+p4mup3*F(2319)+(p1mup3+p2mup3)*F(2346))-p1mup3*p1m
     -   up5*F(2353)-p1mup5*p2mup3*F(2353)-p1mup3*p2mup5*F(2353)-p2mu
     -   p3*p2mup5*F(2353))-p5mup3*p5mup5*F(2631)-p4mup3*p5mup5*F(263
     -   2)-p5mup5*p6mup3*F(2635)+p3mup3*p5mup5*F(2639)+p1mup3*p5mup5
     -   *F(2642)+p2mup3*p5mup5*F(2642)-4*(p4mup3*p6mup5*F(2633)+p5mu
     -   p3*p6mup5*F(2634)+p6mup3*p6mup5*F(2636)+p4mup5*(2*(p5mup3*F(
     -   2307)+p4mup3*F(2308)+(p1mup3+p2mup3)*F(2309)-p3mup3*F(2310))
     -   +p6mup3*F(2637))+mup3mup5*F(2638)-p3mup3*p6mup5*F(2640)-p1mu
     -   p3*p6mup5*F(2643)-p2mup3*p6mup5*F(2643)+p3mup5*p6mup3*F(2645
     -   )+(p1mup5+p2mup5)*p6mup3*F(2646)))
       Fa22(506) = p1mup5*p2mup2*F(2648)+p2mup2*p2mup5*F(2648)+p2mup2*
     -   p5mup5*F(2649)+p2mup2*p6mup5*F(2650)+p2mup2*p3mup5*F(2651)+p
     -   2mup2*p4mup5*F(2652)+p1mup2*((p1mup5+p2mup5)*F(2648)+p5mup5*
     -   F(2649)+p6mup5*F(2650)+p3mup5*F(2651)+p4mup5*F(2652))-2*p5mu
     -   p5*p6mup2*s12*F(2653)+p6mup2*p6mup5*s12*F(2654)+4*(p1mup5*p6
     -   mup2*s12*F(2655)+p2mup5*p6mup2*s12*F(2655)+p3mup5*p6mup2*s12
     -   *F(2656)+p4mup5*p6mup2*s12*F(2657))
       Fa22(507) = (p1mup2+p2mup2)*((p1mup5+p2mup5)*F(2658)+p5mup5*F(2
     -   659)+p3mup5*F(2661)-2*(p6mup5*F(2660)-p4mup5*F(2662)))+p5mup
     -   2*s12*(2*p5mup5*F(2653)-4*((p1mup5+p2mup5)*F(2655)+p3mup5*F(
     -   2656)+p4mup5*F(2657))+p6mup5*F(2663))
       Fa22(508) = (p1mup2+p2mup2)*(p4mup5*p6mup6*F(2588)-4*(mup5mup6*
     -   F(2581)-(p1mup6+p2mup6)*(p5mup5*F(2583)-4*(p3mup5*F(2584)+p4
     -   mup5*F(2586)))+p5mup5*p6mup6*F(2593))-8*(p2mup5*p6mup6*F(258
     -   2)+p4mup6*p5mup5*F(2590)-p5mup5*p5mup6*F(2596)-p1mup5*(-(p6m
     -   up6*F(2582))+2*(p4mup6*F(2589)+p5mup6*F(2592)+(p1mup6+p2mup6
     -   )*F(2597)))+p6mup5*p6mup6*F(9579)-p3mup5*p6mup6*F(9588)-p3mu
     -   p6*(p5mup5*F(2218)+2*((p1mup5+p2mup5)*F(2585)-p4mup5*F(2599)
     -   +p6mup5*F(9588)-p3mup5*F(9590))))+16*(p4mup5*p5mup6*F(2587)+
     -   p2mup5*p4mup6*F(2589)-p4mup6*p6mup5*F(2591)+p2mup5*p5mup6*F(
     -   2592)-p1mup6*p6mup5*F(2594)-p2mup6*p6mup5*F(2594)-p4mup5*p4m
     -   up6*F(2595)+p1mup6*p2mup5*F(2597)+p2mup5*p2mup6*F(2597)+p3mu
     -   p5*p4mup6*F(2598)+p5mup6*p6mup5*F(9583)-p3mup5*p5mup6*F(9594
     -   )))
       Fa22(509) = (p1mup2+p2mup2)*(p6mup3*p6mup6*F(2608)-16*(-(p5mup6
     -   *p6mup3*F(2606))+p3mup3*p5mup6*F(2610)+p1mup3*p5mup6*F(2613)
     -   +p2mup3*p5mup6*F(2613))+4*(-(mup3mup6*F(2602))+p4mup3*(8*(p4
     -   mup6*F(1967)-p3mup6*F(2214)-(p1mup6+p2mup6)*F(2228))+4*p5mup
     -   6*F(2603)+p6mup6*F(2604))+p3mup3*p6mup6*F(2612)+p5mup3*(p6mu
     -   p6*F(2607)+4*(p5mup6*F(2605)+p4mup6*F(2609)+p3mup6*F(2611)+(
     -   p1mup6+p2mup6)*F(2614)))+p1mup3*p6mup6*F(2615)+p2mup3*p6mup6
     -   *F(2615))-32*(p3mup6*p6mup3*F(9581)+p1mup6*p6mup3*F(9584)+p2
     -   mup6*p6mup3*F(9584)+p4mup6*p6mup3*F(9587)-p1mup3*p4mup6*F(95
     -   98)-p2mup3*p4mup6*F(9598)-p3mup3*p3mup6*F(9599)-p1mup3*p1mup
     -   6*F(9603)-p1mup6*p2mup3*F(9603)-p1mup3*p2mup6*F(9603)-p2mup3
     -   *p2mup6*F(9603)-p1mup6*p3mup3*F(9604)-p2mup6*p3mup3*F(9604)-
     -   p1mup3*p3mup6*F(9604)-p2mup3*p3mup6*F(9604)-p3mup3*p4mup6*F(
     -   9606)))
       Fa22(510) = (p1mup2+p2mup2)*(32*p3mup5*(p5mup3*F(2276)+p3mup3*F
     -   (2279)+p4mup3*F(2319)+(p1mup3+p2mup3)*F(2346))-16*p4mup5*(2*
     -   (p5mup3*F(2307)+p4mup3*F(2308)+(p1mup3+p2mup3)*F(2309)-p3mup
     -   3*F(2310))+p6mup3*F(2637))+p3mup5*p6mup3*F(2641)+4*(8*(p1mup
     -   5*(p4mup3*F(2320)+p5mup3*F(2341)+p3mup3*F(2346)-(p1mup3+p2mu
     -   p3)*F(2353))+p2mup5*(p4mup3*F(2320)+p5mup3*F(2341)+p3mup3*F(
     -   2346)-(p1mup3+p2mup3)*F(2353)))-p5mup3*p5mup5*F(2631)-p4mup3
     -   *p5mup5*F(2632)-p5mup5*p6mup3*F(2635)+p3mup3*p5mup5*F(2639)+
     -   p1mup3*p5mup5*F(2642)+p2mup3*p5mup5*F(2642)-4*(p4mup3*p6mup5
     -   *F(2633)+p5mup3*p6mup5*F(2634)+p6mup3*p6mup5*F(2636)+mup3mup
     -   5*F(2638)-p3mup3*p6mup5*F(2640)-(p1mup3+p2mup3)*p6mup5*F(264
     -   3)))+(p1mup5+p2mup5)*p6mup3*F(2644))
       Fa22(511) = -2*(p2mup6*p5mup5*F(2583)-p5mup5*p6mup6*F(2593)-p1m
     -   up6*(-(p5mup5*F(2583))+4*(p3mup5*F(2584)+p4mup5*F(2586)+p6mu
     -   p5*F(2594)-(p1mup5+p2mup5)*F(2597))))+mup5mup6*F(2664)+p4mup
     -   5*p6mup6*F(2666)+4*((p1mup5+p2mup5)*p6mup6*F(2582)+p4mup6*p5
     -   mup5*F(2590)-p5mup5*p5mup6*F(2596)+p6mup5*p6mup6*F(9579)-p3m
     -   up5*p6mup6*F(9588)-p3mup6*(p5mup5*F(2218)+2*((p1mup5+p2mup5)
     -   *F(2585)-p4mup5*F(2599)+p6mup5*F(9588)-p3mup5*F(9590))))+8*(
     -   p2mup6*p3mup5*F(2584)+p2mup6*p4mup5*F(2586)+p4mup6*p6mup5*F(
     -   2591)-(p1mup5+p2mup5)*(p4mup6*F(2589)+p5mup6*F(2592))+p2mup6
     -   *p6mup5*F(2594)+p4mup5*p4mup6*F(2595)-p1mup5*p2mup6*F(2597)-
     -   p2mup5*p2mup6*F(2597)+p4mup5*p5mup6*F(2665)+p3mup5*p4mup6*F(
     -   2667)-p5mup6*p6mup5*F(9583)+p3mup5*p5mup6*F(9594))
       Fa22(512) = (p1mup5+p2mup5)*F(2668)+p5mup5*F(2669)+p6mup5*F(267
     -   0)+p3mup5*F(2671)+p4mup5*F(2672)
       Fa22(513) = (p1mup5+p2mup5)*F(2673)+p5mup5*F(2674)+p6mup5*F(267
     -   5)+p3mup5*F(2676)+p4mup5*F(2677)
       Fa22(514) = (p1mup6+p2mup6)*F(2678)+p5mup6*F(2679)+p6mup6*F(268
     -   0)+p3mup6*F(2681)+p4mup6*F(2682)
       Fa22(515) = (p1mup6+p2mup6)*F(2683)+p5mup6*F(2684)+p6mup6*F(268
     -   5)+p3mup6*F(2686)+p4mup6*F(2687)
       Fa22(516) = 2*(mup3mup6*F(2602)-p4mup3*(8*(p4mup6*F(1967)-p3mup
     -   6*F(2214)-(p1mup6+p2mup6)*F(2228))+4*p5mup6*F(2603)+p6mup6*F
     -   (2604))+p6mup3*p6mup6*F(2616)-p5mup3*(p6mup6*F(2607)+4*(p5mu
     -   p6*F(2605)+p4mup6*F(2609)-p3mup6*F(2618))))+8*(p3mup3*p5mup6
     -   *F(2610)+p1mup3*p5mup6*F(2613)+p2mup3*p5mup6*F(2613)+(p1mup6
     -   +p2mup6)*p5mup3*F(2619)+p5mup6*p6mup3*F(2689))+p3mup3*p6mup6
     -   *F(2690)+p1mup3*p6mup6*F(2691)+p2mup3*p6mup6*F(2691)+16*(p3m
     -   up6*p6mup3*F(9581)+p1mup6*p6mup3*F(9584)+p2mup6*p6mup3*F(958
     -   4)+p4mup6*p6mup3*F(9587)-p1mup3*p4mup6*F(9598)-p2mup3*p4mup6
     -   *F(9598)-p3mup3*p3mup6*F(9599)-p1mup3*p1mup6*F(9603)-p1mup6*
     -   p2mup3*F(9603)-p1mup3*p2mup6*F(9603)-p2mup3*p2mup6*F(9603)-p
     -   1mup6*p3mup3*F(9604)-p2mup6*p3mup3*F(9604)-p1mup3*p3mup6*F(9
     -   604)-p2mup3*p3mup6*F(9604)-p3mup3*p4mup6*F(9606))
       Fa22(517) = -16*(p1mup5*p4mup3*F(2320)+p2mup5*p4mup3*F(2320)+p1
     -   mup5*p5mup3*F(2341)+p2mup5*p5mup3*F(2341)+p1mup5*p3mup3*F(23
     -   46)+p2mup5*p3mup3*F(2346)-p1mup3*p1mup5*F(2353)-p1mup5*p2mup
     -   3*F(2353)-p1mup3*p2mup5*F(2353)-p2mup3*p2mup5*F(2353))+2*(p5
     -   mup3*p5mup5*F(2631)+p4mup3*p5mup5*F(2632))+8*(p4mup3*p6mup5*
     -   F(2633)+p5mup3*p6mup5*F(2634)+p6mup3*p6mup5*F(2636)+p4mup5*(
     -   2*(p5mup3*F(2307)+p4mup3*F(2308)+(p1mup3+p2mup3)*F(2309)-p3m
     -   up3*F(2310))+p6mup3*F(2637))+mup3mup5*F(2638)-p3mup3*p6mup5*
     -   F(2640)-p1mup3*p6mup5*F(2643)-p2mup3*p6mup5*F(2643)+p3mup5*(
     -   -2*(p5mup3*F(2276)+p3mup3*F(2279)+p4mup3*F(2319)+(p1mup3+p2m
     -   up3)*F(2346))+p6mup3*F(2645))+p1mup5*p6mup3*F(2646)+p2mup5*p
     -   6mup3*F(2646))+p5mup5*p6mup3*F(2692)+p3mup3*p5mup5*F(2693)+(
     -   p1mup3+p2mup3)*p5mup5*F(2694)
       Fa22(518) = -4*mup3mup5*(-(p6mup6*F(1911))+2*((p1mup6+p2mup6)*F
     -   (1875)+p3mup6*F(1879)+p5mup6*F(2391)+p4mup6*F(2397)))+mup5mu
     -   p6*(-4*(p6mup3*F(2376)+p4mup3*F(2381)+p5mup3*F(2382))+(p1mup
     -   3+p2mup3)*F(2695)+p3mup3*F(2696))+mup3mup6*(-2*p5mup5*F(2390
     -   )+8*((p1mup5+p2mup5)*F(1874)+p3mup5*F(1878)+p4mup5*F(2394)+p
     -   6mup5*F(9579)))
       Fa22(519) = p6mup6*F(2557)+(p1mup6+p2mup6)*F(2558)+p3mup6*F(255
     -   9)+p4mup6*F(2560)+p5mup6*F(2697)
       Fa22(520) = p5mup5*F(2653)+p6mup5*F(2698)+(p1mup5+p2mup5)*F(269
     -   9)+p3mup5*F(2700)+p4mup5*F(2701)
       Fa22(521) = mup5mup6*(p2mup3*F(2702)+p5mup3*F(2703)+p4mup3*F(27
     -   04)+2*p3mup3*F(2705))
       Fa22(522) = mup3mup6*(p1mup5*F(2707)+p2mup5*F(2709)+p5mup5*F(27
     -   11)+p6mup5*F(2712))+mup3mup5*(p6mup6*F(2713)+4*(p1mup6*F(270
     -   8)+p2mup6*F(2710)+p5mup6*F(2714)))
       Fa22(523) = mup3mup6*p4mup5*F(2715)+mup5mup6*p6mup3*F(2716)
       Fa22(524) = p4mup6*(mup3mup5*F(2718)+p2mup3*p5mup5*F(2721)+p3mu
     -   p3*p5mup5*F(2723)+p2mup3*p4mup5*F(2726)+p4mup3*(-(p5mup5*F(2
     -   719))+2*(p3mup5*F(2725)-p2mup5*F(2727))+p4mup5*F(2728))-p5mu
     -   p3*(p5mup5*F(2720)-2*(p2mup5*F(2722)+p3mup5*F(2724))+p4mup5*
     -   F(2729))+p3mup3*p4mup5*F(2732)-2*(p2mup3*p2mup5*F(2730)-p2mu
     -   p3*p3mup5*F(2731)-p2mup5*p3mup3*F(2733)+p3mup3*p3mup5*F(2734
     -   )))
       Fa22(525) = 4*p4mup6*(2*p1mup5*(-(p4mup3*F(2736))-p5mup3*F(2738
     -   )+p2mup3*F(2739)+p3mup3*F(2741))+p1mup3*(-(p5mup5*F(2737))+p
     -   4mup5*F(2742)+2*(p1mup5*F(2735)+p3mup5*F(2740)+p2mup5*F(2743
     -   ))))
       Fa22(526) = p2mup6*(p4mup3*(p4mup5*F(2746)-4*(p5mup5*F(2745)-p2
     -   mup5*F(2758)-p3mup5*F(2768)))+p4mup5*(p5mup3*F(2749)+4*(p1mu
     -   p3*F(2751)+p3mup3*F(2764)+p2mup3*F(2770))))+4*(p4mup3*(p4mup
     -   5*p5mup6*F(2747)-p1mup6*p4mup5*F(2752)+p5mup5*p5mup6*F(2756)
     -   +p3mup6*p4mup5*F(2763)-2*(p3mup6*p5mup5*F(2755)-p2mup5*p5mup
     -   6*F(2760)-p2mup5*p3mup6*F(2767))+4*(p3mup5*p3mup6*F(2769)+p3
     -   mup5*p5mup6*F(2772)))+p4mup5*(p5mup3*p5mup6*F(2748)+p1mup3*p
     -   1mup6*F(2750)+p1mup3*p5mup6*F(2753)-p1mup6*p5mup3*F(2754)+p1
     -   mup6*p2mup3*F(2757)+p2mup3*p5mup6*F(2759)-p3mup3*p5mup6*F(27
     -   61)+p3mup6*p5mup3*F(2762)+p1mup3*p3mup6*F(2765)+p1mup6*p3mup
     -   3*F(2766)+p2mup3*p3mup6*F(2771)-p3mup3*p3mup6*F(2773)))
       Fa22(527) = p1mup3*p5mup5*p5mup6*F(2781)+p2mup3*p5mup5*p5mup6*F
     -   (2785)+p2mup3*p2mup6*p5mup5*F(2786)+p1mup3*p2mup6*p5mup5*F(2
     -   789)+p5mup3*p5mup5*p5mup6*F(2795)+2*(p2mup3*p3mup6*p5mup5*F(
     -   2793)-p1mup3*p3mup6*p5mup5*F(2796))+p1mup6*(p5mup3*p5mup5*F(
     -   2782)+2*(p1mup3*p5mup5*F(2778)+p4mup3*(p5mup5*F(2775)+p2mup5
     -   *F(2784)+2*p3mup5*F(2792))+p2mup3*(p5mup5*F(2776)+2*(p2mup5*
     -   F(2788)+p3mup5*F(2794)))+p3mup3*p5mup5*F(2798))-4*(p2mup5*p5
     -   mup3*F(2780)-p3mup5*p5mup3*F(2808)-p1mup3*p2mup5*F(2810)-p1m
     -   up3*p3mup5*F(2811)-p2mup5*p3mup3*F(2817)-p3mup3*p3mup5*F(282
     -   7)))+p1mup5*(p2mup6*(p4mup3*F(2783)+4*(p2mup3*F(2787)+p5mup3
     -   *F(2790)+p1mup3*F(2809)+p3mup3*F(2816)))+4*(p1mup6*(p4mup3*F
     -   (2774)+p2mup3*F(2777)-p5mup3*F(2779)+p1mup3*F(2825)+p3mup3*F
     -   (2826))+p3mup6*(p4mup3*F(2792)+p2mup3*F(2794)+p5mup3*F(2797)
     -   +p1mup3*F(2811)+p3mup3*F(2827))+p5mup6*(p4mup3*F(2799)+p2mup
     -   3*F(2800)+p5mup3*F(2801)+p1mup3*F(2813)+p3mup3*F(2828))))+4*
     -   (p2mup5*p2mup6*p5mup3*F(2791)+p2mup3*p2mup5*p2mup6*F(2802)+p
     -   2mup3*p3mup5*p3mup6*F(2803)+p2mup3*p2mup5*p3mup6*F(2804)+p2m
     -   up3*p2mup6*p3mup5*F(2805)+p2mup3*p2mup5*p5mup6*F(2806)+p2mup
     -   3*p3mup5*p5mup6*F(2807)+p1mup3*p2mup5*p2mup6*F(2812)+p1mup3*
     -   p2mup5*p3mup6*F(2814)+p1mup3*p2mup6*p3mup5*F(2815)+p1mup3*p2
     -   mup5*p5mup6*F(2818)+p1mup3*p3mup5*p3mup6*F(2819)+p2mup5*p2mu
     -   p6*p3mup3*F(2820)+p1mup3*p3mup5*p5mup6*F(2821)+p2mup5*p3mup3
     -   *p3mup6*F(2822)+p2mup6*p3mup3*p3mup5*F(2823)+p2mup5*p3mup3*p
     -   5mup6*F(2824)+p3mup3*p3mup5*p3mup6*F(2829))
       Fa22(528) = p3mup3*p5mup5*p5mup6*F(2832)+2*(p3mup3*p3mup6*p5mup
     -   5*F(2834)+p3mup6*p5mup3*p5mup5*F(2835))+p2mup6*(-(p5mup3*p5m
     -   up5*F(2830))+p3mup3*p5mup5*F(2836)+4*p3mup5*p5mup3*F(2837))+
     -   4*(p2mup5*p5mup3*(p3mup6*F(2831)+p5mup6*F(2833))+p3mup5*p5mu
     -   p3*p5mup6*F(2838)+p3mup5*p3mup6*p5mup3*F(2839)+p3mup3*p3mup5
     -   *p5mup6*F(2840))
       Return
       End
