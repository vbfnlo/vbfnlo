c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenEvenGSim9(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4
     -   ,p3mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5m
     -   up3,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mu
     -   p3,mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa9)
       IMPLICIT NONE
      Complex*16   p1mup1, p1mup2, p1mup3, p1mup4, p1mup5, p2mup1, 
     -          p2mup2, p2mup3, p2mup4, p2mup5, p3mup1, p3mup2, 
     -          p3mup3, p3mup4, p3mup5, p4mup1, p4mup2, p4mup3, 
     -          p4mup4, p4mup5, p5mup1, p5mup2, p5mup3, p5mup4, p5mup5
       Complex*16   mup1mup2, mup1mup3, mup1mup4, mup1mup5, mup2mup3, 
     -          mup2mup4, mup2mup5, mup3mup4, mup3mup5, mup4mup5
        common/FPenEvenGSim/F
      Complex*16 Fa(123),F(3806)
      Real*8 P(1320) 
        COMMON/PPenEvenGSim/P
       Complex*16 Fa9(65:72)
       COMMON/FaPenEvenGSim/Fa
       Fa9(65) = mup1mup3*mup4mup5*F(2485)+mup1mup4*mup3mup5*F(2486)+m
     -   up1mup5*mup3mup4*F(2487)+4*(mup4mup5*p3mup3*p4mup1*F(2506)+m
     -   up3mup5*p3mup4*p4mup1*F(2506)+mup3mup4*p3mup5*p4mup1*F(2507)
     -   +mup4mup5*p3mup1*p4mup3*F(2509)+mup3mup5*p3mup1*p4mup4*F(250
     -   9)+mup3mup4*p3mup1*p4mup5*F(2510))-2*(mup1mup5*p2mup4*p3mup3
     -   *F(2203)-mup1mup4*p2mup5*p3mup3*F(2203)+mup1mup5*p2mup3*p3mu
     -   p4*F(2203)-mup1mup3*p2mup5*p3mup4*F(2203)-mup1mup4*p2mup3*p3
     -   mup5*F(2203)+mup1mup3*p2mup4*p3mup5*F(2203)-mup4mup5*p1mup1*
     -   p3mup3*F(2488)-mup3mup5*p1mup1*p3mup4*F(2488)-mup3mup4*p1mup
     -   1*p3mup5*F(2489)-mup4mup5*p1mup3*p3mup1*F(2490)-mup1mup5*(p1
     -   mup4*p3mup3+p1mup3*p3mup4)*F(2491)-mup1mup3*p1mup4*p3mup5*F(
     -   2491)-mup1mup3*p1mup5*p3mup4*F(2492)-mup1mup4*(p1mup5*p3mup3
     -   +p1mup3*p3mup5)*F(2492)-mup3mup5*p1mup4*p3mup1*F(2493)-mup3m
     -   up4*p1mup5*p3mup1*F(2494)+mup4mup5*p2mup1*p3mup3*F(2495)+mup
     -   3mup5*p2mup1*p3mup4*F(2495)-mup3mup4*p2mup1*p3mup5*F(2495)+m
     -   up4mup5*p2mup3*p3mup1*F(2496)-mup3mup5*p2mup4*p3mup1*F(2497)
     -   +mup3mup4*p2mup5*p3mup1*F(2497)-mup3mup5*p1mup4*p4mup1*F(249
     -   8)-mup4mup5*p1mup1*p4mup3*F(2498)-mup3mup5*p1mup1*p4mup4*F(2
     -   498)-mup3mup4*p1mup5*p4mup1*F(2499)-mup3mup4*p1mup1*p4mup5*F
     -   (2499)-mup4mup5*p1mup3*p4mup1*F(2500)-mup1mup5*(p1mup4*p4mup
     -   3+p1mup3*p4mup4)*F(2501)-mup1mup3*p1mup4*p4mup5*F(2501)-mup1
     -   mup3*p1mup5*p4mup4*F(2502)-mup1mup4*(p1mup5*p4mup3+p1mup3*p4
     -   mup5)*F(2502)-mup4mup5*p2mup3*p4mup1*F(2503)-mup3mup5*p2mup4
     -   *p4mup1*F(2504)-mup3mup4*p2mup5*p4mup1*F(2505)+mup4mup5*p2mu
     -   p1*p4mup3*F(2508)+mup3mup5*p2mup1*p4mup4*F(2508)-mup3mup4*p2
     -   mup1*p4mup5*F(2508)-(mup1mup4*p2mup5*p4mup3+mup1mup3*p2mup5*
     -   p4mup4-mup1mup5*(p2mup4*p4mup3+p2mup3*p4mup4)+mup1mup4*p2mup
     -   3*p4mup5-mup1mup3*p2mup4*p4mup5)*F(3789))
       Fa9(66) = 4*(mup3mup5*p3mup4*p4mup1*F(2510)+mup3mup5*p3mup1*p4m
     -   up4*F(2510))+mup1mup4*mup3mup5*F(2512)+mup4mup5*(mup1mup3*F(
     -   2511)+2*(p1mup3*p4mup1*F(2498)+p1mup1*p4mup3*F(2498)-p2mup3*
     -   p4mup1*F(2508)-p2mup1*p4mup3*F(2508)+2*(p3mup3*p4mup1*F(2509
     -   )+p3mup1*p4mup3*F(2509))+p1mup3*p3mup1*F(2513)+p1mup1*p3mup3
     -   *F(2513)-(p2mup3*p3mup1+p2mup1*p3mup3)*F(2517)))+mup3mup4*(m
     -   up1mup5*F(2511)+2*(p1mup5*p4mup1*F(2498)+p1mup1*p4mup5*F(249
     -   8)-p2mup5*p4mup1*F(2508)-p2mup1*p4mup5*F(2508)+2*(p3mup5*p4m
     -   up1*F(2509)+p3mup1*p4mup5*F(2509))+p1mup5*p3mup1*F(2513)+p1m
     -   up1*p3mup5*F(2513)-(p2mup5*p3mup1+p2mup1*p3mup5)*F(2517)))+2
     -   *(mup3mup5*p1mup4*p4mup1*F(2499)+mup3mup5*p1mup1*p4mup4*F(24
     -   99)+mup1mup5*p1mup4*p4mup3*F(2501)+mup1mup3*p1mup5*p4mup4*F(
     -   2501)+mup1mup4*p1mup3*p4mup5*F(2501)+mup1mup4*p1mup5*p4mup3*
     -   F(2502)+mup1mup5*p1mup3*p4mup4*F(2502)+mup1mup3*p1mup4*p4mup
     -   5*F(2502)+mup3mup5*p2mup4*p4mup1*F(2508)+mup3mup5*p2mup1*p4m
     -   up4*F(2508)+mup3mup5*p1mup4*p3mup1*F(2514)+mup3mup5*p1mup1*p
     -   3mup4*F(2514)+mup1mup4*p1mup5*p3mup3*F(2515)+mup1mup5*p1mup3
     -   *p3mup4*F(2515)+mup1mup3*p1mup4*p3mup5*F(2515)+mup1mup5*p1mu
     -   p4*p3mup3*F(2516)+mup1mup3*p1mup5*p3mup4*F(2516)+mup1mup4*p1
     -   mup3*p3mup5*F(2516)+mup3mup5*p2mup4*p3mup1*F(2517)+mup3mup5*
     -   p2mup1*p3mup4*F(2517)-mup1mup5*p2mup4*p3mup3*F(3788)+mup1mup
     -   4*p2mup5*p3mup3*F(3788)+mup1mup5*p2mup3*p3mup4*F(3788)-mup1m
     -   up3*p2mup5*p3mup4*F(3788)-mup1mup4*p2mup3*p3mup5*F(3788)+mup
     -   1mup3*p2mup4*p3mup5*F(3788)-(mup1mup5*p2mup4*p4mup3-mup1mup4
     -   *p2mup5*p4mup3-mup1mup5*p2mup3*p4mup4+mup1mup3*p2mup5*p4mup4
     -   +mup1mup4*p2mup3*p4mup5-mup1mup3*p2mup4*p4mup5)*F(3789))
       Fa9(67) = mup1mup3*mup4mup5*F(2518)+mup1mup4*(mup3mup5*F(2519)-
     -   2*((p2mup5*p3mup3+p2mup3*p3mup5)*F(2203)+(p2mup5*p4mup3+p2mu
     -   p3*p4mup5)*F(2212)-(p1mup5*p3mup3+p1mup3*p3mup5)*F(2491)-(p1
     -   mup5*p4mup3+p1mup3*p4mup5)*F(2525)))+mup1mup5*(mup3mup4*F(25
     -   20)+2*(p2mup4*p3mup3*F(2203)-p2mup3*p3mup4*F(2203)+p2mup4*p4
     -   mup3*F(2212)-p2mup3*p4mup4*F(2212)+p1mup3*p3mup4*F(2491)+p1m
     -   up4*p3mup3*F(2492)+p1mup3*p4mup4*F(2525)+p1mup4*p4mup3*F(252
     -   7)))+2*(mup1mup3*p2mup5*p3mup4*F(2203)+mup1mup3*p2mup4*p3mup
     -   5*F(2203)+mup1mup3*p2mup5*p4mup4*F(2212)+mup1mup3*p2mup4*p4m
     -   up5*F(2212)+mup3mup5*p1mup1*p3mup4*F(2488)+mup3mup4*p1mup1*p
     -   3mup5*F(2488)+mup4mup5*p1mup1*p3mup3*F(2489)+mup3mup5*p1mup4
     -   *p3mup1*F(2490)+mup1mup3*p1mup5*p3mup4*F(2492)+mup1mup3*p1mu
     -   p4*p3mup5*F(2492)+mup3mup4*p1mup5*p3mup1*F(2493)+mup4mup5*p2
     -   mup1*p3mup3*F(2495)-mup3mup5*p2mup1*p3mup4*F(2495)-mup3mup4*
     -   p2mup1*p3mup5*F(2495)+mup4mup5*p2mup3*p3mup1*F(2496)-mup3mup
     -   5*p2mup4*p3mup1*F(2496)+mup3mup4*p2mup5*p3mup1*F(2497)+mup4m
     -   up5*p1mup3*p3mup1*F(2521)+mup4mup5*p1mup1*p4mup3*F(2522)+mup
     -   3mup5*p1mup1*p4mup4*F(2523)+mup3mup4*p1mup1*p4mup5*F(2523)+m
     -   up4mup5*p1mup3*p4mup1*F(2524)+mup3mup5*p1mup4*p4mup1*F(2526)
     -   +mup1mup3*p1mup5*p4mup4*F(2527)+mup1mup3*p1mup4*p4mup5*F(252
     -   7)+mup3mup4*p1mup5*p4mup1*F(2528)+mup4mup5*p2mup1*p4mup3*F(2
     -   529)-mup3mup5*p2mup1*p4mup4*F(2529)-mup3mup4*p2mup1*p4mup5*F
     -   (2529)+mup4mup5*p2mup3*p4mup1*F(2530)-mup3mup5*p2mup4*p4mup1
     -   *F(2530)+mup3mup4*p2mup5*p4mup1*F(2531))+4*(mup3mup5*p3mup4*
     -   p4mup1*F(2506)+mup3mup4*p3mup5*p4mup1*F(2506)+mup4mup5*p3mup
     -   3*p4mup1*F(2507)+mup4mup5*p3mup1*p4mup3*F(2532)+p3mup1*(mup3
     -   mup5*p4mup4+mup3mup4*p4mup5)*F(2533))
       Fa9(68) = -(p1mup3*p1mup4*p2mup1*p4mup5*F(2550))-p1mup3*p1mup5*
     -   p2mup1*p4mup4*F(2551)-p1mup3*p2mup1*p2mup4*p4mup5*F(2552)-p1
     -   mup3*p2mup1*p2mup5*p4mup4*F(2553)-p1mup4*p2mup1*p2mup3*p4mup
     -   5*F(2554)-p1mup5*p2mup1*p2mup3*p4mup4*F(2555)-p1mup3*p2mup1*
     -   p3mup4*p4mup5*F(2570)-p1mup3*p2mup1*p3mup5*p4mup4*F(2571)-p1
     -   mup4*p2mup1*p3mup3*p4mup5*F(2572)+p1mup4*p2mup1*p3mup5*p4mup
     -   3*F(2573)-p1mup5*p2mup1*p3mup3*p4mup4*F(2574)-p1mup5*p2mup1*
     -   p3mup4*p4mup3*F(2575)-p2mup1*p2mup3*p3mup4*p4mup5*F(2576)-p2
     -   mup1*p2mup3*p3mup5*p4mup4*F(2577)-p2mup1*p2mup4*p3mup3*p4mup
     -   5*F(2578)+p2mup1*p2mup4*p3mup5*p4mup3*F(2579)-p2mup1*p2mup5*
     -   p3mup3*p4mup4*F(2580)-p2mup1*p2mup5*p3mup4*p4mup3*F(2581)+p1
     -   mup3*p1mup4*p2mup5*p4mup1*F(2582)-p1mup3*p1mup5*p2mup4*p4mup
     -   1*F(2583)-p1mup4*p2mup1*p2mup5*p4mup3*F(2584)+p1mup4*p2mup3*
     -   p2mup5*p4mup1*F(2585)+p1mup5*p2mup1*p2mup4*p4mup3*F(2586)-p1
     -   mup5*p2mup3*p2mup4*p4mup1*F(2587)+p1mup3*p1mup4*p3mup1*p4mup
     -   5*F(2593)+p1mup3*p1mup4*p3mup5*p4mup1*F(2594)+p1mup3*p1mup5*
     -   p3mup1*p4mup4*F(2595)-p1mup3*p1mup5*p3mup4*p4mup1*F(2596)+p1
     -   mup3*p2mup4*p3mup1*p4mup5*F(2597)+p1mup3*p2mup5*p3mup1*p4mup
     -   4*F(2598)+p1mup4*p2mup3*p3mup5*p4mup1*F(2599)-p1mup4*p2mup5*
     -   p3mup1*p4mup3*F(2600)+p1mup4*p2mup5*p3mup3*p4mup1*F(2601)-p1
     -   mup5*p2mup3*p3mup4*p4mup1*F(2602)+p1mup5*p2mup4*p3mup1*p4mup
     -   3*F(2603)-p1mup5*p2mup4*p3mup3*p4mup1*F(2604)-p1mup3*p2mup4*
     -   p3mup5*p4mup1*F(2610)+p1mup3*p2mup5*p3mup4*p4mup1*F(2611)-p1
     -   mup4*p2mup3*p3mup1*p4mup5*F(2612)-p1mup5*p2mup3*p3mup1*p4mup
     -   4*F(2613)-p2mup3*p2mup4*p3mup1*p4mup5*F(2614)-p2mup3*p2mup4*
     -   p3mup5*p4mup1*F(2615)-p2mup3*p2mup5*p3mup1*p4mup4*F(2616)+p2
     -   mup3*p2mup5*p3mup4*p4mup1*F(2617)+p1mup3*p2mup4*p2mup5*p4mup
     -   1*F(2620)-p1mup4*p1mup5*p2mup1*p4mup3*F(2621)+p1mup4*p1mup5*
     -   p2mup3*p4mup1*F(2622)-p1mup4*p1mup5*p3mup1*p4mup3*F(2623)+p1
     -   mup4*p1mup5*p3mup3*p4mup1*F(2624)+p1mup1*(p1mup5*p2mup3*p4mu
     -   p4*F(2538)+p1mup5*p2mup4*p4mup3*F(2539)+p1mup3*(p2mup4*p4mup
     -   5*F(2534)+p2mup5*p4mup4*F(2535)+2*p3mup4*p3mup5*F(2540)+p3mu
     -   p4*p4mup5*F(2542)+p3mup5*p4mup4*F(2543))+p1mup5*p3mup3*p4mup
     -   4*F(2546)+p1mup5*p3mup4*p4mup3*F(2547)+p2mup3*p2mup4*p4mup5*
     -   F(2548)-p2mup3*p2mup5*p4mup4*F(2549)-p2mup3*p3mup4*p4mup5*F(
     -   2564)-p2mup3*p3mup5*p4mup4*F(2565)-p2mup4*p3mup3*p4mup5*F(25
     -   66)+p2mup4*p3mup5*p4mup3*F(2567)-p2mup5*p3mup3*p4mup4*F(2568
     -   )-p2mup5*p3mup4*p4mup3*F(2569)+2*(p1mup5*p3mup3*p3mup4*F(254
     -   1)-p2mup3*p3mup4*p3mup5*F(2558)-p2mup5*p3mup3*p3mup4*F(2559)
     -   +p2mup4*p3mup3*p3mup5*F(2589))+p2mup4*p2mup5*p4mup3*F(2618)+
     -   p1mup4*(p2mup3*p4mup5*F(2536)+p3mup3*p4mup5*F(2544)+2*p3mup3
     -   *p3mup5*F(2588)-p4mup3*(p2mup5*F(2537)+p3mup5*F(2545)-4*p1mu
     -   p5*F(2630))))-2*(p2mup1*p2mup3*p2mup4*p4mup5*F(2556)+p2mup1*
     -   p2mup3*p2mup5*p4mup4*F(2557)+p1mup3*p2mup1*p3mup4*p3mup5*F(2
     -   560)+p1mup5*p2mup1*p3mup3*p3mup4*F(2561)+p2mup1*p2mup3*p3mup
     -   4*p3mup5*F(2562)+p2mup1*p2mup5*p3mup3*p3mup4*F(2563)-p1mup3*
     -   p1mup5*p3mup1*p3mup4*F(2588)-p1mup3*p1mup4*p3mup1*p3mup5*F(2
     -   590)-p1mup3*p2mup4*p3mup1*p3mup5*F(2591)-p1mup3*p2mup5*p3mup
     -   1*p3mup4*F(2592)+p1mup4*p2mup1*p3mup3*p3mup5*F(2605)+p1mup4*
     -   p2mup3*p3mup1*p3mup5*F(2606)+p1mup5*p2mup3*p3mup1*p3mup4*F(2
     -   607)+p2mup3*p2mup5*p3mup1*p3mup4*F(2608)+p2mup1*p2mup4*p3mup
     -   3*p3mup5*F(2608)+p2mup3*p2mup4*p3mup1*p3mup5*F(2609)-p1mup3*
     -   p1mup4*p1mup5*p4mup1*F(2619)+p2mup3*p2mup4*p2mup5*p4mup1*F(2
     -   625)-p1mup4*p1mup5*p3mup1*p3mup3*F(2626)-p1mup4*p2mup5*p3mup
     -   1*p3mup3*F(2627)-p1mup5*p2mup4*p3mup1*p3mup3*F(2628)-p2mup4*
     -   p2mup5*p3mup1*p3mup3*F(2629)+p2mup1*p2mup4*p2mup5*p4mup3*F(2
     -   631))
       Fa9(69) = -(p1mup4*p2mup1*p4mup3*p4mup5*F(2642))+p1mup4*p2mup3*
     -   p4mup1*p4mup5*F(2643)+p2mup4*p2mup5*p3mup1*p4mup3*F(2644)-p2
     -   mup4*p2mup5*p3mup3*p4mup1*F(2645)+4*(p1mup4*p3mup1*p3mup3*p3
     -   mup5*F(2648)+p2mup4*p3mup1*p3mup3*p3mup5*F(2650))+p1mup3*(p2
     -   mup1*p4mup4*p4mup5*F(2640)+p2mup4*p4mup1*p4mup5*F(2641)+2*(p
     -   3mup1*p3mup4*p3mup5*F(2633)+p1mup4*p4mup1*p4mup5*F(2639)+p1m
     -   up5*p4mup1*p4mup4*F(2652))+p2mup5*p4mup1*p4mup4*F(2654))+p1m
     -   up4*p2mup5*p4mup1*p4mup3*F(2656)-p1mup5*p2mup1*p4mup3*p4mup4
     -   *F(2657)+p1mup5*p2mup3*p4mup1*p4mup4*F(2658)+p1mup5*p2mup4*p
     -   4mup1*p4mup3*F(2659)+p1mup1*(2*(p3mup3*p3mup4*p3mup5*F(2632)
     -   +p1mup3*p4mup4*p4mup5*F(2636))+p2mup3*p4mup4*p4mup5*F(2637)+
     -   p4mup3*(p2mup4*p4mup5*F(2638)+p2mup5*p4mup4*F(2653)+2*(p1mup
     -   4*p4mup5*F(2652)+p1mup5*p4mup4*F(2662))))-2*(p2mup1*p3mup3*p
     -   3mup4*p3mup5*F(2634)+p2mup3*p3mup1*p3mup4*p3mup5*F(2635)+p2m
     -   up1*p2mup3*p4mup4*p4mup5*F(2646)+p2mup3*p2mup4*p4mup1*p4mup5
     -   *F(2647)-p1mup5*p3mup1*p3mup3*p3mup4*F(2649)-p2mup5*p3mup1*p
     -   3mup3*p3mup4*F(2651)-p1mup4*p1mup5*p4mup1*p4mup3*F(2655)-p2m
     -   up4*p2mup5*p4mup1*p4mup3*F(2660)+p2mup3*p2mup5*p4mup1*p4mup4
     -   *F(2661)+p2mup1*p2mup4*p4mup3*p4mup5*F(2661)+p2mup1*p2mup5*p
     -   4mup3*p4mup4*F(2663))
       Fa9(70) = p1mup4*p2mup2*p3mup5*p4mup1*F(2686)+p1mup4*p2mup5*p3m
     -   up2*p4mup1*F(2688)+p1mup5*p2mup1*p3mup4*p4mup2*F(2691)+p1mup
     -   5*p2mup2*p3mup4*p4mup1*F(2692)+p1mup5*p2mup4*p3mup1*p4mup2*F
     -   (2693)+p1mup5*p3mup2*p3mup4*p4mup1*F(2696)+2*(p1mup4*p2mup1*
     -   p3mup5*p4mup2*F(2685)+p1mup4*p2mup5*p3mup1*p4mup2*F(2687)+p1
     -   mup5*p2mup4*p3mup2*p4mup1*F(2694)-p1mup4*p2mup1*p3mup2*p4mup
     -   5*F(2701)+p1mup4*p2mup2*p3mup1*p4mup5*F(2702))+p1mup5*p2mup1
     -   *p3mup2*p4mup4*F(2722)+p1mup5*p2mup2*p3mup1*p4mup4*F(2723)+p
     -   1mup5*p3mup1*p3mup2*p4mup4*F(2724)+p1mup1*(16*p1mup2*p1mup4*
     -   p1mup5*F(2664)+p2mup5*p3mup2*p4mup4*F(2676)+p3mup2*p3mup4*p4
     -   mup5*F(2677)+p3mup2*p3mup5*p4mup4*F(2678)+p2mup2*(p3mup4*p4m
     -   up5*F(2673)+p3mup5*p4mup4*F(2674)+4*(p2mup4*p4mup5*F(2665)+p
     -   2mup5*p4mup4*F(2666)+p4mup4*p4mup5*F(2703)))+p2mup4*(-2*p3mu
     -   p2*p4mup5*F(2675)+p4mup2*(p3mup5*F(2698)+4*(p2mup5*F(2697)+p
     -   4mup5*F(2704))))+p2mup5*p3mup4*p4mup2*F(2709)+p3mup4*p4mup2*
     -   p4mup5*F(2714)-4*(p3mup4*p3mup5*p4mup2*F(2710)+p3mup2*p4mup4
     -   *p4mup5*F(2713)-p2mup5*p4mup2*p4mup4*F(2725))+p3mup5*p4mup2*
     -   p4mup4*F(2726))+p1mup2*(p2mup5*p3mup1*p4mup4*F(2682)+p2mup1*
     -   (p3mup4*p4mup5*F(2679)+p3mup5*p4mup4*F(2680)+4*(p2mup4*p4mup
     -   5*F(2667)+p2mup5*p4mup4*F(2668)+p4mup4*p4mup5*F(2705)))-2*p2
     -   mup4*(p3mup1*p4mup5*F(2681)+p3mup5*p4mup1*F(2700)-2*(p2mup5*
     -   p4mup1*F(2699)+p4mup1*p4mup5*F(2706)))+p2mup5*p3mup4*p4mup1*
     -   F(2711)-4*(p3mup1*p3mup4*p4mup5*F(2683)+p3mup1*p3mup5*p4mup4
     -   *F(2684)+p3mup4*p3mup5*p4mup1*F(2712)+p3mup1*p4mup4*p4mup5*F
     -   (2715)+p3mup4*p4mup1*p4mup5*F(2716)-p2mup5*p4mup1*p4mup4*F(2
     -   727)+p3mup5*p4mup1*p4mup4*F(2728)))-4*(p1mup4*p2mup1*p2mup5*
     -   p4mup2*F(2669)-p1mup4*p2mup2*p2mup5*p4mup1*F(2670)-p1mup5*p2
     -   mup1*p2mup4*p4mup2*F(2671)+p1mup5*p2mup2*p2mup4*p4mup1*F(267
     -   2)-p1mup4*p3mup1*p3mup5*p4mup2*F(2689)+p1mup4*p3mup2*p3mup5*
     -   p4mup1*F(2690)+p1mup5*p3mup1*p3mup4*p4mup2*F(2695)+p1mup4*p2
     -   mup1*p4mup2*p4mup5*F(2707)-p1mup4*p2mup2*p4mup1*p4mup5*F(270
     -   8)-p1mup4*p3mup1*p4mup2*p4mup5*F(2717)+p1mup4*p3mup2*p4mup1*
     -   p4mup5*F(2718)-p1mup4*p2mup1*p2mup2*p4mup5*F(2719)-p1mup5*p2
     -   mup1*p2mup2*p4mup4*F(2720)+p1mup4*p3mup1*p3mup2*p4mup5*F(272
     -   1)-p1mup4*p2mup5*p4mup1*p4mup2*F(2729))
       Fa9(71) = p1mup1*p2mup5*p3mup4*p4mup2*F(1312)-p1mup1*p2mup4*p3m
     -   up5*p4mup2*F(1312)-p1mup1*p2mup2*p3mup5*p4mup4*F(1312)+p1mup
     -   1*p2mup2*p3mup4*p4mup5*F(1312)+p1mup1*p2mup5*p3mup2*p4mup4*F
     -   (1314)-p1mup1*p2mup4*p3mup2*p4mup5*F(1314)-p1mup2*p2mup1*p3m
     -   up5*p4mup4*F(1315)+p1mup2*p2mup1*p3mup4*p4mup5*F(1315)+p1mup
     -   2*p2mup5*p3mup4*p4mup1*F(1317)-p1mup2*p2mup4*p3mup5*p4mup1*F
     -   (1317)-2*(p1mup1*p3mup2*p3mup5*p4mup4*F(1324)-p1mup1*p3mup2*
     -   p3mup4*p4mup5*F(1324)+p1mup1*p3mup5*p4mup2*p4mup4*F(1458)-p1
     -   mup1*p3mup4*p4mup2*p4mup5*F(1458))-4*(p1mup2*p3mup1*p3mup5*p
     -   4mup4*F(1456)-p1mup2*p3mup1*p3mup4*p4mup5*F(1456)+p1mup2*p3m
     -   up5*p4mup1*p4mup4*F(1459)-p1mup2*p3mup4*p4mup1*p4mup5*F(1459
     -   ))+p1mup2*p3mup1*(p2mup5*p4mup4-p2mup4*p4mup5)*F(3791)+p1mup
     -   5*(p2mup2*p3mup4*p4mup1*F(1312)-p2mup4*p3mup2*p4mup1*F(1314)
     -   +p2mup4*p3mup1*p4mup2*F(1320)+p2mup2*p3mup1*p4mup4*F(1320)-p
     -   2mup1*p3mup4*p4mup2*F(1322)+4*p3mup1*p3mup4*p4mup2*F(1456)-p
     -   2mup1*p3mup2*p4mup4*F(3791)+2*(p3mup2*p3mup4*p4mup1*F(1324)+
     -   p3mup1*p3mup2*p4mup4*F(3792)))-p1mup4*(p2mup2*p3mup5*p4mup1*
     -   F(1312)-p2mup5*p3mup2*p4mup1*F(1314)+p2mup5*p3mup1*p4mup2*F(
     -   1320)+p2mup2*p3mup1*p4mup5*F(1320)-p2mup1*p3mup5*p4mup2*F(13
     -   22)+4*p3mup1*p3mup5*p4mup2*F(1456)-p2mup1*p3mup2*p4mup5*F(37
     -   91)+2*(p3mup2*p3mup5*p4mup1*F(1324)+p3mup1*p3mup2*p4mup5*F(3
     -   792)+p3mup1*p4mup2*p4mup5*F(3794)))
       Fa9(72) = p1mup5*p2mup1*p3mup2*p4mup4*F(1333)+p1mup4*p2mup1*p3m
     -   up2*p4mup5*F(1333)+p1mup4*p2mup2*p3mup5*p4mup1*F(1336)-p1mup
     -   5*p2mup2*p3mup1*p4mup4*F(1337)-p1mup4*p2mup2*p3mup1*p4mup5*F
     -   (1337)-p1mup4*p2mup5*p3mup2*p4mup1*F(1338)-p1mup5*p2mup4*p3m
     -   up1*p4mup2*F(1339)+p1mup4*p2mup1*p3mup5*p4mup2*F(1341)+p1mup
     -   4*p2mup5*p3mup1*p4mup2*F(1342)+p1mup5*p2mup1*p3mup4*p4mup2*F
     -   (1343)+p1mup5*p2mup2*p3mup4*p4mup1*F(1344)+p1mup5*p2mup4*p3m
     -   up2*p4mup1*F(1346)+p1mup1*(p2mup2*(p3mup5*p4mup4+p3mup4*p4mu
     -   p5)*F(1327)+p2mup4*p3mup5*p4mup2*F(1328)-p2mup5*p3mup2*p4mup
     -   4*F(1329)-p2mup4*p3mup2*p4mup5*F(1329)+p2mup5*p3mup4*p4mup2*
     -   F(1330)+4*p3mup4*p3mup5*p4mup2*F(1450)+2*(p3mup2*p3mup5*p4mu
     -   p4*F(1347)+p3mup2*p3mup4*p4mup5*F(1347)+p4mup2*(p3mup5*p4mup
     -   4+p3mup4*p4mup5)*F(1466)))+p1mup2*(p2mup1*(p3mup5*p4mup4+p3m
     -   up4*p4mup5)*F(1331)+p2mup4*p3mup5*p4mup1*F(1332)-p2mup5*p3mu
     -   p1*p4mup4*F(1333)-p2mup4*p3mup1*p4mup5*F(1333)+p2mup5*p3mup4
     -   *p4mup1*F(1334)+8*p3mup4*p3mup5*p4mup1*F(1451)+4*(p3mup1*p3m
     -   up5*p4mup4*F(1462)+p3mup1*p3mup4*p4mup5*F(1462)+p4mup1*(p3mu
     -   p5*p4mup4+p3mup4*p4mup5)*F(1467)))-2*(-(p1mup4*p3mup2*p3mup5
     -   *p4mup1*F(1463))+p1mup5*p3mup1*p3mup2*p4mup4*F(1464)+p1mup4*
     -   p3mup1*p3mup2*p4mup5*F(1464)-p1mup5*p3mup2*p3mup4*p4mup1*F(1
     -   465)+p1mup4*p3mup1*p4mup2*p4mup5*F(1469))+4*(p1mup4*p3mup1*p
     -   3mup5*p4mup2*F(1456)+p1mup5*p3mup1*p3mup4*p4mup2*F(2730))
       Return
       End
