c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenEvenG8(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup
     -   1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3
     -   mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3
     -   ,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,
     -   mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa8)
       IMPLICIT NONE
      Complex*16   p1mup1, p1mup2, p1mup3, p1mup4, p1mup5, p2mup1, 
     -          p2mup2, p2mup3, p2mup4, p2mup5, p3mup1, p3mup2, 
     -          p3mup3, p3mup4, p3mup5, p4mup1, p4mup2, p4mup3, 
     -          p4mup4, p4mup5, p5mup1, p5mup2, p5mup3, p5mup4, p5mup5
       Complex*16   mup1mup2, mup1mup3, mup1mup4, mup1mup5, mup2mup3, 
     -          mup2mup4, mup2mup5, mup3mup4, mup3mup5, mup4mup5
        common/FPenEvenG/F
      Complex*16 Fa(150),F(4300)
      Real*8 P(2563) 
        COMMON/PPenEvenG/P
       Complex*16 Fa8(71:80)
       COMMON/FaPenEvenG/Fa
       Fa8(71) = p1mup1*p2mup5*p3mup4*p4mup2*F(1312)-p1mup1*p2mup4*p3m
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
     -   ))+p1mup2*p3mup1*(p2mup5*p4mup4-p2mup4*p4mup5)*F(4285)+p1mup
     -   5*(p2mup2*p3mup4*p4mup1*F(1312)-p2mup4*p3mup2*p4mup1*F(1314)
     -   +p2mup4*p3mup1*p4mup2*F(1320)+p2mup2*p3mup1*p4mup4*F(1320)-p
     -   2mup1*p3mup4*p4mup2*F(1322)+4*p3mup1*p3mup4*p4mup2*F(1456)-p
     -   2mup1*p3mup2*p4mup4*F(4285)+2*(p3mup2*p3mup4*p4mup1*F(1324)+
     -   p3mup1*p3mup2*p4mup4*F(4286)))-p1mup4*(p2mup2*p3mup5*p4mup1*
     -   F(1312)-p2mup5*p3mup2*p4mup1*F(1314)+p2mup5*p3mup1*p4mup2*F(
     -   1320)+p2mup2*p3mup1*p4mup5*F(1320)-p2mup1*p3mup5*p4mup2*F(13
     -   22)+4*p3mup1*p3mup5*p4mup2*F(1456)-p2mup1*p3mup2*p4mup5*F(42
     -   85)+2*(p3mup2*p3mup5*p4mup1*F(1324)+p3mup1*p3mup2*p4mup5*F(4
     -   286)+p3mup1*p4mup2*p4mup5*F(4288)))
       Fa8(72) = p1mup5*p2mup1*p3mup2*p4mup4*F(1333)+p1mup4*p2mup1*p3m
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
       Fa8(73) = mup2mup5*p1mup1*p1mup4*F(2734)+mup2mup5*p2mup1*p3mup4
     -   *F(2738)+mup2mup4*p2mup1*p3mup5*F(2739)+mup1mup4*p2mup2*p3mu
     -   p5*F(2743)+mup2mup5*p2mup4*p3mup1*F(2745)+mup1mup5*(-2*(p1mu
     -   p2*p1mup4*F(2732)-p2mup2*p2mup4*F(2740))+p2mup2*p3mup4*F(274
     -   2)+p2mup4*p3mup2*F(2746))+mup1mup2*p2mup4*p3mup5*F(2747)+mup
     -   2mup4*p2mup5*p3mup1*F(2748)+mup1mup4*p2mup5*p3mup2*F(2749)+m
     -   up1mup2*p2mup5*p3mup4*F(2750)+mup4mup5*(p1mup1*p1mup2*F(2731
     -   )+p2mup1*p3mup2*F(2737)+p2mup2*p3mup1*F(2741)+2*(p2mup1*p2mu
     -   p2*F(2751)+p3mup1*p3mup2*F(2753)))+2*(mup1mup4*p1mup2*p1mup5
     -   *F(2733)-mup1mup2*p1mup4*p1mup5*F(2735)+mup2mup5*p2mup1*p2mu
     -   p4*F(2736)+mup1mup2*p2mup4*p2mup5*F(2744)+mup2mup4*p2mup1*p2
     -   mup5*F(2752)+mup2mup4*p1mup1*p1mup5*F(2754)+mup1mup4*p2mup2*
     -   p2mup5*F(2755))+p2mup1*p2mup2*p3mup4*p4mup5*F(2758)+p2mup1*p
     -   2mup2*p3mup5*p4mup4*F(2759)+p2mup1*p2mup5*p3mup2*p4mup4*F(27
     -   61)+p2mup2*p3mup1*p3mup4*p4mup5*F(2766)+p2mup2*p3mup1*p3mup5
     -   *p4mup4*F(2767)+p2mup5*p3mup1*p3mup2*p4mup4*F(2769)+p2mup2*p
     -   3mup4*p4mup1*p4mup5*F(2787)+p2mup2*p3mup5*p4mup1*p4mup4*F(27
     -   88)+p1mup5*p3mup1*p4mup2*p4mup4*F(2795)+p1mup5*p3mup4*p4mup1
     -   *p4mup2*F(2797)+p2mup5*p3mup4*p4mup1*p4mup2*F(2798)-4*(p2mup
     -   1*p2mup4*p3mup2*p4mup5*F(2760)+p2mup4*p3mup1*p3mup2*p4mup5*F
     -   (2768)+p1mup4*p3mup5*p4mup1*p4mup2*F(2772)+p1mup5*p2mup1*p4m
     -   up2*p4mup4*F(2773)-p1mup5*p2mup2*p4mup1*p4mup4*F(2774)-p1mup
     -   5*p2mup4*p4mup1*p4mup2*F(2775)+p2mup4*p3mup1*p3mup5*p4mup2*F
     -   (2783)-p2mup4*p3mup2*p4mup1*p4mup5*F(2791)+p2mup4*p3mup5*p4m
     -   up1*p4mup2*F(2792)+p1mup5*p3mup2*p4mup1*p4mup4*F(2796)+p2mup
     -   1*p2mup4*p3mup5*p4mup2*F(2804))+p2mup1*p2mup5*p3mup4*p4mup2*
     -   F(2809)+p2mup5*p3mup1*p3mup4*p4mup2*F(2811)+p2mup5*p3mup2*p4
     -   mup1*p4mup4*F(2817)+8*(p2mup1*p2mup2*p2mup4*p4mup5*F(2756)+p
     -   2mup1*p2mup2*p2mup5*p4mup4*F(2757)-p2mup1*p3mup2*p3mup4*p4mu
     -   p5*F(2762)-p2mup1*p3mup2*p3mup5*p4mup4*F(2763)-p2mup2*p2mup4
     -   *p3mup1*p4mup5*F(2764)-p2mup2*p2mup5*p3mup1*p4mup4*F(2765)-p
     -   3mup1*p3mup2*p3mup4*p4mup5*F(2770)-p3mup1*p3mup2*p3mup5*p4mu
     -   p4*F(2771)+p1mup1*p4mup2*p4mup4*p4mup5*F(2776)+p1mup4*p4mup1
     -   *p4mup2*p4mup5*F(2778)+p2mup2*p2mup4*p3mup5*p4mup1*F(2779)+p
     -   2mup2*p3mup4*p3mup5*p4mup1*F(2780)+p2mup4*p2mup5*p3mup1*p4mu
     -   p2*F(2781)+p2mup4*p2mup5*p3mup2*p4mup1*F(2782)+p2mup4*p3mup2
     -   *p3mup5*p4mup1*F(2784)+p2mup2*p2mup4*p4mup1*p4mup5*F(2785)+p
     -   2mup2*p3mup1*p4mup4*p4mup5*F(2786)+p2mup4*p2mup5*p4mup1*p4mu
     -   p2*F(2789)+p2mup4*p3mup1*p4mup2*p4mup5*F(2790)+p2mup2*p4mup1
     -   *p4mup4*p4mup5*F(2793)+p2mup4*p4mup1*p4mup2*p4mup5*F(2794)-p
     -   3mup2*p3mup4*p4mup1*p4mup5*F(2799)-p3mup4*p3mup5*p4mup1*p4mu
     -   p2*F(2800)-p3mup2*p4mup1*p4mup4*p4mup5*F(2801)-p3mup4*p4mup1
     -   *p4mup2*p4mup5*F(2802)+p2mup1*p2mup4*p2mup5*p4mup2*F(2803)+p
     -   2mup2*p2mup5*p3mup4*p4mup1*F(2806)+p2mup1*p2mup2*p4mup4*p4mu
     -   p5*F(2807)+p2mup2*p2mup5*p4mup1*p4mup4*F(2808)+p2mup1*p2mup4
     -   *p4mup2*p4mup5*F(2808)-p2mup1*p3mup4*p3mup5*p4mup2*F(2810)+p
     -   2mup5*p3mup2*p3mup4*p4mup1*F(2812)-p3mup1*p3mup4*p3mup5*p4mu
     -   p2*F(2813)-p2mup1*p3mup2*p4mup4*p4mup5*F(2815)-p2mup1*p3mup4
     -   *p4mup2*p4mup5*F(2816)-p3mup1*p3mup2*p4mup4*p4mup5*F(2818)-p
     -   3mup2*p3mup5*p4mup1*p4mup4*F(2819)-p3mup1*p3mup4*p4mup2*p4mu
     -   p5*F(2819)+p1mup5*p4mup1*p4mup2*p4mup4*F(2820)+p2mup5*p3mup1
     -   *p4mup2*p4mup4*F(2821)-p3mup1*p3mup5*p4mup2*p4mup4*F(2822)+p
     -   2mup5*p4mup1*p4mup2*p4mup4*F(2823)-p3mup5*p4mup1*p4mup2*p4mu
     -   p4*F(2825)+p2mup1*p2mup5*p4mup2*p4mup4*F(2826)-p2mup1*p3mup5
     -   *p4mup2*p4mup4*F(2827))+16*(p1mup2*p4mup1*p4mup4*p4mup5*F(27
     -   77)+p2mup2*p2mup4*p2mup5*p4mup1*F(2805)-p3mup2*p3mup4*p3mup5
     -   *p4mup1*F(2814)+p3mup1*p4mup2*p4mup4*p4mup5*F(2824)+p2mup1*p
     -   4mup2*p4mup4*p4mup5*F(2828))
       Fa8(74) = mup4mup5*p2mup1*p3mup2*F(2832)-mup2mup5*p2mup1*p3mup4
     -   *F(2833)+mup2mup4*p2mup1*p3mup5*F(2834)+mup4mup5*p2mup2*p3mu
     -   p1*F(2835)+mup1mup5*p2mup2*p3mup4*F(2836)+mup1mup4*p2mup2*p3
     -   mup5*F(2837)+mup2mup5*p2mup4*p3mup1*F(2838)+mup1mup5*p2mup4*
     -   p3mup2*F(2839)+mup1mup2*p2mup4*p3mup5*F(2840)+mup2mup4*p2mup
     -   5*p3mup1*F(2841)+mup1mup4*p2mup5*p3mup2*F(2842)+mup1mup2*p2m
     -   up5*p3mup4*F(2843)-4*(p2mup5*p3mup1*p3mup4*p4mup2*F(1582)-p2
     -   mup4*p3mup1*p3mup5*p4mup2*F(1582)-p2mup2*p3mup1*p3mup5*p4mup
     -   4*F(1582)+p2mup2*p3mup1*p3mup4*p4mup5*F(1582)+p2mup1*p2mup5*
     -   p3mup4*p4mup2*F(1585)-p2mup1*p2mup4*p3mup5*p4mup2*F(1585)-p2
     -   mup1*p2mup2*p3mup5*p4mup4*F(1585)+p2mup1*p2mup2*p3mup4*p4mup
     -   5*F(1585)+p2mup5*p3mup4*p4mup1*p4mup2*F(1701)-p2mup4*p3mup5*
     -   p4mup1*p4mup2*F(1701)-p2mup2*p3mup5*p4mup1*p4mup4*F(1701)+p2
     -   mup2*p3mup4*p4mup1*p4mup5*F(1701)+p2mup5*p3mup1*p3mup2*p4mup
     -   4*F(4286)-p2mup4*p3mup1*p3mup2*p4mup5*F(4286)+p2mup5*p3mup2*
     -   p4mup1*p4mup4*F(4288)-p2mup4*p3mup2*p4mup1*p4mup5*F(4288)+p2
     -   mup1*p2mup5*p3mup2*p4mup4*F(4290)-p2mup1*p2mup4*p3mup2*p4mup
     -   5*F(4290))-8*(p2mup2*p2mup5*p3mup4*p4mup1*F(1587)-p2mup2*p2m
     -   up4*p3mup5*p4mup1*F(1587)+p2mup5*p3mup2*p3mup4*p4mup1*F(2848
     -   )-p2mup4*p3mup2*p3mup5*p4mup1*F(2848)-p3mup1*p3mup2*p3mup5*p
     -   4mup4*F(4293)+p3mup1*p3mup2*p3mup4*p4mup5*F(4293)-p2mup1*p3m
     -   up5*p4mup2*p4mup4*F(4295)+p2mup1*p3mup4*p4mup2*p4mup5*F(4295
     -   )-p3mup2*p3mup5*p4mup1*p4mup4*F(4296)-p3mup1*p3mup5*p4mup2*p
     -   4mup4*F(4296)+p3mup2*p3mup4*p4mup1*p4mup5*F(4296)+p3mup1*p3m
     -   up4*p4mup2*p4mup5*F(4296)-p2mup1*p3mup2*p3mup5*p4mup4*F(4298
     -   )+p2mup1*p3mup2*p3mup4*p4mup5*F(4298)-p3mup5*p4mup1*p4mup2*p
     -   4mup4*F(4300)+p3mup4*p4mup1*p4mup2*p4mup5*F(4300))+2*(mup4mu
     -   p5*p1mup1*p1mup2*F(2829)+mup4mup5*p2mup1*p2mup2*F(2844)+mup2
     -   mup5*p2mup1*p2mup4*F(2845)+mup2mup4*p2mup1*p2mup5*F(2846)+mu
     -   p4mup5*p3mup1*p3mup2*F(2847)+mup1mup5*p2mup2*p2mup4*F(932)*P
     -   (27)+p1mup5*(mup1mup4*p1mup2*F(2831)-2*p4mup2*(p3mup4*p4mup1
     -   *F(1458)+p3mup1*p4mup4*F(4288))+mup2mup4*p1mup1*F(1047)*P(41
     -   ))+p1mup4*(2*p3mup5*p4mup1*p4mup2*F(1458)+mup1mup5*p1mup2*F(
     -   2830)-mup2mup5*p1mup1*F(1047)*P(41))+mup1mup4*p2mup2*p2mup5*
     -   F(932)*P(386))
       Fa8(75) = mup1mup4*p2mup2*p3mup5*F(2836)+mup1mup5*p2mup2*p3mup4
     -   *F(2837)+mup1mup4*p2mup5*p3mup2*F(2839)+mup1mup2*p2mup5*p3mu
     -   p4*F(2840)+mup1mup5*p2mup4*p3mup2*F(2842)+mup1mup2*p2mup4*p3
     -   mup5*F(2843)+mup4mup5*p2mup1*p3mup2*F(2852)+mup2mup5*p2mup1*
     -   p3mup4*F(2853)+mup2mup4*p2mup1*p3mup5*F(2854)-mup4mup5*p2mup
     -   2*p3mup1*F(2855)+mup2mup5*p2mup4*p3mup1*F(2856)+mup2mup4*p2m
     -   up5*p3mup1*F(2857)+4*(p2mup5*p3mup1*p3mup4*p4mup2*F(1582)-p2
     -   mup4*p3mup1*p3mup5*p4mup2*F(1582)-p2mup2*p3mup1*p3mup5*p4mup
     -   4*F(1582)+p2mup2*p3mup1*p3mup4*p4mup5*F(1582)+p2mup1*p2mup5*
     -   p3mup4*p4mup2*F(1585)-p2mup1*p2mup4*p3mup5*p4mup2*F(1585)-p2
     -   mup1*p2mup2*p3mup5*p4mup4*F(1585)+p2mup1*p2mup2*p3mup4*p4mup
     -   5*F(1585)+p2mup5*p3mup4*p4mup1*p4mup2*F(1701)-p2mup4*p3mup5*
     -   p4mup1*p4mup2*F(1701)-p2mup2*p3mup5*p4mup1*p4mup4*F(1701)+p2
     -   mup2*p3mup4*p4mup1*p4mup5*F(1701)+p2mup5*p3mup1*p3mup2*p4mup
     -   4*F(4286)-p2mup4*p3mup1*p3mup2*p4mup5*F(4286)+p2mup5*p3mup2*
     -   p4mup1*p4mup4*F(4288)-p2mup4*p3mup2*p4mup1*p4mup5*F(4288)+p2
     -   mup1*p2mup5*p3mup2*p4mup4*F(4290)-p2mup1*p2mup4*p3mup2*p4mup
     -   5*F(4290))+8*(p2mup2*p2mup5*p3mup4*p4mup1*F(1587)-p2mup2*p2m
     -   up4*p3mup5*p4mup1*F(1587)+p2mup5*p3mup2*p3mup4*p4mup1*F(2848
     -   )-p2mup4*p3mup2*p3mup5*p4mup1*F(2848)-p3mup1*p3mup2*p3mup5*p
     -   4mup4*F(4293)+p3mup1*p3mup2*p3mup4*p4mup5*F(4293)-p2mup1*p3m
     -   up5*p4mup2*p4mup4*F(4295)+p2mup1*p3mup4*p4mup2*p4mup5*F(4295
     -   )-p3mup2*p3mup5*p4mup1*p4mup4*F(4296)-p3mup1*p3mup5*p4mup2*p
     -   4mup4*F(4296)+p3mup2*p3mup4*p4mup1*p4mup5*F(4296)+p3mup1*p3m
     -   up4*p4mup2*p4mup5*F(4296)-p2mup1*p3mup2*p3mup5*p4mup4*F(4298
     -   )+p2mup1*p3mup2*p3mup4*p4mup5*F(4298)-p3mup5*p4mup1*p4mup2*p
     -   4mup4*F(4300)+p3mup4*p4mup1*p4mup2*p4mup5*F(4300))+2*(mup4mu
     -   p5*p1mup1*p1mup2*F(2849)+mup4mup5*p2mup1*p2mup2*F(2858)+mup2
     -   mup5*p2mup1*p2mup4*F(2859)+mup2mup4*p2mup1*p2mup5*F(2860)-mu
     -   p4mup5*p3mup1*p3mup2*F(2861)+mup1mup4*p2mup2*p2mup5*F(932)*P
     -   (27)+p1mup5*(mup1mup4*p1mup2*F(2851)+2*(p3mup4*p4mup1*p4mup2
     -   *F(1458)+p3mup1*p4mup2*p4mup4*F(4288))-mup2mup4*p1mup1*F(104
     -   7)*P(41))+p1mup4*(-2*p3mup5*p4mup1*p4mup2*F(1458)+mup1mup5*p
     -   1mup2*F(2850)+mup2mup5*p1mup1*F(1047)*P(41))+mup1mup5*p2mup2
     -   *p2mup4*F(932)*P(386))
       Fa8(76) = -8*p2mup4*p2mup5*p3mup2*p4mup1*F(1580)-4*(p2mup5*p3mu
     -   p1*p3mup2*p4mup4*F(1464)+p2mup4*p3mup1*p3mup2*p4mup5*F(1464)
     -   +p2mup5*p3mup2*p4mup1*p4mup4*F(1469)+p1mup5*p3mup1*p4mup2*p4
     -   mup4*F(1469)+p2mup4*p3mup2*p4mup1*p4mup5*F(1469)+p2mup1*p2mu
     -   p5*p3mup2*p4mup4*F(1592)+p2mup1*p2mup4*p3mup2*p4mup5*F(1592)
     -   )+mup4mup5*p2mup1*p3mup2*F(2865)+mup4mup5*p2mup2*p3mup1*F(28
     -   68)+mup2mup5*p2mup4*p3mup1*F(2871)+mup1mup5*p2mup4*p3mup2*F(
     -   2872)+p3mup5*(16*p3mup1*p3mup4*p4mup2*F(1700)+8*(p3mup1*p3mu
     -   p2*p4mup4*F(1708)+p3mup1*p4mup2*p4mup4*F(1916))+4*(p2mup2*p3
     -   mup1*p4mup4*F(1594)+p2mup4*p3mup1*p4mup2*F(1703)+p4mup1*(p2m
     -   up2*(2*(p2mup4*F(1579)+p3mup4*F(1581))+p4mup4*F(1705))+p4mup
     -   2*(2*p4mup4*F(1059)+p1mup4*F(1468)+p2mup4*F(1706)+4*p3mup4*F
     -   (1848))+2*p3mup2*(2*p3mup4*F(1700)+p4mup4*F(1916)+p2mup4*F(2
     -   848))))+p2mup1*(8*p3mup2*p4mup4*F(1595)+4*(p2mup2*p4mup4*F(1
     -   590)+p4mup2*(4*p3mup4*F(1588)+p2mup4*F(1591)+2*p4mup4*F(1709
     -   )))+mup2mup4*F(2867))+mup1mup4*p2mup2*F(2870)+mup1mup2*p2mup
     -   4*F(2873))+mup2mup4*p2mup5*p3mup1*F(2874)+mup1mup4*p2mup5*p3
     -   mup2*F(2875)+p3mup4*(8*(p3mup1*p3mup2*p4mup5*F(1708)+p3mup1*
     -   p4mup2*p4mup5*F(1916))+p2mup1*(4*(p2mup2*p4mup5*F(1590)+p2mu
     -   p5*p4mup2*F(1593))+8*(p3mup2*p4mup5*F(1595)+p4mup2*p4mup5*F(
     -   1709))+mup2mup5*F(2866))+mup1mup5*p2mup2*F(2869)+mup1mup2*p2
     -   mup5*F(2876)+4*(p2mup2*p3mup1*p4mup5*F(1594)+p2mup5*p3mup1*p
     -   4mup2*F(1704)+p4mup1*(p2mup2*p4mup5*F(1705)+p4mup2*(2*p4mup5
     -   *F(1059)+p1mup5*F(1470)+p2mup5*F(1707))+2*(p2mup2*p2mup5*F(1
     -   587)+p3mup2*p4mup5*F(1916)+p2mup5*p3mup2*F(2881)))))+2*(mup4
     -   mup5*p1mup1*p1mup2*F(2862)+mup1mup5*p1mup2*p1mup4*F(2863)+mu
     -   p1mup4*p1mup2*p1mup5*F(2864)+mup4mup5*p2mup1*p2mup2*F(2877)+
     -   mup2mup5*p2mup1*p2mup4*F(2878)+mup2mup4*p2mup1*p2mup5*F(2879
     -   )+mup4mup5*p3mup1*p3mup2*F(2880)+mup1mup4*p2mup2*p2mup5*F(93
     -   2)*P(27)+mup1mup2*p2mup4*p2mup5*F(1163)*P(27)+mup1mup5*p2mup
     -   2*p2mup4*F(1165)*P(27)-(mup2mup5*p1mup1*p1mup4*F(929)-mup1mu
     -   p2*p1mup4*p1mup5*F(1039)+mup2mup4*p1mup1*p1mup5*F(1047))*P(4
     -   1))
       Fa8(77) = (mup1mup5*mup2mup4*F(2884))/2.+mup4mup5*p1mup1*p2mup2
     -   *F(2885)+mup2mup5*p1mup1*p2mup4*F(2886)+mup2mup4*p1mup1*p2mu
     -   p5*F(2887)+mup4mup5*p1mup2*p2mup1*F(2888)+mup1mup5*p1mup2*p2
     -   mup4*F(2889)+mup2mup5*p1mup4*p2mup1*F(2891)+mup1mup5*p1mup4*
     -   p2mup2*F(2892)+mup2mup4*p1mup5*p2mup1*F(2894)+mup4mup5*p1mup
     -   1*p3mup2*F(2897)+mup2mup5*p1mup1*p3mup4*F(2898)+mup2mup4*p1m
     -   up1*p3mup5*F(2899)+mup4mup5*p1mup2*p3mup1*F(2900)+mup1mup5*p
     -   1mup2*p3mup4*F(2901)+mup2mup5*p1mup4*p3mup1*F(2903)+mup1mup5
     -   *p1mup4*p3mup2*F(2904)+mup2mup4*p1mup5*p3mup1*F(2906)+mup4mu
     -   p5*p1mup1*p4mup2*F(2909)+mup2mup5*p1mup1*p4mup4*F(2910)+mup2
     -   mup4*p1mup1*p4mup5*F(2911)+mup4mup5*p1mup2*p4mup1*F(2912)+mu
     -   p1mup5*p1mup2*p4mup4*F(2913)+mup2mup5*p1mup4*p4mup1*F(2915)+
     -   mup1mup5*p1mup4*p4mup2*F(2916)+mup2mup4*p1mup5*p4mup1*F(2918
     -   )+mup4mup5*p2mup1*p4mup2*F(2921)+mup2mup5*p2mup1*p4mup4*F(29
     -   22)+mup2mup4*p2mup1*p4mup5*F(2923)+mup4mup5*p2mup2*p4mup1*F(
     -   2924)+mup1mup5*p2mup2*p4mup4*F(2925)+mup2mup5*p2mup4*p4mup1*
     -   F(2927)+mup1mup5*p2mup4*p4mup2*F(2928)+mup2mup4*p2mup5*p4mup
     -   1*F(2930)+mup2mup5*p3mup1*p3mup4*F(2933)+mup1mup5*p3mup2*p3m
     -   up4*F(2934)+mup4mup5*p3mup1*p4mup2*F(2936)+mup2mup5*p3mup1*p
     -   4mup4*F(2937)+mup2mup4*p3mup1*p4mup5*F(2938)+mup4mup5*p3mup2
     -   *p4mup1*F(2939)+mup1mup5*p3mup2*p4mup4*F(2940)+mup2mup5*p3mu
     -   p4*p4mup1*F(2942)+mup1mup5*p3mup4*p4mup2*F(2943)+mup2mup4*p3
     -   mup5*p4mup1*F(2945)+mup1mup4*((mup2mup5*F(2883))/2.+p1mup2*p
     -   2mup5*F(2890)+p1mup5*p2mup2*F(2895)+p1mup2*p3mup5*F(2902)+p1
     -   mup5*p3mup2*F(2907)+p1mup2*p4mup5*F(2914)+p1mup5*p4mup2*F(29
     -   19)+p2mup2*p4mup5*F(2926)+p2mup5*p4mup2*F(2931)+p3mup2*p4mup
     -   5*F(2941)+p3mup5*p4mup2*F(2946)+2*(p3mup2*p3mup5*F(2949)+p4m
     -   up2*p4mup5*F(2953)))+mup1mup2*((mup4mup5*F(2882))/2.+p1mup4*
     -   p2mup5*F(2893)+p1mup5*p2mup4*F(2896)+p1mup4*p3mup5*F(2905)+p
     -   1mup5*p3mup4*F(2908)+p1mup4*p4mup5*F(2917)+p1mup5*p4mup4*F(2
     -   920)+p2mup4*p4mup5*F(2929)+p2mup5*p4mup4*F(2932)+p3mup4*p3mu
     -   p5*F(2935)+p3mup4*p4mup5*F(2944)+p3mup5*p4mup4*F(2947)-2*p4m
     -   up4*p4mup5*F(2954))-2*(mup2mup4*p3mup1*p3mup5*F(2948)-mup2mu
     -   p5*p4mup1*p4mup4*F(2950)-mup2mup4*p4mup1*p4mup5*F(2951)-mup1
     -   mup5*p4mup2*p4mup4*F(2952)-mup4mup5*p4mup1*p4mup2*F(2955))+p
     -   1mup1*p1mup2*p1mup4*p3mup5*F(2957)+8*(p1mup1*p1mup2*p1mup4*p
     -   2mup5*F(2956)+p1mup1*p1mup2*p1mup5*p3mup4*F(2959)+p1mup1*p1m
     -   up4*p1mup5*p2mup2*F(2960)+p1mup1*p1mup2*p1mup4*p4mup5*F(2961
     -   )-p1mup1*p1mup4*p1mup5*p3mup2*F(2962)+p1mup2*p1mup4*p1mup5*p
     -   2mup1*F(2963)+p1mup2*p1mup4*p1mup5*p3mup1*F(2964)+p2mup1*p2m
     -   up2*p2mup4*p3mup5*F(2966)+p2mup1*p2mup2*p2mup5*p3mup4*F(2967
     -   )-p2mup1*p2mup2*p3mup4*p3mup5*F(2968)+p2mup1*p2mup4*p2mup5*p
     -   3mup2*F(2969)-p2mup2*p2mup5*p3mup1*p3mup4*F(2970)-p2mup1*p2m
     -   up4*p3mup2*p3mup5*F(2970)+p2mup1*p2mup5*p3mup2*p3mup4*F(2971
     -   )-p2mup2*p2mup4*p3mup1*p3mup5*F(2974)-p2mup2*p3mup1*p3mup4*p
     -   3mup5*F(2975)+p2mup4*p2mup5*p3mup1*p3mup2*F(2976)-p2mup4*p3m
     -   up1*p3mup2*p3mup5*F(2977)+p2mup5*p3mup1*p3mup2*p3mup4*F(2978
     -   ))+16*(p1mup1*p1mup2*p1mup5*p2mup4*F(2958)+p2mup1*p2mup2*p2m
     -   up4*p2mup5*F(2965)-p2mup1*p3mup2*p3mup4*p3mup5*F(2972)-p2mup
     -   2*p2mup4*p2mup5*p3mup1*F(2973)-p3mup1*p3mup2*p3mup4*p3mup5*F
     -   (2979)+p4mup1*p4mup2*p4mup4*p4mup5*F(2980))
       Fa8(78) = mup1mup2*mup4mup5*F(2981)+mup1mup4*mup2mup5*F(2982)+m
     -   up1mup5*mup2mup4*F(2983)-16*(p2mup2*p2mup5*p3mup1*p3mup4*F(1
     -   583)-p2mup2*p2mup4*p3mup1*p3mup5*F(1583)+p2mup5*p3mup1*p3mup
     -   2*p3mup4*F(3034)-p2mup4*p3mup1*p3mup2*p3mup5*F(3034))-4*(mup
     -   2mup5*p3mup1*p3mup4*F(3008)-mup2mup4*p3mup1*p3mup5*F(3008)-m
     -   up1mup5*p3mup2*p3mup4*F(3009)-mup1mup4*p3mup2*p3mup5*F(3010)
     -   -mup4mup5*p4mup1*p4mup2*F(3028)-mup1mup5*p4mup2*p4mup4*F(428
     -   8)*P(27)+mup2mup5*p4mup1*p4mup4*F(4288)*P(41)-mup2mup4*p4mup
     -   1*p4mup5*F(4288)*P(41))+2*(mup1mup2*p3mup5*p4mup4*F(1842)+mu
     -   p1mup2*p3mup4*p4mup5*F(1874)+mup1mup5*p1mup2*p2mup4*F(2988)+
     -   mup1mup4*p1mup2*p2mup5*F(2989)+mup1mup5*p1mup4*p2mup2*F(2991
     -   )+mup1mup2*p1mup4*p2mup5*F(2992)+mup1mup4*p1mup5*p2mup2*F(29
     -   94)+mup1mup2*p1mup5*p2mup4*F(2995)+mup4mup5*p1mup2*p3mup1*F(
     -   2999)-mup1mup5*p1mup2*p3mup4*F(3000)+mup1mup4*p1mup2*p3mup5*
     -   F(3001)+mup2mup5*p1mup4*p3mup1*F(3002)+mup1mup5*p1mup4*p3mup
     -   2*F(3003)-mup1mup2*p1mup4*p3mup5*F(3004)+mup2mup4*p1mup5*p3m
     -   up1*F(3005)+mup1mup4*p1mup5*p3mup2*F(3006)+mup1mup2*p1mup5*p
     -   3mup4*F(3007)+mup4mup5*p1mup2*p4mup1*F(3012)+mup1mup5*p1mup4
     -   *p4mup2*F(3013)+mup1mup4*p1mup5*p4mup2*F(3014)+mup4mup5*p2mu
     -   p2*p4mup1*F(3015)+mup2mup5*p2mup4*p4mup1*F(3016)+mup2mup4*p2
     -   mup5*p4mup1*F(3017)+mup4mup5*p3mup2*p4mup1*F(3018)+mup2mup5*
     -   p3mup4*p4mup1*F(3019)+mup2mup4*p3mup5*p4mup1*F(3020)+mup2mup
     -   5*p3mup1*p4mup4*F(3021)+mup2mup4*p3mup1*p4mup5*F(3022)+mup1m
     -   up5*p3mup2*p4mup4*F(3023)+mup1mup4*p3mup2*p4mup5*F(3024)+mup
     -   4mup5*p3mup1*p4mup2*F(3025)+mup1mup5*p3mup4*p4mup2*F(3026)+m
     -   up1mup4*p3mup5*p4mup2*F(3027)+mup1mup5*p1mup2*p4mup4*F(3030)
     -   +mup1mup4*p1mup2*p4mup5*F(3031)+p2mup1*(-8*(p2mup5*p3mup4-p2
     -   mup4*p3mup5)*(p2mup2*F(933)+p3mup2*F(1583))+mup4mup5*p1mup2*
     -   F(2987)+mup2mup5*p1mup4*F(2990)+mup2mup4*p1mup5*F(2993)+mup4
     -   mup5*p4mup2*F(3029)+mup2mup5*p4mup4*F(3032)+mup2mup4*p4mup5*
     -   F(3033))+mup1mup5*p2mup4*p4mup2*F(1320)*P(27)+mup1mup5*p2mup
     -   2*p4mup4*F(1320)*P(27)+mup1mup2*p2mup5*p4mup4*F(4285)*P(27)-
     -   mup2mup5*p1mup4*p4mup1*F(1314)*P(41)+mup2mup4*p1mup5*p4mup1*
     -   F(1314)*P(41)+mup1mup2*p1mup5*p4mup4*F(4285)*P(41)-mup1mup2*
     -   p1mup4*p4mup5*F(4285)*P(41)+p1mup1*(-8*p1mup2*(p1mup5*p3mup4
     -   -p1mup4*p3mup5)*F(930)+mup4mup5*p2mup2*F(2984)+mup2mup5*p2mu
     -   p4*F(2985)+mup2mup4*p2mup5*F(2986)+mup4mup5*p3mup2*F(2996)-m
     -   up2mup5*p3mup4*F(2997)+mup2mup4*p3mup5*F(2998)+mup4mup5*p4mu
     -   p2*F(3011)-(mup2mup5*p4mup4-mup2mup4*p4mup5)*F(1314)*P(41))+
     -   (mup1mup4*(p2mup5*p4mup2+p2mup2*p4mup5)*F(1320)+mup1mup2*p2m
     -   up4*p4mup5*F(4285)+2*mup1mup4*p4mup2*p4mup5*F(4288))*P(386))
       Fa8(79) = 16*(p2mup2*p2mup5*p3mup1*p3mup4*F(1583)-p2mup2*p2mup4
     -   *p3mup1*p3mup5*F(1583)+p2mup5*p3mup1*p3mup2*p3mup4*F(3034)-p
     -   2mup4*p3mup1*p3mup2*p3mup5*F(3034))+mup1mup2*mup4mup5*F(3035
     -   )+mup1mup4*mup2mup5*F(3036)+mup1mup5*mup2mup4*F(3037)+4*(mup
     -   1mup4*p3mup2*p3mup5*F(3009)+mup1mup5*p3mup2*p3mup4*F(3010)-m
     -   up2mup5*p3mup1*p3mup4*F(3062)+mup2mup4*p3mup1*p3mup5*F(3062)
     -   +mup4mup5*p4mup1*p4mup2*F(3074)+mup1mup4*p4mup2*p4mup5*F(428
     -   8)*P(27)+mup2mup5*p4mup1*p4mup4*F(4288)*P(41)-mup2mup4*p4mup
     -   1*p4mup5*F(4288)*P(41))+2*(mup1mup2*p3mup4*p4mup5*F(1842)+mu
     -   p1mup2*p3mup5*p4mup4*F(1874)+mup2mup4*p3mup1*p4mup5*F(3021)+
     -   mup2mup5*p3mup1*p4mup4*F(3022)+mup1mup4*p3mup2*p4mup5*F(3023
     -   )+mup1mup5*p3mup2*p4mup4*F(3024)+mup1mup4*p3mup5*p4mup2*F(30
     -   26)+mup1mup5*p3mup4*p4mup2*F(3027)+mup1mup4*p1mup2*p4mup5*F(
     -   3030)+mup1mup5*p1mup2*p4mup4*F(3031)+mup1mup5*p1mup2*p2mup4*
     -   F(3042)+mup1mup4*p1mup2*p2mup5*F(3043)+mup1mup5*p1mup4*p2mup
     -   2*F(3045)+mup1mup2*p1mup4*p2mup5*F(3046)+mup1mup4*p1mup5*p2m
     -   up2*F(3048)+mup1mup2*p1mup5*p2mup4*F(3049)+mup4mup5*p1mup2*p
     -   3mup1*F(3053)+mup1mup5*p1mup2*p3mup4*F(3054)+mup1mup4*p1mup2
     -   *p3mup5*F(3055)+mup2mup5*p1mup4*p3mup1*F(3056)+mup1mup5*p1mu
     -   p4*p3mup2*F(3057)+mup1mup2*p1mup4*p3mup5*F(3058)+mup2mup4*p1
     -   mup5*p3mup1*F(3059)+mup1mup4*p1mup5*p3mup2*F(3060)+mup1mup2*
     -   p1mup5*p3mup4*F(3061)+mup4mup5*p1mup2*p4mup1*F(3064)+mup1mup
     -   5*p1mup4*p4mup2*F(3065)+mup1mup4*p1mup5*p4mup2*F(3066)+mup4m
     -   up5*p2mup2*p4mup1*F(3067)+mup2mup5*p2mup4*p4mup1*F(3068)+mup
     -   2mup4*p2mup5*p4mup1*F(3069)+mup4mup5*p3mup2*p4mup1*F(3070)+m
     -   up2mup5*p3mup4*p4mup1*F(3071)+mup2mup4*p3mup5*p4mup1*F(3072)
     -   +mup4mup5*p3mup1*p4mup2*F(3073)+p2mup1*(8*(p2mup5*p3mup4-p2m
     -   up4*p3mup5)*(p2mup2*F(933)+p3mup2*F(1583))+mup2mup4*p4mup5*F
     -   (3032)+mup2mup5*p4mup4*F(3033)+mup4mup5*p1mup2*F(3041)+mup2m
     -   up5*p1mup4*F(3044)+mup2mup4*p1mup5*F(3047)+mup4mup5*p4mup2*F
     -   (3075))+mup1mup4*p2mup5*p4mup2*F(1320)*P(27)+mup1mup4*p2mup2
     -   *p4mup5*F(1320)*P(27)+mup1mup2*p2mup4*p4mup5*F(4285)*P(27)+m
     -   up2mup5*p1mup4*p4mup1*F(1314)*P(41)-mup2mup4*p1mup5*p4mup1*F
     -   (1314)*P(41)-mup1mup2*p1mup5*p4mup4*F(4285)*P(41)+mup1mup2*p
     -   1mup4*p4mup5*F(4285)*P(41)+p1mup1*(8*p1mup2*(p1mup5*p3mup4-p
     -   1mup4*p3mup5)*F(930)+mup4mup5*p2mup2*F(3038)+mup2mup5*p2mup4
     -   *F(3039)+mup2mup4*p2mup5*F(3040)+mup4mup5*p3mup2*F(3050)+mup
     -   2mup5*p3mup4*F(3051)+mup2mup4*p3mup5*F(3052)+mup4mup5*p4mup2
     -   *F(3063)+(mup2mup5*p4mup4-mup2mup4*p4mup5)*F(1314)*P(41))+(m
     -   up1mup5*(p2mup4*p4mup2+p2mup2*p4mup4)*F(1320)+mup1mup2*p2mup
     -   5*p4mup4*F(4285)+2*mup1mup5*p4mup2*p4mup4*F(4288))*P(386))
       Fa8(80) = 32*(p3mup1*p3mup2*p3mup4*p3mup5*F(1054)+p2mup1*p3mup2
     -   *p3mup4*p3mup5*F(1584))+mup1mup2*mup4mup5*F(3076)+mup1mup4*m
     -   up2mup5*F(3077)+mup1mup5*mup2mup4*F(3078)+16*(p2mup1*(p2mup2
     -   *p2mup4*p3mup5*F(931)-p2mup4*p2mup5*p3mup2*F(932)+p2mup2*p2m
     -   up5*p3mup4*F(933))+p2mup1*p2mup2*p3mup4*p3mup5*F(1447)+p2mup
     -   2*p2mup4*p3mup1*p3mup5*F(1573)-p2mup4*p2mup5*p3mup1*p3mup2*F
     -   (1574)+p2mup2*p3mup1*p3mup4*p3mup5*F(1575)+p2mup2*p2mup5*p3m
     -   up1*p3mup4*F(1583)+p2mup1*p2mup4*p3mup2*p3mup5*F(1583)+p2mup
     -   4*p3mup1*p3mup2*p3mup5*F(3034)+p2mup1*p2mup5*p3mup2*p3mup4*F
     -   (3127)+p2mup5*p3mup1*p3mup2*p3mup4*F(3128))+4*(mup1mup2*p3mu
     -   p4*p3mup5*F(1835)+mup2mup5*p3mup1*p3mup4*F(3103)+mup1mup5*p3
     -   mup2*p3mup4*F(3104)-mup2mup4*p3mup1*p3mup5*F(3113)+mup1mup4*
     -   p3mup2*p3mup5*F(3114)+mup4mup5*p4mup1*p4mup2*F(3122)+mup1mup
     -   5*p4mup2*p4mup4*F(1469)*P(27)+mup1mup4*p4mup2*p4mup5*F(1469)
     -   *P(27))+2*(mup1mup2*p3mup5*p4mup4*F(1915)+mup1mup2*p3mup4*p4
     -   mup5*F(1915)+mup2mup5*p1mup4*p2mup1*F(3085)+mup1mup5*p1mup4*
     -   p2mup2*F(3086)+mup1mup2*p1mup4*p2mup5*F(3087)+mup2mup4*p1mup
     -   5*p2mup1*F(3088)+mup1mup4*p1mup5*p2mup2*F(3089)+mup1mup2*p1m
     -   up5*p2mup4*F(3090)-mup2mup5*p1mup4*p3mup1*F(3097)+mup1mup5*p
     -   1mup4*p3mup2*F(3098)+mup1mup2*p1mup4*p3mup5*F(3099)+mup2mup4
     -   *p1mup5*p3mup1*F(3100)-mup1mup4*p1mup5*p3mup2*F(3101)+mup1mu
     -   p2*p1mup5*p3mup4*F(3102)+mup4mup5*p3mup1*p4mup2*F(3105)+mup2
     -   mup5*p3mup1*p4mup4*F(3106)+mup2mup4*p3mup1*p4mup5*F(3106)+mu
     -   p4mup5*p3mup2*p4mup1*F(3107)+mup1mup5*p3mup2*p4mup4*F(3108)+
     -   mup1mup4*p3mup2*p4mup5*F(3108)+mup2mup5*p3mup4*p4mup1*F(3109
     -   )+mup1mup5*p3mup4*p4mup2*F(3110)+mup2mup4*p3mup5*p4mup1*F(31
     -   11)+mup1mup4*p3mup5*p4mup2*F(3112)+mup1mup5*p1mup4*p4mup2*F(
     -   3117)+mup1mup4*p1mup5*p4mup2*F(3118)+mup4mup5*p2mup2*p4mup1*
     -   F(3119)+mup2mup5*p2mup4*p4mup1*F(3120)+mup2mup4*p2mup5*p4mup
     -   1*F(3121)+mup4mup5*p2mup1*p4mup2*F(3123)+p1mup2*(-8*p1mup4*p
     -   1mup5*p3mup1*F(929)+mup4mup5*p2mup1*F(3082)+mup1mup5*p2mup4*
     -   F(3083)+mup1mup4*p2mup5*F(3084)+mup4mup5*p3mup1*F(3094)+mup1
     -   mup5*p3mup4*F(3095)+mup1mup4*p3mup5*F(3096)+mup4mup5*p4mup1*
     -   F(3116)+(mup1mup5*p4mup4+mup1mup4*p4mup5)*F(3124))+mup2mup5*
     -   p2mup1*p4mup4*F(3125)+mup2mup4*p2mup1*p4mup5*F(3125)+mup1mup
     -   4*p2mup5*p4mup2*F(3126)+mup1mup2*p2mup5*p4mup4*F(1333)*P(27)
     -   +mup1mup2*p2mup4*p4mup5*F(1333)*P(27)+mup1mup5*p2mup2*p4mup4
     -   *F(1337)*P(27)+mup1mup4*p2mup2*p4mup5*F(1337)*P(27)+mup1mup5
     -   *p2mup4*p4mup2*F(1339)*P(27)+(mup1mup2*(p1mup5*p4mup4+p1mup4
     -   *p4mup5)*F(1333)+p4mup1*(-(mup2mup5*p1mup4*F(1338))+mup2mup4
     -   *p1mup5*F(1346)-2*(mup2mup5*p4mup4+mup2mup4*p4mup5)*F(1469))
     -   )*P(41)+p1mup1*(8*p1mup2*(p1mup4*p3mup5*F(928)+p1mup5*p3mup4
     -   *F(930))+mup4mup5*p2mup2*F(3079)+mup2mup5*p2mup4*F(3080)+mup
     -   2mup4*p2mup5*F(3081)+mup4mup5*p3mup2*F(3091)+mup2mup5*p3mup4
     -   *F(3092)-mup2mup4*p3mup5*F(3093)+mup4mup5*p4mup2*F(3115)-(mu
     -   p2mup5*p4mup4+mup2mup4*p4mup5)*F(1329)*P(41)))
       Return
       End
