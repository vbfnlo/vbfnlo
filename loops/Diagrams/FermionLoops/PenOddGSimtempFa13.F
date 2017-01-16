c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenOddGSim13(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4
     -   ,p3mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5m
     -   up3,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mu
     -   p3,mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa13)
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
       Complex*16 Fa13(1189:1287)
       COMMON/FaPenOddGSim/Fa
       Fa13(1189) = mup1mup4*((p1mup5+p2mup5)*F(2561)+p3mup5*F(2708)+p
     -   4mup5*F(2843))
       Fa13(1190) = mup1mup4*(p1mup5*F(2527)+p4mup5*F(2531)+(p2mup5+p3
     -   mup5)*F(2566))
       Fa13(1191) = mup1mup5*F(2524)+4*(p1mup5*(p2mup1*F(2527)+p3mup1*
     -   F(2528))+p2mup1*p4mup5*F(2531)+p2mup1*(p2mup5+p3mup5)*F(2566
     -   )+p3mup1*(p4mup5*F(2532)-p2mup5*F(2578)+p3mup5*F(2733)))
       Fa13(1192) = mup1mup5*F(2572)+4*p2mup1*((p1mup5+p2mup5)*F(2561)
     -   +p3mup5*F(2708)+p4mup5*F(2843))
       Fa13(1193) = (p1mup3+p2mup3)*p3mup1*F(2559)-p3mup1*p3mup3*F(256
     -   2)
       Fa13(1194) = p2mup1*(p1mup3+p2mup3)*F(2559)-p2mup1*p3mup3*F(256
     -   2)
       Fa13(1195) = mup1mup4*(p1mup5*F(2528)+p4mup5*F(2532)-p2mup5*F(2
     -   578)+p3mup5*F(2733))
       Fa13(1196) = p2mup1*((p1mup4+p2mup4)*F(2560)+p3mup4*F(2842)-p4m
     -   up4*F(2844))
       Fa13(1197) = p3mup1*((p1mup4+p2mup4)*F(2560)+p3mup4*F(2842)-p4m
     -   up4*F(2844))
       Fa13(1198) = -(mup1mup4*(p1mup3+p2mup3)*F(2559))+mup1mup4*p3mup
     -   3*F(2562)+mup1mup3*((p1mup4+p2mup4)*F(2560)+p3mup4*F(2842)-p
     -   4mup4*F(2844))
       Fa13(1199) = p4mup1*(p1mup4*F(2526)-p4mup4*F(2536)-p2mup4*F(257
     -   7)+p3mup4*F(2847))
       Fa13(1200) = p4mup1*(p1mup4*F(2525)-p4mup4*F(2535)+(p2mup4+p3mu
     -   p4)*F(2565))
       Fa13(1201) = p1mup4*(p2mup1*F(2525)+p3mup1*F(2526))-p2mup1*p4mu
     -   p4*F(2535)+p2mup1*(p2mup4+p3mup4)*F(2565)-p3mup1*(p4mup4*F(2
     -   536)+p2mup4*F(2577)-p3mup4*F(2847))
       Fa13(1202) = -(mup1mup4*p1mup3*F(2564))+mup1mup3*(p1mup4*F(2525
     -   )-p4mup4*F(2535)+(p2mup4+p3mup4)*F(2565))+mup1mup4*(p2mup3+p
     -   3mup3)*F(2567)
       Fa13(1203) = mup1mup4*(-(p1mup3*F(2575))+p2mup3*F(2576)+p3mup3*
     -   F(2579))+mup1mup3*(p1mup4*F(2526)-p4mup4*F(2536)-p2mup4*F(25
     -   77)+p3mup4*F(2847))
       Fa13(1204) = -(p2mup1*(p2mup3+p3mup3)*F(2567))+p1mup3*(p2mup1*F
     -   (2564)+p3mup1*F(2575))-p3mup1*(p2mup3*F(2576)+p3mup3*F(2579)
     -   )
       Fa13(1205) = mup1mup5*F(2573)+4*p3mup1*((p1mup5+p2mup5)*F(2561)
     -   +p3mup5*F(2708)+p4mup5*F(2843))
       Fa13(1206) = mup1mup5*F(2582)-4*(p1mup5*(p2mup1*F(2527)+p3mup1*
     -   F(2528))+p2mup1*p4mup5*F(2531)+p2mup1*(p2mup5+p3mup5)*F(2566
     -   )+p3mup1*(p4mup5*F(2532)-p2mup5*F(2578)+p3mup5*F(2733)))
       Fa13(1207) = mup1mup5*F(2574)-4*p2mup1*((p1mup5+p2mup5)*F(2561)
     -   +p3mup5*F(2708)+p4mup5*F(2843))
       Fa13(1208) = p2mup1*(p2mup3+p3mup3)*F(2567)-p1mup3*(p2mup1*F(25
     -   64)+p3mup1*F(2575))+p3mup1*(p2mup3*F(2576)+p3mup3*F(2579))
       Fa13(1209) = p4mup1*(p1mup5*F(2527)+p4mup5*F(2531)+(p2mup5+p3mu
     -   p5)*F(2566))
       Fa13(1210) = p4mup1*((p1mup5+p2mup5)*F(2561)+p3mup5*F(2708)+p4m
     -   up5*F(2843))
       Fa13(1211) = p4mup1*(p1mup5*F(2528)+p4mup5*F(2532)-p2mup5*F(257
     -   8)+p3mup5*F(2733))
       Fa13(1212) = p4mup1*((p1mup3+p2mup3)*F(2559)-p3mup3*F(2562)+p4m
     -   up3*F(2845))
       Fa13(1213) = p4mup1*(-(p1mup3*F(2564))+(p2mup3+p3mup3)*F(2567)-
     -   p4mup3*F(2846))
       Fa13(1214) = p4mup1*(p1mup3*F(2575)-p2mup3*F(2576)-p3mup3*F(257
     -   9)+p4mup3*F(2848))
       Fa13(1215) = p4mup1*(p1mup3*F(2564)-(p2mup3+p3mup3)*F(2567)+p4m
     -   up3*F(2846))
       Fa13(1216) = mup1mup5*((p1mup4+p2mup4)*F(2560)+p3mup4*F(2842)-p
     -   4mup4*F(2844))
       Fa13(1217) = mup1mup5*(p1mup4*F(2525)-p4mup4*F(2535)+(p2mup4+p3
     -   mup4)*F(2565))
       Fa13(1218) = mup1mup5*(p1mup3*F(2564)-(p2mup3+p3mup3)*F(2567)+p
     -   4mup3*F(2846))
       Fa13(1219) = mup1mup5*((p1mup3+p2mup3)*F(2559)-p3mup3*F(2562)+p
     -   4mup3*F(2845))
       Fa13(1220) = mup1mup5*(p1mup3*F(2575)-p2mup3*F(2576)-p3mup3*F(2
     -   579)+p4mup3*F(2848))
       Fa13(1221) = mup1mup5*(p1mup4*F(2526)-p4mup4*F(2536)-p2mup4*F(2
     -   577)+p3mup4*F(2847))
       Fa13(1222) = mup1mup5*(-(p1mup3*F(2564))+(p2mup3+p3mup3)*F(2567
     -   )-p4mup3*F(2846))
       Fa13(1223) = mup4mup5*F(2882)+p1mup5*p3mup4*F(2886)-8*(p1mup5*p
     -   2mup4*F(2885)+p2mup4*p2mup5*F(2892))+p2mup5*p3mup4*F(2894)+p
     -   3mup4*p4mup5*F(2895)+p3mup5*p4mup4*F(2896)+4*(p1mup5*p4mup4*
     -   F(2888)+p2mup4*p4mup5*F(2889)+p2mup5*p4mup4*F(2890)-p1mup4*(
     -   p3mup5*F(2884)+p4mup5*F(2887)+2*(p1mup5*F(2883)+p2mup5*F(289
     -   1)))-p2mup4*p3mup5*F(2893)-p3mup4*p3mup5*F(2897)+p4mup4*p4mu
     -   p5*F(2898))
       Fa13(1224) = mup4mup5*F(2899)+p1mup4*(p2mup5*F(2900)-4*p1mup5*F
     -   (2902)+p3mup5*F(2903)+p4mup5*F(2904))+p3mup5*p4mup4*F(2907)-
     -   2*p1mup5*p4mup4*F(2909)+p2mup5*p4mup4*F(2913)-4*(p1mup5*p2mu
     -   p4*F(2901)+p3mup4*p4mup5*F(2906)-p1mup5*p3mup4*F(2908)+p2mup
     -   4*(p4mup5*F(2905)+2*(p3mup5*F(2910)+p2mup5*F(2914))))+8*(p3m
     -   up4*p3mup5*F(2911)-p4mup4*p4mup5*F(2912)+p2mup5*p3mup4*F(291
     -   5))
       Fa13(1225) = mup4mup5*F(2916)+p1mup5*p3mup4*F(2920)+2*p1mup4*(p
     -   3mup5*F(2919)+2*(p2mup5*F(2917)+p1mup5*F(2921)+p4mup5*F(2922
     -   )))+p3mup4*p4mup5*F(2924)-4*(p1mup5*p2mup4*F(2918)-p3mup5*p4
     -   mup4*F(2925)-p1mup5*p4mup4*F(2926)-p2mup4*p3mup5*F(2927))+8*
     -   (p2mup4*p4mup5*F(2923)+p3mup4*p3mup5*F(2928)+p4mup4*p4mup5*F
     -   (2929)+p2mup5*p4mup4*F(2930)+p2mup4*p2mup5*F(2931))+p2mup5*p
     -   3mup4*F(2932)
       Fa13(1226) = p2mup5*F(2933)+p3mup5*F(2934)+p4mup5*F(2935)+2*p1m
     -   up5*F(2936)
       Fa13(1227) = p1mup5*F(2937)+p3mup5*F(2938)+p4mup5*F(2939)+2*p2m
     -   up5*F(2940)
       Fa13(1228) = p1mup5*F(2941)+p2mup5*F(2942)+p4mup5*F(2943)-2*p3m
     -   up5*F(2944)
       Fa13(1229) = p1mup5*F(2945)+p2mup5*F(2946)+p3mup5*F(2947)+2*p4m
     -   up5*F(2948)
       Fa13(1230) = p1mup4*F(2949)+p2mup4*F(2950)+p3mup4*F(2951)+p4mup
     -   4*F(2952)
       Fa13(1231) = p1mup4*F(2953)+p2mup4*F(2954)+p3mup4*F(2955)+p4mup
     -   4*F(2956)
       Fa13(1232) = 8*p3mup4*p3mup5*F(2438)+mup4mup5*F(2957)+2*(p2mup4
     -   *p4mup5*F(2960)+p2mup5*p4mup4*F(2961))+4*(p3mup4*p4mup5*F(24
     -   35)+p3mup5*p4mup4*F(2436)+p4mup4*p4mup5*F(2439)+p1mup4*(p4mu
     -   p5*F(2433)+2*(p1mup5*F(2431)+p3mup5*F(2437))+p2mup5*F(2958))
     -   +p1mup5*(2*p3mup4*F(2432)+p4mup4*F(2434)-p2mup4*F(2959))+p2m
     -   up4*p2mup5*F(2962)-p2mup4*p3mup5*F(2963)+p2mup5*p3mup4*F(296
     -   4))
       Fa13(1233) = -8*(p4mup4*p4mup5*F(2462)+p3mup5*p4mup4*F(2464)+p2
     -   mup4*p2mup5*F(2465))+mup4mup5*F(2965)-2*(p1mup4*(2*(p3mup5*F
     -   (2436)+p1mup5*F(2456)+p4mup5*F(2457))+p2mup5*F(2966))-p1mup5
     -   *p2mup4*F(2967))-4*(p1mup5*p4mup4*F(2461)-p2mup4*p4mup5*F(29
     -   68)+p2mup5*p4mup4*F(2969)-p2mup4*p3mup5*F(2970)+p3mup4*(2*(p
     -   4mup5*F(2460)+p3mup5*F(2465))-p1mup5*F(2657)+p2mup5*F(2971))
     -   )
       Fa13(1234) = p1mup2*F(2973)+p3mup2*F(2974)+p4mup2*F(2975)
       Fa13(1235) = p1mup2*F(2976)+p3mup2*F(2977)+p4mup2*F(2978)
       Fa13(1236) = p1mup2*F(2979)+p4mup2*F(2980)+p3mup2*F(2981)
       Fa13(1237) = p1mup2*F(2982)+p3mup2*F(2983)+p4mup2*F(2984)
       Fa13(1238) = p2mup4*F(2985)+p3mup4*F(2986)+p1mup4*F(2987)+p4mup
     -   4*F(2988)
       Fa13(1239) = p2mup4*F(2989)+p1mup4*F(2990)+p3mup4*F(2991)+2*p4m
     -   up4*F(2992)
       Fa13(1240) = mup4mup5*F(3000)+p2mup4*p3mup5*F(3005)+p2mup5*p3mu
     -   p4*F(3006)-4*(p1mup4*p3mup5*F(2525)-p1mup5*p3mup4*F(2527)-p3
     -   mup4*p4mup5*F(2531)-p3mup4*p3mup5*F(2533)-p3mup5*p4mup4*F(25
     -   35)-p1mup4*p2mup5*F(3001)-p1mup5*p2mup4*F(3002)-p2mup4*p4mup
     -   5*F(3003)-p2mup4*p2mup5*F(3004)+p2mup5*p4mup4*F(3007))
       Fa13(1241) = p4mup1*F(3000)+p2mup1*F(3008)+p3mup1*F(3009)
       Fa13(1242) = p3mup1*F(3010)+p4mup1*F(3011)
       Fa13(1243) = p2mup1*F(2999)+p4mup1*F(3012)
       Fa13(1244) = p2mup1*F(3013)+p3mup1*F(3014)
       Fa13(1245) = -4*(p1mup5*p3mup4*F(263)-p1mup4*p3mup5*F(263)-(-(p
     -   2mup5*p3mup4)+p2mup4*p3mup5)*F(379))+2*(p3mup5*p4mup4-p3mup4
     -   *p4mup5)*F(383)+mup4mup5*F(3018)
       Fa13(1246) = mup4mup5*F(3019)+2*(p1mup5*(2*(p3mup4*F(273)+p2mup
     -   4*F(379))+p4mup4*F(4464))-p1mup4*(2*(p3mup5*F(273)+p2mup5*F(
     -   379))+p4mup5*F(4464)))+4*((-(p3mup5*p4mup4)+p3mup4*p4mup5)*F
     -   (4465)-(p2mup5*p4mup4-p2mup4*p4mup5)*F(4469))
       Fa13(1247) = 2*(p1mup5*p3mup4*F(383)-p1mup4*p3mup5*F(383))+mup4
     -   mup5*F(3020)+4*((-(p3mup5*p4mup4)+p3mup4*p4mup5)*F(4467)+(p2
     -   mup5*p3mup4-p2mup4*p3mup5)*F(4469))
       Fa13(1248) = p2mup5*F(3021)+p3mup5*F(3022)-(2*p1mup5*F(263)+p4m
     -   up5*F(383))*P(1)
       Fa13(1249) = p1mup5*F(3023)+p3mup5*F(3024)+2*p2mup5*F(3025)+p4m
     -   up5*F(3026)
       Fa13(1250) = p1mup5*F(3027)-p2mup5*F(3028)+2*p3mup5*F(3029)+p4m
     -   up5*F(3030)
       Fa13(1251) = p2mup5*F(3031)+p3mup5*F(3032)-(p1mup5*F(383)+2*p4m
     -   up5*F(4467))*P(1)
       Fa13(1252) = p1mup4*F(3033)+p3mup4*F(3034)+2*p2mup4*F(3035)+p4m
     -   up4*F(3036)
       Fa13(1253) = 2*p3mup4*F(3029)+p4mup4*F(3030)+p1mup4*F(3037)+p2m
     -   up4*F(3038)
       Fa13(1254) = 4*(p1mup5*p2mup4*F(263)-p1mup4*p2mup5*F(263)-(p2mu
     -   p5*p3mup4-p2mup4*p3mup5)*F(273))-2*(p2mup5*p4mup4-p2mup4*p4m
     -   up5)*F(383)+mup4mup5*F(3039)
       Fa13(1255) = 2*p1mup2*F(3041)+p3mup2*F(3042)+p4mup2*F(3043)
       Fa13(1256) = p1mup2*F(3044)+p3mup2*F(3045)+p4mup2*F(3046)
       Fa13(1257) = -(p1mup2*F(3047))-2*p3mup2*F(3048)+p4mup2*F(3049)
       Fa13(1258) = p1mup2*F(3050)+p3mup2*F(3051)+2*p4mup2*F(3052)
       Fa13(1259) = p2mup4*F(3053)+p3mup4*F(3054)+(2*p1mup4*F(263)+p4m
     -   up4*F(383))*P(1)
       Fa13(1260) = p2mup4*F(3055)+p3mup4*F(3056)+(p1mup4*F(383)+2*p4m
     -   up4*F(4467))*P(1)
       Fa13(1261) = -2*(p2mup5*p3mup4-p2mup4*p3mup5)+mup4mup5*P(1)
       Fa13(1262) = p2mup1*F(3018)+p3mup1*F(3039)+p4mup1*F(4464)*P(1)
       Fa13(1263) = p4mup1*F(3020)+p3mup1*F(3060)
       Fa13(1264) = p4mup1*F(390)+p2mup1*F(3059)
       Fa13(1265) = p3mup1*F(413)+p2mup1*F(3061)
       Fa13(1266) = -4*(p1mup5*p3mup4*F(263)-p1mup4*p3mup5*F(263)-(-(p
     -   2mup5*p3mup4)+p2mup4*p3mup5)*F(379))+2*(p3mup5*p4mup4-p3mup4
     -   *p4mup5)*F(383)+mup4mup5*F(3064)
       Fa13(1267) = mup4mup5*F(3065)+2*(p1mup5*(2*(p3mup4*F(273)+p2mup
     -   4*F(379))+p4mup4*F(4464))-p1mup4*(2*(p3mup5*F(273)+p2mup5*F(
     -   379))+p4mup5*F(4464)))+4*((-(p3mup5*p4mup4)+p3mup4*p4mup5)*F
     -   (4465)-(p2mup5*p4mup4-p2mup4*p4mup5)*F(4469))
       Fa13(1268) = 2*(p1mup5*p3mup4*F(383)-p1mup4*p3mup5*F(383))+mup4
     -   mup5*F(3066)+4*((-(p3mup5*p4mup4)+p3mup4*p4mup5)*F(4467)+(p2
     -   mup5*p3mup4-p2mup4*p3mup5)*F(4469))
       Fa13(1269) = p2mup5*F(3067)+p3mup5*F(3068)-(2*p1mup5*F(263)+p4m
     -   up5*F(383))*P(1)
       Fa13(1270) = p1mup5*F(3069)+p3mup5*F(3070)+2*p2mup5*F(3071)+p4m
     -   up5*F(3072)
       Fa13(1271) = -(p1mup5*F(3073))+p2mup5*F(3074)+2*p3mup5*F(3075)+
     -   p4mup5*F(3076)
       Fa13(1272) = p3mup5*F(3032)+p2mup5*F(3077)-(p1mup5*F(383)+2*p4m
     -   up5*F(4467))*P(1)
       Fa13(1273) = p1mup4*F(3078)+p3mup4*F(3079)+2*p2mup4*F(3080)+p4m
     -   up4*F(3081)
       Fa13(1274) = -2*p3mup4*F(3075)+p1mup4*F(3082)+p2mup4*F(3083)+p4
     -   mup4*F(3084)
       Fa13(1275) = 4*(p1mup5*p2mup4*F(263)-p1mup4*p2mup5*F(263)-(p2mu
     -   p5*p3mup4-p2mup4*p3mup5)*F(273))-2*(p2mup5*p4mup4-p2mup4*p4m
     -   up5)*F(383)+mup4mup5*F(3085)
       Fa13(1276) = 2*p1mup2*F(3087)+p3mup2*F(3088)+p4mup2*F(3089)
       Fa13(1277) = p1mup2*F(3090)+p3mup2*F(3091)+p4mup2*F(3092)
       Fa13(1278) = p1mup2*F(3093)-2*p3mup2*F(3094)+p4mup2*F(3095)
       Fa13(1279) = p1mup2*F(3096)+p3mup2*F(3097)+2*p4mup2*F(3098)
       Fa13(1280) = p2mup4*F(3099)+p3mup4*F(3100)+(2*p1mup4*F(263)+p4m
     -   up4*F(383))*P(1)
       Fa13(1281) = p3mup4*F(3056)+p2mup4*F(3101)+(p1mup4*F(383)+2*p4m
     -   up4*F(4467))*P(1)
       Fa13(1282) = p2mup1*F(3064)+p3mup1*F(3085)+p4mup1*F(4464)*P(1)
       Fa13(1283) = p4mup1*F(3066)+p3mup1*F(3105)
       Fa13(1284) = p4mup1*F(390)+p2mup1*F(3104)
       Fa13(1285) = p3mup1*F(413)+p2mup1*F(3106)
       Fa13(1286) = -4*p3mup5*(p1mup4*F(376)+p2mup4*F(377))+2*(p3mup4*
     -   (2*(p1mup5*F(263)+p2mup5*F(379)-p3mup5*F(396))+p4mup5*F(447)
     -   )-p3mup5*p4mup4*F(448))+mup4mup5*F(3109)
       Fa13(1287) = 2*(p1mup5*p4mup4*F(454)+p1mup4*(2*(p3mup5*F(273)+p
     -   2mup5*F(379)+p1mup5*F(451))+p4mup5*F(454)))+8*(p3mup4*p3mup5
     -   *F(382)+p4mup4*p4mup5*F(459)+p2mup4*p3mup5*F(1500))-mup4mup5
     -   *F(3110)+4*(p3mup5*p4mup4*F(445)+p1mup5*p2mup4*F(453)+p2mup4
     -   *p4mup5*F(455)+p2mup5*(p4mup4*F(444)+2*(p3mup4*F(381)+p2mup4
     -   *F(456)))+p1mup5*p3mup4*F(3111)+p3mup4*p4mup5*F(3112))
       Return
       End
