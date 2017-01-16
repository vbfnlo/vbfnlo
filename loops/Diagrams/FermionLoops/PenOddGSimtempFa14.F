c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenOddGSim14(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4
     -   ,p3mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5m
     -   up3,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mu
     -   p3,mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa14)
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
       Complex*16 Fa14(1288:1386)
       COMMON/FaPenOddGSim/Fa
       Fa14(1288) = -2*(p1mup5*p3mup4*F(443)-p1mup4*p3mup5*F(443))-4*(
     -   (p2mup5*p3mup4-p2mup4*p3mup5)*F(444)-(p3mup5*p4mup4-p3mup4*p
     -   4mup5)*F(459))+mup4mup5*F(3113)
       Fa14(1289) = p2mup5*F(3114)+p3mup5*F(3115)+(2*p1mup5*F(263)+p4m
     -   up5*F(447))*P(1)
       Fa14(1290) = p3mup5*F(3116)+p1mup5*F(3117)+p4mup5*F(3118)+2*p2m
     -   up5*F(3119)
       Fa14(1291) = p2mup5*F(3120)+p1mup5*F(3121)+p4mup5*F(3122)+2*p3m
     -   up5*F(3123)
       Fa14(1292) = p2mup5*F(3124)+p3mup5*F(3125)+(p1mup5*F(443)+2*p4m
     -   up5*F(459))*P(1)
       Fa14(1293) = p3mup4*F(3126)+2*p2mup4*F(3127)+p1mup4*F(3128)+p4m
     -   up4*F(3129)
       Fa14(1294) = p2mup4*F(3130)+p1mup4*F(3131)+2*p3mup4*F(3132)+p4m
     -   up4*F(3133)
       Fa14(1295) = 2*p2mup5*p4mup4*F(448)+p2mup4*(-4*(p3mup5*F(273)-p
     -   2mup5*F(449))+2*p4mup5*F(485))+mup4mup5*F(3134)-4*(p1mup5*p2
     -   mup4*F(263)-p1mup4*p2mup5*F(376)-p2mup5*p3mup4*F(3135))
       Fa14(1296) = 2*p1mup2*F(3137)+p3mup2*F(3138)+p4mup2*F(3139)
       Fa14(1297) = p3mup2*F(3140)+p1mup2*F(3141)+p4mup2*F(3142)
       Fa14(1298) = p1mup2*F(3143)+p4mup2*F(3144)+2*p3mup2*F(3145)
       Fa14(1299) = p1mup2*F(3146)+p3mup2*F(3147)+2*p4mup2*F(3148)
       Fa14(1300) = p2mup4*F(3149)+p3mup4*F(3150)-(2*p1mup4*F(376)+p4m
     -   up4*F(448))*P(1)
       Fa14(1301) = p2mup4*F(3151)+p3mup4*F(3152)-(p1mup4*F(443)+2*p4m
     -   up4*F(459))*P(1)
       Fa14(1302) = p2mup1*F(3109)+p3mup1*F(3134)-p4mup1*F(477)*P(1)
       Fa14(1303) = p3mup1*F(3110)+p4mup1*F(3113)
       Fa14(1304) = p4mup1*F(460)+p2mup1*F(3155)
       Fa14(1305) = p3mup1*F(480)+p2mup1*F(3156)
       Fa14(1306) = -4*(p1mup5*p3mup4*F(263)+p2mup5*p3mup4*F(379))+2*(
     -   p3mup5*(2*(p1mup4*F(376)+p2mup4*F(377)+p3mup4*F(396))+p4mup4
     -   *F(448))+p3mup4*p4mup5*F(485))+mup4mup5*F(3157)
       Fa14(1307) = 2*(p1mup5*p4mup4*F(454)+p1mup4*(2*(p3mup5*F(273)+p
     -   2mup5*F(379)+p1mup5*F(451))+p4mup5*F(454)))+8*(p3mup4*p3mup5
     -   *F(382)+p4mup4*p4mup5*F(459)+p2mup4*p3mup5*F(1500))+4*(p3mup
     -   5*p4mup4*F(445)+p1mup5*p2mup4*F(453)+p2mup4*p4mup5*F(455)+p2
     -   mup5*(p4mup4*F(444)+2*(p3mup4*F(381)+p2mup4*F(456)))+p1mup5*
     -   p3mup4*F(3111)+p3mup4*p4mup5*F(3112))+mup4mup5*F(3155)
       Fa14(1308) = 2*(p1mup5*p3mup4*F(443)-p1mup4*p3mup5*F(443))+4*((
     -   p2mup5*p3mup4-p2mup4*p3mup5)*F(444)-(p3mup5*p4mup4-p3mup4*p4
     -   mup5)*F(459))+mup4mup5*F(3156)
       Fa14(1309) = p2mup5*F(3159)+p3mup5*F(3160)+(-2*p1mup5*F(263)+p4
     -   mup5*F(485))*P(1)
       Fa14(1310) = p3mup5*F(3161)+p1mup5*F(3162)+p4mup5*F(3163)+2*p2m
     -   up5*F(3164)
       Fa14(1311) = 2*p3mup5*F(3123)-p2mup5*F(3165)+p1mup5*F(3166)+p4m
     -   up5*F(3167)
       Fa14(1312) = p2mup5*F(3151)+p3mup5*F(3152)-(p1mup5*F(443)+2*p4m
     -   up5*F(459))*P(1)
       Fa14(1313) = p3mup4*F(3168)+2*p2mup4*F(3169)+p1mup4*F(3170)+p4m
     -   up4*F(3171)
       Fa14(1314) = -(p2mup4*F(3172))+p1mup4*F(3173)+2*p3mup4*F(3174)+
     -   p4mup4*F(3175)
       Fa14(1315) = 4*p1mup5*p2mup4*F(263)+2*(p2mup4*(p4mup5*F(447)+2*
     -   (p3mup5*F(273)-p2mup5*F(449)))-p2mup5*(p4mup4*F(448)+2*(p1mu
     -   p4*F(376)+p3mup4*F(3135))))+mup4mup5*F(3158)
       Fa14(1316) = 2*p1mup2*F(3177)+p3mup2*F(3178)+p4mup2*F(3179)
       Fa14(1317) = p3mup2*F(3180)+p1mup2*F(3181)+p4mup2*F(3182)
       Fa14(1318) = -(p1mup2*F(3183))+p4mup2*F(3184)-2*p3mup2*F(3185)
       Fa14(1319) = p1mup2*F(3186)+p3mup2*F(3187)+2*p4mup2*F(3188)
       Fa14(1320) = p2mup4*F(3189)+p3mup4*F(3190)+(2*p1mup4*F(376)+p4m
     -   up4*F(448))*P(1)
       Fa14(1321) = p2mup4*F(3124)+p3mup4*F(3125)+(p1mup4*F(443)+2*p4m
     -   up4*F(459))*P(1)
       Fa14(1322) = p2mup1*F(3157)+p3mup1*F(3158)+p4mup1*F(477)*P(1)
       Fa14(1323) = p3mup1*F(3155)+p4mup1*F(3156)
       Fa14(1324) = p4mup1*F(480)+p2mup1*F(3110)
       Fa14(1325) = p3mup1*F(460)+p2mup1*F(3113)
       Fa14(1326) = 4*(p1mup5*p3mup4*F(263)-p1mup4*p3mup5*F(263)+(p2mu
     -   p5*p3mup4-p2mup4*p3mup5)*F(379))-2*(p3mup5*p4mup4-p3mup4*p4m
     -   up5)*F(383)+mup4mup5*F(3107)
       Fa14(1327) = mup4mup5*F(3193)-2*(p1mup5*(2*(p3mup4*F(273)+p2mup
     -   4*F(379))+p4mup4*F(4464))-p1mup4*(2*(p3mup5*F(273)+p2mup5*F(
     -   379))+p4mup5*F(4464)))+4*((p3mup5*p4mup4-p3mup4*p4mup5)*F(44
     -   65)+(p2mup5*p4mup4-p2mup4*p4mup5)*F(4469))
       Fa14(1328) = -2*(p1mup5*p3mup4*F(383)-p1mup4*p3mup5*F(383))+mup
     -   4mup5*F(3106)+4*((p3mup5*p4mup4-p3mup4*p4mup5)*F(4467)-(p2mu
     -   p5*p3mup4-p2mup4*p3mup5)*F(4469))
       Fa14(1329) = p2mup5*F(3099)+p3mup5*F(3100)+(2*p1mup5*F(263)+p4m
     -   up5*F(383))*P(1)
       Fa14(1330) = p1mup5*F(3078)+p3mup5*F(3079)+2*p2mup5*F(3080)+p4m
     -   up5*F(3081)
       Fa14(1331) = -2*p3mup5*F(3075)+p1mup5*F(3082)+p2mup5*F(3083)+p4
     -   mup5*F(3084)
       Fa14(1332) = p3mup5*F(3056)+p2mup5*F(3101)+(p1mup5*F(383)+2*p4m
     -   up5*F(4467))*P(1)
       Fa14(1333) = p1mup4*F(3069)+p3mup4*F(3070)+2*p2mup4*F(3071)+p4m
     -   up4*F(3072)
       Fa14(1334) = -(p1mup4*F(3073))+p2mup4*F(3074)+2*p3mup4*F(3075)+
     -   p4mup4*F(3076)
       Fa14(1335) = -4*(p1mup5*p2mup4*F(263)-p1mup4*p2mup5*F(263)-(p2m
     -   up5*p3mup4-p2mup4*p3mup5)*F(273))+2*(p2mup5*p4mup4-p2mup4*p4
     -   mup5)*F(383)+mup4mup5*F(3108)
       Fa14(1336) = 2*p1mup2*F(3195)+p3mup2*F(3196)+p4mup2*F(3197)
       Fa14(1337) = p1mup2*F(3198)+p3mup2*F(3199)+p4mup2*F(3200)
       Fa14(1338) = p1mup2*F(3201)+2*p3mup2*F(3202)+p4mup2*F(3203)
       Fa14(1339) = p1mup2*F(3204)+p3mup2*F(3205)+2*p4mup2*F(3206)
       Fa14(1340) = p2mup4*F(3067)+p3mup4*F(3068)-(2*p1mup4*F(263)+p4m
     -   up4*F(383))*P(1)
       Fa14(1341) = p3mup4*F(3032)+p2mup4*F(3077)-(p1mup4*F(383)+2*p4m
     -   up4*F(4467))*P(1)
       Fa14(1342) = p2mup1*F(3107)+p3mup1*F(3108)-p4mup1*F(4464)*P(1)
       Fa14(1343) = p3mup1*F(3104)+p4mup1*F(3106)
       Fa14(1344) = p4mup1*F(413)+p2mup1*F(3105)
       Fa14(1345) = p3mup1*F(390)+p2mup1*F(3066)
       Fa14(1346) = 4*(p1mup5*p3mup4*F(263)-p1mup4*p3mup5*F(263)+(p2mu
     -   p5*p3mup4-p2mup4*p3mup5)*F(379))-2*(p3mup5*p4mup4-p3mup4*p4m
     -   up5)*F(383)+mup4mup5*F(3062)
       Fa14(1347) = mup4mup5*F(3209)-2*(p1mup5*(2*(p3mup4*F(273)+p2mup
     -   4*F(379))+p4mup4*F(4464))-p1mup4*(2*(p3mup5*F(273)+p2mup5*F(
     -   379))+p4mup5*F(4464)))+4*((p3mup5*p4mup4-p3mup4*p4mup5)*F(44
     -   65)+(p2mup5*p4mup4-p2mup4*p4mup5)*F(4469))
       Fa14(1348) = -2*(p1mup5*p3mup4*F(383)-p1mup4*p3mup5*F(383))+mup
     -   4mup5*F(3061)+4*((p3mup5*p4mup4-p3mup4*p4mup5)*F(4467)-(p2mu
     -   p5*p3mup4-p2mup4*p3mup5)*F(4469))
       Fa14(1349) = p2mup5*F(3053)+p3mup5*F(3054)+(2*p1mup5*F(263)+p4m
     -   up5*F(383))*P(1)
       Fa14(1350) = p1mup5*F(3033)+p3mup5*F(3034)+2*p2mup5*F(3035)+p4m
     -   up5*F(3036)
       Fa14(1351) = 2*p3mup5*F(3029)+p4mup5*F(3030)+p1mup5*F(3037)+p2m
     -   up5*F(3038)
       Fa14(1352) = p2mup5*F(3055)+p3mup5*F(3056)+(p1mup5*F(383)+2*p4m
     -   up5*F(4467))*P(1)
       Fa14(1353) = p1mup4*F(3023)+p3mup4*F(3024)+2*p2mup4*F(3025)+p4m
     -   up4*F(3026)
       Fa14(1354) = p1mup4*F(3027)-p2mup4*F(3028)+2*p3mup4*F(3029)+p4m
     -   up4*F(3030)
       Fa14(1355) = -4*(p1mup5*p2mup4*F(263)-p1mup4*p2mup5*F(263)-(p2m
     -   up5*p3mup4-p2mup4*p3mup5)*F(273))+2*(p2mup5*p4mup4-p2mup4*p4
     -   mup5)*F(383)+mup4mup5*F(3063)
       Fa14(1356) = 2*p1mup2*F(3211)+p3mup2*F(3212)+p4mup2*F(3213)
       Fa14(1357) = p1mup2*F(3214)+p3mup2*F(3215)+p4mup2*F(3216)
       Fa14(1358) = p3mup2*F(3218)+M1*M4*(p1mup2*F(3217)+p4mup2*F(3219
     -   ))
       Fa14(1359) = p1mup2*F(3220)+p3mup2*F(3221)+2*p4mup2*F(3222)
       Fa14(1360) = p2mup4*F(3021)+p3mup4*F(3022)-(2*p1mup4*F(263)+p4m
     -   up4*F(383))*P(1)
       Fa14(1361) = p2mup4*F(3031)+p3mup4*F(3032)-(p1mup4*F(383)+2*p4m
     -   up4*F(4467))*P(1)
       Fa14(1362) = p2mup1*F(3062)+p3mup1*F(3063)-p4mup1*F(4464)*P(1)
       Fa14(1363) = p3mup1*F(3059)+p4mup1*F(3061)
       Fa14(1364) = p4mup1*F(413)+p2mup1*F(3060)
       Fa14(1365) = p3mup1*F(390)+p2mup1*F(3020)
       Fa14(1366) = 8*p2mup4*p2mup5*F(2892)+mup4mup5*F(3225)+p1mup5*(8
     -   *p2mup4*F(2885)+p3mup4*F(3226)+4*p4mup4*F(3227))+p2mup5*p3mu
     -   p4*F(3229)+p3mup4*p4mup5*F(3230)+2*p3mup5*p4mup4*F(3231)+4*(
     -   p1mup4*(p3mup5*F(2884)+p4mup5*F(2887)+2*(p1mup5*F(2883)+p2mu
     -   p5*F(2891)))+p2mup4*p3mup5*F(2893)+p3mup4*p3mup5*F(2897)+p2m
     -   up5*p4mup4*F(2917)+p2mup4*p4mup5*F(3228)+p4mup4*p4mup5*F(323
     -   2))
       Fa14(1367) = 8*p4mup4*p4mup5*F(2912)+4*p2mup4*(p4mup5*F(2905)+2
     -   *(p3mup5*F(2910)+p2mup5*F(2914)))+mup4mup5*F(3233)+p1mup4*p2
     -   mup5*F(3234)+p1mup4*p3mup5*F(3235)+p1mup4*p4mup5*F(3236)+p3m
     -   up4*p4mup5*F(3237)+p3mup5*p4mup4*F(3238)+p1mup5*(4*(p2mup4*F
     -   (2901)+p1mup4*F(2902))+2*p4mup4*F(2909)+p3mup4*F(3239))+p3mu
     -   p4*p3mup5*F(3240)+p2mup5*p4mup4*F(3241)+p2mup5*p3mup4*F(3242
     -   )
       Fa14(1368) = -4*p2mup4*p3mup5*F(2927)+mup4mup5*F(3243)+p1mup4*(
     -   -2*p3mup5*F(2919)+4*(p2mup5*F(2890)+p1mup5*F(3245)+p4mup5*F(
     -   3246)))+p3mup4*p4mup5*F(3248)+p3mup5*p4mup4*F(3249)+p1mup5*(
     -   p3mup4*F(3244)+4*(p2mup4*F(2918)+p4mup4*F(3250)))+8*(p2mup4*
     -   p4mup5*F(3247)+p3mup4*p3mup5*F(3251)+p4mup4*p4mup5*F(3252)+p
     -   2mup5*p4mup4*F(3253)+p2mup4*p2mup5*F(3254))+p2mup5*p3mup4*F(
     -   3255)
       Fa14(1369) = -2*p1mup5*F(2936)+p2mup5*F(3256)+p3mup5*F(3257)+p4
     -   mup5*F(3258)
       Fa14(1370) = -2*p2mup5*F(2940)+p1mup5*F(3259)+p3mup5*F(3260)+p4
     -   mup5*F(3261)
       Fa14(1371) = 2*p3mup5*F(2944)+p1mup5*F(3262)+p2mup5*F(3263)+p4m
     -   up5*F(3264)
       Fa14(1372) = -2*p4mup5*F(2948)+p1mup5*F(3265)+p2mup5*F(3266)+p3
     -   mup5*F(3267)
       Fa14(1373) = p1mup4*F(3268)+p2mup4*F(3269)+p3mup4*F(3270)+p4mup
     -   4*F(3271)
       Fa14(1374) = p1mup4*F(3272)+p2mup4*F(3273)+p3mup4*F(3274)+p4mup
     -   4*F(3275)
       Fa14(1375) = -2*(4*p3mup4*p3mup5*F(2438)+p2mup4*p4mup5*F(2960)+
     -   p2mup5*p4mup4*F(2961)+2*(p3mup4*p4mup5*F(2435)+p3mup5*p4mup4
     -   *F(2436)+p4mup4*p4mup5*F(2439)+p1mup4*(p4mup5*F(2433)+2*(p1m
     -   up5*F(2431)+p3mup5*F(2437))+p2mup5*F(2958))+p1mup5*(2*p3mup4
     -   *F(2432)+p4mup4*F(2434)-p2mup4*F(2959))+p2mup4*p2mup5*F(2962
     -   )-p2mup4*p3mup5*F(2963)+p2mup5*p3mup4*F(2964)))+mup4mup5*F(3
     -   276)
       Fa14(1376) = 8*(p3mup4*p4mup5*F(2460)+p4mup4*p4mup5*F(2462)+p3m
     -   up5*p4mup4*F(2464)+p2mup4*p2mup5*F(2465)+p3mup4*p3mup5*F(246
     -   5))+2*p1mup4*(2*(p3mup5*F(2436)+p1mup5*F(2456)+p4mup5*F(2457
     -   ))+p2mup5*F(2966))+p1mup5*(4*(p3mup4*F(2455)+p4mup4*F(2461))
     -   -2*p2mup4*F(2967))-4*(p2mup4*p4mup5*F(2968)-p2mup5*p4mup4*F(
     -   2969)+p2mup4*p3mup5*F(2970)-p2mup5*p3mup4*F(2971))+mup4mup5*
     -   F(3277)
       Fa14(1377) = p1mup2*F(3279)+p3mup2*F(3280)+p4mup2*F(3281)
       Fa14(1378) = p1mup2*F(3282)+p3mup2*F(3283)+p4mup2*F(3284)
       Fa14(1379) = p1mup2*F(3285)+p4mup2*F(3286)+p3mup2*F(3287)
       Fa14(1380) = p1mup2*F(3288)+p3mup2*F(3289)+p4mup2*F(3290)
       Fa14(1381) = p2mup4*F(3291)+p3mup4*F(3292)+p1mup4*F(3293)+p4mup
     -   4*F(3294)
       Fa14(1382) = p2mup4*F(3295)+p1mup4*F(3296)+p3mup4*F(3297)-2*p4m
     -   up4*F(3298)
       Fa14(1383) = mup4mup5*F(3017)+4*(p1mup4*p3mup5*F(2525)-p1mup5*p
     -   3mup4*F(2527)-p3mup4*p4mup5*F(2531)-p3mup4*p3mup5*F(2533)-p3
     -   mup5*p4mup4*F(2535)-p1mup4*p2mup5*F(3001)+p2mup5*p4mup4*F(30
     -   07)+p1mup5*p2mup4*F(3305)+p2mup4*p4mup5*F(3306)+p2mup4*p2mup
     -   5*F(3307))+2*p2mup4*p3mup5*F(3308)+p2mup5*p3mup4*F(3309)
       Fa14(1384) = p2mup1*F(3015)+p3mup1*F(3016)+p4mup1*F(3017)
       Fa14(1385) = p3mup1*F(2999)+p4mup1*F(3013)
       Fa14(1386) = p2mup1*F(3010)+p4mup1*F(3014)
       Return
       End
