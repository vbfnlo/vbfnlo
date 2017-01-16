c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenOddGSim15(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4
     -   ,p3mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5m
     -   up3,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mu
     -   p3,mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa15)
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
       Complex*16 Fa15(1387:1485)
       COMMON/FaPenOddGSim/Fa
       Fa15(1387) = p2mup1*F(3011)+p3mup1*F(3012)
       Fa15(1388) = mup2mup5*F(3015)+2*p3mup5*p4mup2*F(3315)+p3mup2*p4
     -   mup5*F(3317)+4*(p1mup2*p3mup5*F(3310)+p1mup2*p4mup5*F(3311)+
     -   p1mup5*(p4mup2*F(3305)+p3mup2*F(3312))+p2mup5*p3mup2*F(3313)
     -   +p2mup5*p4mup2*F(3314)+p3mup2*p3mup5*F(3316)+p4mup2*p4mup5*F
     -   (3318))
       Fa15(1389) = mup2mup5*F(3319)-2*(p1mup2*(2*(p1mup5*F(3320)+p2mu
     -   p5*F(3321)+p3mup5*F(3322))+p4mup5*F(3323))+2*p3mup2*(-(p4mup
     -   5*F(3325))+p1mup5*F(3327)+2*(p3mup5*F(3324)+p2mup5*F(3329)))
     -   -p4mup2*(-4*p4mup5*F(3326)+p1mup5*F(3328)+2*(p3mup5*F(3325)+
     -   p2mup5*F(3330))))
       Fa15(1390) = mup2mup5*F(3331)+p1mup2*(-2*p3mup5*F(3334)+4*(p1mu
     -   p5*F(3332)+p2mup5*F(3333)+p4mup5*F(3335)))+p3mup2*(p2mup5*F(
     -   3336)-2*p1mup5*F(3337)-8*p3mup5*F(3338)+p4mup5*F(3339))+p4mu
     -   p2*(8*(p4mup5*F(3247)+p2mup5*F(3254))+p3mup5*F(3339)+4*p1mup
     -   5*F(3340))
       Fa15(1391) = mup2mup5*F(3341)-2*p2mup5*(2*(p3mup2*F(3327)+p1mup
     -   2*F(3343))+p4mup2*F(3346))-4*(-(p1mup2*p4mup5*F(3344))+p1mup
     -   5*p4mup2*F(3345)+p3mup5*p4mup2*F(3347)+p4mup2*p4mup5*F(3349)
     -   -p3mup2*p4mup5*F(3350)+2*(p1mup2*p1mup5*F(3342)+p1mup2*p3mup
     -   5*F(3348)+p3mup2*p3mup5*F(3351)+p1mup5*p3mup2*F(3352)))
       Fa15(1392) = mup2mup5*F(3353)+2*p2mup5*(2*p4mup2*F(3340)+p3mup2
     -   *F(3357))-4*(-(p1mup5*p3mup2*F(3356))-p3mup2*p3mup5*F(3358)+
     -   p1mup2*(p2mup5*F(3354)-p3mup5*F(3355)+2*(p1mup5*F(2883)+p4mu
     -   p5*F(3359)))-p3mup2*p4mup5*F(3360)-p3mup5*p4mup2*F(3361)-2*(
     -   p4mup2*p4mup5*F(3362)+p1mup5*p4mup2*F(3363)))
       Fa15(1393) = mup2mup5*F(3364)+2*(p1mup2*(p2mup5*F(3366)+2*(p1mu
     -   p5*F(3365)+p3mup5*F(3367)+p4mup5*F(3368)))+2*(p2mup5*p3mup2*
     -   F(3369)+p1mup5*p3mup2*F(3370)+p1mup5*p4mup2*F(3373)+p2mup5*p
     -   4mup2*F(3374)+2*(p3mup2*p3mup5*F(3371)+p4mup2*p4mup5*F(3372)
     -   +(p3mup5*p4mup2+p3mup2*p4mup5)*F(3375))))
       Fa15(1394) = 4*(p1mup5*p3mup2*F(263)+p2mup5*p3mup2*F(379)+p3mup
     -   2*p3mup5*F(1268)-p1mup2*p3mup5*F(1608))+2*(p3mup2*p4mup5*F(3
     -   83)+p3mup5*p4mup2*F(1609))-mup2mup5*F(3062)
       Fa15(1395) = -2*(p1mup5*p4mup2*F(1002)-p1mup2*(p4mup5*F(1001)+2
     -   *(p1mup5*F(678)+p2mup5*F(680)+p3mup5*F(1292))))+4*(p2mup5*p4
     -   mup2*F(1004)+p3mup5*p4mup2*F(1010)+p3mup2*(p4mup5*F(1010)+p1
     -   mup5*F(1302)+2*(p3mup5*F(311)+p2mup5*F(1303))))+mup2mup5*F(3
     -   224)+8*p4mup2*p4mup5*F(4472)
       Fa15(1396) = 2*(p1mup5*p3mup2*F(383)+p1mup2*p3mup5*F(383))+mup2
     -   mup5*F(3061)+8*p3mup2*p3mup5*F(4465)+4*(p3mup5*p4mup2*F(4467
     -   )+p3mup2*p4mup5*F(4467)+p2mup5*p3mup2*F(4469))
       Fa15(1397) = 2*p2mup5*p4mup2*F(993)+8*(p3mup5*p4mup2*F(945)+p1m
     -   up5*(p4mup2*F(364)+p3mup2*F(947))+p3mup2*p3mup5*F(1025))+4*(
     -   p1mup2*(2*(p3mup5*F(266)+p1mup5*F(634))+p4mup5*F(684)+p2mup5
     -   *F(940))+p3mup2*p4mup5*F(995)+p4mup2*p4mup5*F(996)+p2mup5*p3
     -   mup2*F(1302))+mup2mup5*F(3223)
       Fa15(1398) = 4*p3mup5*(p1mup2*F(357)+p3mup2*F(4465)+p4mup2*F(44
     -   67))-F(4464)*(2*p2mup5*p3mup2+mup2mup5*P(1))
       Fa15(1399) = 4*(p1mup5*p3mup2*F(263)+p2mup5*p3mup2*F(379)+p3mup
     -   2*p3mup5*F(1268)-p1mup2*p3mup5*F(1608))+2*(p3mup2*p4mup5*F(3
     -   83)+p3mup5*p4mup2*F(1609))-mup2mup5*F(3107)
       Fa15(1400) = -2*(p1mup5*p4mup2*F(1002)-p1mup2*(p4mup5*F(1001)+2
     -   *(p1mup5*F(678)+p2mup5*F(680)+p3mup5*F(1292))))+4*(p2mup5*p4
     -   mup2*F(1004)+p3mup5*p4mup2*F(1010)+p3mup2*(p4mup5*F(1010)+p1
     -   mup5*F(1302)+2*(p3mup5*F(311)+p2mup5*F(1303))))+mup2mup5*F(3
     -   208)+8*p4mup2*p4mup5*F(4472)
       Fa15(1401) = 2*(p1mup5*p3mup2*F(383)+p1mup2*p3mup5*F(383))+mup2
     -   mup5*F(3106)+8*p3mup2*p3mup5*F(4465)+4*(p3mup5*p4mup2*F(4467
     -   )+p3mup2*p4mup5*F(4467)+p2mup5*p3mup2*F(4469))
       Fa15(1402) = 2*p2mup5*p4mup2*F(993)+8*(p3mup5*p4mup2*F(945)+p1m
     -   up5*(p4mup2*F(364)+p3mup2*F(947))+p3mup2*p3mup5*F(1025))+4*(
     -   p1mup2*(2*(p3mup5*F(266)+p1mup5*F(634))+p4mup5*F(684)+p2mup5
     -   *F(940))+p3mup2*p4mup5*F(995)+p4mup2*p4mup5*F(996)+p2mup5*p3
     -   mup2*F(1302))+mup2mup5*F(3207)
       Fa15(1403) = 4*p1mup2*p3mup5*F(1262)+2*(p3mup2*(p4mup5*F(447)+2
     -   *(p1mup5*F(263)+p2mup5*F(379)+p3mup5*F(1268)))+p3mup5*p4mup2
     -   *F(1623))+mup2mup5*F(3157)
       Fa15(1404) = p1mup2*(-4*(p1mup5*F(678)+p2mup5*F(680))+2*p4mup5*
     -   F(1035))+2*p1mup5*p4mup2*F(1036)+8*p4mup2*p4mup5*F(1041)+mup
     -   2mup5*F(3376)+4*(p2mup5*p4mup2*F(1038)-p1mup2*p3mup5*F(1292)
     -   -p3mup2*(p1mup5*F(1302)+2*(p3mup5*F(311)+p2mup5*F(1303)))+(p
     -   3mup5*p4mup2+p3mup2*p4mup5)*F(3377))
       Fa15(1405) = 2*(p1mup5*p3mup2*F(443)+p1mup2*p3mup5*F(443))+8*p3
     -   mup2*p3mup5*F(445)+4*(p2mup5*p3mup2*F(444)+(p3mup5*p4mup2+p3
     -   mup2*p4mup5)*F(459))-mup2mup5*F(3156)
       Fa15(1406) = 8*(p3mup5*p4mup2*F(945)+p1mup5*(p4mup2*F(364)+p3mu
     -   p2*F(947))+p3mup2*p3mup5*F(1025))+2*p2mup5*p4mup2*F(1027)+4*
     -   (p1mup2*(2*(p3mup5*F(266)+p1mup5*F(634))+p2mup5*F(940)-p4mup
     -   5*F(1028))-p3mup2*p4mup5*F(1030)-p4mup2*p4mup5*F(1032)+p2mup
     -   5*p3mup2*F(1302))-mup2mup5*F(3378)
       Fa15(1407) = p3mup2*(-4*p3mup5*F(445)+2*p2mup5*F(477))-4*p3mup5
     -   *(p4mup2*F(459)+p1mup2*F(1629))+mup2mup5*F(477)*P(1)
       Fa15(1408) = -4*(p1mup5*p3mup2*F(263)+p2mup5*p3mup2*F(379)+p3mu
     -   p2*p3mup5*F(1268)-p1mup2*p3mup5*F(1608))+2*(p3mup2*p4mup5*F(
     -   485)-p3mup5*p4mup2*F(1623))+mup2mup5*F(3109)
       Fa15(1409) = -2*(p1mup5*p4mup2*F(1036)-p1mup2*p4mup5*F(1062))+8
     -   *(p3mup2*p3mup5*F(311)+p4mup2*p4mup5*F(1064)+p2mup5*p3mup2*F
     -   (1303))+mup2mup5*F(3379)+4*(p1mup2*p1mup5*F(678)+p1mup2*p2mu
     -   p5*F(680)+p2mup5*p4mup2*F(1063)+p1mup2*p3mup5*F(1292)+p1mup5
     -   *p3mup2*F(1302)+(p3mup5*p4mup2+p3mup2*p4mup5)*F(3380))
       Fa15(1410) = 2*(p1mup5*p3mup2*F(443)+p1mup2*p3mup5*F(443))+8*p3
     -   mup2*p3mup5*F(445)+4*(p2mup5*p3mup2*F(444)+(p3mup5*p4mup2+p3
     -   mup2*p4mup5)*F(459))+mup2mup5*F(3113)
       Fa15(1411) = 8*(p3mup5*p4mup2*F(945)+p1mup5*(p4mup2*F(364)+p3mu
     -   p2*F(947))+p3mup2*p3mup5*F(1025))+2*p2mup5*p4mup2*F(1027)+4*
     -   (p1mup2*(2*(p3mup5*F(266)+p1mup5*F(634))+p2mup5*F(940)-p4mup
     -   5*F(1028))+p3mup2*p4mup5*F(1059)+p4mup2*p4mup5*F(1060)+p2mup
     -   5*p3mup2*F(1302))+mup2mup5*F(3381)
       Fa15(1412) = -2*p2mup5*p3mup2*F(477)+4*(p3mup2*p3mup5*F(445)+p3
     -   mup5*p4mup2*F(459)+p1mup2*p3mup5*F(1629))-mup2mup5*F(477)*P(
     -   1)
       Fa15(1413) = 4*p1mup2*p3mup5*F(1262)+2*(p3mup2*(p4mup5*F(383)+2
     -   *(p1mup5*F(263)+p2mup5*F(379)+p3mup5*F(1268)))+p3mup5*p4mup2
     -   *F(1609))+mup2mup5*F(3064)
       Fa15(1414) = -2*(p1mup5*p4mup2*F(1002)-p1mup2*(p4mup5*F(1001)+2
     -   *(p1mup5*F(678)+p2mup5*F(680)+p3mup5*F(1292))))+4*(p2mup5*p4
     -   mup2*F(1004)+p3mup5*p4mup2*F(1010)+p3mup2*(p4mup5*F(1010)+p1
     -   mup5*F(1302)+2*(p3mup5*F(311)+p2mup5*F(1303))))-mup2mup5*F(3
     -   103)+8*p4mup2*p4mup5*F(4472)
       Fa15(1415) = 2*(p1mup5*p3mup2*F(383)+p1mup2*p3mup5*F(383))-mup2
     -   mup5*F(3066)+8*p3mup2*p3mup5*F(4465)+4*(p3mup5*p4mup2*F(4467
     -   )+p3mup2*p4mup5*F(4467)+p2mup5*p3mup2*F(4469))
       Fa15(1416) = 2*p2mup5*p4mup2*F(993)+8*(p3mup5*p4mup2*F(945)+p1m
     -   up5*(p4mup2*F(364)+p3mup2*F(947))+p3mup2*p3mup5*F(1025))+4*(
     -   p1mup2*(2*(p3mup5*F(266)+p1mup5*F(634))+p4mup5*F(684)+p2mup5
     -   *F(940))+p3mup2*p4mup5*F(995)+p4mup2*p4mup5*F(996)+p2mup5*p3
     -   mup2*F(1302))-mup2mup5*F(3102)
       Fa15(1417) = -4*p3mup5*(p1mup2*F(357)+p3mup2*F(4465)+p4mup2*F(4
     -   467))+F(4464)*(2*p2mup5*p3mup2+mup2mup5*P(1))
       Fa15(1418) = 4*p1mup2*p3mup5*F(1262)+2*(p3mup2*(p4mup5*F(383)+2
     -   *(p1mup5*F(263)+p2mup5*F(379)+p3mup5*F(1268)))+p3mup5*p4mup2
     -   *F(1609))+mup2mup5*F(3018)
       Fa15(1419) = -2*(p1mup5*p4mup2*F(1002)-p1mup2*(p4mup5*F(1001)+2
     -   *(p1mup5*F(678)+p2mup5*F(680)+p3mup5*F(1292))))+4*(p2mup5*p4
     -   mup2*F(1004)+p3mup5*p4mup2*F(1010)+p3mup2*(p4mup5*F(1010)+p1
     -   mup5*F(1302)+2*(p3mup5*F(311)+p2mup5*F(1303))))-mup2mup5*F(3
     -   058)+8*p4mup2*p4mup5*F(4472)
       Fa15(1420) = 2*(p1mup5*p3mup2*F(383)+p1mup2*p3mup5*F(383))-mup2
     -   mup5*F(3020)+8*p3mup2*p3mup5*F(4465)+4*(p3mup5*p4mup2*F(4467
     -   )+p3mup2*p4mup5*F(4467)+p2mup5*p3mup2*F(4469))
       Fa15(1421) = 2*p2mup5*p4mup2*F(993)+8*(p3mup5*p4mup2*F(945)+p1m
     -   up5*(p4mup2*F(364)+p3mup2*F(947))+p3mup2*p3mup5*F(1025))+4*(
     -   p1mup2*(2*(p3mup5*F(266)+p1mup5*F(634))+p4mup5*F(684)+p2mup5
     -   *F(940))+p3mup2*p4mup5*F(995)+p4mup2*p4mup5*F(996)+p2mup5*p3
     -   mup2*F(1302))-mup2mup5*F(3057)
       Fa15(1422) = mup2mup5*F(3008)+4*(p1mup5*p4mup2*F(3002)-p1mup2*(
     -   p3mup5*F(3310)+p4mup5*F(3311))-p1mup5*p3mup2*F(3312)-p2mup5*
     -   p3mup2*F(3313)-p3mup2*p3mup5*F(3316)-p4mup2*p4mup5*F(3318)+p
     -   2mup5*p4mup2*F(3382))+p3mup5*p4mup2*F(3383)+2*p3mup2*p4mup5*
     -   F(3384)
       Fa15(1423) = mup2mup5*F(3385)+p1mup2*(4*(p1mup5*F(3320)+p2mup5*
     -   F(3321)+p3mup5*F(3322))+p4mup5*F(3386))+p3mup2*(4*p1mup5*F(3
     -   327)+8*(p3mup5*F(3324)+p2mup5*F(3329))+p4mup5*F(3387))+p4mup
     -   2*(8*p4mup5*F(3326)-2*p1mup5*F(3328)+p3mup5*F(3387)+p2mup5*F
     -   (3388))
       Fa15(1424) = mup2mup5*F(3389)+2*(p1mup2*(p3mup5*F(3334)+2*(p1mu
     -   p5*F(3390)+p2mup5*F(3391)+p4mup5*F(3392)))+p3mup2*(p1mup5*F(
     -   3337)+4*p3mup5*F(3338)+2*(p2mup5*F(3393)+p4mup5*F(3394))))+4
     -   *p4mup2*(2*(p4mup5*F(2923)+p2mup5*F(2931))+p3mup5*F(3394)+p1
     -   mup5*F(3395))
       Fa15(1425) = 2*p2mup5*(2*(p3mup2*F(3327)+p1mup2*F(3343))+p4mup2
     -   *F(3346))+p1mup2*(-4*p4mup5*F(3344)+8*(p1mup5*F(3342)+p3mup5
     -   *F(3348)))+mup2mup5*F(3396)+4*(p1mup5*p4mup2*F(3345)+p3mup5*
     -   p4mup2*F(3347)+p4mup2*p4mup5*F(3349)+2*(p3mup2*p3mup5*F(3351
     -   )+p1mup5*p3mup2*F(3352))+p3mup2*p4mup5*F(3397))
       Fa15(1426) = -2*p2mup5*p3mup2*F(3357)+mup2mup5*F(3398)-4*(p1mup
     -   5*p3mup2*F(3356)+p3mup2*p3mup5*F(3358)-p1mup2*(p2mup5*F(3354
     -   )-p3mup5*F(3355)+2*(p1mup5*F(2883)+p4mup5*F(3359)))+p3mup2*p
     -   4mup5*F(3360)-p2mup5*p4mup2*F(3395)-p3mup5*p4mup2*F(3399))+8
     -   *(p4mup2*p4mup5*F(3400)+p1mup5*p4mup2*F(3401))
       Fa15(1427) = -2*p1mup2*(p2mup5*F(3366)+2*(p1mup5*F(3365)+p3mup5
     -   *F(3367)+p4mup5*F(3368)))-4*(p2mup5*p3mup2*F(3369)+p1mup5*p3
     -   mup2*F(3370)+2*(p3mup2*p3mup5*F(3371)+p4mup2*p4mup5*F(3372))
     -   +p1mup5*p4mup2*F(3373)+p2mup5*p4mup2*F(3374))-8*(p3mup5*p4mu
     -   p2+p3mup2*p4mup5)*F(3375)+mup2mup5*F(3402)
       Fa15(1428) = p1mup2*(8*(p1mup4*F(3403)+p2mup4*F(3404))+p3mup4*F
     -   (3405)+4*p4mup4*F(3407))+p3mup4*p4mup2*F(3411)+p3mup2*p3mup4
     -   *F(3412)-2*p3mup2*p4mup4*F(3413)-4*(p1mup4*(p3mup2*F(3406)+p
     -   4mup2*F(3408))+p2mup4*p3mup2*F(3409)+p2mup4*p4mup2*F(3410)-p
     -   4mup2*p4mup4*F(3414))
       Fa15(1429) = p1mup2*(p2mup4*F(3416)+8*(p1mup4*F(3415)+p3mup4*F(
     -   3417))+p4mup4*F(3418))+2*p2mup4*p4mup2*F(3420)+8*(p3mup2*p3m
     -   up4*F(3422)+p1mup4*p3mup2*F(3423))+4*(p1mup4*p4mup2*F(3419)+
     -   p3mup4*p4mup2*F(3421)+p2mup4*p3mup2*F(3424)+p3mup2*p4mup4*F(
     -   3425)+p4mup2*p4mup4*F(3426))
       Fa15(1430) = p3mup2*p3mup4*F(3431)-2*p2mup4*p3mup2*F(3432)+p1mu
     -   p2*(4*p2mup4*F(3428)+p3mup4*F(3429)+8*(p1mup4*F(3427)-p4mup4
     -   *F(3433)))-4*(p1mup4*(p3mup2*F(3430)+2*p4mup2*F(3433))+p2mup
     -   4*p4mup2*F(3434)-p3mup2*p4mup4*F(3435))+p3mup4*p4mup2*F(3436
     -   )-8*p4mup2*p4mup4*F(3437)
       Fa15(1431) = p1mup2*(4*(p1mup4*F(3438)+p2mup4*F(3439))+p3mup4*F
     -   (3440)+2*p4mup4*F(3442))+p1mup4*p4mup2*F(3443)+p2mup4*p4mup2
     -   *F(3445)+8*(p2mup4*p3mup2*F(3444)+p3mup2*p3mup4*F(3446)+p4mu
     -   p2*p4mup4*F(3447))+4*(p1mup4*p3mup2*F(3441)+p3mup2*p4mup4*F(
     -   3448))+p3mup4*p4mup2*F(3449)
       Fa15(1432) = p1mup2*(p3mup4*F(3452)+4*(p1mup4*F(3450)+p2mup4*F(
     -   3451)+p4mup4*F(3456)))+p4mup2*(-4*p1mup4*F(3457)+p3mup4*F(34
     -   59)-8*(p2mup4*F(3458)+p4mup4*F(3460)))+2*p3mup2*(4*p3mup4*F(
     -   3453)+p1mup4*F(3454)+2*(p2mup4*F(3455)+p4mup4*F(3461)))
       Fa15(1433) = p1mup2*(p2mup4*F(3463)+2*(p1mup4*F(3462)+p3mup4*F(
     -   3464)+p4mup4*F(3468)))+2*(p4mup2*(p1mup4*F(3469)+p2mup4*F(34
     -   70)+2*(p3mup4*F(3471)+p4mup4*F(3472)))+p3mup2*(p1mup4*F(3466
     -   )+p2mup4*F(3467)+2*(p3mup4*F(3465)+p4mup4*F(3473))))
       Fa15(1434) = -2*mup2mup4*(p1mup5*F(3312)+p2mup5*F(3313)+p3mup5*
     -   F(3476))+mup2mup4*p4mup5*F(3477)+mup2mup5*(2*(p1mup4*F(3474)
     -   +p2mup4*F(3475)+p3mup4*F(3478))+p4mup4*F(3479))-mup4mup5*(2*
     -   (p1mup2*F(3310)-p3mup2*F(3480))+p4mup2*F(3481))
       Fa15(1435) = mup4mup5*p1mup2*F(3482)-mup2mup5*p1mup4*F(3482)+mu
     -   p2mup4*p1mup5*F(3482)-2*(mup2mup5*p2mup4*F(3483)-mup2mup4*p2
     -   mup5*F(3483)-mup4mup5*p3mup2*F(3484)+mup2mup5*p3mup4*F(3484)
     -   -mup2mup4*p3mup5*F(3484)-(mup4mup5*p4mup2-mup2mup5*p4mup4+mu
     -   p2mup4*p4mup5)*F(3485))
       Fa15(1436) = mup2mup5*(4*(p1mup4*F(3001)-p4mup4*F(3007)+p2mup4*
     -   F(3486))+p3mup4*F(3487))+mup2mup4*(4*(p1mup5*F(3002)+p4mup5*
     -   F(3003)+p2mup5*F(3382))+p3mup5*F(3488))+2*mup4mup5*(p3mup2*F
     -   (3489)-2*(p1mup2*F(3311)+p4mup2*F(3490)))
       Fa15(1437) = -2*(mup2mup5*(2*(p1mup4*F(2525)-p4mup4*F(2535)+p3m
     -   up4*F(2565))+p2mup4*F(3492))-mup2mup4*(2*(p1mup5*F(2527)+p4m
     -   up5*F(2531)+p3mup5*F(2566))+p2mup5*F(3493)))-4*mup4mup5*(p1m
     -   up2*F(3491)+p3mup2*F(3494)+p4mup2*F(3495))
       Fa15(1438) = 2*p1mup2*p3mup4*F(1262)+p3mup2*(p4mup4*F(383)+2*(p
     -   1mup4*F(263)+p2mup4*F(379)+p3mup4*F(1268)))+p3mup4*p4mup2*F(
     -   1609)
       Fa15(1439) = p2mup4*p4mup2*F(993)+4*(p3mup4*p4mup2*F(945)+p1mup
     -   4*(p4mup2*F(364)+p3mup2*F(947))+p3mup2*p3mup4*F(1025))+2*(p1
     -   mup2*(2*(p3mup4*F(266)+p1mup4*F(634))+p4mup4*F(684)+p2mup4*F
     -   (940))+p3mup2*p4mup4*F(995)+p4mup2*p4mup4*F(996)+p2mup4*p3mu
     -   p2*F(1302))
       Fa15(1440) = p2mup4*p3mup2*F(4464)-2*p3mup4*(p1mup2*F(357)+p3mu
     -   p2*F(4465)+p4mup2*F(4467))
       Fa15(1441) = p1mup2*(p4mup4*F(1001)+2*(p1mup4*F(678)+p2mup4*F(6
     -   80)+p3mup4*F(1292)))+2*p3mup2*(p4mup4*F(1010)+p1mup4*F(1302)
     -   +2*(p3mup4*F(311)+p2mup4*F(1303)))+p4mup2*(-(p1mup4*F(1002))
     -   +2*(p2mup4*F(1004)+p3mup4*F(1010))+4*p4mup4*F(4472))
       Fa15(1442) = p1mup4*p3mup2*F(383)+p1mup2*p3mup4*F(383)+2*(2*p3m
     -   up2*p3mup4*F(4465)+p3mup4*p4mup2*F(4467)+p3mup2*p4mup4*F(446
     -   7)+p2mup4*p3mup2*F(4469))
       Fa15(1443) = -(mup2mup5*(2*(p1mup4*F(263)+p3mup4*F(273)+p2mup4*
     -   F(379))+p4mup4*F(383)))+mup2mup4*(2*(p1mup5*F(263)+p3mup5*F(
     -   273)+p2mup5*F(379))+p4mup5*F(383))+mup4mup5*(p4mup2*F(1609)+
     -   2*(p1mup2*F(1262)+p3mup2*F(2368)))
       Fa15(1444) = mup4mup5*p3mup2+mup2mup5*p3mup4-mup2mup4*p3mup5
       Fa15(1445) = -2*p1mup2*p3mup4*F(1608)+p3mup4*p4mup2*F(1623)+p3m
     -   up2*(p4mup4*F(448)+2*(p1mup4*F(376)+p2mup4*F(377)+p3mup4*F(3
     -   496)))
       Fa15(1446) = p2mup4*p4mup2*F(1634)+2*(p4mup2*p4mup4*F(1630)+p3m
     -   up2*p4mup4*F(1635)+p1mup2*(p2mup4*F(1613)+p4mup4*F(1633)+2*p
     -   3mup4*F(3497))+p2mup4*p3mup2*F(3499))+4*(p1mup4*(p3mup2*F(26
     -   6)+p4mup2*F(354)+p1mup2*F(1611))+p3mup2*p3mup4*F(3498)+p3mup
     -   4*p4mup2*F(3500))
       Fa15(1447) = -(p2mup4*p3mup2*F(477))+2*(p3mup2*p3mup4*F(445)+p3
     -   mup4*p4mup2*F(459)+p1mup2*p3mup4*F(1629))
       Fa15(1448) = p1mup2*p4mup4*F(1624)+4*(p2mup4*p3mup2*F(883)+p3mu
     -   p2*p3mup4*F(3502))+2*(p1mup2*p1mup4*F(1617)+p1mup2*p2mup4*F(
     -   1618)+p1mup2*p3mup4*F(3501)+p1mup4*p3mup2*F(3503)+p3mup2*p4m
     -   up4*F(3504))+p4mup2*(4*p4mup4*F(314)+p1mup4*F(1625)+2*(p2mup
     -   4*F(1626)+p3mup4*F(3505)))
       Fa15(1449) = p1mup4*p3mup2*F(443)+p1mup2*p3mup4*F(443)+4*p3mup2
     -   *p3mup4*F(445)+2*(p2mup4*p3mup2*F(444)+(p3mup4*p4mup2+p3mup2
     -   *p4mup4)*F(459))
       Fa15(1450) = -2*mup2mup4*(p1mup5*F(263)+p3mup5*F(273)+p2mup5*F(
     -   379))+mup2mup4*p4mup5*F(485)+mup2mup5*(p4mup4*F(448)+2*(p1mu
     -   p4*F(376)+p2mup4*F(377)+p3mup4*F(3135)))+mup4mup5*(-(p4mup2*
     -   F(1623))+2*(p1mup2*F(1608)+p3mup2*F(3506)))
       Fa15(1451) = p3mup2*p4mup4*F(448)+p3mup4*p4mup2*F(1623)+2*(p1mu
     -   p4*p3mup2*F(376)+p2mup4*p3mup2*F(377)+p1mup2*p3mup4*F(1262)+
     -   p3mup2*p3mup4*F(3496))
       Fa15(1452) = p1mup2*p4mup4*F(1624)+4*(p2mup4*p3mup2*F(883)+p3mu
     -   p2*p3mup4*F(3502))+2*(p1mup2*p1mup4*F(1617)+p1mup2*p2mup4*F(
     -   1618)+p1mup2*p3mup4*F(3501)+p1mup4*p3mup2*F(3503)+p3mup2*p4m
     -   up4*F(3504))+p4mup2*(4*p4mup4*F(314)-p1mup4*F(1637)-2*(p2mup
     -   4*F(1638)+p3mup4*F(3507)))
       Fa15(1453) = mup2mup4*(2*(p1mup5*F(263)+p3mup5*F(273)+p2mup5*F(
     -   379))+p4mup5*F(447))+mup4mup5*(p4mup2*F(1623)+2*(p1mup2*F(12
     -   62)+p3mup2*F(2368)))-mup2mup5*(p4mup4*F(448)+2*(p1mup4*F(376
     -   )+p2mup4*F(377)+p3mup4*F(3135)))
       Fa15(1454) = p3mup2*p4mup4*F(383)+2*(p1mup4*p3mup2*F(263)+p2mup
     -   4*p3mup2*F(379)+p3mup2*p3mup4*F(1268)-p1mup2*p3mup4*F(1608))
     -   +p3mup4*p4mup2*F(1609)
       Fa15(1455) = -(p2mup4*p3mup2*F(4464))+2*p3mup4*(p1mup2*F(357)+p
     -   3mup2*F(4465)+p4mup2*F(4467))
       Fa15(1456) = 2*(mup2mup5*(p1mup4*F(263)+p3mup4*F(273)+p2mup4*F(
     -   379))-mup2mup4*(p1mup5*F(263)+p3mup5*F(273)+p2mup5*F(379)))+
     -   (mup2mup5*p4mup4-mup2mup4*p4mup5)*F(383)+mup4mup5*(-(p4mup2*
     -   F(1609))+2*(p1mup2*F(1608)+p3mup2*F(3506)))
       Fa15(1457) = 2*p3mup2*p4mup4*F(3413)+4*(p1mup4*(p3mup2*F(3406)+
     -   p4mup2*F(3408))+p2mup4*p3mup2*F(3409)+p2mup4*p4mup2*F(3410)-
     -   p4mup2*p4mup4*F(3414))+p1mup2*(-8*(p1mup4*F(3403)+p2mup4*F(3
     -   404))-4*p4mup4*F(3407)+p3mup4*F(3508))+p3mup4*p4mup2*F(3509)
     -   +p3mup2*p3mup4*F(3510)
       Fa15(1458) = p2mup4*p4mup2*F(3420)+4*(p3mup2*p3mup4*F(3422)+p1m
     -   up4*p3mup2*F(3423))+2*(p1mup4*p4mup2*F(3419)+p3mup4*p4mup2*F
     -   (3421)+p2mup4*p3mup2*F(3424)+p3mup2*p4mup4*F(3425)+p4mup2*p4
     -   mup4*F(3426)+p1mup2*(2*(p1mup4*F(3415)+p3mup4*F(3417))+p2mup
     -   4*F(3511)+p4mup4*F(3512)))
       Fa15(1459) = 2*p2mup4*p3mup2*F(3432)+8*(p1mup4*p4mup2*F(3433)+p
     -   4mup2*p4mup4*F(3437))+p1mup2*(-4*p2mup4*F(3428)-8*(p1mup4*F(
     -   3427)-p4mup4*F(3433))+p3mup4*F(3513))+p3mup2*p3mup4*F(3514)+
     -   4*(p1mup4*p3mup2*F(3430)+p2mup4*p4mup2*F(3434)+p3mup2*p4mup4
     -   *F(3515))+p3mup4*p4mup2*F(3516)
       Fa15(1460) = -8*(p2mup4*p3mup2*F(3444)+p3mup2*p3mup4*F(3446)+p4
     -   mup2*p4mup4*F(3447))-2*p1mup2*(p4mup4*F(3442)+2*(p1mup4*F(34
     -   38)+p2mup4*F(3439)+p3mup4*F(3517)))+p1mup4*p4mup2*F(3518)-4*
     -   (p1mup4*p3mup2*F(3441)+p3mup2*p4mup4*F(3448)-p2mup4*p4mup2*F
     -   (3519)-p3mup4*p4mup2*F(3520))
       Fa15(1461) = -2*p3mup2*(4*p3mup4*F(3453)+p1mup4*F(3454)+2*(p2mu
     -   p4*F(3455)+p4mup4*F(3461)))+p1mup2*(p3mup4*F(3521)-4*(p1mup4
     -   *F(3450)+p2mup4*F(3451)-p4mup4*F(3522)))+p4mup2*(4*p1mup4*F(
     -   3457)+8*(p2mup4*F(3458)+p4mup4*F(3460))+p3mup4*F(3523))
       Fa15(1462) = 4*(p4mup2*(p1mup4*F(3469)+p2mup4*F(3470)+2*(p3mup4
     -   *F(3471)+p4mup4*F(3472)))+p3mup2*(p1mup4*F(3466)+p2mup4*F(34
     -   67)+2*(p3mup4*F(3465)+p4mup4*F(3473))))+p1mup2*(2*p2mup4*F(3
     -   463)+4*(p1mup4*F(3462)+p3mup4*F(3464))+p4mup4*F(3524))
       Fa15(1463) = mup2mup4*(4*(p1mup5*F(3312)+p2mup5*F(3313)+p3mup5*
     -   F(3476))+p4mup5*F(3525))-2*(mup2mup5*(2*(p1mup4*F(3474)+p2mu
     -   p4*F(3475)+p3mup4*F(3478))+p4mup4*F(3479))-mup4mup5*(p4mup2*
     -   F(3481)+2*(p1mup2*F(3310)+p3mup2*F(3526))))
       Fa15(1464) = mup2mup5*(-4*(p1mup4*F(3001)-p4mup4*F(3007)+p2mup4
     -   *F(3486))+p3mup4*F(3527))+2*(mup4mup5*(-(p3mup2*F(3489))+2*(
     -   p1mup2*F(3311)+p4mup2*F(3490)))+mup2mup4*(2*(p1mup5*F(3305)+
     -   p4mup5*F(3306)+p2mup5*F(3314))+p3mup5*F(3528)))
       Fa15(1465) = 2*(mup2mup5*(2*(p1mup4*F(2525)-p4mup4*F(2535)+p3mu
     -   p4*F(2565))+p2mup4*F(3492))-mup2mup4*(2*(p1mup5*F(2527)+p4mu
     -   p5*F(2531)+p3mup5*F(2566))+p2mup5*F(3493)))+4*mup4mup5*(p1mu
     -   p2*F(3491)+p3mup2*F(3494)+p4mup2*F(3495))
       Fa15(1466) = -4*mup1mup2*(p1mup5*F(3312)+p2mup5*F(3313)+p3mup5*
     -   F(3476))+2*mup1mup2*p4mup5*F(3477)
       Fa15(1467) = mup1mup2*(4*(p1mup5*F(3002)+p4mup5*F(3003)+p2mup5*
     -   F(3382))+p3mup5*F(3488))
       Fa15(1468) = mup1mup2*(p1mup5*F(3482)+2*(p2mup5*F(3483)+p3mup5*
     -   F(3484)+p4mup5*F(3485)))
       Fa15(1469) = 2*mup1mup2*(2*(p1mup5*F(2527)+p4mup5*F(2531)+p3mup
     -   5*F(2566))+p2mup5*F(3493))
       Fa15(1470) = mup1mup2*(2*(p1mup5*F(263)+p3mup5*F(273)+p2mup5*F(
     -   379))+p4mup5*F(383))
       Fa15(1471) = mup1mup2*(2*(p1mup5*F(263)+p3mup5*F(273)+p2mup5*F(
     -   379))-p4mup5*F(485))
       Fa15(1472) = mup1mup2*(2*(p1mup5*F(263)+p3mup5*F(273)+p2mup5*F(
     -   379))+p4mup5*F(447))
       Fa15(1473) = mup1mup2*(4*(p1mup5*F(3312)+p2mup5*F(3313)+p3mup5*
     -   F(3476))+p4mup5*F(3525))
       Fa15(1474) = 2*mup1mup2*(2*(p1mup5*F(3305)+p4mup5*F(3306)+p2mup
     -   5*F(3314))+p3mup5*F(3528))
       Fa15(1475) = 2*mup1mup2*(p1mup5*F(3482)+2*(p2mup5*F(3483)+p3mup
     -   5*F(3484)+p4mup5*F(3485)))
       Fa15(1476) = mup1mup2*(2*(p1mup5*F(2527)+p4mup5*F(2531)+p3mup5*
     -   F(2566))+p2mup5*F(3493))
       Fa15(1477) = mup1mup5*F(3008)-4*p3mup1*(p1mup5*F(3312)+p2mup5*F
     -   (3313)+p3mup5*F(3476))+2*p3mup1*p4mup5*F(3477)
       Fa15(1478) = mup1mup5*F(3011)+2*p3mup1*(p1mup5*F(3482)+2*(p2mup
     -   5*F(3483)+p3mup5*F(3484)+p4mup5*F(3485)))
       Fa15(1479) = 2*p4mup1*(p1mup4*F(3482)+2*(p2mup4*F(3483)+p3mup4*
     -   F(3484)+p4mup4*F(3485)))
       Fa15(1480) = mup1mup4*(4*(p1mup5*F(3312)+p2mup5*F(3313)+p3mup5*
     -   F(3476))+p4mup5*F(3525))
       Fa15(1481) = 2*mup1mup4*(p1mup5*F(3482)+2*(p2mup5*F(3483)+p3mup
     -   5*F(3484)+p4mup5*F(3485)))
       Fa15(1482) = 2*mup1mup4*(2*(p1mup5*F(3305)+p4mup5*F(3306)+p2mup
     -   5*F(3314))+p3mup5*F(3528))
       Fa15(1483) = mup1mup5*F(3009)+p2mup1*(4*(p1mup5*F(3312)+p2mup5*
     -   F(3313)+p3mup5*F(3476))+p4mup5*F(3525))
       Fa15(1484) = mup1mup5*F(3000)-4*(p1mup5*(p3mup1*F(2527)-p2mup1*
     -   F(3305))-p2mup1*p4mup5*F(3306)-p2mup1*p2mup5*F(3314))-2*(p3m
     -   up1*(2*(p4mup5*F(2531)+p3mup5*F(2566))+p2mup5*F(3493))-p2mup
     -   1*p3mup5*F(3528))
       Fa15(1485) = mup1mup5*F(3012)-2*p2mup1*(p1mup5*F(3482)+2*(p2mup
     -   5*F(3483)+p3mup5*F(3484)+p4mup5*F(3485)))
       Return
       End
