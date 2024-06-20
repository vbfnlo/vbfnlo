c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenEvenGSim11(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p
     -   2mup1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup
     -   4,p3mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5
     -   mup3,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2m
     -   up3,mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa11)
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
       Complex*16 Fa11(81:88)
       COMMON/FaPenEvenGSim/Fa
       Fa11(81) = p1mup4*p2mup1*p2mup5*p3mup2*F(3168)+p1mup5*p2mup2*p2
     -   mup4*p3mup1*F(3171)+p1mup4*p2mup2*p3mup1*p3mup5*F(3197)+4*(p
     -   1mup4*p2mup2*p2mup5*p3mup1*F(3169)-p1mup5*p2mup1*p2mup4*p3mu
     -   p2*F(3170)-p1mup4*p1mup5*p2mup1*p3mup2*F(3190)+p1mup4*p1mup5
     -   *p2mup2*p3mup1*F(3191)-p1mup4*p1mup5*p2mup1*p4mup2*F(3192)+p
     -   1mup4*p1mup5*p2mup2*p4mup1*F(3193)-p1mup4*p2mup1*p3mup2*p3mu
     -   p5*F(3196)+p1mup4*p1mup5*p3mup1*p4mup2*F(3198)-p1mup4*p1mup5
     -   *p3mup2*p4mup1*F(3199))+16*p1mup4*p2mup1*p2mup2*p2mup5*F(320
     -   1)+p1mup4*p2mup1*p2mup2*p3mup5*F(3202)+p1mup5*p2mup1*p2mup2*
     -   p3mup4*F(3204)+p1mup4*p2mup5*p3mup1*p3mup2*F(3206)+p1mup5*p2
     -   mup1*p3mup2*p3mup4*F(3208)+p1mup5*p2mup2*p3mup1*p3mup4*F(320
     -   9)+p1mup5*p2mup4*p3mup1*p3mup2*F(3210)+p1mup2*(p2mup1*p2mup5
     -   *p3mup4*F(3167)+p2mup1*p3mup4*p3mup5*F(3184)+4*(p2mup1*p2mup
     -   4*p3mup5*F(3166)+p2mup4*p2mup5*p3mup1*F(3185)-p2mup4*p3mup1*
     -   p3mup5*F(3186))+p2mup5*p3mup1*p3mup4*F(3187)+8*(p2mup1*p2mup
     -   4*p2mup5*F(3165)-p3mup1*p3mup4*p3mup5*F(3188))+p1mup4*(p2mup
     -   1*p3mup5*F(3154)+4*(p2mup5*p3mup1*F(3155)+p2mup1*p4mup5*F(31
     -   58)+p2mup5*p4mup1*F(3159)-p3mup1*p4mup5*F(3175))+p3mup5*p4mu
     -   p1*F(3176)+8*(p2mup1*p2mup5*F(3153)-p3mup1*p3mup5*F(3173)+p4
     -   mup1*p4mup5*F(3195)))+p1mup5*(p2mup4*p3mup1*F(3157)+p2mup1*(
     -   8*p2mup4*F(3152)+p3mup4*F(3156)+4*p4mup4*F(3160))-4*p2mup4*p
     -   4mup1*F(3161)+p3mup1*p4mup4*F(3177)+p3mup4*p4mup1*F(3178)+8*
     -   (p3mup1*p3mup4*F(3174)+p1mup4*p4mup1*F(3189)+p4mup1*p4mup4*F
     -   (3212))))+8*(p1mup4*p1mup5*p2mup1*p2mup2*F(3200)+p1mup5*p2mu
     -   p1*p2mup2*p2mup4*F(3203)-p1mup4*p1mup5*p3mup1*p3mup2*F(3205)
     -   -p1mup4*p3mup1*p3mup2*p3mup5*F(3207)-p1mup5*p3mup1*p3mup2*p3
     -   mup4*F(3211)+p1mup4*p1mup5*p4mup1*p4mup2*F(3213))+p1mup1*(p1
     -   mup5*p2mup2*p3mup4*F(3136)+p2mup2*p2mup4*p3mup5*F(3163)+p2mu
     -   p2*p2mup5*p3mup4*F(3164)+p2mup2*p3mup4*p3mup5*F(3179)+4*(p1m
     -   up5*p2mup4*p3mup2*F(3137)+p1mup5*p2mup2*p4mup4*F(3142)+p1mup
     -   5*p2mup4*p4mup2*F(3143)-p1mup5*p3mup2*p4mup4*F(3150)-p1mup5*
     -   p3mup4*p4mup2*F(3151)+p2mup4*p2mup5*p3mup2*F(3180))+p2mup4*p
     -   3mup2*p3mup5*F(3181)+p2mup5*p3mup2*p3mup4*F(3182)+p1mup2*(p2
     -   mup4*(8*p2mup5*F(3130)+p3mup5*F(3131)+4*p4mup5*F(3138))+p3mu
     -   p4*(p2mup5*F(3132)-8*p3mup5*F(3144)+p4mup5*F(3146))+p4mup4*(
     -   4*p2mup5*F(3139)+p3mup5*F(3147)+8*(p1mup5*F(3129)+p4mup5*F(3
     -   194))))+p1mup4*(p2mup2*p3mup5*F(3133)+p2mup5*p3mup2*F(3134)+
     -   4*(p2mup2*p4mup5*F(3140)-p2mup5*p4mup2*F(3141)-p3mup2*p4mup5
     -   *F(3148)+p3mup5*p4mup2*F(3149))+8*(p2mup2*p2mup5*F(3152)-p3m
     -   up2*p3mup5*F(3172)+p4mup2*p4mup5*F(3212))+16*p1mup5*p4mup2*F
     -   (3214))+8*(p1mup5*p2mup2*p2mup4*F(3135)-p1mup5*p3mup2*p3mup4
     -   *F(3145)+p2mup2*p2mup4*p2mup5*F(3162)-p3mup2*p3mup4*p3mup5*F
     -   (3183)+p1mup5*p4mup2*p4mup4*F(3215)))
       Fa11(82) = -(p1mup1*p1mup2*p2mup5*p3mup4*F(1045))+p1mup1*p1mup2
     -   *p2mup4*p3mup5*F(1045)+p1mup1*p1mup2*p3mup5*p4mup4*F(1052)-p
     -   1mup1*p1mup2*p3mup4*p4mup5*F(1052)-p1mup1*p2mup2*p2mup5*p3mu
     -   p4*F(1188)+p1mup1*p2mup2*p2mup4*p3mup5*F(1188)-p1mup2*p2mup1
     -   *p2mup5*p3mup4*F(1189)+p1mup2*p2mup1*p2mup4*p3mup5*F(1189)-p
     -   1mup1*p2mup5*p3mup2*p3mup4*F(3216)+p1mup1*p2mup4*p3mup2*p3mu
     -   p5*F(3216)-p1mup2*p3mup1*(p2mup5*p3mup4-p2mup4*p3mup5)*F(321
     -   7)+p1mup4*(-(p1mup1*p2mup5*p3mup2*F(1047))+p1mup2*p2mup1*p3m
     -   up5*F(1048)+p1mup2*p2mup5*p3mup1*F(1049)-p2mup1*p2mup5*p3mup
     -   2*F(1049)+p1mup2*p3mup5*p4mup1*F(1052)+p1mup2*p3mup1*p4mup5*
     -   F(1192)+2*(p1mup2*p3mup1*p3mup5*F(1050)+p1mup1*p3mup2*p3mup5
     -   *F(1055)+p1mup1*p3mup5*p4mup2*F(1056)+p3mup1*p3mup2*p3mup5*F
     -   (1326))+p2mup2*(p2mup5*p3mup1*F(932)+p1mup1*p3mup5*F(1046)+p
     -   2mup1*p3mup5*F(1193)+p3mup1*p3mup5*F(3217))+p2mup1*p3mup2*p3
     -   mup5*F(3218)+p2mup5*p3mup1*p3mup2*F(3220))-p1mup5*(-(p1mup1*
     -   p2mup4*p3mup2*F(1047))+p1mup2*p2mup1*p3mup4*F(1048)+p1mup2*p
     -   2mup4*p3mup1*F(1049)-p2mup1*p2mup4*p3mup2*F(1049)+p1mup2*p3m
     -   up4*p4mup1*F(1052)+p1mup2*p3mup1*p4mup4*F(1192)+2*(p1mup2*p3
     -   mup1*p3mup4*F(1050)+p1mup1*p3mup2*p3mup4*F(1055)+p1mup1*p3mu
     -   p4*p4mup2*F(1056)+p3mup1*p3mup2*p3mup4*F(1326))+p2mup2*(p2mu
     -   p4*p3mup1*F(932)+p1mup1*p3mup4*F(1046)+p2mup1*p3mup4*F(1193)
     -   +p3mup1*p3mup4*F(3217))-p2mup1*p3mup2*p3mup4*F(3219)-p2mup4*
     -   p3mup1*p3mup2*F(3221))
       Fa11(83) = -(p1mup1*p1mup2*p2mup5*p3mup4*F(1045))+p1mup1*p1mup2
     -   *p2mup4*p3mup5*F(1045)+p1mup1*p1mup2*p3mup5*p4mup4*F(1052)-p
     -   1mup1*p1mup2*p3mup4*p4mup5*F(1052)-p1mup1*p2mup2*p2mup5*p3mu
     -   p4*F(1188)+p1mup1*p2mup2*p2mup4*p3mup5*F(1188)-p1mup2*p2mup1
     -   *p2mup5*p3mup4*F(1189)+p1mup2*p2mup1*p2mup4*p3mup5*F(1189)-p
     -   1mup1*p2mup5*p3mup2*p3mup4*F(3216)+p1mup1*p2mup4*p3mup2*p3mu
     -   p5*F(3216)-p1mup2*p3mup1*(p2mup5*p3mup4-p2mup4*p3mup5)*F(321
     -   7)-p1mup5*(-(p1mup1*p2mup4*p3mup2*F(1047))+p1mup2*p2mup1*p3m
     -   up4*F(1048)+p1mup2*p2mup4*p3mup1*F(1049)-p2mup1*p2mup4*p3mup
     -   2*F(1049)+p1mup2*p3mup4*p4mup1*F(1052)+p1mup2*p3mup1*p4mup4*
     -   F(1192)+2*(p1mup2*p3mup1*p3mup4*F(1050)+p1mup1*p3mup2*p3mup4
     -   *F(1055)+p1mup1*p3mup4*p4mup2*F(1056)+p3mup1*p3mup2*p3mup4*F
     -   (1326))+p2mup2*(p2mup4*p3mup1*F(932)+p1mup1*p3mup4*F(1046)+p
     -   2mup1*p3mup4*F(1193)+p3mup1*p3mup4*F(3217))+p2mup1*p3mup2*p3
     -   mup4*F(3218)+p2mup4*p3mup1*p3mup2*F(3220))+p1mup4*(-(p1mup1*
     -   p2mup5*p3mup2*F(1047))+p1mup2*p2mup1*p3mup5*F(1048)+p1mup2*p
     -   2mup5*p3mup1*F(1049)-p2mup1*p2mup5*p3mup2*F(1049)+p1mup2*p3m
     -   up5*p4mup1*F(1052)+p1mup2*p3mup1*p4mup5*F(1192)+2*(p1mup2*p3
     -   mup1*p3mup5*F(1050)+p1mup1*p3mup2*p3mup5*F(1055)+p1mup1*p3mu
     -   p5*p4mup2*F(1056)+p3mup1*p3mup2*p3mup5*F(1326))+p2mup2*(p2mu
     -   p5*p3mup1*F(932)+p1mup1*p3mup5*F(1046)+p2mup1*p3mup5*F(1193)
     -   +p3mup1*p3mup5*F(3217))-p2mup1*p3mup2*p3mup5*F(3219)-p2mup5*
     -   p3mup1*p3mup2*F(3221))
       Fa11(84) = -(p1mup2*p1mup5*p2mup4*p3mup1*F(1038))+p1mup2*p1mup5
     -   *p2mup1*p3mup4*F(1048)+p1mup2*p2mup1*p2mup4*p3mup5*F(1162)-p
     -   1mup2*p2mup4*p2mup5*p3mup1*F(1163)-p1mup5*p2mup2*p2mup4*p3mu
     -   p1*F(1165)+p1mup5*p2mup1*p2mup4*p3mup2*F(1166)+p1mup2*p2mup1
     -   *p3mup4*p3mup5*F(1172)+p1mup2*p2mup1*p2mup5*p3mup4*F(1189)+p
     -   1mup5*p2mup1*p2mup2*p3mup4*F(1193)-p1mup2*p1mup5*p3mup1*p4mu
     -   p4*F(1201)-p1mup2*p1mup5*p3mup4*p4mup1*F(1202)+4*p1mup2*p3mu
     -   p1*p3mup4*p3mup5*F(1307)+p1mup1*(p1mup5*p2mup2*p3mup4*F(1046
     -   )-p1mup5*p2mup4*p3mup2*F(1047)+p1mup4*(-(p2mup5*p3mup2*F(929
     -   ))+p3mup5*(p2mup2*F(1035)+2*(p3mup2*F(1050)+p4mup2*F(1053)))
     -   )+p1mup2*(p2mup4*p3mup5*F(1034)+2*p3mup4*p3mup5*F(1041)+p2mu
     -   p5*p3mup4*F(1045)+(p3mup5*p4mup4+p3mup4*p4mup5)*F(1058))+p2m
     -   up2*p2mup4*p3mup5*F(1160)-p2mup4*p2mup5*p3mup2*F(1161)+p2mup
     -   2*p3mup4*p3mup5*F(1167)+p2mup2*p2mup5*p3mup4*F(1188)+2*(p1mu
     -   p5*p3mup2*p3mup4*F(1055)+p1mup5*p3mup4*p4mup2*F(1056)+p3mup2
     -   *p3mup4*p3mup5*F(1306))+p2mup4*p3mup2*p3mup5*F(3216)+p2mup5*
     -   p3mup2*p3mup4*F(3223))+p1mup2*p2mup5*p3mup1*p3mup4*F(3224)+p
     -   1mup5*p2mup1*p3mup2*p3mup4*F(3225)+p1mup5*p2mup2*p3mup1*p3mu
     -   p4*F(3226)+p1mup2*p2mup4*p3mup1*p3mup5*F(3227)-p1mup5*p2mup4
     -   *p3mup1*p3mup2*F(3229)+p1mup4*(p1mup5*p2mup1*p3mup2*F(1039)-
     -   p2mup1*p2mup5*p3mup2*F(1049)+p1mup2*(p2mup1*p3mup5*F(1036)-p
     -   2mup5*p3mup1*F(1037)+2*p3mup1*p3mup5*F(1050)+p3mup5*p4mup1*F
     -   (1200)-p3mup1*p4mup5*F(1201))-2*(p1mup5*p3mup1*p3mup2*F(1051
     -   )+p1mup5*p3mup1*p4mup2*F(1192)-p3mup1*p3mup2*p3mup5*F(1326))
     -   +p2mup2*(-(p2mup5*p3mup1*F(932))-p1mup5*p3mup1*F(1040)+p2mup
     -   1*p3mup5*F(1190)+p3mup1*p3mup5*F(3217))+p2mup1*p3mup2*p3mup5
     -   *F(3228)-p2mup5*p3mup1*p3mup2*F(3229))+2*(p1mup2*p1mup5*p3mu
     -   p1*p3mup4*F(3222)+p1mup5*p3mup1*p3mup2*p3mup4*F(3230))
       Fa11(85) = p1mup2*p2mup1*(p2mup3*p3mup5*F(3234)+p2mup5*p3mup3*F
     -   (3235)-p2mup3*p4mup5*F(3238)-p2mup5*p4mup3*F(3239))+p1mup1*(
     -   4*p1mup2*p1mup3*p1mup5*F(3231)+p2mup2*(p2mup3*p3mup5*F(3232)
     -   +p2mup5*p3mup3*F(3233)-p2mup3*p4mup5*F(3236)-p2mup5*p4mup3*F
     -   (3237))+p2mup3*p2mup5*(p3mup2*F(3240)-p4mup2*F(3241)))
       Fa11(86) = p1mup1*(p2mup2*p2mup3*p4mup5*F(1160)-p2mup3*p2mup5*p
     -   4mup2*F(1161)-p2mup2*p2mup5*p4mup3*F(1188))+p1mup2*p2mup1*(p
     -   2mup3*p4mup5*F(1162)-p2mup5*p4mup3*F(1189))
       Fa11(87) = (p2mup5*p4mup3+p2mup3*p4mup5)*(p1mup1*p2mup2*F(1188)
     -   +p1mup2*p2mup1*F(1189))
       Fa11(88) = p1mup1*(-(p2mup2*p2mup3*p4mup5*F(1160))+p2mup3*p2mup
     -   5*p4mup2*F(1161)+p2mup2*p2mup5*p4mup3*F(1188))+p1mup2*p2mup1
     -   *(-(p2mup3*p4mup5*F(1162))+p2mup5*p4mup3*F(1189))
       Return
       End
