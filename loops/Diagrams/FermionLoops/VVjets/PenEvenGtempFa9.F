c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenEvenG9(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup
     -   1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3
     -   mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3
     -   ,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,
     -   mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa9)
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
       Complex*16 Fa9(81:90)
       COMMON/FaPenEvenG/Fa
       Fa9(81) = p1mup4*p2mup1*p2mup5*p3mup2*F(3168)+p1mup5*p2mup2*p2m
     -   up4*p3mup1*F(3171)+p1mup4*p2mup2*p3mup1*p3mup5*F(3197)+4*(p1
     -   mup4*p2mup2*p2mup5*p3mup1*F(3169)-p1mup5*p2mup1*p2mup4*p3mup
     -   2*F(3170)-p1mup4*p1mup5*p2mup1*p3mup2*F(3190)+p1mup4*p1mup5*
     -   p2mup2*p3mup1*F(3191)-p1mup4*p1mup5*p2mup1*p4mup2*F(3192)+p1
     -   mup4*p1mup5*p2mup2*p4mup1*F(3193)-p1mup4*p2mup1*p3mup2*p3mup
     -   5*F(3196)+p1mup4*p1mup5*p3mup1*p4mup2*F(3198)-p1mup4*p1mup5*
     -   p3mup2*p4mup1*F(3199))+16*p1mup4*p2mup1*p2mup2*p2mup5*F(3201
     -   )+p1mup4*p2mup1*p2mup2*p3mup5*F(3202)+p1mup5*p2mup1*p2mup2*p
     -   3mup4*F(3204)+p1mup4*p2mup5*p3mup1*p3mup2*F(3206)+p1mup5*p2m
     -   up1*p3mup2*p3mup4*F(3208)+p1mup5*p2mup2*p3mup1*p3mup4*F(3209
     -   )+p1mup5*p2mup4*p3mup1*p3mup2*F(3210)+p1mup2*(p2mup1*p2mup5*
     -   p3mup4*F(3167)+p2mup1*p3mup4*p3mup5*F(3184)+4*(p2mup1*p2mup4
     -   *p3mup5*F(3166)+p2mup4*p2mup5*p3mup1*F(3185)-p2mup4*p3mup1*p
     -   3mup5*F(3186))+p2mup5*p3mup1*p3mup4*F(3187)+8*(p2mup1*p2mup4
     -   *p2mup5*F(3165)-p3mup1*p3mup4*p3mup5*F(3188))+p1mup4*(p2mup1
     -   *p3mup5*F(3154)+4*(p2mup5*p3mup1*F(3155)+p2mup1*p4mup5*F(315
     -   8)+p2mup5*p4mup1*F(3159)-p3mup1*p4mup5*F(3175))+p3mup5*p4mup
     -   1*F(3176)+8*(p2mup1*p2mup5*F(3153)-p3mup1*p3mup5*F(3173)+p4m
     -   up1*p4mup5*F(3195)))+p1mup5*(p2mup4*p3mup1*F(3157)+p2mup1*(8
     -   *p2mup4*F(3152)+p3mup4*F(3156)+4*p4mup4*F(3160))-4*p2mup4*p4
     -   mup1*F(3161)+p3mup1*p4mup4*F(3177)+p3mup4*p4mup1*F(3178)+8*(
     -   p3mup1*p3mup4*F(3174)+p1mup4*p4mup1*F(3189)+p4mup1*p4mup4*F(
     -   3212))))+8*(p1mup4*p1mup5*p2mup1*p2mup2*F(3200)+p1mup5*p2mup
     -   1*p2mup2*p2mup4*F(3203)-p1mup4*p1mup5*p3mup1*p3mup2*F(3205)-
     -   p1mup4*p3mup1*p3mup2*p3mup5*F(3207)-p1mup5*p3mup1*p3mup2*p3m
     -   up4*F(3211)+p1mup4*p1mup5*p4mup1*p4mup2*F(3213))+p1mup1*(p1m
     -   up5*p2mup2*p3mup4*F(3136)+p2mup2*p2mup4*p3mup5*F(3163)+p2mup
     -   2*p2mup5*p3mup4*F(3164)+p2mup2*p3mup4*p3mup5*F(3179)+4*(p1mu
     -   p5*p2mup4*p3mup2*F(3137)+p1mup5*p2mup2*p4mup4*F(3142)+p1mup5
     -   *p2mup4*p4mup2*F(3143)-p1mup5*p3mup2*p4mup4*F(3150)-p1mup5*p
     -   3mup4*p4mup2*F(3151)+p2mup4*p2mup5*p3mup2*F(3180))+p2mup4*p3
     -   mup2*p3mup5*F(3181)+p2mup5*p3mup2*p3mup4*F(3182)+p1mup2*(p2m
     -   up4*(8*p2mup5*F(3130)+p3mup5*F(3131)+4*p4mup5*F(3138))+p3mup
     -   4*(p2mup5*F(3132)-8*p3mup5*F(3144)+p4mup5*F(3146))+p4mup4*(4
     -   *p2mup5*F(3139)+p3mup5*F(3147)+8*(p1mup5*F(3129)+p4mup5*F(31
     -   94))))+p1mup4*(p2mup2*p3mup5*F(3133)+p2mup5*p3mup2*F(3134)+4
     -   *(p2mup2*p4mup5*F(3140)-p2mup5*p4mup2*F(3141)-p3mup2*p4mup5*
     -   F(3148)+p3mup5*p4mup2*F(3149))+8*(p2mup2*p2mup5*F(3152)-p3mu
     -   p2*p3mup5*F(3172)+p4mup2*p4mup5*F(3212))+16*p1mup5*p4mup2*F(
     -   3214))+8*(p1mup5*p2mup2*p2mup4*F(3135)-p1mup5*p3mup2*p3mup4*
     -   F(3145)+p2mup2*p2mup4*p2mup5*F(3162)-p3mup2*p3mup4*p3mup5*F(
     -   3183)+p1mup5*p4mup2*p4mup4*F(3215)))
       Fa9(82) = -(p1mup1*p1mup2*p2mup5*p3mup4*F(1045))+p1mup1*p1mup2*
     -   p2mup4*p3mup5*F(1045)+p1mup1*p1mup2*p3mup5*p4mup4*F(1052)-p1
     -   mup1*p1mup2*p3mup4*p4mup5*F(1052)-p1mup1*p2mup2*p2mup5*p3mup
     -   4*F(1188)+p1mup1*p2mup2*p2mup4*p3mup5*F(1188)-p1mup2*p2mup1*
     -   p2mup5*p3mup4*F(1189)+p1mup2*p2mup1*p2mup4*p3mup5*F(1189)-p1
     -   mup1*p2mup5*p3mup2*p3mup4*F(3216)+p1mup1*p2mup4*p3mup2*p3mup
     -   5*F(3216)-p1mup2*p3mup1*(p2mup5*p3mup4-p2mup4*p3mup5)*F(3217
     -   )+p1mup4*(-(p1mup1*p2mup5*p3mup2*F(1047))+p1mup2*p2mup1*p3mu
     -   p5*F(1048)+p1mup2*p2mup5*p3mup1*F(1049)-p2mup1*p2mup5*p3mup2
     -   *F(1049)+p1mup2*p3mup5*p4mup1*F(1052)+p1mup2*p3mup1*p4mup5*F
     -   (1192)+2*(p1mup2*p3mup1*p3mup5*F(1050)+p1mup1*p3mup2*p3mup5*
     -   F(1055)+p1mup1*p3mup5*p4mup2*F(1056)+p3mup1*p3mup2*p3mup5*F(
     -   1326))+p2mup2*(p2mup5*p3mup1*F(932)+p1mup1*p3mup5*F(1046)+p2
     -   mup1*p3mup5*F(1193)+p3mup1*p3mup5*F(3217))+p2mup1*p3mup2*p3m
     -   up5*F(3218)+p2mup5*p3mup1*p3mup2*F(3220))-p1mup5*(-(p1mup1*p
     -   2mup4*p3mup2*F(1047))+p1mup2*p2mup1*p3mup4*F(1048)+p1mup2*p2
     -   mup4*p3mup1*F(1049)-p2mup1*p2mup4*p3mup2*F(1049)+p1mup2*p3mu
     -   p4*p4mup1*F(1052)+p1mup2*p3mup1*p4mup4*F(1192)+2*(p1mup2*p3m
     -   up1*p3mup4*F(1050)+p1mup1*p3mup2*p3mup4*F(1055)+p1mup1*p3mup
     -   4*p4mup2*F(1056)+p3mup1*p3mup2*p3mup4*F(1326))+p2mup2*(p2mup
     -   4*p3mup1*F(932)+p1mup1*p3mup4*F(1046)+p2mup1*p3mup4*F(1193)+
     -   p3mup1*p3mup4*F(3217))-p2mup1*p3mup2*p3mup4*F(3219)-p2mup4*p
     -   3mup1*p3mup2*F(3221))
       Fa9(83) = -(p1mup1*p1mup2*p2mup5*p3mup4*F(1045))+p1mup1*p1mup2*
     -   p2mup4*p3mup5*F(1045)+p1mup1*p1mup2*p3mup5*p4mup4*F(1052)-p1
     -   mup1*p1mup2*p3mup4*p4mup5*F(1052)-p1mup1*p2mup2*p2mup5*p3mup
     -   4*F(1188)+p1mup1*p2mup2*p2mup4*p3mup5*F(1188)-p1mup2*p2mup1*
     -   p2mup5*p3mup4*F(1189)+p1mup2*p2mup1*p2mup4*p3mup5*F(1189)-p1
     -   mup1*p2mup5*p3mup2*p3mup4*F(3216)+p1mup1*p2mup4*p3mup2*p3mup
     -   5*F(3216)-p1mup2*p3mup1*(p2mup5*p3mup4-p2mup4*p3mup5)*F(3217
     -   )-p1mup5*(-(p1mup1*p2mup4*p3mup2*F(1047))+p1mup2*p2mup1*p3mu
     -   p4*F(1048)+p1mup2*p2mup4*p3mup1*F(1049)-p2mup1*p2mup4*p3mup2
     -   *F(1049)+p1mup2*p3mup4*p4mup1*F(1052)+p1mup2*p3mup1*p4mup4*F
     -   (1192)+2*(p1mup2*p3mup1*p3mup4*F(1050)+p1mup1*p3mup2*p3mup4*
     -   F(1055)+p1mup1*p3mup4*p4mup2*F(1056)+p3mup1*p3mup2*p3mup4*F(
     -   1326))+p2mup2*(p2mup4*p3mup1*F(932)+p1mup1*p3mup4*F(1046)+p2
     -   mup1*p3mup4*F(1193)+p3mup1*p3mup4*F(3217))+p2mup1*p3mup2*p3m
     -   up4*F(3218)+p2mup4*p3mup1*p3mup2*F(3220))+p1mup4*(-(p1mup1*p
     -   2mup5*p3mup2*F(1047))+p1mup2*p2mup1*p3mup5*F(1048)+p1mup2*p2
     -   mup5*p3mup1*F(1049)-p2mup1*p2mup5*p3mup2*F(1049)+p1mup2*p3mu
     -   p5*p4mup1*F(1052)+p1mup2*p3mup1*p4mup5*F(1192)+2*(p1mup2*p3m
     -   up1*p3mup5*F(1050)+p1mup1*p3mup2*p3mup5*F(1055)+p1mup1*p3mup
     -   5*p4mup2*F(1056)+p3mup1*p3mup2*p3mup5*F(1326))+p2mup2*(p2mup
     -   5*p3mup1*F(932)+p1mup1*p3mup5*F(1046)+p2mup1*p3mup5*F(1193)+
     -   p3mup1*p3mup5*F(3217))-p2mup1*p3mup2*p3mup5*F(3219)-p2mup5*p
     -   3mup1*p3mup2*F(3221))
       Fa9(84) = -(p1mup2*p1mup5*p2mup4*p3mup1*F(1038))+p1mup2*p1mup5*
     -   p2mup1*p3mup4*F(1048)+p1mup2*p2mup1*p2mup4*p3mup5*F(1162)-p1
     -   mup2*p2mup4*p2mup5*p3mup1*F(1163)-p1mup5*p2mup2*p2mup4*p3mup
     -   1*F(1165)+p1mup5*p2mup1*p2mup4*p3mup2*F(1166)+p1mup2*p2mup1*
     -   p3mup4*p3mup5*F(1172)+p1mup2*p2mup1*p2mup5*p3mup4*F(1189)+p1
     -   mup5*p2mup1*p2mup2*p3mup4*F(1193)-p1mup2*p1mup5*p3mup1*p4mup
     -   4*F(1201)-p1mup2*p1mup5*p3mup4*p4mup1*F(1202)+4*p1mup2*p3mup
     -   1*p3mup4*p3mup5*F(1307)+p1mup1*(p1mup5*p2mup2*p3mup4*F(1046)
     -   -p1mup5*p2mup4*p3mup2*F(1047)+p1mup4*(-(p2mup5*p3mup2*F(929)
     -   )+p3mup5*(p2mup2*F(1035)+2*(p3mup2*F(1050)+p4mup2*F(1053))))
     -   +p1mup2*(p2mup4*p3mup5*F(1034)+2*p3mup4*p3mup5*F(1041)+p2mup
     -   5*p3mup4*F(1045)+(p3mup5*p4mup4+p3mup4*p4mup5)*F(1058))+p2mu
     -   p2*p2mup4*p3mup5*F(1160)-p2mup4*p2mup5*p3mup2*F(1161)+p2mup2
     -   *p3mup4*p3mup5*F(1167)+p2mup2*p2mup5*p3mup4*F(1188)+2*(p1mup
     -   5*p3mup2*p3mup4*F(1055)+p1mup5*p3mup4*p4mup2*F(1056)+p3mup2*
     -   p3mup4*p3mup5*F(1306))+p2mup4*p3mup2*p3mup5*F(3216)+p2mup5*p
     -   3mup2*p3mup4*F(3223))+p1mup2*p2mup5*p3mup1*p3mup4*F(3224)+p1
     -   mup5*p2mup1*p3mup2*p3mup4*F(3225)+p1mup5*p2mup2*p3mup1*p3mup
     -   4*F(3226)+p1mup2*p2mup4*p3mup1*p3mup5*F(3227)-p1mup5*p2mup4*
     -   p3mup1*p3mup2*F(3229)+p1mup4*(p1mup5*p2mup1*p3mup2*F(1039)-p
     -   2mup1*p2mup5*p3mup2*F(1049)+p1mup2*(p2mup1*p3mup5*F(1036)-p2
     -   mup5*p3mup1*F(1037)+2*p3mup1*p3mup5*F(1050)+p3mup5*p4mup1*F(
     -   1200)-p3mup1*p4mup5*F(1201))-2*(p1mup5*p3mup1*p3mup2*F(1051)
     -   +p1mup5*p3mup1*p4mup2*F(1192)-p3mup1*p3mup2*p3mup5*F(1326))+
     -   p2mup2*(-(p2mup5*p3mup1*F(932))-p1mup5*p3mup1*F(1040)+p2mup1
     -   *p3mup5*F(1190)+p3mup1*p3mup5*F(3217))+p2mup1*p3mup2*p3mup5*
     -   F(3228)-p2mup5*p3mup1*p3mup2*F(3229))+2*(p1mup2*p1mup5*p3mup
     -   1*p3mup4*F(3222)+p1mup5*p3mup1*p3mup2*p3mup4*F(3230))
       Fa9(85) = p1mup2*p2mup1*(p2mup3*p3mup5*F(3234)+p2mup5*p3mup3*F(
     -   3235)-p2mup3*p4mup5*F(3238)-p2mup5*p4mup3*F(3239))+p1mup1*(4
     -   *p1mup2*p1mup3*p1mup5*F(3231)+p2mup2*(p2mup3*p3mup5*F(3232)+
     -   p2mup5*p3mup3*F(3233)-p2mup3*p4mup5*F(3236)-p2mup5*p4mup3*F(
     -   3237))+p2mup3*p2mup5*(p3mup2*F(3240)-p4mup2*F(3241)))
       Fa9(86) = p1mup1*(p2mup2*p2mup3*p4mup5*F(1160)-p2mup3*p2mup5*p4
     -   mup2*F(1161)-p2mup2*p2mup5*p4mup3*F(1188))+p1mup2*p2mup1*(p2
     -   mup3*p4mup5*F(1162)-p2mup5*p4mup3*F(1189))
       Fa9(87) = (p2mup5*p4mup3+p2mup3*p4mup5)*(p1mup1*p2mup2*F(1188)+
     -   p1mup2*p2mup1*F(1189))
       Fa9(88) = p1mup1*(-(p2mup2*p2mup3*p4mup5*F(1160))+p2mup3*p2mup5
     -   *p4mup2*F(1161)+p2mup2*p2mup5*p4mup3*F(1188))+p1mup2*p2mup1*
     -   (-(p2mup3*p4mup5*F(1162))+p2mup5*p4mup3*F(1189))
       Fa9(89) = p1mup1*p2mup2*p3mup5*p4mup3*F(3258)+p1mup1*p2mup2*p3m
     -   up3*p4mup5*F(3258)+p1mup1*p2mup3*p3mup2*p4mup5*F(3260)+p1mup
     -   1*p2mup3*p3mup5*p4mup2*F(3261)+p1mup2*p2mup1*p3mup5*p4mup3*F
     -   (3263)+p1mup2*p2mup1*p3mup3*p4mup5*F(3263)+p1mup2*p2mup3*p3m
     -   up1*p4mup5*F(3265)+p1mup2*p2mup3*p3mup5*p4mup1*F(3266)+p1mup
     -   1*p2mup5*p3mup2*p4mup3*F(3291)+p1mup1*p2mup5*p3mup3*p4mup2*F
     -   (3292)+p1mup2*p2mup5*p3mup1*p4mup3*F(3294)+p1mup2*p2mup5*p3m
     -   up3*p4mup1*F(3295)-2*(p1mup2*p2mup3*p2mup5*p3mup1*F(3250)+p1
     -   mup1*p2mup2*p3mup3*p3mup5*F(3251)+p1mup1*p2mup3*p3mup2*p3mup
     -   5*F(3252)+p1mup2*p2mup1*p3mup3*p3mup5*F(3253)+p1mup2*p2mup3*
     -   p3mup1*p3mup5*F(3254)-p1mup2*p2mup3*p2mup5*p4mup1*F(3257)-p1
     -   mup1*p2mup2*p4mup3*p4mup5*F(3259)-p1mup1*p2mup3*p4mup2*p4mup
     -   5*F(3262)-p1mup2*p2mup1*p4mup3*p4mup5*F(3264)-p1mup2*p2mup3*
     -   p4mup1*p4mup5*F(3267)+p1mup1*p2mup5*p3mup2*p3mup3*F(3282)+p1
     -   mup2*p2mup5*p3mup1*p3mup3*F(3283)-p1mup1*p2mup5*p4mup2*p4mup
     -   3*F(3293)-p1mup2*p2mup5*p4mup1*p4mup3*F(3296)-p1mup1*p3mup2*
     -   p3mup5*p4mup3*F(3309)-p1mup1*p3mup2*p3mup3*p4mup5*F(3309)-p1
     -   mup1*p3mup5*p4mup2*p4mup3*F(3312)-p1mup1*p3mup3*p4mup2*p4mup
     -   5*F(3312))-8*(p1mup2*p3mup1*p3mup3*p3mup5*F(3289)-p1mup2*p3m
     -   up1*p4mup3*p4mup5*F(3315)-p1mup2*p3mup3*p3mup5*p4mup1*F(3316
     -   )-p1mup2*p4mup1*p4mup3*p4mup5*F(3318))+p1mup3*(p2mup2*(p3mup
     -   1*p4mup5*F(3271)+p3mup5*p4mup1*F(3272)-2*(p2mup5*p3mup1*F(32
     -   43)-p2mup5*p4mup1*F(3247)+p3mup1*p3mup5*F(3256)-p4mup1*p4mup
     -   5*F(3273)))+p2mup1*(-(p3mup2*p4mup5*F(3268))-p3mup5*p4mup2*F
     -   (3269)+2*(p2mup5*p3mup2*F(3242)-p2mup5*p4mup2*F(3246)+p3mup2
     -   *p3mup5*F(3255)-p4mup2*p4mup5*F(3270)-p2mup2*p3mup5*F(3274)+
     -   p2mup2*p4mup5*F(3278)))+p2mup5*p3mup1*p4mup2*F(3297)+p2mup5*
     -   p3mup2*p4mup1*F(3298)-2*(p2mup5*p3mup1*p3mup2*F(3284)-p2mup5
     -   *p4mup1*p4mup2*F(3299)-p3mup1*p3mup2*p4mup5*F(3319)-p3mup1*p
     -   4mup2*p4mup5*F(3321)-p3mup2*p3mup5*p4mup1*F(3322)-p3mup5*p4m
     -   up1*p4mup2*F(3324))-4*(p3mup1*p3mup2*p3mup5*F(3290)-p3mup1*p
     -   3mup5*p4mup2*F(3320)-p3mup2*p4mup1*p4mup5*F(3323)-p4mup1*p4m
     -   up2*p4mup5*F(3325)))-4*(p2mup1*p2mup2*p2mup3*p3mup5*F(3276)+
     -   p2mup1*p2mup2*p2mup5*p3mup3*F(3277)-p2mup1*p2mup2*p2mup3*p4m
     -   up5*F(3280)-p2mup1*p2mup2*p2mup5*p4mup3*F(3281)+p1mup1*p3mup
     -   2*p3mup3*p3mup5*F(3288)-p1mup1*p3mup2*p4mup3*p4mup5*F(3310)-
     -   p1mup1*p3mup3*p3mup5*p4mup2*F(3311)-p1mup1*p4mup2*p4mup3*p4m
     -   up5*F(3313)-p1mup2*p3mup1*p3mup5*p4mup3*F(3314)-p1mup2*p3mup
     -   1*p3mup3*p4mup5*F(3314)-p1mup2*p3mup5*p4mup1*p4mup3*F(3317)-
     -   p1mup2*p3mup3*p4mup1*p4mup5*F(3317)+p2mup1*p2mup3*p2mup5*p3m
     -   up2*F(3326)-p2mup1*p2mup3*p2mup5*p4mup2*F(3327))+p1mup5*(p2m
     -   up1*(-(p3mup2*p4mup3*F(3300))+p3mup3*p4mup2*F(3301)-2*(p2mup
     -   3*p3mup2*F(3244)-p2mup3*p4mup2*F(3248)+p2mup2*p3mup3*F(3275)
     -   -p2mup2*p4mup3*F(3279)-p3mup2*p3mup3*F(3285)-p4mup2*p4mup3*F
     -   (3302)))+p2mup2*(p3mup1*p4mup3*F(3303)+p3mup3*p4mup1*F(3304)
     -   +2*(p2mup3*p3mup1*F(3245)-p2mup3*p4mup1*F(3249)-p3mup1*p3mup
     -   3*F(3286)+p4mup1*p4mup3*F(3305)))+p2mup3*p3mup1*p4mup2*F(330
     -   6)+p2mup3*p3mup2*p4mup1*F(3307)-2*(p2mup3*p3mup1*p3mup2*F(32
     -   87)-p2mup3*p4mup1*p4mup2*F(3308)-p3mup1*p3mup2*p4mup3*F(3329
     -   )-p3mup1*p4mup2*p4mup3*F(3331)-p3mup2*p3mup3*p4mup1*F(3332)-
     -   p3mup3*p4mup1*p4mup2*F(3334))-4*(p3mup1*p3mup2*p3mup3*F(3328
     -   )-p3mup1*p3mup3*p4mup2*F(3330)-p3mup2*p4mup1*p4mup3*F(3333)-
     -   p4mup1*p4mup2*p4mup3*F(3335)))
       Fa9(90) = -(p1mup1*p2mup2*p3mup5*p4mup3*F(1167))+p1mup1*p2mup2*
     -   p3mup3*p4mup5*F(1167)+p1mup1*p2mup3*p3mup2*p4mup5*F(1168)-p1
     -   mup1*p2mup5*p3mup3*p4mup2*F(1169)-p1mup1*p2mup3*p3mup5*p4mup
     -   2*F(1169)+p1mup1*p2mup5*p3mup2*p4mup3*F(1171)+p1mup2*p2mup3*
     -   p3mup1*p4mup5*F(1173)-p1mup2*p2mup5*p3mup3*p4mup1*F(1174)-p1
     -   mup2*p2mup3*p3mup5*p4mup1*F(1174)+p1mup2*p2mup5*p3mup1*p4mup
     -   3*F(1176)-p1mup5*p2mup2*p3mup3*p4mup1*F(1179)-p1mup5*p2mup3*
     -   p3mup2*p4mup1*F(1186)+p1mup5*p2mup3*p3mup1*p4mup2*F(1187)+p1
     -   mup5*p2mup2*p3mup1*p4mup3*F(1305)-4*(p1mup2*p3mup1*p3mup5*p4
     -   mup3*F(1307)-p1mup2*p3mup1*p3mup3*p4mup5*F(1307)+p1mup2*p3mu
     -   p5*p4mup1*p4mup3*F(1451)-p1mup2*p3mup3*p4mup1*p4mup5*F(1451)
     -   +p1mup5*p3mup2*p4mup1*p4mup3*F(1456)+p1mup5*p4mup1*p4mup2*p4
     -   mup3*F(1459))-2*(p1mup2*p2mup3*p2mup5*p4mup1*F(1163)+p1mup5*
     -   p2mup2*p2mup3*p4mup1*F(1165)+p1mup1*p3mup2*p3mup5*p4mup3*F(1
     -   306)-p1mup1*p3mup2*p3mup3*p4mup5*F(1306)+p1mup5*p3mup2*p3mup
     -   3*p4mup1*F(1309)-p1mup5*p3mup1*p3mup2*p4mup3*F(1311)+p1mup1*
     -   p3mup5*p4mup2*p4mup3*F(1450)-p1mup1*p3mup3*p4mup2*p4mup5*F(1
     -   450)+p1mup5*p3mup3*p4mup1*p4mup2*F(1453)-p1mup5*p3mup1*p4mup
     -   2*p4mup3*F(1455)+p1mup1*p2mup5*p4mup2*p4mup3*F(3336)-p1mup1*
     -   p2mup3*p4mup2*p4mup5*F(3336)-p1mup2*p2mup3*p4mup1*p4mup5*F(3
     -   337)+p1mup5*p2mup2*p4mup1*p4mup3*F(3338)+p1mup2*p2mup5*p4mup
     -   1*p4mup3*F(3338)+p1mup5*p2mup3*p4mup1*p4mup2*F(3340))+p1mup3
     -   *(-(p2mup5*p3mup1*p4mup2*F(1180))+p2mup5*p3mup2*p4mup1*F(118
     -   2)+4*(p3mup2*p4mup1*p4mup5*F(1456)+p4mup1*p4mup2*p4mup5*F(14
     -   59))+p2mup2*(p3mup1*p4mup5*F(1178)-p3mup5*p4mup1*F(1179)-2*(
     -   p2mup5*p4mup1*F(932)-p4mup1*p4mup5*F(3338)))+2*(p3mup1*p3mup
     -   2*p4mup5*F(1308)-p3mup2*p3mup5*p4mup1*F(1309)+p3mup1*p4mup2*
     -   p4mup5*F(1452)-p3mup5*p4mup1*p4mup2*F(1453)-p2mup5*p4mup1*p4
     -   mup2*F(3340)))+p2mup1*(-4*p2mup3*p2mup5*p4mup2*F(932)+p1mup2
     -   *(-(p3mup5*p4mup3)+p3mup3*p4mup5)*F(1172)+p1mup5*p3mup3*p4mu
     -   p2*F(1174)+p1mup3*p3mup5*p4mup2*F(1174)+p1mup3*p3mup2*p4mup5
     -   *F(1177)-p1mup5*p3mup2*p4mup3*F(1183)-2*(p1mup3*p2mup5*p4mup
     -   2*F(1049)-p1mup5*p2mup3*p4mup2*F(1166)-p2mup2*(2*(p2mup3*p4m
     -   up5*F(931)-p2mup5*p4mup3*F(933))+p1mup3*p4mup5*F(1190)-p1mup
     -   5*p4mup3*F(1193))-p1mup3*p4mup2*p4mup5*F(3339)-p1mup5*p4mup2
     -   *p4mup3*F(3341)))
       Return
       End
