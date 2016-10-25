c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenEvenGSim13(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p
     -   2mup1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup
     -   4,p3mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5
     -   mup3,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2m
     -   up3,mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa13)
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
       Complex*16 Fa13(97:103)
       COMMON/FaPenEvenGSim/Fa
       Fa13(97) = (mup1mup5*mup2mup3*F(3522))/2.+mup2mup5*p1mup3*p2mup
     -   1*F(3523)+mup1mup5*p1mup3*p2mup2*F(3524)+mup2mup3*p1mup5*p2m
     -   up1*F(3526)+mup3mup5*p1mup1*p3mup2*F(3529)+mup2mup5*p1mup1*p
     -   3mup3*F(3530)+mup2mup3*p1mup1*p3mup5*F(3531)+mup3mup5*p1mup2
     -   *p3mup1*F(3532)+mup1mup5*p1mup2*p3mup3*F(3533)+mup2mup5*p1mu
     -   p3*p3mup1*F(3535)+mup1mup5*p1mup3*p3mup2*F(3536)+mup2mup3*p1
     -   mup5*p3mup1*F(3538)+mup3mup5*p1mup1*p4mup2*F(3541)+mup2mup5*
     -   p1mup1*p4mup3*F(3542)+mup2mup3*p1mup1*p4mup5*F(3543)+mup3mup
     -   5*p1mup2*p4mup1*F(3544)+mup1mup5*p1mup2*p4mup3*F(3545)+mup2m
     -   up5*p1mup3*p4mup1*F(3547)+mup1mup5*p1mup3*p4mup2*F(3548)+mup
     -   2mup3*p1mup5*p4mup1*F(3550)+mup3mup5*p2mup1*p3mup2*F(3553)+m
     -   up2mup5*p2mup1*p3mup3*F(3554)+mup2mup3*p2mup1*p3mup5*F(3555)
     -   +mup3mup5*p2mup2*p3mup1*F(3556)+mup1mup5*p2mup2*p3mup3*F(355
     -   7)+mup2mup5*p2mup3*p3mup1*F(3559)+mup1mup5*p2mup3*p3mup2*F(3
     -   560)+mup2mup3*p2mup5*p3mup1*F(3562)+mup3mup5*p2mup1*p4mup2*F
     -   (3565)+mup2mup5*p2mup1*p4mup3*F(3566)+mup2mup3*p2mup1*p4mup5
     -   *F(3567)+mup3mup5*p2mup2*p4mup1*F(3568)+mup1mup5*p2mup2*p4mu
     -   p3*F(3569)+mup2mup5*p2mup3*p4mup1*F(3571)+mup1mup5*p2mup3*p4
     -   mup2*F(3572)+mup2mup3*p2mup5*p4mup1*F(3574)+mup1mup3*((mup2m
     -   up5*F(3521))/2.+p1mup5*p2mup2*F(3527)+p1mup2*p3mup5*F(3534)+
     -   p1mup5*p3mup2*F(3539)+p1mup2*p4mup5*F(3546)+p1mup5*p4mup2*F(
     -   3551)+p2mup2*p3mup5*F(3558)+p2mup5*p3mup2*F(3563)+p2mup2*p4m
     -   up5*F(3570)+p2mup5*p4mup2*F(3575))+mup1mup2*((mup3mup5*F(352
     -   0))/2.+p1mup3*p2mup5*F(3525)+p1mup5*p2mup3*F(3528)+p1mup3*p3
     -   mup5*F(3537)+p1mup5*p3mup3*F(3540)+p1mup3*p4mup5*F(3549)+p1m
     -   up5*p4mup3*F(3552)+p2mup3*p3mup5*F(3561)+p2mup5*p3mup3*F(356
     -   4)+p2mup3*p4mup5*F(3573)+p2mup5*p4mup3*F(3576))+4*(p1mup1*p1
     -   mup2*p2mup3*p3mup5*F(3584)+p1mup1*p1mup2*p2mup5*p3mup3*F(358
     -   5)+p1mup1*p1mup3*p2mup2*p3mup5*F(3586)-p1mup1*p1mup3*p2mup5*
     -   p3mup2*F(3587)-p1mup1*p1mup2*p2mup3*p4mup5*F(3590)-p1mup1*p1
     -   mup2*p2mup5*p4mup3*F(3591)-p1mup1*p1mup3*p2mup2*p4mup5*F(359
     -   2)+p1mup1*p1mup3*p2mup5*p4mup2*F(3593))+8*(p1mup1*p1mup2*p1m
     -   up3*p2mup5*F(3577)+p1mup1*p1mup3*p1mup5*p2mup2*F(3579)+p1mup
     -   1*p1mup2*p1mup3*p3mup5*F(3580)+p1mup1*p1mup2*p1mup5*p3mup3*F
     -   (3581)+p1mup1*p1mup2*p2mup3*p2mup5*F(3582)+p1mup1*p1mup5*p2m
     -   up2*p2mup3*F(3583)-p1mup1*p1mup2*p1mup3*p4mup5*F(3588)-p1mup
     -   1*p1mup2*p1mup5*p4mup3*F(3589)+p1mup2*p1mup3*p1mup5*p2mup1*F
     -   (3594)+p1mup2*p1mup5*p2mup1*p2mup3*F(3595)+p1mup1*p1mup3*p2m
     -   up2*p2mup5*F(3595)+p1mup2*p1mup3*p2mup1*p2mup5*F(3596)+p1mup
     -   2*p1mup3*p1mup5*p3mup1*F(3597)-p1mup2*p1mup3*p1mup5*p4mup1*F
     -   (3598)+p1mup3*p1mup5*p2mup1*p2mup2*F(3600)-p3mup1*p3mup2*p3m
     -   up5*p4mup3*F(3604)-p3mup1*p3mup2*p3mup3*p4mup5*F(3604)-p3mup
     -   2*p3mup5*p4mup1*p4mup3*F(3607)-p3mup1*p3mup5*p4mup2*p4mup3*F
     -   (3607)-p3mup2*p3mup3*p4mup1*p4mup5*F(3607)-p3mup1*p3mup3*p4m
     -   up2*p4mup5*F(3607)-p3mup5*p4mup1*p4mup2*p4mup3*F(3610)-p3mup
     -   3*p4mup1*p4mup2*p4mup5*F(3610))+16*(p1mup1*p1mup2*p1mup5*p2m
     -   up3*F(3578)+p2mup1*p2mup2*p2mup3*p2mup5*F(3599)+p1mup1*p1mup
     -   3*p1mup5*p3mup2*F(3601)-p1mup1*p1mup3*p1mup5*p4mup2*F(3602)+
     -   p3mup1*p3mup2*p3mup3*p3mup5*F(3603)-p3mup1*p3mup2*p4mup3*p4m
     -   up5*F(3605)-p3mup2*p3mup3*p3mup5*p4mup1*F(3606)-p3mup1*p3mup
     -   3*p3mup5*p4mup2*F(3606)-p3mup2*p4mup1*p4mup3*p4mup5*F(3608)-
     -   p3mup1*p4mup2*p4mup3*p4mup5*F(3608)-p3mup3*p3mup5*p4mup1*p4m
     -   up2*F(3609)-p4mup1*p4mup2*p4mup3*p4mup5*F(3611))
       Fa13(98) = mup1mup2*mup3mup5*F(3612)+mup1mup3*mup2mup5*F(3613)+
     -   mup1mup5*mup2mup3*F(3614)+2*(-8*(p3mup5*p4mup1*p4mup2*p4mup3
     -   *F(1848)-p3mup3*p4mup1*p4mup2*p4mup5*F(1848))+mup2mup5*p1mup
     -   3*p2mup1*F(3615)+mup1mup5*p1mup3*p2mup2*F(3616)+mup1mup2*p1m
     -   up3*p2mup5*F(3617)+mup2mup3*p1mup5*p2mup1*F(3618)+mup1mup3*p
     -   1mup5*p2mup2*F(3619)+mup1mup2*p1mup5*p2mup3*F(3620)-mup2mup5
     -   *p1mup3*p4mup1*F(3627)+mup1mup5*p1mup3*p4mup2*F(3628)+mup1mu
     -   p2*p1mup3*p4mup5*F(3629)+mup2mup3*p1mup5*p4mup1*F(3630)+mup1
     -   mup3*p1mup5*p4mup2*F(3631)+mup1mup2*p1mup5*p4mup3*F(3632)+mu
     -   p3mup5*p2mup1*p4mup2*F(3633)+mup2mup5*p2mup1*p4mup3*F(3634)+
     -   mup2mup3*p2mup1*p4mup5*F(3635)-mup3mup5*p2mup2*p4mup1*F(3636
     -   )+mup1mup5*p2mup2*p4mup3*F(3637)+mup1mup3*p2mup2*p4mup5*F(36
     -   38)+mup2mup5*p2mup3*p4mup1*F(3639)+mup1mup5*p2mup3*p4mup2*F(
     -   3640)+mup1mup2*p2mup3*p4mup5*F(3641)+mup2mup3*p2mup5*p4mup1*
     -   F(3642)+mup1mup3*p2mup5*p4mup2*F(3643)+mup1mup2*p2mup5*p4mup
     -   3*F(3644)+p1mup2*(-8*p1mup3*p1mup5*p4mup1*F(929)+mup3mup5*p4
     -   mup1*F(3624)+mup1mup5*p4mup3*F(3625)+mup1mup3*p4mup5*F(3626)
     -   +mup3mup5*p3mup1*F(3646)+(mup1mup5*p3mup3+mup1mup3*p3mup5)*F
     -   (3653))+mup2mup5*p2mup1*p3mup3*F(3655)+mup2mup3*p2mup1*p3mup
     -   5*F(3655)+mup1mup2*p2mup5*p3mup3*F(1174)*P(16)+mup1mup2*p2mu
     -   p3*p3mup5*F(1174)*P(16)+mup1mup5*p2mup2*p3mup3*F(1179)*P(16)
     -   +mup1mup3*p2mup2*p3mup5*F(1179)*P(16)+p3mup2*(8*p4mup1*(-(p3
     -   mup5*p4mup3)+p3mup3*p4mup5)*F(1700)+mup1mup5*p1mup3*F(3647)+
     -   mup1mup3*p1mup5*F(3648)+mup3mup5*p2mup1*F(3652)+mup1mup3*p2m
     -   up5*F(3656)+mup1mup5*p2mup3*F(1186)*P(16))+p1mup1*(-4*(p1mup
     -   3*p2mup5*p4mup2*F(929)-p1mup3*p2mup2*p4mup5*F(1035)-p1mup2*(
     -   2*(p1mup3*p4mup5*F(928)-p1mup5*p4mup3*F(930))+p2mup3*p4mup5*
     -   F(1034)-p2mup5*p4mup3*F(1045)))+mup3mup5*p4mup2*F(3621)-mup2
     -   mup5*p4mup3*F(3622)+mup2mup3*p4mup5*F(3623)+mup3mup5*p3mup2*
     -   F(3645)+(mup2mup5*p3mup3+mup2mup3*p3mup5)*F(1169)*P(45))+p3m
     -   up1*(-8*(p3mup5*p4mup3-p3mup3*p4mup5)*(p3mup2*F(1054)+p4mup2
     -   *F(1700))+mup3mup5*p2mup2*F(3649)+mup2mup5*p2mup3*F(3650)+mu
     -   p2mup3*p2mup5*F(3651)+mup2mup3*p1mup5*F(3654)+mup2mup5*p1mup
     -   3*F(1180)*P(45))+mup1mup2*(p1mup5*p3mup3+p1mup3*p3mup5)*F(11
     -   74)*P(247))
       Fa13(99) = mup1mup3*mup2mup5*F(3658)+mup1mup5*mup2mup3*F(3659)+
     -   2*(mup2mup5*p1mup3*p2mup1*F(3660)+mup1mup5*p1mup3*p2mup2*F(3
     -   661)+mup2mup3*p1mup5*p2mup1*F(3663)+mup1mup3*p1mup5*p2mup2*F
     -   (3664)+mup1mup3*p1mup2*p4mup5*F(3671)+mup2mup5*p1mup3*p4mup1
     -   *F(3672)+mup2mup3*p1mup5*p4mup1*F(3675)+mup1mup3*p1mup5*p4mu
     -   p2*F(3676)+mup2mup5*p2mup1*p4mup3*F(3679)+mup2mup3*p2mup1*p4
     -   mup5*F(3680)+mup1mup3*p2mup2*p4mup5*F(3683)+mup2mup5*p2mup3*
     -   p4mup1*F(3684)+mup2mup3*p2mup5*p4mup1*F(3687)+mup1mup3*p2mup
     -   5*p4mup2*F(3688)+mup1mup3*p1mup5*p3mup2*F(3693)+mup2mup5*p2m
     -   up3*p3mup1*F(3695)+mup2mup3*p2mup5*p3mup1*F(3696)+mup3mup5*(
     -   p1mup2*p4mup1*F(3669)+p2mup1*p4mup2*F(3678)+p2mup2*p4mup1*F(
     -   3681)+p1mup2*p3mup1*F(3691)+p2mup2*p3mup1*F(3694)+p2mup1*p3m
     -   up2*F(3697))+mup1mup3*p1mup2*p3mup5*F(3699)+mup2mup5*p2mup1*
     -   p3mup3*F(3700)+mup2mup3*p2mup1*p3mup5*F(3701)+8*(p3mup2*p3mu
     -   p5*p4mup1*p4mup3*F(3799)+p3mup1*p3mup5*p4mup2*p4mup3*F(3799)
     -   +p3mup2*p3mup3*p4mup1*p4mup5*F(3799)+p3mup1*p3mup3*p4mup2*p4
     -   mup5*F(3799)+p3mup5*p4mup1*p4mup2*p4mup3*F(3802)+p3mup3*p4mu
     -   p1*p4mup2*p4mup5*F(3802)+p3mup1*p3mup2*p3mup5*p4mup3*F(3805)
     -   +p3mup1*p3mup2*p3mup3*p4mup5*F(3805))+16*(p3mup1*p3mup2*p4mu
     -   p3*p4mup5*F(3799)+p3mup2*p4mup1*p4mup3*p4mup5*F(3802)+p3mup1
     -   *p4mup2*p4mup3*p4mup5*F(3802)+p4mup1*p4mup2*p4mup3*p4mup5*F(
     -   3806))+mup1mup5*(p1mup2*p4mup3*F(3670)+p1mup3*p4mup2*F(3673)
     -   +p2mup2*p4mup3*F(3682)+p2mup3*p4mup2*F(3685)+p1mup3*p3mup2*F
     -   (3692)+p1mup2*p3mup3*F(3698)+(p2mup3*p3mup2+p2mup2*p3mup3)*F
     -   (1199)*P(16))+mup2mup5*p1mup3*p3mup1*F(1195)*P(45)+mup2mup3*
     -   p1mup5*p3mup1*F(1195)*P(247)+p1mup1*(8*p1mup2*(p1mup5*p4mup3
     -   +p1mup3*p4mup5)*F(930)+4*(p1mup2*(p2mup5*p4mup3+p2mup3*p4mup
     -   5)*F(1045)+p1mup3*p2mup2*p4mup5*F(1046)-p1mup3*p2mup5*p4mup2
     -   *F(1047))+mup3mup5*p4mup2*F(3666)+mup2mup5*p4mup3*F(3667)+mu
     -   p2mup3*p4mup5*F(3668)+mup3mup5*p3mup2*F(3690)+mup2mup5*p3mup
     -   3*F(1195)*P(45)+mup2mup3*p3mup5*F(1195)*P(247))+mup1mup3*(p2
     -   mup5*p3mup2+p2mup2*p3mup5)*F(1199)*P(263))+mup1mup2*(mup3mup
     -   5*F(3657)+2*(p1mup5*p2mup3*F(3665)+p1mup5*p4mup3*F(3677)+p2m
     -   up3*p4mup5*F(3686)+p2mup5*p4mup3*F(3689)+p2mup5*p3mup3*F(379
     -   0)*P(16)+p1mup3*(p2mup5*F(3662)+p4mup5*F(3674)+p3mup5*F(3790
     -   )*P(45))+p1mup5*p3mup3*F(3790)*P(247)+p2mup3*p3mup5*F(3790)*
     -   P(263)))
       Fa13(100) = mup1mup2*mup3mup5*F(3702)+mup1mup3*mup2mup5*F(3703)
     -   +mup1mup5*mup2mup3*F(3704)+2*(8*(p3mup5*p4mup1*p4mup2*p4mup3
     -   *F(1848)-p3mup3*p4mup1*p4mup2*p4mup5*F(1848))+mup2mup5*p1mup
     -   3*p2mup1*F(3705)+mup1mup5*p1mup3*p2mup2*F(3706)+mup1mup2*p1m
     -   up3*p2mup5*F(3707)+mup2mup3*p1mup5*p2mup1*F(3708)+mup1mup3*p
     -   1mup5*p2mup2*F(3709)+mup1mup2*p1mup5*p2mup3*F(3710)+mup2mup5
     -   *p1mup3*p4mup1*F(3717)-mup1mup5*p1mup3*p4mup2*F(3718)-mup1mu
     -   p2*p1mup3*p4mup5*F(3719)+mup2mup3*p1mup5*p4mup1*F(3720)+mup1
     -   mup3*p1mup5*p4mup2*F(3721)+mup1mup2*p1mup5*p4mup3*F(3722)+mu
     -   p3mup5*p2mup1*p4mup2*F(3723)+mup2mup5*p2mup1*p4mup3*F(3724)+
     -   mup2mup3*p2mup1*p4mup5*F(3725)+mup3mup5*p2mup2*p4mup1*F(3726
     -   )+mup1mup5*p2mup2*p4mup3*F(3727)-mup1mup3*p2mup2*p4mup5*F(37
     -   28)+mup2mup5*p2mup3*p4mup1*F(3729)+mup1mup5*p2mup3*p4mup2*F(
     -   3730)+mup1mup2*p2mup3*p4mup5*F(3731)+mup2mup3*p2mup5*p4mup1*
     -   F(3732)+mup1mup3*p2mup5*p4mup2*F(3733)+mup1mup2*p2mup5*p4mup
     -   3*F(3734)+p1mup2*(8*p1mup3*p1mup5*p4mup1*F(929)+mup3mup5*p4m
     -   up1*F(3714)+mup1mup5*p4mup3*F(3715)+mup1mup3*p4mup5*F(3716)+
     -   mup3mup5*p3mup1*F(3736)+(mup1mup5*p3mup3+mup1mup3*p3mup5)*F(
     -   3743))+mup2mup5*p2mup1*p3mup3*F(3745)+mup2mup3*p2mup1*p3mup5
     -   *F(3745)+mup1mup2*p1mup5*p3mup3*F(1174)*P(45)+mup1mup2*p1mup
     -   3*p3mup5*F(1174)*P(45)+p1mup1*(4*(p1mup3*p2mup5*p4mup2*F(929
     -   )-p1mup3*p2mup2*p4mup5*F(1035))+p1mup2*(-8*(p1mup3*p4mup5*F(
     -   928)-p1mup5*p4mup3*F(930))-4*(p2mup3*p4mup5*F(1034)-p2mup5*p
     -   4mup3*F(1045)))-mup3mup5*p4mup2*F(3711)+mup2mup5*p4mup3*F(37
     -   12)+mup2mup3*p4mup5*F(3713)+mup3mup5*p3mup2*F(3735)+(mup2mup
     -   5*p3mup3+mup2mup3*p3mup5)*F(1169)*P(247))+p3mup1*(8*(p3mup5*
     -   p4mup3-p3mup3*p4mup5)*(p3mup2*F(1054)+p4mup2*F(1700))+mup3mu
     -   p5*p2mup2*F(3739)+mup2mup5*p2mup3*F(3740)+mup2mup3*p2mup5*F(
     -   3741)+mup2mup3*p1mup5*F(3744)+mup2mup5*p1mup3*F(1180)*P(247)
     -   )+(mup1mup2*(p2mup5*p3mup3+p2mup3*p3mup5)*F(1174)+p2mup2*(mu
     -   p1mup5*p3mup3+mup1mup3*p3mup5)*F(1179))*P(263)+p3mup2*(8*p4m
     -   up1*(p3mup5*p4mup3-p3mup3*p4mup5)*F(1700)+mup1mup5*p1mup3*F(
     -   3737)+mup1mup3*p1mup5*F(3738)+mup3mup5*p2mup1*F(3742)+mup1mu
     -   p3*p2mup5*F(3746)+mup1mup5*p2mup3*F(1186)*P(263)))
       Fa13(101) = -(p1mup3*p1mup5*p2mup1*p3mup2*F(3761))+p1mup3*p1mup
     -   5*p2mup2*p3mup1*F(3762)+p1mup3*p1mup5*p2mup1*p4mup2*F(3765)-
     -   p1mup3*p1mup5*p2mup2*p4mup1*F(3766)+4*p1mup3*p2mup1*p2mup2*p
     -   2mup5*F(3772)+p1mup2*(2*p2mup1*p2mup3*p2mup5*F(3756)+p1mup3*
     -   (p2mup1*p3mup5*F(3751)+p2mup5*p3mup1*F(3752)-p2mup1*p4mup5*F
     -   (3757)-p2mup5*p4mup1*F(3758)-p3mup1*p4mup5*F(3769)-p3mup5*p4
     -   mup1*F(3770)+2*(p3mup1*p3mup5*F(3764)-p4mup1*p4mup5*F(3771))
     -   )+p1mup5*(p2mup1*p3mup3*F(3753)-p2mup3*p3mup1*F(3754)-p2mup1
     -   *p4mup3*F(3759)+p2mup3*p4mup1*F(3760)+2*(p3mup1*p3mup3*F(377
     -   4)-p4mup1*p4mup3*F(3778))-p3mup1*p4mup3*F(3779)-p3mup3*p4mup
     -   1*F(3780)))+2*(p1mup5*p2mup1*p2mup2*p2mup3*F(3773)+p1mup3*p1
     -   mup5*p3mup1*p3mup2*F(3775)-p1mup3*p1mup5*p3mup2*p4mup1*F(378
     -   1)-p1mup3*p1mup5*p3mup1*p4mup2*F(3781)-p1mup3*p1mup5*p4mup1*
     -   p4mup2*F(3782))+p1mup1*(p1mup2*(-(p3mup5*p4mup3*F(3767))-p3m
     -   up3*p4mup5*F(3767)+2*(p3mup3*p3mup5*F(3763)-p4mup3*p4mup5*F(
     -   3768)))+2*(p2mup2*p2mup3*p2mup5*F(3755)+p1mup3*(p3mup2*(p3mu
     -   p5*F(3774)-p4mup5*F(3776))-p4mup2*(p3mup5*F(3777)+p4mup5*F(3
     -   778))))+p1mup5*(p2mup2*p3mup3*F(3747)+p2mup3*p3mup2*F(3748)-
     -   p2mup2*p4mup3*F(3749)-p2mup3*p4mup2*F(3750)+2*(p3mup2*p3mup3
     -   *F(3783)-p3mup2*p4mup3*F(3784)-p3mup3*p4mup2*F(3785)-p4mup2*
     -   p4mup3*F(3786))))
       Fa13(102) = p1mup2*(p1mup1*(p3mup5*p4mup3-p3mup3*p4mup5)*F(1041
     -   )+p1mup5*(p2mup3*p4mup1*F(1038)+p3mup3*p4mup1*F(1043)-p3mup1
     -   *p4mup3*F(1044)+p2mup1*p4mup3*F(1048)+2*p4mup1*p4mup3*F(1053
     -   ))+p1mup3*(p4mup1*(p2mup5*F(1037)+p3mup5*F(1043))-p4mup5*(p2
     -   mup1*F(1036)+p3mup1*F(1042)+2*p4mup1*F(1053))))+p1mup1*p1mup
     -   5*(p2mup2*p4mup3*F(1046)+p2mup3*p4mup2*F(1047)+2*(p3mup2*p4m
     -   up3*F(1055)+p4mup2*p4mup3*F(1056)))+p1mup3*(-2*p1mup1*p4mup5
     -   *(p3mup2*F(1050)+p4mup2*F(1053))+p1mup5*(-(p2mup1*p4mup2*F(1
     -   039))+p2mup2*p4mup1*F(1040)+2*(p3mup2*p4mup1*F(1051)+p4mup1*
     -   p4mup2*F(1192))))
       Fa13(103) = p1mup1*(2*p4mup5*(p1mup2*p4mup3*F(1052)+p1mup3*p3mu
     -   p2*F(1055)+p1mup3*p4mup2*F(1056))+p1mup5*(p2mup2*p4mup3*F(10
     -   46)+p2mup3*p4mup2*F(1047)+2*(p3mup2*p4mup3*F(1055)+p4mup2*p4
     -   mup3*F(1056)))+p1mup2*(p3mup5*p4mup3+p3mup3*p4mup5)*F(1057))
     -   +p1mup2*(p1mup3*(p2mup1*p4mup5*F(1048)+p2mup5*p4mup1*F(1049)
     -   +p3mup5*p4mup1*F(1051)+2*p4mup1*p4mup5*F(1053)+p3mup1*p4mup5
     -   *F(1057))+p1mup5*(-(p4mup1*(p2mup3*F(1049)+p3mup3*F(1051)))+
     -   p4mup3*(p2mup1*F(1048)+2*p4mup1*F(1056)+p3mup1*F(1057))))
       Return
       End
