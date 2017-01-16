c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenOddGSim18(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4
     -   ,p3mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5m
     -   up3,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mu
     -   p3,mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa18)
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
       Complex*16 Fa18(1684:1782)
       COMMON/FaPenOddGSim/Fa
       Fa18(1684) = p1mup4*F(3691)+p2mup4*F(3759)+2*p3mup4*F(3760)+p4m
     -   up4*F(3761)
       Fa18(1685) = p1mup4*F(3693)+p2mup4*F(3762)+p3mup4*F(3763)+2*p4m
     -   up4*F(3764)
       Fa18(1686) = p3mup1*F(3661)+p4mup1*F(3662)-p2mup1*F(3722)
       Fa18(1687) = p3mup1*F(3723)+p4mup1*F(3724)
       Fa18(1688) = p4mup1*F(3663)+p2mup1*F(3720)
       Fa18(1689) = p3mup1*F(3664)+p2mup1*F(3721)
       Fa18(1690) = 2*p1mup2*F(3732)+p3mup2*F(3766)+p4mup2*F(3767)
       Fa18(1691) = p1mup2*F(3735)+p3mup2*F(3768)+p4mup2*F(3769)
       Fa18(1692) = p3mup2*F(3771)+M1*M3*(p1mup2*F(3770)+p4mup2*F(3772
     -   ))
       Fa18(1693) = M1*M3*(p1mup2*F(3773)+p3mup2*F(3774))+p4mup2*F(377
     -   5)
       Fa18(1694) = 2*p1mup3*F(3688)+p2mup3*F(3752)
       Fa18(1695) = M1*M3*p1mup3*F(3755)+p2mup3*F(3756)
       Fa18(1696) = p1mup3*F(3691)+p2mup3*F(3759)
       Fa18(1697) = p1mup3*F(3693)+p2mup3*F(3762)
       Fa18(1698) = 2*p1mup4*F(3696)+p3mup4*F(3697)+p4mup4*F(3698)+p2m
     -   up4*F(3747)
       Fa18(1699) = p1mup4*F(3748)+2*p2mup4*F(3749)+p3mup4*F(3779)+p4m
     -   up4*F(3780)
       Fa18(1700) = p1mup4*F(3703)+2*p3mup4*F(3705)+p4mup4*F(3706)+p2m
     -   up4*F(3750)
       Fa18(1701) = p1mup4*F(3707)+p3mup4*F(3709)+2*p4mup4*F(3710)+p2m
     -   up4*F(3751)
       Fa18(1702) = p2mup1*F(3658)+p3mup1*F(3717)+p4mup1*F(3718)
       Fa18(1703) = p3mup1*F(3659)+p4mup1*F(3660)
       Fa18(1704) = p2mup1*F(3656)+p4mup1*F(3664)
       Fa18(1705) = p2mup1*F(3657)+p3mup1*F(3663)
       Fa18(1706) = 2*p1mup2*F(3672)+p3mup2*F(3782)+p4mup2*F(3783)
       Fa18(1707) = p1mup2*F(3784)+p3mup2*F(3785)+p4mup2*F(3786)
       Fa18(1708) = p1mup2*F(3787)-2*p3mup2*F(3788)+p4mup2*F(3789)
       Fa18(1709) = p1mup2*F(3790)+p3mup2*F(3791)-2*p4mup2*F(3792)
       Fa18(1710) = p2mup3*F(3695)+2*p1mup3*F(3696)
       Fa18(1711) = p1mup3*F(3699)-2*p2mup3*F(3700)
       Fa18(1712) = p1mup3*F(3703)+p2mup3*F(3704)
       Fa18(1713) = p1mup3*F(3707)+p2mup3*F(3708)
       Fa18(1714) = -(p2mup4*F(3687))+2*p1mup4*F(3688)+p3mup4*F(3753)+
     -   p4mup4*F(3754)
       Fa18(1715) = p1mup4*F(3689)+2*p2mup4*F(3690)+p3mup4*F(3796)+p4m
     -   up4*F(3797)
       Fa18(1716) = p1mup4*F(3691)+p2mup4*F(3692)+2*p3mup4*F(3760)+p4m
     -   up4*F(3761)
       Fa18(1717) = p1mup4*F(3693)+p2mup4*F(3694)+p3mup4*F(3763)+2*p4m
     -   up4*F(3764)
       Fa18(1718) = p2mup1*F(3597)+p3mup1*F(3653)+p4mup1*F(3654)
       Fa18(1719) = p3mup1*F(3598)+p4mup1*F(3599)
       Fa18(1720) = p2mup1*F(3595)+p4mup1*F(3603)
       Fa18(1721) = p2mup1*F(3596)+p3mup1*F(3602)
       Fa18(1722) = 2*p1mup2*F(3803)+p3mup2*F(3804)+p4mup2*F(3805)
       Fa18(1723) = p1mup2*F(3806)+p3mup2*F(3807)+p4mup2*F(3808)
       Fa18(1724) = p1mup2*F(3809)+p4mup2*F(3810)-2*p3mup2*F(3811)
       Fa18(1725) = -(p1mup2*F(3812))+p3mup2*F(3813)-2*p4mup2*F(3814)
       Fa18(1726) = p2mup3*F(3815)+2*p1mup3*F(3816)
       Fa18(1727) = 2*p2mup3*F(3817)+p1mup3*F(3818)
       Fa18(1728) = p2mup3*F(3819)+p1mup3*F(3820)
       Fa18(1729) = -(p2mup3*F(3821))+p1mup3*F(3822)
       Fa18(1730) = p2mup4*F(3823)+2*p1mup4*F(3824)+p3mup4*F(3825)+p4m
     -   up4*F(3826)
       Fa18(1731) = p1mup4*F(3827)+p3mup4*F(3828)+2*p2mup4*F(3829)+p4m
     -   up4*F(3830)
       Fa18(1732) = p2mup4*F(3831)+p1mup4*F(3832)+2*p3mup4*F(3833)+p4m
     -   up4*F(3834)
       Fa18(1733) = p2mup4*F(3835)+p1mup4*F(3836)+p3mup4*F(3837)+2*p4m
     -   up4*F(3838)
       Fa18(1734) = p2mup1*F(3533)+p3mup1*F(3592)+p4mup1*F(3593)
       Fa18(1735) = p3mup1*F(3534)+p4mup1*F(3535)
       Fa18(1736) = p2mup1*F(3531)+p4mup1*F(3539)
       Fa18(1737) = p2mup1*F(3532)+p3mup1*F(3538)
       Fa18(1738) = p1mup2*F(3850)+p3mup2*F(3851)+p4mup2*F(3852)
       Fa18(1739) = p1mup2*F(3853)+p3mup2*F(3854)+p4mup2*F(3855)
       Fa18(1740) = p1mup2*F(3856)+p4mup2*F(3857)+p3mup2*F(3858)
       Fa18(1741) = p3mup2*F(3859)+p1mup2*F(3860)+2*p4mup2*F(3861)
       Fa18(1742) = p2mup3*F(3865)+p1mup3*F(3866)
       Fa18(1743) = p1mup3*F(3867)+p2mup3*F(3868)
       Fa18(1744) = p1mup3*F(3869)+p2mup3*F(3870)
       Fa18(1745) = p2mup3*F(3871)+p1mup3*F(3872)
       Fa18(1746) = p2mup4*F(3873)+p3mup4*F(3874)+p1mup4*F(3875)+p4mup
     -   4*F(3876)
       Fa18(1747) = p1mup4*F(3877)+p3mup4*F(3878)+p2mup4*F(3879)+p4mup
     -   4*F(3880)
       Fa18(1748) = p1mup4*F(3881)+p2mup4*F(3882)+p3mup4*F(3883)+p4mup
     -   4*F(3884)
       Fa18(1749) = p2mup4*F(3885)+p3mup4*F(3886)+2*p4mup4*F(3887)+p1m
     -   up4*F(3888)
       Fa18(1750) = 8*p2mup3*p2mup4*F(3899)+p2mup3*p3mup4*F(3901)+p2mu
     -   p4*p3mup3*F(3902)+4*p3mup3*p3mup4*F(3903)+p1mup3*(p1mup4*F(3
     -   895)+p3mup4*F(3897)-8*p2mup4*F(3900)+p4mup4*F(3904))+p2mup3*
     -   p4mup4*F(3905)+p3mup3*p4mup4*F(3906)+p3mup4*p4mup3*F(3907)+p
     -   1mup4*(8*p2mup3*F(3896)+p3mup3*F(3898)+p4mup3*F(3908))+p2mup
     -   4*p4mup3*F(3909)+p4mup3*p4mup4*F(3910)
       Fa18(1751) = p3mup3*p3mup4*F(3914)+8*(p2mup3*p2mup4*F(3913)+p2m
     -   up3*p3mup4*F(3915)+p2mup4*p3mup3*F(3916))+p3mup3*p4mup4*F(39
     -   18)+p2mup4*p4mup3*F(3921)+p3mup4*p4mup3*F(3922)+4*(p2mup3*p4
     -   mup4*F(3917)+p4mup3*p4mup4*F(3924))+p1mup3*(p2mup4*F(3911)+4
     -   *p1mup4*F(3919)+p3mup4*F(3920)+p4mup4*F(3925))+p1mup4*(p2mup
     -   3*F(3912)+p3mup3*F(3923)+p4mup3*F(3926))
       Fa18(1752) = p3mup4*p4mup3*F(3933)+8*(p4mup3*p4mup4*F(3934)+p2m
     -   up4*p4mup3*F(3935))+p2mup3*(p3mup4*F(3928)+8*(p2mup4*F(3929)
     -   +p4mup4*F(3931))+p1mup4*F(3937))+p1mup3*p3mup4*F(3938)+p3mup
     -   3*(4*(p3mup4*F(3927)+p2mup4*F(3930)+p4mup4*F(3932))+p1mup4*F
     -   (3939))+4*(p1mup3*p1mup4*F(3936)+p1mup3*p2mup4*F(3940)+p1mup
     -   3*p4mup4*F(3941))+p1mup4*p4mup3*F(3942)
       Fa18(1753) = p3mup3*F(3943)+p4mup3*F(3944)
       Fa18(1754) = p3mup3*F(3945)+p4mup3*F(3946)
       Fa18(1755) = p4mup3*F(3947)+2*p3mup3*F(3948)
       Fa18(1756) = p3mup3*F(3949)+2*p4mup3*F(3950)
       Fa18(1757) = p2mup3*p3mup4*F(3954)-4*(p1mup3*p2mup4*F(3951)-p1m
     -   up4*p2mup3*F(3952)+p2mup3*p2mup4*F(3953)-p2mup4*p3mup3*F(395
     -   5))+p2mup3*p4mup4*F(3956)+2*p2mup4*p4mup3*F(3957)
       Fa18(1758) = -2*p2mup4*(2*(p3mup3*F(3959)+p4mup3*F(3962))+p1mup
     -   3*F(3963))+p2mup3*(-4*(p2mup4*F(3958)-p3mup4*F(3960)-p4mup4*
     -   F(3961))+2*p1mup4*F(3963))
       Fa18(1759) = p2mup3*p3mup4*F(3966)+p2mup4*p3mup3*F(3967)-4*(p2m
     -   up3*p2mup4*F(3964)+p1mup3*p2mup4*F(3965)+p2mup4*p4mup3*F(396
     -   8)-m1sq*p2mup3*(p1mup4*F(346)+p4mup4*F(3969)))
       Fa18(1760) = p3mup3*p4mup4*F(3979)+p3mup4*p4mup3*F(3980)+2*(p3m
     -   up3*p3mup4*F(3971)+p1mup4*p3mup3*F(3974)+p1mup3*(2*(p1mup4+p
     -   2mup4)*F(376)+p3mup4*F(3973)+p4mup4*F(3975))+p2mup4*p3mup3*F
     -   (3977)+p2mup3*(2*(p2mup4*F(3970)+p1mup4*F(3972))+p3mup4*F(39
     -   76)+p4mup4*F(3978))-p1mup4*p4mup3*F(3981)-p2mup4*p4mup3*F(39
     -   82)+p4mup3*p4mup4*F(3983))
       Fa18(1761) = p1mup3*p4mup4*F(3986)-p1mup4*p4mup3*F(3987)+2*(p1m
     -   up3*p1mup4*F(451)+p1mup3*p3mup4*F(3985)+p1mup4*p2mup3*F(3988
     -   )+p2mup3*p4mup4*F(3990)+p2mup4*(2*(p3mup3*F(382)+p2mup3*F(45
     -   6))+p1mup3*F(3984)-p4mup3*F(3991))+p3mup3*(2*p3mup4*F(382)+p
     -   1mup4*F(3992)+p4mup4*F(3993))-p3mup4*p4mup3*F(3994))+4*(p2mu
     -   p3*p3mup4*F(3989)+p4mup3*p4mup4*F(4467))
       Fa18(1762) = p1mup3*(p3mup4*F(3996)+2*(p1mup4*F(451)+p2mup4*F(3
     -   998)-p4mup4*F(4001)))+p3mup3*(4*p3mup4*F(382)+p1mup4*F(3997)
     -   +2*(p2mup4*F(4000)+p4mup4*F(4004)))+2*(p2mup3*(p1mup4*F(3995
     -   )+p3mup4*F(3999)+2*(p2mup4*F(456)+p4mup4*F(4003)))+p4mup3*(p
     -   1mup4*F(4002)+p3mup4*F(4005)+2*(p2mup4+p4mup4)*F(4467)))
       Fa18(1763) = p3mup3*F(4006)+p4mup3*F(4007)
       Fa18(1764) = p3mup3*F(4008)+p4mup3*F(4009)
       Fa18(1765) = p4mup3*F(4010)-2*p3mup3*F(4011)
       Fa18(1766) = p3mup3*F(4012)+2*p4mup3*F(4013)
       Fa18(1767) = 2*p1mup3*p2mup4*F(258)+p2mup3*(2*(p2mup4*F(449)+p1
     -   mup4*F(4014)+p3mup4*F(4015))+p4mup4*F(4016))+p2mup4*(2*p3mup
     -   3*F(600)-p4mup3*F(4017))
       Fa18(1768) = -(p1mup4*p2mup3*F(4018))+p1mup3*p2mup4*F(4018)+2*(
     -   (p2mup4*p3mup3-p2mup3*p3mup4)*F(4019)+(p2mup4*p4mup3-p2mup3*
     -   p4mup4)*F(4020))
       Fa18(1769) = p2mup3*(p3mup4*F(4022)+2*(p4mup4*F(362)+p2mup4*F(4
     -   49)+p1mup4*F(4024)))+p2mup4*(p3mup3*F(4023)+2*(p1mup3*F(4021
     -   )+p4mup3*F(4025)))
       Fa18(1770) = 2*(p1mup3*(p3mup4*F(258)+p4mup4*F(346)+2*p2mup4*F(
     -   1262))+p1mup4*(p3mup3*F(264)+2*p2mup3*F(1608)+p4mup3*F(4024)
     -   )+p2mup3*p3mup4*F(4026)+p2mup3*p4mup4*F(4027)-p2mup4*(p3mup3
     -   *F(4026)+p4mup3*F(4027)))+p3mup3*p4mup4*F(4028)+p3mup4*p4mup
     -   3*F(4029)
       Fa18(1771) = p1mup3*(2*(p3mup4*F(282)+p2mup4*F(4026))+p4mup4*F(
     -   4030))+p1mup4*(2*(p3mup3*F(600)-p2mup3*F(4026))+p4mup3*F(403
     -   1))+2*(p2mup3*p4mup4*F(4032)+p2mup4*p4mup3*F(4033)+p3mup3*p4
     -   mup4*F(4034)+p3mup4*p4mup3*F(4035)+2*(p2mup3*p3mup4*F(4036)+
     -   p2mup4*p3mup3*F(4037)))
       Fa18(1772) = p1mup3*(2*(p4mup4*F(362)+p2mup4*F(4027))+p3mup4*F(
     -   4038))+p1mup4*(2*(p4mup3*F(3969)-p2mup3*F(4027))+p3mup3*F(40
     -   39))+2*(p2mup3*p3mup4*F(4032)+p2mup4*p3mup3*F(4033)+p3mup3*p
     -   4mup4*F(4040)+p3mup4*p4mup3*F(4041)+2*(p2mup3*p4mup4*F(4042)
     -   +p2mup4*p4mup3*F(4043)))
       Fa18(1773) = p3mup3*F(3753)+p4mup3*F(3754)
       Fa18(1774) = p3mup3*F(3796)+p4mup3*F(3797)
       Fa18(1775) = 2*p3mup3*F(3760)+p4mup3*F(3761)
       Fa18(1776) = p3mup3*F(3763)+2*p4mup3*F(3764)
       Fa18(1777) = 2*(p1mup3*p2mup4*F(258)+p1mup4*p2mup3*F(264)+p2mup
     -   3*p3mup4*F(282)+p2mup4*p3mup3*F(600))+p2mup3*p4mup4*F(4038)+
     -   p2mup4*p4mup3*F(4039)
       Fa18(1778) = p1mup3*p2mup4*F(4028)+p1mup4*p2mup3*F(4029)+2*(p2m
     -   up3*p3mup4*F(4034)+p2mup4*p3mup3*F(4035)+p2mup3*p4mup4*F(404
     -   0)+p2mup4*p4mup3*F(4041))
       Fa18(1779) = 2*(p1mup3*p2mup4*F(346)+p2mup3*p4mup4*F(362)+p2mup
     -   4*p4mup3*F(3969)+p1mup4*p2mup3*F(4024))+p2mup3*p3mup4*F(4030
     -   )+p2mup4*p3mup3*F(4031)
       Fa18(1780) = 2*(p1mup4*(p3mup3*F(258)+p4mup3*F(346)+2*p2mup3*F(
     -   1262))+p1mup3*(p3mup4*F(264)+2*p2mup4*F(1608)+p4mup4*F(4024)
     -   )+p2mup4*p3mup3*F(4026)-p2mup3*p3mup4*F(4026)+p2mup4*p4mup3*
     -   F(4027)-p2mup3*p4mup4*F(4027))+p3mup4*p4mup3*F(4028)+p3mup3*
     -   p4mup4*F(4029)
       Fa18(1781) = p1mup4*(2*(p3mup3*F(282)+p2mup3*F(4026))+p4mup3*F(
     -   4030))+p1mup3*(2*(p3mup4*F(600)-p2mup4*F(4026))+p4mup4*F(403
     -   1))+2*(p2mup4*p4mup3*F(4032)+p2mup3*p4mup4*F(4033)+p3mup4*p4
     -   mup3*F(4034)+p3mup3*p4mup4*F(4035)+2*(p2mup4*p3mup3*F(4036)+
     -   p2mup3*p3mup4*F(4037)))
       Fa18(1782) = p1mup4*(2*(p4mup3*F(362)+p2mup3*F(4027))+p3mup3*F(
     -   4038))+p1mup3*(2*(p4mup4*F(3969)-p2mup4*F(4027))+p3mup4*F(40
     -   39))+2*(p2mup4*p3mup3*F(4032)+p2mup3*p3mup4*F(4033)+p3mup4*p
     -   4mup3*F(4040)+p3mup3*p4mup4*F(4041)+2*(p2mup4*p4mup3*F(4042)
     -   +p2mup3*p4mup4*F(4043)))
       Return
       End
