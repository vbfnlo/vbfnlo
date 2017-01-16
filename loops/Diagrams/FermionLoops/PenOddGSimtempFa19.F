c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenOddGSim19(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4
     -   ,p3mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5m
     -   up3,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mu
     -   p3,mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa19)
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
       Complex*16 Fa19(1783:1881)
       COMMON/FaPenOddGSim/Fa
       Fa19(1783) = p3mup3*F(3697)+p4mup3*F(3698)
       Fa19(1784) = p3mup3*F(3779)+p4mup3*F(3780)
       Fa19(1785) = 2*p3mup3*F(3705)+p4mup3*F(3706)
       Fa19(1786) = p3mup3*F(3709)+2*p4mup3*F(3710)
       Fa19(1787) = 2*(p1mup4*p2mup3*F(258)+p1mup3*p2mup4*F(264)+p2mup
     -   4*p3mup3*F(282)+p2mup3*p3mup4*F(600))+p2mup4*p4mup3*F(4038)+
     -   p2mup3*p4mup4*F(4039)
       Fa19(1788) = p1mup4*p2mup3*F(4028)+p1mup3*p2mup4*F(4029)+2*(p2m
     -   up4*p3mup3*F(4034)+p2mup3*p3mup4*F(4035)+p2mup4*p4mup3*F(404
     -   0)+p2mup3*p4mup4*F(4041))
       Fa19(1789) = 2*(p1mup4*p2mup3*F(346)+p2mup4*p4mup3*F(362)+p2mup
     -   3*p4mup4*F(3969)+p1mup3*p2mup4*F(4024))+p2mup4*p3mup3*F(4030
     -   )+p2mup3*p3mup4*F(4031)
       Fa19(1790) = p3mup3*F(3757)+p4mup3*F(3758)
       Fa19(1791) = p3mup3*F(3701)+p4mup3*F(3702)
       Fa19(1792) = p3mup3*(4*p3mup4*F(382)+p1mup4*F(3997)+2*(p2mup4*F
     -   (4000)+p4mup4*F(4004)))+p1mup3*(p3mup4*F(3996)+2*(p1mup4*F(4
     -   51)+p2mup4*F(3998)+p4mup4*F(4044)))+2*(p2mup3*(p1mup4*F(3995
     -   )+p3mup4*F(3999)+2*(p2mup4*F(456)+p4mup4*F(4003)))+p4mup3*(p
     -   1mup4*F(4002)+p3mup4*F(4005)+2*(p2mup4+p4mup4)*F(4467)))
       Fa19(1793) = p3mup3*F(4045)+p4mup3*F(4046)
       Fa19(1794) = p3mup3*F(4047)-p4mup3*F(4048)
       Fa19(1795) = p4mup3*F(4049)+2*p3mup3*F(4050)
       Fa19(1796) = p3mup3*F(4051)+2*p4mup3*F(4052)
       Fa19(1797) = 2*p1mup3*p2mup4*F(264)-p2mup3*(2*(p1mup4*F(4014)+p
     -   3mup4*F(4015))+p4mup4*F(4016))+p2mup4*(2*(p3mup3*F(282)-p2mu
     -   p3*F(449))+p4mup3*F(4017))
       Fa19(1798) = p1mup4*p2mup3*F(4018)-p1mup3*p2mup4*F(4018)-2*((p2
     -   mup4*p3mup3-p2mup3*p3mup4)*F(4019)+(p2mup4*p4mup3-p2mup3*p4m
     -   up4)*F(4020))
       Fa19(1799) = -2*p1mup4*p2mup3*F(346)+p2mup3*(2*(p2mup4*F(449)-p
     -   4mup4*F(3969))+p3mup4*F(4022))+p2mup4*(p3mup3*F(4023)+2*(p1m
     -   up3*F(4021)+p4mup3*F(4025)))
       Fa19(1800) = p2mup3*p2mup4*F(4057)+p2mup3*p3mup4*F(4058)+p2mup4
     -   *p3mup3*F(4059)+p3mup3*p3mup4*F(4060)+p1mup3*(8*p2mup4*F(390
     -   0)+p1mup4*F(4053)+p3mup4*F(4055)+p4mup4*F(4061))+p2mup3*p4mu
     -   p4*F(4062)+p3mup3*p4mup4*F(4063)+p3mup4*p4mup3*F(4064)+p1mup
     -   4*(p2mup3*F(4054)+p3mup3*F(4056)+p4mup3*F(4065))+p2mup4*p4mu
     -   p3*F(4066)+p4mup3*p4mup4*F(4067)
       Fa19(1801) = p2mup3*(-8*(p2mup4*F(3913)+p3mup4*F(3915))-4*p4mup
     -   4*F(3917)+p1mup4*F(4069))-8*p3mup3*p3mup4*F(4070)+p3mup3*p4m
     -   up4*F(4071)+p1mup3*p3mup4*F(4072)+p2mup4*(-8*p3mup3*F(3916)+
     -   p1mup3*F(4068)+p4mup3*F(4073))-4*(p1mup3*p1mup4*F(3919)+p4mu
     -   p3*p4mup4*F(3924)-p3mup4*p4mup3*F(4074))+p1mup4*p3mup3*F(407
     -   5)+p1mup3*p4mup4*F(4076)+p1mup4*p4mup3*F(4077)
       Fa19(1802) = -8*(p4mup3*p4mup4*F(3934)+p2mup4*p4mup3*F(3935))-4
     -   *(p3mup3*(p3mup4*F(3927)+p2mup4*F(3930)+p4mup4*F(3932))+p1mu
     -   p3*p1mup4*F(3936))+p3mup4*p4mup3*F(4079)+p2mup3*(-8*(p2mup4*
     -   F(3929)+p4mup4*F(3931))+p3mup4*F(4078)+p1mup4*F(4080))+p1mup
     -   3*p3mup4*F(4081)+p1mup4*p3mup3*F(4082)+p1mup3*p2mup4*F(4083)
     -   +p1mup3*p4mup4*F(4084)+p1mup4*p4mup3*F(4085)
       Fa19(1803) = p3mup3*F(4086)+p4mup3*F(4087)
       Fa19(1804) = p3mup3*F(4088)+p4mup3*F(4089)
       Fa19(1805) = p4mup3*F(4090)-2*p3mup3*F(4091)
       Fa19(1806) = p3mup3*F(4092)+2*p4mup3*F(4093)
       Fa19(1807) = -2*p2mup4*p4mup3*F(3957)+4*(p1mup3*p2mup4*F(3951)+
     -   p2mup4*p3mup3*F(4096))+p2mup3*(4*(p2mup4*F(3953)+p1mup4*F(40
     -   94))+p3mup4*F(4095)+p4mup4*F(4097))
       Fa19(1808) = 4*(p2mup3*p2mup4*F(3958)+p2mup4*p3mup3*F(3959)-p2m
     -   up3*p3mup4*F(3960)-p2mup3*p4mup4*F(3961)+p2mup4*p4mup3*F(396
     -   2))+2*(-(p1mup4*p2mup3)+p1mup3*p2mup4)*F(3963)
       Fa19(1809) = 4*(p2mup3*p2mup4*F(3964)+p1mup3*p2mup4*F(3965)+p2m
     -   up4*p4mup3*F(3968)+m1sq*p2mup3*(p4mup4*F(362)+p1mup4*F(4024)
     -   ))+p2mup3*p3mup4*F(4098)+p2mup4*p3mup3*F(4099)
       Fa19(1810) = -(p1mup2*(2*(p1mup3*F(4100)+p2mup3*F(4101))+p3mup3
     -   *F(4102)))+p1mup3*p3mup2*F(4103)+p2mup3*p3mup2*F(4104)+p3mup
     -   2*p3mup3*F(4105)+p4mup2*p4mup3*F(4106)
       Fa19(1811) = p1mup3*p3mup2*F(4111)-4*(p1mup2*(p1mup3*F(4107)+p2
     -   mup3*F(4108)-p3mup3*F(4110))-p4mup2*p4mup3*F(4112))+8*(p2mup
     -   3*p3mup2*F(4109)+p3mup2*p3mup3*F(4113))
       Fa19(1812) = -4*(p2mup3*p3mup2*F(4114)+p3mup2*p3mup3*F(4115)+2*
     -   p4mup2*p4mup3*F(4116)+p1mup2*p1mup3*F(4117)+p1mup2*p2mup3*F(
     -   4118))-2*p1mup2*p3mup3*F(4119)+p1mup3*p3mup2*F(4120)
       Fa19(1813) = p1mup2*p3mup4*F(4121)+4*(m1sq*(p1mup2*(2*(p1mup4*F
     -   (633)+p2mup4*F(635))+p4mup4*F(1367))-p4mup2*(p2mup4*F(640)+p
     -   1mup4*F(1368)-p4mup4*F(1370)))-p3mup2*(p1mup4*F(4122)+p2mup4
     -   *F(4123)))+p3mup2*p3mup4*F(4124)+p3mup2*p4mup4*F(4125)+p3mup
     -   4*p4mup2*F(4126)
       Fa19(1814) = p3mup2*(-8*(p3mup4*F(4129)+p2mup4*F(4130))+p4mup4*
     -   F(4131)+p1mup4*F(4136))+p1mup2*(4*p1mup4*F(4127)+p2mup4*F(41
     -   28)+p3mup4*F(4132)+p4mup4*F(4137))+p4mup2*(p2mup4*F(4133)-4*
     -   (p3mup4*F(4134)+p4mup4*F(4135))+2*p1mup4*F(4138))
       Fa19(1815) = p3mup4*(p4mup2*F(4142)+p1mup2*F(4143))-2*p3mup2*(2
     -   *(p3mup4*F(4139)+p2mup4*F(4140)+p4mup4*F(4141))+p1mup4*F(414
     -   4))+4*m1sq*(p1mup2*(p1mup4*F(677)+p2mup4*F(679)-p4mup4*F(140
     -   0))+p4mup2*(p1mup4*F(1401)+2*(p4mup4*F(1402)+p2mup4*F(4145))
     -   ))
       Fa19(1816) = p1mup2*(p2mup4*F(4147)+8*(p1mup4*F(4146)+p3mup4*F(
     -   4148))+p4mup4*F(4152))+p3mup2*(8*(p3mup4*F(4149)+p1mup4*F(41
     -   50))+p2mup4*F(4151)+p4mup4*F(4154))+p4mup2*(p2mup4*F(4153)+p
     -   4mup4*F(4155)-8*(p1mup4*F(4156)+p3mup4*F(4157)))
       Fa19(1817) = 2*p2mup4*(2*(p3mup2*F(4158)+p4mup2*F(4159))+p1mup2
     -   *F(4165))+4*(p1mup2*p1mup4*F(4164)+p1mup4*p3mup2*F(4167)+2*(
     -   p3mup2*p3mup4*F(4162)+p3mup2*p4mup4*F(4163)+p1mup2*p3mup4*F(
     -   4166)+p1mup2*p4mup4*F(4168))+p4mup2*(2*(p4mup4*F(4160)+p3mup
     -   4*F(4161))+p1mup4*F(4169)))
       Fa19(1818) = p1mup2*p3mup4*F(4170)+2*p2mup4*p3mup2*F(4171)+p3mu
     -   p2*p3mup4*F(4172)+p3mup4*p4mup2*F(4173)+4*m1sq*(p1mup2*p2mup
     -   4*F(939)+p4mup2*(p2mup4*F(946)+2*(p4mup4*F(365)+p1mup4*F(417
     -   5)))+2*(p1mup2*p1mup4*F(938)+p1mup4*p3mup2*F(4174)+p1mup2*p4
     -   mup4*F(4175)+p3mup2*p4mup4*F(4176)))
       Fa19(1819) = p2mup3*p3mup2*F(4179)+8*(p1mup3*p3mup2*F(4180)+p3m
     -   up2*p3mup3*F(4181))-4*(p1mup2*(p2mup3*F(4178)+2*(p1mup3*F(41
     -   77)-p3mup3*F(4180)))-p4mup2*p4mup3*F(4182))
       Fa19(1820) = p2mup3*p3mup2*F(4184)+2*(p1mup2*(2*p1mup3*F(4183)+
     -   p2mup3*F(4185)+p3mup3*F(4186))+p3mup2*p3mup3*F(4188))+4*(p1m
     -   up3*p3mup2*F(4187)+p4mup2*p4mup3*F(4189))
       Fa19(1821) = p1mup2*(2*p1mup3*F(4164)+p2mup3*F(4165)+4*p3mup3*F
     -   (4166))+2*p3mup2*(p1mup3*F(4167)+p2mup3*F(4190)+2*p3mup3*F(4
     -   191))+4*p4mup2*p4mup3*F(4192)
       Fa19(1822) = p1mup2*(2*(p1mup3*F(633)+p2mup3*F(635))+p3mup3*F(1
     -   261))-p1mup3*p3mup2*F(1263)-p2mup3*p3mup2*F(1265)+p3mup2*p3m
     -   up3*F(1274)-p4mup2*p4mup3*F(4193)
       Fa19(1823) = p1mup2*(p1mup3*F(1289)+p2mup3*F(1290)-p3mup3*F(130
     -   0))+p1mup3*p3mup2*F(1301)+2*(p3mup2*p3mup3*F(1306)+p4mup2*p4
     -   mup3*F(1308)+p2mup3*p3mup2*F(4194))
       Fa19(1824) = p1mup2*(2*(p1mup3*F(4195)+p2mup3*F(4196))+p3mup3*F
     -   (4197))+p1mup3*p3mup2*F(4198)+2*p2mup3*p3mup2*F(4199)+4*(p3m
     -   up2*p3mup3*F(4200)+p4mup2*p4mup3*F(4201))
       Fa19(1825) = p3mup2*p4mup4*F(4205)+p3mup4*p4mup2*F(4206)+p1mup2
     -   *(4*p1mup4*F(633)-2*(2*p2mup4*F(4207)+p3mup4*F(4208)+p4mup4*
     -   F(4209)))-2*(p2mup4*p4mup2*F(640)+p1mup4*p4mup2*F(1368)-p3mu
     -   p2*p3mup4*F(4202)-p1mup4*p3mup2*F(4203)-p2mup4*p3mup2*F(4204
     -   )+p4mup2*p4mup4*F(4210))
       Fa19(1826) = p1mup2*(2*(p1mup4*F(4211)+p2mup4*F(4212)+p3mup4*F(
     -   4213))+p4mup4*F(4214))-p1mup4*p4mup2*F(4215)-2*(p2mup4*p4mup
     -   2*F(4216)-p1mup4*p3mup2*F(4218)-p3mup2*p4mup4*F(4220)+p3mup4
     -   *p4mup2*F(4221))+4*(p3mup2*p3mup4*F(4217)+p2mup4*p3mup2*F(42
     -   19)+p4mup2*p4mup4*F(4222))
       Fa19(1827) = p1mup4*p3mup2*F(4224)+p1mup2*(p3mup4*F(4223)+2*(p4
     -   mup4*F(1400)+p1mup4*F(4225)))+2*(p2mup4*(2*p4mup2*F(689)+p1m
     -   up2*F(4226)+p3mup2*F(4227))+p3mup2*p4mup4*F(4229)+p3mup4*p4m
     -   up2*F(4230)+p1mup4*p4mup2*F(4231)+2*(p3mup2*p3mup4*F(4228)+p
     -   4mup2*p4mup4*F(4232)))
       Fa19(1828) = p2mup4*p4mup2*F(4237)+2*(p1mup2*(p2mup4*F(4234)+2*
     -   (p1mup4*F(4233)+p3mup4*F(4235))+p4mup4*F(4236))+p2mup4*p3mup
     -   2*F(4240)+p3mup2*p4mup4*F(4241)+p4mup2*p4mup4*F(4242))+4*(p3
     -   mup2*p3mup4*F(4238)+p1mup4*p3mup2*F(4239)+p1mup4*p4mup2*F(42
     -   43)+p3mup4*p4mup2*F(4244))
       Fa19(1829) = p1mup2*(2*p1mup4*F(4245)+p2mup4*F(4246)+4*(p3mup4*
     -   F(4247)+p4mup4*F(4250)))+2*(p1mup4*p3mup2*F(4248)+p2mup4*p3m
     -   up2*F(4249)+p1mup4*p4mup2*F(4251)+p2mup4*p4mup2*F(4252)+2*(p
     -   3mup4*p4mup2*F(4253)+p3mup2*p4mup4*F(4253)+p3mup2*p3mup4*F(4
     -   254)+p4mup2*p4mup4*F(4255)))
       Fa19(1830) = -(p2mup4*p3mup2*F(4257))+2*(p1mup2*(2*p4mup4*F(359
     -   )-p2mup4*F(939)+p3mup4*F(4256))+p3mup2*p3mup4*F(4258)+p3mup4
     -   *p4mup2*F(4259)+p2mup4*p4mup2*F(4260))+4*(p1mup4*(p4mup2*F(3
     -   59)-p1mup2*F(938)+p3mup2*F(944))+p3mup2*p4mup4*F(949)+p4mup2
     -   *p4mup4*F(4261))
       Fa19(1831) = p1mup2*p2mup3*F(1499)+p2mup3*p3mup2*F(1505)-p4mup2
     -   *p4mup3*F(4262)+2*(p1mup2*p1mup3*F(882)+p3mup2*p3mup3*F(1507
     -   )+(p1mup3*p3mup2+p1mup2*p3mup3)*F(4263))
       Fa19(1832) = p2mup3*p3mup2*F(4265)+2*(p1mup2*(2*p1mup3*F(4264)+
     -   p2mup3*F(4266)+p3mup3*F(4267))+p3mup2*p3mup3*F(4269))+4*(p1m
     -   up3*p3mup2*F(4268)+p4mup2*p4mup3*F(4270))
       Fa19(1833) = p1mup2*(2*p1mup3*F(4245)+p2mup3*F(4246)+4*p3mup3*F
     -   (4247))+2*p3mup2*(p1mup3*F(4248)+p2mup3*F(4249)+2*p3mup3*F(4
     -   254))+4*p4mup2*p4mup3*F(4255)
       Fa19(1834) = p1mup2*(2*(p1mup3*F(633)+p2mup3*F(635))+p3mup3*F(1
     -   261))-p1mup3*p3mup2*F(1263)-p2mup3*p3mup2*F(1265)+p3mup2*p3m
     -   up3*F(1274)+p4mup2*p4mup3*F(1370)
       Fa19(1835) = 4*p4mup2*p4mup3*F(1402)+p1mup2*(2*(p1mup3*F(677)+p
     -   2mup3*F(679))+p3mup3*F(4271))+p3mup2*(p1mup3*F(4272)+2*p2mup
     -   3*F(4273)+4*p3mup3*F(4274))
       Fa19(1836) = p1mup2*(-4*(p1mup4*F(633)-p2mup4*F(4207))+2*(p3mup
     -   4*F(4208)+p4mup4*F(4209)))+p3mup2*p4mup4*F(4275)+p3mup4*p4mu
     -   p2*F(4276)+2*(p1mup4*p3mup2*F(1263)+p2mup4*(p4mup2*F(640)+p3
     -   mup2*F(1265))+p1mup4*p4mup2*F(1368)+p4mup2*p4mup4*F(4210)+p3
     -   mup2*p3mup4*F(4277))
       Fa19(1837) = 4*p4mup2*p4mup4*F(1661)+p1mup2*p4mup4*F(4280)-p1mu
     -   p4*p4mup2*F(4281)+2*(p1mup2*p3mup4*F(1300)+p1mup2*p1mup4*F(4
     -   278)+p2mup4*(2*p3mup2*F(1299)+p1mup2*F(4279)+p4mup2*F(4282))
     -   +p3mup4*p4mup2*F(4284)+p3mup2*(p4mup4*F(4283)+p1mup4*F(4285)
     -   +2*p3mup4*F(4286)))
       Fa19(1838) = p1mup2*(2*(p4mup4*F(1400)+p1mup4*F(4225))+p3mup4*F
     -   (4287))+p1mup4*p3mup2*F(4288)+2*(p1mup4*p4mup2*F(4231)+p2mup
     -   4*(2*p4mup2*F(689)+p1mup2*F(4226)+p3mup2*F(4289))+p3mup2*p4m
     -   up4*F(4290)+p3mup4*p4mup2*F(4291))+4*(p4mup2*p4mup4*F(4232)+
     -   p3mup2*p3mup4*F(4292))
       Fa19(1839) = p2mup4*p4mup2*F(4294)+2*(p4mup2*p4mup4*F(598)+p1mu
     -   p2*(2*p3mup4*F(279)-p2mup4*F(1499)+p4mup4*F(4293))+p3mup2*p4
     -   mup4*F(4295)+p2mup4*p3mup2*F(4296))+4*(p1mup4*(p3mup2*F(279)
     -   -p1mup2*F(882)+p4mup2*F(1503))+p3mup4*p4mup2*F(1508)+p3mup2*
     -   p3mup4*F(4297))
       Fa19(1840) = p1mup2*(2*p1mup4*F(4298)+p2mup4*F(4299)+4*(p3mup4*
     -   F(4300)+p4mup4*F(4302)))+2*(p1mup4*p3mup2*F(4301)+p1mup4*p4m
     -   up2*F(4303)+p2mup4*p3mup2*F(4304)+p2mup4*p4mup2*F(4305)+2*(p
     -   3mup2*p3mup4*F(4306)+p3mup4*p4mup2*F(4307)+p3mup2*p4mup4*F(4
     -   307)+p4mup2*p4mup4*F(4308)))
       Fa19(1841) = 4*(p1mup4*(p4mup2*F(359)-p1mup2*F(938)+p3mup2*F(94
     -   4))+p3mup2*p4mup4*F(949)+p4mup2*p4mup4*F(4261))-p2mup4*p3mup
     -   2*F(4310)+2*(p2mup4*p4mup2*F(4260)+p1mup2*(2*p4mup4*F(359)-p
     -   2mup4*F(939)+p3mup4*F(4309))+p3mup2*p3mup4*F(4311)+p3mup4*p4
     -   mup2*F(4312))
       Fa19(1842) = p1mup2*p2mup3*F(1499)+p2mup3*p3mup2*F(1505)+2*(p1m
     -   up2*p1mup3*F(882)+p3mup2*p3mup3*F(1507)+(p1mup3*p3mup2+p1mup
     -   2*p3mup3)*F(4263))+p4mup2*p4mup3*F(4313)
       Fa19(1843) = 4*p4mup2*p4mup3*F(365)+2*p1mup2*(2*p1mup3*F(938)+p
     -   2mup3*F(939)+p3mup3*F(4314))+p3mup2*(4*p1mup3*F(4174)+p2mup3
     -   *F(4310)+2*p3mup3*F(4315))
       Fa19(1844) = 4*p4mup2*p4mup3*F(318)+p1mup2*(2*p1mup3*F(889)+p2m
     -   up3*F(4316)+4*p3mup3*F(4318))+2*p3mup2*(p2mup3*F(4317)+p1mup
     -   3*F(4319)+2*p3mup3*F(4320))
       Fa19(1845) = p1mup3*p3mup2*F(1263)+p2mup3*p3mup2*F(1265)+p1mup2
     -   *(-2*(p1mup3*F(633)-p2mup3*F(4207))+p3mup3*F(4208))+p4mup2*p
     -   4mup3*F(4210)+p3mup2*p3mup3*F(4277)
       Fa19(1846) = p1mup2*p3mup3*F(1300)+p1mup2*p1mup3*F(4278)+p1mup2
     -   *p2mup3*F(4279)+p1mup3*p3mup2*F(4285)+2*(p2mup3*p3mup2*F(129
     -   9)+p4mup2*p4mup3*F(1661)+p3mup2*p3mup3*F(4286))
       Fa19(1847) = 4*p4mup2*p4mup3*F(4232)+p1mup2*(2*(p1mup3*F(4225)+
     -   p2mup3*F(4226))+p3mup3*F(4287))+p3mup2*(p1mup3*F(4288)+2*p2m
     -   up3*F(4289)+4*p3mup3*F(4292))
       Fa19(1848) = -2*(p1mup4*p3mup2*F(1263)+p2mup4*(p4mup2*F(640)+p3
     -   mup2*F(1265))-p3mup2*p3mup4*F(1274)-p1mup2*(2*(p1mup4*F(633)
     -   +p2mup4*F(635))+p3mup4*F(1261)+p4mup4*F(1367))+p1mup4*p4mup2
     -   *F(1368)-p4mup2*p4mup4*F(1370))+p3mup2*p4mup4*F(4321)+p3mup4
     -   *p4mup2*F(4322)
       Fa19(1849) = 4*(p3mup2*p3mup4*F(1306)+p4mup2*p4mup4*F(1308)+p2m
     -   up4*p3mup2*F(4194))+p1mup2*(2*(p1mup4*F(1289)+p2mup4*F(1290)
     -   -p3mup4*F(1300))-p4mup4*F(4280))+p1mup4*p4mup2*F(4281)+2*(p1
     -   mup4*p3mup2*F(1301)+p2mup4*p4mup2*F(4323)+p3mup2*p4mup4*F(43
     -   24)+p3mup4*p4mup2*F(4325))
       Fa19(1850) = p1mup2*(2*(p1mup4*F(677)+p2mup4*F(679)-p4mup4*F(14
     -   00))+p3mup4*F(4271))+p1mup4*p3mup2*F(4272)+4*(p4mup2*p4mup4*
     -   F(1402)+p2mup4*p4mup2*F(4145)+p3mup2*p3mup4*F(4274))+2*(p1mu
     -   p4*p4mup2*F(1401)+p2mup4*p3mup2*F(4273)+p3mup2*p4mup4*F(4326
     -   )+p3mup4*p4mup2*F(4327))
       Fa19(1851) = p2mup4*p4mup2*F(4328)+2*(p2mup4*p3mup2*F(1505)+p1m
     -   up2*(p2mup4*F(1499)+2*(p1mup4*F(882)+p3mup4*F(4263))-p4mup4*
     -   F(4293))+p4mup2*p4mup4*F(4313)+p3mup2*p4mup4*F(4330))+4*(p3m
     -   up2*p3mup4*F(1507)+p1mup4*p3mup2*F(4263)+p1mup4*p4mup2*F(432
     -   9)+p3mup4*p4mup2*F(4331))
       Fa19(1852) = p1mup2*p2mup4*F(4316)+4*(p4mup2*p4mup4*F(318)+p3mu
     -   p4*(p4mup2*F(317)+p1mup2*F(4318)+p3mup2*F(4320))+p1mup2*p4mu
     -   p4*F(4332))+2*(p1mup2*p1mup4*F(889)+p3mup2*(2*p4mup4*F(317)+
     -   p2mup4*F(4317)+p1mup4*F(4319))+p1mup4*p4mup2*F(4333)+p2mup4*
     -   p4mup2*F(4334))
       Fa19(1853) = p3mup2*(4*(p1mup4*F(4174)+p4mup4*F(4176))+p2mup4*F
     -   (4310)+2*p3mup4*F(4315))+2*(p1mup2*(p2mup4*F(939)+2*(p1mup4*
     -   F(938)+p4mup4*F(4175))+p3mup4*F(4314))+p4mup2*(p2mup4*F(946)
     -   +2*(p4mup4*F(365)+p1mup4*F(4175))+p3mup4*F(4335)))
       Fa19(1854) = p4mup2*p4mup3*F(598)-p1mup2*p2mup3*F(1499)+p2mup3*
     -   p3mup2*F(4296)+2*(p1mup3*p3mup2*F(279)+p1mup2*p3mup3*F(279)-
     -   p1mup2*p1mup3*F(882)+p3mup2*p3mup3*F(4297))
       Fa19(1855) = 4*(p1mup3*p3mup2*F(944)+p4mup2*p4mup3*F(4261))-p2m
     -   up3*p3mup2*F(4310)-2*(p1mup2*(2*p1mup3*F(938)+p2mup3*F(939)-
     -   p3mup3*F(4309))-p3mup2*p3mup3*F(4311))
       Fa19(1856) = p1mup2*(2*p1mup3*F(4298)+p2mup3*F(4299)+4*p3mup3*F
     -   (4300))+2*p3mup2*(p1mup3*F(4301)+p2mup3*F(4304)+2*p3mup3*F(4
     -   306))+4*p4mup2*p4mup3*F(4308)
       Fa19(1857) = p1mup3*p3mup2*F(1263)+p2mup3*p3mup2*F(1265)+p4mup2
     -   *p4mup3*F(4193)+p1mup2*(-2*(p1mup3*F(633)-p2mup3*F(4207))+p3
     -   mup3*F(4208))+p3mup2*p3mup3*F(4277)
       Fa19(1858) = p1mup2*(2*(p1mup3*F(4195)+p2mup3*F(4196))+p3mup3*F
     -   (4197))+4*p4mup2*p4mup3*F(4201)+p3mup2*(4*p3mup3*F(4200)-p1m
     -   up3*F(4336)-2*p2mup3*F(4337))
       Fa19(1859) = -2*(p2mup4*p4mup2*F(640)-p1mup2*(2*(p1mup4*F(633)+
     -   p2mup4*F(635))+p3mup4*F(1261)+p4mup4*F(1367))+p1mup4*p4mup2*
     -   F(1368)-p4mup2*p4mup4*F(1370)-p3mup2*p3mup4*F(4202)-p1mup4*p
     -   3mup2*F(4203)-p2mup4*p3mup2*F(4204))+p3mup2*p4mup4*F(4205)+p
     -   3mup4*p4mup2*F(4206)
       Fa19(1860) = p1mup2*(2*(p1mup4*F(677)+p2mup4*F(679)-p4mup4*F(14
     -   00))-p3mup4*F(4223))-p1mup4*p3mup2*F(4224)+4*(p4mup2*p4mup4*
     -   F(1402)+p2mup4*p4mup2*F(4145)-p3mup2*p3mup4*F(4228))+2*(p1mu
     -   p4*p4mup2*F(1401)-p3mup2*p4mup4*F(4229)-p3mup4*p4mup2*F(4230
     -   )+p2mup4*p3mup2*F(4338))
       Fa19(1861) = p3mup2*(4*(p1mup4*F(4174)+p4mup4*F(4176))+p2mup4*F
     -   (4257)-2*p3mup4*F(4258))+2*(p1mup2*(p2mup4*F(939)+2*(p1mup4*
     -   F(938)+p4mup4*F(4175))-p3mup4*F(4256))+p4mup2*(p2mup4*F(946)
     -   +2*(p4mup4*F(365)+p1mup4*F(4175))-p3mup4*F(4259)))
       Fa19(1862) = -(p1mup2*p2mup3*F(1499))+p4mup2*p4mup3*F(4262)+p2m
     -   up3*p3mup2*F(4296)+2*(p1mup3*p3mup2*F(279)+p1mup2*p3mup3*F(2
     -   79)-p1mup2*p1mup3*F(882)+p3mup2*p3mup3*F(4297))
       Fa19(1863) = p1mup2*(2*(p1mup3*F(4100)+p2mup3*F(4101))+p3mup3*F
     -   (4102))-p1mup3*p3mup2*F(4103)-p2mup3*p3mup2*F(4104)-p4mup2*p
     -   4mup3*F(4106)+p3mup2*p3mup3*F(4339)
       Fa19(1864) = p1mup2*(p1mup3*F(4107)+p2mup3*F(4108)-p3mup3*F(411
     -   0))+p1mup3*p3mup2*F(4340)+p4mup2*p4mup3*F(4341)-2*(p2mup3*p3
     -   mup2*F(4109)-p3mup2*p3mup3*F(4342))
       Fa19(1865) = 8*p4mup2*p4mup3*F(4116)+4*(p3mup2*p3mup3*F(4115)+p
     -   1mup2*p1mup3*F(4117)+p2mup3*(p3mup2*F(4114)+p1mup2*F(4118)))
     -   +p1mup2*p3mup3*F(4343)+p1mup3*p3mup2*F(4344)
       Fa19(1866) = 4*m1sq*(p1mup2*(-2*(p1mup4*F(633)-p2mup4*F(4207))+
     -   p4mup4*F(4209))+p4mup2*(p2mup4*F(640)+p1mup4*F(1368)+p4mup4*
     -   F(4210)))+p1mup2*p3mup4*F(4345)+p3mup2*(4*(p1mup4*F(4122)+p2
     -   mup4*F(4123))+p3mup4*F(4346)+p4mup4*F(4347))+p3mup4*p4mup2*F
     -   (4348)
       Fa19(1867) = 4*p1mup2*p1mup4*F(4349)+p3mup2*p3mup4*F(4351)+p3mu
     -   p2*p4mup4*F(4352)+p1mup2*p3mup4*F(4353)+p2mup4*(8*p3mup2*F(4
     -   130)+p1mup2*F(4350)+p4mup2*F(4354))+p4mup2*(4*p4mup4*F(4135)
     -   -2*p1mup4*F(4138)+p3mup4*F(4355))+p1mup4*p3mup2*F(4356)+p1mu
     -   p2*p4mup4*F(4357)
       Fa19(1868) = 2*p3mup2*(2*(p3mup4*F(4139)+p2mup4*F(4140)+p4mup4*
     -   F(4141))+p1mup4*F(4144))+4*m1sq*(p1mup2*p4mup4*F(1400)+p1mup
     -   2*p1mup4*F(4225)+p1mup2*p2mup4*F(4226)+p1mup4*p4mup2*F(4231)
     -   +2*(p2mup4*p4mup2*F(689)+p4mup2*p4mup4*F(4232)))+p3mup4*(p4m
     -   up2*F(4358)+p1mup2*F(4359))
       Fa19(1869) = p3mup2*p3mup4*F(4363)+p1mup4*(p1mup2*F(4360)+8*(p4
     -   mup2*F(4156)+p3mup2*F(4364)))+p2mup4*p3mup2*F(4365)+p1mup2*(
     -   p2mup4*F(4361)+p3mup4*F(4362)+p4mup4*F(4366))+p2mup4*p4mup2*
     -   F(4367)+p3mup2*p4mup4*F(4368)+p4mup2*p4mup4*F(4369)+p3mup4*p
     -   4mup2*F(4370)
       Fa19(1870) = p2mup4*(2*(p3mup2*F(4158)+p4mup2*F(4159))+p1mup2*F
     -   (4165))+2*(p1mup2*p1mup4*F(4164)+p1mup4*p3mup2*F(4167)+2*(p3
     -   mup2*p3mup4*F(4162)+p3mup2*p4mup4*F(4163)+p1mup2*p3mup4*F(41
     -   66)+p1mup2*p4mup4*F(4168))+p4mup2*(2*(p4mup4*F(4160)+p3mup4*
     -   F(4161))+p1mup4*F(4169)))
       Fa19(1871) = -2*p2mup4*p3mup2*F(4171)+4*m1sq*(-(p1mup2*p2mup4*F
     -   (939))+p2mup4*p4mup2*F(4260)+2*(p1mup2*p4mup4*F(359)+p1mup4*
     -   (p4mup2*F(359)-p1mup2*F(938)+p3mup2*F(944))+p3mup2*p4mup4*F(
     -   949)+p4mup2*p4mup4*F(4261)))+p3mup4*(p1mup2*F(4371)+p3mup2*F
     -   (4372)+p4mup2*F(4373))
       Fa19(1872) = p1mup2*(p2mup3*F(4178)+2*(p1mup3*F(4177)-p3mup3*F(
     -   4180)))-p4mup2*p4mup3*F(4182)+p2mup3*p3mup2*F(4374)-2*(p1mup
     -   3*p3mup2*F(4180)-p3mup2*p3mup3*F(4375))
       Fa19(1873) = 2*p2mup3*p3mup2*F(4184)+4*(p1mup2*(2*p1mup3*F(4183
     -   )+p2mup3*F(4185)+p3mup3*F(4186))+p3mup2*p3mup3*F(4188))+8*(p
     -   1mup3*p3mup2*F(4187)+p4mup2*p4mup3*F(4189))
       Fa19(1874) = 2*p1mup2*(2*p1mup3*F(4164)+p2mup3*F(4165)+4*p3mup3
     -   *F(4166))+4*p3mup2*(p1mup3*F(4167)+p2mup3*F(4190)+2*p3mup3*F
     -   (4191))+8*p4mup2*p4mup3*F(4192)
       Fa19(1875) = 4*(p1mup3*p4mup2*F(4376)+p2mup3*p4mup2*F(4377))-2*
     -   p3mup2*p4mup3*F(4378)+p3mup3*p4mup2*F(4379)+p1mup2*p4mup3*F(
     -   4380)
       Fa19(1876) = 2*(p4mup3*(2*p3mup2*F(4382)+p1mup2*F(4384))-p4mup2
     -   *(2*(p2mup3*F(4381)+p3mup3*F(4383))+p1mup3*F(4385)))
       Fa19(1877) = 2*p2mup3*p4mup2*F(4386)-p3mup2*p4mup3*F(4387)+p3mu
     -   p3*p4mup2*F(4388)+p1mup2*p4mup3*F(4389)+p1mup3*p4mup2*F(4390
     -   )
       Fa19(1878) = 8*m4sq*p4mup2*(p1mup3*F(1503)+p3mup3*F(1508))+2*p2
     -   mup3*p4mup2*F(4392)+4*(p1mup2*p4mup3*F(4391)+p3mup2*p4mup3*F
     -   (4393))
       Fa19(1879) = p2mup3*p4mup2*F(4395)+p3mup3*p4mup2*F(4396)+2*(p1m
     -   up3*p4mup2*F(4394)+p1mup2*p4mup3*F(4394)+p3mup2*p4mup3*F(439
     -   7))
       Fa19(1880) = p1mup3*p4mup2*F(4169)+p2mup3*p4mup2*F(4398)+2*(p1m
     -   up2*p4mup3*F(4168)+(p3mup3*p4mup2+p3mup2*p4mup3)*F(4399))
       Fa19(1881) = mup2mup4*(2*p4mup3*F(3957)-4*(p1mup3*F(3951)-p3mup
     -   3*F(3955)+p2mup3*F(4401)))+mup2mup3*(p3mup4*F(3954)+p4mup4*F
     -   (3956)+4*(p1mup4*F(3952)+p2mup4*F(4402)))+2*mup3mup4*(2*(p1m
     -   up2*F(4400)+p3mup2*F(4403))-p4mup2*F(4404))
       Return
       End
