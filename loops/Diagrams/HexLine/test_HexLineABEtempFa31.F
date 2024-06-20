c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineAbe31(p1mup2,p1mup3,p1mup4,p1m
     -   up5,p2mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,
     -   p4mup2,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mu
     -   p2,p6mup3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,
     -   mup3mup5,mup4mup5,Fa31)
       IMPLICIT NONE
      Complex*16   p1mup2, p1mup3, p1mup4, p1mup5, p2mup2, p2mup3, 
     -          p2mup4, p2mup5, p3mup2, p3mup3, p3mup4, p3mup5, 
     -          p4mup2, p4mup3, p4mup4, p4mup5, p5mup2, p5mup3, 
     -          p5mup4, p5mup5, p6mup2, p6mup3, p6mup4, p6mup5
       Complex*16   mup2mup3, mup2mup4, mup2mup5, mup3mup4, mup3mup5, 
     -          mup4mup5
        common/test_FhlFunctionsAbe/F
       COMMON/Invariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12,s23,s34,s45
     -   ,s56,s16,s123,s234,s345
       Real*8   p1sq 
       Real*8   p2sq 
       Real*8   p3sq 
       Real*8   p4sq 
       Real*8   p5sq 
       Real*8   p6sq 
       Real*8   s12, s16 
       Real*8   s23 
       Real*8   s34 
       Real*8   s45 
       Real*8   s56 
       Real*8   s123,s234,s345
      Complex*16 Fa(445),F(5930)
      Real*8 P(11674) 
        COMMON/test_PFunctionsAbe/P
       Complex*16 Fa31(331:341)
       COMMON/test_FaFunctionsAbe/Fa
       Fa31(331) = p1mup3*p3mup5*p4mup4*F(3996)+p3mup3*p3mup5*p4mup4*F
     -   (4010)+p3mup3*p3mup4*p4mup5*F(4012)+p1mup3*p4mup4*p4mup5*F(4
     -   015)+p3mup3*p4mup4*p4mup5*F(4016)+p3mup4*p4mup3*p4mup5*F(401
     -   9)+p1mup5*(p2mup4*p3mup3*F(3985)+p3mup3*p3mup4*F(3986)+16*(p
     -   3mup3*p4mup4*F(3987)+p4mup3*p4mup4*F(4020)))+p1mup3*p1mup4*p
     -   3mup5*F(4021)+p1mup3*p1mup4*p4mup5*F(4022)+p1mup4*p4mup3*p4m
     -   up5*F(4029)+p2mup4*p3mup3*p3mup5*F(4035)+p2mup3*(p2mup4*(p2m
     -   up5*F(3988)+16*(p6mup5*F(3999)+p3mup5*F(4005)+p4mup5*F(4008)
     -   ))+16*(p2mup5*p3mup4*F(3989)-p2mup5*p4mup4*F(3990)-p1mup4*p4
     -   mup5*F(3997)+p3mup4*p6mup5*F(4000)-p4mup4*p6mup5*F(4001)+p3m
     -   up4*p3mup5*F(4006)+p3mup5*p4mup4*F(4007)+p3mup4*p4mup5*F(401
     -   1)+p1mup4*p6mup5*F(4025)+p1mup4*p3mup5*F(4030))+p4mup4*p4mup
     -   5*F(4037))+16*(p2mup4*p2mup5*p3mup3*F(3991)-p2mup5*p3mup3*p3
     -   mup4*F(3992)+p2mup5*p3mup3*p4mup4*F(3993)-p1mup3*p2mup4*p3mu
     -   p5*F(3994)-p1mup3*p3mup4*p3mup5*F(3995)-p1mup4*p3mup3*p4mup5
     -   *F(3998)-p2mup4*p3mup3*p6mup5*F(4002)+p3mup3*p3mup4*p6mup5*F
     -   (4003)-p3mup3*p4mup4*p6mup5*F(4004)+p2mup4*p3mup3*p4mup5*F(4
     -   009)-p1mup3*p2mup4*p4mup5*F(4013)-p1mup3*p3mup4*p4mup5*F(401
     -   4)+p3mup4*p3mup5*p4mup3*F(4017)+p2mup5*p4mup3*p4mup4*F(4018)
     -   -p1mup3*p2mup4*p6mup5*F(4023)-p1mup3*p4mup4*p6mup5*F(4024)+p
     -   1mup4*p3mup3*p6mup5*F(4026)-p2mup4*p4mup3*p4mup5*F(4027)-p3m
     -   up5*p4mup3*p4mup4*F(4028)+p3mup4*p3mup5*p6mup3*F(4031)-p1mup
     -   3*p1mup4*p6mup5*F(4032)-p1mup3*p3mup4*p6mup5*F(4033)-p1mup4*
     -   p3mup3*p3mup5*F(4034)-p3mup3*p3mup4*p3mup5*F(4036)+p3mup3*p3
     -   mup5*p6mup4*F(4038)-p4mup3*p4mup4*p4mup5*F(4039))
       Fa31(332) = -8*mup3mup5*(p2mup4*F(4040)-p3mup4*F(4041)-p4mup4*F
     -   (4042)+p1mup4*F(4043))+p1mup5*p2mup4*p4mup3*F(4047)+p2mup4*p
     -   2mup5*p4mup3*F(4051)+p2mup3*p2mup5*p6mup4*F(4058)-16*(p1mup5
     -   *p3mup4*p4mup3*F(4045)+p2mup4*p3mup5*p4mup3*F(4048)-p2mup5*p
     -   3mup4*p4mup3*F(4049)-p4mup3*p4mup4*p6mup5*F(4050)+p2mup4*p4m
     -   up3*p6mup5*F(4052)-p3mup4*p4mup3*p6mup5*F(4053)-p1mup3*p2mup
     -   5*p6mup4*F(4054)+p1mup5*p2mup3*p6mup4*F(4056)-p1mup5*p2mup4*
     -   p6mup3*F(4057)-p2mup3*p6mup4*p6mup5*F(4059)+p2mup4*p2mup5*p6
     -   mup3*F(4060)-p2mup4*p6mup3*p6mup5*F(4061)-p2mup5*p3mup3*p6mu
     -   p4*F(4062)+p2mup5*p3mup4*p6mup3*F(4063)+p2mup5*p4mup4*p6mup3
     -   *F(4064)-p1mup3*p3mup5*p6mup4*F(4065)+p1mup5*p3mup3*p6mup4*F
     -   (4066)-p1mup5*p3mup4*p6mup3*F(4067)-p3mup3*p6mup4*p6mup5*F(4
     -   068)-p3mup4*p6mup3*p6mup5*F(4069)-p1mup3*p4mup5*p6mup4*F(407
     -   0)+p1mup5*p4mup3*p6mup4*F(4072)-p1mup5*p4mup4*p6mup3*F(4073)
     -   -p4mup4*p6mup3*p6mup5*F(4074)-p4mup3*p4mup5*p6mup4*F(4075)+p
     -   4mup4*p4mup5*p6mup3*F(4076)-p2mup3*p4mup5*p6mup4*F(4077)-p2m
     -   up4*p4mup5*p6mup3*F(4078)-p2mup5*p4mup3*p6mup4*F(4079)-p3mup
     -   3*p4mup5*p6mup4*F(4080)-p3mup4*p4mup5*p6mup3*F(4081)+p3mup5*
     -   p4mup3*p6mup4*F(4082)-p4mup3*p6mup4*p6mup5*F(4084)-p2mup4*p3
     -   mup5*p6mup3*F(4085)-p2mup5*p6mup3*p6mup4*F(4086)-p3mup5*p4mu
     -   p4*p6mup3*F(4087)+p1mup5*p6mup3*p6mup4*F(4090)+p6mup3*p6mup4
     -   *p6mup5*F(4091)+p2mup3*p3mup5*p6mup4*F(4092)-p3mup5*p6mup3*p
     -   6mup4*F(4094)-p1mup4*(p3mup5*p4mup3*F(4044)-p2mup5*p4mup3*F(
     -   4046)-p2mup5*p6mup3*F(4055)+p4mup3*p6mup5*F(4071)-p1mup5*p4m
     -   up3*F(4083)-p1mup5*p6mup3*F(4088)-p6mup3*p6mup5*F(4089)+p3mu
     -   p5*p6mup3*F(4093)+p4mup5*p6mup3*F(4095))+p4mup5*p6mup3*p6mup
     -   4*F(4096)-p1mup3*p1mup5*p6mup4*F(4097)-p1mup3*p6mup4*p6mup5*
     -   F(4098))
       Fa31(333) = mup4mup5*(p4mup2*F(4099)-p1mup2*F(4100)+p2mup2*F(41
     -   01)+p3mup2*F(4102))
       Fa31(334) = mup4mup5*p6mup2*F(4103)-mup2mup5*p6mup4*F(4104)+mup
     -   2mup4*p6mup5*F(4105)
       Fa31(335) = mup2mup4*(p1mup5*F(4106)-p2mup5*F(4107)-p3mup5*F(41
     -   08)-p4mup5*F(4109))+2*(p2mup2*(p3mup4*p3mup5*F(4111)-p1mup5*
     -   p6mup4*F(4124)+p1mup4*(p2mup5*F(4117)+p3mup5*F(4118)-p4mup5*
     -   F(4122)+p6mup5*F(4136)+p1mup5*F(4150))+p6mup4*(p3mup5*F(4125
     -   )-p4mup5*F(4138)-p2mup5*F(4143)-p6mup5*F(4151))+p1mup5*p2mup
     -   4*F(4155)+p4mup4*p4mup5*F(4157))+p1mup2*(p3mup5*p4mup4*F(411
     -   3)-p2mup5*p4mup4*F(4116)-p4mup4*p4mup5*F(4119)+p2mup5*p6mup4
     -   *F(4123)+p3mup5*p6mup4*F(4127)+p4mup5*p6mup4*F(4130)+p1mup4*
     -   p3mup5*F(4132)-p4mup4*p6mup5*F(4135)+p1mup5*p4mup4*F(4141)+p
     -   1mup5*p6mup4*F(4146)+p6mup4*p6mup5*F(4147)+p2mup4*(p3mup5*F(
     -   4112)-p2mup5*F(4114)-p4mup5*F(4120)-p6mup5*F(4133)+p1mup5*F(
     -   4148))-p3mup4*(p3mup5*F(4110)+p2mup5*F(4115)+p4mup5*F(4121)+
     -   p6mup5*F(4134)-p1mup5*F(4149))+p1mup4*p1mup5*F(4153)+p1mup4*
     -   p2mup5*F(4154)+p1mup4*p6mup5*F(4156)+p1mup4*p4mup5*F(4158))+
     -   p6mup4*(p3mup2*p6mup5*F(4129)+p4mup2*p4mup5*F(4137)+p3mup2*p
     -   4mup5*F(4140)-p3mup2*p3mup5*F(4142)-p3mup5*p4mup2*F(4144)-p3
     -   mup5*p6mup2*F(4145)-p4mup2*p6mup5*F(4152)-p2mup5*(p3mup2*F(4
     -   126)+p4mup2*F(4139)-p6mup2*F(4159))+p4mup5*p6mup2*F(4160)+p1
     -   mup5*(p3mup2*F(4128)-p4mup2*F(4131)+p6mup2*F(4161))+p6mup2*p
     -   6mup5*F(4162)))
       Fa31(336) = p2mup2*(p3mup5*p4mup4*F(4164)+p1mup5*p3mup4*F(4165)
     -   +p1mup5*p4mup4*F(4166)+p4mup4*p6mup5*F(4167)+p3mup4*p6mup5*F
     -   (4168)+p2mup5*p3mup4*F(4170)+p2mup5*p4mup4*F(4171)+p2mup4*(p
     -   3mup5*F(4163)+p2mup5*F(4169)+p6mup5*F(4172)-p4mup5*F(4173))-
     -   p3mup4*p4mup5*F(4174))
       Fa31(337) = p3mup2*(p2mup4*p3mup5*F(4186)+16*(-(p2mup5*p3mup4*F
     -   (4179))+p2mup5*p4mup4*F(4180)+p1mup4*p4mup5*F(4181)+p3mup4*p
     -   6mup5*F(4183)+p4mup4*p6mup5*F(4184)+p2mup4*(p2mup5*F(4178)+p
     -   6mup5*F(4182)+p4mup5*F(4185))+p3mup4*p3mup5*F(4187)+p3mup5*p
     -   4mup4*F(4188)-p3mup4*p4mup5*F(4189)+p4mup4*p4mup5*F(4190)-p1
     -   mup5*(p2mup4*F(4175)+p3mup4*F(4176)-p4mup4*F(4177)-p1mup4*F(
     -   4191))+p1mup4*p2mup5*F(4192)+p1mup4*p6mup5*F(4193)+p1mup4*p3
     -   mup5*F(4194)))
       Fa31(338) = mup2mup5*(p2mup4*F(4195)+p3mup4*F(4196)+p4mup4*F(41
     -   97))+2*(p2mup4*p3mup5*p4mup2*F(4202)-p2mup5*p3mup4*p4mup2*F(
     -   4203)+p3mup4*p4mup2*p4mup5*F(4204)-p3mup5*p4mup2*p4mup4*F(42
     -   05)-p3mup4*p3mup5*p6mup2*F(4207)-p2mup4*p3mup5*p6mup2*F(4210
     -   )-p3mup5*p4mup4*p6mup2*F(4212)+p4mup4*p4mup5*p6mup2*F(4216)-
     -   p2mup4*p4mup2*p6mup5*F(4217)-p2mup4*p4mup5*p6mup2*F(4218)+p2
     -   mup5*p4mup4*p6mup2*F(4219)-p3mup4*p4mup2*p6mup5*F(4220)-p3mu
     -   p4*p4mup5*p6mup2*F(4221)+p2mup4*p2mup5*p6mup2*F(4223)+p2mup5
     -   *p3mup4*p6mup2*F(4224)+p1mup4*(-(p3mup5*p4mup2*F(4198))+p2mu
     -   p5*p4mup2*F(4200)+p4mup2*p4mup5*F(4206)+p2mup5*p6mup2*F(4208
     -   )-p4mup2*p6mup5*F(4213)+p4mup5*p6mup2*F(4214)+p1mup5*p4mup2*
     -   F(4222)+p3mup5*p6mup2*F(4225))+p2mup4*p6mup2*p6mup5*F(4226)+
     -   p3mup4*p6mup2*p6mup5*F(4227)+p4mup4*p6mup2*p6mup5*F(4228)-p3
     -   mup4*p3mup5*p4mup2*F(4229)-p2mup5*p4mup2*p4mup4*F(4230)+p1mu
     -   p5*(p3mup4*p4mup2*F(4199)+p2mup4*p4mup2*F(4201)+p2mup4*p6mup
     -   2*F(4209)+p3mup4*p6mup2*F(4211)+p4mup4*p6mup2*F(4215)-p4mup2
     -   *p4mup4*F(4231))+p2mup4*p2mup5*p4mup2*F(4232)+p2mup4*p4mup2*
     -   p4mup5*F(4233)-p4mup2*p4mup4*p6mup5*F(4234)-p4mup2*p4mup4*p4
     -   mup5*F(4235))
       Fa31(339) = p1mup4*(mup2mup5*F(4236)+2*p6mup2*(p1mup5*F(4237)+p
     -   6mup5*F(4238)))
       Fa31(340) = 4*mup4mup5*F(4239)+p3mup5*p4mup4*F(4242)+p1mup5*p4m
     -   up4*F(4246)-8*(p1mup4*p4mup5*F(4243)-p4mup4*p6mup5*F(4249)-p
     -   1mup4*p3mup5*F(4250)-p1mup4*p6mup5*F(4251)+p1mup5*p6mup4*F(4
     -   252)-p3mup5*p6mup4*F(4253)-p2mup4*(p3mup5*F(4241)+p1mup5*F(4
     -   244)+p6mup5*F(4247)-p2mup5*F(4254)-p4mup5*F(4255))-p2mup5*p4
     -   mup4*F(4257)-p4mup4*p4mup5*F(4258)-p3mup4*(p3mup5*F(4240)+p1
     -   mup5*F(4245)-p6mup5*F(4248)-p2mup5*F(4256)-p4mup5*F(4259))-p
     -   4mup5*p6mup4*F(4260)-p1mup4*p1mup5*F(4261)-p1mup4*p2mup5*F(4
     -   262)+p2mup5*p6mup4*F(4263)+p6mup4*p6mup5*F(4264))
       Fa31(341) = mup4mup5*(p4mup2*F(4265)+p1mup2*F(4266)-p3mup2*F(42
     -   67)-p2mup2*F(4268))
       Return
       End
