c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineAbe32(p1mup2,p1mup3,p1mup4,p1m
     -   up5,p2mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,
     -   p4mup2,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mu
     -   p2,p6mup3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,
     -   mup3mup5,mup4mup5,Fa32)
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
       Complex*16 Fa32(342:352)
       COMMON/test_FaFunctionsAbe/Fa
       Fa32(342) = -(mup2mup5*p6mup4*F(4269))+mup4mup5*p6mup2*F(4270)+
     -   mup2mup4*p6mup5*F(4271)
       Fa32(343) = mup2mup4*(p1mup5*F(4272)+p4mup5*F(4273)-p3mup5*F(42
     -   74)+p2mup5*F(4275))+2*(-(p2mup2*p3mup4*p3mup5*F(4297))-p2mup
     -   2*p4mup4*p4mup5*F(4299)+p2mup2*p6mup4*p6mup5*F(4302)-p2mup5*
     -   p4mup2*p6mup4*F(4303)+p1mup2*(p2mup4*p2mup5*F(4282)-p2mup4*p
     -   6mup5*F(4283)+p2mup5*p3mup4*F(4284)+p2mup5*p4mup4*F(4285)-p3
     -   mup4*p6mup5*F(4286)-p4mup4*p6mup5*F(4287)+p3mup4*p3mup5*F(42
     -   89)+p2mup4*p3mup5*F(4291)+p3mup5*p4mup4*F(4292)+p4mup4*p4mup
     -   5*F(4293)+p1mup4*(p1mup5*F(4276)+p2mup5*F(4277)-p6mup5*F(427
     -   8)+p3mup5*F(4290)+p4mup5*F(4294))+p2mup4*p4mup5*F(4295)+p3mu
     -   p4*p4mup5*F(4296)-p2mup5*p6mup4*F(4301)-p4mup5*p6mup4*F(4304
     -   )+p1mup5*(p2mup4*F(4279)+p3mup4*F(4280)+p4mup4*F(4281)-p6mup
     -   4*F(4305))+p6mup4*p6mup5*F(4306)-p3mup5*p6mup4*F(4308))+p3mu
     -   p2*p6mup4*p6mup5*F(4309)-p3mup5*p4mup2*p6mup4*F(4310)-p4mup2
     -   *p4mup5*p6mup4*F(4311)+p1mup5*p4mup2*p6mup4*F(4312)-p1mup4*p
     -   2mup2*(p6mup5*F(4288)+p3mup5*F(4298)+p4mup5*F(4300)-p2mup5*F
     -   (4307)-p1mup5*F(4313))+p1mup5*p2mup2*p2mup4*F(4314)+p1mup5*p
     -   2mup2*p6mup4*F(4315)+p1mup5*p3mup2*p6mup4*F(4316)-p1mup5*p6m
     -   up2*p6mup4*F(4317)+p2mup2*p2mup5*p6mup4*F(4318)-p2mup2*p4mup
     -   5*p6mup4*F(4319)+p2mup5*p3mup2*p6mup4*F(4320)+p2mup2*p3mup5*
     -   p6mup4*F(4320)-p3mup2*p4mup5*p6mup4*F(4321)+p4mup2*p6mup4*p6
     -   mup5*F(4322)-p2mup5*p6mup2*p6mup4*F(4323)+p4mup5*p6mup2*p6mu
     -   p4*F(4324)-p6mup2*p6mup4*p6mup5*F(4325)+p3mup2*p3mup5*p6mup4
     -   *F(4326)-p3mup5*p6mup2*p6mup4*F(4327))
       Fa32(344) = p2mup2*(-(p2mup4*(p6mup5*F(4328)+p3mup5*F(4331)+p4m
     -   up5*F(4334)-p2mup5*F(4336)))+p4mup4*(-(p6mup5*F(4330))-p3mup
     -   5*F(4333)+p1mup5*F(4335)+p2mup5*F(4338))+p3mup4*(-(p6mup5*F(
     -   4329))+p4mup5*F(4332)+p2mup5*F(4337)+p1mup5*F(4339)))
       Fa32(345) = p3mup2*(p2mup5*(p3mup4*F(4297)+p1mup4*F(4298)+p2mup
     -   4*F(4331)+p4mup4*F(4333))+p2mup4*p6mup5*F(4341)+p3mup4*p6mup
     -   5*F(4342)+p4mup4*p6mup5*F(4343)-p1mup5*p2mup4*F(4345)-p1mup5
     -   *p3mup4*F(4346)-p1mup5*p4mup4*F(4347)-p2mup4*p4mup5*F(4348)+
     -   p3mup4*p4mup5*F(4350)+p4mup4*p4mup5*F(4351)+p1mup4*(p6mup5*F
     -   (4340)-p1mup5*F(4344)+p4mup5*F(4349)-p3mup5*F(4352))-p2mup4*
     -   p3mup5*F(4353)-p3mup4*p3mup5*F(4354)-p3mup5*p4mup4*F(4355))
       Fa32(346) = mup2mup5*(p2mup4*F(4356)+p3mup4*F(4357)+p4mup4*F(43
     -   58))-2*(p1mup4*p2mup5*p4mup2*F(4359)-p2mup4*p4mup2*p4mup5*F(
     -   4365)+p4mup2*p4mup4*p6mup5*F(4367)-p2mup4*p2mup5*p4mup2*F(43
     -   68)+p2mup4*p4mup2*p6mup5*F(4369)-p2mup5*p3mup4*p4mup2*F(4370
     -   )+p3mup4*p4mup2*p6mup5*F(4371)-p2mup4*p2mup5*p6mup2*F(4375)+
     -   p2mup4*p4mup5*p6mup2*F(4376)-p2mup4*p6mup2*p6mup5*F(4377)-p2
     -   mup5*p3mup4*p6mup2*F(4378)-p2mup5*p4mup4*p6mup2*F(4379)-p3mu
     -   p4*p6mup2*p6mup5*F(4381)-p1mup5*(p2mup4*p4mup2*F(4360)+p3mup
     -   4*p4mup2*F(4362)-p1mup4*p4mup2*F(4364)+p4mup2*p4mup4*F(4366)
     -   -p2mup4*p6mup2*F(4374)-p3mup4*p6mup2*F(4380)-p4mup4*p6mup2*F
     -   (4382))-p4mup4*p6mup2*p6mup5*F(4383)-p3mup4*p3mup5*p4mup2*F(
     -   4384)-p3mup4*p4mup2*p4mup5*F(4385)-p2mup4*p3mup5*p4mup2*F(43
     -   86)-p3mup4*p3mup5*p6mup2*F(4387)-p2mup4*p3mup5*p6mup2*F(4388
     -   )+p3mup4*p4mup5*p6mup2*F(4390)-p3mup5*p4mup4*p6mup2*F(4391)-
     -   p4mup4*p4mup5*p6mup2*F(4392)-p1mup4*(-(p3mup5*p4mup2*F(4361)
     -   )-p4mup2*p4mup5*F(4363)-p4mup2*p6mup5*F(4372)+p2mup5*p6mup2*
     -   F(4373)+p3mup5*p6mup2*F(4389)+p4mup5*p6mup2*F(4393))+p2mup5*
     -   p4mup2*p4mup4*F(4394)+p3mup5*p4mup2*p4mup4*F(4395)+p4mup2*p4
     -   mup4*p4mup5*F(4396))
       Fa32(347) = p1mup4*(mup2mup5*F(4397)+2*(p1mup5*p6mup2*F(4398)-p
     -   6mup2*p6mup5*F(4399)))
       Fa32(348) = mup4mup5*F(4400)+2*(-(p2mup4*p2mup5*F(4405))+p2mup4
     -   *p6mup5*F(4406)-p2mup5*p3mup4*F(4407)-p2mup5*p4mup4*F(4408)+
     -   p3mup4*p6mup5*F(4409)+p4mup4*p6mup5*F(4410)-p2mup4*p4mup5*F(
     -   4411)+p3mup4*p3mup5*F(4412)+p2mup4*p3mup5*F(4414)-p3mup4*p4m
     -   up5*F(4415)+p3mup5*p4mup4*F(4416)+p4mup4*p4mup5*F(4417)-p1mu
     -   p5*(-(p2mup4*F(4402))+p3mup4*F(4403)+p4mup4*F(4404)+p6mup4*F
     -   (4420))+p2mup5*p6mup4*F(4421)-p6mup4*p6mup5*F(4422)-p4mup5*p
     -   6mup4*F(4423)+p1mup4*(p4mup5*F(4401)+p2mup5*F(4413)+p1mup5*F
     -   (4418)-p6mup5*F(4419)+p3mup5*F(4424))-p3mup5*p6mup4*F(4425))
       Fa32(349) = mup4mup5*(p4mup3*F(4426)+p1mup3*F(4427)+p3mup3*F(44
     -   28)+p2mup3*F(4429))
       Fa32(350) = mup3mup5*p6mup4*F(4430)+mup4mup5*p6mup3*F(4431)+mup
     -   3mup4*p6mup5*F(4432)
       Fa32(351) = mup3mup4*(p1mup5*F(4433)-p4mup5*F(4434)+p3mup5*F(44
     -   35)+p2mup5*F(4436))
       Fa32(352) = -(p1mup5*p2mup3*(p2mup4*F(4442)+p3mup4*F(4443)+p4mu
     -   p4*F(4444)))+p1mup3*(p2mup4*p2mup5*F(4438)+p2mup5*p3mup4*F(4
     -   439)+p2mup5*p4mup4*F(4440)+p1mup5*p2mup4*F(4446)+p1mup5*p3mu
     -   p4*F(4447)+p1mup5*p4mup4*F(4448))+p1mup4*(p1mup3*p2mup5*F(44
     -   37)-p2mup3*(p2mup5*F(4298)+p1mup5*F(4441))+p1mup3*p1mup5*F(4
     -   445)-p2mup5*p3mup3*F(4449)+p1mup5*p3mup3*F(4450))
       Return
       End
