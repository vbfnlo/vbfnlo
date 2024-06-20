c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineAbe33(p1mup2,p1mup3,p1mup4,p1m
     -   up5,p2mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,
     -   p4mup2,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mu
     -   p2,p6mup3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,
     -   mup3mup5,mup4mup5,Fa33)
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
       Complex*16 Fa33(353:363)
       COMMON/test_FaFunctionsAbe/Fa
       Fa33(353) = -(p1mup3*p1mup4*p3mup5*F(4451))-p1mup3*p3mup4*p3mup
     -   5*F(4453)-p1mup3*p3mup5*p4mup4*F(4454)+p1mup3*p1mup4*p6mup5*
     -   F(4455)+p1mup3*p3mup4*p6mup5*F(4457)+p1mup3*p4mup4*p6mup5*F(
     -   4458)+p1mup4*p3mup3*p6mup5*F(4460)+p3mup3*p3mup4*p6mup5*F(44
     -   65)+p3mup3*p4mup4*p6mup5*F(4466)-p1mup3*p1mup4*p4mup5*F(4467
     -   )-p1mup3*p3mup4*p4mup5*F(4469)-p1mup3*p4mup4*p4mup5*F(4470)+
     -   p1mup4*p4mup3*p4mup5*F(4471)+p2mup5*p3mup3*p3mup4*F(4472)+p2
     -   mup5*p3mup3*p4mup4*F(4473)+p2mup3*(p2mup5*p3mup4*F(4297)+p2m
     -   up5*p4mup4*F(4333)-p1mup4*p3mup5*F(4352)-p3mup4*p3mup5*F(435
     -   4)-p3mup5*p4mup4*F(4355)+p1mup4*p6mup5*F(4459)+p3mup4*p6mup5
     -   *F(4462)+p4mup4*p6mup5*F(4463)+p1mup4*p4mup5*F(4476)+p3mup4*
     -   p4mup5*F(4477)+p4mup4*p4mup5*F(4478))+p1mup4*p3mup3*p4mup5*F
     -   (4480)+p3mup3*p3mup4*p4mup5*F(4481)+p3mup3*p4mup4*p4mup5*F(4
     -   482)-p3mup4*p3mup5*p4mup3*F(4483)-p3mup4*p4mup3*p4mup5*F(448
     -   4)-p1mup5*p4mup3*p4mup4*F(4485)-p1mup5*p3mup3*p4mup4*F(4486)
     -   -p1mup4*p3mup3*p3mup5*F(4487)-p3mup3*p3mup4*p3mup5*F(4488)-p
     -   3mup3*p3mup5*p4mup4*F(4489)-p3mup4*p3mup5*p6mup3*F(4490)+p2m
     -   up4*(p2mup5*p3mup3*F(4297)+p2mup3*p2mup5*F(4331)-p2mup3*p3mu
     -   p5*F(4353)-p3mup3*p3mup5*F(4354)-p1mup3*p3mup5*F(4452)+p1mup
     -   3*p6mup5*F(4456)+p2mup3*p6mup5*F(4461)+p3mup3*p6mup5*F(4464)
     -   -p1mup3*p4mup5*F(4468)-p2mup3*p4mup5*F(4474)-p3mup3*p4mup5*F
     -   (4475)-p4mup3*p4mup5*F(4479)-p1mup5*p3mup3*F(4491))-p1mup5*p
     -   3mup3*p3mup4*F(4492)+p2mup5*p4mup3*p4mup4*F(4493)+p3mup5*p4m
     -   up3*p4mup4*F(4494)-p3mup3*p3mup5*p6mup4*F(4495)+p4mup3*p4mup
     -   4*p4mup5*F(4496)
       Fa33(354) = mup3mup5*(p3mup4*F(4498)+8*(p2mup4*F(4497)+p4mup4*F
     -   (4499)+p1mup4*F(4500)))+16*(-(p1mup5*p2mup4*p4mup3*F(4502))-
     -   p1mup5*p3mup4*p4mup3*F(4504)-p2mup4*p2mup5*p4mup3*F(4506)-p2
     -   mup5*p3mup4*p4mup3*F(4507)+p1mup3*p2mup5*p6mup4*F(4508)-p2mu
     -   p4*p2mup5*p6mup3*F(4512)-p2mup4*p4mup3*p6mup5*F(4513)+p2mup4
     -   *p4mup5*p6mup3*F(4514)-p2mup5*p3mup4*p6mup3*F(4515)+p2mup5*p
     -   4mup3*p6mup4*F(4516)-p2mup5*p4mup4*p6mup3*F(4517)+p4mup3*p4m
     -   up4*p6mup5*F(4518)-p2mup4*p3mup5*p4mup3*F(4519)+p3mup4*p4mup
     -   3*p6mup5*F(4520)+p1mup5*p2mup4*p6mup3*F(4522)-p2mup3*p6mup4*
     -   (p2mup5*F(4320)+p3mup5*F(4326)+p1mup5*F(4510)+p4mup5*F(4511)
     -   +p6mup5*F(4523))-p2mup4*p3mup5*p6mup3*F(4524)-p2mup4*p6mup3*
     -   p6mup5*F(4525)+p1mup3*p3mup5*p6mup4*F(4526)+p1mup5*p3mup4*p6
     -   mup3*F(4528)-p3mup3*p6mup4*p6mup5*F(4529)+p3mup4*p4mup5*p6mu
     -   p3*F(4530)-p3mup4*p6mup3*p6mup5*F(4531)+p3mup5*p4mup3*p6mup4
     -   *F(4532)-p3mup5*p4mup4*p6mup3*F(4533)+p1mup3*p4mup5*p6mup4*F
     -   (4534)+p1mup5*p4mup4*p6mup3*F(4535)-p4mup4*p6mup3*p6mup5*F(4
     -   536)+p1mup3*p1mup5*p6mup4*F(4537)-p1mup3*p6mup4*p6mup5*F(453
     -   8)+p4mup3*p4mup5*p6mup4*F(4541)-p4mup4*p4mup5*p6mup3*F(4542)
     -   +p1mup4*(p2mup5*p4mup3*F(4501)+p3mup5*p4mup3*F(4503)+p1mup5*
     -   p4mup3*F(4505)-p2mup5*p6mup3*F(4509)+p4mup3*p6mup5*F(4521)-p
     -   3mup5*p6mup3*F(4527)+p1mup5*p6mup3*F(4539)-p6mup3*p6mup5*F(4
     -   540)-p4mup5*p6mup3*F(4543))-p1mup5*p4mup3*p6mup4*F(4544)-p1m
     -   up5*p3mup3*p6mup4*F(4545)+p1mup5*p6mup3*p6mup4*F(4546)-p2mup
     -   5*p3mup3*p6mup4*F(4547)+p2mup5*p6mup3*p6mup4*F(4548)+p3mup3*
     -   p4mup5*p6mup4*F(4549)-p4mup3*p6mup4*p6mup5*F(4550)+p3mup5*p6
     -   mup3*p6mup4*F(4551)-p4mup5*p6mup3*p6mup4*F(4552)+p6mup3*p6mu
     -   p4*p6mup5*F(4553))
       Fa33(355) = 4*mup4mup5*F(4554)+8*(p2mup5*p3mup4*F(4560)+p3mup4*
     -   p3mup5*F(4561)+p3mup4*p4mup5*F(4564)-p3mup4*p6mup5*F(4565)-p
     -   4mup4*p6mup5*F(4566)+p1mup5*(-(p3mup4*F(4556))+p2mup4*F(4558
     -   )+p4mup4*F(4562)-p6mup4*F(4568))+p4mup4*p4mup5*F(4569)+p1mup
     -   4*(p3mup5*F(4555)+p2mup5*F(4557)+p6mup5*F(4567)+p1mup5*F(457
     -   0)+p4mup5*F(4571))-p2mup4*(p3mup5*F(4559)+p6mup5*F(4563)-p2m
     -   up5*F(4572)+p4mup5*F(4573))+p2mup5*p6mup4*F(4574)-p3mup5*p6m
     -   up4*F(4575)+p4mup5*p6mup4*F(4576)+p6mup4*p6mup5*F(4577)+p2mu
     -   p5*p4mup4*F(4578)+p3mup5*p4mup4*F(4579))
       Fa33(356) = mup4mup5*(p3mup2*F(4580)+p1mup2*F(4581)+p2mup2*F(45
     -   84)+p4mup2*F(4587)+p6mup2*F(4588))+mup2mup4*(-(p1mup5*F(4582
     -   ))+p2mup5*F(4585)+p6mup5*F(4590)+p3mup5*F(4591)+p4mup5*F(459
     -   2))+mup2mup5*(p1mup4*F(4583)+p2mup4*F(4586)+p6mup4*F(4589)+p
     -   3mup4*F(4593)-p4mup4*F(4594))+2*(-(p2mup4*p3mup2*p6mup5*F(46
     -   02))-p3mup2*p3mup4*p6mup5*F(4603)-p3mup2*p4mup4*p6mup5*F(460
     -   4)+p2mup5*p3mup2*p3mup4*F(4605)-p2mup4*p3mup2*p3mup5*F(4606)
     -   +p2mup2*p3mup4*p3mup5*F(4607)+p1mup2*p3mup4*p4mup5*F(4608)-p
     -   1mup2*p3mup5*p4mup4*F(4609)+p2mup2*p2mup4*p3mup5*F(4611)+p2m
     -   up5*p3mup4*p4mup2*F(4612)+p2mup2*p3mup4*p4mup5*F(4612)-p2mup
     -   4*p3mup5*p4mup2*F(4614)+p2mup5*p3mup2*p6mup4*F(4615)+p1mup2*
     -   p3mup5*p6mup4*F(4616)+p1mup2*p2mup4*p3mup5*F(4619)-p1mup2*p2
     -   mup5*p3mup4*F(4620)-p1mup2*p3mup4*p6mup5*F(4621)-p2mup2*p2mu
     -   p5*p3mup4*F(4623)+p2mup4*p3mup5*p6mup2*F(4624)-p2mup5*p3mup4
     -   *p6mup2*F(4625)-p3mup2*p3mup4*p3mup5*F(4626)-p3mup2*p3mup5*p
     -   4mup4*F(4627)+p3mup2*p3mup4*p4mup5*F(4628)+p3mup2*p3mup5*p6m
     -   up4*F(4629)+p3mup2*p4mup5*p6mup4*F(4631)+p3mup2*p6mup4*p6mup
     -   5*F(4632)-p1mup2*p3mup4*p3mup5*F(4633)-p3mup4*p4mup2*p4mup5*
     -   F(4635)-2*(p2mup4*p2mup5*p3mup2*F(4612)+p2mup5*p3mup2*p4mup4
     -   *F(4613)-p2mup4*p3mup2*p4mup5*F(4635)-p3mup2*p4mup4*p4mup5*F
     -   (4636))-p3mup4*p3mup5*p4mup2*F(4637)-p3mup4*p3mup5*p6mup2*F(
     -   4638)+p1mup4*(p2mup5*p3mup2*F(4595)+p3mup2*(p3mup5*F(4598)+p
     -   1mup5*F(4599)+p4mup5*F(4600)-p6mup5*F(4601))+p3mup5*(p2mup2*
     -   F(4610)+p1mup2*F(4617)-p4mup2*F(4639)-p6mup2*F(4640)))-p3mup
     -   4*p4mup5*p6mup2*F(4642)-p1mup5*(p2mup4*p3mup2*F(4596)+p3mup2
     -   *(p3mup4*F(4597)+p6mup4*F(4630)-p4mup4*F(4634))+p3mup4*(p1mu
     -   p2*F(4618)+p2mup2*F(4622)+p4mup2*F(4641)-p6mup2*F(4643)))-p2
     -   mup2*p3mup5*p4mup4*F(4644)+p2mup2*p3mup5*p6mup4*F(4645)-p2mu
     -   p2*p3mup4*p6mup5*F(4646)-p3mup5*p4mup2*p4mup4*F(4647)+p3mup5
     -   *p4mup4*p6mup2*F(4648)+p3mup5*p4mup2*p6mup4*F(4648)-p3mup4*p
     -   4mup2*p6mup5*F(4649)+p3mup5*p6mup2*p6mup4*F(4650)+p3mup4*p6m
     -   up2*p6mup5*F(4651))
       Fa33(357) = mup4mup5*(p2mup3*F(4652)-p3mup3*F(4653)-p1mup3*F(46
     -   54)+p4mup3*F(4655))
       Fa33(358) = mup4mup5*p6mup3*F(4656)+mup3mup5*p6mup4*F(4657)+mup
     -   3mup4*p6mup5*F(4658)
       Fa33(359) = mup3mup4*(p1mup5*F(4659)+p2mup5*F(4660)+p3mup5*F(46
     -   61)+p4mup5*F(4662))
       Fa33(360) = p1mup4*(p2mup3*p2mup5*F(4665)+p2mup5*p3mup3*F(4666)
     -   +p1mup5*p2mup3*F(4671)+p1mup5*p3mup3*F(4672))+p1mup3*(p1mup4
     -   *p2mup5*F(4663)-p1mup5*(p2mup4*F(4664)+p3mup4*F(4668))+p1mup
     -   4*p1mup5*F(4670)+p2mup5*p3mup4*F(4673)-p2mup4*p2mup5*F(4674)
     -   -p2mup5*p4mup4*F(4675)+p1mup5*p4mup4*F(4676))-p1mup5*p2mup3*
     -   (p2mup4*F(4667)+p3mup4*F(4669)-p4mup4*F(4677))
       Fa33(361) = p1mup5*p3mup3*(p2mup4*F(4678)+p3mup4*F(4679))-p1mup
     -   4*(p2mup3*p3mup5*F(4681)+p3mup3*p3mup5*F(4682)+p2mup3*p4mup5
     -   *F(4684)+p3mup3*p4mup5*F(4685))+p1mup3*(-(p1mup4*p3mup5*F(46
     -   80))-p1mup4*p4mup5*F(4683)+p1mup4*p6mup5*F(4686)+p2mup4*p6mu
     -   p5*F(4687)+p3mup4*p6mup5*F(4688))+p1mup4*p2mup3*p6mup5*F(469
     -   0)+p1mup4*p3mup3*p6mup5*F(4691)+p2mup3*p2mup4*p6mup5*F(4692)
     -   +p2mup3*p3mup4*p6mup5*F(4693)+p2mup3*p4mup4*p6mup5*F(4694)+p
     -   2mup4*p3mup3*p6mup5*F(4695)+p3mup3*p3mup4*p6mup5*F(4696)+p3m
     -   up3*p4mup4*p6mup5*F(4697)+p2mup3*p2mup4*p3mup5*F(4701)+p2mup
     -   4*p3mup3*p3mup5*F(4702)+p2mup5*(p4mup3*p4mup4*F(4396)+p2mup3
     -   *(-(p3mup4*F(4698))+p2mup4*F(4703)+p4mup4*F(4704))+p3mup3*(-
     -   (p3mup4*F(4699))+p2mup4*F(4705)+p4mup4*F(4706)))-p2mup3*p3mu
     -   p4*p3mup5*F(4707)-p2mup3*p3mup5*p4mup4*F(4708)+p3mup3*p3mup4
     -   *p3mup5*F(4711)+p3mup3*p3mup5*p4mup4*F(4712)-p2mup3*p3mup4*p
     -   4mup5*F(4714)-p3mup3*p3mup4*p4mup5*F(4715)-p3mup3*p3mup5*p6m
     -   up4*F(4716)-p1mup5*p3mup3*p4mup4*F(4717)+p1mup3*(p2mup4*p3mu
     -   p5*F(4700)+p3mup4*p3mup5*F(4709)-p3mup4*p4mup5*F(4713)-p2mup
     -   4*p4mup5*F(4718)+p4mup4*(p6mup5*F(4689)+p3mup5*F(4710)-p4mup
     -   5*F(4719)))-p2mup3*p2mup4*p4mup5*F(4720)-p2mup3*p4mup4*p4mup
     -   5*F(4721)-p2mup4*p3mup3*p4mup5*F(4722)-p3mup3*p4mup4*p4mup5*
     -   F(4723)+p3mup4*p3mup5*p4mup3*F(4724)+p3mup4*p3mup5*p6mup3*F(
     -   4725)+p3mup4*p4mup3*p4mup5*F(4726)+p1mup4*p4mup3*p4mup5*F(47
     -   27)+p1mup5*p4mup3*p4mup4*F(4728)+p3mup5*p4mup3*p4mup4*F(4729
     -   )-p2mup4*p4mup3*p4mup5*F(4730)-p4mup3*p4mup4*p4mup5*F(4731)
       Fa33(362) = -(mup3mup5*(p1mup4*F(4732)+p2mup4*F(4733)+p3mup4*F(
     -   4734)+p4mup4*F(4735)))+2*(-(p2mup3*p3mup5*p6mup4*F(4744))+p1
     -   mup4*p3mup5*p4mup3*F(4747)+p1mup3*p3mup5*p6mup4*F(4749)-p1mu
     -   p4*p3mup5*p6mup3*F(4750)+p4mup3*p4mup4*p6mup5*F(4752)+p2mup3
     -   *p4mup5*p6mup4*F(4754)+p2mup3*p6mup4*p6mup5*F(4755)+p3mup3*p
     -   4mup5*p6mup4*F(4759)+p3mup3*p6mup4*p6mup5*F(4760)+p3mup4*p4m
     -   up3*p6mup5*F(4761)-p3mup4*p4mup5*p6mup3*F(4762)+p3mup4*p6mup
     -   3*p6mup5*F(4763)+p1mup3*p4mup5*p6mup4*F(4764)+p1mup4*p4mup3*
     -   p6mup5*F(4766)-p1mup4*p4mup5*p6mup3*F(4767)+p4mup4*p6mup3*p6
     -   mup5*F(4768)+p1mup3*p6mup4*p6mup5*F(4770)+p1mup4*p6mup3*p6mu
     -   p5*F(4772)-p1mup5*(p2mup4*p4mup3*F(4737)-p2mup4*p6mup3*F(474
     -   2)+p3mup4*p4mup3*F(4748)-p3mup4*p6mup3*F(4751)+p2mup3*p6mup4
     -   *F(4753)+p3mup3*p6mup4*F(4758)+p1mup4*p4mup3*F(4765)+p1mup3*
     -   p6mup4*F(4769)+p1mup4*p6mup3*F(4771)-p4mup3*p6mup4*F(4773)-p
     -   4mup4*p6mup3*F(4774)+p6mup3*p6mup4*F(4775))+p2mup5*(p1mup4*p
     -   4mup3*F(4736)+p3mup4*p4mup3*F(4739)-p1mup4*p6mup3*F(4741)-p3
     -   mup4*p6mup3*F(4746)-p2mup4*p4mup3*F(4776)+p2mup4*p6mup3*F(47
     -   77)+p4mup4*p6mup3*F(4779)+p6mup4*(p3mup3*F(4517)+p1mup3*F(47
     -   40)+p2mup3*F(4743)+p4mup3*F(4778)-p6mup3*F(4780)))+p3mup5*p4
     -   mup3*p6mup4*F(4781)+p3mup5*p4mup4*p6mup3*F(4782)-p3mup5*p6mu
     -   p3*p6mup4*F(4783)+p4mup3*p4mup5*p6mup4*F(4784)+p4mup4*p4mup5
     -   *p6mup3*F(4785)+p2mup4*(-(p3mup5*p4mup3*F(4738))+p3mup5*p6mu
     -   p3*F(4745)+p4mup3*p6mup5*F(4756)+p6mup3*p6mup5*F(4757)+p4mup
     -   5*p6mup3*F(4786))+p4mup3*p6mup4*p6mup5*F(4787)-p4mup5*p6mup3
     -   *p6mup4*F(4788)-p6mup3*p6mup4*p6mup5*F(4789))
       Fa33(363) = mup4mup5*(p1mup2*F(4790)+p2mup2*F(4791)-p3mup2*F(47
     -   92)-p4mup2*F(4793))
       Return
       End
