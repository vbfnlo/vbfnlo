c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineNoAbe1g34(p1mup2,p1mup3,p1mup4
     -   ,p1mup6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3m
     -   up6,p4mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,
     -   p6mup2,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3m
     -   up4,mup3mup6,mup4mup6,Fa34)
       IMPLICIT NONE
      Complex*16   p1mup2, p1mup3, p1mup4, p1mup6, p2mup2, p2mup3, 
     -          p2mup4, p2mup6, p3mup2, p3mup3, p3mup4, p3mup6, 
     -          p4mup2, p4mup3, p4mup4, p4mup6, p5mup2, p5mup3, 
     -          p5mup4, p5mup6, p6mup2, p6mup3, p6mup4, p6mup6
       Complex*16   mup2mup3, mup2mup4, mup2mup6, mup3mup4, mup3mup6, 
     -          mup4mup6
        common/test_FhlFunctionsNoAbe1g/F
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
      Complex*16 Fa(587),F(6395)
      Real*8 P(7129) 
        COMMON/test_PFunctionsNoAbe1g/P
       Complex*16 Fa34(463:476)
       COMMON/test_FaFunctionsNoAbe1g/Fa
       Fa34(463) = 16*p3mup6*(p3mup4*(p4mup3*p5mup2+p4mup2*p5mup3)*F(4
     -   483)-p3mup3*p5mup4*(p4mup2*F(4482)+p5mup2*F(4484))+p3mup4*p5
     -   mup2*p5mup3*F(4485))+p3mup2*p5mup4*(-16*p3mup3*(p4mup6*F(448
     -   0)+p5mup6*F(4481))+p3mup6*(p4mup3*F(4486)+p5mup3*F(4487)))
       Fa34(464) = p3mup3*p4mup2*F(4490)+p3mup3*p5mup2*F(4492)+p2mup3*
     -   p4mup2*F(4494)+p2mup3*p5mup2*F(4495)+p1mup2*p2mup3*F(4496)+p
     -   1mup2*p3mup3*F(4497)+p3mup2*(p4mup3*F(4489)+p5mup3*F(4491)+p
     -   1mup3*F(4498)+p2mup3*F(4499)+p3mup3*F(4500))+p4mup2*p5mup3*F
     -   (4502)+p2mup2*(p2mup3*F(4488)+p3mup3*F(4493)+4*(p4mup3*F(450
     -   3)+p5mup3*F(4504))+p1mup3*F(4505))+p1mup3*p4mup2*F(4506)+p5m
     -   up2*p5mup3*F(4507)+p1mup3*p5mup2*F(4509)+p1mup2*p1mup3*F(451
     -   1)+p1mup2*p4mup3*F(4512)+4*(p4mup2*p4mup3*F(4501)+p4mup3*p5m
     -   up2*F(4508)+mup2mup3*F(4510)+p1mup2*p5mup3*F(4513))
       Fa34(465) = p5mup4*(p3mup2*p4mup3*F(4515)+p2mup2*p4mup3*F(4525)
     -   +p2mup2*p5mup3*F(4526)+8*(p4mup2*p4mup3*F(4517)+p3mup2*p5mup
     -   3*F(4519)+p4mup2*p5mup3*F(4521)+p4mup3*p5mup2*F(4524)-p5mup2
     -   *p5mup3*F(4527))+p2mup3*(8*(p4mup2*F(4514)-p5mup2*F(4518))+p
     -   2mup2*F(4523)+p3mup2*F(4528))+p3mup3*(8*(p4mup2*F(4516)-p5mu
     -   p2*F(4520))+p2mup2*F(4522)+p3mup2*F(4529))+mup2mup3*F(4530))
       Fa34(466) = p2mup2*(p1mup3*p5mup4*F(4536)+8*p4mup4*(p2mup3*F(45
     -   33)+p3mup3*F(4534)+p1mup3*F(4542)))+p1mup3*p5mup4*(p3mup2*F(
     -   4538)-8*(p4mup2*F(4540)+p5mup2*F(4544)))+p1mup2*(p2mup3*p5mu
     -   p4*F(4535)+p3mup3*p5mup4*F(4537)+p4mup3*p5mup4*F(4539)+8*(p2
     -   mup3*p4mup4*F(4531)+p3mup3*p4mup4*F(4532)+p1mup3*p4mup4*F(45
     -   41))+p1mup3*p5mup4*F(4543)+p5mup3*p5mup4*F(4545))
       Fa34(467) = p2mup2*(p3mup3*p3mup4*F(4547)+p2mup3*p3mup4*F(4549)
     -   +p4mup3*p4mup4*F(4551)+p3mup4*p4mup3*F(4553)+p3mup4*p5mup3*F
     -   (4554)+p4mup4*p5mup3*F(4555)+p1mup4*p5mup3*F(4556)+p2mup4*p5
     -   mup3*F(4557)-mup3mup4*F(4558)+p1mup3*p3mup4*F(4560))+p1mup2*
     -   (p3mup3*p3mup4*F(4546)+p2mup3*p3mup4*F(4548)+p4mup3*p4mup4*F
     -   (4550)+p3mup4*p4mup3*F(4552)+p1mup3*p3mup4*F(4559)-p3mup4*p5
     -   mup3*F(4561)-p4mup4*p5mup3*F(4562)-p1mup4*p5mup3*F(4563)-p2m
     -   up4*p5mup3*F(4564)+mup3mup4*F(4565))
       Fa34(468) = p3mup2*p3mup4*p5mup3*F(4567)+p3mup2*p4mup4*p5mup3*F
     -   (4569)+p2mup3*p3mup2*p4mup4*F(4577)+p2mup3*p3mup2*p3mup4*F(4
     -   582)+p2mup3*p3mup4*p4mup2*F(4583)+p3mup2*p3mup3*p3mup4*F(458
     -   4)+p3mup2*p3mup3*p4mup4*F(4586)+p3mup2*p3mup4*p4mup3*F(4587)
     -   +p3mup3*p3mup4*p4mup2*F(4588)+p2mup3*p4mup2*p4mup4*F(4589)+p
     -   3mup3*p4mup2*p4mup4*F(4590)+p4mup2*p4mup4*p5mup3*F(4602)+p3m
     -   up4*p4mup2*p5mup3*F(4604)+mup3mup4*p3mup2*F(4605)-mup3mup4*p
     -   4mup2*F(4606)+p3mup2*p4mup3*p4mup4*F(4616)+p3mup4*p4mup2*p4m
     -   up3*F(4617)+p4mup2*p4mup3*p4mup4*F(4618)+p2mup4*(p3mup3*p4mu
     -   p2*F(4579)+p3mup2*(p5mup3*F(4568)+p4mup3*F(4578)+p2mup3*F(45
     -   80)+p3mup3*F(4585))-2*(p5mup2*(p2mup3*F(4571)+p3mup3*F(4573)
     -   )+p5mup2*(p4mup3*F(4592)-p1mup3*F(4595)+p5mup3*F(4599))+mup2
     -   mup3*F(4608))+p4mup2*(p2mup3*F(4581)+p5mup3*F(4603)+p4mup3*F
     -   (4615))+p1mup2*p1mup3*F(4621))+p1mup2*p1mup4*p2mup3*F(4622)+
     -   p2mup2*(p2mup3*p2mup4*F(4566)+p2mup4*p3mup3*F(4570)+p2mup4*p
     -   4mup3*F(4591)+p1mup3*p1mup4*F(4623))+p1mup2*p1mup4*p3mup3*F(
     -   4624)+p1mup3*p1mup4*p3mup2*F(4628)+p1mup2*p1mup4*p4mup3*F(46
     -   29)+p1mup3*p1mup4*p4mup2*F(4630)-2*(p2mup3*p3mup4*p5mup2*F(4
     -   572)+p2mup3*p4mup4*p5mup2*F(4574)+p3mup3*p3mup4*p5mup2*F(457
     -   5)+p3mup3*p4mup4*p5mup2*F(4576)+p3mup4*p4mup3*p5mup2*F(4593)
     -   +p4mup3*p4mup4*p5mup2*F(4594)-p1mup3*p3mup4*p5mup2*F(4596)-p
     -   1mup3*p4mup4*p5mup2*F(4597)+p1mup4*p4mup3*p5mup2*F(4598)+p3m
     -   up4*p5mup2*p5mup3*F(4600)+p4mup4*p5mup2*p5mup3*F(4601)-mup2m
     -   up4*p2mup3*F(4607)-mup2mup4*p3mup3*F(4609)+mup2mup3*p3mup4*F
     -   (4610)-mup2mup4*p4mup3*F(4611)+mup2mup3*p4mup4*F(4612)+mup2m
     -   up4*p5mup3*F(4613)-mup2mup4*p1mup3*F(4614)-p1mup4*p5mup2*p5m
     -   up3*F(4619)-mup3mup4*p5mup2*F(4620)+p1mup3*p1mup4*p5mup2*F(4
     -   625)+p1mup4*p2mup3*p5mup2*F(4626)+p1mup4*p3mup3*p5mup2*F(462
     -   7)+mup2mup3*p1mup4*F(4631))
       Fa34(469) = p1mup2*(p2mup3*p2mup4*F(4632)+p2mup4*p3mup3*F(4633)
     -   +p2mup4*p4mup3*F(4640)+p1mup3*p1mup4*F(4647))+p1mup4*(p2mup3
     -   *p4mup2*F(4635)+p3mup2*p4mup3*F(4638)+p3mup3*p4mup2*F(4639)+
     -   p4mup2*p4mup3*F(4641)+p4mup2*p5mup3*F(4642)+p3mup2*p5mup3*F(
     -   4643)+p2mup2*p2mup3*F(4645)+p2mup2*p3mup3*F(4646)+p2mup3*p3m
     -   up2*F(4648)+p3mup2*p3mup3*F(4649)+p2mup2*p4mup3*F(4651))+p1m
     -   up3*(p2mup4*p3mup2*F(4634)+p3mup2*p3mup4*F(4636)-p3mup2*p4mu
     -   p4*F(4637)+p2mup2*p2mup4*F(4644)+p2mup4*p4mup2*F(4650)+p3mup
     -   4*p4mup2*F(4652)+p4mup2*p4mup4*F(4653))
       Fa34(470) = mup4mup6*(p4mup3*F(4654)+p1mup3*F(4656)+4*(p2mup3*F
     -   (4655)+p3mup3*F(4657)+p5mup3*F(4658)))+mup3mup6*p5mup4*F(465
     -   9)
       Fa34(471) = mup3mup6*(p3mup4*F(4660)-p4mup4*F(4661))
       Fa34(472) = p5mup6*(mup3mup4*F(4662)+2*(p3mup3*p3mup4*F(4668)+p
     -   2mup4*p3mup3*F(4669)+p2mup3*p3mup4*F(4669)+p2mup3*p4mup4*F(4
     -   670)+p3mup3*p4mup4*F(4671)+p2mup3*p2mup4*F(4672)-p2mup4*p4mu
     -   p3*F(4675)-p3mup4*p4mup3*F(4676)+p2mup4*p5mup3*F(4677)+p3mup
     -   4*p5mup3*F(4678)-p4mup3*p5mup4*F(4679)+p4mup4*p5mup3*F(4680)
     -   -p1mup3*(p1mup4*F(4663)+p2mup4*F(4664)+p3mup4*F(4665)+p4mup4
     -   *F(4666)+p5mup4*F(4681))+p1mup4*(p4mup3*F(4667)+p2mup3*F(467
     -   3)+p3mup3*F(4674)+p5mup3*F(4682))+p4mup3*p4mup4*F(4683)+p2mu
     -   p3*p5mup4*F(4684)+p3mup3*p5mup4*F(4685)+p5mup3*p5mup4*F(4686
     -   )))
       Fa34(473) = mup3mup6*(p2mup4*F(4687)+p1mup4*F(4688))+p5mup4*(-(
     -   p1mup3*p4mup6*F(4690))-p4mup3*(p4mup6*F(4689)+p2mup6*F(4691)
     -   +p3mup6*F(4692)+p1mup6*F(4693))+p2mup3*p4mup6*F(4694)+p3mup3
     -   *p4mup6*F(4695)+p4mup6*p5mup3*F(4696)+p2mup6*p5mup3*F(4697)+
     -   p3mup6*p5mup3*F(4698))
       Fa34(474) = mup3mup4*(p3mup6*F(4699)+p4mup6*F(4700)-p1mup6*F(47
     -   01)+p2mup6*F(4702))+2*(-(p1mup3*p1mup6*(p1mup4*F(4703)+p2mup
     -   4*F(4704)+p3mup4*F(4705)+p4mup4*F(4706)))+p2mup3*p2mup4*p3mu
     -   p6*F(4711)+p2mup4*p3mup3*p3mup6*F(4712)+p2mup3*p3mup4*p3mup6
     -   *F(4712)+p2mup3*p3mup6*p4mup4*F(4713)+p2mup4*p2mup6*p3mup3*F
     -   (4714)+p2mup3*p2mup6*p3mup4*F(4714)+p2mup6*p3mup3*p3mup4*F(4
     -   715)+p2mup6*p3mup3*p4mup4*F(4716)+p2mup3*p2mup4*p4mup6*F(471
     -   8)+p2mup3*p3mup4*p4mup6*F(4719)+p2mup3*p4mup4*p4mup6*F(4720)
     -   +p2mup4*p3mup3*p4mup6*F(4722)+p3mup3*p3mup4*p4mup6*F(4723)+p
     -   3mup3*p4mup4*p4mup6*F(4724)+p2mup3*p2mup4*p2mup6*F(4726)+p2m
     -   up3*p2mup6*p4mup4*F(4727)+p3mup3*p3mup4*p3mup6*F(4729)-p2mup
     -   4*p2mup6*p4mup3*F(4730)-p2mup6*p3mup4*p4mup3*F(4731)-p2mup4*
     -   p3mup6*p4mup3*F(4732)-p3mup4*p3mup6*p4mup3*F(4733)+p2mup4*p2
     -   mup6*p5mup3*F(4735)+p2mup6*p3mup4*p5mup3*F(4736)+p2mup6*p4mu
     -   p4*p5mup3*F(4737)+p2mup4*p3mup6*p5mup3*F(4739)+p3mup4*p3mup6
     -   *p5mup3*F(4740)+p3mup6*p4mup4*p5mup3*F(4741)+p1mup4*(p2mup6*
     -   (p4mup3*F(4707)+p3mup3*F(4710)+p2mup3*F(4725)+p5mup3*F(4734)
     -   )+p3mup6*(p4mup3*F(4708)+p2mup3*F(4709)+p3mup3*F(4728)+p5mup
     -   3*F(4738))+p4mup6*(p2mup3*F(4717)+p3mup3*F(4721)+p5mup3*F(47
     -   42)))+p2mup4*p4mup6*p5mup3*F(4743)+p3mup4*p4mup6*p5mup3*F(47
     -   44)+p4mup4*p4mup6*p5mup3*F(4745)+p2mup6*p4mup3*p4mup4*F(4746
     -   )+p3mup6*p4mup3*p4mup4*F(4747)+p1mup6*p5mup3*p5mup4*F(4748))
       Fa34(475) = -(p1mup3*(p2mup4*p2mup6*F(4751)+p2mup6*p3mup4*F(475
     -   2)+p2mup6*p4mup4*F(4753)+p2mup4*p3mup6*F(4755)+p3mup4*p3mup6
     -   *F(4756)+p3mup6*p4mup4*F(4757)+p2mup4*p4mup6*F(4759)+p3mup4*
     -   p4mup6*F(4760)+p4mup4*p4mup6*F(4761)))+p1mup6*p2mup3*p2mup4*
     -   F(4764)+p1mup6*p2mup3*p3mup4*F(4765)+p1mup6*p2mup3*p4mup4*F(
     -   4766)-p1mup6*p2mup4*p4mup3*F(4767)+p3mup3*p3mup6*p4mup4*F(47
     -   68)+p1mup6*p2mup4*p3mup3*F(4770)+p1mup6*p3mup3*p3mup4*F(4771
     -   )+p1mup6*p3mup3*p4mup4*F(4772)-p1mup6*p3mup4*p4mup3*F(4773)-
     -   p2mup4*p4mup3*p4mup6*F(4774)-p3mup4*p4mup3*p4mup6*F(4775)+p1
     -   mup4*(-(p1mup3*(p2mup6*F(4750)+p3mup6*F(4754)+p4mup6*F(4758)
     -   ))+p4mup3*(p4mup6*F(4749)+p1mup6*F(4762))+p1mup6*(p2mup3*F(4
     -   763)+p3mup3*F(4769)+p5mup3*F(4776)))+p1mup6*p2mup4*p5mup3*F(
     -   4777)+p1mup6*p3mup4*p5mup3*F(4778)+p1mup6*p4mup4*p5mup3*F(47
     -   79)+p4mup3*p4mup4*p4mup6*F(4780)+p1mup6*p4mup3*p4mup4*F(4781
     -   )
       Fa34(476) = p5mup4*(-(p1mup3*(p1mup6*F(4782)+p2mup6*F(4783)+p3m
     -   up6*F(4784)))+p1mup6*p2mup3*F(4785)+p3mup3*p3mup6*F(4786)+p2
     -   mup3*p3mup6*F(4787)+p2mup6*p3mup3*F(4788)+p2mup3*p2mup6*F(47
     -   89)+p1mup6*p3mup3*F(4790))
       Return
       End
