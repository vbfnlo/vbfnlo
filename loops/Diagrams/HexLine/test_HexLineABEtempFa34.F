c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineAbe34(p1mup2,p1mup3,p1mup4,p1m
     -   up5,p2mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,
     -   p4mup2,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mu
     -   p2,p6mup3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,
     -   mup3mup5,mup4mup5,Fa34)
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
       Complex*16 Fa34(364:374)
       COMMON/test_FaFunctionsAbe/Fa
       Fa34(364) = mup4mup5*p6mup2*F(4794)-mup2mup5*p6mup4*F(4795)+mup
     -   2mup4*p6mup5*F(4796)
       Fa34(365) = mup2mup4*p1mup5*F(4797)-mup2mup4*p2mup5*F(4798)-mup
     -   2mup4*p3mup5*F(4799)-mup2mup4*p4mup5*F(4800)+2*(p1mup5*p4mup
     -   2*p6mup4*F(4773)-p1mup5*p6mup2*p6mup4*F(4775)+p4mup2*p4mup5*
     -   p6mup4*F(4784)+p4mup2*p6mup4*p6mup5*F(4787)-p4mup5*p6mup2*p6
     -   mup4*F(4788)-p6mup2*p6mup4*p6mup5*F(4789)+p2mup2*(-(p1mup5*p
     -   2mup4*F(4804))+p1mup4*(p2mup5*F(4803)+p1mup5*F(4810)+p4mup5*
     -   F(4811)-p6mup5*F(4812))+p3mup5*(p6mup4*F(4391)+p1mup4*F(4807
     -   )-p3mup4*F(4821))+p6mup4*(p2mup5*F(4819)-p1mup5*F(4824)+p4mu
     -   p5*F(4825)+p6mup5*F(4826))+p4mup4*p4mup5*F(4827))+p1mup2*(-(
     -   p4mup4*p6mup5*F(4689))+p4mup4*p4mup5*F(4719)+p4mup5*p6mup4*F
     -   (4764)+p6mup4*p6mup5*F(4770)+p1mup4*(p1mup5*F(4670)+p4mup5*F
     -   (4683)-p6mup5*F(4686)+p2mup5*F(4801)+p3mup5*F(4805))-p1mup5*
     -   (-(p4mup4*F(4676))+p6mup4*F(4769)+p2mup4*F(4802)+p3mup4*F(48
     -   06))-p2mup4*p6mup5*F(4808)-p3mup4*p6mup5*F(4809)+p2mup4*p3mu
     -   p5*F(4813)-p2mup5*p3mup4*F(4814)-p2mup4*p2mup5*F(4815)+p2mup
     -   4*p4mup5*F(4816)-p2mup5*p4mup4*F(4817)+p2mup5*p6mup4*F(4818)
     -   -p3mup4*p3mup5*F(4820)-p3mup5*p4mup4*F(4822)+p3mup5*p6mup4*F
     -   (4823)+p3mup4*p4mup5*F(4828))+p1mup5*p3mup2*p6mup4*F(4829)+p
     -   2mup5*p3mup2*p6mup4*F(4830)+p2mup5*p4mup2*p6mup4*F(4831)-p2m
     -   up5*p6mup2*p6mup4*F(4832)+p3mup2*p3mup5*p6mup4*F(4833)+p3mup
     -   5*p4mup2*p6mup4*F(4834)+p3mup2*p4mup5*p6mup4*F(4834)+p3mup2*
     -   p6mup4*p6mup5*F(4835)-p3mup5*p6mup2*p6mup4*F(4836))
       Fa34(366) = p2mup2*(p2mup4*(-(p3mup5*F(4348))+p2mup5*F(4703)+p6
     -   mup5*F(4838)-p4mup5*F(4842))+p4mup4*(p3mup5*F(4351)+p6mup5*F
     -   (4840)+p2mup5*F(4843)-p1mup5*F(4844))+p3mup4*(p1mup5*F(4837)
     -   +p6mup5*F(4839)+p2mup5*F(4841)-p4mup5*F(4845)))
       Fa34(367) = p3mup2*(p2mup4*p2mup5*F(4386)-p2mup5*p4mup4*F(4395)
     -   -p4mup4*p4mup5*F(4496)+p1mup4*p2mup5*F(4846)-p1mup5*p2mup4*F
     -   (4847)+p1mup4*p1mup5*F(4848)+p1mup5*p3mup4*F(4849)+p2mup5*p3
     -   mup4*F(4850)+p1mup4*p4mup5*F(4851)+p2mup4*p4mup5*F(4852)-p1m
     -   up4*p6mup5*F(4853)-p2mup4*p6mup5*F(4854)-p3mup4*p6mup5*F(485
     -   5)-p4mup4*p6mup5*F(4856)-p1mup4*p3mup5*F(4857)-p2mup4*p3mup5
     -   *F(4858)-p1mup5*p4mup4*F(4859)-p3mup4*p3mup5*F(4860)-p3mup5*
     -   p4mup4*F(4861)-p3mup4*p4mup5*F(4862))
       Fa34(368) = mup2mup5*p2mup4*F(4863)-mup2mup5*p3mup4*F(4864)+mup
     -   2mup5*p4mup4*F(4865)+2*(p1mup5*p4mup2*p4mup4*F(4728)-p4mup2*
     -   p4mup4*p4mup5*F(4731)-p4mup2*p4mup4*p6mup5*F(4752)-p4mup4*p6
     -   mup2*p6mup5*F(4768)-p1mup5*p4mup4*p6mup2*F(4774)-p4mup4*p4mu
     -   p5*p6mup2*F(4785)+p1mup5*p2mup4*p4mup2*F(4867)+p1mup4*(p4mup
     -   2*p4mup5*F(4727)+p1mup5*p4mup2*F(4765)-p4mup2*p6mup5*F(4766)
     -   +p4mup5*p6mup2*F(4767)-p2mup5*p4mup2*F(4866)+p2mup5*p6mup2*F
     -   (4868))-p1mup5*p2mup4*p6mup2*F(4869)+p2mup5*p3mup4*p4mup2*F(
     -   4870)+p1mup5*p3mup4*p4mup2*F(4871)-p2mup5*p3mup4*p6mup2*F(48
     -   73)-p1mup5*p3mup4*p6mup2*F(4875)+p2mup4*p4mup2*p4mup5*F(4876
     -   )+p2mup4*p2mup5*p4mup2*F(4877)+p2mup4*p2mup5*p6mup2*F(4878)-
     -   p2mup4*p4mup2*p6mup5*F(4879)+p2mup4*p4mup5*p6mup2*F(4880)-p2
     -   mup4*p6mup2*p6mup5*F(4881)-p3mup4*p4mup2*p6mup5*F(4882)-p3mu
     -   p4*p6mup2*p6mup5*F(4883)+p2mup5*p4mup2*p4mup4*F(4884)-p2mup5
     -   *p4mup4*p6mup2*F(4885)+p3mup5*(p4mup2*(p4mup4*F(4496)-p1mup4
     -   *F(4851)-p2mup4*F(4852)+p3mup4*F(4862))+p6mup2*(p2mup4*F(487
     -   2)+p1mup4*F(4874)-p3mup4*F(4886)-p4mup4*F(4887)))-p3mup4*p4m
     -   up2*p4mup5*F(4888)-p3mup4*p4mup5*p6mup2*F(4889))
       Fa34(369) = p1mup4*(2*(p1mup5*p6mup2*F(4771)-p6mup2*p6mup5*F(47
     -   72))+mup2mup5*F(4890))
       Fa34(370) = mup4mup5*F(4891)+p2mup5*(6*(p3mup4*F(4607)+p2mup4*F
     -   (4612)+p4mup4*F(4613))-2*(p1mup4*F(4892)+p6mup4*F(4902)))-2*
     -   (p1mup4*p3mup5*F(4894)+p1mup4*p1mup5*F(4896)+p1mup4*p4mup5*F
     -   (4897)-p1mup4*p6mup5*F(4898)-p3mup4*(3*(p3mup5*F(4626)-p4mup
     -   5*F(4637))+p1mup5*F(4895)+p6mup5*F(4900))+p2mup4*(-(p1mup5*F
     -   (4893))-p6mup5*F(4899)+3*(p4mup5*F(4635)-p3mup5*F(4903)))+p3
     -   mup5*p6mup4*F(4905)+p4mup4*(-(p6mup5*F(4901))+3*(p4mup5*F(46
     -   36)-p3mup5*F(4904))+p1mup5*F(4906))-p1mup5*p6mup4*F(4907)+p4
     -   mup5*p6mup4*F(4908)+p6mup4*p6mup5*F(4909))
       Fa34(371) = 4*(mup3mup5*F(4910)+2*(-(p1mup3*p4mup5*F(4911))-p2m
     -   up3*p4mup5*F(4915)+p1mup3*p2mup5*F(4916)+p1mup5*p2mup3*F(491
     -   7)+p4mup3*(p6mup5*F(4912)-p1mup5*F(4913)+p2mup5*F(4914)+p3mu
     -   p5*F(4918))+p4mup3*p4mup5*F(4919)+p2mup3*p3mup5*F(4920)-p2mu
     -   p5*p3mup3*F(4921)+p2mup3*p6mup5*F(4922)+p2mup5*p6mup3*F(4923
     -   )+p3mup3*p4mup5*F(4924)-p4mup5*p6mup3*F(4925)-p2mup3*p2mup5*
     -   F(4926)-(p1mup3*(p1mup5*F(4927)+p3mup5*F(4928)-p6mup5*F(4929
     -   ))+p1mup5*p3mup3*F(4930)-p3mup3*p3mup5*F(4931)+p3mup5*p6mup3
     -   *F(4932)-p1mup5*p6mup3*F(4933)+p3mup3*p6mup5*F(4934)-p6mup3*
     -   p6mup5*F(4935))*P(283)))
       Fa34(372) = -(mup3mup5*F(4936))+2*(p2mup3*p2mup5*F(4940)+p2mup3
     -   *p6mup5*F(4941)-p2mup5*p3mup3*F(4942)+p3mup3*p6mup5*F(4943)-
     -   p2mup3*p3mup5*F(4945)+p3mup3*p3mup5*F(4946)+p2mup3*p4mup5*F(
     -   4947)+p3mup3*p4mup5*F(4948)-p2mup5*p4mup3*F(4949)-p6mup3*p6m
     -   up5*F(4952)+p2mup5*p6mup3*F(4954)+p1mup3*(p3mup5*F(4937)+p2m
     -   up5*F(4944)+p1mup5*F(4950)+p6mup5*F(4951)+p4mup5*F(4955))+p4
     -   mup5*p6mup3*F(4956)-p3mup5*p6mup3*F(4957)-p4mup3*p6mup5*F(49
     -   58)-p1mup5*(p2mup3*F(4938)+p3mup3*F(4939)+p6mup3*F(4953)-p4m
     -   up3*F(4959))+p3mup5*p4mup3*F(4960)-p4mup3*p4mup5*F(4961))
       Fa34(373) = 4*mup2mup5*F(4962)+8*(p3mup2*p6mup5*F(4965)+p1mup2*
     -   p3mup5*F(4966)+p2mup2*p3mup5*F(4967)+p3mup2*p3mup5*F(4970)+p
     -   3mup2*p4mup5*F(4971)-p2mup2*p6mup5*F(4972)+p3mup5*p4mup2*F(4
     -   974)+p4mup2*p6mup5*F(4975)+p1mup2*p6mup5*F(4976)-p3mup5*p6mu
     -   p2*F(4978)+p1mup5*(p3mup2*F(4963)-p2mup2*F(4969)+p4mup2*F(49
     -   73)+p6mup2*F(4977)-p1mup2*F(4979))-p1mup2*p4mup5*F(4980)+p2m
     -   up2*p4mup5*F(4982)+p4mup2*p4mup5*F(4983)+p4mup5*p6mup2*F(498
     -   5)+p6mup2*p6mup5*F(4986)+p2mup5*(p3mup2*F(4964)+p1mup2*F(496
     -   8)+p2mup2*F(4981)+p6mup2*F(4984)+p4mup2*F(4987)))
       Fa34(374) = -4*mup3mup4*F(4988)+p4mup4*p6mup3*F(5000)-8*(p4mup3
     -   *p6mup4*F(4995)-p1mup4*p4mup3*F(4996)-p2mup4*p4mup3*F(4997)-
     -   p2mup4*p6mup3*F(4998)+p3mup4*p6mup3*F(4999)+p3mup4*p4mup3*F(
     -   5002)-p1mup3*(p2mup4*F(4989)+p3mup4*F(4990)+p4mup4*F(4991)+p
     -   1mup4*F(5003)-p6mup4*F(5004))-p2mup3*(-(p2mup4*F(4992))-p3mu
     -   p4*F(4993)+p4mup4*F(4994)+p1mup4*F(5001)+p6mup4*F(5005))+p6m
     -   up3*p6mup4*F(5006)-p1mup4*p6mup3*F(5007)+p4mup3*p4mup4*F(500
     -   8)-p2mup4*p3mup3*F(5009)+p1mup4*p3mup3*F(5010)+p3mup3*p3mup4
     -   *F(5011)+p3mup3*p4mup4*F(5012)-p3mup3*p6mup4*F(5013))
       Return
       End
