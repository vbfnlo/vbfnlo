c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineNoAbe1g35(p1mup2,p1mup3,p1mup4
     -   ,p1mup6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3m
     -   up6,p4mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,
     -   p6mup2,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3m
     -   up4,mup3mup6,mup4mup6,Fa35)
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
       Complex*16 Fa35(477:490)
       COMMON/test_FaFunctionsNoAbe1g/Fa
       Fa35(477) = p4mup3*F(4791)+p5mup3*F(4792)+p2mup3*F(4793)+p3mup3
     -   *F(4794)+p1mup3*F(4795)
       Fa35(478) = mup4mup6*(-4*(p2mup3*F(4796)+p3mup3*F(4797)+p1mup3*
     -   F(4798))+p4mup3*F(4799)+p5mup3*F(4800))+8*mup3mup6*p5mup4*F(
     -   4801)
       Fa35(479) = mup3mup6*(p4mup4*F(4802)+p3mup4*F(4803))
       Fa35(480) = p5mup6*(mup3mup4*F(4804)-2*(p3mup3*p3mup4*F(4808)+p
     -   2mup3*p3mup4*F(4809)-p2mup4*p3mup3*F(4811)-p2mup3*p4mup4*F(4
     -   812)-p3mup3*p4mup4*F(4813)+p3mup4*p4mup3*F(4814)-p2mup3*p2mu
     -   p4*F(4815)-p2mup3*p5mup4*F(4819)-p3mup3*p5mup4*F(4820)+p3mup
     -   4*p5mup3*F(4821)+p1mup3*(p1mup4*F(4805)+p3mup4*F(4810)-p2mup
     -   4*F(4816)-p4mup4*F(4817)-p5mup4*F(4822))+p1mup4*(p2mup3*F(48
     -   06)+p3mup3*F(4807)+p4mup3*F(4818)+p5mup3*F(4823))+p4mup3*p4m
     -   up4*F(4824)-p2mup4*p5mup3*F(4825)-p4mup3*p5mup4*F(4826)-p4mu
     -   p4*p5mup3*F(4827)+p2mup4*p4mup3*F(4828)-p5mup3*p5mup4*F(4829
     -   )))
       Fa35(481) = mup3mup6*(p1mup4*F(4830)+p2mup4*F(4831))-p5mup4*(p2
     -   mup3*p4mup6*F(4832)+p3mup3*p4mup6*F(4833)+p1mup3*p4mup6*F(48
     -   34)+p4mup3*p4mup6*F(4835)+p2mup6*p4mup3*F(4836)+p3mup6*p4mup
     -   3*F(4837)+p1mup6*p4mup3*F(4838)+p4mup6*p5mup3*F(4839)+p2mup6
     -   *p5mup3*F(4840)+p3mup6*p5mup3*F(4841))
       Fa35(482) = mup3mup4*(p4mup6*F(4843)+4*(p3mup6*F(4842)-p1mup6*F
     -   (4844)-p2mup6*F(4845)))+8*(p2mup3*p3mup4*p3mup6*F(4855)+p2mu
     -   p6*p3mup3*p3mup4*F(4856)+p2mup3*p3mup4*p4mup6*F(4857)+p3mup3
     -   *p3mup4*p4mup6*F(4858)+p2mup3*p2mup6*p3mup4*F(4859)+p3mup3*p
     -   3mup4*p3mup6*F(4860)+p1mup3*p1mup6*(p1mup4*F(4846)-p2mup4*F(
     -   4853)-p3mup4*F(4854)-p4mup4*F(4861))-p2mup3*p2mup4*p3mup6*F(
     -   4862)-p2mup3*p3mup6*p4mup4*F(4863)-p2mup4*p2mup6*p3mup3*F(48
     -   64)-p2mup6*p3mup3*p4mup4*F(4865)-p2mup3*p2mup4*p4mup6*F(4867
     -   )-p2mup3*p4mup4*p4mup6*F(4868)+p2mup6*p3mup4*p4mup3*F(4869)-
     -   p2mup4*p3mup3*p4mup6*F(4871)-p3mup3*p4mup4*p4mup6*F(4872)+p3
     -   mup4*p3mup6*p4mup3*F(4873)-p2mup3*p2mup4*p2mup6*F(4874)-p2mu
     -   p3*p2mup6*p4mup4*F(4875)-p2mup4*p3mup3*p3mup6*F(4876)+p2mup6
     -   *p3mup4*p5mup3*F(4878)+p3mup4*p3mup6*p5mup3*F(4880)+p1mup4*(
     -   p2mup3*p3mup6*F(4847)+p2mup6*p3mup3*F(4848)+p2mup3*p4mup6*F(
     -   4849)+p3mup3*p4mup6*F(4850)+p2mup3*p2mup6*F(4851)+p3mup3*p3m
     -   up6*F(4852)+p2mup6*p4mup3*F(4866)+p3mup6*p4mup3*F(4870)+p2mu
     -   p6*p5mup3*F(4877)+p3mup6*p5mup3*F(4879)+p4mup6*p5mup3*F(4881
     -   ))+p3mup4*p4mup6*p5mup3*F(4882)-p2mup4*p3mup6*p4mup3*F(4883)
     -   -p2mup4*p3mup6*p5mup3*F(4884)-p2mup4*p2mup6*p5mup3*F(4885)-p
     -   2mup6*p4mup4*p5mup3*F(4886)-p3mup6*p4mup4*p5mup3*F(4887)-p2m
     -   up4*p4mup6*p5mup3*F(4888)-p4mup4*p4mup6*p5mup3*F(4889)-p2mup
     -   4*p2mup6*p4mup3*F(4890)-p2mup6*p4mup3*p4mup4*F(4891)-p3mup6*
     -   p4mup3*p4mup4*F(4892)-p1mup6*p5mup3*p5mup4*F(4893))
       Fa35(483) = p1mup6*p2mup3*p2mup4*F(4899)+p1mup6*p2mup4*p3mup3*F
     -   (4900)+p1mup6*p2mup3*p3mup4*F(4902)+p1mup6*p3mup3*p3mup4*F(4
     -   904)+p3mup3*p3mup6*p4mup4*F(4906)-p3mup4*p4mup3*p4mup6*F(490
     -   8)+p1mup6*p2mup3*p4mup4*F(4911)+p1mup6*p3mup3*p4mup4*F(4914)
     -   +p1mup3*(-(p1mup4*(p2mup6*F(4894)+p3mup6*F(4896)+p4mup6*F(48
     -   98)))-p2mup6*p3mup4*F(4901)-p3mup4*p3mup6*F(4903)-p3mup4*p4m
     -   up6*F(4905)+p2mup4*p2mup6*F(4909)+p2mup6*p4mup4*F(4910)+p2mu
     -   p4*p3mup6*F(4912)+p3mup6*p4mup4*F(4913)+p2mup4*p4mup6*F(4915
     -   )+p4mup4*p4mup6*F(4916))-p1mup6*p3mup4*p4mup3*F(4918)-p1mup4
     -   *(p4mup3*p4mup6*F(4907)+p1mup6*(p2mup3*F(4895)+p3mup3*F(4897
     -   )+p4mup3*F(4917)+p5mup3*F(4919)))-p1mup6*p3mup4*p5mup3*F(492
     -   0)+p1mup6*p2mup4*p4mup3*F(4921)+p1mup6*p2mup4*p5mup3*F(4922)
     -   +p1mup6*p4mup4*p5mup3*F(4923)+p2mup4*p4mup3*p4mup6*F(4924)+p
     -   4mup3*p4mup4*p4mup6*F(4925)+p1mup6*p4mup3*p4mup4*F(4926)
       Fa35(484) = p5mup4*(p3mup3*p3mup6*F(4927)+p2mup3*p3mup6*F(4928)
     -   +p1mup3*p1mup6*F(4930)+p2mup6*(p3mup3*F(4737)+p2mup3*F(4929)
     -   +p1mup3*F(4931))+p1mup6*p2mup3*F(4932)+p1mup3*p3mup6*F(4933)
     -   +p1mup6*p3mup3*F(4934))
       Fa35(485) = p2mup3*F(4935)+p3mup3*F(4936)+p5mup3*F(4937)+p4mup3
     -   *F(4938)+p1mup3*F(4939)
       Fa35(486) = p4mup3*p5mup4*F(4946)+p1mup3*p2mup4*F(4952)+p1mup3*
     -   p5mup4*F(4958)+p1mup4*p5mup3*F(4960)+p5mup3*p5mup4*F(4962)+4
     -   *(p1mup3*p4mup4*F(4947)-p3mup4*p5mup3*F(4948)-p4mup4*p5mup3*
     -   F(4949)-p2mup4*p5mup3*F(4950)+p1mup3*p3mup4*F(4951)-p1mup4*p
     -   4mup3*F(4953)-p4mup3*p4mup4*F(4954)-p3mup4*p4mup3*F(4955)-p2
     -   mup4*p4mup3*F(4956)+mup3mup4*F(4961)+p1mup3*p1mup4*F(4963))-
     -   p2mup3*(-(p5mup4*F(4959))+4*(p4mup4*F(4940)+p3mup4*F(4943)+p
     -   2mup4*F(4944)+p1mup4*F(4964)))+p3mup3*(p2mup4*F(4945)+p5mup4
     -   *F(4957)-4*(p4mup4*F(4941)+p3mup4*F(4942)+p1mup4*F(4965)))
       Fa35(487) = mup4mup6*(p3mup2*F(4966)+8*(p1mup2*F(4967)-p2mup2*F
     -   (4974)+p5mup2*F(4976)+p4mup2*F(4978)))+4*(mup2mup4*(p1mup6*F
     -   (4968)+p2mup6*F(4969)+p3mup6*F(4970)+p4mup6*F(4971)+p5mup6*F
     -   (4972))-2*(mup2mup6*(p1mup4*F(4973)+p2mup4*F(4975)-p5mup4*F(
     -   4977)-p4mup4*F(4979)+p3mup4*F(4980))+p2mup2*p3mup4*p5mup6*F(
     -   4987)+p1mup6*p2mup2*p3mup4*F(4992)+p2mup2*p2mup6*p3mup4*F(49
     -   93)+p2mup2*p3mup4*p3mup6*F(4994)+p2mup2*p3mup4*p4mup6*F(4995
     -   )+p2mup6*p3mup2*p3mup4*F(4997)+p3mup2*p3mup4*p4mup6*F(4998)+
     -   p3mup2*p3mup4*p5mup6*F(4999)+p3mup2*p3mup4*p3mup6*F(5000)-p1
     -   mup6*p3mup2*p3mup4*F(5001)+p2mup6*p3mup4*p4mup2*F(5002)+p3mu
     -   p4*p4mup2*p4mup6*F(5004)+p3mup4*p3mup6*p4mup2*F(5006)+p3mup4
     -   *p4mup2*p5mup6*F(5009)+p1mup6*p3mup4*p4mup2*F(5013)-p3mup2*p
     -   3mup6*p5mup4*F(5014)-p2mup6*p3mup2*p5mup4*F(5015)-p1mup6*p3m
     -   up2*p5mup4*F(5016)-p3mup2*p4mup6*p5mup4*F(5017)-p3mup2*p5mup
     -   4*p5mup6*F(5018)+p2mup6*p3mup4*p5mup2*F(5019)+p3mup4*p3mup6*
     -   p5mup2*F(5020)+p3mup4*p4mup6*p5mup2*F(5021)+p3mup4*p5mup2*p5
     -   mup6*F(5022)+p1mup6*p3mup4*p5mup2*F(5023)+p1mup2*(p3mup4*(p5
     -   mup6*F(4986)+p1mup6*F(4988)+p2mup6*F(4989)+p3mup6*F(4990)+p4
     -   mup6*F(4991))-2*p3mup6*(p4mup4*F(5024)+p5mup4*F(5028)+p2mup4
     -   *F(5031)))+p1mup4*(p2mup6*p3mup2*F(4981)+p3mup2*(p4mup6*F(49
     -   82)+p5mup6*F(4983)+p3mup6*F(4984)+p1mup6*F(4985))+2*p3mup6*(
     -   p2mup2*F(5025)+p5mup2*F(5030)+p4mup2*F(5032)-p1mup2*F(5033))
     -   )-2*(p1mup6*p2mup4*p3mup2*F(4996)+p2mup4*p2mup6*p3mup2*F(500
     -   2)+p2mup6*p3mup2*p4mup4*F(5003)+p2mup4*p3mup2*p4mup6*F(5004)
     -   +p3mup2*p4mup4*p4mup6*F(5005)+p2mup4*p3mup2*p5mup6*F(5007)+p
     -   3mup2*p4mup4*p5mup6*F(5008)+p2mup4*p3mup2*p3mup6*F(5010)+p3m
     -   up2*p3mup6*p4mup4*F(5011)+p1mup6*p3mup2*p4mup4*F(5012)-p2mup
     -   2*p2mup4*p3mup6*F(5026)+p2mup4*p3mup6*p4mup2*F(5027)+p2mup4*
     -   p3mup6*p5mup2*F(5029)+p2mup2*p3mup6*p4mup4*F(5034)+p2mup2*p3
     -   mup6*p5mup4*F(5035)+p3mup6*p4mup2*p4mup4*F(5036)+p3mup6*p4mu
     -   p4*p5mup2*F(5037)+p3mup6*p4mup2*p5mup4*F(5037)+p3mup6*p5mup2
     -   *p5mup4*F(5038))))
       Fa35(488) = p3mup2*F(5039)+4*(p2mup2*F(5040)+p4mup2*F(5041)+p5m
     -   up2*F(5042))+p1mup2*F(5043)
       Fa35(489) = mup4mup6*(p1mup2*F(5044)+p2mup2*F(5045)+p3mup2*F(50
     -   46)+p4mup2*F(5047)+p5mup2*F(5048))-2*mup2mup6*p5mup4*F(5049)
       Fa35(490) = p5mup6*(p1mup2*(p4mup4*F(4817)+p5mup4*F(4822)+p3mup
     -   4*F(5050))+p2mup2*(p3mup4*F(5051)+p4mup4*F(5052)+p5mup4*F(50
     -   53)))
       Return
       End
