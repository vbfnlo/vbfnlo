c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineAbe35(p1mup2,p1mup3,p1mup4,p1m
     -   up5,p2mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,
     -   p4mup2,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mu
     -   p2,p6mup3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,
     -   mup3mup5,mup4mup5,Fa35)
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
       Complex*16 Fa35(375:385)
       COMMON/test_FaFunctionsAbe/Fa
       Fa35(375) = 4*mup3mup4*F(5014)+8*(-(p2mup3*p2mup4*F(5018))+p2mu
     -   p4*p3mup3*F(5019)+p2mup3*p6mup4*F(5020)+p3mup3*p6mup4*F(5021
     -   )-p2mup3*p3mup4*F(5023)+p3mup3*p3mup4*F(5024)-p2mup4*p4mup3*
     -   F(5026)+p2mup4*p6mup3*F(5027)+p3mup4*p6mup3*F(5028)-p4mup4*p
     -   6mup3*F(5029)+p1mup3*(p3mup4*F(5015)+p2mup4*F(5022)-p4mup4*F
     -   (5025)+p1mup4*F(5030)+p6mup4*F(5031))-p6mup3*p6mup4*F(5032)+
     -   p2mup3*p4mup4*F(5034)+p3mup3*p4mup4*F(5035)-p4mup3*p6mup4*F(
     -   5036)+p1mup4*(-(p2mup3*F(5016))+p3mup3*F(5017)-p6mup3*F(5033
     -   )+p4mup3*F(5037))+p3mup4*p4mup3*F(5038)-p4mup3*p4mup4*F(5039
     -   ))
       Fa35(376) = -(mup2mup4*F(5040))+2*(p1mup2*p3mup4*F(5043)-p2mup2
     -   *p3mup4*F(5044)+p3mup2*p6mup4*F(5045)+p3mup2*p3mup4*F(5048)+
     -   p1mup2*p4mup4*F(5049)-p2mup2*p6mup4*F(5050)+p3mup4*p4mup2*F(
     -   5053)+p3mup4*p6mup2*F(5054)-p4mup2*p6mup4*F(5055)+p4mup4*p6m
     -   up2*F(5056)+p1mup2*p6mup4*F(5057)+p3mup2*p4mup4*F(5059)+p2mu
     -   p2*p4mup4*F(5060)+p1mup4*(p3mup2*F(5041)-p2mup2*F(5047)+p4mu
     -   p2*F(5051)+p6mup2*F(5058)-p1mup2*F(5062))+p2mup4*(p3mup2*F(5
     -   042)+p1mup2*F(5046)+p6mup2*F(5052)-p4mup2*F(5061)-p2mup2*F(5
     -   063))+p6mup2*p6mup4*F(5064)+p4mup2*p4mup4*F(5065))
       Fa35(377) = p3mup4*F(5066)+p1mup4*F(5067)+p6mup4*F(5068)+p2mup4
     -   *F(5069)-p4mup4*F(5070)
       Fa35(378) = p3mup5*F(5071)+p1mup5*F(5072)-p6mup5*F(5073)+p2mup5
     -   *F(5074)-p4mup5*F(5075)
       Fa35(379) = mup2mup5*(p4mup3*F(5076)-p1mup3*F(5078)-p3mup3*F(50
     -   82)+p2mup3*F(5086)+p6mup3*F(5087))+mup3mup5*(p4mup2*F(5077)+
     -   p1mup2*F(5079)+p3mup2*F(5083)+p6mup2*F(5088))+mup2mup3*(p1mu
     -   p5*F(5080)-p4mup5*F(5081)+p3mup5*F(5084)+p6mup5*F(5085)+p2mu
     -   p5*F(5090))+p2mup2*(-4*p2mup5*p4mup3*F(4623)+mup3mup5*F(5089
     -   )+2*(p4mup3*(p6mup5*F(5091)+2*p3mup5*F(5109)+p4mup5*F(5114)-
     -   p1mup5*F(5122))+p4mup5*(p1mup3*F(5107)+p2mup3*F(5108)-p3mup3
     -   *F(5137)+p6mup3*F(5139))))+2*(-(p1mup5*p2mup3*p4mup2*F(5096)
     -   )-p2mup3*p4mup2*p6mup5*F(5097)+p3mup2*p4mup3*p6mup5*F(5098)-
     -   p3mup3*p4mup2*p6mup5*F(5099)-p1mup5*p3mup2*p4mup3*F(5110)-p2
     -   mup3*p3mup2*p4mup5*F(5111)+p4mup2*p4mup3*p6mup5*F(5113)-p1mu
     -   p5*p4mup2*p4mup3*F(5115)-p2mup3*p4mup2*p4mup5*F(5116)-p4mup3
     -   *p6mup2*p6mup5*F(5119)+p2mup5*p4mup2*p6mup3*F(5120)+p1mup2*(
     -   p4mup3*(p6mup5*F(5092)-p1mup5*F(5093)-p2mup5*F(5094)-p3mup5*
     -   F(5101)-p4mup5*F(5104))+p4mup5*(-(p3mup3*F(5106))+p1mup3*F(5
     -   117)+p2mup3*F(5118)+p6mup3*F(5121)))+p1mup5*p4mup3*p6mup2*F(
     -   5123)+p2mup3*p4mup5*p6mup2*F(5124)-p2mup5*p4mup3*p6mup2*F(51
     -   25)-p3mup3*p4mup2*p4mup5*F(5126)+p3mup2*p4mup3*p4mup5*F(5126
     -   )+p4mup2*p6mup3*p6mup5*F(5127)+p3mup5*p4mup2*p6mup3*F(5128)-
     -   p1mup5*p4mup2*p6mup3*F(5129)+p1mup5*p3mup3*p4mup2*F(5131)+p1
     -   mup3*(p2mup5*p4mup2*F(5095)+p4mup2*(-(p6mup5*F(5100))+p1mup5
     -   *F(5102)+p3mup5*F(5103)+p4mup5*F(5105))-p4mup5*(p3mup2*F(513
     -   0)+p6mup2*F(5132)))-p3mup5*p4mup3*p6mup2*F(5135)-p4mup3*p4mu
     -   p5*p6mup2*F(5136)+p4mup2*p4mup5*p6mup3*F(5136)-p3mup2*p3mup3
     -   *p4mup5*F(5140)-2*(p2mup3*p2mup5*p4mup2*F(5109)-p2mup5*p3mup
     -   2*p4mup3*F(5109)+p2mup5*p3mup3*p4mup2*F(5112)-p2mup3*p3mup5*
     -   p4mup2*F(5133)+p3mup2*p3mup5*p4mup3*F(5133)-p3mup3*p3mup5*p4
     -   mup2*F(5134)-p2mup5*p4mup2*p4mup3*F(5138)-p3mup5*p4mup2*p4mu
     -   p3*F(5141))+p3mup3*p4mup5*p6mup2*F(5142)+p3mup2*p4mup5*p6mup
     -   3*F(5142)+p4mup2*p4mup3*p4mup5*F(5143)+p4mup5*p6mup2*p6mup3*
     -   F(5144))
       Fa35(380) = p1mup2*(p1mup5*p2mup3*F(4279)+p2mup3*p2mup5*F(4282)
     -   -p2mup3*p6mup5*F(4283)+p1mup5*p3mup3*F(5145)+p2mup5*p3mup3*F
     -   (5146)-p3mup3*p6mup5*F(5147)+p2mup3*p3mup5*F(5149)+p3mup3*p3
     -   mup5*F(5150)+p1mup3*(p1mup5*F(4276)+p2mup5*F(4277)-p6mup5*F(
     -   4278)+p3mup5*F(5148)+p4mup5*F(5151))+p2mup3*p4mup5*F(5152)+p
     -   3mup3*p4mup5*F(5153))+p2mup2*(p1mup5*p2mup3*F(4314)+p2mup3*p
     -   2mup5*F(4336)-p1mup3*(p6mup5*F(4288)-p2mup5*F(4307)-p1mup5*F
     -   (4313)+p4mup5*F(4674)+p3mup5*F(5154))-p2mup3*p3mup5*F(5155)+
     -   p1mup5*p3mup3*F(5156))
       Fa35(381) = p2mup2*(p2mup3*(p6mup5*F(4328)+p4mup5*F(4703))+p3mu
     -   p3*(p4mup5*F(4705)+p6mup5*F(5157)+p3mup5*F(5158)-p2mup5*F(51
     -   59)))
       Fa35(382) = mup2mup5*p2mup3*F(5160)+mup2mup5*p3mup3*F(5161)-mup
     -   2mup5*p1mup3*F(5162)+mup2mup5*p4mup3*F(5164)-mup2mup5*p6mup3
     -   *F(5165)+mup2mup3*p6mup5*F(5167)+mup3mup5*p6mup2*F(5168)+p2m
     -   up2*(-(mup3mup5*F(5166))+2*(p6mup3*p6mup5*F(4302)-p4mup3*p6m
     -   up5*F(5176)))+p1mup2*(mup3mup5*F(5163)+2*(p6mup3*p6mup5*F(43
     -   06)-p4mup3*p6mup5*F(5177)))-2*(p6mup2*p6mup3*p6mup5*F(4325)-
     -   p2mup3*p6mup2*p6mup5*F(4377)-p1mup3*p6mup2*p6mup5*F(4399)+p2
     -   mup3*p4mup2*p6mup5*F(5169)-p3mup2*p4mup3*p6mup5*F(5170)+p3mu
     -   p3*p4mup2*p6mup5*F(5171)+p1mup3*p4mup2*p6mup5*F(5172)+p3mup2
     -   *p3mup3*p6mup5*F(5173)+p2mup3*p3mup2*p6mup5*F(5174)+p1mup3*p
     -   3mup2*p6mup5*F(5175)-p3mup3*p6mup2*p6mup5*F(5178)-p4mup2*p4m
     -   up3*p6mup5*F(5179)-p4mup2*p6mup3*p6mup5*F(5180)-p3mup2*p6mup
     -   3*p6mup5*F(5181)-p4mup3*p6mup2*p6mup5*F(5182))
       Fa35(383) = -(mup3mup5*p3mup2*F(5183))+mup3mup5*p4mup2*F(5184)-
     -   2*(p2mup2*p2mup5*p6mup3*F(4318)-p2mup5*p6mup2*p6mup3*F(4323)
     -   +p2mup5*p4mup2*p6mup3*F(4777)+p2mup2*p4mup5*p6mup3*F(4777)+p
     -   4mup2*p4mup5*p6mup3*F(4786)-p2mup5*p3mup2*p6mup3*F(5185)-p3m
     -   up2*p3mup5*p6mup3*F(5187)-p3mup2*p4mup5*p6mup3*F(5189)-p1mup
     -   2*p6mup3*(p2mup5*F(4301)+p1mup5*F(4305)+p3mup5*F(5186)+p4mup
     -   5*F(5190))+p1mup5*p6mup3*(p2mup2*F(4315)-p6mup2*F(4317)+p3mu
     -   p2*F(5188)+p4mup2*F(5191))-p2mup2*p3mup5*p6mup3*F(5192)-p3mu
     -   p5*p4mup2*p6mup3*F(5193)+p3mup5*p6mup2*p6mup3*F(5194)-p4mup5
     -   *p6mup2*p6mup3*F(5195))
       Fa35(384) = mup2mup3*(p1mup5*F(5196)+p2mup5*F(5197)+p3mup5*F(51
     -   98)+p4mup5*F(5199))+2*(p1mup2*p4mup3*(p1mup5*F(5200)+p2mup5*
     -   F(5201)+p3mup5*F(5202)-p4mup5*F(5203))-p2mup2*p4mup3*(p4mup5
     -   *F(4776)+p1mup5*F(5204)+p2mup5*F(5205)-p3mup5*F(5206)))
       Fa35(385) = p2mup3*p2mup5*p6mup2*F(4375)+p2mup3*p3mup2*p4mup5*F
     -   (4386)-p2mup3*p2mup5*p4mup2*F(4703)-p2mup5*p3mup3*p4mup2*F(4
     -   705)+p2mup3*p4mup2*p4mup5*F(4720)+p3mup3*p4mup2*p4mup5*F(472
     -   2)+p4mup2*p4mup3*p4mup5*F(4730)-p2mup5*p4mup2*p4mup3*F(4776)
     -   -p3mup3*p3mup5*p4mup2*F(4821)+p2mup3*p2mup5*p3mup2*F(5213)+p
     -   2mup3*p3mup2*p3mup5*F(5214)+p2mup3*p3mup5*p4mup2*F(5218)+p2m
     -   up5*p3mup2*p4mup3*F(5219)-p2mup3*p3mup5*p6mup2*F(5221)+p2mup
     -   5*p3mup3*p6mup2*F(5222)+p3mup2*p3mup5*p4mup3*F(5225)+p3mup3*
     -   p3mup5*p6mup2*F(5227)+p1mup3*(p2mup5*p6mup2*F(4373)-p1mup5*p
     -   6mup2*F(4398)-p2mup5*p4mup2*F(4674)+p4mup2*p4mup5*F(4718)-p2
     -   mup5*p3mup2*F(5207)-p1mup5*p3mup2*F(5209)-p3mup2*p3mup5*F(52
     -   10)-p3mup2*p4mup5*F(5211)+p1mup5*p4mup2*F(5215)-p3mup5*p4mup
     -   2*F(5223)+p3mup5*p6mup2*F(5226)+p4mup5*p6mup2*F(5228))+p2mup
     -   3*p4mup5*p6mup2*F(5229)-p3mup2*p4mup3*p4mup5*F(5230)+p3mup3*
     -   p4mup5*p6mup2*F(5231)+p1mup5*(-(p2mup3*p6mup2*F(4374))+p2mup
     -   3*p3mup2*F(5208)+p3mup2*p3mup3*F(5212)+p2mup3*p4mup2*F(5216)
     -   +p3mup3*p4mup2*F(5217)-p3mup3*p6mup2*F(5220)-p3mup2*p4mup3*F
     -   (5224)-p4mup2*p4mup3*F(5232)+p4mup3*p6mup2*F(5233))-p2mup5*p
     -   3mup2*p3mup3*F(5234)-p3mup2*p3mup3*p3mup5*F(5235)-p3mup2*p3m
     -   up3*p4mup5*F(5236)+p3mup5*p4mup2*p4mup3*F(5237)+p2mup5*p4mup
     -   3*p6mup2*F(5238)-p3mup5*p4mup3*p6mup2*F(5239)+p4mup3*p4mup5*
     -   p6mup2*F(5240)
       Return
       End
