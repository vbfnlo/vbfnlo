c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineNoAbe1g36(p1mup2,p1mup3,p1mup4
     -   ,p1mup6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3m
     -   up6,p4mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,
     -   p6mup2,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3m
     -   up4,mup3mup6,mup4mup6,Fa36)
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
       Complex*16 Fa36(491:504)
       COMMON/test_FaFunctionsNoAbe1g/Fa
       Fa36(491) = p1mup2*(p1mup4*p5mup6*F(4805)-p4mup6*p5mup4*F(4834)
     -   -p1mup6*p4mup4*F(4861)-p4mup4*p4mup6*F(4916)-p1mup6*p5mup4*F
     -   (4930)+p2mup4*p5mup6*F(5055)-p1mup6*p3mup4*F(5057)-p2mup6*p3
     -   mup4*F(5059)-p3mup4*p3mup6*F(5060)-p3mup4*p4mup6*F(5061)-p2m
     -   up6*p4mup4*F(5065)-p3mup6*p4mup4*F(5066)-p2mup6*p5mup4*F(507
     -   1)-p3mup6*p5mup4*F(5072))+p2mup2*(p1mup4*p5mup6*F(5054)+p2mu
     -   p4*p5mup6*F(5056)-p1mup6*p3mup4*F(5058)-p2mup6*p3mup4*F(5062
     -   )-p3mup4*p3mup6*F(5063)-p3mup4*p4mup6*F(5064)-p1mup6*p4mup4*
     -   F(5067)-p2mup6*p4mup4*F(5068)-p3mup6*p4mup4*F(5069)-p4mup4*p
     -   4mup6*F(5070)-p1mup6*p5mup4*F(5073)-p2mup6*p5mup4*F(5074)-p3
     -   mup6*p5mup4*F(5075)-p4mup6*p5mup4*F(5076))
       Fa36(492) = p1mup4*p2mup2*p2mup6*F(5080)+p1mup4*p2mup2*p3mup6*F
     -   (5081)+p1mup4*p2mup2*p4mup6*F(5082)+p1mup2*(p1mup4*p4mup6*F(
     -   4898)+p1mup4*p2mup6*F(5077)+p1mup4*p3mup6*F(5078)+p2mup4*p2m
     -   up6*F(5084)+p2mup4*p3mup6*F(5085)+p2mup4*p4mup6*F(5086))+p2m
     -   up2*p2mup4*p2mup6*F(5088)+p2mup2*p2mup4*p3mup6*F(5089)+p2mup
     -   2*p2mup4*p4mup6*F(5090)-p3mup2*p3mup6*p5mup4*F(5091)-p2mup6*
     -   p3mup2*p5mup4*F(5092)-p1mup6*(p4mup2*p5mup4*F(4838)-p1mup2*p
     -   1mup4*F(4846)+p5mup2*p5mup4*F(4893)-p1mup4*p2mup2*F(5079)+p1
     -   mup2*p2mup4*F(5083)+p2mup2*p2mup4*F(5087)+p3mup2*p5mup4*F(50
     -   93))-p3mup2*p4mup6*p5mup4*F(5094)-p3mup2*p5mup4*p5mup6*F(509
     -   5)-p2mup6*p4mup2*p5mup4*F(5096)-p3mup6*p4mup2*p5mup4*F(5097)
     -   -p2mup6*p5mup2*p5mup4*F(5098)-p3mup6*p5mup2*p5mup4*F(5099)
       Fa36(493) = 2*(p4mup2*p5mup4*(p5mup6*F(4826)+p4mup6*F(4835))+p5
     -   mup2*p5mup4*(p5mup6*F(4829)+p4mup6*F(4839)))+mup2mup4*(p1mup
     -   6*F(5100)+p2mup6*F(5101)+p3mup6*F(5102)+p4mup6*F(5103)+p5mup
     -   6*F(5104))
       Fa36(494) = -(p4mup4*p5mup2*p5mup6*F(4827))-p4mup4*p4mup6*p5mup
     -   2*F(4889)-p1mup6*p4mup4*p5mup2*F(4923)+mup2mup6*p2mup4*F(510
     -   6)+mup2mup6*p3mup4*F(5107)-mup2mup6*p4mup4*F(5108)+p2mup4*p2
     -   mup6*p3mup2*F(5110)+p2mup4*p3mup2*p4mup6*F(5112)+p2mup4*p3mu
     -   p2*p5mup6*F(5114)+p2mup4*p3mup2*p3mup6*F(5116)+p1mup6*p2mup4
     -   *p3mup2*F(5118)+p2mup4*p2mup6*p5mup2*F(5127)+p1mup4*(-(mup2m
     -   up6*F(5105))+p3mup2*(p2mup6*F(5109)+p4mup6*F(5111)+p5mup6*F(
     -   5113)+p3mup6*F(5115)+p1mup6*F(5117))+p4mup2*(p5mup6*F(4818)+
     -   p4mup6*F(4907)+p1mup6*F(4917)+p2mup6*F(5119)+p3mup6*F(5121))
     -   +p5mup2*(p5mup6*F(4823)+p4mup6*F(4881)+p1mup6*F(4919)+p2mup6
     -   *F(5126)+p3mup6*F(5128)))+p2mup4*p3mup6*p5mup2*F(5129)+p2mup
     -   4*p4mup6*p5mup2*F(5130)+p2mup4*p5mup2*p5mup6*F(5131)+p1mup6*
     -   p2mup4*p5mup2*F(5132)+p2mup6*p3mup2*p3mup4*F(5133)-p2mup6*p3
     -   mup2*p4mup4*F(5134)-p3mup2*p3mup4*p4mup6*F(5136)-p3mup2*p4mu
     -   p4*p4mup6*F(5137)-p3mup2*p3mup4*p5mup6*F(5138)-p3mup2*p4mup4
     -   *p5mup6*F(5139)-p3mup2*p3mup4*p3mup6*F(5140)-p3mup2*p3mup6*p
     -   4mup4*F(5141)+p1mup6*p3mup2*p3mup4*F(5142)-p1mup6*p3mup2*p4m
     -   up4*F(5143)-p2mup6*p3mup4*p5mup2*F(5145)-p1mup6*p3mup4*p5mup
     -   2*F(5146)-p2mup6*p4mup4*p5mup2*F(5147)-p3mup4*p3mup6*p5mup2*
     -   F(5148)-p3mup6*p4mup4*p5mup2*F(5149)-p3mup4*p4mup6*p5mup2*F(
     -   5151)-p3mup4*p5mup2*p5mup6*F(5152)-p4mup2*(-(p2mup4*(p2mup6*
     -   F(5120)+p3mup6*F(5122)+p5mup6*F(5123)+p4mup6*F(5124)+p1mup6*
     -   F(5125)))+p4mup4*(-(p5mup6*F(4824))+p4mup6*F(4925)+p1mup6*F(
     -   4926)+p2mup6*F(5153)+p3mup6*F(5155))+p3mup4*(p2mup6*F(5135)+
     -   p1mup6*F(5144)-p5mup6*F(5150)+p3mup6*F(5154)+p4mup6*F(5156))
     -   )
       Fa36(495) = p2mup2*F(5157)+p4mup2*F(5158)+p5mup2*F(5159)+p3mup2
     -   *F(5160)+p1mup2*F(5161)
       Fa36(496) = p3mup2*(p3mup4*F(5163)+p5mup4*F(5175)+4*(p4mup4*F(5
     -   162)+p2mup4*F(5164)+p1mup4*F(5183)))-p1mup2*(-(p5mup4*F(5176
     -   ))+8*(p4mup4*F(5165)+p3mup4*F(5168)+p2mup4*F(5170)+p1mup4*F(
     -   5184)))+8*(-(mup2mup4*F(5181))+p2mup2*(-(p5mup4*F(5166))+p1m
     -   up4*F(5171)+p4mup4*F(5178)+p3mup4*F(5179)-p2mup4*F(5185))-p4
     -   mup2*(p5mup4*F(5167)-p1mup4*F(5169)+p2mup4*F(5180)-p3mup4*F(
     -   5182)+p4mup4*F(5186))+p5mup2*(p2mup4*F(5172)+p3mup4*F(5173)+
     -   p4mup4*F(5174)+p1mup4*F(5177)-p5mup4*F(5187)))
       Fa36(497) = p2mup2*p5mup4*F(5203)+4*(p2mup2*p4mup4*F(5190)+p4mu
     -   p2*p5mup4*F(5192)-p2mup2*p3mup4*F(5194)+p1mup4*p4mup2*F(5195
     -   )-p2mup2*p2mup4*F(5198)-p2mup4*p4mup2*F(5202)-p3mup4*p4mup2*
     -   F(5204)+p1mup4*p2mup2*F(5205)-p4mup2*p4mup4*F(5206)-p5mup2*p
     -   5mup4*F(5209))+p3mup2*(p5mup4*F(5207)+8*(p4mup4*F(5188)+p3mu
     -   p4*F(5191)+p2mup4*F(5196)+p1mup4*F(5211)))+8*(p2mup4*p5mup2*
     -   F(5199)+p3mup4*p5mup2*F(5200)+p4mup4*p5mup2*F(5201)-mup2mup4
     -   *F(5210)+p1mup4*p5mup2*F(5212))+p1mup2*(p5mup4*F(5208)+4*(p4
     -   mup4*F(5189)-p3mup4*F(5193)+p2mup4*F(5197)+p1mup4*F(5213)))
       Fa36(498) = mup4mup6*(p4mup2*F(5214)+p2mup2*F(5215)+p1mup2*F(52
     -   16)+p3mup2*F(5217)+p5mup2*F(5218))+4*mup2mup6*p5mup4*F(5219)
       Fa36(499) = p5mup6*(-(p1mup2*(p3mup4*F(5220)+p4mup4*F(5221)+p5m
     -   up4*F(5224)))+p2mup2*(p3mup4*F(5222)+p4mup4*F(5223)+p5mup4*F
     -   (5225)))
       Fa36(500) = p1mup2*(p1mup4*p5mup6*F(5226)+p2mup4*p5mup6*F(5227)
     -   +p1mup6*p3mup4*F(5228)+p2mup6*p3mup4*F(5229)+p3mup4*p3mup6*F
     -   (5230)+p3mup4*p4mup6*F(5231)+p1mup6*p4mup4*F(5232)+p2mup6*p4
     -   mup4*F(5233)+p3mup6*p4mup4*F(5234)+p4mup4*p4mup6*F(5235)+p1m
     -   up6*p5mup4*F(5246)+p2mup6*p5mup4*F(5247)+p3mup6*p5mup4*F(524
     -   8)+p4mup6*p5mup4*F(5249))-p2mup2*(p1mup4*p5mup6*F(5236)+p2mu
     -   p4*p5mup6*F(5237)+p1mup6*p3mup4*F(5238)+p2mup6*p3mup4*F(5239
     -   )+p3mup4*p3mup6*F(5240)+p3mup4*p4mup6*F(5241)+p1mup6*p4mup4*
     -   F(5242)+p2mup6*p4mup4*F(5243)+p3mup6*p4mup4*F(5244)+p4mup4*p
     -   4mup6*F(5245)+p1mup6*p5mup4*F(5250)+p2mup6*p5mup4*F(5251)+p3
     -   mup6*p5mup4*F(5252)+p4mup6*p5mup4*F(5253))
       Fa36(501) = -(p1mup2*(p1mup4*(p1mup6*F(5254)+p2mup6*F(5255)+p3m
     -   up6*F(5256)+p4mup6*F(5257))+p2mup4*(p1mup6*F(5258)+p2mup6*F(
     -   5259)+p3mup6*F(5260)+p4mup6*F(5261))))+p1mup4*p1mup6*p2mup2*
     -   F(5262)+p1mup4*p2mup2*p2mup6*F(5263)+p1mup4*p2mup2*p3mup6*F(
     -   5264)+p1mup4*p2mup2*p4mup6*F(5265)+p1mup6*p2mup2*p2mup4*F(52
     -   66)+p2mup2*p2mup4*p2mup6*F(5267)+p2mup2*p2mup4*p3mup6*F(5268
     -   )+p2mup2*p2mup4*p4mup6*F(5269)-p2mup6*p4mup2*p5mup4*F(5270)-
     -   p3mup6*p4mup2*p5mup4*F(5271)-p1mup6*p4mup2*p5mup4*F(5272)+p3
     -   mup2*p5mup4*(p3mup6*F(4787)+p2mup6*F(4789)+p1mup6*F(5273)+p4
     -   mup6*F(5274)+p5mup6*F(5275))+p2mup6*p5mup2*p5mup4*F(5276)+p3
     -   mup6*p5mup2*p5mup4*F(5277)+p1mup6*p5mup2*p5mup4*F(5278)
       Fa36(502) = mup2mup4*(p1mup6*F(5279)+4*(p2mup6*F(5280)+p3mup6*F
     -   (5281)+p4mup6*F(5282)+p5mup6*F(5283)))+8*p5mup4*(-(p4mup2*(p
     -   4mup6*F(5284)+p5mup6*F(5285)))+p5mup2*(p4mup6*F(5286)+p5mup6
     -   *F(5287)))
       Fa36(503) = mup2mup6*p3mup4*F(5289)+8*(p3mup2*p3mup4*p3mup6*F(4
     -   712)+p3mup2*p3mup6*p4mup4*F(4713)+p2mup6*p3mup2*p3mup4*F(471
     -   4)+p2mup6*p3mup2*p4mup4*F(4727)+p3mup6*p4mup4*p5mup2*F(5075)
     -   +mup2mup6*p4mup4*F(5290)+p3mup2*p3mup4*p4mup6*F(5300)+p3mup2
     -   *p4mup4*p4mup6*F(5301)+p3mup2*p3mup4*p5mup6*F(5304)+p3mup2*p
     -   4mup4*p5mup6*F(5305)+p1mup6*p3mup2*p3mup4*F(5308)+p1mup6*p3m
     -   up2*p4mup4*F(5309)-p1mup6*p3mup4*p4mup2*F(5310)-p2mup6*p3mup
     -   4*p4mup2*F(5312)-p3mup4*p3mup6*p4mup2*F(5314)-p3mup4*p4mup2*
     -   p5mup6*F(5316)-p3mup4*p4mup2*p4mup6*F(5318)+p2mup6*p3mup4*p5
     -   mup2*F(5321)+p2mup6*p4mup4*p5mup2*F(5322)+p3mup4*p3mup6*p5mu
     -   p2*F(5325)+p3mup4*p4mup6*p5mup2*F(5328)+p4mup4*p4mup6*p5mup2
     -   *F(5329)+p3mup4*p5mup2*p5mup6*F(5332)+p4mup4*p5mup2*p5mup6*F
     -   (5333)+p1mup4*(mup2mup6*F(5291)+p2mup6*p4mup2*F(5292)+p3mup6
     -   *p4mup2*F(5293)+p4mup2*p5mup6*F(5294)+p4mup2*p4mup6*F(5295)+
     -   p1mup6*p4mup2*F(5296)+p3mup2*(p3mup6*F(4709)+p2mup6*F(4725)+
     -   p4mup6*F(5298)+p5mup6*F(5302)+p1mup6*F(5306))+p2mup6*p5mup2*
     -   F(5319)+p3mup6*p5mup2*F(5323)+p4mup6*p5mup2*F(5326)+p5mup2*p
     -   5mup6*F(5330)+p1mup6*p5mup2*F(5334))+p2mup4*(mup2mup6*F(5288
     -   )-p1mup6*p4mup2*F(5297)+p3mup2*(p3mup6*F(4711)+p2mup6*F(4726
     -   )+p4mup6*F(5299)+p5mup6*F(5303)+p1mup6*F(5307))-p2mup6*p4mup
     -   2*F(5311)-p3mup6*p4mup2*F(5313)-p4mup2*p5mup6*F(5315)-p4mup2
     -   *p4mup6*F(5317)+p2mup6*p5mup2*F(5320)+p3mup6*p5mup2*F(5324)+
     -   p4mup6*p5mup2*F(5327)+p5mup2*p5mup6*F(5331)+p1mup6*p5mup2*F(
     -   5335))+p1mup6*p3mup4*p5mup2*F(5336)+p1mup6*p4mup4*p5mup2*F(5
     -   337)+p2mup6*p4mup2*p4mup4*F(5338)+p3mup6*p4mup2*p4mup4*F(533
     -   9)+p4mup2*p4mup4*p5mup6*F(5340)+p4mup2*p4mup4*p4mup6*F(5341)
     -   +p1mup6*p4mup2*p4mup4*F(5342))
       Fa36(504) = mup4mup6*F(5343)-4*(p1mup4*(p5mup6*F(5344)+p1mup6*F
     -   (5345))+p1mup4*p3mup6*F(5347)+p1mup4*p4mup6*F(5348)-p1mup6*p
     -   3mup4*F(5350)-p3mup4*p5mup6*F(5351)-p4mup4*p5mup6*F(5352)-p1
     -   mup6*p4mup4*F(5357)-3*(p4mup4*p4mup6*F(5005)+p3mup4*p3mup6*F
     -   (5355)+p3mup4*p4mup6*F(5356)+p3mup6*p4mup4*F(5358))-p2mup4*(
     -   p1mup6*F(5349)+p5mup6*F(5353)+3*(p2mup6*F(5002)+p4mup6*F(500
     -   4)+p3mup6*F(5359)))-p5mup4*p5mup6*F(5360)-p1mup6*p5mup4*F(53
     -   61)-p2mup6*(-(p1mup4*F(5346))+3*(p4mup4*F(5003)+p3mup4*F(535
     -   4))+p5mup4*F(5362))-p3mup6*p5mup4*F(5363)-p4mup6*p5mup4*F(53
     -   64))
       Return
       End
