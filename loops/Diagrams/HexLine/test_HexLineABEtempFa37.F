c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineAbe37(p1mup2,p1mup3,p1mup4,p1m
     -   up5,p2mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,
     -   p4mup2,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mu
     -   p2,p6mup3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,
     -   mup3mup5,mup4mup5,Fa37)
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
       Complex*16 Fa37(397:407)
       COMMON/test_FaFunctionsAbe/Fa
       Fa37(397) = p6mup4*(p3mup2*p4mup3*F(5447)-p4mup2*(p2mup3*F(5446
     -   )+p3mup3*F(5448))+p4mup2*p4mup3*F(5449)+p3mup2*p6mup3*F(5454
     -   )-p6mup2*(p2mup3*F(5453)+p3mup3*F(5455))+p1mup2*(p4mup3*F(54
     -   45)-p1mup3*F(5450)-p2mup3*F(5451)-p3mup3*F(5457))+p1mup3*p3m
     -   up2*F(5458)+p2mup3*p3mup2*F(5460)+p2mup2*(p4mup3*F(5444)+p6m
     -   up3*F(5452)+p1mup3*F(5456)+p2mup3*F(5459)+p3mup3*F(5461))+p3
     -   mup2*p3mup3*F(5462))
       Fa37(398) = mup2mup3*p6mup4*F(5466)+mup2mup4*(p4mup3*F(5463)+p1
     -   mup3*F(5464)+p3mup3*F(5467)+p2mup3*F(5468)+p6mup3*F(5469))+m
     -   up3mup4*(p1mup2*F(5465)+p6mup2*F(5470)-p2mup2*F(5471))+2*(p2
     -   mup2*p2mup4*p6mup3*F(5479)+p2mup2*p3mup4*p6mup3*F(5480)+p2mu
     -   p2*p4mup4*p6mup3*F(5481)+p1mup3*p4mup2*p6mup4*F(5482)-p2mup2
     -   *p3mup4*p4mup3*F(5484)-p1mup3*p3mup2*p3mup4*F(5485)+p2mup3*p
     -   3mup2*p3mup4*F(5486)+p2mup2*p2mup4*p4mup3*F(5488)+p4mup3*p6m
     -   up2*p6mup4*F(5489)+p1mup2*(p1mup4*p4mup3*F(5472)+p1mup4*p6mu
     -   p3*F(5473)+p2mup4*p4mup3*F(5474)+p2mup4*p6mup3*F(5475)+p3mup
     -   4*p6mup3*F(5476)+p4mup4*p6mup3*F(5477)+p3mup4*p4mup3*F(5483)
     -   -p6mup3*p6mup4*F(5490)+p4mup3*p4mup4*F(5491))-p4mup2*p6mup3*
     -   p6mup4*F(5492)+p1mup3*p6mup2*p6mup4*F(5493)+p1mup4*(p2mup2*p
     -   6mup3*F(5478)+p2mup2*p4mup3*F(5487)+p3mup2*p3mup3*F(5494))+p
     -   2mup4*p3mup2*p3mup3*F(5495)+p3mup2*p3mup3*p3mup4*F(5496)-p2m
     -   up2*p4mup3*p4mup4*F(5497)+p6mup2*p6mup3*p6mup4*F(5498))
       Fa37(399) = mup2mup3*(p3mup4*F(5500)+8*(p2mup4*F(5499)+p4mup4*F
     -   (5502)+p1mup4*F(5503)))+8*mup3mup4*(p4mup2*F(5501)-p3mup2*F(
     -   5504))-16*(-(p3mup2*p3mup4*p6mup3*F(5507))+p3mup3*p3mup4*p6m
     -   up2*F(5508)+p2mup3*p2mup4*p6mup2*F(5511)+p2mup3*p3mup4*p6mup
     -   2*F(5512)+p2mup3*p4mup4*p6mup2*F(5513)-p2mup4*p3mup2*p6mup3*
     -   F(5514)+p2mup4*p3mup3*p6mup2*F(5515)-p3mup2*p4mup4*p6mup3*F(
     -   5520)+p3mup3*p4mup4*p6mup2*F(5521)-p2mup3*p2mup4*p3mup2*F(55
     -   27)-p2mup3*p2mup4*p4mup2*F(5528)-p2mup3*p3mup2*p4mup4*F(5529
     -   )+p2mup4*p3mup2*p4mup3*F(5530)-p2mup4*p3mup3*p4mup2*F(5531)+
     -   p1mup3*(p3mup4*p4mup2*F(5516)+p3mup4*p6mup2*F(5517)+p1mup4*p
     -   4mup2*F(5522)+p4mup4*p6mup2*F(5523)+p1mup4*p6mup2*F(5524)+p2
     -   mup4*(p4mup2*F(5505)+p6mup2*F(5509)+p3mup2*F(5525))+p1mup4*p
     -   3mup2*F(5532)+p3mup2*p4mup4*F(5533)+p4mup2*p4mup4*F(5535))-p
     -   4mup2*p4mup4*p6mup3*F(5537)+p2mup3*p3mup4*p4mup2*F(5538)-p2m
     -   up4*p4mup2*p6mup3*F(5539)-p3mup4*p4mup2*p6mup3*F(5540)-p2mup
     -   4*p4mup3*p6mup2*F(5542)+p2mup4*p6mup2*p6mup3*F(5543)+p3mup4*
     -   p6mup2*p6mup3*F(5544)+p4mup4*p6mup2*p6mup3*F(5546)-p3mup2*p3
     -   mup4*p4mup3*F(5548)+p3mup3*p3mup4*p4mup2*F(5549)+p4mup3*p4mu
     -   p4*p6mup2*F(5550)-p1mup4*(p3mup2*p6mup3*F(5518)+p2mup3*(p4mu
     -   p2*F(5506)-p6mup2*F(5510)+p3mup2*F(5526))+p3mup2*p4mup3*F(55
     -   34)+p4mup2*p4mup3*F(5536)+p4mup2*p6mup3*F(5541)-p4mup3*p6mup
     -   2*F(5545)-p6mup2*p6mup3*F(5547)-p3mup3*(p6mup2*F(5519)+p4mup
     -   2*F(5551)))+p3mup4*p4mup3*p6mup2*F(5552)-p2mup3*p4mup2*p4mup
     -   4*F(5553)+p2mup4*p4mup2*p4mup3*F(5554)-p3mup2*p3mup3*p4mup4*
     -   F(5555)-p3mup3*p4mup2*p4mup4*F(5556)+p3mup4*p4mup2*p4mup3*F(
     -   5557)+p3mup2*p4mup3*p4mup4*F(5558)-p4mup2*p4mup3*p4mup4*F(55
     -   59))
       Fa37(400) = mup4mup5*(-(p2mup2*p2mup3*F(5108))+p2mup2*p3mup3*F(
     -   5137)+p2mup2*p4mup3*F(5138)+p3mup3*p4mup2*F(5141)+p3mup2*p4m
     -   up3*F(5141)+p2mup3*p4mup2*F(5560)+p1mup2*(p3mup3*F(5106)-p1m
     -   up3*F(5117)-p2mup3*F(5118)+p4mup3*F(5561))-p1mup3*(p2mup2*F(
     -   5107)+p4mup2*F(5562)))
       Fa37(401) = mup4mup5*(p2mup3*(p3mup2*F(5111)-p6mup2*F(5124))+p1
     -   mup3*p3mup2*F(5130)+p1mup3*p6mup2*F(5132)+p3mup2*p3mup3*F(51
     -   40)-p3mup3*p6mup2*F(5142)+p4mup2*p4mup3*F(5143)+mup2mup3*F(5
     -   563)-p4mup3*p6mup2*F(5564)-p6mup3*(p1mup2*F(5121)+p2mup2*F(5
     -   139)+p3mup2*F(5142)+p4mup2*F(5564)-p6mup2*F(5565)))
       Fa37(402) = mup3mup4*(p2mup2*p2mup5*F(4623)+p4mup2*p4mup5*F(463
     -   5)+p3mup2*p3mup5*F(5133)+p1mup2*(-(p4mup5*F(4608))+p1mup5*F(
     -   4618)+p2mup5*F(4620)-p3mup5*F(5566)))
       Fa37(403) = mup3mup4*(p1mup2*p6mup5*F(4621)+p1mup5*p4mup2*F(464
     -   1)+p4mup5*p6mup2*F(4642)-p1mup5*p6mup2*F(4643)+p4mup2*p6mup5
     -   *F(4649)-p6mup2*p6mup5*F(4651)-p3mup5*p4mup2*F(4903)-p3mup2*
     -   p4mup5*F(4903)-p2mup5*(p4mup2*F(4612)-p6mup2*F(4625)+p3mup2*
     -   F(5109))+p2mup2*(-(p4mup5*F(4612))+p1mup5*F(4622)+p6mup5*F(4
     -   646)-p3mup5*F(5109))+p1mup5*p3mup2*F(5395)+p3mup5*p6mup2*F(5
     -   401)-mup2mup5*F(5567)+p3mup2*p6mup5*F(5568))
       Fa37(404) = mup2mup3*(p3mup4*p3mup5*F(4626)-p1mup5*p3mup4*F(463
     -   3)-p2mup4*p4mup5*F(4635)-p4mup4*p4mup5*F(4636)-p3mup4*p4mup5
     -   *F(4637)-p1mup5*p2mup4*F(4641)+p2mup4*p3mup5*F(4903)+p3mup5*
     -   p4mup4*F(4904)+p1mup4*p3mup5*F(5389)+p3mup4*p6mup5*F(5399)-p
     -   3mup5*p6mup4*F(5402)-p1mup4*p1mup5*F(5404)-p1mup5*p4mup4*F(5
     -   405)-p1mup4*p4mup5*F(5406)-p1mup4*p6mup5*F(5407)-p4mup5*p6mu
     -   p4*F(5410)+p1mup5*p6mup4*F(5411)-p4mup4*p6mup5*F(5418)+p2mup
     -   4*p6mup5*F(5569)+p2mup5*(p3mup4*F(4607)+p1mup4*F(4608)+p2mup
     -   4*F(4612)+p4mup4*F(4613)-p6mup4*F(5570))+p6mup4*p6mup5*F(557
     -   1))
       Fa37(405) = -(mup2mup5*(p3mup3*p3mup4*F(5388)-p3mup4*p4mup3*F(5
     -   392)+p3mup4*p6mup3*F(5403)-p3mup3*p4mup4*F(5412)+p3mup3*p6mu
     -   p4*F(5413)-p4mup3*p4mup4*F(5414)+p4mup4*p6mup3*F(5416)+p4mup
     -   3*p6mup4*F(5416)-p2mup3*p2mup4*F(5572)+p1mup3*(-(p4mup4*F(53
     -   86))+p1mup4*F(5393)+p3mup4*F(5394)+p6mup4*F(5400)-p2mup4*F(5
     -   574))-p2mup4*p3mup3*F(5577)+p1mup4*(p3mup3*F(5387)-p4mup3*F(
     -   5408)-p6mup3*F(5409)+p2mup3*F(5578))+p2mup3*p3mup4*F(5579)+p
     -   2mup3*p4mup4*F(5580)-p2mup4*p4mup3*F(5581)-p2mup3*p6mup4*F(5
     -   582)+p2mup4*p6mup3*F(5584)))+mup2mup4*(-(p3mup3*p4mup5*F(462
     -   6))+p1mup5*p4mup3*F(4633)+p4mup3*p4mup5*F(4637)+p4mup5*p6mup
     -   3*F(4638)-p2mup3*p4mup5*F(4903)+p1mup3*p1mup5*F(4927)+p1mup5
     -   *p3mup3*F(4930)-p1mup5*p6mup3*F(4933)+p3mup3*p6mup5*F(4934)+
     -   p3mup5*(-(p4mup3*F(4626))+p1mup3*F(4928)+p6mup3*F(4932)+p2mu
     -   p3*F(5133)+p3mup3*F(5134))-p1mup3*p4mup5*F(5389)+p1mup5*p2mu
     -   p3*F(5395)+p1mup3*p6mup5*F(5397)+p4mup3*p6mup5*F(5415)-p2mup
     -   3*p6mup5*F(5583)-p2mup5*(p4mup3*F(4607)+p2mup3*F(5109)+p3mup
     -   3*F(5112)+p1mup3*F(5566)-p6mup3*F(5585)))+mup3mup5*(p2mup2*p
     -   2mup4*F(4611)-p2mup4*p4mup2*F(4614)+p2mup4*p6mup2*F(4624)-p2
     -   mup2*p4mup4*F(4644)+p2mup2*p6mup4*F(4645)-p4mup2*p4mup4*F(46
     -   47)+p4mup4*p6mup2*F(4648)+p4mup2*p6mup4*F(4648)+p1mup2*(-(p4
     -   mup4*F(4609))+p6mup4*F(4616)+p1mup4*F(4617)+p2mup4*F(4619)-p
     -   3mup4*F(5573))-p2mup2*p3mup4*F(5575)-p2mup4*p3mup2*F(5576)+p
     -   1mup4*(p2mup2*F(4610)-p4mup2*F(4639)-p6mup2*F(4640)+p3mup2*F
     -   (5586))+p3mup2*p3mup4*F(5587)-p3mup4*p4mup2*F(5588)+p3mup4*p
     -   6mup2*F(5589)-p3mup2*p4mup4*F(5590)+p3mup2*p6mup4*F(5591))
       Fa37(406) = mup3mup5*p6mup2*p6mup4*F(4650)-mup2mup4*mup3mup5*F(
     -   5567)+mup2mup5*p6mup3*p6mup4*F(5592)+mup2mup4*p6mup3*p6mup5*
     -   F(5593)
       Fa37(407) = p1mup2*(p1mup4*p2mup3*F(5596)+p1mup4*p3mup3*F(5597)
     -   +p2mup3*p2mup4*F(5598)+p2mup4*p3mup3*F(5599)+p2mup3*p3mup4*F
     -   (5601)+p3mup3*p3mup4*F(5602)+p1mup3*(p1mup4*F(5594)+p2mup4*F
     -   (5595)+p3mup4*F(5600)-p4mup4*F(5606))-p2mup3*p4mup4*F(5607)-
     -   p3mup3*p4mup4*F(5608))+p2mup2*(p1mup3*(-(p3mup4*F(5603))+p1m
     -   up4*F(5609)+p2mup4*F(5610)-p4mup4*F(5615))+p2mup3*(-(p3mup4*
     -   F(5604))+p1mup4*F(5611)+p2mup4*F(5613)-p4mup4*F(5616))+p3mup
     -   3*(-(p3mup4*F(5605))+p1mup4*F(5612)+p2mup4*F(5614)-p4mup4*F(
     -   5617)))
       Return
       End
