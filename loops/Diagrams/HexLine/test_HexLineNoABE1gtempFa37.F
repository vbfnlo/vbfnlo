c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineNoAbe1g37(p1mup2,p1mup3,p1mup4
     -   ,p1mup6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3m
     -   up6,p4mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,
     -   p6mup2,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3m
     -   up4,mup3mup6,mup4mup6,Fa37)
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
       Complex*16 Fa37(505:518)
       COMMON/test_FaFunctionsNoAbe1g/Fa
       Fa37(505) = p4mup6*F(5365)+2*p2mup6*F(5366)+p5mup6*F(5367)+p1mu
     -   p6*F(5368)+p3mup6*F(5369)
       Fa37(506) = p4mup2*F(5370)+p2mup2*F(5371)+p5mup2*F(5372)+p3mup2
     -   *F(5373)+p1mup2*F(5374)
       Fa37(507) = p4mup4*F(5376)+p2mup4*F(5377)+p3mup4*F(5378)+p5mup4
     -   *F(5379)+p1mup4*F(5380)
       Fa37(508) = p3mup2*(p4mup4*F(5381)+p5mup4*F(5382)+p3mup4*F(5383
     -   )-p1mup4*F(5384)-p2mup4*F(5385))-2*(p4mup2*p5mup4*F(5388)-p1
     -   mup4*p4mup2*F(5390)-p5mup2*(p2mup4*F(5393)+p3mup4*F(5394)+p4
     -   mup4*F(5395))-p1mup4*p5mup2*F(5397)+mup2mup4*F(5398)-p3mup4*
     -   p4mup2*F(5400)+p2mup4*p4mup2*F(5401)+p1mup2*(p4mup4*F(5386)+
     -   p3mup4*F(5389)+p2mup4*F(5391)+p5mup4*F(5396)+p1mup4*F(5403))
     -   +p2mup2*(p5mup4*F(5387)-p1mup4*F(5392)+p4mup4*F(5399)+p3mup4
     -   *F(5402)+p2mup4*F(5404))+p4mup2*p4mup4*F(5405)+p5mup2*p5mup4
     -   *F(5406))
       Fa37(509) = p5mup4*F(5407)-p4mup4*F(5408)+p2mup4*F(5409)+p1mup4
     -   *F(5410)-p3mup4*F(5411)
       Fa37(510) = p3mup4*p5mup3*F(5417)+p4mup4*p5mup3*F(5418)+p2mup4*
     -   p5mup3*F(5419)+p2mup4*p3mup3*F(5423)+p1mup3*(p3mup4*F(5420)-
     -   4*(p4mup4*F(5413)+p2mup4*F(5424)))+mup3mup4*F(5433)+p5mup3*(
     -   p5mup4*F(5429)+p1mup4*F(5434))-4*(p3mup3*p4mup4*F(5415)+p3mu
     -   p3*p3mup4*F(5416)+p2mup3*(p4mup4*F(5414)+p3mup4*F(5421)+p2mu
     -   p4*F(5422))-p4mup3*(p5mup4*F(5412)-p1mup4*F(5425)+p2mup4*F(5
     -   426)+p3mup4*F(5427)-p4mup4*F(5428))+p3mup3*p5mup4*F(5430)+p1
     -   mup3*p5mup4*F(5431)+p2mup3*p5mup4*F(5432)+p1mup3*p1mup4*F(54
     -   35)+p1mup4*p2mup3*F(5436)+p1mup4*p3mup3*F(5437))
       Fa37(511) = p2mup3*p3mup4*F(5443)+p1mup4*p2mup3*F(5446)+p4mup4*
     -   p5mup3*F(5452)+p5mup3*p5mup4*F(5457)-4*(p2mup3*p5mup4*F(5439
     -   )+p2mup3*p4mup4*F(5440)+p2mup3*p2mup4*F(5444)+p3mup3*(p5mup4
     -   *F(5438)+p4mup4*F(5441)-p3mup4*F(5442)+p2mup4*F(5445)-p1mup4
     -   *F(5447))+p1mup3*p4mup4*F(5448)+p4mup3*p4mup4*F(5449)+p3mup4
     -   *p4mup3*F(5450)+p3mup4*p5mup3*F(5451)+p2mup4*p5mup3*F(5453)-
     -   p1mup3*p3mup4*F(5454)+p2mup4*p4mup3*F(5455)+p1mup3*p2mup4*F(
     -   5456)+p4mup3*p5mup4*F(5458)+p1mup3*p5mup4*F(5459)-mup3mup4*F
     -   (5460)+p1mup3*p1mup4*F(5461)+p1mup4*p5mup3*F(5462)+p1mup4*p4
     -   mup3*F(5463))
       Fa37(512) = mup2mup6*p3mup3*F(5465)+mup2mup6*p1mup3*F(5468)+4*(
     -   mup3mup6*p3mup2*F(5466)+mup2mup3*p3mup6*F(5470))+mup2mup3*p4
     -   mup6*F(5471)+mup2mup6*p4mup3*F(5472)+mup2mup6*p5mup3*F(5473)
     -   +p2mup3*(mup2mup6*F(5464)+8*(p4mup2*p4mup6*F(4867)+p2mup2*p3
     -   mup6*F(5478)+p2mup2*p4mup6*F(5480)-p3mup2*p4mup6*F(5482)-p3m
     -   up2*p3mup6*F(5484)+p3mup6*p4mup2*F(5486)+p3mup6*p5mup2*F(549
     -   2)+p4mup6*p5mup2*F(5495)))+p2mup6*(-4*mup2mup3*F(5469)+8*(p4
     -   mup2*p5mup3*F(4885)+p4mup2*p4mup3*F(4890)+p2mup2*p5mup3*F(52
     -   51)+p5mup2*p5mup3*F(5276)+p2mup3*(p4mup2*F(4874)+p2mup2*F(52
     -   67)+p5mup2*F(5320)-p3mup2*F(5483))-p3mup2*p4mup3*F(5487)-p3m
     -   up2*p5mup3*F(5493)+p3mup3*(p4mup2*F(4864)+p2mup2*F(5479)+p5m
     -   up2*F(5494)+p3mup2*F(5499))+p2mup2*p4mup3*F(5503)+p4mup3*p5m
     -   up2*F(5507)))+8*(p3mup3*p4mup2*p4mup6*F(4871)+p4mup2*p4mup6*
     -   p5mup3*F(4888)+p4mup2*p4mup3*p4mup6*F(4924)+mup3mup6*p1mup2*
     -   F(5467)+mup3mup6*p2mup2*F(5474)+mup3mup6*p5mup2*F(5475)+mup3
     -   mup6*p4mup2*F(5476)+p2mup2*p3mup3*p4mup6*F(5477)+p2mup2*p3mu
     -   p3*p3mup6*F(5481)-p3mup2*p4mup3*p4mup6*F(5485)-p3mup2*p3mup6
     -   *p4mup3*F(5488)+p3mup3*p3mup6*p4mup2*F(5489)-p3mup2*p3mup6*p
     -   5mup3*F(5490)+p3mup3*p3mup6*p5mup2*F(5491)-p3mup2*p4mup6*p5m
     -   up3*F(5496)+p3mup3*p4mup6*p5mup2*F(5497)+p3mup2*p3mup3*p3mup
     -   6*F(5498)+p3mup2*p3mup3*p4mup6*F(5500)+p2mup2*p4mup3*p4mup6*
     -   F(5501)+p2mup2*p3mup6*p4mup3*F(5502)+p2mup2*p3mup6*p5mup3*F(
     -   5504)+p2mup2*p4mup6*p5mup3*F(5505)+p4mup3*p4mup6*p5mup2*F(55
     -   06)+p3mup6*p4mup2*p5mup3*F(5508)+p3mup6*p4mup3*p5mup2*F(5509
     -   )+p3mup6*p4mup2*p4mup3*F(5510)+p3mup6*p5mup2*p5mup3*F(5511)+
     -   p4mup6*p5mup2*p5mup3*F(5512))
       Fa37(513) = 4*(-2*(p1mup2*(p2mup3*p3mup6*F(5516)+p2mup3*p4mup6*
     -   F(5521)+p3mup3*p3mup6*F(5522)+p3mup3*p4mup6*F(5524)+p1mup6*(
     -   p5mup3*F(5246)+p2mup3*F(5258)+p3mup3*F(5515)+p4mup3*F(5533))
     -   +p2mup6*(p5mup3*F(5247)+p2mup3*F(5259)+p3mup3*F(5517)+p4mup3
     -   *F(5534))+p3mup6*p4mup3*F(5536)+p4mup3*p4mup6*F(5538)+p3mup6
     -   *p5mup3*F(5541)+p4mup6*p5mup3*F(5543))-p1mup3*(p1mup6*p2mup2
     -   *F(5262)+p1mup6*p5mup2*F(5334)-p1mup2*(p1mup6*F(5254)+p3mup6
     -   *F(5514))+p1mup6*p3mup2*F(5518)+p2mup6*(p4mup2*F(4909)-p1mup
     -   2*F(5255)+p2mup2*F(5263)+p5mup2*F(5319)+p3mup2*F(5519))+p3mu
     -   p2*p3mup6*F(5523)+p2mup2*p3mup6*F(5528)+p4mup2*(p4mup6*F(491
     -   5)+p1mup6*F(5535)+p3mup6*F(5537))+p3mup6*p5mup2*F(5542)+p4mu
     -   p6*(-(p1mup2*F(5520))+p3mup2*F(5525)+p2mup2*F(5529)+p5mup2*F
     -   (5544))))+p1mup6*(-(mup2mup3*F(5513))+2*(-(p2mup3*p3mup2*F(5
     -   527))+p2mup3*p4mup2*F(5530)-p3mup2*p4mup3*F(5531)+p3mup3*p4m
     -   up2*F(5532)-p3mup2*p5mup3*F(5539)+p2mup2*(p5mup3*F(5250)+p2m
     -   up3*F(5266)+p3mup3*F(5526)+p4mup3*F(5545))+p4mup2*p5mup3*F(5
     -   546)+p5mup2*(p5mup3*F(5278)+p2mup3*F(5335)+p3mup3*F(5540)+p4
     -   mup3*F(5547))+p3mup2*p3mup3*F(5548)+p4mup2*p4mup3*F(5549))))
       Fa37(514) = p5mup6*(-(p1mup2*(p1mup3*F(5226)+p2mup3*F(5227)+p3m
     -   up3*F(5550)))+p2mup2*(p1mup3*F(5236)+p2mup3*F(5237)+p3mup3*F
     -   (5551)))
       Fa37(515) = p5mup6*(-(p4mup3*(p1mup2*F(5555)+p3mup2*F(5556)-p2m
     -   up2*F(5557)))+p4mup2*(p2mup3*F(5552)+p3mup3*F(5553)+p1mup3*F
     -   (5554)+p4mup3*F(5558)+p5mup3*F(5559)))
       Fa37(516) = 4*p5mup6*(mup2mup3*F(5560)-2*(p1mup2*p5mup3*F(5224)
     -   -p2mup2*p5mup3*F(5225)-p3mup2*(p1mup3*F(5561)-p2mup3*F(5562)
     -   -p5mup3*F(5563)+p3mup3*F(5565))-p5mup2*(p5mup3*F(5287)+p1mup
     -   3*F(5330)+p2mup3*F(5331)+p3mup3*F(5564)+p4mup3*F(5566))))
       Fa37(517) = mup3mup6*p1mup2*F(5567)+mup3mup6*p3mup2*F(5574)+mup
     -   3mup6*p4mup2*F(5575)+mup3mup6*p5mup2*F(5576)+p2mup2*(mup3mup
     -   6*F(5568)-8*(-(p4mup3*p4mup6*F(5064))-p4mup6*p5mup3*F(5497)+
     -   p2mup3*p2mup6*F(5577)+p2mup3*p4mup6*F(5579)-p3mup3*p4mup6*F(
     -   5580)-p2mup6*p3mup3*F(5582)-p2mup6*p4mup3*F(5587)-p2mup6*p5m
     -   up3*F(5592)))+4*(mup2mup3*p2mup6*F(5570)+mup2mup3*p4mup6*F(5
     -   572)+p3mup6*(mup2mup3*F(5571)+2*(p2mup2*p4mup3*F(5063)+p3mup
     -   3*p4mup2*F(5140)+p3mup2*p4mup3*F(5140)+p4mup2*p4mup3*F(5154)
     -   +p2mup2*p3mup3*F(5584)-p2mup3*(p4mup2*F(4732)+p2mup2*F(5578)
     -   +p3mup2*F(5585))+p2mup2*p5mup3*F(5589)+p5mup2*(p5mup3*F(4698
     -   )+p3mup3*F(4786)+p4mup3*F(5148)-p2mup3*F(5590))+p3mup2*p5mup
     -   3*F(5598)+p4mup2*p5mup3*F(5606)+p3mup2*p3mup3*F(5607))))+8*(
     -   p3mup3*p4mup2*p4mup6*F(5136)+p3mup2*p4mup3*p4mup6*F(5136)+p4
     -   mup3*p4mup6*p5mup2*F(5151)+p4mup2*p4mup3*p4mup6*F(5156)+mup2
     -   mup6*p1mup3*F(5569)-mup2mup6*p2mup3*F(5573)-p2mup3*p3mup2*p4
     -   mup6*F(5581)-p2mup3*p2mup6*p3mup2*F(5583)-p2mup3*p4mup2*p4mu
     -   p6*F(5586)-p2mup3*p2mup6*p4mup2*F(5588)-p2mup3*p4mup6*p5mup2
     -   *F(5591)-p2mup3*p2mup6*p5mup2*F(5593)+mup2mup6*p3mup3*F(5594
     -   )+mup2mup6*p4mup3*F(5595)+mup2mup6*p5mup3*F(5596)+p2mup6*p3m
     -   up3*p4mup2*F(5597)+p2mup6*p3mup2*p4mup3*F(5597)+p2mup6*p3mup
     -   2*p5mup3*F(5599)+p2mup6*p3mup3*p5mup2*F(5600)+p3mup2*p4mup6*
     -   p5mup3*F(5601)+p3mup3*p4mup6*p5mup2*F(5602)+p4mup2*p4mup6*p5
     -   mup3*F(5603)+p2mup6*p4mup2*p5mup3*F(5604)+p2mup6*p4mup3*p5mu
     -   p2*F(5605)+p2mup6*p3mup2*p3mup3*F(5608)+p3mup2*p3mup3*p4mup6
     -   *F(5609)+p2mup6*p4mup2*p4mup3*F(5610)+p2mup6*p5mup2*p5mup3*F
     -   (5611)+p4mup6*p5mup2*p5mup3*F(5612))
       Fa37(518) = 8*(p1mup2*(p3mup6*p5mup3*F(4738)+p1mup6*p5mup3*F(47
     -   76)+p3mup6*p4mup3*F(5060)+p4mup3*p4mup6*F(5061)-p1mup3*(p1mu
     -   p6*F(4703)+p3mup6*F(4754)+p2mup6*F(5614)+p4mup6*F(5618))+p1m
     -   up6*p2mup3*F(5621)-p2mup3*p2mup6*F(5622)-p2mup3*p3mup6*F(562
     -   3)-p2mup3*p4mup6*F(5624)+p3mup3*(p3mup6*F(4728)+p1mup6*F(476
     -   9)+p2mup6*F(5627)+p4mup6*F(5631))+p1mup6*p4mup3*F(5635)+p2mu
     -   p6*p4mup3*F(5636)+p2mup6*p5mup3*F(5643)+p4mup6*p5mup3*F(5645
     -   ))-p1mup3*(p2mup6*(p2mup2*F(5616)+p3mup2*F(5629)+p4mup2*F(56
     -   38))+p3mup6*(p5mup2*F(4784)+p2mup2*F(5617)+p3mup2*F(5630)+p4
     -   mup2*F(5639))+p2mup6*p5mup2*F(5644)+p4mup6*(p2mup2*F(5619)+p
     -   3mup2*F(5632)+p4mup2*F(5640)+p5mup2*F(5646))))+p1mup6*(mup2m
     -   up3*F(5613)+8*(-(p2mup3*p3mup2*F(5626))-p2mup3*p4mup2*F(5634
     -   )-p1mup3*(p2mup2*F(5615)+p3mup2*F(5628)+p4mup2*F(5637))+p2mu
     -   p2*(p2mup3*F(5620)+p3mup3*F(5625)+p4mup3*F(5633)+p5mup3*F(56
     -   41))+p3mup3*p4mup2*F(5647)+p3mup2*p4mup3*F(5647)+p3mup2*p5mu
     -   p3*F(5648)+p4mup2*p5mup3*F(5649)+p5mup2*(p5mup3*F(4748)-p1mu
     -   p3*F(4782)+p3mup3*F(4790)-p2mup3*F(5642)+p4mup3*F(5650))+p3m
     -   up2*p3mup3*F(5651)+p4mup2*p4mup3*F(5652)))
       Return
       End
