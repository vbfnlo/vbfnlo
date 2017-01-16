c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineAbe38(p1mup2,p1mup3,p1mup4,p1m
     -   up5,p2mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,
     -   p4mup2,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mu
     -   p2,p6mup3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,
     -   mup3mup5,mup4mup5,Fa38)
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
       Complex*16 Fa38(408:418)
       COMMON/test_FaFunctionsAbe/Fa
       Fa38(408) = p6mup4*(p1mup3*p2mup2*F(5621)+p2mup2*p2mup3*F(5622)
     -   +p2mup2*p3mup3*F(5623)+p2mup3*p4mup2*F(5624)-p3mup2*p4mup3*F
     -   (5625)+p3mup3*p4mup2*F(5626)+p3mup2*p3mup3*F(5627)+p2mup3*p3
     -   mup2*F(5628)+p1mup3*p3mup2*F(5629)+p2mup2*p4mup3*F(5630)+p1m
     -   up2*(p1mup3*F(5618)+p2mup3*F(5619)+p3mup3*F(5620)+p4mup3*F(5
     -   631))-p2mup2*p6mup3*F(5632)-p2mup3*p6mup2*F(5633)-p3mup3*p6m
     -   up2*F(5634)-p4mup2*p4mup3*F(5635)-p3mup2*p6mup3*F(5636))
       Fa38(409) = mup2mup4*(-(p2mup3*F(5637))+p3mup3*F(5638)+p1mup3*F
     -   (5639))+mup3mup4*p1mup2*F(5640)+mup2mup4*(-(p4mup3*F(5641))+
     -   p6mup3*F(5642))+mup2mup3*p6mup4*F(5644)+mup3mup4*p6mup2*F(56
     -   45)+2*(p1mup3*p3mup2*p3mup4*F(5646)+p1mup3*p4mup2*p6mup4*F(5
     -   647)-p1mup4*p3mup2*p3mup3*F(5648)-p2mup3*p3mup2*p3mup4*F(564
     -   9)-p1mup2*(p1mup4*p4mup3*F(5650)+p2mup4*p4mup3*F(5651)+p3mup
     -   4*p4mup3*F(5652)+p6mup3*p6mup4*F(5653)-p1mup4*p6mup3*F(5655)
     -   -p2mup4*p6mup3*F(5656)-p3mup4*p6mup3*F(5657)+p4mup3*p4mup4*F
     -   (5661)-p4mup4*p6mup3*F(5662))+p3mup2*p3mup3*(p2mup4*F(5663)+
     -   p3mup4*F(5667))-p6mup4*(p1mup3*p6mup2*F(5654)+p4mup2*p6mup3*
     -   F(5658)+p4mup3*p6mup2*F(5668)-p6mup2*p6mup3*F(5669)))+p2mup2
     -   *(mup3mup4*F(5643)+2*(p2mup4*p4mup3*F(5616)+p1mup4*p4mup3*F(
     -   5659)-p1mup4*p6mup3*F(5660)-p3mup4*p4mup3*F(5664)-p2mup4*p6m
     -   up3*F(5665)+p3mup4*p6mup3*F(5666)+p4mup3*p4mup4*F(5670)-p4mu
     -   p4*p6mup3*F(5671)))
       Fa38(410) = mup3mup4*(p3mup2*F(5672)+p4mup2*F(5676))+mup2mup3*(
     -   p1mup4*F(5673)+p2mup4*F(5674)+p3mup4*F(5675)+p4mup4*F(5677))
     -   +2*(p3mup2*p3mup4*p6mup3*F(5681)-p3mup3*p3mup4*p6mup2*F(5682
     -   )-p2mup3*p2mup4*p6mup2*F(5685)-p2mup3*p3mup4*p6mup2*F(5686)-
     -   p2mup3*p4mup4*p6mup2*F(5687)+p2mup4*p3mup2*p6mup3*F(5688)-p2
     -   mup4*p3mup3*p6mup2*F(5689)+p3mup2*p4mup4*p6mup3*F(5693)-p3mu
     -   p3*p4mup4*p6mup2*F(5694)-p2mup3*p2mup4*p3mup2*F(5697)+p2mup3
     -   *p2mup4*p4mup2*F(5700)+p2mup3*p3mup2*p4mup4*F(5701)-p2mup3*p
     -   3mup4*p4mup2*F(5702)-p2mup4*p3mup2*p4mup3*F(5703)+p2mup4*p3m
     -   up3*p4mup2*F(5704)-p3mup2*p3mup4*p4mup3*F(5707)+p3mup3*p3mup
     -   4*p4mup2*F(5708)+p4mup3*p4mup4*p6mup2*F(5712)+p2mup4*p4mup3*
     -   p6mup2*F(5713)+p2mup4*p6mup2*p6mup3*F(5714)+p3mup4*p4mup3*p6
     -   mup2*F(5715)+p3mup4*p6mup2*p6mup3*F(5716)+p4mup4*p6mup2*p6mu
     -   p3*F(5718)+p1mup4*(p3mup2*p6mup3*F(5691)-p2mup3*(p3mup2*F(56
     -   79)+p6mup2*F(5684)+p4mup2*F(5699))-p3mup3*(p6mup2*F(5692)+p4
     -   mup2*F(5705))+p3mup2*p4mup3*F(5711)+p4mup3*p6mup2*F(5717)+p6
     -   mup2*p6mup3*F(5719)+p4mup2*p4mup3*F(5720)-p4mup2*p6mup3*F(57
     -   21))+p3mup2*p3mup3*p4mup4*F(5722)+p2mup4*p4mup2*p4mup3*F(572
     -   3)+p2mup3*p4mup2*p4mup4*F(5723)+p3mup3*p4mup2*p4mup4*F(5724)
     -   -p3mup4*p4mup2*p4mup3*F(5725)+p1mup3*(p1mup4*p3mup2*F(5680)-
     -   p6mup2*(p3mup4*F(5690)+p4mup4*F(5695))-p1mup4*p6mup2*F(5696)
     -   +p2mup4*(p3mup2*F(5678)-p6mup2*F(5683)+p4mup2*F(5698))-p1mup
     -   4*p4mup2*F(5706)-p3mup2*p4mup4*F(5709)+p3mup4*p4mup2*F(5710)
     -   +p4mup2*p4mup4*F(5726))-p2mup4*p4mup2*p6mup3*F(5727)+p3mup4*
     -   p4mup2*p6mup3*F(5728)-p3mup2*p4mup3*p4mup4*F(5729)-p4mup2*p4
     -   mup3*p4mup4*F(5730)-p4mup2*p4mup4*p6mup3*F(5731))
       Fa38(411) = p1mup2*(-(p1mup4*p2mup3*F(5733))+p2mup3*(-(p2mup4*F
     -   (5740))+p3mup4*F(5741)+p4mup4*F(5746))+p1mup3*(p2mup4*F(5732
     -   )+p1mup4*F(5736)+p3mup4*F(5737)-p4mup4*F(5748))+p3mup3*(-(p2
     -   mup4*F(5742))+p1mup4*F(5749)+p3mup4*F(5750)-p4mup4*F(5754)))
     -   +p2mup2*(-(p1mup4*p2mup3*F(5735))+p2mup3*(-(p2mup4*F(5743))+
     -   p3mup4*F(5744)+p4mup4*F(5747))+p1mup3*(p2mup4*F(5734)+p1mup4
     -   *F(5738)+p3mup4*F(5739)-p4mup4*F(5751))+p3mup3*(-(p2mup4*F(5
     -   745))+p1mup4*F(5752)+p3mup4*F(5753)-p4mup4*F(5755)))
       Fa38(412) = p6mup4*(p1mup3*p2mup2*F(5759)+p2mup2*p2mup3*F(5760)
     -   +p2mup2*p3mup3*F(5761)+p2mup2*p4mup3*F(5762)+p1mup2*(p1mup3*
     -   F(5756)+p2mup3*F(5757)+p3mup3*F(5758)+p4mup3*F(5763))+p2mup3
     -   *p4mup2*F(5764)+p3mup3*p4mup2*F(5765)+p3mup2*p4mup3*F(5765)-
     -   p3mup2*p3mup3*F(5766)-p2mup2*p6mup3*F(5767)-p2mup3*p3mup2*F(
     -   5768)-p2mup3*p6mup2*F(5769)-p1mup3*p3mup2*F(5770)-p3mup3*p6m
     -   up2*F(5771)-p4mup2*p4mup3*F(5772)-p3mup2*p6mup3*F(5773))
       Fa38(413) = mup3mup4*(p1mup2*F(5775)+p2mup2*F(5776)+p6mup2*F(57
     -   78))+mup2mup3*p6mup4*F(5779)+mup2mup4*(p1mup3*F(5774)-p2mup3
     -   *F(5777)+p3mup3*F(5788)+p4mup3*F(5789)-p6mup3*F(5792))-2*(p6
     -   mup3*(p1mup2*(p1mup4*F(5780)+p2mup4*F(5781)+p3mup4*F(5782)+p
     -   4mup4*F(5783))+p2mup2*(p1mup4*F(5784)+p2mup4*F(5785)+p3mup4*
     -   F(5786)+p4mup4*F(5787)))-p1mup3*p4mup2*p6mup4*F(5793)-p1mup3
     -   *p3mup2*p3mup4*F(5798)-p2mup3*p3mup2*p3mup4*F(5799)+p1mup3*p
     -   6mup2*p6mup4*F(5801)+p4mup3*p6mup2*p6mup4*F(5802)-p1mup4*p3m
     -   up2*p3mup3*F(5803)-p2mup4*p3mup2*p3mup3*F(5804)+p3mup2*p3mup
     -   3*p3mup4*F(5805)-p1mup2*(-(p6mup3*p6mup4*F(5800))+p4mup3*(p2
     -   mup4*F(5790)-p1mup4*F(5794)+p3mup4*F(5795)+p4mup4*F(5806)))-
     -   p2mup2*p4mup3*(p2mup4*F(5791)-p1mup4*F(5796)+p3mup4*F(5797)+
     -   p4mup4*F(5807))+p4mup2*p6mup3*p6mup4*F(5808)-p6mup2*p6mup3*p
     -   6mup4*F(5809))
       Fa38(414) = -(mup3mup4*(p3mup2*F(5812)+p4mup2*F(5814)))+mup2mup
     -   3*(p1mup4*F(5810)+p2mup4*F(5811)+p3mup4*F(5813)+p4mup4*F(581
     -   5))+2*(-(p3mup2*p3mup4*p6mup3*F(5822))+p2mup3*p2mup4*p3mup2*
     -   F(5823)-p2mup3*p2mup4*p4mup2*F(5824)-p2mup3*p3mup4*p4mup2*F(
     -   5825)-p2mup4*p3mup2*p6mup3*F(5826)-p3mup2*p4mup4*p6mup3*F(58
     -   30)+p2mup3*p2mup4*p6mup2*F(5832)+p2mup3*p3mup4*p6mup2*F(5833
     -   )+p4mup2*p4mup4*p6mup3*F(5836)+p2mup4*p4mup2*p6mup3*F(5837)+
     -   p2mup4*p6mup2*p6mup3*F(5838)+p3mup4*p4mup2*p6mup3*F(5839)+p3
     -   mup4*p6mup2*p6mup3*F(5840)+p4mup4*p6mup2*p6mup3*F(5842)-p1mu
     -   p3*(p2mup4*(p3mup2*F(5816)+p4mup2*F(5818)-p6mup2*F(5820))-p1
     -   mup4*p3mup2*F(5827)+p3mup4*p4mup2*F(5828)+p1mup4*p4mup2*F(58
     -   31)-p3mup4*p6mup2*F(5834)-p1mup4*p6mup2*F(5835)+p3mup2*p4mup
     -   4*F(5844)+p4mup2*p4mup4*F(5847)-p4mup4*p6mup2*F(5849))+p1mup
     -   4*(p2mup3*p3mup2*F(5817)+p2mup3*p4mup2*F(5819)-p2mup3*p6mup2
     -   *F(5821)-p3mup2*p6mup3*F(5829)+p4mup2*p6mup3*F(5841)+p6mup2*
     -   p6mup3*F(5843)+p3mup3*p4mup2*F(5845)+p3mup2*p4mup3*F(5845)-p
     -   3mup3*p6mup2*F(5846)+p4mup2*p4mup3*F(5848)-p4mup3*p6mup2*F(5
     -   850))-p2mup3*p3mup2*p4mup4*F(5851)+p2mup4*p3mup3*p4mup2*F(58
     -   52)+p2mup4*p3mup2*p4mup3*F(5852)-p2mup4*p3mup3*p6mup2*F(5853
     -   )-p2mup3*p4mup2*p4mup4*F(5854)+p2mup4*p4mup2*p4mup3*F(5855)+
     -   p2mup3*p4mup4*p6mup2*F(5856)-p2mup4*p4mup3*p6mup2*F(5857)+p3
     -   mup3*p3mup4*p4mup2*F(5858)+p3mup2*p3mup4*p4mup3*F(5858)-p3mu
     -   p3*p3mup4*p6mup2*F(5859)+p3mup4*p4mup2*p4mup3*F(5860)-p3mup4
     -   *p4mup3*p6mup2*F(5861)-p4mup2*p4mup3*p4mup4*F(5862)+p3mup2*p
     -   3mup3*p4mup4*F(5863)+p3mup3*p4mup2*p4mup4*F(5864)+p3mup2*p4m
     -   up3*p4mup4*F(5864)-p3mup3*p4mup4*p6mup2*F(5865)-p4mup3*p4mup
     -   4*p6mup2*F(5866))
       Fa38(415) = -6*p2mup3*p2mup4*F(5572)-mup3mup4*F(5867)+2*(-3*(p4
     -   mup3*p4mup4*F(5414)-p2mup3*p4mup4*F(5580)+p2mup4*p4mup3*F(55
     -   81))+p1mup3*p2mup4*F(5868)-p1mup4*p2mup3*F(5869)+p1mup3*p1mu
     -   p4*F(5870)-p1mup3*p6mup4*F(5872)-p2mup3*p6mup4*F(5873)-p4mup
     -   3*p6mup4*F(5875)+p1mup4*p6mup3*F(5876)+p2mup4*p6mup3*F(5877)
     -   +p3mup4*(-3*(p4mup3*F(5392)-p2mup3*F(5579))+p1mup3*F(5871)+p
     -   6mup3*F(5878))+p4mup4*p6mup3*F(5879)-p1mup3*p4mup4*F(5880)+p
     -   3mup3*(3*(p3mup4*F(5388)-p4mup4*F(5412)-p2mup4*F(5577))-p6mu
     -   p4*F(5874)+p1mup4*F(5881))+p1mup4*p4mup3*F(5882)+p6mup3*p6mu
     -   p4*F(5883))
       Fa38(416) = p1mup3*F(5884)-p2mup3*F(5885)+p6mup3*F(5886)-p3mup3
     -   *F(5887)-p4mup3*F(5888)
       Fa38(417) = -(mup2mup3*F(5889))+2*(-3*(p3mup3*p4mup2*F(5141)+p4
     -   mup2*p4mup3*F(5143))+p1mup2*p1mup3*F(5890)-p1mup2*p6mup3*F(5
     -   893)+p3mup3*p6mup2*F(5895)+p2mup3*(-3*(p3mup2*F(5111)+p4mup2
     -   *F(5560))+p1mup2*F(5891)+p6mup2*F(5896))+p1mup3*p6mup2*F(589
     -   7)-p4mup2*p6mup3*F(5898)+p4mup3*p6mup2*F(5899)+p2mup2*(3*(p2
     -   mup3*F(5108)-p3mup3*F(5137)-p4mup3*F(5138))-p6mup3*F(5892)+p
     -   1mup3*F(5900))-p1mup2*p3mup3*F(5901)-p3mup2*(3*(p3mup3*F(514
     -   0)+p4mup3*F(5141))+p6mup3*F(5894)-p1mup3*F(5902))-p1mup2*p4m
     -   up3*F(5903)+p1mup3*p4mup2*F(5904)+p6mup2*p6mup3*F(5905))
       Fa38(418) = p1mup4*F(5906)-p2mup4*F(5907)+p6mup4*F(5908)+p3mup4
     -   *F(5909)-p4mup4*F(5910)
       Return
       End
