c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineNoAbe1g38(p1mup2,p1mup3,p1mup4
     -   ,p1mup6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3m
     -   up6,p4mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,
     -   p6mup2,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3m
     -   up4,mup3mup6,mup4mup6,Fa38)
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
       Complex*16 Fa38(519:532)
       COMMON/test_FaFunctionsNoAbe1g/Fa
       Fa38(519) = p5mup6*(p1mup2*(p1mup3*F(4663)-p3mup3*F(4674)+p2mup
     -   3*F(5654))+p2mup2*(p1mup3*F(5653)+p2mup3*F(5655)-p3mup3*F(56
     -   56)))
       Fa38(520) = p5mup6*(p4mup3*(p1mup2*F(5658)+p2mup2*F(5659)+p3mup
     -   2*F(5661))+p4mup2*(-(p2mup3*F(4675))-p1mup3*F(5657)+p3mup3*F
     -   (5660)+p4mup3*F(5662)+p5mup3*F(5663)))
       Fa38(521) = 4*p5mup6*(mup2mup3*F(5664)+2*(p1mup2*p5mup3*F(4682)
     -   +p3mup3*p5mup2*F(4685)+p5mup2*p5mup3*F(4686)-p2mup3*p3mup2*F
     -   (5665)-p1mup3*(p5mup2*F(4681)+p3mup2*F(5666))+p2mup2*p5mup3*
     -   F(5667)-p2mup3*p5mup2*F(5668)+p3mup2*p5mup3*F(5669)+p3mup2*p
     -   3mup3*F(5670)+p4mup3*p5mup2*F(5671)))
       Fa38(522) = mup3mup4*(p1mup6*F(5677)+4*(p5mup6*F(5673)+p3mup6*F
     -   (5675)+p4mup6*F(5676)+p2mup6*F(5678)))+8*(mup3mup6*(p1mup4*F
     -   (5674)-p3mup4*F(5681)+p5mup4*F(5682)+p4mup4*F(5683))-p2mup4*
     -   (p3mup6*p5mup3*F(1724)+p1mup6*p5mup3*F(1726)+p5mup3*p5mup6*F
     -   (1727)+p3mup6*p4mup3*F(5355)+p4mup3*p4mup6*F(5356)+mup3mup6*
     -   F(5680)+p2mup3*p5mup6*F(5684)+p2mup3*p3mup6*F(5685)+p2mup3*p
     -   4mup6*F(5687)+p3mup3*(p3mup6*F(1717)+p1mup6*F(1721)+p5mup6*F
     -   (1723)+p2mup6*F(5686)+p4mup6*F(5688))+p2mup3*p2mup6*F(5689)+
     -   p1mup6*p2mup3*F(5691)+p1mup3*(p1mup6*F(1719)+p3mup6*F(1720)+
     -   p5mup6*F(1722)+p2mup6*F(5690)+p4mup6*F(5692))+p4mup3*p5mup6*
     -   F(5693)+p2mup6*p4mup3*F(5698)+p1mup6*p4mup3*F(5703)+p2mup6*p
     -   5mup3*F(5707)+p4mup6*p5mup3*F(5708))+p2mup3*(mup4mup6*F(5679
     -   )+p1mup4*p5mup6*F(5695)+p1mup4*p3mup6*F(5696)+p1mup4*p4mup6*
     -   F(5697)+p1mup4*p2mup6*F(5699)+p1mup4*p1mup6*F(5702)+p3mup4*(
     -   p3mup6*F(5355)+p4mup6*F(5356)+p5mup6*F(5693)+p2mup6*F(5700)+
     -   p1mup6*F(5703))+p4mup4*(p4mup6*F(5005)+p3mup6*F(5358)+p5mup6
     -   *F(5694)+p2mup6*F(5701)+p1mup6*F(5704))+p5mup4*p5mup6*F(5705
     -   )+p4mup6*p5mup4*F(5706)+p3mup6*p5mup4*F(5709)+p2mup6*p5mup4*
     -   F(5710)+p1mup6*p5mup4*F(5711))-mup4mup6*(p1mup3*F(5672)-p3mu
     -   p3*F(5712)-p4mup3*F(5713)-p5mup3*F(5714))+2*p2mup6*(-(p3mup3
     -   *p3mup4*F(5716))+p3mup4*p4mup3*F(5717)+p1mup3*p4mup4*F(5718)
     -   +p3mup4*p5mup3*F(5719)+p1mup3*p5mup4*F(5720)+p1mup3*p3mup4*F
     -   (5722)-p1mup4*(p3mup3*F(5715)+p5mup3*F(5721)+p4mup3*F(5723)-
     -   p1mup3*F(5724))+p3mup3*p4mup4*F(5725)+p3mup3*p5mup4*F(5726)+
     -   p4mup3*p4mup4*F(5727)+p4mup4*p5mup3*F(5728)+p4mup3*p5mup4*F(
     -   5728)+p5mup3*p5mup4*F(5729)))
       Fa38(523) = p5mup3*F(5730)+p2mup3*F(5731)+p4mup3*F(5732)+p3mup3
     -   *F(5733)-2*p1mup3*F(5734)
       Fa38(524) = -16*p2mup2*p2mup6*p4mup3*F(4993)+mup2mup6*(p4mup3*F
     -   (5735)+p2mup3*F(5739)+p1mup3*F(5740))+mup3mup6*(p4mup2*F(573
     -   6)+4*(p2mup2*F(5737)+p1mup2*F(5738)+p3mup2*F(5743)+p5mup2*F(
     -   5748)))+mup2mup6*p5mup3*F(5749)+p3mup3*(mup2mup6*F(5741)+16*
     -   (p3mup6*p4mup2*F(1717)+p4mup2*(p4mup6*F(5767)+p2mup6*F(5770)
     -   +p1mup6*F(5772)+p5mup6*F(5774))+p4mup6*(p1mup2*F(5789)+p2mup
     -   2*F(5801)+p3mup2*F(5803)+p5mup2*F(5804))))-4*(mup2mup3*p2mup
     -   6*F(5742)+mup2mup3*(p3mup6*F(5744)+p4mup6*F(5745)-p1mup6*F(5
     -   746)+p5mup6*F(5747))+2*(-(p4mup2*p4mup6*p5mup3*F(5776))+p4mu
     -   p3*p4mup6*p5mup2*F(5777)-p2mup6*p4mup2*p5mup3*F(5778)+p2mup6
     -   *p4mup3*p5mup2*F(5779)-p3mup6*p4mup2*p5mup3*F(5780)+p3mup6*p
     -   4mup3*p5mup2*F(5781)-p1mup6*p4mup2*p5mup3*F(5782)+p1mup6*p4m
     -   up3*p5mup2*F(5783)-p4mup2*p5mup3*p5mup6*F(5784)+p4mup3*p5mup
     -   2*p5mup6*F(5785)+p1mup3*(p1mup6*p4mup2*F(5752)+p2mup6*p4mup2
     -   *F(5753)+p3mup6*p4mup2*F(5755)+p4mup2*p4mup6*F(5757)+p4mup2*
     -   p5mup6*F(5758)+2*(p2mup2*p4mup6*F(5788)+p4mup6*p5mup2*F(5796
     -   )+p3mup2*p4mup6*F(5798)))-p1mup2*(p1mup6*p4mup3*F(5750)+p2mu
     -   p6*p4mup3*F(5751)+p3mup6*p4mup3*F(5754)+p4mup3*p4mup6*F(5756
     -   )+p4mup3*p5mup6*F(5759)+2*(p4mup6*p5mup3*F(5795)+p2mup3*p4mu
     -   p6*F(5797)+p1mup3*p4mup6*F(5799)))+4*(p4mup2*p4mup3*p5mup6*F
     -   (5793)-p1mup6*p4mup2*p4mup3*F(5800))+2*(p2mup2*p4mup3*p4mup6
     -   *F(5760)+p2mup2*p3mup6*p4mup3*F(5761)+p1mup6*p2mup2*p4mup3*F
     -   (5762)-p1mup6*p2mup3*p4mup2*F(5763)+p2mup2*p4mup3*p5mup6*F(5
     -   764)-p2mup3*p4mup2*p4mup6*F(5765)+p3mup2*p4mup3*p4mup6*F(576
     -   6)-p2mup3*p3mup6*p4mup2*F(5768)+p3mup2*p3mup6*p4mup3*F(5768)
     -   -p2mup3*p2mup6*p4mup2*F(5769)+p2mup6*p3mup2*p4mup3*F(5769)+p
     -   1mup6*p3mup2*p4mup3*F(5771)-p2mup3*p4mup2*p5mup6*F(5773)+p3m
     -   up2*p4mup3*p5mup6*F(5775)-p2mup3*p3mup2*p4mup6*F(5786)+p2mup
     -   2*p2mup3*p4mup6*F(5787)+p2mup6*p4mup2*p4mup3*F(5790)+p3mup6*
     -   p4mup2*p4mup3*F(5791)+p4mup2*p4mup3*p4mup6*F(5792)-p2mup3*p4
     -   mup6*p5mup2*F(5794)-p2mup2*p4mup6*p5mup3*F(5802)-p3mup2*p4mu
     -   p6*p5mup3*F(5804)-p4mup6*p5mup2*p5mup3*F(5805))))
       Fa38(525) = mup3mup6*F(5806)+4*(p1mup3*p1mup6*F(5807)+p1mup3*p4
     -   mup6*F(5809)+p1mup3*p5mup6*F(5811)-p1mup6*p2mup3*F(5812)-3*(
     -   p4mup3*p4mup6*F(5356)+p2mup3*p4mup6*F(5813))-p2mup3*p5mup6*F
     -   (5815)-p3mup3*(3*(p3mup6*F(1717)+p4mup6*F(5688)+p2mup6*F(577
     -   0))+p1mup6*F(5814)+p5mup6*F(5816))-p1mup6*p4mup3*F(5817)-p4m
     -   up3*p5mup6*F(5818)-p3mup6*(3*(p4mup3*F(5355)+p2mup3*F(5768))
     -   -p1mup3*F(5810)+p5mup3*F(5819))-p4mup6*p5mup3*F(5820)-p1mup6
     -   *p5mup3*F(5821)-p2mup6*(3*(p4mup3*F(5354)+p2mup3*F(5769))-p1
     -   mup3*F(5808)+p5mup3*F(5822))-p5mup3*p5mup6*F(5823))
       Fa38(526) = p3mup6*F(5824)+p4mup6*F(5826)+2*(p2mup6*F(5825)-p5m
     -   up6*F(5827))+p1mup6*F(5828)
       Fa38(527) = p5mup4*F(5829)+p4mup4*F(5830)+p2mup4*F(5831)+p1mup4
     -   *F(5832)+p3mup4*F(5833)
       Fa38(528) = p2mup6*F(5834)+p4mup6*F(5835)+p5mup6*F(5836)+p1mup6
     -   *F(5837)+p3mup6*F(5838)
       Fa38(529) = mup2mup6*F(5839)+4*(-(p1mup2*p4mup6*F(5840))-p1mup2
     -   *p5mup6*F(5841)-p1mup2*p1mup6*F(5842)-p1mup2*p3mup6*F(5844)+
     -   p2mup2*(3*(p2mup6*F(4993)+p4mup6*F(4995)+p3mup6*F(5761))+p5m
     -   up6*F(5845)+p1mup6*F(5846))+p3mup2*(3*(p3mup6*F(5768)+p4mup6
     -   *F(5813))+p5mup6*F(5847)+p1mup6*F(5848))+p4mup2*(3*(p4mup6*F
     -   (5004)+p3mup6*F(5359))+p1mup6*F(5849)+p5mup6*F(5850))+p4mup6
     -   *p5mup2*F(5851)+p2mup6*(3*(p4mup2*F(5002)+p3mup2*F(5769))-p1
     -   mup2*F(5843)+p5mup2*F(5852))+p3mup6*p5mup2*F(5853)+p1mup6*p5
     -   mup2*F(5854)+p5mup2*p5mup6*F(5855))
       Fa38(530) = p5mup4*(p2mup3*p4mup2*F(5857)+p3mup3*p4mup2*F(5859)
     -   -p4mup2*p4mup3*F(5860)-p4mup3*p5mup2*F(5862)+p2mup2*(-(p4mup
     -   3*F(5856))+p5mup3*F(5861)+p3mup3*F(5863)+p2mup3*F(5864))-p2m
     -   up3*p5mup2*F(5865)-p3mup3*p5mup2*F(5867)+p4mup2*p5mup3*F(586
     -   8)+p3mup2*(-(p4mup3*F(5858))+p5mup3*F(5866)+p3mup3*F(5869)+2
     -   *p2mup3*F(5870))+p5mup2*p5mup3*F(5871)+mup2mup3*F(5872))
       Fa38(531) = p2mup2*p4mup4*(p2mup3*F(5878)+p3mup3*F(5879))+p1mup
     -   2*(p4mup4*(p1mup3*F(5874)+p2mup3*F(5875)+p3mup3*F(5876))+p5m
     -   up4*(-(p4mup3*F(5873))+p1mup3*F(5881)+p2mup3*F(5882)+p3mup3*
     -   F(5884)+p5mup3*F(5885)))+p1mup3*(p2mup2*p4mup4*F(5877)+p4mup
     -   2*p5mup4*F(5880)+p2mup2*p5mup4*F(5883)-p5mup2*p5mup4*F(5886)
     -   +p3mup2*p5mup4*F(5887))
       Fa38(532) = p1mup2*(p3mup3*p3mup4*F(5888)+p4mup4*p5mup3*F(5892)
     -   +p1mup4*p5mup3*F(5894)+p2mup4*p5mup3*F(5895)+p3mup4*(p5mup3*
     -   F(5890)+p1mup3*F(5898)+p2mup3*F(5899)-p4mup3*F(5902))+p4mup3
     -   *p4mup4*F(5904)-mup3mup4*F(5906))+p2mup2*(p3mup3*p3mup4*F(58
     -   89)+p4mup4*p5mup3*F(5893)+p1mup4*p5mup3*F(5896)+p2mup4*p5mup
     -   3*F(5897)+p3mup4*(p5mup3*F(5891)+p1mup3*F(5900)+p2mup3*F(590
     -   1)-p4mup3*F(5903))+p4mup3*p4mup4*F(5905)-mup3mup4*F(5907))
       Return
       End
