c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineNoAbe1g39(p1mup2,p1mup3,p1mup4
     -   ,p1mup6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3m
     -   up6,p4mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,
     -   p6mup2,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3m
     -   up4,mup3mup6,mup4mup6,Fa39)
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
       Complex*16 Fa39(533:546)
       COMMON/test_FaFunctionsNoAbe1g/Fa
       Fa39(533) = p3mup2*p3mup4*p5mup3*F(5914)+p2mup4*p3mup2*p5mup3*F
     -   (5915)+p3mup2*p4mup4*p5mup3*F(5916)-p1mup3*p1mup4*p5mup2*F(5
     -   917)-p1mup3*p2mup4*p5mup2*F(5918)-p1mup4*p2mup3*p5mup2*F(591
     -   9)-p1mup3*p3mup4*p5mup2*F(5920)-p1mup4*p3mup3*p5mup2*F(5921)
     -   -p1mup3*p4mup4*p5mup2*F(5922)-p2mup3*p2mup4*p5mup2*F(5923)-p
     -   2mup3*p3mup4*p5mup2*F(5924)-p2mup4*p3mup3*p5mup2*F(5925)-p2m
     -   up3*p4mup4*p5mup2*F(5926)-p3mup3*p3mup4*p5mup2*F(5927)-p3mup
     -   3*p4mup4*p5mup2*F(5928)-p2mup4*p3mup2*p4mup3*F(5930)+p2mup4*
     -   p3mup3*p4mup2*F(5931)+p3mup2*p3mup3*p3mup4*F(5932)+p2mup3*p2
     -   mup4*p4mup2*F(5934)+2*(p2mup3*p3mup2*p4mup4*F(5929)+p2mup3*p
     -   2mup4*p3mup2*F(5933)+p2mup3*p3mup2*p3mup4*F(5935))+p2mup4*p3
     -   mup2*p3mup3*F(5936)+p3mup2*p3mup3*p4mup4*F(5937)+p1mup3*p1mu
     -   p4*p3mup2*F(5938)+p2mup3*p4mup2*p4mup4*F(5939)+p3mup3*p4mup2
     -   *p4mup4*F(5940)+p4mup2*p4mup4*p5mup3*F(5941)+p2mup3*p3mup4*p
     -   4mup2*F(5942)+p2mup4*p4mup2*p5mup3*F(5943)+p3mup4*p4mup2*p5m
     -   up3*F(5944)+p2mup2*(p2mup3*p2mup4*F(5908)+p2mup4*p3mup3*F(59
     -   09)+p1mup3*p1mup4*F(5912)+p2mup4*p4mup3*F(5945))+p3mup3*p3mu
     -   p4*p4mup2*F(5946)+p1mup2*(p1mup3*p2mup4*F(5910)+p1mup4*p2mup
     -   3*F(5911)+p1mup4*p3mup3*F(5913)+p1mup4*p4mup3*F(5947))+p1mup
     -   3*p1mup4*p4mup2*F(5948)+p1mup4*p5mup2*p5mup3*F(5949)-p2mup4*
     -   p4mup3*p5mup2*F(5950)+p2mup4*p5mup2*p5mup3*F(5951)+p3mup4*p5
     -   mup2*p5mup3*F(5952)+p4mup4*p5mup2*p5mup3*F(5953)+p3mup2*p3mu
     -   p4*p4mup3*F(5954)-p1mup4*p4mup3*p5mup2*F(5955)-p3mup4*p4mup3
     -   *p5mup2*F(5956)-p4mup3*p4mup4*p5mup2*F(5957)-mup3mup4*p3mup2
     -   *F(5958)+mup3mup4*p4mup2*F(5959)+mup2mup4*p2mup3*F(5960)+mup
     -   2mup3*p2mup4*F(5961)+mup2mup4*p3mup3*F(5962)+mup2mup3*p3mup4
     -   *F(5963)-mup2mup4*p4mup3*F(5964)+mup2mup3*p4mup4*F(5965)-mup
     -   3mup4*p5mup2*F(5966)+mup2mup4*p5mup3*F(5967)+mup2mup4*p1mup3
     -   *F(5968)+mup2mup3*p1mup4*F(5969)-p2mup4*p4mup2*p4mup3*F(5970
     -   )+p3mup2*p4mup3*p4mup4*F(5971)-p3mup4*p4mup2*p4mup3*F(5972)+
     -   p4mup2*p4mup3*p4mup4*F(5973)
       Fa39(534) = p1mup2*(p2mup3*p2mup4*F(5974)+p2mup4*p3mup3*F(5977)
     -   +p1mup3*p1mup4*F(5979)+p2mup4*p4mup3*F(5987))+p1mup3*(p2mup2
     -   *p2mup4*F(5975)+p2mup4*p3mup2*F(5981)+p2mup4*p4mup2*F(5983)+
     -   p3mup2*p3mup4*F(5984)+p3mup2*p4mup4*F(5986)+p3mup4*p4mup2*F(
     -   5990)+p4mup2*p4mup4*F(5992))+p1mup4*(p2mup2*(p2mup3*F(5976)+
     -   p3mup3*F(5978)+p4mup3*F(5988))+p3mup2*(p5mup3*F(5980)+2*p2mu
     -   p3*F(5982)+p3mup3*F(5985)+p4mup3*F(5994))+p4mup2*(p2mup3*F(5
     -   989)+p3mup3*F(5991)+p5mup3*F(5993)-p4mup3*F(5995)))
       Fa39(535) = mup4mup6*(p2mup2*p2mup3*F(5787)+p3mup3*p4mup2*F(599
     -   6)+p3mup2*p4mup3*F(5996)+p3mup3*p5mup2*F(5997)+p3mup2*p5mup3
     -   *F(5997)+p4mup2*p4mup3*F(5998)+p4mup2*p5mup3*F(5999)+p5mup2*
     -   p5mup3*F(6000))
       Fa39(536) = mup4mup6*(p1mup2*p3mup3*F(5789)+p1mup2*p1mup3*F(579
     -   9)-p4mup3*p5mup2*F(5999)+p2mup3*(p3mup2*F(5786)+p5mup2*F(579
     -   4)+p1mup2*F(5797)+p4mup2*F(6001))+mup2mup3*F(6002)-p2mup2*(p
     -   1mup3*F(5788)+p3mup3*F(6003))-p2mup2*p4mup3*F(6004)-p2mup2*p
     -   5mup3*F(6005)-p3mup2*p3mup3*F(6006))
       Fa39(537) = mup4mup6*(-(p1mup2*(p5mup3*F(5795)+p4mup3*F(6007)))
     -   +p1mup3*(p5mup2*F(5796)+p3mup2*F(5798)+p4mup2*F(6008)))
       Fa39(538) = mup3mup6*(p4mup2*p4mup4*F(5036)+p4mup2*p5mup4*F(503
     -   7)+p1mup2*(p4mup4*F(5024)+p5mup4*F(5028)+p3mup4*F(6009))+p2m
     -   up2*(p4mup4*F(5034)+p5mup4*F(5035)+p3mup4*F(6010))-p3mup2*p3
     -   mup4*F(6011)+p3mup4*p4mup2*F(6012)+p3mup2*p4mup4*F(6013)+p3m
     -   up2*p5mup4*F(6014))
       Fa39(539) = mup3mup6*(p5mup2*p5mup4*F(5038)+p2mup4*(-(p2mup2*F(
     -   5026))+p4mup2*F(5027)+p1mup2*F(5031)+p3mup2*F(6015))-p1mup4*
     -   (p2mup2*F(5025)+p4mup2*F(5032)-p1mup2*F(5033)+p3mup2*F(6016)
     -   ))+mup2mup6*p5mup4*(p4mup3*F(6017)+p5mup3*F(6018))
       Fa39(540) = mup3mup4*(p2mup6*p4mup2*F(5002)+p4mup2*p4mup6*F(500
     -   4)+p4mup2*p5mup6*F(5009)+p1mup6*p4mup2*F(5013)+p2mup6*p5mup2
     -   *F(5019)+p4mup6*p5mup2*F(5021)+p5mup2*p5mup6*F(5022)+p1mup6*
     -   p5mup2*F(5023)+p3mup6*p4mup2*F(5359)+p2mup2*(p5mup6*F(4987)+
     -   p1mup6*F(4992)+p2mup6*F(4993)+p4mup6*F(4995)+p3mup6*F(5761))
     -   +p3mup2*p3mup6*F(5768)+p2mup6*p3mup2*F(5769)+p3mup2*p4mup6*F
     -   (5813)+p1mup2*(p5mup6*F(4986)+p1mup6*F(4988)+p2mup6*F(4989)+
     -   p4mup6*F(4991)+p3mup6*F(6019))+p3mup2*p5mup6*F(6020)+p1mup6*
     -   p3mup2*F(6021)+p3mup6*p5mup2*F(6022)+2*mup2mup6*F(6023))
       Fa39(541) = mup2mup6*(-(p1mup3*p4mup4*F(5718))-p1mup3*p3mup4*F(
     -   5722)-p1mup3*p2mup4*F(6024)-p2mup3*p2mup4*F(6025)+p1mup4*(p3
     -   mup3*F(5715)-p1mup3*F(5724)+p2mup3*F(6027))+p2mup3*p3mup4*F(
     -   6028)+p2mup3*p4mup4*F(6029)+p3mup3*(p3mup4*F(5716)-p2mup4*F(
     -   6026)+p4mup4*F(6030)))
       Fa39(542) = mup2mup3*(p4mup4*p4mup6*F(5005)+p3mup4*p3mup6*F(535
     -   5)+p3mup4*p4mup6*F(5356)+p3mup6*p4mup4*F(5358)+p2mup4*(p2mup
     -   6*F(5002)+p4mup6*F(5004)+p5mup6*F(5009)+p1mup6*F(5013)+p3mup
     -   6*F(5359))+p3mup4*p5mup6*F(5693)+p4mup4*p5mup6*F(5694)+p1mup
     -   4*p5mup6*F(5695)+p1mup4*p3mup6*F(5696)+p1mup4*p4mup6*F(5697)
     -   +p1mup4*p1mup6*F(5702)+p1mup6*p3mup4*F(5703)+p1mup6*p4mup4*F
     -   (5704)+p5mup4*p5mup6*F(5705)+p4mup6*p5mup4*F(5706)+p3mup6*p5
     -   mup4*F(5709)+p1mup6*p5mup4*F(5711)+p2mup6*(p4mup4*F(5003)+p3
     -   mup4*F(5354)+p1mup4*F(6031)+p5mup4*F(6032)))
       Fa39(543) = mup2mup4*(p3mup6*p5mup3*F(1724)+p1mup6*p5mup3*F(172
     -   6)+p5mup3*p5mup6*F(1727)+p2mup6*p4mup3*F(5354)+p3mup6*p4mup3
     -   *F(5355)+p4mup3*p4mup6*F(5356)+p4mup3*p5mup6*F(5693)+p1mup6*
     -   p4mup3*F(5703)+p4mup6*p5mup3*F(5708)+p2mup3*p3mup6*F(5768)+p
     -   2mup3*p2mup6*F(5769)+p3mup3*(p3mup6*F(1717)+p1mup6*F(1721)+p
     -   5mup6*F(1723)+p4mup6*F(5688)+p2mup6*F(5770))+p2mup3*p4mup6*F
     -   (5813)+p2mup3*p5mup6*F(6020)+p1mup6*p2mup3*F(6021)+2*mup3mup
     -   6*F(6023)+p1mup3*(p1mup6*F(1719)+p3mup6*F(1720)+p5mup6*F(172
     -   2)+p4mup6*F(5692)+p2mup6*F(6033))+p2mup6*p5mup3*F(6034))+2*(
     -   mup3mup6*p5mup2*(p2mup4*F(5029)-p1mup4*F(5030)+p4mup4*F(5037
     -   )+p3mup4*F(6038))+mup2mup6*(-(p3mup4*(p4mup3*F(5717)+p5mup3*
     -   F(5719)))-p1mup3*p5mup4*F(5720)+p1mup4*p5mup3*F(5721)+p1mup4
     -   *p4mup3*F(5723)+p4mup4*p5mup3*F(6017)-p2mup4*p4mup3*F(6035)+
     -   p2mup3*p5mup4*F(6036)-p2mup4*p5mup3*F(6037)+p3mup3*p5mup4*F(
     -   6039)+p4mup3*p4mup4*F(6040)))
       Fa39(544) = p1mup2*p1mup3*F(6045)+p2mup2*(p3mup3*F(6042)+p4mup3
     -   *F(6043)+p5mup3*F(6044)+p1mup3*F(6047))+p1mup2*p3mup3*F(6048
     -   )+p1mup2*p4mup3*F(6049)+p1mup2*p5mup3*F(6050)+p3mup2*p4mup3*
     -   F(6051)-p3mup3*p4mup2*F(6052)+p3mup2*p5mup3*F(6053)-p3mup3*p
     -   5mup2*F(6054)+p2mup3*(2*p3mup2*F(5786)+p2mup2*F(6041)+p1mup2
     -   *F(6046)-p4mup2*F(6055)-p5mup2*F(6056))-p1mup3*p4mup2*F(6057
     -   )-p1mup3*p5mup2*F(6058)+p4mup2*p5mup3*F(6059)-p4mup3*p5mup2*
     -   F(6060)+p3mup2*p3mup3*F(6061)+p1mup3*p3mup2*F(6062)+p4mup2*p
     -   4mup3*F(6063)+p5mup2*p5mup3*F(6064)+mup2mup3*F(6065)
       Fa39(545) = p5mup4*(-(p3mup2*p4mup3*F(6069))+p3mup3*p4mup2*F(60
     -   70)+p3mup2*p3mup3*F(6071)+p2mup2*(p3mup3*F(6066)+p2mup3*F(60
     -   67)+p4mup3*F(6073)+p5mup3*F(6074))+p2mup3*(p4mup2*F(6068)+p3
     -   mup2*F(6072)+p5mup2*F(6075))+p3mup3*p5mup2*F(6076)-p4mup2*p4
     -   mup3*F(6077)+p3mup2*p5mup3*F(6078)+p4mup2*p5mup3*F(6079)+mup
     -   2mup3*F(6080)+p4mup3*p5mup2*F(6081)+p5mup2*p5mup3*F(6082))
       Fa39(546) = p2mup2*p4mup4*(p2mup3*F(6091)+p3mup3*F(6092))+p1mup
     -   2*(p1mup3*p5mup4*F(6083)+p2mup3*p5mup4*F(6084)+p3mup3*p5mup4
     -   *F(6086)+p1mup3*p4mup4*F(6087)+p2mup3*p4mup4*F(6088)+p3mup3*
     -   p4mup4*F(6089)+p4mup3*p5mup4*F(6095)+p5mup3*p5mup4*F(6096))+
     -   p1mup3*(p2mup2*p5mup4*F(6085)+p2mup2*p4mup4*F(6090)+p4mup2*p
     -   5mup4*F(6093)+p3mup2*p5mup4*F(6094)+p5mup2*p5mup4*F(6097))
       Return
       End
