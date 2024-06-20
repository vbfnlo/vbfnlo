c************** Calling the Fa functions*************************
       subroutine HexBoxCross1Fa11(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,
     &   p3mup4,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa
     &   11)
       IMPLICIT NONE
      Complex*16   p1mup3, p1mup4, p2mup3, 
     -          p2mup4, p3mup3, p3mup4, 
     -          p4mup3, p4mup4, p5mup3, 
     -          p5mup4, p6mup3, p6mup4
       Complex*16   mup3mup4

        common/HexBoxCross1FhlFunctions/F
       COMMON/HexBoxCross1Invariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12
     &   ,s23,s34,s45,s56,s16,s123,s234,s345
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
      Complex*16 Fa(272),F(1461)
      Real*8 P(1825) 
        COMMON/HexBoxCross1PFunctions/P
       Complex*16 Fa11(131:143)
       COMMON/HexBoxCross1FaFunctions/Fa
       Fa11(131) = 2*(p1mup3*F(293)+p2mup3*F(295)+p6mup3*F(296))+(p4mu
     &   p3+p5mup3)*F(464)
       Fa11(132) = 2*(p4mup4+p5mup4)*(-(p1mup3*F(465))-p2mup3*F(466)+(
     &   p4mup3+p5mup3)*F(467)+p6mup3*F(1458))+mup3mup4*F(1455)*P(52)
       Fa11(133) = -2*mup3mup4*s45*F(274)+4*(p4mup4+p5mup4)*(p2mup3*F(
     &   360)+p1mup3*F(362)+(p4mup3+p5mup3)*F(364)+p6mup3*F(365))
       Fa11(134) = 2*(p4mup4+p5mup4)*(-(p1mup3*F(334))-p2mup3*F(336)+(
     &   p4mup3+p5mup3)*F(338)+p6mup3*F(339))+mup3mup4*s45*F(468)
       Fa11(135) = 2*(p4mup4+p5mup4)*(-(p1mup3*F(334))-p2mup3*F(336)+(
     &   p4mup3+p5mup3)*F(338)+p6mup3*F(339))+mup3mup4*s45*F(367)
       Fa11(136) = 4*(p1mup4*(p2mup3*F(360)+p1mup3*F(362)+(p4mup3+p5mu
     &   p3)*F(364)+p6mup3*F(365))-p6mup4*(p1mup3*F(376)-(p4mup3+p5mu
     &   p3)*F(377)+p6mup3*F(378)+p2mup3*F(393)))-mup3mup4*F(469)
       Fa11(137) = mup3mup4*F(313)+2*(p1mup3*(p6mup4*F(333)-p1mup4*F(3
     &   34))+p1mup4*(-(p2mup3*F(336))+(p4mup3+p5mup3)*F(338)+p6mup3*
     &   F(339))+p6mup4*(p6mup3*F(337)+p2mup3*F(340)-(p4mup3+p5mup3)*
     &   F(392)))
       Fa11(138) = -2*mup3mup4*F(370)-(p4mup4+p5mup4)*(4*(p1mup3*F(376
     &   )-(p4mup3+p5mup3)*F(377)+p6mup3*F(378))-p2mup3*F(470))
       Fa11(139) = -2*mup3mup4*F(458)+(p4mup4+p5mup4)*(4*(p1mup3*F(333
     &   )+p6mup3*F(337)+p2mup3*F(340))+(p4mup3+p5mup3)*F(471))
       Fa11(140) = -2*mup3mup4*F(384)+(p4mup4+p5mup4)*(4*(p1mup3*F(333
     &   )+p6mup3*F(337)+p2mup3*F(340))+(p4mup3+p5mup3)*F(471))
       Fa11(141) = -(p1mup3*F(341))+(p4mup3+p5mup3)*F(342)-p2mup3*F(34
     &   3)+p6mup3*F(478)
       Fa11(142) = p1mup3*F(333)+p6mup3*F(337)+p2mup3*F(340)-(p4mup3+p
     &   5mup3)*F(392)
       Fa11(143) = 2*p6mup4*F(275)+(p4mup4+p5mup4)*F(431)
       Return
       End
