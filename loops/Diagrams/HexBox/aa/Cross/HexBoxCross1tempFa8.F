c************** Calling the Fa functions*************************
       subroutine HexBoxCross1Fa8(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p
     &   3mup4,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa8
     &   )
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
       Complex*16 Fa8(92:104)
       COMMON/HexBoxCross1FaFunctions/Fa
       Fa8(92) = -(p2mup3*F(358))+(p4mup3+p5mup3)*F(359)-p1mup3*F(361)
     &   +p6mup3*F(363)
       Fa8(93) = mup3mup4*F(389)-8*(-(p4mup4*p5mup3*F(338))-p5mup3*p5m
     &   up4*F(338)-p4mup4*p6mup3*F(339)-p5mup4*p6mup3*F(339)-p5mup3*
     &   p6mup4*F(359)+p1mup3*((p4mup4+p5mup4)*F(334)+p6mup4*F(361)-p
     &   2mup4*F(362))-p6mup3*p6mup4*F(363)-p2mup4*p5mup3*F(364)-p4mu
     &   p3*((p4mup4+p5mup4)*F(338)+p6mup4*F(359)+p2mup4*F(364))-p2mu
     &   p4*p6mup3*F(365)+p2mup3*((p4mup4+p5mup4)*F(336)-p2mup4*F(360
     &   )+p6mup4*F(390)))
       Fa8(94) = mup3mup4*F(391)-2*(p4mup4*(p6mup3*F(337)+p2mup3*F(340
     &   ))+p5mup4*(p6mup3*F(337)+p2mup3*F(340))-p1mup4*p6mup3*F(363)
     &   +p1mup3*((p4mup4+p5mup4)*F(333)+p1mup4*F(361)-p2mup4*F(376))
     &   -(p4mup3+p5mup3)*(p1mup4*F(359)-p2mup4*F(377))-p2mup4*p6mup3
     &   *F(378)+p1mup4*p2mup3*F(390)-p4mup4*(p4mup3+p5mup3)*F(392)-(
     &   p4mup3+p5mup3)*p5mup4*F(392)-p2mup3*p2mup4*F(393))
       Fa8(95) = p2mup3*F(400)+p4mup3*F(401)+p5mup3*F(402)+p1mup3*F(40
     &   3)+p6mup3*s45*F(404)
       Fa8(96) = 4*(p2mup3*F(263)-p6mup3*F(268))+p1mup3*F(405)+(p4mup3
     &   +p5mup3)*F(406)
       Fa8(97) = (p4mup4+p5mup4)*F(398)-p1mup4*F(408)-p6mup4*F(409)
       Fa8(98) = p1mup3*F(410)+p2mup3*F(411)+(p4mup3+p5mup3)*F(412)+4*
     &   p6mup3*F(413)
       Fa8(99) = mup3mup4*F(415)-16*(p4mup4+p5mup4)*(p1mup3*F(416)-(p4
     &   mup3+p5mup3)*F(417)+p2mup3*F(418)+p6mup3*F(419))
       Fa8(100) = p2mup3*F(421)+p6mup3*F(422)+(p4mup3+p5mup3)*F(423)+p
     &   1mup3*F(424)
       Fa8(101) = 2*((p4mup3+p5mup3)*F(335)+p6mup3*F(337))+p1mup3*F(35
     &   0)-p2mup3*F(352)
       Fa8(102) = -4*p2mup4*s45*F(274)+p4mup4*F(426)+p5mup4*F(427)+p6m
     &   up4*F(428)
       Fa8(103) = -4*p2mup4*F(370)+(p4mup4+p5mup4)*F(429)+p1mup4*F(430
     &   )
       Fa8(104) = p2mup4*F(431)+p6mup4*F(432)
       Return
       End
