c************** Calling the Fa functions*************************
       subroutine HexBoxFa9(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,
     &   p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa9)
       IMPLICIT NONE
      Complex*16   p1mup3, p1mup4, p2mup3, 
     -          p2mup4, p3mup3, p3mup4, 
     -          p4mup3, p4mup4, p5mup3, 
     -          p5mup4, p6mup3, p6mup4
       Complex*16   mup3mup4

        common/HexBoxFhlFunctions/F
       COMMON/HexBoxInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12,s23,s
     &   34,s45,s56,s16,s123,s234,s345
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
      Complex*16 Fa(259),F(1314)
      Real*8 P(1362) 
        COMMON/HexBoxPFunctions/P
       Complex*16 Fa9(97:108)
       COMMON/HexBoxFaFunctions/Fa
       Fa9(97) = 2*(p1mup4*s45*F(288)+p6mup4*F(346))+(p4mup4+p5mup4)*F
     &   (404)
       Fa9(98) = -2*mup3mup4*s45*F(288)+4*(p4mup4+p5mup4)*(p2mup3*F(31
     &   1)+p1mup3*F(312)+(p4mup3+p5mup3)*F(315)+p6mup3*F(316))
       Fa9(99) = 4*(2*(-(p4mup4*p5mup3*F(279))-p5mup3*p5mup4*F(279)+p2
     &   mup3*p4mup4*F(282)+p2mup3*p5mup4*F(282)+p4mup4*p6mup3*F(284)
     &   +p5mup4*p6mup3*F(284)-p1mup4*p2mup3*F(311)-p1mup4*p5mup3*F(3
     &   15)-p1mup4*p6mup3*F(316)+p1mup3*((p4mup4+p5mup4)*F(273)-p1mu
     &   p4*F(312)+p6mup4*F(337))-p5mup3*p6mup4*F(338)-p4mup3*((p4mup
     &   4+p5mup4)*F(279)+p1mup4*F(315)+p6mup4*F(338))+p6mup3*p6mup4*
     &   F(339)+p2mup3*p6mup4*F(380))+mup3mup4*F(405))
       Fa9(100) = 4*(p1mup4*(p2mup3*F(311)+p1mup3*F(312)+(p4mup3+p5mup
     &   3)*F(315)+p6mup3*F(316))-p6mup4*(p1mup3*F(337)-(p4mup3+p5mup
     &   3)*F(338)+p6mup3*F(339)+p2mup3*F(380)))-mup3mup4*F(406)
       Fa9(101) = -4*p1mup3*(p4mup4+p5mup4)*F(337)-2*mup3mup4*F(346)+(
     &   p4mup4+p5mup4)*(4*((p4mup3+p5mup3)*F(338)-p6mup3*F(339))+p2m
     &   up3*F(407))
       Fa9(102) = (p4mup4+p5mup4)*F(287)+p1mup4*s45*F(399)+p6mup4*F(41
     &   1)
       Fa9(103) = 2*(p1mup3*F(412)+p2mup3*F(413)+p6mup3*F(414))+(p4mup
     &   3+p5mup3)*F(415)
       Fa9(104) = 8*(p4mup3*p4mup4*F(293)+p4mup4*p5mup3*F(293)+p4mup3*
     &   p5mup4*F(293)+p5mup3*p5mup4*F(293)+p4mup4*p6mup3*F(294)+p5mu
     &   p4*p6mup3*F(294)+p2mup4*p4mup3*F(315)+p2mup4*p5mup3*F(315)+p
     &   2mup4*p6mup3*F(316)+p4mup3*p6mup4*F(329)+p5mup3*p6mup4*F(329
     &   )+p1mup3*(-((p4mup4+p5mup4)*F(292))+p2mup4*F(312)+p6mup4*F(3
     &   30))-p6mup3*p6mup4*F(332)+p2mup3*(-((p4mup4+p5mup4)*F(291))+
     &   p2mup4*F(311)+p6mup4*F(377)))+mup3mup4*F(416)
       Fa9(105) = 2*(p4mup4+p5mup4)*(-(p2mup3*F(291))-p1mup3*F(292)+(p
     &   4mup3+p5mup3)*F(293)+p6mup3*F(294))+mup3mup4*s45*F(417)
       Fa9(106) = mup3mup4*F(287)-2*(p1mup4*p2mup3*F(291)-p1mup3*(p6mu
     &   p4*F(290)-p1mup4*F(292))-p6mup3*(p6mup4*F(281)+p1mup4*F(294)
     &   )-p2mup3*p6mup4*F(295)-(p4mup3+p5mup3)*(p1mup4*F(293)-p6mup4
     &   *F(379)))
       Fa9(107) = -2*mup3mup4*F(411)+(p4mup4+p5mup4)*(4*(p6mup3*F(281)
     &   +p1mup3*F(290)+p2mup3*F(295))+(p4mup3+p5mup3)*F(418))
       Fa9(108) = (p4mup4+p5mup4)*F(287)-p1mup4*s45*F(396)+p6mup4*F(39
     &   7)
       Return
       End
