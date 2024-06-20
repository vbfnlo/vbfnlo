c************** Calling the Fa functions*************************
       subroutine HexBoxCross1Fa9(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p
     &   3mup4,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa9
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
       Complex*16 Fa9(105:117)
       COMMON/HexBoxCross1FaFunctions/Fa
       Fa9(105) = p2mup4*F(382)-p1mup4*F(395)
       Fa9(106) = mup3mup4*F(433)+16*(-(p4mup4*p5mup3*F(342))-p5mup3*p
     &   5mup4*F(342)+p2mup3*p4mup4*F(343)+p2mup3*p5mup4*F(343)+p4mup
     &   4*p6mup3*F(344)+p5mup4*p6mup3*F(344)-p1mup4*p2mup3*F(360)-p1
     &   mup4*p5mup3*F(434)+p1mup3*((p4mup4+p5mup4)*F(341)-p1mup4*F(3
     &   62)+p6mup4*F(435))-p1mup4*p6mup3*F(436)-p5mup3*p6mup4*F(437)
     &   -p4mup3*((p4mup4+p5mup4)*F(342)+p1mup4*F(434)+p6mup4*F(437))
     &   +p2mup3*p6mup4*F(438)+p6mup3*p6mup4*F(439))
       Fa9(107) = 2*(p2mup3*F(360)+p1mup3*F(362))+(p4mup3+p5mup3)*F(44
     &   0)+p6mup3*F(441)
       Fa9(108) = p1mup3*F(435)-(p4mup3+p5mup3)*F(437)+p2mup3*F(438)+p
     &   6mup3*F(439)
       Fa9(109) = p2mup3*F(360)+p1mup3*F(362)+(p4mup3+p5mup3)*F(364)+p
     &   6mup3*F(365)
       Fa9(110) = p1mup3*F(341)-(p4mup3+p5mup3)*F(342)+p2mup3*F(343)+p
     &   6mup3*F(344)
       Fa9(111) = -(p1mup3*F(376))+(p4mup3+p5mup3)*F(377)-p6mup3*F(378
     &   )+p2mup3*F(379)
       Fa9(112) = (p4mup4+p5mup4)*F(431)+p6mup4*F(442)
       Fa9(113) = (p4mup4+p5mup4)*F(382)-p1mup4*F(443)
       Fa9(114) = p6mup4*F(395)+p2mup4*F(431)
       Fa9(115) = p2mup4*F(382)+p1mup4*F(444)
       Fa9(116) = p6mup4*F(382)+p1mup4*F(431)
       Fa9(117) = (p4mup4+p5mup4)*F(274)+p6mup4*F(275)
       Return
       End
