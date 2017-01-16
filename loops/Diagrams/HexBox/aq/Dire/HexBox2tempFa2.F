c************** Calling the Fa functions*************************
       subroutine HexBox2Fa2(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6
     &   ,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa2)
       IMPLICIT NONE
      Complex*16   p1mup3, p1mup6, p2mup3, p2mup6, p3mup3, p3mup6, 
     -          p4mup3, p4mup6, p5mup3, p5mup6, p6mup3, p6mup6
       Complex*16   mup3mup6
        common/HexBox2FhlFunctions/F
       COMMON/HexBox2Invariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12,s23,
     &   s34,s45,s56,s16,s123,s234,s345
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
      Complex*16 Fa(385),F(1646)
      Real*8 P(1185) 
        COMMON/HexBox2PFunctions/P
       Complex*16 Fa2(20:38)
       COMMON/HexBox2FaFunctions/Fa
       Fa2(20) = (p5mup6+p6mup6)*F(64)+p1mup6*s56*F(65)
       Fa2(21) = p5mup6+p6mup6
       Fa2(22) = p4mup6*F(68)+(p5mup6+p6mup6)*F(69)+p1mup6*F(70)
       Fa2(23) = p5mup6*F(76)+p6mup6*F(77)-p4mup6*s56*F(78)
       Fa2(24) = -4*((p5mup6+p6mup6)*F(57)-p1mup6*F(59))+p4mup6*F(79)
       Fa2(25) = p4mup6*F(83)-6*(p1mup6+p5mup6+p6mup6)*F(1632)
       Fa2(26) = p4mup6*F(87)+6*(p1mup6+p5mup6+p6mup6)*F(1632)
       Fa2(27) = p5mup6*F(95)+p6mup6*s12*F(96)-p2mup6*F(97)-p1mup6*F(9
     &   8)
       Fa2(28) = 16*(-2*p1mup6*p2mup3*F(99)+p1mup3*p2mup6*F(100))
       Fa2(29) = p2mup3*F(102)-p1mup3*F(103)
       Fa2(30) = (p5mup6+p6mup6)*F(104)+p1mup6*s56*F(105)
       Fa2(31) = p1mup3*(p5mup6+p6mup6)
       Fa2(32) = (p5mup6+p6mup6)*F(108)-p2mup6*F(109)+p1mup6*F(110)
       Fa2(33) = 8*p2mup3*s56*F(105)+(p5mup3+p6mup3)*F(112)
       Fa2(34) = 8*((p5mup3+p6mup3)*F(108)+p1mup3*F(113))+p2mup3*F(114
     &   )
       Fa2(35) = p6mup3*p6mup6*s12*F(96)+p1mup3*p2mup6*s56*F(105)+p1mu
     &   p3*p6mup6*F(116)-p2mup3*p6mup6*F(117)+p5mup6*(p6mup3*s12*F(9
     &   6)+p1mup3*F(116)-p2mup3*F(117))-p2mup6*p6mup3*F(118)+p1mup6*
     &   p6mup3*F(119)+p5mup3*((p5mup6+p6mup6)*s12*F(96)-p2mup6*F(118
     &   )+p1mup6*F(119))
       Fa2(36) = p5mup6*F(120)+p6mup6*F(121)+p2mup6*s56*F(122)
       Fa2(37) = -4*((p5mup6+p6mup6)*s12*F(96)-p1mup6*F(98))+p2mup6*F(
     &   123)
       Fa2(38) = 8*p2mup3*(p5mup6+p6mup6)*F(99)+4*mup3mup6*F(119)
       Return
       End
