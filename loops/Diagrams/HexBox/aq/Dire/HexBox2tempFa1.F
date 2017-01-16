c************** Calling the Fa functions*************************
       subroutine HexBox2Fa1(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6
     &   ,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa1)
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
       Complex*16 Fa1(1:19)
       COMMON/HexBox2FaFunctions/Fa
       Fa1(1) = p1mup6*F(5)+4*p6mup6*s16*F(6)
       Fa1(2) = p5mup6*s16*F(6)+p1mup6*s45*F(13)
       Fa1(3) = 4*p6mup6*s16*F(23)+p1mup6*F(24)
       Fa1(4) = 8*(p1mup3*s16*F(23)+p2mup3*s16*F(26))+p6mup3*F(27)
       Fa1(5) = p1mup6*(p1mup3+p6mup3)
       Fa1(6) = p1mup3+p6mup3
       Fa1(7) = p2mup3*s16*F(31)+p6mup3*F(32)
       Fa1(8) = p1mup6*((p1mup3+p6mup3)*F(33)+p5mup3*F(34)-p2mup3*F(35
     &   ))
       Fa1(9) = p5mup6*s16*F(23)-p1mup6*F(36)
       Fa1(10) = p1mup3*s16*F(23)+p6mup3*F(37)
       Fa1(11) = 4*mup3mup6*s16*F(23)+8*p1mup6*p2mup3*F(38)
       Fa1(12) = mup3mup6*s16*F(23)+2*p1mup6*p2mup3*F(38)
       Fa1(13) = -4*((p1mup3+p6mup3)*F(33)-p2mup3*F(35))+p5mup3*F(39)
       Fa1(14) = p5mup3*s16*F(23)-6*p6mup3*F(1633)
       Fa1(15) = -2*p5mup3*s16*F(23)+12*p6mup3*F(1633)
       Fa1(16) = p6mup3*F(42)-4*s16*(p1mup3*F(31)+p5mup3*F(43))
       Fa1(17) = p1mup6*(p1mup3+p5mup3+p6mup3)
       Fa1(18) = p1mup3+p5mup3+p6mup3
       Fa1(19) = p5mup6*F(56)+p6mup6*F(57)+p4mup6*F(58)-p1mup6*F(59)
       Return
       End
