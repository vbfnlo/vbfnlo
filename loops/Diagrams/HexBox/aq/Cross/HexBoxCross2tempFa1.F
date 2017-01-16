c************** Calling the Fa functions*************************
       subroutine HexBoxCross2Fa1(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p
     &   3mup6,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa1
     &   )
       IMPLICIT NONE
      Complex*16   p1mup3, p1mup6, p2mup3, p2mup6, p3mup3, p3mup6, 
     -          p4mup3, p4mup6, p5mup3, p5mup6, p6mup3, p6mup6
       Complex*16   mup3mup6
        common/HexBoxCross2FhlFunctions/F
       COMMON/HexBoxCross2Invariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12
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
      Complex*16 Fa(392),F(1853)
      Real*8 P(1355) 
        COMMON/HexBoxCross2PFunctions/P
       Complex*16 Fa1(1:19)
       COMMON/HexBoxCross2FaFunctions/Fa
       Fa1(1) = (p1mup6+p6mup6)*F(1)+p5mup6*s16*F(2)
       Fa1(2) = 8*((p1mup6+p6mup6)*F(3)+p4mup6*F(4))+p5mup6*F(5)
       Fa1(3) = p1mup6*F(23)-4*(p5mup6*F(21)+p6mup6*F(24)+p4mup6*F(25)
     &   )
       Fa1(4) = p1mup6*F(26)+p6mup6*F(27)-4*p4mup6*F(28)
       Fa1(5) = -4*(p5mup6*F(21)+(p1mup6+p6mup6)*F(24))+p4mup6*F(29)
       Fa1(6) = p1mup6+p6mup6
       Fa1(7) = p4mup6*F(22)+6*(p1mup6+p5mup6+p6mup6)*F(1842)
       Fa1(8) = p5mup6*F(43)-8*((p1mup6+p6mup6)*F(42)-p2mup6*F(44))
       Fa1(9) = p2mup6*p5mup3*F(45)-p2mup3*p5mup6*F(46)
       Fa1(10) = (p1mup6+p6mup6)*F(48)-p5mup6*s16*F(49)
       Fa1(11) = p5mup3*s16*F(49)+(p1mup3+p6mup3)*F(50)
       Fa1(12) = 8*(p1mup3+p6mup3)*F(42)+p5mup3*F(51)+16*p2mup3*F(52)
       Fa1(13) = p5mup3*(p1mup6+p6mup6)
       Fa1(14) = -2*p2mup3*(p1mup6+p6mup6)*F(46)+p2mup6*(p1mup3+p6mup3
     &   )*F(54)
       Fa1(15) = -2*p2mup3*p5mup6*F(46)+p2mup6*p5mup3*F(54)
       Fa1(16) = (p1mup3+p6mup3)*(p1mup6+p6mup6)
       Fa1(17) = -2*p5mup6*(p1mup3+p6mup3)*F(46)+p5mup3*(p1mup6+p6mup6
     &   )*F(54)
       Fa1(18) = p2mup6*p5mup3*s16*F(49)+(p1mup3+p6mup3)*(p1mup6+p6mup
     &   6)*F(60)+p2mup6*(p1mup3+p6mup3)*F(62)-p5mup6*(p1mup3+p6mup3)
     &   *F(63)+(p1mup6+p6mup6)*(p5mup3*F(61)-p2mup3*F(64))
       Fa1(19) = -(p2mup3*F(50))+2*p5mup3*F(66)
       Return
       End
