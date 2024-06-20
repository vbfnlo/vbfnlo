c************** Calling the Fa functions*************************
       subroutine HexBoxCross2Fa2(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p
     &   3mup6,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa2
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
       Complex*16 Fa2(20:38)
       COMMON/HexBoxCross2FaFunctions/Fa
       Fa2(20) = p2mup3*s16*F(49)+(p1mup3+p6mup3)*F(67)
       Fa2(21) = -8*((p1mup3+p6mup3)*F(42)+p2mup3*F(57))+p5mup3*F(68)
       Fa2(22) = p2mup6*F(50)-2*p5mup6*F(63)+(p1mup6+p6mup6)*F(71)
       Fa2(23) = 8*p2mup6*s16*F(49)+(p1mup6+p6mup6)*F(72)
       Fa2(24) = (p1mup6+p6mup6)*F(42)+p2mup6*F(73)-2*p5mup6*F(74)
       Fa2(25) = 2*p2mup3*(p1mup6+p6mup6)*F(46)+mup3mup6*F(63)
       Fa2(26) = 2*(-(p2mup6*(p1mup3+p6mup3))+p2mup3*(p1mup6+p6mup6))*
     &   F(46)+mup3mup6*F(74)
       Fa2(27) = p6mup6*F(60)+p2mup6*F(62)-p5mup6*F(63)+p1mup6*F(76)
       Fa2(28) = p1mup6*F(77)+p6mup6*F(78)+4*p2mup6*F(79)
       Fa2(29) = 4*((p1mup6+p6mup6)*F(60)-p5mup6*F(63))+p2mup6*F(80)
       Fa2(30) = 4*(p1mup3*F(62)+p5mup3*F(79))+p6mup3*F(80)
       Fa2(31) = 4*(p2mup3*F(52)+(p1mup3+p6mup3)*F(60))+p5mup3*F(82)
       Fa2(32) = p2mup6*p5mup3*F(46)+3*mup3mup6*F(1841)
       Fa2(33) = p1mup3+p6mup3
       Fa2(34) = p5mup6*(p1mup3+p6mup3)-p5mup3*(p1mup6+p6mup6)
       Fa2(35) = p5mup3*F(74)+6*(p1mup3+p6mup3)*F(1841)
       Fa2(36) = p2mup6*F(64)+6*(p1mup6+p5mup6+p6mup6)*F(1841)
       Fa2(37) = p1mup3+p5mup3+p6mup3
       Fa2(38) = p2mup6*(p1mup3+p6mup3)
       Return
       End
