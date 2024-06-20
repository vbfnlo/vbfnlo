c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineNoAbe1g4(p1mup2,p1mup3,p1mup4,
     -   p1mup6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mu
     -   p6,p4mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p
     -   6mup2,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mu
     -   p4,mup3mup6,mup4mup6,Fa4)
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
       Complex*16 Fa4(43:56)
       COMMON/test_FaFunctionsNoAbe1g/Fa
       Fa4(43) = mup4mup6*(p1mup2+p2mup2)*s45*F(10)+2*p5mup4*(p1mup2*(
     -   p3mup6*F(11)-(p1mup6+p2mup6)*s345*F(12)+(p4mup6+p5mup6)*s345
     -   *F(13))+p2mup2*(p3mup6*F(11)-(p1mup6+p2mup6)*s345*F(12)+(p4m
     -   up6+p5mup6)*s345*F(13))+p3mup2*s12*(-((p3mup6+p4mup6+p5mup6)
     -   *F(12))+(p1mup6+p2mup6)*F(14)))-3*mup4mup6*p4mup2*s12*s345*F
     -   (16)
       Fa4(44) = p1mup2*(mup3mup6*s45*F(10)-2*(p4mup6*(p5mup3*F(10)+p4
     -   mup3*F(11))-p4mup3*(-(p3mup6*F(11))+(p1mup6+p2mup6)*s345*F(1
     -   2)-p5mup6*s345*F(13))))+p2mup2*(mup3mup6*s45*F(10)-2*(p4mup6
     -   *(p5mup3*F(10)+p4mup3*F(11))-p4mup3*(-(p3mup6*F(11))+(p1mup6
     -   +p2mup6)*s345*F(12)-p5mup6*s345*F(13))))+2*(p3mup2*p4mup3+p4
     -   mup2*(p4mup3+p5mup3))*s12*((p3mup6+p4mup6+p5mup6)*F(12)-(p1m
     -   up6+p2mup6)*F(14))-3*mup3mup6*p4mup2*s12*s345*F(16)
       Fa4(45) = (p1mup2+p2mup2)*s45*F(10)-3*p4mup2*s12*s345*F(16)
       Fa4(46) = mup4mup6*(p3mup2+p4mup2)*p4mup3+mup4mup6*p4mup2*p5mup
     -   3+mup3mup6*p3mup2*p5mup4
       Fa4(47) = (p1mup2+p2mup2)*(mup4mup6*p4mup3+mup3mup6*p5mup4)
       Fa4(48) = (p1mup2+p2mup2)*p4mup3
       Fa4(49) = mup4mup6*(p1mup2+p2mup2)*(p4mup3+p5mup3)
       Fa4(50) = (p1mup2+p2mup2)*(p4mup3+p5mup3)
       Fa4(51) = (p1mup2+p2mup2)*(p4mup3+p5mup3)*p5mup4
       Fa4(52) = mup4mup6*(p1mup2+p2mup2)
       Fa4(53) = p1mup2+p2mup2
       Fa4(54) = (p1mup2+p2mup2)*p5mup4
       Fa4(55) = mup3mup6*(p1mup2+p2mup2)
       Fa4(56) = -(p3mup6*F(11))+(p1mup6+p2mup6)*s345*F(12)-(p4mup6+p5
     -   mup6)*s345*F(13)
       Return
       End
