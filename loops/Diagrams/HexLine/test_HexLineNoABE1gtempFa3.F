c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineNoAbe1g3(p1mup2,p1mup3,p1mup4,
     -   p1mup6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mu
     -   p6,p4mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p
     -   6mup2,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mu
     -   p4,mup3mup6,mup4mup6,Fa3)
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
       Complex*16 Fa3(29:42)
       COMMON/test_FaFunctionsNoAbe1g/Fa
       Fa3(29) = mup3mup6*p3mup2+mup2mup6*(p4mup3+p5mup3)
       Fa3(30) = mup3mup6*(p3mup2+p4mup2+p5mup2)
       Fa3(31) = -(mup3mup6*p4mup2)+mup2mup6*p4mup3
       Fa3(32) = 2*p2mup6*F(3)-p1mup6*s16*F(5)+(p3mup6+p4mup6+p5mup6)*
     -   s16*F(6)
       Fa3(33) = mup4mup6*s45*F(10)+2*p5mup4*(p3mup6*F(11)-(p1mup6+p2m
     -   up6)*s345*F(12)+(p4mup6+p5mup6)*s345*F(13))
       Fa3(34) = 2*mup3mup6*s45*F(10)-4*(p4mup6*(p5mup3*F(10)+p4mup3*F
     -   (11))-p4mup3*(-(p3mup6*F(11))+(p1mup6+p2mup6)*s345*F(12)-p5m
     -   up6*s345*F(13)))
       Fa3(35) = 2*((p3mup6+p5mup6)*F(12)-(p1mup6+p2mup6)*F(14))+p4mup
     -   6*s345*F(15)
       Fa3(36) = -(p3mup2*p4mup6*F(10))+p4mup2*(p3mup6*F(11)-(p1mup6+p
     -   2mup6)*s345*F(12)+(p4mup6+p5mup6)*s345*F(13))
       Fa3(37) = p5mup4*(p3mup6+p4mup6+p5mup6)*F(12)-(p1mup6+p2mup6)*p
     -   5mup4*F(14)
       Fa3(38) = p4mup3*(p3mup6+p4mup6+p5mup6)*F(12)-(p1mup6+p2mup6)*p
     -   4mup3*F(14)
       Fa3(39) = p1mup2*p4mup6*F(10)+p2mup2*p4mup6*F(10)+p4mup2*s12*(-
     -   ((p3mup6+p4mup6+p5mup6)*F(12))+(p1mup6+p2mup6)*F(14))
       Fa3(40) = 4*(p4mup3+p5mup3)*((p3mup6+p4mup6+p5mup6)*F(12)-(p1mu
     -   p6+p2mup6)*F(14))-6*mup3mup6*s345*F(16)
       Fa3(41) = p1mup2*(p3mup6*F(11)-(p1mup6+p2mup6)*s345*F(12)+(p4mu
     -   p6+p5mup6)*s345*F(13))+p2mup2*(p3mup6*F(11)-(p1mup6+p2mup6)*
     -   s345*F(12)+(p4mup6+p5mup6)*s345*F(13))+p3mup2*s12*(-((p3mup6
     -   +p4mup6+p5mup6)*F(12))+(p1mup6+p2mup6)*F(14))
       Fa3(42) = 2*(p4mup3+p5mup3)*p5mup4*(-((p3mup6+p4mup6+p5mup6)*F(
     -   12))+(p1mup6+p2mup6)*F(14))+3*(mup4mup6*p4mup3+mup3mup6*p5mu
     -   p4)*s345*F(16)
       Return
       End
