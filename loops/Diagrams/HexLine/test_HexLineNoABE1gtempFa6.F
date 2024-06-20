c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineNoAbe1g6(p1mup2,p1mup3,p1mup4,
     -   p1mup6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mu
     -   p6,p4mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p
     -   6mup2,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mu
     -   p4,mup3mup6,mup4mup6,Fa6)
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
       Complex*16 Fa6(71:84)
       COMMON/test_FaFunctionsNoAbe1g/Fa
       Fa6(71) = (p1mup2+p2mup2)*(p1mup3+p2mup3+p3mup3)
       Fa6(72) = (p1mup2+p2mup2)*(p1mup3+p2mup3+p3mup3)*p5mup4
       Fa6(73) = mup4mup6*(p1mup3+p2mup3+p3mup3)
       Fa6(74) = p1mup3+p2mup3+p3mup3
       Fa6(75) = (p1mup3+p2mup3+p3mup3)*p5mup4
       Fa6(76) = -((p4mup6+p5mup6)*F(18))+(p1mup6+p2mup6+p3mup6)*F(19)
       Fa6(77) = (p1mup3+p2mup3+p3mup3)*((p1mup6+p2mup6+p3mup6+p4mup6)
     -   *F(24)-p5mup6*F(25))
       Fa6(78) = (p1mup6+p2mup6+p3mup6+p4mup6)*F(24)-p5mup6*F(25)
       Fa6(79) = p3mup2*(p1mup6+p2mup6+p3mup6+p4mup6)*F(24)-p3mup2*p5m
     -   up6*F(25)
       Fa6(80) = (p1mup2+p2mup2)*((p1mup6+p2mup6+p3mup6+p4mup6)*F(24)-
     -   p5mup6*F(25))
       Fa6(81) = (p1mup2+p2mup2)*(p1mup3+p2mup3+p3mup3)*((p1mup6+p2mup
     -   6+p3mup6+p4mup6)*F(24)-p5mup6*F(25))
       Fa6(82) = p1mup4+p2mup4+p3mup4+p4mup4
       Fa6(83) = p4mup2*s12+(p1mup2+p2mup2)*s123
       Fa6(84) = -4*p3mup2*p4mup3*F(27)+(p1mup3+p2mup3+p3mup3)*p4mup2*
     -   F(29)
       Return
       End
