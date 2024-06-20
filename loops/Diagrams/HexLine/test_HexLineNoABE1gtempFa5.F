c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineNoAbe1g5(p1mup2,p1mup3,p1mup4,
     -   p1mup6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mu
     -   p6,p4mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p
     -   6mup2,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mu
     -   p4,mup3mup6,mup4mup6,Fa5)
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
       Complex*16 Fa5(57:70)
       COMMON/test_FaFunctionsNoAbe1g/Fa
       Fa5(57) = (p3mup6+p4mup6+p5mup6)*F(12)-(p1mup6+p2mup6)*F(14)
       Fa5(58) = 6*mup4mup6*s45*F(17)-4*(p5mup4*(p4mup6+p5mup6)*F(18)-
     -   (p1mup6+p2mup6+p3mup6)*p5mup4*F(19))
       Fa5(59) = p4mup3*s123*((p4mup6+p5mup6)*F(18)-(p1mup6+p2mup6+p3m
     -   up6)*F(19))+(p1mup3+p2mup3+p3mup3)*((p1mup6+p2mup6+p3mup6)*s
     -   45*F(18)-p4mup6*F(20)-p5mup6*s45*F(21))
       Fa5(60) = -2*((p1mup6+p2mup6+p3mup6)*s45*F(18)-p5mup6*s45*F(21)
     -   )+p4mup6*F(22)
       Fa5(61) = p4mup2*s123*((p4mup6+p5mup6)*F(18)-(p1mup6+p2mup6+p3m
     -   up6)*F(19))+p3mup2*((p1mup6+p2mup6+p3mup6)*s45*F(18)-p4mup6*
     -   F(20)-p5mup6*s45*F(21))
       Fa5(62) = -((p1mup2+p2mup2)*((p1mup6+p2mup6+p3mup6)*s45*F(18)-p
     -   4mup6*F(20)-p5mup6*s45*F(21)))
       Fa5(63) = (p1mup2+p2mup2)*(-((p4mup6+p5mup6)*F(18))+(p1mup6+p2m
     -   up6+p3mup6)*F(19))
       Fa5(64) = (p1mup2+p2mup2)*s123*(-3*mup4mup6*s45*F(17)+2*(p5mup4
     -   *(p4mup6+p5mup6)*F(18)-(p1mup6+p2mup6+p3mup6)*p5mup4*F(19)))
     -   +mup4mup6*p4mup2*s12*F(23)
       Fa5(65) = (p1mup2+p2mup2)*(p4mup3*s123*((p4mup6+p5mup6)*F(18)-(
     -   p1mup6+p2mup6+p3mup6)*F(19))+(p1mup3+p2mup3+p3mup3)*((p1mup6
     -   +p2mup6+p3mup6)*s45*F(18)-p4mup6*F(20)-p5mup6*s45*F(21)))
       Fa5(66) = 6*(p1mup2+p2mup2)*s123*s45*F(17)-2*p4mup2*s12*F(23)
       Fa5(67) = -(mup4mup6*(p1mup3+p2mup3+p3mup3)*p4mup2)+mup4mup6*p3
     -   mup2*p4mup3
       Fa5(68) = (p1mup3+p2mup3+p3mup3)*p4mup2-p3mup2*p4mup3
       Fa5(69) = mup4mup6*(p1mup2+p2mup2)*p4mup3
       Fa5(70) = mup4mup6*(p1mup2+p2mup2)*(p1mup3+p2mup3+p3mup3)
       Return
       End
