c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineNoAbe1g2(p1mup2,p1mup3,p1mup4,
     -   p1mup6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mu
     -   p6,p4mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p
     -   6mup2,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mu
     -   p4,mup3mup6,mup4mup6,Fa2)
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
       Complex*16 Fa2(15:28)
       COMMON/test_FaFunctionsNoAbe1g/Fa
       Fa2(15) = p3mup2*p4mup3+p4mup2*(p4mup3+p5mup3)
       Fa2(16) = (p3mup2+p4mup2+p5mup2)*(mup4mup6*p4mup3+mup3mup6*p5mu
     -   p4)
       Fa2(17) = p4mup3*(p3mup2+p4mup2+p5mup2)
       Fa2(18) = mup4mup6*(p3mup2+p4mup2+p5mup2)*(p4mup3+p5mup3)
       Fa2(19) = (p3mup2+p4mup2+p5mup2)*(p4mup3+p5mup3)
       Fa2(20) = (p3mup2+p4mup2+p5mup2)*(p4mup3+p5mup3)*p5mup4
       Fa2(21) = mup4mup6*p4mup3+mup3mup6*p5mup4
       Fa2(22) = mup4mup6*(p4mup3+p5mup3)
       Fa2(23) = p4mup3+p5mup3
       Fa2(24) = (p4mup3+p5mup3)*p5mup4
       Fa2(25) = mup4mup6*(p3mup2+p4mup2+p5mup2)
       Fa2(26) = p3mup2+p4mup2+p5mup2
       Fa2(27) = (p3mup2+p4mup2+p5mup2)*p5mup4
       Fa2(28) = mup4mup6*p4mup2+mup2mup6*p5mup4
       Return
       End
