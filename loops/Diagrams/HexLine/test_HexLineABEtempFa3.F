c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineAbe3(p1mup2,p1mup3,p1mup4,p1mu
     -   p5,p2mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p
     -   4mup2,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup
     -   2,p6mup3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,m
     -   up3mup5,mup4mup5,Fa3)
       IMPLICIT NONE
      Complex*16   p1mup2, p1mup3, p1mup4, p1mup5, p2mup2, p2mup3, 
     -          p2mup4, p2mup5, p3mup2, p3mup3, p3mup4, p3mup5, 
     -          p4mup2, p4mup3, p4mup4, p4mup5, p5mup2, p5mup3, 
     -          p5mup4, p5mup5, p6mup2, p6mup3, p6mup4, p6mup5
       Complex*16   mup2mup3, mup2mup4, mup2mup5, mup3mup4, mup3mup5, 
     -          mup4mup5
        common/test_FhlFunctionsAbe/F
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
      Complex*16 Fa(445),F(5930)
      Real*8 P(11674) 
        COMMON/test_PFunctionsAbe/P
       Complex*16 Fa3(23:33)
       COMMON/test_FaFunctionsAbe/Fa
       Fa3(23) = -(p4mup2*s12*((p1mup3+p2mup3)*F(14)+p3mup3*F(15)))+p1
     -   mup2*(p4mup3*F(12)+(p1mup3+p2mup3)*s123*F(13)+p3mup3*s123*F(
     -   16))+p2mup2*(p4mup3*F(12)+(p1mup3+p2mup3)*s123*F(13)+p3mup3*
     -   s123*F(16))
       Fa3(24) = p1mup3*p3mup2*s123*F(13)+p2mup3*p3mup2*s123*F(13)+p3m
     -   up2*p3mup3*s123*F(16)+p3mup2*p4mup3*F(17)+p1mup3*p4mup2*F(18
     -   )+p2mup3*p4mup2*F(18)-p3mup3*p4mup2*F(19)
       Fa3(25) = p4mup3*F(12)+(p1mup3+p2mup3)*s123*F(13)+p3mup3*s123*F
     -   (16)
       Fa3(26) = (p1mup4+p2mup4+p3mup4+p4mup4)*(p1mup3*p3mup2*s12*F(14
     -   )+p2mup3*p3mup2*s12*F(14)+p3mup2*p3mup3*s12*F(15)+p1mup3*(p1
     -   mup2+p2mup2)*F(18)+(p1mup2+p2mup2)*p2mup3*F(18)-(p1mup2+p2mu
     -   p2)*p3mup3*F(19))
       Fa3(27) = p1mup3*p3mup2*s12*F(14)+p2mup3*p3mup2*s12*F(14)+p3mup
     -   2*p3mup3*s12*F(15)+p1mup3*(p1mup2+p2mup2)*F(18)+(p1mup2+p2mu
     -   p2)*p2mup3*F(18)+(p1mup2+p2mup2)*p3mup3*F(20)
       Fa3(28) = (p1mup4+p2mup4+p3mup4+p4mup4)*((p1mup3+p2mup3)*F(18)+
     -   p3mup3*F(20))
       Fa3(29) = (p1mup3+p2mup3)*F(18)+p3mup3*F(20)
       Fa3(30) = (p1mup4+p2mup4+p3mup4+p4mup4)*((p1mup3+p2mup3)*F(14)+
     -   p3mup3*F(15))
       Fa3(31) = (p1mup3+p2mup3)*F(14)+p3mup3*F(15)
       Fa3(32) = p4mup2*s12*F(21)-(p1mup2+p2mup2)*s123*F(22)
       Fa3(33) = (p1mup3+p2mup3+p3mup3)*p4mup2*F(24)-p3mup2*p4mup3*F(2
     -   5)
       Return
       End
