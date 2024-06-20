c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineNoAbe1g1(p1mup2,p1mup3,p1mup4,
     -   p1mup6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mu
     -   p6,p4mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p
     -   6mup2,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mu
     -   p4,mup3mup6,mup4mup6,Fa1)
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
       Complex*16 Fa1(1:14)
       COMMON/test_FaFunctionsNoAbe1g/Fa
       Fa1(1) = 2*p3mup6*p5mup4+mup4mup6*s45
       Fa1(2) = -4*(p3mup6*p4mup3+p4mup6*(p4mup3+p5mup3))*F(1)+mup3mup
     -   6*s45*F(2)
       Fa1(3) = 4*p2mup6*s45*F(3)+p4mup6*F(4)+2*s16*s45*(-(p1mup6*F(5)
     -   )+(p3mup6+p5mup6)*F(6))
       Fa1(4) = -2*(p3mup6*p4mup2-p3mup2*p4mup6)+mup2mup6*s45
       Fa1(5) = p5mup4*(2*p2mup6*F(3)-p1mup6*s16*F(5)+(p3mup6+p4mup6+p
     -   5mup6)*s16*F(6))
       Fa1(6) = p4mup3*(2*p2mup6*F(3)-p1mup6*s16*F(5)+(p3mup6+p4mup6+p
     -   5mup6)*s16*F(6))
       Fa1(7) = p3mup2*p4mup6*F(1)+p4mup6*p5mup2*F(1)+2*(p2mup6*p4mup2
     -   *F(3)+p4mup2*p4mup6*F(3))-p1mup6*p4mup2*s16*F(5)+p4mup2*(p3m
     -   up6+p5mup6)*s16*F(6)
       Fa1(8) = mup3mup6*s345*F(2)+4*(p4mup3+p5mup3)*(2*p2mup6*F(3)-p1
     -   mup6*s16*F(5)+(p3mup6+p4mup6+p5mup6)*s16*F(6))
       Fa1(9) = -2*mup2mup6*s345*F(1)+8*p3mup2*p3mup6*F(3)+4*(p3mup6*(
     -   p4mup2+p5mup2)*F(1)+p3mup2*(2*p2mup6*F(3)-p1mup6*s16*F(5)+(p
     -   4mup6+p5mup6)*s16*F(6)))
       Fa1(10) = mup4mup6*p4mup3*s345*F(1)+p5mup4*(mup3mup6*s345*F(1)+
     -   2*(p4mup3+p5mup3)*(2*p2mup6*F(3)-p1mup6*s16*F(5)+(p3mup6+p4m
     -   up6+p5mup6)*s16*F(6)))
       Fa1(11) = 4*(-(mup2mup6*p5mup4*s345*F(1))+mup4mup6*p3mup2*s45*F
     -   (1)+mup4mup6*p5mup2*s45*F(1)+4*p2mup6*p3mup2*p5mup4*F(3)+2*(
     -   p3mup6*p5mup4*((p4mup2+p5mup2)*F(1)+2*p3mup2*F(3))-p1mup6*p3
     -   mup2*p5mup4*s16*F(5)+p3mup2*p5mup4*(p4mup6+p5mup6)*s16*F(6))
     -   )-mup4mup6*p4mup2*F(7)
       Fa1(12) = -(mup3mup6*p3mup2*s45*F(1))-mup3mup6*p5mup2*s45*F(1)+
     -   4*(p2mup6*p4mup2*p5mup3*F(3)+p4mup2*p4mup6*p5mup3*F(3))+2*(p
     -   3mup2*p4mup6*p5mup3*F(1)+p4mup6*p5mup2*p5mup3*F(1)-p1mup6*p4
     -   mup2*p5mup3*s16*F(5)+p4mup2*p5mup3*p5mup6*s16*F(6)+p3mup6*(p
     -   4mup3*p5mup2*F(1)+2*(p3mup2+p4mup2)*p4mup3*F(3)+p4mup2*p5mup
     -   3*s16*F(6)))+p4mup3*(-(mup2mup6*s345*F(1))+4*(p3mup2+p4mup2)
     -   *p4mup6*F(3)+2*(p4mup6*p5mup2*F(1)+(p3mup2+p4mup2)*(2*p2mup6
     -   *F(3)-p1mup6*s16*F(5)+p5mup6*s16*F(6))))+mup3mup6*p4mup2*F(8
     -   )
       Fa1(13) = (p3mup2+p5mup2)*s45*F(2)-p4mup2*F(9)
       Fa1(14) = mup4mup6*(p3mup2+p4mup2)*p4mup3+mup4mup6*p4mup2*p5mup
     -   3+(mup3mup6*p3mup2+mup2mup6*(p4mup3+p5mup3))*p5mup4
       Return
       End
