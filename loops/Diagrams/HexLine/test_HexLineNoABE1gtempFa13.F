c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineNoAbe1g13(p1mup2,p1mup3,p1mup4
     -   ,p1mup6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3m
     -   up6,p4mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,
     -   p6mup2,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3m
     -   up4,mup3mup6,mup4mup6,Fa13)
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
       Complex*16 Fa13(169:182)
       COMMON/test_FaFunctionsNoAbe1g/Fa
       Fa13(169) = 2*p4mup2*s12*F(254)+(p1mup2+p2mup2)*s123*F(266)
       Fa13(170) = mup4mup6*(8*p3mup2*p4mup3*F(267)+p3mup3*p4mup2*F(26
     -   8))
       Fa13(171) = p3mup2*p4mup3*(p3mup4*p5mup6*F(245)+p2mup6*p3mup4*F
     -   (246)+p3mup4*p3mup6*F(246)+p3mup4*p4mup6*F(247)+p1mup4*(p5mu
     -   p6*F(245)+(p2mup6+p3mup6)*F(246)+p4mup6*F(247))+p2mup4*(p5mu
     -   p6*F(245)+(p2mup6+p3mup6)*F(246)+p4mup6*F(247))-p4mup4*p5mup
     -   6*F(248)+p2mup6*p4mup4*F(249)+p3mup6*p4mup4*F(249)+p4mup4*p4
     -   mup6*F(250))
       Fa13(172) = p1mup6*p3mup2*p4mup3*((p1mup4+p2mup4+p3mup4)*F(246)
     -   +p4mup4*F(249))
       Fa13(173) = p1mup3*p2mup6*(p1mup4+p2mup4+p3mup4)*p4mup2
       Fa13(174) = p4mup2*(p1mup3*(p3mup4*p5mup6*F(245)+p3mup4*p3mup6*
     -   F(246)+p3mup4*p4mup6*F(247)+p1mup4*(p5mup6*F(245)+p3mup6*F(2
     -   46)+p4mup6*F(247))+p2mup4*(p5mup6*F(245)+p3mup6*F(246)+p4mup
     -   6*F(247))-p4mup4*p5mup6*F(248)+p2mup6*p4mup4*F(249)+p3mup6*p
     -   4mup4*F(249)+p4mup4*p4mup6*F(250))+(p2mup3+p3mup3)*(p3mup4*p
     -   5mup6*F(245)+p2mup6*p3mup4*F(246)+p3mup4*p3mup6*F(246)+p3mup
     -   4*p4mup6*F(247)+p1mup4*(p5mup6*F(245)+(p2mup6+p3mup6)*F(246)
     -   +p4mup6*F(247))+p2mup4*(p5mup6*F(245)+(p2mup6+p3mup6)*F(246)
     -   +p4mup6*F(247))-p4mup4*p5mup6*F(248)+p2mup6*p4mup4*F(249)+p3
     -   mup6*p4mup4*F(249)+p4mup4*p4mup6*F(250)))
       Fa13(175) = p1mup6*p2mup3*p2mup4*p4mup2*F(246)+p1mup6*p2mup4*p3
     -   mup3*p4mup2*F(246)+p1mup4*p1mup6*(p2mup3+p3mup3)*p4mup2*F(24
     -   6)+p1mup6*p2mup3*p3mup4*p4mup2*F(246)+p1mup6*p3mup3*p3mup4*p
     -   4mup2*F(246)+p1mup6*p2mup3*p4mup2*p4mup4*F(249)+p1mup6*p3mup
     -   3*p4mup2*p4mup4*F(249)+p1mup3*p1mup6*p4mup2*((p1mup4+p2mup4+
     -   p3mup4)*F(246)+p4mup4*F(249))+p3mup3*p4mup2*p4mup6*p5mup4*F(
     -   253)-p3mup2*p4mup3*p4mup6*p5mup4*F(253)+(p3mup3*p4mup2-p3mup
     -   2*p4mup3)*p5mup4*p5mup6*F(256)
       Fa13(176) = (p1mup6*p2mup3+p1mup3*(p1mup6+p2mup6+p3mup6))*p4mup
     -   2*p5mup4
       Fa13(177) = p5mup4*(16*p4mup2*(-((p1mup3+p2mup3)*(p4mup6*F(253)
     -   +p5mup6*F(256)))+((p1mup6+p2mup6)*p3mup3+p2mup3*(p2mup6+p3mu
     -   p6))*F(257))+(p1mup6+p2mup6)*p3mup2*p4mup3*F(269))
       Fa13(178) = p3mup6*p5mup4*(16*p3mup3*p4mup2*F(257)+p3mup2*p4mup
     -   3*F(269))
       Fa13(179) = -4*p3mup2*p4mup3*F(254)+(p1mup3+p2mup3+p3mup3)*p4mu
     -   p2*F(270)
       Fa13(180) = (p1mup2+p2mup2)*p4mup3*(16*(p3mup4*p5mup6*F(245)+p1
     -   mup6*p3mup4*F(246)+p2mup6*p3mup4*F(246)+p3mup4*p3mup6*F(246)
     -   +p3mup4*p4mup6*F(247)+p1mup4*(p5mup6*F(245)+(p1mup6+p2mup6+p
     -   3mup6)*F(246)+p4mup6*F(247))+p2mup4*(p5mup6*F(245)+(p1mup6+p
     -   2mup6+p3mup6)*F(246)+p4mup6*F(247))-p4mup4*p5mup6*F(248)+p1m
     -   up6*p4mup4*F(249)+p2mup6*p4mup4*F(249)+p3mup6*p4mup4*F(249)+
     -   p4mup4*p4mup6*F(250)+p4mup6*p5mup4*F(253)+p5mup4*p5mup6*F(25
     -   6)-(p1mup6+p2mup6+p3mup6)*p5mup4*F(257))-mup4mup6*F(268))
       Fa13(181) = mup4mup6*(p2mup2*(p1mup3+p3mup3)+p1mup2*(p1mup3+p2m
     -   up3+p3mup3))
       Fa13(182) = p2mup2*(p1mup3*p1mup6*p2mup4*F(246)+p2mup4*p2mup6*p
     -   3mup3*F(246)+p1mup3*p1mup6*p3mup4*F(246)+p2mup6*p3mup3*p3mup
     -   4*F(246)+p2mup4*p3mup3*p4mup6*F(247)+p1mup4*((p1mup3*p1mup6+
     -   p2mup6*p3mup3)*F(246)+p2mup3*(p5mup6*F(245)+p3mup6*F(246))+(
     -   p2mup3+p3mup3)*p4mup6*F(247))+(p1mup3*p1mup6+p2mup6*p3mup3)*
     -   p4mup4*F(249)+p2mup3*(p2mup4*p5mup6*F(245)+p3mup4*p5mup6*F(2
     -   45)+p2mup4*p3mup6*F(246)+p3mup4*p3mup6*F(246)+p2mup4*p4mup6*
     -   F(247)+p3mup4*p4mup6*F(247)-p4mup4*p5mup6*F(248)+p3mup6*p4mu
     -   p4*F(249)+p4mup4*p4mup6*F(250)))+p1mup2*(p1mup3*p1mup6*p2mup
     -   4*F(246)+p2mup4*p2mup6*p3mup3*F(246)+p1mup3*p1mup6*p3mup4*F(
     -   246)+p2mup6*p3mup3*p3mup4*F(246)+p2mup4*p3mup3*p4mup6*F(247)
     -   +p3mup3*p3mup4*p4mup6*F(247)+p1mup4*((p1mup3*p1mup6+p2mup6*p
     -   3mup3)*F(246)+p2mup3*(p5mup6*F(245)+p3mup6*F(246))+(p2mup3+p
     -   3mup3)*p4mup6*F(247))+p1mup3*p1mup6*p4mup4*F(249)+p2mup6*p3m
     -   up3*p4mup4*F(249)+p3mup3*p4mup4*p4mup6*F(250)+p2mup3*(p2mup4
     -   *p5mup6*F(245)+p3mup4*p5mup6*F(245)+p2mup4*p3mup6*F(246)+p3m
     -   up4*p3mup6*F(246)+p2mup4*p4mup6*F(247)+p3mup4*p4mup6*F(247)-
     -   p4mup4*p5mup6*F(248)+p3mup6*p4mup4*F(249)+p4mup4*p4mup6*F(25
     -   0)))
       Return
       End
