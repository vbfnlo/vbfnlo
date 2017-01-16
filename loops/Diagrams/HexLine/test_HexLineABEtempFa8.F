c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineAbe8(p1mup2,p1mup3,p1mup4,p1mu
     -   p5,p2mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p
     -   4mup2,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup
     -   2,p6mup3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,m
     -   up3mup5,mup4mup5,Fa8)
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
       Complex*16 Fa8(78:88)
       COMMON/test_FaFunctionsAbe/Fa
       Fa8(78) = 8*(p1mup3*p3mup4*F(141)+p2mup3*p3mup4*F(141)+p3mup3*p
     -   3mup4*F(142)+p3mup3*p4mup4*F(143)+p1mup3*p4mup4*F(149)+p2mup
     -   3*p4mup4*F(149)+p1mup4*(p3mup3*F(141)+(p1mup3+p2mup3)*F(147)
     -   +p4mup3*F(152))+p2mup4*(p3mup3*F(141)+(p1mup3+p2mup3)*F(147)
     -   +p4mup3*F(152))+p3mup4*p4mup3*F(153)+p4mup3*p4mup4*F(155))+m
     -   up3mup4*F(183)
       Fa8(79) = 8*((p1mup3+p2mup3)*(p4mup4*F(145)+(p1mup4+p2mup4)*F(1
     -   48))+p1mup3*p3mup4*F(150)+p2mup3*p3mup4*F(150)+p3mup3*(p3mup
     -   4*F(144)+p4mup4*F(146)+(p1mup4+p2mup4)*F(151))+(p1mup4+p2mup
     -   4)*p4mup3*F(154)+p3mup4*p4mup3*F(156)+p4mup3*p4mup4*F(157))+
     -   mup3mup4*F(184)
       Fa8(80) = p1mup2*(p4mup4*F(163)+(p1mup4+p2mup4)*F(164)+p3mup4*F
     -   (165))+p2mup2*(p4mup4*F(163)+(p1mup4+p2mup4)*F(164)+p3mup4*F
     -   (165))+p3mup2*s12*(p4mup4*F(174)+(p1mup4+p2mup4)*F(175)+p3mu
     -   p4*F(176))
       Fa8(81) = p4mup4*F(163)+(p1mup4+p2mup4)*F(164)+p3mup4*F(165)
       Fa8(82) = p4mup2*(mup3mup4*F(185)+2*((p1mup3+p2mup3)*(p4mup4*F(
     -   187)+(p1mup4+p2mup4)*F(189))+p1mup3*p3mup4*F(190)+p2mup3*p3m
     -   up4*F(190)+p3mup3*(p3mup4*F(186)+p4mup4*F(188)+(p1mup4+p2mup
     -   4)*F(191))-(p1mup4+p2mup4)*p4mup3*F(192)-p3mup4*p4mup3*F(193
     -   )-p4mup3*p4mup4*F(194)))
       Fa8(83) = p3mup2*(mup3mup4*F(185)+2*((p1mup3+p2mup3)*(p4mup4*F(
     -   187)+(p1mup4+p2mup4)*F(189))+p1mup3*p3mup4*F(190)+p2mup3*p3m
     -   up4*F(190)+p3mup3*(p3mup4*F(186)+p4mup4*F(188)+(p1mup4+p2mup
     -   4)*F(191))-(p1mup4+p2mup4)*p4mup3*F(192)-p3mup4*p4mup3*F(193
     -   )-p4mup3*p4mup4*F(194)))
       Fa8(84) = (p1mup2+p2mup2)*(mup3mup4*F(185)+2*((p1mup3+p2mup3)*(
     -   p4mup4*F(187)+(p1mup4+p2mup4)*F(189))+p1mup3*p3mup4*F(190)+p
     -   2mup3*p3mup4*F(190)+p3mup3*(p3mup4*F(186)+p4mup4*F(188)+(p1m
     -   up4+p2mup4)*F(191))-(p1mup4+p2mup4)*p4mup3*F(192)-p3mup4*p4m
     -   up3*F(193)-p4mup3*p4mup4*F(194)))
       Fa8(85) = 4*mup3mup4*F(185)+8*((p1mup3+p2mup3)*(p4mup4*F(187)+(
     -   p1mup4+p2mup4)*F(189))+p1mup3*p3mup4*F(190)+p2mup3*p3mup4*F(
     -   190)+p3mup3*(p3mup4*F(186)+p4mup4*F(188)+(p1mup4+p2mup4)*F(1
     -   91))-(p1mup4+p2mup4)*p4mup3*F(192)-p3mup4*p4mup3*F(193)-p4mu
     -   p3*p4mup4*F(194))
       Fa8(86) = (p1mup3+p2mup3)*F(122)+p3mup3*F(123)+p4mup3*F(124)
       Fa8(87) = p4mup4*F(174)+(p1mup4+p2mup4)*F(175)+p3mup4*F(176)
       Fa8(88) = p4mup2*s12*F(195)+4*(p1mup2+p2mup2)*s123*F(196)
       Return
       End
