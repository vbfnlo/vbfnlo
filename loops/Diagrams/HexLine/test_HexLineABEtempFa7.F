c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineAbe7(p1mup2,p1mup3,p1mup4,p1mu
     -   p5,p2mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p
     -   4mup2,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup
     -   2,p6mup3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,m
     -   up3mup5,mup4mup5,Fa7)
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
       Complex*16 Fa7(67:77)
       COMMON/test_FaFunctionsAbe/Fa
       Fa7(67) = p1mup2*p1mup3*F(120)+p1mup3*p2mup2*F(120)+p1mup2*p2mu
     -   p3*F(120)+p2mup2*p2mup3*F(120)-p1mup2*p4mup3*F(121)-p2mup2*p
     -   4mup3*F(121)+p1mup3*p4mup2*s12*F(122)+p2mup3*p4mup2*s12*F(12
     -   2)+p3mup3*p4mup2*s12*F(123)+p4mup2*p4mup3*s12*F(124)-(p1mup2
     -   +p2mup2)*p3mup3*F(125)
       Fa7(68) = p1mup3*p3mup2*F(120)+p2mup3*p3mup2*F(120)-p3mup2*p4mu
     -   p3*F(121)-p3mup2*p3mup3*F(125)-p1mup3*p4mup2*F(126)-p2mup3*p
     -   4mup2*F(126)-p3mup3*p4mup2*F(127)+p4mup2*p4mup3*F(128)
       Fa7(69) = (p1mup3+p2mup3)*F(120)-p4mup3*F(121)-p3mup3*F(125)
       Fa7(70) = -(mup3mup4*p4mup2*s12*F(133))-mup3mup4*p3mup2*s12*F(1
     -   36)-8*(mup3mup4*p2mup2*F(131)-p1mup2*(-(mup3mup4*F(131))+2*(
     -   p3mup3*(-(p3mup4*F(130))+(p1mup4+p2mup4)*F(134)-p4mup4*F(135
     -   ))-p4mup3*(p3mup4*F(132)+p4mup4*F(137)+(p1mup4+p2mup4)*F(138
     -   ))+p1mup3*(p3mup4*F(129)+(p1mup4+p2mup4)*F(139)+p4mup4*F(140
     -   ))+p2mup3*(p3mup4*F(129)+(p1mup4+p2mup4)*F(139)+p4mup4*F(140
     -   )))))+16*(p2mup2*p2mup3*p3mup4*F(129)-p2mup2*p3mup3*p3mup4*F
     -   (130)-p2mup2*p3mup4*p4mup3*F(132)+p1mup4*p2mup2*p3mup3*F(134
     -   )+p2mup2*p2mup4*p3mup3*F(134)-p2mup2*p3mup3*p4mup4*F(135)-p2
     -   mup2*p4mup3*p4mup4*F(137)-p1mup4*p2mup2*p4mup3*F(138)-p2mup2
     -   *p2mup4*p4mup3*F(138)+p1mup4*p2mup2*p2mup3*F(139)+p2mup2*p2m
     -   up3*p2mup4*F(139)+p2mup2*p2mup3*p4mup4*F(140)+p1mup4*p3mup2*
     -   p3mup3*s12*F(141)+p2mup4*p3mup2*p3mup3*s12*F(141)+p2mup3*p3m
     -   up2*p3mup4*s12*F(141)+p3mup2*p3mup3*p3mup4*s12*F(142)+p3mup2
     -   *p3mup3*p4mup4*s12*F(143)+p3mup3*p3mup4*p4mup2*s12*F(144)+p2
     -   mup3*p4mup2*p4mup4*s12*F(145)+p3mup3*p4mup2*p4mup4*s12*F(146
     -   )+p1mup4*p2mup3*p3mup2*s12*F(147)+p2mup3*p2mup4*p3mup2*s12*F
     -   (147)+p1mup4*p2mup3*p4mup2*s12*F(148)+p2mup3*p2mup4*p4mup2*s
     -   12*F(148)+p2mup3*p3mup2*p4mup4*s12*F(149)+p2mup3*p3mup4*p4mu
     -   p2*s12*F(150)+p1mup3*(p2mup2*(p3mup4*F(129)+(p1mup4+p2mup4)*
     -   F(139)+p4mup4*F(140))+s12*(p3mup2*(p3mup4*F(141)+(p1mup4+p2m
     -   up4)*F(147)+p4mup4*F(149))+p4mup2*(p4mup4*F(145)+(p1mup4+p2m
     -   up4)*F(148)+p3mup4*F(150))))+p1mup4*p3mup3*p4mup2*s12*F(151)
     -   +p2mup4*p3mup3*p4mup2*s12*F(151)+p1mup4*p3mup2*p4mup3*s12*F(
     -   152)+p2mup4*p3mup2*p4mup3*s12*F(152)+p3mup2*p3mup4*p4mup3*s1
     -   2*F(153)+p1mup4*p4mup2*p4mup3*s12*F(154)+p2mup4*p4mup2*p4mup
     -   3*s12*F(154)+p3mup2*p4mup3*p4mup4*s12*F(155)+p3mup4*p4mup2*p
     -   4mup3*s12*F(156)+p4mup2*p4mup3*p4mup4*s12*F(157))
       Fa7(71) = p1mup3*p3mup2*s12*F(122)+p2mup3*p3mup2*s12*F(122)+p3m
     -   up2*p3mup3*s12*F(123)+p3mup2*p4mup3*s12*F(124)+p1mup3*(p1mup
     -   2+p2mup2)*F(126)+(p1mup2+p2mup2)*p2mup3*F(126)+p1mup2*p3mup3
     -   *F(127)+p2mup2*p3mup3*F(127)+(p1mup2+p2mup2)*p4mup3*F(158)
       Fa7(72) = p3mup4*F(159)+4*(p1mup4+p2mup4)*F(160)+p4mup4*F(161)
       Fa7(73) = 8*p4mup2*(p4mup4*F(163)+(p1mup4+p2mup4)*F(164)+p3mup4
     -   *F(165))+p3mup2*((p1mup4+p2mup4)*F(166)+8*(p3mup4*F(162)+p4m
     -   up4*F(167)))
       Fa7(74) = mup3mup4*F(168)+8*(-(p4mup3*(p3mup4*F(132)+p4mup4*F(1
     -   37)))+p3mup3*(-(p3mup4*F(130))+(p1mup4+p2mup4)*F(134)+p4mup4
     -   *F(169))+(p1mup4+p2mup4)*p4mup3*F(170)+p1mup3*(p3mup4*F(129)
     -   +(p1mup4+p2mup4)*F(139)+p4mup4*F(171))+p2mup3*(p3mup4*F(129)
     -   +(p1mup4+p2mup4)*F(139)+p4mup4*F(171)))
       Fa7(75) = p1mup2*(p4mup4*F(167)+p3mup4*F(172)+(p1mup4+p2mup4)*F
     -   (173))+p2mup2*(p4mup4*F(167)+p3mup4*F(172)+(p1mup4+p2mup4)*F
     -   (173))-p4mup2*s12*(p4mup4*F(174)+(p1mup4+p2mup4)*F(175)+p3mu
     -   p4*F(176))
       Fa7(76) = (p1mup3+p2mup3)*F(126)+p3mup3*F(127)+p4mup3*F(158)
       Fa7(77) = 4*(p4mup4*F(167)+p3mup4*F(172))+(p1mup4+p2mup4)*F(182
     -   )
       Return
       End
