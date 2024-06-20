c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineAbe39(p1mup2,p1mup3,p1mup4,p1m
     -   up5,p2mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,
     -   p4mup2,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mu
     -   p2,p6mup3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,
     -   mup3mup5,mup4mup5,Fa39)
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
       Complex*16 Fa39(419:429)
       COMMON/test_FaFunctionsAbe/Fa
       Fa39(419) = -(mup2mup4*F(5911))+2*(-3*(p4mup2*p4mup4*F(4647)-p3
     -   mup2*p3mup4*F(5587)+p3mup4*p4mup2*F(5588)+p3mup2*p4mup4*F(55
     -   90))-p1mup2*p6mup4*F(5914)+p1mup2*p3mup4*F(5916)-p3mup2*p6mu
     -   p4*F(5917)+p1mup4*(p1mup2*F(5912)+p6mup2*F(5918))+p2mup4*(-3
     -   *(p4mup2*F(4614)+p3mup2*F(5576))+p1mup2*F(5913)+p6mup2*F(591
     -   9))+p3mup4*p6mup2*F(5920)-p4mup2*p6mup4*F(5921)+p4mup4*p6mup
     -   2*F(5922)+p1mup4*p3mup2*F(5923)-p1mup2*p4mup4*F(5924)+p2mup2
     -   *(3*(p2mup4*F(4611)-p4mup4*F(4644)-p3mup4*F(5575))-p6mup4*F(
     -   5915)+p1mup4*F(5925))+p1mup4*p4mup2*F(5926)+p6mup2*p6mup4*F(
     -   5927))
       Fa39(420) = F(225)*Fa(102)+2*(Fa(103)+Fa(104)+Fa(105)+Fa(106)+F
     -   a(107))
       Fa39(421) = F(225)*Fa(109)+2*(Fa(110)+Fa(111)+Fa(112)+Fa(113)+F
     -   a(114)+Fa(115)+Fa(116))
       Fa39(422) = Fa(162)+Fa(163)+Fa(164)-8*s12*Fa(165)
       Fa39(423) = Fa(179)+16*Fa(180)+Fa(181)+Fa(182)+32*(Fa(183)+Fa(1
     -   84)+Fa(185)+Fa(186)+Fa(187)-Fa(188))+Fa(189)
       Fa39(424) = 16*(Fa(190)+Fa(191)+Fa(192)-Fa(193)+2*(Fa(194)+Fa(1
     -   95)+Fa(196)+Fa(197)+Fa(198)+Fa(199)))-Fa(200)
       Fa39(425) = Fa(201)+16*Fa(202)+Fa(203)+Fa(204)-32*(Fa(205)-Fa(2
     -   06)+Fa(207)+Fa(208)-Fa(209)-Fa(210))+Fa(211)
       Fa39(426) = Fa(212)+Fa(213)+Fa(214)-2*(Fa(215)+Fa(216))-Fa(217)
       Fa39(427) = 8*(Fa(222)+Fa(223)+Fa(224)-2*(Fa(225)+Fa(226)))+Fa(
     -   227)
       Fa39(428) = Fa(230)+Fa(231)+Fa(232)-16*(Fa(233)+Fa(234))+Fa(235
     -   )
       Fa39(429) = Fa(263)+Fa(264)+Fa(265)+Fa(266)
       Return
       End
