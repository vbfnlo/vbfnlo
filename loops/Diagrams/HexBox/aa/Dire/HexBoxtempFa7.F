c************** Calling the Fa functions*************************
       subroutine HexBoxFa7(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,
     &   p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa7)
       IMPLICIT NONE
      Complex*16   p1mup3, p1mup4, p2mup3, 
     -          p2mup4, p3mup3, p3mup4, 
     -          p4mup3, p4mup4, p5mup3, 
     -          p5mup4, p6mup3, p6mup4
       Complex*16   mup3mup4

        common/HexBoxFhlFunctions/F
       COMMON/HexBoxInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12,s23,s
     &   34,s45,s56,s16,s123,s234,s345
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
      Complex*16 Fa(259),F(1314)
      Real*8 P(1362) 
        COMMON/HexBoxPFunctions/P
       Complex*16 Fa7(73:84)
       COMMON/HexBoxFaFunctions/Fa
       Fa7(73) = mup3mup4*F(302)+16*(p4mup4+p5mup4)*((p4mup3+p5mup3)*F
     &   (295)+p1mup3*F(303)+p2mup3*F(304)+p6mup3*F(305))
       Fa7(74) = mup3mup4*F(308)+2*(-(p4mup4*p5mup3*F(279))-p5mup3*p5m
     &   up4*F(279)+p2mup3*p4mup4*F(282)+p2mup3*p5mup4*F(282)+p4mup4*
     &   p6mup3*F(284)+p5mup4*p6mup3*F(284)-p5mup3*p6mup4*F(310)+p1mu
     &   p3*((p4mup4+p5mup4)*F(273)+p6mup4*F(309)-2*p1mup4*F(312))-p2
     &   mup3*p6mup4*F(313)+p6mup3*p6mup4*F(314)-p4mup3*((p4mup4+p5mu
     &   p4)*F(279)+p6mup4*F(310)+2*p1mup4*F(315))-2*(p1mup4*p2mup3*F
     &   (311)+p1mup4*p5mup3*F(315)+p1mup4*p6mup3*F(316)))
       Fa7(75) = p1mup4*F(318)-p2mup4*F(319)-p6mup4*F(320)+(p4mup4+p5m
     &   up4)*F(321)
       Fa7(76) = p1mup3*F(322)+p2mup3*F(323)-p6mup3*F(324)-(p4mup3+p5m
     &   up3)*F(325)
       Fa7(77) = mup3mup4*F(326)-16*(p2mup3*((p4mup4+p5mup4)*F(273)+p6
     &   mup4*F(309)-2*p2mup4*F(311))-2*p2mup4*((p4mup3+p5mup3)*F(315
     &   )+p6mup3*F(316))-p4mup4*p5mup3*F(328)-p5mup3*p5mup4*F(328)-p
     &   4mup3*(p4mup4+p5mup4)*F(328)-p4mup3*p6mup4*F(329)-p5mup3*p6m
     &   up4*F(329)-p1mup3*(2*p2mup4*F(312)-(p4mup4+p5mup4)*F(327)+p6
     &   mup4*F(330))+p4mup4*p6mup3*F(331)+p5mup4*p6mup3*F(331)+p6mup
     &   3*p6mup4*F(332))
       Fa7(78) = mup3mup4*F(333)+2*(p2mup4*((p4mup3+p5mup3)*F(279)-p2m
     &   up3*F(282)-p6mup3*F(284))+p1mup4*(p2mup3*F(273)+p1mup3*F(327
     &   )-(p4mup3+p5mup3)*F(328)+p6mup3*F(331))+p2mup3*p6mup4*F(334)
     &   -p1mup3*(p2mup4*F(273)+p6mup4*F(335)))
       Fa7(79) = mup3mup4*F(336)+2*(p4mup4*p5mup3*F(280)+p5mup3*p5mup4
     &   *F(280)+p4mup4*p6mup3*F(281)+p5mup4*p6mup3*F(281)-p1mup4*p5m
     &   up3*F(329)+p1mup4*p6mup3*F(332)+p1mup3*((p4mup4+p5mup4)*F(27
     &   5)-p1mup4*F(330)-2*p2mup4*F(337))+p4mup3*((p4mup4+p5mup4)*F(
     &   280)-p1mup4*F(329)+2*p2mup4*F(338))+2*(p2mup4*p5mup3*F(338)-
     &   p2mup4*p6mup3*F(339))-p2mup3*((p4mup4+p5mup4)*F(278)-p1mup4*
     &   F(309)-2*p2mup4*F(340)))
       Fa7(80) = mup3mup4*F(262)-2*(p4mup4+p5mup4)*(p2mup3*F(334)-p1mu
     &   p3*F(335))
       Fa7(81) = p1mup4*F(267)+p6mup4*F(343)
       Fa7(82) = p1mup3*F(273)-(p4mup3+p5mup3)*F(279)+p2mup3*F(282)+p6
     &   mup3*F(284)
       Fa7(83) = (p4mup4+p5mup4)*F(344)+p1mup4*F(345)-8*p2mup4*F(346)
       Fa7(84) = p1mup3*F(349)+p2mup3*F(350)+p6mup3*F(351)+(p4mup3+p5m
     &   up3)*F(352)
       Return
       End
