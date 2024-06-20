c************** Calling the Fa functions*************************
       subroutine HexBoxCross1Fa12(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,
     &   p3mup4,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa
     &   12)
       IMPLICIT NONE
      Complex*16   p1mup3, p1mup4, p2mup3, 
     -          p2mup4, p3mup3, p3mup4, 
     -          p4mup3, p4mup4, p5mup3, 
     -          p5mup4, p6mup3, p6mup4
       Complex*16   mup3mup4

        common/HexBoxCross1FhlFunctions/F
       COMMON/HexBoxCross1Invariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12
     &   ,s23,s34,s45,s56,s16,s123,s234,s345
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
      Complex*16 Fa(272),F(1461)
      Real*8 P(1825) 
        COMMON/HexBoxCross1PFunctions/P
       Complex*16 Fa12(144:156)
       COMMON/HexBoxCross1FaFunctions/Fa
       Fa12(144) = -(p1mup4*F(275))+(p4mup4+p5mup4)*F(382)
       Fa12(145) = -4*(p1mup3*F(376)-(p4mup3+p5mup3)*F(377)+p6mup3*F(3
     &   78))+p2mup3*F(470)
       Fa12(146) = -(p1mup3*F(465))-p2mup3*F(466)+(p4mup3+p5mup3)*F(46
     &   7)+p6mup3*F(1458)
       Fa12(147) = -8*(p1mup3*F(376)-(p4mup3+p5mup3)*F(377)+p6mup3*F(3
     &   78))+p2mup3*F(480)
       Fa12(148) = (p4mup3+p5mup3)*F(359)-p1mup3*F(361)+p6mup3*F(363)-
     &   p2mup3*F(390)
       Fa12(149) = p1mup4*F(396)+p6mup4*F(397)
       Fa12(150) = -2*(p1mup3*F(376)-(p4mup3+p5mup3)*F(377)+p6mup3*F(3
     &   78))+p2mup3*F(481)
       Fa12(151) = 2*(p1mup3*F(333)+p6mup3*F(337)+p2mup3*F(340))+(p4mu
     &   p3+p5mup3)*F(482)
       Fa12(152) = mup3mup4*F(486)+p1mup3*(8*(p1mup4*F(487)-p5mup4*F(4
     &   88)-p6mup4*F(489)+p4mup4*F(490))+p2mup4*F(491))+p2mup4*p6mup
     &   3*F(509)+8*(p4mup3*p4mup4*F(496)+p4mup3*p5mup4*F(497)+p5mup3
     &   *p5mup4*F(498)+p5mup3*p6mup4*F(499)+p4mup4*p5mup3*F(500)+p4m
     &   up3*p6mup4*F(501)+p1mup4*p6mup3*F(502)+p5mup4*p6mup3*F(503)+
     &   p4mup4*p6mup3*F(504)-p6mup3*p6mup4*F(505)-p2mup3*(p4mup4*F(4
     &   92)+p1mup4*F(493)+p6mup4*F(494)+p5mup4*F(495)-p2mup4*F(506))
     &   +p1mup4*p4mup3*F(507)+p1mup4*p5mup3*F(508)-p2mup4*p4mup3*F(5
     &   10)+p2mup4*p5mup3*F(511))
       Fa12(153) = mup3mup4*F(513)+2*(p5mup4*p6mup3*F(517)+p4mup4*p6mu
     &   p3*F(518)+p6mup3*p6mup4*F(519)+p1mup4*p4mup3*F(520)+p4mup3*p
     &   4mup4*F(521)+p4mup3*p5mup4*F(522)+p2mup4*p4mup3*F(523)-p1mup
     &   3*p4mup4*F(525)+p4mup3*p6mup4*F(526)+p1mup4*p6mup3*F(527)+p5
     &   mup3*(p5mup4*F(514)+p6mup4*F(515)+p4mup4*F(516)+p1mup4*F(524
     &   )+p2mup4*F(528))+p2mup4*p6mup3*F(529)-p2mup3*p4mup4*F(530))
       Fa12(154) = -(p1mup4*F(532))+p5mup4*F(533)-p2mup4*F(534)+p6mup4
     &   *F(535)+p4mup4*F(536)
       Fa12(155) = p4mup4*F(537)+p1mup4*F(538)+p5mup4*F(539)+p6mup4*F(
     &   540)-p2mup4*F(541)
       Fa12(156) = p1mup3*F(542)+p5mup3*F(543)+p2mup3*F(544)+p6mup3*F(
     &   545)+p4mup3*F(546)
       Return
       End
