c************** Calling the Fa functions*************************
       subroutine HexBoxCross1Fa10(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,
     &   p3mup4,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa
     &   10)
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
       Complex*16 Fa10(118:130)
       COMMON/HexBoxCross1FaFunctions/Fa
       Fa10(118) = p1mup4*F(274)+p6mup4*F(448)
       Fa10(119) = -(p1mup4*F(275))+(p4mup4+p5mup4)*F(448)
       Fa10(120) = p2mup3*F(349)+p6mup3*F(351)+(p1mup3+p4mup3+p5mup3)*
     &   F(353)
       Fa10(121) = (p4mup3+p5mup3)*F(335)+p6mup3*F(337)+p1mup3*F(346)-
     &   p2mup3*F(347)
       Fa10(122) = 4*(2*(-(p4mup4*p5mup3*F(342))-p5mup3*p5mup4*F(342)+
     &   p2mup3*p4mup4*F(343)+p2mup3*p5mup4*F(343)+p4mup4*p6mup3*F(34
     &   4)+p5mup4*p6mup3*F(344)-p1mup4*p2mup3*F(360)-p1mup4*p5mup3*F
     &   (364)-p1mup4*p6mup3*F(365)+p1mup3*((p4mup4+p5mup4)*F(341)-p1
     &   mup4*F(362)+p6mup4*F(376))-p5mup3*p6mup4*F(377)-p4mup3*((p4m
     &   up4+p5mup4)*F(342)+p1mup4*F(364)+p6mup4*F(377))+p6mup3*p6mup
     &   4*F(378)+p2mup3*p6mup4*F(393))+mup3mup4*F(450))
       Fa10(123) = -(p1mup3*F(334))-p2mup3*F(336)+(p4mup3+p5mup3)*F(33
     &   8)+p6mup3*F(339)
       Fa10(124) = p1mup3*F(333)+(p4mup3+p5mup3)*F(335)+p6mup3*F(337)+
     &   p2mup3*F(340)
       Fa10(125) = -2*p2mup4*F(275)+(p4mup4+p5mup4)*F(395)
       Fa10(126) = -8*(p4mup4+p5mup4)*((p4mup3+p5mup3)*F(359)-p1mup3*F
     &   (361)+p6mup3*F(363)-p2mup3*F(390))+mup3mup4*F(454)
       Fa10(127) = 2*(p1mup4*s45*F(274)+p6mup4*F(370))+(p4mup4+p5mup4)
     &   *F(457)
       Fa10(128) = (p4mup4+p5mup4)*F(313)+p1mup4*s45*F(453)+p6mup4*F(4
     &   58)
       Fa10(129) = (p4mup4+p5mup4)*F(313)-p1mup4*s45*F(367)+p6mup4*F(3
     &   84)
       Fa10(130) = (p4mup3+p5mup3)*F(461)+2*(p1mup3*F(460)+p2mup3*F(46
     &   2)+p6mup3*F(463))
       Return
       End
