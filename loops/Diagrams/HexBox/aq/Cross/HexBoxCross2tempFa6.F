c************** Calling the Fa functions*************************
       subroutine HexBoxCross2Fa6(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p
     &   3mup6,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa6
     &   )
       IMPLICIT NONE
      Complex*16   p1mup3, p1mup6, p2mup3, p2mup6, p3mup3, p3mup6, 
     -          p4mup3, p4mup6, p5mup3, p5mup6, p6mup3, p6mup6
       Complex*16   mup3mup6
        common/HexBoxCross2FhlFunctions/F
       COMMON/HexBoxCross2Invariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12
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
      Complex*16 Fa(392),F(1853)
      Real*8 P(1355) 
        COMMON/HexBoxCross2PFunctions/P
       Complex*16 Fa6(96:114)
       COMMON/HexBoxCross2FaFunctions/Fa
       Fa6(96) = p1mup6*F(387)+p6mup6*F(388)+p2mup6*F(389)+p5mup6*F(39
     &   0)
       Fa6(97) = 4*mup3mup6*F(393)-8*(p6mup3*(p6mup6*F(395)+p5mup6*F(3
     &   96)-p1mup6*F(397)-p2mup6*F(398))+p2mup3*(-(p6mup6*F(394))-p5
     &   mup6*F(399)+p2mup6*F(400)+p1mup6*F(401)))
       Fa6(98) = mup3mup6*F(402)+2*(p1mup3*(p6mup6*F(403)-p1mup6*F(404
     &   )-p2mup6*F(405)-p5mup6*F(407))+p6mup3*(p6mup6*F(317)+p5mup6*
     &   F(318)+p1mup6*F(409)+p2mup6*F(410))+p2mup3*(p6mup6*F(406)+p2
     &   mup6*F(408)+p5mup6*F(411)+p1mup6*F(412)))
       Fa6(99) = p1mup6*F(415)+p6mup6*F(416)+p2mup6*F(417)+p5mup6*F(41
     &   8)
       Fa6(100) = p1mup6*F(419)+p6mup6*F(420)+p2mup6*F(421)+p5mup6*F(4
     &   22)
       Fa6(101) = p1mup6*F(423)+p6mup6*F(424)+p2mup6*F(425)+p5mup6*F(4
     &   26)
       Fa6(102) = p1mup3*F(430)+p5mup3*F(431)+p6mup3*F(432)
       Fa6(103) = p1mup3*F(433)+p5mup3*F(434)+p6mup3*F(435)+4*p2mup3*F
     &   (436)
       Fa6(104) = p5mup6*F(446)+p6mup6*F(447)+p1mup6*F(448)+p2mup6*F(4
     &   49)
       Fa6(105) = p2mup3*F(450)-p6mup3*F(451)
       Fa6(106) = mup3mup6*F(453)-2*(-(p5mup3*p5mup6*F(407))+p6mup6*(p
     &   5mup3*F(403)+p2mup3*F(452))-p6mup3*p6mup6*F(454)-p5mup6*p6mu
     &   p3*F(455)+p2mup3*p5mup6*F(456)-p2mup6*(p5mup3*F(405)+p2mup3*
     &   F(457)+p6mup3*F(458))-p1mup6*(p6mup3*F(325)+p5mup3*F(404)+p2
     &   mup3*F(459)))
       Fa6(107) = -(p6mup6*F(461))+p5mup6*F(462)+p2mup6*F(463)+p1mup6*
     &   F(464)
       Fa6(108) = p5mup6*F(279)+p6mup6*F(280)-p1mup6*F(281)-p2mup6*F(2
     &   82)
       Fa6(109) = p1mup6*F(309)+p5mup6*F(310)+p6mup6*F(311)-p2mup6*F(3
     &   19)
       Fa6(110) = p1mup6*F(321)+p5mup6*F(322)+p6mup6*F(323)+2*p2mup6*F
     &   (327)
       Fa6(111) = p5mup3*F(465)-p2mup3*F(466)+p6mup3*F(1848)
       Fa6(112) = p1mup3*F(465)+p2mup3*F(466)+p6mup3*F(467)
       Fa6(113) = -(p2mup3*F(466))+p1mup3*F(468)+p5mup3*F(469)+p6mup3*
     &   F(470)
       Fa6(114) = -(p2mup3*F(450))+p6mup3*F(472)
       Return
       End
