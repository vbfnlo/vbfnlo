c************** Calling the Fa functions*************************
       subroutine NoAbe2gaHexLineFa5(p1mup2,p1mup3,p1mup5,p1mup6,p2mup
     -   2,p2mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4
     -   mup3,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3
     -   ,p6mup5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,
     -   mup5mup6,Fa5)
       IMPLICIT NONE
      Complex*16   p1mup2, p1mup3, p1mup5, p1mup6, p2mup2, p2mup3, 
     -          p2mup5, p2mup6, p3mup2, p3mup3, p3mup5, p3mup6, 
     -          p4mup2, p4mup3, p4mup5, p4mup6, p5mup2, p5mup3, 
     -          p5mup5, p5mup6, p6mup2, p6mup3, p6mup5, p6mup6
       Complex*16   mup2mup3, mup2mup5, mup2mup6, mup3mup5, mup3mup6, 
     -          mup5mup6
        common/NoAbe2gaHexLineFhlFunctions/F
       COMMON/NoAbe2gaHexLineInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,
     -   s12,s23,s34,s45,s56,s16,s123,s234,s345
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
      Complex*16 Fa(872),F(9706)
      Real*8 P(9687) 
        COMMON/NoAbe2gaHexLinePFunctions/P
       Complex*16 Fa5(97:120)
       COMMON/NoAbe2gaHexLineFaFunctions/Fa
       Fa5(97) = p1mup5*p2mup6*p6mup2*F(103)+(2*p2mup5*((p1mup2+p2mup2
     -   )*p3mup6+(p3mup6+p4mup6)*p5mup2)+mup2mup5*(p3mup6+p4mup6)*s2
     -   34)*F(109)-4*p1mup5*(p3mup2+p4mup2)*(p3mup6+p4mup6)*F(331)-2
     -   *(p3mup2+p4mup2)*(p3mup6+p4mup6)*(p5mup5*F(332)+p2mup5*F(333
     -   )+2*p6mup5*F(334))
       Fa5(98) = 2*p2mup5*(p3mup6+p4mup6)*p6mup2*F(109)-p2mup6*p6mup2*
     -   p6mup5*F(294)-4*p4mup2*(p3mup6*p4mup5+(p3mup5+p4mup5)*p4mup6
     -   )*F(9575)
       Fa5(99) = p3mup5*p3mup6*p4mup2+p3mup2*(p3mup5+p4mup5)*(p3mup6+p
     -   4mup6)
       Fa5(100) = 8*(p3mup2*p4mup3*p5mup6*F(100)+p4mup2*p4mup3*p5mup6*
     -   F(100)+p4mup3*p5mup2*p5mup6*F(100)+p4mup3*p5mup6*p6mup2*F(10
     -   0)+p1mup6*p3mup2*p4mup3*F(104)+p1mup6*p4mup2*p4mup3*F(104)+p
     -   1mup6*p4mup3*p5mup2*F(104)+p1mup6*p4mup3*p6mup2*F(104)+p2mup
     -   6*p3mup2*p4mup3*F(109)+p3mup2*p3mup6*p4mup3*F(109)+p2mup6*p4
     -   mup2*p4mup3*F(109)+p3mup6*p4mup2*p4mup3*F(109)+p3mup2*p4mup3
     -   *p4mup6*F(109)+p4mup2*p4mup3*p4mup6*F(109)+p2mup6*p4mup3*p5m
     -   up2*F(109)+p3mup6*p4mup3*p5mup2*F(109)+p4mup3*p4mup6*p5mup2*
     -   F(109)+p2mup6*p4mup3*p6mup2*F(109)+p3mup6*p4mup3*p6mup2*F(10
     -   9)+p4mup3*p4mup6*p6mup2*F(109))+p3mup2*p4mup3*p6mup6*F(323)+
     -   p4mup2*p4mup3*p6mup6*F(323)+p4mup3*p5mup2*p6mup6*F(323)+p4mu
     -   p3*p6mup2*p6mup6*F(323)+p1mup2*p4mup3*(8*(p5mup6*F(100)+p1mu
     -   p6*F(104)+(p2mup6+p3mup6+p4mup6)*F(109))+p6mup6*F(323))+p2mu
     -   p2*p4mup3*(8*(p5mup6*F(100)+p1mup6*F(104)+(p2mup6+p3mup6+p4m
     -   up6)*F(109))+p6mup6*F(323))+4*mup3mup6*(p3mup2+p4mup2)*F(335
     -   )
       Fa5(101) = -(p3mup2*p4mup3*p5mup5*F(99))-p4mup2*p4mup3*p5mup5*F
     -   (99)-p4mup3*p5mup2*p5mup5*F(99)-p4mup3*p5mup5*p6mup2*F(99)+p
     -   1mup5*p3mup2*p4mup3*F(103)+p1mup5*p4mup2*p4mup3*F(103)+p1mup
     -   5*p4mup3*p5mup2*F(103)+p1mup5*p4mup3*p6mup2*F(103)+p2mup5*p3
     -   mup2*p4mup3*F(107)+p3mup2*p3mup5*p4mup3*F(107)+p2mup5*p4mup2
     -   *p4mup3*F(107)+p3mup5*p4mup2*p4mup3*F(107)+p3mup2*p4mup3*p4m
     -   up5*F(107)+p4mup2*p4mup3*p4mup5*F(107)+p2mup5*p4mup3*p5mup2*
     -   F(107)+p3mup5*p4mup3*p5mup2*F(107)+p4mup3*p4mup5*p5mup2*F(10
     -   7)+p2mup5*p4mup3*p6mup2*F(107)+p3mup5*p4mup3*p6mup2*F(107)+p
     -   4mup3*p4mup5*p6mup2*F(107)-p3mup2*p4mup3*p6mup5*F(294)-p4mup
     -   2*p4mup3*p6mup5*F(294)-p4mup3*p5mup2*p6mup5*F(294)-p4mup3*p6
     -   mup2*p6mup5*F(294)+p1mup2*p4mup3*(-(p5mup5*F(99))+p1mup5*F(1
     -   03)+(p2mup5+p3mup5+p4mup5)*F(107)-p6mup5*F(294))+p2mup2*p4mu
     -   p3*(-(p5mup5*F(99))+p1mup5*F(103)+(p2mup5+p3mup5+p4mup5)*F(1
     -   07)-p6mup5*F(294))-mup3mup5*(p3mup2+p4mup2)*F(335)
       Fa5(102) = p3mup2+p4mup2
       Fa5(103) = 2*(p5mup5*(2*p5mup6*F(316)+p6mup6*F(317))-p6mup5*p6m
     -   up6*F(318)-p1mup6*(p5mup5*F(307)+4*(p6mup5*F(311)+p1mup5*F(3
     -   20))))+mup5mup6*F(336)+p1mup5*p6mup6*F(337)+p2mup5*p6mup6*F(
     -   338)+p3mup5*p6mup6*F(338)+p4mup5*p6mup6*F(338)+(p2mup6+p3mup
     -   6+p4mup6)*p5mup5*F(339)+8*(p1mup5*p5mup6*F(309)+p5mup6*p6mup
     -   5*F(319)+p1mup6*(p2mup5+p3mup5+p4mup5)*F(327)+p2mup5*p5mup6*
     -   F(328)+p3mup5*p5mup6*F(328)+p4mup5*p5mup6*F(328)-p1mup5*(p2m
     -   up6+p3mup6+p4mup6)*F(331)-p2mup6*p6mup5*F(334)-p3mup6*p6mup5
     -   *F(334)-p4mup6*p6mup5*F(334)-(p2mup5+p3mup5+p4mup5)*(p2mup6+
     -   p3mup6+p4mup6)*F(9575))
       Fa5(104) = p1mup5*F(272)+p5mup5*F(273)+(p2mup5+p3mup5+p4mup5)*F
     -   (275)+p6mup5*F(340)
       Fa5(105) = p1mup5*F(341)+p5mup5*F(342)+p6mup5*F(343)+(p2mup5+p3
     -   mup5+p4mup5)*F(344)
       Fa5(106) = p5mup6*F(181)+(p2mup6+p3mup6+p4mup6)*F(185)-p1mup6*F
     -   (345)-p6mup6*F(346)
       Fa5(107) = p1mup6*F(347)+p5mup6*F(348)+p6mup6*F(349)+(p2mup6+p3
     -   mup6+p4mup6)*F(350)
       Fa5(108) = mup3mup6*(p2mup5*F(352)+s234*(2*(p5mup5*F(99)-p1mup5
     -   *F(103)-(p3mup5+p4mup5)*F(107))+p6mup5*F(353)))+mup3mup5*(-4
     -   *p2mup6*F(325)+s234*(-4*p5mup6*F(100)+p6mup6*F(354)+(p3mup6+
     -   p4mup6)*F(355)+p1mup6*F(356)))
       Fa5(109) = 4*mup3mup6*((p1mup2+p2mup2+p3mup2+p4mup2+p5mup2+p6mu
     -   p2)*(p5mup5*F(99)-p1mup5*F(103)-(p2mup5+p3mup5+p4mup5)*F(107
     -   )+p6mup5*F(294))-mup2mup5*F(335))+mup3mup5*(4*mup2mup6*F(335
     -   )-(p1mup2+p2mup2+p3mup2+p4mup2+p5mup2+p6mup2)*(8*p5mup6*F(10
     -   0)-p6mup6*F(357)-(p2mup6+p3mup6+p4mup6)*F(358)-p1mup6*F(359)
     -   ))
       Fa5(110) = 2*p2mup6*F(325)+s234*(2*p5mup6*F(100)+p6mup6*F(102)+
     -   (p3mup6+p4mup6)*F(360)+p1mup6*F(361))
       Fa5(111) = (p1mup2+p2mup2+p3mup2+p4mup2+p5mup2+p6mup2)*(4*(p5mu
     -   p6*F(100)+p1mup6*F(104)+(p2mup6+p3mup6+p4mup6)*F(109))+p6mup
     -   6*F(326))-2*mup2mup6*F(335)
       Fa5(112) = p2mup5*F(362)+s234*(p6mup5*F(101)+p1mup5*F(103)+(p3m
     -   up5+p4mup5)*F(107)+p5mup5*F(363))
       Fa5(113) = -((p1mup2+p2mup2+p3mup2+p4mup2+p5mup2+p6mup2)*(2*(p5
     -   mup5*F(99)-p1mup5*F(103)-(p2mup5+p3mup5+p4mup5)*F(107))-p6mu
     -   p5*F(324)))+2*mup2mup5*F(335)
       Fa5(114) = mup3mup6*((p1mup5+p2mup5)*F(364)+p5mup5*F(366)+2*(p6
     -   mup5*F(368)-(p3mup5+p4mup5)*F(370)))+mup3mup5*((p1mup6+p2mup
     -   6)*F(365)-p5mup6*F(367)+p6mup6*F(369)+(p3mup6+p4mup6)*F(371)
     -   )
       Fa5(115) = mup3mup6*(4*(p1mup2+p2mup2)*F(372)+s12*(p5mup2*F(373
     -   )+4*p6mup2*F(374)))
       Fa5(116) = mup5mup6*(-(p2mup2*F(375))+p5mup2*s12*F(376))
       Fa5(117) = p2mup2*p2mup6*(p2mup5*F(377)+4*(p5mup5*F(378)-p6mup5
     -   *F(379)))+p1mup5*(p1mup6*p2mup2*F(377)+p2mup2*p2mup6*F(377)+
     -   4*(p2mup2*p5mup6*F(380)+p5mup2*p5mup6*s12*F(381)))+4*p5mup2*
     -   p5mup6*s12*(p2mup5*F(381)+p5mup5*F(382)+p6mup5*F(383))
       Fa5(118) = mup5mup6*p6mup2*s12*F(388)+p2mup2*(4*(p2mup5*(p5mup6
     -   *F(380)-p6mup6*F(384))+p5mup6*p6mup5*F(386)-p6mup6*(p1mup5*F
     -   (384)+p5mup5*F(385)+p6mup5*F(387)))+p5mup5*p5mup6*F(389))+p1
     -   mup6*(p2mup2*(p2mup5*F(377)+4*(p5mup5*F(378)-p6mup5*F(379)))
     -   +4*p5mup2*s12*((p1mup5+p2mup5)*F(392)+p5mup5*F(393)+p6mup5*F
     -   (394)))+4*p5mup2*s12*(-2*(p1mup5*p6mup6*F(390)+p2mup5*p6mup6
     -   *F(390)+p6mup5*p6mup6*F(391))+p1mup5*p2mup6*F(392)+p2mup5*p2
     -   mup6*F(392)+p2mup6*p5mup5*F(393)+p2mup6*p6mup5*F(394)+p5mup5
     -   *p6mup6*F(395))
       Fa5(119) = 4*p2mup2*p4mup5*(p5mup6*F(396)+p6mup6*F(397)+4*p4mup
     -   6*F(398))+p1mup2*(p1mup6*p2mup5*F(377)+p2mup5*p2mup6*F(377)+
     -   p1mup5*((p1mup6+p2mup6)*F(377)+4*(p5mup6*F(380)-p6mup6*F(384
     -   )))+p5mup5*p5mup6*F(389)+4*(p1mup6*p5mup5*F(378)+p2mup6*p5mu
     -   p5*F(378)-p1mup6*p6mup5*F(379)-p2mup6*p6mup5*F(379)+p2mup5*p
     -   5mup6*F(380)-p2mup5*p6mup6*F(384)-p5mup5*p6mup6*F(385)+p5mup
     -   6*p6mup5*F(386)-p6mup5*p6mup6*F(387)+p4mup5*p5mup6*F(396)+p4
     -   mup5*p6mup6*F(397))+16*p4mup5*p4mup6*F(398))
       Fa5(120) = p1mup2*p3mup5*(p5mup6*F(396)+p6mup6*F(397))+p2mup2*p
     -   3mup5*(p5mup6*F(396)+p6mup6*F(397))+p1mup2*(p1mup6+p2mup6)*(
     -   p3mup5+p4mup5)*F(399)+p2mup2*(p1mup6+p2mup6)*(p3mup5+p4mup5)
     -   *F(399)+s12*(p2mup5*p6mup2*(p5mup6*F(400)+p6mup6*F(401)+p1mu
     -   p6*F(402))+2*(p3mup5+p4mup5)*p5mup2*(-(p6mup6*F(403))+(p1mup
     -   6+p2mup6)*F(404)+p5mup6*F(405)))
       Return
       End
