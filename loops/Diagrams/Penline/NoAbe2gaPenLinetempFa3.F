c************** Calling the Fa functions*************************
       subroutine NoAbe2gaPenLineFa3(p1mup2,p1mup4,p1mup5,p2mup2,p2mup
     &   4,p2mup5,p3mup2,p3mup4,p3mup5,p4mup2,p4mup4,p4mup5,p5mup2,p5
     &   mup4,p5mup5,mup2mup4,mup2mup5,mup4mup5,Fa3)
       IMPLICIT NONE
      Complex*16   p1mup2, p1mup4, p1mup5, p2mup2, p2mup4, p2mup5, 
     -          p3mup2, p3mup4, p3mup5, p4mup2, p4mup4, p4mup5, 
     -          p5mup2, p5mup4, p5mup5
       Complex*16   mup2mup4, mup2mup5, mup4mup5
        common/NoAbe2gaPenLineFhlFunctions/F
       COMMON/NoAbe2gaPenLineInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,s12,s
     &   23,s34,s45,s15
       Real*8   p1sq 
       Real*8   p2sq 
       Real*8   p3sq 
       Real*8   p4sq 
       Real*8   p5sq 
       Real*8   s12, s15 
       Real*8   s23 
       Real*8   s34 
       Real*8   s45
      Complex*16 Fa(84),F(426)
      Real*8 P(767) 
        COMMON/NoAbe2gaPenLinePFunctions/P
       Complex*16 Fa3(57:84)
       COMMON/NoAbe2gaPenLineFaFunctions/Fa
       Fa3(57) = p3mup5*F(283)+p1mup5*F(284)+p5mup5*F(285)
       Fa3(58) = p3mup4*F(286)+p5mup4*F(287)
       Fa3(59) = p5mup4*F(290)+p3mup4*F(291)
       Fa3(60) = p3mup5*F(292)+p5mup5*F(293)
       Fa3(61) = p3mup4*F(294)+p5mup4*F(295)
       Fa3(62) = p2mup4*F(297)+p5mup4*F(298)+p3mup4*F(299)
       Fa3(63) = p5mup2*F(300)+p1mup2*F(301)
       Fa3(64) = p3mup5*F(302)+p5mup5*F(303)
       Fa3(65) = p1mup2*F(304)+p3mup2*F(305)+p5mup2*F(306)
       Fa3(66) = p3mup2*F(307)+p1mup2*F(308)+p5mup2*F(309)
       Fa3(67) = p2mup5*F(313)+p3mup5*F(314)+p5mup5*F(315)
       Fa3(68) = p2mup5*F(316)+p1mup5*F(317)+p3mup5*F(318)+p5mup5*F(31
     &   9)
       Fa3(69) = p1mup2*F(320)+p5mup2*F(321)
       Fa3(70) = p1mup2*F(322)+p5mup2*F(323)
       Fa3(71) = p1mup2*F(324)+p5mup2*F(325)
       Fa3(72) = p1mup2*F(326)+p3mup2*F(327)+p5mup2*F(328)
       Fa3(73) = p1mup2*F(329)+p5mup2*F(330)
       Fa3(74) = p1mup2*F(335)+p5mup2*F(336)
       Fa3(75) = p1mup2*p2mup5*p5mup4*F(274)+p3mup2*p3mup5*p5mup4*F(27
     &   7)+p3mup2*p5mup4*p5mup5*F(280)+p3mup5*p5mup2*p5mup4*F(281)+p
     &   1mup2*p1mup4*p3mup5*F(282)+p1mup4*p3mup5*p5mup2*F(288)+p1mup
     &   2*p3mup5*p5mup4*F(289)+p2mup2*p3mup5*p5mup4*F(296)+p1mup2*p3
     &   mup4*p5mup5*F(310)+p1mup2*p2mup4*p5mup5*F(311)+p2mup5*p5mup2
     &   *p5mup4*F(312)+p2mup4*p5mup2*p5mup5*F(331)+p3mup4*p5mup2*p5m
     &   up5*F(332)+p1mup2*p5mup4*p5mup5*F(333)+p2mup2*p5mup4*p5mup5*
     &   F(334)+p5mup2*p5mup4*p5mup5*F(337)+mup2mup4*p1mup5*F(338)+mu
     &   p2mup4*p3mup5*F(339)+mup2mup4*p2mup5*F(340)+mup2mup4*p5mup5*
     &   F(341)+mup2mup5*p1mup4*F(342)+mup2mup5*p3mup4*F(343)+mup2mup
     &   5*p2mup4*F(344)+mup2mup5*p5mup4*F(345)+mup4mup5*p2mup2*F(346
     &   )+mup4mup5*p3mup2*F(347)+mup4mup5*p1mup2*F(348)+mup4mup5*p5m
     &   up2*F(349)+p3mup2*p3mup4*Fa(55)+p3mup4*p3mup5*Fa(56)+p1mup4*
     &   p3mup2*Fa3(57)+p1mup5*p3mup2*Fa3(58)+p2mup2*p2mup5*Fa3(59)+p
     &   2mup4*p3mup2*Fa3(60)+p2mup5*p3mup2*Fa3(61)+p1mup5*p2mup2*Fa3
     &   (62)+p2mup4*p3mup5*Fa3(63)+p2mup2*p3mup4*Fa3(64)+p1mup4*p2mu
     &   p5*Fa3(65)+p1mup5*p2mup4*Fa3(66)+p2mup2*p2mup4*Fa3(67)+p1mup
     &   4*p2mup2*Fa3(68)+p1mup4*p1mup5*Fa3(69)+p1mup5*p5mup4*Fa3(70)
     &   +p1mup5*p3mup4*Fa3(71)+p2mup4*p2mup5*Fa3(72)+p2mup5*p3mup4*F
     &   a3(73)+p1mup4*p5mup5*Fa3(74)
       Fa3(76) = p2mup5*F(350)+p3mup5*F(351)+p1mup5*F(352)+p5mup5*F(35
     &   3)
       Fa3(77) = p1mup5*F(354)+p3mup5*F(355)+p5mup5*F(356)+p2mup5*F(35
     &   7)
       Fa3(78) = p1mup5*p5mup2*F(358)+p1mup2*p1mup5*F(359)+16*p2mup5*p
     &   5mup2*F(360)+p1mup2*p2mup5*F(361)+p1mup5*p2mup2*F(362)+p2mup
     &   2*p2mup5*F(363)+p3mup5*p5mup2*F(364)+p1mup5*p3mup2*F(365)+p1
     &   mup2*p3mup5*F(366)+p2mup5*p3mup2*F(367)+p2mup2*p3mup5*F(368)
     &   +p3mup2*p3mup5*F(369)+p1mup2*p5mup5*F(370)+p2mup2*p5mup5*F(3
     &   71)+p3mup2*p5mup5*F(372)+p5mup2*p5mup5*F(373)+mup2mup5*F(374
     &   )
       Fa3(79) = p1mup2*F(376)+p2mup2*F(377)+p3mup2*F(378)+p5mup2*F(37
     &   9)
       Fa3(80) = p2mup2*F(380)+p3mup2*F(381)+p1mup2*F(382)+p5mup2*F(38
     &   3)
       Fa3(81) = p2mup4*F(384)+p3mup4*F(385)+p1mup4*F(386)+p5mup4*F(38
     &   7)
       Fa3(82) = p1mup4*F(388)+p3mup4*F(389)+p2mup4*F(390)+p5mup4*F(39
     &   1)
       Fa3(83) = p1mup2*p5mup4*F(392)+p1mup2*p1mup4*F(393)+p2mup2*p5mu
     &   p4*F(394)+p1mup2*p2mup4*F(395)+p1mup4*p2mup2*F(396)+p2mup2*p
     &   2mup4*F(397)+p3mup2*p5mup4*F(398)+p1mup4*p3mup2*F(399)+p1mup
     &   2*p3mup4*F(400)+p2mup4*p3mup2*F(401)+p2mup2*p3mup4*F(402)+p3
     &   mup2*p3mup4*F(403)+p1mup4*p5mup2*F(404)+p2mup4*p5mup2*F(405)
     &   +p3mup4*p5mup2*F(406)+p5mup2*p5mup4*F(407)+mup2mup4*F(408)
       Fa3(84) = p2mup5*p5mup4*F(409)+p1mup4*p2mup5*F(410)+p2mup4*p2mu
     &   p5*F(411)+p3mup5*p5mup4*F(412)+p1mup4*p3mup5*F(413)+p2mup5*p
     &   3mup4*F(414)+p2mup4*p3mup5*F(415)+p3mup4*p3mup5*F(416)+p1mup
     &   5*p5mup4*F(417)+p1mup4*p1mup5*F(418)+p1mup5*p2mup4*F(419)+p1
     &   mup5*p3mup4*F(420)+p2mup4*p5mup5*F(421)+p3mup4*p5mup5*F(422)
     &   +p1mup4*p5mup5*F(423)+p5mup4*p5mup5*F(424)+mup4mup5*F(425)
       Return
       End
