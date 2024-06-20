c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineNoAbe1g15(p1mup2,p1mup3,p1mup4
     -   ,p1mup6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3m
     -   up6,p4mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,
     -   p6mup2,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3m
     -   up4,mup3mup6,mup4mup6,Fa15)
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
       Complex*16 Fa15(197:210)
       COMMON/test_FaFunctionsNoAbe1g/Fa
       Fa15(197) = -2*p5mup6*F(234)+(p1mup6+p2mup6+p3mup6)*F(238)+p4mu
     -   p6*F(239)
       Fa15(198) = p4mup4*F(263)-2*((p1mup4+p2mup4+p3mup4)*F(264)-p5mu
     -   p4*F(271))
       Fa15(199) = (p1mup2+p2mup2)*(p5mup4*F(273)+p4mup4*F(274)+(p1mup
     -   4+p2mup4+p3mup4)*F(275))
       Fa15(200) = 2*(p5mup4*F(273)+p4mup4*F(274))+(p1mup4+p2mup4+p3mu
     -   p4)*F(280)
       Fa15(201) = p4mup6*F(235)+(p1mup6+p2mup6+p3mup6)*F(236)+p5mup6*
     -   F(237)
       Fa15(202) = -(mup4mup6*s45*F(282))+p5mup4*(4*(p3mup6*F(281)+(p4
     -   mup6+p5mup6)*F(283)-p1mup6*F(284))+p2mup6*F(285))
       Fa15(203) = -(mup3mup6*s45*F(289))+p4mup3*p4mup6*F(290)-p4mup6*
     -   p5mup3*F(291)-p1mup3*p4mup6*F(292)+p2mup6*p4mup3*F(294)+p4mu
     -   p3*p5mup6*F(295)+p2mup3*(-(p4mup6*F(288))+4*s45*(p1mup6*F(29
     -   7)+p3mup6*F(301)+p2mup6*F(304)+p5mup6*F(308)))+p3mup3*(-(p4m
     -   up6*F(286))+4*s45*(p1mup6*F(298)+p3mup6*F(299)+p2mup6*F(306)
     -   +p5mup6*F(309)))+4*(p3mup6*p4mup3*F(287)-p1mup6*p4mup3*F(293
     -   )+p1mup3*p1mup6*s45*F(296)+p1mup3*p3mup6*s45*F(300)+p3mup6*p
     -   5mup3*s45*F(302)+p1mup3*p2mup6*s45*F(303)+p1mup6*p5mup3*s45*
     -   F(305)+p1mup3*p5mup6*s45*F(307)+p2mup6*p5mup3*s45*F(310)+p5m
     -   up3*p5mup6*s45*F(311))
       Fa15(204) = p4mup6*F(312)-2*s45*(p3mup6*F(313)+p2mup6*F(314)+p5
     -   mup6*F(315)+p1mup6*F(316))
       Fa15(205) = p3mup2*p4mup6*F(317)+2*mup2mup6*s45*F(319)+p1mup2*p
     -   4mup6*F(320)+p2mup2*p4mup6*F(321)+p4mup6*p5mup2*F(323)+p4mup
     -   2*(p3mup6*F(318)+p4mup6*F(322)+p2mup6*F(324)+4*(p5mup6*F(283
     -   )-p1mup6*F(325)))+4*s45*(p1mup6*p3mup2*F(326)+p1mup2*p3mup6*
     -   F(327)+p2mup2*p3mup6*F(328)+p3mup2*p3mup6*F(329)+2*(p3mup2*p
     -   5mup6*F(330)+p2mup6*p3mup2*F(331))+p3mup6*p5mup2*F(332))
       Fa15(206) = -(mup4mup6*s45*F(334))+p5mup4*((p4mup6+p5mup6)*F(33
     -   5)+p2mup6*F(337)-4*(p1mup6*F(336)-p3mup6*F(338)))
       Fa15(207) = 4*p4mup6*(p3mup3*F(341)+p5mup3*F(342)+p1mup3*F(343)
     -   )-p4mup3*(p5mup6*F(335)-4*(p1mup6*F(336)-p3mup6*F(338)-p4mup
     -   6*F(340))+p2mup6*F(345))+2*s45*(-(mup3mup6*F(339))+2*p2mup6*
     -   (p1mup3*F(346)+p3mup3*F(348)+p5mup3*F(350)))+p2mup3*(p4mup6*
     -   F(344)+4*s45*(p5mup6*F(330)+p2mup6*F(347)+p1mup6*F(349)+p3mu
     -   p6*F(351)))
       Fa15(208) = 2*mup2mup6*s45*F(352)+p4mup2*p4mup6*F(356)+p4mup6*p
     -   5mup2*F(357)+p2mup6*p4mup2*F(358)+p4mup2*p5mup6*F(361)+p1mup
     -   6*(p4mup2*F(360)-4*s45*(p5mup2*F(305)+p2mup2*F(362)+p3mup2*F
     -   (365)))-4*(p5mup2*p5mup6*s45*F(311)+p3mup2*p4mup6*F(353)+p2m
     -   up2*p4mup6*F(355)-p3mup6*p4mup2*F(359)+p2mup2*p2mup6*s45*F(3
     -   64)+p2mup6*p3mup2*s45*F(366)+p2mup2*p5mup6*s45*F(367)+p1mup2
     -   *(p4mup6*F(354)+s45*(p1mup6*F(296)+p5mup6*F(307)+p2mup6*F(36
     -   3)+p3mup6*F(368)))+p2mup2*p3mup6*s45*F(369)+p2mup6*p5mup2*s4
     -   5*F(370)+p3mup2*p5mup6*s45*F(371)+p3mup2*p3mup6*s45*F(372)+p
     -   3mup6*p5mup2*s45*F(373))
       Fa15(209) = mup3mup6*F(376)+p1mup6*p2mup3*F(378)+p2mup3*p2mup6*
     -   F(379)+p2mup3*(p4mup6+p5mup6)*F(380)+(p4mup3+p5mup3)*(4*(p4m
     -   up6+p5mup6)*F(386)-p3mup6*F(387)+p1mup6*F(388))+p3mup3*((p4m
     -   up6+p5mup6)*F(377)+4*(p1mup6*F(382)+p2mup6*F(384)+p3mup6*F(3
     -   90)))+4*(p1mup3*p3mup6*F(381)+p2mup3*p3mup6*F(383)-p2mup6*(p
     -   4mup3+p5mup3)*F(389)+p1mup3*(p2mup6*F(385)+p1mup6*F(391)))+p
     -   1mup3*(p4mup6+p5mup6)*F(392)
       Fa15(210) = mup2mup6*F(393)+p3mup6*(p1mup2*F(397)-p2mup2*F(398)
     -   )-p1mup2*p2mup6*F(400)-p3mup2*(4*((p4mup6+p5mup6)*F(394)-p1m
     -   up6*F(395))+p2mup6*F(396)+p3mup6*F(403))+p2mup6*p4mup2*F(404
     -   )+p2mup6*p5mup2*F(404)+p3mup6*p4mup2*F(405)+p3mup6*p5mup2*F(
     -   405)+p1mup6*p4mup2*F(406)+p1mup6*p5mup2*F(406)-4*(p2mup2*((p
     -   4mup6+p5mup6)*F(399)-p1mup6*F(401)+p2mup6*F(402))+p1mup2*(p4
     -   mup6+p5mup6)*F(407)+p1mup2*p1mup6*F(408)-(p4mup2+p5mup2)*(p4
     -   mup6+p5mup6)*F(409))
       Return
       End
