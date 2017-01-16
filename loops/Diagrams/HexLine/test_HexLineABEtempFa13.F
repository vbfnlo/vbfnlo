c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineAbe13(p1mup2,p1mup3,p1mup4,p1m
     -   up5,p2mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,
     -   p4mup2,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mu
     -   p2,p6mup3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,
     -   mup3mup5,mup4mup5,Fa13)
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
       Complex*16 Fa13(133:143)
       COMMON/test_FaFunctionsAbe/Fa
       Fa13(133) = p4mup3*F(288)+p1mup3*F(289)+p2mup3*F(290)+p3mup3*F(
     -   291)
       Fa13(134) = mup2mup4*(p1mup3*F(292)+p2mup3*F(293)+p3mup3*F(294)
     -   +p4mup3*F(295))+mup3mup4*(p1mup2*F(296)+p2mup2*F(297)+p3mup2
     -   *F(299)+p4mup2*F(301))+mup2mup3*(p2mup4*F(298)+p3mup4*F(300)
     -   +p4mup4*F(302)+p1mup4*F(303))+p2mup2*p2mup3*p2mup4*F(307)+p2
     -   mup2*p2mup3*p3mup4*F(308)+p2mup2*p3mup3*p4mup4*F(309)+p2mup2
     -   *p3mup4*p4mup3*F(311)+p2mup3*p3mup2*p3mup4*F(312)+p2mup3*p2m
     -   up4*p3mup2*F(313)+p2mup3*p2mup4*p4mup2*F(314)+p2mup3*p3mup2*
     -   p4mup4*F(315)+p2mup4*p3mup2*p4mup3*F(316)+p2mup4*p3mup3*p4mu
     -   p2*F(317)-p1mup3*p3mup2*p4mup4*F(318)+p1mup3*p3mup4*p4mup2*F
     -   (319)-p1mup4*p3mup2*p4mup3*F(320)+p1mup4*p3mup3*p4mup2*F(321
     -   )+p1mup3*p2mup2*p3mup4*F(323)+p1mup3*p2mup2*p4mup4*F(324)+p1
     -   mup4*p2mup2*p2mup3*F(325)+p2mup2*p2mup3*p4mup4*F(326)+p2mup2
     -   *p3mup3*p3mup4*F(327)-p1mup4*p2mup2*p4mup3*F(329)-p1mup3*p2m
     -   up4*p4mup2*F(330)+p1mup4*p2mup3*p4mup2*F(331)+p2mup3*p3mup4*
     -   p4mup2*F(332)+p1mup4*p2mup2*p3mup3*F(333)+p2mup2*p2mup4*p3mu
     -   p3*F(334)+p1mup3*p3mup2*p3mup4*F(335)+p2mup4*p3mup2*p3mup3*F
     -   (336)+p3mup2*p3mup3*p3mup4*F(337)+p3mup2*p3mup4*p4mup3*F(338
     -   )-p1mup3*p2mup4*p3mup2*F(339)+p2mup2*p2mup4*p4mup3*F(340)+p3
     -   mup3*p3mup4*p4mup2*F(341)+p2mup3*p4mup2*p4mup4*F(342)+p2mup4
     -   *p4mup2*p4mup3*F(343)+p3mup4*p4mup2*p4mup3*F(344)+p1mup4*p3m
     -   up2*p3mup3*F(345)+p3mup2*p3mup3*p4mup4*F(346)+p1mup4*p4mup2*
     -   p4mup3*F(347)+p4mup2*p4mup3*p4mup4*F(348)-p1mup4*p2mup3*p3mu
     -   p2*F(349)+p2mup2*p4mup3*p4mup4*F(350)+p3mup2*p4mup3*p4mup4*F
     -   (351)+p3mup3*p4mup2*p4mup4*F(352)-p1mup3*p1mup4*p2mup2*F(361
     -   )+p1mup3*p2mup2*p2mup4*F(362)-p1mup3*p1mup4*p3mup2*F(364)-p1
     -   mup3*p1mup4*p4mup2*F(365)+p1mup2*(p3mup4*p4mup3*F(310)+p2mup
     -   4*p4mup3*F(328)-p1mup3*p1mup4*F(353)+p1mup3*p2mup4*F(354)+p1
     -   mup3*p3mup4*F(355)+p1mup3*p4mup4*F(356)+p2mup3*(p2mup4*F(304
     -   )+p3mup4*F(305)+p4mup4*F(322)-p1mup4*F(357))+p3mup3*(p4mup4*
     -   F(306)-p1mup4*F(358)+p2mup4*F(359)+p3mup4*F(360))-p1mup4*p4m
     -   up3*F(363)+p4mup3*p4mup4*F(366))+p1mup3*p4mup2*p4mup4*F(367)
       Fa13(135) = mup2mup3*F(368)+p1mup3*p3mup2*F(373)+p2mup3*p3mup2*
     -   F(374)+p3mup2*p4mup3*F(375)+p2mup3*p4mup2*F(376)+p3mup3*p4mu
     -   p2*F(377)+p3mup2*p3mup3*F(379)+p2mup2*(p2mup3*F(370)+p3mup3*
     -   F(372)+p1mup3*F(378)+p4mup3*F(380))+p1mup3*p4mup2*F(381)+p4m
     -   up2*p4mup3*F(382)+p1mup2*(p2mup3*F(369)+p3mup3*F(371)+p1mup3
     -   *F(383)+p4mup3*F(384))
       Fa13(136) = 4*p1mup4*F(385)+p2mup4*F(386)+p3mup4*F(387)+p4mup4*
     -   F(388)
       Fa13(137) = mup2mup4*F(389)+p3mup4*p4mup2*F(391)+p1mup4*p2mup2*
     -   F(397)+p2mup2*p3mup4*F(398)+p1mup4*p4mup2*F(399)+p3mup2*(p4m
     -   up4*F(390)+p2mup4*F(392)+p3mup4*F(393)+p1mup4*F(400))+p2mup2
     -   *p4mup4*F(401)+p1mup2*(p3mup4*F(395)+p4mup4*F(396)+8*(p2mup4
     -   *F(394)+p1mup4*F(403)))+8*(p2mup4*p4mup2*F(402)-p2mup2*p2mup
     -   4*F(404)-p4mup2*p4mup4*F(405))
       Fa13(138) = mup3mup4*F(406)+p3mup4*p4mup3*F(410)+p1mup3*p2mup4*
     -   F(411)+p1mup3*p3mup4*F(412)+p1mup3*p4mup4*F(413)+p1mup4*p4mu
     -   p3*F(417)+p3mup3*(p4mup4*F(409)+p2mup4*F(415)+p3mup4*F(416)+
     -   p1mup4*F(418))+p2mup4*p4mup3*F(419)+p4mup3*p4mup4*F(420)+p1m
     -   up3*p1mup4*F(421)+p2mup3*(p2mup4*F(407)+p3mup4*F(408)+p4mup4
     -   *F(414)-8*p1mup4*F(422))
       Fa13(139) = mup2mup4*F(423)+p3mup2*p3mup4*F(425)+p3mup2*p4mup4*
     -   F(426)+p1mup2*p3mup4*F(428)+p1mup2*p4mup4*F(429)+p2mup2*p3mu
     -   p4*F(431)+p2mup2*p4mup4*F(432)+p1mup4*p4mup2*F(433)+p2mup4*(
     -   p3mup2*F(424)+p1mup2*F(427)+p2mup2*F(430)+p4mup2*F(434))+p1m
     -   up4*p2mup2*F(435)+p3mup4*p4mup2*F(436)+p4mup2*p4mup4*F(437)+
     -   p1mup4*p3mup2*F(438)-8*p1mup2*p1mup4*F(439)
       Fa13(140) = p4mup2*F(441)+p2mup2*F(442)+p3mup2*F(443)+p1mup2*F(
     -   444)
       Fa13(141) = p3mup2*F(447)+p4mup2*F(448)+p1mup2*F(449)+p2mup2*F(
     -   450)
       Fa13(142) = p2mup3*F(451)+p3mup3*F(452)+p1mup3*F(453)+p4mup3*F(
     -   454)
       Fa13(143) = p3mup2*F(455)+p4mup2*F(456)+p1mup2*F(457)+p2mup2*F(
     -   458)
       Return
       End
