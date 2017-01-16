c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineAbe40(p1mup2,p1mup3,p1mup4,p1m
     -   up5,p2mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,
     -   p4mup2,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mu
     -   p2,p6mup3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,
     -   mup3mup5,mup4mup5,Fa40)
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
       Complex*16 Fa40(430:445)
       COMMON/test_FaFunctionsAbe/Fa
       Fa40(430) = Fa(277)-16*Fa(278)+Fa(279)+8*Fa(280)+Fa(281)+Fa(282
     -   )
       Fa40(431) = Fa(290)+Fa(292)-32*(Fa(293)-Fa(294)+Fa(295)-Fa(296)
     -   +Fa(297)-Fa(298))+Fa(299)+16*(Fa(288)+Fa(289)-Fa(291)+Fa(300
     -   ))
       Fa40(432) = Fa(302)+Fa(303)-32*(Fa(306)+Fa(307)-Fa(308)+Fa(309)
     -   +Fa(310)-Fa(311))+Fa(312)-16*(Fa(301)-Fa(304)-Fa(305)+Fa(313
     -   ))
       Fa40(433) = Fa(314)+Fa(315)+Fa(316)+Fa(317)-Fa(318)-2*(Fa(319)+
     -   Fa(320)+Fa(321)+Fa(322)+Fa(323)+Fa(324))+Fa(325)+Fa(326)
       Fa40(434) = Fa(327)+8*(Fa(328)+Fa(329))+Fa(330)+Fa(331)+Fa(332)
       Fa40(435) = -16*Fa(336)+Fa(337)+8*(Fa(333)+Fa(334)-Fa(335)+Fa(3
     -   38)+Fa(339))
       Fa40(436) = Fa(341)+Fa(342)+Fa(343)+2*(Fa(344)-Fa(345))+Fa(346)
     -   -Fa(347)
       Fa40(437) = 8*(Fa(349)+Fa(350)+Fa(351)-2*(Fa(352)-Fa(353)))+Fa(
     -   354)
       Fa40(438) = Fa(357)-Fa(358)-Fa(359)+2*(Fa(360)-Fa(361))+Fa(362)
       Fa40(439) = Fa(363)+Fa(364)-Fa(365)+2*(Fa(366)-Fa(367))+Fa(368)
     -   +Fa(369)
       Fa40(440) = Fa(382)-Fa(383)+Fa(384)+2*(Fa(380)-Fa(381)+Fa(385))
       Fa40(441) = Fa(388)+Fa(389)-Fa(390)+2*(Fa(386)-Fa(387)-Fa(391))
       Fa40(442) = -16*(Fa(396)+Fa(397))+8*Fa(398)+Fa(399)
       Fa40(443) = Fa(400)+Fa(401)-Fa(402)-Fa(403)+Fa(404)+Fa(405)+Fa(
     -   406)
       Fa40(444) = -2*(Fa(407)-Fa(408))+Fa(409)+Fa(410)
       Fa40(445) = -2*(Fa(411)-Fa(412))+Fa(413)+Fa(414)
       Return
       End
