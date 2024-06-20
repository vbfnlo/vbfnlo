c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineNoAbe1g16(p1mup2,p1mup3,p1mup4
     -   ,p1mup6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3m
     -   up6,p4mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,
     -   p6mup2,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3m
     -   up4,mup3mup6,mup4mup6,Fa16)
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
       Complex*16 Fa16(211:224)
       COMMON/test_FaFunctionsNoAbe1g/Fa
       Fa16(211) = -(mup3mup6*p5mup4*F(411))+p2mup3*(-4*mup4mup6*s45*F
     -   (413)+8*p5mup4*(p3mup6*F(383)+(p4mup6+p5mup6)*F(415)+p2mup6*
     -   F(418)+p1mup6*F(419)))-mup4mup6*p5mup3*s45*F(423)-4*(mup4mup
     -   6*p3mup3*s45*F(414)+mup4mup6*p1mup3*s45*F(424))+p4mup3*(-(mu
     -   p4mup6*F(412))+8*p5mup4*((p4mup6+p5mup6)*F(386)-p2mup6*F(389
     -   )-p3mup6*F(425)-p1mup6*F(426)))+8*(p4mup6*p5mup3*p5mup4*F(38
     -   6)+p5mup3*p5mup4*p5mup6*F(386)-p2mup6*p5mup3*p5mup4*F(389)+p
     -   3mup3*p3mup6*p5mup4*F(390)+p1mup3*p1mup6*p5mup4*F(391)+p3mup
     -   3*p4mup6*p5mup4*F(416)+p3mup3*p5mup4*p5mup6*F(416)+p2mup6*p3
     -   mup3*p5mup4*F(417)+p1mup3*p3mup6*p5mup4*F(420)+p1mup6*p3mup3
     -   *p5mup4*F(421)+p1mup3*p2mup6*p5mup4*F(422)-p3mup6*p5mup3*p5m
     -   up4*F(425)-p1mup6*p5mup3*p5mup4*F(426)+p1mup3*p5mup4*(p4mup6
     -   +p5mup6)*F(427))
       Fa16(212) = p1mup2*p3mup6*p5mup4*F(431)+p1mup6*p3mup2*p5mup4*F(
     -   434)+p3mup2*p3mup6*p5mup4*F(440)+4*(mup2mup6*p5mup4*F(428)+m
     -   up4mup6*p3mup2*s45*F(430)+mup4mup6*p5mup2*s45*F(441)+mup4mup
     -   6*p1mup2*s45*F(442))+mup4mup6*p2mup2*s45*F(443)+p2mup6*p5mup
     -   2*p5mup4*F(444)+p1mup6*p5mup2*p5mup4*F(446)+p4mup2*(mup4mup6
     -   *F(429)+p5mup4*(p2mup6*F(444)-8*((p4mup6+p5mup6)*F(409)+p3mu
     -   p6*F(445))+p1mup6*F(446)))-8*(p5mup2*p5mup4*p5mup6*F(409)-p2
     -   mup2*p3mup6*p5mup4*F(432)+p2mup6*p3mup2*p5mup4*F(433)-p3mup2
     -   *p5mup4*p5mup6*F(435)-p2mup2*p5mup4*p5mup6*F(436)-p1mup2*p2m
     -   up6*p5mup4*F(437)-p1mup6*p2mup2*p5mup4*F(438)-p2mup2*p2mup6*
     -   p5mup4*F(439)+p3mup6*p5mup2*p5mup4*F(445)-p1mup2*p5mup4*p5mu
     -   p6*F(447)-p4mup6*p5mup4*(-(p5mup2*F(409))+p3mup2*F(435)+p2mu
     -   p2*F(436)+p1mup2*F(447)))+p1mup2*p1mup6*p5mup4*F(448)
       Fa16(213) = p2mup2*p3mup3*p4mup6*F(454)-p3mup2*p4mup6*p5mup3*F(
     -   460)+p3mup3*p4mup6*p5mup2*F(461)+p2mup2*p4mup3*p4mup6*F(465)
     -   +p3mup2*p4mup3*p4mup6*F(468)+p3mup2*p3mup3*p4mup6*F(469)+mup
     -   3mup6*p5mup2*s45*F(471)-p3mup3*p4mup2*p4mup6*F(476)+p4mup6*p
     -   5mup2*p5mup3*F(480)+p2mup2*p4mup6*p5mup3*F(481)+p2mup6*p4mup
     -   2*p5mup3*F(482)+p3mup6*p4mup2*p5mup3*F(484)+4*(mup3mup6*p4mu
     -   p2*F(450)+mup3mup6*p3mup2*s45*F(462)+mup3mup6*p1mup2*s45*F(4
     -   72)+mup3mup6*p2mup2*s45*F(473)+mup2mup3*p2mup6*s45*F(486)+mu
     -   p2mup3*p3mup6*s45*F(487))+mup2mup3*p4mup6*F(488)+mup2mup6*(p
     -   4mup3*F(449)-4*s45*(p2mup3*F(451)+p3mup3*F(452)+p1mup3*F(474
     -   )-p5mup3*F(495)))+p2mup3*(p4mup6*(p5mup2*F(459)+p2mup2*F(463
     -   )-p3mup2*F(466))-8*(p4mup2*(p4mup6*F(453)+p3mup6*F(456)+p2mu
     -   p6*F(464))-s45*(p2mup2*p2mup6*F(491)+p2mup2*p3mup6*F(493)-p2
     -   mup6*p3mup2*F(499)+p2mup6*p5mup2*F(500)-p3mup2*p3mup6*F(501)
     -   +p3mup6*p5mup2*F(503))))+8*(p2mup2*p3mup6*p4mup3*F(455)+p2mu
     -   p6*p3mup2*p4mup3*F(457)-p2mup6*p3mup3*p4mup2*F(458)+p2mup2*p
     -   2mup6*p4mup3*F(467)+p3mup2*p3mup6*p4mup3*F(470)+p2mup6*p4mup
     -   2*p4mup3*F(475)+p3mup6*p4mup2*p4mup3*F(477)+p4mup2*p4mup6*p5
     -   mup3*F(478)-p4mup3*p4mup6*p5mup2*F(479)-p2mup6*p4mup3*p5mup2
     -   *F(483)-p3mup6*p4mup3*p5mup2*F(485)+p4mup2*p4mup3*p4mup6*F(4
     -   89)-p3mup3*p3mup6*p4mup2*F(490)+p2mup2*p2mup6*p3mup3*s45*F(4
     -   92)+p2mup2*p3mup3*p3mup6*s45*F(494)+p2mup2*p3mup6*p5mup3*s45
     -   *F(496)-p2mup6*p3mup2*p5mup3*s45*F(497)+p2mup6*p3mup3*p5mup2
     -   *s45*F(498)+p3mup3*p3mup6*p5mup2*s45*F(502)+p3mup2*p3mup6*p5
     -   mup3*s45*F(502)+p3mup2*p3mup3*p3mup6*s45*F(504)+p2mup6*p5mup
     -   2*p5mup3*s45*F(505)+p3mup6*p5mup2*p5mup3*s45*F(506)+p2mup6*p
     -   3mup2*p3mup3*s45*F(507)+p2mup2*p2mup6*p5mup3*s45*F(508))
       Fa16(214) = p1mup3*(p4mup6*p5mup2*F(525)+p2mup2*p4mup6*F(535)+p
     -   3mup2*(p4mup6*F(515)-8*s45*(p2mup6*F(556)+p3mup6*F(557)))-8*
     -   (p2mup6*p4mup2*F(519)+p4mup2*p4mup6*F(524)-p2mup2*p2mup6*s45
     -   *F(539)-p2mup6*p5mup2*s45*F(553)-p3mup6*(-(p4mup2*F(516))+p2
     -   mup2*s45*F(546)+p5mup2*s45*F(558))))+p1mup6*(p3mup2*p4mup3*F
     -   (510)+p4mup2*p5mup3*F(520)+4*mup2mup3*s45*F(522)+p4mup2*p4mu
     -   p3*F(523)+p1mup2*p4mup3*F(533)-p1mup3*p4mup2*F(534)+p3mup3*(
     -   p4mup2*F(511)+8*s45*(p2mup2*F(527)+p1mup2*F(531)+p3mup2*F(54
     -   1)+p5mup2*F(551)))+8*(p2mup2*p4mup3*F(517)-p4mup3*p5mup2*F(5
     -   21)+p1mup3*s45*(p1mup2*F(528)+p2mup2*F(530))+p1mup3*p3mup2*s
     -   45*F(545)+p2mup2*p5mup3*s45*F(549)+p2mup3*(-(p4mup2*F(509))+
     -   s45*(p2mup2*F(526)+p1mup2*F(529)+p3mup2*F(542)+p5mup2*F(550)
     -   ))+p1mup3*p5mup2*s45*F(552)+p1mup2*p5mup3*s45*F(552)+p3mup2*
     -   p5mup3*s45*F(555)+p5mup2*p5mup3*s45*F(559)))+p1mup2*(p1mup3*
     -   p4mup6*F(532)+p4mup3*p4mup6*F(536)+p2mup3*(p4mup6*F(512)+8*s
     -   45*(p2mup6*F(538)+p3mup6*F(544)))+p3mup3*(p4mup6*F(513)+8*(p
     -   2mup6*s45*F(540)+p3mup6*s45*F(547)))+p4mup6*p5mup3*F(548)+8*
     -   (p3mup6*p4mup3*F(514)+p2mup6*p4mup3*F(518)+p1mup3*p2mup6*s45
     -   *F(537)+p1mup3*p3mup6*s45*F(543)+p3mup6*p5mup3*s45*F(554)-p2
     -   mup6*p5mup3*s45*F(560)))
       Fa16(215) = p5mup6*(p1mup2*(p1mup3*F(561)+p2mup3*F(562)+p3mup3*
     -   F(564))+p2mup2*(-(p3mup3*F(563))+p1mup3*F(565)+p2mup3*F(566)
     -   ))
       Fa16(216) = p5mup6*(p3mup2*p4mup3*F(569)+p2mup2*p4mup3*F(570)+p
     -   4mup2*p4mup3*F(571)+p4mup2*p5mup3*F(572)-8*(p2mup3*p4mup2*F(
     -   567)+p3mup3*p4mup2*F(568)+p1mup3*p4mup2*F(573))+p1mup2*p4mup
     -   3*F(574))
       Fa16(217) = p5mup6*(2*p4mup3*p5mup2*F(576)+s45*(-(mup2mup3*F(57
     -   5))+2*(p1mup3*p3mup2*F(564)+p2mup3*p3mup2*F(577)+p1mup3*p5mu
     -   p2*F(578)+p1mup2*p5mup3*F(578)+p2mup3*p5mup2*F(579)+p3mup2*p
     -   3mup3*F(580)+p2mup2*p5mup3*F(581)+p3mup3*p5mup2*F(582)+p3mup
     -   2*p5mup3*F(582)+p5mup2*p5mup3*F(583))))
       Fa16(218) = (p4mup3*F(585))/s45+p2mup3*F(586)+p3mup3*F(587)-2*p
     -   5mup3*F(588)+p1mup3*F(589)
       Fa16(219) = (p4mup2*F(590))/s45+p3mup2*F(591)+2*(p2mup2*F(592)+
     -   p5mup2*F(593)+p1mup2*F(594))
       Fa16(220) = mup4mup6*(p3mup3*p4mup2*F(596)+p4mup2*p5mup3*F(597)
     -   +p4mup2*p4mup3*F(598)-8*(p2mup2*p2mup3*s45*F(599)-p3mup2*(p4
     -   mup3*F(595)+p5mup3*s45*F(600))+p5mup2*s45*(p3mup3*F(601)+p5m
     -   up3*F(602))))
       Fa16(221) = mup4mup6*(p2mup3*(p4mup2*F(603)-8*s45*(p1mup2*F(608
     -   )+p5mup2*F(611)-p3mup2*F(612)))-8*(-(p4mup3*(p5mup2*F(604)+p
     -   2mup2*F(605)))+s45*(p2mup2*p3mup3*F(606)+p1mup2*p1mup3*F(607
     -   )+p1mup3*p2mup2*F(609)+p1mup2*p3mup3*F(610)-p3mup2*p3mup3*F(
     -   613)+mup2mup3*F(614)-p2mup2*p5mup3*F(615))))
       Fa16(222) = mup4mup6*(p1mup3*(p4mup2*F(617)-8*(p5mup2*s45*F(618
     -   )-p3mup2*s45*F(619)))+8*p1mup2*(p4mup3*F(616)+p5mup3*s45*F(6
     -   20)))
       Fa16(223) = mup3mup6*p5mup4*(p3mup2*F(621)+p1mup2*F(622)+p2mup2
     -   *F(623)+p4mup2*F(624))
       Fa16(224) = p5mup4*(mup3mup6*p5mup2*F(624)+8*mup2mup6*(p4mup3+p
     -   5mup3)*F(625))
       Return
       End
