c************** Calling the Fa functions*************************
       subroutine NoAbe2gaHexLineFa6(p1mup2,p1mup3,p1mup5,p1mup6,p2mup
     -   2,p2mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4
     -   mup3,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3
     -   ,p6mup5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,
     -   mup5mup6,Fa6)
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
       Complex*16 Fa6(121:144)
       COMMON/NoAbe2gaHexLineFaFunctions/Fa
       Fa6(121) = (p1mup2+p2mup2)*p4mup6*(4*p3mup5*F(398)+p5mup5*F(406
     -   )+(p1mup5+p2mup5)*F(407))
       Fa6(122) = p6mup2*(p2mup5*p2mup6*F(402)+p5mup5*(p5mup6*F(408)-p
     -   2mup6*F(409))+2*(p3mup5+p4mup5)*((p1mup6+p2mup6)*F(410)+p5mu
     -   p6*F(411)+p6mup6*F(412)))
       Fa6(123) = p1mup5*p6mup2*s12*(p5mup6*F(400)+p6mup6*F(401)+(p1mu
     -   p6+p2mup6)*F(402))+p1mup2*p3mup6*(4*p4mup5*F(398)+p5mup5*F(4
     -   06)+(p1mup5+p2mup5)*F(407))+p2mup2*p3mup6*(4*p4mup5*F(398)+p
     -   5mup5*F(406)+(p1mup5+p2mup5)*F(407))+p1mup2*(p3mup6+p4mup6)*
     -   p6mup5*F(413)+p2mup2*(p3mup6+p4mup6)*p6mup5*F(413)+s12*(-(p1
     -   mup6*p5mup5*p6mup2*F(409))+2*p5mup5*p6mup2*p6mup6*F(414)+(p3
     -   mup6+p4mup6)*p5mup2*p5mup5*F(416)+(p3mup6+p4mup6)*p5mup2*((p
     -   1mup5+p2mup5)*F(415)+p6mup5*F(417)))
       Fa6(124) = -(p1mup6*p6mup2*p6mup5*F(418))-p2mup6*p6mup2*p6mup5*
     -   F(418)-p3mup6*p5mup5*p6mup2*F(419)-p4mup6*p5mup5*p6mup2*F(41
     -   9)-p3mup6*p6mup2*p6mup5*F(420)-p4mup6*p6mup2*p6mup5*F(420)+p
     -   1mup5*p3mup6*p6mup2*F(421)+p2mup5*p3mup6*p6mup2*F(421)+p1mup
     -   5*p4mup6*p6mup2*F(421)+p2mup5*p4mup6*p6mup2*F(421)+4*p3mup5*
     -   p3mup6*p5mup2*F(422)+p5mup6*p6mup2*p6mup5*F(423)
       Fa6(125) = p1mup2*p3mup5*p3mup6*F(398)+p2mup2*p3mup5*p3mup6*F(3
     -   98)+(p3mup6*p4mup5+(p3mup5+p4mup5)*p4mup6)*p5mup2*s12*F(422)
     -   +(p3mup5+p4mup5)*(p3mup6+p4mup6)*p6mup2*s12*F(424)
       Fa6(126) = 2*(p1mup2+p2mup2)*F(372)+s12*(2*p6mup2*F(374)+p5mup2
     -   *F(426))
       Fa6(127) = 4*mup3mup6*(p1mup5*p5mup2*F(428)+p2mup5*p5mup2*F(428
     -   )-p5mup2*p5mup5*F(430)-p5mup2*p6mup5*F(431)+p1mup5*p6mup2*F(
     -   433)+p2mup5*p6mup2*F(433)+p5mup5*p6mup2*F(434)-p6mup2*p6mup5
     -   *F(435)+p3mup5*p5mup2*F(439)+p4mup5*p5mup2*F(439)+(p3mup5+p4
     -   mup5)*p6mup2*F(440))+mup3mup5*(p5mup2*((p1mup6+p2mup6)*F(429
     -   )+4*(p5mup6*F(427)-p6mup6*F(432)+(p3mup6+p4mup6)*F(441)))+p6
     -   mup2*(p5mup6*F(436)+(p1mup6+p2mup6)*F(437)-4*(p6mup6*F(438)-
     -   (p3mup6+p4mup6)*F(442))))
       Fa6(128) = mup3mup5*((p1mup2+p2mup2)*F(443)+s12*(p5mup2*F(444)+
     -   p6mup2*F(445)))
       Fa6(129) = 4*mup3mup6*(-(p1mup5*p2mup2*F(428))-p2mup2*p2mup5*F(
     -   428)+p2mup2*p5mup5*F(430)+p2mup2*p6mup5*F(431)-p2mup2*p3mup5
     -   *F(439)-p2mup2*p4mup5*F(439)-p1mup2*((p1mup5+p2mup5)*F(428)-
     -   p5mup5*F(430)-p6mup5*F(431)+(p3mup5+p4mup5)*F(439))+p5mup5*p
     -   6mup2*s12*F(448)+2*(p6mup2*p6mup5*s12*F(449)-(p1mup5+p2mup5)
     -   *p6mup2*s12*F(452)-(p3mup5+p4mup5)*p6mup2*s12*F(454)))+mup3m
     -   up5*(p1mup2*(-4*(p5mup6*F(427)+(p3mup6+p4mup6)*F(441))+(p1mu
     -   p6+p2mup6)*F(446)+p6mup6*F(447))+p2mup2*(-4*(p5mup6*F(427)+(
     -   p3mup6+p4mup6)*F(441))+(p1mup6+p2mup6)*F(446)+p6mup6*F(447))
     -   +4*p6mup2*s12*(-(p5mup6*F(450))-p6mup6*F(451)+(p1mup6+p2mup6
     -   )*F(453)+(p3mup6+p4mup6)*F(455)))
       Fa6(130) = 4*mup3mup6*(p5mup2*s12*(-(p5mup5*F(448))+2*(-(p6mup5
     -   *F(449))+(p1mup5+p2mup5)*F(452)+(p3mup5+p4mup5)*F(454)))+(p1
     -   mup2+p2mup2)*((p1mup5+p2mup5)*F(456)+p5mup5*F(458)+p6mup5*F(
     -   460)-(p3mup5+p4mup5)*F(462)))+mup3mup5*(4*p5mup2*s12*(p5mup6
     -   *F(450)+p6mup6*F(451)-(p1mup6+p2mup6)*F(453)-(p3mup6+p4mup6)
     -   *F(455))+(p1mup2+p2mup2)*((p1mup6+p2mup6)*F(457)-p5mup6*F(45
     -   9)+4*(p6mup6*F(461)-(p3mup6+p4mup6)*F(463))))
       Fa6(131) = mup5mup6*F(464)-(p1mup5+p2mup5)*p5mup6*F(466)+p5mup5
     -   *p6mup6*F(468)+p5mup6*p6mup5*F(469)-p5mup5*p5mup6*F(471)+p1m
     -   up6*((p1mup5+p2mup5)*F(465)+p5mup5*F(467)+2*(p6mup5*F(379)-(
     -   p3mup5+p4mup5)*F(472)))+p2mup6*((p1mup5+p2mup5)*F(465)+p5mup
     -   5*F(467)+2*(p6mup5*F(379)-(p3mup5+p4mup5)*F(472)))+2*(p6mup5
     -   *p6mup6*F(387)+(p3mup5+p4mup5)*(-(p5mup6*F(396))-4*(p3mup6+p
     -   4mup6)*F(398))+(p1mup5+p2mup5)*(p6mup6*F(384)-(p3mup6+p4mup6
     -   )*F(407))+p3mup5*p6mup6*F(470)+p4mup5*p6mup6*F(470)+p3mup6*p
     -   5mup5*F(473)+p4mup6*p5mup5*F(473)+(p3mup6+p4mup6)*p6mup5*F(4
     -   74))
       Fa6(132) = mup3mup6*(-2*((p1mup5+p2mup5)*F(428)-p5mup5*F(430))+
     -   p6mup5*F(478)+(p3mup5+p4mup5)*F(480))+mup3mup5*((p1mup6+p2mu
     -   p6)*F(476)+p5mup6*F(477)+p6mup6*F(479)+(p3mup6+p4mup6)*F(481
     -   ))
       Fa6(133) = mup3mup6*(2*(p5mup5*F(458)+p6mup5*F(460))+(p1mup5+p2
     -   mup5)*F(482)+(p3mup5+p4mup5)*F(486))+mup3mup5*((p1mup6+p2mup
     -   6)*F(483)+p5mup6*F(484)+p6mup6*F(485)+(p3mup6+p4mup6)*F(487)
     -   )
       Fa6(134) = -4*(p6mup5*p6mup6*F(391)+p3mup5*p6mup6*F(403)+p4mup5
     -   *p6mup6*F(403)-p1mup6*p3mup5*F(404)-p2mup6*p3mup5*F(404)-p1m
     -   up6*p4mup5*F(404)-p2mup6*p4mup5*F(404)-p3mup5*p5mup6*F(405)-
     -   p4mup5*p5mup6*F(405))+8*(p3mup5*p3mup6*F(422)+p3mup6*p4mup5*
     -   F(422)+p3mup5*p4mup6*F(422)+p4mup5*p4mup6*F(422))-mup5mup6*F
     -   (488)+p1mup6*p6mup5*F(490)+p2mup6*p6mup5*F(490)+2*(p5mup5*p5
     -   mup6*F(382)+p1mup6*p5mup5*F(393)+p2mup6*p5mup5*F(393)+p5mup5
     -   *p6mup6*F(395)+p3mup6*p5mup5*F(416)+p4mup6*p5mup5*F(416)+p1m
     -   up5*(p5mup6*F(381)-2*p6mup6*F(390)-(p1mup6+p2mup6)*F(489)-(p
     -   3mup6+p4mup6)*F(491))+p2mup5*(p5mup6*F(381)-2*p6mup6*F(390)-
     -   (p1mup6+p2mup6)*F(489)-(p3mup6+p4mup6)*F(491)))+p3mup6*p6mup
     -   5*F(492)+p4mup6*p6mup5*F(492)+p5mup6*p6mup5*F(493)
       Fa6(135) = 4*(p1mup6*p3mup5*F(410)+p2mup6*p3mup5*F(410)+p1mup6*
     -   p4mup5*F(410)+p2mup6*p4mup5*F(410)+p3mup5*p5mup6*F(411)+p4mu
     -   p5*p5mup6*F(411)+p3mup5*p6mup6*F(412)+p4mup5*p6mup6*F(412)+p
     -   5mup5*p6mup6*F(414))+8*(p3mup5*p3mup6*F(424)+p3mup6*p4mup5*F
     -   (424)+p3mup5*p4mup6*F(424)+p4mup5*p4mup6*F(424))-2*(-(p5mup5
     -   *p5mup6*F(408))+p1mup6*p5mup5*F(409)+p2mup6*p5mup5*F(409)+p3
     -   mup6*p5mup5*F(419)+p4mup6*p5mup5*F(419)-p1mup5*(p5mup6*F(400
     -   )+p6mup6*F(401)+(p1mup6+p2mup6)*F(402)+(p3mup6+p4mup6)*F(421
     -   ))-p2mup5*(p5mup6*F(400)+p6mup6*F(401)+(p1mup6+p2mup6)*F(402
     -   )+(p3mup6+p4mup6)*F(421))+p6mup5*p6mup6*F(425))-mup5mup6*F(4
     -   94)-p1mup6*p6mup5*F(495)-p2mup6*p6mup5*F(495)-p3mup6*p6mup5*
     -   F(496)-p4mup6*p6mup5*F(496)-p5mup6*p6mup5*F(497)
       Fa6(136) = (p1mup6+p2mup6)*F(499)+p5mup6*F(500)-p6mup6*F(501)+(
     -   p3mup6+p4mup6)*F(502)
       Fa6(137) = p1mup6*p5mup2*F(476)+p2mup6*p5mup2*F(476)+p5mup2*p5m
     -   up6*F(477)+p5mup2*p6mup6*F(479)+p3mup6*p5mup2*F(481)+p4mup6*
     -   p5mup2*F(481)+p1mup6*p6mup2*F(483)+p2mup6*p6mup2*F(483)+p5mu
     -   p6*p6mup2*F(484)+p6mup2*p6mup6*F(485)+(p3mup6+p4mup6)*p6mup2
     -   *F(487)
       Fa6(138) = 2*(p5mup6*p6mup2*s12*F(450)+p6mup2*s12*(p6mup6*F(451
     -   )-(p1mup6+p2mup6)*F(453)-(p3mup6+p4mup6)*F(455)))+(p1mup2+p2
     -   mup2)*p5mup6*F(504)+(p1mup2+p2mup2)*((p1mup6+p2mup6)*F(503)+
     -   p6mup6*F(505)+(p3mup6+p4mup6)*F(506))
       Fa6(139) = 2*p5mup2*s12*(-(p5mup6*F(450))-p6mup6*F(451)+(p1mup6
     -   +p2mup6)*F(453)+(p3mup6+p4mup6)*F(455))+(p1mup2+p2mup2)*((p1
     -   mup6+p2mup6)*F(507)+p5mup6*F(508)+p6mup6*F(509)+(p3mup6+p4mu
     -   p6)*F(510))
       Fa6(140) = (p3mup5+p4mup5)*F(370)+(p1mup5+p2mup5)*F(511)+p5mup5
     -   *F(512)+p6mup5*F(513)
       Fa6(141) = p6mup2*(-(mup5mup6*F(516))+2*(-(p5mup5*p6mup6*F(522)
     -   )+(p1mup6+p2mup6)*(p5mup5*F(518)-4*((p3mup5+p4mup5)*F(520)-p
     -   6mup5*F(524)))+p6mup5*p6mup6*F(525)+2*(p5mup5*p5mup6*F(521)+
     -   p3mup5*p6mup6*F(527)+p4mup5*p6mup6*F(527)+p3mup6*p5mup5*F(52
     -   8)+p4mup6*p5mup5*F(528))+p1mup5*(-(p6mup6*F(517))+4*(-(p5mup
     -   6*F(519))+(p3mup6+p4mup6)*F(523)+(p1mup6+p2mup6)*F(531)))+p2
     -   mup5*(-(p6mup6*F(517))+4*(-(p5mup6*F(519))+(p3mup6+p4mup6)*F
     -   (523)+(p1mup6+p2mup6)*F(531)))-4*(p3mup5*p5mup6*F(526)+p4mup
     -   5*p5mup6*F(526)-p3mup6*p6mup5*F(529)-p4mup6*p6mup5*F(529)+p5
     -   mup6*p6mup5*F(530)-(p3mup5+p4mup5)*(p3mup6+p4mup6)*F(9576)))
     -   )
       Fa6(142) = p5mup2*(4*(p5mup5*p6mup6*F(522)-(p1mup6+p2mup6)*(p5m
     -   up5*F(518)-4*((p3mup5+p4mup5)*F(520)-p6mup5*F(524)))-p6mup5*
     -   p6mup6*F(525)+p1mup5*(p6mup6*F(517)-4*(-(p5mup6*F(519))+(p3m
     -   up6+p4mup6)*F(523)+(p1mup6+p2mup6)*F(531)))+p2mup5*(p6mup6*F
     -   (517)-4*(-(p5mup6*F(519))+(p3mup6+p4mup6)*F(523)+(p1mup6+p2m
     -   up6)*F(531))))+mup5mup6*F(532)-8*(p3mup6*p5mup5*F(528)+p4mup
     -   6*p5mup5*F(528)-p5mup5*p5mup6*F(533))+p3mup5*p6mup6*F(534)+p
     -   4mup5*p6mup6*F(534)-16*(-((p3mup5+p4mup5)*p5mup6*F(526))+p3m
     -   up6*p6mup5*F(529)+p4mup6*p6mup5*F(529)-p5mup6*p6mup5*F(530)+
     -   (p3mup5+p4mup5)*(p3mup6+p4mup6)*F(9576)))
       Fa6(143) = -2*(p1mup5*p5mup2*F(428)+p2mup5*p5mup2*F(428)-p5mup2
     -   *p5mup5*F(430)-p5mup5*p6mup2*F(458)-p6mup2*p6mup5*F(460))+p5
     -   mup2*p6mup5*F(478)+p3mup5*p5mup2*F(480)+p4mup5*p5mup2*F(480)
     -   +p1mup5*p6mup2*F(482)+p2mup5*p6mup2*F(482)+(p3mup5+p4mup5)*p
     -   6mup2*F(486)
       Fa6(144) = 2*p6mup2*s12*(-(p5mup5*F(448))+2*(-(p6mup5*F(449))+(
     -   p1mup5+p2mup5)*F(452)+(p3mup5+p4mup5)*F(454)))+p1mup2*(-2*p5
     -   mup5*F(430)+(p1mup5+p2mup5)*F(536)+p6mup5*F(537)+(p3mup5+p4m
     -   up5)*F(538))+p2mup2*(-2*p5mup5*F(430)+(p1mup5+p2mup5)*F(536)
     -   +p6mup5*F(537)+(p3mup5+p4mup5)*F(538))
       Return
       End
