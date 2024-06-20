c************** Calling the Fa functions*************************
       subroutine HexBoxFa13(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4
     &   ,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa13)
       IMPLICIT NONE
      Complex*16   p1mup3, p1mup4, p2mup3, 
     -          p2mup4, p3mup3, p3mup4, 
     -          p4mup3, p4mup4, p5mup3, 
     -          p5mup4, p6mup3, p6mup4
       Complex*16   mup3mup4

        common/HexBoxFhlFunctions/F
       COMMON/HexBoxInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12,s23,s
     &   34,s45,s56,s16,s123,s234,s345
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
      Complex*16 Fa(259),F(1314)
      Real*8 P(1362) 
        COMMON/HexBoxPFunctions/P
       Complex*16 Fa13(145:156)
       COMMON/HexBoxFaFunctions/Fa
       Fa13(145) = 2*(p6mup3*F(281)+p1mup3*F(290)+p2mup3*F(295))+(p4mu
     &   p3+p5mup3)*F(456)
       Fa13(146) = p1mup3*F(337)-(p4mup3+p5mup3)*F(338)+p6mup3*F(339)-
     &   p2mup3*F(340)
       Fa13(147) = -(p1mup3*F(401))-p2mup3*F(402)+(p4mup3+p5mup3)*F(40
     &   3)+p6mup3*F(1311)
       Fa13(148) = -4*mup3mup4*F(459)+p1mup3*(p4mup4*F(460)+p5mup4*F(4
     &   62)+8*(p1mup4*F(461)+p6mup4*F(463)-p2mup4*F(464)))+p1mup4*(p
     &   4mup3*F(465)+8*(p5mup3*F(466)+p6mup3*F(467)-p2mup3*F(468)))+
     &   p2mup3*p4mup4*F(469)+p4mup3*p4mup4*F(472)+p4mup3*p6mup4*F(47
     &   4)+p4mup4*p6mup3*F(476)+p2mup4*p4mup3*F(482)+8*(p2mup3*p5mup
     &   4*F(470)+p2mup3*p6mup4*F(471)+p4mup3*p5mup4*F(473)+p4mup4*p5
     &   mup3*F(475)+p5mup3*p5mup4*F(477)+p5mup3*p6mup4*F(478)+p5mup4
     &   *p6mup3*F(479)+p6mup3*p6mup4*F(480)-p2mup3*p2mup4*F(481)+p2m
     &   up4*p5mup3*F(483)+p2mup4*p6mup3*F(484))
       Fa13(149) = p1mup4*F(486)+p2mup4*F(487)-p6mup4*F(488)-p4mup4*F(
     &   489)-p5mup4*F(490)
       Fa13(150) = p1mup3*F(491)-p2mup3*F(492)-p6mup3*F(493)-p4mup3*F(
     &   494)-p5mup3*F(495)
       Fa13(151) = mup3mup4*F(496)-2*(p2mup3*(p4mup4*F(497)+p5mup4*F(4
     &   98)+p6mup4*F(499)+p1mup4*F(500)+p2mup4*F(501))-p4mup3*p4mup4
     &   *F(502)-p4mup3*p5mup4*F(503)-p4mup3*p6mup4*F(504)-p4mup4*p5m
     &   up3*F(505)-p5mup3*p5mup4*F(506)-p5mup3*p6mup4*F(507)+p4mup4*
     &   p6mup3*F(508)+p5mup4*p6mup3*F(509)+p6mup3*p6mup4*F(510)+p2mu
     &   p4*p4mup3*F(511)+p2mup4*p5mup3*F(512)+p2mup4*p6mup3*F(513)-p
     &   1mup3*p4mup4*F(514)-p1mup4*p4mup3*F(515)-p1mup4*p5mup3*F(516
     &   )-p1mup3*p5mup4*F(516)-p1mup4*p6mup3*F(517)-p1mup3*p6mup4*F(
     &   518)-p1mup3*p2mup4*F(519)+p1mup3*p1mup4*F(520))
       Fa13(152) = mup3mup4*F(521)+2*(p2mup4*p4mup3*F(527)+p2mup4*p5mu
     &   p3*F(528)+p2mup4*p6mup3*F(529)-p2mup3*(p5mup4*F(530)+p6mup4*
     &   F(531))+p1mup3*(-(p2mup4*F(522))+p4mup4*F(524)+p5mup4*F(525)
     &   +p6mup4*F(526)-p1mup4*F(532))+p4mup3*p4mup4*F(534)+p4mup3*p5
     &   mup4*F(535)+p4mup3*p6mup4*F(536)+p4mup4*p5mup3*F(537)+p4mup4
     &   *p6mup3*F(538)+p5mup3*p5mup4*F(539)+p5mup3*p6mup4*F(540)+p5m
     &   up4*p6mup3*F(541)+p6mup3*p6mup4*F(542)+p1mup4*(p2mup3*F(523)
     &   -p4mup3*F(533)-p5mup3*F(543)+p6mup3*F(544))-p2mup3*p2mup4*F(
     &   545)+p2mup3*p4mup4*F(546))
       Fa13(153) = mup3mup4*F(553)+8*(p1mup3*(p4mup4*F(554)+p5mup4*F(5
     &   55)+p6mup4*F(556)+p2mup4*F(557)+p1mup4*F(558))-p4mup3*p4mup4
     &   *F(562)-p4mup3*p5mup4*F(563)-p4mup3*p6mup4*F(564)-p4mup4*p5m
     &   up3*F(565)-p4mup4*p6mup3*F(566)-p5mup3*p5mup4*F(567)-p5mup3*
     &   p6mup4*F(568)-p5mup4*p6mup3*F(569)-p6mup3*p6mup4*F(570)+p1mu
     &   p4*p4mup3*F(572)+p1mup4*p5mup3*F(573)+p1mup4*p6mup3*F(574)+p
     &   2mup4*p4mup3*F(575)+p2mup4*p5mup3*F(576)+p2mup4*p6mup3*F(577
     &   )+p2mup3*(p4mup4*F(559)+p5mup4*F(560)+p6mup4*F(561)-p1mup4*F
     &   (571)-p2mup4*F(578)))
       Fa13(154) = mup3mup4*F(579)+2*(p5mup3*p5mup4*F(539)+p5mup3*p6mu
     &   p4*F(540)+p5mup4*p6mup3*F(541)+p6mup3*p6mup4*F(542)+p1mup3*p
     &   5mup4*F(581)+p1mup3*p6mup4*F(582)-p1mup3*p2mup4*F(583)-p1mup
     &   3*p1mup4*F(584)+p4mup3*(p4mup4*F(534)+p5mup4*F(535)+p6mup4*F
     &   (536)+p2mup4*F(546)+p1mup4*F(585))+p1mup4*p5mup3*F(586)+p1mu
     &   p4*p6mup3*F(587)+p4mup4*(p5mup3*F(537)+p6mup3*F(538)+p1mup3*
     &   F(580)+p2mup3*F(588))+p2mup4*p5mup3*F(589)+p2mup3*p5mup4*F(5
     &   89)+p2mup3*p6mup4*F(590)+p2mup4*p6mup3*F(591)-p1mup4*p2mup3*
     &   F(592)-p2mup3*p2mup4*F(593))
       Fa13(155) = p4mup4*F(598)+p5mup4*F(599)+8*(p1mup4*F(596)+p6mup4
     &   *F(597)+p2mup4*F(600))
       Fa13(156) = p1mup3*F(601)+p2mup3*F(602)+p6mup3*F(603)+p4mup3*F(
     &   604)+p5mup3*F(605)
       Return
       End
