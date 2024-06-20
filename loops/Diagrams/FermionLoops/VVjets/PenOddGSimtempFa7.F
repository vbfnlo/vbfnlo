c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenOddGSim7(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2m
     -   up1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,
     -   p3mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mu
     -   p3,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup
     -   3,mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa7)
       IMPLICIT NONE
      Complex*16   p1mup1, p1mup2, p1mup3, p1mup4, p1mup5, p2mup1, 
     -          p2mup2, p2mup3, p2mup4, p2mup5, p3mup1, p3mup2, 
     -          p3mup3, p3mup4, p3mup5, p4mup1, p4mup2, p4mup3, 
     -          p4mup4, p4mup5, p5mup1, p5mup2, p5mup3, p5mup4, p5mup5
       Complex*16   mup1mup2, mup1mup3, mup1mup4, mup1mup5, mup2mup3, 
     -          mup2mup4, mup2mup5, mup3mup4, mup3mup5, mup4mup5
        common/FPenOddGSim/F
      Complex*16 EpsF(126)
       common/EPSFPenOddGSim/EPSF
       COMMON/InvariantsPenOddGSim/p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s3
     -   4,s45,s234,m1,m1sq,m3,m3sq,m4,m4sq
       Real*8   p1sq 
        Real*8   p2sq 
       Real*8   p3sq 
       Real*8   p4sq 
       Real*8   p5sq 
       Real*8   s12 
       Real*8   s23 
       Real*8   s34 
       Real*8   s45 
       Real*8   s234,m1,m1sq,m3,m3sq,m4,m4sq
      Complex*16 Fa(3196),F(4476)
      Real*8 P(1383) 
        COMMON/PPenOddGSim/P
       Complex*16 Fa7(595:693)
       COMMON/FaPenOddGSim/Fa
       Fa7(595) = p4mup1*(p1mup5*F(383)+2*(p3mup5*F(4465)+p4mup5*F(446
     -   7)+p2mup5*F(4469)))
       Fa7(596) = mup1mup4*(p1mup2*F(380)+2*p3mup2*F(382))
       Fa7(597) = p2mup1*p4mup3*F(1012)+p2mup3*p4mup1*F(4464)
       Fa7(598) = p4mup1*(p1mup3*F(383)+2*(p3mup3*F(4465)+p4mup3*F(446
     -   7)+p2mup3*F(4469)))
       Fa7(599) = p2mup1*p3mup2*F(1607)+2*(p1mup2*p2mup1*F(1262)-p3mup
     -   1*(p3mup2*F(382)+p1mup2*F(1610)))
       Fa7(600) = p3mup1*(p1mup2*F(380)+2*p3mup2*F(382))
       Fa7(601) = p2mup1*(p1mup2*F(380)+2*p3mup2*F(382))
       Fa7(602) = p2mup1*(p1mup2*F(383)+2*p3mup2*F(4465))
       Fa7(603) = -(mup1mup3*(p1mup2*F(383)+2*p3mup2*F(4465)))+2*mup1m
     -   up2*(p3mup3*F(4465)+p4mup3*F(4467)+p2mup3*F(4469))
       Fa7(604) = -(mup1mup3*(p1mup2*F(380)+2*p3mup2*F(382)))+2*mup1mu
     -   p2*(p4mup3*F(356)+p2mup3*F(381)+p3mup3*F(382))
       Fa7(605) = mup1mup4*(p1mup2*F(383)+2*p3mup2*F(4465))
       Fa7(606) = -(mup3mup4*p3mup1*F(396))+mup1mup4*p4mup3*F(1012)
       Fa7(607) = mup2mup3*p3mup1*F(396)+mup1mup2*(2*p2mup3*F(377)+p3m
     -   up3*F(378)+p4mup3*F(1012))-mup1mup3*(2*p1mup2*F(1262)+p3mup2
     -   *F(1607))
       Fa7(608) = (-(mup2mup3*p2mup1)+mup1mup2*p2mup3)*F(396)+2*mup1mu
     -   p3*(p3mup2*F(382)+p1mup2*F(1610))
       Fa7(609) = -(mup2mup4*p2mup1*F(396))+2*mup1mup4*(p3mup2*F(382)+
     -   p1mup2*F(1610))
       Fa7(610) = mup2mup4*p3mup1*F(396)-mup1mup4*p3mup2*F(1607)+2*mup
     -   1mup4*p1mup2*F(1608)
       Fa7(611) = -(mup2mup4*p2mup1*F(4464))+2*mup1mup4*(p1mup2*F(357)
     -   +p3mup2*F(4465))
       Fa7(612) = (-(mup2mup3*p2mup1)+mup1mup2*p2mup3)*F(4464)+2*mup1m
     -   up3*(p1mup2*F(357)+p3mup2*F(4465))
       Fa7(613) = mup1mup4*(p1mup2*F(418)+2*p3mup2*F(4466))
       Fa7(614) = p2mup1*p4mup3*F(383)+p2mup3*p4mup1*F(4464)
       Fa7(615) = p2mup1*p3mup2*F(1622)-2*(p1mup2*p3mup1*F(273)-p1mup2
     -   *p2mup1*F(1262)+p3mup1*p3mup2*F(4466))
       Fa7(616) = p3mup1*(p1mup2*F(418)+2*p3mup2*F(4466))
       Fa7(617) = p2mup1*(p1mup2*F(418)+2*p3mup2*F(4466))
       Fa7(618) = -(mup1mup3*(p1mup2*F(418)+2*p3mup2*F(4466)))+2*mup1m
     -   up2*(p4mup3*F(4465)+p3mup3*F(4466)+p2mup3*F(4468))
       Fa7(619) = mup1mup4*p4mup3*F(383)-mup3mup4*p3mup1*F(4463)
       Fa7(620) = mup1mup2*(2*p2mup3*F(379)+p4mup3*F(383)+p3mup3*F(418
     -   ))-mup1mup3*(2*p1mup2*F(1262)+p3mup2*F(1622))+mup2mup3*p3mup
     -   1*F(4463)
       Fa7(621) = (-(mup2mup3*p2mup1)+mup1mup2*p2mup3)*F(4463)+2*mup1m
     -   up3*(p1mup2*F(273)+p3mup2*F(4466))
       Fa7(622) = -(mup2mup4*p2mup1*F(4463))+2*mup1mup4*(p1mup2*F(273)
     -   +p3mup2*F(4466))
       Fa7(623) = 2*mup1mup4*p1mup2*F(1608)-mup1mup4*p3mup2*F(1622)+mu
     -   p2mup4*p3mup1*F(4463)
       Fa7(624) = p4mup1*(p1mup5*F(443)+2*(p2mup5*F(444)+p3mup5*F(445)
     -   +p4mup5*F(459)))
       Fa7(625) = p2mup1*p4mup3*F(448)+p2mup3*p4mup1*F(477)
       Fa7(626) = p4mup1*(p1mup3*F(443)+2*(p2mup3*F(444)+p3mup3*F(445)
     -   +p4mup3*F(459)))
       Fa7(627) = -(p2mup1*p3mup2*F(1607))+2*(p3mup1*p3mup2*F(382)+p1m
     -   up2*p2mup1*F(1608)+p1mup2*p3mup1*F(1610))
       Fa7(628) = p2mup1*(p1mup2*F(443)+2*p3mup2*F(445))
       Fa7(629) = mup1mup3*(p1mup2*F(443)+2*p3mup2*F(445))-2*mup1mup2*
     -   (p2mup3*F(444)+p3mup3*F(445)+p4mup3*F(459))
       Fa7(630) = mup1mup3*(p1mup2*F(380)+2*p3mup2*F(382))-2*mup1mup2*
     -   (p4mup3*F(356)+p2mup3*F(381)+p3mup3*F(382))
       Fa7(631) = mup1mup4*(p1mup2*F(443)+2*p3mup2*F(445))
       Fa7(632) = mup3mup4*p3mup1*F(396)-mup1mup4*p4mup3*F(448)
       Fa7(633) = mup2mup3*p3mup1*F(396)+mup1mup2*(2*p2mup3*F(377)+p3m
     -   up3*F(378)+p4mup3*F(448))-mup1mup3*p3mup2*F(1607)+2*mup1mup3
     -   *p1mup2*F(1608)
       Fa7(634) = mup2mup4*p2mup1*F(396)-2*mup1mup4*(p3mup2*F(382)+p1m
     -   up2*F(1610))
       Fa7(635) = -(mup2mup4*p3mup1*F(396))+2*mup1mup4*p1mup2*F(1262)+
     -   mup1mup4*p3mup2*F(1607)
       Fa7(636) = mup2mup4*p2mup1*F(477)-2*mup1mup4*(p3mup2*F(445)+p1m
     -   up2*F(1629))
       Fa7(637) = (-(mup2mup3*p2mup1)+mup1mup2*p2mup3)*F(477)+2*mup1mu
     -   p3*(p3mup2*F(445)+p1mup2*F(1629))
       Fa7(638) = -(mup1mup3*(p1mup2*F(443)+2*p3mup2*F(445)))+2*mup1mu
     -   p2*(p2mup3*F(444)+p3mup3*F(445)+p4mup3*F(459))
       Fa7(639) = -(mup3mup4*p3mup1*F(396))+mup1mup4*p4mup3*F(448)
       Fa7(640) = mup2mup3*p3mup1*F(396)+mup1mup2*(2*p2mup3*F(377)+p3m
     -   up3*F(378)+p4mup3*F(448))-mup1mup3*(2*p1mup2*F(1262)+p3mup2*
     -   F(1607))
       Fa7(641) = -(mup2mup4*p2mup1*F(477))+2*mup1mup4*(p3mup2*F(445)+
     -   p1mup2*F(1629))
       Fa7(642) = -(p2mup1*p3mup2*F(1622))+2*(p1mup2*(p3mup1*F(273)+p2
     -   mup1*F(1608))+p3mup1*p3mup2*F(4466))
       Fa7(643) = mup1mup3*(p1mup2*F(383)+2*p3mup2*F(4465))-2*mup1mup2
     -   *(p3mup3*F(4465)+p4mup3*F(4467)+p2mup3*F(4469))
       Fa7(644) = mup1mup3*(p1mup2*F(418)+2*p3mup2*F(4466))-2*mup1mup2
     -   *(p4mup3*F(4465)+p3mup3*F(4466)+p2mup3*F(4468))
       Fa7(645) = -(mup1mup4*p4mup3*F(383))+mup3mup4*p3mup1*F(4463)
       Fa7(646) = mup1mup2*(2*p2mup3*F(379)+p4mup3*F(383)+p3mup3*F(418
     -   ))+2*mup1mup3*p1mup2*F(1608)-mup1mup3*p3mup2*F(1622)+mup2mup
     -   3*p3mup1*F(4463)
       Fa7(647) = mup2mup4*p2mup1*F(4463)-2*mup1mup4*(p1mup2*F(273)+p3
     -   mup2*F(4466))
       Fa7(648) = 2*mup1mup4*p1mup2*F(1262)+mup1mup4*p3mup2*F(1622)-mu
     -   p2mup4*p3mup1*F(4463)
       Fa7(649) = mup2mup4*p2mup1*F(4464)-2*mup1mup4*(p1mup2*F(357)+p3
     -   mup2*F(4465))
       Fa7(650) = mup3mup4*p3mup1*F(396)-mup1mup4*p4mup3*F(1012)
       Fa7(651) = mup2mup3*p3mup1*F(396)+mup1mup2*(2*p2mup3*F(377)+p3m
     -   up3*F(378)+p4mup3*F(1012))-mup1mup3*p3mup2*F(1607)+2*mup1mup
     -   3*p1mup2*F(1608)
       Fa7(652) = 4*p4mup1*(p3mup3*p4mup5*F(297)+p3mup5*p4mup3*F(302)+
     -   p1mup3*p3mup5*F(344)+p3mup3*p3mup5*F(351)+p1mup5*(p3mup3*F(3
     -   48)+2*(p2mup3*F(255)+p4mup3*F(364)+p1mup3*F(634)))+p2mup3*p3
     -   mup5*F(1266)+p2mup5*p3mup3*F(1267)+2*(p2mup3*p4mup5*F(353)+p
     -   1mup3*p4mup5*F(354)+p4mup3*p4mup5*F(361)+p1mup3*p2mup5*F(638
     -   )+p2mup3*p2mup5*F(1264)+p2mup5*p4mup3*F(1277)))+mup3mup5*p4m
     -   up1*F(2343)
       Fa7(653) = -4*p4mup1*(p1mup3*(p3mup5*F(293)+p4mup5*F(297)+p1mup
     -   5*F(828)+p2mup5*F(1267))+p1mup5*(p4mup3*F(302)+p3mup3*F(310)
     -   -p2mup3*F(1291))+2*(p2mup3*p4mup5*F(304)+p3mup5*p4mup3*F(306
     -   )+p3mup3*p4mup5*F(306)+p2mup3*p3mup5*F(311)+p3mup3*p3mup5*F(
     -   312)+p4mup3*p4mup5*F(316)+p2mup3*p2mup5*F(1303)+p2mup5*p3mup
     -   3*F(1304)+p2mup5*p4mup3*F(1831)))+mup3mup5*p4mup1*F(2342)
       Fa7(654) = p4mup1*(p1mup5*F(631)+p3mup5*F(1247)+p4mup5*F(1248)+
     -   p2mup5*F(1464))
       Fa7(655) = p4mup1*(p1mup5*F(1390)+2*(p2mup5*F(1395)+p4mup5*F(13
     -   98))+p3mup5*F(1828))
       Fa7(656) = p4mup1*(p1mup5*F(603)-2*(p3mup5*F(885)+p4mup5*F(888)
     -   )+p2mup5*F(2189))
       Fa7(657) = -4*mup2mup4*(p1mup3*p4mup5*F(258)+p1mup5*p4mup3*F(26
     -   3)+p3mup5*p4mup3*F(273)+p4mup3*p4mup5*F(277)-p3mup3*p4mup5*F
     -   (282)+p2mup5*p4mup3*F(379)+p2mup3*p4mup5*F(2338))+mup2mup4*m
     -   up3mup5*F(2363)
       Fa7(658) = -4*mup2mup4*(p1mup3*(p3mup5*F(293)+p4mup5*F(297)+p1m
     -   up5*F(828)+p2mup5*F(1267))+p1mup5*(p4mup3*F(302)+p3mup3*F(31
     -   0)-p2mup3*F(1291))+2*(p2mup3*p4mup5*F(304)+p3mup5*p4mup3*F(3
     -   06)+p3mup3*p4mup5*F(306)+p2mup3*p3mup5*F(311)+p3mup3*p3mup5*
     -   F(312)+p4mup3*p4mup5*F(316)+p2mup3*p2mup5*F(1303)+p2mup5*p3m
     -   up3*F(1304)+p2mup5*p4mup3*F(1831)))+mup2mup4*mup3mup5*F(2366
     -   )
       Fa7(659) = 4*mup1mup4*(p1mup3*p4mup5*F(258)+p1mup5*p4mup3*F(263
     -   )+p3mup5*p4mup3*F(273)+p4mup3*p4mup5*F(277)-p3mup3*p4mup5*F(
     -   282)+p2mup5*p4mup3*F(379)+p2mup3*p4mup5*F(2338))+mup1mup4*mu
     -   p3mup5*F(2362)
       Fa7(660) = mup1mup4*(4*(p1mup3*(p3mup5*F(293)+p4mup5*F(297)+p1m
     -   up5*F(828)+p2mup5*F(1267))+p1mup5*(p4mup3*F(302)+p3mup3*F(31
     -   0)+p2mup3*F(1657)))+8*(p2mup3*p4mup5*F(304)+p3mup5*p4mup3*F(
     -   306)+p3mup3*p4mup5*F(306)+p2mup3*p3mup5*F(311)+p3mup3*p3mup5
     -   *F(312)+p4mup3*p4mup5*F(316)+p2mup3*p2mup5*F(1303)+p2mup5*p3
     -   mup3*F(1304)+p2mup5*p4mup3*F(1831))+mup3mup5*F(2336))
       Fa7(661) = mup1mup4*(p1mup2*F(1647)+2*p3mup2*F(1655))
       Fa7(662) = 4*mup2mup4*(p3mup3*p4mup5*F(297)+p3mup5*p4mup3*F(302
     -   )+p1mup3*p3mup5*F(344)+p3mup3*p3mup5*F(351)+p1mup5*(p3mup3*F
     -   (348)+2*(p2mup3*F(255)+p4mup3*F(364)+p1mup3*F(634)))+p2mup3*
     -   p3mup5*F(1266)+p2mup5*p3mup3*F(1267)+2*(p2mup3*p4mup5*F(353)
     -   +p1mup3*p4mup5*F(354)+p4mup3*p4mup5*F(361)+p1mup3*p2mup5*F(6
     -   38)+p2mup3*p2mup5*F(1264)+p2mup5*p4mup3*F(1277)))+mup2mup4*m
     -   up3mup5*F(2343)
       Fa7(663) = -4*mup1mup4*(p3mup3*p4mup5*F(297)+p3mup5*p4mup3*F(30
     -   2)+p1mup3*p3mup5*F(344)+p3mup3*p3mup5*F(351)+p1mup5*(p3mup3*
     -   F(348)+2*(p2mup3*F(255)+p4mup3*F(364)+p1mup3*F(634)))+p2mup3
     -   *p3mup5*F(1266)+p2mup5*p3mup3*F(1267)+2*(p2mup3*p4mup5*F(353
     -   )+p1mup3*p4mup5*F(354)+p4mup3*p4mup5*F(361)+p1mup3*p2mup5*F(
     -   638)+p2mup3*p2mup5*F(1264)+p2mup5*p4mup3*F(1277)))+mup1mup4*
     -   mup3mup5*F(2335)
       Fa7(664) = p4mup1*(p2mup3*F(625)+p3mup3*F(629)+2*(p1mup3*F(628)
     -   +p4mup3*F(986)))+p2mup1*p4mup3*F(1256)
       Fa7(665) = p4mup1*(p1mup3*F(2323)+p2mup3*F(2324)+p3mup3*F(2326)
     -   -2*p4mup3*F(2344))
       Fa7(666) = p4mup1*(p1mup3*F(601)+2*(p3mup3*F(290)+p4mup3*F(886)
     -   )+p2mup3*F(2325))+p2mup1*p4mup3*F(2352)
       Fa7(667) = p2mup1*p3mup2*F(1252)+2*(p1mup2*p2mup1*F(1258)+p1mup
     -   2*p3mup1*F(1494)+p3mup1*p3mup2*F(1498))
       Fa7(668) = p3mup1*(p1mup2*F(1278)+2*p3mup2*F(1286))
       Fa7(669) = p2mup1*(p1mup2*F(1647)+2*p3mup2*F(1655))
       Fa7(670) = mup3mup5*p4mup1*F(2305)+4*(p1mup3*(-(p2mup5*p4mup1)+
     -   p2mup1*p4mup5)*F(258)+p1mup5*(p2mup3*p4mup1+p2mup1*p4mup3)*F
     -   (263)+p2mup3*p3mup5*p4mup1*F(273)+p2mup1*p3mup5*p4mup3*F(273
     -   )+p2mup1*p4mup3*p4mup5*F(277)+p2mup5*p3mup3*p4mup1*F(282)-p2
     -   mup1*p3mup3*p4mup5*F(282)+p2mup3*p4mup1*p4mup5*F(357)+p2mup1
     -   *p2mup5*p4mup3*F(379)+p2mup3*p2mup5*p4mup1*F(1500)+p2mup1*p2
     -   mup3*p4mup5*F(2338)+p2mup5*p4mup1*p4mup3*F(2345))+mup3mup5*p
     -   2mup1*F(2362)
       Fa7(671) = -4*p2mup1*(p3mup3*p4mup5*F(297)+p3mup5*p4mup3*F(302)
     -   +p1mup3*p3mup5*F(344)+p3mup3*p3mup5*F(351)+p1mup5*(p3mup3*F(
     -   348)+2*(p2mup3*F(255)+p4mup3*F(364)+p1mup3*F(634)))+p2mup3*p
     -   3mup5*F(1266)+p2mup5*p3mup3*F(1267)+2*(p2mup3*p4mup5*F(353)+
     -   p1mup3*p4mup5*F(354)+p4mup3*p4mup5*F(361)+p1mup3*p2mup5*F(63
     -   8)+p2mup3*p2mup5*F(1264)+p2mup5*p4mup3*F(1277)))+mup3mup5*p2
     -   mup1*F(2335)
       Fa7(672) = -4*p2mup1*(p1mup3*(p3mup5*F(293)+p4mup5*F(297)+p1mup
     -   5*F(828)+p2mup5*F(1267))+p1mup5*(p4mup3*F(302)+p3mup3*F(310)
     -   -p2mup3*F(1291))+2*(p2mup3*p4mup5*F(304)+p3mup5*p4mup3*F(306
     -   )+p3mup3*p4mup5*F(306)+p2mup3*p3mup5*F(311)+p3mup3*p3mup5*F(
     -   312)+p4mup3*p4mup5*F(316)+p2mup3*p2mup5*F(1303)+p2mup5*p3mup
     -   3*F(1304)+p2mup5*p4mup3*F(1831)))+mup3mup5*p2mup1*F(2342)
       Fa7(673) = p1mup2*p2mup1*F(1388)+p2mup1*p3mup2*F(2350)+p1mup2*p
     -   3mup1*F(2355)+p3mup1*p3mup2*F(2357)
       Fa7(674) = 2*mup1mup4*p4mup3*F(2344)+mup3mup4*p3mup1*F(2354)
       Fa7(675) = mup1mup3*p1mup2*F(1697)+mup1mup2*p2mup3*F(2324)+mup1
     -   mup2*p3mup3*F(2326)-2*mup1mup2*p4mup3*F(2344)+mup1mup3*p3mup
     -   2*F(2347)+mup2mup3*p3mup1*F(2356)
       Fa7(676) = mup1mup2*(p1mup3*p1mup5*F(828)+2*(p2mup3*p3mup5*F(31
     -   1)+p4mup3*p4mup5*F(316)+p3mup3*(p3mup5*F(312)+p2mup5*F(1304)
     -   )))
       Fa7(677) = mup1mup3*p1mup2*F(1278)+2*(mup1mup2*p2mup3*F(1284)+m
     -   up1mup3*p3mup2*F(1286)+mup1mup2*p3mup3*F(1287))+mup1mup2*p4m
     -   up3*F(2340)
       Fa7(678) = mup2mup4*(4*(p2mup3*(p1mup5*F(263)+p3mup5*F(273)+p4m
     -   up5*F(357))+p2mup5*(-(p1mup3*F(258))+p3mup3*F(282)+p2mup3*F(
     -   1500)+p4mup3*F(2345)))+mup3mup5*F(2365))
       Fa7(679) = -2*mup1mup4*p4mup3*F(886)+mup3mup4*p2mup1*F(2349)
       Fa7(680) = mup1mup4*(mup3mup5*F(2187)+4*(p1mup3*p2mup5*F(258)-p
     -   2mup5*p3mup3*F(282)-p2mup3*(p1mup5*F(263)+p3mup5*F(273)+p4mu
     -   p5*F(357)+p2mup5*F(1500))-p2mup5*p4mup3*F(2345)))
       Fa7(681) = mup1mup2*(2*(p3mup3*F(290)+p4mup3*F(886))+p2mup3*F(2
     -   325))+mup1mup3*p1mup2*F(2346)+mup1mup3*p3mup2*F(2348)+mup2mu
     -   p3*p2mup1*F(2351)
       Fa7(682) = mup1mup4*p1mup2*F(2355)+mup1mup4*p3mup2*F(2357)+mup2
     -   mup4*p2mup1*F(2359)
       Fa7(683) = mup1mup4*p1mup2*F(1388)+mup1mup4*p3mup2*F(2350)+mup2
     -   mup4*p3mup1*F(2358)
       Fa7(684) = mup1mup4*p4mup3*F(1644)+mup3mup4*p3mup1*F(2361)
       Fa7(685) = mup1mup3*p3mup2*F(1252)+mup1mup2*(-2*p2mup3*F(1250)+
     -   p3mup3*F(1253)+p4mup3*F(1256))+2*mup1mup3*p1mup2*F(1258)+mup
     -   2mup3*p3mup1*F(2362)
       Fa7(686) = mup1mup4*p4mup3*F(879)+mup3mup4*p2mup1*F(2360)
       Fa7(687) = 2*mup1mup3*(p1mup2*F(1494)+p3mup2*F(1498))+mup2mup3*
     -   p2mup1*F(2341)+mup1mup2*(p4mup3*F(1215)+p2mup3*F(1492)+2*p3m
     -   up3*F(600)*P(7))
       Fa7(688) = mup1mup4*p1mup2*F(1773)+2*mup1mup4*p3mup2*F(1777)+mu
     -   p2mup4*p2mup1*F(2337)
       Fa7(689) = 2*mup1mup4*p1mup2*F(1364)+mup1mup4*p3mup2*F(1640)+mu
     -   p2mup4*p3mup1*F(2363)
       Fa7(690) = -2*mup1mup4*p4mup3*F(986)+mup3mup4*p3mup1*F(2305)+mu
     -   p3mup4*p2mup1*F(2367)
       Fa7(691) = 2*mup1mup4*p1mup2*F(1465)+mup1mup4*p3mup2*F(1466)+mu
     -   p2mup4*p2mup1*F(2319)+mup2mup4*p3mup1*F(2353)
       Fa7(692) = mup1mup2*(p2mup3*F(625)+p3mup3*F(629)+2*p4mup3*F(986
     -   ))+2*mup1mup3*p1mup2*F(1745)+mup1mup3*p3mup2*F(1746)+mup2mup
     -   3*p2mup1*F(2317)+mup2mup3*p3mup1*F(2328)
       Fa7(693) = 4*p1mup5*(p3mup2*p4mup1-p3mup1*p4mup2)*F(263)-p3mup1
     -   *(-(mup2mup5*F(2341))+4*(p3mup5*p4mup2*F(273)+p2mup5*p4mup2*
     -   F(379)+p4mup5*(p1mup2*F(1262)+p4mup2*F(1369)+p3mup2*F(2368))
     -   ))+p4mup1*(mup2mup5*F(2320)+4*(p1mup2*p3mup5*F(1262)+p3mup2*
     -   (p4mup5*F(357)+p2mup5*F(379)+p3mup5*F(1268))+p3mup5*p4mup2*F
     -   (2369)))
       Return
       End
