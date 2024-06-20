c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenOddGSim8(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2m
     -   up1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,
     -   p3mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mu
     -   p3,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup
     -   3,mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa8)
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
       Complex*16 Fa8(694:792)
       COMMON/FaPenOddGSim/Fa
       Fa8(694) = -4*p4mup1*(p1mup2*p1mup5*F(678)+p1mup2*p2mup5*F(680)
     -   +2*(p4mup2*p4mup5*F(314)+p2mup5*p4mup2*F(692))+p1mup5*p4mup2
     -   *F(948)+p3mup5*(2*(p4mup2*F(304)+p3mup2*F(311))+p1mup2*F(129
     -   2))+p3mup2*(p1mup5*F(1302)+2*(p4mup5*F(304)+p2mup5*F(1303)))
     -   +p1mup2*p4mup5*F(1399))+mup2mup5*p4mup1*F(2370)
       Fa8(695) = 4*p3mup1*(p1mup2*p1mup5*F(678)+p1mup2*p2mup5*F(680)+
     -   2*(p4mup2*p4mup5*F(314)+p2mup5*p4mup2*F(692))+p1mup5*p4mup2*
     -   F(948)+p3mup5*(2*(p4mup2*F(304)+p3mup2*F(311))+p1mup2*F(1292
     -   ))+p3mup2*(p1mup5*F(1302)+2*(p4mup5*F(304)+p2mup5*F(1303)))+
     -   p1mup2*p4mup5*F(1399))+mup2mup5*p3mup1*F(2371)
       Fa8(696) = -(mup3mup4*(-(mup2mup5*F(2360))+4*(p1mup2*p4mup5*F(1
     -   262)+p4mup2*(p1mup5*F(263)+p3mup5*F(273)+p2mup5*F(379)+p4mup
     -   5*F(1369))+p3mup2*p4mup5*F(2368))))
       Fa8(697) = -4*mup3mup4*(p1mup2*p1mup5*F(678)+p1mup2*p2mup5*F(68
     -   0)+2*(p4mup2*p4mup5*F(314)+p2mup5*p4mup2*F(692))+p1mup5*p4mu
     -   p2*F(948)+p3mup5*(2*(p4mup2*F(304)+p3mup2*F(311))+p1mup2*F(1
     -   292))+p3mup2*(p1mup5*F(1302)+2*(p4mup5*F(304)+p2mup5*F(1303)
     -   ))+p1mup2*p4mup5*F(1399))+mup2mup5*mup3mup4*F(2366)
       Fa8(698) = mup2mup4*p4mup1*F(2358)+mup1mup4*p4mup2*F(2372)
       Fa8(699) = mup1mup4*(mup2mup5*F(2362)+4*(p1mup2*p4mup5*F(1262)+
     -   p4mup2*(p1mup5*F(263)+p3mup5*F(273)+p2mup5*F(379)+p4mup5*F(1
     -   369))+p3mup2*p4mup5*F(2368)))
       Fa8(700) = 4*mup1mup4*(p1mup2*p1mup5*F(678)+p1mup2*p2mup5*F(680
     -   )+2*(p4mup2*p4mup5*F(314)+p2mup5*p4mup2*F(692))+p1mup5*p4mup
     -   2*F(948)+p3mup5*(2*(p4mup2*F(304)+p3mup2*F(311))+p1mup2*F(12
     -   92))+p3mup2*(p1mup5*F(1302)+2*(p4mup5*F(304)+p2mup5*F(1303))
     -   )+p1mup2*p4mup5*F(1399))+mup1mup4*mup2mup5*F(2371)
       Fa8(701) = mup3mup4*(mup2mup5*F(2367)+4*(p1mup2*p3mup5*F(1262)+
     -   p3mup2*(p1mup5*F(263)+p4mup5*F(357)+p2mup5*F(379)+p3mup5*F(1
     -   268))+p3mup5*p4mup2*F(2369)))
       Fa8(702) = -(mup1mup4*(-(mup2mup5*F(2319))+4*(p1mup2*p3mup5*F(1
     -   262)+p3mup2*(p1mup5*F(263)+p4mup5*F(357)+p2mup5*F(379)+p3mup
     -   5*F(1268))+p3mup5*p4mup2*F(2369))))
       Fa8(703) = p2mup1*p4mup2*F(1361)+p3mup2*p4mup1*F(1466)+2*(p1mup
     -   2*p4mup1*F(1465)+p4mup1*p4mup2*F(1469))+p3mup1*p4mup2*F(1774
     -   )
       Fa8(704) = p1mup2*p4mup1*F(1388)+2*p4mup1*p4mup2*F(1396)+p3mup2
     -   *p4mup1*F(2350)+p3mup1*p4mup2*F(2373)
       Fa8(705) = p1mup2*p4mup1*F(2355)+p3mup2*p4mup1*F(2357)+p2mup1*p
     -   4mup2*F(2372)+p4mup1*p4mup2*F(2374)
       Fa8(706) = -8*p1mup2*p4mup1*(p3mup5*F(266)+p4mup5*F(354))+4*(p1
     -   mup5*(p2mup1*p4mup2*F(263)-2*p4mup1*(p4mup2*F(364)+p1mup2*F(
     -   634)+p3mup2*F(947)))-p4mup1*(p1mup2*p2mup5*F(940)+p4mup2*(2*
     -   (p4mup5*F(361)+p3mup5*F(945))+p2mup5*F(948))+p3mup2*(2*(p4mu
     -   p5*F(945)+p3mup5*F(1025))+p2mup5*F(1302)))+p2mup1*(p3mup5*p4
     -   mup2*F(273)+p2mup5*p4mup2*F(379)+p4mup2*p4mup5*F(1369)+p3mup
     -   2*p4mup5*F(2368)))+mup2mup5*p4mup1*F(2375)
       Fa8(707) = -4*(mup1mup3*p4mup2*(p1mup5*F(263)+p3mup5*F(273)+p2m
     -   up5*F(379))-mup1mup2*(p1mup3*p4mup5*F(258)+p1mup5*p4mup3*F(2
     -   63)+p3mup5*p4mup3*F(273)-p3mup3*p4mup5*F(282)+p2mup5*p4mup3*
     -   F(379)+p2mup3*p4mup5*F(2338)))+mup1mup2*mup3mup5*F(2362)+mup
     -   1mup3*(mup2mup5*F(2341)-4*p4mup5*(p1mup2*F(1262)+p4mup2*F(13
     -   69)+p3mup2*F(2368)))
       Fa8(708) = p4mup2*(2*p2mup1*F(1704)+p3mup1*F(2376))
       Fa8(709) = 8*p1mup2*p3mup1*(p3mup5*F(266)+p4mup5*F(354))+p1mup5
     -   *(-4*p2mup1*p3mup2*F(263)+8*p3mup1*(p4mup2*F(364)+p1mup2*F(6
     -   34)+p3mup2*F(947)))-4*p2mup1*(p3mup2*(p4mup5*F(357)+p2mup5*F
     -   (379)+p3mup5*F(1268))+p3mup5*p4mup2*F(2369))+p3mup1*(8*(p3mu
     -   p2*p4mup5*F(945)+p3mup2*p3mup5*F(1025))+4*(p1mup2*p2mup5*F(9
     -   40)+p4mup2*(2*(p4mup5*F(361)+p3mup5*F(945))+p2mup5*F(948))+p
     -   2mup5*p3mup2*F(1302))+mup2mup5*F(2377))
       Fa8(710) = p2mup1*(2*p3mup5*(p4mup2*F(304)+p3mup2*F(311))+p4mup
     -   2*(2*(p4mup5*F(314)+p2mup5*F(692))+p1mup5*F(948))+p3mup2*(p1
     -   mup5*F(1302)+2*(p4mup5*F(304)+p2mup5*F(1303))))
       Fa8(711) = -4*mup1mup2*(p3mup3*p4mup5*F(297)+p3mup5*p4mup3*F(30
     -   2)+p1mup3*p3mup5*F(344)+p1mup5*(p3mup3*F(348)+2*(p2mup3*F(25
     -   5)+p4mup3*F(364)))+2*(p2mup3*p4mup5*F(353)+p1mup3*p4mup5*F(3
     -   54)+p1mup3*p2mup5*F(638)+p2mup3*p2mup5*F(1264)+p2mup5*p4mup3
     -   *F(1277)))+mup1mup2*mup3mup5*F(2335)+mup1mup3*(mup2mup5*F(23
     -   20)+4*(p1mup2*p3mup5*F(1262)+p3mup2*(p1mup5*F(263)+p4mup5*F(
     -   357)+p2mup5*F(379)+p3mup5*F(1268))+p3mup5*p4mup2*F(2369)))
       Fa8(712) = -4*mup1mup3*(p1mup2*p1mup5*F(678)+p1mup2*p2mup5*F(68
     -   0)+2*(p4mup2*p4mup5*F(314)+p2mup5*p4mup2*F(692))+p1mup5*p4mu
     -   p2*F(948)+p3mup5*(2*(p4mup2*F(304)+p3mup2*F(311))+p1mup2*F(1
     -   292))+p3mup2*(p1mup5*F(1302)+2*(p4mup5*F(304)+p2mup5*F(1303)
     -   ))+p1mup2*p4mup5*F(1399))+mup1mup2*(4*(p1mup3*(p3mup5*F(293)
     -   +p4mup5*F(297)+p2mup5*F(1267))+p1mup5*(p4mup3*F(302)+p3mup3*
     -   F(310)+p2mup3*F(1657)))+8*((p3mup5*p4mup3+p3mup3*p4mup5)*F(3
     -   06)+p2mup3*(p4mup5*F(304)+p2mup5*F(1303))+p2mup5*p4mup3*F(18
     -   31))+mup3mup5*F(2336))+mup1mup3*mup2mup5*F(2370)
       Fa8(713) = mup2mup3*p4mup1*F(2356)+mup1mup3*p4mup2*F(2373)
       Fa8(714) = -4*mup3mup4*(p1mup2*(2*(p3mup5*F(266)+p4mup5*F(354)+
     -   p1mup5*F(634))+p2mup5*F(940))+p4mup2*(2*(p4mup5*F(361)+p1mup
     -   5*F(364)+p3mup5*F(945))+p2mup5*F(948))+p3mup2*(2*(p4mup5*F(9
     -   45)+p1mup5*F(947)+p3mup5*F(1025))+p2mup5*F(1302)))+mup2mup5*
     -   mup3mup4*F(2375)
       Fa8(715) = -4*mup1mup3*(p1mup2*(2*(p3mup5*F(266)+p4mup5*F(354)+
     -   p1mup5*F(634))+p2mup5*F(940))+p4mup2*(2*(p4mup5*F(361)+p1mup
     -   5*F(364)+p3mup5*F(945))+p2mup5*F(948))+p3mup2*(2*(p4mup5*F(9
     -   45)+p1mup5*F(947)+p3mup5*F(1025))+p2mup5*F(1302)))+mup1mup2*
     -   (mup3mup5*F(2187)+4*(p1mup3*p2mup5*F(258)-p1mup5*p2mup3*F(26
     -   3)-p2mup3*p4mup5*F(357)-p2mup3*p2mup5*F(1500)-p2mup5*p4mup3*
     -   F(2345)))+mup1mup3*mup2mup5*F(2375)
       Fa8(716) = mup1mup4*(4*(p1mup2*(2*(p3mup5*F(266)+p4mup5*F(354)+
     -   p1mup5*F(634))+p2mup5*F(940))+p4mup2*(2*(p4mup5*F(361)+p1mup
     -   5*F(364)+p3mup5*F(945))+p2mup5*F(948))+p3mup2*(2*(p4mup5*F(9
     -   45)+p1mup5*F(947)+p3mup5*F(1025))+p2mup5*F(1302)))+mup2mup5*
     -   F(2377))
       Fa8(717) = mup1mup2*mup3mup4*F(2354)+mup1mup4*mup2mup3*F(2356)+
     -   mup1mup3*mup2mup4*F(2358)
       Fa8(718) = mup1mup3*p4mup2*F(1361)+mup2mup3*p4mup1*F(2317)
       Fa8(719) = mup1mup3*p4mup2*F(1774)+mup2mup3*p4mup1*F(2328)
       Fa8(720) = mup1mup4*p4mup2*F(1495)+mup2mup4*p4mup1*F(2353)
       Fa8(721) = mup1mup4*p4mup2*F(1682)+mup2mup4*p4mup1*F(2319)
       Fa8(722) = mup1mup3*mup2mup4*F(2337)+mup1mup4*mup2mup3*F(2341)+
     -   mup1mup2*mup3mup4*F(2360)
       Fa8(723) = mup1mup4*mup2mup3*F(2317)+mup1mup3*mup2mup4*F(2319)+
     -   mup1mup2*mup3mup4*F(2367)
       Fa8(724) = mup1mup2*mup3mup4*F(2305)+mup1mup4*mup2mup3*F(2328)+
     -   mup1mup3*mup2mup4*F(2353)
       Fa8(725) = mup2mup5*(p3mup1*F(396)+p4mup1*F(4464))
       Fa8(726) = -(p2mup1*p4mup2*F(1609))+2*(p3mup1*p4mup2*F(356)+p4m
     -   up1*(p1mup2*F(357)+p3mup2*F(4465)+p4mup2*F(4467)))
       Fa8(727) = 2*p3mup1*p4mup2*F(356)+p4mup1*(p1mup2*F(383)+2*(p3mu
     -   p2*F(4465)+p4mup2*F(4467)))
       Fa8(728) = -(mup1mup3*mup2mup5)+mup1mup2*mup3mup5
       Fa8(729) = mup1mup3*p4mup2*F(1609)-mup2mup3*p4mup1*F(4464)
       Fa8(730) = -(mup1mup4*p4mup2*F(1609))+mup2mup4*p4mup1*F(4464)
       Fa8(731) = mup1mup4*mup2mup3-mup1mup3*mup2mup4+mup1mup2*mup3mup
     -   4
       Fa8(732) = mup2mup5*(p3mup1*F(4463)+p4mup1*F(4464))
       Fa8(733) = -(p2mup1*p4mup2*F(1609))+2*(p1mup2*p4mup1*F(357)+(p3
     -   mup2*p4mup1+p3mup1*p4mup2)*F(4465)+p4mup1*p4mup2*F(4467))
       Fa8(734) = p1mup2*p4mup1*F(383)+2*((p3mup2*p4mup1+p3mup1*p4mup2
     -   )*F(4465)+p4mup1*p4mup2*F(4467))
       Fa8(735) = mup2mup5*(p3mup1*F(396)+p4mup1*F(477))
       Fa8(736) = -(p2mup1*p4mup2*F(1623))+2*(p3mup1*p4mup2*F(356)+p3m
     -   up2*p4mup1*F(445)+p4mup1*p4mup2*F(459)+p1mup2*p4mup1*F(1629)
     -   )
       Fa8(737) = 2*p3mup1*p4mup2*F(356)+p4mup1*(p1mup2*F(443)+2*(p3mu
     -   p2*F(445)+p4mup2*F(459)))
       Fa8(738) = mup1mup3*mup2mup5-mup1mup2*mup3mup5
       Fa8(739) = mup2mup3*p4mup1*F(477)-mup1mup3*p4mup2*F(1623)
       Fa8(740) = -(mup2mup4*p4mup1*F(477))+mup1mup4*p4mup2*F(1623)
       Fa8(741) = -(mup2mup3*p4mup1*F(477))+mup1mup3*p4mup2*F(1623)
       Fa8(742) = mup2mup4*p4mup1*F(477)-mup1mup4*p4mup2*F(1623)
       Fa8(743) = -(mup1mup3*p4mup2*F(1609))+mup2mup3*p4mup1*F(4464)
       Fa8(744) = mup1mup4*p4mup2*F(1609)-mup2mup4*p4mup1*F(4464)
       Fa8(745) = 4*p1mup5*(-(p3mup2*p4mup1)+p3mup1*p4mup2)*F(263)+p3m
     -   up1*(mup2mup5*F(2362)+4*(p3mup5*p4mup2*F(273)+p2mup5*p4mup2*
     -   F(379)+p4mup5*(p1mup2*F(1262)+p4mup2*F(1369)+p3mup2*F(2368))
     -   ))-p4mup1*(-(mup2mup5*F(2319))+4*(p1mup2*p3mup5*F(1262)+p3mu
     -   p2*(p4mup5*F(357)+p2mup5*F(379)+p3mup5*F(1268))+p3mup5*p4mup
     -   2*F(2369)))
       Fa8(746) = 4*p4mup1*(p1mup2*p1mup5*F(678)+p1mup2*p2mup5*F(680)+
     -   2*(p4mup2*p4mup5*F(314)+p2mup5*p4mup2*F(692))+p1mup5*p4mup2*
     -   F(948)+p3mup5*(2*(p4mup2*F(304)+p3mup2*F(311))+p1mup2*F(1292
     -   ))+p3mup2*(p1mup5*F(1302)+2*(p4mup5*F(304)+p2mup5*F(1303)))+
     -   p1mup2*p4mup5*F(1399))+mup2mup5*p4mup1*F(2371)
       Fa8(747) = -4*p3mup1*(p1mup2*p1mup5*F(678)+p1mup2*p2mup5*F(680)
     -   +2*(p4mup2*p4mup5*F(314)+p2mup5*p4mup2*F(692))+p1mup5*p4mup2
     -   *F(948)+p3mup5*(2*(p4mup2*F(304)+p3mup2*F(311))+p1mup2*F(129
     -   2))+p3mup2*(p1mup5*F(1302)+2*(p4mup5*F(304)+p2mup5*F(1303)))
     -   +p1mup2*p4mup5*F(1399))+mup2mup5*p3mup1*F(2370)
       Fa8(748) = mup3mup4*(mup2mup5*F(2361)+4*(p1mup2*p4mup5*F(1262)+
     -   p4mup2*(p1mup5*F(263)+p3mup5*F(273)+p2mup5*F(379)+p4mup5*F(1
     -   369))+p3mup2*p4mup5*F(2368)))
       Fa8(749) = 4*mup3mup4*(p1mup2*p1mup5*F(678)+p1mup2*p2mup5*F(680
     -   )+2*(p4mup2*p4mup5*F(314)+p2mup5*p4mup2*F(692))+p1mup5*p4mup
     -   2*F(948)+p3mup5*(2*(p4mup2*F(304)+p3mup2*F(311))+p1mup2*F(12
     -   92))+p3mup2*(p1mup5*F(1302)+2*(p4mup5*F(304)+p2mup5*F(1303))
     -   )+p1mup2*p4mup5*F(1399))+mup2mup5*mup3mup4*F(2339)
       Fa8(750) = mup2mup4*p4mup1*F(2359)+mup1mup4*p4mup2*F(2373)
       Fa8(751) = -(mup1mup4*(-(mup2mup5*F(2341))+4*(p1mup2*p4mup5*F(1
     -   262)+p4mup2*(p1mup5*F(263)+p3mup5*F(273)+p2mup5*F(379)+p4mup
     -   5*F(1369))+p3mup2*p4mup5*F(2368))))
       Fa8(752) = -4*mup1mup4*(p1mup2*p1mup5*F(678)+p1mup2*p2mup5*F(68
     -   0)+2*(p4mup2*p4mup5*F(314)+p2mup5*p4mup2*F(692))+p1mup5*p4mu
     -   p2*F(948)+p3mup5*(2*(p4mup2*F(304)+p3mup2*F(311))+p1mup2*F(1
     -   292))+p3mup2*(p1mup5*F(1302)+2*(p4mup5*F(304)+p2mup5*F(1303)
     -   ))+p1mup2*p4mup5*F(1399))+mup1mup4*mup2mup5*F(2370)
       Fa8(753) = -(mup3mup4*(-(mup2mup5*F(2364))+4*(p1mup2*p3mup5*F(1
     -   262)+p3mup2*(p1mup5*F(263)+p4mup5*F(357)+p2mup5*F(379)+p3mup
     -   5*F(1268))+p3mup5*p4mup2*F(2369))))
       Fa8(754) = mup1mup4*(mup2mup5*F(2320)+4*(p1mup2*p3mup5*F(1262)+
     -   p3mup2*(p1mup5*F(263)+p4mup5*F(357)+p2mup5*F(379)+p3mup5*F(1
     -   268))+p3mup5*p4mup2*F(2369)))
       Fa8(755) = p3mup1*p4mup2*F(1495)+p2mup1*p4mup2*F(1682)+p4mup1*(
     -   p3mup2*F(1746)+2*(p1mup2*F(1745)+p4mup2*F(1749)))
       Fa8(756) = p1mup2*p4mup1*F(1697)+2*p4mup1*p4mup2*F(1704)+p3mup2
     -   *p4mup1*F(2347)+p3mup1*p4mup2*F(2372)
       Fa8(757) = p1mup2*p4mup1*F(2346)+p3mup2*p4mup1*F(2348)+p2mup1*p
     -   4mup2*F(2373)+p4mup1*p4mup2*F(2376)
       Fa8(758) = 8*p1mup2*p4mup1*(p3mup5*F(266)+p4mup5*F(354))+p1mup5
     -   *(-4*p2mup1*p4mup2*F(263)+8*p4mup1*(p4mup2*F(364)+p1mup2*F(6
     -   34)+p3mup2*F(947)))-4*p2mup1*(p3mup5*p4mup2*F(273)+p2mup5*p4
     -   mup2*F(379)+p4mup2*p4mup5*F(1369)+p3mup2*p4mup5*F(2368))+p4m
     -   up1*(8*(p3mup2*p4mup5*F(945)+p3mup2*p3mup5*F(1025))+4*(p1mup
     -   2*p2mup5*F(940)+p4mup2*(2*(p4mup5*F(361)+p3mup5*F(945))+p2mu
     -   p5*F(948))+p2mup5*p3mup2*F(1302))+mup2mup5*F(2377))
       Fa8(759) = -4*mup1mup2*(p1mup3*p4mup5*F(258)+p1mup5*p4mup3*F(26
     -   3)+p3mup5*p4mup3*F(273)-p3mup3*p4mup5*F(282)+p2mup5*p4mup3*F
     -   (379)+p2mup3*p4mup5*F(2338))+mup1mup2*mup3mup5*F(2341)+mup1m
     -   up3*(mup2mup5*F(2362)+4*(p1mup2*p4mup5*F(1262)+p4mup2*(p1mup
     -   5*F(263)+p3mup5*F(273)+p2mup5*F(379)+p4mup5*F(1369))+p3mup2*
     -   p4mup5*F(2368)))
       Fa8(760) = p4mup2*(2*p2mup1*F(1396)+p3mup1*F(2374))
       Fa8(761) = -8*p1mup2*p3mup1*(p3mup5*F(266)+p4mup5*F(354))+4*(p1
     -   mup5*(p2mup1*p3mup2*F(263)-2*p3mup1*(p4mup2*F(364)+p1mup2*F(
     -   634)+p3mup2*F(947)))-p3mup1*(p1mup2*p2mup5*F(940)+p4mup2*(2*
     -   (p4mup5*F(361)+p3mup5*F(945))+p2mup5*F(948))+p3mup2*(2*(p4mu
     -   p5*F(945)+p3mup5*F(1025))+p2mup5*F(1302)))+p2mup1*(p3mup2*(p
     -   4mup5*F(357)+p2mup5*F(379)+p3mup5*F(1268))+p3mup5*p4mup2*F(2
     -   369)))+mup2mup5*p3mup1*F(2375)
       Fa8(762) = 4*mup1mup2*(p3mup3*p4mup5*F(297)+p3mup5*p4mup3*F(302
     -   )+p1mup3*p3mup5*F(344)+p1mup5*(p3mup3*F(348)+2*(p2mup3*F(255
     -   )+p4mup3*F(364)))+2*(p2mup3*p4mup5*F(353)+p1mup3*p4mup5*F(35
     -   4)+p1mup3*p2mup5*F(638)+p2mup3*p2mup5*F(1264)+p2mup5*p4mup3*
     -   F(1277)))+mup1mup2*mup3mup5*F(2343)-mup1mup3*(-(mup2mup5*F(2
     -   319))+4*(p1mup2*p3mup5*F(1262)+p3mup2*(p1mup5*F(263)+p4mup5*
     -   F(357)+p2mup5*F(379)+p3mup5*F(1268))+p3mup5*p4mup2*F(2369)))
       Fa8(763) = 4*(mup1mup3*(p1mup2*p1mup5*F(678)+p1mup2*p2mup5*F(68
     -   0)+2*(p4mup2*p4mup5*F(314)+p2mup5*p4mup2*F(692))+p1mup5*p4mu
     -   p2*F(948)+p3mup5*(2*(p4mup2*F(304)+p3mup2*F(311))+p1mup2*F(1
     -   292))+p3mup2*(p1mup5*F(1302)+2*(p4mup5*F(304)+p2mup5*F(1303)
     -   ))+p1mup2*p4mup5*F(1399))-mup1mup2*(p1mup3*(p3mup5*F(293)+p4
     -   mup5*F(297)+p2mup5*F(1267))+p1mup5*(p4mup3*F(302)+p3mup3*F(3
     -   10)-p2mup3*F(1291))+2*((p3mup5*p4mup3+p3mup3*p4mup5)*F(306)+
     -   p2mup3*(p4mup5*F(304)+p2mup5*F(1303))+p2mup5*p4mup3*F(1831))
     -   ))+mup1mup2*mup3mup5*F(2342)+mup1mup3*mup2mup5*F(2371)
       Fa8(764) = mup2mup3*p4mup1*F(2351)+mup1mup3*p4mup2*F(2372)
       Fa8(765) = mup3mup4*(4*(p1mup2*(2*(p3mup5*F(266)+p4mup5*F(354)+
     -   p1mup5*F(634))+p2mup5*F(940))+p4mup2*(2*(p4mup5*F(361)+p1mup
     -   5*F(364)+p3mup5*F(945))+p2mup5*F(948))+p3mup2*(2*(p4mup5*F(9
     -   45)+p1mup5*F(947)+p3mup5*F(1025))+p2mup5*F(1302)))+mup2mup5*
     -   F(2377))
       Fa8(766) = mup1mup2*(mup3mup5*F(2305)-4*(p1mup3*p2mup5*F(258)-p
     -   1mup5*p2mup3*F(263)-p2mup3*p4mup5*F(357)-p2mup3*p2mup5*F(150
     -   0)-p2mup5*p4mup3*F(2345)))+mup1mup3*(4*(p1mup2*(2*(p3mup5*F(
     -   266)+p4mup5*F(354)+p1mup5*F(634))+p2mup5*F(940))+p4mup2*(2*(
     -   p4mup5*F(361)+p1mup5*F(364)+p3mup5*F(945))+p2mup5*F(948))+p3
     -   mup2*(2*(p4mup5*F(945)+p1mup5*F(947)+p3mup5*F(1025))+p2mup5*
     -   F(1302)))+mup2mup5*F(2377))
       Fa8(767) = -4*mup1mup4*(p1mup2*(2*(p3mup5*F(266)+p4mup5*F(354)+
     -   p1mup5*F(634))+p2mup5*F(940))+p4mup2*(2*(p4mup5*F(361)+p1mup
     -   5*F(364)+p3mup5*F(945))+p2mup5*F(948))+p3mup2*(2*(p4mup5*F(9
     -   45)+p1mup5*F(947)+p3mup5*F(1025))+p2mup5*F(1302)))+mup1mup4*
     -   mup2mup5*F(2375)
       Fa8(768) = mup1mup2*mup3mup4*F(2349)+mup1mup4*mup2mup3*F(2351)+
     -   mup1mup3*mup2mup4*F(2359)
       Fa8(769) = mup1mup3*p4mup2*F(1682)+mup2mup3*p4mup1*F(2330)
       Fa8(770) = mup1mup3*p4mup2*F(1495)+mup2mup3*p4mup1*F(2334)
       Fa8(771) = mup1mup4*p4mup2*F(1774)+mup2mup4*p4mup1*F(2365)
       Fa8(772) = mup1mup4*p4mup2*F(1361)+mup2mup4*p4mup1*F(2320)
       Fa8(773) = mup1mup2*mup3mup4*F(2361)+mup1mup4*mup2mup3*F(2362)+
     -   mup1mup3*mup2mup4*F(2363)
       Fa8(774) = mup1mup3*mup2mup4*F(2320)+mup1mup4*mup2mup3*F(2330)+
     -   mup1mup2*mup3mup4*F(2364)
       Fa8(775) = mup1mup2*mup3mup4*F(2187)+mup1mup4*mup2mup3*F(2334)+
     -   mup1mup3*mup2mup4*F(2365)
       Fa8(776) = mup1mup5*(2*p1mup2*F(1364)+p3mup2*F(1640))
       Fa8(777) = mup1mup5*(p1mup2*F(1647)+2*p3mup2*F(1655))
       Fa8(778) = mup1mup5*(p1mup2*F(1388)+p3mup2*F(2350))
       Fa8(779) = mup1mup5*(2*(p2mup3*F(1250)+p1mup3*F(1251))+p3mup3*F
     -   (1641)+p4mup3*F(1644))+4*p3mup1*(p3mup5*p4mup3*F(273)+p4mup3
     -   *p4mup5*F(277)+p2mup5*p4mup3*F(379)+p2mup3*p4mup5*F(2338))+m
     -   up3mup5*p3mup1*F(2362)
       Fa8(780) = mup1mup5*(p1mup3*F(1648)-2*(p3mup3*F(1287)-p2mup3*F(
     -   1653))+p4mup3*F(2352))
       Fa8(781) = -4*p2mup3*p3mup1*(-(p1mup5*F(1291))+2*(p4mup5*F(304)
     -   +p3mup5*F(311)+p2mup5*F(1303)))-8*p3mup1*p4mup3*(p3mup5*F(30
     -   6)+p4mup5*F(316)+p2mup5*F(1831))+mup1mup5*p2mup3*F(2331)+mup
     -   3mup5*p3mup1*F(2342)+mup1mup5*(p1mup3*F(2327)+p3mup3*F(2333)
     -   +2*p4mup3*F(2344))
       Fa8(782) = mup1mup5*p3mup4*F(1678)+2*(mup1mup5*p2mup4*F(1680)+m
     -   up1mup5*p1mup4*F(1681))+mup1mup5*p4mup4*F(1683)-4*(p1mup5*p3
     -   mup4*p4mup1*F(263)+p3mup4*p4mup1*(p3mup5*F(356)+p4mup5*F(357
     -   )+p2mup5*F(379))-p2mup4*p3mup5*p4mup1*F(2318))
       Fa8(783) = -4*(p2mup4*p4mup1*(-(p1mup5*F(681))+2*(p3mup5*F(304)
     -   +p4mup5*F(314)+p2mup5*F(692)))+p3mup4*p4mup1*(-(p1mup5*F(607
     -   ))+2*(p3mup5*F(306)+p4mup5*F(316)+p2mup5*F(1831))))+mup1mup5
     -   *p1mup4*F(2307)+mup1mup5*p2mup4*F(2308)+mup1mup5*p3mup4*F(23
     -   09)+mup1mup5*p4mup4*F(2310)
       Fa8(784) = mup1mup5*(p1mup4*F(1698)+2*(p2mup4*F(1394)+p4mup4*F(
     -   1705))+p3mup4*F(2311))
       Fa8(785) = mup2mup5*(p4mup3*p4mup4*F(278)+p3mup3*p4mup4*F(283)-
     -   4*(p3mup4*p4mup3*F(274)-p1mup3*(p4mup4*F(259)+2*(p2mup4*F(25
     -   4)+p3mup4*F(256)+p1mup4*F(882)))+p1mup4*(p4mup3*F(261)+2*(p3
     -   mup3*F(279)-p2mup3*F(2379)))-p2mup3*p4mup4*F(2380)+p2mup4*p4
     -   mup3*F(2381))-8*(p2mup4*p3mup3*F(280)+p3mup3*p3mup4*F(281)-p
     -   2mup3*p3mup4*F(2378)-p2mup3*p2mup4*F(2382)))
       Fa8(786) = -(mup2mup5*(4*(p1mup4*p4mup3*F(301)+p1mup3*(p3mup4*F
     -   (294)-p4mup4*F(609)+p1mup4*F(889)+p2mup4*F(2383))-p1mup4*p2m
     -   up3*F(2384))-p3mup3*(p1mup4*F(608)+8*(p3mup4*F(307)-p4mup4*F
     -   (317)+p2mup4*F(2385)))+8*(p3mup4*p4mup3*F(309)+p4mup3*p4mup4
     -   *F(318)+p2mup3*p4mup4*F(2386)+p2mup3*p3mup4*F(2387)-p2mup4*p
     -   4mup3*F(2388)+p2mup3*p2mup4*F(2389))))
       Fa8(787) = mup1mup5*(p2mup4*F(251)-2*(p1mup4*F(253)+p3mup4*F(87
     -   8))+p4mup4*F(880))+4*(p2mup4*p4mup1*(p1mup5*F(263)+p3mup5*F(
     -   273)+p4mup5*F(357)+p2mup5*F(1470))-p2mup5*p3mup4*p4mup1*F(21
     -   88))
       Fa8(788) = mup1mup5*(p1mup4*F(285)+2*(p3mup4*F(884)+p4mup4*F(88
     -   7))+p2mup4*F(2316))
       Fa8(789) = -4*p1mup3*p3mup4*p4mup1*F(294)+p1mup4*p3mup3*p4mup1*
     -   F(608)
       Fa8(790) = 4*mup3mup4*p4mup1*(p1mup5*F(263)+p3mup5*F(273)+p4mup
     -   5*F(357)+p2mup5*F(379))+mup1mup5*mup3mup4*F(2361)
       Fa8(791) = -4*mup2mup4*p4mup1*(p1mup5*F(263)+p3mup5*F(273)+p4mu
     -   p5*F(357)+p2mup5*F(379))+mup1mup5*mup2mup4*F(2363)
       Fa8(792) = mup2mup5*(p3mup3*p4mup4*F(363)+4*(p1mup3*p3mup4*F(34
     -   5)+p3mup3*p3mup4*F(352)+p3mup4*p4mup3*F(358)+p1mup4*(p3mup3*
     -   F(347)+2*(p2mup3*F(343)-p4mup3*F(359)+p1mup3*F(938)))+p2mup3
     -   *p3mup4*F(2390)+p2mup4*p3mup3*F(2391))-8*(p4mup4*(p1mup3*F(3
     -   59)+p2mup3*F(360))-p4mup3*p4mup4*F(365)-p1mup3*p2mup4*F(2392
     -   )-p2mup3*p2mup4*F(2393)+p2mup4*p4mup3*F(2394)))
       Return
       End
