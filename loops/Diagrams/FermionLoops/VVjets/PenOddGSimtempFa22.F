c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenOddGSim22(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4
     -   ,p3mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5m
     -   up3,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mu
     -   p3,mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa22)
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
       Complex*16 Fa22(2080:2178)
       COMMON/FaPenOddGSim/Fa
       Fa22(2080) = -(p2mup1*p3mup2*F(4440))+2*(p1mup2*p2mup1*F(1608)+
     -   p3mup1*p3mup2*F(4019)+p1mup2*p3mup1*F(4443))
       Fa22(2081) = mup1mup2*(p3mup3*F(4023)+2*(p4mup3*F(4025)+p2mup3*
     -   F(4439)))-mup1mup3*(2*p1mup2*F(1262)+p3mup2*F(4440))+mup2mup
     -   3*p3mup1*F(4444)
       Fa22(2082) = -(mup1mup3*(p1mup2*F(4018)+2*p3mup2*F(4019)))+2*mu
     -   p1mup2*(p3mup3*F(4019)+p4mup3*F(4020)+p2mup3*F(4438))
       Fa22(2083) = 2*mup1mup4*(p3mup2*F(4019)+p1mup2*F(4443))-mup2mup
     -   4*p2mup1*F(4444)
       Fa22(2084) = p3mup1*(p1mup2*F(4400)-p3mup2*F(4459))
       Fa22(2085) = p3mup1*(p1mup2*F(3963)+2*p3mup2*F(4406))
       Fa22(2086) = p4mup1*(p3mup3*F(4099)+4*(p1mup3*F(3965)+p4mup3*F(
     -   3968)+p2mup3*F(4409)))
       Fa22(2087) = p4mup1*(p1mup3*F(3963)+2*(p3mup3*F(3959)+p4mup3*F(
     -   3962)+p2mup3*F(4405)))
       Fa22(2088) = mup1mup4*(p1mup2*F(4400)+p3mup2*F(4403))
       Fa22(2089) = mup1mup4*(p1mup2*F(3963)+2*p3mup2*F(4406))
       Fa22(2090) = -4*mup1mup4*p4mup3*F(3968)+2*mup3mup4*p3mup1*F(441
     -   2)
       Fa22(2091) = 4*m1sq*mup1mup4*p1mup2*F(1262)+2*(mup1mup4*p3mup2*
     -   F(4410)-mup2mup4*p3mup1*F(4412))
       Fa22(2092) = p2mup1*(p1mup2*F(4400)+p3mup2*F(4403))
       Fa22(2093) = -2*mup1mup2*p4mup3*F(3957)+4*(mup1mup2*p3mup3*F(40
     -   96)-mup1mup3*p1mup2*F(4400)+mup1mup2*p2mup3*F(4401)+mup1mup3
     -   *p3mup2*F(4459))
       Fa22(2094) = 2*p2mup1*p3mup2*F(4410)+4*(m1sq*p1mup2*p2mup1*F(12
     -   62)-p3mup1*(p1mup2*F(4411)+p3mup2*F(4413)))
       Fa22(2095) = 2*p2mup1*(p1mup2*F(3963)+2*p3mup2*F(4406))
       Fa22(2096) = 4*m1sq*mup1mup3*p1mup2*F(1608)+mup1mup2*(p3mup3*F(
     -   4099)+4*(p4mup3*F(3968)+p2mup3*F(4409)))-2*(mup1mup3*p3mup2*
     -   F(4410)-mup2mup3*p3mup1*F(4412))
       Fa22(2097) = -4*mup1mup2*(p3mup3*F(3959)+p4mup3*F(3962)+p2mup3*
     -   F(4405))+2*mup1mup3*(p1mup2*F(3963)+2*p3mup2*F(4406))
       Fa22(2098) = 2*(-(mup2mup3*p2mup1)+mup1mup2*p2mup3)*F(4412)+4*m
     -   up1mup3*(p1mup2*F(4411)+p3mup2*F(4413))
       Fa22(2099) = 2*mup2mup4*p2mup1*F(4412)-4*mup1mup4*(p1mup2*F(441
     -   1)+p3mup2*F(4413))
       Fa22(2100) = 4*m1sq*p4mup1*(p1mup2*F(1262)+p4mup2*F(2369))-2*(p
     -   3mup1*p4mup2*F(4404)-p3mup2*p4mup1*F(4410))
       Fa22(2101) = p4mup1*(p1mup2*F(3963)+2*(p3mup2*F(4406)+p4mup2*F(
     -   4408)))
       Fa22(2102) = mup1mup4*p4mup2*F(4404)-mup2mup4*p4mup1*F(4412)
       Fa22(2103) = 2*p2mup1*p4mup2*F(4404)-4*p4mup1*(p1mup2*F(4411)+p
     -   3mup2*F(4413)+p4mup2*F(4414))
       Fa22(2104) = mup1mup3*p4mup2*F(4404)-mup2mup3*p4mup1*F(4412)
       Fa22(2105) = p4mup2*(p3mup1*F(4414)+m1sq*p2mup1*F(4442))
       Fa22(2106) = p3mup1*p4mup2*F(4437)-p3mup2*p4mup1*F(4440)+2*(p1m
     -   up2*p4mup1*F(1608)+p4mup1*p4mup2*F(4442))
       Fa22(2107) = p4mup1*(p1mup2*F(4018)+2*(p3mup2*F(4019)+p4mup2*F(
     -   4020)))
       Fa22(2108) = -(mup1mup4*p4mup2*F(4437))+mup2mup4*p4mup1*F(4444)
       Fa22(2109) = -(p2mup1*p4mup2*F(4437))+2*(p3mup2*p4mup1*F(4019)+
     -   p4mup1*(p4mup2*F(4020)+p1mup2*F(4443)))
       Fa22(2110) = mup1mup3*p4mup2*F(4437)-mup2mup3*p4mup1*F(4444)
       Fa22(2111) = p4mup2*(p2mup1*F(2369)-p3mup1*F(4020))
       Fa22(2112) = 2*(p1mup2*p4mup1*F(1608)+p4mup1*p4mup2*F(4442))+p3
     -   mup1*p4mup2*F(4446)+p3mup2*p4mup1*F(4449)
       Fa22(2113) = p4mup1*(p1mup2*F(4028)+2*(p3mup2*F(4035)+p4mup2*F(
     -   4041)))
       Fa22(2114) = mup2mup4*p4mup1*F(1057)+mup1mup4*p4mup2*F(4451)
       Fa22(2115) = 2*p4mup1*(p3mup2*F(4035)+p4mup2*F(4041)+p1mup2*F(4
     -   450))+p2mup1*p4mup2*F(4451)
       Fa22(2116) = mup2mup3*p4mup1*F(1056)+mup1mup3*p4mup2*F(4446)
       Fa22(2117) = p4mup2*(p2mup1*F(2369)+p3mup1*F(4040))
       Fa22(2118) = (mup1mup4*mup2mup3+mup1mup2*mup3mup4)*F(1056)+mup1
     -   mup3*mup2mup4*F(1057)
       Fa22(2119) = 2*(p1mup2*p4mup1*F(1262)+p4mup1*p4mup2*F(2369))+p3
     -   mup1*p4mup2*F(4451)+p3mup2*p4mup1*F(4452)
       Fa22(2120) = p4mup1*(p1mup2*F(4029)+2*(p3mup2*F(4034)+p4mup2*F(
     -   4040)))
       Fa22(2121) = mup2mup4*p4mup1*F(1056)+mup1mup4*p4mup2*F(4446)
       Fa22(2122) = p2mup1*p4mup2*F(4446)+2*(p3mup2*p4mup1*F(4034)+p4m
     -   up1*(p4mup2*F(4040)+p1mup2*F(4453)))
       Fa22(2123) = mup2mup3*p4mup1*F(1057)+mup1mup3*p4mup2*F(4451)
       Fa22(2124) = p4mup2*(p3mup1*F(4041)+p2mup1*F(4442))
       Fa22(2125) = mup1mup3*mup2mup4*F(1056)+(mup1mup4*mup2mup3+mup1m
     -   up2*mup3mup4)*F(1057)
       Fa22(2126) = 2*(p1mup2*p4mup1*F(1262)+p4mup1*p4mup2*F(2369))-p3
     -   mup1*p4mup2*F(4437)+p3mup2*p4mup1*F(4440)
       Fa22(2127) = mup1mup4*p4mup2*F(4437)-mup2mup4*p4mup1*F(4444)
       Fa22(2128) = -(mup1mup3*p4mup2*F(4437))+mup2mup3*p4mup1*F(4444)
       Fa22(2129) = p4mup2*(p3mup1*F(4020)+p2mup1*F(4442))
       Fa22(2130) = 2*(p3mup1*p4mup2*F(4404)-p3mup2*p4mup1*F(4410))+4*
     -   m1sq*p4mup1*(p1mup2*F(1608)+p4mup2*F(4442))
       Fa22(2131) = 2*p4mup1*(p1mup2*F(3963)+2*(p3mup2*F(4406)+p4mup2*
     -   F(4408)))
       Fa22(2132) = -2*p2mup1*p4mup2*F(4404)+4*p4mup1*(p1mup2*F(4411)+
     -   p3mup2*F(4413)+p4mup2*F(4414))
       Fa22(2133) = p4mup2*(m1sq*p2mup1*F(2369)-p3mup1*F(4414))
       Fa22(2134) = Fa(2)+Fa(96)
       Fa22(2135) = Fa(3)+Fa(97)
       Fa22(2136) = Fa(4)+Fa(98)
       Fa22(2137) = Fa(10)+Fa(105)
       Fa22(2138) = Fa(11)+2*Fa(106)
       Fa22(2139) = Fa(13)+Fa(110)
       Fa22(2140) = Fa(14)+Fa(113)
       Fa22(2141) = 4*Fa(89)+Fa(216)
       Fa22(2142) = -4*Fa(90)+Fa(217)
       Fa22(2143) = -4*Fa(91)+Fa(218)
       Fa22(2144) = (mup3mup4*F(248))/2.+Fa(103)
       Fa22(2145) = (mup3mup4*F(249))/2.+Fa(104)
       Fa22(2146) = -4*Fa(107)+Fa(220)
       Fa22(2147) = 4*Fa(108)+Fa(221)
       Fa22(2148) = -4*Fa(109)+Fa(223)
       Fa22(2149) = (mup3mup4*F(330))/2.+Fa(111)
       Fa22(2150) = (mup3mup4*F(331))/2.+Fa(112)
       Fa22(2151) = (mup2mup4*F(322))/2.+Fa(219)
       Fa22(2152) = (mup2mup4*F(323))/2.+Fa(222)
       Fa22(2153) = (mup2mup3*F(326))/2.+Fa(224)+Fa(312)
       Fa22(2154) = (mup2mup3*F(327))/2.+Fa(225)+Fa(313)
       Fa22(2155) = (mup2mup3*F(328))/2.+Fa(226)+Fa(314)
       Fa22(2156) = (mup2mup3*F(329))/2.+Fa(228)+Fa(315)
       Fa22(2157) = (mup2mup4*F(332))/2.+Fa(231)
       Fa22(2158) = (mup2mup4*F(333))/2.+Fa(232)
       Fa22(2159) = mup1mup5*p4mup2*F(1361)+Fa(304)+Fa(693)+Fa(776)
       Fa22(2160) = mup1mup5*p4mup2*F(2373)+Fa(305)+Fa(694)+Fa(777)
       Fa22(2161) = 2*mup1mup5*p4mup2*F(1396)+Fa(306)+Fa(695)+Fa(778)
       Fa22(2162) = Fa(307)+Fa(552)+Fa(779)+4*(p1mup5*p3mup1*p4mup3*F(
     -   263)+Fa(976))
       Fa22(2163) = Fa(308)+Fa(553)+Fa(780)
       Fa22(2164) = Fa(309)+Fa(781)-4*(p1mup5*p3mup1*p4mup3*F(302)+Fa(
     -   977)+Fa(1007))
       Fa22(2165) = Fa(310)+Fa(400)+Fa(782)+Fa(978)
       Fa22(2166) = Fa(311)+Fa(402)+Fa(784)
       Fa22(2167) = 4*Fa(322)+Fa(411)+Fa(787)+Fa(980)
       Fa22(2168) = -4*Fa(323)+Fa(412)+Fa(788)
       Fa22(2169) = Fa(324)+Fa(420)+Fa(793)
       Fa22(2170) = Fa(325)+Fa(570)+Fa(796)
       Fa22(2171) = 4*p1mup2*p2mup1*p4mup5*F(1262)+mup1mup5*p4mup2*F(1
     -   774)+Fa(326)+Fa(706)+Fa(797)
       Fa22(2172) = Fa(327)+Fa(571)+Fa(799)-4*Fa(986)
       Fa22(2173) = Fa(328)+Fa(572)+Fa(800)
       Fa22(2174) = -4*p1mup2*p2mup1*p3mup5*F(1262)+2*mup1mup5*p4mup2*
     -   F(1469)+Fa(329)+Fa(709)+Fa(801)
       Fa22(2175) = mup1mup5*p4mup2*F(2374)+Fa(330)+Fa(802)-4*(Fa(710)
     -   +Fa(1012))
       Fa22(2176) = Fa(401)+Fa(783)+Fa(979)-4*Fa(1008)
       Fa22(2177) = Fa(403)+Fa(554)
       Fa22(2178) = Fa(404)+Fa(555)
       Return
       End
