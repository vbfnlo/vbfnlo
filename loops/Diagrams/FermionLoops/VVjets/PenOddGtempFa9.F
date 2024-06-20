c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenOddG9(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1
     -   ,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3m
     -   up5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,
     -   p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,m
     -   up2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa9)
       IMPLICIT NONE
      Complex*16   p1mup1, p1mup2, p1mup3, p1mup4, p1mup5, p2mup1, 
     -          p2mup2, p2mup3, p2mup4, p2mup5, p3mup1, p3mup2, 
     -          p3mup3, p3mup4, p3mup5, p4mup1, p4mup2, p4mup3, 
     -          p4mup4, p4mup5, p5mup1, p5mup2, p5mup3, p5mup4, p5mup5
       Complex*16   mup1mup2, mup1mup3, mup1mup4, mup1mup5, mup2mup3, 
     -          mup2mup4, mup2mup5, mup3mup4, mup3mup5, mup4mup5
        common/FPenOddG/F
      Complex*16 EpsF(126)
       common/EPSFPenOddG/EPSF
       COMMON/InvariantsPenOddG/p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,s
     -   45,s234,m1,m1sq,m3,m3sq,m4,m4sq
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
      Complex*16 Fa(4390),F(7459)
      Real*8 P(1989) 
        COMMON/PPenOddG/P
       Complex*16 Fa9(2337:2628)
       COMMON/FaPenOddG/Fa
       Fa9(2337) = p1mup1*(p1mup3*F(5794)-p3mup3*F(5799)+p2mup3*F(6342
     -   ))
       Fa9(2338) = p1mup1*(p1mup2*F(6168)-p2mup2*F(6174)+p3mup2*F(6343
     -   ))
       Fa9(2339) = mup1mup3*p1mup2*F(5801)+mup1mup2*p1mup3*F(5802)-4*(
     -   mup1mup3*p3mup2*F(5803)-mup1mup2*p3mup3*F(5803)+(mup1mup3*p2
     -   mup2-mup1mup2*p2mup3)*F(6345))-2*mup2mup3*p1mup1*F(6350)
       Fa9(2340) = -2*p4mup1*(-(p3mup5*F(370))+2*(p1mup5*F(276)+p2mup5
     -   *F(359)+p4mup5*F(1012)))+mup1mup5*F(5847)
       Fa9(2341) = 2*p4mup1*(2*(p4mup5*F(335)+p2mup5*F(361)+p3mup5*F(3
     -   62))+p1mup5*F(367))+mup1mup5*F(5848)
       Fa9(2342) = -2*p2mup5*p4mup1+mup1mup5*P(7)
       Fa9(2343) = 2*p4mup1*(p3mup5*F(402)+2*(p1mup5*F(276)+p2mup5*F(3
     -   59)+p4mup5*F(1012)))-mup1mup5*F(5901)
       Fa9(2344) = 2*p1mup5*p4mup1*F(402)+mup1mup5*F(5902)+4*p4mup1*(p
     -   4mup5*F(7442)+p3mup5*F(7443)+p2mup5*F(7447))
       Fa9(2345) = 2*p4mup1*(p3mup5*F(350)+2*(p1mup5*F(276)+p2mup5*F(3
     -   59)+p4mup5*F(1012)))+mup1mup5*F(5955)
       Fa9(2346) = 2*p4mup1*(2*(p4mup5*F(335)+p2mup5*F(361)+p3mup5*F(3
     -   62))+p1mup5*F(367))-mup1mup5*F(5956)
       Fa9(2347) = 2*p2mup5*p4mup1+mup1mup5*P(23)
       Fa9(2348) = -2*p4mup1*(-(p3mup5*F(370))+2*(p1mup5*F(276)+p2mup5
     -   *F(359)+p4mup5*F(1012)))+mup1mup5*F(5952)
       Fa9(2349) = 2*p4mup1*(2*(p4mup5*F(335)+p2mup5*F(361)+p3mup5*F(3
     -   62))+p1mup5*F(367))+mup1mup5*F(5953)
       Fa9(2350) = 2*p4mup1*(p3mup5*F(402)+2*(p1mup5*F(276)+p2mup5*F(3
     -   59)+p4mup5*F(1012)))+mup1mup5*F(5898)
       Fa9(2351) = 2*p1mup5*p4mup1*F(402)-mup1mup5*F(5899)+4*p4mup1*(p
     -   4mup5*F(7442)+p3mup5*F(7443)+p2mup5*F(7447))
       Fa9(2352) = 2*p4mup1*(p3mup5*F(350)+2*(p1mup5*F(276)+p2mup5*F(3
     -   59)+p4mup5*F(1012)))+mup1mup5*F(5844)
       Fa9(2353) = 2*p4mup1*(2*(p4mup5*F(335)+p2mup5*F(361)+p3mup5*F(3
     -   62))+p1mup5*F(367))-mup1mup5*F(5845)
       Fa9(2354) = mup1mup5*F(5709)-2*p4mup1*(p3mup5*F(5840)+2*(p1mup5
     -   *F(5454)+p2mup5*F(5455)+p4mup5*F(5841)))
       Fa9(2355) = mup1mup5*F(5780)+p4mup1*(p1mup5*F(5802)+4*(p3mup5*F
     -   (5803)+p4mup5*F(5804)+p2mup5*F(6345)))
       Fa9(2356) = p2mup1*p3mup3*F(5842)-4*(p1mup3*p2mup1*F(5838)+p2mu
     -   p1*p2mup3*F(6346))-2*p2mup3*p3mup1*F(6350)
       Fa9(2357) = p1mup3*p3mup1*F(5802)+p1mup1*p3mup3*F(6166)+4*(p3mu
     -   p1*p3mup3*F(5803)+p1mup1*p1mup3*F(5838)+p2mup3*p3mup1*F(6345
     -   )+p1mup1*p2mup3*F(6346))
       Fa9(2358) = p1mup3*p2mup1*F(5801)-4*p2mup1*(p3mup3*F(5803)+p2mu
     -   p3*F(6345))+2*p1mup1*p2mup3*F(6350)
       Fa9(2359) = p2mup1*p3mup2*F(6347)+p2mup2*p3mup1*F(6349)+4*(p1mu
     -   p2*p2mup1*F(5452)+p2mup1*p2mup2*F(5461)-p3mup1*p3mup2*F(5803
     -   )-p1mup2*p3mup1*F(6351))
       Fa9(2360) = p3mup1*(p1mup2*F(5801)-4*(p3mup2*F(5803)+p2mup2*F(6
     -   345)))+p1mup1*(-4*(p1mup2*F(5452)-p2mup2*F(5532))+p3mup2*F(6
     -   383))
       Fa9(2361) = p1mup2*p2mup1*F(5802)+4*(p1mup1*p3mup2*F(5803)+p2mu
     -   p1*p3mup2*F(5803)+p2mup1*p2mup2*F(6345)+p1mup1*p1mup2*F(6351
     -   ))+p1mup1*p2mup2*F(6385)
       Fa9(2362) = mup1mup5*F(5779)-2*p2mup5*p4mup1*F(6350)
       Fa9(2363) = -(mup1mup2*(p1mup3*F(5794)-p3mup3*F(5799)+p2mup3*F(
     -   6342)))+mup1mup3*(p1mup2*F(6168)-p2mup2*F(6174)+p3mup2*F(634
     -   3))
       Fa9(2364) = mup1mup2*(p3mup3*F(6166)+4*(p1mup3*F(5838)+p2mup3*F
     -   (6346)))+mup1mup3*(-4*(p1mup2*F(5452)-p2mup2*F(5532))+p3mup2
     -   *F(6383))
       Fa9(2365) = 2*mup1mup2*p2mup3*F(6350)+mup1mup3*(4*(p3mup2*F(580
     -   3)+p1mup2*F(6351))+p2mup2*F(6385))
       Fa9(2366) = mup1mup3*(p1mup2*F(5802)+4*(p3mup2*F(5803)+p2mup2*F
     -   (6345)))+mup1mup2*(p1mup3*F(5801)-4*(p3mup3*F(5803)+p2mup3*F
     -   (6345)))+2*mup2mup3*p1mup1*F(6350)
       Fa9(2367) = p4mup3*(p3mup1*F(5804)-p1mup1*F(5843))
       Fa9(2368) = p3mup1*p4mup3*F(5800)+p4mup1*(p3mup3*F(5842)-4*(p1m
     -   up3*F(5838)-p4mup3*F(5843)+p2mup3*F(6346)))
       Fa9(2369) = p4mup1*(p2mup2*F(6349)-4*(p3mup2*F(5803)+p1mup2*F(6
     -   351)))
       Fa9(2370) = p2mup1*p4mup3*F(6134)-2*p2mup3*p4mup1*F(6350)
       Fa9(2371) = 4*(p1mup2*p4mup1*F(5452)+p2mup2*p4mup1*F(5461))+p3m
     -   up1*p4mup2*F(6344)+p3mup2*p4mup1*F(6347)
       Fa9(2372) = 4*mup1mup2*p4mup3*F(5843)-2*mup2mup3*p3mup1*F(6350)
       Fa9(2373) = p1mup1*p4mup3*F(5800)+p4mup1*(p1mup3*F(5802)+4*(p3m
     -   up3*F(5803)+p4mup3*F(5804)+p2mup3*F(6345)))
       Fa9(2374) = p4mup1*(p1mup2*F(5801)-4*(p3mup2*F(5803)+p2mup2*F(6
     -   345)))
       Fa9(2375) = p4mup3*(p3mup1*F(335)+p1mup1*F(1012))
       Fa9(2376) = p4mup1*(p3mup3*F(351)+2*(p1mup3*F(349)+p4mup3*F(101
     -   2)+p2mup3*F(1050)))-p3mup1*p4mup3*F(1051)
       Fa9(2377) = p4mup1*(p2mup2*F(1624)+2*(p3mup2*F(362)+p1mup2*F(16
     -   25)))
       Fa9(2378) = p2mup3*p4mup1+p2mup1*p4mup3
       Fa9(2379) = p3mup1*p4mup2*F(1051)+p4mup1*(-(p3mup2*F(1621))+2*(
     -   p2mup2*F(1335)+p1mup2*F(1646)))
       Fa9(2380) = 2*mup1mup2*p4mup3*F(1012)+mup2mup3*p3mup1*F(1051)
       Fa9(2381) = p4mup1*(2*(p4mup3*F(335)+p2mup3*F(361)+p3mup3*F(362
     -   ))+p1mup3*F(367))+p1mup1*p4mup3*F(1051)
       Fa9(2382) = p4mup1*(2*(p2mup2*F(361)+p3mup2*F(362))+p1mup2*F(36
     -   7))
       Fa9(2383) = p4mup3*(p3mup1*F(7442)+p1mup1*F(7446))
       Fa9(2384) = -(p3mup1*p4mup3*F(7440))+p4mup1*(p3mup3*F(402)+2*(p
     -   1mup3*F(276)+p2mup3*F(359)+p4mup3*F(7446)))
       Fa9(2385) = p4mup1*(p2mup2*F(1644)+2*(p1mup2*F(888)+p3mup2*F(74
     -   43)))
       Fa9(2386) = -(p3mup2*p4mup1*F(1643))+2*(p2mup2*p4mup1*F(1335)+p
     -   1mup2*p4mup1*F(1646))+p3mup1*p4mup2*F(7440)
       Fa9(2387) = mup2mup3*p3mup1*F(7440)+2*mup1mup2*p4mup3*F(7446)
       Fa9(2388) = p1mup3*p4mup1*F(402)+p1mup1*p4mup3*F(7440)+2*p4mup1
     -   *(p4mup3*F(7442)+p3mup3*F(7443)+p2mup3*F(7447))
       Fa9(2389) = p4mup1*(p1mup2*F(402)+2*(p3mup2*F(7443)+p2mup2*F(74
     -   47)))
       Fa9(2390) = -(p3mup1*p4mup2*F(1051))+p4mup1*(p3mup2*F(1621)+2*(
     -   p1mup2*F(1325)+p2mup2*F(1623)))
       Fa9(2391) = 2*(p1mup2*p4mup1*F(1325)+p2mup2*p4mup1*F(1623))+p3m
     -   up2*p4mup1*F(1643)-p3mup1*p4mup2*F(7440)
       Fa9(2392) = p4mup3*(p3mup1*F(5804)+p1mup1*F(6167))
       Fa9(2393) = p3mup1*p4mup3*F(6134)+p4mup1*(p3mup3*F(6166)+4*(p1m
     -   up3*F(5838)+p4mup3*F(6167)+p2mup3*F(6346)))
       Fa9(2394) = p4mup1*(4*(p3mup2*F(5803)+p1mup2*F(6351))+p2mup2*F(
     -   6385))
       Fa9(2395) = p2mup1*p4mup3*F(5800)+2*p2mup3*p4mup1*F(6350)
       Fa9(2396) = -4*(p1mup2*p4mup1*F(5452)-p2mup2*p4mup1*F(5532))+p3
     -   mup1*p4mup2*F(6382)+p3mup2*p4mup1*F(6383)
       Fa9(2397) = 4*mup1mup2*p4mup3*F(6167)+2*mup2mup3*p3mup1*F(6350)
       Fa9(2398) = p1mup3*p4mup1*F(5801)+p1mup1*p4mup3*F(6134)-4*(p3mu
     -   p3*p4mup1*F(5803)+p4mup1*p4mup3*F(5804)+p2mup3*p4mup1*F(6345
     -   ))
       Fa9(2399) = p4mup1*(p1mup2*F(5802)+4*(p3mup2*F(5803)+p2mup2*F(6
     -   345)))
       Fa9(2400) = -4*p4mup1*p4mup2*F(5804)+p2mup1*p4mup2*F(6382)
       Fa9(2401) = mup1mup3*p4mup2*F(6344)+2*mup2mup3*p4mup1*F(6350)
       Fa9(2402) = -4*p4mup1*p4mup2*F(5804)+p1mup1*p4mup2*F(6344)
       Fa9(2403) = p4mup2*(2*p4mup1*F(335)-p2mup1*F(1051))
       Fa9(2404) = -(mup2mup3*p4mup1)+mup1mup3*p4mup2
       Fa9(2405) = p4mup2*(2*p4mup1*F(335)+p1mup1*F(1051))
       Fa9(2406) = p4mup2*(p2mup1*F(7440)-2*p4mup1*F(7442))
       Fa9(2407) = p4mup2*(p1mup1*F(7440)+2*p4mup1*F(7442))
       Fa9(2408) = -2*p4mup1*p4mup2*F(335)+p2mup1*p4mup2*F(1051)
       Fa9(2409) = mup2mup3*p4mup1-mup1mup3*p4mup2
       Fa9(2410) = p4mup2*(4*p4mup1*F(5804)+p2mup1*F(6344))
       Fa9(2411) = -2*mup2mup3*p4mup1*F(6350)+mup1mup3*p4mup2*F(6382)
       Fa9(2412) = p4mup2*(4*p4mup1*F(5804)+p1mup1*F(6382))
       Fa9(2413) = p1mup1*(p3mup5*F(5797)+p2mup5*F(6232))
       Fa9(2414) = 2*p1mup1*(2*p2mup5*F(5455)+p3mup5*F(5840))
       Fa9(2415) = mup1mup5*(p4mup3*F(6134)-4*(p1mup3*F(5794)-p3mup3*F
     -   (5799)+p2mup3*F(6342)))
       Fa9(2416) = mup1mup5*(p3mup3*F(6166)+4*(p1mup3*F(5838)+p4mup3*F
     -   (6167)+p2mup3*F(6346)))-2*mup3mup5*p3mup1*F(6350)
       Fa9(2417) = mup1mup5*(p1mup2*F(6168)-p2mup2*F(6174)+p3mup2*F(63
     -   43))
       Fa9(2418) = -2*mup2mup5*p2mup1*F(6350)+mup1mup5*(4*(p3mup2*F(58
     -   03)+p1mup2*F(6351))+p2mup2*F(6385))
       Fa9(2419) = mup1mup5*(p1mup3*F(5801)-4*(p3mup3*F(5803)+p4mup3*F
     -   (5804)+p2mup3*F(6345)))-2*mup3mup5*p1mup1*F(6350)
       Fa9(2420) = mup1mup5*(p1mup2*F(5802)+4*(p3mup2*F(5803)+p2mup2*F
     -   (6345)))+2*mup2mup5*p1mup1*F(6350)
       Fa9(2421) = mup1mup5*(-4*(p1mup2*F(5452)-p2mup2*F(5532))+p3mup2
     -   *F(6383))
       Fa9(2422) = mup1mup5*(p3mup3*F(351)+2*(p1mup3*F(349)+p4mup3*F(1
     -   012)+p2mup3*F(1050)))-mup3mup5*p3mup1*F(1051)
       Fa9(2423) = -(mup3mup5*p3mup1*F(7440))+mup1mup5*(p3mup3*F(402)+
     -   2*(p1mup3*F(276)+p2mup3*F(359)+p4mup3*F(7446)))
       Fa9(2424) = p1mup1*(p3mup5*F(6133)+p2mup5*F(6169))
       Fa9(2425) = p1mup1*(2*p2mup5*F(5455)+p3mup5*F(5840))
       Fa9(2426) = mup1mup5*(p4mup3*F(5800)+4*(p1mup3*F(5794)-p3mup3*F
     -   (5799)+p2mup3*F(6342)))
       Fa9(2427) = mup1mup5*(p3mup3*F(5842)-4*(p1mup3*F(5838)-p4mup3*F
     -   (5843)+p2mup3*F(6346)))+2*mup3mup5*p3mup1*F(6350)
       Fa9(2428) = 2*mup2mup5*p2mup1*F(6350)+mup1mup5*(p2mup2*F(6349)-
     -   4*(p3mup2*F(5803)+p1mup2*F(6351)))
       Fa9(2429) = mup1mup5*(p1mup3*F(5802)+4*(p3mup3*F(5803)+p4mup3*F
     -   (5804)+p2mup3*F(6345)))+2*mup3mup5*p1mup1*F(6350)
       Fa9(2430) = mup1mup5*(p1mup2*F(5801)-4*(p3mup2*F(5803)+p2mup2*F
     -   (6345)))-2*mup2mup5*p1mup1*F(6350)
       Fa9(2431) = mup1mup5*(4*(p1mup2*F(5452)+p2mup2*F(5461))+p3mup2*
     -   F(6347))
       Fa9(2432) = p2mup1*F(6386)+p3mup1*F(6392)+p4mup1*F(6393)
       Fa9(2433) = p3mup1*F(6387)+p4mup1*F(6388)+p1mup1*F(6389)
       Fa9(2434) = p2mup1*F(6390)+p1mup1*F(6394)+p4mup1*F(6395)
       Fa9(2435) = p2mup1*F(6391)+p3mup1*F(6396)+p1mup1*F(6397)
       Fa9(2436) = -2*p1mup2*F(6405)+p2mup2*F(6406)+p3mup2*F(6407)+p4m
     -   up2*F(6408)
       Fa9(2437) = p1mup2*F(6409)+p3mup2*F(6410)+p4mup2*F(6411)-2*p2mu
     -   p2*F(6412)
       Fa9(2438) = p1mup2*F(6413)+p2mup2*F(6414)+p4mup2*F(6415)+p3mup2
     -   *F(6416)
       Fa9(2439) = p2mup2*F(6417)+p3mup2*F(6418)+p1mup2*F(6419)-2*p4mu
     -   p2*F(6420)
       Fa9(2440) = p2mup3*F(6424)+p1mup3*F(6425)+p3mup3*F(6426)+p4mup3
     -   *F(6427)
       Fa9(2441) = p1mup3*F(6428)+p2mup3*F(6429)+p3mup3*F(6430)+p4mup3
     -   *F(6431)
       Fa9(2442) = p1mup3*F(6432)+p2mup3*F(6433)+p4mup3*F(6434)+2*p3mu
     -   p3*F(6435)
       Fa9(2443) = p2mup3*F(6436)+p3mup3*F(6437)+2*p4mup3*F(6438)+p1mu
     -   p3*F(6439)
       Fa9(2444) = p2mup4*F(6440)+p3mup4*F(6441)+p1mup4*F(6442)+p4mup4
     -   *F(6443)
       Fa9(2445) = p1mup4*F(6444)+p3mup4*F(6445)+p2mup4*F(6446)+p4mup4
     -   *F(6447)
       Fa9(2446) = p1mup4*F(6448)+p2mup4*F(6449)+p3mup4*F(6450)+p4mup4
     -   *F(6451)
       Fa9(2447) = p2mup4*F(6452)+p3mup4*F(6453)+2*p4mup4*F(6454)+p1mu
     -   p4*F(6455)
       Fa9(2448) = p2mup1*F(6462)+p3mup1*F(6468)+p4mup1*F(6469)
       Fa9(2449) = p3mup1*F(6463)+p4mup1*F(6464)+p1mup1*F(6465)
       Fa9(2450) = p2mup1*F(6466)+p1mup1*F(6470)+p4mup1*F(6471)
       Fa9(2451) = p2mup1*F(6467)+p1mup1*F(6472)+p3mup1*F(6473)
       Fa9(2452) = 2*p1mup2*F(6481)+p2mup2*F(6482)+p3mup2*F(6483)+p4mu
     -   p2*F(6484)
       Fa9(2453) = p1mup2*F(6485)+p3mup2*F(6486)+p4mup2*F(6487)+2*p2mu
     -   p2*F(6488)
       Fa9(2454) = p1mup2*F(6489)+p2mup2*F(6490)+p4mup2*F(6491)+2*p3mu
     -   p2*F(6492)
       Fa9(2455) = p1mup2*F(6493)+p2mup2*F(6494)+p3mup2*F(6495)+2*p4mu
     -   p2*F(6496)
       Fa9(2456) = -(p2mup3*F(6497))-2*p1mup3*F(6498)+p4mup3*F(6499)+p
     -   3mup3*F(6500)
       Fa9(2457) = 2*p2mup3*F(6501)+p1mup3*F(6502)+p3mup3*F(6503)+p4mu
     -   p3*F(6504)
       Fa9(2458) = p2mup3*F(6505)+p1mup3*F(6506)+p4mup3*F(6507)-2*p3mu
     -   p3*F(6508)
       Fa9(2459) = p2mup3*F(6509)+p1mup3*F(6510)+p3mup3*F(6511)+2*p4mu
     -   p3*F(6512)
       Fa9(2460) = p2mup4*F(6513)+2*p1mup4*F(6514)+p3mup4*F(6515)+p4mu
     -   p4*F(6516)
       Fa9(2461) = p1mup4*F(6517)+p3mup4*F(6518)-2*p2mup4*F(6519)+p4mu
     -   p4*F(6520)
       Fa9(2462) = p2mup4*F(6521)+p1mup4*F(6522)+2*p3mup4*F(6523)+p4mu
     -   p4*F(6524)
       Fa9(2463) = p2mup4*F(6525)+p1mup4*F(6526)+p3mup4*F(6527)+2*p4mu
     -   p4*F(6528)
       Fa9(2464) = p2mup1*F(6535)+p3mup1*F(6541)+p4mup1*F(6542)
       Fa9(2465) = p3mup1*F(6536)+p4mup1*F(6537)+p1mup1*F(6538)
       Fa9(2466) = p2mup1*F(6539)+p1mup1*F(6543)+p4mup1*F(6544)
       Fa9(2467) = p2mup1*F(6540)+p1mup1*F(6545)+p3mup1*F(6546)
       Fa9(2468) = -2*p1mup2*F(6554)+p2mup2*F(6555)+p3mup2*F(6556)+p4m
     -   up2*F(6557)
       Fa9(2469) = p1mup2*F(6558)+p3mup2*F(6559)+p4mup2*F(6560)+2*p2mu
     -   p2*F(6561)
       Fa9(2470) = p1mup2*F(6562)+p2mup2*F(6563)+2*p3mup2*F(6564)+p4mu
     -   p2*F(6565)
       Fa9(2471) = p1mup2*F(6566)+p2mup2*F(6567)+p3mup2*F(6568)+2*p4mu
     -   p2*F(6569)
       Fa9(2472) = -(p2mup3*F(6573))+2*p1mup3*F(6574)+p4mup3*F(6575)+p
     -   3mup3*F(6576)
       Fa9(2473) = p1mup3*F(6577)+2*p2mup3*F(6578)+p4mup3*F(6579)+p3mu
     -   p3*F(6580)
       Fa9(2474) = p1mup3*F(6581)+p2mup3*F(6582)+p4mup3*F(6583)+2*p3mu
     -   p3*F(6584)
       Fa9(2475) = p1mup3*F(6585)+p2mup3*F(6586)+p3mup3*F(6587)+2*p4mu
     -   p3*F(6588)
       Fa9(2476) = p2mup4*F(6589)+2*p1mup4*F(6590)+p3mup4*F(6591)+p4mu
     -   p4*F(6592)
       Fa9(2477) = p1mup4*F(6593)-2*p2mup4*F(6594)+p3mup4*F(6595)+p4mu
     -   p4*F(6596)
       Fa9(2478) = p1mup4*F(6597)+p2mup4*F(6598)+2*p3mup4*F(6599)+p4mu
     -   p4*F(6600)
       Fa9(2479) = p1mup4*F(6601)+p2mup4*F(6602)+p3mup4*F(6603)+2*p4mu
     -   p4*F(6604)
       Fa9(2480) = p3mup1*F(6543)+p4mup1*F(6545)+p2mup1*F(6611)
       Fa9(2481) = p3mup1*F(6612)+p4mup1*F(6613)-p1mup1*F(6614)
       Fa9(2482) = p1mup1*F(6541)+p4mup1*F(6546)+p2mup1*F(6615)
       Fa9(2483) = p1mup1*F(6542)+p3mup1*F(6544)+p2mup1*F(6616)
       Fa9(2484) = -2*p1mup2*F(6624)-p2mup2*F(6625)+p3mup2*F(6626)+p4m
     -   up2*F(6627)
       Fa9(2485) = -(p1mup2*F(6628))+p3mup2*F(6629)+p4mup2*F(6630)-2*p
     -   2mup2*F(6631)
       Fa9(2486) = -(p1mup2*F(6632))+p2mup2*F(6633)-2*p3mup2*F(6634)+p
     -   4mup2*F(6635)
       Fa9(2487) = -(p1mup2*F(6636))+p2mup2*F(6637)+p3mup2*F(6638)-2*p
     -   4mup2*F(6639)
       Fa9(2488) = 2*p1mup3*F(6590)+p3mup3*F(6591)+p4mup3*F(6592)+p2mu
     -   p3*F(6643)
       Fa9(2489) = p1mup3*F(6644)+2*p2mup3*F(6645)+p4mup3*F(6646)+p3mu
     -   p3*F(6647)
       Fa9(2490) = p1mup3*F(6597)+2*p3mup3*F(6599)+p4mup3*F(6600)+p2mu
     -   p3*F(6648)
       Fa9(2491) = p1mup3*F(6601)+p3mup3*F(6603)+2*p4mup3*F(6604)+p2mu
     -   p3*F(6649)
       Fa9(2492) = 2*p1mup4*F(6574)+p4mup4*F(6575)+p3mup4*F(6576)+p2mu
     -   p4*F(6650)
       Fa9(2493) = p2mup4*F(6652)+M1*M3*(p1mup4*F(6651)+p3mup4*F(6653)
     -   +p4mup4*F(6654))
       Fa9(2494) = p1mup4*F(6581)+p4mup4*F(6583)+2*p3mup4*F(6584)+p2mu
     -   p4*F(6655)
       Fa9(2495) = p1mup4*F(6585)+p3mup4*F(6587)+2*p4mup4*F(6588)+p2mu
     -   p4*F(6656)
       Fa9(2496) = p3mup1*F(6541)+p4mup1*F(6542)-p2mup1*F(6614)
       Fa9(2497) = p1mup1*F(6611)+p3mup1*F(6615)+p4mup1*F(6616)
       Fa9(2498) = p1mup1*F(6543)+p4mup1*F(6544)+p2mup1*F(6612)
       Fa9(2499) = p1mup1*F(6545)+p3mup1*F(6546)+p2mup1*F(6613)
       Fa9(2500) = 2*p1mup2*F(6624)+p2mup2*F(6658)+p3mup2*F(6659)+p4mu
     -   p2*F(6660)
       Fa9(2501) = p1mup2*F(6628)+2*p2mup2*F(6631)+p3mup2*F(6661)+p4mu
     -   p2*F(6662)
       Fa9(2502) = p3mup2*F(6665)+M1*M3*(p1mup2*F(6663)+p2mup2*F(6664)
     -   +p4mup2*F(6666))
       Fa9(2503) = M1*M3*(p1mup2*F(6667)+p2mup2*F(6668)+p3mup2*F(6669)
     -   )+p4mup2*F(6670)
       Fa9(2504) = 2*p1mup3*F(6574)+p4mup3*F(6575)+p3mup3*F(6576)+p2mu
     -   p3*F(6650)
       Fa9(2505) = p2mup3*F(6652)+M1*M3*(p1mup3*F(6651)+p3mup3*F(6653)
     -   +p4mup3*F(6654))
       Fa9(2506) = p1mup3*F(6581)+p4mup3*F(6583)+2*p3mup3*F(6584)+p2mu
     -   p3*F(6655)
       Fa9(2507) = p1mup3*F(6585)+p3mup3*F(6587)+2*p4mup3*F(6588)+p2mu
     -   p3*F(6656)
       Fa9(2508) = 2*p1mup4*F(6590)+p3mup4*F(6591)+p4mup4*F(6592)+p2mu
     -   p4*F(6643)
       Fa9(2509) = p1mup4*F(6644)+2*p2mup4*F(6645)+p4mup4*F(6646)+p3mu
     -   p4*F(6647)
       Fa9(2510) = p1mup4*F(6597)+2*p3mup4*F(6599)+p4mup4*F(6600)+p2mu
     -   p4*F(6648)
       Fa9(2511) = p1mup4*F(6601)+p3mup4*F(6603)+2*p4mup4*F(6604)+p2mu
     -   p4*F(6649)
       Fa9(2512) = p2mup1*F(6538)+p3mup1*F(6543)+p4mup1*F(6545)
       Fa9(2513) = p1mup1*F(6535)+p3mup1*F(6539)+p4mup1*F(6540)
       Fa9(2514) = p2mup1*F(6536)+p1mup1*F(6541)+p4mup1*F(6546)
       Fa9(2515) = p2mup1*F(6537)+p1mup1*F(6542)+p3mup1*F(6544)
       Fa9(2516) = 2*p1mup2*F(6554)+p2mup2*F(6675)+p3mup2*F(6676)+p4mu
     -   p2*F(6677)
       Fa9(2517) = p1mup2*F(6678)+p3mup2*F(6679)+p4mup2*F(6680)+2*p2mu
     -   p2*F(6681)
       Fa9(2518) = p1mup2*F(6682)+p2mup2*F(6683)-2*p3mup2*F(6684)+p4mu
     -   p2*F(6685)
       Fa9(2519) = p1mup2*F(6686)-p2mup2*F(6687)+p3mup2*F(6688)-2*p4mu
     -   p2*F(6689)
       Fa9(2520) = p2mup3*F(6589)+2*p1mup3*F(6590)+p3mup3*F(6591)+p4mu
     -   p3*F(6592)
       Fa9(2521) = p1mup3*F(6593)-2*p2mup3*F(6594)+p3mup3*F(6595)+p4mu
     -   p3*F(6596)
       Fa9(2522) = p1mup3*F(6597)+p2mup3*F(6598)+2*p3mup3*F(6599)+p4mu
     -   p3*F(6600)
       Fa9(2523) = p1mup3*F(6601)+p2mup3*F(6602)+p3mup3*F(6603)+2*p4mu
     -   p3*F(6604)
       Fa9(2524) = -(p2mup4*F(6573))+2*p1mup4*F(6574)+p4mup4*F(6575)+p
     -   3mup4*F(6576)
       Fa9(2525) = p1mup4*F(6577)+2*p2mup4*F(6578)+p4mup4*F(6579)+p3mu
     -   p4*F(6580)
       Fa9(2526) = p1mup4*F(6581)+p2mup4*F(6582)+p4mup4*F(6583)+2*p3mu
     -   p4*F(6584)
       Fa9(2527) = p1mup4*F(6585)+p2mup4*F(6586)+p3mup4*F(6587)+2*p4mu
     -   p4*F(6588)
       Fa9(2528) = p2mup1*F(6465)+p3mup1*F(6470)+p4mup1*F(6472)
       Fa9(2529) = p1mup1*F(6462)+p3mup1*F(6466)+p4mup1*F(6467)
       Fa9(2530) = p2mup1*F(6463)+p1mup1*F(6468)+p4mup1*F(6473)
       Fa9(2531) = p2mup1*F(6464)+p1mup1*F(6469)+p3mup1*F(6471)
       Fa9(2532) = 2*p1mup2*F(6698)+p2mup2*F(6699)+p3mup2*F(6700)+p4mu
     -   p2*F(6701)
       Fa9(2533) = p1mup2*F(6702)+p3mup2*F(6703)+p4mup2*F(6704)-2*p2mu
     -   p2*F(6705)
       Fa9(2534) = p1mup2*F(6706)+p2mup2*F(6707)+p4mup2*F(6708)-2*p3mu
     -   p2*F(6709)
       Fa9(2535) = p1mup2*F(6710)+p2mup2*F(6711)-p3mup2*F(6712)+2*p4mu
     -   p2*F(6713)
       Fa9(2536) = p2mup3*F(6714)+2*p1mup3*F(6715)+p4mup3*F(6716)+p3mu
     -   p3*F(6717)
       Fa9(2537) = 2*p2mup3*F(6718)+p1mup3*F(6719)+p3mup3*F(6720)-p4mu
     -   p3*F(6721)
       Fa9(2538) = p2mup3*F(6722)+p1mup3*F(6723)+p4mup3*F(6724)+2*p3mu
     -   p3*F(6725)
       Fa9(2539) = -(p2mup3*F(6726))+p1mup3*F(6727)+p3mup3*F(6728)+2*p
     -   4mup3*F(6729)
       Fa9(2540) = p2mup4*F(6730)+2*p1mup4*F(6731)+p3mup4*F(6732)+p4mu
     -   p4*F(6733)
       Fa9(2541) = p1mup4*F(6734)+p3mup4*F(6735)+2*p2mup4*F(6736)+p4mu
     -   p4*F(6737)
       Fa9(2542) = p2mup4*F(6738)+p1mup4*F(6739)+2*p3mup4*F(6740)+p4mu
     -   p4*F(6741)
       Fa9(2543) = p2mup4*F(6742)+p1mup4*F(6743)+p3mup4*F(6744)+2*p4mu
     -   p4*F(6745)
       Fa9(2544) = p2mup1*F(6389)+p3mup1*F(6394)+p4mup1*F(6397)
       Fa9(2545) = p1mup1*F(6386)+p3mup1*F(6390)+p4mup1*F(6391)
       Fa9(2546) = p2mup1*F(6387)+p1mup1*F(6392)+p4mup1*F(6396)
       Fa9(2547) = p2mup1*F(6388)+p1mup1*F(6393)+p3mup1*F(6395)
       Fa9(2548) = p1mup2*F(6757)+p2mup2*F(6758)+p3mup2*F(6759)+p4mup2
     -   *F(6760)
       Fa9(2549) = p1mup2*F(6761)+p3mup2*F(6762)+p4mup2*F(6763)+2*p2mu
     -   p2*F(6764)
       Fa9(2550) = p1mup2*F(6765)+p2mup2*F(6766)+p4mup2*F(6767)+p3mup2
     -   *F(6768)
       Fa9(2551) = p2mup2*F(6769)+p3mup2*F(6770)+p1mup2*F(6771)+2*p4mu
     -   p2*F(6772)
       Fa9(2552) = p2mup3*F(6776)+p1mup3*F(6777)+p3mup3*F(6778)+p4mup3
     -   *F(6779)
       Fa9(2553) = p1mup3*F(6780)+p2mup3*F(6781)+p3mup3*F(6782)+p4mup3
     -   *F(6783)
       Fa9(2554) = p1mup3*F(6784)+p2mup3*F(6785)+p4mup3*F(6786)-2*p3mu
     -   p3*F(6787)
       Fa9(2555) = p2mup3*F(6788)+p3mup3*F(6789)+2*p4mup3*F(6790)+p1mu
     -   p3*F(6791)
       Fa9(2556) = p2mup4*F(6792)+p3mup4*F(6793)+p1mup4*F(6794)+p4mup4
     -   *F(6795)
       Fa9(2557) = p1mup4*F(6796)+p3mup4*F(6797)+p2mup4*F(6798)+p4mup4
     -   *F(6799)
       Fa9(2558) = p1mup4*F(6800)+p2mup4*F(6801)+p3mup4*F(6802)+p4mup4
     -   *F(6803)
       Fa9(2559) = p2mup4*F(6804)+p3mup4*F(6805)+2*p4mup4*F(6806)+p1mu
     -   p4*F(6807)
       Fa9(2560) = 8*p2mup3*p2mup4*F(6818)+p2mup3*p3mup4*F(6820)+p2mup
     -   4*p3mup3*F(6821)+4*p3mup3*p3mup4*F(6822)+p1mup3*(p1mup4*F(68
     -   14)+p3mup4*F(6816)-8*p2mup4*F(6819)+p4mup4*F(6823))+p2mup3*p
     -   4mup4*F(6824)+p3mup3*p4mup4*F(6825)+p3mup4*p4mup3*F(6826)+p1
     -   mup4*(8*p2mup3*F(6815)+p3mup3*F(6817)+p4mup3*F(6827))+p2mup4
     -   *p4mup3*F(6828)+p4mup3*p4mup4*F(6829)
       Fa9(2561) = p3mup3*p3mup4*F(6833)+8*(p2mup3*p2mup4*F(6832)+p2mu
     -   p3*p3mup4*F(6834)+p2mup4*p3mup3*F(6835))+p3mup3*p4mup4*F(683
     -   7)+p2mup4*p4mup3*F(6840)+p3mup4*p4mup3*F(6841)+4*(p2mup3*p4m
     -   up4*F(6836)+p4mup3*p4mup4*F(6843))+p1mup3*(p2mup4*F(6830)+4*
     -   p1mup4*F(6838)+p3mup4*F(6839)+p4mup4*F(6844))+p1mup4*(p2mup3
     -   *F(6831)+p3mup3*F(6842)+p4mup3*F(6845))
       Fa9(2562) = p3mup4*p4mup3*F(6852)+8*(p4mup3*p4mup4*F(6853)+p2mu
     -   p4*p4mup3*F(6854))+p2mup3*(p3mup4*F(6847)+8*(p2mup4*F(6848)+
     -   p4mup4*F(6850))+p1mup4*F(6856))+p1mup3*p3mup4*F(6857)+p3mup3
     -   *(4*(p3mup4*F(6846)+p2mup4*F(6849)+p4mup4*F(6851))+p1mup4*F(
     -   6858))+4*(p1mup3*p1mup4*F(6855)+p1mup3*p2mup4*F(6859)+p1mup3
     -   *p4mup4*F(6860))+p1mup4*p4mup3*F(6861)
       Fa9(2563) = p2mup3*p3mup4*F(6865)-4*(p1mup3*p2mup4*F(6862)-p1mu
     -   p4*p2mup3*F(6863)+p2mup3*p2mup4*F(6864)-p2mup4*p3mup3*F(6866
     -   ))+p2mup3*p4mup4*F(6867)+2*p2mup4*p4mup3*F(6868)
       Fa9(2564) = -2*p2mup4*(2*(p3mup3*F(6870)+p4mup3*F(6873))+p1mup3
     -   *F(6874))+p2mup3*(-4*(p2mup4*F(6869)-p3mup4*F(6871)-p4mup4*F
     -   (6872))+2*p1mup4*F(6874))
       Fa9(2565) = p2mup3*p3mup4*F(6877)+p2mup4*p3mup3*F(6878)-4*(p2mu
     -   p3*p2mup4*F(6875)+p1mup3*p2mup4*F(6876)+p2mup4*p4mup3*F(6879
     -   )-m1sq*p2mup3*(p1mup4*F(327)+p4mup4*F(6880)))
       Fa9(2566) = p3mup3*p4mup4*F(6890)+p3mup4*p4mup3*F(6891)+2*(p3mu
     -   p3*p3mup4*F(6882)+p1mup4*p3mup3*F(6885)+p1mup3*(2*(p1mup4+p2
     -   mup4)*F(349)+p3mup4*F(6884)+p4mup4*F(6886))+p2mup4*p3mup3*F(
     -   6888)+p2mup3*(2*(p2mup4*F(6881)+p1mup4*F(6883))+p3mup4*F(688
     -   7)+p4mup4*F(6889))-p1mup4*p4mup3*F(6892)-p2mup4*p4mup3*F(689
     -   3)+p4mup3*p4mup4*F(6894))
       Fa9(2567) = p1mup3*(2*(p1mup4*F(355)+p2mup4*F(6895)+p3mup4*F(68
     -   96))+p4mup4*F(6897))-p1mup4*p4mup3*F(6898)+2*(p1mup4*p2mup3*
     -   F(6899)+p2mup3*p4mup4*F(6901)-p2mup4*p4mup3*F(6902)+p1mup4*p
     -   3mup3*F(6903)+p3mup3*p4mup4*F(6904)-p3mup4*p4mup3*F(6905))+4
     -   *(p2mup3*p2mup4*F(358)+p3mup3*(p2mup4+p3mup4)*F(362)+p2mup3*
     -   p3mup4*F(6900)+p4mup3*p4mup4*F(7444))
       Fa9(2568) = p1mup4*p3mup3*F(6908)+p1mup3*(p3mup4*F(6907)+2*(p1m
     -   up4*F(355)+p2mup4*F(6909)-p4mup4*F(6912)))+2*(p2mup4*p3mup3*
     -   F(6911)+p1mup4*p4mup3*F(6913)+p2mup3*(p1mup4*F(6906)+p3mup4*
     -   F(6910)+2*(p2mup4*F(358)+p4mup4*F(6914)))+p3mup3*p4mup4*F(69
     -   15)+p3mup4*p4mup3*F(6916))+4*(p3mup3*p3mup4*F(362)+p4mup3*(p
     -   2mup4+p4mup4)*F(7444))
       Fa9(2569) = 2*p1mup3*p2mup4*F(271)+p2mup3*(2*(p2mup4*F(348)+p1m
     -   up4*F(6917)+p3mup4*F(6918))+p4mup4*F(6919))+p2mup4*(2*p3mup3
     -   *F(1235)-p4mup3*F(6920))
       Fa9(2570) = -(p1mup4*p2mup3*F(6921))+p1mup3*p2mup4*F(6921)+2*((
     -   p2mup4*p3mup3-p2mup3*p3mup4)*F(6922)+(p2mup4*p4mup3-p2mup3*p
     -   4mup4)*F(6923))
       Fa9(2571) = p2mup3*(p3mup4*F(6925)+2*(p2mup4*F(348)+p4mup4*F(10
     -   16)+p1mup4*F(6927)))+p2mup4*(p3mup3*F(6926)+2*(p1mup3*F(6924
     -   )+p4mup3*F(6928)))
       Fa9(2572) = 2*(p1mup3*(p3mup4*F(271)+p4mup4*F(327)+2*p2mup4*F(1
     -   325))+p1mup4*(p3mup3*F(885)+2*p2mup3*F(1646)+p4mup3*F(6927))
     -   +p2mup3*p3mup4*F(6929)+p2mup3*p4mup4*F(6930)-p2mup4*(p3mup3*
     -   F(6929)+p4mup3*F(6930)))+p3mup3*p4mup4*F(6931)+p3mup4*p4mup3
     -   *F(6932)
       Fa9(2573) = p1mup3*(2*(p3mup4*F(894)+p2mup4*F(6929))+p4mup4*F(6
     -   933))+p1mup4*(2*(p3mup3*F(1235)-p2mup3*F(6929))+p4mup3*F(693
     -   4))+2*(p2mup3*p4mup4*F(6935)+p2mup4*p4mup3*F(6936)+p3mup3*p4
     -   mup4*F(6937)+p3mup4*p4mup3*F(6938)+2*(p2mup3*p3mup4*F(6939)+
     -   p2mup4*p3mup3*F(6940)))
       Fa9(2574) = p1mup3*(2*(p4mup4*F(1016)+p2mup4*F(6930))+p3mup4*F(
     -   6941))+p1mup4*(2*(p4mup3*F(6880)-p2mup3*F(6930))+p3mup3*F(69
     -   42))+2*(p2mup3*p3mup4*F(6935)+p2mup4*p3mup3*F(6936)+p3mup3*p
     -   4mup4*F(6943)+p3mup4*p4mup3*F(6944)+2*(p2mup3*p4mup4*F(6945)
     -   +p2mup4*p4mup3*F(6946)))
       Fa9(2575) = 2*(p1mup3*p2mup4*F(271)+p1mup4*p2mup3*F(885)+p2mup3
     -   *p3mup4*F(894)+p2mup4*p3mup3*F(1235))+p2mup3*p4mup4*F(6941)+
     -   p2mup4*p4mup3*F(6942)
       Fa9(2576) = p1mup3*p2mup4*F(6931)+p1mup4*p2mup3*F(6932)+2*(p2mu
     -   p3*p3mup4*F(6937)+p2mup4*p3mup3*F(6938)+p2mup3*p4mup4*F(6943
     -   )+p2mup4*p4mup3*F(6944))
       Fa9(2577) = 2*(p1mup3*p2mup4*F(327)+p2mup3*p4mup4*F(1016)+p2mup
     -   4*p4mup3*F(6880)+p1mup4*p2mup3*F(6927))+p2mup3*p3mup4*F(6933
     -   )+p2mup4*p3mup3*F(6934)
       Fa9(2578) = 2*(p1mup4*(p3mup3*F(271)+p4mup3*F(327)+2*p2mup3*F(1
     -   325))+p1mup3*(p3mup4*F(885)+2*p2mup4*F(1646)+p4mup4*F(6927))
     -   +p2mup4*p3mup3*F(6929)-p2mup3*p3mup4*F(6929)+p2mup4*p4mup3*F
     -   (6930)-p2mup3*p4mup4*F(6930))+p3mup4*p4mup3*F(6931)+p3mup3*p
     -   4mup4*F(6932)
       Fa9(2579) = p1mup4*(2*(p3mup3*F(894)+p2mup3*F(6929))+p4mup3*F(6
     -   933))+p1mup3*(2*(p3mup4*F(1235)-p2mup4*F(6929))+p4mup4*F(693
     -   4))+2*(p2mup4*p4mup3*F(6935)+p2mup3*p4mup4*F(6936)+p3mup4*p4
     -   mup3*F(6937)+p3mup3*p4mup4*F(6938)+2*(p2mup4*p3mup3*F(6939)+
     -   p2mup3*p3mup4*F(6940)))
       Fa9(2580) = p1mup4*(2*(p4mup3*F(1016)+p2mup3*F(6930))+p3mup3*F(
     -   6941))+p1mup3*(2*(p4mup4*F(6880)-p2mup4*F(6930))+p3mup4*F(69
     -   42))+2*(p2mup4*p3mup3*F(6935)+p2mup3*p3mup4*F(6936)+p3mup4*p
     -   4mup3*F(6943)+p3mup3*p4mup4*F(6944)+2*(p2mup4*p4mup3*F(6945)
     -   +p2mup3*p4mup4*F(6946)))
       Fa9(2581) = 2*(p1mup4*p2mup3*F(271)+p1mup3*p2mup4*F(885)+p2mup4
     -   *p3mup3*F(894)+p2mup3*p3mup4*F(1235))+p2mup4*p4mup3*F(6941)+
     -   p2mup3*p4mup4*F(6942)
       Fa9(2582) = p1mup4*p2mup3*F(6931)+p1mup3*p2mup4*F(6932)+2*(p2mu
     -   p4*p3mup3*F(6937)+p2mup3*p3mup4*F(6938)+p2mup4*p4mup3*F(6943
     -   )+p2mup3*p4mup4*F(6944))
       Fa9(2583) = 2*(p1mup4*p2mup3*F(327)+p2mup4*p4mup3*F(1016)+p2mup
     -   3*p4mup4*F(6880)+p1mup3*p2mup4*F(6927))+p2mup4*p3mup3*F(6933
     -   )+p2mup3*p3mup4*F(6934)
       Fa9(2584) = p1mup4*p3mup3*F(6908)+2*(p2mup4*p3mup3*F(6911)+p1mu
     -   p4*p4mup3*F(6913)+p2mup3*(p1mup4*F(6906)+p3mup4*F(6910)+2*(p
     -   2mup4*F(358)+p4mup4*F(6914)))+p3mup3*p4mup4*F(6915)+p3mup4*p
     -   4mup3*F(6916))+p1mup3*(p3mup4*F(6907)+2*(p1mup4*F(355)+p2mup
     -   4*F(6909)+p4mup4*F(6947)))+4*(p3mup3*p3mup4*F(362)+p4mup3*(p
     -   2mup4+p4mup4)*F(7444))
       Fa9(2585) = p2mup3*(2*(p2mup4*F(348)+p1mup4*F(6917)+p3mup4*F(69
     -   18))+p4mup4*F(6919))-p2mup4*(2*(p1mup3*F(885)+p3mup3*F(894))
     -   +p4mup3*F(6920))
       Fa9(2586) = p1mup4*p2mup3*F(6921)-p1mup3*p2mup4*F(6921)-2*((p2m
     -   up4*p3mup3-p2mup3*p3mup4)*F(6922)+(p2mup4*p4mup3-p2mup3*p4mu
     -   p4)*F(6923))
       Fa9(2587) = -2*p1mup4*p2mup3*F(327)+p2mup3*(2*(p2mup4*F(348)-p4
     -   mup4*F(6880))+p3mup4*F(6925))+p2mup4*(p3mup3*F(6926)+2*(p1mu
     -   p3*F(6924)+p4mup3*F(6928)))
       Fa9(2588) = p2mup3*p2mup4*F(6952)+p2mup3*p3mup4*F(6953)+p2mup4*
     -   p3mup3*F(6954)+p3mup3*p3mup4*F(6955)+p1mup3*(8*p2mup4*F(6819
     -   )+p1mup4*F(6948)+p3mup4*F(6950)+p4mup4*F(6956))+p2mup3*p4mup
     -   4*F(6957)+p3mup3*p4mup4*F(6958)+p3mup4*p4mup3*F(6959)+p1mup4
     -   *(p2mup3*F(6949)+p3mup3*F(6951)+p4mup3*F(6960))+p2mup4*p4mup
     -   3*F(6961)+p4mup3*p4mup4*F(6962)
       Fa9(2589) = p2mup3*(-8*(p2mup4*F(6832)+p3mup4*F(6834))-4*p4mup4
     -   *F(6836)+p1mup4*F(6964))-8*p3mup3*p3mup4*F(6965)+p3mup3*p4mu
     -   p4*F(6966)+p1mup3*p3mup4*F(6967)+p2mup4*(-8*p3mup3*F(6835)+p
     -   1mup3*F(6963)+p4mup3*F(6968))-4*(p1mup3*p1mup4*F(6838)+p4mup
     -   3*p4mup4*F(6843)-p3mup4*p4mup3*F(6969))+p1mup4*p3mup3*F(6970
     -   )+p1mup3*p4mup4*F(6971)+p1mup4*p4mup3*F(6972)
       Fa9(2590) = -8*(p4mup3*p4mup4*F(6853)+p2mup4*p4mup3*F(6854))-4*
     -   (p3mup3*(p3mup4*F(6846)+p2mup4*F(6849)+p4mup4*F(6851))+p1mup
     -   3*p1mup4*F(6855))+p3mup4*p4mup3*F(6974)+p2mup3*(-8*(p2mup4*F
     -   (6848)+p4mup4*F(6850))+p3mup4*F(6973)+p1mup4*F(6975))+p1mup3
     -   *p3mup4*F(6976)+p1mup4*p3mup3*F(6977)+p1mup3*p2mup4*F(6978)+
     -   p1mup3*p4mup4*F(6979)+p1mup4*p4mup3*F(6980)
       Fa9(2591) = -2*p2mup4*p4mup3*F(6868)+4*(p1mup3*p2mup4*F(6862)+p
     -   2mup4*p3mup3*F(6983))+p2mup3*(4*(p2mup4*F(6864)+p1mup4*F(698
     -   1))+p3mup4*F(6982)+p4mup4*F(6984))
       Fa9(2592) = 4*(p2mup3*p2mup4*F(6869)+p2mup4*p3mup3*F(6870)-p2mu
     -   p3*p3mup4*F(6871)-p2mup3*p4mup4*F(6872)+p2mup4*p4mup3*F(6873
     -   ))+2*(-(p1mup4*p2mup3)+p1mup3*p2mup4)*F(6874)
       Fa9(2593) = 4*(p2mup3*p2mup4*F(6875)+p1mup3*p2mup4*F(6876)+p2mu
     -   p4*p4mup3*F(6879)+m1sq*p2mup3*(p4mup4*F(1016)+p1mup4*F(6927)
     -   ))+p2mup3*p3mup4*F(6985)+p2mup4*p3mup3*F(6986)
       Fa9(2594) = -(p1mup2*(2*(p1mup3*F(6987)+p2mup3*F(6988))+p3mup3*
     -   F(6989)))+p1mup3*p3mup2*F(6990)+p2mup3*p3mup2*F(6991)+p2mup2
     -   *p3mup3*F(6992)+2*(p1mup3*p2mup2*F(6993)+p2mup2*p2mup3*F(699
     -   4))+p3mup2*p3mup3*F(6995)
       Fa9(2595) = -4*(p1mup3*p2mup2*F(6998)+p1mup2*(p1mup3*F(6996)+p2
     -   mup3*F(6997)-p3mup3*F(7001)))+p1mup3*p3mup2*F(7002)-8*(p2mup
     -   2*p2mup3*F(6999)-p2mup3*p3mup2*F(7000)-p3mup2*p3mup3*F(7003)
     -   -p2mup2*p3mup3*F(7004))
       Fa9(2596) = -2*p1mup2*p3mup3*F(7011)+p1mup3*p3mup2*F(7012)-4*(p
     -   2mup3*p3mup2*F(7007)+p3mup2*p3mup3*F(7008)+p1mup2*p1mup3*F(7
     -   009)+p1mup2*p2mup3*F(7010)+p2mup2*(p3mup3*F(7005)+2*p2mup3*F
     -   (7006)+p1mup3*F(7013)))
       Fa9(2597) = p1mup2*p3mup4*F(7014)+p3mup2*p3mup4*F(7017)+p2mup2*
     -   p3mup4*F(7018)+p3mup2*p4mup4*F(7019)+p3mup4*p4mup2*F(7020)-4
     -   *(p3mup2*(p1mup4*F(7015)+p2mup4*F(7016))-m1sq*(p1mup2*(2*(p1
     -   mup4*F(634)+p2mup4*F(635))+p4mup4*F(1326))-p1mup4*p4mup2*F(1
     -   327)-p2mup4*p4mup2*F(1329)+p4mup2*p4mup4*F(1337)+2*(p2mup2*p
     -   2mup4*F(7021)+p1mup4*p2mup2*F(7022))+p2mup2*p4mup4*F(7023)))
       Fa9(2598) = -8*(p3mup2*p3mup4*F(7029)+p2mup4*p3mup2*F(7030))+p3
     -   mup2*p4mup4*F(7031)+p2mup4*p4mup2*F(7033)+p2mup2*(-8*(p2mup4
     -   *F(7024)+p3mup4*F(7025))+p4mup4*F(7026)+p1mup4*F(7035))-4*(p
     -   3mup4*p4mup2*F(7034)+p4mup2*p4mup4*F(7036))+p1mup4*p3mup2*F(
     -   7037)+p1mup2*(4*p1mup4*F(7027)+p2mup4*F(7028)+p3mup4*F(7032)
     -   +p4mup4*F(7038))+2*p1mup4*p4mup2*F(7039)
       Fa9(2599) = p2mup2*p3mup4*F(7040)+4*(m1sq*(p1mup4*p2mup2*F(675)
     -   +p1mup2*(p1mup4*F(667)+p2mup4*F(669)-p4mup4*F(1392))+p1mup4*
     -   p4mup2*F(1393)+2*(p2mup2*p2mup4*F(1390)+p2mup2*p4mup4*F(1396
     -   )+p4mup2*p4mup4*F(1397)+p2mup4*p4mup2*F(4112)))-p3mup2*p3mup
     -   4*F(7041)-p2mup4*p3mup2*F(7042)-p3mup2*p4mup4*F(7043))+p3mup
     -   4*p4mup2*F(7044)+p1mup2*p3mup4*F(7045)-2*p1mup4*p3mup2*F(704
     -   6)
       Fa9(2600) = p2mup2*p2mup4*F(7050)+p2mup2*p3mup4*F(7052)+p2mup4*
     -   p3mup2*F(7055)+p1mup2*(p2mup4*F(7048)+8*(p1mup4*F(7047)+p3mu
     -   p4*F(7051))+p4mup4*F(7056))+p2mup2*p4mup4*F(7057)+p2mup4*p4m
     -   up2*F(7058)+p3mup2*p4mup4*F(7059)+p4mup2*p4mup4*F(7060)-4*p1
     -   mup4*(p2mup2*F(7049)-2*(p3mup2*F(7054)-p4mup2*F(7061)))+8*(p
     -   3mup2*p3mup4*F(7053)-p3mup4*p4mup2*F(7062))
       Fa9(2601) = p2mup4*(2*(p3mup2*F(7064)+p4mup2*F(7067))+p1mup2*F(
     -   7073))+p2mup2*(2*(p2mup4*F(7063)+p3mup4*F(7065)+p4mup4*F(706
     -   6))+p1mup4*F(7074))+2*(p1mup2*p1mup4*F(7072)+p1mup4*p3mup2*F
     -   (7076)+2*(p4mup2*p4mup4*F(7068)+p3mup4*p4mup2*F(7069)+p3mup2
     -   *p3mup4*F(7070)+p3mup2*p4mup4*F(7071)+p1mup2*p3mup4*F(7075)+
     -   p1mup2*p4mup4*F(7077))+p1mup4*p4mup2*F(7078))
       Fa9(2602) = 2*p2mup4*p3mup2*F(7081)+p3mup4*(p1mup2*F(7079)+p2mu
     -   p2*F(7080)+p3mup2*F(7082)+p4mup2*F(7083))+4*m1sq*(p1mup4*p2m
     -   up2*F(963)+p2mup2*p2mup4*F(1518)+p2mup4*p4mup2*F(1522)+p1mup
     -   2*(p2mup4*F(961)+2*(p1mup4*F(960)+p4mup4*F(7085)))+p2mup2*p4
     -   mup4*F(7086)+2*(p4mup2*p4mup4*F(1018)+p1mup4*p3mup2*F(7084)+
     -   p1mup4*p4mup2*F(7085)+p3mup2*p4mup4*F(7087)))
       Fa9(2603) = p2mup3*p3mup2*F(7093)-4*(p1mup3*p2mup2*F(7090)+p2mu
     -   p2*p2mup3*F(7091)-p2mup2*p3mup3*F(7092)+p1mup2*(p2mup3*F(708
     -   9)+2*(p1mup3*F(7088)-p3mup3*F(7094))))+8*(p1mup3*p3mup2*F(70
     -   94)+p3mup2*p3mup3*F(7095))
       Fa9(2604) = p2mup2*p3mup3*F(7097)+p2mup3*p3mup2*F(7098)+2*(p1mu
     -   p3*p2mup2*F(7100)+p1mup2*(2*p1mup3*F(7096)+p2mup3*F(7099)+p3
     -   mup3*F(7101))+p2mup2*p2mup3*F(7102)+2*p1mup3*p3mup2*F(7103)+
     -   p3mup2*p3mup3*F(7104))
       Fa9(2605) = p1mup3*p2mup2*F(7074)+p1mup2*(2*p1mup3*F(7072)+p2mu
     -   p3*F(7073)+4*p3mup3*F(7075))+2*(p1mup3*p3mup2*F(7076)+p2mup2
     -   *p2mup3*F(7105)+p2mup2*p3mup3*F(7106)+p2mup3*p3mup2*F(7107))
     -   +4*p3mup2*p3mup3*F(7108)
       Fa9(2606) = p1mup2*(2*(p1mup3*F(634)+p2mup3*F(635))+p3mup3*F(18
     -   38))-p1mup3*p3mup2*F(1839)-p2mup3*p3mup2*F(1840)+p3mup2*p3mu
     -   p3*F(1849)+2*(p2mup2*p2mup3*F(7021)+p1mup3*p2mup2*F(7022))+p
     -   2mup2*p3mup3*F(7109)
       Fa9(2607) = p1mup2*(p1mup3*F(1285)+p2mup3*F(1286)-p3mup3*F(1871
     -   ))+p3mup2*(p1mup3*F(1872)+2*(p3mup3*F(1876)+p2mup3*F(4199)))
     -   +p2mup2*(p1mup3*F(1289)+2*(p2mup3*F(1288)+p3mup3*F(4200)))
       Fa9(2608) = p1mup2*(2*(p1mup3*F(7110)+p2mup3*F(7111))+p3mup3*F(
     -   7112))+p1mup3*p3mup2*F(7113)+2*(p2mup2*p3mup3*F(7114)+p2mup3
     -   *p3mup2*F(7115)+p1mup3*p2mup2*F(7116))+4*(p2mup2*p2mup3*F(71
     -   17)+p3mup2*p3mup3*F(7118))
       Fa9(2609) = 4*p2mup2*p2mup4*F(1281)-p3mup2*p4mup4*F(7122)-p3mup
     -   4*p4mup2*F(7123)+p1mup2*(-4*(p1mup4*F(634)-p2mup4*F(7124))+2
     -   *(p3mup4*F(7125)+p4mup4*F(7126)))+2*(p2mup4*p4mup2*F(1329)+p
     -   2mup2*p4mup4*F(1680)+p2mup2*p3mup4*F(2285)-p3mup2*p3mup4*F(7
     -   119)+p1mup4*(2*p2mup2*F(641)+p4mup2*F(1327)-p3mup2*F(7120))-
     -   p2mup4*p3mup2*F(7121)+p4mup2*p4mup4*F(7127))
       Fa9(2610) = p1mup2*(2*(p1mup4*F(7128)+p2mup4*F(7129)+p3mup4*F(7
     -   130))+p4mup4*F(7131))-p1mup4*p4mup2*F(7132)+2*(p1mup4*p2mup2
     -   *F(7134)+p2mup2*p4mup4*F(7136)-p2mup4*p4mup2*F(7137)+p1mup4*
     -   p3mup2*F(7139)+p3mup2*p4mup4*F(7141)-p3mup4*p4mup2*F(7142))+
     -   4*(p2mup2*p2mup4*F(7133)+p2mup2*p3mup4*F(7135)+p3mup2*p3mup4
     -   *F(7138)+p2mup4*p3mup2*F(7140)+p4mup2*p4mup4*F(7143))
       Fa9(2611) = p1mup2*p3mup4*F(7144)+4*(p2mup2*p4mup4*F(3879)+p4mu
     -   p2*p4mup4*F(3880)+p3mup2*p3mup4*F(7150))+2*(p1mup2*p4mup4*F(
     -   1392)+p2mup2*p3mup4*F(7148)+p2mup4*(2*(p4mup2*F(1391)+p2mup2
     -   *F(3878))+p1mup2*F(7147)+p3mup2*F(7149))+p3mup2*p4mup4*F(715
     -   1)+p3mup4*p4mup2*F(7152))+p1mup4*(p3mup2*F(7145)-2*(p2mup2*F
     -   (675)-p1mup2*F(7146)-p4mup2*F(7153)))
       Fa9(2612) = p2mup2*(2*(p2mup4*F(7155)+p1mup4*F(7157)+p3mup4*F(7
     -   160))+p4mup4*F(7161))+p2mup4*p4mup2*F(7162)+2*(p1mup2*(p2mup
     -   4*F(7156)+2*(p1mup4*F(7154)+p3mup4*F(7158))+p4mup4*F(7159))+
     -   p2mup4*p3mup2*F(7165)+p3mup2*p4mup4*F(7166)+p4mup2*p4mup4*F(
     -   7167))+4*(p3mup2*p3mup4*F(7163)+p1mup4*p3mup2*F(7164)+p1mup4
     -   *p4mup2*F(7168)+p3mup4*p4mup2*F(7169))
       Fa9(2613) = p1mup2*(2*p1mup4*F(7170)+p2mup4*F(7171)+4*(p3mup4*F
     -   (7174)+p4mup4*F(7177)))+p1mup4*(p2mup2*F(7172)+2*(p3mup2*F(7
     -   175)+p4mup2*F(7178)))+2*(p2mup4*p3mup2*F(7176)+p2mup2*p3mup4
     -   *F(7176)+p2mup4*p4mup2*F(7179)+p2mup2*p4mup4*F(7179)+2*(p2mu
     -   p2*p2mup4*F(7173)+p3mup4*p4mup2*F(7180)+p3mup2*p4mup4*F(7180
     -   )+p3mup2*p3mup4*F(7181)+p4mup2*p4mup4*F(7182)))
       Fa9(2614) = p2mup2*p3mup4*F(7184)-p2mup4*p3mup2*F(7185)-2*(p2mu
     -   p2*p2mup4*F(1518)-p2mup2*p4mup4*F(1777)-p1mup2*(2*p4mup4*F(3
     -   36)-p2mup4*F(961)+p3mup4*F(7183))-p1mup4*(2*(p4mup2*F(336)-p
     -   1mup2*F(960)+p3mup2*F(1520))+p2mup2*F(7186))-p3mup2*p3mup4*F
     -   (7187)-p3mup4*p4mup2*F(7188)-p2mup4*p4mup2*F(7189))+4*(p3mup
     -   2*p4mup4*F(1523)+p4mup2*p4mup4*F(7190))
       Fa9(2615) = p1mup2*p2mup3*F(1551)+p1mup3*p2mup2*F(1552)+p2mup2*
     -   p2mup3*F(1553)+p2mup3*p3mup2*F(1940)+p2mup2*p3mup3*F(7191)+2
     -   *(p1mup2*p1mup3*F(886)+p3mup2*p3mup3*F(1941)+(p1mup3*p3mup2+
     -   p1mup2*p3mup3)*F(7192))
       Fa9(2616) = p2mup2*p3mup3*F(7194)+p2mup3*p3mup2*F(7195)+2*(p1mu
     -   p3*p2mup2*F(7197)+p1mup2*(2*p1mup3*F(7193)+p2mup3*F(7196)+p3
     -   mup3*F(7198))+p2mup2*p2mup3*F(7199)+2*p1mup3*p3mup2*F(7200)+
     -   p3mup2*p3mup3*F(7201))
       Fa9(2617) = p1mup3*p2mup2*F(7172)+p1mup2*(2*p1mup3*F(7170)+p2mu
     -   p3*F(7171)+4*p3mup3*F(7174))+2*(p1mup3*p3mup2*F(7175)+p2mup3
     -   *p3mup2*F(7176)+p2mup2*p3mup3*F(7176))+4*(p2mup2*p2mup3*F(71
     -   73)+p3mup2*p3mup3*F(7181))
       Fa9(2618) = p1mup2*(2*(p1mup3*F(667)+p2mup3*F(669))+p3mup3*F(72
     -   02))+p1mup3*p3mup2*F(7203)+2*(p1mup3*p2mup2*F(675)+p2mup3*p3
     -   mup2*F(7204)+p2mup2*p3mup3*F(7205))+4*(p2mup2*p2mup3*F(1390)
     -   +p3mup2*p3mup3*F(7206))
       Fa9(2619) = 4*p2mup2*p2mup4*F(1281)+p1mup2*(-4*(p1mup4*F(634)-p
     -   2mup4*F(7124))+2*(p3mup4*F(7125)+p4mup4*F(7126)))+p3mup2*p4m
     -   up4*F(7207)+p3mup4*p4mup2*F(7208)+2*(p2mup4*p4mup2*F(1329)+p
     -   2mup2*p4mup4*F(1680)+p1mup4*(2*p2mup2*F(641)+p4mup2*F(1327)+
     -   p3mup2*F(1839))+p2mup4*p3mup2*F(1840)+p2mup2*p3mup4*F(2285)+
     -   p4mup2*p4mup4*F(7127)+p3mup2*p3mup4*F(7209))
       Fa9(2620) = 4*(p2mup4*p3mup2*F(1870)+p2mup2*p3mup4*F(1875)+p2mu
     -   p2*p2mup4*F(4087)+p3mup2*p3mup4*F(4088))+p1mup2*p4mup4*F(721
     -   2)+p1mup4*(-(p4mup2*F(7213))-2*(p2mup2*F(1289)-p1mup2*F(7210
     -   )-p3mup2*F(7217)))+2*(p1mup2*p3mup4*F(1871)+p1mup2*p2mup4*F(
     -   7211)+p3mup2*p4mup4*F(7215)+p4mup2*(2*p4mup4*F(1667)+p2mup4*
     -   F(7214)+p3mup4*F(7216))+p2mup2*p4mup4*F(7218))
       Fa9(2621) = p1mup2*p3mup4*F(7219)+p1mup4*(-2*(p2mup2*F(675)-p1m
     -   up2*F(7146)-p4mup2*F(7153))+p3mup2*F(7220))+2*(p1mup2*p4mup4
     -   *F(1392)+p2mup4*(2*(p4mup2*F(1391)+p2mup2*F(3878))+p1mup2*F(
     -   7147)+p3mup2*F(7221))+p3mup2*p4mup4*F(7222)+p3mup4*p4mup2*F(
     -   7223)+p2mup2*p3mup4*F(7224))+4*(p2mup2*p4mup4*F(3879)+p4mup2
     -   *p4mup4*F(3880)+p3mup2*p3mup4*F(7225))
       Fa9(2622) = p2mup2*p4mup4*F(7228)+p2mup4*p4mup2*F(7229)+2*(p4mu
     -   p2*p4mup4*F(1233)-p2mup2*p2mup4*F(1553)+p2mup2*p3mup4*F(2346
     -   )+p1mup4*(2*(p3mup2*F(284)-p1mup2*F(886)+p4mup2*F(1937))+p2m
     -   up2*F(7226))+p1mup2*(2*p3mup4*F(284)-p2mup4*F(1551)+p4mup4*F
     -   (7227))+p3mup2*p4mup4*F(7230)+p2mup4*p3mup2*F(7231))+4*(p3mu
     -   p4*p4mup2*F(1942)+p3mup2*p3mup4*F(7232))
       Fa9(2623) = p1mup4*p2mup2*F(7235)+p1mup2*(2*p1mup4*F(7233)+p2mu
     -   p4*F(7234)+4*(p3mup4*F(3685)+p4mup4*F(7237)))+2*(p1mup4*p3mu
     -   p2*F(7236)+p1mup4*p4mup2*F(7238)+p2mup4*p3mup2*F(7239)+p2mup
     -   4*p4mup2*F(7240)+p2mup2*p4mup4*F(7240)+2*(p3mup2*p4mup4*F(37
     -   72)+p4mup2*p4mup4*F(3773)+p2mup2*p2mup4*F(7241))+p3mup4*(p2m
     -   up2*F(7239)+2*(p4mup2*F(3772)+p3mup2*F(7242))))
       Fa9(2624) = 4*(p3mup2*p4mup4*F(1523)+p4mup2*p4mup4*F(7190))+p2m
     -   up2*p3mup4*F(7244)-p2mup4*p3mup2*F(7245)-2*(p2mup2*p2mup4*F(
     -   1518)-p2mup2*p4mup4*F(1777)-p1mup4*(2*(p4mup2*F(336)-p1mup2*
     -   F(960)+p3mup2*F(1520))+p2mup2*F(7186))-p2mup4*p4mup2*F(7189)
     -   -p1mup2*(2*p4mup4*F(336)-p2mup4*F(961)+p3mup4*F(7243))-p3mup
     -   2*p3mup4*F(7246)-p3mup4*p4mup2*F(7247))
       Fa9(2625) = 4*p1mup3*p3mup2*F(7084)+p2mup3*p3mup2*F(7245)+p2mup
     -   2*p3mup3*F(7248)+2*(p1mup3*p2mup2*F(963)+p2mup2*p2mup3*F(151
     -   8)+p1mup2*(2*p1mup3*F(960)+p2mup3*F(961)+p3mup3*F(7249))+p3m
     -   up2*p3mup3*F(7250))
       Fa9(2626) = 2*(p2mup3*p3mup2+p2mup2*p3mup3)*F(7253)+p1mup2*(2*p
     -   1mup3*F(908)+p2mup3*F(7251)+4*p3mup3*F(7254))+p1mup3*(p2mup2
     -   *F(7252)+2*p3mup2*F(7255))+4*(p2mup2*p2mup3*F(2731)+p3mup2*p
     -   3mup3*F(7256))
       Fa9(2627) = 2*(p1mup3*p2mup2*F(641)+p2mup2*p2mup3*F(1281))+p1mu
     -   p3*p3mup2*F(1839)+p2mup3*p3mup2*F(1840)+p2mup2*p3mup3*F(2285
     -   )+p1mup2*(-2*(p1mup3*F(634)-p2mup3*F(7124))+p3mup3*F(7125))+
     -   p3mup2*p3mup3*F(7209)
       Fa9(2628) = p3mup3*(p1mup2*F(1871)+2*(p2mup2*F(1875)+p3mup2*F(4
     -   088)))+p2mup3*(2*(p3mup2*F(1870)+p2mup2*F(4087))+p1mup2*F(72
     -   11))+p1mup3*(-(p2mup2*F(1289))+p1mup2*F(7210)+p3mup2*F(7217)
     -   )
       Return
       End
