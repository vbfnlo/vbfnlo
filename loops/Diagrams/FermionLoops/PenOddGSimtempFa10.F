c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenOddGSim10(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4
     -   ,p3mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5m
     -   up3,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mu
     -   p3,mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa10)
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
       Complex*16 Fa10(892:990)
       COMMON/FaPenOddGSim/Fa
       Fa10(892) = p3mup1*p4mup3*p4mup4*F(278)+p3mup1*p3mup3*p4mup4*F(
     -   283)+p3mup3*p4mup1*p4mup4*F(363)+p1mup4*(-4*p3mup1*p4mup3*F(
     -   261)+8*(p2mup3*p4mup1*F(343)-p4mup1*p4mup3*F(359)+p2mup3*p3m
     -   up1*F(2379)))-4*(p3mup1*p3mup4*p4mup3*F(274)-p3mup3*p3mup4*p
     -   4mup1*F(352)-p3mup4*p4mup1*p4mup3*F(358)-p2mup3*p3mup1*p4mup
     -   4*F(2380)+p2mup4*p3mup1*p4mup3*F(2381)-p2mup3*p3mup4*p4mup1*
     -   F(2390)-p2mup4*p3mup3*p4mup1*F(2391)-p1mup3*(p3mup1*p4mup4*F
     -   (259)-2*(p4mup1*p4mup4*F(359)-p1mup4*p3mup1*F(882)-p1mup4*p4
     -   mup1*F(938)-p2mup4*(p3mup1*F(254)+p4mup1*F(2392)))))-8*(p2mu
     -   p4*p3mup1*p3mup3*F(280)+p3mup1*p3mup3*p3mup4*F(281)+p2mup3*p
     -   4mup1*p4mup4*F(360)-p4mup1*p4mup3*p4mup4*F(365)-p2mup3*p3mup
     -   1*p3mup4*F(2378)-p2mup3*p2mup4*p3mup1*F(2382)-p2mup3*p2mup4*
     -   p4mup1*F(2393)+p2mup4*p4mup1*p4mup3*F(2394))
       Fa10(893) = p2mup1*p3mup4*p4mup3*F(274)+p2mup4*p3mup3*p4mup1*F(
     -   282)-p2mup3*p4mup1*p4mup4*F(362)+p2mup1*p4mup3*p4mup4*F(598)
     -   +p2mup1*p3mup3*p4mup4*F(599)-p1mup3*(p2mup4*p4mup1*F(258)+p2
     -   mup1*(p4mup4*F(259)+2*(p2mup4*F(254)+p3mup4*F(256)+p1mup4*F(
     -   882))))+p2mup3*p3mup4*p4mup1*F(2188)+p2mup4*p4mup1*p4mup3*F(
     -   2345)+p1mup4*(p2mup1*p4mup3*F(261)+p2mup3*p4mup1*F(346)+2*(p
     -   2mup1*p3mup3*F(279)-p2mup1*p2mup3*F(2379)))+2*(p2mup1*p2mup4
     -   *p3mup3*F(280)+p2mup1*p3mup3*p3mup4*F(281)-p2mup1*p2mup3*p3m
     -   up4*F(2378)-p2mup1*p2mup3*p2mup4*F(2382))+p2mup3*p2mup4*p4mu
     -   p1*F(2395)
       Fa10(894) = p2mup1*(-(p2mup4*p4mup2*F(640))+p4mup2*p4mup4*F(137
     -   0)+2*(p1mup2*p3mup4*F(636)+p3mup2*p3mup4*F(2397)))+p4mup1*(p
     -   2mup4*p4mup2*F(946)+2*(p4mup2*p4mup4*F(365)-p3mup2*p3mup4*F(
     -   2409)))
       Fa10(895) = -8*p1mup2*p3mup1*p3mup4*F(636)+p3mup1*(p4mup2*p4mup
     -   4*F(1685)+4*p3mup4*p4mup2*F(2396)-8*p3mup2*(p1mup4*F(2398)+p
     -   2mup4*F(2399))+p3mup2*p4mup4*F(2400))+4*p3mup2*p3mup4*p4mup1
     -   *F(2407)
       Fa10(896) = mup1mup2*(p4mup3*p4mup4*F(278)+p3mup3*p4mup4*F(283)
     -   -4*(p3mup4*p4mup3*F(274)-p1mup3*(p4mup4*F(259)+2*(p2mup4*F(2
     -   54)+p3mup4*F(256)+p1mup4*F(882)))+p1mup4*(p4mup3*F(261)+2*(p
     -   3mup3*F(279)-p2mup3*F(2379)))-p2mup3*p4mup4*F(2380)+p2mup4*p
     -   4mup3*F(2381))-8*(p2mup4*p3mup3*F(280)+p3mup3*p3mup4*F(281)-
     -   p2mup3*p3mup4*F(2378)-p2mup3*p2mup4*F(2382)))
       Fa10(897) = p1mup4*p2mup3*p3mup1*F(346)+p1mup4*p2mup1*p3mup3*F(
     -   347)+p2mup1*p3mup3*p3mup4*F(352)+p2mup1*p3mup4*p4mup3*F(358)
     -   -p2mup3*p3mup1*p4mup4*F(362)+2*(p1mup4*p2mup1*p2mup3*F(343)-
     -   p1mup4*p2mup1*p4mup3*F(359)+p2mup1*p4mup3*p4mup4*F(365))-p2m
     -   up1*p3mup3*p4mup4*F(632)+p2mup3*p3mup1*p3mup4*F(2188)+p2mup1
     -   *p2mup3*p3mup4*F(2390)+p1mup3*(-(p2mup4*p3mup1*F(258))+p2mup
     -   1*p3mup4*F(345)-2*(p2mup1*p4mup4*F(359)-p1mup4*p2mup1*F(938)
     -   -p2mup1*p2mup4*F(2392)))+p2mup4*(p3mup1*p3mup3*F(282)+p3mup1
     -   *p4mup3*F(2345)+p2mup1*p3mup3*F(2391)+2*p2mup1*p2mup3*F(2393
     -   )+p2mup3*p3mup1*F(2395))
       Fa10(898) = 4*p2mup1*p3mup4*(p1mup2*F(1262)+p3mup2*F(2413))+p3m
     -   up1*(8*(p1mup4*p3mup2*F(944)+p3mup2*p4mup4*F(949)-p1mup2*p3m
     -   up4*F(2408)-p4mup2*(p4mup4*F(365)+p3mup4*F(2411)))+p2mup4*p3
     -   mup2*F(2414))
       Fa10(899) = mup1mup2*(p3mup3*p4mup4*F(363)+4*(p1mup3*p3mup4*F(3
     -   45)+p3mup3*p3mup4*F(352)+p3mup4*p4mup3*F(358)+p1mup4*(p3mup3
     -   *F(347)+2*(p2mup3*F(343)-p4mup3*F(359)+p1mup3*F(938)))+p2mup
     -   3*p3mup4*F(2390)+p2mup4*p3mup3*F(2391))-8*(p4mup4*(p1mup3*F(
     -   359)+p2mup3*F(360))-p4mup3*p4mup4*F(365)-p1mup3*p2mup4*F(239
     -   2)-p2mup3*p2mup4*F(2393)+p2mup4*p4mup3*F(2394)))
       Fa10(900) = mup1mup2*(p2mup3*(p1mup4*F(346)-p4mup4*F(362)+p3mup
     -   4*F(2188))+p2mup4*(-(p1mup3*F(258))+p3mup3*F(282)+p4mup3*F(2
     -   345)+p2mup3*F(2395)))
       Fa10(901) = mup3mup5*(-(p2mup4*p4mup2*F(640))+p1mup2*(2*(p1mup4
     -   *F(633)+p2mup4*F(635)+p3mup4*F(636))+p4mup4*F(1367))-p1mup4*
     -   p4mup2*F(1368)+p4mup2*p4mup4*F(1370)-p3mup4*p4mup2*F(2396)+2
     -   *(p3mup2*p3mup4*F(2397)+p1mup4*p3mup2*F(2398)+p2mup4*p3mup2*
     -   F(2399))+p3mup2*p4mup4*F(2415))
       Fa10(902) = mup3mup5*(p1mup4*p4mup2*F(1707)-4*p1mup2*(p1mup4*F(
     -   677)+p2mup4*F(679)-p4mup4*F(1400)+p3mup4*F(2401))-8*(p4mup2*
     -   p4mup4*F(1402)-p3mup2*p3mup4*F(2402)-p2mup4*(p4mup2*F(689)+p
     -   3mup2*F(2404))+p3mup4*p4mup2*F(2405)-p3mup2*p4mup4*F(2406))+
     -   p1mup4*p3mup2*F(2412))
       Fa10(903) = p4mup1*(p1mup3*(p4mup4*F(300)+4*(p1mup4*F(889)+p2mu
     -   p4*F(2383)))+4*p1mup4*(p4mup3*F(301)-p2mup3*F(2384))+8*(-(p3
     -   mup3*p3mup4*F(307))+p3mup4*p4mup3*F(309)+p3mup3*p4mup4*F(317
     -   )+p4mup3*p4mup4*F(318)-p2mup4*p3mup3*F(2385)+p2mup3*p4mup4*F
     -   (2386)+p2mup3*p3mup4*F(2387)-p2mup4*p4mup3*F(2388)+p2mup3*p2
     -   mup4*F(2389)))
       Fa10(904) = mup1mup5*(p4mup3*p4mup4*F(278)+p3mup3*p4mup4*F(283)
     -   -4*(p3mup4*p4mup3*F(274)-p1mup3*(p4mup4*F(259)+2*(p2mup4*F(2
     -   54)+p3mup4*F(256)+p1mup4*F(882)))+p1mup4*(p4mup3*F(261)+2*(p
     -   3mup3*F(279)-p2mup3*F(2379)))-p2mup3*p4mup4*F(2380)+p2mup4*p
     -   4mup3*F(2381))-8*(p2mup4*p3mup3*F(280)+p3mup3*p3mup4*F(281)-
     -   p2mup3*p3mup4*F(2378)-p2mup3*p2mup4*F(2382)))
       Fa10(905) = -(mup1mup5*(4*(p1mup4*p4mup3*F(301)+p1mup3*(p3mup4*
     -   F(294)-p4mup4*F(609)+p1mup4*F(889)+p2mup4*F(2383))-p1mup4*p2
     -   mup3*F(2384))-p3mup3*(p1mup4*F(608)+8*(p3mup4*F(307)-p4mup4*
     -   F(317)+p2mup4*F(2385)))+8*(p3mup4*p4mup3*F(309)+p4mup3*p4mup
     -   4*F(318)+p2mup3*p4mup4*F(2386)+p2mup3*p3mup4*F(2387)-p2mup4*
     -   p4mup3*F(2388)+p2mup3*p2mup4*F(2389))))
       Fa10(906) = mup3mup5*(p1mup4*p3mup2*F(346)-p3mup2*p4mup4*F(362)
     -   +p1mup2*p3mup4*F(1262)+p2mup4*p3mup2*F(2318)+p3mup4*p4mup2*F
     -   (2369)+p3mup2*p3mup4*F(2413))
       Fa10(907) = mup1mup5*(p3mup3*p4mup4*F(363)+4*(p1mup3*p3mup4*F(3
     -   45)+p3mup3*p3mup4*F(352)+p3mup4*p4mup3*F(358)+p1mup4*(p3mup3
     -   *F(347)+2*(p2mup3*F(343)-p4mup3*F(359)+p1mup3*F(938)))+p2mup
     -   3*p3mup4*F(2390)+p2mup4*p3mup3*F(2391))-8*(p4mup4*(p1mup3*F(
     -   359)+p2mup3*F(360))-p4mup3*p4mup4*F(365)-p1mup3*p2mup4*F(239
     -   2)-p2mup3*p2mup4*F(2393)+p2mup4*p4mup3*F(2394)))
       Fa10(908) = p3mup1*(-(p1mup3*p4mup4*F(609))+p1mup3*p2mup4*F(238
     -   3)+p1mup4*(p4mup3*F(301)+p1mup3*F(889)-p2mup3*F(2384))+2*(p3
     -   mup4*p4mup3*F(309)+p4mup3*p4mup4*F(318)-p3mup3*(p3mup4*F(307
     -   )-p4mup4*F(317)+p2mup4*F(2385))+p2mup3*p4mup4*F(2386)+p2mup3
     -   *p3mup4*F(2387)-p2mup4*p4mup3*F(2388)+p2mup3*p2mup4*F(2389))
     -   )
       Fa10(909) = mup3mup5*(p2mup4*p4mup2*F(1221)+p1mup2*(-4*p2mup4*F
     -   (939)+8*(p4mup4*F(359)-p3mup4*F(2408)))-8*(p4mup2*p4mup4*F(3
     -   65)-p1mup4*(p4mup2*F(359)-p1mup2*F(938)+p3mup2*F(944))-p3mup
     -   2*p4mup4*F(949)-p3mup2*p3mup4*F(2409)+p3mup4*p4mup2*F(2411))
     -   +p2mup4*p3mup2*F(2414))
       Fa10(910) = p2mup1*(p1mup3*(p4mup4*F(300)+4*(p3mup4*F(294)+p1mu
     -   p4*F(889)+p2mup4*F(2383)))+4*(p1mup4*(p3mup3*F(299)+p4mup3*F
     -   (301))-p1mup4*p2mup3*F(2384))-8*(-(p3mup4*p4mup3*F(309))-p4m
     -   up3*p4mup4*F(318)+p3mup3*(p3mup4*F(307)-p4mup4*F(317)+p2mup4
     -   *F(2385))-p2mup3*p3mup4*F(2387)-p2mup3*p2mup4*F(2389)))
       Fa10(911) = 4*p2mup1*(-(p1mup2*p3mup4*F(2401))+2*(p2mup4*p4mup2
     -   *F(689)-p4mup2*p4mup4*F(1402)+p3mup2*p3mup4*F(2402)))
       Fa10(912) = p3mup1*(p1mup2*p3mup4*F(2401)+2*p4mup2*(p4mup4*F(14
     -   02)+p3mup4*F(2405))+p3mup2*(p1mup4*F(2403)-2*(p2mup4*F(2404)
     -   +p4mup4*F(2406))))
       Fa10(913) = mup1mup2*(p1mup3*(p4mup4*F(300)+4*(p3mup4*F(294)+p1
     -   mup4*F(889)+p2mup4*F(2383)))+4*(p1mup4*(p3mup3*F(299)+p4mup3
     -   *F(301))-p1mup4*p2mup3*F(2384))-8*(-(p3mup4*p4mup3*F(309))-p
     -   4mup3*p4mup4*F(318)+p3mup3*(p3mup4*F(307)-p4mup4*F(317)+p2mu
     -   p4*F(2385))-p2mup3*p4mup4*F(2386)-p2mup3*p3mup4*F(2387)+p2mu
     -   p4*p4mup3*F(2388)-p2mup3*p2mup4*F(2389)))
       Fa10(914) = -(p3mup1*p3mup4*p4mup3*F(274))+p3mup3*p3mup4*p4mup1
     -   *F(352)+p3mup4*p4mup1*p4mup3*F(358)-p3mup1*p4mup3*p4mup4*F(5
     -   98)-p3mup1*p3mup3*p4mup4*F(599)-p3mup3*p4mup1*p4mup4*F(632)+
     -   p1mup4*(-(p3mup1*p4mup3*F(261))+2*(p2mup3*p4mup1*F(343)-p4mu
     -   p1*p4mup3*F(359)+p2mup3*p3mup1*F(2379)))+p2mup3*p3mup1*p4mup
     -   4*F(2380)-p2mup4*p3mup1*p4mup3*F(2381)+p2mup3*p3mup4*p4mup1*
     -   F(2390)+p2mup4*p3mup3*p4mup1*F(2391)+p1mup3*(p3mup1*p4mup4*F
     -   (259)-2*(p4mup1*p4mup4*F(359)-p1mup4*p3mup1*F(882)-p1mup4*p4
     -   mup1*F(938)-p2mup4*(p3mup1*F(254)+p4mup1*F(2392))))-2*(p2mup
     -   4*p3mup1*p3mup3*F(280)+p3mup1*p3mup3*p3mup4*F(281)+p2mup3*p4
     -   mup1*p4mup4*F(360)-p4mup1*p4mup3*p4mup4*F(365)-p2mup3*p3mup1
     -   *p3mup4*F(2378)-p2mup3*p2mup4*p3mup1*F(2382)-p2mup3*p2mup4*p
     -   4mup1*F(2393)+p2mup4*p4mup1*p4mup3*F(2394))
       Fa10(915) = p2mup1*p4mup3*p4mup4*F(278)+p2mup1*p3mup3*p4mup4*F(
     -   283)-8*(p2mup1*p2mup4*p3mup3*F(280)+p2mup1*p3mup3*p3mup4*F(2
     -   81)-p2mup1*p2mup3*p3mup4*F(2378)-p2mup1*p2mup3*p2mup4*F(2382
     -   ))-4*(p2mup1*p3mup4*p4mup3*F(274)+p2mup4*p3mup3*p4mup1*F(282
     -   )-p2mup3*p4mup1*p4mup4*F(362)-p1mup3*(p2mup4*p4mup1*F(258)+p
     -   2mup1*(p4mup4*F(259)+2*(p2mup4*F(254)+p3mup4*F(256)+p1mup4*F
     -   (882))))+p2mup3*p3mup4*p4mup1*F(2188)+p2mup4*p4mup1*p4mup3*F
     -   (2345)+p1mup4*(p2mup1*p4mup3*F(261)+p2mup3*p4mup1*F(346)+2*(
     -   p2mup1*p3mup3*F(279)-p2mup1*p2mup3*F(2379)))+p2mup3*p2mup4*p
     -   4mup1*F(2395))
       Fa10(916) = p2mup1*(4*p2mup4*p4mup2*F(640)+p4mup2*p4mup4*F(1685
     -   )-8*(p1mup2*p3mup4*F(636)+p3mup2*p3mup4*F(2397)))+p4mup1*(p2
     -   mup4*p4mup2*F(1221)-8*(p4mup2*p4mup4*F(365)-p3mup2*p3mup4*F(
     -   2409)))
       Fa10(917) = 2*p1mup2*p3mup1*p3mup4*F(636)+p3mup2*p3mup4*p4mup1*
     -   F(2413)+p3mup1*(p4mup2*p4mup4*F(1370)-p3mup4*p4mup2*F(2396)+
     -   2*(p1mup4*p3mup2*F(2398)+p2mup4*p3mup2*F(2399))+p3mup2*p4mup
     -   4*F(2415))
       Fa10(918) = mup1mup2*(p3mup4*p4mup3*F(274)+p4mup3*p4mup4*F(598)
     -   +p3mup3*p4mup4*F(599)-p1mup3*(p4mup4*F(259)+2*(p2mup4*F(254)
     -   +p3mup4*F(256)+p1mup4*F(882)))+p1mup4*(p4mup3*F(261)+2*(p3mu
     -   p3*F(279)-p2mup3*F(2379)))-p2mup3*p4mup4*F(2380)+p2mup4*p4mu
     -   p3*F(2381)+2*(p2mup4*p3mup3*F(280)+p3mup3*p3mup4*F(281)-p2mu
     -   p3*p3mup4*F(2378)-p2mup3*p2mup4*F(2382)))
       Fa10(919) = p2mup1*p3mup3*p4mup4*F(363)+8*(p1mup4*p2mup1*p2mup3
     -   *F(343)-p1mup4*p2mup1*p4mup3*F(359)+p2mup1*p4mup3*p4mup4*F(3
     -   65))+p1mup3*(-4*(p2mup4*p3mup1*F(258)-p2mup1*p3mup4*F(345))-
     -   8*(p2mup1*p4mup4*F(359)-p1mup4*p2mup1*F(938)-p2mup1*p2mup4*F
     -   (2392)))+4*(p1mup4*p2mup3*p3mup1*F(346)+p1mup4*p2mup1*p3mup3
     -   *F(347)+p2mup1*p3mup3*p3mup4*F(352)+p2mup1*p3mup4*p4mup3*F(3
     -   58)-p2mup3*p3mup1*p4mup4*F(362)+p2mup3*p3mup1*p3mup4*F(2188)
     -   +p2mup1*p2mup3*p3mup4*F(2390)+p2mup4*(p3mup1*p3mup3*F(282)+p
     -   3mup1*p4mup3*F(2345)+p2mup1*p3mup3*F(2391)+2*p2mup1*p2mup3*F
     -   (2393)+p2mup3*p3mup1*F(2395)))
       Fa10(920) = p2mup1*p3mup4*(-(p1mup2*F(1262))+p3mup2*F(2407))+p3
     -   mup1*(p2mup4*p3mup2*F(2410)-2*(p1mup4*p3mup2*F(944)+p3mup2*p
     -   4mup4*F(949)-p1mup2*p3mup4*F(2408)-p4mup2*(p4mup4*F(365)+p3m
     -   up4*F(2411))))
       Fa10(921) = mup1mup2*(p1mup3*p3mup4*F(345)+p3mup4*p4mup3*F(358)
     -   +p1mup4*(p3mup3*F(347)+2*(p2mup3*F(343)-p4mup3*F(359)+p1mup3
     -   *F(938)))+p2mup3*p3mup4*F(2390)+p3mup3*(p3mup4*F(352)-p4mup4
     -   *F(632)+p2mup4*F(2391))-2*(p4mup4*(p1mup3*F(359)+p2mup3*F(36
     -   0))-p4mup3*p4mup4*F(365)-p1mup3*p2mup4*F(2392)-p2mup3*p2mup4
     -   *F(2393)+p2mup4*p4mup3*F(2394)))
       Fa10(922) = mup4mup5*(-(p1mup3*p4mup2*F(258))+p3mup3*p4mup2*F(2
     -   82)+p1mup2*p4mup3*F(1262)-p2mup3*p4mup2*F(2338)+p3mup2*p4mup
     -   3*F(2368)+p4mup2*p4mup3*F(2416))
       Fa10(923) = mup4mup5*(p1mup3*p3mup2*F(1660)+8*(-(p3mup2*p3mup3*
     -   F(1306))+p3mup3*p4mup2*F(1307)+p3mup2*p4mup3*F(1307)-p4mup2*
     -   p4mup3*F(1308)+p2mup3*(p3mup2*F(1299)+p4mup2*F(2417)))-4*p1m
     -   up2*(p1mup3*F(1289)+p2mup3*F(1290)-p3mup3*F(1300)-p4mup3*F(2
     -   418))+p1mup3*p4mup2*F(2419))
       Fa10(924) = p4mup1*(p1mup4*p4mup2*F(1401)+p1mup2*(p1mup4*F(677)
     -   +p2mup4*F(679)-p4mup4*F(1400)+p3mup4*F(2401))+p1mup4*p3mup2*
     -   F(2403)-2*(p2mup4*p3mup2*F(2404)-p3mup4*p4mup2*F(2405)+p3mup
     -   2*p4mup4*F(2406)))
       Fa10(925) = mup1mup5*(-(p2mup4*p4mup2*F(640))+p1mup2*(2*(p1mup4
     -   *F(633)+p2mup4*F(635)+p3mup4*F(636))+p4mup4*F(1367))-p1mup4*
     -   p4mup2*F(1368)+p4mup2*p4mup4*F(1370)-p3mup4*p4mup2*F(2396)+2
     -   *(p3mup2*p3mup4*F(2397)+p1mup4*p3mup2*F(2398)+p2mup4*p3mup2*
     -   F(2399))+p3mup2*p4mup4*F(2415))
       Fa10(926) = mup1mup5*(p1mup4*p4mup2*F(1707)-4*p1mup2*(p1mup4*F(
     -   677)+p2mup4*F(679)-p4mup4*F(1400)+p3mup4*F(2401))-8*(p4mup2*
     -   p4mup4*F(1402)-p3mup2*p3mup4*F(2402)-p2mup4*(p4mup2*F(689)+p
     -   3mup2*F(2404))+p3mup4*p4mup2*F(2405)-p3mup2*p4mup4*F(2406))+
     -   p1mup4*p3mup2*F(2412))
       Fa10(927) = mup4mup5*(4*(p1mup3*p3mup2*F(1263)+p2mup3*p3mup2*F(
     -   1265)-p1mup2*(p3mup3*F(1261)+2*(p1mup3*F(633)+p2mup3*F(635)+
     -   p4mup3*F(1276))))+p3mup2*p3mup3*F(1646)+p3mup2*p4mup3*F(2420
     -   )+p3mup3*p4mup2*F(2422)-8*(p1mup3*p4mup2*F(1276)+p2mup3*p4mu
     -   p2*F(2421)+p4mup2*p4mup3*F(2423)))
       Fa10(928) = mup1mup5*(p1mup4*p3mup2*F(346)-p3mup2*p4mup4*F(362)
     -   +p1mup2*p3mup4*F(1262)+p2mup4*p3mup2*F(2318)+p3mup4*p4mup2*F
     -   (2369)+p3mup2*p3mup4*F(2413))
       Fa10(929) = mup1mup5*p4mup2*(p2mup3*F(2338)-p4mup3*F(2416))
       Fa10(930) = p3mup1*(p1mup3*p3mup2*F(1301)-2*(-(p3mup2*p3mup3*F(
     -   1306))+p3mup3*p4mup2*F(1307)+p3mup2*p4mup3*F(1307)-p4mup2*p4
     -   mup3*F(1308)+p2mup3*(p3mup2*F(1299)+p4mup2*F(2417)))+p1mup2*
     -   (p1mup3*F(1289)+p2mup3*F(1290)-p3mup3*F(1300)-p4mup3*F(2418)
     -   )+p1mup3*p4mup2*F(2424))
       Fa10(931) = mup1mup5*(p1mup2*p3mup3*F(1261)+p3mup2*p3mup3*F(127
     -   4)+2*(p1mup2*p1mup3*F(633)+p2mup3*p4mup2*F(2421)+p4mup2*p4mu
     -   p3*F(2423)))
       Fa10(932) = 4*mup1mup5*(p1mup2*p1mup3*F(1289)-p1mup2*p3mup3*F(1
     -   300)+2*(p3mup2*p3mup3*F(1306)+p4mup2*p4mup3*F(1308)-p2mup3*p
     -   4mup2*F(2417)))
       Fa10(933) = p4mup1*(p1mup3*p3mup2*F(1660)+8*(-(p3mup2*p3mup3*F(
     -   1306))+p3mup3*p4mup2*F(1307)+p3mup2*p4mup3*F(1307)-p4mup2*p4
     -   mup3*F(1308)+p2mup3*(p3mup2*F(1299)+p4mup2*F(2417)))-4*p1mup
     -   2*(p1mup3*F(1289)+p2mup3*F(1290)-p3mup3*F(1300)-p4mup3*F(241
     -   8))+p1mup3*p4mup2*F(2419))
       Fa10(934) = mup4mup5*(p1mup2*(-4*p2mup3*F(1499)+8*(p3mup3*F(279
     -   )+p4mup3*F(1503)))+p2mup3*p3mup2*F(1778)+p2mup3*p4mup2*F(242
     -   5)+8*(p1mup3*(p3mup2*F(279)-p1mup2*F(882)+p4mup2*F(1503))-p3
     -   mup2*p3mup3*F(1507)+p3mup3*p4mup2*F(1508)+p3mup2*p4mup3*F(15
     -   08)-p4mup2*p4mup3*F(2426)))
       Fa10(935) = p2mup1*(p1mup3*p3mup2*F(1660)+8*(-(p3mup2*p3mup3*F(
     -   1306))+p3mup3*p4mup2*F(1307)+p3mup2*p4mup3*F(1307)-p4mup2*p4
     -   mup3*F(1308)+p2mup3*(p3mup2*F(1299)+p4mup2*F(2417)))-4*p1mup
     -   2*(p1mup3*F(1289)+p2mup3*F(1290)-p3mup3*F(1300)-p4mup3*F(241
     -   8))+p1mup3*p4mup2*F(2419))
       Fa10(936) = 4*mup1mup5*(-2*(p1mup2*p3mup3*F(279)-p1mup2*p1mup3*
     -   F(882)-p3mup2*p3mup3*F(1507)-p4mup2*p4mup3*F(2426))+p2mup3*p
     -   4mup2*F(2427))
       Fa10(937) = mup1mup5*(p2mup4*p4mup2*F(1221)+p1mup2*(-4*p2mup4*F
     -   (939)+8*(p4mup4*F(359)-p3mup4*F(2408)))-8*(p4mup2*p4mup4*F(3
     -   65)-p1mup4*(p4mup2*F(359)-p1mup2*F(938)+p3mup2*F(944))-p3mup
     -   2*p4mup4*F(949)-p3mup2*p3mup4*F(2409)+p3mup4*p4mup2*F(2411))
     -   +p2mup4*p3mup2*F(2414))
       Fa10(938) = p2mup1*(p1mup2*(p1mup4*F(677)+p2mup4*F(679)-p4mup4*
     -   F(1400))+p1mup4*p4mup2*F(1401)+p1mup4*p3mup2*F(2403)-2*(p2mu
     -   p4*p3mup2*F(2404)-p3mup4*p4mup2*F(2405)+p3mup2*p4mup4*F(2406
     -   )))
       Fa10(939) = p3mup1*(-4*p1mup2*(p1mup4*F(677)+p2mup4*F(679)-p4mu
     -   p4*F(1400))+p1mup4*p4mup2*F(1707)+8*(p2mup4*p4mup2*F(689)+p3
     -   mup2*p3mup4*F(2402)))
       Fa10(940) = mup1mup3*(p1mup4*p4mup2*F(1401)+p1mup2*(p1mup4*F(67
     -   7)+p2mup4*F(679)-p4mup4*F(1400)+p3mup4*F(2401))+p1mup4*p3mup
     -   2*F(2403)+2*(p4mup2*p4mup4*F(1402)-p3mup2*p3mup4*F(2402)-p2m
     -   up4*(p4mup2*F(689)+p3mup2*F(2404))+p3mup4*p4mup2*F(2405)-p3m
     -   up2*p4mup4*F(2406)))
       Fa10(941) = -4*(p1mup3*(p3mup1*p4mup2*F(258)+p4mup1*(-(p3mup2*F
     -   (1263))+2*(p1mup2*F(633)+p4mup2*F(1276))))-p3mup1*(p3mup3*p4
     -   mup2*F(282)+p1mup2*p4mup3*F(1262)-p2mup3*p4mup2*F(2338)+p3mu
     -   p2*p4mup3*F(2368)+p4mup2*p4mup3*F(2416)))+p4mup1*(4*(p2mup3*
     -   p3mup2*F(1265)-p1mup2*(p3mup3*F(1261)+2*(p2mup3*F(635)+p4mup
     -   3*F(1276))))+p3mup2*p3mup3*F(1646)+p3mup2*p4mup3*F(2420)+p3m
     -   up3*p4mup2*F(2422)-8*(p2mup3*p4mup2*F(2421)+p4mup2*p4mup3*F(
     -   2423)))
       Fa10(942) = p1mup2*(-4*(p2mup1*p4mup3*F(1262)+p2mup3*p4mup1*F(1
     -   499))+8*(p3mup3*p4mup1*F(279)+p4mup1*p4mup3*F(1503)))+p2mup3
     -   *p3mup2*p4mup1*F(1778)-4*(p2mup1*p3mup3*p4mup2*F(282)-p1mup3
     -   *(p2mup1*p4mup2*F(258)+2*p4mup1*(p3mup2*F(279)-p1mup2*F(882)
     -   +p4mup2*F(1503)))-p2mup1*p2mup3*p4mup2*F(2338)+p2mup1*p3mup2
     -   *p4mup3*F(2368)+p2mup1*p4mup2*p4mup3*F(2416))+p2mup3*p4mup1*
     -   p4mup2*F(2425)+8*(-(p3mup2*p3mup3*p4mup1*F(1507))+p3mup3*p4m
     -   up1*p4mup2*F(1508)+p3mup2*p4mup1*p4mup3*F(1508)-p4mup1*p4mup
     -   2*p4mup3*F(2426))
       Fa10(943) = -(p2mup1*p3mup4*p4mup2*F(2396))+p1mup4*(-2*p4mup1*(
     -   p4mup2*F(359)-p1mup2*F(938)+p3mup2*F(944))+p2mup1*(-(p4mup2*
     -   F(1368))+2*(p1mup2*F(633)+p3mup2*F(2398))))+p1mup2*(p2mup4*p
     -   4mup1*F(939)+p2mup1*p4mup4*F(1367)-2*(p4mup1*p4mup4*F(359)-p
     -   2mup1*p2mup4*F(635)-p3mup4*p4mup1*F(2408)))+p2mup4*p3mup2*p4
     -   mup1*F(2410)-2*(p3mup2*p4mup1*p4mup4*F(949)-p2mup1*p2mup4*p3
     -   mup2*F(2399)-p3mup4*p4mup1*p4mup2*F(2411))+p2mup1*p3mup2*p4m
     -   up4*F(2415)
       Fa10(944) = 2*p1mup2*p2mup4*p3mup1*F(635)-p2mup4*p3mup1*p4mup2*
     -   F(640)+p1mup2*p3mup4*p4mup1*F(1262)+p1mup2*p3mup1*p4mup4*F(1
     -   367)+p1mup4*(p3mup2*p4mup1*F(346)+2*p1mup2*p3mup1*F(633)-p3m
     -   up1*p4mup2*F(1368))+p3mup4*p4mup1*p4mup2*F(2369)+p3mup2*(-(p
     -   4mup1*p4mup4*F(362))+p2mup4*p4mup1*F(2318)+2*p3mup1*p3mup4*F
     -   (2397))
       Fa10(945) = mup1mup3*(p4mup2*p4mup4*F(1685)+4*(p2mup4*p4mup2*F(
     -   640)-p1mup2*(2*(p1mup4*F(633)+p2mup4*F(635)+p3mup4*F(636))+p
     -   4mup4*F(1367))+p1mup4*p4mup2*F(1368)+p3mup4*p4mup2*F(2396))-
     -   8*(p3mup2*p3mup4*F(2397)+p1mup4*p3mup2*F(2398)+p2mup4*p3mup2
     -   *F(2399))+p3mup2*p4mup4*F(2400))
       Fa10(946) = p2mup4*p3mup1*p4mup2*F(1221)-4*(p2mup1*p3mup2*p4mup
     -   4*F(362)-p1mup4*(p2mup1*p3mup2*F(346)+2*(p3mup1*p4mup2*F(359
     -   )-p1mup2*p3mup1*F(938)))+p1mup2*p2mup4*p3mup1*F(939)-p2mup1*
     -   p2mup4*p3mup2*F(2318)-p2mup1*p3mup4*p4mup2*F(2369))+8*(p1mup
     -   2*p3mup1*p4mup4*F(359)+p3mup1*p3mup2*p3mup4*F(2409))
       Fa10(947) = mup1mup3*(p1mup4*p3mup2*F(346)-p3mup2*p4mup4*F(362)
     -   +p1mup2*p3mup4*F(1262)+p2mup4*p3mup2*F(2318)+p3mup4*p4mup2*F
     -   (2369)-p3mup2*p3mup4*F(2407))
       Fa10(948) = mup1mup3*(p2mup4*p4mup2*F(946)+p1mup2*(p2mup4*F(939
     -   )-2*(p4mup4*F(359)-p3mup4*F(2408)))+p2mup4*p3mup2*F(2410)+2*
     -   (p4mup2*p4mup4*F(365)-p1mup4*(p4mup2*F(359)-p1mup2*F(938)+p3
     -   mup2*F(944))-p3mup2*p4mup4*F(949)-p3mup2*p3mup4*F(2409)+p3mu
     -   p4*p4mup2*F(2411)))
       Fa10(949) = -(p2mup1*p2mup3*p3mup2*F(1265))+p2mup1*p3mup2*p3mup
     -   3*F(1274)+p1mup3*(p2mup1*(-(p3mup2*F(1263))+2*(p1mup2*F(633)
     -   +p4mup2*F(1276)))-2*p3mup1*(p3mup2*F(279)-p1mup2*F(882)+p4mu
     -   p2*F(1503)))+p1mup2*(p2mup1*p3mup3*F(1261)+p2mup3*p3mup1*F(1
     -   499)-2*(p3mup1*p3mup3*F(279)-p2mup1*p2mup3*F(635)-p2mup1*p4m
     -   up3*F(1276)+p3mup1*p4mup3*F(1503)))+p2mup3*p3mup1*p3mup2*F(1
     -   505)-2*(-(p3mup1*p3mup2*p3mup3*F(1507))+p3mup1*p3mup3*p4mup2
     -   *F(1508)+p3mup1*p3mup2*p4mup3*F(1508)-p2mup1*p2mup3*p4mup2*F
     -   (2421)-p2mup1*p4mup2*p4mup3*F(2423)-p3mup1*p4mup2*p4mup3*F(2
     -   426))+p2mup3*p3mup1*p4mup2*F(2427)+p2mup1*p3mup2*p4mup3*F(24
     -   28)+p2mup1*p3mup3*p4mup2*F(2429)
       Fa10(950) = -(mup4mup5*(-(p1mup3*p4mup2*F(258))+p3mup3*p4mup2*F
     -   (282)+p1mup2*p4mup3*F(1262)-p2mup3*p4mup2*F(2338)+p3mup2*p4m
     -   up3*F(2368)+p4mup2*p4mup3*F(2416)))
       Fa10(951) = mup4mup5*(p1mup3*p3mup2*F(1301)-2*(-(p3mup2*p3mup3*
     -   F(1306))+p3mup3*p4mup2*F(1307)+p3mup2*p4mup3*F(1307)-p4mup2*
     -   p4mup3*F(1308)+p2mup3*(p3mup2*F(1299)+p4mup2*F(2417)))+p1mup
     -   2*(p1mup3*F(1289)+p2mup3*F(1290)-p3mup3*F(1300)-p4mup3*F(241
     -   8))+p1mup3*p4mup2*F(2424))
       Fa10(952) = p4mup1*(p1mup4*p4mup2*F(1707)-4*p1mup2*(p1mup4*F(67
     -   7)+p2mup4*F(679)-p4mup4*F(1400)+p3mup4*F(2401))+8*(p2mup4*p3
     -   mup2*F(2404)-p3mup4*p4mup2*F(2405)+p3mup2*p4mup4*F(2406))+p1
     -   mup4*p3mup2*F(2412))
       Fa10(953) = mup1mup5*(p4mup2*p4mup4*F(1685)+4*(p2mup4*p4mup2*F(
     -   640)-p1mup2*(2*(p1mup4*F(633)+p2mup4*F(635)+p3mup4*F(636))+p
     -   4mup4*F(1367))+p1mup4*p4mup2*F(1368)+p3mup4*p4mup2*F(2396))-
     -   8*(p3mup2*p3mup4*F(2397)+p1mup4*p3mup2*F(2398)+p2mup4*p3mup2
     -   *F(2399))+p3mup2*p4mup4*F(2400))
       Fa10(954) = mup1mup5*(p1mup4*p4mup2*F(1401)+p1mup2*(p1mup4*F(67
     -   7)+p2mup4*F(679)-p4mup4*F(1400)+p3mup4*F(2401))+p1mup4*p3mup
     -   2*F(2403)+2*(p4mup2*p4mup4*F(1402)-p3mup2*p3mup4*F(2402)-p2m
     -   up4*(p4mup2*F(689)+p3mup2*F(2404))+p3mup4*p4mup2*F(2405)-p3m
     -   up2*p4mup4*F(2406)))
       Fa10(955) = mup4mup5*(-(p1mup3*p3mup2*F(1263))-p2mup3*p3mup2*F(
     -   1265)+p3mup2*p3mup3*F(1274)+p1mup2*(p3mup3*F(1261)+2*(p1mup3
     -   *F(633)+p2mup3*F(635)+p4mup3*F(1276)))+2*(p1mup3*p4mup2*F(12
     -   76)+p2mup3*p4mup2*F(2421)+p4mup2*p4mup3*F(2423))+p3mup2*p4mu
     -   p3*F(2428)+p3mup3*p4mup2*F(2429))
       Fa10(956) = mup1mup5*(p1mup4*p3mup2*F(346)-p3mup2*p4mup4*F(362)
     -   +p1mup2*p3mup4*F(1262)+p2mup4*p3mup2*F(2318)+p3mup4*p4mup2*F
     -   (2369)-p3mup2*p3mup4*F(2407))
       Fa10(957) = p3mup1*(p1mup3*p3mup2*F(1660)+8*(-(p3mup2*p3mup3*F(
     -   1306))+p3mup3*p4mup2*F(1307)+p3mup2*p4mup3*F(1307)-p4mup2*p4
     -   mup3*F(1308)+p2mup3*(p3mup2*F(1299)+p4mup2*F(2417)))-4*p1mup
     -   2*(p1mup3*F(1289)+p2mup3*F(1290)-p3mup3*F(1300)-p4mup3*F(241
     -   8))+p1mup3*p4mup2*F(2419))
       Fa10(958) = -(mup1mup5*(4*p1mup2*p3mup3*F(1261)-p3mup2*p3mup3*F
     -   (1646)+8*(p1mup2*p1mup3*F(633)+p2mup3*p4mup2*F(2421)+p4mup2*
     -   p4mup3*F(2423))))
       Fa10(959) = -4*mup1mup5*(p1mup2*p1mup3*F(1289)-p1mup2*p3mup3*F(
     -   1300)+2*(p3mup2*p3mup3*F(1306)+p4mup2*p4mup3*F(1308)-p2mup3*
     -   p4mup2*F(2417)))
       Fa10(960) = p4mup1*(p1mup3*p3mup2*F(1301)-2*(-(p3mup2*p3mup3*F(
     -   1306))+p3mup3*p4mup2*F(1307)+p3mup2*p4mup3*F(1307)-p4mup2*p4
     -   mup3*F(1308)+p2mup3*(p3mup2*F(1299)+p4mup2*F(2417)))+p1mup2*
     -   (p1mup3*F(1289)+p2mup3*F(1290)-p3mup3*F(1300)-p4mup3*F(2418)
     -   )+p1mup3*p4mup2*F(2424))
       Fa10(961) = mup4mup5*(p1mup2*(p2mup3*F(1499)-2*(p3mup3*F(279)+p
     -   4mup3*F(1503)))+p2mup3*p3mup2*F(1505)-2*(p1mup3*(p3mup2*F(27
     -   9)-p1mup2*F(882)+p4mup2*F(1503))-p3mup2*p3mup3*F(1507)+p3mup
     -   3*p4mup2*F(1508)+p3mup2*p4mup3*F(1508)-p4mup2*p4mup3*F(2426)
     -   )+p2mup3*p4mup2*F(2427))
       Fa10(962) = p2mup1*(p1mup3*p3mup2*F(1301)-2*(-(p3mup2*p3mup3*F(
     -   1306))+p3mup3*p4mup2*F(1307)+p3mup2*p4mup3*F(1307)-p4mup2*p4
     -   mup3*F(1308)+p2mup3*(p3mup2*F(1299)+p4mup2*F(2417)))+p1mup2*
     -   (p1mup3*F(1289)+p2mup3*F(1290)-p3mup3*F(1300)-p4mup3*F(2418)
     -   )+p1mup3*p4mup2*F(2424))
       Fa10(963) = mup1mup5*(p2mup3*p4mup2*F(2425)+8*(p1mup2*p3mup3*F(
     -   279)-p1mup2*p1mup3*F(882)-p3mup2*p3mup3*F(1507)-p4mup2*p4mup
     -   3*F(2426)))
       Fa10(964) = mup1mup5*(p2mup4*p4mup2*F(946)+p1mup2*(p2mup4*F(939
     -   )-2*(p4mup4*F(359)-p3mup4*F(2408)))+p2mup4*p3mup2*F(2410)+2*
     -   (p4mup2*p4mup4*F(365)-p1mup4*(p4mup2*F(359)-p1mup2*F(938)+p3
     -   mup2*F(944))-p3mup2*p4mup4*F(949)-p3mup2*p3mup4*F(2409)+p3mu
     -   p4*p4mup2*F(2411)))
       Fa10(965) = p2mup1*(-4*p1mup2*(p1mup4*F(677)+p2mup4*F(679)-p4mu
     -   p4*F(1400))+p1mup4*p4mup2*F(1707)+8*(p2mup4*p3mup2*F(2404)-p
     -   3mup4*p4mup2*F(2405)+p3mup2*p4mup4*F(2406))+p1mup4*p3mup2*F(
     -   2412))
       Fa10(966) = p3mup1*(p1mup2*(p1mup4*F(677)+p2mup4*F(679)-p4mup4*
     -   F(1400))+p1mup4*p4mup2*F(1401)-2*(p2mup4*p4mup2*F(689)+p3mup
     -   2*p3mup4*F(2402)))
       Fa10(967) = mup1mup3*(p1mup4*p4mup2*F(1707)-4*p1mup2*(p1mup4*F(
     -   677)+p2mup4*F(679)-p4mup4*F(1400)+p3mup4*F(2401))-8*(p4mup2*
     -   p4mup4*F(1402)-p3mup2*p3mup4*F(2402)-p2mup4*(p4mup2*F(689)+p
     -   3mup2*F(2404))+p3mup4*p4mup2*F(2405)-p3mup2*p4mup4*F(2406))+
     -   p1mup4*p3mup2*F(2412))
       Fa10(968) = p1mup3*(p3mup1*p4mup2*F(258)+p4mup1*(-(p3mup2*F(126
     -   3))+2*(p1mup2*F(633)+p4mup2*F(1276))))-p3mup1*(p3mup3*p4mup2
     -   *F(282)+p1mup2*p4mup3*F(1262)-p2mup3*p4mup2*F(2338)+p3mup2*p
     -   4mup3*F(2368)+p4mup2*p4mup3*F(2416))+p4mup1*(p1mup2*p3mup3*F
     -   (1261)-p2mup3*p3mup2*F(1265)+p3mup2*p3mup3*F(1274)+2*(p1mup2
     -   *p2mup3*F(635)+p1mup2*p4mup3*F(1276)+p2mup3*p4mup2*F(2421)+p
     -   4mup2*p4mup3*F(2423))+p3mup2*p4mup3*F(2428)+p3mup3*p4mup2*F(
     -   2429))
       Fa10(969) = p2mup1*p3mup3*p4mup2*F(282)-p1mup3*(p2mup1*p4mup2*F
     -   (258)+2*p4mup1*(p3mup2*F(279)-p1mup2*F(882)+p4mup2*F(1503)))
     -   +p1mup2*(p2mup1*p4mup3*F(1262)+p2mup3*p4mup1*F(1499)-2*(p3mu
     -   p3*p4mup1*F(279)+p4mup1*p4mup3*F(1503)))+p2mup3*p3mup2*p4mup
     -   1*F(1505)-p2mup1*p2mup3*p4mup2*F(2338)+p2mup1*p3mup2*p4mup3*
     -   F(2368)+p2mup1*p4mup2*p4mup3*F(2416)-2*(-(p3mup2*p3mup3*p4mu
     -   p1*F(1507))+p3mup3*p4mup1*p4mup2*F(1508)+p3mup2*p4mup1*p4mup
     -   3*F(1508)-p4mup1*p4mup2*p4mup3*F(2426))+p2mup3*p4mup1*p4mup2
     -   *F(2427)
       Fa10(970) = p2mup1*p3mup2*p4mup4*F(2400)+4*(p2mup1*p3mup4*p4mup
     -   2*F(2396)+p1mup4*(2*p4mup1*(p4mup2*F(359)-p1mup2*F(938)+p3mu
     -   p2*F(944))+p2mup1*(p4mup2*F(1368)-2*(p1mup2*F(633)+p3mup2*F(
     -   2398))))-p1mup2*(p2mup4*p4mup1*F(939)+p2mup1*p4mup4*F(1367)-
     -   2*(p4mup1*p4mup4*F(359)-p2mup1*p2mup4*F(635)-p3mup4*p4mup1*F
     -   (2408))))+8*(p3mup2*p4mup1*p4mup4*F(949)-p2mup1*p2mup4*p3mup
     -   2*F(2399)-p3mup4*p4mup1*p4mup2*F(2411))+p2mup4*p3mup2*p4mup1
     -   *F(2414)
       Fa10(971) = mup1mup3*(-(p2mup4*p4mup2*F(640))+p1mup2*(2*(p1mup4
     -   *F(633)+p2mup4*F(635)+p3mup4*F(636))+p4mup4*F(1367))-p1mup4*
     -   p4mup2*F(1368)+p4mup2*p4mup4*F(1370)-p3mup4*p4mup2*F(2396)+2
     -   *(p3mup2*p3mup4*F(2397)+p1mup4*p3mup2*F(2398)+p2mup4*p3mup2*
     -   F(2399))+p3mup2*p4mup4*F(2415))
       Fa10(972) = -(p2mup1*p3mup2*p4mup4*F(362))+p1mup4*(p2mup1*p3mup
     -   2*F(346)+2*(p3mup1*p4mup2*F(359)-p1mup2*p3mup1*F(938)))+p1mu
     -   p2*p3mup1*(2*p4mup4*F(359)-p2mup4*F(939))-p2mup4*p3mup1*p4mu
     -   p2*F(946)+p2mup1*p2mup4*p3mup2*F(2318)+p2mup1*p3mup4*p4mup2*
     -   F(2369)+2*p3mup1*p3mup2*p3mup4*F(2409)
       Fa10(973) = mup1mup3*(p1mup4*p3mup2*F(346)-p3mup2*p4mup4*F(362)
     -   +p1mup2*p3mup4*F(1262)+p2mup4*p3mup2*F(2318)+p3mup4*p4mup2*F
     -   (2369)+p3mup2*p3mup4*F(2413))
       Fa10(974) = mup1mup3*(p2mup4*p4mup2*F(1221)+p1mup2*(-4*p2mup4*F
     -   (939)+8*(p4mup4*F(359)-p3mup4*F(2408)))-8*(p4mup2*p4mup4*F(3
     -   65)-p1mup4*(p4mup2*F(359)-p1mup2*F(938)+p3mup2*F(944))-p3mup
     -   2*p4mup4*F(949)-p3mup2*p3mup4*F(2409)+p3mup4*p4mup2*F(2411))
     -   +p2mup4*p3mup2*F(2414))
       Fa10(975) = 4*(p2mup1*p2mup3*p3mup2*F(1265)+p1mup3*(p2mup1*(p3m
     -   up2*F(1263)-2*(p1mup2*F(633)+p4mup2*F(1276)))+2*p3mup1*(p3mu
     -   p2*F(279)-p1mup2*F(882)+p4mup2*F(1503)))-p1mup2*(p2mup1*p3mu
     -   p3*F(1261)+p2mup3*p3mup1*F(1499)-2*(p3mup1*p3mup3*F(279)-p2m
     -   up1*p2mup3*F(635)-p2mup1*p4mup3*F(1276)+p3mup1*p4mup3*F(1503
     -   ))))+p2mup1*p3mup2*p3mup3*F(1646)+p2mup3*p3mup1*p3mup2*F(177
     -   8)+p2mup1*p3mup2*p4mup3*F(2420)+p2mup1*p3mup3*p4mup2*F(2422)
     -   +p2mup3*p3mup1*p4mup2*F(2425)+8*(-(p3mup1*p3mup2*p3mup3*F(15
     -   07))+p3mup1*p3mup3*p4mup2*F(1508)+p3mup1*p3mup2*p4mup3*F(150
     -   8)-p2mup1*p2mup3*p4mup2*F(2421)-p2mup1*p4mup2*p4mup3*F(2423)
     -   -p3mup1*p4mup2*p4mup3*F(2426))
       Fa10(976) = p3mup1*p4mup5*(p1mup3*F(258)-p3mup3*F(282))
       Fa10(977) = p3mup1*(p1mup3*p4mup5*F(297)+p1mup3*p2mup5*F(1267)+
     -   2*(p3mup3*p4mup5*F(306)+p2mup5*p3mup3*F(1304)))
       Fa10(978) = p4mup1*(4*(p1mup4*p3mup5*F(346)-p3mup5*p4mup4*F(362
     -   ))+mup4mup5*F(2320))
       Fa10(979) = -4*p4mup1*(p1mup4*p3mup5*F(272)+p1mup4*p2mup5*F(641
     -   )+2*(p3mup5*p4mup4*F(308)+p2mup5*p4mup4*F(691)))+mup4mup5*p4
     -   mup1*F(2200)
       Fa10(980) = p4mup1*(-4*(p1mup4*p2mup5*F(346)-p2mup5*p4mup4*F(36
     -   2))+mup4mup5*F(2305))
       Fa10(981) = p4mup1*(mup3mup5*(p1mup4*F(346)-p4mup4*F(362)+p3mup
     -   4*F(2188)+p2mup4*F(2318))+mup4mup5*(-(p1mup3*F(258))+p3mup3*
     -   F(282)-p2mup3*F(2338)+p4mup3*F(2345)))
       Fa10(982) = mup2mup5*p4mup1*(p1mup4*F(346)-p4mup4*F(362)+p3mup4
     -   *F(2188)+p2mup4*F(2318))
       Fa10(983) = p3mup1*(mup3mup5*(p1mup4*F(346)-p4mup4*F(362)+p3mup
     -   4*F(2188)+p2mup4*F(2318))+mup4mup5*(-(p1mup3*F(258))+p3mup3*
     -   F(282)-p2mup3*F(2338)+p4mup3*F(2345)))
       Fa10(984) = mup2mup5*p3mup1*(p1mup4*F(346)-p4mup4*F(362)+p3mup4
     -   *F(2188)+p2mup4*F(2318))
       Fa10(985) = mup2mup5*p4mup1*(p1mup3*F(258)-p3mup3*F(282)+p2mup3
     -   *F(2338)-p4mup3*F(2345))+mup1mup5*(-(p1mup3*p4mup2*F(258))+p
     -   3mup3*p4mup2*F(282)+p1mup2*p4mup3*F(1262)+p3mup2*p4mup3*F(23
     -   68))
       Fa10(986) = p2mup5*p3mup1*(p1mup3*F(258)-p3mup3*F(282))
       Fa10(987) = mup2mup5*p3mup1*(p1mup3*F(258)-p3mup3*F(282)+p2mup3
     -   *F(2338)-p4mup3*F(2345))+mup1mup5*(-(p1mup3*p3mup2*F(1263))-
     -   p2mup3*p3mup2*F(1265)+2*(p1mup3*p4mup2*F(1276)+p1mup2*(p2mup
     -   3*F(635)+p4mup3*F(1276)))+p3mup2*p4mup3*F(2428)+p3mup3*p4mup
     -   2*F(2429))
       Fa10(988) = mup1mup5*(p1mup2*p2mup3*F(1290)+p1mup3*p3mup2*F(130
     -   1)-2*(p2mup3*p3mup2*F(1299)+p3mup3*p4mup2*F(1307)+p3mup2*p4m
     -   up3*F(1307))-p1mup2*p4mup3*F(2418)+p1mup3*p4mup2*F(2424))
       Fa10(989) = p2mup1*(mup3mup5*(p1mup4*F(346)-p4mup4*F(362)+p3mup
     -   4*F(2188)+p2mup4*F(2318))+mup4mup5*(-(p1mup3*F(258))+p3mup3*
     -   F(282)-p2mup3*F(2338)+p4mup3*F(2345)))
       Fa10(990) = mup1mup5*(-(p1mup2*p2mup3*F(1499))-p2mup3*p3mup2*F(
     -   1505)+2*(p1mup2*p4mup3*F(1503)+p1mup3*(p3mup2*F(279)+p4mup2*
     -   F(1503))+(p3mup3*p4mup2+p3mup2*p4mup3)*F(1508)))+mup2mup5*p2
     -   mup1*(p1mup3*F(258)-p3mup3*F(282)+p2mup3*F(2338)-p4mup3*F(23
     -   45))
       Return
       End
