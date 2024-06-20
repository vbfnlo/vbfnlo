c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenOddG14(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup
     -   1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3
     -   mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3
     -   ,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,
     -   mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa14)
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
       Complex*16 Fa14(3797:4088)
       COMMON/FaPenOddG/Fa
       Fa14(3797) = mup2mup4*F(5164)+2*Fa(1850)
       Fa14(3798) = mup2mup4*F(5124)+2*Fa(1851)
       Fa14(3799) = mup2mup4*F(5125)+2*Fa(1853)
       Fa14(3800) = -2*Fa(1889)+Fa(1982)
       Fa14(3801) = 2*p1mup1*p4mup5*F(352)+Fa(1983)+4*(p1mup1*p1mup5*F
     -   (276)+Fa(2034))
       Fa14(3802) = mup1mup5*F(5125)-2*(p1mup1*p3mup5*F(7441)-Fa(1890)
     -   )
       Fa14(3803) = 2*Fa(1896)+Fa(1984)
       Fa14(3804) = mup1mup4*F(5128)+2*(p1mup1*p3mup4*F(7441)-Fa(1893)
     -   )
       Fa14(3805) = mup1mup2*F(5125)+4*p3mup1*p4mup2*F(7444)+2*Fa(1949
     -   )
       Fa14(3806) = mup1mup4*F(5126)+2*(p3mup4*p4mup1*F(7441)+Fa(1894)
     -   )
       Fa14(3807) = mup1mup2*F(5123)+4*Fa(1966)+2*Fa(2013)
       Fa14(3808) = mup1mup2*F(5143)+4*Fa(1967)+2*(Fa(730)-Fa(2007))
       Fa14(3809) = mup1mup4*F(5141)-2*(p2mup4*p4mup1*F(7441)+Fa(1895)
     -   )
       Fa14(3810) = mup1mup4*F(5142)+2*Fa(1955)
       Fa14(3811) = mup1mup2*F(5144)+4*Fa(1971)-2*(Fa(731)+Fa(2008))
       Fa14(3812) = mup2mup5*F(5048)+Fa(1834)
       Fa14(3813) = mup2mup5*F(5358)+Fa(1835)
       Fa14(3814) = mup2mup5*F(5359)+2*Fa(1836)
       Fa14(3815) = mup2mup5*F(5426)+Fa(1837)
       Fa14(3816) = mup2mup5*F(5431)+Fa(1838)
       Fa14(3817) = mup2mup5*F(5432)-2*Fa(1819)
       Fa14(3818) = mup2mup4*F(5437)+Fa(1867)
       Fa14(3819) = mup2mup4*F(5438)-2*Fa(1868)
       Fa14(3820) = mup2mup4*F(5439)+Fa(1869)
       Fa14(3821) = mup2mup4*F(5444)+Fa(1870)
       Fa14(3822) = mup2mup4*F(5445)+Fa(1871)
       Fa14(3823) = mup2mup4*F(5446)+Fa(1872)
       Fa14(3824) = Fa(1910)+Fa(1985)
       Fa14(3825) = p1mup1*p4mup5*F(5704)+Fa(1986)+4*(p1mup1*p1mup5*F(
     -   5454)+Fa(2023))
       Fa14(3826) = mup1mup5*F(5056)+4*(p1mup1*p1mup5*F(5447)+p1mup1*p
     -   4mup5*F(5448))-2*Fa(1911)+Fa(2039)
       Fa14(3827) = mup1mup5*F(5057)+Fa(1912)
       Fa14(3828) = -4*p1mup1*(p1mup5*F(4528)+p4mup5*F(4668))+mup1mup5
     -   *F(5058)+Fa(1913)-2*Fa(2040)
       Fa14(3829) = -2*mup4mup5*p4mup1*F(5708)+Fa(1991)+Fa(2041)
       Fa14(3830) = 2*mup4mup5*p3mup1*F(5708)+Fa(1992)+Fa(2042)
       Fa14(3831) = 2*(mup4mup5*p1mup1*F(5708)-Fa(1998))+Fa(2045)
       Fa14(3832) = Fa(1921)+Fa(1999)
       Fa14(3833) = Fa(2000)-2*(mup4mup5*p2mup1*F(5708)+Fa(2046))
       Fa14(3834) = mup1mup4*F(5050)+Fa(1914)
       Fa14(3835) = mup1mup4*F(5052)+Fa(1915)+Fa(1987)
       Fa14(3836) = mup1mup4*F(5055)-2*Fa(1916)+Fa(1988)
       Fa14(3837) = mup1mup2*F(5057)-4*p4mup2*(p2mup1*F(5662)+p3mup1*F
     -   (5668))-2*Fa(1989)
       Fa14(3838) = mup1mup2*F(5056)-4*(p3mup1*p4mup2*F(5657)-p1mup1*p
     -   4mup2*F(5662))-2*Fa(1929)
       Fa14(3839) = mup1mup2*F(5058)+2*Fa(1930)+4*Fa(2017)
       Fa14(3840) = mup1mup4*F(5048)-2*Fa(1917)+Fa(1990)
       Fa14(3841) = mup1mup2*F(4978)+4*(p4mup1*p4mup2*F(5662)+Fa(1993)
     -   )+Fa(2014)
       Fa14(3842) = mup1mup2*F(5053)+4*Fa(1994)+Fa(2015)+Fa(2021)
       Fa14(3843) = mup1mup4*F(5049)+Fa(1920)+Fa(1995)
       Fa14(3844) = mup1mup4*F(5051)-2*Fa(1940)
       Fa14(3845) = mup1mup2*F(5054)+4*Fa(2001)+Fa(2016)+Fa(2022)
       Fa14(3846) = -4*mup1mup5*p4mup2*F(5662)-2*mup2mup5*p3mup1*F(570
     -   8)+Fa(1918)+Fa(2043)
       Fa14(3847) = 2*mup1mup5*p4mup2*F(5668)+Fa(1919)+Fa(2044)
       Fa14(3848) = -2*(mup1mup5*p4mup2*F(5653)-mup2mup5*p4mup1*F(5708
     -   ))+Fa(1922)+4*Fa(2047)
       Fa14(3849) = Fa(1923)-2*(2*mup1mup5*p4mup2*F(5657)+Fa(2033))
       Fa14(3850) = Fa(1996)+2*Fa(2019)
       Fa14(3851) = 4*mup1mup4*p4mup2*F(5662)+Fa(1997)
       Fa14(3852) = 4*mup1mup4*p4mup2*F(5657)+Fa(2002)
       Fa14(3853) = mup2mup5*F(5709)-2*Fa(2212)
       Fa14(3854) = mup2mup5*F(5710)+Fa(2213)
       Fa14(3855) = mup2mup5*F(5711)+Fa(2214)
       Fa14(3856) = mup2mup5*F(5809)+Fa(2215)
       Fa14(3857) = mup2mup5*F(5810)+Fa(2216)
       Fa14(3858) = mup2mup5*F(5819)+Fa(2217)
       Fa14(3859) = mup2mup3*F(5828)+4*Fa(2238)+Fa(2250)
       Fa14(3860) = mup2mup3*F(5829)+4*Fa(2239)+Fa(2251)
       Fa14(3861) = mup2mup3*F(5830)+4*Fa(2240)+Fa(2252)
       Fa14(3862) = mup2mup3*F(5835)+Fa(2241)+Fa(2253)
       Fa14(3863) = mup2mup3*F(5836)+Fa(2242)+Fa(2254)
       Fa14(3864) = mup2mup3*F(5837)+Fa(2243)+Fa(2255)
       Fa14(3865) = Fa(2291)+Fa(2321)
       Fa14(3866) = mup1mup5*F(5778)+4*(p1mup1*p1mup5*F(5454)+p1mup1*p
     -   4mup5*F(5841))+Fa(2292)+Fa(2414)
       Fa14(3867) = mup1mup5*F(5776)-2*Fa(2293)
       Fa14(3868) = mup1mup5*F(5781)+2*p1mup1*p2mup5*F(6350)+Fa(2294)
       Fa14(3869) = Fa(2295)+Fa(2329)
       Fa14(3870) = 2*mup3mup5*p4mup1*F(6350)+Fa(2296)+Fa(2415)
       Fa14(3871) = -2*Fa(2297)+Fa(2416)
       Fa14(3872) = Fa(2299)+Fa(2419)
       Fa14(3873) = 2*p1mup1*p4mup5*F(5798)+Fa(2322)+4*(p1mup1*p1mup5*
     -   F(5795)+Fa(2413))
       Fa14(3874) = mup1mup3*F(5782)+4*p2mup1*p4mup3*F(6167)+Fa(2323)
       Fa14(3875) = mup1mup3*F(5783)+Fa(2324)-4*Fa(2367)
       Fa14(3876) = mup1mup3*F(5784)+4*p2mup1*p4mup3*F(5804)+Fa(2325)
       Fa14(3877) = mup1mup2*F(5776)+4*p4mup2*(p3mup1*F(5804)+p2mup1*F
     -   (6384))+Fa(2326)
       Fa14(3878) = mup1mup2*F(5778)+4*p4mup2*(p3mup1*F(5804)+p1mup1*F
     -   (6348))+Fa(2327)
       Fa14(3879) = mup1mup2*F(5781)-4*F(5804)*Fa(732)+Fa(2328)
       Fa14(3880) = mup1mup3*F(5709)+4*Fa(2330)+Fa(2368)
       Fa14(3881) = mup1mup2*F(5775)+4*Fa(2331)+Fa(2369)+Fa(2400)
       Fa14(3882) = -2*mup2mup5*p4mup1*F(6350)+mup1mup5*p4mup2*F(6382)
     -   +Fa(2298)+4*Fa(2417)
       Fa14(3883) = 4*mup1mup5*p4mup2*F(5804)+2*mup1mup2*p2mup5*F(6350
     -   )+Fa(2418)
       Fa14(3884) = mup1mup3*F(5779)-4*Fa(2332)+Fa(2370)
       Fa14(3885) = mup1mup2*F(5774)+4*(p4mup1*p4mup2*F(6348)-Fa(2333)
     -   )+Fa(2371)
       Fa14(3886) = mup1mup3*F(5780)+4*Fa(2337)+Fa(2373)
       Fa14(3887) = mup1mup2*F(5777)-4*Fa(2338)+Fa(2374)+Fa(2402)
       Fa14(3888) = 4*mup1mup5*p4mup2*F(5804)+Fa(2300)+Fa(2420)
       Fa14(3889) = 2*mup2mup5*p3mup1*F(6350)+4*mup1mup5*p4mup2*F(6384
     -   )+Fa(2301)+Fa(2421)
       Fa14(3890) = mup1mup2*p4mup3*F(5800)+4*Fa(2334)+Fa(2401)
       Fa14(3891) = 4*mup1mup3*p4mup2*F(6348)+Fa(2335)+Fa(2372)
       Fa14(3892) = -4*mup1mup3*p4mup2*F(5804)+2*mup2mup3*p2mup1*F(635
     -   0)+Fa(2336)
       Fa14(3893) = 4*(-(mup1mup3*p4mup2)+mup1mup2*p4mup3)*F(5804)+Fa(
     -   2339)
       Fa14(3894) = mup2mup5*F(5844)+2*Fa(2218)
       Fa14(3895) = mup2mup5*F(5845)-2*Fa(2219)
       Fa14(3896) = mup2mup5*F(5846)+2*Fa(2220)
       Fa14(3897) = -2*Fa(2221)+mup2mup5*F(1051)*P(23)
       Fa14(3898) = mup2mup5*F(5885)-2*Fa(2222)
       Fa14(3899) = mup2mup3*F(5847)-2*Fa(2260)
       Fa14(3900) = -2*Fa(2261)+mup2mup3*F(1051)*P(7)
       Fa14(3901) = mup2mup3*F(5848)+2*Fa(2262)
       Fa14(3902) = mup2mup3*F(5896)+2*Fa(216)+4*Fa(2263)
       Fa14(3903) = mup2mup3*F(5897)+2*Fa(217)+4*Fa(2264)
       Fa14(3904) = -2*p3mup1*p4mup5*F(1051)+Fa(2340)
       Fa14(3905) = -4*p1mup1*(p1mup5*F(276)+p4mup5*F(1012))+mup1mup5*
     -   F(5849)-2*(Fa(358)-Fa(840))
       Fa14(3906) = mup1mup5*F(5863)+2*Fa(2302)
       Fa14(3907) = 2*p2mup1*p4mup5+Fa(2342)
       Fa14(3908) = -(mup3mup5*p4mup1)+mup1mup5*p4mup3+mup1mup3*p4mup5
       Fa14(3909) = Fa(2303)-Fa(2422)
       Fa14(3910) = -2*p1mup1*p4mup5*F(1051)+Fa(2341)
       Fa14(3911) = -4*p2mup1*p4mup3*F(1012)+mup1mup3*F(5865)-2*Fa(494
     -   )
       Fa14(3912) = mup1mup3*F(5864)+2*Fa(495)+4*Fa(2375)
       Fa14(3913) = -4*(p3mup1*p4mup2*F(335)-p2mup1*p4mup2*F(2718))+mu
     -   p1mup2*F(5863)+2*Fa(500)
       Fa14(3914) = mup1mup2*F(5849)-4*(p3mup1*p4mup2*F(335)-p1mup1*p4
     -   mup2*F(6355))-2*Fa(501)
       Fa14(3915) = mup1mup3*F(5844)+2*Fa(2376)
       Fa14(3916) = 2*(Fa(2377)+Fa(2403))+mup1mup2*F(1051)*P(7)
       Fa14(3917) = -(mup2mup5*p4mup1)+mup1mup5*p4mup2+mup1mup2*p4mup5
       Fa14(3918) = 2*Fa(2378)+mup1mup3*P(23)
       Fa14(3919) = mup1mup2*F(5847)+4*p4mup1*p4mup2*F(6355)+2*Fa(2379
     -   )
       Fa14(3920) = mup1mup3*F(5845)-2*Fa(2381)
       Fa14(3921) = mup1mup2*F(5848)+2*(Fa(2382)+Fa(2405))
       Fa14(3922) = -(mup2mup5*p3mup1*F(1051))+2*mup1mup5*p4mup2*F(271
     -   8)+Fa(814)-Fa(2304)
       Fa14(3923) = -(mup1mup2*p4mup3)+Fa(2404)
       Fa14(3924) = 2*mup1mup3*p4mup2*F(6355)+Fa(551)+Fa(2380)
       Fa14(3925) = mup2mup5*F(5898)+2*Fa(2223)
       Fa14(3926) = mup2mup5*F(5899)-2*Fa(2224)
       Fa14(3927) = mup2mup5*F(5900)-2*Fa(2225)
       Fa14(3928) = -2*Fa(2226)+mup2mup5*F(7440)*P(23)
       Fa14(3929) = mup2mup5*F(5939)-2*Fa(2227)
       Fa14(3930) = mup2mup3*F(5901)-2*Fa(2267)
       Fa14(3931) = -2*Fa(2268)+mup2mup3*F(7440)*P(7)
       Fa14(3932) = mup2mup3*F(5902)+2*Fa(2269)
       Fa14(3933) = mup2mup3*F(5950)+2*Fa(226)+4*Fa(2270)
       Fa14(3934) = mup2mup3*F(5951)+2*Fa(227)+4*Fa(2271)
       Fa14(3935) = 2*p3mup1*p4mup5*F(7440)+Fa(2343)
       Fa14(3936) = 4*(p1mup1*p1mup5*F(276)+p1mup1*p4mup5*F(1012))-mup
     -   1mup5*F(5917)+2*(Fa(375)+Fa(832))
       Fa14(3937) = mup1mup5*F(5916)+2*Fa(2305)
       Fa14(3938) = Fa(2306)-Fa(2423)
       Fa14(3939) = -2*p1mup1*p4mup5*F(7440)+Fa(2344)
       Fa14(3940) = mup1mup3*F(5918)-4*p2mup1*p4mup3*F(7446)-2*Fa(521)
       Fa14(3941) = mup1mup3*F(5919)+2*Fa(522)+4*Fa(2383)
       Fa14(3942) = mup1mup2*F(5916)+4*(p2mup1*p4mup2*F(2718)-p3mup1*p
     -   4mup2*F(7442))+2*Fa(524)
       Fa14(3943) = mup1mup2*F(5917)+4*(p1mup1*p4mup2*F(6355)-p3mup1*p
     -   4mup2*F(7442))-2*Fa(525)
       Fa14(3944) = mup1mup3*F(5898)+2*Fa(2384)
       Fa14(3945) = 2*(Fa(2385)-Fa(2406))+mup1mup2*F(7440)*P(7)
       Fa14(3946) = mup1mup2*F(5901)+4*p4mup1*p4mup2*F(6355)+2*Fa(2386
     -   )
       Fa14(3947) = mup1mup3*F(5899)-2*Fa(2388)
       Fa14(3948) = mup1mup2*F(5902)+2*(Fa(2389)+Fa(2407))
       Fa14(3949) = -2*mup1mup5*p4mup2*F(2718)+mup2mup5*p3mup1*F(7440)
     -   -Fa(826)+Fa(2307)
       Fa14(3950) = 2*mup1mup3*p4mup2*F(6355)+Fa(569)+Fa(2387)
       Fa14(3951) = mup2mup5*F(5952)-2*Fa(2228)
       Fa14(3952) = mup2mup5*F(5953)+2*Fa(2219)
       Fa14(3953) = mup2mup5*F(5954)+2*Fa(2229)
       Fa14(3954) = 2*Fa(2221)+mup2mup5*F(1051)*P(7)
       Fa14(3955) = mup2mup5*F(5990)+2*Fa(2230)
       Fa14(3956) = mup2mup3*F(5955)+2*Fa(2273)
       Fa14(3957) = 2*Fa(2261)+mup2mup3*F(1051)*P(23)
       Fa14(3958) = mup2mup3*F(5956)+2*Fa(2274)
       Fa14(3959) = mup2mup3*F(6000)-2*Fa(216)-4*Fa(2263)
       Fa14(3960) = mup2mup3*F(6001)+2*Fa(235)+4*Fa(2275)
       Fa14(3961) = 2*p3mup1*p4mup5*F(1051)+Fa(2345)
       Fa14(3962) = 4*(p1mup1*p1mup5*F(276)+p1mup1*p4mup5*F(1012))+mup
     -   1mup5*F(5957)+2*(Fa(358)+Fa(823))
       Fa14(3963) = mup1mup5*F(5865)-2*Fa(2308)
       Fa14(3964) = -2*p2mup1*p4mup5+Fa(2347)
       Fa14(3965) = -Fa(2309)+Fa(2422)
       Fa14(3966) = -2*p1mup1*p4mup5*F(1051)+Fa(2346)
       Fa14(3967) = 4*p2mup1*p4mup3*F(1012)+mup1mup3*F(5863)+2*Fa(494)
       Fa14(3968) = mup1mup3*F(5970)-2*Fa(495)-4*Fa(2375)
       Fa14(3969) = mup1mup2*F(5865)+4*(p3mup1*p4mup2*F(335)+p2mup1*p4
     -   mup2*F(6355))+2*Fa(543)
       Fa14(3970) = 4*(p3mup1*p4mup2*F(335)+p1mup1*p4mup2*F(2718))+mup
     -   1mup2*F(5957)+2*Fa(544)
       Fa14(3971) = mup1mup3*F(5952)-2*Fa(2376)
       Fa14(3972) = -2*(Fa(2377)-Fa(2408))+mup1mup2*F(1051)*P(23)
       Fa14(3973) = -2*Fa(2378)+mup1mup3*P(7)
       Fa14(3974) = 4*p4mup1*p4mup2*F(2718)+mup1mup2*F(5955)+2*Fa(2390
     -   )
       Fa14(3975) = mup1mup3*F(5953)+2*Fa(2381)
       Fa14(3976) = mup1mup2*F(5956)-2*(Fa(2382)+Fa(2405))
       Fa14(3977) = mup2mup5*p3mup1*F(1051)+2*mup1mup5*p4mup2*F(6355)+
     -   Fa(834)+Fa(2310)
       Fa14(3978) = mup1mup2*p4mup3+Fa(2409)
       Fa14(3979) = 2*mup1mup3*p4mup2*F(2718)+Fa(510)-Fa(2380)
       Fa14(3980) = mup2mup5*F(5955)+2*Fa(2218)
       Fa14(3981) = mup2mup5*F(5956)-2*Fa(2219)
       Fa14(3982) = mup2mup5*F(6002)+2*Fa(2220)
       Fa14(3983) = mup2mup3*F(5952)-2*Fa(2260)
       Fa14(3984) = mup2mup3*F(5953)+2*Fa(2262)
       Fa14(3985) = mup2mup3*F(6032)+2*Fa(217)+4*Fa(2264)
       Fa14(3986) = -2*p3mup1*p4mup5*F(1051)+Fa(2348)
       Fa14(3987) = -4*p1mup1*(p1mup5*F(276)+p4mup5*F(1012))+mup1mup5*
     -   F(5970)-2*(Fa(358)-Fa(840))
       Fa14(3988) = -2*p1mup1*p4mup5*F(1051)+Fa(2349)
       Fa14(3989) = mup1mup3*F(5957)+2*Fa(495)+4*Fa(2375)
       Fa14(3990) = mup1mup2*F(5970)-4*(p3mup1*p4mup2*F(335)-p1mup1*p4
     -   mup2*F(6355))-2*Fa(501)
       Fa14(3991) = mup1mup3*F(5955)+2*Fa(2376)
       Fa14(3992) = mup1mup2*F(5952)+4*p4mup1*p4mup2*F(6355)+2*Fa(2379
     -   )
       Fa14(3993) = mup1mup3*F(5956)-2*Fa(2381)
       Fa14(3994) = mup1mup2*F(5953)+2*(Fa(2382)+Fa(2405))
       Fa14(3995) = mup2mup5*F(5901)-2*Fa(2231)
       Fa14(3996) = mup2mup5*F(5902)+2*Fa(2224)
       Fa14(3997) = mup2mup5*F(5951)+2*Fa(2225)
       Fa14(3998) = 2*Fa(2226)+mup2mup5*F(7440)*P(7)
       Fa14(3999) = mup2mup5*F(5950)+2*Fa(2227)
       Fa14(4000) = mup2mup3*F(5898)+2*Fa(2278)
       Fa14(4001) = 2*Fa(2268)+mup2mup3*F(7440)*P(23)
       Fa14(4002) = mup2mup3*F(5899)+2*Fa(2279)
       Fa14(4003) = mup2mup3*F(5939)-2*Fa(226)-4*Fa(2270)
       Fa14(4004) = mup2mup3*F(5900)-2*Fa(227)-4*Fa(2271)
       Fa14(4005) = 2*p3mup1*p4mup5*F(7440)+Fa(2350)
       Fa14(4006) = 4*(p1mup1*p1mup5*F(276)+p1mup1*p4mup5*F(1012))+mup
     -   1mup5*F(5919)+2*(Fa(375)+Fa(832))
       Fa14(4007) = mup1mup5*F(5918)-2*Fa(2305)
       Fa14(4008) = -Fa(2306)+Fa(2423)
       Fa14(4009) = 2*p1mup1*p4mup5*F(7440)-Fa(2351)
       Fa14(4010) = mup1mup3*F(5916)+4*p2mup1*p4mup3*F(7446)+2*Fa(521)
       Fa14(4011) = mup1mup3*F(5917)-2*Fa(522)-4*Fa(2383)
       Fa14(4012) = mup1mup2*F(5918)+4*p4mup2*(p2mup1*F(6355)+p3mup1*F
     -   (7442))+2*Fa(567)
       Fa14(4013) = mup1mup2*F(5919)+4*(p1mup1*p4mup2*F(2718)+p3mup1*p
     -   4mup2*F(7442))+2*Fa(568)
       Fa14(4014) = mup1mup3*F(5901)-2*Fa(2384)
       Fa14(4015) = -2*(Fa(2385)-Fa(2406))+mup1mup2*F(7440)*P(23)
       Fa14(4016) = 4*p4mup1*p4mup2*F(2718)+mup1mup2*F(5898)+2*Fa(2391
     -   )
       Fa14(4017) = mup1mup3*F(5902)+2*Fa(2388)
       Fa14(4018) = mup1mup2*F(5899)-2*(Fa(2389)+Fa(2407))
       Fa14(4019) = 2*mup1mup5*p4mup2*F(6355)+mup2mup5*p3mup1*F(7440)+
     -   Fa(842)+Fa(2307)
       Fa14(4020) = 2*mup1mup3*p4mup2*F(2718)+Fa(529)-Fa(2387)
       Fa14(4021) = mup2mup5*F(5847)-2*Fa(2228)
       Fa14(4022) = mup2mup5*F(5848)+2*Fa(2219)
       Fa14(4023) = mup2mup5*F(6059)+2*Fa(2229)
       Fa14(4024) = mup2mup3*F(5844)+2*Fa(2273)
       Fa14(4025) = mup2mup3*F(5845)+2*Fa(2274)
       Fa14(4026) = mup2mup3*F(6089)+2*Fa(235)+4*Fa(2275)
       Fa14(4027) = 2*p3mup1*p4mup5*F(1051)+Fa(2352)
       Fa14(4028) = 4*(p1mup1*p1mup5*F(276)+p1mup1*p4mup5*F(1012))+mup
     -   1mup5*F(5864)+2*(Fa(358)+Fa(823))
       Fa14(4029) = -2*p1mup1*p4mup5*F(1051)+Fa(2353)
       Fa14(4030) = mup1mup3*F(5849)-2*Fa(495)-4*Fa(2375)
       Fa14(4031) = 4*(p3mup1*p4mup2*F(335)+p1mup1*p4mup2*F(2718))+mup
     -   1mup2*F(5864)+2*Fa(544)
       Fa14(4032) = mup1mup3*F(5847)-2*Fa(2376)
       Fa14(4033) = 4*p4mup1*p4mup2*F(2718)+mup1mup2*F(5844)+2*Fa(2390
     -   )
       Fa14(4034) = mup1mup3*F(5848)+2*Fa(2381)
       Fa14(4035) = mup1mup2*F(5845)-2*(Fa(2382)+Fa(2405))
       Fa14(4036) = mup2mup5*F(5774)+2*Fa(2232)
       Fa14(4037) = mup2mup5*F(6090)+Fa(2233)
       Fa14(4038) = mup2mup5*F(6091)-2*Fa(2234)
       Fa14(4039) = mup2mup5*F(6139)+Fa(2235)
       Fa14(4040) = mup2mup5*F(5837)+Fa(2236)
       Fa14(4041) = mup2mup5*F(6148)-2*Fa(2237)
       Fa14(4042) = mup2mup3*F(6157)+4*Fa(2244)+Fa(2281)
       Fa14(4043) = mup2mup3*F(6158)+Fa(2245)+Fa(2282)
       Fa14(4044) = mup2mup3*F(6159)+4*Fa(2246)+Fa(2283)
       Fa14(4045) = mup2mup3*F(6164)+Fa(2247)+Fa(2284)
       Fa14(4046) = mup2mup3*F(6165)+Fa(2248)+Fa(2285)
       Fa14(4047) = mup2mup3*F(5810)+Fa(2249)+Fa(2286)
       Fa14(4048) = Fa(2311)+Fa(2354)
       Fa14(4049) = mup1mup5*F(5783)-4*p1mup1*(p1mup5*F(5454)+p4mup5*F
     -   (5841))+Fa(2312)-2*Fa(2425)
       Fa14(4050) = mup1mup5*F(5782)+2*Fa(2293)
       Fa14(4051) = mup1mup5*F(5784)-2*p1mup1*p2mup5*F(6350)+Fa(2313)
       Fa14(4052) = Fa(2314)+Fa(2362)
       Fa14(4053) = -2*mup3mup5*p4mup1*F(6350)+Fa(2315)+Fa(2426)
       Fa14(4054) = Fa(2316)+Fa(2427)
       Fa14(4055) = Fa(2318)+Fa(2429)
       Fa14(4056) = -2*p1mup1*p4mup5*F(5798)+Fa(2355)+4*(p1mup1*p1mup5
     -   *F(6131)+Fa(2424))
       Fa14(4057) = mup1mup3*F(5776)+4*p2mup1*p4mup3*F(5843)+Fa(2356)
       Fa14(4058) = mup1mup3*F(5778)+Fa(2357)+4*Fa(2392)
       Fa14(4059) = mup1mup3*F(5781)-4*p2mup1*p4mup3*F(5804)+Fa(2358)
       Fa14(4060) = mup1mup2*F(5782)-4*(p3mup1*p4mup2*F(5804)-p2mup1*p
     -   4mup2*F(6348))+Fa(2359)
       Fa14(4061) = mup1mup2*F(5783)-4*(p3mup1*p4mup2*F(5804)-p1mup1*p
     -   4mup2*F(6384))+Fa(2360)
       Fa14(4062) = mup1mup2*F(5784)+4*F(5804)*Fa(732)+Fa(2361)
       Fa14(4063) = mup1mup3*F(5774)-4*Fa(2330)+Fa(2393)
       Fa14(4064) = mup1mup2*F(5779)-4*Fa(2331)+Fa(2394)+Fa(2410)
       Fa14(4065) = mup1mup5*p4mup2*F(6344)+2*mup2mup5*p4mup1*F(6350)+
     -   Fa(2317)-4*Fa(2417)
       Fa14(4066) = -4*mup1mup5*p4mup2*F(5804)-2*mup1mup2*p2mup5*F(635
     -   0)+Fa(2428)
       Fa14(4067) = mup1mup3*F(5775)+4*Fa(2332)+Fa(2395)
       Fa14(4068) = mup1mup2*F(5709)+4*(p4mup1*p4mup2*F(6384)+Fa(2333)
     -   )+Fa(2396)
       Fa14(4069) = mup1mup3*F(5777)-4*Fa(2337)+Fa(2398)
       Fa14(4070) = mup1mup2*F(5780)+4*Fa(2338)+Fa(2399)+Fa(2412)
       Fa14(4071) = -4*mup1mup5*p4mup2*F(5804)+Fa(2319)+Fa(2430)
       Fa14(4072) = 4*mup1mup5*p4mup2*F(6348)-2*(mup2mup5*p3mup1*F(635
     -   0)+Fa(2320))+Fa(2431)
       Fa14(4073) = mup1mup2*p4mup3*F(6134)+4*Fa(2363)+Fa(2411)
       Fa14(4074) = 4*mup1mup3*p4mup2*F(6384)+Fa(2364)+Fa(2397)
       Fa14(4075) = 4*mup1mup3*p4mup2*F(5804)-2*mup2mup3*p2mup1*F(6350
     -   )+Fa(2365)
       Fa14(4076) = 4*(mup1mup3*p4mup2-mup1mup2*p4mup3)*F(5804)+Fa(236
     -   6)
       Fa14(4077) = mup3mup4*F(6421)+Fa(2560)
       Fa14(4078) = mup3mup4*F(6422)+Fa(2561)
       Fa14(4079) = mup3mup4*F(6423)+Fa(2562)
       Fa14(4080) = mup3mup4*F(6392)+Fa(2563)
       Fa14(4081) = mup3mup4*F(6395)+Fa(2564)
       Fa14(4082) = mup3mup4*F(6393)+Fa(2565)
       Fa14(4083) = mup2mup4*F(6401)+Fa(2597)
       Fa14(4084) = mup2mup4*F(6402)+Fa(2598)
       Fa14(4085) = mup2mup4*F(6403)+Fa(2599)
       Fa14(4086) = mup2mup4*F(6456)+Fa(2600)
       Fa14(4087) = mup2mup4*F(6457)+2*Fa(2601)
       Fa14(4088) = mup2mup4*F(6458)+Fa(2602)
       Return
       End
