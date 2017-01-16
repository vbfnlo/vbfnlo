c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenOddGSim20(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4
     -   ,p3mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5m
     -   up3,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mu
     -   p3,mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa20)
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
       Complex*16 Fa20(1882:1980)
       COMMON/FaPenOddGSim/Fa
       Fa20(1882) = -(mup2mup4*(p1mup3*F(3963)+2*(p3mup3*F(3959)+p4mup
     -   3*F(3962)+p2mup3*F(4405))))+mup2mup3*(p1mup4*F(3963)+2*(p3mu
     -   p4*F(3960)+p4mup4*F(3961)+p2mup4*F(4407)))+mup3mup4*(p1mup2*
     -   F(3963)+2*(p3mup2*F(4406)+p4mup2*F(4408)))
       Fa20(1883) = mup2mup3*p3mup4*F(3966)+4*m1sq*(mup2mup3*p1mup4*F(
     -   346)+mup3mup4*p1mup2*F(1262)+mup2mup3*p2mup4*F(2318)+mup3mup
     -   4*p4mup2*F(2369)+mup2mup3*p4mup4*F(3969))+mup2mup4*(p3mup3*F
     -   (3967)-4*(p1mup3*F(3965)+p4mup3*F(3968)+p2mup3*F(4409)))+2*m
     -   up3mup4*p3mup2*F(4410)
       Fa20(1884) = (mup2mup4*p2mup3-mup2mup3*p2mup4)*F(4412)+2*mup3mu
     -   p4*(p1mup2*F(4411)+p3mup2*F(4413)+p4mup2*F(4414))
       Fa20(1885) = 2*p1mup2*p4mup3*F(1367)+p3mup2*p4mup3*F(4417)-p4mu
     -   p2*(2*(p1mup3*F(4415)+p2mup3*F(4416))+p3mup3*F(4418))
       Fa20(1886) = -(p1mup2*p4mup3*F(4419))+p1mup3*p4mup2*F(4420)+2*(
     -   p2mup3*p4mup2*F(4421)-p3mup2*p4mup3*F(4422)+p3mup3*p4mup2*F(
     -   4423))
       Fa20(1887) = p1mup2*p4mup3*F(4424)+p1mup3*p4mup2*F(4425)+2*p2mu
     -   p3*p4mup2*F(4426)+p3mup2*p4mup3*F(4427)+p3mup3*p4mup2*F(4428
     -   )
       Fa20(1888) = p4mup2*(4*(p1mup3*F(4329)+p3mup3*F(4331))-p2mup3*F
     -   (4430))-2*p4mup3*(p1mup2*F(4429)+p3mup2*F(4431))
       Fa20(1889) = p2mup3*p4mup2*F(4433)+p3mup3*p4mup2*F(4434)+2*(p1m
     -   up3*p4mup2*F(4432)+p1mup2*p4mup3*F(4432)+p3mup2*p4mup3*F(443
     -   5))
       Fa20(1890) = p1mup3*p4mup2*F(4251)+p2mup3*p4mup2*F(4252)+2*(p1m
     -   up2*p4mup3*F(4250)+(p3mup3*p4mup2+p3mup2*p4mup3)*F(4253))
       Fa20(1891) = mup2mup4*(2*(p1mup3*F(258)+p3mup3*F(600)+p2mup3*F(
     -   2338))-p4mup3*F(4017))+mup2mup3*(p4mup4*F(4016)+2*(p1mup4*F(
     -   4014)+p3mup4*F(4015)+p2mup4*F(4436)))+mup3mup4*(2*(p1mup2*F(
     -   1608)+p3mup2*F(3506))+p4mup2*F(4437))
       Fa20(1892) = mup3mup4*(p1mup2*F(4018)+2*(p3mup2*F(4019)+p4mup2*
     -   F(4020)))-mup2mup4*(p1mup3*F(4018)+2*(p3mup3*F(4019)+p4mup3*
     -   F(4020)+p2mup3*F(4438)))+mup2mup3*(p1mup4*F(4018)+2*(p3mup4*
     -   F(4019)+p4mup4*F(4020)+p2mup4*F(4438)))
       Fa20(1893) = mup2mup4*(p3mup3*F(4023)+2*(p1mup3*F(4021)+p4mup3*
     -   F(4025)+p2mup3*F(4439)))+mup2mup3*(p3mup4*F(4022)+2*(p4mup4*
     -   F(362)+p1mup4*F(4024)+p2mup4*F(4441)))+mup3mup4*(-(p3mup2*F(
     -   4440))+2*(p1mup2*F(1608)+p4mup2*F(4442)))
       Fa20(1894) = 2*mup3mup4*(p3mup2*F(4019)+p4mup2*F(4020)+p1mup2*F
     -   (4443))+(mup2mup4*p2mup3-mup2mup3*p2mup4)*F(4444)
       Fa20(1895) = -2*(p2mup3*p4mup2*F(640)-p1mup2*p4mup3*F(1367)+p1m
     -   up3*p4mup2*F(1368))+p3mup2*p4mup3*F(4321)+p3mup3*p4mup2*F(43
     -   22)
       Fa20(1896) = -(p1mup2*p4mup3*F(4280))+p1mup3*p4mup2*F(4281)+2*(
     -   p2mup3*p4mup2*F(4323)+p3mup2*p4mup3*F(4324)+p3mup3*p4mup2*F(
     -   4325))
       Fa20(1897) = -(p1mup2*p4mup3*F(1400))+p1mup3*p4mup2*F(1401)+2*p
     -   2mup3*p4mup2*F(4145)+p3mup2*p4mup3*F(4326)+p3mup3*p4mup2*F(4
     -   327)
       Fa20(1898) = p2mup3*p4mup2*F(4328)-2*(p1mup2*p4mup3*F(4293)-p3m
     -   up2*p4mup3*F(4330))+4*(p1mup3*p4mup2*F(4329)+p3mup3*p4mup2*F
     -   (4331))
       Fa20(1899) = p2mup3*p4mup2*F(946)+2*(p1mup3*p4mup2*F(4175)+p1mu
     -   p2*p4mup3*F(4175)+p3mup2*p4mup3*F(4176))+p3mup3*p4mup2*F(433
     -   5)
       Fa20(1900) = 2*(p3mup3*p4mup2*F(317)+p3mup2*p4mup3*F(317)+p1mup
     -   2*p4mup3*F(4332))+p1mup3*p4mup2*F(4333)+p2mup3*p4mup2*F(4334
     -   )
       Fa20(1901) = mup2mup4*(2*(p1mup3*F(258)+p3mup3*F(600)+p2mup3*F(
     -   2338))+p4mup3*F(4039))+mup2mup3*(p4mup4*F(4038)+2*(p1mup4*F(
     -   264)+p3mup4*F(282)+p2mup4*F(4445)))+mup3mup4*(2*(p1mup2*F(16
     -   08)+p3mup2*F(3506))+p4mup2*F(4446))
       Fa20(1902) = mup3mup4*(p1mup2*F(4028)+2*(p3mup2*F(4035)+p4mup2*
     -   F(4041)))+mup2mup4*(p1mup3*F(4029)+2*(p3mup3*F(4034)+p4mup3*
     -   F(4040)+p2mup3*F(4447)))+mup2mup3*(p1mup4*F(4028)+2*(p3mup4*
     -   F(4035)+p4mup4*F(4041)+p2mup4*F(4448)))
       Fa20(1903) = mup2mup4*(2*(p1mup3*F(346)+p2mup3*F(2318)+p4mup3*F
     -   (3969))+p3mup3*F(4031))+mup2mup3*(p3mup4*F(4030)+2*(p4mup4*F
     -   (362)+p1mup4*F(4024)+p2mup4*F(4441)))+mup3mup4*(2*(p1mup2*F(
     -   1608)+p4mup2*F(4442))+p3mup2*F(4449))
       Fa20(1904) = mup2mup3*p2mup4*F(1056)+mup2mup4*p2mup3*F(1057)+2*
     -   mup3mup4*(p3mup2*F(4035)+p4mup2*F(4041)+p1mup2*F(4450))
       Fa20(1905) = 2*(p2mup3*p4mup2*F(640)+p1mup3*p4mup2*F(1368)+p1mu
     -   p2*p4mup3*F(4209))+p3mup2*p4mup3*F(4275)+p3mup3*p4mup2*F(427
     -   6)
       Fa20(1906) = p1mup2*p4mup3*F(4280)-p1mup3*p4mup2*F(4281)+2*(p2m
     -   up3*p4mup2*F(4282)+p3mup2*p4mup3*F(4283)+p3mup3*p4mup2*F(428
     -   4))
       Fa20(1907) = 2*p2mup3*p4mup2*F(689)+p1mup2*p4mup3*F(1400)+p1mup
     -   3*p4mup2*F(4231)+p3mup2*p4mup3*F(4290)+p3mup3*p4mup2*F(4291)
       Fa20(1908) = 4*(p1mup3*p4mup2*F(1503)+p3mup3*p4mup2*F(1508))+p2
     -   mup3*p4mup2*F(4294)+2*(p1mup2*p4mup3*F(4293)+p3mup2*p4mup3*F
     -   (4295))
       Fa20(1909) = 2*(p1mup3*p4mup2*F(359)+p1mup2*p4mup3*F(359)+p3mup
     -   2*p4mup3*F(949))+p2mup3*p4mup2*F(4260)+p3mup3*p4mup2*F(4312)
       Fa20(1910) = p1mup3*p4mup2*F(4303)+p2mup3*p4mup2*F(4305)+2*(p1m
     -   up2*p4mup3*F(4302)+(p3mup3*p4mup2+p3mup2*p4mup3)*F(4307))
       Fa20(1911) = mup2mup3*(2*(p1mup4*F(258)+p3mup4*F(600)+p2mup4*F(
     -   2338))+p4mup4*F(4039))+mup2mup4*(p4mup3*F(4038)+2*(p1mup3*F(
     -   264)+p3mup3*F(282)+p2mup3*F(4445)))+mup3mup4*(2*(p1mup2*F(12
     -   62)+p3mup2*F(2368))+p4mup2*F(4451))
       Fa20(1912) = mup3mup4*(p1mup2*F(4029)+2*(p3mup2*F(4034)+p4mup2*
     -   F(4040)))+mup2mup3*(p1mup4*F(4029)+2*(p3mup4*F(4034)+p4mup4*
     -   F(4040)+p2mup4*F(4447)))+mup2mup4*(p1mup3*F(4028)+2*(p3mup3*
     -   F(4035)+p4mup3*F(4041)+p2mup3*F(4448)))
       Fa20(1913) = mup2mup3*(2*(p1mup4*F(346)+p2mup4*F(2318)+p4mup4*F
     -   (3969))+p3mup4*F(4031))+mup2mup4*(p3mup3*F(4030)+2*(p4mup3*F
     -   (362)+p1mup3*F(4024)+p2mup3*F(4441)))+mup3mup4*(2*(p1mup2*F(
     -   1262)+p4mup2*F(2369))+p3mup2*F(4452))
       Fa20(1914) = mup2mup4*p2mup3*F(1056)+mup2mup3*p2mup4*F(1057)+2*
     -   mup3mup4*(p3mup2*F(4034)+p4mup2*F(4040)+p1mup2*F(4453))
       Fa20(1915) = 2*(p1mup2*p4mup3*F(4209)+p1mup3*p4mup2*F(4415)+p2m
     -   up3*p4mup2*F(4416))-p3mup2*p4mup3*F(4417)+p3mup3*p4mup2*F(44
     -   18)
       Fa20(1916) = p1mup2*p4mup3*F(4419)+2*p3mup2*p4mup3*F(4422)-p4mu
     -   p2*(p1mup3*F(4420)+2*(p2mup3*F(4421)+p3mup3*F(4423)))
       Fa20(1917) = p1mup2*p4mup3*F(4424)+p1mup3*p4mup2*F(4425)+2*p2mu
     -   p3*p4mup2*F(4426)+p3mup3*p4mup2*F(4428)-p3mup2*p4mup3*F(4454
     -   )
       Fa20(1918) = 4*(p1mup3*p4mup2*F(1503)+p3mup3*p4mup2*F(1508))+p2
     -   mup3*p4mup2*F(4430)+2*(p1mup2*p4mup3*F(4429)+p3mup2*p4mup3*F
     -   (4431))
       Fa20(1919) = -(mup2mup3*(p4mup4*F(4016)+2*(p1mup4*F(4014)+p3mup
     -   4*F(4015)+p2mup4*F(4436))))+mup3mup4*(2*(p1mup2*F(1262)+p3mu
     -   p2*F(2368))-p4mup2*F(4437))+mup2mup4*(p4mup3*F(4017)+2*(p1mu
     -   p3*F(264)+p3mup3*F(282)+p2mup3*F(4445)))
       Fa20(1920) = mup2mup3*(2*(p1mup4*F(346)+p2mup4*F(2318)+p4mup4*F
     -   (3969))-p3mup4*F(4022))-mup2mup4*(p3mup3*F(4023)+2*(p1mup3*F
     -   (4021)+p4mup3*F(4025)+p2mup3*F(4439)))+mup3mup4*(2*(p1mup2*F
     -   (1262)+p4mup2*F(2369))+p3mup2*F(4440))
       Fa20(1921) = 2*p3mup2*p4mup3*F(4378)+p3mup3*p4mup2*F(4455)-4*(p
     -   1mup3*p4mup2*F(4376)+p2mup3*p4mup2*F(4377)-p1mup2*p4mup3*F(4
     -   456))
       Fa20(1922) = 4*(p2mup3*p4mup2*F(4381)-p3mup2*p4mup3*F(4382)+p3m
     -   up3*p4mup2*F(4383))-2*(p1mup2*p4mup3*F(4384)-p1mup3*p4mup2*F
     -   (4385))
       Fa20(1923) = 2*p2mup3*p4mup2*F(4386)+p3mup3*p4mup2*F(4388)+p1mu
     -   p2*p4mup3*F(4389)+p1mup3*p4mup2*F(4390)+p3mup2*p4mup3*F(4457
     -   )
       Fa20(1924) = 8*m4sq*p4mup2*(p1mup3*F(4329)+p3mup3*F(4331))-2*p2
     -   mup3*p4mup2*F(4392)-4*p4mup3*(p1mup2*F(4391)+p3mup2*F(4393))
       Fa20(1925) = mup2mup4*(-2*p4mup3*F(3957)+4*(p1mup3*F(3951)+p3mu
     -   p3*F(4096)+p2mup3*F(4401)))+mup2mup3*(p3mup4*F(4095)+p4mup4*
     -   F(4097)+4*(p1mup4*F(4094)+p2mup4*F(4458)))+2*mup3mup4*(p4mup
     -   2*F(4404)-2*(p1mup2*F(4400)-p3mup2*F(4459)))
       Fa20(1926) = mup2mup3*p3mup4*F(4098)+mup2mup4*(p3mup3*F(4099)+4
     -   *(p1mup3*F(3965)+p4mup3*F(3968)+p2mup3*F(4409)))-2*mup3mup4*
     -   p3mup2*F(4410)+4*m1sq*(mup2mup3*p4mup4*F(362)+mup3mup4*p1mup
     -   2*F(1608)+mup2mup3*p1mup4*F(4024)+mup2mup3*p2mup4*F(4441)+mu
     -   p3mup4*p4mup2*F(4442))
       Fa20(1927) = 2*(mup2mup4*p2mup3-mup2mup3*p2mup4)*F(4412)+4*mup3
     -   mup4*(p1mup2*F(4411)+p3mup2*F(4413)+p4mup2*F(4414))
       Fa20(1928) = p3mup1*(p1mup3*F(3951)-p3mup3*F(3955)+p2mup3*F(440
     -   1))
       Fa20(1929) = p3mup1*(p1mup3*F(3963)+2*(p3mup3*F(3959)+p2mup3*F(
     -   4405)))
       Fa20(1930) = p3mup1*p3mup4*F(3954)+p3mup1*p4mup4*F(3956)+p3mup4
     -   *p4mup1*F(3966)+4*(p1mup4*p3mup1*F(3952)+m1sq*p4mup1*(p1mup4
     -   *F(346)+p2mup4*F(2318)+p4mup4*F(3969))+p2mup4*p3mup1*F(4402)
     -   )
       Fa20(1931) = p4mup1*(p1mup4*F(3963)+2*(p3mup4*F(3960)+p4mup4*F(
     -   3961)+p2mup4*F(4407)))
       Fa20(1932) = 2*p3mup1*(p1mup4*F(3963)+2*(p3mup4*F(3960)+p4mup4*
     -   F(3961)+p2mup4*F(4407)))
       Fa20(1933) = 2*p2mup4*p4mup1*F(4412)+p2mup1*(p3mup4*F(4095)+p4m
     -   up4*F(4097)+4*(p1mup4*F(4094)+p2mup4*F(4458)))
       Fa20(1934) = p2mup1*(p1mup4*F(3963)+2*(p3mup4*F(3960)+p4mup4*F(
     -   3961)+p2mup4*F(4407)))
       Fa20(1935) = 4*mup1mup4*(p1mup3*F(3951)+p3mup3*F(4096)+p2mup3*F
     -   (4401))+mup1mup3*(p3mup4*F(4095)+p4mup4*F(4097)+4*(p1mup4*F(
     -   4094)+p2mup4*F(4458)))
       Fa20(1936) = -2*(mup1mup4*(p1mup3*F(3963)+2*(p3mup3*F(3959)+p2m
     -   up3*F(4405)))-mup1mup3*(p1mup4*F(3963)+2*(p3mup4*F(3960)+p4m
     -   up4*F(3961)+p2mup4*F(4407))))
       Fa20(1937) = mup1mup2*(p3mup4*F(3954)+p4mup4*F(3956)+4*(p1mup4*
     -   F(3952)+p2mup4*F(4402)))
       Fa20(1938) = mup1mup2*(p1mup4*F(3963)+2*(p3mup4*F(3960)+p4mup4*
     -   F(3961)+p2mup4*F(4407)))
       Fa20(1939) = p2mup1*p3mup4*F(4098)-2*p2mup4*p3mup1*F(4412)+4*m1
     -   sq*p2mup1*(p4mup4*F(362)+p1mup4*F(4024)+p2mup4*F(4441))
       Fa20(1940) = mup1mup3*p3mup4*F(4098)+mup1mup4*(p3mup3*F(4099)+4
     -   *(p1mup3*F(3965)+p2mup3*F(4409)))+4*m1sq*mup1mup3*(p4mup4*F(
     -   362)+p1mup4*F(4024)+p2mup4*F(4441))
       Fa20(1941) = mup1mup2*(p3mup4*F(3966)+4*m1sq*(p1mup4*F(346)+p2m
     -   up4*F(2318)+p4mup4*F(3969)))
       Fa20(1942) = p2mup1*(p1mup3*F(3951)+p3mup3*F(4096)+p2mup3*F(440
     -   1))
       Fa20(1943) = p2mup1*p3mup3*F(4099)+4*(p1mup3*p2mup1*F(3965)+p2m
     -   up1*p2mup3*F(4409))+2*p2mup3*p3mup1*F(4412)
       Fa20(1944) = 2*p2mup1*(p1mup3*F(3963)+2*(p3mup3*F(3959)+p2mup3*
     -   F(4405)))
       Fa20(1945) = p3mup1*(p1mup3*F(258)+p3mup3*F(600)+p2mup3*F(2338)
     -   )
       Fa20(1946) = p3mup1*(p1mup3*F(4018)+2*(p3mup3*F(4019)+p2mup3*F(
     -   4438)))
       Fa20(1947) = p3mup1*(p4mup4*F(4016)+2*(p1mup4*F(4014)+p3mup4*F(
     -   4015)+p2mup4*F(4436)))+p4mup1*(p3mup4*F(4022)+2*(p4mup4*F(36
     -   2)+p1mup4*F(4024)+p2mup4*F(4441)))
       Fa20(1948) = p4mup1*(p1mup4*F(4018)+2*(p3mup4*F(4019)+p4mup4*F(
     -   4020)+p2mup4*F(4438)))
       Fa20(1949) = p3mup1*(p1mup4*F(4018)+2*(p3mup4*F(4019)+p4mup4*F(
     -   4020)+p2mup4*F(4438)))
       Fa20(1950) = p2mup1*(p4mup4*F(4016)+2*(p1mup4*F(4014)+p3mup4*F(
     -   4015)+p2mup4*F(4436)))+p2mup4*p4mup1*F(4444)
       Fa20(1951) = p2mup1*(p1mup4*F(4018)+2*(p3mup4*F(4019)+p4mup4*F(
     -   4020)+p2mup4*F(4438)))
       Fa20(1952) = -(mup1mup3*(p4mup4*F(4016)+2*(p1mup4*F(4014)+p3mup
     -   4*F(4015)+p2mup4*F(4436))))+2*mup1mup4*(p1mup3*F(264)+p3mup3
     -   *F(282)+p2mup3*F(4445))
       Fa20(1953) = -(mup1mup4*(p1mup3*F(4018)+2*(p3mup3*F(4019)+p2mup
     -   3*F(4438))))+mup1mup3*(p1mup4*F(4018)+2*(p3mup4*F(4019)+p4mu
     -   p4*F(4020)+p2mup4*F(4438)))
       Fa20(1954) = mup1mup2*(p4mup4*F(4016)+2*(p1mup4*F(4014)+p3mup4*
     -   F(4015)+p2mup4*F(4436)))
       Fa20(1955) = mup1mup2*(p1mup4*F(4018)+2*(p3mup4*F(4019)+p4mup4*
     -   F(4020)+p2mup4*F(4438)))
       Fa20(1956) = p2mup1*(2*(p1mup4*F(346)+p2mup4*F(2318)+p4mup4*F(3
     -   969))-p3mup4*F(4022))+p2mup4*p3mup1*F(4444)
       Fa20(1957) = mup1mup3*(2*(p1mup4*F(346)+p2mup4*F(2318)+p4mup4*F
     -   (3969))-p3mup4*F(4022))-mup1mup4*(p3mup3*F(4023)+2*(p1mup3*F
     -   (4021)+p2mup3*F(4439)))
       Fa20(1958) = mup1mup2*(p3mup4*F(4022)+2*(p4mup4*F(362)+p1mup4*F
     -   (4024)+p2mup4*F(4441)))
       Fa20(1959) = p2mup1*(p1mup3*F(264)+p3mup3*F(282)+p2mup3*F(4445)
     -   )
       Fa20(1960) = p2mup1*p3mup3*F(4023)+2*(p1mup3*p2mup1*F(4021)+p2m
     -   up1*p2mup3*F(4439))+p2mup3*p3mup1*F(4444)
       Fa20(1961) = p2mup1*(p1mup3*F(4018)+2*(p3mup3*F(4019)+p2mup3*F(
     -   4438)))
       Fa20(1962) = p3mup1*(p1mup3*F(4028)+2*(p3mup3*F(4035)+p2mup3*F(
     -   4448)))
       Fa20(1963) = 2*p1mup4*(p3mup1*F(264)+p4mup1*F(4024))+p4mup1*(p3
     -   mup4*F(4030)+2*(p4mup4*F(362)+p2mup4*F(4441)))+p3mup1*(p4mup
     -   4*F(4038)+2*(p3mup4*F(282)+p2mup4*F(4445)))
       Fa20(1964) = p4mup1*(p1mup4*F(4028)+2*(p3mup4*F(4035)+p4mup4*F(
     -   4041)+p2mup4*F(4448)))
       Fa20(1965) = p3mup1*(p1mup4*F(4029)+2*(p3mup4*F(4034)+p4mup4*F(
     -   4040)+p2mup4*F(4447)))
       Fa20(1966) = p2mup4*p4mup1*F(1056)+p2mup1*(2*(p1mup4*F(258)+p3m
     -   up4*F(600)+p2mup4*F(2338))+p4mup4*F(4039))
       Fa20(1967) = p2mup1*(p1mup4*F(4028)+2*(p3mup4*F(4035)+p4mup4*F(
     -   4041)+p2mup4*F(4448)))
       Fa20(1968) = mup1mup3*(2*(p1mup4*F(258)+p3mup4*F(600)+p2mup4*F(
     -   2338))+p4mup4*F(4039))+2*mup1mup4*(p1mup3*F(264)+p3mup3*F(28
     -   2)+p2mup3*F(4445))
       Fa20(1969) = mup1mup3*(p1mup4*F(4029)+2*(p3mup4*F(4034)+p4mup4*
     -   F(4040)+p2mup4*F(4447)))+mup1mup4*(p1mup3*F(4028)+2*(p3mup3*
     -   F(4035)+p2mup3*F(4448)))
       Fa20(1970) = mup1mup2*(p4mup4*F(4038)+2*(p1mup4*F(264)+p3mup4*F
     -   (282)+p2mup4*F(4445)))
       Fa20(1971) = mup1mup2*(p1mup4*F(4028)+2*(p3mup4*F(4035)+p4mup4*
     -   F(4041)+p2mup4*F(4448)))
       Fa20(1972) = p2mup4*p3mup1*F(1057)+p2mup1*(2*(p1mup4*F(346)+p2m
     -   up4*F(2318)+p4mup4*F(3969))+p3mup4*F(4031))
       Fa20(1973) = mup1mup3*(2*(p1mup4*F(346)+p2mup4*F(2318)+p4mup4*F
     -   (3969))+p3mup4*F(4031))+mup1mup4*(p3mup3*F(4030)+2*(p1mup3*F
     -   (4024)+p2mup3*F(4441)))
       Fa20(1974) = mup1mup2*(p3mup4*F(4030)+2*(p4mup4*F(362)+p1mup4*F
     -   (4024)+p2mup4*F(4441)))
       Fa20(1975) = p2mup3*p3mup1*F(1056)+p2mup1*p3mup3*F(4030)+2*(p1m
     -   up3*p2mup1*F(4024)+p2mup1*p2mup3*F(4441))
       Fa20(1976) = p2mup1*(p1mup3*F(4029)+2*(p3mup3*F(4034)+p2mup3*F(
     -   4447)))
       Fa20(1977) = mup1mup4*p2mup3*F(1056)+mup1mup3*p2mup4*F(1057)
       Fa20(1978) = p3mup1*(p1mup3*F(264)+p3mup3*F(282)+p2mup3*F(4445)
     -   )
       Fa20(1979) = p3mup1*(p1mup3*F(4029)+2*(p3mup3*F(4034)+p2mup3*F(
     -   4447)))
       Fa20(1980) = 2*p1mup4*(p3mup1*F(258)+p4mup1*F(346))+p4mup1*(2*(
     -   p2mup4*F(2318)+p4mup4*F(3969))+p3mup4*F(4031))+p3mup1*(2*(p3
     -   mup4*F(600)+p2mup4*F(2338))+p4mup4*F(4039))
       Return
       End
