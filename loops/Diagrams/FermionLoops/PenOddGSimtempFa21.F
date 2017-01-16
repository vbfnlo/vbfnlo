c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenOddGSim21(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4
     -   ,p3mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5m
     -   up3,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mu
     -   p3,mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa21)
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
       Complex*16 Fa21(1981:2079)
       COMMON/FaPenOddGSim/Fa
       Fa21(1981) = p4mup1*(p1mup4*F(4029)+2*(p3mup4*F(4034)+p4mup4*F(
     -   4040)+p2mup4*F(4447)))
       Fa21(1982) = p3mup1*(p1mup4*F(4028)+2*(p3mup4*F(4035)+p4mup4*F(
     -   4041)+p2mup4*F(4448)))
       Fa21(1983) = p2mup4*p4mup1*F(1057)+p2mup1*(p4mup4*F(4038)+2*(p1
     -   mup4*F(264)+p3mup4*F(282)+p2mup4*F(4445)))
       Fa21(1984) = p2mup1*(p1mup4*F(4029)+2*(p3mup4*F(4034)+p4mup4*F(
     -   4040)+p2mup4*F(4447)))
       Fa21(1985) = 2*mup1mup4*(p1mup3*F(258)+p3mup3*F(600)+p2mup3*F(2
     -   338))+mup1mup3*(p4mup4*F(4038)+2*(p1mup4*F(264)+p3mup4*F(282
     -   )+p2mup4*F(4445)))
       Fa21(1986) = mup1mup4*(p1mup3*F(4029)+2*(p3mup3*F(4034)+p2mup3*
     -   F(4447)))+mup1mup3*(p1mup4*F(4028)+2*(p3mup4*F(4035)+p4mup4*
     -   F(4041)+p2mup4*F(4448)))
       Fa21(1987) = mup1mup2*(2*(p1mup4*F(258)+p3mup4*F(600)+p2mup4*F(
     -   2338))+p4mup4*F(4039))
       Fa21(1988) = mup1mup2*(p1mup4*F(4029)+2*(p3mup4*F(4034)+p4mup4*
     -   F(4040)+p2mup4*F(4447)))
       Fa21(1989) = p2mup4*p3mup1*F(1056)+p2mup1*(p3mup4*F(4030)+2*(p4
     -   mup4*F(362)+p1mup4*F(4024)+p2mup4*F(4441)))
       Fa21(1990) = mup1mup4*(2*(p1mup3*F(346)+p2mup3*F(2318))+p3mup3*
     -   F(4031))+mup1mup3*(p3mup4*F(4030)+2*(p4mup4*F(362)+p1mup4*F(
     -   4024)+p2mup4*F(4441)))
       Fa21(1991) = mup1mup2*(2*(p1mup4*F(346)+p2mup4*F(2318)+p4mup4*F
     -   (3969))+p3mup4*F(4031))
       Fa21(1992) = p2mup1*(p1mup3*F(258)+p3mup3*F(600)+p2mup3*F(2338)
     -   )
       Fa21(1993) = p2mup3*p3mup1*F(1057)+2*(p1mup3*p2mup1*F(346)+p2mu
     -   p1*p2mup3*F(2318))+p2mup1*p3mup3*F(4031)
       Fa21(1994) = p2mup1*(p1mup3*F(4028)+2*(p3mup3*F(4035)+p2mup3*F(
     -   4448)))
       Fa21(1995) = mup1mup3*p2mup4*F(1056)+mup1mup4*p2mup3*F(1057)
       Fa21(1996) = p3mup1*p4mup4*F(4016)+p3mup4*p4mup1*F(4022)-2*(p1m
     -   up4*p4mup1*F(346)+p2mup4*p4mup1*F(2318)+p4mup1*p4mup4*F(3969
     -   )-p1mup4*p3mup1*F(4014)-p3mup1*p3mup4*F(4015)-p2mup4*p3mup1*
     -   F(4436))
       Fa21(1997) = 2*mup1mup4*(p1mup3*F(258)+p3mup3*F(600)+p2mup3*F(2
     -   338))+mup1mup3*(p4mup4*F(4016)+2*(p1mup4*F(4014)+p3mup4*F(40
     -   15)+p2mup4*F(4436)))
       Fa21(1998) = p2mup1*(p3mup4*F(4022)+2*(p4mup4*F(362)+p1mup4*F(4
     -   024)+p2mup4*F(4441)))-p2mup4*p3mup1*F(4444)
       Fa21(1999) = mup1mup4*(p3mup3*F(4023)+2*(p1mup3*F(4021)+p2mup3*
     -   F(4439)))+mup1mup3*(p3mup4*F(4022)+2*(p4mup4*F(362)+p1mup4*F
     -   (4024)+p2mup4*F(4441)))
       Fa21(2000) = mup1mup2*(2*(p1mup4*F(346)+p2mup4*F(2318)+p4mup4*F
     -   (3969))-p3mup4*F(4022))
       Fa21(2001) = p3mup1*(p1mup3*F(3951)+p3mup3*F(4096)+p2mup3*F(440
     -   1))
       Fa21(2002) = 2*p3mup1*(p1mup3*F(3963)+2*(p3mup3*F(3959)+p2mup3*
     -   F(4405)))
       Fa21(2003) = p3mup1*p3mup4*F(4095)+p3mup1*p4mup4*F(4097)+p3mup4
     -   *p4mup1*F(4098)+4*(p1mup4*p3mup1*F(4094)+m1sq*p4mup1*(p4mup4
     -   *F(362)+p1mup4*F(4024)+p2mup4*F(4441))+p2mup4*p3mup1*F(4458)
     -   )
       Fa21(2004) = 2*p4mup1*(p1mup4*F(3963)+2*(p3mup4*F(3960)+p4mup4*
     -   F(3961)+p2mup4*F(4407)))
       Fa21(2005) = p3mup1*(p1mup4*F(3963)+2*(p3mup4*F(3960)+p4mup4*F(
     -   3961)+p2mup4*F(4407)))
       Fa21(2006) = p2mup1*(p3mup4*F(3954)+p4mup4*F(3956)+4*(p1mup4*F(
     -   3952)+p2mup4*F(4402)))-2*p2mup4*p4mup1*F(4412)
       Fa21(2007) = 2*p2mup1*(p1mup4*F(3963)+2*(p3mup4*F(3960)+p4mup4*
     -   F(3961)+p2mup4*F(4407)))
       Fa21(2008) = -4*mup1mup4*(p1mup3*F(3951)-p3mup3*F(3955)+p2mup3*
     -   F(4401))+mup1mup3*(p3mup4*F(3954)+p4mup4*F(3956)+4*(p1mup4*F
     -   (3952)+p2mup4*F(4402)))
       Fa21(2009) = 2*(mup1mup4*(p1mup3*F(3963)+2*(p3mup3*F(3959)+p2mu
     -   p3*F(4405)))-mup1mup3*(p1mup4*F(3963)+2*(p3mup4*F(3960)+p4mu
     -   p4*F(3961)+p2mup4*F(4407))))
       Fa21(2010) = mup1mup2*(p3mup4*F(4095)+p4mup4*F(4097)+4*(p1mup4*
     -   F(4094)+p2mup4*F(4458)))
       Fa21(2011) = 2*mup1mup2*(p1mup4*F(3963)+2*(p3mup4*F(3960)+p4mup
     -   4*F(3961)+p2mup4*F(4407)))
       Fa21(2012) = p2mup1*p3mup4*F(3966)+4*m1sq*p2mup1*(p1mup4*F(346)
     -   +p2mup4*F(2318)+p4mup4*F(3969))+2*p2mup4*p3mup1*F(4412)
       Fa21(2013) = mup1mup3*p3mup4*F(3966)+4*m1sq*mup1mup3*(p1mup4*F(
     -   346)+p2mup4*F(2318)+p4mup4*F(3969))+mup1mup4*(p3mup3*F(3967)
     -   -4*(p1mup3*F(3965)+p2mup3*F(4409)))
       Fa21(2014) = mup1mup2*(p3mup4*F(4098)+4*m1sq*(p4mup4*F(362)+p1m
     -   up4*F(4024)+p2mup4*F(4441)))
       Fa21(2015) = p2mup1*(p1mup3*F(3951)-p3mup3*F(3955)+p2mup3*F(440
     -   1))
       Fa21(2016) = p2mup1*p3mup3*F(3967)-4*(p1mup3*p2mup1*F(3965)+p2m
     -   up1*p2mup3*F(4409))-2*p2mup3*p3mup1*F(4412)
       Fa21(2017) = p2mup1*(p1mup3*F(3963)+2*(p3mup3*F(3959)+p2mup3*F(
     -   4405)))
       Fa21(2018) = p3mup1*(p1mup2*F(4400)+p3mup2*F(4403))
       Fa21(2019) = 2*p3mup1*(p1mup2*F(3963)+2*p3mup2*F(4406))
       Fa21(2020) = p4mup1*(p3mup3*F(3967)-4*(p1mup3*F(3965)+p4mup3*F(
     -   3968)+p2mup3*F(4409)))
       Fa21(2021) = 2*p4mup1*(p1mup3*F(3963)+2*(p3mup3*F(3959)+p4mup3*
     -   F(3962)+p2mup3*F(4405)))
       Fa21(2022) = mup1mup4*(p1mup2*F(4400)-p3mup2*F(4459))
       Fa21(2023) = 2*mup1mup4*(p1mup2*F(3963)+2*p3mup2*F(4406))
       Fa21(2024) = 4*mup1mup4*p4mup3*F(3968)-2*mup3mup4*p3mup1*F(4412
     -   )
       Fa21(2025) = 4*m1sq*mup1mup4*p1mup2*F(1608)-2*(mup1mup4*p3mup2*
     -   F(4410)-mup2mup4*p3mup1*F(4412))
       Fa21(2026) = p2mup1*p4mup3*F(3957)+p2mup3*p4mup1*F(4412)
       Fa21(2027) = p2mup1*(p1mup2*F(4400)-p3mup2*F(4459))
       Fa21(2028) = 2*mup1mup2*(p4mup3*F(3957)+2*(p3mup3*F(3955)-p2mup
     -   3*F(4401)))+4*mup1mup3*(p1mup2*F(4400)+p3mup2*F(4403))
       Fa21(2029) = -2*p2mup1*p3mup2*F(4410)+4*(m1sq*p1mup2*p2mup1*F(1
     -   608)+p3mup1*(p1mup2*F(4411)+p3mup2*F(4413)))
       Fa21(2030) = p2mup1*(p1mup2*F(3963)+2*p3mup2*F(4406))
       Fa21(2031) = 4*m1sq*mup1mup3*p1mup2*F(1262)+mup1mup2*(p3mup3*F(
     -   3967)-4*(p4mup3*F(3968)+p2mup3*F(4409)))+2*(mup1mup3*p3mup2*
     -   F(4410)-mup2mup3*p3mup1*F(4412))
       Fa21(2032) = 4*mup1mup2*(p3mup3*F(3959)+p4mup3*F(3962)+p2mup3*F
     -   (4405))-2*mup1mup3*(p1mup2*F(3963)+2*p3mup2*F(4406))
       Fa21(2033) = (-(mup2mup3*p2mup1)+mup1mup2*p2mup3)*F(4412)+2*mup
     -   1mup3*(p1mup2*F(4411)+p3mup2*F(4413))
       Fa21(2034) = -2*mup2mup4*p2mup1*F(4412)+4*mup1mup4*(p1mup2*F(44
     -   11)+p3mup2*F(4413))
       Fa21(2035) = p3mup1*(p1mup2*F(4018)+2*p3mup2*F(4019))
       Fa21(2036) = p4mup1*(p3mup3*F(4023)+2*(p1mup3*F(4021)+p4mup3*F(
     -   4025)+p2mup3*F(4439)))
       Fa21(2037) = p4mup1*(p1mup3*F(4018)+2*(p3mup3*F(4019)+p4mup3*F(
     -   4020)+p2mup3*F(4438)))
       Fa21(2038) = mup1mup4*(p1mup2*F(4018)+2*p3mup2*F(4019))
       Fa21(2039) = -2*mup1mup4*p4mup3*F(4025)+mup3mup4*p3mup1*F(4444)
       Fa21(2040) = 2*mup1mup4*p1mup2*F(1262)+mup1mup4*p3mup2*F(4440)-
     -   mup2mup4*p3mup1*F(4444)
       Fa21(2041) = p2mup1*p4mup3*F(4017)+p2mup3*p4mup1*F(4444)
       Fa21(2042) = 2*(mup1mup2*p3mup3*F(600)+mup1mup3*p1mup2*F(1608)+
     -   mup1mup2*p2mup3*F(2338)+mup1mup3*p3mup2*F(3506))-mup1mup2*p4
     -   mup3*F(4017)
       Fa21(2043) = p2mup1*p3mup2*F(4440)+2*(p1mup2*p2mup1*F(1262)-p3m
     -   up1*p3mup2*F(4019)-p1mup2*p3mup1*F(4443))
       Fa21(2044) = p2mup1*(p1mup2*F(4018)+2*p3mup2*F(4019))
       Fa21(2045) = mup1mup2*p3mup3*F(4023)+2*(mup1mup3*p1mup2*F(1608)
     -   +mup1mup2*p4mup3*F(4025)+mup1mup2*p2mup3*F(4439))-mup1mup3*p
     -   3mup2*F(4440)+mup2mup3*p3mup1*F(4444)
       Fa21(2046) = mup1mup3*(p1mup2*F(4018)+2*p3mup2*F(4019))-2*mup1m
     -   up2*(p3mup3*F(4019)+p4mup3*F(4020)+p2mup3*F(4438))
       Fa21(2047) = 2*mup1mup3*(p3mup2*F(4019)+p1mup2*F(4443))+(-(mup2
     -   mup3*p2mup1)+mup1mup2*p2mup3)*F(4444)
       Fa21(2048) = -2*mup1mup4*(p3mup2*F(4019)+p1mup2*F(4443))+mup2mu
     -   p4*p2mup1*F(4444)
       Fa21(2049) = p3mup1*(p1mup2*F(4029)+2*p3mup2*F(4034))
       Fa21(2050) = p4mup1*(2*(p1mup3*F(346)+p2mup3*F(2318)+p4mup3*F(3
     -   969))+p3mup3*F(4031))
       Fa21(2051) = p4mup1*(p1mup3*F(4029)+2*(p3mup3*F(4034)+p4mup3*F(
     -   4040)+p2mup3*F(4447)))
       Fa21(2052) = mup1mup4*(p1mup2*F(4029)+2*p3mup2*F(4034))
       Fa21(2053) = 2*mup1mup4*p4mup3*F(362)+mup3mup4*p3mup1*F(1057)
       Fa21(2054) = mup2mup4*p3mup1*F(1056)+2*mup1mup4*p1mup2*F(1262)+
     -   mup1mup4*p3mup2*F(4452)
       Fa21(2055) = p2mup3*p4mup1*F(1057)+p2mup1*p4mup3*F(4038)
       Fa21(2056) = 2*(mup1mup2*p3mup3*F(600)+mup1mup3*p1mup2*F(1608)+
     -   mup1mup2*p2mup3*F(2338)+mup1mup3*p3mup2*F(3506))+mup1mup2*p4
     -   mup3*F(4039)
       Fa21(2057) = p2mup1*p3mup2*F(4452)+2*(p1mup2*p2mup1*F(1262)+p3m
     -   up1*p3mup2*F(4034)+p1mup2*p3mup1*F(4453))
       Fa21(2058) = p2mup1*(p1mup2*F(4028)+2*p3mup2*F(4035))
       Fa21(2059) = mup2mup3*p3mup1*F(1057)+2*(mup1mup3*p1mup2*F(1608)
     -   +mup1mup2*p2mup3*F(2318)+mup1mup2*p4mup3*F(3969))+mup1mup2*p
     -   3mup3*F(4031)+mup1mup3*p3mup2*F(4449)
       Fa21(2060) = mup1mup3*(p1mup2*F(4028)+2*p3mup2*F(4035))+2*mup1m
     -   up2*(p3mup3*F(4034)+p4mup3*F(4040)+p2mup3*F(4447))
       Fa21(2061) = mup2mup3*p2mup1*F(1056)+mup1mup2*p2mup3*F(1057)+2*
     -   mup1mup3*(p3mup2*F(4035)+p1mup2*F(4450))
       Fa21(2062) = mup2mup4*p2mup1*F(1057)+2*mup1mup4*(p3mup2*F(4034)
     -   +p1mup2*F(4453))
       Fa21(2063) = p3mup1*(p1mup2*F(4028)+2*p3mup2*F(4035))
       Fa21(2064) = p4mup1*(p3mup3*F(4030)+2*(p4mup3*F(362)+p1mup3*F(4
     -   024)+p2mup3*F(4441)))
       Fa21(2065) = p4mup1*(p1mup3*F(4028)+2*(p3mup3*F(4035)+p4mup3*F(
     -   4041)+p2mup3*F(4448)))
       Fa21(2066) = mup1mup4*(p1mup2*F(4028)+2*p3mup2*F(4035))
       Fa21(2067) = mup3mup4*p3mup1*F(1056)+2*mup1mup4*p4mup3*F(3969)
       Fa21(2068) = mup2mup4*p3mup1*F(1057)+2*mup1mup4*p1mup2*F(1608)+
     -   mup1mup4*p3mup2*F(4449)
       Fa21(2069) = p2mup3*p4mup1*F(1056)+p2mup1*p4mup3*F(4039)
       Fa21(2070) = 2*mup1mup3*(p1mup2*F(1262)+p3mup2*F(2368))+mup1mup
     -   2*(p4mup3*F(4038)+2*(p3mup3*F(282)+p2mup3*F(4445)))
       Fa21(2071) = p2mup1*p3mup2*F(4449)+2*(p1mup2*p2mup1*F(1608)+p3m
     -   up1*p3mup2*F(4035)+p1mup2*p3mup1*F(4450))
       Fa21(2072) = p2mup1*(p1mup2*F(4029)+2*p3mup2*F(4034))
       Fa21(2073) = mup2mup3*p3mup1*F(1056)+2*mup1mup3*p1mup2*F(1262)+
     -   mup1mup2*(p3mup3*F(4030)+2*(p4mup3*F(362)+p2mup3*F(4441)))+m
     -   up1mup3*p3mup2*F(4452)
       Fa21(2074) = mup1mup3*(p1mup2*F(4029)+2*p3mup2*F(4034))+2*mup1m
     -   up2*(p3mup3*F(4035)+p4mup3*F(4041)+p2mup3*F(4448))
       Fa21(2075) = mup1mup2*p2mup3*F(1056)+mup2mup3*p2mup1*F(1057)+2*
     -   mup1mup3*(p3mup2*F(4034)+p1mup2*F(4453))
       Fa21(2076) = mup2mup4*p2mup1*F(1056)+2*mup1mup4*(p3mup2*F(4035)
     -   +p1mup2*F(4450))
       Fa21(2077) = 2*mup1mup4*p4mup3*F(4025)-mup3mup4*p3mup1*F(4444)
       Fa21(2078) = 2*mup1mup4*p1mup2*F(1608)-mup1mup4*p3mup2*F(4440)+
     -   mup2mup4*p3mup1*F(4444)
       Fa21(2079) = 2*mup1mup3*(p1mup2*F(1262)+p3mup2*F(2368))+mup1mup
     -   2*(p4mup3*F(4017)+2*(p3mup3*F(282)+p2mup3*F(4445)))
       Return
       End
