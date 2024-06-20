c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenOddG6(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1
     -   ,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3m
     -   up5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,
     -   p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,m
     -   up2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa6)
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
       Complex*16 Fa6(1461:1752)
       COMMON/FaPenOddG/Fa
       Fa6(1461) = mup4mup5*F(4575)+4*(p1mup5*p3mup4*F(4577)+p2mup5*p3
     -   mup4*F(4579)+p1mup4*(p2mup5*F(4561)+p3mup5*F(4576)+p1mup5*F(
     -   4581)+p4mup5*F(4582))+p3mup4*p4mup5*F(4583)+p3mup4*p3mup5*F(
     -   4584)+p1mup5*p4mup4*F(4586)+p3mup5*p4mup4*F(4588)+2*(p4mup4*
     -   p4mup5*F(4587)+p2mup5*p4mup4*F(4589))+p2mup4*(p3mup5*F(4578)
     -   +p1mup5*F(4580)+2*(p4mup5*F(4585)+p2mup5*F(4590))))
       Fa6(1462) = p2mup5*F(4591)+p3mup5*F(4592)+p4mup5*F(4593)+2*p1mu
     -   p5*F(4594)
       Fa6(1463) = p1mup5*F(4595)+p3mup5*F(4596)+p4mup5*F(4597)+2*p2mu
     -   p5*F(4598)
       Fa6(1464) = p1mup5*F(4599)+p2mup5*F(4600)+p4mup5*F(4601)+2*p3mu
     -   p5*F(4602)
       Fa6(1465) = p1mup5*F(4603)+p2mup5*F(4604)+p3mup5*F(4605)+2*p4mu
     -   p5*F(4606)
       Fa6(1466) = p2mup1*F(4607)+p3mup1*F(4608)+p4mup1*F(4609)
       Fa6(1467) = p1mup1*F(4526)+p3mup1*F(4610)+p4mup1*F(4611)
       Fa6(1468) = p2mup1*F(4532)+p1mup1*F(4612)+p4mup1*F(4613)
       Fa6(1469) = p1mup1*F(4614)+p2mup1*F(4615)+p3mup1*F(4616)
       Fa6(1470) = p1mup4*F(4618)+p2mup4*F(4619)+p3mup4*F(4620)+p4mup4
     -   *F(4621)
       Fa6(1471) = mup4mup5*F(4622)-4*(p1mup5*p4mup4*F(4627)-p2mup4*p4
     -   mup5*F(4628)-p2mup5*p4mup4*F(4629)+p3mup4*p4mup5*F(4633)+p3m
     -   up5*p4mup4*F(4634)+p1mup4*(p4mup5*F(4626)+2*(p1mup5*F(4623)-
     -   p2mup5*F(4630)+p3mup5*F(4635)))+2*(p1mup5*p2mup4*F(4624)+p1m
     -   up5*p3mup4*F(4625)-p2mup4*p2mup5*F(4631)-p2mup5*p3mup4*F(463
     -   2)+p2mup4*p3mup5*F(4636)+p3mup4*p3mup5*F(4637))+p4mup4*p4mup
     -   5*F(4638))
       Fa6(1472) = mup4mup5*F(4639)-4*(p1mup5*(p2mup4*F(4641)+p3mup4*F
     -   (4642))-p2mup5*p3mup4*F(4644)-p1mup4*(p3mup5*F(4634)+p2mup5*
     -   F(4640)+p1mup5*F(4645)+p4mup5*F(4646))+p2mup4*(p3mup5*F(4643
     -   )+p4mup5*F(4647)+p2mup5*F(4648))-p1mup5*p4mup4*F(4650)+p2mup
     -   5*p4mup4*F(4652)-2*(p3mup4*p4mup5*F(4649)+p4mup4*p4mup5*F(46
     -   51)+p3mup5*p4mup4*F(4653)+p3mup4*p3mup5*F(4654)))
       Fa6(1473) = p1mup4*F(4655)+p3mup4*F(4656)+p2mup4*F(4657)+p4mup4
     -   *F(4658)
       Fa6(1474) = p1mup4*F(4659)+p2mup4*F(4660)+p3mup4*F(4661)+2*p4mu
     -   p4*F(4662)
       Fa6(1475) = mup4mup5*F(4609)+4*(-(p1mup5*p2mup4*F(4528))+p1mup4
     -   *(p2mup5*F(4663)+p3mup5*F(4664))-p1mup5*p3mup4*F(4665)-p2mup
     -   4*p3mup5*F(4666)+p2mup5*p3mup4*F(4667)-p2mup4*p4mup5*F(4668)
     -   -p3mup4*p4mup5*F(4669)-p2mup4*p2mup5*F(4670)-p3mup4*p3mup5*F
     -   (4671)-p2mup5*p4mup4*F(4672)-p3mup5*p4mup4*F(4673))
       Fa6(1476) = p2mup3*F(4674)-2*p1mup3*F(4675)+p3mup3*F(4676)+p4mu
     -   p3*F(4677)
       Fa6(1477) = p1mup3*F(4678)+2*p2mup3*F(4679)+p3mup3*F(4680)+p4mu
     -   p3*F(4681)
       Fa6(1478) = p1mup3*F(4682)+p2mup3*F(4683)+p4mup3*F(4684)+2*p3mu
     -   p3*F(4685)
       Fa6(1479) = mup3mup5*F(4612)+4*(p1mup5*p4mup3*F(4665)-p2mup5*p4
     -   mup3*F(4688)+p3mup5*p4mup3*F(4689)+p4mup5*(p1mup3*F(4686)-p2
     -   mup3*F(4687)+p4mup3*F(4690)-p3mup3*F(4691)))
       Fa6(1480) = p1mup3*F(4692)+p2mup3*F(4693)+p3mup3*F(4694)+2*p4mu
     -   p3*F(4695)
       Fa6(1481) = mup3mup5*F(4696)+p2mup3*(p2mup5*F(4701)-4*(p4mup5*F
     -   (4548)+p1mup5*F(4698)-p3mup5*F(4702)))-4*(p1mup3*p3mup5*F(46
     -   99)+p1mup5*p3mup3*F(4700)+p2mup5*(p4mup3*F(4552)+p1mup3*F(46
     -   97)-p3mup3*F(4703))+p3mup3*p3mup5*F(4704)+p3mup3*p4mup5*F(47
     -   05)+p3mup5*p4mup3*F(4708)+2*(p1mup3*p1mup5*F(4623)+p1mup3*p4
     -   mup5*F(4706)+p4mup3*p4mup5*F(4707)+p1mup5*p4mup3*F(4709)))
       Fa6(1482) = mup3mup5*F(4710)+4*(-(p1mup3*p4mup5*F(4705))-p1mup5
     -   *p4mup3*F(4708)-p1mup3*p1mup5*F(4711)-p1mup3*p3mup5*F(4713)+
     -   p2mup3*(p4mup5*F(4549)+p2mup5*F(4714)+p1mup5*F(4715)+p3mup5*
     -   F(4716))-2*(p4mup3*p4mup5*F(4649)+p3mup5*p4mup3*F(4654)+p3mu
     -   p3*p4mup5*F(4654)+p3mup3*p3mup5*F(4717))-p1mup5*p3mup3*F(471
     -   8)+p2mup5*(p4mup3*F(4553)-p1mup3*F(4712)+p3mup3*F(4719)))
       Fa6(1483) = p2mup4*F(4720)+p3mup4*F(4721)-2*p1mup4*F(4722)+p4mu
     -   p4*F(4723)
       Fa6(1484) = 2*p1mup5*F(4665)-p2mup5*F(4732)+p3mup5*F(4733)+p4mu
     -   p5*F(4734)
       Fa6(1485) = p1mup5*F(4735)+p3mup5*F(4736)+2*p2mup5*F(4737)+p4mu
     -   p5*F(4738)
       Fa6(1486) = p1mup5*F(4739)-p2mup5*F(4740)+2*(p4mup5*F(4671)+p3m
     -   up5*F(4741))
       Fa6(1487) = p1mup5*F(4734)-p2mup5*F(4742)+2*(p4mup5*F(4743)+p3m
     -   up5*F(4744))
       Fa6(1488) = -(p2mup1*F(4730))+p3mup1*F(4745)+p4mup1*F(4746)
       Fa6(1489) = p1mup1*F(4730)+p3mup1*F(4731)+p4mup1*F(7439)
       Fa6(1490) = -(p2mup1*F(4731))+p1mup1*F(4747)
       Fa6(1491) = p1mup1*F(4748)-p2mup1*F(7439)
       Fa6(1492) = p2mup4*F(4740)+p1mup4*F(4750)+2*(p3mup4*F(4751)+p4m
     -   up4*F(4752))
       Fa6(1493) = p1mup4*F(4735)+p3mup4*F(4736)-2*p2mup4*F(4753)-p4mu
     -   p4*F(4754)
       Fa6(1494) = p2mup4*F(4742)+p1mup4*F(4755)+2*(p4mup4*F(4756)+p3m
     -   up4*F(4757))
       Fa6(1495) = p2mup3*F(4758)+2*p1mup3*F(4759)+p3mup3*F(4760)+p4mu
     -   p3*F(4761)
       Fa6(1496) = p1mup3*F(4762)-2*p2mup3*F(4763)-p3mup3*F(4764)+p4mu
     -   p3*F(4765)
       Fa6(1497) = p1mup3*F(4739)-p2mup3*F(4740)+2*(p4mup3*F(4671)+p3m
     -   up3*F(4741))
       Fa6(1498) = p1mup3*F(4734)-p2mup3*F(4742)+2*(p4mup3*F(4743)+p3m
     -   up3*F(4744))
       Fa6(1499) = p2mup4*F(4732)+p4mup4*F(4755)+p3mup4*F(4766)+2*p1mu
     -   p4*F(4767)
       Fa6(1500) = 2*p1mup5*F(4665)+p4mup5*F(4734)+p2mup5*F(4768)+p3mu
     -   p5*F(4769)
       Fa6(1501) = p4mup5*F(4742)-p1mup5*F(4768)+2*p2mup5*F(4770)+p3mu
     -   p5*F(4771)
       Fa6(1502) = p1mup5*F(4769)-p2mup5*F(4771)+2*(p4mup5*F(4744)+p3m
     -   up5*F(4772))
       Fa6(1503) = p4mup1*F(4746)+p3mup1*F(4773)-p2mup1*F(7437)
       Fa6(1504) = p1mup1*F(7437)+p3mup1*F(7438)+p4mup1*F(7439)
       Fa6(1505) = p1mup1*F(4774)-p2mup1*F(7438)
       Fa6(1506) = p2mup4*F(4771)+p1mup4*F(4776)+2*(p4mup4*F(4757)+p3m
     -   up4*F(4777))
       Fa6(1507) = p4mup4*F(4742)+2*p2mup4*F(4770)+p3mup4*F(4771)+p1mu
     -   p4*F(4778)
       Fa6(1508) = 2*p1mup3*F(4665)+p4mup3*F(4734)+p2mup3*F(4768)+p3mu
     -   p3*F(4769)
       Fa6(1509) = p4mup3*F(4742)-p1mup3*F(4768)+2*p2mup3*F(4770)+p3mu
     -   p3*F(4771)
       Fa6(1510) = p1mup3*F(4769)-p2mup3*F(4771)+2*(p4mup3*F(4744)+p3m
     -   up3*F(4772))
       Fa6(1511) = p4mup4*F(4755)+2*p1mup4*F(4767)+p3mup4*F(4776)+p2mu
     -   p4*F(4778)
       Fa6(1512) = p2mup5*F(4732)+p3mup5*F(4766)+2*p1mup5*F(4767)+p4mu
     -   p5*F(4780)
       Fa6(1513) = p1mup5*F(4735)+p3mup5*F(4736)-2*p2mup5*F(4753)+p4mu
     -   p5*F(4781)
       Fa6(1514) = p2mup5*F(4740)+p1mup5*F(4750)+2*(p3mup5*F(4751)+p4m
     -   up5*F(4752))
       Fa6(1515) = p1mup5*F(4782)+p2mup5*F(4783)+2*(p3mup5*F(4784)+p4m
     -   up5*F(4785))
       Fa6(1516) = p2mup1*F(4730)+p3mup1*F(4747)+p4mup1*F(4786)
       Fa6(1517) = p1mup1*F(4730)+p3mup1*F(4731)+p4mup1*F(4779)
       Fa6(1518) = p2mup1*F(4731)+p1mup1*F(4745)
       Fa6(1519) = p2mup1*F(4779)+p1mup1*F(4787)
       Fa6(1520) = p1mup4*F(4739)-p2mup4*F(4740)+2*(p4mup4*F(4671)+p3m
     -   up4*F(4741))
       Fa6(1521) = -(p1mup4*F(4762))+2*p2mup4*F(4763)+p3mup4*F(4764)+p
     -   4mup4*F(4789)
       Fa6(1522) = -(p2mup4*F(4783))+p1mup4*F(4790)+2*(p3mup4*F(4791)+
     -   p4mup4*F(4792))
       Fa6(1523) = -(p2mup3*F(4758))+2*p1mup3*F(4793)+p4mup3*F(4794)+p
     -   3mup3*F(4795)
       Fa6(1524) = -(p1mup3*F(4762))+2*p2mup3*F(4763)+p3mup3*F(4764)+p
     -   4mup3*F(4789)
       Fa6(1525) = p2mup3*F(4740)+p1mup3*F(4750)+2*(p3mup3*F(4751)+p4m
     -   up3*F(4752))
       Fa6(1526) = p1mup3*F(4782)+p2mup3*F(4783)+2*(p3mup3*F(4784)+p4m
     -   up3*F(4785))
       Fa6(1527) = p2mup4*F(4758)+2*p1mup4*F(4759)+p3mup4*F(4760)+p4mu
     -   p4*F(4796)
       Fa6(1528) = 2*p1mup5*F(4665)-p2mup5*F(4732)+p3mup5*F(4733)+p4mu
     -   p5*F(4797)
       Fa6(1529) = p1mup5*F(4735)+p3mup5*F(4736)+2*p2mup5*F(4737)+p4mu
     -   p5*F(4781)
       Fa6(1530) = -(p2mup5*F(4783))+p1mup5*F(4790)+2*(p3mup5*F(4791)+
     -   p4mup5*F(4792))
       Fa6(1531) = -(p2mup1*F(4730))+p3mup1*F(4745)+p4mup1*F(4787)
       Fa6(1532) = -(p2mup1*F(4779))+p1mup1*F(4786)
       Fa6(1533) = p1mup4*F(4782)+p2mup4*F(4783)+2*(p3mup4*F(4784)+p4m
     -   up4*F(4785))
       Fa6(1534) = p2mup3*F(4758)+2*p1mup3*F(4759)+p3mup3*F(4760)+p4mu
     -   p3*F(4796)
       Fa6(1535) = -(p2mup3*F(4783))+p1mup3*F(4790)+2*(p3mup3*F(4791)+
     -   p4mup3*F(4792))
       Fa6(1536) = -(p2mup4*F(4758))+2*p1mup4*F(4793)+p4mup4*F(4794)+p
     -   3mup4*F(4795)
       Fa6(1537) = p4mup5*F(4755)+2*p1mup5*F(4767)+p3mup5*F(4776)+p2mu
     -   p5*F(4778)
       Fa6(1538) = p4mup5*F(4742)+2*p2mup5*F(4770)+p3mup5*F(4771)+p1mu
     -   p5*F(4778)
       Fa6(1539) = p2mup5*F(4771)+p1mup5*F(4776)+2*(p4mup5*F(4757)+p3m
     -   up5*F(4777))
       Fa6(1540) = p2mup5*F(4742)+p1mup5*F(4755)+2*(p4mup5*F(4756)+p3m
     -   up5*F(4757))
       Fa6(1541) = p4mup1*F(4748)+p3mup1*F(4774)+p2mup1*F(7437)
       Fa6(1542) = p1mup1*F(4773)+p2mup1*F(7438)
       Fa6(1543) = p1mup1*F(4746)+p2mup1*F(7439)
       Fa6(1544) = p1mup4*F(4769)-p2mup4*F(4771)+2*(p4mup4*F(4744)+p3m
     -   up4*F(4772))
       Fa6(1545) = p4mup4*F(4742)-p1mup4*F(4768)+2*p2mup4*F(4770)+p3mu
     -   p4*F(4771)
       Fa6(1546) = p1mup4*F(4734)-p2mup4*F(4742)+2*(p4mup4*F(4743)+p3m
     -   up4*F(4744))
       Fa6(1547) = p4mup3*F(4755)+2*p1mup3*F(4767)+p3mup3*F(4776)+p2mu
     -   p3*F(4778)
       Fa6(1548) = p4mup3*F(4742)+2*p2mup3*F(4770)+p3mup3*F(4771)+p1mu
     -   p3*F(4778)
       Fa6(1549) = p2mup3*F(4771)+p1mup3*F(4776)+2*(p4mup3*F(4757)+p3m
     -   up3*F(4777))
       Fa6(1550) = p2mup3*F(4742)+p1mup3*F(4755)+2*(p4mup3*F(4756)+p3m
     -   up3*F(4757))
       Fa6(1551) = 2*p1mup4*F(4665)+p4mup4*F(4734)+p2mup4*F(4768)+p3mu
     -   p4*F(4769)
       Fa6(1552) = p2mup5*F(4732)+p4mup5*F(4755)+p3mup5*F(4766)+2*p1mu
     -   p5*F(4767)
       Fa6(1553) = p1mup5*F(4735)+p3mup5*F(4736)-2*p2mup5*F(4753)-p4mu
     -   p5*F(4754)
       Fa6(1554) = p2mup1*F(4730)+p3mup1*F(4747)+p4mup1*F(4748)
       Fa6(1555) = p1mup4*F(4735)+p3mup4*F(4736)+2*p2mup4*F(4737)+p4mu
     -   p4*F(4738)
       Fa6(1556) = -(p2mup3*F(4758))+2*p1mup3*F(4793)+p3mup3*F(4795)+p
     -   4mup3*F(4801)
       Fa6(1557) = -(p1mup3*F(4762))+2*p2mup3*F(4763)+p3mup3*F(4764)+p
     -   4mup3*F(4802)
       Fa6(1558) = 2*p1mup4*F(4665)-p2mup4*F(4732)+p3mup4*F(4733)+p4mu
     -   p4*F(4734)
       Fa6(1559) = -4*(p1mup3*p4mup5*F(4527)+p4mup3*(p1mup5*F(4528)+(p
     -   2mup5+p3mup5)*F(4529)+p4mup5*F(4530))-(p2mup3+p3mup3)*p4mup5
     -   *F(4531))+mup3mup5*F(4607)
       Fa6(1560) = 4*((p1mup5+p2mup5)*p4mup3*F(4534)+p3mup5*p4mup3*F(4
     -   535)+p4mup5*((p1mup3+p2mup3)*F(4533)+p4mup3*F(4536)-p3mup3*F
     -   (4537)))+mup3mup5*F(4610)
       Fa6(1561) = mup3mup5*F(4803)-4*(p1mup5*p3mup3*F(4544)+p2mup5*p3
     -   mup3*F(4546)+p3mup3*p3mup5*F(4547)+p3mup3*p4mup5*F(4549)-p1m
     -   up5*p4mup3*F(4552)+p3mup5*p4mup3*F(4553)-2*(p4mup3*p4mup5*F(
     -   4551)-p2mup5*p4mup3*F(4554))-p1mup3*(-(p3mup5*F(4543))+p4mup
     -   5*F(4548)+p1mup5*F(4804)+p2mup5*F(4805))+p2mup3*(p3mup5*F(45
     -   45)+2*(p2mup5*F(4540)+p4mup5*F(4550))-p1mup5*F(4806)))
       Fa6(1562) = 4*(2*p1mup5*(p2mup4+p3mup4)*F(4557)+p1mup5*p4mup4*F
     -   (4559)+(p2mup5+p3mup5)*p4mup4*F(4561)+p1mup4*(p4mup5*F(4558)
     -   +2*(p1mup5*F(4556)+(p2mup5+p3mup5)*F(4562)))+(p2mup4+p3mup4)
     -   *(p4mup5*F(4560)+2*(p2mup5+p3mup5)*F(4563))+p4mup4*p4mup5*F(
     -   4564))+mup4mup5*F(4807)
       Fa6(1563) = -4*(p2mup4*p4mup5*F(4567)+p2mup5*p4mup4*F(4568)+p1m
     -   up5*(p4mup4*F(4568)+2*(p2mup4*F(4566)+p3mup4*F(4569)))+p3mup
     -   4*p4mup5*F(4570)+p3mup5*p4mup4*F(4571)+p1mup4*(p4mup5*F(4567
     -   )+2*((p1mup5+p2mup5)*F(4566)+p3mup5*F(4572)))+2*(p2mup4*p2mu
     -   p5*F(4566)+p2mup5*p3mup4*F(4569)+p2mup4*p3mup5*F(4572)+p3mup
     -   4*p3mup5*F(4573))+p4mup4*p4mup5*F(4574))+mup4mup5*F(4808)
       Fa6(1564) = mup4mup5*F(4809)-4*(p2mup5*p3mup4*F(4579)+p3mup4*p4
     -   mup5*F(4583)+p3mup4*p3mup5*F(4584)+p3mup5*p4mup4*F(4588)+2*(
     -   p2mup4*p4mup5*F(4585)+p4mup4*p4mup5*F(4587)+p2mup5*p4mup4*F(
     -   4589)+p2mup4*p2mup5*F(4590))+p1mup4*(p2mup5*F(4561)+p1mup5*F
     -   (4581)+p4mup5*F(4582)-p3mup5*F(4810))-p2mup4*p3mup5*F(4811)+
     -   p1mup5*(p3mup4*F(4577)+p4mup4*F(4586)-p2mup4*F(4812)))
       Fa6(1565) = p2mup5*F(4813)+p3mup5*F(4814)+p4mup5*F(4815)+2*p1mu
     -   p5*F(4816)
       Fa6(1566) = p1mup5*F(4817)+p3mup5*F(4818)+p4mup5*F(4819)+2*p2mu
     -   p5*F(4820)
       Fa6(1567) = p1mup5*F(4821)+p2mup5*F(4822)+p4mup5*F(4823)+2*p3mu
     -   p5*F(4824)
       Fa6(1568) = p1mup5*F(4825)+p2mup5*F(4826)+p3mup5*F(4827)+2*p4mu
     -   p5*F(4828)
       Fa6(1569) = p2mup1*F(4526)+p3mup1*F(4612)+p4mup1*F(4614)
       Fa6(1570) = p3mup1*F(4532)+p1mup1*F(4607)+p4mup1*F(4615)
       Fa6(1571) = p1mup1*F(4608)+p2mup1*F(4610)+p4mup1*F(4616)
       Fa6(1572) = p1mup1*F(4609)+p2mup1*F(4611)+p3mup1*F(4613)
       Fa6(1573) = p1mup4*F(4830)+p2mup4*F(4831)+p3mup4*F(4832)+p4mup4
     -   *F(4833)
       Fa6(1574) = 4*(p1mup5*(2*(p2mup4*F(4624)+p3mup4*F(4625))+p4mup4
     -   *F(4627))-p2mup4*p4mup5*F(4628)-p2mup5*p4mup4*F(4629)+p3mup4
     -   *p4mup5*F(4633)+p3mup5*p4mup4*F(4634)+p1mup4*(p4mup5*F(4626)
     -   +2*(p1mup5*F(4623)-p2mup5*F(4630)+p3mup5*F(4635)))-2*(p2mup4
     -   *p2mup5*F(4631)+p2mup5*p3mup4*F(4632)-p2mup4*p3mup5*F(4636)-
     -   p3mup4*p3mup5*F(4637))+p4mup4*p4mup5*F(4638))+mup4mup5*F(483
     -   4)
       Fa6(1575) = 4*(p1mup5*(p2mup4*F(4641)+p3mup4*F(4642))-p2mup5*p3
     -   mup4*F(4644)-p1mup4*(p3mup5*F(4634)+p2mup5*F(4640)+p1mup5*F(
     -   4645)+p4mup5*F(4646))+p2mup4*(p3mup5*F(4643)+p4mup5*F(4647)+
     -   p2mup5*F(4648))-p1mup5*p4mup4*F(4650)+p2mup5*p4mup4*F(4652)-
     -   2*(p3mup4*p4mup5*F(4649)+p4mup4*p4mup5*F(4651)+p3mup5*p4mup4
     -   *F(4653)+p3mup4*p3mup5*F(4654)))+mup4mup5*F(4835)
       Fa6(1576) = p1mup4*F(4836)+p3mup4*F(4837)+p2mup4*F(4838)+p4mup4
     -   *F(4839)
       Fa6(1577) = p1mup4*F(4840)+p2mup4*F(4841)+p3mup4*F(4842)+2*p4mu
     -   p4*F(4843)
       Fa6(1578) = mup4mup5*F(4614)-4*(-(p1mup5*p2mup4*F(4528))+p1mup4
     -   *(p2mup5*F(4663)+p3mup5*F(4664))-p1mup5*p3mup4*F(4665)-p2mup
     -   4*p3mup5*F(4666)+p2mup5*p3mup4*F(4667)-p2mup4*p4mup5*F(4668)
     -   -p3mup4*p4mup5*F(4669)-p2mup4*p2mup5*F(4670)-p3mup4*p3mup5*F
     -   (4671)-p2mup5*p4mup4*F(4672)-p3mup5*p4mup4*F(4673))
       Fa6(1579) = p2mup3*F(4844)+2*p1mup3*F(4845)+p3mup3*F(4846)+p4mu
     -   p3*F(4847)
       Fa6(1580) = p1mup3*F(4848)+2*p2mup3*F(4849)+p3mup3*F(4850)+p4mu
     -   p3*F(4851)
       Fa6(1581) = p1mup3*F(4852)+p2mup3*F(4853)+p4mup3*F(4854)+2*p3mu
     -   p3*F(4855)
       Fa6(1582) = mup3mup5*F(4608)-4*(p1mup5*p4mup3*F(4665)-p2mup5*p4
     -   mup3*F(4688)+p3mup5*p4mup3*F(4689)+p4mup5*(p1mup3*F(4686)-p2
     -   mup3*F(4687)+p4mup3*F(4690)-p3mup3*F(4691)))
       Fa6(1583) = p1mup3*F(4856)+p2mup3*F(4857)+p3mup3*F(4858)-2*p4mu
     -   p3*F(4859)
       Fa6(1584) = mup3mup5*F(4860)+4*(p1mup3*p3mup5*F(4699)+p1mup5*p3
     -   mup3*F(4700)+p2mup5*(p4mup3*F(4552)+p1mup3*F(4697)-p3mup3*F(
     -   4703))+p3mup3*p3mup5*F(4704)+p3mup3*p4mup5*F(4705)+p3mup5*p4
     -   mup3*F(4708)+2*(p1mup3*p1mup5*F(4623)+p1mup3*p4mup5*F(4706)+
     -   p4mup3*p4mup5*F(4707)+p1mup5*p4mup3*F(4709))+p2mup3*(p4mup5*
     -   F(4548)+p1mup5*F(4698)-p3mup5*F(4702)+p2mup5*F(4861)))
       Fa6(1585) = mup3mup5*F(4862)-4*(p2mup5*p4mup3*F(4553)-p4mup3*(2
     -   *(p4mup5*F(4649)+p3mup5*F(4654))+p1mup5*F(4708))-p1mup3*(p4m
     -   up5*F(4705)+p1mup5*F(4711)+p2mup5*F(4712)+p3mup5*F(4713))-p3
     -   mup3*(2*(p4mup5*F(4654)+p3mup5*F(4717))+p1mup5*F(4718))+p2mu
     -   p3*(p4mup5*F(4549)+p2mup5*F(4714)+p3mup5*F(4716)-p1mup5*F(48
     -   63)))+p2mup5*p3mup3*F(4864)
       Fa6(1586) = p2mup4*F(4865)+p3mup4*F(4866)+2*p1mup4*F(4867)+p4mu
     -   p4*F(4868)
       Fa6(1587) = -8*(p2mup3*p2mup4*F(4880)+p2mup4*p3mup3*F(4880)+p2m
     -   up3*p3mup4*F(4880)+p3mup3*p3mup4*F(4880)+p1mup4*(p2mup3+p3mu
     -   p3)*F(4881))+p2mup3*p4mup4*F(4882)+p3mup3*p4mup4*F(4882)-4*(
     -   -(p1mup3*(2*(p1mup4*F(4875)+(p2mup4+p3mup4)*F(4876))+p4mup4*
     -   F(4877)))+p1mup4*p4mup3*F(4878)+p2mup4*p4mup3*F(4879)+p3mup4
     -   *p4mup3*F(4879)-p4mup3*p4mup4*F(4883))
       Fa6(1588) = p1mup3*(2*(p1mup4*F(4884)+p2mup4*F(4885)+p3mup4*F(4
     -   886))+p4mup4*F(4887))-p1mup4*p4mup3*F(4888)-p2mup3*p4mup4*F(
     -   4892)+p2mup4*p4mup3*F(4893)-p3mup4*p4mup3*F(4894)-2*(p2mup3*
     -   p2mup4*F(4889)+p1mup4*p2mup3*F(4890)+p2mup3*p3mup4*F(4891)+p
     -   3mup3*p3mup4*F(4895)+p1mup4*p3mup3*F(4896)+p2mup4*p3mup3*F(4
     -   897))-p3mup3*p4mup4*F(4898)+p4mup3*p4mup4*F(4899)
       Fa6(1589) = -8*(p2mup3*p2mup4*F(4900)+p2mup3*p3mup4*F(4901)-p3m
     -   up3*p3mup4*F(4905)-p2mup4*p3mup3*F(4906))+p1mup4*(4*p4mup3*F
     -   (4903)-8*(p2mup3*F(4900)-p3mup3*F(4906)))-4*(p2mup3*p4mup4*F
     -   (4902)+p1mup3*(2*((p1mup4+p2mup4)*F(4900)+p3mup4*F(4901))+p4
     -   mup4*F(4902))-p2mup4*p4mup3*F(4903)-p3mup4*p4mup3*F(4904)-p3
     -   mup3*p4mup4*F(4907))+p4mup3*p4mup4*F(4908)
       Fa6(1590) = p2mup3*p4mup4*F(4919)-4*(p2mup3*p2mup4*F(4914)+p2mu
     -   p3*p3mup4*F(4915)+p2mup4*p3mup3*F(4916)-p3mup3*p3mup4*F(4917
     -   )-p1mup4*(p2mup3*F(4911)+p3mup3*F(4913)-2*p4mup3*F(4918))-p1
     -   mup3*(p2mup4*F(4910)+p3mup4*F(4912)+2*(p1mup4*F(4909)-p4mup4
     -   *F(4918)))-p2mup4*p4mup3*F(4920)+p3mup3*p4mup4*F(4921)+p3mup
     -   4*p4mup3*F(4922))-8*p4mup3*p4mup4*F(4923)
       Fa6(1591) = p2mup3*p3mup4*F(4930)+p2mup4*p3mup3*F(4931)+p3mup3*
     -   p3mup4*F(4932)+p1mup3*(p1mup4*F(4924)+p2mup4*F(4926)+p3mup4*
     -   F(4928)+p4mup4*F(4933))-p1mup4*(p2mup3*F(4927)+p3mup3*F(4929
     -   )+p4mup3*F(4934))+p3mup4*p4mup3*F(4936)-p3mup3*p4mup4*F(4938
     -   )+2*(p2mup3*p2mup4*F(4925)+p2mup4*p4mup3*F(4935)-p2mup3*p4mu
     -   p4*F(4937)+p4mup3*p4mup4*F(4939))
       Fa6(1592) = p2mup4*p3mup3*F(4948)+4*(p1mup4*p3mup3*F(4947)+p1mu
     -   p3*(p1mup4*F(4940)+p2mup4*F(4941)+p3mup4*F(4942)-p4mup4*F(49
     -   49))-p1mup4*p4mup3*F(4950)-p2mup4*p4mup3*F(4951)-p2mup3*(p2m
     -   up4*F(4943)-p1mup4*F(4944)+p3mup4*F(4945)-p4mup4*F(4953)))+8
     -   *(p3mup3*p3mup4*F(4946)-p3mup4*p4mup3*F(4952)-p4mup3*p4mup4*
     -   F(4954)+p3mup3*p4mup4*F(4955))
       Fa6(1593) = mup3mup4*((p1mup5+p2mup5)*F(4534)+p3mup5*F(4535)+p4
     -   mup5*F(4958))-mup3mup5*((p1mup4+p2mup4)*F(4956)+p3mup4*F(495
     -   7)-p4mup4*F(4959))+mup4mup5*((p1mup3+p2mup3)*F(4533)-p3mup3*
     -   F(4537)+p4mup3*F(4960))
       Fa6(1594) = mup3mup4*(p1mup5*F(4528)+(p2mup5+p3mup5)*F(4529)+p4
     -   mup5*F(4668))+mup3mup5*(-(p1mup4*F(4663))+p4mup4*F(4672)-(p2
     -   mup4+p3mup4)*F(4961))+mup4mup5*(p1mup3*F(4527)-(p2mup3+p3mup
     -   3)*F(4531)+p4mup3*F(4962))
       Fa6(1595) = mup3mup4*(p1mup5*F(4665)+p4mup5*F(4669)-p2mup5*F(46
     -   88)+p3mup5*F(4689))+mup3mup5*(-(p1mup4*F(4664))+p4mup4*F(467
     -   3)+p2mup4*F(4963)-p3mup4*F(4964))+mup4mup5*(p1mup3*F(4686)-p
     -   2mup3*F(4687)-p3mup3*F(4691)+p4mup3*F(4965))
       Fa6(1596) = 8*p1mup4*(p2mup3+p3mup3)*F(4881)-4*(p1mup3*(2*(p1mu
     -   p4*F(4875)+(p2mup4+p3mup4)*F(4876))+p4mup4*F(4877))-p1mup4*p
     -   4mup3*F(4878)-(p2mup4+p3mup4)*(p4mup3*F(4879)+2*(p2mup3+p3mu
     -   p3)*F(4880))-(p2mup3+p3mup3)*p4mup4*F(4966))+p4mup3*p4mup4*F
     -   (4967)
       Fa6(1597) = -(p1mup3*(2*(p1mup4*F(4884)+p2mup4*F(4885)+p3mup4*F
     -   (4886))+p4mup4*F(4887)))+p2mup3*p4mup4*F(4892)-p2mup4*p4mup3
     -   *F(4893)+p3mup4*p4mup3*F(4894)+p1mup4*(p4mup3*F(4888)+2*(p2m
     -   up3*F(4890)+p3mup3*F(4896)))+2*(p2mup3*p2mup4*F(4889)+p2mup3
     -   *p3mup4*F(4891)+p3mup3*p3mup4*F(4895)+p2mup4*p3mup3*F(4897))
     -   +p3mup3*p4mup4*F(4898)-p4mup3*p4mup4*F(4899)
       Fa6(1598) = 8*(p2mup3*p2mup4*F(4900)+p2mup3*p3mup4*F(4901)-p3mu
     -   p3*p3mup4*F(4905)-p2mup4*p3mup3*F(4906))+p1mup4*(-4*p4mup3*F
     -   (4903)+8*(p2mup3*F(4900)-p3mup3*F(4906)))+p3mup3*p4mup4*F(49
     -   68)+4*(p2mup3*p4mup4*F(4902)+p1mup3*(2*((p1mup4+p2mup4)*F(49
     -   00)+p3mup4*F(4901))+p4mup4*F(4902))-p2mup4*p4mup3*F(4903)-p3
     -   mup4*p4mup3*F(4904)+p4mup3*p4mup4*F(4969))
       Fa6(1599) = -4*(-(p2mup3*p2mup4*F(4914))-p2mup3*p3mup4*F(4915)-
     -   p2mup4*p3mup3*F(4916)+p3mup3*p3mup4*F(4917)+p1mup4*(p2mup3*F
     -   (4911)+p3mup3*F(4913)-2*p4mup3*F(4918))+p1mup3*(p2mup4*F(491
     -   0)+p3mup4*F(4912)+2*(p1mup4*F(4909)-p4mup4*F(4918)))-p3mup3*
     -   p4mup4*F(4921)-p3mup4*p4mup3*F(4922)-2*p4mup3*p4mup4*F(4923)
     -   -p2mup3*p4mup4*F(4970))+p2mup4*p4mup3*F(4971)
       Fa6(1600) = 8*(p2mup4*p4mup3*F(4935)+p4mup3*p4mup4*F(4939))+p1m
     -   up4*p4mup3*F(4974)+4*(p2mup4*p3mup3*F(4931)+p3mup3*p3mup4*F(
     -   4932)+p1mup3*(p1mup4*F(4924)+p2mup4*F(4926)+p3mup4*F(4928)+p
     -   4mup4*F(4933))+p3mup4*p4mup3*F(4936)+p2mup3*(p3mup4*F(4930)+
     -   2*(p2mup4*F(4925)-p4mup4*F(4937))+p1mup4*F(4972))+p1mup4*p3m
     -   up3*F(4973)+p3mup3*p4mup4*F(4975))
       Fa6(1601) = -(p1mup4*p3mup3*F(4947))-p1mup3*(p1mup4*F(4940)+p2m
     -   up4*F(4941)+p3mup4*F(4942)-p4mup4*F(4949))+p1mup4*p4mup3*F(4
     -   950)+p2mup4*p4mup3*F(4951)-2*(p3mup3*p3mup4*F(4946)-p3mup4*p
     -   4mup3*F(4952)-p4mup3*p4mup4*F(4954)+p3mup3*p4mup4*F(4955))+p
     -   2mup3*(p2mup4*F(4943)+p3mup4*F(4945)-p4mup4*F(4953)+p1mup4*F
     -   (4976))+p2mup4*p3mup3*F(4977)
       Fa6(1602) = p3mup1*((p1mup5+p2mup5)*F(4534)+p3mup5*F(4535)+p4mu
     -   p5*F(4958))
       Fa6(1603) = p2mup1*(p2mup5+p3mup5)*F(4529)+p1mup5*(p2mup1*F(452
     -   8)+p3mup1*F(4665))+p2mup1*p4mup5*F(4668)+p3mup1*(p4mup5*F(46
     -   69)-p2mup5*F(4688)+p3mup5*F(4689))
       Fa6(1604) = p2mup1*((p1mup5+p2mup5)*F(4534)+p3mup5*F(4535)+p4mu
     -   p5*F(4958))
       Fa6(1605) = p1mup4*(p2mup1*F(4663)+p3mup1*F(4664))-p2mup1*p4mup
     -   4*F(4672)+p2mup1*(p2mup4+p3mup4)*F(4961)-p3mup1*(p4mup4*F(46
     -   73)+p2mup4*F(4963)-p3mup4*F(4964))
       Fa6(1606) = p3mup1*((p1mup4+p2mup4)*F(4956)+p3mup4*F(4957)-p4mu
     -   p4*F(4959))
       Fa6(1607) = p2mup1*((p1mup4+p2mup4)*F(4956)+p3mup4*F(4957)-p4mu
     -   p4*F(4959))
       Fa6(1608) = mup1mup3*(p1mup5*F(4528)+(p2mup5+p3mup5)*F(4529)+p4
     -   mup5*F(4668))
       Fa6(1609) = mup1mup3*(p1mup5*F(4665)+p4mup5*F(4669)-p2mup5*F(46
     -   88)+p3mup5*F(4689))
       Fa6(1610) = mup1mup3*((p1mup5+p2mup5)*F(4534)+p3mup5*F(4535)+p4
     -   mup5*F(4958))
       Fa6(1611) = mup1mup5*F(4607)-4*p4mup1*(p1mup5*F(4528)+(p2mup5+p
     -   3mup5)*F(4529)+p4mup5*F(4668))
       Fa6(1612) = mup1mup5*F(4610)+4*p4mup1*((p1mup5+p2mup5)*F(4534)+
     -   p3mup5*F(4535)+p4mup5*F(4958))
       Fa6(1613) = p1mup1*(p1mup4*F(4663)-p4mup4*F(4672)+(p2mup4+p3mup
     -   4)*F(4961))
       Fa6(1614) = p1mup1*(p1mup4*F(4664)-p4mup4*F(4673)-p2mup4*F(4963
     -   )+p3mup4*F(4964))
       Fa6(1615) = -(p2mup1*(p2mup3+p3mup3)*F(4531))+p1mup3*(p2mup1*F(
     -   4527)+p3mup1*F(4686))-p3mup1*(p2mup3*F(4687)+p3mup3*F(4691))
       Fa6(1616) = -(p1mup1*p1mup3*F(4527))+p1mup1*(p2mup3+p3mup3)*F(4
     -   531)-(p1mup3+p2mup3)*p3mup1*F(4533)+p3mup1*p3mup3*F(4537)
       Fa6(1617) = p1mup3*p2mup1*F(4533)+p2mup1*p2mup3*F(4533)-p2mup1*
     -   p3mup3*F(4537)-p1mup1*p1mup3*F(4686)+p1mup1*p2mup3*F(4687)+p
     -   1mup1*p3mup3*F(4691)
       Fa6(1618) = mup1mup5*F(4608)-4*p4mup1*(p1mup5*F(4665)+p4mup5*F(
     -   4669)-p2mup5*F(4688)+p3mup5*F(4689))
       Fa6(1619) = p4mup1*(p1mup4*F(4664)-p4mup4*F(4673)-p2mup4*F(4963
     -   )+p3mup4*F(4964))
       Fa6(1620) = mup1mup4*(p1mup5*F(4665)+p4mup5*F(4669)-p2mup5*F(46
     -   88)+p3mup5*F(4689))
       Fa6(1621) = p4mup1*(p1mup4*F(4663)-p4mup4*F(4672)+(p2mup4+p3mup
     -   4)*F(4961))
       Fa6(1622) = -(mup1mup4*p1mup3*F(4527))+mup1mup4*(p2mup3+p3mup3)
     -   *F(4531)+mup1mup3*(p1mup4*F(4663)-p4mup4*F(4672)+(p2mup4+p3m
     -   up4)*F(4961))
       Fa6(1623) = mup1mup4*(-(p1mup3*F(4686))+p2mup3*F(4687)+p3mup3*F
     -   (4691))+mup1mup3*(p1mup4*F(4664)-p4mup4*F(4673)-p2mup4*F(496
     -   3)+p3mup4*F(4964))
       Fa6(1624) = p4mup1*((p1mup4+p2mup4)*F(4956)+p3mup4*F(4957)-p4mu
     -   p4*F(4959))
       Fa6(1625) = mup1mup4*((p1mup5+p2mup5)*F(4534)+p3mup5*F(4535)+p4
     -   mup5*F(4958))
       Fa6(1626) = mup1mup4*(p1mup5*F(4528)+(p2mup5+p3mup5)*F(4529)+p4
     -   mup5*F(4668))
       Fa6(1627) = -(mup1mup4*(p1mup3+p2mup3)*F(4533))+mup1mup4*p3mup3
     -   *F(4537)+mup1mup3*((p1mup4+p2mup4)*F(4956)+p3mup4*F(4957)-p4
     -   mup4*F(4959))
       Fa6(1628) = mup1mup5*F(4526)+4*p4mup1*(p1mup5*F(4528)+(p2mup5+p
     -   3mup5)*F(4529)+p4mup5*F(4668))
       Fa6(1629) = mup1mup5*F(4532)-4*p4mup1*((p1mup5+p2mup5)*F(4534)+
     -   p3mup5*F(4535)+p4mup5*F(4958))
       Fa6(1630) = p2mup1*(p2mup3+p3mup3)*F(4531)-p1mup3*(p2mup1*F(452
     -   7)+p3mup1*F(4686))+p3mup1*(p2mup3*F(4687)+p3mup3*F(4691))
       Fa6(1631) = p1mup1*p1mup3*F(4527)-p1mup1*(p2mup3+p3mup3)*F(4531
     -   )+(p1mup3+p2mup3)*p3mup1*F(4533)-p3mup1*p3mup3*F(4537)
       Fa6(1632) = mup1mup5*F(4612)+4*p4mup1*(p1mup5*F(4665)+p4mup5*F(
     -   4669)-p2mup5*F(4688)+p3mup5*F(4689))
       Fa6(1633) = p4mup3*(p2mup1*F(4962)+p3mup1*F(4965))
       Fa6(1634) = p4mup3*(p3mup1*F(4960)+p1mup1*F(4962))
       Fa6(1635) = p4mup3*(p2mup1*F(4960)-p1mup1*F(4965))
       Fa6(1636) = p4mup1*(-(p1mup3*F(4527))+(p2mup3+p3mup3)*F(4531)-p
     -   4mup3*F(4962))
       Fa6(1637) = p4mup1*(p1mup3*F(4686)-p2mup3*F(4687)-p3mup3*F(4691
     -   )+p4mup3*F(4965))
       Fa6(1638) = p4mup1*((p1mup3+p2mup3)*F(4533)-p3mup3*F(4537)+p4mu
     -   p3*F(4960))
       Fa6(1639) = p4mup1*(p1mup3*F(4527)-(p2mup3+p3mup3)*F(4531)+p4mu
     -   p3*F(4962))
       Fa6(1640) = p1mup1*(p2mup5+p3mup5)
       Fa6(1641) = p1mup1*(p2mup5*F(4688)-p3mup5*F(4689))
       Fa6(1642) = mup1mup5*(p1mup4*F(4664)-p4mup4*F(4673)-p2mup4*F(49
     -   63)+p3mup4*F(4964))
       Fa6(1643) = mup1mup5*(p1mup3*F(4527)-(p2mup3+p3mup3)*F(4531)+p4
     -   mup3*F(4962))
       Fa6(1644) = mup1mup5*(p1mup3*F(4686)-p2mup3*F(4687)-p3mup3*F(46
     -   91)+p4mup3*F(4965))
       Fa6(1645) = mup1mup5*((p1mup4+p2mup4)*F(4956)+p3mup4*F(4957)-p4
     -   mup4*F(4959))
       Fa6(1646) = mup1mup5*(p1mup4*F(4663)-p4mup4*F(4672)+(p2mup4+p3m
     -   up4)*F(4961))
       Fa6(1647) = mup1mup5*((p1mup3+p2mup3)*F(4533)-p3mup3*F(4537)+p4
     -   mup3*F(4960))
       Fa6(1648) = mup1mup5*(-(p1mup3*F(4527))+(p2mup3+p3mup3)*F(4531)
     -   -p4mup3*F(4962))
       Fa6(1649) = mup4mup5*F(4981)+p1mup5*p3mup4*F(4985)-8*(p1mup5*p2
     -   mup4*F(4984)+p2mup4*p2mup5*F(4991))+p2mup5*p3mup4*F(4993)+p3
     -   mup4*p4mup5*F(4994)+p3mup5*p4mup4*F(4995)+4*(p1mup5*p4mup4*F
     -   (4987)+p2mup4*p4mup5*F(4988)+p2mup5*p4mup4*F(4989)-p1mup4*(p
     -   3mup5*F(4983)+p4mup5*F(4986)+2*(p1mup5*F(4982)+p2mup5*F(4990
     -   )))-p2mup4*p3mup5*F(4992)-p3mup4*p3mup5*F(4996)+p4mup4*p4mup
     -   5*F(4997))
       Fa6(1650) = mup4mup5*F(4998)+p1mup4*(p2mup5*F(4999)-4*p1mup5*F(
     -   5001)+p3mup5*F(5002)+p4mup5*F(5003))+p3mup5*p4mup4*F(5006)-2
     -   *p1mup5*p4mup4*F(5008)+p2mup5*p4mup4*F(5012)-4*(p1mup5*p2mup
     -   4*F(5000)+p3mup4*p4mup5*F(5005)-p1mup5*p3mup4*F(5007)+p2mup4
     -   *(p4mup5*F(5004)+2*(p3mup5*F(5009)+p2mup5*F(5013))))+8*(p3mu
     -   p4*p3mup5*F(5010)-p4mup4*p4mup5*F(5011)+p2mup5*p3mup4*F(5014
     -   ))
       Fa6(1651) = mup4mup5*F(5015)+p1mup5*p3mup4*F(5019)+2*p1mup4*(p3
     -   mup5*F(5018)+2*(p2mup5*F(5016)+p1mup5*F(5020)+p4mup5*F(5021)
     -   ))+p3mup4*p4mup5*F(5023)-4*(p1mup5*p2mup4*F(5017)-p3mup5*p4m
     -   up4*F(5024)-p1mup5*p4mup4*F(5025)-p2mup4*p3mup5*F(5026))+8*(
     -   p2mup4*p4mup5*F(5022)+p3mup4*p3mup5*F(5027)+p4mup4*p4mup5*F(
     -   5028)+p2mup5*p4mup4*F(5029)+p2mup4*p2mup5*F(5030))+p2mup5*p3
     -   mup4*F(5031)
       Fa6(1652) = p2mup5*F(5032)+p3mup5*F(5033)+p4mup5*F(5034)+2*p1mu
     -   p5*F(5035)
       Fa6(1653) = p1mup5*F(5036)+p4mup5*F(5037)+p3mup5*F(5038)+2*p2mu
     -   p5*F(5039)
       Fa6(1654) = p1mup5*F(5040)+p2mup5*F(5041)+p4mup5*F(5042)-2*p3mu
     -   p5*F(5043)
       Fa6(1655) = p1mup5*F(5044)+p2mup5*F(5045)+p3mup5*F(5046)+2*p4mu
     -   p5*F(5047)
       Fa6(1656) = p2mup1*F(5048)+p3mup1*F(5049)+p4mup1*F(5050)
       Fa6(1657) = p1mup1*F(4978)+p3mup1*F(5051)+p4mup1*F(5052)
       Fa6(1658) = p1mup1*F(5053)+p2mup1*F(5054)+p4mup1*F(5055)
       Fa6(1659) = p2mup1*F(5056)+p1mup1*F(5057)+p3mup1*F(5058)
       Fa6(1660) = p1mup4*F(5059)+p2mup4*F(5060)+p3mup4*F(5061)+p4mup4
     -   *F(5062)
       Fa6(1661) = p1mup4*F(5063)+p2mup4*F(5064)+p3mup4*F(5065)+p4mup4
     -   *F(5066)
       Fa6(1662) = 8*p3mup4*p3mup5*F(4563)+mup4mup5*F(5067)+2*(p2mup4*
     -   p4mup5*F(5070)+p2mup5*p4mup4*F(5071))+4*(p3mup4*p4mup5*F(456
     -   0)+p3mup5*p4mup4*F(4561)+p4mup4*p4mup5*F(4564)+p1mup4*(p4mup
     -   5*F(4558)+2*(p1mup5*F(4556)+p3mup5*F(4562))+p2mup5*F(5068))+
     -   p1mup5*(2*p3mup4*F(4557)+p4mup4*F(4559)-p2mup4*F(5069))+p2mu
     -   p4*p2mup5*F(5072)-p2mup4*p3mup5*F(5073)+p2mup5*p3mup4*F(5074
     -   ))
       Fa6(1663) = -8*(p4mup4*p4mup5*F(4587)+p3mup5*p4mup4*F(4589)+p2m
     -   up4*p2mup5*F(4590))+mup4mup5*F(5075)-2*(p1mup4*(2*(p3mup5*F(
     -   4561)+p1mup5*F(4581)+p4mup5*F(4582))+p2mup5*F(5076))-p1mup5*
     -   p2mup4*F(5077))-4*(p1mup5*p4mup4*F(4586)-p2mup4*p4mup5*F(507
     -   8)+p2mup5*p4mup4*F(5079)-p2mup4*p3mup5*F(5080)+p3mup4*(2*(p4
     -   mup5*F(4585)+p3mup5*F(4590))-p1mup5*F(4812)+p2mup5*F(5081)))
       Fa6(1664) = p1mup2*F(5083)+p2mup2*F(5084)+p3mup2*F(5085)+p4mup2
     -   *F(5086)
       Fa6(1665) = p1mup2*F(5087)+p3mup2*F(5088)+p4mup2*F(5089)+2*p2mu
     -   p2*F(5090)
       Fa6(1666) = p1mup2*F(5091)+p2mup2*F(5092)+p4mup2*F(5093)+p3mup2
     -   *F(5094)
       Fa6(1667) = p1mup2*F(5096)+p2mup2*F(5097)+p3mup2*F(5098)+p4mup2
     -   *F(5099)
       Fa6(1668) = p2mup4*F(5102)+p3mup4*F(5103)+p1mup4*F(5104)+p4mup4
     -   *F(5105)
       Fa6(1669) = p2mup4*F(5109)+p1mup4*F(5110)+p3mup4*F(5111)+2*p4mu
     -   p4*F(5112)
       Fa6(1670) = mup4mup5*F(5050)+p2mup4*p3mup5*F(5120)+p2mup5*p3mup
     -   4*F(5121)-4*(p1mup4*p3mup5*F(4663)-p3mup4*p4mup5*F(4668)-p3m
     -   up4*p3mup5*F(4670)-p3mup5*p4mup4*F(4672)-p1mup4*p2mup5*F(511
     -   6)-p1mup5*(p3mup4*F(4528)+p2mup4*F(5117))-p2mup4*p4mup5*F(51
     -   18)-p2mup4*p2mup5*F(5119)+p2mup5*p4mup4*F(5122))
       Fa6(1671) = 2*(p3mup5*p4mup4-p3mup4*p4mup5)*F(352)-4*(p1mup5*p3
     -   mup4*F(276)-p1mup4*p3mup5*F(276)+(p2mup5*p3mup4-p2mup4*p3mup
     -   5)*F(359))+mup4mup5*F(5126)
       Fa6(1672) = mup4mup5*F(5127)+2*(p1mup5*(2*(p2mup4*F(359)+p3mup4
     -   *F(888))+p4mup4*F(7441))-p1mup4*(2*(p2mup5*F(359)+p3mup5*F(8
     -   88))+p4mup5*F(7441)))+4*((-(p3mup5*p4mup4)+p3mup4*p4mup5)*F(
     -   7442)-(p2mup5*p4mup4-p2mup4*p4mup5)*F(7448))
       Fa6(1673) = 2*(p1mup5*p3mup4*F(352)-p1mup4*p3mup5*F(352))+mup4m
     -   up5*F(5128)+4*((-(p3mup5*p4mup4)+p3mup4*p4mup5)*F(7444)+(p2m
     -   up5*p3mup4-p2mup4*p3mup5)*F(7448))
       Fa6(1674) = p2mup5*F(5129)+p3mup5*F(5130)-(2*p1mup5*F(276)+p4mu
     -   p5*F(352))*P(1)
       Fa6(1675) = p1mup5*F(5131)+p3mup5*F(5132)+2*p2mup5*F(5133)+p4mu
     -   p5*F(5134)
       Fa6(1676) = p1mup5*F(5135)-p2mup5*F(5136)+2*p3mup5*F(5137)+p4mu
     -   p5*F(5138)
       Fa6(1677) = p2mup5*F(5139)+p3mup5*F(5140)-(p1mup5*F(352)+2*p4mu
     -   p5*F(7444))*P(1)
       Fa6(1678) = p2mup1*F(5126)+p3mup1*F(5141)+p4mup1*F(7441)*P(1)
       Fa6(1679) = p1mup1*F(5123)+p4mup1*F(5128)+p3mup1*F(5142)
       Fa6(1680) = p4mup1*F(374)+p1mup1*F(5143)+p2mup1*F(5144)
       Fa6(1681) = p3mup1*F(368)+p2mup1*F(5125)-p1mup1*F(7441)*P(1)
       Fa6(1682) = p1mup4*F(5145)+p3mup4*F(5146)+2*p2mup4*F(5147)+p4mu
     -   p4*F(5148)
       Fa6(1683) = 2*p3mup4*F(5137)+p4mup4*F(5138)+p1mup4*F(5149)+p2mu
     -   p4*F(5150)
       Fa6(1684) = -2*(p2mup5*p4mup4-p2mup4*p4mup5)*F(352)+4*(p1mup5*p
     -   2mup4*F(276)-p1mup4*p2mup5*F(276)-(p2mup5*p3mup4-p2mup4*p3mu
     -   p5)*F(888))+mup4mup5*F(5141)
       Fa6(1685) = 2*p1mup2*F(5152)+p2mup2*F(5153)+p3mup2*F(5154)+p4mu
     -   p2*F(5155)
       Fa6(1686) = p1mup2*F(5156)+p3mup2*F(5157)+2*p2mup2*F(5158)+p4mu
     -   p2*F(5159)
       Fa6(1687) = -(p1mup2*F(5160))+p2mup2*F(5161)-2*p3mup2*F(5162)+p
     -   4mup2*F(5163)
       Fa6(1688) = p1mup2*F(5165)+p2mup2*F(5166)+p3mup2*F(5167)+2*p4mu
     -   p2*F(5168)
       Fa6(1689) = p2mup4*F(5169)+p3mup4*F(5170)+(2*p1mup4*F(276)+p4mu
     -   p4*F(352))*P(1)
       Fa6(1690) = p2mup4*F(5173)+p3mup4*F(5174)+(p1mup4*F(352)+2*p4mu
     -   p4*F(7444))*P(1)
       Fa6(1691) = -2*(p2mup5*p3mup4-p2mup4*p3mup5)+mup4mup5*P(1)
       Fa6(1692) = 2*(p3mup5*p4mup4-p3mup4*p4mup5)*F(352)-4*(p1mup5*p3
     -   mup4*F(276)-p1mup4*p3mup5*F(276)+(p2mup5*p3mup4-p2mup4*p3mup
     -   5)*F(359))+mup4mup5*F(5178)
       Fa6(1693) = mup4mup5*F(5179)+2*(p1mup5*(2*(p2mup4*F(359)+p3mup4
     -   *F(888))+p4mup4*F(7441))-p1mup4*(2*(p2mup5*F(359)+p3mup5*F(8
     -   88))+p4mup5*F(7441)))+4*((-(p3mup5*p4mup4)+p3mup4*p4mup5)*F(
     -   7442)-(p2mup5*p4mup4-p2mup4*p4mup5)*F(7448))
       Fa6(1694) = 2*(p1mup5*p3mup4*F(352)-p1mup4*p3mup5*F(352))+mup4m
     -   up5*F(5180)+4*((-(p3mup5*p4mup4)+p3mup4*p4mup5)*F(7444)+(p2m
     -   up5*p3mup4-p2mup4*p3mup5)*F(7448))
       Fa6(1695) = p2mup5*F(5181)+p3mup5*F(5182)-(2*p1mup5*F(276)+p4mu
     -   p5*F(352))*P(1)
       Fa6(1696) = p1mup5*F(5183)+p3mup5*F(5184)+2*p2mup5*F(5185)+p4mu
     -   p5*F(5186)
       Fa6(1697) = -(p1mup5*F(5187))+p2mup5*F(5188)+2*p3mup5*F(5189)+p
     -   4mup5*F(5190)
       Fa6(1698) = p3mup5*F(5140)+p2mup5*F(5191)-(p1mup5*F(352)+2*p4mu
     -   p5*F(7444))*P(1)
       Fa6(1699) = p2mup1*F(5178)+p3mup1*F(5192)+p4mup1*F(7441)*P(1)
       Fa6(1700) = p1mup1*F(5175)+p4mup1*F(5180)+p3mup1*F(5193)
       Fa6(1701) = p4mup1*F(374)+p1mup1*F(5194)+p2mup1*F(5195)
       Fa6(1702) = p3mup1*F(368)+p2mup1*F(5177)-p1mup1*F(7441)*P(1)
       Fa6(1703) = p1mup4*F(5196)+p3mup4*F(5197)+2*p2mup4*F(5198)+p4mu
     -   p4*F(5199)
       Fa6(1704) = -2*p3mup4*F(5189)+p1mup4*F(5200)+p2mup4*F(5201)+p4m
     -   up4*F(5202)
       Fa6(1705) = -2*(p2mup5*p4mup4-p2mup4*p4mup5)*F(352)+4*(p1mup5*p
     -   2mup4*F(276)-p1mup4*p2mup5*F(276)-(p2mup5*p3mup4-p2mup4*p3mu
     -   p5)*F(888))+mup4mup5*F(5192)
       Fa6(1706) = 2*p1mup2*F(5204)+p2mup2*F(5205)+p3mup2*F(5206)+p4mu
     -   p2*F(5207)
       Fa6(1707) = p1mup2*F(5208)+p3mup2*F(5209)+2*p2mup2*F(5210)+p4mu
     -   p2*F(5211)
       Fa6(1708) = p1mup2*F(5212)+p2mup2*F(5213)-2*p3mup2*F(5214)+p4mu
     -   p2*F(5215)
       Fa6(1709) = p1mup2*F(5217)+p2mup2*F(5218)+p3mup2*F(5219)+2*p4mu
     -   p2*F(5220)
       Fa6(1710) = p2mup4*F(5221)+p3mup4*F(5222)+(2*p1mup4*F(276)+p4mu
     -   p4*F(352))*P(1)
       Fa6(1711) = p3mup4*F(5174)+p2mup4*F(5225)+(p1mup4*F(352)+2*p4mu
     -   p4*F(7444))*P(1)
       Fa6(1712) = 4*(p1mup5*p3mup4*F(276)+p2mup5*p3mup4*F(359))+2*(p3
     -   mup4*p4mup5*F(441)-p3mup5*(p4mup4*F(430)+2*(p1mup4*F(349)+p2
     -   mup4*F(1050)+p3mup4*F(1051))))+mup4mup5*F(5229)
       Fa6(1713) = -(mup4mup5*F(5230))+2*(p1mup4*(p4mup5*F(432)+2*(p1m
     -   up5*F(355)+p2mup5*F(359)+p3mup5*F(888)))+p1mup5*(p4mup4*F(43
     -   2)+2*(p2mup4*F(356)+p3mup4*F(5231))))+4*(p2mup4*p4mup5*F(433
     -   )+p2mup5*p4mup4*F(434)+p3mup5*p4mup4*F(439)+2*(p2mup4*p2mup5
     -   *F(358)+p2mup5*p3mup4*F(361)+p3mup4*p3mup5*F(362)+p4mup4*p4m
     -   up5*F(436)+p2mup4*p3mup5*F(1554))+p3mup4*p4mup5*F(5232))
       Fa6(1714) = -4*(p2mup5*p3mup4*F(434)-p2mup4*p3mup5*F(434)-(p3mu
     -   p5*p4mup4-p3mup4*p4mup5)*F(436))-2*(p1mup5*p3mup4-p1mup4*p3m
     -   up5)*F(438)+mup4mup5*F(5233)
       Fa6(1715) = p2mup5*F(5234)+p3mup5*F(5235)+(2*p1mup5*F(276)+p4mu
     -   p5*F(441))*P(1)
       Fa6(1716) = p3mup5*F(5236)+p4mup5*F(5237)+p1mup5*F(5238)+2*p2mu
     -   p5*F(5239)
       Fa6(1717) = p2mup5*F(5240)+p4mup5*F(5241)+p1mup5*F(5242)+2*p3mu
     -   p5*F(5243)
       Fa6(1718) = p2mup5*F(5244)+p3mup5*F(5245)+(2*p4mup5*F(436)+p1mu
     -   p5*F(438))*P(1)
       Fa6(1719) = p2mup1*F(5229)+p3mup1*F(5246)-p4mup1*F(1095)*P(1)
       Fa6(1720) = p1mup1*F(5226)+p3mup1*F(5230)+p4mup1*F(5233)
       Fa6(1721) = p4mup1*F(443)+p2mup1*F(5247)+p1mup1*F(5248)
       Fa6(1722) = p3mup1*F(437)+p2mup1*F(5228)+p1mup1*F(1095)*P(1)
       Fa6(1723) = p3mup4*F(5249)+2*p2mup4*F(5250)+p1mup4*F(5251)+p4mu
     -   p4*F(5252)
       Fa6(1724) = p2mup4*F(5253)+p1mup4*F(5254)+2*p3mup4*F(5255)+p4mu
     -   p4*F(5256)
       Fa6(1725) = 2*(p2mup5*p4mup4*F(430)+p2mup4*(p4mup5*F(429)+2*(p2
     -   mup5*F(348)-p3mup5*F(888))))+mup4mup5*F(5246)-4*(p1mup5*p2mu
     -   p4*F(276)-p1mup4*p2mup5*F(349)-p2mup5*p3mup4*F(5257))
       Fa6(1726) = 2*p1mup2*F(5259)+p2mup2*F(5260)+p3mup2*F(5261)+p4mu
     -   p2*F(5262)
       Fa6(1727) = p3mup2*F(5263)+p1mup2*F(5264)+p4mup2*F(5265)+2*p2mu
     -   p2*F(5266)
       Fa6(1728) = p1mup2*F(5267)+p2mup2*F(5268)+p4mup2*F(5269)+2*p3mu
     -   p2*F(5270)
       Fa6(1729) = p1mup2*F(5272)+p2mup2*F(5273)+p3mup2*F(5274)+2*p4mu
     -   p2*F(5275)
       Fa6(1730) = p2mup4*F(5276)+p3mup4*F(5277)-(2*p1mup4*F(349)+p4mu
     -   p4*F(430))*P(1)
       Fa6(1731) = p2mup4*F(5280)+p3mup4*F(5281)-(2*p4mup4*F(436)+p1mu
     -   p4*F(438))*P(1)
       Fa6(1732) = -4*(p1mup5*p3mup4*F(276)+p2mup5*p3mup4*F(359))+2*(p
     -   3mup4*p4mup5*F(429)+p3mup5*(p4mup4*F(430)+2*(p1mup4*F(349)+p
     -   2mup4*F(1050)+p3mup4*F(1051))))+mup4mup5*F(5226)
       Fa6(1733) = 2*(p1mup4*(p4mup5*F(432)+2*(p1mup5*F(355)+p2mup5*F(
     -   359)+p3mup5*F(888)))+p1mup5*(p4mup4*F(432)+2*(p2mup4*F(356)+
     -   p3mup4*F(5231))))+4*(p2mup4*p4mup5*F(433)+p2mup5*p4mup4*F(43
     -   4)+p3mup5*p4mup4*F(439)+2*(p2mup4*p2mup5*F(358)+p2mup5*p3mup
     -   4*F(361)+p3mup4*p3mup5*F(362)+p4mup4*p4mup5*F(436)+p2mup4*p3
     -   mup5*F(1554))+p3mup4*p4mup5*F(5232))+mup4mup5*F(5247)
       Fa6(1734) = 4*((p2mup5*p3mup4-p2mup4*p3mup5)*F(434)-(p3mup5*p4m
     -   up4-p3mup4*p4mup5)*F(436))+2*(p1mup5*p3mup4-p1mup4*p3mup5)*F
     -   (438)+mup4mup5*F(5228)
       Fa6(1735) = p2mup5*F(5283)+p3mup5*F(5284)+(-2*p1mup5*F(276)+p4m
     -   up5*F(429))*P(1)
       Fa6(1736) = p3mup5*F(5285)+p4mup5*F(5286)+p1mup5*F(5287)+2*p2mu
     -   p5*F(5288)
       Fa6(1737) = 2*p3mup5*F(5243)-p2mup5*F(5289)+p4mup5*F(5290)+p1mu
     -   p5*F(5291)
       Fa6(1738) = p2mup5*F(5280)+p3mup5*F(5281)-(2*p4mup5*F(436)+p1mu
     -   p5*F(438))*P(1)
       Fa6(1739) = p2mup1*F(5226)+p3mup1*F(5248)+p4mup1*F(1095)*P(1)
       Fa6(1740) = p4mup1*F(5228)+p1mup1*F(5229)+p3mup1*F(5247)
       Fa6(1741) = p4mup1*F(437)+p2mup1*F(5230)+p1mup1*F(5246)
       Fa6(1742) = p3mup1*F(443)+p2mup1*F(5233)-p1mup1*F(1095)*P(1)
       Fa6(1743) = p3mup4*F(5292)+2*p2mup4*F(5293)+p1mup4*F(5294)+p4mu
     -   p4*F(5295)
       Fa6(1744) = -(p2mup4*F(5296))+p1mup4*F(5297)+2*p3mup4*F(5298)+p
     -   4mup4*F(5299)
       Fa6(1745) = 4*p1mup5*p2mup4*F(276)+p2mup4*(2*p4mup5*F(441)-4*(p
     -   2mup5*F(348)-p3mup5*F(888)))+mup4mup5*F(5248)-2*p2mup5*(p4mu
     -   p4*F(430)+2*(p1mup4*F(349)+p3mup4*F(5257)))
       Fa6(1746) = 2*p1mup2*F(5301)+p2mup2*F(5302)+p3mup2*F(5303)+p4mu
     -   p2*F(5304)
       Fa6(1747) = p3mup2*F(5305)+p1mup2*F(5306)+p4mup2*F(5307)+2*p2mu
     -   p2*F(5308)
       Fa6(1748) = -(p1mup2*F(5309))+p2mup2*F(5310)+p4mup2*F(5311)-2*p
     -   3mup2*F(5312)
       Fa6(1749) = p1mup2*F(5314)+p2mup2*F(5315)+p3mup2*F(5316)+2*p4mu
     -   p2*F(5317)
       Fa6(1750) = p2mup4*F(5318)+p3mup4*F(5319)+(2*p1mup4*F(349)+p4mu
     -   p4*F(430))*P(1)
       Fa6(1751) = p2mup4*F(5244)+p3mup4*F(5245)+(2*p4mup4*F(436)+p1mu
     -   p4*F(438))*P(1)
       Fa6(1752) = -2*(p3mup5*p4mup4-p3mup4*p4mup5)*F(352)+4*(p1mup5*p
     -   3mup4*F(276)-p1mup4*p3mup5*F(276)+(p2mup5*p3mup4-p2mup4*p3mu
     -   p5)*F(359))+mup4mup5*F(5175)
       Return
       End
