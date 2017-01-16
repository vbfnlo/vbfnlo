c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenOddGSim17(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4
     -   ,p3mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5m
     -   up3,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mu
     -   p3,mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa17)
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
       Complex*16 Fa17(1585:1683)
       COMMON/FaPenOddGSim/Fa
       Fa17(1585) = -4*p1mup2*(p2mup1*F(3311)+p3mup1*F(3491))+2*p3mup2
     -   *(p2mup1*F(3489)-2*p3mup1*F(3494))
       Fa17(1586) = p4mup1*(p1mup2*F(3482)+2*(p3mup2*F(3484)+p4mup2*F(
     -   3485)))
       Fa17(1587) = p3mup1*p4mup2*F(3481)-p3mup2*p4mup1*F(3489)+2*(p1m
     -   up2*p4mup1*F(3311)+p4mup1*p4mup2*F(3490))
       Fa17(1588) = 2*p2mup1*p4mup2*F(3481)-4*p4mup1*(p1mup2*F(3491)+p
     -   3mup2*F(3494)+p4mup2*F(3495))
       Fa17(1589) = mup1mup4*p4mup2*F(3481)-mup2mup4*p4mup1*F(3529)
       Fa17(1590) = p4mup2*(p2mup1*F(3490)+p3mup1*F(3495))
       Fa17(1591) = p4mup1*(p1mup2*F(383)+2*(p3mup2*F(4465)+p4mup2*F(4
     -   467)))
       Fa17(1592) = p3mup1*p4mup2*F(1609)-p3mup2*p4mup1*F(4464)
       Fa17(1593) = -(p2mup1*p4mup2*F(1609))+2*p4mup1*(p1mup2*F(357)+p
     -   3mup2*F(4465)+p4mup2*F(4467))
       Fa17(1594) = p4mup1*(p1mup2*F(443)+2*(p3mup2*F(445)+p4mup2*F(45
     -   9)))
       Fa17(1595) = p3mup2*p4mup1*F(477)-p3mup1*p4mup2*F(1623)
       Fa17(1596) = -(p2mup1*p4mup2*F(1623))+2*(p3mup2*p4mup1*F(445)+p
     -   4mup1*(p4mup2*F(459)+p1mup2*F(1629)))
       Fa17(1597) = -(p3mup2*p4mup1*F(477))+p3mup1*p4mup2*F(1623)
       Fa17(1598) = -(p3mup1*p4mup2*F(1609))+p3mup2*p4mup1*F(4464)
       Fa17(1599) = 2*p4mup1*(p1mup2*F(3482)+2*(p3mup2*F(3484)+p4mup2*
     -   F(3485)))
       Fa17(1600) = 2*(p3mup1*p4mup2*F(3481)-p3mup2*p4mup1*F(3489))+4*
     -   (p1mup2*p4mup1*F(3311)+p4mup1*p4mup2*F(3490))
       Fa17(1601) = -2*p2mup1*p4mup2*F(3481)+4*p4mup1*(p1mup2*F(3491)+
     -   p3mup2*F(3494)+p4mup2*F(3495))
       Fa17(1602) = mup1mup5*(2*(p1mup4*F(3474)+p2mup4*F(3475)+p3mup4*
     -   F(3478))+p4mup4*F(3479))
       Fa17(1603) = mup1mup5*(p1mup4*F(3482)+2*(p2mup4*F(3483)+p3mup4*
     -   F(3484)+p4mup4*F(3485)))
       Fa17(1604) = mup1mup5*(-4*(p1mup4*F(3001)-p4mup4*F(3007)+p2mup4
     -   *F(3486))+p3mup4*F(3527))
       Fa17(1605) = mup1mup5*(p1mup2*F(3310)+p3mup2*F(3526))
       Fa17(1606) = mup1mup5*(4*p1mup2*F(3311)-2*p3mup2*F(3489))
       Fa17(1607) = 2*mup1mup5*(p1mup2*F(3482)+2*p3mup2*F(3484))
       Fa17(1608) = 2*mup1mup5*(2*(p1mup4*F(2525)-p4mup4*F(2535)+p3mup
     -   4*F(2565))+p2mup4*F(3492))
       Fa17(1609) = 4*mup1mup5*(p1mup2*F(3491)+p3mup2*F(3494))-2*mup2m
     -   up5*p2mup1*F(3529)
       Fa17(1610) = mup1mup5*(2*(p1mup4*F(263)+p3mup4*F(273)+p2mup4*F(
     -   379))+p4mup4*F(383))
       Fa17(1611) = mup1mup5*(p1mup2*F(1608)+p3mup2*F(3506))
       Fa17(1612) = mup1mup5*(p4mup4*F(448)+2*(p1mup4*F(376)+p2mup4*F(
     -   377)+p3mup4*F(3135)))
       Fa17(1613) = mup1mup5*(p1mup2*F(1262)+p3mup2*F(2368))
       Fa17(1614) = 2*mup1mup5*(2*(p1mup4*F(3474)+p2mup4*F(3475)+p3mup
     -   4*F(3478))+p4mup4*F(3479))
       Fa17(1615) = 2*mup1mup5*(p1mup4*F(3482)+2*(p2mup4*F(3483)+p3mup
     -   4*F(3484)+p4mup4*F(3485)))
       Fa17(1616) = mup1mup5*(4*(p1mup4*F(3001)-p4mup4*F(3007)+p2mup4*
     -   F(3486))+p3mup4*F(3487))
       Fa17(1617) = mup1mup5*(p1mup2*F(3310)-p3mup2*F(3480))
       Fa17(1618) = mup1mup5*(-4*p1mup2*F(3311)+2*p3mup2*F(3489))
       Fa17(1619) = mup1mup5*(p1mup2*F(3482)+2*p3mup2*F(3484))
       Fa17(1620) = mup1mup5*(2*(p1mup4*F(2525)-p4mup4*F(2535)+p3mup4*
     -   F(2565))+p2mup4*F(3492))
       Fa17(1621) = -4*mup1mup5*(p1mup2*F(3491)+p3mup2*F(3494))+2*mup2
     -   mup5*p2mup1*F(3529)
       Fa17(1622) = p2mup1*F(3530)+p3mup1*F(3536)+p4mup1*F(3537)
       Fa17(1623) = p3mup1*F(3531)+p4mup1*F(3532)
       Fa17(1624) = p2mup1*F(3534)+p4mup1*F(3538)
       Fa17(1625) = p2mup1*F(3535)+p3mup1*F(3539)
       Fa17(1626) = -2*p1mup2*F(3547)+p3mup2*F(3548)+p4mup2*F(3549)
       Fa17(1627) = p1mup2*F(3550)+p3mup2*F(3551)+p4mup2*F(3552)
       Fa17(1628) = p1mup2*F(3553)+p4mup2*F(3554)+p3mup2*F(3555)
       Fa17(1629) = p3mup2*F(3556)+p1mup2*F(3557)-2*p4mup2*F(3558)
       Fa17(1630) = p2mup3*F(3562)+p1mup3*F(3563)
       Fa17(1631) = p1mup3*F(3564)+p2mup3*F(3565)
       Fa17(1632) = p1mup3*F(3566)+p2mup3*F(3567)
       Fa17(1633) = p2mup3*F(3568)+p1mup3*F(3569)
       Fa17(1634) = p2mup4*F(3570)+p3mup4*F(3571)+p1mup4*F(3572)+p4mup
     -   4*F(3573)
       Fa17(1635) = p1mup4*F(3574)+p3mup4*F(3575)+p2mup4*F(3576)+p4mup
     -   4*F(3577)
       Fa17(1636) = p1mup4*F(3578)+p2mup4*F(3579)+p3mup4*F(3580)+p4mup
     -   4*F(3581)
       Fa17(1637) = p2mup4*F(3582)+p3mup4*F(3583)+2*p4mup4*F(3584)+p1m
     -   up4*F(3585)
       Fa17(1638) = p2mup1*F(3594)+p3mup1*F(3600)+p4mup1*F(3601)
       Fa17(1639) = p3mup1*F(3595)+p4mup1*F(3596)
       Fa17(1640) = p2mup1*F(3598)+p4mup1*F(3602)
       Fa17(1641) = p2mup1*F(3599)+p3mup1*F(3603)
       Fa17(1642) = 2*p1mup2*F(3611)+p3mup2*F(3612)+p4mup2*F(3613)
       Fa17(1643) = p1mup2*F(3614)+p3mup2*F(3615)+p4mup2*F(3616)
       Fa17(1644) = p1mup2*F(3617)+p4mup2*F(3618)+2*p3mup2*F(3619)
       Fa17(1645) = p1mup2*F(3620)+p3mup2*F(3621)+2*p4mup2*F(3622)
       Fa17(1646) = p2mup3*F(3623)+2*p1mup3*F(3624)
       Fa17(1647) = 2*p2mup3*F(3625)+p1mup3*F(3626)
       Fa17(1648) = p2mup3*F(3627)+p1mup3*F(3628)
       Fa17(1649) = p2mup3*F(3629)+p1mup3*F(3630)
       Fa17(1650) = p2mup4*F(3631)+2*p1mup4*F(3632)+p3mup4*F(3633)+p4m
     -   up4*F(3634)
       Fa17(1651) = p1mup4*F(3635)+p3mup4*F(3636)-2*p2mup4*F(3637)+p4m
     -   up4*F(3638)
       Fa17(1652) = p2mup4*F(3639)+p1mup4*F(3640)+2*p3mup4*F(3641)+p4m
     -   up4*F(3642)
       Fa17(1653) = p2mup4*F(3643)+p1mup4*F(3644)+p3mup4*F(3645)+2*p4m
     -   up4*F(3646)
       Fa17(1654) = p2mup1*F(3655)+p3mup1*F(3661)+p4mup1*F(3662)
       Fa17(1655) = p3mup1*F(3656)+p4mup1*F(3657)
       Fa17(1656) = p2mup1*F(3659)+p4mup1*F(3663)
       Fa17(1657) = p2mup1*F(3660)+p3mup1*F(3664)
       Fa17(1658) = -2*p1mup2*F(3672)+p3mup2*F(3673)+p4mup2*F(3674)
       Fa17(1659) = p1mup2*F(3675)+p3mup2*F(3676)+p4mup2*F(3677)
       Fa17(1660) = p1mup2*F(3678)+2*p3mup2*F(3679)+p4mup2*F(3680)
       Fa17(1661) = p1mup2*F(3681)+p3mup2*F(3682)+2*p4mup2*F(3683)
       Fa17(1662) = -(p2mup3*F(3687))+2*p1mup3*F(3688)
       Fa17(1663) = p1mup3*F(3689)+2*p2mup3*F(3690)
       Fa17(1664) = p1mup3*F(3691)+p2mup3*F(3692)
       Fa17(1665) = p1mup3*F(3693)+p2mup3*F(3694)
       Fa17(1666) = p2mup4*F(3695)+2*p1mup4*F(3696)+p3mup4*F(3697)+p4m
     -   up4*F(3698)
       Fa17(1667) = p1mup4*F(3699)-2*p2mup4*F(3700)+p3mup4*F(3701)+p4m
     -   up4*F(3702)
       Fa17(1668) = p1mup4*F(3703)+p2mup4*F(3704)+2*p3mup4*F(3705)+p4m
     -   up4*F(3706)
       Fa17(1669) = p1mup4*F(3707)+p2mup4*F(3708)+p3mup4*F(3709)+2*p4m
     -   up4*F(3710)
       Fa17(1670) = p3mup1*F(3717)+p4mup1*F(3718)+p2mup1*F(3719)
       Fa17(1671) = p3mup1*F(3720)+p4mup1*F(3721)
       Fa17(1672) = p4mup1*F(3664)+p2mup1*F(3723)
       Fa17(1673) = p3mup1*F(3663)+p2mup1*F(3724)
       Fa17(1674) = -2*p1mup2*F(3732)+p3mup2*F(3733)+p4mup2*F(3734)
       Fa17(1675) = -(p1mup2*F(3735))+p3mup2*F(3736)+p4mup2*F(3737)
       Fa17(1676) = -(p1mup2*F(3738))-2*p3mup2*F(3739)+p4mup2*F(3740)
       Fa17(1677) = -(p1mup2*F(3741))+p3mup2*F(3742)-2*p4mup2*F(3743)
       Fa17(1678) = 2*p1mup3*F(3696)+p2mup3*F(3747)
       Fa17(1679) = p1mup3*F(3748)+2*p2mup3*F(3749)
       Fa17(1680) = p1mup3*F(3703)+p2mup3*F(3750)
       Fa17(1681) = p1mup3*F(3707)+p2mup3*F(3751)
       Fa17(1682) = 2*p1mup4*F(3688)+p2mup4*F(3752)+p3mup4*F(3753)+p4m
     -   up4*F(3754)
       Fa17(1683) = p2mup4*F(3756)+M1*M3*(p1mup4*F(3755)+p3mup4*F(3757
     -   )+p4mup4*F(3758))
       Return
       End
