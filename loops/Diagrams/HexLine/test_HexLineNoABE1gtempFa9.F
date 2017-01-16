c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineNoAbe1g9(p1mup2,p1mup3,p1mup4,
     -   p1mup6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mu
     -   p6,p4mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p
     -   6mup2,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mu
     -   p4,mup3mup6,mup4mup6,Fa9)
       IMPLICIT NONE
      Complex*16   p1mup2, p1mup3, p1mup4, p1mup6, p2mup2, p2mup3, 
     -          p2mup4, p2mup6, p3mup2, p3mup3, p3mup4, p3mup6, 
     -          p4mup2, p4mup3, p4mup4, p4mup6, p5mup2, p5mup3, 
     -          p5mup4, p5mup6, p6mup2, p6mup3, p6mup4, p6mup6
       Complex*16   mup2mup3, mup2mup4, mup2mup6, mup3mup4, mup3mup6, 
     -          mup4mup6
        common/test_FhlFunctionsNoAbe1g/F
       COMMON/Invariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12,s23,s34,s45
     -   ,s56,s16,s123,s234,s345
       Real*8   p1sq 
       Real*8   p2sq 
       Real*8   p3sq 
       Real*8   p4sq 
       Real*8   p5sq 
       Real*8   p6sq 
       Real*8   s12, s16 
       Real*8   s23 
       Real*8   s34 
       Real*8   s45 
       Real*8   s56 
       Real*8   s123,s234,s345
      Complex*16 Fa(587),F(6395)
      Real*8 P(7129) 
        COMMON/test_PFunctionsNoAbe1g/P
       Complex*16 Fa9(113:126)
       COMMON/test_FaFunctionsNoAbe1g/Fa
       Fa9(113) = p5mup4*(p3mup2*F(104)+s345*(p2mup2*F(106)+p1mup2*F(1
     -   07)+(p4mup2+p5mup2)*F(108)))
       Fa9(114) = p5mup4*((p3mup2+p4mup2+p5mup2)*F(64)+p1mup2*F(65)+p2
     -   mup2*F(66))
       Fa9(115) = mup4mup6*p4mup2*F(37)+p5mup4*(mup2mup6*F(105)+2*(p1m
     -   up2*((p3mup6+p4mup6+p5mup6)*F(46)+p1mup6*F(47)+p2mup6*F(48))
     -   -p2mup2*((p3mup6+p4mup6+p5mup6)*F(51)-p1mup6*F(52)-p2mup6*F(
     -   53))+(p3mup2+p4mup2+p5mup2)*(p2mup6*F(60)+(p3mup6+p4mup6+p5m
     -   up6)*F(109)+p1mup6*F(114))))
       Fa9(116) = mup3mup6*p3mup2*F(104)+mup2mup6*p4mup3*F(105)+mup2mu
     -   p6*p5mup3*F(105)+mup3mup6*p2mup2*s345*F(106)+mup3mup6*p1mup2
     -   *s345*F(107)+mup3mup6*p4mup2*s345*F(108)+mup3mup6*p5mup2*s34
     -   5*F(108)+2*(p1mup2*(p4mup3+p5mup3)*((p3mup6+p4mup6+p5mup6)*F
     -   (46)+p1mup6*F(47)+p2mup6*F(48))-p2mup2*(p4mup3+p5mup3)*((p3m
     -   up6+p4mup6+p5mup6)*F(51)-p1mup6*F(52)-p2mup6*F(53))+p2mup6*p
     -   3mup2*p4mup3*F(60)+p2mup6*p4mup2*p4mup3*F(60)+p2mup6*p4mup3*
     -   p5mup2*F(60)+p2mup6*p3mup2*p5mup3*F(60)+p2mup6*p4mup2*p5mup3
     -   *F(60)+p2mup6*p5mup2*p5mup3*F(60)+p3mup2*p3mup6*p4mup3*F(109
     -   )+p3mup6*p4mup2*p4mup3*F(109)+p3mup2*p4mup3*p4mup6*F(109)+p4
     -   mup2*p4mup3*p4mup6*F(109)+p3mup6*p4mup3*p5mup2*F(109)+p4mup3
     -   *p4mup6*p5mup2*F(109)+p3mup2*p3mup6*p5mup3*F(109)+p3mup6*p4m
     -   up2*p5mup3*F(109)+p3mup2*p4mup6*p5mup3*F(109)+p4mup2*p4mup6*
     -   p5mup3*F(109)+p3mup6*p5mup2*p5mup3*F(109)+p4mup6*p5mup2*p5mu
     -   p3*F(109)+p3mup2*p4mup3*p5mup6*F(109)+p4mup2*p4mup3*p5mup6*F
     -   (109)+p4mup3*p5mup2*p5mup6*F(109)+p3mup2*p5mup3*p5mup6*F(109
     -   )+p4mup2*p5mup3*p5mup6*F(109)+p5mup2*p5mup3*p5mup6*F(109)+p1
     -   mup6*(p3mup2+p4mup2+p5mup2)*(p4mup3+p5mup3)*F(114))
       Fa9(117) = mup3mup6*(-(p1mup2*F(65))-p2mup2*F(66)+(p3mup2+p4mup
     -   2+p5mup2)*F(111))
       Fa9(118) = -(mup3mup6*p4mup2*F(37))+p4mup3*(mup2mup6*F(105)+2*(
     -   p1mup2*((p3mup6+p4mup6+p5mup6)*F(46)+p1mup6*F(47)+p2mup6*F(4
     -   8))-p2mup2*((p3mup6+p4mup6+p5mup6)*F(51)-p1mup6*F(52)-p2mup6
     -   *F(53))+(p3mup2+p4mup2+p5mup2)*(p2mup6*F(60)+(p3mup6+p4mup6+
     -   p5mup6)*F(109)+p1mup6*F(114))))
       Fa9(119) = p1mup6*s345*F(33)+(p4mup6+p5mup6)*s345*F(34)+p2mup6*
     -   s345*F(35)+p3mup6*F(115)
       Fa9(120) = -2*((p3mup6+p4mup6+p5mup6)*F(40)-p2mup6*F(63))+p1mup
     -   6*F(116)
       Fa9(121) = -4*(p1mup2*((p3mup6+p4mup6+p5mup6)*F(46)+p1mup6*F(47
     -   )+p2mup6*F(48))-p2mup2*((p3mup6+p4mup6+p5mup6)*F(51)-p1mup6*
     -   F(52)-p2mup6*F(53)))+(p3mup2+p4mup2+p5mup2)*((p3mup6+p4mup6+
     -   p5mup6)*F(58)+p1mup6*F(59)-4*p2mup6*F(60))+2*mup2mup6*F(117)
       Fa9(122) = mup4mup6*s45*F(118)+p5mup4*((p1mup6+p2mup6)*F(120)+4
     -   *((p4mup6+p5mup6)*F(119)-p3mup6*F(121)))
       Fa9(123) = p4mup6*(p3mup3*F(124)-4*(p5mup3*F(126)-(p1mup3+p2mup
     -   3)*F(127)))+p4mup3*((p1mup6+p2mup6)*F(128)+4*(p4mup6*F(123)+
     -   p3mup6*F(125)+p5mup6*F(129)))+s45*(mup3mup6*F(122)-4*(p1mup6
     -   *p2mup3*F(130)+p2mup3*p2mup6*F(130)+p1mup6*p3mup3*F(131)+p2m
     -   up6*p3mup3*F(131)+p3mup3*p3mup6*F(132)+p2mup3*p3mup6*F(133)+
     -   p3mup6*p5mup3*F(134)+p1mup6*p5mup3*F(135)+p2mup6*p5mup3*F(13
     -   5)-p2mup3*p5mup6*F(136)+p1mup3*((p1mup6+p2mup6)*F(130)+p3mup
     -   6*F(133)-p5mup6*F(136))-p5mup3*p5mup6*F(137)-p3mup3*p5mup6*F
     -   (138)))
       Fa9(124) = (p4mup6*F(139))/s45+2*p3mup6*F(140)
       Fa9(125) = p4mup2*((p1mup6+p2mup6)*F(142)-4*((p4mup6+p5mup6)*F(
     -   119)-p3mup6*F(143)))+4*p3mup2*(-(p4mup6*F(141))+s45*((p1mup6
     -   +p2mup6)*F(144)+p5mup6*F(145)+p3mup6*F(146)))
       Fa9(126) = -6*mup4mup6*s45*F(148)+4*p5mup4*((p4mup6+p5mup6)*F(1
     -   49)+(p1mup6+p2mup6)*F(150)+p3mup6*F(151))
       Return
       End
