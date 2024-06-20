c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineNoAbe1g8(p1mup2,p1mup3,p1mup4,
     -   p1mup6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mu
     -   p6,p4mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p
     -   6mup2,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mu
     -   p4,mup3mup6,mup4mup6,Fa8)
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
       Complex*16 Fa8(99:112)
       COMMON/test_FaFunctionsNoAbe1g/Fa
       Fa8(99) = 8*p3mup2*p5mup4*((p4mup6+p5mup6)*F(72)+p2mup6*F(74)+p
     -   3mup6*F(75))+mup4mup6*s45*(4*(p1mup2*F(65)+p2mup2*F(66))+(p3
     -   mup2+p5mup2)*F(90)+p4mup2*s345*F(92))+4*p5mup4*(mup2mup6*s34
     -   5*F(91)-2*(-(p1mup2*(p3mup6*F(76)-s345*(p1mup6*F(79)+p2mup6*
     -   F(81))+(p4mup6+p5mup6)*s345*F(83)))+p2mup6*p4mup2*s345*F(84)
     -   +p2mup6*p5mup2*s345*F(84)+p1mup6*p4mup2*s345*F(85)+p1mup6*p5
     -   mup2*s345*F(85)-p4mup2*p4mup6*s345*F(86)-p4mup6*p5mup2*s345*
     -   F(86)-p4mup2*p5mup6*s345*F(86)-p5mup2*p5mup6*s345*F(86)-p2mu
     -   p2*(p3mup6*F(77)-s345*(p1mup6*F(80)+p2mup6*F(82))+(p4mup6+p5
     -   mup6)*s345*F(87))-p3mup6*p4mup2*F(93)-p3mup6*p5mup2*F(93)-p1
     -   mup6*p3mup2*F(94)))
       Fa8(100) = 4*mup2mup6*p4mup3*s345*F(91)+mup3mup6*p3mup2*s45*F(9
     -   5)+mup3mup6*p5mup2*s45*F(95)-mup3mup6*p1mup2*s45*F(96)-mup3m
     -   up6*p2mup2*s45*F(97)+8*(p3mup2*p4mup6*p5mup3*F(64)+p4mup6*p5
     -   mup2*p5mup3*F(64)+p1mup2*p4mup6*p5mup3*F(65)+p2mup2*p4mup6*p
     -   5mup3*F(66)+p3mup2*p4mup3*p5mup6*F(72)+p2mup6*p3mup2*p4mup3*
     -   F(74)+p3mup2*p4mup3*p4mup6*F(75)+p1mup2*p4mup3*p4mup6*F(76)+
     -   p2mup2*p4mup3*p4mup6*F(77)-p1mup2*p1mup6*p4mup3*s345*F(79)-p
     -   1mup6*p2mup2*p4mup3*s345*F(80)-p1mup2*p2mup6*p4mup3*s345*F(8
     -   1)-p2mup2*p2mup6*p4mup3*s345*F(82)+p1mup2*p4mup3*p5mup6*s345
     -   *F(83)-p2mup6*p4mup3*p5mup2*s345*F(84)-p1mup6*p4mup3*p5mup2*
     -   s345*F(85)+p4mup3*p5mup2*p5mup6*s345*F(86)+p2mup2*p4mup3*p5m
     -   up6*s345*F(87)+p4mup3*p4mup6*p5mup2*F(93)+p3mup6*(p4mup2*(p5
     -   mup3*F(40)+p4mup3*F(75))+p4mup3*(p3mup2*F(75)+p1mup2*F(76)+p
     -   2mup2*F(77)+p5mup2*F(93)))+p1mup6*p3mup2*p4mup3*F(94)+p4mup2
     -   *(p4mup3*(p5mup6*F(72)+p2mup6*F(74)+p4mup6*F(75)+p1mup6*F(94
     -   ))+p5mup3*(p5mup6*F(40)-p2mup6*F(63)+p1mup6*F(89)-p4mup6*F(9
     -   8))))+mup3mup6*p4mup2*s345*s45*F(99)
       Fa8(101) = -(p2mup2*F(100))+p4mup2*s345*F(101)+(p3mup2+p5mup2)*
     -   F(102)-p1mup2*F(103)
       Fa8(102) = 16*p1mup2*p4mup3*(p3mup6+p4mup6)*p5mup4*F(46)+8*mup4
     -   mup6*(p3mup2*p4mup3*F(104)+p4mup2*(p5mup3*F(37)+p4mup3*F(104
     -   ))+p4mup3*s345*(p2mup2*F(106)+p1mup2*F(107)+p5mup2*F(108)))+
     -   p5mup4*(8*(mup3mup6*p3mup2*F(104)+mup2mup6*p4mup3*F(105)+mup
     -   2mup6*p5mup3*F(105)+p2mup2*(-2*(p4mup3+p5mup3)*((p3mup6+p4mu
     -   p6+p5mup6)*F(51)-p1mup6*F(52)-p2mup6*F(53))+mup3mup6*s345*F(
     -   106))+p1mup2*(2*(p3mup6*p5mup3*F(46)+p4mup6*p5mup3*F(46)+(p4
     -   mup3+p5mup3)*(p5mup6*F(46)+p1mup6*F(47)+p2mup6*F(48)))+mup3m
     -   up6*s345*F(107))+mup3mup6*p4mup2*s345*F(108)+mup3mup6*p5mup2
     -   *s345*F(108))+16*(p2mup6*p3mup2*p4mup3*F(60)+p2mup6*p4mup2*p
     -   4mup3*F(60)+p2mup6*p4mup3*p5mup2*F(60)+p2mup6*p3mup2*p5mup3*
     -   F(60)+p2mup6*p4mup2*p5mup3*F(60)+p2mup6*p5mup2*p5mup3*F(60)+
     -   p3mup2*p3mup6*p4mup3*F(109)+p3mup6*p4mup2*p4mup3*F(109)+p3mu
     -   p2*p4mup3*p4mup6*F(109)+p4mup2*p4mup3*p4mup6*F(109)+p3mup6*p
     -   4mup3*p5mup2*F(109)+p4mup3*p4mup6*p5mup2*F(109)+p3mup2*p3mup
     -   6*p5mup3*F(109)+p3mup6*p4mup2*p5mup3*F(109)+p3mup2*p4mup6*p5
     -   mup3*F(109)+p4mup2*p4mup6*p5mup3*F(109)+p3mup6*p5mup2*p5mup3
     -   *F(109)+p4mup6*p5mup2*p5mup3*F(109)+p3mup2*p4mup3*p5mup6*F(1
     -   09)+p4mup2*p4mup3*p5mup6*F(109)+p4mup3*p5mup2*p5mup6*F(109)+
     -   p3mup2*p5mup3*p5mup6*F(109)+p4mup2*p5mup3*p5mup6*F(109)+p5mu
     -   p2*p5mup3*p5mup6*F(109))+p1mup6*(p3mup2+p4mup2+p5mup2)*(p4mu
     -   p3+p5mup3)*F(110))
       Fa8(103) = p3mup2*p4mup3*F(104)+p4mup2*(p5mup3*F(37)+p4mup3*F(1
     -   04))+p4mup3*s345*(p2mup2*F(106)+p1mup2*F(107)+p5mup2*F(108))
       Fa8(104) = -((mup4mup6*p4mup3+mup3mup6*p5mup4)*(p1mup2*F(65)+p2
     -   mup2*F(66)-(p3mup2+p4mup2+p5mup2)*F(111)))
       Fa8(105) = p4mup3*((p3mup2+p4mup2+p5mup2)*F(64)+p1mup2*F(65)+p2
     -   mup2*F(66))
       Fa8(106) = -(mup4mup6*(p4mup3+p5mup3)*(p1mup2*F(65)+p2mup2*F(66
     -   )-(p3mup2+p4mup2+p5mup2)*F(111)))
       Fa8(107) = (p4mup3+p5mup3)*((p3mup2+p4mup2+p5mup2)*F(64)+p1mup2
     -   *F(65)+p2mup2*F(66))
       Fa8(108) = (p4mup3+p5mup3)*p5mup4*((p3mup2+p4mup2+p5mup2)*F(64)
     -   +p1mup2*F(65)+p2mup2*F(66))
       Fa8(109) = mup4mup6*(p3mup2*F(104)+s345*(p2mup2*F(106)+p1mup2*F
     -   (107)+(p4mup2+p5mup2)*F(108)))
       Fa8(110) = 2*p3mup2*F(104)+s345*(2*(p4mup2+p5mup2)*F(108)+p2mup
     -   2*F(112)+p1mup2*F(113))
       Fa8(111) = mup4mup6*((p3mup2+p4mup2+p5mup2)*F(64)+p1mup2*F(65)+
     -   p2mup2*F(66))
       Fa8(112) = -(p1mup2*F(65))-p2mup2*F(66)+(p3mup2+p4mup2+p5mup2)*
     -   F(111)
       Return
       End
