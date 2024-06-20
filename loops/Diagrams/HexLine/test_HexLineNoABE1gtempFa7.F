c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineNoAbe1g7(p1mup2,p1mup3,p1mup4,
     -   p1mup6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mu
     -   p6,p4mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p
     -   6mup2,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mu
     -   p4,mup3mup6,mup4mup6,Fa7)
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
       Complex*16 Fa7(85:98)
       COMMON/test_FaFunctionsNoAbe1g/Fa
       Fa7(85) = (p1mup2+p2mup2)*(p1mup3+p2mup3+p3mup3)*(p1mup4+p2mup4
     -   +p3mup4+p4mup4)
       Fa7(86) = (p1mup3+p2mup3+p3mup3)*(p1mup4+p2mup4+p3mup4+p4mup4)
       Fa7(87) = p3mup2*(p1mup4+p2mup4+p3mup4+p4mup4)
       Fa7(88) = (p1mup2+p2mup2)*(p1mup4+p2mup4+p3mup4+p4mup4)
       Fa7(89) = 2*mup4mup6*s45*F(31)+4*p5mup4*(-(p3mup6*F(32))+s345*(
     -   p1mup6*F(33)+(p4mup6+p5mup6)*F(34)+p2mup6*F(35)))
       Fa7(90) = mup3mup6*s45*F(36)-4*(p3mup6*p4mup3*F(32)-p4mup3*(-(p
     -   4mup6*F(32))+p1mup6*s345*F(33)+p5mup6*s345*F(34)+p2mup6*s345
     -   *F(35))+p4mup6*p5mup3*F(37))
       Fa7(91) = p2mup6*F(38)+p4mup6*s345*F(39)-2*(p3mup6+p5mup6)*F(40
     -   )+p1mup6*F(41)
       Fa7(92) = mup2mup6*s45*F(42)-p3mup2*p4mup6*F(43)+p4mup2*p4mup6*
     -   s345*s45*F(54)+p4mup6*p5mup2*s345*s45*F(55)+p2mup6*p4mup2*s3
     -   45*s45*F(56)+p3mup2*p3mup6*s45*F(58)+p3mup6*p5mup2*s45*F(58)
     -   +p3mup2*p5mup6*s45*F(58)+p5mup2*p5mup6*s45*F(58)+p1mup6*p3mu
     -   p2*s45*F(59)+p1mup6*p5mup2*s45*F(59)+4*(p3mup6*p4mup2*F(44)-
     -   p1mup2*s45*(p4mup6*F(45)+(p3mup6+p5mup6)*F(46)+p1mup6*F(47)+
     -   p2mup6*F(48))-p1mup6*p4mup2*F(49)+p2mup2*s45*(p4mup6*s345*F(
     -   50)+(p3mup6+p5mup6)*F(51))-p1mup6*p2mup2*s45*F(52)-p2mup2*p2
     -   mup6*s45*F(53)-p4mup2*p5mup6*F(57)-p2mup6*(p3mup2+p5mup2)*s4
     -   5*F(60))
       Fa7(93) = p5mup4*(p1mup6*F(62)-4*((p3mup6+p4mup6+p5mup6)*F(40)-
     -   p2mup6*F(63)))
       Fa7(94) = p4mup3*(p1mup6*F(62)-4*((p3mup6+p4mup6+p5mup6)*F(40)-
     -   p2mup6*F(63)))
       Fa7(95) = 4*p4mup6*((p3mup2+p5mup2)*F(64)+p1mup2*F(65)+p2mup2*F
     -   (66))+p4mup2*(4*((p3mup6+p5mup6)*F(40)-p2mup6*F(63)+p4mup6*F
     -   (67))+p1mup6*F(68))
       Fa7(96) = -((p4mup3+p5mup3)*(-(p1mup6*F(62))+4*((p3mup6+p4mup6+
     -   p5mup6)*F(40)-p2mup6*F(63))))+mup3mup6*s345*F(70)
       Fa7(97) = 2*mup2mup6*s345*F(71)+p3mup2*(p1mup6*F(73)+4*((p4mup6
     -   +p5mup6)*F(72)+p2mup6*F(74)+p3mup6*F(75)))+4*(p1mup2*(p3mup6
     -   *F(76)-s345*(p1mup6*F(79)+p2mup6*F(81))+(p4mup6+p5mup6)*s345
     -   *F(83))+(p4mup2+p5mup2)*(p3mup6*F(78)-s345*(p2mup6*F(84)+p1m
     -   up6*F(85))+(p4mup6+p5mup6)*s345*F(86))+p2mup2*(p3mup6*F(77)-
     -   s345*(p1mup6*F(80)+p2mup6*F(82))+(p4mup6+p5mup6)*s345*F(87))
     -   )
       Fa7(98) = 8*p3mup6*(p4mup3+p5mup3)*p5mup4*F(40)+4*p5mup4*(-(mup
     -   3mup6*s345*F(88))+2*p5mup3*((p4mup6+p5mup6)*F(40)-p2mup6*F(6
     -   3)+p1mup6*F(89)))+p4mup3*(-4*mup4mup6*s345*F(88)+8*p5mup4*((
     -   p4mup6+p5mup6)*F(40)-p2mup6*F(63)+p1mup6*F(89)))
       Return
       End
