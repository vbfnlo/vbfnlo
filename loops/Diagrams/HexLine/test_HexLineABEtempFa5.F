c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineAbe5(p1mup2,p1mup3,p1mup4,p1mu
     -   p5,p2mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p
     -   4mup2,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup
     -   2,p6mup3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,m
     -   up3mup5,mup4mup5,Fa5)
       IMPLICIT NONE
      Complex*16   p1mup2, p1mup3, p1mup4, p1mup5, p2mup2, p2mup3, 
     -          p2mup4, p2mup5, p3mup2, p3mup3, p3mup4, p3mup5, 
     -          p4mup2, p4mup3, p4mup4, p4mup5, p5mup2, p5mup3, 
     -          p5mup4, p5mup5, p6mup2, p6mup3, p6mup4, p6mup5
       Complex*16   mup2mup3, mup2mup4, mup2mup5, mup3mup4, mup3mup5, 
     -          mup4mup5
        common/test_FhlFunctionsAbe/F
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
      Complex*16 Fa(445),F(5930)
      Real*8 P(11674) 
        COMMON/test_PFunctionsAbe/P
       Complex*16 Fa5(45:55)
       COMMON/test_FaFunctionsAbe/Fa
       Fa5(45) = -((p1mup3+p2mup3+p3mup3)*(p6mup5*F(35)-(p1mup5+p2mup5
     -   +p3mup5+p4mup5)*F(36)))
       Fa5(46) = -((p1mup2+p2mup2)*(p6mup5*F(35)-(p1mup5+p2mup5+p3mup5
     -   +p4mup5)*F(36)))
       Fa5(47) = p4mup2*F(37)+s123*(p2mup2*F(38)+p3mup2*F(39)-p1mup2*F
     -   (40))
       Fa5(48) = p4mup3*F(42)+4*s123*(p2mup3*F(41)+p1mup3*F(43)+p3mup3
     -   *F(44))
       Fa5(49) = -4*mup2mup3*s123*F(51)+8*(p4mup2*(p1mup3*F(48)-p2mup3
     -   *F(49)+p3mup3*F(50))+p1mup2*(-(p4mup3*F(45))+s123*(p1mup3*F(
     -   52)+p2mup3*F(54)+p3mup3*F(56)))+p3mup2*(p4mup3*F(47)+p1mup3*
     -   s123*F(56)+p2mup3*s123*F(57)+p3mup3*s123*F(58))+p2mup2*(-(p4
     -   mup3*F(46))+p1mup3*s123*F(53)+p2mup3*s123*F(55)+p3mup3*s123*
     -   F(59)))
       Fa5(50) = p1mup2*p4mup3*F(61)+4*mup2mup3*s123*F(65)-8*(p3mup2*p
     -   4mup3*F(62)-p1mup3*(p4mup2*F(63)+p1mup2*s123*F(66)+p3mup2*s1
     -   23*F(69)+p2mup2*s123*F(70))+p2mup2*p4mup3*F(72)-p2mup3*(p4mu
     -   p2*F(60)+s123*(p1mup2*F(67)+p2mup2*F(71)+p3mup2*F(73)))-p3mu
     -   p3*(p4mup2*F(64)+p1mup2*s123*F(68)+p3mup2*s123*F(74)+p2mup2*
     -   s123*F(75)))
       Fa5(51) = p4mup3*F(76)+s123*(-(p1mup3*F(77))+p2mup3*F(78)+p3mup
     -   3*F(79))
       Fa5(52) = (p1mup4+p2mup4+p3mup4+p4mup4)*(mup2mup3*F(83)+16*(p2m
     -   up2*(p2mup3*F(81)-p1mup3*F(84)+p3mup3*F(85))+p3mup2*(p2mup3*
     -   F(82)-p1mup3*F(86)+p3mup3*F(87))+p1mup2*(p2mup3*F(80)-p1mup3
     -   *F(88)-p3mup3*F(89))))
       Fa5(53) = mup2mup3*F(90)+2*(p2mup2*(p2mup3*F(81)-p1mup3*F(84)+p
     -   3mup3*F(85))+p3mup2*(p2mup3*F(82)-p1mup3*F(86)+p3mup3*F(87))
     -   +p1mup2*(p2mup3*F(80)-p1mup3*F(88)-p3mup3*F(91)))
       Fa5(54) = (p1mup4+p2mup4+p3mup4+p4mup4)*(p3mup2*F(62)+p2mup2*F(
     -   72)+p1mup2*F(92))
       Fa5(55) = -((p1mup4+p2mup4+p3mup4+p4mup4)*(8*(p2mup3*F(60)+p3mu
     -   p3*F(64))-p1mup3*F(93)))
       Return
       End
