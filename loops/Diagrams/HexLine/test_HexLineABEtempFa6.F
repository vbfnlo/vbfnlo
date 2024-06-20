c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineAbe6(p1mup2,p1mup3,p1mup4,p1mu
     -   p5,p2mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p
     -   4mup2,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup
     -   2,p6mup3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,m
     -   up3mup5,mup4mup5,Fa6)
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
       Complex*16 Fa6(56:66)
       COMMON/test_FaFunctionsAbe/Fa
       Fa6(56) = (p1mup4+p2mup4+p3mup4+p4mup4)*(p3mup2*F(94)+p1mup2*F(
     -   95)+p2mup2*F(96))
       Fa6(57) = p4mup2*F(76)+s123*(p3mup2*F(78)-p1mup2*F(98)+p2mup2*F
     -   (99))
       Fa6(58) = 4*(p3mup2*F(62)+p2mup2*F(72))+p1mup2*F(102)
       Fa6(59) = -4*(p2mup3*F(60)+p3mup3*F(64))+p1mup3*F(103)
       Fa6(60) = p3mup2*F(94)+p1mup2*F(95)+p2mup2*F(96)
       Fa6(61) = (p1mup4+p2mup4+p3mup4+p4mup4)*(p1mup3*F(48)-p2mup3*F(
     -   49)+p3mup3*F(50))
       Fa6(62) = (p1mup4+p2mup4+p3mup4+p4mup4)*(mup2mup3*F(104)+2*(p1m
     -   up2*(p1mup3*F(105)+p2mup3*F(106)+p3mup3*F(107))+p1mup3*p2mup
     -   2*F(108)+p2mup2*p2mup3*F(109)+p2mup2*p3mup3*F(110)-p1mup3*p3
     -   mup2*F(111)-p3mup2*p3mup3*F(112)-p2mup3*p3mup2*F(113)))
       Fa6(63) = p1mup3*F(48)-p2mup3*F(49)+p3mup3*F(50)
       Fa6(64) = 4*mup2mup3*F(104)+8*(p1mup2*(p1mup3*F(105)+p2mup3*F(1
     -   06)+p3mup3*F(107))+p1mup3*p2mup2*F(108)+p2mup2*p2mup3*F(109)
     -   +p2mup2*p3mup3*F(110)-p1mup3*p3mup2*F(111)-p3mup2*p3mup3*F(1
     -   12)-p2mup3*p3mup2*F(113))
       Fa6(65) = p4mup2*s12*F(114)+4*(p1mup2+p2mup2)*F(115)+p3mup2*s12
     -   *F(116)
       Fa6(66) = p3mup3*F(117)+p4mup3*F(118)+4*(p1mup3+p2mup3)*F(119)
       Return
       End
