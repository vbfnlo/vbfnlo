c************** Calling the Fa functions*************************
       subroutine NoAbe4vertHexLineFa5(p1mup2,p1mup3,p1mup5,p1mup6,p2m
     -   up2,p2mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,
     -   p4mup3,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mu
     -   p3,p6mup5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup
     -   6,mup5mup6,Fa5)
       IMPLICIT NONE
      Complex*16   p1mup2, p1mup3, p1mup5, p1mup6, p2mup2, p2mup3, 
     -          p2mup5, p2mup6, p3mup2, p3mup3, p3mup5, p3mup6, 
     -          p4mup2, p4mup3, p4mup5, p4mup6, p5mup2, p5mup3, 
     -          p5mup5, p5mup6, p6mup2, p6mup3, p6mup5, p6mup6
       Complex*16   mup2mup3, mup2mup5, mup2mup6, mup3mup5, mup3mup6, 
     -          mup5mup6
        common/NoAbe4vertHexLineFhlFunctions/F
       COMMON/NoAbe4vertHexLineInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6s
     -   q,s12,s23,s34,s45,s56,s16,s123,s234,s345
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
      Complex*16 Fa(199),F(140)
      Real*8 P(39) 
        COMMON/NoAbe4vertHexLinePFunctions/P
       Complex*16 Fa5(77:95)
       COMMON/NoAbe4vertHexLineFaFunctions/Fa
       Fa5(77) = -4*(p3mup5*((p3mup2+p4mup2)*F(96)+p2mup2*F(100)+p1mup
     -   2*F(101))+p1mup5*p3mup2*F(103))+2*mup2mup5*s34*F(126)
       Fa5(78) = -(p3mup2*p4mup5*F(126))+p3mup5*(p4mup2*F(96)+p2mup2*F
     -   (100)+p1mup2*F(101)+p3mup2*F(130))
       Fa5(79) = -2*mup2mup5*F(92)+4*(p1mup5*((p3mup2+p4mup2)*F(95)+p2
     -   mup2*F(98)+p1mup2*F(99))+p2mup5*((p3mup2+p4mup2)*F(96)+p2mup
     -   2*F(100)+p1mup2*F(101)))
       Fa5(80) = (p3mup5+p4mup5)*(p1mup2*F(95)+p2mup2*F(96)+(p3mup2+p4
     -   mup2)*F(128))
       Fa5(81) = mup3mup6*p1mup5*F(103)+mup3mup5*(p3mup6+p4mup6)*F(126
     -   )
       Fa5(82) = mup3mup5*p1mup6*F(103)+mup3mup6*(p3mup5+p4mup5)*F(126
     -   )
       Fa5(83) = mup2mup5*(p3mup6+p4mup6)
       Fa5(84) = mup5mup6*(p1mup2*F(98)+p2mup2*F(100)+(p3mup2+p4mup2)*
     -   F(130))
       Fa5(85) = mup5mup6*((p3mup2+p4mup2)*F(96)+p2mup2*F(100)+p1mup2*
     -   F(101))-mup2mup6*(p3mup5+p4mup5)*F(126)
       Fa5(86) = mup3mup6*((p3mup2+p4mup2)*F(96)+p2mup2*F(100)+p1mup2*
     -   F(101))
       Fa5(87) = (p3mup2+p4mup2)*F(96)+p2mup2*F(100)+p1mup2*F(101)
       Fa5(88) = mup3mup5*((p3mup2+p4mup2)*F(96)+p2mup2*F(100)+p1mup2*
     -   F(101))
       Fa5(89) = p3mup6+p4mup6
       Fa5(90) = p3mup5+p4mup5
       Fa5(91) = mup5mup6*((p1mup3+p2mup3)*F(104)+p3mup3*F(105)+p4mup3
     -   *F(106))
       Fa5(92) = mup5mup6*((p1mup3+p2mup3)*F(104)+p4mup3*F(106)+p3mup3
     -   *F(107))
       Fa5(93) = (p1mup2+p2mup2)*p4mup6*((p1mup3+p2mup3)*F(106)+p3mup3
     -   *F(108)+p4mup3*F(129))
       Fa5(94) = -(mup3mup6*p2mup2*F(109))-4*mup3mup6*p3mup2*s12*F(110
     -   )+8*(p1mup6*p2mup2*p2mup3*F(104)+p2mup2*p2mup3*p2mup6*F(104)
     -   +p2mup2*p2mup3*p3mup6*F(105)+p1mup3*p2mup2*((p1mup6+p2mup6)*
     -   F(104)+p3mup6*F(105))+p1mup6*p2mup2*p4mup3*F(106)+p2mup2*p2m
     -   up6*p4mup3*F(106)+p1mup6*p2mup2*p3mup3*F(107)+p2mup2*p2mup6*
     -   p3mup3*F(107)+p2mup2*p3mup6*p4mup3*F(108)+p2mup2*p3mup3*p3mu
     -   p6*F(111))+p1mup2*(-(mup3mup6*F(109))+8*(p1mup6*p2mup3*F(104
     -   )+p2mup3*p2mup6*F(104)+p1mup3*(p1mup6+p2mup6)*F(104)+p1mup3*
     -   p3mup6*F(105)+p2mup3*p3mup6*F(105)+p1mup6*p4mup3*F(106)+p2mu
     -   p6*p4mup3*F(106)+p1mup6*p3mup3*F(107)+p2mup6*p3mup3*F(107)+p
     -   3mup6*p4mup3*F(108)+p3mup3*p3mup6*F(111)))
       Fa5(95) = (p1mup3+p2mup3)*F(104)+p4mup3*F(106)+p3mup3*F(107)
       Return
       End
