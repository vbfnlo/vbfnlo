c************** Calling the Fa functions*************************
       subroutine NoAbe4vertHexLineFa4(p1mup2,p1mup3,p1mup5,p1mup6,p2m
     -   up2,p2mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,
     -   p4mup3,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mu
     -   p3,p6mup5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup
     -   6,mup5mup6,Fa4)
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
       Complex*16 Fa4(58:76)
       COMMON/NoAbe4vertHexLineFaFunctions/Fa
       Fa4(58) = -(mup2mup5*mup3mup6*F(92))+2*(mup3mup6*(p1mup5*((p3mu
     -   p2+p4mup2)*F(95)+p2mup2*F(98)+p1mup2*F(99))+p2mup5*((p3mup2+
     -   p4mup2)*F(96)+p2mup2*F(100)+p1mup2*F(101)))+mup3mup5*(p3mup6
     -   +p4mup6)*(p1mup2*F(95)+p2mup2*F(96)+(p3mup2+p4mup2)*F(128)))
       Fa4(59) = -(mup2mup6*mup3mup5*F(92))+2*(mup3mup5*(p1mup6*((p3mu
     -   p2+p4mup2)*F(95)+p2mup2*F(98)+p1mup2*F(99))+p2mup6*((p3mup2+
     -   p4mup2)*F(96)+p2mup2*F(100)+p1mup2*F(101)))+mup3mup6*(p3mup5
     -   +p4mup5)*(p1mup2*F(95)+p2mup2*F(96)+(p3mup2+p4mup2)*F(128)))
       Fa4(60) = -((p3mup2*F(94))/s34)+p1mup2*F(95)+p2mup2*F(96)+p4mup
     -   2*F(128)
       Fa4(61) = p1mup2*F(95)+p2mup2*F(96)+(p3mup2*F(97))/s34+p4mup2*F
     -   (128)
       Fa4(62) = p4mup3*((p3mup2+p4mup2)*F(96)+p2mup2*F(100)+p1mup2*F(
     -   101))
       Fa4(63) = p1mup5*p3mup6*F(103)+p3mup5*(p3mup6+p4mup6)*F(126)
       Fa4(64) = -4*p1mup6*p3mup5*F(103)+2*(-2*p3mup6*(p3mup5+p4mup5)+
     -   mup5mup6*s34)*F(126)
       Fa4(65) = p2mup5*(p3mup6+p4mup6)
       Fa4(66) = -2*mup5mup6*F(92)+4*p2mup6*(p3mup5+p4mup5)*F(126)
       Fa4(67) = -2*p4mup3*(p3mup6+p4mup6)+mup3mup6*s34
       Fa4(68) = (-4*p1mup6*p4mup3*F(103))/s34-2*mup3mup6*F(126)
       Fa4(69) = -((p3mup6*F(92))/s34)+p2mup6*F(126)
       Fa4(70) = -(p3mup2*p4mup6*F(126))+p3mup6*(p4mup2*F(96)+p2mup2*F
     -   (100)+p1mup2*F(101)+p3mup2*F(130))
       Fa4(71) = -4*(p3mup6*((p3mup2+p4mup2)*F(96)+p2mup2*F(100)+p1mup
     -   2*F(101))+p1mup6*p3mup2*F(103))+2*mup2mup6*s34*F(126)
       Fa4(72) = (p3mup6+p4mup6)*(p1mup2*F(95)+p2mup2*F(96)+(p3mup2+p4
     -   mup2)*F(128))
       Fa4(73) = -2*mup2mup6*F(92)+4*(p1mup6*((p3mup2+p4mup2)*F(95)+p2
     -   mup2*F(98)+p1mup2*F(99))+p2mup6*((p3mup2+p4mup2)*F(96)+p2mup
     -   2*F(100)+p1mup2*F(101)))
       Fa4(74) = -((p3mup5*F(92))/s34)+p2mup5*F(126)
       Fa4(75) = (-4*p1mup5*p4mup3*F(103))/s34-2*mup3mup5*F(126)
       Fa4(76) = -2*p4mup3*(p3mup5+p4mup5)+mup3mup5*s34
       Return
       End
