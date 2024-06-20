c************** Calling the Fa functions*************************
       subroutine NoAbe4vertHexLineFa10(p1mup2,p1mup3,p1mup5,p1mup6,p2
     -   mup2,p2mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2
     -   ,p4mup3,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6m
     -   up3,p6mup5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mu
     -   p6,mup5mup6,Fa10)
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
       Complex*16 Fa10(172:199)
       COMMON/NoAbe4vertHexLineFaFunctions/Fa
       Fa10(172) = (p1mup2+p2mup2)*(mup3mup5*(p1mup6+p2mup6)*F(121)+mu
     -   p3mup6*(p3mup5+p4mup5)*F(124))
       Fa10(173) = p3mup2*s12*F(121)+(p1mup2+p2mup2)*s34*F(124)
       Fa10(174) = (p1mup5+p2mup5)*p3mup6*F(121)+p3mup5*(p3mup6+p4mup6
     -   )*F(124)
       Fa10(175) = -4*(p1mup6+p2mup6)*p3mup5*F(121)-2*(2*p3mup6*(p3mup
     -   5+p4mup5)-mup5mup6*s34)*F(124)
       Fa10(176) = 2*(p1mup6+p2mup6)*p4mup3*F(121)+mup3mup6*s34*F(124)
       Fa10(177) = p3mup2*(p3mup6+p4mup6)
       Fa10(178) = (p1mup2+p2mup2)*(p3mup6+p4mup6)
       Fa10(179) = 2*(p1mup5+p2mup5)*p4mup3*F(121)+mup3mup5*s34*F(124)
       Fa10(180) = p3mup2*(p3mup5+p4mup5)
       Fa10(181) = (p1mup2+p2mup2)*(p3mup5+p4mup5)
       Fa10(182) = mup3mup6*(p1mup5+p2mup5)*F(121)+mup3mup5*(p3mup6+p4
     -   mup6)*F(124)
       Fa10(183) = mup3mup5*(p1mup6+p2mup6)*F(121)+mup3mup6*(p3mup5+p4
     -   mup5)*F(124)
       Fa10(184) = (p1mup2+p2mup2)*(p1mup3+p2mup3+p3mup3)*p4mup6
       Fa10(185) = (p1mup2+p2mup2)*(p1mup3+p2mup3+p3mup3)*(p1mup6+p2mu
     -   p6+p3mup6)
       Fa10(186) = (p1mup2+p2mup2)*(p1mup3+p2mup3+p3mup3)*(p1mup5+p2mu
     -   p5+p3mup5)
       Fa10(187) = (p1mup2+p2mup2)*(p1mup3+p2mup3+p3mup3)*p4mup5
       Fa10(188) = mup5mup6*(p1mup2+p2mup2)
       Fa10(189) = p1mup2+p2mup2
       Fa10(190) = (p1mup3+p2mup3+p3mup3)*p4mup6
       Fa10(191) = (p1mup3+p2mup3+p3mup3)*(p1mup6+p2mup6+p3mup6)
       Fa10(192) = p1mup6+p2mup6+p3mup6
       Fa10(193) = p3mup2*(p1mup6+p2mup6+p3mup6)
       Fa10(194) = (p1mup2+p2mup2)*(p1mup6+p2mup6+p3mup6)
       Fa10(195) = p1mup5+p2mup5+p3mup5
       Fa10(196) = (p1mup3+p2mup3+p3mup3)*(p1mup5+p2mup5+p3mup5)
       Fa10(197) = (p1mup3+p2mup3+p3mup3)*p4mup5
       Fa10(198) = p3mup2*(p1mup5+p2mup5+p3mup5)
       Fa10(199) = (p1mup2+p2mup2)*(p1mup5+p2mup5+p3mup5)
       Return
       End
