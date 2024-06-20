c************** Calling the Fa functions*************************
       subroutine NoAbe4vertHexLineFa8(p1mup2,p1mup3,p1mup5,p1mup6,p2m
     -   up2,p2mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,
     -   p4mup3,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mu
     -   p3,p6mup5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup
     -   6,mup5mup6,Fa8)
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
       Complex*16 Fa8(134:152)
       COMMON/NoAbe4vertHexLineFaFunctions/Fa
       Fa8(134) = p1mup5*(mup3mup6*p3mup2+mup2mup6*p4mup3)*F(120)+(mup
     -   3mup5*p3mup2+mup2mup5*p4mup3)*(p2mup6+p3mup6+p4mup6)*F(123)
       Fa8(135) = p1mup6*(mup3mup5*p3mup2+mup2mup5*p4mup3)*F(120)+(mup
     -   3mup6*p3mup2+mup2mup6*p4mup3)*(p2mup5+p3mup5+p4mup5)*F(123)
       Fa8(136) = 8*mup3mup6*p1mup5*(p3mup2+p4mup2)*F(120)+4*(2*mup3mu
     -   p5*(p3mup2+p4mup2)*(p2mup6+p3mup6+p4mup6)-mup2mup6*mup3mup5*
     -   s234+mup2mup5*mup3mup6*s234)*F(123)
       Fa8(137) = 8*mup3mup5*p1mup6*(p3mup2+p4mup2)*F(120)+4*(2*mup3mu
     -   p6*(p3mup2+p4mup2)*(p2mup5+p3mup5+p4mup5)+mup2mup6*mup3mup5*
     -   s234-mup2mup5*mup3mup6*s234)*F(123)
       Fa8(138) = p1mup5*p3mup6*F(120)+p3mup5*(p2mup6+p3mup6+p4mup6)*F
     -   (123)
       Fa8(139) = p1mup6*p3mup5*F(120)+p3mup6*(p2mup5+p3mup5+p4mup5)*F
     -   (123)
       Fa8(140) = p1mup5*p2mup6*F(120)+p2mup5*(p2mup6+p3mup6+p4mup6)*F
     -   (123)
       Fa8(141) = 4*p1mup6*p2mup5*F(120)+2*(2*p2mup6*(p2mup5+p3mup5+p4
     -   mup5)-mup5mup6*s234)*F(123)
       Fa8(142) = p4mup3*(p2mup6+p3mup6+p4mup6)
       Fa8(143) = p2mup6+p4mup6+p3mup6*s234*P(39)
       Fa8(144) = (p1mup6*F(120))/s234+(p3mup6*F(123))/s34
       Fa8(145) = p3mup2*(p2mup6+p3mup6+p4mup6)
       Fa8(146) = -2*(p3mup2+p4mup2)*(p2mup6+p3mup6+p4mup6)+mup2mup6*s
     -   234
       Fa8(147) = 2*p1mup6*(p3mup2+p4mup2)*F(120)+mup2mup6*s234*F(123)
       Fa8(148) = (p1mup5*F(120))/s234+(p3mup5*F(123))/s34
       Fa8(149) = p2mup5+p4mup5+p3mup5*s234*P(39)
       Fa8(150) = p4mup3*(p2mup5+p3mup5+p4mup5)
       Fa8(151) = p3mup2*(p2mup5+p3mup5+p4mup5)
       Fa8(152) = 2*p1mup5*(p3mup2+p4mup2)*F(120)+mup2mup5*s234*F(123)
       Return
       End
