c************** Calling the Fa functions*************************
       subroutine NoAbe4vertHexLineFa9(p1mup2,p1mup3,p1mup5,p1mup6,p2m
     -   up2,p2mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,
     -   p4mup3,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mu
     -   p3,p6mup5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup
     -   6,mup5mup6,Fa9)
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
       Complex*16 Fa9(153:171)
       COMMON/NoAbe4vertHexLineFaFunctions/Fa
       Fa9(153) = -2*(p3mup2+p4mup2)*(p2mup5+p3mup5+p4mup5)+mup2mup5*s
     -   234
       Fa9(154) = mup3mup6*p1mup5*F(120)+mup3mup5*(p2mup6+p3mup6+p4mup
     -   6)*F(123)
       Fa9(155) = mup3mup5*p1mup6*F(120)+mup3mup6*(p2mup5+p3mup5+p4mup
     -   5)*F(123)
       Fa9(156) = mup2mup6*p1mup5*F(120)+mup2mup5*(p2mup6+p3mup6+p4mup
     -   6)*F(123)
       Fa9(157) = mup2mup5*p1mup6*F(120)+mup2mup6*(p2mup5+p3mup5+p4mup
     -   5)*F(123)
       Fa9(158) = p2mup6+p3mup6+p4mup6
       Fa9(159) = p2mup5+p3mup5+p4mup5
       Fa9(160) = mup3mup6*p3mup5-mup3mup5*p3mup6
       Fa9(161) = mup3mup6*p3mup5-mup3mup5*p3mup6+mup5mup6*p4mup3
       Fa9(162) = -4*(mup3mup6*p3mup2*s12*F(121)-(p1mup2+p2mup2)*(2*p4
     -   mup3*(p3mup6+p4mup6)-mup3mup6*s34)*F(124))
       Fa9(163) = (2*(p1mup2+p2mup2)*(p1mup6+p2mup6)*p4mup3+mup3mup6*p
     -   3mup2*s12)*F(121)+mup3mup6*(p1mup2+p2mup2)*s34*F(124)
       Fa9(164) = (2*(p1mup2+p2mup2)*(p1mup5+p2mup5)*p4mup3+mup3mup5*p
     -   3mup2*s12)*F(121)+mup3mup5*(p1mup2+p2mup2)*s34*F(124)
       Fa9(165) = -4*(mup3mup5*p3mup2*s12*F(121)-(p1mup2+p2mup2)*(2*p4
     -   mup3*(p3mup5+p4mup5)-mup3mup5*s34)*F(124))
       Fa9(166) = (p1mup2+p2mup2)*((p1mup5+p2mup5)*p3mup6*F(121)+p3mup
     -   5*(p3mup6+p4mup6)*F(124))
       Fa9(167) = mup5mup6*(p3mup2*s12*F(121)+(p1mup2+p2mup2)*s34*F(12
     -   4))
       Fa9(168) = 4*((2*(p1mup2+p2mup2)*(p1mup6+p2mup6)*p3mup5-mup5mup
     -   6*p3mup2*s12)*F(121)+(p1mup2+p2mup2)*(2*p3mup6*(p3mup5+p4mup
     -   5)-mup5mup6*s34)*F(124))
       Fa9(169) = p3mup2*(mup3mup6*(p1mup5+p2mup5)*F(121)+mup3mup5*(p3
     -   mup6+p4mup6)*F(124))
       Fa9(170) = p3mup2*(mup3mup5*(p1mup6+p2mup6)*F(121)+mup3mup6*(p3
     -   mup5+p4mup5)*F(124))
       Fa9(171) = (p1mup2+p2mup2)*(mup3mup6*(p1mup5+p2mup5)*F(121)+mup
     -   3mup5*(p3mup6+p4mup6)*F(124))
       Return
       End
