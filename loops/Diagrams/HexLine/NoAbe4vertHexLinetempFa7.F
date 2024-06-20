c************** Calling the Fa functions*************************
       subroutine NoAbe4vertHexLineFa7(p1mup2,p1mup3,p1mup5,p1mup6,p2m
     -   up2,p2mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,
     -   p4mup3,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mu
     -   p3,p6mup5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup
     -   6,mup5mup6,Fa7)
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
       Complex*16 Fa7(115:133)
       COMMON/NoAbe4vertHexLineFaFunctions/Fa
       Fa7(115) = (p1mup2+p2mup2)*(p1mup6+p2mup6)
       Fa7(116) = 4*(p2mup3*p2mup5*F(104)+p1mup3*(p1mup5+p2mup5)*F(104
     -   )+p1mup3*p3mup5*F(105)+p2mup3*p3mup5*F(105)+p2mup5*p4mup3*F(
     -   106)+p2mup5*p3mup3*F(107)+p1mup5*(p2mup3*F(104)+p4mup3*F(106
     -   )+p3mup3*F(107))+p3mup5*p4mup3*F(108)+p3mup3*p3mup5*F(111))+
     -   mup3mup5*F(114)
       Fa7(117) = p4mup5*((p1mup3+p2mup3)*F(106)+p3mup3*F(108)+p4mup3*
     -   F(129))
       Fa7(118) = (p1mup5+p2mup5)*p3mup2
       Fa7(119) = (p1mup2+p2mup2)*(p1mup5+p2mup5)
       Fa7(120) = (p1mup2+p2mup2)*p4mup5
       Fa7(121) = mup5mup6*((p1mup3+p2mup3)*F(107)+p3mup3*F(111)+p4mup
     -   3*F(131))
       Fa7(122) = mup5mup6*((p1mup3+p2mup3)*F(105)+p4mup3*F(108)+p3mup
     -   3*F(111))-mup3mup6*p4mup5*F(127)
       Fa7(123) = (p1mup3+p2mup3)*F(105)+p4mup3*F(108)+p3mup3*F(111)
       Fa7(124) = p1mup6+p2mup6
       Fa7(125) = p1mup5+p2mup5
       Fa7(126) = p1mup5*(2*p2mup6*p4mup3+mup3mup6*s34)*F(120)+F(123)*
     -   (2*p2mup5*p4mup3*(p2mup6+p3mup6+p4mup6)+mup3mup6*p3mup5*s234
     -   +mup3mup5*s34*(p2mup6+p4mup6+p3mup6*s234*P(39)))
       Fa7(127) = p1mup6*(2*p2mup5*p4mup3+mup3mup5*s34)*F(120)+F(123)*
     -   (2*p2mup6*p4mup3*(p3mup5+p4mup5)+mup3mup5*p3mup6*s234-mup5mu
     -   p6*p4mup3*s234+mup3mup6*p4mup5*s34+p2mup5*(2*p2mup6*p4mup3+m
     -   up3mup6*s34)+mup3mup6*p3mup5*s234*s34*P(39))
       Fa7(128) = -2*(p3mup2+p4mup2)*p4mup3*(p2mup6+p3mup6+p4mup6)+(mu
     -   p3mup6*p3mup2+mup2mup6*p4mup3)*s234
       Fa7(129) = 2*p1mup6*(p3mup2+p4mup2)*p4mup3*F(120)+(mup3mup6*p3m
     -   up2+mup2mup6*p4mup3)*s234*F(123)
       Fa7(130) = 2*p1mup5*(p3mup2+p4mup2)*p4mup3*F(120)+(mup3mup5*p3m
     -   up2+mup2mup5*p4mup3)*s234*F(123)
       Fa7(131) = -2*(p3mup2+p4mup2)*p4mup3*(p2mup5+p3mup5+p4mup5)+(mu
     -   p3mup5*p3mup2+mup2mup5*p4mup3)*s234
       Fa7(132) = p1mup5*(-2*(p2mup6*p3mup2+p3mup6*(p3mup2+p4mup2))+mu
     -   p2mup6*s34)*F(120)+F(123)*(-2*(p3mup2*p3mup5*p3mup6+p3mup5*p
     -   3mup6*p4mup2+p2mup6*p3mup5*(p3mup2+p4mup2)+p3mup2*p3mup5*p4m
     -   up6+p3mup5*p4mup2*p4mup6+p2mup5*p3mup2*(p2mup6+p3mup6+p4mup6
     -   ))+mup2mup6*p3mup5*s234+mup2mup5*p2mup6*s34+mup2mup5*p4mup6*
     -   s34+mup2mup5*p3mup6*s234*s34*P(39))
       Fa7(133) = p1mup6*(-2*(p2mup5*p3mup2+p3mup5*(p3mup2+p4mup2))+mu
     -   p2mup5*s34)*F(120)+F(123)*(-2*(p3mup2*p3mup5*p3mup6+p3mup5*p
     -   3mup6*p4mup2+p2mup5*(p2mup6*p3mup2+p3mup6*(p3mup2+p4mup2))+p
     -   3mup2*p3mup6*p4mup5+p3mup6*p4mup2*p4mup5+p2mup6*p3mup2*(p3mu
     -   p5+p4mup5))+mup5mup6*p3mup2*s234+mup2mup5*p3mup6*s234+mup2mu
     -   p6*p2mup5*s34+mup2mup6*p4mup5*s34+mup2mup6*p3mup5*s234*s34*P
     -   (39))
       Return
       End
