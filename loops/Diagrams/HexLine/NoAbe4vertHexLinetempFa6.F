c************** Calling the Fa functions*************************
       subroutine NoAbe4vertHexLineFa6(p1mup2,p1mup3,p1mup5,p1mup6,p2m
     -   up2,p2mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,
     -   p4mup3,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mu
     -   p3,p6mup5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup
     -   6,mup5mup6,Fa6)
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
       Complex*16 Fa6(96:114)
       COMMON/NoAbe4vertHexLineFaFunctions/Fa
       Fa6(96) = -(mup3mup5*p2mup2*F(109))-4*mup3mup5*p3mup2*s12*F(110
     -   )+8*(p1mup5*p2mup2*p2mup3*F(104)+p2mup2*p2mup3*p2mup5*F(104)
     -   +p2mup2*p2mup3*p3mup5*F(105)+p1mup3*p2mup2*((p1mup5+p2mup5)*
     -   F(104)+p3mup5*F(105))+p1mup5*p2mup2*p4mup3*F(106)+p2mup2*p2m
     -   up5*p4mup3*F(106)+p1mup5*p2mup2*p3mup3*F(107)+p2mup2*p2mup5*
     -   p3mup3*F(107)+p2mup2*p3mup5*p4mup3*F(108)+p2mup2*p3mup3*p3mu
     -   p5*F(111))+p1mup2*(-(mup3mup5*F(109))+8*(p1mup5*p2mup3*F(104
     -   )+p2mup3*p2mup5*F(104)+p1mup3*(p1mup5+p2mup5)*F(104)+p1mup3*
     -   p3mup5*F(105)+p2mup3*p3mup5*F(105)+p1mup5*p4mup3*F(106)+p2mu
     -   p5*p4mup3*F(106)+p1mup5*p3mup3*F(107)+p2mup5*p3mup3*F(107)+p
     -   3mup5*p4mup3*F(108)+p3mup3*p3mup5*F(111)))
       Fa6(97) = (p1mup2+p2mup2)*p4mup5*((p1mup3+p2mup3)*F(106)+p3mup3
     -   *F(108)+p4mup3*F(129))
       Fa6(98) = (p1mup2+p2mup2)*p3mup5*p4mup6
       Fa6(99) = mup5mup6*(4*p3mup2*s12*F(110)+(p1mup2+p2mup2)*F(112))
       Fa6(100) = -4*mup5mup6*p3mup2*s12*F(110)+mup5mup6*(p1mup2+p2mup
     -   2)*F(113)+8*(p1mup2+p2mup2)*p3mup6*p4mup5*F(127)
       Fa6(101) = mup5mup6*p3mup2*((p1mup3+p2mup3)*F(107)+p3mup3*F(111
     -   )+p4mup3*F(131))
       Fa6(102) = p3mup2*(mup5mup6*((p1mup3+p2mup3)*F(105)+p4mup3*F(10
     -   8)+p3mup3*F(111))-mup3mup6*p4mup5*F(127))
       Fa6(103) = mup3mup5*(p1mup2+p2mup2)*p4mup6
       Fa6(104) = mup5mup6*(p1mup2+p2mup2)*((p1mup3+p2mup3)*F(107)+p3m
     -   up3*F(111)+p4mup3*F(131))
       Fa6(105) = (p1mup2+p2mup2)*(mup5mup6*((p1mup3+p2mup3)*F(105)+p4
     -   mup3*F(108)+p3mup3*F(111))-mup3mup6*p4mup5*F(127))
       Fa6(106) = -2*p3mup2*s12*F(110)+(p1mup2+p2mup2)*F(114)
       Fa6(107) = 2*p3mup2*s12*F(110)+(p1mup2+p2mup2)*F(115)
       Fa6(108) = p3mup2*((p1mup3+p2mup3)*F(105)+p4mup3*F(108)+p3mup3*
     -   F(111))
       Fa6(109) = (p1mup2+p2mup2)*((p1mup3+p2mup3)*F(105)+p4mup3*F(108
     -   )+p3mup3*F(111))
       Fa6(110) = mup5mup6*F(115)-4*p3mup6*p4mup5*F(127)
       Fa6(111) = p4mup6*((p1mup3+p2mup3)*F(106)+p3mup3*F(108)+p4mup3*
     -   F(129))
       Fa6(112) = 4*(p2mup3*p2mup6*F(104)+p1mup3*(p1mup6+p2mup6)*F(104
     -   )+p1mup3*p3mup6*F(105)+p2mup3*p3mup6*F(105)+p2mup6*p4mup3*F(
     -   106)+p2mup6*p3mup3*F(107)+p1mup6*(p2mup3*F(104)+p4mup3*F(106
     -   )+p3mup3*F(107))+p3mup6*p4mup3*F(108)+p3mup3*p3mup6*F(111))+
     -   mup3mup6*F(114)
       Fa6(113) = (p1mup6+p2mup6)*p3mup2
       Fa6(114) = (p1mup2+p2mup2)*p4mup6
       Return
       End
