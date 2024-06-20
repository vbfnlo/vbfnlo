c************** Calling the Fa functions*************************
       subroutine NoAbe4vertHexLineFa3(p1mup2,p1mup3,p1mup5,p1mup6,p2m
     -   up2,p2mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,
     -   p4mup3,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mu
     -   p3,p6mup5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup
     -   6,mup5mup6,Fa3)
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
       Complex*16 Fa3(39:57)
       COMMON/NoAbe4vertHexLineFaFunctions/Fa
       Fa3(39) = mup5mup6*(p4mup3*F(2)+p3mup3*F(77)+p1mup3*F(78)+p2mup
     -   3*F(79))-mup3mup6*p4mup5*F(132)
       Fa3(40) = p4mup3*F(2)+p3mup3*F(77)+p1mup3*F(78)+p2mup3*F(79)
       Fa3(41) = mup5mup6*(p3mup2*F(79)+p2mup2*F(82)-p1mup2*F(90)+p4mu
     -   p2*F(91))
       Fa3(42) = mup5mup6*(p4mup2*F(4)+p3mup2*F(79)+p2mup2*F(82)+p1mup
     -   2*F(83))-mup2mup6*p4mup5*F(132)
       Fa3(43) = mup2mup6*(p4mup3*F(2)+p3mup3*F(77)+p1mup3*F(78)+p2mup
     -   3*F(79))-mup3mup6*(p4mup2*F(4)+p3mup2*F(79)+p2mup2*F(82)+p1m
     -   up2*F(83))
       Fa3(44) = p4mup2*F(4)+p3mup2*F(79)+p2mup2*F(82)+p1mup2*F(83)
       Fa3(45) = mup2mup5*(p4mup3*F(2)+p3mup3*F(77)+p1mup3*F(78)+p2mup
     -   3*F(79))-mup3mup5*(p4mup2*F(4)+p3mup2*F(79)+p2mup2*F(82)+p1m
     -   up2*F(83))
       Fa3(46) = mup3mup5*p3mup6*F(92)-2*p2mup5*p4mup3*(p3mup6+p4mup6)
     -   *F(126)+mup3mup6*(-(p3mup5*F(92))+p2mup5*s34*F(126))
       Fa3(47) = (mup3mup6*p3mup5-mup3mup5*p3mup6+mup5mup6*p4mup3)*F(9
     -   2)-p2mup6*(2*p4mup3*(p3mup5+p4mup5)-mup3mup5*s34)*F(126)
       Fa3(48) = 8*p4mup3*(p3mup6+p4mup6)*(p1mup2*F(95)+p2mup2*F(96)+(
     -   p3mup2+p4mup2)*F(128))+4*mup3mup6*(p3mup2*F(94)-s34*(p1mup2*
     -   F(95)+p2mup2*F(96)+p4mup2*F(128)))
       Fa3(49) = -(mup2mup6*p4mup3*F(92))+mup3mup6*p1mup2*s34*F(95)+mu
     -   p3mup6*p2mup2*s34*F(96)+mup3mup6*p3mup2*F(97)+2*(p2mup6*(p3m
     -   up2+p4mup2)*p4mup3*F(96)+p1mup6*p4mup3*((p3mup2+p4mup2)*F(95
     -   )+p2mup2*F(98)+p1mup2*F(99))+p2mup2*p2mup6*p4mup3*F(100)+p1m
     -   up2*p2mup6*p4mup3*F(101))+mup3mup6*p4mup2*s34*F(128)
       Fa3(50) = -(mup2mup5*p4mup3*F(92))+mup3mup5*p1mup2*s34*F(95)+mu
     -   p3mup5*p2mup2*s34*F(96)+mup3mup5*p3mup2*F(97)+2*(p2mup5*(p3m
     -   up2+p4mup2)*p4mup3*F(96)+p1mup5*p4mup3*((p3mup2+p4mup2)*F(95
     -   )+p2mup2*F(98)+p1mup2*F(99))+p2mup2*p2mup5*p4mup3*F(100)+p1m
     -   up2*p2mup5*p4mup3*F(101))+mup3mup5*p4mup2*s34*F(128)
       Fa3(51) = 8*p4mup3*(p3mup5+p4mup5)*(p1mup2*F(95)+p2mup2*F(96)+(
     -   p3mup2+p4mup2)*F(128))+4*mup3mup5*(p3mup2*F(94)-s34*(p1mup2*
     -   F(95)+p2mup2*F(96)+p4mup2*F(128)))
       Fa3(52) = -4*mup2mup5*p3mup6*F(92)+8*(p1mup2*p3mup5*(p3mup6+p4m
     -   up6)*F(95)+p2mup2*p3mup5*p3mup6*F(96)+p2mup5*p3mup6*p4mup2*F
     -   (96)+p2mup2*p3mup5*p4mup6*F(96)+p1mup5*p3mup6*((p3mup2+p4mup
     -   2)*F(95)+p2mup2*F(98)+p1mup2*F(99))+p2mup2*p2mup5*p3mup6*F(1
     -   00)+p1mup2*p2mup5*p3mup6*F(101)-p2mup5*p3mup2*p4mup6*F(126)+
     -   p3mup2*p3mup5*p3mup6*F(128)+p3mup5*p3mup6*p4mup2*F(128)+p3mu
     -   p2*p3mup5*p4mup6*F(128)+p3mup5*p4mup2*p4mup6*F(128)+p2mup5*p
     -   3mup2*p3mup6*F(130))
       Fa3(53) = mup5mup6*(p1mup2*F(95)-(p3mup2*F(102))/s34+p4mup2*F(1
     -   28)+p2mup2*F(130))
       Fa3(54) = -4*(mup2mup6*p3mup5*F(92)-mup5mup6*(p3mup2*F(94)-s34*
     -   (p1mup2*F(95)+p2mup2*F(96)+p4mup2*F(128))))+8*(p1mup2*p3mup6
     -   *(p3mup5+p4mup5)*F(95)+p2mup2*p3mup5*p3mup6*F(96)+p2mup6*p3m
     -   up5*p4mup2*F(96)+p2mup2*p3mup6*p4mup5*F(96)+p1mup6*p3mup5*((
     -   p3mup2+p4mup2)*F(95)+p2mup2*F(98)+p1mup2*F(99))+p2mup2*p2mup
     -   6*p3mup5*F(100)+p1mup2*p2mup6*p3mup5*F(101)-p2mup6*p3mup2*p4
     -   mup5*F(126)+p3mup2*p3mup5*p3mup6*F(128)+p3mup5*p3mup6*p4mup2
     -   *F(128)+p3mup2*p3mup6*p4mup5*F(128)+p3mup6*p4mup2*p4mup5*F(1
     -   28)+p2mup6*p3mup2*p3mup5*F(130))
       Fa3(55) = -4*mup2mup5*mup3mup6*s34*F(126)+8*(mup3mup6*(p3mup5*(
     -   (p3mup2+p4mup2)*F(96)+p2mup2*F(100)+p1mup2*F(101))+p1mup5*p3
     -   mup2*F(103))+mup2mup5*p4mup3*(p3mup6+p4mup6)*F(126)-mup3mup5
     -   *(-(p3mup2*p4mup6*F(126))+p3mup6*(p4mup2*F(96)+p2mup2*F(100)
     -   +p1mup2*F(101)+p3mup2*F(130))))
       Fa3(56) = mup5mup6*p4mup3*(p1mup2*F(98)+p2mup2*F(100)+(p3mup2+p
     -   4mup2)*F(130))
       Fa3(57) = -4*mup2mup6*mup3mup5*s34*F(126)+8*(-(mup5mup6*p4mup3*
     -   ((p3mup2+p4mup2)*F(96)+p2mup2*F(100)+p1mup2*F(101)))+mup3mup
     -   5*(p3mup6*((p3mup2+p4mup2)*F(96)+p2mup2*F(100)+p1mup2*F(101)
     -   )+p1mup6*p3mup2*F(103))+mup2mup6*p4mup3*(p3mup5+p4mup5)*F(12
     -   6)-mup3mup6*(-(p3mup2*p4mup5*F(126))+p3mup5*(p4mup2*F(96)+p2
     -   mup2*F(100)+p1mup2*F(101)+p3mup2*F(130))))
       Return
       End
