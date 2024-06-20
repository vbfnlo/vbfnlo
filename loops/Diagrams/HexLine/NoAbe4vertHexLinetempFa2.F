c************** Calling the Fa functions*************************
       subroutine NoAbe4vertHexLineFa2(p1mup2,p1mup3,p1mup5,p1mup6,p2m
     -   up2,p2mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,
     -   p4mup3,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mu
     -   p3,p6mup5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup
     -   6,mup5mup6,Fa2)
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
       Complex*16 Fa2(20:38)
       COMMON/NoAbe4vertHexLineFaFunctions/Fa
       Fa2(20) = p3mup3*p4mup2*F(16)+p4mup2*p4mup3*F(22)+p1mup2*p3mup3
     -   *F(31)+p2mup2*p3mup3*F(35)-p3mup2*p4mup3*F(43)+mup2mup3*F(45
     -   )+p3mup2*p3mup3*F(49)+p1mup3*(p4mup2*F(14)+p1mup2*F(33)+p2mu
     -   p2*F(34)+p3mup2*F(50))+p2mup3*(p4mup2*F(13)+p2mup2*F(32)+p1m
     -   up2*F(34)+p3mup2*F(51))+p2mup2*p4mup3*F(52)+p1mup2*p4mup3*F(
     -   53)
       Fa2(21) = p3mup3*p4mup2*F(18)+p3mup2*p4mup3*F(18)+p1mup3*p4mup2
     -   *F(19)+p2mup3*p4mup2*F(20)+p4mup2*p4mup3*F(23)+p1mup2*(p4mup
     -   3*F(15)+p3mup3*F(36)+p2mup3*F(38)+p1mup3*F(41))+p2mup2*(p4mu
     -   p3*F(16)+p3mup3*F(37)+p1mup3*F(38)+p2mup3*F(42))+p1mup3*p3mu
     -   p2*F(46)+p2mup3*p3mup2*F(47)+p3mup2*p3mup3*F(48)+mup2mup3*F(
     -   138)
       Fa2(22) = mup5mup6*F(71)+2*p3mup6*p4mup5*F(132)
       Fa2(23) = -2*mup5mup6*F(73)+4*p2mup6*p4mup5*F(132)
       Fa2(24) = p4mup6*(p3mup3*F(2)+p1mup3*F(3)+p2mup3*F(4)+p4mup3*F(
     -   134))
       Fa2(25) = 2*mup3mup6*F(71)+4*(p1mup6*(p4mup3*F(3)+p1mup3*F(74)+
     -   p2mup3*F(75)+p3mup3*F(76))+p3mup6*(p4mup3*F(2)+p3mup3*F(77)+
     -   p1mup3*F(78)+p2mup3*F(79)))
       Fa2(26) = 2*(p2mup6*F(71)+p3mup6*F(73))+p1mup6*F(80)
       Fa2(27) = mup2mup6*F(71)+2*(p1mup6*p3mup2*F(81)+p3mup6*(p4mup2*
     -   F(4)+p3mup2*F(79)+p2mup2*F(82)+p1mup2*F(83)))
       Fa2(28) = -(mup3mup6*F(73))+2*p2mup6*(p4mup3*F(2)+p3mup3*F(77)+
     -   p1mup3*F(78)+p2mup3*F(79))
       Fa2(29) = p4mup6*(p1mup2*F(3)+p2mup2*F(4)+p3mup2*F(133)+p4mup2*
     -   F(134))
       Fa2(30) = -2*mup2mup6*F(73)+4*(p1mup6*(p4mup2*F(3)+p1mup2*F(74)
     -   +p2mup2*F(75)+p3mup2*F(78))+p2mup6*(p4mup2*F(4)+p3mup2*F(79)
     -   +p2mup2*F(82)+p1mup2*F(83)))
       Fa2(31) = 2*(p2mup5*F(71)+p3mup5*F(73))+p1mup5*F(80)
       Fa2(32) = 2*mup3mup5*F(71)+4*(p1mup5*(p4mup3*F(3)+p1mup3*F(74)+
     -   p2mup3*F(75)+p3mup3*F(76))+p3mup5*(p4mup3*F(2)+p3mup3*F(77)+
     -   p1mup3*F(78)+p2mup3*F(79)))
       Fa2(33) = p4mup5*(p3mup3*F(2)+p1mup3*F(3)+p2mup3*F(4)+p4mup3*F(
     -   134))
       Fa2(34) = mup2mup5*F(71)+2*(p1mup5*p3mup2*F(81)+p3mup5*(p4mup2*
     -   F(4)+p3mup2*F(79)+p2mup2*F(82)+p1mup2*F(83)))
       Fa2(35) = -(mup3mup5*F(73))+2*p2mup5*(p4mup3*F(2)+p3mup3*F(77)+
     -   p1mup3*F(78)+p2mup3*F(79))
       Fa2(36) = -2*mup2mup5*F(73)+4*(p1mup5*(p4mup2*F(3)+p1mup2*F(74)
     -   +p2mup2*F(75)+p3mup2*F(78))+p2mup5*(p4mup2*F(4)+p3mup2*F(79)
     -   +p2mup2*F(82)+p1mup2*F(83)))
       Fa2(37) = p4mup5*(p1mup2*F(3)+p2mup2*F(4)+p3mup2*F(133)+p4mup2*
     -   F(134))
       Fa2(38) = mup5mup6*(p3mup3*F(77)+p2mup3*F(79)-p1mup3*F(88)+p4mu
     -   p3*F(89))
       Return
       End
