c************** Calling the Fa functions*************************
       subroutine NoAbe4vertBoxLineFa1(p1mup2,p1mup4,p1mup5,p2mup2,p2m
     &   up4,p2mup5,p3mup2,p3mup4,p3mup5,p4mup2,p4mup4,p4mup5,p5mup2,
     &   p5mup4,p5mup5,mup2mup4,mup2mup5,mup4mup5,Fa1)
       IMPLICIT NONE
      Complex*16   p1mup2, p1mup4, p1mup5, p2mup2, p2mup4, p2mup5, 
     -          p3mup2, p3mup4, p3mup5, p4mup2, p4mup4, p4mup5, 
     -          p5mup2, p5mup4, p5mup5
       Complex*16   mup2mup4, mup2mup5, mup4mup5
        common/NoAbe4vertBoxLineFhlFunctions/F
       COMMON/NoAbe4vertBoxLineInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,s12
     &   ,s23,s34,s45,s15
       Real*8   p1sq 
       Real*8   p2sq 
       Real*8   p3sq 
       Real*8   p4sq 
       Real*8   p5sq 
       Real*8   s12, s15 
       Real*8   s23 
       Real*8   s34 
       Real*8   s45
      Complex*16 Fa(27),F(21)
      Real*8 P(3) 
        COMMON/NoAbe4vertBoxLinePFunctions/P
       Complex*16 Fa1(1:27)
       COMMON/NoAbe4vertBoxLineFaFunctions/Fa
       Fa1(1) = p2mup2*F(6)+p1mup2*F(13)+p3mup2*F(17)
       Fa1(2) = p1mup4*p2mup2+p1mup2*p2mup4
       Fa1(3) = p1mup4*F(13)+p2mup4*F(14)
       Fa1(4) = p2mup2*p2mup4*F(5)+F(12)*Fa1(2)+p3mup2*Fa1(3)
       Fa1(5) = p1mup2*p1mup4*F(11)+mup2mup4*F(15)+4*Fa1(4)
       Fa1(6) = p2mup2*F(5)+p1mup2*F(12)+p3mup2*F(14)
       Fa1(7) = p1mup2*F(4)+p2mup2*F(5)+p3mup2*F(6)
       Fa1(8) = -4*p2mup5*p3mup4*F(16)+mup4mup5*F(20)
       Fa1(9) = p1mup5*p2mup2+p1mup2*p2mup5
       Fa1(10) = p1mup5*F(13)+p2mup5*F(14)
       Fa1(11) = p2mup2*p2mup5*F(5)+F(12)*Fa1(9)+p3mup2*Fa1(10)
       Fa1(12) = p1mup2*p1mup5*F(11)+mup2mup5*F(15)+4*Fa1(11)
       Fa1(13) = -(mup2mup5*p3mup4*F(16))+mup4mup5*Fa1(6)
       Fa1(14) = p2mup5+p3mup5
       Fa1(15) = p1mup4*p2mup5*s23*F(7)+p2mup4*F(1)*Fa1(14)
       Fa1(16) = mup2mup5*p1mup4*F(7)+(mup2mup4*F(1)*Fa1(14))/s23
       Fa1(17) = p2mup2+p3mup2
       Fa1(18) = -(mup2mup5*s23)+2*Fa1(14)*Fa1(17)
       Fa1(19) = mup2mup4*F(1)+2*p1mup4*F(7)*Fa1(17)
       Fa1(20) = p2mup4+p3mup4
       Fa1(21) = 2*mup4mup5*s23-4*p2mup5*Fa1(20)
       Fa1(22) = -4*p1mup5*p2mup4*s23*F(7)+F(1)*Fa1(21)
       Fa1(23) = mup2mup4*p1mup5*F(7)+(mup2mup5*F(1)*Fa1(20))/s23
       Fa1(24) = mup2mup5*F(1)+2*p1mup5*F(7)*Fa1(17)
       Fa1(25) = -(mup2mup4*s23)+2*Fa1(17)*Fa1(20)
       Fa1(26) = p1mup4+p2mup4
       Fa1(27) = p1mup5+p2mup5
       Return
       End
