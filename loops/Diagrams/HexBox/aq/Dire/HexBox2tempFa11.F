c************** Calling the Fa functions*************************
       subroutine HexBox2Fa11(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup
     &   6,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa11)
       IMPLICIT NONE
      Complex*16   p1mup3, p1mup6, p2mup3, p2mup6, p3mup3, p3mup6, 
     -          p4mup3, p4mup6, p5mup3, p5mup6, p6mup3, p6mup6
       Complex*16   mup3mup6
        common/HexBox2FhlFunctions/F
       COMMON/HexBox2Invariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12,s23,
     &   s34,s45,s56,s16,s123,s234,s345
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
      Complex*16 Fa(385),F(1646)
      Real*8 P(1185) 
        COMMON/HexBox2PFunctions/P
       Complex*16 Fa11(191:209)
       COMMON/HexBox2FaFunctions/Fa
       Fa11(191) = 4*(p5mup6*F(633)+p6mup6*F(634))+p2mup6*F(635)+p1mup
     &   6*F(636)+p4mup6*F(637)
       Fa11(192) = mup3mup6*F(638)+2*(p5mup6+p6mup6)*(p2mup3*F(639)-p1
     &   mup3*F(640))
       Fa11(193) = mup3mup6*F(642)+2*(-(p1mup3*p2mup6*F(644))+p1mup3*p
     &   4mup6*F(645)+p2mup6*p4mup3*F(647)+p4mup3*p4mup6*F(649)+p4mup
     &   6*(p5mup3+p6mup3)*F(650)+p2mup3*(-(p4mup6*F(646))+p2mup6*F(6
     &   51))-p2mup6*(p5mup3+p6mup3)*F(652)+p1mup6*(p2mup3*F(643)-p4m
     &   up3*F(648)-p1mup3*F(653)+(p5mup3+p6mup3)*F(654)))
       Fa11(194) = p2mup3*F(657)+p1mup3*F(658)+(p5mup3+p6mup3)*F(659)+
     &   p4mup3*F(660)
       Fa11(195) = p4mup6*F(661)+(p5mup6+p6mup6)*s56*F(662)-8*p1mup6*F
     &   (663)
       Fa11(196) = mup3mup6*F(664)+2*(p5mup6+p6mup6)*(p2mup3*F(665)+p1
     &   mup3*F(666)+p4mup3*F(667)+(p5mup3+p6mup3)*F(668))
       Fa11(197) = 16*p1mup6*F(670)-8*p4mup6*F(671)
       Fa11(198) = (p5mup6+p6mup6)*F(672)+p2mup6*F(673)+p1mup6*F(674)
       Fa11(199) = p2mup6*F(671)-2*p1mup6*F(676)
       Fa11(200) = p1mup6*F(677)+p2mup6*F(678)+p4mup6*F(679)+8*(p5mup6
     &   +p6mup6)*F(680)
       Fa11(201) = p4mup6*F(640)-p2mup6*F(681)
       Fa11(202) = p1mup3*F(644)-p4mup3*F(647)-p2mup3*F(651)+(p5mup3+p
     &   6mup3)*F(652)
       Fa11(203) = 2*(p4mup3*F(649)+(p5mup3+p6mup3)*F(650))+p1mup3*F(6
     &   85)-p2mup3*F(686)
       Fa11(204) = p2mup3*F(688)+p4mup3*F(689)+p1mup3*F(690)+(p5mup3+p
     &   6mup3)*F(691)
       Fa11(205) = p1mup3*F(692)+p2mup3*F(693)+(p5mup3+p6mup3)*F(694)+
     &   p4mup3*F(695)
       Fa11(206) = mup3mup6*F(698)+p2mup3*(p2mup6*F(701)+(p5mup6+p6mup
     &   6)*F(702)+8*p4mup6*F(703)+p1mup6*F(704))+p2mup6*p4mup3*F(705
     &   )+p4mup3*p5mup6*F(706)+p4mup3*p6mup6*F(706)+p5mup3*p5mup6*F(
     &   707)+p5mup6*p6mup3*F(707)+p5mup3*p6mup6*F(707)+p6mup3*p6mup6
     &   *F(707)+p1mup3*((p5mup6+p6mup6)*F(697)+p2mup6*F(699)-8*p4mup
     &   6*F(700)+p1mup6*F(708))+8*(-(p4mup3*p4mup6*F(709))+p4mup6*p5
     &   mup3*F(710)+p4mup6*p6mup3*F(710))+p2mup6*p5mup3*F(711)+p2mup
     &   6*p6mup3*F(711)+p1mup6*p5mup3*F(712)+p1mup6*p6mup3*F(712)+p1
     &   mup6*p4mup3*F(713)
       Fa11(207) = p5mup6*F(714)+p6mup6*F(715)+p4mup6*F(716)+p2mup6*F(
     &   717)
       Fa11(208) = 4*(p5mup6+p6mup6)*F(634)+p1mup6*F(636)+p2mup6*F(718
     &   )+p4mup6*F(719)
       Fa11(209) = p1mup3*F(720)+(p5mup3+p6mup3)*F(721)+4*p4mup3*F(722
     &   )+p2mup3*F(723)
       Return
       End
