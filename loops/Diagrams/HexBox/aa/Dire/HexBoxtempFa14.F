c************** Calling the Fa functions*************************
       subroutine HexBoxFa14(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4
     &   ,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa14)
       IMPLICIT NONE
      Complex*16   p1mup3, p1mup4, p2mup3, 
     -          p2mup4, p3mup3, p3mup4, 
     -          p4mup3, p4mup4, p5mup3, 
     -          p5mup4, p6mup3, p6mup4
       Complex*16   mup3mup4

        common/HexBoxFhlFunctions/F
       COMMON/HexBoxInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12,s23,s
     &   34,s45,s56,s16,s123,s234,s345
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
      Complex*16 Fa(259),F(1314)
      Real*8 P(1362) 
        COMMON/HexBoxPFunctions/P
       Complex*16 Fa14(157:168)
       COMMON/HexBoxFaFunctions/Fa
       Fa14(157) = mup3mup4*F(606)-2*(-(p1mup3*(p4mup4*F(607)+p5mup4*F
     &   (611)+p6mup4*F(612)-p1mup4*F(614)-p2mup4*F(615)))-p2mup3*(p4
     &   mup4*F(608)+p5mup4*F(609)+p6mup4*F(610)-p1mup4*F(613)-p2mup4
     &   *F(616))+p4mup4*p5mup3*F(619)+p4mup4*p6mup3*F(620)+p5mup3*p5
     &   mup4*F(621)+p5mup3*p6mup4*F(622)+p5mup4*p6mup3*F(623)+p6mup3
     &   *p6mup4*F(624)-p1mup4*p5mup3*F(626)-p1mup4*p6mup3*F(627)+p4m
     &   up3*(p5mup4*F(591)+p4mup4*F(617)+p6mup4*F(618)-p1mup4*F(625)
     &   -p2mup4*F(628))-p2mup4*p5mup3*F(629)-p2mup4*p6mup3*F(630))
       Fa14(158) = p1mup4*F(636)+p6mup4*F(637)+p4mup4*F(638)+p5mup4*F(
     &   639)+p2mup4*F(640)
       Fa14(159) = p1mup3*F(641)+p2mup3*F(642)+p6mup3*F(643)+p4mup3*F(
     &   644)+p5mup3*F(645)
       Fa14(160) = mup3mup4*F(648)+2*(-(p2mup4*(p4mup3*F(617)+p5mup3*F
     &   (619)+p6mup3*F(620)-p2mup3*F(628)))+p1mup3*(p2mup4*F(607)+p4
     &   mup4*F(650)+p5mup4*F(651)+p6mup4*F(652)+p1mup4*F(653))+p2mup
     &   3*p5mup4*F(654)+p2mup3*p6mup4*F(655)-p4mup3*p4mup4*F(656)-p4
     &   mup3*p5mup4*F(657)-p4mup3*p6mup4*F(658)-p4mup4*p5mup3*F(659)
     &   -p4mup4*p6mup3*F(660)-p5mup3*p5mup4*F(661)-p5mup3*p6mup4*F(6
     &   62)-p5mup4*p6mup3*F(663)-p6mup3*p6mup4*F(664)-p1mup4*(p2mup3
     &   *F(649)-p4mup3*F(665)-p5mup3*F(666)-p6mup3*F(667))+p2mup3*p4
     &   mup4*F(668))
       Fa14(161) = mup3mup4*F(671)-2*(p1mup3*(p4mup4*F(675)+p5mup4*F(6
     &   76)+p6mup4*F(677)+p2mup4*F(678)+p1mup4*F(679))+p2mup3*(p4mup
     &   4*F(672)+p5mup4*F(673)+p6mup4*F(674)-p1mup4*F(678)-p2mup4*F(
     &   680)))
       Fa14(162) = p1mup4*F(682)-p6mup4*F(683)-p4mup4*F(684)-p5mup4*F(
     &   685)+p2mup4*F(686)
       Fa14(163) = -(p1mup3*F(687))+p6mup3*F(688)+p4mup3*F(689)+p5mup3
     &   *F(690)-p2mup3*F(691)
       Fa14(164) = 8*(p1mup4*F(695)-p6mup4*F(696)-p4mup4*F(697)+p5mup4
     &   *F(698))+p2mup4*F(699)
       Fa14(165) = p1mup3*F(702)+p6mup3*F(703)+p4mup3*F(704)+p5mup3*F(
     &   705)-p2mup3*F(706)
       Fa14(166) = p4mup4*F(711)+4*(p1mup4*F(709)-p6mup4*F(710)+p5mup4
     &   *F(712)+p2mup4*F(713))
       Fa14(167) = 4*(p1mup4*F(714)+p2mup4*F(715)-p6mup4*F(716)-p5mup4
     &   *F(717))+p4mup4*F(718)
       Fa14(168) = p6mup4*F(719)+p5mup4*F(721)+4*(p4mup4*F(720)+p1mup4
     &   *F(722)+p2mup4*F(723))
       Return
       End
