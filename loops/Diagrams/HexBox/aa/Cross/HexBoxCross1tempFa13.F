c************** Calling the Fa functions*************************
       subroutine HexBoxCross1Fa13(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,
     &   p3mup4,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa
     &   13)
       IMPLICIT NONE
      Complex*16   p1mup3, p1mup4, p2mup3, 
     -          p2mup4, p3mup3, p3mup4, 
     -          p4mup3, p4mup4, p5mup3, 
     -          p5mup4, p6mup3, p6mup4
       Complex*16   mup3mup4

        common/HexBoxCross1FhlFunctions/F
       COMMON/HexBoxCross1Invariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12
     &   ,s23,s34,s45,s56,s16,s123,s234,s345
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
      Complex*16 Fa(272),F(1461)
      Real*8 P(1825) 
        COMMON/HexBoxCross1PFunctions/P
       Complex*16 Fa13(157:169)
       COMMON/HexBoxCross1FaFunctions/Fa
       Fa13(157) = p4mup3*F(547)+p1mup3*F(548)+p2mup3*F(549)+p6mup3*F(
     &   550)-p5mup3*F(551)
       Fa13(158) = -(p1mup3*F(552))+p2mup3*F(553)+p4mup3*F(554)+p5mup3
     &   *F(555)+p6mup3*F(1460)
       Fa13(159) = -(p1mup3*F(556))+p4mup3*F(557)+p5mup3*F(558)+p6mup3
     &   *F(559)+p2mup3*F(560)
       Fa13(160) = p2mup3*F(561)+8*(p1mup3*F(562)+p5mup3*F(563)-p6mup3
     &   *F(564)-p4mup3*F(565))
       Fa13(161) = p1mup4*F(566)+p5mup4*F(567)+p6mup4*F(568)-p2mup4*F(
     &   569)+p4mup4*F(570)
       Fa13(162) = p4mup3*F(571)+p1mup3*F(572)-p5mup3*F(573)-p6mup3*F(
     &   574)+p2mup3*F(575)
       Fa13(163) = p1mup4*F(578)-4*(p2mup4*F(579)-p4mup4*F(580)-p5mup4
     &   *F(581))+p6mup4*F(582)
       Fa13(164) = p4mup4*F(587)-p1mup4*F(588)+p5mup4*F(589)+p6mup4*F(
     &   590)-p2mup4*F(591)
       Fa13(165) = p5mup3*F(594)+p6mup3*F(595)+8*(p4mup3*F(592)+p1mup3
     &   *F(593)+p2mup3*F(596))
       Fa13(166) = p5mup3*F(560)+p1mup3*F(597)+p2mup3*F(598)+p4mup3*F(
     &   599)+p6mup3*F(600)
       Fa13(167) = mup3mup4*F(603)-2*(p5mup4*p6mup3*F(610)+p4mup4*p6mu
     &   p3*F(611)+p6mup3*p6mup4*F(612)-p2mup4*p6mup3*F(613)+p2mup3*(
     &   p4mup4*F(604)+p6mup4*F(605)+p5mup4*F(606)+p1mup4*F(614)+p2mu
     &   p4*F(615))-p4mup3*(p4mup4*F(607)+p5mup4*F(608)+p6mup4*F(609)
     &   +p2mup4*F(616)-p1mup4*F(617))+p5mup3*p5mup4*F(618)+p5mup3*p6
     &   mup4*F(619)+p4mup4*p5mup3*F(620)+p1mup3*p4mup4*F(621)-p1mup4
     &   *p6mup3*F(622)-p2mup4*p5mup3*F(623)+p1mup4*p5mup3*F(624)+p1m
     &   up3*p5mup4*F(624)+p1mup3*p6mup4*F(625)+p1mup3*p2mup4*F(626)-
     &   p1mup3*p1mup4*F(627))
       Fa13(168) = mup3mup4*F(631)+2*(p1mup3*(p5mup4*F(632)+p6mup4*F(6
     &   33)+p4mup4*F(634)-p2mup4*F(635)-p1mup4*F(638))+p5mup3*p5mup4
     &   *F(639)+p5mup3*p6mup4*F(640)+p4mup4*p5mup3*F(641)+p5mup4*p6m
     &   up3*F(643)+p4mup4*p6mup3*F(644)+p6mup3*p6mup4*F(645)+p1mup4*
     &   p5mup3*F(647)+p1mup4*p6mup3*F(648)+p4mup3*(p4mup4*F(636)+p5m
     &   up4*F(637)+p6mup4*F(642)+p1mup4*F(646)+p2mup4*F(649))+p2mup4
     &   *p5mup3*F(650)+p2mup3*p5mup4*F(650)+p2mup4*p6mup3*F(651)+p2m
     &   up3*p4mup4*F(652)+p2mup3*p6mup4*F(653)-p1mup4*p2mup3*F(654)-
     &   p2mup3*p2mup4*F(655))
       Fa13(169) = mup3mup4*F(659)+8*(-(p5mup3*p5mup4*F(666))-p5mup3*p
     &   6mup4*F(667)-p4mup4*p5mup3*F(668)-p5mup4*p6mup3*F(670)+p1mup
     &   3*(p5mup4*F(660)+p6mup4*F(661)+p4mup4*F(662)+p1mup4*F(665)+p
     &   2mup4*F(671))-p4mup4*p6mup3*F(672)+p2mup3*p4mup4*F(673)+p2mu
     &   p3*p6mup4*F(674)+p2mup3*p5mup4*F(675)-p6mup3*p6mup4*F(676)+p
     &   1mup4*p5mup3*F(678)+p1mup4*p6mup3*F(679)-p1mup4*p2mup3*F(680
     &   )-p4mup3*(p4mup4*F(663)+p5mup4*F(664)+p6mup4*F(669)-p1mup4*F
     &   (677)-p2mup4*F(681))+p2mup4*p5mup3*F(682)+p2mup4*p6mup3*F(68
     &   3)-p2mup3*p2mup4*F(684))
       Return
       End
