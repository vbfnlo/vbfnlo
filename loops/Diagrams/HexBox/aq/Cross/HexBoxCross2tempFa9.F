c************** Calling the Fa functions*************************
       subroutine HexBoxCross2Fa9(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p
     &   3mup6,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa9
     &   )
       IMPLICIT NONE
      Complex*16   p1mup3, p1mup6, p2mup3, p2mup6, p3mup3, p3mup6, 
     -          p4mup3, p4mup6, p5mup3, p5mup6, p6mup3, p6mup6
       Complex*16   mup3mup6
        common/HexBoxCross2FhlFunctions/F
       COMMON/HexBoxCross2Invariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12
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
      Complex*16 Fa(392),F(1853)
      Real*8 P(1355) 
        COMMON/HexBoxCross2PFunctions/P
       Complex*16 Fa9(153:171)
       COMMON/HexBoxCross2FaFunctions/Fa
       Fa9(153) = p6mup6*F(406)+p2mup6*F(408)+p5mup6*F(411)+p1mup6*F(4
     &   12)
       Fa9(154) = p6mup6*F(317)+p5mup6*F(318)+p1mup6*F(409)+p2mup6*F(4
     &   10)
       Fa9(155) = (p1mup3+p5mup3+p6mup3)*(p2mup6*F(466)-p6mup6*F(513)+
     &   p1mup6*F(514)-p5mup6*F(515))
       Fa9(156) = (p1mup6+p6mup6)*F(517)+8*(p2mup6*F(518)-p5mup6*F(519
     &   ))
       Fa9(157) = (p1mup6+p6mup6)*F(520)+p4mup6*F(521)+p5mup6*F(522)-8
     &   *p2mup6*F(523)
       Fa9(158) = mup3mup6*F(524)+2*(p1mup6+p6mup6)*((p1mup3+p6mup3)*F
     &   (525)+p5mup3*F(526)+p4mup3*F(527)-p2mup3*F(528))
       Fa9(159) = mup3mup6*F(529)-16*(p4mup3*p4mup6*F(530)+p4mup6*p6mu
     &   p3*F(530)+p4mup6*p5mup3*F(531)-p2mup3*p4mup6*F(532)+p5mup6*p
     &   6mup3*F(533)+p4mup3*p5mup6*F(534)+p5mup3*p5mup6*F(535)+p1mup
     &   3*(p4mup6*F(530)+p5mup6*F(533)-2*p2mup6*F(536))-p2mup6*p5mup
     &   3*F(537)-p2mup6*p4mup3*F(538)+p2mup3*p5mup6*F(539)-2*(p2mup6
     &   *p6mup3*F(536)+p2mup3*p2mup6*F(540)))
       Fa9(160) = -(p4mup6*F(518))+(p1mup6+p6mup6)*F(543)-p5mup6*F(544
     &   )
       Fa9(161) = p2mup3*F(545)+(p1mup3+p6mup3)*F(546)+p4mup3*F(547)+p
     &   5mup3*F(548)
       Fa9(162) = 8*((p1mup3+p6mup3)*F(549)-p4mup3*F(550)-p5mup3*F(551
     &   ))+p2mup3*F(552)
       Fa9(163) = mup3mup6*F(553)-2*(p1mup6+p6mup6)*(p4mup3*F(554)+p5m
     &   up3*F(555))
       Fa9(164) = (p1mup3+p6mup3)*F(557)+p4mup3*F(558)+p5mup3*F(559)-p
     &   2mup3*F(560)
       Fa9(165) = 2*(p1mup6+p6mup6)*(p2mup3*F(562)+(p1mup3+p6mup3)*F(5
     &   63))-mup3mup6*F(561)*P(124)
       Fa9(166) = mup3mup6*F(564)+16*(p6mup3*p6mup6*F(533)+p4mup3*p6mu
     &   p6*F(534)+p5mup3*p6mup6*F(535)+p2mup3*p6mup6*F(539)+p1mup6*(
     &   p6mup3*F(533)+p4mup3*F(534)+p5mup3*F(535)+p2mup3*F(539))+p2m
     &   up3*p4mup6*F(565)+p4mup6*p6mup3*F(566)+p4mup6*p5mup3*F(567)+
     &   p2mup6*p6mup3*F(568)+p1mup3*((p1mup6+p6mup6)*F(533)+p4mup6*F
     &   (566)+p2mup6*F(568))+p2mup6*p5mup3*F(569)+p4mup3*p4mup6*F(57
     &   0)+p2mup6*p4mup3*F(571)+p2mup3*p2mup6*F(572))
       Fa9(167) = mup3mup6*F(573)+2*(p4mup3*p6mup6*F(530)+p6mup3*p6mup
     &   6*F(530)+p5mup3*p6mup6*F(531)-p2mup3*p6mup6*F(532)+p1mup6*((
     &   p4mup3+p6mup3)*F(530)+p5mup3*F(531)-p2mup3*F(532))+p2mup6*p6
     &   mup3*F(574)+p2mup6*p5mup3*F(575)+p2mup6*p4mup3*F(576)+p1mup3
     &   *((p1mup6+p6mup6)*F(530)+p2mup6*F(574)-2*p5mup6*F(577))-p2mu
     &   p3*p2mup6*F(578)-2*(p5mup3*p5mup6*F(567)+p4mup3*p5mup6*F(570
     &   )+p5mup6*p6mup3*F(577)-p2mup3*p5mup6*F(579)))
       Fa9(168) = mup3mup6*F(580)-16*(p4mup6*(p6mup3*F(525)+p5mup3*F(5
     &   26)+p4mup3*F(527)-p2mup3*F(528))+p1mup3*(p4mup6*F(525)+p5mup
     &   6*F(533)-p2mup6*F(536))+p5mup6*(p6mup3*F(533)+p4mup3*F(534)+
     &   p5mup3*F(535)+p2mup3*F(539))-p2mup6*(p6mup3*F(536)+p2mup3*F(
     &   540)+p5mup3*F(581)+p4mup3*F(582)))
       Fa9(169) = mup3mup6*F(586)-2*(p1mup6+p6mup6)*(p2mup3*F(587)+p4m
     &   up3*F(588)-(p1mup3+p6mup3)*F(589)+p5mup3*F(590))
       Fa9(170) = -(mup3mup6*F(591))+2*(p6mup3*p6mup6*F(536)+p2mup3*p6
     &   mup6*F(540)+p4mup6*p6mup3*F(574)+p4mup6*p5mup3*F(575)+p4mup3
     &   *p4mup6*F(576)-p2mup3*p4mup6*F(578)+p5mup3*p6mup6*F(581)+p4m
     &   up3*p6mup6*F(582)+p1mup6*(p6mup3*F(536)+p2mup3*F(540)+p5mup3
     &   *F(581)+p4mup3*F(582))+p1mup3*((p1mup6+p6mup6)*F(536)+p4mup6
     &   *F(574)+2*p5mup6*F(592))+2*(p5mup6*p6mup3*F(592)+p4mup3*p5mu
     &   p6*F(593)+p5mup3*p5mup6*F(594)-p2mup3*p5mup6*F(595)))
       Fa9(171) = -4*mup3mup6*F(599)+p2mup3*p6mup6*F(601)+p1mup6*(p2mu
     &   p3*F(601)-8*(p6mup3*F(600)-p5mup3*F(602)-p4mup3*F(603)))-8*(
     &   p6mup3*p6mup6*F(600)-p5mup3*p6mup6*F(602)-p4mup3*p6mup6*F(60
     &   3)+p2mup3*p2mup6*F(604)-p4mup6*p6mup3*F(605)-p4mup6*p5mup3*F
     &   (606)-p4mup3*p4mup6*F(607)-p5mup6*p6mup3*F(608)-p4mup3*p5mup
     &   6*F(609)-p5mup3*p5mup6*F(610)-p2mup6*p6mup3*F(611)-p1mup3*(-
     &   ((p1mup6+p6mup6)*F(600))+p4mup6*F(605)+p5mup6*F(608)+p2mup6*
     &   F(611))-p2mup6*p5mup3*F(612)-p2mup6*p4mup3*F(613)-p2mup3*p5m
     &   up6*F(614)-p2mup3*p4mup6*F(615))
       Return
       End
