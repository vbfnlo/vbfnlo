c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineNoAbe1g17(p1mup2,p1mup3,p1mup4
     -   ,p1mup6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3m
     -   up6,p4mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,
     -   p6mup2,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3m
     -   up4,mup3mup6,mup4mup6,Fa17)
       IMPLICIT NONE
      Complex*16   p1mup2, p1mup3, p1mup4, p1mup6, p2mup2, p2mup3, 
     -          p2mup4, p2mup6, p3mup2, p3mup3, p3mup4, p3mup6, 
     -          p4mup2, p4mup3, p4mup4, p4mup6, p5mup2, p5mup3, 
     -          p5mup4, p5mup6, p6mup2, p6mup3, p6mup4, p6mup6
       Complex*16   mup2mup3, mup2mup4, mup2mup6, mup3mup4, mup3mup6, 
     -          mup4mup6
        common/test_FhlFunctionsNoAbe1g/F
       COMMON/Invariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12,s23,s34,s45
     -   ,s56,s16,s123,s234,s345
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
      Complex*16 Fa(587),F(6395)
      Real*8 P(7129) 
        COMMON/test_PFunctionsNoAbe1g/P
       Complex*16 Fa17(225:238)
       COMMON/test_FaFunctionsNoAbe1g/Fa
       Fa17(225) = p5mup4*(mup2mup3*((p4mup6+p5mup6)*F(626)+p3mup6*F(6
     -   27)+p2mup6*F(628)-p1mup6*F(629))-2*(p3mup2*(p4mup3+p5mup3)*(
     -   p4mup6+p5mup6)*F(630)-(p4mup2+p5mup2)*(p3mup3*(p4mup6+p5mup6
     -   )*F(631)-p3mup6*(p4mup3+p5mup3)*F(632))-(p4mup3*(p4mup2+p5mu
     -   p2)+p4mup2*p5mup3)*(p4mup6+p5mup6)*F(633)))
       Fa17(226) = p5mup4*(-(mup2mup6*(p2mup3*F(634)+8*(p3mup3*F(635)+
     -   p1mup3*F(636))))+16*(p1mup2*(p1mup3*p1mup6*F(637)+p1mup3*p2m
     -   up6*F(638)+p1mup6*p2mup3*F(639)+p1mup3*p3mup6*F(640)+p1mup6*
     -   p3mup3*F(641))-p1mup3*p1mup6*p3mup2*F(642)-p1mup3*p2mup6*p3m
     -   up2*F(647)-p1mup3*p3mup2*p3mup6*F(650)+p1mup3*p1mup6*p2mup2*
     -   F(653)+p1mup3*p2mup2*p2mup6*F(654)+p1mup6*p2mup2*p2mup3*F(65
     -   5)+p1mup3*p2mup2*p3mup6*F(656)+p1mup6*p2mup3*p3mup2*F(657)-p
     -   1mup3*p2mup2*p5mup6*F(658)-p4mup6*(p5mup2*p5mup3*F(633)-p1mu
     -   p2*(p1mup3*F(643)+p2mup3*F(648)+p3mup3*F(651))+p1mup3*(p3mup
     -   2*F(652)+p2mup2*F(658)))+p1mup2*(p2mup3*(p2mup6*F(644)+p3mup
     -   6*F(645))+p2mup6*p3mup3*F(646)+p5mup6*(p1mup3*F(643)+p2mup3*
     -   F(648))+p1mup6*(p4mup3+p5mup3)*F(659)+p2mup6*(p4mup3+p5mup3)
     -   *F(660)+p3mup6*(p3mup3*F(649)+p4mup3*F(661)))-p1mup3*p3mup6*
     -   p4mup2*F(662)+p1mup6*p2mup2*p3mup3*F(663)-p1mup6*p3mup2*p3mu
     -   p3*F(664)-p1mup3*p1mup6*p4mup2*F(665)-p1mup6*p2mup3*p4mup2*F
     -   (666)+p1mup3*p2mup6*p4mup2*F(667)-p5mup2*(p5mup3*p5mup6*F(63
     -   3)+p1mup3*p1mup6*F(665)+p1mup6*p2mup3*F(666)-p1mup3*p2mup6*F
     -   (667))+p1mup6*p2mup2*(p4mup3+p5mup3)*F(668)))
       Fa17(227) = p5mup4*(-(p1mup3*(p3mup2*p5mup6*F(652)+p3mup6*p5mup
     -   2*F(662)))+p2mup2*p2mup3*p2mup6*F(669)+p2mup2*p2mup3*p3mup6*
     -   F(670)-p2mup2*p2mup3*p4mup6*F(671)-p2mup2*p2mup3*p5mup6*F(67
     -   1)-p2mup2*p3mup3*p4mup6*F(672)-p2mup2*p3mup3*p5mup6*F(672)+p
     -   2mup3*p2mup6*p3mup2*F(673)+p2mup3*p3mup2*p3mup6*F(674)+p2mup
     -   3*p3mup2*p4mup6*F(675)+p2mup3*p3mup2*p5mup6*F(675)+p1mup6*p3
     -   mup2*p4mup3*F(676)+p1mup6*p3mup2*p5mup3*F(676)-p1mup6*p3mup3
     -   *p4mup2*F(677)-p1mup6*p3mup3*p5mup2*F(677)+p1mup2*(p3mup3*p5
     -   mup6*F(651)+p3mup6*p5mup3*F(661)+(p4mup3+p5mup3)*(p4mup6+p5m
     -   up6)*F(678))+p2mup6*p3mup2*p4mup3*F(679)+p2mup6*p3mup2*p5mup
     -   3*F(679)-p2mup6*p3mup3*p4mup2*F(680)-p2mup6*p3mup3*p5mup2*F(
     -   680)-p2mup3*p3mup6*p4mup2*F(681)-p2mup3*p3mup6*p5mup2*F(681)
     -   +p2mup6*p3mup2*p3mup3*F(682)+p3mup2*p3mup3*p3mup6*F(683)+p1m
     -   up3*(p4mup2+p5mup2)*(p4mup6+p5mup6)*F(684)+p2mup3*p4mup2*p4m
     -   up6*F(685)+p2mup2*p4mup3*p4mup6*F(685)+p2mup3*p4mup6*p5mup2*
     -   F(685)+p2mup2*p4mup6*p5mup3*F(685)+p2mup3*p4mup2*p5mup6*F(68
     -   5)+p2mup2*p4mup3*p5mup6*F(685)+p2mup3*p5mup2*p5mup6*F(685)+p
     -   2mup2*p5mup3*p5mup6*F(685)+p2mup2*p2mup6*p3mup3*F(686)+p2mup
     -   2*p3mup3*p3mup6*F(687)+p2mup3*p2mup6*p4mup2*F(688)+p2mup2*p2
     -   mup6*p4mup3*F(688)+p2mup3*p2mup6*p5mup2*F(688)+p2mup2*p2mup6
     -   *p5mup3*F(688)+p2mup2*p3mup6*p4mup3*F(689)+p2mup2*p3mup6*p5m
     -   up3*F(689)-p1mup6*p4mup2*p4mup3*F(690)-p1mup6*p4mup3*p5mup2*
     -   F(690)-p1mup6*p4mup2*p5mup3*F(690)-p1mup6*p5mup2*p5mup3*F(69
     -   0)+p2mup6*(p4mup2+p5mup2)*(p4mup3+p5mup3)*F(691))
       Fa17(228) = p5mup4*(p3mup2*p3mup6*(p4mup3+p5mup3)*F(692)-p3mup3
     -   *p3mup6*(p4mup2+p5mup2)*F(693)+p3mup2*p3mup3*(p4mup6+p5mup6)
     -   *F(694))
       Fa17(229) = p3mup3*p4mup2*F(696)+p2mup3*p4mup2*F(697)+p4mup2*p5
     -   mup3*F(698)+p4mup3*p5mup2*F(699)+p1mup2*p4mup3*F(701)+p1mup3
     -   *p4mup2*F(702)+p4mup2*p4mup3*F(703)+p3mup2*(p4mup3*F(695)+4*
     -   s45*(-(p5mup3*F(600))+p2mup3*F(704)+p3mup3*F(705)+p1mup3*F(7
     -   06)))+p2mup2*(p4mup3*F(700)+4*s45*(p2mup3*F(599)+p3mup3*F(60
     -   6)+p1mup3*F(609)+p5mup3*F(707)))+4*(p3mup3*p5mup2*s45*F(601)
     -   +p5mup2*p5mup3*s45*F(602)+p1mup2*p1mup3*s45*F(607)+p1mup2*p2
     -   mup3*s45*F(608)+p1mup2*p3mup3*s45*F(610)+p2mup3*p5mup2*s45*F
     -   (611)+mup2mup3*s45*F(614)+p1mup3*p5mup2*s45*F(618)+p1mup2*p5
     -   mup3*s45*F(708))
       Fa17(230) = mup4mup6*(p4mup2*p5mup3*F(711)-8*(p3mup2*p4mup3*F(7
     -   09)-p3mup3*p4mup2*F(710)+p4mup2*p4mup3*F(712)-p2mup2*p2mup3*
     -   s45*F(713)-p3mup2*p5mup3*s45*F(714)-p3mup3*p5mup2*s45*F(715)
     -   -p5mup2*p5mup3*s45*F(716)))
       Fa17(231) = (8*mup4mup6*(p4mup3*(p5mup2*F(717)+p2mup2*F(718))-p
     -   2mup3*p4mup2*F(719)-s45*(p1mup2*(p1mup3*F(721)+p2mup3*F(722)
     -   +p3mup3*F(724))+p2mup2*(p3mup3*F(720)+p1mup3*F(723)+p5mup3*F
     -   (725))-p2mup3*p3mup2*F(726)+p2mup3*p5mup2*F(727)-p3mup2*p3mu
     -   p3*F(728)))+mup2mup3*mup4mup6*s45*F(729))/s45
       Fa17(232) = mup4mup6*(-(p1mup2*p4mup3*F(730))+p1mup3*p4mup2*F(7
     -   31)+p1mup2*p5mup3*s45*F(732)-p1mup3*p3mup2*s45*F(733)+p1mup3
     -   *p5mup2*s45*F(734))
       Fa17(233) = mup3mup6*p5mup4*(p1mup2*F(735)+p2mup2*F(736)+p3mup2
     -   *F(737)+p4mup2*F(738))
       Fa17(234) = p5mup4*(mup3mup6*p5mup2*F(738)-16*mup2mup6*(p4mup3+
     -   p5mup3)*F(739))
       Fa17(235) = p5mup4*(mup2mup3*((p4mup6+p5mup6)*F(740)+p3mup6*F(7
     -   41)+p2mup6*F(742)-p1mup6*F(743))+2*(p3mup3*(p4mup2+p5mup2)*(
     -   p4mup6+p5mup6)*F(744)+(p4mup3+p5mup3)*(p3mup2*(p4mup6+p5mup6
     -   )*F(745)-p3mup6*(p4mup2+p5mup2)*F(746))+(p4mup3*(p4mup2+p5mu
     -   p2)+p4mup2*p5mup3)*(p4mup6+p5mup6)*F(747)))
       Fa17(236) = p5mup4*(p5mup2*p5mup3*p5mup6*F(747)-mup2mup6*p1mup3
     -   *F(748)-mup2mup6*p2mup3*F(749)-p1mup2*p1mup3*p1mup6*F(750)-p
     -   1mup2*p1mup3*p2mup6*F(751)-p1mup3*p1mup6*p2mup2*F(752)-p1mup
     -   2*p1mup3*p3mup6*F(753)-p1mup2*p1mup3*p5mup6*F(754)-p1mup3*p2
     -   mup2*p2mup6*F(755)-p1mup3*p2mup2*p3mup6*F(756)-p1mup3*p2mup2
     -   *p5mup6*F(757)+p1mup2*p1mup6*p2mup3*F(758)-p1mup2*p2mup3*p2m
     -   up6*F(759)+p1mup6*p2mup2*p2mup3*F(760)-p1mup2*p2mup3*p3mup6*
     -   F(761)-p1mup2*p2mup3*p5mup6*F(762)-p1mup3*p1mup6*p3mup2*F(76
     -   3)-p1mup2*p1mup6*p3mup3*F(763)-p1mup2*p2mup6*p3mup3*F(764)-p
     -   1mup3*p2mup6*p3mup2*F(765)-p1mup6*p2mup3*p3mup2*F(766)-p1mup
     -   3*p3mup2*p3mup6*F(767)-p1mup2*p3mup3*p3mup6*F(767)+p4mup6*(p
     -   5mup2*p5mup3*F(747)-p1mup2*p1mup3*F(754)-p1mup3*p2mup2*F(757
     -   )-p1mup2*p2mup3*F(762)+p1mup2*p3mup3*F(768)-p1mup3*p3mup2*F(
     -   769))-p1mup3*p1mup6*p4mup2*F(770)-p1mup3*p1mup6*p5mup2*F(770
     -   )-p1mup3*p2mup6*p4mup2*F(771)-p1mup3*p2mup6*p5mup2*F(771)+p1
     -   mup6*p2mup2*p4mup3*F(772)+p1mup6*p2mup2*p5mup3*F(772)-p1mup6
     -   *p2mup3*p4mup2*F(773)-p1mup6*p2mup3*p5mup2*F(773)-p1mup3*p3m
     -   up6*p4mup2*F(774)-p1mup2*p1mup6*p4mup3*F(775)-p1mup2*p1mup6*
     -   p5mup3*F(775)-p1mup2*p2mup6*p4mup3*F(776)-p1mup2*p2mup6*p5mu
     -   p3*F(776)-p1mup2*p3mup6*p4mup3*F(777)-mup2mup6*p3mup3*F(778)
     -   -p1mup6*p2mup2*p3mup3*F(779)-p1mup6*p3mup2*p3mup3*F(780))
       Fa17(237) = p5mup4*(p2mup2*p2mup3*p2mup6*F(781)+p2mup2*p2mup3*p
     -   3mup6*F(782)+p2mup2*p2mup3*p4mup6*F(783)+p2mup2*p2mup3*p5mup
     -   6*F(783)+p2mup3*p2mup6*p3mup2*F(784)+p2mup3*p3mup2*p3mup6*F(
     -   785)-p2mup2*p3mup3*p4mup6*F(786)-p2mup2*p3mup3*p5mup6*F(786)
     -   +p2mup3*p3mup2*p4mup6*F(787)+p2mup3*p3mup2*p5mup6*F(787)+p1m
     -   up3*(p3mup2*p5mup6*F(769)+p3mup6*p5mup2*F(774)+(p4mup2+p5mup
     -   2)*(p4mup6+p5mup6)*F(788))-p2mup2*p2mup6*p4mup3*F(789)-p2mup
     -   2*p2mup6*p5mup3*F(789)+p2mup3*p2mup6*p4mup2*F(790)+p2mup3*p2
     -   mup6*p5mup2*F(790)-p2mup2*p3mup6*p4mup3*F(791)-p2mup2*p3mup6
     -   *p5mup3*F(791)+p2mup3*p3mup6*p4mup2*F(792)+p2mup3*p3mup6*p5m
     -   up2*F(792)-p2mup2*p4mup3*p4mup6*F(793)-p2mup2*p4mup6*p5mup3*
     -   F(793)-p2mup2*p4mup3*p5mup6*F(793)-p2mup2*p5mup3*p5mup6*F(79
     -   3)+p2mup3*p4mup2*p4mup6*F(794)+p2mup3*p4mup6*p5mup2*F(794)+p
     -   2mup3*p4mup2*p5mup6*F(794)+p2mup3*p5mup2*p5mup6*F(794)+p1mup
     -   6*p3mup3*p4mup2*F(795)+p1mup6*p3mup3*p5mup2*F(795)+p1mup2*(-
     -   (p3mup3*p5mup6*F(768))+p3mup6*p5mup3*F(777)-(p4mup3+p5mup3)*
     -   (p4mup6+p5mup6)*F(796))+p2mup6*p3mup3*p4mup2*F(797)+p2mup6*p
     -   3mup3*p5mup2*F(797)+p2mup6*p3mup2*p3mup3*F(798)+p3mup2*p3mup
     -   3*p3mup6*F(799)+p2mup2*p2mup6*p3mup3*F(800)+p2mup2*p3mup3*p3
     -   mup6*F(801)+p1mup6*p3mup2*p4mup3*F(802)+p1mup6*p3mup2*p5mup3
     -   *F(802)+p2mup6*p3mup2*p4mup3*F(803)+p2mup6*p3mup2*p5mup3*F(8
     -   03)-p1mup6*p4mup2*p4mup3*F(804)-p1mup6*p4mup3*p5mup2*F(804)-
     -   p1mup6*p4mup2*p5mup3*F(804)-p1mup6*p5mup2*p5mup3*F(804)-p2mu
     -   p6*(p4mup2+p5mup2)*(p4mup3+p5mup3)*F(805))
       Fa17(238) = p5mup4*(p3mup2*p3mup3*(p4mup6+p5mup6)*F(806)+p3mup3
     -   *p3mup6*(p4mup2+p5mup2)*F(807)+p3mup2*p3mup6*(p4mup3+p5mup3)
     -   *F(808))
       Return
       End
