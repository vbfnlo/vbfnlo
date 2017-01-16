c************** Calling the Fa functions*************************
       subroutine HexBox2Fa12(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup
     &   6,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa12)
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
       Complex*16 Fa12(210:228)
       COMMON/HexBox2FaFunctions/Fa
       Fa12(210) = (p5mup6+p6mup6)*(p4mup3*F(648)+p1mup3*F(653)-(p5mup
     &   3+p6mup3)*F(654)+p2mup3*F(724))
       Fa12(211) = -(p4mup6*((p5mup3+p6mup3)*F(650)+p2mup3*F(667)+p1mu
     &   p3*F(726)-p4mup3*F(728)))+mup3mup6*F(729)+p2mup6*(p1mup3*F(7
     &   25)-p4mup3*F(727)+(p5mup3+p6mup3)*F(730)+p2mup3*F(731))
       Fa12(212) = -(p2mup6*F(670))+p4mup6*F(676)
       Fa12(213) = p1mup3*F(733)+p2mup3*F(734)-(p5mup3+p6mup3)*F(735)+
     &   p4mup3*F(736)
       Fa12(214) = -2*(p5mup6+p6mup6)*(p2mup3*F(737)-p4mup3*F(738)+(p5
     &   mup3+p6mup3)*F(739)-p1mup3*F(741))+mup3mup6*F(740)*P(88)
       Fa12(215) = mup3mup6*F(742)-16*(p2mup3*((p5mup6+p6mup6)*F(646)-
     &   p2mup6*F(744)+2*p1mup6*F(745))-p1mup3*((p5mup6+p6mup6)*F(645
     &   )+p2mup6*F(743)-2*p1mup6*F(747))-(p5mup3+p6mup3)*((p5mup6+p6
     &   mup6)*F(650)+p2mup6*F(748)+2*p1mup6*F(749))-p4mup3*((p5mup6+
     &   p6mup6)*F(649)-p2mup6*F(746)+2*p1mup6*F(750)))
       Fa12(216) = p4mup3*F(649)+(p5mup3+p6mup3)*F(650)+p2mup3*F(667)+
     &   p1mup3*F(726)
       Fa12(217) = p1mup3*F(753)+p2mup3*F(754)+(p5mup3+p6mup3)*F(755)+
     &   p4mup3*F(756)
       Fa12(218) = mup3mup6*F(757)+16*(p2mup3*((p5mup6+p6mup6)*F(643)+
     &   p4mup6*F(743)-2*p2mup6*F(760))-p1mup3*((p5mup6+p6mup6)*F(653
     &   )-p4mup6*F(758)+2*p2mup6*F(761))-p4mup3*((p5mup6+p6mup6)*F(6
     &   48)+p4mup6*F(759)+2*p2mup6*F(763))+(p5mup3+p6mup3)*((p5mup6+
     &   p6mup6)*F(654)+p4mup6*F(762)-2*p2mup6*F(764)))
       Fa12(219) = mup3mup6*F(765)+2*(-(p1mup6*p4mup3*F(648))+p1mup3*(
     &   p4mup6*F(640)-p2mup6*F(643)-p1mup6*F(653))+p2mup6*p4mup3*F(7
     &   66)+p2mup3*(-(p4mup6*F(639))+p1mup6*F(643)+p2mup6*F(767))+(p
     &   5mup3+p6mup3)*(p1mup6*F(654)-p2mup6*F(768)))
       Fa12(220) = p2mup3*F(643)-p4mup3*F(648)-p1mup3*F(653)+(p5mup3+p
     &   6mup3)*F(654)
       Fa12(221) = p6mup3*F(770)+p1mup3*F(771)+p2mup3*F(772)+p4mup3*F(
     &   773)+p5mup3*F(774)
       Fa12(222) = p1mup3*F(775)+p2mup3*F(776)+p4mup3*F(777)+(p5mup3+p
     &   6mup3)*F(778)
       Fa12(223) = 2*p4mup3*F(722)+p1mup3*F(779)+(p5mup3+p6mup3)*F(780
     &   )+p2mup3*F(781)
       Fa12(224) = p4mup3*F(648)+p1mup3*F(653)-(p5mup3+p6mup3)*F(654)+
     &   p2mup3*F(724)
       Fa12(225) = p2mup3*F(785)+p1mup3*F(786)+(p5mup3+p6mup3)*F(787)+
     &   p4mup3*F(788)
       Fa12(226) = -(p5mup3*p5mup6*F(654))-p5mup6*p6mup3*F(654)-p5mup3
     &   *p6mup6*F(654)-p6mup3*p6mup6*F(654)+p2mup3*p5mup6*F(724)+p2m
     &   up3*p6mup6*F(724)-p2mup3*p4mup6*F(745)+p4mup6*p5mup3*F(749)+
     &   p4mup6*p6mup3*F(749)+p2mup3*p2mup6*F(760)+p1mup3*((p5mup6+p6
     &   mup6)*F(653)-p4mup6*F(747)+p2mup6*F(761))+p4mup3*((p5mup6+p6
     &   mup6)*F(648)+p4mup6*F(750)+p2mup6*F(763))+p2mup6*(p5mup3+p6m
     &   up3)*F(764)-mup3mup6*F(789)
       Fa12(227) = (p5mup6+p6mup6)*F(687)+p2mup6*F(790)
       Fa12(228) = p2mup3*F(745)+p1mup3*F(747)-(p5mup3+p6mup3)*F(749)-
     &   p4mup3*F(750)
       Return
       End
