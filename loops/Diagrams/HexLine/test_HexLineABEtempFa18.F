c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineAbe18(p1mup2,p1mup3,p1mup4,p1m
     -   up5,p2mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,
     -   p4mup2,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mu
     -   p2,p6mup3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,
     -   mup3mup5,mup4mup5,Fa18)
       IMPLICIT NONE
      Complex*16   p1mup2, p1mup3, p1mup4, p1mup5, p2mup2, p2mup3, 
     -          p2mup4, p2mup5, p3mup2, p3mup3, p3mup4, p3mup5, 
     -          p4mup2, p4mup3, p4mup4, p4mup5, p5mup2, p5mup3, 
     -          p5mup4, p5mup5, p6mup2, p6mup3, p6mup4, p6mup5
       Complex*16   mup2mup3, mup2mup4, mup2mup5, mup3mup4, mup3mup5, 
     -          mup4mup5
        common/test_FhlFunctionsAbe/F
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
      Complex*16 Fa(445),F(5930)
      Real*8 P(11674) 
        COMMON/test_PFunctionsAbe/P
       Complex*16 Fa18(188:198)
       COMMON/test_FaFunctionsAbe/Fa
       Fa18(188) = p4mup2*(p2mup3*(-(p2mup4*p2mup5*F(1112))+p1mup4*p6m
     -   up5*F(1113)+p2mup4*p6mup5*F(1113)+p2mup5*p6mup4*F(1117)-p6mu
     -   p4*p6mup5*F(1118))+p1mup5*(-(p2mup3*p2mup4*F(1112))+p2mup3*p
     -   6mup4*F(1117)+p2mup4*p6mup3*F(1119)+p6mup3*p6mup4*F(1125))+p
     -   6mup3*(p1mup4*p2mup5*F(1119)+p2mup4*p2mup5*F(1119)-p1mup4*p6
     -   mup5*F(1124)-p2mup4*p6mup5*F(1124)+p2mup5*p6mup4*F(1125)+p6m
     -   up4*p6mup5*F(1126)))
       Fa18(189) = mup3mup5*(p1mup4*p4mup2*F(1127)+p2mup4*p4mup2*F(112
     -   7)-p1mup4*p3mup2*F(1130)-p2mup4*p3mup2*F(1130)+p3mup4*p4mup2
     -   *F(1131)+16*(p3mup2*p3mup4*F(1128)-p3mup2*p4mup4*F(1129)-p3m
     -   up2*p6mup4*F(1132))+p4mup2*p4mup4*F(1133)+p4mup2*p6mup4*F(11
     -   34))
       Fa18(190) = mup4mup5*(p1mup3*p4mup2*s12*F(814)+p2mup3*p4mup2*s1
     -   2*F(814)-p3mup3*p4mup2*s12*F(815)+(p1mup2+p2mup2)*p4mup3*F(9
     -   96)+p1mup3*(p1mup2+p2mup2)*F(998)+(p1mup2+p2mup2)*p2mup3*F(9
     -   98)-(p1mup2+p2mup2)*p3mup3*F(1001))
       Fa18(191) = mup4mup5*(p4mup2*s12*(p4mup3*F(816)-p6mup3*F(968))+
     -   (p1mup2+p2mup2)*p6mup3*F(1002))
       Fa18(192) = mup3mup4*(-(p1mup2*p3mup5*F(1004))-(p2mup2*p2mup5+p
     -   1mup2*(p1mup5+p2mup5))*F(1006)+p4mup5*(p4mup2*s12*F(991)+p1m
     -   up2*F(1008)))
       Fa18(193) = mup3mup4*(-(p1mup5*p4mup2*s12*F(950))-p2mup5*p4mup2
     -   *s12*F(950)-p3mup5*p4mup2*s12*F(954)+p2mup2*p3mup5*F(1004)+p
     -   1mup5*p2mup2*F(1006)-p2mup2*p4mup5*F(1008)+p1mup2*p6mup5*F(1
     -   009)+p2mup2*p6mup5*F(1009)+p4mup2*p6mup5*s12*F(1135))
       Fa18(194) = p1mup2*(p1mup3*p3mup5*p4mup4*F(1014)+p2mup3*p3mup5*
     -   p4mup4*F(1014)+p1mup3*p3mup4*p4mup5*F(1021)+p2mup3*p3mup4*p4
     -   mup5*F(1021)-p3mup4*p4mup3*p4mup5*F(1032)-p3mup3*p4mup4*p4mu
     -   p5*F(1043)+p3mup5*p4mup3*p4mup4*F(1045)-p1mup4*p3mup5*p4mup3
     -   *F(1049)-p2mup4*p3mup5*p4mup3*F(1049)+p1mup5*p3mup4*p4mup3*F
     -   (1051)+p2mup5*p3mup4*p4mup3*F(1051)-p3mup3*p4mup4*p6mup5*F(1
     -   053)-p3mup4*p3mup5*p4mup3*F(1061)+p3mup4*p4mup5*p6mup3*F(106
     -   2)+p1mup5*p3mup3*p4mup4*F(1066)+p2mup5*p3mup3*p4mup4*F(1066)
     -   +p3mup3*p3mup4*p4mup5*F(1067)+p1mup4*p3mup3*p4mup5*F(1068)+p
     -   2mup4*p3mup3*p4mup5*F(1068)-p3mup4*p4mup3*p6mup5*F(1071)+p3m
     -   up5*p4mup3*p6mup4*F(1072)-p3mup5*p4mup4*p6mup3*F(1073)-p3mup
     -   3*p3mup5*p4mup4*F(1076)-p3mup3*p4mup5*p6mup4*F(1077)+(p1mup4
     -   *p1mup5*p2mup3+p1mup3*p1mup5*p2mup4+p1mup3*p1mup4*p2mup5)*F(
     -   1081))+p2mup2*(-(p3mup4*p4mup3*p4mup5*F(1032))-p3mup3*p4mup4
     -   *p4mup5*F(1043)+p3mup5*p4mup3*p4mup4*F(1045)-p1mup4*p3mup5*p
     -   4mup3*F(1049)+p1mup4*p3mup3*p4mup5*F(1068)+p1mup3*(p3mup5*p4
     -   mup4*F(1014)+p3mup4*p4mup5*F(1021)+p1mup4*p1mup5*F(1081)))
       Fa18(195) = p2mup2*(p2mup3*p3mup5*p4mup4*F(1014)+p2mup3*p3mup4*
     -   p4mup5*F(1021)-p2mup4*p3mup5*p4mup3*F(1049)+p1mup5*p3mup4*p4
     -   mup3*F(1051)+p2mup5*p3mup4*p4mup3*F(1051)-p3mup3*p4mup4*p6mu
     -   p5*F(1053)-p3mup4*p3mup5*p4mup3*F(1061)+p3mup4*p4mup5*p6mup3
     -   *F(1062)+p1mup5*p3mup3*p4mup4*F(1066)+p2mup5*p3mup3*p4mup4*F
     -   (1066)+p3mup3*p3mup4*p4mup5*F(1067)+p2mup4*p3mup3*p4mup5*F(1
     -   068)-p3mup4*p4mup3*p6mup5*F(1071)+p3mup5*p4mup3*p6mup4*F(107
     -   2)-p3mup5*p4mup4*p6mup3*F(1073)-p3mup3*p3mup5*p4mup4*F(1076)
     -   -p3mup3*p4mup5*p6mup4*F(1077))
       Fa18(196) = p2mup2*(p1mup3*(p1mup5*p3mup4*F(1082)+p2mup5*p3mup4
     -   *F(1082)+p1mup4*p3mup5*F(1082)+p2mup4*p3mup5*F(1082)+p3mup4*
     -   p3mup5*F(1084)-p3mup4*p6mup5*F(1085)+p3mup5*p6mup4*F(1090))+
     -   p1mup4*(p2mup3*p3mup5*F(1082)-p3mup3*p6mup5*F(1086)+p2mup5*p
     -   3mup3*F(1087)+p3mup3*p3mup5*F(1088)-p3mup5*p6mup3*F(1094))+p
     -   1mup5*(p2mup3*p3mup4*F(1082)+(p1mup4+p2mup4)*p3mup3*F(1087)+
     -   p3mup3*p3mup4*F(1088)-p3mup3*p6mup4*F(1094)+p3mup4*p6mup3*F(
     -   1102)))+p1mup2*(p2mup3*p2mup5*p3mup4*F(1082)+p1mup4*p2mup3*p
     -   3mup5*F(1082)+p2mup3*p2mup4*p3mup5*F(1082)+p2mup3*p3mup4*p3m
     -   up5*F(1084)-p2mup3*p3mup4*p6mup5*F(1085)-p1mup4*p3mup3*p6mup
     -   5*F(1086)-p2mup4*p3mup3*p6mup5*F(1086)+p1mup4*p2mup5*p3mup3*
     -   F(1087)+p2mup4*p2mup5*p3mup3*F(1087)+p2mup5*p3mup3*p3mup4*F(
     -   1088)+p1mup4*p3mup3*p3mup5*F(1088)+p2mup4*p3mup3*p3mup5*F(10
     -   88)+p2mup3*p3mup5*p6mup4*F(1090)+p1mup3*(p1mup5*p3mup4*F(108
     -   2)+p2mup5*p3mup4*F(1082)+p1mup4*p3mup5*F(1082)+p2mup4*p3mup5
     -   *F(1082)+p3mup4*p3mup5*F(1084)-p3mup4*p6mup5*F(1085)+p3mup5*
     -   p6mup4*F(1090))-p1mup4*p3mup5*p6mup3*F(1094)-p2mup4*p3mup5*p
     -   6mup3*F(1094)-p2mup5*p3mup3*p6mup4*F(1094)-p3mup3*p3mup4*p6m
     -   up5*F(1095)+p3mup3*p3mup4*p3mup5*F(1100)-p3mup4*p6mup3*p6mup
     -   5*F(1101)-p3mup3*p6mup4*p6mup5*F(1101)+p2mup5*p3mup4*p6mup3*
     -   F(1102)+p1mup5*(p2mup3*p3mup4*F(1082)+(p1mup4+p2mup4)*p3mup3
     -   *F(1087)+p3mup3*p3mup4*F(1088)-p3mup3*p6mup4*F(1094)+p3mup4*
     -   p6mup3*F(1102))-p3mup4*p3mup5*p6mup3*F(1103)-p3mup3*p3mup5*p
     -   6mup4*F(1103)-p3mup5*p6mup3*p6mup4*F(1136))
       Fa18(197) = p1mup2*(p1mup3*p4mup4*p4mup5*F(1016)+p2mup3*p4mup4*
     -   p4mup5*F(1016)-p1mup4*p4mup3*p4mup5*F(1017)-p2mup4*p4mup3*p4
     -   mup5*F(1017)+p1mup5*p4mup3*p4mup4*F(1029)+p2mup5*p4mup3*p4mu
     -   p4*F(1029)+p4mup3*p4mup4*p6mup5*F(1033)-p4mup3*p4mup5*p6mup4
     -   *F(1034)+p4mup4*p4mup5*p6mup3*F(1035)+p4mup3*p4mup4*p4mup5*F
     -   (1044))+p2mup2*(p4mup3*p4mup4*p4mup5*F(1044)-p2mup4*p3mup3*p
     -   6mup5*F(1086)+p2mup4*p2mup5*p3mup3*F(1087)+p2mup5*p3mup3*p3m
     -   up4*F(1088)+p2mup4*p3mup3*p3mup5*F(1088)+p2mup3*(p2mup5*p3mu
     -   p4*F(1082)+p2mup4*p3mup5*F(1082)+p3mup4*p3mup5*F(1084)-p3mup
     -   4*p6mup5*F(1085)+p3mup5*p6mup4*F(1090))-p2mup4*p3mup5*p6mup3
     -   *F(1094)-p2mup5*p3mup3*p6mup4*F(1094)-p3mup3*p3mup4*p6mup5*F
     -   (1095)+p3mup3*p3mup4*p3mup5*F(1100)-p3mup4*p6mup3*p6mup5*F(1
     -   101)-p3mup3*p6mup4*p6mup5*F(1101)+p2mup5*p3mup4*p6mup3*F(110
     -   2)-p3mup4*p3mup5*p6mup3*F(1103)-p3mup3*p3mup5*p6mup4*F(1103)
     -   -p3mup5*p6mup3*p6mup4*F(1136))
       Fa18(198) = p2mup2*(-(p1mup4*p4mup3*p4mup5*F(1017))-p2mup4*p4mu
     -   p3*p4mup5*F(1017)+p2mup3*p4mup5*(p4mup4*F(1016)+p1mup4*F(102
     -   0))-p1mup4*p1mup5*p4mup3*F(1023)-p1mup4*p2mup5*p4mup3*F(1023
     -   )+p1mup5*p4mup3*p4mup4*F(1029)+p2mup5*p4mup3*p4mup4*F(1029)+
     -   p4mup3*p4mup4*p6mup5*F(1033)-p4mup3*p4mup5*p6mup4*F(1034)+p4
     -   mup4*p4mup5*p6mup3*F(1035)+p1mup3*(p1mup5*p4mup4*F(1012)+p2m
     -   up5*p4mup4*F(1012)-p4mup4*p6mup5*F(1015)+p4mup4*p4mup5*F(101
     -   6)+p1mup4*p4mup5*F(1020)+p2mup4*p4mup5*F(1020)-p4mup5*p6mup4
     -   *F(1040))+p1mup4*p4mup5*p6mup3*F(1056)-p1mup4*p4mup3*p6mup5*
     -   F(1069))+p1mup2*(p2mup3*p2mup5*p4mup4*F(1012)-p2mup3*p4mup4*
     -   p6mup5*F(1015)+p1mup4*p2mup3*p4mup5*F(1020)+p2mup3*p2mup4*p4
     -   mup5*F(1020)-p1mup4*p2mup5*p4mup3*F(1023)-p2mup4*p2mup5*p4mu
     -   p3*F(1023)-p2mup3*p4mup5*p6mup4*F(1040)+p1mup3*((p1mup5+p2mu
     -   p5)*p4mup4*F(1012)-p4mup4*p6mup5*F(1015)+(p1mup4+p2mup4)*p4m
     -   up5*F(1020)-p4mup5*p6mup4*F(1040))+p1mup4*p4mup5*p6mup3*F(10
     -   56)+p2mup4*p4mup5*p6mup3*F(1056)+p2mup5*p4mup3*p6mup4*F(1058
     -   )-p1mup4*p4mup3*p6mup5*F(1069)-p2mup4*p4mup3*p6mup5*F(1069)-
     -   p2mup5*p4mup4*p6mup3*F(1070)+p1mup5*(p2mup3*p4mup4*F(1012)-(
     -   p1mup4+p2mup4)*p4mup3*F(1023)+p4mup3*p6mup4*F(1058)-p4mup4*p
     -   6mup3*F(1070))-p4mup4*p6mup3*p6mup5*F(1074)-p4mup3*p6mup4*p6
     -   mup5*F(1078)+p4mup5*p6mup3*p6mup4*F(1079))
       Return
       End
