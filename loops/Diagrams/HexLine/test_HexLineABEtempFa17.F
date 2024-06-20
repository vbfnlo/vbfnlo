c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineAbe17(p1mup2,p1mup3,p1mup4,p1m
     -   up5,p2mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,
     -   p4mup2,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mu
     -   p2,p6mup3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,
     -   mup3mup5,mup4mup5,Fa17)
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
       Complex*16 Fa17(177:187)
       COMMON/test_FaFunctionsAbe/Fa
       Fa17(177) = p1mup5*(p1mup2+p2mup2)*F(987)+(p1mup2+p2mup2)*p2mup
     -   5*F(987)+(p1mup2+p2mup2)*(p4mup5*F(986)-p6mup5*F(988))-p1mup
     -   2*p3mup5*F(989)-p2mup2*p3mup5*F(989)+p4mup2*p6mup5*s12*F(990
     -   )+3*(p1mup5*p4mup2*s12*F(950)+p2mup5*p4mup2*s12*F(950)+p3mup
     -   5*p4mup2*s12*F(954)+p4mup2*p4mup5*s12*F(991))
       Fa17(178) = -(p6mup5*F(988))+p4mup5*F(992)+(p1mup5+p2mup5)*F(99
     -   3)-p3mup5*F(994)
       Fa17(179) = mup4mup5*((p1mup3+p2mup3)*p4mup2*F(995)+16*p3mup2*p
     -   4mup3*F(996)+p3mup3*p4mup2*F(997))
       Fa17(180) = mup4mup5*(p4mup2*(p4mup3*F(999)-p6mup3*F(1000))+p3m
     -   up2*((p1mup3+p2mup3)*F(998)-p3mup3*F(1001)+p6mup3*F(1002)))
       Fa17(181) = mup3mup4*(p4mup2*p4mup5*F(1003)-16*p3mup2*p3mup5*F(
     -   1004))
       Fa17(182) = mup3mup4*(p1mup5*p4mup2*F(1005)+p2mup5*p4mup2*F(100
     -   5)+p3mup5*p4mup2*F(1007)-16*(p1mup5*p3mup2*F(1006)+p2mup5*p3
     -   mup2*F(1006)-p3mup2*p4mup5*F(1008)+p3mup2*p6mup5*F(1009)+p4m
     -   up2*p6mup5*F(1010)))
       Fa17(183) = p2mup3*p3mup2*p4mup4*p4mup5*F(1016)-p2mup4*p3mup2*p
     -   4mup3*p4mup5*F(1017)+p2mup4*p3mup3*p4mup2*p4mup5*F(1018)-p3m
     -   up3*p4mup2*p4mup4*p6mup5*F(1019)-p3mup4*p4mup2*p4mup3*p4mup5
     -   *F(1025)-p2mup3*p3mup4*p4mup2*p4mup5*F(1026)-p2mup3*p3mup5*p
     -   4mup2*p4mup4*F(1027)-p2mup4*p3mup5*p4mup2*p4mup3*F(1028)+p1m
     -   up4*(p2mup3*(p3mup5*p4mup2*F(1011)+p3mup2*p4mup5*F(1020))-p3
     -   mup2*p4mup3*(p4mup5*F(1017)+(p1mup5+p2mup5)*F(1023))+p3mup3*
     -   p4mup2*(p4mup5*F(1018)+p1mup5*F(1024))-p3mup5*p4mup2*p4mup3*
     -   F(1028))+p1mup5*p3mup2*p4mup3*p4mup4*F(1029)+p2mup5*p3mup2*p
     -   4mup3*p4mup4*F(1029)-p1mup5*p3mup3*p4mup2*p4mup4*F(1030)-p2m
     -   up5*p3mup3*p4mup2*p4mup4*F(1030)-p1mup5*p3mup4*p4mup2*p4mup3
     -   *F(1031)-p2mup5*p3mup4*p4mup2*p4mup3*F(1031)-p3mup3*p3mup4*p
     -   4mup2*p4mup5*F(1032)-p3mup2*p3mup4*p4mup3*p4mup5*F(1032)+p3m
     -   up2*p4mup3*p4mup4*p6mup5*F(1033)-p3mup2*p4mup3*p4mup5*p6mup4
     -   *F(1034)+p3mup2*p4mup4*p4mup5*p6mup3*F(1035)+p3mup3*p4mup2*p
     -   4mup5*p6mup4*F(1036)-p3mup4*p3mup5*p4mup2*p4mup3*F(1037)+p3m
     -   up4*p4mup2*p4mup3*p6mup5*F(1038)-p3mup4*p4mup2*p4mup5*p6mup3
     -   *F(1039)+p1mup3*(p1mup4*p3mup5*p4mup2*F(1011)+p2mup4*p3mup5*
     -   p4mup2*F(1011)+(p1mup5+p2mup5)*p3mup2*p4mup4*F(1012)-p1mup5*
     -   p3mup4*p4mup2*F(1013)-p2mup5*p3mup4*p4mup2*F(1013)+p3mup2*p3
     -   mup5*p4mup4*F(1014)-p3mup2*p4mup4*p6mup5*F(1015)+p3mup2*p4mu
     -   p4*p4mup5*F(1016)+p1mup4*p3mup2*p4mup5*F(1020)+p2mup4*p3mup2
     -   *p4mup5*F(1020)+p3mup2*p3mup4*p4mup5*F(1021)-p3mup4*p4mup2*p
     -   6mup5*F(1022)-p3mup4*p4mup2*p4mup5*F(1026)-p3mup5*p4mup2*p4m
     -   up4*F(1027)-p3mup2*p4mup5*p6mup4*F(1040)-p3mup4*p3mup5*p4mup
     -   2*F(1041)-p3mup5*p4mup2*p6mup4*F(1042))-p3mup2*p3mup3*p4mup4
     -   *p4mup5*F(1043)+p3mup3*p4mup2*p4mup4*p4mup5*F(1044)+p3mup2*p
     -   4mup3*p4mup4*p4mup5*F(1044)+p3mup3*p3mup5*p4mup2*p4mup4*F(10
     -   45)+p3mup2*p3mup5*p4mup3*p4mup4*F(1045)-p3mup5*p4mup2*p4mup3
     -   *p4mup4*F(1046)-p3mup5*p4mup2*p4mup3*p6mup4*F(1047)-p3mup5*p
     -   4mup2*p4mup4*p6mup3*F(1048)
       Fa17(184) = -(p2mup4*p2mup5*p3mup2*p4mup3*F(1023))+p1mup4*p2mup
     -   5*p3mup3*p4mup2*F(1024)+p2mup4*p2mup5*p3mup3*p4mup2*F(1024)-
     -   p2mup3*(-((p1mup5+p2mup5)*p3mup2*p4mup4*F(1012))+p1mup5*p3mu
     -   p4*p4mup2*F(1013)+p2mup5*p3mup4*p4mup2*F(1013)-p3mup2*p3mup5
     -   *p4mup4*F(1014)+p3mup2*p4mup4*p6mup5*F(1015)-p2mup4*(p3mup5*
     -   p4mup2*F(1011)+p3mup2*p4mup5*F(1020))-p3mup2*p3mup4*p4mup5*F
     -   (1021)+p3mup4*p4mup2*p6mup5*F(1022)+p3mup2*p4mup5*p6mup4*F(1
     -   040)+p3mup4*p3mup5*p4mup2*F(1041)+p3mup5*p4mup2*p6mup4*F(104
     -   2))-p1mup4*p3mup2*p3mup5*p4mup3*F(1049)-p2mup4*p3mup2*p3mup5
     -   *p4mup3*F(1049)+p1mup4*p3mup3*p3mup5*p4mup2*F(1050)+p2mup4*p
     -   3mup3*p3mup5*p4mup2*F(1050)+p2mup5*p3mup2*p3mup4*p4mup3*F(10
     -   51)-p2mup5*p3mup3*p3mup4*p4mup2*F(1052)-p3mup2*p3mup3*p4mup4
     -   *p6mup5*F(1053)-p1mup4*p3mup3*p4mup2*p6mup5*F(1054)-p2mup4*p
     -   3mup3*p4mup2*p6mup5*F(1054)-p3mup3*p3mup4*p4mup2*p6mup5*F(10
     -   55)+p1mup4*p3mup2*p4mup5*p6mup3*F(1056)+p2mup4*p3mup2*p4mup5
     -   *p6mup3*F(1056)-p1mup4*p3mup5*p4mup2*p6mup3*F(1057)-p2mup4*p
     -   3mup5*p4mup2*p6mup3*F(1057)+p2mup5*p3mup2*p4mup3*p6mup4*F(10
     -   58)-p2mup5*p3mup3*p4mup2*p6mup4*F(1059)-p2mup5*p3mup4*p4mup2
     -   *p6mup3*F(1060)-p3mup3*p3mup4*p3mup5*p4mup2*F(1061)-p3mup2*p
     -   3mup4*p3mup5*p4mup3*F(1061)+p3mup2*p3mup4*p4mup5*p6mup3*F(10
     -   62)+p3mup3*p4mup2*p6mup4*p6mup5*F(1063)-p3mup4*p3mup5*p4mup2
     -   *p6mup3*F(1064)+p3mup4*p4mup2*p6mup3*p6mup5*F(1065)+p2mup5*p
     -   3mup2*p3mup3*p4mup4*F(1066)+p3mup2*p3mup3*p3mup4*p4mup5*F(10
     -   67)+p1mup4*p3mup2*p3mup3*p4mup5*F(1068)+p2mup4*p3mup2*p3mup3
     -   *p4mup5*F(1068)-p1mup4*p3mup2*p4mup3*p6mup5*F(1069)-p2mup4*p
     -   3mup2*p4mup3*p6mup5*F(1069)-p2mup5*p3mup2*p4mup4*p6mup3*F(10
     -   70)-p1mup5*(p2mup4*(p3mup2*p4mup3*F(1023)-p3mup3*p4mup2*F(10
     -   24))+p4mup2*(p3mup3*p3mup4*F(1052)+p3mup3*p6mup4*F(1059)+p3m
     -   up4*p6mup3*F(1060))-p3mup2*(p3mup4*p4mup3*F(1051)+p4mup3*p6m
     -   up4*F(1058)+p3mup3*p4mup4*F(1066)-p4mup4*p6mup3*F(1070)))-p3
     -   mup2*p3mup4*p4mup3*p6mup5*F(1071)+p3mup2*p3mup5*p4mup3*p6mup
     -   4*F(1072)-p3mup2*p3mup5*p4mup4*p6mup3*F(1073)-p3mup2*p4mup4*
     -   p6mup3*p6mup5*F(1074)-p3mup3*p3mup5*p4mup2*p6mup4*F(1075)-p3
     -   mup2*p3mup3*p3mup5*p4mup4*F(1076)-p3mup2*p3mup3*p4mup5*p6mup
     -   4*F(1077)-p3mup2*p4mup3*p6mup4*p6mup5*F(1078)+p3mup2*p4mup5*
     -   p6mup3*p6mup4*F(1079)-p3mup5*p4mup2*p6mup3*p6mup4*F(1080)
       Fa17(185) = p3mup2*(p1mup5*p2mup3*(p2mup4*F(1081)+p3mup4*F(1082
     -   )+p6mup4*F(1089))+p1mup3*(p1mup5*p2mup4*F(1081)+p2mup4*p2mup
     -   5*F(1081)+p1mup5*p3mup4*F(1082)+p2mup5*p3mup4*F(1082)+p2mup4
     -   *p3mup5*F(1082)-p2mup4*p6mup5*F(1083)+p1mup4*((p1mup5+p2mup5
     -   )*F(1081)+p3mup5*F(1082)-p6mup5*F(1083))+p3mup4*p3mup5*F(108
     -   4)-p3mup4*p6mup5*F(1085)+p1mup5*p6mup4*F(1089)+p2mup5*p6mup4
     -   *F(1089)+p3mup5*p6mup4*F(1090)-p6mup4*p6mup5*F(1091))+p1mup4
     -   *(p1mup5*p2mup3*F(1081)+p2mup3*p2mup5*F(1081)+p2mup3*p3mup5*
     -   F(1082)-p2mup3*p6mup5*F(1083)-p3mup3*p6mup5*F(1086)+p1mup5*p
     -   3mup3*F(1087)+p2mup5*p3mup3*F(1087)+p3mup3*p3mup5*F(1088)-p6
     -   mup3*p6mup5*F(1092)+p1mup5*p6mup3*F(1093)+p2mup5*p6mup3*F(10
     -   93)-p3mup5*p6mup3*F(1094)))
       Fa17(186) = p1mup5*p3mup2*p3mup3*p3mup4*F(1088)+p2mup5*p3mup2*p
     -   3mup3*p3mup4*F(1088)+p2mup3*p3mup2*(p2mup5*p3mup4*F(1082)+p2
     -   mup4*(p2mup5*F(1081)+p3mup5*F(1082)-p6mup5*F(1083))+p3mup4*p
     -   3mup5*F(1084)-p3mup4*p6mup5*F(1085)+p2mup5*p6mup4*F(1089)+p3
     -   mup5*p6mup4*F(1090)-p6mup4*p6mup5*F(1091))-p1mup5*p3mup2*p3m
     -   up3*p6mup4*F(1094)-p2mup5*p3mup2*p3mup3*p6mup4*F(1094)-p2mup
     -   4*p3mup2*(p3mup3*(p6mup5*F(1086)-(p1mup5+p2mup5)*F(1087)-p3m
     -   up5*F(1088))+p6mup3*(p6mup5*F(1092)-(p1mup5+p2mup5)*F(1093)+
     -   p3mup5*F(1094)))-p3mup2*p3mup3*p3mup4*p6mup5*F(1095)-p2mup4*
     -   p4mup2*p4mup3*p4mup5*F(1097)-p2mup3*p4mup2*p4mup4*p4mup5*F(1
     -   097)+p4mup2*p4mup5*(p1mup3*p1mup4*F(1096)-(p1mup4*p4mup3+p1m
     -   up3*p4mup4)*F(1097))+p4mup2*p4mup3*p4mup4*p6mup5*F(1098)-p1m
     -   up3*p1mup5*p4mup2*p4mup4*F(1099)+p3mup2*p3mup3*p3mup4*p3mup5
     -   *F(1100)-p3mup2*p3mup4*p6mup3*p6mup5*F(1101)-p3mup2*p3mup3*p
     -   6mup4*p6mup5*F(1101)+p1mup5*p3mup2*p3mup4*p6mup3*F(1102)+p2m
     -   up5*p3mup2*p3mup4*p6mup3*F(1102)-p3mup2*p3mup4*p3mup5*p6mup3
     -   *F(1103)-p3mup2*p3mup3*p3mup5*p6mup4*F(1103)-p4mup2*p4mup3*p
     -   4mup4*p4mup5*F(1104)-p1mup5*p4mup2*p4mup3*p4mup4*F(1105)-p2m
     -   up5*p4mup2*p4mup3*p4mup4*F(1105)-p4mup2*p4mup4*p4mup5*p6mup3
     -   *F(1106)-p1mup5*p3mup2*p6mup3*p6mup4*F(1107)-p2mup5*p3mup2*p
     -   6mup3*p6mup4*F(1107)-p3mup2*p6mup3*p6mup4*p6mup5*F(1108)-p4m
     -   up2*p4mup3*p4mup5*p6mup4*F(1109)+p3mup2*p3mup5*p6mup3*p6mup4
     -   *F(1110)
       Fa17(187) = p4mup2*(p2mup3*p2mup4*p4mup5*F(1096)-p1mup5*p2mup4*
     -   p4mup3*F(1099)-p2mup4*p2mup5*p4mup3*F(1099)-p1mup5*p2mup3*p4
     -   mup4*F(1099)-p2mup3*p2mup5*p4mup4*F(1099)-p2mup3*p4mup4*p6mu
     -   p5*F(1111)+p2mup4*p4mup3*p6mup5*F(1114)-p2mup4*p4mup5*p6mup3
     -   *F(1115)+p4mup4*p6mup3*p6mup5*F(1116)+p1mup4*(p2mup3*p4mup5*
     -   F(1096)-(p1mup5+p2mup5)*p4mup3*F(1099)+p1mup3*p1mup5*F(1112)
     -   +p1mup3*p2mup5*F(1112)+p2mup3*(p1mup5+p2mup5)*F(1112)-p1mup3
     -   *p6mup5*F(1113)+p4mup3*p6mup5*F(1114)-p6mup3*(p4mup5*F(1115)
     -   +p1mup5*F(1119)))-p1mup5*p4mup4*p6mup3*F(1120)-p2mup5*p4mup4
     -   *p6mup3*F(1120)-p2mup3*p4mup5*p6mup4*F(1120)+p1mup3*(-(p4mup
     -   4*(p2mup5*F(1099)+p6mup5*F(1111)))+p2mup4*(p4mup5*F(1096)+(p
     -   1mup5+p2mup5)*F(1112)-p6mup5*F(1113))-p6mup4*((p1mup5+p2mup5
     -   )*F(1117)-p6mup5*F(1118)+p4mup5*F(1120)))-p1mup5*p4mup3*p6mu
     -   p4*F(1121)-p2mup5*p4mup3*p6mup4*F(1121)+p4mup3*p6mup4*p6mup5
     -   *F(1122)+p4mup5*p6mup3*p6mup4*F(1123))
       Return
       End
