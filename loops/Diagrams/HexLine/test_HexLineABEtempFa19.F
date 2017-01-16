c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineAbe19(p1mup2,p1mup3,p1mup4,p1m
     -   up5,p2mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,
     -   p4mup2,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mu
     -   p2,p6mup3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,
     -   mup3mup5,mup4mup5,Fa19)
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
       Complex*16 Fa19(199:209)
       COMMON/test_FaFunctionsAbe/Fa
       Fa19(199) = p2mup2*(p2mup5*p4mup3*p6mup4*F(1058)-p2mup5*p4mup4*
     -   p6mup3*F(1070)-p4mup4*p6mup3*p6mup5*F(1074)-p4mup3*p6mup4*p6
     -   mup5*F(1078)+p4mup5*p6mup3*p6mup4*F(1079)+p1mup3*p1mup4*p2mu
     -   p5*F(1081)-p1mup3*p1mup4*p6mup5*F(1083)-p2mup4*(p2mup5*p4mup
     -   3*F(1023)-p4mup5*p6mup3*F(1056)+p4mup3*p6mup5*F(1069)-p1mup3
     -   *p2mup5*F(1081)+p1mup3*p6mup5*F(1083))+p2mup3*(p2mup5*p4mup4
     -   *F(1012)-p4mup4*p6mup5*F(1015)+p2mup4*p4mup5*F(1020)-p4mup5*
     -   p6mup4*F(1040)+p1mup4*p2mup5*F(1081)-p1mup4*p6mup5*F(1083))+
     -   p1mup3*p2mup5*p6mup4*F(1089)-p1mup3*p6mup4*p6mup5*F(1091)-p1
     -   mup4*p6mup3*p6mup5*F(1092)+p1mup4*p2mup5*p6mup3*F(1093))+p1m
     -   up5*(p2mup2*(p2mup3*p4mup4*F(1012)-p2mup4*p4mup3*F(1023)+p4m
     -   up3*p6mup4*F(1058)-p4mup4*p6mup3*F(1070)+p1mup4*p2mup3*F(108
     -   1)+p1mup3*p2mup4*F(1081)+p1mup3*p6mup4*F(1089)+p1mup4*p6mup3
     -   *F(1093))+p1mup2*(p1mup3*(p1mup4*F(1081)+p6mup4*F(1089))+p2m
     -   up3*(p2mup4*F(1081)+p6mup4*F(1089))+(p1mup4+p2mup4)*p6mup3*F
     -   (1093)-p6mup3*p6mup4*F(1107)))+p1mup2*(p2mup3*p2mup4*p2mup5*
     -   F(1081)-p2mup3*p2mup4*p6mup5*F(1083)+p2mup3*p2mup5*p6mup4*F(
     -   1089)-p2mup3*p6mup4*p6mup5*F(1091)+p1mup3*(p2mup4*p2mup5*F(1
     -   081)-p2mup4*p6mup5*F(1083)+p2mup5*p6mup4*F(1089)-p6mup4*p6mu
     -   p5*F(1091))-p2mup4*p6mup3*p6mup5*F(1092)+p2mup4*p2mup5*p6mup
     -   3*F(1093)+p1mup4*(p2mup3*(p2mup5*F(1081)-p6mup5*F(1083))-p6m
     -   up5*(p1mup3*F(1083)+p6mup3*F(1092))+p2mup5*p6mup3*F(1093))-p
     -   2mup5*p6mup3*p6mup4*F(1107)-p6mup3*p6mup4*p6mup5*F(1108))
       Fa19(200) = 32*p2mup2*(p2mup3*(-(p2mup4*p2mup5*F(1081))+p2mup4*
     -   p6mup5*F(1083)-p2mup5*p6mup4*F(1089)+p6mup4*p6mup5*F(1091))-
     -   p1mup5*(p2mup3*p2mup4*F(1081)+p2mup3*p6mup4*F(1089)+p2mup4*p
     -   6mup3*F(1093)-p6mup3*p6mup4*F(1107))+p6mup3*(p2mup4*p6mup5*F
     -   (1092)-p2mup4*p2mup5*F(1093)+p2mup5*p6mup4*F(1107)+p6mup4*p6
     -   mup5*F(1108)))+mup3mup5*(p1mup4*(p1mup2+p2mup2)*F(1138)+(p1m
     -   up2+p2mup2)*p2mup4*F(1138)+16*(p1mup4*p4mup2*s12*F(913)+p2mu
     -   p4*p4mup2*s12*F(913)+p3mup4*p4mup2*s12*F(914)+(p1mup2+p2mup2
     -   )*(p4mup4*F(1129)+p6mup4*F(1132))+(p1mup2+p2mup2)*p3mup4*F(1
     -   137)+p4mup2*s12*(p4mup4*F(915)+p6mup4*F(1139))))
       Fa19(201) = mup4mup5*(16*p3mup2*p4mup3*s12*F(816)-(p1mup2+p2mup
     -   2)*(16*p4mup3*F(999)-(p1mup3+p2mup3)*F(1140)-p3mup3*F(1141))
     -   )
       Fa19(202) = mup4mup5*((p1mup2+p2mup2)*p6mup3*F(1000)+p3mup2*s12
     -   *((p1mup3+p2mup3)*F(814)-p3mup3*F(815)+p6mup3*F(1142)))
       Fa19(203) = mup3mup4*(16*p3mup2*p3mup5*s12*F(954)+(p2mup2*p2mup
     -   5+p1mup2*(p1mup5+p2mup5))*F(1143)+p1mup2*p3mup5*F(1144)+p1mu
     -   p2*p4mup5*F(1145))
       Fa19(204) = mup3mup4*(p1mup5*p2mup2*F(1143)+p2mup2*p3mup5*F(114
     -   4)+p2mup2*p4mup5*F(1145)+16*(p1mup5*p3mup2*s12*F(950)+p2mup5
     -   *p3mup2*s12*F(950)+p3mup2*p4mup5*s12*F(991)+(p1mup2+p2mup2)*
     -   p6mup5*F(1010)+p3mup2*p6mup5*s12*F(1146)))
       Fa19(205) = p1mup2*(p1mup4*p3mup3*p4mup5*F(1018)+p2mup4*p3mup3*
     -   p4mup5*F(1018)-p3mup3*p4mup4*p6mup5*F(1019)-p3mup4*p4mup3*p4
     -   mup5*F(1025)-p1mup3*p3mup4*p4mup5*F(1026)-p2mup3*p3mup4*p4mu
     -   p5*F(1026)-p1mup3*p3mup5*p4mup4*F(1027)-p2mup3*p3mup5*p4mup4
     -   *F(1027)-p1mup4*p3mup5*p4mup3*F(1028)-p2mup4*p3mup5*p4mup3*F
     -   (1028)-p1mup5*p3mup3*p4mup4*F(1030)-p2mup5*p3mup3*p4mup4*F(1
     -   030)-p1mup5*p3mup4*p4mup3*F(1031)-p2mup5*p3mup4*p4mup3*F(103
     -   1)-p3mup3*p3mup4*p4mup5*F(1032)+p3mup3*p4mup5*p6mup4*F(1036)
     -   -p3mup4*p3mup5*p4mup3*F(1037)+p3mup4*p4mup3*p6mup5*F(1038)-p
     -   3mup4*p4mup5*p6mup3*F(1039)+p3mup3*p4mup4*p4mup5*F(1044)+p3m
     -   up3*p3mup5*p4mup4*F(1045)-p3mup5*p4mup3*p4mup4*F(1046)-p3mup
     -   5*p4mup3*p6mup4*F(1047)-p3mup5*p4mup4*p6mup3*F(1048)+(p1mup4
     -   *p1mup5*p2mup3+p1mup3*p1mup5*p2mup4+p1mup3*p1mup4*p2mup5)*F(
     -   1112))-p2mup2*(p3mup4*p4mup5*(p4mup3*F(1025)+p1mup3*F(1026))
     -   +p4mup4*(p1mup3*p3mup5*F(1027)-p3mup3*p4mup5*F(1044)+p3mup5*
     -   p4mup3*F(1046))-p1mup4*(p3mup3*p4mup5*F(1018)-p3mup5*p4mup3*
     -   F(1028)+p1mup3*p1mup5*F(1112)))
       Fa19(206) = p2mup2*(p2mup3*p3mup4*p4mup5*F(1026)+p2mup3*p3mup5*
     -   p4mup4*F(1027)+p2mup4*(-(p3mup3*p4mup5*F(1018))+p3mup5*p4mup
     -   3*F(1028))+(p1mup5+p2mup5)*p3mup4*p4mup3*F(1031)+p3mup4*p3mu
     -   p5*p4mup3*F(1037)-p3mup4*p4mup3*p6mup5*F(1038)+p3mup4*p4mup5
     -   *p6mup3*F(1039)+p3mup3*(p4mup5*(p3mup4*F(1032)-p6mup4*F(1036
     -   ))+p4mup4*(p6mup5*F(1019)+(p1mup5+p2mup5)*F(1030)-p3mup5*F(1
     -   045)))+p3mup5*p4mup3*p6mup4*F(1047)+p3mup5*p4mup4*p6mup3*F(1
     -   048))
       Fa19(207) = p2mup2*(p1mup3*(-((p1mup5+p2mup5)*p3mup4*F(1013))-p
     -   3mup4*p6mup5*F(1022)+p3mup5*((p1mup4+p2mup4)*F(1011)-p3mup4*
     -   F(1041)-p6mup4*F(1042)))+p1mup4*(p2mup3*p3mup5*F(1011)+(p1mu
     -   p5+p2mup5)*p3mup3*F(1024)+p3mup3*p3mup5*F(1050)-p3mup3*p6mup
     -   5*F(1054)-p3mup5*p6mup3*F(1057))-p1mup5*(p2mup3*p3mup4*F(101
     -   3)-p2mup4*p3mup3*F(1024)+p3mup3*p3mup4*F(1052)+p3mup3*p6mup4
     -   *F(1059)+p3mup4*p6mup3*F(1060)))+p1mup2*(p2mup3*p2mup4*p3mup
     -   5*F(1011)-p1mup5*p2mup3*p3mup4*F(1013)-p2mup3*p2mup5*p3mup4*
     -   F(1013)-p2mup3*p3mup4*p6mup5*F(1022)+p1mup5*p2mup4*p3mup3*F(
     -   1024)+p2mup4*p2mup5*p3mup3*F(1024)-p2mup3*p3mup4*p3mup5*F(10
     -   41)-p2mup3*p3mup5*p6mup4*F(1042)+p1mup3*(-((p1mup5+p2mup5)*p
     -   3mup4*F(1013))-p3mup4*p6mup5*F(1022)+p3mup5*((p1mup4+p2mup4)
     -   *F(1011)-p3mup4*F(1041)-p6mup4*F(1042)))+p2mup4*p3mup3*p3mup
     -   5*F(1050)-p1mup5*p3mup3*p3mup4*F(1052)-p2mup5*p3mup3*p3mup4*
     -   F(1052)-p2mup4*p3mup3*p6mup5*F(1054)-p3mup3*p3mup4*p6mup5*F(
     -   1055)-p2mup4*p3mup5*p6mup3*F(1057)+p1mup4*(p2mup3*p3mup5*F(1
     -   011)+(p1mup5+p2mup5)*p3mup3*F(1024)+p3mup3*p3mup5*F(1050)-p3
     -   mup3*p6mup5*F(1054)-p3mup5*p6mup3*F(1057))-p1mup5*p3mup3*p6m
     -   up4*F(1059)-p2mup5*p3mup3*p6mup4*F(1059)-p1mup5*p3mup4*p6mup
     -   3*F(1060)-p2mup5*p3mup4*p6mup3*F(1060)-p3mup3*p3mup4*p3mup5*
     -   F(1061)+p3mup3*p6mup4*p6mup5*F(1063)-p3mup4*p3mup5*p6mup3*F(
     -   1064)+p3mup4*p6mup3*p6mup5*F(1065)-p3mup3*p3mup5*p6mup4*F(10
     -   75)-p3mup5*p6mup3*p6mup4*F(1080))
       Fa19(208) = p2mup2*(p2mup4*p2mup5*p3mup3*F(1024)+p2mup3*(p2mup4
     -   *p3mup5*F(1011)-p2mup5*p3mup4*F(1013)-p3mup4*p6mup5*F(1022)-
     -   p3mup4*p3mup5*F(1041)-p3mup5*p6mup4*F(1042))+p2mup4*p3mup3*p
     -   3mup5*F(1050)-p2mup5*p3mup3*p3mup4*F(1052)-p2mup4*p3mup3*p6m
     -   up5*F(1054)-p3mup3*p3mup4*p6mup5*F(1055)-p2mup4*p3mup5*p6mup
     -   3*F(1057)-p2mup5*p3mup3*p6mup4*F(1059)-p2mup5*p3mup4*p6mup3*
     -   F(1060)-p3mup3*p3mup4*p3mup5*F(1061)+p3mup3*p6mup4*p6mup5*F(
     -   1063)-p3mup4*p3mup5*p6mup3*F(1064)+p3mup4*p6mup3*p6mup5*F(10
     -   65)-p3mup3*p3mup5*p6mup4*F(1075)-p3mup5*p6mup3*p6mup4*F(1080
     -   )-p4mup3*p4mup4*p4mup5*F(1104))-p1mup2*(p1mup4*p4mup3*p4mup5
     -   *F(1097)+p2mup4*p4mup3*p4mup5*F(1097)+p4mup4*(p1mup3*p4mup5*
     -   F(1097)+p2mup3*p4mup5*F(1097)-p4mup3*p6mup5*F(1098)+p4mup3*p
     -   4mup5*F(1104)+p1mup5*p4mup3*F(1105)+p2mup5*p4mup3*F(1105)+p4
     -   mup5*p6mup3*F(1106))+p4mup3*p4mup5*p6mup4*F(1109))
       Fa19(209) = p2mup2*(p2mup4*p4mup3*p4mup5*F(1097)+p2mup3*p4mup4*
     -   p4mup5*F(1097)-p4mup3*p4mup4*p6mup5*F(1098)+p1mup5*p4mup3*p4
     -   mup4*F(1105)+p2mup5*p4mup3*p4mup4*F(1105)+p4mup4*p4mup5*p6mu
     -   p3*F(1106)+p4mup3*p4mup5*p6mup4*F(1109)+p1mup4*(-(p2mup3*p4m
     -   up5*F(1096))+p4mup3*(p4mup5*F(1097)+(p1mup5+p2mup5)*F(1099)-
     -   p6mup5*F(1114))+p4mup5*p6mup3*F(1115))+p1mup3*(-(p1mup4*p4mu
     -   p5*F(1096))-p2mup4*p4mup5*F(1096)+p4mup4*(p4mup5*F(1097)+(p1
     -   mup5+p2mup5)*F(1099)+p6mup5*F(1111))+p4mup5*p6mup4*F(1120)))
     -   +p1mup2*(-(p2mup3*p2mup4*p4mup5*F(1096))+p1mup5*p2mup4*p4mup
     -   3*F(1099)+p2mup4*p2mup5*p4mup3*F(1099)+p1mup5*p2mup3*p4mup4*
     -   F(1099)+p2mup3*p2mup5*p4mup4*F(1099)+p2mup3*p4mup4*p6mup5*F(
     -   1111)-p2mup4*p4mup3*p6mup5*F(1114)+p2mup4*p4mup5*p6mup3*F(11
     -   15)+p1mup4*(-(p2mup3*p4mup5*F(1096))+(p1mup5+p2mup5)*p4mup3*
     -   F(1099)-p4mup3*p6mup5*F(1114)+p4mup5*p6mup3*F(1115))-p4mup4*
     -   p6mup3*p6mup5*F(1116)+p1mup5*p4mup4*p6mup3*F(1120)+p2mup5*p4
     -   mup4*p6mup3*F(1120)+p2mup3*p4mup5*p6mup4*F(1120)+p1mup3*(-((
     -   p1mup4+p2mup4)*p4mup5*F(1096))+(p1mup5+p2mup5)*p4mup4*F(1099
     -   )+p4mup4*p6mup5*F(1111)+p4mup5*p6mup4*F(1120))+p1mup5*p4mup3
     -   *p6mup4*F(1121)+p2mup5*p4mup3*p6mup4*F(1121)-p4mup3*p6mup4*p
     -   6mup5*F(1122)+p4mup5*p6mup3*p6mup4*F(1147))
       Return
       End
