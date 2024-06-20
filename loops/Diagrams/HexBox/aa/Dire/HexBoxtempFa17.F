c************** Calling the Fa functions*************************
       subroutine HexBoxFa17(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4
     &   ,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa17)
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
       Complex*16 Fa17(193:204)
       COMMON/HexBoxFaFunctions/Fa
       Fa17(193) = 2*(p1mup4*F(976)+p6mup4*F(977)+p4mup4*F(978)-p5mup4
     &   *F(979))+p2mup4*F(980)
       Fa17(194) = p5mup3*F(984)+2*(p4mup3*F(981)+p1mup3*F(982)+p6mup3
     &   *F(983)-p2mup3*F(985))
       Fa17(195) = mup3mup4*F(986)+4*(p4mup3*(p5mup4*F(586)-p2mup4*F(5
     &   92)+p4mup4*F(931)+p6mup4*F(933)-p1mup4*F(943))-p1mup3*(p4mup
     &   4*F(990)+p5mup4*F(991)+p6mup4*F(992))-p2mup3*(p4mup4*F(987)+
     &   p5mup4*F(988)+p6mup4*F(989)-p1mup4*F(993))+p5mup3*(p4mup4*F(
     &   995)+p5mup4*F(997)+p6mup4*F(998))+p6mup3*(p4mup4*F(996)+p5mu
     &   p4*F(999)+p6mup4*F(1000))-p1mup4*(p1mup3*F(994)+p5mup3*F(100
     &   2)+p6mup3*F(1003))-p2mup4*(p6mup3*F(626)-p1mup3*F(993)-p2mup
     &   3*F(1001)+p5mup3*F(1004)))
       Fa17(196) = -2*(mup3mup4*F(1005)-2*(p6mup3*p6mup4*F(662)-p1mup4
     &   *p6mup3*F(666)+p4mup3*p6mup4*F(958)-p1mup4*p4mup3*F(961)+p1m
     &   up4*p2mup3*F(1006)-p1mup3*p6mup4*F(1009)-p1mup3*p1mup4*F(101
     &   0)-p2mup3*p6mup4*F(1013)+p5mup4*(p6mup3*F(661)+p4mup3*F(957)
     &   -p1mup3*F(1008)-p2mup3*F(1012)+p5mup3*F(1015))+p5mup3*p6mup4
     &   *F(1016)-p1mup4*p5mup3*F(1017)))+4*(p2mup4*(p4mup3*F(588)-p2
     &   mup3*F(593)+p6mup3*F(619)-p1mup3*F(987)+p5mup3*F(1011))+p4mu
     &   p4*(p6mup3*F(659)+p4mup3*F(956)-p1mup3*F(1007)+p5mup3*F(1014
     &   )-p2mup3*F(1018)))
       Fa17(197) = 2*mup3mup4*F(1019)+p5mup3*p6mup4*F(1026)+4*(p2mup4*
     &   p6mup3*F(622)+p2mup3*(p4mup4*F(619)+p5mup4*F(621)+p6mup4*F(6
     &   22)-p1mup4*F(626)-p2mup4*F(629))+p4mup3*(p5mup4*F(540)+p2mup
     &   4*F(590)+p1mup4*F(933)+p4mup4*F(968)+p6mup4*F(970))-p1mup3*p
     &   2mup4*F(989)-p1mup3*p1mup4*F(992)+p1mup4*p5mup3*F(998)+p1mup
     &   4*p6mup3*F(1000)+p1mup3*p4mup4*F(1020)+p1mup3*p5mup4*F(1021)
     &   +p1mup3*p6mup4*F(1022)+p4mup4*p5mup3*F(1023)+p4mup4*p6mup3*F
     &   (1024)+p5mup3*p5mup4*F(1025)+p5mup4*p6mup3*F(1027)+p6mup3*p6
     &   mup4*F(1028)+p2mup4*p5mup3*F(1029))
       Fa17(198) = p1mup4*F(1033)-p2mup4*F(1034)+p6mup4*F(1035)+p5mup4
     &   *F(1036)-p4mup4*F(1037)
       Fa17(199) = p1mup3*F(1038)-p2mup3*F(1039)+p6mup3*F(1040)+p4mup3
     &   *F(1041)+p5mup3*F(1042)
       Fa17(200) = p6mup4*F(1045)+p4mup4*F(1046)+p5mup4*F(1047)+p1mup4
     &   *F(1048)+p2mup4*F(1049)
       Fa17(201) = p2mup3*F(1050)+p1mup3*F(1051)+p6mup3*F(1052)+p4mup3
     &   *F(1053)+p5mup3*F(1054)
       Fa17(202) = -(p2mup4*F(1055))+p6mup4*F(1056)+p4mup4*F(1057)+p5m
     &   up4*F(1058)-p1mup4*F(1059)
       Fa17(203) = p6mup4*F(1060)+p4mup4*F(1061)+p5mup4*F(1062)-p2mup4
     &   *F(1063)-p1mup4*F(1064)
       Fa17(204) = p1mup4*F(1065)-p2mup4*F(1066)-p6mup4*F(1067)-p5mup4
     &   *F(1068)+p4mup4*F(1069)
       Return
       End
