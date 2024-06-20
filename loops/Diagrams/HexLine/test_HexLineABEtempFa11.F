c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineAbe11(p1mup2,p1mup3,p1mup4,p1m
     -   up5,p2mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,
     -   p4mup2,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mu
     -   p2,p6mup3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,
     -   mup3mup5,mup4mup5,Fa11)
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
       Complex*16 Fa11(111:121)
       COMMON/test_FaFunctionsAbe/Fa
       Fa11(111) = p2mup2*(p2mup4*p3mup3*p4mup5*F(226)+p2mup3*(p3mup4*
     -   p4mup5*F(226)+p3mup5*p4mup4*F(231))+p3mup3*(p3mup4*p4mup5*F(
     -   226)-p4mup4*p6mup5*F(227)+(p1mup5+p2mup5+p3mup5)*p4mup4*F(23
     -   1)-p4mup5*p6mup4*F(233)))
       Fa11(112) = p2mup2*(p1mup5*(p2mup4*p3mup3+(p1mup3+p2mup3+p3mup3
     -   )*p3mup4)*F(228)+p1mup4*((p1mup5+p2mup5)*p3mup3+(p1mup3+p2mu
     -   p3+p3mup3)*p3mup5)*F(228)-p1mup4*p3mup3*p6mup5*F(229)+p1mup5
     -   *p3mup3*p6mup4*F(232)+p1mup3*(p2mup5*p3mup4*F(228)+p2mup4*p3
     -   mup5*F(228)+p3mup4*p3mup5*F(228)-p3mup4*p6mup5*F(229)+p3mup5
     -   *p6mup4*F(232)))+p1mup2*(p2mup4*p2mup5*p3mup3*F(228)+p1mup3*
     -   p2mup5*p3mup4*F(228)+p2mup3*p2mup5*p3mup4*F(228)+p2mup5*p3mu
     -   p3*p3mup4*F(228)+p1mup5*(p2mup4*p3mup3+(p1mup3+p2mup3+p3mup3
     -   )*p3mup4)*F(228)+p1mup3*p2mup4*p3mup5*F(228)+p2mup3*p2mup4*p
     -   3mup5*F(228)+p2mup4*p3mup3*p3mup5*F(228)+p1mup3*p3mup4*p3mup
     -   5*F(228)+p2mup3*p3mup4*p3mup5*F(228)+p3mup3*p3mup4*p3mup5*F(
     -   228)+p1mup4*((p1mup5+p2mup5)*p3mup3+(p1mup3+p2mup3+p3mup3)*p
     -   3mup5)*F(228)-p1mup4*p3mup3*p6mup5*F(229)-p2mup4*p3mup3*p6mu
     -   p5*F(229)-p1mup3*p3mup4*p6mup5*F(229)-p2mup3*p3mup4*p6mup5*F
     -   (229)-p3mup3*p3mup4*p6mup5*F(229)+p1mup5*p3mup3*p6mup4*F(232
     -   )+p2mup5*p3mup3*p6mup4*F(232)+p1mup3*p3mup5*p6mup4*F(232)+p2
     -   mup3*p3mup5*p6mup4*F(232)+p3mup3*p3mup5*p6mup4*F(232)-p3mup3
     -   *p6mup4*p6mup5*F(234))
       Fa11(113) = p1mup2*(p1mup3+p2mup3)*p4mup4*p4mup5*F(230)+p2mup2*
     -   (p2mup4*(p2mup5*p3mup3+(p2mup3+p3mup3)*p3mup5)*F(228)-p2mup4
     -   *p3mup3*p6mup5*F(229)+p2mup3*(p2mup5*p3mup4*F(228)+p3mup4*p3
     -   mup5*F(228)-p3mup4*p6mup5*F(229)+p3mup5*p6mup4*F(232))+p3mup
     -   3*(p2mup5*p3mup4*F(228)+p3mup4*p3mup5*F(228)-p3mup4*p6mup5*F
     -   (229)+p2mup5*p6mup4*F(232)+p3mup5*p6mup4*F(232)-p6mup4*p6mup
     -   5*F(234)))
       Fa11(114) = p1mup2*(p1mup3+p2mup3)*((p1mup4+p2mup4)*p4mup5*F(22
     -   6)-p4mup4*p6mup5*F(227)+(p1mup5+p2mup5)*p4mup4*F(231)-p4mup5
     -   *p6mup4*F(233))+p2mup2*(p2mup3*p4mup5*(p1mup4*F(226)+p4mup4*
     -   F(230))+p1mup3*(p1mup4*p4mup5*F(226)+p2mup4*p4mup5*F(226)-p4
     -   mup4*p6mup5*F(227)+p4mup4*p4mup5*F(230)+p1mup5*p4mup4*F(231)
     -   +p2mup5*p4mup4*F(231)-p4mup5*p6mup4*F(233)))
       Fa11(115) = p1mup2*p1mup3*p1mup4*p1mup5*F(228)+p1mup2*(p1mup5*p
     -   2mup3*p2mup4*F(228)+p1mup4*p2mup3*p2mup5*F(228)-p1mup4*(p1mu
     -   p3+p2mup3)*p6mup5*F(229)+p1mup5*(p1mup3+p2mup3)*p6mup4*F(232
     -   )+(p1mup3+p2mup3)*(p2mup4*p2mup5*F(228)-p2mup4*p6mup5*F(229)
     -   +p2mup5*p6mup4*F(232)-p6mup4*p6mup5*F(234)))+p2mup2*(p2mup3*
     -   (p2mup4*p4mup5*F(226)+p1mup4*(p1mup5+p2mup5)*F(228)-p6mup5*(
     -   p4mup4*F(227)+p1mup4*F(229))+(p1mup5+p2mup5)*p4mup4*F(231)-p
     -   4mup5*p6mup4*F(233))+p1mup3*(p1mup5*p2mup4*F(228)+p1mup4*p2m
     -   up5*F(228)+p2mup4*p2mup5*F(228)-p1mup4*p6mup5*F(229)-p2mup4*
     -   p6mup5*F(229)+p1mup5*p6mup4*F(232)+p2mup5*p6mup4*F(232)-p6mu
     -   p4*p6mup5*F(234)))
       Fa11(116) = p2mup2*p2mup3*(p1mup5*p2mup4*F(228)+p2mup4*p2mup5*F
     -   (228)-p2mup4*p6mup5*F(229)+p1mup5*p6mup4*F(232)+p2mup5*p6mup
     -   4*F(232)-p6mup4*p6mup5*F(234))
       Fa11(117) = p4mup3*(mup4mup5*F(225)+2*(p3mup4*p4mup5*F(226)-p4m
     -   up4*p6mup5*F(227)+p1mup5*p3mup4*F(228)+p2mup5*p3mup4*F(228)+
     -   p3mup4*p3mup5*F(228)-p3mup4*p6mup5*F(229)+p1mup4*(p4mup5*F(2
     -   26)+(p1mup5+p2mup5+p3mup5)*F(228)-p6mup5*F(229))+p2mup4*(p4m
     -   up5*F(226)+(p1mup5+p2mup5+p3mup5)*F(228)-p6mup5*F(229))+p4mu
     -   p4*p4mup5*F(230)+p1mup5*p4mup4*F(231)+p2mup5*p4mup4*F(231)+p
     -   3mup5*p4mup4*F(231)+p1mup5*p6mup4*F(232)+p2mup5*p6mup4*F(232
     -   )+p3mup5*p6mup4*F(232)-p4mup5*p6mup4*F(233)-p6mup4*p6mup5*F(
     -   234)))
       Fa11(118) = mup4mup5*p4mup2*s12*F(225)-mup4mup5*(p1mup2+p2mup2)
     -   *s123*F(235)+2*(p3mup4*p4mup2*p4mup5*s12*F(226)-p4mup2*p4mup
     -   4*p6mup5*s12*F(227)+p1mup5*p3mup4*p4mup2*s12*F(228)+p2mup5*p
     -   3mup4*p4mup2*s12*F(228)+p3mup4*p3mup5*p4mup2*s12*F(228)-p3mu
     -   p4*p4mup2*p6mup5*s12*F(229)+p4mup2*p4mup4*p4mup5*s12*F(230)+
     -   p1mup5*p4mup2*p4mup4*s12*F(231)+p2mup5*p4mup2*p4mup4*s12*F(2
     -   31)+p3mup5*p4mup2*p4mup4*s12*F(231)+p1mup5*p4mup2*p6mup4*s12
     -   *F(232)+p2mup5*p4mup2*p6mup4*s12*F(232)+p3mup5*p4mup2*p6mup4
     -   *s12*F(232)-p4mup2*p4mup5*p6mup4*s12*F(233)-p4mup2*p6mup4*p6
     -   mup5*s12*F(234)-p1mup2*p3mup4*p6mup5*s123*F(236)-p2mup2*p3mu
     -   p4*p6mup5*s123*F(236)-p1mup2*p4mup4*p6mup5*s123*F(237)-p2mup
     -   2*p4mup4*p6mup5*s123*F(237)+p1mup2*p1mup5*p3mup4*s123*F(238)
     -   +p1mup5*p2mup2*p3mup4*s123*F(238)+p1mup2*p2mup5*p3mup4*s123*
     -   F(238)+p2mup2*p2mup5*p3mup4*s123*F(238)+p1mup2*p3mup4*p3mup5
     -   *s123*F(238)+p2mup2*p3mup4*p3mup5*s123*F(238)+p1mup2*p3mup4*
     -   p4mup5*s123*F(239)+p2mup2*p3mup4*p4mup5*s123*F(239)+p1mup4*(
     -   p4mup2*s12*(p4mup5*F(226)+(p1mup5+p2mup5+p3mup5)*F(228)-p6mu
     -   p5*F(229))+(p1mup2+p2mup2)*s123*(-(p6mup5*F(236))+(p1mup5+p2
     -   mup5+p3mup5)*F(238)+p4mup5*F(239)))+p2mup4*(p4mup2*s12*(p4mu
     -   p5*F(226)+(p1mup5+p2mup5+p3mup5)*F(228)-p6mup5*F(229))+(p1mu
     -   p2+p2mup2)*s123*(-(p6mup5*F(236))+(p1mup5+p2mup5+p3mup5)*F(2
     -   38)+p4mup5*F(239)))+p1mup2*p6mup4*p6mup5*s123*F(240)+p2mup2*
     -   p6mup4*p6mup5*s123*F(240)+p1mup2*p1mup5*p6mup4*s123*F(241)+p
     -   1mup5*p2mup2*p6mup4*s123*F(241)+p1mup2*p2mup5*p6mup4*s123*F(
     -   241)+p2mup2*p2mup5*p6mup4*s123*F(241)+p1mup2*p3mup5*p6mup4*s
     -   123*F(241)+p2mup2*p3mup5*p6mup4*s123*F(241)+p1mup2*p4mup5*p6
     -   mup4*s123*F(242)+p2mup2*p4mup5*p6mup4*s123*F(242)+p1mup2*p4m
     -   up4*p4mup5*s123*F(243)+p2mup2*p4mup4*p4mup5*s123*F(243)+(p1m
     -   up2+p2mup2)*(p1mup5+p2mup5+p3mup5)*p4mup4*s123*F(244))
       Fa11(119) = 2*(p1mup4*(-(p6mup5*F(236))+(p1mup5+p2mup5+p3mup5)*
     -   F(238)+p4mup5*F(239))+p2mup4*(-(p6mup5*F(236))+(p1mup5+p2mup
     -   5+p3mup5)*F(238)+p4mup5*F(239))+p3mup4*(-(p6mup5*F(236))+(p1
     -   mup5+p2mup5+p3mup5)*F(238)+p4mup5*F(239))+p6mup4*(p6mup5*F(2
     -   40)+(p1mup5+p2mup5+p3mup5)*F(241)+p4mup5*F(242))+p4mup4*(-(p
     -   6mup5*F(237))+p4mup5*F(243)+(p1mup5+p2mup5+p3mup5)*F(244)))+
     -   mup4mup5*F(245)
       Fa11(120) = p4mup2*(mup4mup5*F(225)+2*(p3mup4*p4mup5*F(226)-p4m
     -   up4*p6mup5*F(227)+p1mup5*p3mup4*F(228)+p2mup5*p3mup4*F(228)+
     -   p3mup4*p3mup5*F(228)-p3mup4*p6mup5*F(229)+p1mup4*(p4mup5*F(2
     -   26)+(p1mup5+p2mup5+p3mup5)*F(228)-p6mup5*F(229))+p2mup4*(p4m
     -   up5*F(226)+(p1mup5+p2mup5+p3mup5)*F(228)-p6mup5*F(229))+p4mu
     -   p4*p4mup5*F(230)+p1mup5*p4mup4*F(231)+p2mup5*p4mup4*F(231)+p
     -   3mup5*p4mup4*F(231)+p1mup5*p6mup4*F(232)+p2mup5*p6mup4*F(232
     -   )+p3mup5*p6mup4*F(232)-p4mup5*p6mup4*F(233)-p6mup4*p6mup5*F(
     -   234)))
       Fa11(121) = 4*mup4mup5*F(225)+8*(p3mup4*p4mup5*F(226)-p4mup4*p6
     -   mup5*F(227)+p1mup5*p3mup4*F(228)+p2mup5*p3mup4*F(228)+p3mup4
     -   *p3mup5*F(228)-p3mup4*p6mup5*F(229)+p1mup4*(p4mup5*F(226)+(p
     -   1mup5+p2mup5+p3mup5)*F(228)-p6mup5*F(229))+p2mup4*(p4mup5*F(
     -   226)+(p1mup5+p2mup5+p3mup5)*F(228)-p6mup5*F(229))+p4mup4*p4m
     -   up5*F(230)+p1mup5*p4mup4*F(231)+p2mup5*p4mup4*F(231)+p3mup5*
     -   p4mup4*F(231)+p1mup5*p6mup4*F(232)+p2mup5*p6mup4*F(232)+p3mu
     -   p5*p6mup4*F(232)-p4mup5*p6mup4*F(233)-p6mup4*p6mup5*F(234))
       Return
       End
