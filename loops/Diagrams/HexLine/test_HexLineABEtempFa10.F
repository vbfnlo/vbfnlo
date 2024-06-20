c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineAbe10(p1mup2,p1mup3,p1mup4,p1m
     -   up5,p2mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,
     -   p4mup2,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mu
     -   p2,p6mup3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,
     -   mup3mup5,mup4mup5,Fa10)
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
       Complex*16 Fa10(100:110)
       COMMON/test_FaFunctionsAbe/Fa
       Fa10(100) = -((p1mup2+p2mup2)*(p4mup5*F(217)-p6mup5*F(218)-(p1m
     -   up5+p2mup5+p3mup5)*F(221)))
       Fa10(101) = -(p4mup5*F(217))+p6mup5*F(218)+(p1mup5+p2mup5+p3mup
     -   5)*F(221)
       Fa10(102) = mup4mup5*(p1mup3+p2mup3+p3mup3)*p4mup2-mup4mup5*p3m
     -   up2*p4mup3
       Fa10(103) = p1mup3*p3mup4*p4mup2*p4mup5*F(226)+p2mup3*p3mup4*p4
     -   mup2*p4mup5*F(226)+p3mup3*p3mup4*p4mup2*p4mup5*F(226)-p3mup2
     -   *p3mup4*p4mup3*p4mup5*F(226)+p1mup4*(p3mup3*p4mup2-p3mup2*p4
     -   mup3)*p4mup5*F(226)+p2mup4*(p3mup3*p4mup2-p3mup2*p4mup3)*p4m
     -   up5*F(226)-p3mup3*p4mup2*p4mup4*p6mup5*F(227)+p3mup2*p4mup3*
     -   p4mup4*p6mup5*F(227)+p1mup3*p1mup5*p3mup4*p4mup2*F(228)+p1mu
     -   p3*p2mup5*p3mup4*p4mup2*F(228)+p1mup3*p2mup4*p3mup5*p4mup2*F
     -   (228)+p1mup3*p3mup4*p3mup5*p4mup2*F(228)+p1mup4*((p1mup5*p3m
     -   up3+(p1mup3+p2mup3)*p3mup5)*p4mup2-(p1mup5+p2mup5)*p3mup2*p4
     -   mup3)*F(228)-p1mup3*p3mup4*p4mup2*p6mup5*F(229)+p3mup3*p4mup
     -   2*p4mup4*p4mup5*F(230)-p3mup2*p4mup3*p4mup4*p4mup5*F(230)+p1
     -   mup5*p3mup3*p4mup2*p4mup4*F(231)+p2mup5*p3mup3*p4mup2*p4mup4
     -   *F(231)+p1mup3*p3mup5*p4mup2*p4mup4*F(231)+p2mup3*p3mup5*p4m
     -   up2*p4mup4*F(231)+p3mup3*p3mup5*p4mup2*p4mup4*F(231)-p1mup5*
     -   p3mup2*p4mup3*p4mup4*F(231)-p2mup5*p3mup2*p4mup3*p4mup4*F(23
     -   1)-p3mup2*p3mup5*p4mup3*p4mup4*F(231)+p1mup3*p3mup5*p4mup2*p
     -   6mup4*F(232)-(p3mup3*p4mup2-p3mup2*p4mup3)*p4mup5*p6mup4*F(2
     -   33)
       Fa10(104) = p2mup4*p2mup5*p3mup3*p4mup2*F(228)+p2mup3*p2mup5*p3
     -   mup4*p4mup2*F(228)+p2mup5*p3mup3*p3mup4*p4mup2*F(228)+p2mup3
     -   *p2mup4*p3mup5*p4mup2*F(228)+p2mup4*p3mup3*p3mup5*p4mup2*F(2
     -   28)+p2mup3*p3mup4*p3mup5*p4mup2*F(228)+p3mup3*p3mup4*p3mup5*
     -   p4mup2*F(228)-p2mup4*p2mup5*p3mup2*p4mup3*F(228)-p2mup5*p3mu
     -   p2*p3mup4*p4mup3*F(228)-p2mup4*p3mup2*p3mup5*p4mup3*F(228)-p
     -   3mup2*p3mup4*p3mup5*p4mup3*F(228)+p1mup5*((p2mup4*p3mup3+(p2
     -   mup3+p3mup3)*p3mup4)*p4mup2-p3mup2*(p2mup4+p3mup4)*p4mup3)*F
     -   (228)-p2mup4*p3mup3*p4mup2*p6mup5*F(229)-p2mup3*p3mup4*p4mup
     -   2*p6mup5*F(229)-p3mup3*p3mup4*p4mup2*p6mup5*F(229)+p2mup4*p3
     -   mup2*p4mup3*p6mup5*F(229)+p3mup2*p3mup4*p4mup3*p6mup5*F(229)
     -   +p1mup4*(p2mup5*p3mup3*p4mup2*F(228)+(p3mup3*p4mup2-p3mup2*p
     -   4mup3)*(p3mup5*F(228)-p6mup5*F(229)))+p2mup5*p3mup3*p4mup2*p
     -   6mup4*F(232)+p2mup3*p3mup5*p4mup2*p6mup4*F(232)+p3mup3*p3mup
     -   5*p4mup2*p6mup4*F(232)-p2mup5*p3mup2*p4mup3*p6mup4*F(232)-p3
     -   mup2*p3mup5*p4mup3*p6mup4*F(232)+p1mup5*(p3mup3*p4mup2-p3mup
     -   2*p4mup3)*p6mup4*F(232)-(p3mup3*p4mup2-p3mup2*p4mup3)*p6mup4
     -   *p6mup5*F(234)
       Fa10(105) = p4mup2*(p2mup3*p4mup4*p4mup5*F(230)+p1mup3*(p1mup4*
     -   p4mup5*F(226)+p4mup4*p4mup5*F(230)+p1mup5*p4mup4*F(231)))
       Fa10(106) = p4mup2*(p1mup4*(p2mup3*p4mup5*F(226)+p1mup3*(p1mup5
     -   +p2mup5)*F(228)+p2mup3*(p1mup5+p2mup5)*F(228)-p1mup3*p6mup5*
     -   F(229))+p2mup3*(p2mup4*p4mup5*F(226)-p4mup4*p6mup5*F(227)+(p
     -   1mup5+p2mup5)*p4mup4*F(231)-p4mup5*p6mup4*F(233))+p1mup3*(p2
     -   mup4*(p4mup5*F(226)+(p1mup5+p2mup5)*F(228)-p6mup5*F(229))+p4
     -   mup4*(-(p6mup5*F(227))+p2mup5*F(231))+p6mup4*((p1mup5+p2mup5
     -   )*F(232)-p4mup5*F(233)-p6mup5*F(234))))
       Fa10(107) = p2mup3*p4mup2*(p1mup5*p2mup4*F(228)+p2mup4*p2mup5*F
     -   (228)-p1mup4*p6mup5*F(229)-p2mup4*p6mup5*F(229)+p1mup5*p6mup
     -   4*F(232)+p2mup5*p6mup4*F(232)-p6mup4*p6mup5*F(234))
       Fa10(108) = (p1mup2+p2mup2)*p4mup3*(mup4mup5*F(225)+2*(p3mup4*p
     -   4mup5*F(226)-p4mup4*p6mup5*F(227)+p1mup5*p3mup4*F(228)+p2mup
     -   5*p3mup4*F(228)+p3mup4*p3mup5*F(228)-p3mup4*p6mup5*F(229)+p1
     -   mup4*(p4mup5*F(226)+(p1mup5+p2mup5+p3mup5)*F(228)-p6mup5*F(2
     -   29))+p2mup4*(p4mup5*F(226)+(p1mup5+p2mup5+p3mup5)*F(228)-p6m
     -   up5*F(229))+p4mup4*p4mup5*F(230)+p1mup5*p4mup4*F(231)+p2mup5
     -   *p4mup4*F(231)+p3mup5*p4mup4*F(231)+p1mup5*p6mup4*F(232)+p2m
     -   up5*p6mup4*F(232)+p3mup5*p6mup4*F(232)-p4mup5*p6mup4*F(233)-
     -   p6mup4*p6mup5*F(234)))
       Fa10(109) = mup4mup5*(p1mup2+p2mup2)*(p1mup3+p2mup3+p3mup3)
       Fa10(110) = p2mup2*(p1mup4*p3mup3*p4mup5*F(226)+p1mup3*p3mup4*p
     -   4mup5*F(226)+p1mup3*p1mup4*p1mup5*F(228)+p3mup3*p4mup4*p4mup
     -   5*F(230)+p1mup3*p3mup5*p4mup4*F(231))+p1mup2*(p1mup4*p3mup3*
     -   p4mup5*F(226)+p2mup4*p3mup3*p4mup5*F(226)+p1mup3*p3mup4*p4mu
     -   p5*F(226)+p2mup3*p3mup4*p4mup5*F(226)+p3mup3*p3mup4*p4mup5*F
     -   (226)-p3mup3*p4mup4*p6mup5*F(227)+p1mup3*p1mup5*p2mup4*F(228
     -   )+p1mup4*(p1mup5*p2mup3+p1mup3*p2mup5)*F(228)+p3mup3*p4mup4*
     -   p4mup5*F(230)+p1mup5*p3mup3*p4mup4*F(231)+p2mup5*p3mup3*p4mu
     -   p4*F(231)+p1mup3*p3mup5*p4mup4*F(231)+p2mup3*p3mup5*p4mup4*F
     -   (231)+p3mup3*p3mup5*p4mup4*F(231)-p3mup3*p4mup5*p6mup4*F(233
     -   ))
       Return
       End
