c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineAbe12(p1mup2,p1mup3,p1mup4,p1m
     -   up5,p2mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,
     -   p4mup2,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mu
     -   p2,p6mup3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,
     -   mup3mup5,mup4mup5,Fa12)
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
       Complex*16 Fa12(122:132)
       COMMON/test_FaFunctionsAbe/Fa
       Fa12(122) = p3mup2*(mup4mup5*F(225)+2*(p3mup4*p4mup5*F(226)-p4m
     -   up4*p6mup5*F(227)+p1mup5*p3mup4*F(228)+p2mup5*p3mup4*F(228)+
     -   p3mup4*p3mup5*F(228)-p3mup4*p6mup5*F(229)+p1mup4*(p4mup5*F(2
     -   26)+(p1mup5+p2mup5+p3mup5)*F(228)-p6mup5*F(229))+p2mup4*(p4m
     -   up5*F(226)+(p1mup5+p2mup5+p3mup5)*F(228)-p6mup5*F(229))+p4mu
     -   p4*p4mup5*F(230)+p1mup5*p4mup4*F(231)+p2mup5*p4mup4*F(231)+p
     -   3mup5*p4mup4*F(231)+p1mup5*p6mup4*F(232)+p2mup5*p6mup4*F(232
     -   )+p3mup5*p6mup4*F(232)-p4mup5*p6mup4*F(233)-p6mup4*p6mup5*F(
     -   234)))
       Fa12(123) = (p1mup3+p2mup3+p3mup3)*(mup4mup5*F(225)+2*(p3mup4*p
     -   4mup5*F(226)-p4mup4*p6mup5*F(227)+p1mup5*p3mup4*F(228)+p2mup
     -   5*p3mup4*F(228)+p3mup4*p3mup5*F(228)-p3mup4*p6mup5*F(229)+p1
     -   mup4*(p4mup5*F(226)+(p1mup5+p2mup5+p3mup5)*F(228)-p6mup5*F(2
     -   29))+p2mup4*(p4mup5*F(226)+(p1mup5+p2mup5+p3mup5)*F(228)-p6m
     -   up5*F(229))+p4mup4*p4mup5*F(230)+p1mup5*p4mup4*F(231)+p2mup5
     -   *p4mup4*F(231)+p3mup5*p4mup4*F(231)+p1mup5*p6mup4*F(232)+p2m
     -   up5*p6mup4*F(232)+p3mup5*p6mup4*F(232)-p4mup5*p6mup4*F(233)-
     -   p6mup4*p6mup5*F(234)))
       Fa12(124) = (p1mup2+p2mup2)*(mup4mup5*F(225)+2*(p3mup4*p4mup5*F
     -   (226)-p4mup4*p6mup5*F(227)+p1mup5*p3mup4*F(228)+p2mup5*p3mup
     -   4*F(228)+p3mup4*p3mup5*F(228)-p3mup4*p6mup5*F(229)+p1mup4*(p
     -   4mup5*F(226)+(p1mup5+p2mup5+p3mup5)*F(228)-p6mup5*F(229))+p2
     -   mup4*(p4mup5*F(226)+(p1mup5+p2mup5+p3mup5)*F(228)-p6mup5*F(2
     -   29))+p4mup4*p4mup5*F(230)+p1mup5*p4mup4*F(231)+p2mup5*p4mup4
     -   *F(231)+p3mup5*p4mup4*F(231)+p1mup5*p6mup4*F(232)+p2mup5*p6m
     -   up4*F(232)+p3mup5*p6mup4*F(232)-p4mup5*p6mup4*F(233)-p6mup4*
     -   p6mup5*F(234)))
       Fa12(125) = (p1mup2+p2mup2)*(p6mup5*F(219)+(p1mup5+p2mup5+p3mup
     -   5)*F(220)+p4mup5*F(222))
       Fa12(126) = (p1mup2+p2mup2)*(p6mup4*F(202)+(p1mup4+p2mup4+p3mup
     -   4)*F(203)-p4mup4*F(204))
       Fa12(127) = p6mup4*F(202)+(p1mup4+p2mup4+p3mup4)*F(203)-p4mup4*
     -   F(204)
       Fa12(128) = p6mup5*F(219)+(p1mup5+p2mup5+p3mup5)*F(220)+p4mup5*
     -   F(222)
       Fa12(129) = p3mup2*F(246)+p4mup2*F(247)+p1mup2*F(248)+p2mup2*F(
     -   249)
       Fa12(130) = p1mup3*F(250)+p2mup3*F(251)+p3mup3*F(252)+4*p4mup3*
     -   F(253)
       Fa12(131) = mup2mup3*F(254)+p3mup2*p4mup3*F(257)+p1mup3*p4mup2*
     -   F(258)+p2mup3*p4mup2*F(259)+p3mup3*p4mup2*F(260)+p2mup2*(p4m
     -   up3*F(256)+p2mup3*F(261)+p1mup3*F(262)+p3mup3*F(264))+p2mup3
     -   *p3mup2*F(265)+p1mup3*p3mup2*F(266)+p4mup2*p4mup3*F(267)+p1m
     -   up2*(p4mup3*F(255)+p2mup3*F(263)+p1mup3*F(268)+p3mup3*F(269)
     -   )+p3mup2*p3mup3*F(270)
       Fa12(132) = mup2mup3*F(271)+p3mup3*p4mup2*F(273)+p1mup3*p2mup2*
     -   F(274)+p2mup2*p3mup3*F(275)+p1mup3*p4mup2*F(277)+p2mup3*p4mu
     -   p2*F(278)+p2mup2*p4mup3*F(282)+p3mup2*(p4mup3*F(272)+p1mup3*
     -   F(280)+p2mup3*F(281)+p3mup3*F(283))+p4mup2*p4mup3*F(284)+p1m
     -   up2*(p4mup3*F(276)+p3mup3*F(279)+p1mup3*F(285)+p2mup3*F(286)
     -   )+p2mup2*p2mup3*F(287)
       Return
       End
