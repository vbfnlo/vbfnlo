c************** Calling the Fa functions*************************
       subroutine HexBoxCross2Fa7(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p
     &   3mup6,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa7
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
       Complex*16 Fa7(115:133)
       COMMON/HexBoxCross2FaFunctions/Fa
       Fa7(115) = p2mup3*F(300)+p5mup3*F(472)
       Fa7(116) = p2mup3*F(471)+p1mup3*F(472)
       Fa7(117) = p1mup3*p2mup6*F(475)-p1mup3*p6mup6*F(476)+p1mup6*(p6
     &   mup3*F(321)+p1mup3*F(474)+p5mup3*F(477))-p1mup3*p5mup6*F(478
     &   )+p5mup3*p6mup6*F(479)+p6mup3*(p6mup6*F(323)+2*p2mup6*F(327)
     &   +p5mup6*F(480))+p5mup3*p5mup6*F(481)-p2mup6*p5mup3*F(482)
       Fa7(118) = p5mup3*F(413)+p6mup3*F(484)
       Fa7(119) = p1mup3*F(393)+p6mup3*F(485)
       Fa7(120) = p5mup3*F(402)+p1mup3*F(453)-6*p6mup3*F(1845)
       Fa7(121) = p5mup6*F(489)+p6mup6*F(490)+p2mup6*F(491)+p1mup6*F(4
     &   92)
       Fa7(122) = p5mup6*F(493)+p6mup6*F(494)-p2mup6*F(495)+p1mup6*F(4
     &   96)
       Fa7(123) = 2*p5mup3*F(498)+p6mup3*F(499)
       Fa7(124) = p1mup3*F(500)-p6mup3*F(501)
       Fa7(125) = p5mup3*F(501)+p1mup3*F(502)
       Fa7(126) = -(p6mup6*F(301))-2*(p5mup6*F(304)-p2mup6*F(305)-p1mu
     &   p6*F(306))
       Fa7(127) = p5mup6*F(503)+p6mup6*F(504)-2*(p2mup6*F(505)+p1mup6*
     &   F(506))
       Fa7(128) = p5mup6*F(340)+p6mup6*F(341)+p1mup6*F(342)+p2mup6*F(3
     &   43)
       Fa7(129) = p2mup6*F(277)+p6mup6*F(283)+2*(p5mup6*F(286)+p1mup6*
     &   F(287))
       Fa7(130) = p2mup6*F(278)+p6mup6*F(284)+p5mup6*F(285)
       Fa7(131) = p1mup6*F(315)+p2mup6*F(316)-p6mup6*F(317)-p5mup6*F(3
     &   18)
       Fa7(132) = p1mup6*F(321)+p6mup6*F(323)+2*p2mup6*F(327)+p5mup6*F
     &   (480)
       Fa7(133) = p6mup3*F(300)+p5mup3*F(450)
       Return
       End
