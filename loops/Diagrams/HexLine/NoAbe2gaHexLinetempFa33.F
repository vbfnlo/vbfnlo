c************** Calling the Fa functions*************************
       subroutine NoAbe2gaHexLineFa33(p1mup2,p1mup3,p1mup5,p1mup6,p2mu
     -   p2,p2mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p
     -   4mup3,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup
     -   3,p6mup5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6
     -   ,mup5mup6,Fa33)
       IMPLICIT NONE
      Complex*16   p1mup2, p1mup3, p1mup5, p1mup6, p2mup2, p2mup3, 
     -          p2mup5, p2mup6, p3mup2, p3mup3, p3mup5, p3mup6, 
     -          p4mup2, p4mup3, p4mup5, p4mup6, p5mup2, p5mup3, 
     -          p5mup5, p5mup6, p6mup2, p6mup3, p6mup5, p6mup6
       Complex*16   mup2mup3, mup2mup5, mup2mup6, mup3mup5, mup3mup6, 
     -          mup5mup6
        common/NoAbe2gaHexLineFhlFunctions/F
       COMMON/NoAbe2gaHexLineInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,
     -   s12,s23,s34,s45,s56,s16,s123,s234,s345
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
      Complex*16 Fa(872),F(9706)
      Real*8 P(9687) 
        COMMON/NoAbe2gaHexLinePFunctions/P
       Complex*16 Fa33(769:792)
       COMMON/NoAbe2gaHexLineFaFunctions/Fa
       Fa33(769) = p1mup2*p3mup5*p3mup6*F(5624)+p2mup2*p3mup5*p3mup6*F
     -   (9118)-p3mup5*p4mup6*p5mup2*F(9119)+p3mup6*p4mup5*p5mup2*F(9
     -   120)+p4mup5*p4mup6*p5mup2*F(9122)+p3mup2*(-(p3mup5*p4mup6*F(
     -   9116))+p3mup6*p4mup5*F(9117)+p4mup5*p4mup6*F(9121)+p3mup5*p3
     -   mup6*F(9683))+p3mup5*p3mup6*p4mup2*F(9685)
       Fa33(770) = -(mup3mup6*(p4mup2*F(9124)+p5mup2*F(9125)))-mup3mup
     -   6*p3mup2*F(9126)+p4mup2*p4mup3*p6mup6*F(9135)+p4mup2*p5mup3*
     -   p6mup6*F(9138)+p3mup2*p4mup3*p6mup6*F(9139)+p3mup2*p5mup3*p6
     -   mup6*F(9141)+p5mup2*p5mup3*p6mup6*F(9143)+p4mup3*p5mup2*p6mu
     -   p6*F(9146)+p2mup2*(-(mup3mup6*F(9123))+p4mup3*p6mup6*F(9129)
     -   +p5mup3*p6mup6*F(9130)+p2mup3*p6mup6*F(9133)+4*(p5mup3*p5mup
     -   6*F(9127)+p4mup3*p5mup6*F(9128)+p2mup6*(2*(p3mup3*F(6898)+p2
     -   mup3*F(7046)-p4mup3*F(7143))+p5mup3*F(9131))-p2mup3*p5mup6*F
     -   (9132)+p3mup3*p5mup6*F(9145))+p3mup3*p6mup6*F(9147))+p2mup3*
     -   p3mup2*p6mup6*F(9157)+p2mup3*p4mup2*p6mup6*F(9160)+p2mup3*p5
     -   mup2*p6mup6*F(9164)+p3mup3*p4mup2*p6mup6*F(9167)+p3mup2*p3mu
     -   p3*p6mup6*F(9169)+p3mup3*p5mup2*p6mup6*F(9171)+mup2mup3*p6mu
     -   p6*F(9175)+4*(p4mup2*p4mup3*p5mup6*F(9134)+p4mup2*p5mup3*p5m
     -   up6*F(9136)-p3mup2*p4mup3*p5mup6*F(9137)-p3mup2*p5mup3*p5mup
     -   6*F(9140)+p5mup2*p5mup3*p5mup6*F(9142)+p4mup3*p5mup2*p5mup6*
     -   F(9144)+p2mup6*p3mup2*p5mup3*F(9155)+p2mup3*p3mup2*p5mup6*F(
     -   9156)+p2mup3*p4mup2*p5mup6*F(9159)+p2mup6*p4mup2*p5mup3*F(91
     -   61)+p2mup3*p5mup2*p5mup6*F(9162)+p2mup6*p5mup2*p5mup3*F(9163
     -   )+p3mup3*p4mup2*p5mup6*F(9166)+p3mup2*p3mup3*p5mup6*F(9168)+
     -   p3mup3*p5mup2*p5mup6*F(9170)+mup2mup6*p5mup3*F(9173)+mup2mup
     -   3*p5mup6*F(9174)+p1mup6*(p2mup2*(p5mup3*F(9148)+2*(p3mup3*F(
     -   6801)-p4mup3*F(7155)+p2mup3*F(9149)))+p4mup2*p5mup3*F(9152)+
     -   p5mup2*p5mup3*F(9153)+p3mup2*(p5mup3*F(9150)+2*(p4mup3*F(716
     -   2)+p2mup3*F(9663)+p3mup3*F(9667)))+2*(p4mup2*p4mup3*F(9151)+
     -   p4mup3*p5mup2*F(9154)+p2mup3*p4mup2*F(9664)+p2mup3*p5mup2*F(
     -   9665)+p3mup3*p4mup2*F(9668)+p3mup3*p5mup2*F(9669)+mup2mup3*F
     -   (9694))))+8*(p2mup6*p3mup2*p4mup3*F(7207)+p2mup6*p4mup2*p4mu
     -   p3*F(9158)+p2mup6*p4mup3*p5mup2*F(9165)+mup2mup6*p4mup3*F(91
     -   72)+p2mup3*p2mup6*p3mup2*F(9673)+p2mup3*p2mup6*p4mup2*F(9674
     -   )+p2mup3*p2mup6*p5mup2*F(9675)+p2mup6*p3mup2*p3mup3*F(9676)+
     -   p2mup6*p3mup3*p4mup2*F(9678)+p2mup6*p3mup3*p5mup2*F(9679)+mu
     -   p2mup6*p3mup3*F(9693)+(mup2mup6*p2mup3+mup2mup3*p2mup6)*F(96
     -   95))
       Fa33(771) = p1mup2*(p4mup3*p6mup6*F(9177)+p5mup3*p6mup6*F(9179)
     -   +mup3mup6*F(9180)+p2mup3*p6mup6*F(9182)+p3mup3*p6mup6*F(9185
     -   )+p1mup3*p6mup6*F(9187)+4*(p4mup3*p5mup6*F(9176)+p5mup3*p5mu
     -   p6*F(9178)-p2mup3*p5mup6*F(9181)+p3mup3*p5mup6*F(9184)-p1mup
     -   3*p5mup6*F(9186)+p2mup6*(p5mup3*F(9183)+2*(p3mup3*F(5500)+p1
     -   mup3*F(5598)-p4mup3*F(7140)+p2mup3*F(9197)))+p1mup6*(p5mup3*
     -   F(9190)+2*(p3mup3*F(5490)+p2mup3*F(5598)-p4mup3*F(7153)+p1mu
     -   p3*F(9198)))))+p1mup3*(p4mup2*p6mup6*F(9194)+4*(p4mup2*p5mup
     -   6*F(9193)+p5mup2*p5mup6*F(9195))+p5mup2*p6mup6*F(9196)+p2mup
     -   2*(-4*p5mup6*F(9188)+p6mup6*F(9189)+8*(p2mup6*F(9149)+p1mup6
     -   *F(9199)))+p3mup2*(4*p5mup6*F(9191)+p6mup6*F(9192)+8*(p1mup6
     -   *F(9659)+p2mup6*F(9663)))+8*(p1mup6*p4mup2*F(9660)+p1mup6*p5
     -   mup2*F(9661)+p2mup6*p4mup2*F(9664)+p2mup6*p5mup2*F(9665)+mup
     -   2mup6*F(9694)))
       Fa33(772) = p1mup2*(-8*p4mup6*p6mup3*F(4052)+p6mup3*p6mup6*F(92
     -   02)+4*(p5mup6*p6mup3*F(9200)+p4mup6*p5mup3*F(9204)))+p2mup2*
     -   (-8*p4mup6*p6mup3*F(7106)+p6mup3*p6mup6*F(9203)+4*(p5mup6*p6
     -   mup3*F(9201)+p4mup6*p5mup3*F(9205)))
       Fa33(773) = p1mup2*(2*p4mup6*(p1mup3*F(5601)+p2mup3*F(5612)-p4m
     -   up3*F(9207))+p3mup6*(-2*(p6mup3*F(3977)-p1mup3*F(5490)-p2mup
     -   3*F(5500)-p4mup3*F(7007))+p5mup3*F(9208)))+p2mup2*p3mup6*p5m
     -   up3*F(9209)+2*(p1mup3*p2mup2*(p3mup6*F(6801)+p4mup6*F(6893))
     -   +p2mup2*p2mup3*p3mup6*F(6898)+p2mup2*p2mup3*p4mup6*F(6901)-p
     -   2mup2*p3mup6*p6mup3*F(6937)+p2mup2*p3mup6*p4mup3*F(7016)+p3m
     -   up2*p3mup6*p4mup3*F(7078)+p3mup6*p4mup3*p5mup2*F(7087)-p2mup
     -   2*p4mup3*p4mup6*F(9206)+p3mup6*p4mup2*p4mup3*F(9210)+p3mup2*
     -   p4mup3*p4mup6*F(9210)+p4mup2*p4mup3*p4mup6*F(9211)+p4mup3*p4
     -   mup6*p5mup2*F(9212))
       Fa33(774) = p3mup6*p4mup2*p5mup3*F(9213)+p3mup2*p4mup6*p5mup3*F
     -   (9213)+p3mup2*p3mup6*p5mup3*F(9214)+p3mup6*p5mup2*p5mup3*F(9
     -   215)+p4mup2*p4mup6*p5mup3*F(9216)+p4mup6*p5mup2*p5mup3*F(921
     -   7)+2*(-(p1mup2*(p2mup6*p6mup3*F(3975)+p1mup6*p6mup3*F(3978)-
     -   p3mup3*(p4mup6*F(5618)+p3mup6*F(5695))))+p2mup2*(p3mup3*p4mu
     -   p6*F(6905)-p2mup6*p6mup3*F(6924)-p1mup6*p6mup3*F(6944)+p3mup
     -   3*p3mup6*F(6951))+p3mup2*p3mup3*p3mup6*F(9666)+p1mup3*p3mup2
     -   *p3mup6*F(9667)+p1mup3*p3mup6*p4mup2*F(9668)+p1mup3*p3mup2*p
     -   4mup6*F(9668)+p1mup3*p3mup6*p5mup2*F(9669)+p1mup3*p4mup2*p4m
     -   up6*F(9670)+p1mup3*p4mup6*p5mup2*F(9671)+p2mup3*p3mup2*p3mup
     -   6*F(9676)+p2mup3*p3mup6*p4mup2*F(9678)+p2mup3*p3mup2*p4mup6*
     -   F(9678)+p2mup3*p3mup6*p5mup2*F(9679)+p2mup3*p4mup2*p4mup6*F(
     -   9680)+p2mup3*p4mup6*p5mup2*F(9681)+p3mup3*p3mup6*p4mup2*F(96
     -   82)+p3mup2*p3mup3*p4mup6*F(9682)+p3mup3*p3mup6*p5mup2*F(9683
     -   )+p3mup3*p4mup2*p4mup6*F(9684)+p3mup3*p4mup6*p5mup2*F(9685)+
     -   mup2mup3*p3mup6*F(9693)+mup2mup3*p4mup6*F(9696))
       Fa33(775) = p6mup3*(p3mup2*(p5mup6*F(9220)+2*(p1mup6*F(9630)+p2
     -   mup6*F(9637)+p4mup6*F(9642)+p3mup6*F(9643)))+p4mup2*(-(p5mup
     -   6*F(9218))+2*(p1mup6*F(9631)+p2mup6*F(9638)+p3mup6*F(9642)+p
     -   4mup6*F(9645)))+p5mup2*(-(p5mup6*F(9219))+2*(p1mup6*F(9632)+
     -   p2mup6*F(9639)+p3mup6*F(9644)+p4mup6*F(9646)))+2*mup2mup6*F(
     -   9652))
       Fa33(776) = p6mup3*p6mup6*(p4mup2*F(9221)-4*p3mup2*F(9222)+p5mu
     -   p2*F(9223))
       Fa33(777) = p4mup5*(p2mup2*(p4mup3*F(9224)+p5mup3*F(9225)-p2mup
     -   3*F(9226)-p3mup3*F(9228))+p4mup2*p4mup3*F(9229)+p4mup2*p5mup
     -   3*F(9230)+p4mup3*p5mup2*F(9231)+p5mup2*p5mup3*F(9232)+p2mup3
     -   *p4mup2*F(9234)+p2mup3*p5mup2*F(9235)+p3mup2*(p5mup3*F(7397)
     -   +p4mup3*F(9227)+p2mup3*F(9233)+p3mup3*F(9236))+p3mup3*p4mup2
     -   *F(9237)+p3mup3*p5mup2*F(9238)+mup2mup3*F(9239))
       Fa33(778) = p4mup5*(p1mup2*(p4mup3*F(9240)+p5mup3*F(9241)-p1mup
     -   3*F(9242)-p2mup3*F(9244)-p3mup3*F(9245))+p1mup3*(-(p2mup2*F(
     -   9243))+p3mup2*F(9246)+p4mup2*F(9247)+p5mup2*F(9248)))
       Fa33(779) = mup3mup5*(p2mup2*F(9249)+p1mup2*F(9250)+p3mup2*F(92
     -   51)+p4mup2*F(9252)+p5mup2*F(9253))
       Fa33(780) = p4mup5*p6mup3*(p1mup2*F(9254)+p2mup2*F(9255))
       Fa33(781) = p4mup3*(p3mup2*p5mup5*F(9259)-p4mup2*p5mup5*F(9260)
     -   -p5mup2*p5mup5*F(9262)+p2mup2*(-(p5mup5*F(9256))+8*(p2mup5*F
     -   (9257)+p1mup5*F(9264)))+8*(p1mup2*p2mup5*F(9258)+p1mup2*p1mu
     -   p5*F(9263)+p2mup5*p3mup2*F(9265)+p2mup5*p4mup2*F(9266)+p2mup
     -   5*p5mup2*F(9267)+p3mup5*(p4mup2*F(9227)-p2mup2*F(9261)-p3mup
     -   2*F(9268)+p5mup2*F(9269))+mup2mup5*F(9270)))
       Fa33(782) = p1mup2*p4mup3*p5mup5*F(9271)+p1mup2*p5mup3*p5mup5*F
     -   (9272)+p1mup2*p5mup5*p6mup3*F(9274)-p1mup2*p2mup3*p5mup5*F(9
     -   276)+p1mup2*p3mup3*p5mup5*F(9285)+p5mup2*p5mup3*p5mup5*F(928
     -   7)-p1mup2*p1mup3*p5mup5*F(9288)-p2mup3*p5mup2*p5mup5*F(9298)
     -   +p2mup2*(p5mup3*p5mup5*F(9273)+p5mup5*p6mup3*F(9275)+4*p3mup
     -   5*p6mup3*F(9284)+p3mup3*p5mup5*F(9286)-p1mup3*p5mup5*F(9289)
     -   +p2mup3*(-(p5mup5*F(9277))+8*(p2mup5*F(5962)+p3mup5*F(9297))
     -   )-8*(p2mup5*p5mup3*F(9279)-p3mup5*p5mup3*F(9283)-p1mup3*p3mu
     -   p5*F(9292)-p2mup5*p3mup3*F(9297)-p3mup3*p3mup5*F(9300)))-p3m
     -   up3*p5mup2*p5mup5*F(9301)-8*(p2mup5*p5mup2*p5mup3*F(7495)+p1
     -   mup2*p2mup5*p5mup3*F(9278)+p1mup5*p4mup3*(p3mup2*F(9290)+p4m
     -   up2*F(9293))+p1mup5*p4mup3*p5mup2*F(9294)-p1mup2*p2mup3*p2mu
     -   p5*F(9295)+p2mup5*p3mup3*p5mup2*F(9302))+p3mup5*(4*p1mup2*(p
     -   6mup3*F(9282)+2*(p4mup3*F(9280)+p5mup3*F(9281)+p1mup3*F(9291
     -   )+p2mup3*F(9296)+p3mup3*F(9299)))-8*p5mup2*(p5mup3*F(6028)+p
     -   2mup3*F(9302)+p3mup3*F(9303)))
       Fa33(783) = -(p4mup2*p5mup3*p5mup5*F(9306))+p3mup2*p5mup3*p5mup
     -   5*F(9307)+p1mup3*p3mup2*p5mup5*F(9314)+p1mup3*p4mup2*p5mup5*
     -   F(9315)+p1mup3*p5mup2*p5mup5*F(9317)+p2mup3*p3mup2*p5mup5*F(
     -   9319)+p2mup3*p4mup2*p5mup5*F(9320)+p3mup3*p4mup2*p5mup5*F(93
     -   21)+p3mup2*p3mup3*p5mup5*F(9322)+mup2mup3*p5mup5*F(9324)+8*(
     -   p1mup5*p3mup2*p5mup3*F(7474)+p1mup5*p3mup3*p4mup2*F(9246)-p1
     -   mup2*p1mup5*p3mup3*F(9291)-p1mup5*p2mup2*p3mup3*F(9292)+p1mu
     -   p2*p1mup5*p5mup3*F(9308)+p1mup5*p2mup2*p5mup3*F(9309)-p1mup2
     -   *p1mup5*p2mup3*F(9312)-p1mup5*p2mup2*p2mup3*F(9313)+p1mup5*p
     -   4mup2*p5mup3*F(9316)+p1mup5*p5mup2*p5mup3*F(9318)+mup2mup5*p
     -   5mup3*F(9323)-p1mup3*p1mup5*p2mup2*F(9325)-p1mup2*p1mup3*p1m
     -   up5*F(9326)+p1mup3*p1mup5*p3mup2*F(9327)+p1mup3*p1mup5*p4mup
     -   2*F(9328)+p1mup3*p1mup5*p5mup2*F(9329)+p1mup5*p2mup3*p3mup2*
     -   F(9330)+p1mup5*p2mup3*p4mup2*F(9331)+p1mup5*p2mup3*p5mup2*F(
     -   9332)+p1mup5*p3mup2*p3mup3*F(9334)+p1mup5*p3mup3*p5mup2*F(93
     -   35)+mup2mup5*p3mup3*F(9340)+p3mup5*(p4mup2*(p5mup3*F(7397)+p
     -   2mup3*F(9233)+p3mup3*F(9236)+p1mup3*F(9246))+p1mup3*p5mup2*F
     -   (9335)+p3mup2*(p5mup3*F(7497)+p3mup3*F(9333)+p1mup3*F(9334)+
     -   p2mup3*F(9339))+mup2mup3*F(9340))+mup2mup5*p1mup3*F(9341)+mu
     -   p2mup3*p1mup5*F(9341)+mup2mup5*p2mup3*F(9342))+4*(p1mup2*p1m
     -   up5*p6mup3*F(9310)+p1mup5*p2mup2*p6mup3*F(9311)+p2mup5*(p1mu
     -   p2*p6mup3*F(9304)+p2mup2*p6mup3*F(9305)-2*(p1mup2*p3mup3*F(9
     -   296)+p1mup2*p1mup3*F(9312)+p1mup3*p2mup2*F(9313)-p1mup3*p5mu
     -   p2*F(9332)-p4mup2*(p5mup3*F(7395)+p3mup3*F(9233)+p1mup3*F(93
     -   31)+p2mup3*F(9337))-p2mup3*p5mup2*F(9338)-p3mup2*(p5mup3*F(7
     -   494)+p1mup3*F(9330)+p2mup3*F(9336)+p3mup3*F(9339))-mup2mup3*
     -   F(9342))))
       Fa33(784) = p6mup5*(-(p1mup2*p4mup3*F(9343))-p1mup2*p5mup3*F(93
     -   44)-p2mup2*p4mup3*F(9345)-p2mup2*p5mup3*F(9346)+p2mup3*(p3mu
     -   p2*F(7511)+p4mup2*F(7512)+p5mup2*F(7532)+p2mup2*F(9347)+p1mu
     -   p2*F(9348))-p4mup2*p4mup3*F(9349)+p3mup2*p4mup3*F(9350)-p4mu
     -   p2*p5mup3*F(9351)-p5mup2*p5mup3*F(9352)+p3mup3*(p4mup2*F(751
     -   3)+p3mup2*F(7534)+p5mup2*F(7536)+p1mup2*F(9353)+p2mup2*F(935
     -   4))+p3mup2*p5mup3*F(9355)-p4mup3*p5mup2*F(9356)+p1mup2*p1mup
     -   3*F(9357)+p1mup3*p2mup2*F(9358)+p1mup3*p3mup2*F(9359)+p1mup3
     -   *p4mup2*F(9360)+p1mup3*p5mup2*F(9361)+mup2mup3*F(9362))
       Fa33(785) = p6mup3*(-16*p6mup5*(p1mup2*F(9363)+p2mup2*F(9364))+
     -   4*p5mup5*(p4mup2*F(9365)+p3mup2*F(9366)+p5mup2*F(9367))+p2mu
     -   p5*p5mup2*F(9368))
       Fa33(786) = p6mup3*(p1mup5*(p3mup2*F(9370)+p4mup2*F(9371)+p5mup
     -   2*F(9372))+p2mup5*p3mup2*F(9373)+p2mup5*p4mup2*F(9374)+p3mup
     -   5*p4mup2*F(9375)+p3mup2*p4mup5*F(9375)+p3mup2*p3mup5*F(9376)
     -   +p3mup5*p5mup2*F(9377)+p4mup2*p4mup5*F(9378)+p4mup5*p5mup2*F
     -   (9379))
       Fa33(787) = p6mup3*p6mup5*(p4mup2*F(4314)+p3mup2*F(4315)+p5mup2
     -   *F(8318))
       Fa33(788) = p2mup2*F(9380)+p4mup2*F(9381)+p5mup2*F(9382)+p3mup2
     -   *F(9383)-p1mup2*F(9384)
       Fa33(789) = mup5mup6*F(9385)-4*p3mup6*(3*p5mup5*F(5684)+2*(p1mu
     -   p5*F(9391)-p4mup5*F(9405))+6*(p6mup5*F(9621)-p2mup5*F(9639)-
     -   p3mup5*F(9644)))-2*(-(p4mup5*p6mup6*F(9393))-4*(p4mup5*p5mup
     -   6*F(9392)+p4mup6*p6mup5*F(9397))+p5mup5*(p1mup6*F(9387)-2*p4
     -   mup6*F(9395)-3*p6mup6*F(9400)+6*(p5mup6*F(5688)+p2mup6*F(940
     -   3)))-2*(p1mup5*(p6mup6*F(9386)-2*(p2mup6*F(9388)+p4mup6*F(93
     -   96)+p5mup6*F(9398)+p1mup6*F(9407)))+3*(p6mup5*p6mup6*F(9617)
     -   -p2mup5*p6mup6*F(9619)-p3mup5*p6mup6*F(9621))+2*(p1mup6*(p2m
     -   up5*F(9389)+p3mup5*F(9390)+p4mup5*F(9394)+p6mup5*F(9399))+p2
     -   mup5*p4mup6*F(9402)+p3mup5*p4mup6*F(9404)+p4mup5*p4mup6*F(94
     -   06)+p2mup6*(p4mup5*F(9401)-3*(p6mup5*F(9619)-p2mup5*F(9636)-
     -   p3mup5*F(9639))))-6*(p5mup6*p6mup5*F(9625)-p2mup5*p5mup6*F(9
     -   641)-p3mup5*p5mup6*F(9647))))
       Fa33(790) = p1mup5*F(9408)+p5mup5*F(9409)+p6mup5*F(9410)+p2mup5
     -   *F(9411)+p3mup5*F(9412)+p4mup5*F(9413)
       Fa33(791) = p1mup5*F(9414)+p5mup5*F(9415)+p6mup5*F(9416)+p2mup5
     -   *F(9417)+p3mup5*F(9418)+p4mup5*F(9419)
       Fa33(792) = p1mup6*F(9420)+p5mup6*F(9421)+p6mup6*F(9422)+p2mup6
     -   *F(9423)+p3mup6*F(9424)+p4mup6*F(9425)
       Return
       End
