c************** Calling the Fa functions*************************
       subroutine NoAbe2gaHexLineFa28(p1mup2,p1mup3,p1mup5,p1mup6,p2mu
     -   p2,p2mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p
     -   4mup3,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup
     -   3,p6mup5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6
     -   ,mup5mup6,Fa28)
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
       Complex*16 Fa28(649:672)
       COMMON/NoAbe2gaHexLineFaFunctions/Fa
       Fa28(649) = mup5mup6*F(7297)+p1mup6*p5mup5*F(7301)+p1mup6*p6mup
     -   5*F(7302)+p5mup5*p5mup6*F(7303)+p5mup5*p6mup6*F(7304)+p5mup6
     -   *p6mup5*F(7305)+p6mup5*p6mup6*F(7306)+2*p4mup6*p5mup5*F(7308
     -   )+p4mup6*p6mup5*F(7310)+p1mup6*p2mup5*F(7312)+p2mup5*p6mup6*
     -   F(7313)+p2mup6*p6mup5*F(7315)+p1mup6*p3mup5*F(7319)+p3mup5*p
     -   6mup6*F(7321)+p1mup5*(p1mup6*F(7299)+p6mup6*F(7300)+p5mup6*F
     -   (7318)-4*(p4mup6*F(7309)-p2mup6*F(7311)-p3mup6*F(7322)))+4*(
     -   p2mup6*p5mup5*F(7314)+p2mup5*p4mup6*F(7317)+p3mup5*p4mup6*F(
     -   7320)+p3mup6*p5mup5*F(7324))+p3mup6*p6mup5*F(7325)+p3mup5*p5
     -   mup6*F(7328)+p2mup5*p5mup6*F(7330)-8*(p2mup5*p2mup6*F(7327)+
     -   p3mup5*p3mup6*F(7329)+(p2mup6*p3mup5+p2mup5*p3mup6)*F(7331))
     -   +p4mup5*(p5mup6*F(7298)+p1mup6*F(7307)+p2mup6*F(7316)+p3mup6
     -   *F(7323)+p4mup6*F(7326)+p6mup6*F(7332))
       Fa28(650) = mup3mup6*F(7333)+p5mup6*p6mup3*F(7339)+p1mup6*p2mup
     -   3*F(7343)+p2mup3*p5mup6*F(7344)+p1mup6*p3mup3*F(7346)+p1mup6
     -   *p4mup3*F(7347)+p4mup3*p5mup6*F(7348)+p4mup3*p6mup6*F(7349)+
     -   p6mup3*p6mup6*F(7351)+p5mup3*(p5mup6*F(7334)+p1mup6*F(7335)+
     -   p6mup6*F(7336)+p2mup6*F(7341)+p3mup6*F(7353)+p4mup6*F(7355))
     -   -2*(p1mup6*p6mup3*F(7340)-p2mup6*p6mup3*F(7345)-p3mup6*p6mup
     -   3*F(7354)-p4mup6*p6mup3*F(7356))+p1mup3*(p1mup6*F(7337)+p5mu
     -   p6*F(7338)+p2mup6*F(7342)+p6mup6*F(7350)+p3mup6*F(7352)+p4mu
     -   p6*F(7357))+p2mup3*p2mup6*F(7358)+p3mup3*p6mup6*F(7359)+p2mu
     -   p3*p3mup6*F(7360)+p3mup6*p4mup3*F(7361)+p4mup3*p4mup6*F(7362
     -   )+p3mup3*p5mup6*F(7363)+p2mup6*p4mup3*F(7364)+p2mup3*p6mup6*
     -   F(7365)+p2mup3*p4mup6*F(7366)+p3mup3*p4mup6*F(7367)+p3mup3*p
     -   3mup6*F(7368)+p2mup6*p3mup3*F(7369)
       Fa28(651) = mup5mup6*(p4mup3*F(7372)+p5mup3*F(7373)+2*(p2mup3*F
     -   (7371)-p3mup3*F(7374)))
       Fa28(652) = mup3mup5*(p1mup6*F(7376)+p6mup6*F(7379)+p2mup6*F(73
     -   80)+p5mup6*F(7382))+mup3mup6*(p5mup5*F(7377)+p6mup5*F(7378)+
     -   p1mup5*F(7381)+p2mup5*F(7383))
       Fa28(653) = 2*mup5mup6*p6mup3*F(7384)+mup3mup6*p4mup5*F(7385)
       Fa28(654) = p4mup6*(mup3mup5*F(7387)+8*(p5mup3*p5mup5*F(7391)-p
     -   2mup3*p5mup5*F(7394)-p3mup3*p5mup5*F(7396)+p4mup3*(p5mup5*F(
     -   7389)-2*(p4mup5*F(7388)+p2mup5*F(7393)+p3mup5*F(7399)))-2*(p
     -   2mup5*p5mup3*F(7395)+p3mup5*p5mup3*F(7397)-p4mup5*(-(p5mup3*
     -   F(7390))+p2mup3*F(7392)+p3mup3*F(7398))-p2mup3*p2mup5*F(7400
     -   )-p2mup5*p3mup3*F(7401)-p2mup3*p3mup5*F(7401)-p3mup3*p3mup5*
     -   F(7402))))
       Fa28(655) = p4mup6*(p1mup3*(-(p5mup5*F(7403))+2*(p2mup5*F(7404)
     -   +p1mup5*F(7406)+p3mup5*F(7407)+p4mup5*F(7408)))+2*p1mup5*(-(
     -   p2mup3*F(7405))-p3mup3*F(7409)+p4mup3*F(7410)+p5mup3*F(7411)
     -   ))
       Fa28(656) = -(p4mup3*p5mup5*p5mup6*F(7423))+p2mup6*(2*p4mup5*(p
     -   5mup3*F(7414)+p1mup3*F(7416)-p2mup3*F(7430)-p3mup3*F(7432))+
     -   p4mup3*(-(p5mup5*F(7415))+2*(p4mup5*F(7413)+p2mup5*F(7431)+p
     -   3mup5*F(7433))))+p3mup6*(-2*p4mup5*(p5mup3*F(7417)+p4mup3*F(
     -   7418)-p1mup3*F(7419)+p2mup3*F(7432)-p3mup3*F(7436))+p4mup3*(
     -   p5mup5*F(7420)+2*(p2mup5*F(7433)-p3mup5*F(7437))))+2*(p4mup3
     -   *p4mup5*p5mup6*F(7421)+p4mup5*p5mup3*p5mup6*F(7422)+p1mup3*p
     -   4mup5*p5mup6*F(7424)+p1mup6*p4mup3*p4mup5*F(7425)+p1mup6*p4m
     -   up5*p5mup3*F(7426)+p1mup3*p1mup6*p4mup5*F(7427)-p1mup6*p2mup
     -   3*p4mup5*F(7428)-p1mup6*p3mup3*p4mup5*F(7429)-p2mup3*p4mup5*
     -   p5mup6*F(7434)+p2mup5*p4mup3*p5mup6*F(7435)-p3mup3*p4mup5*p5
     -   mup6*F(7438)+p3mup5*p4mup3*p5mup6*F(7439))
       Fa28(657) = p1mup3*p5mup5*(p1mup6*F(7440)+2*p2mup6*F(7442))-p2m
     -   up3*(p1mup6*p5mup5*F(7441)+4*p1mup5*p2mup6*F(7445))+p1mup6*p
     -   4mup3*p5mup5*F(7451)+p1mup6*p5mup3*p5mup5*F(7455)+2*(p5mup3*
     -   p5mup5*p5mup6*F(7452)+p1mup3*(p5mup5*p5mup6*F(7444)+p3mup6*(
     -   p5mup5*F(7443)-2*(p2mup5*F(7448)+p1mup5*F(7457)+p3mup5*F(745
     -   8)))-2*(p1mup5*p2mup6*F(7447)+p2mup6*p3mup5*F(7448)+p2mup5*p
     -   5mup6*F(7450)+p1mup5*p5mup6*F(7453)+p3mup5*p5mup6*F(7454)+p2
     -   mup5*p2mup6*F(7456)+p1mup6*p2mup5*F(7461)+p1mup5*p1mup6*F(74
     -   62)+p1mup6*p3mup5*F(7463)))-p2mup3*p2mup6*p5mup5*F(7480)-p2m
     -   up3*p3mup6*p5mup5*F(7482)-p2mup3*p5mup5*p5mup6*F(7483))+4*(p
     -   1mup5*p2mup3*p3mup6*F(7446)+p1mup5*p2mup3*p5mup6*F(7449)+p1m
     -   up5*p1mup6*p2mup3*F(7460)-p1mup5*p2mup6*p4mup3*F(7464)-p1mup
     -   6*p2mup5*p4mup3*F(7465)-p1mup5*p2mup6*p5mup3*F(7467)-p1mup6*
     -   p2mup5*p5mup3*F(7468)-p1mup5*p3mup6*p4mup3*F(7469)-p1mup6*p3
     -   mup5*p4mup3*F(7470)-p1mup5*p3mup6*p5mup3*F(7473)-p1mup6*p3mu
     -   p5*p5mup3*F(7474)-p1mup5*p1mup6*p4mup3*F(7475)-p1mup5*p4mup3
     -   *p5mup6*F(7476)-p1mup5*p1mup6*p5mup3*F(7477)-p1mup5*p5mup3*p
     -   5mup6*F(7479)-p2mup5*p2mup6*p5mup3*F(7481)+p1mup6*p2mup3*p2m
     -   up5*F(7484)+p1mup6*p2mup3*p3mup5*F(7485)+p2mup3*p2mup5*p2mup
     -   6*F(7487)+p2mup3*p2mup6*p3mup5*F(7488)+p2mup3*p2mup5*p3mup6*
     -   F(7488)+p2mup3*p2mup5*p5mup6*F(7489)+p2mup3*p3mup5*p3mup6*F(
     -   7490)+p2mup3*p3mup5*p5mup6*F(7491))+p3mup3*(p1mup6*(p5mup5*F
     -   (7459)-4*(p1mup5*F(7471)-p2mup5*F(7485)))+4*(p1mup5*p2mup6*F
     -   (7466)+p1mup5*p3mup6*F(7472)+p1mup5*p5mup6*F(7478)+p2mup5*p2
     -   mup6*F(7488)+p2mup5*p3mup6*F(7490)+p3mup5*(p3mup6*F(5966)+p1
     -   mup6*F(7486)+p2mup6*F(7490))+p2mup5*p5mup6*F(7491)))
       Fa28(658) = p2mup6*p5mup5*(p3mup3*F(7482)-p5mup3*F(7492))+p3mup
     -   6*p5mup3*p5mup5*F(7493)+p3mup3*p3mup6*p5mup5*F(7496)+p3mup3*
     -   p5mup5*p5mup6*F(7498)+2*(p2mup5*p3mup6*p5mup3*F(7494)+p2mup5
     -   *p5mup3*p5mup6*F(7495)+p3mup5*(p5mup3*(p5mup6*F(6028)+p2mup6
     -   *F(7494)+p3mup6*F(7497))-p3mup3*p5mup6*F(7499)))
       Fa28(659) = p6mup3*(-(p5mup5*p5mup6*F(7508))-2*(p1mup6*(p2mup5*
     -   F(7501)+p3mup5*F(7503)+p4mup5*F(7505))-p2mup6*p4mup5*F(7506)
     -   +p3mup6*p4mup5*F(7507)-p1mup5*(p2mup6*F(7500)+p3mup6*F(7502)
     -   +p4mup6*F(7504)+p1mup6*F(7509))-p2mup5*p2mup6*F(7510)-p2mup5
     -   *p3mup6*F(7511)-p2mup5*p4mup6*F(7512)-p3mup5*p4mup6*F(7513))
     -   )
       Fa28(660) = -(p1mup6*p5mup5*p6mup3*F(7514))+p2mup6*(-4*(p3mup5*
     -   p6mup3*F(7511)-p1mup3*p6mup5*F(7517))+2*p5mup5*p6mup3*F(7520
     -   ))+2*(p4mup6*p5mup5*p6mup3*F(7524)+p3mup6*p5mup5*p6mup3*F(75
     -   26))+p2mup3*(p6mup6*(2*(p2mup5*F(4085)+p3mup5*F(4087)-p1mup5
     -   *F(7515)-p4mup5*F(7521))+p5mup5*F(7522))+4*p6mup5*(-(p1mup6*
     -   F(7516))+p2mup6*F(7529)+p3mup6*F(7530)+p4mup6*F(7531)+p5mup6
     -   *F(7533)))+4*(p1mup3*p1mup6*p6mup5*F(7518)-p1mup5*p5mup6*p6m
     -   up3*F(7519)-p4mup5*p4mup6*p6mup3*F(7523)+p4mup3*p4mup6*p6mup
     -   5*F(7525)-p4mup5*p5mup6*p6mup3*F(7527)+p5mup3*p5mup6*p6mup5*
     -   F(7528)-p2mup5*p5mup6*p6mup3*F(7532)-p3mup5*p3mup6*p6mup3*F(
     -   7534)+p3mup3*p3mup6*p6mup5*F(7535)-p3mup5*p5mup6*p6mup3*F(75
     -   36))
       Fa28(661) = p6mup5*(p1mup3*(p3mup6*F(7537)+p4mup6*F(7538)+p5mup
     -   6*F(7539))+p2mup6*(p3mup3*F(7530)+p4mup3*F(7540)+p5mup3*F(75
     -   41))+p1mup6*p4mup3*F(7542)+p3mup6*p4mup3*F(7543)+p4mup6*p5mu
     -   p3*F(7544)+p1mup6*p5mup3*F(7545)+p3mup6*p5mup3*F(7546)+p4mup
     -   3*p5mup6*F(7547)+p1mup6*p3mup3*F(7548)+p3mup3*p4mup6*F(7549)
     -   +p3mup3*p5mup6*F(7550))
       Fa28(662) = p6mup3*p6mup5*(p1mup6*F(7551)+2*(p2mup6*F(4301)+p4m
     -   up6*F(4310)+p3mup6*F(4316)+p5mup6*F(7552)))
       Fa28(663) = p6mup6*(p2mup3*p6mup5*F(7553)-p2mup5*p5mup3*F(7554)
     -   )
       Fa28(664) = p6mup6*(p5mup3*p5mup5*F(7566)+8*(p4mup5*p5mup3*F(75
     -   57)+p4mup3*p4mup5*F(7558)+p1mup3*(p2mup5*F(7556)+p3mup5*F(75
     -   59)+p1mup5*F(7560)+p4mup5*F(7562))-p2mup5*p4mup3*F(7563)+p1m
     -   up5*p4mup3*F(7565)+p3mup5*p4mup3*F(7567)+p3mup5*p5mup3*F(756
     -   9)+p1mup5*p5mup3*F(7570))-4*(p1mup3*p5mup5*F(7555)-p4mup3*p5
     -   mup5*F(7561)-p3mup3*(p5mup5*F(7568)+2*(p3mup5*F(4003)+p2mup5
     -   *F(4087)+p4mup5*F(7564)+p1mup5*F(7571)))))
       Fa28(665) = p6mup3*p6mup6*(p1mup5*F(7572)+p4mup5*F(7573)+p5mup5
     -   *F(7574)-p2mup5*F(7575)-p3mup5*F(7576))
       Fa28(666) = p6mup5*p6mup6*(-(p1mup3*F(7577))-p4mup3*F(7578)+p5m
     -   up3*F(7579)+2*p6mup3*F(7580)+p3mup3*F(7581))
       Fa28(667) = mup5mup6*(p5mup3*F(7582)+p4mup3*F(7584)-2*(p2mup3*F
     -   (7583)+p3mup3*F(7585)))
       Fa28(668) = mup3mup6*(p1mup5*F(7587)+p5mup5*F(7588)+p6mup5*F(75
     -   89)+p2mup5*F(7590))+2*mup3mup5*(-(p6mup6*F(7592))+2*(p1mup6*
     -   F(7591)+p2mup6*F(7593)+p5mup6*F(7594)))
       Fa28(669) = mup5mup6*p6mup3*F(7595)+mup3mup6*p4mup5*F(7596)
       Fa28(670) = p4mup6*(mup3mup5*F(7598)+2*(p2mup3*p5mup5*F(7602)+p
     -   3mup3*p5mup5*F(7604)+p4mup3*(p5mup5*F(7599)+2*(p2mup5*F(7601
     -   )+p3mup5*F(7606)+p4mup5*F(7607)))+p5mup3*(p5mup5*F(7600)-2*(
     -   p2mup5*F(7603)-p3mup5*F(7605)+p4mup5*F(7608)))+2*(p2mup3*p2m
     -   up5*F(7609)+p2mup5*p3mup3*F(7610)+p2mup3*p3mup5*F(7610)+p2mu
     -   p3*p4mup5*F(7611)+p3mup3*p4mup5*F(7612)+p3mup3*p3mup5*F(7613
     -   ))))
       Fa28(671) = p4mup6*(p1mup3*(p5mup5*F(7614)-2*(p2mup5*F(7616)+p3
     -   mup5*F(7617)+p4mup5*F(7618)+p1mup5*F(7620)))+2*p1mup5*(p5mup
     -   3*F(7615)+p4mup3*F(7619)+p2mup3*F(7621)+p3mup3*F(7622)))
       Fa28(672) = p4mup3*(p5mup5*(p5mup6*F(7626)+2*(p2mup6*F(7630)+p3
     -   mup6*F(7638)))+4*(-(p2mup5*p5mup6*F(7629))+p3mup5*p5mup6*F(7
     -   636)-p2mup5*p2mup6*F(7641)-p2mup6*p3mup5*F(7642)-p2mup5*p3mu
     -   p6*F(7642)+p4mup5*(p5mup6*F(7624)+p1mup6*F(7640)-p2mup6*F(76
     -   43)-p3mup6*F(7646))-p3mup5*p3mup6*F(7647)))+4*p4mup5*(p1mup3
     -   *p2mup6*F(7631)+p3mup3*p5mup6*F(7633)+p1mup3*p1mup6*F(7634)+
     -   p1mup3*p3mup6*F(7635)+p1mup3*p5mup6*F(7637)-p1mup6*p3mup3*F(
     -   7639)+p5mup3*(p5mup6*F(7625)+p1mup6*F(7627)+p2mup6*F(7644)+p
     -   3mup6*F(7645))-p2mup6*p3mup3*F(9678)-p2mup3*(p5mup6*F(7628)+
     -   p1mup6*F(7632)+p2mup6*F(9674)+p3mup6*F(9678))-p3mup3*p3mup6*
     -   F(9682))
       Return
       End
