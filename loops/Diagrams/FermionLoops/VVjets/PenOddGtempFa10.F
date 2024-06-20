c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenOddG10(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup
     -   1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3
     -   mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3
     -   ,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,
     -   mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa10)
       IMPLICIT NONE
      Complex*16   p1mup1, p1mup2, p1mup3, p1mup4, p1mup5, p2mup1, 
     -          p2mup2, p2mup3, p2mup4, p2mup5, p3mup1, p3mup2, 
     -          p3mup3, p3mup4, p3mup5, p4mup1, p4mup2, p4mup3, 
     -          p4mup4, p4mup5, p5mup1, p5mup2, p5mup3, p5mup4, p5mup5
       Complex*16   mup1mup2, mup1mup3, mup1mup4, mup1mup5, mup2mup3, 
     -          mup2mup4, mup2mup5, mup3mup4, mup3mup5, mup4mup5
        common/FPenOddG/F
      Complex*16 EpsF(126)
       common/EPSFPenOddG/EPSF
       COMMON/InvariantsPenOddG/p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,s
     -   45,s234,m1,m1sq,m3,m3sq,m4,m4sq
       Real*8   p1sq 
        Real*8   p2sq 
       Real*8   p3sq 
       Real*8   p4sq 
       Real*8   p5sq 
       Real*8   s12 
       Real*8   s23 
       Real*8   s34 
       Real*8   s45 
       Real*8   s234,m1,m1sq,m3,m3sq,m4,m4sq
      Complex*16 Fa(4390),F(7459)
      Real*8 P(1989) 
        COMMON/PPenOddG/P
       Complex*16 Fa10(2629:2920)
       COMMON/FaPenOddG/Fa
       Fa10(2629) = p1mup2*p3mup3*F(7219)+p1mup3*(-2*(p2mup2*F(675)-p1
     -   mup2*F(7146))+p3mup2*F(7220))+2*(p1mup2*p2mup3*F(7147)+p2mup
     -   3*p3mup2*F(7221)+p2mup2*p3mup3*F(7224))+4*(p2mup2*p2mup3*F(3
     -   878)+p3mup2*p3mup3*F(7225))
       Fa10(2630) = 4*p2mup2*p2mup4*F(7021)-2*(p2mup4*p4mup2*F(1329)-p
     -   4mup2*p4mup4*F(1337)-p1mup2*(2*(p1mup4*F(634)+p2mup4*F(635))
     -   +p4mup4*F(1326)+p3mup4*F(1838))+p2mup4*p3mup2*F(1840)-p3mup2
     -   *p3mup4*F(1849)+p1mup4*(p4mup2*F(1327)+p3mup2*F(1839)-2*p2mu
     -   p2*F(7022))-p2mup2*p4mup4*F(7023)-p2mup2*p3mup4*F(7109))+p3m
     -   up2*p4mup4*F(7257)+p3mup4*p4mup2*F(7258)
       Fa10(2631) = 4*(p4mup2*p4mup4*F(1292)+p3mup2*p3mup4*F(1876)+p2m
     -   up4*p3mup2*F(4199))+p1mup2*(2*(p1mup4*F(1285)+p2mup4*F(1286)
     -   -p3mup4*F(1871))-p4mup4*F(7212))+p1mup4*p4mup2*F(7213)+2*(p1
     -   mup4*p3mup2*F(1872)+p2mup4*p4mup2*F(7259)+p3mup2*p4mup4*F(72
     -   60)+p3mup4*p4mup2*F(7261)+p2mup2*(p1mup4*F(1289)+2*(p2mup4*F
     -   (1288)+p3mup4*F(4200))+p4mup4*F(7262)))
       Fa10(2632) = p1mup2*(2*(p1mup4*F(667)+p2mup4*F(669)-p4mup4*F(13
     -   92))+p3mup4*F(7202))+p1mup4*(2*(p2mup2*F(675)+p4mup2*F(1393)
     -   )+p3mup2*F(7203))+2*(p2mup4*p3mup2*F(7204)+p2mup2*p3mup4*F(7
     -   205)+2*(p2mup2*p2mup4*F(1390)+p2mup2*p4mup4*F(1396)+p4mup2*p
     -   4mup4*F(1397)+p2mup4*p4mup2*F(4112)+p3mup2*p3mup4*F(7206))+p
     -   3mup2*p4mup4*F(7263)+p3mup4*p4mup2*F(7264))
       Fa10(2633) = p2mup2*p4mup4*F(7265)+p2mup4*p4mup2*F(7266)+2*(p2m
     -   up2*p2mup4*F(1553)+p2mup4*p3mup2*F(1940)+p2mup2*p3mup4*F(719
     -   1)+p1mup2*(p2mup4*F(1551)+2*(p1mup4*F(886)+p3mup4*F(7192))-p
     -   4mup4*F(7227))+p4mup2*p4mup4*F(7267)+p1mup4*(p2mup2*F(1552)+
     -   2*(p3mup2*F(7192)+p4mup2*F(7268)))+p3mup2*p4mup4*F(7269))+4*
     -   (p3mup2*p3mup4*F(1941)+p3mup4*p4mup2*F(7270))
       Fa10(2634) = 4*((p3mup4*p4mup2+p3mup2*p4mup4)*F(913)+p4mup2*p4m
     -   up4*F(914)+p2mup2*p2mup4*F(2731))+p1mup2*(2*p1mup4*F(908)+p2
     -   mup4*F(7251)+4*(p4mup4*F(3686)+p3mup4*F(7254)))+p1mup4*(p2mu
     -   p2*F(7252)+2*(p3mup2*F(7255)+p4mup2*F(7271)))+2*(p2mup4*p3mu
     -   p2*F(7253)+p2mup2*p3mup4*F(7253)+2*p3mup2*p3mup4*F(7256)+(p2
     -   mup4*p4mup2+p2mup2*p4mup4)*F(7272))
       Fa10(2635) = 4*(p4mup2*p4mup4*F(1018)+p3mup2*p4mup4*F(7087))+p2
     -   mup4*p3mup2*F(7245)+p2mup2*p3mup4*F(7248)+2*(p2mup2*p2mup4*F
     -   (1518)+p2mup4*p4mup2*F(1522)+p1mup4*(p2mup2*F(963)+2*(p3mup2
     -   *F(7084)+p4mup2*F(7085)))+p2mup2*p4mup4*F(7086)+p1mup2*(p2mu
     -   p4*F(961)+2*(p1mup4*F(960)+p4mup4*F(7085))+p3mup4*F(7249))+p
     -   3mup2*p3mup4*F(7250)+p3mup4*p4mup2*F(7273))
       Fa10(2636) = -(p1mup2*p2mup3*F(1551))-p2mup2*p2mup3*F(1553)+p2m
     -   up2*p3mup3*F(2346)+p1mup3*(2*(p3mup2*F(284)-p1mup2*F(886))+p
     -   2mup2*F(7226))+p2mup3*p3mup2*F(7231)+2*(p1mup2*p3mup3*F(284)
     -   +p3mup2*p3mup3*F(7232))
       Fa10(2637) = -2*p1mup2*(2*p1mup3*F(960)+p2mup3*F(961)-p3mup3*F(
     -   7243))+p2mup2*(-2*(p2mup3*F(1518)-p1mup3*F(7186))+p3mup3*F(7
     -   244))+p3mup2*(4*p1mup3*F(1520)-p2mup3*F(7245)+2*p3mup3*F(724
     -   6))
       Fa10(2638) = p1mup2*(4*p3mup3*F(3685)+2*p1mup3*F(7233)+p2mup3*F
     -   (7234))+p1mup3*p2mup2*F(7235)+2*(p1mup3*p3mup2*F(7236)+p2mup
     -   3*p3mup2*F(7239)+p2mup2*p3mup3*F(7239))+4*(p2mup2*p2mup3*F(7
     -   241)+p3mup2*p3mup3*F(7242))
       Fa10(2639) = p1mup2*(2*(p1mup3*F(7110)+p2mup3*F(7111))+p3mup3*F
     -   (7112))+2*p2mup2*(p3mup3*F(7114)+p1mup3*F(7116)+2*p2mup3*F(7
     -   117))+p3mup2*(4*p3mup3*F(7118)-p1mup3*F(7274)-2*p2mup3*F(727
     -   5))
       Fa10(2640) = 4*p2mup2*p2mup4*F(7021)+p1mup4*(4*p2mup2*F(7022)-2
     -   *(p4mup2*F(1327)-p3mup2*F(7120)))-2*(p2mup4*p4mup2*F(1329)-p
     -   4mup2*p4mup4*F(1337)-p1mup2*(2*(p1mup4*F(634)+p2mup4*F(635))
     -   +p4mup4*F(1326)+p3mup4*F(1838))-p2mup2*p4mup4*F(7023)-p2mup2
     -   *p3mup4*F(7109)-p3mup2*p3mup4*F(7119)-p2mup4*p3mup2*F(7121))
     -   +p3mup2*p4mup4*F(7122)+p3mup4*p4mup2*F(7123)
       Fa10(2641) = p1mup2*(2*(p1mup4*F(667)+p2mup4*F(669)-p4mup4*F(13
     -   92))-p3mup4*F(7144))+p1mup4*(2*(p2mup2*F(675)+p4mup2*F(1393)
     -   )-p3mup2*F(7145))+4*(p2mup2*p2mup4*F(1390)+p2mup2*p4mup4*F(1
     -   396)+p4mup2*p4mup4*F(1397)+p2mup4*p4mup2*F(4112)-p3mup2*p3mu
     -   p4*F(7150))-2*(p2mup2*p3mup4*F(7148)+p3mup2*p4mup4*F(7151)+p
     -   3mup4*p4mup2*F(7152)-p2mup4*p3mup2*F(7276))
       Fa10(2642) = 4*(p4mup2*p4mup4*F(1018)+p3mup2*p4mup4*F(7087))-p2
     -   mup2*p3mup4*F(7184)+p2mup4*p3mup2*F(7185)+2*(p2mup2*p2mup4*F
     -   (1518)+p2mup4*p4mup2*F(1522)+p1mup4*(p2mup2*F(963)+2*(p3mup2
     -   *F(7084)+p4mup2*F(7085)))+p2mup2*p4mup4*F(7086)+p1mup2*(p2mu
     -   p4*F(961)+2*(p1mup4*F(960)+p4mup4*F(7085))-p3mup4*F(7183))-p
     -   3mup2*p3mup4*F(7187)-p3mup4*p4mup2*F(7188))
       Fa10(2643) = p1mup2*(2*(p1mup3*F(6987)+p2mup3*F(6988))+p3mup3*F
     -   (6989))-p2mup3*p3mup2*F(6991)-p1mup3*(p3mup2*F(6990)+2*p2mup
     -   2*F(6993))-2*p2mup2*p2mup3*F(6994)+p2mup2*p3mup3*F(7277)+p3m
     -   up2*p3mup3*F(7278)
       Fa10(2644) = p1mup3*p2mup2*F(6998)+p1mup2*(p1mup3*F(6996)+p2mup
     -   3*F(6997)-p3mup3*F(7001))+p1mup3*p3mup2*F(7279)+2*(p2mup2*p2
     -   mup3*F(6999)-p2mup3*p3mup2*F(7000)+p3mup2*p3mup3*F(7280)+p2m
     -   up2*p3mup3*F(7281))
       Fa10(2645) = 4*(p3mup2*p3mup3*F(7008)+p1mup2*p1mup3*F(7009)+p2m
     -   up3*(p3mup2*F(7007)+p1mup2*F(7010))+p2mup2*(p3mup3*F(7005)+2
     -   *p2mup3*F(7006)+p1mup3*F(7013)))+p1mup2*p3mup3*F(7282)+p1mup
     -   3*p3mup2*F(7283)
       Fa10(2646) = 4*(p1mup4*p3mup2*F(7015)+p2mup4*p3mup2*F(7016)+m1s
     -   q*(2*(p1mup4*p2mup2*F(641)+p2mup2*p2mup4*F(1281))+p1mup4*p4m
     -   up2*F(1327)+p2mup4*p4mup2*F(1329)+p2mup2*p4mup4*F(1680)+p1mu
     -   p2*(-2*(p1mup4*F(634)-p2mup4*F(7124))+p4mup4*F(7126))+p4mup2
     -   *p4mup4*F(7127)))+p1mup2*p3mup4*F(7284)+p3mup2*p3mup4*F(7285
     -   )+p2mup2*p3mup4*F(7286)+p3mup2*p4mup4*F(7287)+p3mup4*p4mup2*
     -   F(7288)
       Fa10(2647) = -2*p1mup4*p4mup2*F(7039)+4*(p4mup2*p4mup4*F(7036)+
     -   p1mup2*p1mup4*F(7291))+p3mup2*p3mup4*F(7293)+p3mup2*p4mup4*F
     -   (7294)+p1mup2*p3mup4*F(7295)+p2mup4*(8*p3mup2*F(7030)+p1mup2
     -   *F(7292)+p4mup2*F(7296))+p3mup4*p4mup2*F(7297)+p2mup2*(8*p2m
     -   up4*F(7024)+p3mup4*F(7289)+p4mup4*F(7290)+p1mup4*F(7298))+p1
     -   mup4*p3mup2*F(7299)+p1mup2*p4mup4*F(7300)
       Fa10(2648) = 2*p3mup2*(2*(p3mup4*F(7041)+p2mup4*F(7042)+p4mup4*
     -   F(7043))+p1mup4*F(7046))+4*m1sq*(p1mup2*p4mup4*F(1392)+2*(p2
     -   mup4*p4mup2*F(1391)+p2mup2*p2mup4*F(3878)+p2mup2*p4mup4*F(38
     -   79)+p4mup2*p4mup4*F(3880))+p1mup2*p2mup4*F(7147)+p1mup4*(-(p
     -   2mup2*F(675))+p1mup2*F(7146)+p4mup2*F(7153)))+p3mup4*(p2mup2
     -   *F(7301)+p4mup2*F(7302)+p1mup2*F(7303))
       Fa10(2649) = p2mup2*p2mup4*F(7306)+p2mup2*p3mup4*F(7308)+p3mup2
     -   *p3mup4*F(7309)+p1mup4*(4*p2mup2*F(7049)+p1mup2*F(7304)+8*(p
     -   4mup2*F(7061)+p3mup2*F(7310)))+p2mup4*p3mup2*F(7311)+p1mup2*
     -   (p2mup4*F(7305)+p3mup4*F(7307)+p4mup4*F(7312))+p2mup2*p4mup4
     -   *F(7313)+p2mup4*p4mup2*F(7314)+p3mup2*p4mup4*F(7315)+p4mup2*
     -   p4mup4*F(7316)+p3mup4*p4mup2*F(7317)
       Fa10(2650) = -2*p2mup4*p3mup2*F(7081)+4*m1sq*(-(p1mup2*p2mup4*F
     -   (961))-p2mup2*p2mup4*F(1518)+p2mup2*p4mup4*F(1777)+p1mup4*(2
     -   *(p4mup2*F(336)-p1mup2*F(960)+p3mup2*F(1520))+p2mup2*F(7186)
     -   )+p2mup4*p4mup2*F(7189)+2*(p1mup2*p4mup4*F(336)+p3mup2*p4mup
     -   4*F(1523)+p4mup2*p4mup4*F(7190)))+p3mup4*(p1mup2*F(7318)+p2m
     -   up2*F(7319)+p3mup2*F(7320)+p4mup2*F(7321))
       Fa10(2651) = p1mup3*p2mup2*F(7090)+p2mup2*p2mup3*F(7091)+p1mup2
     -   *(p2mup3*F(7089)+2*(p1mup3*F(7088)-p3mup3*F(7094)))+p2mup2*p
     -   3mup3*F(7322)+p2mup3*p3mup2*F(7323)-2*(p1mup3*p3mup2*F(7094)
     -   -p3mup2*p3mup3*F(7324))
       Fa10(2652) = 4*(p1mup3*p4mup2*F(7325)+p2mup3*p4mup2*F(7326))+p3
     -   mup3*p4mup2*F(7328)+p4mup3*(-2*p3mup2*F(7327)+p1mup2*F(7329)
     -   +4*p4mup2*F(7330)+p2mup2*F(7331))
       Fa10(2653) = -4*p4mup2*(p2mup3*F(7333)+p3mup3*F(7335))+2*(p4mup
     -   3*(2*(p2mup2*F(7332)+p3mup2*F(7334)+p4mup2*F(7336))+p1mup2*F
     -   (7337))-p1mup3*p4mup2*F(7338))
       Fa10(2654) = -(p3mup2*p4mup3*F(7340))+p3mup3*p4mup2*F(7341)+2*(
     -   p2mup3*p4mup2*F(7339)+p4mup2*p4mup3*F(7342)+p2mup2*p4mup3*F(
     -   7343))+p1mup2*p4mup3*F(7344)+p1mup3*p4mup2*F(7345)
       Fa10(2655) = 4*m4sq*p4mup2*(p1mup3*F(1937)+p3mup3*F(1942))+p2mu
     -   p3*p4mup2*F(7348)+p4mup3*(p2mup2*F(7347)+2*(p1mup2*F(7346)+p
     -   3mup2*F(7349)+p4mup2*F(7350)))
       Fa10(2656) = p2mup2*p4mup3*F(7352)+p2mup3*p4mup2*F(7353)+p3mup3
     -   *p4mup2*F(7354)+2*(p1mup3*p4mup2*F(7351)+p1mup2*p4mup3*F(735
     -   1)+p3mup2*p4mup3*F(7355)+p4mup2*p4mup3*F(7356))
       Fa10(2657) = p1mup3*p4mup2*F(7078)+p2mup2*p4mup3*F(7357)+p2mup3
     -   *p4mup2*F(7358)+2*(p1mup2*p4mup3*F(7077)+p4mup2*p4mup3*F(735
     -   9)+(p3mup3*p4mup2+p3mup2*p4mup3)*F(7360))
       Fa10(2658) = 2*(p1mup2*p4mup3*F(1326)+p2mup2*p4mup3*F(7023))+p3
     -   mup2*p4mup3*F(7363)-p4mup2*(p3mup3*F(7364)+2*(p1mup3*F(7361)
     -   +p2mup3*F(7362)+p4mup3*F(7365)))
       Fa10(2659) = -(p4mup3*(p1mup2*F(7366)+2*(p2mup2*F(7368)+p3mup2*
     -   F(7370))))+p4mup2*(4*p4mup3*F(1292)+p1mup3*F(7367)+2*(p2mup3
     -   *F(7369)+p3mup3*F(7371)))
       Fa10(2660) = p1mup2*p4mup3*F(7372)+p3mup2*p4mup3*F(7375)+p4mup2
     -   *(p1mup3*F(7373)+p3mup3*F(7376)+2*(p2mup3*F(7374)+p4mup3*F(7
     -   377)))+2*p2mup2*p4mup3*F(7378)
       Fa10(2661) = -4*(p1mup3*p4mup2*F(7268)+p3mup3*p4mup2*F(7270))+p
     -   2mup3*p4mup2*F(7381)+p4mup3*(p2mup2*F(7380)+2*(p1mup2*F(7379
     -   )+p3mup2*F(7382)+p4mup2*F(7383)))
       Fa10(2662) = p2mup2*p4mup3*F(7385)+p2mup3*p4mup2*F(7386)+p3mup3
     -   *p4mup2*F(7387)+2*(p1mup3*p4mup2*F(7384)+p1mup2*p4mup3*F(738
     -   4)+p3mup2*p4mup3*F(7388)+p4mup2*p4mup3*F(7389))
       Fa10(2663) = p4mup2*(p1mup3*F(7178)+p2mup3*F(7179)+2*p3mup3*F(7
     -   180))+p4mup3*(p2mup2*F(7179)+2*(p1mup2*F(7177)+p3mup2*F(7180
     -   )+p4mup2*F(7182)))
       Fa10(2664) = -2*p4mup2*(p1mup3*F(1327)+p2mup3*F(1329))+p4mup3*(
     -   2*(p1mup2*F(1326)+p4mup2*F(1337)+p2mup2*F(7023))+p3mup2*F(72
     -   57))+p3mup3*p4mup2*F(7258)
       Fa10(2665) = p4mup2*(4*p4mup3*F(1292)+p1mup3*F(7213)+2*(p2mup3*
     -   F(7259)+p3mup3*F(7261)))+p4mup3*(-(p1mup2*F(7212))+2*(p3mup2
     -   *F(7260)+p2mup2*F(7262)))
       Fa10(2666) = -(p1mup2*p4mup3*F(1392))+p1mup3*p4mup2*F(1393)+2*(
     -   p2mup2*p4mup3*F(1396)+p4mup2*p4mup3*F(1397)+p2mup3*p4mup2*F(
     -   4112))+p3mup2*p4mup3*F(7263)+p3mup3*p4mup2*F(7264)
       Fa10(2667) = p2mup2*p4mup3*F(7265)+p2mup3*p4mup2*F(7266)-2*(p1m
     -   up2*p4mup3*F(7227)-p4mup2*p4mup3*F(7267)-p3mup2*p4mup3*F(726
     -   9))+4*(p1mup3*p4mup2*F(7268)+p3mup3*p4mup2*F(7270))
       Fa10(2668) = p4mup3*(p2mup2*F(7086)+2*(p1mup2*F(7085)+p3mup2*F(
     -   7087)))+p4mup2*(p2mup3*F(1522)+2*(p4mup3*F(1018)+p1mup3*F(70
     -   85))+p3mup3*F(7273))
       Fa10(2669) = 2*(p3mup3*p4mup2*F(913)+p3mup2*p4mup3*F(913)+p4mup
     -   2*p4mup3*F(914)+p1mup2*p4mup3*F(3686))+p1mup3*p4mup2*F(7271)
     -   +(p2mup3*p4mup2+p2mup2*p4mup3)*F(7272)
       Fa10(2670) = 2*(p1mup3*p4mup2*F(1327)+p2mup3*p4mup2*F(1329))+p4
     -   mup3*(2*(p2mup2*F(1680)+p1mup2*F(7126)+p4mup2*F(7127))+p3mup
     -   2*F(7207))+p3mup3*p4mup2*F(7208)
       Fa10(2671) = p4mup2*(4*p4mup3*F(1667)-p1mup3*F(7213)+2*(p2mup3*
     -   F(7214)+p3mup3*F(7216)))+p4mup3*(p1mup2*F(7212)+2*(p3mup2*F(
     -   7215)+p2mup2*F(7218)))
       Fa10(2672) = p1mup2*p4mup3*F(1392)+2*(p2mup3*p4mup2*F(1391)+p2m
     -   up2*p4mup3*F(3879)+p4mup2*p4mup3*F(3880))+p1mup3*p4mup2*F(71
     -   53)+p3mup2*p4mup3*F(7222)+p3mup3*p4mup2*F(7223)
       Fa10(2673) = p4mup2*(2*p4mup3*F(1233)+4*(p1mup3*F(1937)+p3mup3*
     -   F(1942))+p2mup3*F(7229))+p4mup3*(p2mup2*F(7228)+2*(p1mup2*F(
     -   7227)+p3mup2*F(7230)))
       Fa10(2674) = 2*p1mup3*p4mup2*F(336)+p2mup3*p4mup2*F(7189)+p4mup
     -   3*(p2mup2*F(1777)+2*(p1mup2*F(336)+p3mup2*F(1523)+p4mup2*F(7
     -   190)))+p3mup3*p4mup2*F(7247)
       Fa10(2675) = 2*(p3mup3*p4mup2*F(3772)+p3mup2*p4mup3*F(3772)+p4m
     -   up2*p4mup3*F(3773)+p1mup2*p4mup3*F(7237))+p1mup3*p4mup2*F(72
     -   38)+(p2mup3*p4mup2+p2mup2*p4mup3)*F(7240)
       Fa10(2676) = 2*(p2mup2*p4mup3*F(1680)+p1mup2*p4mup3*F(7126))-p3
     -   mup2*p4mup3*F(7363)+p4mup2*(p3mup3*F(7364)+2*(p1mup3*F(7361)
     -   +p2mup3*F(7362)+p4mup3*F(7365)))
       Fa10(2677) = p4mup3*(p1mup2*F(7366)+2*(p2mup2*F(7368)+p3mup2*F(
     -   7370)))+p4mup2*(4*p4mup3*F(1667)-p1mup3*F(7367)-2*(p2mup3*F(
     -   7369)+p3mup3*F(7371)))
       Fa10(2678) = p4mup2*(p1mup3*F(7373)+2*p2mup3*F(7374)+p3mup3*F(7
     -   376))+p4mup3*(p1mup2*F(7372)+2*(p4mup2*F(7377)+p2mup2*F(7378
     -   ))-p3mup2*F(7390))
       Fa10(2679) = 4*(p1mup3*p4mup2*F(1937)+p3mup3*p4mup2*F(1942))+p2
     -   mup3*p4mup2*F(7381)+p4mup3*(p2mup2*F(7380)+2*(p1mup2*F(7379)
     -   +p3mup2*F(7382)+p4mup2*F(7383)))
       Fa10(2680) = 2*p3mup2*p4mup3*F(7327)+p3mup3*p4mup2*F(7391)-4*(p
     -   1mup3*p4mup2*F(7325)+p2mup3*p4mup2*F(7326)+p4mup2*p4mup3*F(7
     -   330)-p1mup2*p4mup3*F(7392))+p2mup2*p4mup3*F(7393)
       Fa10(2681) = 2*(p4mup2*(2*(p2mup3*F(7333)+p3mup3*F(7335))+p1mup
     -   3*F(7338))-p4mup3*(p1mup2*F(7337)+2*(p2mup2*F(7332)+p3mup2*F
     -   (7334)-p4mup2*F(7394))))
       Fa10(2682) = p3mup3*p4mup2*F(7341)+2*(p2mup3*p4mup2*F(7339)+p4m
     -   up2*p4mup3*F(7342)+p2mup2*p4mup3*F(7343))+p1mup2*p4mup3*F(73
     -   44)+p1mup3*p4mup2*F(7345)+p3mup2*p4mup3*F(7395)
       Fa10(2683) = 8*m4sq*p4mup2*(p1mup3*F(7268)+p3mup3*F(7270))-2*(p
     -   2mup3*p4mup2*F(7348)+p4mup3*(p2mup2*F(7347)+2*(p1mup2*F(7346
     -   )+p3mup2*F(7349)+p4mup2*F(7350))))
       Fa10(2684) = p3mup1*(p3mup4*F(6865)+p4mup4*F(6867)+4*(p1mup4*F(
     -   6863)+p2mup4*F(7396)))
       Fa10(2685) = 2*p3mup1*(p1mup4*F(6874)+2*(p3mup4*F(6871)+p4mup4*
     -   F(6872)+p2mup4*F(7397)))
       Fa10(2686) = p2mup1*p3mup4*F(6985)-2*p2mup4*p3mup1*F(7398)+4*m1
     -   sq*p2mup1*(p4mup4*F(1016)+p1mup4*F(6927)+p2mup4*F(7399))
       Fa10(2687) = p2mup1*(p1mup4*F(6874)+2*(p3mup4*F(6871)+p4mup4*F(
     -   6872)+p2mup4*F(7397)))
       Fa10(2688) = p2mup1*(p3mup4*F(6982)+p4mup4*F(6984)+4*(p1mup4*F(
     -   6981)+p2mup4*F(7400)))
       Fa10(2689) = mup2mup3*(p3mup4*F(6865)+p4mup4*F(6867)+4*(p1mup4*
     -   F(6863)+p2mup4*F(7396)))+mup2mup4*(2*p4mup3*F(6868)-4*(p1mup
     -   3*F(6862)-p3mup3*F(6866)+p2mup3*F(7402)))+4*mup3mup4*(p1mup2
     -   *F(7401)+p2mup2*F(7403)+p3mup2*F(7404))-2*mup3mup4*p4mup2*F(
     -   7405)
       Fa10(2690) = mup2mup3*(p1mup4*F(6874)+2*(p3mup4*F(6871)+p4mup4*
     -   F(6872)+p2mup4*F(7397)))-mup2mup4*(p1mup3*F(6874)+2*(p3mup3*
     -   F(6870)+p4mup3*F(6873)+p2mup3*F(7407)))+mup3mup4*(p1mup2*F(6
     -   874)+2*(p2mup2*F(7406)+p3mup2*F(7408)+p4mup2*F(7409)))
       Fa10(2691) = mup2mup3*p3mup4*F(6877)+4*m1sq*(mup3mup4*(p1mup2*F
     -   (1325)+p2mup2*F(1623)+p4mup2*F(2718))+mup2mup3*(p1mup4*F(327
     -   )+p2mup4*F(2414)+p4mup4*F(6880)))+mup2mup4*(p3mup3*F(6878)-4
     -   *(p1mup3*F(6876)+p4mup3*F(6879)+p2mup3*F(7410)))+2*mup3mup4*
     -   p3mup2*F(7411)
       Fa10(2692) = (mup2mup4*p2mup3-mup2mup3*p2mup4)*F(7398)+mup3mup4
     -   *(p2mup2*F(7413)+2*(p1mup2*F(7412)+p3mup2*F(7414)+p4mup2*F(7
     -   415)))
       Fa10(2693) = p3mup1*(p4mup4*F(6919)+2*(p1mup4*F(6917)+p3mup4*F(
     -   6918)+p2mup4*F(7416)))
       Fa10(2694) = p3mup1*(p1mup4*F(6921)+2*(p3mup4*F(6922)+p4mup4*F(
     -   6923)+p2mup4*F(7417)))
       Fa10(2695) = p2mup1*(2*(p1mup4*F(327)+p2mup4*F(2414)+p4mup4*F(6
     -   880))-p3mup4*F(6925))+p2mup4*p3mup1*F(7418)
       Fa10(2696) = p2mup1*(p1mup4*F(6921)+2*(p3mup4*F(6922)+p4mup4*F(
     -   6923)+p2mup4*F(7417)))
       Fa10(2697) = p2mup1*(p4mup4*F(6919)+2*(p1mup4*F(6917)+p3mup4*F(
     -   6918)+p2mup4*F(7416)))
       Fa10(2698) = mup2mup4*(2*(p1mup3*F(271)+p3mup3*F(1235)+p2mup3*F
     -   (2710))-p4mup3*F(6920))+mup2mup3*(p4mup4*F(6919)+2*(p1mup4*F
     -   (6917)+p3mup4*F(6918)+p2mup4*F(7416)))+mup3mup4*(2*(p2mup2*F
     -   (1335)+p1mup2*F(1646)+p3mup2*F(5681))+p4mup2*F(7419))
       Fa10(2699) = mup3mup4*p1mup2*F(6921)-mup2mup4*p1mup3*F(6921)+mu
     -   p2mup3*p1mup4*F(6921)+2*(mup3mup4*p3mup2*F(6922)-mup2mup4*p3
     -   mup3*F(6922)+mup2mup3*p3mup4*F(6922)+mup3mup4*p4mup2*F(6923)
     -   -mup2mup4*p4mup3*F(6923)+mup2mup3*p4mup4*F(6923)+(mup3mup4*p
     -   2mup2-mup2mup4*p2mup3+mup2mup3*p2mup4)*F(7417))
       Fa10(2700) = mup2mup3*(p3mup4*F(6925)+2*(p4mup4*F(1016)+p1mup4*
     -   F(6927)+p2mup4*F(7399)))+mup2mup4*(p3mup3*F(6926)+2*(p1mup3*
     -   F(6924)+p4mup3*F(6928)+p2mup3*F(7420)))+mup3mup4*(2*(p2mup2*
     -   F(1335)+p1mup2*F(1646)+p4mup2*F(6355))-p3mup2*F(7421))
       Fa10(2701) = (mup2mup4*p2mup3-mup2mup3*p2mup4)*F(7418)+mup3mup4
     -   *(2*(p3mup2*F(6922)+p4mup2*F(6923)+p1mup2*F(7422))+p2mup2*F(
     -   7423))
       Fa10(2702) = p3mup1*(p4mup4*F(6941)+2*(p1mup4*F(885)+p3mup4*F(8
     -   94)+p2mup4*F(7424)))
       Fa10(2703) = p3mup1*(p1mup4*F(6932)+2*(p3mup4*F(6937)+p4mup4*F(
     -   6943)+p2mup4*F(7425)))
       Fa10(2704) = p2mup4*p3mup1*F(1094)+p2mup1*(2*(p1mup4*F(327)+p2m
     -   up4*F(2414)+p4mup4*F(6880))+p3mup4*F(6934))
       Fa10(2705) = p2mup1*(p1mup4*F(6931)+2*(p3mup4*F(6938)+p4mup4*F(
     -   6944)+p2mup4*F(7426)))
       Fa10(2706) = p2mup1*(2*(p1mup4*F(271)+p3mup4*F(1235)+p2mup4*F(2
     -   710))+p4mup4*F(6942))
       Fa10(2707) = mup2mup4*(2*(p1mup3*F(271)+p3mup3*F(1235)+p2mup3*F
     -   (2710))+p4mup3*F(6942))+mup2mup3*(p4mup4*F(6941)+2*(p1mup4*F
     -   (885)+p3mup4*F(894)+p2mup4*F(7424)))+mup3mup4*(2*(p2mup2*F(1
     -   335)+p1mup2*F(1646)+p3mup2*F(5681))+p4mup2*F(7427))
       Fa10(2708) = mup3mup4*p1mup2*F(6931)+mup2mup3*p1mup4*F(6931)+mu
     -   p2mup4*(p1mup3*F(6932)+2*(p3mup3*F(6937)+p4mup3*F(6943)+p2mu
     -   p3*F(7425)))+2*(mup3mup4*p3mup2*F(6938)+mup2mup3*p3mup4*F(69
     -   38)+mup3mup4*p4mup2*F(6944)+mup2mup3*p4mup4*F(6944)+(mup3mup
     -   4*p2mup2+mup2mup3*p2mup4)*F(7426))
       Fa10(2709) = mup2mup4*(2*(p1mup3*F(327)+p2mup3*F(2414)+p4mup3*F
     -   (6880))+p3mup3*F(6934))+mup2mup3*(p3mup4*F(6933)+2*(p4mup4*F
     -   (1016)+p1mup4*F(6927)+p2mup4*F(7399)))+mup3mup4*(2*(p2mup2*F
     -   (1335)+p1mup2*F(1646)+p4mup2*F(6355))+p3mup2*F(7428))
       Fa10(2710) = mup2mup3*p2mup4*F(1093)+mup2mup4*p2mup3*F(1094)+mu
     -   p3mup4*(2*(p3mup2*F(6938)+p4mup2*F(6944)+p1mup2*F(7429))+p2m
     -   up2*F(7430))
       Fa10(2711) = p3mup1*(2*(p1mup4*F(271)+p3mup4*F(1235)+p2mup4*F(2
     -   710))+p4mup4*F(6942))
       Fa10(2712) = p3mup1*(p1mup4*F(6931)+2*(p3mup4*F(6938)+p4mup4*F(
     -   6944)+p2mup4*F(7426)))
       Fa10(2713) = p2mup4*p3mup1*F(1093)+p2mup1*(p3mup4*F(6933)+2*(p4
     -   mup4*F(1016)+p1mup4*F(6927)+p2mup4*F(7399)))
       Fa10(2714) = p2mup1*(p1mup4*F(6932)+2*(p3mup4*F(6937)+p4mup4*F(
     -   6943)+p2mup4*F(7425)))
       Fa10(2715) = p2mup1*(p4mup4*F(6941)+2*(p1mup4*F(885)+p3mup4*F(8
     -   94)+p2mup4*F(7424)))
       Fa10(2716) = mup2mup3*(2*(p1mup4*F(271)+p3mup4*F(1235)+p2mup4*F
     -   (2710))+p4mup4*F(6942))+mup2mup4*(p4mup3*F(6941)+2*(p1mup3*F
     -   (885)+p3mup3*F(894)+p2mup3*F(7424)))+mup3mup4*(2*(p1mup2*F(1
     -   325)+p2mup2*F(1623)+p3mup2*F(2717))+p4mup2*F(7431))
       Fa10(2717) = mup3mup4*p1mup2*F(6932)+mup2mup3*p1mup4*F(6932)+2*
     -   (mup3mup4*p3mup2*F(6937)+mup2mup3*p3mup4*F(6937)+mup3mup4*p4
     -   mup2*F(6943)+mup2mup3*p4mup4*F(6943)+(mup3mup4*p2mup2+mup2mu
     -   p3*p2mup4)*F(7425))+mup2mup4*(p1mup3*F(6931)+2*(p3mup3*F(693
     -   8)+p4mup3*F(6944)+p2mup3*F(7426)))
       Fa10(2718) = mup2mup3*(2*(p1mup4*F(327)+p2mup4*F(2414)+p4mup4*F
     -   (6880))+p3mup4*F(6934))+mup2mup4*(p3mup3*F(6933)+2*(p4mup3*F
     -   (1016)+p1mup3*F(6927)+p2mup3*F(7399)))+mup3mup4*(2*(p1mup2*F
     -   (1325)+p2mup2*F(1623)+p4mup2*F(2718))+p3mup2*F(7432))
       Fa10(2719) = mup2mup4*p2mup3*F(1093)+mup2mup3*p2mup4*F(1094)+mu
     -   p3mup4*(2*(p3mup2*F(6937)+p4mup2*F(6943)+p1mup2*F(7433))+p2m
     -   up2*F(7434))
       Fa10(2720) = p2mup1*(p3mup4*F(6925)+2*(p4mup4*F(1016)+p1mup4*F(
     -   6927)+p2mup4*F(7399)))-p2mup4*p3mup1*F(7418)
       Fa10(2721) = -(mup2mup3*(p4mup4*F(6919)+2*(p1mup4*F(6917)+p3mup
     -   4*F(6918)+p2mup4*F(7416))))+mup3mup4*(2*(p1mup2*F(1325)+p2mu
     -   p2*F(1623)+p3mup2*F(2717))-p4mup2*F(7419))+mup2mup4*(p4mup3*
     -   F(6920)+2*(p1mup3*F(885)+p3mup3*F(894)+p2mup3*F(7424)))
       Fa10(2722) = mup2mup3*(2*(p1mup4*F(327)+p2mup4*F(2414)+p4mup4*F
     -   (6880))-p3mup4*F(6925))-mup2mup4*(p3mup3*F(6926)+2*(p1mup3*F
     -   (6924)+p4mup3*F(6928)+p2mup3*F(7420)))+mup3mup4*(2*(p1mup2*F
     -   (1325)+p2mup2*F(1623)+p4mup2*F(2718))+p3mup2*F(7421))
       Fa10(2723) = p3mup1*(p3mup4*F(6982)+p4mup4*F(6984)+4*(p1mup4*F(
     -   6981)+p2mup4*F(7400)))
       Fa10(2724) = p3mup1*(p1mup4*F(6874)+2*(p3mup4*F(6871)+p4mup4*F(
     -   6872)+p2mup4*F(7397)))
       Fa10(2725) = p2mup1*p3mup4*F(6877)+4*m1sq*p2mup1*(p1mup4*F(327)
     -   +p2mup4*F(2414)+p4mup4*F(6880))+2*p2mup4*p3mup1*F(7398)
       Fa10(2726) = 2*p2mup1*(p1mup4*F(6874)+2*(p3mup4*F(6871)+p4mup4*
     -   F(6872)+p2mup4*F(7397)))
       Fa10(2727) = p2mup1*(p3mup4*F(6865)+p4mup4*F(6867)+4*(p1mup4*F(
     -   6863)+p2mup4*F(7396)))
       Fa10(2728) = mup2mup3*(p3mup4*F(6982)+p4mup4*F(6984)+4*(p1mup4*
     -   F(6981)+p2mup4*F(7400)))+mup2mup4*(-2*p4mup3*F(6868)+4*(p1mu
     -   p3*F(6862)+p3mup3*F(6983)+p2mup3*F(7402)))+2*mup3mup4*(p4mup
     -   2*F(7405)-2*(p1mup2*F(7401)-p2mup2*F(7435)-p3mup2*F(7436)))
       Fa10(2729) = mup2mup3*p3mup4*F(6985)+4*m1sq*(mup3mup4*(p2mup2*F
     -   (1335)+p1mup2*F(1646)+p4mup2*F(6355))+mup2mup3*(p4mup4*F(101
     -   6)+p1mup4*F(6927)+p2mup4*F(7399)))+mup2mup4*(p3mup3*F(6986)+
     -   4*(p1mup3*F(6876)+p4mup3*F(6879)+p2mup3*F(7410)))-2*mup3mup4
     -   *p3mup2*F(7411)
       Fa10(2730) = p3mup1*(p1mup2*F(7401)+p2mup2*F(7403)+p3mup2*F(740
     -   4))
       Fa10(2731) = p1mup1*(p1mup2*F(7401)+p2mup2*F(7403)+p3mup2*F(740
     -   4))
       Fa10(2732) = p1mup2*p3mup1*F(6874)+2*(m1sq*p1mup1*(p1mup2*F(132
     -   5)+p2mup2*F(1623))+p2mup2*p3mup1*F(7406)+p3mup1*p3mup2*F(740
     -   8))+p1mup1*p3mup2*F(7411)
       Fa10(2733) = p3mup1*(p1mup3*F(6862)-p3mup3*F(6866)+p2mup3*F(740
     -   2))
       Fa10(2734) = p1mup1*(p1mup3*F(6862)-p3mup3*F(6866)+p2mup3*F(740
     -   2))
       Fa10(2735) = -2*p3mup1*(p1mup3*F(6874)+2*(p3mup3*F(6870)+p2mup3
     -   *F(7407)))+p1mup1*(p3mup3*F(6878)-4*(p1mup3*F(6876)+p2mup3*F
     -   (7410)))
       Fa10(2736) = p4mup1*(p3mup4*F(6877)+4*m1sq*(p1mup4*F(327)+p2mup
     -   4*F(2414)+p4mup4*F(6880)))
       Fa10(2737) = p1mup1*(p3mup4*F(6865)+p4mup4*F(6867)+4*(p1mup4*F(
     -   6863)+p2mup4*F(7396)))-2*p4mup1*(p1mup4*F(6874)+2*(p3mup4*F(
     -   6871)+p4mup4*F(6872)+p2mup4*F(7397)))
       Fa10(2738) = p1mup1*(p3mup4*F(6877)+4*m1sq*(p1mup4*F(327)+p2mup
     -   4*F(2414)+p4mup4*F(6880)))
       Fa10(2739) = p2mup1*p3mup3*F(6986)+2*p2mup3*p3mup1*F(7398)+4*(p
     -   1mup3*p2mup1*F(6876)+p2mup1*p2mup3*F(7410))
       Fa10(2740) = -2*(p1mup1*p2mup3*F(7398)-p2mup1*(p1mup3*F(6874)+2
     -   *(p3mup3*F(6870)+p2mup3*F(7407))))
       Fa10(2741) = 4*m1sq*p2mup1*(p2mup2*F(1335)+p1mup2*F(1646))-2*(p
     -   2mup1*p3mup2*F(7411)-p3mup1*(p2mup2*F(7413)+2*(p1mup2*F(7412
     -   )+p3mup2*F(7414))))
       Fa10(2742) = p1mup2*p2mup1*F(6874)+p1mup1*p2mup2*F(7413)+2*(p2m
     -   up1*p2mup2*F(7406)+p2mup1*p3mup2*F(7408)+p1mup1*p1mup2*F(741
     -   2)+p1mup1*p3mup2*F(7414))
       Fa10(2743) = mup1mup3*(p3mup4*F(6982)+p4mup4*F(6984)+4*(p1mup4*
     -   F(6981)+p2mup4*F(7400)))+4*mup1mup4*(p1mup3*F(6862)+p3mup3*F
     -   (6983)+p2mup3*F(7402))
       Fa10(2744) = mup1mup3*p3mup4*F(6985)-2*mup3mup4*p3mup1*F(7398)+
     -   4*m1sq*mup1mup3*(p4mup4*F(1016)+p1mup4*F(6927)+p2mup4*F(7399
     -   ))+mup1mup4*(p3mup3*F(6986)+4*(p1mup3*F(6876)+p2mup3*F(7410)
     -   ))
       Fa10(2745) = p2mup1*(-(p1mup2*F(7401))+p2mup2*F(7435)+p3mup2*F(
     -   7436))
       Fa10(2746) = mup1mup2*(p3mup4*F(6865)+p4mup4*F(6867)+4*(p1mup4*
     -   F(6863)+p2mup4*F(7396)))+4*mup1mup4*(-(p1mup2*F(7401))+p2mup
     -   2*F(7435)+p3mup2*F(7436))
       Fa10(2747) = (-(mup2mup4*p2mup1)+mup1mup2*p2mup4)*F(7398)+mup1m
     -   up4*(p2mup2*F(7413)+2*(p1mup2*F(7412)+p3mup2*F(7414)))
       Fa10(2748) = 4*m1sq*mup1mup3*(p1mup2*F(1325)+p2mup2*F(1623))+mu
     -   p1mup2*(p3mup3*F(6878)-4*(p1mup3*F(6876)+p2mup3*F(7410)))+2*
     -   mup1mup3*p3mup2*F(7411)
       Fa10(2749) = mup1mup2*p2mup3*F(7398)+mup1mup3*(p2mup2*F(7413)+2
     -   *(p1mup2*F(7412)+p3mup2*F(7414)))
       Fa10(2750) = 2*(mup1mup3*(p1mup4*F(6874)+2*(p3mup4*F(6871)+p4mu
     -   p4*F(6872)+p2mup4*F(7397)))-mup3mup4*p1mup1*F(7398)-mup1mup4
     -   *(p1mup3*F(6874)+2*(p3mup3*F(6870)+p2mup3*F(7407))))
       Fa10(2751) = 2*(-(mup1mup2*(p1mup4*F(6874)+2*(p3mup4*F(6871)+p4
     -   mup4*F(6872)+p2mup4*F(7397))))+mup2mup4*p1mup1*F(7398)+mup1m
     -   up4*(p1mup2*F(6874)+2*(p2mup2*F(7406)+p3mup2*F(7408))))
       Fa10(2752) = p2mup1*(p1mup3*F(6862)+p3mup3*F(6983)+p2mup3*F(740
     -   2))
       Fa10(2753) = mup1mup2*p3mup4*F(6877)+4*m1sq*(mup1mup2*p1mup4*F(
     -   327)+mup1mup4*p2mup2*F(1335)+mup1mup4*p1mup2*F(1646)+mup1mup
     -   2*p2mup4*F(2414)+mup1mup2*p4mup4*F(6880))+2*(mup2mup4*p3mup1
     -   *F(7398)-mup1mup4*p3mup2*F(7411))
       Fa10(2754) = -(mup1mup2*(p1mup3*F(6862)-p3mup3*F(6866)+p2mup3*F
     -   (7402)))+mup1mup3*(p1mup2*F(7401)+p2mup2*F(7403)+p3mup2*F(74
     -   04))
       Fa10(2755) = 2*(-(mup2mup3*p1mup1*F(7398))+mup1mup2*(p1mup3*F(6
     -   874)+2*(p3mup3*F(6870)+p2mup3*F(7407)))-mup1mup3*(p1mup2*F(6
     -   874)+2*(p2mup2*F(7406)+p3mup2*F(7408))))
       Fa10(2756) = p3mup1*(p1mup2*F(6921)+2*(p3mup2*F(6922)+p2mup2*F(
     -   7417)))+p1mup1*(-2*(p2mup2*F(1335)+p1mup2*F(1646))+p3mup2*F(
     -   7421))
       Fa10(2757) = p3mup1*(p1mup3*F(271)+p3mup3*F(1235)+p2mup3*F(2710
     -   ))
       Fa10(2758) = p1mup1*(p1mup3*F(271)+p3mup3*F(1235)+p2mup3*F(2710
     -   ))
       Fa10(2759) = p1mup3*p3mup1*F(6921)+p1mup1*p3mup3*F(6926)+2*(p3m
     -   up1*p3mup3*F(6922)+p1mup1*p1mup3*F(6924)+p2mup3*p3mup1*F(741
     -   7)+p1mup1*p2mup3*F(7420))
       Fa10(2760) = p4mup1*(p3mup4*F(6925)+2*(p4mup4*F(1016)+p1mup4*F(
     -   6927)+p2mup4*F(7399)))
       Fa10(2761) = p1mup1*(p4mup4*F(6919)+2*(p1mup4*F(6917)+p3mup4*F(
     -   6918)+p2mup4*F(7416)))+p4mup1*(p1mup4*F(6921)+2*(p3mup4*F(69
     -   22)+p4mup4*F(6923)+p2mup4*F(7417)))
       Fa10(2762) = p1mup1*(p3mup4*F(6925)+2*(p4mup4*F(1016)+p1mup4*F(
     -   6927)+p2mup4*F(7399)))
       Fa10(2763) = p2mup1*p3mup3*F(6926)+p2mup3*p3mup1*F(7418)+2*(p1m
     -   up3*p2mup1*F(6924)+p2mup1*p2mup3*F(7420))
       Fa10(2764) = p1mup3*p2mup1*F(6921)+2*p2mup1*(p3mup3*F(6922)+p2m
     -   up3*F(7417))-p1mup1*p2mup3*F(7418)
       Fa10(2765) = p2mup1*p3mup2*F(7421)+2*(p1mup2*p2mup1*F(1325)+p2m
     -   up1*p2mup2*F(1623)-p3mup1*p3mup2*F(6922)-p1mup2*p3mup1*F(742
     -   2))-p2mup2*p3mup1*F(7423)
       Fa10(2766) = p1mup2*p2mup1*F(6921)+2*(p1mup1*p3mup2*F(6922)+p2m
     -   up1*p3mup2*F(6922)+p2mup1*p2mup2*F(7417)+p1mup1*p1mup2*F(742
     -   2))+p1mup1*p2mup2*F(7423)
       Fa10(2767) = -(mup1mup3*(p4mup4*F(6919)+2*(p1mup4*F(6917)+p3mup
     -   4*F(6918)+p2mup4*F(7416))))+2*mup1mup4*(p1mup3*F(885)+p3mup3
     -   *F(894)+p2mup3*F(7424))
       Fa10(2768) = mup1mup3*(2*(p1mup4*F(327)+p2mup4*F(2414)+p4mup4*F
     -   (6880))-p3mup4*F(6925))+mup3mup4*p3mup1*F(7418)-mup1mup4*(p3
     -   mup3*F(6926)+2*(p1mup3*F(6924)+p2mup3*F(7420)))
       Fa10(2769) = 2*mup1mup4*(p1mup2*F(1325)+p2mup2*F(1623)+p3mup2*F
     -   (2717))+mup1mup2*(p4mup4*F(6919)+2*(p1mup4*F(6917)+p3mup4*F(
     -   6918)+p2mup4*F(7416)))
       Fa10(2770) = (-(mup2mup4*p2mup1)+mup1mup2*p2mup4)*F(7418)+mup1m
     -   up4*(2*(p3mup2*F(6922)+p1mup2*F(7422))+p2mup2*F(7423))
       Fa10(2771) = mup1mup2*(p3mup3*F(6926)+2*(p1mup3*F(6924)+p2mup3*
     -   F(7420)))+mup1mup3*(2*(p2mup2*F(1335)+p1mup2*F(1646))-p3mup2
     -   *F(7421))
       Fa10(2772) = mup1mup2*p2mup3*F(7418)+mup1mup3*(2*(p3mup2*F(6922
     -   )+p1mup2*F(7422))+p2mup2*F(7423))
       Fa10(2773) = mup1mup4*(p1mup3*F(6921)+2*(p3mup3*F(6922)+p2mup3*
     -   F(7417)))-mup1mup3*(p1mup4*F(6921)+2*(p3mup4*F(6922)+p4mup4*
     -   F(6923)+p2mup4*F(7417)))+mup3mup4*p1mup1*F(7418)
       Fa10(2774) = -(mup1mup4*(p1mup2*F(6921)+2*(p3mup2*F(6922)+p2mup
     -   2*F(7417))))+mup1mup2*(p1mup4*F(6921)+2*(p3mup4*F(6922)+p4mu
     -   p4*F(6923)+p2mup4*F(7417)))-mup2mup4*p1mup1*F(7418)
       Fa10(2775) = p2mup1*(p1mup3*F(885)+p3mup3*F(894)+p2mup3*F(7424)
     -   )
       Fa10(2776) = mup1mup2*(p3mup4*F(6925)+2*(p4mup4*F(1016)+p1mup4*
     -   F(6927)+p2mup4*F(7399)))-mup2mup4*p3mup1*F(7418)+mup1mup4*(2
     -   *(p1mup2*F(1325)+p2mup2*F(1623))+p3mup2*F(7421))
       Fa10(2777) = mup1mup2*(p1mup3*F(271)+p3mup3*F(1235)+p2mup3*F(27
     -   10))+mup1mup3*(p2mup2*F(1335)+p1mup2*F(1646)+p3mup2*F(5681))
       Fa10(2778) = mup1mup3*p1mup2*F(6921)-mup1mup2*p1mup3*F(6921)+2*
     -   (mup1mup3*p3mup2*F(6922)-mup1mup2*p3mup3*F(6922)+(mup1mup3*p
     -   2mup2-mup1mup2*p2mup3)*F(7417))+mup2mup3*p1mup1*F(7418)
       Fa10(2779) = p3mup1*(p1mup2*F(6932)+2*(p3mup2*F(6937)+p2mup2*F(
     -   7425)))+p1mup1*(2*(p2mup2*F(1335)+p1mup2*F(1646))+p3mup2*F(7
     -   428))
       Fa10(2780) = p1mup1*(2*(p1mup3*F(327)+p2mup3*F(2414))+p3mup3*F(
     -   6934))+p3mup1*(p1mup3*F(6931)+2*(p3mup3*F(6938)+p2mup3*F(742
     -   6)))
       Fa10(2781) = p4mup1*(p3mup4*F(6933)+2*(p4mup4*F(1016)+p1mup4*F(
     -   6927)+p2mup4*F(7399)))
       Fa10(2782) = p1mup1*(p4mup4*F(6941)+2*(p1mup4*F(885)+p3mup4*F(8
     -   94)+p2mup4*F(7424)))+p4mup1*(p1mup4*F(6931)+2*(p3mup4*F(6938
     -   )+p4mup4*F(6944)+p2mup4*F(7426)))
       Fa10(2783) = p1mup1*(p3mup4*F(6933)+2*(p4mup4*F(1016)+p1mup4*F(
     -   6927)+p2mup4*F(7399)))
       Fa10(2784) = p2mup3*p3mup1*F(1093)+p2mup1*p3mup3*F(6933)+2*(p1m
     -   up3*p2mup1*F(6927)+p2mup1*p2mup3*F(7399))
       Fa10(2785) = p1mup1*p2mup3*F(1094)+p2mup1*(p1mup3*F(6932)+2*(p3
     -   mup3*F(6937)+p2mup3*F(7425)))
       Fa10(2786) = p2mup1*p3mup2*F(7432)+2*(p2mup1*p2mup2*F(1623)+p3m
     -   up1*p3mup2*F(6937)+p1mup2*(p2mup1*F(1325)+p3mup1*F(7433)))+p
     -   2mup2*p3mup1*F(7434)
       Fa10(2787) = p1mup2*p2mup1*F(6931)+2*(p1mup1*p3mup2*F(6938)+p2m
     -   up1*p3mup2*F(6938)+p2mup1*p2mup2*F(7426)+p1mup1*p1mup2*F(742
     -   9))+p1mup1*p2mup2*F(7430)
       Fa10(2788) = mup1mup3*(2*(p1mup4*F(271)+p3mup4*F(1235)+p2mup4*F
     -   (2710))+p4mup4*F(6942))+2*mup1mup4*(p1mup3*F(885)+p3mup3*F(8
     -   94)+p2mup3*F(7424))
       Fa10(2789) = mup3mup4*p3mup1*F(1094)+mup1mup3*(2*(p1mup4*F(327)
     -   +p2mup4*F(2414)+p4mup4*F(6880))+p3mup4*F(6934))+mup1mup4*(p3
     -   mup3*F(6933)+2*(p1mup3*F(6927)+p2mup3*F(7399)))
       Fa10(2790) = 2*mup1mup4*(p1mup2*F(1325)+p2mup2*F(1623)+p3mup2*F
     -   (2717))+mup1mup2*(p4mup4*F(6941)+2*(p1mup4*F(885)+p3mup4*F(8
     -   94)+p2mup4*F(7424)))
       Fa10(2791) = mup1mup2*p2mup4*F(1093)+mup2mup4*p2mup1*F(1094)+mu
     -   p1mup4*(2*(p3mup2*F(6937)+p1mup2*F(7433))+p2mup2*F(7434))
       Fa10(2792) = mup1mup2*(2*(p1mup3*F(327)+p2mup3*F(2414))+p3mup3*
     -   F(6934))+mup1mup3*(2*(p2mup2*F(1335)+p1mup2*F(1646))+p3mup2*
     -   F(7428))
       Fa10(2793) = mup1mup2*p2mup3*F(1094)+mup1mup3*(2*(p3mup2*F(6938
     -   )+p1mup2*F(7429))+p2mup2*F(7430))
       Fa10(2794) = mup3mup4*p1mup1*F(1094)+mup1mup3*(p1mup4*F(6932)+2
     -   *(p3mup4*F(6937)+p4mup4*F(6943)+p2mup4*F(7425)))+mup1mup4*(p
     -   1mup3*F(6931)+2*(p3mup3*F(6938)+p2mup3*F(7426)))
       Fa10(2795) = mup2mup4*p1mup1*F(1093)+mup1mup4*(p1mup2*F(6932)+2
     -   *(p3mup2*F(6937)+p2mup2*F(7425)))+mup1mup2*(p1mup4*F(6931)+2
     -   *(p3mup4*F(6938)+p4mup4*F(6944)+p2mup4*F(7426)))
       Fa10(2796) = (mup3mup4*p2mup1+mup1mup4*p2mup3)*F(1093)+mup1mup3
     -   *p2mup4*F(1094)
       Fa10(2797) = mup2mup4*p3mup1*F(1093)+mup1mup2*(p3mup4*F(6933)+2
     -   *(p4mup4*F(1016)+p1mup4*F(6927)+p2mup4*F(7399)))+mup1mup4*(2
     -   *(p1mup2*F(1325)+p2mup2*F(1623))+p3mup2*F(7432))
       Fa10(2798) = mup2mup3*p1mup1*F(1094)+mup1mup2*(p1mup3*F(6932)+2
     -   *(p3mup3*F(6937)+p2mup3*F(7425)))+mup1mup3*(p1mup2*F(6931)+2
     -   *(p3mup2*F(6938)+p2mup2*F(7426)))
       Fa10(2799) = p3mup1*(p1mup2*F(6931)+2*(p3mup2*F(6938)+p2mup2*F(
     -   7426)))+p1mup1*(2*(p1mup2*F(1325)+p2mup2*F(1623))+p3mup2*F(7
     -   432))
       Fa10(2800) = p3mup1*(p1mup3*F(885)+p3mup3*F(894)+p2mup3*F(7424)
     -   )
       Fa10(2801) = p1mup1*(p1mup3*F(885)+p3mup3*F(894)+p2mup3*F(7424)
     -   )
       Fa10(2802) = p1mup1*(p3mup3*F(6933)+2*(p1mup3*F(6927)+p2mup3*F(
     -   7399)))+p3mup1*(p1mup3*F(6932)+2*(p3mup3*F(6937)+p2mup3*F(74
     -   25)))
       Fa10(2803) = p4mup1*(2*(p1mup4*F(327)+p2mup4*F(2414)+p4mup4*F(6
     -   880))+p3mup4*F(6934))
       Fa10(2804) = p1mup1*(2*(p1mup4*F(271)+p3mup4*F(1235)+p2mup4*F(2
     -   710))+p4mup4*F(6942))+p4mup1*(p1mup4*F(6932)+2*(p3mup4*F(693
     -   7)+p4mup4*F(6943)+p2mup4*F(7425)))
       Fa10(2805) = p1mup1*(2*(p1mup4*F(327)+p2mup4*F(2414)+p4mup4*F(6
     -   880))+p3mup4*F(6934))
       Fa10(2806) = p2mup3*p3mup1*F(1094)+2*(p1mup3*p2mup1*F(327)+p2mu
     -   p1*p2mup3*F(2414))+p2mup1*p3mup3*F(6934)
       Fa10(2807) = p1mup1*p2mup3*F(1093)+p2mup1*(p1mup3*F(6931)+2*(p3
     -   mup3*F(6938)+p2mup3*F(7426)))
       Fa10(2808) = p2mup1*(2*(p2mup2*F(1335)+p1mup2*F(1646))+p3mup2*F
     -   (7428))+p3mup1*(2*(p3mup2*F(6938)+p1mup2*F(7429))+p2mup2*F(7
     -   430))
       Fa10(2809) = p1mup2*p2mup1*F(6932)+2*(p1mup1*p3mup2*F(6937)+p2m
     -   up1*p3mup2*F(6937)+p2mup1*p2mup2*F(7425)+p1mup1*p1mup2*F(743
     -   3))+p1mup1*p2mup2*F(7434)
       Fa10(2810) = 2*mup1mup4*(p1mup3*F(271)+p3mup3*F(1235)+p2mup3*F(
     -   2710))+mup1mup3*(p4mup4*F(6941)+2*(p1mup4*F(885)+p3mup4*F(89
     -   4)+p2mup4*F(7424)))
       Fa10(2811) = mup3mup4*p3mup1*F(1093)+mup1mup4*(2*(p1mup3*F(327)
     -   +p2mup3*F(2414))+p3mup3*F(6934))+mup1mup3*(p3mup4*F(6933)+2*
     -   (p4mup4*F(1016)+p1mup4*F(6927)+p2mup4*F(7399)))
       Fa10(2812) = 2*mup1mup4*(p2mup2*F(1335)+p1mup2*F(1646)+p3mup2*F
     -   (5681))+mup1mup2*(2*(p1mup4*F(271)+p3mup4*F(1235)+p2mup4*F(2
     -   710))+p4mup4*F(6942))
       Fa10(2813) = mup2mup4*p2mup1*F(1093)+mup1mup2*p2mup4*F(1094)+mu
     -   p1mup4*(2*(p3mup2*F(6938)+p1mup2*F(7429))+p2mup2*F(7430))
       Fa10(2814) = mup1mup2*(p3mup3*F(6933)+2*(p1mup3*F(6927)+p2mup3*
     -   F(7399)))+mup1mup3*(2*(p1mup2*F(1325)+p2mup2*F(1623))+p3mup2
     -   *F(7432))
       Fa10(2815) = mup1mup2*p2mup3*F(1093)+mup1mup3*(2*(p3mup2*F(6937
     -   )+p1mup2*F(7433))+p2mup2*F(7434))
       Fa10(2816) = mup3mup4*p1mup1*F(1093)+mup1mup4*(p1mup3*F(6932)+2
     -   *(p3mup3*F(6937)+p2mup3*F(7425)))+mup1mup3*(p1mup4*F(6931)+2
     -   *(p3mup4*F(6938)+p4mup4*F(6944)+p2mup4*F(7426)))
       Fa10(2817) = mup2mup4*p1mup1*F(1094)+mup1mup2*(p1mup4*F(6932)+2
     -   *(p3mup4*F(6937)+p4mup4*F(6943)+p2mup4*F(7425)))+mup1mup4*(p
     -   1mup2*F(6931)+2*(p3mup2*F(6938)+p2mup2*F(7426)))
       Fa10(2818) = p2mup1*(p1mup3*F(271)+p3mup3*F(1235)+p2mup3*F(2710
     -   ))
       Fa10(2819) = mup1mup3*p2mup4*F(1093)+(mup3mup4*p2mup1+mup1mup4*
     -   p2mup3)*F(1094)
       Fa10(2820) = mup2mup4*p3mup1*F(1094)+mup1mup2*(2*(p1mup4*F(327)
     -   +p2mup4*F(2414)+p4mup4*F(6880))+p3mup4*F(6934))+mup1mup4*(2*
     -   (p2mup2*F(1335)+p1mup2*F(1646))+p3mup2*F(7428))
       Fa10(2821) = mup1mup3*(p1mup2*F(1325)+p2mup2*F(1623)+p3mup2*F(2
     -   717))+mup1mup2*(p1mup3*F(885)+p3mup3*F(894)+p2mup3*F(7424))
       Fa10(2822) = mup2mup3*p1mup1*F(1093)+mup1mup3*(p1mup2*F(6932)+2
     -   *(p3mup2*F(6937)+p2mup2*F(7425)))+mup1mup2*(p1mup3*F(6931)+2
     -   *(p3mup3*F(6938)+p2mup3*F(7426)))
       Fa10(2823) = p3mup1*(p1mup2*F(6921)+2*(p3mup2*F(6922)+p2mup2*F(
     -   7417)))+p1mup1*(2*(p1mup2*F(1325)+p2mup2*F(1623))+p3mup2*F(7
     -   421))
       Fa10(2824) = p4mup1*(2*(p1mup4*F(327)+p2mup4*F(2414)+p4mup4*F(6
     -   880))-p3mup4*F(6925))
       Fa10(2825) = p1mup1*(2*(p1mup4*F(327)+p2mup4*F(2414)+p4mup4*F(6
     -   880))-p3mup4*F(6925))
       Fa10(2826) = p2mup1*(2*(p2mup2*F(1335)+p1mup2*F(1646))-p3mup2*F
     -   (7421))+p3mup1*(2*(p3mup2*F(6922)+p1mup2*F(7422))+p2mup2*F(7
     -   423))
       Fa10(2827) = 2*mup1mup4*(p1mup3*F(271)+p3mup3*F(1235)+p2mup3*F(
     -   2710))+mup1mup3*(p4mup4*F(6919)+2*(p1mup4*F(6917)+p3mup4*F(6
     -   918)+p2mup4*F(7416)))
       Fa10(2828) = mup1mup3*(p3mup4*F(6925)+2*(p4mup4*F(1016)+p1mup4*
     -   F(6927)+p2mup4*F(7399)))-mup3mup4*p3mup1*F(7418)+mup1mup4*(p
     -   3mup3*F(6926)+2*(p1mup3*F(6924)+p2mup3*F(7420)))
       Fa10(2829) = 2*mup1mup4*(p2mup2*F(1335)+p1mup2*F(1646)+p3mup2*F
     -   (5681))-mup1mup2*(p4mup4*F(6919)+2*(p1mup4*F(6917)+p3mup4*F(
     -   6918)+p2mup4*F(7416)))
       Fa10(2830) = -(mup1mup2*(p3mup3*F(6926)+2*(p1mup3*F(6924)+p2mup
     -   3*F(7420))))+mup1mup3*(2*(p1mup2*F(1325)+p2mup2*F(1623))+p3m
     -   up2*F(7421))
       Fa10(2831) = -(mup1mup4*(p1mup3*F(6921)+2*(p3mup3*F(6922)+p2mup
     -   3*F(7417))))+mup1mup3*(p1mup4*F(6921)+2*(p3mup4*F(6922)+p4mu
     -   p4*F(6923)+p2mup4*F(7417)))-mup3mup4*p1mup1*F(7418)
       Fa10(2832) = mup1mup4*(p1mup2*F(6921)+2*(p3mup2*F(6922)+p2mup2*
     -   F(7417)))-mup1mup2*(p1mup4*F(6921)+2*(p3mup4*F(6922)+p4mup4*
     -   F(6923)+p2mup4*F(7417)))+mup2mup4*p1mup1*F(7418)
       Fa10(2833) = mup1mup2*(2*(p1mup4*F(327)+p2mup4*F(2414)+p4mup4*F
     -   (6880))-p3mup4*F(6925))+mup2mup4*p3mup1*F(7418)+mup1mup4*(2*
     -   (p2mup2*F(1335)+p1mup2*F(1646))-p3mup2*F(7421))
       Fa10(2834) = -(mup1mup3*(p1mup2*F(6921)+2*(p3mup2*F(6922)+p2mup
     -   2*F(7417))))+mup1mup2*(p1mup3*F(6921)+2*(p3mup3*F(6922)+p2mu
     -   p3*F(7417)))-mup2mup3*p1mup1*F(7418)
       Fa10(2835) = p3mup1*(p1mup2*F(7401)-p2mup2*F(7435)-p3mup2*F(743
     -   6))
       Fa10(2836) = p1mup1*(-(p1mup2*F(7401))+p2mup2*F(7435)+p3mup2*F(
     -   7436))
       Fa10(2837) = 4*m1sq*p1mup1*(p2mup2*F(1335)+p1mup2*F(1646))-2*(p
     -   1mup2*p3mup1*F(6874)+2*(p2mup2*p3mup1*F(7406)+p3mup1*p3mup2*
     -   F(7408))+p1mup1*p3mup2*F(7411))
       Fa10(2838) = p3mup1*(p1mup3*F(6862)+p3mup3*F(6983)+p2mup3*F(740
     -   2))
       Fa10(2839) = p1mup1*(p1mup3*F(6862)+p3mup3*F(6983)+p2mup3*F(740
     -   2))
       Fa10(2840) = 2*p3mup1*(p1mup3*F(6874)+2*(p3mup3*F(6870)+p2mup3*
     -   F(7407)))+p1mup1*(p3mup3*F(6986)+4*(p1mup3*F(6876)+p2mup3*F(
     -   7410)))
       Fa10(2841) = p4mup1*(p3mup4*F(6985)+4*m1sq*(p4mup4*F(1016)+p1mu
     -   p4*F(6927)+p2mup4*F(7399)))
       Fa10(2842) = 2*p1mup4*p4mup1*F(6874)+p1mup1*p3mup4*F(6982)+p1mu
     -   p1*p4mup4*F(6984)+4*(p3mup4*p4mup1*F(6871)+p4mup1*p4mup4*F(6
     -   872)+p1mup1*p1mup4*F(6981)+p2mup4*p4mup1*F(7397)+p1mup1*p2mu
     -   p4*F(7400))
       Fa10(2843) = p1mup1*(p3mup4*F(6985)+4*m1sq*(p4mup4*F(1016)+p1mu
     -   p4*F(6927)+p2mup4*F(7399)))
       Fa10(2844) = p2mup1*p3mup3*F(6878)-2*p2mup3*p3mup1*F(7398)-4*(p
     -   1mup3*p2mup1*F(6876)+p2mup1*p2mup3*F(7410))
       Fa10(2845) = -(p1mup1*p2mup3*F(7398))+p2mup1*(p1mup3*F(6874)+2*
     -   (p3mup3*F(6870)+p2mup3*F(7407)))
       Fa10(2846) = 2*(2*m1sq*p2mup1*(p1mup2*F(1325)+p2mup2*F(1623))+p
     -   2mup1*p3mup2*F(7411)-p3mup1*(p2mup2*F(7413)+2*(p1mup2*F(7412
     -   )+p3mup2*F(7414))))
       Fa10(2847) = mup1mup3*(p3mup4*F(6865)+p4mup4*F(6867)+4*(p1mup4*
     -   F(6863)+p2mup4*F(7396)))-4*mup1mup4*(p1mup3*F(6862)-p3mup3*F
     -   (6866)+p2mup3*F(7402))
       Fa10(2848) = mup1mup3*p3mup4*F(6877)+4*m1sq*mup1mup3*(p1mup4*F(
     -   327)+p2mup4*F(2414)+p4mup4*F(6880))+2*mup3mup4*p3mup1*F(7398
     -   )+mup1mup4*(p3mup3*F(6878)-4*(p1mup3*F(6876)+p2mup3*F(7410))
     -   )
       Fa10(2849) = p2mup1*(p1mup2*F(7401)+p2mup2*F(7403)+p3mup2*F(740
     -   4))
       Fa10(2850) = mup1mup2*(p3mup4*F(6982)+p4mup4*F(6984)+4*(p1mup4*
     -   F(6981)+p2mup4*F(7400)))+4*mup1mup4*(p1mup2*F(7401)+p2mup2*F
     -   (7403)+p3mup2*F(7404))
       Fa10(2851) = 4*m1sq*mup1mup3*(p2mup2*F(1335)+p1mup2*F(1646))+mu
     -   p1mup2*(p3mup3*F(6986)+4*(p1mup3*F(6876)+p2mup3*F(7410)))-2*
     -   mup1mup3*p3mup2*F(7411)
       Fa10(2852) = -(mup1mup3*(p1mup4*F(6874)+2*(p3mup4*F(6871)+p4mup
     -   4*F(6872)+p2mup4*F(7397))))+mup3mup4*p1mup1*F(7398)+mup1mup4
     -   *(p1mup3*F(6874)+2*(p3mup3*F(6870)+p2mup3*F(7407)))
       Fa10(2853) = 2*(mup1mup2*(p1mup4*F(6874)+2*(p3mup4*F(6871)+p4mu
     -   p4*F(6872)+p2mup4*F(7397)))-mup2mup4*p1mup1*F(7398)-mup1mup4
     -   *(p1mup2*F(6874)+2*(p2mup2*F(7406)+p3mup2*F(7408))))
       Fa10(2854) = p2mup1*(p1mup3*F(6862)-p3mup3*F(6866)+p2mup3*F(740
     -   2))
       Fa10(2855) = mup1mup2*p3mup4*F(6985)+4*m1sq*(mup1mup2*p4mup4*F(
     -   1016)+mup1mup4*p1mup2*F(1325)+mup1mup4*p2mup2*F(1623)+mup1mu
     -   p2*p1mup4*F(6927)+mup1mup2*p2mup4*F(7399))-2*(mup2mup4*p3mup
     -   1*F(7398)-mup1mup4*p3mup2*F(7411))
       Fa10(2856) = mup1mup2*(p1mup3*F(6862)+p3mup3*F(6983)+p2mup3*F(7
     -   402))+mup1mup3*(-(p1mup2*F(7401))+p2mup2*F(7435)+p3mup2*F(74
     -   36))
       Fa10(2857) = 2*(mup2mup3*p1mup1*F(7398)-mup1mup2*(p1mup3*F(6874
     -   )+2*(p3mup3*F(6870)+p2mup3*F(7407)))+mup1mup3*(p1mup2*F(6874
     -   )+2*(p2mup2*F(7406)+p3mup2*F(7408))))
       Fa10(2858) = 4*m1sq*p4mup1*(p1mup2*F(1325)+p2mup2*F(1623))-2*(p
     -   3mup1*p4mup2*F(7405)-p3mup2*p4mup1*F(7411))
       Fa10(2859) = p4mup1*(p1mup2*F(6874)+2*(p2mup2*F(7406)+p3mup2*F(
     -   7408)))
       Fa10(2860) = 2*p3mup1*p4mup3*F(6868)+p4mup1*(p3mup3*F(6878)-4*(
     -   p1mup3*F(6876)+p4mup3*F(6879)+p2mup3*F(7410)))
       Fa10(2861) = p1mup1*p4mup3*F(6868)+p4mup1*(p1mup3*F(6874)+2*(p3
     -   mup3*F(6870)+p4mup3*F(6873)+p2mup3*F(7407)))
       Fa10(2862) = p4mup3*(p3mup1*F(6873)+p1mup1*F(6879))
       Fa10(2863) = p4mup1*(p2mup2*F(7413)+2*(p1mup2*F(7412)+p3mup2*F(
     -   7414)))
       Fa10(2864) = 2*mup1mup2*p4mup3*F(6879)+mup2mup3*p3mup1*F(7398)
       Fa10(2865) = p2mup1*p4mup3*F(6868)+p2mup3*p4mup1*F(7398)
       Fa10(2866) = 2*(p2mup2*p4mup1*F(1335)+p1mup2*p4mup1*F(1646))+p3
     -   mup1*p4mup2*F(7419)-p3mup2*p4mup1*F(7421)
       Fa10(2867) = p4mup1*(p1mup2*F(6921)+2*(p3mup2*F(6922)+p2mup2*F(
     -   7417)))
       Fa10(2868) = -(p3mup1*p4mup3*F(6920))+p4mup1*(p3mup3*F(6926)+2*
     -   (p1mup3*F(6924)+p4mup3*F(6928)+p2mup3*F(7420)))
       Fa10(2869) = p1mup1*p4mup3*F(6920)+p4mup1*(p1mup3*F(6921)+2*(p3
     -   mup3*F(6922)+p4mup3*F(6923)+p2mup3*F(7417)))
       Fa10(2870) = p4mup3*(p3mup1*F(6923)+p1mup1*F(6928))
       Fa10(2871) = p4mup1*(2*(p3mup2*F(6922)+p1mup2*F(7422))+p2mup2*F
     -   (7423))
       Fa10(2872) = 2*mup1mup2*p4mup3*F(6928)+mup2mup3*p3mup1*F(7418)
       Fa10(2873) = p2mup1*p4mup3*F(6920)+p2mup3*p4mup1*F(7418)
       Fa10(2874) = 2*(p2mup2*p4mup1*F(1335)+p1mup2*p4mup1*F(1646))+p3
     -   mup1*p4mup2*F(7427)+p3mup2*p4mup1*F(7428)
       Fa10(2875) = p4mup1*(p1mup2*F(6931)+2*(p3mup2*F(6938)+p2mup2*F(
     -   7426)))
       Fa10(2876) = p4mup1*(2*(p1mup3*F(327)+p2mup3*F(2414)+p4mup3*F(6
     -   880))+p3mup3*F(6934))+p3mup1*p4mup3*F(6942)
       Fa10(2877) = p1mup3*p4mup1*F(6932)+p1mup1*p4mup3*F(6942)+2*p4mu
     -   p1*(p3mup3*F(6937)+p4mup3*F(6943)+p2mup3*F(7425))
       Fa10(2878) = p4mup3*(p1mup1*F(6880)+p3mup1*F(6944))
       Fa10(2879) = p4mup1*(2*(p3mup2*F(6938)+p1mup2*F(7429))+p2mup2*F
     -   (7430))
       Fa10(2880) = mup2mup3*p3mup1*F(1094)+2*mup1mup2*p4mup3*F(6880)
       Fa10(2881) = p2mup3*p4mup1*F(1094)+p2mup1*p4mup3*F(6941)
       Fa10(2882) = 2*(p1mup2*p4mup1*F(1325)+p2mup2*p4mup1*F(1623))+p3
     -   mup1*p4mup2*F(7431)+p3mup2*p4mup1*F(7432)
       Fa10(2883) = p4mup1*(p1mup2*F(6932)+2*(p3mup2*F(6937)+p2mup2*F(
     -   7425)))
       Fa10(2884) = p3mup1*p4mup3*F(6941)+p4mup1*(p3mup3*F(6933)+2*(p4
     -   mup3*F(1016)+p1mup3*F(6927)+p2mup3*F(7399)))
       Fa10(2885) = p1mup3*p4mup1*F(6931)+p1mup1*p4mup3*F(6941)+2*p4mu
     -   p1*(p3mup3*F(6938)+p4mup3*F(6944)+p2mup3*F(7426))
       Fa10(2886) = p4mup3*(p1mup1*F(1016)+p3mup1*F(6943))
       Fa10(2887) = p4mup1*(2*(p3mup2*F(6937)+p1mup2*F(7433))+p2mup2*F
     -   (7434))
       Fa10(2888) = 2*mup1mup2*p4mup3*F(1016)+mup2mup3*p3mup1*F(1093)
       Fa10(2889) = p2mup3*p4mup1*F(1093)+p2mup1*p4mup3*F(6942)
       Fa10(2890) = 2*(p1mup2*p4mup1*F(1325)+p2mup2*p4mup1*F(1623))-p3
     -   mup1*p4mup2*F(7419)+p3mup2*p4mup1*F(7421)
       Fa10(2891) = 4*m1sq*p4mup1*(p2mup2*F(1335)+p1mup2*F(1646))+2*(p
     -   3mup1*p4mup2*F(7405)-p3mup2*p4mup1*F(7411))
       Fa10(2892) = 2*p4mup1*(p1mup2*F(6874)+2*(p2mup2*F(7406)+p3mup2*
     -   F(7408)))
       Fa10(2893) = -2*p3mup1*p4mup3*F(6868)+p4mup1*(p3mup3*F(6986)+4*
     -   (p1mup3*F(6876)+p4mup3*F(6879)+p2mup3*F(7410)))
       Fa10(2894) = 2*p4mup1*(p2mup2*F(7413)+2*(p1mup2*F(7412)+p3mup2*
     -   F(7414)))
       Fa10(2895) = 4*mup1mup2*p4mup3*F(6879)+2*mup2mup3*p3mup1*F(7398
     -   )
       Fa10(2896) = p4mup2*(p1mup1*F(7405)+2*p4mup1*F(7409))
       Fa10(2897) = p4mup2*(p2mup1*F(7409)+p1mup1*F(7415))
       Fa10(2898) = 2*p4mup2*(p2mup1*F(7405)-2*p4mup1*F(7415))
       Fa10(2899) = mup2mup4*p4mup1*F(7398)-mup1mup4*p4mup2*F(7405)
       Fa10(2900) = mup2mup3*p4mup1*F(7398)-mup1mup3*p4mup2*F(7405)
       Fa10(2901) = p4mup2*(2*p4mup1*F(6923)+p1mup1*F(7419))
       Fa10(2902) = p4mup2*(2*p4mup1*F(6923)-p2mup1*F(7419))
       Fa10(2903) = mup2mup4*p4mup1*F(7418)-mup1mup4*p4mup2*F(7419)
       Fa10(2904) = -(mup2mup3*p4mup1*F(7418))+mup1mup3*p4mup2*F(7419)
       Fa10(2905) = p4mup2*(2*p4mup1*F(6944)+p1mup1*F(7427))
       Fa10(2906) = p4mup2*(2*p4mup1*F(6944)+p2mup1*F(7431))
       Fa10(2907) = mup2mup4*p4mup1*F(1094)+mup1mup4*p4mup2*F(7431)
       Fa10(2908) = mup2mup3*p4mup1*F(1093)+mup1mup3*p4mup2*F(7427)
       Fa10(2909) = (mup1mup4*mup2mup3+mup1mup2*mup3mup4)*F(1093)+mup1
     -   mup3*mup2mup4*F(1094)
       Fa10(2910) = p4mup2*(2*p4mup1*F(6943)+p1mup1*F(7431))
       Fa10(2911) = p4mup2*(2*p4mup1*F(6943)+p2mup1*F(7427))
       Fa10(2912) = mup2mup4*p4mup1*F(1093)+mup1mup4*p4mup2*F(7427)
       Fa10(2913) = mup2mup3*p4mup1*F(1094)+mup1mup3*p4mup2*F(7431)
       Fa10(2914) = mup1mup3*mup2mup4*F(1093)+(mup1mup4*mup2mup3+mup1m
     -   up2*mup3mup4)*F(1094)
       Fa10(2915) = -2*p4mup1*p4mup2*F(6923)+p2mup1*p4mup2*F(7419)
       Fa10(2916) = -(mup2mup4*p4mup1*F(7418))+mup1mup4*p4mup2*F(7419)
       Fa10(2917) = mup2mup3*p4mup1*F(7418)-mup1mup3*p4mup2*F(7419)
       Fa10(2918) = 2*p4mup2*(p1mup1*F(7405)+2*p4mup1*F(7409))
       Fa10(2919) = -2*p2mup1*p4mup2*F(7405)+4*p4mup1*p4mup2*F(7415)
       Fa10(2920) = -4*Fa(2)+Fa(115)
       Return
       End
