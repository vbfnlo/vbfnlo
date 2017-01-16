c************** Calling the Fa functions*************************
       subroutine NoAbe2gaHexLineFa30(p1mup2,p1mup3,p1mup5,p1mup6,p2mu
     -   p2,p2mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p
     -   4mup3,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup
     -   3,p6mup5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6
     -   ,mup5mup6,Fa30)
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
       Complex*16 Fa30(697:720)
       COMMON/NoAbe2gaHexLineFaFunctions/Fa
       Fa30(697) = mup5mup6*(p2mup2*F(8139)+p3mup2*F(8140)+p4mup2*F(81
     -   41)+p5mup2*F(8142))
       Fa30(698) = p2mup6*p3mup2*p6mup5*F(8178)-4*(p1mup6*p3mup2*p5mup
     -   5*F(8149)+p3mup2*p5mup5*p6mup6*F(8158)-p4mup2*(-(p5mup5*(p6m
     -   up6*F(8145)+p1mup6*F(8156)+2*(p2mup6*F(8152)+p5mup6*F(8159))
     -   ))+2*(p2mup5*p6mup6*F(8151)-p1mup5*p6mup6*F(8155)+p6mup5*(p6
     -   mup6*F(8144)-2*(p2mup6*F(8153)+p1mup6*F(8157)+p5mup6*F(8160)
     -   )))+4*(p1mup5*p2mup6*F(8167)+p1mup6*p2mup5*F(8168)+p1mup5*p1
     -   mup6*F(8169)+p1mup5*p5mup6*F(8172)+p2mup5*p2mup6*F(8175)+p2m
     -   up5*p5mup6*F(8179))))-8*(p1mup5*p3mup2*p6mup6*F(8148)+p3mup2
     -   *p6mup5*p6mup6*F(8154)+p5mup2*p5mup5*p5mup6*F(8161)+p2mup2*(
     -   -2*(p1mup5*p1mup6*F(8164)-p1mup5*p5mup6*F(8170))+p2mup6*(p5m
     -   up5*F(8146)+2*(p6mup5*F(8147)+p1mup5*F(8163)+p2mup5*F(8174))
     -   ))+p2mup5*p3mup2*p6mup6*F(8176)-p2mup6*p3mup2*p5mup5*F(8177)
     -   -p3mup2*p5mup5*p5mup6*F(8180))+p3mup2*p5mup6*p6mup5*F(8181)-
     -   16*(p1mup6*p3mup2*p6mup5*F(8150)+p5mup2*p5mup6*p6mup5*F(8162
     -   )-p1mup5*p1mup6*p3mup2*F(8165)+p1mup5*p2mup6*p3mup2*F(8166)+
     -   p1mup5*p3mup2*p5mup6*F(8171)-p1mup5*p5mup2*p5mup6*F(8173)-p2
     -   mup5*p5mup2*p5mup6*F(8182)-p1mup6*p2mup5*p3mup2*F(8183)-p2mu
     -   p5*p2mup6*p3mup2*F(8184)-p2mup5*p3mup2*p5mup6*F(8185))
       Fa30(699) = -2*mup5mup6*p6mup2*F(7384)+mup2mup6*(2*p5mup5*F(818
     -   6)+p6mup5*F(8187)+p1mup5*F(8188)+p2mup5*F(8189))+p5mup2*(-4*
     -   p5mup5*(2*p2mup6*F(8192)+p1mup6*F(8199))+p5mup5*p6mup6*F(820
     -   1)+16*(p1mup5*p2mup6*F(8206)+p1mup5*p1mup6*F(8208))-8*(p1mup
     -   5*p6mup6*F(8198)-p6mup5*(p6mup6*F(8191)-2*(p2mup6*F(8193)+p1
     -   mup6*F(8200)))-p2mup5*(p6mup6*F(8194)+2*(p1mup6*F(8207)+p2mu
     -   p6*F(8209)))))+p2mup2*(p5mup5*p6mup6*F(8196)-4*p1mup6*p5mup5
     -   *F(8203)+8*(p5mup5*p5mup6*F(8195)-p1mup5*p6mup6*F(8202)-p6mu
     -   p5*(p6mup6*F(8190)-2*(p5mup6*F(8197)-p1mup6*F(8204)))+p2mup5
     -   *(p6mup6*F(4109)-2*(p1mup6*F(8205)+p5mup6*F(8210)))))
       Fa30(700) = mup2mup6*p4mup5*F(8211)+4*(2*p2mup2*p4mup5*(p6mup6*
     -   F(8221)-2*(p4mup6*F(8225)+p5mup6*F(8229)))+p1mup2*(p5mup5*(p
     -   6mup6*F(8212)-p1mup6*F(8214)-2*(p2mup6*F(8218)+p5mup6*F(8222
     -   )))+2*(-(p2mup5*p6mup6*F(8217))+p4mup5*p6mup6*F(8220)+p6mup5
     -   *(p6mup6*F(8213)-2*(p1mup6*F(8215)+p2mup6*F(8219)+p5mup6*F(8
     -   223)))+p1mup5*(-(p6mup6*F(8216))+2*(p5mup6*F(8226)+p2mup6*F(
     -   8230)+p1mup6*F(8232)))-2*(p4mup5*p4mup6*F(8224)-p2mup5*p5mup
     -   6*F(8227)+p4mup5*p5mup6*F(8228)-p2mup5*p2mup6*F(8231)-p1mup6
     -   *p2mup5*F(8233)))))
       Fa30(701) = mup2mup5*p1mup6*F(8234)+mup2mup5*p6mup6*F(8235)+mup
     -   2mup5*p2mup6*F(8236)+mup2mup5*p5mup6*F(8237)+8*(p2mup5*p6mup
     -   2*p6mup6*F(7575)-p4mup2*p4mup5*p6mup6*F(8239)-p4mup5*p5mup2*
     -   p6mup6*F(8240)+p3mup2*p4mup5*p6mup6*F(8245)-p2mup2*(-2*(p2mu
     -   p6*p4mup5*F(8243)-p1mup6*p4mup5*F(8252))+p3mup5*(p6mup6*F(41
     -   11)+2*(p1mup6*F(8257)+p2mup6*F(8262)+p5mup6*F(8265)))))-16*(
     -   p1mup6*p2mup5*p6mup2*F(7501)-p2mup5*p5mup6*p6mup2*F(7532)-p2
     -   mup6*p4mup2*p4mup5*F(8242)-p2mup6*p4mup5*p5mup2*F(8244)-p1mu
     -   p2*p2mup6*p4mup5*F(8247)-p4mup2*p4mup5*p5mup6*F(8248)-p4mup5
     -   *p5mup2*p5mup6*F(8249)+p1mup2*p1mup6*p4mup5*F(8253)-p1mup6*p
     -   4mup2*p4mup5*F(8254)-p1mup6*p4mup5*p5mup2*F(8255)+p1mup6*p3m
     -   up2*p4mup5*F(8258)+p2mup6*p3mup2*p4mup5*F(8261)+p3mup2*p4mup
     -   5*p5mup6*F(8266))+p3mup5*(mup2mup6*F(8238)-8*(-(p1mup2*(-(p6
     -   mup6*F(8241))+2*(p2mup6*F(8246)+p5mup6*F(8251)+p1mup6*F(8256
     -   ))))+p5mup2*(p6mup6*F(5681)-2*(p1mup6*F(8259)+p2mup6*F(8263)
     -   +p5mup6*F(8267)))-p4mup2*(-(p6mup6*F(8250))+2*(p1mup6*F(8260
     -   )+p2mup6*F(8264)+p5mup6*F(8268)))+p3mup2*(p6mup6*F(6075)-2*(
     -   p1mup6*F(8269)+p2mup6*F(8270)+p5mup6*F(8271)))))
       Fa30(702) = p4mup6*(p1mup2*(-(p5mup5*F(8272))+2*(p1mup5*F(8274)
     -   +p2mup5*F(8275)+p3mup5*F(8276)))+p2mup2*(p5mup5*F(8273)-2*(p
     -   1mup5*F(8277)+p2mup5*F(8278)+p3mup5*F(8279))))
       Fa30(703) = p6mup2*(-(p5mup5*p5mup6*F(7508))+p2mup6*(2*(p4mup5*
     -   F(7506)+p2mup5*F(7510)+p3mup5*F(7511))-p5mup5*F(7520))-2*(p1
     -   mup6*(p3mup5*F(7503)+p4mup5*F(7505))-p4mup5*p5mup6*F(7527)-p
     -   3mup5*p5mup6*F(7536))+p3mup5*p6mup6*F(7576)+p4mup5*p6mup6*F(
     -   8280))
       Fa30(704) = 4*p1mup6*p5mup5*p6mup2*F(7514)+mup2mup5*p3mup6*F(82
     -   81)+mup2mup5*p4mup6*F(8282)+p3mup2*p3mup6*p6mup5*F(8309)-8*(
     -   p1mup2*p3mup6*p5mup5*F(8283)-p5mup5*p6mup2*p6mup6*F(8286)-p2
     -   mup2*p3mup6*p5mup5*F(8287)+p4mup2*p4mup6*p5mup5*F(8293)+p4mu
     -   p6*p5mup2*p5mup5*F(8295)-p3mup6*p4mup2*p5mup5*F(8296)-p3mup6
     -   *p5mup2*p5mup5*F(8299)-p1mup5*(p6mup2*(2*(p2mup6*F(7500)+p1m
     -   up6*F(7509)+p5mup6*F(7519))-p6mup6*F(7572))+2*(p1mup2*p3mup6
     -   *F(8289)-p2mup2*p3mup6*F(8301)+p3mup6*p5mup2*F(8302)-p3mup2*
     -   p3mup6*F(8303)+p3mup6*p4mup2*F(8304)+p4mup2*p4mup6*F(8305)+p
     -   4mup6*p5mup2*F(8306)-p3mup2*p4mup6*F(8307)))-p3mup2*p3mup6*p
     -   5mup5*F(8308)-p3mup2*p4mup6*p5mup5*F(8312))+p3mup2*p4mup6*p6
     -   mup5*F(8313)+16*(p2mup5*p3mup6*(p1mup2*F(8246)-p2mup2*F(8262
     -   )+p5mup2*F(8263)+p4mup2*F(8264)+p3mup2*F(8270))-p1mup2*p3mup
     -   6*p6mup5*F(8284)-p1mup2*p4mup6*p6mup5*F(8285)+p2mup2*p3mup6*
     -   p6mup5*F(8288)-p1mup2*p3mup6*p4mup5*F(8290)-p2mup2*p3mup6*p4
     -   mup5*F(8291)+p2mup2*p4mup6*p6mup5*F(8292)-p4mup2*p4mup6*p6mu
     -   p5*F(8294)-p4mup6*p5mup2*p6mup5*F(8297)+p3mup6*p4mup2*p6mup5
     -   *F(8298)+p3mup6*p5mup2*p6mup5*F(8300)+p2mup5*p4mup2*p4mup6*F
     -   (8310)+p2mup5*p4mup6*p5mup2*F(8311)+p2mup5*p3mup2*p4mup6*F(8
     -   314))
       Fa30(705) = p3mup6*p5mup5*p6mup2*F(8317)+4*(p2mup6*p6mup2*p6mup
     -   5*F(4299)+p3mup6*p6mup2*p6mup5*F(4315)+p5mup6*p6mup2*p6mup5*
     -   F(8318))+p4mup6*(p6mup2*(4*p6mup5*F(4314)+2*(p1mup5*F(7504)+
     -   p2mup5*F(7512))-p5mup5*F(7524))+2*p4mup2*(p4mup5*F(8316)+p3m
     -   up5*F(8319)))+2*(p1mup5*p3mup6*p6mup2*F(7502)+p2mup5*p3mup6*
     -   p6mup2*F(7511)-p1mup6*p6mup2*p6mup5*F(7551)-p3mup6*p4mup2*p4
     -   mup5*F(8315)+p3mup5*p3mup6*p5mup2*F(8320))
       Fa30(706) = -(p3mup6*(p4mup5*(p6mup2*F(7507)+p5mup2*F(8322)+p3m
     -   up2*F(8325))-p3mup5*(p6mup2*F(7534)+p1mup2*F(8323)-p2mup2*F(
     -   8324)-p4mup2*F(8327)+p3mup2*F(8329))))+p4mup6*(p4mup5*(p6mup
     -   2*F(7523)+p5mup2*F(8321)-p3mup2*F(8326))+p3mup5*(p6mup2*F(75
     -   13)+p5mup2*F(8328)+p3mup2*F(8330)))
       Fa30(707) = mup5mup6*(p5mup2*F(8332)+p2mup2*F(8333)+2*p3mup2*F(
     -   8334)+p4mup2*F(8335))
       Fa30(708) = p5mup2*p5mup5*p5mup6*F(8340)+2*p2mup2*p2mup6*p5mup5
     -   *F(8342)-4*(p1mup5*p5mup2*p5mup6*F(8341)-p2mup2*p2mup6*p6mup
     -   5*F(8343)+p2mup5*p5mup2*p5mup6*F(8350)-p5mup2*p5mup6*p6mup5*
     -   F(8354)-p1mup5*p2mup2*p5mup6*F(8355)+p1mup5*p2mup2*p2mup6*F(
     -   8358)-p1mup5*p1mup6*p2mup2*F(8359)+p2mup2*p2mup5*p2mup6*F(83
     -   65))+p4mup2*(p5mup5*(p5mup6*F(8337)+p6mup6*F(8338)+2*(p2mup6
     -   *F(8348)-p1mup6*F(8351)))+2*(p6mup5*(p6mup6*F(8339)+2*(p2mup
     -   6*F(8347)-p1mup6*F(8352)+p5mup6*F(8353)))+p1mup5*(p6mup6*F(8
     -   357)-2*(p5mup6*F(8356)-p2mup6*F(8361)-p1mup6*F(8364)))+p2mup
     -   5*(p6mup6*F(8347)-2*(p5mup6*F(8346)-p1mup6*F(8363)-p2mup6*F(
     -   8366)))))+p3mup2*(p5mup5*p6mup6*F(7712)+p1mup6*(-2*p5mup5*F(
     -   8344)-4*(p2mup5*F(6945)-p6mup5*F(8345)-p1mup5*F(8360)))+p5mu
     -   p6*(-(p5mup5*F(8349))+4*(p6mup5*F(6931)+p1mup5*F(8362)+p2mup
     -   5*F(8367)))-2*(p2mup6*p5mup5*F(7689)+p6mup5*p6mup6*F(9618)-p
     -   1mup5*p6mup6*F(9628)-p2mup5*p6mup6*F(9635))+4*(p2mup6*p6mup5
     -   *F(9635)-p1mup5*p2mup6*F(9662)-p2mup5*p2mup6*F(9672)))
       Fa30(709) = mup2mup6*(p1mup5*F(8368)+p5mup5*F(8369)+p6mup5*F(83
     -   70)+p2mup5*F(8371))+mup5mup6*p6mup2*F(8372)-4*(p2mup2*(p5mup
     -   5*p5mup6*F(8375)-p5mup5*p6mup6*F(8376)-4*(p5mup6*p6mup5*F(83
     -   86)-p1mup6*p6mup5*F(8388))+2*(p6mup5*p6mup6*F(8377)+p1mup6*p
     -   5mup5*F(8387)-p1mup5*p6mup6*F(8389))+p2mup5*(-2*p6mup6*F(834
     -   3)+4*(p5mup6*F(8383)+p1mup6*F(8391))))-p5mup2*(-(p1mup6*p5mu
     -   p5*F(8374))+p5mup5*p6mup6*F(8380)+2*(p1mup5*p6mup6*F(8373)+p
     -   6mup5*p6mup6*F(8381)+p2mup6*p5mup5*F(8384)-p2mup5*p6mup6*F(8
     -   385))-4*(p1mup5*p2mup6*F(8378)-p1mup6*p2mup5*F(8379)-p1mup6*
     -   p6mup5*F(8382)+p2mup6*p6mup5*F(8385)+p1mup5*p1mup6*F(8390)-p
     -   2mup5*p2mup6*F(8392))))
       Fa30(710) = mup2mup6*p4mup5*F(8393)+4*(p1mup2*(-2*(-(p6mup5*p6m
     -   up6*F(8395))+p1mup5*p6mup6*F(8396)+2*p1mup6*p6mup5*F(8398)+(
     -   2*p2mup6*p6mup5+p2mup5*p6mup6)*F(8399))+p5mup5*(p6mup6*F(839
     -   4)-2*(p1mup6*F(8397)+p2mup6*F(8400))-p5mup6*F(8401)))-p2mup2
     -   *p4mup5*(p6mup6*F(8405)+4*(p5mup6*F(8403)+p4mup6*F(8409)))-p
     -   1mup2*(p4mup5*(p6mup6*F(8404)+4*(p5mup6*F(8402)+p4mup6*F(840
     -   8)))-4*(-(p5mup6*p6mup5*F(8406))-p2mup5*p5mup6*F(8407)+p1mup
     -   6*p2mup5*F(8410)+p2mup5*p2mup6*F(8411)+p1mup5*p1mup6*F(8412)
     -   +p1mup5*p2mup6*F(8413)+p1mup5*p5mup6*F(8414))))
       Fa30(711) = 2*mup2mup5*p6mup6*F(8416)+mup2mup6*p3mup5*F(8417)+p
     -   3mup2*p4mup5*p6mup6*F(8430)-4*(mup2mup5*p2mup6*F(8418)+mup2m
     -   up5*p5mup6*F(8419)-p4mup2*p4mup5*p6mup6*F(8422)-p4mup5*p5mup
     -   2*p6mup6*F(8424)-p1mup6*(mup2mup5*F(8415)-4*(p3mup2*(p3mup5*
     -   F(6803)+p4mup5*F(6895))-p3mup5*p5mup2*F(8420)-p4mup5*p5mup2*
     -   F(8426)-p1mup2*p3mup5*F(8431)-p1mup2*p4mup5*F(8434)-p3mup5*p
     -   4mup2*F(8442)-p4mup2*p4mup5*F(8443)-p2mup5*p6mup2*F(9628))))
     -   -8*(p1mup2*p3mup5*p6mup6*F(8425)-p3mup5*p4mup2*p6mup6*F(8437
     -   )+p3mup5*p5mup2*p6mup6*F(8439)-p2mup2*(p3mup5*(p6mup6*F(5664
     -   )+2*(p5mup6*F(8428)-p1mup6*F(8440)-p2mup6*F(8444)))-2*p4mup5
     -   *(p1mup6*F(8441)+p2mup6*F(8445)))+p2mup5*p6mup2*p6mup6*F(961
     -   8)-p3mup2*p3mup5*p6mup6*F(9637))+16*(p4mup2*p4mup5*p5mup6*F(
     -   8421)+p4mup5*p5mup2*p5mup6*F(8423)+p2mup5*p5mup6*p6mup2*F(84
     -   27)-p3mup2*p4mup5*p5mup6*F(8429)+p1mup2*p2mup6*p3mup5*F(8432
     -   )+p1mup2*p3mup5*p5mup6*F(8433)+p1mup2*p2mup6*p4mup5*F(8435)+
     -   p3mup5*p4mup2*p5mup6*F(8436)+p3mup5*p5mup2*p5mup6*F(8438)+p2
     -   mup6*p3mup5*p5mup2*F(8446)+p3mup2*p3mup5*p5mup6*F(8447)+p2mu
     -   p6*p3mup5*p4mup2*F(8448)+p2mup6*p4mup2*p4mup5*F(8449)+p2mup6
     -   *p4mup5*p5mup2*F(8450)-p2mup6*p3mup2*p3mup5*F(9673)-p2mup6*p
     -   3mup2*p4mup5*F(9674))
       Fa30(712) = p4mup6*(p1mup2*(p5mup5*F(8451)-2*(p3mup5*F(8453)-p2
     -   mup5*F(8455)+p1mup5*F(8457)))+p2mup2*(p5mup5*F(8452)-2*(p3mu
     -   p5*F(8454)-p2mup5*F(8456)+p1mup5*F(8458))))
       Fa30(713) = p6mup2*(p5mup5*p5mup6*F(7701)+2*p2mup6*p5mup5*F(771
     -   2)+4*p2mup5*p2mup6*F(9635)+p3mup5*(-2*p6mup6*F(9624)+4*(p5mu
     -   p6*F(6934)+p1mup6*F(9630)+p2mup6*F(9637)))+p4mup5*(p6mup6*F(
     -   4474)+4*(p5mup6*F(7707)+p1mup6*F(9631)+p2mup6*F(9638))))
       Fa30(714) = -(p5mup5*p6mup2*p6mup6*F(7744))-mup2mup5*p3mup6*F(8
     -   459)-mup2mup5*p4mup6*F(8460)-2*(p3mup2*p3mup6*p5mup5*F(7692)
     -   -p1mup6*p5mup5*p6mup2*F(7719)+p4mup2*p4mup6*p5mup5*F(8464)+p
     -   4mup6*p5mup2*p5mup5*F(8466)+p1mup2*p3mup6*p5mup5*F(8468)-p3m
     -   up2*p4mup6*p5mup5*F(8473)+p3mup6*p5mup2*p5mup5*F(8476)-p3mup
     -   6*p4mup2*p5mup5*F(8477)+p2mup2*(-2*p4mup6*p6mup5*F(8461)+p3m
     -   up6*(p5mup5*F(8471)-2*(p6mup5*F(5664)-p2mup5*F(8444)-p1mup5*
     -   F(8480)-p4mup5*F(8484)))))+p1mup5*(2*p6mup2*(-(p6mup6*F(9626
     -   ))+2*(p5mup6*F(6938)+p2mup6*F(9628)+p1mup6*F(9656)))-4*(p4mu
     -   p6*(p5mup2*F(8469)-p4mup2*F(8478)+p3mup2*F(8481))+p3mup6*(p5
     -   mup2*F(8463)-p1mup2*F(8479)-p4mup2*F(8482)+p3mup2*F(9663))))
     -   -4*(p1mup2*p3mup6*p6mup5*F(8425)-p1mup2*p2mup5*p3mup6*F(8432
     -   )-p3mup6*p4mup2*p6mup5*F(8437)+p3mup6*p5mup2*p6mup5*F(8439)-
     -   p2mup5*p3mup6*p5mup2*F(8446)-p2mup5*p3mup6*p4mup2*F(8448)-p1
     -   mup2*p4mup6*p6mup5*F(8462)+p4mup2*p4mup6*p6mup5*F(8465)+p4mu
     -   p6*p5mup2*p6mup5*F(8467)-p2mup5*p4mup2*p4mup6*F(8470)-p2mup5
     -   *p4mup6*p5mup2*F(8472)-p3mup2*p4mup6*p6mup5*F(8474)-p1mup2*p
     -   3mup6*p4mup5*F(8475)+p2mup5*p3mup2*p4mup6*F(8483)-p3mup2*p3m
     -   up6*p6mup5*F(9637)+p2mup5*p3mup2*p3mup6*F(9673))
       Fa30(715) = p3mup6*p5mup5*p6mup2*F(7717)+p4mup6*p5mup5*p6mup2*F
     -   (8485)-2*(p4mup6*p6mup2*p6mup5*F(7733)+p3mup5*p4mup2*p4mup6*
     -   F(8486)-p4mup2*p4mup5*p4mup6*F(8487)-p3mup6*p4mup2*p4mup5*F(
     -   8488)-p3mup5*p3mup6*p5mup2*F(8489)+p2mup6*p6mup2*p6mup5*F(96
     -   18)+p5mup6*p6mup2*p6mup5*F(9621)+p3mup6*p6mup2*p6mup5*F(9624
     -   )+p1mup6*p6mup2*p6mup5*F(9626)-p1mup5*p6mup2*(p4mup6*F(7704)
     -   +p3mup6*F(9630))-p2mup5*p6mup2*(p4mup6*F(7702)+p3mup6*F(9637
     -   )))
       Fa30(716) = -(p3mup5*(p4mup6*(p5mup2*F(8491)-p6mup2*F(8492)+p3m
     -   up2*F(8495))+p3mup6*(-(p1mup2*F(8493))+p2mup2*F(8494)+p4mup2
     -   *F(8497)-p6mup2*F(9643)+p3mup2*F(9676))))+p4mup5*(p4mup6*(p6
     -   mup2*F(7714)+p5mup2*F(8490)-p3mup2*F(8496))+p3mup6*(p5mup2*F
     -   (8498)+p6mup2*F(9642)-p3mup2*F(9678)))
       Fa30(717) = mup3mup6*p4mup2*F(8502)+mup3mup6*p5mup2*F(8503)+p2m
     -   up2*(mup3mup6*F(8499)+8*(p5mup3*p6mup6*F(8507)-2*(p4mup3*p6m
     -   up6*F(8511)-p2mup6*p5mup3*F(8513)+p5mup3*p5mup6*F(8517)+p1mu
     -   p6*p5mup3*F(8524)-p2mup3*(p6mup6*F(4110)+2*(p1mup6*F(4208)+p
     -   2mup6*F(4220)+p5mup6*F(8538))))+4*(p1mup6*p4mup3*F(4212)+p2m
     -   up6*p4mup3*F(4224)+p4mup3*p5mup6*F(8543))+p3mup3*(-2*p6mup6*
     -   F(4112)+4*(p1mup6*F(4210)+p2mup6*F(4223)+p5mup6*F(8545)))))-
     -   16*(p3mup3*p4mup2*p6mup6*F(5732)-p1mup6*(-(p5mup2*(p4mup3*F(
     -   8526)+p2mup3*F(8529)+p5mup3*F(8530)+p3mup3*F(8531)))+mup2mup
     -   3*F(8551)+p3mup2*(-(p5mup3*F(8527))+2*(p4mup3*F(4218)+p2mup3
     -   *F(8554)+p3mup3*F(8555)))+p4mup2*(p5mup3*F(8525)-2*(p2mup3*F
     -   (8528)+p3mup3*F(8532)-p4mup3*F(8556)))))+4*(mup2mup6*p4mup3*
     -   F(8500)+mup2mup6*p5mup3*F(8501)+mup2mup6*p3mup3*F(8504)+mup2
     -   mup6*p2mup3*F(8505)-2*(p4mup2*p5mup3*p6mup6*F(8508)-p4mup3*p
     -   5mup2*p6mup6*F(8509)-p2mup3*p5mup2*p6mup6*F(8514)-p5mup2*p5m
     -   up3*p6mup6*F(8518)+p3mup3*p5mup2*p6mup6*F(8520)+mup2mup3*p6m
     -   up6*F(8523))+4*(p2mup3*p4mup2*p6mup6*F(8512)+p4mup2*p5mup3*p
     -   5mup6*F(8521)-p4mup3*p5mup2*p5mup6*F(8522)-p2mup3*p5mup2*p5m
     -   up6*F(8542)+p4mup2*p4mup3*p6mup6*F(8546)-p5mup2*p5mup3*p5mup
     -   6*F(8547)-p3mup3*p5mup2*p5mup6*F(8548)+mup2mup3*p5mup6*F(855
     -   3)+p2mup6*(-(p5mup2*(p4mup3*F(8516)+p2mup3*F(8535)+p5mup3*F(
     -   8537)+p3mup3*F(8539)))+mup2mup3*F(8552)+p3mup2*(-(p5mup3*F(8
     -   533))+2*(p4mup3*F(4234)+p2mup3*F(8557)+p3mup3*F(8558)))+p4mu
     -   p2*(p5mup3*F(8515)-2*(p2mup3*F(8534)+p3mup3*F(8540)-p4mup3*F
     -   (8560)))))+p3mup2*(-(mup3mup6*F(8506))+2*p5mup3*p6mup6*F(851
     -   0)+4*(p3mup3*p6mup6*F(5893)+p4mup3*p6mup6*F(8519)+p2mup3*p6m
     -   up6*F(8536)-p5mup3*p5mup6*F(8541))+8*(p4mup3*p5mup6*F(8549)+
     -   p2mup3*p5mup6*F(8559)+p3mup3*p5mup6*F(8561)))-8*(p2mup3*p4mu
     -   p2*p5mup6*F(8544)+p3mup3*p4mup2*p5mup6*F(8550)-p4mup2*p4mup3
     -   *p5mup6*F(8562)))
       Fa30(718) = p1mup2*(mup3mup6*F(8563)+2*p5mup3*p6mup6*F(8565)-4*
     -   (p3mup3*p6mup6*F(4480)-p2mup3*p6mup6*F(8566)-p2mup6*p5mup3*F
     -   (8567)+p4mup3*p6mup6*F(8568)+p5mup3*p5mup6*F(8569)+p1mup3*p6
     -   mup6*F(8570))+8*(p1mup3*p2mup6*F(4207)+p2mup6*p4mup3*F(4211)
     -   +p2mup3*p2mup6*F(4221)+p2mup6*p3mup3*F(4222)+p1mup3*p5mup6*F
     -   (8576)+p4mup3*p5mup6*F(8578)+p2mup3*p5mup6*F(8581)+p3mup3*p5
     -   mup6*F(8582))+p1mup6*(-4*p5mup3*F(8572)+8*(p2mup3*F(4207)+p3
     -   mup3*F(4209)+p4mup3*F(4215)+p1mup3*F(8583))))+p1mup3*(mup2mu
     -   p6*F(8564)+2*p5mup2*p6mup6*F(8574)+p2mup2*(-4*p6mup6*F(8571)
     -   +8*(p2mup6*F(4208)+p5mup6*F(8577)+p1mup6*F(8584)))-4*(p2mup6
     -   *p5mup2*F(8529)-p4mup2*p6mup6*F(8573)-p3mup2*p6mup6*F(8575)+
     -   p5mup2*p5mup6*F(8580)+p1mup6*p5mup2*F(8586))-8*(p2mup6*p4mup
     -   2*F(8528)-p2mup6*p3mup2*F(8554)+p4mup2*p5mup6*F(8579)+p1mup6
     -   *p4mup2*F(8585)-p1mup6*p3mup2*F(8587)-p3mup2*p5mup6*F(8588))
     -   )
       Fa30(719) = p1mup2*(p4mup6*p5mup3*F(8591)+p6mup3*(-(p6mup6*F(85
     -   89))+2*(p4mup6*F(8593)+p5mup6*F(8594))))+p2mup2*(-(p6mup3*p6
     -   mup6*F(8590))+p4mup6*p5mup3*F(8592)+2*(p4mup6*p6mup3*F(4115)
     -   +p5mup6*p6mup3*F(8595)))
       Fa30(720) = p4mup3*p6mup2*(p6mup6*F(8596)+2*(p5mup6*F(8597)+p1m
     -   up6*F(8598)))
       Return
       End
