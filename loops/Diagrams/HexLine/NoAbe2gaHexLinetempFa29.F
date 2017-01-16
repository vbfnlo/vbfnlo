c************** Calling the Fa functions*************************
       subroutine NoAbe2gaHexLineFa29(p1mup2,p1mup3,p1mup5,p1mup6,p2mu
     -   p2,p2mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p
     -   4mup3,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup
     -   3,p6mup5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6
     -   ,mup5mup6,Fa29)
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
       Complex*16 Fa29(673:696)
       COMMON/NoAbe2gaHexLineFaFunctions/Fa
       Fa29(673) = p5mup5*(p1mup6*p5mup3*F(7648)-2*(p1mup6*p2mup3*F(76
     -   49)-p1mup3*p2mup6*F(7650)))-p5mup3*p5mup5*p5mup6*F(7653)+p1m
     -   up3*p5mup5*p5mup6*F(7658)+p2mup3*p5mup5*p5mup6*F(7660)+2*(p1
     -   mup3*p3mup6*p5mup5*F(7655)+p1mup6*(p5mup5*(p1mup3*F(7654)+p4
     -   mup3*F(7666)+p3mup3*F(7671))-2*(p2mup5*(p5mup3*F(7652)+p2mup
     -   3*F(7661)+p1mup3*F(7663)+p4mup3*F(7680)-p3mup3*F(7682))+p3mu
     -   p5*(p5mup3*F(7657)-p2mup3*F(7662)+p1mup3*F(7669)+p4mup3*F(76
     -   84)-p3mup3*F(7686))))+p2mup3*p2mup6*p5mup5*F(7689)+p2mup3*p3
     -   mup6*p5mup5*F(7692))-4*(p1mup3*p2mup6*p3mup5*F(7664)+p1mup3*
     -   p2mup5*p3mup6*F(7664)+p1mup3*p2mup5*p5mup6*F(7665)+p1mup3*p3
     -   mup5*p5mup6*F(7667)+p1mup3*p2mup5*p2mup6*F(7668)+p1mup3*p3mu
     -   p5*p3mup6*F(7670)-p2mup3*p2mup5*p5mup6*F(7688)+p2mup5*p2mup6
     -   *p5mup3*F(7690)-p2mup5*p3mup3*p5mup6*F(7691)-p2mup3*p3mup5*p
     -   5mup6*F(7691)-p3mup3*p3mup5*p3mup6*F(9666)-p1mup5*(p3mup6*p5
     -   mup3*F(7656)+p5mup3*p5mup6*F(7659)+p1mup6*p2mup3*F(7672)+p4m
     -   up3*p5mup6*F(7674)-p1mup3*p1mup6*F(7675)-p1mup3*p3mup6*F(767
     -   6)-p1mup3*p5mup6*F(7677)+p1mup6*p5mup3*F(7678)+p2mup3*p5mup6
     -   *F(7681)+p3mup6*p4mup3*F(7683)-p1mup6*p4mup3*F(7687)+p2mup6*
     -   (p5mup3*F(7651)-p1mup3*F(7673)+p4mup3*F(7679)+p2mup3*F(9662)
     -   )+p2mup3*p3mup6*F(9663)+p3mup3*(p5mup6*F(7163)+p1mup6*F(7685
     -   )+p2mup6*F(9663)+p3mup6*F(9667)))-p2mup3*p2mup5*p2mup6*F(967
     -   2)-p2mup5*p2mup6*p3mup3*F(9673)-p2mup3*p2mup6*p3mup5*F(9673)
     -   -p2mup3*p2mup5*p3mup6*F(9673)-(p2mup6*p3mup3*p3mup5+p2mup5*p
     -   3mup3*p3mup6+p2mup3*p3mup5*p3mup6)*F(9676))
       Fa29(674) = p3mup3*(4*p3mup5*p5mup6*F(7084)+p5mup5*(p5mup6*F(76
     -   95)+2*(p2mup6*F(7692)+p3mup6*F(7699))))+2*p5mup3*(-(p2mup6*p
     -   5mup5*F(7694))+p3mup6*p5mup5*F(7697)-2*(-(p2mup5*p5mup6*F(76
     -   93))+p3mup5*p5mup6*F(7696)+p2mup6*p3mup5*F(7698)+p2mup5*p3mu
     -   p6*F(7698)+p3mup5*p3mup6*F(7700)))
       Fa29(675) = p6mup3*(p5mup5*p5mup6*F(7701)+4*(-(p3mup5*p4mup6*F(
     -   7703))+p1mup5*p4mup6*F(7704)+p1mup5*p2mup6*F(9628)+p1mup6*p3
     -   mup5*F(9630)+p1mup5*p3mup6*F(9630)+p1mup6*p4mup5*F(9631)+p2m
     -   up5*(p4mup6*F(7702)+p1mup6*F(9628)+p2mup6*F(9635)+p3mup6*F(9
     -   637))+p2mup6*p4mup5*F(9638)+p3mup6*p4mup5*F(9642)+p1mup5*p1m
     -   up6*F(9656)))
       Fa29(676) = p2mup3*(4*p6mup5*(p4mup6*F(7702)+p1mup6*F(7705)+p5m
     -   up6*F(7713)+p2mup6*F(9635)+p3mup6*F(9637))+p6mup6*(-(p4mup5*
     -   F(4477))+p5mup5*F(7712)+2*(p1mup5*F(9628)+p2mup5*F(9635)+p3m
     -   up5*F(9637))))+2*(p6mup3*(-(p4mup6*p5mup5*F(7708))+2*(p4mup5
     -   *p5mup6*F(7707)-p2mup5*p5mup6*F(7711)))+p6mup3*(p2mup6*p5mup
     -   5*F(7712)+p3mup6*p5mup5*F(7717)+p1mup6*p5mup5*F(7719)+2*(p4m
     -   up5*p4mup6*F(7714)-p3mup5*p5mup6*F(7716)-p1mup5*p5mup6*F(771
     -   8)+p2mup6*p3mup5*F(9637)))-2*(p1mup3*p6mup5*(p2mup6*F(7706)+
     -   p1mup6*F(7709))+p4mup3*p4mup6*p6mup5*F(7710)+p5mup3*p5mup6*p
     -   6mup5*F(7715)-p3mup6*(p3mup5*p6mup3+p3mup3*p6mup5)*F(9643)))
       Fa29(677) = p6mup5*(p4mup6*p5mup3*F(7721)+p1mup3*(p3mup6*F(7720
     -   )+p4mup6*F(7722)+p5mup6*F(7723))-p1mup6*p5mup3*F(7724)+p2mup
     -   6*p4mup3*F(7725)+p2mup6*p5mup3*F(7726)+p1mup6*p4mup3*F(7727)
     -   -p3mup6*p4mup3*F(7728)+p3mup6*p5mup3*F(7730)+p4mup3*p5mup6*F
     -   (7731)+p3mup3*(p4mup6*F(7703)-p5mup6*F(7729)+p1mup6*F(7732)-
     -   p2mup6*F(9637)))
       Fa29(678) = p6mup3*p6mup5*(p4mup6*F(7733)+p2mup6*F(9618)+p5mup6
     -   *F(9621)+p3mup6*F(9624)+p1mup6*F(9626))
       Fa29(679) = p6mup6*(p2mup5*p5mup3*F(7726)+p2mup3*p6mup5*F(9618)
     -   )
       Fa29(680) = p6mup6*(p4mup3*p5mup5*F(7737)+p4mup5*p5mup3*F(7738)
     -   -p4mup3*p4mup5*F(7739)+p5mup3*p5mup5*F(7740)+p1mup3*(-(p5mup
     -   5*F(7734))+2*(p2mup5*F(7706)+p3mup5*F(7720)+p1mup5*F(7736))+
     -   p4mup5*F(7741))+2*(p2mup5*p4mup3*F(7725)-p3mup5*p4mup3*F(772
     -   8)+p3mup5*p5mup3*F(7730)-p1mup5*p5mup3*F(7735)+p1mup5*p4mup3
     -   *F(7742))-p3mup3*(p4mup5*F(4492)+p5mup5*F(7717)+2*(p1mup5*F(
     -   9630)+p2mup5*F(9637)+p3mup5*F(9643))))
       Fa29(681) = p6mup3*p6mup6*(p4mup5*F(7743)+4*p5mup5*F(7744)+8*(p
     -   2mup5*F(9618)+p3mup5*F(9624)+p1mup5*F(9626)))
       Fa29(682) = p6mup5*p6mup6*(p1mup3*F(7745)-p4mup3*F(7746)+p5mup3
     -   *F(7747)+p6mup3*F(9616)-p3mup3*F(9624))
       Fa29(683) = p5mup5*p6mup3*F(7750)+p1mup5*p6mup3*F(7753)+p2mup3*
     -   p5mup5*F(7757)+p1mup5*p2mup3*F(7760)+p2mup5*p6mup3*F(7761)+p
     -   4mup3*p4mup5*F(7763)+p4mup5*p6mup3*F(7766)+p1mup3*(p5mup5*F(
     -   7749)+p1mup5*F(7751)-4*p6mup5*F(7754)+p2mup5*F(7759)+p4mup5*
     -   F(7765)+p3mup5*F(7767))+p5mup3*(p5mup5*F(7748)+p1mup5*F(7752
     -   )-4*p6mup5*F(7755)+p2mup5*F(7758)+p4mup5*F(7764)+p3mup5*F(77
     -   68))+p3mup5*p6mup3*F(7769)+p4mup3*p5mup5*F(7770)+p1mup5*p4mu
     -   p3*F(7771)+p4mup3*p6mup5*F(7772)+p2mup3*p2mup5*F(7773)+mup3m
     -   up5*F(7774)+p3mup3*p4mup5*F(7775)+p2mup3*p3mup5*F(7776)+p3mu
     -   p5*p4mup3*F(7777)+p3mup3*p5mup5*F(7778)+p1mup5*p3mup3*F(7779
     -   )-4*(p6mup3*p6mup5*F(7756)-p2mup3*p6mup5*F(7762)+p3mup3*p6mu
     -   p5*F(7780))+p2mup3*p4mup5*F(7781)+p2mup5*p4mup3*F(7782)+p3mu
     -   p3*p3mup5*F(7783)+p2mup5*p3mup3*F(7784)
       Fa29(684) = p4mup3*F(7785)+p5mup3*F(7786)+2*(p1mup3*F(7787)-p6m
     -   up3*F(7788)+p3mup3*F(7789)+p2mup3*F(7790))
       Fa29(685) = p5mup3*F(7792)+2*(p4mup3*F(7791)-p1mup3*F(7793)-p6m
     -   up3*F(7794)+p3mup3*F(7795)+p2mup3*F(7796))
       Fa29(686) = mup5mup6*F(7797)+p1mup6*p2mup5*F(7800)+p2mup5*p6mup
     -   6*F(7803)+2*(p1mup6*p5mup5*F(7801)+p2mup6*p5mup5*F(7804))+p2
     -   mup6*p6mup5*F(7805)+p2mup5*p2mup6*F(7809)+p2mup5*p3mup6*F(78
     -   14)+p4mup5*(p6mup6*F(7798)+p1mup6*F(7807)+p2mup6*F(7810)+p3m
     -   up6*F(7815)+p4mup6*F(7818)+p5mup6*F(7819))+p4mup6*p5mup5*F(7
     -   820)+p2mup5*p4mup6*F(7822)+p1mup5*(p2mup6*F(7808)-8*(p1mup6*
     -   F(7799)-p3mup6*F(7813)-p5mup6*F(7824))+4*(p4mup6*F(7821)+p6m
     -   up6*F(7825)))+p2mup5*p5mup6*F(7826)+p5mup5*p6mup6*F(7828)+p2
     -   mup6*p3mup5*F(7830)-4*(p6mup5*p6mup6*F(7806)+p3mup5*p6mup6*F
     -   (7812)-p3mup6*p5mup5*F(7816)-p4mup6*p6mup5*F(7823)-p5mup5*p5
     -   mup6*F(7827)-p3mup5*p4mup6*F(7832))+8*(p1mup6*p6mup5*F(7802)
     -   +p1mup6*p3mup5*F(7811)+p3mup6*p6mup5*F(7817)+p5mup6*p6mup5*F
     -   (7829)-p3mup5*p3mup6*F(7831)-p3mup5*p5mup6*F(7833))
       Fa29(687) = mup5mup6*F(7834)+p2mup6*p6mup5*F(7843)+p5mup5*(p2mu
     -   p6*F(7842)-2*(p5mup6*F(7844)-p6mup6*F(7845))+p4mup6*F(7847))
     -   +p2mup5*p2mup6*F(7851)+2*p2mup5*p6mup6*F(7853)+p2mup5*p4mup6
     -   *F(7854)+p2mup6*p3mup5*F(7856)+p3mup5*p4mup6*F(7857)+p1mup5*
     -   (p2mup6*F(7838)+4*p6mup6*F(7840)+p4mup6*F(7848)-8*(p1mup6*F(
     -   7837)+p5mup6*F(7839)+p3mup6*F(7860)))-4*(p6mup5*p6mup6*F(784
     -   6)-p4mup6*p6mup5*F(7849)+p1mup6*p2mup5*F(7850)-p2mup5*p5mup6
     -   *F(7852)-p3mup5*p6mup6*F(7859)+p3mup6*p5mup5*F(7861)-p1mup6*
     -   p5mup5*F(7865))+p2mup5*p3mup6*F(7867)+8*(p1mup6*p6mup5*F(784
     -   1)-p1mup6*p3mup5*F(7855)+p3mup5*p5mup6*F(7858)+p3mup6*p6mup5
     -   *F(7859)+p5mup6*p6mup5*F(7866)-p3mup5*p3mup6*F(7868))+p4mup5
     -   *(p4mup6*F(7862)+p2mup6*F(7864)-2*(p5mup6*F(7835)+p6mup6*F(7
     -   836)+p1mup6*F(7863)+p3mup6*F(7869)))
       Fa29(688) = p1mup3*p6mup6*F(7874)+p6mup3*p6mup6*F(7875)+p2mup3*
     -   p6mup6*F(7880)+p3mup3*p6mup6*F(7883)+p1mup6*(p5mup3*F(7871)+
     -   2*(p1mup3*F(7873)+p6mup3*F(7876)+p2mup3*F(7879)+p3mup3*F(788
     -   2)-p4mup3*F(7885)))+p4mup3*p5mup6*F(7887)+p4mup3*p6mup6*F(78
     -   88)+p5mup3*(p6mup6*F(7872)+p2mup6*F(7877)+p5mup6*F(7884)+p3m
     -   up6*F(7891)+p4mup6*F(7893))+2*(mup3mup6*F(7870)+p1mup3*p2mup
     -   6*F(7878)-p2mup6*p6mup3*F(7881)+p1mup3*p3mup6*F(7882)+p2mup6
     -   *p4mup3*F(7886)-p1mup3*p5mup6*F(7889)-p5mup6*p6mup3*F(7890)-
     -   p3mup6*p6mup3*F(7892)-p4mup6*p6mup3*F(7894)+p1mup3*p4mup6*F(
     -   7895)+p2mup3*p2mup6*F(7896)+p2mup6*p3mup3*F(7897)+p3mup3*p5m
     -   up6*F(7898)+p3mup6*p4mup3*F(7899)+p3mup3*p3mup6*F(7900)+p4mu
     -   p3*p4mup6*F(7901)+p2mup3*p5mup6*F(7902)+p2mup3*p3mup6*F(7903
     -   )+p2mup3*p4mup6*F(7904)+p3mup3*p4mup6*F(7905))
       Fa29(689) = mup3mup6*F(7906)+p5mup3*(p5mup6*F(7907)+p1mup6*F(79
     -   08)+2*p6mup6*F(7909)+p2mup6*F(7914)+p3mup6*F(7925)+p4mup6*F(
     -   7928))+2*(-(p6mup3*p6mup6*F(7913))+p2mup3*p6mup6*F(7916)+p3m
     -   up3*p6mup6*F(7918)+p1mup6*p4mup3*F(7919)+p2mup6*p4mup3*F(792
     -   0)+p5mup6*(-(p6mup3*F(7912))+p2mup3*F(7915)+p3mup3*F(7917)+p
     -   4mup3*F(7921))+p4mup3*p6mup6*F(7922)+p1mup6*p6mup3*F(7926)-p
     -   2mup6*p6mup3*F(7927)-p1mup6*p3mup3*F(7929)+p2mup6*p3mup3*F(7
     -   930)-p3mup6*p4mup3*F(7931)-p3mup6*p6mup3*F(7932)+p4mup3*p4mu
     -   p6*F(7933)+p1mup6*p2mup3*F(7934)+p2mup3*p2mup6*F(7935)-p4mup
     -   6*p6mup3*F(7936)+p1mup3*(p5mup6*F(7910)+p6mup6*F(7911)+p1mup
     -   6*F(7923)-p2mup6*F(7924)-p3mup6*F(7929)-p4mup6*F(7937))+p2mu
     -   p3*p4mup6*F(7938)+p3mup3*p3mup6*F(7939)+p2mup3*p3mup6*F(7940
     -   )+p3mup3*p4mup6*F(7941))
       Fa29(690) = p1mup6*F(7942)+p2mup6*F(7943)+p5mup6*F(7944)+p6mup6
     -   *F(7945)+p3mup6*F(7946)+p4mup6*F(7947)
       Fa29(691) = mup2mup6*F(7948)+p5mup6*p6mup2*F(7954)+p1mup6*p6mup
     -   2*F(7955)+p2mup2*p6mup6*F(7957)+p2mup6*p6mup2*F(7960)+p4mup2
     -   *p6mup6*F(7965)+p6mup2*p6mup6*F(7967)+p3mup6*p6mup2*F(7970)+
     -   p5mup2*(p5mup6*F(7949)+p1mup6*F(7950)+p6mup6*F(7951)+p2mup6*
     -   F(7958)+p3mup6*F(7968)+p4mup6*F(7971))+p1mup2*(p5mup6*F(7953
     -   )+p6mup6*F(7966)+2*(p1mup6*F(7952)-p2mup6*F(7959)-p3mup6*F(7
     -   969)-p4mup6*F(7972)))+p4mup6*p6mup2*F(7973)+p2mup2*p2mup6*F(
     -   7974)+p3mup2*p6mup6*F(7975)+p3mup6*p4mup2*F(7976)+p2mup2*p3m
     -   up6*F(7977)+p4mup2*p4mup6*F(7978)+p2mup6*p4mup2*F(7980)+p2mu
     -   p2*p4mup6*F(7981)+2*(p2mup2*p5mup6*F(7956)+p1mup6*p2mup2*F(7
     -   961)+p1mup6*p3mup2*F(7962)+p1mup6*p4mup2*F(7963)+p4mup2*p5mu
     -   p6*F(7964)+p3mup2*p5mup6*F(7979)+p3mup2*p4mup6*F(7982)+p3mup
     -   2*p3mup6*F(7983)+p2mup6*p3mup2*F(7984))
       Fa29(692) = mup2mup6*F(7985)+p5mup2*p6mup6*F(7990)+p1mup2*p6mup
     -   6*F(7992)+p6mup2*p6mup6*F(7993)+p2mup6*p5mup2*F(7996)+p1mup2
     -   *p2mup6*F(7997)+p1mup6*(p2mup2*F(7986)+p5mup2*F(7989)+p1mup2
     -   *F(7991)-2*(p6mup2*F(7994)-p3mup2*F(7999))+p4mup2*F(8000))+p
     -   4mup2*p6mup6*F(8001)+p2mup6*p4mup2*F(8002)+p3mup6*p5mup2*F(8
     -   003)+p1mup2*p3mup6*F(8005)+p4mup6*p5mup2*F(8007)+p1mup2*p4mu
     -   p6*F(8008)+p2mup2*(p5mup6*F(7987)+p6mup6*F(7988)+p2mup6*F(79
     -   95)+p3mup6*F(8004)+p4mup6*F(8009))+p5mup2*p5mup6*F(8011)+p3m
     -   up2*p6mup6*F(8015)+p2mup6*p3mup2*F(8018)+2*(p2mup6*p6mup2*F(
     -   7998)+p3mup6*p6mup2*F(8006)+p4mup6*p6mup2*F(8010)+p1mup2*p5m
     -   up6*F(8012)+p5mup6*p6mup2*F(8013)+p4mup2*p5mup6*F(8014)-p3mu
     -   p6*p4mup2*F(8016)+p4mup2*p4mup6*F(8017)+p3mup2*p5mup6*F(8019
     -   )+p3mup2*p3mup6*F(8020)+p3mup2*p4mup6*F(8021))
       Fa29(693) = mup2mup6*F(8022)+p5mup2*p5mup6*F(8026)+p5mup2*p6mup
     -   6*F(8028)+p1mup2*p5mup6*F(8029)+p1mup2*p6mup6*F(8030)+p2mup6
     -   *p5mup2*F(8034)+p4mup2*p5mup6*F(8037)+p4mup2*p6mup6*F(8038)+
     -   p1mup2*p2mup6*F(8041)+p1mup6*(p2mup2*F(8023)+p5mup2*F(8027)+
     -   p3mup2*F(8035)+p1mup2*F(8040)+2*(p4mup2*F(8036)+p6mup2*F(804
     -   2)))+p3mup6*p5mup2*F(8044)+p4mup6*p5mup2*F(8046)+p2mup2*(p5m
     -   up6*F(8024)+p6mup6*F(8025)+p2mup6*F(8033)+p3mup6*F(8045)+p4m
     -   up6*F(8047))+p3mup2*p6mup6*F(8048)+p2mup6*p3mup2*F(8049)+p3m
     -   up2*p5mup6*F(8050)+p3mup2*p3mup6*F(8052)+p3mup2*p4mup6*F(805
     -   3)+p1mup2*p3mup6*F(8055)+p1mup2*p4mup6*F(8056)+2*(p5mup6*p6m
     -   up2*F(8031)+p6mup2*p6mup6*F(8032)-p2mup6*p4mup2*F(8039)+p2mu
     -   p6*p6mup2*F(8043)-p4mup2*p4mup6*F(8051)+p3mup6*p4mup2*F(8054
     -   )+p4mup6*p6mup2*F(8057)+p3mup6*p6mup2*F(8058))
       Fa29(694) = p4mup3*p6mup5*F(8072)+p5mup3*p6mup5*F(8073)+p5mup5*
     -   (p6mup3*F(8060)-2*(p2mup3*F(8063)+p3mup3*F(8068))+p4mup3*F(8
     -   076)+p5mup3*F(8077))+2*p1mup3*(p5mup5*F(8059)+2*(p6mup5*F(80
     -   62)+p1mup5*F(8069)+p4mup5*F(8079)-p3mup5*F(8081)-p2mup5*F(80
     -   87)))+4*(p1mup5*p6mup3*F(8061)+p2mup5*p6mup3*F(8064)+p3mup5*
     -   p6mup3*F(8066)+p1mup5*p4mup3*F(8067)+p1mup5*p5mup3*F(8070)-p
     -   6mup5*(p2mup3*F(8065)+p3mup3*F(8071))-p6mup3*p6mup5*F(8074)+
     -   p4mup5*p6mup3*F(8075)+p4mup5*p5mup3*F(8078)+mup3mup5*F(8080)
     -   -p1mup5*p3mup3*F(8081)-p3mup5*p4mup3*F(8082)-p3mup5*p5mup3*F
     -   (8083)+p4mup3*p4mup5*F(8084)-p2mup5*p4mup3*F(8085)-p2mup5*p5
     -   mup3*F(8086)-p1mup5*p2mup3*F(8088)-p3mup3*p3mup5*F(8089)+p2m
     -   up3*p4mup5*F(8090)-p2mup3*p3mup5*F(8091)-p2mup5*p3mup3*F(809
     -   2)-p3mup3*p4mup5*F(8093)+p2mup3*p2mup5*F(8094))
       Fa29(695) = p4mup3*p5mup5*F(8113)+p4mup3*p6mup5*F(8117)+p5mup3*
     -   (p5mup5*F(8095)+4*(p1mup5*F(8099)+p6mup5*F(8102)-p2mup5*F(81
     -   05)-p3mup5*F(8111)-p4mup5*F(8118)))+2*(p5mup5*p6mup3*F(8097)
     -   -p2mup3*p5mup5*F(8104)-p3mup3*p5mup5*F(8115)+p1mup3*(p5mup5*
     -   F(8096)+2*(p1mup5*F(8098)+p6mup5*F(8101)+p2mup5*F(8106)+p3mu
     -   p5*F(8110)+p4mup5*F(8119))))+4*(p1mup5*p6mup3*F(8100)+p6mup3
     -   *p6mup5*F(8103)+p1mup5*p2mup3*F(8107)-p2mup5*p6mup3*F(8108)-
     -   p2mup3*p6mup5*F(8109)+p1mup5*p3mup3*F(8110)+p3mup5*p6mup3*F(
     -   8112)-p1mup5*p4mup3*F(8114)-p3mup3*p6mup5*F(8116)+p4mup5*p6m
     -   up3*F(8120)+p2mup3*p2mup5*F(8121)+mup3mup5*F(8122)+p3mup5*p4
     -   mup3*F(8123)+p3mup3*p3mup5*F(8124)+p4mup3*p4mup5*F(8125)+p2m
     -   up3*p3mup5*F(8126)+p2mup5*p4mup3*F(8127)+p2mup5*p3mup3*F(812
     -   8)+p2mup3*p4mup5*F(8129)+p3mup3*p4mup5*F(8130))
       Fa29(696) = p1mup5*F(8131)+p2mup5*F(8132)+p5mup5*F(8133)+p6mup5
     -   *F(8134)+p4mup5*F(8135)+p3mup5*F(8136)
       Return
       End
