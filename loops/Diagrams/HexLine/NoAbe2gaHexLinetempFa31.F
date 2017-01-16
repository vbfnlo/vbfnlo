c************** Calling the Fa functions*************************
       subroutine NoAbe2gaHexLineFa31(p1mup2,p1mup3,p1mup5,p1mup6,p2mu
     -   p2,p2mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p
     -   4mup3,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup
     -   3,p6mup5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6
     -   ,mup5mup6,Fa31)
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
       Complex*16 Fa31(721:744)
       COMMON/NoAbe2gaHexLineFaFunctions/Fa
       Fa31(721) = p1mup2*(-(p3mup6*p5mup3*F(8599))+2*(p3mup6*p4mup3*F
     -   (4217)+p4mup3*p4mup6*F(4219)+p2mup3*p3mup6*F(4222)+p2mup3*p4
     -   mup6*F(4229)-p3mup6*p6mup3*F(4480)+p1mup3*(p3mup6*F(4209)+p4
     -   mup6*F(8603))))+p2mup2*(p3mup6*(-2*(p6mup3*F(4112)-p1mup3*F(
     -   4210)-p2mup3*F(4223)-p4mup3*F(4233))-p5mup3*F(8600))+2*p4mup
     -   6*(p2mup3*F(4230)+p4mup3*F(4274)+p1mup3*F(8604)))+p4mup3*(p3
     -   mup6*p5mup2*F(8601)-p4mup6*p5mup2*F(8602)-2*(p3mup2*p3mup6*F
     -   (4283)-p3mup2*p4mup6*F(4284)-p3mup6*p4mup2*F(8605)-p4mup2*p4
     -   mup6*F(8606)))
       Fa31(722) = p2mup6*p5mup3*p6mup2*F(4299)+p4mup6*p5mup3*p6mup2*F
     -   (4314)+p3mup6*p5mup3*p6mup2*F(4315)+p5mup3*p5mup6*p6mup2*F(8
     -   318)-p1mup3*p3mup6*p5mup2*F(8531)-p2mup3*p3mup6*p5mup2*F(853
     -   9)-p3mup6*p4mup2*p5mup3*F(8607)+p4mup2*p4mup6*p5mup3*F(8608)
     -   -p1mup3*p4mup6*p5mup2*F(8610)-p3mup2*p3mup6*p5mup3*F(8611)-p
     -   2mup3*p4mup6*p5mup2*F(8612)-p3mup2*p4mup6*p5mup3*F(8614)-p3m
     -   up6*p5mup2*p5mup3*F(8615)-p4mup6*p5mup2*p5mup3*F(8616)-p3mup
     -   3*p4mup6*p5mup2*F(8617)-p3mup3*p3mup6*p5mup2*F(8620)+mup2mup
     -   3*p3mup6*F(8621)+mup2mup3*p4mup6*F(8622)+2*(p2mup6*p4mup3*p6
     -   mup2*F(4300)+p4mup3*p4mup6*p6mup2*F(4307)+p3mup6*p4mup3*p6mu
     -   p2*F(4311)-p1mup3*p3mup6*p4mup2*F(8532)-p2mup3*p3mup6*p4mup2
     -   *F(8540)+p1mup3*p3mup2*p3mup6*F(8555)+p2mup3*p3mup2*p3mup6*F
     -   (8558)+p1mup2*(p3mup3*p3mup6*F(4231)+p3mup3*p4mup6*F(4236)+p
     -   2mup6*p6mup3*F(8566)-p1mup6*p6mup3*F(8570))+p2mup2*(p2mup6*p
     -   6mup3*F(4110)+p3mup3*p3mup6*F(4232)+p3mup3*p4mup6*F(4237)-p1
     -   mup6*p6mup3*F(8571))-p1mup3*p4mup2*p4mup6*F(8609)-p2mup3*p4m
     -   up2*p4mup6*F(8613)+p3mup3*p3mup6*p4mup2*F(8618)-p3mup3*p4mup
     -   2*p4mup6*F(8619)+p1mup3*p3mup2*p4mup6*F(8623)+p3mup2*p3mup3*
     -   p3mup6*F(8624)+p2mup3*p3mup2*p4mup6*F(8625)+p3mup2*p3mup3*p4
     -   mup6*F(8626))
       Fa31(723) = p6mup3*(-(mup2mup6*F(8627))-4*(p3mup6*(2*(p4mup2*F(
     -   5732)-p3mup2*F(5893))+p5mup2*F(8520))-2*(p4mup2*(p4mup6*F(57
     -   34)+p5mup6*F(5737)+p2mup6*F(8512)+p1mup6*F(8573))+p3mup2*(p4
     -   mup6*F(5895)+p5mup6*F(5896)+p2mup6*F(8536)+p1mup6*F(8575)))-
     -   p5mup2*(p2mup6*F(8514)+p1mup6*F(8574)+p4mup6*F(8628)+p5mup6*
     -   F(8629))))
       Fa31(724) = p6mup3*p6mup6*(p5mup2*F(8632)+2*(p4mup2*F(8631)-p3m
     -   up2*F(8633)))
       Fa31(725) = mup2mup6*p5mup3*F(8636)+4*(p3mup2*p3mup3*p6mup6*F(4
     -   324)+p3mup3*p4mup2*p6mup6*F(5732)+mup2mup6*p3mup3*F(7275)-p3
     -   mup2*p4mup3*p6mup6*F(8519)+mup2mup6*p4mup3*F(8635)-mup3mup6*
     -   p4mup2*F(8637)+mup3mup6*p5mup2*F(8638)-mup3mup6*p3mup2*F(864
     -   0)-p2mup2*(p2mup3*(4*p5mup6+p6mup6)*F(4110)-p3mup3*(4*p5mup6
     -   +p6mup6)*F(4112)-p4mup3*(4*p5mup6+p6mup6)*F(8511)+mup3mup6*F
     -   (8634)+p5mup3*(-4*(p2mup6*F(8343)+p1mup6*F(8389)-p5mup6*F(86
     -   41))+p6mup6*F(8642)))+p4mup2*p5mup3*p6mup6*F(8645)-p4mup3*p5
     -   mup2*p6mup6*F(8646)+p3mup2*p5mup3*p6mup6*F(8648)-p5mup2*p5mu
     -   p3*p6mup6*F(8652)-mup2mup3*p6mup6*F(8654)+p3mup3*p5mup2*p6mu
     -   p6*F(8655)+p4mup2*p4mup3*p6mup6*F(8662)+p2mup3*(p3mup2*(4*p5
     -   mup6+p6mup6)*F(4322)-p4mup2*(4*p5mup6+p6mup6)*F(8512)+mup2mu
     -   p6*F(8639)+p5mup2*(p6mup6*F(8649)+4*(p5mup6*F(8659)-p1mup6*F
     -   (9628)-p2mup6*F(9635)))))+16*(p3mup2*p3mup3*p5mup6*F(4324)+p
     -   3mup3*p4mup2*p5mup6*F(5732)-p5mup2*p5mup3*p5mup6*F(8518)-p3m
     -   up2*p4mup3*p5mup6*F(8519)+p4mup2*p5mup3*p5mup6*F(8643)-p4mup
     -   3*p5mup2*p5mup6*F(8644)-p3mup2*p5mup3*p5mup6*F(8647)+p2mup6*
     -   p4mup2*p5mup3*F(8650)+p2mup6*p4mup3*p5mup2*F(8651)+mup2mup3*
     -   p5mup6*F(8653)+p1mup6*p4mup2*p5mup3*F(8656)+p1mup6*p4mup3*p5
     -   mup2*F(8657)+p1mup6*p5mup2*p5mup3*F(8658)+p2mup6*p5mup2*p5mu
     -   p3*F(8660)+p3mup3*p5mup2*p5mup6*F(8661)+p4mup2*p4mup3*p5mup6
     -   *F(8662)+p1mup6*p3mup2*p5mup3*F(9628)-p1mup6*p3mup3*p5mup2*F
     -   (9630)+p2mup6*p3mup2*p5mup3*F(9635)-p2mup6*p3mup3*p5mup2*F(9
     -   637)+mup2mup3*p1mup6*F(9650)+mup2mup3*p2mup6*F(9651))
       Fa31(726) = p1mup2*(p3mup3*(4*p5mup6+p6mup6)*F(4480)-p2mup3*(4*
     -   p5mup6+p6mup6)*F(8566)+p4mup3*p6mup6*F(8568)+p1mup3*p6mup6*F
     -   (8570)+mup3mup6*F(8663)-p5mup3*p6mup6*F(8666)+4*(p4mup3*p5mu
     -   p6*F(8568)+p1mup3*p5mup6*F(8570)-p5mup3*p5mup6*F(8665)+p2mup
     -   6*p5mup3*F(8667)+p1mup6*p5mup3*F(8668)))+p1mup3*(p2mup2*(4*p
     -   5mup6+p6mup6)*F(8571)-p4mup2*(4*p5mup6+p6mup6)*F(8573)-mup2m
     -   up6*F(8664)+p5mup2*p6mup6*F(8669)+p3mup2*p6mup6*F(8670)+4*(p
     -   3mup2*p5mup6*F(8670)+p5mup2*p5mup6*F(8671)-p2mup6*p5mup2*F(9
     -   628)-p1mup6*p5mup2*F(9656)))
       Fa31(727) = p1mup2*p6mup3*(4*p5mup6+p6mup6)*F(8589)+p2mup2*p6mu
     -   p3*(4*p5mup6+p6mup6)*F(8590)-4*(p1mup2*p4mup6*p5mup3*F(8673)
     -   +p2mup2*p4mup6*p5mup3*F(8674))
       Fa31(728) = -4*p4mup3*p6mup2*(p6mup6*F(8596)-4*p5mup6*F(8675))
       Fa31(729) = p2mup2*p3mup6*p5mup3*F(5664)-p3mup6*p4mup3*p5mup2*F
     -   (7109)+p1mup2*p3mup6*p5mup3*F(8676)+p4mup3*p4mup6*p5mup2*F(8
     -   677)
       Fa31(730) = -(p4mup2*p4mup6*p5mup3*F(8679))-p4mup6*p5mup2*p5mup
     -   3*F(8681)+p2mup6*p5mup3*p6mup2*F(9618)+p4mup6*p5mup3*p6mup2*
     -   F(9620)+p5mup3*p5mup6*p6mup2*F(9621)+p1mup3*p4mup6*p5mup2*F(
     -   9631)+p2mup3*p4mup6*p5mup2*F(9638)-p3mup2*p4mup6*p5mup3*F(96
     -   38)+p3mup3*p4mup6*p5mup2*F(9642)+p3mup6*(-(p4mup2*p5mup3*F(8
     -   678))+p5mup3*p6mup2*F(9624)-p3mup2*p5mup3*F(9637)+p5mup2*(-(
     -   p5mup3*F(8680))+p1mup3*F(9630)+p2mup3*F(9637)+p3mup3*F(9643)
     -   )-mup2mup3*F(9652))-mup2mup3*p4mup6*F(9653)
       Fa31(731) = p6mup3*(-(mup2mup6*F(8672))+4*(p4mup2*p5mup6*F(8631
     -   )+p3mup2*p5mup6*F(8682)+p5mup2*(-(p5mup6*F(8683))+p2mup6*F(9
     -   618)+p4mup6*F(9620)+p3mup6*F(9624)+p1mup6*F(9626))))
       Fa31(732) = p6mup3*p6mup6*(-(p4mup2*F(8631))+p3mup2*F(8633)+p5m
     -   up2*F(8684))
       Fa31(733) = p1mup2*p5mup5*F(8689)+p1mup2*p6mup5*F(8690)+p5mup5*
     -   p6mup2*F(8692)+p6mup2*p6mup5*F(8693)+p2mup2*p5mup5*F(8695)+p
     -   2mup2*p6mup5*F(8696)+p1mup2*p2mup5*F(8698)+p2mup2*p4mup5*F(8
     -   699)+p4mup2*p5mup5*F(8702)+p4mup2*p6mup5*F(8703)+p1mup2*p4mu
     -   p5*F(8706)+p5mup2*(p5mup5*F(8686)+p6mup5*F(8687)+p2mup5*F(86
     -   97)+p4mup5*F(8704)+p3mup5*F(8707))+p4mup2*p4mup5*F(8708)+p1m
     -   up2*p3mup5*F(8709)-4*(p2mup5*p6mup2*F(8700)-p4mup5*p6mup2*F(
     -   8705)+p3mup5*p6mup2*F(8710))+p2mup2*p2mup5*F(8711)+p1mup5*(p
     -   5mup2*F(8685)+p1mup2*F(8688)+4*p6mup2*F(8691)+p2mup2*F(8694)
     -   +p4mup2*F(8701)+2*p3mup2*F(8713))+p3mup2*p5mup5*F(8714)+p3mu
     -   p2*p6mup5*F(8715)+p3mup5*p4mup2*F(8716)+p2mup2*p3mup5*F(8717
     -   )+p2mup5*p4mup2*F(8718)-2*(mup2mup5*F(8712)-p3mup2*p4mup5*F(
     -   8719)-p3mup2*p3mup5*F(8720)-p2mup5*p3mup2*F(8721))
       Fa31(734) = mup3mup5*(-(p6mup2*F(8722))+p1mup2*F(9650)+p2mup2*F
     -   (9651)+p3mup2*F(9652)+p4mup2*F(9653)+p5mup2*F(9654))
       Fa31(735) = p4mup5*p6mup3*(p2mup2*F(4114)+p1mup2*F(8723))
       Fa31(736) = p4mup3*(p5mup5*(-(p2mup2*F(8511))+p3mup2*F(8519)-p4
     -   mup2*F(8662)+p5mup2*F(8724))-4*mup2mup5*F(9653))
       Fa31(737) = p3mup3*p5mup5*(2*p6mup2*F(4060)+p2mup2*F(4112)+p1mu
     -   p2*F(4480)+p5mup2*F(5680))+p1mup2*p4mup3*p5mup5*F(8568)+p1mu
     -   p2*p1mup3*p5mup5*F(8570)+p1mup3*p2mup2*p5mup5*F(8571)+p1mup2
     -   *p5mup5*p6mup3*F(8725)+p2mup2*p5mup5*p6mup3*F(8726)+p5mup3*p
     -   5mup5*p6mup2*F(8727)+p4mup3*p5mup5*p6mup2*F(8728)+p2mup2*p5m
     -   up3*p5mup5*F(8730)+p1mup2*p5mup3*p5mup5*F(8732)-p2mup3*(p5mu
     -   p5*(2*p6mup2*F(4044)+p2mup2*F(4110)+p1mup2*F(8566)+p5mup2*F(
     -   8731))+4*p4mup5*p6mup2*F(8735))+p5mup2*p5mup3*p5mup5*F(8736)
     -   -4*(p2mup2*p3mup5*p6mup3*F(4111)+p2mup5*p4mup3*p6mup2*F(7521
     -   )+p2mup5*p5mup3*p6mup2*F(7554)-p3mup5*p5mup3*p6mup2*F(7569)+
     -   p1mup2*p3mup5*p6mup3*F(8729)+p4mup3*p4mup5*p6mup2*F(8733)+p4
     -   mup5*p5mup3*p6mup2*F(8734)+p3mup3*p3mup5*p6mup2*F(8737))
       Fa31(738) = -(p5mup5*(p2mup3*p3mup2*F(4322)+p3mup2*p3mup3*F(432
     -   4)+p3mup3*p4mup2*F(5732)-p2mup3*p4mup2*F(8512)-p1mup3*p4mup2
     -   *F(8573)+p1mup3*p3mup2*F(8670)+p4mup2*p5mup3*F(8739)-p3mup2*
     -   p5mup3*F(8740)+mup2mup3*F(8741)-p1mup3*p5mup2*F(8743)))-4*(p
     -   2mup2*p2mup5*p6mup3*F(4109)-p1mup5*p2mup2*p6mup3*F(4122)+p3m
     -   up5*p4mup3*p6mup2*F(7564)+p1mup2*p2mup5*p6mup3*F(8738)-p1mup
     -   2*p1mup5*p6mup3*F(8742)-p2mup3*p2mup5*p6mup2*F(8744)-p2mup5*
     -   p3mup3*p6mup2*F(8745)-p2mup3*p3mup5*p6mup2*F(8745)-p3mup3*p4
     -   mup5*p6mup2*F(8746)+mup2mup5*p1mup3*F(9650)-mup2mup3*p1mup5*
     -   F(9650)+mup2mup5*p2mup3*F(9651)-mup2mup3*p2mup5*F(9651)+mup2
     -   mup5*p3mup3*F(9652)-mup2mup3*p3mup5*F(9652)+mup2mup5*p5mup3*
     -   F(9654))
       Fa31(739) = p6mup5*(-(p3mup3*p5mup2*F(5680))-p3mup3*p4mup2*F(57
     -   32)+p3mup2*p3mup3*F(5893)+p2mup3*p4mup2*F(8512)+p3mup2*p4mup
     -   3*F(8519)+p2mup3*p3mup2*F(8536)+p4mup2*p4mup3*F(8546)+p1mup3
     -   *p4mup2*F(8573)+p1mup3*p3mup2*F(8575)+p4mup3*p5mup2*F(8724)+
     -   p2mup2*(p2mup3*F(4110)-p3mup3*F(4112)-p4mup3*F(8511)-p1mup3*
     -   F(8571)-p5mup3*F(8730))+p2mup3*p5mup2*F(8731)-p1mup2*(p3mup3
     -   *F(4480)-p2mup3*F(8566)+p4mup3*F(8568)+p1mup3*F(8570)+p5mup3
     -   *F(8732))-p4mup2*p5mup3*F(8739)+p3mup2*p5mup3*F(8740)+p1mup3
     -   *p5mup2*F(8743)-mup2mup3*F(8749)+p5mup2*p5mup3*F(8750))
       Fa31(740) = 4*(p2mup2*p6mup3*p6mup5*F(4299)+p2mup3*p6mup2*p6mup
     -   5*F(4301))+p6mup3*(p5mup2*p5mup5*F(5646)+p4mup2*p5mup5*F(572
     -   5)-p3mup2*p5mup5*F(5870)-4*(p2mup5*p5mup2*F(8194)-p1mup2*p6m
     -   up5*F(8630)))
       Fa31(741) = 2*mup2mup5*p6mup3*F(8751)+p6mup2*(p1mup3*(2*p2mup5*
     -   F(8748)+p5mup5*F(8752))+2*(p4mup3*p6mup5*F(4310)+p3mup3*p6mu
     -   p5*F(4316)+p5mup3*p6mup5*F(7552)+p1mup5*p5mup3*F(8753)+p1mup
     -   5*p3mup3*F(8754)+p1mup3*p3mup5*F(8754)+p1mup3*p4mup5*F(8755)
     -   +p1mup3*p1mup5*F(8756)))
       Fa31(742) = p6mup3*(p3mup5*(p5mup2*F(5681)+p4mup2*F(5733)+p3mup
     -   2*F(6075))+p2mup5*p3mup2*F(8176)-p2mup5*p4mup2*F(8757)-p4mup
     -   2*p4mup5*F(8758)-p4mup5*p5mup2*F(8759)-p1mup5*p4mup2*F(8760)
     -   -p1mup5*p5mup2*F(8761)+p1mup5*p3mup2*F(8762)+p3mup2*p4mup5*F
     -   (8763))
       Fa31(743) = p6mup5*(p4mup2*p6mup3*F(4314)+p3mup2*p6mup3*F(4315)
     -   +p5mup2*p6mup3*F(8318)+p1mup3*p6mup2*F(8764))
       Fa31(744) = p4mup5*(mup2mup3*F(8622)+2*(p2mup2*(p2mup3*F(4230)+
     -   p3mup3*F(4237)+p4mup3*F(4274)-p5mup3*F(8765))+p4mup2*(p4mup3
     -   *F(8606)-p2mup3*F(8613)-p3mup3*F(8619)+p5mup3*F(8766))+p3mup
     -   2*(p4mup3*F(4284)+p2mup3*F(8625)+p3mup3*F(8626)-p5mup3*F(876
     -   8))+p5mup2*(-(p4mup3*F(8767))-p5mup3*F(8769)+p2mup3*F(8770)+
     -   p3mup3*F(8771))))
       Return
       End
