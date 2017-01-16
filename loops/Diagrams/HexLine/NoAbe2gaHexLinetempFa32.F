c************** Calling the Fa functions*************************
       subroutine NoAbe2gaHexLineFa32(p1mup2,p1mup3,p1mup5,p1mup6,p2mu
     -   p2,p2mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p
     -   4mup3,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup
     -   3,p6mup5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6
     -   ,mup5mup6,Fa32)
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
       Complex*16 Fa32(745:768)
       COMMON/NoAbe2gaHexLineFaFunctions/Fa
       Fa32(745) = p4mup5*(p1mup2*(p4mup3*F(4219)+p2mup3*F(4229)+p3mup
     -   3*F(4236)+p1mup3*F(8603)-p5mup3*F(8772))+p1mup3*(p2mup2*F(86
     -   04)-p4mup2*F(8609)+p3mup2*F(8623)+p5mup2*F(8773)))
       Fa32(746) = mup3mup5*(p5mup2*F(8774)+p6mup2*F(8775)-2*(p1mup2*F
     -   (8776)+p2mup2*F(8777)+p3mup2*F(8778)+p4mup2*F(8779)))
       Fa32(747) = p4mup5*p6mup3*(p1mup2*F(8780)+p2mup2*F(8781))
       Fa32(748) = p4mup3*(p5mup2*p5mup5*F(8782)-4*(p2mup2*p5mup5*F(87
     -   83)-p3mup2*p5mup5*F(8785)+p4mup2*p5mup5*F(8787))+8*(p1mup5*p
     -   2mup2*F(4212)+p1mup2*(p2mup5*F(4211)+p1mup5*F(4215))+p2mup2*
     -   p2mup5*F(4224)+p2mup2*p3mup5*F(4233)+p2mup5*p3mup2*F(4234)-p
     -   3mup2*p3mup5*F(4283)+p2mup5*p4mup2*F(8560)+p3mup5*p4mup2*F(8
     -   605)-p2mup5*p5mup2*F(8784)+p3mup5*p5mup2*F(8786)+mup2mup5*F(
     -   8788)))
       Fa32(749) = -(p5mup3*p5mup5*p6mup2*F(8793))+2*(p6mup2*(p2mup3*p
     -   5mup5*F(4408)+p3mup3*p5mup5*F(4426)+2*p4mup5*p5mup3*F(7231))
     -   +p4mup3*p5mup5*p6mup2*F(8794))+p2mup3*p5mup2*p5mup5*F(8802)+
     -   p5mup2*p5mup3*p5mup5*F(8805)+p3mup3*p5mup2*p5mup5*F(8806)+p1
     -   mup2*(8*p2mup3*p2mup5*F(4221)+p3mup5*(-8*(p5mup3*F(4184)-p1m
     -   up3*F(4209)-p4mup3*F(4217)-p2mup3*F(4222)-p3mup3*F(4231))-4*
     -   p6mup3*F(8795))+p5mup3*(p5mup5*F(8789)+8*p2mup5*F(8803))+2*p
     -   5mup5*(p6mup3*F(8791)-2*(-(p2mup3*F(8798))+p3mup3*F(8800)+p4
     -   mup3*F(8804)+p1mup3*F(8807))))+p2mup2*(p5mup3*p5mup5*F(8790)
     -   +2*p5mup5*p6mup3*F(8792)+8*(p1mup3*p3mup5*F(4210)+p2mup3*p3m
     -   up5*F(4223)+p2mup5*(p5mup3*F(4187)+p2mup3*F(4220)+p3mup3*F(4
     -   223))+p3mup3*p3mup5*F(4232)-p3mup5*p5mup3*F(8796))-4*(p3mup5
     -   *p6mup3*F(8797)-p2mup3*p5mup5*F(8799)+p3mup3*p5mup5*F(8801)+
     -   p1mup3*p5mup5*F(8808)))+8*(p1mup5*p3mup2*p4mup3*F(4218)+p1mu
     -   p5*p4mup2*p4mup3*F(8556)+p2mup5*p3mup3*p5mup2*F(8559)+p2mup3
     -   *p3mup5*p5mup2*F(8559)+p3mup3*p3mup5*p5mup2*F(8561)-p1mup5*p
     -   4mup3*p5mup2*F(8809)-p2mup5*p5mup2*p5mup3*F(8810)-p3mup5*p5m
     -   up2*p5mup3*F(8811))-4*(p2mup5*p4mup3*p6mup2*F(7702)-p2mup5*p
     -   5mup3*p6mup2*F(7711)+p4mup3*p4mup5*p6mup2*F(7714)-p3mup5*p5m
     -   up3*p6mup2*F(7716)+p2mup3*p4mup5*p6mup2*F(9638)+p3mup3*p3mup
     -   5*p6mup2*F(9643))
       Fa32(750) = -(p4mup2*p5mup3*p5mup5*F(8812))+p3mup2*p5mup3*p5mup
     -   5*F(8815)-2*mup2mup3*p5mup5*F(8821)+p1mup3*p5mup2*p5mup5*F(8
     -   828)+8*(p1mup3*p2mup2*p2mup5*F(4208)-p3mup2*p3mup5*p5mup3*F(
     -   4272)-p1mup3*p2mup5*p4mup2*F(8528)-p1mup3*p3mup5*p4mup2*F(85
     -   32)-p2mup3*p2mup5*p4mup2*F(8534)-p2mup5*p3mup3*p4mup2*F(8540
     -   )-p2mup3*p3mup5*p4mup2*F(8540)+p1mup3*p2mup5*p3mup2*F(8554)+
     -   p1mup3*p3mup2*p3mup5*F(8555)+p2mup3*p2mup5*p3mup2*F(8557)+p2
     -   mup5*p3mup2*p3mup3*F(8558)+p2mup3*p3mup2*p3mup5*F(8558)+p1mu
     -   p3*p3mup5*p5mup2*F(8588)+p3mup3*p3mup5*p4mup2*F(8618)+p3mup2
     -   *p3mup3*p3mup5*F(8624)+p2mup5*p4mup2*p5mup3*F(8814)-p3mup5*p
     -   4mup2*p5mup3*F(8819)-p2mup5*p3mup2*p5mup3*F(8833)+mup2mup5*p
     -   1mup3*F(8835)+mup2mup5*p2mup3*F(8836)+mup2mup5*p3mup3*F(8837
     -   )+p1mup3*p2mup5*p5mup2*F(8839)+p2mup3*p2mup5*p5mup2*F(8840))
     -   +4*(p3mup5*p4mup3*p6mup2*F(7703)+mup2mup3*p2mup5*F(8552)+mup
     -   2mup3*p3mup5*F(8621)+p2mup3*p4mup2*p5mup5*F(8813)+p2mup2*p2m
     -   up5*p6mup3*F(8817)-p3mup3*p4mup2*p5mup5*F(8818)-mup2mup5*p5m
     -   up3*F(8820)+p1mup2*(p2mup5*(2*(p1mup3*F(4207)+p3mup3*F(4222)
     -   )+p6mup3*F(8816))+p1mup5*(2*(p2mup3*F(4207)+p3mup3*F(4209)+p
     -   1mup3*F(8583)-p5mup3*F(8822))-p6mup3*F(8824)))+p1mup3*p4mup2
     -   *p5mup5*F(8826)-p1mup3*p3mup2*p5mup5*F(8829)-p2mup3*p3mup2*p
     -   5mup5*F(8832)-p3mup2*p3mup3*p5mup5*F(8834)+p1mup5*(mup2mup3*
     -   F(8551)+p2mup2*(2*(p2mup3*F(4208)+p3mup3*F(4210)+p1mup3*F(85
     -   84)-p5mup3*F(8823))-p6mup3*F(8825))+2*(-(p2mup3*p4mup2*F(852
     -   8))-p3mup3*p4mup2*F(8532)+p2mup3*p3mup2*F(8554)+p3mup2*p3mup
     -   3*F(8555)-p1mup3*p4mup2*F(8585)+p1mup3*p3mup2*F(8587)+p3mup3
     -   *p5mup2*F(8588)+p4mup2*p5mup3*F(8827)-p3mup2*p5mup3*F(8830)-
     -   p5mup2*p5mup3*F(8831)+p1mup3*p5mup2*F(8838)+p2mup3*p5mup2*F(
     -   8839)))-p2mup3*p2mup5*p6mup2*F(9635)-p2mup5*p3mup3*p6mup2*F(
     -   9637)-p2mup3*p3mup5*p6mup2*F(9637)-p3mup3*p4mup5*p6mup2*F(96
     -   42))
       Fa32(751) = p6mup5*(p5mup2*p5mup3*F(8518)-mup2mup3*F(8523)-2*(p
     -   3mup3*p4mup2*F(5732)-p3mup2*p3mup3*F(5893)-p2mup3*p4mup2*F(8
     -   512)-p3mup2*p4mup3*F(8519)-p2mup3*p3mup2*F(8536)-p4mup2*p4mu
     -   p3*F(8546)+p1mup2*(p3mup3*F(4480)-p2mup3*F(8566)+p4mup3*F(85
     -   68)+p1mup3*F(8570))-p1mup3*p4mup2*F(8573)-p1mup3*p3mup2*F(85
     -   75))+p1mup2*p5mup3*F(8841)+p2mup2*(2*(p2mup3*F(4110)-p3mup3*
     -   F(4112)-p4mup3*F(8511)-p1mup3*F(8571))+p5mup3*F(8842))-p4mup
     -   2*p5mup3*F(8843)+p4mup3*p5mup2*F(8844)+p2mup3*p5mup2*F(8845)
     -   +p3mup2*p5mup3*F(8846)+p3mup3*p5mup2*F(8847)+p1mup3*p5mup2*F
     -   (8848))
       Fa32(752) = p6mup3*(p5mup2*p5mup5*F(8852)-2*(p4mup2*p5mup5*F(88
     -   51)-p3mup2*p5mup5*F(8853))+4*(p1mup2*p6mup5*F(8849)+p2mup2*p
     -   6mup5*F(8850)-p2mup5*p5mup2*F(8854)))+4*p2mup3*p6mup2*p6mup5
     -   *F(9618)
       Fa32(753) = -2*mup2mup5*p6mup3*F(8855)+p6mup2*(p1mup3*(p5mup5*F
     -   (8856)-2*(p2mup5*F(9628)+p3mup5*F(9630)+p4mup5*F(9631)))+2*(
     -   p6mup5*(p4mup3*F(7733)+p5mup3*F(9621)+p3mup3*F(9624))+p1mup5
     -   *(p5mup3*F(7718)-p3mup3*F(9630)-p1mup3*F(9656))))
       Fa32(754) = p6mup3*(-(p2mup5*p4mup2*F(8857))-p4mup2*p4mup5*F(88
     -   58)+p3mup5*p4mup2*F(8859)+p4mup5*p5mup2*F(8860)+p3mup5*p5mup
     -   2*F(8861)-p1mup5*p4mup2*F(8862)+p1mup5*p5mup2*F(8863)+p1mup5
     -   *p3mup2*F(8864)+p2mup5*p3mup2*F(8865)+p3mup2*p3mup5*F(8866)+
     -   p3mup2*p4mup5*F(8867))
       Fa32(755) = p6mup5*(p5mup2*p6mup3*F(8683)-p4mup2*p6mup3*F(8868)
     -   +p3mup2*p6mup3*F(8869)+p1mup3*p6mup2*F(9626))
       Fa32(756) = p2mup2*F(8870)+p4mup2*F(8871)+p5mup2*F(8872)+2*(p6m
     -   up2*F(7788)+p1mup2*F(8873)-p3mup2*F(8874))
       Fa32(757) = p6mup2*F(7794)-p2mup2*F(8875)-p4mup2*F(8876)+p5mup2
     -   *F(8877)+p1mup2*F(8878)-p3mup2*F(8879)
       Fa32(758) = p1mup2*p5mup5*F(8883)+p1mup2*p6mup5*F(8884)+p5mup5*
     -   p6mup2*F(8886)+2*p4mup2*p5mup5*F(8889)+p4mup2*p6mup5*F(8890)
     -   +p5mup2*p5mup5*F(8891)+p5mup2*p6mup5*F(8892)-4*(p2mup5*p6mup
     -   2*F(8888)-p6mup2*p6mup5*F(8893)+p3mup5*p6mup2*F(8895)+p4mup5
     -   *p6mup2*F(8896))+p2mup2*(p5mup5*F(8881)+p6mup5*F(8882)+p2mup
     -   5*F(8887)+p3mup5*F(8894)+p4mup5*F(8898))+p4mup5*p5mup2*F(890
     -   0)+p1mup2*p4mup5*F(8901)+p4mup2*p4mup5*F(8902)+mup2mup5*F(89
     -   03)+p3mup2*p5mup5*F(8905)+p3mup2*p6mup5*F(8906)+p1mup5*(p2mu
     -   p2*F(8880)-4*p6mup2*F(8885)+p1mup2*F(8897)+p4mup2*F(8899)+p3
     -   mup2*F(8904)+p5mup2*F(8907))+p1mup2*p2mup5*F(8908)+p3mup2*p4
     -   mup5*F(8909)+p1mup2*p3mup5*F(8910)+p3mup5*p5mup2*F(8911)+p3m
     -   up2*p3mup5*F(8912)+p2mup5*p5mup2*F(8913)+p2mup5*p4mup2*F(891
     -   4)+p2mup5*p3mup2*F(8915)+p3mup5*p4mup2*F(8916)
       Fa32(759) = p5mup2*p5mup5*F(8921)+p5mup2*p6mup5*F(8922)+p1mup2*
     -   p5mup5*F(8924)+p1mup2*p6mup5*F(8925)+p1mup2*p2mup5*F(8931)+p
     -   4mup2*p5mup5*F(8934)+p1mup2*p3mup5*F(8938)+4*(p6mup2*p6mup5*
     -   F(8928)+p2mup5*p6mup2*F(8932)+p3mup5*p6mup2*F(8939)+p4mup5*p
     -   6mup2*F(8941))+p1mup2*p4mup5*F(8942)+p2mup2*(p5mup5*F(8918)+
     -   p6mup5*F(8919)+p2mup5*F(8929)+p3mup5*F(8937)+p4mup5*F(8943))
     -   +mup2mup5*F(8944)+p1mup5*(p2mup2*F(8917)+p1mup2*F(8923)-4*p6
     -   mup2*F(8926)+2*(p5mup2*F(8920)+p4mup2*F(8933)+p3mup2*F(8945)
     -   ))+p3mup2*p6mup5*F(8947)+p4mup2*p4mup5*F(8949)-2*(p5mup5*p6m
     -   up2*F(8927)-p2mup5*p5mup2*F(8930)-p4mup2*p6mup5*F(8935)-p3mu
     -   p5*p5mup2*F(8936)-p4mup5*p5mup2*F(8940)-p3mup2*p5mup5*F(8946
     -   )+p3mup5*p4mup2*F(8948)+p2mup5*p4mup2*F(8950)+p3mup2*p4mup5*
     -   F(8951)+p3mup2*p3mup5*F(8952)+p2mup5*p3mup2*F(8953))
       Fa32(760) = mup5mup6*(p4mup2*F(8955)+p5mup2*F(8956)-4*(p2mup2*F
     -   (8954)-p3mup2*F(8957)))
       Fa32(761) = 4*(p2mup2*p2mup6*p5mup5*F(8960)-p1mup6*p4mup2*p5mup
     -   5*F(8962)-p2mup6*p4mup2*p5mup5*F(8964))+p4mup2*p5mup5*p6mup6
     -   *F(8980)+8*(p2mup5*p4mup2*p6mup6*F(8963)+p1mup5*(-2*(p5mup2*
     -   p5mup6*F(8976)-p2mup2*(p2mup6*F(8959)+p5mup6*F(8965)+p1mup6*
     -   F(8981)))+p4mup2*(p6mup6*F(8961)-2*(p2mup6*F(8966)+p5mup6*F(
     -   8970)-p1mup6*F(8985))))+p4mup2*p5mup5*p5mup6*F(8990)+p5mup2*
     -   p5mup5*p5mup6*F(8991)-p4mup2*p6mup5*p6mup6*F(9622))+p3mup2*(
     -   p5mup5*(8*p5mup6*F(5539)+4*(p2mup6*F(8968)+p1mup6*F(8971))+p
     -   6mup6*F(8975))-8*(p6mup6*(p2mup5*F(7711)+p1mup5*F(7718))+p6m
     -   up5*p6mup6*F(9621)-2*(p1mup6*(p6mup5*F(8972)+p1mup5*F(8982)+
     -   p2mup5*F(8986))+p2mup6*(p6mup5*F(8974)+p1mup5*F(8983)+p2mup5
     -   *F(8987))+p5mup6*(p1mup5*F(8984)+p2mup5*F(8988)+p6mup5*F(964
     -   7)))))-16*(p2mup2*p2mup5*p2mup6*F(5812)-p1mup6*p2mup5*p4mup2
     -   *F(8967)-p2mup2*p2mup6*p6mup5*F(8969)+p2mup5*p4mup2*p5mup6*F
     -   (8973)+p2mup5*p5mup2*p5mup6*F(8977)-p1mup6*p4mup2*p6mup5*F(8
     -   978)-p2mup6*p4mup2*p6mup5*F(8979)-p2mup5*p2mup6*p4mup2*F(898
     -   9)-p4mup2*p5mup6*p6mup5*F(9649)-p5mup2*p5mup6*p6mup5*F(9655)
     -   )
       Fa32(762) = mup2mup6*(-2*p5mup5*F(8993)+8*(p1mup5*F(8992)+p2mup
     -   5*F(8995)+p6mup5*F(8996)))-4*(mup5mup6*p6mup2*F(8994)-p1mup6
     -   *p2mup2*p5mup5*F(8998)+p2mup6*p5mup2*p5mup5*F(9001)+p1mup6*p
     -   5mup2*p5mup5*F(9003)+p2mup2*p5mup5*p6mup6*F(9006))+p5mup2*p5
     -   mup5*p6mup6*F(9015)-16*(p1mup6*p2mup2*p2mup5*F(8999)+p1mup6*
     -   p2mup2*p6mup5*F(9005)-p1mup6*p2mup5*p5mup2*F(9009)-p2mup2*p2
     -   mup5*p5mup6*F(9010)-p2mup2*p5mup6*p6mup5*F(9012)-p2mup6*p5mu
     -   p2*p6mup5*F(9013)-p1mup6*p5mup2*p6mup5*F(9014)-p2mup5*p2mup6
     -   *p5mup2*F(9017))-8*(p2mup2*p2mup5*p6mup6*F(9000)-p2mup5*p5mu
     -   p2*p6mup6*F(9004)+p2mup2*p6mup5*p6mup6*F(9007)-p2mup2*p5mup5
     -   *p5mup6*F(9011)+p1mup5*(p2mup2*p6mup6*F(8997)-p5mup2*(p6mup6
     -   *F(9002)-2*(p2mup6*F(9008)-p1mup6*F(9016))))+p5mup2*p6mup5*p
     -   6mup6*F(9625))
       Fa32(763) = 2*mup2mup6*p4mup5*F(9018)+p2mup2*p4mup5*(p6mup6*F(9
     -   028)+4*(p5mup6*F(9026)-p4mup6*F(9036)))+p1mup2*(p2mup6*p5mup
     -   5*F(9024)+p4mup5*p6mup6*F(9027)-p5mup5*p6mup6*F(9031)-2*(p1m
     -   up5*p6mup6*F(9020)-p2mup5*p6mup6*F(9023)+p6mup5*p6mup6*F(903
     -   2)-p5mup5*p5mup6*F(9037))+p1mup6*(p5mup5*F(9019)-4*(p2mup5*F
     -   (9021)-p6mup5*F(9029)-p1mup5*F(9039)))+4*(p1mup5*p2mup6*F(90
     -   22)+p4mup5*p5mup6*F(9025)+p1mup5*p5mup6*F(9030)+p2mup5*p5mup
     -   6*F(9033)+p2mup6*p6mup5*F(9034)-p4mup5*p4mup6*F(9035)+p5mup6
     -   *p6mup5*F(9038)+p2mup5*p2mup6*F(9040)))
       Fa32(764) = 4*mup2mup5*p6mup6*F(9042)+p3mup2*p4mup5*p6mup6*F(90
     -   61)+p4mup2*p4mup5*p6mup6*F(9066)+p4mup5*p5mup2*p6mup6*F(9068
     -   )-8*(mup2mup5*p2mup6*F(9043)-mup2mup6*p3mup5*F(9044)+mup2mup
     -   5*p5mup6*F(9045)-p2mup2*p3mup5*p6mup6*F(9064)-p1mup6*(-(mup2
     -   mup5*F(9041))+2*(p2mup5*p6mup2*F(4304)-p2mup2*p3mup5*F(9046)
     -   -p1mup2*p3mup5*F(9047)-p2mup2*p4mup5*F(9048)-p1mup2*p4mup5*F
     -   (9049)+p4mup2*p4mup5*F(9055)+p4mup5*p5mup2*F(9058)-p3mup5*p4
     -   mup2*F(9060)-p3mup5*p5mup2*F(9062)-p3mup2*p3mup5*F(9069)-p3m
     -   up2*p4mup5*F(9070)))-p1mup2*(p3mup5*p6mup6*F(5877)-2*(p2mup6
     -   *p3mup5*F(9051)+p2mup6*p4mup5*F(9053)+p3mup5*p5mup6*F(9074))
     -   )-p2mup5*p6mup2*p6mup6*F(9616)-p3mup2*p3mup5*p6mup6*F(9644)-
     -   p3mup5*p4mup2*p6mup6*F(9646)-p3mup5*p5mup2*p6mup6*F(9647))-1
     -   6*(p2mup2*p2mup6*p3mup5*F(9050)+p2mup2*p2mup6*p4mup5*F(9052)
     -   -p2mup6*p4mup2*p4mup5*F(9054)+p2mup6*p3mup5*p4mup2*F(9056)-p
     -   2mup6*p4mup5*p5mup2*F(9057)-p3mup2*p4mup5*p5mup6*F(9059)+p2m
     -   up6*p3mup5*p5mup2*F(9063)+p4mup2*p4mup5*p5mup6*F(9065)-p4mup
     -   5*p5mup2*p5mup6*F(9067)+p2mup6*p3mup2*p3mup5*F(9071)+p2mup6*
     -   p3mup2*p4mup5*F(9072)+p2mup2*p3mup5*p5mup6*F(9073)+p2mup5*p5
     -   mup6*p6mup2*F(9621)+p3mup2*p3mup5*p5mup6*F(9686)+p3mup5*p4mu
     -   p2*p5mup6*F(9689)+p3mup5*p5mup2*p5mup6*F(9690))
       Fa32(765) = p4mup6*(p1mup2*(-(p5mup5*F(9075))+2*(p1mup5*F(9077)
     -   +p2mup5*F(9079)+p3mup5*F(9081)))+p2mup2*(-(p5mup5*F(9076))+2
     -   *(p1mup5*F(9078)+p2mup5*F(9080)+p3mup5*F(9082))))
       Fa32(766) = p6mup2*(-(p2mup6*p5mup5*F(5723))+4*(p1mup6*p4mup5*F
     -   (9620)+p2mup6*p4mup5*F(9620)+(p1mup6+p2mup6)*p3mup5*F(9624))
     -   )
       Fa32(767) = -(p1mup6*p5mup5*p6mup2*F(5723))-2*(p4mup6*(mup2mup5
     -   *F(9084)+2*(p2mup2*p6mup5*F(9087)+p1mup2*p6mup5*F(9088)+p2mu
     -   p5*p4mup2*F(9092)))-p4mup6*(-(p3mup2*p5mup5*F(9097))+p4mup2*
     -   p5mup5*F(9105)-p5mup2*p5mup5*F(9107)-2*(p2mup5*p5mup2*F(9095
     -   )-p3mup2*p6mup5*F(9099)-p4mup2*p6mup5*F(9106)-p5mup2*p6mup5*
     -   F(9108)-p2mup5*p3mup2*F(9112)))-p1mup5*(-2*p4mup6*(p4mup2*F(
     -   9093)+p5mup2*F(9096)-p3mup2*F(9110))+p6mup2*(p6mup6*F(9616)+
     -   2*(p2mup6*F(9089)-p5mup6*F(9621))))-p3mup6*(p3mup2*p5mup5*F(
     -   5623)+p4mup2*p5mup5*F(5631)-mup2mup5*F(9083)+p1mup2*(p5mup5*
     -   F(5586)+2*(p6mup5*F(5877)+p1mup5*F(9086)+p2mup5*F(9091)-p4mu
     -   p5*F(9100)))+p2mup2*p5mup5*F(9104)+2*(p2mup5*p4mup2*F(9094)+
     -   p1mup5*p4mup2*F(9098)+p2mup2*(p6mup5*F(9064)+p1mup5*F(9085)+
     -   p2mup5*F(9090)-p4mup5*F(9101))+p1mup5*p3mup2*F(9109)+p2mup5*
     -   p3mup2*F(9111)+p3mup2*p6mup5*F(9644)+p4mup2*p6mup5*F(9646))+
     -   p5mup2*(p5mup5*F(5539)+2*(p1mup5*F(9102)+p2mup5*F(9103)+p6mu
     -   p5*F(9647)))))
       Fa32(768) = p3mup6*p4mup2*p4mup5*F(9114)+p4mup2*p4mup5*p4mup6*F
     -   (9115)+(p1mup5+p2mup5)*p4mup6*p6mup2*F(9620)+(p1mup5+p2mup5)
     -   *p3mup6*p6mup2*F(9624)+p3mup5*(-(p4mup2*p4mup6*F(9113))+p3mu
     -   p6*p5mup2*F(9686))
       Return
       End
