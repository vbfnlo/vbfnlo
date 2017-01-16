      subroutine ten_red3_Gram1(p1sq,p2sq,s12,musq,B023,B013,B012,C30,Cij)
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c************************************************************************************
      Implicit none
      integer k,l,jj
      Real*8 p1sq,p2sq,s12,p1p2,p1sq2,p1sq3,p2sq2,p2sq3,s122
      real*8 ZMax,TX1,TX2,det3,detABS,ratio,line
      Complex*16 C30,Cij(30,9)
      Complex*16 S2000000
      Real*8 r10,r21
      Real*8 Inv00,Inv001,Inv002,Inv0000,Inv0011,Inv0021,Inv0022
      Real*8 Inv00001,Inv00002,Inv00111,Inv00211,Inv00221,Inv00222
      Real*8 Inv000000,Inv000011,Inv000021,Inv000022
      Real*8 Inv002222,Inv002221,Inv002211,Inv002111,Inv001111
      Real*8 Inv0000001,Inv0000002
      Real*8 Inv0000222,Inv0000221,Inv0000211,Inv0000111
      Real*8 Inv0022222,Inv0022221,Inv0022211,Inv0022111,Inv0021111,Inv0011111
      Real*8 Inv00000000
      Real*8 Inv00000022,Inv00000021,Inv00000011
      Real*8 Inv00002222,Inv00002221,Inv00002211
      Real*8 Inv00002111,Inv00001111
      Real*8 Inv00222222,Inv00222221,Inv00222211
      Real*8 Inv00222111,Inv00221111,Inv00211111
      Real*8 Inv00111111
      Real*8 Inv000000001,Inv000000002
      Real*8 Inv000000222,Inv000000221,Inv000000211,Inv000000111
      Real*8 Inv000022222,Inv000022221,Inv000022211,Inv000022111
      Real*8 Inv000021111,Inv000011111
      Real*8 Inv002222222,Inv002222221,Inv002222211,Inv002222111
      Real*8 Inv002221111,Inv002211111,Inv002111111,Inv001111111
      Real*8 Z(2,2),IX,I6Z,I10Z,I12Z,I14Z,I18Z,I22Z,I26Z,I30Z
      Real*8 I34Z,I38Z,I42Z,I46Z,ZZ(2,2,2,2),F(5),F1
      Complex*16 B012,B013,B023
      Complex*16 B12(6,11),B13(6,11),B23(6,11)
      Complex*16 S200,S21(2),S2001(2)
      Complex*16 S2h001(2)
      Complex*16 S20000,S200001(2)
      Complex*16 S2h00001(2)
      Complex*16 S2h0000001(2),S20000001(2)
      Complex*16 auxC30,tempC30
      Complex*16 auxC300,tempC300
      Complex*16 auxC30000,tempC30000
      Complex*16 auxC3000000,tempC3000000
      Complex*16 aux1(2),temp1(2)
      Complex*16 aux2(2,2),temp2(2,2)
      Complex*16 aux3(2,2,2),temp3(2,2,2)
      Complex*16 aux001(2),temp001(2)
      Complex*16 aux002(2,2),temp002(2,2)
      Complex*16 aux003(2,2,2),temp003(2,2,2)
      Complex*16 aux00001(2),temp00001(2)
      Complex*16 aux00002(2,2),temp00002(2,2)
      Complex*16 aux00003(2,2,2),temp00003(2,2,2)
      Complex*16 aux0000001(2),temp0000001(2)
      include 'CijGram.inc'             
      real*8 tempjj1,tempjj2,tempkl1,tempkl2,tempkl3,tempkl4
      real*8 IXtemp,Zmaxtemp
      Integer jjtemp,ktemp,ltemp,jjinit,kinit,linit,cont
      real*8 tempjj,tempkl,IXinit,Zmaxinit
      Common/Decide/tempjj,tempkl,IX,ZMax,jjinit
      Save/Decide/
      integer order
      integer it,jt
      real*8 musq
      real*8 accuracyC(0:4,5),AccuracyD(0:5,4)
      real*8 accuracyCR(9,0:4,5)
      Common/Accuracy/AccuracyC,AccuracyD
      integer i1,i2,index(0:4),ac
       logical printmy
       common/mprint/printmy
      Common/musqInv/musqcp
      real*8 musqcp
      musqcp=musq

      index(0)=1
      index(1)=2
      index(2)=4
      index(3)=6
      index(4)=9

      cont=0
!      print*, "order", order
      order=12

      do jt=1,9
         do it=1,30
            Cij(it,jt)=0d0
         enddo
      enddo
      call ten_red2_forGram(p1sq,B012,B12)
      call ten_red2_forGram(s12,B013,B13)
      call ten_red2_forGram(p2sq,B023,B23)


cFC      call ten_red2_forGram_QUAD_D(p1sq,musq,B012,B12)
cFC      call ten_red2_forGram_QUAD_D(s12,musq,B013,B13)
cFC      call ten_red2_forGram_QUAD_D(p2sq,musq,B023,B23)

       Z(1,1)=2*p2sq
       Z(2,1)=p1sq+p2sq-s12
       Z(1,2)=p1sq+p2sq-s12
       Z(2,2)=2*p1sq
       If(abs(Z(1,1)).ge.abs(Z(1,2))) then
        If(abs(Z(1,1)).ge.abs(Z(2,2))) then
         k=1
         l=1
         ZMax=Z(1,1)
        else
         k=2
         l=2
         ZMax=Z(2,2)
        endif
      else
        If(abs(Z(1,2)).ge.abs(Z(2,2))) then
         k=2
         l=1
         ZMax=Z(1,2)
        else
         k=2
         l=2
         ZMax=Z(2,2)
        endif
      endif
 
           
      r10=p1sq
      r21=s12-p1sq
      p1p2=(s12-p1sq-p2sq)/2d0
      p1sq2=p1sq*p1sq
      p1sq3=p1sq*p1sq*p1sq
      p2sq2=p2sq*p2sq
      p2sq3=p2sq2*p2sq
      s122=s12*s12
      Inv00=-1d0/2d0
      Inv002=1d0/6d0
      Inv001=1d0/3d0
      Inv0000=1d0/48d0*(p1sq+p2sq+s12)
      Inv0022=-1d0/12d0
      Inv0021=-1d0/8d0
      Inv0011=-1d0/4d0
      Inv00002=-1d0/240d0*(2d0*(p2sq+s12)+p1sq)
      Inv00001=-1d0/240d0*(3d0*(p1sq+s12)+4d0*p2sq)
      Inv00222=1d0/20d0
      Inv00221=1d0/15d0
      Inv00211=1d0/10d0
      Inv00111=1d0/5d0
      Inv000000=-1d0/1440d0*(p1sq*p1sq+p2sq*p2sq+p1sq*(s12+p2sq)+s12*(s12+p2sq)            )
      Inv000022=1d0/720d0*(3d0*(p2sq+s12)+p1sq)
      Inv000021=1d0/720d0*(2d0*(p1sq+2d0*s12)+5d0*p2sq)
      Inv000011=1d0/360d0*(3d0*(p1sq+s12)+5d0*p2sq)
      Inv002222=-1d0/30d0
      Inv002221=-1d0/24d0
      Inv002211=-1d0/18d0
      Inv002111=-1d0/12d0
      Inv001111=-1d0/6d0
      Inv0000001=(1d0/10080d0)*(p1sq*(4d0*(p1sq+s12)+5d0*p2sq)+s12*(4d0*s12+5d0*p2sq)+6d0*p2sq*p2sq)
      Inv0000002=(1d0/10080d0)*(p1sq*p1sq+2d0*(p2sq+s12)*p1sq+3d0*(p2sq*p2sq+s12*(s12+p2sq)))
      Inv0000222=-(1d0/1680d0)*(p1sq+4d0*(s12+p2sq))
      Inv0000221=-(1d0/5040d0)*(5d0*(p1sq+3d0*s12)+18d0*p2sq)
      Inv0000211=-(1d0/504d0)*(p1sq+2d0*s12+3d0*p2sq)
      Inv0000111=-(1d0/168d0)*(p1sq+s12+2d0*p2sq)
      Inv0022222=(1d0/42d0)
      Inv0022221=(1d0/35d0)
      Inv0022211=(1d0/28d0)
      Inv0022111=(1d0/21d0)
      Inv0021111=(1d0/14d0)
      Inv0011111=(1d0/7d0)
      Inv00000000=(3*p1sq3 + 3*p1sq2*(p2sq + s12) + 3*(p2sq + s12)*(p2sq2 + s122) +p1sq*(3*p2sq2 + 4*p2sq*s12 + 3*s122))/161280d0
      Inv00000022=(-p1sq2 - 3*p1sq*(p2sq + s12) - 6*(p2sq2 + p2sq*s12 + s122))/40320d0
      Inv00000021=(-5*p1sq2 - 12*p1sq*p2sq - 21*p2sq2 - 10*p1sq*s12 - 18*p2sq*s12 - 15*s122)/80640d0
      Inv00000011=(-10*p1sq2 - 15*p1sq*p2sq - 21*p2sq2 - 10*p1sq*s12 - 15*p2sq*s12 - 10*s122)/40320d0
      Inv00002222=(p1sq + 5*(p2sq + s12))/3360d0
      Inv00002221=(3*p1sq + 14*p2sq + 12*s12)/6720d0
      Inv00002211=(5*p1sq + 21*p2sq + 15*s12)/6720d0
      Inv00002111=(2*p1sq + 7*p2sq + 4*s12)/1344d0
      Inv00001111=(3*p1sq + 7*p2sq + 3*s12)/672d0
      Inv00222222=-1d0/56d0
      Inv00222221=-1d0/48d0
      Inv00222211=-1d0/40d0
      Inv00222111=-1d0/32d0
      Inv00221111=-1d0/24d0
      Inv00211111=-1d0/16d0
      Inv00111111=-1d0/8d0
      Inv000000002=(-p1sq3 - 2*p1sq2*(p2sq + s12) - 4*(p2sq + s12)*(p2sq2 + s122) -
     -    p1sq*(3*p2sq2 + 4*p2sq*s12 + 3*s122))/483840d0
      Inv000000001=(-5*p1sq3 - 6*p1sq2*p2sq - 7*p1sq*p2sq2 - 8*p2sq3 - 5*p1sq2*s12 - 8*p1sq*p2sq*s12 -
     -    7*p2sq2*s12 - 5*s12**3 - 5*p1sq*s122 - 6*p2sq*s122)/483840d0
      Inv000000222=(p1sq2 + 4*p1sq*(p2sq + s12) + 10*(p2sq2 + p2sq*s12 + s122))/120960d0
      Inv000000221=(2*p1sq2 + 7*p1sq*p2sq + 16*p2sq2 + 6*p1sq*s12 + 14*p2sq*s12 + 12*s122)/120960d0
      Inv000000211=(5*p1sq2 + 28*p2sq2 + 21*p2sq*s12 + 2*p1sq*(7*p2sq + 5*s12) + 15*s122)/120960d0
      Inv000000111=(20*p1sq2 + 56*p2sq2 + 35*p2sq*s12 + 5*p1sq*(7*p2sq + 4*s12) + 20*s122)/120960d0
      Inv000022222=(-p1sq - 6*(p2sq + s12))/6048d0
      Inv000022221=(-7*p1sq - 5*(8*p2sq + 7*s12))/30240d0
      Inv000022211=(-3*p1sq - 4*(4*p2sq + 3*s12))/8640d0
      Inv000022111=(-5*p1sq - 3*(8*p2sq + 5*s12))/8640d0
      Inv000021111=(-p1sq - 4*p2sq - 2*s12)/864d0
      Inv000011111=(-3*p1sq - 8*p2sq - 3*s12)/864d0
      Inv002222222=1/72d0
      Inv002222221=1/63d0
      Inv002222211=1/54d0
      Inv002222111=1/45d0
      Inv002221111=1/36d0
      Inv002211111=1/27d0
      Inv002111111=1/18d0
      Inv001111111=1/9d0
      Inv100=(2*p1sq**4+2*p1sq**3*(p2sq+s12)+p1sq**2*(2*p2sq**2+3*p2sq*s12+2*s
     &  12**2)+p1sq*(2*p2sq**3+3*p2sq**2*s12+3*p2sq*s12**2+2*s12**3)+2*(
     &  p2sq**4+p2sq**3*s12+p2sq**2*s12**2+p2sq*s12**3+s12**4))/2.4192d6
      Inv8011=(-36*p2sq**3-28*p2sq**2*s12-21*p2sq*s12**2-3*p1sq**2*(7*p2
     &  sq+5*s12)-15*(p1sq**3+s12**3)-p1sq*(15*s12**2+28*p2sq*(p2sq+s12)))
     &  /1.2096d6
      Inv8021=(-3*p1sq**3-p1sq**2*(7*p2sq+6*s12)-p1sq*(12*p2sq**2+14*p2sq*
     &  s12+9*s12**2)-2*(9*p2sq**3+8*p2sq**2*s12+7*p2sq*s12**2+6*s12**3)
     &  )/1.2096d6
      Inv8022=(-p1sq**3-3*p1sq**2*(p2sq+s12)-10*(p2sq**3+p2sq**2*s12+p2sq
     &  *s12**2+s12**3)-p1sq*(8*p2sq*s12+6*(p2sq**2+s12**2)))/1.2096d6
      Inv601111=(18*p2sq**2+10*p2sq*s12+5*(p1sq**2+s12**2+p1sq*(2*p2sq+
     &  s12)))/21600.d0
      Inv602111=(5*p1sq**2+16*p1sq*p2sq+36*p2sq**2+10*p1sq*s12+24*p2sq*
     &  s12+15*s12**2)/86400.d0
      Inv602211=(7*p1sq**2+28*p1sq*p2sq+72*p2sq**2+21*p1sq*s12+56*p2sq*
     &  s12+42*s12**2)/302400.d0
      Inv602221=(7*p1sq**2+32*p1sq*p2sq+90*p2sq**2+28*p1sq*s12+80*p2sq*
     &  s12+70*s12**2)/604800.d0
      Inv602222=(p1sq**2+5*p1sq*(p2sq+s12)+15*(p2sq**2+p2sq*s12+s12**2)
     &  )/151200.d0
      Inv40111111=(-p1sq-3*p2sq-s12)/180.d0
      Inv40211111=(-2*p1sq-9*p2sq-4*s12)/1080.d0
      Inv40221111=(-5*p1sq-27*p2sq-15*s12)/5400.d0
      Inv40222111=(-p1sq-6*p2sq-4*s12)/1800.d0
      Inv40222211=(-7*p1sq-45*p2sq-35*s12)/18900.d0
      Inv40222221=(-4*p1sq-27*p2sq-24*s12)/15120.d0
      Inv40222222=(-p1sq-7*(p2sq+s12))/5040.d0
      Inv2011111111=1d0/5d0
      Inv2021111111=1d0/10d0
      Inv2022111111=1d0/15d0
      Inv2022211111=1d0/20d0
      Inv2022221111=1d0/25d0
      Inv2022222111=1d0/30d0
      Inv2022222211=1d0/35d0
      Inv2022222221=1d0/40d0
      Inv2022222222=1d0/45d0
      Inv1001=(-12*p1sq**4-p1sq**2*(16*p2sq**2+21*p2sq*s12+12*s12**2)-3*
     &  p1sq*(6*p2sq**3+8*p2sq**2*s12+7*p2sq*s12**2+4*s12**3)-2*(10*p2sq
     &  **4+9*p2sq**3*s12+8*p2sq**2*s12**2+7*p2sq*s12**3+6*s12**4+p1sq**
     &  3*(7*p2sq+6*s12)))/2.66112d7
      Inv1002=(-2*p1sq**4-10*(p2sq**4+p2sq**3*s12+p2sq**2*s12**2+p2sq*s1
     &  2**3+s12**4)-4*p1sq*(p2sq+s12)*(p1sq**2+p2sq*s12+2*(p2sq**2+s12*
     &  *2))-p1sq**2*(9*p2sq*s12+6*(p2sq**2+s12**2)))/2.66112d7
      Inv80111=(105*p1sq**3+21*p1sq**2*(8*p2sq+5*s12)+7*p1sq*(36*p2sq**
     &  2+32*p2sq*s12+15*s12**2)+3*(120*p2sq**3+84*p2sq**2*s12+56*p2sq*s
     &  12**2+35*s12**3))/1.33056d7
      Inv80211=(21*p1sq**3+14*p1sq**2*(4*p2sq+3*s12)+p1sq*(108*p2sq**2+
     &  112*p2sq*s12+63*s12**2)+4*(45*p2sq**3+36*p2sq**2*s12+28*p2sq*s12
     &  **2+21*s12**3))/1.33056d7
      Inv80221=(7*p1sq**3+3*p1sq**2*(8*p2sq+7*s12)+p1sq*(54*p2sq**2+64*
     &  p2sq*s12+42*s12**2)+10*(10*p2sq**3+9*p2sq**2*s12+8*p2sq*s12**2+7
     &  *s12**3))/1.33056d7
      Inv80222=(3*p1sq**3+12*p1sq**2*(p2sq+s12)+60*(p2sq**3+p2sq**2*s12
     &  +p2sq*s12**2+s12**3)+10*p1sq*(4*p2sq*s12+3*(p2sq**2+s12**2)))/1.
     &  33056d7
      Inv6011111=(-18*p2sq**2-9*p2sq*(p1sq+s12)-4*(p1sq**2+p1sq*s12+s12
     &  **2))/23760.d0
      Inv6021111=(-5*p1sq**2-18*p1sq*p2sq-45*p2sq**2-10*p1sq*s12-27*p2s
     &  q*s12-15*s12**2)/118800.d0
      Inv6022111=(-14*p1sq**2-63*p1sq*p2sq-180*p2sq**2-42*(p1sq+3*p2sq)
     &  *s12-84*s12**2)/831600.d0
      Inv6022211=(-14*p1sq**2-72*p1sq*p2sq-225*p2sq**2-56*p1sq*s12-180*
     &  p2sq*s12-140*s12**2)/1.6632d6
      Inv6022221=(-8*p1sq**2-5*p1sq*(9*p2sq+8*s12)-15*(10*p2sq**2+9*p2s
     &  q*s12+8*s12**2))/1.6632d6
      Inv6022222=(-p1sq**2-6*p1sq*(p2sq+s12)-21*(p2sq**2+p2sq*s12+s12**
     &  2))/332640.d0
      Inv401111111=(10*p2sq+3*(p1sq+s12))/660.d0
      Inv402111111=(p1sq+5*p2sq+2*s12)/660.d0
      Inv402211111=(p1sq+6*p2sq+3*s12)/1320.d0
      Inv402221111=(3*p1sq+20*p2sq+12*s12)/6600.d0
      Inv402222111=(7*p1sq+50*p2sq+35*s12)/23100.d0
      Inv402222211=(2*p1sq+15*p2sq+12*s12)/9240.d0
      Inv402222221=(9*p1sq+70*p2sq+63*s12)/55440.d0
      Inv402222222=(p1sq+8*(p2sq+s12))/7920.d0
      Inv20111111111=-1d0/5.5d0
      Inv20211111111=-1d0/11d0
      Inv20221111111=-1d0/16.5d0
      Inv20222111111=-1d0/22d0
      Inv20222211111=-1d0/27.5d0
      Inv20222221111=-1d0/33d0
      Inv20222222111=-1d0/38.5d0
      Inv20222222211=-1d0/44d0
      Inv20222222221=-1d0/49.5d0
      Inv20222222222=-1d0/55d0
      Inv120=(-5*p1sq**5-5*p1sq**4*(p2sq+s12)-p1sq**3*(5*p2sq**2+8*p2sq*s12+5*s12**2
     &  )-p1sq**2*(5*p2sq**3+9*p2sq**2*s12+9*p2sq*s12**2+5*s12**3)-p1sq*(5*p2
     &  sq**4+8*p2sq**3*s12+9*p2sq**2*s12**2+8*p2sq*s12**3+5*s12**4)-5*(p2sq*
     &  *5+p2sq**4*s12+p2sq**3*s12**2+p2sq**2*s12**3+p2sq*s12**4+s12**5))/3.
     &  193344d8
         Inv10011=(55*p2sq**4+45*p2sq**3*s12+36*p2sq**2*s12**2+28*p2sq*s12
     &  **3+7*p1sq**3*(4*p2sq+3*s12)+21*(p1sq**4+s12**4)+3*p1sq*(15*p2sq**3
     &  +18*p2sq**2*s12+14*p2sq*s12**2+7*s12**3+p1sq*(12*p2sq**2+14*p2sq*
     &  s12+7*s12**2)))/7.98336d7
      Inv10021=(7*p1sq**4+2*p1sq**3*(8*p2sq+7*s12)+3*p1sq**2*(9*p2sq**2+12
     &  *p2sq*s12+7*s12**2)+p1sq*(40*p2sq**3+54*p2sq**2*s12+48*p2sq*s12**
     &  2+28*s12**3)+5*(11*p2sq**4+10*p2sq**3*s12+9*p2sq**2*s12**2+8*p2s
     &  q*s12**3+7*s12**4))/1.596672d8
      Inv10022=(p1sq**4+3*p1sq**3*(p2sq+s12)+5*p1sq*(2*p2sq**3+3*p2sq**2*s
     &  12+3*p2sq*s12**2+2*s12**3)+15*(p2sq**4+p2sq**3*s12+p2sq**2*s12**
     &  2+p2sq*s12**3+s12**4)+p1sq**2*(9*p2sq*s12+6*(p2sq**2+s12**2)))/7.
     &  98336d7
      Inv801111=(-165*p2sq**3-105*p2sq**2*s12-63*p2sq*s12**2-7*p1sq*(9*
     &  p1sq*p2sq+15*p2sq**2+5*p1sq*s12+12*p2sq*s12+5*s12**2)-35*(p1sq**
     &  3+s12**3))/6.6528d6
      Inv802111=(-165*p2sq**3-120*p2sq**2*s12-84*p2sq*s12**2-56*s12**3-
     &  14*p1sq**2*(p1sq+3*p2sq+2*s12)-6*p1sq*(15*p2sq**2+14*p2sq*s12+7*
     &  s12**2))/1.33056d7
      Inv802211=(-14*p1sq**3-6*p1sq**2*(9*p2sq+7*s12)-3*p1sq*(45*p2sq**
     &  2+48*p2sq*s12+28*s12**2)-5*(55*p2sq**3+45*p2sq**2*s12+36*p2sq*s1
     &  2**2+28*s12**3))/3.99168d7
      Inv802221=(-2*p1sq**3-p1sq**2*(9*p2sq+8*s12)-5*(11*p2sq**3+10*p2s
     &  q**2*s12+9*p2sq*s12**2+8*s12**3+p1sq*(5*p2sq**2+6*p2sq*s12+4*s12
     &  **2)))/1.33056d7
      Inv802222=(-p1sq**3-35*(p2sq**3+p2sq**2*s12+p2sq*s12**2+s12**3)-5
     &  *p1sq*(4*p2sq*s12+p1sq*(p2sq+s12)+3*(p2sq**2+s12**2)))/1.33056d7
      Inv60111111=(11*p2sq**2+5*p2sq*(p1sq+s12)+2*(p1sq**2+p1sq*s12+s12
     &  **2))/15840.d0
      Inv60211111=(p1sq**2+4*p1sq*p2sq+11*p2sq**2+2*p1sq*s12+6*p2sq*s12
     &  +3*s12**2)/31680.d0
      Inv60221111=(7*p1sq**2+35*p1sq*p2sq+110*p2sq**2+21*p1sq*s12+70*p2
     &  sq*s12+42*s12**2)/554400.d0
      Inv60222111=(14*p1sq**2+80*p1sq*p2sq+275*p2sq**2+56*p1sq*s12+200*
     &  p2sq*s12+140*s12**2)/2.2176d6
      Inv60222211=(12*p1sq**2+75*p1sq*p2sq+275*p2sq**2+60*p1sq*s12+225*
     &  p2sq*s12+180*s12**2)/3.3264d6
      Inv60222221=(3*p1sq**2+20*p1sq*p2sq+77*p2sq**2+18*p1sq*s12+70*p2s
     &  q*s12+63*s12**2)/1.33056d6
      Inv60222222=(p1sq**2+7*p1sq*(p2sq+s12)+28*(p2sq**2+p2sq*s12+s12**
     &  2))/665280.d0
      Inv4011111111=(-11*p2sq-3*(p1sq+s12))/792.d0
      Inv4021111111=(-2*p1sq-11*p2sq-4*s12)/1584.d0
      Inv4022111111=(-5*p1sq-33*p2sq-15*s12)/7920.d0
      Inv4022211111=(-3*p1sq-22*p2sq-12*s12)/7920.d0
      Inv4022221111=(-7*p1sq-55*p2sq-35*s12)/27720.d0
      Inv4022222111=(-4*p1sq-33*p2sq-24*s12)/22176.d0
      Inv4022222211=(-9*p1sq-77*p2sq-63*s12)/66528.d0
      Inv4022222221=(-5*p1sq-44*p2sq-40*s12)/47520.d0
      Inv4022222222=(-p1sq-9*(p2sq+s12))/11880.d0
      Inv201111111111=1d0/6d0
      Inv202111111111=1d0/12d0
      Inv202211111111=1d0/18d0
      Inv202221111111=1d0/24d0
      Inv202222111111=1d0/30d0
      Inv202222211111=1d0/36d0
      Inv202222221111=1d0/42d0
      Inv202222222111=1d0/48d0
      Inv202222222211=1d0/54d0
      Inv202222222221=1d0/60d0
      Inv202222222222=1d0/66d0
                       
      TX1=-2*p2sq*r10 + r21*(-p1sq - p2sq + s12)
      TX2=-2*p1sq*r21 + r10*(-p1sq - p2sq + s12)

      If(abs(TX1).gt.abs(TX2))then
         jj=1
         jjinit=1
         IX=1d0/TX1
        else
         jjinit=2
         jj=2
         IX=1d0/TX2
      endif

       kinit=k
       linit=l
       Zmaxinit=Zmax
       IXinit=IX


c FC %       if(jj.eq.1)  IX=1/TX1
c FC %       if(jj.eq.2)  IX=1/TX2
c FC %
c FC %       if((k+l).eq.2) Zmax=Z(1,1)
c FC %       if((k+l).eq.3) Zmax=Z(2,1)
c FC %       if((k*l).eq.4) Zmax=Z(2,2)

c FC %       print*, 'k',k
c FC %       print*, 'l',l
c FC %       print*, 'jj',jj
c FC %       print*, 'cont',cont

       ZZ(1,1,1,1)=0d0
       ZZ(2,1,1,1)=0d0
       ZZ(1,2,1,1)=0d0
       ZZ(2,2,1,1)=0d0
       ZZ(1,1,2,1)=0d0
       ZZ(2,1,2,1)=-1d0
       ZZ(1,2,2,1)=1d0
       ZZ(2,2,2,1)=0d0
       ZZ(1,1,1,2)=0d0
       ZZ(2,1,1,2)=1d0
       ZZ(1,2,1,2)=-1d0
       ZZ(2,2,1,2)=0d0
       ZZ(1,1,2,2)=0d0
       ZZ(2,1,2,2)=0d0
       ZZ(1,2,2,2)=0d0
       ZZ(2,2,2,2)=0d0

 100  I6Z=1d0/(ZMax*6d0)
      I10Z=1d0/(ZMax*10d0)
      I12Z=1d0/(ZMax*12d0)
      I14Z=1d0/(ZMax*14d0)
      I18Z=1d0/(ZMax*18d0)
      I22Z=1d0/(ZMax*22d0)
      I26Z=1d0/(ZMax*26d0)
      I30Z=1d0/(ZMax*30d0)
      I34Z=1d0/(ZMax*34d0)
      I38Z=1d0/(ZMax*38d0)
      I42Z=1d0/(ZMax*42d0)
      I46Z=1d0/(ZMax*46d0)

       F(1)=r10*ZZ(k,1,l,1)+r21*ZZ(k,2,l,1)
       F(2)=r10*ZZ(k,1,l,2)+r21*ZZ(k,2,l,2)
       F(3)=r10*r10*ZZ(k,1,l,1)+r21*(r10*(ZZ(k,1,l,2)
     &  +ZZ(k,2,l,1))+r21*ZZ(k,2,l,2))
       F(4)=2*r10*ZZ(k,1,l,1)+r21*(ZZ(k,1,l,2)+ZZ(k,2,l,1))
       F(5)=2*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+4*r21*ZZ(k,2,l,2)

       det3=4d0*(p1sq*p2sq-p1p2*p1p2)
       detABS=4d0*(abs(p1sq*p2sq)+abs(p1p2*p1p2))
       ratio= abs(det3/detABS)
       line= (0.0025d0-0.001d0)/(0.0118d0-0.007d0)*(ratio) 
     &  +0.0025d0-(0.0025d0-0.001d0)/(0.0118d0-0.007d0)*(0.0118d0)
cFC %       print*, "F", Real(F(1))
cFC %       print*, "F", Real(F(2))
cFC %       print*, "F", Real(F(3))
cFC %       print*, "F", Real(F(4))
cFC %       print*, "F", Real(F(5))

c                Iteration0
c                Step1
       S21(1)=B013-B023
       S21(2)=B012-B013
       auxC30=-(S21(1)*Z(jj,1))-S21(2)*Z(jj,2)
       tempC30=auxC30*IX
       if(order.eq.0) goto 500
c                Iteration1
c                Step1
       S221(1)=B023+B13(1,1)
       S221(2)=B13(1,1)-B23(1,1)
       S222(1)=B12(1,1)-B13(1,1)
       S222(2)=-B13(1,1)
       S200=2*B023
       auxC300=-(F(1)*S21(1))-F(2)*S21(2)+S221(k)*Z(1,l)+S222(k)*Z(2,l)-
     &  (2*Inv00-S200+S221(1)+S222(2))*Z(k,l)
       tempC300=I6Z*(auxC300+tempC30*F(3))
       aux1(1)=-(S221(1)*Z(jj,1))-S222(1)*Z(jj,2)
       aux1(2)=-(S221(2)*Z(jj,1))-S222(2)*Z(jj,2)
       temp1(1)=IX*(aux1(1)+2*tempC300*Z(jj,1))
       temp1(2)=IX*(aux1(2)+2*tempC300*Z(jj,2))
c                Step2
       tempC30=IX*(auxC30+det3*temp1(jj))
       if(order.eq.1) goto 500
c                Iteration2
c                Step1
       S2311(1)=-B023+B13(1,2)
       S2311(2)=B13(1,2)+B23(1,1)
       S2312(1)=B13(1,2)+B23(1,1)
       S2312(2)=B13(1,2)-B23(1,2)
       S2321(1)=B12(1,2)-B13(1,2)
       S2321(2)=-B13(1,2)
       S2322(1)=-B13(1,2)
       S2322(2)=-B13(1,2)
       S2h001(1)=B13(2,2)-B23(2,2)
       S2h001(2)=B12(2,2)-B13(2,2)
       S2001(1)=-2*B023
       S2001(2)=2*B23(1,1)
       aux001(1)=-(F(1)*S221(1))-F(2)*S222(1)+S2311(k)*Z(1,l)+S2321(k)*Z
     &  (2,l)+(-2*Inv001+S2001(1)-S2311(1)-S2322(1))*Z(k,l)-2*(S2h001(1)
     &  *ZZ(k,1,l,1)+S2h001(2)*ZZ(k,1,l,2))
       aux001(2)=-(F(1)*S221(2))-F(2)*S222(2)+S2312(k)*Z(1,l)+S2322(k)*Z
     &  (2,l)+(-2*Inv002+S2001(2)-S2312(1)-S2322(2))*Z(k,l)-2*(S2h001(1)
     &  *ZZ(k,2,l,1)+S2h001(2)*ZZ(k,2,l,2))
       temp001(1)=I10Z*(aux001(1)+2*tempC300*F(4)+F(3)*temp1(1))
       temp001(2)=I10Z*(aux001(2)+tempC300*F(5)+F(3)*temp1(2))
       aux2(1,1)=-(S2311(1)*Z(jj,1))-S2321(1)*Z(jj,2)
       aux2(2,1)=-(S2312(1)*Z(jj,1))-S2322(1)*Z(jj,2)
       aux2(2,2)=-(S2312(2)*Z(jj,1))-S2322(2)*Z(jj,2)
       temp2(1,1)=IX*(aux2(1,1)+4*temp001(1)*Z(jj,1))
       temp2(2,1)=IX*(aux2(2,1)+2*temp001(2)*Z(jj,1)+2*temp001(1)*Z(jj,2
     &  ))
       temp2(2,2)=IX*(aux2(2,2)+4*temp001(2)*Z(jj,2))
       temp2(1,2)=temp2(2,1)
c                Step2
       tempC300=I6Z*(auxC300+tempC30*F(3)-det3*temp2(k,l))
       temp1(1)=IX*(aux1(1)+det3*temp2(1,jj)+2*tempC300*Z(jj,1))
       temp1(2)=IX*(aux1(2)+det3*temp2(2,jj)+2*tempC300*Z(jj,2))
c                Step3
       tempC30=IX*(auxC30+det3*temp1(jj))
!        C30=tempC30
!            Cij(1,1)=temp1(1)
!            Cij(2,1)=temp1(2)
!            Cij(1,2)=temp2(1,1)
!            Cij(2,2)=temp2(2,2)
!            Cij(3,2)=temp2(2,1)
!            Cij(4,2)=tempC300
!            Cij(1,3)=temp3(1,1,1)
!            Cij(2,3)=temp3(2,2,2)
!            Cij(3,3)=temp3(2,1,1)
!            Cij(4,3)=temp3(2,2,1)
!            Cij(5,3)=temp001(1)
!            Cij(6,3)=temp001(2)
          tempjj1=abs(det3*temp1(1)/TX1)
          tempjj2=abs(det3*temp1(2)/TX2)
          If(tempjj1.lt.tempjj2)then
             jjtemp=1
             IXtemp=1d0/TX1
             tempjj=tempjj1/abs(auxC30)*abs(TX1)
          else
             jjtemp=2
             IXtemp=1d0/TX2
             tempjj=tempjj2/abs(auxC30)*abs(TX2)
          endif
          tempkl1=abs(det3*temp2(1,1)/(Z(1,1)*6d0))
          tempkl2=abs(det3*temp2(2,1)/(Z(2,1)*6d0))
          tempkl3=abs(det3*temp2(1,2)/(Z(1,2)*6d0))
          tempkl4=abs(det3*temp2(2,2)/(Z(2,2)*6d0)) 
          tempkl4=tempkl4*0.62d0


cFC %          print*, "Ktemp", k
cFC %          print*, "ltemp", l,Zmax,Z(k,l),IX
cFC %          PRINT*,  "JJ", JJ
cFC %           print*,'tempjj1',tempjj1
cFC %          print*,'tempjj2',tempjj2
cFC %          print*,'tempkl1',tempkl1
cFC %          print*,'tempkl2',tempkl2
cFC %          print*,'tempkl3',tempkl3
cFC %          print*,'tempkl4',tempkl4
cFC %          print*,'temp11',temp2(1,1)
cFC %          print*,'temp21',temp2(2,1)
cFC %          print*,'temp12',temp2(1,2)
cFC %          print*,'temp22',temp2(2,2)
cFC %          print*,'temp1',temp1(1)
cFC %          print*,'temp2',temp1(2)

          if(tempkl1.lt.tempkl2) then
             if(tempkl1.lt.tempkl4) then
                  ktemp=1
                  ltemp=1
                  Zmaxtemp=Z(1,1)
          F1=r10*r10*ZZ(ktemp,1,ltemp,1)+r21*(r10*(ZZ(ktemp,1,ltemp,2)
     &  +ZZ(ktemp,2,ltemp,1))+r21*ZZ(ktemp,2,ltemp,2))
                  tempkl=tempkl1/abs(auxC300+tempC30*F1)*abs((Z(1,1)*6d0))
             else
                  ktemp=2
                  ltemp=2
                 Zmaxtemp=Z(2,2)
          F1=r10*r10*ZZ(ktemp,1,ltemp,1)+r21*(r10*(ZZ(ktemp,1,ltemp,2)
     &  +ZZ(ktemp,2,ltemp,1))+r21*ZZ(ktemp,2,ltemp,2))
                  tempkl=tempkl4/abs(auxC300+tempC30*F1)*abs((Z(2,2)*6d0))
             endif
          else
             if(tempkl2.lt.tempkl4) then
                  ktemp=2
                  ltemp=1
                 Zmaxtemp=Z(2,1)
          F1=r10*r10*ZZ(ktemp,1,ltemp,1)+r21*(r10*(ZZ(ktemp,1,ltemp,2)
     &  +ZZ(ktemp,2,ltemp,1))+r21*ZZ(ktemp,2,ltemp,2))
                  tempkl=tempkl2/abs(auxC300+tempC30*F1)*abs((Z(2,1)*6d0))
             else
                  ktemp=2
                  ltemp=2
                 Zmaxtemp=Z(2,2)
          F1=r10*r10*ZZ(ktemp,1,ltemp,1)+r21*(r10*(ZZ(ktemp,1,ltemp,2)
     &  +ZZ(ktemp,2,ltemp,1))+r21*ZZ(ktemp,2,ltemp,2))
                  tempkl=tempkl4/abs(auxC300+tempC30*F1)*abs((Z(2,2)*6d0))
             endif
          endif
  
      if(cont.eq.4)  goto 200
      if(cont.lt.4) then
         if(ktemp.eq.k) then
             if (ltemp.eq.l) then
                  if (jjtemp.eq.jj) then
                    goto 200
                  else
                   jj=jjtemp
                   IX=IXtemp
                   k=ktemp
                   l=ltemp
                   Zmax=Zmaxtemp
                   cont=cont+1
                 
                   goto 100
                 endif
              else
                   jj=jjtemp
                   IX=IXtemp
                   k=ktemp
                   l=ltemp
                   Zmax=Zmaxtemp
                   cont=cont+1
                   goto 100 
              endif
            else
                    jj=jjtemp
                   IX=IXtemp
                   k=ktemp
                   l=ltemp
                   Zmax=Zmaxtemp
                   cont=cont+1
                   goto 100
             endif
             else
                   jj=jjinit
                   IX=IXinit
                   k=kinit
                   l=linit
                   Zmax=Zmaxinit
                   goto 100
                   cont=4
             endif

 200   If((abs(tempkl).ge.2d-3).and.(abs(tempjj).ge.1d-2)) then
       !return
       endif
       
!       If(tempjj.ge.1d-2) then
!       return
!       endif
c FC %       
       if(abs(IX*det3).ge.3.2d-1) then
!.or.(abs(IX*det3).gt.abs(20d0*line))
          if(printmy)then
         print*, "CAUTION SMALL CAYLEY,IX*det3 in C",abs(IX*det3)
         endif
c         print*, "LINE should be larger than IX*det3",abs(line)
c          print*
       !return
       endif     
c FC %       print*, 'k',k
c FC %       print*, 'l',l
c FC %       print*, 'jj',jj
c FC %       print*, 'cont',cont
       if(order.eq.2) goto 500
c                Iteration3
c                Step1
       S2h0021(1)=B13(2,3)+B23(2,2)
       S2h0021(2)=B13(2,3)-B23(2,3)
       S2h0022(1)=B12(2,3)-B13(2,3)
       S2h0022(2)=-B13(2,3)
       S20000=2*B23(2,2)
       auxC30000=-(F(1)*S2h001(1))-F(2)*S2h001(2)+S2h0021(k)*Z(1,l)+S2h0
     &  022(k)*Z(2,l)+(-2*Inv0000+S20000-S2h0021(1)-S2h0022(2))*Z(k,l)
       tempC30000=I10Z*(auxC30000+tempC300*F(3))
       S20021(1)=2*B023
       S20021(2)=-2*B23(1,1)
       S20022(1)=-2*B23(1,1)
       S20022(2)=2*B23(1,2)
       S24111(1)=B023+B13(1,3)
       S24111(2)=B13(1,3)-B23(1,1)
       S24121(1)=B13(1,3)-B23(1,1)
       S24121(2)=B13(1,3)+B23(1,2)
       S24122(1)=B13(1,3)+B23(1,2)
       S24122(2)=B13(1,3)-B23(1,3)
       S24211(1)=B12(1,3)-B13(1,3)
       S24211(2)=-B13(1,3)
       S24221(1)=-B13(1,3)
       S24221(2)=-B13(1,3)
       S24222(1)=-B13(1,3)
       S24222(2)=-B13(1,3)
       aux002(1,1)=-(F(1)*S2311(1))-F(2)*S2321(1)+S24111(k)*Z(1,l)+S2421
     &  1(k)*Z(2,l)+(-2*Inv0011+S20021(1)-S24111(1)-S24221(1))*Z(k,l)-4*
     &  (S2h0021(1)*ZZ(k,1,l,1)+S2h0022(1)*ZZ(k,1,l,2))
       aux002(2,1)=-(F(1)*S2312(1))-F(2)*S2322(1)+S24121(k)*Z(1,l)+S2422
     &  1(k)*Z(2,l)+(-2*Inv0021+S20022(1)-S24121(1)-S24222(1))*Z(k,l)-2*
     &  (S2h0021(2)*ZZ(k,1,l,1)+S2h0022(2)*ZZ(k,1,l,2))-2*S2h0021(1)*ZZ(
     &  k,2,l,1)-2*S2h0022(1)*ZZ(k,2,l,2)
       aux002(2,2)=-(F(1)*S2312(2))-F(2)*S2322(2)+S24122(k)*Z(1,l)+S2422
     &  2(k)*Z(2,l)+(-2*Inv0022+S20022(2)-S24122(1)-S24222(2))*Z(k,l)-4*
     &  (S2h0021(2)*ZZ(k,2,l,1)+S2h0022(2)*ZZ(k,2,l,2))
       temp002(1,1)=I14Z*(aux002(1,1)+4*F(4)*temp001(1)+F(3)*temp2(1,1)+
     &  8*tempC30000*ZZ(k,1,l,1))
       temp002(2,1)=I14Z*(aux002(2,1)+F(5)*temp001(1)+2*F(4)*temp001(2)+
     &  F(3)*temp2(2,1)+4*tempC30000*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp002(2,2)=I14Z*(aux002(2,2)+2*F(5)*temp001(2)+F(3)*temp2(2,2)+
     &  8*tempC30000*ZZ(k,2,l,2))
       temp002(1,2)=temp002(2,1)
       aux3(1,1,1)=-(S24111(1)*Z(jj,1))-S24211(1)*Z(jj,2)
       aux3(2,1,1)=-(S24121(1)*Z(jj,1))-S24221(1)*Z(jj,2)
       aux3(2,2,1)=-(S24122(1)*Z(jj,1))-S24222(1)*Z(jj,2)
       aux3(2,2,2)=-(S24122(2)*Z(jj,1))-S24222(2)*Z(jj,2)
       temp3(1,1,1)=IX*(aux3(1,1,1)+6*temp002(1,1)*Z(jj,1))
       temp3(2,1,1)=IX*(aux3(2,1,1)+4*temp002(2,1)*Z(jj,1)+2*temp002(1,1
     &  )*Z(jj,2))
       temp3(2,2,1)=IX*(aux3(2,2,1)+2*temp002(2,2)*Z(jj,1)+4*temp002(2,1
     &  )*Z(jj,2))
       temp3(2,2,2)=IX*(aux3(2,2,2)+6*temp002(2,2)*Z(jj,2))
       temp3(1,1,2)=temp3(2,1,1)
       temp3(1,2,1)=temp3(2,1,1)
       temp3(1,2,2)=temp3(2,2,1)
       temp3(2,1,2)=temp3(2,2,1)
c                Step2
       temp001(1)=I10Z*(aux001(1)+2*tempC300*F(4)+F(3)*temp1(1)-det3*tem
     &  p3(1,k,l))
       temp001(2)=I10Z*(aux001(2)+tempC300*F(5)+F(3)*temp1(2)-det3*temp3
     &  (2,k,l))
       temp2(1,1)=IX*(aux2(1,1)+det3*temp3(1,1,jj)+4*temp001(1)*Z(jj,1))
       temp2(2,1)=IX*(aux2(2,1)+det3*temp3(2,1,jj)+2*(temp001(2)*Z(jj,1)
     &  +temp001(1)*Z(jj,2)))
       temp2(2,2)=IX*(aux2(2,2)+det3*temp3(2,2,jj)+4*temp001(2)*Z(jj,2))
       temp2(1,2)=temp2(2,1)
c                Step3
       tempC300=I6Z*(auxC300+tempC30*F(3)-det3*temp2(k,l))
       temp1(1)=IX*(aux1(1)+det3*temp2(1,jj)+2*tempC300*Z(jj,1))
       temp1(2)=IX*(aux1(2)+det3*temp2(2,jj)+2*tempC300*Z(jj,2))
c                Step4
       tempC30=IX*(auxC30+det3*temp1(jj))
       if(order.eq.3) goto 500
c                Iteration4
c                Step1
       S2h00311(1)=B13(2,4)-B23(2,2)
       S2h00311(2)=B13(2,4)+B23(2,3)
       S2h00312(1)=B13(2,4)+B23(2,3)
       S2h00312(2)=B13(2,4)-B23(2,4)
       S2h00321(1)=B12(2,4)-B13(2,4)
       S2h00321(2)=-B13(2,4)
       S2h00322(1)=-B13(2,4)
       S2h00322(2)=-B13(2,4)
       S200001(1)=-2*B23(2,2)
       S200001(2)=2*B23(2,3)
       S2h00001(1)=B13(3,4)-B23(3,4)
       S2h00001(2)=B12(3,4)-B13(3,4)
       aux00001(1)=-(F(1)*S2h0021(1))-F(2)*S2h0022(1)+S2h00311(k)*Z(1,l)
     &  +S2h00321(k)*Z(2,l)+(-2*Inv00001+S200001(1)-S2h00311(1)-S2h00322
     &  (1))*Z(k,l)-2*(S2h00001(1)*ZZ(k,1,l,1)+S2h00001(2)*ZZ(k,1,l,2))
       aux00001(2)=-(F(1)*S2h0021(2))-F(2)*S2h0022(2)+S2h00312(k)*Z(1,l)
     &  +S2h00322(k)*Z(2,l)+(-2*Inv00002+S200001(2)-S2h00312(1)-S2h00322
     &  (2))*Z(k,l)-2*(S2h00001(1)*ZZ(k,2,l,1)+S2h00001(2)*ZZ(k,2,l,2))
       temp00001(1)=I14Z*(aux00001(1)+2*tempC30000*F(4)+F(3)*temp001(1))
       temp00001(2)=I14Z*(aux00001(2)+tempC30000*F(5)+F(3)*temp001(2))
       S200311(1)=-2*B023
       S200311(2)=2*B23(1,1)
       S200312(1)=2*B23(1,1)
       S200312(2)=-2*B23(1,2)
       S200321(1)=2*B23(1,1)
       S200321(2)=-2*B23(1,2)
       S200322(1)=-2*B23(1,2)
       S200322(2)=2*B23(1,3)
       S251111(1)=-B023+B13(1,4)
       S251111(2)=B13(1,4)+B23(1,1)
       S251211(1)=B13(1,4)+B23(1,1)
       S251211(2)=B13(1,4)-B23(1,2)
       S251221(1)=B13(1,4)-B23(1,2)
       S251221(2)=B13(1,4)+B23(1,3)
       S251222(1)=B13(1,4)+B23(1,3)
       S251222(2)=B13(1,4)-B23(1,4)
       S252111(1)=B12(1,4)-B13(1,4)
       S252111(2)=-B13(1,4)
       S252211(1)=-B13(1,4)
       S252211(2)=-B13(1,4)
       S252221(1)=-B13(1,4)
       S252221(2)=-B13(1,4)
       S252222(1)=-B13(1,4)
       S252222(2)=-B13(1,4)
       aux003(1,1,1)=-(F(1)*S24111(1))-F(2)*S24211(1)+S251111(k)*Z(1,l)+
     &  S252111(k)*Z(2,l)+(-2*Inv00111+S200311(1)-S251111(1)-S252211(1))
     &  *Z(k,l)-6*(S2h00311(1)*ZZ(k,1,l,1)+S2h00321(1)*ZZ(k,1,l,2))
       aux003(2,1,1)=-(F(1)*S24121(1))-F(2)*S24221(1)+S251211(k)*Z(1,l)+
     &  S252211(k)*Z(2,l)+(-2*Inv00211+S200321(1)-S251211(1)-S252221(1))
     &  *Z(k,l)-4*(S2h00312(1)*ZZ(k,1,l,1)+S2h00322(1)*ZZ(k,1,l,2))-2*S2
     &  h00311(1)*ZZ(k,2,l,1)-2*S2h00321(1)*ZZ(k,2,l,2)
       aux003(2,2,1)=-(F(1)*S24122(1))-F(2)*S24222(1)+S251221(k)*Z(1,l)+
     &  S252221(k)*Z(2,l)+(-2*Inv00221+S200322(1)-S251221(1)-S252222(1))
     &  *Z(k,l)-2*(S2h00312(2)*ZZ(k,1,l,1)+S2h00322(2)*ZZ(k,1,l,2))-4*S2
     &  h00312(1)*ZZ(k,2,l,1)-4*S2h00322(1)*ZZ(k,2,l,2)
       aux003(2,2,2)=-(F(1)*S24122(2))-F(2)*S24222(2)+S251222(k)*Z(1,l)+
     &  S252222(k)*Z(2,l)+(-2*Inv00222+S200322(2)-S251222(1)-S252222(2))
     &  *Z(k,l)-6*(S2h00312(2)*ZZ(k,2,l,1)+S2h00322(2)*ZZ(k,2,l,2))
       temp003(1,1,1)=I18Z*(aux003(1,1,1)+6*F(4)*temp002(1,1)+F(3)*temp3
     &  (1,1,1)+24*temp00001(1)*ZZ(k,1,l,1))
       temp003(2,1,1)=I18Z*(aux003(2,1,1)+F(5)*temp002(1,1)+4*F(4)*temp0
     &  02(2,1)+F(3)*temp3(2,1,1)+8*(temp00001(2)*ZZ(k,1,l,1)+temp00001(
     &  1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))))
       temp003(2,2,1)=I18Z*(aux003(2,2,1)+2*F(5)*temp002(2,1)+2*F(4)*tem
     &  p002(2,2)+F(3)*temp3(2,2,1)+8*(temp00001(2)*(ZZ(k,1,l,2)+ZZ(k,2,
     &  l,1))+temp00001(1)*ZZ(k,2,l,2)))
       temp003(2,2,2)=I18Z*(aux003(2,2,2)+F(3)*temp3(2,2,2)+24*temp00001
     &  (2)*ZZ(k,2,l,2)+temp002(2,2)*(6*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+12
     &  *r21*ZZ(k,2,l,2)))
       temp003(1,1,2)=temp003(2,1,1)
       temp003(1,2,1)=temp003(2,1,1)
       temp003(1,2,2)=temp003(2,2,1)
       temp003(2,1,2)=temp003(2,2,1)
       aux41(1,1,1)=-(S251111(1)*Z(jj,1))-S252111(1)*Z(jj,2)
       aux42(1,1,1)=-(S251211(1)*Z(jj,1))-S252211(1)*Z(jj,2)
       aux42(2,1,1)=-(S251221(1)*Z(jj,1))-S252221(1)*Z(jj,2)
       aux42(2,2,1)=-(S251222(1)*Z(jj,1))-S252222(1)*Z(jj,2)
       aux42(2,2,2)=-(S251222(2)*Z(jj,1))-S252222(2)*Z(jj,2)
       temp41(1,1,1)=IX*(aux41(1,1,1)+8*temp003(1,1,1)*Z(jj,1))
       temp42(1,1,1)=IX*(aux42(1,1,1)+6*temp003(2,1,1)*Z(jj,1)+2*temp003
     &  (1,1,1)*Z(jj,2))
       temp42(2,1,1)=IX*(aux42(2,1,1)+4*(temp003(2,2,1)*Z(jj,1)+temp003(
     &  2,1,1)*Z(jj,2)))
       temp42(2,2,1)=IX*(aux42(2,2,1)+2*temp003(2,2,2)*Z(jj,1)+6*temp003
     &  (2,2,1)*Z(jj,2))
       temp42(2,2,2)=IX*(aux42(2,2,2)+8*temp003(2,2,2)*Z(jj,2))
       temp41(1,1,2)=temp42(1,1,1)
       temp41(1,2,1)=temp42(1,1,1)
       temp41(1,2,2)=temp42(2,1,1)
       temp42(1,1,2)=temp42(2,1,1)
       temp42(1,2,1)=temp42(2,1,1)
       temp42(1,2,2)=temp42(2,2,1)
       temp42(2,1,2)=temp42(2,2,1)
c                Step2
       tempC30000=I10Z*(auxC30000+tempC300*F(3)-det3*temp002(k,l))
       temp002(1,1)=I14Z*(aux002(1,1)+4*F(4)*temp001(1)+F(3)*temp2(1,1)-
     &  det3*temp41(1,k,l)+8*tempC30000*ZZ(k,1,l,1))
       temp002(2,1)=I14Z*(aux002(2,1)+F(5)*temp001(1)+2*F(4)*temp001(2)+
     &  F(3)*temp2(2,1)-det3*temp42(1,k,l)+4*tempC30000*(ZZ(k,1,l,2)+ZZ(
     &  k,2,l,1)))
       temp002(2,2)=I14Z*(aux002(2,2)+2*F(5)*temp001(2)+F(3)*temp2(2,2)-
     &  det3*temp42(2,k,l)+8*tempC30000*ZZ(k,2,l,2))
       temp002(1,2)=temp002(2,1)
       temp002(1,2)=temp002(2,1)
       temp3(1,1,1)=IX*(aux3(1,1,1)+det3*temp41(1,1,jj)+6*temp002(1,1)*Z
     &  (jj,1))
       temp3(2,1,1)=IX*(aux3(2,1,1)+det3*temp42(1,1,jj)+4*temp002(2,1)*Z
     &  (jj,1)+2*temp002(1,1)*Z(jj,2))
       temp3(2,2,1)=IX*(aux3(2,2,1)+det3*temp42(2,1,jj)+2*temp002(2,2)*Z
     &  (jj,1)+4*temp002(2,1)*Z(jj,2))
       temp3(2,2,2)=IX*(aux3(2,2,2)+det3*temp42(2,2,jj)+6*temp002(2,2)*Z
     &  (jj,2))
       temp3(1,1,2)=temp3(2,1,1)
       temp3(1,2,1)=temp3(2,1,1)
       temp3(1,2,2)=temp3(2,2,1)
       temp3(2,1,2)=temp3(2,2,1)
c                Step3
       temp001(1)=I10Z*(aux001(1)+2*tempC300*F(4)+F(3)*temp1(1)-det3*tem
     &  p3(1,k,l))
       temp001(2)=I10Z*(aux001(2)+tempC300*F(5)+F(3)*temp1(2)-det3*temp3
     &  (2,k,l))
       temp2(1,1)=IX*(aux2(1,1)+det3*temp3(1,1,jj)+4*temp001(1)*Z(jj,1))
       temp2(2,1)=IX*(aux2(2,1)+det3*temp3(2,1,jj)+2*(temp001(2)*Z(jj,1)
     &  +temp001(1)*Z(jj,2)))
       temp2(2,2)=IX*(aux2(2,2)+det3*temp3(2,2,jj)+4*temp001(2)*Z(jj,2))
       temp2(1,2)=temp2(2,1)
c                Step4
       tempC300=I6Z*(auxC300+tempC30*F(3)-det3*temp2(k,l))
       temp1(1)=IX*(aux1(1)+det3*temp2(1,jj)+2*tempC300*Z(jj,1))
       temp1(2)=IX*(aux1(2)+det3*temp2(2,jj)+2*tempC300*Z(jj,2))
c                Step5
       tempC30=IX*(auxC30+det3*temp1(jj))
       if(order.eq.4) goto 500
c                Iteration5
c                Step1
       S2000000=2*B23(3,4)
       S2000021(1)=2*B23(2,2)
       S2000021(2)=-2*B23(2,3)
       S2000022(1)=-2*B23(2,3)
       S2000022(2)=2*B23(2,4)
       S2h000021(1)=B13(3,5)+B23(3,4)
       S2h000021(2)=B13(3,5)-B23(3,5)
       S2h000022(1)=B12(3,5)-B13(3,5)
       S2h000022(2)=-B13(3,5)
       auxC3000000=-(F(1)*S2h00001(1))-F(2)*S2h00001(2)+S2h000021(k)*Z(1
     &  ,l)+S2h000022(k)*Z(2,l)+(-2*Inv000000+S2000000-S2h000021(1)-S2h0
     &  00022(2))*Z(k,l)
       tempC3000000=I14Z*(auxC3000000+tempC30000*F(3))
       S2h004111(1)=B13(2,5)+B23(2,2)
       S2h004111(2)=B13(2,5)-B23(2,3)
       S2h004121(1)=B13(2,5)-B23(2,3)
       S2h004121(2)=B13(2,5)+B23(2,4)
       S2h004122(1)=B13(2,5)+B23(2,4)
       S2h004122(2)=B13(2,5)-B23(2,5)
       S2h004211(1)=B12(2,5)-B13(2,5)
       S2h004211(2)=-B13(2,5)
       S2h004221(1)=-B13(2,5)
       S2h004221(2)=-B13(2,5)
       S2h004222(1)=-B13(2,5)
       S2h004222(2)=-B13(2,5)
       aux00002(1,1)=-(F(1)*S2h00311(1))-F(2)*S2h00321(1)+S2h004111(k)*Z
     &  (1,l)+S2h004211(k)*Z(2,l)+(-2*Inv000011+S2000021(1)-S2h004111(1)
     &  -S2h004221(1))*Z(k,l)-4*(S2h000021(1)*ZZ(k,1,l,1)+S2h000022(1)*Z
     &  Z(k,1,l,2))
       aux00002(2,1)=-(F(1)*S2h00312(1))-F(2)*S2h00322(1)+S2h004121(k)*Z
     &  (1,l)+S2h004221(k)*Z(2,l)+(-2*Inv000021+S2000022(1)-S2h004121(1)
     &  -S2h004222(1))*Z(k,l)-2*(S2h000021(2)*ZZ(k,1,l,1)+S2h000022(2)*Z
     &  Z(k,1,l,2))-2*S2h000021(1)*ZZ(k,2,l,1)-2*S2h000022(1)*ZZ(k,2,l,2
     &  )
       aux00002(2,2)=-(F(1)*S2h00312(2))-F(2)*S2h00322(2)+S2h004122(k)*Z
     &  (1,l)+S2h004222(k)*Z(2,l)+(-2*Inv000022+S2000022(2)-S2h004122(1)
     &  -S2h004222(2))*Z(k,l)-4*(S2h000021(2)*ZZ(k,2,l,1)+S2h000022(2)*Z
     &  Z(k,2,l,2))
       temp00002(1,1)=I18Z*(aux00002(1,1)+4*F(4)*temp00001(1)+F(3)*temp0
     &  02(1,1)+8*tempC3000000*ZZ(k,1,l,1))
       temp00002(2,1)=I18Z*(aux00002(2,1)+F(5)*temp00001(1)+2*F(4)*temp0
     &  0001(2)+F(3)*temp002(2,1)+4*tempC3000000*(ZZ(k,1,l,2)+ZZ(k,2,l,1
     &  )))
       temp00002(2,2)=I18Z*(aux00002(2,2)+2*F(5)*temp00001(2)+F(3)*temp0
     &  02(2,2)+8*tempC3000000*ZZ(k,2,l,2))
       temp00002(1,2)=temp00002(2,1)
       S2004111(1)=2*B023
       S2004111(2)=-2*B23(1,1)
       S2004121(1)=-2*B23(1,1)
       S2004121(2)=2*B23(1,2)
       S2004122(1)=2*B23(1,2)
       S2004122(2)=-2*B23(1,3)
       S2004211(1)=-2*B23(1,1)
       S2004211(2)=2*B23(1,2)
       S2004221(1)=2*B23(1,2)
       S2004221(2)=-2*B23(1,3)
       S2004222(1)=-2*B23(1,3)
       S2004222(2)=2*B23(1,4)
       S2611111(1)=B023+B13(1,5)
       S2611111(2)=B13(1,5)-B23(1,1)
       S2612111(1)=B13(1,5)-B23(1,1)
       S2612111(2)=B13(1,5)+B23(1,2)
       S2612211(1)=B13(1,5)+B23(1,2)
       S2612211(2)=B13(1,5)-B23(1,3)
       S2612221(1)=B13(1,5)-B23(1,3)
       S2612221(2)=B13(1,5)+B23(1,4)
       S2612222(1)=B13(1,5)+B23(1,4)
       S2612222(2)=B13(1,5)-B23(1,5)
       S2621111(1)=B12(1,5)-B13(1,5)
       S2621111(2)=-B13(1,5)
       S2622111(1)=-B13(1,5)
       S2622111(2)=-B13(1,5)
       S2622211(1)=-B13(1,5)
       S2622211(2)=-B13(1,5)
       S2622221(1)=-B13(1,5)
       S2622221(2)=-B13(1,5)
       S2622222(1)=-B13(1,5)
       S2622222(2)=-B13(1,5)
       aux0041(1,1,1)=-(F(1)*S251111(1))-F(2)*S252111(1)+S2611111(k)*Z(1
     &  ,l)+S2621111(k)*Z(2,l)+(-2*Inv001111+S2004111(1)-S2611111(1)-S26
     &  22111(1))*Z(k,l)-8*(S2h004111(1)*ZZ(k,1,l,1)+S2h004211(1)*ZZ(k,1
     &  ,l,2))
       aux0042(1,1,1)=-(F(1)*S251211(1))-F(2)*S252211(1)+S2612111(k)*Z(1
     &  ,l)+S2622111(k)*Z(2,l)+(-2*Inv002111+S2004211(1)-S2612111(1)-S26
     &  22211(1))*Z(k,l)-6*(S2h004121(1)*ZZ(k,1,l,1)+S2h004221(1)*ZZ(k,1
     &  ,l,2))-2*S2h004111(1)*ZZ(k,2,l,1)-2*S2h004211(1)*ZZ(k,2,l,2)
       aux0042(2,1,1)=-(F(1)*S251221(1))-F(2)*S252221(1)+S2612211(k)*Z(1
     &  ,l)+S2622211(k)*Z(2,l)+(-2*Inv002211+S2004221(1)-S2612211(1)-S26
     &  22221(1))*Z(k,l)-4*(S2h004122(1)*ZZ(k,1,l,1)+S2h004222(1)*ZZ(k,1
     &  ,l,2))-4*S2h004121(1)*ZZ(k,2,l,1)-4*S2h004221(1)*ZZ(k,2,l,2)
       aux0042(2,2,1)=-(F(1)*S251222(1))-F(2)*S252222(1)+S2612221(k)*Z(1
     &  ,l)+S2622221(k)*Z(2,l)+(-2*Inv002221+S2004222(1)-S2612221(1)-S26
     &  22222(1))*Z(k,l)-2*(S2h004122(2)*ZZ(k,1,l,1)+S2h004222(2)*ZZ(k,1
     &  ,l,2))-6*S2h004122(1)*ZZ(k,2,l,1)-6*S2h004222(1)*ZZ(k,2,l,2)
       aux0042(2,2,2)=-(F(1)*S251222(2))-F(2)*S252222(2)+S2612222(k)*Z(1
     &  ,l)+S2622222(k)*Z(2,l)+(-2*Inv002222+S2004222(2)-S2612222(1)-S26
     &  22222(2))*Z(k,l)-8*(S2h004122(2)*ZZ(k,2,l,1)+S2h004222(2)*ZZ(k,2
     &  ,l,2))
       temp0041(1,1,1)=I22Z*(aux0041(1,1,1)+8*F(4)*temp003(1,1,1)+F(3)*t
     &  emp41(1,1,1)+48*temp00002(1,1)*ZZ(k,1,l,1))
       temp0042(1,1,1)=I22Z*(aux0042(1,1,1)+F(5)*temp003(1,1,1)+6*F(4)*t
     &  emp003(2,1,1)+F(3)*temp42(1,1,1)+24*temp00002(2,1)*ZZ(k,1,l,1)+1
     &  2*temp00002(1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp0042(2,1,1)=I22Z*(aux0042(2,1,1)+2*F(5)*temp003(2,1,1)+4*F(4)
     &  *temp003(2,2,1)+F(3)*temp42(2,1,1)+16*temp00002(2,1)*(ZZ(k,1,l,2
     &  )+ZZ(k,2,l,1))+8*(temp00002(2,2)*ZZ(k,1,l,1)+temp00002(1,1)*ZZ(k
     &  ,2,l,2)))
       temp0042(2,2,1)=I22Z*(aux0042(2,2,1)+2*F(4)*temp003(2,2,2)+F(3)*t
     &  emp42(2,2,1)+12*temp00002(2,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+24*temp
     &  00002(2,1)*ZZ(k,2,l,2)+temp003(2,2,1)*(6*r10*(ZZ(k,1,l,2)+ZZ(k,2
     &  ,l,1))+12*r21*ZZ(k,2,l,2)))
       temp0042(2,2,2)=I22Z*(aux0042(2,2,2)+4*F(5)*temp003(2,2,2)+F(3)*t
     &  emp42(2,2,2)+48*temp00002(2,2)*ZZ(k,2,l,2))
       temp0041(1,1,2)=temp0042(1,1,1)
       temp0041(1,2,1)=temp0042(1,1,1)
       temp0041(1,2,2)=temp0042(2,1,1)
       temp0042(1,1,2)=temp0042(2,1,1)
       temp0042(1,2,1)=temp0042(2,1,1)
       temp0042(1,2,2)=temp0042(2,2,1)
       temp0042(2,1,2)=temp0042(2,2,1)
       aux511(1,1,1)=-(S2611111(1)*Z(jj,1))-S2621111(1)*Z(jj,2)
       aux521(1,1,1)=-(S2612111(1)*Z(jj,1))-S2622111(1)*Z(jj,2)
       aux522(1,1,1)=-(S2612211(1)*Z(jj,1))-S2622211(1)*Z(jj,2)
       aux522(2,1,1)=-(S2612221(1)*Z(jj,1))-S2622221(1)*Z(jj,2)
       aux522(2,2,1)=-(S2612222(1)*Z(jj,1))-S2622222(1)*Z(jj,2)
       aux522(2,2,2)=-(S2612222(2)*Z(jj,1))-S2622222(2)*Z(jj,2)
       temp511(1,1,1)=IX*(aux511(1,1,1)+10*temp0041(1,1,1)*Z(jj,1))
       temp521(1,1,1)=IX*(aux521(1,1,1)+8*temp0042(1,1,1)*Z(jj,1)+2*temp
     &  0041(1,1,1)*Z(jj,2))
       temp522(1,1,1)=IX*(aux522(1,1,1)+6*temp0042(2,1,1)*Z(jj,1)+4*temp
     &  0042(1,1,1)*Z(jj,2))
       temp522(2,1,1)=IX*(aux522(2,1,1)+4*temp0042(2,2,1)*Z(jj,1)+6*temp
     &  0042(2,1,1)*Z(jj,2))
       temp522(2,2,1)=IX*(aux522(2,2,1)+2*temp0042(2,2,2)*Z(jj,1)+8*temp
     &  0042(2,2,1)*Z(jj,2))
       temp522(2,2,2)=IX*(aux522(2,2,2)+10*temp0042(2,2,2)*Z(jj,2))
       temp511(1,1,2)=temp521(1,1,1)
       temp511(1,2,1)=temp521(1,1,1)
       temp511(1,2,2)=temp522(1,1,1)
       temp521(1,1,2)=temp522(1,1,1)
       temp521(1,2,1)=temp522(1,1,1)
       temp521(1,2,2)=temp522(2,1,1)
       temp522(1,1,2)=temp522(2,1,1)
       temp522(1,2,1)=temp522(2,1,1)
       temp522(1,2,2)=temp522(2,2,1)
       temp522(2,1,2)=temp522(2,2,1)
c                Step2
       temp00001(1)=I14Z*(aux00001(1)+2*tempC30000*F(4)+F(3)*temp001(1)-
     &  det3*temp003(1,k,l))
       temp00001(2)=I14Z*(aux00001(2)+tempC30000*F(5)+F(3)*temp001(2)-de
     &  t3*temp003(2,k,l))
       temp003(1,1,1)=I18Z*(aux003(1,1,1)+6*F(4)*temp002(1,1)+F(3)*temp3
     &  (1,1,1)-det3*temp511(1,k,l)+24*temp00001(1)*ZZ(k,1,l,1))
       temp003(2,1,1)=I18Z*(aux003(2,1,1)+F(5)*temp002(1,1)+4*F(4)*temp0
     &  02(2,1)+F(3)*temp3(2,1,1)-det3*temp521(1,k,l)+8*(temp00001(2)*ZZ
     &  (k,1,l,1)+temp00001(1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))))
       temp003(2,2,1)=I18Z*(aux003(2,2,1)+2*F(5)*temp002(2,1)+2*F(4)*tem
     &  p002(2,2)+F(3)*temp3(2,2,1)-det3*temp522(1,k,l)+8*(temp00001(2)*
     &  (ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp00001(1)*ZZ(k,2,l,2)))
       temp003(2,2,2)=I18Z*(aux003(2,2,2)+F(3)*temp3(2,2,2)-det3*temp522
     &  (2,k,l)+24*temp00001(2)*ZZ(k,2,l,2)+temp002(2,2)*(6*r10*(ZZ(k,1,
     &  l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2)))
       temp003(1,1,2)=temp003(2,1,1)
       temp003(1,2,1)=temp003(2,1,1)
       temp003(1,2,2)=temp003(2,2,1)
       temp003(2,1,2)=temp003(2,2,1)
       temp41(1,1,1)=IX*(aux41(1,1,1)+det3*temp511(1,1,jj)+8*temp003(1,1
     &  ,1)*Z(jj,1))
       temp42(1,1,1)=IX*(aux42(1,1,1)+det3*temp521(1,1,jj)+6*temp003(2,1
     &  ,1)*Z(jj,1)+2*temp003(1,1,1)*Z(jj,2))
       temp42(2,1,1)=IX*(aux42(2,1,1)+det3*temp522(1,1,jj)+4*(temp003(2,
     &  2,1)*Z(jj,1)+temp003(2,1,1)*Z(jj,2)))
       temp42(2,2,1)=IX*(aux42(2,2,1)+det3*temp522(2,1,jj)+2*temp003(2,2
     &  ,2)*Z(jj,1)+6*temp003(2,2,1)*Z(jj,2))
       temp42(2,2,2)=IX*(aux42(2,2,2)+det3*temp522(2,2,jj)+8*temp003(2,2
     &  ,2)*Z(jj,2))
       temp41(1,1,2)=temp42(1,1,1)
       temp41(1,2,1)=temp42(1,1,1)
       temp41(1,2,2)=temp42(2,1,1)
       temp42(1,1,2)=temp42(2,1,1)
       temp42(1,2,1)=temp42(2,1,1)
       temp42(1,2,2)=temp42(2,2,1)
       temp42(2,1,2)=temp42(2,2,1)
c                Step3
       tempC30000=I10Z*(auxC30000+tempC300*F(3)-det3*temp002(k,l))
       temp002(1,1)=I14Z*(aux002(1,1)+4*F(4)*temp001(1)+F(3)*temp2(1,1)-
     &  det3*temp41(1,k,l)+8*tempC30000*ZZ(k,1,l,1))
       temp002(2,1)=I14Z*(aux002(2,1)+F(5)*temp001(1)+2*F(4)*temp001(2)+
     &  F(3)*temp2(2,1)-det3*temp42(1,k,l)+4*tempC30000*(ZZ(k,1,l,2)+ZZ(
     &  k,2,l,1)))
       temp002(2,2)=I14Z*(aux002(2,2)+2*F(5)*temp001(2)+F(3)*temp2(2,2)-
     &  det3*temp42(2,k,l)+8*tempC30000*ZZ(k,2,l,2))
       temp002(1,2)=temp002(2,1)
       temp3(1,1,1)=IX*(aux3(1,1,1)+det3*temp41(1,1,jj)+6*temp002(1,1)*Z
     &  (jj,1))
       temp3(2,1,1)=IX*(aux3(2,1,1)+det3*temp42(1,1,jj)+4*temp002(2,1)*Z
     &  (jj,1)+2*temp002(1,1)*Z(jj,2))
       temp3(2,2,1)=IX*(aux3(2,2,1)+det3*temp42(2,1,jj)+2*temp002(2,2)*Z
     &  (jj,1)+4*temp002(2,1)*Z(jj,2))
       temp3(2,2,2)=IX*(aux3(2,2,2)+det3*temp42(2,2,jj)+6*temp002(2,2)*Z
     &  (jj,2))
       temp3(1,1,2)=temp3(2,1,1)
       temp3(1,2,1)=temp3(2,1,1)
       temp3(1,2,2)=temp3(2,2,1)
       temp3(2,1,2)=temp3(2,2,1)
c                Step4
       temp001(1)=I10Z*(aux001(1)+2*tempC300*F(4)+F(3)*temp1(1)-det3*tem
     &  p3(1,k,l))
       temp001(2)=I10Z*(aux001(2)+tempC300*F(5)+F(3)*temp1(2)-det3*temp3
     &  (2,k,l))
       temp2(1,1)=IX*(aux2(1,1)+det3*temp3(1,1,jj)+4*temp001(1)*Z(jj,1))
       temp2(2,1)=IX*(aux2(2,1)+det3*temp3(2,1,jj)+2*(temp001(2)*Z(jj,1)
     &  +temp001(1)*Z(jj,2)))
       temp2(2,2)=IX*(aux2(2,2)+det3*temp3(2,2,jj)+4*temp001(2)*Z(jj,2))
       temp2(1,2)=temp2(2,1)
c                Step5
       tempC300=I6Z*(auxC300+tempC30*F(3)-det3*temp2(k,l))
       temp1(1)=IX*(aux1(1)+det3*temp2(1,jj)+2*tempC300*Z(jj,1))
       temp1(2)=IX*(aux1(2)+det3*temp2(2,jj)+2*tempC300*Z(jj,2))
c                Step6
       tempC30=IX*(auxC30+det3*temp1(jj))
       if(order.eq.5) goto 500
c                Iteration6
c                Step1
       S20000001(1)=-2*B23(3,4)
       S20000001(2)=2*B23(3,5)
       S2h0000001(1)=B13(4,6)-B23(4,6)
       S2h0000001(2)=B12(4,6)-B13(4,6)
       S2h0000311(1)=B13(3,6)-B23(3,4)
       S2h0000311(2)=B13(3,6)+B23(3,5)
       S2h0000312(1)=B13(3,6)+B23(3,5)
       S2h0000312(2)=B13(3,6)-B23(3,6)
       S2h0000321(1)=B12(3,6)-B13(3,6)
       S2h0000321(2)=-B13(3,6)
       S2h0000322(1)=-B13(3,6)
       S2h0000322(2)=-B13(3,6)
       aux0000001(1)=-(F(1)*S2h000021(1))-F(2)*S2h000022(1)+S2h0000311(k
     &  )*Z(1,l)+S2h0000321(k)*Z(2,l)+(-2*Inv0000001+S20000001(1)-S2h000
     &  0311(1)-S2h0000322(1))*Z(k,l)-2*(S2h0000001(1)*ZZ(k,1,l,1)+S2h00
     &  00001(2)*ZZ(k,1,l,2))
       aux0000001(2)=-(F(1)*S2h000021(2))-F(2)*S2h000022(2)+S2h0000312(k
     &  )*Z(1,l)+S2h0000322(k)*Z(2,l)+(-2*Inv0000002+S20000001(2)-S2h000
     &  0312(1)-S2h0000322(2))*Z(k,l)-2*(S2h0000001(1)*ZZ(k,2,l,1)+S2h00
     &  00001(2)*ZZ(k,2,l,2))
       temp0000001(1)=I18Z*(aux0000001(1)+2*tempC3000000*F(4)+F(3)*temp0
     &  0001(1))
       temp0000001(2)=I18Z*(aux0000001(2)+tempC3000000*F(5)+F(3)*temp000
     &  01(2))
       S20000311(1)=-2*B23(2,2)
       S20000311(2)=2*B23(2,3)
       S20000312(1)=2*B23(2,3)
       S20000312(2)=-2*B23(2,4)
       S20000321(1)=2*B23(2,3)
       S20000321(2)=-2*B23(2,4)
       S20000322(1)=-2*B23(2,4)
       S20000322(2)=2*B23(2,5)
       S2h0051111(1)=B13(2,6)-B23(2,2)
       S2h0051111(2)=B13(2,6)+B23(2,3)
       S2h0051211(1)=B13(2,6)+B23(2,3)
       S2h0051211(2)=B13(2,6)-B23(2,4)
       S2h0051221(1)=B13(2,6)-B23(2,4)
       S2h0051221(2)=B13(2,6)+B23(2,5)
       S2h0051222(1)=B13(2,6)+B23(2,5)
       S2h0051222(2)=B13(2,6)-B23(2,6)
       S2h0052111(1)=B12(2,6)-B13(2,6)
       S2h0052111(2)=-B13(2,6)
       S2h0052211(1)=-B13(2,6)
       S2h0052211(2)=-B13(2,6)
       S2h0052221(1)=-B13(2,6)
       S2h0052221(2)=-B13(2,6)
       S2h0052222(1)=-B13(2,6)
       S2h0052222(2)=-B13(2,6)
       aux00003(1,1,1)=-(F(1)*S2h004111(1))-F(2)*S2h004211(1)+S2h0051111
     &  (k)*Z(1,l)+S2h0052111(k)*Z(2,l)+(-2*Inv0000111+S20000311(1)-S2h0
     &  051111(1)-S2h0052211(1))*Z(k,l)-6*(S2h0000311(1)*ZZ(k,1,l,1)+S2h
     &  0000321(1)*ZZ(k,1,l,2))
       aux00003(2,1,1)=-(F(1)*S2h004121(1))-F(2)*S2h004221(1)+S2h0051211
     &  (k)*Z(1,l)+S2h0052211(k)*Z(2,l)+(-2*Inv0000211+S20000321(1)-S2h0
     &  051211(1)-S2h0052221(1))*Z(k,l)-4*(S2h0000312(1)*ZZ(k,1,l,1)+S2h
     &  0000322(1)*ZZ(k,1,l,2))-2*S2h0000311(1)*ZZ(k,2,l,1)-2*S2h0000321
     &  (1)*ZZ(k,2,l,2)
       aux00003(2,2,1)=-(F(1)*S2h004122(1))-F(2)*S2h004222(1)+S2h0051221
     &  (k)*Z(1,l)+S2h0052221(k)*Z(2,l)+(-2*Inv0000221+S20000322(1)-S2h0
     &  051221(1)-S2h0052222(1))*Z(k,l)-2*(S2h0000312(2)*ZZ(k,1,l,1)+S2h
     &  0000322(2)*ZZ(k,1,l,2))-4*S2h0000312(1)*ZZ(k,2,l,1)-4*S2h0000322
     &  (1)*ZZ(k,2,l,2)
       aux00003(2,2,2)=-(F(1)*S2h004122(2))-F(2)*S2h004222(2)+S2h0051222
     &  (k)*Z(1,l)+S2h0052222(k)*Z(2,l)+(-2*Inv0000222+S20000322(2)-S2h0
     &  051222(1)-S2h0052222(2))*Z(k,l)-6*(S2h0000312(2)*ZZ(k,2,l,1)+S2h
     &  0000322(2)*ZZ(k,2,l,2))
       temp00003(1,1,1)=I22Z*(aux00003(1,1,1)+6*F(4)*temp00002(1,1)+F(3)
     &  *temp003(1,1,1)+24*temp0000001(1)*ZZ(k,1,l,1))
       temp00003(2,1,1)=I22Z*(aux00003(2,1,1)+F(5)*temp00002(1,1)+4*F(4)
     &  *temp00002(2,1)+F(3)*temp003(2,1,1)+8*(temp0000001(2)*ZZ(k,1,l,1
     &  )+temp0000001(1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))))
       temp00003(2,2,1)=I22Z*(aux00003(2,2,1)+2*F(5)*temp00002(2,1)+2*F(
     &  4)*temp00002(2,2)+F(3)*temp003(2,2,1)+8*(temp0000001(2)*(ZZ(k,1,
     &  l,2)+ZZ(k,2,l,1))+temp0000001(1)*ZZ(k,2,l,2)))
       temp00003(2,2,2)=I22Z*(aux00003(2,2,2)+F(3)*temp003(2,2,2)+24*tem
     &  p0000001(2)*ZZ(k,2,l,2)+temp00002(2,2)*(6*r10*(ZZ(k,1,l,2)+ZZ(k,
     &  2,l,1))+12*r21*ZZ(k,2,l,2)))
       temp00003(1,1,2)=temp00003(2,1,1)
       temp00003(1,2,1)=temp00003(2,1,1)
       temp00003(1,2,2)=temp00003(2,2,1)
       temp00003(2,1,2)=temp00003(2,2,1)
       S20051111(1)=-2*B023
       S20051111(2)=2*B23(1,1)
       S20051211(1)=2*B23(1,1)
       S20051211(2)=-2*B23(1,2)
       S20051221(1)=-2*B23(1,2)
       S20051221(2)=2*B23(1,3)
       S20051222(1)=2*B23(1,3)
       S20051222(2)=-2*B23(1,4)
       S20052111(1)=2*B23(1,1)
       S20052111(2)=-2*B23(1,2)
       S20052211(1)=-2*B23(1,2)
       S20052211(2)=2*B23(1,3)
       S20052221(1)=2*B23(1,3)
       S20052221(2)=-2*B23(1,4)
       S20052222(1)=-2*B23(1,4)
       S20052222(2)=2*B23(1,5)
       S27111111(1)=-B023+B13(1,6)
       S27111111(2)=B13(1,6)+B23(1,1)
       S27121111(1)=B13(1,6)+B23(1,1)
       S27121111(2)=B13(1,6)-B23(1,2)
       S27122111(1)=B13(1,6)-B23(1,2)
       S27122111(2)=B13(1,6)+B23(1,3)
       S27122211(1)=B13(1,6)+B23(1,3)
       S27122211(2)=B13(1,6)-B23(1,4)
       S27122221(1)=B13(1,6)-B23(1,4)
       S27122221(2)=B13(1,6)+B23(1,5)
       S27122222(1)=B13(1,6)+B23(1,5)
       S27122222(2)=B13(1,6)-B23(1,6)
       S27211111(1)=B12(1,6)-B13(1,6)
       S27211111(2)=-B13(1,6)
       S27221111(1)=-B13(1,6)
       S27221111(2)=-B13(1,6)
       S27222111(1)=-B13(1,6)
       S27222111(2)=-B13(1,6)
       S27222211(1)=-B13(1,6)
       S27222211(2)=-B13(1,6)
       S27222221(1)=-B13(1,6)
       S27222221(2)=-B13(1,6)
       S27222222(1)=-B13(1,6)
       S27222222(2)=-B13(1,6)
       aux00511(1,1,1)=-(F(1)*S2611111(1))-F(2)*S2621111(1)+S27111111(k)
     &  *Z(1,l)+S27211111(k)*Z(2,l)+(-2*Inv0011111+S20051111(1)-S2711111
     &  1(1)-S27221111(1))*Z(k,l)-10*(S2h0051111(1)*ZZ(k,1,l,1)+S2h00521
     &  11(1)*ZZ(k,1,l,2))
       aux00521(1,1,1)=-(F(1)*S2612111(1))-F(2)*S2622111(1)+S27121111(k)
     &  *Z(1,l)+S27221111(k)*Z(2,l)+(-2*Inv0021111+S20052111(1)-S2712111
     &  1(1)-S27222111(1))*Z(k,l)-8*(S2h0051211(1)*ZZ(k,1,l,1)+S2h005221
     &  1(1)*ZZ(k,1,l,2))-2*S2h0051111(1)*ZZ(k,2,l,1)-2*S2h0052111(1)*ZZ
     &  (k,2,l,2)
       aux00522(1,1,1)=-(F(1)*S2612211(1))-F(2)*S2622211(1)+S27122111(k)
     &  *Z(1,l)+S27222111(k)*Z(2,l)+(-2*Inv0022111+S20052211(1)-S2712211
     &  1(1)-S27222211(1))*Z(k,l)-6*(S2h0051221(1)*ZZ(k,1,l,1)+S2h005222
     &  1(1)*ZZ(k,1,l,2))-4*S2h0051211(1)*ZZ(k,2,l,1)-4*S2h0052211(1)*ZZ
     &  (k,2,l,2)
       aux00522(2,1,1)=-(F(1)*S2612221(1))-F(2)*S2622221(1)+S27122211(k)
     &  *Z(1,l)+S27222211(k)*Z(2,l)+(-2*Inv0022211+S20052221(1)-S2712221
     &  1(1)-S27222221(1))*Z(k,l)-4*(S2h0051222(1)*ZZ(k,1,l,1)+S2h005222
     &  2(1)*ZZ(k,1,l,2))-6*S2h0051221(1)*ZZ(k,2,l,1)-6*S2h0052221(1)*ZZ
     &  (k,2,l,2)
       aux00522(2,2,1)=-(F(1)*S2612222(1))-F(2)*S2622222(1)+S27122221(k)
     &  *Z(1,l)+S27222221(k)*Z(2,l)+(-2*Inv0022221+S20052222(1)-S2712222
     &  1(1)-S27222222(1))*Z(k,l)-2*(S2h0051222(2)*ZZ(k,1,l,1)+S2h005222
     &  2(2)*ZZ(k,1,l,2))-8*S2h0051222(1)*ZZ(k,2,l,1)-8*S2h0052222(1)*ZZ
     &  (k,2,l,2)
       aux00522(2,2,2)=-(F(1)*S2612222(2))-F(2)*S2622222(2)+S27122222(k)
     &  *Z(1,l)+S27222222(k)*Z(2,l)+(-2*Inv0022222+S20052222(2)-S2712222
     &  2(1)-S27222222(2))*Z(k,l)-10*(S2h0051222(2)*ZZ(k,2,l,1)+S2h00522
     &  22(2)*ZZ(k,2,l,2))
       temp00511(1,1,1)=I26Z*(aux00511(1,1,1)+10*F(4)*temp0041(1,1,1)+F(
     &  3)*temp511(1,1,1)+80*temp00003(1,1,1)*ZZ(k,1,l,1))
       temp00521(1,1,1)=I26Z*(aux00521(1,1,1)+F(5)*temp0041(1,1,1)+8*F(4
     &  )*temp0042(1,1,1)+F(3)*temp521(1,1,1)+48*temp00003(2,1,1)*ZZ(k,1
     &  ,l,1)+16*temp00003(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp00522(1,1,1)=I26Z*(aux00522(1,1,1)+2*F(5)*temp0042(1,1,1)+6*F
     &  (4)*temp0042(2,1,1)+F(3)*temp522(1,1,1)+24*(temp00003(2,2,1)*ZZ(
     &  k,1,l,1)+temp00003(2,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))+8*temp00003
     &  (1,1,1)*ZZ(k,2,l,2))
       temp00522(2,1,1)=I26Z*(aux00522(2,1,1)+4*F(4)*temp0042(2,2,1)+F(3
     &  )*temp522(2,1,1)+8*temp00003(2,2,2)*ZZ(k,1,l,1)+temp0042(2,1,1)*
     &  (6*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2))+24*(temp000
     &  03(2,2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp00003(2,1,1)*ZZ(k,2,l,2)
     &  ))
       temp00522(2,2,1)=I26Z*(aux00522(2,2,1)+4*F(5)*temp0042(2,2,1)+2*F
     &  (4)*temp0042(2,2,2)+F(3)*temp522(2,2,1)+16*temp00003(2,2,2)*(ZZ(
     &  k,1,l,2)+ZZ(k,2,l,1))+48*temp00003(2,2,1)*ZZ(k,2,l,2))
       temp00522(2,2,2)=I26Z*(aux00522(2,2,2)+F(3)*temp522(2,2,2)+80*tem
     &  p00003(2,2,2)*ZZ(k,2,l,2)+temp0042(2,2,2)*(10*r10*(ZZ(k,1,l,2)+Z
     &  Z(k,2,l,1))+20*r21*ZZ(k,2,l,2)))
       temp00511(1,1,2)=temp00521(1,1,1)
       temp00511(1,2,1)=temp00521(1,1,1)
       temp00511(1,2,2)=temp00522(1,1,1)
       temp00521(1,1,2)=temp00522(1,1,1)
       temp00521(1,2,1)=temp00522(1,1,1)
       temp00521(1,2,2)=temp00522(2,1,1)
       temp00522(1,1,2)=temp00522(2,1,1)
       temp00522(1,2,1)=temp00522(2,1,1)
       temp00522(1,2,2)=temp00522(2,2,1)
       temp00522(2,1,2)=temp00522(2,2,1)
       aux6111(1,1,1)=-(S27111111(1)*Z(jj,1))-S27211111(1)*Z(jj,2)
       aux6211(1,1,1)=-(S27121111(1)*Z(jj,1))-S27221111(1)*Z(jj,2)
       aux6221(1,1,1)=-(S27122111(1)*Z(jj,1))-S27222111(1)*Z(jj,2)
       aux6222(1,1,1)=-(S27122211(1)*Z(jj,1))-S27222211(1)*Z(jj,2)
       aux6222(2,1,1)=-(S27122221(1)*Z(jj,1))-S27222221(1)*Z(jj,2)
       aux6222(2,2,1)=-(S27122222(1)*Z(jj,1))-S27222222(1)*Z(jj,2)
       aux6222(2,2,2)=-(S27122222(2)*Z(jj,1))-S27222222(2)*Z(jj,2)
       temp6111(1,1,1)=IX*(aux6111(1,1,1)+12*temp00511(1,1,1)*Z(jj,1))
       temp6211(1,1,1)=IX*(aux6211(1,1,1)+10*temp00521(1,1,1)*Z(jj,1)+2*
     &  temp00511(1,1,1)*Z(jj,2))
       temp6221(1,1,1)=IX*(aux6221(1,1,1)+8*temp00522(1,1,1)*Z(jj,1)+4*t
     &  emp00521(1,1,1)*Z(jj,2))
       temp6222(1,1,1)=IX*(aux6222(1,1,1)+6*(temp00522(2,1,1)*Z(jj,1)+te
     &  mp00522(1,1,1)*Z(jj,2)))
       temp6222(2,1,1)=IX*(aux6222(2,1,1)+4*temp00522(2,2,1)*Z(jj,1)+8*t
     &  emp00522(2,1,1)*Z(jj,2))
       temp6222(2,2,1)=IX*(aux6222(2,2,1)+2*temp00522(2,2,2)*Z(jj,1)+10*
     &  temp00522(2,2,1)*Z(jj,2))
       temp6222(2,2,2)=IX*(aux6222(2,2,2)+12*temp00522(2,2,2)*Z(jj,2))
       temp6111(1,1,2)=temp6211(1,1,1)
       temp6111(1,2,1)=temp6211(1,1,1)
       temp6111(1,2,2)=temp6221(1,1,1)
       temp6211(1,1,2)=temp6221(1,1,1)
       temp6211(1,2,1)=temp6221(1,1,1)
       temp6211(1,2,2)=temp6222(1,1,1)
       temp6221(1,1,2)=temp6222(1,1,1)
       temp6221(1,2,1)=temp6222(1,1,1)
       temp6221(1,2,2)=temp6222(2,1,1)
       temp6222(1,1,2)=temp6222(2,1,1)
       temp6222(1,2,1)=temp6222(2,1,1)
       temp6222(1,2,2)=temp6222(2,2,1)
       temp6222(2,1,2)=temp6222(2,2,1)
c                Step2
       tempC3000000=I14Z*(auxC3000000+tempC30000*F(3)-det3*temp00002(k,l
     &  ))
       temp00002(1,1)=I18Z*(aux00002(1,1)+4*F(4)*temp00001(1)+F(3)*temp0
     &  02(1,1)-det3*temp0041(1,k,l)+8*tempC3000000*ZZ(k,1,l,1))
       temp00002(2,1)=I18Z*(aux00002(2,1)+F(5)*temp00001(1)+2*F(4)*temp0
     &  0001(2)+F(3)*temp002(2,1)-det3*temp0042(1,k,l)+4*tempC3000000*(Z
     &  Z(k,1,l,2)+ZZ(k,2,l,1)))
       temp00002(2,2)=I18Z*(aux00002(2,2)+2*F(5)*temp00001(2)+F(3)*temp0
     &  02(2,2)-det3*temp0042(2,k,l)+8*tempC3000000*ZZ(k,2,l,2))
       temp00002(1,2)=temp00002(2,1)
       temp0041(1,1,1)=I22Z*(aux0041(1,1,1)+8*F(4)*temp003(1,1,1)+F(3)*t
     &  emp41(1,1,1)-det3*temp6111(1,k,l)+48*temp00002(1,1)*ZZ(k,1,l,1))
       temp0042(1,1,1)=I22Z*(aux0042(1,1,1)+F(5)*temp003(1,1,1)+6*F(4)*t
     &  emp003(2,1,1)+F(3)*temp42(1,1,1)-det3*temp6211(1,k,l)+24*temp000
     &  02(2,1)*ZZ(k,1,l,1)+12*temp00002(1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp0042(2,1,1)=I22Z*(aux0042(2,1,1)+2*F(5)*temp003(2,1,1)+4*F(4)
     &  *temp003(2,2,1)+F(3)*temp42(2,1,1)-det3*temp6221(1,k,l)+16*temp0
     &  0002(2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+8*(temp00002(2,2)*ZZ(k,1,l,1
     &  )+temp00002(1,1)*ZZ(k,2,l,2)))
       temp0042(2,2,1)=I22Z*(aux0042(2,2,1)+2*F(4)*temp003(2,2,2)+F(3)*t
     &  emp42(2,2,1)-det3*temp6222(1,k,l)+12*temp00002(2,2)*(ZZ(k,1,l,2)
     &  +ZZ(k,2,l,1))+24*temp00002(2,1)*ZZ(k,2,l,2)+temp003(2,2,1)*(6*r1
     &  0*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2)))
       temp0042(2,2,2)=I22Z*(aux0042(2,2,2)+4*F(5)*temp003(2,2,2)+F(3)*t
     &  emp42(2,2,2)-det3*temp6222(2,k,l)+48*temp00002(2,2)*ZZ(k,2,l,2))
       temp0041(1,1,2)=temp0042(1,1,1)
       temp0041(1,2,1)=temp0042(1,1,1)
       temp0041(1,2,2)=temp0042(2,1,1)
       temp0042(1,1,2)=temp0042(2,1,1)
       temp0042(1,2,1)=temp0042(2,1,1)
       temp0042(1,2,2)=temp0042(2,2,1)
       temp0042(2,1,2)=temp0042(2,2,1)
       temp511(1,1,1)=IX*(aux511(1,1,1)+det3*temp6111(1,1,jj)+10*temp004
     &  1(1,1,1)*Z(jj,1))
       temp521(1,1,1)=IX*(aux521(1,1,1)+det3*temp6211(1,1,jj)+8*temp0042
     &  (1,1,1)*Z(jj,1)+2*temp0041(1,1,1)*Z(jj,2))
       temp522(1,1,1)=IX*(aux522(1,1,1)+det3*temp6221(1,1,jj)+6*temp0042
     &  (2,1,1)*Z(jj,1)+4*temp0042(1,1,1)*Z(jj,2))
       temp522(2,1,1)=IX*(aux522(2,1,1)+det3*temp6222(1,1,jj)+4*temp0042
     &  (2,2,1)*Z(jj,1)+6*temp0042(2,1,1)*Z(jj,2))
       temp522(2,2,1)=IX*(aux522(2,2,1)+det3*temp6222(2,1,jj)+2*temp0042
     &  (2,2,2)*Z(jj,1)+8*temp0042(2,2,1)*Z(jj,2))
       temp522(2,2,2)=IX*(aux522(2,2,2)+det3*temp6222(2,2,jj)+10*temp004
     &  2(2,2,2)*Z(jj,2))
       temp511(1,1,2)=temp521(1,1,1)
       temp511(1,2,1)=temp521(1,1,1)
       temp511(1,2,2)=temp522(1,1,1)
       temp521(1,1,2)=temp522(1,1,1)
       temp521(1,2,1)=temp522(1,1,1)
       temp521(1,2,2)=temp522(2,1,1)
       temp522(1,1,2)=temp522(2,1,1)
       temp522(1,2,1)=temp522(2,1,1)
       temp522(1,2,2)=temp522(2,2,1)
       temp522(2,1,2)=temp522(2,2,1)
c                Step3
       temp00001(1)=I14Z*(aux00001(1)+2*tempC30000*F(4)+F(3)*temp001(1)-
     &  det3*temp003(1,k,l))
       temp00001(2)=I14Z*(aux00001(2)+tempC30000*F(5)+F(3)*temp001(2)-de
     &  t3*temp003(2,k,l))
       temp003(1,1,1)=I18Z*(aux003(1,1,1)+6*F(4)*temp002(1,1)+F(3)*temp3
     &  (1,1,1)-det3*temp511(1,k,l)+24*temp00001(1)*ZZ(k,1,l,1))
       temp003(2,1,1)=I18Z*(aux003(2,1,1)+F(5)*temp002(1,1)+4*F(4)*temp0
     &  02(2,1)+F(3)*temp3(2,1,1)-det3*temp521(1,k,l)+8*(temp00001(2)*ZZ
     &  (k,1,l,1)+temp00001(1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))))
       temp003(2,2,1)=I18Z*(aux003(2,2,1)+2*F(5)*temp002(2,1)+2*F(4)*tem
     &  p002(2,2)+F(3)*temp3(2,2,1)-det3*temp522(1,k,l)+8*(temp00001(2)*
     &  (ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp00001(1)*ZZ(k,2,l,2)))
       temp003(2,2,2)=I18Z*(aux003(2,2,2)+F(3)*temp3(2,2,2)-det3*temp522
     &  (2,k,l)+24*temp00001(2)*ZZ(k,2,l,2)+temp002(2,2)*(6*r10*(ZZ(k,1,
     &  l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2)))
       temp003(1,1,2)=temp003(2,1,1)
       temp003(1,2,1)=temp003(2,1,1)
       temp003(1,2,2)=temp003(2,2,1)
       temp003(2,1,2)=temp003(2,2,1)
       temp41(1,1,1)=IX*(aux41(1,1,1)+det3*temp511(1,1,jj)+8*temp003(1,1
     &  ,1)*Z(jj,1))
       temp42(1,1,1)=IX*(aux42(1,1,1)+det3*temp521(1,1,jj)+6*temp003(2,1
     &  ,1)*Z(jj,1)+2*temp003(1,1,1)*Z(jj,2))
       temp42(2,1,1)=IX*(aux42(2,1,1)+det3*temp522(1,1,jj)+4*(temp003(2,
     &  2,1)*Z(jj,1)+temp003(2,1,1)*Z(jj,2)))
       temp42(2,2,1)=IX*(aux42(2,2,1)+det3*temp522(2,1,jj)+2*temp003(2,2
     &  ,2)*Z(jj,1)+6*temp003(2,2,1)*Z(jj,2))
       temp42(2,2,2)=IX*(aux42(2,2,2)+det3*temp522(2,2,jj)+8*temp003(2,2
     &  ,2)*Z(jj,2))
       temp41(1,1,2)=temp42(1,1,1)
       temp41(1,2,1)=temp42(1,1,1)
       temp41(1,2,2)=temp42(2,1,1)
       temp42(1,1,2)=temp42(2,1,1)
       temp42(1,2,1)=temp42(2,1,1)
       temp42(1,2,2)=temp42(2,2,1)
       temp42(2,1,2)=temp42(2,2,1)
c                Step4
       tempC30000=I10Z*(auxC30000+tempC300*F(3)-det3*temp002(k,l))
       temp002(1,1)=I14Z*(aux002(1,1)+4*F(4)*temp001(1)+F(3)*temp2(1,1)-
     &  det3*temp41(1,k,l)+8*tempC30000*ZZ(k,1,l,1))
       temp002(2,1)=I14Z*(aux002(2,1)+F(5)*temp001(1)+2*F(4)*temp001(2)+
     &  F(3)*temp2(2,1)-det3*temp42(1,k,l)+4*tempC30000*(ZZ(k,1,l,2)+ZZ(
     &  k,2,l,1)))
       temp002(2,2)=I14Z*(aux002(2,2)+2*F(5)*temp001(2)+F(3)*temp2(2,2)-
     &  det3*temp42(2,k,l)+8*tempC30000*ZZ(k,2,l,2))
       temp002(1,2)=temp002(2,1)
       temp3(1,1,1)=IX*(aux3(1,1,1)+det3*temp41(1,1,jj)+6*temp002(1,1)*Z
     &  (jj,1))
       temp3(2,1,1)=IX*(aux3(2,1,1)+det3*temp42(1,1,jj)+4*temp002(2,1)*Z
     &  (jj,1)+2*temp002(1,1)*Z(jj,2))
       temp3(2,2,1)=IX*(aux3(2,2,1)+det3*temp42(2,1,jj)+2*temp002(2,2)*Z
     &  (jj,1)+4*temp002(2,1)*Z(jj,2))
       temp3(2,2,2)=IX*(aux3(2,2,2)+det3*temp42(2,2,jj)+6*temp002(2,2)*Z
     &  (jj,2))
       temp3(1,1,2)=temp3(2,1,1)
       temp3(1,2,1)=temp3(2,1,1)
       temp3(1,2,2)=temp3(2,2,1)
       temp3(2,1,2)=temp3(2,2,1)
c                Step5
       temp001(1)=I10Z*(aux001(1)+2*tempC300*F(4)+F(3)*temp1(1)-det3*tem
     &  p3(1,k,l))
       temp001(2)=I10Z*(aux001(2)+tempC300*F(5)+F(3)*temp1(2)-det3*temp3
     &  (2,k,l))
       temp2(1,1)=IX*(aux2(1,1)+det3*temp3(1,1,jj)+4*temp001(1)*Z(jj,1))
       temp2(2,1)=IX*(aux2(2,1)+det3*temp3(2,1,jj)+2*(temp001(2)*Z(jj,1)
     &  +temp001(1)*Z(jj,2)))
       temp2(2,2)=IX*(aux2(2,2)+det3*temp3(2,2,jj)+4*temp001(2)*Z(jj,2))
       temp2(1,2)=temp2(2,1)
c                Step6
       tempC300=I6Z*(auxC300+tempC30*F(3)-det3*temp2(k,l))
       temp1(1)=IX*(aux1(1)+det3*temp2(1,jj)+2*tempC300*Z(jj,1))
       temp1(2)=IX*(aux1(2)+det3*temp2(2,jj)+2*tempC300*Z(jj,2))
c                Step7
       tempC30=IX*(auxC30+det3*temp1(jj))
c       Cij(1,4)= temp41(1,1,1)
           C30=tempC30
           Cij(1,1)=temp1(1)
           Cij(2,1)=temp1(2)
           Cij(1,2)=temp2(1,1)
           Cij(2,2)=temp2(2,2)
           Cij(3,2)=temp2(2,1)
           Cij(4,2)=tempC300
           Cij(1,3)=temp3(1,1,1)
           Cij(2,3)=temp3(2,2,2)
           Cij(3,3)=temp3(2,1,1)
           Cij(4,3)=temp3(2,2,1)
           Cij(5,3)=temp001(1)
           Cij(6,3)=temp001(2)
           Cij(1,4)=temp41(1,1,1)
           Cij(2,4)=temp42(2,2,2)
           Cij(3,4)=temp42(1,1,1)
           Cij(4,4)=temp42(2,1,1)
           Cij(5,4)=temp42(2,2,1)
           Cij(6,4)=temp002(1,1)
           Cij(7,4)=temp002(2,2)
           Cij(8,4)=temp002(2,1)
           Cij(9,4)=tempC30000

       if(order.eq.6) goto 500
c                Iteration7
c                Step1
       S200000000=2*B23(4,6)
       S200000021(1)=2*B23(3,4)
       S200000021(2)=-2*B23(3,5)
       S200000022(1)=-2*B23(3,5)
       S200000022(2)=2*B23(3,6)
       S2h00000021(1)=B13(4,7)+B23(4,6)
       S2h00000021(2)=B13(4,7)-B23(4,7)
       S2h00000022(1)=B12(4,7)-B13(4,7)
       S2h00000022(2)=-B13(4,7)
       auxC300000000=-(F(1)*S2h0000001(1))-F(2)*S2h0000001(2)+S2h0000002
     &  1(k)*Z(1,l)+S2h00000022(k)*Z(2,l)+(-2*Inv00000000+S200000000-S2h
     &  00000021(1)-S2h00000022(2))*Z(k,l)
       tempC300000000=I18Z*(auxC300000000+tempC3000000*F(3))
       S2h00004111(1)=B13(3,7)+B23(3,4)
       S2h00004111(2)=B13(3,7)-B23(3,5)
       S2h00004121(1)=B13(3,7)-B23(3,5)
       S2h00004121(2)=B13(3,7)+B23(3,6)
       S2h00004122(1)=B13(3,7)+B23(3,6)
       S2h00004122(2)=B13(3,7)-B23(3,7)
       S2h00004211(1)=B12(3,7)-B13(3,7)
       S2h00004211(2)=-B13(3,7)
       S2h00004221(1)=-B13(3,7)
       S2h00004221(2)=-B13(3,7)
       S2h00004222(1)=-B13(3,7)
       S2h00004222(2)=-B13(3,7)
       aux0000002(1,1)=-(F(1)*S2h0000311(1))-F(2)*S2h0000321(1)+S2h00004
     &  111(k)*Z(1,l)+S2h00004211(k)*Z(2,l)+(-2*Inv00000011+S200000021(1
     &  )-S2h00004111(1)-S2h00004221(1))*Z(k,l)-4*(S2h00000021(1)*ZZ(k,1
     &  ,l,1)+S2h00000022(1)*ZZ(k,1,l,2))
       aux0000002(2,1)=-(F(1)*S2h0000312(1))-F(2)*S2h0000322(1)+S2h00004
     &  121(k)*Z(1,l)+S2h00004221(k)*Z(2,l)+(-2*Inv00000021+S200000022(1
     &  )-S2h00004121(1)-S2h00004222(1))*Z(k,l)-2*(S2h00000021(2)*ZZ(k,1
     &  ,l,1)+S2h00000022(2)*ZZ(k,1,l,2))-2*S2h00000021(1)*ZZ(k,2,l,1)-2
     &  *S2h00000022(1)*ZZ(k,2,l,2)
       aux0000002(2,2)=-(F(1)*S2h0000312(2))-F(2)*S2h0000322(2)+S2h00004
     &  122(k)*Z(1,l)+S2h00004222(k)*Z(2,l)+(-2*Inv00000022+S200000022(2
     &  )-S2h00004122(1)-S2h00004222(2))*Z(k,l)-4*(S2h00000021(2)*ZZ(k,2
     &  ,l,1)+S2h00000022(2)*ZZ(k,2,l,2))
       temp0000002(1,1)=I22Z*(aux0000002(1,1)+4*F(4)*temp0000001(1)+F(3)
     &  *temp00002(1,1)+8*tempC300000000*ZZ(k,1,l,1))
       temp0000002(2,1)=I22Z*(aux0000002(2,1)+F(5)*temp0000001(1)+2*F(4)
     &  *temp0000001(2)+F(3)*temp00002(2,1)+4*tempC300000000*(ZZ(k,1,l,2
     &  )+ZZ(k,2,l,1)))
       temp0000002(2,2)=I22Z*(aux0000002(2,2)+2*F(5)*temp0000001(2)+F(3)
     &  *temp00002(2,2)+8*tempC300000000*ZZ(k,2,l,2))
       temp0000002(1,2)=temp0000002(2,1)
       S200004111(1)=2*B23(2,2)
       S200004111(2)=-2*B23(2,3)
       S200004121(1)=-2*B23(2,3)
       S200004121(2)=2*B23(2,4)
       S200004122(1)=2*B23(2,4)
       S200004122(2)=-2*B23(2,5)
       S200004211(1)=-2*B23(2,3)
       S200004211(2)=2*B23(2,4)
       S200004221(1)=2*B23(2,4)
       S200004221(2)=-2*B23(2,5)
       S200004222(1)=-2*B23(2,5)
       S200004222(2)=2*B23(2,6)
       S200611111(1)=2*B023
       S200611111(2)=-2*B23(1,1)
       S200612111(1)=-2*B23(1,1)
       S200612111(2)=2*B23(1,2)
       S200612211(1)=2*B23(1,2)
       S200612211(2)=-2*B23(1,3)
       S200612221(1)=-2*B23(1,3)
       S200612221(2)=2*B23(1,4)
       S200612222(1)=2*B23(1,4)
       S200612222(2)=-2*B23(1,5)
       S200621111(1)=-2*B23(1,1)
       S200621111(2)=2*B23(1,2)
       S200622111(1)=2*B23(1,2)
       S200622111(2)=-2*B23(1,3)
       S200622211(1)=-2*B23(1,3)
       S200622211(2)=2*B23(1,4)
       S200622221(1)=2*B23(1,4)
       S200622221(2)=-2*B23(1,5)
       S200622222(1)=-2*B23(1,5)
       S200622222(2)=2*B23(1,6)
       S2h00611111(1)=B13(2,7)+B23(2,2)
       S2h00611111(2)=B13(2,7)-B23(2,3)
       S2h00612111(1)=B13(2,7)-B23(2,3)
       S2h00612111(2)=B13(2,7)+B23(2,4)
       S2h00612211(1)=B13(2,7)+B23(2,4)
       S2h00612211(2)=B13(2,7)-B23(2,5)
       S2h00612221(1)=B13(2,7)-B23(2,5)
       S2h00612221(2)=B13(2,7)+B23(2,6)
       S2h00612222(1)=B13(2,7)+B23(2,6)
       S2h00612222(2)=B13(2,7)-B23(2,7)
       S2h00621111(1)=B12(2,7)-B13(2,7)
       S2h00621111(2)=-B13(2,7)
       S2h00622111(1)=-B13(2,7)
       S2h00622111(2)=-B13(2,7)
       S2h00622211(1)=-B13(2,7)
       S2h00622211(2)=-B13(2,7)
       S2h00622221(1)=-B13(2,7)
       S2h00622221(2)=-B13(2,7)
       S2h00622222(1)=-B13(2,7)
       S2h00622222(2)=-B13(2,7)
       aux000041(1,1,1)=-(F(1)*S2h0051111(1))-F(2)*S2h0052111(1)+S2h0061
     &  1111(k)*Z(1,l)+S2h00621111(k)*Z(2,l)+(-2*Inv00001111+S200004111(
     &  1)-S2h00611111(1)-S2h00622111(1))*Z(k,l)-8*(S2h00004111(1)*ZZ(k,
     &  1,l,1)+S2h00004211(1)*ZZ(k,1,l,2))
       aux000042(1,1,1)=-(F(1)*S2h0051211(1))-F(2)*S2h0052211(1)+S2h0061
     &  2111(k)*Z(1,l)+S2h00622111(k)*Z(2,l)+(-2*Inv00002111+S200004211(
     &  1)-S2h00612111(1)-S2h00622211(1))*Z(k,l)-6*(S2h00004121(1)*ZZ(k,
     &  1,l,1)+S2h00004221(1)*ZZ(k,1,l,2))-2*S2h00004111(1)*ZZ(k,2,l,1)-
     &  2*S2h00004211(1)*ZZ(k,2,l,2)
       aux000042(2,1,1)=-(F(1)*S2h0051221(1))-F(2)*S2h0052221(1)+S2h0061
     &  2211(k)*Z(1,l)+S2h00622211(k)*Z(2,l)+(-2*Inv00002211+S200004221(
     &  1)-S2h00612211(1)-S2h00622221(1))*Z(k,l)-4*(S2h00004122(1)*ZZ(k,
     &  1,l,1)+S2h00004222(1)*ZZ(k,1,l,2))-4*S2h00004121(1)*ZZ(k,2,l,1)-
     &  4*S2h00004221(1)*ZZ(k,2,l,2)
       aux000042(2,2,1)=-(F(1)*S2h0051222(1))-F(2)*S2h0052222(1)+S2h0061
     &  2221(k)*Z(1,l)+S2h00622221(k)*Z(2,l)+(-2*Inv00002221+S200004222(
     &  1)-S2h00612221(1)-S2h00622222(1))*Z(k,l)-2*(S2h00004122(2)*ZZ(k,
     &  1,l,1)+S2h00004222(2)*ZZ(k,1,l,2))-6*S2h00004122(1)*ZZ(k,2,l,1)-
     &  6*S2h00004222(1)*ZZ(k,2,l,2)
       aux000042(2,2,2)=-(F(1)*S2h0051222(2))-F(2)*S2h0052222(2)+S2h0061
     &  2222(k)*Z(1,l)+S2h00622222(k)*Z(2,l)+(-2*Inv00002222+S200004222(
     &  2)-S2h00612222(1)-S2h00622222(2))*Z(k,l)-8*(S2h00004122(2)*ZZ(k,
     &  2,l,1)+S2h00004222(2)*ZZ(k,2,l,2))
       temp000041(1,1,1)=I26Z*(aux000041(1,1,1)+8*F(4)*temp00003(1,1,1)+
     &  F(3)*temp0041(1,1,1)+48*temp0000002(1,1)*ZZ(k,1,l,1))
       temp000042(1,1,1)=I26Z*(aux000042(1,1,1)+F(5)*temp00003(1,1,1)+6*
     &  F(4)*temp00003(2,1,1)+F(3)*temp0042(1,1,1)+24*temp0000002(2,1)*Z
     &  Z(k,1,l,1)+12*temp0000002(1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp000042(2,1,1)=I26Z*(aux000042(2,1,1)+2*F(5)*temp00003(2,1,1)+
     &  4*F(4)*temp00003(2,2,1)+F(3)*temp0042(2,1,1)+16*temp0000002(2,1)
     &  *(ZZ(k,1,l,2)+ZZ(k,2,l,1))+8*(temp0000002(2,2)*ZZ(k,1,l,1)+temp0
     &  000002(1,1)*ZZ(k,2,l,2)))
       temp000042(2,2,1)=I26Z*(aux000042(2,2,1)+2*F(4)*temp00003(2,2,2)+
     &  F(3)*temp0042(2,2,1)+12*temp0000002(2,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1
     &  ))+24*temp0000002(2,1)*ZZ(k,2,l,2)+temp00003(2,2,1)*(6*r10*(ZZ(k
     &  ,1,l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2)))
       temp000042(2,2,2)=I26Z*(aux000042(2,2,2)+4*F(5)*temp00003(2,2,2)+
     &  F(3)*temp0042(2,2,2)+48*temp0000002(2,2)*ZZ(k,2,l,2))
       temp000041(1,1,2)=temp000042(1,1,1)
       temp000041(1,2,1)=temp000042(1,1,1)
       temp000041(1,2,2)=temp000042(2,1,1)
       temp000042(1,1,2)=temp000042(2,1,1)
       temp000042(1,2,1)=temp000042(2,1,1)
       temp000042(1,2,2)=temp000042(2,2,1)
       temp000042(2,1,2)=temp000042(2,2,1)
       S281111111(1)=B023+B13(1,7)
       S281111111(2)=B13(1,7)-B23(1,1)
       S281211111(1)=B13(1,7)-B23(1,1)
       S281211111(2)=B13(1,7)+B23(1,2)
       S281221111(1)=B13(1,7)+B23(1,2)
       S281221111(2)=B13(1,7)-B23(1,3)
       S281222111(1)=B13(1,7)-B23(1,3)
       S281222111(2)=B13(1,7)+B23(1,4)
       S281222211(1)=B13(1,7)+B23(1,4)
       S281222211(2)=B13(1,7)-B23(1,5)
       S281222221(1)=B13(1,7)-B23(1,5)
       S281222221(2)=B13(1,7)+B23(1,6)
       S281222222(1)=B13(1,7)+B23(1,6)
       S281222222(2)=B13(1,7)-B23(1,7)
       S282111111(1)=B12(1,7)-B13(1,7)
       S282111111(2)=-B13(1,7)
       S282211111(1)=-B13(1,7)
       S282211111(2)=-B13(1,7)
       S282221111(1)=-B13(1,7)
       S282221111(2)=-B13(1,7)
       S282222111(1)=-B13(1,7)
       S282222111(2)=-B13(1,7)
       S282222211(1)=-B13(1,7)
       S282222211(2)=-B13(1,7)
       S282222221(1)=-B13(1,7)
       S282222221(2)=-B13(1,7)
       S282222222(1)=-B13(1,7)
       S282222222(2)=-B13(1,7)
       aux006111(1,1,1)=-(F(1)*S27111111(1))-F(2)*S27211111(1)+S28111111
     &  1(k)*Z(1,l)+S282111111(k)*Z(2,l)+(-2*Inv00111111+S200611111(1)-S
     &  281111111(1)-S282211111(1))*Z(k,l)-12*(S2h00611111(1)*ZZ(k,1,l,1
     &  )+S2h00621111(1)*ZZ(k,1,l,2))
       aux006211(1,1,1)=-(F(1)*S27121111(1))-F(2)*S27221111(1)+S28121111
     &  1(k)*Z(1,l)+S282211111(k)*Z(2,l)+(-2*Inv00211111+S200621111(1)-S
     &  281211111(1)-S282221111(1))*Z(k,l)-10*(S2h00612111(1)*ZZ(k,1,l,1
     &  )+S2h00622111(1)*ZZ(k,1,l,2))-2*S2h00611111(1)*ZZ(k,2,l,1)-2*S2h
     &  00621111(1)*ZZ(k,2,l,2)
       aux006221(1,1,1)=-(F(1)*S27122111(1))-F(2)*S27222111(1)+S28122111
     &  1(k)*Z(1,l)+S282221111(k)*Z(2,l)+(-2*Inv00221111+S200622111(1)-S
     &  281221111(1)-S282222111(1))*Z(k,l)-8*(S2h00612211(1)*ZZ(k,1,l,1)
     &  +S2h00622211(1)*ZZ(k,1,l,2))-4*S2h00612111(1)*ZZ(k,2,l,1)-4*S2h0
     &  0622111(1)*ZZ(k,2,l,2)
       aux006222(1,1,1)=-(F(1)*S27122211(1))-F(2)*S27222211(1)+S28122211
     &  1(k)*Z(1,l)+S282222111(k)*Z(2,l)+(-2*Inv00222111+S200622211(1)-S
     &  281222111(1)-S282222211(1))*Z(k,l)-6*(S2h00612221(1)*ZZ(k,1,l,1)
     &  +S2h00622221(1)*ZZ(k,1,l,2))-6*S2h00612211(1)*ZZ(k,2,l,1)-6*S2h0
     &  0622211(1)*ZZ(k,2,l,2)
       aux006222(2,1,1)=-(F(1)*S27122221(1))-F(2)*S27222221(1)+S28122221
     &  1(k)*Z(1,l)+S282222211(k)*Z(2,l)+(-2*Inv00222211+S200622221(1)-S
     &  281222211(1)-S282222221(1))*Z(k,l)-4*(S2h00612222(1)*ZZ(k,1,l,1)
     &  +S2h00622222(1)*ZZ(k,1,l,2))-8*S2h00612221(1)*ZZ(k,2,l,1)-8*S2h0
     &  0622221(1)*ZZ(k,2,l,2)
       aux006222(2,2,1)=-(F(1)*S27122222(1))-F(2)*S27222222(1)+S28122222
     &  1(k)*Z(1,l)+S282222221(k)*Z(2,l)+(-2*Inv00222221+S200622222(1)-S
     &  281222221(1)-S282222222(1))*Z(k,l)-2*(S2h00612222(2)*ZZ(k,1,l,1)
     &  +S2h00622222(2)*ZZ(k,1,l,2))-10*S2h00612222(1)*ZZ(k,2,l,1)-10*S2
     &  h00622222(1)*ZZ(k,2,l,2)
       aux006222(2,2,2)=-(F(1)*S27122222(2))-F(2)*S27222222(2)+S28122222
     &  2(k)*Z(1,l)+S282222222(k)*Z(2,l)+(-2*Inv00222222+S200622222(2)-S
     &  281222222(1)-S282222222(2))*Z(k,l)-12*(S2h00612222(2)*ZZ(k,2,l,1
     &  )+S2h00622222(2)*ZZ(k,2,l,2))
       temp006111(1,1,1)=I30Z*(aux006111(1,1,1)+12*F(4)*temp00511(1,1,1)
     &  +F(3)*temp6111(1,1,1)+120*temp000041(1,1,1)*ZZ(k,1,l,1))
       temp006211(1,1,1)=I30Z*(aux006211(1,1,1)+F(5)*temp00511(1,1,1)+10
     &  *F(4)*temp00521(1,1,1)+F(3)*temp6211(1,1,1)+80*temp000042(1,1,1)
     &  *ZZ(k,1,l,1)+20*temp000041(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp006221(1,1,1)=I30Z*(aux006221(1,1,1)+2*F(5)*temp00521(1,1,1)+
     &  F(3)*temp6221(1,1,1)+48*temp000042(2,1,1)*ZZ(k,1,l,1)+32*temp000
     &  042(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+8*(F(4)*temp00522(1,1,1)+te
     &  mp000041(1,1,1)*ZZ(k,2,l,2)))
       temp006222(1,1,1)=I30Z*(aux006222(1,1,1)+6*F(4)*temp00522(2,1,1)+
     &  F(3)*temp6222(1,1,1)+36*temp000042(2,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,
     &  1))+temp00522(1,1,1)*(6*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+12*r21*ZZ(
     &  k,2,l,2))+24*(temp000042(2,2,1)*ZZ(k,1,l,1)+temp000042(1,1,1)*ZZ
     &  (k,2,l,2)))
       temp006222(2,1,1)=I30Z*(aux006222(2,1,1)+4*F(5)*temp00522(2,1,1)+
     &  4*F(4)*temp00522(2,2,1)+F(3)*temp6222(2,1,1)+8*temp000042(2,2,2)
     &  *ZZ(k,1,l,1)+32*temp000042(2,2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+48*t
     &  emp000042(2,1,1)*ZZ(k,2,l,2))
       temp006222(2,2,1)=I30Z*(aux006222(2,2,1)+2*F(4)*temp00522(2,2,2)+
     &  F(3)*temp6222(2,2,1)+20*temp000042(2,2,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,
     &  1))+80*temp000042(2,2,1)*ZZ(k,2,l,2)+temp00522(2,2,1)*(10*r10*(Z
     &  Z(k,1,l,2)+ZZ(k,2,l,1))+20*r21*ZZ(k,2,l,2)))
       temp006222(2,2,2)=I30Z*(aux006222(2,2,2)+F(3)*temp6222(2,2,2)+120
     &  *temp000042(2,2,2)*ZZ(k,2,l,2)+temp00522(2,2,2)*(12*r10*(ZZ(k,1,
     &  l,2)+ZZ(k,2,l,1))+24*r21*ZZ(k,2,l,2)))
       temp006111(1,1,2)=temp006211(1,1,1)
       temp006111(1,2,1)=temp006211(1,1,1)
       temp006111(1,2,2)=temp006221(1,1,1)
       temp006211(1,1,2)=temp006221(1,1,1)
       temp006211(1,2,1)=temp006221(1,1,1)
       temp006211(1,2,2)=temp006222(1,1,1)
       temp006221(1,1,2)=temp006222(1,1,1)
       temp006221(1,2,1)=temp006222(1,1,1)
       temp006221(1,2,2)=temp006222(2,1,1)
       temp006222(1,1,2)=temp006222(2,1,1)
       temp006222(1,2,1)=temp006222(2,1,1)
       temp006222(1,2,2)=temp006222(2,2,1)
       temp006222(2,1,2)=temp006222(2,2,1)
       aux71111(1,1,1)=-(S281111111(1)*Z(jj,1))-S282111111(1)*Z(jj,2)
       aux72111(1,1,1)=-(S281211111(1)*Z(jj,1))-S282211111(1)*Z(jj,2)
       aux72211(1,1,1)=-(S281221111(1)*Z(jj,1))-S282221111(1)*Z(jj,2)
       aux72221(1,1,1)=-(S281222111(1)*Z(jj,1))-S282222111(1)*Z(jj,2)
       aux72222(1,1,1)=-(S281222211(1)*Z(jj,1))-S282222211(1)*Z(jj,2)
       aux72222(2,1,1)=-(S281222221(1)*Z(jj,1))-S282222221(1)*Z(jj,2)
       aux72222(2,2,1)=-(S281222222(1)*Z(jj,1))-S282222222(1)*Z(jj,2)
       aux72222(2,2,2)=-(S281222222(2)*Z(jj,1))-S282222222(2)*Z(jj,2)
       temp71111(1,1,1)=IX*(aux71111(1,1,1)+14*temp006111(1,1,1)*Z(jj,1)
     &  )
       temp72111(1,1,1)=IX*(aux72111(1,1,1)+12*temp006211(1,1,1)*Z(jj,1)
     &  +2*temp006111(1,1,1)*Z(jj,2))
       temp72211(1,1,1)=IX*(aux72211(1,1,1)+10*temp006221(1,1,1)*Z(jj,1)
     &  +4*temp006211(1,1,1)*Z(jj,2))
       temp72221(1,1,1)=IX*(aux72221(1,1,1)+8*temp006222(1,1,1)*Z(jj,1)+
     &  6*temp006221(1,1,1)*Z(jj,2))
       temp72222(1,1,1)=IX*(aux72222(1,1,1)+6*temp006222(2,1,1)*Z(jj,1)+
     &  8*temp006222(1,1,1)*Z(jj,2))
       temp72222(2,1,1)=IX*(aux72222(2,1,1)+4*temp006222(2,2,1)*Z(jj,1)+
     &  10*temp006222(2,1,1)*Z(jj,2))
       temp72222(2,2,1)=IX*(aux72222(2,2,1)+2*temp006222(2,2,2)*Z(jj,1)+
     &  12*temp006222(2,2,1)*Z(jj,2))
       temp72222(2,2,2)=IX*(aux72222(2,2,2)+14*temp006222(2,2,2)*Z(jj,2)
     &  )
       temp71111(1,1,2)=temp72111(1,1,1)
       temp71111(1,2,1)=temp72111(1,1,1)
       temp71111(1,2,2)=temp72211(1,1,1)
       temp72111(1,1,2)=temp72211(1,1,1)
       temp72111(1,2,1)=temp72211(1,1,1)
       temp72111(1,2,2)=temp72221(1,1,1)
       temp72211(1,1,2)=temp72221(1,1,1)
       temp72211(1,2,1)=temp72221(1,1,1)
       temp72211(1,2,2)=temp72222(1,1,1)
       temp72221(1,1,2)=temp72222(1,1,1)
       temp72221(1,2,1)=temp72222(1,1,1)
       temp72221(1,2,2)=temp72222(2,1,1)
       temp72222(1,1,2)=temp72222(2,1,1)
       temp72222(1,2,1)=temp72222(2,1,1)
       temp72222(1,2,2)=temp72222(2,2,1)
       temp72222(2,1,2)=temp72222(2,2,1)
c                Step2
       temp0000001(1)=I18Z*(aux0000001(1)+2*tempC3000000*F(4)+F(3)*temp0
     &  0001(1)-det3*temp00003(1,k,l))
       temp0000001(2)=I18Z*(aux0000001(2)+tempC3000000*F(5)+F(3)*temp000
     &  01(2)-det3*temp00003(2,k,l))
       temp00003(1,1,1)=I22Z*(aux00003(1,1,1)+6*F(4)*temp00002(1,1)+F(3)
     &  *temp003(1,1,1)-det3*temp00511(1,k,l)+24*temp0000001(1)*ZZ(k,1,l
     &  ,1))
       temp00003(2,1,1)=I22Z*(aux00003(2,1,1)+F(5)*temp00002(1,1)+4*F(4)
     &  *temp00002(2,1)+F(3)*temp003(2,1,1)-det3*temp00521(1,k,l)+8*(tem
     &  p0000001(2)*ZZ(k,1,l,1)+temp0000001(1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))
     &  ))
       temp00003(2,2,1)=I22Z*(aux00003(2,2,1)+2*F(5)*temp00002(2,1)+2*F(
     &  4)*temp00002(2,2)+F(3)*temp003(2,2,1)-det3*temp00522(1,k,l)+8*(t
     &  emp0000001(2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp0000001(1)*ZZ(k,2,l,
     &  2)))
       temp00003(2,2,2)=I22Z*(aux00003(2,2,2)+F(3)*temp003(2,2,2)-det3*t
     &  emp00522(2,k,l)+24*temp0000001(2)*ZZ(k,2,l,2)+temp00002(2,2)*(6*
     &  r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2)))
       temp00003(1,1,2)=temp00003(2,1,1)
       temp00003(1,2,1)=temp00003(2,1,1)
       temp00003(1,2,2)=temp00003(2,2,1)
       temp00003(2,1,2)=temp00003(2,2,1)
       temp00511(1,1,1)=I26Z*(aux00511(1,1,1)+10*F(4)*temp0041(1,1,1)+F(
     &  3)*temp511(1,1,1)-det3*temp71111(1,k,l)+80*temp00003(1,1,1)*ZZ(k
     &  ,1,l,1))
       temp00521(1,1,1)=I26Z*(aux00521(1,1,1)+F(5)*temp0041(1,1,1)+8*F(4
     &  )*temp0042(1,1,1)+F(3)*temp521(1,1,1)-det3*temp72111(1,k,l)+48*t
     &  emp00003(2,1,1)*ZZ(k,1,l,1)+16*temp00003(1,1,1)*(ZZ(k,1,l,2)+ZZ(
     &  k,2,l,1)))
       temp00522(1,1,1)=I26Z*(aux00522(1,1,1)+2*F(5)*temp0042(1,1,1)+6*F
     &  (4)*temp0042(2,1,1)+F(3)*temp522(1,1,1)-det3*temp72211(1,k,l)+24
     &  *(temp00003(2,2,1)*ZZ(k,1,l,1)+temp00003(2,1,1)*(ZZ(k,1,l,2)+ZZ(
     &  k,2,l,1)))+8*temp00003(1,1,1)*ZZ(k,2,l,2))
       temp00522(2,1,1)=I26Z*(aux00522(2,1,1)+4*F(4)*temp0042(2,2,1)+F(3
     &  )*temp522(2,1,1)-det3*temp72221(1,k,l)+8*temp00003(2,2,2)*ZZ(k,1
     &  ,l,1)+temp0042(2,1,1)*(6*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+12*r21*ZZ
     &  (k,2,l,2))+24*(temp00003(2,2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp00
     &  003(2,1,1)*ZZ(k,2,l,2)))
       temp00522(2,2,1)=I26Z*(aux00522(2,2,1)+4*F(5)*temp0042(2,2,1)+2*F
     &  (4)*temp0042(2,2,2)+F(3)*temp522(2,2,1)-det3*temp72222(1,k,l)+16
     &  *temp00003(2,2,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+48*temp00003(2,2,1)*
     &  ZZ(k,2,l,2))
       temp00522(2,2,2)=I26Z*(aux00522(2,2,2)+F(3)*temp522(2,2,2)-det3*t
     &  emp72222(2,k,l)+80*temp00003(2,2,2)*ZZ(k,2,l,2)+temp0042(2,2,2)*
     &  (10*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+20*r21*ZZ(k,2,l,2)))
       temp00511(1,1,2)=temp00521(1,1,1)
       temp00511(1,2,1)=temp00521(1,1,1)
       temp00511(1,2,2)=temp00522(1,1,1)
       temp00521(1,1,2)=temp00522(1,1,1)
       temp00521(1,2,1)=temp00522(1,1,1)
       temp00521(1,2,2)=temp00522(2,1,1)
       temp00522(1,1,2)=temp00522(2,1,1)
       temp00522(1,2,1)=temp00522(2,1,1)
       temp00522(1,2,2)=temp00522(2,2,1)
       temp00522(2,1,2)=temp00522(2,2,1)
       temp6111(1,1,1)=IX*(aux6111(1,1,1)+det3*temp71111(1,1,jj)+12*temp
     &  00511(1,1,1)*Z(jj,1))
       temp6211(1,1,1)=IX*(aux6211(1,1,1)+det3*temp72111(1,1,jj)+10*temp
     &  00521(1,1,1)*Z(jj,1)+2*temp00511(1,1,1)*Z(jj,2))
       temp6221(1,1,1)=IX*(aux6221(1,1,1)+det3*temp72211(1,1,jj)+8*temp0
     &  0522(1,1,1)*Z(jj,1)+4*temp00521(1,1,1)*Z(jj,2))
       temp6222(1,1,1)=IX*(aux6222(1,1,1)+det3*temp72221(1,1,jj)+6*(temp
     &  00522(2,1,1)*Z(jj,1)+temp00522(1,1,1)*Z(jj,2)))
       temp6222(2,1,1)=IX*(aux6222(2,1,1)+det3*temp72222(1,1,jj)+4*temp0
     &  0522(2,2,1)*Z(jj,1)+8*temp00522(2,1,1)*Z(jj,2))
       temp6222(2,2,1)=IX*(aux6222(2,2,1)+det3*temp72222(2,1,jj)+2*temp0
     &  0522(2,2,2)*Z(jj,1)+10*temp00522(2,2,1)*Z(jj,2))
       temp6222(2,2,2)=IX*(aux6222(2,2,2)+det3*temp72222(2,2,jj)+12*temp
     &  00522(2,2,2)*Z(jj,2))
       temp6111(1,1,2)=temp6211(1,1,1)
       temp6111(1,2,1)=temp6211(1,1,1)
       temp6111(1,2,2)=temp6221(1,1,1)
       temp6211(1,1,2)=temp6221(1,1,1)
       temp6211(1,2,1)=temp6221(1,1,1)
       temp6211(1,2,2)=temp6222(1,1,1)
       temp6221(1,1,2)=temp6222(1,1,1)
       temp6221(1,2,1)=temp6222(1,1,1)
       temp6221(1,2,2)=temp6222(2,1,1)
       temp6222(1,1,2)=temp6222(2,1,1)
       temp6222(1,2,1)=temp6222(2,1,1)
       temp6222(1,2,2)=temp6222(2,2,1)
       temp6222(2,1,2)=temp6222(2,2,1)
c                Step3
       tempC3000000=I14Z*(auxC3000000+tempC30000*F(3)-det3*temp00002(k,l
     &  ))
       temp00002(1,1)=I18Z*(aux00002(1,1)+4*F(4)*temp00001(1)+F(3)*temp0
     &  02(1,1)-det3*temp0041(1,k,l)+8*tempC3000000*ZZ(k,1,l,1))
       temp00002(2,1)=I18Z*(aux00002(2,1)+F(5)*temp00001(1)+2*F(4)*temp0
     &  0001(2)+F(3)*temp002(2,1)-det3*temp0042(1,k,l)+4*tempC3000000*(Z
     &  Z(k,1,l,2)+ZZ(k,2,l,1)))
       temp00002(2,2)=I18Z*(aux00002(2,2)+2*F(5)*temp00001(2)+F(3)*temp0
     &  02(2,2)-det3*temp0042(2,k,l)+8*tempC3000000*ZZ(k,2,l,2))
       temp00002(1,2)=temp00002(2,1)
       temp0041(1,1,1)=I22Z*(aux0041(1,1,1)+8*F(4)*temp003(1,1,1)+F(3)*t
     &  emp41(1,1,1)-det3*temp6111(1,k,l)+48*temp00002(1,1)*ZZ(k,1,l,1))
       temp0042(1,1,1)=I22Z*(aux0042(1,1,1)+F(5)*temp003(1,1,1)+6*F(4)*t
     &  emp003(2,1,1)+F(3)*temp42(1,1,1)-det3*temp6211(1,k,l)+24*temp000
     &  02(2,1)*ZZ(k,1,l,1)+12*temp00002(1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp0042(2,1,1)=I22Z*(aux0042(2,1,1)+2*F(5)*temp003(2,1,1)+4*F(4)
     &  *temp003(2,2,1)+F(3)*temp42(2,1,1)-det3*temp6221(1,k,l)+16*temp0
     &  0002(2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+8*(temp00002(2,2)*ZZ(k,1,l,1
     &  )+temp00002(1,1)*ZZ(k,2,l,2)))
       temp0042(2,2,1)=I22Z*(aux0042(2,2,1)+2*F(4)*temp003(2,2,2)+F(3)*t
     &  emp42(2,2,1)-det3*temp6222(1,k,l)+12*temp00002(2,2)*(ZZ(k,1,l,2)
     &  +ZZ(k,2,l,1))+24*temp00002(2,1)*ZZ(k,2,l,2)+temp003(2,2,1)*(6*r1
     &  0*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2)))
       temp0042(2,2,2)=I22Z*(aux0042(2,2,2)+4*F(5)*temp003(2,2,2)+F(3)*t
     &  emp42(2,2,2)-det3*temp6222(2,k,l)+48*temp00002(2,2)*ZZ(k,2,l,2))
       temp0041(1,1,2)=temp0042(1,1,1)
       temp0041(1,2,1)=temp0042(1,1,1)
       temp0041(1,2,2)=temp0042(2,1,1)
       temp0042(1,1,2)=temp0042(2,1,1)
       temp0042(1,2,1)=temp0042(2,1,1)
       temp0042(1,2,2)=temp0042(2,2,1)
       temp0042(2,1,2)=temp0042(2,2,1)
       temp511(1,1,1)=IX*(aux511(1,1,1)+det3*temp6111(1,1,jj)+10*temp004
     &  1(1,1,1)*Z(jj,1))
       temp521(1,1,1)=IX*(aux521(1,1,1)+det3*temp6211(1,1,jj)+8*temp0042
     &  (1,1,1)*Z(jj,1)+2*temp0041(1,1,1)*Z(jj,2))
       temp522(1,1,1)=IX*(aux522(1,1,1)+det3*temp6221(1,1,jj)+6*temp0042
     &  (2,1,1)*Z(jj,1)+4*temp0042(1,1,1)*Z(jj,2))
       temp522(2,1,1)=IX*(aux522(2,1,1)+det3*temp6222(1,1,jj)+4*temp0042
     &  (2,2,1)*Z(jj,1)+6*temp0042(2,1,1)*Z(jj,2))
       temp522(2,2,1)=IX*(aux522(2,2,1)+det3*temp6222(2,1,jj)+2*temp0042
     &  (2,2,2)*Z(jj,1)+8*temp0042(2,2,1)*Z(jj,2))
       temp522(2,2,2)=IX*(aux522(2,2,2)+det3*temp6222(2,2,jj)+10*temp004
     &  2(2,2,2)*Z(jj,2))
       temp511(1,1,2)=temp521(1,1,1)
       temp511(1,2,1)=temp521(1,1,1)
       temp511(1,2,2)=temp522(1,1,1)
       temp521(1,1,2)=temp522(1,1,1)
       temp521(1,2,1)=temp522(1,1,1)
       temp521(1,2,2)=temp522(2,1,1)
       temp522(1,1,2)=temp522(2,1,1)
       temp522(1,2,1)=temp522(2,1,1)
       temp522(1,2,2)=temp522(2,2,1)
       temp522(2,1,2)=temp522(2,2,1)
c                Step4
       temp00001(1)=I14Z*(aux00001(1)+2*tempC30000*F(4)+F(3)*temp001(1)-
     &  det3*temp003(1,k,l))
       temp00001(2)=I14Z*(aux00001(2)+tempC30000*F(5)+F(3)*temp001(2)-de
     &  t3*temp003(2,k,l))
       temp003(1,1,1)=I18Z*(aux003(1,1,1)+6*F(4)*temp002(1,1)+F(3)*temp3
     &  (1,1,1)-det3*temp511(1,k,l)+24*temp00001(1)*ZZ(k,1,l,1))
       temp003(2,1,1)=I18Z*(aux003(2,1,1)+F(5)*temp002(1,1)+4*F(4)*temp0
     &  02(2,1)+F(3)*temp3(2,1,1)-det3*temp521(1,k,l)+8*(temp00001(2)*ZZ
     &  (k,1,l,1)+temp00001(1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))))
       temp003(2,2,1)=I18Z*(aux003(2,2,1)+2*F(5)*temp002(2,1)+2*F(4)*tem
     &  p002(2,2)+F(3)*temp3(2,2,1)-det3*temp522(1,k,l)+8*(temp00001(2)*
     &  (ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp00001(1)*ZZ(k,2,l,2)))
       temp003(2,2,2)=I18Z*(aux003(2,2,2)+F(3)*temp3(2,2,2)-det3*temp522
     &  (2,k,l)+24*temp00001(2)*ZZ(k,2,l,2)+temp002(2,2)*(6*r10*(ZZ(k,1,
     &  l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2)))
       temp003(1,1,2)=temp003(2,1,1)
       temp003(1,2,1)=temp003(2,1,1)
       temp003(1,2,2)=temp003(2,2,1)
       temp003(2,1,2)=temp003(2,2,1)
       temp41(1,1,1)=IX*(aux41(1,1,1)+det3*temp511(1,1,jj)+8*temp003(1,1
     &  ,1)*Z(jj,1))
        temp42(1,1,1)=IX*(aux42(1,1,1)+det3*temp521(1,1,jj)+6*temp003(2,1
     &  ,1)*Z(jj,1)+2*temp003(1,1,1)*Z(jj,2))
       temp42(2,1,1)=IX*(aux42(2,1,1)+det3*temp522(1,1,jj)+4*(temp003(2,
     &  2,1)*Z(jj,1)+temp003(2,1,1)*Z(jj,2)))
       temp42(2,2,1)=IX*(aux42(2,2,1)+det3*temp522(2,1,jj)+2*temp003(2,2
     &  ,2)*Z(jj,1)+6*temp003(2,2,1)*Z(jj,2))
       temp42(2,2,2)=IX*(aux42(2,2,2)+det3*temp522(2,2,jj)+8*temp003(2,2
     &  ,2)*Z(jj,2))
       temp41(1,1,2)=temp42(1,1,1)
       temp41(1,2,1)=temp42(1,1,1)
       temp41(1,2,2)=temp42(2,1,1)
       temp42(1,1,2)=temp42(2,1,1)
       temp42(1,2,1)=temp42(2,1,1)
       temp42(1,2,2)=temp42(2,2,1)
       temp42(2,1,2)=temp42(2,2,1)
c                Step5
       tempC30000=I10Z*(auxC30000+tempC300*F(3)-det3*temp002(k,l))
       temp002(1,1)=I14Z*(aux002(1,1)+4*F(4)*temp001(1)+F(3)*temp2(1,1)-
     &  det3*temp41(1,k,l)+8*tempC30000*ZZ(k,1,l,1))
       temp002(2,1)=I14Z*(aux002(2,1)+F(5)*temp001(1)+2*F(4)*temp001(2)+
     &  F(3)*temp2(2,1)-det3*temp42(1,k,l)+4*tempC30000*(ZZ(k,1,l,2)+ZZ(
     &  k,2,l,1)))
       temp002(2,2)=I14Z*(aux002(2,2)+2*F(5)*temp001(2)+F(3)*temp2(2,2)-
     &  det3*temp42(2,k,l)+8*tempC30000*ZZ(k,2,l,2))
       temp002(1,2)=temp002(2,1)
       temp3(1,1,1)=IX*(aux3(1,1,1)+det3*temp41(1,1,jj)+6*temp002(1,1)*Z
     &  (jj,1))
       temp3(2,1,1)=IX*(aux3(2,1,1)+det3*temp42(1,1,jj)+4*temp002(2,1)*Z
     &  (jj,1)+2*temp002(1,1)*Z(jj,2))
       temp3(2,2,1)=IX*(aux3(2,2,1)+det3*temp42(2,1,jj)+2*temp002(2,2)*Z
     &  (jj,1)+4*temp002(2,1)*Z(jj,2))
       temp3(2,2,2)=IX*(aux3(2,2,2)+det3*temp42(2,2,jj)+6*temp002(2,2)*Z
     &  (jj,2))
       temp3(1,1,2)=temp3(2,1,1)
       temp3(1,2,1)=temp3(2,1,1)
       temp3(1,2,2)=temp3(2,2,1)
       temp3(2,1,2)=temp3(2,2,1)
c                Step6
       temp001(1)=I10Z*(aux001(1)+2*tempC300*F(4)+F(3)*temp1(1)-det3*tem
     &  p3(1,k,l))
       temp001(2)=I10Z*(aux001(2)+tempC300*F(5)+F(3)*temp1(2)-det3*temp3
     &  (2,k,l))
       temp2(1,1)=IX*(aux2(1,1)+det3*temp3(1,1,jj)+4*temp001(1)*Z(jj,1))
       temp2(2,1)=IX*(aux2(2,1)+det3*temp3(2,1,jj)+2*(temp001(2)*Z(jj,1)
     &  +temp001(1)*Z(jj,2)))
       temp2(2,2)=IX*(aux2(2,2)+det3*temp3(2,2,jj)+4*temp001(2)*Z(jj,2))
       temp2(1,2)=temp2(2,1)
c                Step7
       tempC300=I6Z*(auxC300+tempC30*F(3)-det3*temp2(k,l))
       temp1(1)=IX*(aux1(1)+det3*temp2(1,jj)+2*tempC300*Z(jj,1))
       temp1(2)=IX*(aux1(2)+det3*temp2(2,jj)+2*tempC300*Z(jj,2))
c                Step8
       tempC30=IX*(auxC30+det3*temp1(jj))

         ac=1
         accuracyCR(1,0,ac) = abs(C30     /tempC30      -1d0) 
         accuracyCR(1,1,ac) = abs(Cij(1,1)/temp1(1)      -1d0) 
         accuracyCR(2,1,ac) = abs(Cij(2,1)/temp1(2)      -1d0)
         accuracyCR(1,2,ac) = abs(Cij(1,2)/temp2(1,1)    -1d0)
         accuracyCR(2,2,ac) = abs(Cij(2,2)/temp2(2,2)    -1d0)
         accuracyCR(3,2,ac) = abs(Cij(3,2)/temp2(2,1)    -1d0)
         accuracyCR(4,2,ac) = abs(Cij(4,2)/tempC300      -1d0)
         accuracyCR(1,3,ac) = abs(Cij(1,3)/temp3(1,1,1)  -1d0)
         accuracyCR(2,3,ac) = abs(Cij(2,3)/temp3(2,2,2)  -1d0)
         accuracyCR(3,3,ac) = abs(Cij(3,3)/temp3(2,1,1)  -1d0)
         accuracyCR(4,3,ac) = abs(Cij(4,3)/temp3(2,2,1)  -1d0)
         accuracyCR(5,3,ac) = abs(Cij(5,3)/temp001(1)    -1d0)
         accuracyCR(6,3,ac) = abs(Cij(6,3)/temp001(2)    -1d0)
         accuracyCR(1,4,ac) = abs(Cij(1,4)/temp41(1,1,1) -1d0)
         accuracyCR(2,4,ac) = abs(Cij(2,4)/temp42(2,2,2) -1d0)
         accuracyCR(3,4,ac) = abs(Cij(3,4)/temp42(1,1,1) -1d0)
         accuracyCR(4,4,ac) = abs(Cij(4,4)/temp42(2,1,1) -1d0)
         accuracyCR(5,4,ac) = abs(Cij(5,4)/temp42(2,2,1) -1d0)
         accuracyCR(6,4,ac) = abs(Cij(6,4)/temp002(1,1)  -1d0)
         accuracyCR(7,4,ac) = abs(Cij(7,4)/temp002(2,2)  -1d0)
         accuracyCR(8,4,ac) = abs(Cij(8,4)/temp002(2,1)  -1d0)
         accuracyCR(9,4,ac) = abs(Cij(9,4)/tempC30000    -1d0)


      DO I1=0,4
           accuracyC(i1,ac)=accuracyCR(1,i1,ac)
            DO I2=2,INDEX(I1)  
       if(accuracyCR(i2,i1,ac).gt.accuracyC(i1,ac)) then
          accuracyC(i1,ac)=accuracyCR(i2,i1,ac)
       endif
          enddo
        enddo

c          if(accuracyC(4,ac).lt.1d-16) goto 500

           C30=tempC30
           Cij(1,1)=temp1(1)
           Cij(2,1)=temp1(2)
           Cij(1,2)=temp2(1,1)
           Cij(2,2)=temp2(2,2)
           Cij(3,2)=temp2(2,1)
           Cij(4,2)=tempC300
           Cij(1,3)=temp3(1,1,1)
           Cij(2,3)=temp3(2,2,2)
           Cij(3,3)=temp3(2,1,1)
           Cij(4,3)=temp3(2,2,1)
           Cij(5,3)=temp001(1)
           Cij(6,3)=temp001(2)
           Cij(1,4)=temp41(1,1,1)
           Cij(2,4)=temp42(2,2,2)
           Cij(3,4)=temp42(1,1,1)
           Cij(4,4)=temp42(2,1,1)
           Cij(5,4)=temp42(2,2,1)
           Cij(6,4)=temp002(1,1)
           Cij(7,4)=temp002(2,2)
           Cij(8,4)=temp002(2,1)
           Cij(9,4)=tempC30000

       if(order.eq.7) goto 500
c                Iteration8
c                Step1
       S2000000001(1)=-2*B23(4,6)
       S2000000001(2)=2*B23(4,7)
       S2h000000001(1)=B13(5,8)-B23(5,8)
       S2h000000001(2)=B12(5,8)-B13(5,8)
       S2000000311(1)=-2*B23(3,4)
       S2000000311(2)=2*B23(3,5)
       S2000000312(1)=2*B23(3,5)
       S2000000312(2)=-2*B23(3,6)
       S2000000321(1)=2*B23(3,5)
       S2000000321(2)=-2*B23(3,6)
       S2000000322(1)=-2*B23(3,6)
       S2000000322(2)=2*B23(3,7)
       S2h000000311(1)=B13(4,8)-B23(4,6)
       S2h000000311(2)=B13(4,8)+B23(4,7)
       S2h000000312(1)=B13(4,8)+B23(4,7)
       S2h000000312(2)=B13(4,8)-B23(4,8)
       S2h000000321(1)=B12(4,8)-B13(4,8)
       S2h000000321(2)=-B13(4,8)
       S2h000000322(1)=-B13(4,8)
       S2h000000322(2)=-B13(4,8)
       aux000000001(1)=-(F(1)*S2h00000021(1))-F(2)*S2h00000022(1)+S2h000
     &  000311(k)*Z(1,l)+S2h000000321(k)*Z(2,l)+(-2*Inv000000001+S200000
     &  0001(1)-S2h000000311(1)-S2h000000322(1))*Z(k,l)-2*S2h000000001(1
     &  )*ZZ(k,1,l,1)-2*S2h000000001(2)*ZZ(k,1,l,2)
       aux000000001(2)=-(F(1)*S2h00000021(2))-F(2)*S2h00000022(2)+S2h000
     &  000312(k)*Z(1,l)+S2h000000322(k)*Z(2,l)+(-2*Inv000000002+S200000
     &  0001(2)-S2h000000312(1)-S2h000000322(2))*Z(k,l)-2*S2h000000001(1
     &  )*ZZ(k,2,l,1)-2*S2h000000001(2)*ZZ(k,2,l,2)
       temp000000001(1)=I22Z*(aux000000001(1)+2*tempC300000000*F(4)+F(3)
     &  *temp0000001(1))
       temp000000001(2)=I22Z*(aux000000001(2)+tempC300000000*F(5)+F(3)*t
     &  emp0000001(2))
       S2000051111(1)=-2*B23(2,2)
       S2000051111(2)=2*B23(2,3)
       S2000051211(1)=2*B23(2,3)
       S2000051211(2)=-2*B23(2,4)
       S2000051221(1)=-2*B23(2,4)
       S2000051221(2)=2*B23(2,5)
       S2000051222(1)=2*B23(2,5)
       S2000051222(2)=-2*B23(2,6)
       S2000052111(1)=2*B23(2,3)
       S2000052111(2)=-2*B23(2,4)
       S2000052211(1)=-2*B23(2,4)
       S2000052211(2)=2*B23(2,5)
       S2000052221(1)=2*B23(2,5)
       S2000052221(2)=-2*B23(2,6)
       S2000052222(1)=-2*B23(2,6)
       S2000052222(2)=2*B23(2,7)
       S2h000051111(1)=B13(3,8)-B23(3,4)
       S2h000051111(2)=B13(3,8)+B23(3,5)
       S2h000051211(1)=B13(3,8)+B23(3,5)
       S2h000051211(2)=B13(3,8)-B23(3,6)
       S2h000051221(1)=B13(3,8)-B23(3,6)
       S2h000051221(2)=B13(3,8)+B23(3,7)
       S2h000051222(1)=B13(3,8)+B23(3,7)
       S2h000051222(2)=B13(3,8)-B23(3,8)
       S2h000052111(1)=B12(3,8)-B13(3,8)
       S2h000052111(2)=-B13(3,8)
       S2h000052211(1)=-B13(3,8)
       S2h000052211(2)=-B13(3,8)
       S2h000052221(1)=-B13(3,8)
       S2h000052221(2)=-B13(3,8)
       S2h000052222(1)=-B13(3,8)
       S2h000052222(2)=-B13(3,8)
       aux0000003(1,1,1)=-(F(1)*S2h00004111(1))-F(2)*S2h00004211(1)+S2h0
     &  00051111(k)*Z(1,l)+S2h000052111(k)*Z(2,l)+(-2*Inv000000111+S2000
     &  000311(1)-S2h000051111(1)-S2h000052211(1))*Z(k,l)-6*(S2h00000031
     &  1(1)*ZZ(k,1,l,1)+S2h000000321(1)*ZZ(k,1,l,2))
       aux0000003(2,1,1)=-(F(1)*S2h00004121(1))-F(2)*S2h00004221(1)+S2h0
     &  00051211(k)*Z(1,l)+S2h000052211(k)*Z(2,l)+(-2*Inv000000211+S2000
     &  000321(1)-S2h000051211(1)-S2h000052221(1))*Z(k,l)-4*(S2h00000031
     &  2(1)*ZZ(k,1,l,1)+S2h000000322(1)*ZZ(k,1,l,2))-2*S2h000000311(1)*
     &  ZZ(k,2,l,1)-2*S2h000000321(1)*ZZ(k,2,l,2)
       aux0000003(2,2,1)=-(F(1)*S2h00004122(1))-F(2)*S2h00004222(1)+S2h0
     &  00051221(k)*Z(1,l)+S2h000052221(k)*Z(2,l)+(-2*Inv000000221+S2000
     &  000322(1)-S2h000051221(1)-S2h000052222(1))*Z(k,l)-2*(S2h00000031
     &  2(2)*ZZ(k,1,l,1)+S2h000000322(2)*ZZ(k,1,l,2))-4*S2h000000312(1)*
     &  ZZ(k,2,l,1)-4*S2h000000322(1)*ZZ(k,2,l,2)
       aux0000003(2,2,2)=-(F(1)*S2h00004122(2))-F(2)*S2h00004222(2)+S2h0
     &  00051222(k)*Z(1,l)+S2h000052222(k)*Z(2,l)+(-2*Inv000000222+S2000
     &  000322(2)-S2h000051222(1)-S2h000052222(2))*Z(k,l)-6*(S2h00000031
     &  2(2)*ZZ(k,2,l,1)+S2h000000322(2)*ZZ(k,2,l,2))
       temp0000003(1,1,1)=I26Z*(aux0000003(1,1,1)+6*F(4)*temp0000002(1,1
     &  )+F(3)*temp00003(1,1,1)+24*temp000000001(1)*ZZ(k,1,l,1))
       temp0000003(2,1,1)=I26Z*(aux0000003(2,1,1)+F(5)*temp0000002(1,1)+
     &  4*F(4)*temp0000002(2,1)+F(3)*temp00003(2,1,1)+8*(temp000000001(2
     &  )*ZZ(k,1,l,1)+temp000000001(1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))))
       temp0000003(2,2,1)=I26Z*(aux0000003(2,2,1)+2*F(5)*temp0000002(2,1
     &  )+2*F(4)*temp0000002(2,2)+F(3)*temp00003(2,2,1)+8*(temp000000001
     &  (2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp000000001(1)*ZZ(k,2,l,2)))
       temp0000003(2,2,2)=I26Z*(aux0000003(2,2,2)+F(3)*temp00003(2,2,2)+
     &  24*temp000000001(2)*ZZ(k,2,l,2)+temp0000002(2,2)*(6*r10*(ZZ(k,1,
     &  l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2)))
       temp0000003(1,1,2)=temp0000003(2,1,1)
       temp0000003(1,2,1)=temp0000003(2,1,1)
       temp0000003(1,2,2)=temp0000003(2,2,1)
       temp0000003(2,1,2)=temp0000003(2,2,1)
       S2007111111(1)=-2*B023
       S2007111111(2)=2*B23(1,1)
       S2007121111(1)=2*B23(1,1)
       S2007121111(2)=-2*B23(1,2)
       S2007122111(1)=-2*B23(1,2)
       S2007122111(2)=2*B23(1,3)
       S2007122211(1)=2*B23(1,3)
       S2007122211(2)=-2*B23(1,4)
       S2007122221(1)=-2*B23(1,4)
       S2007122221(2)=2*B23(1,5)
       S2007122222(1)=2*B23(1,5)
       S2007122222(2)=-2*B23(1,6)
       S2007211111(1)=2*B23(1,1)
       S2007211111(2)=-2*B23(1,2)
       S2007221111(1)=-2*B23(1,2)
       S2007221111(2)=2*B23(1,3)
       S2007222111(1)=2*B23(1,3)
       S2007222111(2)=-2*B23(1,4)
       S2007222211(1)=-2*B23(1,4)
       S2007222211(2)=2*B23(1,5)
       S2007222221(1)=2*B23(1,5)
       S2007222221(2)=-2*B23(1,6)
       S2007222222(1)=-2*B23(1,6)
       S2007222222(2)=2*B23(1,7)
       S2h007111111(1)=B13(2,8)-B23(2,2)
       S2h007111111(2)=B13(2,8)+B23(2,3)
       S2h007121111(1)=B13(2,8)+B23(2,3)
       S2h007121111(2)=B13(2,8)-B23(2,4)
       S2h007122111(1)=B13(2,8)-B23(2,4)
       S2h007122111(2)=B13(2,8)+B23(2,5)
       S2h007122211(1)=B13(2,8)+B23(2,5)
       S2h007122211(2)=B13(2,8)-B23(2,6)
       S2h007122221(1)=B13(2,8)-B23(2,6)
       S2h007122221(2)=B13(2,8)+B23(2,7)
       S2h007122222(1)=B13(2,8)+B23(2,7)
       S2h007122222(2)=B13(2,8)-B23(2,8)
       S2h007211111(1)=B12(2,8)-B13(2,8)
       S2h007211111(2)=-B13(2,8)
       S2h007221111(1)=-B13(2,8)
       S2h007221111(2)=-B13(2,8)
       S2h007222111(1)=-B13(2,8)
       S2h007222111(2)=-B13(2,8)
       S2h007222211(1)=-B13(2,8)
       S2h007222211(2)=-B13(2,8)
       S2h007222221(1)=-B13(2,8)
       S2h007222221(2)=-B13(2,8)
       S2h007222222(1)=-B13(2,8)
       S2h007222222(2)=-B13(2,8)
       aux0000511(1,1,1)=-(F(1)*S2h00611111(1))-F(2)*S2h00621111(1)+S2h0
     &  07111111(k)*Z(1,l)+S2h007211111(k)*Z(2,l)+(-2*Inv000011111+S2000
     &  051111(1)-S2h007111111(1)-S2h007221111(1))*Z(k,l)-10*(S2h0000511
     &  11(1)*ZZ(k,1,l,1)+S2h000052111(1)*ZZ(k,1,l,2))
       aux0000521(1,1,1)=-(F(1)*S2h00612111(1))-F(2)*S2h00622111(1)+S2h0
     &  07121111(k)*Z(1,l)+S2h007221111(k)*Z(2,l)+(-2*Inv000021111+S2000
     &  052111(1)-S2h007121111(1)-S2h007222111(1))*Z(k,l)-8*(S2h00005121
     &  1(1)*ZZ(k,1,l,1)+S2h000052211(1)*ZZ(k,1,l,2))-2*S2h000051111(1)*
     &  ZZ(k,2,l,1)-2*S2h000052111(1)*ZZ(k,2,l,2)
       aux0000522(1,1,1)=-(F(1)*S2h00612211(1))-F(2)*S2h00622211(1)+S2h0
     &  07122111(k)*Z(1,l)+S2h007222111(k)*Z(2,l)+(-2*Inv000022111+S2000
     &  052211(1)-S2h007122111(1)-S2h007222211(1))*Z(k,l)-6*(S2h00005122
     &  1(1)*ZZ(k,1,l,1)+S2h000052221(1)*ZZ(k,1,l,2))-4*S2h000051211(1)*
     &  ZZ(k,2,l,1)-4*S2h000052211(1)*ZZ(k,2,l,2)
       aux0000522(2,1,1)=-(F(1)*S2h00612221(1))-F(2)*S2h00622221(1)+S2h0
     &  07122211(k)*Z(1,l)+S2h007222211(k)*Z(2,l)+(-2*Inv000022211+S2000
     &  052221(1)-S2h007122211(1)-S2h007222221(1))*Z(k,l)-4*(S2h00005122
     &  2(1)*ZZ(k,1,l,1)+S2h000052222(1)*ZZ(k,1,l,2))-6*S2h000051221(1)*
     &  ZZ(k,2,l,1)-6*S2h000052221(1)*ZZ(k,2,l,2)
       aux0000522(2,2,1)=-(F(1)*S2h00612222(1))-F(2)*S2h00622222(1)+S2h0
     &  07122221(k)*Z(1,l)+S2h007222221(k)*Z(2,l)+(-2*Inv000022221+S2000
     &  052222(1)-S2h007122221(1)-S2h007222222(1))*Z(k,l)-2*(S2h00005122
     &  2(2)*ZZ(k,1,l,1)+S2h000052222(2)*ZZ(k,1,l,2))-8*S2h000051222(1)*
     &  ZZ(k,2,l,1)-8*S2h000052222(1)*ZZ(k,2,l,2)
       aux0000522(2,2,2)=-(F(1)*S2h00612222(2))-F(2)*S2h00622222(2)+S2h0
     &  07122222(k)*Z(1,l)+S2h007222222(k)*Z(2,l)+(-2*Inv000022222+S2000
     &  052222(2)-S2h007122222(1)-S2h007222222(2))*Z(k,l)-10*(S2h0000512
     &  22(2)*ZZ(k,2,l,1)+S2h000052222(2)*ZZ(k,2,l,2))
       temp0000511(1,1,1)=I30Z*(aux0000511(1,1,1)+10*F(4)*temp000041(1,1
     &  ,1)+F(3)*temp00511(1,1,1)+80*temp0000003(1,1,1)*ZZ(k,1,l,1))
       temp0000521(1,1,1)=I30Z*(aux0000521(1,1,1)+F(5)*temp000041(1,1,1)
     &  +8*F(4)*temp000042(1,1,1)+F(3)*temp00521(1,1,1)+48*temp0000003(2
     &  ,1,1)*ZZ(k,1,l,1)+16*temp0000003(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)
     &  ))
       temp0000522(1,1,1)=I30Z*(aux0000522(1,1,1)+2*F(5)*temp000042(1,1,
     &  1)+6*F(4)*temp000042(2,1,1)+F(3)*temp00522(1,1,1)+24*(temp000000
     &  3(2,2,1)*ZZ(k,1,l,1)+temp0000003(2,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)
     &  ))+8*temp0000003(1,1,1)*ZZ(k,2,l,2))
       temp0000522(2,1,1)=I30Z*(aux0000522(2,1,1)+4*F(4)*temp000042(2,2,
     &  1)+F(3)*temp00522(2,1,1)+8*temp0000003(2,2,2)*ZZ(k,1,l,1)+temp00
     &  0042(2,1,1)*(6*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2))
     &  +24*(temp0000003(2,2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp0000003(2,
     &  1,1)*ZZ(k,2,l,2)))
       temp0000522(2,2,1)=I30Z*(aux0000522(2,2,1)+4*F(5)*temp000042(2,2,
     &  1)+2*F(4)*temp000042(2,2,2)+F(3)*temp00522(2,2,1)+16*temp0000003
     &  (2,2,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+48*temp0000003(2,2,1)*ZZ(k,2,l
     &  ,2))
       temp0000522(2,2,2)=I30Z*(aux0000522(2,2,2)+F(3)*temp00522(2,2,2)+
     &  80*temp0000003(2,2,2)*ZZ(k,2,l,2)+temp000042(2,2,2)*(10*r10*(ZZ(
     &  k,1,l,2)+ZZ(k,2,l,1))+20*r21*ZZ(k,2,l,2)))
       temp0000511(1,1,2)=temp0000521(1,1,1)
       temp0000511(1,2,1)=temp0000521(1,1,1)
       temp0000511(1,2,2)=temp0000522(1,1,1)
       temp0000521(1,1,2)=temp0000522(1,1,1)
       temp0000521(1,2,1)=temp0000522(1,1,1)
       temp0000521(1,2,2)=temp0000522(2,1,1)
       temp0000522(1,1,2)=temp0000522(2,1,1)
       temp0000522(1,2,1)=temp0000522(2,1,1)
       temp0000522(1,2,2)=temp0000522(2,2,1)
       temp0000522(2,1,2)=temp0000522(2,2,1)
       S2911111111(1)=-B023+B13(1,8)
       S2911111111(2)=B13(1,8)+B23(1,1)
       S2912111111(1)=B13(1,8)+B23(1,1)
       S2912111111(2)=B13(1,8)-B23(1,2)
       S2912211111(1)=B13(1,8)-B23(1,2)
       S2912211111(2)=B13(1,8)+B23(1,3)
       S2912221111(1)=B13(1,8)+B23(1,3)
       S2912221111(2)=B13(1,8)-B23(1,4)
       S2912222111(1)=B13(1,8)-B23(1,4)
       S2912222111(2)=B13(1,8)+B23(1,5)
       S2912222211(1)=B13(1,8)+B23(1,5)
       S2912222211(2)=B13(1,8)-B23(1,6)
       S2912222221(1)=B13(1,8)-B23(1,6)
       S2912222221(2)=B13(1,8)+B23(1,7)
       S2912222222(1)=B13(1,8)+B23(1,7)
       S2912222222(2)=B13(1,8)-B23(1,8)
       S2921111111(1)=B12(1,8)-B13(1,8)
       S2921111111(2)=-B13(1,8)
       S2922111111(1)=-B13(1,8)
       S2922111111(2)=-B13(1,8)
       S2922211111(1)=-B13(1,8)
       S2922211111(2)=-B13(1,8)
       S2922221111(1)=-B13(1,8)
       S2922221111(2)=-B13(1,8)
       S2922222111(1)=-B13(1,8)
       S2922222111(2)=-B13(1,8)
       S2922222211(1)=-B13(1,8)
       S2922222211(2)=-B13(1,8)
       S2922222221(1)=-B13(1,8)
       S2922222221(2)=-B13(1,8)
       S2922222222(1)=-B13(1,8)
       S2922222222(2)=-B13(1,8)
       aux0071111(1,1,1)=-(F(1)*S281111111(1))-F(2)*S282111111(1)+S29111
     &  11111(k)*Z(1,l)+S2921111111(k)*Z(2,l)+(-2*Inv001111111+S20071111
     &  11(1)-S2911111111(1)-S2922111111(1))*Z(k,l)-14*(S2h007111111(1)*
     &  ZZ(k,1,l,1)+S2h007211111(1)*ZZ(k,1,l,2))
       aux0072111(1,1,1)=-(F(1)*S281211111(1))-F(2)*S282211111(1)+S29121
     &  11111(k)*Z(1,l)+S2922111111(k)*Z(2,l)+(-2*Inv002111111+S20072111
     &  11(1)-S2912111111(1)-S2922211111(1))*Z(k,l)-12*(S2h007121111(1)*
     &  ZZ(k,1,l,1)+S2h007221111(1)*ZZ(k,1,l,2))-2*S2h007111111(1)*ZZ(k,
     &  2,l,1)-2*S2h007211111(1)*ZZ(k,2,l,2)
       aux0072211(1,1,1)=-(F(1)*S281221111(1))-F(2)*S282221111(1)+S29122
     &  11111(k)*Z(1,l)+S2922211111(k)*Z(2,l)+(-2*Inv002211111+S20072211
     &  11(1)-S2912211111(1)-S2922221111(1))*Z(k,l)-10*(S2h007122111(1)*
     &  ZZ(k,1,l,1)+S2h007222111(1)*ZZ(k,1,l,2))-4*S2h007121111(1)*ZZ(k,
     &  2,l,1)-4*S2h007221111(1)*ZZ(k,2,l,2)
       aux0072221(1,1,1)=-(F(1)*S281222111(1))-F(2)*S282222111(1)+S29122
     &  21111(k)*Z(1,l)+S2922221111(k)*Z(2,l)+(-2*Inv002221111+S20072221
     &  11(1)-S2912221111(1)-S2922222111(1))*Z(k,l)-8*(S2h007122211(1)*Z
     &  Z(k,1,l,1)+S2h007222211(1)*ZZ(k,1,l,2))-6*S2h007122111(1)*ZZ(k,2
     &  ,l,1)-6*S2h007222111(1)*ZZ(k,2,l,2)
       aux0072222(1,1,1)=-(F(1)*S281222211(1))-F(2)*S282222211(1)+S29122
     &  22111(k)*Z(1,l)+S2922222111(k)*Z(2,l)+(-2*Inv002222111+S20072222
     &  11(1)-S2912222111(1)-S2922222211(1))*Z(k,l)-6*(S2h007122221(1)*Z
     &  Z(k,1,l,1)+S2h007222221(1)*ZZ(k,1,l,2))-8*S2h007122211(1)*ZZ(k,2
     &  ,l,1)-8*S2h007222211(1)*ZZ(k,2,l,2)
       aux0072222(2,1,1)=-(F(1)*S281222221(1))-F(2)*S282222221(1)+S29122
     &  22211(k)*Z(1,l)+S2922222211(k)*Z(2,l)+(-2*Inv002222211+S20072222
     &  21(1)-S2912222211(1)-S2922222221(1))*Z(k,l)-4*(S2h007122222(1)*Z
     &  Z(k,1,l,1)+S2h007222222(1)*ZZ(k,1,l,2))-10*S2h007122221(1)*ZZ(k,
     &  2,l,1)-10*S2h007222221(1)*ZZ(k,2,l,2)
       aux0072222(2,2,1)=-(F(1)*S281222222(1))-F(2)*S282222222(1)+S29122
     &  22221(k)*Z(1,l)+S2922222221(k)*Z(2,l)+(-2*Inv002222221+S20072222
     &  22(1)-S2912222221(1)-S2922222222(1))*Z(k,l)-2*(S2h007122222(2)*Z
     &  Z(k,1,l,1)+S2h007222222(2)*ZZ(k,1,l,2))-12*S2h007122222(1)*ZZ(k,
     &  2,l,1)-12*S2h007222222(1)*ZZ(k,2,l,2)
       aux0072222(2,2,2)=-(F(1)*S281222222(2))-F(2)*S282222222(2)+S29122
     &  22222(k)*Z(1,l)+S2922222222(k)*Z(2,l)+(-2*Inv002222222+S20072222
     &  22(2)-S2912222222(1)-S2922222222(2))*Z(k,l)-14*(S2h007122222(2)*
     &  ZZ(k,2,l,1)+S2h007222222(2)*ZZ(k,2,l,2))
       temp0071111(1,1,1)=I34Z*(aux0071111(1,1,1)+14*F(4)*temp006111(1,1
     &  ,1)+F(3)*temp71111(1,1,1)+168*temp0000511(1,1,1)*ZZ(k,1,l,1))
       temp0072111(1,1,1)=I34Z*(aux0072111(1,1,1)+F(5)*temp006111(1,1,1)
     &  +12*F(4)*temp006211(1,1,1)+F(3)*temp72111(1,1,1)+120*temp0000521
     &  (1,1,1)*ZZ(k,1,l,1)+24*temp0000511(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,
     &  1)))
       temp0072211(1,1,1)=I34Z*(aux0072211(1,1,1)+2*F(5)*temp006211(1,1,
     &  1)+10*F(4)*temp006221(1,1,1)+F(3)*temp72211(1,1,1)+80*temp000052
     &  2(1,1,1)*ZZ(k,1,l,1)+40*temp0000521(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l
     &  ,1))+8*temp0000511(1,1,1)*ZZ(k,2,l,2))
       temp0072221(1,1,1)=I34Z*(aux0072221(1,1,1)+8*F(4)*temp006222(1,1,
     &  1)+F(3)*temp72221(1,1,1)+48*(temp0000522(2,1,1)*ZZ(k,1,l,1)+temp
     &  0000522(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))+24*temp0000521(1,1,1)*
     &  ZZ(k,2,l,2)+temp006221(1,1,1)*(6*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+1
     &  2*r21*ZZ(k,2,l,2)))
       temp0072222(1,1,1)=I34Z*(aux0072222(1,1,1)+4*F(5)*temp006222(1,1,
     &  1)+6*F(4)*temp006222(2,1,1)+F(3)*temp72222(1,1,1)+24*temp0000522
     &  (2,2,1)*ZZ(k,1,l,1)+48*(temp0000522(2,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l
     &  ,1))+temp0000522(1,1,1)*ZZ(k,2,l,2)))
       temp0072222(2,1,1)=I34Z*(aux0072222(2,1,1)+4*F(4)*temp006222(2,2,
     &  1)+F(3)*temp72222(2,1,1)+8*temp0000522(2,2,2)*ZZ(k,1,l,1)+40*tem
     &  p0000522(2,2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+80*temp0000522(2,1,1)*
     &  ZZ(k,2,l,2)+temp006222(2,1,1)*(10*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+
     &  20*r21*ZZ(k,2,l,2)))
       temp0072222(2,2,1)=I34Z*(aux0072222(2,2,1)+2*F(4)*temp006222(2,2,
     &  2)+F(3)*temp72222(2,2,1)+24*temp0000522(2,2,2)*(ZZ(k,1,l,2)+ZZ(k
     &  ,2,l,1))+120*temp0000522(2,2,1)*ZZ(k,2,l,2)+temp006222(2,2,1)*(1
     &  2*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+24*r21*ZZ(k,2,l,2)))
       temp0072222(2,2,2)=I34Z*(aux0072222(2,2,2)+F(3)*temp72222(2,2,2)+
     &  168*temp0000522(2,2,2)*ZZ(k,2,l,2)+temp006222(2,2,2)*(14*r10*(ZZ
     &  (k,1,l,2)+ZZ(k,2,l,1))+28*r21*ZZ(k,2,l,2)))
       temp0071111(1,1,2)=temp0072111(1,1,1)
       temp0071111(1,2,1)=temp0072111(1,1,1)
       temp0071111(1,2,2)=temp0072211(1,1,1)
       temp0072111(1,1,2)=temp0072211(1,1,1)
       temp0072111(1,2,1)=temp0072211(1,1,1)
       temp0072111(1,2,2)=temp0072221(1,1,1)
       temp0072211(1,1,2)=temp0072221(1,1,1)
       temp0072211(1,2,1)=temp0072221(1,1,1)
       temp0072211(1,2,2)=temp0072222(1,1,1)
       temp0072221(1,1,2)=temp0072222(1,1,1)
       temp0072221(1,2,1)=temp0072222(1,1,1)
       temp0072221(1,2,2)=temp0072222(2,1,1)
       temp0072222(1,1,2)=temp0072222(2,1,1)
       temp0072222(1,2,1)=temp0072222(2,1,1)
       temp0072222(1,2,2)=temp0072222(2,2,1)
       temp0072222(2,1,2)=temp0072222(2,2,1)
       aux811111(1,1,1)=-(S2911111111(1)*Z(jj,1))-S2921111111(1)*Z(jj,2)
       aux821111(1,1,1)=-(S2912111111(1)*Z(jj,1))-S2922111111(1)*Z(jj,2)
       aux822111(1,1,1)=-(S2912211111(1)*Z(jj,1))-S2922211111(1)*Z(jj,2)
       aux822211(1,1,1)=-(S2912221111(1)*Z(jj,1))-S2922221111(1)*Z(jj,2)
       aux822221(1,1,1)=-(S2912222111(1)*Z(jj,1))-S2922222111(1)*Z(jj,2)
       aux822222(1,1,1)=-(S2912222211(1)*Z(jj,1))-S2922222211(1)*Z(jj,2)
       aux822222(2,1,1)=-(S2912222221(1)*Z(jj,1))-S2922222221(1)*Z(jj,2)
       aux822222(2,2,1)=-(S2912222222(1)*Z(jj,1))-S2922222222(1)*Z(jj,2)
       aux822222(2,2,2)=-(S2912222222(2)*Z(jj,1))-S2922222222(2)*Z(jj,2)
       temp811111(1,1,1)=IX*(aux811111(1,1,1)+16*temp0071111(1,1,1)*Z(jj
     &  ,1))
       temp821111(1,1,1)=IX*(aux821111(1,1,1)+14*temp0072111(1,1,1)*Z(jj
     &  ,1)+2*temp0071111(1,1,1)*Z(jj,2))
       temp822111(1,1,1)=IX*(aux822111(1,1,1)+12*temp0072211(1,1,1)*Z(jj
     &  ,1)+4*temp0072111(1,1,1)*Z(jj,2))
       temp822211(1,1,1)=IX*(aux822211(1,1,1)+10*temp0072221(1,1,1)*Z(jj
     &  ,1)+6*temp0072211(1,1,1)*Z(jj,2))
       temp822221(1,1,1)=IX*(aux822221(1,1,1)+8*(temp0072222(1,1,1)*Z(jj
     &  ,1)+temp0072221(1,1,1)*Z(jj,2)))
       temp822222(1,1,1)=IX*(aux822222(1,1,1)+6*temp0072222(2,1,1)*Z(jj,
     &  1)+10*temp0072222(1,1,1)*Z(jj,2))
       temp822222(2,1,1)=IX*(aux822222(2,1,1)+4*temp0072222(2,2,1)*Z(jj,
     &  1)+12*temp0072222(2,1,1)*Z(jj,2))
       temp822222(2,2,1)=IX*(aux822222(2,2,1)+2*temp0072222(2,2,2)*Z(jj,
     &  1)+14*temp0072222(2,2,1)*Z(jj,2))
       temp822222(2,2,2)=IX*(aux822222(2,2,2)+16*temp0072222(2,2,2)*Z(jj
     &  ,2))
       temp811111(1,1,2)=temp821111(1,1,1)
       temp811111(1,2,1)=temp821111(1,1,1)
       temp811111(1,2,2)=temp822111(1,1,1)
       temp821111(1,1,2)=temp822111(1,1,1)
       temp821111(1,2,1)=temp822111(1,1,1)
       temp821111(1,2,2)=temp822211(1,1,1)
       temp822111(1,1,2)=temp822211(1,1,1)
       temp822111(1,2,1)=temp822211(1,1,1)
       temp822111(1,2,2)=temp822221(1,1,1)
       temp822211(1,1,2)=temp822221(1,1,1)
       temp822211(1,2,1)=temp822221(1,1,1)
       temp822211(1,2,2)=temp822222(1,1,1)
       temp822221(1,1,2)=temp822222(1,1,1)
       temp822221(1,2,1)=temp822222(1,1,1)
       temp822221(1,2,2)=temp822222(2,1,1)
       temp822222(1,1,2)=temp822222(2,1,1)
       temp822222(1,2,1)=temp822222(2,1,1)
       temp822222(1,2,2)=temp822222(2,2,1)
       temp822222(2,1,2)=temp822222(2,2,1)
c                Step2
       tempC300000000=I18Z*(auxC300000000+tempC3000000*F(3)-det3*temp000
     &  0002(k,l))
       temp0000002(1,1)=I22Z*(aux0000002(1,1)+4*F(4)*temp0000001(1)+F(3)
     &  *temp00002(1,1)-det3*temp000041(1,k,l)+8*tempC300000000*ZZ(k,1,l
     &  ,1))
       temp0000002(2,1)=I22Z*(aux0000002(2,1)+F(5)*temp0000001(1)+2*F(4)
     &  *temp0000001(2)+F(3)*temp00002(2,1)-det3*temp000042(1,k,l)+4*tem
     &  pC300000000*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp0000002(2,2)=I22Z*(aux0000002(2,2)+2*F(5)*temp0000001(2)+F(3)
     &  *temp00002(2,2)-det3*temp000042(2,k,l)+8*tempC300000000*ZZ(k,2,l
     &  ,2))
       temp0000002(1,2)=temp0000002(2,1)
       temp000041(1,1,1)=I26Z*(aux000041(1,1,1)+8*F(4)*temp00003(1,1,1)+
     &  F(3)*temp0041(1,1,1)-det3*temp006111(1,k,l)+48*temp0000002(1,1)*
     &  ZZ(k,1,l,1))
       temp000042(1,1,1)=I26Z*(aux000042(1,1,1)+F(5)*temp00003(1,1,1)+6*
     &  F(4)*temp00003(2,1,1)+F(3)*temp0042(1,1,1)-det3*temp006211(1,k,l
     &  )+24*temp0000002(2,1)*ZZ(k,1,l,1)+12*temp0000002(1,1)*(ZZ(k,1,l,
     &  2)+ZZ(k,2,l,1)))
       temp000042(2,1,1)=I26Z*(aux000042(2,1,1)+2*F(5)*temp00003(2,1,1)+
     &  4*F(4)*temp00003(2,2,1)+F(3)*temp0042(2,1,1)-det3*temp006221(1,k
     &  ,l)+16*temp0000002(2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+8*(temp0000002
     &  (2,2)*ZZ(k,1,l,1)+temp0000002(1,1)*ZZ(k,2,l,2)))
       temp000042(2,2,1)=I26Z*(aux000042(2,2,1)+2*F(4)*temp00003(2,2,2)+
     &  F(3)*temp0042(2,2,1)-det3*temp006222(1,k,l)+12*temp0000002(2,2)*
     &  (ZZ(k,1,l,2)+ZZ(k,2,l,1))+24*temp0000002(2,1)*ZZ(k,2,l,2)+temp00
     &  003(2,2,1)*(6*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2)))
       temp000042(2,2,2)=I26Z*(aux000042(2,2,2)+4*F(5)*temp00003(2,2,2)+
     &  F(3)*temp0042(2,2,2)-det3*temp006222(2,k,l)+48*temp0000002(2,2)*
     &  ZZ(k,2,l,2))
       temp000041(1,1,2)=temp000042(1,1,1)
       temp000041(1,2,1)=temp000042(1,1,1)
       temp000041(1,2,2)=temp000042(2,1,1)
       temp000042(1,1,2)=temp000042(2,1,1)
       temp000042(1,2,1)=temp000042(2,1,1)
       temp000042(1,2,2)=temp000042(2,2,1)
       temp000042(2,1,2)=temp000042(2,2,1)
       temp006111(1,1,1)=I30Z*(aux006111(1,1,1)+12*F(4)*temp00511(1,1,1)
     &  +F(3)*temp6111(1,1,1)-det3*temp811111(1,k,l)+120*temp000041(1,1,
     &  1)*ZZ(k,1,l,1))
       temp006211(1,1,1)=I30Z*(aux006211(1,1,1)+F(5)*temp00511(1,1,1)+10
     &  *F(4)*temp00521(1,1,1)+F(3)*temp6211(1,1,1)-det3*temp821111(1,k,
     &  l)+80*temp000042(1,1,1)*ZZ(k,1,l,1)+20*temp000041(1,1,1)*(ZZ(k,1
     &  ,l,2)+ZZ(k,2,l,1)))
       temp006221(1,1,1)=I30Z*(aux006221(1,1,1)+2*F(5)*temp00521(1,1,1)+
     &  F(3)*temp6221(1,1,1)-det3*temp822111(1,k,l)+48*temp000042(2,1,1)
     &  *ZZ(k,1,l,1)+32*temp000042(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+8*(F
     &  (4)*temp00522(1,1,1)+temp000041(1,1,1)*ZZ(k,2,l,2)))
       temp006222(1,1,1)=I30Z*(aux006222(1,1,1)+6*F(4)*temp00522(2,1,1)+
     &  F(3)*temp6222(1,1,1)-det3*temp822211(1,k,l)+36*temp000042(2,1,1)
     &  *(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp00522(1,1,1)*(6*r10*(ZZ(k,1,l,2)+
     &  ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2))+24*(temp000042(2,2,1)*ZZ(k,1,l,
     &  1)+temp000042(1,1,1)*ZZ(k,2,l,2)))
       temp006222(2,1,1)=I30Z*(aux006222(2,1,1)+4*F(5)*temp00522(2,1,1)+
     &  4*F(4)*temp00522(2,2,1)+F(3)*temp6222(2,1,1)-det3*temp822221(1,k
     &  ,l)+8*temp000042(2,2,2)*ZZ(k,1,l,1)+32*temp000042(2,2,1)*(ZZ(k,1
     &  ,l,2)+ZZ(k,2,l,1))+48*temp000042(2,1,1)*ZZ(k,2,l,2))
       temp006222(2,2,1)=I30Z*(aux006222(2,2,1)+2*F(4)*temp00522(2,2,2)+
     &  F(3)*temp6222(2,2,1)-det3*temp822222(1,k,l)+20*temp000042(2,2,2)
     &  *(ZZ(k,1,l,2)+ZZ(k,2,l,1))+80*temp000042(2,2,1)*ZZ(k,2,l,2)+temp
     &  00522(2,2,1)*(10*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+20*r21*ZZ(k,2,l,2
     &  )))
       temp006222(2,2,2)=I30Z*(aux006222(2,2,2)+F(3)*temp6222(2,2,2)-det
     &  3*temp822222(2,k,l)+120*temp000042(2,2,2)*ZZ(k,2,l,2)+temp00522(
     &  2,2,2)*(12*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+24*r21*ZZ(k,2,l,2)))
       temp006111(1,1,2)=temp006211(1,1,1)
       temp006111(1,2,1)=temp006211(1,1,1)
       temp006111(1,2,2)=temp006221(1,1,1)
       temp006211(1,1,2)=temp006221(1,1,1)
       temp006211(1,2,1)=temp006221(1,1,1)
       temp006211(1,2,2)=temp006222(1,1,1)
       temp006221(1,1,2)=temp006222(1,1,1)
       temp006221(1,2,1)=temp006222(1,1,1)
       temp006221(1,2,2)=temp006222(2,1,1)
       temp006222(1,1,2)=temp006222(2,1,1)
       temp006222(1,2,1)=temp006222(2,1,1)
       temp006222(1,2,2)=temp006222(2,2,1)
       temp006222(2,1,2)=temp006222(2,2,1)
       temp71111(1,1,1)=IX*(aux71111(1,1,1)+det3*temp811111(1,1,jj)+14*t
     &  emp006111(1,1,1)*Z(jj,1))
       temp72111(1,1,1)=IX*(aux72111(1,1,1)+det3*temp821111(1,1,jj)+12*t
     &  emp006211(1,1,1)*Z(jj,1)+2*temp006111(1,1,1)*Z(jj,2))
       temp72211(1,1,1)=IX*(aux72211(1,1,1)+det3*temp822111(1,1,jj)+10*t
     &  emp006221(1,1,1)*Z(jj,1)+4*temp006211(1,1,1)*Z(jj,2))
       temp72221(1,1,1)=IX*(aux72221(1,1,1)+det3*temp822211(1,1,jj)+8*te
     &  mp006222(1,1,1)*Z(jj,1)+6*temp006221(1,1,1)*Z(jj,2))
       temp72222(1,1,1)=IX*(aux72222(1,1,1)+det3*temp822221(1,1,jj)+6*te
     &  mp006222(2,1,1)*Z(jj,1)+8*temp006222(1,1,1)*Z(jj,2))
       temp72222(2,1,1)=IX*(aux72222(2,1,1)+det3*temp822222(1,1,jj)+4*te
     &  mp006222(2,2,1)*Z(jj,1)+10*temp006222(2,1,1)*Z(jj,2))
       temp72222(2,2,1)=IX*(aux72222(2,2,1)+det3*temp822222(2,1,jj)+2*te
     &  mp006222(2,2,2)*Z(jj,1)+12*temp006222(2,2,1)*Z(jj,2))
       temp72222(2,2,2)=IX*(aux72222(2,2,2)+det3*temp822222(2,2,jj)+14*t
     &  emp006222(2,2,2)*Z(jj,2))
       temp71111(1,1,2)=temp72111(1,1,1)
       temp71111(1,2,1)=temp72111(1,1,1)
       temp71111(1,2,2)=temp72211(1,1,1)
       temp72111(1,1,2)=temp72211(1,1,1)
       temp72111(1,2,1)=temp72211(1,1,1)
       temp72111(1,2,2)=temp72221(1,1,1)
       temp72211(1,1,2)=temp72221(1,1,1)
       temp72211(1,2,1)=temp72221(1,1,1)
       temp72211(1,2,2)=temp72222(1,1,1)
       temp72221(1,1,2)=temp72222(1,1,1)
       temp72221(1,2,1)=temp72222(1,1,1)
       temp72221(1,2,2)=temp72222(2,1,1)
       temp72222(1,1,2)=temp72222(2,1,1)
       temp72222(1,2,1)=temp72222(2,1,1)
       temp72222(1,2,2)=temp72222(2,2,1)
       temp72222(2,1,2)=temp72222(2,2,1)
c                Step3
       temp0000001(1)=I18Z*(aux0000001(1)+2*tempC3000000*F(4)+F(3)*temp0
     &  0001(1)-det3*temp00003(1,k,l))
       temp0000001(2)=I18Z*(aux0000001(2)+tempC3000000*F(5)+F(3)*temp000
     &  01(2)-det3*temp00003(2,k,l))
       temp00003(1,1,1)=I22Z*(aux00003(1,1,1)+6*F(4)*temp00002(1,1)+F(3)
     &  *temp003(1,1,1)-det3*temp00511(1,k,l)+24*temp0000001(1)*ZZ(k,1,l
     &  ,1))
       temp00003(2,1,1)=I22Z*(aux00003(2,1,1)+F(5)*temp00002(1,1)+4*F(4)
     &  *temp00002(2,1)+F(3)*temp003(2,1,1)-det3*temp00521(1,k,l)+8*(tem
     &  p0000001(2)*ZZ(k,1,l,1)+temp0000001(1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))
     &  ))
       temp00003(2,2,1)=I22Z*(aux00003(2,2,1)+2*F(5)*temp00002(2,1)+2*F(
     &  4)*temp00002(2,2)+F(3)*temp003(2,2,1)-det3*temp00522(1,k,l)+8*(t
     &  emp0000001(2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp0000001(1)*ZZ(k,2,l,
     &  2)))
       temp00003(2,2,2)=I22Z*(aux00003(2,2,2)+F(3)*temp003(2,2,2)-det3*t
     &  emp00522(2,k,l)+24*temp0000001(2)*ZZ(k,2,l,2)+temp00002(2,2)*(6*
     &  r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2)))
       temp00003(1,1,2)=temp00003(2,1,1)
       temp00003(1,2,1)=temp00003(2,1,1)
       temp00003(1,2,2)=temp00003(2,2,1)
       temp00003(2,1,2)=temp00003(2,2,1)
       temp00511(1,1,1)=I26Z*(aux00511(1,1,1)+10*F(4)*temp0041(1,1,1)+F(
     &  3)*temp511(1,1,1)-det3*temp71111(1,k,l)+80*temp00003(1,1,1)*ZZ(k
     &  ,1,l,1))
       temp00521(1,1,1)=I26Z*(aux00521(1,1,1)+F(5)*temp0041(1,1,1)+8*F(4
     &  )*temp0042(1,1,1)+F(3)*temp521(1,1,1)-det3*temp72111(1,k,l)+48*t
     &  emp00003(2,1,1)*ZZ(k,1,l,1)+16*temp00003(1,1,1)*(ZZ(k,1,l,2)+ZZ(
     &  k,2,l,1)))
       temp00522(1,1,1)=I26Z*(aux00522(1,1,1)+2*F(5)*temp0042(1,1,1)+6*F
     &  (4)*temp0042(2,1,1)+F(3)*temp522(1,1,1)-det3*temp72211(1,k,l)+24
     &  *(temp00003(2,2,1)*ZZ(k,1,l,1)+temp00003(2,1,1)*(ZZ(k,1,l,2)+ZZ(
     &  k,2,l,1)))+8*temp00003(1,1,1)*ZZ(k,2,l,2))
       temp00522(2,1,1)=I26Z*(aux00522(2,1,1)+4*F(4)*temp0042(2,2,1)+F(3
     &  )*temp522(2,1,1)-det3*temp72221(1,k,l)+8*temp00003(2,2,2)*ZZ(k,1
     &  ,l,1)+temp0042(2,1,1)*(6*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+12*r21*ZZ
     &  (k,2,l,2))+24*(temp00003(2,2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp00
     &  003(2,1,1)*ZZ(k,2,l,2)))
       temp00522(2,2,1)=I26Z*(aux00522(2,2,1)+4*F(5)*temp0042(2,2,1)+2*F
     &  (4)*temp0042(2,2,2)+F(3)*temp522(2,2,1)-det3*temp72222(1,k,l)+16
     &  *temp00003(2,2,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+48*temp00003(2,2,1)*
     &  ZZ(k,2,l,2))
       temp00522(2,2,2)=I26Z*(aux00522(2,2,2)+F(3)*temp522(2,2,2)-det3*t
     &  emp72222(2,k,l)+80*temp00003(2,2,2)*ZZ(k,2,l,2)+temp0042(2,2,2)*
     &  (10*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+20*r21*ZZ(k,2,l,2)))
       temp00511(1,1,2)=temp00521(1,1,1)
       temp00511(1,2,1)=temp00521(1,1,1)
       temp00511(1,2,2)=temp00522(1,1,1)
       temp00521(1,1,2)=temp00522(1,1,1)
       temp00521(1,2,1)=temp00522(1,1,1)
       temp00521(1,2,2)=temp00522(2,1,1)
       temp00522(1,1,2)=temp00522(2,1,1)
       temp00522(1,2,1)=temp00522(2,1,1)
       temp00522(1,2,2)=temp00522(2,2,1)
       temp00522(2,1,2)=temp00522(2,2,1)
       temp6111(1,1,1)=IX*(aux6111(1,1,1)+det3*temp71111(1,1,jj)+12*temp
     &  00511(1,1,1)*Z(jj,1))
       temp6211(1,1,1)=IX*(aux6211(1,1,1)+det3*temp72111(1,1,jj)+10*temp
     &  00521(1,1,1)*Z(jj,1)+2*temp00511(1,1,1)*Z(jj,2))
       temp6221(1,1,1)=IX*(aux6221(1,1,1)+det3*temp72211(1,1,jj)+8*temp0
     &  0522(1,1,1)*Z(jj,1)+4*temp00521(1,1,1)*Z(jj,2))
       temp6222(1,1,1)=IX*(aux6222(1,1,1)+det3*temp72221(1,1,jj)+6*(temp
     &  00522(2,1,1)*Z(jj,1)+temp00522(1,1,1)*Z(jj,2)))
       temp6222(2,1,1)=IX*(aux6222(2,1,1)+det3*temp72222(1,1,jj)+4*temp0
     &  0522(2,2,1)*Z(jj,1)+8*temp00522(2,1,1)*Z(jj,2))
       temp6222(2,2,1)=IX*(aux6222(2,2,1)+det3*temp72222(2,1,jj)+2*temp0
     &  0522(2,2,2)*Z(jj,1)+10*temp00522(2,2,1)*Z(jj,2))
       temp6222(2,2,2)=IX*(aux6222(2,2,2)+det3*temp72222(2,2,jj)+12*temp
     &  00522(2,2,2)*Z(jj,2))
       temp6111(1,1,2)=temp6211(1,1,1)
       temp6111(1,2,1)=temp6211(1,1,1)
       temp6111(1,2,2)=temp6221(1,1,1)
       temp6211(1,1,2)=temp6221(1,1,1)
       temp6211(1,2,1)=temp6221(1,1,1)
       temp6211(1,2,2)=temp6222(1,1,1)
       temp6221(1,1,2)=temp6222(1,1,1)
       temp6221(1,2,1)=temp6222(1,1,1)
       temp6221(1,2,2)=temp6222(2,1,1)
       temp6222(1,1,2)=temp6222(2,1,1)
       temp6222(1,2,1)=temp6222(2,1,1)
       temp6222(1,2,2)=temp6222(2,2,1)
       temp6222(2,1,2)=temp6222(2,2,1)
c                Step4
       tempC3000000=I14Z*(auxC3000000+tempC30000*F(3)-det3*temp00002(k,l
     &  ))
       temp00002(1,1)=I18Z*(aux00002(1,1)+4*F(4)*temp00001(1)+F(3)*temp0
     &  02(1,1)-det3*temp0041(1,k,l)+8*tempC3000000*ZZ(k,1,l,1))
       temp00002(2,1)=I18Z*(aux00002(2,1)+F(5)*temp00001(1)+2*F(4)*temp0
     &  0001(2)+F(3)*temp002(2,1)-det3*temp0042(1,k,l)+4*tempC3000000*(Z
     &  Z(k,1,l,2)+ZZ(k,2,l,1)))
       temp00002(2,2)=I18Z*(aux00002(2,2)+2*F(5)*temp00001(2)+F(3)*temp0
     &  02(2,2)-det3*temp0042(2,k,l)+8*tempC3000000*ZZ(k,2,l,2))
       temp00002(1,2)=temp00002(2,1)
       temp0041(1,1,1)=I22Z*(aux0041(1,1,1)+8*F(4)*temp003(1,1,1)+F(3)*t
     &  emp41(1,1,1)-det3*temp6111(1,k,l)+48*temp00002(1,1)*ZZ(k,1,l,1))
       temp0042(1,1,1)=I22Z*(aux0042(1,1,1)+F(5)*temp003(1,1,1)+6*F(4)*t
     &  emp003(2,1,1)+F(3)*temp42(1,1,1)-det3*temp6211(1,k,l)+24*temp000
     &  02(2,1)*ZZ(k,1,l,1)+12*temp00002(1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp0042(2,1,1)=I22Z*(aux0042(2,1,1)+2*F(5)*temp003(2,1,1)+4*F(4)
     &  *temp003(2,2,1)+F(3)*temp42(2,1,1)-det3*temp6221(1,k,l)+16*temp0
     &  0002(2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+8*(temp00002(2,2)*ZZ(k,1,l,1
     &  )+temp00002(1,1)*ZZ(k,2,l,2)))
       temp0042(2,2,1)=I22Z*(aux0042(2,2,1)+2*F(4)*temp003(2,2,2)+F(3)*t
     &  emp42(2,2,1)-det3*temp6222(1,k,l)+12*temp00002(2,2)*(ZZ(k,1,l,2)
     &  +ZZ(k,2,l,1))+24*temp00002(2,1)*ZZ(k,2,l,2)+temp003(2,2,1)*(6*r1
     &  0*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2)))
       temp0042(2,2,2)=I22Z*(aux0042(2,2,2)+4*F(5)*temp003(2,2,2)+F(3)*t
     &  emp42(2,2,2)-det3*temp6222(2,k,l)+48*temp00002(2,2)*ZZ(k,2,l,2))
       temp0041(1,1,2)=temp0042(1,1,1)
       temp0041(1,2,1)=temp0042(1,1,1)
       temp0041(1,2,2)=temp0042(2,1,1)
       temp0042(1,1,2)=temp0042(2,1,1)
       temp0042(1,2,1)=temp0042(2,1,1)
       temp0042(1,2,2)=temp0042(2,2,1)
       temp0042(2,1,2)=temp0042(2,2,1)
       temp511(1,1,1)=IX*(aux511(1,1,1)+det3*temp6111(1,1,jj)+10*temp004
     &  1(1,1,1)*Z(jj,1))
       temp521(1,1,1)=IX*(aux521(1,1,1)+det3*temp6211(1,1,jj)+8*temp0042
     &  (1,1,1)*Z(jj,1)+2*temp0041(1,1,1)*Z(jj,2))
       temp522(1,1,1)=IX*(aux522(1,1,1)+det3*temp6221(1,1,jj)+6*temp0042
     &  (2,1,1)*Z(jj,1)+4*temp0042(1,1,1)*Z(jj,2))
       temp522(2,1,1)=IX*(aux522(2,1,1)+det3*temp6222(1,1,jj)+4*temp0042
     &  (2,2,1)*Z(jj,1)+6*temp0042(2,1,1)*Z(jj,2))
       temp522(2,2,1)=IX*(aux522(2,2,1)+det3*temp6222(2,1,jj)+2*temp0042
     &  (2,2,2)*Z(jj,1)+8*temp0042(2,2,1)*Z(jj,2))
       temp522(2,2,2)=IX*(aux522(2,2,2)+det3*temp6222(2,2,jj)+10*temp004
     &  2(2,2,2)*Z(jj,2))
       temp511(1,1,2)=temp521(1,1,1)
       temp511(1,2,1)=temp521(1,1,1)
       temp511(1,2,2)=temp522(1,1,1)
       temp521(1,1,2)=temp522(1,1,1)
       temp521(1,2,1)=temp522(1,1,1)
       temp521(1,2,2)=temp522(2,1,1)
       temp522(1,1,2)=temp522(2,1,1)
       temp522(1,2,1)=temp522(2,1,1)
       temp522(1,2,2)=temp522(2,2,1)
       temp522(2,1,2)=temp522(2,2,1)
c                Step5
       temp00001(1)=I14Z*(aux00001(1)+2*tempC30000*F(4)+F(3)*temp001(1)-
     &  det3*temp003(1,k,l))
       temp00001(2)=I14Z*(aux00001(2)+tempC30000*F(5)+F(3)*temp001(2)-de
     &  t3*temp003(2,k,l))
       temp003(1,1,1)=I18Z*(aux003(1,1,1)+6*F(4)*temp002(1,1)+F(3)*temp3
     &  (1,1,1)-det3*temp511(1,k,l)+24*temp00001(1)*ZZ(k,1,l,1))
       temp003(2,1,1)=I18Z*(aux003(2,1,1)+F(5)*temp002(1,1)+4*F(4)*temp0
     &  02(2,1)+F(3)*temp3(2,1,1)-det3*temp521(1,k,l)+8*(temp00001(2)*ZZ
     &  (k,1,l,1)+temp00001(1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))))
       temp003(2,2,1)=I18Z*(aux003(2,2,1)+2*F(5)*temp002(2,1)+2*F(4)*tem
     &  p002(2,2)+F(3)*temp3(2,2,1)-det3*temp522(1,k,l)+8*(temp00001(2)*
     &  (ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp00001(1)*ZZ(k,2,l,2)))
       temp003(2,2,2)=I18Z*(aux003(2,2,2)+F(3)*temp3(2,2,2)-det3*temp522
     &  (2,k,l)+24*temp00001(2)*ZZ(k,2,l,2)+temp002(2,2)*(6*r10*(ZZ(k,1,
     &  l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2)))
       temp003(1,1,2)=temp003(2,1,1)
       temp003(1,2,1)=temp003(2,1,1)
       temp003(1,2,2)=temp003(2,2,1)
       temp003(2,1,2)=temp003(2,2,1)
       temp41(1,1,1)=IX*(aux41(1,1,1)+det3*temp511(1,1,jj)+8*temp003(1,1
     &  ,1)*Z(jj,1))
       temp42(1,1,1)=IX*(aux42(1,1,1)+det3*temp521(1,1,jj)+6*temp003(2,1
     &  ,1)*Z(jj,1)+2*temp003(1,1,1)*Z(jj,2))
       temp42(2,1,1)=IX*(aux42(2,1,1)+det3*temp522(1,1,jj)+4*(temp003(2,
     &  2,1)*Z(jj,1)+temp003(2,1,1)*Z(jj,2)))
       temp42(2,2,1)=IX*(aux42(2,2,1)+det3*temp522(2,1,jj)+2*temp003(2,2
     &  ,2)*Z(jj,1)+6*temp003(2,2,1)*Z(jj,2))
       temp42(2,2,2)=IX*(aux42(2,2,2)+det3*temp522(2,2,jj)+8*temp003(2,2
     &  ,2)*Z(jj,2))
       temp41(1,1,2)=temp42(1,1,1)
       temp41(1,2,1)=temp42(1,1,1)
       temp41(1,2,2)=temp42(2,1,1)
       temp42(1,1,2)=temp42(2,1,1)
       temp42(1,2,1)=temp42(2,1,1)
       temp42(1,2,2)=temp42(2,2,1)
       temp42(2,1,2)=temp42(2,2,1)
c                Step6
       tempC30000=I10Z*(auxC30000+tempC300*F(3)-det3*temp002(k,l))
       temp002(1,1)=I14Z*(aux002(1,1)+4*F(4)*temp001(1)+F(3)*temp2(1,1)-
     &  det3*temp41(1,k,l)+8*tempC30000*ZZ(k,1,l,1))
       temp002(2,1)=I14Z*(aux002(2,1)+F(5)*temp001(1)+2*F(4)*temp001(2)+
     &  F(3)*temp2(2,1)-det3*temp42(1,k,l)+4*tempC30000*(ZZ(k,1,l,2)+ZZ(
     &  k,2,l,1)))
       temp002(2,2)=I14Z*(aux002(2,2)+2*F(5)*temp001(2)+F(3)*temp2(2,2)-
     &  det3*temp42(2,k,l)+8*tempC30000*ZZ(k,2,l,2))
       temp002(1,2)=temp002(2,1)
       temp3(1,1,1)=IX*(aux3(1,1,1)+det3*temp41(1,1,jj)+6*temp002(1,1)*Z
     &  (jj,1))
       temp3(2,1,1)=IX*(aux3(2,1,1)+det3*temp42(1,1,jj)+4*temp002(2,1)*Z
     &  (jj,1)+2*temp002(1,1)*Z(jj,2))
       temp3(2,2,1)=IX*(aux3(2,2,1)+det3*temp42(2,1,jj)+2*temp002(2,2)*Z
     &  (jj,1)+4*temp002(2,1)*Z(jj,2))
       temp3(2,2,2)=IX*(aux3(2,2,2)+det3*temp42(2,2,jj)+6*temp002(2,2)*Z
     &  (jj,2))
       temp3(1,1,2)=temp3(2,1,1)
       temp3(1,2,1)=temp3(2,1,1)
       temp3(1,2,2)=temp3(2,2,1)
       temp3(2,1,2)=temp3(2,2,1)
c                Step7
       temp001(1)=I10Z*(aux001(1)+2*tempC300*F(4)+F(3)*temp1(1)-det3*tem
     &  p3(1,k,l))
       temp001(2)=I10Z*(aux001(2)+tempC300*F(5)+F(3)*temp1(2)-det3*temp3
     &  (2,k,l))
       temp2(1,1)=IX*(aux2(1,1)+det3*temp3(1,1,jj)+4*temp001(1)*Z(jj,1))
       temp2(2,1)=IX*(aux2(2,1)+det3*temp3(2,1,jj)+2*(temp001(2)*Z(jj,1)
     &  +temp001(1)*Z(jj,2)))
       temp2(2,2)=IX*(aux2(2,2)+det3*temp3(2,2,jj)+4*temp001(2)*Z(jj,2))
       temp2(1,2)=temp2(2,1)
c                Step8
       tempC300=I6Z*(auxC300+tempC30*F(3)-det3*temp2(k,l))
       temp1(1)=IX*(aux1(1)+det3*temp2(1,jj)+2*tempC300*Z(jj,1))
       temp1(2)=IX*(aux1(2)+det3*temp2(2,jj)+2*tempC300*Z(jj,2))
c                Step9
       tempC30=IX*(auxC30+det3*temp1(jj))

         ac=2
         accuracyCR(1,0,ac) = abs(C30     /tempC30      -1d0) 
         accuracyCR(1,1,ac) = abs(Cij(1,1)/temp1(1)      -1d0) 
         accuracyCR(2,1,ac) = abs(Cij(2,1)/temp1(2)      -1d0)
         accuracyCR(1,2,ac) = abs(Cij(1,2)/temp2(1,1)    -1d0)
         accuracyCR(2,2,ac) = abs(Cij(2,2)/temp2(2,2)    -1d0)
         accuracyCR(3,2,ac) = abs(Cij(3,2)/temp2(2,1)    -1d0)
         accuracyCR(4,2,ac) = abs(Cij(4,2)/tempC300      -1d0)
         accuracyCR(1,3,ac) = abs(Cij(1,3)/temp3(1,1,1)  -1d0)
         accuracyCR(2,3,ac) = abs(Cij(2,3)/temp3(2,2,2)  -1d0)
         accuracyCR(3,3,ac) = abs(Cij(3,3)/temp3(2,1,1)  -1d0)
         accuracyCR(4,3,ac) = abs(Cij(4,3)/temp3(2,2,1)  -1d0)
         accuracyCR(5,3,ac) = abs(Cij(5,3)/temp001(1)    -1d0)
         accuracyCR(6,3,ac) = abs(Cij(6,3)/temp001(2)    -1d0)
         accuracyCR(1,4,ac) = abs(Cij(1,4)/temp41(1,1,1) -1d0)
         accuracyCR(2,4,ac) = abs(Cij(2,4)/temp42(2,2,2) -1d0)
         accuracyCR(3,4,ac) = abs(Cij(3,4)/temp42(1,1,1) -1d0)
         accuracyCR(4,4,ac) = abs(Cij(4,4)/temp42(2,1,1) -1d0)
         accuracyCR(5,4,ac) = abs(Cij(5,4)/temp42(2,2,1) -1d0)
         accuracyCR(6,4,ac) = abs(Cij(6,4)/temp002(1,1)  -1d0)
         accuracyCR(7,4,ac) = abs(Cij(7,4)/temp002(2,2)  -1d0)
         accuracyCR(8,4,ac) = abs(Cij(8,4)/temp002(2,1)  -1d0)
         accuracyCR(9,4,ac) = abs(Cij(9,4)/tempC30000    -1d0)


      DO I1=0,4
           accuracyC(i1,ac)=accuracyCR(1,i1,ac)
            DO I2=2,INDEX(I1)  
       if(accuracyCR(i2,i1,ac).gt.accuracyC(i1,ac)) then
          accuracyC(i1,ac)=accuracyCR(i2,i1,ac)
       endif
          enddo
        enddo

c          if(accuracyC(4,ac).lt.1d-16) goto 500

           C30=tempC30
           Cij(1,1)=temp1(1)
           Cij(2,1)=temp1(2)
           Cij(1,2)=temp2(1,1)
           Cij(2,2)=temp2(2,2)
           Cij(3,2)=temp2(2,1)
           Cij(4,2)=tempC300
           Cij(1,3)=temp3(1,1,1)
           Cij(2,3)=temp3(2,2,2)
           Cij(3,3)=temp3(2,1,1)
           Cij(4,3)=temp3(2,2,1)
           Cij(5,3)=temp001(1)
           Cij(6,3)=temp001(2)
           Cij(1,4)=temp41(1,1,1)
           Cij(2,4)=temp42(2,2,2)
           Cij(3,4)=temp42(1,1,1)
           Cij(4,4)=temp42(2,1,1)
           Cij(5,4)=temp42(2,2,1)
           Cij(6,4)=temp002(1,1)
           Cij(7,4)=temp002(2,2)
           Cij(8,4)=temp002(2,1)
           Cij(9,4)=tempC30000


       if(order.eq.8) goto 500
c                Iteration9
c                Step1
       S20000000000=2*B23(5,8)
       S20000000021(1)=2*B23(4,6)
       S20000000021(2)=-2*B23(4,7)
       S20000000022(1)=-2*B23(4,7)
       S20000000022(2)=2*B23(4,8)
       S2h0000000021(1)=B13(5,9)+B23(5,8)
       S2h0000000021(2)=B13(5,9)-B23(5,9)
       S2h0000000022(1)=B12(5,9)-B13(5,9)
       S2h0000000022(2)=-B13(5,9)
       auxC30000000000=-(F(1)*S2h000000001(1))-F(2)*S2h000000001(2)+S2h0
     &  000000021(k)*Z(1,l)+S2h0000000022(k)*Z(2,l)+(Inv100+S20000000000
     &  -S2h0000000021(1)-S2h0000000022(2))*Z(k,l)
       tempC30000000000=I22Z*(auxC30000000000+tempC300000000*F(3))
       S2h0000004111(1)=B13(4,9)+B23(4,6)
       S2h0000004111(2)=B13(4,9)-B23(4,7)
       S2h0000004121(1)=B13(4,9)-B23(4,7)
       S2h0000004121(2)=B13(4,9)+B23(4,8)
       S2h0000004122(1)=B13(4,9)+B23(4,8)
       S2h0000004122(2)=B13(4,9)-B23(4,9)
       S2h0000004211(1)=B12(4,9)-B13(4,9)
       S2h0000004211(2)=-B13(4,9)
       S2h0000004221(1)=-B13(4,9)
       S2h0000004221(2)=-B13(4,9)
       S2h0000004222(1)=-B13(4,9)
       S2h0000004222(2)=-B13(4,9)
       aux000000002(1,1)=-(F(1)*S2h000000311(1))-F(2)*S2h000000321(1)+S2
     &  h0000004111(k)*Z(1,l)+S2h0000004211(k)*Z(2,l)+(Inv8011+S20000000
     &  021(1)-S2h0000004111(1)-S2h0000004221(1))*Z(k,l)-4*(S2h000000002
     &  1(1)*ZZ(k,1,l,1)+S2h0000000022(1)*ZZ(k,1,l,2))
       aux000000002(2,1)=-(F(1)*S2h000000312(1))-F(2)*S2h000000322(1)+S2
     &  h0000004121(k)*Z(1,l)+S2h0000004221(k)*Z(2,l)+(Inv8021+S20000000
     &  022(1)-S2h0000004121(1)-S2h0000004222(1))*Z(k,l)-2*(S2h000000002
     &  1(2)*ZZ(k,1,l,1)+S2h0000000022(2)*ZZ(k,1,l,2))-2*S2h0000000021(1
     &  )*ZZ(k,2,l,1)-2*S2h0000000022(1)*ZZ(k,2,l,2)
       aux000000002(2,2)=-(F(1)*S2h000000312(2))-F(2)*S2h000000322(2)+S2
     &  h0000004122(k)*Z(1,l)+S2h0000004222(k)*Z(2,l)+(Inv8022+S20000000
     &  022(2)-S2h0000004122(1)-S2h0000004222(2))*Z(k,l)-4*(S2h000000002
     &  1(2)*ZZ(k,2,l,1)+S2h0000000022(2)*ZZ(k,2,l,2))
       temp000000002(1,1)=I26Z*(aux000000002(1,1)+4*F(4)*temp000000001(1
     &  )+F(3)*temp0000002(1,1)+8*tempC30000000000*ZZ(k,1,l,1))
       temp000000002(2,1)=I26Z*(aux000000002(2,1)+F(5)*temp000000001(1)+
     &  2*F(4)*temp000000001(2)+F(3)*temp0000002(2,1)+4*tempC30000000000
     &  *(ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp000000002(2,2)=I26Z*(aux000000002(2,2)+2*F(5)*temp000000001(2
     &  )+F(3)*temp0000002(2,2)+8*tempC30000000000*ZZ(k,2,l,2))
       temp000000002(1,2)=temp000000002(2,1)
       S20000004111(1)=2*B23(3,4)
       S20000004111(2)=-2*B23(3,5)
       S20000004121(1)=-2*B23(3,5)
       S20000004121(2)=2*B23(3,6)
       S20000004122(1)=2*B23(3,6)
       S20000004122(2)=-2*B23(3,7)
       S20000004211(1)=-2*B23(3,5)
       S20000004211(2)=2*B23(3,6)
       S20000004221(1)=2*B23(3,6)
       S20000004221(2)=-2*B23(3,7)
       S20000004222(1)=-2*B23(3,7)
       S20000004222(2)=2*B23(3,8)
       S20000611111(1)=2*B23(2,2)
       S20000611111(2)=-2*B23(2,3)
       S20000612111(1)=-2*B23(2,3)
       S20000612111(2)=2*B23(2,4)
       S20000612211(1)=2*B23(2,4)
       S20000612211(2)=-2*B23(2,5)
       S20000612221(1)=-2*B23(2,5)
       S20000612221(2)=2*B23(2,6)
       S20000612222(1)=2*B23(2,6)
       S20000612222(2)=-2*B23(2,7)
       S20000621111(1)=-2*B23(2,3)
       S20000621111(2)=2*B23(2,4)
       S20000622111(1)=2*B23(2,4)
       S20000622111(2)=-2*B23(2,5)
       S20000622211(1)=-2*B23(2,5)
       S20000622211(2)=2*B23(2,6)
       S20000622221(1)=2*B23(2,6)
       S20000622221(2)=-2*B23(2,7)
       S20000622222(1)=-2*B23(2,7)
       S20000622222(2)=2*B23(2,8)
       S2h0000611111(1)=B13(3,9)+B23(3,4)
       S2h0000611111(2)=B13(3,9)-B23(3,5)
       S2h0000612111(1)=B13(3,9)-B23(3,5)
       S2h0000612111(2)=B13(3,9)+B23(3,6)
       S2h0000612211(1)=B13(3,9)+B23(3,6)
       S2h0000612211(2)=B13(3,9)-B23(3,7)
       S2h0000612221(1)=B13(3,9)-B23(3,7)
       S2h0000612221(2)=B13(3,9)+B23(3,8)
       S2h0000612222(1)=B13(3,9)+B23(3,8)
       S2h0000612222(2)=B13(3,9)-B23(3,9)
       S2h0000621111(1)=B12(3,9)-B13(3,9)
       S2h0000621111(2)=-B13(3,9)
       S2h0000622111(1)=-B13(3,9)
       S2h0000622111(2)=-B13(3,9)
       S2h0000622211(1)=-B13(3,9)
       S2h0000622211(2)=-B13(3,9)
       S2h0000622221(1)=-B13(3,9)
       S2h0000622221(2)=-B13(3,9)
       S2h0000622222(1)=-B13(3,9)
       S2h0000622222(2)=-B13(3,9)
       aux00000041(1,1,1)=-(F(1)*S2h000051111(1))-F(2)*S2h000052111(1)+S
     &  2h0000611111(k)*Z(1,l)+S2h0000621111(k)*Z(2,l)+(Inv601111+S20000
     &  004111(1)-S2h0000611111(1)-S2h0000622111(1))*Z(k,l)-8*(S2h000000
     &  4111(1)*ZZ(k,1,l,1)+S2h0000004211(1)*ZZ(k,1,l,2))
       aux00000042(1,1,1)=-(F(1)*S2h000051211(1))-F(2)*S2h000052211(1)+S
     &  2h0000612111(k)*Z(1,l)+S2h0000622111(k)*Z(2,l)+(Inv602111+S20000
     &  004211(1)-S2h0000612111(1)-S2h0000622211(1))*Z(k,l)-6*(S2h000000
     &  4121(1)*ZZ(k,1,l,1)+S2h0000004221(1)*ZZ(k,1,l,2))-2*S2h000000411
     &  1(1)*ZZ(k,2,l,1)-2*S2h0000004211(1)*ZZ(k,2,l,2)
       aux00000042(2,1,1)=-(F(1)*S2h000051221(1))-F(2)*S2h000052221(1)+S
     &  2h0000612211(k)*Z(1,l)+S2h0000622211(k)*Z(2,l)+(Inv602211+S20000
     &  004221(1)-S2h0000612211(1)-S2h0000622221(1))*Z(k,l)-4*(S2h000000
     &  4122(1)*ZZ(k,1,l,1)+S2h0000004222(1)*ZZ(k,1,l,2))-4*S2h000000412
     &  1(1)*ZZ(k,2,l,1)-4*S2h0000004221(1)*ZZ(k,2,l,2)
       aux00000042(2,2,1)=-(F(1)*S2h000051222(1))-F(2)*S2h000052222(1)+S
     &  2h0000612221(k)*Z(1,l)+S2h0000622221(k)*Z(2,l)+(Inv602221+S20000
     &  004222(1)-S2h0000612221(1)-S2h0000622222(1))*Z(k,l)-2*(S2h000000
     &  4122(2)*ZZ(k,1,l,1)+S2h0000004222(2)*ZZ(k,1,l,2))-6*S2h000000412
     &  2(1)*ZZ(k,2,l,1)-6*S2h0000004222(1)*ZZ(k,2,l,2)
       aux00000042(2,2,2)=-(F(1)*S2h000051222(2))-F(2)*S2h000052222(2)+S
     &  2h0000612222(k)*Z(1,l)+S2h0000622222(k)*Z(2,l)+(Inv602222+S20000
     &  004222(2)-S2h0000612222(1)-S2h0000622222(2))*Z(k,l)-8*(S2h000000
     &  4122(2)*ZZ(k,2,l,1)+S2h0000004222(2)*ZZ(k,2,l,2))
       temp00000041(1,1,1)=I30Z*(aux00000041(1,1,1)+8*F(4)*temp0000003(1
     &  ,1,1)+F(3)*temp000041(1,1,1)+48*temp000000002(1,1)*ZZ(k,1,l,1))
       temp00000042(1,1,1)=I30Z*(aux00000042(1,1,1)+F(5)*temp0000003(1,1
     &  ,1)+6*F(4)*temp0000003(2,1,1)+F(3)*temp000042(1,1,1)+24*temp0000
     &  00002(2,1)*ZZ(k,1,l,1)+12*temp000000002(1,1)*(ZZ(k,1,l,2)+ZZ(k,2
     &  ,l,1)))
       temp00000042(2,1,1)=I30Z*(aux00000042(2,1,1)+2*F(5)*temp0000003(2
     &  ,1,1)+4*F(4)*temp0000003(2,2,1)+F(3)*temp000042(2,1,1)+16*temp00
     &  0000002(2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+8*(temp000000002(2,2)*ZZ(
     &  k,1,l,1)+temp000000002(1,1)*ZZ(k,2,l,2)))
       temp00000042(2,2,1)=I30Z*(aux00000042(2,2,1)+2*F(4)*temp0000003(2
     &  ,2,2)+F(3)*temp000042(2,2,1)+12*temp000000002(2,2)*(ZZ(k,1,l,2)+
     &  ZZ(k,2,l,1))+24*temp000000002(2,1)*ZZ(k,2,l,2)+temp0000003(2,2,1
     &  )*(6*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2)))
       temp00000042(2,2,2)=I30Z*(aux00000042(2,2,2)+4*F(5)*temp0000003(2
     &  ,2,2)+F(3)*temp000042(2,2,2)+48*temp000000002(2,2)*ZZ(k,2,l,2))
       temp00000041(1,1,2)=temp00000042(1,1,1)
       temp00000041(1,2,1)=temp00000042(1,1,1)
       temp00000041(1,2,2)=temp00000042(2,1,1)
       temp00000042(1,1,2)=temp00000042(2,1,1)
       temp00000042(1,2,1)=temp00000042(2,1,1)
       temp00000042(1,2,2)=temp00000042(2,2,1)
       temp00000042(2,1,2)=temp00000042(2,2,1)
       S2h0081111111(1)=B13(2,9)+B23(2,2)
       S2h0081111111(2)=B13(2,9)-B23(2,3)
       S2h0081211111(1)=B13(2,9)-B23(2,3)
       S2h0081211111(2)=B13(2,9)+B23(2,4)
       S2h0081221111(1)=B13(2,9)+B23(2,4)
       S2h0081221111(2)=B13(2,9)-B23(2,5)
       S2h0081222111(1)=B13(2,9)-B23(2,5)
       S2h0081222111(2)=B13(2,9)+B23(2,6)
       S2h0081222211(1)=B13(2,9)+B23(2,6)
       S2h0081222211(2)=B13(2,9)-B23(2,7)
       S2h0081222221(1)=B13(2,9)-B23(2,7)
       S2h0081222221(2)=B13(2,9)+B23(2,8)
       S2h0081222222(1)=B13(2,9)+B23(2,8)
       S2h0081222222(2)=B13(2,9)-B23(2,9)
       S2h0082111111(1)=B12(2,9)-B13(2,9)
       S2h0082111111(2)=-B13(2,9)
       S2h0082211111(1)=-B13(2,9)
       S2h0082211111(2)=-B13(2,9)
       S2h0082221111(1)=-B13(2,9)
       S2h0082221111(2)=-B13(2,9)
       S2h0082222111(1)=-B13(2,9)
       S2h0082222111(2)=-B13(2,9)
       S2h0082222211(1)=-B13(2,9)
       S2h0082222211(2)=-B13(2,9)
       S2h0082222221(1)=-B13(2,9)
       S2h0082222221(2)=-B13(2,9)
       S2h0082222222(1)=-B13(2,9)
       S2h0082222222(2)=-B13(2,9)
       aux00006111(1,1,1)=-(F(1)*S2h007111111(1))-F(2)*S2h007211111(1)+S
     &  2h0081111111(k)*Z(1,l)+S2h0082111111(k)*Z(2,l)+(Inv40111111+S2
     &  0000611111(1)-S2h0081111111(1)-S2h0082211111(1))*Z(k,l)-12*(S2h0
     &  000611111(1)*ZZ(k,1,l,1)+S2h0000621111(1)*ZZ(k,1,l,2))
       aux00006211(1,1,1)=-(F(1)*S2h007121111(1))-F(2)*S2h007221111(1)+S
     &  2h0081211111(k)*Z(1,l)+S2h0082211111(k)*Z(2,l)+(Inv40211111+S2
     &  0000621111(1)-S2h0081211111(1)-S2h0082221111(1))*Z(k,l)-10*(S2h0
     &  000612111(1)*ZZ(k,1,l,1)+S2h0000622111(1)*ZZ(k,1,l,2))-2*S2h0000
     &  611111(1)*ZZ(k,2,l,1)-2*S2h0000621111(1)*ZZ(k,2,l,2)
       aux00006221(1,1,1)=-(F(1)*S2h007122111(1))-F(2)*S2h007222111(1)+S
     &  2h0081221111(k)*Z(1,l)+S2h0082221111(k)*Z(2,l)+(Inv40221111+S2
     &  0000622111(1)-S2h0081221111(1)-S2h0082222111(1))*Z(k,l)-8*(S2h00
     &  00612211(1)*ZZ(k,1,l,1)+S2h0000622211(1)*ZZ(k,1,l,2))-4*S2h00006
     &  12111(1)*ZZ(k,2,l,1)-4*S2h0000622111(1)*ZZ(k,2,l,2)
       aux00006222(1,1,1)=-(F(1)*S2h007122211(1))-F(2)*S2h007222211(1)+S
     &  2h0081222111(k)*Z(1,l)+S2h0082222111(k)*Z(2,l)+(Inv40222111+S2
     &  0000622211(1)-S2h0081222111(1)-S2h0082222211(1))*Z(k,l)-6*(S2h00
     &  00612221(1)*ZZ(k,1,l,1)+S2h0000622221(1)*ZZ(k,1,l,2))-6*S2h00006
     &  12211(1)*ZZ(k,2,l,1)-6*S2h0000622211(1)*ZZ(k,2,l,2)
       aux00006222(2,1,1)=-(F(1)*S2h007122221(1))-F(2)*S2h007222221(1)+S
     &  2h0081222211(k)*Z(1,l)+S2h0082222211(k)*Z(2,l)+(Inv40222211+S2
     &  0000622221(1)-S2h0081222211(1)-S2h0082222221(1))*Z(k,l)-4*(S2h00
     &  00612222(1)*ZZ(k,1,l,1)+S2h0000622222(1)*ZZ(k,1,l,2))-8*S2h00006
     &  12221(1)*ZZ(k,2,l,1)-8*S2h0000622221(1)*ZZ(k,2,l,2)
       aux00006222(2,2,1)=-(F(1)*S2h007122222(1))-F(2)*S2h007222222(1)+S
     &  2h0081222221(k)*Z(1,l)+S2h0082222221(k)*Z(2,l)+(Inv40222221+S2
     &  0000622222(1)-S2h0081222221(1)-S2h0082222222(1))*Z(k,l)-2*(S2h00
     &  00612222(2)*ZZ(k,1,l,1)+S2h0000622222(2)*ZZ(k,1,l,2))-10*S2h0000
     &  612222(1)*ZZ(k,2,l,1)-10*S2h0000622222(1)*ZZ(k,2,l,2)
       aux00006222(2,2,2)=-(F(1)*S2h007122222(2))-F(2)*S2h007222222(2)+S
     &  2h0081222222(k)*Z(1,l)+S2h0082222222(k)*Z(2,l)+(Inv40222222+S2
     &  0000622222(2)-S2h0081222222(1)-S2h0082222222(2))*Z(k,l)-12*(S2h0
     &  000612222(2)*ZZ(k,2,l,1)+S2h0000622222(2)*ZZ(k,2,l,2))
       temp00006111(1,1,1)=I34Z*(aux00006111(1,1,1)+12*F(4)*temp0000511(
     &  1,1,1)+F(3)*temp006111(1,1,1)+120*temp00000041(1,1,1)*ZZ(k,1,l,1
     &  ))
       temp00006211(1,1,1)=I34Z*(aux00006211(1,1,1)+F(5)*temp0000511(1,1
     &  ,1)+10*F(4)*temp0000521(1,1,1)+F(3)*temp006211(1,1,1)+80*temp000
     &  00042(1,1,1)*ZZ(k,1,l,1)+20*temp00000041(1,1,1)*(ZZ(k,1,l,2)+ZZ(
     &  k,2,l,1)))
       temp00006221(1,1,1)=I34Z*(aux00006221(1,1,1)+2*F(5)*temp0000521(1
     &  ,1,1)+F(3)*temp006221(1,1,1)+48*temp00000042(2,1,1)*ZZ(k,1,l,1)+
     &  32*temp00000042(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+8*(F(4)*temp000
     &  0522(1,1,1)+temp00000041(1,1,1)*ZZ(k,2,l,2)))
       temp00006222(1,1,1)=I34Z*(aux00006222(1,1,1)+6*F(4)*temp0000522(2
     &  ,1,1)+F(3)*temp006222(1,1,1)+36*temp00000042(2,1,1)*(ZZ(k,1,l,2)
     &  +ZZ(k,2,l,1))+temp0000522(1,1,1)*(6*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1)
     &  )+12*r21*ZZ(k,2,l,2))+24*(temp00000042(2,2,1)*ZZ(k,1,l,1)+temp00
     &  000042(1,1,1)*ZZ(k,2,l,2)))
       temp00006222(2,1,1)=I34Z*(aux00006222(2,1,1)+4*F(5)*temp0000522(2
     &  ,1,1)+4*F(4)*temp0000522(2,2,1)+F(3)*temp006222(2,1,1)+8*temp000
     &  00042(2,2,2)*ZZ(k,1,l,1)+32*temp00000042(2,2,1)*(ZZ(k,1,l,2)+ZZ(
     &  k,2,l,1))+48*temp00000042(2,1,1)*ZZ(k,2,l,2))
       temp00006222(2,2,1)=I34Z*(aux00006222(2,2,1)+2*F(4)*temp0000522(2
     &  ,2,2)+F(3)*temp006222(2,2,1)+20*temp00000042(2,2,2)*(ZZ(k,1,l,2)
     &  +ZZ(k,2,l,1))+80*temp00000042(2,2,1)*ZZ(k,2,l,2)+temp0000522(2,2
     &  ,1)*(10*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+20*r21*ZZ(k,2,l,2)))
       temp00006222(2,2,2)=I34Z*(aux00006222(2,2,2)+F(3)*temp006222(2,2,
     &  2)+120*temp00000042(2,2,2)*ZZ(k,2,l,2)+temp0000522(2,2,2)*(12*r1
     &  0*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+24*r21*ZZ(k,2,l,2)))
       temp00006111(1,1,2)=temp00006211(1,1,1)
       temp00006111(1,2,1)=temp00006211(1,1,1)
       temp00006111(1,2,2)=temp00006221(1,1,1)
       temp00006211(1,1,2)=temp00006221(1,1,1)
       temp00006211(1,2,1)=temp00006221(1,1,1)
       temp00006211(1,2,2)=temp00006222(1,1,1)
       temp00006221(1,1,2)=temp00006222(1,1,1)
       temp00006221(1,2,1)=temp00006222(1,1,1)
       temp00006221(1,2,2)=temp00006222(2,1,1)
       temp00006222(1,1,2)=temp00006222(2,1,1)
       temp00006222(1,2,1)=temp00006222(2,1,1)
       temp00006222(1,2,2)=temp00006222(2,2,1)
       temp00006222(2,1,2)=temp00006222(2,2,1)
       S20081111111(1)=2*B023
       S20081111111(2)=-2*B23(1,1)
       S20081211111(1)=-2*B23(1,1)
       S20081211111(2)=2*B23(1,2)
       S20081221111(1)=2*B23(1,2)
       S20081221111(2)=-2*B23(1,3)
       S20081222111(1)=-2*B23(1,3)
       S20081222111(2)=2*B23(1,4)
       S20081222211(1)=2*B23(1,4)
       S20081222211(2)=-2*B23(1,5)
       S20081222221(1)=-2*B23(1,5)
       S20081222221(2)=2*B23(1,6)
       S20081222222(1)=2*B23(1,6)
       S20081222222(2)=-2*B23(1,7)
       S20082111111(1)=-2*B23(1,1)
       S20082111111(2)=2*B23(1,2)
       S20082211111(1)=2*B23(1,2)
       S20082211111(2)=-2*B23(1,3)
       S20082221111(1)=-2*B23(1,3)
       S20082221111(2)=2*B23(1,4)
       S20082222111(1)=2*B23(1,4)
       S20082222111(2)=-2*B23(1,5)
       S20082222211(1)=-2*B23(1,5)
       S20082222211(2)=2*B23(1,6)
       S20082222221(1)=2*B23(1,6)
       S20082222221(2)=-2*B23(1,7)
       S20082222222(1)=-2*B23(1,7)
       S20082222222(2)=2*B23(1,8)
       S210111111111(1)=B023+B13(1,9)
       S210111111111(2)=B13(1,9)-B23(1,1)
       S210121111111(1)=B13(1,9)-B23(1,1)
       S210121111111(2)=B13(1,9)+B23(1,2)
       S210122111111(1)=B13(1,9)+B23(1,2)
       S210122111111(2)=B13(1,9)-B23(1,3)
       S210122211111(1)=B13(1,9)-B23(1,3)
       S210122211111(2)=B13(1,9)+B23(1,4)
       S210122221111(1)=B13(1,9)+B23(1,4)
       S210122221111(2)=B13(1,9)-B23(1,5)
       S210122222111(1)=B13(1,9)-B23(1,5)
       S210122222111(2)=B13(1,9)+B23(1,6)
       S210122222211(1)=B13(1,9)+B23(1,6)
       S210122222211(2)=B13(1,9)-B23(1,7)
       S210122222221(1)=B13(1,9)-B23(1,7)
       S210122222221(2)=B13(1,9)+B23(1,8)
       S210122222222(1)=B13(1,9)+B23(1,8)
       S210122222222(2)=B13(1,9)-B23(1,9)
       S210211111111(1)=B12(1,9)-B13(1,9)
       S210211111111(2)=-B13(1,9)
       S210221111111(1)=-B13(1,9)
       S210221111111(2)=-B13(1,9)
       S210222111111(1)=-B13(1,9)
       S210222111111(2)=-B13(1,9)
       S210222211111(1)=-B13(1,9)
       S210222211111(2)=-B13(1,9)
       S210222221111(1)=-B13(1,9)
       S210222221111(2)=-B13(1,9)
       S210222222111(1)=-B13(1,9)
       S210222222111(2)=-B13(1,9)
       S210222222211(1)=-B13(1,9)
       S210222222211(2)=-B13(1,9)
       S210222222221(1)=-B13(1,9)
       S210222222221(2)=-B13(1,9)
       S210222222222(1)=-B13(1,9)
       S210222222222(2)=-B13(1,9)
       aux00811111(1,1,1)=-(F(1)*S2911111111(1))-F(2)*S2921111111(1)+S21
     &  0111111111(k)*Z(1,l)+S210211111111(k)*Z(2,l)+(Inv2011111111+S200
     &  81111111(1)-S210111111111(1)-S210221111111(1))*Z(k,l)-16*(S2h008
     &  1111111(1)*ZZ(k,1,l,1)+S2h0082111111(1)*ZZ(k,1,l,2))
       aux00821111(1,1,1)=-(F(1)*S2912111111(1))-F(2)*S2922111111(1)+S21
     &  0121111111(k)*Z(1,l)+S210221111111(k)*Z(2,l)+(Inv2021111111+S200
     &  82111111(1)-S210121111111(1)-S210222111111(1))*Z(k,l)-14*(S2h008
     &  1211111(1)*ZZ(k,1,l,1)+S2h0082211111(1)*ZZ(k,1,l,2))-2*(S2h00811
     &  11111(1)*ZZ(k,2,l,1)+S2h0082111111(1)*ZZ(k,2,l,2))
       aux00822111(1,1,1)=-(F(1)*S2912211111(1))-F(2)*S2922211111(1)+S21
     &  0122111111(k)*Z(1,l)+S210222111111(k)*Z(2,l)+(Inv2022111111+S200
     &  82211111(1)-S210122111111(1)-S210222211111(1))*Z(k,l)-12*(S2h008
     &  1221111(1)*ZZ(k,1,l,1)+S2h0082221111(1)*ZZ(k,1,l,2))-4*(S2h00812
     &  11111(1)*ZZ(k,2,l,1)+S2h0082211111(1)*ZZ(k,2,l,2))
       aux00822211(1,1,1)=-(F(1)*S2912221111(1))-F(2)*S2922221111(1)+S21
     &  0122211111(k)*Z(1,l)+S210222211111(k)*Z(2,l)+(Inv2022211111+S200
     &  82221111(1)-S210122211111(1)-S210222221111(1))*Z(k,l)-10*(S2h008
     &  1222111(1)*ZZ(k,1,l,1)+S2h0082222111(1)*ZZ(k,1,l,2))-6*(S2h00812
     &  21111(1)*ZZ(k,2,l,1)+S2h0082221111(1)*ZZ(k,2,l,2))
       aux00822221(1,1,1)=-(F(1)*S2912222111(1))-F(2)*S2922222111(1)+S21
     &  0122221111(k)*Z(1,l)+S210222221111(k)*Z(2,l)+(Inv2022221111+S200
     &  82222111(1)-S210122221111(1)-S210222222111(1))*Z(k,l)-8*(S2h0081
     &  222211(1)*ZZ(k,1,l,1)+S2h0082222211(1)*ZZ(k,1,l,2)+S2h0081222111
     &  (1)*ZZ(k,2,l,1))-8*S2h0082222111(1)*ZZ(k,2,l,2)
       aux00822222(1,1,1)=-(F(1)*S2912222211(1))-F(2)*S2922222211(1)+S21
     &  0122222111(k)*Z(1,l)+S210222222111(k)*Z(2,l)+(Inv2022222111+S200
     &  82222211(1)-S210122222111(1)-S210222222211(1))*Z(k,l)-6*(S2h0081
     &  222221(1)*ZZ(k,1,l,1)+S2h0082222221(1)*ZZ(k,1,l,2))-10*(S2h00812
     &  22211(1)*ZZ(k,2,l,1)+S2h0082222211(1)*ZZ(k,2,l,2))
       aux00822222(2,1,1)=-(F(1)*S2912222221(1))-F(2)*S2922222221(1)+S21
     &  0122222211(k)*Z(1,l)+S210222222211(k)*Z(2,l)+(Inv2022222211+S200
     &  82222221(1)-S210122222211(1)-S210222222221(1))*Z(k,l)-4*(S2h0081
     &  222222(1)*ZZ(k,1,l,1)+S2h0082222222(1)*ZZ(k,1,l,2))-12*(S2h00812
     &  22221(1)*ZZ(k,2,l,1)+S2h0082222221(1)*ZZ(k,2,l,2))
       aux00822222(2,2,1)=-(F(1)*S2912222222(1))-F(2)*S2922222222(1)+S21
     &  0122222221(k)*Z(1,l)+S210222222221(k)*Z(2,l)+(Inv2022222221+S200
     &  82222222(1)-S210122222221(1)-S210222222222(1))*Z(k,l)-2*(S2h0081
     &  222222(2)*ZZ(k,1,l,1)+S2h0082222222(2)*ZZ(k,1,l,2))-14*(S2h00812
     &  22222(1)*ZZ(k,2,l,1)+S2h0082222222(1)*ZZ(k,2,l,2))
       aux00822222(2,2,2)=-(F(1)*S2912222222(2))-F(2)*S2922222222(2)+S21
     &  0122222222(k)*Z(1,l)+S210222222222(k)*Z(2,l)+(Inv2022222222+S200
     &  82222222(2)-S210122222222(1)-S210222222222(2))*Z(k,l)-16*(S2h008
     &  1222222(2)*ZZ(k,2,l,1)+S2h0082222222(2)*ZZ(k,2,l,2))
       temp00811111(1,1,1)=I38Z*(aux00811111(1,1,1)+16*F(4)*temp0071111(
     &  1,1,1)+F(3)*temp811111(1,1,1)+224*temp00006111(1,1,1)*ZZ(k,1,l,1
     &  ))
       temp00821111(1,1,1)=I38Z*(aux00821111(1,1,1)+F(5)*temp0071111(1,1
     &  ,1)+14*F(4)*temp0072111(1,1,1)+F(3)*temp821111(1,1,1)+168*temp00
     &  006211(1,1,1)*ZZ(k,1,l,1)+28*temp00006111(1,1,1)*(ZZ(k,1,l,2)+ZZ
     &  (k,2,l,1)))
       temp00822111(1,1,1)=I38Z*(aux00822111(1,1,1)+2*F(5)*temp0072111(1
     &  ,1,1)+12*F(4)*temp0072211(1,1,1)+F(3)*temp822111(1,1,1)+120*temp
     &  00006221(1,1,1)*ZZ(k,1,l,1)+48*temp00006211(1,1,1)*(ZZ(k,1,l,2)+
     &  ZZ(k,2,l,1))+8*temp00006111(1,1,1)*ZZ(k,2,l,2))
       temp00822211(1,1,1)=I38Z*(aux00822211(1,1,1)+10*F(4)*temp0072221(
     &  1,1,1)+F(3)*temp822211(1,1,1)+80*temp00006222(1,1,1)*ZZ(k,1,l,1)
     &  +60*temp00006221(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+24*temp0000621
     &  1(1,1,1)*ZZ(k,2,l,2)+temp0072211(1,1,1)*(6*r10*(ZZ(k,1,l,2)+ZZ(k
     &  ,2,l,1))+12*r21*ZZ(k,2,l,2)))
       temp00822221(1,1,1)=I38Z*(aux00822221(1,1,1)+4*F(5)*temp0072221(1
     &  ,1,1)+8*F(4)*temp0072222(1,1,1)+F(3)*temp822221(1,1,1)+64*temp00
     &  006222(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+48*(temp00006222(2,1,1)*
     &  ZZ(k,1,l,1)+temp00006221(1,1,1)*ZZ(k,2,l,2)))
       temp00822222(1,1,1)=I38Z*(aux00822222(1,1,1)+6*F(4)*temp0072222(2
     &  ,1,1)+F(3)*temp822222(1,1,1)+24*temp00006222(2,2,1)*ZZ(k,1,l,1)+
     &  60*temp00006222(2,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+80*temp00006222
     &  (1,1,1)*ZZ(k,2,l,2)+temp0072222(1,1,1)*(10*r10*(ZZ(k,1,l,2)+ZZ(k
     &  ,2,l,1))+20*r21*ZZ(k,2,l,2)))
       temp00822222(2,1,1)=I38Z*(aux00822222(2,1,1)+4*F(4)*temp0072222(2
     &  ,2,1)+F(3)*temp822222(2,1,1)+8*temp00006222(2,2,2)*ZZ(k,1,l,1)+4
     &  8*temp00006222(2,2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+120*temp00006222
     &  (2,1,1)*ZZ(k,2,l,2)+temp0072222(2,1,1)*(12*r10*(ZZ(k,1,l,2)+ZZ(k
     &  ,2,l,1))+24*r21*ZZ(k,2,l,2)))
       temp00822222(2,2,1)=I38Z*(aux00822222(2,2,1)+2*F(4)*temp0072222(2
     &  ,2,2)+F(3)*temp822222(2,2,1)+28*temp00006222(2,2,2)*(ZZ(k,1,l,2)
     &  +ZZ(k,2,l,1))+168*temp00006222(2,2,1)*ZZ(k,2,l,2)+temp0072222(2,
     &  2,1)*(14*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+28*r21*ZZ(k,2,l,2)))
       temp00822222(2,2,2)=I38Z*(aux00822222(2,2,2)+8*F(5)*temp0072222(2
     &  ,2,2)+F(3)*temp822222(2,2,2)+224*temp00006222(2,2,2)*ZZ(k,2,l,2)
     &  )
       temp00811111(1,1,2)=temp00821111(1,1,1)
       temp00811111(1,2,1)=temp00821111(1,1,1)
       temp00811111(1,2,2)=temp00822111(1,1,1)
       temp00821111(1,1,2)=temp00822111(1,1,1)
       temp00821111(1,2,1)=temp00822111(1,1,1)
       temp00821111(1,2,2)=temp00822211(1,1,1)
       temp00822111(1,1,2)=temp00822211(1,1,1)
       temp00822111(1,2,1)=temp00822211(1,1,1)
       temp00822111(1,2,2)=temp00822221(1,1,1)
       temp00822211(1,1,2)=temp00822221(1,1,1)
       temp00822211(1,2,1)=temp00822221(1,1,1)
       temp00822211(1,2,2)=temp00822222(1,1,1)
       temp00822221(1,1,2)=temp00822222(1,1,1)
       temp00822221(1,2,1)=temp00822222(1,1,1)
       temp00822221(1,2,2)=temp00822222(2,1,1)
       temp00822222(1,1,2)=temp00822222(2,1,1)
       temp00822222(1,2,1)=temp00822222(2,1,1)
       temp00822222(1,2,2)=temp00822222(2,2,1)
       temp00822222(2,1,2)=temp00822222(2,2,1)
       aux9111111(1,1,1)=-(S210111111111(1)*Z(jj,1))-S210211111111(1)*Z(
     &  jj,2)
       aux9211111(1,1,1)=-(S210121111111(1)*Z(jj,1))-S210221111111(1)*Z(
     &  jj,2)
       aux9221111(1,1,1)=-(S210122111111(1)*Z(jj,1))-S210222111111(1)*Z(
     &  jj,2)
       aux9222111(1,1,1)=-(S210122211111(1)*Z(jj,1))-S210222211111(1)*Z(
     &  jj,2)
       aux9222211(1,1,1)=-(S210122221111(1)*Z(jj,1))-S210222221111(1)*Z(
     &  jj,2)
       aux9222221(1,1,1)=-(S210122222111(1)*Z(jj,1))-S210222222111(1)*Z(
     &  jj,2)
       aux9222222(1,1,1)=-(S210122222211(1)*Z(jj,1))-S210222222211(1)*Z(
     &  jj,2)
       aux9222222(2,1,1)=-(S210122222221(1)*Z(jj,1))-S210222222221(1)*Z(
     &  jj,2)
       aux9222222(2,2,1)=-(S210122222222(1)*Z(jj,1))-S210222222222(1)*Z(
     &  jj,2)
       aux9222222(2,2,2)=-(S210122222222(2)*Z(jj,1))-S210222222222(2)*Z(
     &  jj,2)
       temp9111111(1,1,1)=IX*(aux9111111(1,1,1)+18*temp00811111(1,1,1)*Z
     &  (jj,1))
       temp9211111(1,1,1)=IX*(aux9211111(1,1,1)+16*temp00821111(1,1,1)*Z
     &  (jj,1)+2*temp00811111(1,1,1)*Z(jj,2))
       temp9221111(1,1,1)=IX*(aux9221111(1,1,1)+14*temp00822111(1,1,1)*Z
     &  (jj,1)+4*temp00821111(1,1,1)*Z(jj,2))
       temp9222111(1,1,1)=IX*(aux9222111(1,1,1)+12*temp00822211(1,1,1)*Z
     &  (jj,1)+6*temp00822111(1,1,1)*Z(jj,2))
       temp9222211(1,1,1)=IX*(aux9222211(1,1,1)+10*temp00822221(1,1,1)*Z
     &  (jj,1)+8*temp00822211(1,1,1)*Z(jj,2))
       temp9222221(1,1,1)=IX*(aux9222221(1,1,1)+8*temp00822222(1,1,1)*Z(
     &  jj,1)+10*temp00822221(1,1,1)*Z(jj,2))
       temp9222222(1,1,1)=IX*(aux9222222(1,1,1)+6*temp00822222(2,1,1)*Z(
     &  jj,1)+12*temp00822222(1,1,1)*Z(jj,2))
       temp9222222(2,1,1)=IX*(aux9222222(2,1,1)+4*temp00822222(2,2,1)*Z(
     &  jj,1)+14*temp00822222(2,1,1)*Z(jj,2))
       temp9222222(2,2,1)=IX*(aux9222222(2,2,1)+2*temp00822222(2,2,2)*Z(
     &  jj,1)+16*temp00822222(2,2,1)*Z(jj,2))
       temp9222222(2,2,2)=IX*(aux9222222(2,2,2)+18*temp00822222(2,2,2)*Z
     &  (jj,2))
       temp9111111(1,1,2)=temp9211111(1,1,1)
       temp9111111(1,2,1)=temp9211111(1,1,1)
       temp9111111(1,2,2)=temp9221111(1,1,1)
       temp9211111(1,1,2)=temp9221111(1,1,1)
       temp9211111(1,2,1)=temp9221111(1,1,1)
       temp9211111(1,2,2)=temp9222111(1,1,1)
       temp9221111(1,1,2)=temp9222111(1,1,1)
       temp9221111(1,2,1)=temp9222111(1,1,1)
       temp9221111(1,2,2)=temp9222211(1,1,1)
       temp9222111(1,1,2)=temp9222211(1,1,1)
       temp9222111(1,2,1)=temp9222211(1,1,1)
       temp9222111(1,2,2)=temp9222221(1,1,1)
       temp9222211(1,1,2)=temp9222221(1,1,1)
       temp9222211(1,2,1)=temp9222221(1,1,1)
       temp9222211(1,2,2)=temp9222222(1,1,1)
       temp9222221(1,1,2)=temp9222222(1,1,1)
       temp9222221(1,2,1)=temp9222222(1,1,1)
       temp9222221(1,2,2)=temp9222222(2,1,1)
       temp9222222(1,1,2)=temp9222222(2,1,1)
       temp9222222(1,2,1)=temp9222222(2,1,1)
       temp9222222(1,2,2)=temp9222222(2,2,1)
       temp9222222(2,1,2)=temp9222222(2,2,1)
c                Step2
       temp000000001(1)=I22Z*(aux000000001(1)+2*tempC300000000*F(4)+F(3)
     &  *temp0000001(1)-det3*temp0000003(1,k,l))
       temp000000001(2)=I22Z*(aux000000001(2)+tempC300000000*F(5)+F(3)*t
     &  emp0000001(2)-det3*temp0000003(2,k,l))
       temp0000003(1,1,1)=I26Z*(aux0000003(1,1,1)+6*F(4)*temp0000002(1,1
     &  )+F(3)*temp00003(1,1,1)-det3*temp0000511(1,k,l)+24*temp000000001
     &  (1)*ZZ(k,1,l,1))
       temp0000003(2,1,1)=I26Z*(aux0000003(2,1,1)+F(5)*temp0000002(1,1)+
     &  4*F(4)*temp0000002(2,1)+F(3)*temp00003(2,1,1)-det3*temp0000521(1
     &  ,k,l)+8*(temp000000001(2)*ZZ(k,1,l,1)+temp000000001(1)*(ZZ(k,1,l
     &  ,2)+ZZ(k,2,l,1))))
       temp0000003(2,2,1)=I26Z*(aux0000003(2,2,1)+2*F(5)*temp0000002(2,1
     &  )+2*F(4)*temp0000002(2,2)+F(3)*temp00003(2,2,1)-det3*temp0000522
     &  (1,k,l)+8*(temp000000001(2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp000000
     &  001(1)*ZZ(k,2,l,2)))
       temp0000003(2,2,2)=I26Z*(aux0000003(2,2,2)+F(3)*temp00003(2,2,2)-
     &  det3*temp0000522(2,k,l)+24*temp000000001(2)*ZZ(k,2,l,2)+temp0000
     &  002(2,2)*(6*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2)))
       temp0000003(1,1,2)=temp0000003(2,1,1)
       temp0000003(1,2,1)=temp0000003(2,1,1)
       temp0000003(1,2,2)=temp0000003(2,2,1)
       temp0000003(2,1,2)=temp0000003(2,2,1)
       temp0000511(1,1,1)=I30Z*(aux0000511(1,1,1)+10*F(4)*temp000041(1,1
     &  ,1)+F(3)*temp00511(1,1,1)-det3*temp0071111(1,k,l)+80*temp0000003
     &  (1,1,1)*ZZ(k,1,l,1))
       temp0000521(1,1,1)=I30Z*(aux0000521(1,1,1)+F(5)*temp000041(1,1,1)
     &  +8*F(4)*temp000042(1,1,1)+F(3)*temp00521(1,1,1)-det3*temp0072111
     &  (1,k,l)+48*temp0000003(2,1,1)*ZZ(k,1,l,1)+16*temp0000003(1,1,1)*
     &  (ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp0000522(1,1,1)=I30Z*(aux0000522(1,1,1)+2*F(5)*temp000042(1,1,
     &  1)+6*F(4)*temp000042(2,1,1)+F(3)*temp00522(1,1,1)-det3*temp00722
     &  11(1,k,l)+24*(temp0000003(2,2,1)*ZZ(k,1,l,1)+temp0000003(2,1,1)*
     &  (ZZ(k,1,l,2)+ZZ(k,2,l,1)))+8*temp0000003(1,1,1)*ZZ(k,2,l,2))
       temp0000522(2,1,1)=I30Z*(aux0000522(2,1,1)+4*F(4)*temp000042(2,2,
     &  1)+F(3)*temp00522(2,1,1)-det3*temp0072221(1,k,l)+8*temp0000003(2
     &  ,2,2)*ZZ(k,1,l,1)+temp000042(2,1,1)*(6*r10*(ZZ(k,1,l,2)+ZZ(k,2,l
     &  ,1))+12*r21*ZZ(k,2,l,2))+24*(temp0000003(2,2,1)*(ZZ(k,1,l,2)+ZZ(
     &  k,2,l,1))+temp0000003(2,1,1)*ZZ(k,2,l,2)))
       temp0000522(2,2,1)=I30Z*(aux0000522(2,2,1)+4*F(5)*temp000042(2,2,
     &  1)+2*F(4)*temp000042(2,2,2)+F(3)*temp00522(2,2,1)-det3*temp00722
     &  22(1,k,l)+16*temp0000003(2,2,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+48*tem
     &  p0000003(2,2,1)*ZZ(k,2,l,2))
       temp0000522(2,2,2)=I30Z*(aux0000522(2,2,2)+F(3)*temp00522(2,2,2)-
     &  det3*temp0072222(2,k,l)+80*temp0000003(2,2,2)*ZZ(k,2,l,2)+temp00
     &  0042(2,2,2)*(10*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+20*r21*ZZ(k,2,l,2)
     &  ))
       temp0000511(1,1,2)=temp0000521(1,1,1)
       temp0000511(1,2,1)=temp0000521(1,1,1)
       temp0000511(1,2,2)=temp0000522(1,1,1)
       temp0000521(1,1,2)=temp0000522(1,1,1)
       temp0000521(1,2,1)=temp0000522(1,1,1)
       temp0000521(1,2,2)=temp0000522(2,1,1)
       temp0000522(1,1,2)=temp0000522(2,1,1)
       temp0000522(1,2,1)=temp0000522(2,1,1)
       temp0000522(1,2,2)=temp0000522(2,2,1)
       temp0000522(2,1,2)=temp0000522(2,2,1)
       temp0071111(1,1,1)=I34Z*(aux0071111(1,1,1)+14*F(4)*temp006111(1,1
     &  ,1)+F(3)*temp71111(1,1,1)-det3*temp9111111(1,k,l)+168*temp000051
     &  1(1,1,1)*ZZ(k,1,l,1))
       temp0072111(1,1,1)=I34Z*(aux0072111(1,1,1)+F(5)*temp006111(1,1,1)
     &  +12*F(4)*temp006211(1,1,1)+F(3)*temp72111(1,1,1)-det3*temp921111
     &  1(1,k,l)+120*temp0000521(1,1,1)*ZZ(k,1,l,1)+24*temp0000511(1,1,1
     &  )*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp0072211(1,1,1)=I34Z*(aux0072211(1,1,1)+2*F(5)*temp006211(1,1,
     &  1)+10*F(4)*temp006221(1,1,1)+F(3)*temp72211(1,1,1)-det3*temp9221
     &  111(1,k,l)+80*temp0000522(1,1,1)*ZZ(k,1,l,1)+40*temp0000521(1,1,
     &  1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+8*temp0000511(1,1,1)*ZZ(k,2,l,2))
       temp0072221(1,1,1)=I34Z*(aux0072221(1,1,1)+8*F(4)*temp006222(1,1,
     &  1)+F(3)*temp72221(1,1,1)-det3*temp9222111(1,k,l)+48*(temp0000522
     &  (2,1,1)*ZZ(k,1,l,1)+temp0000522(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))
     &  )+24*temp0000521(1,1,1)*ZZ(k,2,l,2)+temp006221(1,1,1)*(6*r10*(ZZ
     &  (k,1,l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2)))
       temp0072222(1,1,1)=I34Z*(aux0072222(1,1,1)+4*F(5)*temp006222(1,1,
     &  1)+6*F(4)*temp006222(2,1,1)+F(3)*temp72222(1,1,1)-det3*temp92222
     &  11(1,k,l)+24*temp0000522(2,2,1)*ZZ(k,1,l,1)+48*(temp0000522(2,1,
     &  1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp0000522(1,1,1)*ZZ(k,2,l,2)))
       temp0072222(2,1,1)=I34Z*(aux0072222(2,1,1)+4*F(4)*temp006222(2,2,
     &  1)+F(3)*temp72222(2,1,1)-det3*temp9222221(1,k,l)+8*temp0000522(2
     &  ,2,2)*ZZ(k,1,l,1)+40*temp0000522(2,2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)
     &  )+80*temp0000522(2,1,1)*ZZ(k,2,l,2)+temp006222(2,1,1)*(10*r10*(Z
     &  Z(k,1,l,2)+ZZ(k,2,l,1))+20*r21*ZZ(k,2,l,2)))
       temp0072222(2,2,1)=I34Z*(aux0072222(2,2,1)+2*F(4)*temp006222(2,2,
     &  2)+F(3)*temp72222(2,2,1)-det3*temp9222222(1,k,l)+24*temp0000522(
     &  2,2,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+120*temp0000522(2,2,1)*ZZ(k,2,l
     &  ,2)+temp006222(2,2,1)*(12*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+24*r21*Z
     &  Z(k,2,l,2)))
       temp0072222(2,2,2)=I34Z*(aux0072222(2,2,2)+F(3)*temp72222(2,2,2)-
     &  det3*temp9222222(2,k,l)+168*temp0000522(2,2,2)*ZZ(k,2,l,2)+temp0
     &  06222(2,2,2)*(14*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+28*r21*ZZ(k,2,l,2
     &  )))
       temp0071111(1,1,2)=temp0072111(1,1,1)
       temp0071111(1,2,1)=temp0072111(1,1,1)
       temp0071111(1,2,2)=temp0072211(1,1,1)
       temp0072111(1,1,2)=temp0072211(1,1,1)
       temp0072111(1,2,1)=temp0072211(1,1,1)
       temp0072111(1,2,2)=temp0072221(1,1,1)
       temp0072211(1,1,2)=temp0072221(1,1,1)
       temp0072211(1,2,1)=temp0072221(1,1,1)
       temp0072211(1,2,2)=temp0072222(1,1,1)
       temp0072221(1,1,2)=temp0072222(1,1,1)
       temp0072221(1,2,1)=temp0072222(1,1,1)
       temp0072221(1,2,2)=temp0072222(2,1,1)
       temp0072222(1,1,2)=temp0072222(2,1,1)
       temp0072222(1,2,1)=temp0072222(2,1,1)
       temp0072222(1,2,2)=temp0072222(2,2,1)
       temp0072222(2,1,2)=temp0072222(2,2,1)
       temp811111(1,1,1)=IX*(aux811111(1,1,1)+det3*temp9111111(1,1,jj)+1
     &  6*temp0071111(1,1,1)*Z(jj,1))
       temp821111(1,1,1)=IX*(aux821111(1,1,1)+det3*temp9211111(1,1,jj)+1
     &  4*temp0072111(1,1,1)*Z(jj,1)+2*temp0071111(1,1,1)*Z(jj,2))
       temp822111(1,1,1)=IX*(aux822111(1,1,1)+det3*temp9221111(1,1,jj)+1
     &  2*temp0072211(1,1,1)*Z(jj,1)+4*temp0072111(1,1,1)*Z(jj,2))
       temp822211(1,1,1)=IX*(aux822211(1,1,1)+det3*temp9222111(1,1,jj)+1
     &  0*temp0072221(1,1,1)*Z(jj,1)+6*temp0072211(1,1,1)*Z(jj,2))
       temp822221(1,1,1)=IX*(aux822221(1,1,1)+det3*temp9222211(1,1,jj)+8
     &  *(temp0072222(1,1,1)*Z(jj,1)+temp0072221(1,1,1)*Z(jj,2)))
       temp822222(1,1,1)=IX*(aux822222(1,1,1)+det3*temp9222221(1,1,jj)+6
     &  *temp0072222(2,1,1)*Z(jj,1)+10*temp0072222(1,1,1)*Z(jj,2))
       temp822222(2,1,1)=IX*(aux822222(2,1,1)+det3*temp9222222(1,1,jj)+4
     &  *temp0072222(2,2,1)*Z(jj,1)+12*temp0072222(2,1,1)*Z(jj,2))
       temp822222(2,2,1)=IX*(aux822222(2,2,1)+det3*temp9222222(2,1,jj)+2
     &  *temp0072222(2,2,2)*Z(jj,1)+14*temp0072222(2,2,1)*Z(jj,2))
       temp822222(2,2,2)=IX*(aux822222(2,2,2)+det3*temp9222222(2,2,jj)+1
     &  6*temp0072222(2,2,2)*Z(jj,2))
       temp811111(1,1,2)=temp821111(1,1,1)
       temp811111(1,2,1)=temp821111(1,1,1)
       temp811111(1,2,2)=temp822111(1,1,1)
       temp821111(1,1,2)=temp822111(1,1,1)
       temp821111(1,2,1)=temp822111(1,1,1)
       temp821111(1,2,2)=temp822211(1,1,1)
       temp822111(1,1,2)=temp822211(1,1,1)
       temp822111(1,2,1)=temp822211(1,1,1)
       temp822111(1,2,2)=temp822221(1,1,1)
       temp822211(1,1,2)=temp822221(1,1,1)
       temp822211(1,2,1)=temp822221(1,1,1)
       temp822211(1,2,2)=temp822222(1,1,1)
       temp822221(1,1,2)=temp822222(1,1,1)
       temp822221(1,2,1)=temp822222(1,1,1)
       temp822221(1,2,2)=temp822222(2,1,1)
       temp822222(1,1,2)=temp822222(2,1,1)
       temp822222(1,2,1)=temp822222(2,1,1)
       temp822222(1,2,2)=temp822222(2,2,1)
       temp822222(2,1,2)=temp822222(2,2,1)
c                Step3
       tempC300000000=I18Z*(auxC300000000+tempC3000000*F(3)-det3*temp000
     &  0002(k,l))
       temp0000002(1,1)=I22Z*(aux0000002(1,1)+4*F(4)*temp0000001(1)+F(3)
     &  *temp00002(1,1)-det3*temp000041(1,k,l)+8*tempC300000000*ZZ(k,1,l
     &  ,1))
       temp0000002(2,1)=I22Z*(aux0000002(2,1)+F(5)*temp0000001(1)+2*F(4)
     &  *temp0000001(2)+F(3)*temp00002(2,1)-det3*temp000042(1,k,l)+4*tem
     &  pC300000000*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp0000002(2,2)=I22Z*(aux0000002(2,2)+2*F(5)*temp0000001(2)+F(3)
     &  *temp00002(2,2)-det3*temp000042(2,k,l)+8*tempC300000000*ZZ(k,2,l
     &  ,2))
       temp0000002(1,2)=temp0000002(2,1)
       temp000041(1,1,1)=I26Z*(aux000041(1,1,1)+8*F(4)*temp00003(1,1,1)+
     &  F(3)*temp0041(1,1,1)-det3*temp006111(1,k,l)+48*temp0000002(1,1)*
     &  ZZ(k,1,l,1))
       temp000042(1,1,1)=I26Z*(aux000042(1,1,1)+F(5)*temp00003(1,1,1)+6*
     &  F(4)*temp00003(2,1,1)+F(3)*temp0042(1,1,1)-det3*temp006211(1,k,l
     &  )+24*temp0000002(2,1)*ZZ(k,1,l,1)+12*temp0000002(1,1)*(ZZ(k,1,l,
     &  2)+ZZ(k,2,l,1)))
       temp000042(2,1,1)=I26Z*(aux000042(2,1,1)+2*F(5)*temp00003(2,1,1)+
     &  4*F(4)*temp00003(2,2,1)+F(3)*temp0042(2,1,1)-det3*temp006221(1,k
     &  ,l)+16*temp0000002(2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+8*(temp0000002
     &  (2,2)*ZZ(k,1,l,1)+temp0000002(1,1)*ZZ(k,2,l,2)))
       temp000042(2,2,1)=I26Z*(aux000042(2,2,1)+2*F(4)*temp00003(2,2,2)+
     &  F(3)*temp0042(2,2,1)-det3*temp006222(1,k,l)+12*temp0000002(2,2)*
     &  (ZZ(k,1,l,2)+ZZ(k,2,l,1))+24*temp0000002(2,1)*ZZ(k,2,l,2)+temp00
     &  003(2,2,1)*(6*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2)))
       temp000042(2,2,2)=I26Z*(aux000042(2,2,2)+4*F(5)*temp00003(2,2,2)+
     &  F(3)*temp0042(2,2,2)-det3*temp006222(2,k,l)+48*temp0000002(2,2)*
     &  ZZ(k,2,l,2))
       temp000041(1,1,2)=temp000042(1,1,1)
       temp000041(1,2,1)=temp000042(1,1,1)
       temp000041(1,2,2)=temp000042(2,1,1)
       temp000042(1,1,2)=temp000042(2,1,1)
       temp000042(1,2,1)=temp000042(2,1,1)
       temp000042(1,2,2)=temp000042(2,2,1)
       temp000042(2,1,2)=temp000042(2,2,1)
       temp006111(1,1,1)=I30Z*(aux006111(1,1,1)+12*F(4)*temp00511(1,1,1)
     &  +F(3)*temp6111(1,1,1)-det3*temp811111(1,k,l)+120*temp000041(1,1,
     &  1)*ZZ(k,1,l,1))
       temp006211(1,1,1)=I30Z*(aux006211(1,1,1)+F(5)*temp00511(1,1,1)+10
     &  *F(4)*temp00521(1,1,1)+F(3)*temp6211(1,1,1)-det3*temp821111(1,k,
     &  l)+80*temp000042(1,1,1)*ZZ(k,1,l,1)+20*temp000041(1,1,1)*(ZZ(k,1
     &  ,l,2)+ZZ(k,2,l,1)))
       temp006221(1,1,1)=I30Z*(aux006221(1,1,1)+2*F(5)*temp00521(1,1,1)+
     &  F(3)*temp6221(1,1,1)-det3*temp822111(1,k,l)+48*temp000042(2,1,1)
     &  *ZZ(k,1,l,1)+32*temp000042(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+8*(F
     &  (4)*temp00522(1,1,1)+temp000041(1,1,1)*ZZ(k,2,l,2)))
       temp006222(1,1,1)=I30Z*(aux006222(1,1,1)+6*F(4)*temp00522(2,1,1)+
     &  F(3)*temp6222(1,1,1)-det3*temp822211(1,k,l)+36*temp000042(2,1,1)
     &  *(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp00522(1,1,1)*(6*r10*(ZZ(k,1,l,2)+
     &  ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2))+24*(temp000042(2,2,1)*ZZ(k,1,l,
     &  1)+temp000042(1,1,1)*ZZ(k,2,l,2)))
       temp006222(2,1,1)=I30Z*(aux006222(2,1,1)+4*F(5)*temp00522(2,1,1)+
     &  4*F(4)*temp00522(2,2,1)+F(3)*temp6222(2,1,1)-det3*temp822221(1,k
     &  ,l)+8*temp000042(2,2,2)*ZZ(k,1,l,1)+32*temp000042(2,2,1)*(ZZ(k,1
     &  ,l,2)+ZZ(k,2,l,1))+48*temp000042(2,1,1)*ZZ(k,2,l,2))
       temp006222(2,2,1)=I30Z*(aux006222(2,2,1)+2*F(4)*temp00522(2,2,2)+
     &  F(3)*temp6222(2,2,1)-det3*temp822222(1,k,l)+20*temp000042(2,2,2)
     &  *(ZZ(k,1,l,2)+ZZ(k,2,l,1))+80*temp000042(2,2,1)*ZZ(k,2,l,2)+temp
     &  00522(2,2,1)*(10*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+20*r21*ZZ(k,2,l,2
     &  )))
       temp006222(2,2,2)=I30Z*(aux006222(2,2,2)+F(3)*temp6222(2,2,2)-det
     &  3*temp822222(2,k,l)+120*temp000042(2,2,2)*ZZ(k,2,l,2)+temp00522(
     &  2,2,2)*(12*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+24*r21*ZZ(k,2,l,2)))
       temp006111(1,1,2)=temp006211(1,1,1)
       temp006111(1,2,1)=temp006211(1,1,1)
       temp006111(1,2,2)=temp006221(1,1,1)
       temp006211(1,1,2)=temp006221(1,1,1)
       temp006211(1,2,1)=temp006221(1,1,1)
       temp006211(1,2,2)=temp006222(1,1,1)
       temp006221(1,1,2)=temp006222(1,1,1)
       temp006221(1,2,1)=temp006222(1,1,1)
       temp006221(1,2,2)=temp006222(2,1,1)
       temp006222(1,1,2)=temp006222(2,1,1)
       temp006222(1,2,1)=temp006222(2,1,1)
       temp006222(1,2,2)=temp006222(2,2,1)
       temp006222(2,1,2)=temp006222(2,2,1)
       temp71111(1,1,1)=IX*(aux71111(1,1,1)+det3*temp811111(1,1,jj)+14*t
     &  emp006111(1,1,1)*Z(jj,1))
       temp72111(1,1,1)=IX*(aux72111(1,1,1)+det3*temp821111(1,1,jj)+12*t
     &  emp006211(1,1,1)*Z(jj,1)+2*temp006111(1,1,1)*Z(jj,2))
       temp72211(1,1,1)=IX*(aux72211(1,1,1)+det3*temp822111(1,1,jj)+10*t
     &  emp006221(1,1,1)*Z(jj,1)+4*temp006211(1,1,1)*Z(jj,2))
       temp72221(1,1,1)=IX*(aux72221(1,1,1)+det3*temp822211(1,1,jj)+8*te
     &  mp006222(1,1,1)*Z(jj,1)+6*temp006221(1,1,1)*Z(jj,2))
       temp72222(1,1,1)=IX*(aux72222(1,1,1)+det3*temp822221(1,1,jj)+6*te
     &  mp006222(2,1,1)*Z(jj,1)+8*temp006222(1,1,1)*Z(jj,2))
       temp72222(2,1,1)=IX*(aux72222(2,1,1)+det3*temp822222(1,1,jj)+4*te
     &  mp006222(2,2,1)*Z(jj,1)+10*temp006222(2,1,1)*Z(jj,2))
       temp72222(2,2,1)=IX*(aux72222(2,2,1)+det3*temp822222(2,1,jj)+2*te
     &  mp006222(2,2,2)*Z(jj,1)+12*temp006222(2,2,1)*Z(jj,2))
       temp72222(2,2,2)=IX*(aux72222(2,2,2)+det3*temp822222(2,2,jj)+14*t
     &  emp006222(2,2,2)*Z(jj,2))
       temp71111(1,1,2)=temp72111(1,1,1)
       temp71111(1,2,1)=temp72111(1,1,1)
       temp71111(1,2,2)=temp72211(1,1,1)
       temp72111(1,1,2)=temp72211(1,1,1)
       temp72111(1,2,1)=temp72211(1,1,1)
       temp72111(1,2,2)=temp72221(1,1,1)
       temp72211(1,1,2)=temp72221(1,1,1)
       temp72211(1,2,1)=temp72221(1,1,1)
       temp72211(1,2,2)=temp72222(1,1,1)
       temp72221(1,1,2)=temp72222(1,1,1)
       temp72221(1,2,1)=temp72222(1,1,1)
       temp72221(1,2,2)=temp72222(2,1,1)
       temp72222(1,1,2)=temp72222(2,1,1)
       temp72222(1,2,1)=temp72222(2,1,1)
       temp72222(1,2,2)=temp72222(2,2,1)
       temp72222(2,1,2)=temp72222(2,2,1)
c                Step4
       temp0000001(1)=I18Z*(aux0000001(1)+2*tempC3000000*F(4)+F(3)*temp0
     &  0001(1)-det3*temp00003(1,k,l))
       temp0000001(2)=I18Z*(aux0000001(2)+tempC3000000*F(5)+F(3)*temp000
     &  01(2)-det3*temp00003(2,k,l))
       temp00003(1,1,1)=I22Z*(aux00003(1,1,1)+6*F(4)*temp00002(1,1)+F(3)
     &  *temp003(1,1,1)-det3*temp00511(1,k,l)+24*temp0000001(1)*ZZ(k,1,l
     &  ,1))
       temp00003(2,1,1)=I22Z*(aux00003(2,1,1)+F(5)*temp00002(1,1)+4*F(4)
     &  *temp00002(2,1)+F(3)*temp003(2,1,1)-det3*temp00521(1,k,l)+8*(tem
     &  p0000001(2)*ZZ(k,1,l,1)+temp0000001(1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))
     &  ))
       temp00003(2,2,1)=I22Z*(aux00003(2,2,1)+2*F(5)*temp00002(2,1)+2*F(
     &  4)*temp00002(2,2)+F(3)*temp003(2,2,1)-det3*temp00522(1,k,l)+8*(t
     &  emp0000001(2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp0000001(1)*ZZ(k,2,l,
     &  2)))
       temp00003(2,2,2)=I22Z*(aux00003(2,2,2)+F(3)*temp003(2,2,2)-det3*t
     &  emp00522(2,k,l)+24*temp0000001(2)*ZZ(k,2,l,2)+temp00002(2,2)*(6*
     &  r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2)))
       temp00003(1,1,2)=temp00003(2,1,1)
       temp00003(1,2,1)=temp00003(2,1,1)
       temp00003(1,2,2)=temp00003(2,2,1)
       temp00003(2,1,2)=temp00003(2,2,1)
       temp00511(1,1,1)=I26Z*(aux00511(1,1,1)+10*F(4)*temp0041(1,1,1)+F(
     &  3)*temp511(1,1,1)-det3*temp71111(1,k,l)+80*temp00003(1,1,1)*ZZ(k
     &  ,1,l,1))
       temp00521(1,1,1)=I26Z*(aux00521(1,1,1)+F(5)*temp0041(1,1,1)+8*F(4
     &  )*temp0042(1,1,1)+F(3)*temp521(1,1,1)-det3*temp72111(1,k,l)+48*t
     &  emp00003(2,1,1)*ZZ(k,1,l,1)+16*temp00003(1,1,1)*(ZZ(k,1,l,2)+ZZ(
     &  k,2,l,1)))
       temp00522(1,1,1)=I26Z*(aux00522(1,1,1)+2*F(5)*temp0042(1,1,1)+6*F
     &  (4)*temp0042(2,1,1)+F(3)*temp522(1,1,1)-det3*temp72211(1,k,l)+24
     &  *(temp00003(2,2,1)*ZZ(k,1,l,1)+temp00003(2,1,1)*(ZZ(k,1,l,2)+ZZ(
     &  k,2,l,1)))+8*temp00003(1,1,1)*ZZ(k,2,l,2))
       temp00522(2,1,1)=I26Z*(aux00522(2,1,1)+4*F(4)*temp0042(2,2,1)+F(3
     &  )*temp522(2,1,1)-det3*temp72221(1,k,l)+8*temp00003(2,2,2)*ZZ(k,1
     &  ,l,1)+temp0042(2,1,1)*(6*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+12*r21*ZZ
     &  (k,2,l,2))+24*(temp00003(2,2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp00
     &  003(2,1,1)*ZZ(k,2,l,2)))
       temp00522(2,2,1)=I26Z*(aux00522(2,2,1)+4*F(5)*temp0042(2,2,1)+2*F
     &  (4)*temp0042(2,2,2)+F(3)*temp522(2,2,1)-det3*temp72222(1,k,l)+16
     &  *temp00003(2,2,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+48*temp00003(2,2,1)*
     &  ZZ(k,2,l,2))
       temp00522(2,2,2)=I26Z*(aux00522(2,2,2)+F(3)*temp522(2,2,2)-det3*t
     &  emp72222(2,k,l)+80*temp00003(2,2,2)*ZZ(k,2,l,2)+temp0042(2,2,2)*
     &  (10*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+20*r21*ZZ(k,2,l,2)))
       temp00511(1,1,2)=temp00521(1,1,1)
       temp00511(1,2,1)=temp00521(1,1,1)
       temp00511(1,2,2)=temp00522(1,1,1)
       temp00521(1,1,2)=temp00522(1,1,1)
       temp00521(1,2,1)=temp00522(1,1,1)
       temp00521(1,2,2)=temp00522(2,1,1)
       temp00522(1,1,2)=temp00522(2,1,1)
       temp00522(1,2,1)=temp00522(2,1,1)
       temp00522(1,2,2)=temp00522(2,2,1)
       temp00522(2,1,2)=temp00522(2,2,1)
       temp6111(1,1,1)=IX*(aux6111(1,1,1)+det3*temp71111(1,1,jj)+12*temp
     &  00511(1,1,1)*Z(jj,1))
       temp6211(1,1,1)=IX*(aux6211(1,1,1)+det3*temp72111(1,1,jj)+10*temp
     &  00521(1,1,1)*Z(jj,1)+2*temp00511(1,1,1)*Z(jj,2))
       temp6221(1,1,1)=IX*(aux6221(1,1,1)+det3*temp72211(1,1,jj)+8*temp0
     &  0522(1,1,1)*Z(jj,1)+4*temp00521(1,1,1)*Z(jj,2))
       temp6222(1,1,1)=IX*(aux6222(1,1,1)+det3*temp72221(1,1,jj)+6*(temp
     &  00522(2,1,1)*Z(jj,1)+temp00522(1,1,1)*Z(jj,2)))
       temp6222(2,1,1)=IX*(aux6222(2,1,1)+det3*temp72222(1,1,jj)+4*temp0
     &  0522(2,2,1)*Z(jj,1)+8*temp00522(2,1,1)*Z(jj,2))
       temp6222(2,2,1)=IX*(aux6222(2,2,1)+det3*temp72222(2,1,jj)+2*temp0
     &  0522(2,2,2)*Z(jj,1)+10*temp00522(2,2,1)*Z(jj,2))
       temp6222(2,2,2)=IX*(aux6222(2,2,2)+det3*temp72222(2,2,jj)+12*temp
     &  00522(2,2,2)*Z(jj,2))
       temp6111(1,1,2)=temp6211(1,1,1)
       temp6111(1,2,1)=temp6211(1,1,1)
       temp6111(1,2,2)=temp6221(1,1,1)
       temp6211(1,1,2)=temp6221(1,1,1)
       temp6211(1,2,1)=temp6221(1,1,1)
       temp6211(1,2,2)=temp6222(1,1,1)
       temp6221(1,1,2)=temp6222(1,1,1)
       temp6221(1,2,1)=temp6222(1,1,1)
       temp6221(1,2,2)=temp6222(2,1,1)
       temp6222(1,1,2)=temp6222(2,1,1)
       temp6222(1,2,1)=temp6222(2,1,1)
       temp6222(1,2,2)=temp6222(2,2,1)
       temp6222(2,1,2)=temp6222(2,2,1)
c                Step5
       tempC3000000=I14Z*(auxC3000000+tempC30000*F(3)-det3*temp00002(k,l
     &  ))
       temp00002(1,1)=I18Z*(aux00002(1,1)+4*F(4)*temp00001(1)+F(3)*temp0
     &  02(1,1)-det3*temp0041(1,k,l)+8*tempC3000000*ZZ(k,1,l,1))
       temp00002(2,1)=I18Z*(aux00002(2,1)+F(5)*temp00001(1)+2*F(4)*temp0
     &  0001(2)+F(3)*temp002(2,1)-det3*temp0042(1,k,l)+4*tempC3000000*(Z
     &  Z(k,1,l,2)+ZZ(k,2,l,1)))
       temp00002(2,2)=I18Z*(aux00002(2,2)+2*F(5)*temp00001(2)+F(3)*temp0
     &  02(2,2)-det3*temp0042(2,k,l)+8*tempC3000000*ZZ(k,2,l,2))
       temp00002(1,2)=temp00002(2,1)
       temp0041(1,1,1)=I22Z*(aux0041(1,1,1)+8*F(4)*temp003(1,1,1)+F(3)*t
     &  emp41(1,1,1)-det3*temp6111(1,k,l)+48*temp00002(1,1)*ZZ(k,1,l,1))
       temp0042(1,1,1)=I22Z*(aux0042(1,1,1)+F(5)*temp003(1,1,1)+6*F(4)*t
     &  emp003(2,1,1)+F(3)*temp42(1,1,1)-det3*temp6211(1,k,l)+24*temp000
     &  02(2,1)*ZZ(k,1,l,1)+12*temp00002(1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp0042(2,1,1)=I22Z*(aux0042(2,1,1)+2*F(5)*temp003(2,1,1)+4*F(4)
     &  *temp003(2,2,1)+F(3)*temp42(2,1,1)-det3*temp6221(1,k,l)+16*temp0
     &  0002(2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+8*(temp00002(2,2)*ZZ(k,1,l,1
     &  )+temp00002(1,1)*ZZ(k,2,l,2)))
       temp0042(2,2,1)=I22Z*(aux0042(2,2,1)+2*F(4)*temp003(2,2,2)+F(3)*t
     &  emp42(2,2,1)-det3*temp6222(1,k,l)+12*temp00002(2,2)*(ZZ(k,1,l,2)
     &  +ZZ(k,2,l,1))+24*temp00002(2,1)*ZZ(k,2,l,2)+temp003(2,2,1)*(6*r1
     &  0*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2)))
       temp0042(2,2,2)=I22Z*(aux0042(2,2,2)+4*F(5)*temp003(2,2,2)+F(3)*t
     &  emp42(2,2,2)-det3*temp6222(2,k,l)+48*temp00002(2,2)*ZZ(k,2,l,2))
       temp0041(1,1,2)=temp0042(1,1,1)
       temp0041(1,2,1)=temp0042(1,1,1)
       temp0041(1,2,2)=temp0042(2,1,1)
       temp0042(1,1,2)=temp0042(2,1,1)
       temp0042(1,2,1)=temp0042(2,1,1)
       temp0042(1,2,2)=temp0042(2,2,1)
       temp0042(2,1,2)=temp0042(2,2,1)
       temp511(1,1,1)=IX*(aux511(1,1,1)+det3*temp6111(1,1,jj)+10*temp004
     &  1(1,1,1)*Z(jj,1))
       temp521(1,1,1)=IX*(aux521(1,1,1)+det3*temp6211(1,1,jj)+8*temp0042
     &  (1,1,1)*Z(jj,1)+2*temp0041(1,1,1)*Z(jj,2))
       temp522(1,1,1)=IX*(aux522(1,1,1)+det3*temp6221(1,1,jj)+6*temp0042
     &  (2,1,1)*Z(jj,1)+4*temp0042(1,1,1)*Z(jj,2))
       temp522(2,1,1)=IX*(aux522(2,1,1)+det3*temp6222(1,1,jj)+4*temp0042
     &  (2,2,1)*Z(jj,1)+6*temp0042(2,1,1)*Z(jj,2))
       temp522(2,2,1)=IX*(aux522(2,2,1)+det3*temp6222(2,1,jj)+2*temp0042
     &  (2,2,2)*Z(jj,1)+8*temp0042(2,2,1)*Z(jj,2))
       temp522(2,2,2)=IX*(aux522(2,2,2)+det3*temp6222(2,2,jj)+10*temp004
     &  2(2,2,2)*Z(jj,2))
       temp511(1,1,2)=temp521(1,1,1)
       temp511(1,2,1)=temp521(1,1,1)
       temp511(1,2,2)=temp522(1,1,1)
       temp521(1,1,2)=temp522(1,1,1)
       temp521(1,2,1)=temp522(1,1,1)
       temp521(1,2,2)=temp522(2,1,1)
       temp522(1,1,2)=temp522(2,1,1)
       temp522(1,2,1)=temp522(2,1,1)
       temp522(1,2,2)=temp522(2,2,1)
       temp522(2,1,2)=temp522(2,2,1)
c                Step6
       temp00001(1)=I14Z*(aux00001(1)+2*tempC30000*F(4)+F(3)*temp001(1)-
     &  det3*temp003(1,k,l))
       temp00001(2)=I14Z*(aux00001(2)+tempC30000*F(5)+F(3)*temp001(2)-de
     &  t3*temp003(2,k,l))
       temp003(1,1,1)=I18Z*(aux003(1,1,1)+6*F(4)*temp002(1,1)+F(3)*temp3
     &  (1,1,1)-det3*temp511(1,k,l)+24*temp00001(1)*ZZ(k,1,l,1))
       temp003(2,1,1)=I18Z*(aux003(2,1,1)+F(5)*temp002(1,1)+4*F(4)*temp0
     &  02(2,1)+F(3)*temp3(2,1,1)-det3*temp521(1,k,l)+8*(temp00001(2)*ZZ
     &  (k,1,l,1)+temp00001(1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))))
       temp003(2,2,1)=I18Z*(aux003(2,2,1)+2*F(5)*temp002(2,1)+2*F(4)*tem
     &  p002(2,2)+F(3)*temp3(2,2,1)-det3*temp522(1,k,l)+8*(temp00001(2)*
     &  (ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp00001(1)*ZZ(k,2,l,2)))
       temp003(2,2,2)=I18Z*(aux003(2,2,2)+F(3)*temp3(2,2,2)-det3*temp522
     &  (2,k,l)+24*temp00001(2)*ZZ(k,2,l,2)+temp002(2,2)*(6*r10*(ZZ(k,1,
     &  l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2)))
       temp003(1,1,2)=temp003(2,1,1)
       temp003(1,2,1)=temp003(2,1,1)
       temp003(1,2,2)=temp003(2,2,1)
       temp003(2,1,2)=temp003(2,2,1)

       temp41(1,1,1)=IX*(aux41(1,1,1)+det3*temp511(1,1,jj)+8*temp003(1,1
     &  ,1)*Z(jj,1))
       temp42(1,1,1)=IX*(aux42(1,1,1)+det3*temp521(1,1,jj)+6*temp003(2,1
     &  ,1)*Z(jj,1)+2*temp003(1,1,1)*Z(jj,2))
       temp42(2,1,1)=IX*(aux42(2,1,1)+det3*temp522(1,1,jj)+4*(temp003(2,
     &  2,1)*Z(jj,1)+temp003(2,1,1)*Z(jj,2)))
       temp42(2,2,1)=IX*(aux42(2,2,1)+det3*temp522(2,1,jj)+2*temp003(2,2
     &  ,2)*Z(jj,1)+6*temp003(2,2,1)*Z(jj,2))
       temp42(2,2,2)=IX*(aux42(2,2,2)+det3*temp522(2,2,jj)+8*temp003(2,2
     &  ,2)*Z(jj,2))
       temp41(1,1,2)=temp42(1,1,1)
       temp41(1,2,1)=temp42(1,1,1)
       temp41(1,2,2)=temp42(2,1,1)
       temp42(1,1,2)=temp42(2,1,1)
       temp42(1,2,1)=temp42(2,1,1)
       temp42(1,2,2)=temp42(2,2,1)
       temp42(2,1,2)=temp42(2,2,1)
c                Step7
       tempC30000=I10Z*(auxC30000+tempC300*F(3)-det3*temp002(k,l))
       temp002(1,1)=I14Z*(aux002(1,1)+4*F(4)*temp001(1)+F(3)*temp2(1,1)-
     &  det3*temp41(1,k,l)+8*tempC30000*ZZ(k,1,l,1))
       temp002(2,1)=I14Z*(aux002(2,1)+F(5)*temp001(1)+2*F(4)*temp001(2)+
     &  F(3)*temp2(2,1)-det3*temp42(1,k,l)+4*tempC30000*(ZZ(k,1,l,2)+ZZ(
     &  k,2,l,1)))
       temp002(2,2)=I14Z*(aux002(2,2)+2*F(5)*temp001(2)+F(3)*temp2(2,2)-
     &  det3*temp42(2,k,l)+8*tempC30000*ZZ(k,2,l,2))
       temp002(1,2)=temp002(2,1)
       temp3(1,1,1)=IX*(aux3(1,1,1)+det3*temp41(1,1,jj)+6*temp002(1,1)*Z
     &  (jj,1))
       temp3(2,1,1)=IX*(aux3(2,1,1)+det3*temp42(1,1,jj)+4*temp002(2,1)*Z
     &  (jj,1)+2*temp002(1,1)*Z(jj,2))
       temp3(2,2,1)=IX*(aux3(2,2,1)+det3*temp42(2,1,jj)+2*temp002(2,2)*Z
     &  (jj,1)+4*temp002(2,1)*Z(jj,2))
       temp3(2,2,2)=IX*(aux3(2,2,2)+det3*temp42(2,2,jj)+6*temp002(2,2)*Z
     &  (jj,2))
       temp3(1,1,2)=temp3(2,1,1)
       temp3(1,2,1)=temp3(2,1,1)
       temp3(1,2,2)=temp3(2,2,1)
       temp3(2,1,2)=temp3(2,2,1)
c                Step8
       temp001(1)=I10Z*(aux001(1)+2*tempC300*F(4)+F(3)*temp1(1)-det3*tem
     &  p3(1,k,l))
       temp001(2)=I10Z*(aux001(2)+tempC300*F(5)+F(3)*temp1(2)-det3*temp3
     &  (2,k,l))
       temp2(1,1)=IX*(aux2(1,1)+det3*temp3(1,1,jj)+4*temp001(1)*Z(jj,1))
       temp2(2,1)=IX*(aux2(2,1)+det3*temp3(2,1,jj)+2*(temp001(2)*Z(jj,1)
     &  +temp001(1)*Z(jj,2)))
       temp2(2,2)=IX*(aux2(2,2)+det3*temp3(2,2,jj)+4*temp001(2)*Z(jj,2))
       temp2(1,2)=temp2(2,1)
c                Step9
       tempC300=I6Z*(auxC300+tempC30*F(3)-det3*temp2(k,l))
       temp1(1)=IX*(aux1(1)+det3*temp2(1,jj)+2*tempC300*Z(jj,1))
       temp1(2)=IX*(aux1(2)+det3*temp2(2,jj)+2*tempC300*Z(jj,2))
c                Step10
       tempC30=IX*(auxC30+det3*temp1(jj))

         ac=3
         accuracyCR(1,0,ac) = abs(C30     /tempC30      -1d0) 
         accuracyCR(1,1,ac) = abs(Cij(1,1)/temp1(1)      -1d0) 
         accuracyCR(2,1,ac) = abs(Cij(2,1)/temp1(2)      -1d0)
         accuracyCR(1,2,ac) = abs(Cij(1,2)/temp2(1,1)    -1d0)
         accuracyCR(2,2,ac) = abs(Cij(2,2)/temp2(2,2)    -1d0)
         accuracyCR(3,2,ac) = abs(Cij(3,2)/temp2(2,1)    -1d0)
         accuracyCR(4,2,ac) = abs(Cij(4,2)/tempC300      -1d0)
         accuracyCR(1,3,ac) = abs(Cij(1,3)/temp3(1,1,1)  -1d0)
         accuracyCR(2,3,ac) = abs(Cij(2,3)/temp3(2,2,2)  -1d0)
         accuracyCR(3,3,ac) = abs(Cij(3,3)/temp3(2,1,1)  -1d0)
         accuracyCR(4,3,ac) = abs(Cij(4,3)/temp3(2,2,1)  -1d0)
         accuracyCR(5,3,ac) = abs(Cij(5,3)/temp001(1)    -1d0)
         accuracyCR(6,3,ac) = abs(Cij(6,3)/temp001(2)    -1d0)
         accuracyCR(1,4,ac) = abs(Cij(1,4)/temp41(1,1,1) -1d0)
         accuracyCR(2,4,ac) = abs(Cij(2,4)/temp42(2,2,2) -1d0)
         accuracyCR(3,4,ac) = abs(Cij(3,4)/temp42(1,1,1) -1d0)
         accuracyCR(4,4,ac) = abs(Cij(4,4)/temp42(2,1,1) -1d0)
         accuracyCR(5,4,ac) = abs(Cij(5,4)/temp42(2,2,1) -1d0)
         accuracyCR(6,4,ac) = abs(Cij(6,4)/temp002(1,1)  -1d0)
         accuracyCR(7,4,ac) = abs(Cij(7,4)/temp002(2,2)  -1d0)
         accuracyCR(8,4,ac) = abs(Cij(8,4)/temp002(2,1)  -1d0)
         accuracyCR(9,4,ac) = abs(Cij(9,4)/tempC30000    -1d0)


      DO I1=0,4
           accuracyC(i1,ac)=accuracyCR(1,i1,ac)
            DO I2=2,INDEX(I1)  
       if(accuracyCR(i2,i1,ac).gt.accuracyC(i1,ac)) then
          accuracyC(i1,ac)=accuracyCR(i2,i1,ac)
       endif
          enddo
        enddo

c          if(accuracyC(4,ac).lt.1d-16) goto 500

           C30=tempC30
           Cij(1,1)=temp1(1)
           Cij(2,1)=temp1(2)
           Cij(1,2)=temp2(1,1)
           Cij(2,2)=temp2(2,2)
           Cij(3,2)=temp2(2,1)
           Cij(4,2)=tempC300
           Cij(1,3)=temp3(1,1,1)
           Cij(2,3)=temp3(2,2,2)
           Cij(3,3)=temp3(2,1,1)
           Cij(4,3)=temp3(2,2,1)
           Cij(5,3)=temp001(1)
           Cij(6,3)=temp001(2)
           Cij(1,4)=temp41(1,1,1)
           Cij(2,4)=temp42(2,2,2)
           Cij(3,4)=temp42(1,1,1)
           Cij(4,4)=temp42(2,1,1)
           Cij(5,4)=temp42(2,2,1)
           Cij(6,4)=temp002(1,1)
           Cij(7,4)=temp002(2,2)
           Cij(8,4)=temp002(2,1)
           Cij(9,4)=tempC30000


       if(order.eq.9) goto 500
c       Print*, "here WRONG AFTER order 9"
c                Iteration10
c                Step1
       S200000000001(1)=-2*B23(5,8)
       S200000000001(2)=2*B23(5,9)
       S2h00000000001(1)=B13(6,10)-B23(6,10)
       S2h00000000001(2)=B12(6,10)-B13(6,10)
       S200000000311(1)=-2*B23(4,6)
       S200000000311(2)=2*B23(4,7)
       S200000000312(1)=2*B23(4,7)
       S200000000312(2)=-2*B23(4,8)
       S200000000321(1)=2*B23(4,7)
       S200000000321(2)=-2*B23(4,8)
       S200000000322(1)=-2*B23(4,8)
       S200000000322(2)=2*B23(4,9)
       S2h00000000311(1)=B13(5,10)-B23(5,8)
       S2h00000000311(2)=B13(5,10)+B23(5,9)
       S2h00000000312(1)=B13(5,10)+B23(5,9)
       S2h00000000312(2)=B13(5,10)-B23(5,10)
       S2h00000000321(1)=B12(5,10)-B13(5,10)
       S2h00000000321(2)=-B13(5,10)
       S2h00000000322(1)=-B13(5,10)
       S2h00000000322(2)=-B13(5,10)
       aux00000000001(1)=-(F(1)*S2h0000000021(1))-F(2)*S2h0000000022(1)+
     &  S2h00000000311(k)*Z(1,l)+S2h00000000321(k)*Z(2,l)+(Inv1001+S20000
     &  0000001(1)-S2h00000000311(1)-S2h00000000322(1))*Z(k,l)-2*S2h0000
     &  0000001(1)*ZZ(k,1,l,1)-2*S2h00000000001(2)*ZZ(k,1,l,2)
       aux00000000001(2)=-(F(1)*S2h0000000021(2))-F(2)*S2h0000000022(2)+
     &  S2h00000000312(k)*Z(1,l)+S2h00000000322(k)*Z(2,l)+(Inv1002+S20000
     &  0000001(2)-S2h00000000311(2)-S2h00000000322(2))*Z(k,l)-2*S2h0000
     &  0000001(1)*ZZ(k,2,l,1)-2*S2h00000000001(2)*ZZ(k,2,l,2)
       temp00000000001(1)=I26Z*(aux00000000001(1)+2*tempC30000000000*F(4
     &  )+F(3)*temp000000001(1))
       temp00000000001(2)=I26Z*(aux00000000001(2)+tempC30000000000*F(5)+
     &  F(3)*temp000000001(2))
       S200000051111(1)=-2*B23(3,4)
       S200000051111(2)=2*B23(3,5)
       S200000051211(1)=2*B23(3,5)
       S200000051211(2)=-2*B23(3,6)
       S200000051221(1)=-2*B23(3,6)
       S200000051221(2)=2*B23(3,7)
       S200000051222(1)=2*B23(3,7)
       S200000051222(2)=-2*B23(3,8)
       S200000052111(1)=2*B23(3,5)
       S200000052111(2)=-2*B23(3,6)
       S200000052211(1)=-2*B23(3,6)
       S200000052211(2)=2*B23(3,7)
       S200000052221(1)=2*B23(3,7)
       S200000052221(2)=-2*B23(3,8)
       S200000052222(1)=-2*B23(3,8)
       S200000052222(2)=2*B23(3,9)
       S2h00000051111(1)=B13(4,10)-B23(4,6)
       S2h00000051111(2)=B13(4,10)+B23(4,7)
       S2h00000051211(1)=B13(4,10)+B23(4,7)
       S2h00000051211(2)=B13(4,10)-B23(4,8)
       S2h00000051221(1)=B13(4,10)-B23(4,8)
       S2h00000051221(2)=B13(4,10)+B23(4,9)
       S2h00000051222(1)=B13(4,10)+B23(4,9)
       S2h00000051222(2)=B13(4,10)-B23(4,10)
       S2h00000052111(1)=B12(4,10)-B13(4,10)
       S2h00000052111(2)=-B13(4,10)
       S2h00000052211(1)=-B13(4,10)
       S2h00000052211(2)=-B13(4,10)
       S2h00000052221(1)=-B13(4,10)
       S2h00000052221(2)=-B13(4,10)
       S2h00000052222(1)=-B13(4,10)
       S2h00000052222(2)=-B13(4,10)
       aux000000003(1,1,1)=-(F(1)*S2h0000004111(1))-F(2)*S2h0000004211(1
     &  )+S2h00000051111(k)*Z(1,l)+S2h00000052111(k)*Z(2,l)+(Inv80111+S2
     &  00000000311(1)-S2h00000051111(1)-S2h00000052211(1))*Z(k,l)-6*(S2
     &  h00000000311(1)*ZZ(k,1,l,1)+S2h00000000321(1)*ZZ(k,1,l,2))
       aux000000003(2,1,1)=-(F(1)*S2h0000004121(1))-F(2)*S2h0000004221(1
     &  )+S2h00000051211(k)*Z(1,l)+S2h00000052211(k)*Z(2,l)+(Inv80211+S2
     &  00000000321(1)-S2h00000051211(1)-S2h00000052221(1))*Z(k,l)-4*(S2
     &  h00000000312(1)*ZZ(k,1,l,1)+S2h00000000322(1)*ZZ(k,1,l,2))-2*S2h
     &  00000000311(1)*ZZ(k,2,l,1)-2*S2h00000000321(1)*ZZ(k,2,l,2)
       aux000000003(2,2,1)=-(F(1)*S2h0000004122(1))-F(2)*S2h0000004222(1
     &  )+S2h00000051221(k)*Z(1,l)+S2h00000052221(k)*Z(2,l)+(Inv80221+S2
     &  00000000322(1)-S2h00000051221(1)-S2h00000052222(1))*Z(k,l)-2*(S2
     &  h00000000312(2)*ZZ(k,1,l,1)+S2h00000000322(2)*ZZ(k,1,l,2))-4*S2h
     &  00000000312(1)*ZZ(k,2,l,1)-4*S2h00000000322(1)*ZZ(k,2,l,2)
       aux000000003(2,2,2)=-(F(1)*S2h0000004122(2))-F(2)*S2h0000004222(2
     &  )+S2h00000051222(k)*Z(1,l)+S2h00000052222(k)*Z(2,l)+(Inv80222+S2
     &  00000000322(2)-S2h00000051222(1)-S2h00000052222(2))*Z(k,l)-6*(S2
     &  h00000000312(2)*ZZ(k,2,l,1)+S2h00000000322(2)*ZZ(k,2,l,2))
       temp000000003(1,1,1)=I30Z*(aux000000003(1,1,1)+6*F(4)*temp0000000
     &  02(1,1)+F(3)*temp0000003(1,1,1)+24*temp00000000001(1)*ZZ(k,1,l,1
     &  ))
       temp000000003(2,1,1)=I30Z*(aux000000003(2,1,1)+F(5)*temp000000002
     &  (1,1)+4*F(4)*temp000000002(2,1)+F(3)*temp0000003(2,1,1)+8*(temp0
     &  0000000001(2)*ZZ(k,1,l,1)+temp00000000001(1)*(ZZ(k,1,l,2)+ZZ(k,2
     &  ,l,1))))
       temp000000003(2,2,1)=I30Z*(aux000000003(2,2,1)+2*F(5)*temp0000000
     &  02(2,1)+2*F(4)*temp000000002(2,2)+F(3)*temp0000003(2,2,1)+8*(tem
     &  p00000000001(2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp00000000001(1)*ZZ(
     &  k,2,l,2)))
       temp000000003(2,2,2)=I30Z*(aux000000003(2,2,2)+F(3)*temp0000003(2
     &  ,2,2)+24*temp00000000001(2)*ZZ(k,2,l,2)+temp000000002(2,2)*(6*r1
     &  0*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2)))
       temp000000003(1,1,2)=temp000000003(2,1,1)
       temp000000003(1,2,1)=temp000000003(2,1,1)
       temp000000003(1,2,2)=temp000000003(2,2,1)
       temp000000003(2,1,2)=temp000000003(2,2,1)
       S200007111111(1)=-2*B23(2,2)
       S200007111111(2)=2*B23(2,3)
       S200007121111(1)=2*B23(2,3)
       S200007121111(2)=-2*B23(2,4)
       S200007122111(1)=-2*B23(2,4)
       S200007122111(2)=2*B23(2,5)
       S200007122211(1)=2*B23(2,5)
       S200007122211(2)=-2*B23(2,6)
       S200007122221(1)=-2*B23(2,6)
       S200007122221(2)=2*B23(2,7)
       S200007122222(1)=2*B23(2,7)
       S200007122222(2)=-2*B23(2,8)
       S200007211111(1)=2*B23(2,3)
       S200007211111(2)=-2*B23(2,4)
       S200007221111(1)=-2*B23(2,4)
       S200007221111(2)=2*B23(2,5)
       S200007222111(1)=2*B23(2,5)
       S200007222111(2)=-2*B23(2,6)
       S200007222211(1)=-2*B23(2,6)
       S200007222211(2)=2*B23(2,7)
       S200007222221(1)=2*B23(2,7)
       S200007222221(2)=-2*B23(2,8)
       S200007222222(1)=-2*B23(2,8)
       S200007222222(2)=2*B23(2,9)
       S2h00007111111(1)=B13(3,10)-B23(3,4)
       S2h00007111111(2)=B13(3,10)+B23(3,5)
       S2h00007121111(1)=B13(3,10)+B23(3,5)
       S2h00007121111(2)=B13(3,10)-B23(3,6)
       S2h00007122111(1)=B13(3,10)-B23(3,6)
       S2h00007122111(2)=B13(3,10)+B23(3,7)
       S2h00007122211(1)=B13(3,10)+B23(3,7)
       S2h00007122211(2)=B13(3,10)-B23(3,8)
       S2h00007122221(1)=B13(3,10)-B23(3,8)
       S2h00007122221(2)=B13(3,10)+B23(3,9)
       S2h00007122222(1)=B13(3,10)+B23(3,9)
       S2h00007122222(2)=B13(3,10)-B23(3,10)
       S2h00007211111(1)=B12(3,10)-B13(3,10)
       S2h00007211111(2)=-B13(3,10)
       S2h00007221111(1)=-B13(3,10)
       S2h00007221111(2)=-B13(3,10)
       S2h00007222111(1)=-B13(3,10)
       S2h00007222111(2)=-B13(3,10)
       S2h00007222211(1)=-B13(3,10)
       S2h00007222211(2)=-B13(3,10)
       S2h00007222221(1)=-B13(3,10)
       S2h00007222221(2)=-B13(3,10)
       S2h00007222222(1)=-B13(3,10)
       S2h00007222222(2)=-B13(3,10)
       aux000000511(1,1,1)=-(F(1)*S2h0000611111(1))-F(2)*S2h0000621111(1
     &  )+S2h00007111111(k)*Z(1,l)+S2h00007211111(k)*Z(2,l)+(Inv6011111+
     &  S200000051111(1)-S2h00007111111(1)-S2h00007221111(1))*Z(k,l)-10*
     &  (S2h00000051111(1)*ZZ(k,1,l,1)+S2h00000052111(1)*ZZ(k,1,l,2))
       aux000000521(1,1,1)=-(F(1)*S2h0000612111(1))-F(2)*S2h0000622111(1
     &  )+S2h00007121111(k)*Z(1,l)+S2h00007221111(k)*Z(2,l)+(Inv6021111+
     &  S200000052111(1)-S2h00007121111(1)-S2h00007222111(1))*Z(k,l)-8*(
     &  S2h00000051211(1)*ZZ(k,1,l,1)+S2h00000052211(1)*ZZ(k,1,l,2))-2*S
     &  2h00000051111(1)*ZZ(k,2,l,1)-2*S2h00000052111(1)*ZZ(k,2,l,2)
       aux000000522(1,1,1)=-(F(1)*S2h0000612211(1))-F(2)*S2h0000622211(1
     &  )+S2h00007122111(k)*Z(1,l)+S2h00007222111(k)*Z(2,l)+(Inv6022111+
     &  S200000052211(1)-S2h00007122111(1)-S2h00007222211(1))*Z(k,l)-6*(
     &  S2h00000051221(1)*ZZ(k,1,l,1)+S2h00000052221(1)*ZZ(k,1,l,2))-4*S
     &  2h00000051211(1)*ZZ(k,2,l,1)-4*S2h00000052211(1)*ZZ(k,2,l,2)
       aux000000522(2,1,1)=-(F(1)*S2h0000612221(1))-F(2)*S2h0000622221(1
     &  )+S2h00007122211(k)*Z(1,l)+S2h00007222211(k)*Z(2,l)+(Inv6022211+
     &  S200000052221(1)-S2h00007122211(1)-S2h00007222221(1))*Z(k,l)-4*(
     &  S2h00000051222(1)*ZZ(k,1,l,1)+S2h00000052222(1)*ZZ(k,1,l,2))-6*S
     &  2h00000051221(1)*ZZ(k,2,l,1)-6*S2h00000052221(1)*ZZ(k,2,l,2)
       aux000000522(2,2,1)=-(F(1)*S2h0000612222(1))-F(2)*S2h0000622222(1
     &  )+S2h00007122221(k)*Z(1,l)+S2h00007222221(k)*Z(2,l)+(Inv6022221+
     &  S200000052222(1)-S2h00007122221(1)-S2h00007222222(1))*Z(k,l)-2*(
     &  S2h00000051222(2)*ZZ(k,1,l,1)+S2h00000052222(2)*ZZ(k,1,l,2))-8*S
     &  2h00000051222(1)*ZZ(k,2,l,1)-8*S2h00000052222(1)*ZZ(k,2,l,2)
       aux000000522(2,2,2)=-(F(1)*S2h0000612222(2))-F(2)*S2h0000622222(2
     &  )+S2h00007122222(k)*Z(1,l)+S2h00007222222(k)*Z(2,l)+(Inv6022222+
     &  S200000052222(2)-S2h00007122222(1)-S2h00007222222(2))*Z(k,l)-10*
     &  (S2h00000051222(2)*ZZ(k,2,l,1)+S2h00000052222(2)*ZZ(k,2,l,2))
       temp000000511(1,1,1)=I34Z*(aux000000511(1,1,1)+10*F(4)*temp000000
     &  41(1,1,1)+F(3)*temp0000511(1,1,1)+80*temp000000003(1,1,1)*ZZ(k,1
     &  ,l,1))
       temp000000521(1,1,1)=I34Z*(aux000000521(1,1,1)+F(5)*temp00000041(
     &  1,1,1)+8*F(4)*temp00000042(1,1,1)+F(3)*temp0000521(1,1,1)+48*tem
     &  p000000003(2,1,1)*ZZ(k,1,l,1)+16*temp000000003(1,1,1)*(ZZ(k,1,l,
     &  2)+ZZ(k,2,l,1)))
       temp000000522(1,1,1)=I34Z*(aux000000522(1,1,1)+2*F(5)*temp0000004
     &  2(1,1,1)+6*F(4)*temp00000042(2,1,1)+F(3)*temp0000522(1,1,1)+24*(
     &  temp000000003(2,2,1)*ZZ(k,1,l,1)+temp000000003(2,1,1)*(ZZ(k,1,l,
     &  2)+ZZ(k,2,l,1)))+8*temp000000003(1,1,1)*ZZ(k,2,l,2))
       temp000000522(2,1,1)=I34Z*(aux000000522(2,1,1)+4*F(4)*temp0000004
     &  2(2,2,1)+F(3)*temp0000522(2,1,1)+8*temp000000003(2,2,2)*ZZ(k,1,l
     &  ,1)+temp00000042(2,1,1)*(6*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+12*r21*
     &  ZZ(k,2,l,2))+24*(temp000000003(2,2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+
     &  temp000000003(2,1,1)*ZZ(k,2,l,2)))
       temp000000522(2,2,1)=I34Z*(aux000000522(2,2,1)+4*F(5)*temp0000004
     &  2(2,2,1)+2*F(4)*temp00000042(2,2,2)+F(3)*temp0000522(2,2,1)+16*t
     &  emp000000003(2,2,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+48*temp000000003(2
     &  ,2,1)*ZZ(k,2,l,2))
       temp000000522(2,2,2)=I34Z*(aux000000522(2,2,2)+F(3)*temp0000522(2
     &  ,2,2)+80*temp000000003(2,2,2)*ZZ(k,2,l,2)+temp00000042(2,2,2)*(1
     &  0*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+20*r21*ZZ(k,2,l,2)))
       temp000000511(1,1,2)=temp000000521(1,1,1)
       temp000000511(1,2,1)=temp000000521(1,1,1)
       temp000000511(1,2,2)=temp000000522(1,1,1)
       temp000000521(1,1,2)=temp000000522(1,1,1)
       temp000000521(1,2,1)=temp000000522(1,1,1)
       temp000000521(1,2,2)=temp000000522(2,1,1)
       temp000000522(1,1,2)=temp000000522(2,1,1)
       temp000000522(1,2,1)=temp000000522(2,1,1)
       temp000000522(1,2,2)=temp000000522(2,2,1)
       temp000000522(2,1,2)=temp000000522(2,2,1)
       S2h00911111111(1)=B13(2,10)-B23(2,2)
       S2h00911111111(2)=B13(2,10)+B23(2,3)
       S2h00912111111(1)=B13(2,10)+B23(2,3)
       S2h00912111111(2)=B13(2,10)-B23(2,4)
       S2h00912211111(1)=B13(2,10)-B23(2,4)
       S2h00912211111(2)=B13(2,10)+B23(2,5)
       S2h00912221111(1)=B13(2,10)+B23(2,5)
       S2h00912221111(2)=B13(2,10)-B23(2,6)
       S2h00912222111(1)=B13(2,10)-B23(2,6)
       S2h00912222111(2)=B13(2,10)+B23(2,7)
       S2h00912222211(1)=B13(2,10)+B23(2,7)
       S2h00912222211(2)=B13(2,10)-B23(2,8)
       S2h00912222221(1)=B13(2,10)-B23(2,8)
       S2h00912222221(2)=B13(2,10)+B23(2,9)
       S2h00912222222(1)=B13(2,10)+B23(2,9)
       S2h00912222222(2)=B13(2,10)-B23(2,10)
       S2h00921111111(1)=B12(2,10)-B13(2,10)
       S2h00921111111(2)=-B13(2,10)
       S2h00922111111(1)=-B13(2,10)
       S2h00922111111(2)=-B13(2,10)
       S2h00922211111(1)=-B13(2,10)
       S2h00922211111(2)=-B13(2,10)
       S2h00922221111(1)=-B13(2,10)
       S2h00922221111(2)=-B13(2,10)
       S2h00922222111(1)=-B13(2,10)
       S2h00922222111(2)=-B13(2,10)
       S2h00922222211(1)=-B13(2,10)
       S2h00922222211(2)=-B13(2,10)
       S2h00922222221(1)=-B13(2,10)
       S2h00922222221(2)=-B13(2,10)
       S2h00922222222(1)=-B13(2,10)
       S2h00922222222(2)=-B13(2,10)
       aux000071111(1,1,1)=-(F(1)*S2h0081111111(1))-F(2)*S2h0082111111(1
     &  )+S2h00911111111(k)*Z(1,l)+S2h00921111111(k)*Z(2,l)+(Inv401111
     &  111+S200007111111(1)-S2h00911111111(1)-S2h00922111111(1))*Z(k,l)
     &  -14*(S2h00007111111(1)*ZZ(k,1,l,1)+S2h00007211111(1)*ZZ(k,1,l,2)
     &  )
       aux000072111(1,1,1)=-(F(1)*S2h0081211111(1))-F(2)*S2h0082211111(1
     &  )+S2h00912111111(k)*Z(1,l)+S2h00922111111(k)*Z(2,l)+(Inv402111
     &  111+S200007211111(1)-S2h00912111111(1)-S2h00922211111(1))*Z(k,l)
     &  -12*(S2h00007121111(1)*ZZ(k,1,l,1)+S2h00007221111(1)*ZZ(k,1,l,2)
     &  )-2*S2h00007111111(1)*ZZ(k,2,l,1)-2*S2h00007211111(1)*ZZ(k,2,l,2
     &  )
       aux000072211(1,1,1)=-(F(1)*S2h0081221111(1))-F(2)*S2h0082221111(1
     &  )+S2h00912211111(k)*Z(1,l)+S2h00922211111(k)*Z(2,l)+(Inv402211
     &  111+S200007221111(1)-S2h00912211111(1)-S2h00922221111(1))*Z(k,l)
     &  -10*(S2h00007122111(1)*ZZ(k,1,l,1)+S2h00007222111(1)*ZZ(k,1,l,2)
     &  )-4*S2h00007121111(1)*ZZ(k,2,l,1)-4*S2h00007221111(1)*ZZ(k,2,l,2
     &  )
       aux000072221(1,1,1)=-(F(1)*S2h0081222111(1))-F(2)*S2h0082222111(1
     &  )+S2h00912221111(k)*Z(1,l)+S2h00922221111(k)*Z(2,l)+(Inv402221
     &  111+S200007222111(1)-S2h00912221111(1)-S2h00922222111(1))*Z(k,l)
     &  -8*(S2h00007122211(1)*ZZ(k,1,l,1)+S2h00007222211(1)*ZZ(k,1,l,2))
     &  -6*S2h00007122111(1)*ZZ(k,2,l,1)-6*S2h00007222111(1)*ZZ(k,2,l,2)
       aux000072222(1,1,1)=-(F(1)*S2h0081222211(1))-F(2)*S2h0082222211(1
     &  )+S2h00912222111(k)*Z(1,l)+S2h00922222111(k)*Z(2,l)+(Inv402222
     &  111+S200007222211(1)-S2h00912222111(1)-S2h00922222211(1))*Z(k,l)
     &  -6*(S2h00007122221(1)*ZZ(k,1,l,1)+S2h00007222221(1)*ZZ(k,1,l,2))
     &  -8*S2h00007122211(1)*ZZ(k,2,l,1)-8*S2h00007222211(1)*ZZ(k,2,l,2)
       aux000072222(2,1,1)=-(F(1)*S2h0081222221(1))-F(2)*S2h0082222221(1
     &  )+S2h00912222211(k)*Z(1,l)+S2h00922222211(k)*Z(2,l)+(Inv402222
     &  211+S200007222221(1)-S2h00912222211(1)-S2h00922222221(1))*Z(k,l)
     &  -4*(S2h00007122222(1)*ZZ(k,1,l,1)+S2h00007222222(1)*ZZ(k,1,l,2))
     &  -10*S2h00007122221(1)*ZZ(k,2,l,1)-10*S2h00007222221(1)*ZZ(k,2,l,
     &  2)
       aux000072222(2,2,1)=-(F(1)*S2h0081222222(1))-F(2)*S2h0082222222(1
     &  )+S2h00912222221(k)*Z(1,l)+S2h00922222221(k)*Z(2,l)+(Inv402222
     &  221+S200007222222(1)-S2h00912222221(1)-S2h00922222222(1))*Z(k,l)
     &  -2*(S2h00007122222(2)*ZZ(k,1,l,1)+S2h00007222222(2)*ZZ(k,1,l,2))
     &  -12*S2h00007122222(1)*ZZ(k,2,l,1)-12*S2h00007222222(1)*ZZ(k,2,l,
     &  2)
       aux000072222(2,2,2)=-(F(1)*S2h0081222222(2))-F(2)*S2h0082222222(2
     &  )+S2h00912222222(k)*Z(1,l)+S2h00922222222(k)*Z(2,l)+(Inv402222
     &  222+S200007222222(2)-S2h00912222222(1)-S2h00922222222(2))*Z(k,l)
     &  -14*(S2h00007122222(2)*ZZ(k,2,l,1)+S2h00007222222(2)*ZZ(k,2,l,2)
     &  )
       temp000071111(1,1,1)=I38Z*(aux000071111(1,1,1)+14*F(4)*temp000061
     &  11(1,1,1)+F(3)*temp0071111(1,1,1)+168*temp000000511(1,1,1)*ZZ(k,
     &  1,l,1))
       temp000072111(1,1,1)=I38Z*(aux000072111(1,1,1)+F(3)*temp0072111(1
     &  ,1,1)+120*temp000000521(1,1,1)*ZZ(k,1,l,1)+2*(r10*temp00006111(1
     &  ,1,1)+6*r21*temp00006211(1,1,1))*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+24*(r
     &  10*temp00006211(1,1,1)*ZZ(k,1,l,1)+temp000000511(1,1,1)*(ZZ(k,1,
     &  l,2)+ZZ(k,2,l,1)))+4*r21*temp00006111(1,1,1)*ZZ(k,2,l,2))
       temp000072211(1,1,1)=I38Z*(aux000072211(1,1,1)+F(3)*temp0072211(1
     &  ,1,1)+80*temp000000522(1,1,1)*ZZ(k,1,l,1)+20*r10*temp00006221(1,
     &  1,1)*ZZ(k,1,l,1)+40*temp000000521(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1
     &  ))+2*(2*r10*temp00006211(1,1,1)+5*r21*temp00006221(1,1,1))*(ZZ(k
     &  ,1,l,2)+ZZ(k,2,l,1))+8*(temp000000511(1,1,1)*ZZ(k,2,l,2)+r21*tem
     &  p00006211(1,1,1)*ZZ(k,2,l,2)))
       temp000072221(1,1,1)=I38Z*(aux000072221(1,1,1)+F(3)*temp0072221(1
     &  ,1,1)+16*r10*temp00006222(1,1,1)*ZZ(k,1,l,1)+2*(3*r10*temp000062
     &  21(1,1,1)+4*r21*temp00006222(1,1,1))*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+4
     &  8*(temp000000522(2,1,1)*ZZ(k,1,l,1)+temp000000522(1,1,1)*(ZZ(k,1
     &  ,l,2)+ZZ(k,2,l,1)))+24*temp000000521(1,1,1)*ZZ(k,2,l,2)+12*r21*t
     &  emp00006221(1,1,1)*ZZ(k,2,l,2))
       temp000072222(1,1,1)=I38Z*(aux000072222(1,1,1)+F(3)*temp0072222(1
     &  ,1,1)+24*temp000000522(2,2,1)*ZZ(k,1,l,1)+12*r10*temp00006222(2,
     &  1,1)*ZZ(k,1,l,1)+2*(4*r10*temp00006222(1,1,1)+3*r21*temp00006222
     &  (2,1,1))*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+16*r21*temp00006222(1,1,1)*ZZ
     &  (k,2,l,2)+48*(temp000000522(2,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+tem
     &  p000000522(1,1,1)*ZZ(k,2,l,2)))
       temp000072222(2,1,1)=I38Z*(aux000072222(2,1,1)+F(3)*temp0072222(2
     &  ,1,1)+8*(temp000000522(2,2,2)*ZZ(k,1,l,1)+r10*temp00006222(2,2,1
     &  )*ZZ(k,1,l,1))+40*temp000000522(2,2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))
     &  +2*(5*r10*temp00006222(2,1,1)+2*r21*temp00006222(2,2,1))*(ZZ(k,1
     &  ,l,2)+ZZ(k,2,l,1))+80*temp000000522(2,1,1)*ZZ(k,2,l,2)+20*r21*te
     &  mp00006222(2,1,1)*ZZ(k,2,l,2))
       temp000072222(2,2,1)=I38Z*(aux000072222(2,2,1)+F(3)*temp0072222(2
     &  ,2,1)+4*r10*temp00006222(2,2,2)*ZZ(k,1,l,1)+2*(6*r10*temp0000622
     &  2(2,2,1)+r21*temp00006222(2,2,2))*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+120*
     &  temp000000522(2,2,1)*ZZ(k,2,l,2)+24*(temp000000522(2,2,2)*(ZZ(k,
     &  1,l,2)+ZZ(k,2,l,1))+r21*temp00006222(2,2,1)*ZZ(k,2,l,2)))
       temp000072222(2,2,2)=I38Z*(aux000072222(2,2,2)+F(3)*temp0072222(2
     &  ,2,2)+168*temp000000522(2,2,2)*ZZ(k,2,l,2)+14*temp00006222(2,2,2
     &  )*(r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+2*r21*ZZ(k,2,l,2)))
       temp000071111(1,1,2)=temp000072111(1,1,1)
       temp000071111(1,2,1)=temp000072111(1,1,1)
       temp000071111(1,2,2)=temp000072211(1,1,1)
       temp000072111(1,1,2)=temp000072211(1,1,1)
       temp000072111(1,2,1)=temp000072211(1,1,1)
       temp000072111(1,2,2)=temp000072221(1,1,1)
       temp000072211(1,1,2)=temp000072221(1,1,1)
       temp000072211(1,2,1)=temp000072221(1,1,1)
       temp000072211(1,2,2)=temp000072222(1,1,1)
       temp000072221(1,1,2)=temp000072222(1,1,1)
       temp000072221(1,2,1)=temp000072222(1,1,1)
       temp000072221(1,2,2)=temp000072222(2,1,1)
       temp000072222(1,1,2)=temp000072222(2,1,1)
       temp000072222(1,2,1)=temp000072222(2,1,1)
       temp000072222(1,2,2)=temp000072222(2,2,1)
       temp000072222(2,1,2)=temp000072222(2,2,1)
       S200911111111(1)=-2*B023
       S200911111111(2)=2*B23(1,1)
       S200912111111(1)=2*B23(1,1)
       S200912111111(2)=-2*B23(1,2)
       S200912211111(1)=-2*B23(1,2)
       S200912211111(2)=2*B23(1,3)
       S200912221111(1)=2*B23(1,3)
       S200912221111(2)=-2*B23(1,4)
       S200912222111(1)=-2*B23(1,4)
       S200912222111(2)=2*B23(1,5)
       S200912222211(1)=2*B23(1,5)
       S200912222211(2)=-2*B23(1,6)
       S200912222221(1)=-2*B23(1,6)
       S200912222221(2)=2*B23(1,7)
       S200912222222(1)=2*B23(1,7)
       S200912222222(2)=-2*B23(1,8)
       S200921111111(1)=2*B23(1,1)
       S200921111111(2)=-2*B23(1,2)
       S200922111111(1)=-2*B23(1,2)
       S200922111111(2)=2*B23(1,3)
       S200922211111(1)=2*B23(1,3)
       S200922211111(2)=-2*B23(1,4)
       S200922221111(1)=-2*B23(1,4)
       S200922221111(2)=2*B23(1,5)
       S200922222111(1)=2*B23(1,5)
       S200922222111(2)=-2*B23(1,6)
       S200922222211(1)=-2*B23(1,6)
       S200922222211(2)=2*B23(1,7)
       S200922222221(1)=2*B23(1,7)
       S200922222221(2)=-2*B23(1,8)
       S200922222222(1)=-2*B23(1,8)
       S200922222222(2)=2*B23(1,9)
       S2111111111111(1)=-B023+B13(1,10)
       S2111111111111(2)=B13(1,10)+B23(1,1)
       S2111211111111(1)=B13(1,10)+B23(1,1)
       S2111211111111(2)=B13(1,10)-B23(1,2)
       S2111221111111(1)=B13(1,10)-B23(1,2)
       S2111221111111(2)=B13(1,10)+B23(1,3)
       S2111222111111(1)=B13(1,10)+B23(1,3)
       S2111222111111(2)=B13(1,10)-B23(1,4)
       S2111222211111(1)=B13(1,10)-B23(1,4)
       S2111222211111(2)=B13(1,10)+B23(1,5)
       S2111222221111(1)=B13(1,10)+B23(1,5)
       S2111222221111(2)=B13(1,10)-B23(1,6)
       S2111222222111(1)=B13(1,10)-B23(1,6)
       S2111222222111(2)=B13(1,10)+B23(1,7)
       S2111222222211(1)=B13(1,10)+B23(1,7)
       S2111222222211(2)=B13(1,10)-B23(1,8)
       S2111222222221(1)=B13(1,10)-B23(1,8)
       S2111222222221(2)=B13(1,10)+B23(1,9)
       S2111222222222(1)=B13(1,10)+B23(1,9)
       S2111222222222(2)=B13(1,10)-B23(1,10)
       S2112111111111(1)=B12(1,10)-B13(1,10)
       S2112111111111(2)=-B13(1,10)
       S2112211111111(1)=-B13(1,10)
       S2112211111111(2)=-B13(1,10)
       S2112221111111(1)=-B13(1,10)
       S2112221111111(2)=-B13(1,10)
       S2112222111111(1)=-B13(1,10)
       S2112222111111(2)=-B13(1,10)
       S2112222211111(1)=-B13(1,10)
       S2112222211111(2)=-B13(1,10)
       S2112222221111(1)=-B13(1,10)
       S2112222221111(2)=-B13(1,10)
       S2112222222111(1)=-B13(1,10)
       S2112222222111(2)=-B13(1,10)
       S2112222222211(1)=-B13(1,10)
       S2112222222211(2)=-B13(1,10)
       S2112222222221(1)=-B13(1,10)
       S2112222222221(2)=-B13(1,10)
       S2112222222222(1)=-B13(1,10)
       S2112222222222(2)=-B13(1,10)
       aux009111111(1,1,1)=-(F(1)*S210111111111(1))-F(2)*S210211111111(1
     &  )+S2111111111111(k)*Z(1,l)+S2112111111111(k)*Z(2,l)+(Inv20111111
     &  111+S200911111111(1)-S2111111111111(1)-S2112211111111(1))*Z(k,l)
     &  -18*(S2h00911111111(1)*ZZ(k,1,l,1)+S2h00921111111(1)*ZZ(k,1,l,2)
     &  )
       aux009211111(1,1,1)=-(F(1)*S210121111111(1))-F(2)*S210221111111(1
     &  )+S2111211111111(k)*Z(1,l)+S2112211111111(k)*Z(2,l)+(Inv20211111
     &  111+S200921111111(1)-S2111211111111(1)-S2112221111111(1))*Z(k,l)
     &  -16*(S2h00912111111(1)*ZZ(k,1,l,1)+S2h00922111111(1)*ZZ(k,1,l,2)
     &  )-2*(S2h00911111111(1)*ZZ(k,2,l,1)+S2h00921111111(1)*ZZ(k,2,l,2)
     &  )
       aux009221111(1,1,1)=-(F(1)*S210122111111(1))-F(2)*S210222111111(1
     &  )+S2111221111111(k)*Z(1,l)+S2112221111111(k)*Z(2,l)+(Inv20221111
     &  111+S200922111111(1)-S2111221111111(1)-S2112222111111(1))*Z(k,l)
     &  -14*(S2h00912211111(1)*ZZ(k,1,l,1)+S2h00922211111(1)*ZZ(k,1,l,2)
     &  )-4*(S2h00912111111(1)*ZZ(k,2,l,1)+S2h00922111111(1)*ZZ(k,2,l,2)
     &  )
       aux009222111(1,1,1)=-(F(1)*S210122211111(1))-F(2)*S210222211111(1
     &  )+S2111222111111(k)*Z(1,l)+S2112222111111(k)*Z(2,l)+(Inv20222111
     &  111+S200922211111(1)-S2111222111111(1)-S2112222211111(1))*Z(k,l)
     &  -12*(S2h00912221111(1)*ZZ(k,1,l,1)+S2h00922221111(1)*ZZ(k,1,l,2)
     &  )-6*(S2h00912211111(1)*ZZ(k,2,l,1)+S2h00922211111(1)*ZZ(k,2,l,2)
     &  )
       aux009222211(1,1,1)=-(F(1)*S210122221111(1))-F(2)*S210222221111(1
     &  )+S2111222211111(k)*Z(1,l)+S2112222211111(k)*Z(2,l)+(Inv20222211
     &  111+S200922221111(1)-S2111222211111(1)-S2112222221111(1))*Z(k,l)
     &  -10*(S2h00912222111(1)*ZZ(k,1,l,1)+S2h00922222111(1)*ZZ(k,1,l,2)
     &  )-8*(S2h00912221111(1)*ZZ(k,2,l,1)+S2h00922221111(1)*ZZ(k,2,l,2)
     &  )
       aux009222221(1,1,1)=-(F(1)*S210122222111(1))-F(2)*S210222222111(1
     &  )+S2111222221111(k)*Z(1,l)+S2112222221111(k)*Z(2,l)+(Inv20222221
     &  111+S200922222111(1)-S2111222221111(1)-S2112222222111(1))*Z(k,l)
     &  -8*(S2h00912222211(1)*ZZ(k,1,l,1)+S2h00922222211(1)*ZZ(k,1,l,2))
     &  -10*(S2h00912222111(1)*ZZ(k,2,l,1)+S2h00922222111(1)*ZZ(k,2,l,2)
     &  )
       aux009222222(1,1,1)=-(F(1)*S210122222211(1))-F(2)*S210222222211(1
     &  )+S2111222222111(k)*Z(1,l)+S2112222222111(k)*Z(2,l)+(Inv20222222
     &  111+S200922222211(1)-S2111222222111(1)-S2112222222211(1))*Z(k,l)
     &  -6*(S2h00912222221(1)*ZZ(k,1,l,1)+S2h00922222221(1)*ZZ(k,1,l,2))
     &  -12*(S2h00912222211(1)*ZZ(k,2,l,1)+S2h00922222211(1)*ZZ(k,2,l,2)
     &  )
       aux009222222(2,1,1)=-(F(1)*S210122222221(1))-F(2)*S210222222221(1
     &  )+S2111222222211(k)*Z(1,l)+S2112222222211(k)*Z(2,l)+(Inv20222222
     &  211+S200922222221(1)-S2111222222211(1)-S2112222222221(1))*Z(k,l)
     &  -4*(S2h00912222222(1)*ZZ(k,1,l,1)+S2h00922222222(1)*ZZ(k,1,l,2))
     &  -14*(S2h00912222221(1)*ZZ(k,2,l,1)+S2h00922222221(1)*ZZ(k,2,l,2)
     &  )
       aux009222222(2,2,1)=-(F(1)*S210122222222(1))-F(2)*S210222222222(1
     &  )+S2111222222221(k)*Z(1,l)+S2112222222221(k)*Z(2,l)+(Inv20222222
     &  221+S200922222222(1)-S2111222222221(1)-S2112222222222(1))*Z(k,l)
     &  -2*(S2h00912222222(2)*ZZ(k,1,l,1)+S2h00922222222(2)*ZZ(k,1,l,2))
     &  -16*(S2h00912222222(1)*ZZ(k,2,l,1)+S2h00922222222(1)*ZZ(k,2,l,2)
     &  )
       aux009222222(2,2,2)=-(F(1)*S210122222222(2))-F(2)*S210222222222(2
     &  )+S2111222222222(k)*Z(1,l)+S2112222222222(k)*Z(2,l)+(Inv20222222
     &  222+S200922222222(2)-S2111222222222(1)-S2112222222222(2))*Z(k,l)
     &  -18*(S2h00912222222(2)*ZZ(k,2,l,1)+S2h00922222222(2)*ZZ(k,2,l,2)
     &  )
       temp009111111(1,1,1)=I42Z*(aux009111111(1,1,1)+18*F(4)*temp008111
     &  11(1,1,1)+F(3)*temp9111111(1,1,1)+288*temp000071111(1,1,1)*ZZ(k,
     &  1,l,1))
       temp009211111(1,1,1)=I42Z*(aux009211111(1,1,1)+F(5)*temp00811111(
     &  1,1,1)+16*F(4)*temp00821111(1,1,1)+F(3)*temp9211111(1,1,1)+224*t
     &  emp000072111(1,1,1)*ZZ(k,1,l,1)+32*temp000071111(1,1,1)*(ZZ(k,1,
     &  l,2)+ZZ(k,2,l,1)))
       temp009221111(1,1,1)=I42Z*(aux009221111(1,1,1)+2*F(5)*temp0082111
     &  1(1,1,1)+14*F(4)*temp00822111(1,1,1)+F(3)*temp9221111(1,1,1)+168
     &  *temp000072211(1,1,1)*ZZ(k,1,l,1)+56*temp000072111(1,1,1)*(ZZ(k,
     &  1,l,2)+ZZ(k,2,l,1))+8*temp000071111(1,1,1)*ZZ(k,2,l,2))
       temp009222111(1,1,1)=I42Z*(aux009222111(1,1,1)+12*F(4)*temp008222
     &  11(1,1,1)+F(3)*temp9222111(1,1,1)+120*temp000072221(1,1,1)*ZZ(k,
     &  1,l,1)+72*temp000072211(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+24*temp
     &  000072111(1,1,1)*ZZ(k,2,l,2)+temp00822111(1,1,1)*(6*r10*(ZZ(k,1,
     &  l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2)))
       temp009222211(1,1,1)=I42Z*(aux009222211(1,1,1)+4*F(5)*temp0082221
     &  1(1,1,1)+10*F(4)*temp00822221(1,1,1)+F(3)*temp9222211(1,1,1)+80*
     &  (temp000072222(1,1,1)*ZZ(k,1,l,1)+temp000072221(1,1,1)*(ZZ(k,1,l
     &  ,2)+ZZ(k,2,l,1)))+48*temp000072211(1,1,1)*ZZ(k,2,l,2))
       temp009222221(1,1,1)=I42Z*(aux009222221(1,1,1)+8*F(4)*temp0082222
     &  2(1,1,1)+F(3)*temp9222221(1,1,1)+48*temp000072222(2,1,1)*ZZ(k,1,
     &  l,1)+temp00822221(1,1,1)*(10*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+20*r2
     &  1*ZZ(k,2,l,2))+80*(temp000072222(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)
     &  )+temp000072221(1,1,1)*ZZ(k,2,l,2)))
       temp009222222(1,1,1)=I42Z*(aux009222222(1,1,1)+6*F(4)*temp0082222
     &  2(2,1,1)+F(3)*temp9222222(1,1,1)+24*temp000072222(2,2,1)*ZZ(k,1,
     &  l,1)+72*temp000072222(2,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+120*temp0
     &  00072222(1,1,1)*ZZ(k,2,l,2)+temp00822222(1,1,1)*(12*r10*(ZZ(k,1,
     &  l,2)+ZZ(k,2,l,1))+24*r21*ZZ(k,2,l,2)))
       temp009222222(2,1,1)=I42Z*(aux009222222(2,1,1)+4*F(4)*temp0082222
     &  2(2,2,1)+F(3)*temp9222222(2,1,1)+8*temp000072222(2,2,2)*ZZ(k,1,l
     &  ,1)+56*temp000072222(2,2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+168*temp00
     &  0072222(2,1,1)*ZZ(k,2,l,2)+temp00822222(2,1,1)*(14*r10*(ZZ(k,1,l
     &  ,2)+ZZ(k,2,l,1))+28*r21*ZZ(k,2,l,2)))
       temp009222222(2,2,1)=I42Z*(aux009222222(2,2,1)+8*F(5)*temp0082222
     &  2(2,2,1)+2*F(4)*temp00822222(2,2,2)+F(3)*temp9222222(2,2,1)+32*t
     &  emp000072222(2,2,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+224*temp000072222(
     &  2,2,1)*ZZ(k,2,l,2))
       temp009222222(2,2,2)=I42Z*(aux009222222(2,2,2)+F(3)*temp9222222(2
     &  ,2,2)+288*temp000072222(2,2,2)*ZZ(k,2,l,2)+temp00822222(2,2,2)*(
     &  18*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+36*r21*ZZ(k,2,l,2)))
       temp009111111(1,1,2)=temp009211111(1,1,1)
       temp009111111(1,2,1)=temp009211111(1,1,1)
       temp009111111(1,2,2)=temp009221111(1,1,1)
       temp009211111(1,1,2)=temp009221111(1,1,1)
       temp009211111(1,2,1)=temp009221111(1,1,1)
       temp009211111(1,2,2)=temp009222111(1,1,1)
       temp009221111(1,1,2)=temp009222111(1,1,1)
       temp009221111(1,2,1)=temp009222111(1,1,1)
       temp009221111(1,2,2)=temp009222211(1,1,1)
       temp009222111(1,1,2)=temp009222211(1,1,1)
       temp009222111(1,2,1)=temp009222211(1,1,1)
       temp009222111(1,2,2)=temp009222221(1,1,1)
       temp009222211(1,1,2)=temp009222221(1,1,1)
       temp009222211(1,2,1)=temp009222221(1,1,1)
       temp009222211(1,2,2)=temp009222222(1,1,1)
       temp009222221(1,1,2)=temp009222222(1,1,1)
       temp009222221(1,2,1)=temp009222222(1,1,1)
       temp009222221(1,2,2)=temp009222222(2,1,1)
       temp009222222(1,1,2)=temp009222222(2,1,1)
       temp009222222(1,2,1)=temp009222222(2,1,1)
       temp009222222(1,2,2)=temp009222222(2,2,1)
       temp009222222(2,1,2)=temp009222222(2,2,1)
       aux101111111(1,1,1)=-(S2111111111111(1)*Z(jj,1))-S2112111111111(1
     &  )*Z(jj,2)
       aux102111111(1,1,1)=-(S2111211111111(1)*Z(jj,1))-S2112211111111(1
     &  )*Z(jj,2)
       aux102211111(1,1,1)=-(S2111221111111(1)*Z(jj,1))-S2112221111111(1
     &  )*Z(jj,2)
       aux102221111(1,1,1)=-(S2111222111111(1)*Z(jj,1))-S2112222111111(1
     &  )*Z(jj,2)
       aux102222111(1,1,1)=-(S2111222211111(1)*Z(jj,1))-S2112222211111(1
     &  )*Z(jj,2)
       aux102222211(1,1,1)=-(S2111222221111(1)*Z(jj,1))-S2112222221111(1
     &  )*Z(jj,2)
       aux102222221(1,1,1)=-(S2111222222111(1)*Z(jj,1))-S2112222222111(1
     &  )*Z(jj,2)
       aux102222222(1,1,1)=-(S2111222222211(1)*Z(jj,1))-S2112222222211(1
     &  )*Z(jj,2)
       aux102222222(2,1,1)=-(S2111222222221(1)*Z(jj,1))-S2112222222221(1
     &  )*Z(jj,2)
       aux102222222(2,2,1)=-(S2111222222222(1)*Z(jj,1))-S2112222222222(1
     &  )*Z(jj,2)
       aux102222222(2,2,2)=-(S2111222222222(2)*Z(jj,1))-S2112222222222(2
     &  )*Z(jj,2)
       temp101111111(1,1,1)=IX*(aux101111111(1,1,1)+20*temp009111111(1,1
     &  ,1)*Z(jj,1))
       temp102111111(1,1,1)=IX*(aux102111111(1,1,1)+18*temp009211111(1,1
     &  ,1)*Z(jj,1)+2*temp009111111(1,1,1)*Z(jj,2))
       temp102211111(1,1,1)=IX*(aux102211111(1,1,1)+16*temp009221111(1,1
     &  ,1)*Z(jj,1)+4*temp009211111(1,1,1)*Z(jj,2))
       temp102221111(1,1,1)=IX*(aux102221111(1,1,1)+14*temp009222111(1,1
     &  ,1)*Z(jj,1)+6*temp009221111(1,1,1)*Z(jj,2))
       temp102222111(1,1,1)=IX*(aux102222111(1,1,1)+12*temp009222211(1,1
     &  ,1)*Z(jj,1)+8*temp009222111(1,1,1)*Z(jj,2))
       temp102222211(1,1,1)=IX*(aux102222211(1,1,1)+10*(temp009222221(1,
     &  1,1)*Z(jj,1)+temp009222211(1,1,1)*Z(jj,2)))
       temp102222221(1,1,1)=IX*(aux102222221(1,1,1)+8*temp009222222(1,1,
     &  1)*Z(jj,1)+12*temp009222221(1,1,1)*Z(jj,2))
       temp102222222(1,1,1)=IX*(aux102222222(1,1,1)+6*temp009222222(2,1,
     &  1)*Z(jj,1)+14*temp009222222(1,1,1)*Z(jj,2))
       temp102222222(2,1,1)=IX*(aux102222222(2,1,1)+4*temp009222222(2,2,
     &  1)*Z(jj,1)+16*temp009222222(2,1,1)*Z(jj,2))
       temp102222222(2,2,1)=IX*(aux102222222(2,2,1)+2*temp009222222(2,2,
     &  2)*Z(jj,1)+18*temp009222222(2,2,1)*Z(jj,2))
       temp102222222(2,2,2)=IX*(aux102222222(2,2,2)+20*temp009222222(2,2
     &  ,2)*Z(jj,2))
       temp101111111(1,1,2)=temp102111111(1,1,1)
       temp101111111(1,2,1)=temp102111111(1,1,1)
       temp101111111(1,2,2)=temp102211111(1,1,1)
       temp102111111(1,1,2)=temp102211111(1,1,1)
       temp102111111(1,2,1)=temp102211111(1,1,1)
       temp102111111(1,2,2)=temp102221111(1,1,1)
       temp102211111(1,1,2)=temp102221111(1,1,1)
       temp102211111(1,2,1)=temp102221111(1,1,1)
       temp102211111(1,2,2)=temp102222111(1,1,1)
       temp102221111(1,1,2)=temp102222111(1,1,1)
       temp102221111(1,2,1)=temp102222111(1,1,1)
       temp102221111(1,2,2)=temp102222211(1,1,1)
       temp102222111(1,1,2)=temp102222211(1,1,1)
       temp102222111(1,2,1)=temp102222211(1,1,1)
       temp102222111(1,2,2)=temp102222221(1,1,1)
       temp102222211(1,1,2)=temp102222221(1,1,1)
       temp102222211(1,2,1)=temp102222221(1,1,1)
       temp102222211(1,2,2)=temp102222222(1,1,1)
       temp102222221(1,1,2)=temp102222222(1,1,1)
       temp102222221(1,2,1)=temp102222222(1,1,1)
       temp102222221(1,2,2)=temp102222222(2,1,1)
       temp102222222(1,1,2)=temp102222222(2,1,1)
       temp102222222(1,2,1)=temp102222222(2,1,1)
       temp102222222(1,2,2)=temp102222222(2,2,1)
       temp102222222(2,1,2)=temp102222222(2,2,1)
c                Step2
       tempC30000000000=I22Z*(auxC30000000000+tempC300000000*F(3)-det3*t
     &  emp000000002(k,l))
       temp000000002(1,1)=I26Z*(aux000000002(1,1)+4*F(4)*temp000000001(1
     &  )+F(3)*temp0000002(1,1)-det3*temp00000041(1,k,l)+8*tempC30000000
     &  000*ZZ(k,1,l,1))
       temp000000002(2,1)=I26Z*(aux000000002(2,1)+F(5)*temp000000001(1)+
     &  2*F(4)*temp000000001(2)+F(3)*temp0000002(2,1)-det3*temp00000042(
     &  1,k,l)+4*tempC30000000000*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp000000002(2,2)=I26Z*(aux000000002(2,2)+2*F(5)*temp000000001(2
     &  )+F(3)*temp0000002(2,2)-det3*temp00000042(2,k,l)+8*tempC30000000
     &  000*ZZ(k,2,l,2))
       temp000000002(1,2)=temp000000002(2,1)
       temp00000041(1,1,1)=I30Z*(aux00000041(1,1,1)+8*F(4)*temp0000003(1
     &  ,1,1)+F(3)*temp000041(1,1,1)-det3*temp00006111(1,k,l)+48*temp000
     &  000002(1,1)*ZZ(k,1,l,1))
       temp00000042(1,1,1)=I30Z*(aux00000042(1,1,1)+F(5)*temp0000003(1,1
     &  ,1)+6*F(4)*temp0000003(2,1,1)+F(3)*temp000042(1,1,1)-det3*temp00
     &  006211(1,k,l)+24*temp000000002(2,1)*ZZ(k,1,l,1)+12*temp000000002
     &  (1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp00000042(2,1,1)=I30Z*(aux00000042(2,1,1)+2*F(5)*temp0000003(2
     &  ,1,1)+4*F(4)*temp0000003(2,2,1)+F(3)*temp000042(2,1,1)-det3*temp
     &  00006221(1,k,l)+16*temp000000002(2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+
     &  8*(temp000000002(2,2)*ZZ(k,1,l,1)+temp000000002(1,1)*ZZ(k,2,l,2)
     &  ))
       temp00000042(2,2,1)=I30Z*(aux00000042(2,2,1)+2*F(4)*temp0000003(2
     &  ,2,2)+F(3)*temp000042(2,2,1)-det3*temp00006222(1,k,l)+12*temp000
     &  000002(2,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+24*temp000000002(2,1)*ZZ(k
     &  ,2,l,2)+temp0000003(2,2,1)*(6*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+12*r
     &  21*ZZ(k,2,l,2)))
       temp00000042(2,2,2)=I30Z*(aux00000042(2,2,2)+4*F(5)*temp0000003(2
     &  ,2,2)+F(3)*temp000042(2,2,2)-det3*temp00006222(2,k,l)+48*temp000
     &  000002(2,2)*ZZ(k,2,l,2))
       temp00000041(1,1,2)=temp00000042(1,1,1)
       temp00000041(1,2,1)=temp00000042(1,1,1)
       temp00000041(1,2,2)=temp00000042(2,1,1)
       temp00000042(1,1,2)=temp00000042(2,1,1)
       temp00000042(1,2,1)=temp00000042(2,1,1)
       temp00000042(1,2,2)=temp00000042(2,2,1)
       temp00000042(2,1,2)=temp00000042(2,2,1)
       temp00006111(1,1,1)=I34Z*(aux00006111(1,1,1)+12*F(4)*temp0000511(
     &  1,1,1)+F(3)*temp006111(1,1,1)-det3*temp00811111(1,k,l)+120*temp0
     &  0000041(1,1,1)*ZZ(k,1,l,1))
       temp00006211(1,1,1)=I34Z*(aux00006211(1,1,1)+F(5)*temp0000511(1,1
     &  ,1)+10*F(4)*temp0000521(1,1,1)+F(3)*temp006211(1,1,1)-det3*temp0
     &  0821111(1,k,l)+80*temp00000042(1,1,1)*ZZ(k,1,l,1)+20*temp0000004
     &  1(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp00006221(1,1,1)=I34Z*(aux00006221(1,1,1)+2*F(5)*temp0000521(1
     &  ,1,1)+F(3)*temp006221(1,1,1)-det3*temp00822111(1,k,l)+48*temp000
     &  00042(2,1,1)*ZZ(k,1,l,1)+32*temp00000042(1,1,1)*(ZZ(k,1,l,2)+ZZ(
     &  k,2,l,1))+8*(F(4)*temp0000522(1,1,1)+temp00000041(1,1,1)*ZZ(k,2,
     &  l,2)))
       temp00006222(1,1,1)=I34Z*(aux00006222(1,1,1)+6*F(4)*temp0000522(2
     &  ,1,1)+F(3)*temp006222(1,1,1)-det3*temp00822211(1,k,l)+36*temp000
     &  00042(2,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp0000522(1,1,1)*(6*r10
     &  *(ZZ(k,1,l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2))+24*(temp00000042(
     &  2,2,1)*ZZ(k,1,l,1)+temp00000042(1,1,1)*ZZ(k,2,l,2)))
       temp00006222(2,1,1)=I34Z*(aux00006222(2,1,1)+4*F(5)*temp0000522(2
     &  ,1,1)+4*F(4)*temp0000522(2,2,1)+F(3)*temp006222(2,1,1)-det3*temp
     &  00822221(1,k,l)+8*temp00000042(2,2,2)*ZZ(k,1,l,1)+32*temp0000004
     &  2(2,2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+48*temp00000042(2,1,1)*ZZ(k,2
     &  ,l,2))
       temp00006222(2,2,1)=I34Z*(aux00006222(2,2,1)+2*F(4)*temp0000522(2
     &  ,2,2)+F(3)*temp006222(2,2,1)-det3*temp00822222(1,k,l)+20*temp000
     &  00042(2,2,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+80*temp00000042(2,2,1)*ZZ
     &  (k,2,l,2)+temp0000522(2,2,1)*(10*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+2
     &  0*r21*ZZ(k,2,l,2)))
       temp00006222(2,2,2)=I34Z*(aux00006222(2,2,2)+F(3)*temp006222(2,2,
     &  2)-det3*temp00822222(2,k,l)+120*temp00000042(2,2,2)*ZZ(k,2,l,2)+
     &  temp0000522(2,2,2)*(12*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+24*r21*ZZ(k
     &  ,2,l,2)))
       temp00006111(1,1,2)=temp00006211(1,1,1)
       temp00006111(1,2,1)=temp00006211(1,1,1)
       temp00006111(1,2,2)=temp00006221(1,1,1)
       temp00006211(1,1,2)=temp00006221(1,1,1)
       temp00006211(1,2,1)=temp00006221(1,1,1)
       temp00006211(1,2,2)=temp00006222(1,1,1)
       temp00006221(1,1,2)=temp00006222(1,1,1)
       temp00006221(1,2,1)=temp00006222(1,1,1)
       temp00006221(1,2,2)=temp00006222(2,1,1)
       temp00006222(1,1,2)=temp00006222(2,1,1)
       temp00006222(1,2,1)=temp00006222(2,1,1)
       temp00006222(1,2,2)=temp00006222(2,2,1)
       temp00006222(2,1,2)=temp00006222(2,2,1)
       temp00811111(1,1,1)=I38Z*(aux00811111(1,1,1)+16*F(4)*temp0071111(
     &  1,1,1)-det3*temp101111111(1,k,l)+F(3)*temp811111(1,1,1)+224*temp
     &  00006111(1,1,1)*ZZ(k,1,l,1))
       temp00821111(1,1,1)=I38Z*(aux00821111(1,1,1)+F(5)*temp0071111(1,1
     &  ,1)+14*F(4)*temp0072111(1,1,1)-det3*temp102111111(1,k,l)+F(3)*te
     &  mp821111(1,1,1)+168*temp00006211(1,1,1)*ZZ(k,1,l,1)+28*temp00006
     &  111(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp00822111(1,1,1)=I38Z*(aux00822111(1,1,1)+2*F(5)*temp0072111(1
     &  ,1,1)+12*F(4)*temp0072211(1,1,1)-det3*temp102211111(1,k,l)+F(3)*
     &  temp822111(1,1,1)+120*temp00006221(1,1,1)*ZZ(k,1,l,1)+48*temp000
     &  06211(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+8*temp00006111(1,1,1)*ZZ(
     &  k,2,l,2))
       temp00822211(1,1,1)=I38Z*(aux00822211(1,1,1)+10*F(4)*temp0072221(
     &  1,1,1)-det3*temp102221111(1,k,l)+F(3)*temp822211(1,1,1)+80*temp0
     &  0006222(1,1,1)*ZZ(k,1,l,1)+60*temp00006221(1,1,1)*(ZZ(k,1,l,2)+Z
     &  Z(k,2,l,1))+24*temp00006211(1,1,1)*ZZ(k,2,l,2)+temp0072211(1,1,1
     &  )*(6*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2)))
       temp00822221(1,1,1)=I38Z*(aux00822221(1,1,1)+4*F(5)*temp0072221(1
     &  ,1,1)+8*F(4)*temp0072222(1,1,1)-det3*temp102222111(1,k,l)+F(3)*t
     &  emp822221(1,1,1)+64*temp00006222(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)
     &  )+48*(temp00006222(2,1,1)*ZZ(k,1,l,1)+temp00006221(1,1,1)*ZZ(k,2
     &  ,l,2)))
       temp00822222(1,1,1)=I38Z*(aux00822222(1,1,1)+6*F(4)*temp0072222(2
     &  ,1,1)-det3*temp102222211(1,k,l)+F(3)*temp822222(1,1,1)+24*temp00
     &  006222(2,2,1)*ZZ(k,1,l,1)+60*temp00006222(2,1,1)*(ZZ(k,1,l,2)+ZZ
     &  (k,2,l,1))+80*temp00006222(1,1,1)*ZZ(k,2,l,2)+temp0072222(1,1,1)
     &  *(10*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+20*r21*ZZ(k,2,l,2)))
       temp00822222(2,1,1)=I38Z*(aux00822222(2,1,1)+4*F(4)*temp0072222(2
     &  ,2,1)-det3*temp102222221(1,k,l)+F(3)*temp822222(2,1,1)+8*temp000
     &  06222(2,2,2)*ZZ(k,1,l,1)+48*temp00006222(2,2,1)*(ZZ(k,1,l,2)+ZZ(
     &  k,2,l,1))+120*temp00006222(2,1,1)*ZZ(k,2,l,2)+temp0072222(2,1,1)
     &  *(12*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+24*r21*ZZ(k,2,l,2)))
       temp00822222(2,2,1)=I38Z*(aux00822222(2,2,1)+2*F(4)*temp0072222(2
     &  ,2,2)-det3*temp102222222(1,k,l)+F(3)*temp822222(2,2,1)+28*temp00
     &  006222(2,2,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+168*temp00006222(2,2,1)*
     &  ZZ(k,2,l,2)+temp0072222(2,2,1)*(14*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))
     &  +28*r21*ZZ(k,2,l,2)))
       temp00822222(2,2,2)=I38Z*(aux00822222(2,2,2)+8*F(5)*temp0072222(2
     &  ,2,2)-det3*temp102222222(2,k,l)+F(3)*temp822222(2,2,2)+224*temp0
     &  0006222(2,2,2)*ZZ(k,2,l,2))
       temp00811111(1,1,2)=temp00821111(1,1,1)
       temp00811111(1,2,1)=temp00821111(1,1,1)
       temp00811111(1,2,2)=temp00822111(1,1,1)
       temp00821111(1,1,2)=temp00822111(1,1,1)
       temp00821111(1,2,1)=temp00822111(1,1,1)
       temp00821111(1,2,2)=temp00822211(1,1,1)
       temp00822111(1,1,2)=temp00822211(1,1,1)
       temp00822111(1,2,1)=temp00822211(1,1,1)
       temp00822111(1,2,2)=temp00822221(1,1,1)
       temp00822211(1,1,2)=temp00822221(1,1,1)
       temp00822211(1,2,1)=temp00822221(1,1,1)
       temp00822211(1,2,2)=temp00822222(1,1,1)
       temp00822221(1,1,2)=temp00822222(1,1,1)
       temp00822221(1,2,1)=temp00822222(1,1,1)
       temp00822221(1,2,2)=temp00822222(2,1,1)
       temp00822222(1,1,2)=temp00822222(2,1,1)
       temp00822222(1,2,1)=temp00822222(2,1,1)
       temp00822222(1,2,2)=temp00822222(2,2,1)
       temp00822222(2,1,2)=temp00822222(2,2,1)
       temp9111111(1,1,1)=IX*(aux9111111(1,1,1)+det3*temp101111111(1,1,j
     &  j)+18*temp00811111(1,1,1)*Z(jj,1))
       temp9211111(1,1,1)=IX*(aux9211111(1,1,1)+det3*temp102111111(1,1,j
     &  j)+16*temp00821111(1,1,1)*Z(jj,1)+2*temp00811111(1,1,1)*Z(jj,2))
       temp9221111(1,1,1)=IX*(aux9221111(1,1,1)+det3*temp102211111(1,1,j
     &  j)+14*temp00822111(1,1,1)*Z(jj,1)+4*temp00821111(1,1,1)*Z(jj,2))
       temp9222111(1,1,1)=IX*(aux9222111(1,1,1)+det3*temp102221111(1,1,j
     &  j)+12*temp00822211(1,1,1)*Z(jj,1)+6*temp00822111(1,1,1)*Z(jj,2))
       temp9222211(1,1,1)=IX*(aux9222211(1,1,1)+det3*temp102222111(1,1,j
     &  j)+10*temp00822221(1,1,1)*Z(jj,1)+8*temp00822211(1,1,1)*Z(jj,2))
       temp9222221(1,1,1)=IX*(aux9222221(1,1,1)+det3*temp102222211(1,1,j
     &  j)+8*temp00822222(1,1,1)*Z(jj,1)+10*temp00822221(1,1,1)*Z(jj,2))
       temp9222222(1,1,1)=IX*(aux9222222(1,1,1)+det3*temp102222221(1,1,j
     &  j)+6*temp00822222(2,1,1)*Z(jj,1)+12*temp00822222(1,1,1)*Z(jj,2))
       temp9222222(2,1,1)=IX*(aux9222222(2,1,1)+det3*temp102222222(1,1,j
     &  j)+4*temp00822222(2,2,1)*Z(jj,1)+14*temp00822222(2,1,1)*Z(jj,2))
       temp9222222(2,2,1)=IX*(aux9222222(2,2,1)+det3*temp102222222(2,1,j
     &  j)+2*temp00822222(2,2,2)*Z(jj,1)+16*temp00822222(2,2,1)*Z(jj,2))
       temp9222222(2,2,2)=IX*(aux9222222(2,2,2)+det3*temp102222222(2,2,j
     &  j)+18*temp00822222(2,2,2)*Z(jj,2))
       temp9111111(1,1,2)=temp9211111(1,1,1)
       temp9111111(1,2,1)=temp9211111(1,1,1)
       temp9111111(1,2,2)=temp9221111(1,1,1)
       temp9211111(1,1,2)=temp9221111(1,1,1)
       temp9211111(1,2,1)=temp9221111(1,1,1)
       temp9211111(1,2,2)=temp9222111(1,1,1)
       temp9221111(1,1,2)=temp9222111(1,1,1)
       temp9221111(1,2,1)=temp9222111(1,1,1)
       temp9221111(1,2,2)=temp9222211(1,1,1)
       temp9222111(1,1,2)=temp9222211(1,1,1)
       temp9222111(1,2,1)=temp9222211(1,1,1)
       temp9222111(1,2,2)=temp9222221(1,1,1)
       temp9222211(1,1,2)=temp9222221(1,1,1)
       temp9222211(1,2,1)=temp9222221(1,1,1)
       temp9222211(1,2,2)=temp9222222(1,1,1)
       temp9222221(1,1,2)=temp9222222(1,1,1)
       temp9222221(1,2,1)=temp9222222(1,1,1)
       temp9222221(1,2,2)=temp9222222(2,1,1)
       temp9222222(1,1,2)=temp9222222(2,1,1)
       temp9222222(1,2,1)=temp9222222(2,1,1)
       temp9222222(1,2,2)=temp9222222(2,2,1)
       temp9222222(2,1,2)=temp9222222(2,2,1)
c                Step3
       temp000000001(1)=I22Z*(aux000000001(1)+2*tempC300000000*F(4)+F(3)
     &  *temp0000001(1)-det3*temp0000003(1,k,l))
       temp000000001(2)=I22Z*(aux000000001(2)+tempC300000000*F(5)+F(3)*t
     &  emp0000001(2)-det3*temp0000003(2,k,l))
       temp0000003(1,1,1)=I26Z*(aux0000003(1,1,1)+6*F(4)*temp0000002(1,1
     &  )+F(3)*temp00003(1,1,1)-det3*temp0000511(1,k,l)+24*temp000000001
     &  (1)*ZZ(k,1,l,1))
       temp0000003(2,1,1)=I26Z*(aux0000003(2,1,1)+F(5)*temp0000002(1,1)+
     &  4*F(4)*temp0000002(2,1)+F(3)*temp00003(2,1,1)-det3*temp0000521(1
     &  ,k,l)+8*(temp000000001(2)*ZZ(k,1,l,1)+temp000000001(1)*(ZZ(k,1,l
     &  ,2)+ZZ(k,2,l,1))))
       temp0000003(2,2,1)=I26Z*(aux0000003(2,2,1)+2*F(5)*temp0000002(2,1
     &  )+2*F(4)*temp0000002(2,2)+F(3)*temp00003(2,2,1)-det3*temp0000522
     &  (1,k,l)+8*(temp000000001(2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp000000
     &  001(1)*ZZ(k,2,l,2)))
       temp0000003(2,2,2)=I26Z*(aux0000003(2,2,2)+F(3)*temp00003(2,2,2)-
     &  det3*temp0000522(2,k,l)+24*temp000000001(2)*ZZ(k,2,l,2)+temp0000
     &  002(2,2)*(6*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2)))
       temp0000003(1,1,2)=temp0000003(2,1,1)
       temp0000003(1,2,1)=temp0000003(2,1,1)
       temp0000003(1,2,2)=temp0000003(2,2,1)
       temp0000003(2,1,2)=temp0000003(2,2,1)
       temp0000511(1,1,1)=I30Z*(aux0000511(1,1,1)+10*F(4)*temp000041(1,1
     &  ,1)+F(3)*temp00511(1,1,1)-det3*temp0071111(1,k,l)+80*temp0000003
     &  (1,1,1)*ZZ(k,1,l,1))
       temp0000521(1,1,1)=I30Z*(aux0000521(1,1,1)+F(5)*temp000041(1,1,1)
     &  +8*F(4)*temp000042(1,1,1)+F(3)*temp00521(1,1,1)-det3*temp0072111
     &  (1,k,l)+48*temp0000003(2,1,1)*ZZ(k,1,l,1)+16*temp0000003(1,1,1)*
     &  (ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp0000522(1,1,1)=I30Z*(aux0000522(1,1,1)+2*F(5)*temp000042(1,1,
     &  1)+6*F(4)*temp000042(2,1,1)+F(3)*temp00522(1,1,1)-det3*temp00722
     &  11(1,k,l)+24*(temp0000003(2,2,1)*ZZ(k,1,l,1)+temp0000003(2,1,1)*
     &  (ZZ(k,1,l,2)+ZZ(k,2,l,1)))+8*temp0000003(1,1,1)*ZZ(k,2,l,2))
       temp0000522(2,1,1)=I30Z*(aux0000522(2,1,1)+4*F(4)*temp000042(2,2,
     &  1)+F(3)*temp00522(2,1,1)-det3*temp0072221(1,k,l)+8*temp0000003(2
     &  ,2,2)*ZZ(k,1,l,1)+temp000042(2,1,1)*(6*r10*(ZZ(k,1,l,2)+ZZ(k,2,l
     &  ,1))+12*r21*ZZ(k,2,l,2))+24*(temp0000003(2,2,1)*(ZZ(k,1,l,2)+ZZ(
     &  k,2,l,1))+temp0000003(2,1,1)*ZZ(k,2,l,2)))
       temp0000522(2,2,1)=I30Z*(aux0000522(2,2,1)+4*F(5)*temp000042(2,2,
     &  1)+2*F(4)*temp000042(2,2,2)+F(3)*temp00522(2,2,1)-det3*temp00722
     &  22(1,k,l)+16*temp0000003(2,2,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+48*tem
     &  p0000003(2,2,1)*ZZ(k,2,l,2))
       temp0000522(2,2,2)=I30Z*(aux0000522(2,2,2)+F(3)*temp00522(2,2,2)-
     &  det3*temp0072222(2,k,l)+80*temp0000003(2,2,2)*ZZ(k,2,l,2)+temp00
     &  0042(2,2,2)*(10*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+20*r21*ZZ(k,2,l,2)
     &  ))
       temp0000511(1,1,2)=temp0000521(1,1,1)
       temp0000511(1,2,1)=temp0000521(1,1,1)
       temp0000511(1,2,2)=temp0000522(1,1,1)
       temp0000521(1,1,2)=temp0000522(1,1,1)
       temp0000521(1,2,1)=temp0000522(1,1,1)
       temp0000521(1,2,2)=temp0000522(2,1,1)
       temp0000522(1,1,2)=temp0000522(2,1,1)
       temp0000522(1,2,1)=temp0000522(2,1,1)
       temp0000522(1,2,2)=temp0000522(2,2,1)
       temp0000522(2,1,2)=temp0000522(2,2,1)
       temp0071111(1,1,1)=I34Z*(aux0071111(1,1,1)+14*F(4)*temp006111(1,1
     &  ,1)+F(3)*temp71111(1,1,1)-det3*temp9111111(1,k,l)+168*temp000051
     &  1(1,1,1)*ZZ(k,1,l,1))
       temp0072111(1,1,1)=I34Z*(aux0072111(1,1,1)+F(5)*temp006111(1,1,1)
     &  +12*F(4)*temp006211(1,1,1)+F(3)*temp72111(1,1,1)-det3*temp921111
     &  1(1,k,l)+120*temp0000521(1,1,1)*ZZ(k,1,l,1)+24*temp0000511(1,1,1
     &  )*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp0072211(1,1,1)=I34Z*(aux0072211(1,1,1)+2*F(5)*temp006211(1,1,
     &  1)+10*F(4)*temp006221(1,1,1)+F(3)*temp72211(1,1,1)-det3*temp9221
     &  111(1,k,l)+80*temp0000522(1,1,1)*ZZ(k,1,l,1)+40*temp0000521(1,1,
     &  1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+8*temp0000511(1,1,1)*ZZ(k,2,l,2))
       temp0072221(1,1,1)=I34Z*(aux0072221(1,1,1)+8*F(4)*temp006222(1,1,
     &  1)+F(3)*temp72221(1,1,1)-det3*temp9222111(1,k,l)+48*(temp0000522
     &  (2,1,1)*ZZ(k,1,l,1)+temp0000522(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))
     &  )+24*temp0000521(1,1,1)*ZZ(k,2,l,2)+temp006221(1,1,1)*(6*r10*(ZZ
     &  (k,1,l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2)))
       temp0072222(1,1,1)=I34Z*(aux0072222(1,1,1)+4*F(5)*temp006222(1,1,
     &  1)+6*F(4)*temp006222(2,1,1)+F(3)*temp72222(1,1,1)-det3*temp92222
     &  11(1,k,l)+24*temp0000522(2,2,1)*ZZ(k,1,l,1)+48*(temp0000522(2,1,
     &  1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp0000522(1,1,1)*ZZ(k,2,l,2)))
       temp0072222(2,1,1)=I34Z*(aux0072222(2,1,1)+4*F(4)*temp006222(2,2,
     &  1)+F(3)*temp72222(2,1,1)-det3*temp9222221(1,k,l)+8*temp0000522(2
     &  ,2,2)*ZZ(k,1,l,1)+40*temp0000522(2,2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)
     &  )+80*temp0000522(2,1,1)*ZZ(k,2,l,2)+temp006222(2,1,1)*(10*r10*(Z
     &  Z(k,1,l,2)+ZZ(k,2,l,1))+20*r21*ZZ(k,2,l,2)))
       temp0072222(2,2,1)=I34Z*(aux0072222(2,2,1)+2*F(4)*temp006222(2,2,
     &  2)+F(3)*temp72222(2,2,1)-det3*temp9222222(1,k,l)+24*temp0000522(
     &  2,2,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+120*temp0000522(2,2,1)*ZZ(k,2,l
     &  ,2)+temp006222(2,2,1)*(12*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+24*r21*Z
     &  Z(k,2,l,2)))
       temp0072222(2,2,2)=I34Z*(aux0072222(2,2,2)+F(3)*temp72222(2,2,2)-
     &  det3*temp9222222(2,k,l)+168*temp0000522(2,2,2)*ZZ(k,2,l,2)+temp0
     &  06222(2,2,2)*(14*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+28*r21*ZZ(k,2,l,2
     &  )))
       temp0071111(1,1,2)=temp0072111(1,1,1)
       temp0071111(1,2,1)=temp0072111(1,1,1)
       temp0071111(1,2,2)=temp0072211(1,1,1)
       temp0072111(1,1,2)=temp0072211(1,1,1)
       temp0072111(1,2,1)=temp0072211(1,1,1)
       temp0072111(1,2,2)=temp0072221(1,1,1)
       temp0072211(1,1,2)=temp0072221(1,1,1)
       temp0072211(1,2,1)=temp0072221(1,1,1)
       temp0072211(1,2,2)=temp0072222(1,1,1)
       temp0072221(1,1,2)=temp0072222(1,1,1)
       temp0072221(1,2,1)=temp0072222(1,1,1)
       temp0072221(1,2,2)=temp0072222(2,1,1)
       temp0072222(1,1,2)=temp0072222(2,1,1)
       temp0072222(1,2,1)=temp0072222(2,1,1)
       temp0072222(1,2,2)=temp0072222(2,2,1)
       temp0072222(2,1,2)=temp0072222(2,2,1)
       temp811111(1,1,1)=IX*(aux811111(1,1,1)+det3*temp9111111(1,1,jj)+1
     &  6*temp0071111(1,1,1)*Z(jj,1))
       temp821111(1,1,1)=IX*(aux821111(1,1,1)+det3*temp9211111(1,1,jj)+1
     &  4*temp0072111(1,1,1)*Z(jj,1)+2*temp0071111(1,1,1)*Z(jj,2))
       temp822111(1,1,1)=IX*(aux822111(1,1,1)+det3*temp9221111(1,1,jj)+1
     &  2*temp0072211(1,1,1)*Z(jj,1)+4*temp0072111(1,1,1)*Z(jj,2))
       temp822211(1,1,1)=IX*(aux822211(1,1,1)+det3*temp9222111(1,1,jj)+1
     &  0*temp0072221(1,1,1)*Z(jj,1)+6*temp0072211(1,1,1)*Z(jj,2))
       temp822221(1,1,1)=IX*(aux822221(1,1,1)+det3*temp9222211(1,1,jj)+8
     &  *(temp0072222(1,1,1)*Z(jj,1)+temp0072221(1,1,1)*Z(jj,2)))
       temp822222(1,1,1)=IX*(aux822222(1,1,1)+det3*temp9222221(1,1,jj)+6
     &  *temp0072222(2,1,1)*Z(jj,1)+10*temp0072222(1,1,1)*Z(jj,2))
       temp822222(2,1,1)=IX*(aux822222(2,1,1)+det3*temp9222222(1,1,jj)+4
     &  *temp0072222(2,2,1)*Z(jj,1)+12*temp0072222(2,1,1)*Z(jj,2))
       temp822222(2,2,1)=IX*(aux822222(2,2,1)+det3*temp9222222(2,1,jj)+2
     &  *temp0072222(2,2,2)*Z(jj,1)+14*temp0072222(2,2,1)*Z(jj,2))
       temp822222(2,2,2)=IX*(aux822222(2,2,2)+det3*temp9222222(2,2,jj)+1
     &  6*temp0072222(2,2,2)*Z(jj,2))
       temp811111(1,1,2)=temp821111(1,1,1)
       temp811111(1,2,1)=temp821111(1,1,1)
       temp811111(1,2,2)=temp822111(1,1,1)
       temp821111(1,1,2)=temp822111(1,1,1)
       temp821111(1,2,1)=temp822111(1,1,1)
       temp821111(1,2,2)=temp822211(1,1,1)
       temp822111(1,1,2)=temp822211(1,1,1)
       temp822111(1,2,1)=temp822211(1,1,1)
       temp822111(1,2,2)=temp822221(1,1,1)
       temp822211(1,1,2)=temp822221(1,1,1)
       temp822211(1,2,1)=temp822221(1,1,1)
       temp822211(1,2,2)=temp822222(1,1,1)
       temp822221(1,1,2)=temp822222(1,1,1)
       temp822221(1,2,1)=temp822222(1,1,1)
       temp822221(1,2,2)=temp822222(2,1,1)
       temp822222(1,1,2)=temp822222(2,1,1)
       temp822222(1,2,1)=temp822222(2,1,1)
       temp822222(1,2,2)=temp822222(2,2,1)
       temp822222(2,1,2)=temp822222(2,2,1)
c                Step4
       tempC300000000=I18Z*(auxC300000000+tempC3000000*F(3)-det3*temp000
     &  0002(k,l))
       temp0000002(1,1)=I22Z*(aux0000002(1,1)+4*F(4)*temp0000001(1)+F(3)
     &  *temp00002(1,1)-det3*temp000041(1,k,l)+8*tempC300000000*ZZ(k,1,l
     &  ,1))
       temp0000002(2,1)=I22Z*(aux0000002(2,1)+F(5)*temp0000001(1)+2*F(4)
     &  *temp0000001(2)+F(3)*temp00002(2,1)-det3*temp000042(1,k,l)+4*tem
     &  pC300000000*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp0000002(2,2)=I22Z*(aux0000002(2,2)+2*F(5)*temp0000001(2)+F(3)
     &  *temp00002(2,2)-det3*temp000042(2,k,l)+8*tempC300000000*ZZ(k,2,l
     &  ,2))
       temp0000002(1,2)=temp0000002(2,1)
       temp000041(1,1,1)=I26Z*(aux000041(1,1,1)+8*F(4)*temp00003(1,1,1)+
     &  F(3)*temp0041(1,1,1)-det3*temp006111(1,k,l)+48*temp0000002(1,1)*
     &  ZZ(k,1,l,1))
       temp000042(1,1,1)=I26Z*(aux000042(1,1,1)+F(5)*temp00003(1,1,1)+6*
     &  F(4)*temp00003(2,1,1)+F(3)*temp0042(1,1,1)-det3*temp006211(1,k,l
     &  )+24*temp0000002(2,1)*ZZ(k,1,l,1)+12*temp0000002(1,1)*(ZZ(k,1,l,
     &  2)+ZZ(k,2,l,1)))
       temp000042(2,1,1)=I26Z*(aux000042(2,1,1)+2*F(5)*temp00003(2,1,1)+
     &  4*F(4)*temp00003(2,2,1)+F(3)*temp0042(2,1,1)-det3*temp006221(1,k
     &  ,l)+16*temp0000002(2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+8*(temp0000002
     &  (2,2)*ZZ(k,1,l,1)+temp0000002(1,1)*ZZ(k,2,l,2)))
       temp000042(2,2,1)=I26Z*(aux000042(2,2,1)+2*F(4)*temp00003(2,2,2)+
     &  F(3)*temp0042(2,2,1)-det3*temp006222(1,k,l)+12*temp0000002(2,2)*
     &  (ZZ(k,1,l,2)+ZZ(k,2,l,1))+24*temp0000002(2,1)*ZZ(k,2,l,2)+temp00
     &  003(2,2,1)*(6*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2)))
       temp000042(2,2,2)=I26Z*(aux000042(2,2,2)+4*F(5)*temp00003(2,2,2)+
     &  F(3)*temp0042(2,2,2)-det3*temp006222(2,k,l)+48*temp0000002(2,2)*
     &  ZZ(k,2,l,2))
       temp000041(1,1,2)=temp000042(1,1,1)
       temp000041(1,2,1)=temp000042(1,1,1)
       temp000041(1,2,2)=temp000042(2,1,1)
       temp000042(1,1,2)=temp000042(2,1,1)
       temp000042(1,2,1)=temp000042(2,1,1)
       temp000042(1,2,2)=temp000042(2,2,1)
       temp000042(2,1,2)=temp000042(2,2,1)
       temp006111(1,1,1)=I30Z*(aux006111(1,1,1)+12*F(4)*temp00511(1,1,1)
     &  +F(3)*temp6111(1,1,1)-det3*temp811111(1,k,l)+120*temp000041(1,1,
     &  1)*ZZ(k,1,l,1))
       temp006211(1,1,1)=I30Z*(aux006211(1,1,1)+F(5)*temp00511(1,1,1)+10
     &  *F(4)*temp00521(1,1,1)+F(3)*temp6211(1,1,1)-det3*temp821111(1,k,
     &  l)+80*temp000042(1,1,1)*ZZ(k,1,l,1)+20*temp000041(1,1,1)*(ZZ(k,1
     &  ,l,2)+ZZ(k,2,l,1)))
       temp006221(1,1,1)=I30Z*(aux006221(1,1,1)+2*F(5)*temp00521(1,1,1)+
     &  F(3)*temp6221(1,1,1)-det3*temp822111(1,k,l)+48*temp000042(2,1,1)
     &  *ZZ(k,1,l,1)+32*temp000042(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+8*(F
     &  (4)*temp00522(1,1,1)+temp000041(1,1,1)*ZZ(k,2,l,2)))
       temp006222(1,1,1)=I30Z*(aux006222(1,1,1)+6*F(4)*temp00522(2,1,1)+
     &  F(3)*temp6222(1,1,1)-det3*temp822211(1,k,l)+36*temp000042(2,1,1)
     &  *(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp00522(1,1,1)*(6*r10*(ZZ(k,1,l,2)+
     &  ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2))+24*(temp000042(2,2,1)*ZZ(k,1,l,
     &  1)+temp000042(1,1,1)*ZZ(k,2,l,2)))
       temp006222(2,1,1)=I30Z*(aux006222(2,1,1)+4*F(5)*temp00522(2,1,1)+
     &  4*F(4)*temp00522(2,2,1)+F(3)*temp6222(2,1,1)-det3*temp822221(1,k
     &  ,l)+8*temp000042(2,2,2)*ZZ(k,1,l,1)+32*temp000042(2,2,1)*(ZZ(k,1
     &  ,l,2)+ZZ(k,2,l,1))+48*temp000042(2,1,1)*ZZ(k,2,l,2))
       temp006222(2,2,1)=I30Z*(aux006222(2,2,1)+2*F(4)*temp00522(2,2,2)+
     &  F(3)*temp6222(2,2,1)-det3*temp822222(1,k,l)+20*temp000042(2,2,2)
     &  *(ZZ(k,1,l,2)+ZZ(k,2,l,1))+80*temp000042(2,2,1)*ZZ(k,2,l,2)+temp
     &  00522(2,2,1)*(10*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+20*r21*ZZ(k,2,l,2
     &  )))
       temp006222(2,2,2)=I30Z*(aux006222(2,2,2)+F(3)*temp6222(2,2,2)-det
     &  3*temp822222(2,k,l)+120*temp000042(2,2,2)*ZZ(k,2,l,2)+temp00522(
     &  2,2,2)*(12*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+24*r21*ZZ(k,2,l,2)))
       temp006111(1,1,2)=temp006211(1,1,1)
       temp006111(1,2,1)=temp006211(1,1,1)
       temp006111(1,2,2)=temp006221(1,1,1)
       temp006211(1,1,2)=temp006221(1,1,1)
       temp006211(1,2,1)=temp006221(1,1,1)
       temp006211(1,2,2)=temp006222(1,1,1)
       temp006221(1,1,2)=temp006222(1,1,1)
       temp006221(1,2,1)=temp006222(1,1,1)
       temp006221(1,2,2)=temp006222(2,1,1)
       temp006222(1,1,2)=temp006222(2,1,1)
       temp006222(1,2,1)=temp006222(2,1,1)
       temp006222(1,2,2)=temp006222(2,2,1)
       temp006222(2,1,2)=temp006222(2,2,1)
       temp71111(1,1,1)=IX*(aux71111(1,1,1)+det3*temp811111(1,1,jj)+14*t
     &  emp006111(1,1,1)*Z(jj,1))
       temp72111(1,1,1)=IX*(aux72111(1,1,1)+det3*temp821111(1,1,jj)+12*t
     &  emp006211(1,1,1)*Z(jj,1)+2*temp006111(1,1,1)*Z(jj,2))
       temp72211(1,1,1)=IX*(aux72211(1,1,1)+det3*temp822111(1,1,jj)+10*t
     &  emp006221(1,1,1)*Z(jj,1)+4*temp006211(1,1,1)*Z(jj,2))
       temp72221(1,1,1)=IX*(aux72221(1,1,1)+det3*temp822211(1,1,jj)+8*te
     &  mp006222(1,1,1)*Z(jj,1)+6*temp006221(1,1,1)*Z(jj,2))
       temp72222(1,1,1)=IX*(aux72222(1,1,1)+det3*temp822221(1,1,jj)+6*te
     &  mp006222(2,1,1)*Z(jj,1)+8*temp006222(1,1,1)*Z(jj,2))
       temp72222(2,1,1)=IX*(aux72222(2,1,1)+det3*temp822222(1,1,jj)+4*te
     &  mp006222(2,2,1)*Z(jj,1)+10*temp006222(2,1,1)*Z(jj,2))
       temp72222(2,2,1)=IX*(aux72222(2,2,1)+det3*temp822222(2,1,jj)+2*te
     &  mp006222(2,2,2)*Z(jj,1)+12*temp006222(2,2,1)*Z(jj,2))
       temp72222(2,2,2)=IX*(aux72222(2,2,2)+det3*temp822222(2,2,jj)+14*t
     &  emp006222(2,2,2)*Z(jj,2))
       temp71111(1,1,2)=temp72111(1,1,1)
       temp71111(1,2,1)=temp72111(1,1,1)
       temp71111(1,2,2)=temp72211(1,1,1)
       temp72111(1,1,2)=temp72211(1,1,1)
       temp72111(1,2,1)=temp72211(1,1,1)
       temp72111(1,2,2)=temp72221(1,1,1)
       temp72211(1,1,2)=temp72221(1,1,1)
       temp72211(1,2,1)=temp72221(1,1,1)
       temp72211(1,2,2)=temp72222(1,1,1)
       temp72221(1,1,2)=temp72222(1,1,1)
       temp72221(1,2,1)=temp72222(1,1,1)
       temp72221(1,2,2)=temp72222(2,1,1)
       temp72222(1,1,2)=temp72222(2,1,1)
       temp72222(1,2,1)=temp72222(2,1,1)
       temp72222(1,2,2)=temp72222(2,2,1)
       temp72222(2,1,2)=temp72222(2,2,1)
c                Step5
       temp0000001(1)=I18Z*(aux0000001(1)+2*tempC3000000*F(4)+F(3)*temp0
     &  0001(1)-det3*temp00003(1,k,l))
       temp0000001(2)=I18Z*(aux0000001(2)+tempC3000000*F(5)+F(3)*temp000
     &  01(2)-det3*temp00003(2,k,l))
       temp00003(1,1,1)=I22Z*(aux00003(1,1,1)+6*F(4)*temp00002(1,1)+F(3)
     &  *temp003(1,1,1)-det3*temp00511(1,k,l)+24*temp0000001(1)*ZZ(k,1,l
     &  ,1))
       temp00003(2,1,1)=I22Z*(aux00003(2,1,1)+F(5)*temp00002(1,1)+4*F(4)
     &  *temp00002(2,1)+F(3)*temp003(2,1,1)-det3*temp00521(1,k,l)+8*(tem
     &  p0000001(2)*ZZ(k,1,l,1)+temp0000001(1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))
     &  ))
       temp00003(2,2,1)=I22Z*(aux00003(2,2,1)+2*F(5)*temp00002(2,1)+2*F(
     &  4)*temp00002(2,2)+F(3)*temp003(2,2,1)-det3*temp00522(1,k,l)+8*(t
     &  emp0000001(2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp0000001(1)*ZZ(k,2,l,
     &  2)))
       temp00003(2,2,2)=I22Z*(aux00003(2,2,2)+F(3)*temp003(2,2,2)-det3*t
     &  emp00522(2,k,l)+24*temp0000001(2)*ZZ(k,2,l,2)+temp00002(2,2)*(6*
     &  r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2)))
       temp00003(1,1,2)=temp00003(2,1,1)
       temp00003(1,2,1)=temp00003(2,1,1)
       temp00003(1,2,2)=temp00003(2,2,1)
       temp00003(2,1,2)=temp00003(2,2,1)
       temp00511(1,1,1)=I26Z*(aux00511(1,1,1)+10*F(4)*temp0041(1,1,1)+F(
     &  3)*temp511(1,1,1)-det3*temp71111(1,k,l)+80*temp00003(1,1,1)*ZZ(k
     &  ,1,l,1))
       temp00521(1,1,1)=I26Z*(aux00521(1,1,1)+F(5)*temp0041(1,1,1)+8*F(4
     &  )*temp0042(1,1,1)+F(3)*temp521(1,1,1)-det3*temp72111(1,k,l)+48*t
     &  emp00003(2,1,1)*ZZ(k,1,l,1)+16*temp00003(1,1,1)*(ZZ(k,1,l,2)+ZZ(
     &  k,2,l,1)))
       temp00522(1,1,1)=I26Z*(aux00522(1,1,1)+2*F(5)*temp0042(1,1,1)+6*F
     &  (4)*temp0042(2,1,1)+F(3)*temp522(1,1,1)-det3*temp72211(1,k,l)+24
     &  *(temp00003(2,2,1)*ZZ(k,1,l,1)+temp00003(2,1,1)*(ZZ(k,1,l,2)+ZZ(
     &  k,2,l,1)))+8*temp00003(1,1,1)*ZZ(k,2,l,2))
       temp00522(2,1,1)=I26Z*(aux00522(2,1,1)+4*F(4)*temp0042(2,2,1)+F(3
     &  )*temp522(2,1,1)-det3*temp72221(1,k,l)+8*temp00003(2,2,2)*ZZ(k,1
     &  ,l,1)+temp0042(2,1,1)*(6*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+12*r21*ZZ
     &  (k,2,l,2))+24*(temp00003(2,2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp00
     &  003(2,1,1)*ZZ(k,2,l,2)))
       temp00522(2,2,1)=I26Z*(aux00522(2,2,1)+4*F(5)*temp0042(2,2,1)+2*F
     &  (4)*temp0042(2,2,2)+F(3)*temp522(2,2,1)-det3*temp72222(1,k,l)+16
     &  *temp00003(2,2,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+48*temp00003(2,2,1)*
     &  ZZ(k,2,l,2))
       temp00522(2,2,2)=I26Z*(aux00522(2,2,2)+F(3)*temp522(2,2,2)-det3*t
     &  emp72222(2,k,l)+80*temp00003(2,2,2)*ZZ(k,2,l,2)+temp0042(2,2,2)*
     &  (10*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+20*r21*ZZ(k,2,l,2)))
       temp00511(1,1,2)=temp00521(1,1,1)
       temp00511(1,2,1)=temp00521(1,1,1)
       temp00511(1,2,2)=temp00522(1,1,1)
       temp00521(1,1,2)=temp00522(1,1,1)
       temp00521(1,2,1)=temp00522(1,1,1)
       temp00521(1,2,2)=temp00522(2,1,1)
       temp00522(1,1,2)=temp00522(2,1,1)
       temp00522(1,2,1)=temp00522(2,1,1)
       temp00522(1,2,2)=temp00522(2,2,1)
       temp00522(2,1,2)=temp00522(2,2,1)
       temp6111(1,1,1)=IX*(aux6111(1,1,1)+det3*temp71111(1,1,jj)+12*temp
     &  00511(1,1,1)*Z(jj,1))
       temp6211(1,1,1)=IX*(aux6211(1,1,1)+det3*temp72111(1,1,jj)+10*temp
     &  00521(1,1,1)*Z(jj,1)+2*temp00511(1,1,1)*Z(jj,2))
       temp6221(1,1,1)=IX*(aux6221(1,1,1)+det3*temp72211(1,1,jj)+8*temp0
     &  0522(1,1,1)*Z(jj,1)+4*temp00521(1,1,1)*Z(jj,2))
       temp6222(1,1,1)=IX*(aux6222(1,1,1)+det3*temp72221(1,1,jj)+6*(temp
     &  00522(2,1,1)*Z(jj,1)+temp00522(1,1,1)*Z(jj,2)))
       temp6222(2,1,1)=IX*(aux6222(2,1,1)+det3*temp72222(1,1,jj)+4*temp0
     &  0522(2,2,1)*Z(jj,1)+8*temp00522(2,1,1)*Z(jj,2))
       temp6222(2,2,1)=IX*(aux6222(2,2,1)+det3*temp72222(2,1,jj)+2*temp0
     &  0522(2,2,2)*Z(jj,1)+10*temp00522(2,2,1)*Z(jj,2))
       temp6222(2,2,2)=IX*(aux6222(2,2,2)+det3*temp72222(2,2,jj)+12*temp
     &  00522(2,2,2)*Z(jj,2))
       temp6111(1,1,2)=temp6211(1,1,1)
       temp6111(1,2,1)=temp6211(1,1,1)
       temp6111(1,2,2)=temp6221(1,1,1)
       temp6211(1,1,2)=temp6221(1,1,1)
       temp6211(1,2,1)=temp6221(1,1,1)
       temp6211(1,2,2)=temp6222(1,1,1)
       temp6221(1,1,2)=temp6222(1,1,1)
       temp6221(1,2,1)=temp6222(1,1,1)
       temp6221(1,2,2)=temp6222(2,1,1)
       temp6222(1,1,2)=temp6222(2,1,1)
       temp6222(1,2,1)=temp6222(2,1,1)
       temp6222(1,2,2)=temp6222(2,2,1)
       temp6222(2,1,2)=temp6222(2,2,1)
c                Step6
       tempC3000000=I14Z*(auxC3000000+tempC30000*F(3)-det3*temp00002(k,l
     &  ))
       temp00002(1,1)=I18Z*(aux00002(1,1)+4*F(4)*temp00001(1)+F(3)*temp0
     &  02(1,1)-det3*temp0041(1,k,l)+8*tempC3000000*ZZ(k,1,l,1))
       temp00002(2,1)=I18Z*(aux00002(2,1)+F(5)*temp00001(1)+2*F(4)*temp0
     &  0001(2)+F(3)*temp002(2,1)-det3*temp0042(1,k,l)+4*tempC3000000*(Z
     &  Z(k,1,l,2)+ZZ(k,2,l,1)))
       temp00002(2,2)=I18Z*(aux00002(2,2)+2*F(5)*temp00001(2)+F(3)*temp0
     &  02(2,2)-det3*temp0042(2,k,l)+8*tempC3000000*ZZ(k,2,l,2))
       temp00002(1,2)=temp00002(2,1)
       temp0041(1,1,1)=I22Z*(aux0041(1,1,1)+8*F(4)*temp003(1,1,1)+F(3)*t
     &  emp41(1,1,1)-det3*temp6111(1,k,l)+48*temp00002(1,1)*ZZ(k,1,l,1))
       temp0042(1,1,1)=I22Z*(aux0042(1,1,1)+F(5)*temp003(1,1,1)+6*F(4)*t
     &  emp003(2,1,1)+F(3)*temp42(1,1,1)-det3*temp6211(1,k,l)+24*temp000
     &  02(2,1)*ZZ(k,1,l,1)+12*temp00002(1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp0042(2,1,1)=I22Z*(aux0042(2,1,1)+2*F(5)*temp003(2,1,1)+4*F(4)
     &  *temp003(2,2,1)+F(3)*temp42(2,1,1)-det3*temp6221(1,k,l)+16*temp0
     &  0002(2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+8*(temp00002(2,2)*ZZ(k,1,l,1
     &  )+temp00002(1,1)*ZZ(k,2,l,2)))
       temp0042(2,2,1)=I22Z*(aux0042(2,2,1)+2*F(4)*temp003(2,2,2)+F(3)*t
     &  emp42(2,2,1)-det3*temp6222(1,k,l)+12*temp00002(2,2)*(ZZ(k,1,l,2)
     &  +ZZ(k,2,l,1))+24*temp00002(2,1)*ZZ(k,2,l,2)+temp003(2,2,1)*(6*r1
     &  0*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2)))
       temp0042(2,2,2)=I22Z*(aux0042(2,2,2)+4*F(5)*temp003(2,2,2)+F(3)*t
     &  emp42(2,2,2)-det3*temp6222(2,k,l)+48*temp00002(2,2)*ZZ(k,2,l,2))
       temp0041(1,1,2)=temp0042(1,1,1)
       temp0041(1,2,1)=temp0042(1,1,1)
       temp0041(1,2,2)=temp0042(2,1,1)
       temp0042(1,1,2)=temp0042(2,1,1)
       temp0042(1,2,1)=temp0042(2,1,1)
       temp0042(1,2,2)=temp0042(2,2,1)
       temp0042(2,1,2)=temp0042(2,2,1)
       temp511(1,1,1)=IX*(aux511(1,1,1)+det3*temp6111(1,1,jj)+10*temp004
     &  1(1,1,1)*Z(jj,1))
       temp521(1,1,1)=IX*(aux521(1,1,1)+det3*temp6211(1,1,jj)+8*temp0042
     &  (1,1,1)*Z(jj,1)+2*temp0041(1,1,1)*Z(jj,2))
       temp522(1,1,1)=IX*(aux522(1,1,1)+det3*temp6221(1,1,jj)+6*temp0042
     &  (2,1,1)*Z(jj,1)+4*temp0042(1,1,1)*Z(jj,2))
       temp522(2,1,1)=IX*(aux522(2,1,1)+det3*temp6222(1,1,jj)+4*temp0042
     &  (2,2,1)*Z(jj,1)+6*temp0042(2,1,1)*Z(jj,2))
       temp522(2,2,1)=IX*(aux522(2,2,1)+det3*temp6222(2,1,jj)+2*temp0042
     &  (2,2,2)*Z(jj,1)+8*temp0042(2,2,1)*Z(jj,2))
       temp522(2,2,2)=IX*(aux522(2,2,2)+det3*temp6222(2,2,jj)+10*temp004
     &  2(2,2,2)*Z(jj,2))
       temp511(1,1,2)=temp521(1,1,1)
       temp511(1,2,1)=temp521(1,1,1)
       temp511(1,2,2)=temp522(1,1,1)
       temp521(1,1,2)=temp522(1,1,1)
       temp521(1,2,1)=temp522(1,1,1)
       temp521(1,2,2)=temp522(2,1,1)
       temp522(1,1,2)=temp522(2,1,1)
       temp522(1,2,1)=temp522(2,1,1)
       temp522(1,2,2)=temp522(2,2,1)
       temp522(2,1,2)=temp522(2,2,1)
c                Step7
       temp00001(1)=I14Z*(aux00001(1)+2*tempC30000*F(4)+F(3)*temp001(1)-
     &  det3*temp003(1,k,l))
       temp00001(2)=I14Z*(aux00001(2)+tempC30000*F(5)+F(3)*temp001(2)-de
     &  t3*temp003(2,k,l))
       temp003(1,1,1)=I18Z*(aux003(1,1,1)+6*F(4)*temp002(1,1)+F(3)*temp3
     &  (1,1,1)-det3*temp511(1,k,l)+24*temp00001(1)*ZZ(k,1,l,1))
       temp003(2,1,1)=I18Z*(aux003(2,1,1)+F(5)*temp002(1,1)+4*F(4)*temp0
     &  02(2,1)+F(3)*temp3(2,1,1)-det3*temp521(1,k,l)+8*(temp00001(2)*ZZ
     &  (k,1,l,1)+temp00001(1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))))
       temp003(2,2,1)=I18Z*(aux003(2,2,1)+2*F(5)*temp002(2,1)+2*F(4)*tem
     &  p002(2,2)+F(3)*temp3(2,2,1)-det3*temp522(1,k,l)+8*(temp00001(2)*
     &  (ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp00001(1)*ZZ(k,2,l,2)))
       temp003(2,2,2)=I18Z*(aux003(2,2,2)+F(3)*temp3(2,2,2)-det3*temp522
     &  (2,k,l)+24*temp00001(2)*ZZ(k,2,l,2)+temp002(2,2)*(6*r10*(ZZ(k,1,
     &  l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2)))
       temp003(1,1,2)=temp003(2,1,1)
       temp003(1,2,1)=temp003(2,1,1)
       temp003(1,2,2)=temp003(2,2,1)
       temp003(2,1,2)=temp003(2,2,1)
       temp41(1,1,1)=IX*(aux41(1,1,1)+det3*temp511(1,1,jj)+8*temp003(1,1
     &  ,1)*Z(jj,1))
       temp42(1,1,1)=IX*(aux42(1,1,1)+det3*temp521(1,1,jj)+6*temp003(2,1
     &  ,1)*Z(jj,1)+2*temp003(1,1,1)*Z(jj,2))
       temp42(2,1,1)=IX*(aux42(2,1,1)+det3*temp522(1,1,jj)+4*(temp003(2,
     &  2,1)*Z(jj,1)+temp003(2,1,1)*Z(jj,2)))
       temp42(2,2,1)=IX*(aux42(2,2,1)+det3*temp522(2,1,jj)+2*temp003(2,2
     &  ,2)*Z(jj,1)+6*temp003(2,2,1)*Z(jj,2))
       temp42(2,2,2)=IX*(aux42(2,2,2)+det3*temp522(2,2,jj)+8*temp003(2,2
     &  ,2)*Z(jj,2))
       temp41(1,1,2)=temp42(1,1,1)
       temp41(1,2,1)=temp42(1,1,1)
       temp41(1,2,2)=temp42(2,1,1)
       temp42(1,1,2)=temp42(2,1,1)
       temp42(1,2,1)=temp42(2,1,1)
       temp42(1,2,2)=temp42(2,2,1)
       temp42(2,1,2)=temp42(2,2,1)
c                Step8
       tempC30000=I10Z*(auxC30000+tempC300*F(3)-det3*temp002(k,l))
       temp002(1,1)=I14Z*(aux002(1,1)+4*F(4)*temp001(1)+F(3)*temp2(1,1)-
     &  det3*temp41(1,k,l)+8*tempC30000*ZZ(k,1,l,1))
       temp002(2,1)=I14Z*(aux002(2,1)+F(5)*temp001(1)+2*F(4)*temp001(2)+
     &  F(3)*temp2(2,1)-det3*temp42(1,k,l)+4*tempC30000*(ZZ(k,1,l,2)+ZZ(
     &  k,2,l,1)))
       temp002(2,2)=I14Z*(aux002(2,2)+2*F(5)*temp001(2)+F(3)*temp2(2,2)-
     &  det3*temp42(2,k,l)+8*tempC30000*ZZ(k,2,l,2))
       temp002(1,2)=temp002(2,1)
       temp3(1,1,1)=IX*(aux3(1,1,1)+det3*temp41(1,1,jj)+6*temp002(1,1)*Z
     &  (jj,1))
       temp3(2,1,1)=IX*(aux3(2,1,1)+det3*temp42(1,1,jj)+4*temp002(2,1)*Z
     &  (jj,1)+2*temp002(1,1)*Z(jj,2))
       temp3(2,2,1)=IX*(aux3(2,2,1)+det3*temp42(2,1,jj)+2*temp002(2,2)*Z
     &  (jj,1)+4*temp002(2,1)*Z(jj,2))
       temp3(2,2,2)=IX*(aux3(2,2,2)+det3*temp42(2,2,jj)+6*temp002(2,2)*Z
     &  (jj,2))
       temp3(1,1,2)=temp3(2,1,1)
       temp3(1,2,1)=temp3(2,1,1)
       temp3(1,2,2)=temp3(2,2,1)
       temp3(2,1,2)=temp3(2,2,1)
c                Step9
       temp001(1)=I10Z*(aux001(1)+2*tempC300*F(4)+F(3)*temp1(1)-det3*tem
     &  p3(1,k,l))
       temp001(2)=I10Z*(aux001(2)+tempC300*F(5)+F(3)*temp1(2)-det3*temp3
     &  (2,k,l))
       temp2(1,1)=IX*(aux2(1,1)+det3*temp3(1,1,jj)+4*temp001(1)*Z(jj,1))
       temp2(2,1)=IX*(aux2(2,1)+det3*temp3(2,1,jj)+2*(temp001(2)*Z(jj,1)
     &  +temp001(1)*Z(jj,2)))
       temp2(2,2)=IX*(aux2(2,2)+det3*temp3(2,2,jj)+4*temp001(2)*Z(jj,2))
       temp2(1,2)=temp2(2,1)
c                Step10
       tempC300=I6Z*(auxC300+tempC30*F(3)-det3*temp2(k,l))
       temp1(1)=IX*(aux1(1)+det3*temp2(1,jj)+2*tempC300*Z(jj,1))
       temp1(2)=IX*(aux1(2)+det3*temp2(2,jj)+2*tempC300*Z(jj,2))
c                Step11
       tempC30=IX*(auxC30+det3*temp1(jj))

         ac=4
         accuracyCR(1,0,ac) = abs(C30     /tempC30      -1d0) 
         accuracyCR(1,1,ac) = abs(Cij(1,1)/temp1(1)      -1d0) 
         accuracyCR(2,1,ac) = abs(Cij(2,1)/temp1(2)      -1d0)
         accuracyCR(1,2,ac) = abs(Cij(1,2)/temp2(1,1)    -1d0)
         accuracyCR(2,2,ac) = abs(Cij(2,2)/temp2(2,2)    -1d0)
         accuracyCR(3,2,ac) = abs(Cij(3,2)/temp2(2,1)    -1d0)
         accuracyCR(4,2,ac) = abs(Cij(4,2)/tempC300      -1d0)
         accuracyCR(1,3,ac) = abs(Cij(1,3)/temp3(1,1,1)  -1d0)
         accuracyCR(2,3,ac) = abs(Cij(2,3)/temp3(2,2,2)  -1d0)
         accuracyCR(3,3,ac) = abs(Cij(3,3)/temp3(2,1,1)  -1d0)
         accuracyCR(4,3,ac) = abs(Cij(4,3)/temp3(2,2,1)  -1d0)
         accuracyCR(5,3,ac) = abs(Cij(5,3)/temp001(1)    -1d0)
         accuracyCR(6,3,ac) = abs(Cij(6,3)/temp001(2)    -1d0)
         accuracyCR(1,4,ac) = abs(Cij(1,4)/temp41(1,1,1) -1d0)
         accuracyCR(2,4,ac) = abs(Cij(2,4)/temp42(2,2,2) -1d0)
         accuracyCR(3,4,ac) = abs(Cij(3,4)/temp42(1,1,1) -1d0)
         accuracyCR(4,4,ac) = abs(Cij(4,4)/temp42(2,1,1) -1d0)
         accuracyCR(5,4,ac) = abs(Cij(5,4)/temp42(2,2,1) -1d0)
         accuracyCR(6,4,ac) = abs(Cij(6,4)/temp002(1,1)  -1d0)
         accuracyCR(7,4,ac) = abs(Cij(7,4)/temp002(2,2)  -1d0)
         accuracyCR(8,4,ac) = abs(Cij(8,4)/temp002(2,1)  -1d0)
         accuracyCR(9,4,ac) = abs(Cij(9,4)/tempC30000    -1d0)


      DO I1=0,4
           accuracyC(i1,ac)=accuracyCR(1,i1,ac)
            DO I2=2,INDEX(I1)  
       if(accuracyCR(i2,i1,ac).gt.accuracyC(i1,ac)) then
          accuracyC(i1,ac)=accuracyCR(i2,i1,ac)
       endif
          enddo
        enddo

c          if(accuracyC(4,ac).lt.1d-16) goto 500

           C30=tempC30
           Cij(1,1)=temp1(1)
           Cij(2,1)=temp1(2)
           Cij(1,2)=temp2(1,1)
           Cij(2,2)=temp2(2,2)
           Cij(3,2)=temp2(2,1)
           Cij(4,2)=tempC300
           Cij(1,3)=temp3(1,1,1)
           Cij(2,3)=temp3(2,2,2)
           Cij(3,3)=temp3(2,1,1)
           Cij(4,3)=temp3(2,2,1)
           Cij(5,3)=temp001(1)
           Cij(6,3)=temp001(2)
           Cij(1,4)=temp41(1,1,1)
           Cij(2,4)=temp42(2,2,2)
           Cij(3,4)=temp42(1,1,1)
           Cij(4,4)=temp42(2,1,1)
           Cij(5,4)=temp42(2,2,1)
           Cij(6,4)=temp002(1,1)
           Cij(7,4)=temp002(2,2)
           Cij(8,4)=temp002(2,1)
           Cij(9,4)=tempC30000

       if(order.eq.10) goto 500
c                Iteration11
c                Step1
       S2000000000000=2*B23(6,10)
       S2000000000021(1)=2*B23(5,8)
       S2000000000021(2)=-2*B23(5,9)
       S2000000000022(1)=-2*B23(5,9)
       S2000000000022(2)=2*B23(5,10)
       S2h000000000021(1)=B13(6,11)+B23(6,10)
       S2h000000000021(2)=B13(6,11)-B23(6,11)
       S2h000000000022(1)=B12(6,11)-B13(6,11)
       S2h000000000022(2)=-B13(6,11)
       auxC3000000000000=-(F(1)*S2h00000000001(1))-F(2)*S2h00000000001(2
     &  )+S2h000000000021(k)*Z(1,l)+S2h000000000022(k)*Z(2,l)+(Inv120+S2
     &  000000000000-S2h000000000021(1)-S2h000000000022(2))*Z(k,l)
       tempC3000000000000=I26Z*(auxC3000000000000+tempC30000000000*F(3))
       S2h000000004111(1)=B13(5,11)+B23(5,8)
       S2h000000004111(2)=B13(5,11)-B23(5,9)
       S2h000000004121(1)=B13(5,11)-B23(5,9)
       S2h000000004121(2)=B13(5,11)+B23(5,10)
       S2h000000004122(1)=B13(5,11)+B23(5,10)
       S2h000000004122(2)=B13(5,11)-B23(5,11)
       S2h000000004211(1)=B12(5,11)-B13(5,11)
       S2h000000004211(2)=-B13(5,11)
       S2h000000004221(1)=-B13(5,11)
       S2h000000004221(2)=-B13(5,11)
       S2h000000004222(1)=-B13(5,11)
       S2h000000004222(2)=-B13(5,11)
       aux00000000002(1,1)=-(F(1)*S2h00000000311(1))-F(2)*S2h00000000321
     &  (1)+S2h000000004111(k)*Z(1,l)+S2h000000004211(k)*Z(2,l)+(Inv1001
     &  1+S2000000000021(1)-S2h000000004111(1)-S2h000000004221(1))*Z(k,l
     &  )-4*(S2h000000000021(1)*ZZ(k,1,l,1)+S2h000000000022(1)*ZZ(k,1,l,
     &  2))
       aux00000000002(2,1)=-(F(1)*S2h00000000312(1))-F(2)*S2h00000000322
     &  (1)+S2h000000004121(k)*Z(1,l)+S2h000000004221(k)*Z(2,l)+(Inv1002
     &  1+S2000000000022(1)-S2h000000004121(1)-S2h000000004222(1))*Z(k,l
     &  )-2*(S2h000000000021(2)*ZZ(k,1,l,1)+S2h000000000022(2)*ZZ(k,1,l,
     &  2))-2*S2h000000000021(1)*ZZ(k,2,l,1)-2*S2h000000000022(1)*ZZ(k,2
     &  ,l,2)
       aux00000000002(2,2)=-(F(1)*S2h00000000312(2))-F(2)*S2h00000000322
     &  (2)+S2h000000004122(k)*Z(1,l)+S2h000000004222(k)*Z(2,l)+(Inv1002
     &  2+S2000000000022(2)-S2h000000004122(1)-S2h000000004222(2))*Z(k,l
     &  )-4*(S2h000000000021(2)*ZZ(k,2,l,1)+S2h000000000022(2)*ZZ(k,2,l,
     &  2))
       temp00000000002(1,1)=I30Z*(aux00000000002(1,1)+4*F(4)*temp0000000
     &  0001(1)+F(3)*temp000000002(1,1)+8*tempC3000000000000*ZZ(k,1,l,1)
     &  )
       temp00000000002(2,1)=I30Z*(aux00000000002(2,1)+F(5)*temp000000000
     &  01(1)+2*F(4)*temp00000000001(2)+F(3)*temp000000002(2,1)+4*tempC3
     &  000000000000*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp00000000002(2,2)=I30Z*(aux00000000002(2,2)+2*F(5)*temp0000000
     &  0001(2)+F(3)*temp000000002(2,2)+8*tempC3000000000000*ZZ(k,2,l,2)
     &  )
       temp00000000002(1,2)=temp00000000002(2,1)
       S2000000004111(1)=2*B23(4,6)
       S2000000004111(2)=-2*B23(4,7)
       S2000000004121(1)=-2*B23(4,7)
       S2000000004121(2)=2*B23(4,8)
       S2000000004122(1)=2*B23(4,8)
       S2000000004122(2)=-2*B23(4,9)
       S2000000004211(1)=-2*B23(4,7)
       S2000000004211(2)=2*B23(4,8)
       S2000000004221(1)=2*B23(4,8)
       S2000000004221(2)=-2*B23(4,9)
       S2000000004222(1)=-2*B23(4,9)
       S2000000004222(2)=2*B23(4,10)
       S2000000611111(1)=2*B23(3,4)
       S2000000611111(2)=-2*B23(3,5)
       S2000000612111(1)=-2*B23(3,5)
       S2000000612111(2)=2*B23(3,6)
       S2000000612211(1)=2*B23(3,6)
       S2000000612211(2)=-2*B23(3,7)
       S2000000612221(1)=-2*B23(3,7)
       S2000000612221(2)=2*B23(3,8)
       S2000000612222(1)=2*B23(3,8)
       S2000000612222(2)=-2*B23(3,9)
       S2000000621111(1)=-2*B23(3,5)
       S2000000621111(2)=2*B23(3,6)
       S2000000622111(1)=2*B23(3,6)
       S2000000622111(2)=-2*B23(3,7)
       S2000000622211(1)=-2*B23(3,7)
       S2000000622211(2)=2*B23(3,8)
       S2000000622221(1)=2*B23(3,8)
       S2000000622221(2)=-2*B23(3,9)
       S2000000622222(1)=-2*B23(3,9)
       S2000000622222(2)=2*B23(3,10)
       S2h000000611111(1)=B13(4,11)+B23(4,6)
       S2h000000611111(2)=B13(4,11)-B23(4,7)
       S2h000000612111(1)=B13(4,11)-B23(4,7)
       S2h000000612111(2)=B13(4,11)+B23(4,8)
       S2h000000612211(1)=B13(4,11)+B23(4,8)
       S2h000000612211(2)=B13(4,11)-B23(4,9)
       S2h000000612221(1)=B13(4,11)-B23(4,9)
       S2h000000612221(2)=B13(4,11)+B23(4,10)
       S2h000000612222(1)=B13(4,11)+B23(4,10)
       S2h000000612222(2)=B13(4,11)-B23(4,11)
       S2h000000621111(1)=B12(4,11)-B13(4,11)
       S2h000000621111(2)=-B13(4,11)
       S2h000000622111(1)=-B13(4,11)
       S2h000000622111(2)=-B13(4,11)
       S2h000000622211(1)=-B13(4,11)
       S2h000000622211(2)=-B13(4,11)
       S2h000000622221(1)=-B13(4,11)
       S2h000000622221(2)=-B13(4,11)
       S2h000000622222(1)=-B13(4,11)
       S2h000000622222(2)=-B13(4,11)
       aux0000000041(1,1,1)=-(F(1)*S2h00000051111(1))-F(2)*S2h0000005211
     &  1(1)+S2h000000611111(k)*Z(1,l)+S2h000000621111(k)*Z(2,l)+(Inv801
     &  111+S2000000004111(1)-S2h000000611111(1)-S2h000000622111(1))*Z(k
     &  ,l)-8*(S2h000000004111(1)*ZZ(k,1,l,1)+S2h000000004211(1)*ZZ(k,1,
     &  l,2))
       aux0000000042(1,1,1)=-(F(1)*S2h00000051211(1))-F(2)*S2h0000005221
     &  1(1)+S2h000000612111(k)*Z(1,l)+S2h000000622111(k)*Z(2,l)+(Inv802
     &  111+S2000000004211(1)-S2h000000612111(1)-S2h000000622211(1))*Z(k
     &  ,l)-2*(3*(S2h000000004121(1)*ZZ(k,1,l,1)+S2h000000004221(1)*ZZ(k
     &  ,1,l,2))+S2h000000004111(1)*ZZ(k,2,l,1)+S2h000000004211(1)*ZZ(k,
     &  2,l,2))
       aux0000000042(2,1,1)=-(F(1)*S2h00000051221(1))-F(2)*S2h0000005222
     &  1(1)+S2h000000612211(k)*Z(1,l)+S2h000000622211(k)*Z(2,l)+(Inv802
     &  211+S2000000004221(1)-S2h000000612211(1)-S2h000000622221(1))*Z(k
     &  ,l)-4*(S2h000000004122(1)*ZZ(k,1,l,1)+S2h000000004222(1)*ZZ(k,1,
     &  l,2)+S2h000000004121(1)*ZZ(k,2,l,1)+S2h000000004221(1)*ZZ(k,2,l,
     &  2))
       aux0000000042(2,2,1)=-(F(1)*S2h00000051222(1))-F(2)*S2h0000005222
     &  2(1)+S2h000000612221(k)*Z(1,l)+S2h000000622221(k)*Z(2,l)+(Inv802
     &  221+S2000000004222(1)-S2h000000612221(1)-S2h000000622222(1))*Z(k
     &  ,l)-2*(S2h000000004122(2)*ZZ(k,1,l,1)+S2h000000004222(2)*ZZ(k,1,
     &  l,2)+3*(S2h000000004122(1)*ZZ(k,2,l,1)+S2h000000004222(1)*ZZ(k,2
     &  ,l,2)))
       aux0000000042(2,2,2)=-(F(1)*S2h00000051222(2))-F(2)*S2h0000005222
     &  2(2)+S2h000000612222(k)*Z(1,l)+S2h000000622222(k)*Z(2,l)+(Inv802
     &  222+S2000000004222(2)-S2h000000612222(1)-S2h000000622222(2))*Z(k
     &  ,l)-8*(S2h000000004122(2)*ZZ(k,2,l,1)+S2h000000004222(2)*ZZ(k,2,
     &  l,2))
       temp0000000041(1,1,1)=I34Z*(aux0000000041(1,1,1)+8*F(4)*temp00000
     &  0003(1,1,1)+F(3)*temp00000041(1,1,1)+48*temp00000000002(1,1)*ZZ(
     &  k,1,l,1))
       temp0000000042(1,1,1)=I34Z*(aux0000000042(1,1,1)+F(5)*temp0000000
     &  03(1,1,1)+6*F(4)*temp000000003(2,1,1)+F(3)*temp00000042(1,1,1)+2
     &  4*temp00000000002(2,1)*ZZ(k,1,l,1)+12*temp00000000002(1,1)*(ZZ(k
     &  ,1,l,2)+ZZ(k,2,l,1)))
       temp0000000042(2,1,1)=I34Z*(aux0000000042(2,1,1)+2*F(5)*temp00000
     &  0003(2,1,1)+4*F(4)*temp000000003(2,2,1)+F(3)*temp00000042(2,1,1)
     &  +16*temp00000000002(2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+8*(temp000000
     &  00002(2,2)*ZZ(k,1,l,1)+temp00000000002(1,1)*ZZ(k,2,l,2)))
       temp0000000042(2,2,1)=I34Z*(aux0000000042(2,2,1)+2*F(4)*temp00000
     &  0003(2,2,2)+F(3)*temp00000042(2,2,1)+12*temp00000000002(2,2)*(ZZ
     &  (k,1,l,2)+ZZ(k,2,l,1))+24*temp00000000002(2,1)*ZZ(k,2,l,2)+temp0
     &  00000003(2,2,1)*(6*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l
     &  ,2)))
       temp0000000042(2,2,2)=I34Z*(aux0000000042(2,2,2)+4*F(5)*temp00000
     &  0003(2,2,2)+F(3)*temp00000042(2,2,2)+48*temp00000000002(2,2)*ZZ(
     &  k,2,l,2))
       temp0000000041(1,1,2)=temp0000000042(1,1,1)
       temp0000000041(1,2,1)=temp0000000042(1,1,1)
       temp0000000041(1,2,2)=temp0000000042(2,1,1)
       temp0000000042(1,1,2)=temp0000000042(2,1,1)
       temp0000000042(1,2,1)=temp0000000042(2,1,1)
       temp0000000042(1,2,2)=temp0000000042(2,2,1)
       temp0000000042(2,1,2)=temp0000000042(2,2,1)
       S2h000081111111(1)=B13(3,11)+B23(3,4)
       S2h000081111111(2)=B13(3,11)-B23(3,5)
       S2h000081211111(1)=B13(3,11)-B23(3,5)
       S2h000081211111(2)=B13(3,11)+B23(3,6)
       S2h000081221111(1)=B13(3,11)+B23(3,6)
       S2h000081221111(2)=B13(3,11)-B23(3,7)
       S2h000081222111(1)=B13(3,11)-B23(3,7)
       S2h000081222111(2)=B13(3,11)+B23(3,8)
       S2h000081222211(1)=B13(3,11)+B23(3,8)
       S2h000081222211(2)=B13(3,11)-B23(3,9)
       S2h000081222221(1)=B13(3,11)-B23(3,9)
       S2h000081222221(2)=B13(3,11)+B23(3,10)
       S2h000081222222(1)=B13(3,11)+B23(3,10)
       S2h000081222222(2)=B13(3,11)-B23(3,11)
       S2h000082111111(1)=B12(3,11)-B13(3,11)
       S2h000082111111(2)=-B13(3,11)
       S2h000082211111(1)=-B13(3,11)
       S2h000082211111(2)=-B13(3,11)
       S2h000082221111(1)=-B13(3,11)
       S2h000082221111(2)=-B13(3,11)
       S2h000082222111(1)=-B13(3,11)
       S2h000082222111(2)=-B13(3,11)
       S2h000082222211(1)=-B13(3,11)
       S2h000082222211(2)=-B13(3,11)
       S2h000082222221(1)=-B13(3,11)
       S2h000082222221(2)=-B13(3,11)
       S2h000082222222(1)=-B13(3,11)
       S2h000082222222(2)=-B13(3,11)
       aux0000006111(1,1,1)=-(F(1)*S2h00007111111(1))-F(2)*S2h0000721111
     &  1(1)+S2h000081111111(k)*Z(1,l)+S2h000082111111(k)*Z(2,l)+(Inv601
     &  11111+S2000000611111(1)-S2h000081111111(1)-S2h000082211111(1))*Z
     &  (k,l)-12*(S2h000000611111(1)*ZZ(k,1,l,1)+S2h000000621111(1)*ZZ(k
     &  ,1,l,2))
       aux0000006211(1,1,1)=-(F(1)*S2h00007121111(1))-F(2)*S2h0000722111
     &  1(1)+S2h000081211111(k)*Z(1,l)+S2h000082211111(k)*Z(2,l)+(Inv602
     &  11111+S2000000621111(1)-S2h000081211111(1)-S2h000082221111(1))*Z
     &  (k,l)-2*(5*(S2h000000612111(1)*ZZ(k,1,l,1)+S2h000000622111(1)*ZZ
     &  (k,1,l,2))+S2h000000611111(1)*ZZ(k,2,l,1)+S2h000000621111(1)*ZZ(
     &  k,2,l,2))
       aux0000006221(1,1,1)=-(F(1)*S2h00007122111(1))-F(2)*S2h0000722211
     &  1(1)+S2h000081221111(k)*Z(1,l)+S2h000082221111(k)*Z(2,l)+(Inv602
     &  21111+S2000000622111(1)-S2h000081221111(1)-S2h000082222111(1))*Z
     &  (k,l)-4*(2*(S2h000000612211(1)*ZZ(k,1,l,1)+S2h000000622211(1)*ZZ
     &  (k,1,l,2))+S2h000000612111(1)*ZZ(k,2,l,1)+S2h000000622111(1)*ZZ(
     &  k,2,l,2))
       aux0000006222(1,1,1)=-(F(1)*S2h00007122211(1))-F(2)*S2h0000722221
     &  1(1)+S2h000081222111(k)*Z(1,l)+S2h000082222111(k)*Z(2,l)+(Inv602
     &  22111+S2000000622211(1)-S2h000081222111(1)-S2h000082222211(1))*Z
     &  (k,l)-6*(S2h000000612221(1)*ZZ(k,1,l,1)+S2h000000622221(1)*ZZ(k,
     &  1,l,2)+S2h000000612211(1)*ZZ(k,2,l,1)+S2h000000622211(1)*ZZ(k,2,
     &  l,2))
       aux0000006222(2,1,1)=-(F(1)*S2h00007122221(1))-F(2)*S2h0000722222
     &  1(1)+S2h000081222211(k)*Z(1,l)+S2h000082222211(k)*Z(2,l)+(Inv602
     &  22211+S2000000622221(1)-S2h000081222211(1)-S2h000082222221(1))*Z
     &  (k,l)-4*(S2h000000612222(1)*ZZ(k,1,l,1)+S2h000000622222(1)*ZZ(k,
     &  1,l,2)+2*(S2h000000612221(1)*ZZ(k,2,l,1)+S2h000000622221(1)*ZZ(k
     &  ,2,l,2)))
       aux0000006222(2,2,1)=-(F(1)*S2h00007122222(1))-F(2)*S2h0000722222
     &  2(1)+S2h000081222221(k)*Z(1,l)+S2h000082222221(k)*Z(2,l)+(Inv602
     &  22221+S2000000622222(1)-S2h000081222221(1)-S2h000082222222(1))*Z
     &  (k,l)-2*(S2h000000612222(2)*ZZ(k,1,l,1)+S2h000000622222(2)*ZZ(k,
     &  1,l,2)+5*(S2h000000612222(1)*ZZ(k,2,l,1)+S2h000000622222(1)*ZZ(k
     &  ,2,l,2)))
       aux0000006222(2,2,2)=-(F(1)*S2h00007122222(2))-F(2)*S2h0000722222
     &  2(2)+S2h000081222222(k)*Z(1,l)+S2h000082222222(k)*Z(2,l)+(Inv602
     &  22222+S2000000622222(2)-S2h000081222222(1)-S2h000082222222(2))*Z
     &  (k,l)-12*(S2h000000612222(2)*ZZ(k,2,l,1)+S2h000000622222(2)*ZZ(k
     &  ,2,l,2))
       temp0000006111(1,1,1)=I38Z*(aux0000006111(1,1,1)+12*F(4)*temp0000
     &  00511(1,1,1)+F(3)*temp00006111(1,1,1)+120*temp0000000041(1,1,1)*
     &  ZZ(k,1,l,1))
       temp0000006211(1,1,1)=I38Z*(aux0000006211(1,1,1)+F(5)*temp0000005
     &  11(1,1,1)+10*F(4)*temp000000521(1,1,1)+F(3)*temp00006211(1,1,1)+
     &  80*temp0000000042(1,1,1)*ZZ(k,1,l,1)+20*temp0000000041(1,1,1)*(Z
     &  Z(k,1,l,2)+ZZ(k,2,l,1)))
       temp0000006221(1,1,1)=I38Z*(aux0000006221(1,1,1)+2*F(5)*temp00000
     &  0521(1,1,1)+F(3)*temp00006221(1,1,1)+48*temp0000000042(2,1,1)*ZZ
     &  (k,1,l,1)+32*temp0000000042(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+8*(
     &  F(4)*temp000000522(1,1,1)+temp0000000041(1,1,1)*ZZ(k,2,l,2)))
       temp0000006222(1,1,1)=I38Z*(aux0000006222(1,1,1)+6*F(4)*temp00000
     &  0522(2,1,1)+F(3)*temp00006222(1,1,1)+36*temp0000000042(2,1,1)*(Z
     &  Z(k,1,l,2)+ZZ(k,2,l,1))+temp000000522(1,1,1)*(6*r10*(ZZ(k,1,l,2)
     &  +ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2))+24*(temp0000000042(2,2,1)*ZZ(k
     &  ,1,l,1)+temp0000000042(1,1,1)*ZZ(k,2,l,2)))
       temp0000006222(2,1,1)=I38Z*(aux0000006222(2,1,1)+4*F(5)*temp00000
     &  0522(2,1,1)+4*F(4)*temp000000522(2,2,1)+F(3)*temp00006222(2,1,1)
     &  +8*temp0000000042(2,2,2)*ZZ(k,1,l,1)+32*temp0000000042(2,2,1)*(Z
     &  Z(k,1,l,2)+ZZ(k,2,l,1))+48*temp0000000042(2,1,1)*ZZ(k,2,l,2))
       temp0000006222(2,2,1)=I38Z*(aux0000006222(2,2,1)+2*F(4)*temp00000
     &  0522(2,2,2)+F(3)*temp00006222(2,2,1)+20*temp0000000042(2,2,2)*(Z
     &  Z(k,1,l,2)+ZZ(k,2,l,1))+80*temp0000000042(2,2,1)*ZZ(k,2,l,2)+tem
     &  p000000522(2,2,1)*(10*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+20*r21*ZZ(k,
     &  2,l,2)))
       temp0000006222(2,2,2)=I38Z*(aux0000006222(2,2,2)+F(3)*temp0000622
     &  2(2,2,2)+120*temp0000000042(2,2,2)*ZZ(k,2,l,2)+temp000000522(2,2
     &  ,2)*(12*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+24*r21*ZZ(k,2,l,2)))
       temp0000006111(1,1,2)=temp0000006211(1,1,1)
       temp0000006111(1,2,1)=temp0000006211(1,1,1)
       temp0000006111(1,2,2)=temp0000006221(1,1,1)
       temp0000006211(1,1,2)=temp0000006221(1,1,1)
       temp0000006211(1,2,1)=temp0000006221(1,1,1)
       temp0000006211(1,2,2)=temp0000006222(1,1,1)
       temp0000006221(1,1,2)=temp0000006222(1,1,1)
       temp0000006221(1,2,1)=temp0000006222(1,1,1)
       temp0000006221(1,2,2)=temp0000006222(2,1,1)
       temp0000006222(1,1,2)=temp0000006222(2,1,1)
       temp0000006222(1,2,1)=temp0000006222(2,1,1)
       temp0000006222(1,2,2)=temp0000006222(2,2,1)
       temp0000006222(2,1,2)=temp0000006222(2,2,1)
       S2000081111111(1)=2*B23(2,2)
       S2000081111111(2)=-2*B23(2,3)
       S2000081211111(1)=-2*B23(2,3)
       S2000081211111(2)=2*B23(2,4)
       S2000081221111(1)=2*B23(2,4)
       S2000081221111(2)=-2*B23(2,5)
       S2000081222111(1)=-2*B23(2,5)
       S2000081222111(2)=2*B23(2,6)
       S2000081222211(1)=2*B23(2,6)
       S2000081222211(2)=-2*B23(2,7)
       S2000081222221(1)=-2*B23(2,7)
       S2000081222221(2)=2*B23(2,8)
       S2000081222222(1)=2*B23(2,8)
       S2000081222222(2)=-2*B23(2,9)
       S2000082111111(1)=-2*B23(2,3)
       S2000082111111(2)=2*B23(2,4)
       S2000082211111(1)=2*B23(2,4)
       S2000082211111(2)=-2*B23(2,5)
       S2000082221111(1)=-2*B23(2,5)
       S2000082221111(2)=2*B23(2,6)
       S2000082222111(1)=2*B23(2,6)
       S2000082222111(2)=-2*B23(2,7)
       S2000082222211(1)=-2*B23(2,7)
       S2000082222211(2)=2*B23(2,8)
       S2000082222221(1)=2*B23(2,8)
       S2000082222221(2)=-2*B23(2,9)
       S2000082222222(1)=-2*B23(2,9)
       S2000082222222(2)=2*B23(2,10)
       S2h0010111111111(1)=B13(2,11)+B23(2,2)
       S2h0010111111111(2)=B13(2,11)-B23(2,3)
       S2h0010121111111(1)=B13(2,11)-B23(2,3)
       S2h0010121111111(2)=B13(2,11)+B23(2,4)
       S2h0010122111111(1)=B13(2,11)+B23(2,4)
       S2h0010122111111(2)=B13(2,11)-B23(2,5)
       S2h0010122211111(1)=B13(2,11)-B23(2,5)
       S2h0010122211111(2)=B13(2,11)+B23(2,6)
       S2h0010122221111(1)=B13(2,11)+B23(2,6)
       S2h0010122221111(2)=B13(2,11)-B23(2,7)
       S2h0010122222111(1)=B13(2,11)-B23(2,7)
       S2h0010122222111(2)=B13(2,11)+B23(2,8)
       S2h0010122222211(1)=B13(2,11)+B23(2,8)
       S2h0010122222211(2)=B13(2,11)-B23(2,9)
       S2h0010122222221(1)=B13(2,11)-B23(2,9)
       S2h0010122222221(2)=B13(2,11)+B23(2,10)
       S2h0010122222222(1)=B13(2,11)+B23(2,10)
       S2h0010122222222(2)=B13(2,11)-B23(2,11)
       S2h0010211111111(1)=B12(2,11)-B13(2,11)
       S2h0010211111111(2)=-B13(2,11)
       S2h0010221111111(1)=-B13(2,11)
       S2h0010221111111(2)=-B13(2,11)
       S2h0010222111111(1)=-B13(2,11)
       S2h0010222111111(2)=-B13(2,11)
       S2h0010222211111(1)=-B13(2,11)
       S2h0010222211111(2)=-B13(2,11)
       S2h0010222221111(1)=-B13(2,11)
       S2h0010222221111(2)=-B13(2,11)
       S2h0010222222111(1)=-B13(2,11)
       S2h0010222222111(2)=-B13(2,11)
       S2h0010222222211(1)=-B13(2,11)
       S2h0010222222211(2)=-B13(2,11)
       S2h0010222222221(1)=-B13(2,11)
       S2h0010222222221(2)=-B13(2,11)
       S2h0010222222222(1)=-B13(2,11)
       S2h0010222222222(2)=-B13(2,11)
       aux0000811111(1,1,1)=-(F(1)*S2h00911111111(1))-F(2)*S2h0092111111
     &  1(1)+S2h0010111111111(k)*Z(1,l)+S2h0010211111111(k)*Z(2,l)+(Inv40
     &  11111111+S2000081111111(1)-S2h0010111111111(1)-S2h00102211111
     &  11(1))*Z(k,l)-16*(S2h000081111111(1)*ZZ(k,1,l,1)+S2h000082111111
     &  (1)*ZZ(k,1,l,2))
       aux0000821111(1,1,1)=-(F(1)*S2h00912111111(1))-F(2)*S2h0092211111
     &  1(1)+S2h0010121111111(k)*Z(1,l)+S2h0010221111111(k)*Z(2,l)+(Inv40
     &  21111111+S2000082111111(1)-S2h0010121111111(1)-S2h00102221111
     &  11(1))*Z(k,l)-14*(S2h000081211111(1)*ZZ(k,1,l,1)+S2h000082211111
     &  (1)*ZZ(k,1,l,2))-2*(S2h000081111111(1)*ZZ(k,2,l,1)+S2h0000821111
     &  11(1)*ZZ(k,2,l,2))
       aux0000822111(1,1,1)=-(F(1)*S2h00912211111(1))-F(2)*S2h0092221111
     &  1(1)+S2h0010122111111(k)*Z(1,l)+S2h0010222111111(k)*Z(2,l)+(Inv40
     &  22111111+S2000082211111(1)-S2h0010122111111(1)-S2h00102222111
     &  11(1))*Z(k,l)-12*(S2h000081221111(1)*ZZ(k,1,l,1)+S2h000082221111
     &  (1)*ZZ(k,1,l,2))-4*(S2h000081211111(1)*ZZ(k,2,l,1)+S2h0000822111
     &  11(1)*ZZ(k,2,l,2))
       aux0000822211(1,1,1)=-(F(1)*S2h00912221111(1))-F(2)*S2h0092222111
     &  1(1)+S2h0010122211111(k)*Z(1,l)+S2h0010222211111(k)*Z(2,l)+(Inv40
     &  22211111+S2000082221111(1)-S2h0010122211111(1)-S2h00102222211
     &  11(1))*Z(k,l)-10*(S2h000081222111(1)*ZZ(k,1,l,1)+S2h000082222111
     &  (1)*ZZ(k,1,l,2))-6*(S2h000081221111(1)*ZZ(k,2,l,1)+S2h0000822211
     &  11(1)*ZZ(k,2,l,2))
       aux0000822221(1,1,1)=-(F(1)*S2h00912222111(1))-F(2)*S2h0092222211
     &  1(1)+S2h0010122221111(k)*Z(1,l)+S2h0010222221111(k)*Z(2,l)+(Inv40
     &  22221111+S2000082222111(1)-S2h0010122221111(1)-S2h00102222221
     &  11(1))*Z(k,l)-8*(S2h000081222211(1)*ZZ(k,1,l,1)+S2h000082222211(
     &  1)*ZZ(k,1,l,2)+S2h000081222111(1)*ZZ(k,2,l,1))-8*S2h000082222111
     &  (1)*ZZ(k,2,l,2)
       aux0000822222(1,1,1)=-(F(1)*S2h00912222211(1))-F(2)*S2h0092222221
     &  1(1)+S2h0010122222111(k)*Z(1,l)+S2h0010222222111(k)*Z(2,l)+(Inv40
     &  22222111+S2000082222211(1)-S2h0010122222111(1)-S2h00102222222
     &  11(1))*Z(k,l)-6*(S2h000081222221(1)*ZZ(k,1,l,1)+S2h000082222221(
     &  1)*ZZ(k,1,l,2))-10*(S2h000081222211(1)*ZZ(k,2,l,1)+S2h0000822222
     &  11(1)*ZZ(k,2,l,2))
       aux0000822222(2,1,1)=-(F(1)*S2h00912222221(1))-F(2)*S2h0092222222
     &  1(1)+S2h0010122222211(k)*Z(1,l)+S2h0010222222211(k)*Z(2,l)+(Inv40
     &  22222211+S2000082222221(1)-S2h0010122222211(1)-S2h00102222222
     &  21(1))*Z(k,l)-4*(S2h000081222222(1)*ZZ(k,1,l,1)+S2h000082222222(
     &  1)*ZZ(k,1,l,2))-12*(S2h000081222221(1)*ZZ(k,2,l,1)+S2h0000822222
     &  21(1)*ZZ(k,2,l,2))
       aux0000822222(2,2,1)=-(F(1)*S2h00912222222(1))-F(2)*S2h0092222222
     &  2(1)+S2h0010122222221(k)*Z(1,l)+S2h0010222222221(k)*Z(2,l)+(Inv40
     &  22222221+S2000082222222(1)-S2h0010122222221(1)-S2h00102222222
     &  22(1))*Z(k,l)-2*(S2h000081222222(2)*ZZ(k,1,l,1)+S2h000082222222(
     &  2)*ZZ(k,1,l,2))-14*(S2h000081222222(1)*ZZ(k,2,l,1)+S2h0000822222
     &  22(1)*ZZ(k,2,l,2))
       aux0000822222(2,2,2)=-(F(1)*S2h00912222222(2))-F(2)*S2h0092222222
     &  2(2)+S2h0010122222222(k)*Z(1,l)+S2h0010222222222(k)*Z(2,l)+(Inv40
     &  22222222+S2000082222222(2)-S2h0010122222222(1)-S2h00102222222
     &  22(2))*Z(k,l)-16*(S2h000081222222(2)*ZZ(k,2,l,1)+S2h000082222222
     &  (2)*ZZ(k,2,l,2))
       temp0000811111(1,1,1)=I42Z*(aux0000811111(1,1,1)+16*F(4)*temp0000
     &  71111(1,1,1)+F(3)*temp00811111(1,1,1)+224*temp0000006111(1,1,1)*
     &  ZZ(k,1,l,1))
       temp0000821111(1,1,1)=I42Z*(aux0000821111(1,1,1)+F(5)*temp0000711
     &  11(1,1,1)+14*F(4)*temp000072111(1,1,1)+F(3)*temp00821111(1,1,1)+
     &  168*temp0000006211(1,1,1)*ZZ(k,1,l,1)+28*temp0000006111(1,1,1)*(
     &  ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp0000822111(1,1,1)=I42Z*(aux0000822111(1,1,1)+2*F(5)*temp00007
     &  2111(1,1,1)+12*F(4)*temp000072211(1,1,1)+F(3)*temp00822111(1,1,1
     &  )+120*temp0000006221(1,1,1)*ZZ(k,1,l,1)+48*temp0000006211(1,1,1)
     &  *(ZZ(k,1,l,2)+ZZ(k,2,l,1))+8*temp0000006111(1,1,1)*ZZ(k,2,l,2))
       temp0000822211(1,1,1)=I42Z*(aux0000822211(1,1,1)+10*F(4)*temp0000
     &  72221(1,1,1)+F(3)*temp00822211(1,1,1)+80*temp0000006222(1,1,1)*Z
     &  Z(k,1,l,1)+60*temp0000006221(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+24
     &  *temp0000006211(1,1,1)*ZZ(k,2,l,2)+temp000072211(1,1,1)*(6*r10*(
     &  ZZ(k,1,l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2)))
       temp0000822221(1,1,1)=I42Z*(aux0000822221(1,1,1)+4*F(5)*temp00007
     &  2221(1,1,1)+8*F(4)*temp000072222(1,1,1)+F(3)*temp00822221(1,1,1)
     &  +64*temp0000006222(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+48*(temp0000
     &  006222(2,1,1)*ZZ(k,1,l,1)+temp0000006221(1,1,1)*ZZ(k,2,l,2)))
       temp0000822222(1,1,1)=I42Z*(aux0000822222(1,1,1)+6*F(4)*temp00007
     &  2222(2,1,1)+F(3)*temp00822222(1,1,1)+24*temp0000006222(2,2,1)*ZZ
     &  (k,1,l,1)+60*temp0000006222(2,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+80*
     &  temp0000006222(1,1,1)*ZZ(k,2,l,2)+temp000072222(1,1,1)*(10*r10*(
     &  ZZ(k,1,l,2)+ZZ(k,2,l,1))+20*r21*ZZ(k,2,l,2)))
       temp0000822222(2,1,1)=I42Z*(aux0000822222(2,1,1)+4*F(4)*temp00007
     &  2222(2,2,1)+F(3)*temp00822222(2,1,1)+8*temp0000006222(2,2,2)*ZZ(
     &  k,1,l,1)+48*temp0000006222(2,2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+120*
     &  temp0000006222(2,1,1)*ZZ(k,2,l,2)+temp000072222(2,1,1)*(12*r10*(
     &  ZZ(k,1,l,2)+ZZ(k,2,l,1))+24*r21*ZZ(k,2,l,2)))
       temp0000822222(2,2,1)=I42Z*(aux0000822222(2,2,1)+2*F(4)*temp00007
     &  2222(2,2,2)+F(3)*temp00822222(2,2,1)+28*temp0000006222(2,2,2)*(Z
     &  Z(k,1,l,2)+ZZ(k,2,l,1))+168*temp0000006222(2,2,1)*ZZ(k,2,l,2)+te
     &  mp000072222(2,2,1)*(14*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+28*r21*ZZ(k
     &  ,2,l,2)))
       temp0000822222(2,2,2)=I42Z*(aux0000822222(2,2,2)+8*F(5)*temp00007
     &  2222(2,2,2)+F(3)*temp00822222(2,2,2)+224*temp0000006222(2,2,2)*Z
     &  Z(k,2,l,2))
       temp0000811111(1,1,2)=temp0000821111(1,1,1)
       temp0000811111(1,2,1)=temp0000821111(1,1,1)
       temp0000811111(1,2,2)=temp0000822111(1,1,1)
       temp0000821111(1,1,2)=temp0000822111(1,1,1)
       temp0000821111(1,2,1)=temp0000822111(1,1,1)
       temp0000821111(1,2,2)=temp0000822211(1,1,1)
       temp0000822111(1,1,2)=temp0000822211(1,1,1)
       temp0000822111(1,2,1)=temp0000822211(1,1,1)
       temp0000822111(1,2,2)=temp0000822221(1,1,1)
       temp0000822211(1,1,2)=temp0000822221(1,1,1)
       temp0000822211(1,2,1)=temp0000822221(1,1,1)
       temp0000822211(1,2,2)=temp0000822222(1,1,1)
       temp0000822221(1,1,2)=temp0000822222(1,1,1)
       temp0000822221(1,2,1)=temp0000822222(1,1,1)
       temp0000822221(1,2,2)=temp0000822222(2,1,1)
       temp0000822222(1,1,2)=temp0000822222(2,1,1)
       temp0000822222(1,2,1)=temp0000822222(2,1,1)
       temp0000822222(1,2,2)=temp0000822222(2,2,1)
       temp0000822222(2,1,2)=temp0000822222(2,2,1)
       S20010111111111(1)=2*B023
       S20010111111111(2)=-2*B23(1,1)
       S20010121111111(1)=-2*B23(1,1)
       S20010121111111(2)=2*B23(1,2)
       S20010122111111(1)=2*B23(1,2)
       S20010122111111(2)=-2*B23(1,3)
       S20010122211111(1)=-2*B23(1,3)
       S20010122211111(2)=2*B23(1,4)
       S20010122221111(1)=2*B23(1,4)
       S20010122221111(2)=-2*B23(1,5)
       S20010122222111(1)=-2*B23(1,5)
       S20010122222111(2)=2*B23(1,6)
       S20010122222211(1)=2*B23(1,6)
       S20010122222211(2)=-2*B23(1,7)
       S20010122222221(1)=-2*B23(1,7)
       S20010122222221(2)=2*B23(1,8)
       S20010122222222(1)=2*B23(1,8)
       S20010122222222(2)=-2*B23(1,9)
       S20010211111111(1)=-2*B23(1,1)
       S20010211111111(2)=2*B23(1,2)
       S20010221111111(1)=2*B23(1,2)
       S20010221111111(2)=-2*B23(1,3)
       S20010222111111(1)=-2*B23(1,3)
       S20010222111111(2)=2*B23(1,4)
       S20010222211111(1)=2*B23(1,4)
       S20010222211111(2)=-2*B23(1,5)
       S20010222221111(1)=-2*B23(1,5)
       S20010222221111(2)=2*B23(1,6)
       S20010222222111(1)=2*B23(1,6)
       S20010222222111(2)=-2*B23(1,7)
       S20010222222211(1)=-2*B23(1,7)
       S20010222222211(2)=2*B23(1,8)
       S20010222222221(1)=2*B23(1,8)
       S20010222222221(2)=-2*B23(1,9)
       S20010222222222(1)=-2*B23(1,9)
       S20010222222222(2)=2*B23(1,10)
       S21211111111111(1)=B023+B13(1,11)
       S21211111111111(2)=B13(1,11)-B23(1,1)
       S21212111111111(1)=B13(1,11)-B23(1,1)
       S21212111111111(2)=B13(1,11)+B23(1,2)
       S21212211111111(1)=B13(1,11)+B23(1,2)
       S21212211111111(2)=B13(1,11)-B23(1,3)
       S21212221111111(1)=B13(1,11)-B23(1,3)
       S21212221111111(2)=B13(1,11)+B23(1,4)
       S21212222111111(1)=B13(1,11)+B23(1,4)
       S21212222111111(2)=B13(1,11)-B23(1,5)
       S21212222211111(1)=B13(1,11)-B23(1,5)
       S21212222211111(2)=B13(1,11)+B23(1,6)
       S21212222221111(1)=B13(1,11)+B23(1,6)
       S21212222221111(2)=B13(1,11)-B23(1,7)
       S21212222222111(1)=B13(1,11)-B23(1,7)
       S21212222222111(2)=B13(1,11)+B23(1,8)
       S21212222222211(1)=B13(1,11)+B23(1,8)
       S21212222222211(2)=B13(1,11)-B23(1,9)
       S21212222222221(1)=B13(1,11)-B23(1,9)
       S21212222222221(2)=B13(1,11)+B23(1,10)
       S21212222222222(1)=B13(1,11)+B23(1,10)
       S21212222222222(2)=B13(1,11)-B23(1,11)
       S21221111111111(1)=B12(1,11)-B13(1,11)
       S21221111111111(2)=-B13(1,11)
       S21222111111111(1)=-B13(1,11)
       S21222111111111(2)=-B13(1,11)
       S21222211111111(1)=-B13(1,11)
       S21222211111111(2)=-B13(1,11)
       S21222221111111(1)=-B13(1,11)
       S21222221111111(2)=-B13(1,11)
       S21222222111111(1)=-B13(1,11)
       S21222222111111(2)=-B13(1,11)
       S21222222211111(1)=-B13(1,11)
       S21222222211111(2)=-B13(1,11)
       S21222222221111(1)=-B13(1,11)
       S21222222221111(2)=-B13(1,11)
       S21222222222111(1)=-B13(1,11)
       S21222222222111(2)=-B13(1,11)
       S21222222222211(1)=-B13(1,11)
       S21222222222211(2)=-B13(1,11)
       S21222222222221(1)=-B13(1,11)
       S21222222222221(2)=-B13(1,11)
       S21222222222222(1)=-B13(1,11)
       S21222222222222(2)=-B13(1,11)
       aux00101111111(1,1,1)=-(F(1)*S2111111111111(1))-F(2)*S21121111111
     &  11(1)+S21211111111111(k)*Z(1,l)+S21221111111111(k)*Z(2,l)+(Inv20
     &  1111111111+S20010111111111(1)-S21211111111111(1)-S21222111111111
     &  (1))*Z(k,l)-20*(S2h0010111111111(1)*ZZ(k,1,l,1)+S2h0010211111111
     &  (1)*ZZ(k,1,l,2))
       aux00102111111(1,1,1)=-(F(1)*S2111211111111(1))-F(2)*S21122111111
     &  11(1)+S21212111111111(k)*Z(1,l)+S21222111111111(k)*Z(2,l)+(Inv20
     &  2111111111+S20010211111111(1)-S21212111111111(1)-S21222211111111
     &  (1))*Z(k,l)-18*(S2h0010121111111(1)*ZZ(k,1,l,1)+S2h0010221111111
     &  (1)*ZZ(k,1,l,2))-2*(S2h0010111111111(1)*ZZ(k,2,l,1)+S2h001021111
     &  1111(1)*ZZ(k,2,l,2))
       aux00102211111(1,1,1)=-(F(1)*S2111221111111(1))-F(2)*S21122211111
     &  11(1)+S21212211111111(k)*Z(1,l)+S21222211111111(k)*Z(2,l)+(Inv20
     &  2211111111+S20010221111111(1)-S21212211111111(1)-S21222221111111
     &  (1))*Z(k,l)-16*(S2h0010122111111(1)*ZZ(k,1,l,1)+S2h0010222111111
     &  (1)*ZZ(k,1,l,2))-4*(S2h0010121111111(1)*ZZ(k,2,l,1)+S2h001022111
     &  1111(1)*ZZ(k,2,l,2))
       aux00102221111(1,1,1)=-(F(1)*S2111222111111(1))-F(2)*S21122221111
     &  11(1)+S21212221111111(k)*Z(1,l)+S21222221111111(k)*Z(2,l)+(Inv20
     &  2221111111+S20010222111111(1)-S21212221111111(1)-S21222222111111
     &  (1))*Z(k,l)-14*(S2h0010122211111(1)*ZZ(k,1,l,1)+S2h0010222211111
     &  (1)*ZZ(k,1,l,2))-6*(S2h0010122111111(1)*ZZ(k,2,l,1)+S2h001022211
     &  1111(1)*ZZ(k,2,l,2))
       aux00102222111(1,1,1)=-(F(1)*S2111222211111(1))-F(2)*S21122222111
     &  11(1)+S21212222111111(k)*Z(1,l)+S21222222111111(k)*Z(2,l)+(Inv20
     &  2222111111+S20010222211111(1)-S21212222111111(1)-S21222222211111
     &  (1))*Z(k,l)-12*(S2h0010122221111(1)*ZZ(k,1,l,1)+S2h0010222221111
     &  (1)*ZZ(k,1,l,2))-8*(S2h0010122211111(1)*ZZ(k,2,l,1)+S2h001022221
     &  1111(1)*ZZ(k,2,l,2))
       aux00102222211(1,1,1)=-(F(1)*S2111222221111(1))-F(2)*S21122222211
     &  11(1)+S21212222211111(k)*Z(1,l)+S21222222211111(k)*Z(2,l)+(Inv20
     &  2222211111+S20010222221111(1)-S21212222211111(1)-S21222222221111
     &  (1))*Z(k,l)-10*(S2h0010122222111(1)*ZZ(k,1,l,1)+S2h0010222222111
     &  (1)*ZZ(k,1,l,2)+S2h0010122221111(1)*ZZ(k,2,l,1))-10*S2h001022222
     &  1111(1)*ZZ(k,2,l,2)
       aux00102222221(1,1,1)=-(F(1)*S2111222222111(1))-F(2)*S21122222221
     &  11(1)+S21212222221111(k)*Z(1,l)+S21222222221111(k)*Z(2,l)+(Inv20
     &  2222221111+S20010222222111(1)-S21212222221111(1)-S21222222222111
     &  (1))*Z(k,l)-8*(S2h0010122222211(1)*ZZ(k,1,l,1)+S2h0010222222211(
     &  1)*ZZ(k,1,l,2))-12*(S2h0010122222111(1)*ZZ(k,2,l,1)+S2h001022222
     &  2111(1)*ZZ(k,2,l,2))
       aux00102222222(1,1,1)=-(F(1)*S2111222222211(1))-F(2)*S21122222222
     &  11(1)+S21212222222111(k)*Z(1,l)+S21222222222111(k)*Z(2,l)+(Inv20
     &  2222222111+S20010222222211(1)-S21212222222111(1)-S21222222222211
     &  (1))*Z(k,l)-6*(S2h0010122222221(1)*ZZ(k,1,l,1)+S2h0010222222221(
     &  1)*ZZ(k,1,l,2))-14*(S2h0010122222211(1)*ZZ(k,2,l,1)+S2h001022222
     &  2211(1)*ZZ(k,2,l,2))
       aux00102222222(2,1,1)=-(F(1)*S2111222222221(1))-F(2)*S21122222222
     &  21(1)+S21212222222211(k)*Z(1,l)+S21222222222211(k)*Z(2,l)+(Inv20
     &  2222222211+S20010222222221(1)-S21212222222211(1)-S21222222222221
     &  (1))*Z(k,l)-4*(S2h0010122222222(1)*ZZ(k,1,l,1)+S2h0010222222222(
     &  1)*ZZ(k,1,l,2))-16*(S2h0010122222221(1)*ZZ(k,2,l,1)+S2h001022222
     &  2221(1)*ZZ(k,2,l,2))
       aux00102222222(2,2,1)=-(F(1)*S2111222222222(1))-F(2)*S21122222222
     &  22(1)+S21212222222221(k)*Z(1,l)+S21222222222221(k)*Z(2,l)+(Inv20
     &  2222222221+S20010222222222(1)-S21212222222221(1)-S21222222222222
     &  (1))*Z(k,l)-2*(S2h0010122222222(2)*ZZ(k,1,l,1)+S2h0010222222222(
     &  2)*ZZ(k,1,l,2))-18*(S2h0010122222222(1)*ZZ(k,2,l,1)+S2h001022222
     &  2222(1)*ZZ(k,2,l,2))
       aux00102222222(2,2,2)=-(F(1)*S2111222222222(2))-F(2)*S21122222222
     &  22(2)+S21212222222222(k)*Z(1,l)+S21222222222222(k)*Z(2,l)+(Inv20
     &  2222222222+S20010222222222(2)-S21212222222222(1)-S21222222222222
     &  (2))*Z(k,l)-20*(S2h0010122222222(2)*ZZ(k,2,l,1)+S2h0010222222222
     &  (2)*ZZ(k,2,l,2))
       temp00101111111(1,1,1)=I46Z*(aux00101111111(1,1,1)+20*F(4)*temp00
     &  9111111(1,1,1)+F(3)*temp101111111(1,1,1)+360*temp0000811111(1,1,
     &  1)*ZZ(k,1,l,1))
       temp00102111111(1,1,1)=I46Z*(aux00102111111(1,1,1)+F(5)*temp00911
     &  1111(1,1,1)+18*F(4)*temp009211111(1,1,1)+F(3)*temp102111111(1,1,
     &  1)+288*temp0000821111(1,1,1)*ZZ(k,1,l,1)+36*temp0000811111(1,1,1
     &  )*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp00102211111(1,1,1)=I46Z*(aux00102211111(1,1,1)+2*F(5)*temp009
     &  211111(1,1,1)+16*F(4)*temp009221111(1,1,1)+F(3)*temp102211111(1,
     &  1,1)+224*temp0000822111(1,1,1)*ZZ(k,1,l,1)+64*temp0000821111(1,1
     &  ,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+8*temp0000811111(1,1,1)*ZZ(k,2,l,2
     &  ))
       temp00102221111(1,1,1)=I46Z*(aux00102221111(1,1,1)+14*F(4)*temp00
     &  9222111(1,1,1)+F(3)*temp102221111(1,1,1)+168*temp0000822211(1,1,
     &  1)*ZZ(k,1,l,1)+84*temp0000822111(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)
     &  )+24*temp0000821111(1,1,1)*ZZ(k,2,l,2)+temp009221111(1,1,1)*(6*r
     &  10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2)))
       temp00102222111(1,1,1)=I46Z*(aux00102222111(1,1,1)+4*F(5)*temp009
     &  222111(1,1,1)+12*F(4)*temp009222211(1,1,1)+F(3)*temp102222111(1,
     &  1,1)+120*temp0000822221(1,1,1)*ZZ(k,1,l,1)+96*temp0000822211(1,1
     &  ,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+48*temp0000822111(1,1,1)*ZZ(k,2,l,
     &  2))
       temp00102222211(1,1,1)=I46Z*(aux00102222211(1,1,1)+10*F(4)*temp00
     &  9222221(1,1,1)+F(3)*temp102222211(1,1,1)+100*temp0000822221(1,1,
     &  1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp009222211(1,1,1)*(10*r10*(ZZ(k,
     &  1,l,2)+ZZ(k,2,l,1))+20*r21*ZZ(k,2,l,2))+80*(temp0000822222(1,1,1
     &  )*ZZ(k,1,l,1)+temp0000822211(1,1,1)*ZZ(k,2,l,2)))
       temp00102222221(1,1,1)=I46Z*(aux00102222221(1,1,1)+8*F(4)*temp009
     &  222222(1,1,1)+F(3)*temp102222221(1,1,1)+48*temp0000822222(2,1,1)
     &  *ZZ(k,1,l,1)+96*temp0000822222(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+
     &  120*temp0000822221(1,1,1)*ZZ(k,2,l,2)+temp009222221(1,1,1)*(12*r
     &  10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+24*r21*ZZ(k,2,l,2)))
       temp00102222222(1,1,1)=I46Z*(aux00102222222(1,1,1)+6*F(4)*temp009
     &  222222(2,1,1)+F(3)*temp102222222(1,1,1)+24*temp0000822222(2,2,1)
     &  *ZZ(k,1,l,1)+84*temp0000822222(2,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+
     &  168*temp0000822222(1,1,1)*ZZ(k,2,l,2)+temp009222222(1,1,1)*(14*r
     &  10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+28*r21*ZZ(k,2,l,2)))
       temp00102222222(2,1,1)=I46Z*(aux00102222222(2,1,1)+4*F(4)*temp009
     &  222222(2,2,1)+F(3)*temp102222222(2,1,1)+8*(F(5)*temp009222222(2,
     &  1,1)+temp0000822222(2,2,2)*ZZ(k,1,l,1))+64*temp0000822222(2,2,1)
     &  *(ZZ(k,1,l,2)+ZZ(k,2,l,1))+224*temp0000822222(2,1,1)*ZZ(k,2,l,2)
     &  )
       temp00102222222(2,2,1)=I46Z*(aux00102222222(2,2,1)+2*F(4)*temp009
     &  222222(2,2,2)+F(3)*temp102222222(2,2,1)+36*temp0000822222(2,2,2)
     &  *(ZZ(k,1,l,2)+ZZ(k,2,l,1))+288*temp0000822222(2,2,1)*ZZ(k,2,l,2)
     &  +temp009222222(2,2,1)*(18*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+36*r21*Z
     &  Z(k,2,l,2)))
       temp00102222222(2,2,2)=I46Z*(aux00102222222(2,2,2)+F(3)*temp10222
     &  2222(2,2,2)+360*temp0000822222(2,2,2)*ZZ(k,2,l,2)+temp009222222(
     &  2,2,2)*(20*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+40*r21*ZZ(k,2,l,2)))
       temp00101111111(1,1,2)=temp00102111111(1,1,1)
       temp00101111111(1,2,1)=temp00102111111(1,1,1)
       temp00101111111(1,2,2)=temp00102211111(1,1,1)
       temp00102111111(1,1,2)=temp00102211111(1,1,1)
       temp00102111111(1,2,1)=temp00102211111(1,1,1)
       temp00102111111(1,2,2)=temp00102221111(1,1,1)
       temp00102211111(1,1,2)=temp00102221111(1,1,1)
       temp00102211111(1,2,1)=temp00102221111(1,1,1)
       temp00102211111(1,2,2)=temp00102222111(1,1,1)
       temp00102221111(1,1,2)=temp00102222111(1,1,1)
       temp00102221111(1,2,1)=temp00102222111(1,1,1)
       temp00102221111(1,2,2)=temp00102222211(1,1,1)
       temp00102222111(1,1,2)=temp00102222211(1,1,1)
       temp00102222111(1,2,1)=temp00102222211(1,1,1)
       temp00102222111(1,2,2)=temp00102222221(1,1,1)
       temp00102222211(1,1,2)=temp00102222221(1,1,1)
       temp00102222211(1,2,1)=temp00102222221(1,1,1)
       temp00102222211(1,2,2)=temp00102222222(1,1,1)
       temp00102222221(1,1,2)=temp00102222222(1,1,1)
       temp00102222221(1,2,1)=temp00102222222(1,1,1)
       temp00102222221(1,2,2)=temp00102222222(2,1,1)
       temp00102222222(1,1,2)=temp00102222222(2,1,1)
       temp00102222222(1,2,1)=temp00102222222(2,1,1)
       temp00102222222(1,2,2)=temp00102222222(2,2,1)
       temp00102222222(2,1,2)=temp00102222222(2,2,1)
       aux1111111111(1,1,1)=-(S21211111111111(1)*Z(jj,1))-S2122111111111
     &  1(1)*Z(jj,2)
       aux1121111111(1,1,1)=-(S21212111111111(1)*Z(jj,1))-S2122211111111
     &  1(1)*Z(jj,2)
       aux1122111111(1,1,1)=-(S21212211111111(1)*Z(jj,1))-S2122221111111
     &  1(1)*Z(jj,2)
       aux1122211111(1,1,1)=-(S21212221111111(1)*Z(jj,1))-S2122222111111
     &  1(1)*Z(jj,2)
       aux1122221111(1,1,1)=-(S21212222111111(1)*Z(jj,1))-S2122222211111
     &  1(1)*Z(jj,2)
       aux1122222111(1,1,1)=-(S21212222211111(1)*Z(jj,1))-S2122222221111
     &  1(1)*Z(jj,2)
       aux1122222211(1,1,1)=-(S21212222221111(1)*Z(jj,1))-S2122222222111
     &  1(1)*Z(jj,2)
       aux1122222221(1,1,1)=-(S21212222222111(1)*Z(jj,1))-S2122222222211
     &  1(1)*Z(jj,2)
       aux1122222222(1,1,1)=-(S21212222222211(1)*Z(jj,1))-S2122222222221
     &  1(1)*Z(jj,2)
       aux1122222222(2,1,1)=-(S21212222222221(1)*Z(jj,1))-S2122222222222
     &  1(1)*Z(jj,2)
       aux1122222222(2,2,1)=-(S21212222222222(1)*Z(jj,1))-S2122222222222
     &  2(1)*Z(jj,2)
       aux1122222222(2,2,2)=-(S21212222222222(2)*Z(jj,1))-S2122222222222
     &  2(2)*Z(jj,2)
       temp1111111111(1,1,1)=IX*(aux1111111111(1,1,1)+22*temp00101111111
     &  (1,1,1)*Z(jj,1))
       temp1121111111(1,1,1)=IX*(aux1121111111(1,1,1)+20*temp00102111111
     &  (1,1,1)*Z(jj,1)+2*temp00101111111(1,1,1)*Z(jj,2))
       temp1122111111(1,1,1)=IX*(aux1122111111(1,1,1)+18*temp00102211111
     &  (1,1,1)*Z(jj,1)+4*temp00102111111(1,1,1)*Z(jj,2))
       temp1122211111(1,1,1)=IX*(aux1122211111(1,1,1)+16*temp00102221111
     &  (1,1,1)*Z(jj,1)+6*temp00102211111(1,1,1)*Z(jj,2))
       temp1122221111(1,1,1)=IX*(aux1122221111(1,1,1)+14*temp00102222111
     &  (1,1,1)*Z(jj,1)+8*temp00102221111(1,1,1)*Z(jj,2))
       temp1122222111(1,1,1)=IX*(aux1122222111(1,1,1)+12*temp00102222211
     &  (1,1,1)*Z(jj,1)+10*temp00102222111(1,1,1)*Z(jj,2))
       temp1122222211(1,1,1)=IX*(aux1122222211(1,1,1)+10*temp00102222221
     &  (1,1,1)*Z(jj,1)+12*temp00102222211(1,1,1)*Z(jj,2))
       temp1122222221(1,1,1)=IX*(aux1122222221(1,1,1)+8*temp00102222222(
     &  1,1,1)*Z(jj,1)+14*temp00102222221(1,1,1)*Z(jj,2))
       temp1122222222(1,1,1)=IX*(aux1122222222(1,1,1)+6*temp00102222222(
     &  2,1,1)*Z(jj,1)+16*temp00102222222(1,1,1)*Z(jj,2))
       temp1122222222(2,1,1)=IX*(aux1122222222(2,1,1)+4*temp00102222222(
     &  2,2,1)*Z(jj,1)+18*temp00102222222(2,1,1)*Z(jj,2))
       temp1122222222(2,2,1)=IX*(aux1122222222(2,2,1)+2*temp00102222222(
     &  2,2,2)*Z(jj,1)+20*temp00102222222(2,2,1)*Z(jj,2))
       temp1122222222(2,2,2)=IX*(aux1122222222(2,2,2)+22*temp00102222222
     &  (2,2,2)*Z(jj,2))
       temp1111111111(1,1,2)=temp1121111111(1,1,1)
       temp1111111111(1,2,1)=temp1121111111(1,1,1)
       temp1111111111(1,2,2)=temp1122111111(1,1,1)
       temp1121111111(1,1,2)=temp1122111111(1,1,1)
       temp1121111111(1,2,1)=temp1122111111(1,1,1)
       temp1121111111(1,2,2)=temp1122211111(1,1,1)
       temp1122111111(1,1,2)=temp1122211111(1,1,1)
       temp1122111111(1,2,1)=temp1122211111(1,1,1)
       temp1122111111(1,2,2)=temp1122221111(1,1,1)
       temp1122211111(1,1,2)=temp1122221111(1,1,1)
       temp1122211111(1,2,1)=temp1122221111(1,1,1)
       temp1122211111(1,2,2)=temp1122222111(1,1,1)
       temp1122221111(1,1,2)=temp1122222111(1,1,1)
       temp1122221111(1,2,1)=temp1122222111(1,1,1)
       temp1122221111(1,2,2)=temp1122222211(1,1,1)
       temp1122222111(1,1,2)=temp1122222211(1,1,1)
       temp1122222111(1,2,1)=temp1122222211(1,1,1)
       temp1122222111(1,2,2)=temp1122222221(1,1,1)
       temp1122222211(1,1,2)=temp1122222221(1,1,1)
       temp1122222211(1,2,1)=temp1122222221(1,1,1)
       temp1122222211(1,2,2)=temp1122222222(1,1,1)
       temp1122222221(1,1,2)=temp1122222222(1,1,1)
       temp1122222221(1,2,1)=temp1122222222(1,1,1)
       temp1122222221(1,2,2)=temp1122222222(2,1,1)
       temp1122222222(1,1,2)=temp1122222222(2,1,1)
       temp1122222222(1,2,1)=temp1122222222(2,1,1)
       temp1122222222(1,2,2)=temp1122222222(2,2,1)
       temp1122222222(2,1,2)=temp1122222222(2,2,1)
c                Step2
       temp00000000001(1)=I26Z*(aux00000000001(1)+2*tempC30000000000*F(4
     &  )+F(3)*temp000000001(1)-det3*temp000000003(1,k,l))
       temp00000000001(2)=I26Z*(aux00000000001(2)+tempC30000000000*F(5)+
     &  F(3)*temp000000001(2)-det3*temp000000003(2,k,l))
       temp000000003(1,1,1)=I30Z*(aux000000003(1,1,1)+6*F(4)*temp0000000
     &  02(1,1)+F(3)*temp0000003(1,1,1)-det3*temp000000511(1,k,l)+24*tem
     &  p00000000001(1)*ZZ(k,1,l,1))
       temp000000003(2,1,1)=I30Z*(aux000000003(2,1,1)+F(5)*temp000000002
     &  (1,1)+4*F(4)*temp000000002(2,1)+F(3)*temp0000003(2,1,1)-det3*tem
     &  p000000521(1,k,l)+8*(temp00000000001(2)*ZZ(k,1,l,1)+temp00000000
     &  001(1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))))
       temp000000003(2,2,1)=I30Z*(aux000000003(2,2,1)+2*F(5)*temp0000000
     &  02(2,1)+2*F(4)*temp000000002(2,2)+F(3)*temp0000003(2,2,1)-det3*t
     &  emp000000522(1,k,l)+8*(temp00000000001(2)*(ZZ(k,1,l,2)+ZZ(k,2,l,
     &  1))+temp00000000001(1)*ZZ(k,2,l,2)))
       temp000000003(2,2,2)=I30Z*(aux000000003(2,2,2)+F(3)*temp0000003(2
     &  ,2,2)-det3*temp000000522(2,k,l)+24*temp00000000001(2)*ZZ(k,2,l,2
     &  )+temp000000002(2,2)*(6*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+12*r21*ZZ(
     &  k,2,l,2)))
       temp000000003(1,1,2)=temp000000003(2,1,1)
       temp000000003(1,2,1)=temp000000003(2,1,1)
       temp000000003(1,2,2)=temp000000003(2,2,1)
       temp000000003(2,1,2)=temp000000003(2,2,1)
       temp000000511(1,1,1)=I34Z*(aux000000511(1,1,1)+10*F(4)*temp000000
     &  41(1,1,1)+F(3)*temp0000511(1,1,1)-det3*temp000071111(1,k,l)+80*t
     &  emp000000003(1,1,1)*ZZ(k,1,l,1))
       temp000000521(1,1,1)=I34Z*(aux000000521(1,1,1)+F(5)*temp00000041(
     &  1,1,1)+8*F(4)*temp00000042(1,1,1)+F(3)*temp0000521(1,1,1)-det3*t
     &  emp000072111(1,k,l)+48*temp000000003(2,1,1)*ZZ(k,1,l,1)+16*temp0
     &  00000003(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp000000522(1,1,1)=I34Z*(aux000000522(1,1,1)+2*F(5)*temp0000004
     &  2(1,1,1)+6*F(4)*temp00000042(2,1,1)+F(3)*temp0000522(1,1,1)-det3
     &  *temp000072211(1,k,l)+24*(temp000000003(2,2,1)*ZZ(k,1,l,1)+temp0
     &  00000003(2,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))+8*temp000000003(1,1,1
     &  )*ZZ(k,2,l,2))
       temp000000522(2,1,1)=I34Z*(aux000000522(2,1,1)+4*F(4)*temp0000004
     &  2(2,2,1)+F(3)*temp0000522(2,1,1)-det3*temp000072221(1,k,l)+8*tem
     &  p000000003(2,2,2)*ZZ(k,1,l,1)+temp00000042(2,1,1)*(6*r10*(ZZ(k,1
     &  ,l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2))+24*(temp000000003(2,2,1)*
     &  (ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp000000003(2,1,1)*ZZ(k,2,l,2)))
       temp000000522(2,2,1)=I34Z*(aux000000522(2,2,1)+4*F(5)*temp0000004
     &  2(2,2,1)+2*F(4)*temp00000042(2,2,2)+F(3)*temp0000522(2,2,1)-det3
     &  *temp000072222(1,k,l)+16*temp000000003(2,2,2)*(ZZ(k,1,l,2)+ZZ(k,
     &  2,l,1))+48*temp000000003(2,2,1)*ZZ(k,2,l,2))
       temp000000522(2,2,2)=I34Z*(aux000000522(2,2,2)+F(3)*temp0000522(2
     &  ,2,2)-det3*temp000072222(2,k,l)+80*temp000000003(2,2,2)*ZZ(k,2,l
     &  ,2)+temp00000042(2,2,2)*(10*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+20*r21
     &  *ZZ(k,2,l,2)))
       temp000000511(1,1,2)=temp000000521(1,1,1)
       temp000000511(1,2,1)=temp000000521(1,1,1)
       temp000000511(1,2,2)=temp000000522(1,1,1)
       temp000000521(1,1,2)=temp000000522(1,1,1)
       temp000000521(1,2,1)=temp000000522(1,1,1)
       temp000000521(1,2,2)=temp000000522(2,1,1)
       temp000000522(1,1,2)=temp000000522(2,1,1)
       temp000000522(1,2,1)=temp000000522(2,1,1)
       temp000000522(1,2,2)=temp000000522(2,2,1)
       temp000000522(2,1,2)=temp000000522(2,2,1)
       temp000071111(1,1,1)=I38Z*(aux000071111(1,1,1)+14*F(4)*temp000061
     &  11(1,1,1)+F(3)*temp0071111(1,1,1)-det3*temp009111111(1,k,l)+168*
     &  temp000000511(1,1,1)*ZZ(k,1,l,1))
       temp000072111(1,1,1)=I38Z*(aux000072111(1,1,1)+F(3)*temp0072111(1
     &  ,1,1)-det3*temp009211111(1,k,l)+120*temp000000521(1,1,1)*ZZ(k,1,
     &  l,1)+2*(r10*temp00006111(1,1,1)+6*r21*temp00006211(1,1,1))*(ZZ(k
     &  ,1,l,2)+ZZ(k,2,l,1))+24*(r10*temp00006211(1,1,1)*ZZ(k,1,l,1)+tem
     &  p000000511(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))+4*r21*temp00006111(
     &  1,1,1)*ZZ(k,2,l,2))
       temp000072211(1,1,1)=I38Z*(aux000072211(1,1,1)+F(3)*temp0072211(1
     &  ,1,1)-det3*temp009221111(1,k,l)+80*temp000000522(1,1,1)*ZZ(k,1,l
     &  ,1)+20*r10*temp00006221(1,1,1)*ZZ(k,1,l,1)+40*temp000000521(1,1,
     &  1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+2*(2*r10*temp00006211(1,1,1)+5*r21*
     &  temp00006221(1,1,1))*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+8*(temp000000511(
     &  1,1,1)*ZZ(k,2,l,2)+r21*temp00006211(1,1,1)*ZZ(k,2,l,2)))
       temp000072221(1,1,1)=I38Z*(aux000072221(1,1,1)+F(3)*temp0072221(1
     &  ,1,1)-det3*temp009222111(1,k,l)+16*r10*temp00006222(1,1,1)*ZZ(k,
     &  1,l,1)+2*(3*r10*temp00006221(1,1,1)+4*r21*temp00006222(1,1,1))*(
     &  ZZ(k,1,l,2)+ZZ(k,2,l,1))+48*(temp000000522(2,1,1)*ZZ(k,1,l,1)+te
     &  mp000000522(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))+24*temp000000521(1
     &  ,1,1)*ZZ(k,2,l,2)+12*r21*temp00006221(1,1,1)*ZZ(k,2,l,2))
       temp000072222(1,1,1)=I38Z*(aux000072222(1,1,1)+F(3)*temp0072222(1
     &  ,1,1)-det3*temp009222211(1,k,l)+24*temp000000522(2,2,1)*ZZ(k,1,l
     &  ,1)+12*r10*temp00006222(2,1,1)*ZZ(k,1,l,1)+2*(4*r10*temp00006222
     &  (1,1,1)+3*r21*temp00006222(2,1,1))*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+16*
     &  r21*temp00006222(1,1,1)*ZZ(k,2,l,2)+48*(temp000000522(2,1,1)*(ZZ
     &  (k,1,l,2)+ZZ(k,2,l,1))+temp000000522(1,1,1)*ZZ(k,2,l,2)))
       temp000072222(2,1,1)=I38Z*(aux000072222(2,1,1)+F(3)*temp0072222(2
     &  ,1,1)-det3*temp009222221(1,k,l)+8*(temp000000522(2,2,2)*ZZ(k,1,l
     &  ,1)+r10*temp00006222(2,2,1)*ZZ(k,1,l,1))+40*temp000000522(2,2,1)
     &  *(ZZ(k,1,l,2)+ZZ(k,2,l,1))+2*(5*r10*temp00006222(2,1,1)+2*r21*te
     &  mp00006222(2,2,1))*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+80*temp000000522(2,
     &  1,1)*ZZ(k,2,l,2)+20*r21*temp00006222(2,1,1)*ZZ(k,2,l,2))
       temp000072222(2,2,1)=I38Z*(aux000072222(2,2,1)+F(3)*temp0072222(2
     &  ,2,1)-det3*temp009222222(1,k,l)+4*r10*temp00006222(2,2,2)*ZZ(k,1
     &  ,l,1)+2*(6*r10*temp00006222(2,2,1)+r21*temp00006222(2,2,2))*(ZZ(
     &  k,1,l,2)+ZZ(k,2,l,1))+120*temp000000522(2,2,1)*ZZ(k,2,l,2)+24*(t
     &  emp000000522(2,2,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+r21*temp00006222(2
     &  ,2,1)*ZZ(k,2,l,2)))
       temp000072222(2,2,2)=I38Z*(aux000072222(2,2,2)+F(3)*temp0072222(2
     &  ,2,2)-det3*temp009222222(2,k,l)+168*temp000000522(2,2,2)*ZZ(k,2,
     &  l,2)+14*temp00006222(2,2,2)*(r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+2*r21
     &  *ZZ(k,2,l,2)))
       temp000071111(1,1,2)=temp000072111(1,1,1)
       temp000071111(1,2,1)=temp000072111(1,1,1)
       temp000071111(1,2,2)=temp000072211(1,1,1)
       temp000072111(1,1,2)=temp000072211(1,1,1)
       temp000072111(1,2,1)=temp000072211(1,1,1)
       temp000072111(1,2,2)=temp000072221(1,1,1)
       temp000072211(1,1,2)=temp000072221(1,1,1)
       temp000072211(1,2,1)=temp000072221(1,1,1)
       temp000072211(1,2,2)=temp000072222(1,1,1)
       temp000072221(1,1,2)=temp000072222(1,1,1)
       temp000072221(1,2,1)=temp000072222(1,1,1)
       temp000072221(1,2,2)=temp000072222(2,1,1)
       temp000072222(1,1,2)=temp000072222(2,1,1)
       temp000072222(1,2,1)=temp000072222(2,1,1)
       temp000072222(1,2,2)=temp000072222(2,2,1)
       temp000072222(2,1,2)=temp000072222(2,2,1)
       temp009111111(1,1,1)=I42Z*(aux009111111(1,1,1)+18*F(4)*temp008111
     &  11(1,1,1)-det3*temp1111111111(1,k,l)+F(3)*temp9111111(1,1,1)+288
     &  *temp000071111(1,1,1)*ZZ(k,1,l,1))
       temp009211111(1,1,1)=I42Z*(aux009211111(1,1,1)+F(5)*temp00811111(
     &  1,1,1)+16*F(4)*temp00821111(1,1,1)-det3*temp1121111111(1,k,l)+F(
     &  3)*temp9211111(1,1,1)+224*temp000072111(1,1,1)*ZZ(k,1,l,1)+32*te
     &  mp000071111(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp009221111(1,1,1)=I42Z*(aux009221111(1,1,1)+2*F(5)*temp0082111
     &  1(1,1,1)+14*F(4)*temp00822111(1,1,1)-det3*temp1122111111(1,k,l)+
     &  F(3)*temp9221111(1,1,1)+168*temp000072211(1,1,1)*ZZ(k,1,l,1)+56*
     &  temp000072111(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+8*temp000071111(1
     &  ,1,1)*ZZ(k,2,l,2))
       temp009222111(1,1,1)=I42Z*(aux009222111(1,1,1)+12*F(4)*temp008222
     &  11(1,1,1)-det3*temp1122211111(1,k,l)+F(3)*temp9222111(1,1,1)+120
     &  *temp000072221(1,1,1)*ZZ(k,1,l,1)+72*temp000072211(1,1,1)*(ZZ(k,
     &  1,l,2)+ZZ(k,2,l,1))+24*temp000072111(1,1,1)*ZZ(k,2,l,2)+temp0082
     &  2111(1,1,1)*(6*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2))
     &  )
       temp009222211(1,1,1)=I42Z*(aux009222211(1,1,1)+4*F(5)*temp0082221
     &  1(1,1,1)+10*F(4)*temp00822221(1,1,1)-det3*temp1122221111(1,k,l)+
     &  F(3)*temp9222211(1,1,1)+80*(temp000072222(1,1,1)*ZZ(k,1,l,1)+tem
     &  p000072221(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))+48*temp000072211(1,
     &  1,1)*ZZ(k,2,l,2))
       temp009222221(1,1,1)=I42Z*(aux009222221(1,1,1)+8*F(4)*temp0082222
     &  2(1,1,1)-det3*temp1122222111(1,k,l)+F(3)*temp9222221(1,1,1)+48*t
     &  emp000072222(2,1,1)*ZZ(k,1,l,1)+temp00822221(1,1,1)*(10*r10*(ZZ(
     &  k,1,l,2)+ZZ(k,2,l,1))+20*r21*ZZ(k,2,l,2))+80*(temp000072222(1,1,
     &  1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp000072221(1,1,1)*ZZ(k,2,l,2)))
       temp009222222(1,1,1)=I42Z*(aux009222222(1,1,1)+6*F(4)*temp0082222
     &  2(2,1,1)-det3*temp1122222211(1,k,l)+F(3)*temp9222222(1,1,1)+24*t
     &  emp000072222(2,2,1)*ZZ(k,1,l,1)+72*temp000072222(2,1,1)*(ZZ(k,1,
     &  l,2)+ZZ(k,2,l,1))+120*temp000072222(1,1,1)*ZZ(k,2,l,2)+temp00822
     &  222(1,1,1)*(12*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+24*r21*ZZ(k,2,l,2))
     &  )
       temp009222222(2,1,1)=I42Z*(aux009222222(2,1,1)+4*F(4)*temp0082222
     &  2(2,2,1)-det3*temp1122222221(1,k,l)+F(3)*temp9222222(2,1,1)+8*te
     &  mp000072222(2,2,2)*ZZ(k,1,l,1)+56*temp000072222(2,2,1)*(ZZ(k,1,l
     &  ,2)+ZZ(k,2,l,1))+168*temp000072222(2,1,1)*ZZ(k,2,l,2)+temp008222
     &  22(2,1,1)*(14*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+28*r21*ZZ(k,2,l,2)))
       temp009222222(2,2,1)=I42Z*(aux009222222(2,2,1)+8*F(5)*temp0082222
     &  2(2,2,1)+2*F(4)*temp00822222(2,2,2)-det3*temp1122222222(1,k,l)+F
     &  (3)*temp9222222(2,2,1)+32*temp000072222(2,2,2)*(ZZ(k,1,l,2)+ZZ(k
     &  ,2,l,1))+224*temp000072222(2,2,1)*ZZ(k,2,l,2))
       temp009222222(2,2,2)=I42Z*(aux009222222(2,2,2)-det3*temp112222222
     &  2(2,k,l)+F(3)*temp9222222(2,2,2)+288*temp000072222(2,2,2)*ZZ(k,2
     &  ,l,2)+temp00822222(2,2,2)*(18*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+36*r
     &  21*ZZ(k,2,l,2)))
       temp009111111(1,1,2)=temp009211111(1,1,1)
       temp009111111(1,2,1)=temp009211111(1,1,1)
       temp009111111(1,2,2)=temp009221111(1,1,1)
       temp009211111(1,1,2)=temp009221111(1,1,1)
       temp009211111(1,2,1)=temp009221111(1,1,1)
       temp009211111(1,2,2)=temp009222111(1,1,1)
       temp009221111(1,1,2)=temp009222111(1,1,1)
       temp009221111(1,2,1)=temp009222111(1,1,1)
       temp009221111(1,2,2)=temp009222211(1,1,1)
       temp009222111(1,1,2)=temp009222211(1,1,1)
       temp009222111(1,2,1)=temp009222211(1,1,1)
       temp009222111(1,2,2)=temp009222221(1,1,1)
       temp009222211(1,1,2)=temp009222221(1,1,1)
       temp009222211(1,2,1)=temp009222221(1,1,1)
       temp009222211(1,2,2)=temp009222222(1,1,1)
       temp009222221(1,1,2)=temp009222222(1,1,1)
       temp009222221(1,2,1)=temp009222222(1,1,1)
       temp009222221(1,2,2)=temp009222222(2,1,1)
       temp009222222(1,1,2)=temp009222222(2,1,1)
       temp009222222(1,2,1)=temp009222222(2,1,1)
       temp009222222(1,2,2)=temp009222222(2,2,1)
       temp009222222(2,1,2)=temp009222222(2,2,1)
       temp101111111(1,1,1)=IX*(aux101111111(1,1,1)+det3*temp1111111111(
     &  1,1,jj)+20*temp009111111(1,1,1)*Z(jj,1))
       temp102111111(1,1,1)=IX*(aux102111111(1,1,1)+det3*temp1121111111(
     &  1,1,jj)+18*temp009211111(1,1,1)*Z(jj,1)+2*temp009111111(1,1,1)*Z
     &  (jj,2))
       temp102211111(1,1,1)=IX*(aux102211111(1,1,1)+det3*temp1122111111(
     &  1,1,jj)+16*temp009221111(1,1,1)*Z(jj,1)+4*temp009211111(1,1,1)*Z
     &  (jj,2))
       temp102221111(1,1,1)=IX*(aux102221111(1,1,1)+det3*temp1122211111(
     &  1,1,jj)+14*temp009222111(1,1,1)*Z(jj,1)+6*temp009221111(1,1,1)*Z
     &  (jj,2))
       temp102222111(1,1,1)=IX*(aux102222111(1,1,1)+det3*temp1122221111(
     &  1,1,jj)+12*temp009222211(1,1,1)*Z(jj,1)+8*temp009222111(1,1,1)*Z
     &  (jj,2))
       temp102222211(1,1,1)=IX*(aux102222211(1,1,1)+det3*temp1122222111(
     &  1,1,jj)+10*(temp009222221(1,1,1)*Z(jj,1)+temp009222211(1,1,1)*Z(
     &  jj,2)))
       temp102222221(1,1,1)=IX*(aux102222221(1,1,1)+det3*temp1122222211(
     &  1,1,jj)+8*temp009222222(1,1,1)*Z(jj,1)+12*temp009222221(1,1,1)*Z
     &  (jj,2))
       temp102222222(1,1,1)=IX*(aux102222222(1,1,1)+det3*temp1122222221(
     &  1,1,jj)+6*temp009222222(2,1,1)*Z(jj,1)+14*temp009222222(1,1,1)*Z
     &  (jj,2))
       temp102222222(2,1,1)=IX*(aux102222222(2,1,1)+det3*temp1122222222(
     &  1,1,jj)+4*temp009222222(2,2,1)*Z(jj,1)+16*temp009222222(2,1,1)*Z
     &  (jj,2))
       temp102222222(2,2,1)=IX*(aux102222222(2,2,1)+det3*temp1122222222(
     &  2,1,jj)+2*temp009222222(2,2,2)*Z(jj,1)+18*temp009222222(2,2,1)*Z
     &  (jj,2))
       temp102222222(2,2,2)=IX*(aux102222222(2,2,2)+det3*temp1122222222(
     &  2,2,jj)+20*temp009222222(2,2,2)*Z(jj,2))
       temp101111111(1,1,2)=temp102111111(1,1,1)
       temp101111111(1,2,1)=temp102111111(1,1,1)
       temp101111111(1,2,2)=temp102211111(1,1,1)
       temp102111111(1,1,2)=temp102211111(1,1,1)
       temp102111111(1,2,1)=temp102211111(1,1,1)
       temp102111111(1,2,2)=temp102221111(1,1,1)
       temp102211111(1,1,2)=temp102221111(1,1,1)
       temp102211111(1,2,1)=temp102221111(1,1,1)
       temp102211111(1,2,2)=temp102222111(1,1,1)
       temp102221111(1,1,2)=temp102222111(1,1,1)
       temp102221111(1,2,1)=temp102222111(1,1,1)
       temp102221111(1,2,2)=temp102222211(1,1,1)
       temp102222111(1,1,2)=temp102222211(1,1,1)
       temp102222111(1,2,1)=temp102222211(1,1,1)
       temp102222111(1,2,2)=temp102222221(1,1,1)
       temp102222211(1,1,2)=temp102222221(1,1,1)
       temp102222211(1,2,1)=temp102222221(1,1,1)
       temp102222211(1,2,2)=temp102222222(1,1,1)
       temp102222221(1,1,2)=temp102222222(1,1,1)
       temp102222221(1,2,1)=temp102222222(1,1,1)
       temp102222221(1,2,2)=temp102222222(2,1,1)
       temp102222222(1,1,2)=temp102222222(2,1,1)
       temp102222222(1,2,1)=temp102222222(2,1,1)
       temp102222222(1,2,2)=temp102222222(2,2,1)
       temp102222222(2,1,2)=temp102222222(2,2,1)
c                Step3
       tempC30000000000=I22Z*(auxC30000000000+tempC300000000*F(3)-det3*t
     &  emp000000002(k,l))
       temp000000002(1,1)=I26Z*(aux000000002(1,1)+4*F(4)*temp000000001(1
     &  )+F(3)*temp0000002(1,1)-det3*temp00000041(1,k,l)+8*tempC30000000
     &  000*ZZ(k,1,l,1))
       temp000000002(2,1)=I26Z*(aux000000002(2,1)+F(5)*temp000000001(1)+
     &  2*F(4)*temp000000001(2)+F(3)*temp0000002(2,1)-det3*temp00000042(
     &  1,k,l)+4*tempC30000000000*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp000000002(2,2)=I26Z*(aux000000002(2,2)+2*F(5)*temp000000001(2
     &  )+F(3)*temp0000002(2,2)-det3*temp00000042(2,k,l)+8*tempC30000000
     &  000*ZZ(k,2,l,2))
       temp000000002(1,2)=temp000000002(2,1)
       temp00000041(1,1,1)=I30Z*(aux00000041(1,1,1)+8*F(4)*temp0000003(1
     &  ,1,1)+F(3)*temp000041(1,1,1)-det3*temp00006111(1,k,l)+48*temp000
     &  000002(1,1)*ZZ(k,1,l,1))
       temp00000042(1,1,1)=I30Z*(aux00000042(1,1,1)+F(5)*temp0000003(1,1
     &  ,1)+6*F(4)*temp0000003(2,1,1)+F(3)*temp000042(1,1,1)-det3*temp00
     &  006211(1,k,l)+24*temp000000002(2,1)*ZZ(k,1,l,1)+12*temp000000002
     &  (1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp00000042(2,1,1)=I30Z*(aux00000042(2,1,1)+2*F(5)*temp0000003(2
     &  ,1,1)+4*F(4)*temp0000003(2,2,1)+F(3)*temp000042(2,1,1)-det3*temp
     &  00006221(1,k,l)+16*temp000000002(2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+
     &  8*(temp000000002(2,2)*ZZ(k,1,l,1)+temp000000002(1,1)*ZZ(k,2,l,2)
     &  ))
       temp00000042(2,2,1)=I30Z*(aux00000042(2,2,1)+2*F(4)*temp0000003(2
     &  ,2,2)+F(3)*temp000042(2,2,1)-det3*temp00006222(1,k,l)+12*temp000
     &  000002(2,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+24*temp000000002(2,1)*ZZ(k
     &  ,2,l,2)+temp0000003(2,2,1)*(6*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+12*r
     &  21*ZZ(k,2,l,2)))
       temp00000042(2,2,2)=I30Z*(aux00000042(2,2,2)+4*F(5)*temp0000003(2
     &  ,2,2)+F(3)*temp000042(2,2,2)-det3*temp00006222(2,k,l)+48*temp000
     &  000002(2,2)*ZZ(k,2,l,2))
       temp00000041(1,1,2)=temp00000042(1,1,1)
       temp00000041(1,2,1)=temp00000042(1,1,1)
       temp00000041(1,2,2)=temp00000042(2,1,1)
       temp00000042(1,1,2)=temp00000042(2,1,1)
       temp00000042(1,2,1)=temp00000042(2,1,1)
       temp00000042(1,2,2)=temp00000042(2,2,1)
       temp00000042(2,1,2)=temp00000042(2,2,1)
       temp00006111(1,1,1)=I34Z*(aux00006111(1,1,1)+12*F(4)*temp0000511(
     &  1,1,1)+F(3)*temp006111(1,1,1)-det3*temp00811111(1,k,l)+120*temp0
     &  0000041(1,1,1)*ZZ(k,1,l,1))
       temp00006211(1,1,1)=I34Z*(aux00006211(1,1,1)+F(5)*temp0000511(1,1
     &  ,1)+10*F(4)*temp0000521(1,1,1)+F(3)*temp006211(1,1,1)-det3*temp0
     &  0821111(1,k,l)+80*temp00000042(1,1,1)*ZZ(k,1,l,1)+20*temp0000004
     &  1(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp00006221(1,1,1)=I34Z*(aux00006221(1,1,1)+2*F(5)*temp0000521(1
     &  ,1,1)+F(3)*temp006221(1,1,1)-det3*temp00822111(1,k,l)+48*temp000
     &  00042(2,1,1)*ZZ(k,1,l,1)+32*temp00000042(1,1,1)*(ZZ(k,1,l,2)+ZZ(
     &  k,2,l,1))+8*(F(4)*temp0000522(1,1,1)+temp00000041(1,1,1)*ZZ(k,2,
     &  l,2)))
       temp00006222(1,1,1)=I34Z*(aux00006222(1,1,1)+6*F(4)*temp0000522(2
     &  ,1,1)+F(3)*temp006222(1,1,1)-det3*temp00822211(1,k,l)+36*temp000
     &  00042(2,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp0000522(1,1,1)*(6*r10
     &  *(ZZ(k,1,l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2))+24*(temp00000042(
     &  2,2,1)*ZZ(k,1,l,1)+temp00000042(1,1,1)*ZZ(k,2,l,2)))
       temp00006222(2,1,1)=I34Z*(aux00006222(2,1,1)+4*F(5)*temp0000522(2
     &  ,1,1)+4*F(4)*temp0000522(2,2,1)+F(3)*temp006222(2,1,1)-det3*temp
     &  00822221(1,k,l)+8*temp00000042(2,2,2)*ZZ(k,1,l,1)+32*temp0000004
     &  2(2,2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+48*temp00000042(2,1,1)*ZZ(k,2
     &  ,l,2))
       temp00006222(2,2,1)=I34Z*(aux00006222(2,2,1)+2*F(4)*temp0000522(2
     &  ,2,2)+F(3)*temp006222(2,2,1)-det3*temp00822222(1,k,l)+20*temp000
     &  00042(2,2,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+80*temp00000042(2,2,1)*ZZ
     &  (k,2,l,2)+temp0000522(2,2,1)*(10*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+2
     &  0*r21*ZZ(k,2,l,2)))
       temp00006222(2,2,2)=I34Z*(aux00006222(2,2,2)+F(3)*temp006222(2,2,
     &  2)-det3*temp00822222(2,k,l)+120*temp00000042(2,2,2)*ZZ(k,2,l,2)+
     &  temp0000522(2,2,2)*(12*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+24*r21*ZZ(k
     &  ,2,l,2)))
       temp00006111(1,1,2)=temp00006211(1,1,1)
       temp00006111(1,2,1)=temp00006211(1,1,1)
       temp00006111(1,2,2)=temp00006221(1,1,1)
       temp00006211(1,1,2)=temp00006221(1,1,1)
       temp00006211(1,2,1)=temp00006221(1,1,1)
       temp00006211(1,2,2)=temp00006222(1,1,1)
       temp00006221(1,1,2)=temp00006222(1,1,1)
       temp00006221(1,2,1)=temp00006222(1,1,1)
       temp00006221(1,2,2)=temp00006222(2,1,1)
       temp00006222(1,1,2)=temp00006222(2,1,1)
       temp00006222(1,2,1)=temp00006222(2,1,1)
       temp00006222(1,2,2)=temp00006222(2,2,1)
       temp00006222(2,1,2)=temp00006222(2,2,1)
       temp00811111(1,1,1)=I38Z*(aux00811111(1,1,1)+16*F(4)*temp0071111(
     &  1,1,1)-det3*temp101111111(1,k,l)+F(3)*temp811111(1,1,1)+224*temp
     &  00006111(1,1,1)*ZZ(k,1,l,1))
       temp00821111(1,1,1)=I38Z*(aux00821111(1,1,1)+F(5)*temp0071111(1,1
     &  ,1)+14*F(4)*temp0072111(1,1,1)-det3*temp102111111(1,k,l)+F(3)*te
     &  mp821111(1,1,1)+168*temp00006211(1,1,1)*ZZ(k,1,l,1)+28*temp00006
     &  111(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp00822111(1,1,1)=I38Z*(aux00822111(1,1,1)+2*F(5)*temp0072111(1
     &  ,1,1)+12*F(4)*temp0072211(1,1,1)-det3*temp102211111(1,k,l)+F(3)*
     &  temp822111(1,1,1)+120*temp00006221(1,1,1)*ZZ(k,1,l,1)+48*temp000
     &  06211(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+8*temp00006111(1,1,1)*ZZ(
     &  k,2,l,2))
       temp00822211(1,1,1)=I38Z*(aux00822211(1,1,1)+10*F(4)*temp0072221(
     &  1,1,1)-det3*temp102221111(1,k,l)+F(3)*temp822211(1,1,1)+80*temp0
     &  0006222(1,1,1)*ZZ(k,1,l,1)+60*temp00006221(1,1,1)*(ZZ(k,1,l,2)+Z
     &  Z(k,2,l,1))+24*temp00006211(1,1,1)*ZZ(k,2,l,2)+temp0072211(1,1,1
     &  )*(6*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2)))
       temp00822221(1,1,1)=I38Z*(aux00822221(1,1,1)+4*F(5)*temp0072221(1
     &  ,1,1)+8*F(4)*temp0072222(1,1,1)-det3*temp102222111(1,k,l)+F(3)*t
     &  emp822221(1,1,1)+64*temp00006222(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)
     &  )+48*(temp00006222(2,1,1)*ZZ(k,1,l,1)+temp00006221(1,1,1)*ZZ(k,2
     &  ,l,2)))
       temp00822222(1,1,1)=I38Z*(aux00822222(1,1,1)+6*F(4)*temp0072222(2
     &  ,1,1)-det3*temp102222211(1,k,l)+F(3)*temp822222(1,1,1)+24*temp00
     &  006222(2,2,1)*ZZ(k,1,l,1)+60*temp00006222(2,1,1)*(ZZ(k,1,l,2)+ZZ
     &  (k,2,l,1))+80*temp00006222(1,1,1)*ZZ(k,2,l,2)+temp0072222(1,1,1)
     &  *(10*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+20*r21*ZZ(k,2,l,2)))
       temp00822222(2,1,1)=I38Z*(aux00822222(2,1,1)+4*F(4)*temp0072222(2
     &  ,2,1)-det3*temp102222221(1,k,l)+F(3)*temp822222(2,1,1)+8*temp000
     &  06222(2,2,2)*ZZ(k,1,l,1)+48*temp00006222(2,2,1)*(ZZ(k,1,l,2)+ZZ(
     &  k,2,l,1))+120*temp00006222(2,1,1)*ZZ(k,2,l,2)+temp0072222(2,1,1)
     &  *(12*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+24*r21*ZZ(k,2,l,2)))
       temp00822222(2,2,1)=I38Z*(aux00822222(2,2,1)+2*F(4)*temp0072222(2
     &  ,2,2)-det3*temp102222222(1,k,l)+F(3)*temp822222(2,2,1)+28*temp00
     &  006222(2,2,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+168*temp00006222(2,2,1)*
     &  ZZ(k,2,l,2)+temp0072222(2,2,1)*(14*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))
     &  +28*r21*ZZ(k,2,l,2)))
       temp00822222(2,2,2)=I38Z*(aux00822222(2,2,2)+8*F(5)*temp0072222(2
     &  ,2,2)-det3*temp102222222(2,k,l)+F(3)*temp822222(2,2,2)+224*temp0
     &  0006222(2,2,2)*ZZ(k,2,l,2))
       temp00811111(1,1,2)=temp00821111(1,1,1)
       temp00811111(1,2,1)=temp00821111(1,1,1)
       temp00811111(1,2,2)=temp00822111(1,1,1)
       temp00821111(1,1,2)=temp00822111(1,1,1)
       temp00821111(1,2,1)=temp00822111(1,1,1)
       temp00821111(1,2,2)=temp00822211(1,1,1)
       temp00822111(1,1,2)=temp00822211(1,1,1)
       temp00822111(1,2,1)=temp00822211(1,1,1)
       temp00822111(1,2,2)=temp00822221(1,1,1)
       temp00822211(1,1,2)=temp00822221(1,1,1)
       temp00822211(1,2,1)=temp00822221(1,1,1)
       temp00822211(1,2,2)=temp00822222(1,1,1)
       temp00822221(1,1,2)=temp00822222(1,1,1)
       temp00822221(1,2,1)=temp00822222(1,1,1)
       temp00822221(1,2,2)=temp00822222(2,1,1)
       temp00822222(1,1,2)=temp00822222(2,1,1)
       temp00822222(1,2,1)=temp00822222(2,1,1)
       temp00822222(1,2,2)=temp00822222(2,2,1)
       temp00822222(2,1,2)=temp00822222(2,2,1)
       temp9111111(1,1,1)=IX*(aux9111111(1,1,1)+det3*temp101111111(1,1,j
     &  j)+18*temp00811111(1,1,1)*Z(jj,1))
       temp9211111(1,1,1)=IX*(aux9211111(1,1,1)+det3*temp102111111(1,1,j
     &  j)+16*temp00821111(1,1,1)*Z(jj,1)+2*temp00811111(1,1,1)*Z(jj,2))
       temp9221111(1,1,1)=IX*(aux9221111(1,1,1)+det3*temp102211111(1,1,j
     &  j)+14*temp00822111(1,1,1)*Z(jj,1)+4*temp00821111(1,1,1)*Z(jj,2))
       temp9222111(1,1,1)=IX*(aux9222111(1,1,1)+det3*temp102221111(1,1,j
     &  j)+12*temp00822211(1,1,1)*Z(jj,1)+6*temp00822111(1,1,1)*Z(jj,2))
       temp9222211(1,1,1)=IX*(aux9222211(1,1,1)+det3*temp102222111(1,1,j
     &  j)+10*temp00822221(1,1,1)*Z(jj,1)+8*temp00822211(1,1,1)*Z(jj,2))
       temp9222221(1,1,1)=IX*(aux9222221(1,1,1)+det3*temp102222211(1,1,j
     &  j)+8*temp00822222(1,1,1)*Z(jj,1)+10*temp00822221(1,1,1)*Z(jj,2))
       temp9222222(1,1,1)=IX*(aux9222222(1,1,1)+det3*temp102222221(1,1,j
     &  j)+6*temp00822222(2,1,1)*Z(jj,1)+12*temp00822222(1,1,1)*Z(jj,2))
       temp9222222(2,1,1)=IX*(aux9222222(2,1,1)+det3*temp102222222(1,1,j
     &  j)+4*temp00822222(2,2,1)*Z(jj,1)+14*temp00822222(2,1,1)*Z(jj,2))
       temp9222222(2,2,1)=IX*(aux9222222(2,2,1)+det3*temp102222222(2,1,j
     &  j)+2*temp00822222(2,2,2)*Z(jj,1)+16*temp00822222(2,2,1)*Z(jj,2))
       temp9222222(2,2,2)=IX*(aux9222222(2,2,2)+det3*temp102222222(2,2,j
     &  j)+18*temp00822222(2,2,2)*Z(jj,2))
       temp9111111(1,1,2)=temp9211111(1,1,1)
       temp9111111(1,2,1)=temp9211111(1,1,1)
       temp9111111(1,2,2)=temp9221111(1,1,1)
       temp9211111(1,1,2)=temp9221111(1,1,1)
       temp9211111(1,2,1)=temp9221111(1,1,1)
       temp9211111(1,2,2)=temp9222111(1,1,1)
       temp9221111(1,1,2)=temp9222111(1,1,1)
       temp9221111(1,2,1)=temp9222111(1,1,1)
       temp9221111(1,2,2)=temp9222211(1,1,1)
       temp9222111(1,1,2)=temp9222211(1,1,1)
       temp9222111(1,2,1)=temp9222211(1,1,1)
       temp9222111(1,2,2)=temp9222221(1,1,1)
       temp9222211(1,1,2)=temp9222221(1,1,1)
       temp9222211(1,2,1)=temp9222221(1,1,1)
       temp9222211(1,2,2)=temp9222222(1,1,1)
       temp9222221(1,1,2)=temp9222222(1,1,1)
       temp9222221(1,2,1)=temp9222222(1,1,1)
       temp9222221(1,2,2)=temp9222222(2,1,1)
       temp9222222(1,1,2)=temp9222222(2,1,1)
       temp9222222(1,2,1)=temp9222222(2,1,1)
       temp9222222(1,2,2)=temp9222222(2,2,1)
       temp9222222(2,1,2)=temp9222222(2,2,1)
c                Step4
       temp000000001(1)=I22Z*(aux000000001(1)+2*tempC300000000*F(4)+F(3)
     &  *temp0000001(1)-det3*temp0000003(1,k,l))
       temp000000001(2)=I22Z*(aux000000001(2)+tempC300000000*F(5)+F(3)*t
     &  emp0000001(2)-det3*temp0000003(2,k,l))
       temp0000003(1,1,1)=I26Z*(aux0000003(1,1,1)+6*F(4)*temp0000002(1,1
     &  )+F(3)*temp00003(1,1,1)-det3*temp0000511(1,k,l)+24*temp000000001
     &  (1)*ZZ(k,1,l,1))
       temp0000003(2,1,1)=I26Z*(aux0000003(2,1,1)+F(5)*temp0000002(1,1)+
     &  4*F(4)*temp0000002(2,1)+F(3)*temp00003(2,1,1)-det3*temp0000521(1
     &  ,k,l)+8*(temp000000001(2)*ZZ(k,1,l,1)+temp000000001(1)*(ZZ(k,1,l
     &  ,2)+ZZ(k,2,l,1))))
       temp0000003(2,2,1)=I26Z*(aux0000003(2,2,1)+2*F(5)*temp0000002(2,1
     &  )+2*F(4)*temp0000002(2,2)+F(3)*temp00003(2,2,1)-det3*temp0000522
     &  (1,k,l)+8*(temp000000001(2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp000000
     &  001(1)*ZZ(k,2,l,2)))
       temp0000003(2,2,2)=I26Z*(aux0000003(2,2,2)+F(3)*temp00003(2,2,2)-
     &  det3*temp0000522(2,k,l)+24*temp000000001(2)*ZZ(k,2,l,2)+temp0000
     &  002(2,2)*(6*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2)))
       temp0000003(1,1,2)=temp0000003(2,1,1)
       temp0000003(1,2,1)=temp0000003(2,1,1)
       temp0000003(1,2,2)=temp0000003(2,2,1)
       temp0000003(2,1,2)=temp0000003(2,2,1)
       temp0000511(1,1,1)=I30Z*(aux0000511(1,1,1)+10*F(4)*temp000041(1,1
     &  ,1)+F(3)*temp00511(1,1,1)-det3*temp0071111(1,k,l)+80*temp0000003
     &  (1,1,1)*ZZ(k,1,l,1))
       temp0000521(1,1,1)=I30Z*(aux0000521(1,1,1)+F(5)*temp000041(1,1,1)
     &  +8*F(4)*temp000042(1,1,1)+F(3)*temp00521(1,1,1)-det3*temp0072111
     &  (1,k,l)+48*temp0000003(2,1,1)*ZZ(k,1,l,1)+16*temp0000003(1,1,1)*
     &  (ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp0000522(1,1,1)=I30Z*(aux0000522(1,1,1)+2*F(5)*temp000042(1,1,
     &  1)+6*F(4)*temp000042(2,1,1)+F(3)*temp00522(1,1,1)-det3*temp00722
     &  11(1,k,l)+24*(temp0000003(2,2,1)*ZZ(k,1,l,1)+temp0000003(2,1,1)*
     &  (ZZ(k,1,l,2)+ZZ(k,2,l,1)))+8*temp0000003(1,1,1)*ZZ(k,2,l,2))
       temp0000522(2,1,1)=I30Z*(aux0000522(2,1,1)+4*F(4)*temp000042(2,2,
     &  1)+F(3)*temp00522(2,1,1)-det3*temp0072221(1,k,l)+8*temp0000003(2
     &  ,2,2)*ZZ(k,1,l,1)+temp000042(2,1,1)*(6*r10*(ZZ(k,1,l,2)+ZZ(k,2,l
     &  ,1))+12*r21*ZZ(k,2,l,2))+24*(temp0000003(2,2,1)*(ZZ(k,1,l,2)+ZZ(
     &  k,2,l,1))+temp0000003(2,1,1)*ZZ(k,2,l,2)))
       temp0000522(2,2,1)=I30Z*(aux0000522(2,2,1)+4*F(5)*temp000042(2,2,
     &  1)+2*F(4)*temp000042(2,2,2)+F(3)*temp00522(2,2,1)-det3*temp00722
     &  22(1,k,l)+16*temp0000003(2,2,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+48*tem
     &  p0000003(2,2,1)*ZZ(k,2,l,2))
       temp0000522(2,2,2)=I30Z*(aux0000522(2,2,2)+F(3)*temp00522(2,2,2)-
     &  det3*temp0072222(2,k,l)+80*temp0000003(2,2,2)*ZZ(k,2,l,2)+temp00
     &  0042(2,2,2)*(10*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+20*r21*ZZ(k,2,l,2)
     &  ))
       temp0000511(1,1,2)=temp0000521(1,1,1)
       temp0000511(1,2,1)=temp0000521(1,1,1)
       temp0000511(1,2,2)=temp0000522(1,1,1)
       temp0000521(1,1,2)=temp0000522(1,1,1)
       temp0000521(1,2,1)=temp0000522(1,1,1)
       temp0000521(1,2,2)=temp0000522(2,1,1)
       temp0000522(1,1,2)=temp0000522(2,1,1)
       temp0000522(1,2,1)=temp0000522(2,1,1)
       temp0000522(1,2,2)=temp0000522(2,2,1)
       temp0000522(2,1,2)=temp0000522(2,2,1)
       temp0071111(1,1,1)=I34Z*(aux0071111(1,1,1)+14*F(4)*temp006111(1,1
     &  ,1)+F(3)*temp71111(1,1,1)-det3*temp9111111(1,k,l)+168*temp000051
     &  1(1,1,1)*ZZ(k,1,l,1))
       temp0072111(1,1,1)=I34Z*(aux0072111(1,1,1)+F(5)*temp006111(1,1,1)
     &  +12*F(4)*temp006211(1,1,1)+F(3)*temp72111(1,1,1)-det3*temp921111
     &  1(1,k,l)+120*temp0000521(1,1,1)*ZZ(k,1,l,1)+24*temp0000511(1,1,1
     &  )*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp0072211(1,1,1)=I34Z*(aux0072211(1,1,1)+2*F(5)*temp006211(1,1,
     &  1)+10*F(4)*temp006221(1,1,1)+F(3)*temp72211(1,1,1)-det3*temp9221
     &  111(1,k,l)+80*temp0000522(1,1,1)*ZZ(k,1,l,1)+40*temp0000521(1,1,
     &  1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+8*temp0000511(1,1,1)*ZZ(k,2,l,2))
       temp0072221(1,1,1)=I34Z*(aux0072221(1,1,1)+8*F(4)*temp006222(1,1,
     &  1)+F(3)*temp72221(1,1,1)-det3*temp9222111(1,k,l)+48*(temp0000522
     &  (2,1,1)*ZZ(k,1,l,1)+temp0000522(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))
     &  )+24*temp0000521(1,1,1)*ZZ(k,2,l,2)+temp006221(1,1,1)*(6*r10*(ZZ
     &  (k,1,l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2)))
       temp0072222(1,1,1)=I34Z*(aux0072222(1,1,1)+4*F(5)*temp006222(1,1,
     &  1)+6*F(4)*temp006222(2,1,1)+F(3)*temp72222(1,1,1)-det3*temp92222
     &  11(1,k,l)+24*temp0000522(2,2,1)*ZZ(k,1,l,1)+48*(temp0000522(2,1,
     &  1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp0000522(1,1,1)*ZZ(k,2,l,2)))
       temp0072222(2,1,1)=I34Z*(aux0072222(2,1,1)+4*F(4)*temp006222(2,2,
     &  1)+F(3)*temp72222(2,1,1)-det3*temp9222221(1,k,l)+8*temp0000522(2
     &  ,2,2)*ZZ(k,1,l,1)+40*temp0000522(2,2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)
     &  )+80*temp0000522(2,1,1)*ZZ(k,2,l,2)+temp006222(2,1,1)*(10*r10*(Z
     &  Z(k,1,l,2)+ZZ(k,2,l,1))+20*r21*ZZ(k,2,l,2)))
       temp0072222(2,2,1)=I34Z*(aux0072222(2,2,1)+2*F(4)*temp006222(2,2,
     &  2)+F(3)*temp72222(2,2,1)-det3*temp9222222(1,k,l)+24*temp0000522(
     &  2,2,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+120*temp0000522(2,2,1)*ZZ(k,2,l
     &  ,2)+temp006222(2,2,1)*(12*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+24*r21*Z
     &  Z(k,2,l,2)))
       temp0072222(2,2,2)=I34Z*(aux0072222(2,2,2)+F(3)*temp72222(2,2,2)-
     &  det3*temp9222222(2,k,l)+168*temp0000522(2,2,2)*ZZ(k,2,l,2)+temp0
     &  06222(2,2,2)*(14*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+28*r21*ZZ(k,2,l,2
     &  )))
       temp0071111(1,1,2)=temp0072111(1,1,1)
       temp0071111(1,2,1)=temp0072111(1,1,1)
       temp0071111(1,2,2)=temp0072211(1,1,1)
       temp0072111(1,1,2)=temp0072211(1,1,1)
       temp0072111(1,2,1)=temp0072211(1,1,1)
       temp0072111(1,2,2)=temp0072221(1,1,1)
       temp0072211(1,1,2)=temp0072221(1,1,1)
       temp0072211(1,2,1)=temp0072221(1,1,1)
       temp0072211(1,2,2)=temp0072222(1,1,1)
       temp0072221(1,1,2)=temp0072222(1,1,1)
       temp0072221(1,2,1)=temp0072222(1,1,1)
       temp0072221(1,2,2)=temp0072222(2,1,1)
       temp0072222(1,1,2)=temp0072222(2,1,1)
       temp0072222(1,2,1)=temp0072222(2,1,1)
       temp0072222(1,2,2)=temp0072222(2,2,1)
       temp0072222(2,1,2)=temp0072222(2,2,1)
       temp811111(1,1,1)=IX*(aux811111(1,1,1)+det3*temp9111111(1,1,jj)+1
     &  6*temp0071111(1,1,1)*Z(jj,1))
       temp821111(1,1,1)=IX*(aux821111(1,1,1)+det3*temp9211111(1,1,jj)+1
     &  4*temp0072111(1,1,1)*Z(jj,1)+2*temp0071111(1,1,1)*Z(jj,2))
       temp822111(1,1,1)=IX*(aux822111(1,1,1)+det3*temp9221111(1,1,jj)+1
     &  2*temp0072211(1,1,1)*Z(jj,1)+4*temp0072111(1,1,1)*Z(jj,2))
       temp822211(1,1,1)=IX*(aux822211(1,1,1)+det3*temp9222111(1,1,jj)+1
     &  0*temp0072221(1,1,1)*Z(jj,1)+6*temp0072211(1,1,1)*Z(jj,2))
       temp822221(1,1,1)=IX*(aux822221(1,1,1)+det3*temp9222211(1,1,jj)+8
     &  *(temp0072222(1,1,1)*Z(jj,1)+temp0072221(1,1,1)*Z(jj,2)))
       temp822222(1,1,1)=IX*(aux822222(1,1,1)+det3*temp9222221(1,1,jj)+6
     &  *temp0072222(2,1,1)*Z(jj,1)+10*temp0072222(1,1,1)*Z(jj,2))
       temp822222(2,1,1)=IX*(aux822222(2,1,1)+det3*temp9222222(1,1,jj)+4
     &  *temp0072222(2,2,1)*Z(jj,1)+12*temp0072222(2,1,1)*Z(jj,2))
       temp822222(2,2,1)=IX*(aux822222(2,2,1)+det3*temp9222222(2,1,jj)+2
     &  *temp0072222(2,2,2)*Z(jj,1)+14*temp0072222(2,2,1)*Z(jj,2))
       temp822222(2,2,2)=IX*(aux822222(2,2,2)+det3*temp9222222(2,2,jj)+1
     &  6*temp0072222(2,2,2)*Z(jj,2))
       temp811111(1,1,2)=temp821111(1,1,1)
       temp811111(1,2,1)=temp821111(1,1,1)
       temp811111(1,2,2)=temp822111(1,1,1)
       temp821111(1,1,2)=temp822111(1,1,1)
       temp821111(1,2,1)=temp822111(1,1,1)
       temp821111(1,2,2)=temp822211(1,1,1)
       temp822111(1,1,2)=temp822211(1,1,1)
       temp822111(1,2,1)=temp822211(1,1,1)
       temp822111(1,2,2)=temp822221(1,1,1)
       temp822211(1,1,2)=temp822221(1,1,1)
       temp822211(1,2,1)=temp822221(1,1,1)
       temp822211(1,2,2)=temp822222(1,1,1)
       temp822221(1,1,2)=temp822222(1,1,1)
       temp822221(1,2,1)=temp822222(1,1,1)
       temp822221(1,2,2)=temp822222(2,1,1)
       temp822222(1,1,2)=temp822222(2,1,1)
       temp822222(1,2,1)=temp822222(2,1,1)
       temp822222(1,2,2)=temp822222(2,2,1)
       temp822222(2,1,2)=temp822222(2,2,1)
c                Step5
       tempC300000000=I18Z*(auxC300000000+tempC3000000*F(3)-det3*temp000
     &  0002(k,l))
       temp0000002(1,1)=I22Z*(aux0000002(1,1)+4*F(4)*temp0000001(1)+F(3)
     &  *temp00002(1,1)-det3*temp000041(1,k,l)+8*tempC300000000*ZZ(k,1,l
     &  ,1))
       temp0000002(2,1)=I22Z*(aux0000002(2,1)+F(5)*temp0000001(1)+2*F(4)
     &  *temp0000001(2)+F(3)*temp00002(2,1)-det3*temp000042(1,k,l)+4*tem
     &  pC300000000*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp0000002(2,2)=I22Z*(aux0000002(2,2)+2*F(5)*temp0000001(2)+F(3)
     &  *temp00002(2,2)-det3*temp000042(2,k,l)+8*tempC300000000*ZZ(k,2,l
     &  ,2))
       temp0000002(1,2)=temp0000002(2,1)
       temp000041(1,1,1)=I26Z*(aux000041(1,1,1)+8*F(4)*temp00003(1,1,1)+
     &  F(3)*temp0041(1,1,1)-det3*temp006111(1,k,l)+48*temp0000002(1,1)*
     &  ZZ(k,1,l,1))
       temp000042(1,1,1)=I26Z*(aux000042(1,1,1)+F(5)*temp00003(1,1,1)+6*
     &  F(4)*temp00003(2,1,1)+F(3)*temp0042(1,1,1)-det3*temp006211(1,k,l
     &  )+24*temp0000002(2,1)*ZZ(k,1,l,1)+12*temp0000002(1,1)*(ZZ(k,1,l,
     &  2)+ZZ(k,2,l,1)))
       temp000042(2,1,1)=I26Z*(aux000042(2,1,1)+2*F(5)*temp00003(2,1,1)+
     &  4*F(4)*temp00003(2,2,1)+F(3)*temp0042(2,1,1)-det3*temp006221(1,k
     &  ,l)+16*temp0000002(2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+8*(temp0000002
     &  (2,2)*ZZ(k,1,l,1)+temp0000002(1,1)*ZZ(k,2,l,2)))
       temp000042(2,2,1)=I26Z*(aux000042(2,2,1)+2*F(4)*temp00003(2,2,2)+
     &  F(3)*temp0042(2,2,1)-det3*temp006222(1,k,l)+12*temp0000002(2,2)*
     &  (ZZ(k,1,l,2)+ZZ(k,2,l,1))+24*temp0000002(2,1)*ZZ(k,2,l,2)+temp00
     &  003(2,2,1)*(6*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2)))
       temp000042(2,2,2)=I26Z*(aux000042(2,2,2)+4*F(5)*temp00003(2,2,2)+
     &  F(3)*temp0042(2,2,2)-det3*temp006222(2,k,l)+48*temp0000002(2,2)*
     &  ZZ(k,2,l,2))
       temp000041(1,1,2)=temp000042(1,1,1)
       temp000041(1,2,1)=temp000042(1,1,1)
       temp000041(1,2,2)=temp000042(2,1,1)
       temp000042(1,1,2)=temp000042(2,1,1)
       temp000042(1,2,1)=temp000042(2,1,1)
       temp000042(1,2,2)=temp000042(2,2,1)
       temp000042(2,1,2)=temp000042(2,2,1)
       temp006111(1,1,1)=I30Z*(aux006111(1,1,1)+12*F(4)*temp00511(1,1,1)
     &  +F(3)*temp6111(1,1,1)-det3*temp811111(1,k,l)+120*temp000041(1,1,
     &  1)*ZZ(k,1,l,1))
       temp006211(1,1,1)=I30Z*(aux006211(1,1,1)+F(5)*temp00511(1,1,1)+10
     &  *F(4)*temp00521(1,1,1)+F(3)*temp6211(1,1,1)-det3*temp821111(1,k,
     &  l)+80*temp000042(1,1,1)*ZZ(k,1,l,1)+20*temp000041(1,1,1)*(ZZ(k,1
     &  ,l,2)+ZZ(k,2,l,1)))
       temp006221(1,1,1)=I30Z*(aux006221(1,1,1)+2*F(5)*temp00521(1,1,1)+
     &  F(3)*temp6221(1,1,1)-det3*temp822111(1,k,l)+48*temp000042(2,1,1)
     &  *ZZ(k,1,l,1)+32*temp000042(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+8*(F
     &  (4)*temp00522(1,1,1)+temp000041(1,1,1)*ZZ(k,2,l,2)))
       temp006222(1,1,1)=I30Z*(aux006222(1,1,1)+6*F(4)*temp00522(2,1,1)+
     &  F(3)*temp6222(1,1,1)-det3*temp822211(1,k,l)+36*temp000042(2,1,1)
     &  *(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp00522(1,1,1)*(6*r10*(ZZ(k,1,l,2)+
     &  ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2))+24*(temp000042(2,2,1)*ZZ(k,1,l,
     &  1)+temp000042(1,1,1)*ZZ(k,2,l,2)))
       temp006222(2,1,1)=I30Z*(aux006222(2,1,1)+4*F(5)*temp00522(2,1,1)+
     &  4*F(4)*temp00522(2,2,1)+F(3)*temp6222(2,1,1)-det3*temp822221(1,k
     &  ,l)+8*temp000042(2,2,2)*ZZ(k,1,l,1)+32*temp000042(2,2,1)*(ZZ(k,1
     &  ,l,2)+ZZ(k,2,l,1))+48*temp000042(2,1,1)*ZZ(k,2,l,2))
       temp006222(2,2,1)=I30Z*(aux006222(2,2,1)+2*F(4)*temp00522(2,2,2)+
     &  F(3)*temp6222(2,2,1)-det3*temp822222(1,k,l)+20*temp000042(2,2,2)
     &  *(ZZ(k,1,l,2)+ZZ(k,2,l,1))+80*temp000042(2,2,1)*ZZ(k,2,l,2)+temp
     &  00522(2,2,1)*(10*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+20*r21*ZZ(k,2,l,2
     &  )))
       temp006222(2,2,2)=I30Z*(aux006222(2,2,2)+F(3)*temp6222(2,2,2)-det
     &  3*temp822222(2,k,l)+120*temp000042(2,2,2)*ZZ(k,2,l,2)+temp00522(
     &  2,2,2)*(12*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+24*r21*ZZ(k,2,l,2)))
       temp006111(1,1,2)=temp006211(1,1,1)
       temp006111(1,2,1)=temp006211(1,1,1)
       temp006111(1,2,2)=temp006221(1,1,1)
       temp006211(1,1,2)=temp006221(1,1,1)
       temp006211(1,2,1)=temp006221(1,1,1)
       temp006211(1,2,2)=temp006222(1,1,1)
       temp006221(1,1,2)=temp006222(1,1,1)
       temp006221(1,2,1)=temp006222(1,1,1)
       temp006221(1,2,2)=temp006222(2,1,1)
       temp006222(1,1,2)=temp006222(2,1,1)
       temp006222(1,2,1)=temp006222(2,1,1)
       temp006222(1,2,2)=temp006222(2,2,1)
       temp006222(2,1,2)=temp006222(2,2,1)
       temp71111(1,1,1)=IX*(aux71111(1,1,1)+det3*temp811111(1,1,jj)+14*t
     &  emp006111(1,1,1)*Z(jj,1))
       temp72111(1,1,1)=IX*(aux72111(1,1,1)+det3*temp821111(1,1,jj)+12*t
     &  emp006211(1,1,1)*Z(jj,1)+2*temp006111(1,1,1)*Z(jj,2))
       temp72211(1,1,1)=IX*(aux72211(1,1,1)+det3*temp822111(1,1,jj)+10*t
     &  emp006221(1,1,1)*Z(jj,1)+4*temp006211(1,1,1)*Z(jj,2))
       temp72221(1,1,1)=IX*(aux72221(1,1,1)+det3*temp822211(1,1,jj)+8*te
     &  mp006222(1,1,1)*Z(jj,1)+6*temp006221(1,1,1)*Z(jj,2))
       temp72222(1,1,1)=IX*(aux72222(1,1,1)+det3*temp822221(1,1,jj)+6*te
     &  mp006222(2,1,1)*Z(jj,1)+8*temp006222(1,1,1)*Z(jj,2))
       temp72222(2,1,1)=IX*(aux72222(2,1,1)+det3*temp822222(1,1,jj)+4*te
     &  mp006222(2,2,1)*Z(jj,1)+10*temp006222(2,1,1)*Z(jj,2))
       temp72222(2,2,1)=IX*(aux72222(2,2,1)+det3*temp822222(2,1,jj)+2*te
     &  mp006222(2,2,2)*Z(jj,1)+12*temp006222(2,2,1)*Z(jj,2))
       temp72222(2,2,2)=IX*(aux72222(2,2,2)+det3*temp822222(2,2,jj)+14*t
     &  emp006222(2,2,2)*Z(jj,2))
       temp71111(1,1,2)=temp72111(1,1,1)
       temp71111(1,2,1)=temp72111(1,1,1)
       temp71111(1,2,2)=temp72211(1,1,1)
       temp72111(1,1,2)=temp72211(1,1,1)
       temp72111(1,2,1)=temp72211(1,1,1)
       temp72111(1,2,2)=temp72221(1,1,1)
       temp72211(1,1,2)=temp72221(1,1,1)
       temp72211(1,2,1)=temp72221(1,1,1)
       temp72211(1,2,2)=temp72222(1,1,1)
       temp72221(1,1,2)=temp72222(1,1,1)
       temp72221(1,2,1)=temp72222(1,1,1)
       temp72221(1,2,2)=temp72222(2,1,1)
       temp72222(1,1,2)=temp72222(2,1,1)
       temp72222(1,2,1)=temp72222(2,1,1)
       temp72222(1,2,2)=temp72222(2,2,1)
       temp72222(2,1,2)=temp72222(2,2,1)
c                Step6
       temp0000001(1)=I18Z*(aux0000001(1)+2*tempC3000000*F(4)+F(3)*temp0
     &  0001(1)-det3*temp00003(1,k,l))
       temp0000001(2)=I18Z*(aux0000001(2)+tempC3000000*F(5)+F(3)*temp000
     &  01(2)-det3*temp00003(2,k,l))
       temp00003(1,1,1)=I22Z*(aux00003(1,1,1)+6*F(4)*temp00002(1,1)+F(3)
     &  *temp003(1,1,1)-det3*temp00511(1,k,l)+24*temp0000001(1)*ZZ(k,1,l
     &  ,1))
       temp00003(2,1,1)=I22Z*(aux00003(2,1,1)+F(5)*temp00002(1,1)+4*F(4)
     &  *temp00002(2,1)+F(3)*temp003(2,1,1)-det3*temp00521(1,k,l)+8*(tem
     &  p0000001(2)*ZZ(k,1,l,1)+temp0000001(1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))
     &  ))
       temp00003(2,2,1)=I22Z*(aux00003(2,2,1)+2*F(5)*temp00002(2,1)+2*F(
     &  4)*temp00002(2,2)+F(3)*temp003(2,2,1)-det3*temp00522(1,k,l)+8*(t
     &  emp0000001(2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp0000001(1)*ZZ(k,2,l,
     &  2)))
       temp00003(2,2,2)=I22Z*(aux00003(2,2,2)+F(3)*temp003(2,2,2)-det3*t
     &  emp00522(2,k,l)+24*temp0000001(2)*ZZ(k,2,l,2)+temp00002(2,2)*(6*
     &  r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2)))
       temp00003(1,1,2)=temp00003(2,1,1)
       temp00003(1,2,1)=temp00003(2,1,1)
       temp00003(1,2,2)=temp00003(2,2,1)
       temp00003(2,1,2)=temp00003(2,2,1)
       temp00511(1,1,1)=I26Z*(aux00511(1,1,1)+10*F(4)*temp0041(1,1,1)+F(
     &  3)*temp511(1,1,1)-det3*temp71111(1,k,l)+80*temp00003(1,1,1)*ZZ(k
     &  ,1,l,1))
       temp00521(1,1,1)=I26Z*(aux00521(1,1,1)+F(5)*temp0041(1,1,1)+8*F(4
     &  )*temp0042(1,1,1)+F(3)*temp521(1,1,1)-det3*temp72111(1,k,l)+48*t
     &  emp00003(2,1,1)*ZZ(k,1,l,1)+16*temp00003(1,1,1)*(ZZ(k,1,l,2)+ZZ(
     &  k,2,l,1)))
       temp00522(1,1,1)=I26Z*(aux00522(1,1,1)+2*F(5)*temp0042(1,1,1)+6*F
     &  (4)*temp0042(2,1,1)+F(3)*temp522(1,1,1)-det3*temp72211(1,k,l)+24
     &  *(temp00003(2,2,1)*ZZ(k,1,l,1)+temp00003(2,1,1)*(ZZ(k,1,l,2)+ZZ(
     &  k,2,l,1)))+8*temp00003(1,1,1)*ZZ(k,2,l,2))
       temp00522(2,1,1)=I26Z*(aux00522(2,1,1)+4*F(4)*temp0042(2,2,1)+F(3
     &  )*temp522(2,1,1)-det3*temp72221(1,k,l)+8*temp00003(2,2,2)*ZZ(k,1
     &  ,l,1)+temp0042(2,1,1)*(6*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+12*r21*ZZ
     &  (k,2,l,2))+24*(temp00003(2,2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp00
     &  003(2,1,1)*ZZ(k,2,l,2)))
       temp00522(2,2,1)=I26Z*(aux00522(2,2,1)+4*F(5)*temp0042(2,2,1)+2*F
     &  (4)*temp0042(2,2,2)+F(3)*temp522(2,2,1)-det3*temp72222(1,k,l)+16
     &  *temp00003(2,2,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+48*temp00003(2,2,1)*
     &  ZZ(k,2,l,2))
       temp00522(2,2,2)=I26Z*(aux00522(2,2,2)+F(3)*temp522(2,2,2)-det3*t
     &  emp72222(2,k,l)+80*temp00003(2,2,2)*ZZ(k,2,l,2)+temp0042(2,2,2)*
     &  (10*r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+20*r21*ZZ(k,2,l,2)))
       temp00511(1,1,2)=temp00521(1,1,1)
       temp00511(1,2,1)=temp00521(1,1,1)
       temp00511(1,2,2)=temp00522(1,1,1)
       temp00521(1,1,2)=temp00522(1,1,1)
       temp00521(1,2,1)=temp00522(1,1,1)
       temp00521(1,2,2)=temp00522(2,1,1)
       temp00522(1,1,2)=temp00522(2,1,1)
       temp00522(1,2,1)=temp00522(2,1,1)
       temp00522(1,2,2)=temp00522(2,2,1)
       temp00522(2,1,2)=temp00522(2,2,1)
       temp6111(1,1,1)=IX*(aux6111(1,1,1)+det3*temp71111(1,1,jj)+12*temp
     &  00511(1,1,1)*Z(jj,1))
       temp6211(1,1,1)=IX*(aux6211(1,1,1)+det3*temp72111(1,1,jj)+10*temp
     &  00521(1,1,1)*Z(jj,1)+2*temp00511(1,1,1)*Z(jj,2))
       temp6221(1,1,1)=IX*(aux6221(1,1,1)+det3*temp72211(1,1,jj)+8*temp0
     &  0522(1,1,1)*Z(jj,1)+4*temp00521(1,1,1)*Z(jj,2))
       temp6222(1,1,1)=IX*(aux6222(1,1,1)+det3*temp72221(1,1,jj)+6*(temp
     &  00522(2,1,1)*Z(jj,1)+temp00522(1,1,1)*Z(jj,2)))
       temp6222(2,1,1)=IX*(aux6222(2,1,1)+det3*temp72222(1,1,jj)+4*temp0
     &  0522(2,2,1)*Z(jj,1)+8*temp00522(2,1,1)*Z(jj,2))
       temp6222(2,2,1)=IX*(aux6222(2,2,1)+det3*temp72222(2,1,jj)+2*temp0
     &  0522(2,2,2)*Z(jj,1)+10*temp00522(2,2,1)*Z(jj,2))
       temp6222(2,2,2)=IX*(aux6222(2,2,2)+det3*temp72222(2,2,jj)+12*temp
     &  00522(2,2,2)*Z(jj,2))
       temp6111(1,1,2)=temp6211(1,1,1)
       temp6111(1,2,1)=temp6211(1,1,1)
       temp6111(1,2,2)=temp6221(1,1,1)
       temp6211(1,1,2)=temp6221(1,1,1)
       temp6211(1,2,1)=temp6221(1,1,1)
       temp6211(1,2,2)=temp6222(1,1,1)
       temp6221(1,1,2)=temp6222(1,1,1)
       temp6221(1,2,1)=temp6222(1,1,1)
       temp6221(1,2,2)=temp6222(2,1,1)
       temp6222(1,1,2)=temp6222(2,1,1)
       temp6222(1,2,1)=temp6222(2,1,1)
       temp6222(1,2,2)=temp6222(2,2,1)
       temp6222(2,1,2)=temp6222(2,2,1)
c                Step7
       tempC3000000=I14Z*(auxC3000000+tempC30000*F(3)-det3*temp00002(k,l
     &  ))
       temp00002(1,1)=I18Z*(aux00002(1,1)+4*F(4)*temp00001(1)+F(3)*temp0
     &  02(1,1)-det3*temp0041(1,k,l)+8*tempC3000000*ZZ(k,1,l,1))
       temp00002(2,1)=I18Z*(aux00002(2,1)+F(5)*temp00001(1)+2*F(4)*temp0
     &  0001(2)+F(3)*temp002(2,1)-det3*temp0042(1,k,l)+4*tempC3000000*(Z
     &  Z(k,1,l,2)+ZZ(k,2,l,1)))
       temp00002(2,2)=I18Z*(aux00002(2,2)+2*F(5)*temp00001(2)+F(3)*temp0
     &  02(2,2)-det3*temp0042(2,k,l)+8*tempC3000000*ZZ(k,2,l,2))
       temp00002(1,2)=temp00002(2,1)
       temp0041(1,1,1)=I22Z*(aux0041(1,1,1)+8*F(4)*temp003(1,1,1)+F(3)*t
     &  emp41(1,1,1)-det3*temp6111(1,k,l)+48*temp00002(1,1)*ZZ(k,1,l,1))
       temp0042(1,1,1)=I22Z*(aux0042(1,1,1)+F(5)*temp003(1,1,1)+6*F(4)*t
     &  emp003(2,1,1)+F(3)*temp42(1,1,1)-det3*temp6211(1,k,l)+24*temp000
     &  02(2,1)*ZZ(k,1,l,1)+12*temp00002(1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp0042(2,1,1)=I22Z*(aux0042(2,1,1)+2*F(5)*temp003(2,1,1)+4*F(4)
     &  *temp003(2,2,1)+F(3)*temp42(2,1,1)-det3*temp6221(1,k,l)+16*temp0
     &  0002(2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+8*(temp00002(2,2)*ZZ(k,1,l,1
     &  )+temp00002(1,1)*ZZ(k,2,l,2)))
       temp0042(2,2,1)=I22Z*(aux0042(2,2,1)+2*F(4)*temp003(2,2,2)+F(3)*t
     &  emp42(2,2,1)-det3*temp6222(1,k,l)+12*temp00002(2,2)*(ZZ(k,1,l,2)
     &  +ZZ(k,2,l,1))+24*temp00002(2,1)*ZZ(k,2,l,2)+temp003(2,2,1)*(6*r1
     &  0*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2)))
       temp0042(2,2,2)=I22Z*(aux0042(2,2,2)+4*F(5)*temp003(2,2,2)+F(3)*t
     &  emp42(2,2,2)-det3*temp6222(2,k,l)+48*temp00002(2,2)*ZZ(k,2,l,2))
       temp0041(1,1,2)=temp0042(1,1,1)
       temp0041(1,2,1)=temp0042(1,1,1)
       temp0041(1,2,2)=temp0042(2,1,1)
       temp0042(1,1,2)=temp0042(2,1,1)
       temp0042(1,2,1)=temp0042(2,1,1)
       temp0042(1,2,2)=temp0042(2,2,1)
       temp0042(2,1,2)=temp0042(2,2,1)
       temp511(1,1,1)=IX*(aux511(1,1,1)+det3*temp6111(1,1,jj)+10*temp004
     &  1(1,1,1)*Z(jj,1))
       temp521(1,1,1)=IX*(aux521(1,1,1)+det3*temp6211(1,1,jj)+8*temp0042
     &  (1,1,1)*Z(jj,1)+2*temp0041(1,1,1)*Z(jj,2))
       temp522(1,1,1)=IX*(aux522(1,1,1)+det3*temp6221(1,1,jj)+6*temp0042
     &  (2,1,1)*Z(jj,1)+4*temp0042(1,1,1)*Z(jj,2))
       temp522(2,1,1)=IX*(aux522(2,1,1)+det3*temp6222(1,1,jj)+4*temp0042
     &  (2,2,1)*Z(jj,1)+6*temp0042(2,1,1)*Z(jj,2))
       temp522(2,2,1)=IX*(aux522(2,2,1)+det3*temp6222(2,1,jj)+2*temp0042
     &  (2,2,2)*Z(jj,1)+8*temp0042(2,2,1)*Z(jj,2))
       temp522(2,2,2)=IX*(aux522(2,2,2)+det3*temp6222(2,2,jj)+10*temp004
     &  2(2,2,2)*Z(jj,2))
       temp511(1,1,2)=temp521(1,1,1)
       temp511(1,2,1)=temp521(1,1,1)
       temp511(1,2,2)=temp522(1,1,1)
       temp521(1,1,2)=temp522(1,1,1)
       temp521(1,2,1)=temp522(1,1,1)
       temp521(1,2,2)=temp522(2,1,1)
       temp522(1,1,2)=temp522(2,1,1)
       temp522(1,2,1)=temp522(2,1,1)
       temp522(1,2,2)=temp522(2,2,1)
       temp522(2,1,2)=temp522(2,2,1)
c                Step8
       temp00001(1)=I14Z*(aux00001(1)+2*tempC30000*F(4)+F(3)*temp001(1)-
     &  det3*temp003(1,k,l))
       temp00001(2)=I14Z*(aux00001(2)+tempC30000*F(5)+F(3)*temp001(2)-de
     &  t3*temp003(2,k,l))
       temp003(1,1,1)=I18Z*(aux003(1,1,1)+6*F(4)*temp002(1,1)+F(3)*temp3
     &  (1,1,1)-det3*temp511(1,k,l)+24*temp00001(1)*ZZ(k,1,l,1))
       temp003(2,1,1)=I18Z*(aux003(2,1,1)+F(5)*temp002(1,1)+4*F(4)*temp0
     &  02(2,1)+F(3)*temp3(2,1,1)-det3*temp521(1,k,l)+8*(temp00001(2)*ZZ
     &  (k,1,l,1)+temp00001(1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))))
       temp003(2,2,1)=I18Z*(aux003(2,2,1)+2*F(5)*temp002(2,1)+2*F(4)*tem
     &  p002(2,2)+F(3)*temp3(2,2,1)-det3*temp522(1,k,l)+8*(temp00001(2)*
     &  (ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp00001(1)*ZZ(k,2,l,2)))
       temp003(2,2,2)=I18Z*(aux003(2,2,2)+F(3)*temp3(2,2,2)-det3*temp522
     &  (2,k,l)+24*temp00001(2)*ZZ(k,2,l,2)+temp002(2,2)*(6*r10*(ZZ(k,1,
     &  l,2)+ZZ(k,2,l,1))+12*r21*ZZ(k,2,l,2)))
       temp003(1,1,2)=temp003(2,1,1)
       temp003(1,2,1)=temp003(2,1,1)
       temp003(1,2,2)=temp003(2,2,1)
       temp003(2,1,2)=temp003(2,2,1)
       temp41(1,1,1)=IX*(aux41(1,1,1)+det3*temp511(1,1,jj)+8*temp003(1,1
     &  ,1)*Z(jj,1))
       temp42(1,1,1)=IX*(aux42(1,1,1)+det3*temp521(1,1,jj)+6*temp003(2,1
     &  ,1)*Z(jj,1)+2*temp003(1,1,1)*Z(jj,2))
       temp42(2,1,1)=IX*(aux42(2,1,1)+det3*temp522(1,1,jj)+4*(temp003(2,
     &  2,1)*Z(jj,1)+temp003(2,1,1)*Z(jj,2)))
       temp42(2,2,1)=IX*(aux42(2,2,1)+det3*temp522(2,1,jj)+2*temp003(2,2
     &  ,2)*Z(jj,1)+6*temp003(2,2,1)*Z(jj,2))
       temp42(2,2,2)=IX*(aux42(2,2,2)+det3*temp522(2,2,jj)+8*temp003(2,2
     &  ,2)*Z(jj,2))
       temp41(1,1,2)=temp42(1,1,1)
       temp41(1,2,1)=temp42(1,1,1)
       temp41(1,2,2)=temp42(2,1,1)
       temp42(1,1,2)=temp42(2,1,1)
       temp42(1,2,1)=temp42(2,1,1)
       temp42(1,2,2)=temp42(2,2,1)
       temp42(2,1,2)=temp42(2,2,1)
c                Step9
       tempC30000=I10Z*(auxC30000+tempC300*F(3)-det3*temp002(k,l))
       temp002(1,1)=I14Z*(aux002(1,1)+4*F(4)*temp001(1)+F(3)*temp2(1,1)-
     &  det3*temp41(1,k,l)+8*tempC30000*ZZ(k,1,l,1))
       temp002(2,1)=I14Z*(aux002(2,1)+F(5)*temp001(1)+2*F(4)*temp001(2)+
     &  F(3)*temp2(2,1)-det3*temp42(1,k,l)+4*tempC30000*(ZZ(k,1,l,2)+ZZ(
     &  k,2,l,1)))
       temp002(2,2)=I14Z*(aux002(2,2)+2*F(5)*temp001(2)+F(3)*temp2(2,2)-
     &  det3*temp42(2,k,l)+8*tempC30000*ZZ(k,2,l,2))
       temp002(1,2)=temp002(2,1)
       temp3(1,1,1)=IX*(aux3(1,1,1)+det3*temp41(1,1,jj)+6*temp002(1,1)*Z
     &  (jj,1))
       temp3(2,1,1)=IX*(aux3(2,1,1)+det3*temp42(1,1,jj)+4*temp002(2,1)*Z
     &  (jj,1)+2*temp002(1,1)*Z(jj,2))
       temp3(2,2,1)=IX*(aux3(2,2,1)+det3*temp42(2,1,jj)+2*temp002(2,2)*Z
     &  (jj,1)+4*temp002(2,1)*Z(jj,2))
       temp3(2,2,2)=IX*(aux3(2,2,2)+det3*temp42(2,2,jj)+6*temp002(2,2)*Z
     &  (jj,2))
       temp3(1,1,2)=temp3(2,1,1)
       temp3(1,2,1)=temp3(2,1,1)
       temp3(1,2,2)=temp3(2,2,1)
       temp3(2,1,2)=temp3(2,2,1)
c                Step10
       temp001(1)=I10Z*(aux001(1)+2*tempC300*F(4)+F(3)*temp1(1)-det3*tem
     &  p3(1,k,l))
       temp001(2)=I10Z*(aux001(2)+tempC300*F(5)+F(3)*temp1(2)-det3*temp3
     &  (2,k,l))
       temp2(1,1)=IX*(aux2(1,1)+det3*temp3(1,1,jj)+4*temp001(1)*Z(jj,1))
       temp2(2,1)=IX*(aux2(2,1)+det3*temp3(2,1,jj)+2*(temp001(2)*Z(jj,1)
     &  +temp001(1)*Z(jj,2)))
       temp2(2,2)=IX*(aux2(2,2)+det3*temp3(2,2,jj)+4*temp001(2)*Z(jj,2))
       temp2(1,2)=temp2(2,1)
c                Step11
       tempC300=I6Z*(auxC300+tempC30*F(3)-det3*temp2(k,l))
       temp1(1)=IX*(aux1(1)+det3*temp2(1,jj)+2*tempC300*Z(jj,1))
       temp1(2)=IX*(aux1(2)+det3*temp2(2,jj)+2*tempC300*Z(jj,2))
c                Step12
       tempC30=IX*(auxC30+det3*temp1(jj))

         ac=5
         accuracyCR(1,0,ac) = abs(C30     /tempC30      -1d0) 
         accuracyCR(1,1,ac) = abs(Cij(1,1)/temp1(1)      -1d0) 
         accuracyCR(2,1,ac) = abs(Cij(2,1)/temp1(2)      -1d0)
         accuracyCR(1,2,ac) = abs(Cij(1,2)/temp2(1,1)    -1d0)
         accuracyCR(2,2,ac) = abs(Cij(2,2)/temp2(2,2)    -1d0)
         accuracyCR(3,2,ac) = abs(Cij(3,2)/temp2(2,1)    -1d0)
         accuracyCR(4,2,ac) = abs(Cij(4,2)/tempC300      -1d0)
         accuracyCR(1,3,ac) = abs(Cij(1,3)/temp3(1,1,1)  -1d0)
         accuracyCR(2,3,ac) = abs(Cij(2,3)/temp3(2,2,2)  -1d0)
         accuracyCR(3,3,ac) = abs(Cij(3,3)/temp3(2,1,1)  -1d0)
         accuracyCR(4,3,ac) = abs(Cij(4,3)/temp3(2,2,1)  -1d0)
         accuracyCR(5,3,ac) = abs(Cij(5,3)/temp001(1)    -1d0)
         accuracyCR(6,3,ac) = abs(Cij(6,3)/temp001(2)    -1d0)
         accuracyCR(1,4,ac) = abs(Cij(1,4)/temp41(1,1,1) -1d0)
         accuracyCR(2,4,ac) = abs(Cij(2,4)/temp42(2,2,2) -1d0)
         accuracyCR(3,4,ac) = abs(Cij(3,4)/temp42(1,1,1) -1d0)
         accuracyCR(4,4,ac) = abs(Cij(4,4)/temp42(2,1,1) -1d0)
         accuracyCR(5,4,ac) = abs(Cij(5,4)/temp42(2,2,1) -1d0)
         accuracyCR(6,4,ac) = abs(Cij(6,4)/temp002(1,1)  -1d0)
         accuracyCR(7,4,ac) = abs(Cij(7,4)/temp002(2,2)  -1d0)
         accuracyCR(8,4,ac) = abs(Cij(8,4)/temp002(2,1)  -1d0)
         accuracyCR(9,4,ac) = abs(Cij(9,4)/tempC30000    -1d0)


      DO I1=0,4
           accuracyC(i1,ac)=accuracyCR(1,i1,ac)
            DO I2=2,INDEX(I1)  
       if(accuracyCR(i2,i1,ac).gt.accuracyC(i1,ac)) then
          accuracyC(i1,ac)=accuracyCR(i2,i1,ac)
       endif
          enddo
        enddo

 500     C30=tempC30
           Cij(1,1)=temp1(1)
           Cij(2,1)=temp1(2)
           Cij(1,2)=temp2(1,1)
           Cij(2,2)=temp2(2,2)
           Cij(3,2)=temp2(2,1)
           Cij(4,2)=tempC300
           Cij(1,3)=temp3(1,1,1)
           Cij(2,3)=temp3(2,2,2)
           Cij(3,3)=temp3(2,1,1)
           Cij(4,3)=temp3(2,2,1)
           Cij(5,3)=temp001(1)
           Cij(6,3)=temp001(2)
           Cij(1,4)=temp41(1,1,1)
           Cij(2,4)=temp42(2,2,2)
           Cij(3,4)=temp42(1,1,1)
           Cij(4,4)=temp42(2,1,1)
           Cij(5,4)=temp42(2,2,1)
           Cij(6,4)=temp002(1,1)
           Cij(7,4)=temp002(2,2)
           Cij(8,4)=temp002(2,1)
           Cij(9,4)=tempC30000
           Cij(1,5)=temp511(1,1,1)
           Cij(2,5)=temp522(2,2,2)
           Cij(3,5)=temp521(1,1,1)
           Cij(4,5)=temp522(1,1,1)
           Cij(5,5)=temp522(2,1,1)
           Cij(6,5)=temp522(2,2,1)
           Cij(7,5)=temp003(1,1,1)
           Cij(8,5)=temp003(2,2,2)
           Cij(9,5)=temp003(2,1,1)
           Cij(10,5)=temp003(2,2,1)
           Cij(11,5)=temp00001(1)
           Cij(12,5)=temp00001(2)
           Cij(1,5)=temp511(1,1,1)
           Cij(2,5)=temp522(2,2,2)
           Cij(3,5)=temp521(1,1,1)
           Cij(4,5)=temp522(1,1,1)
           Cij(5,5)=temp522(2,1,1)
           Cij(6,5)=temp522(2,2,1)
           Cij(7,5)=temp003(1,1,1)
           Cij(8,5)=temp003(2,2,2)
           Cij(9,5)=temp003(2,1,1)
           Cij(10,5)=temp003(2,2,1)
           Cij(11,5)=temp00001(1)
           Cij(12,5)=temp00001(2)
           Cij(1,6)=temp6111(1,1,1)
           Cij(2,6)=temp6222(2,2,2)
           Cij(3,6)=temp6211(1,1,1)
           Cij(4,6)=temp6221(1,1,1)
           Cij(5,6)=temp6222(1,1,1)
           Cij(6,6)=temp6222(2,1,1)
           Cij(7,6)=temp6222(2,2,1)
           Cij(8,6)=temp0041(1,1,1)
           Cij(9,6)=temp0042(2,2,2)
           Cij(10,6)=temp0042(1,1,1)
           Cij(11,6)=temp0042(2,1,1)
           Cij(12,6)=temp0042(2,2,1)
           Cij(13,6)=temp00002(1,1)
           Cij(14,6)=temp00002(2,2)
           Cij(15,6)=temp00002(2,1)
           Cij(16,6)=tempC3000000
           Cij(1,7)=temp71111(1,1,1)
           Cij(2,7)=temp72222(2,2,2)
           Cij(3,7)=temp72111(1,1,1)
           Cij(4,7)=temp72211(1,1,1)
           Cij(5,7)=temp72221(1,1,1)
           Cij(6,7)=temp72222(1,1,1)
           Cij(7,7)=temp72222(2,1,1)
           Cij(8,7)=temp72222(2,2,1)
           Cij(9,7)=temp00511(1,1,1)
           Cij(10,7)=temp00522(2,2,2)
           Cij(11,7)=temp00521(1,1,1)
           Cij(12,7)=temp00522(1,1,1)
           Cij(13,7)=temp00522(2,1,1)
           Cij(14,7)=temp00522(2,2,1)
           Cij(15,7)=temp00003(1,1,1)
           Cij(16,7)=temp00003(2,2,2)
           Cij(17,7)=temp00003(2,1,1)
           Cij(18,7)=temp00003(2,2,1)
           Cij(19,7)=temp0000001(1)
           Cij(20,7)=temp0000001(2)
           Cij(1,8)=temp811111(1,1,1)
           Cij(2,8)=temp822222(2,2,2)
           Cij(3,8)=temp821111(1,1,1)
           Cij(4,8)=temp822111(1,1,1)
           Cij(5,8)=temp822211(1,1,1)
           Cij(6,8)=temp822221(1,1,1)
           Cij(7,8)=temp822222(1,1,1)
           Cij(8,8)=temp822222(2,1,1)
           Cij(9,8)=temp822222(2,2,1)
           Cij(10,8)=temp006111(1,1,1)
           Cij(11,8)=temp006222(2,2,2)
           Cij(12,8)=temp006211(1,1,1)
           Cij(13,8)=temp006221(1,1,1)
           Cij(14,8)=temp006222(1,1,1)
           Cij(15,8)=temp006222(2,1,1)
           Cij(16,8)=temp006222(2,2,1)
           Cij(17,8)=temp000041(1,1,1)
           Cij(18,8)=temp000042(2,2,2)
           Cij(19,8)=temp000042(1,1,1)
           Cij(20,8)=temp000042(2,1,1)
           Cij(21,8)=temp000042(2,2,1)
           Cij(22,8)=temp0000002(1,1)
           Cij(23,8)=temp0000002(2,2)
           Cij(24,8)=temp0000002(2,1)
           Cij(25,8)=tempC300000000
           Cij(11,9)=temp0071111(1,1,1)
           Cij(12,9)=temp0072222(2,2,2)
           Cij(13,9)=temp0072111(1,1,1)
           Cij(14,9)=temp0072211(1,1,1)
           Cij(15,9)=temp0072221(1,1,1)
           Cij(16,9)=temp0072222(1,1,1)
           Cij(17,9)=temp0072222(2,1,1)
           Cij(18,9)=temp0072222(2,2,1)
           Cij(19,9)=temp0000511(1,1,1)
           Cij(20,9)=temp0000522(2,2,2)
           Cij(21,9)=temp0000521(1,1,1)
           Cij(22,9)=temp0000522(1,1,1)
           Cij(23,9)=temp0000522(2,1,1)
           Cij(24,9)=temp0000522(2,2,1)
           Cij(25,9)=temp0000003(1,1,1)
           Cij(26,9)=temp0000003(2,2,2)
           Cij(27,9)=temp0000003(2,1,1)
           Cij(28,9)=temp0000003(2,2,1)
           Cij(29,9)=temp000000001(1)
           Cij(30,9)=temp000000001(2)


           Cij(1,9)=temp9111111(1,1,1)
           Cij(2,9)=temp9222222(2,2,2)
           Cij(3,9)=temp9211111(1,1,1)
           Cij(4,9)=temp9221111(1,1,1)
           Cij(5,9)=temp9222111(1,1,1)
           Cij(6,9)=temp9222211(1,1,1)
           Cij(7,9)=temp9222221(1,1,1)
           Cij(8,9)=temp9222222(1,1,1)
           Cij(9,9)=temp9222222(2,1,1)        
           Cij(10,9)=temp9222222(2,2,1)
c FC %
c FC %           Cij(1,10)=temp101111111(1,1,1)
c FC %           Cij(2,10)=temp102222222(2,2,2)
c FC %           Cij(3,10)=temp102111111(1,1,1)
c FC %           Cij(4,10)=temp102211111(1,1,1)
c FC %           Cij(5,10)=temp102221111(1,1,1)
c FC %           Cij(6,10)=temp102222111(1,1,1)
c FC %           Cij(7,10)=temp102222211(1,1,1)
c FC %           Cij(8,10)=temp102222221(1,1,1)
c FC %           Cij(9,10)=temp102222222(1,1,1)        
c FC %          Cij(10,10)=temp102222222(2,1,1)
c FC %          Cij(11,10)=temp102222222(2,2,1)
c FC %
c FC %           Cij(12,10)=temp00811111(1,1,1)
c FC %           Cij(13,10)=temp00822222(2,2,2)
c FC %           Cij(14,10)=temp00821111(1,1,1)
c FC %           Cij(15,10)=temp00822111(1,1,1)
c FC %           Cij(16,10)=temp00822211(1,1,1)
c FC %           Cij(17,10)=temp00822221(1,1,1)
c FC %           Cij(18,10)=temp00822222(1,1,1)
c FC %           Cij(19,10)=temp00822222(2,1,1)
c FC %           Cij(20,10)=temp00822222(2,2,1)
c FC %           Cij(21,10)=temp00006111(1,1,1)
c FC %           Cij(22,10)=temp00006222(2,2,2)
c FC %           Cij(23,10)=temp00006211(1,1,1)
c FC %           Cij(24,10)=temp00006221(1,1,1)
c FC %           Cij(25,10)=temp00006222(1,1,1)
c FC %           Cij(26,10)=temp00006222(2,1,1)
c FC %           Cij(27,10)=temp00006222(2,2,1)
c FC %           Cij(28,10)=temp00000041(1,1,1)
c FC %           Cij(29,10)=temp00000042(2,2,2)
c FC %           Cij(30,10)=temp00000042(1,1,1)
c FC %           Cij(31,10)=temp00000042(2,1,1)
c FC %           Cij(32,10)=temp00000042(2,2,1)
c FC %           Cij(33,10)=temp000000002(1,1)
c FC %           Cij(34,10)=temp000000002(2,2)
c FC %           Cij(35,10)=temp000000002(2,1)
c FC %           Cij(36,10)=tempC30000000000
c FC %
c FC %
c FC %           Cij(1,11)=temp1111111111(1,1,1)
c FC %           Cij(2,11)=temp1122222222(2,2,2)
c FC %           Cij(3,11)=temp1121111111(1,1,1)
c FC %           Cij(4,11)=temp1122111111(1,1,1)
c FC %           Cij(5,11)=temp1122211111(1,1,1)
c FC %           Cij(6,11)=temp1122221111(1,1,1)
c FC %           Cij(7,11)=temp1122222111(1,1,1)
c FC %           Cij(8,11)=temp1122222211(1,1,1)
c FC %           Cij(9,11)=temp1122222221(1,1,1)        
c FC %          Cij(10,11)=temp1122222222(1,1,1)
c FC %          Cij(11,11)=temp1122222222(2,1,1)
c FC %          Cij(12,11)=temp1122222222(2,2,1)
c FC %
c FC %           Cij(13,11)=temp009111111(1,1,1)
c FC %           Cij(14,11)=temp009222222(2,2,2)
c FC %           Cij(15,11)=temp009211111(1,1,1)
c FC %           Cij(16,11)=temp009221111(1,1,1)
c FC %           Cij(17,11)=temp009222111(1,1,1)
c FC %           Cij(18,11)=temp009222211(1,1,1)
c FC %           Cij(19,11)=temp009222221(1,1,1)
c FC %           Cij(20,11)=temp009222222(1,1,1)
c FC %           Cij(21,11)=temp009222222(2,1,1)        
c FC %           Cij(22,11)=temp009222222(2,2,1)
c FC %
c FC %           Cij(23,11)=temp000071111(1,1,1)
c FC %           Cij(24,11)=temp000072222(2,2,2)
c FC %           Cij(25,11)=temp000072111(1,1,1)
c FC %           Cij(26,11)=temp000072211(1,1,1)
c FC %           Cij(27,11)=temp000072221(1,1,1)
c FC %           Cij(28,11)=temp000072222(1,1,1)
c FC %           Cij(29,11)=temp000072222(2,1,1)
c FC %           Cij(30,11)=temp000072222(2,2,1)
c FC %           Cij(31,11)=temp000000511(1,1,1)
c FC %           Cij(32,11)=temp000000522(2,2,2)
c FC %           Cij(33,11)=temp000000521(1,1,1)
c FC %           Cij(34,11)=temp000000522(1,1,1)
c FC %           Cij(35,11)=temp000000522(2,1,1)
c FC %           Cij(36,11)=temp000000522(2,2,1)
c FC %           Cij(37,11)=temp000000003(1,1,1)
c FC %           Cij(38,11)=temp000000003(2,2,2)
c FC %           Cij(39,11)=temp000000003(2,1,1)
c FC %           Cij(40,11)=temp000000003(2,2,1)
c FC %           Cij(41,11)=temp00000000001(1)
c FC %           Cij(42,11)=temp00000000001(2)
c FC %
c FC %
c FC %c FC %           Cij(1,12)=temp12111111111(1,1,1)
c FC %c FC %           Cij(2,12)=temp12222222222(2,2,2)
c FC %c FC %           Cij(3,12)=temp12211111111(1,1,1)
c FC %c FC %           Cij(4,12)=temp12221111111(1,1,1)
c FC %c FC %           Cij(5,12)=temp12222111111(1,1,1)
c FC %c FC %           Cij(6,12)=temp12222211111(1,1,1)
c FC %c FC %           Cij(7,12)=temp12222221111(1,1,1)
c FC %c FC %           Cij(8,12)=temp12222222111(1,1,1)
c FC %c FC %           Cij(9,12)=temp12222222211(1,1,1)        
c FC %c FC %          Cij(10,12)=temp12222222221(1,1,1)
c FC %c FC %          Cij(11,12)=temp12222222222(1,1,1)
c FC %c FC %          Cij(12,12)=temp12222222222(2,1,1)
c FC %c FC %          Cij(13,12)=temp12222222222(2,1,1)
c FC %
c FC %           Cij(14,12)=temp00101111111(1,1,1)
c FC %           Cij(15,12)=temp00102222222(2,2,2)
c FC %           Cij(16,12)=temp00102111111(1,1,1)
c FC %           Cij(17,12)=temp00102211111(1,1,1)
c FC %           Cij(18,12)=temp00102221111(1,1,1)
c FC %           Cij(19,12)=temp00102222111(1,1,1)
c FC %           Cij(20,12)=temp00102222211(1,1,1)
c FC %           Cij(21,12)=temp00102222221(1,1,1)
c FC %           Cij(22,12)=temp00102222222(1,1,1)        
c FC %           Cij(23,12)=temp00102222222(2,1,1)
c FC %           Cij(24,12)=temp00102222222(2,2,1)
c FC %
c FC %           Cij(25,12)=temp0000811111(1,1,1)
c FC %           Cij(26,12)=temp0000822222(2,2,2)
c FC %           Cij(27,12)=temp0000821111(1,1,1)
c FC %           Cij(28,12)=temp0000822111(1,1,1)
c FC %           Cij(29,12)=temp0000822211(1,1,1)
c FC %           Cij(30,12)=temp0000822221(1,1,1)
c FC %           Cij(31,12)=temp0000822222(1,1,1)
c FC %           Cij(32,12)=temp0000822222(2,1,1)
c FC %           Cij(33,12)=temp0000822222(2,2,1)
c FC %           Cij(34,12)=temp0000006111(1,1,1)
c FC %           Cij(35,12)=temp0000006222(2,2,2)
c FC %           Cij(36,12)=temp0000006211(1,1,1)
c FC %           Cij(37,12)=temp0000006221(1,1,1)
c FC %           Cij(38,12)=temp0000006222(1,1,1)
c FC %           Cij(39,12)=temp0000006222(2,1,1)
c FC %           Cij(40,12)=temp0000006222(2,2,1)
c FC %           Cij(41,12)=temp0000000041(1,1,1)
c FC %           Cij(42,12)=temp0000000042(2,2,2)
c FC %           Cij(43,12)=temp0000000042(1,1,1)
c FC %           Cij(44,12)=temp0000000042(2,1,1)
c FC %           Cij(45,12)=temp0000000042(2,2,1)
c FC %           Cij(46,12)=temp00000000002(1,1)
c FC %           Cij(47,12)=temp00000000002(2,2)
c FC %           Cij(48,12)=temp00000000002(2,1)
c FC %           Cij(49,12)=tempC3000000000000
c FC %

      return
      End
