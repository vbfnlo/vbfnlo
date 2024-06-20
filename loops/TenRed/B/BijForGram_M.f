       subroutine ten_red2_forGram_M(m0,p1sq,A2,A1,B012,B12)
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c************************************************************************************
       Implicit none
       Real*8 p1sq,m0,m0sq
       Complex*16 B012,A2,A1,A202,A402,A602,A802,A1002
       Complex*16 B12(6,11)
        Real*8 Inv2,Inv5,Inv6,Inv8,Inv12,Inv18,Inv20,Inv24,Inv36,Inv48,Inv60
        Real*8 Inv72,Inv192,Inv300,Inv360,Inv600,Inv1344,Inv1440,Inv1470,Inv1920
        Real*8 Inv3600,Inv4536,Inv5880,Inv7920,Inv10890,Inv23040,Inv23520,Inv26880
        Real*8 Inv45360,Inv53760,Inv100800,Inv190080,Inv201600,Inv544320,Inv609840
        Real*8 Inv1088640,Inv2419200,Inv2439360,Inv7983360,Inv14636160,Inv58544640
        Real*8 Inv63866880,Inv127733760
        Real*8 p1sq2,p1sq3,p1sq4,p1sq5,m0sq2,m0sq3,m0sq4,m0sq5
        Real*8 Invp1sq
        real*8 inv11,inv14,inv22,inv35,inv63,inv80,inv99
        real*8 inv294,inv308,inv648,inv900,inv1210,inv1584,inv21780

c FC %       index(1)=1
c FC %       index(2)=2
c FC %       index(3)=2
c FC %       index(4)=3
c FC %       index(5)=3
c FC %       index(6)=4
c FC %       index(7)=4
c FC %       index(8)=5
c FC %       index(9)=5
c FC %       index(10)=6
c FC %       index(11)=6
c       index(12)=49
       If(abs(p1sq).gt.1d-8) then
        Inv2=1d0/2d0
        Inv5=1d0/5d0
        Inv6=1d0/6d0
        Inv8=1d0/8d0
        Inv11=1d0/11d0
        Inv12=1d0/12d0
        inv14=1d0/14d0
        Inv18=1d0/18d0
        Inv20=1d0/20d0
        Inv22=1d0/22d0
        Inv24=1d0/24d0
        INV35=1D0/35D0
        Inv36=1d0/36d0
        Inv48=1d0/48d0
        Inv60=1d0/60d0
        Inv63=1D0/63D0
        Inv72=1d0/72d0
        iNV80=1D0/80D0
        Inv99=1D0/99D0
        Inv192=1d0/192d0
        Inv294=1D0/294D0
        Inv300=1d0/300d0
        Inv308=1d0/308d0
        Inv360=1d0/360d0
        Inv600=1d0/600d0
        Inv648=1D0/648D0
        INV900=1D0/900D0
        Inv1210=1D0/1210D0
        Inv1344=1d0/1344d0
        Inv1440=1d0/1440d0
        Inv1470=1d0/1470d0
        Inv1584=1D0/1584D0
        Inv1920=1d0/1920d0
        Inv3600=1d0/3600d0
        Inv4536=1d0/4536d0
        Inv5880=1d0/5880d0
        Inv7920=1d0/7920d0
        Inv10890=1d0/10890d0
        Inv21780=1D0/21780D0
        Inv23040=1d0/23040d0
        Inv23520=1d0/23520d0
        Inv26880=1d0/26880d0
        Inv45360=1d0/45360d0
        Inv53760=1d0/53760d0
        Inv100800=1d0/100800d0
        Inv190080=1d0/190080d0
        Inv201600=1d0/201600d0
        Inv544320=1d0/544320d0
        Inv609840=1d0/609840d0
        Inv1088640=1d0/1088640d0
        Inv2419200=1d0/2419200d0
        Inv2439360=1d0/2439360d0
        Inv7983360=1d0/7983360d0
        Inv14636160=1d0/14636160d0
        Inv58544640=1d0/58544640d0
        Inv63866880=1d0/63866880d0
        Inv127733760=1d0/127733760d0
        p1sq2=p1sq*p1sq
        p1sq3=p1sq2*p1sq
        p1sq4=p1sq3*p1sq
        p1sq5=p1sq4*p1sq
        m0sq=m0*m0
        m0sq2=m0sq*m0sq
        m0sq3=m0sq2*m0sq
        m0sq4=m0sq3*m0sq
        m0sq5=m0sq4*m0sq
        Invp1sq=1d0/p1sq        
       A202=Inv8*m0sq*(2*A1+m0sq)
       A402=Inv24*(A1+5*Inv6*m0sq)*m0sq2
       A602=Inv192*(A1+13*Inv12*m0sq)*m0sq3
       A802=Inv1920*(A1+77*Inv60*m0sq)*m0sq4
       A1002=Inv23040*(A1+29*Inv20*m0sq)*m0sq5
       B12(1,1)=-(B012*Inv2)
       B12(2,2)=Inv18*(3*A1+6*(1+B012)*m0sq+p1sq*(-1+3*B12(1,1)))
       B12(1,2)=Inv2*Invp1sq*(A1-p1sq*B12(1,1)-2*B12(2,2))
       B12(2,3)=Inv48*(-6*A1+p1sq+6*m0sq*(-1+2*B12(1,1))+6*p1sq*B12(1,2)
     &  )
       B12(1,3)=-(Inv2*Invp1sq*(A1+p1sq*B12(1,2)+4*B12(2,3)))
       B12(3,4)=Inv600*(60*A202+30*m0sq2-10*m0sq*(p1sq-12*B12(2,2))+p1sq
     &  *(p1sq+60*B12(2,3)))
       B12(2,4)=Inv300*(30*A1+20*m0sq-3*p1sq+60*m0sq*B12(1,2)+30*p1sq*B1
     &  2(1,3))
       B12(1,4)=Inv2*Invp1sq*(A1-p1sq*B12(1,3)-6*B12(2,4))
       B12(3,5)=Inv1440*(-120*A202-30*m0sq2+10*m0sq*(p1sq+24*B12(2,3))-p
     &  1sq*(p1sq-120*B12(2,4)))
       B12(2,5)=Inv360*(-30*A1-15*m0sq+2*p1sq+60*m0sq*B12(1,3)+30*p1sq*B
     &  12(1,4))
       B12(1,5)=-(Inv2*Invp1sq*(A1+p1sq*B12(1,4)+8*B12(2,5)))
       B12(4,6)=Inv23520*(1680*A402+140*m0sq3-70*m0sq2*p1sq-p1sq3+14*m0s
     &  q*(p1sq2+240*B12(3,4))+1680*p1sq*B12(3,5))
       B12(3,6)=Inv5880*(70*m0sq2-21*m0sq*p1sq+2*p1sq2+420*(A202+2*m0sq*
     &  B12(2,4)+p1sq*B12(2,5)))
       B12(2,6)=Inv35*m0sq-Inv294*p1sq+Inv14*(A1+2*m0sq*B12(1,4)+p1sq*B1
     &  2(1,5))
       B12(1,6)=Inv2*Invp1sq*(A1-p1sq*B12(1,5)-10*B12(2,6))
       B12(4,7)=Inv53760*(-3360*A402-140*m0sq3+70*m0sq2*p1sq+p1sq3-14*m0
     &  sq*(p1sq2-480*B12(3,5))+3360*p1sq*B12(3,6))
       B12(3,7)=Inv26880*(-210*m0sq2+56*m0sq*p1sq-5*p1sq2+1680*(-A202+2*
     &  m0sq*B12(2,5)+p1sq*B12(2,6)))
       B12(2,7)=Inv1344*(-28*m0sq+3*p1sq+84*(-A1+2*m0sq*B12(1,5)+p1sq*B1
     &  2(1,6)))
       B12(1,7)=-(Inv2*Invp1sq*(A1+p1sq*B12(1,6)+12*B12(2,7)))
       B12(5,8)=Inv1088640*(60480*A602+630*m0sq4-420*m0sq3*p1sq+126*m0sq
     &  2*p1sq2+p1sq4-18*m0sq*(p1sq3-6720*B12(4,6))+60480*p1sq*B12(4,7))
       B12(4,8)=Inv544320*(840*m0sq3-378*m0sq2*p1sq+72*m0sq*p1sq2-5*p1sq
     &  3+30240*(A402+2*m0sq*B12(3,6)+p1sq*B12(3,7)))
       B12(3,8)=Inv45360*(252*m0sq2-60*m0sq*p1sq+5*p1sq2+2520*(A202+2*m0
     &  sq*B12(2,6)+p1sq*B12(2,7)))
       B12(2,8)=Inv63*m0sq-Inv648*p1sq+Inv18*(A1+2*m0sq*B12(1,6)+p1sq*B1
     &  2(1,7))
       B12(1,8)=Inv2*Invp1sq*(A1-p1sq*B12(1,7)-14*B12(2,8))
       B12(5,9)=Inv2419200*(-120960*A602-630*m0sq4+420*m0sq3*p1sq-126*m0
     &  sq2*p1sq2-p1sq4+18*m0sq*(p1sq3+13440*B12(4,7))+120960*p1sq*B12(4
     &  ,8))
       B12(4,9)=Inv201600*(-10080*A402-210*m0sq3+84*m0sq2*p1sq+p1sq3-15*
     &  m0sq*(p1sq2-1344*B12(3,7))+10080*p1sq*B12(3,8))
       B12(3,9)=Inv100800*(-420*m0sq2+90*m0sq*p1sq-7*p1sq2+5040*(-A202+2
     &  *m0sq*B12(2,7)+p1sq*B12(2,8)))
       B12(2,9)=-(Inv80*m0sq)+Inv900*p1sq+Inv20*(-A1+2*m0sq*B12(1,7)+p1s
     &  q*B12(1,8))
       B12(1,9)=-(Inv2*Invp1sq*(A1+p1sq*B12(1,8)+16*B12(2,9)))
       B12(6,10)=A802*Inv22+Inv58544640*(2772*m0sq5-2310*m0sq4*p1sq+924*
     &  m0sq3*p1sq2-198*m0sq2*p1sq3+22*m0sq*p1sq4-p1sq5)+Inv11*m0sq*B12(
     &  5,8)+Inv22*p1sq*B12(5,9)
       B12(5,10)=Inv14636160*(2310*m0sq4-1386*m0sq3*p1sq+396*m0sq2*p1sq2
     &  -55*m0sq*p1sq3+3*p1sq4+665280*(A602+2*m0sq*B12(4,8)+p1sq*B12(4,9
     &  )))
       B12(4,10)=Inv2439360*(1848*m0sq3-660*m0sq2*p1sq+110*m0sq*p1sq2-7*
     &  p1sq3+110880*(A402+2*m0sq*B12(3,8)+p1sq*B12(3,9)))
       B12(3,10)=Inv308*m0sq2-Inv1584*m0sq*p1sq+Inv21780*p1sq2+Inv22*(A2
     &  02+2*m0sq*B12(2,8)+p1sq*B12(2,9))
       B12(2,10)=Inv99*m0sq-Inv1210*p1sq+Inv22*(A1+2*m0sq*B12(1,8)+p1sq*
     &  B12(1,9))
       B12(1,10)=Inv2*Invp1sq*(A1-p1sq*B12(1,9)-18*B12(2,10))
       B12(6,11)=-(A802*Inv24)+Inv127733760*(-2772*m0sq5+2310*m0sq4*p1sq
     &  -924*m0sq3*p1sq2+198*m0sq2*p1sq3-22*m0sq*p1sq4+p1sq5)+Inv12*m0sq
     &  *B12(5,9)+Inv24*p1sq*B12(5,10)
       B12(5,11)=Inv63866880*(-6930*m0sq4+3696*m0sq3*p1sq-990*m0sq2*p1sq
     &  2+132*m0sq*p1sq3-7*p1sq4+2661120*(-A602+2*m0sq*B12(4,9)+p1sq*B12
     &  (4,10)))
       B12(4,11)=Inv7983360*(-4620*m0sq3+1485*m0sq2*p1sq-231*m0sq*p1sq2+
     &  14*p1sq3+332640*(-A402+2*m0sq*B12(3,9)+p1sq*B12(3,10)))
       B12(3,11)=Inv190080*(-495*m0sq2+88*m0sq*p1sq-6*p1sq2+7920*(-A202+
     &  2*m0sq*B12(2,9)+p1sq*B12(2,10)))
       B12(2,11)=Inv7920*(-66*m0sq+5*p1sq+330*(-A1+2*m0sq*B12(1,9)+p1sq*
     &  B12(1,10)))
       B12(1,11)=-(Inv2*Invp1sq*(A1+p1sq*B12(1,10)+20*B12(2,11)))
c FC %       DO I1=1,11
c FC %          print*,''
c FC %          print*, 'Rank=', i1
c FC %              DO I2=1,INDEX(I1)
c FC %                 print*,'B12(I2,I1)',B12(I2,I1)
c FC %              enddo
c FC %       enddo       
       else
       print*,''
       print*,'no yet implemented'
       print*,''
       B12(1,1)=0d0
       B12(2,2)=0d0
       B12(1,2)=0d0
       B12(2,3)=0d0
       B12(1,3)=0d0
       B12(3,4)=0d0
       B12(2,4)=0d0
       B12(1,4)=0d0
       B12(3,5)=0d0
       B12(2,5)=0d0
       B12(1,5)=0d0
       B12(4,6)=0d0
       B12(3,6)=0d0
       B12(2,6)=0d0
       B12(1,6)=0d0
       B12(4,7)=0d0
       B12(3,7)=0d0
       B12(2,7)=0d0
       B12(1,7)=0d0
       B12(5,8)=0d0
       B12(4,8)=0d0
       B12(3,8)=0d0
       B12(2,8)=0d0
       B12(1,8)=0d0
       B12(5,9)=0d0
       B12(4,9)=0d0
       B12(3,9)=0d0
       B12(2,9)=0d0
       B12(1,9)=0d0
       B12(6,10)=0d0
       B12(5,10)=0d0
       B12(4,10)=0d0
       B12(3,10)=0d0
       B12(2,10)=0d0
       B12(1,10)=0d0
       B12(6,11)=0d0
       B12(5,11)=0d0
       B12(4,11)=0d0
       B12(3,11)=0d0
       B12(2,11)=0d0
       B12(1,11)=0d0
       Endif
       return
       End
