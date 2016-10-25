       subroutine ten_red2_forGram_G(m0,m1,p1sq,A2,A1,B012,B12)
       Implicit none
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c************************************************************************************
       Real*8 p1sq,m0,m1,m0sq,m1sq
       Complex*16 B012,A2,A1,A202,A402,A602,A802,A1002,A1202
       Complex*16 B12(6,11)
        Real*8 Inv2,Inv5,Inv6,Inv8,Inv12,Inv18,Inv20,Inv24,Inv36,Inv48,Inv60
        Real*8 Inv72,Inv192,Inv300,Inv360,Inv600,Inv1344,Inv1440,Inv1470,Inv1920
        Real*8 Inv3600,Inv4536,Inv5880,Inv7920,Inv10890,Inv23040,Inv23520,Inv26880
        Real*8 Inv45360,Inv53760,Inv100800,Inv190080,Inv201600,Inv544320,Inv609840
        Real*8 Inv1088640,Inv2419200,Inv2439360,Inv7983360,Inv14636160,Inv58544640
        Real*8 Inv63866880,Inv127733760
        Real*8 p1sq2,p1sq3,p1sq4,p1sq5,m0sq2,m0sq3,m0sq4,m0sq5,m1sq2,m1sq3,m1sq4,m1sq5
        Real*8 Invp1sq
       If(abs(p1sq).gt.1d-8) then
        Inv2=1d0/2d0
        Inv5=1d0/5d0
        Inv6=1d0/6d0
        Inv8=1d0/8d0
        Inv12=1d0/12d0
        Inv18=1d0/18d0
        Inv20=1d0/20d0
        Inv24=1d0/24d0
        Inv36=1d0/36d0
        Inv48=1d0/48d0
        Inv60=1d0/60d0
        Inv72=1d0/72d0
        Inv192=1d0/192d0
        Inv300=1d0/300d0
        Inv360=1d0/360d0
        Inv600=1d0/600d0
        Inv1344=1d0/1344d0
        Inv1440=1d0/1440d0
        Inv1470=1d0/1470d0
        Inv1920=1d0/1920d0
        Inv3600=1d0/3600d0
        Inv4536=1d0/4536d0
        Inv5880=1d0/5880d0
        Inv7920=1d0/7920d0
        Inv10890=1d0/10890d0
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
        m1sq=m1*m1
        m1sq2=m1sq*m1sq
        m1sq3=m1sq2*m1sq
        m1sq4=m1sq3*m1sq
        m1sq5=m1sq4*m1sq
        m0sq=m0*m0
        m0sq2=m0sq*m0sq
        m0sq3=m0sq2*m0sq
        m0sq4=m0sq3*m0sq
        m0sq5=m0sq4*m0sq
        Invp1sq=1d0/p1sq        
       A202=Inv8*m1sq*(2*A2+m1sq)
       A402=Inv24*(A2+5*Inv6*m1sq)*m1sq2
       A602=Inv192*(A2+13*Inv12*m1sq)*m1sq3
       A802=Inv1920*(A2+77*Inv60*m1sq)*m1sq4
       A1002=Inv23040*(A2+29*Inv20*m1sq)*m1sq5
       B12(1,1)=-(Inv2*Invp1sq*(-A1+A2+B012*(m0sq-m1sq+p1sq)))
       B12(2,2)=Inv18*(3*A2+6*B012*m0sq+3*(m0sq+m1sq)-p1sq+3*(m0sq-m1sq+
     &  p1sq)*B12(1,1))
       B12(1,2)=Inv2*Invp1sq*(A2-(m0sq-m1sq+p1sq)*B12(1,1)-2*B12(2,2))
       B12(2,3)=Inv48*(-6*A2-2*m0sq-4*m1sq+p1sq+12*m0sq*B12(1,1)+6*(m0sq
     &  -m1sq+p1sq)*B12(1,2))
       B12(1,3)=-(Inv2*Invp1sq*(A2+(m0sq-m1sq+p1sq)*B12(1,2)+4*B12(2,3))
     &  )
       B12(3,4)=Inv600*(10*(m0sq2+m0sq*m1sq+m1sq2)-5*(m0sq+m1sq)*p1sq+p1
     &  sq2+60*(A202+2*m0sq*B12(2,2)+(m0sq-m1sq+p1sq)*B12(2,3)))
       B12(2,4)=Inv300*(5*m0sq+15*m1sq-3*p1sq+30*(A2+2*m0sq*B12(1,2)+(m0
     &  sq-m1sq+p1sq)*B12(1,3)))
       B12(1,4)=Inv2*Invp1sq*(A2-(m0sq-m1sq+p1sq)*B12(1,3)-6*B12(2,4))
       B12(3,5)=Inv1440*(-5*(m0sq2+2*m0sq*m1sq+3*m1sq2)+4*m0sq*p1sq+6*m1
     &  sq*p1sq-p1sq2+120*(-A202+2*m0sq*B12(2,3)+(m0sq-m1sq+p1sq)*B12(2,
     &  4)))
       B12(2,5)=Inv360*(-30*A2-3*(m0sq+4*m1sq)+2*p1sq+60*m0sq*B12(1,3)+3
     &  0*(m0sq-m1sq+p1sq)*B12(1,4))
       B12(1,5)=-(Inv2*Invp1sq*(A2+(m0sq-m1sq+p1sq)*B12(1,4)+8*B12(2,5))
     &  )
       B12(4,6)=Inv23520*(1680*A402+35*(m0sq+m1sq)*(m0sq2+m1sq2)-7*(3*m0
     &  sq2+4*m0sq*m1sq+3*m1sq2)*p1sq+7*(m0sq+m1sq)*p1sq2-p1sq3+1680*(2*
     &  m0sq*B12(3,4)+(m0sq-m1sq+p1sq)*B12(3,5)))
       B12(3,6)=Inv5880*(7*(m0sq2+3*m0sq*m1sq+6*m1sq2)-7*(m0sq+2*m1sq)*p
     &  1sq+2*p1sq2+420*(A202+2*m0sq*B12(2,4)+(m0sq-m1sq+p1sq)*B12(2,5))
     &  )
       B12(2,6)=Inv1470*(7*(m0sq+5*m1sq)-5*p1sq+105*(A2+2*m0sq*B12(1,4)+
     &  (m0sq-m1sq+p1sq)*B12(1,5)))
       B12(1,6)=Inv2*Invp1sq*(A2-(m0sq-m1sq+p1sq)*B12(1,5)-10*B12(2,6))
       B12(4,7)=Inv53760*(-14*(m0sq3+2*m0sq2*m1sq+3*m0sq*m1sq2+4*m1sq3)+
     &  14*(m0sq2+2*m0sq*m1sq+2*m1sq2)*p1sq-2*(3*m0sq+4*m1sq)*p1sq2+p1sq
     &  3+3360*(-A402+2*m0sq*B12(3,5)+(m0sq-m1sq+p1sq)*B12(3,6)))
       B12(3,7)=Inv26880*(-14*(m0sq2+4*m0sq*m1sq+10*m1sq2)+8*(2*m0sq+5*m
     &  1sq)*p1sq-5*p1sq2+1680*(-A202+2*m0sq*B12(2,5)+(m0sq-m1sq+p1sq)*B
     &  12(2,6)))
       B12(2,7)=Inv1344*(-4*(m0sq+6*m1sq)+3*p1sq+84*(-A2+2*m0sq*B12(1,5)
     &  +(m0sq-m1sq+p1sq)*B12(1,6)))
       B12(1,7)=-(Inv2*Invp1sq*(A2+(m0sq-m1sq+p1sq)*B12(1,6)+12*B12(2,7)
     &  ))
       B12(5,8)=Inv1088640*(126*(m0sq4+m0sq3*m1sq+m0sq2*m1sq2+m0sq*m1sq3
     &  +m1sq4)-42*(m0sq+m1sq)*(2*m0sq2+m0sq*m1sq+2*m1sq2)*p1sq+18*(2*m0
     &  sq2+3*m0sq*m1sq+2*m1sq2)*p1sq2-9*(m0sq+m1sq)*p1sq3+p1sq4+60480*(
     &  A602+2*m0sq*B12(4,6)+(m0sq-m1sq+p1sq)*B12(4,7)))
       B12(4,8)=Inv544320*(42*(m0sq3+3*m0sq2*m1sq+6*m0sq*m1sq2+10*m1sq3)
     &  -18*(3*m0sq2+8*m0sq*m1sq+10*m1sq2)*p1sq+9*(3*m0sq+5*m1sq)*p1sq2-
     &  5*p1sq3+30240*(A402+2*m0sq*B12(3,6)+(m0sq-m1sq+p1sq)*B12(3,7)))
       B12(3,8)=Inv45360*(12*(m0sq2+5*m0sq*m1sq+15*m1sq2)-15*(m0sq+3*m1s
     &  q)*p1sq+5*p1sq2+2520*(A202+2*m0sq*B12(2,6)+(m0sq-m1sq+p1sq)*B12(
     &  2,7)))
       B12(2,8)=Inv4536*(9*(m0sq+7*m1sq)-7*p1sq+252*(A2+2*m0sq*B12(1,6)+
     &  (m0sq-m1sq+p1sq)*B12(1,7)))
       B12(1,8)=Inv2*Invp1sq*(A2-(m0sq-m1sq+p1sq)*B12(1,7)-14*B12(2,8))
       B12(5,9)=Inv2419200*(-42*(m0sq4+2*m0sq3*m1sq+3*m0sq2*m1sq2+4*m0sq
     &  *m1sq3+5*m1sq4)+12*(4*m0sq3+9*m0sq2*m1sq+12*m0sq*m1sq2+10*m1sq3)
     &  *p1sq-9*(3*m0sq2+6*m0sq*m1sq+5*m1sq2)*p1sq2+2*(4*m0sq+5*m1sq)*p1
     &  sq3-p1sq4+120960*(-A602+2*m0sq*B12(4,7)+(m0sq-m1sq+p1sq)*B12(4,8
     &  )))
       B12(4,9)=Inv201600*(-6*(m0sq3+4*m0sq2*m1sq+10*m0sq*m1sq2+20*m1sq3
     &  )+3*(3*m0sq2+10*m0sq*m1sq+15*m1sq2)*p1sq-5*(m0sq+2*m1sq)*p1sq2+p
     &  1sq3+10080*(-A402+2*m0sq*B12(3,7)+(m0sq-m1sq+p1sq)*B12(3,8)))
       B12(3,9)=Inv100800*(-15*(m0sq2+6*m0sq*m1sq+21*m1sq2)+10*(2*m0sq+7
     &  *m1sq)*p1sq-7*p1sq2+5040*(-A202+2*m0sq*B12(2,7)+(m0sq-m1sq+p1sq)
     &  *B12(2,8)))
       B12(2,9)=Inv3600*(-5*(m0sq+8*m1sq)+4*p1sq+180*(-A2+2*m0sq*B12(1,7
     &  )+(m0sq-m1sq+p1sq)*B12(1,8)))
       B12(1,9)=-(Inv2*Invp1sq*(A2+(m0sq-m1sq+p1sq)*B12(1,8)+16*B12(2,9)
     &  ))
       B12(6,10)=Inv58544640*(2661120*A802+462*(m0sq+m1sq)*(m0sq2-m0sq*m
     &  1sq+m1sq2)*(m0sq2+m0sq*m1sq+m1sq2)-66*(5*m0sq4+8*m0sq3*m1sq+9*m0
     &  sq2*m1sq2+8*m0sq*m1sq3+5*m1sq4)*p1sq+33*(m0sq+m1sq)*(5*m0sq2+4*m
     &  0sq*m1sq+5*m1sq2)*p1sq2-11*(5*m0sq2+8*m0sq*m1sq+5*m1sq2)*p1sq3+1
     &  1*(m0sq+m1sq)*p1sq4-p1sq5+2661120*(2*m0sq*B12(5,8)+(m0sq-m1sq+p1
     &  sq)*B12(5,9)))
       B12(5,10)=Inv14636160*(66*(m0sq4+3*m0sq3*m1sq+6*m0sq2*m1sq2+10*m0
     &  sq*m1sq3+15*m1sq4)-99*(m0sq3+3*m0sq2*m1sq+5*m0sq*m1sq2+5*m1sq3)*
     &  p1sq+33*(2*m0sq2+5*m0sq*m1sq+5*m1sq2)*p1sq2-11*(2*m0sq+3*m1sq)*p
     &  1sq3+3*p1sq4+665280*(A602+2*m0sq*B12(4,8)+(m0sq-m1sq+p1sq)*B12(4
     &  ,9)))
       B12(4,10)=Inv2439360*(33*(m0sq3+5*m0sq2*m1sq+15*m0sq*m1sq2+35*m1s
     &  q3)-55*(m0sq2+4*m0sq*m1sq+7*m1sq2)*p1sq+11*(3*m0sq+7*m1sq)*p1sq2
     &  -7*p1sq3+110880*(A402+2*m0sq*B12(3,8)+(m0sq-m1sq+p1sq)*B12(3,9))
     &  )
       B12(3,10)=Inv609840*(55*(m0sq2+7*m0sq*m1sq+28*m1sq2)-77*(m0sq+4*m
     &  1sq)*p1sq+28*p1sq2+27720*(A202+2*m0sq*B12(2,8)+(m0sq-m1sq+p1sq)*
     &  B12(2,9)))
       B12(2,10)=Inv10890*(11*(m0sq+9*m1sq)-9*p1sq+495*(A2+2*m0sq*B12(1,
     &  8)+(m0sq-m1sq+p1sq)*B12(1,9)))
       B12(1,10)=Inv2*Invp1sq*(A2-(m0sq-m1sq+p1sq)*B12(1,9)-18*B12(2,10)
     &  )
       B12(6,11)=Inv127733760*(-132*(m0sq5+2*m0sq4*m1sq+3*m0sq3*m1sq2+4*
     &  m0sq2*m1sq3+5*m0sq*m1sq4+6*m1sq5)+33*(5*m0sq4+12*m0sq3*m1sq+18*m
     &  0sq2*m1sq2+20*m0sq*m1sq3+15*m1sq4)*p1sq-22*(5*m0sq3+12*m0sq2*m1s
     &  q+15*m0sq*m1sq2+10*m1sq3)*p1sq2+22*(2*m0sq2+4*m0sq*m1sq+3*m1sq2)
     &  *p1sq3-2*(5*m0sq+6*m1sq)*p1sq4+p1sq5+5322240*(-A802+2*m0sq*B12(5
     &  ,9)+(m0sq-m1sq+p1sq)*B12(5,10)))
       B12(5,11)=Inv63866880*(-99*(m0sq4+4*m0sq3*m1sq+10*m0sq2*m1sq2+20*
     &  m0sq*m1sq3+35*m1sq4)+44*(4*m0sq3+15*m0sq2*m1sq+30*m0sq*m1sq2+35*
     &  m1sq3)*p1sq-66*(2*m0sq2+6*m0sq*m1sq+7*m1sq2)*p1sq2+12*(4*m0sq+7*
     &  m1sq)*p1sq3-7*p1sq4+2661120*(-A602+2*m0sq*B12(4,9)+(m0sq-m1sq+p1
     &  sq)*B12(4,10)))
       B12(4,11)=Inv7983360*(-55*(m0sq3+6*m0sq2*m1sq+21*m0sq*m1sq2+56*m1
     &  sq3)+33*(3*m0sq2+14*m0sq*m1sq+28*m1sq2)*p1sq-21*(3*m0sq+8*m1sq)*
     &  p1sq2+14*p1sq3+332640*(-A402+2*m0sq*B12(3,9)+(m0sq-m1sq+p1sq)*B1
     &  2(3,10)))
       B12(3,11)=Inv190080*(-11*(m0sq2+8*m0sq*m1sq+36*m1sq2)+8*(2*m0sq+9
     &  *m1sq)*p1sq-6*p1sq2+7920*(-A202+2*m0sq*B12(2,9)+(m0sq-m1sq+p1sq)
     &  *B12(2,10)))
       B12(2,11)=Inv7920*(-6*(m0sq+10*m1sq)+5*p1sq+330*(-A2+2*m0sq*B12(1
     &  ,9)+(m0sq-m1sq+p1sq)*B12(1,10)))
       B12(1,11)=-(Inv2*Invp1sq*(A2+(m0sq-m1sq+p1sq)*B12(1,10)+20*B12(2,
     &  11)))
       else
       print*,'no yet implemented'
       stop
       Endif
       return
       End
