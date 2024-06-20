       subroutine tens_red2_R3_G_QUAD(m0,m1,p1sq,
     &  A2,A2R,A2I,A1,A1R,A1I,B012,B012R,B012I,BijR,BijI)
       implicit none
       real*16 p1sq,m0,m1,m0sq,m1sq,msq
       complex*32 B012,A2,A1
       real*16 B012R,B012I, A0R,A0I,A2R,A2I,A1R,A1I
       real*16 BijR(2,3)
       real*16 BijI(2,3)
       complex*32 B12(2,3)
       real*16 Inv2,Inv3,Inv4,Inv6,Inv18,Inv36,Inv48,Invp1sq
       parameter(Inv2=1q0/2q0,Inv3=1q0/3q0, Inv6=Inv2*Inv3)
       parameter(Inv36=Inv6*Inv6,Inv4=1q0/4q0, Inv48=1q0/48q0)
       real*16 m0m1,m0m12,m0m13,m0m14,Logm1m0
       real*16 DerB0R, Der2B0R,DerB1R
     
         B012R=REAL(B012)
         B012I=IMAG(B012)

c  m1=m2
       If((abs(m0-m1).lt.1q-6).or.(abs(m0/m1-1q0).lt.1q-6)) then
          msq=m0*m0
          A0R=A2R
          A0I=A2I
 

c   p1sq>0
          If(abs(p1sq).gt.1q-7) then
c         
          if(msq.gt.1q-7) then   

          Inv18=Inv3*Inv6/p1sq
   
          BijR(1,1)=-(B012R*Inv2)
          BijR(1,2)= Inv18*(p1sq - 6q0*(msq - A0R + (msq - p1sq)*B012R))
          BijR(2,2)= Inv36*(12q0*msq - 2q0*p1sq + 6q0*A0R + 3q0*(4q0*msq - p1sq)*B012R)
     
          BijI(1,1)=-(B012I*Inv2)
          BijI(1,2)= Inv18*(-6q0*(-A0I + (msq - p1sq)*B012I))
          BijI(2,2)= Inv36*(+ 6q0*A0I + 3q0*(4q0*msq - p1sq)*B012I)
   
         B12(2,3)=Inv48*(-6*A2-2*msq-4*msq+p1sq
     & +12*msq*CMPLX(BijR(1,1),BijI(1,1),16)
     &  +6*(p1sq)*CMPLX(BijR(1,2),BijI(1,2),16))

          BijR(2,3)= REAL(B12(2,3))
          BijI(2,3)= IMAG(B12(2,3))

          return
          else
          Inv18=Inv3*Inv6
          
   
          BijR(1,1)=-(B012R*Inv2)
          BijR(1,2)= Inv18*(1q0 + 6q0*B012R)
          BijR(2,2)=-Inv36*p1sq*(2q0+3q0*B012R)
     
          BijI(1,1)=-(B012I*Inv2)
          BijI(1,2)= Inv3*(B012I)
          BijI(2,2)=-Inv36*3q0*p1sq*B012I

         B12(2,3)=Inv48*p1sq*(1q0+6*CMPLX(BijR(1,2),BijI(1,2),16))

          BijR(2,3)= REAL(B12(2,3))
          BijI(2,3)= IMAG(B12(2,3))


          endif
c p1sq=0
          else
   
          BijR(1,1)=-(B012R*Inv2)
          BijR(1,2)= Inv3*B012R
          BijR(2,2)= Inv6*(2q0*(B012R +1q0)*msq + A0R)
   
   
          BijI(1,1)=-(B012I*Inv2)
          BijI(1,2)= Inv3*B012I
          BijI(2,2)= Inv6*(2q0*(B012I)*msq + A0I)

         B12(2,3)=Inv48*6*(-A2-msq+2*msq*CMPLX(BijR(1,1),BijI(1,1),16))

          BijR(2,3)= REAL(B12(2,3))
          BijI(2,3)= IMAG(B12(2,3))

          return
          Endif
c  m1=!=m2
       else

          m0sq=m0*m0
          m1sq=m1*m1
c p1sq>0
          If(abs(p1sq).gt.1q-7) then
         
          inv18=1q0/(18q0)
          invp1sq=1q0/p1sq

          B12(1,1)=-(Inv2*Invp1sq*(-A1+A2+B012*(m0sq-m1sq+p1sq)))
          B12(2,2)=Inv18*(3q0*A2+6q0*B012*m0sq+3q0*(m0sq+m1sq)
     &             -p1sq+3q0*(m0sq-m1sq+p1sq)*B12(1,1))
          B12(1,2)=Inv2*Invp1sq*(A2-(m0sq-m1sq+p1sq)*B12(1,1)-2q0*B12(2,2))

       B12(2,3)=Inv48*(-6*A2-2*m0sq-4*m1sq+p1sq
     & +12*m0sq*B12(1,1)
     & +6*(m0sq-m1sq+p1sq)*B12(1,2))


          BijR(1,1)= REAL(B12(1,1))
          BijR(1,2)= REAL(B12(1,2))
          BijR(2,2)= REAL(B12(2,2))

          BijI(1,1)= IMAG(B12(1,1))
          BijI(1,2)= IMAG(B12(1,2))
          BijI(2,2)= IMAG(B12(2,2))

          BijR(2,3)= REAL(B12(2,3))
          BijI(2,3)= IMAG(B12(2,3))


          return
c p1sq=0
          else
          
           inv18=1q0/(18q0)
c m1 & m2> 0
                If((m0sq.gt.1q-7).and.(m1sq.gt.1q-7)) then  ! m0sq>0 & m1sq>0
                  m0m1=m0sq-m1sq
                  m0m12=m0m1*m0m1
                  m0m13=m0m12*m0m1
                  m0m14=m0m12*m0m12
                  Logm1m0=log(m1sq/m0sq)
  

                  DerB0R=(m0sq+m1sq)/(2q0*m0m12)+m0sq*m1sq/m0m13*logm1m0

                  Der2B0R=(m0sq*m0sq+10q0*m0sq*m1sq+m1sq*m1sq)/(3q0*m0m14) 
     &                    +2q0*m0sq*m1sq*(m0sq+m1sq)/(m0m14*m0m1)*logm1m0

                  DerB1R= -Inv2*DerB0R- m0m1*Inv4*Der2B0R

c m1 or m2 = 0
                else   ! m0sq=0 or m1sq=0
                  m0m1=m0sq-m1sq
                  m0m12=m0m1*m0m1
                  m0m14=m0m12*m0m12
                  DerB0R=(m0sq+m1sq)/(2q0*m0m12)

                  Der2B0R=(m0sq*m0sq+10*m0sq*m1sq+m1sq*m1sq)/(3q0*m0m14) 

                  DerB1R= -Inv2*DerB0R- m0m1*Inv4*Der2B0R

                endif

               BijR(1,1)=-Inv2*(B012r + m0m1*DerB0r)
               BijI(1,1)=-Inv2*B012I
               B12(1,1)=CMPLX(BijR(1,1),BijI(1,1),16)

               B12(2,2)=Inv18*(3q0*A2+6*B012*m0sq+3q0*(m0sq+m1sq)
     &                  -p1sq+3q0*(m0sq-m1sq+ p1sq)*B12(1,1))
               B12(1,2)=Inv18*(-6q0*DerB0R*m0sq-12q0*b12(1,1)
     &                         -12q0*m0m1*(m0sq+m1sq)*DerB1R+1q0)

       B12(2,3)=Inv48*(-6*A2-2*m0sq-4*m1sq+p1sq
     & +12*m0sq*B12(1,1)
     & +6*(m0sq-m1sq+p1sq)*B12(1,2))


          BijR(1,2)= REAL(B12(1,2))
          BijR(2,2)= REAL(B12(2,2))


          BijI(1,2)= IMAG(B12(1,2))
          BijI(2,2)= IMAG(B12(2,2))

          BijR(2,3)= REAL(B12(2,3))
          BijI(2,3)= IMAG(B12(2,3))


         endif


       endif

       return
       End

