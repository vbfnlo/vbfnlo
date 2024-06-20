       subroutine tens_red2_R2_GDiv(m0,m1,p1sq,
     &  A2,A2R,A2I,A1,A1R,A1I,B012,B012R,B012I,BijR,BijI)
       implicit none
       Real*8 p1sq,m0,m1,m0sq,m1sq,msq
       Complex*16 B012,A2,A1
       Real*8 B012R,B012I, A0R,A0I,A2R,A2I,A1R,A1I
       Real*8 BijR(2,2)
       Real*8 BijI(2,2)
       complex*16 B12(2,2)
       Real*8 Inv2,Inv3,Inv4,Inv6,Inv18,Inv36,Invp1sq
       parameter(Inv2=1d0/2d0,Inv3=1d0/3d0, Inv6=Inv2*Inv3)
       parameter(Inv36=Inv6*Inv6,Inv4=1d0/4d0)
       real*8 DerB0R, Der2B0R,DerB1R
       
         B012R=Dble(B012)
         B012I=DImag(B012)

c  m1=m2
       If((abs(m0-m1).lt.1d-6).or.(m0/m1-1d0.lt.1d-6)) then
   
          msq=m0*m0
          A0R=A2R
          A0I=A2I
 

c   p1sq>0
          If(abs(p1sq).gt.1d-7) then
c         
          if(msq.gt.1d-7) then   

          Inv18=Inv3*Inv6/p1sq
   
          BijR(1,1)=-(B012R*Inv2)
          BijR(1,2)= Inv18*( - 6d0*( - A0R + (msq - p1sq)*B012R))
          BijR(2,2)= Inv36*(  6d0*A0R + 3d0*(4d0*msq - p1sq)*B012R)
     
          BijI(1,1)=-(B012I*Inv2)
          BijI(1,2)= Inv18*(-6d0*(-A0I + (msq - p1sq)*B012I))
          BijI(2,2)= Inv36*(+ 6d0*A0I + 3d0*(4d0*msq - p1sq)*B012I)
   
          return
          else
          Inv18=Inv3*Inv6
          
   
          BijR(1,1)=-(B012R*Inv2)
          BijR(1,2)= Inv18*(+ 6d0*B012R)
          BijR(2,2)=-Inv36*p1sq*(2d0+3d0*B012R)
     
          BijI(1,1)=-(B012I*Inv2)
          BijI(1,2)= Inv3*(B012I)
          BijI(2,2)=-Inv36*3d0*p1sq*B012I


          endif
c p1sq=0
          else
   
          BijR(1,1)=-(B012R*Inv2)
          BijR(1,2)= Inv3*B012R
          BijR(2,2)= Inv6*(2d0*(B012R)*msq + A0R)
   
   
          BijI(1,1)=-(B012I*Inv2)
          BijI(1,2)= Inv3*B012I
          BijI(2,2)= Inv6*(2d0*(B012I)*msq + A0I)
          return
          Endif
c  m1=!=m2
       else
          m0sq=m0*m0
          m1sq=m1*m1
c p1sq>0
          If(abs(p1sq).gt.1d-7) then
         
          inv18=1d0/(18d0)
          invp1sq=1d0/p1sq

          B12(1,1)=-(Inv2*Invp1sq*(-A1+A2+B012*(m0sq-m1sq+p1sq)))
          B12(2,2)=Inv18*(3d0*A2+6d0*B012*m0sq
     &                    +3d0*(m0sq-m1sq+p1sq)*B12(1,1))
          B12(1,2)=Inv2*Invp1sq*(A2-(m0sq-m1sq+p1sq)*B12(1,1)-2d0*B12(2,2))

          BijR(1,1)= Dble(B12(1,1))
          BijR(1,2)= Dble(B12(1,2))
          BijR(2,2)= Dble(B12(2,2))

          BijI(1,1)= DImag(B12(1,1))
          BijI(1,2)= DImag(B12(1,2))
          BijI(2,2)= DImag(B12(2,2))
          return
c p1sq=0
          else
          
           inv18=1d0/(18d0)
c m1 & m2> 0
                If((m0sq.gt.1d-7).and.(m1sq.gt.1d-7)) then  ! m0sq>0 & m1sq>0
                    DerB0R=0
                  Der2B0R=0
                  DerB1R= 0

c m1 or m2 = 0
                else   ! m0sq=0 or m1sq=0
                    DerB0R=0
                   Der2B0R=0
                   DerB1R=0

                endif

               BijR(1,1)=-Inv2*(B012r)
               BijI(1,1)=-Inv2*B012I
               B12(1,1)=DCMPLX(BijR(1,1),BijI(1,1))

               B12(2,2)=Inv18*(3d0*A2+6*B012*m0sq
     &                  +3d0*(m0sq-m1sq+ p1sq)*B12(1,1))
               B12(1,2)=Inv18*(-12d0*b12(1,1)
     &                         )


          BijR(1,2)= Dble(B12(1,2))
          BijR(2,2)= Dble(B12(2,2))


          BijI(1,2)= DImag(B12(1,2))
          BijI(2,2)= DImag(B12(2,2))


         endif


       endif

       return
       End
