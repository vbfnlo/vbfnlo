       subroutine tens_red2_Re_1MDiv(m,p1sq,
     &  A0R,A0I,B012,B012R,B012I,BijR,BijI)
       implicit none
       Real*8 p1sq
       Complex*16 B012
       Real*8 B012R,B012I
       Real*8 BijR(2,2)
       Real*8 BijI(2,2)
       Real*8 Inv2,Inv3,Inv6,m,msq,Inv18,Inv36
       Real*8 A0R,A0I

       msq=m*m
       Inv2=1d0/2d0
       Inv3=1d0/3d0
       Inv6=Inv2*Inv3
       

       B012R=Dble(B012)
       B012I=DImag(B012)

       If(abs(p1sq).gt.1d-8) then
        
       Inv18=Inv3*Inv6/p1sq
       Inv36=Inv6*Inv6

       BijR(1,1)=-(B012R*Inv2)
       BijR(1,2)=Inv18*(-6d0*(-A0R + (msq - p1sq)*B012R))
       BijR(2,2)=Inv36*(6d0*A0R + 3d0*(4d0*msq - p1sq)*B012R)


       BijI(1,1)=-(B012I*Inv2)
       BijI(1,2)=Inv18*(-6d0*(-A0I + (msq - p1sq)*B012I))
       BijI(2,2)=Inv36*(6d0*A0I + 3d0*(4d0*msq - p1sq)*B012I)

      return
       else

       BijR(1,1)=-(B012R*Inv2)
       BijR(1,2)= Inv3*B012R
       BijR(2,2)=Inv6*(2d0*(B012R)*msq + A0R)


       BijI(1,1)=-(B012I*Inv2)
       BijI(1,2)= Inv3*B012I
       BijI(2,2)= Inv6*(2d0*(B012I)*msq + A0I)

       Endif

       return
       End
