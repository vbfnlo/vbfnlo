       subroutine tens_red2_Re(p1sq,B012,B012R,B012I,BijR,BijI)
       implicit none
       Real*8 p1sq
       Complex*16 B012
       Real*8 B012R,B012I
       Real*8 BijR(2,2)
       Real*8 BijI(2,2)
       Real*8 Inv2,Inv18,Inv36

       B012R=Dble(B012)
       B012I=DImag(B012)

       If(abs(p1sq).gt.1d-7) then

        Inv2=1d0/2d0
        Inv18=1d0/18d0
        Inv36=1d0/36d0

       BijR(1,1)=-(B012R*Inv2)
       BijR(1,2)=(1d0+6d0*B012R)*Inv18
       BijR(2,2)=-((2d0+3d0*B012R)*Inv36*p1sq)

       BijI(1,1)=-(B012I*Inv2)
       BijI(1,2)=(6d0*B012I)*Inv18
       BijI(2,2)=-((3d0*B012I)*Inv36*p1sq)

      return
       else
       BijR(1,1)=0d0
       BijR(1,2)=0d0
       BijR(2,2)=0d0

       BijI(1,1)=0d0
       BijI(1,2)=0d0
       BijI(2,2)=0d0

       Endif

       return
       End
