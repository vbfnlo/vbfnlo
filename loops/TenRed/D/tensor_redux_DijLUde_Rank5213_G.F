c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables
c************************************************************************************
c***********************************************************************************
       subroutine ten4_213_R5_G(Div,m0sq,m1sq,m2sq,m3sq,p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,
     &  r1,r2r1,r3r2,
     &  z11,z12,z13,z21,z22,z23,z31,z32,z33,iz11,iz22,iz33,
     &                     C0234R, C0134R, C0124R, C0123R,
     &                     C234R, C134R, C124R, C123R,
     &                     C0234I, C0134I, C0124I, C0123I,
     &                     C234I, C134I, C124I, C123I,
     &                     D0R,D0I,DijR,DijI)

       implicit none
       real*8 p1sq, p2sq, p3sq, p1p2, p1p3,p2p3,p4sq,s12,s23
       real*8 r1, r2r1, r3r2,Iv
       real*8 z11,z12,z13,z21,z22,z23,z31,z32,z33
       real*8 iz11,iz22,iz33
       real*8 PRr(1:3), PRI(1:3)
       real*8 Rr(1:3), RI(1:3)
       real*8 C0234R, C0134R, C0124R, C0123R
       real*8 C234R(9,4), C134R(9,4), C124R(9,4), C123R(9,4)
       real*8 C0234I, C0134I, C0124I, C0123I
       real*8 C234I(9,4), C134I(9,4), C124I(9,4), C123I(9,4)
c       real*8 D0r, Dijr(70,7)
c       real*8 D0I, DijI(70,7)
       real*8 D0r, Dijr(34,5)
       real*8 D0I, DijI(34,5)
       real*8 r10,r21,r32
       real*8 Inv6, Inv36,Inv2400,Inv10,Inv14,Inv282240

       real*8 m0sq,m1sq,m2sq,m3sq
       integer Div
       Inv6=1d0/6d0
       Inv10=1d-1
       Inv36=1d0/36d0
       Inv2400=1d0/24d2
       inv14=1d0/14d0
       Inv282240=1d0/282240d0


       r10=r1
       r21=r2r1
       r32=r3r2


       s12 = (p1sq +p2sq+ 2*p1p2)
       s23 = (p2sq +p3sq+ 2*p2p3)

       p4sq=p1sq+p2sq+p3sq+2d0*(p1p2+p1p3+p2p3)


       PRr(1) = (C0134R - C0234R - D0r*r1   )
       PRr(2) = (C0124R - C0134R - D0r*r2r1 )
       PRr(3) = (C0123R - C0124R - D0r*r3r2 )


           Rr(3)=(PRr(3)-z32*(PRr(1)-z21*PRr(2))-z31*PRr(2))*iz33
           Rr(2)=(PRr(1)-z21*PRr(2) - z23*Rr(3))*iz22
           Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11


       Dijr(1,1) = Rr(1)
       Dijr(2,1) = Rr(2)
       Dijr(3,1) = Rr(3)
c
       PRI(1) = (C0134I - C0234I - D0I*r1   )
       PRI(2) = (C0124I - C0134I - D0I*r2r1 )
       PRI(3) = (C0123I - C0124I - D0I*r3r2 )
c
           RI(3)=(PRI(3)-z32*(PRI(1)-z21*PRI(2)) -z31*PRI(2))*iz33
           RI(2)=(PRI(1)-z21*PRI(2) - z23*RI(3))*iz22
           RI(1)=(PRI(2)-  z12*RI(2)-  z13*RI(3))*iz11

       DijI(1,1) = RI(1)
       DijI(2,1) = RI(2)
       DijI(3,1) = RI(3)
CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c D00
CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       Dijr(7,2) = (C0234R + 2d0*m0sq*D0r+ Dijr(1,1)*r1 + Dijr(2,1)*r2r1
     &      +Dijr(3,1)*r3r2)*0.5d0

       DijI(7,2) = (C0234I+2d0*m0sq*D0I + DijI(1,1)*r1 + DijI(2,1)*r2r1
     &      +DijI(3,1)*r3r2)*0.5d0


CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c 30-3 1-32

c
       PRr(1) = (C134R(1,1) + C0234R - Dijr(1,1)*r1 - Dijr(7,2)*2.d0)
       PRr(2) = (C124R(1,1) - C134R(1,1) - Dijr(1,1)*r2r1        )
       PRr(3) = (C123R(1,1) - C124R(1,1) - Dijr(1,1)*r3r2        )
c
           Rr(3)=(PRr(3)-z32*(PRr(1)-z21*PRr(2))-z31*PRr(2))*iz33
           Rr(2)=(PRr(1)-z21*PRr(2) - z23*Rr(3))*iz22
           Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11

       Dijr(1,2) = Rr(1)
       Dijr(4,2) = Rr(2)
       Dijr(5,2) = Rr(3)
c
       PRI(1) = (C134I(1,1) + C0234I  - DijI(1,1)*r1 - DijI(7,2)*2.d0)
       PRI(2) = (C124I(1,1) - C134I(1,1) - DijI(1,1)*r2r1        )
       PRI(3) = (C123I(1,1) - C124I(1,1) - DijI(1,1)*r3r2        )
c
       RI(3)=(PRI(3)-z32*(PRI(1)-z21*PRI(2)) -z31*PRI(2))*iz33
           RI(2)=(PRI(1)-z21*PRI(2) - z23*RI(3))*iz22
           RI(1)=(PRI(2)-  z12*RI(2)-  z13*RI(3))*iz11


       DijI(1,2) = RI(1)
       DijI(4,2) = RI(2)
       DijI(5,2) = RI(3)

c 33-3 4-35
       PRr(1) = (C134R(1,1) -C234R(1,1) - Dijr(2,1)*r1          )
       PRr(2)=(C124R(2,1)-C134R(1,1)-Dijr(2,1)*r2r1-Dijr(7,2)*2.d0)
       PRr(3) = (C123R(2,1) -C124R(2,1) - Dijr(2,1)*r3r2        )

       Rr(3)=(PRr(3)-z32*(PRr(1)-z21*PRr(2))-z31*PRr(2))*iz33
           Rr(2)=(PRr(1)-z21*PRr(2) - z23*Rr(3))*iz22
c          Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11

c       Dijr(4,2) = Rr(1)
       Dijr(2,2) = Rr(2)
       Dijr(6,2) = Rr(3)

       PRI(1) = (C134I(1,1) -C234I(1,1) - DijI(2,1)*r1          )
       PRI(2)=(C124I(2,1)-C134I(1,1)-DijI(2,1)*r2r1-DijI(7,2)*2.d0)
       PRI(3) = (C123I(2,1) -C124I(2,1) - DijI(2,1)*r3r2        )

       RI(3)=(PRI(3)-z32*(PRI(1)-z21*PRI(2)) -z31*PRI(2))*iz33
           RI(2)=(PRI(1)-z21*PRI(2) - z23*RI(3))*iz22
c          RI(1)=(PRI(2)-  z12*RI(2)-  z13*RI(3))*iz11

c       DijI(4,2) = RI(1)
       DijI(2,2) = RI(2)
       DijI(6,2) = RI(3)

c 36-3 7-38

       PRr(1) = (C134R(2,1) - C234R(2,1) - Dijr(3,1)*r1      )
       PRr(2) = (C124R(2,1) - C134R(2,1) - Dijr(3,1)*r2r1    )
       PRr(3) = (       - C124R(2,1) - Dijr(3,1)*r3r2 -Dijr(7,2)*2.d0)

       Rr(3)=(PRr(3)-z32*(PRr(1)-z21*PRr(2))-z31*PRr(2))*iz33
c          Rr(2)=(PRr(1)-z21*PRr(2) - z23*Rr(3))*iz22
c          Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11

c       Dijr(5,2) = Rr(1)
c       Dijr(6,2) = Rr(2)
       Dijr(3,2) = Rr(3)

       PRI(1) = (C134I(2,1) - C234I(2,1) - DijI(3,1)*r1      )
       PRI(2) = (C124I(2,1) - C134I(2,1) - DijI(3,1)*r2r1    )
       PRI(3) = (       - C124I(2,1) - DijI(3,1)*r3r2 -DijI(7,2)*2.d0)


       RI(3)=(PRI(3)-z32*(PRI(1)-z21*PRI(2)) -z31*PRI(2))*iz33
c          RI(2)=(PRI(1)-z21*PRI(2) - z23*RI(3))*iz22
c          RI(1)=(PRI(2)-  z12*RI(2)-  z13*RI(3))*iz11

c       DijI(5,2) = RI(1)
c       DijI(6,2) = RI(2)
       DijI(3,2) = RI(3)


CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       Dijr(11,3)=0.25d0*(-C0234R+2d0*m0sq*Dijr(1,1)+Dijr(1,2)*r1
     & +Dijr(4,2)*r2r1 +Dijr(5,2)*r3r2)

       Dijr(12,3)=0.25d0*(C234R(1,1)+2d0*m0sq*Dijr(2,1)+Dijr(4,2)*r1
     &   +Dijr(2,2)*r2r1 +Dijr(6,2)*r3r2)
       Dijr(13,3)=0.25d0*(C234R(2,1)+2d0*m0sq*Dijr(3,1)+Dijr(5,2)*r1
     &    + Dijr(6,2)*r2r1 + Dijr(3,2)*r3r2)

        DijI(11,3)=0.25d0*(-C0234I+2d0*m0sq*DijI(1,1)+DijI(1,2)*r1
     &   +DijI(4,2)*r2r1 +DijI(5,2)*r3r2)

       DijI(12,3)=0.25d0*(C234I(1,1)+2d0*m0sq*DijI(2,1)+DijI(4,2)*r1
     &   +DijI(2,2)*r2r1 +DijI(6,2)*r3r2)
       DijI(13,3)=0.25d0*(C234I(2,1)+2d0*m0sq*DijI(3,1)+DijI(5,2)*r1
     &  + DijI(6,2)*r2r1 + DijI(3,2)*r3r2)

c 41-4 2-43

       PRr(1) =(C134R(1,2) - C0234R - Dijr(1,2)*r1  - Dijr(11,3)*4.d0)
       PRr(2) =(C124R(1,2) - C134R(1,2) - Dijr(1,2)*r2r1       )
       PRr(3) =(C123R(1,2) - C124R(1,2) - Dijr(1,2)*r3r2       )

       Rr(3)=(PRr(3)-z32*(PRr(1)-z21*PRr(2))-z31*PRr(2))*iz33
           Rr(2)=(PRr(1)-z21*PRr(2) - z23*Rr(3))*iz22
           Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11

       Dijr(1,3) = Rr(1)
       Dijr(4,3) = Rr(2)
       Dijr(5,3) = Rr(3)

       PRI(1) =(C134I(1,2) - C0234I - DijI(1,2)*r1  - DijI(11,3)*4.d0)
       PRI(2) =(C124I(1,2) - C134I(1,2) - DijI(1,2)*r2r1       )
       PRI(3) =(C123I(1,2) - C124I(1,2) - DijI(1,2)*r3r2       )

       RI(3)=(PRI(3)-z32*(PRI(1)-z21*PRI(2)) -z31*PRI(2))*iz33
           RI(2)=(PRI(1)-z21*PRI(2) - z23*RI(3))*iz22
           RI(1)=(PRI(2)-  z12*RI(2)-  z13*RI(3))*iz11

       DijI(1,3) = RI(1)
       DijI(4,3) = RI(2)
       DijI(5,3) = RI(3)

c 50-5 1-52

       PRr(1) =  (C134R(1,2) - C234R(1,2) - Dijr(2,2)*r1)
       PRr(2) =  (C124R(2,2) - C134R(1,2) - Dijr(2,2)*r2r1
     &       -Dijr(12,3)*4.d0)
       PRr(3) =  (C123R(2,2) - C124R(2,2) - Dijr(2,2)*r3r2)

       Rr(3)=(PRr(3)-z32*(PRr(1)-z21*PRr(2))-z31*PRr(2))*iz33
           Rr(2)=(PRr(1)-z21*PRr(2) - z23*Rr(3))*iz22
           Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11

       Dijr(6,3) = Rr(1)
       Dijr(2,3) = Rr(2)
       Dijr(8,3) = Rr(3)

       PRI(1) =  (C134I(1,2) - C234I(1,2) - DijI(2,2)*r1)
       PRI(2) =  (C124I(2,2) - C134I(1,2) - DijI(2,2)*r2r1
     &      -DijI(12,3)*4.d0)
       PRI(3) =  (C123I(2,2) - C124I(2,2) - DijI(2,2)*r3r2)


       RI(3)=(PRI(3)-z32*(PRI(1)-z21*PRI(2)) -z31*PRI(2))*iz33
           RI(2)=(PRI(1)-z21*PRI(2) - z23*RI(3))*iz22
           RI(1)=(PRI(2)-  z12*RI(2)-  z13*RI(3))*iz11

       DijI(6,3) = RI(1)
       DijI(2,3) = RI(2)
       DijI(8,3) = RI(3)

c
c 56-5 7-58
c
c
       PRr(1) = (C134R(2,2) - C234R(2,2) - Dijr(3,2)*r1)
       PRr(2) = (C124R(2,2) - C134R(2,2) - Dijr(3,2)*r2r1)
       PRr(3) = (     -C124R(2,2) - Dijr(3,2)*r3r2 - Dijr(13,3)*4.d0)
c

       Rr(3)=(PRr(3)-z32*(PRr(1)-z21*PRr(2))-z31*PRr(2))*iz33
           Rr(2)=(PRr(1)-z21*PRr(2) - z23*Rr(3))*iz22
           Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11

       Dijr(7,3) = Rr(1)
       Dijr(9,3) = Rr(2)
       Dijr(3,3) = Rr(3)
c
       PRI(1) = (C134I(2,2) - C234I(2,2) - DijI(3,2)*r1)
       PRI(2) = (C124I(2,2) - C134I(2,2) - DijI(3,2)*r2r1)
       PRI(3) = (     -C124I(2,2) - DijI(3,2)*r3r2 - DijI(13,3)*4.d0)
c

       RI(3)=(PRI(3)-z32*(PRI(1)-z21*PRI(2)) -z31*PRI(2))*iz33
           RI(2)=(PRI(1)-z21*PRI(2) - z23*RI(3))*iz22
           RI(1)=(PRI(2)-  z12*RI(2)-  z13*RI(3))*iz11

       DijI(7,3) = RI(1)
       DijI(9,3) = RI(2)
       DijI(3,3) = RI(3)
c 44-4 5-46
       PRr(1) = (C134R(1,2) + C234R(1,1) - Dijr(4,2)*r1
     &      -Dijr(12,3)*2.d0)
       PRr(2) = (C124R(3,2) - C134R(1,2) - Dijr(4,2)*r2r1
     &      - Dijr(11,3)*2.d0)
       PRr(3) = (C123R(3,2) - C124R(3,2)  -Dijr(4,2)*r3r2)
c

       Rr(3)=(PRr(3)-z32*(PRr(1)-z21*PRr(2))-z31*PRr(2))*iz33
c          Rr(2)=(PRr(1)-z21*PRr(2) - z23*Rr(3))*iz22
c          Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11

C       Dijr(4,3) = Rr(1)
C       Dijr(6,3) = Rr(2)
       Dijr(10,3) = Rr(3)
c
c
       PRI(1) = (C134I(1,2) + C234I(1,1) - DijI(4,2)*r1
     &      -DijI(12,3)*2.d0)
       PRI(2) = (C124I(3,2) - C134I(1,2) - DijI(4,2)*r2r1
     &      - DijI(11,3)*2.d0)
       PRI(3) = (C123I(3,2) - C124I(3,2)  -DijI(4,2)*r3r2)


       RI(3)=(PRI(3)-z32*(PRI(1)-z21*PRI(2)) -z31*PRI(2))*iz33
c          RI(2)=(PRI(1)-z21*PRI(2) - z23*RI(3))*iz22
c          RI(1)=(PRI(2)-  z12*RI(2)-  z13*RI(3))*iz11

C       DijI(4,3) = RI(1)
C       DijI(6,3) = RI(2)
       DijI(10,3) = RI(3)
c
c to d efine D00ij and D0000 functions!
c In P V notation we have:
c Dij( 1,7)=D416, Dij(1,8)=D417, Dij(1,9)=D418
c Dij( 1,10)=D419, Dij(1,11)=D420, Dij(1,12)=D421, D(1,13)=D422
c
       Iv=1.d0/6.d0
c

c
       Dijr(7,1)=(C0234R+2*m0sq*Dijr(1,2)+Dijr(1,3)*r1+Dijr(4,3)*r2r1+Dijr(5,3)*r3r2 )*Iv
c
       Dijr(8,1)=(C234R(1,2)+2*m0sq*Dijr(2,2)+Dijr(6,3)*r1+Dijr(2,3)*r2r1
     &  +Dijr(8,3)*r3r2)*Iv
c
       Dijr(9,1)=(C234R(2,2)+2*m0sq*Dijr(3,2)+Dijr(7,3)*r1+Dijr(9,3)*r2r1
     &  +Dijr(3,3)*r3r2)*Iv
c
       Dijr(10,1)=(-C234R(1,1)+2*m0sq*Dijr(4,2)+Dijr(4,3)*r1+Dijr(6,3)*r2r1
     &  +Dijr(10,3)*r3r2)*Iv
c
       Dijr(11,1)=(-C234R(2,1)+2*m0sq*Dijr(5,2)+Dijr(5,3)*r1+Dijr(10,3)*r2r1
     &  +Dijr(7,3)*r3r2)*Iv
c
       Dijr(12,1)=(C234R(3,2)+2*m0sq*Dijr(6,2)+Dijr(10,3)*r1+Dijr(8,3)*r2r1
     &  +Dijr(9,3)*r3r2)*Iv
c

       if(Div.eq.0) then
       Dijr(13,1)=(C234R(4,2)+2*m0sq*Dijr(7,2)+Dijr(11,3)*r1+Dijr(12,3)*r2r1
     &       + Dijr(13,3)*r3r2+Iv)*Iv
       else
       Dijr(13,1)=(C234R(4,2)+2*m0sq*Dijr(7,2)+Dijr(11,3)*r1+Dijr(12,3)*r2r1
     &       + Dijr(13,3)*r3r2)*Iv
       endif

       DijI(7,1)=(C0234I+2*m0sq*DijI(1,2)+DijI(1,3)*r1+DijI(4,3)*r2r1+DijI(5,3)*r3r2 )*Iv
c
       DijI(8,1)=(C234I(1,2)+2*m0sq*DijI(2,2)+DijI(6,3)*r1+DijI(2,3)*r2r1
     &  +DijI(8,3)*r3r2)*Iv
c
       DijI(9,1)=(C234I(2,2)+2*m0sq*DijI(3,2)+DijI(7,3)*r1+DijI(9,3)*r2r1
     & +DijI(3,3)*r3r2)*Iv
c
       DijI(10,1)=(-C234I(1,1)+2*m0sq*DijI(4,2)+DijI(4,3)*r1+DijI(6,3)*r2r1
     &  +DijI(10,3)*r3r2)*Iv
c
       DijI(11,1)=(-C234I(2,1)+2*m0sq*DijI(5,2)+DijI(5,3)*r1+DijI(10,3)*r2r1
     &  +DijI(7,3)*r3r2)*Iv
c
       DijI(12,1)=(C234I(3,2)+2*m0sq*DijI(6,2)+DijI(10,3)*r1+DijI(8,3)*r2r1
     &  +DijI(9,3)*r3r2)*Iv
c
       DijI(13,1)=(C234I(4,2)+2*m0sq*DijI(7,2)+DijI(11,3)*r1+DijI(12,3)*r2r1
     &       + DijI(13,3)*r3r2)*Iv

c to d efine D00ij and D0000 functions!
c In P V notation we have:
c Dij( 1,7)=D416, Dij(1,8)=D417, Dij(1,9)=D418
c Dij( 1,10)=D419, Dij(1,11)=D420, Dij(1,12)=D421, D(1,13)=D422


c                {0, 0, 0, 0}
       DijR(22,4)=DijR(13,1)
       DijI(22,4)=DijI(13,1)
c                {0, 0, ii1, 1}

       DijR(16,4)=DijR(7,1)
       DijR(17,4)=DijR(8,1)
       DijR(18,4)=DijR(9,1)
       DijR(19,4)=DijR(10,1)
       DijR(20,4)=DijR(11,1)
       DijR(21,4)=DijR(12,1)

       DijI(16,4)=DijI(7,1)
       DijI(17,4)=DijI(8,1)
       DijI(18,4)=DijI(9,1)
       DijI(19,4)=DijI(10,1)
       DijI(20,4)=DijI(11,1)
       DijI(21,4)=DijI(12,1)


c FC %c                {0, 0, 0, 0}
c FC %       DijR(22,4)=Inv36+Inv6*(C234R(4,2)+r10*DijR(11,3)+r21*DijR(1
c FC %     &   2,3)+r32*DijR(13,3))
c FC %       DijI(22,4)=Inv6*(C234I(4,2)+r10*DijI(11,3)+r21*DijI(12,3)+r
c FC %     &   32*DijI(13,3))
c FC %c                {0, 0, ii1, 1}
c FC %       PRr(1)=C134R(5,3)+C234R(4,2)-r10*DijR(11,3)-2*DijR(22,4)
c FC %       PRr(2)=C124R(5,3)-C134R(5,3)-r21*DijR(11,3)
c FC %       PRr(3)=C123R(5,3)-C124R(5,3)-r32*DijR(11,3)
c FC %       Rr(3)=(PRr(3)-z32*(PRr(1)-z21*PRr(2))-z31*PRr(2))*iz33
c FC %       Rr(2)=(PRr(1)-z21*PRr(2)-z23*Rr(3))*iz22
c FC %       Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11
c FC %       DijR(16,4)=Rr(1)
c FC %       DijR(19,4)=Rr(2)
c FC %       DijR(20,4)=Rr(3)
c FC %       PRI(1)=C134I(5,3)+C234I(4,2)-r10*DijI(11,3)-2*DijI(22,4)
c FC %       PRI(2)=C124I(5,3)-C134I(5,3)-r21*DijI(11,3)
c FC %       PRI(3)=C123I(5,3)-C124I(5,3)-r32*DijI(11,3)
c FC %       RI(3)=(PRI(3)-z32*(PRI(1)-z21*PRI(2)) -z31*PRI(2))*iz33
c FC %       RI(2)=(PRI(1)-z21*PRI(2)-z23*RI(3))*iz22
c FC %       RI(1)=(PRI(2)-  z12*RI(2)-  z13*RI(3))*iz11
c FC %       DijI(16,4)=RI(1)
c FC %       DijI(19,4)=RI(2)
c FC %       DijI(20,4)=RI(3)
c FC %c                {0, 0, ii1, 2}
c FC %       PRr(1)=C134R(5,3)-C234R(5,3)-r10*DijR(12,3)
c FC %       PRr(2)=C124R(6,3)-C134R(5,3)-r21*DijR(12,3)-2*DijR(22,4)
c FC %       PRr(3)=C123R(6,3)-C124R(6,3)-r32*DijR(12,3)
c FC %       Rr(3)=(PRr(3)-z32*(PRr(1)-z21*PRr(2))-z31*PRr(2))*iz33
c FC %       Rr(2)=(PRr(1)-z21*PRr(2)-z23*Rr(3))*iz22
c FC %       Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11
c FC %       DijR(19,4)=Rr(1)
c FC %       DijR(17,4)=Rr(2)
c FC %       DijR(21,4)=Rr(3)
c FC %       PRI(1)=C134I(5,3)-C234I(5,3)-r10*DijI(12,3)
c FC %       PRI(2)=C124I(6,3)-C134I(5,3)-r21*DijI(12,3)-2*DijI(22,4)
c FC %       PRI(3)=C123I(6,3)-C124I(6,3)-r32*DijI(12,3)
c FC %       RI(3)=(PRI(3)-z32*(PRI(1)-z21*PRI(2)) -z31*PRI(2))*iz33
c FC %       RI(2)=(PRI(1)-z21*PRI(2)-z23*RI(3))*iz22
c FC %       RI(1)=(PRI(2)-  z12*RI(2)-  z13*RI(3))*iz11
c FC %       DijI(19,4)=RI(1)
c FC %       DijI(17,4)=RI(2)
c FC %       DijI(21,4)=RI(3)
c FC %c                {0, 0, ii1, 3}
c FC %       PRr(1)=C134R(6,3)-C234R(6,3)-r10*DijR(13,3)
c FC %       PRr(2)=C124R(6,3)-C134R(6,3)-r21*DijR(13,3)
c FC %       PRr(3)=-C124R(6,3)-r32*DijR(13,3)-2*DijR(22,4)
c FC %       Rr(3)=(PRr(3)-z32*(PRr(1)-z21*PRr(2))-z31*PRr(2))*iz33
c FC %       Rr(2)=(PRr(1)-z21*PRr(2)-z23*Rr(3))*iz22
c FC %       Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11
c FC %       DijR(20,4)=Rr(1)
c FC %       DijR(21,4)=Rr(2)
c FC %       DijR(18,4)=Rr(3)
c FC %       PRI(1)=C134I(6,3)-C234I(6,3)-r10*DijI(13,3)
c FC %       PRI(2)=C124I(6,3)-C134I(6,3)-r21*DijI(13,3)
c FC %       PRI(3)=-C124I(6,3)-r32*DijI(13,3)-2*DijI(22,4)
c FC %       RI(3)=(PRI(3)-z32*(PRI(1)-z21*PRI(2)) -z31*PRI(2))*iz33
c FC %       RI(2)=(PRI(1)-z21*PRI(2)-z23*RI(3))*iz22
c FC %       RI(1)=(PRI(2)-  z12*RI(2)-  z13*RI(3))*iz11
c FC %       DijI(20,4)=RI(1)
c FC %       DijI(21,4)=RI(2)
c FC %       DijI(18,4)=RI(3)
c                {ii1, 1, 1, 1}
       PRr(1)=C0234R+C134R(1,3)-r10*DijR(1,3)-6*DijR(16,4)
       PRr(2)=C124R(1,3)-C134R(1,3)-r21*DijR(1,3)
       PRr(3)=C123R(1,3)-C124R(1,3)-r32*DijR(1,3)
       Rr(3)=(PRr(3)-z32*(PRr(1)-z21*PRr(2))-z31*PRr(2))*iz33
       Rr(2)=(PRr(1)-z21*PRr(2)-z23*Rr(3))*iz22
       Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11
       DijR(1,4)=Rr(1)
       DijR(4,4)=Rr(2)
       DijR(5,4)=Rr(3)
       PRI(1)=C0234I+C134I(1,3)-r10*DijI(1,3)-6*DijI(16,4)
       PRI(2)=C124I(1,3)-C134I(1,3)-r21*DijI(1,3)
       PRI(3)=C123I(1,3)-C124I(1,3)-r32*DijI(1,3)
       RI(3)=(PRI(3)-z32*(PRI(1)-z21*PRI(2)) -z31*PRI(2))*iz33
       RI(2)=(PRI(1)-z21*PRI(2)-z23*RI(3))*iz22
       RI(1)=(PRI(2)-  z12*RI(2)-  z13*RI(3))*iz11
       DijI(1,4)=RI(1)
       DijI(4,4)=RI(2)
       DijI(5,4)=RI(3)
c                {ii1, 2, 1, 1}
       PRr(1)=C134R(1,3)-C234R(1,1)-r10*DijR(4,3)-4*DijR(19,4)
       PRr(2)=C124R(3,3)-C134R(1,3)-r21*DijR(4,3)-2*DijR(16,4)
       PRr(3)=C123R(3,3)-C124R(3,3)-r32*DijR(4,3)
       Rr(3)=(PRr(3)-z32*(PRr(1)-z21*PRr(2))-z31*PRr(2))*iz33
       Rr(2)=(PRr(1)-z21*PRr(2)-z23*Rr(3))*iz22
       Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11
       DijR(4,4)=Rr(1)
       DijR(6,4)=Rr(2)
       DijR(7,4)=Rr(3)
       PRI(1)=C134I(1,3)-C234I(1,1)-r10*DijI(4,3)-4*DijI(19,4)
       PRI(2)=C124I(3,3)-C134I(1,3)-r21*DijI(4,3)-2*DijI(16,4)
       PRI(3)=C123I(3,3)-C124I(3,3)-r32*DijI(4,3)
       RI(3)=(PRI(3)-z32*(PRI(1)-z21*PRI(2)) -z31*PRI(2))*iz33
       RI(2)=(PRI(1)-z21*PRI(2)-z23*RI(3))*iz22
       RI(1)=(PRI(2)-  z12*RI(2)-  z13*RI(3))*iz11
       DijI(4,4)=RI(1)
       DijI(6,4)=RI(2)
       DijI(7,4)=RI(3)
c                {ii1, 2, 2, 1}
       PRr(1)=C134R(1,3)+C234R(1,2)-r10*DijR(6,3)-2*DijR(17,4)
       PRr(2)=C124R(4,3)-C134R(1,3)-r21*DijR(6,3)-4*DijR(19,4)
       PRr(3)=C123R(4,3)-C124R(4,3)-r32*DijR(6,3)
       Rr(3)=(PRr(3)-z32*(PRr(1)-z21*PRr(2))-z31*PRr(2))*iz33
       Rr(2)=(PRr(1)-z21*PRr(2)-z23*Rr(3))*iz22
       Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11
       DijR(6,4)=Rr(1)
       DijR(9,4)=Rr(2)
       DijR(10,4)=Rr(3)
       PRI(1)=C134I(1,3)+C234I(1,2)-r10*DijI(6,3)-2*DijI(17,4)
       PRI(2)=C124I(4,3)-C134I(1,3)-r21*DijI(6,3)-4*DijI(19,4)
       PRI(3)=C123I(4,3)-C124I(4,3)-r32*DijI(6,3)
       RI(3)=(PRI(3)-z32*(PRI(1)-z21*PRI(2)) -z31*PRI(2))*iz33
       RI(2)=(PRI(1)-z21*PRI(2)-z23*RI(3))*iz22
       RI(1)=(PRI(2)-  z12*RI(2)-  z13*RI(3))*iz11
       DijI(6,4)=RI(1)
       DijI(9,4)=RI(2)
       DijI(10,4)=RI(3)
c                {ii1, 2, 2, 2}
       PRr(1)=C134R(1,3)-C234R(1,3)-r10*DijR(2,3)
       PRr(2)=C124R(2,3)-C134R(1,3)-r21*DijR(2,3)-6*DijR(17,4)
       PRr(3)=C123R(2,3)-C124R(2,3)-r32*DijR(2,3)
       Rr(3)=(PRr(3)-z32*(PRr(1)-z21*PRr(2))-z31*PRr(2))*iz33
       Rr(2)=(PRr(1)-z21*PRr(2)-z23*Rr(3))*iz22
       Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11
       DijR(9,4)=Rr(1)
       DijR(2,4)=Rr(2)
       DijR(13,4)=Rr(3)
       PRI(1)=C134I(1,3)-C234I(1,3)-r10*DijI(2,3)
       PRI(2)=C124I(2,3)-C134I(1,3)-r21*DijI(2,3)-6*DijI(17,4)
       PRI(3)=C123I(2,3)-C124I(2,3)-r32*DijI(2,3)
       RI(3)=(PRI(3)-z32*(PRI(1)-z21*PRI(2)) -z31*PRI(2))*iz33
       RI(2)=(PRI(1)-z21*PRI(2)-z23*RI(3))*iz22
       RI(1)=(PRI(2)-  z12*RI(2)-  z13*RI(3))*iz11
       DijI(9,4)=RI(1)
       DijI(2,4)=RI(2)
       DijI(13,4)=RI(3)
c                {ii1, 3, 1, 1}
       PRr(1)=C134R(3,3)-C234R(2,1)-r10*DijR(5,3)-4*DijR(20,4)
       PRr(2)=C124R(3,3)-C134R(3,3)-r21*DijR(5,3)
       PRr(3)=-C124R(3,3)-r32*DijR(5,3)-2*DijR(16,4)
       Rr(3)=(PRr(3)-z32*(PRr(1)-z21*PRr(2))-z31*PRr(2))*iz33
       Rr(2)=(PRr(1)-z21*PRr(2)-z23*Rr(3))*iz22
       Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11
       DijR(5,4)=Rr(1)
       DijR(7,4)=Rr(2)
       DijR(8,4)=Rr(3)
       PRI(1)=C134I(3,3)-C234I(2,1)-r10*DijI(5,3)-4*DijI(20,4)
       PRI(2)=C124I(3,3)-C134I(3,3)-r21*DijI(5,3)
       PRI(3)=-C124I(3,3)-r32*DijI(5,3)-2*DijI(16,4)
       RI(3)=(PRI(3)-z32*(PRI(1)-z21*PRI(2)) -z31*PRI(2))*iz33
       RI(2)=(PRI(1)-z21*PRI(2)-z23*RI(3))*iz22
       RI(1)=(PRI(2)-  z12*RI(2)-  z13*RI(3))*iz11
       DijI(5,4)=RI(1)
       DijI(7,4)=RI(2)
       DijI(8,4)=RI(3)
c                {ii1, 3, 2, 1}
       PRr(1)=C134R(3,3)+C234R(3,2)-r10*DijR(10,3)-2*DijR(21,4)
       PRr(2)=C124R(4,3)-C134R(3,3)-r21*DijR(10,3)-2*DijR(20,4)
       PRr(3)=-C124R(4,3)-r32*DijR(10,3)-2*DijR(19,4)
       Rr(3)=(PRr(3)-z32*(PRr(1)-z21*PRr(2))-z31*PRr(2))*iz33
       Rr(2)=(PRr(1)-z21*PRr(2)-z23*Rr(3))*iz22
       Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11
       DijR(7,4)=Rr(1)
       DijR(10,4)=Rr(2)
       DijR(11,4)=Rr(3)
       PRI(1)=C134I(3,3)+C234I(3,2)-r10*DijI(10,3)-2*DijI(21,4)
       PRI(2)=C124I(4,3)-C134I(3,3)-r21*DijI(10,3)-2*DijI(20,4)
       PRI(3)=-C124I(4,3)-r32*DijI(10,3)-2*DijI(19,4)
       RI(3)=(PRI(3)-z32*(PRI(1)-z21*PRI(2)) -z31*PRI(2))*iz33
       RI(2)=(PRI(1)-z21*PRI(2)-z23*RI(3))*iz22
       RI(1)=(PRI(2)-  z12*RI(2)-  z13*RI(3))*iz11
       DijI(7,4)=RI(1)
       DijI(10,4)=RI(2)
       DijI(11,4)=RI(3)
c                {ii1, 3, 2, 2}
       PRr(1)=C134R(3,3)-C234R(3,3)-r10*DijR(8,3)
       PRr(2)=C124R(2,3)-C134R(3,3)-r21*DijR(8,3)-4*DijR(21,4)
       PRr(3)=-C124R(2,3)-r32*DijR(8,3)-2*DijR(17,4)
       Rr(3)=(PRr(3)-z32*(PRr(1)-z21*PRr(2))-z31*PRr(2))*iz33
       Rr(2)=(PRr(1)-z21*PRr(2)-z23*Rr(3))*iz22
       Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11
       DijR(10,4)=Rr(1)
       DijR(13,4)=Rr(2)
       DijR(14,4)=Rr(3)
       PRI(1)=C134I(3,3)-C234I(3,3)-r10*DijI(8,3)
       PRI(2)=C124I(2,3)-C134I(3,3)-r21*DijI(8,3)-4*DijI(21,4)
       PRI(3)=-C124I(2,3)-r32*DijI(8,3)-2*DijI(17,4)
       RI(3)=(PRI(3)-z32*(PRI(1)-z21*PRI(2)) -z31*PRI(2))*iz33
       RI(2)=(PRI(1)-z21*PRI(2)-z23*RI(3))*iz22
       RI(1)=(PRI(2)-  z12*RI(2)-  z13*RI(3))*iz11
       DijI(10,4)=RI(1)
       DijI(13,4)=RI(2)
       DijI(14,4)=RI(3)
c                {ii1, 3, 3, 1}
       PRr(1)=C134R(4,3)+C234R(2,2)-r10*DijR(7,3)-2*DijR(18,4)
       PRr(2)=C124R(4,3)-C134R(4,3)-r21*DijR(7,3)
       PRr(3)=-C124R(4,3)-r32*DijR(7,3)-4*DijR(20,4)
       Rr(3)=(PRr(3)-z32*(PRr(1)-z21*PRr(2))-z31*PRr(2))*iz33
       Rr(2)=(PRr(1)-z21*PRr(2)-z23*Rr(3))*iz22
       Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11
       DijR(8,4)=Rr(1)
       DijR(11,4)=Rr(2)
       DijR(12,4)=Rr(3)
       PRI(1)=C134I(4,3)+C234I(2,2)-r10*DijI(7,3)-2*DijI(18,4)
       PRI(2)=C124I(4,3)-C134I(4,3)-r21*DijI(7,3)
       PRI(3)=-C124I(4,3)-r32*DijI(7,3)-4*DijI(20,4)
       RI(3)=(PRI(3)-z32*(PRI(1)-z21*PRI(2)) -z31*PRI(2))*iz33
       RI(2)=(PRI(1)-z21*PRI(2)-z23*RI(3))*iz22
       RI(1)=(PRI(2)-  z12*RI(2)-  z13*RI(3))*iz11
       DijI(8,4)=RI(1)
       DijI(11,4)=RI(2)
       DijI(12,4)=RI(3)
c                {ii1, 3, 3, 2}
       PRr(1)=C134R(4,3)-C234R(4,3)-r10*DijR(9,3)
       PRr(2)=C124R(2,3)-C134R(4,3)-r21*DijR(9,3)-2*DijR(18,4)
       PRr(3)=-C124R(2,3)-r32*DijR(9,3)-4*DijR(21,4)
       Rr(3)=(PRr(3)-z32*(PRr(1)-z21*PRr(2))-z31*PRr(2))*iz33
       Rr(2)=(PRr(1)-z21*PRr(2)-z23*Rr(3))*iz22
       Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11
       DijR(11,4)=Rr(1)
       DijR(14,4)=Rr(2)
       DijR(15,4)=Rr(3)
       PRI(1)=C134I(4,3)-C234I(4,3)-r10*DijI(9,3)
       PRI(2)=C124I(2,3)-C134I(4,3)-r21*DijI(9,3)-2*DijI(18,4)
       PRI(3)=-C124I(2,3)-r32*DijI(9,3)-4*DijI(21,4)
       RI(3)=(PRI(3)-z32*(PRI(1)-z21*PRI(2)) -z31*PRI(2))*iz33
       RI(2)=(PRI(1)-z21*PRI(2)-z23*RI(3))*iz22
       RI(1)=(PRI(2)-  z12*RI(2)-  z13*RI(3))*iz11
       DijI(11,4)=RI(1)
       DijI(14,4)=RI(2)
       DijI(15,4)=RI(3)
c                {ii1, 3, 3, 3}
       PRr(1)=C134R(2,3)-C234R(2,3)-r10*DijR(3,3)
       PRr(2)=C124R(2,3)-C134R(2,3)-r21*DijR(3,3)
       PRr(3)=-C124R(2,3)-r32*DijR(3,3)-6*DijR(18,4)
       Rr(3)=(PRr(3)-z32*(PRr(1)-z21*PRr(2))-z31*PRr(2))*iz33
       Rr(2)=(PRr(1)-z21*PRr(2)-z23*Rr(3))*iz22
       Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11
       DijR(12,4)=Rr(1)
       DijR(15,4)=Rr(2)
       DijR(3,4)=Rr(3)
       PRI(1)=C134I(2,3)-C234I(2,3)-r10*DijI(3,3)
       PRI(2)=C124I(2,3)-C134I(2,3)-r21*DijI(3,3)
       PRI(3)=-C124I(2,3)-r32*DijI(3,3)-6*DijI(18,4)
       RI(3)=(PRI(3)-z32*(PRI(1)-z21*PRI(2)) -z31*PRI(2))*iz33
       RI(2)=(PRI(1)-z21*PRI(2)-z23*RI(3))*iz22
       RI(1)=(PRI(2)-  z12*RI(2)-  z13*RI(3))*iz11
       DijI(12,4)=RI(1)
       DijI(15,4)=RI(2)
       DijI(3,4)=RI(3)
c                {0, 0, 0, 0, ii1}
       PRr(1)=C134R(9,4)-C234R(9,4)-r10*DijR(22,4)
       PRr(2)=C124R(9,4)-C134R(9,4)-r21*DijR(22,4)
       PRr(3)=C123R(9,4)-C124R(9,4)-r32*DijR(22,4)
       Rr(3)=(PRr(3)-z32*(PRr(1)-z21*PRr(2))-z31*PRr(2))*iz33
       Rr(2)=(PRr(1)-z21*PRr(2)-z23*Rr(3))*iz22
       Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11
       DijR(32,5)=Rr(1)
       DijR(33,5)=Rr(2)
       DijR(34,5)=Rr(3)
       PRI(1)=C134I(9,4)-C234I(9,4)-r10*DijI(22,4)
       PRI(2)=C124I(9,4)-C134I(9,4)-r21*DijI(22,4)
       PRI(3)=C123I(9,4)-C124I(9,4)-r32*DijI(22,4)
       RI(3)=(PRI(3)-z32*(PRI(1)-z21*PRI(2)) -z31*PRI(2))*iz33
       RI(2)=(PRI(1)-z21*PRI(2)-z23*RI(3))*iz22
       RI(1)=(PRI(2)-  z12*RI(2)-  z13*RI(3))*iz11
       DijI(32,5)=RI(1)
       DijI(33,5)=RI(2)
       DijI(34,5)=RI(3)
c                {0, 0, ii1, 1, 1}
       PRr(1)=C134R(6,4)-C234R(4,2)-r10*DijR(16,4)-4*DijR(32,5)
       PRr(2)=C124R(6,4)-C134R(6,4)-r21*DijR(16,4)
       PRr(3)=C123R(6,4)-C124R(6,4)-r32*DijR(16,4)
       Rr(3)=(PRr(3)-z32*(PRr(1)-z21*PRr(2))-z31*PRr(2))*iz33
       Rr(2)=(PRr(1)-z21*PRr(2)-z23*Rr(3))*iz22
       Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11
       DijR(22,5)=Rr(1)
       DijR(25,5)=Rr(2)
       DijR(26,5)=Rr(3)
       PRI(1)=C134I(6,4)-C234I(4,2)-r10*DijI(16,4)-4*DijI(32,5)
       PRI(2)=C124I(6,4)-C134I(6,4)-r21*DijI(16,4)
       PRI(3)=C123I(6,4)-C124I(6,4)-r32*DijI(16,4)
       RI(3)=(PRI(3)-z32*(PRI(1)-z21*PRI(2)) -z31*PRI(2))*iz33
       RI(2)=(PRI(1)-z21*PRI(2)-z23*RI(3))*iz22
       RI(1)=(PRI(2)-  z12*RI(2)-  z13*RI(3))*iz11
       DijI(22,5)=RI(1)
       DijI(25,5)=RI(2)
       DijI(26,5)=RI(3)
c                {0, 0, ii1, 2, 1}
       PRr(1)=C134R(6,4)+C234R(5,3)-r10*DijR(19,4)-2*DijR(33,5)
       PRr(2)=C124R(8,4)-C134R(6,4)-r21*DijR(19,4)-2*DijR(32,5)
       PRr(3)=C123R(8,4)-C124R(8,4)-r32*DijR(19,4)
       Rr(3)=(PRr(3)-z32*(PRr(1)-z21*PRr(2))-z31*PRr(2))*iz33
       Rr(2)=(PRr(1)-z21*PRr(2)-z23*Rr(3))*iz22
       Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11
       DijR(25,5)=Rr(1)
       DijR(27,5)=Rr(2)
       DijR(28,5)=Rr(3)
       PRI(1)=C134I(6,4)+C234I(5,3)-r10*DijI(19,4)-2*DijI(33,5)
       PRI(2)=C124I(8,4)-C134I(6,4)-r21*DijI(19,4)-2*DijI(32,5)
       PRI(3)=C123I(8,4)-C124I(8,4)-r32*DijI(19,4)
       RI(3)=(PRI(3)-z32*(PRI(1)-z21*PRI(2)) -z31*PRI(2))*iz33
       RI(2)=(PRI(1)-z21*PRI(2)-z23*RI(3))*iz22
       RI(1)=(PRI(2)-  z12*RI(2)-  z13*RI(3))*iz11
       DijI(25,5)=RI(1)
       DijI(27,5)=RI(2)
       DijI(28,5)=RI(3)
c                {0, 0, ii1, 2, 2}
       PRr(1)=C134R(6,4)-C234R(6,4)-r10*DijR(17,4)
       PRr(2)=C124R(7,4)-C134R(6,4)-r21*DijR(17,4)-4*DijR(33,5)
       PRr(3)=C123R(7,4)-C124R(7,4)-r32*DijR(17,4)
       Rr(3)=(PRr(3)-z32*(PRr(1)-z21*PRr(2))-z31*PRr(2))*iz33
       Rr(2)=(PRr(1)-z21*PRr(2)-z23*Rr(3))*iz22
       Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11
       DijR(27,5)=Rr(1)
       DijR(23,5)=Rr(2)
       DijR(30,5)=Rr(3)
       PRI(1)=C134I(6,4)-C234I(6,4)-r10*DijI(17,4)
       PRI(2)=C124I(7,4)-C134I(6,4)-r21*DijI(17,4)-4*DijI(33,5)
       PRI(3)=C123I(7,4)-C124I(7,4)-r32*DijI(17,4)
       RI(3)=(PRI(3)-z32*(PRI(1)-z21*PRI(2)) -z31*PRI(2))*iz33
       RI(2)=(PRI(1)-z21*PRI(2)-z23*RI(3))*iz22
       RI(1)=(PRI(2)-  z12*RI(2)-  z13*RI(3))*iz11
       DijI(27,5)=RI(1)
       DijI(23,5)=RI(2)
       DijI(30,5)=RI(3)
c                {0, 0, ii1, 3, 1}
       PRr(1)=C134R(8,4)+C234R(6,3)-r10*DijR(20,4)-2*DijR(34,5)
       PRr(2)=C124R(8,4)-C134R(8,4)-r21*DijR(20,4)
       PRr(3)=-C124R(8,4)-r32*DijR(20,4)-2*DijR(32,5)
       Rr(3)=(PRr(3)-z32*(PRr(1)-z21*PRr(2))-z31*PRr(2))*iz33
       Rr(2)=(PRr(1)-z21*PRr(2)-z23*Rr(3))*iz22
       Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11
       DijR(26,5)=Rr(1)
       DijR(28,5)=Rr(2)
       DijR(29,5)=Rr(3)
       PRI(1)=C134I(8,4)+C234I(6,3)-r10*DijI(20,4)-2*DijI(34,5)
       PRI(2)=C124I(8,4)-C134I(8,4)-r21*DijI(20,4)
       PRI(3)=-C124I(8,4)-r32*DijI(20,4)-2*DijI(32,5)
       RI(3)=(PRI(3)-z32*(PRI(1)-z21*PRI(2)) -z31*PRI(2))*iz33
       RI(2)=(PRI(1)-z21*PRI(2)-z23*RI(3))*iz22
       RI(1)=(PRI(2)-  z12*RI(2)-  z13*RI(3))*iz11
       DijI(26,5)=RI(1)
       DijI(28,5)=RI(2)
       DijI(29,5)=RI(3)
c                {0, 0, ii1, 3, 2}
       PRr(1)=C134R(8,4)-C234R(8,4)-r10*DijR(21,4)
       PRr(2)=C124R(7,4)-C134R(8,4)-r21*DijR(21,4)-2*DijR(34,5)
       PRr(3)=-C124R(7,4)-r32*DijR(21,4)-2*DijR(33,5)
       Rr(3)=(PRr(3)-z32*(PRr(1)-z21*PRr(2))-z31*PRr(2))*iz33
       Rr(2)=(PRr(1)-z21*PRr(2)-z23*Rr(3))*iz22
       Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11
       DijR(28,5)=Rr(1)
       DijR(30,5)=Rr(2)
       DijR(31,5)=Rr(3)
       PRI(1)=C134I(8,4)-C234I(8,4)-r10*DijI(21,4)
       PRI(2)=C124I(7,4)-C134I(8,4)-r21*DijI(21,4)-2*DijI(34,5)
       PRI(3)=-C124I(7,4)-r32*DijI(21,4)-2*DijI(33,5)
       RI(3)=(PRI(3)-z32*(PRI(1)-z21*PRI(2)) -z31*PRI(2))*iz33
       RI(2)=(PRI(1)-z21*PRI(2)-z23*RI(3))*iz22
       RI(1)=(PRI(2)-  z12*RI(2)-  z13*RI(3))*iz11
       DijI(28,5)=RI(1)
       DijI(30,5)=RI(2)
       DijI(31,5)=RI(3)
c                {0, 0, ii1, 3, 3}
       PRr(1)=C134R(7,4)-C234R(7,4)-r10*DijR(18,4)
       PRr(2)=C124R(7,4)-C134R(7,4)-r21*DijR(18,4)
       PRr(3)=-C124R(7,4)-r32*DijR(18,4)-4*DijR(34,5)
       Rr(3)=(PRr(3)-z32*(PRr(1)-z21*PRr(2))-z31*PRr(2))*iz33
       Rr(2)=(PRr(1)-z21*PRr(2)-z23*Rr(3))*iz22
       Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11
       DijR(29,5)=Rr(1)
       DijR(31,5)=Rr(2)
       DijR(24,5)=Rr(3)
       PRI(1)=C134I(7,4)-C234I(7,4)-r10*DijI(18,4)
       PRI(2)=C124I(7,4)-C134I(7,4)-r21*DijI(18,4)
       PRI(3)=-C124I(7,4)-r32*DijI(18,4)-4*DijI(34,5)
       RI(3)=(PRI(3)-z32*(PRI(1)-z21*PRI(2)) -z31*PRI(2))*iz33
       RI(2)=(PRI(1)-z21*PRI(2)-z23*RI(3))*iz22
       RI(1)=(PRI(2)-  z12*RI(2)-  z13*RI(3))*iz11
       DijI(29,5)=RI(1)
       DijI(31,5)=RI(2)
       DijI(24,5)=RI(3)


	RETURN

        end



