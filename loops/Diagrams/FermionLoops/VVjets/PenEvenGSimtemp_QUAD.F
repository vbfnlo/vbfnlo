
       subroutine PenEvenGSim_QUAD(cpl,MM1,MM3,MM4,p1,p2,p3,p4,p5, mup1, mu
     -   p2,mup3,mup4,mup5,musq,comp,result)
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 1/02/2012
c Modified:17/7/2012
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c************************************************************************************
       IMPLICIT NONE

#include "VBFNLO/utilities/process.inc"

       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16   mup1(0:3),mup2(0:3),mup3(0:3),mup4(0:3),mup5(0:3)
       Complex*16   p1mup1, p1mup2, p1mup3, p1mup4, p1mup5, p2mup1, 
     -          p2mup2, p2mup3, p2mup4, p2mup5, p3mup1, p3mup2, 
     -          p3mup3, p3mup4, p3mup5, p4mup1, p4mup2, p4mup3, 
     -          p4mup4, p4mup5, p5mup1, p5mup2, p5mup3, p5mup4, p5mup5
       Complex*16   mup1mup2, mup1mup3, mup1mup4, mup1mup5, mup2mup3, 
     -          mup2mup4, mup2mup5, mup3mup4, mup3mup5, mup4mup5
       Real*8 dotrr
       Complex*16 B0finG,C0finG,D0finG,E0finG
       EXTERNAL dotrr,B0finG,C0finG,D0finG,E0finG
        Real*8   p1sq, p1p2, p1p3, p1p4, p1p5 
       Real*8   p2sq, p2p3, p2p4, p2p5 
       Real*8   p3sq, p3p4, p3p5 
       Real*8   p4sq, p4p5 
       Real*8   p5sq 
       Real*8   s12, s13, s14, s15 
       Real*8   s23, s24, s25 
       Real*8   s34, s35 
       Real*8   s45 
       Real*8   s234


       Real*8  B012R,B013R,B014R,B015R 
       Real*8  B023R,B024R,B025R 
       Real*8  B034R,B035R 
       Real*8  B045R 
       Real*8  B012I,B013I,B014I,B015I 
       Real*8  B023I,B024I,B025I 
       Real*8  B034I,B035I 
       Real*8  B045I 
       Real*8  B12R(3,3),B13R(3,3),B14R(3,3),B15R(3,3) 
       Real*8  B23R(3,3),B24R(3,3),B25R(3,3) 
       Real*8  B34R(3,3),B35R(3,3) 
       Real*8  B45R(3,3) 
       Real*8  B12I(3,3),B13I(3,3),B14I(3,3),B15I(3,3) 
       Real*8  B23I(3,3),B24I(3,3),B25I(3,3) 
       Real*8  B34I(3,3),B35I(3,3) 
       Real*8  B45I(3,3)


       Real*8 C0123R,C0124R,C0125R 
       Real*8 C0134R,C0135R 
       Real*8 C0145R 
       Real*8 C0234R,C0235R 
       Real*8 C0245R 
       Real*8 C0345R 
       Real*8 C0123I,C0124I,C0125I 
       Real*8 C0134I,C0135I 
       Real*8 C0145I 
       Real*8 C0234I,C0235I 
       Real*8 C0245I 
       Real*8 C0345I  
       Real*8 C123R(9,4),C124R(9,4),C125R(9,4) 
       Real*8 C134R(9,4),C135R(9,4) 
       Real*8 C145R(9,4) 
       Real*8 C234R(9,4),C235R(9,4) 
       Real*8 C245R(9,4) 
       Real*8 C345R(9,4) 
       Real*8 C123I(9,4),C124I(9,4),C125I(9,4) 
       Real*8 C134I(9,4),C135I(9,4) 
       Real*8 C145I(9,4) 
       Real*8 C234I(9,4),C235I(9,4) 
       Real*8 C245I(9,4) 
       Real*8 C345I(9,4)


       Real*8 D01234R,D01235R 
       Real*8 D01245R 
       Real*8 D01345R 
       Real*8 D02345R  
       Real*8 D01234I,D01235I 
       Real*8 D01245I 
       Real*8 D01345I 
       Real*8 D02345I 
       Real*8 D1234R(34,5),D1235R(34,5) 
       Real*8 D1245R(34,5) 
       Real*8 D1345R(34,5) 
       Real*8 D2345R(34,5) 
       Real*8 D1234I(34,5),D1235I(34,5) 
       Real*8 D1245I(34,5) 
       Real*8 D1345I(34,5) 
       Real*8 D2345I(34,5) 

       Real*8 E012345R 
       Real*8 E012345I 
       Real*8 E12345R(80,5) 
       Real*8 E12345I(80,5) 
       Logical PrintB,PrintC,PrintD,PrintE 
       Integer comp 
       Real*8  MM1,MM3,MM4,M1,M3,M4,m1sq,m3sq,m4sq 
       Real*8 cpl(10) 
      Complex*16 F(3806),Fa(123)
      Real*8 P(1320)
      Complex*16 dotrc,dotcc,result(1,6)
       Real*8 musq
      EXTERNAL   dotrc,dotcc

       COMMON/FaPenEvenGSim/Fa
       COMMON/FPenEvenGSim/F
       COMMON/PPenEvenGSim/P
      Save/FPenEvenGSim/
      Save/PPenEvenGSim/
       COMMON/InvariantsPenEvenGSim/p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s
     -   34,s45,s234,m1,m1sq,m3,m3sq,m4,m4sq
       COMMON/EVALFunctionsPenEvenGSim/ E012345R,E12345R,E012345I,E123
     -   45I
       COMMON/DVALFunctionsPenEvenGSim/ D01234R,D1234R,D01234I,D1234I,
     -    D01235R,D1235R,D01235I,D1235I, D01245R,D1245R,D01245I,D1245
     -   I, D01345R,D1345R,D01345I,D1345I, D02345R,D2345R,D02345I,D23
     -   45I
       COMMON/CVALFunctionsPenEvenGSim/ C0123R,C123R,C0123I,C123I, C01
     -   24R,C124R,C0124I,C124I, C0125R,C125R,C0125I,C125I, C0134R,C1
     -   34R,C0134I,C134I, C0135R,C135R,C0135I,C135I, C0145R,C145R,C0
     -   145I,C145I, C0234R,C234R,C0234I,C234I, C0235R,C235R,C0235I,C
     -   235I, C0245R,C245R,C0245I,C245I, C0345R,C345R,C0345I,C345I 
       COMMON/BVALFunctionsPenEvenGSim/ B012R,B012I, B013R,B013I, B014
     -   R,B014I, B015R,B015I, B023R,B023I, B024R,B024I, B025R,B025I,
     -    B034R,B034I, B035R,B035I, B045R,B045I, B12R,B12I, B13R,B13I
     -   , B14R,B14I, B15R,B15I, B23R,B23I, B24R,B24I, B25R,B25I, B34
     -   R,B34I, B35R,B35I, B45R,B45I 

ccccccccccccccccccccc For QUAD precission ccccccccccccccccccccccccc
       Real*16 musq_QUAD, dotrr_QUAD 
       External dotrr_QUAD 

       Real*16   p1sq_QUAD, p1p2_QUAD, p1p3_QUAD, p1p4_QUAD, p1p5_QUAD
       Real*16   p2sq_QUAD, p2p3_QUAD, p2p4_QUAD, p2p5_QUAD 
       Real*16   p3sq_QUAD, p3p4_QUAD, p3p5_QUAD 
       Real*16   p4sq_QUAD, p4p5_QUAD 
       Real*16   p5sq_QUAD

       Real*16   s12_QUAD, s13_QUAD, s14_QUAD, s15_QUAD
       Real*16   s23_QUAD, s24_QUAD, s25_QUAD
       Real*16   s34_QUAD, s35_QUAD
       Real*16   s45_QUAD

       Real*16  p1_QUAD(0:3),p2_QUAD(0:3),p3_QUAD(0:3),p4_QUAD(0:3),
     -          p5_QUAD(0:3)
       Real*16 ptemp(0:3,5)
       Real*8 psq(5)
       Integer i,jtemp
       Real*16   M1_QUAD,M3_QUAD
       Real*16   M4_QUAD
       integer j1,zaxis,j2,j3,j4,j5
       real*16 c,d,e03,e04

      double precision     hmass, wmass, zmass, amass,
     &                     tmass, bmass, lmass, cmass
      common /COUPL_MASS/  hmass, wmass, zmass, amass,
     &                     tmass, bmass, lmass, cmass
       integer numberZ, numberw
       double precision vmass

       M1_QUAD=MM1
       M3_QUAD=MM3
       M4_QUAD=MM4


c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c       Definition of the scalar products. Not inlcueded the contraction of the
c       moments with the external currents  
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       do i=0,3 
       ptemp(i,1)=p1(i)
       ptemp(i,2)=p2(i)
       ptemp(i,3)=p3(i)
       ptemp(i,4)=p4(i)
       ptemp(i,5)=p5(i)
       enddo

cfc       print*, "QUAD_INSIDE"
cfc             do i=0,3
cfc            write(*,'(A3,I1,A2,E25.16,E25.16)') "p1(",i,")=",ptemp(i,1)
cfc            enddo
cfc             do i=0,3
cfc            write(*,'(A3,I1,A2,E25.16,E25.16)') "p2(",i,")=",ptemp(i,2)
cfc            enddo
cfc             do i=0,3
cfc            write(*,'(A3,I1,A2,E25.16,E25.16)') "p3(",i,")=",ptemp(i,3)
cfc            enddo
cfc             do i=0,3
cfc            write(*,'(A3,I1,A2,E25.16,E25.16)') "p4(",i,")=",ptemp(i,4)
cfc             enddo
cfc             do i=0,3
cfc            write(*,'(A3,I1,A2,E25.16,E25.16)') "p5(",i,")=",ptemp(i,5)
cfc             enddo


       psq(1)=abs(dotrr(p1,p1))
       psq(2)=abs(dotrr(p2,p2))
       psq(3)=abs(dotrr(p3,p3))
       psq(4)=abs(dotrr(p4,p4))
       psq(5)=abs(dotrr(p5,p5))

       if (narrowwidth.eq.1) then
c find on-shell bosons
         numberZ=0
         numberW=0
cccccccc Identify the Z cccccccccccccccccccccccc
         j1=0
         do i=1,5
           if(j1.eq.0) then
             if(abs(psq(i)/(ZMASS*ZMASS)-1d0).lt.1d-6) then
               j1=i
               numberZ=numberZ+1
             else if(abs(psq(i)/(WMASS*WMASS)-1d0).lt.1d-6) then
               j1=i
               numberW=numberW+1
             endif
           else
             if(abs(psq(i)/(ZMASS*ZMASS)-1d0).lt.1d-6) then
               j2=i
               numberZ=numberZ+1
             else if(abs(psq(i)/(WMASS*WMASS)-1d0).lt.1d-6) then
               j2=i
               numberW=numberW+1
             endif
           endif
         enddo
         if (numberz.eq.2) vmass=zmass
         if (numberw.eq.2) vmass=wmass
       endif
       
       if ((narrowwidth.eq.1) .and. ((numberw.eq.2) .or. (NUMBERZ.EQ.2))) THEN
ccccccccccccccccccccccccccccccccccc
c first Z

c          do i=0,3
c            write(*,'(A3,I1,A2,E25.16,E25.16)') "p1(",i,")=",ptemp(i,j1)
c            enddo

       jtemp=j1
       if (ptemp(0,jtemp).ge.0d0) then
       ptemp(0,jtemp)=SQRT(VMASS*VMASS
     &                 +ptemp(1,jtemp)*ptemp(1,jtemp)+
     -                  ptemp(2,jtemp)*ptemp(2,jtemp)+
     -                  ptemp(3,jtemp)*ptemp(3,jtemp))
       else
       ptemp(0,jtemp)=-SQRT(VMASS*VMASS
     &                 +ptemp(1,jtemp)*ptemp(1,jtemp)+
     -                  ptemp(2,jtemp)*ptemp(2,jtemp)+
     -                  ptemp(3,jtemp)*ptemp(3,jtemp))
       endif
c          do i=0,3
c           write(*,'(A3,I1,A2,E25.16,E25.16)') "p1(",i,")=",ptemp(i,j1)
c           enddo
c           stop

c second Z
       jtemp=j2
c          do i=0,3
c            write(*,'(A3,I1,A2,E25.16,E25.16)') "p1(",i,")=",ptemp(i,j2)
c            enddo
       if (ptemp(0,jtemp).ge.0d0) then
       ptemp(0,jtemp)=SQRT(VMASS*VMASS
     &                 +ptemp(1,jtemp)*ptemp(1,jtemp)+
     -                  ptemp(2,jtemp)*ptemp(2,jtemp)+
     -                  ptemp(3,jtemp)*ptemp(3,jtemp))
       else
       ptemp(0,jtemp)=-SQRT(VMASS*VMASS
     &                 +ptemp(1,jtemp)*ptemp(1,jtemp)+
     -                  ptemp(2,jtemp)*ptemp(2,jtemp)+
     -                  ptemp(3,jtemp)*ptemp(3,jtemp))
       endif
c          do i=0,3
c            write(*,'(A3,I1,A2,E25.16,E25.16)') "p1(",i,")=",ptemp(i,j2)
c            enddo
c            stop
ccccccccccccccccccccccccccccccccccccccccccccccccccc
c  identify gluons in the axis and out of the axis
       j3=0
       j4=0
       zaxis=0
       do jtemp=1,5
       if(jtemp.ne.j1.and.jtemp.ne.j2) then
c          print*, "j",jtemp
       if(abs(ptemp(1,jtemp)).le.10d-7.and.
     &    abs(ptemp(2,jtemp)).le.10d-7)then
c p_0 >0 axes
       if(j3.eq.0) then
       j3=jtemp
       zaxis=zaxis+1
       else
       j4=jtemp
       zaxis=zaxis + 1
       endif
       else
       j5=jtemp
       endif
       endif
       enddo
c   Impose x, y conservation gluon out of axes with
c Z bosons

c          do i=0,3
c            write(*,'(A3,I1,A2,E25.18,E25.18)') "p1(",i,")=",ptemp(i,j5)
c            enddo
       jtemp=j5
c       print*,"j5",j5
c       print*,"x",ptemp(1,j5)+ptemp(1,j1)+ptemp(1,j2)
c       print*,"y",ptemp(2,j5)+ptemp(2,j1)+ptemp(2,j2)
       ptemp(1,jtemp)=-ptemp(1,j1)-ptemp(1,j2)
       ptemp(2,jtemp)=-ptemp(2,j1)-ptemp(2,j2)
c       print*,"x",ptemp(1,j5)+ptemp(1,j1)+ptemp(1,j2)
c       print*,"y",ptemp(2,j5)+ptemp(2,j1)+ptemp(2,j2)
       if (ptemp(0,jtemp).ge.0d0) then
       ptemp(0,jtemp)=SQRT(
     &                 +ptemp(1,jtemp)*ptemp(1,jtemp)+
     -                  ptemp(2,jtemp)*ptemp(2,jtemp)+
     -                  ptemp(3,jtemp)*ptemp(3,jtemp))
       else
       ptemp(0,jtemp)=-SQRT(
     &                 +ptemp(1,jtemp)*ptemp(1,jtemp)+
     -                  ptemp(2,jtemp)*ptemp(2,jtemp)+
     -                  ptemp(3,jtemp)*ptemp(3,jtemp))
       endif

c          do i=0,3
c            write(*,'(A3,I1,A2,E25.18,E25.18)') "p1(",i,")=",ptemp(i,j5)
c            enddo
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    Now we have to reconstuct the pgluons of the z axis
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       print*, "j1",j1
c       print*, "j2",j2
c       print*, "j3",j3
c       print*, "j4",j4
c       print*, "j5",j5

      jtemp=j3
      c= ptemp(0,j1)+ ptemp(0,j2)+ ptemp(0,j5)
      d= ptemp(3,j1)+ ptemp(3,j2)+ ptemp(3,j5)

      if (ptemp(0,jtemp).ge.0d0) then   
         If(ptemp(3,jtemp).lt.0d0) then
            if (ptemp(0,j4).ge.0d0) then   
               If(ptemp(3,j4).ge.0d0) then
c     +- ++
                   e04=(-c-d)/2q0
                   e03=(d-c)/2q0
c j3
                   ptemp(0,jtemp)=e03
                   ptemp(1,jtemp)=0q0
                   ptemp(2,jtemp)=0q0
                   ptemp(3,jtemp)=-e03
c j4
                   ptemp(0,j4)=e04
                   ptemp(1,j4)=0q0
                   ptemp(2,j4)=0q0
                   ptemp(3,j4)=e04
                   else
c   +-   +-
                    print*,"PenQUAD_queck momenta"
                   endif
            else  
c + -   - +
                  If(ptemp(3,j4).ge.0d0) then
                    print*,"PenQUAD_queck momenta"
                  else
c +-     - -        
                   e04=(c+d)/2q0
                   e03=(d-c)/2q0
                   ptemp(0,jtemp)=e03
                   ptemp(1,jtemp)=0q0
                   ptemp(2,jtemp)=0q0
                   ptemp(3,jtemp)=-e03
c j4
                   ptemp(0,j4)=-e04
                   ptemp(1,j4)=0q0
                   ptemp(2,j4)=0q0
                   ptemp(3,j4)=-e04
                 endif   
            endif
         else  
c   + +
            if (ptemp(0,j4).ge.0d0) then   
               If(ptemp(3,j4).ge.0d0) then
c     ++     ++
                print*,"PenQUAD_queck momenta"
                   else
c   ++    +-
                e04=(d-c)/2q0
                e03=(-c-d)/2q0
                   ptemp(0,jtemp)=e03
                   ptemp(1,jtemp)=0q0
                   ptemp(2,jtemp)=0q0
                   ptemp(3,jtemp)=e03
c j4
                   ptemp(0,j4)=e04
                   ptemp(1,j4)=0q0
                   ptemp(2,j4)=0q0
                   ptemp(3,j4)=-e04
                   endif
            else  
c + +   - +
                  If(ptemp(3,j4).ge.0d0) then
                e04=(-d-c)/2q0
                e03=(d-c)/2q0
                   ptemp(0,jtemp)=e03
                   ptemp(1,jtemp)=0q0
                   ptemp(2,jtemp)=0q0
                   ptemp(3,jtemp)=e03
c j4
                   ptemp(0,j4)=-e04
                   ptemp(1,j4)=0q0
                   ptemp(2,j4)=0q0
                   ptemp(3,j4)=e04
                  else
c ++     - -        
                print*,"PenQUAD_queck momenta"
                 endif   
            endif
        endif
      else
c -
         If(ptemp(3,jtemp).lt.0d0) then
            if (ptemp(0,j4).ge.0d0) then   
               If(ptemp(3,j4).ge.0d0) then
c     -- ++
                    print*,"PenQUAD_queck momenta"
                   else
c   --   +-
                   e04=(d-c)/2q0
                   e03=(c+d)/2q0
                   ptemp(0,jtemp)=-e03
                   ptemp(1,jtemp)=0q0
                   ptemp(2,jtemp)=0q0
                   ptemp(3,jtemp)=-e03
c j4
                   ptemp(0,j4)=e04
                   ptemp(1,j4)=0q0
                   ptemp(2,j4)=0q0
                   ptemp(3,j4)=-e04

                   endif
            else  
c  - -   - +
                  If(ptemp(3,j4).ge.0d0) then
                   e04=(c-d)/2q0
                   e03=(c+d)/2q0
                   ptemp(0,jtemp)=-e03
                   ptemp(1,jtemp)=0q0
                   ptemp(2,jtemp)=0q0
                   ptemp(3,jtemp)=-e03
c j4
                   ptemp(0,j4)=-e04
                   ptemp(1,j4)=0q0
                   ptemp(2,j4)=0q0
                   ptemp(3,j4)=e04
                  else
c --     - -        
                 print*,"PenQUAD_queck momenta"
                 endif   
            endif
         else  
c   - +
            if (ptemp(0,j4).ge.0d0) then   
               If(ptemp(3,j4).ge.0d0) then
c     -+     ++
                   e04=(-c-d)/2q0
                   e03=(c-d)/2q0
                   ptemp(0,jtemp)=-e03
                   ptemp(1,jtemp)=0q0
                   ptemp(2,jtemp)=0q0
                   ptemp(3,jtemp)=+e03
c j4
                   ptemp(0,j4)=e04
                   ptemp(1,j4)=0q0
                   ptemp(2,j4)=0q0
                   ptemp(3,j4)=e04
                   else
c   -+    +-
                print*,"PenQUAD_queck momenta"
                   endif
            else  
c - +   - +
                  If(ptemp(3,j4).ge.0d0) then
                print*,"PenQUAD_queck momenta"
                  else
c -+     - -        
                   e04=(c+d)/2q0
                   e03=(c-d)/2q0
                   ptemp(0,jtemp)=-e03
                   ptemp(1,jtemp)=0q0
                   ptemp(2,jtemp)=0q0
                   ptemp(3,jtemp)=+e03
c j4
                   ptemp(0,j4)=-e04
                   ptemp(1,j4)=0q0
                   ptemp(2,j4)=0q0
                   ptemp(3,j4)=-e04
                 endif   
            endif
        endif
       endif 


       else ! off-shell or no two W/Z found

c find particle with largest |q^2|
c put everything below m=10^-3 massless
         j1=0
         j2=0
         j3=0
         j4=0
         j5=0
         do jtemp=1,5
           if(psq(jtemp).le.1d-6) THEN 
             if (ptemp(0,jtemp).ge.0d0) then
               ptemp(0,jtemp)=SQRT(
     &                         +ptemp(1,jtemp)*ptemp(1,jtemp)+
     -                          ptemp(2,jtemp)*ptemp(2,jtemp)+
     -                          ptemp(3,jtemp)*ptemp(3,jtemp))
             else
               ptemp(0,jtemp)=-SQRT(
     &                         +ptemp(1,jtemp)*ptemp(1,jtemp)+
     -                          ptemp(2,jtemp)*ptemp(2,jtemp)+
     -                          ptemp(3,jtemp)*ptemp(3,jtemp))
             endif
           endif
           if (j1.eq.0) then
             j1=jtemp
           else if (j2.eq.0) then
             if (abs(psq(j1)) .gt. abs(psq(jtemp))) then
               j2 = j1
               j1 = jtemp
             else
               j2 = jtemp
             endif
           else if (j3.eq.0) then
             if (abs(psq(j2)) .gt. abs(psq(jtemp))) then
               j3 = j2
               j2 = jtemp
             else
               j3 = jtemp
             endif
           else if (j4.eq.0) then
             if (abs(psq(j3)) .gt. abs(psq(jtemp))) then
               j4 = j3
               j3 = jtemp
             else
               j4 = jtemp
             endif
           else
             if (abs(psq(j4)) .gt. abs(psq(jtemp))) then
               j5 = j4
               j4 = jtemp
             else
               j5 = jtemp
             endif
           endif
         enddo

cc off-shell 
         do i=0,3
           ptemp(i,j5)=-ptemp(i,j1)-ptemp(i,j2)-ptemp(i,j3)-ptemp(i,j4)
         enddo
       endif


          do i=0,3
            if(abs(p1(i)-ptemp(i,1)).gt.1d-8) then
            write(*,'(A3,I1,A2,E25.16,E25.16)') "p1(",i,")=",p1(i)-ptemp(i,1)
!             stop
            endif
            enddo
             do i=0,3
                if(abs(p2(i)-ptemp(i,2)).gt.1d-8)then
            write(*,'(A3,I1,A2,E25.16,E25.16)') "p2(",i,")=",p2(i)-ptemp(i,2)
!             stop
            endif
            enddo
             do i=0,3
                if(abs(p3(i)-ptemp(i,3)).gt.1d-8)then
            write(*,'(A3,I1,A2,E25.16,E25.16)') "p3(",i,")=",p3(i)-ptemp(i,3)
!             stop
            endif
            enddo
             do i=0,3
                if(abs(p4(i)-ptemp(i,4)).gt.1d-8)then
            write(*,'(A3,I1,A2,E25.16,E25.16)') "p4(",i,")=",p4(i)-ptemp(i,4)
!             stop
            endif
             enddo
             do i=0,3
                if(abs(p5(i)-ptemp(i,5)).gt.1d-8)then
            write(*,'(A3,I1,A2,E25.16,E25.16)') "p5(",i,")=",p5(i)-ptemp(i,5)
!             stop
            endif
           enddo

c       pause
cfC
cfC       jtemp=j1
cfC       if(ptemp(3,jtemp).lt.0) then
cfC       if(m(jtemp).gt.1d-7) then 
cfC       ptemp(3,jtemp)=-sqrt(ptemp(0,jtemp)**2-m(jtemp)*m(jtemp))
cfC       else
cfC       ptemp(3,jtemp)=-ptemp(0,jtemp)
cfC       endif
cfC       ptemp(0,jtemp)=ptemp(0,jtemp)
cfC       ptemp(1,jtemp)=0q0
cfC       ptemp(2,jtemp)=0q0
cfC
cfC       jtemp=j2
cfC
cfC       do jtemp=1,5
cfC       if (jtemp.ne.j2) then
cfC       ptemp(0,j2)=-ptemp(i,jtemp)+ptemp(i,j2)
cfC       endif
cfC       enddo
cfC       enddo
cfC       jtemp=j2
cfC       if(m(jtemp).gt.1d-7) then 
cfC       ptemp(3,jtemp)=sqrt(ptemp(0,jtemp)**2-m(jtemp)*m(jtemp))
cfC       else
cfC       ptemp(3,jtemp)=ptemp(0,jtemp)
cfC       endif   
cfC       ptemp(0,jtemp)=ptemp(0,jtemp)
cfC       ptemp(1,jtemp)=0q0
cfC       ptemp(2,jtemp)=0q0
cfC 
cfC       else
cfC
cfC       if(m(jtemp).gt.1d-7) then 
cfC       ptemp(3,jtemp)=sqrt(ptemp(0,jtemp)**2-m(jtemp)*m(jtemp))
cfC       else
cfC       ptemp(3,jtemp)=ptemp(0,jtemp)
cfC       endif
cfC       ptemp(0,jtemp)=ptemp(0,jtemp)
cfC       ptemp(1,jtemp)=0q0
cfC       ptemp(2,jtemp)=0q0
cfC
cfC       jtemp=j2
cfC
cfC       do jtemp=1,5
cfC       if (jtemp.ne.j2) then
cfC       ptemp(0,j2)=-ptemp(i,jtemp)+ptemp(i,j2)
cfC       endif
cfC       enddo
cfC       enddo
cfC       jtemp=j2
cfC       if(m(jtemp).gt.1d-7) then 
cfC       ptemp(3,jtemp)=-sqrt(ptemp(0,jtemp)**2-m(jtemp)*m(jtemp))
cfC       else
cfC       ptemp(3,jtemp)=-ptemp(0,jtemp)
cfC       endif   
cfC       ptemp(0,jtemp)=ptemp(0,jtemp)
cfC       ptemp(1,jtemp)=0q0
cfC       ptemp(2,jtemp)=0q0
cfC 
cfC       endif
cfC
cfC       psq(1)=abs(dotrr(p1,p1))
cfC       psq(2)=abs(dotrr(p2,p2))
cfC       psq(3)=abs(dotrr(p3,p3))
cfC       psq(4)=abs(dotrr(p4,p4))
cfC       psq(5)=abs(dotrr(p5,p5))
cfC       if (psq(1).gt.10d-7) then
cfC       j=1
cfC       elseif (psq(2).gt.10d-7)then
cfC       j=2
cfC       elseif (psq(3).gt.10d-7)then
cfC       j=3
cfC       elseif (psq(4).gt.10d-7)then
cfC       j=4
cfC       elseif (psq(5).gt.10d-7)then
cfC       j=5
cfC       else
cfC       j=7
cfCc$$$       write(*,*) 'No massive invariant'
cfCc$$$       Return
cfC       endif
cfC       if(j.ne.7) then
cfC       do jtemp=1,5
cfCc massless 
cfC       if (abs(psq(jtemp)).le.10d-7) then
cfCc x y axes
cfC       if(abs(ptemp(1,jtemp)).le.10d-7.and.
cfC     &    abs(ptemp(2,jtemp)).le.10d-7)then
cfCc p_0 >0 axes
cfC       if (ptemp(0,jtemp).ge.0d0) then
cfC       ptemp(0,jtemp)=abs(ptemp(3,jtemp))
cfC       else
cfC       ptemp(0,jtemp)=-abs(ptemp(3,jtemp))
cfC       endif ! p_0
cfCc  no x y axis
cfC       else
cfC       if (ptemp(0,jtemp).ge.0d0) then
cfC       ptemp(0,jtemp)=SQRT(ptemp(1,jtemp)*ptemp(1,jtemp)+
cfC     -                  ptemp(2,jtemp)*ptemp(2,jtemp)+
cfC     -                  ptemp(3,jtemp)*ptemp(3,jtemp))
cfC       else
cfCc$$$       print*, 'Minus jtemp', jtemp
cfC      ptemp(0,jtemp)=-SQRT(ptemp(1,jtemp)*ptemp(1,jtemp)+
cfC     -                  ptemp(2,jtemp)*ptemp(2,jtemp)+
cfC     -                  ptemp(3,jtemp)*ptemp(3,jtemp))
cfC      endif
cfC      endif
cfCc  off-shell  or massive
cfC      else
cfC      if (ptemp(0,jtemp).ge.0d0) then
cfC      ptemp(0,jtemp)=SQRT(m(jtemp)*m(jtemp)+
cfC     &                  ptemp(1,jtemp)*ptemp(1,jtemp)+
cfC     -                  ptemp(2,jtemp)*ptemp(2,jtemp)+
cfC     -                  ptemp(3,jtemp)*ptemp(3,jtemp))
cfC      else
cfCc$$$       print*, 'Minus jtemp', jtemp
cfC      ptemp(0,jtemp)=-SQRT(m(jtemp)*m(jtemp)+
cfC     &                   ptemp(1,jtemp)*ptemp(1,jtemp)+
cfC     -                  ptemp(2,jtemp)*ptemp(2,jtemp)+
cfC     -                  ptemp(3,jtemp)*ptemp(3,jtemp))
cfC      endif
cfC      endif
cfC      endif
cfC      enddo
cfC
cfCcc Impose conservation of momenta for the first off-shell particle
cfC       do i=0,3
cfC       ptemp(i,j)=0q0
cfC       do jtemp=1,5
cfC       if (jtemp.ne.j) then
cfC       ptemp(i,j)=-ptemp(i,jtemp)+ptemp(i,j)
cfC       endif
cfC       enddo
cfC       enddo
cfC
cfC       elseif(j.eq.7) then
cfC      do i=0,3 
cfC      ptemp(i,1)=p1(i)
cfC      ptemp(i,2)=p2(i)
cfC      ptemp(i,3)=p3(i)
cfC      ptemp(i,4)=p4(i)
cfC      ptemp(i,5)=p5(i)
cfC
cfC      enddo 
cfC      endif
cfC
      do i=0,3 
      p1_QUAD(i)=ptemp(i,1)
      p2_QUAD(i)=ptemp(i,2)
      p3_QUAD(i)=ptemp(i,3)
      p4_QUAD(i)=ptemp(i,4)
      p5_QUAD(i)=ptemp(i,5)

      enddo
cfc      print*, 'p1QUAD',p1_QUAD
cfc      print*, 'p2QUAD',p2_QUAD
cfc      print*, 'p3QUAD',p3_QUAD
cfc      print*, 'p4QUAD',p4_QUAD
cfc      print*, 'p5QUAD',p5_QUAD
cfcc$$$      print*, 'p6QUAD',p6_QUAD
cfc      Print*, 'Conser',(p1_QUAD(i)+p2_QUAD(i)+p3_QUAD(i)+
cfc     -   p4_QUAD(i)+p5_QUAD(i),i=0,3)
cfc      Print*, 'Conser',(p1(i)+p2(i)+p3(i)+
cfc     -   p4(i)+p5(i),i=0,3)

       p1sq_QUAD = dotrr_QUAD(p1_QUAD,p1_QUAD)
       p1p2_QUAD = dotrr_QUAD(p1_QUAD,p2_QUAD)
       p1p3_QUAD = dotrr_QUAD(p1_QUAD,p3_QUAD)
       p1p4_QUAD = dotrr_QUAD(p1_QUAD,p4_QUAD)
       p1p5_QUAD = dotrr_QUAD(p1_QUAD,p5_QUAD)

       p2sq_QUAD = dotrr_QUAD(p2_QUAD,p2_QUAD)
       p2p3_QUAD = dotrr_QUAD(p2_QUAD,p3_QUAD)
       p2p4_QUAD = dotrr_QUAD(p2_QUAD,p4_QUAD)
       p2p5_QUAD = dotrr_QUAD(p2_QUAD,p5_QUAD)

       p3sq_QUAD = dotrr_QUAD(p3_QUAD,p3_QUAD)
       p3p4_QUAD = dotrr_QUAD(p3_QUAD,p4_QUAD)
       p3p5_QUAD = dotrr_QUAD(p3_QUAD,p5_QUAD)

       p4sq_QUAD = dotrr_QUAD(p4_QUAD,p4_QUAD)
       p4p5_QUAD = dotrr_QUAD(p4_QUAD,p5_QUAD)

       p5sq_QUAD = dotrr_QUAD(p5_QUAD,p5_QUAD)


       s12_QUAD = (p1sq_QUAD +p2sq_QUAD+ 2*p1p2_QUAD) 
       s13_QUAD = (p1sq_QUAD +p3sq_QUAD+ 2*p1p3_QUAD) 
       s14_QUAD = (p1sq_QUAD +p4sq_QUAD+ 2*p1p4_QUAD) 
       s15_QUAD = (p1sq_QUAD +p5sq_QUAD+ 2*p1p5_QUAD) 

       s23_QUAD = (p2sq_QUAD +p3sq_QUAD+ 2*p2p3_QUAD) 
       s24_QUAD = (p2sq_QUAD +p4sq_QUAD+ 2*p2p4_QUAD) 
       s25_QUAD = (p2sq_QUAD +p5sq_QUAD+ 2*p2p5_QUAD) 

       s34_QUAD = (p3sq_QUAD +p4sq_QUAD+ 2*p3p4_QUAD) 
       s35_QUAD = (p3sq_QUAD +p5sq_QUAD+ 2*p3p5_QUAD) 

       s45_QUAD = (p4sq_QUAD +p5sq_QUAD+ 2*p4p5_QUAD) 

c       Write(*,'(a5,E44.36)')," p1sq_QUAD ", p1sq_QUAD 
c       Write(*,'(a5,E44.36)')," p1p2_QUAD ", p1p2_QUAD
c       Write(*,'(a5,E44.36)')," p1p3_QUAD ", p1p3_QUAD
c       Write(*,'(a5,E44.36)')," p1p4_QUAD ", p1p4_QUAD
c       Write(*,'(a5,E44.36)')," p1p5_QUAD ", p1p5_QUAD
c       Write(*,'(a5,E44.36)')," p1p6_QUAD ", p1p6_QUAD
c       Write(*,'(a5,E44.36)')," p2sq_QUAD ", p2sq_QUAD 
c       Write(*,'(a5,E44.36)')," p2p3_QUAD ", p2p3_QUAD
c       Write(*,'(a5,E44.36)')," p2p4_QUAD ", p2p4_QUAD
c       Write(*,'(a5,E44.36)')," p2p5_QUAD ", p2p5_QUAD
c       Write(*,'(a5,E44.36)')," p2p6_QUAD ", p2p6_QUAD
c       Write(*,'(a5,E44.36)')," p3sq_QUAD ", p3sq_QUAD 
c       Write(*,'(a5,E44.36)')," p3p4_QUAD ", p3p4_QUAD
c       Write(*,'(a5,E44.36)')," p3p5_QUAD ", p3p5_QUAD
c       Write(*,'(a5,E44.36)')," p3p6_QUAD ", p3p6_QUAD
c       Write(*,'(a5,E44.36)')," p4sq_QUAD ", p4sq_QUAD 
c       Write(*,'(a5,E44.36)')," p4p5_QUAD ", p4p5_QUAD
c       Write(*,'(a5,E44.36)')," p4p6_QUAD ", p4p6_QUAD
c       Write(*,'(a5,E44.36)')," p5sq_QUAD ", p5sq_QUAD 
c       Write(*,'(a5,E44.36)')," p5p6_QUAD ", p5p6_QUAD
c       Write(*,'(a5,E44.36)')," p6sq_QUAD ", p6sq_QUAD 

       p1sq=p1sq_QUAD 
       p1p2=p1p2_QUAD
       p1p3=p1p3_QUAD
       p1p4=p1p4_QUAD
       p1p5=p1p5_QUAD

       p2sq=p2sq_QUAD 
       p2p3=p2p3_QUAD
       p2p4=p2p4_QUAD
       p2p5=p2p5_QUAD

       p3sq=p3sq_QUAD 
       p3p4=p3p4_QUAD
       p3p5=p3p5_QUAD

       p4sq=p4sq_QUAD 
       p4p5=p4p5_QUAD

       p5sq=p5sq_QUAD 
 
       s12=s12_QUAD 
       s13=s13_QUAD 
       s14=s14_QUAD 
       s15=s15_QUAD 

       s23=s23_QUAD 
       s24=s24_QUAD 
       s25=s25_QUAD 

       s34=s34_QUAD 
       s35=s35_QUAD 

       s45=s45_QUAD


        s234=s15
        m1=mm1
        m3=mm3
        m4=mm4
        m1sq=m1*m1
        m3sq=m3*m3
        m4sq=m4*m4


       do i=0,3
         p1(i)=p1_QUAD(i)
         p2(i)=p2_QUAD(i)
         p3(i)=p3_QUAD(i)
         p4(i)=p4_QUAD(i)
         p5(i)=p5_QUAD(i)

         enddo

          If(abs(musq-s12).lt.10d-7) then
          musq_QUAD=s12_QUAD
          elseIf(abs(musq-s23).lt.10d-7) then
          musq_QUAD=s23_QUAD
          elseIf(abs(musq-s34).lt.10d-7) then
          musq_QUAD=s34_QUAD
          elseIf(abs(musq-s45).lt.10d-7) then
          musq_QUAD=s45_QUAD
          elseIf(abs(musq-s15).lt.10d-7) then
          musq_QUAD=s15_QUAD
          else
          musq_QUAD=musq
          endif


cFC       p1sq = dotrr(p1,p1)
cFC       p1p2 = dotrr(p1,p2)
cFC       p1p3 = dotrr(p1,p3)
cFC       p1p4 = dotrr(p1,p4)
cFC       p1p5 = dotrr(p1,p5)
cFC       p2sq = dotrr(p2,p2)
cFC       p2p3 = dotrr(p2,p3)
cFC       p2p4 = dotrr(p2,p4)
cFC       p2p5 = dotrr(p2,p5)
cFC       p3sq = dotrr(p3,p3)
cFC       p3p4 = dotrr(p3,p4)
cFC       p3p5 = dotrr(p3,p5)
cFC       p4sq = dotrr(p4,p4)
cFC       p4p5 = dotrr(p4,p5)
cFC       p5sq = dotrr(p5,p5)
cFC       s12 = (p1sq +p2sq+ 2*p1p2) 
cFC       s13 = (p1sq +p3sq+ 2*p1p3) 
cFC       s14 = (p1sq +p4sq+ 2*p1p4) 
cFC       s15 = (p1sq +p5sq+ 2*p1p5) 
cFC       s23 = (p2sq +p3sq+ 2*p2p3) 
cFC       s24 = (p2sq +p4sq+ 2*p2p4) 
cFC       s25 = (p2sq +p5sq+ 2*p2p5) 
cFC       s34 = (p3sq +p4sq+ 2*p3p4) 
cFC       s35 = (p3sq +p5sq+ 2*p3p5) 
cFC       s45 = (p4sq +p5sq+ 2*p4p5) 
c       Write(*,'(a5,F20.10)')," p1sq ", p1sq 
c       Write(*,'(a5,F20.10)')," p1p2 ", p1p2
c       Write(*,'(a5,F20.10)')," p1p3 ", p1p3
c       Write(*,'(a5,F20.10)')," p1p4 ", p1p4
c       Write(*,'(a5,F20.10)')," p1p5 ", p1p5
c       Write(*,'(a5,F20.10)')," p2sq ", p2sq 
c       Write(*,'(a5,F20.10)')," p2p3 ", p2p3
c       Write(*,'(a5,F20.10)')," p2p4 ", p2p4
c       Write(*,'(a5,F20.10)')," p2p5 ", p2p5
c       Write(*,'(a5,F20.10)')," p3sq ", p3sq 
c       Write(*,'(a5,F20.10)')," p3p4 ", p3p4
c       Write(*,'(a5,F20.10)')," p3p5 ", p3p5
c       Write(*,'(a5,F20.10)')," p4sq ", p4sq 
c       Write(*,'(a5,F20.10)')," p4p5 ", p4p5
c       Write(*,'(a5,F20.10)')," p5sq ", p5sq


      PrintB=.False. 
      PrintC=.False. 
      PrintD=.False. 
      PrintE=.False.
c************************************************************************************
c************************************************************************************
        if (comp.gt.0) then    
c    Calling C_ij,D_ij,E_ij,F_ij Functions    
c************************************************************************************
c************************************************************************************
      call  XEgetSim_R5_QUAD(M1_QUAD,M1_QUAD,M1_QUAD,M3_QUAD,M4_QUAD,
     &   p1sq_QUAD,p2sq_QUAD,p3sq_QUAD,p4sq_QUAD,p5sq_QUAD,
     &   p1p2_QUAD,p1p3_QUAD,p1p4_QUAD,p2p3_QUAD,p2p4_QUAD,p3p4_QUAD,
     &   p1p5_QUAD,p4p5_QUAD,musq_QUAD)

c      print*, "here_XGETEM1M3 end"
c************************************************************************************
c************************************************************************************
c************************************************************************************
c       Definition of the F,P functions:Independent of the currents    
c************************************************************************************
c************************************************************************************
c************************************************************************************
       call FPenEvenGSim1(F(1))
       call FPenEvenGSim2(F(191))
       call FPenEvenGSim3(F(381))
       call FPenEvenGSim4(F(571))
       call FPenEvenGSim5(F(761))
       call FPenEvenGSim6(F(951))
       call FPenEvenGSim7(F(1141))
       call FPenEvenGSim8(F(1331))
       call FPenEvenGSim9(F(1521))
       call FPenEvenGSim10(F(1711))
       call FPenEvenGSim11(F(1901))
       call FPenEvenGSim12(F(2091))
       call FPenEvenGSim13(F(2281))
       call FPenEvenGSim14(F(2471))
       call FPenEvenGSim15(F(2661))
       call FPenEvenGSim16(F(2851))
       call FPenEvenGSim17(F(3041))
       call FPenEvenGSim18(F(3231))
       call FPenEvenGSim19(F(3421))
       call FPenEvenGSim20(F(3611))
cfc       do j1=1,1320
cfc          print*,"P",j1,P(j1)
cfc       enddo
cfc       do j1=1,3806
cfc          print*,"F",j1,Dble(F(j1))
cfc          print*,"F",j1,DImag(F(j1))
cfc       enddo

c************************************************************************************
c************************************************************************************
c************************************************************************************
       endif  
c               PART THAT DEPENDS ON THE EXTERNAL CURRENT
c************************************************************************************
c************************************************************************************
c************************************************************************************
       p1mup1 = dotrc(p1,mup1)
       p1mup2 = dotrc(p1,mup2)
       p1mup3 = dotrc(p1,mup3)
       p1mup4 = dotrc(p1,mup4)
       p1mup5 = dotrc(p1,mup5)
       p2mup1 = dotrc(p2,mup1)
       p2mup2 = dotrc(p2,mup2)
       p2mup3 = dotrc(p2,mup3)
       p2mup4 = dotrc(p2,mup4)
       p2mup5 = dotrc(p2,mup5)
       p3mup1 = dotrc(p3,mup1)
       p3mup2 = dotrc(p3,mup2)
       p3mup3 = dotrc(p3,mup3)
       p3mup4 = dotrc(p3,mup4)
       p3mup5 = dotrc(p3,mup5)
       p4mup1 = dotrc(p4,mup1)
       p4mup2 = dotrc(p4,mup2)
       p4mup3 = dotrc(p4,mup3)
       p4mup4 = dotrc(p4,mup4)
       p4mup5 = dotrc(p4,mup5)
       p5mup1 = dotrc(p5,mup1)
       p5mup2 = dotrc(p5,mup2)
       p5mup3 = dotrc(p5,mup3)
       p5mup4 = dotrc(p5,mup4)
       p5mup5 = dotrc(p5,mup5)
       mup1mup2 = dotcc(mup1,mup2)
       mup1mup3 = dotcc(mup1,mup3)
       mup1mup4 = dotcc(mup1,mup4)
       mup1mup5 = dotcc(mup1,mup5)
       mup2mup3 = dotcc(mup2,mup3)
       mup2mup4 = dotcc(mup2,mup4)
       mup2mup5 = dotcc(mup2,mup5)
       mup3mup4 = dotcc(mup3,mup4)
       mup3mup5 = dotcc(mup3,mup5)
       mup4mup5 = dotcc(mup4,mup5)


cfc       print*, "p1mup1", p1mup1
cfc       print*, "p1mup2", p1mup2
cfc       print*, "p1mup3", p1mup3
cfc       print*, "p1mup4", p1mup4
cfc       print*, "p1mup5", p1mup5
cfc       print*, "p2mup1", p2mup1
cfc       print*, "p2mup2", p2mup2
cfc       print*, "p2mup3", p2mup3
cfc       print*, "p2mup4", p2mup4
cfc       print*, "p2mup5", p2mup5
cfc       print*, "p3mup1", p3mup1
cfc       print*, "p3mup2", p3mup2
cfc       print*, "p3mup3", p3mup3
cfc       print*, "p3mup4", p3mup4
cfc       print*, "p3mup5", p3mup5
cfc       print*, "p4mup1", p4mup1
cfc       print*, "p4mup2", p4mup2
cfc       print*, "p4mup3", p4mup3
cfc       print*, "p4mup4", p4mup4
cfc       print*, "p4mup5", p4mup5
cfc       print*, "p5mup1", p5mup1
cfc       print*, "p5mup2", p5mup2
cfc       print*, "p5mup3", p5mup3
cfc       print*, "p5mup4", p5mup4
cfc       print*, "p5mup5", p5mup5
cfc       print*, "mup1mup2",mup1mup2
cfc       print*, "mup1mup3",mup1mup3
cfc       print*, "mup1mup4",mup1mup4
cfc       print*, "mup1mup5",mup1mup5
cfc       print*, "mup2mup3",mup2mup3
cfc       print*, "mup2mup4",mup2mup4
cfc       print*, "mup2mup5",mup2mup5
cfc       print*, "mup3mup4",mup3mup4
cfc       print*, "mup3mup5",mup3mup5
cfc       print*, "mup4mup5",mup4mup5


c************** Calling the Fa functions**********************************************************************
c************************************************************************************
c************************************************************************************
       call FaPenEvenGSim1(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p
     -   2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup
     -   5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5
     -   mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup
     -   2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(1))
       call FaPenEvenGSim2(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p
     -   2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup
     -   5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5
     -   mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup
     -   2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(9))
       call FaPenEvenGSim3(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p
     -   2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup
     -   5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5
     -   mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup
     -   2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(17))
       call FaPenEvenGSim4(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p
     -   2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup
     -   5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5
     -   mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup
     -   2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(25))
       call FaPenEvenGSim5(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p
     -   2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup
     -   5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5
     -   mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup
     -   2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(33))
       call FaPenEvenGSim6(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p
     -   2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup
     -   5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5
     -   mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup
     -   2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(41))
       if(comp.gt.0) then
       call FaPenEvenGSim7(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p
     -   2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup
     -   5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5
     -   mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup
     -   2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(49))
       call FaPenEvenGSim8(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p
     -   2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup
     -   5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5
     -   mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup
     -   2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(57))
       call FaPenEvenGSim9(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p
     -   2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup
     -   5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5
     -   mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup
     -   2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(65))
       call FaPenEvenGSim10(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,
     -   p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mu
     -   p5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p
     -   5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mu
     -   p2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(73))
       call FaPenEvenGSim11(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,
     -   p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mu
     -   p5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p
     -   5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mu
     -   p2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(81))
       call FaPenEvenGSim12(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,
     -   p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mu
     -   p5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p
     -   5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mu
     -   p2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(89))
        call FaPenEvenGSim13(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,
     -   p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mu
     -   p5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p
     -   5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mu
     -   p2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(97))
         endif
       call FaPenEvenGSim14(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,
     -   p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mu
     -   p5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p
     -   5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mu
     -   p2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(104))
       if(comp.gt.0) then
       call FaPenEvenGSim15(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,
     -   p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mu
     -   p5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p
     -   5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mu
     -   p2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(113))
       endif


cfc      do j1=1,123
cfc          print*,"Fa",j1,Fa(j1)
cfc       enddo

c************************************************************************************
c************************************************************************************
c************************************************************************************
c       Definition of the Matrix Element  
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c       Amplitude                         
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************


c The Finite virtual piece should be multiplied for (-1)  since 
c I have multiplied by (-I) to get the F's and k's without (I) factor
c . The factorization from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps Gamma[1+Eps] (musq)^(-Eps))
c  c So, I*I=(-1)!!!
       result(1,1) = Cpl(8)*(Cpl(9)*(Cpl(10)*Fa(104)+M1*M4*Cpl(5)*Fa(1
     -   05))+M3*Cpl(4)*(M4*Cpl(10)*Fa(106)+M1*Cpl(5)*Fa(107)))+Cpl(3
     -   )*(Cpl(4)*(Cpl(5)*Fa(104)+M1*M4*Cpl(10)*Fa(105))+M3*Cpl(9)*(
     -   M4*Cpl(5)*Fa(106)+M1*Cpl(10)*Fa(107)))
       if(comp.gt.0)then
       result(1,2) = Cpl(8)*(Cpl(9)*(Cpl(10)*Fa(108)+M1*M4*Cpl(5)*Fa(1
     -   09))+M3*Cpl(4)*(M4*Cpl(10)*Fa(110)+M1*Cpl(5)*Fa(111)))+Cpl(3
     -   )*(Cpl(4)*(Cpl(5)*Fa(108)+M1*M4*Cpl(10)*Fa(109))+M3*Cpl(9)*(
     -   M4*Cpl(5)*Fa(110)+M1*Cpl(10)*Fa(111)))
       result(1,3) = Cpl(8)*(Cpl(9)*(Cpl(10)*Fa(112)+M1*M4*Cpl(5)*Fa(1
     -   13))+M3*Cpl(4)*(M4*Cpl(10)*Fa(114)+M1*Cpl(5)*Fa(115)))+Cpl(3
     -   )*(Cpl(4)*(Cpl(5)*Fa(112)+M1*M4*Cpl(10)*Fa(113))+M3*Cpl(9)*(
     -   M4*Cpl(5)*Fa(114)+M1*Cpl(10)*Fa(115)))
       result(1,4) = (Cpl(8)*(Cpl(9)*(2*Cpl(10)*Fa(116)+M1*M4*Cpl(5)*F
     -   a(117))+M3*Cpl(4)*(M4*Cpl(10)*Fa(118)+M1*Cpl(5)*Fa(119)))+Cp
     -   l(3)*(Cpl(4)*(2*Cpl(5)*Fa(116)+M1*M4*Cpl(10)*Fa(117))+M3*Cpl
     -   (9)*(M4*Cpl(5)*Fa(118)+M1*Cpl(10)*Fa(119))))/2.
       result(1,5) = (Cpl(8)*(Cpl(9)*(2*Cpl(10)*Fa(120)+M1*M4*Cpl(5)*F
     -   a(121))+M3*Cpl(4)*(M4*Cpl(10)*Fa(122)+M1*Cpl(5)*Fa(123)))+Cp
     -   l(3)*(Cpl(4)*(2*Cpl(5)*Fa(120)+M1*M4*Cpl(10)*Fa(121))+M3*Cpl
     -   (9)*(M4*Cpl(5)*Fa(122)+M1*Cpl(10)*Fa(123))))/2.
       endif
       result(1,6) = 0
       result(1,1) =-result(1,1)
       result(1,2) =-result(1,2)
       result(1,3) =-result(1,3)
       result(1,4) =-result(1,4)
       result(1,5) =-result(1,5)
       result(1,6) =-result(1,6)
c************************************************************************************
c************************************************************************************
       Return
       End


cMK       function dotrr_QUAD(p1,p2)
cMK       implicit none      
cMK       real*16 p1(0:3),p2(0:3),dotrr_QUAD
cMK       dotrr_QUAD = p1(0)*p2(0) - p1(1)*p2(1) - p1(2)*p2(2) - p1(3)*p2(3)
cMK       end
cMK      
cMK      function dotrc_QUAD(p1,p2)
cMK      implicit none      
cMK      real*16 p1(0:3)
cMK      complex*32  p2(0:3),dotrc_QUAD
cMK      dotrc_QUAD = p1(0)*p2(0) - p1(1)*p2(1) - p1(2)*p2(2) - p1(3)*p2(3)
cMK      end
cMK
cMK      function dotcc_QUAD(p1,p2)
cMK      implicit none      
cMK      complex*32 dotcc_QUAD,p1(0:3),p2(0:3)
cMK      dotcc_QUAD = p1(0)*p2(0) - p1(1)*p2(1) - p1(2)*p2(2) - 
cMK     #     p1(3)*p2(3)
cMK      end


