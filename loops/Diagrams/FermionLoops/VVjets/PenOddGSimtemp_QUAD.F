       subroutine PenOddGSim_QUAD(cpl,MM1,MM3,MM4,p1,p2,p3,p4,p5, mup1, mup
     -   2,mup3,mup4,mup5,musq,comp,result)
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 1/02/2012
c Modified:16/7/2012
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
      Complex*16 F(4476),Fa(3196)
      Real*8 P(1383)
      Complex*16 dotrc,dotcc,result(1,6)
       Real*8 musq
      EXTERNAL   dotrc,dotcc

      Complex*16 EpsF(126)
       common/EPSFPenOddGSim/EPSF
      Complex*16 EPSCCCC,EPSCCCR,EPSCCRR,EPSCRRR
       Real*8 EPSRRRR
      EXTERNAL   EPSCCCC,EPSCCCR,EPSCCRR,EPSCRRR,EPSRRRR
       COMMON/FaPenOddGSim/Fa
       COMMON/FPenOddGSim/F
       COMMON/PPenOddGSim/P
      Save/FPenOddGSim/
      Save/PPenOddGSim/
       COMMON/InvariantsPenOddGSim/p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s3
     -   4,s45,s234,m1,m1sq,m3,m3sq,m4,m4sq
       COMMON/EVALFunctionsPenOddGSim/ E012345R,E12345R,E012345I,E1234
     -   5I
       COMMON/DVALFunctionsPenOddGSim/ D01234R,D1234R,D01234I,D1234I, 
     -   D01235R,D1235R,D01235I,D1235I, D01245R,D1245R,D01245I,D1245I
     -   , D01345R,D1345R,D01345I,D1345I, D02345R,D2345R,D02345I,D234
     -   5I
       COMMON/CVALFunctionsPenOddGSim/ C0123R,C123R,C0123I,C123I, C012
     -   4R,C124R,C0124I,C124I, C0125R,C125R,C0125I,C125I, C0134R,C13
     -   4R,C0134I,C134I, C0135R,C135R,C0135I,C135I, C0145R,C145R,C01
     -   45I,C145I, C0234R,C234R,C0234I,C234I, C0235R,C235R,C0235I,C2
     -   35I, C0245R,C245R,C0245I,C245I, C0345R,C345R,C0345I,C345I 
       COMMON/BVALFunctionsPenOddGSim/ B012R,B012I, B013R,B013I, B014R
     -   ,B014I, B015R,B015I, B023R,B023I, B024R,B024I, B025R,B025I, 
     -   B034R,B034I, B035R,B035I, B045R,B045I, B12R,B12I, B13R,B13I,
     -    B14R,B14I, B15R,B15I, B23R,B23I, B24R,B24I, B25R,B25I, B34R
     -   ,B34I, B35R,B35I, B45R,B45I 

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
c$$$      print*, 'p1QUAD',p1_QUAD
c$$$      print*, 'p2QUAD',p2_QUAD
c$$$      print*, 'p3QUAD',p3_QUAD
c$$$      print*, 'p4QUAD',p4_QUAD
c$$$      print*, 'p5QUAD',p5_QUAD
c$$$c$$$      print*, 'p6QUAD',p6_QUAD
c$$$      Print*, 'Conser',(p1_QUAD(i)+p2_QUAD(i)+p3_QUAD(i)+
c$$$     -   p4_QUAD(i)+p5_QUAD(i),i=0,3)
c$$$      Print*, 'Conser',(p1(i)+p2(i)+p3(i)+
c$$$     -   p4(i)+p5(i)+p6(i),i=0,3)

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
c      s234=s15
c       m1=mm1
c        m3=mm3
c        m4=mm4
c       m1sq=m1*m1
c      m3sq=m3*m3
c      m4sq=m4*m4

      PrintB=.False. 
      PrintC=.False. 
      PrintD=.False. 
      PrintE=.False.


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

cfc
cfc
cfc       p1sq = dotrr(p1,p1)
cfc       p1p2 = dotrr(p1,p2)
cfc       p1p3 = dotrr(p1,p3)
cfc       p1p4 = dotrr(p1,p4)
cfc       p1p5 = dotrr(p1,p5)
cfc       p2sq = dotrr(p2,p2)
cfc       p2p3 = dotrr(p2,p3)
cfc       p2p4 = dotrr(p2,p4)
cfc       p2p5 = dotrr(p2,p5)
cfc       p3sq = dotrr(p3,p3)
cfc       p3p4 = dotrr(p3,p4)
cfc       p3p5 = dotrr(p3,p5)
cfc       p4sq = dotrr(p4,p4)
cfc       p4p5 = dotrr(p4,p5)
cfc       p5sq = dotrr(p5,p5)
cfc       s12 = (p1sq +p2sq+ 2*p1p2) 
cfc       s13 = (p1sq +p3sq+ 2*p1p3) 
cfc       s14 = (p1sq +p4sq+ 2*p1p4) 
cfc       s15 = (p1sq +p5sq+ 2*p1p5) 
cfc       s23 = (p2sq +p3sq+ 2*p2p3) 
cfc       s24 = (p2sq +p4sq+ 2*p2p4) 
cfc       s25 = (p2sq +p5sq+ 2*p2p5) 
cfc       s34 = (p3sq +p4sq+ 2*p3p4) 
cfc       s35 = (p3sq +p5sq+ 2*p3p5) 
cfc       s45 = (p4sq +p5sq+ 2*p4p5) 
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

c************************************************************************************
c************************************************************************************
        if (comp.gt.0) then    
c    Calling C_ij,D_ij,E_ij,F_ij Functions    
c************************************************************************************
c************************************************************************************
       call  XEgetPenOddGSim_R5_QUAD(M1_QUAD,M1_QUAD,M1_QUAD,M3_QUAD,M4_QUAD,
     &   p1sq_QUAD,p2sq_QUAD,p3sq_QUAD,p4sq_QUAD,p5sq_QUAD,
     &   p1p2_QUAD,p1p3_QUAD,p1p4_QUAD,p2p3_QUAD,p2p4_QUAD,p3p4_QUAD,
     &   p1p5_QUAD,p4p5_QUAD,musq_QUAD)
c************************************************************************************
c************************************************************************************
c************************************************************************************
c       Definition of the F,P functions:Independent of the currents    
c************************************************************************************
c************************************************************************************
c************************************************************************************
       call FPenOddGSim1(F(1))
       call FPenOddGSim2(F(224))
       call FPenOddGSim3(F(447))
       call FPenOddGSim4(F(670))
       call FPenOddGSim5(F(893))
       call FPenOddGSim6(F(1116))
       call FPenOddGSim7(F(1339))
       call FPenOddGSim8(F(1562))
       call FPenOddGSim9(F(1785))
       call FPenOddGSim10(F(2008))
       call FPenOddGSim11(F(2231))
       call FPenOddGSim12(F(2454))
       call FPenOddGSim13(F(2677))
       call FPenOddGSim14(F(2900))
       call FPenOddGSim15(F(3123))
       call FPenOddGSim16(F(3346))
       call FPenOddGSim17(F(3569))
       call FPenOddGSim18(F(3792))
       call FPenOddGSim19(F(4015))
       call FPenOddGSim20(F(4238))
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
c************** Calling the Fa functions**********************************************************************
       EpsF(1)=EPSCCRR(mup3,mup4,p1,p2)
       EpsF(2)=EPSCCRR(mup3,mup4,p2,p3)
       EpsF(3)=EPSCCRR(mup3,mup4,p2,p4)
       EpsF(4)=EPSCCRR(mup2,mup4,p1,p2)
       EpsF(5)=EPSCCRR(mup2,mup4,p2,p3)
       EpsF(6)=EPSCCRR(mup2,mup4,p2,p4)
       EpsF(7)=EPSCCRR(mup2,mup3,p1,p2)
       EpsF(8)=EPSCCRR(mup2,mup3,p2,p3)
       EpsF(9)=EPSCCRR(mup2,mup3,p2,p4)
       EpsF(10)=EPSCCCR(mup2,mup3,mup4,p1)
       EpsF(11)=EPSCCCR(mup2,mup3,mup4,p2)
       EpsF(12)=EPSCCCR(mup2,mup3,mup4,p3)
       EpsF(13)=EPSCCCR(mup2,mup3,mup4,p4)
       EpsF(14)=EPSCCRR(mup1,mup4,p1,p2)
       EpsF(15)=EPSCCRR(mup1,mup4,p2,p3)
       EpsF(16)=EPSCCRR(mup1,mup4,p2,p4)
       EpsF(17)=EPSCCRR(mup1,mup3,p1,p2)
       EpsF(18)=EPSCCRR(mup1,mup3,p2,p3)
       EpsF(19)=EPSCCRR(mup1,mup3,p2,p4)
       EpsF(20)=EPSCCCC(mup1,mup2,mup3,mup4)
       EpsF(21)=EPSCCCR(mup1,mup3,mup4,p1)
       EpsF(22)=EPSCCCR(mup1,mup3,mup4,p2)
       EpsF(23)=EPSCCCR(mup1,mup3,mup4,p3)
       EpsF(24)=EPSCCCR(mup1,mup3,mup4,p4)
       EpsF(25)=EPSCCRR(mup1,mup2,p1,p2)
       EpsF(26)=EPSCCRR(mup1,mup2,p2,p3)
       EpsF(27)=EPSCCRR(mup1,mup2,p2,p4)
       EpsF(28)=EPSCCCR(mup1,mup2,mup4,p1)
       EpsF(29)=EPSCCCR(mup1,mup2,mup4,p2)
       EpsF(30)=EPSCCCR(mup1,mup2,mup4,p3)
       EpsF(31)=EPSCCCR(mup1,mup2,mup4,p4)
       EpsF(32)=EPSCCCR(mup1,mup2,mup3,p1)
       EpsF(33)=EPSCCCR(mup1,mup2,mup3,p2)
       EpsF(34)=EPSCCCR(mup1,mup2,mup3,p3)
       EpsF(35)=EPSCCCR(mup1,mup2,mup3,p4)
       EpsF(36)=EPSCCCC(mup1,mup2,mup3,mup5)
       EpsF(37)=EPSCCCR(mup1,mup2,mup5,p2)
       EpsF(38)=EPSCCCR(mup1,mup2,mup5,p3)
       EpsF(39)=EPSCCRR(mup1,mup2,p1,p3)
       EpsF(40)=EPSCCRR(mup1,mup2,p3,p4)
       EpsF(41)=EPSCCCC(mup1,mup3,mup4,mup5)
       EpsF(42)=EPSCCCR(mup1,mup3,mup5,p3)
       EpsF(43)=EPSCCRR(mup1,mup3,p1,p3)
       EpsF(44)=EPSCCRR(mup1,mup3,p3,p4)
       EpsF(45)=EPSCCCR(mup1,mup3,mup5,p2)
       EpsF(46)=EPSCCRR(mup1,mup3,p1,p4)
       EpsF(47)=EPSCCCC(mup1,mup2,mup4,mup5)
       EpsF(48)=EPSCCCR(mup1,mup4,mup5,p1)
       EpsF(49)=EPSCCCR(mup1,mup4,mup5,p2)
       EpsF(50)=EPSCCCR(mup1,mup4,mup5,p3)
       EpsF(51)=EPSCCRR(mup1,mup4,p1,p3)
       EpsF(52)=EPSCCCR(mup1,mup4,mup5,p4)
       EpsF(53)=EPSCCRR(mup1,mup4,p1,p4)
       EpsF(54)=EPSCCRR(mup1,mup4,p3,p4)
       EpsF(55)=EPSCCCR(mup1,mup2,mup5,p1)
       EpsF(56)=EPSCCCR(mup1,mup2,mup5,p4)
       EpsF(57)=EPSCCCR(mup1,mup3,mup5,p1)
       EpsF(58)=EPSCCCR(mup1,mup3,mup5,p4)
       EpsF(59)=EPSCCRR(mup1,mup5,p1,p2)
       EpsF(60)=EPSCCRR(mup1,mup5,p1,p3)
       EpsF(61)=EPSCCRR(mup1,mup5,p1,p4)
       EpsF(62)=EPSCCRR(mup1,mup5,p2,p3)
       EpsF(63)=EPSCCRR(mup1,mup5,p2,p4)
       EpsF(64)=EPSCCRR(mup1,mup5,p3,p4)
       EpsF(65)=EPSCRRR(mup1,p1,p2,p3)
       EpsF(66)=EPSCRRR(mup1,p2,p3,p4)
       EpsF(67)=EPSCCCR(mup2,mup3,mup5,p2)
       EpsF(68)=EPSCCCR(mup2,mup3,mup5,p3)
       EpsF(69)=EPSCCRR(mup2,mup3,p1,p3)
       EpsF(70)=EPSCCRR(mup2,mup3,p3,p4)
       EpsF(71)=EPSCCRR(mup2,mup5,p2,p3)
       EpsF(72)=EPSCRRR(mup2,p1,p2,p3)
       EpsF(73)=EPSCRRR(mup2,p2,p3,p4)
       EpsF(74)=EPSCCRR(mup3,mup5,p2,p3)
       EpsF(75)=EPSCRRR(mup3,p1,p2,p3)
       EpsF(76)=EPSCRRR(mup3,p2,p3,p4)
       EpsF(77)=EPSCCRR(mup1,mup2,p1,p4)
       EpsF(78)=EPSCRRR(mup1,p1,p2,p4)
       EpsF(79)=EPSCCRR(mup2,mup3,p1,p4)
       EpsF(80)=EPSCRRR(mup2,p1,p2,p4)
       EpsF(81)=EPSCRRR(mup3,p1,p2,p4)
       EpsF(82)=EPSCCCR(mup2,mup3,mup5,p1)
       EpsF(83)=EPSCCCR(mup2,mup4,mup5,p1)
       EpsF(84)=EPSCCCR(mup2,mup4,mup5,p2)
       EpsF(85)=EPSCCCR(mup2,mup4,mup5,p3)
       EpsF(86)=EPSCCCR(mup3,mup4,mup5,p1)
       EpsF(87)=EPSCCCR(mup3,mup4,mup5,p2)
       EpsF(88)=EPSCCCR(mup3,mup4,mup5,p3)
       EpsF(89)=EPSCCCC(mup2,mup3,mup4,mup5)
       EpsF(90)=EPSCCRR(mup2,mup4,p1,p3)
       EpsF(91)=EPSCCRR(mup3,mup4,p1,p3)
       EpsF(92)=EPSCRRR(mup4,p1,p2,p3)
       EpsF(93)=EPSCCCR(mup2,mup3,mup5,p4)
       EpsF(94)=EPSCCCR(mup2,mup4,mup5,p4)
       EpsF(95)=EPSCCRR(mup2,mup4,p1,p4)
       EpsF(96)=EPSCCRR(mup2,mup4,p3,p4)
       EpsF(97)=EPSCCCR(mup3,mup4,mup5,p4)
       EpsF(98)=EPSCCRR(mup3,mup4,p1,p4)
       EpsF(99)=EPSCCRR(mup3,mup4,p3,p4)
       EpsF(100)=EPSCCRR(mup2,mup5,p2,p4)
       EpsF(101)=EPSCCRR(mup4,mup5,p2,p4)
       EpsF(102)=EPSCRRR(mup4,p1,p2,p4)
       EpsF(103)=EPSCRRR(mup4,p2,p3,p4)
       EpsF(104)=EPSCCRR(mup4,mup5,p2,p3)
       EpsF(105)=EPSCRRR(mup1,p1,p3,p4)
       EpsF(106)=EPSCCRR(mup2,mup5,p3,p4)
       EpsF(107)=EPSCRRR(mup2,p1,p3,p4)
       EpsF(108)=EPSCCRR(mup4,mup5,p3,p4)
       EpsF(109)=EPSCRRR(mup4,p1,p3,p4)
       EpsF(110)=EPSCRRR(mup3,p1,p3,p4)
       EpsF(111)=EPSCCRR(mup3,mup5,p3,p4)
       EpsF(112)=EPSCCRR(mup3,mup5,p2,p4)
       EpsF(113)=EPSRRRR(p1,p2,p3,p4)
       EpsF(114)=EPSCRRR(mup5,p2,p3,p4)
       EpsF(115)=EPSCCRR(mup2,mup5,p1,p2)
       EpsF(116)=EPSCCRR(mup4,mup5,p1,p2)
       EpsF(117)=EPSCCRR(mup2,mup5,p1,p3)
       EpsF(118)=EPSCCRR(mup4,mup5,p1,p3)
       EpsF(119)=EPSCCRR(mup3,mup5,p1,p3)
       EpsF(120)=EPSCCRR(mup3,mup5,p1,p2)
       EpsF(121)=EPSCRRR(mup5,p1,p2,p3)
       EpsF(122)=EPSCCRR(mup2,mup5,p1,p4)
       EpsF(123)=EPSCCRR(mup3,mup5,p1,p4)
       EpsF(124)=EPSCRRR(mup5,p1,p2,p4)
       EpsF(125)=EPSCRRR(mup5,p1,p3,p4)
       EpsF(126)=EPSCCRR(mup4,mup5,p1,p4)
c************************************************************************************
c************************************************************************************
       call FaPenOddGSim1(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p2
     -   mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup5
     -   ,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5m
     -   up4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup2
     -   mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(1))
       call FaPenOddGSim2(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p2
     -   mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup5
     -   ,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5m
     -   up4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup2
     -   mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(100))
       call FaPenOddGSim3(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p2
     -   mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup5
     -   ,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5m
     -   up4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup2
     -   mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(199))
       call FaPenOddGSim4(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p2
     -   mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup5
     -   ,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5m
     -   up4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup2
     -   mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(298))
       call FaPenOddGSim5(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p2
     -   mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup5
     -   ,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5m
     -   up4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup2
     -   mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(397))
       call FaPenOddGSim6(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p2
     -   mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup5
     -   ,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5m
     -   up4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup2
     -   mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(496))
       call FaPenOddGSim7(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p2
     -   mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup5
     -   ,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5m
     -   up4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup2
     -   mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(595))
       call FaPenOddGSim8(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p2
     -   mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup5
     -   ,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5m
     -   up4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup2
     -   mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(694))
       call FaPenOddGSim9(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p2
     -   mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup5
     -   ,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5m
     -   up4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup2
     -   mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(793))
       call FaPenOddGSim10(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p
     -   2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup
     -   5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5
     -   mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup
     -   2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(892))
       call FaPenOddGSim11(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p
     -   2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup
     -   5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5
     -   mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup
     -   2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(991))

       if (comp.gt.0) then
       call FaPenOddGSim12(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p
     -   2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup
     -   5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5
     -   mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup
     -   2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(1090))
       call FaPenOddGSim13(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p
     -   2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup
     -   5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5
     -   mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup
     -   2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(1189))
       call FaPenOddGSim14(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p
     -   2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup
     -   5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5
     -   mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup
     -   2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(1288))
       call FaPenOddGSim15(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p
     -   2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup
     -   5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5
     -   mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup
     -   2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(1387))
       call FaPenOddGSim16(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p
     -   2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup
     -   5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5
     -   mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup
     -   2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(1486))
       call FaPenOddGSim17(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p
     -   2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup
     -   5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5
     -   mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup
     -   2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(1585))
       call FaPenOddGSim18(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p
     -   2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup
     -   5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5
     -   mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup
     -   2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(1684))
       call FaPenOddGSim19(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p
     -   2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup
     -   5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5
     -   mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup
     -   2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(1783))
       call FaPenOddGSim20(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p
     -   2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup
     -   5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5
     -   mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup
     -   2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(1882))
       call FaPenOddGSim21(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p
     -   2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup
     -   5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5
     -   mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup
     -   2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(1981))
       endif

       call FaPenOddGSim22(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p
     -   2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup
     -   5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5
     -   mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup
     -   2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(2080))
       call FaPenOddGSim23(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p
     -   2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup
     -   5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5
     -   mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup
     -   2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(2179))
       call FaPenOddGSim24(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p
     -   2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup
     -   5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5
     -   mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup
     -   2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(2278))
       call FaPenOddGSim25(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p
     -   2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup
     -   5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5
     -   mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup
     -   2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(2377))
       call FaPenOddGSim26(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p
     -   2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup
     -   5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5
     -   mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup
     -   2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(2476))
       call FaPenOddGSim27(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p
     -   2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup
     -   5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5
     -   mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup
     -   2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(2575))
       if(comp.gt.0) then
       call FaPenOddGSim28(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p
     -   2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup
     -   5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5
     -   mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup
     -   2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(2674))
       call FaPenOddGSim29(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p
     -   2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup
     -   5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5
     -   mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup
     -   2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(2773))
       call FaPenOddGSim30(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p
     -   2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup
     -   5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5
     -   mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup
     -   2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(2872))
       call FaPenOddGSim31(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p
     -   2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup
     -   5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5
     -   mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup
     -   2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(2971))
       call FaPenOddGSim32(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p
     -   2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup
     -   5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5
     -   mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup
     -   2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(3070))
       endif
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
       result(1,1) = (0,-0.5)*(Cpl(3)*Cpl(4)*Cpl(5)*(-(EpsF(20)*Fa(1))
     -   -EpsF(36)*Fa(9)-EpsF(41)*Fa(12)-EpsF(89)*Fa(15)-EpsF(47)*Fa(
     -   2139)-2*(EpsF(32)*Fa(5)+EpsF(33)*Fa(6)+EpsF(34)*Fa(7)+EpsF(3
     -   5)*Fa(8)+EpsF(21)*Fa(92)+EpsF(22)*Fa(93)+EpsF(23)*Fa(94)+Eps
     -   F(24)*Fa(95)+EpsF(28)*Fa(99)+EpsF(29)*Fa(100)+EpsF(30)*Fa(10
     -   1)+EpsF(31)*Fa(102)+EpsF(14)*Fa(213)+EpsF(15)*Fa(214)+EpsF(1
     -   6)*Fa(215)+EpsF(51)*Fa(227)+EpsF(53)*Fa(229)+EpsF(54)*Fa(230
     -   )+EpsF(59)*Fa(316)+EpsF(60)*Fa(317)+EpsF(61)*Fa(318)+EpsF(62
     -   )*Fa(319)+EpsF(63)*Fa(320)+EpsF(64)*Fa(321)+EpsF(13)*Fa(406)
     -   +EpsF(25)*Fa(2134)+EpsF(26)*Fa(2135)+EpsF(27)*Fa(2136)+EpsF(
     -   39)*Fa(2137)+EpsF(40)*Fa(2138)+EpsF(77)*Fa(2140)+EpsF(17)*Fa
     -   (2141)+EpsF(18)*Fa(2142)+EpsF(19)*Fa(2143)+EpsF(37)*Fa(2144)
     -   +EpsF(38)*Fa(2145)+EpsF(43)*Fa(2146)+EpsF(44)*Fa(2147)+EpsF(
     -   46)*Fa(2148)+EpsF(55)*Fa(2149)+EpsF(56)*Fa(2150)+EpsF(42)*Fa
     -   (2151)+EpsF(45)*Fa(2152)+EpsF(48)*Fa(2153)+EpsF(49)*Fa(2154)
     -   +EpsF(50)*Fa(2155)+EpsF(52)*Fa(2156)+EpsF(57)*Fa(2157)+EpsF(
     -   58)*Fa(2158)+EpsF(1)*Fa(2159)+EpsF(2)*Fa(2160)+EpsF(3)*Fa(21
     -   61)+EpsF(4)*Fa(2162)+EpsF(5)*Fa(2163)+EpsF(6)*Fa(2164)+EpsF(
     -   7)*Fa(2165)+EpsF(9)*Fa(2166)+EpsF(69)*Fa(2167)+EpsF(70)*Fa(2
     -   168)+EpsF(79)*Fa(2169)+EpsF(90)*Fa(2170)+EpsF(91)*Fa(2171)+E
     -   psF(95)*Fa(2172)+EpsF(96)*Fa(2173)+EpsF(98)*Fa(2174)+EpsF(99
     -   )*Fa(2175)+EpsF(8)*Fa(2176)+EpsF(10)*Fa(2177)+EpsF(11)*Fa(21
     -   78)+EpsF(12)*Fa(2179)+EpsF(65)*Fa(2180)+EpsF(66)*Fa(2181)+Ep
     -   sF(67)*Fa(2182)+EpsF(68)*Fa(2183)+EpsF(71)*Fa(2184)+EpsF(72)
     -   *Fa(2185)+EpsF(73)*Fa(2186)+EpsF(74)*Fa(2187)+EpsF(75)*Fa(21
     -   88)+EpsF(76)*Fa(2189)+EpsF(78)*Fa(2190)+EpsF(80)*Fa(2191)+Ep
     -   sF(81)*Fa(2192)+EpsF(82)*Fa(2193)+EpsF(83)*Fa(2194)+EpsF(84)
     -   *Fa(2195)+EpsF(85)*Fa(2196)+EpsF(93)*Fa(2197)+EpsF(94)*Fa(21
     -   98)+EpsF(100)*Fa(2199)+EpsF(105)*Fa(2200)+EpsF(106)*Fa(2201)
     -   +EpsF(107)*Fa(2202)+EpsF(110)*Fa(2203)+EpsF(111)*Fa(2204)+Ep
     -   sF(112)*Fa(2205)+EpsF(115)*Fa(2206)+EpsF(117)*Fa(2207)+EpsF(
     -   119)*Fa(2208)+EpsF(120)*Fa(2209)+EpsF(122)*Fa(2210)+EpsF(123
     -   )*Fa(2211)+EpsF(86)*Fa(2212)+EpsF(87)*Fa(2213)+EpsF(88)*Fa(2
     -   214)+EpsF(97)*Fa(2215)+EpsF(101)*Fa(2216)+EpsF(102)*Fa(2217)
     -   +EpsF(103)*Fa(2218)+EpsF(104)*Fa(2219)+EpsF(108)*Fa(2220)+Ep
     -   sF(109)*Fa(2221)+EpsF(116)*Fa(2222)+EpsF(118)*Fa(2223)+EpsF(
     -   126)*Fa(2224)+EpsF(92)*Fa(2225)+EpsF(114)*Fa(2226)+EpsF(121)
     -   *Fa(2227)+EpsF(124)*Fa(2228)+EpsF(125)*Fa(2229))-8*EpsF(113)
     -   *Fa(2230))+Cpl(8)*Cpl(9)*Cpl(10)*(-(EpsF(20)*Fa(74))-EpsF(36
     -   )*Fa(82)-EpsF(41)*Fa(85)-EpsF(89)*Fa(88)+8*EpsF(113)*Fa(2230
     -   )-EpsF(47)*Fa(2492)-2*(EpsF(32)*Fa(78)+EpsF(33)*Fa(79)+EpsF(
     -   34)*Fa(80)+EpsF(35)*Fa(81)+EpsF(21)*Fa(192)+EpsF(22)*Fa(193)
     -   +EpsF(23)*Fa(194)+EpsF(24)*Fa(195)+EpsF(28)*Fa(199)+EpsF(29)
     -   *Fa(200)+EpsF(30)*Fa(201)+EpsF(31)*Fa(202)+EpsF(14)*Fa(284)+
     -   EpsF(15)*Fa(285)+EpsF(16)*Fa(286)+EpsF(51)*Fa(298)+EpsF(53)*
     -   Fa(300)+EpsF(54)*Fa(301)+EpsF(59)*Fa(386)+EpsF(60)*Fa(387)+E
     -   psF(61)*Fa(388)+EpsF(62)*Fa(389)+EpsF(63)*Fa(390)+EpsF(64)*F
     -   a(391)+EpsF(13)*Fa(516)+EpsF(25)*Fa(2487)+EpsF(26)*Fa(2488)+
     -   EpsF(27)*Fa(2489)+EpsF(39)*Fa(2490)+EpsF(40)*Fa(2491)+EpsF(7
     -   7)*Fa(2493)+EpsF(17)*Fa(2494)+EpsF(18)*Fa(2495)+EpsF(19)*Fa(
     -   2496)+EpsF(37)*Fa(2497)+EpsF(38)*Fa(2498)+EpsF(43)*Fa(2499)+
     -   EpsF(44)*Fa(2500)+EpsF(46)*Fa(2501)+EpsF(55)*Fa(2502)+EpsF(5
     -   6)*Fa(2503)+EpsF(42)*Fa(2504)+EpsF(45)*Fa(2505)+EpsF(48)*Fa(
     -   2506)+EpsF(49)*Fa(2507)+EpsF(50)*Fa(2508)+EpsF(52)*Fa(2509)+
     -   EpsF(57)*Fa(2510)+EpsF(58)*Fa(2511)+EpsF(1)*Fa(2512)+EpsF(2)
     -   *Fa(2513)+EpsF(3)*Fa(2514)+EpsF(4)*Fa(2515)+EpsF(5)*Fa(2516)
     -   +EpsF(6)*Fa(2517)+EpsF(7)*Fa(2518)+EpsF(9)*Fa(2519)+EpsF(69)
     -   *Fa(2520)+EpsF(70)*Fa(2521)+EpsF(79)*Fa(2522)+EpsF(90)*Fa(25
     -   23)+EpsF(91)*Fa(2524)+EpsF(95)*Fa(2525)+EpsF(96)*Fa(2526)+Ep
     -   sF(98)*Fa(2527)+EpsF(99)*Fa(2528)+EpsF(8)*Fa(2529)+EpsF(10)*
     -   Fa(2530)+EpsF(11)*Fa(2531)+EpsF(12)*Fa(2532)+EpsF(65)*Fa(253
     -   3)+EpsF(66)*Fa(2534)+EpsF(67)*Fa(2535)+EpsF(68)*Fa(2536)+Eps
     -   F(71)*Fa(2537)+EpsF(72)*Fa(2538)+EpsF(73)*Fa(2539)+EpsF(74)*
     -   Fa(2540)+EpsF(75)*Fa(2541)+EpsF(76)*Fa(2542)+EpsF(78)*Fa(254
     -   3)+EpsF(80)*Fa(2544)+EpsF(81)*Fa(2545)+EpsF(82)*Fa(2546)+Eps
     -   F(83)*Fa(2547)+EpsF(84)*Fa(2548)+EpsF(85)*Fa(2549)+EpsF(93)*
     -   Fa(2550)+EpsF(94)*Fa(2551)+EpsF(100)*Fa(2552)+EpsF(105)*Fa(2
     -   553)+EpsF(106)*Fa(2554)+EpsF(107)*Fa(2555)+EpsF(110)*Fa(2556
     -   )+EpsF(111)*Fa(2557)+EpsF(112)*Fa(2558)+EpsF(115)*Fa(2559)+E
     -   psF(117)*Fa(2560)+EpsF(119)*Fa(2561)+EpsF(120)*Fa(2562)+EpsF
     -   (122)*Fa(2563)+EpsF(123)*Fa(2564)+EpsF(86)*Fa(2565)+EpsF(87)
     -   *Fa(2566)+EpsF(88)*Fa(2567)+EpsF(97)*Fa(2568)+EpsF(101)*Fa(2
     -   569)+EpsF(102)*Fa(2570)+EpsF(103)*Fa(2571)+EpsF(104)*Fa(2572
     -   )+EpsF(108)*Fa(2573)+EpsF(109)*Fa(2574)+EpsF(116)*Fa(2575)+E
     -   psF(118)*Fa(2576)+EpsF(126)*Fa(2577)+EpsF(92)*Fa(2578)+EpsF(
     -   114)*Fa(2579)+EpsF(121)*Fa(2580)+EpsF(124)*Fa(2581)+EpsF(125
     -   )*Fa(2582)))+M3*M4*(Cpl(3)*Cpl(5)*Cpl(9)*(-2*(EpsF(35)*Fa(23
     -   )+EpsF(20)*Fa(29)+EpsF(32)*Fa(32)+EpsF(33)*Fa(33)+EpsF(34)*F
     -   a(34)+EpsF(36)*Fa(35)+EpsF(41)*Fa(36)+EpsF(89)*Fa(38)+EpsF(2
     -   4)*Fa(117)+EpsF(31)*Fa(124)+EpsF(21)*Fa(130)+EpsF(22)*Fa(131
     -   )+EpsF(23)*Fa(132)+EpsF(28)*Fa(135)+EpsF(29)*Fa(136)+EpsF(30
     -   )*Fa(137)-EpsF(13)*Fa(448)-EpsF(12)*Fa(469)+EpsF(56)*Fa(2238
     -   )+EpsF(52)*Fa(2244)+EpsF(58)*Fa(2246)+EpsF(93)*Fa(2273)+EpsF
     -   (94)*Fa(2274)+EpsF(97)*Fa(2285)+EpsF(47)*Fa(2294)+EpsF(37)*F
     -   a(2295)+EpsF(38)*Fa(2296)+EpsF(55)*Fa(2297)+EpsF(42)*Fa(2298
     -   )+EpsF(45)*Fa(2299)+EpsF(48)*Fa(2300)+EpsF(49)*Fa(2301)+EpsF
     -   (50)*Fa(2302)+EpsF(57)*Fa(2303)-EpsF(10)*Fa(2311)+EpsF(11)*F
     -   a(2312)+EpsF(67)*Fa(2313)+EpsF(68)*Fa(2314)+EpsF(82)*Fa(2317
     -   )+EpsF(83)*Fa(2318)+EpsF(84)*Fa(2319)+EpsF(85)*Fa(2320)+EpsF
     -   (86)*Fa(2324)+EpsF(87)*Fa(2325)+EpsF(88)*Fa(2326))-4*(EpsF(1
     -   6)*Fa(235)-EpsF(19)*Fa(238)-EpsF(46)*Fa(242)+EpsF(53)*Fa(248
     -   )+EpsF(14)*Fa(251)+EpsF(15)*Fa(252)+EpsF(17)*Fa(253)-EpsF(18
     -   )*Fa(254)-EpsF(43)*Fa(256)+EpsF(51)*Fa(261)-EpsF(61)*Fa(343)
     -   -EpsF(63)*Fa(345)+EpsF(59)*Fa(353)-EpsF(60)*Fa(354)-EpsF(62)
     -   *Fa(355)+EpsF(27)*Fa(2233)-EpsF(3)*Fa(2249)+EpsF(6)*Fa(2252)
     -   +EpsF(9)*Fa(2254)+EpsF(98)*Fa(2258)+EpsF(100)*Fa(2275)+EpsF(
     -   112)*Fa(2276)+EpsF(123)*Fa(2281)+EpsF(101)*Fa(2286)-EpsF(126
     -   )*Fa(2291)-EpsF(25)*Fa(2292)+EpsF(26)*Fa(2293)-EpsF(1)*Fa(23
     -   04)-EpsF(2)*Fa(2305)+EpsF(4)*Fa(2306)+EpsF(5)*Fa(2307)-EpsF(
     -   7)*Fa(2308)+EpsF(91)*Fa(2309)+EpsF(8)*Fa(2310)+EpsF(71)*Fa(2
     -   315)+EpsF(74)*Fa(2316)-EpsF(115)*Fa(2321)+EpsF(119)*Fa(2322)
     -   +EpsF(120)*Fa(2323)+EpsF(104)*Fa(2327)-EpsF(116)*Fa(2328)-Ep
     -   sF(118)*Fa(2329)-EpsF(39)*Fa(25)*F(4463)-EpsF(121)*Fa(731)*F
     -   (4463)+EpsF(69)*Fa(2255)*F(4463)+EpsF(90)*Fa(2256)*F(4463)-E
     -   psF(65)*Fa(2262)*F(4463)+EpsF(72)*Fa(2266)*F(4463)-EpsF(75)*
     -   Fa(2268)*F(4463)+EpsF(117)*Fa(2278)*F(4463)+EpsF(92)*Fa(2287
     -   )*F(4463)-EpsF(77)*Fa(25)*F(4464)-EpsF(124)*Fa(731)*F(4464)+
     -   EpsF(79)*Fa(2255)*F(4464)+EpsF(95)*Fa(2256)*F(4464)-EpsF(78)
     -   *Fa(2262)*F(4464)+EpsF(80)*Fa(2266)*F(4464)-EpsF(81)*Fa(2268
     -   )*F(4464)+EpsF(122)*Fa(2278)*F(4464)+EpsF(102)*Fa(2287)*F(44
     -   64)))+Cpl(4)*Cpl(8)*Cpl(10)*(-2*(EpsF(20)*Fa(58)+EpsF(32)*Fa
     -   (59)+EpsF(33)*Fa(60)+EpsF(34)*Fa(61)+EpsF(35)*Fa(62)+EpsF(36
     -   )*Fa(63)+EpsF(41)*Fa(64)+EpsF(89)*Fa(66)-EpsF(21)*Fa(169)-Ep
     -   sF(22)*Fa(170)-EpsF(23)*Fa(171)-EpsF(24)*Fa(172)+EpsF(28)*Fa
     -   (174)+EpsF(29)*Fa(175)+EpsF(30)*Fa(176)+EpsF(31)*Fa(177)+Eps
     -   F(12)*Fa(502)+EpsF(13)*Fa(503)+EpsF(47)*Fa(2422)+EpsF(37)*Fa
     -   (2423)+EpsF(38)*Fa(2424)+EpsF(55)*Fa(2425)+EpsF(56)*Fa(2426)
     -   +EpsF(42)*Fa(2427)+EpsF(45)*Fa(2428)+EpsF(48)*Fa(2429)+EpsF(
     -   49)*Fa(2430)+EpsF(50)*Fa(2431)+EpsF(52)*Fa(2432)+EpsF(57)*Fa
     -   (2433)+EpsF(58)*Fa(2434)+EpsF(10)*Fa(2442)-EpsF(11)*Fa(2443)
     -   +EpsF(67)*Fa(2444)+EpsF(68)*Fa(2445)+EpsF(82)*Fa(2447)+EpsF(
     -   83)*Fa(2448)+EpsF(84)*Fa(2449)+EpsF(85)*Fa(2450)+EpsF(93)*Fa
     -   (2451)+EpsF(94)*Fa(2452)+EpsF(86)*Fa(2458)+EpsF(87)*Fa(2459)
     -   +EpsF(88)*Fa(2460)+EpsF(97)*Fa(2461))+4*(EpsF(16)*Fa(235)-Ep
     -   sF(19)*Fa(238)-EpsF(46)*Fa(242)+EpsF(53)*Fa(248)+EpsF(15)*Fa
     -   (252)-EpsF(18)*Fa(254)-EpsF(43)*Fa(256)+EpsF(51)*Fa(261)-Eps
     -   F(14)*Fa(278)-EpsF(17)*Fa(279)-EpsF(61)*Fa(343)-EpsF(63)*Fa(
     -   345)-EpsF(60)*Fa(354)-EpsF(62)*Fa(355)-EpsF(59)*Fa(370)-EpsF
     -   (3)*Fa(2249)+EpsF(112)*Fa(2276)-EpsF(126)*Fa(2291)-EpsF(25)*
     -   Fa(2292)+EpsF(26)*Fa(2293)-EpsF(2)*Fa(2305)+EpsF(4)*Fa(2306)
     -   -EpsF(7)*Fa(2308)+EpsF(74)*Fa(2316)-EpsF(118)*Fa(2329)+EpsF(
     -   27)*Fa(2421)-EpsF(1)*Fa(2435)-EpsF(5)*Fa(2436)-EpsF(6)*Fa(24
     -   37)-EpsF(9)*Fa(2438)-EpsF(91)*Fa(2439)-EpsF(98)*Fa(2440)-Eps
     -   F(8)*Fa(2441)-EpsF(71)*Fa(2446)-EpsF(100)*Fa(2453)-EpsF(115)
     -   *Fa(2454)-EpsF(119)*Fa(2455)-EpsF(120)*Fa(2456)-EpsF(123)*Fa
     -   (2457)-EpsF(101)*Fa(2462)-EpsF(104)*Fa(2463)-EpsF(116)*Fa(24
     -   64)-EpsF(39)*Fa(25)*F(4463)-EpsF(121)*Fa(731)*F(4463)-EpsF(6
     -   5)*Fa(2262)*F(4463)+EpsF(72)*Fa(2266)*F(4463)-EpsF(75)*Fa(22
     -   68)*F(4463)-EpsF(69)*Fa(2353)*F(4463)-EpsF(90)*Fa(2354)*F(44
     -   63)-EpsF(117)*Fa(2372)*F(4463)-EpsF(92)*Fa(2381)*F(4463)-Eps
     -   F(77)*Fa(25)*F(4464)-EpsF(124)*Fa(731)*F(4464)-EpsF(78)*Fa(2
     -   262)*F(4464)+EpsF(80)*Fa(2266)*F(4464)-EpsF(81)*Fa(2268)*F(4
     -   464)-EpsF(79)*Fa(2353)*F(4464)-EpsF(95)*Fa(2354)*F(4464)-Eps
     -   F(122)*Fa(2372)*F(4464)-EpsF(102)*Fa(2381)*F(4464))))+M1*(M3
     -   *(Cpl(3)*Cpl(9)*Cpl(10)*(-2*(EpsF(20)*Fa(39)+EpsF(32)*Fa(42)
     -   +EpsF(33)*Fa(43)+EpsF(34)*Fa(44)+EpsF(35)*Fa(45)+EpsF(36)*Fa
     -   (46)+EpsF(41)*Fa(47)+EpsF(89)*Fa(49)+EpsF(21)*Fa(142)+EpsF(2
     -   2)*Fa(143)-EpsF(23)*Fa(144)-EpsF(24)*Fa(145)+EpsF(28)*Fa(148
     -   )-EpsF(29)*Fa(149)+EpsF(30)*Fa(150)+EpsF(31)*Fa(151)+EpsF(12
     -   )*Fa(484)+EpsF(13)*Fa(485)+EpsF(47)*Fa(2332)+EpsF(37)*Fa(233
     -   3)+EpsF(38)*Fa(2334)+EpsF(55)*Fa(2335)+EpsF(56)*Fa(2336)+Eps
     -   F(42)*Fa(2337)+EpsF(45)*Fa(2338)+EpsF(48)*Fa(2339)-EpsF(49)*
     -   Fa(2340)+EpsF(50)*Fa(2341)+EpsF(52)*Fa(2342)+EpsF(57)*Fa(234
     -   3)+EpsF(58)*Fa(2344)+EpsF(10)*Fa(2358)-EpsF(11)*Fa(2359)+Eps
     -   F(67)*Fa(2360)+EpsF(68)*Fa(2361)+EpsF(82)*Fa(2363)+EpsF(83)*
     -   Fa(2364)+EpsF(84)*Fa(2365)+EpsF(85)*Fa(2366)+EpsF(93)*Fa(236
     -   7)+EpsF(94)*Fa(2368)+EpsF(86)*Fa(2376)+EpsF(87)*Fa(2377)+Eps
     -   F(88)*Fa(2378)+EpsF(97)*Fa(2379))+4*(-(EpsF(39)*F(396)*Fa(25
     -   ))-EpsF(77)*F(477)*Fa(25)+EpsF(15)*Fa(234)-EpsF(18)*Fa(237)-
     -   EpsF(43)*Fa(240)+EpsF(51)*Fa(246)-EpsF(14)*Fa(263)+EpsF(16)*
     -   Fa(264)-EpsF(17)*Fa(265)-EpsF(19)*Fa(266)-EpsF(46)*Fa(268)+E
     -   psF(53)*Fa(272)-EpsF(60)*Fa(342)-EpsF(62)*Fa(344)-EpsF(59)*F
     -   a(363)-EpsF(61)*Fa(364)-EpsF(63)*Fa(365)-EpsF(121)*F(396)*Fa
     -   (731)-EpsF(124)*F(477)*Fa(731)+EpsF(26)*Fa(2232)-EpsF(65)*F(
     -   396)*Fa(2262)-EpsF(78)*F(477)*Fa(2262)+EpsF(71)*Fa(2265)+Eps
     -   F(72)*F(396)*Fa(2266)+EpsF(80)*F(477)*Fa(2266)-EpsF(75)*F(39
     -   6)*Fa(2268)-EpsF(81)*F(477)*Fa(2268)-EpsF(118)*Fa(2290)-EpsF
     -   (25)*Fa(2330)+EpsF(27)*Fa(2331)-EpsF(1)*Fa(2345)-EpsF(2)*Fa(
     -   2346)-EpsF(3)*Fa(2347)-EpsF(4)*Fa(2348)-EpsF(5)*Fa(2349)-Eps
     -   F(6)*Fa(2350)-EpsF(7)*Fa(2351)-EpsF(9)*Fa(2352)-EpsF(69)*F(3
     -   96)*Fa(2353)-EpsF(79)*F(477)*Fa(2353)-EpsF(90)*F(396)*Fa(235
     -   4)-EpsF(95)*F(477)*Fa(2354)-EpsF(91)*Fa(2355)-EpsF(98)*Fa(23
     -   56)-EpsF(8)*Fa(2357)-EpsF(74)*Fa(2362)+EpsF(100)*Fa(2369)-Ep
     -   sF(112)*Fa(2370)-EpsF(115)*Fa(2371)-EpsF(117)*F(396)*Fa(2372
     -   )-EpsF(122)*F(477)*Fa(2372)-EpsF(119)*Fa(2373)-EpsF(120)*Fa(
     -   2374)-EpsF(123)*Fa(2375)-EpsF(101)*Fa(2380)-EpsF(92)*F(396)*
     -   Fa(2381)-EpsF(102)*F(477)*Fa(2381)-EpsF(104)*Fa(2382)-EpsF(1
     -   16)*Fa(2383)-EpsF(126)*Fa(2384)))+Cpl(4)*Cpl(5)*Cpl(8)*(-2*(
     -   EpsF(34)*Fa(22)+EpsF(20)*Fa(50)+EpsF(32)*Fa(51)-EpsF(33)*Fa(
     -   52)+EpsF(35)*Fa(53)+EpsF(36)*Fa(54)+EpsF(41)*Fa(55)+EpsF(89)
     -   *Fa(57)+EpsF(23)*Fa(116)+EpsF(30)*Fa(123)+EpsF(21)*Fa(157)+E
     -   psF(22)*Fa(158)+EpsF(24)*Fa(159)+EpsF(28)*Fa(162)+EpsF(29)*F
     -   a(163)+EpsF(31)*Fa(164)-EpsF(12)*Fa(447)-EpsF(13)*Fa(499)+Ep
     -   sF(38)*Fa(2236)+EpsF(42)*Fa(2239)+EpsF(50)*Fa(2243)+EpsF(68)
     -   *Fa(2264)+EpsF(85)*Fa(2272)+EpsF(88)*Fa(2284)+EpsF(47)*Fa(23
     -   87)+EpsF(37)*Fa(2388)+EpsF(55)*Fa(2389)+EpsF(56)*Fa(2390)-Ep
     -   sF(45)*Fa(2391)+EpsF(48)*Fa(2392)+EpsF(49)*Fa(2393)+EpsF(52)
     -   *Fa(2394)+EpsF(57)*Fa(2395)+EpsF(58)*Fa(2396)+EpsF(10)*Fa(24
     -   04)+EpsF(11)*Fa(2405)+EpsF(67)*Fa(2406)+EpsF(82)*Fa(2407)+Ep
     -   sF(83)*Fa(2408)+EpsF(84)*Fa(2409)+EpsF(93)*Fa(2410)+EpsF(94)
     -   *Fa(2411)+EpsF(86)*Fa(2416)+EpsF(87)*Fa(2417)+EpsF(97)*Fa(24
     -   18))-4*(-(EpsF(39)*F(396)*Fa(25))-EpsF(77)*F(477)*Fa(25)+Eps
     -   F(15)*Fa(234)-EpsF(18)*Fa(237)-EpsF(43)*Fa(240)+EpsF(51)*Fa(
     -   246)+EpsF(16)*Fa(264)-EpsF(19)*Fa(266)-EpsF(46)*Fa(268)+EpsF
     -   (53)*Fa(272)+EpsF(14)*Fa(275)+EpsF(17)*Fa(276)-EpsF(60)*Fa(3
     -   42)-EpsF(62)*Fa(344)-EpsF(61)*Fa(364)-EpsF(63)*Fa(365)+EpsF(
     -   59)*Fa(369)-EpsF(121)*F(396)*Fa(731)-EpsF(124)*F(477)*Fa(731
     -   )+EpsF(26)*Fa(2232)+EpsF(2)*Fa(2248)+EpsF(5)*Fa(2251)+EpsF(6
     -   9)*F(396)*Fa(2255)+EpsF(79)*F(477)*Fa(2255)+EpsF(90)*F(396)*
     -   Fa(2256)+EpsF(95)*F(477)*Fa(2256)+EpsF(91)*Fa(2257)+EpsF(8)*
     -   Fa(2259)-EpsF(65)*F(396)*Fa(2262)-EpsF(78)*F(477)*Fa(2262)+E
     -   psF(71)*Fa(2265)+EpsF(72)*F(396)*Fa(2266)+EpsF(80)*F(477)*Fa
     -   (2266)+EpsF(74)*Fa(2267)-EpsF(75)*F(396)*Fa(2268)-EpsF(81)*F
     -   (477)*Fa(2268)+EpsF(117)*F(396)*Fa(2278)+EpsF(122)*F(477)*Fa
     -   (2278)+EpsF(119)*Fa(2279)+EpsF(92)*F(396)*Fa(2287)+EpsF(102)
     -   *F(477)*Fa(2287)+EpsF(104)*Fa(2288)-EpsF(118)*Fa(2290)+EpsF(
     -   100)*Fa(2369)-EpsF(126)*Fa(2384)+EpsF(25)*Fa(2385)+EpsF(27)*
     -   Fa(2386)+EpsF(1)*Fa(2397)+EpsF(3)*Fa(2398)+EpsF(4)*Fa(2399)+
     -   EpsF(6)*Fa(2400)+EpsF(7)*Fa(2401)+EpsF(9)*Fa(2402)+EpsF(98)*
     -   Fa(2403)+EpsF(112)*Fa(2412)-EpsF(115)*Fa(2413)+EpsF(120)*Fa(
     -   2414)+EpsF(123)*Fa(2415)+EpsF(101)*Fa(2419)-EpsF(116)*Fa(242
     -   0))))+M4*(Cpl(3)*Cpl(4)*Cpl(10)*(-2*(EpsF(20)*Fa(16)+EpsF(32
     -   )*Fa(20)+EpsF(33)*Fa(21)+EpsF(34)*Fa(22)+EpsF(35)*Fa(23)+Eps
     -   F(36)*Fa(24)+EpsF(41)*Fa(26)+EpsF(89)*Fa(28)+EpsF(21)*Fa(114
     -   )+EpsF(22)*Fa(115)+EpsF(23)*Fa(116)+EpsF(24)*Fa(117)+EpsF(28
     -   )*Fa(121)+EpsF(29)*Fa(122)+EpsF(30)*Fa(123)+EpsF(31)*Fa(124)
     -   -EpsF(12)*Fa(447)-EpsF(13)*Fa(448)+EpsF(47)*Fa(2234)+EpsF(37
     -   )*Fa(2235)+EpsF(38)*Fa(2236)+EpsF(55)*Fa(2237)+EpsF(56)*Fa(2
     -   238)+EpsF(42)*Fa(2239)+EpsF(45)*Fa(2240)+EpsF(48)*Fa(2241)+E
     -   psF(49)*Fa(2242)+EpsF(50)*Fa(2243)+EpsF(52)*Fa(2244)+EpsF(57
     -   )*Fa(2245)+EpsF(58)*Fa(2246)+EpsF(10)*Fa(2260)+EpsF(11)*Fa(2
     -   261)+EpsF(67)*Fa(2263)+EpsF(68)*Fa(2264)+EpsF(82)*Fa(2269)+E
     -   psF(83)*Fa(2270)+EpsF(84)*Fa(2271)+EpsF(85)*Fa(2272)+EpsF(93
     -   )*Fa(2273)+EpsF(94)*Fa(2274)+EpsF(86)*Fa(2282)+EpsF(87)*Fa(2
     -   283)+EpsF(88)*Fa(2284)+EpsF(97)*Fa(2285))-4*(-(EpsF(39)*F(39
     -   6)*Fa(25))+EpsF(14)*Fa(233)+EpsF(15)*Fa(234)+EpsF(16)*Fa(235
     -   )+EpsF(17)*Fa(236)-EpsF(18)*Fa(237)-EpsF(19)*Fa(238)-EpsF(43
     -   )*Fa(240)-EpsF(46)*Fa(242)+EpsF(51)*Fa(246)+EpsF(53)*Fa(248)
     -   +EpsF(59)*Fa(341)-EpsF(60)*Fa(342)-EpsF(61)*Fa(343)-EpsF(62)
     -   *Fa(344)-EpsF(63)*Fa(345)-EpsF(121)*F(396)*Fa(731)+EpsF(25)*
     -   Fa(2231)+EpsF(26)*Fa(2232)+EpsF(27)*Fa(2233)-EpsF(1)*Fa(2247
     -   )+EpsF(2)*Fa(2248)-EpsF(3)*Fa(2249)+EpsF(4)*Fa(2250)+EpsF(5)
     -   *Fa(2251)+EpsF(6)*Fa(2252)-EpsF(7)*Fa(2253)+EpsF(9)*Fa(2254)
     -   +EpsF(69)*F(396)*Fa(2255)+EpsF(90)*F(396)*Fa(2256)+EpsF(91)*
     -   Fa(2257)+EpsF(98)*Fa(2258)+EpsF(8)*Fa(2259)-EpsF(65)*F(396)*
     -   Fa(2262)+EpsF(71)*Fa(2265)+EpsF(72)*F(396)*Fa(2266)+EpsF(74)
     -   *Fa(2267)-EpsF(75)*F(396)*Fa(2268)+EpsF(100)*Fa(2275)+EpsF(1
     -   12)*Fa(2276)+EpsF(115)*Fa(2277)+EpsF(117)*F(396)*Fa(2278)+Ep
     -   sF(119)*Fa(2279)+EpsF(120)*Fa(2280)+EpsF(123)*Fa(2281)+EpsF(
     -   101)*Fa(2286)+EpsF(92)*F(396)*Fa(2287)+EpsF(104)*Fa(2288)-Ep
     -   sF(116)*Fa(2289)-EpsF(118)*Fa(2290)-EpsF(126)*Fa(2291)-EpsF(
     -   77)*Fa(25)*F(4464)-EpsF(124)*Fa(731)*F(4464)+EpsF(79)*Fa(225
     -   5)*F(4464)+EpsF(95)*Fa(2256)*F(4464)-EpsF(78)*Fa(2262)*F(446
     -   4)+EpsF(80)*Fa(2266)*F(4464)-EpsF(81)*Fa(2268)*F(4464)+EpsF(
     -   122)*Fa(2278)*F(4464)+EpsF(102)*Fa(2287)*F(4464)))+Cpl(5)*Cp
     -   l(8)*Cpl(9)*(-2*(EpsF(34)*Fa(44)+EpsF(35)*Fa(62)+EpsF(20)*Fa
     -   (67)+EpsF(32)*Fa(68)+EpsF(33)*Fa(69)+EpsF(36)*Fa(70)+EpsF(41
     -   )*Fa(71)+EpsF(89)*Fa(73)-EpsF(23)*Fa(144)+EpsF(30)*Fa(150)-E
     -   psF(24)*Fa(172)+EpsF(31)*Fa(177)-EpsF(21)*Fa(183)+EpsF(22)*F
     -   a(184)+EpsF(28)*Fa(186)-EpsF(29)*Fa(187)+EpsF(12)*Fa(484)+Ep
     -   sF(13)*Fa(503)+EpsF(38)*Fa(2334)+EpsF(42)*Fa(2337)+EpsF(50)*
     -   Fa(2341)+EpsF(68)*Fa(2361)+EpsF(85)*Fa(2366)+EpsF(88)*Fa(237
     -   8)+EpsF(56)*Fa(2426)+EpsF(52)*Fa(2432)+EpsF(58)*Fa(2434)+Eps
     -   F(93)*Fa(2451)+EpsF(94)*Fa(2452)+EpsF(97)*Fa(2461)+EpsF(47)*
     -   Fa(2466)+EpsF(37)*Fa(2467)+EpsF(55)*Fa(2468)+EpsF(45)*Fa(246
     -   9)+EpsF(48)*Fa(2470)-EpsF(49)*Fa(2471)+EpsF(57)*Fa(2472)+Eps
     -   F(10)*Fa(2476)-EpsF(11)*Fa(2477)+EpsF(67)*Fa(2478)+EpsF(82)*
     -   Fa(2479)+EpsF(83)*Fa(2480)+EpsF(84)*Fa(2481)+EpsF(86)*Fa(248
     -   4)+EpsF(87)*Fa(2485))+4*(-(EpsF(39)*F(396)*Fa(25))+EpsF(15)*
     -   Fa(234)+EpsF(16)*Fa(235)-EpsF(18)*Fa(237)-EpsF(19)*Fa(238)-E
     -   psF(43)*Fa(240)-EpsF(46)*Fa(242)+EpsF(51)*Fa(246)+EpsF(53)*F
     -   a(248)-EpsF(14)*Fa(280)-EpsF(17)*Fa(281)-EpsF(60)*Fa(342)-Ep
     -   sF(61)*Fa(343)-EpsF(62)*Fa(344)-EpsF(63)*Fa(345)-EpsF(59)*Fa
     -   (373)-EpsF(121)*F(396)*Fa(731)+EpsF(26)*Fa(2232)-EpsF(3)*Fa(
     -   2249)-EpsF(65)*F(396)*Fa(2262)+EpsF(71)*Fa(2265)+EpsF(72)*F(
     -   396)*Fa(2266)-EpsF(75)*F(396)*Fa(2268)+EpsF(112)*Fa(2276)-Ep
     -   sF(118)*Fa(2290)-EpsF(126)*Fa(2291)-EpsF(2)*Fa(2346)-EpsF(5)
     -   *Fa(2349)-EpsF(69)*F(396)*Fa(2353)-EpsF(90)*F(396)*Fa(2354)-
     -   EpsF(91)*Fa(2355)-EpsF(8)*Fa(2357)-EpsF(74)*Fa(2362)-EpsF(11
     -   7)*F(396)*Fa(2372)-EpsF(119)*Fa(2373)-EpsF(92)*F(396)*Fa(238
     -   1)-EpsF(104)*Fa(2382)+EpsF(27)*Fa(2421)-EpsF(6)*Fa(2437)-Eps
     -   F(9)*Fa(2438)-EpsF(98)*Fa(2440)-EpsF(100)*Fa(2453)-EpsF(123)
     -   *Fa(2457)-EpsF(101)*Fa(2462)-EpsF(25)*Fa(2465)-EpsF(1)*Fa(24
     -   73)+EpsF(4)*Fa(2474)-EpsF(7)*Fa(2475)-EpsF(115)*Fa(2482)-Eps
     -   F(120)*Fa(2483)-EpsF(116)*Fa(2486)-EpsF(77)*Fa(25)*F(4464)-E
     -   psF(124)*Fa(731)*F(4464)-EpsF(78)*Fa(2262)*F(4464)+EpsF(80)*
     -   Fa(2266)*F(4464)-EpsF(81)*Fa(2268)*F(4464)-EpsF(79)*Fa(2353)
     -   *F(4464)-EpsF(95)*Fa(2354)*F(4464)-EpsF(122)*Fa(2372)*F(4464
     -   )-EpsF(102)*Fa(2381)*F(4464))))))
       if(comp.gt.0) then
       result(1,2) = 0
       result(1,3) = (0,-0.5)*(Cpl(3)*Cpl(4)*Cpl(5)*(-(EpsF(41)*F(2482
     -   ))-8*(EpsF(66)*Fa(2587)+EpsF(78)*Fa(2588)+EpsF(105)*Fa(2589)
     -   -EpsF(102)*Fa(2596)+EpsF(103)*Fa(2597)-EpsF(109)*Fa(2598)-Ep
     -   sF(76)*Fa(2600)+EpsF(81)*Fa(2601)+EpsF(110)*Fa(2604)-EpsF(11
     -   4)*Fa(2607)+EpsF(124)*Fa(2611)+EpsF(125)*Fa(2612))-2*(EpsF(1
     -   7)*Fa(1021)+EpsF(18)*Fa(1022)+EpsF(19)*Fa(1023)+EpsF(21)*Fa(
     -   1024)+EpsF(22)*Fa(1025)+EpsF(23)*Fa(1026)+EpsF(24)*Fa(1027)+
     -   EpsF(42)*Fa(1028)+EpsF(43)*Fa(1029)+EpsF(44)*Fa(1030)+EpsF(4
     -   5)*Fa(1031)+EpsF(46)*Fa(1032)+EpsF(57)*Fa(1037)+EpsF(58)*Fa(
     -   1038)+EpsF(86)*Fa(1041)+EpsF(87)*Fa(1042)+EpsF(88)*Fa(1043)+
     -   EpsF(97)*Fa(1044)+EpsF(14)*Fa(1131)+EpsF(15)*Fa(1132)+EpsF(1
     -   6)*Fa(1133)+EpsF(51)*Fa(1137)+EpsF(53)*Fa(1139)+EpsF(54)*Fa(
     -   1140)+EpsF(3)*Fa(1187)+EpsF(98)*Fa(1191)+EpsF(99)*Fa(1192)+E
     -   psF(48)*Fa(2583)+EpsF(49)*Fa(2584)+EpsF(50)*Fa(2585)+EpsF(52
     -   )*Fa(2586)+EpsF(59)*Fa(2590)+EpsF(60)*Fa(2591)+EpsF(61)*Fa(2
     -   592)+EpsF(62)*Fa(2593)+EpsF(63)*Fa(2594)+EpsF(64)*Fa(2595)+E
     -   psF(74)*Fa(2599)+EpsF(101)*Fa(2602)+EpsF(108)*Fa(2603)+EpsF(
     -   111)*Fa(2605)+EpsF(112)*Fa(2606)+EpsF(119)*Fa(2608)+EpsF(120
     -   )*Fa(2609)+EpsF(123)*Fa(2610)+EpsF(126)*Fa(2613)+EpsF(1)*Fa(
     -   2614)+EpsF(2)*Fa(2615)+EpsF(91)*Fa(2616)+EpsF(104)*Fa(2617)+
     -   EpsF(116)*Fa(2618)+EpsF(118)*Fa(2619)))+Cpl(8)*Cpl(9)*Cpl(10
     -   )*(-(EpsF(41)*F(2674))+8*(EpsF(66)*Fa(2587)+EpsF(78)*Fa(2588
     -   )+EpsF(105)*Fa(2589)+EpsF(103)*Fa(2597)-EpsF(109)*Fa(2598)-E
     -   psF(76)*Fa(2600)+EpsF(81)*Fa(2601)+EpsF(110)*Fa(2604)-EpsF(1
     -   14)*Fa(2607)+EpsF(124)*Fa(2611)+EpsF(125)*Fa(2612)-EpsF(102)
     -   *Fa(2650))-2*(EpsF(17)*Fa(1109)+EpsF(18)*Fa(1110)+EpsF(19)*F
     -   a(1111)+EpsF(21)*Fa(1112)+EpsF(22)*Fa(1113)+EpsF(23)*Fa(1114
     -   )+EpsF(24)*Fa(1115)+EpsF(42)*Fa(1116)+EpsF(43)*Fa(1117)+EpsF
     -   (44)*Fa(1118)+EpsF(45)*Fa(1119)+EpsF(46)*Fa(1120)+EpsF(57)*F
     -   a(1125)+EpsF(58)*Fa(1126)+EpsF(86)*Fa(1127)+EpsF(87)*Fa(1128
     -   )+EpsF(88)*Fa(1129)+EpsF(97)*Fa(1130)+EpsF(14)*Fa(1168)+EpsF
     -   (15)*Fa(1169)+EpsF(16)*Fa(1170)+EpsF(51)*Fa(1174)+EpsF(53)*F
     -   a(1176)+EpsF(54)*Fa(1177)+EpsF(3)*Fa(1205)+EpsF(98)*Fa(1206)
     -   +EpsF(99)*Fa(1207)+EpsF(48)*Fa(2640)+EpsF(49)*Fa(2641)+EpsF(
     -   50)*Fa(2642)+EpsF(52)*Fa(2643)+EpsF(59)*Fa(2644)+EpsF(60)*Fa
     -   (2645)+EpsF(61)*Fa(2646)+EpsF(62)*Fa(2647)+EpsF(63)*Fa(2648)
     -   +EpsF(64)*Fa(2649)+EpsF(74)*Fa(2651)+EpsF(101)*Fa(2652)+EpsF
     -   (108)*Fa(2653)+EpsF(111)*Fa(2654)+EpsF(112)*Fa(2655)+EpsF(11
     -   9)*Fa(2656)+EpsF(120)*Fa(2657)+EpsF(123)*Fa(2658)+EpsF(126)*
     -   Fa(2659)+EpsF(1)*Fa(2660)+EpsF(2)*Fa(2661)+EpsF(91)*Fa(2662)
     -   +EpsF(104)*Fa(2663)+EpsF(116)*Fa(2664)+EpsF(118)*Fa(2665)))+
     -   M3*M4*(Cpl(3)*Cpl(5)*Cpl(9)*(-2*EpsF(41)*F(2622)-4*(mup4mup5
     -   *EpsF(46)*F(2605)+mup3mup4*EpsF(61)*F(2605)+mup1mup5*EpsF(98
     -   )*F(2605)+mup1mup3*EpsF(126)*F(2605)+mup3mup5*EpsF(53)*F(261
     -   6)+mup1mup4*EpsF(123)*F(2616)+mup4mup5*EpsF(43)*F(2625)+mup3
     -   mup4*EpsF(60)*F(2625)+mup1mup5*EpsF(91)*F(2625)+mup1mup3*Eps
     -   F(118)*F(2625)+mup3mup5*EpsF(51)*F(2627)+mup1mup4*EpsF(119)*
     -   F(2627)+EpsF(24)*Fa(1049)+EpsF(58)*Fa(1057)+EpsF(21)*Fa(1060
     -   )-EpsF(22)*Fa(1061)+EpsF(23)*Fa(1062)+EpsF(42)*Fa(1063)+EpsF
     -   (45)*Fa(1064)+EpsF(57)*Fa(1068)+EpsF(86)*Fa(1069)+EpsF(87)*F
     -   a(1070)+EpsF(52)*Fa(2623)+EpsF(48)*Fa(2624)+EpsF(49)*Fa(2625
     -   )+EpsF(50)*Fa(2626)-mup1mup5*EpsF(1)*F(4460)+mup3mup5*EpsF(1
     -   4)*F(4460)-mup4mup5*EpsF(17)*F(4460)-mup3mup4*EpsF(59)*F(446
     -   0)-mup1mup3*EpsF(116)*F(4460)+mup1mup4*EpsF(120)*F(4460)+mup
     -   1mup5*EpsF(2)*F(4461)-mup3mup5*EpsF(15)*F(4461)+mup4mup5*Eps
     -   F(18)*F(4461)+mup3mup4*EpsF(62)*F(4461)-mup1mup4*EpsF(74)*F(
     -   4461)-p2mup1*EpsF(88)*F(4461)+mup1mup3*EpsF(104)*F(4461)+mup
     -   1mup5*EpsF(3)*F(4462)-mup3mup5*EpsF(16)*F(4462)+mup4mup5*Eps
     -   F(19)*F(4462)+mup3mup4*EpsF(63)*F(4462)-p2mup1*EpsF(97)*F(44
     -   62)+mup1mup3*EpsF(101)*F(4462)-mup1mup4*EpsF(112)*F(4462)))+
     -   Cpl(4)*Cpl(8)*Cpl(10)*(-2*EpsF(41)*F(2650)-4*(mup3mup5*EpsF(
     -   53)*F(2605)+mup1mup4*EpsF(123)*F(2605)+mup4mup5*EpsF(46)*F(2
     -   616)+mup3mup4*EpsF(61)*F(2616)+mup1mup5*EpsF(98)*F(2616)+mup
     -   1mup3*EpsF(126)*F(2616)+mup3mup5*EpsF(51)*F(2625)+mup1mup4*E
     -   psF(119)*F(2625)+mup4mup5*EpsF(43)*F(2627)+mup3mup4*EpsF(60)
     -   *F(2627)+mup1mup5*EpsF(91)*F(2627)+mup1mup3*EpsF(118)*F(2627
     -   )-EpsF(87)*Fa(1070)+EpsF(21)*Fa(1091)+EpsF(22)*Fa(1092)+EpsF
     -   (23)*Fa(1093)+EpsF(24)*Fa(1094)+EpsF(42)*Fa(1095)-EpsF(45)*F
     -   a(1096)+EpsF(57)*Fa(1101)+EpsF(58)*Fa(1102)+EpsF(86)*Fa(1103
     -   )+EpsF(48)*Fa(2634)+EpsF(49)*Fa(2635)+EpsF(50)*Fa(2636)+EpsF
     -   (52)*Fa(2637)+mup1mup5*EpsF(1)*F(4460)-mup3mup5*EpsF(14)*F(4
     -   460)+mup4mup5*EpsF(17)*F(4460)+mup3mup4*EpsF(59)*F(4460)+mup
     -   1mup3*EpsF(116)*F(4460)-mup1mup4*EpsF(120)*F(4460)-mup1mup5*
     -   EpsF(2)*F(4461)+mup3mup5*EpsF(15)*F(4461)-mup4mup5*EpsF(18)*
     -   F(4461)-mup3mup4*EpsF(62)*F(4461)+mup1mup4*EpsF(74)*F(4461)+
     -   p2mup1*EpsF(88)*F(4461)-mup1mup3*EpsF(104)*F(4461)-mup1mup5*
     -   EpsF(3)*F(4462)+mup3mup5*EpsF(16)*F(4462)-mup4mup5*EpsF(19)*
     -   F(4462)-mup3mup4*EpsF(63)*F(4462)+p2mup1*EpsF(97)*F(4462)-mu
     -   p1mup3*EpsF(101)*F(4462)+mup1mup4*EpsF(112)*F(4462))))+M1*(M
     -   3*(Cpl(3)*Cpl(9)*Cpl(10)*(-2*EpsF(41)*F(2635)-4*(mup1mup5*Ep
     -   sF(1)*F(2583)-mup3mup5*EpsF(14)*F(2583)+mup4mup5*EpsF(17)*F(
     -   2583)+mup3mup4*EpsF(59)*F(2583)+mup1mup3*EpsF(116)*F(2583)-m
     -   up1mup4*EpsF(120)*F(2583)-mup1mup5*EpsF(2)*F(2584)+mup3mup5*
     -   EpsF(15)*F(2584)-mup4mup5*EpsF(18)*F(2584)-mup3mup4*EpsF(62)
     -   *F(2584)+mup1mup4*EpsF(74)*F(2584)+p2mup1*EpsF(88)*F(2584)-m
     -   up1mup3*EpsF(104)*F(2584)+mup3mup5*EpsF(51)*F(2602)+mup1mup4
     -   *EpsF(119)*F(2602)+mup4mup5*EpsF(43)*F(2615)+mup3mup4*EpsF(6
     -   0)*F(2615)+mup1mup5*EpsF(91)*F(2615)+mup1mup3*EpsF(118)*F(26
     -   15)-mup1mup5*EpsF(3)*F(2628)+mup3mup5*EpsF(16)*F(2628)-mup4m
     -   up5*EpsF(19)*F(2628)-mup3mup4*EpsF(63)*F(2628)+p2mup1*EpsF(9
     -   7)*F(2628)-mup1mup3*EpsF(101)*F(2628)+mup1mup4*EpsF(112)*F(2
     -   628)+mup4mup5*EpsF(46)*F(2638)+mup3mup4*EpsF(61)*F(2638)+mup
     -   1mup5*EpsF(98)*F(2638)+mup1mup3*EpsF(126)*F(2638)+mup3mup5*E
     -   psF(53)*F(2645)+mup1mup4*EpsF(123)*F(2645)+EpsF(21)*Fa(1071)
     -   -EpsF(22)*Fa(1072)+EpsF(23)*Fa(1073)+EpsF(24)*Fa(1074)+EpsF(
     -   42)*Fa(1075)-EpsF(45)*Fa(1076)+EpsF(57)*Fa(1080)+EpsF(58)*Fa
     -   (1081)+EpsF(86)*Fa(1082)-EpsF(87)*Fa(1083)+EpsF(48)*Fa(2627)
     -   +EpsF(49)*Fa(2628)+EpsF(50)*Fa(2629)+EpsF(52)*Fa(2630)))+Cpl
     -   (4)*Cpl(5)*Cpl(8)*(-2*EpsF(41)*F(2647)-4*(-(mup1mup5*EpsF(1)
     -   *F(2583))+mup3mup5*EpsF(14)*F(2583)-mup4mup5*EpsF(17)*F(2583
     -   )-mup3mup4*EpsF(59)*F(2583)-mup1mup3*EpsF(116)*F(2583)+mup1m
     -   up4*EpsF(120)*F(2583)+mup1mup5*EpsF(2)*F(2584)-mup3mup5*EpsF
     -   (15)*F(2584)+mup4mup5*EpsF(18)*F(2584)+mup3mup4*EpsF(62)*F(2
     -   584)-mup1mup4*EpsF(74)*F(2584)-p2mup1*EpsF(88)*F(2584)+mup1m
     -   up3*EpsF(104)*F(2584)+mup4mup5*EpsF(43)*F(2602)+mup3mup4*Eps
     -   F(60)*F(2602)+mup1mup5*EpsF(91)*F(2602)+mup1mup3*EpsF(118)*F
     -   (2602)+mup3mup5*EpsF(51)*F(2615)+mup1mup4*EpsF(119)*F(2615)+
     -   mup1mup5*EpsF(3)*F(2628)-mup3mup5*EpsF(16)*F(2628)+mup4mup5*
     -   EpsF(19)*F(2628)+mup3mup4*EpsF(63)*F(2628)-p2mup1*EpsF(97)*F
     -   (2628)+mup1mup3*EpsF(101)*F(2628)-mup1mup4*EpsF(112)*F(2628)
     -   +mup3mup5*EpsF(53)*F(2638)+mup1mup4*EpsF(123)*F(2638)+mup4mu
     -   p5*EpsF(46)*F(2645)+mup3mup4*EpsF(61)*F(2645)+mup1mup5*EpsF(
     -   98)*F(2645)+mup1mup3*EpsF(126)*F(2645)+EpsF(23)*Fa(1048)+Eps
     -   F(42)*Fa(1050)+EpsF(45)*Fa(1076)+EpsF(87)*Fa(1083)+EpsF(21)*
     -   Fa(1084)+EpsF(22)*Fa(1085)+EpsF(24)*Fa(1086)+EpsF(57)*Fa(108
     -   8)+EpsF(58)*Fa(1089)+EpsF(86)*Fa(1090)+EpsF(50)*Fa(2622)+Eps
     -   F(48)*Fa(2631)+EpsF(49)*Fa(2632)+EpsF(52)*Fa(2633))))+M4*(Cp
     -   l(3)*Cpl(4)*Cpl(10)*(-2*EpsF(41)*F(2598)-4*(-(mup1mup5*EpsF(
     -   1)*F(2583))+mup3mup5*EpsF(14)*F(2583)-mup4mup5*EpsF(17)*F(25
     -   83)-mup3mup4*EpsF(59)*F(2583)-mup1mup3*EpsF(116)*F(2583)+mup
     -   1mup4*EpsF(120)*F(2583)+mup1mup5*EpsF(2)*F(2584)-mup3mup5*Ep
     -   sF(15)*F(2584)+mup4mup5*EpsF(18)*F(2584)+mup3mup4*EpsF(62)*F
     -   (2584)-mup1mup4*EpsF(74)*F(2584)-p2mup1*EpsF(88)*F(2584)+mup
     -   1mup3*EpsF(104)*F(2584)+mup4mup5*EpsF(43)*F(2602)+mup3mup4*E
     -   psF(60)*F(2602)+mup1mup5*EpsF(91)*F(2602)+mup1mup3*EpsF(118)
     -   *F(2602)+mup4mup5*EpsF(46)*F(2605)+mup3mup4*EpsF(61)*F(2605)
     -   +mup1mup5*EpsF(98)*F(2605)+mup1mup3*EpsF(126)*F(2605)+mup3mu
     -   p5*EpsF(51)*F(2615)+mup1mup4*EpsF(119)*F(2615)+mup3mup5*EpsF
     -   (53)*F(2616)+mup1mup4*EpsF(123)*F(2616)+EpsF(21)*Fa(1046)+Ep
     -   sF(22)*Fa(1047)+EpsF(23)*Fa(1048)+EpsF(24)*Fa(1049)+EpsF(42)
     -   *Fa(1050)-EpsF(45)*Fa(1051)+EpsF(57)*Fa(1056)+EpsF(58)*Fa(10
     -   57)+EpsF(86)*Fa(1058)+EpsF(87)*Fa(1059)+EpsF(48)*Fa(2620)+Ep
     -   sF(49)*Fa(2621)+EpsF(50)*Fa(2622)+EpsF(52)*Fa(2623)+mup1mup5
     -   *EpsF(3)*F(4462)-mup3mup5*EpsF(16)*F(4462)+mup4mup5*EpsF(19)
     -   *F(4462)+mup3mup4*EpsF(63)*F(4462)-p2mup1*EpsF(97)*F(4462)+m
     -   up1mup3*EpsF(101)*F(4462)-mup1mup4*EpsF(112)*F(4462)))+Cpl(5
     -   )*Cpl(8)*Cpl(9)*(-2*EpsF(41)*F(2651)+4*(-(mup1mup5*EpsF(1)*F
     -   (2583))+mup3mup5*EpsF(14)*F(2583)-mup4mup5*EpsF(17)*F(2583)-
     -   mup3mup4*EpsF(59)*F(2583)-mup1mup3*EpsF(116)*F(2583)+mup1mup
     -   4*EpsF(120)*F(2583)+mup1mup5*EpsF(2)*F(2584)-mup3mup5*EpsF(1
     -   5)*F(2584)+mup4mup5*EpsF(18)*F(2584)+mup3mup4*EpsF(62)*F(258
     -   4)-mup1mup4*EpsF(74)*F(2584)-p2mup1*EpsF(88)*F(2584)+mup1mup
     -   3*EpsF(104)*F(2584)-mup3mup5*EpsF(51)*F(2602)-mup1mup4*EpsF(
     -   119)*F(2602)-mup3mup5*EpsF(53)*F(2605)-mup1mup4*EpsF(123)*F(
     -   2605)-mup4mup5*EpsF(43)*F(2615)-mup3mup4*EpsF(60)*F(2615)-mu
     -   p1mup5*EpsF(91)*F(2615)-mup1mup3*EpsF(118)*F(2615)-mup4mup5*
     -   EpsF(46)*F(2616)-mup3mup4*EpsF(61)*F(2616)-mup1mup5*EpsF(98)
     -   *F(2616)-mup1mup3*EpsF(126)*F(2616)+EpsF(87)*Fa(1059)-EpsF(2
     -   3)*Fa(1073)-EpsF(42)*Fa(1075)-EpsF(24)*Fa(1094)-EpsF(58)*Fa(
     -   1102)-EpsF(21)*Fa(1104)+EpsF(22)*Fa(1105)-EpsF(45)*Fa(1106)-
     -   EpsF(57)*Fa(1107)-EpsF(86)*Fa(1108)-EpsF(50)*Fa(2629)-EpsF(5
     -   2)*Fa(2637)-EpsF(48)*Fa(2638)-EpsF(49)*Fa(2639)+mup1mup5*Eps
     -   F(3)*F(4462)-mup3mup5*EpsF(16)*F(4462)+mup4mup5*EpsF(19)*F(4
     -   462)+mup3mup4*EpsF(63)*F(4462)-p2mup1*EpsF(97)*F(4462)+mup1m
     -   up3*EpsF(101)*F(4462)-mup1mup4*EpsF(112)*F(4462))))))
       result(1,4) = (0,-0.5)*(-2*Cpl(5)*Cpl(8)*Cpl(9)*EpsF(50)*Fa(135
     -   8)+Cpl(3)*Cpl(4)*Cpl(5)*(-(EpsF(47)*F(2972))-4*(EpsF(65)*Fa(
     -   1434)-EpsF(66)*Fa(1435)+EpsF(113)*F(3529)*Fa(2268))-2*(EpsF(
     -   25)*Fa(1223)+EpsF(26)*Fa(1224)+EpsF(27)*Fa(1225)+EpsF(28)*Fa
     -   (1226)+EpsF(29)*Fa(1227)+EpsF(30)*Fa(1228)+EpsF(31)*Fa(1229)
     -   +EpsF(37)*Fa(1230)+EpsF(38)*Fa(1231)+EpsF(39)*Fa(1232)+EpsF(
     -   40)*Fa(1233)+EpsF(48)*Fa(1234)+EpsF(49)*Fa(1235)+EpsF(50)*Fa
     -   (1236)+EpsF(52)*Fa(1237)+EpsF(55)*Fa(1238)+EpsF(56)*Fa(1239)
     -   +EpsF(77)*Fa(1240)+EpsF(83)*Fa(1241)+EpsF(84)*Fa(1242)+EpsF(
     -   85)*Fa(1243)+EpsF(94)*Fa(1244)+EpsF(14)*Fa(1388)+EpsF(15)*Fa
     -   (1389)+EpsF(16)*Fa(1390)+EpsF(51)*Fa(1391)+EpsF(53)*Fa(1392)
     -   +EpsF(54)*Fa(1393)+EpsF(78)*Fa(1436)+EpsF(105)*Fa(1437)+EpsF
     -   (6)*Fa(1478)+EpsF(95)*Fa(1484)+EpsF(96)*Fa(1485)+EpsF(59)*Fa
     -   (2666)+EpsF(60)*Fa(2667)+EpsF(61)*Fa(2668)+EpsF(62)*Fa(2669)
     -   +EpsF(63)*Fa(2670)+EpsF(64)*Fa(2671)+EpsF(92)*Fa(2672)+EpsF(
     -   102)*Fa(2673)+EpsF(103)*Fa(2674)+EpsF(109)*Fa(2675)+EpsF(4)*
     -   Fa(2676)+EpsF(71)*Fa(2677)+EpsF(72)*Fa(2678)+EpsF(73)*Fa(267
     -   9)+EpsF(80)*Fa(2680)+EpsF(90)*Fa(2681)+EpsF(100)*Fa(2682)+Ep
     -   sF(106)*Fa(2683)+EpsF(107)*Fa(2684)+EpsF(114)*Fa(2685)+EpsF(
     -   115)*Fa(2686)+EpsF(117)*Fa(2687)+EpsF(121)*Fa(2688)+EpsF(122
     -   )*Fa(2689)+EpsF(124)*Fa(2690)+EpsF(125)*Fa(2691)+EpsF(5)*Fa(
     -   2692)+EpsF(101)*Fa(2693)+EpsF(108)*Fa(2694)+EpsF(116)*Fa(269
     -   5)+EpsF(118)*Fa(2696)+EpsF(126)*Fa(2697)+EpsF(104)*Fa(2698))
     -   )+Cpl(8)*Cpl(9)*Cpl(10)*(-(EpsF(47)*F(3278))-4*(EpsF(66)*Fa(
     -   1435)-EpsF(113)*F(3529)*Fa(2268))-2*(EpsF(25)*Fa(1366)+EpsF(
     -   26)*Fa(1367)+EpsF(27)*Fa(1368)+EpsF(28)*Fa(1369)+EpsF(29)*Fa
     -   (1370)+EpsF(30)*Fa(1371)+EpsF(31)*Fa(1372)+EpsF(37)*Fa(1373)
     -   +EpsF(38)*Fa(1374)+EpsF(39)*Fa(1375)+EpsF(40)*Fa(1376)+EpsF(
     -   48)*Fa(1377)+EpsF(49)*Fa(1378)+EpsF(50)*Fa(1379)+EpsF(52)*Fa
     -   (1380)+EpsF(55)*Fa(1381)+EpsF(56)*Fa(1382)+EpsF(77)*Fa(1383)
     -   +EpsF(83)*Fa(1384)+EpsF(84)*Fa(1385)+EpsF(85)*Fa(1386)+EpsF(
     -   94)*Fa(1387)+EpsF(14)*Fa(1422)+EpsF(15)*Fa(1423)+EpsF(16)*Fa
     -   (1424)+EpsF(51)*Fa(1425)+EpsF(53)*Fa(1426)+EpsF(54)*Fa(1427)
     -   +EpsF(65)*Fa(1463)+EpsF(78)*Fa(1464)+EpsF(105)*Fa(1465)+EpsF
     -   (6)*Fa(1531)+EpsF(95)*Fa(1537)+EpsF(96)*Fa(1538)+EpsF(59)*Fa
     -   (2814)+EpsF(60)*Fa(2815)+EpsF(61)*Fa(2816)+EpsF(62)*Fa(2817)
     -   +EpsF(63)*Fa(2818)+EpsF(64)*Fa(2819)+EpsF(92)*Fa(2820)+EpsF(
     -   102)*Fa(2821)+EpsF(103)*Fa(2822)+EpsF(109)*Fa(2823)+EpsF(4)*
     -   Fa(2824)+EpsF(71)*Fa(2825)+EpsF(72)*Fa(2826)+EpsF(73)*Fa(282
     -   7)+EpsF(80)*Fa(2828)+EpsF(90)*Fa(2829)+EpsF(100)*Fa(2830)+Ep
     -   sF(106)*Fa(2831)+EpsF(107)*Fa(2832)+EpsF(114)*Fa(2833)+EpsF(
     -   115)*Fa(2834)+EpsF(117)*Fa(2835)+EpsF(121)*Fa(2836)+EpsF(122
     -   )*Fa(2837)+EpsF(124)*Fa(2838)+EpsF(125)*Fa(2839)+EpsF(5)*Fa(
     -   2840)+EpsF(101)*Fa(2841)+EpsF(108)*Fa(2842)+EpsF(116)*Fa(284
     -   3)+EpsF(118)*Fa(2844)+EpsF(126)*Fa(2845)+EpsF(104)*Fa(2846))
     -   )+M3*M4*(Cpl(3)*Cpl(5)*Cpl(9)*(-(EpsF(47)*F(3086))-2*(EpsF(4
     -   0)*Fa(23)-EpsF(54)*Fa(172)+EpsF(96)*Fa(503)+EpsF(25)*Fa(1266
     -   )+EpsF(26)*Fa(1267)+EpsF(27)*Fa(1268)+EpsF(28)*Fa(1269)+EpsF
     -   (29)*Fa(1270)+EpsF(30)*Fa(1271)+EpsF(31)*Fa(1272)+EpsF(37)*F
     -   a(1273)+EpsF(38)*Fa(1274)+EpsF(39)*Fa(1275)+EpsF(48)*Fa(1276
     -   )+EpsF(49)*Fa(1277)+EpsF(50)*Fa(1278)+EpsF(52)*Fa(1279)+EpsF
     -   (55)*Fa(1280)+EpsF(56)*Fa(1281)+EpsF(83)*Fa(1282)+EpsF(84)*F
     -   a(1283)+EpsF(85)*Fa(1284)+EpsF(94)*Fa(1285)+EpsF(53)*Fa(1398
     -   )-EpsF(14)*Fa(1399)+EpsF(15)*Fa(1400)+EpsF(16)*Fa(1401)+EpsF
     -   (51)*Fa(1402)-EpsF(6)*Fa(1508)+EpsF(108)*Fa(2285)+EpsF(64)*F
     -   a(2434)+EpsF(106)*Fa(2451)+EpsF(61)*Fa(2701)+EpsF(126)*Fa(27
     -   21)+EpsF(59)*Fa(2723)+EpsF(60)*Fa(2724)+EpsF(62)*Fa(2725)+Ep
     -   sF(63)*Fa(2726)+EpsF(4)*Fa(2727)+EpsF(71)*Fa(2728)-EpsF(90)*
     -   Fa(2729)+EpsF(100)*Fa(2730)+EpsF(115)*Fa(2731)+EpsF(117)*Fa(
     -   2732)+EpsF(5)*Fa(2733)+EpsF(101)*Fa(2734)+EpsF(116)*Fa(2735)
     -   +EpsF(118)*Fa(2736)+EpsF(104)*Fa(2737)+EpsF(77)*Fa(1261)*F(4
     -   464)+EpsF(95)*Fa(1501)*F(4464)+EpsF(122)*Fa(2715)*F(4464))-4
     -   *(EpsF(66)*Fa(238)+EpsF(105)*Fa(242)+EpsF(65)*Fa(1443)-EpsF(
     -   103)*Fa(2249)+EpsF(109)*Fa(2258)+EpsF(114)*Fa(2276)+EpsF(125
     -   )*Fa(2281)+EpsF(73)*Fa(2438)+EpsF(92)*Fa(2704)-EpsF(72)*Fa(2
     -   708)+EpsF(121)*Fa(2714)-EpsF(78)*Fa(1444)*F(4464)-EpsF(113)*
     -   Fa(2268)*F(4464)+EpsF(107)*Fa(2353)*F(4464)+EpsF(102)*Fa(270
     -   5)*F(4464)+EpsF(80)*Fa(2709)*F(4464)+EpsF(124)*Fa(2716)*F(44
     -   64)))+Cpl(4)*Cpl(8)*Cpl(10)*(-(EpsF(47)*F(3194))-2*(EpsF(40)
     -   *Fa(62)+EpsF(54)*Fa(117)-EpsF(96)*Fa(448)+EpsF(25)*Fa(1326)+
     -   EpsF(26)*Fa(1327)+EpsF(27)*Fa(1328)+EpsF(28)*Fa(1329)+EpsF(2
     -   9)*Fa(1330)+EpsF(30)*Fa(1331)+EpsF(31)*Fa(1332)+EpsF(37)*Fa(
     -   1333)+EpsF(38)*Fa(1334)+EpsF(39)*Fa(1335)+EpsF(48)*Fa(1336)+
     -   EpsF(49)*Fa(1337)+EpsF(50)*Fa(1338)+EpsF(52)*Fa(1339)+EpsF(5
     -   5)*Fa(1340)+EpsF(56)*Fa(1341)+EpsF(83)*Fa(1342)+EpsF(84)*Fa(
     -   1343)+EpsF(85)*Fa(1344)+EpsF(94)*Fa(1345)+EpsF(14)*Fa(1413)-
     -   EpsF(15)*Fa(1414)-EpsF(16)*Fa(1415)-EpsF(51)*Fa(1416)+EpsF(5
     -   3)*Fa(1417)+EpsF(6)*Fa(1525)+EpsF(64)*Fa(2246)+EpsF(106)*Fa(
     -   2273)+EpsF(108)*Fa(2461)+EpsF(59)*Fa(2780)+EpsF(60)*Fa(2781)
     -   +EpsF(61)*Fa(2782)+EpsF(62)*Fa(2783)+EpsF(63)*Fa(2784)-EpsF(
     -   4)*Fa(2786)+EpsF(71)*Fa(2787)+EpsF(90)*Fa(2788)+EpsF(100)*Fa
     -   (2789)+EpsF(115)*Fa(2790)+EpsF(117)*Fa(2791)+EpsF(5)*Fa(2793
     -   )+EpsF(101)*Fa(2794)+EpsF(116)*Fa(2795)+EpsF(118)*Fa(2796)+E
     -   psF(126)*Fa(2797)+EpsF(104)*Fa(2798)-EpsF(77)*Fa(1261)*F(446
     -   4)-EpsF(95)*Fa(1501)*F(4464)+EpsF(122)*Fa(2752)*F(4464))+4*(
     -   EpsF(66)*Fa(238)+EpsF(105)*Fa(242)-EpsF(65)*Fa(1456)-EpsF(10
     -   3)*Fa(2249)-EpsF(73)*Fa(2254)+EpsF(114)*Fa(2276)-EpsF(109)*F
     -   a(2440)-EpsF(125)*Fa(2457)-EpsF(72)*Fa(2708)+EpsF(92)*Fa(278
     -   5)-EpsF(121)*Fa(2792)-EpsF(78)*Fa(1444)*F(4464)-EpsF(107)*Fa
     -   (2255)*F(4464)-EpsF(113)*Fa(2268)*F(4464)+EpsF(102)*Fa(2705)
     -   *F(4464)+EpsF(80)*Fa(2709)*F(4464)-EpsF(124)*Fa(2753)*F(4464
     -   ))))+M1*(M3*(Cpl(3)*Cpl(9)*Cpl(10)*(-(EpsF(47)*F(3136))+4*(E
     -   psF(66)*Fa(266)+EpsF(105)*Fa(268)-EpsF(78)*F(477)*Fa(1444)-E
     -   psF(65)*Fa(1450)-EpsF(107)*F(477)*Fa(2255)-EpsF(113)*F(477)*
     -   Fa(2268)-EpsF(103)*Fa(2347)-EpsF(109)*Fa(2356)-EpsF(114)*Fa(
     -   2370)-EpsF(125)*Fa(2375)-EpsF(73)*Fa(2402)+EpsF(102)*F(477)*
     -   Fa(2705)+EpsF(80)*F(477)*Fa(2709)+EpsF(92)*Fa(2743)-EpsF(72)
     -   *Fa(2746)-EpsF(121)*Fa(2751)-EpsF(124)*F(477)*Fa(2753))-2*(E
     -   psF(40)*Fa(45)+EpsF(54)*Fa(159)-EpsF(96)*Fa(499)-EpsF(77)*F(
     -   477)*Fa(1261)+EpsF(25)*Fa(1286)-EpsF(26)*Fa(1287)+EpsF(27)*F
     -   a(1288)+EpsF(28)*Fa(1289)+EpsF(29)*Fa(1290)+EpsF(30)*Fa(1291
     -   )+EpsF(31)*Fa(1292)+EpsF(37)*Fa(1293)+EpsF(38)*Fa(1294)+EpsF
     -   (39)*Fa(1295)+EpsF(48)*Fa(1296)+EpsF(49)*Fa(1297)+EpsF(50)*F
     -   a(1298)+EpsF(52)*Fa(1299)+EpsF(55)*Fa(1300)+EpsF(56)*Fa(1301
     -   )+EpsF(83)*Fa(1302)+EpsF(84)*Fa(1303)+EpsF(85)*Fa(1304)+EpsF
     -   (94)*Fa(1305)+EpsF(14)*Fa(1403)+EpsF(15)*Fa(1404)-EpsF(16)*F
     -   a(1405)-EpsF(51)*Fa(1406)+EpsF(53)*Fa(1407)-EpsF(95)*F(477)*
     -   Fa(1501)+EpsF(6)*Fa(1511)+EpsF(108)*Fa(2379)+EpsF(64)*Fa(239
     -   6)+EpsF(106)*Fa(2410)+EpsF(59)*Fa(2738)+EpsF(60)*Fa(2739)+Ep
     -   sF(61)*Fa(2740)+EpsF(62)*Fa(2741)+EpsF(63)*Fa(2742)+EpsF(4)*
     -   Fa(2744)+EpsF(71)*Fa(2745)+EpsF(90)*Fa(2747)+EpsF(100)*Fa(27
     -   48)+EpsF(115)*Fa(2749)+EpsF(117)*Fa(2750)+EpsF(122)*F(477)*F
     -   a(2752)+EpsF(5)*Fa(2754)+EpsF(101)*Fa(2755)+EpsF(116)*Fa(275
     -   6)+EpsF(118)*Fa(2757)+EpsF(126)*Fa(2758)+EpsF(104)*Fa(2759))
     -   )+Cpl(4)*Cpl(5)*Cpl(8)*(-(EpsF(47)*F(3176))-4*(EpsF(66)*Fa(2
     -   66)+EpsF(105)*Fa(268)-EpsF(78)*F(477)*Fa(1444)+EpsF(65)*Fa(1
     -   453)-EpsF(113)*F(477)*Fa(2268)+EpsF(73)*Fa(2352)+EpsF(107)*F
     -   (477)*Fa(2353)+EpsF(103)*Fa(2398)+EpsF(109)*Fa(2403)+EpsF(11
     -   4)*Fa(2412)+EpsF(125)*Fa(2415)+EpsF(102)*F(477)*Fa(2705)+Eps
     -   F(80)*F(477)*Fa(2709)+EpsF(124)*F(477)*Fa(2716)+EpsF(92)*Fa(
     -   2765)-EpsF(72)*Fa(2768)+EpsF(121)*Fa(2773))-2*(EpsF(40)*Fa(5
     -   3)-EpsF(54)*Fa(145)+EpsF(96)*Fa(485)+EpsF(77)*F(477)*Fa(1261
     -   )+EpsF(25)*Fa(1306)+EpsF(26)*Fa(1307)+EpsF(27)*Fa(1308)+EpsF
     -   (28)*Fa(1309)+EpsF(29)*Fa(1310)-EpsF(30)*Fa(1311)+EpsF(31)*F
     -   a(1312)+EpsF(37)*Fa(1313)+EpsF(38)*Fa(1314)+EpsF(39)*Fa(1315
     -   )+EpsF(48)*Fa(1316)+EpsF(49)*Fa(1317)+EpsF(50)*Fa(1318)+EpsF
     -   (52)*Fa(1319)+EpsF(55)*Fa(1320)+EpsF(56)*Fa(1321)+EpsF(83)*F
     -   a(1322)+EpsF(84)*Fa(1323)+EpsF(85)*Fa(1324)+EpsF(94)*Fa(1325
     -   )+EpsF(14)*Fa(1408)+EpsF(15)*Fa(1409)+EpsF(16)*Fa(1410)+EpsF
     -   (51)*Fa(1411)+EpsF(53)*Fa(1412)+EpsF(95)*F(477)*Fa(1501)-Eps
     -   F(6)*Fa(1521)+EpsF(64)*Fa(2344)+EpsF(106)*Fa(2367)+EpsF(108)
     -   *Fa(2418)+EpsF(122)*F(477)*Fa(2715)+EpsF(59)*Fa(2760)+EpsF(6
     -   0)*Fa(2761)+EpsF(61)*Fa(2762)+EpsF(62)*Fa(2763)+EpsF(63)*Fa(
     -   2764)+EpsF(4)*Fa(2766)+EpsF(71)*Fa(2767)+EpsF(90)*Fa(2769)+E
     -   psF(100)*Fa(2770)+EpsF(115)*Fa(2771)+EpsF(117)*Fa(2772)+EpsF
     -   (5)*Fa(2774)+EpsF(101)*Fa(2775)+EpsF(116)*Fa(2776)+EpsF(118)
     -   *Fa(2777)+EpsF(126)*Fa(2778)+EpsF(104)*Fa(2779))))+M4*(Cpl(3
     -   )*Cpl(4)*Cpl(10)*(-(EpsF(47)*F(3040))-2*(EpsF(40)*Fa(23)-Eps
     -   F(54)*Fa(172)+EpsF(96)*Fa(503)+EpsF(25)*Fa(1245)+EpsF(26)*Fa
     -   (1246)+EpsF(27)*Fa(1247)+EpsF(28)*Fa(1248)+EpsF(29)*Fa(1249)
     -   -EpsF(30)*Fa(1250)+EpsF(31)*Fa(1251)+EpsF(37)*Fa(1252)+EpsF(
     -   38)*Fa(1253)+EpsF(39)*Fa(1254)+EpsF(48)*Fa(1255)+EpsF(49)*Fa
     -   (1256)+EpsF(50)*Fa(1257)+EpsF(52)*Fa(1258)+EpsF(55)*Fa(1259)
     -   +EpsF(56)*Fa(1260)+EpsF(83)*Fa(1262)+EpsF(84)*Fa(1263)+EpsF(
     -   85)*Fa(1264)+EpsF(94)*Fa(1265)-EpsF(14)*Fa(1394)+EpsF(15)*Fa
     -   (1395)+EpsF(16)*Fa(1396)+EpsF(51)*Fa(1397)+EpsF(53)*Fa(1398)
     -   -EpsF(6)*Fa(1497)+EpsF(108)*Fa(2285)+EpsF(64)*Fa(2434)+EpsF(
     -   106)*Fa(2451)+EpsF(59)*Fa(2699)+EpsF(60)*Fa(2700)+EpsF(61)*F
     -   a(2701)+EpsF(62)*Fa(2702)+EpsF(63)*Fa(2703)+EpsF(4)*Fa(2706)
     -   +EpsF(71)*Fa(2707)-EpsF(90)*Fa(2710)+EpsF(100)*Fa(2711)+EpsF
     -   (115)*Fa(2712)+EpsF(117)*Fa(2713)+EpsF(5)*Fa(2717)+EpsF(101)
     -   *Fa(2718)+EpsF(116)*Fa(2719)+EpsF(118)*Fa(2720)+EpsF(126)*Fa
     -   (2721)+EpsF(104)*Fa(2722)+EpsF(77)*Fa(1261)*F(4464)+EpsF(95)
     -   *Fa(1501)*F(4464)+EpsF(122)*Fa(2715)*F(4464))-4*(EpsF(66)*Fa
     -   (238)+EpsF(105)*Fa(242)+EpsF(65)*Fa(1443)-EpsF(103)*Fa(2249)
     -   +EpsF(109)*Fa(2258)+EpsF(114)*Fa(2276)+EpsF(125)*Fa(2281)+Ep
     -   sF(73)*Fa(2438)+EpsF(92)*Fa(2704)-EpsF(72)*Fa(2708)+EpsF(121
     -   )*Fa(2714)-EpsF(78)*Fa(1444)*F(4464)-EpsF(113)*Fa(2268)*F(44
     -   64)+EpsF(107)*Fa(2353)*F(4464)+EpsF(102)*Fa(2705)*F(4464)+Ep
     -   sF(80)*Fa(2709)*F(4464)+EpsF(124)*Fa(2716)*F(4464)))+Cpl(5)*
     -   Cpl(8)*Cpl(9)*(-(EpsF(47)*F(3210))-2*(EpsF(40)*Fa(62)+EpsF(5
     -   4)*Fa(117)-EpsF(96)*Fa(448)+EpsF(25)*Fa(1346)+EpsF(26)*Fa(13
     -   47)+EpsF(27)*Fa(1348)+EpsF(28)*Fa(1349)+EpsF(29)*Fa(1350)+Ep
     -   sF(30)*Fa(1351)+EpsF(31)*Fa(1352)+EpsF(37)*Fa(1353)-EpsF(38)
     -   *Fa(1354)+EpsF(39)*Fa(1355)+EpsF(48)*Fa(1356)+EpsF(49)*Fa(13
     -   57)+EpsF(52)*Fa(1359)+EpsF(55)*Fa(1360)+EpsF(56)*Fa(1361)+Ep
     -   sF(83)*Fa(1362)+EpsF(84)*Fa(1363)+EpsF(85)*Fa(1364)+EpsF(94)
     -   *Fa(1365)+EpsF(53)*Fa(1417)+EpsF(14)*Fa(1418)-EpsF(15)*Fa(14
     -   19)-EpsF(16)*Fa(1420)-EpsF(51)*Fa(1421)+EpsF(6)*Fa(1528)+Eps
     -   F(64)*Fa(2246)+EpsF(106)*Fa(2273)+EpsF(108)*Fa(2461)+EpsF(61
     -   )*Fa(2782)+EpsF(126)*Fa(2797)+EpsF(59)*Fa(2799)+EpsF(60)*Fa(
     -   2800)+EpsF(62)*Fa(2801)+EpsF(63)*Fa(2802)-EpsF(4)*Fa(2803)+E
     -   psF(71)*Fa(2804)+EpsF(90)*Fa(2805)+EpsF(100)*Fa(2806)+EpsF(1
     -   15)*Fa(2807)+EpsF(117)*Fa(2808)+EpsF(5)*Fa(2809)+EpsF(101)*F
     -   a(2810)+EpsF(116)*Fa(2811)+EpsF(118)*Fa(2812)+EpsF(104)*Fa(2
     -   813)-EpsF(77)*Fa(1261)*F(4464)-EpsF(95)*Fa(1501)*F(4464)+Eps
     -   F(122)*Fa(2752)*F(4464))+4*(EpsF(66)*Fa(238)+EpsF(105)*Fa(24
     -   2)-EpsF(65)*Fa(1456)-EpsF(103)*Fa(2249)-EpsF(73)*Fa(2254)+Ep
     -   sF(114)*Fa(2276)-EpsF(109)*Fa(2440)-EpsF(125)*Fa(2457)-EpsF(
     -   72)*Fa(2708)+EpsF(92)*Fa(2785)-EpsF(121)*Fa(2792)-EpsF(78)*F
     -   a(1444)*F(4464)-EpsF(107)*Fa(2255)*F(4464)-EpsF(113)*Fa(2268
     -   )*F(4464)+EpsF(102)*Fa(2705)*F(4464)+EpsF(80)*Fa(2709)*F(446
     -   4)-EpsF(124)*Fa(2753)*F(4464))))))
       result(1,5) = 0
       result(1,6) = (0,-0.5)*(Cpl(3)*(-2*Cpl(9)*Cpl(10)*EpsF(33)*Fa(1
     -   683)+Cpl(4)*Cpl(5)*(-(EpsF(20)*F(3546))+4*(-(EpsF(113)*F(441
     -   2)*Fa(731))+EpsF(66)*Fa(1882)+EpsF(105)*Fa(1884)-EpsF(107)*F
     -   (4412)*Fa(2278)+EpsF(109)*Fa(2895))-2*(EpsF(10)*Fa(1622)+Eps
     -   F(11)*Fa(1623)+EpsF(12)*Fa(1624)+EpsF(13)*Fa(1625)+EpsF(21)*
     -   Fa(1626)+EpsF(22)*Fa(1627)+EpsF(23)*Fa(1628)+EpsF(24)*Fa(162
     -   9)+EpsF(32)*Fa(1634)+EpsF(33)*Fa(1635)+EpsF(34)*Fa(1636)+Eps
     -   F(35)*Fa(1637)+EpsF(65)*Fa(1881)+EpsF(78)*Fa(1883)+EpsF(28)*
     -   Fa(2847)+EpsF(29)*Fa(2848)+EpsF(30)*Fa(2849)+EpsF(31)*Fa(285
     -   0)+EpsF(25)*Fa(2851)+EpsF(26)*Fa(2852)+EpsF(27)*Fa(2853)+Eps
     -   F(39)*Fa(2854)+EpsF(40)*Fa(2855)+EpsF(77)*Fa(2856)+EpsF(14)*
     -   Fa(2857)+EpsF(15)*Fa(2858)+EpsF(16)*Fa(2859)+EpsF(17)*Fa(286
     -   0)+EpsF(18)*Fa(2861)+EpsF(19)*Fa(2862)+EpsF(43)*Fa(2863)+Eps
     -   F(44)*Fa(2864)+EpsF(46)*Fa(2865)+EpsF(51)*Fa(2866)+EpsF(53)*
     -   Fa(2867)+EpsF(54)*Fa(2868)+EpsF(4)*Fa(2869)+EpsF(6)*Fa(2870)
     -   +EpsF(7)*Fa(2871)+EpsF(8)*Fa(2872)+EpsF(9)*Fa(2873)+EpsF(69)
     -   *Fa(2874)+EpsF(70)*Fa(2875)+EpsF(72)*Fa(2876)+EpsF(73)*Fa(28
     -   77)+EpsF(75)*Fa(2878)+EpsF(76)*Fa(2879)+EpsF(79)*Fa(2880)+Ep
     -   sF(80)*Fa(2881)+EpsF(81)*Fa(2882)+EpsF(90)*Fa(2883)+EpsF(95)
     -   *Fa(2884)+EpsF(96)*Fa(2885)+EpsF(1)*Fa(2886)+EpsF(3)*Fa(2887
     -   )+EpsF(5)*Fa(2888)+EpsF(91)*Fa(2889)+EpsF(92)*Fa(2890)+EpsF(
     -   98)*Fa(2891)+EpsF(99)*Fa(2892)+EpsF(102)*Fa(2893)+EpsF(103)*
     -   Fa(2894)+EpsF(110)*Fa(2896)+EpsF(2)*Fa(2897))))-2*Cpl(4)*Cpl
     -   (5)*Cpl(8)*(EpsF(23)*Fa(1692)+EpsF(24)*Fa(1693)+EpsF(29)*Fa(
     -   3054))+M3*M4*(Cpl(3)*Cpl(5)*Cpl(9)*(-(EpsF(20)*F(3671))-2*(E
     -   psF(10)*Fa(1654)+EpsF(11)*Fa(1655)+EpsF(12)*Fa(1656)+EpsF(13
     -   )*Fa(1657)+EpsF(21)*Fa(1658)+EpsF(22)*Fa(1659)+EpsF(23)*Fa(1
     -   660)+EpsF(24)*Fa(1661)+EpsF(32)*Fa(1666)+EpsF(33)*Fa(1667)+E
     -   psF(34)*Fa(1668)+EpsF(35)*Fa(1669)+EpsF(28)*Fa(2949)+EpsF(29
     -   )*Fa(2950)+EpsF(30)*Fa(2951)+EpsF(31)*Fa(2952)+EpsF(25)*Fa(2
     -   953)+EpsF(26)*Fa(2954)+EpsF(27)*Fa(2955)+EpsF(39)*Fa(2956)+E
     -   psF(40)*Fa(2957)+EpsF(77)*Fa(2958)+EpsF(14)*Fa(2959)+EpsF(15
     -   )*Fa(2960)+EpsF(16)*Fa(2961)+EpsF(17)*Fa(2962)+EpsF(18)*Fa(2
     -   963)+EpsF(19)*Fa(2964)+EpsF(43)*Fa(2965)+EpsF(44)*Fa(2966)+E
     -   psF(46)*Fa(2967)+EpsF(51)*Fa(2968)+EpsF(53)*Fa(2969)+EpsF(54
     -   )*Fa(2970)+EpsF(4)*Fa(2971)+EpsF(6)*Fa(2972)+EpsF(7)*Fa(2973
     -   )+EpsF(8)*Fa(2974)+EpsF(9)*Fa(2975)+EpsF(69)*Fa(2976)+EpsF(7
     -   0)*Fa(2977)+EpsF(79)*Fa(2982)+EpsF(90)*Fa(2985)+EpsF(95)*Fa(
     -   2986)+EpsF(96)*Fa(2987)+EpsF(1)*Fa(2988)+EpsF(3)*Fa(2989)+Ep
     -   sF(5)*Fa(2990)+EpsF(91)*Fa(2991)+EpsF(98)*Fa(2993)+EpsF(99)*
     -   Fa(2994)+EpsF(2)*Fa(2999))-4*(EpsF(65)*Fa(1901)+EpsF(66)*Fa(
     -   1902)+EpsF(78)*Fa(1903)+EpsF(105)*Fa(1904)+EpsF(113)*Fa(2118
     -   )+EpsF(72)*Fa(2978)+EpsF(73)*Fa(2979)+EpsF(75)*Fa(2980)+EpsF
     -   (76)*Fa(2981)+EpsF(80)*Fa(2983)+EpsF(81)*Fa(2984)+EpsF(92)*F
     -   a(2992)+EpsF(102)*Fa(2995)+EpsF(103)*Fa(2996)+EpsF(109)*Fa(2
     -   997)+EpsF(110)*Fa(2998)+EpsF(107)*Fa(3000)))+Cpl(4)*Cpl(8)*C
     -   pl(10)*(-(EpsF(20)*F(3781))-4*(EpsF(65)*Fa(1911)+EpsF(66)*Fa
     -   (1912)+EpsF(78)*Fa(1913)+EpsF(105)*Fa(1914)+EpsF(113)*Fa(212
     -   5)+EpsF(72)*Fa(3030)+EpsF(73)*Fa(3031)+EpsF(75)*Fa(3032)+Eps
     -   F(76)*Fa(3033)+EpsF(80)*Fa(3035)+EpsF(81)*Fa(3036)+EpsF(92)*
     -   Fa(3044)+EpsF(102)*Fa(3047)+EpsF(103)*Fa(3048)+EpsF(109)*Fa(
     -   3049)+EpsF(110)*Fa(3050)+EpsF(107)*Fa(3052))-2*(EpsF(10)*Fa(
     -   1702)+EpsF(11)*Fa(1703)+EpsF(12)*Fa(1704)+EpsF(13)*Fa(1705)+
     -   EpsF(21)*Fa(1706)+EpsF(22)*Fa(1707)+EpsF(23)*Fa(1708)+EpsF(2
     -   4)*Fa(1709)+EpsF(32)*Fa(1714)+EpsF(33)*Fa(1715)+EpsF(34)*Fa(
     -   1716)+EpsF(35)*Fa(1717)+EpsF(39)*Fa(3008)+EpsF(40)*Fa(3009)+
     -   EpsF(77)*Fa(3010)+EpsF(43)*Fa(3017)+EpsF(44)*Fa(3018)+EpsF(4
     -   6)*Fa(3019)+EpsF(51)*Fa(3020)+EpsF(53)*Fa(3021)+EpsF(54)*Fa(
     -   3022)+EpsF(69)*Fa(3028)+EpsF(70)*Fa(3029)+EpsF(79)*Fa(3034)+
     -   EpsF(90)*Fa(3037)+EpsF(95)*Fa(3038)+EpsF(96)*Fa(3039)+EpsF(9
     -   1)*Fa(3043)+EpsF(98)*Fa(3045)+EpsF(99)*Fa(3046)+EpsF(28)*Fa(
     -   3075)+EpsF(29)*Fa(3076)+EpsF(30)*Fa(3077)+EpsF(31)*Fa(3078)+
     -   EpsF(25)*Fa(3079)+EpsF(26)*Fa(3080)+EpsF(27)*Fa(3081)+EpsF(1
     -   4)*Fa(3082)+EpsF(15)*Fa(3083)+EpsF(16)*Fa(3084)+EpsF(17)*Fa(
     -   3085)+EpsF(18)*Fa(3086)+EpsF(19)*Fa(3087)+EpsF(4)*Fa(3088)+E
     -   psF(6)*Fa(3089)+EpsF(7)*Fa(3090)+EpsF(8)*Fa(3091)+EpsF(9)*Fa
     -   (3092)+EpsF(1)*Fa(3093)+EpsF(3)*Fa(3094)+EpsF(5)*Fa(3095)+Ep
     -   sF(2)*Fa(3096))))+M1*(M3*(Cpl(3)*Cpl(9)*Cpl(10)*(-(EpsF(20)*
     -   F(3731))-2*(EpsF(10)*Fa(1670)+EpsF(11)*Fa(1671)+EpsF(12)*Fa(
     -   1672)+EpsF(13)*Fa(1673)+EpsF(21)*Fa(1674)+EpsF(22)*Fa(1675)+
     -   EpsF(23)*Fa(1676)+EpsF(24)*Fa(1677)+EpsF(32)*Fa(1682)+EpsF(3
     -   4)*Fa(1684)+EpsF(35)*Fa(1685)+EpsF(28)*Fa(3001)+EpsF(29)*Fa(
     -   3002)+EpsF(30)*Fa(3003)+EpsF(31)*Fa(3004)+EpsF(25)*Fa(3005)+
     -   EpsF(26)*Fa(3006)+EpsF(27)*Fa(3007)+EpsF(39)*Fa(3008)+EpsF(4
     -   0)*Fa(3009)+EpsF(77)*Fa(3010)+EpsF(14)*Fa(3011)+EpsF(15)*Fa(
     -   3012)+EpsF(16)*Fa(3013)+EpsF(17)*Fa(3014)+EpsF(18)*Fa(3015)+
     -   EpsF(19)*Fa(3016)+EpsF(43)*Fa(3017)+EpsF(44)*Fa(3018)+EpsF(4
     -   6)*Fa(3019)+EpsF(51)*Fa(3020)+EpsF(53)*Fa(3021)+EpsF(54)*Fa(
     -   3022)+EpsF(4)*Fa(3023)+EpsF(6)*Fa(3024)+EpsF(7)*Fa(3025)+Eps
     -   F(8)*Fa(3026)+EpsF(9)*Fa(3027)+EpsF(69)*Fa(3028)+EpsF(70)*Fa
     -   (3029)+EpsF(79)*Fa(3034)+EpsF(90)*Fa(3037)+EpsF(95)*Fa(3038)
     -   +EpsF(96)*Fa(3039)+EpsF(1)*Fa(3040)+EpsF(3)*Fa(3041)+EpsF(5)
     -   *Fa(3042)+EpsF(91)*Fa(3043)+EpsF(98)*Fa(3045)+EpsF(99)*Fa(30
     -   46)+EpsF(2)*Fa(3051))-4*(EpsF(65)*Fa(1911)+EpsF(66)*Fa(1912)
     -   +EpsF(78)*Fa(1913)+EpsF(105)*Fa(1914)+EpsF(113)*Fa(2125)+Eps
     -   F(72)*Fa(3030)+EpsF(73)*Fa(3031)+EpsF(75)*Fa(3032)+EpsF(76)*
     -   Fa(3033)+EpsF(80)*Fa(3035)+EpsF(81)*Fa(3036)+EpsF(92)*Fa(304
     -   4)+EpsF(102)*Fa(3047)+EpsF(103)*Fa(3048)+EpsF(109)*Fa(3049)+
     -   EpsF(110)*Fa(3050)+EpsF(107)*Fa(3052)))+Cpl(4)*Cpl(5)*Cpl(8)
     -   *(-(EpsF(20)*F(3765))-4*(EpsF(65)*Fa(1901)+EpsF(66)*Fa(1902)
     -   +EpsF(78)*Fa(1903)+EpsF(105)*Fa(1904)+EpsF(113)*Fa(2118)+Eps
     -   F(72)*Fa(2978)+EpsF(73)*Fa(2979)+EpsF(75)*Fa(2980)+EpsF(76)*
     -   Fa(2981)+EpsF(80)*Fa(2983)+EpsF(81)*Fa(2984)+EpsF(92)*Fa(299
     -   2)+EpsF(102)*Fa(2995)+EpsF(103)*Fa(2996)+EpsF(109)*Fa(2997)+
     -   EpsF(110)*Fa(2998)+EpsF(107)*Fa(3000))-2*(EpsF(10)*Fa(1686)+
     -   EpsF(11)*Fa(1687)+EpsF(12)*Fa(1688)+EpsF(13)*Fa(1689)+EpsF(2
     -   1)*Fa(1690)+EpsF(22)*Fa(1691)+EpsF(32)*Fa(1698)+EpsF(33)*Fa(
     -   1699)+EpsF(34)*Fa(1700)+EpsF(35)*Fa(1701)+EpsF(39)*Fa(2956)+
     -   EpsF(40)*Fa(2957)+EpsF(77)*Fa(2958)+EpsF(43)*Fa(2965)+EpsF(4
     -   4)*Fa(2966)+EpsF(46)*Fa(2967)+EpsF(51)*Fa(2968)+EpsF(53)*Fa(
     -   2969)+EpsF(54)*Fa(2970)+EpsF(69)*Fa(2976)+EpsF(70)*Fa(2977)+
     -   EpsF(79)*Fa(2982)+EpsF(90)*Fa(2985)+EpsF(95)*Fa(2986)+EpsF(9
     -   6)*Fa(2987)+EpsF(91)*Fa(2991)+EpsF(98)*Fa(2993)+EpsF(99)*Fa(
     -   2994)+EpsF(28)*Fa(3053)+EpsF(30)*Fa(3055)+EpsF(31)*Fa(3056)+
     -   EpsF(25)*Fa(3057)+EpsF(26)*Fa(3058)+EpsF(27)*Fa(3059)+EpsF(1
     -   4)*Fa(3060)+EpsF(15)*Fa(3061)+EpsF(16)*Fa(3062)+EpsF(17)*Fa(
     -   3063)+EpsF(18)*Fa(3064)+EpsF(19)*Fa(3065)+EpsF(4)*Fa(3066)+E
     -   psF(6)*Fa(3067)+EpsF(7)*Fa(3068)+EpsF(8)*Fa(3069)+EpsF(9)*Fa
     -   (3070)+EpsF(1)*Fa(3071)+EpsF(3)*Fa(3072)+EpsF(5)*Fa(3073)+Ep
     -   sF(2)*Fa(3074))))+M4*(Cpl(3)*Cpl(4)*Cpl(10)*(-(EpsF(20)*F(36
     -   10))-4*(-(EpsF(113)*F(4444)*Fa(731))+EpsF(65)*Fa(1891)+EpsF(
     -   66)*Fa(1892)+EpsF(78)*Fa(1893)+EpsF(105)*Fa(1894)+EpsF(107)*
     -   F(4444)*Fa(2372)+EpsF(72)*Fa(2927)-EpsF(73)*Fa(2928)+EpsF(75
     -   )*Fa(2929)+EpsF(76)*Fa(2930)+EpsF(80)*Fa(2932)+EpsF(81)*Fa(2
     -   933)+EpsF(92)*Fa(2941)+EpsF(102)*Fa(2944)+EpsF(103)*Fa(2945)
     -   +EpsF(109)*Fa(2946)+EpsF(110)*Fa(2947))-2*(EpsF(10)*Fa(1638)
     -   +EpsF(11)*Fa(1639)+EpsF(12)*Fa(1640)+EpsF(13)*Fa(1641)+EpsF(
     -   21)*Fa(1642)+EpsF(22)*Fa(1643)+EpsF(23)*Fa(1644)+EpsF(24)*Fa
     -   (1645)+EpsF(32)*Fa(1650)+EpsF(33)*Fa(1651)+EpsF(34)*Fa(1652)
     -   +EpsF(35)*Fa(1653)+EpsF(28)*Fa(2898)+EpsF(29)*Fa(2899)+EpsF(
     -   30)*Fa(2900)+EpsF(31)*Fa(2901)+EpsF(25)*Fa(2902)+EpsF(26)*Fa
     -   (2903)+EpsF(27)*Fa(2904)+EpsF(39)*Fa(2905)+EpsF(40)*Fa(2906)
     -   +EpsF(77)*Fa(2907)+EpsF(14)*Fa(2908)+EpsF(15)*Fa(2909)+EpsF(
     -   16)*Fa(2910)+EpsF(17)*Fa(2911)+EpsF(18)*Fa(2912)+EpsF(19)*Fa
     -   (2913)+EpsF(43)*Fa(2914)+EpsF(44)*Fa(2915)+EpsF(46)*Fa(2916)
     -   +EpsF(51)*Fa(2917)+EpsF(53)*Fa(2918)+EpsF(54)*Fa(2919)+EpsF(
     -   4)*Fa(2920)+EpsF(6)*Fa(2921)+EpsF(7)*Fa(2922)+EpsF(8)*Fa(292
     -   3)+EpsF(9)*Fa(2924)+EpsF(69)*Fa(2925)+EpsF(70)*Fa(2926)+EpsF
     -   (79)*Fa(2931)+EpsF(90)*Fa(2934)+EpsF(95)*Fa(2935)+EpsF(96)*F
     -   a(2936)+EpsF(1)*Fa(2937)+EpsF(3)*Fa(2938)+EpsF(5)*Fa(2939)+E
     -   psF(91)*Fa(2940)+EpsF(98)*Fa(2942)+EpsF(99)*Fa(2943)+EpsF(2)
     -   *Fa(2948)))+Cpl(5)*Cpl(8)*Cpl(9)*(-(EpsF(20)*F(3802))+4*(-(E
     -   psF(113)*F(4444)*Fa(731))+EpsF(66)*Fa(1892)+EpsF(105)*Fa(189
     -   4)-EpsF(65)*Fa(1919)-EpsF(78)*Fa(1920)-EpsF(107)*F(4444)*Fa(
     -   2278)-EpsF(73)*Fa(2928)+EpsF(109)*Fa(2946)-EpsF(72)*Fa(3126)
     -   -EpsF(75)*Fa(3127)-EpsF(76)*Fa(3128)-EpsF(80)*Fa(3130)-EpsF(
     -   81)*Fa(3131)-EpsF(92)*Fa(3139)+EpsF(102)*Fa(3142)-EpsF(103)*
     -   Fa(3143)-EpsF(110)*Fa(3144))-2*(EpsF(10)*Fa(1718)+EpsF(11)*F
     -   a(1719)+EpsF(12)*Fa(1720)+EpsF(13)*Fa(1721)+EpsF(21)*Fa(1722
     -   )+EpsF(22)*Fa(1723)+EpsF(23)*Fa(1724)+EpsF(24)*Fa(1725)+EpsF
     -   (32)*Fa(1730)+EpsF(33)*Fa(1731)+EpsF(34)*Fa(1732)+EpsF(35)*F
     -   a(1733)+EpsF(28)*Fa(3097)+EpsF(29)*Fa(3098)+EpsF(30)*Fa(3099
     -   )+EpsF(31)*Fa(3100)+EpsF(25)*Fa(3101)+EpsF(26)*Fa(3102)+EpsF
     -   (27)*Fa(3103)+EpsF(39)*Fa(3104)+EpsF(40)*Fa(3105)+EpsF(77)*F
     -   a(3106)+EpsF(14)*Fa(3107)+EpsF(15)*Fa(3108)-EpsF(16)*Fa(3109
     -   )+EpsF(17)*Fa(3110)+EpsF(18)*Fa(3111)+EpsF(19)*Fa(3112)+EpsF
     -   (43)*Fa(3113)+EpsF(44)*Fa(3114)+EpsF(46)*Fa(3115)+EpsF(51)*F
     -   a(3116)+EpsF(53)*Fa(3117)+EpsF(54)*Fa(3118)+EpsF(4)*Fa(3119)
     -   +EpsF(6)*Fa(3120)+EpsF(7)*Fa(3121)+EpsF(8)*Fa(3122)+EpsF(9)*
     -   Fa(3123)+EpsF(69)*Fa(3124)+EpsF(70)*Fa(3125)+EpsF(79)*Fa(312
     -   9)+EpsF(90)*Fa(3132)+EpsF(95)*Fa(3133)+EpsF(96)*Fa(3134)+Eps
     -   F(1)*Fa(3135)+EpsF(3)*Fa(3136)+EpsF(5)*Fa(3137)+EpsF(91)*Fa(
     -   3138)+EpsF(98)*Fa(3140)+EpsF(99)*Fa(3141)+EpsF(2)*Fa(3145)))
     -   ))+Cpl(8)*Cpl(9)*Cpl(10)*(-(EpsF(20)*F(3849))-4*(-(EpsF(113)
     -   *F(4412)*Fa(731))+EpsF(66)*Fa(1882)+EpsF(107)*F(4412)*Fa(237
     -   2))-2*(EpsF(10)*Fa(1734)+EpsF(11)*Fa(1735)+EpsF(12)*Fa(1736)
     -   +EpsF(13)*Fa(1737)+EpsF(21)*Fa(1738)+EpsF(22)*Fa(1739)+EpsF(
     -   23)*Fa(1740)+EpsF(24)*Fa(1741)+EpsF(32)*Fa(1746)+EpsF(33)*Fa
     -   (1747)+EpsF(34)*Fa(1748)+EpsF(35)*Fa(1749)+EpsF(65)*Fa(1925)
     -   +EpsF(78)*Fa(1926)+EpsF(105)*Fa(1927)+EpsF(28)*Fa(3146)+EpsF
     -   (29)*Fa(3147)+EpsF(30)*Fa(3148)+EpsF(31)*Fa(3149)+EpsF(25)*F
     -   a(3150)+EpsF(26)*Fa(3151)+EpsF(27)*Fa(3152)+EpsF(39)*Fa(3153
     -   )+EpsF(40)*Fa(3154)+EpsF(77)*Fa(3155)+EpsF(14)*Fa(3156)+EpsF
     -   (15)*Fa(3157)+EpsF(16)*Fa(3158)+EpsF(17)*Fa(3159)+EpsF(18)*F
     -   a(3160)+EpsF(19)*Fa(3161)+EpsF(43)*Fa(3162)+EpsF(44)*Fa(3163
     -   )+EpsF(46)*Fa(3164)+EpsF(51)*Fa(3165)+EpsF(53)*Fa(3166)+EpsF
     -   (54)*Fa(3167)+EpsF(4)*Fa(3168)+EpsF(6)*Fa(3169)+EpsF(7)*Fa(3
     -   170)+EpsF(8)*Fa(3171)+EpsF(9)*Fa(3172)+EpsF(69)*Fa(3173)+Eps
     -   F(70)*Fa(3174)+EpsF(72)*Fa(3175)+EpsF(73)*Fa(3176)+EpsF(75)*
     -   Fa(3177)+EpsF(76)*Fa(3178)+EpsF(79)*Fa(3179)+EpsF(80)*Fa(318
     -   0)+EpsF(81)*Fa(3181)+EpsF(90)*Fa(3182)+EpsF(95)*Fa(3183)+Eps
     -   F(96)*Fa(3184)+EpsF(1)*Fa(3185)+EpsF(3)*Fa(3186)+EpsF(5)*Fa(
     -   3187)+EpsF(91)*Fa(3188)+EpsF(92)*Fa(3189)+EpsF(98)*Fa(3190)+
     -   EpsF(99)*Fa(3191)+EpsF(102)*Fa(3192)+EpsF(103)*Fa(3193)+EpsF
     -   (109)*Fa(3194)+EpsF(110)*Fa(3195)+EpsF(2)*Fa(3196))))
       endif
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
