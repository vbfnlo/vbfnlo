       subroutine FermionBoxOddCouplings_QUAD(za1,za2,za3,za4,zb1,zb2,zb3,z
     -   b4,m,p1,p2,p3,p4,mup1,mup2,mup3,mup4,musq,comp,result,result
     -   n)
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************

c ************************************************************************************
c determine the fermion looop
c Note: No couplings are assigned
c Note: To make it shorter in the promgram: mu_p2,...->mup2,... 
c Notation of External momenta: p1+p2+p3+p4=0 
c mu_p1,mu_p2,mu_p3,mu_p4 should be think as external current 
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       IMPLICIT NONE

#include "VBFNLO/utilities/process.inc"

      Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3)
       Real*8   p1sq, p1p2, p1p3, p1p4, p2sq, p2p3, p2p4, p3sq, 
     -          p3p4, p4sq
       Real*8   s12, s13, s14, s23, s24, s34
       Complex*16   p1mup1, p1mup2, p1mup3, p1mup4, p2mup1, p2mup2, 
     -          p2mup3, p2mup4, p3mup1, p3mup2, p3mup3, p3mup4, 
     -          p4mup1, p4mup2, p4mup3, p4mup4
       Complex*16   mup1mup2, mup1mup3, mup1mup4, mup2mup3, mup2mup4, 
     -          mup3mup4
       Real*8 dotrr
       Complex*16 B0fin,C0fin,D0fin
       Complex*16 A0fin,B0t1,C0fin1M,D04
       EXTERNAL dotrr,A0fin,B0t1,C0fin1M,D04,B0fin,C0fin,D0fin

      Real*8 B012R,B013R,B014R,B023R,B024R,B034R
      Real*8 B012I,B013I,B014I,B023I,B024I,B034I
      Real*8 B12R(2,2),B13R(2,2),B14R(2,2),B23R(2,2),B24R(2,2),B34R(2,2)
      Real*8 B12I(2,2),B13I(2,2),B14I(2,2),B23I(2,2),B24I(2,2),B34I(2,2)

       Real*8 C0123R,C0124R,C0134R,C0234R
       Real*8 C0123I,C0124I,C0134I,C0234I
       Real*8 C123R(6,3),C124R(6,3),C134R(6,3),C234R(6,3)
       Real*8 C123I(6,3),C124I(6,3),C134I(6,3),C234I(6,3)

        Real*8 D01234R
       Real*8 D01234I
       Real*8 D1234R(22,4)
        Real*8 D1234I(22,4)
      Complex*16 F(242)
      Real*8 FI(242),FR(242),P(321)
      Real*8 za1,za2,za3,za4,zb1,zb2,zb3,zb4
      Complex*16 mup1(0:3),mup2(0:3),mup3(0:3),mup4(0:3)
      Complex*16 dotrc,dotcc,result,resultn(5)
       Real*8 musq,m,msq
      EXTERNAL   dotrc,dotcc
      Integer comp
       common/FerBoxOddFunctions/F
       common/PFerBoxOddFunctions/P
      Save/FerBoxOddFunctions/ 
      Save/PFerBoxOddFunctions/ 
      Complex*16 Fa(35)
      Common/musqInv/musqcp
      real*8 musqcp
       logical printmy
      common/mprint/printmy

      COMMON/DVALFunctionsGSim/ D01234R,D1234R,D01234I,
     &                D1234I
       COMMON/CVALFunctionsGSim/ C0123R,C123R,C0123I,C123I, 
     &             C0124R,C124R,C0124I,C124I, 
     &             C0134R,C134R,C0134I,C134I, 
     -             C0234R,C234R,C0234I,C234I
       COMMON/BVALFunctionsGSim/ 
     &           B012R,B012I,B013R,B013I, 
     &           B014R,B014I,B023R,B023I,
     -           B024R,B024I,B034R,B034I,
     -           B12R,B12I, B13R,B13I,
     &           B14R,B14I,B23R,B23I,
     -           B24R,B24I, B34R, B34I

ccccccccccccccccccccc For QUAD precission ccccccccccccccccccccccccc
       Real*16 musq_QUAD, dotrr_QUAD 
       External dotrr_QUAD 

       Real*16   p1sq_QUAD, p1p2_QUAD, p1p3_QUAD, p1p4_QUAD
       Real*16   p2sq_QUAD, p2p3_QUAD, p2p4_QUAD
       Real*16   p3sq_QUAD, p3p4_QUAD
       Real*16   p4sq_QUAD


       Real*16   s12_QUAD, s13_QUAD, s14_QUAD
       Real*16   s23_QUAD, s24_QUAD
       Real*16   s34_QUAD


       Real*16  p1_QUAD(0:3),p2_QUAD(0:3),p3_QUAD(0:3),p4_QUAD(0:3)
       Real*16 ptemp(0:3,4)
       Real*8 psq(4)
       Integer i,jtemp
       Real*16   M_QUAD,msq_QUAD
       integer j1,j2,j3,j4

      double precision     hmass, wmass, zmass, amass,
     &                     tmass, bmass, lmass, cmass
      common /COUPL_MASS/  hmass, wmass, zmass, amass,
     &                     tmass, bmass, lmass, cmass
       integer numberZ, numberw
       double precision vmass
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
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
c************************************************************************************
c************************************************************************************
c************************************************************************************

      musqcp=musq


       m_quad=m
       msq_QUAD=dotrr(m_quad,m_quad)
       msq=msq_quad
c`       msq=m*m


       do i=0,3 
       ptemp(i,1)=p1(i)
       ptemp(i,2)=p2(i)
       ptemp(i,3)=p3(i)
       ptemp(i,4)=p4(i)
       enddo


       psq(1)=abs(dotrr(p1,p1))
       psq(2)=abs(dotrr(p2,p2))
       psq(3)=abs(dotrr(p3,p3))
       psq(4)=abs(dotrr(p4,p4))
 

       if (narrowwidth.eq.1) then
c find on-shell bosons
         numberZ=0
         numberW=0
cccccccc Identify the Z cccccccccccccccccccccccc
         j1=0
         do i=1,4
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

         jtemp=j1
         if (ptemp(0,jtemp).ge.0d0) then
           ptemp(0,jtemp)=SQRT(VMASS*VMASS
     &                     +ptemp(1,jtemp)*ptemp(1,jtemp)+
     -                      ptemp(2,jtemp)*ptemp(2,jtemp)+
     -                      ptemp(3,jtemp)*ptemp(3,jtemp))
         else
           ptemp(0,jtemp)=-SQRT(VMASS*VMASS
     &                     +ptemp(1,jtemp)*ptemp(1,jtemp)+
     -                      ptemp(2,jtemp)*ptemp(2,jtemp)+
     -                      ptemp(3,jtemp)*ptemp(3,jtemp))
         endif

c second Z
         jtemp=j2
         if (ptemp(0,jtemp).ge.0d0) then
           ptemp(0,jtemp)=SQRT(VMASS*VMASS
     &                     +ptemp(1,jtemp)*ptemp(1,jtemp)+
     -                      ptemp(2,jtemp)*ptemp(2,jtemp)+
     -                      ptemp(3,jtemp)*ptemp(3,jtemp))
         else
           ptemp(0,jtemp)=-SQRT(VMASS*VMASS
     &                     +ptemp(1,jtemp)*ptemp(1,jtemp)+
     -                      ptemp(2,jtemp)*ptemp(2,jtemp)+
     -                      ptemp(3,jtemp)*ptemp(3,jtemp))
         endif


ccccccccccccccccccccccccccccccccccccccccccccccccccc
c  identify GLUON ONSHELL -- take the one with smaller |q^2|
         j3=0
         do jtemp=1,4
           if (jtemp.eq.j1) cycle
           if (jtemp.eq.j2) cycle
           if (j3.eq.0) then
             j3=jtemp
           else
             if (abs(psq(j3)) .gt. abs(psq(jtemp))) then
               j4 = j3
               j3 = jtemp
             else
               j4 = jtemp
             endif
           endif
         enddo
         if (ptemp(0,j3).ge.0d0) then
           ptemp(0,j3)=SQRT(
     &                     +ptemp(1,j3)*ptemp(1,j3)+
     -                      ptemp(2,j3)*ptemp(2,j3)+
     -                      ptemp(3,j3)*ptemp(3,j3))
         else
           ptemp(0,j3)=-SQRT(
     &                     +ptemp(1,j3)*ptemp(1,j3)+
     -                      ptemp(2,j3)*ptemp(2,j3)+
     -                      ptemp(3,j3)*ptemp(3,j3))
         endif
c Off shell gluon
         do i=0,3
            ptemp(i,j4)=-(ptemp(i,j1)+ptemp(i,j2)+ptemp(i,j3))
         enddo

       else ! off-shell or no two W/Z found

c find particle with largest |q^2|
c put everything below m=10^-3 massless
         j1=0
         j2=0
         j3=0
         j4=0
         do jtemp=1,4
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
           else
             if (abs(psq(j3)) .gt. abs(psq(jtemp))) then
               j4 = j3
               j3 = jtemp
             else
               j4 = jtemp
             endif
           endif
         enddo

cc off-shell 
         do i=0,3
           ptemp(i,j4)=-ptemp(i,j1)-ptemp(i,j2)-ptemp(i,j3)
         enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cMK       else
cMK          print*,"D:Wrong momenta recusntraciont"
cMK          print*, "number of Z", NUMBERZ
cMK          do i=1,4
cMK          print*, "psq particle ", i,"is",  psq(i)
cMK          enddo
cMKc          stop
cMK
cMK      do i=0,3
cMK      print*, "conv",ptemp(i,1)+ptemp(i,2)+ptemp(i,3)+ptemp(i,4)
cMK            write(*,'(A3,I1,A2,E25.16,E25.16)') "p1(",i,")=",p1(i)-ptemp(i,1)
cMK            write(*,'(A3,I1,A2,E25.16,E25.16)') "p2(",i,")=",p2(i)-ptemp(i,2)
cMK            write(*,'(A3,I1,A2,E25.16,E25.16)') "p3(",i,")=",p3(i)-ptemp(i,3)
cMK            write(*,'(A3,I1,A2,E25.16,E25.16)') "p4(",i,")=",p4(i)-ptemp(i,4)
cMK             enddo
cMK
       endif


c      do i=0,3
c      if(abs(ptemp(i,1)+ptemp(i,2)+ptemp(i,3)+ptemp(i,4)).gt.1d-14) then
c      print*, "conv",ptemp(i,1)+ptemp(i,2)+ptemp(i,3)+ptemp(i,4)
!       stop
c      endif
c      enddo


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
      p1_QUAD(i)=ptemp(i,1)
      p2_QUAD(i)=ptemp(i,2)
      p3_QUAD(i)=ptemp(i,3)
      p4_QUAD(i)=ptemp(i,4)
  

      enddo

       p1sq_QUAD = dotrr_QUAD(p1_QUAD,p1_QUAD)
       p1p2_QUAD = dotrr_QUAD(p1_QUAD,p2_QUAD)
       p1p3_QUAD = dotrr_QUAD(p1_QUAD,p3_QUAD)
       p1p4_QUAD = dotrr_QUAD(p1_QUAD,p4_QUAD)
       p2sq_QUAD = dotrr_QUAD(p2_QUAD,p2_QUAD)
       p2p3_QUAD = dotrr_QUAD(p2_QUAD,p3_QUAD)
       p2p4_QUAD = dotrr_QUAD(p2_QUAD,p4_QUAD)
       p3sq_QUAD = dotrr_QUAD(p3_QUAD,p3_QUAD)
       p3p4_QUAD = dotrr_QUAD(p3_QUAD,p4_QUAD)
       p4sq_QUAD = dotrr_QUAD(p4_QUAD,p4_QUAD)
       s12_QUAD = (p1sq_QUAD +p2sq_QUAD+ 2*p1p2_QUAD) 
       s13_QUAD = (p1sq_QUAD +p3sq_QUAD+ 2*p1p3_QUAD) 
       s14_QUAD = (p1sq_QUAD +p4sq_QUAD+ 2*p1p4_QUAD) 
       s23_QUAD = (p2sq_QUAD +p3sq_QUAD+ 2*p2p3_QUAD) 
       s24_QUAD = (p2sq_QUAD +p4sq_QUAD+ 2*p2p4_QUAD) 
       s34_QUAD = (p3sq_QUAD +p4sq_QUAD+ 2*p3p4_QUAD)


       p1sq=p1sq_QUAD 
       p1p2=p1p2_QUAD
       p1p3=p1p3_QUAD
       p1p4=p1p4_QUAD
       p2sq=p2sq_QUAD 
       p2p3=p2p3_QUAD
       p2p4=p2p4_QUAD
       p3sq=p3sq_QUAD 
       p3p4=p3p4_QUAD
       p4sq=p4sq_QUAD 
       s12=s12_QUAD 
       s13=s13_QUAD 
       s14=s14_QUAD 
       s23=s23_QUAD 
       s24=s24_QUAD 
       s34=s34_QUAD


        do i=0,3
         p1(i)=p1_QUAD(i)
         p2(i)=p2_QUAD(i)
         p3(i)=p3_QUAD(i)
         p4(i)=p4_QUAD(i)
           enddo


          If(abs(musq-s12).lt.10d-7) then
          musq_QUAD=s12_QUAD
          elseIf(abs(musq-s23).lt.10d-7) then
          musq_QUAD=s23_QUAD
          else
          musq_QUAD=musq
          endif


cfc       msq=m*m
cfc       p1sq = dotrr(p1,p1)
cfc       p1p2 = dotrr(p1,p2)
cfc       p1p3 = dotrr(p1,p3)
cfc       p1p4 = dotrr(p1,p4)
cfc       p2sq = dotrr(p2,p2)
cfc       p2p3 = dotrr(p2,p3)
cfc       p2p4 = dotrr(p2,p4)
cfc       p3sq = dotrr(p3,p3)
cfc       p3p4 = dotrr(p3,p4)
cfc       p4sq = dotrr(p4,p4)
cfc       s12 = (p1sq +p2sq+ 2*p1p2) 
cfc       s13 = (p1sq +p3sq+ 2*p1p3) 
cfc       s14 = (p1sq +p4sq+ 2*p1p4) 
cfc       s23 = (p2sq +p3sq+ 2*p2p3) 
cfc       s24 = (p2sq +p4sq+ 2*p2p4) 
cfc       s34 = (p3sq +p4sq+ 2*p3p4) 
c       Write(*,'(a5,F20.10)')," p1sq ", p1sq 
c       Write(*,'(a5,F20.10)')," p1p2 ", p1p2
c       Write(*,'(a5,F20.10)')," p1p3 ", p1p3
c       Write(*,'(a5,F20.10)')," p1p4 ", p1p4
c       Write(*,'(a5,F20.10)')," p2sq ", p2sq 
c       Write(*,'(a5,F20.10)')," p2p3 ", p2p3
c       Write(*,'(a5,F20.10)')," p2p4 ", p2p4
c       Write(*,'(a5,F20.10)')," p3sq ", p3sq 
c       Write(*,'(a5,F20.10)')," p3p4 ", p3p4
c       Write(*,'(a5,F20.10)')," p4sq ", p4sq 
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
        if (comp.gt.0) then    
   
       call XDget_1M_QUAD(M_QUAD,p1sq_QUAD,p2sq_QUAD,p3sq_QUAD,
     &   p1p2_QUAD,p1p3_QUAD,p2p3_QUAD,musq_QUAD)


         endif
c       Definition of the F,P functions:Independent of the currents    
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       P(1) = p1sq*za1*za2*za3+p3sq*za1*za2*za3-s12*za1*za2*za3-s23*za
     -   1*za2*za3-p4sq*za2*za3*zb1+p4sq*za3*zb1*zb2+p4sq*za1*za2*zb3
     -   -p4sq*za1*zb2*zb3+p1sq*zb1*zb2*zb3+p3sq*zb1*zb2*zb3-s12*zb1*
     -   zb2*zb3-s23*zb1*zb2*zb3
       P(2) = za4-zb4
       P(3) = za1*za2-zb1*zb2
       P(4) = za3-zb3
       P(5) = za1-zb1
       P(6) = za2*za3-zb2*zb3
       P(7) = -(s12*za1*za2*za3)-s23*za1*za2*za3-p4sq*za2*za3*zb1+s12*
     -   za2*za3*zb1+p4sq*za3*zb1*zb2+s23*za3*zb1*zb2+p4sq*za1*za2*zb
     -   3+s23*za1*za2*zb3-p4sq*za1*zb2*zb3+s12*za1*zb2*zb3-s12*zb1*z
     -   b2*zb3-s23*zb1*zb2*zb3+p1sq*P(3)*P(4)+p3sq*P(5)*P(6)
       P(8) = za2-zb2
       P(9) = za3*zb1+za1*zb3
       P(10) = za4+zb4
       P(11) = za1*za2*za3*za4-zb1*zb2*zb3*zb4
       P(12) = -(za1*za2*za3*za4)+zb1*zb2*zb3*zb4
       P(13) = M**2*P(8)*P(9)*P(10)+s12*P(11)+p1sq*P(12)
       P(14) = 2*zb3*zb4+za3*P(2)
       P(15) = -za4+zb4
       P(16) = 2*za3*za4+zb3*P(15)
       P(17) = zb1*P(14)+za1*P(16)
       P(18) = p2sq*P(11)-M**2*P(8)*P(17)
       P(19) = za1*za4-zb1*zb4
       P(20) = p1sq-p2sq-s12
       P(21) = -2*M**2*P(4)*P(8)*P(19)-P(11)*P(20)
       P(22) = p2sq+p3sq-s23
       P(23) = -2*M**2*P(4)*P(8)*P(19)+P(11)*P(22)
       P(24) = p1sq*za1*za2*za3-s23*za1*za2*za3-p4sq*za2*za3*zb1-p4sq*
     -   za1*zb2*zb3+p1sq*zb1*zb2*zb3-s23*zb1*zb2*zb3
       P(25) = za1+zb1
       P(26) = p1sq*P(12)+M**2*P(6)*P(10)*P(25)
       P(27) = s23*P(11)-M**2*P(6)*P(10)*P(25)
       P(28) = 2*M**2*P(4)*P(8)*P(19)-P(11)*P(22)
       P(29) = -(p2sq*za1*za2*za3)-s12*za2*za3*zb1+p1sq*za3*zb1*zb2+p1
     -   sq*za1*za2*zb3-s12*za1*zb2*zb3-p2sq*zb1*zb2*zb3+2*M**2*P(4)*
     -   P(8)*P(25)
       P(30) = za3-2*zb3
       P(31) = zb2*zb3+za2*P(30)
       P(32) = -2*za3+zb3
       P(33) = za2*za3+zb2*P(32)
       P(34) = za1*P(31)+zb1*P(33)
       P(35) = za2*za3+zb2*zb3
       P(36) = za1*za2*za3+zb1*zb2*zb3
       P(37) = za3*zb1-za1*zb3
       P(38) = za3+zb3
       P(39) = -2*zb2*zb3+za2*P(38)
       P(40) = -2*za2*za3+zb2*P(38)
       P(41) = za1*P(39)+zb1*P(40)
       P(42) = 2*zb2*zb3+za2*P(4)
       P(43) = -za3+zb3
       P(44) = 2*za2*za3+zb2*P(43)
       P(45) = za1*P(42)+zb1*P(44)
       P(46) = za2*za3*zb1+za1*zb2*zb3
       P(47) = zb1*zb2*P(30)+za1*za2*P(32)
       P(48) = 2*za2*za3-za2*zb3+zb2*zb3
       P(49) = za2*za3-za3*zb2+2*zb2*zb3
       P(50) = za1*P(48)+zb1*P(49)
       P(51) = za3-3*zb3
       P(52) = -3*za3+zb3
       P(53) = zb1*zb2*P(51)+za1*za2*P(52)
       P(54) = -(za2*za3*zb1)-za3*zb1*zb2-za1*za2*zb3-za1*zb2*zb3+2*P(
     -   36)
       P(55) = za1*za2*za4+zb1*zb2*zb4
       P(56) = -(p4sq*za1*za2*za4)+s12*za1*za2*za4+p4sq*za1*za4*zb2-p3
     -   sq*za4*zb1*zb2+s23*za4*zb1*zb2-p3sq*za1*za2*zb4+s23*za1*za2*
     -   zb4+p4sq*za2*zb1*zb4-p4sq*zb1*zb2*zb4+s12*zb1*zb2*zb4+2*M**2
     -   *P(2)*P(8)*P(25)-p1sq*P(55)
       P(57) = s12*za1*za2*za4-s23*za1*za2*za4+p3sq*za1*za4*zb2-s12*za
     -   1*za4*zb2-p3sq*za4*zb1*zb2+s23*za4*zb1*zb2-p3sq*za1*za2*zb4+
     -   s23*za1*za2*zb4+p3sq*za2*zb1*zb4-s12*za2*zb1*zb4+s12*zb1*zb2
     -   *zb4-s23*zb1*zb2*zb4+p2sq*P(2)*P(3)-p4sq*P(8)*P(19)
       P(58) = 2*za2*za4-za4*zb2-za2*zb4
       P(59) = za2-2*zb2
       P(60) = za4*zb2+zb4*P(59)
       P(61) = za1*P(58)-zb1*P(60)
       P(62) = p1sq-s12
       P(63) = 2*za3*za4-za3*zb4+zb3*zb4
       P(64) = za3*za4-za4*zb3+2*zb3*zb4
       P(65) = zb1*P(63)+za1*P(64)
       P(66) = P(11)*P(62)-M**2*P(8)*P(65)
       P(67) = za1*za2+zb1*zb2
       P(68) = p4sq*P(11)-M**2*P(2)*P(38)*P(67)
       P(69) = M**2*P(2)*P(8)*P(37)
       P(70) = p1sq*za1*za2*za3*za4-p2sq*za1*za2*za3*za4-s12*za1*za2*z
     -   a3*za4-p1sq*zb1*zb2*zb3*zb4+p2sq*zb1*zb2*zb3*zb4+s12*zb1*zb2
     -   *zb3*zb4-2*P(69)
       P(71) = -(p1sq*za1*za2*za3*za4)+p4sq*za1*za2*za3*za4+s23*za1*za
     -   2*za3*za4+p1sq*zb1*zb2*zb3*zb4-p4sq*zb1*zb2*zb3*zb4-s23*zb1*
     -   zb2*zb3*zb4+2*P(69)
       P(72) = p1sq*za1*za2*za3*za4-p4sq*za1*za2*za3*za4-s23*za1*za2*z
     -   a3*za4-p1sq*zb1*zb2*zb3*zb4+p4sq*zb1*zb2*zb3*zb4+s23*zb1*zb2
     -   *zb3*zb4-2*P(69)
       P(73) = p4sq*za1*za4*zb2+s23*za4*zb1*zb2+s23*za1*za2*zb4+p4sq*z
     -   a2*zb1*zb4+2*M**2*P(2)*P(8)*P(25)-p1sq*P(55)
       P(74) = -2*za4+zb4
       P(75) = -(zb3*zb4)+za3*P(74)
       P(76) = zb2*P(64)+za2*P(75)
       P(77) = p1sq*P(11)+M**2*P(25)*P(76)
       P(78) = p2sq-s23
       P(79) = -(M**2*P(4)*P(8)*P(19))+P(11)*P(78)
       P(80) = s12*za1*za4*zb2+p2sq*za4*zb1*zb2+p2sq*za1*za2*zb4+s12*z
     -   a2*zb1*zb4+2*M**2*P(2)*P(8)*P(25)-p1sq*P(55)
       P(81) = za4*zb1*zb2+za1*za2*zb4
       P(82) = za4*zb2-za2*zb4
       P(83) = za1*za4*zb2+za2*zb1*zb4
       P(84) = za4*zb2+za2*zb4
       P(85) = -(za2*za4)+za4*zb2+2*za2*zb4
       P(86) = 2*za4*zb2+zb4*P(8)
       P(87) = za1*P(85)+zb1*P(86)
       P(88) = za4-2*zb4
       P(89) = zb1*zb2*P(74)+za1*za2*P(88)
       P(90) = za4*zb2+za2*P(88)
       P(91) = za2+zb2
       P(92) = -2*za4*zb2+zb4*P(91)
       P(93) = za1*P(90)+zb1*P(92)
       P(94) = za1*za4*zb3+za3*zb1*zb4
       P(95) = za3*za4*zb1+za1*zb3*zb4
       P(96) = P(62)*P(94)+p3sq*P(95)-s23*P(95)
       P(97) = za3*za4-zb3*zb4
       P(98) = -(s12*za1*za3*za4)-s23*za1*za3*za4-p2sq*za3*za4*zb1+s23
     -   *za3*za4*zb1+p2sq*za1*za4*zb3+s12*za1*za4*zb3+p2sq*za3*zb1*z
     -   b4+s12*za3*zb1*zb4-p2sq*za1*zb3*zb4+s23*za1*zb3*zb4-s12*zb1*
     -   zb3*zb4-s23*zb1*zb3*zb4+p1sq*P(4)*P(19)+p3sq*P(5)*P(97)
       P(99) = p4sq*P(12)+M**2*P(2)*P(50)
       P(100) = p4sq*P(11)-M**2*P(2)*P(50)
       P(101) = p1sq-p4sq-s23
       P(102) = 2*M**2*P(2)*P(3)*P(4)+P(11)*P(101)
       P(103) = p3sq+p4sq-s12
       P(104) = -2*M**2*P(2)*P(3)*P(4)+P(11)*P(103)
       P(105) = -2*M**2*P(2)*P(3)*P(4)-P(11)*P(101)
       P(106) = 2*M**2*P(2)*P(3)*P(4)-P(11)*P(103)
       P(107) = -(p4sq*za1*za3*za4)-s23*za3*za4*zb1+p1sq*za1*za4*zb3+p
     -   1sq*za3*zb1*zb4-s23*za1*zb3*zb4-p4sq*zb1*zb3*zb4+2*M**2*P(2)
     -   *P(4)*P(25)
       P(108) = zb3*zb4+za3*P(88)
       P(109) = za3*za4+zb3*P(74)
       P(110) = zb1*P(108)+za1*P(109)
       P(111) = za3*za4+zb3*zb4
       P(112) = p1sq*za1*za3*za4-s12*za1*za3*za4-p2sq*za3*za4*zb1-p2sq
     -   *za1*zb3*zb4+p1sq*zb1*zb3*zb4-s12*zb1*zb3*zb4
       P(113) = p1sq*P(12)+M**2*P(25)*P(91)*P(97)
       P(114) = s12*P(11)-M**2*P(25)*P(91)*P(97)
       P(115) = za1*za4+zb1*zb4
       P(116) = za1*za3*za4+zb1*zb3*zb4
       P(117) = p1sq*p3sq+p2sq*p4sq-s12*s23
       P(118) = p3sq*za1*za2*za3*za4+p4sq*za1*za2*za3*za4-s12*za1*za2*
     -   za3*za4-s23*za1*za2*za3*za4+p3sq*za2*za3*za4*zb1-p4sq*za1*za
     -   3*za4*zb2-s23*za3*za4*zb1*zb2+s12*za1*za4*zb2*zb3-s12*za2*za
     -   3*zb1*zb4+s23*za1*za2*zb3*zb4+p4sq*za2*zb1*zb3*zb4-p3sq*za1*
     -   zb2*zb3*zb4-p3sq*zb1*zb2*zb3*zb4-p4sq*zb1*zb2*zb3*zb4+s12*zb
     -   1*zb2*zb3*zb4+s23*zb1*zb2*zb3*zb4+p2sq*P(2)*P(36)+p1sq*P(4)*
     -   P(55)
       P(119) = -4*M**4*P(2)*P(4)*P(8)*P(25)-P(11)*P(117)+2*M**2*P(118
     -   )
       P(120) = s12+s23
       P(121) = p1sq+p2sq+3*p3sq+p4sq-2*P(120)
       P(122) = za3*za4-2*za4*zb3+zb3*zb4
       P(123) = 3*za4*zb3+za3*zb4-2*P(111)
       P(124) = zb2*P(122)+za2*P(123)
       P(125) = za3*za4-2*za3*zb4+zb3*zb4
       P(126) = za4*zb3+3*za3*zb4-2*P(111)
       P(127) = za2*P(125)+zb2*P(126)
       P(128) = -(za1*P(124))+zb1*P(127)
       P(129) = -(P(11)*P(121))+2*M**2*P(128)
       P(130) = p1sq+p2sq-s12
       P(131) = -2*M**2*P(8)*P(25)*P(111)+P(11)*P(130)
       P(132) = p1sq+p4sq-s23
       P(133) = -2*M**2*P(2)*P(25)*P(35)+P(11)*P(132)
       FR(1) = (8*(C123R(4,2)+C124R(4,2))*P(11)+2*(B013R+B024R)*P(12)+
     -   C0234R*P(23)+C0134R*P(104)+D01234R*P(119)+4*D1234R(7,2)*P(12
     -   9)+C0123R*P(131)+C0124R*P(133))/2.
       FI(1) = (8*(C123I(4,2)+C124I(4,2))*P(11)+2*(B013I+B024I)*P(12)+
     -   C0234I*P(23)+C0134I*P(104)+D01234I*P(119)+4*D1234I(7,2)*P(12
     -   9)+C0123I*P(131)+C0124I*P(133))/2.
       F(1)=DCMPLX(FR(1),FI(1))
       P(134) = s12*P(11)+p1sq*P(12)+2*M**2*P(8)*P(95)
       P(135) = -(P(11)*P(20))-2*M**2*P(5)*P(8)*P(97)
       P(136) = p2sq-p3sq+s23
       P(137) = 2*M**2*P(8)*P(38)*P(115)-P(11)*P(136)
       FR(2) = (C0134R-2*(C124R(1,1)-C134R(2,1)+2*(-C124R(2,2)+C124R(3
     -   ,2)+D1234R(7,2)+2*D1234R(13,3))))*P(11)+C0124R*P(12)+2*(D123
     -   4R(5,2)*P(23)+D1234R(3,2)*P(28))+D01234R*P(134)+D1234R(1,1)*
     -   P(135)+D1234R(3,1)*P(137)
       FI(2) = (C0134I-2*(C124I(1,1)-C134I(2,1)+2*(-C124I(2,2)+C124I(3
     -   ,2)+D1234I(7,2)+2*D1234I(13,3))))*P(11)+C0124I*P(12)+2*(D123
     -   4I(5,2)*P(23)+D1234I(3,2)*P(28))+D01234I*P(134)+D1234I(1,1)*
     -   P(135)+D1234I(3,1)*P(137)
       F(2)=DCMPLX(FR(2),FI(2))
       P(138) = p3sq-s23
       P(139) = P(11)*P(138)
       P(140) = 2*M**2*P(8)*P(94)+P(139)
       P(141) = p3sq*P(11)+s23*P(12)+2*M**2*P(8)*P(94)
       FR(3) = 2*((C124R(1,1)-C124R(2,1)-C134R(1,1)+C134R(2,1)+2*(C124
     -   R(1,2)+C124R(2,2)-2*(C124R(3,2)-D1234R(11,3)+D1234R(13,3))))
     -   *P(11)+2*D1234R(5,2)*P(23)+(D1234R(1,2)+D1234R(3,2))*P(28)-D
     -   1234R(1,1)*P(140)+D1234R(3,1)*P(141))
       FI(3) = 2*((C124I(1,1)-C124I(2,1)-C134I(1,1)+C134I(2,1)+2*(C124
     -   I(1,2)+C124I(2,2)-2*(C124I(3,2)-D1234I(11,3)+D1234I(13,3))))
     -   *P(11)+2*D1234I(5,2)*P(23)+(D1234I(1,2)+D1234I(3,2))*P(28)-D
     -   1234I(1,1)*P(140)+D1234I(3,1)*P(141))
       F(3)=DCMPLX(FR(3),FI(3))
       P(142) = p4sq*P(11)-2*M**2*P(2)*P(36)
       P(143) = -2*M**2*P(2)*P(5)*P(6)-P(11)*P(101)
       P(144) = M**2*P(8)*P(9)*P(10)
       P(145) = p1sq*za1*za2*za3*za4-p3sq*za1*za2*za3*za4-s12*za1*za2*
     -   za3*za4+s23*za1*za2*za3*za4-p1sq*zb1*zb2*zb3*zb4+p3sq*zb1*zb
     -   2*zb3*zb4+s12*zb1*zb2*zb3*zb4-s23*zb1*zb2*zb3*zb4-2*P(144)
       P(146) = p3sq-p4sq+s12-2*s23
       P(147) = za4*zb2-2*za2*zb4+zb2*zb4
       P(148) = zb2*zb3*P(15)+za3*P(147)
       P(149) = za3*P(2)+zb3*P(10)
       P(150) = -2*za4*zb2*zb3+za2*P(149)
       P(151) = -(zb1*P(148))+za1*P(150)
       P(152) = P(11)*P(146)+2*M**2*P(151)
       FR(4) = (C0123R+C0234R)*P(12)-2*((C123R(1,1)+C134R(1,1)-C134R(2
     -   ,1)-4*(D1234R(12,3)-D1234R(13,3)))*P(11)-(D1234R(5,2)+D1234R
     -   (6,2))*P(23)-(D1234R(3,2)+D1234R(4,2))*P(28))+D01234R*P(142)
     -   +D1234R(1,1)*P(143)+D1234R(2,1)*P(145)+D1234R(3,1)*P(152)
       FI(4) = (C0123I+C0234I)*P(12)-2*((C123I(1,1)+C134I(1,1)-C134I(2
     -   ,1)-4*(D1234I(12,3)-D1234I(13,3)))*P(11)-(D1234I(5,2)+D1234I
     -   (6,2))*P(23)-(D1234I(3,2)+D1234I(4,2))*P(28))+D01234I*P(142)
     -   +D1234I(1,1)*P(143)+D1234I(2,1)*P(145)+D1234I(3,1)*P(152)
       F(4)=DCMPLX(FR(4),FI(4))
       P(153) = za2*za3*za4-zb2*zb3*zb4
       P(154) = p1sq*P(12)+2*M**2*P(25)*P(153)
       P(155) = p2sq-p4sq+s12+s23
       P(156) = za3*za4*zb2-za4*zb2*zb3+za2*za3*zb4+za2*zb3*zb4-2*zb2*
     -   zb3*zb4
       P(157) = 2*za3*za4+za3*zb4-zb3*zb4
       P(158) = -(za4*zb2*P(38))+za2*P(157)
       P(159) = zb1*P(156)+za1*P(158)
       P(160) = -(P(11)*P(155))+2*M**2*P(159)
       FR(5) = (C0134R+C0234R+2*(C134R(2,1)-4*D1234R(7,2)-8*D1234R(13,
     -   3)))*P(11)+2*(D1234R(6,2)*P(23)+D1234R(3,2)*P(28))+D1234R(2,
     -   1)*P(135)+D01234R*P(154)+D1234R(3,1)*P(160)
       FI(5) = (C0134I+C0234I+2*(C134I(2,1)-4*D1234I(7,2)-8*D1234I(13,
     -   3)))*P(11)+2*(D1234I(6,2)*P(23)+D1234I(3,2)*P(28))+D1234I(2,
     -   1)*P(135)+D01234I*P(154)+D1234I(3,1)*P(160)
       F(5)=DCMPLX(FR(5),FI(5))
       P(161) = p4sq*P(12)+2*M**2*P(2)*P(36)
       P(162) = p1sq-p4sq+s23
       P(163) = -2*M**2*P(6)*P(10)*P(25)+P(11)*P(162)
       P(164) = -(p1sq*za1*za2*za3*za4)-p3sq*za1*za2*za3*za4+s12*za1*z
     -   a2*za3*za4+s23*za1*za2*za3*za4+p1sq*zb1*zb2*zb3*zb4+p3sq*zb1
     -   *zb2*zb3*zb4-s12*zb1*zb2*zb3*zb4-s23*zb1*zb2*zb3*zb4+2*P(69)
       P(165) = p3sq+p4sq-s12-2*s23
       P(166) = za4*zb2+2*za2*zb4-zb2*zb4
       P(167) = -(zb2*zb3*P(10))+za3*P(166)
       P(168) = zb3*P(2)+za3*P(10)
       P(169) = -2*za4*zb2*zb3+za2*P(168)
       P(170) = zb1*P(167)+za1*P(169)
       P(171) = P(11)*P(165)+2*M**2*P(170)
       FR(6) = (C0123R+C0234R-2*(C134R(1,1)-C134R(2,1)-8*(D1234R(11,3)
     -   -D1234R(13,3))))*P(11)+2*((D1234R(5,2)+D1234R(6,2))*P(23)+(D
     -   1234R(3,2)+D1234R(4,2))*P(28))+D01234R*P(161)+D1234R(1,1)*P(
     -   163)+D1234R(2,1)*P(164)+D1234R(3,1)*P(171)
       FI(6) = (C0123I+C0234I-2*(C134I(1,1)-C134I(2,1)-8*(D1234I(11,3)
     -   -D1234I(13,3))))*P(11)+2*((D1234I(5,2)+D1234I(6,2))*P(23)+(D
     -   1234I(3,2)+D1234I(4,2))*P(28))+D01234I*P(161)+D1234I(1,1)*P(
     -   163)+D1234I(2,1)*P(164)+D1234I(3,1)*P(171)
       F(6)=DCMPLX(FR(6),FI(6))
       P(172) = s23*P(11)-2*M**2*P(6)*P(115)
       FR(7) = 2*(-((C123R(2,1)+C134R(1,1)-C134R(2,1)-8*(D1234R(12,3)-
     -   D1234R(13,3)))*P(11))+2*D1234R(6,2)*P(23)+(D1234R(2,2)+D1234
     -   R(3,2))*P(28)+(D1234R(2,1)-D1234R(3,1))*P(172))
       FI(7) = 2*(-((C123I(2,1)+C134I(1,1)-C134I(2,1)-8*(D1234I(12,3)-
     -   D1234I(13,3)))*P(11))+2*D1234I(6,2)*P(23)+(D1234I(2,2)+D1234
     -   I(3,2))*P(28)+(D1234I(2,1)-D1234I(3,1))*P(172))
       F(7)=DCMPLX(FR(7),FI(7))
       P(173) = p2sq*P(11)-2*M**2*P(8)*P(116)
       FR(8) = 2*((C124R(2,1)+2*C124R(2,2)+C134R(2,1)-4*D1234R(13,3))*
     -   P(11)+D1234R(3,2)*P(28)-D1234R(3,1)*P(173))
       FI(8) = 2*((C124I(2,1)+2*C124I(2,2)+C134I(2,1)-4*D1234I(13,3))*
     -   P(11)+D1234I(3,2)*P(28)-D1234I(3,1)*P(173))
       F(8)=DCMPLX(FR(8),FI(8))
       P(174) = p1sq*P(11)+s12*P(12)-2*M**2*P(8)*P(95)
       FR(9) = (C0124R+2*(-C134R(1,1)+C134R(2,1)+2*(C124R(2,2)-C124R(3
     -   ,2)+D1234R(7,2)+2*(D1234R(11,3)-D1234R(13,3)))))*P(11)+C0134
     -   R*P(12)+2*(D1234R(5,2)*P(23)+D1234R(3,2)*P(28))+D1234R(1,1)*
     -   P(131)+D1234R(3,1)*P(137)+D01234R*P(174)
       FI(9) = (C0124I+2*(-C134I(1,1)+C134I(2,1)+2*(C124I(2,2)-C124I(3
     -   ,2)+D1234I(7,2)+2*(D1234I(11,3)-D1234I(13,3)))))*P(11)+C0134
     -   I*P(12)+2*(D1234I(5,2)*P(23)+D1234I(3,2)*P(28))+D1234I(1,1)*
     -   P(131)+D1234I(3,1)*P(137)+D01234I*P(174)
       F(9)=DCMPLX(FR(9),FI(9))
       P(175) = p1sq*P(11)-2*M**2*P(25)*P(153)
       P(176) = p2sq+p4sq-s12+s23
       P(177) = -(za3*za4*zb2)+za4*zb2*zb3+za2*za3*zb4+za2*zb3*zb4-2*z
     -   b2*zb3*zb4
       P(178) = -(za4*zb2*P(38))+za2*P(63)
       P(179) = zb1*P(177)+za1*P(178)
       P(180) = -(P(11)*P(176))+2*M**2*P(179)
       FR(10) = (C0134R+C0234R)*P(12)+2*((-C134R(1,1)+C134R(2,1)+4*(D1
     -   234R(7,2)+D1234R(12,3)-D1234R(13,3)))*P(11)+D1234R(6,2)*P(23
     -   )+D1234R(3,2)*P(28))+D1234R(2,1)*P(131)+D01234R*P(175)+D1234
     -   R(3,1)*P(180)
       FI(10) = (C0134I+C0234I)*P(12)+2*((-C134I(1,1)+C134I(2,1)+4*(D1
     -   234I(7,2)+D1234I(12,3)-D1234I(13,3)))*P(11)+D1234I(6,2)*P(23
     -   )+D1234I(3,2)*P(28))+D1234I(2,1)*P(131)+D01234I*P(175)+D1234
     -   I(3,1)*P(180)
       F(10)=DCMPLX(FR(10),FI(10))
       P(181) = p1sq*za1*za2*za3*za4+p3sq*za1*za2*za3*za4-s12*za1*za2*
     -   za3*za4-s23*za1*za2*za3*za4-p1sq*zb1*zb2*zb3*zb4-p3sq*zb1*zb
     -   2*zb3*zb4+s12*zb1*zb2*zb3*zb4+s23*zb1*zb2*zb3*zb4-2*P(69)
       FR(11) = 2*(D1234R(1,1)+D1234R(1,2)-D1234R(3,1)+D1234R(3,2)-2*D
     -   1234R(5,2))*P(181)
       FI(11) = 2*(D1234I(1,1)+D1234I(1,2)-D1234I(3,1)+D1234I(3,2)-2*D
     -   1234I(5,2))*P(181)
       F(11)=DCMPLX(FR(11),FI(11))
       FR(12) = (C0134R-2*(C124R(1,1)-C124R(2,1))-4*D1234R(7,2))*P(11)
     -   +C0124R*P(12)+D1234R(3,1)*P(28)+D01234R*P(134)+D1234R(1,1)*P
     -   (135)+2*(D1234R(3,2)-D1234R(5,2))*P(181)
       FI(12) = (C0134I-2*(C124I(1,1)-C124I(2,1))-4*D1234I(7,2))*P(11)
     -   +C0124I*P(12)+D1234I(3,1)*P(28)+D01234I*P(134)+D1234I(1,1)*P
     -   (135)+2*(D1234I(3,2)-D1234I(5,2))*P(181)
       F(12)=DCMPLX(FR(12),FI(12))
       P(182) = 2*M**2*P(2)*P(5)*P(6)+P(11)*P(101)
       P(183) = p1sq-s23
       P(184) = p3sq-p4sq-s12+2*P(183)
       P(185) = -2*M**2*P(2)*P(41)-P(11)*P(184)
       FR(13) = (C0123R+C0234R+2*C123R(1,1)-4*D1234R(7,2))*P(11)+D0123
     -   4R*P(161)+(D1234R(2,1)+2*(D1234R(3,2)+D1234R(4,2)-D1234R(5,2
     -   )-D1234R(6,2)))*P(181)+D1234R(1,1)*P(182)+D1234R(3,1)*P(185)
       FI(13) = (C0123I+C0234I+2*C123I(1,1)-4*D1234I(7,2))*P(11)+D0123
     -   4I*P(161)+(D1234I(2,1)+2*(D1234I(3,2)+D1234I(4,2)-D1234I(5,2
     -   )-D1234I(6,2)))*P(181)+D1234I(1,1)*P(182)+D1234I(3,1)*P(185)
       F(13)=DCMPLX(FR(13),FI(13))
       P(186) = p1sq*P(11)+s23*P(12)-2*M**2*P(2)*P(46)
       P(187) = s23*P(11)+p1sq*P(12)+2*M**2*P(2)*P(46)
       FR(14) = 2*((C123R(2,1)-C234R(1,1)+C234R(2,1)+4*(D1234R(12,3)-D
     -   1234R(13,3)))*P(11)+(D1234R(2,2)+D1234R(3,2)-2*D1234R(6,2))*
     -   P(181)+D1234R(2,1)*P(186)+D1234R(3,1)*P(187))
       FI(14) = 2*((C123I(2,1)-C234I(1,1)+C234I(2,1)+4*(D1234I(12,3)-D
     -   1234I(13,3)))*P(11)+(D1234I(2,2)+D1234I(3,2)-2*D1234I(6,2))*
     -   P(181)+D1234I(2,1)*P(186)+D1234I(3,1)*P(187))
       F(14)=DCMPLX(FR(14),FI(14))
       P(188) = p2sq+p4sq+s12-s23
       P(189) = za3*za4*zb2+za4*zb2*zb3-za2*za3*zb4+za2*zb3*zb4-2*zb2*
     -   zb3*zb4
       P(190) = -2*za3*za4+za3*zb4+zb3*zb4
       P(191) = za4*zb2*P(4)+za2*P(190)
       P(192) = zb1*P(189)-za1*P(191)
       P(193) = -(P(11)*P(188))+2*M**2*P(192)
       FR(15) = (C0134R+C0234R+2*(C234R(2,1)-4*(D1234R(7,2)+D1234R(13,
     -   3))))*P(11)+D1234R(2,1)*P(135)+D01234R*P(154)+2*(D1234R(3,2)
     -   -D1234R(6,2))*P(181)+D1234R(3,1)*P(193)
       FI(15) = (C0134I+C0234I+2*(C234I(2,1)-4*(D1234I(7,2)+D1234I(13,
     -   3))))*P(11)+D1234I(2,1)*P(135)+D01234I*P(154)+2*(D1234I(3,2)
     -   -D1234I(6,2))*P(181)+D1234I(3,1)*P(193)
       F(15)=DCMPLX(FR(15),FI(15))
       P(194) = p3sq+p4sq-s12+2*P(183)
       P(195) = 2*M**2*P(2)*P(45)-P(11)*P(194)
       FR(16) = (C0234R+2*C234R(2,1)+4*D1234R(7,2)+8*(D1234R(11,3)-D12
     -   34R(13,3)))*P(11)+C0123R*P(12)+D1234R(1,1)*P(133)+D01234R*P(
     -   142)+(D1234R(2,1)+2*(D1234R(3,2)+D1234R(4,2)-D1234R(5,2)-D12
     -   34R(6,2)))*P(181)+D1234R(3,1)*P(195)
       FI(16) = (C0234I+2*C234I(2,1)+4*D1234I(7,2)+8*(D1234I(11,3)-D12
     -   34I(13,3)))*P(11)+C0123I*P(12)+D1234I(1,1)*P(133)+D01234I*P(
     -   142)+(D1234I(2,1)+2*(D1234I(3,2)+D1234I(4,2)-D1234I(5,2)-D12
     -   34I(6,2)))*P(181)+D1234I(3,1)*P(195)
       F(16)=DCMPLX(FR(16),FI(16))
       FR(17) = 2*((C124R(2,1)-C134R(2,1)+4*D1234R(13,3))*P(11)+D1234R
     -   (3,1)*P(174)+D1234R(3,2)*P(181))
       FI(17) = 2*((C124I(2,1)-C134I(2,1)+4*D1234I(13,3))*P(11)+D1234I
     -   (3,1)*P(174)+D1234I(3,2)*P(181))
       F(17)=DCMPLX(FR(17),FI(17))
       P(196) = 2*M**2*P(8)*P(25)*P(111)-P(11)*P(130)
       P(197) = -2*M**2*P(8)*P(38)*P(115)+P(11)*P(136)
       FR(18) = (C0134R+2*(C134R(1,1)-C134R(2,1))-4*D1234R(7,2)-8*(D12
     -   34R(11,3)-D1234R(13,3)))*P(11)+C0124R*P(12)+D01234R*P(134)+2
     -   *(D1234R(3,2)-D1234R(5,2))*P(181)+D1234R(1,1)*P(196)+D1234R(
     -   3,1)*P(197)
       FI(18) = (C0134I+2*(C134I(1,1)-C134I(2,1))-4*D1234I(7,2)-8*(D12
     -   34I(11,3)-D1234I(13,3)))*P(11)+C0124I*P(12)+D01234I*P(134)+2
     -   *(D1234I(3,2)-D1234I(5,2))*P(181)+D1234I(1,1)*P(196)+D1234I(
     -   3,1)*P(197)
       F(18)=DCMPLX(FR(18),FI(18))
       P(198) = P(11)*P(176)-2*M**2*P(179)
       FR(19) = (C0134R+C0234R+2*(C134R(1,1)-C134R(2,1)-4*(D1234R(7,2)
     -   +D1234R(12,3)-D1234R(13,3))))*P(11)+D01234R*P(154)+2*(D1234R
     -   (3,2)-D1234R(6,2))*P(181)+D1234R(2,1)*P(196)+D1234R(3,1)*P(1
     -   98)
       FI(19) = (C0134I+C0234I+2*(C134I(1,1)-C134I(2,1)-4*(D1234I(7,2)
     -   +D1234I(12,3)-D1234I(13,3))))*P(11)+D01234I*P(154)+2*(D1234I
     -   (3,2)-D1234I(6,2))*P(181)+D1234I(2,1)*P(196)+D1234I(3,1)*P(1
     -   98)
       F(19)=DCMPLX(FR(19),FI(19))
       FR(20) = -2*(C123R(1,1)-2*D1234R(7,2)-4*(D1234R(11,3)-D1234R(13
     -   ,3)))*P(11)+(C0123R+C0234R)*P(12)+D1234R(3,1)*P(106)+D01234R
     -   *P(142)+D1234R(1,1)*P(143)+D1234R(2,1)*P(181)
       FI(20) = -2*(C123I(1,1)-2*D1234I(7,2)-4*(D1234I(11,3)-D1234I(13
     -   ,3)))*P(11)+(C0123I+C0234I)*P(12)+D1234I(3,1)*P(106)+D01234I
     -   *P(142)+D1234I(1,1)*P(143)+D1234I(2,1)*P(181)
       F(20)=DCMPLX(FR(20),FI(20))
       FR(21) = 8*(D1234R(1,2)+D1234R(1,3)+D1234R(3,2)-D1234R(3,3)-2*D
     -   1234R(5,2)-3*(D1234R(5,3)-D1234R(7,3)))*P(11)
       FI(21) = 8*(D1234I(1,2)+D1234I(1,3)+D1234I(3,2)-D1234I(3,3)-2*D
     -   1234I(5,2)-3*(D1234I(5,3)-D1234I(7,3)))*P(11)
       F(21)=DCMPLX(FR(21),FI(21))
       FR(22) = -4*(D1234R(1,1)+D1234R(1,2)-D1234R(3,1)-D1234R(5,2)+2*
     -   (D1234R(3,3)+D1234R(5,3)-2*D1234R(7,3)))*P(11)
       FI(22) = -4*(D1234I(1,1)+D1234I(1,2)-D1234I(3,1)-D1234I(5,2)+2*
     -   (D1234I(3,3)+D1234I(5,3)-2*D1234I(7,3)))*P(11)
       F(22)=DCMPLX(FR(22),FI(22))
       FR(23) = -4*(D1234R(1,1)+D1234R(1,2)-D1234R(3,1)-D1234R(3,2)+2*
     -   (D1234R(3,3)+D1234R(5,3)-2*D1234R(7,3)))*P(11)
       FI(23) = -4*(D1234I(1,1)+D1234I(1,2)-D1234I(3,1)-D1234I(3,2)+2*
     -   (D1234I(3,3)+D1234I(5,3)-2*D1234I(7,3)))*P(11)
       F(23)=DCMPLX(FR(23),FI(23))
       FR(24) = 4*(D1234R(3,1)-D1234R(3,2)+2*(-D1234R(3,3)+D1234R(5,2)
     -   +D1234R(7,3)))*P(11)
       FI(24) = 4*(D1234I(3,1)-D1234I(3,2)+2*(-D1234I(3,3)+D1234I(5,2)
     -   +D1234I(7,3)))*P(11)
       F(24)=DCMPLX(FR(24),FI(24))
       FR(25) = 8*(D1234R(2,2)+D1234R(3,2)-D1234R(3,3)+D1234R(6,3)+D12
     -   34R(7,3)-D1234R(8,3)-2*(D1234R(6,2)-D1234R(9,3)+D1234R(10,3)
     -   ))*P(11)
       FI(25) = 8*(D1234I(2,2)+D1234I(3,2)-D1234I(3,3)+D1234I(6,3)+D12
     -   34I(7,3)-D1234I(8,3)-2*(D1234I(6,2)-D1234I(9,3)+D1234I(10,3)
     -   ))*P(11)
       F(25)=DCMPLX(FR(25),FI(25))
       FR(26) = 8*(D1234R(3,2)-D1234R(3,3)+D1234R(4,2)+D1234R(4,3)-D12
     -   34R(5,2)-D1234R(5,3)-D1234R(6,2)+D1234R(9,3)-2*(-D1234R(7,3)
     -   +D1234R(10,3)))*P(11)
       FI(26) = 8*(D1234I(3,2)-D1234I(3,3)+D1234I(4,2)+D1234I(4,3)-D12
     -   34I(5,2)-D1234I(5,3)-D1234I(6,2)+D1234I(9,3)-2*(-D1234I(7,3)
     -   +D1234I(10,3)))*P(11)
       F(26)=DCMPLX(FR(26),FI(26))
       FR(27) = 4*(-D1234R(2,1)+D1234R(3,1)-D1234R(4,2)+D1234R(5,2)-2*
     -   (D1234R(3,3)-D1234R(7,3)-D1234R(9,3)+D1234R(10,3)))*P(11)
       FI(27) = 4*(-D1234I(2,1)+D1234I(3,1)-D1234I(4,2)+D1234I(5,2)-2*
     -   (D1234I(3,3)-D1234I(7,3)-D1234I(9,3)+D1234I(10,3)))*P(11)
       F(27)=DCMPLX(FR(27),FI(27))
       FR(28) = 4*(-D1234R(2,1)+D1234R(3,1)+D1234R(3,2)-D1234R(4,2)+D1
     -   234R(5,2)-D1234R(6,2)-2*(D1234R(3,3)-D1234R(7,3)-D1234R(9,3)
     -   +D1234R(10,3)))*P(11)
       FI(28) = 4*(-D1234I(2,1)+D1234I(3,1)+D1234I(3,2)-D1234I(4,2)+D1
     -   234I(5,2)-D1234I(6,2)-2*(D1234I(3,3)-D1234I(7,3)-D1234I(9,3)
     -   +D1234I(10,3)))*P(11)
       F(28)=DCMPLX(FR(28),FI(28))
       FR(29) = (C0134R-2*(C124R(1,1)-C124R(2,1))-4*D1234R(7,2))*P(11)
     -   +C0124R*P(12)+D1234R(3,1)*P(28)+D01234R*P(134)+D1234R(1,1)*P
     -   (135)
       FI(29) = (C0134I-2*(C124I(1,1)-C124I(2,1))-4*D1234I(7,2))*P(11)
     -   +C0124I*P(12)+D1234I(3,1)*P(28)+D01234I*P(134)+D1234I(1,1)*P
     -   (135)
       F(29)=DCMPLX(FR(29),FI(29))
       FR(30) = 4*(D1234R(3,1)+D1234R(5,2))*P(11)
       FI(30) = 4*(D1234I(3,1)+D1234I(5,2))*P(11)
       F(30)=DCMPLX(FR(30),FI(30))
       FR(31) = 4*(-D1234R(1,1)-D1234R(1,2)+D1234R(3,1)+D1234R(5,2))*P
     -   (11)
       FI(31) = 4*(-D1234I(1,1)-D1234I(1,2)+D1234I(3,1)+D1234I(5,2))*P
     -   (11)
       F(31)=DCMPLX(FR(31),FI(31))
       FR(32) = 4*(D1234R(3,2)-D1234R(5,2))*P(11)
       FI(32) = 4*(D1234I(3,2)-D1234I(5,2))*P(11)
       F(32)=DCMPLX(FR(32),FI(32))
       FR(33) = 4*(D1234R(1,1)+D1234R(1,2)-D1234R(3,1)+D1234R(3,2)-2*D
     -   1234R(5,2))*P(11)
       FI(33) = 4*(D1234I(1,1)+D1234I(1,2)-D1234I(3,1)+D1234I(3,2)-2*D
     -   1234I(5,2))*P(11)
       F(33)=DCMPLX(FR(33),FI(33))
       FR(34) = 4*(D1234R(3,2)-D1234R(6,2))*P(11)
       FI(34) = 4*(D1234I(3,2)-D1234I(6,2))*P(11)
       F(34)=DCMPLX(FR(34),FI(34))
       P(199) = M**2*P(4)*P(5)*P(82)
       P(200) = -(p2sq*za1*za2*za3*za4)+p4sq*za1*za2*za3*za4-s12*za1*z
     -   a2*za3*za4+s23*za1*za2*za3*za4+p2sq*zb1*zb2*zb3*zb4-p4sq*zb1
     -   *zb2*zb3*zb4+s12*zb1*zb2*zb3*zb4-s23*zb1*zb2*zb3*zb4-2*P(199
     -   )
       FR(35) = (C0134R+C0234R-4*D1234R(7,2)+8*D1234R(13,3))*P(11)+D12
     -   34R(2,1)*P(135)+D01234R*P(154)+D1234R(3,1)*P(200)
       FI(35) = (C0134I+C0234I-4*D1234I(7,2)+8*D1234I(13,3))*P(11)+D12
     -   34I(2,1)*P(135)+D01234I*P(154)+D1234I(3,1)*P(200)
       F(35)=DCMPLX(FR(35),FI(35))
       FR(36) = 8*(D1234R(3,2)+D1234R(3,3))*P(11)
       FI(36) = 8*(D1234I(3,2)+D1234I(3,3))*P(11)
       F(36)=DCMPLX(FR(36),FI(36))
       FR(37) = -4*(D1234R(3,1)+D1234R(5,2)-2*(D1234R(3,3)-D1234R(7,3)
     -   ))*P(11)
       FI(37) = -4*(D1234I(3,1)+D1234I(5,2)-2*(D1234I(3,3)-D1234I(7,3)
     -   ))*P(11)
       F(37)=DCMPLX(FR(37),FI(37))
       FR(38) = 4*(D1234R(3,2)-D1234R(5,2)+2*(D1234R(3,3)-D1234R(7,3))
     -   )*P(11)
       FI(38) = 4*(D1234I(3,2)-D1234I(5,2)+2*(D1234I(3,3)-D1234I(7,3))
     -   )*P(11)
       F(38)=DCMPLX(FR(38),FI(38))
       FR(39) = 4*(-D1234R(3,2)+D1234R(5,2)+2*(D1234R(3,3)+D1234R(5,3)
     -   -2*D1234R(7,3)))*P(11)
       FI(39) = 4*(-D1234I(3,2)+D1234I(5,2)+2*(D1234I(3,3)+D1234I(5,3)
     -   -2*D1234I(7,3)))*P(11)
       F(39)=DCMPLX(FR(39),FI(39))
       FR(40) = 8*(D1234R(3,3)-D1234R(9,3))*P(11)
       FI(40) = 8*(D1234I(3,3)-D1234I(9,3))*P(11)
       F(40)=DCMPLX(FR(40),FI(40))
       FR(41) = 8*(D1234R(3,2)+D1234R(3,3)-D1234R(6,2)-D1234R(9,3))*P(
     -   11)
       FI(41) = 8*(D1234I(3,2)+D1234I(3,3)-D1234I(6,2)-D1234I(9,3))*P(
     -   11)
       F(41)=DCMPLX(FR(41),FI(41))
       FR(42) = 4*(D1234R(2,1)-D1234R(3,1)+D1234R(4,2)-D1234R(5,2)+2*(
     -   D1234R(3,3)-D1234R(7,3)-D1234R(9,3)+D1234R(10,3)))*P(11)
       FI(42) = 4*(D1234I(2,1)-D1234I(3,1)+D1234I(4,2)-D1234I(5,2)+2*(
     -   D1234I(3,3)-D1234I(7,3)-D1234I(9,3)+D1234I(10,3)))*P(11)
       F(42)=DCMPLX(FR(42),FI(42))
       FR(43) = 4*(-D1234R(2,1)+D1234R(3,1)-D1234R(3,2)-D1234R(4,2)+D1
     -   234R(5,2)+D1234R(6,2)+2*(D1234R(3,3)-D1234R(7,3)-D1234R(9,3)
     -   +D1234R(10,3)))*P(11)
       FI(43) = 4*(-D1234I(2,1)+D1234I(3,1)-D1234I(3,2)-D1234I(4,2)+D1
     -   234I(5,2)+D1234I(6,2)+2*(D1234I(3,3)-D1234I(7,3)-D1234I(9,3)
     -   +D1234I(10,3)))*P(11)
       F(43)=DCMPLX(FR(43),FI(43))
       FR(44) = 8*(D1234R(3,3)+D1234R(8,3)-2*D1234R(9,3))*P(11)
       FI(44) = 8*(D1234I(3,3)+D1234I(8,3)-2*D1234I(9,3))*P(11)
       F(44)=DCMPLX(FR(44),FI(44))
       P(201) = P(11)*P(20)+2*M**2*P(5)*P(8)*P(97)
       P(202) = p2sq+p3sq+2*p4sq-s23
       P(203) = -(za2*za3)+za3*zb2+za2*zb3-3*zb2*zb3
       P(204) = 2*za4*zb2*zb3+zb4*P(203)
       P(205) = 3*za3*za4-za4*zb3-2*za3*zb4
       P(206) = za4*zb2*P(43)+za2*P(205)
       P(207) = zb1*P(204)+za1*P(206)
       P(208) = P(11)*P(202)-2*M**2*P(207)
       FR(45) = (C0124R+2*(C124R(1,1)-C124R(2,1)-C234R(2,1))+8*D1234R(
     -   13,3))*P(11)+C0134R*P(12)+2*(D1234R(5,2)*P(104)+D1234R(3,2)*
     -   P(106))+D01234R*P(174)+D1234R(1,1)*P(201)+D1234R(3,1)*P(208)
       FI(45) = (C0124I+2*(C124I(1,1)-C124I(2,1)-C234I(2,1))+8*D1234I(
     -   13,3))*P(11)+C0134I*P(12)+2*(D1234I(5,2)*P(104)+D1234I(3,2)*
     -   P(106))+D01234I*P(174)+D1234I(1,1)*P(201)+D1234I(3,1)*P(208)
       F(45)=DCMPLX(FR(45),FI(45))
       P(209) = p3sq+2*p4sq-s12
       P(210) = -2*M**2*P(2)*P(47)-P(11)*P(209)
       P(211) = 2*M**2*P(2)*P(47)+P(11)*P(209)
       FR(46) = 2*((C0123R+3*C123R(1,1)+2*C123R(1,2)-C234R(2,1)-4*(D12
     -   34R(7,2)+D1234R(11,3)-D1234R(13,3)))*P(11)+2*D1234R(5,2)*P(1
     -   04)+(D1234R(1,2)+D1234R(3,2))*P(106)-D01234R*P(142)+D1234R(1
     -   ,1)*P(210)+D1234R(3,1)*P(211))
       FI(46) = 2*((C0123I+3*C123I(1,1)+2*C123I(1,2)-C234I(2,1)-4*(D12
     -   34I(7,2)+D1234I(11,3)-D1234I(13,3)))*P(11)+2*D1234I(5,2)*P(1
     -   04)+(D1234I(1,2)+D1234I(3,2))*P(106)-D01234I*P(142)+D1234I(1
     -   ,1)*P(210)+D1234I(3,1)*P(211))
       F(46)=DCMPLX(FR(46),FI(46))
       P(212) = p1sq+p3sq+2*p4sq-s12-s23
       P(213) = 2*M**2*P(2)*P(50)-P(11)*P(212)
       P(214) = p3sq+3*p4sq-s12
       P(215) = 2*M**2*P(2)*P(53)+P(11)*P(214)
       FR(47) = (C0123R+C0234R+2*(C123R(1,1)+C123R(2,1)+C234R(1,1)-C23
     -   4R(2,1)+2*(C123R(3,2)-D1234R(7,2))-4*(D1234R(12,3)-D1234R(13
     -   ,3))))*P(11)+2*((D1234R(5,2)+D1234R(6,2))*P(104)+(D1234R(3,2
     -   )+D1234R(4,2))*P(106))+D01234R*P(161)+D1234R(1,1)*P(182)+D12
     -   34R(2,1)*P(213)+D1234R(3,1)*P(215)
       FI(47) = (C0123I+C0234I+2*(C123I(1,1)+C123I(2,1)+C234I(1,1)-C23
     -   4I(2,1)+2*(C123I(3,2)-D1234I(7,2))-4*(D1234I(12,3)-D1234I(13
     -   ,3))))*P(11)+2*((D1234I(5,2)+D1234I(6,2))*P(104)+(D1234I(3,2
     -   )+D1234I(4,2))*P(106))+D01234I*P(161)+D1234I(1,1)*P(182)+D12
     -   34I(2,1)*P(213)+D1234I(3,1)*P(215)
       F(47)=DCMPLX(FR(47),FI(47))
       P(216) = P(11)*P(188)-2*M**2*P(192)
       FR(48) = (C0134R+C0234R)*P(12)-2*((C234R(2,1)-4*(D1234R(7,2)+D1
     -   234R(13,3)))*P(11)-D1234R(6,2)*P(104)-D1234R(3,2)*P(106))+D0
     -   1234R*P(175)+D1234R(2,1)*P(201)+D1234R(3,1)*P(216)
       FI(48) = (C0134I+C0234I)*P(12)-2*((C234I(2,1)-4*(D1234I(7,2)+D1
     -   234I(13,3)))*P(11)-D1234I(6,2)*P(104)-D1234I(3,2)*P(106))+D0
     -   1234I*P(175)+D1234I(2,1)*P(201)+D1234I(3,1)*P(216)
       F(48)=DCMPLX(FR(48),FI(48))
       P(217) = 2*M**2*P(2)*P(25)*P(35)-P(11)*P(132)
       P(218) = p1sq-p3sq-2*p4sq+s12-s23
       P(219) = 2*M**2*P(2)*P(54)+P(11)*P(218)
       FR(49) = (C0123R-2*C234R(2,1)+4*(C123R(2,1)+C123R(3,2)-D1234R(7
     -   ,2))-8*(D1234R(11,3)-D1234R(13,3)))*P(11)+C0234R*P(12)+2*((D
     -   1234R(5,2)+D1234R(6,2))*P(104)+(D1234R(3,2)+D1234R(4,2))*P(1
     -   06))+D01234R*P(161)+D1234R(3,1)*P(215)+D1234R(1,1)*P(217)+D1
     -   234R(2,1)*P(219)
       FI(49) = (C0123I-2*C234I(2,1)+4*(C123I(2,1)+C123I(3,2)-D1234I(7
     -   ,2))-8*(D1234I(11,3)-D1234I(13,3)))*P(11)+C0234I*P(12)+2*((D
     -   1234I(5,2)+D1234I(6,2))*P(104)+(D1234I(3,2)+D1234I(4,2))*P(1
     -   06))+D01234I*P(161)+D1234I(3,1)*P(215)+D1234I(1,1)*P(217)+D1
     -   234I(2,1)*P(219)
       F(49)=DCMPLX(FR(49),FI(49))
       FR(50) = 2*((C123R(2,1)+2*C123R(2,2)+C234R(1,1)-C234R(2,1)-4*(D
     -   1234R(12,3)-D1234R(13,3)))*P(11)+2*D1234R(6,2)*P(104)+(D1234
     -   R(2,2)+D1234R(3,2))*P(106)+(-D1234R(2,1)+D1234R(3,1))*P(142)
     -   )
       FI(50) = 2*((C123I(2,1)+2*C123I(2,2)+C234I(1,1)-C234I(2,1)-4*(D
     -   1234I(12,3)-D1234I(13,3)))*P(11)+2*D1234I(6,2)*P(104)+(D1234
     -   I(2,2)+D1234I(3,2))*P(106)+(-D1234I(2,1)+D1234I(3,1))*P(142)
     -   )
       F(50)=DCMPLX(FR(50),FI(50))
       P(220) = s12*P(11)-2*M**2*P(67)*P(97)
       FR(51) = 2*(-((C124R(2,1)+C234R(2,1)-8*D1234R(13,3))*P(11))+D12
     -   34R(3,2)*P(106)+D1234R(3,1)*P(220))
       FI(51) = 2*(-((C124I(2,1)+C234I(2,1)-8*D1234I(13,3))*P(11))+D12
     -   34I(3,2)*P(106)+D1234I(3,1)*P(220))
       F(51)=DCMPLX(FR(51),FI(51))
       P(221) = p1sq+s12
       P(222) = za2*za3*za4+za3*za4*zb2-2*zb2*zb3*zb4
       P(223) = 2*za2*za3*za4-za2*zb3*zb4-zb2*zb3*zb4
       P(224) = zb1*P(222)+za1*P(223)
       P(225) = -(P(11)*P(221))+2*M**2*P(224)
       P(226) = p1sq-p2sq+s12
       P(227) = 2*M**2*P(25)*P(91)*P(97)-P(11)*P(226)
       P(228) = p2sq-p3sq-2*p4sq-s23
       P(229) = za2*za3+za3*zb2-za2*zb3-3*zb2*zb3
       P(230) = 2*za4*zb2*zb3+zb4*P(229)
       P(231) = za4*zb2*P(4)+za2*P(205)
       P(232) = zb1*P(230)+za1*P(231)
       P(233) = -(P(11)*P(228))-2*M**2*P(232)
       FR(52) = (C0124R+C0134R-2*(C234R(2,1)+8*(D1234R(7,2)+D1234R(11,
     -   3)-D1234R(13,3))))*P(11)+2*(D1234R(5,2)*P(104)+D1234R(3,2)*P
     -   (106))+D01234R*P(225)+D1234R(1,1)*P(227)+D1234R(3,1)*P(233)
       FI(52) = (C0124I+C0134I-2*(C234I(2,1)+8*(D1234I(7,2)+D1234I(11,
     -   3)-D1234I(13,3))))*P(11)+2*(D1234I(5,2)*P(104)+D1234I(3,2)*P
     -   (106))+D01234I*P(225)+D1234I(1,1)*P(227)+D1234I(3,1)*P(233)
       F(52)=DCMPLX(FR(52),FI(52))
       P(234) = p2sq-p4sq-s12-s23
       P(235) = za3*za4*zb2+za4*zb2*zb3+za2*za3*zb4-za2*zb3*zb4-2*zb2*
     -   zb3*zb4
       P(236) = 2*za3*za4-za3*zb4-zb3*zb4
       P(237) = za4*zb2*P(4)+za2*P(236)
       P(238) = zb1*P(235)+za1*P(237)
       P(239) = -(P(11)*P(234))-2*M**2*P(238)
       FR(53) = (C0134R+C0234R+2*(C234R(1,1)-C234R(2,1)-4*D1234R(7,2)-
     -   8*(D1234R(12,3)-D1234R(13,3))))*P(11)+2*(D1234R(6,2)*P(104)+
     -   D1234R(3,2)*P(106))+D01234R*P(154)+D1234R(2,1)*P(227)+D1234R
     -   (3,1)*P(239)
       FI(53) = (C0134I+C0234I+2*(C234I(1,1)-C234I(2,1)-4*D1234I(7,2)-
     -   8*(D1234I(12,3)-D1234I(13,3))))*P(11)+2*(D1234I(6,2)*P(104)+
     -   D1234I(3,2)*P(106))+D01234I*P(154)+D1234I(2,1)*P(227)+D1234I
     -   (3,1)*P(239)
       F(53)=DCMPLX(FR(53),FI(53))
       FR(54) = (C0123R+C0234R+2*C123R(1,1)-4*D1234R(7,2)-8*(D1234R(11
     -   ,3)-D1234R(12,3)))*P(11)+D1234R(3,1)*P(104)+D01234R*P(161)+D
     -   1234R(2,1)*P(164)+D1234R(1,1)*P(182)
       FI(54) = (C0123I+C0234I+2*C123I(1,1)-4*D1234I(7,2)-8*(D1234I(11
     -   ,3)-D1234I(12,3)))*P(11)+D1234I(3,1)*P(104)+D01234I*P(161)+D
     -   1234I(2,1)*P(164)+D1234I(1,1)*P(182)
       F(54)=DCMPLX(FR(54),FI(54))
       FR(55) = -4*(D1234R(3,1)-D1234R(6,2)+2*(D1234R(5,2)+D1234R(7,3)
     -   -D1234R(9,3)))*P(11)
       FI(55) = -4*(D1234I(3,1)-D1234I(6,2)+2*(D1234I(5,2)+D1234I(7,3)
     -   -D1234I(9,3)))*P(11)
       F(55)=DCMPLX(FR(55),FI(55))
       FR(56) = 4*(D1234R(1,1)+D1234R(1,2)-D1234R(3,1)-D1234R(5,2)-2*(
     -   -D1234R(5,3)+D1234R(7,3)-D1234R(9,3)+D1234R(10,3)))*P(11)
       FI(56) = 4*(D1234I(1,1)+D1234I(1,2)-D1234I(3,1)-D1234I(5,2)-2*(
     -   -D1234I(5,3)+D1234I(7,3)-D1234I(9,3)+D1234I(10,3)))*P(11)
       F(56)=DCMPLX(FR(56),FI(56))
       FR(57) = 4*(D1234R(1,1)+D1234R(1,2)-D1234R(2,1)-D1234R(4,2)+D12
     -   34R(5,2)-D1234R(6,2)-2*(-D1234R(5,3)+D1234R(7,3)-D1234R(9,3)
     -   +D1234R(10,3)))*P(11)
       FI(57) = 4*(D1234I(1,1)+D1234I(1,2)-D1234I(2,1)-D1234I(4,2)+D12
     -   34I(5,2)-D1234I(6,2)-2*(-D1234I(5,3)+D1234I(7,3)-D1234I(9,3)
     -   +D1234I(10,3)))*P(11)
       F(57)=DCMPLX(FR(57),FI(57))
       FR(58) = -8*(D1234R(1,2)+D1234R(1,3)-D1234R(4,2)-D1234R(4,3)-D1
     -   234R(5,2)+D1234R(6,2)+D1234R(7,3)-D1234R(9,3)+2*(-D1234R(5,3
     -   )+D1234R(10,3)))*P(11)
       FI(58) = -8*(D1234I(1,2)+D1234I(1,3)-D1234I(4,2)-D1234I(4,3)-D1
     -   234I(5,2)+D1234I(6,2)+D1234I(7,3)-D1234I(9,3)+2*(-D1234I(5,3
     -   )+D1234I(10,3)))*P(11)
       F(58)=DCMPLX(FR(58),FI(58))
       FR(59) = 8*(-D1234R(4,3)+D1234R(5,3)+D1234R(6,3)-D1234R(7,3)-D1
     -   234R(8,3)+D1234R(9,3))*P(11)
       FI(59) = 8*(-D1234I(4,3)+D1234I(5,3)+D1234I(6,3)-D1234I(7,3)-D1
     -   234I(8,3)+D1234I(9,3))*P(11)
       F(59)=DCMPLX(FR(59),FI(59))
       FR(60) = 8*(D1234R(2,2)-D1234R(4,2)-D1234R(4,3)+D1234R(5,2)+D12
     -   34R(5,3)-D1234R(6,2)+D1234R(6,3)-D1234R(7,3)-D1234R(8,3)+D12
     -   34R(9,3))*P(11)
       FI(60) = 8*(D1234I(2,2)-D1234I(4,2)-D1234I(4,3)+D1234I(5,2)+D12
     -   34I(5,3)-D1234I(6,2)+D1234I(6,3)-D1234I(7,3)-D1234I(8,3)+D12
     -   34I(9,3))*P(11)
       F(60)=DCMPLX(FR(60),FI(60))
       FR(61) = 4*(D1234R(2,1)-D1234R(3,1)+D1234R(4,2)-D1234R(5,2)+2*(
     -   -D1234R(7,3)-D1234R(8,3)+D1234R(9,3)+D1234R(10,3)))*P(11)
       FI(61) = 4*(D1234I(2,1)-D1234I(3,1)+D1234I(4,2)-D1234I(5,2)+2*(
     -   -D1234I(7,3)-D1234I(8,3)+D1234I(9,3)+D1234I(10,3)))*P(11)
       F(61)=DCMPLX(FR(61),FI(61))
       FR(62) = 4*(-D1234R(2,2)+D1234R(4,2)-D1234R(5,2)+D1234R(6,2)+2*
     -   (-D1234R(7,3)-D1234R(8,3)+D1234R(9,3)+D1234R(10,3)))*P(11)
       FI(62) = 4*(-D1234I(2,2)+D1234I(4,2)-D1234I(5,2)+D1234I(6,2)+2*
     -   (-D1234I(7,3)-D1234I(8,3)+D1234I(9,3)+D1234I(10,3)))*P(11)
       F(62)=DCMPLX(FR(62),FI(62))
       FR(63) = 8*(D1234R(2,3)-D1234R(6,3)-D1234R(7,3)+D1234R(9,3)+2*(
     -   -D1234R(8,3)+D1234R(10,3)))*P(11)
       FI(63) = 8*(D1234I(2,3)-D1234I(6,3)-D1234I(7,3)+D1234I(9,3)+2*(
     -   -D1234I(8,3)+D1234I(10,3)))*P(11)
       F(63)=DCMPLX(FR(63),FI(63))
       FR(64) = (C0124R+2*(C124R(1,1)-C124R(2,1))-4*D1234R(7,2)-8*D123
     -   4R(11,3))*P(11)+C0134R*P(12)+D1234R(3,1)*P(23)+D01234R*P(174
     -   )+D1234R(1,1)*P(201)
       FI(64) = (C0124I+2*(C124I(1,1)-C124I(2,1))-4*D1234I(7,2)-8*D123
     -   4I(11,3))*P(11)+C0134I*P(12)+D1234I(3,1)*P(23)+D01234I*P(174
     -   )+D1234I(1,1)*P(201)
       F(64)=DCMPLX(FR(64),FI(64))
       FR(65) = -8*(D1234R(3,2)+D1234R(7,3))*P(11)
       FI(65) = -8*(D1234I(3,2)+D1234I(7,3))*P(11)
       F(65)=DCMPLX(FR(65),FI(65))
       FR(66) = -8*(D1234R(1,1)+D1234R(1,3)-D1234R(3,1)+D1234R(3,2)-3*
     -   D1234R(5,2)+2*(D1234R(1,2)-D1234R(5,3))+D1234R(7,3))*P(11)
       FI(66) = -8*(D1234I(1,1)+D1234I(1,3)-D1234I(3,1)+D1234I(3,2)-3*
     -   D1234I(5,2)+2*(D1234I(1,2)-D1234I(5,3))+D1234I(7,3))*P(11)
       F(66)=DCMPLX(FR(66),FI(66))
       FR(67) = 8*(D1234R(3,1)-D1234R(3,2)+2*D1234R(5,2)+D1234R(5,3)-D
     -   1234R(7,3))*P(11)
       FI(67) = 8*(D1234I(3,1)-D1234I(3,2)+2*D1234I(5,2)+D1234I(5,3)-D
     -   1234I(7,3))*P(11)
       F(67)=DCMPLX(FR(67),FI(67))
       FR(68) = 8*(-D1234R(3,2)+D1234R(5,2)+D1234R(5,3)-D1234R(7,3))*P
     -   (11)
       FI(68) = 8*(-D1234I(3,2)+D1234I(5,2)+D1234I(5,3)-D1234I(7,3))*P
     -   (11)
       F(68)=DCMPLX(FR(68),FI(68))
       FR(69) = -4*(D1234R(2,1)+D1234R(2,2)-D1234R(3,1)-3*D1234R(6,2)+
     -   2*(D1234R(3,2)+D1234R(4,2)-D1234R(5,2)+D1234R(6,3)+D1234R(7,
     -   3)-2*D1234R(10,3)))*P(11)
       FI(69) = -4*(D1234I(2,1)+D1234I(2,2)-D1234I(3,1)-3*D1234I(6,2)+
     -   2*(D1234I(3,2)+D1234I(4,2)-D1234I(5,2)+D1234I(6,3)+D1234I(7,
     -   3)-2*D1234I(10,3)))*P(11)
       F(69)=DCMPLX(FR(69),FI(69))
       FR(70) = 4*(D1234R(3,1)+D1234R(5,2)-2*(D1234R(3,2)-D1234R(6,2)+
     -   D1234R(7,3)-D1234R(10,3)))*P(11)
       FI(70) = 4*(D1234I(3,1)+D1234I(5,2)-2*(D1234I(3,2)-D1234I(6,2)+
     -   D1234I(7,3)-D1234I(10,3)))*P(11)
       F(70)=DCMPLX(FR(70),FI(70))
       FR(71) = 4*(D1234R(5,2)+D1234R(6,2)-2*(D1234R(3,2)+D1234R(7,3)-
     -   D1234R(10,3)))*P(11)
       FI(71) = 4*(D1234I(5,2)+D1234I(6,2)-2*(D1234I(3,2)+D1234I(7,3)-
     -   D1234I(10,3)))*P(11)
       F(71)=DCMPLX(FR(71),FI(71))
       FR(72) = -4*(D1234R(1,1)+D1234R(1,2)-D1234R(3,1)-3*D1234R(5,2)+
     -   2*(D1234R(3,2)+D1234R(4,2)+D1234R(4,3)-D1234R(5,3)-D1234R(6,
     -   2)+D1234R(7,3)-D1234R(10,3)))*P(11)
       FI(72) = -4*(D1234I(1,1)+D1234I(1,2)-D1234I(3,1)-3*D1234I(5,2)+
     -   2*(D1234I(3,2)+D1234I(4,2)+D1234I(4,3)-D1234I(5,3)-D1234I(6,
     -   2)+D1234I(7,3)-D1234I(10,3)))*P(11)
       F(72)=DCMPLX(FR(72),FI(72))
       FR(73) = -4*(D1234R(1,1)+D1234R(1,2)+D1234R(2,1)+3*D1234R(4,2)-
     -   5*D1234R(5,2)-D1234R(6,2)-2*(D1234R(3,1)-D1234R(3,2)-D1234R(
     -   4,3)+D1234R(5,3)-D1234R(7,3)+D1234R(10,3)))*P(11)
       FI(73) = -4*(D1234I(1,1)+D1234I(1,2)+D1234I(2,1)+3*D1234I(4,2)-
     -   5*D1234I(5,2)-D1234I(6,2)-2*(D1234I(3,1)-D1234I(3,2)-D1234I(
     -   4,3)+D1234I(5,3)-D1234I(7,3)+D1234I(10,3)))*P(11)
       F(73)=DCMPLX(FR(73),FI(73))
       P(240) = p2sq*za1*za2*za3*za4-p4sq*za1*za2*za3*za4+s12*za1*za2*
     -   za3*za4-s23*za1*za2*za3*za4-p2sq*zb1*zb2*zb3*zb4+p4sq*zb1*zb
     -   2*zb3*zb4-s12*zb1*zb2*zb3*zb4+s23*zb1*zb2*zb3*zb4+2*P(199)
       FR(74) = 4*(D1234R(7,2)-2*D1234R(12,3))*P(11)+(C0134R+C0234R)*P
     -   (12)+D01234R*P(175)+D1234R(2,1)*P(201)+D1234R(3,1)*P(240)
       FI(74) = 4*(D1234I(7,2)-2*D1234I(12,3))*P(11)+(C0134I+C0234I)*P
     -   (12)+D01234I*P(175)+D1234I(2,1)*P(201)+D1234I(3,1)*P(240)
       F(74)=DCMPLX(FR(74),FI(74))
       FR(75) = -8*(D1234R(3,2)+D1234R(9,3))*P(11)
       FI(75) = -8*(D1234I(3,2)+D1234I(9,3))*P(11)
       F(75)=DCMPLX(FR(75),FI(75))
       FR(76) = -4*(D1234R(2,1)-D1234R(3,1)+3*(D1234R(4,2)-D1234R(6,2)
     -   )+2*(D1234R(3,2)+D1234R(4,3)-D1234R(5,2)+D1234R(9,3)-2*D1234
     -   R(10,3)))*P(11)
       FI(76) = -4*(D1234I(2,1)-D1234I(3,1)+3*(D1234I(4,2)-D1234I(6,2)
     -   )+2*(D1234I(3,2)+D1234I(4,3)-D1234I(5,2)+D1234I(9,3)-2*D1234
     -   I(10,3)))*P(11)
       F(76)=DCMPLX(FR(76),FI(76))
       FR(77) = 4*(D1234R(3,1)+D1234R(5,2)-2*(D1234R(3,2)-D1234R(6,2)+
     -   D1234R(9,3)-D1234R(10,3)))*P(11)
       FI(77) = 4*(D1234I(3,1)+D1234I(5,2)-2*(D1234I(3,2)-D1234I(6,2)+
     -   D1234I(9,3)-D1234I(10,3)))*P(11)
       F(77)=DCMPLX(FR(77),FI(77))
       FR(78) = 4*(D1234R(5,2)+D1234R(6,2)-2*(D1234R(3,2)+D1234R(9,3)-
     -   D1234R(10,3)))*P(11)
       FI(78) = 4*(D1234I(5,2)+D1234I(6,2)-2*(D1234I(3,2)+D1234I(9,3)-
     -   D1234I(10,3)))*P(11)
       F(78)=DCMPLX(FR(78),FI(78))
       FR(79) = -8*(D1234R(2,2)+D1234R(2,3)+D1234R(3,2)-2*(D1234R(6,2)
     -   +D1234R(8,3))+D1234R(9,3))*P(11)
       FI(79) = -8*(D1234I(2,2)+D1234I(2,3)+D1234I(3,2)-2*(D1234I(6,2)
     -   +D1234I(8,3))+D1234I(9,3))*P(11)
       F(79)=DCMPLX(FR(79),FI(79))
       FR(80) = 8*(-D1234R(3,2)+D1234R(6,2)+D1234R(8,3)-D1234R(9,3))*P
     -   (11)
       FI(80) = 8*(-D1234I(3,2)+D1234I(6,2)+D1234I(8,3)-D1234I(9,3))*P
     -   (11)
       F(80)=DCMPLX(FR(80),FI(80))
       FR(81) = -4*(D1234R(2,1)-D1234R(3,1)+D1234R(4,2)-D1234R(5,2)+2*
     -   (D1234R(2,2)+D1234R(3,2)-2*D1234R(6,2)+D1234R(6,3)-D1234R(8,
     -   3)+D1234R(9,3)-D1234R(10,3)))*P(11)
       FI(81) = -4*(D1234I(2,1)-D1234I(3,1)+D1234I(4,2)-D1234I(5,2)+2*
     -   (D1234I(2,2)+D1234I(3,2)-2*D1234I(6,2)+D1234I(6,3)-D1234I(8,
     -   3)+D1234I(9,3)-D1234I(10,3)))*P(11)
       F(81)=DCMPLX(FR(81),FI(81))
       FR(82) = -4*(D1234R(2,2)+D1234R(4,2)-D1234R(5,2)-3*D1234R(6,2)+
     -   2*(D1234R(3,2)+D1234R(6,3)-D1234R(8,3)+D1234R(9,3)-D1234R(10
     -   ,3)))*P(11)
       FI(82) = -4*(D1234I(2,2)+D1234I(4,2)-D1234I(5,2)-3*D1234I(6,2)+
     -   2*(D1234I(3,2)+D1234I(6,3)-D1234I(8,3)+D1234I(9,3)-D1234I(10
     -   ,3)))*P(11)
       F(82)=DCMPLX(FR(82),FI(82))
       FR(83) = -2*(C123R(1,1)-2*D1234R(7,2))*P(11)+(C0123R+C0234R)*P(
     -   12)+D1234R(3,1)*P(106)+D01234R*P(142)+D1234R(1,1)*P(143)+D12
     -   34R(2,1)*P(181)
       FI(83) = -2*(C123I(1,1)-2*D1234I(7,2))*P(11)+(C0123I+C0234I)*P(
     -   12)+D1234I(3,1)*P(106)+D01234I*P(142)+D1234I(1,1)*P(143)+D12
     -   34I(2,1)*P(181)
       F(83)=DCMPLX(FR(83),FI(83))
       FR(84) = 4*(D1234R(1,1)+D1234R(1,2)-D1234R(2,1)-D1234R(4,2)-D12
     -   34R(5,2)+D1234R(6,2))*P(11)
       FI(84) = 4*(D1234I(1,1)+D1234I(1,2)-D1234I(2,1)-D1234I(4,2)-D12
     -   34I(5,2)+D1234I(6,2))*P(11)
       F(84)=DCMPLX(FR(84),FI(84))
       FR(85) = 4*(-D1234R(2,2)+D1234R(4,2)-D1234R(5,2)+D1234R(6,2))*P
     -   (11)
       FI(85) = 4*(-D1234I(2,2)+D1234I(4,2)-D1234I(5,2)+D1234I(6,2))*P
     -   (11)
       F(85)=DCMPLX(FR(85),FI(85))
       FR(86) = -4*(D1234R(1,1)+D1234R(1,2)-D1234R(3,1)+D1234R(3,2)-2*
     -   D1234R(5,2))*P(11)
       FI(86) = -4*(D1234I(1,1)+D1234I(1,2)-D1234I(3,1)+D1234I(3,2)-2*
     -   D1234I(5,2))*P(11)
       F(86)=DCMPLX(FR(86),FI(86))
       FR(87) = 4*(-D1234R(3,2)+D1234R(6,2))*P(11)
       FI(87) = 4*(-D1234I(3,2)+D1234I(6,2))*P(11)
       F(87)=DCMPLX(FR(87),FI(87))
       FR(88) = 4*(-D1234R(2,1)+D1234R(3,1)-D1234R(3,2)-D1234R(4,2)+D1
     -   234R(5,2)+D1234R(6,2))*P(11)
       FI(88) = 4*(-D1234I(2,1)+D1234I(3,1)-D1234I(3,2)-D1234I(4,2)+D1
     -   234I(5,2)+D1234I(6,2))*P(11)
       F(88)=DCMPLX(FR(88),FI(88))
       FR(89) = (C0134R-2*(C124R(1,1)-C124R(2,1))+4*D1234R(7,2)+8*(D12
     -   34R(11,3)-D1234R(13,3)))*P(11)+C0124R*P(12)+D1234R(3,1)*P(28
     -   )+D01234R*P(134)+D1234R(1,1)*P(135)
       FI(89) = (C0134I-2*(C124I(1,1)-C124I(2,1))+4*D1234I(7,2)+8*(D12
     -   34I(11,3)-D1234I(13,3)))*P(11)+C0124I*P(12)+D1234I(3,1)*P(28
     -   )+D01234I*P(134)+D1234I(1,1)*P(135)
       F(89)=DCMPLX(FR(89),FI(89))
       FR(90) = 8*(D1234R(1,1)+D1234R(1,3)-D1234R(3,1)+2*(D1234R(1,2)+
     -   D1234R(3,2))-D1234R(3,3)-4*D1234R(5,2)-3*(D1234R(5,3)-D1234R
     -   (7,3)))*P(11)
       FI(90) = 8*(D1234I(1,1)+D1234I(1,3)-D1234I(3,1)+2*(D1234I(1,2)+
     -   D1234I(3,2))-D1234I(3,3)-4*D1234I(5,2)-3*(D1234I(5,3)-D1234I
     -   (7,3)))*P(11)
       F(90)=DCMPLX(FR(90),FI(90))
       FR(91) = 8*(D1234R(3,2)-D1234R(3,3)-D1234R(5,2)-D1234R(5,3)+2*D
     -   1234R(7,3))*P(11)
       FI(91) = 8*(D1234I(3,2)-D1234I(3,3)-D1234I(5,2)-D1234I(5,3)+2*D
     -   1234I(7,3))*P(11)
       F(91)=DCMPLX(FR(91),FI(91))
       FR(92) = 8*(-D1234R(3,3)+D1234R(7,3))*P(11)
       FI(92) = 8*(-D1234I(3,3)+D1234I(7,3))*P(11)
       F(92)=DCMPLX(FR(92),FI(92))
       FR(93) = 4*(D1234R(2,1)+D1234R(2,2)-D1234R(3,1)+3*D1234R(3,2)+2
     -   *(-D1234R(3,3)+D1234R(4,2)-D1234R(5,2)+D1234R(6,3)+D1234R(7,
     -   3)-D1234R(8,3)-2*(D1234R(6,2)-D1234R(9,3)+D1234R(10,3))))*P(
     -   11)
       FI(93) = 4*(D1234I(2,1)+D1234I(2,2)-D1234I(3,1)+3*D1234I(3,2)+2
     -   *(-D1234I(3,3)+D1234I(4,2)-D1234I(5,2)+D1234I(6,3)+D1234I(7,
     -   3)-D1234I(8,3)-2*(D1234I(6,2)-D1234I(9,3)+D1234I(10,3))))*P(
     -   11)
       F(93)=DCMPLX(FR(93),FI(93))
       FR(94) = 4*(D1234R(1,1)+D1234R(1,2)+D1234R(2,1)-5*D1234R(5,2)+3
     -   *(D1234R(4,2)-D1234R(6,2))-2*(D1234R(3,1)+D1234R(3,3)-D1234R
     -   (4,3)+D1234R(5,3)-D1234R(9,3))+4*(D1234R(3,2)+D1234R(7,3)-D1
     -   234R(10,3)))*P(11)
       FI(94) = 4*(D1234I(1,1)+D1234I(1,2)+D1234I(2,1)-5*D1234I(5,2)+3
     -   *(D1234I(4,2)-D1234I(6,2))-2*(D1234I(3,1)+D1234I(3,3)-D1234I
     -   (4,3)+D1234I(5,3)-D1234I(9,3))+4*(D1234I(3,2)+D1234I(7,3)-D1
     -   234I(10,3)))*P(11)
       F(94)=DCMPLX(FR(94),FI(94))
       FR(95) = 4*(D1234R(1,1)+D1234R(1,2)-D1234R(3,1)+3*D1234R(3,2)+2
     -   *(-D1234R(3,3)+D1234R(4,2)+D1234R(4,3)-D1234R(5,3)-D1234R(6,
     -   2)+D1234R(9,3)-2*(D1234R(5,2)-D1234R(7,3)+D1234R(10,3))))*P(
     -   11)
       FI(95) = 4*(D1234I(1,1)+D1234I(1,2)-D1234I(3,1)+3*D1234I(3,2)+2
     -   *(-D1234I(3,3)+D1234I(4,2)+D1234I(4,3)-D1234I(5,3)-D1234I(6,
     -   2)+D1234I(9,3)-2*(D1234I(5,2)-D1234I(7,3)+D1234I(10,3))))*P(
     -   11)
       F(95)=DCMPLX(FR(95),FI(95))
       FR(96) = 4*(-D1234R(5,2)-D1234R(6,2)+2*(D1234R(3,2)-D1234R(3,3)
     -   +D1234R(7,3)+D1234R(9,3)-D1234R(10,3)))*P(11)
       FI(96) = 4*(-D1234I(5,2)-D1234I(6,2)+2*(D1234I(3,2)-D1234I(3,3)
     -   +D1234I(7,3)+D1234I(9,3)-D1234I(10,3)))*P(11)
       F(96)=DCMPLX(FR(96),FI(96))
       FR(97) = 4*(D1234R(3,2)-D1234R(5,2)-2*(D1234R(3,3)-D1234R(7,3)-
     -   D1234R(9,3)+D1234R(10,3)))*P(11)
       FI(97) = 4*(D1234I(3,2)-D1234I(5,2)-2*(D1234I(3,3)-D1234I(7,3)-
     -   D1234I(9,3)+D1234I(10,3)))*P(11)
       F(97)=DCMPLX(FR(97),FI(97))
       FR(98) = (C0134R+C0234R-4*D1234R(7,2)+8*(D1234R(12,3)-D1234R(13
     -   ,3)))*P(11)+D1234R(2,1)*P(135)+D01234R*P(154)+D1234R(3,1)*P(
     -   200)
       FI(98) = (C0134I+C0234I-4*D1234I(7,2)+8*(D1234I(12,3)-D1234I(13
     -   ,3)))*P(11)+D1234I(2,1)*P(135)+D01234I*P(154)+D1234I(3,1)*P(
     -   200)
       F(98)=DCMPLX(FR(98),FI(98))
       FR(99) = 8*(D1234R(2,2)+D1234R(2,3)+D1234R(3,2)-D1234R(3,3)-2*D
     -   1234R(6,2)-3*(D1234R(8,3)-D1234R(9,3)))*P(11)
       FI(99) = 8*(D1234I(2,2)+D1234I(2,3)+D1234I(3,2)-D1234I(3,3)-2*D
     -   1234I(6,2)-3*(D1234I(8,3)-D1234I(9,3)))*P(11)
       F(99)=DCMPLX(FR(99),FI(99))
       FR(100) = -8*(D1234R(3,3)+D1234R(8,3)-2*D1234R(9,3))*P(11)
       FI(100) = -8*(D1234I(3,3)+D1234I(8,3)-2*D1234I(9,3))*P(11)
       F(100)=DCMPLX(FR(100),FI(100))
       FR(101) = 8*(D1234R(3,2)-D1234R(3,3)-D1234R(6,2)-D1234R(8,3)+2*
     -   D1234R(9,3))*P(11)
       FI(101) = 8*(D1234I(3,2)-D1234I(3,3)-D1234I(6,2)-D1234I(8,3)+2*
     -   D1234I(9,3))*P(11)
       F(101)=DCMPLX(FR(101),FI(101))
       FR(102) = 4*(D1234R(2,1)-D1234R(3,1)+3*D1234R(3,2)+D1234R(4,2)-
     -   D1234R(5,2)-5*D1234R(6,2)+2*(D1234R(2,2)-D1234R(3,3)+D1234R(
     -   6,3)+D1234R(7,3)-D1234R(8,3))-4*(-D1234R(9,3)+D1234R(10,3)))
     -   *P(11)
       FI(102) = 4*(D1234I(2,1)-D1234I(3,1)+3*D1234I(3,2)+D1234I(4,2)-
     -   D1234I(5,2)-5*D1234I(6,2)+2*(D1234I(2,2)-D1234I(3,3)+D1234I(
     -   6,3)+D1234I(7,3)-D1234I(8,3))-4*(-D1234I(9,3)+D1234I(10,3)))
     -   *P(11)
       F(102)=DCMPLX(FR(102),FI(102))
       FR(103) = 4*(D1234R(2,2)+D1234R(4,2)-D1234R(5,2)-3*D1234R(6,2)+
     -   2*(D1234R(3,2)-D1234R(3,3)+D1234R(6,3)+D1234R(7,3)-D1234R(8,
     -   3))-4*(-D1234R(9,3)+D1234R(10,3)))*P(11)
       FI(103) = 4*(D1234I(2,2)+D1234I(4,2)-D1234I(5,2)-3*D1234I(6,2)+
     -   2*(D1234I(3,2)-D1234I(3,3)+D1234I(6,3)+D1234I(7,3)-D1234I(8,
     -   3))-4*(-D1234I(9,3)+D1234I(10,3)))*P(11)
       F(103)=DCMPLX(FR(103),FI(103))
       FR(104) = 8*(-D1234R(3,3)+D1234R(9,3))*P(11)
       FI(104) = 8*(-D1234I(3,3)+D1234I(9,3))*P(11)
       F(104)=DCMPLX(FR(104),FI(104))
       FR(105) = 4*(D1234R(2,1)-D1234R(3,1)+3*(D1234R(3,2)+D1234R(4,2)
     -   -D1234R(5,2)-D1234R(6,2))-2*(D1234R(3,3)-D1234R(4,3)+D1234R(
     -   5,3)-D1234R(9,3))-4*(-D1234R(7,3)+D1234R(10,3)))*P(11)
       FI(105) = 4*(D1234I(2,1)-D1234I(3,1)+3*(D1234I(3,2)+D1234I(4,2)
     -   -D1234I(5,2)-D1234I(6,2))-2*(D1234I(3,3)-D1234I(4,3)+D1234I(
     -   5,3)-D1234I(9,3))-4*(-D1234I(7,3)+D1234I(10,3)))*P(11)
       F(105)=DCMPLX(FR(105),FI(105))
       FR(106) = 4*(D1234R(3,2)+D1234R(5,2)-2*(D1234R(3,3)+D1234R(6,2)
     -   -D1234R(7,3)-D1234R(9,3)+D1234R(10,3)))*P(11)
       FI(106) = 4*(D1234I(3,2)+D1234I(5,2)-2*(D1234I(3,3)+D1234I(6,2)
     -   -D1234I(7,3)-D1234I(9,3)+D1234I(10,3)))*P(11)
       F(106)=DCMPLX(FR(106),FI(106))
       FR(107) = 4*(-D1234R(2,1)+D1234R(3,1)-D1234R(4,2)+D1234R(5,2))*
     -   P(11)
       FI(107) = 4*(-D1234I(2,1)+D1234I(3,1)-D1234I(4,2)+D1234I(5,2))*
     -   P(11)
       F(107)=DCMPLX(FR(107),FI(107))
       FR(108) = 4*(D1234R(5,2)-D1234R(6,2))*P(11)
       FI(108) = 4*(D1234I(5,2)-D1234I(6,2))*P(11)
       F(108)=DCMPLX(FR(108),FI(108))
       FR(109) = 4*(-D1234R(1,1)-D1234R(1,2)+D1234R(2,1)+D1234R(4,2)+D
     -   1234R(5,2)-D1234R(6,2))*P(11)
       FI(109) = 4*(-D1234I(1,1)-D1234I(1,2)+D1234I(2,1)+D1234I(4,2)+D
     -   1234I(5,2)-D1234I(6,2))*P(11)
       F(109)=DCMPLX(FR(109),FI(109))
       FR(110) = 4*(D1234R(2,2)-D1234R(4,2)+D1234R(5,2)-D1234R(6,2))*P
     -   (11)
       FI(110) = 4*(D1234I(2,2)-D1234I(4,2)+D1234I(5,2)-D1234I(6,2))*P
     -   (11)
       F(110)=DCMPLX(FR(110),FI(110))
       FR(111) = 4*(D1234R(2,1)-D1234R(3,1)+D1234R(3,2)+D1234R(4,2)-D1
     -   234R(5,2)-D1234R(6,2))*P(11)
       FI(111) = 4*(D1234I(2,1)-D1234I(3,1)+D1234I(3,2)+D1234I(4,2)-D1
     -   234I(5,2)-D1234I(6,2))*P(11)
       F(111)=DCMPLX(FR(111),FI(111))
       FR(112) = (C0123R+C0234R+2*C123R(1,1)-4*D1234R(7,2))*P(11)+D123
     -   4R(3,1)*P(104)+D01234R*P(161)+D1234R(2,1)*P(164)+D1234R(1,1)
     -   *P(182)
       FI(112) = (C0123I+C0234I+2*C123I(1,1)-4*D1234I(7,2))*P(11)+D123
     -   4I(3,1)*P(104)+D01234I*P(161)+D1234I(2,1)*P(164)+D1234I(1,1)
     -   *P(182)
       F(112)=DCMPLX(FR(112),FI(112))
       FR(113) = (C0134R+C0234R-12*D1234R(7,2))*P(11)+D1234R(2,1)*P(13
     -   5)+D01234R*P(154)+D1234R(3,1)*P(200)
       FI(113) = (C0134I+C0234I-12*D1234I(7,2))*P(11)+D1234I(2,1)*P(13
     -   5)+D01234I*P(154)+D1234I(3,1)*P(200)
       F(113)=DCMPLX(FR(113),FI(113))
       FR(114) = 12*D1234R(7,2)*P(11)+(C0134R+C0234R)*P(12)+D01234R*P(
     -   175)+D1234R(2,1)*P(201)+D1234R(3,1)*P(240)
       FI(114) = 12*D1234I(7,2)*P(11)+(C0134I+C0234I)*P(12)+D01234I*P(
     -   175)+D1234I(2,1)*P(201)+D1234I(3,1)*P(240)
       F(114)=DCMPLX(FR(114),FI(114))
       FR(115) = (C0124R+2*(C124R(1,1)-C124R(2,1))+4*D1234R(7,2))*P(11
     -   )+C0134R*P(12)+D1234R(3,1)*P(23)+D01234R*P(174)+D1234R(1,1)*
     -   P(201)
       FI(115) = (C0124I+2*(C124I(1,1)-C124I(2,1))+4*D1234I(7,2))*P(11
     -   )+C0134I*P(12)+D1234I(3,1)*P(23)+D01234I*P(174)+D1234I(1,1)*
     -   P(201)
       F(115)=DCMPLX(FR(115),FI(115))
       FR(116) = -4*(D1234R(3,1)+D1234R(5,2))*P(11)
       FI(116) = -4*(D1234I(3,1)+D1234I(5,2))*P(11)
       F(116)=DCMPLX(FR(116),FI(116))
       FR(117) = 4*(D1234R(1,1)+D1234R(1,2)-D1234R(3,1)-D1234R(5,2))*P
     -   (11)
       FI(117) = 4*(D1234I(1,1)+D1234I(1,2)-D1234I(3,1)-D1234I(5,2))*P
     -   (11)
       F(117)=DCMPLX(FR(117),FI(117))
       FR(118) = 4*(D1234R(2,1)-D1234R(3,1)+D1234R(4,2)-D1234R(5,2))*P
     -   (11)
       FI(118) = 4*(D1234I(2,1)-D1234I(3,1)+D1234I(4,2)-D1234I(5,2))*P
     -   (11)
       F(118)=DCMPLX(FR(118),FI(118))
       FR(119) = 4*(-D1234R(3,2)+D1234R(5,2))*P(11)
       FI(119) = 4*(-D1234I(3,2)+D1234I(5,2))*P(11)
       F(119)=DCMPLX(FR(119),FI(119))
       FR(120) = 4*(-D1234R(5,2)+D1234R(6,2))*P(11)
       FI(120) = 4*(-D1234I(5,2)+D1234I(6,2))*P(11)
       F(120)=DCMPLX(FR(120),FI(120))
       P(241) = -(p3sq*za2*za3*za4)-p4sq*za2*za3*za4+s23*za2*za3*za4+p
     -   4sq*za3*za4*zb2+p1sq*za4*zb2*zb3-s12*za4*zb2*zb3-p1sq*za2*za
     -   3*zb4+s12*za2*za3*zb4-p4sq*za2*zb3*zb4+p3sq*zb2*zb3*zb4+p4sq
     -   *zb2*zb3*zb4-s23*zb2*zb3*zb4+2*M**2*P(2)*P(4)*P(8)
       P(242) = -(s12*za2*za3*za4)-s23*za2*za3*za4+s23*za3*za4*zb2+s12
     -   *za4*zb2*zb3+s12*za2*za3*zb4+s23*za2*zb3*zb4-s12*zb2*zb3*zb4
     -   -s23*zb2*zb3*zb4+p2sq*P(2)*P(6)+p4sq*P(8)*P(97)
       P(243) = p1sq*P(4)*P(61)+P(5)*P(242)
       P(244) = -(P(11)*P(117))+M**2*P(243)
       P(245) = -(p1sq*P(4)*P(61))-P(5)*P(242)
       P(246) = P(11)*P(117)+M**2*P(245)
       P(247) = p1sq+2*p4sq-s12+3*P(138)
       P(248) = 2*M**2*P(128)-P(11)*P(247)
       P(249) = M**2*P(2)*P(5)*P(6)+p4sq*P(12)
       P(250) = M**2*P(8)*P(110)-P(11)*P(138)
       FR(121) = P(11)*(B013R+B024R-4*(C123R(4,2)+C124R(4,2))+C123R(1,
     -   1)*P(130)+(C124R(1,1)-C124R(2,1))*P(132))+M**2*P(25)*(C0124R
     -   *P(2)*P(35)+C0123R*P(8)*P(111)+D01234R*P(241))+D1234R(1,1)*P
     -   (244)+D1234R(3,1)*P(246)-2*D1234R(7,2)*P(248)+C0134R*P(249)+
     -   C0234R*P(250)
       FI(121) = P(11)*(B013I+B024I-4*(C123I(4,2)+C124I(4,2))+C123I(1,
     -   1)*P(130)+(C124I(1,1)-C124I(2,1))*P(132))+M**2*P(25)*(C0124I
     -   *P(2)*P(35)+C0123I*P(8)*P(111)+D01234I*P(241))+D1234I(1,1)*P
     -   (244)+D1234I(3,1)*P(246)-2*D1234I(7,2)*P(248)+C0134I*P(249)+
     -   C0234I*P(250)
       F(121)=DCMPLX(FR(121),FI(121))
       P(251) = p4sq*za2*za3*za4-s23*za2*za3*za4-p1sq*za4*zb2*zb3+p1sq
     -   *za2*za3*zb4-p4sq*zb2*zb3*zb4+s23*zb2*zb3*zb4
       P(252) = p4sq-s23
       P(253) = p1sq+3*P(252)
       P(254) = -(M**2*P(2)*P(5)*P(6))+s23*P(11)
       FR(122) = B023R*P(12)+D1234R(2,1)*P(244)+D1234R(3,1)*P(246)+C01
     -   34R*P(249)-D01234R*M**2*P(25)*P(251)+P(11)*(B013R+C123R(2,1)
     -   *P(130)+2*D1234R(7,2)*P(253))+C0234R*P(254)
       FI(122) = B023I*P(12)+D1234I(2,1)*P(244)+D1234I(3,1)*P(246)+C01
     -   34I*P(249)-D01234I*M**2*P(25)*P(251)+P(11)*(B013I+C123I(2,1)
     -   *P(130)+2*D1234I(7,2)*P(253))+C0234I*P(254)
       F(122)=DCMPLX(FR(122),FI(122))
       P(255) = p2sq*za2*za3*za4-s12*za2*za3*za4+p1sq*za3*za4*zb2-p1sq
     -   *za2*zb3*zb4-p2sq*zb2*zb3*zb4+s12*zb2*zb3*zb4
       P(256) = p2sq-s12
       P(257) = p1sq+3*P(256)
       P(258) = s12*P(12)+M**2*P(5)*P(8)*P(97)
       P(259) = p2sq*P(11)-M**2*P(5)*P(8)*P(97)
       FR(123) = B024R*P(12)+D1234R(3,1)*P(246)+D01234R*M**2*P(25)*P(2
     -   55)+P(11)*(B014R-C124R(2,1)*P(132)-2*D1234R(7,2)*P(257))+C01
     -   34R*P(258)+C0234R*P(259)
       FI(123) = B024I*P(12)+D1234I(3,1)*P(246)+D01234I*M**2*P(25)*P(2
     -   55)+P(11)*(B014I-C124I(2,1)*P(132)-2*D1234I(7,2)*P(257))+C01
     -   34I*P(258)+C0234I*P(259)
       F(123)=DCMPLX(FR(123),FI(123))
       P(260) = p4sq*P(11)-2*M**2*P(6)*P(115)
       P(261) = p1sq-p3sq-s12+s23
       P(262) = 2*M**2*P(4)*P(8)*P(19)+P(11)*P(261)
       P(263) = p3sq-p4sq+s12
       P(264) = 2*M**2*P(4)*P(8)*P(19)-P(11)*P(263)
       P(265) = p3sq+s12
       P(266) = p1sq+p4sq+s23-2*P(265)
       P(267) = -4*M**2*P(4)*P(8)*P(19)-P(11)*P(266)
       FR(124) = -2*(P(11)*(C123R(2,1)+C234R(1,1)-C234R(2,1)-8*(D1234R
     -   (12,3)-D1234R(13,3))+(D1234R(4,2)-D1234R(5,2))*P(162))-D1234
     -   R(2,1)*P(260)+D1234R(3,1)*P(260)-D1234R(2,2)*P(262)-D1234R(3
     -   ,2)*P(264)-D1234R(6,2)*P(267))
       FI(124) = -2*(P(11)*(C123I(2,1)+C234I(1,1)-C234I(2,1)-8*(D1234I
     -   (12,3)-D1234I(13,3))+(D1234I(4,2)-D1234I(5,2))*P(162))-D1234
     -   I(2,1)*P(260)+D1234I(3,1)*P(260)-D1234I(2,2)*P(262)-D1234I(3
     -   ,2)*P(264)-D1234I(6,2)*P(267))
       F(124)=DCMPLX(FR(124),FI(124))
       P(268) = P(11)*P(62)+M**2*P(5)*P(8)*P(97)
       P(269) = s12*P(12)+M**2*P(159)
       FR(125) = 2*(D01234R*M**2*P(25)*P(153)+P(11)*(C0234R+C234R(2,1)
     -   +p1sq*D1234R(1,1)-4*D1234R(7,2)-8*D1234R(13,3)+D1234R(5,2)*P
     -   (162))-D1234R(6,2)*P(262)+D1234R(3,2)*P(264)-D1234R(2,1)*P(2
     -   68)+D1234R(3,1)*P(269))
       FI(125) = 2*(D01234I*M**2*P(25)*P(153)+P(11)*(C0234I+C234I(2,1)
     -   +p1sq*D1234I(1,1)-4*D1234I(7,2)-8*D1234I(13,3)+D1234I(5,2)*P
     -   (162))-D1234I(6,2)*P(262)+D1234I(3,2)*P(264)-D1234I(2,1)*P(2
     -   68)+D1234I(3,1)*P(269))
       F(125)=DCMPLX(FR(125),FI(125))
       P(270) = p1sq+p3sq-2*p4sq+s12+s23
       P(271) = 2*M**2*P(4)*P(8)*P(19)-P(11)*P(270)
       FR(126) = 2*(P(11)*(C0234R-C123R(1,1)+C234R(2,1)+8*(D1234R(11,3
     -   )-D1234R(13,3))-D1234R(1,2)*P(162))+(D1234R(1,1)-D1234R(3,1)
     -   )*P(260)+(D1234R(4,2)-D1234R(6,2))*P(262)+D1234R(3,2)*P(264)
     -   -D1234R(5,2)*P(271))
       FI(126) = 2*(P(11)*(C0234I-C123I(1,1)+C234I(2,1)+8*(D1234I(11,3
     -   )-D1234I(13,3))-D1234I(1,2)*P(162))+(D1234I(1,1)-D1234I(3,1)
     -   )*P(260)+(D1234I(4,2)-D1234I(6,2))*P(262)+D1234I(3,2)*P(264)
     -   -D1234I(5,2)*P(271))
       F(126)=DCMPLX(FR(126),FI(126))
       P(272) = P(11)*P(62)+2*M**2*P(8)*P(116)
       FR(127) = 2*(D1234R(3,2)*P(28)+P(11)*(C124R(2,1)+2*C124R(2,2)+C
     -   134R(2,1)-4*D1234R(13,3)+D1234R(5,2)*P(130))+D1234R(3,1)*P(2
     -   72))
       FI(127) = 2*(D1234I(3,2)*P(28)+P(11)*(C124I(2,1)+2*C124I(2,2)+C
     -   134I(2,1)-4*D1234I(13,3)+D1234I(5,2)*P(130))+D1234I(3,1)*P(2
     -   72))
       F(127)=DCMPLX(FR(127),FI(127))
       P(273) = p1sq+2*p2sq+p3sq-s12-s23
       P(274) = 2*M**2*P(4)*P(8)*P(19)-P(11)*P(273)
       FR(128) = -2*(-(D1234R(3,2)*P(28))+P(11)*(C124R(1,1)-C124R(2,1)
     -   -2*(C124R(2,2)-C124R(3,2))+C134R(1,1)-C134R(2,1)-4*(D1234R(1
     -   1,3)-D1234R(13,3))+D1234R(1,2)*P(130))+(D1234R(1,1)-D1234R(3
     -   ,1))*P(272)+D1234R(5,2)*P(274))
       FI(128) = -2*(-(D1234I(3,2)*P(28))+P(11)*(C124I(1,1)-C124I(2,1)
     -   -2*(C124I(2,2)-C124I(3,2))+C134I(1,1)-C134I(2,1)-4*(D1234I(1
     -   1,3)-D1234I(13,3))+D1234I(1,2)*P(130))+(D1234I(1,1)-D1234I(3
     -   ,1))*P(272)+D1234I(5,2)*P(274))
       F(128)=DCMPLX(FR(128),FI(128))
       P(275) = p1sq-p4sq
       P(276) = M**2*P(179)+P(11)*P(275)
       FR(129) = -2*(-(D1234R(6,2)*P(23))-D1234R(3,2)*P(28)+P(11)*(C02
     -   34R+C134R(1,1)-C134R(2,1)+p1sq*D1234R(1,1)-4*(D1234R(7,2)+D1
     -   234R(12,3)-D1234R(13,3))+(D1234R(4,2)-D1234R(5,2))*P(130))+M
     -   **2*P(25)*(D1234R(2,1)*P(8)*P(111)+D01234R*P(153))-D1234R(3,
     -   1)*P(276))
       FI(129) = -2*(-(D1234I(6,2)*P(23))-D1234I(3,2)*P(28)+P(11)*(C02
     -   34I+C134I(1,1)-C134I(2,1)+p1sq*D1234I(1,1)-4*(D1234I(7,2)+D1
     -   234I(12,3)-D1234I(13,3))+(D1234I(4,2)-D1234I(5,2))*P(130))+M
     -   **2*P(25)*(D1234I(2,1)*P(8)*P(111)+D01234I*P(153))-D1234I(3,
     -   1)*P(276))
       F(129)=DCMPLX(FR(129),FI(129))
       P(277) = p2sq+p4sq-s12-s23
       FR(130) = 2*P(11)*(C134R(2,1)-C234R(2,1)+p1sq*D1234R(3,1)+D1234
     -   R(6,2)*P(130)-D1234R(3,2)*P(277))
       FI(130) = 2*P(11)*(C134I(2,1)-C234I(2,1)+p1sq*D1234I(3,1)+D1234
     -   I(6,2)*P(130)-D1234I(3,2)*P(277))
       F(130)=DCMPLX(FR(130),FI(130))
       P(278) = p1sq*za1*za2*za3*za4+p4sq*za1*za2*za3*za4-s12*za1*za2*
     -   za3*za4-M**2*za1*za3*za4*zb2+M**2*za3*za4*zb1*zb2+M**2*za1*z
     -   a4*zb2*zb3-M**2*za4*zb1*zb2*zb3+M**2*za1*za2*za3*zb4-M**2*za
     -   2*za3*zb1*zb4-M**2*za1*za2*zb3*zb4+M**2*za2*zb1*zb3*zb4-p1sq
     -   *zb1*zb2*zb3*zb4-p4sq*zb1*zb2*zb3*zb4+s12*zb1*zb2*zb3*zb4
       FR(131) = 2*(D01234R*M**2*P(25)*P(153)-D1234R(2,1)*P(268)-P(11)
     -   *(C134R(1,1)-C134R(2,1)+C234R(2,1)+4*D1234R(7,2)+(D1234R(4,2
     -   )-D1234R(6,2))*P(130)+(D1234R(3,2)-D1234R(5,2))*P(277))+D123
     -   4R(3,1)*P(278))
       FI(131) = 2*(D01234I*M**2*P(25)*P(153)-D1234I(2,1)*P(268)-P(11)
     -   *(C134I(1,1)-C134I(2,1)+C234I(2,1)+4*D1234I(7,2)+(D1234I(4,2
     -   )-D1234I(6,2))*P(130)+(D1234I(3,2)-D1234I(5,2))*P(277))+D123
     -   4I(3,1)*P(278))
       F(131)=DCMPLX(FR(131),FI(131))
       P(279) = p1sq+p4sq-s23+2*P(256)
       FR(132) = -2*P(11)*(C134R(1,1)-C134R(2,1)-C234R(1,1)+C234R(2,1)
     -   +p1sq*D1234R(2,1)-p1sq*D1234R(3,1)+D1234R(2,2)*P(130)+D1234R
     -   (3,2)*P(277)-D1234R(6,2)*P(279))
       FI(132) = -2*P(11)*(C134I(1,1)-C134I(2,1)-C234I(1,1)+C234I(2,1)
     -   +p1sq*D1234I(2,1)-p1sq*D1234I(3,1)+D1234I(2,2)*P(130)+D1234I
     -   (3,2)*P(277)-D1234I(6,2)*P(279))
       F(132)=DCMPLX(FR(132),FI(132))
       P(280) = p4sq*P(11)+2*M**2*P(2)*P(46)
       P(281) = p3sq*za1*za2*za3*za4-p4sq*za1*za2*za3*za4-s12*za1*za2*
     -   za3*za4-p3sq*zb1*zb2*zb3*zb4+p4sq*zb1*zb2*zb3*zb4+s12*zb1*zb
     -   2*zb3*zb4-2*P(69)
       P(282) = p3sq-s12
       P(283) = P(11)*P(282)
       P(284) = -(p1sq*za1*za2*za3*za4)+p4sq*za1*za2*za3*za4+s23*za1*z
     -   a2*za3*za4+p1sq*zb1*zb2*zb3*zb4-p4sq*zb1*zb2*zb3*zb4-s23*zb1
     -   *zb2*zb3*zb4+4*P(69)-2*P(283)
       FR(133) = 2*(P(11)*(C123R(2,1)-C234R(1,1)+C234R(2,1)+4*(D1234R(
     -   12,3)-D1234R(13,3))+(-D1234R(4,2)+D1234R(5,2))*P(132))+D1234
     -   R(2,2)*P(181)+(-D1234R(2,1)+D1234R(3,1))*P(280)+D1234R(3,2)*
     -   P(281)+D1234R(6,2)*P(284))
       FI(133) = 2*(P(11)*(C123I(2,1)-C234I(1,1)+C234I(2,1)+4*(D1234I(
     -   12,3)-D1234I(13,3))+(-D1234I(4,2)+D1234I(5,2))*P(132))+D1234
     -   I(2,2)*P(181)+(-D1234I(2,1)+D1234I(3,1))*P(280)+D1234I(3,2)*
     -   P(281)+D1234I(6,2)*P(284))
       F(133)=DCMPLX(FR(133),FI(133))
       P(285) = s12*P(12)+M**2*P(192)
       FR(134) = 2*(P(11)*(C0234R+C234R(2,1)+p1sq*D1234R(1,1)-4*(D1234
     -   R(7,2)+D1234R(13,3))+D1234R(5,2)*P(132))+D01234R*M**2*P(25)*
     -   P(153)-D1234R(6,2)*P(181)-D1234R(2,1)*P(268)+D1234R(3,2)*P(2
     -   81)+D1234R(3,1)*P(285))
       FI(134) = 2*(P(11)*(C0234I+C234I(2,1)+p1sq*D1234I(1,1)-4*(D1234
     -   I(7,2)+D1234I(13,3))+D1234I(5,2)*P(132))+D01234I*M**2*P(25)*
     -   P(153)-D1234I(6,2)*P(181)-D1234I(2,1)*P(268)+D1234I(3,2)*P(2
     -   81)+D1234I(3,1)*P(285))
       F(134)=DCMPLX(FR(134),FI(134))
       P(286) = p4sq*za1*za2*za3*za4+M**2*za2*za3*za4*zb1-M**2*za3*za4
     -   *zb1*zb2-M**2*za1*za2*za4*zb3+M**2*za1*za4*zb2*zb3-M**2*za2*
     -   za3*zb1*zb4+M**2*za3*zb1*zb2*zb4+M**2*za1*za2*zb3*zb4-M**2*z
     -   a1*zb2*zb3*zb4-p4sq*zb1*zb2*zb3*zb4
       P(287) = p1sq*za1*za2*za3*za4-p3sq*za1*za2*za3*za4+s12*za1*za2*
     -   za3*za4-s23*za1*za2*za3*za4-p1sq*zb1*zb2*zb3*zb4+p3sq*zb1*zb
     -   2*zb3*zb4-s12*zb1*zb2*zb3*zb4+s23*zb1*zb2*zb3*zb4+2*P(286)
       FR(135) = 2*(P(11)*(C0234R+C123R(1,1)+C234R(2,1)+4*(D1234R(11,3
     -   )-D1234R(13,3))-D1234R(1,2)*P(132))+(D1234R(4,2)-D1234R(6,2)
     -   )*P(181)+(-D1234R(1,1)+D1234R(3,1))*P(280)+D1234R(3,2)*P(281
     -   )+D1234R(5,2)*P(287))
       FI(135) = 2*(P(11)*(C0234I+C123I(1,1)+C234I(2,1)+4*(D1234I(11,3
     -   )-D1234I(13,3))-D1234I(1,2)*P(132))+(D1234I(4,2)-D1234I(6,2)
     -   )*P(181)+(-D1234I(1,1)+D1234I(3,1))*P(280)+D1234I(3,2)*P(281
     -   )+D1234I(5,2)*P(287))
       F(135)=DCMPLX(FR(135),FI(135))
       P(288) = p1sq*za1*za2*za3*za4-s12*za1*za2*za3*za4-M**2*za2*za3*
     -   za4*zb1+M**2*za3*za4*zb1*zb2+M**2*za1*za2*za4*zb3-M**2*za1*z
     -   a4*zb2*zb3+M**2*za2*za3*zb1*zb4-M**2*za3*zb1*zb2*zb4-M**2*za
     -   1*za2*zb3*zb4+M**2*za1*zb2*zb3*zb4-p1sq*zb1*zb2*zb3*zb4+s12*
     -   zb1*zb2*zb3*zb4
       P(289) = p2sq*za1*za2*za3*za4+p3sq*za1*za2*za3*za4-s23*za1*za2*
     -   za3*za4-p2sq*zb1*zb2*zb3*zb4-p3sq*zb1*zb2*zb3*zb4+s23*zb1*zb
     -   2*zb3*zb4+2*P(288)
       FR(136) = 2*(P(11)*(C124R(2,1)-C134R(2,1)+4*D1234R(13,3)-D1234R
     -   (5,2)*P(130))+D1234R(3,1)*P(174)+D1234R(3,2)*P(289))
       FI(136) = 2*(P(11)*(C124I(2,1)-C134I(2,1)+4*D1234I(13,3)-D1234I
     -   (5,2)*P(130))+D1234I(3,1)*P(174)+D1234I(3,2)*P(289))
       F(136)=DCMPLX(FR(136),FI(136))
       P(290) = p2sq*za1*za2*za3*za4-M**2*za2*za3*za4*zb1+M**2*za3*za4
     -   *zb1*zb2+M**2*za1*za2*za4*zb3-M**2*za1*za4*zb2*zb3+M**2*za2*
     -   za3*zb1*zb4-M**2*za3*zb1*zb2*zb4-M**2*za1*za2*zb3*zb4+M**2*z
     -   a1*zb2*zb3*zb4-p2sq*zb1*zb2*zb3*zb4
       P(291) = P(11)*P(62)
       P(292) = P(139)+2*P(290)+3*P(291)
       FR(137) = 2*(P(11)*(C124R(1,1)-C124R(2,1)+C134R(1,1)-C134R(2,1)
     -   -4*(D1234R(11,3)-D1234R(13,3))+D1234R(1,2)*P(130))+(D1234R(1
     -   ,1)-D1234R(3,1))*P(272)+D1234R(3,2)*P(289)-D1234R(5,2)*P(292
     -   ))
       FI(137) = 2*(P(11)*(C124I(1,1)-C124I(2,1)+C134I(1,1)-C134I(2,1)
     -   -4*(D1234I(11,3)-D1234I(13,3))+D1234I(1,2)*P(130))+(D1234I(1
     -   ,1)-D1234I(3,1))*P(272)+D1234I(3,2)*P(289)-D1234I(5,2)*P(292
     -   ))
       F(137)=DCMPLX(FR(137),FI(137))
       P(293) = -(p2sq*za1*za2*za3*za4)-p3sq*za1*za2*za3*za4+s23*za1*z
     -   a2*za3*za4+p2sq*zb1*zb2*zb3*zb4+p3sq*zb1*zb2*zb3*zb4-s23*zb1
     -   *zb2*zb3*zb4-2*P(288)
       FR(138) = 2*(P(11)*(C0234R+C134R(1,1)-C134R(2,1)+p1sq*D1234R(1,
     -   1)-4*(D1234R(7,2)+D1234R(12,3)-D1234R(13,3))+(D1234R(4,2)-D1
     -   234R(5,2))*P(130))+M**2*P(25)*(D1234R(2,1)*P(8)*P(111)+D0123
     -   4R*P(153))-D1234R(3,1)*P(276)+D1234R(3,2)*P(289)+D1234R(6,2)
     -   *P(293))
       FI(138) = 2*(P(11)*(C0234I+C134I(1,1)-C134I(2,1)+p1sq*D1234I(1,
     -   1)-4*(D1234I(7,2)+D1234I(12,3)-D1234I(13,3))+(D1234I(4,2)-D1
     -   234I(5,2))*P(130))+M**2*P(25)*(D1234I(2,1)*P(8)*P(111)+D0123
     -   4I*P(153))-D1234I(3,1)*P(276)+D1234I(3,2)*P(289)+D1234I(6,2)
     -   *P(293))
       F(138)=DCMPLX(FR(138),FI(138))
       P(294) = 2*p1sq+p2sq+p4sq-s12-s23
       FR(139) = 2*P(11)*(-C134R(2,1)+C234R(2,1)+p1sq*D1234R(3,1)-D123
     -   4R(6,2)*P(130)+D1234R(3,2)*P(294))
       FI(139) = 2*P(11)*(-C134I(2,1)+C234I(2,1)+p1sq*D1234I(3,1)-D123
     -   4I(6,2)*P(130)+D1234I(3,2)*P(294))
       F(139)=DCMPLX(FR(139),FI(139))
       FR(140) = 2*(-(D01234R*M**2*P(25)*P(153))+D1234R(2,1)*P(268)-D1
     -   234R(3,1)*P(278)+P(11)*(C134R(1,1)-C134R(2,1)+C234R(2,1)+4*D
     -   1234R(7,2)+(D1234R(4,2)-D1234R(6,2))*P(130)+(D1234R(3,2)-D12
     -   34R(5,2))*P(294)))
       FI(140) = 2*(-(D01234I*M**2*P(25)*P(153))+D1234I(2,1)*P(268)-D1
     -   234I(3,1)*P(278)+P(11)*(C134I(1,1)-C134I(2,1)+C234I(2,1)+4*D
     -   1234I(7,2)+(D1234I(4,2)-D1234I(6,2))*P(130)+(D1234I(3,2)-D12
     -   34I(5,2))*P(294)))
       F(140)=DCMPLX(FR(140),FI(140))
       P(295) = 3*p1sq+p4sq-s23+2*P(256)
       FR(141) = 2*P(11)*(C134R(1,1)-C134R(2,1)-C234R(1,1)+C234R(2,1)+
     -   p1sq*D1234R(2,1)-p1sq*D1234R(3,1)+D1234R(2,2)*P(130)+D1234R(
     -   3,2)*P(294)-D1234R(6,2)*P(295))
       FI(141) = 2*P(11)*(C134I(1,1)-C134I(2,1)-C234I(1,1)+C234I(2,1)+
     -   p1sq*D1234I(2,1)-p1sq*D1234I(3,1)+D1234I(2,2)*P(130)+D1234I(
     -   3,2)*P(294)-D1234I(6,2)*P(295))
       F(141)=DCMPLX(FR(141),FI(141))
       P(296) = p4sq*P(11)+s12*P(12)-M**2*P(4)*P(5)*P(82)
       FR(142) = -2*((C0234R+p1sq*D1234R(1,1)-2*D1234R(7,2)+4*D1234R(1
     -   3,3))*P(11)+D01234R*M**2*P(25)*P(153)-D1234R(2,1)*P(268)+D12
     -   34R(3,1)*P(296))
       FI(142) = -2*((C0234I+p1sq*D1234I(1,1)-2*D1234I(7,2)+4*D1234I(1
     -   3,3))*P(11)+D01234I*M**2*P(25)*P(153)-D1234I(2,1)*P(268)+D12
     -   34I(3,1)*P(296))
       F(142)=DCMPLX(FR(142),FI(142))
       FR(143) = -8*(D1234R(3,2)+D1234R(3,3))*P(11)
       FI(143) = -8*(D1234I(3,2)+D1234I(3,3))*P(11)
       F(143)=DCMPLX(FR(143),FI(143))
       FR(144) = -8*(D1234R(3,3)+D1234R(5,3)-2*D1234R(7,3))*P(11)
       FI(144) = -8*(D1234I(3,3)+D1234I(5,3)-2*D1234I(7,3))*P(11)
       F(144)=DCMPLX(FR(144),FI(144))
       FR(145) = 8*(-D1234R(3,2)-D1234R(3,3)+D1234R(5,2)+D1234R(7,3))*
     -   P(11)
       FI(145) = 8*(-D1234I(3,2)-D1234I(3,3)+D1234I(5,2)+D1234I(7,3))*
     -   P(11)
       F(145)=DCMPLX(FR(145),FI(145))
       FR(146) = 8*(-D1234R(3,2)-D1234R(3,3)+D1234R(6,2)+D1234R(9,3))*
     -   P(11)
       FI(146) = 8*(-D1234I(3,2)-D1234I(3,3)+D1234I(6,2)+D1234I(9,3))*
     -   P(11)
       F(146)=DCMPLX(FR(146),FI(146))
       FR(147) = 8*(-D1234R(3,3)+D1234R(7,3)+D1234R(9,3)-D1234R(10,3))
     -   *P(11)
       FI(147) = 8*(-D1234I(3,3)+D1234I(7,3)+D1234I(9,3)-D1234I(10,3))
     -   *P(11)
       F(147)=DCMPLX(FR(147),FI(147))
       P(297) = 2*M**2*P(2)*P(3)*P(4)-P(11)*P(212)
       FR(148) = -2*(-(D1234R(3,2)*P(106))+P(11)*(C0234R+C234R(2,1)+p1
     -   sq*D1234R(1,1)-4*(D1234R(7,2)+D1234R(13,3))+D1234R(5,2)*P(13
     -   2))+D01234R*M**2*P(25)*P(153)-D1234R(2,1)*P(268)+D1234R(3,1)
     -   *P(285)+D1234R(6,2)*P(297))
       FI(148) = -2*(-(D1234I(3,2)*P(106))+P(11)*(C0234I+C234I(2,1)+p1
     -   sq*D1234I(1,1)-4*(D1234I(7,2)+D1234I(13,3))+D1234I(5,2)*P(13
     -   2))+D01234I*M**2*P(25)*P(153)-D1234I(2,1)*P(268)+D1234I(3,1)
     -   *P(285)+D1234I(6,2)*P(297))
       F(148)=DCMPLX(FR(148),FI(148))
       P(298) = p1sq-p3sq+s12-s23
       P(299) = 2*M**2*P(2)*P(3)*P(4)+P(11)*P(298)
       FR(149) = 2*(D1234R(3,2)*P(106)-P(11)*(C0234R+C123R(1,1)-2*(C12
     -   3R(2,1)+C123R(3,2))+C234R(2,1)+4*(D1234R(11,3)-D1234R(13,3))
     -   -D1234R(1,2)*P(132))+D1234R(3,1)*P(142)+2*D1234R(2,1)*P(249)
     -   +D1234R(1,1)*P(280)+(D1234R(4,2)-D1234R(6,2))*P(297)-D1234R(
     -   5,2)*P(299))
       FI(149) = 2*(D1234I(3,2)*P(106)-P(11)*(C0234I+C123I(1,1)-2*(C12
     -   3I(2,1)+C123I(3,2))+C234I(2,1)+4*(D1234I(11,3)-D1234I(13,3))
     -   -D1234I(1,2)*P(132))+D1234I(3,1)*P(142)+2*D1234I(2,1)*P(249)
     -   +D1234I(1,1)*P(280)+(D1234I(4,2)-D1234I(6,2))*P(297)-D1234I(
     -   5,2)*P(299))
       F(149)=DCMPLX(FR(149),FI(149))
       P(300) = p1sq+3*p4sq-s23+2*P(282)
       P(301) = 4*M**2*P(2)*P(3)*P(4)-P(11)*P(300)
       FR(150) = 2*(D1234R(3,2)*P(106)+P(11)*(C123R(2,1)+2*C123R(2,2)+
     -   C234R(1,1)-C234R(2,1)-4*(D1234R(12,3)-D1234R(13,3))+(D1234R(
     -   4,2)-D1234R(5,2))*P(132))+(-D1234R(2,1)+D1234R(3,1))*P(142)+
     -   D1234R(2,2)*P(297)-D1234R(6,2)*P(301))
       FI(150) = 2*(D1234I(3,2)*P(106)+P(11)*(C123I(2,1)+2*C123I(2,2)+
     -   C234I(1,1)-C234I(2,1)-4*(D1234I(12,3)-D1234I(13,3))+(D1234I(
     -   4,2)-D1234I(5,2))*P(132))+(-D1234I(2,1)+D1234I(3,1))*P(142)+
     -   D1234I(2,2)*P(297)-D1234I(6,2)*P(301))
       F(150)=DCMPLX(FR(150),FI(150))
       P(302) = P(11)*P(62)+2*M**2*P(67)*P(97)
       P(303) = p4sq-s12
       P(304) = p2sq+p3sq-s23+2*P(303)
       P(305) = 2*M**2*P(2)*P(3)*P(4)-P(11)*P(304)
       FR(151) = -2*(P(11)*(C124R(2,1)-C134R(2,1)+2*C234R(2,1)-8*D1234
     -   R(13,3)+D1234R(5,2)*P(226))+D1234R(3,1)*P(302)-D1234R(3,2)*P
     -   (305))
       FI(151) = -2*(P(11)*(C124I(2,1)-C134I(2,1)+2*C234I(2,1)-8*D1234
     -   I(13,3)+D1234I(5,2)*P(226))+D1234I(3,1)*P(302)-D1234I(3,2)*P
     -   (305))
       F(151)=DCMPLX(FR(151),FI(151))
       P(306) = za2*zb1+za1*zb2
       P(307) = p1sq*P(11)+s12*P(11)+2*M**2*P(97)*P(306)
       P(308) = p1sq-2*p4sq+s12
       P(309) = za2*za3*P(2)+za4*zb2*P(4)
       P(310) = za2*za3-za2*zb3-zb2*zb3
       P(311) = za4*zb2*zb3+zb4*P(310)
       P(312) = za1*P(309)+zb1*P(311)
       P(313) = P(11)*P(308)+2*M**2*P(312)
       P(314) = p2sq+p4sq
       P(315) = p1sq-p3sq+3*s12+s23-2*P(314)
       P(316) = -2*M**2*P(2)*P(3)*P(4)-P(11)*P(315)
       FR(152) = 2*(2*D01234R*M**2*P(25)*P(153)+P(11)*(-C124R(1,1)+C12
     -   4R(2,1)-C134R(1,1)+C134R(2,1)-2*C234R(2,1)-8*(D1234R(7,2)+D1
     -   234R(11,3)-D1234R(13,3))+D1234R(1,2)*P(226))+D1234R(3,2)*P(3
     -   05)+D1234R(1,1)*P(307)-D1234R(3,1)*P(313)+D1234R(5,2)*P(316)
     -   )
       FI(152) = 2*(2*D01234I*M**2*P(25)*P(153)+P(11)*(-C124I(1,1)+C12
     -   4I(2,1)-C134I(1,1)+C134I(2,1)-2*C234I(2,1)-8*(D1234I(7,2)+D1
     -   234I(11,3)-D1234I(13,3))+D1234I(1,2)*P(226))+D1234I(3,2)*P(3
     -   05)+D1234I(1,1)*P(307)-D1234I(3,1)*P(313)+D1234I(5,2)*P(316)
     -   )
       F(152)=DCMPLX(FR(152),FI(152))
       P(317) = M**2*P(238)+P(11)*P(275)
       P(318) = -2*M**2*P(2)*P(3)*P(4)+P(11)*P(304)
       FR(153) = 2*(M**2*P(25)*(D1234R(2,1)*P(91)*P(97)+D01234R*P(153)
     -   )+P(11)*(C0234R-C134R(1,1)+C134R(2,1)+2*(C234R(1,1)-C234R(2,
     -   1))+p1sq*D1234R(1,1)-4*D1234R(7,2)-8*(D1234R(12,3)-D1234R(13
     -   ,3))+(D1234R(4,2)-D1234R(5,2))*P(226))+D1234R(3,2)*P(305)-D1
     -   234R(3,1)*P(317)+D1234R(6,2)*P(318))
       FI(153) = 2*(M**2*P(25)*(D1234I(2,1)*P(91)*P(97)+D01234I*P(153)
     -   )+P(11)*(C0234I-C134I(1,1)+C134I(2,1)+2*(C234I(1,1)-C234I(2,
     -   1))+p1sq*D1234I(1,1)-4*D1234I(7,2)-8*(D1234I(12,3)-D1234I(13
     -   ,3))+(D1234I(4,2)-D1234I(5,2))*P(226))+D1234I(3,2)*P(305)-D1
     -   234I(3,1)*P(317)+D1234I(6,2)*P(318))
       F(153)=DCMPLX(FR(153),FI(153))
       FR(154) = -2*P(11)*(-C134R(2,1)+C234R(2,1)+p1sq*D1234R(3,1)+D12
     -   34R(6,2)*P(226)+D1234R(3,2)*P(277))
       FI(154) = -2*P(11)*(-C134I(2,1)+C234I(2,1)+p1sq*D1234I(3,1)+D12
     -   34I(6,2)*P(226)+D1234I(3,2)*P(277))
       F(154)=DCMPLX(FR(154),FI(154))
       P(319) = M**2*P(5)*P(8)*P(97)-P(11)*P(221)
       P(320) = -(p1sq*za1*za2*za3*za4)+p4sq*za1*za2*za3*za4-s12*za1*z
     -   a2*za3*za4-M**2*za1*za3*za4*zb2+M**2*za3*za4*zb1*zb2+M**2*za
     -   1*za4*zb2*zb3-M**2*za4*zb1*zb2*zb3+M**2*za1*za2*za3*zb4-M**2
     -   *za2*za3*zb1*zb4-M**2*za1*za2*zb3*zb4+M**2*za2*zb1*zb3*zb4+p
     -   1sq*zb1*zb2*zb3*zb4-p4sq*zb1*zb2*zb3*zb4+s12*zb1*zb2*zb3*zb4
       FR(155) = 2*(D01234R*M**2*P(25)*P(153)-P(11)*(C134R(1,1)-C134R(
     -   2,1)+C234R(2,1)+4*D1234R(7,2)+(-D1234R(4,2)+D1234R(6,2))*P(2
     -   26)+(D1234R(3,2)-D1234R(5,2))*P(277))-D1234R(2,1)*P(319)+D12
     -   34R(3,1)*P(320))
       FI(155) = 2*(D01234I*M**2*P(25)*P(153)-P(11)*(C134I(1,1)-C134I(
     -   2,1)+C234I(2,1)+4*D1234I(7,2)+(-D1234I(4,2)+D1234I(6,2))*P(2
     -   26)+(D1234I(3,2)-D1234I(5,2))*P(277))-D1234I(2,1)*P(319)+D12
     -   34I(3,1)*P(320))
       F(155)=DCMPLX(FR(155),FI(155))
       P(321) = p1sq-p4sq+s23-2*P(256)
       FR(156) = -2*P(11)*(C134R(1,1)-C134R(2,1)-C234R(1,1)+C234R(2,1)
     -   -p1sq*D1234R(2,1)+p1sq*D1234R(3,1)-D1234R(2,2)*P(226)+D1234R
     -   (3,2)*P(277)+D1234R(6,2)*P(321))
       FI(156) = -2*P(11)*(C134I(1,1)-C134I(2,1)-C234I(1,1)+C234I(2,1)
     -   -p1sq*D1234I(2,1)+p1sq*D1234I(3,1)-D1234I(2,2)*P(226)+D1234I
     -   (3,2)*P(277)+D1234I(6,2)*P(321))
       F(156)=DCMPLX(FR(156),FI(156))
       FR(157) = 2*((C0234R+p1sq*D1234R(1,1)-2*D1234R(7,2)+4*D1234R(12
     -   ,3))*P(11)+D01234R*M**2*P(25)*P(153)-D1234R(2,1)*P(268)+D123
     -   4R(3,1)*P(296))
       FI(157) = 2*((C0234I+p1sq*D1234I(1,1)-2*D1234I(7,2)+4*D1234I(12
     -   ,3))*P(11)+D01234I*M**2*P(25)*P(153)-D1234I(2,1)*P(268)+D123
     -   4I(3,1)*P(296))
       F(157)=DCMPLX(FR(157),FI(157))
       FR(158) = 8*(D1234R(3,2)+D1234R(9,3))*P(11)
       FI(158) = 8*(D1234I(3,2)+D1234I(9,3))*P(11)
       F(158)=DCMPLX(FR(158),FI(158))
       FR(159) = 8*(D1234R(3,2)-D1234R(5,2)+D1234R(9,3)-D1234R(10,3))*
     -   P(11)
       FI(159) = 8*(D1234I(3,2)-D1234I(5,2)+D1234I(9,3)-D1234I(10,3))*
     -   P(11)
       F(159)=DCMPLX(FR(159),FI(159))
       FR(160) = 8*(D1234R(3,2)-D1234R(6,2)+D1234R(9,3)-D1234R(10,3))*
     -   P(11)
       FI(160) = 8*(D1234I(3,2)-D1234I(6,2)+D1234I(9,3)-D1234I(10,3))*
     -   P(11)
       F(160)=DCMPLX(FR(160),FI(160))
       FR(161) = 8*(D1234R(3,2)+D1234R(4,2)+D1234R(4,3)-D1234R(5,2)-D1
     -   234R(6,2)+D1234R(9,3)-2*D1234R(10,3))*P(11)
       FI(161) = 8*(D1234I(3,2)+D1234I(4,2)+D1234I(4,3)-D1234I(5,2)-D1
     -   234I(6,2)+D1234I(9,3)-2*D1234I(10,3))*P(11)
       F(161)=DCMPLX(FR(161),FI(161))
       FR(162) = 8*(D1234R(3,2)-D1234R(6,2)-D1234R(8,3)+D1234R(9,3))*P
     -   (11)
       FI(162) = 8*(D1234I(3,2)-D1234I(6,2)-D1234I(8,3)+D1234I(9,3))*P
     -   (11)
       F(162)=DCMPLX(FR(162),FI(162))
       FR(163) = 8*(D1234R(3,2)+D1234R(4,2)-D1234R(5,2)-D1234R(6,2)+D1
     -   234R(6,3)-D1234R(8,3)+D1234R(9,3)-D1234R(10,3))*P(11)
       FI(163) = 8*(D1234I(3,2)+D1234I(4,2)-D1234I(5,2)-D1234I(6,2)+D1
     -   234I(6,3)-D1234I(8,3)+D1234I(9,3)-D1234I(10,3))*P(11)
       F(163)=DCMPLX(FR(163),FI(163))
       FR(164) = 8*(D1234R(2,2)+D1234R(3,2)-2*D1234R(6,2)+D1234R(6,3)-
     -   D1234R(8,3)+D1234R(9,3)-D1234R(10,3))*P(11)
       FI(164) = 8*(D1234I(2,2)+D1234I(3,2)-2*D1234I(6,2)+D1234I(6,3)-
     -   D1234I(8,3)+D1234I(9,3)-D1234I(10,3))*P(11)
       F(164)=DCMPLX(FR(164),FI(164))
       FR(165) = 8*(D1234R(2,2)+D1234R(2,3)+D1234R(3,2)-2*(D1234R(6,2)
     -   +D1234R(8,3))+D1234R(9,3))*P(11)
       FI(165) = 8*(D1234I(2,2)+D1234I(2,3)+D1234I(3,2)-2*(D1234I(6,2)
     -   +D1234I(8,3))+D1234I(9,3))*P(11)
       F(165)=DCMPLX(FR(165),FI(165))
       FR(166) = -2*((C0234R+p1sq*D1234R(1,1)-2*D1234R(7,2)+4*(D1234R(
     -   12,3)-D1234R(13,3)))*P(11)+D01234R*M**2*P(25)*P(153)-D1234R(
     -   2,1)*P(268)+D1234R(3,1)*P(296))
       FI(166) = -2*((C0234I+p1sq*D1234I(1,1)-2*D1234I(7,2)+4*(D1234I(
     -   12,3)-D1234I(13,3)))*P(11)+D01234I*M**2*P(25)*P(153)-D1234I(
     -   2,1)*P(268)+D1234I(3,1)*P(296))
       F(166)=DCMPLX(FR(166),FI(166))
       FR(167) = 8*(D1234R(3,3)-D1234R(5,2)+D1234R(6,2)-D1234R(7,3)-D1
     -   234R(9,3)+D1234R(10,3))*P(11)
       FI(167) = 8*(D1234I(3,3)-D1234I(5,2)+D1234I(6,2)-D1234I(7,3)-D1
     -   234I(9,3)+D1234I(10,3))*P(11)
       F(167)=DCMPLX(FR(167),FI(167))
       FR(168) = 8*(-D1234R(3,2)+D1234R(3,3)+D1234R(5,2)-D1234R(7,3)-D
     -   1234R(9,3)+D1234R(10,3))*P(11)
       FI(168) = 8*(-D1234I(3,2)+D1234I(3,3)+D1234I(5,2)-D1234I(7,3)-D
     -   1234I(9,3)+D1234I(10,3))*P(11)
       F(168)=DCMPLX(FR(168),FI(168))
       FR(169) = -8*(D1234R(3,2)-D1234R(3,3)+D1234R(4,2)+D1234R(4,3)-D
     -   1234R(5,2)-D1234R(5,3)-D1234R(6,2)+D1234R(9,3)-2*(-D1234R(7,
     -   3)+D1234R(10,3)))*P(11)
       FI(169) = -8*(D1234I(3,2)-D1234I(3,3)+D1234I(4,2)+D1234I(4,3)-D
     -   1234I(5,2)-D1234I(5,3)-D1234I(6,2)+D1234I(9,3)-2*(-D1234I(7,
     -   3)+D1234I(10,3)))*P(11)
       F(169)=DCMPLX(FR(169),FI(169))
       FR(170) = 8*(-D1234R(3,2)+D1234R(3,3)+D1234R(6,2)+D1234R(8,3)-2
     -   *D1234R(9,3))*P(11)
       FI(170) = 8*(-D1234I(3,2)+D1234I(3,3)+D1234I(6,2)+D1234I(8,3)-2
     -   *D1234I(9,3))*P(11)
       F(170)=DCMPLX(FR(170),FI(170))
       FR(171) = -8*(D1234R(2,2)+D1234R(3,2)-D1234R(3,3)+D1234R(6,3)+D
     -   1234R(7,3)-D1234R(8,3)-2*(D1234R(6,2)-D1234R(9,3)+D1234R(10,
     -   3)))*P(11)
       FI(171) = -8*(D1234I(2,2)+D1234I(3,2)-D1234I(3,3)+D1234I(6,3)+D
     -   1234I(7,3)-D1234I(8,3)-2*(D1234I(6,2)-D1234I(9,3)+D1234I(10,
     -   3)))*P(11)
       F(171)=DCMPLX(FR(171),FI(171))
       FR(172) = -8*(D1234R(3,2)-D1234R(3,3)+D1234R(4,2)-D1234R(5,2)-D
     -   1234R(6,2)+D1234R(6,3)+D1234R(7,3)-D1234R(8,3)-2*(-D1234R(9,
     -   3)+D1234R(10,3)))*P(11)
       FI(172) = -8*(D1234I(3,2)-D1234I(3,3)+D1234I(4,2)-D1234I(5,2)-D
     -   1234I(6,2)+D1234I(6,3)+D1234I(7,3)-D1234I(8,3)-2*(-D1234I(9,
     -   3)+D1234I(10,3)))*P(11)
       F(172)=DCMPLX(FR(172),FI(172))
       FR(173) = -8*(D1234R(2,2)+D1234R(2,3)+D1234R(3,2)-D1234R(3,3)-2
     -   *D1234R(6,2)-3*(D1234R(8,3)-D1234R(9,3)))*P(11)
       FI(173) = -8*(D1234I(2,2)+D1234I(2,3)+D1234I(3,2)-D1234I(3,3)-2
     -   *D1234I(6,2)-3*(D1234I(8,3)-D1234I(9,3)))*P(11)
       F(173)=DCMPLX(FR(173),FI(173))
       FR(174) = (B012R+B024R+2*(B12R(1,1)-B13R(1,1))-4*(C124R(4,2)-C1
     -   34R(4,2)))*P(11)+(B013R+B034R)*P(12)-M**2*P(8)*(D01234R*P(96
     -   )+(-D1234R(1,1)+D1234R(3,1))*P(98))+C0124R*P(99)+C0134R*P(10
     -   0)+C124R(1,1)*P(102)+C134R(1,1)*P(104)+C124R(2,1)*P(105)+C13
     -   4R(2,1)*P(106)
       FI(174) = (B012I+B024I+2*(B12I(1,1)-B13I(1,1))-4*(C124I(4,2)-C1
     -   34I(4,2)))*P(11)+(B013I+B034I)*P(12)-M**2*P(8)*(D01234I*P(96
     -   )+(-D1234I(1,1)+D1234I(3,1))*P(98))+C0124I*P(99)+C0134I*P(10
     -   0)+C124I(1,1)*P(102)+C134I(1,1)*P(104)+C124I(2,1)*P(105)+C13
     -   4I(2,1)*P(106)
       F(174)=DCMPLX(FR(174),FI(174))
       FR(175) = -2*(B13R(1,1)-2*C134R(4,2))*P(11)+(B013R+B034R)*P(12)
     -   +C0134R*P(100)+C134R(1,1)*P(104)+C134R(2,1)*P(106)-M**2*P(8)
     -   *(C0234R*P(4)*P(19)+(-D1234R(2,1)+D1234R(3,1))*P(98)+D01234R
     -   *P(107)-4*D1234R(7,2)*P(110)+C0123R*P(25)*P(111))
       FI(175) = -2*(B13I(1,1)-2*C134I(4,2))*P(11)+(B013I+B034I)*P(12)
     -   +C0134I*P(100)+C134I(1,1)*P(104)+C134I(2,1)*P(106)-M**2*P(8)
     -   *(C0234I*P(4)*P(19)+(-D1234I(2,1)+D1234I(3,1))*P(98)+D01234I
     -   *P(107)-4*D1234I(7,2)*P(110)+C0123I*P(25)*P(111))
       F(175)=DCMPLX(FR(175),FI(175))
       FR(176) = (B024R-8*(C124R(4,2)-C134R(4,2)))*P(11)+B034R*P(12)+C
     -   124R(2,1)*P(105)+C134R(2,1)*P(106)-M**2*P(8)*(D1234R(3,1)*P(
     -   98)+D01234R*P(112))+C0124R*P(113)+C0134R*P(114)
       FI(176) = (B024I-8*(C124I(4,2)-C134I(4,2)))*P(11)+B034I*P(12)+C
     -   124I(2,1)*P(105)+C134I(2,1)*P(106)-M**2*P(8)*(D1234I(3,1)*P(
     -   98)+D01234I*P(112))+C0124I*P(113)+C0134I*P(114)
       F(176)=DCMPLX(FR(176),FI(176))
       FR(177) = 2*M**2*P(8)*(2*(D1234R(3,2)-D1234R(5,2))*P(4)*P(19)+D
     -   01234R*P(95)-D1234R(1,1)*P(5)*P(97)+D1234R(3,1)*P(38)*P(115)
     -   )
       FI(177) = 2*M**2*P(8)*(2*(D1234I(3,2)-D1234I(5,2))*P(4)*P(19)+D
     -   01234I*P(95)-D1234I(1,1)*P(5)*P(97)+D1234I(3,1)*P(38)*P(115)
     -   )
       F(177)=DCMPLX(FR(177),FI(177))
       FR(178) = 4*M**2*P(8)*((D1234R(1,2)+D1234R(3,2)-2*D1234R(5,2))*
     -   P(4)*P(19)+(-D1234R(1,1)+D1234R(3,1))*P(94))
       FI(178) = 4*M**2*P(8)*((D1234I(1,2)+D1234I(3,2)-2*D1234I(5,2))*
     -   P(4)*P(19)+(-D1234I(1,1)+D1234I(3,1))*P(94))
       F(178)=DCMPLX(FR(178),FI(178))
       FR(179) = 4*M**2*P(8)*((D1234R(3,2)+D1234R(4,2)-D1234R(5,2)-D12
     -   34R(6,2))*P(4)*P(19)+(-D1234R(2,1)+D1234R(3,1))*P(94))
       FI(179) = 4*M**2*P(8)*((D1234I(3,2)+D1234I(4,2)-D1234I(5,2)-D12
     -   34I(6,2))*P(4)*P(19)+(-D1234I(2,1)+D1234I(3,1))*P(94))
       F(179)=DCMPLX(FR(179),FI(179))
       FR(180) = -2*(2*(-C134R(1,1)-C134R(1,2)+C134R(2,1)+C134R(3,2))*
     -   P(11)+M**2*P(8)*(D1234R(3,1)*P(4)*P(19)+D01234R*P(95)-D1234R
     -   (1,1)*P(5)*P(97)))
       FI(180) = -2*(2*(-C134I(1,1)-C134I(1,2)+C134I(2,1)+C134I(3,2))*
     -   P(11)+M**2*P(8)*(D1234I(3,1)*P(4)*P(19)+D01234I*P(95)-D1234I
     -   (1,1)*P(5)*P(97)))
       F(180)=DCMPLX(FR(180),FI(180))
       FR(181) = 4*(C124R(2,1)+C124R(3,2)-C134R(2,1)-C134R(3,2))*P(11)
       FI(181) = 4*(C124I(2,1)+C124I(3,2)-C134I(2,1)-C134I(3,2))*P(11)
       F(181)=DCMPLX(FR(181),FI(181))
       FR(182) = -4*(C124R(1,1)+C124R(1,2)-C124R(2,1)-C124R(3,2)-C134R
     -   (1,1)-C134R(1,2)+C134R(2,1)+C134R(3,2))*P(11)
       FI(182) = -4*(C124I(1,1)+C124I(1,2)-C124I(2,1)-C124I(3,2)-C134I
     -   (1,1)-C134I(1,2)+C134I(2,1)+C134I(3,2))*P(11)
       F(182)=DCMPLX(FR(182),FI(182))
       FR(183) = -4*((C134R(2,1)+C134R(3,2))*P(11)+M**2*P(8)*(D1234R(3
     -   ,2)*P(4)*P(19)+D1234R(3,1)*P(116)))
       FI(183) = -4*((C134I(2,1)+C134I(3,2))*P(11)+M**2*P(8)*(D1234I(3
     -   ,2)*P(4)*P(19)+D1234I(3,1)*P(116)))
       F(183)=DCMPLX(FR(183),FI(183))
       FR(184) = 2*(2*(C134R(1,1)+C134R(1,2)-C134R(2,1)-C134R(3,2))*P(
     -   11)+M**2*P(8)*(2*(-D1234R(3,2)+D1234R(5,2))*P(4)*P(19)+D0123
     -   4R*P(95)+D1234R(1,1)*P(25)*P(111)-D1234R(3,1)*P(38)*P(115)))
       FI(184) = 2*(2*(C134I(1,1)+C134I(1,2)-C134I(2,1)-C134I(3,2))*P(
     -   11)+M**2*P(8)*(2*(-D1234I(3,2)+D1234I(5,2))*P(4)*P(19)+D0123
     -   4I*P(95)+D1234I(1,1)*P(25)*P(111)-D1234I(3,1)*P(38)*P(115)))
       F(184)=DCMPLX(FR(184),FI(184))
       FR(185) = 4*((C134R(1,1)+C134R(1,2)-C134R(2,1)-C134R(3,2))*P(11
     -   )+M**2*P(8)*((-D1234R(3,2)+D1234R(6,2))*P(4)*P(19)+(D1234R(2
     -   ,1)-D1234R(3,1))*P(116)))
       FI(185) = 4*((C134I(1,1)+C134I(1,2)-C134I(2,1)-C134I(3,2))*P(11
     -   )+M**2*P(8)*((-D1234I(3,2)+D1234I(6,2))*P(4)*P(19)+(D1234I(2
     -   ,1)-D1234I(3,1))*P(116)))
       F(185)=DCMPLX(FR(185),FI(185))
       FR(186) = -4*M**2*(D1234R(1,1)+D1234R(1,2)-D1234R(3,1)+D1234R(3
     -   ,2)-2*D1234R(5,2))*P(2)*P(8)*P(37)
       FI(186) = -4*M**2*(D1234I(1,1)+D1234I(1,2)-D1234I(3,1)+D1234I(3
     -   ,2)-2*D1234I(5,2))*P(2)*P(8)*P(37)
       F(186)=DCMPLX(FR(186),FI(186))
       FR(187) = 2*M**2*P(8)*(D1234R(3,1)*P(4)*P(19)+2*(-D1234R(3,2)+D
     -   1234R(5,2))*P(2)*P(37)+D01234R*P(95)-D1234R(1,1)*P(5)*P(97))
       FI(187) = 2*M**2*P(8)*(D1234I(3,1)*P(4)*P(19)+2*(-D1234I(3,2)+D
     -   1234I(5,2))*P(2)*P(37)+D01234I*P(95)-D1234I(1,1)*P(5)*P(97))
       F(187)=DCMPLX(FR(187),FI(187))
       FR(188) = 4*M**2*(-D1234R(2,1)+D1234R(3,1)-D1234R(3,2)-D1234R(4
     -   ,2)+D1234R(5,2)+D1234R(6,2))*P(2)*P(8)*P(37)
       FI(188) = 4*M**2*(-D1234I(2,1)+D1234I(3,1)-D1234I(3,2)-D1234I(4
     -   ,2)+D1234I(5,2)+D1234I(6,2))*P(2)*P(8)*P(37)
       F(188)=DCMPLX(FR(188),FI(188))
       FR(189) = 2*(-2*(C134R(1,1)+C134R(1,2)-C134R(2,1)+C134R(2,2)-2*
     -   C134R(3,2))*P(11)+M**2*P(8)*(D1234R(3,1)*P(4)*P(19)+D01234R*
     -   P(95)-D1234R(1,1)*P(5)*P(97)))
       FI(189) = 2*(-2*(C134I(1,1)+C134I(1,2)-C134I(2,1)+C134I(2,2)-2*
     -   C134I(3,2))*P(11)+M**2*P(8)*(D1234I(3,1)*P(4)*P(19)+D01234I*
     -   P(95)-D1234I(1,1)*P(5)*P(97)))
       F(189)=DCMPLX(FR(189),FI(189))
       FR(190) = 4*(C124R(1,1)+C124R(1,2)-C124R(2,1)+C124R(2,2)-C134R(
     -   1,1)-C134R(1,2)+C134R(2,1)-C134R(2,2)-2*(C124R(3,2)-C134R(3,
     -   2)))*P(11)
       FI(190) = 4*(C124I(1,1)+C124I(1,2)-C124I(2,1)+C124I(2,2)-C134I(
     -   1,1)-C134I(1,2)+C134I(2,1)-C134I(2,2)-2*(C124I(3,2)-C134I(3,
     -   2)))*P(11)
       F(190)=DCMPLX(FR(190),FI(190))
       FR(191) = 4*(C124R(2,2)-C124R(3,2)-C134R(2,2)+C134R(3,2))*P(11)
       FI(191) = 4*(C124I(2,2)-C124I(3,2)-C134I(2,2)+C134I(3,2))*P(11)
       F(191)=DCMPLX(FR(191),FI(191))
       FR(192) = 4*((-C134R(2,2)+C134R(3,2))*P(11)+M**2*P(8)*(D1234R(3
     -   ,2)*P(2)*P(37)+D1234R(3,1)*P(95)))
       FI(192) = 4*((-C134I(2,2)+C134I(3,2))*P(11)+M**2*P(8)*(D1234I(3
     -   ,2)*P(2)*P(37)+D1234I(3,1)*P(95)))
       F(192)=DCMPLX(FR(192),FI(192))
       FR(193) = -2*(2*(C134R(1,1)+C134R(1,2)-C134R(2,1)+C134R(2,2)-2*
     -   C134R(3,2))*P(11)+M**2*P(8)*(2*(-D1234R(3,2)+D1234R(5,2))*P(
     -   2)*P(37)+D01234R*P(95)+D1234R(1,1)*P(25)*P(111)-D1234R(3,1)*
     -   P(38)*P(115)))
       FI(193) = -2*(2*(C134I(1,1)+C134I(1,2)-C134I(2,1)+C134I(2,2)-2*
     -   C134I(3,2))*P(11)+M**2*P(8)*(2*(-D1234I(3,2)+D1234I(5,2))*P(
     -   2)*P(37)+D01234I*P(95)+D1234I(1,1)*P(25)*P(111)-D1234I(3,1)*
     -   P(38)*P(115)))
       F(193)=DCMPLX(FR(193),FI(193))
       FR(194) = -4*((C134R(1,1)+C134R(1,2)-C134R(2,1)+C134R(2,2)-2*C1
     -   34R(3,2))*P(11)-M**2*P(8)*((D1234R(3,2)-D1234R(6,2))*P(2)*P(
     -   37)+(-D1234R(2,1)+D1234R(3,1))*P(116)))
       FI(194) = -4*((C134I(1,1)+C134I(1,2)-C134I(2,1)+C134I(2,2)-2*C1
     -   34I(3,2))*P(11)-M**2*P(8)*((D1234I(3,2)-D1234I(6,2))*P(2)*P(
     -   37)+(-D1234I(2,1)+D1234I(3,1))*P(116)))
       F(194)=DCMPLX(FR(194),FI(194))
       FR(195) = -2*M**2*P(8)*(D1234R(3,1)*P(4)*P(19)+D01234R*P(95)-D1
     -   234R(1,1)*P(5)*P(97))
       FI(195) = -2*M**2*P(8)*(D1234I(3,1)*P(4)*P(19)+D01234I*P(95)-D1
     -   234I(1,1)*P(5)*P(97))
       F(195)=DCMPLX(FR(195),FI(195))
       FR(196) = 2*M**2*P(8)*(D1234R(3,1)*P(4)*P(19)+D01234R*P(95)-D12
     -   34R(1,1)*P(5)*P(97))
       FI(196) = 2*M**2*P(8)*(D1234I(3,1)*P(4)*P(19)+D01234I*P(95)-D12
     -   34I(1,1)*P(5)*P(97))
       F(196)=DCMPLX(FR(196),FI(196))
       FR(197) = 4*(C123R(4,2)+C124R(4,2))*P(11)+(B013R+B024R)*P(12)-M
     -   **2*P(4)*(C0134R*P(2)*P(3)+C0234R*P(8)*P(19)+D01234R*P(56)+(
     -   D1234R(1,1)-D1234R(3,1))*P(57)-4*D1234R(7,2)*P(61))+C0123R*P
     -   (66)+C0124R*P(68)+C123R(1,1)*P(70)+C124R(1,1)*P(71)+C124R(2,
     -   1)*P(72)
       FI(197) = 4*(C123I(4,2)+C124I(4,2))*P(11)+(B013I+B024I)*P(12)-M
     -   **2*P(4)*(C0134I*P(2)*P(3)+C0234I*P(8)*P(19)+D01234I*P(56)+(
     -   D1234I(1,1)-D1234I(3,1))*P(57)-4*D1234I(7,2)*P(61))+C0123I*P
     -   (66)+C0124I*P(68)+C123I(1,1)*P(70)+C124I(1,1)*P(71)+C124I(2,
     -   1)*P(72)
       F(197)=DCMPLX(FR(197),FI(197))
       FR(198) = (B013R+B024R)*P(12)+P(11)*(8*C123R(4,2)-2*p3sq*C234R(
     -   2,1)+C234R(1,1)*P(22))+C123R(2,1)*P(70)-M**2*P(4)*(C0134R*P(
     -   2)*P(3)+(D1234R(2,1)-D1234R(3,1))*P(57)-4*D1234R(7,2)*P(61)+
     -   D01234R*P(73))+C0123R*P(77)+C0234R*P(79)
       FI(198) = (B013I+B024I)*P(12)+P(11)*(8*C123I(4,2)-2*p3sq*C234I(
     -   2,1)+C234I(1,1)*P(22))+C123I(2,1)*P(70)-M**2*P(4)*(C0134I*P(
     -   2)*P(3)+(D1234I(2,1)-D1234I(3,1))*P(57)-4*D1234I(7,2)*P(61)+
     -   D01234I*P(73))+C0123I*P(77)+C0234I*P(79)
       F(198)=DCMPLX(FR(198),FI(198))
       FR(199) = 8*C124R(4,2)*P(11)+(B014R+B024R)*P(12)+C124R(2,1)*P(7
     -   2)+C0124R*P(77)-M**2*P(4)*(C0134R*P(2)*P(3)+C0234R*P(8)*P(19
     -   )-D1234R(3,1)*P(57)-4*D1234R(7,2)*P(61)+D01234R*P(80))
       FI(199) = 8*C124I(4,2)*P(11)+(B014I+B024I)*P(12)+C124I(2,1)*P(7
     -   2)+C0124I*P(77)-M**2*P(4)*(C0134I*P(2)*P(3)+C0234I*P(8)*P(19
     -   )-D1234I(3,1)*P(57)-4*D1234I(7,2)*P(61)+D01234I*P(80))
       F(199)=DCMPLX(FR(199),FI(199))
       FR(200) = 2*M**2*P(4)*((D1234R(2,1)-2*(D1234R(5,2)-D1234R(6,2))
     -   )*P(8)*P(19)-2*D1234R(3,1)*P(55)-D01234R*P(81)+D1234R(1,1)*P
     -   (5)*P(82))
       FI(200) = 2*M**2*P(4)*((D1234I(2,1)-2*(D1234I(5,2)-D1234I(6,2))
     -   )*P(8)*P(19)-2*D1234I(3,1)*P(55)-D01234I*P(81)+D1234I(1,1)*P
     -   (5)*P(82))
       F(200)=DCMPLX(FR(200),FI(200))
       FR(201) = 4*((C123R(1,1)+C123R(1,2))*P(11)+M**2*P(4)*((D1234R(1
     -   ,2)-D1234R(4,2)-D1234R(5,2)+D1234R(6,2))*P(8)*P(19)+(D1234R(
     -   1,1)-D1234R(3,1))*P(55)))
       FI(201) = 4*((C123I(1,1)+C123I(1,2))*P(11)+M**2*P(4)*((D1234I(1
     -   ,2)-D1234I(4,2)-D1234I(5,2)+D1234I(6,2))*P(8)*P(19)+(D1234I(
     -   1,1)-D1234I(3,1))*P(55)))
       F(201)=DCMPLX(FR(201),FI(201))
       FR(202) = 4*((C123R(2,1)+C123R(3,2))*P(11)+M**2*P(4)*((-D1234R(
     -   2,2)+D1234R(4,2)-D1234R(5,2)+D1234R(6,2))*P(8)*P(19)+(D1234R
     -   (2,1)-D1234R(3,1))*P(55)))
       FI(202) = 4*((C123I(2,1)+C123I(3,2))*P(11)+M**2*P(4)*((-D1234I(
     -   2,2)+D1234I(4,2)-D1234I(5,2)+D1234I(6,2))*P(8)*P(19)+(D1234I
     -   (2,1)-D1234I(3,1))*P(55)))
       F(202)=DCMPLX(FR(202),FI(202))
       FR(203) = -4*((C124R(2,1)+C124R(3,2))*P(11)+M**2*(D1234R(3,1)+D
     -   1234R(5,2))*P(4)*P(8)*P(19))
       FI(203) = -4*((C124I(2,1)+C124I(3,2))*P(11)+M**2*(D1234I(3,1)+D
     -   1234I(5,2))*P(4)*P(8)*P(19))
       F(203)=DCMPLX(FR(203),FI(203))
       FR(204) = 4*((C124R(1,1)+C124R(1,2)-C124R(2,1)-C124R(3,2))*P(11
     -   )+M**2*(D1234R(1,1)+D1234R(1,2)-D1234R(3,1)-D1234R(5,2))*P(4
     -   )*P(8)*P(19))
       FI(204) = 4*((C124I(1,1)+C124I(1,2)-C124I(2,1)-C124I(3,2))*P(11
     -   )+M**2*(D1234I(1,1)+D1234I(1,2)-D1234I(3,1)-D1234I(5,2))*P(4
     -   )*P(8)*P(19))
       F(204)=DCMPLX(FR(204),FI(204))
       FR(205) = 2*M**2*P(4)*((D1234R(2,1)-2*(D1234R(3,1)-D1234R(4,2)+
     -   D1234R(5,2)))*P(8)*P(19)+D01234R*P(81)-D1234R(1,1)*P(5)*P(82
     -   ))
       FI(205) = 2*M**2*P(4)*((D1234I(2,1)-2*(D1234I(3,1)-D1234I(4,2)+
     -   D1234I(5,2)))*P(8)*P(19)+D01234I*P(81)-D1234I(1,1)*P(5)*P(82
     -   ))
       F(205)=DCMPLX(FR(205),FI(205))
       FR(206) = 4*M**2*P(4)*(-(D1234R(6,2)*P(8)*P(19))+D1234R(3,1)*P(
     -   83))
       FI(206) = 4*M**2*P(4)*(-(D1234I(6,2)*P(8)*P(19))+D1234I(3,1)*P(
     -   83))
       F(206)=DCMPLX(FR(206),FI(206))
       FR(207) = 2*M**2*P(4)*((D1234R(2,1)+2*(D1234R(4,2)-D1234R(6,2))
     -   )*P(8)*P(19)-D01234R*P(81)+2*D1234R(3,1)*P(83)-D1234R(1,1)*P
     -   (25)*P(84))
       FI(207) = 2*M**2*P(4)*((D1234I(2,1)+2*(D1234I(4,2)-D1234I(6,2))
     -   )*P(8)*P(19)-D01234I*P(81)+2*D1234I(3,1)*P(83)-D1234I(1,1)*P
     -   (25)*P(84))
       F(207)=DCMPLX(FR(207),FI(207))
       FR(208) = 4*M**2*P(4)*((D1234R(2,2)-D1234R(6,2))*P(8)*P(19)+(-D
     -   1234R(2,1)+D1234R(3,1))*P(83))
       FI(208) = 4*M**2*P(4)*((D1234I(2,2)-D1234I(6,2))*P(8)*P(19)+(-D
     -   1234I(2,1)+D1234I(3,1))*P(83))
       F(208)=DCMPLX(FR(208),FI(208))
       FR(209) = 2*M**2*P(4)*(2*(D1234R(5,2)-D1234R(6,2))*P(2)*P(3)+D1
     -   234R(2,1)*P(8)*P(19)-D01234R*P(81)+D1234R(1,1)*P(5)*P(82))
       FI(209) = 2*M**2*P(4)*(2*(D1234I(5,2)-D1234I(6,2))*P(2)*P(3)+D1
     -   234I(2,1)*P(8)*P(19)-D01234I*P(81)+D1234I(1,1)*P(5)*P(82))
       F(209)=DCMPLX(FR(209),FI(209))
       FR(210) = 4*(M**2*(-D1234R(1,1)-D1234R(1,2)+D1234R(2,1)+D1234R(
     -   4,2)+D1234R(5,2)-D1234R(6,2))*P(2)*P(3)*P(4)+(-C123R(1,1)-C1
     -   23R(1,2)+C123R(2,1)+C123R(3,2))*P(11))
       FI(210) = 4*(M**2*(-D1234I(1,1)-D1234I(1,2)+D1234I(2,1)+D1234I(
     -   4,2)+D1234I(5,2)-D1234I(6,2))*P(2)*P(3)*P(4)+(-C123I(1,1)-C1
     -   23I(1,2)+C123I(2,1)+C123I(3,2))*P(11))
       F(210)=DCMPLX(FR(210),FI(210))
       FR(211) = 4*(M**2*(D1234R(2,2)-D1234R(4,2)+D1234R(5,2)-D1234R(6
     -   ,2))*P(2)*P(3)*P(4)+(C123R(2,2)-C123R(3,2))*P(11))
       FI(211) = 4*(M**2*(D1234I(2,2)-D1234I(4,2)+D1234I(5,2)-D1234I(6
     -   ,2))*P(2)*P(3)*P(4)+(C123I(2,2)-C123I(3,2))*P(11))
       F(211)=DCMPLX(FR(211),FI(211))
       FR(212) = 2*M**2*P(4)*(2*(-D1234R(4,2)+D1234R(5,2))*P(2)*P(3)+(
     -   D01234R-2*D1234R(3,1))*P(81)-D1234R(1,1)*P(5)*P(82)+D1234R(2
     -   ,1)*P(87))
       FI(212) = 2*M**2*P(4)*(2*(-D1234I(4,2)+D1234I(5,2))*P(2)*P(3)+(
     -   D01234I-2*D1234I(3,1))*P(81)-D1234I(1,1)*P(5)*P(82)+D1234I(2
     -   ,1)*P(87))
       F(212)=DCMPLX(FR(212),FI(212))
       FR(213) = -4*((C124R(1,1)+C124R(1,2)-C124R(2,1)+C124R(2,2)-2*C1
     -   24R(3,2))*P(11)+M**2*P(4)*((D1234R(1,2)-D1234R(5,2))*P(2)*P(
     -   3)+(-D01234R+D1234R(3,1))*P(81)+D1234R(1,1)*P(89)))
       FI(213) = -4*((C124I(1,1)+C124I(1,2)-C124I(2,1)+C124I(2,2)-2*C1
     -   24I(3,2))*P(11)+M**2*P(4)*((D1234I(1,2)-D1234I(5,2))*P(2)*P(
     -   3)+(-D01234I+D1234I(3,1))*P(81)+D1234I(1,1)*P(89)))
       F(213)=DCMPLX(FR(213),FI(213))
       FR(214) = 4*((-C124R(2,2)+C124R(3,2))*P(11)+M**2*P(4)*(D1234R(5
     -   ,2)*P(2)*P(3)-D1234R(3,1)*P(81)))
       FI(214) = 4*((-C124I(2,2)+C124I(3,2))*P(11)+M**2*P(4)*(D1234I(5
     -   ,2)*P(2)*P(3)-D1234I(3,1)*P(81)))
       F(214)=DCMPLX(FR(214),FI(214))
       FR(215) = 4*M**2*P(4)*(D1234R(6,2)*P(2)*P(3)-D1234R(3,1)*P(81))
       FI(215) = 4*M**2*P(4)*(D1234I(6,2)*P(2)*P(3)-D1234I(3,1)*P(81))
       F(215)=DCMPLX(FR(215),FI(215))
       FR(216) = 2*M**2*P(4)*(2*(-D1234R(4,2)+D1234R(6,2))*P(2)*P(3)+(
     -   D01234R-2*D1234R(3,1))*P(81)+D1234R(1,1)*P(25)*P(84)-D1234R(
     -   2,1)*P(93))
       FI(216) = 2*M**2*P(4)*(2*(-D1234I(4,2)+D1234I(6,2))*P(2)*P(3)+(
     -   D01234I-2*D1234I(3,1))*P(81)+D1234I(1,1)*P(25)*P(84)-D1234I(
     -   2,1)*P(93))
       F(216)=DCMPLX(FR(216),FI(216))
       FR(217) = 4*M**2*P(4)*((-D1234R(2,2)+D1234R(6,2))*P(2)*P(3)+(D1
     -   234R(2,1)-D1234R(3,1))*P(81))
       FI(217) = 4*M**2*P(4)*((-D1234I(2,2)+D1234I(6,2))*P(2)*P(3)+(D1
     -   234I(2,1)-D1234I(3,1))*P(81))
       F(217)=DCMPLX(FR(217),FI(217))
       FR(218) = -2*M**2*P(4)*(D1234R(2,1)*P(8)*P(19)-D01234R*P(81)+D1
     -   234R(1,1)*P(5)*P(82))
       FI(218) = -2*M**2*P(4)*(D1234I(2,1)*P(8)*P(19)-D01234I*P(81)+D1
     -   234I(1,1)*P(5)*P(82))
       F(218)=DCMPLX(FR(218),FI(218))
       FR(219) = 2*M**2*P(4)*(D1234R(2,1)*P(8)*P(19)-D01234R*P(81)+D12
     -   34R(1,1)*P(5)*P(82))
       FI(219) = 2*M**2*P(4)*(D1234I(2,1)*P(8)*P(19)-D01234I*P(81)+D12
     -   34I(1,1)*P(5)*P(82))
       F(219)=DCMPLX(FR(219),FI(219))
       FR(220) = -(M**2*P(2)*(D01234R*P(1)+(D1234R(1,1)-D1234R(3,1))*P
     -   (7)))+(B013R-2*(B12R(1,1)+B24R(1,1)+2*(C123R(4,2)-C234R(4,2)
     -   )))*P(11)+(B012R+B024R+B034R)*P(12)+C0123R*P(13)+C0234R*P(18
     -   )+C123R(1,1)*P(21)+C234R(2,1)*P(23)
       FI(220) = -(M**2*P(2)*(D01234I*P(1)+(D1234I(1,1)-D1234I(3,1))*P
     -   (7)))+(B013I-2*(B12I(1,1)+B24I(1,1)+2*(C123I(4,2)-C234I(4,2)
     -   )))*P(11)+(B012I+B024I+B034I)*P(12)+C0123I*P(13)+C0234I*P(18
     -   )+C123I(1,1)*P(21)+C234I(2,1)*P(23)
       F(220)=DCMPLX(FR(220),FI(220))
       FR(221) = (B013R-8*(C123R(4,2)-C234R(4,2)))*P(11)+B034R*P(12)+C
     -   123R(2,1)*P(21)+C234R(2,1)*P(23)-M**2*P(2)*((D1234R(2,1)-D12
     -   34R(3,1))*P(7)+D01234R*P(24))+C0123R*P(26)+C0234R*P(27)+C234
     -   R(1,1)*P(28)
       FI(221) = (B013I-8*(C123I(4,2)-C234I(4,2)))*P(11)+B034I*P(12)+C
     -   123I(2,1)*P(21)+C234I(2,1)*P(23)-M**2*P(2)*((D1234I(2,1)-D12
     -   34I(3,1))*P(7)+D01234I*P(24))+C0123I*P(26)+C0234I*P(27)+C234
     -   I(1,1)*P(28)
       F(221)=DCMPLX(FR(221),FI(221))
       FR(222) = -2*(B24R(1,1)-2*C234R(4,2))*P(11)+(B024R+B034R)*P(12)
     -   +C0234R*P(18)+C234R(2,1)*P(23)-M**2*P(2)*(C0134R*P(3)*P(4)-D
     -   1234R(3,1)*P(7)+D01234R*P(29)-4*D1234R(7,2)*P(34)+C0124R*P(2
     -   5)*P(35))
       FI(222) = -2*(B24I(1,1)-2*C234I(4,2))*P(11)+(B024I+B034I)*P(12)
     -   +C0234I*P(18)+C234I(2,1)*P(23)-M**2*P(2)*(C0134I*P(3)*P(4)-D
     -   1234I(3,1)*P(7)+D01234I*P(29)-4*D1234I(7,2)*P(34)+C0124I*P(2
     -   5)*P(35))
       F(222)=DCMPLX(FR(222),FI(222))
       FR(223) = -2*(-2*(C234R(2,1)+C234R(2,2))*P(11)+M**2*P(2)*(-(D12
     -   34R(3,1)*P(3)*P(4))+D1234R(1,1)*P(5)*P(6)+D01234R*P(36)+D123
     -   4R(2,1)*P(8)*P(37)))
       FI(223) = -2*(-2*(C234I(2,1)+C234I(2,2))*P(11)+M**2*P(2)*(-(D12
     -   34I(3,1)*P(3)*P(4))+D1234I(1,1)*P(5)*P(6)+D01234I*P(36)+D123
     -   4I(2,1)*P(8)*P(37)))
       F(223)=DCMPLX(FR(223),FI(223))
       FR(224) = 4*(-C123R(1,1)-C123R(1,2)+C234R(2,1)+C234R(2,2))*P(11
     -   )
       FI(224) = 4*(-C123I(1,1)-C123I(1,2)+C234I(2,1)+C234I(2,2))*P(11
     -   )
       F(224)=DCMPLX(FR(224),FI(224))
       FR(225) = -4*(C123R(2,1)+C123R(3,2)-C234R(2,2)+C234R(3,2))*P(11
     -   )
       FI(225) = -4*(C123I(2,1)+C123I(3,2)-C234I(2,2)+C234I(3,2))*P(11
     -   )
       F(225)=DCMPLX(FR(225),FI(225))
       FR(226) = 4*M**2*(D1234R(3,2)-D1234R(5,2))*P(2)*P(8)*P(37)
       FI(226) = 4*M**2*(D1234I(3,2)-D1234I(5,2))*P(2)*P(8)*P(37)
       F(226)=DCMPLX(FR(226),FI(226))
       FR(227) = 4*M**2*(D1234R(1,1)+D1234R(1,2)-D1234R(3,1)+D1234R(3,
     -   2)-2*D1234R(5,2))*P(2)*P(8)*P(37)
       FI(227) = 4*M**2*(D1234I(1,1)+D1234I(1,2)-D1234I(3,1)+D1234I(3,
     -   2)-2*D1234I(5,2))*P(2)*P(8)*P(37)
       F(227)=DCMPLX(FR(227),FI(227))
       FR(228) = 2*M**2*P(2)*(-(D1234R(1,1)*P(5)*P(6))-D01234R*P(36)+(
     -   D1234R(2,1)+2*(D1234R(3,2)+D1234R(4,2)-D1234R(5,2)-D1234R(6,
     -   2)))*P(8)*P(37)+D1234R(3,1)*P(41))
       FI(228) = 2*M**2*P(2)*(-(D1234I(1,1)*P(5)*P(6))-D01234I*P(36)+(
     -   D1234I(2,1)+2*(D1234I(3,2)+D1234I(4,2)-D1234I(5,2)-D1234I(6,
     -   2)))*P(8)*P(37)+D1234I(3,1)*P(41))
       F(228)=DCMPLX(FR(228),FI(228))
       FR(229) = -4*((C234R(2,1)+C234R(2,2))*P(11)+M**2*P(2)*(D1234R(3
     -   ,1)*P(36)+(-D1234R(3,2)+D1234R(6,2))*P(8)*P(37)))
       FI(229) = -4*((C234I(2,1)+C234I(2,2))*P(11)+M**2*P(2)*(D1234I(3
     -   ,1)*P(36)+(-D1234I(3,2)+D1234I(6,2))*P(8)*P(37)))
       F(229)=DCMPLX(FR(229),FI(229))
       FR(230) = 2*(-2*(C234R(2,1)+C234R(2,2))*P(11)+M**2*P(2)*(D1234R
     -   (1,1)*P(25)*P(35)+D01234R*P(36)+(D1234R(2,1)+2*(D1234R(3,2)+
     -   D1234R(4,2)-D1234R(5,2)-D1234R(6,2)))*P(8)*P(37)-D1234R(3,1)
     -   *P(45)))
       FI(230) = 2*(-2*(C234I(2,1)+C234I(2,2))*P(11)+M**2*P(2)*(D1234I
     -   (1,1)*P(25)*P(35)+D01234I*P(36)+(D1234I(2,1)+2*(D1234I(3,2)+
     -   D1234I(4,2)-D1234I(5,2)-D1234I(6,2)))*P(8)*P(37)-D1234I(3,1)
     -   *P(45)))
       F(230)=DCMPLX(FR(230),FI(230))
       FR(231) = 4*((-C234R(2,2)+C234R(3,2))*P(11)+M**2*P(2)*((D1234R(
     -   2,2)+D1234R(3,2)-2*D1234R(6,2))*P(8)*P(37)+(D1234R(2,1)-D123
     -   4R(3,1))*P(46)))
       FI(231) = 4*((-C234I(2,2)+C234I(3,2))*P(11)+M**2*P(2)*((D1234I(
     -   2,2)+D1234I(3,2)-2*D1234I(6,2))*P(8)*P(37)+(D1234I(2,1)-D123
     -   4I(3,1))*P(46)))
       F(231)=DCMPLX(FR(231),FI(231))
       FR(232) = 2*(-2*(C234R(2,1)+C234R(3,2))*P(11)+M**2*P(2)*(-(D123
     -   4R(3,1)*P(3)*P(4))+D1234R(1,1)*P(5)*P(6)+D01234R*P(36)+D1234
     -   R(2,1)*P(8)*P(37)))
       FI(232) = 2*(-2*(C234I(2,1)+C234I(3,2))*P(11)+M**2*P(2)*(-(D123
     -   4I(3,1)*P(3)*P(4))+D1234I(1,1)*P(5)*P(6)+D01234I*P(36)+D1234
     -   I(2,1)*P(8)*P(37)))
       F(232)=DCMPLX(FR(232),FI(232))
       FR(233) = 4*(C123R(1,1)+C123R(1,2)-C123R(2,1)-C123R(3,2)-C234R(
     -   2,1)-C234R(3,2))*P(11)
       FI(233) = 4*(C123I(1,1)+C123I(1,2)-C123I(2,1)-C123I(3,2)-C234I(
     -   2,1)-C234I(3,2))*P(11)
       F(233)=DCMPLX(FR(233),FI(233))
       FR(234) = 4*(-C123R(2,2)+C123R(3,2)+C234R(1,1)+C234R(1,2)-C234R
     -   (2,1)-C234R(3,2))*P(11)
       FI(234) = 4*(-C123I(2,2)+C123I(3,2)+C234I(1,1)+C234I(1,2)-C234I
     -   (2,1)-C234I(3,2))*P(11)
       F(234)=DCMPLX(FR(234),FI(234))
       FR(235) = -4*M**2*P(2)*((D1234R(1,2)+D1234R(3,2)-2*D1234R(5,2))
     -   *P(3)*P(4)+D01234R*P(36)+(-D1234R(1,1)+D1234R(3,1))*P(47))
       FI(235) = -4*M**2*P(2)*((D1234I(1,2)+D1234I(3,2)-2*D1234I(5,2))
     -   *P(3)*P(4)+D01234I*P(36)+(-D1234I(1,1)+D1234I(3,1))*P(47))
       F(235)=DCMPLX(FR(235),FI(235))
       FR(236) = 4*M**2*P(2)*((-D1234R(3,2)+D1234R(5,2))*P(3)*P(4)+D12
     -   34R(3,1)*P(36))
       FI(236) = 4*M**2*P(2)*((-D1234I(3,2)+D1234I(5,2))*P(3)*P(4)+D12
     -   34I(3,1)*P(36))
       F(236)=DCMPLX(FR(236),FI(236))
       FR(237) = -2*M**2*P(2)*(2*(D1234R(3,2)+D1234R(4,2)-D1234R(5,2)-
     -   D1234R(6,2))*P(3)*P(4)+D1234R(1,1)*P(5)*P(6)+D01234R*P(36)+D
     -   1234R(2,1)*P(50)+D1234R(3,1)*P(53))
       FI(237) = -2*M**2*P(2)*(2*(D1234I(3,2)+D1234I(4,2)-D1234I(5,2)-
     -   D1234I(6,2))*P(3)*P(4)+D1234I(1,1)*P(5)*P(6)+D01234I*P(36)+D
     -   1234I(2,1)*P(50)+D1234I(3,1)*P(53))
       F(237)=DCMPLX(FR(237),FI(237))
       FR(238) = 4*((-C234R(1,1)-C234R(1,2)+C234R(2,1)+C234R(3,2))*P(1
     -   1)+M**2*P(2)*(-((D1234R(2,2)+D1234R(3,2)-2*D1234R(6,2))*P(3)
     -   *P(4))+(-D1234R(2,1)+D1234R(3,1))*P(36)))
       FI(238) = 4*((-C234I(1,1)-C234I(1,2)+C234I(2,1)+C234I(3,2))*P(1
     -   1)+M**2*P(2)*(-((D1234I(2,2)+D1234I(3,2)-2*D1234I(6,2))*P(3)
     -   *P(4))+(-D1234I(2,1)+D1234I(3,1))*P(36)))
       F(238)=DCMPLX(FR(238),FI(238))
       FR(239) = 4*((C234R(2,1)+C234R(3,2))*P(11)+M**2*P(2)*((-D1234R(
     -   3,2)+D1234R(6,2))*P(3)*P(4)+D1234R(3,1)*P(36)))
       FI(239) = 4*((C234I(2,1)+C234I(3,2))*P(11)+M**2*P(2)*((-D1234I(
     -   3,2)+D1234I(6,2))*P(3)*P(4)+D1234I(3,1)*P(36)))
       F(239)=DCMPLX(FR(239),FI(239))
       FR(240) = -2*(-2*(C234R(2,1)+C234R(3,2))*P(11)+M**2*P(2)*(2*(D1
     -   234R(3,2)+D1234R(4,2)-D1234R(5,2)-D1234R(6,2))*P(3)*P(4)+D12
     -   34R(1,1)*P(25)*P(35)+D01234R*P(36)+D1234R(3,1)*P(53)+D1234R(
     -   2,1)*P(54)))
       FI(240) = -2*(-2*(C234I(2,1)+C234I(3,2))*P(11)+M**2*P(2)*(2*(D1
     -   234I(3,2)+D1234I(4,2)-D1234I(5,2)-D1234I(6,2))*P(3)*P(4)+D12
     -   34I(1,1)*P(25)*P(35)+D01234I*P(36)+D1234I(3,1)*P(53)+D1234I(
     -   2,1)*P(54)))
       F(240)=DCMPLX(FR(240),FI(240))
       FR(241) = -2*M**2*P(2)*(-(D1234R(3,1)*P(3)*P(4))+D1234R(1,1)*P(
     -   5)*P(6)+D01234R*P(36)+D1234R(2,1)*P(8)*P(37))
       FI(241) = -2*M**2*P(2)*(-(D1234I(3,1)*P(3)*P(4))+D1234I(1,1)*P(
     -   5)*P(6)+D01234I*P(36)+D1234I(2,1)*P(8)*P(37))
       F(241)=DCMPLX(FR(241),FI(241))
       FR(242) = 2*M**2*P(2)*(-(D1234R(3,1)*P(3)*P(4))+D1234R(1,1)*P(5
     -   )*P(6)+D01234R*P(36)+D1234R(2,1)*P(8)*P(37))
       FI(242) = 2*M**2*P(2)*(-(D1234I(3,1)*P(3)*P(4))+D1234I(1,1)*P(5
     -   )*P(6)+D01234I*P(36)+D1234I(2,1)*P(8)*P(37))
       F(242)=DCMPLX(FR(242),FI(242))
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
     
c               PART THAT DEPENDS ON THE EXTERNAL CURRENT
c************************************************************************************
c************************************************************************************
c************************************************************************************
       p1mup1 = dotrc(p1,mup1)
       p1mup2 = dotrc(p1,mup2)
       p1mup3 = dotrc(p1,mup3)
       p1mup4 = dotrc(p1,mup4)
       p2mup1 = dotrc(p2,mup1)
       p2mup2 = dotrc(p2,mup2)
       p2mup3 = dotrc(p2,mup3)
       p2mup4 = dotrc(p2,mup4)
       p3mup1 = dotrc(p3,mup1)
       p3mup2 = dotrc(p3,mup2)
       p3mup3 = dotrc(p3,mup3)
       p3mup4 = dotrc(p3,mup4)
       p4mup1 = dotrc(p4,mup1)
       p4mup2 = dotrc(p4,mup2)
       p4mup3 = dotrc(p4,mup3)
       p4mup4 = dotrc(p4,mup4)
       mup1mup2 = dotcc(mup1,mup2)
       mup1mup3 = dotcc(mup1,mup3)
       mup1mup4 = dotcc(mup1,mup4)
       mup2mup3 = dotcc(mup2,mup3)
       mup2mup4 = dotcc(mup2,mup4)
       mup3mup4 = dotcc(mup3,mup4)
       call boxOddFa(p1,p2,p3,p4,mup1,mup2,mup3,mup4,Fa(1))
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


       resultn(1) = (0,-1)*(F(1)*Fa(1)+(p4mup4*F(2)+p1mup4*F(3)+p2mup4
     -   *F(4))*Fa(2)+(p4mup4*F(5)+p1mup4*F(6)+p2mup4*F(7))*Fa(3)+(p4
     -   mup4*F(8)+p1mup4*F(9)+p2mup4*F(10))*Fa(4)+(p1mup3*F(11)+p4mu
     -   p3*F(12)+p2mup3*F(13))*Fa(5)+(p2mup3*F(14)+p4mup3*F(15)+p1mu
     -   p3*F(16))*Fa(6)+(p4mup3*F(17)+p1mup3*F(18)+p2mup3*F(19))*Fa(
     -   7)+(p4mup2*F(45)+p1mup2*F(46)+p2mup2*F(47))*Fa(8)+(p4mup2*F(
     -   48)+p1mup2*F(49)+p2mup2*F(50))*Fa(9)+(p4mup2*F(51)+p1mup2*F(
     -   52)+p2mup2*F(53))*Fa(10)+(p4mup1*F(29)+p2mup1*F(83))*Fa(11)+
     -   (p1mup1*F(112)+p4mup1*F(113))*Fa(12)+(p2mup1*F(114)+p1mup1*F
     -   (115))*Fa(13)+(mup3mup4*F(20)+p1mup3*p1mup4*F(21)+p1mup4*p4m
     -   up3*F(22)+p1mup3*p4mup4*F(23)+p4mup3*p4mup4*F(24)+p2mup3*p2m
     -   up4*F(25)+p1mup4*p2mup3*F(26)+p1mup3*p2mup4*F(26)+p2mup4*p4m
     -   up3*F(27)+p2mup3*p4mup4*F(28))*Fa(14)+(mup3mup4*F(29)+p2mup3
     -   *p4mup4*F(30)+p1mup4*p2mup3*F(31)+p2mup4*p4mup3*F(32)+p1mup3
     -   *p2mup4*F(33)+p2mup3*p2mup4*F(34))*Fa(15)+(mup3mup4*F(35)+p4
     -   mup3*p4mup4*F(36)+p1mup3*p4mup4*F(37)+p1mup4*p4mup3*F(38)+p1
     -   mup3*p1mup4*F(39)+p2mup3*p4mup4*F(40)+p2mup4*p4mup3*F(41)+p1
     -   mup3*p2mup4*F(42)+p1mup4*p2mup3*F(43)+p2mup3*p2mup4*F(44))*F
     -   a(16)+(mup2mup4*F(54)+p4mup2*p4mup4*F(55)+p1mup4*p4mup2*F(56
     -   )+p1mup2*p4mup4*F(57)+p1mup2*p1mup4*F(58)+p1mup4*p2mup2*F(59
     -   )+p1mup2*p2mup4*F(60)+p2mup4*p4mup2*F(61)+p2mup2*p4mup4*F(62
     -   )+p2mup2*p2mup4*F(63))*Fa(17)+(mup2mup4*F(64)+p4mup2*p4mup4*
     -   F(65)+p1mup2*p1mup4*F(66)+p1mup2*p4mup4*F(67)+p1mup4*p4mup2*
     -   F(68)+p2mup2*p2mup4*F(69)+p2mup2*p4mup4*F(70)+p2mup4*p4mup2*
     -   F(71)+p1mup4*p2mup2*F(72)+p1mup2*p2mup4*F(73))*Fa(18)+(mup2m
     -   up4*F(74)+p4mup2*p4mup4*F(75)+p1mup2*p1mup4*F(76)+p1mup2*p4m
     -   up4*F(77)+p1mup4*p4mup2*F(78)+p2mup2*p2mup4*F(79)+p2mup4*p4m
     -   up2*F(80)+p2mup2*p4mup4*F(80)+p1mup2*p2mup4*F(81)+p1mup4*p2m
     -   up2*F(82))*Fa(19)+(mup2mup3*F(83)+p1mup2*p4mup3*F(84)+p2mup2
     -   *p4mup3*F(85)+p1mup3*p4mup2*F(86)+p4mup2*p4mup3*F(87)+p2mup3
     -   *p4mup2*F(88))*Fa(20)+(mup2mup3*F(89)+p1mup2*p1mup3*F(90)+p1
     -   mup3*p4mup2*F(91)+p1mup2*p4mup3*F(91)+p4mup2*p4mup3*F(92)+p2
     -   mup2*p2mup3*F(93)+p1mup2*p2mup3*F(94)+p1mup3*p2mup2*F(95)+p2
     -   mup3*p4mup2*F(96)+p2mup2*p4mup3*F(97))*Fa(21)+(p1mup3*p4mup2
     -   *F(96)+mup2mup3*F(98)+p2mup2*p2mup3*F(99)+p2mup2*p4mup3*F(10
     -   0)+p2mup3*p4mup2*F(101)+p1mup2*p2mup3*F(102)+p1mup3*p2mup2*F
     -   (103)+p4mup2*p4mup3*F(104)+p1mup2*p1mup3*F(105)+p1mup2*p4mup
     -   3*F(106))*Fa(22)+(p4mup1*p4mup4*F(30)+p1mup4*p4mup1*F(31)+mu
     -   p1mup4*F(83)+p2mup4*p4mup1*F(107))*Fa(23)+(mup1mup4*F(29)+p2
     -   mup1*p4mup4*F(116)+p1mup4*p2mup1*F(117)+p2mup1*p2mup4*F(118)
     -   )*Fa(24)+(p1mup1*p4mup4*F(30)+p1mup1*p1mup4*F(31)+p1mup1*p2m
     -   up4*F(107)+mup1mup4*F(113))*Fa(25)+(p4mup1*p4mup3*F(32)+p1mu
     -   p3*p4mup1*F(33)+p2mup3*p4mup1*F(111)+mup1mup3*F(112))*Fa(26)
     -   +(p1mup3*p2mup1*F(86)+p2mup1*p2mup3*F(88)+mup1mup3*F(115)+p2
     -   mup1*p4mup3*F(119))*Fa(27)+(p1mup1*p4mup3*F(32)+p1mup1*p1mup
     -   3*F(33)+p1mup1*p2mup3*F(111)+mup1mup3*F(114))*Fa(28)+(mup1mu
     -   p2*F(83)+p4mup1*p4mup2*F(108)+p1mup2*p4mup1*F(109)+p2mup2*p4
     -   mup1*F(110))*Fa(29)+(mup1mup2*F(29)+p1mup2*p2mup1*F(84)+p2mu
     -   p1*p2mup2*F(85)+p2mup1*p4mup2*F(120))*Fa(30)+(p1mup1*p4mup2*
     -   F(108)+p1mup1*p1mup2*F(109)+p1mup1*p2mup2*F(110)+mup1mup2*F(
     -   113))*Fa(31)+(mup2mup3*p4mup4*F(30)+mup2mup3*p1mup4*F(31)+mu
     -   p2mup4*p4mup3*F(32)+mup2mup4*p1mup3*F(33)+mup2mup3*p2mup4*F(
     -   107)+mup3mup4*p4mup2*F(108)+mup3mup4*p1mup2*F(109)+mup3mup4*
     -   p2mup2*F(110)+mup2mup4*p2mup3*F(111))*Fa(32)+(mup1mup4*p1mup
     -   3*F(86)+mup1mup4*p2mup3*F(88)+mup1mup3*p4mup4*F(116)+mup1mup
     -   3*p1mup4*F(117)+mup1mup3*p2mup4*F(118)+mup1mup4*p4mup3*F(119
     -   ))*Fa(33)+(mup1mup2*p4mup4*F(30)+mup1mup2*p1mup4*F(31)+mup1m
     -   up4*p1mup2*F(84)+mup1mup4*p2mup2*F(85)+mup1mup2*p2mup4*F(107
     -   )+mup1mup4*p4mup2*F(120))*Fa(34)+(mup1mup2*p4mup3*F(32)+mup1
     -   mup2*p1mup3*F(33)+mup1mup3*p4mup2*F(108)+mup1mup3*p1mup2*F(1
     -   09)+mup1mup3*p2mup2*F(110)+mup1mup2*p2mup3*F(111))*Fa(35))


        if (comp.gt.0) then    


       resultn(2) = (0,-1)*(F(121)*Fa(11)+F(122)*Fa(12)+F(123)*Fa(13)+
     -   (p2mup4*F(124)+p4mup4*F(125)+p1mup4*F(126))*Fa(23)+(p4mup4*F
     -   (127)+p1mup4*F(128)+p2mup4*F(129))*Fa(24)+(p4mup4*F(130)+p1m
     -   up4*F(131)+p2mup4*F(132))*Fa(25)+(p2mup3*F(133)+p4mup3*F(134
     -   )+p1mup3*F(135))*Fa(26)+(p4mup3*F(136)+p1mup3*F(137)+p2mup3*
     -   F(138))*Fa(27)+(p4mup3*F(139)+p1mup3*F(140)+p2mup3*F(141))*F
     -   a(28)+(p4mup2*F(148)+p1mup2*F(149)+p2mup2*F(150))*Fa(29)+(p4
     -   mup2*F(151)+p1mup2*F(152)+p2mup2*F(153))*Fa(30)+(p4mup2*F(15
     -   4)+p1mup2*F(155)+p2mup2*F(156))*Fa(31)+(p1mup3*p4mup4*F(92)+
     -   p2mup3*p2mup4*F(100)+p2mup3*p4mup4*F(104)+mup3mup4*F(142)+p4
     -   mup3*p4mup4*F(143)+p1mup3*p1mup4*F(144)+p1mup4*p4mup3*F(145)
     -   +p2mup4*p4mup3*F(146)+p1mup4*p2mup3*F(147)+p1mup3*p2mup4*F(1
     -   47))*Fa(33)+(mup2mup4*F(157)+p4mup2*p4mup4*F(158)+p1mup4*p4m
     -   up2*F(159)+p1mup2*p4mup4*F(160)+p1mup2*p1mup4*F(161)+p2mup4*
     -   p4mup2*F(162)+p2mup2*p4mup4*F(162)+p1mup4*p2mup2*F(163)+p1mu
     -   p2*p2mup4*F(164)+p2mup2*p2mup4*F(165))*Fa(34)+(p4mup2*p4mup3
     -   *F(40)+p2mup2*p4mup3*F(44)+mup2mup3*F(166)+p1mup2*p4mup3*F(1
     -   67)+p1mup3*p4mup2*F(168)+p1mup2*p1mup3*F(169)+p2mup3*p4mup2*
     -   F(170)+p1mup2*p2mup3*F(171)+p1mup3*p2mup2*F(172)+p2mup2*p2mu
     -   p3*F(173))*Fa(35))
       resultn(3) = (0,-1)*(F(174)*Fa(8)+F(175)*Fa(9)+F(176)*Fa(10)+(p
     -   4mup4*F(177)+p1mup4*F(178)+p2mup4*F(179))*Fa(17)+(p2mup4*F(1
     -   80)+p4mup4*F(181)+p1mup4*F(182))*Fa(18)+(p4mup4*F(183)+p1mup
     -   4*F(184)+p2mup4*F(185))*Fa(19)+(p1mup3*F(186)+p4mup3*F(187)+
     -   p2mup3*F(188))*Fa(20)+(p2mup3*F(189)+p1mup3*F(190)+p4mup3*F(
     -   191))*Fa(21)+(p4mup3*F(192)+p1mup3*F(193)+p2mup3*F(194))*Fa(
     -   22)+p4mup1*F(195)*Fa(29)+p2mup1*F(196)*Fa(30)+p1mup1*F(195)*
     -   Fa(31)+mup3mup4*F(195)*Fa(32)+mup1mup4*F(196)*Fa(34)+mup1mup
     -   3*F(195)*Fa(35))
       resultn(4) = (0,-1)*(F(197)*Fa(5)+F(198)*Fa(6)+F(199)*Fa(7)+(p4
     -   mup4*F(200)+p1mup4*F(201)+p2mup4*F(202))*Fa(14)+(p4mup4*F(20
     -   3)+p1mup4*F(204)+p2mup4*F(205))*Fa(15)+(p4mup4*F(206)+p1mup4
     -   *F(207)+p2mup4*F(208))*Fa(16)+(p4mup2*F(209)+p1mup2*F(210)+p
     -   2mup2*F(211))*Fa(20)+(p2mup2*F(212)+p1mup2*F(213)+p4mup2*F(2
     -   14))*Fa(21)+(p4mup2*F(215)+p1mup2*F(216)+p2mup2*F(217))*Fa(2
     -   2)+p4mup1*F(218)*Fa(26)+p2mup1*F(219)*Fa(27)+p1mup1*F(218)*F
     -   a(28)+mup2mup4*F(218)*Fa(32)+mup1mup4*F(219)*Fa(33)+mup1mup2
     -   *F(218)*Fa(35))
       resultn(5) = (0,-1)*(F(220)*Fa(2)+F(221)*Fa(3)+F(222)*Fa(4)+(p4
     -   mup3*F(223)+p1mup3*F(224)+p2mup3*F(225))*Fa(14)+(p4mup3*F(22
     -   6)+p1mup3*F(227)+p2mup3*F(228))*Fa(15)+(p4mup3*F(229)+p1mup3
     -   *F(230)+p2mup3*F(231))*Fa(16)+(p4mup2*F(232)+p1mup2*F(233)+p
     -   2mup2*F(234))*Fa(17)+(p1mup2*F(235)+p4mup2*F(236)+p2mup2*F(2
     -   37))*Fa(18)+(p2mup2*F(238)+p4mup2*F(239)+p1mup2*F(240))*Fa(1
     -   9)+p4mup1*F(241)*Fa(23)+p2mup1*F(242)*Fa(24)+p1mup1*F(241)*F
     -   a(25)+mup2mup3*F(241)*Fa(32)+mup1mup3*F(242)*Fa(33)+mup1mup2
     -   *F(241)*Fa(34))
       resultn(2) =-resultn(2)
       resultn(3) =-resultn(3)
       resultn(4) =-resultn(4)
       resultn(5) =-resultn(5)


        endif 
c************************************************************************************


c The Finite virtual piece should be multiplied for (-1)  since 
c I have multiplied by (-I) to get the F's and k's without (I) factor
c . The factorization from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps Gamma[1+Eps] (musq)^(-Eps))
c  c So, I*I=(-1)!!!
       result=-resultn(1)
c************************************************************************************
c************************************************************************************
       Return
       End
