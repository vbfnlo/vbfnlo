       subroutine PrintIntPen2gMassLess_QUAD(p1,p2,p3,p4,p5,musq)
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 27/10/2008
c Modified:17/11/2008
c ************************************************************************************
c determine the Hexagon of
c p1(gluon)    ---->------->--------->--------->---   p6(B) p6^2=M_h
c                $         $          $         $     
c                $         $          $         $     
c                $         $          $         $     
c                V         V          V         V     
c                $         $          $         $     
c                $         $          $         $     
c               p2         p3         p4        p5     
c              mu_p2      mu_p3      mu_p4     mu_p5      
c Note: No couplings are assigned
c Note:m the mass of the fermion in the loop  
c Note: To make it shorter in the promgram: mu_p2,...->mup2,... 
c Notation of External momenta: p1+p2+p3+p4+p5+p6=0 
c mu_p1,mu_p2,mu_p3,mu_p4,mu_p5 should be think as external current 
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
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
      Real*8 psq(6)
        Real*8   p1sq
       Real*8   p2sq 
       Real*8   p3sq 
       Real*8   p4sq 
       Real*8   p5sq
       Real*8   s12, s15
       Real*8   s23
       Real*8   s34 
       Real*8   s45
       Complex*16  B012,B013,B014,B015
       Complex*16  B023,B024,B025
       Complex*16  B034,B035
       Complex*16  B045
       Real*8  B012R,B013R,B014R,B015R
       Real*8  B023R,B024R,B025R
       Real*8  B034R,B035R
       Real*8  B045R
       Real*8  B012I,B013I,B014I,B015I
       Real*8  B023I,B024I,B025I
       Real*8  B034I,B035I
       Real*8  B045I
       Complex*16 C0123,C0124,C0125
       Complex*16 C0134,C0135
       Complex*16 C0145
       Complex*16 C0234,C0235
       Complex*16 C0245
       Complex*16 C0345
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
       Real*8 C123R(4,2),C124R(4,2),C125R(4,2)
       Real*8 C134R(4,2),C135R(4,2)
       Real*8 C145R(4,2) 
       Real*8 C234R(4,2),C235R(4,2)
       Real*8 C245R(4,2)
       Real*8 C345R(4,2)
       Real*8 C123I(4,2),C124I(4,2),C125I(4,2)
       Real*8 C134I(4,2),C135I(4,2)
       Real*8 C145I(4,2)
       Real*8 C234I(4,2),C235I(4,2)
       Real*8 C245I(4,2)
       Real*8 C345I(4,2)
       Complex*16  D01234,D01235
       Complex*16 D01245
       Complex*16 D01345
       Complex*16 D02345
       Real*8 D01234R,D01235R
       Real*8 D01245R
       Real*8 D01345R
       Real*8 D02345R
       Real*8 D01234I,D01235I
       Real*8 D01245I
       Real*8 D01345I
       Real*8 D02345I
       Real*8 D1234R(13,3),D1235R(13,3)
       Real*8 D1245R(13,3)
       Real*8 D1345R(13,3)
       Real*8 D2345R(13,3)
       Real*8 D1234I(13,3),D1235I(13,3)
       Real*8 D1245I(13,3)
       Real*8 D1345I(13,3)
       Real*8 D2345I(13,3)
       Complex*16 E012345
       Real*8 E012345R
       Real*8 E012345I
       Real*8 E12345R(46,4)
       Real*8 E12345I(46,4)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  QUAD precision start here
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       Real*16  p1_QUAD(0:3),p2_QUAD(0:3),p3_QUAD(0:3),p4_QUAD(0:3),
     -          p5_QUAD(0:3)
       Real*16 ptemp(0:3,5)
       Real*16   p1sq_QUAD, p1p2_QUAD, p1p3_QUAD, p1p4_QUAD, p1p5_QUAD
       Real*16   p2sq_QUAD, p2p3_QUAD, p2p4_QUAD, p2p5_QUAD 
       Real*16   p3sq_QUAD, p3p4_QUAD, p3p5_QUAD 
       Real*16   p4sq_QUAD, p4p5_QUAD 
       Real*16   p5sq_QUAD 
       Real*16   s12_QUAD, s13_QUAD, s14_QUAD, s15_QUAD
       Real*16   s23_QUAD
       Real*16   s34_QUAD
       Real*16   s45_QUAD

       Complex*32  B012_QUAD,B013_QUAD,B014_QUAD,B015_QUAD
       Complex*32  B023_QUAD,B024_QUAD,B025_QUAD
       Complex*32  B034_QUAD,B035_QUAD
       Complex*32  B045_QUAD
       Real*16  B012R_QUAD,B013R_QUAD,B014R_QUAD,B015R_QUAD
       Real*16  B023R_QUAD,B024R_QUAD,B025R_QUAD
       Real*16  B034R_QUAD,B035R_QUAD      
       Real*16  B045R_QUAD                 
       Real*16  B012I_QUAD,B013I_QUAD,B014I_QUAD,B015I_QUAD
       Real*16  B023I_QUAD,B024I_QUAD,B025I_QUAD
       Real*16  B034I_QUAD,B035I_QUAD
       Real*16  B045I_QUAD
       Complex*32    C0123_QUAD,C0124_QUAD,C0125_QUAD
       Complex*32    C0134_QUAD,C0135_QUAD
       Complex*32    C0145_QUAD      
       Complex*32    C0234_QUAD,C0235_QUAD
       Complex*32    C0245_QUAD
       Complex*32    C0345_QUAD
       Real*16       C0123R_QUAD,C0124R_QUAD,C0125R_QUAD
       Real*16       C0134R_QUAD,C0135R_QUAD
       Real*16       C0145R_QUAD       
       Real*16       C0234R_QUAD,C0235R_QUAD
       Real*16       C0245R_QUAD       
       Real*16       C0345R_QUAD       
       Real*16       C0123I_QUAD,C0124I_QUAD,C0125I_QUAD
       Real*16       C0134I_QUAD,C0135I_QUAD
       Real*16       C0145I_QUAD       
       Real*16       C0234I_QUAD,C0235I_QUAD
       Real*16       C0245I_QUAD
       Real*16       C0345I_QUAD
       Real*16       C123R_QUAD(4,2),C124R_QUAD(4,2),C125R_QUAD(4,2)
       Real*16       C134R_QUAD(4,2),C135R_QUAD(4,2)
       Real*16       C145R_QUAD(4,2)      
       Real*16       C234R_QUAD(4,2),C235R_QUAD(4,2)
       Real*16       C245R_QUAD(4,2)      
       Real*16       C345R_QUAD(4,2)      
       Real*16       C123I_QUAD(4,2),C124I_QUAD(4,2),C125I_QUAD(4,2)
       Real*16       C134I_QUAD(4,2),C135I_QUAD(4,2)
       Real*16       C145I_QUAD(4,2)      
       Real*16       C234I_QUAD(4,2),C235I_QUAD(4,2)
       Real*16       C245I_QUAD(4,2)
       Real*16       C345I_QUAD(4,2)
       Complex*32  D01234_QUAD,D01235_QUAD
       Complex*32  D01245_QUAD
       Complex*32  D01345_QUAD
       Complex*32  D02345_QUAD
       Real*16     D01234R_QUAD,D01235R_QUAD
       Real*16     D01245R_QUAD
       Real*16     D01345R_QUAD
       Real*16     D02345R_QUAD
       Real*16     D01234I_QUAD,D01235I_QUAD
       Real*16     D01245I_QUAD
       Real*16     D01345I_QUAD
       Real*16     D02345I_QUAD
       Real*16     D1234R_QUAD(13,3),D1235R_QUAD(13,3)
       Real*16     D1245R_QUAD(13,3)
       Real*16     D1345R_QUAD(13,3)
       Real*16     D2345R_QUAD(13,3)
       Real*16     D1234I_QUAD(13,3),D1235I_QUAD(13,3)
       Real*16     D1245I_QUAD(13,3)
       Real*16     D1345I_QUAD(13,3)
       Real*16     D2345I_QUAD(13,3)
       Complex*32 E012345_QUAD
       Real*16    E012345R_QUAD
       Real*16    E012345I_QUAD
       Real*16    E12345R_QUAD(46,4)
       Real*16    E12345I_QUAD(46,4)

       Real*16 musq_QUAD, dotrr_QUAD
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      End QUAD precision
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       COMMON/NoAbe2gaPenLineEVALFunctions/ E012345R,E12345R,E012345I,E12345I
       COMMON/NoAbe2gaPenLineDVALFunctions/ D01234R,D1234R,D01234I,D1234I, 
     &    D01235R,D1235R,D01235I,D1235I,
     &    D01245R,D1245R,D01245I,D1245I, 
     &    D01345R,D1345R,D01345I,D1345I, 
     &    D02345R,D2345R,D02345I,D2345I 

       COMMON/NoAbe2gaPenLineCVALFunctions/ C0123R,C123R,C0123I,C123I, 
     &    C0124R,C124R,C0124I,C124I, C0125R,C125R,C0125I,C125I,
     &    C0134R,C134R,C0134I,C134I, C0135R,C135R,C0135I,C135I,
     &    C0145R,C145R,C0145I,C145I, C0234R,C234R,C0234I,C234I, 
     &    C0235R,C235R,C0235I,C235I, C0245R,C245R,C0245I,C245I, 
     &    C0345R,C345R,C0345I,C345I
       COMMON/NoAbe2gaPenLineBVALFunctions/ B012R,B012I, B013R,B013I, B014R,B014I, 
     &    B015R,B015I,  B023R,B023I,  B024R,B024I, B025R,B025I,
     &    B034R,B034I,  B035R,B035I,   B045R,B045I
      real*8 is12,is45,is12s45
       COMMON/InvariantsP/p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,s45,
     &     s15, Is12,Is45,Is12s45
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       Complex*16 B0fin,C0fin,D0fin,E0fin
c       EXTERNAL B0fin,C0fin,D0fin,E0fin
       Logical PrintB,PrintC,PrintD,PrintE

       External dotrr_QUAD ,dotrr
      Complex*16 dotrc,dotcc
       Real*8 dotrr
       Real*8 musq
      EXTERNAL   dotrc,dotcc
       Integer i,j,jtemp

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
      PrintB=.False. 
      PrintC=.False. 
      PrintD=.False. 
      PrintE=.False. 

      do i=0,3 
       ptemp(i,1)=p1(i)
       ptemp(i,2)=p2(i)
       ptemp(i,3)=p3(i)
       ptemp(i,4)=p4(i)
       ptemp(i,5)=p5(i)
 
       enddo
       psq(1)=abs(dotrr(p1,p1))
       psq(2)=abs(dotrr(p2,p2))
       psq(3)=abs(dotrr(p3,p3))
       psq(4)=abs(dotrr(p4,p4))
       psq(5)=abs(dotrr(p5,p5))

       if (psq(1).gt.10d-7) then
       j=1
       elseif (psq(2).gt.10d-7)then
       j=2
       elseif (psq(3).gt.10d-7)then
       j=3
       elseif (psq(4).gt.10d-7)then
       j=4
       elseif (psq(5).gt.10d-7)then
       j=5
       else
       j=7
c$$$       write(*,*) 'No massive invariant'
c$$$       Return
       endif
       if(j.ne.7) then
       do jtemp=1,5
       if (psq(jtemp).le.10d-7) then
       if(abs(ptemp(1,jtemp)).le.10d-7)then
       if(abs(ptemp(2,jtemp)).le.10d-7) then
       ptemp(0,jtemp)=abs(ptemp(3,jtemp))
       else
       if (ptemp(0,jtemp).ge.0d0) then
       ptemp(0,jtemp)=SQRT(ptemp(1,jtemp)*ptemp(1,jtemp)+
     -                  ptemp(2,jtemp)*ptemp(2,jtemp)+
     -                  ptemp(3,jtemp)*ptemp(3,jtemp))
       else
c$$$       print*, 'Minus jtemp', jtemp
      ptemp(0,jtemp)=-SQRT(ptemp(1,jtemp)*ptemp(1,jtemp)+
     -                  ptemp(2,jtemp)*ptemp(2,jtemp)+
     -                  ptemp(3,jtemp)*ptemp(3,jtemp))
      endif
      endif
      else
      if (ptemp(0,jtemp).ge.0d0) then
c$$$       print*, 'jtemp', jtemp
      ptemp(0,jtemp)=SQRT(ptemp(1,jtemp)*ptemp(1,jtemp)+
     -                  ptemp(2,jtemp)*ptemp(2,jtemp)+
     -                  ptemp(3,jtemp)*ptemp(3,jtemp))
      else
c$$$       print*, 'Minus jtemp', jtemp
      ptemp(0,jtemp)=-SQRT(ptemp(1,jtemp)*ptemp(1,jtemp)+
     -                  ptemp(2,jtemp)*ptemp(2,jtemp)+
     -                  ptemp(3,jtemp)*ptemp(3,jtemp))
      endif
      endif
      endif
      enddo
      do i=0,3
      ptemp(i,j)=0.0
      do jtemp=1,5
      if (jtemp.ne.j) then
      ptemp(i,j)=-ptemp(i,jtemp)+ptemp(i,j)
      endif
      enddo
      enddo
      elseif(j.eq.7) then
      do i=0,3 
      ptemp(i,1)=p1(i)
      ptemp(i,2)=p2(i)
      ptemp(i,3)=p3(i)
      ptemp(i,4)=p4(i)
      ptemp(i,5)=p5(i)

      enddo 
      endif


      do i=0,3 
      p1_QUAD(i)=ptemp(i,1)
      p2_QUAD(i)=ptemp(i,2)
      p3_QUAD(i)=ptemp(i,3)
      p4_QUAD(i)=ptemp(i,4)
      p5_QUAD(i)=ptemp(i,5)
       enddo

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
c       s24_QUAD = (p2sq_QUAD +p4sq_QUAD+ 2*p2p4_QUAD) 
c       s25_QUAD = (p2sq_QUAD +p5sq_QUAD+ 2*p2p5_QUAD) 
       s34_QUAD = (p3sq_QUAD +p4sq_QUAD+ 2*p3p4_QUAD) 
c       s35_QUAD = (p3sq_QUAD +p5sq_QUAD+ 2*p3p5_QUAD) 
       s45_QUAD = (p4sq_QUAD +p5sq_QUAD+ 2*p4p5_QUAD) 

cfc       Write(*,'(a5,E44.36)')," p1sq_QUAD ", p1sq_QUAD 
cfc       Write(*,'(a5,E44.36)')," p1p2_QUAD ", p1p2_QUAD
cfc       Write(*,'(a5,E44.36)')," p1p3_QUAD ", p1p3_QUAD
cfc       Write(*,'(a5,E44.36)')," p1p4_QUAD ", p1p4_QUAD
cfc       Write(*,'(a5,E44.36)')," p1p5_QUAD ", p1p5_QUAD
cfc       Write(*,'(a5,E44.36)')," p2sq_QUAD ", p2sq_QUAD 
cfc       Write(*,'(a5,E44.36)')," p2p3_QUAD ", p2p3_QUAD
cfc       Write(*,'(a5,E44.36)')," p2p4_QUAD ", p2p4_QUAD
cfc       Write(*,'(a5,E44.36)')," p2p5_QUAD ", p2p5_QUAD
cfc       Write(*,'(a5,E44.36)')," p3sq_QUAD ", p3sq_QUAD 
cfc       Write(*,'(a5,E44.36)')," p3p4_QUAD ", p3p4_QUAD
cfc       Write(*,'(a5,E44.36)')," p3p5_QUAD ", p3p5_QUAD
cfc       Write(*,'(a5,E44.36)')," p4sq_QUAD ", p4sq_QUAD 
cfc       Write(*,'(a5,E44.36)')," p4p5_QUAD ", p4p5_QUAD
cfc       Write(*,'(a5,E44.36)')," p5sq_QUAD ", p5sq_QUAD 

       p1sq=dble(p1sq_QUAD) 
c       p1p2=dble(p1p2_QUAD)
c       p1p3=dble(p1p3_QUAD)
c       p1p4=dble(p1p4_QUAD)
c       p1p5=dble(p1p5_QUAD)
       p2sq=dble(p2sq_QUAD) 
c       p2p3=dble(p2p3_QUAD)
c       p2p4=dble(p2p4_QUAD)
c       p2p5=dble(p2p5_QUAD)
       p3sq=dble(p3sq_QUAD) 
c       p3p4=dble(p3p4_QUAD)
c       p3p5=dble(p3p5_QUAD)
       p4sq=dble(p4sq_QUAD) 
c       p4p5=dble(p4p5_QUAD)
       p5sq=dble(p5sq_QUAD) 
       s12=dble(s12_QUAD) 
c       s13=dble(s13_QUAD) 
c       s14=dble(s14_QUAD) 
       s15=dble(s15_QUAD) 
       s23=dble(s23_QUAD) 
c       s24=dble(s24_QUAD) 
c       s25=dble(s25_QUAD) 
       s34=dble(s34_QUAD) 
c       s35=dble(s35_QUAD) 
       s45=dble(s45_QUAD) 

       Is12=1d0/s12
       Is45=1d0/s45
       Is12s45=1d0/(s12*s45)

       do i=0,3
         p1(i)=dble(p1_QUAD(i))
         p2(i)=dble(p2_QUAD(i))
         p3(i)=dble(p3_QUAD(i))
         p4(i)=dble(p4_QUAD(i))
         p5(i)=dble(p5_QUAD(i))

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
        
c        print*, "here:Pen2g",j,comp


       call B0fin_QUAD(p1sq_QUAD,musq_QUAD,B012_QUAD)
       B012R_QUAD=real(B012_QUAD)
       B012I_QUAD=imag(B012_QUAD)
       call From_BijR_QUAD_to_Dbl(B012R,B012I,B012R_QUAD,B012I_QUAD)
       B012=DCMPLX(B012R,B012I)

       call B0fin_QUAD(p2sq_QUAD,musq_QUAD,B023_QUAD)
       B023R_QUAD=real(B023_QUAD)
       B023I_QUAD=imag(B023_QUAD)
       call From_BijR_QUAD_to_Dbl(B023R,B023I,B023R_QUAD,B023I_QUAD)
       B023=DCMPLX(B023R,B023I)

       call B0fin_QUAD(p3sq_QUAD,musq_QUAD,B034_QUAD)
       B034R_QUAD=real(B034_QUAD)
       B034I_QUAD=imag(B034_QUAD)
       call From_BijR_QUAD_to_Dbl(B034R,B034I,B034R_QUAD,B034I_QUAD)
       B034=DCMPLX(B034R,B034I)

       call B0fin_QUAD(p4sq_QUAD,musq_QUAD,B045_QUAD)
       B045R_QUAD=real(B045_QUAD)
       B045I_QUAD=imag(B045_QUAD)
       call From_BijR_QUAD_to_Dbl(B045R,B045I,B045R_QUAD,B045I_QUAD)
       B045=DCMPLX(B045R,B045I)


       call B0fin_QUAD(s12_QUAD,musq_QUAD,B013_QUAD)
       B013R_QUAD=real(B013_QUAD)
       B013I_QUAD=imag(B013_QUAD)
       call From_BijR_QUAD_to_Dbl(B013R,B013I,B013R_QUAD,B013I_QUAD)
       B013=DCMPLX(B013R,B013I)

       call B0fin_QUAD(s45_QUAD,musq_QUAD,B014_QUAD)
       B014R_QUAD=real(B014_QUAD)
       B014I_QUAD=imag(B014_QUAD)
       call From_BijR_QUAD_to_Dbl(B014R,B014I,B014R_QUAD,B014I_QUAD)
       B014=DCMPLX(B014R,B014I)


       call B0fin_QUAD(s23_QUAD,musq_QUAD,B024_QUAD)
       B024R_QUAD=real(B024_QUAD)
       B024I_QUAD=imag(B024_QUAD)
       call From_BijR_QUAD_to_Dbl(B024R,B024I,B024R_QUAD,B024I_QUAD)
       B024=DCMPLX(B024R,B024I)

       call B0fin_QUAD(s15_QUAD,musq_QUAD,B025_QUAD)
       B025R_QUAD=real(B025_QUAD)
       B025I_QUAD=imag(B025_QUAD)
       call From_BijR_QUAD_to_Dbl(B025R,B025I,B025R_QUAD,B025I_QUAD)
       B025=DCMPLX(B025R,B025I)

       call B0fin_QUAD(s34_QUAD,musq_QUAD,B035_QUAD)
       B035R_QUAD=real(B035_QUAD)
       B035I_QUAD=imag(B035_QUAD)
       call From_BijR_QUAD_to_Dbl(B035R,B035I,B035R_QUAD,B035I_QUAD)
       B035=DCMPLX(B035R,B035I)


       call B0fin_QUAD(p5sq_QUAD,musq_QUAD,B015_QUAD)
       B015R_QUAD=real(B015_QUAD)
       B015I_QUAD=imag(B015_QUAD)
       call From_BijR_QUAD_to_Dbl(B015R,B015I,B015R_QUAD,B015I_QUAD)
       B015=DCMPLX(B015R,B015I)

       call C0fin_QUAD(p1sq_QUAD,p2sq_QUAD,s12_QUAD,musq_QUAD,C0123_QU
     -   AD)
       call tens_red3_new_Re_Com_QUAD(p1sq_QUAD,p2sq_QUAD,s12_QUAD,B02
     -   3_QUAD,B013_QUAD,B012_QUAD,C0123_QUAD,C0123R_QUAD,C0123I_QUA
     -   D,C123R_QUAD,C123I_QUAD)
       call  From_CijR_QUAD_to_Dbl(C0123R,C0123I,C123R,C123I,C0123R_QU
     -   AD,C0123I_QUAD,C123R_QUAD,C123I_QUAD)
       C0123=DCMPLX(C0123R,C0123I)


       call C0fin_QUAD(p1sq_QUAD,s23_QUAD,s45_QUAD,musq_QUAD,C0124_QU
     -   AD)
       call tens_red3_new_Re_Com_QUAD(p1sq_QUAD,s23_QUAD,s45_QUAD,B02
     -   4_QUAD,B014_QUAD,B012_QUAD,C0124_QUAD,C0124R_QUAD,C0124I_QUA
     -   D,C124R_QUAD,C124I_QUAD)
       call  From_CijR_QUAD_to_Dbl(C0124R,C0124I,C124R,C124I,C0124R_QU
     -   AD,C0124I_QUAD,C124R_QUAD,C124I_QUAD)
       C0124=DCMPLX(C0124R,C0124I)

       call C0fin_QUAD(p1sq_QUAD,s15_QUAD,p5sq_QUAD,musq_QUAD,C0125_QU
     -   AD)
       call tens_red3_new_Re_Com_QUAD(p1sq_QUAD,s15_QUAD,p5sq_QUAD,B02
     -   5_QUAD,B015_QUAD,B012_QUAD,C0125_QUAD,C0125R_QUAD,C0125I_QUA
     -   D,C125R_QUAD,C125I_QUAD)
       call  From_CijR_QUAD_to_Dbl(C0125R,C0125I,C125R,C125I,C0125R_QU
     -   AD,C0125I_QUAD,C125R_QUAD,C125I_QUAD)
       C0125=DCMPLX(C0125R,C0125I)

       call C0fin_QUAD(p2sq_QUAD,p3sq_QUAD,s23_QUAD,musq_QUAD,C0234_QU
     -   AD)
       call tens_red3_new_Re_Com_QUAD(p2sq_QUAD,p3sq_QUAD,s23_QUAD,B03
     -   4_QUAD,B024_QUAD,B023_QUAD,C0234_QUAD,C0234R_QUAD,C0234I_QUA
     -   D,C234R_QUAD,C234I_QUAD)
       call  From_CijR_QUAD_to_Dbl(C0234R,C0234I,C234R,C234I,C0234R_QU
     -   AD,C0234I_QUAD,C234R_QUAD,C234I_QUAD)
       C0234=DCMPLX(C0234R,C0234I)

       call C0fin_QUAD(p2sq_QUAD,s34_QUAD,s15_QUAD,musq_QUAD,C0235_QU
     -   AD)
       call tens_red3_new_Re_Com_QUAD(p2sq_QUAD,s34_QUAD,s15_QUAD,B03
     -   5_QUAD,B025_QUAD,B023_QUAD,C0235_QUAD,C0235R_QUAD,C0235I_QUA
     -   D,C235R_QUAD,C235I_QUAD)
       call  From_CijR_QUAD_to_Dbl(C0235R,C0235I,C235R,C235I,C0235R_QU
     -   AD,C0235I_QUAD,C235R_QUAD,C235I_QUAD)
       C0235=DCMPLX(C0235R,C0235I)


       call C0fin_QUAD(p3sq_QUAD,p4sq_QUAD,s34_QUAD,musq_QUAD,C0345_QU
     -   AD)
       call tens_red3_new_Re_Com_QUAD(p3sq_QUAD,p4sq_QUAD,s34_QUAD,B04
     -   5_QUAD,B035_QUAD,B034_QUAD,C0345_QUAD,C0345R_QUAD,C0345I_QUA
     -   D,C345R_QUAD,C345I_QUAD)
       call  From_CijR_QUAD_to_Dbl(C0345R,C0345I,C345R,C345I,C0345R_QU
     -   AD,C0345I_QUAD,C345R_QUAD,C345I_QUAD)
       C0345=DCMPLX(C0345R,C0345I)


       call C0fin_QUAD(s12_QUAD,p3sq_QUAD,s45_QUAD,musq_QUAD,C0134_QU
     -   AD)
       call tens_red3_new_Re_Com_QUAD(s12_QUAD,p3sq_QUAD,s45_QUAD,B03
     -   4_QUAD,B014_QUAD,B013_QUAD,C0134_QUAD,C0134R_QUAD,C0134I_QUA
     -   D,C134R_QUAD,C134I_QUAD)
       call  From_CijR_QUAD_to_Dbl(C0134R,C0134I,C134R,C134I,C0134R_QU
     -   AD,C0134I_QUAD,C134R_QUAD,C134I_QUAD)
       C0134=DCMPLX(C0134R,C0134I)

       call C0fin_QUAD(s12_QUAD,s34_QUAD,p5sq_QUAD,musq_QUAD,C0135_QUAD
     -   )
       call tens_red3_new_Re_Com_QUAD(s12_QUAD,s34_QUAD,p5sq_QUAD,B035_
     -   QUAD,B015_QUAD,B013_QUAD,C0135_QUAD,C0135R_QUAD,C0135I_QUAD,
     -   C135R_QUAD,C135I_QUAD)
       call  From_CijR_QUAD_to_Dbl(C0135R,C0135I,C135R,C135I,C0135R_QU
     -   AD,C0135I_QUAD,C135R_QUAD,C135I_QUAD)
       C0135=DCMPLX(C0135R,C0135I)


       call C0fin_QUAD(s45_QUAD,p4sq_QUAD,p5sq_QUAD,musq_QUAD,C0145_QU
     -   AD)
       call tens_red3_new_Re_Com_QUAD(s45_QUAD,p4sq_QUAD,p5sq_QUAD,B04
     -   5_QUAD,B015_QUAD,B014_QUAD,C0145_QUAD,C0145R_QUAD,C0145I_QUA
     -   D,C145R_QUAD,C145I_QUAD)
       call  From_CijR_QUAD_to_Dbl(C0145R,C0145I,C145R,C145I,C0145R_QU
     -   AD,C0145I_QUAD,C145R_QUAD,C145I_QUAD)
       C0145=DCMPLX(C0145R,C0145I)


       call C0fin_QUAD(s23_QUAD,p4sq_QUAD,s15_QUAD,musq_QUAD,C0245_QU
     -   AD)
       call tens_red3_new_Re_Com_QUAD(s23_QUAD,p4sq_QUAD,s15_QUAD,B04
     -   5_QUAD,B025_QUAD,B024_QUAD,C0245_QUAD,C0245R_QUAD,C0245I_QUA
     -   D,C245R_QUAD,C245I_QUAD)
       call  From_CijR_QUAD_to_Dbl(C0245R,C0245I,C245R,C245I,C0245R_QU
     -   AD,C0245I_QUAD,C245R_QUAD,C245I_QUAD)
       C0245=DCMPLX(C0245R,C0245I)


       call D0fin_QUAD(s12_QUAD,s23_QUAD,p1sq_QUAD,p2sq_QUAD,p3sq_QUAD
     -   ,s45_QUAD,musq_QUAD,D01234_QUAD)
       call tens_red4_new_Re_Com_QUAD(p1sq_QUAD,p2sq_QUAD,p3sq_QUAD,p1
     -   p2_QUAD,p1p3_QUAD,p2p3_QUAD,C0234R_QUAD,C0134R_QUAD,C0124R_Q
     -   UAD,C0123R_QUAD,C234R_QUAD,C134R_QUAD,C124R_QUAD,C123R_QUAD,
     -   C0234I_QUAD,C0134I_QUAD,C0124I_QUAD,C0123I_QUAD,C234I_QUAD,C
     -   134I_QUAD,C124I_QUAD,C123I_QUAD,D01234_QUAD,D01234R_QUAD,D01
     -   234I_QUAD,D1234R_QUAD,D1234I_QUAD)
       call  From_DijR_QUAD_to_Dbl(D01234R,D01234I,D1234R,D1234I,D0123
     -   4R_QUAD,D01234I_QUAD,D1234R_QUAD,D1234I_QUAD)
       D01234=DCMPLX(D01234R,D01234I)

       call D0fin_QUAD(s12_QUAD,s15_QUAD,p1sq_QUAD,p2sq_QUAD,s34_QUAD
     -   ,p5sq_QUAD,musq_QUAD,D01235_QUAD)
       call tens_red4_new_Re_Com_QUAD(p1sq_QUAD,p2sq_QUAD,s34_QUAD,p1p
     -   2_QUAD,p1p3_QUAD+ p1p4_QUAD,p2p3_QUAD+ p2p4_QUAD,C0235R_QUAD
     -   ,C0135R_QUAD,C0125R_QUAD,C0123R_QUAD,C235R_QUAD,C135R_QUAD,C
     -   125R_QUAD,C123R_QUAD,C0235I_QUAD,C0135I_QUAD,C0125I_QUAD,C01
     -   23I_QUAD,C235I_QUAD,C135I_QUAD,C125I_QUAD,C123I_QUAD,D01235_
     -   QUAD,D01235R_QUAD,D01235I_QUAD,D1235R_QUAD,D1235I_QUAD)
       call  From_DijR_QUAD_to_Dbl(D01235R,D01235I,D1235R,D1235I,D0123
     -   5R_QUAD,D01235I_QUAD,D1235R_QUAD,D1235I_QUAD)
       D01235=DCMPLX(D01235R,D01235I)


       call D0fin_QUAD(s45_QUAD,s15_QUAD,p1sq_QUAD,s23_QUAD,p4sq_QUA
     -   D,p5sq_QUAD,musq_QUAD,D01245_QUAD)
       call tens_red4_new_Re_Com_QUAD(p1sq_QUAD,s23_QUAD,p4sq_QUAD,p1p
     -   2_QUAD+ p1p3_QUAD,p1p4_QUAD,p2p4_QUAD+ p3p4_QUAD,C0245R_QUAD
     -   ,C0145R_QUAD,C0125R_QUAD,C0124R_QUAD,C245R_QUAD,C145R_QUAD,C
     -   125R_QUAD,C124R_QUAD,C0245I_QUAD,C0145I_QUAD,C0125I_QUAD,C01
     -   24I_QUAD,C245I_QUAD,C145I_QUAD,C125I_QUAD,C124I_QUAD,D01245_
     -   QUAD,D01245R_QUAD,D01245I_QUAD,D1245R_QUAD,D1245I_QUAD)
       call  From_DijR_QUAD_to_Dbl(D01245R,D01245I,D1245R,D1245I,D0124
     -   5R_QUAD,D01245I_QUAD,D1245R_QUAD,D1245I_QUAD)
       D01245=DCMPLX(D01245R,D01245I)


       call D0fin_QUAD(s23_QUAD,s34_QUAD,p2sq_QUAD,p3sq_QUAD,p4sq_QUAD
     -   ,s15_QUAD,musq_QUAD,D02345_QUAD)
       call tens_red4_new_Re_Com_QUAD(p2sq_QUAD,p3sq_QUAD,p4sq_QUAD,p2
     -   p3_QUAD,p2p4_QUAD,p3p4_QUAD,C0345R_QUAD,C0245R_QUAD,C0235R_Q
     -   UAD,C0234R_QUAD,C345R_QUAD,C245R_QUAD,C235R_QUAD,C234R_QUAD,
     -   C0345I_QUAD,C0245I_QUAD,C0235I_QUAD,C0234I_QUAD,C345I_QUAD,C
     -   245I_QUAD,C235I_QUAD,C234I_QUAD,D02345_QUAD,D02345R_QUAD,D02
     -   345I_QUAD,D2345R_QUAD,D2345I_QUAD)
       call  From_DijR_QUAD_to_Dbl(D02345R,D02345I,D2345R,D2345I,D0234
     -   5R_QUAD,D02345I_QUAD,D2345R_QUAD,D2345I_QUAD)
       D02345=DCMPLX(D02345R,D02345I)


      call D0fin_QUAD(s45_QUAD,s34_QUAD,s12_QUAD,p3sq_QUAD,p4sq_QUAD
     -   ,p5sq_QUAD,musq_QUAD,D01345_QUAD)
       call tens_red4_new_Re_Com_QUAD(s12_QUAD,p3sq_QUAD,p4sq_QUAD,p1p
     -   3_QUAD+ p2p3_QUAD,p1p4_QUAD+ p2p4_QUAD,p3p4_QUAD,C0345R_QUAD
     -   ,C0145R_QUAD,C0135R_QUAD,C0134R_QUAD,C345R_QUAD,C145R_QUAD,C
     -   135R_QUAD,C134R_QUAD,C0345I_QUAD,C0145I_QUAD,C0135I_QUAD,C01
     -   34I_QUAD,C345I_QUAD,C145I_QUAD,C135I_QUAD,C134I_QUAD,D01345_
     -   QUAD,D01345R_QUAD,D01345I_QUAD,D1345R_QUAD,D1345I_QUAD)
       call  From_DijR_QUAD_to_Dbl(D01345R,D01345I,D1345R,D1345I,D0134
     -   5R_QUAD,D01345I_QUAD,D1345R_QUAD,D1345I_QUAD)
       D01345=DCMPLX(D01345R,D01345I)


       call E0fin_QUAD(p1sq_QUAD,p2sq_QUAD,p3sq_QUAD,p4sq_QUAD,p5sq_QUA
     -   D,s12_QUAD,s23_QUAD,s34_QUAD,s45_QUAD,s15_QUAD,D02345_QUAD
     -   ,D01345_QUAD,D01245_QUAD,D01235_QUAD,D01234_QUAD,E012345_QUA
     -   D)
       E012345R_QUAD=real(E012345_QUAD)
       E012345I_QUAD=imag(E012345_QUAD)
       call tens_red5_new_Re_Com_QUAD(p1sq_QUAD,p2sq_QUAD,p3sq_QUAD,p4
     -   sq_QUAD,p1p2_QUAD,p1p3_QUAD,p1p4_QUAD,p2p3_QUAD,p2p4_QUAD,p3
     -   p4_QUAD,D02345R_QUAD,D01345R_QUAD,D01245R_QUAD,D01235R_QUAD,
     -   D01234R_QUAD,D2345R_QUAD,D1345R_QUAD,D1245R_QUAD,D1235R_QUAD
     -   ,D1234R_QUAD,D02345I_QUAD,D01345I_QUAD,D01245I_QUAD,D01235I_
     -   QUAD,D01234I_QUAD,D2345I_QUAD,D1345I_QUAD,D1245I_QUAD,D1235I
     -   _QUAD,D1234I_QUAD,E12345R_QUAD,E12345I_QUAD)
       call From_EijR_QUAD_to_Dbl(E012345R,E012345I,E12345R,E12345I,E0
     -   12345R_QUAD,E012345I_QUAD,E12345R_QUAD,E12345I_QUAD)
       E012345=DCMPLX(E012345R,E012345I)

c************************************************************************************
c************************************************************************************
c************************************************************************************
c   Printing C_ij,D_ij,E_ij Functions    
c************************************************************************************
c************************************************************************************
       If(PrintB) then
       write(*,*)"B12"
       write(*,*) "B12=", DCMPLX(B012R,B012I)
       write(*,*)"B23"
       write(*,*) "B23=", DCMPLX(B023R,B023I)
       write(*,*)"B34"
       write(*,*) "B34=", DCMPLX(B034R,B034I)
       write(*,*)"B45"
       write(*,*) "B45=", DCMPLX(B045R,B045I)
       write(*,*)"B13"
       write(*,*) "B13=", DCMPLX(B013R,B013I)
       write(*,*)"B14"
       write(*,*) "B14=", DCMPLX(B014R,B014I)
       write(*,*)"B24"
       write(*,*) "B24=", DCMPLX(B024R,B024I)
       write(*,*)"B25"
       write(*,*) "B25=", DCMPLX(B025R,B025I)
       write(*,*)"B35"
       write(*,*) "B35=", DCMPLX(B035R,B035I)
       write(*,*)"B15"
       write(*,*) "B15=", DCMPLX(B015R,B015I)
      Endif    
       If(PrintC) then
       write(*,*)"C123"
       call print_CijR_QUAD(C0123_QUAD,C123R_QUAD,C123I_QUAD)
       call print_CijR_QUAD_D(C0123,C123R,C123I,C0123_QUAD,C123R_QUAD,
     -   C123I_QUAD)
       call print_CijR(C0123,C123R,C123I)
cfc       call tens_Int3_Re_Com(0d0,p1sq,p2sq,s12,B023,B013,B012,C0123R,C
cfc     -   0123I,C123R,C123I)
       write(*,*)"C124"
       call print_CijR_QUAD(C0124_QUAD,C124R_QUAD,C124I_QUAD)
       call print_CijR_QUAD_D(C0124,C124R,C124I,C0124_QUAD,C124R_QUAD,
     -   C124I_QUAD)
       call print_CijR(C0124,C124R,C124I)
cfc       call tens_Int3_Re_Com(0d0,p1sq,s23,s45,B024,B014,B012,C0124R,C
cfc     -   0124I,C124R,C124I)
       write(*,*)"C125"
       call print_CijR_QUAD(C0125_QUAD,C125R_QUAD,C125I_QUAD)
       call print_CijR_QUAD_D(C0125,C125R,C125I,C0125_QUAD,C125R_QUAD,
     -   C125I_QUAD)
       call print_CijR(C0125,C125R,C125I)
cfc       call tens_Int3_Re_Com(0d0,p1sq,s15,p5sq,B025,B015,B012,C0125R,C
cfc     -   0125I,C125R,C125I)
       write(*,*)"C234"
       call print_CijR_QUAD(C0234_QUAD,C234R_QUAD,C234I_QUAD)
       call print_CijR_QUAD_D(C0234,C234R,C234I,C0234_QUAD,C234R_QUAD,
     -   C234I_QUAD)
       call print_CijR(C0234,C234R,C234I)
cfc       call tens_Int3_Re_Com(0d0,p2sq,p3sq,s23,B034,B024,B023,C0234R,C
cfc     -   0234I,C234R,C234I)
       write(*,*)"C235"
       call print_CijR_QUAD(C0235_QUAD,C235R_QUAD,C235I_QUAD)
       call print_CijR_QUAD_D(C0235,C235R,C235I,C0235_QUAD,C235R_QUAD,
     -   C235I_QUAD)
       call print_CijR(C0235,C235R,C235I)
cfc       call tens_Int3_Re_Com(0d0,p2sq,s34,s15,B035,B025,B023,C0235R,C
cfc     -   0235I,C235R,C235I)
       write(*,*)"C345"
       call print_CijR_QUAD(C0345_QUAD,C345R_QUAD,C345I_QUAD)
       call print_CijR_QUAD_D(C0345,C345R,C345I,C0345_QUAD,C345R_QUAD,
     -   C345I_QUAD)
       call print_CijR(C0345,C345R,C345I)
cfc       call tens_Int3_Re_Com(0d0,p3sq,p4sq,s34,B045,B035,B034,C0345R,C
cfc     -   0345I,C345R,C345I)
       write(*,*)"C134"
       call print_CijR_QUAD(C0134_QUAD,C134R_QUAD,C134I_QUAD)
       call print_CijR_QUAD_D(C0134,C134R,C134I,C0134_QUAD,C134R_QUAD,
     -   C134I_QUAD)
       call print_CijR(C0134,C134R,C134I)
cfc       call tens_Int3_Re_Com(0d0,s12,p3sq,s45,B034,B014,B013,C0134R,C
cfc     -   0134I,C134R,C134I)
       write(*,*)"C135"
       call print_CijR_QUAD(C0135_QUAD,C135R_QUAD,C135I_QUAD)
       call print_CijR_QUAD_D(C0135,C135R,C135I,C0135_QUAD,C135R_QUAD,
     -   C135I_QUAD)
       call print_CijR(C0135,C135R,C135I)
cfc       call tens_Int3_Re_Com(0d0,s12,s34,p5sq,B035,B015,B013,C0135R,C01
cfc     -   35I,C135R,C135I)
       write(*,*)"C145"
       call print_CijR_QUAD(C0145_QUAD,C145R_QUAD,C145I_QUAD)
       call print_CijR_QUAD_D(C0145,C145R,C145I,C0145_QUAD,C145R_QUAD,
     -   C145I_QUAD)
       call print_CijR(C0145,C145R,C145I)
cfc       call tens_Int3_Re_Com(0d0,s45,p4sq,p5sq,B045,B015,B014,C0145R,C
cfc     -   0145I,C145R,C145I)
       write(*,*)"C245"
       call print_CijR_QUAD(C0245_QUAD,C245R_QUAD,C245I_QUAD)
       call print_CijR_QUAD_D(C0245,C245R,C245I,C0245_QUAD,C245R_QUAD,
     -   C245I_QUAD)
       call print_CijR(C0245,C245R,C245I)
cfc       call tens_Int3_Re_Com(0d0,s23,p4sq,s15,B045,B025,B024,C0245R,C
cfc     -   0245I,C245R,C245I)
      Endif    

       If(PrintD) then
       write(*,*)"D1234"
       call print_DijR_QUAD(D01234_QUAD,D1234R_QUAD,D1234I_QUAD)
       call print_DijR_QUAD_D(D01234,D1234R,D1234I,D01234_QUAD,D1234R_
     -   QUAD,D1234I_QUAD)
       call print_DijR(D01234,D1234R,D1234I)
cfc       call tens_Int4_Re_Com(0d0,p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,C0234R,
cfc     -   C0134R,C0124R,C0123R,C234R,C134R,C124R,C123R,C0234I,C0134I,C
cfc     -   0124I,C0123I,C234I,C134I,C124I,C123I,D01234R,D01234I,D1234R,
cfc     -   D1234I)
       write(*,*)"D1235"
       call print_DijR_QUAD(D01235_QUAD,D1235R_QUAD,D1235I_QUAD)
       call print_DijR_QUAD_D(D01235,D1235R,D1235I,D01235_QUAD,D1235R_
     -   QUAD,D1235I_QUAD)
       call print_DijR(D01235,D1235R,D1235I)
cfc       call tens_Int4_Re_Com(0d0,p1sq,p2sq,s34,p1p2,p1p3+p1p4,p2p3+p2p
cfc     -   4,C0235R,C0135R,C0125R,C0123R,C235R,C135R,C125R,C123R,C0235I
cfc     -   ,C0135I,C0125I,C0123I,C235I,C135I,C125I,C123I,D01235R,D01235
cfc     -   I,D1235R,D1235I)
       write(*,*)"D1245"
       call print_DijR_QUAD(D01245_QUAD,D1245R_QUAD,D1245I_QUAD)
       call print_DijR_QUAD_D(D01245,D1245R,D1245I,D01245_QUAD,D1245R_
     -   QUAD,D1245I_QUAD)
       call print_DijR(D01245,D1245R,D1245I)
cfc       call tens_Int4_Re_Com(0d0,p1sq,s23,p4sq,p1p2+p1p3,p1p4,p2p4+p3p
cfc     -   4,C0245R,C0145R,C0125R,C0124R,C245R,C145R,C125R,C124R,C0245I
cfc     -   ,C0145I,C0125I,C0124I,C245I,C145I,C125I,C124I,D01245R,D01245
cfc     -   I,D1245R,D1245I)
       write(*,*)"D2345"
       call print_DijR_QUAD(D02345_QUAD,D2345R_QUAD,D2345I_QUAD)
       call print_DijR_QUAD_D(D02345,D2345R,D2345I,D02345_QUAD,D2345R_
     -   QUAD,D2345I_QUAD)
       call print_DijR(D02345,D2345R,D2345I)
cfc       call tens_Int4_Re_Com(0d0,p2sq,p3sq,p4sq,p2p3,p2p4,p3p4,C0345R,
cfc     -   C0245R,C0235R,C0234R,C345R,C245R,C235R,C234R,C0345I,C0245I,C
cfc     -   0235I,C0234I,C345I,C245I,C235I,C234I,D02345R,D02345I,D2345R,
cfc     -   D2345I)
       write(*,*)"D1345"
       call print_DijR_QUAD(D01345_QUAD,D1345R_QUAD,D1345I_QUAD)
       call print_DijR_QUAD_D(D01345,D1345R,D1345I,D01345_QUAD,D1345R_
     -   QUAD,D1345I_QUAD)
       call print_DijR(D01345,D1345R,D1345I)
cfc       call tens_Int4_Re_Com(0d0,s12,p3sq,p4sq,p1p3+p2p3,p1p4+p2p4,p3p
cfc     -   4,C0345R,C0145R,C0135R,C0134R,C345R,C145R,C135R,C134R,C0345I
cfc     -   ,C0145I,C0135I,C0134I,C345I,C145I,C135I,C134I,D01345R,D01345
cfc     -   I,D1345R,D1345I)
      Endif    
       If(PrintE) then
       write(*,*)"E12345"
       call print_EijR_QUAD(E012345_QUAD,E12345R_QUAD,E12345I_QUAD)
       call print_EijR_QUAD_D(E012345,E12345R,E12345I,E012345,E12345R_
     -   QUAD,E12345I_QUAD)
       call print_EijR(E012345,E12345R,E12345I)
cfc       call tens_Int5_Re_Com(0d0,p1sq,p2sq,p3sq,p4sq,p1p2,p1p3,p1p4,p2
cfc     -   p3,p2p4,p3p4,D02345R,D01345R,D01245R,D01235R,D01234R,D2345R,
cfc     -   D1345R,D1245R,D1235R,D1234R,D02345I,D01345I,D01245I,D01235I,
cfc     -   D01234I,D2345I,D1345I,D1245I,D1235I,D1234I,E012345R,E012345I
cfc     -   ,E12345R,E12345I)
      Endif    
c*****************

****************************************************************
       Return
       End
