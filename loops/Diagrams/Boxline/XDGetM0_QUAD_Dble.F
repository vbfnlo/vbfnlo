       subroutine PrintIntBoxMassless_QUAD(p1,p2,p3,p4,musq)
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
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      Real*8 psq(4)
        Real*8   p1sq
       Real*8   p2sq
       Real*8   p3sq
       Real*8   p4sq

       Real*8   s12
       Real*8   s23

       Complex*16  B012,B013,B014
       Complex*16  B023,B024
       Complex*16  B034

       Real*8  B012R,B013R,B014R
       Real*8  B023R,B024R
       Real*8  B034R

       Real*8  B012I,B013I,B014I
       Real*8  B023I,B024I
       Real*8  B034I

       Complex*16 C0123,C0124
       Complex*16 C0134

       Complex*16 C0234

       Real*8 C0123R,C0124R
       Real*8 C0134R

       Real*8 C0234R

       Real*8 C0123I,C0124I
       Real*8 C0134I

       Real*8 C0234I

       Real*8 C123R(4,2),C124R(4,2)
       Real*8 C134R(4,2)

       Real*8 C234R(4,2)

       Real*8 C123I(4,2),C124I(4,2)
       Real*8 C134I(4,2)

       Real*8 C234I(4,2)

       Complex*16  D01234

       Real*8 D01234R
       Real*8 D01234I

       Real*8 D1234R(13,3)
       Real*8 D1234I(13,3)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  QUAD precision start here
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       Real*16  p1_QUAD(0:3),p2_QUAD(0:3),p3_QUAD(0:3),p4_QUAD(0:3)
       Real*16 ptemp(0:3,4)
       Real*16   p1sq_QUAD, p1p2_QUAD, p1p3_QUAD, p1p4_QUAD
       Real*16   p2sq_QUAD, p2p3_QUAD, p2p4_QUAD
       Real*16   p3sq_QUAD, p3p4_QUAD
       Real*16   p4sq_QUAD

       Real*16   s12_QUAD
       Real*16   s23_QUAD

       Complex*32  B012_QUAD,B013_QUAD,B014_QUAD
       Complex*32  B023_QUAD,B024_QUAD
       Complex*32  B034_QUAD

       Real*16  B012R_QUAD,B013R_QUAD,B014R_QUAD
       Real*16  B023R_QUAD,B024R_QUAD
       Real*16  B034R_QUAD

       Real*16  B012I_QUAD,B013I_QUAD,B014I_QUAD
       Real*16  B023I_QUAD,B024I_QUAD
       Real*16  B034I_QUAD

       Complex*32    C0123_QUAD,C0124_QUAD
       Complex*32    C0134_QUAD

       Complex*32    C0234_QUAD
       Real*16       C0123R_QUAD,C0124R_QUAD
       Real*16       C0134R_QUAD

       Real*16       C0234R_QUAD
       Real*16       C0123I_QUAD,C0124I_QUAD
       Real*16       C0134I_QUAD

       Real*16       C0234I_QUAD
       Real*16       C123R_QUAD(4,2),C124R_QUAD(4,2)
       Real*16       C134R_QUAD(4,2)

       Real*16       C234R_QUAD(4,2)
       Real*16       C123I_QUAD(4,2),C124I_QUAD(4,2)
       Real*16       C134I_QUAD(4,2)

       Real*16       C234I_QUAD(4,2)
       Complex*32  D01234_QUAD
       Real*16     D01234R_QUAD
       Real*16     D01234I_QUAD
       Real*16     D1234R_QUAD(13,3)
       Real*16     D1234I_QUAD(13,3)

       Real*16 musq_QUAD, dotrr_QUAD
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      End QUAD precision
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       COMMON/DVALFunctionsD/ D01234R,D1234R,D01234I,D1234I

       COMMON/CVALFunctionsD/ C0123R,C123R,C0123I,C123I, 
     &    C0124R,C124R,C0124I,C124I,
     &    C0134R,C134R,C0134I,C134I,
     &    C0234R,C234R,C0234I,C234I

       COMMON/BVALFunctionsD/ B012R,B012I, B013R,B013I, B014R,B014I, 
     &      B023R,B023I,  B024R,B024I,
     &    B034R,B034I
       COMMON/InvariantsD/p1sq,p2sq,p3sq,p4sq,s12,s23
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       Logical PrintB,PrintC,PrintD
        External dotrr_QUAD ,dotrr
       Complex*16 dotrc,dotcc
       Real*8 dotrr
       Real*8 musq
      EXTERNAL   dotrc,dotcc
       Integer i,j,jtemp
 

c      Write(*,*) 'musq',musq_QUAD
      PrintB=.false. 
      PrintC=.false. 
      PrintD=.false. 
c      PrintF=.false.


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


       if (psq(1).gt.10d-7) then
       j=1
       elseif (psq(2).gt.10d-7)then
       j=2
       elseif (psq(3).gt.10d-7)then
       j=3
       elseif (psq(4).gt.10d-7)then
 
       else
       j=7
c$$$       write(*,*) 'No massive invariant'
c$$$       Return
       endif
       if(j.ne.7) then
       do jtemp=1,4
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
      ptemp(i,j)=0q0
      do jtemp=1,4
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
      enddo 
      endif


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
       s23_QUAD = (p2sq_QUAD +p3sq_QUAD+ 2*p2p3_QUAD) 

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
       s12=dble(s12_QUAD) 
       s23=dble(s23_QUAD) 

       do i=0,3
         p1(i)=dble(p1_QUAD(i))
         p2(i)=dble(p2_QUAD(i))
         p3(i)=dble(p3_QUAD(i))
         p4(i)=dble(p4_QUAD(i))
      enddo

      If(abs(musq-s12).lt.10d-7) then
          musq_QUAD=s12_QUAD
          elseIf(abs(musq-s23).lt.10d-7) then
          musq_QUAD=s23_QUAD
          else
          musq_QUAD=musq
          endif

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


       call B0fin_QUAD(s12_QUAD,musq_QUAD,B013_QUAD)
       B013R_QUAD=real(B013_QUAD)
       B013I_QUAD=imag(B013_QUAD)
       call From_BijR_QUAD_to_Dbl(B013R,B013I,B013R_QUAD,B013I_QUAD)
       B013=DCMPLX(B013R,B013I)

       call B0fin_QUAD(p4sq_QUAD,musq_QUAD,B014_QUAD)
       B014R_QUAD=real(B014_QUAD)
       B014I_QUAD=imag(B014_QUAD)
       call From_BijR_QUAD_to_Dbl(B014R,B014I,B014R_QUAD,B014I_QUAD)
       B014=DCMPLX(B014R,B014I)


       call B0fin_QUAD(s23_QUAD,musq_QUAD,B024_QUAD)
       B024R_QUAD=real(B024_QUAD)
       B024I_QUAD=imag(B024_QUAD)
       call From_BijR_QUAD_to_Dbl(B024R,B024I,B024R_QUAD,B024I_QUAD)
       B024=DCMPLX(B024R,B024I)


       call C0fin_QUAD(p1sq_QUAD,p2sq_QUAD,s12_QUAD,musq_QUAD,C0123_QU
     -   AD)
       call tens_red3_new_Re_Com_QUAD(p1sq_QUAD,p2sq_QUAD,s12_QUAD,B02
     -   3_QUAD,B013_QUAD,B012_QUAD,C0123_QUAD,C0123R_QUAD,C0123I_QUA
     -   D,C123R_QUAD,C123I_QUAD)
       call  From_CijR_QUAD_to_Dbl(C0123R,C0123I,C123R,C123I,C0123R_QU
     -   AD,C0123I_QUAD,C123R_QUAD,C123I_QUAD)
       C0123=DCMPLX(C0123R,C0123I)


       call C0fin_QUAD(p1sq_QUAD,s23_QUAD,p4sq_QUAD,musq_QUAD,C0124_QU
     -   AD)
       call tens_red3_new_Re_Com_QUAD(p1sq_QUAD,s23_QUAD,p4sq_QUAD,B02
     -   4_QUAD,B014_QUAD,B012_QUAD,C0124_QUAD,C0124R_QUAD,C0124I_QUA
     -   D,C124R_QUAD,C124I_QUAD)
       call  From_CijR_QUAD_to_Dbl(C0124R,C0124I,C124R,C124I,C0124R_QU
     -   AD,C0124I_QUAD,C124R_QUAD,C124I_QUAD)
       C0124=DCMPLX(C0124R,C0124I)


       call C0fin_QUAD(p2sq_QUAD,p3sq_QUAD,s23_QUAD,musq_QUAD,C0234_QU
     -   AD)
       call tens_red3_new_Re_Com_QUAD(p2sq_QUAD,p3sq_QUAD,s23_QUAD,B03
     -   4_QUAD,B024_QUAD,B023_QUAD,C0234_QUAD,C0234R_QUAD,C0234I_QUA
     -   D,C234R_QUAD,C234I_QUAD)
       call  From_CijR_QUAD_to_Dbl(C0234R,C0234I,C234R,C234I,C0234R_QU
     -   AD,C0234I_QUAD,C234R_QUAD,C234I_QUAD)
       C0234=DCMPLX(C0234R,C0234I)


       call C0fin_QUAD(s12_QUAD,p3sq_QUAD,p4sq_QUAD,musq_QUAD,C0134_QU
     -   AD)
       call tens_red3_new_Re_Com_QUAD(s12_QUAD,p3sq_QUAD,p4sq_QUAD,B03
     -   4_QUAD,B014_QUAD,B013_QUAD,C0134_QUAD,C0134R_QUAD,C0134I_QUA
     -   D,C134R_QUAD,C134I_QUAD)
       call  From_CijR_QUAD_to_Dbl(C0134R,C0134I,C134R,C134I,C0134R_QU
     -   AD,C0134I_QUAD,C134R_QUAD,C134I_QUAD)
       C0134=DCMPLX(C0134R,C0134I)


       call D0fin_QUAD(s12_QUAD,s23_QUAD,p1sq_QUAD,p2sq_QUAD,p3sq_QUAD
     -   ,p4sq_QUAD,musq_QUAD,D01234_QUAD)
       call tens_red4_new_Re_Com_QUAD(p1sq_QUAD,p2sq_QUAD,p3sq_QUAD,p1
     -   p2_QUAD,p1p3_QUAD,p2p3_QUAD,C0234R_QUAD,C0134R_QUAD,C0124R_Q
     -   UAD,C0123R_QUAD,C234R_QUAD,C134R_QUAD,C124R_QUAD,C123R_QUAD,
     -   C0234I_QUAD,C0134I_QUAD,C0124I_QUAD,C0123I_QUAD,C234I_QUAD,C
     -   134I_QUAD,C124I_QUAD,C123I_QUAD,D01234_QUAD,D01234R_QUAD,D01
     -   234I_QUAD,D1234R_QUAD,D1234I_QUAD)
       call  From_DijR_QUAD_to_Dbl(D01234R,D01234I,D1234R,D1234I,D0123
     -   4R_QUAD,D01234I_QUAD,D1234R_QUAD,D1234I_QUAD)
       D01234=DCMPLX(D01234R,D01234I)


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
       write(*,*)"B13"
       write(*,*) "B13=", DCMPLX(B013R,B013I)
       write(*,*)"B14"
       write(*,*) "B14=", DCMPLX(B014R,B014I)
       write(*,*)"B24"
       write(*,*) "B24=", DCMPLX(B024R,B024I)
      Endif    
       If(PrintC) then
       write(*,*)"C123"
       call print_CijR_QUAD(C0123_QUAD,C123R_QUAD,C123I_QUAD)
       call print_CijR_QUAD_D(C0123,C123R,C123I,C0123_QUAD,C123R_QUAD,
     -   C123I_QUAD)
       call print_CijR(C0123,C123R,C123I)
       write(*,*)"C124"
       call print_CijR_QUAD(C0124_QUAD,C124R_QUAD,C124I_QUAD)
       call print_CijR_QUAD_D(C0124,C124R,C124I,C0124_QUAD,C124R_QUAD,
     -   C124I_QUAD)
       call print_CijR(C0124,C124R,C124I)
       write(*,*)"C234"
       call print_CijR_QUAD(C0234_QUAD,C234R_QUAD,C234I_QUAD)
       call print_CijR_QUAD_D(C0234,C234R,C234I,C0234_QUAD,C234R_QUAD,
     -   C234I_QUAD)
       call print_CijR(C0234,C234R,C234I)
       write(*,*)"C134"
       call print_CijR_QUAD(C0134_QUAD,C134R_QUAD,C134I_QUAD)
       call print_CijR_QUAD_D(C0134,C134R,C134I,C0134_QUAD,C134R_QUAD,
     -   C134I_QUAD)
       call print_CijR(C0134,C134R,C134I)
      Endif    

       If(PrintD) then
       write(*,*)"D1234"
       call print_DijR_QUAD(D01234_QUAD,D1234R_QUAD,D1234I_QUAD)
       call print_DijR_QUAD_D(D01234,D1234R,D1234I,D01234_QUAD,D1234R_
     -   QUAD,D1234I_QUAD)
       call print_DijR(D01234,D1234R,D1234I)
      Endif    
c*****************

****************************************************************
       Return
       End
