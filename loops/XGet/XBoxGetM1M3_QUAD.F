c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
       subroutine XDget_1M_QUAD(M1,p1sq,p2sq,p3sq,p1p2
     & ,p1p3,p2p3,musq)
       Implicit none
       Real*16   p1sq, p1p2, p1p3!, p1p4
       Real*16   p2sq, p2p3!, p2p4
       Real*16   p3sq
       Real*16   p4sq
       Real*16   s12!, s13!, s14
       Real*16   s23!, s24
!       Real*16   s34
       Real*16 M5,M1,M2,M3,M4,musq
       complex*32 A01, A02,A03,A04
       Real*16      A01R, A02R,A03R,A04R
       Real*16      A01I, A02I,A03I,A04I
       complex*32  B012,B013,B014 
       complex*32  B023,B024 
       complex*32  B034 

       Real*16  B012R,B013R,B014R 
       Real*16  B023R,B024R 
       Real*16  B034R 

       Real*16  B012I,B013I,B014I 
       Real*16  B023I,B024I 
       Real*16  B034I 

       Real*16  B12R(3,3),B13R(3,3),B14R(3,3) 
       Real*16  B23R(3,3),B24R(3,3) 
       Real*16  B34R(3,3) 

       Real*16  B12I(3,3),B13I(3,3),B14I(3,3) 
       Real*16  B23I(3,3),B24I(3,3) 
       Real*16  B34I(3,3) 


       complex*32 C0123,C0124 
       complex*32 C0134 

       complex*32 C0234 


       Real*16 C0123R,C0124R 
       Real*16 C0134R 

       Real*16 C0234R 


       Real*16 C0123I,C0124I 
       Real*16 C0134I 

       Real*16 C0234I 


       Real*16 C123R(9,4),C124R(9,4) 
       Real*16 C134R(9,4) 

       Real*16 C234R(9,4) 


       Real*16 C123I(9,4),C124I(9,4) 
       Real*16 C134I(9,4) 

       Real*16 C234I(9,4) 


       complex*32  D01234


       Real*16 D01234R


       Real*16 D01234I


       Real*16 D1234R(34,5)


       Real*16 D1234I(34,5)


       Logical PrintB,PrintC,PrintD 

       COMMON/DVALFunctionsGSim/ D01234R_d,D1234R_d,D01234I_d,
     &                D1234I_d
       COMMON/CVALFunctionsGSim/ C0123R_d,C123R_d,C0123I_d,C123I_d, 
     &             C0124R_d,C124R_d,C0124I_d,
     -             C124I_d, 
     &             C0134R_d,C134R_d,C0134I_d,
     -             C134I_d, 
     -             C0234R_d,C234R_d,C0234I_d,C234I_d
       COMMON/BVALFunctionsGSim/ 
     &           B012R_d,B012I_d, B013R_d,B013I_d, 
     &           B014R_d,B014I_d,
     -           B023R_d,B023I_d, B024R_d,B024I_d, 
     &           B034R_d,B034I_d,
     -           B12R_d,B12I_d, B13R_d,
     &           B13I_d, B14R_d,
     -           B14I_d,  B23R_d,B23I_d, B24R_d,B24I_d, 
     -           B34R_d, B34I_d


       Real*8  B012R_d,B013R_d,B014R_d 
       Real*8  B023R_d,B024R_d 
       Real*8  B034R_d 

       Real*8  B012I_d,B013I_d,B014I_d 
       Real*8  B023I_d,B024I_d 
       Real*8  B034I_d 

       Real*8  B12R_d(2,2),B13R_d(2,2),B14R_d(2,2) 
       Real*8  B23R_d(2,2),B24R_d(2,2) 
       Real*8  B34R_d(2,2) 

       Real*8  B12I_d(2,2),B13I_d(2,2),B14I_d(2,2) 
       Real*8  B23I_d(2,2),B24I_d(2,2) 
       Real*8  B34I_d(2,2)


       Real*8 C0123R_d,C0124R_d 
       Real*8 C0134R_d 

       Real*8 C0234R_d 


       Real*8 C0123I_d,C0124I_d 
       Real*8 C0134I_d 

       Real*8 C0234I_d 


       Real*8 C123R_d(6,3),C124R_d(6,3) 
       Real*8 C134R_d(6,3) 

       Real*8 C234R_d(6,3) 


       Real*8 C123I_d(6,3),C124I_d(6,3) 
       Real*8 C134I_d(6,3) 

       Real*8 C234I_d(6,3)


       Real*8 D01234R_d


       Real*8 D01234I_d


       Real*8 D1234R_d(22,4)


       Real*8 D1234I_d(22,4)


       complex*32 A0finG_QUAD,B0finG_QUAD,C0finG_QUAD,D0finG_QUAD

       external A0finG_QUAD,B0finG_QUAD,C0finG_QUAD,D0finG_QUAD
 

       PRINTb=.TRUE.
       PRINTc=.TRUE.
       PRINTd=.TRUE.

       M5=M1
       M4=M1
       M3=M1
       M2=M1


cfc       PRINTe=.TRUE.

c       print*, "here_XGETEM1M3"

       p4sq=p1sq+p2sq+p3sq+2*(p1p2+p1p3+p2p3)
       s12 = (p1sq +p2sq+ 2*p1p2) 
!       s13 = (p1sq +p3sq+ 2*p1p3) 
!       s14 = (p1sq +p4sq+ 2*p1p4) 
       s23 = (p2sq +p3sq+ 2*p2p3) 
!       s24 = (p2sq +p4sq+ 2*p2p4) 
!       s34 = (p3sq +p4sq+ 2*p3p4) 


        A01=A0finG_QUAD(M5,musq)  
        A01R=real(A01)   
        A01I=imag(A01)   
c
        A02=A0finG_QUAD(M1,musq)  
        A02R=real(A02)   
        A02I=imag(A02)   
c
        A03=A0finG_QUAD(M2,musq)  
        A03R=real(A03)   
        A03I=imag(A03)   
c
        A04=A0finG_QUAD(M3,musq)  
        A04R=real(A04)   
        A04I=imag(A04)  

c
       B012=B0finG_QUAD(M5,M1,p1sq,musq)
       B012R=real(B012)
       B012I=imag(B012)
       call tens_red2_R3_G_QUAD(M5,M1,p1sq,A02,A02R,A02I,A01,A01R,A01I,B012
     -   ,B012R,B012I,B12R,B12I) 
c
       B023=B0finG_QUAD(M1,M2,p2sq,musq)
       B023R=real(B023)
       B023I=imag(B023)
       call tens_red2_R3_G_QUAD(M1,M2,p2sq,A03,A03R,A03I,A02,A02R,A02I,B023
     -   ,B023R,B023I,B23R,B23I) 

       B034=B0finG_QUAD(M2,M3,p3sq,musq)
       B034R=real(B034)
       B034I=imag(B034)
       call tens_red2_R3_G_QUAD(M2,M3,p3sq,A04,A04R,A04I,A03,A03R,A03I,B034
     -   ,B034R,B034I,B34R,B34I) 

        B013=B0finG_QUAD(M5,M2,s12,musq)
       B013R=real(B013)
       B013I=imag(B013)
       call tens_red2_R3_G_QUAD(M5,M2,s12,A03,A03R,A03I,A01,A01R,A01I,B013,
     -   B013R,B013I,B13R,B13I) 

       B014=B0finG_QUAD(M5,M3,p4sq,musq)
       B014R=real(B014)
       B014I=imag(B014)
       call tens_red2_R3_G_QUAD(M5,M3,p4sq,A04,A04R,A04I,A01,A01R,A01I,B014
     -   ,B014R,B014I,B14R,B14I) 

       B024=B0finG_QUAD(M1,M3,s23,musq)
       B024R=real(B024)
       B024I=imag(B024)
       call tens_red2_R3_G_QUAD(M1,M3,s23,A04,A04R,A04I,A02,A02R,A02I,B024,
     -   B024R,B024I,B24R,B24I) 

   
c
c       print*, "here_XGETEM1M3after B"
c
       C0123=C0finG_QUAD(M5,M1,M2,p1sq,p2sq,s12,musq)
       call tens_red3_new_Re_Com_R4_G_QUAD(M5,M1,M2,p1sq,p2sq,s12,
     &  B023R,B013R,B012R,
     &  B23R,B13R,B12R,
     &  B023I,B013I,B012I,
     &  B23I,B13I,B12I,
     &  C0123,C0123R,C0123I,C123R,C123I)

       C0124=C0finG_QUAD(M5,M1,M3,p1sq,s23,p4sq,musq)
       call tens_red3_new_Re_Com_R4_G_QUAD(M5,M1,M3,p1sq,s23,p4sq,
     & B024R,B014R,B012R,
     & B24R,B14R,B12R,
     & B024I,B014I,B012I,
     & B24I,B14I,B12I,
     & C0124,C0124R,C0124I,C124R,C124I)

       C0234=C0finG_QUAD(M1,M2,M3,p2sq,p3sq,s23,musq)
       call tens_red3_new_Re_Com_R4_G_QUAD(M1,M2,M3,p2sq,p3sq,s23,
     & B034R,B024R,B023R,
     & B34R,B24R,B23R,
     & B034I,B024I,B023I,
     & B34I,B24I,B23I,
     & C0234,C0234R,C0234I,C234R,C234I)


       C0134=C0finG_QUAD(M5,M2,M3,s12,p3sq,p4sq,musq)
       call tens_red3_new_Re_Com_R4_G_QUAD(M5,M2,M3,s12,p3sq,p4sq,
     & B034R,B014R,B013R,
     & B34R,B14R,B13R,
     & B034I,B014I,B013I,
     & B34I,B14I,B13I,
     & C0134,C0134R,C0134I,C134R,C134I)


c
c       print*, "here_XGETEM1M3after C"
c
       D01234=D0finG_QUAD(M5,M1,M2,M3,s12,s23,p1sq,p2sq,p3sq,p4sq,musq)
       call tens_red4_new_Re_Com_R5_G_QUAD(M5,M1,M2,M3,p1sq,p2sq,p3sq,p1p2
     - ,p1p3,p2p3,C02
     -   34R,C0134R,C0124R,C0123R,C234R,C134R,C124R,C123R,C0234I,C013
     -   4I,C0124I,C0123I,C234I,C134I,C124I,C123I,D01234,D01234R,D012
     -   34I,D1234R,D1234I)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                    B0
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       call R3_From_BijR_QUAD_to_DblR2(B012R_d,B012I_d,B12R_d,B12I_d,
     &    B012R,B012I,B12R,B12I)
       call R3_From_BijR_QUAD_to_DblR2(B013R_d,B013I_d,B13R_d,B13I_d,
     &    B013R,B013I,B13R,B13I)
       call R3_From_BijR_QUAD_to_DblR2(B014R_d,B014I_d,B14R_d,B14I_d,
     &    B014R,B014I,B14R,B14I)
       call R3_From_BijR_QUAD_to_DblR2(B023R_d,B023I_d,B23R_d,B23I_d,
     &    B023R,B023I,B23R,B23I)
       call R3_From_BijR_QUAD_to_DblR2(B024R_d,B024I_d,B24R_d,B24I_d,
     &    B024R,B024I,B24R,B24I)
       call R3_From_BijR_QUAD_to_DblR2(B034R_d,B034I_d,B34R_d,B34I_d,
     &    B034R,B034I,B34R,B34I)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                    C0
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call  R4_From_CijR_QUAD_to_DblR3(C0123R_d,C0123I_d,C123R_d,C123I_d,
     -              C0123R,C0123I,C123R,C123I)
      call  R4_From_CijR_QUAD_to_DblR3(C0124R_d,C0124I_d,C124R_d,C124I_d,
     -              C0124R,C0124I,C124R,C124I)
      call  R4_From_CijR_QUAD_to_DblR3(C0134R_d,C0134I_d,C134R_d,C134I_d,
     -              C0134R,C0134I,C134R,C134I)
 
      call  R4_From_CijR_QUAD_to_DblR3(C0234R_d,C0234I_d,C234R_d,C234I_d,
     -              C0234R,C0234I,C234R,C234I)
 

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                D0
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call  R5_From_DijR_QUAD_to_DblR4(D01234R_d,D01234I_d,D1234R_d,D1234I_d,
     -            D01234R,D01234I,D1234R,D1234I)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                E0
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  

cfc       If(PrintB) then
cfc       write(*,*)"B12"
cfc       write(*,*) "B12=", DCMPLX(B012)
cfc       endif
cfc       If(PrintB) then
cfc       write(*,*)"B23"
cfc       write(*,*) "B23=", DCMPLX(B023)
cfc       endif
cfc       If(PrintB) then
cfc       write(*,*)"B34"
cfc       write(*,*) "B34=", DCMPLX(B034)
cfc       endif
cfc
cfc       If(PrintB) then
cfc       write(*,*)"B13"
cfc       write(*,*) "B13=", DCMPLX(B013)
cfc       endif
cfc       If(PrintB) then
cfc       write(*,*)"B14"
cfc       write(*,*) "B14=", DCMPLX(B014)
cfc       endif
cfc       If(PrintB) then
cfc       write(*,*)"B24"
cfc       write(*,*) "B24=", DCMPLX(B024)
cfc       endif
cfc
cfc       If(PrintC) then
cfc       write(*,*)"C123",C0123R_d,C0123I_d
cfc       
cfc       call print_CijR3(DCMPLX(C0123),C123R_d,C123I_d)
cfc       Endif
cfc       If(PrintC) then
cfc       write(*,*)"C124",C0124R_d,C0124I_d
cfc       call print_CijR3(DCMPLX(C0124),C124R_d,C124I_d)
cfc       Endif
cfc
cfc       If(PrintC) then
cfc       write(*,*)"C234",C0234R_d,C0234I_d
cfc       call print_CijR3(DCMPLX(C0234),C234R_d,C234I_d)
cfc       Endif
cfc
cfc
cfc       If(PrintC) then
cfc       write(*,*)"C134",C0134R_d,C0134I_d
cfc       call print_CijR3(DCMPLX(C0134),C134R_d,C134I_d)
cfc       Endif
cfc
cfc
cfc       If(PrintD) then
cfc       write(*,*)"D1234",D01234R_d,D01234I_d
cfc       write(*,*)"D1234",D01234R,D01234I
cfc       call print_DijR4(DCMPLX(D01234),D1234R_d,D1234I_d)
cfc       endif


       return
       end


      subroutine R3_From_BijR_QUAD_to_DblR2(B0R,B0I,BijR,BijI,
     &    B0R_QUAD,B0I_QUAD,BijR_QUAD,BijI_QUAD)
      implicit none
      real*16 B0R_QUAD,B0I_QUAD,BijR_QUAD(3,3),BijI_QUAD(3,3)
      real*8 B0R,B0I,BijR(2,2),BijI(2,2)
      integer indx_B(3),i1,i2
cccccccccccccccccccccccccccccccccccccccc
      indx_B(1)=1
      indx_B(2)=2
      indx_B(3)=3
ccccccccccccccccccccccccccccccccccccc
      B0R=B0R_QUAD
      B0I=B0I_QUAD
      do i1=1,2
         do i2=1,indx_B(i1)
           BijR(i2,i1)=BijR_QUAD(i2,i1)
           BijI(i2,i1)=BijI_QUAD(i2,i1)
         enddo
      enddo
      end


      subroutine R4_From_CijR_QUAD_to_DblR3(C0R,C0I,CijR,CijI,
     -              C0R_QUAD,C0I_QUAD,CijR_QUAD,CijI_QUAD)
      implicit none
      real*16 C0R_QUAD,C0I_QUAD
      real*8 C0R,C0I
      real*16  CijR_QUAD(9,4),CijI_QUAD(9,4)
      real*8  CijR(6,3),CijI(6,3)
      integer indx_C(4)
      integer i1,i2
      indx_C(1)=2
      indx_C(2)=4
      indx_C(3)=6
      indx_C(4)=9
      C0R=C0R_QUAD
      C0I=C0I_QUAD
      do i1=1,3
         do i2=1,indx_C(i1)
           CijR(i2,i1)=CijR_QUAD(i2,i1)
           CijI(i2,i1)=CijI_QUAD(i2,i1)
         enddo
      enddo
      end


      subroutine R5_From_DijR_QUAD_to_DblR4(D0R,D0I,DijR,DijI,
     -            D0R_QUAD, D0I_QUAD,DijR_QUAD,DijI_QUAD)
      implicit none
      real*16 D0R_QUAD,D0I_QUAD
      real*8 D0R,D0I
      real*16  DijR_QUAD(34,5),DijI_QUAD(34,5)
      real*8  DijR(22,4),DijI(22,4)
      integer indx_D(5)
      integer i1,i2
      indx_D(1)=13
      indx_D(2)=7
      indx_D(3)=13
      indx_D(4)=22
      indx_D(5)=34


      D0R=D0R_QUAD
      D0I=D0I_QUAD
      do i1=1,4
         do i2=1,indx_D(i1)
           DijR(i2,i1)=DijR_QUAD(i2,i1)
           DijI(i2,i1)=DijI_QUAD(i2,i1)
         enddo
      enddo
      end


      subroutine print_CijR3(C0,CijR,CijI)
      implicit none
      complex*16 C0
      real*8  CijR(6,3),CijI(6,3)
      integer indx_C(3)
c      common/indexes/indx_C,indx_D,indx_E
      integer i1,i2
      indx_C(1)=2
      indx_C(2)=4
      indx_C(3)=6
      write(*,'(A4,E25.16,E25.16,A1)') 'C0=(',C0,')'
      do i1=1,3
         do i2=1,indx_C(i1)
            write(*,'(A4,I1,A1,I1,A3,E25.16,E25.16,A1)') 
     -   'Cij(',i2,',',i1,')=(',DCMPLX(CijR(i2,i1),CijI(i2,i1)),')'
         enddo
      enddo
      end


      subroutine print_DijR4(D0,DijR,DijI)
      implicit none
      complex*16 D0
      real*8  DijR(22,4),DijI(22,4)
      integer indx_D(5)
c      common/indexes/indx_C,indx_D,indx_E
      integer i1,i2
      indx_D(1)=13
      indx_D(2)=7
      indx_D(3)=13
      indx_D(4)=22
      indx_D(5)=34
      write(*,'(A4,E25.16,E25.16,A1)') 'D0=(',D0,')'
      do i1=1,4
         Print*, "********rank",i1
         do i2=1,indx_D(i1)
            write(*,'(A4,I2,A1,I1,A3,E25.16,E25.16,A1)')
     -     'Dij(',i2,',',i1,')=(',DCMPLX(DijR(i2,i1),DijI(i2,i1)),')'
         enddo
         enddo
      end
