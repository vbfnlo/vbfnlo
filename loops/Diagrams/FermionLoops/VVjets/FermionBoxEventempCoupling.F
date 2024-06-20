       subroutine FermionBoxEvenCouplings(za1,za2,za3,za4,zb1,zb2,zb3,
     -   zb4,m,p1,p2,p3,p4,mup1,mup2,mup3,mup4,musq,comp,result,resul
     -   tn)
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
       Complex*16 A0fin,B0t1,C0fin1M,D0finG
       EXTERNAL dotrr,A0fin,B0t1,C0fin1M,D0finG,B0fin,C0fin,D0fin
       Complex*16 A0,B012,B013,B014,B023,B024,B034
      Real*8 A0R,A0I,B012R,B013R,B014R,B023R,B024R,B034R
      Real*8 B012I,B013I,B014I,B023I,B024I,B034I
      Real*8 B12R(2,2),B13R(2,2),B14R(2,2),B23R(2,2),B24R(2,2),B34R(2,2)
      Real*8 B12I(2,2),B13I(2,2),B14I(2,2),B23I(2,2),B24I(2,2),B34I(2,2)
       Complex*16 C0123,C0124,C0134,C0234
       Real*8 C0123R,C0124R,C0134R,C0234R
       Real*8 C0123I,C0124I,C0134I,C0234I
       Real*8 C123R(6,3),C124R(6,3),C134R(6,3),C234R(6,3)
       Real*8 C123I(6,3),C124I(6,3),C134I(6,3),C234I(6,3)
       Complex*16 D01234
        Real*8 D01234R
       Real*8 D01234I
       Real*8 D1234R(22,4)
        Real*8 D1234I(22,4)
      Complex*16 F(279)
      Real*8 FI(279),FR(279),P(452)
      Complex*16 mup1(0:3),mup2(0:3),mup3(0:3),mup4(0:3)
      Complex*16 dotrc,dotcc,result,resultn(5)
       Real*8 musq,m,msq,za1,za2,za3,za4,zb1,zb2,zb3,zb4
      EXTERNAL   dotrc,dotcc
      Integer comp
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
      Common/musqInv/musqcp
      real*8 musqcp
       logical printmy
      common/mprint/printmy
      logical printB,printC,printD
      musqcp=musq
      printB=.true.
      printC=.true.
      printD=.true.
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
       msq=m*m
       p1sq = dotrr(p1,p1)
       p1p2 = dotrr(p1,p2)
       p1p3 = dotrr(p1,p3)
       p1p4 = dotrr(p1,p4)
       p2sq = dotrr(p2,p2)
       p2p3 = dotrr(p2,p3)
       p2p4 = dotrr(p2,p4)
       p3sq = dotrr(p3,p3)
       p3p4 = dotrr(p3,p4)
       p4sq = dotrr(p4,p4)
       s12 = (p1sq +p2sq+ 2*p1p2) 
       s13 = (p1sq +p3sq+ 2*p1p3) 
       s14 = (p1sq +p4sq+ 2*p1p4) 
       s23 = (p2sq +p3sq+ 2*p2p3) 
       s24 = (p2sq +p4sq+ 2*p2p4) 
       s34 = (p3sq +p4sq+ 2*p3p4) 
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
        if(m.gt.1d-6) then
c    Calling C_ij,D_ij Functions    
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
        A0=A0fin(m,musq)  
        A0R=Dble(A0)   
        A0I=DImag(A0)   
       B012=B0t1(m,p1sq,musq)  
       B013=B0t1(m,s12,musq)   
       B014=B0t1(m,p4sq,musq)   
       B023=B0t1(m,p2sq,musq)  
       B024=B0t1(m,s23,musq)   
       B034=B0t1(m,p3sq,musq)  
c************************************************************************************
       call tens_red2_Re_1M(m,p1sq,A0R,A0I,B012,B012R,B012I,B12R,B12I)
     -    
       call tens_red2_Re_1M(m,s12,A0R,A0I,B013,B013R,B013I,B13R,B13I) 
     -   
       call tens_red2_Re_1M(m,p4sq,A0R,A0I,B014,B014R,B014I,B14R,B14I)
     -    
       call tens_red2_Re_1M(m,p2sq,A0R,A0I,B023,B023R,B023I,B23R,B23I)
     -    
       call tens_red2_Re_1M(m,s23,A0R,A0I,B024,B024R,B024I,B24R,B24I) 
     -   
       call tens_red2_Re_1M(m,p3sq,A0R,A0I,B034,B034R,B034I,B34R,B34I)
     -    
c************************************************************************************
      C0123=C0fin1M(m,p1sq,p2sq,s12) 
      C0124=C0fin1M(m,p1sq,s23,p4sq)  
      C0134=C0fin1M(m,s12,p3sq,p4sq)  
      C0234=C0fin1M(m,p2sq,p3sq,s23) 
c************************************************************************************
c************************************************************************************
c************************************************************************************
       call tens_red3_new_Re_Com_1MRank3(m,p1sq,p2sq,s12,B023R,B013R,B
     -   012R,B23R,B13R,B12R,B023I,B013I,B012I,B23I,B13I,B12I,C0123,C
     -   0123R,C0123I,C123R,C123I) 
       call tens_red3_new_Re_Com_1MRank3(m,p1sq,s23,p4sq,B024R,B014R,B
     -   012R,B24R,B14R,B12R,B024I,B014I,B012I,B24I,B14I,B12I,C0124,C
     -   0124R,C0124I,C124R,C124I) 
       call tens_red3_new_Re_Com_1MRank3(m,s12,p3sq,p4sq,B034R,B014R,B
     -   013R,B34R,B14R,B13R,B034I,B014I,B013I,B34I,B14I,B13I,C0134,C
     -   0134R,C0134I,C134R,C134I) 
       call tens_red3_new_Re_Com_1MRank3(m,p2sq,p3sq,s23,B034R,B024R,B
     -   023R,B34R,B24R,B23R,B034I,B024I,B023I,B34I,B24I,B23I,C0234,C
     -   0234R,C0234I,C234R,C234I) 
c************************************************************************************
c************************************************************************************
cMK call D0finG because D04 is buggy !!!
       D01234=D0finG(m,m,m,m,s12,s23,p1sq,p2sq,p3sq,p4sq,musq)
c************************************************************************************
c************************************************************************************
       call tens_red4_new_Re_Com_1MRank4(m,p1sq,p2sq,p3sq,p1p2,p1p3,p2
     -   p3,C0234R,C0134R,C0124R,C0123R,C234R,C134R,C124R,C123R,C0234
     -   I,C0134I,C0124I,C0123I,C234I,C134I,C124I,C123I,D01234,D01234
     -   R,D01234I,D1234R,D1234I)
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
      elseif(m.lt.1d-6) then
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c    Calling C_ij,D_ij Functions    
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       B012=B0fin(p1sq,musq)   
       B013=B0fin(s12,musq)   
       B014=B0fin(p4sq,musq)   
       B023=B0fin(p2sq,musq)  
       B024=B0fin(s23,musq)   
       B034=B0fin(p3sq,musq)  
       call tens_red2_Re(p1sq,B012,B012R,B012I,B12R,B12I) 
       call tens_red2_Re(s12,B013,B013R,B013I,B13R,B13I) 
       call tens_red2_Re(p4sq,B014,B014R,B014I,B14R,B14I) 
       call tens_red2_Re(p2sq,B023,B023R,B023I,B23R,B23I) 
       call tens_red2_Re(s23,B024,B024R,B024I,B24R,B24I) 
       call tens_red2_Re(p3sq,B034,B034R,B034I,B34R,B34I) 
c************************************************************************************
c************************************************************************************
c************************************************************************************
      C0123=C0fin(p1sq,p2sq,s12,musq) 
      C0124=C0fin(p1sq,s23,p4sq,musq)  
      C0134=C0fin(s12,p3sq,p4sq,musq)  
      C0234=C0fin(p2sq,p3sq,s23,musq) 
c************************************************************************************
c************************************************************************************
       call tens_red3_new_Re_ComRank3(p1sq,p2sq,s12,B023R,B013R,B012R,
     -   B23R,B13R,B12R,B023I,B013I,B012I,B23I,B13I,B12I,C0123,C0123R
     -   ,C0123I,C123R,C123I) 
       call tens_red3_new_Re_ComRank3(p1sq,s23,p4sq,B024R,B014R,B012R,
     -   B24R,B14R,B12R,B024I,B014I,B012I,B24I,B14I,B12I,C0124,C0124R
     -   ,C0124I,C124R,C124I) 
       call tens_red3_new_Re_ComRank3(s12,p3sq,p4sq,B034R,B014R,B013R,
     -   B34R,B14R,B13R,B034I,B014I,B013I,B34I,B14I,B13I,C0134,C0134R
     -   ,C0134I,C134R,C134I) 
       call tens_red3_new_Re_ComRank3(p2sq,p3sq,s23,B034R,B024R,B023R,
     -   B34R,B24R,B23R,B034I,B024I,B023I,B34I,B24I,B23I,C0234,C0234R
     -   ,C0234I,C234R,C234I) 
c************************************************************************************
c************************************************************************************
         D01234=D0fin(s12,s23,p1sq,p2sq,p3sq,p4sq,musq)
c       Print*," D0",D0
c************************************************************************************
c************************************************************************************
       call tens_red4_new_Re_ComRank4(p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,C0
     -   234R,C0134R,C0124R,C0123R,C234R,C134R,C124R,C123R,C0234I,C01
     -   34I,C0124I,C0123I,C234I,C134I,C124I,C123I,D01234,D01234R,D01
     -   234I,D1234R,D1234I)
         endif
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       endif  


cfc
cfc       If(PrintB) then
cfc       write(*,*)"B12"
cfc       write(*,*) "B12=", B012
cfc       endif
cfc       If(PrintB) then
cfc       write(*,*)"B23"
cfc       write(*,*) "B23=", B023
cfc       endif
cfc       If(PrintB) then
cfc       write(*,*)"B34"
cfc       write(*,*) "B34=", B034
cfc       endif
cfc
cfc       If(PrintB) then
cfc       write(*,*)"B13"
cfc       write(*,*) "B13=", B013
cfc       endif
cfc       If(PrintB) then
cfc       write(*,*)"B14"
cfc       write(*,*) "B14=", B014
cfc       endif
cfc       If(PrintB) then
cfc       write(*,*)"B24"
cfc       write(*,*) "B24=", B024
cfc       endif
cfc
cfc       If(PrintC) then
cfc       write(*,*)"C123",C0123R,C0123I
cfc       
cfc       call print_CijR3(DCMPLX(C0123),C123R,C123I)
cfc       Endif
cfc       If(PrintC) then
cfc       write(*,*)"C124",C0124R,C0124I
cfc       call print_CijR3(DCMPLX(C0124),C124R,C124I)
cfc       Endif
cfc
cfc       If(PrintC) then
cfc       write(*,*)"C234",C0234R,C0234I
cfc       call print_CijR3(DCMPLX(C0234),C234R,C234I)
cfc       Endif
cfc
cfc
cfc       If(PrintC) then
cfc       write(*,*)"C134",C0134R,C0134I
cfc       call print_CijR3(DCMPLX(C0134),C134R,C134I)
cfc       Endif
cfc
cfc
cfc       If(PrintD) then
cfc       write(*,*)"D1234",D01234R,D01234I
cfc       write(*,*)"D1234",D01234R,D01234I
cfc       call print_DijR4(DCMPLX(D01234),D1234R,D1234I)
cfc       endif


c       Definition of the F,P functions:Independent of the currents    
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       P(1) = za2-zb2
       P(2) = za3-zb3
       P(3) = za1*za4+zb1*zb4
       P(4) = za1*za2*za4-zb1*zb2*zb4
       P(5) = 3*za2-zb2
       P(6) = za2-3*zb2
       P(7) = za1*za4*P(5)+zb1*zb4*P(6)
       P(8) = za1*za2*za3*za4+zb1*zb2*zb3*zb4
       P(9) = za1-zb1
       P(10) = za4*zb2+za2*zb4
       P(11) = za2+2*zb2
       P(12) = za4*zb2+zb4*P(11)
       P(13) = 2*za4+zb4
       P(14) = za4*zb2+za2*P(13)
       P(15) = -(zb1*P(12))+za1*P(14)
       P(16) = za1*za4*zb2-za2*zb1*zb4
       P(17) = za4*zb1*zb2-za1*za2*zb4
       P(18) = za1*za2+zb1*zb2
       P(19) = za4-zb4
       P(20) = za4-3*zb4
       P(21) = 3*za4-zb4
       P(22) = zb1*zb2*P(20)+za1*za2*P(21)
       P(23) = -(p1sq*za1*za2*za4)-p4sq*za1*za2*za4+s12*za1*za2*za4+p4
     -   sq*za1*za4*zb2+p3sq*za4*zb1*zb2-s23*za4*zb1*zb2-p3sq*za1*za2
     -   *zb4+s23*za1*za2*zb4-p4sq*za2*zb1*zb4+p1sq*zb1*zb2*zb4+p4sq*
     -   zb1*zb2*zb4-s12*zb1*zb2*zb4+2*M**2*P(1)*P(9)*P(19)
       P(24) = za1*za2-zb1*zb2
       P(25) = za4+zb4
       P(26) = -(s12*za1*za2*za4)-s23*za1*za2*za4+p3sq*za1*za4*zb2+s12
     -   *za1*za4*zb2-p3sq*za4*zb1*zb2+s23*za4*zb1*zb2+p3sq*za1*za2*z
     -   b4-s23*za1*za2*zb4-p3sq*za2*zb1*zb4-s12*za2*zb1*zb4+s12*zb1*
     -   zb2*zb4+s23*zb1*zb2*zb4+p4sq*P(1)*P(3)+p2sq*P(24)*P(25)
       P(27) = p1sq-s12
       P(28) = za4-2*zb4
       P(29) = -(za3*za4)+zb3*P(28)
       P(30) = 2*za3*za4-za3*zb4+zb3*zb4
       P(31) = za1*P(29)+zb1*P(30)
       P(32) = P(8)*P(27)+M**2*P(1)*P(31)
       P(33) = za3+zb3
       P(34) = p4sq*P(8)-M**2*P(19)*P(24)*P(33)
       P(35) = za3*za4+zb3*zb4
       P(36) = p1sq+p2sq-s12
       P(37) = -2*M**2*P(1)*P(9)*P(35)+P(8)*P(36)
       P(38) = p1sq-p4sq+s23
       P(39) = 2*M**2*P(1)*P(9)*P(35)-P(8)*P(38)
       P(40) = -2*M**2*P(1)*P(9)*P(35)+P(8)*P(38)
       P(41) = -(za1*za2*za3*za4)-zb1*zb2*zb3*zb4
       P(42) = s12*za1*za2*za4-s23*za1*za2*za4+p3sq*za1*za4*zb2-s12*za
     -   1*za4*zb2+p3sq*za4*zb1*zb2-s23*za4*zb1*zb2-p3sq*za1*za2*zb4+
     -   s23*za1*za2*zb4-p3sq*za2*zb1*zb4+s12*za2*zb1*zb4-s12*zb1*zb2
     -   *zb4+s23*zb1*zb2*zb4-p4sq*P(1)*P(3)+p2sq*P(18)*P(19)
       P(43) = za3*za4-za4*zb3+2*zb3*zb4
       P(44) = -2*za4+zb4
       P(45) = -(zb3*zb4)+za3*P(44)
       P(46) = za1*P(43)+zb1*P(45)
       P(47) = -(P(8)*P(27))+M**2*P(1)*P(46)
       P(48) = -(p4sq*P(8))+M**2*P(19)*P(24)*P(33)
       P(49) = za3*zb1+za1*zb3
       P(50) = M**2*P(1)*P(19)*P(49)
       P(51) = p2sq*za1*za2*za3*za4+s12*za1*za2*za3*za4+p2sq*zb1*zb2*z
     -   b3*zb4+s12*zb1*zb2*zb3*zb4-p1sq*P(8)-2*P(50)
       P(52) = p1sq*za1*za2*za3*za4-p4sq*za1*za2*za3*za4-s23*za1*za2*z
     -   a3*za4+p1sq*zb1*zb2*zb3*zb4-p4sq*zb1*zb2*zb3*zb4-s23*zb1*zb2
     -   *zb3*zb4+2*P(50)
       P(53) = p4sq*za1*za2*za3*za4+s23*za1*za2*za3*za4+p4sq*zb1*zb2*z
     -   b3*zb4+s23*zb1*zb2*zb3*zb4-p1sq*P(8)-2*P(50)
       P(54) = -(p1sq*za1*za2*za4)+p4sq*za1*za4*zb2-s23*za4*zb1*zb2+s2
     -   3*za1*za2*zb4-p4sq*za2*zb1*zb4+p1sq*zb1*zb2*zb4+2*M**2*P(1)*
     -   P(9)*P(19)
       P(55) = zb2*P(43)+za2*P(45)
       P(56) = p1sq*P(8)+M**2*P(9)*P(55)
       P(57) = p2sq-s23
       P(58) = -(M**2*P(1)*P(2)*P(3))+P(8)*P(57)
       P(59) = p2sq+p3sq-s23
       P(60) = za2+zb2
       P(61) = za1*za4-zb1*zb4
       P(62) = -(s12*za1*za2*za4)-s23*za1*za2*za4+p3sq*za1*za4*zb2-s12
     -   *za1*za4*zb2-p3sq*za4*zb1*zb2-s23*za4*zb1*zb2+p3sq*za1*za2*z
     -   b4+s23*za1*za2*zb4-p3sq*za2*zb1*zb4+s12*za2*zb1*zb4+s12*zb1*
     -   zb2*zb4+s23*zb1*zb2*zb4+p2sq*P(18)*P(19)+p4sq*P(60)*P(61)
       P(63) = p3sq+3*P(57)
       P(64) = za2*za3+zb2*zb3
       P(65) = p1sq-p2sq+s12
       P(66) = -2*M**2*P(9)*P(19)*P(64)+P(8)*P(65)
       P(67) = p2sq+3*p3sq-s23
       P(68) = zb2*P(29)+za2*P(30)
       P(69) = -(p1sq*P(8))+M**2*P(9)*P(68)
       P(70) = M**2*P(1)*P(2)*P(3)-P(8)*P(57)
       P(71) = za1*za4*zb2-za2*zb1*zb4-2*P(17)
       P(72) = za4*zb2+zb4*P(60)
       P(73) = za4*zb2+za2*P(25)
       P(74) = -(zb1*P(72))+za1*P(73)
       P(75) = 2*za4-zb4
       P(76) = za1*za2*P(28)+zb1*zb2*P(75)
       P(77) = zb1*zb2*P(28)+za1*za2*P(75)
       P(78) = -(za2*za4)+za4*zb2+2*za2*zb4
       P(79) = -za2+zb2
       P(80) = -2*za4*zb2+zb4*P(79)
       P(81) = za1*P(78)+zb1*P(80)
       P(82) = 2*za4*zb2+zb4*P(60)
       P(83) = za4+2*zb4
       P(84) = za4*zb2+za2*P(83)
       P(85) = -(zb1*P(82))+za1*P(84)
       P(86) = za1+zb1
       P(87) = za4*zb2-za2*zb4
       P(88) = 2*za2*za4-za4*zb2-za2*zb4
       P(89) = za2-2*zb2
       P(90) = za4*zb2+zb4*P(89)
       P(91) = za1*P(88)+zb1*P(90)
       P(92) = -(za4*zb2)+za2*P(28)
       P(93) = 2*za4*zb2+zb4*P(1)
       P(94) = za1*P(92)+zb1*P(93)
       P(95) = za4*zb2+za2*P(28)
       P(96) = 2*za4*zb2-zb4*P(60)
       P(97) = za1*P(95)+zb1*P(96)
       P(98) = -(p1sq*za1*za2*za4)+s12*za1*za4*zb2-p2sq*za4*zb1*zb2+p2
     -   sq*za1*za2*zb4-s12*za2*zb1*zb4+p1sq*zb1*zb2*zb4+2*M**2*P(1)*
     -   P(9)*P(19)
       P(99) = -(p1sq*za1*za2*za4)+p4sq*za1*za2*za4-s12*za1*za2*za4+p4
     -   sq*za1*za4*zb2-p3sq*za4*zb1*zb2-s23*za4*zb1*zb2+p3sq*za1*za2
     -   *zb4+s23*za1*za2*zb4-p4sq*za2*zb1*zb4+p1sq*zb1*zb2*zb4-p4sq*
     -   zb1*zb2*zb4+s12*zb1*zb2*zb4+2*M**2*P(1)*P(9)*P(19)
       P(100) = p1sq+s12
       P(101) = za4*P(64)
       P(102) = za2*za3*zb4+za3*zb2*zb4-za2*zb3*zb4+3*zb2*zb3*zb4-2*P(
     -   101)
       P(103) = zb4*P(64)
       P(104) = 3*za2*za3*za4-za3*za4*zb2+za2*za4*zb3+za4*zb2*zb3-2*P(
     -   103)
       P(105) = zb1*P(102)+za1*P(104)
       P(106) = P(8)*P(100)-M**2*P(105)
       P(107) = p1sq+p4sq-s23
       P(108) = 2*M**2*P(9)*P(19)*P(64)-P(8)*P(107)
       P(109) = -2*M**2*P(9)*P(19)*P(64)+P(8)*P(107)
       P(110) = za4*zb2+zb4*P(1)
       P(111) = -za4+zb4
       P(112) = za4*zb2+za2*P(111)
       P(113) = zb1*P(110)-za1*P(112)
       P(114) = za1*za2*za3-zb1*zb2*zb3
       P(115) = -za3+zb3
       P(116) = -2*zb2*zb3+za2*P(115)
       P(117) = 2*za2*za3+zb2*P(115)
       P(118) = za1*P(116)+zb1*P(117)
       P(119) = za2*za3*zb1-za1*zb2*zb3
       P(120) = za3-2*zb3
       P(121) = za2*za3+zb2*P(120)
       P(122) = 2*za2*za3-za2*zb3-zb2*zb3
       P(123) = zb1*P(121)+za1*P(122)
       P(124) = za3-3*zb3
       P(125) = 3*za3-zb3
       P(126) = zb1*zb2*P(124)+za1*za2*P(125)
       P(127) = za2*za3-zb2*zb3
       P(128) = p1sq*za1*za2*za3-s23*za1*za2*za3+p4sq*za2*za3*zb1-p4sq
     -   *za1*zb2*zb3-p1sq*zb1*zb2*zb3+s23*zb1*zb2*zb3
       P(129) = -(s12*za1*za2*za3)-s23*za1*za2*za3+p4sq*za2*za3*zb1-s1
     -   2*za2*za3*zb1-p4sq*za3*zb1*zb2-s23*za3*zb1*zb2+p4sq*za1*za2*
     -   zb3+s23*za1*za2*zb3-p4sq*za1*zb2*zb3+s12*za1*zb2*zb3+s12*zb1
     -   *zb2*zb3+s23*zb1*zb2*zb3+p1sq*P(2)*P(18)+p3sq*P(86)*P(127)
       P(130) = p1sq*P(8)-M**2*P(9)*P(25)*P(127)
       P(131) = -(s23*P(8))+M**2*P(9)*P(25)*P(127)
       P(132) = -2*M**2*P(1)*P(2)*P(3)+P(8)*P(59)
       P(133) = p1sq-p2sq-s12
       P(134) = 2*M**2*P(1)*P(2)*P(3)+P(8)*P(133)
       P(135) = 2*M**2*P(1)*P(2)*P(3)-P(8)*P(59)
       P(136) = s12*za1*za2*za3-s23*za1*za2*za3+p4sq*za2*za3*zb1-s12*z
     -   a2*za3*zb1+p4sq*za3*zb1*zb2-s23*za3*zb1*zb2-p4sq*za1*za2*zb3
     -   +s23*za1*za2*zb3-p4sq*za1*zb2*zb3+s12*za1*zb2*zb3-s12*zb1*zb
     -   2*zb3+s23*zb1*zb2*zb3+p1sq*P(2)*P(18)-p3sq*P(9)*P(64)
       P(137) = -(p1sq*P(8))+M**2*P(9)*P(25)*P(127)
       P(138) = s23*P(8)-M**2*P(9)*P(25)*P(127)
       P(139) = M**2*P(2)*P(9)*P(10)
       P(140) = p3sq*za1*za2*za3*za4+s23*za1*za2*za3*za4+p3sq*zb1*zb2*
     -   zb3*zb4+s23*zb1*zb2*zb3*zb4-p2sq*P(8)-2*P(139)
       P(141) = -(p1sq*za1*za2*za3*za4)+p2sq*za1*za2*za3*za4-s12*za1*z
     -   a2*za3*za4-p1sq*zb1*zb2*zb3*zb4+p2sq*zb1*zb2*zb3*zb4-s12*zb1
     -   *zb2*zb3*zb4+2*P(139)
       P(142) = p2sq*za1*za2*za3*za4-p3sq*za1*za2*za3*za4-s23*za1*za2*
     -   za3*za4+p2sq*zb1*zb2*zb3*zb4-p3sq*zb1*zb2*zb3*zb4-s23*zb1*zb
     -   2*zb3*zb4+2*P(139)
       P(143) = -(s12*za1*za2*za3)-s23*za1*za2*za3+p4sq*za2*za3*zb1+s1
     -   2*za2*za3*zb1-p4sq*za3*zb1*zb2+s23*za3*zb1*zb2+p4sq*za1*za2*
     -   zb3-s23*za1*za2*zb3-p4sq*za1*zb2*zb3-s12*za1*zb2*zb3+s12*zb1
     -   *zb2*zb3+s23*zb1*zb2*zb3+p1sq*P(24)*P(33)+p3sq*P(9)*P(64)
       P(144) = p2sq-p3sq+s23
       P(145) = -2*M**2*P(1)*P(9)*P(35)+P(8)*P(144)
       P(146) = 2*M**2*P(1)*P(9)*P(35)-P(8)*P(36)
       P(147) = 2*M**2*P(1)*P(9)*P(35)-P(8)*P(144)
       P(148) = p1sq*za1*za2*za3-p3sq*za1*za2*za3+s12*za1*za2*za3-s23*
     -   za1*za2*za3+p4sq*za2*za3*zb1+p4sq*za3*zb1*zb2-p4sq*za1*za2*z
     -   b3-p4sq*za1*zb2*zb3-p1sq*zb1*zb2*zb3+p3sq*zb1*zb2*zb3-s12*zb
     -   1*zb2*zb3+s23*zb1*zb2*zb3
       P(149) = za1*za2*za3+zb1*zb2*zb3
       P(150) = -(za2*za3*zb1)-za3*zb1*zb2-za1*za2*zb3-za1*zb2*zb3+2*P
     -   (149)
       P(151) = -(P(8)*P(100))+M**2*P(25)*P(150)
       P(152) = -2*za3*za4+zb3*P(19)
       P(153) = 2*zb3*zb4+za3*P(19)
       P(154) = za1*P(152)+zb1*P(153)
       P(155) = p2sq*P(8)+M**2*P(1)*P(154)
       P(156) = p1sq*za1*za2*za3+p3sq*za1*za2*za3-s12*za1*za2*za3-s23*
     -   za1*za2*za3+p4sq*za2*za3*zb1-p4sq*za3*zb1*zb2+p4sq*za1*za2*z
     -   b3-p4sq*za1*zb2*zb3-p1sq*zb1*zb2*zb3-p3sq*zb1*zb2*zb3+s12*zb
     -   1*zb2*zb3+s23*zb1*zb2*zb3
       P(157) = za3*zb1-za1*zb3
       P(158) = -(p1sq*P(8))+s12*P(8)-M**2*P(1)*P(25)*P(157)
       P(159) = -2*zb3*zb4+za3*P(111)
       P(160) = 2*za3*za4+zb3*P(111)
       P(161) = zb1*P(159)+za1*P(160)
       P(162) = -(p2sq*P(8))+M**2*P(1)*P(161)
       P(163) = p1sq*P(8)-s12*P(8)+M**2*P(1)*P(25)*P(157)
       P(164) = -(p2sq*za1*za2*za3)+s12*za2*za3*zb1-p1sq*za3*zb1*zb2+p
     -   1sq*za1*za2*zb3-s12*za1*zb2*zb3+p2sq*zb1*zb2*zb3+2*M**2*P(1)
     -   *P(2)*P(9)
       P(165) = -2*zb2*zb3+za2*P(33)
       P(166) = 2*za2*za3-zb2*P(33)
       P(167) = za1*P(165)+zb1*P(166)
       P(168) = za2*za3-za3*zb2+2*zb2*zb3
       P(169) = -2*za3+zb3
       P(170) = -(zb2*zb3)+za2*P(169)
       P(171) = zb1*P(168)+za1*P(170)
       P(172) = -2*za2*za3+zb2*P(124)
       P(173) = 3*za2*za3-za2*zb3+2*zb2*zb3
       P(174) = zb1*P(172)+za1*P(173)
       P(175) = 2*za3+zb3
       P(176) = zb2*zb3+za2*P(175)
       P(177) = za3+2*zb3
       P(178) = za2*za3+zb2*P(177)
       P(179) = za1*P(176)-zb1*P(178)
       P(180) = 3*za3+zb3
       P(181) = 2*zb2*zb3+za2*P(180)
       P(182) = za3+3*zb3
       P(183) = 2*za2*za3+zb2*P(182)
       P(184) = za1*P(181)-zb1*P(183)
       P(185) = -(za2*za3)+zb2*P(120)
       P(186) = 2*za2*za3-za2*zb3+zb2*zb3
       P(187) = zb1*P(185)+za1*P(186)
       P(188) = 2*za3-zb3
       P(189) = zb1*zb2*P(120)+za1*za2*P(188)
       P(190) = -(za2*za3*zb1)+za1*zb2*zb3+2*P(114)
       P(191) = p4sq*za2*za3*za4-s23*za2*za3*za4-p1sq*za4*zb2*zb3+p1sq
     -   *za2*za3*zb4-p4sq*zb2*zb3*zb4+s23*zb2*zb3*zb4
       P(192) = -(s12*za2*za3*za4)+s23*za2*za3*za4+p1sq*za3*za4*zb2-s2
     -   3*za3*za4*zb2+p1sq*za4*zb2*zb3-s12*za4*zb2*zb3-p1sq*za2*za3*
     -   zb4+s12*za2*za3*zb4-p1sq*za2*zb3*zb4+s23*za2*zb3*zb4+s12*zb2
     -   *zb3*zb4-s23*zb2*zb3*zb4-p4sq*P(1)*P(35)+p2sq*P(19)*P(64)
       P(193) = p4sq*P(8)-M**2*P(19)*P(86)*P(127)
       P(194) = -(s23*P(8))+M**2*P(19)*P(86)*P(127)
       P(195) = p4sq*za1*za2*za3*za4+s12*za1*za2*za3*za4+p4sq*zb1*zb2*
     -   zb3*zb4+s12*zb1*zb2*zb3*zb4-p3sq*P(8)-2*P(50)
       P(196) = -(p2sq*za1*za2*za3*za4)+p3sq*za1*za2*za3*za4-s23*za1*z
     -   a2*za3*za4-p2sq*zb1*zb2*zb3*zb4+p3sq*zb1*zb2*zb3*zb4-s23*zb1
     -   *zb2*zb3*zb4+2*P(50)
       P(197) = p3sq*za1*za2*za3*za4-p4sq*za1*za2*za3*za4-s12*za1*za2*
     -   za3*za4+p3sq*zb1*zb2*zb3*zb4-p4sq*zb1*zb2*zb3*zb4-s12*zb1*zb
     -   2*zb3*zb4+2*P(50)
       P(198) = s23*za1*za2*za3*za4+s23*zb1*zb2*zb3*zb4+p2sq*P(8)-p3sq
     -   *P(8)-2*P(50)
       P(199) = -(s12*za2*za3*za4)-s23*za2*za3*za4+p1sq*za3*za4*zb2+s2
     -   3*za3*za4*zb2-p1sq*za4*zb2*zb3+s12*za4*zb2*zb3+p1sq*za2*za3*
     -   zb4-s12*za2*za3*zb4-p1sq*za2*zb3*zb4-s23*za2*zb3*zb4+s12*zb2
     -   *zb3*zb4+s23*zb2*zb3*zb4+p4sq*P(1)*P(35)+p2sq*P(25)*P(127)
       P(200) = -(p4sq*P(8))+M**2*P(19)*P(86)*P(127)
       P(201) = s23*P(8)-M**2*P(19)*P(86)*P(127)
       P(202) = p3sq-p4sq+s12
       P(203) = -2*M**2*P(1)*P(2)*P(3)+P(8)*P(202)
       P(204) = 2*M**2*P(1)*P(2)*P(3)-P(8)*P(202)
       P(205) = za4*zb2*zb3-za2*za3*zb4
       P(206) = za2*za3*za4-zb2*zb3*zb4
       P(207) = za3*za4-zb3*zb4
       P(208) = -(s12*za2*za3*za4)-s23*za2*za3*za4+p1sq*za3*za4*zb2-s2
     -   3*za3*za4*zb2-p1sq*za4*zb2*zb3-s12*za4*zb2*zb3+p1sq*za2*za3*
     -   zb4+s12*za2*za3*zb4-p1sq*za2*zb3*zb4+s23*za2*zb3*zb4+s12*zb2
     -   *zb3*zb4+s23*zb2*zb3*zb4+p2sq*P(19)*P(64)+p4sq*P(60)*P(207)
       P(209) = p3sq+p4sq-s12
       P(210) = -2*M**2*P(2)*P(18)*P(19)+P(8)*P(209)
       P(211) = p2sq-p3sq-s23
       P(212) = 2*M**2*P(2)*P(18)*P(19)+P(8)*P(211)
       P(213) = 2*M**2*P(2)*P(18)*P(19)-P(8)*P(209)
       P(214) = -2*M**2*P(2)*P(18)*P(19)-P(8)*P(211)
       P(215) = -(p3sq*za2*za3*za4)-p4sq*za2*za3*za4+s23*za2*za3*za4+p
     -   4sq*za3*za4*zb2+p1sq*za4*zb2*zb3-s12*za4*zb2*zb3-p1sq*za2*za
     -   3*zb4+s12*za2*za3*zb4-p4sq*za2*zb3*zb4+p3sq*zb2*zb3*zb4+p4sq
     -   *zb2*zb3*zb4-s23*zb2*zb3*zb4+2*M**2*P(1)*P(2)*P(19)
       P(216) = p2sq*P(8)-M**2*P(1)*P(86)*P(207)
       P(217) = -(p3sq*za2*za3*za4)+p4sq*za2*za3*za4-s23*za2*za3*za4+p
     -   4sq*za3*za4*zb2-p1sq*za4*zb2*zb3-s12*za4*zb2*zb3+p1sq*za2*za
     -   3*zb4+s12*za2*za3*zb4-p4sq*za2*zb3*zb4+p3sq*zb2*zb3*zb4-p4sq
     -   *zb2*zb3*zb4+s23*zb2*zb3*zb4+2*M**2*P(1)*P(2)*P(19)
       P(218) = -(p2sq*P(8))+M**2*P(1)*P(86)*P(207)
       P(219) = p2sq*za2*za3*za4-s12*za2*za3*za4+p1sq*za3*za4*zb2-p1sq
     -   *za2*zb3*zb4-p2sq*zb2*zb3*zb4+s12*zb2*zb3*zb4
       P(220) = -(s12*P(8))+M**2*P(1)*P(86)*P(207)
       P(221) = s12*P(8)-M**2*P(1)*P(86)*P(207)
       P(222) = zb2*zb3*P(28)+za2*za3*P(75)
       P(223) = za3*za4*zb2-za2*zb3*zb4
       P(224) = za1*za4*zb3-za3*zb1*zb4
       P(225) = za3*za4*zb1-za1*zb3*zb4
       P(226) = -(za3*za4*zb1)+za1*zb3*zb4
       P(227) = -(P(27)*P(224))+p3sq*P(225)+s23*P(226)
       P(228) = -(s12*za1*za3*za4)+s23*za1*za3*za4+p2sq*za3*za4*zb1-s2
     -   3*za3*za4*zb1-p2sq*za1*za4*zb3+s12*za1*za4*zb3+p2sq*za3*zb1*
     -   zb4-s12*za3*zb1*zb4-p2sq*za1*zb3*zb4+s23*za1*zb3*zb4+s12*zb1
     -   *zb3*zb4-s23*zb1*zb3*zb4+p1sq*P(2)*P(3)-p3sq*P(9)*P(35)
       P(229) = -(p4sq*P(8))+M**2*P(19)*P(187)
       P(230) = p4sq*P(8)+M**2*P(19)*P(171)
       P(231) = s23*za1*za2*za3*za4+s23*zb1*zb2*zb3*zb4+p1sq*P(8)-p4sq
     -   *P(8)-2*P(139)
       P(232) = -(p3sq*za1*za2*za3*za4)+p4sq*za1*za2*za3*za4-s12*za1*z
     -   a2*za3*za4-p3sq*zb1*zb2*zb3*zb4+p4sq*zb1*zb2*zb3*zb4-s12*zb1
     -   *zb2*zb3*zb4+2*P(139)
       P(233) = -(p1sq*za1*za2*za3*za4)+p4sq*za1*za2*za3*za4-s23*za1*z
     -   a2*za3*za4-p1sq*zb1*zb2*zb3*zb4+p4sq*zb1*zb2*zb3*zb4-s23*zb1
     -   *zb2*zb3*zb4+2*P(139)
       P(234) = s12*za1*za2*za3*za4+s12*zb1*zb2*zb3*zb4+p3sq*P(8)-p4sq
     -   *P(8)-2*P(139)
       P(235) = -(p4sq*za1*za3*za4)+s23*za3*za4*zb1+p1sq*za1*za4*zb3-p
     -   1sq*za3*zb1*zb4-s23*za1*zb3*zb4+p4sq*zb1*zb3*zb4+2*M**2*P(2)
     -   *P(9)*P(19)
       P(236) = za1*za3*za4-zb1*zb3*zb4
       P(237) = -(s12*za1*za3*za4)-s23*za1*za3*za4+p2sq*za3*za4*zb1+s2
     -   3*za3*za4*zb1+p2sq*za1*za4*zb3-s12*za1*za4*zb3-p2sq*za3*zb1*
     -   zb4+s12*za3*zb1*zb4-p2sq*za1*zb3*zb4-s23*za1*zb3*zb4+s12*zb1
     -   *zb3*zb4+s23*zb1*zb3*zb4+p3sq*P(9)*P(35)+p1sq*P(33)*P(61)
       P(238) = p3sq-p4sq-s12
       P(239) = 2*M**2*P(9)*P(19)*P(64)+P(8)*P(238)
       P(240) = -2*M**2*P(9)*P(19)*P(64)-P(8)*P(238)
       P(241) = -(s12*za1*za3*za4)-s23*za1*za3*za4+p2sq*za3*za4*zb1-s2
     -   3*za3*za4*zb1+p2sq*za1*za4*zb3+s12*za1*za4*zb3-p2sq*za3*zb1*
     -   zb4-s12*za3*zb1*zb4-p2sq*za1*zb3*zb4+s23*za1*zb3*zb4+s12*zb1
     -   *zb3*zb4+s23*zb1*zb3*zb4+p1sq*P(2)*P(3)+p3sq*P(86)*P(207)
       P(242) = P(27)*P(224)+s23*P(225)+p3sq*P(226)
       P(243) = p1sq-p4sq-s23
       P(244) = -2*M**2*P(2)*P(18)*P(19)-P(8)*P(243)
       P(245) = 2*M**2*P(2)*P(18)*P(19)+P(8)*P(243)
       P(246) = p1sq*za1*za3*za4-s12*za1*za3*za4+p2sq*za3*za4*zb1-p2sq
     -   *za1*zb3*zb4-p1sq*zb1*zb3*zb4+s12*zb1*zb3*zb4
       P(247) = p1sq*P(8)-M**2*P(9)*P(60)*P(207)
       P(248) = -(s12*P(8))+M**2*P(9)*P(60)*P(207)
       P(249) = -(p1sq*P(8))+M**2*P(9)*P(60)*P(207)
       P(250) = s12*P(8)-M**2*P(9)*P(60)*P(207)
       P(251) = p1sq*p3sq+p2sq*p4sq-s12*s23
       P(252) = p3sq*za1*za2*za3*za4+p4sq*za1*za2*za3*za4-s12*za1*za2*
     -   za3*za4-s23*za1*za2*za3*za4-p3sq*za2*za3*za4*zb1-p4sq*za1*za
     -   3*za4*zb2+s23*za3*za4*zb1*zb2+s12*za1*za4*zb2*zb3+s12*za2*za
     -   3*zb1*zb4+s23*za1*za2*zb3*zb4-p4sq*za2*zb1*zb3*zb4-p3sq*za1*
     -   zb2*zb3*zb4+p3sq*zb1*zb2*zb3*zb4+p4sq*zb1*zb2*zb3*zb4-s12*zb
     -   1*zb2*zb3*zb4-s23*zb1*zb2*zb3*zb4+p1sq*P(2)*P(4)+p2sq*P(19)*
     -   P(114)
       P(253) = 4*M**4*P(1)*P(2)*P(9)*P(19)+P(8)*P(251)-2*M**2*P(252)
       P(254) = p1sq+p2sq+p3sq+p4sq-2*s23
       P(255) = za3*za4*zb2+za2*za4*zb3+za2*za3*zb4+zb2*zb3*zb4-2*P(10
     -   1)
       P(256) = za2*za3*za4+za4*zb2*zb3+za3*zb2*zb4+za2*zb3*zb4-2*P(10
     -   3)
       P(257) = za1*P(255)+zb1*P(256)
       P(258) = P(8)*P(254)+2*M**2*P(257)
       FR(1) = (2*(B013R+B024R-4*(C124R(4,2)+C234R(4,2)-4*D1234R(22,4)
     -   ))*P(8)+C0124R*P(108)+C0234R*P(135)+C0123R*P(146)+C0134R*P(2
     -   13)+D01234R*P(253)+4*D1234R(7,2)*P(258))/2.
       FI(1) = (2*(B013I+B024I-4*(C124I(4,2)+C234I(4,2)-4*D1234I(22,4)
     -   ))*P(8)+C0124I*P(108)+C0234I*P(135)+C0123I*P(146)+C0134I*P(2
     -   13)+D01234I*P(253)+4*D1234I(7,2)*P(258))/2.
       F(1)=DCMPLX(FR(1),FI(1))
       P(259) = -4*M**4*P(1)*P(2)*P(9)*P(19)-P(8)*P(251)+2*M**2*P(252)
       P(260) = M**2*P(10)*P(49)
       P(261) = s12*za1*za2*za3*za4+s23*za1*za2*za3*za4-M**2*za2*za3*z
     -   a4*zb1-M**2*za1*za3*za4*zb2-M**2*za1*za2*za4*zb3-M**2*za4*zb
     -   1*zb2*zb3-M**2*za1*za2*za3*zb4-M**2*za3*zb1*zb2*zb4-M**2*za2
     -   *zb1*zb3*zb4-M**2*za1*zb2*zb3*zb4+s12*zb1*zb2*zb3*zb4+s23*zb
     -   1*zb2*zb3*zb4
       P(262) = p1sq*za1*za2*za3*za4+p2sq*za1*za2*za3*za4+p3sq*za1*za2
     -   *za3*za4+p4sq*za1*za2*za3*za4+p1sq*zb1*zb2*zb3*zb4+p2sq*zb1*
     -   zb2*zb3*zb4+p3sq*zb1*zb2*zb3*zb4+p4sq*zb1*zb2*zb3*zb4-4*P(26
     -   0)-2*P(261)
       FR(2) = (32*D1234R(22,4)*P(8)+C0123R*P(37)+2*(B013R+B024R)*P(41
     -   )+C0124R*P(109)+C0234R*P(132)+C0134R*P(210)+D01234R*P(259)-4
     -   *D1234R(7,2)*P(262))/2.
       FI(2) = (32*D1234I(22,4)*P(8)+C0123I*P(37)+2*(B013I+B024I)*P(41
     -   )+C0124I*P(109)+C0234I*P(132)+C0134I*P(210)+D01234I*P(259)-4
     -   *D1234I(7,2)*P(262))/2.
       F(2)=DCMPLX(FR(2),FI(2))
       P(263) = p1sq+p2sq+p3sq+p4sq-2*s12
       P(264) = za2*P(35)
       P(265) = za3*za4*zb2+za2*za4*zb3+za2*za3*zb4+zb2*zb3*zb4-2*P(26
     -   4)
       P(266) = zb2*P(35)
       P(267) = za2*za3*za4+za4*zb2*zb3+za3*zb2*zb4+za2*zb3*zb4-2*P(26
     -   6)
       P(268) = za1*P(265)+zb1*P(267)
       P(269) = P(8)*P(263)+2*M**2*P(268)
       FR(3) = (2*(B013R+B024R-4*(C123R(4,2)+C134R(4,2)-4*D1234R(22,4)
     -   ))*P(8)+C0124R*P(108)+C0234R*P(135)+C0123R*P(146)+C0134R*P(2
     -   13)+D01234R*P(253)+4*D1234R(7,2)*P(269))/2.
       FI(3) = (2*(B013I+B024I-4*(C123I(4,2)+C134I(4,2)-4*D1234I(22,4)
     -   ))*P(8)+C0124I*P(108)+C0234I*P(135)+C0123I*P(146)+C0134I*P(2
     -   13)+D01234I*P(253)+4*D1234I(7,2)*P(269))/2.
       F(3)=DCMPLX(FR(3),FI(3))
       P(270) = p1sq*za1*za2*za3*za4+p3sq*za1*za2*za3*za4-s12*za1*za2*
     -   za3*za4-s23*za1*za2*za3*za4+p1sq*zb1*zb2*zb3*zb4+p3sq*zb1*zb
     -   2*zb3*zb4-s12*zb1*zb2*zb3*zb4-s23*zb1*zb2*zb3*zb4+2*P(50)
       FR(4) = 16*(D1234R(11,3)-D1234R(13,3)+D1234R(16,4)+D1234R(18,4)
     -   -2*D1234R(20,4))*P(8)-2*(D1234R(1,1)+D1234R(1,2)-D1234R(3,1)
     -   +D1234R(3,2)-2*D1234R(5,2))*P(270)
       FI(4) = 16*(D1234I(11,3)-D1234I(13,3)+D1234I(16,4)+D1234I(18,4)
     -   -2*D1234I(20,4))*P(8)-2*(D1234I(1,1)+D1234I(1,2)-D1234I(3,1)
     -   +D1234I(3,2)-2*D1234I(5,2))*P(270)
       F(4)=DCMPLX(FR(4),FI(4))
       FR(5) = 2*(-2*(C234R(1,1)+C234R(1,2)-C234R(2,1)+C234R(2,2)-2*(C
     -   234R(3,2)+2*(D1234R(12,3)-D1234R(13,3)+D1234R(17,4)+D1234R(1
     -   8,4)-2*D1234R(21,4))))*P(8)+(D1234R(3,1)+2*D1234R(6,2))*P(10
     -   8)+(D1234R(2,1)+D1234R(2,2)+D1234R(3,2))*P(109))
       FI(5) = 2*(-2*(C234I(1,1)+C234I(1,2)-C234I(2,1)+C234I(2,2)-2*(C
     -   234I(3,2)+2*(D1234I(12,3)-D1234I(13,3)+D1234I(17,4)+D1234I(1
     -   8,4)-2*D1234I(21,4))))*P(8)+(D1234I(3,1)+2*D1234I(6,2))*P(10
     -   8)+(D1234I(2,1)+D1234I(2,2)+D1234I(3,2))*P(109))
       F(5)=DCMPLX(FR(5),FI(5))
       P(271) = p3sq*P(8)-s23*P(8)+2*M**2*P(1)*P(224)
       P(272) = -(p3sq*P(8))+s23*P(8)-2*M**2*P(1)*P(224)
       FR(6) = 2*((-C124R(1,1)+C124R(2,1)+C134R(1,1)-C134R(2,1)-2*(C12
     -   4R(1,2)+C124R(2,2)-2*(C124R(3,2)+D1234R(11,3)-D1234R(13,3)+2
     -   *(D1234R(16,4)+D1234R(18,4)-2*D1234R(20,4)))))*P(8)+(D1234R(
     -   1,2)+D1234R(3,2))*P(132)+2*D1234R(5,2)*P(135)+D1234R(1,1)*P(
     -   271)+D1234R(3,1)*P(272))
       FI(6) = 2*((-C124I(1,1)+C124I(2,1)+C134I(1,1)-C134I(2,1)-2*(C12
     -   4I(1,2)+C124I(2,2)-2*(C124I(3,2)+D1234I(11,3)-D1234I(13,3)+2
     -   *(D1234I(16,4)+D1234I(18,4)-2*D1234I(20,4)))))*P(8)+(D1234I(
     -   1,2)+D1234I(3,2))*P(132)+2*D1234I(5,2)*P(135)+D1234I(1,1)*P(
     -   271)+D1234I(3,1)*P(272))
       F(6)=DCMPLX(FR(6),FI(6))
       P(273) = p4sq*P(8)-2*M**2*P(19)*P(114)
       P(274) = -(P(8)*P(38))+2*M**2*P(9)*P(25)*P(127)
       P(275) = p3sq+p4sq-s12-2*s23
       P(276) = -(za4*zb2)-2*za2*zb4+zb2*zb4
       P(277) = zb2*zb3*P(25)+za3*P(276)
       P(278) = zb3*P(19)+za3*P(25)
       P(279) = -2*za4*zb2*zb3+za2*P(278)
       P(280) = zb1*P(277)+za1*P(279)
       P(281) = -(P(8)*P(275))-2*M**2*P(280)
       FR(7) = (C0123R+C0234R)*P(41)+2*((C134R(1,1)-C134R(2,1)+8*(D123
     -   4R(18,4)+D1234R(19,4)-D1234R(20,4)-D1234R(21,4)))*P(8)+(D123
     -   4R(3,2)+D1234R(4,2))*P(132)+(D1234R(5,2)+D1234R(6,2))*P(135)
     -   )+D1234R(2,1)*P(270)+D01234R*P(273)+D1234R(1,1)*P(274)+D1234
     -   R(3,1)*P(281)
       FI(7) = (C0123I+C0234I)*P(41)+2*((C134I(1,1)-C134I(2,1)+8*(D123
     -   4I(18,4)+D1234I(19,4)-D1234I(20,4)-D1234I(21,4)))*P(8)+(D123
     -   4I(3,2)+D1234I(4,2))*P(132)+(D1234I(5,2)+D1234I(6,2))*P(135)
     -   )+D1234I(2,1)*P(270)+D01234I*P(273)+D1234I(1,1)*P(274)+D1234
     -   I(3,1)*P(281)
       F(7)=DCMPLX(FR(7),FI(7))
       P(282) = p1sq*P(8)-s12*P(8)+2*M**2*P(1)*P(225)
       P(283) = -(p1sq*P(8))+s12*P(8)-2*M**2*P(1)*P(225)
       FR(8) = 2*((C124R(1,1)-C124R(2,1)-C134R(1,1)+C134R(2,1)-2*(C134
     -   R(1,2)+C134R(2,2)-2*(C134R(3,2)+D1234R(11,3)-D1234R(13,3)+2*
     -   (D1234R(16,4)+D1234R(18,4)-2*D1234R(20,4)))))*P(8)+(D1234R(1
     -   ,2)+D1234R(3,2))*P(37)+2*D1234R(5,2)*P(146)+D1234R(1,1)*P(28
     -   2)+D1234R(3,1)*P(283))
       FI(8) = 2*((C124I(1,1)-C124I(2,1)-C134I(1,1)+C134I(2,1)-2*(C134
     -   I(1,2)+C134I(2,2)-2*(C134I(3,2)+D1234I(11,3)-D1234I(13,3)+2*
     -   (D1234I(16,4)+D1234I(18,4)-2*D1234I(20,4)))))*P(8)+(D1234I(1
     -   ,2)+D1234I(3,2))*P(37)+2*D1234I(5,2)*P(146)+D1234I(1,1)*P(28
     -   2)+D1234I(3,1)*P(283))
       F(8)=DCMPLX(FR(8),FI(8))
       FR(9) = 16*(D1234R(1,2)+D1234R(1,4)+D1234R(3,2)+D1234R(3,4)-2*(
     -   -D1234R(1,3)+D1234R(3,3)+D1234R(5,2))-6*(D1234R(5,3)-D1234R(
     -   7,3)-D1234R(8,4))-4*(D1234R(5,4)+D1234R(12,4)))*P(8)
       FI(9) = 16*(D1234I(1,2)+D1234I(1,4)+D1234I(3,2)+D1234I(3,4)-2*(
     -   -D1234I(1,3)+D1234I(3,3)+D1234I(5,2))-6*(D1234I(5,3)-D1234I(
     -   7,3)-D1234I(8,4))-4*(D1234I(5,4)+D1234I(12,4)))*P(8)
       F(9)=DCMPLX(FR(9),FI(9))
       P(284) = -(p4sq*P(8))+2*M**2*P(19)*P(114)
       P(285) = P(8)*P(38)-2*M**2*P(9)*P(25)*P(127)
       P(286) = 2*p1sq+p3sq-p4sq-s12
       P(287) = za3*P(25)+zb3*P(111)
       P(288) = -2*zb2*zb3*zb4+za2*P(287)
       P(289) = za3*P(19)+zb3*P(25)
       P(290) = -2*za2*za3*za4+zb2*P(289)
       P(291) = za1*P(288)+zb1*P(290)
       P(292) = -(P(8)*P(286))+2*M**2*P(291)
       FR(10) = (C0123R+C0234R-2*(C134R(1,1)-C134R(2,1)+2*(C134R(1,2)+
     -   C134R(2,2)-2*(C134R(3,2)+D1234R(11,3)-D1234R(13,3)+2*(D1234R
     -   (18,4)+D1234R(19,4)-D1234R(20,4)-D1234R(21,4))))))*P(8)+2*((
     -   D1234R(3,2)+D1234R(4,2))*P(37)+(D1234R(5,2)+D1234R(6,2))*P(1
     -   46))+D1234R(2,1)*P(270)+D01234R*P(284)+D1234R(1,1)*P(285)+D1
     -   234R(3,1)*P(292)
       FI(10) = (C0123I+C0234I-2*(C134I(1,1)-C134I(2,1)+2*(C134I(1,2)+
     -   C134I(2,2)-2*(C134I(3,2)+D1234I(11,3)-D1234I(13,3)+2*(D1234I
     -   (18,4)+D1234I(19,4)-D1234I(20,4)-D1234I(21,4))))))*P(8)+2*((
     -   D1234I(3,2)+D1234I(4,2))*P(37)+(D1234I(5,2)+D1234I(6,2))*P(1
     -   46))+D1234I(2,1)*P(270)+D01234I*P(284)+D1234I(1,1)*P(285)+D1
     -   234I(3,1)*P(292)
       F(10)=DCMPLX(FR(10),FI(10))
       FR(11) = 8*(D1234R(3,3)-D1234R(7,3)-2*(-D1234R(3,4)+D1234R(12,4
     -   )))*P(8)
       FI(11) = 8*(D1234I(3,3)-D1234I(7,3)-2*(-D1234I(3,4)+D1234I(12,4
     -   )))*P(8)
       F(11)=DCMPLX(FR(11),FI(11))
       P(293) = -(p4sq*P(8))+s23*P(8)+2*M**2*P(9)*P(205)
       P(294) = p4sq*P(8)-s23*P(8)-2*M**2*P(9)*P(205)
       P(295) = p2sq*za1*za2*za3*za4+p4sq*za1*za2*za3*za4-s12*za1*za2*
     -   za3*za4-s23*za1*za2*za3*za4+p2sq*zb1*zb2*zb3*zb4+p4sq*zb1*zb
     -   2*zb3*zb4-s12*zb1*zb2*zb3*zb4-s23*zb1*zb2*zb3*zb4+2*P(139)
       FR(12) = 2*((-C134R(1,1)+C134R(2,1)+C234R(1,1)-C234R(2,1)+4*(D1
     -   234R(12,3)-D1234R(13,3)+2*(D1234R(17,4)+D1234R(18,4)-2*D1234
     -   R(21,4))))*P(8)+D1234R(2,1)*P(293)+D1234R(3,1)*P(294)-(D1234
     -   R(2,2)+D1234R(3,2)-2*D1234R(6,2))*P(295))
       FI(12) = 2*((-C134I(1,1)+C134I(2,1)+C234I(1,1)-C234I(2,1)+4*(D1
     -   234I(12,3)-D1234I(13,3)+2*(D1234I(17,4)+D1234I(18,4)-2*D1234
     -   I(21,4))))*P(8)+D1234I(2,1)*P(293)+D1234I(3,1)*P(294)-(D1234
     -   I(2,2)+D1234I(3,2)-2*D1234I(6,2))*P(295))
       F(12)=DCMPLX(FR(12),FI(12))
       P(296) = p1sq*P(8)-2*M**2*P(9)*P(206)
       FR(13) = 2*((-C134R(1,1)+C134R(2,1)-C234R(1,1)+C234R(2,1)-2*(C1
     -   34R(1,2)+C134R(2,2)-2*(C134R(3,2)+D1234R(12,3)-D1234R(13,3)+
     -   2*(D1234R(17,4)+D1234R(18,4)-2*D1234R(21,4)))))*P(8)+(D1234R
     -   (2,2)+D1234R(3,2))*P(37)+2*D1234R(6,2)*P(146)+(D1234R(2,1)-D
     -   1234R(3,1))*P(296))
       FI(13) = 2*((-C134I(1,1)+C134I(2,1)-C234I(1,1)+C234I(2,1)-2*(C1
     -   34I(1,2)+C134I(2,2)-2*(C134I(3,2)+D1234I(12,3)-D1234I(13,3)+
     -   2*(D1234I(17,4)+D1234I(18,4)-2*D1234I(21,4)))))*P(8)+(D1234I
     -   (2,2)+D1234I(3,2))*P(37)+2*D1234I(6,2)*P(146)+(D1234I(2,1)-D
     -   1234I(3,1))*P(296))
       F(13)=DCMPLX(FR(13),FI(13))
       FR(14) = -16*(D1234R(3,3)-D1234R(3,4)+D1234R(5,3)+D1234R(5,4)-2
     -   *D1234R(7,3)+3*(-D1234R(8,4)+D1234R(12,4)))*P(8)
       FI(14) = -16*(D1234I(3,3)-D1234I(3,4)+D1234I(5,3)+D1234I(5,4)-2
     -   *D1234I(7,3)+3*(-D1234I(8,4)+D1234I(12,4)))*P(8)
       F(14)=DCMPLX(FR(14),FI(14))
       FR(15) = 4*(D1234R(1,1)+D1234R(1,2)-D1234R(3,1)+D1234R(3,2)-2*(
     -   D1234R(5,2)+2*(D1234R(3,3)-D1234R(3,4)-D1234R(7,3)+D1234R(7,
     -   4)-D1234R(8,4)-D1234R(9,3)+D1234R(10,3)-2*(D1234R(11,4)-D123
     -   4R(12,4))+D1234R(15,4))))*P(8)
       FI(15) = 4*(D1234I(1,1)+D1234I(1,2)-D1234I(3,1)+D1234I(3,2)-2*(
     -   D1234I(5,2)+2*(D1234I(3,3)-D1234I(3,4)-D1234I(7,3)+D1234I(7,
     -   4)-D1234I(8,4)-D1234I(9,3)+D1234I(10,3)-2*(D1234I(11,4)-D123
     -   4I(12,4))+D1234I(15,4))))*P(8)
       F(15)=DCMPLX(FR(15),FI(15))
       FR(16) = -8*(D1234R(3,2)+D1234R(3,3)-D1234R(5,2)-D1234R(5,3)-2*
     -   (D1234R(3,4)+D1234R(8,4))+4*D1234R(12,4))*P(8)
       FI(16) = -8*(D1234I(3,2)+D1234I(3,3)-D1234I(5,2)-D1234I(5,3)-2*
     -   (D1234I(3,4)+D1234I(8,4))+4*D1234I(12,4))*P(8)
       F(16)=DCMPLX(FR(16),FI(16))
       FR(17) = 16*(D1234R(3,2)+D1234R(3,4)+D1234R(4,2)+D1234R(4,4)-D1
     -   234R(5,2)-D1234R(5,4)-D1234R(6,2)-2*(D1234R(3,3)-D1234R(4,3)
     -   +D1234R(5,3)-D1234R(9,3))-4*(-D1234R(7,3)+D1234R(10,3))+3*(-
     -   D1234R(7,4)+D1234R(8,4)+D1234R(11,4)-D1234R(12,4))-D1234R(15
     -   ,4))*P(8)
       FI(17) = 16*(D1234I(3,2)+D1234I(3,4)+D1234I(4,2)+D1234I(4,4)-D1
     -   234I(5,2)-D1234I(5,4)-D1234I(6,2)-2*(D1234I(3,3)-D1234I(4,3)
     -   +D1234I(5,3)-D1234I(9,3))-4*(-D1234I(7,3)+D1234I(10,3))+3*(-
     -   D1234I(7,4)+D1234I(8,4)+D1234I(11,4)-D1234I(12,4))-D1234I(15
     -   ,4))*P(8)
       F(17)=DCMPLX(FR(17),FI(17))
       FR(18) = -8*(D1234R(1,1)+D1234R(1,3)-D1234R(3,1)+3*D1234R(3,3)+
     -   D1234R(5,3)+2*(D1234R(1,2)-D1234R(3,4)-D1234R(5,2)+D1234R(5,
     -   4))-5*D1234R(7,3)+6*(-D1234R(8,4)+D1234R(12,4)))*P(8)
       FI(18) = -8*(D1234I(1,1)+D1234I(1,3)-D1234I(3,1)+3*D1234I(3,3)+
     -   D1234I(5,3)+2*(D1234I(1,2)-D1234I(3,4)-D1234I(5,2)+D1234I(5,
     -   4))-5*D1234I(7,3)+6*(-D1234I(8,4)+D1234I(12,4)))*P(8)
       F(18)=DCMPLX(FR(18),FI(18))
       FR(19) = -4*(D1234R(3,2)-D1234R(6,2)-2*(D1234R(3,3)-D1234R(7,3)
     -   -D1234R(9,3)+D1234R(10,3)+2*(D1234R(3,4)+D1234R(11,4)-D1234R
     -   (12,4)-D1234R(15,4))))*P(8)
       FI(19) = -4*(D1234I(3,2)-D1234I(6,2)-2*(D1234I(3,3)-D1234I(7,3)
     -   -D1234I(9,3)+D1234I(10,3)+2*(D1234I(3,4)+D1234I(11,4)-D1234I
     -   (12,4)-D1234I(15,4))))*P(8)
       F(19)=DCMPLX(FR(19),FI(19))
       FR(20) = 8*(D1234R(3,3)+D1234R(8,3)+2*(D1234R(3,4)-D1234R(9,3)+
     -   D1234R(14,4)-2*D1234R(15,4)))*P(8)
       FI(20) = 8*(D1234I(3,3)+D1234I(8,3)+2*(D1234I(3,4)-D1234I(9,3)+
     -   D1234I(14,4)-2*D1234I(15,4)))*P(8)
       F(20)=DCMPLX(FR(20),FI(20))
       FR(21) = 16*(D1234R(3,3)+D1234R(3,4)-D1234R(9,3)-D1234R(15,4))*
     -   P(8)
       FI(21) = 16*(D1234I(3,3)+D1234I(3,4)-D1234I(9,3)-D1234I(15,4))*
     -   P(8)
       F(21)=DCMPLX(FR(21),FI(21))
       FR(22) = 8*(D1234R(1,2)+D1234R(1,3)+D1234R(3,2)-3*D1234R(3,3)-5
     -   *D1234R(5,3)+7*D1234R(7,3)-4*D1234R(10,3)-2*(-D1234R(3,4)-D1
     -   234R(4,3)-D1234R(4,4)+D1234R(5,2)+D1234R(5,4)-D1234R(9,3)-3*
     -   (-D1234R(7,4)+D1234R(8,4)+D1234R(11,4)-D1234R(12,4))+D1234R(
     -   15,4)))*P(8)
       FI(22) = 8*(D1234I(1,2)+D1234I(1,3)+D1234I(3,2)-3*D1234I(3,3)-5
     -   *D1234I(5,3)+7*D1234I(7,3)-4*D1234I(10,3)-2*(-D1234I(3,4)-D1
     -   234I(4,3)-D1234I(4,4)+D1234I(5,2)+D1234I(5,4)-D1234I(9,3)-3*
     -   (-D1234I(7,4)+D1234I(8,4)+D1234I(11,4)-D1234I(12,4))+D1234I(
     -   15,4)))*P(8)
       F(22)=DCMPLX(FR(22),FI(22))
       FR(23) = 8*(D1234R(3,2)-3*D1234R(3,3)+D1234R(4,2)+D1234R(4,3)-D
     -   1234R(5,2)-D1234R(5,3)-D1234R(6,2)+5*D1234R(9,3)-6*D1234R(10
     -   ,3)+4*(D1234R(7,3)-D1234R(10,4))+2*(D1234R(3,4)+D1234R(6,3)+
     -   D1234R(6,4)-D1234R(8,3)+D1234R(8,4)+4*D1234R(11,4)+D1234R(14
     -   ,4)-2*(D1234R(7,4)+D1234R(12,4)+D1234R(15,4))))*P(8)
       FI(23) = 8*(D1234I(3,2)-3*D1234I(3,3)+D1234I(4,2)+D1234I(4,3)-D
     -   1234I(5,2)-D1234I(5,3)-D1234I(6,2)+5*D1234I(9,3)-6*D1234I(10
     -   ,3)+4*(D1234I(7,3)-D1234I(10,4))+2*(D1234I(3,4)+D1234I(6,3)+
     -   D1234I(6,4)-D1234I(8,3)+D1234I(8,4)+4*D1234I(11,4)+D1234I(14
     -   ,4)-2*(D1234I(7,4)+D1234I(12,4)+D1234I(15,4))))*P(8)
       F(23)=DCMPLX(FR(23),FI(23))
       FR(24) = -4*(D1234R(1,1)+D1234R(1,2)-D1234R(3,1)+D1234R(3,2)+2*
     -   (D1234R(4,2)+D1234R(4,3)-D1234R(5,2)-D1234R(6,2)-D1234R(9,3)
     -   -4*(D1234R(11,4)-D1234R(12,4))+2*(D1234R(3,3)-D1234R(3,4)-D1
     -   234R(7,3)+D1234R(7,4)-D1234R(8,4)+D1234R(15,4))))*P(8)
       FI(24) = -4*(D1234I(1,1)+D1234I(1,2)-D1234I(3,1)+D1234I(3,2)+2*
     -   (D1234I(4,2)+D1234I(4,3)-D1234I(5,2)-D1234I(6,2)-D1234I(9,3)
     -   -4*(D1234I(11,4)-D1234I(12,4))+2*(D1234I(3,3)-D1234I(3,4)-D1
     -   234I(7,3)+D1234I(7,4)-D1234I(8,4)+D1234I(15,4))))*P(8)
       F(24)=DCMPLX(FR(24),FI(24))
       FR(25) = 8*(D1234R(3,2)+D1234R(4,2)-D1234R(5,2)-D1234R(6,2)-3*(
     -   D1234R(3,3)-D1234R(4,3)+D1234R(5,3)-D1234R(9,3))-6*(-D1234R(
     -   7,3)+D1234R(7,4)-D1234R(8,4)+D1234R(10,3)-D1234R(11,4)+D1234
     -   R(12,4))-2*(-D1234R(3,4)-D1234R(4,4)+D1234R(5,4)+D1234R(15,4
     -   )))*P(8)
       FI(25) = 8*(D1234I(3,2)+D1234I(4,2)-D1234I(5,2)-D1234I(6,2)-3*(
     -   D1234I(3,3)-D1234I(4,3)+D1234I(5,3)-D1234I(9,3))-6*(-D1234I(
     -   7,3)+D1234I(7,4)-D1234I(8,4)+D1234I(10,3)-D1234I(11,4)+D1234
     -   I(12,4))-2*(-D1234I(3,4)-D1234I(4,4)+D1234I(5,4)+D1234I(15,4
     -   )))*P(8)
       F(25)=DCMPLX(FR(25),FI(25))
       FR(26) = 8*(D1234R(2,2)+D1234R(3,2)-3*(D1234R(3,3)-D1234R(6,3)-
     -   D1234R(7,3)+D1234R(8,3))-6*D1234R(10,3)-2*(-D1234R(3,4)+D123
     -   4R(6,2)-D1234R(6,4)-D1234R(8,4)-3*D1234R(9,3)-4*D1234R(11,4)
     -   -D1234R(14,4)+2*(D1234R(7,4)+D1234R(10,4)+D1234R(12,4)+D1234
     -   R(15,4))))*P(8)
       FI(26) = 8*(D1234I(2,2)+D1234I(3,2)-3*(D1234I(3,3)-D1234I(6,3)-
     -   D1234I(7,3)+D1234I(8,3))-6*D1234I(10,3)-2*(-D1234I(3,4)+D123
     -   4I(6,2)-D1234I(6,4)-D1234I(8,4)-3*D1234I(9,3)-4*D1234I(11,4)
     -   -D1234I(14,4)+2*(D1234I(7,4)+D1234I(10,4)+D1234I(12,4)+D1234
     -   I(15,4))))*P(8)
       F(26)=DCMPLX(FR(26),FI(26))
       FR(27) = -4*(D1234R(2,1)-D1234R(3,1)+D1234R(3,2)+3*(D1234R(4,2)
     -   -D1234R(5,2))-D1234R(6,2)+2*(D1234R(4,3)-D1234R(5,3)-D1234R(
     -   7,3)+D1234R(10,3)-4*(D1234R(11,4)-D1234R(12,4))+2*(D1234R(3,
     -   3)-D1234R(3,4)+D1234R(7,4)-D1234R(8,4)-D1234R(9,3)+D1234R(15
     -   ,4))))*P(8)
       FI(27) = -4*(D1234I(2,1)-D1234I(3,1)+D1234I(3,2)+3*(D1234I(4,2)
     -   -D1234I(5,2))-D1234I(6,2)+2*(D1234I(4,3)-D1234I(5,3)-D1234I(
     -   7,3)+D1234I(10,3)-4*(D1234I(11,4)-D1234I(12,4))+2*(D1234I(3,
     -   3)-D1234I(3,4)+D1234I(7,4)-D1234I(8,4)-D1234I(9,3)+D1234I(15
     -   ,4))))*P(8)
       F(27)=DCMPLX(FR(27),FI(27))
       FR(28) = 8*(D1234R(1,2)+D1234R(1,3)-4*D1234R(3,3)+D1234R(4,2)-D
     -   1234R(6,2)+9*D1234R(7,3)-3*(-D1234R(4,3)+D1234R(5,2)-D1234R(
     -   9,3))-6*(D1234R(5,3)+D1234R(10,3))+2*(D1234R(3,2)+D1234R(3,4
     -   )+D1234R(4,4)-D1234R(5,4)+3*(-D1234R(7,4)+D1234R(8,4)+D1234R
     -   (11,4)-D1234R(12,4))-D1234R(15,4)))*P(8)
       FI(28) = 8*(D1234I(1,2)+D1234I(1,3)-4*D1234I(3,3)+D1234I(4,2)-D
     -   1234I(6,2)+9*D1234I(7,3)-3*(-D1234I(4,3)+D1234I(5,2)-D1234I(
     -   9,3))-6*(D1234I(5,3)+D1234I(10,3))+2*(D1234I(3,2)+D1234I(3,4
     -   )+D1234I(4,4)-D1234I(5,4)+3*(-D1234I(7,4)+D1234I(8,4)+D1234I
     -   (11,4)-D1234I(12,4))-D1234I(15,4)))*P(8)
       F(28)=DCMPLX(FR(28),FI(28))
       FR(29) = 8*(D1234R(2,2)-4*D1234R(3,3)+D1234R(4,2)+D1234R(4,3)-D
     -   1234R(5,2)-D1234R(5,3)+5*D1234R(7,3)-3*(D1234R(6,2)-D1234R(6
     -   ,3)+D1234R(8,3))+7*D1234R(9,3)-8*D1234R(10,3)+2*(D1234R(3,2)
     -   +D1234R(3,4)+D1234R(6,4)+D1234R(8,4)+4*D1234R(11,4)+D1234R(1
     -   4,4)-2*(D1234R(7,4)+D1234R(10,4)+D1234R(12,4)+D1234R(15,4)))
     -   )*P(8)
       FI(29) = 8*(D1234I(2,2)-4*D1234I(3,3)+D1234I(4,2)+D1234I(4,3)-D
     -   1234I(5,2)-D1234I(5,3)+5*D1234I(7,3)-3*(D1234I(6,2)-D1234I(6
     -   ,3)+D1234I(8,3))+7*D1234I(9,3)-8*D1234I(10,3)+2*(D1234I(3,2)
     -   +D1234I(3,4)+D1234I(6,4)+D1234I(8,4)+4*D1234I(11,4)+D1234I(1
     -   4,4)-2*(D1234I(7,4)+D1234I(10,4)+D1234I(12,4)+D1234I(15,4)))
     -   )*P(8)
       F(29)=DCMPLX(FR(29),FI(29))
       FR(30) = -4*(D1234R(1,1)+D1234R(1,2)+D1234R(2,1)+3*(D1234R(4,2)
     -   -D1234R(5,2))-D1234R(6,2)-2*(D1234R(3,1)-D1234R(4,3)-3*(D123
     -   4R(3,3)-D1234R(7,3))-D1234R(10,3)+2*(D1234R(3,4)-D1234R(7,4)
     -   +D1234R(8,4)+D1234R(9,3)+2*(D1234R(11,4)-D1234R(12,4))-D1234
     -   R(15,4))))*P(8)
       FI(30) = -4*(D1234I(1,1)+D1234I(1,2)+D1234I(2,1)+3*(D1234I(4,2)
     -   -D1234I(5,2))-D1234I(6,2)-2*(D1234I(3,1)-D1234I(4,3)-3*(D123
     -   4I(3,3)-D1234I(7,3))-D1234I(10,3)+2*(D1234I(3,4)-D1234I(7,4)
     -   +D1234I(8,4)+D1234I(9,3)+2*(D1234I(11,4)-D1234I(12,4))-D1234
     -   I(15,4))))*P(8)
       F(30)=DCMPLX(FR(30),FI(30))
       FR(31) = -4*(D1234R(2,2)+D1234R(4,2)-D1234R(5,2)-3*D1234R(6,2)+
     -   2*(D1234R(3,2)+D1234R(6,3)-D1234R(8,3)+D1234R(9,3)-D1234R(10
     -   ,3))+4*(-D1234R(3,4)+D1234R(10,4)+D1234R(12,4)-D1234R(14,4))
     -   -8*(D1234R(11,4)-D1234R(15,4)))*P(8)
       FI(31) = -4*(D1234I(2,2)+D1234I(4,2)-D1234I(5,2)-3*D1234I(6,2)+
     -   2*(D1234I(3,2)+D1234I(6,3)-D1234I(8,3)+D1234I(9,3)-D1234I(10
     -   ,3))+4*(-D1234I(3,4)+D1234I(10,4)+D1234I(12,4)-D1234I(14,4))
     -   -8*(D1234I(11,4)-D1234I(15,4)))*P(8)
       F(31)=DCMPLX(FR(31),FI(31))
       FR(32) = -8*(D1234R(2,2)+D1234R(2,3)+D1234R(3,2)+D1234R(9,3)-2*
     -   (D1234R(3,4)+D1234R(6,2)+D1234R(8,3)-D1234R(13,4))-6*(D1234R
     -   (14,4)-D1234R(15,4)))*P(8)
       FI(32) = -8*(D1234I(2,2)+D1234I(2,3)+D1234I(3,2)+D1234I(9,3)-2*
     -   (D1234I(3,4)+D1234I(6,2)+D1234I(8,3)-D1234I(13,4))-6*(D1234I
     -   (14,4)-D1234I(15,4)))*P(8)
       F(32)=DCMPLX(FR(32),FI(32))
       FR(33) = -8*(D1234R(3,2)-D1234R(3,3)-D1234R(6,2)+3*D1234R(9,3)-
     -   2*(D1234R(3,4)+D1234R(8,3)+D1234R(14,4))+4*D1234R(15,4))*P(8
     -   )
       FI(33) = -8*(D1234I(3,2)-D1234I(3,3)-D1234I(6,2)+3*D1234I(9,3)-
     -   2*(D1234I(3,4)+D1234I(8,3)+D1234I(14,4))+4*D1234I(15,4))*P(8
     -   )
       F(33)=DCMPLX(FR(33),FI(33))
       FR(34) = -4*(D1234R(2,1)-D1234R(3,1)+D1234R(3,2)+D1234R(4,2)-D1
     -   234R(5,2)-D1234R(6,2)+2*(D1234R(3,3)-D1234R(7,3)-D1234R(9,3)
     -   +D1234R(10,3))+4*(-D1234R(3,4)+D1234R(10,4)+D1234R(12,4)-D12
     -   34R(14,4))-8*(D1234R(11,4)-D1234R(15,4)))*P(8)
       FI(34) = -4*(D1234I(2,1)-D1234I(3,1)+D1234I(3,2)+D1234I(4,2)-D1
     -   234I(5,2)-D1234I(6,2)+2*(D1234I(3,3)-D1234I(7,3)-D1234I(9,3)
     -   +D1234I(10,3))+4*(-D1234I(3,4)+D1234I(10,4)+D1234I(12,4)-D12
     -   34I(14,4))-8*(D1234I(11,4)-D1234I(15,4)))*P(8)
       F(34)=DCMPLX(FR(34),FI(34))
       FR(35) = -8*(D1234R(3,3)+D1234R(8,3)-2*(D1234R(3,4)+D1234R(9,3)
     -   -D1234R(13,4)+3*(D1234R(14,4)-D1234R(15,4))))*P(8)
       FI(35) = -8*(D1234I(3,3)+D1234I(8,3)-2*(D1234I(3,4)+D1234I(9,3)
     -   -D1234I(13,4)+3*(D1234I(14,4)-D1234I(15,4))))*P(8)
       F(35)=DCMPLX(FR(35),FI(35))
       FR(36) = -8*(D1234R(3,2)-D1234R(6,2)-D1234R(8,3)+D1234R(9,3)-2*
     -   (D1234R(3,4)+D1234R(14,4)-2*D1234R(15,4)))*P(8)
       FI(36) = -8*(D1234I(3,2)-D1234I(6,2)-D1234I(8,3)+D1234I(9,3)-2*
     -   (D1234I(3,4)+D1234I(14,4)-2*D1234I(15,4)))*P(8)
       F(36)=DCMPLX(FR(36),FI(36))
       FR(37) = 4*(D1234R(2,1)+D1234R(2,2)-D1234R(3,1)+D1234R(3,2)-2*(
     -   D1234R(6,2)+2*(D1234R(3,3)-D1234R(3,4)-D1234R(7,3)-D1234R(9,
     -   3)+D1234R(10,3)+D1234R(10,4)+D1234R(12,4)-D1234R(14,4)-2*(D1
     -   234R(11,4)-D1234R(15,4)))))*P(8)
       FI(37) = 4*(D1234I(2,1)+D1234I(2,2)-D1234I(3,1)+D1234I(3,2)-2*(
     -   D1234I(6,2)+2*(D1234I(3,3)-D1234I(3,4)-D1234I(7,3)-D1234I(9,
     -   3)+D1234I(10,3)+D1234I(10,4)+D1234I(12,4)-D1234I(14,4)-2*(D1
     -   234I(11,4)-D1234I(15,4)))))*P(8)
       F(37)=DCMPLX(FR(37),FI(37))
       FR(38) = -16*(D1234R(3,3)-D1234R(3,4)+D1234R(8,3)-2*D1234R(9,3)
     -   +D1234R(13,4)-3*(D1234R(14,4)-D1234R(15,4)))*P(8)
       FI(38) = -16*(D1234I(3,3)-D1234I(3,4)+D1234I(8,3)-2*D1234I(9,3)
     -   +D1234I(13,4)-3*(D1234I(14,4)-D1234I(15,4)))*P(8)
       F(38)=DCMPLX(FR(38),FI(38))
       FR(39) = -8*(D1234R(3,2)+D1234R(3,3)-D1234R(6,2)-D1234R(8,3)-2*
     -   (D1234R(3,4)+D1234R(14,4)-2*D1234R(15,4)))*P(8)
       FI(39) = -8*(D1234I(3,2)+D1234I(3,3)-D1234I(6,2)-D1234I(8,3)-2*
     -   (D1234I(3,4)+D1234I(14,4)-2*D1234I(15,4)))*P(8)
       F(39)=DCMPLX(FR(39),FI(39))
       FR(40) = -8*(-D1234R(4,3)+D1234R(5,3)-D1234R(6,3)+D1234R(8,3)-3
     -   *(D1234R(7,3)+D1234R(9,3))+4*(D1234R(10,3)+D1234R(10,4))+2*(
     -   D1234R(3,3)-D1234R(3,4)-D1234R(6,4)-D1234R(8,4)-4*D1234R(11,
     -   4)-D1234R(14,4)+2*(D1234R(7,4)+D1234R(12,4)+D1234R(15,4))))*
     -   P(8)
       FI(40) = -8*(-D1234I(4,3)+D1234I(5,3)-D1234I(6,3)+D1234I(8,3)-3
     -   *(D1234I(7,3)+D1234I(9,3))+4*(D1234I(10,3)+D1234I(10,4))+2*(
     -   D1234I(3,3)-D1234I(3,4)-D1234I(6,4)-D1234I(8,4)-4*D1234I(11,
     -   4)-D1234I(14,4)+2*(D1234I(7,4)+D1234I(12,4)+D1234I(15,4))))*
     -   P(8)
       F(40)=DCMPLX(FR(40),FI(40))
       FR(41) = -8*(-D1234R(2,3)-D1234R(6,3)-D1234R(7,3)+4*D1234R(8,3)
     -   -5*D1234R(9,3)+2*(D1234R(3,3)-D1234R(3,4)-D1234R(9,4)+D1234R
     -   (10,3)+D1234R(12,4)+D1234R(13,4))+6*(D1234R(10,4)-D1234R(11,
     -   4)-D1234R(14,4)+D1234R(15,4)))*P(8)
       FI(41) = -8*(-D1234I(2,3)-D1234I(6,3)-D1234I(7,3)+4*D1234I(8,3)
     -   -5*D1234I(9,3)+2*(D1234I(3,3)-D1234I(3,4)-D1234I(9,4)+D1234I
     -   (10,3)+D1234I(12,4)+D1234I(13,4))+6*(D1234I(10,4)-D1234I(11,
     -   4)-D1234I(14,4)+D1234I(15,4)))*P(8)
       F(41)=DCMPLX(FR(41),FI(41))
       FR(42) = -4*(D1234R(2,2)+D1234R(4,2)-D1234R(5,2)-3*D1234R(6,2)+
     -   2*(D1234R(3,2)+D1234R(3,3)+D1234R(6,3)-D1234R(9,3)-D1234R(10
     -   ,3))+4*(-D1234R(3,4)+D1234R(10,4)+D1234R(12,4)-D1234R(14,4))
     -   -8*(D1234R(11,4)-D1234R(15,4)))*P(8)
       FI(42) = -4*(D1234I(2,2)+D1234I(4,2)-D1234I(5,2)-3*D1234I(6,2)+
     -   2*(D1234I(3,2)+D1234I(3,3)+D1234I(6,3)-D1234I(9,3)-D1234I(10
     -   ,3))+4*(-D1234I(3,4)+D1234I(10,4)+D1234I(12,4)-D1234I(14,4))
     -   -8*(D1234I(11,4)-D1234I(15,4)))*P(8)
       F(42)=DCMPLX(FR(42),FI(42))
       FR(43) = 8*(D1234R(3,2)-3*D1234R(3,3)+D1234R(4,2)-D1234R(5,2)-D
     -   1234R(6,2)+D1234R(6,3)+5*D1234R(7,3)-D1234R(8,3)-6*D1234R(10
     -   ,3)+2*(D1234R(3,4)+D1234R(4,3)-D1234R(5,3)+D1234R(6,4)+D1234
     -   R(8,4)+4*D1234R(11,4)+D1234R(14,4)+2*(-D1234R(7,4)+D1234R(9,
     -   3)-D1234R(10,4)-D1234R(12,4)-D1234R(15,4))))*P(8)
       FI(43) = 8*(D1234I(3,2)-3*D1234I(3,3)+D1234I(4,2)-D1234I(5,2)-D
     -   1234I(6,2)+D1234I(6,3)+5*D1234I(7,3)-D1234I(8,3)-6*D1234I(10
     -   ,3)+2*(D1234I(3,4)+D1234I(4,3)-D1234I(5,3)+D1234I(6,4)+D1234
     -   I(8,4)+4*D1234I(11,4)+D1234I(14,4)+2*(-D1234I(7,4)+D1234I(9,
     -   3)-D1234I(10,4)-D1234I(12,4)-D1234I(15,4))))*P(8)
       F(43)=DCMPLX(FR(43),FI(43))
       FR(44) = 8*(D1234R(2,2)+D1234R(2,3)+D1234R(3,2)-3*D1234R(3,3)-5
     -   *D1234R(8,3)+7*D1234R(9,3)-4*D1234R(10,3)-2*(-D1234R(3,4)+D1
     -   234R(6,2)-D1234R(6,3)-D1234R(7,3)-D1234R(9,4)+D1234R(12,4)+D
     -   1234R(13,4))-6*(D1234R(10,4)-D1234R(11,4)-D1234R(14,4)+D1234
     -   R(15,4)))*P(8)
       FI(44) = 8*(D1234I(2,2)+D1234I(2,3)+D1234I(3,2)-3*D1234I(3,3)-5
     -   *D1234I(8,3)+7*D1234I(9,3)-4*D1234I(10,3)-2*(-D1234I(3,4)+D1
     -   234I(6,2)-D1234I(6,3)-D1234I(7,3)-D1234I(9,4)+D1234I(12,4)+D
     -   1234I(13,4))-6*(D1234I(10,4)-D1234I(11,4)-D1234I(14,4)+D1234
     -   I(15,4)))*P(8)
       F(44)=DCMPLX(FR(44),FI(44))
       FR(45) = -4*(D1234R(2,1)+D1234R(2,2)-D1234R(3,1)+D1234R(3,2)+2*
     -   (D1234R(4,2)-D1234R(5,2)-D1234R(6,2)+D1234R(6,3)-D1234R(7,3)
     -   +2*(D1234R(3,3)-D1234R(3,4)-D1234R(9,3)+D1234R(10,4)+D1234R(
     -   12,4)-D1234R(14,4)-2*(D1234R(11,4)-D1234R(15,4)))))*P(8)
       FI(45) = -4*(D1234I(2,1)+D1234I(2,2)-D1234I(3,1)+D1234I(3,2)+2*
     -   (D1234I(4,2)-D1234I(5,2)-D1234I(6,2)+D1234I(6,3)-D1234I(7,3)
     -   +2*(D1234I(3,3)-D1234I(3,4)-D1234I(9,3)+D1234I(10,4)+D1234I(
     -   12,4)-D1234I(14,4)-2*(D1234I(11,4)-D1234I(15,4)))))*P(8)
       F(45)=DCMPLX(FR(45),FI(45))
       FR(46) = 8*(D1234R(2,2)+D1234R(3,2)-3*D1234R(3,3)+D1234R(6,3)+D
     -   1234R(7,3)-7*D1234R(8,3)-2*(-D1234R(2,3)-D1234R(3,4)+D1234R(
     -   6,2)-4*D1234R(9,3)-D1234R(9,4)+D1234R(10,3)+D1234R(12,4)+D12
     -   34R(13,4)+3*(D1234R(10,4)-D1234R(11,4)-D1234R(14,4)+D1234R(1
     -   5,4))))*P(8)
       FI(46) = 8*(D1234I(2,2)+D1234I(3,2)-3*D1234I(3,3)+D1234I(6,3)+D
     -   1234I(7,3)-7*D1234I(8,3)-2*(-D1234I(2,3)-D1234I(3,4)+D1234I(
     -   6,2)-4*D1234I(9,3)-D1234I(9,4)+D1234I(10,3)+D1234I(12,4)+D12
     -   34I(13,4)+3*(D1234I(10,4)-D1234I(11,4)-D1234I(14,4)+D1234I(1
     -   5,4))))*P(8)
       F(46)=DCMPLX(FR(46),FI(46))
       FR(47) = -4*(D1234R(2,1)-D1234R(3,1)+D1234R(3,2)+D1234R(4,2)-D1
     -   234R(5,2)-3*D1234R(6,2)+2*(D1234R(2,2)+D1234R(6,3)+D1234R(8,
     -   3)-3*D1234R(9,3)-D1234R(10,3)+2*(D1234R(3,3)-D1234R(3,4)+D12
     -   34R(10,4)+D1234R(12,4)-D1234R(14,4)-2*(D1234R(11,4)-D1234R(1
     -   5,4)))))*P(8)
       FI(47) = -4*(D1234I(2,1)-D1234I(3,1)+D1234I(3,2)+D1234I(4,2)-D1
     -   234I(5,2)-3*D1234I(6,2)+2*(D1234I(2,2)+D1234I(6,3)+D1234I(8,
     -   3)-3*D1234I(9,3)-D1234I(10,3)+2*(D1234I(3,3)-D1234I(3,4)+D12
     -   34I(10,4)+D1234I(12,4)-D1234I(14,4)-2*(D1234I(11,4)-D1234I(1
     -   5,4)))))*P(8)
       F(47)=DCMPLX(FR(47),FI(47))
       FR(48) = -16*(D1234R(3,3)-D1234R(3,4)-D1234R(6,3)-D1234R(7,3)+D
     -   1234R(8,3)-D1234R(9,4)+2*(-D1234R(9,3)+D1234R(10,3))+D1234R(
     -   12,4)+D1234R(13,4)+3*(D1234R(10,4)-D1234R(11,4)-D1234R(14,4)
     -   +D1234R(15,4)))*P(8)
       FI(48) = -16*(D1234I(3,3)-D1234I(3,4)-D1234I(6,3)-D1234I(7,3)+D
     -   1234I(8,3)-D1234I(9,4)+2*(-D1234I(9,3)+D1234I(10,3))+D1234I(
     -   12,4)+D1234I(13,4)+3*(D1234I(10,4)-D1234I(11,4)-D1234I(14,4)
     -   +D1234I(15,4)))*P(8)
       F(48)=DCMPLX(FR(48),FI(48))
       FR(49) = 16*(D1234R(2,3)+D1234R(2,4)-D1234R(3,3)+D1234R(3,4)-3*
     -   (D1234R(8,3)-D1234R(9,3))+6*D1234R(14,4)-4*(D1234R(13,4)+D12
     -   34R(15,4)))*P(8)
       FI(49) = 16*(D1234I(2,3)+D1234I(2,4)-D1234I(3,3)+D1234I(3,4)-3*
     -   (D1234I(8,3)-D1234I(9,3))+6*D1234I(14,4)-4*(D1234I(13,4)+D12
     -   34I(15,4)))*P(8)
       F(49)=DCMPLX(FR(49),FI(49))
       FR(50) = -8*(D1234R(2,2)+D1234R(2,3)+D1234R(3,2)+D1234R(3,3)-D1
     -   234R(8,3)-D1234R(9,3)-2*(D1234R(3,4)+D1234R(6,2)-D1234R(13,4
     -   ))-6*(D1234R(14,4)-D1234R(15,4)))*P(8)
       FI(50) = -8*(D1234I(2,2)+D1234I(2,3)+D1234I(3,2)+D1234I(3,3)-D1
     -   234I(8,3)-D1234I(9,3)-2*(D1234I(3,4)+D1234I(6,2)-D1234I(13,4
     -   ))-6*(D1234I(14,4)-D1234I(15,4)))*P(8)
       F(50)=DCMPLX(FR(50),FI(50))
       P(297) = -2*M**2*P(19)*P(86)*P(127)-P(8)*P(243)
       P(298) = s12*za1*za2*za3*za4+s23*za1*za2*za3*za4+s12*zb1*zb2*zb
     -   3*zb4+s23*zb1*zb2*zb3*zb4-p1sq*P(8)-p3sq*P(8)-2*P(50)
       P(299) = p1sq-s23
       P(300) = p3sq-p4sq-s12+2*P(299)
       P(301) = 2*M**2*P(19)*P(167)+P(8)*P(300)
       FR(51) = (C0123R+C0234R)*P(41)-2*((C123R(1,1)-2*(D1234R(7,2)+2*
     -   (D1234R(11,3)+D1234R(12,3)-2*(D1234R(13,3)-D1234R(18,4)-D123
     -   4R(19,4)+D1234R(20,4)+D1234R(21,4)))))*P(8)-(-D1234R(3,2)-D1
     -   234R(4,2)+D1234R(5,2)+D1234R(6,2))*P(270))+D01234R*P(273)+D1
     -   234R(1,1)*P(297)+D1234R(2,1)*P(298)+D1234R(3,1)*P(301)
       FI(51) = (C0123I+C0234I)*P(41)-2*((C123I(1,1)-2*(D1234I(7,2)+2*
     -   (D1234I(11,3)+D1234I(12,3)-2*(D1234I(13,3)-D1234I(18,4)-D123
     -   4I(19,4)+D1234I(20,4)+D1234I(21,4)))))*P(8)-(-D1234I(3,2)-D1
     -   234I(4,2)+D1234I(5,2)+D1234I(6,2))*P(270))+D01234I*P(273)+D1
     -   234I(1,1)*P(297)+D1234I(2,1)*P(298)+D1234I(3,1)*P(301)
       F(51)=DCMPLX(FR(51),FI(51))
       P(302) = 2*M**2*P(19)*P(86)*P(127)+P(8)*P(243)
       P(303) = M**2*P(1)*P(25)*P(157)
       P(304) = p3sq*za1*za2*za3*za4+s12*za1*za2*za3*za4-s23*za1*za2*z
     -   a3*za4+p3sq*zb1*zb2*zb3*zb4+s12*zb1*zb2*zb3*zb4-s23*zb1*zb2*
     -   zb3*zb4-p1sq*P(8)-2*P(303)
       P(305) = p3sq-p4sq+s12-2*s23
       P(306) = za4*zb2-2*za2*zb4+zb2*zb4
       P(307) = zb2*zb3*P(111)+za3*P(306)
       P(308) = -2*za4*zb2*zb3+za2*P(289)
       P(309) = zb1*P(307)+za1*P(308)
       P(310) = -(P(8)*P(305))-2*M**2*P(309)
       FR(52) = (C0123R+C0234R+2*(C123R(1,1)+C134R(1,1)-C134R(2,1)+4*(
     -   D1234R(12,3)-D1234R(13,3))+8*(D1234R(18,4)+D1234R(19,4)-D123
     -   4R(20,4)-D1234R(21,4))))*P(8)+2*((D1234R(3,2)+D1234R(4,2))*P
     -   (132)+(D1234R(5,2)+D1234R(6,2))*P(135))+D01234R*P(284)+D1234
     -   R(1,1)*P(302)+D1234R(2,1)*P(304)+D1234R(3,1)*P(310)
       FI(52) = (C0123I+C0234I+2*(C123I(1,1)+C134I(1,1)-C134I(2,1)+4*(
     -   D1234I(12,3)-D1234I(13,3))+8*(D1234I(18,4)+D1234I(19,4)-D123
     -   4I(20,4)-D1234I(21,4))))*P(8)+2*((D1234I(3,2)+D1234I(4,2))*P
     -   (132)+(D1234I(5,2)+D1234I(6,2))*P(135))+D01234I*P(284)+D1234
     -   I(1,1)*P(302)+D1234I(2,1)*P(304)+D1234I(3,1)*P(310)
       F(52)=DCMPLX(FR(52),FI(52))
       P(311) = -(p1sq*P(8))+s23*P(8)-2*M**2*P(19)*P(119)
       P(312) = p1sq*P(8)-s23*P(8)+2*M**2*P(19)*P(119)
       FR(53) = 2*(C234R(1,1)*P(8)-(C123R(2,1)+C234R(2,1)-4*(D1234R(12
     -   ,3)-D1234R(13,3)+2*(D1234R(17,4)+D1234R(18,4)-2*D1234R(21,4)
     -   )))*P(8)-(D1234R(2,2)+D1234R(3,2)-2*D1234R(6,2))*P(270)+D123
     -   4R(2,1)*P(311)+D1234R(3,1)*P(312))
       FI(53) = 2*(C234I(1,1)*P(8)-(C123I(2,1)+C234I(2,1)-4*(D1234I(12
     -   ,3)-D1234I(13,3)+2*(D1234I(17,4)+D1234I(18,4)-2*D1234I(21,4)
     -   )))*P(8)-(D1234I(2,2)+D1234I(3,2)-2*D1234I(6,2))*P(270)+D123
     -   4I(2,1)*P(311)+D1234I(3,1)*P(312))
       F(53)=DCMPLX(FR(53),FI(53))
       P(313) = s23*P(8)-2*M**2*P(61)*P(127)
       FR(54) = 2*((C123R(2,1)+C134R(1,1)-C134R(2,1)+8*(D1234R(17,4)+D
     -   1234R(18,4)-2*D1234R(21,4)))*P(8)+(D1234R(2,2)+D1234R(3,2))*
     -   P(132)+2*D1234R(6,2)*P(135)+(-D1234R(2,1)+D1234R(3,1))*P(313
     -   ))
       FI(54) = 2*((C123I(2,1)+C134I(1,1)-C134I(2,1)+8*(D1234I(17,4)+D
     -   1234I(18,4)-2*D1234I(21,4)))*P(8)+(D1234I(2,2)+D1234I(3,2))*
     -   P(132)+2*D1234I(6,2)*P(135)+(-D1234I(2,1)+D1234I(3,1))*P(313
     -   ))
       F(54)=DCMPLX(FR(54),FI(54))
       P(314) = p1sq+p3sq+2*p4sq-s12-s23
       P(315) = 2*M**2*P(19)*P(187)-P(8)*P(314)
       P(316) = p3sq+3*p4sq-s12-2*s23
       P(317) = -3*za4+zb4
       P(318) = za4*zb2+2*za2*zb4-zb2*zb4
       P(319) = zb2*zb3*P(317)+za3*P(318)
       P(320) = za3*za4-za4*zb3-3*za3*zb4+zb3*zb4
       P(321) = 2*za4*zb2*zb3+za2*P(320)
       P(322) = zb1*P(319)+za1*P(321)
       P(323) = P(8)*P(316)-2*M**2*P(322)
       FR(55) = (C0123R+C0234R-2*(-C123R(2,1)+C134R(1,1)-C134R(2,1)-C2
     -   34R(1,1)+C234R(2,1)-4*(D1234R(11,3)-D1234R(13,3))-8*(D1234R(
     -   18,4)+D1234R(19,4)-D1234R(20,4)-D1234R(21,4))))*P(8)+D01234R
     -   *P(284)+D1234R(1,1)*P(285)+2*(-D1234R(3,2)-D1234R(4,2)+D1234
     -   R(5,2)+D1234R(6,2))*P(295)+D1234R(2,1)*P(315)+D1234R(3,1)*P(
     -   323)
       FI(55) = (C0123I+C0234I-2*(-C123I(2,1)+C134I(1,1)-C134I(2,1)-C2
     -   34I(1,1)+C234I(2,1)-4*(D1234I(11,3)-D1234I(13,3))-8*(D1234I(
     -   18,4)+D1234I(19,4)-D1234I(20,4)-D1234I(21,4))))*P(8)+D01234I
     -   *P(284)+D1234I(1,1)*P(285)+2*(-D1234I(3,2)-D1234I(4,2)+D1234
     -   I(5,2)+D1234I(6,2))*P(295)+D1234I(2,1)*P(315)+D1234I(3,1)*P(
     -   323)
       F(55)=DCMPLX(FR(55),FI(55))
       P(324) = p3sq+p4sq-s12+2*P(299)
       P(325) = 2*M**2*P(19)*P(118)+P(8)*P(324)
       FR(56) = (C0123R-2*(C234R(2,1)+2*D1234R(7,2)-4*(D1234R(12,3)-D1
     -   234R(13,3))-8*(D1234R(18,4)+D1234R(19,4)-D1234R(20,4)-D1234R
     -   (21,4))))*P(8)+C0234R*P(41)+D1234R(1,1)*P(108)+2*(-D1234R(3,
     -   2)-D1234R(4,2)+D1234R(5,2)+D1234R(6,2))*P(270)+D01234R*P(284
     -   )+D1234R(2,1)*P(298)+D1234R(3,1)*P(325)
       FI(56) = (C0123I-2*(C234I(2,1)+2*D1234I(7,2)-4*(D1234I(12,3)-D1
     -   234I(13,3))-8*(D1234I(18,4)+D1234I(19,4)-D1234I(20,4)-D1234I
     -   (21,4))))*P(8)+C0234I*P(41)+D1234I(1,1)*P(108)+2*(-D1234I(3,
     -   2)-D1234I(4,2)+D1234I(5,2)+D1234I(6,2))*P(270)+D01234I*P(284
     -   )+D1234I(2,1)*P(298)+D1234I(3,1)*P(325)
       F(56)=DCMPLX(FR(56),FI(56))
       P(326) = p1sq+p3sq-2*p4sq+s12+s23
       P(327) = -2*M**2*P(25)*P(150)+P(8)*P(326)
       P(328) = p3sq-3*p4sq+s12+2*s23
       P(329) = 3*za4+zb4
       P(330) = za4*zb2+2*za2*zb4+zb2*zb4
       P(331) = -(zb2*zb3*P(329))+za3*P(330)
       P(332) = za4+3*zb4
       P(333) = -(zb3*P(25))+za3*P(332)
       P(334) = -2*za4*zb2*zb3+za2*P(333)
       P(335) = -(zb1*P(331))+za1*P(334)
       P(336) = -(P(8)*P(328))+2*M**2*P(335)
       FR(57) = (C0123R-2*(C134R(1,1)-C134R(2,1)+C234R(2,1)-8*(D1234R(
     -   12,3)-D1234R(13,3)+D1234R(18,4)+D1234R(19,4)-D1234R(20,4)-D1
     -   234R(21,4))))*P(8)+C0234R*P(41)+D1234R(1,1)*P(108)+D01234R*P
     -   (284)+2*(-D1234R(3,2)-D1234R(4,2)+D1234R(5,2)+D1234R(6,2))*P
     -   (295)+D1234R(2,1)*P(327)+D1234R(3,1)*P(336)
       FI(57) = (C0123I-2*(C134I(1,1)-C134I(2,1)+C234I(2,1)-8*(D1234I(
     -   12,3)-D1234I(13,3)+D1234I(18,4)+D1234I(19,4)-D1234I(20,4)-D1
     -   234I(21,4))))*P(8)+C0234I*P(41)+D1234I(1,1)*P(108)+D01234I*P
     -   (284)+2*(-D1234I(3,2)-D1234I(4,2)+D1234I(5,2)+D1234I(6,2))*P
     -   (295)+D1234I(2,1)*P(327)+D1234I(3,1)*P(336)
       F(57)=DCMPLX(FR(57),FI(57))
       P(337) = za4*zb1*zb2*zb3+za1*za2*za3*zb4
       P(338) = M**2*P(337)
       P(339) = M**2*za3*za4*zb1*zb2+M**2*za1*za4*zb2*zb3+M**2*za2*za3
     -   *zb1*zb4+M**2*za1*za2*zb3*zb4-p4sq*P(8)
       P(340) = s12*za1*za2*za3*za4+s23*za1*za2*za3*za4+s12*zb1*zb2*zb
     -   3*zb4+s23*zb1*zb2*zb3*zb4-4*P(338)+2*P(339)
       P(341) = za1*za4*zb2-2*za4*zb1*zb2+za1*za2*zb4
       P(342) = -2*za1+zb1
       P(343) = za4*zb1*zb2+za2*zb4*P(342)
       P(344) = zb3*P(341)+za3*P(343)
       P(345) = s12*P(8)+s23*P(8)+2*M**2*P(344)
       P(346) = 2*p4sq*P(8)-P(345)
       FR(58) = 2*((C0123R+C123R(1,1)-C124R(1,1)+C124R(2,1)-C134R(1,1)
     -   +C134R(2,1)-C234R(2,1)+8*(D1234R(11,3)-D1234R(13,3)+D1234R(1
     -   6,4)+D1234R(18,4)-2*D1234R(20,4)))*P(8)-D01234R*P(273)-(D123
     -   4R(1,2)+D1234R(3,2)-2*D1234R(5,2))*P(295)+D1234R(1,1)*P(340)
     -   +D1234R(3,1)*P(346))
       FI(58) = 2*((C0123I+C123I(1,1)-C124I(1,1)+C124I(2,1)-C134I(1,1)
     -   +C134I(2,1)-C234I(2,1)+8*(D1234I(11,3)-D1234I(13,3)+D1234I(1
     -   6,4)+D1234I(18,4)-2*D1234I(20,4)))*P(8)-D01234I*P(273)-(D123
     -   4I(1,2)+D1234I(3,2)-2*D1234I(5,2))*P(295)+D1234I(1,1)*P(340)
     -   +D1234I(3,1)*P(346))
       F(58)=DCMPLX(FR(58),FI(58))
       P(347) = p1sq*za1*za2*za3*za4-p3sq*za1*za2*za3*za4-s12*za1*za2*
     -   za3*za4+s23*za1*za2*za3*za4+p1sq*zb1*zb2*zb3*zb4-p3sq*zb1*zb
     -   2*zb3*zb4-s12*zb1*zb2*zb3*zb4+s23*zb1*zb2*zb3*zb4+2*P(303)
       P(348) = 2*p1sq-p3sq+p4sq-s12
       P(349) = -2*za2*za3*za4+zb2*P(287)
       P(350) = -2*zb2*zb3*zb4+za2*P(289)
       P(351) = zb1*P(349)+za1*P(350)
       P(352) = -(P(8)*P(348))+2*M**2*P(351)
       FR(59) = (C0234R-2*(C134R(1,1)-C134R(2,1)-C234R(2,1)+2*(C134R(1
     -   ,2)+C134R(2,2)-2*(C134R(3,2)+D1234R(12,3)-D1234R(13,3)+2*(D1
     -   234R(18,4)+D1234R(19,4)-D1234R(20,4)-D1234R(21,4))))))*P(8)+
     -   C0123R*P(41)+D1234R(1,1)*P(109)+2*((D1234R(3,2)+D1234R(4,2))
     -   *P(37)+(D1234R(5,2)+D1234R(6,2))*P(146))+D01234R*P(273)+D123
     -   4R(2,1)*P(347)+D1234R(3,1)*P(352)
       FI(59) = (C0234I-2*(C134I(1,1)-C134I(2,1)-C234I(2,1)+2*(C134I(1
     -   ,2)+C134I(2,2)-2*(C134I(3,2)+D1234I(12,3)-D1234I(13,3)+2*(D1
     -   234I(18,4)+D1234I(19,4)-D1234I(20,4)-D1234I(21,4))))))*P(8)+
     -   C0123I*P(41)+D1234I(1,1)*P(109)+2*((D1234I(3,2)+D1234I(4,2))
     -   *P(37)+(D1234I(5,2)+D1234I(6,2))*P(146))+D01234I*P(273)+D123
     -   4I(2,1)*P(347)+D1234I(3,1)*P(352)
       F(59)=DCMPLX(FR(59),FI(59))
       P(353) = p1sq-p3sq+2*p4sq+s12-s23
       P(354) = -2*M**2*P(19)*P(179)+P(8)*P(353)
       P(355) = -p3sq+3*p4sq+s12+2*P(299)
       P(356) = 2*M**2*P(19)*P(184)-P(8)*P(355)
       FR(60) = (C0234R+2*(C234R(2,1)+2*(-C234R(2,2)+C234R(3,2)+D1234R
     -   (7,2)+2*D1234R(11,3)-6*D1234R(13,3)+4*(D1234R(12,3)+D1234R(1
     -   8,4)+D1234R(19,4)-D1234R(20,4)-D1234R(21,4)))))*P(8)+C0123R*
     -   P(41)+2*(D1234R(5,2)+D1234R(6,2))*P(108)+(D1234R(1,1)+2*(D12
     -   34R(3,2)+D1234R(4,2)))*P(109)+D01234R*P(273)+D1234R(2,1)*P(3
     -   54)+D1234R(3,1)*P(356)
       FI(60) = (C0234I+2*(C234I(2,1)+2*(-C234I(2,2)+C234I(3,2)+D1234I
     -   (7,2)+2*D1234I(11,3)-6*D1234I(13,3)+4*(D1234I(12,3)+D1234I(1
     -   8,4)+D1234I(19,4)-D1234I(20,4)-D1234I(21,4)))))*P(8)+C0123I*
     -   P(41)+2*(D1234I(5,2)+D1234I(6,2))*P(108)+(D1234I(1,1)+2*(D12
     -   34I(3,2)+D1234I(4,2)))*P(109)+D01234I*P(273)+D1234I(2,1)*P(3
     -   54)+D1234I(3,1)*P(356)
       F(60)=DCMPLX(FR(60),FI(60))
       P(357) = 2*M**2*P(19)*P(171)+P(8)*P(314)
       P(358) = p3sq+3*p4sq-s12+2*P(299)
       P(359) = 2*M**2*P(19)*P(174)-P(8)*P(358)
       FR(61) = (C0234R+2*(-C123R(2,1)+C234R(1,1)+2*(-C234R(2,2)+C234R
     -   (3,2)+D1234R(7,2)+2*(D1234R(11,3)+D1234R(12,3)-2*(D1234R(13,
     -   3)-D1234R(18,4)-D1234R(19,4)+D1234R(20,4)+D1234R(21,4))))))*
     -   P(8)+C0123R*P(41)+2*(D1234R(5,2)+D1234R(6,2))*P(108)+(D1234R
     -   (1,1)+2*(D1234R(3,2)+D1234R(4,2)))*P(109)+D01234R*P(273)+D12
     -   34R(2,1)*P(357)+D1234R(3,1)*P(359)
       FI(61) = (C0234I+2*(-C123I(2,1)+C234I(1,1)+2*(-C234I(2,2)+C234I
     -   (3,2)+D1234I(7,2)+2*(D1234I(11,3)+D1234I(12,3)-2*(D1234I(13,
     -   3)-D1234I(18,4)-D1234I(19,4)+D1234I(20,4)+D1234I(21,4))))))*
     -   P(8)+C0123I*P(41)+2*(D1234I(5,2)+D1234I(6,2))*P(108)+(D1234I
     -   (1,1)+2*(D1234I(3,2)+D1234I(4,2)))*P(109)+D01234I*P(273)+D12
     -   34I(2,1)*P(357)+D1234I(3,1)*P(359)
       F(61)=DCMPLX(FR(61),FI(61))
       P(360) = P(8)*P(65)-2*M**2*P(9)*P(60)*P(207)
       P(361) = p2sq-p4sq-s12-s23
       P(362) = za3*za4*zb2+za4*zb2*zb3+za2*za3*zb4-za2*zb3*zb4-2*zb2*
     -   zb3*zb4
       P(363) = 2*za3*za4-za3*zb4-zb3*zb4
       P(364) = za4*zb2*P(2)+za2*P(363)
       P(365) = -(zb1*P(362))+za1*P(364)
       P(366) = P(8)*P(361)+2*M**2*P(365)
       FR(62) = (C0134R+C0234R)*P(41)-2*((C234R(1,1)-C234R(2,1)+4*D123
     -   4R(13,3)-8*(D1234R(18,4)-D1234R(21,4)))*P(8)-D1234R(3,2)*P(2
     -   10)-D1234R(6,2)*P(213))+D01234R*P(296)+D1234R(2,1)*P(360)+D1
     -   234R(3,1)*P(366)
       FI(62) = (C0134I+C0234I)*P(41)-2*((C234I(1,1)-C234I(2,1)+4*D123
     -   4I(13,3)-8*(D1234I(18,4)-D1234I(21,4)))*P(8)-D1234I(3,2)*P(2
     -   10)-D1234I(6,2)*P(213))+D01234I*P(296)+D1234I(2,1)*P(360)+D1
     -   234I(3,1)*P(366)
       F(62)=DCMPLX(FR(62),FI(62))
       P(367) = p2sq+p4sq-s12+s23
       P(368) = za4*zb2+za2*zb4-2*zb2*zb4
       P(369) = za3*P(87)-zb3*P(368)
       P(370) = za2*P(30)-za4*zb2*P(33)
       P(371) = zb1*P(369)+za1*P(370)
       P(372) = -(P(8)*P(367))+2*M**2*P(371)
       FR(63) = D1234R(2,1)*P(37)+(C0134R+C0234R)*P(41)-2*((C134R(1,1)
     -   -C134R(2,1)+4*D1234R(13,3)-8*(D1234R(18,4)-D1234R(21,4)))*P(
     -   8)-(-D1234R(3,2)+D1234R(6,2))*P(270))+D01234R*P(296)+D1234R(
     -   3,1)*P(372)
       FI(63) = D1234I(2,1)*P(37)+(C0134I+C0234I)*P(41)-2*((C134I(1,1)
     -   -C134I(2,1)+4*D1234I(13,3)-8*(D1234I(18,4)-D1234I(21,4)))*P(
     -   8)-(-D1234I(3,2)+D1234I(6,2))*P(270))+D01234I*P(296)+D1234I(
     -   3,1)*P(372)
       F(63)=DCMPLX(FR(63),FI(63))
       P(373) = -(p1sq*P(8))+2*M**2*P(9)*P(206)
       FR(64) = (C0134R+C0234R+2*(C134R(1,1)-2*(C134R(2,2)-C134R(3,2)+
     -   D1234R(7,2)+2*D1234R(12,3)-4*(D1234R(18,4)-D1234R(21,4)))))*
     -   P(8)+2*D1234R(3,2)*P(37)+(D1234R(2,1)+2*D1234R(6,2))*P(146)+
     -   D1234R(3,1)*P(295)+D01234R*P(373)
       FI(64) = (C0134I+C0234I+2*(C134I(1,1)-2*(C134I(2,2)-C134I(3,2)+
     -   D1234I(7,2)+2*D1234I(12,3)-4*(D1234I(18,4)-D1234I(21,4)))))*
     -   P(8)+2*D1234I(3,2)*P(37)+(D1234I(2,1)+2*D1234I(6,2))*P(146)+
     -   D1234I(3,1)*P(295)+D01234I*P(373)
       F(64)=DCMPLX(FR(64),FI(64))
       P(374) = -(P(8)*P(65))+2*M**2*P(9)*P(60)*P(207)
       P(375) = M**2*P(9)*P(33)*P(87)
       P(376) = p4sq*za1*za2*za3*za4+s12*za1*za2*za3*za4-s23*za1*za2*z
     -   a3*za4+p4sq*zb1*zb2*zb3*zb4+s12*zb1*zb2*zb3*zb4-s23*zb1*zb2*
     -   zb3*zb4-p2sq*P(8)-2*P(375)
       FR(65) = (C0134R+C0234R+2*(C134R(2,1)+C234R(1,1)-C234R(2,1)-2*D
     -   1234R(7,2)-4*D1234R(12,3)+8*(D1234R(18,4)-D1234R(21,4))))*P(
     -   8)+2*(-D1234R(3,2)+D1234R(6,2))*P(295)+D01234R*P(373)+D1234R
     -   (2,1)*P(374)+D1234R(3,1)*P(376)
       FI(65) = (C0134I+C0234I+2*(C134I(2,1)+C234I(1,1)-C234I(2,1)-2*D
     -   1234I(7,2)-4*D1234I(12,3)+8*(D1234I(18,4)-D1234I(21,4))))*P(
     -   8)+2*(-D1234I(3,2)+D1234I(6,2))*P(295)+D01234I*P(373)+D1234I
     -   (2,1)*P(374)+D1234I(3,1)*P(376)
       F(65)=DCMPLX(FR(65),FI(65))
       P(377) = P(8)*P(133)+2*M**2*P(1)*P(86)*P(207)
       P(378) = p2sq+p4sq+s12-s23
       P(379) = za3*za4*zb2+za4*zb2*zb3-za2*za3*zb4+za2*zb3*zb4-2*zb2*
     -   zb3*zb4
       P(380) = -2*za3*za4+za3*zb4+zb3*zb4
       P(381) = za4*zb2*P(2)+za2*P(380)
       P(382) = zb1*P(379)+za1*P(381)
       P(383) = P(8)*P(378)+2*M**2*P(382)
       FR(66) = (C0134R+C0234R)*P(41)-2*((C234R(2,1)+4*(D1234R(12,3)-D
     -   1234R(13,3))-8*(D1234R(18,4)-D1234R(21,4)))*P(8)-(-D1234R(3,
     -   2)+D1234R(6,2))*P(270))+D01234R*P(296)+D1234R(2,1)*P(377)+D1
     -   234R(3,1)*P(383)
       FI(66) = (C0134I+C0234I)*P(41)-2*((C234I(2,1)+4*(D1234I(12,3)-D
     -   1234I(13,3))-8*(D1234I(18,4)-D1234I(21,4)))*P(8)-(-D1234I(3,
     -   2)+D1234I(6,2))*P(270))+D01234I*P(296)+D1234I(2,1)*P(377)+D1
     -   234I(3,1)*P(383)
       F(66)=DCMPLX(FR(66),FI(66))
       FR(67) = D1234R(2,1)*P(37)+(C0134R+C0234R)*P(41)-2*((C134R(1,1)
     -   -C134R(2,1)+C234R(2,1)-2*D1234R(7,2)-8*(D1234R(18,4)-D1234R(
     -   21,4)))*P(8)-(-D1234R(3,2)+D1234R(6,2))*P(295))+D01234R*P(29
     -   6)+D1234R(3,1)*P(376)
       FI(67) = D1234I(2,1)*P(37)+(C0134I+C0234I)*P(41)-2*((C134I(1,1)
     -   -C134I(2,1)+C234I(2,1)-2*D1234I(7,2)-8*(D1234I(18,4)-D1234I(
     -   21,4)))*P(8)-(-D1234I(3,2)+D1234I(6,2))*P(295))+D01234I*P(29
     -   6)+D1234I(3,1)*P(376)
       F(67)=DCMPLX(FR(67),FI(67))
       FR(68) = (C0124R+2*(C124R(1,1)-C124R(2,1)-2*D1234R(7,2)-4*D1234
     -   R(11,3)+8*(D1234R(18,4)-D1234R(20,4))))*P(8)+C0134R*P(41)+D1
     -   234R(3,1)*P(132)+2*(-D1234R(3,2)+D1234R(5,2))*P(270)+D01234R
     -   *P(282)+D1234R(1,1)*P(377)
       FI(68) = (C0124I+2*(C124I(1,1)-C124I(2,1)-2*D1234I(7,2)-4*D1234
     -   I(11,3)+8*(D1234I(18,4)-D1234I(20,4))))*P(8)+C0134I*P(41)+D1
     -   234I(3,1)*P(132)+2*(-D1234I(3,2)+D1234I(5,2))*P(270)+D01234I
     -   *P(282)+D1234I(1,1)*P(377)
       F(68)=DCMPLX(FR(68),FI(68))
       P(384) = 2*M**2*P(1)*P(33)*P(61)-P(8)*P(144)
       FR(69) = (C0124R-2*(C134R(1,1)-C134R(2,1)+2*D1234R(7,2)+4*D1234
     -   R(13,3)-8*(D1234R(18,4)-D1234R(20,4))))*P(8)+D1234R(1,1)*P(3
     -   7)+C0134R*P(41)+2*(-D1234R(3,2)+D1234R(5,2))*P(270)+D01234R*
     -   P(282)+D1234R(3,1)*P(384)
       FI(69) = (C0124I-2*(C134I(1,1)-C134I(2,1)+2*D1234I(7,2)+4*D1234
     -   I(13,3)-8*(D1234I(18,4)-D1234I(20,4))))*P(8)+D1234I(1,1)*P(3
     -   7)+C0134I*P(41)+2*(-D1234I(3,2)+D1234I(5,2))*P(270)+D01234I*
     -   P(282)+D1234I(3,1)*P(384)
       F(69)=DCMPLX(FR(69),FI(69))
       P(385) = za2*za3*za4+za3*za4*zb2-2*zb2*zb3*zb4
       P(386) = -2*za2*za3*za4+za2*zb3*zb4+zb2*zb3*zb4
       P(387) = zb1*P(385)+za1*P(386)
       P(388) = -(P(8)*P(100))-2*M**2*P(387)
       P(389) = p2sq+p3sq-2*p4sq+s23
       P(390) = zb2*P(2)+za2*P(33)
       P(391) = 2*za4*zb2*zb3-zb4*P(390)
       P(392) = za3*za4-za4*zb3+2*za3*zb4
       P(393) = -(za4*zb2*P(33))+za2*P(392)
       P(394) = zb1*P(391)+za1*P(393)
       P(395) = -(P(8)*P(389))+2*M**2*P(394)
       FR(70) = (C0124R+C0134R+2*(C124R(2,1)+C134R(2,1)-C234R(2,1)-4*(
     -   D1234R(7,2)+D1234R(11,3)+D1234R(13,3)-2*(D1234R(18,4)-D1234R
     -   (20,4)))))*P(8)+2*(-D1234R(3,2)+D1234R(5,2))*P(295)+D1234R(1
     -   ,1)*P(374)+D01234R*P(388)+D1234R(3,1)*P(395)
       FI(70) = (C0124I+C0134I+2*(C124I(2,1)+C134I(2,1)-C234I(2,1)-4*(
     -   D1234I(7,2)+D1234I(11,3)+D1234I(13,3)-2*(D1234I(18,4)-D1234I
     -   (20,4)))))*P(8)+2*(-D1234I(3,2)+D1234I(5,2))*P(295)+D1234I(1
     -   ,1)*P(374)+D01234I*P(388)+D1234I(3,1)*P(395)
       F(70)=DCMPLX(FR(70),FI(70))
       P(396) = p2sq-p3sq-2*p4sq+s23
       P(397) = za2*P(33)+zb2*P(115)
       P(398) = -2*za4*zb2*zb3+zb4*P(397)
       P(399) = za3*za4-za4*zb3-2*za3*zb4
       P(400) = za4*zb2*P(33)+za2*P(399)
       P(401) = zb1*P(398)+za1*P(400)
       P(402) = -(P(8)*P(396))-2*M**2*P(401)
       FR(71) = (C0124R-2*(C134R(1,1)-C134R(2,1)+C234R(2,1)-8*(D1234R(
     -   18,4)-D1234R(20,4))))*P(8)+D1234R(1,1)*P(37)+C0134R*P(41)+D0
     -   1234R*P(282)+2*(-D1234R(3,2)+D1234R(5,2))*P(295)+D1234R(3,1)
     -   *P(402)
       FI(71) = (C0124I-2*(C134I(1,1)-C134I(2,1)+C234I(2,1)-8*(D1234I(
     -   18,4)-D1234I(20,4))))*P(8)+D1234I(1,1)*P(37)+C0134I*P(41)+D0
     -   1234I*P(282)+2*(-D1234I(3,2)+D1234I(5,2))*P(295)+D1234I(3,1)
     -   *P(402)
       F(71)=DCMPLX(FR(71),FI(71))
       FR(72) = (C0134R+2*(-C124R(2,1)+C134R(1,1)-2*(C134R(2,2)-C134R(
     -   3,2)+D1234R(7,2)+2*D1234R(11,3)-4*(D1234R(18,4)-D1234R(20,4)
     -   ))))*P(8)+2*D1234R(3,2)*P(37)+C0124R*P(41)+D1234R(3,1)*P(132
     -   )+(D1234R(1,1)+2*D1234R(5,2))*P(146)+D01234R*P(283)
       FI(72) = (C0134I+2*(-C124I(2,1)+C134I(1,1)-2*(C134I(2,2)-C134I(
     -   3,2)+D1234I(7,2)+2*D1234I(11,3)-4*(D1234I(18,4)-D1234I(20,4)
     -   ))))*P(8)+2*D1234I(3,2)*P(37)+C0124I*P(41)+D1234I(3,1)*P(132
     -   )+(D1234I(1,1)+2*D1234I(5,2))*P(146)+D01234I*P(283)
       F(72)=DCMPLX(FR(72),FI(72))
       P(403) = P(8)*P(100)+2*M**2*P(387)
       P(404) = p2sq-p3sq-2*p4sq-s23
       P(405) = za2*za3+za3*zb2-za2*zb3-3*zb2*zb3
       P(406) = 2*za4*zb2*zb3+zb4*P(405)
       P(407) = 3*za3*za4-za4*zb3-2*za3*zb4
       P(408) = za4*zb2*P(2)+za2*P(407)
       P(409) = -(zb1*P(406))+za1*P(408)
       P(410) = P(8)*P(404)+2*M**2*P(409)
       FR(73) = (C0124R+C0134R)*P(41)+2*((C234R(2,1)-8*(D1234R(13,3)-D
     -   1234R(18,4)+D1234R(20,4)))*P(8)+D1234R(3,2)*P(210)+D1234R(5,
     -   2)*P(213))+D1234R(1,1)*P(360)+D01234R*P(403)+D1234R(3,1)*P(4
     -   10)
       FI(73) = (C0124I+C0134I)*P(41)+2*((C234I(2,1)-8*(D1234I(13,3)-D
     -   1234I(18,4)+D1234I(20,4)))*P(8)+D1234I(3,2)*P(210)+D1234I(5,
     -   2)*P(213))+D1234I(1,1)*P(360)+D01234I*P(403)+D1234I(3,1)*P(4
     -   10)
       F(73)=DCMPLX(FR(73),FI(73))
       P(411) = p2sq*za1*za2*za3*za4-p4sq*za1*za2*za3*za4-s12*za1*za2*
     -   za3*za4+s23*za1*za2*za3*za4+p2sq*zb1*zb2*zb3*zb4-p4sq*zb1*zb
     -   2*zb3*zb4-s12*zb1*zb2*zb3*zb4+s23*zb1*zb2*zb3*zb4+2*P(375)
       FR(74) = (C0134R+C0234R+2*(C134R(1,1)-C134R(2,1)+C234R(2,1)-2*(
     -   C134R(2,2)-C134R(3,2)+D1234R(7,2)+2*(D1234R(12,3)-D1234R(13,
     -   3))-4*(D1234R(18,4)-D1234R(21,4)))))*P(8)+2*D1234R(3,2)*P(37
     -   )+(D1234R(2,1)+2*D1234R(6,2))*P(146)+D01234R*P(373)+D1234R(3
     -   ,1)*P(411)
       FI(74) = (C0134I+C0234I+2*(C134I(1,1)-C134I(2,1)+C234I(2,1)-2*(
     -   C134I(2,2)-C134I(3,2)+D1234I(7,2)+2*(D1234I(12,3)-D1234I(13,
     -   3))-4*(D1234I(18,4)-D1234I(21,4)))))*P(8)+2*D1234I(3,2)*P(37
     -   )+(D1234I(2,1)+2*D1234I(6,2))*P(146)+D01234I*P(373)+D1234I(3
     -   ,1)*P(411)
       F(74)=DCMPLX(FR(74),FI(74))
       P(412) = -(P(8)*P(133))-2*M**2*P(1)*P(86)*P(207)
       P(413) = -(zb1*P(379))-za1*P(381)
       P(414) = -(P(8)*P(378))+2*M**2*P(413)
       FR(75) = (C0134R+C0234R+2*(C234R(2,1)+8*(D1234R(18,4)-D1234R(21
     -   ,4))))*P(8)+2*(D1234R(3,2)*P(210)+D1234R(6,2)*P(213))+D01234
     -   R*P(373)+D1234R(2,1)*P(412)+D1234R(3,1)*P(414)
       FI(75) = (C0134I+C0234I+2*(C234I(2,1)+8*(D1234I(18,4)-D1234I(21
     -   ,4))))*P(8)+2*(D1234I(3,2)*P(210)+D1234I(6,2)*P(213))+D01234
     -   I*P(373)+D1234I(2,1)*P(412)+D1234I(3,1)*P(414)
       F(75)=DCMPLX(FR(75),FI(75))
       P(415) = s12*za1*za2*za3*za4+s23*za1*za2*za3*za4+s12*zb1*zb2*zb
     -   3*zb4+s23*zb1*zb2*zb3*zb4-p2sq*P(8)-p4sq*P(8)-2*P(139)
       FR(76) = (C0134R+C0234R+2*(C234R(1,1)-2*(C234R(2,2)-C234R(3,2)+
     -   D1234R(7,2)+2*D1234R(12,3)-4*(D1234R(18,4)-D1234R(21,4)))))*
     -   P(8)+2*(D1234R(6,2)*P(108)+D1234R(3,2)*P(109))+D01234R*P(373
     -   )+D1234R(2,1)*P(374)+D1234R(3,1)*P(415)
       FI(76) = (C0134I+C0234I+2*(C234I(1,1)-2*(C234I(2,2)-C234I(3,2)+
     -   D1234I(7,2)+2*D1234I(12,3)-4*(D1234I(18,4)-D1234I(21,4)))))*
     -   P(8)+2*(D1234I(6,2)*P(108)+D1234I(3,2)*P(109))+D01234I*P(373
     -   )+D1234I(2,1)*P(374)+D1234I(3,1)*P(415)
       F(76)=DCMPLX(FR(76),FI(76))
       FR(77) = (C0134R+C0234R+2*(C134R(1,1)-C134R(2,1)+C234R(2,1)-2*(
     -   C234R(2,2)-C234R(3,2)+D1234R(7,2)+2*D1234R(13,3)-4*(D1234R(1
     -   8,4)-D1234R(21,4)))))*P(8)+2*(D1234R(6,2)*P(108)+D1234R(3,2)
     -   *P(109))+D1234R(2,1)*P(146)+D01234R*P(373)+D1234R(3,1)*P(411
     -   )
       FI(77) = (C0134I+C0234I+2*(C134I(1,1)-C134I(2,1)+C234I(2,1)-2*(
     -   C234I(2,2)-C234I(3,2)+D1234I(7,2)+2*D1234I(13,3)-4*(D1234I(1
     -   8,4)-D1234I(21,4)))))*P(8)+2*(D1234I(6,2)*P(108)+D1234I(3,2)
     -   *P(109))+D1234I(2,1)*P(146)+D01234I*P(373)+D1234I(3,1)*P(411
     -   )
       F(77)=DCMPLX(FR(77),FI(77))
       P(416) = p2sq-p4sq+s12+s23
       P(417) = za3*za4*zb2-za4*zb2*zb3+za2*za3*zb4+za2*zb3*zb4-2*zb2*
     -   zb3*zb4
       P(418) = zb3*zb4-za3*P(13)
       P(419) = za4*zb2*P(33)+za2*P(418)
       P(420) = zb1*P(417)+za1*P(419)
       P(421) = P(8)*P(416)+2*M**2*P(420)
       FR(78) = (C0134R+C0234R)*P(41)-2*((C134R(2,1)+4*(D1234R(12,3)-D
     -   1234R(13,3))-8*(D1234R(18,4)-D1234R(21,4)))*P(8)-D1234R(3,2)
     -   *P(132)-D1234R(6,2)*P(135))+D01234R*P(296)+D1234R(2,1)*P(377
     -   )+D1234R(3,1)*P(421)
       FI(78) = (C0134I+C0234I)*P(41)-2*((C134I(2,1)+4*(D1234I(12,3)-D
     -   1234I(13,3))-8*(D1234I(18,4)-D1234I(21,4)))*P(8)-D1234I(3,2)
     -   *P(132)-D1234I(6,2)*P(135))+D01234I*P(296)+D1234I(2,1)*P(377
     -   )+D1234I(3,1)*P(421)
       F(78)=DCMPLX(FR(78),FI(78))
       P(422) = -(za3*za4*zb2)+za4*zb2*zb3+za2*za3*zb4+za2*zb3*zb4-2*z
     -   b2*zb3*zb4
       P(423) = -2*za3*za4+za3*zb4-zb3*zb4
       P(424) = za4*zb2*P(33)+za2*P(423)
       P(425) = zb1*P(422)+za1*P(424)
       P(426) = P(8)*P(367)+2*M**2*P(425)
       FR(79) = (C0134R+C0234R+2*(C134R(1,1)-C134R(2,1)+8*(D1234R(18,4
     -   )-D1234R(21,4))))*P(8)+2*(D1234R(3,2)*P(132)+D1234R(6,2)*P(1
     -   35))+D1234R(2,1)*P(146)+D01234R*P(373)+D1234R(3,1)*P(426)
       FI(79) = (C0134I+C0234I+2*(C134I(1,1)-C134I(2,1)+8*(D1234I(18,4
     -   )-D1234I(21,4))))*P(8)+2*(D1234I(3,2)*P(132)+D1234I(6,2)*P(1
     -   35))+D1234I(2,1)*P(146)+D01234I*P(373)+D1234I(3,1)*P(426)
       F(79)=DCMPLX(FR(79),FI(79))
       P(427) = -2*M**2*P(1)*P(33)*P(61)+P(8)*P(144)
       FR(80) = (C0124R+2*(C124R(1,1)-C134R(2,1)-2*(C124R(2,2)-C124R(3
     -   ,2)+D1234R(7,2)+2*D1234R(11,3)-4*(D1234R(18,4)-D1234R(20,4))
     -   )))*P(8)+C0134R*P(41)+2*(D1234R(3,2)*P(132)+D1234R(5,2)*P(13
     -   5))+D01234R*P(282)+D1234R(1,1)*P(377)+D1234R(3,1)*P(427)
       FI(80) = (C0124I+2*(C124I(1,1)-C134I(2,1)-2*(C124I(2,2)-C124I(3
     -   ,2)+D1234I(7,2)+2*D1234I(11,3)-4*(D1234I(18,4)-D1234I(20,4))
     -   )))*P(8)+C0134I*P(41)+2*(D1234I(3,2)*P(132)+D1234I(5,2)*P(13
     -   5))+D01234I*P(282)+D1234I(1,1)*P(377)+D1234I(3,1)*P(427)
       F(80)=DCMPLX(FR(80),FI(80))
       FR(81) = (C0134R+2*(C134R(1,1)-C134R(2,1)-2*(C134R(2,2)-C134R(3
     -   ,2)+D1234R(7,2)+2*(D1234R(11,3)-D1234R(13,3))-4*(D1234R(18,4
     -   )-D1234R(20,4)))))*P(8)+2*D1234R(3,2)*P(37)+C0124R*P(41)+(D1
     -   234R(1,1)+2*D1234R(5,2))*P(146)+D01234R*P(283)+D1234R(3,1)*P
     -   (427)
       FI(81) = (C0134I+2*(C134I(1,1)-C134I(2,1)-2*(C134I(2,2)-C134I(3
     -   ,2)+D1234I(7,2)+2*(D1234I(11,3)-D1234I(13,3))-4*(D1234I(18,4
     -   )-D1234I(20,4)))))*P(8)+2*D1234I(3,2)*P(37)+C0124I*P(41)+(D1
     -   234I(1,1)+2*D1234I(5,2))*P(146)+D01234I*P(283)+D1234I(3,1)*P
     -   (427)
       F(81)=DCMPLX(FR(81),FI(81))
       FR(82) = (C0134R+2*(C134R(1,1)-C134R(2,1)+2*(-C124R(2,2)+C124R(
     -   3,2)+D1234R(7,2)+4*(D1234R(18,4)-D1234R(20,4)))))*P(8)+C0124
     -   R*P(41)+2*(D1234R(3,2)*P(132)+D1234R(5,2)*P(135))+D1234R(1,1
     -   )*P(146)+D01234R*P(283)+D1234R(3,1)*P(427)
       FI(82) = (C0134I+2*(C134I(1,1)-C134I(2,1)+2*(-C124I(2,2)+C124I(
     -   3,2)+D1234I(7,2)+4*(D1234I(18,4)-D1234I(20,4)))))*P(8)+C0124
     -   I*P(41)+2*(D1234I(3,2)*P(132)+D1234I(5,2)*P(135))+D1234I(1,1
     -   )*P(146)+D01234I*P(283)+D1234I(3,1)*P(427)
       F(82)=DCMPLX(FR(82),FI(82))
       P(428) = p2sq+p3sq+2*p4sq-s23
       P(429) = -(za2*za3)+za3*zb2+za2*zb3-3*zb2*zb3
       P(430) = 2*za4*zb2*zb3+zb4*P(429)
       P(431) = za4*zb2*P(115)+za2*P(407)
       P(432) = -(zb1*P(430))+za1*P(431)
       P(433) = -(P(8)*P(428))+2*M**2*P(432)
       FR(83) = C0134R*P(8)+C0124R*P(41)+2*((-C124R(1,1)+C124R(2,1)+C2
     -   34R(2,1)-4*D1234R(13,3)+8*(D1234R(18,4)-D1234R(20,4)))*P(8)+
     -   D1234R(3,2)*P(210)+D1234R(5,2)*P(213))+D01234R*P(283)+D1234R
     -   (1,1)*P(412)+D1234R(3,1)*P(433)
       FI(83) = C0134I*P(8)+C0124I*P(41)+2*((-C124I(1,1)+C124I(2,1)+C2
     -   34I(2,1)-4*D1234I(13,3)+8*(D1234I(18,4)-D1234I(20,4)))*P(8)+
     -   D1234I(3,2)*P(210)+D1234I(5,2)*P(213))+D01234I*P(283)+D1234I
     -   (1,1)*P(412)+D1234I(3,1)*P(433)
       F(83)=DCMPLX(FR(83),FI(83))
       P(434) = p3sq+2*p4sq-s12
       P(435) = 2*M**2*P(19)*P(189)-P(8)*P(434)
       FR(84) = 2*((-C0123R-3*C123R(1,1)+C234R(2,1)-2*(C123R(1,2)-2*(D
     -   1234R(7,2)+3*(D1234R(11,3)-D1234R(13,3))+2*(D1234R(16,4)+D12
     -   34R(18,4)-2*D1234R(20,4)))))*P(8)+(D1234R(1,2)+D1234R(3,2))*
     -   P(210)+2*D1234R(5,2)*P(213)+D01234R*P(273)+(-D1234R(1,1)+D12
     -   34R(3,1))*P(435))
       FI(84) = 2*((-C0123I-3*C123I(1,1)+C234I(2,1)-2*(C123I(1,2)-2*(D
     -   1234I(7,2)+3*(D1234I(11,3)-D1234I(13,3))+2*(D1234I(16,4)+D12
     -   34I(18,4)-2*D1234I(20,4)))))*P(8)+(D1234I(1,2)+D1234I(3,2))*
     -   P(210)+2*D1234I(5,2)*P(213)+D01234I*P(273)+(-D1234I(1,1)+D12
     -   34I(3,1))*P(435))
       F(84)=DCMPLX(FR(84),FI(84))
       P(436) = p3sq+3*p4sq-s12
       P(437) = 2*M**2*P(19)*P(126)-P(8)*P(436)
       FR(85) = (C0123R+C0234R)*P(41)-2*((C123R(1,1)+C123R(2,1)+C234R(
     -   1,1)-C234R(2,1)-2*(-C123R(3,2)+D1234R(7,2)+2*(D1234R(11,3)+D
     -   1234R(12,3)-2*(D1234R(13,3)-D1234R(18,4)-D1234R(19,4)+D1234R
     -   (20,4)+D1234R(21,4)))))*P(8)-(D1234R(3,2)+D1234R(4,2))*P(210
     -   )-(D1234R(5,2)+D1234R(6,2))*P(213))+D01234R*P(273)+D1234R(1,
     -   1)*P(297)+D1234R(2,1)*P(357)+D1234R(3,1)*P(437)
       FI(85) = (C0123I+C0234I)*P(41)-2*((C123I(1,1)+C123I(2,1)+C234I(
     -   1,1)-C234I(2,1)-2*(-C123I(3,2)+D1234I(7,2)+2*(D1234I(11,3)+D
     -   1234I(12,3)-2*(D1234I(13,3)-D1234I(18,4)-D1234I(19,4)+D1234I
     -   (20,4)+D1234I(21,4)))))*P(8)-(D1234I(3,2)+D1234I(4,2))*P(210
     -   )-(D1234I(5,2)+D1234I(6,2))*P(213))+D01234I*P(273)+D1234I(1,
     -   1)*P(297)+D1234I(2,1)*P(357)+D1234I(3,1)*P(437)
       F(85)=DCMPLX(FR(85),FI(85))
       P(438) = p1sq-p3sq-2*p4sq+s12-s23
       P(439) = -2*M**2*P(19)*P(123)-P(8)*P(438)
       FR(86) = (C0234R-2*(-C234R(2,1)+2*(C123R(2,1)+C123R(3,2)-D1234R
     -   (7,2)-4*(D1234R(12,3)-D1234R(13,3)+D1234R(18,4)+D1234R(19,4)
     -   -D1234R(20,4)-D1234R(21,4)))))*P(8)+C0123R*P(41)+D1234R(1,1)
     -   *P(109)+2*((D1234R(3,2)+D1234R(4,2))*P(210)+(D1234R(5,2)+D12
     -   34R(6,2))*P(213))+D01234R*P(273)+D1234R(3,1)*P(437)+D1234R(2
     -   ,1)*P(439)
       FI(86) = (C0234I-2*(-C234I(2,1)+2*(C123I(2,1)+C123I(3,2)-D1234I
     -   (7,2)-4*(D1234I(12,3)-D1234I(13,3)+D1234I(18,4)+D1234I(19,4)
     -   -D1234I(20,4)-D1234I(21,4)))))*P(8)+C0123I*P(41)+D1234I(1,1)
     -   *P(109)+2*((D1234I(3,2)+D1234I(4,2))*P(210)+(D1234I(5,2)+D12
     -   34I(6,2))*P(213))+D01234I*P(273)+D1234I(3,1)*P(437)+D1234I(2
     -   ,1)*P(439)
       F(86)=DCMPLX(FR(86),FI(86))
       FR(87) = 2*(-((C123R(2,1)+2*C123R(2,2)+C234R(1,1)-C234R(2,1)-4*
     -   (D1234R(12,3)-D1234R(13,3))-8*(D1234R(17,4)+D1234R(18,4))+16
     -   *D1234R(21,4))*P(8))+(D1234R(2,2)+D1234R(3,2))*P(210)+2*D123
     -   4R(6,2)*P(213)+(D1234R(2,1)-D1234R(3,1))*P(273))
       FI(87) = 2*(-((C123I(2,1)+2*C123I(2,2)+C234I(1,1)-C234I(2,1)-4*
     -   (D1234I(12,3)-D1234I(13,3))-8*(D1234I(17,4)+D1234I(18,4))+16
     -   *D1234I(21,4))*P(8))+(D1234I(2,2)+D1234I(3,2))*P(210)+2*D123
     -   4I(6,2)*P(213)+(D1234I(2,1)-D1234I(3,1))*P(273))
       F(87)=DCMPLX(FR(87),FI(87))
       P(440) = p1sq+2*p4sq-s23
       P(441) = 2*M**2*P(19)*P(190)-P(8)*P(440)
       FR(88) = 2*(-((C0123R+C123R(1,1)+C234R(2,1)+2*(C234R(2,2)-2*(D1
     -   234R(7,2)+3*(D1234R(11,3)-D1234R(13,3))+2*(D1234R(16,4)+D123
     -   4R(18,4)-2*D1234R(20,4)))))*P(8))+2*D1234R(5,2)*P(108)+(D123
     -   4R(1,2)+D1234R(3,2))*P(109)+D01234R*P(273)+(-D1234R(1,1)+D12
     -   34R(3,1))*P(441))
       FI(88) = 2*(-((C0123I+C123I(1,1)+C234I(2,1)+2*(C234I(2,2)-2*(D1
     -   234I(7,2)+3*(D1234I(11,3)-D1234I(13,3))+2*(D1234I(16,4)+D123
     -   4I(18,4)-2*D1234I(20,4)))))*P(8))+2*D1234I(5,2)*P(108)+(D123
     -   4I(1,2)+D1234I(3,2))*P(109)+D01234I*P(273)+(-D1234I(1,1)+D12
     -   34I(3,1))*P(441))
       F(88)=DCMPLX(FR(88),FI(88))
       P(442) = p2sq-p3sq+2*p4sq-s23
       P(443) = za2*P(2)+zb2*P(33)
       P(444) = -2*za4*zb2*zb3+zb4*P(443)
       P(445) = za3*za4+za4*zb3-2*za3*zb4
       P(446) = za4*zb2*P(115)+za2*P(445)
       P(447) = zb1*P(444)+za1*P(446)
       P(448) = -(P(8)*P(442))+2*M**2*P(447)
       FR(89) = (C0124R+C0134R-2*(C234R(2,1)+2*(C234R(2,2)+2*(D1234R(7
     -   ,2)+D1234R(11,3)+D1234R(13,3)-2*(D1234R(18,4)-D1234R(20,4)))
     -   )))*P(8)+2*(D1234R(5,2)*P(108)+D1234R(3,2)*P(109))+D1234R(1,
     -   1)*P(374)+D01234R*P(388)+D1234R(3,1)*P(448)
       FI(89) = (C0124I+C0134I-2*(C234I(2,1)+2*(C234I(2,2)+2*(D1234I(7
     -   ,2)+D1234I(11,3)+D1234I(13,3)-2*(D1234I(18,4)-D1234I(20,4)))
     -   )))*P(8)+2*(D1234I(5,2)*P(108)+D1234I(3,2)*P(109))+D1234I(1,
     -   1)*P(374)+D01234I*P(388)+D1234I(3,1)*P(448)
       F(89)=DCMPLX(FR(89),FI(89))
       P(449) = P(8)*P(396)+2*M**2*P(401)
       FR(90) = C0134R*P(8)+C0124R*P(41)+2*(C134R(1,1)*P(8)-(C134R(2,1
     -   )+C234R(2,1)+2*C234R(2,2)+4*D1234R(13,3)-8*(D1234R(18,4)-D12
     -   34R(20,4)))*P(8)+D1234R(5,2)*P(108)+D1234R(3,2)*P(109))+D123
     -   4R(1,1)*P(146)+D01234R*P(283)+D1234R(3,1)*P(449)
       FI(90) = C0134I*P(8)+C0124I*P(41)+2*(C134I(1,1)*P(8)-(C134I(2,1
     -   )+C234I(2,1)+2*C234I(2,2)+4*D1234I(13,3)-8*(D1234I(18,4)-D12
     -   34I(20,4)))*P(8)+D1234I(5,2)*P(108)+D1234I(3,2)*P(109))+D123
     -   4I(1,1)*P(146)+D01234I*P(283)+D1234I(3,1)*P(449)
       F(90)=DCMPLX(FR(90),FI(90))
       FR(91) = -4*(C134R(2,1)+C134R(2,2)-4*(D1234R(13,3)+D1234R(18,4)
     -   ))*P(8)+2*(D1234R(3,1)+D1234R(3,2))*P(37)
       FI(91) = -4*(C134I(2,1)+C134I(2,2)-4*(D1234I(13,3)+D1234I(18,4)
     -   ))*P(8)+2*(D1234I(3,1)+D1234I(3,2))*P(37)
       F(91)=DCMPLX(FR(91),FI(91))
       P(450) = s12*P(8)-2*M**2*P(24)*P(207)
       FR(92) = 2*((C124R(2,1)+C234R(2,1)+8*D1234R(18,4))*P(8)+D1234R(
     -   3,2)*P(210)-D1234R(3,1)*P(450))
       FI(92) = 2*((C124I(2,1)+C234I(2,1)+8*D1234I(18,4))*P(8)+D1234I(
     -   3,2)*P(210)-D1234I(3,1)*P(450))
       F(92)=DCMPLX(FR(92),FI(92))
       P(451) = p2sq*P(8)-2*M**2*P(1)*P(236)
       FR(93) = 2*(-((C124R(2,1)+2*C124R(2,2)+C134R(2,1)-4*D1234R(13,3
     -   )-8*D1234R(18,4))*P(8))+D1234R(3,2)*P(132)+D1234R(3,1)*P(451
     -   ))
       FI(93) = 2*(-((C124I(2,1)+2*C124I(2,2)+C134I(2,1)-4*D1234I(13,3
     -   )-8*D1234I(18,4))*P(8))+D1234I(3,2)*P(132)+D1234I(3,1)*P(451
     -   ))
       F(93)=DCMPLX(FR(93),FI(93))
       FR(94) = 2*(-((C134R(2,1)+C234R(2,1)+2*C234R(2,2)-4*D1234R(13,3
     -   )-8*D1234R(18,4))*P(8))+D1234R(3,2)*P(109)+D1234R(3,1)*P(296
     -   ))
       FI(94) = 2*(-((C134I(2,1)+C234I(2,1)+2*C234I(2,2)-4*D1234I(13,3
     -   )-8*D1234I(18,4))*P(8))+D1234I(3,2)*P(109)+D1234I(3,1)*P(296
     -   ))
       F(94)=DCMPLX(FR(94),FI(94))
       FR(95) = 2*((-C124R(2,1)+C134R(2,1)+4*D1234R(13,3)+8*D1234R(18,
     -   4))*P(8)-D1234R(3,2)*P(270)+D1234R(3,1)*P(283))
       FI(95) = 2*((-C124I(2,1)+C134I(2,1)+4*D1234I(13,3)+8*D1234I(18,
     -   4))*P(8)-D1234I(3,2)*P(270)+D1234I(3,1)*P(283))
       F(95)=DCMPLX(FR(95),FI(95))
       P(452) = -(p2sq*P(8))+s12*P(8)-2*M**2*P(9)*P(223)
       FR(96) = 2*((C134R(2,1)-C234R(2,1)+4*D1234R(13,3)+8*D1234R(18,4
     -   ))*P(8)-D1234R(3,2)*P(295)+D1234R(3,1)*P(452))
       FI(96) = 2*((C134I(2,1)-C234I(2,1)+4*D1234I(13,3)+8*D1234I(18,4
     -   ))*P(8)-D1234I(3,2)*P(295)+D1234I(3,1)*P(452))
       F(96)=DCMPLX(FR(96),FI(96))
       FR(97) = -4*(D1234R(3,2)-D1234R(5,2)+2*(D1234R(3,3)-D1234R(7,3)
     -   -D1234R(9,3)+D1234R(10,3)-4*(D1234R(11,4)-D1234R(12,4))+2*(-
     -   D1234R(3,4)+D1234R(7,4)-D1234R(8,4)+D1234R(15,4))))*P(8)
       FI(97) = -4*(D1234I(3,2)-D1234I(5,2)+2*(D1234I(3,3)-D1234I(7,3)
     -   -D1234I(9,3)+D1234I(10,3)-4*(D1234I(11,4)-D1234I(12,4))+2*(-
     -   D1234I(3,4)+D1234I(7,4)-D1234I(8,4)+D1234I(15,4))))*P(8)
       F(97)=DCMPLX(FR(97),FI(97))
       FR(98) = 4*(D1234R(2,1)-D1234R(3,1)+D1234R(4,2)-D1234R(5,2)-2*(
     -   D1234R(3,3)+D1234R(8,3)-2*(D1234R(3,4)+D1234R(9,3)-D1234R(10
     -   ,4)-D1234R(12,4)+D1234R(14,4)+2*(D1234R(11,4)-D1234R(15,4)))
     -   ))*P(8)
       FI(98) = 4*(D1234I(2,1)-D1234I(3,1)+D1234I(4,2)-D1234I(5,2)-2*(
     -   D1234I(3,3)+D1234I(8,3)-2*(D1234I(3,4)+D1234I(9,3)-D1234I(10
     -   ,4)-D1234I(12,4)+D1234I(14,4)+2*(D1234I(11,4)-D1234I(15,4)))
     -   ))*P(8)
       F(98)=DCMPLX(FR(98),FI(98))
       FR(99) = -4*(D1234R(3,1)-D1234R(6,2)+2*(D1234R(3,2)+D1234R(7,3)
     -   -D1234R(10,3))-4*(D1234R(3,4)+D1234R(11,4)-D1234R(12,4)-D123
     -   4R(15,4)))*P(8)
       FI(99) = -4*(D1234I(3,1)-D1234I(6,2)+2*(D1234I(3,2)+D1234I(7,3)
     -   -D1234I(10,3))-4*(D1234I(3,4)+D1234I(11,4)-D1234I(12,4)-D123
     -   4I(15,4)))*P(8)
       F(99)=DCMPLX(FR(99),FI(99))
       FR(100) = -4*(D1234R(1,1)+D1234R(1,2)-D1234R(3,1)-3*D1234R(5,2)
     -   +2*(D1234R(3,2)+D1234R(3,3)+D1234R(4,2)+D1234R(4,3)-D1234R(6
     -   ,2)-D1234R(7,3)-D1234R(10,3)-2*(D1234R(3,4)-D1234R(7,4)+D123
     -   4R(8,4)+2*(D1234R(11,4)-D1234R(12,4))-D1234R(15,4))))*P(8)
       FI(100) = -4*(D1234I(1,1)+D1234I(1,2)-D1234I(3,1)-3*D1234I(5,2)
     -   +2*(D1234I(3,2)+D1234I(3,3)+D1234I(4,2)+D1234I(4,3)-D1234I(6
     -   ,2)-D1234I(7,3)-D1234I(10,3)-2*(D1234I(3,4)-D1234I(7,4)+D123
     -   4I(8,4)+2*(D1234I(11,4)-D1234I(12,4))-D1234I(15,4))))*P(8)
       F(100)=DCMPLX(FR(100),FI(100))
       FR(101) = -4*(D1234R(2,1)+D1234R(2,2)-D1234R(3,1)-3*D1234R(6,2)
     -   +2*(D1234R(3,2)+D1234R(3,3)+D1234R(4,2)-D1234R(5,2)+D1234R(6
     -   ,3)-D1234R(9,3)-D1234R(10,3))+4*(-D1234R(3,4)+D1234R(10,4)+D
     -   1234R(12,4)-D1234R(14,4))-8*(D1234R(11,4)-D1234R(15,4)))*P(8
     -   )
       FI(101) = -4*(D1234I(2,1)+D1234I(2,2)-D1234I(3,1)-3*D1234I(6,2)
     -   +2*(D1234I(3,2)+D1234I(3,3)+D1234I(4,2)-D1234I(5,2)+D1234I(6
     -   ,3)-D1234I(9,3)-D1234I(10,3))+4*(-D1234I(3,4)+D1234I(10,4)+D
     -   1234I(12,4)-D1234I(14,4))-8*(D1234I(11,4)-D1234I(15,4)))*P(8
     -   )
       F(101)=DCMPLX(FR(101),FI(101))
       FR(102) = 4*(D1234R(3,1)-3*D1234R(3,2)+2*(D1234R(5,2)+D1234R(6,
     -   2)-D1234R(7,3)-D1234R(9,3)+2*(D1234R(3,4)+D1234R(10,3)+D1234
     -   R(11,4)-D1234R(12,4)-D1234R(15,4))))*P(8)
       FI(102) = 4*(D1234I(3,1)-3*D1234I(3,2)+2*(D1234I(5,2)+D1234I(6,
     -   2)-D1234I(7,3)-D1234I(9,3)+2*(D1234I(3,4)+D1234I(10,3)+D1234
     -   I(11,4)-D1234I(12,4)-D1234I(15,4))))*P(8)
       F(102)=DCMPLX(FR(102),FI(102))
       FR(103) = -4*(D1234R(2,1)-D1234R(3,1)+3*(D1234R(4,2)-D1234R(6,2
     -   ))+2*(D1234R(3,2)+D1234R(3,3)+D1234R(4,3)-D1234R(5,2)-D1234R
     -   (5,3)-D1234R(9,3)-4*(D1234R(11,4)-D1234R(12,4))+2*(-D1234R(3
     -   ,4)+D1234R(7,4)-D1234R(8,4)+D1234R(15,4))))*P(8)
       FI(103) = -4*(D1234I(2,1)-D1234I(3,1)+3*(D1234I(4,2)-D1234I(6,2
     -   ))+2*(D1234I(3,2)+D1234I(3,3)+D1234I(4,3)-D1234I(5,2)-D1234I
     -   (5,3)-D1234I(9,3)-4*(D1234I(11,4)-D1234I(12,4))+2*(-D1234I(3
     -   ,4)+D1234I(7,4)-D1234I(8,4)+D1234I(15,4))))*P(8)
       F(103)=DCMPLX(FR(103),FI(103))
       FR(104) = -4*(D1234R(2,1)-D1234R(3,1)+D1234R(4,2)-D1234R(5,2)+2
     -   *(D1234R(2,2)+D1234R(3,2)+D1234R(3,3)+D1234R(6,3)+D1234R(7,3
     -   )+D1234R(8,3)-2*(D1234R(3,4)+D1234R(6,2)+D1234R(9,3)+D1234R(
     -   10,3)-D1234R(10,4)-D1234R(12,4)+D1234R(14,4)+2*(D1234R(11,4)
     -   -D1234R(15,4)))))*P(8)
       FI(104) = -4*(D1234I(2,1)-D1234I(3,1)+D1234I(4,2)-D1234I(5,2)+2
     -   *(D1234I(2,2)+D1234I(3,2)+D1234I(3,3)+D1234I(6,3)+D1234I(7,3
     -   )+D1234I(8,3)-2*(D1234I(3,4)+D1234I(6,2)+D1234I(9,3)+D1234I(
     -   10,3)-D1234I(10,4)-D1234I(12,4)+D1234I(14,4)+2*(D1234I(11,4)
     -   -D1234I(15,4)))))*P(8)
       F(104)=DCMPLX(FR(104),FI(104))
       FR(105) = 4*(D1234R(3,1)-3*D1234R(3,2)+4*(D1234R(3,4)+D1234R(6,
     -   2)-D1234R(7,3)+D1234R(10,3)+D1234R(11,4)-D1234R(12,4)-D1234R
     -   (15,4)))*P(8)
       FI(105) = 4*(D1234I(3,1)-3*D1234I(3,2)+4*(D1234I(3,4)+D1234I(6,
     -   2)-D1234I(7,3)+D1234I(10,3)+D1234I(11,4)-D1234I(12,4)-D1234I
     -   (15,4)))*P(8)
       F(105)=DCMPLX(FR(105),FI(105))
       FR(106) = 4*(D1234R(3,2)-D1234R(5,2)-2*(D1234R(5,3)-3*D1234R(7,
     -   3)-D1234R(9,3)+D1234R(10,3)-4*(D1234R(11,4)-D1234R(12,4))+2*
     -   (D1234R(3,3)-D1234R(3,4)+D1234R(7,4)-D1234R(8,4)+D1234R(15,4
     -   ))))*P(8)
       FI(106) = 4*(D1234I(3,2)-D1234I(5,2)-2*(D1234I(5,3)-3*D1234I(7,
     -   3)-D1234I(9,3)+D1234I(10,3)-4*(D1234I(11,4)-D1234I(12,4))+2*
     -   (D1234I(3,3)-D1234I(3,4)+D1234I(7,4)-D1234I(8,4)+D1234I(15,4
     -   ))))*P(8)
       F(106)=DCMPLX(FR(106),FI(106))
       FR(107) = -4*(D1234R(2,1)+D1234R(2,2)-D1234R(3,1)-D1234R(6,2)+2
     -   *(-D1234R(7,3)+D1234R(8,3)-3*D1234R(9,3)+D1234R(10,3)+2*(D12
     -   34R(3,3)-D1234R(3,4)+D1234R(10,4)+D1234R(12,4)-D1234R(14,4)-
     -   2*(D1234R(11,4)-D1234R(15,4)))))*P(8)
       FI(107) = -4*(D1234I(2,1)+D1234I(2,2)-D1234I(3,1)-D1234I(6,2)+2
     -   *(-D1234I(7,3)+D1234I(8,3)-3*D1234I(9,3)+D1234I(10,3)+2*(D12
     -   34I(3,3)-D1234I(3,4)+D1234I(10,4)+D1234I(12,4)-D1234I(14,4)-
     -   2*(D1234I(11,4)-D1234I(15,4)))))*P(8)
       F(107)=DCMPLX(FR(107),FI(107))
       FR(108) = 4*(D1234R(3,1)+D1234R(5,2)-2*(D1234R(3,2)+D1234R(3,3)
     -   -D1234R(6,2)-D1234R(10,3)-2*(D1234R(3,4)+D1234R(11,4)-D1234R
     -   (12,4)-D1234R(15,4))))*P(8)
       FI(108) = 4*(D1234I(3,1)+D1234I(5,2)-2*(D1234I(3,2)+D1234I(3,3)
     -   -D1234I(6,2)-D1234I(10,3)-2*(D1234I(3,4)+D1234I(11,4)-D1234I
     -   (12,4)-D1234I(15,4))))*P(8)
       F(108)=DCMPLX(FR(108),FI(108))
       FR(109) = -4*(D1234R(1,1)+D1234R(1,2)-D1234R(3,1)-D1234R(5,2)+2
     -   *(D1234R(5,3)-3*D1234R(7,3)-D1234R(9,3)+D1234R(10,3)-4*(D123
     -   4R(11,4)-D1234R(12,4))+2*(D1234R(3,3)-D1234R(3,4)+D1234R(7,4
     -   )-D1234R(8,4)+D1234R(15,4))))*P(8)
       FI(109) = -4*(D1234I(1,1)+D1234I(1,2)-D1234I(3,1)-D1234I(5,2)+2
     -   *(D1234I(5,3)-3*D1234I(7,3)-D1234I(9,3)+D1234I(10,3)-4*(D123
     -   4I(11,4)-D1234I(12,4))+2*(D1234I(3,3)-D1234I(3,4)+D1234I(7,4
     -   )-D1234I(8,4)+D1234I(15,4))))*P(8)
       F(109)=DCMPLX(FR(109),FI(109))
       FR(110) = 4*(D1234R(3,2)-D1234R(6,2)-2*(-D1234R(7,3)+D1234R(8,3
     -   )-3*D1234R(9,3)+D1234R(10,3)+2*(D1234R(3,3)-D1234R(3,4)+D123
     -   4R(10,4)+D1234R(12,4)-D1234R(14,4)-2*(D1234R(11,4)-D1234R(15
     -   ,4)))))*P(8)
       FI(110) = 4*(D1234I(3,2)-D1234I(6,2)-2*(-D1234I(7,3)+D1234I(8,3
     -   )-3*D1234I(9,3)+D1234I(10,3)+2*(D1234I(3,3)-D1234I(3,4)+D123
     -   4I(10,4)+D1234I(12,4)-D1234I(14,4)-2*(D1234I(11,4)-D1234I(15
     -   ,4)))))*P(8)
       F(110)=DCMPLX(FR(110),FI(110))
       FR(111) = 4*(D1234R(3,1)+D1234R(6,2)-2*(D1234R(3,2)+D1234R(3,3)
     -   -D1234R(5,2)-D1234R(10,3))+4*(D1234R(3,4)+D1234R(11,4)-D1234
     -   R(12,4)-D1234R(15,4)))*P(8)
       FI(111) = 4*(D1234I(3,1)+D1234I(6,2)-2*(D1234I(3,2)+D1234I(3,3)
     -   -D1234I(5,2)-D1234I(10,3))+4*(D1234I(3,4)+D1234I(11,4)-D1234
     -   I(12,4)-D1234I(15,4)))*P(8)
       F(111)=DCMPLX(FR(111),FI(111))
       FR(112) = 4*(D1234R(2,1)-D1234R(3,1)+D1234R(4,2)-D1234R(6,2)+2*
     -   (D1234R(3,2)-3*D1234R(3,3)-D1234R(5,2)-D1234R(5,3)-2*(-D1234
     -   R(3,4)+D1234R(7,4)-D1234R(8,4)-D1234R(9,3)+D1234R(10,3)-2*(D
     -   1234R(7,3)+D1234R(11,4)-D1234R(12,4))+D1234R(15,4))))*P(8)
       FI(112) = 4*(D1234I(2,1)-D1234I(3,1)+D1234I(4,2)-D1234I(6,2)+2*
     -   (D1234I(3,2)-3*D1234I(3,3)-D1234I(5,2)-D1234I(5,3)-2*(-D1234
     -   I(3,4)+D1234I(7,4)-D1234I(8,4)-D1234I(9,3)+D1234I(10,3)-2*(D
     -   1234I(7,3)+D1234I(11,4)-D1234I(12,4))+D1234I(15,4))))*P(8)
       F(112)=DCMPLX(FR(112),FI(112))
       FR(113) = 4*(D1234R(3,2)-D1234R(6,2)-2*(3*D1234R(3,3)-D1234R(7,
     -   3)-5*D1234R(9,3)+D1234R(10,3)+2*(-D1234R(3,4)+D1234R(8,3)+D1
     -   234R(10,4)+D1234R(12,4)-D1234R(14,4)-2*(D1234R(11,4)-D1234R(
     -   15,4)))))*P(8)
       FI(113) = 4*(D1234I(3,2)-D1234I(6,2)-2*(3*D1234I(3,3)-D1234I(7,
     -   3)-5*D1234I(9,3)+D1234I(10,3)+2*(-D1234I(3,4)+D1234I(8,3)+D1
     -   234I(10,4)+D1234I(12,4)-D1234I(14,4)-2*(D1234I(11,4)-D1234I(
     -   15,4)))))*P(8)
       F(113)=DCMPLX(FR(113),FI(113))
       FR(114) = 4*(D1234R(3,1)+D1234R(5,2)-2*(D1234R(3,2)-D1234R(6,2)
     -   -D1234R(9,3)-D1234R(10,3)+2*(D1234R(3,3)-D1234R(3,4)-D1234R(
     -   11,4)+D1234R(12,4)+D1234R(15,4))))*P(8)
       FI(114) = 4*(D1234I(3,1)+D1234I(5,2)-2*(D1234I(3,2)-D1234I(6,2)
     -   -D1234I(9,3)-D1234I(10,3)+2*(D1234I(3,3)-D1234I(3,4)-D1234I(
     -   11,4)+D1234I(12,4)+D1234I(15,4))))*P(8)
       F(114)=DCMPLX(FR(114),FI(114))
       FR(115) = -8*(D1234R(1,1)+D1234R(1,3)-D1234R(3,1)+D1234R(3,2)+2
     -   *(D1234R(1,2)+D1234R(3,3)-D1234R(3,4)+D1234R(5,4))-3*(D1234R
     -   (5,2)+D1234R(7,3))+6*(-D1234R(8,4)+D1234R(12,4)))*P(8)
       FI(115) = -8*(D1234I(1,1)+D1234I(1,3)-D1234I(3,1)+D1234I(3,2)+2
     -   *(D1234I(1,2)+D1234I(3,3)-D1234I(3,4)+D1234I(5,4))-3*(D1234I
     -   (5,2)+D1234I(7,3))+6*(-D1234I(8,4)+D1234I(12,4)))*P(8)
       F(115)=DCMPLX(FR(115),FI(115))
       FR(116) = -4*(D1234R(1,1)+D1234R(1,2)+D1234R(2,1)+3*D1234R(4,2)
     -   -5*D1234R(5,2)-D1234R(6,2)-2*(D1234R(3,1)-D1234R(3,2)-D1234R
     -   (4,3)+D1234R(5,3)+D1234R(7,3)-D1234R(10,3)+4*(D1234R(11,4)-D
     -   1234R(12,4))-2*(D1234R(3,3)-D1234R(3,4)+D1234R(7,4)-D1234R(8
     -   ,4)-D1234R(9,3)+D1234R(15,4))))*P(8)
       FI(116) = -4*(D1234I(1,1)+D1234I(1,2)+D1234I(2,1)+3*D1234I(4,2)
     -   -5*D1234I(5,2)-D1234I(6,2)-2*(D1234I(3,1)-D1234I(3,2)-D1234I
     -   (4,3)+D1234I(5,3)+D1234I(7,3)-D1234I(10,3)+4*(D1234I(11,4)-D
     -   1234I(12,4))-2*(D1234I(3,3)-D1234I(3,4)+D1234I(7,4)-D1234I(8
     -   ,4)-D1234I(9,3)+D1234I(15,4))))*P(8)
       F(116)=DCMPLX(FR(116),FI(116))
       FR(117) = 8*(D1234R(3,1)-D1234R(3,3)+3*D1234R(5,2)-D1234R(7,3)-
     -   2*(D1234R(3,2)-D1234R(3,4)-D1234R(5,3)-D1234R(8,4))-4*D1234R
     -   (12,4))*P(8)
       FI(117) = 8*(D1234I(3,1)-D1234I(3,3)+3*D1234I(5,2)-D1234I(7,3)-
     -   2*(D1234I(3,2)-D1234I(3,4)-D1234I(5,3)-D1234I(8,4))-4*D1234I
     -   (12,4))*P(8)
       F(117)=DCMPLX(FR(117),FI(117))
       FR(118) = 8*(D1234R(3,2)-D1234R(5,2)-3*(D1234R(3,3)+D1234R(5,3)
     -   )+2*(D1234R(3,4)-D1234R(5,4))+6*(D1234R(7,3)+D1234R(8,4)-D12
     -   34R(12,4)))*P(8)
       FI(118) = 8*(D1234I(3,2)-D1234I(5,2)-3*(D1234I(3,3)+D1234I(5,3)
     -   )+2*(D1234I(3,4)-D1234I(5,4))+6*(D1234I(7,3)+D1234I(8,4)-D12
     -   34I(12,4)))*P(8)
       F(118)=DCMPLX(FR(118),FI(118))
       FR(119) = -4*(D1234R(2,1)-D1234R(3,1)-D1234R(3,2)+D1234R(4,2)-D
     -   1234R(5,2)+D1234R(6,2)+6*(D1234R(3,3)-D1234R(7,3)-D1234R(9,3
     -   )+D1234R(10,3))-8*(D1234R(11,4)-D1234R(12,4))+4*(-D1234R(3,4
     -   )+D1234R(7,4)-D1234R(8,4)+D1234R(15,4)))*P(8)
       FI(119) = -4*(D1234I(2,1)-D1234I(3,1)-D1234I(3,2)+D1234I(4,2)-D
     -   1234I(5,2)+D1234I(6,2)+6*(D1234I(3,3)-D1234I(7,3)-D1234I(9,3
     -   )+D1234I(10,3))-8*(D1234I(11,4)-D1234I(12,4))+4*(-D1234I(3,4
     -   )+D1234I(7,4)-D1234I(8,4)+D1234I(15,4)))*P(8)
       F(119)=DCMPLX(FR(119),FI(119))
       FR(120) = 8*(D1234R(3,1)-D1234R(3,2)+D1234R(5,3)+D1234R(7,3)+2*
     -   (-D1234R(3,3)+D1234R(3,4)+D1234R(5,2)+D1234R(8,4)-2*D1234R(1
     -   2,4)))*P(8)
       FI(120) = 8*(D1234I(3,1)-D1234I(3,2)+D1234I(5,3)+D1234I(7,3)+2*
     -   (-D1234I(3,3)+D1234I(3,4)+D1234I(5,2)+D1234I(8,4)-2*D1234I(1
     -   2,4)))*P(8)
       F(120)=DCMPLX(FR(120),FI(120))
       FR(121) = -4*(3*D1234R(3,2)-D1234R(5,2)-2*(-D1234R(3,3)+D1234R(
     -   6,2)-D1234R(7,3)+D1234R(9,3)+D1234R(10,3)+2*(D1234R(3,4)+D12
     -   34R(11,4)-D1234R(12,4)-D1234R(15,4))))*P(8)
       FI(121) = -4*(3*D1234I(3,2)-D1234I(5,2)-2*(-D1234I(3,3)+D1234I(
     -   6,2)-D1234I(7,3)+D1234I(9,3)+D1234I(10,3)+2*(D1234I(3,4)+D12
     -   34I(11,4)-D1234I(12,4)-D1234I(15,4))))*P(8)
       F(121)=DCMPLX(FR(121),FI(121))
       FR(122) = -16*(D1234R(3,2)-D1234R(3,4)+D1234R(7,3)+D1234R(12,4)
     -   )*P(8)
       FI(122) = -16*(D1234I(3,2)-D1234I(3,4)+D1234I(7,3)+D1234I(12,4)
     -   )*P(8)
       F(122)=DCMPLX(FR(122),FI(122))
       FR(123) = -8*(D1234R(3,2)-D1234R(5,2)-D1234R(5,3)+D1234R(7,3)-2
     -   *(D1234R(3,4)+D1234R(8,4))+4*D1234R(12,4))*P(8)
       FI(123) = -8*(D1234I(3,2)-D1234I(5,2)-D1234I(5,3)+D1234I(7,3)-2
     -   *(D1234I(3,4)+D1234I(8,4))+4*D1234I(12,4))*P(8)
       F(123)=DCMPLX(FR(123),FI(123))
       FR(124) = -4*(D1234R(3,2)-D1234R(6,2)-2*(-D1234R(7,3)+D1234R(10
     -   ,3))-4*(D1234R(3,4)+D1234R(11,4)-D1234R(12,4)-D1234R(15,4)))
     -   *P(8)
       FI(124) = -4*(D1234I(3,2)-D1234I(6,2)-2*(-D1234I(7,3)+D1234I(10
     -   ,3))-4*(D1234I(3,4)+D1234I(11,4)-D1234I(12,4)-D1234I(15,4)))
     -   *P(8)
       F(124)=DCMPLX(FR(124),FI(124))
       FR(125) = -8*(D1234R(3,2)-D1234R(3,3)+2*(-D1234R(3,4)+D1234R(7,
     -   3)+D1234R(12,4)))*P(8)
       FI(125) = -8*(D1234I(3,2)-D1234I(3,3)+2*(-D1234I(3,4)+D1234I(7,
     -   3)+D1234I(12,4)))*P(8)
       F(125)=DCMPLX(FR(125),FI(125))
       FR(126) = 8*(-D1234R(3,3)+D1234R(7,3)+2*(D1234R(3,4)+D1234R(8,4
     -   )-2*D1234R(12,4)))*P(8)
       FI(126) = 8*(-D1234I(3,3)+D1234I(7,3)+2*(D1234I(3,4)+D1234I(8,4
     -   )-2*D1234I(12,4)))*P(8)
       F(126)=DCMPLX(FR(126),FI(126))
       FR(127) = -4*(D1234R(5,2)-D1234R(6,2)+2*(D1234R(3,3)-D1234R(9,3
     -   ))-4*(D1234R(3,4)+D1234R(11,4)-D1234R(12,4)-D1234R(15,4)))*P
     -   (8)
       FI(127) = -4*(D1234I(5,2)-D1234I(6,2)+2*(D1234I(3,3)-D1234I(9,3
     -   ))-4*(D1234I(3,4)+D1234I(11,4)-D1234I(12,4)-D1234I(15,4)))*P
     -   (8)
       F(127)=DCMPLX(FR(127),FI(127))
       FR(128) = -8*(D1234R(3,2)+D1234R(7,3)+2*(-D1234R(3,4)+D1234R(12
     -   ,4)))*P(8)
       FI(128) = -8*(D1234I(3,2)+D1234I(7,3)+2*(-D1234I(3,4)+D1234I(12
     -   ,4)))*P(8)
       F(128)=DCMPLX(FR(128),FI(128))
       FR(129) = -4*(D1234R(3,2)-D1234R(5,2)-2*(-D1234R(9,3)+D1234R(10
     -   ,3))-4*(D1234R(3,4)+D1234R(11,4)-D1234R(12,4)-D1234R(15,4)))
     -   *P(8)
       FI(129) = -4*(D1234I(3,2)-D1234I(5,2)-2*(-D1234I(9,3)+D1234I(10
     -   ,3))-4*(D1234I(3,4)+D1234I(11,4)-D1234I(12,4)-D1234I(15,4)))
     -   *P(8)
       F(129)=DCMPLX(FR(129),FI(129))
       FR(130) = -8*(D1234R(3,2)-D1234R(3,3)+2*(-D1234R(3,4)+D1234R(9,
     -   3)+D1234R(15,4)))*P(8)
       FI(130) = -8*(D1234I(3,2)-D1234I(3,3)+2*(-D1234I(3,4)+D1234I(9,
     -   3)+D1234I(15,4)))*P(8)
       F(130)=DCMPLX(FR(130),FI(130))
       FR(131) = 4*(D1234R(5,2)-D1234R(6,2)-2*(D1234R(3,3)-D1234R(7,3)
     -   )+4*(D1234R(3,4)+D1234R(11,4)-D1234R(12,4)-D1234R(15,4)))*P(
     -   8)
       FI(131) = 4*(D1234I(5,2)-D1234I(6,2)-2*(D1234I(3,3)-D1234I(7,3)
     -   )+4*(D1234I(3,4)+D1234I(11,4)-D1234I(12,4)-D1234I(15,4)))*P(
     -   8)
       F(131)=DCMPLX(FR(131),FI(131))
       FR(132) = 8*(-D1234R(3,3)+D1234R(9,3)+2*(D1234R(3,4)+D1234R(14,
     -   4)-2*D1234R(15,4)))*P(8)
       FI(132) = 8*(-D1234I(3,3)+D1234I(9,3)+2*(D1234I(3,4)+D1234I(14,
     -   4)-2*D1234I(15,4)))*P(8)
       F(132)=DCMPLX(FR(132),FI(132))
       FR(133) = -8*(D1234R(3,2)+D1234R(9,3)+2*(-D1234R(3,4)+D1234R(15
     -   ,4)))*P(8)
       FI(133) = -8*(D1234I(3,2)+D1234I(9,3)+2*(-D1234I(3,4)+D1234I(15
     -   ,4)))*P(8)
       F(133)=DCMPLX(FR(133),FI(133))
       FR(134) = 8*(D1234R(3,3)-D1234R(9,3)-2*(-D1234R(3,4)+D1234R(15,
     -   4)))*P(8)
       FI(134) = 8*(D1234I(3,3)-D1234I(9,3)-2*(-D1234I(3,4)+D1234I(15,
     -   4)))*P(8)
       F(134)=DCMPLX(FR(134),FI(134))
       FR(135) = 16*(D1234R(3,3)+D1234R(3,4))*P(8)
       FI(135) = 16*(D1234I(3,3)+D1234I(3,4))*P(8)
       F(135)=DCMPLX(FR(135),FI(135))
       FR(136) = (B023R+4*(C134R(4,2)-C234R(4,2)+2*(C134R(5,3)-C134R(6
     -   ,3)-C234R(5,3)+C234R(6,3))))*P(8)+B013R*P(41)+M**2*P(9)*(D01
     -   234R*P(191)+(-D1234R(2,1)+D1234R(3,1))*P(192))+C0134R*P(193)
     -   +C0234R*P(194)+C134R(1,1)*P(195)+C234R(1,1)*P(196)+C134R(2,1
     -   )*P(197)+C234R(2,1)*P(198)
       FI(136) = (B023I+4*(C134I(4,2)-C234I(4,2)+2*(C134I(5,3)-C134I(6
     -   ,3)-C234I(5,3)+C234I(6,3))))*P(8)+B013I*P(41)+M**2*P(9)*(D01
     -   234I*P(191)+(-D1234I(2,1)+D1234I(3,1))*P(192))+C0134I*P(193)
     -   +C0234I*P(194)+C134I(1,1)*P(195)+C234I(1,1)*P(196)+C134I(2,1
     -   )*P(197)+C234I(2,1)*P(198)
       F(136)=DCMPLX(FR(136),FI(136))
       FR(137) = 4*((C234R(1,1)+C234R(1,2)-C234R(2,1)-C234R(2,2)+2*(C1
     -   34R(1,2)+C134R(1,3)+C134R(2,2)-C134R(2,3)-3*(C134R(3,3)-C134
     -   R(4,3))+C234R(2,3)+C234R(3,3)-2*(C134R(3,2)+C234R(4,3))))*P(
     -   8)-M**2*(D1234R(2,1)+D1234R(2,2)-D1234R(3,1)+D1234R(3,2)-2*D
     -   1234R(6,2))*P(9)*P(19)*P(64))
       FI(137) = 4*((C234I(1,1)+C234I(1,2)-C234I(2,1)-C234I(2,2)+2*(C1
     -   34I(1,2)+C134I(1,3)+C134I(2,2)-C134I(2,3)-3*(C134I(3,3)-C134
     -   I(4,3))+C234I(2,3)+C234I(3,3)-2*(C134I(3,2)+C234I(4,3))))*P(
     -   8)-M**2*(D1234I(2,1)+D1234I(2,2)-D1234I(3,1)+D1234I(3,2)-2*D
     -   1234I(6,2))*P(9)*P(19)*P(64))
       F(137)=DCMPLX(FR(137),FI(137))
       FR(138) = (B013R+8*(C134R(5,3)-C134R(6,3)-C234R(5,3)+C234R(6,3)
     -   ))*P(8)+B023R*P(41)+C234R(2,1)*P(132)+C234R(1,1)*P(135)-M**2
     -   *P(9)*(D01234R*P(191)+(D1234R(2,1)-D1234R(3,1))*P(199))+C013
     -   4R*P(200)+C0234R*P(201)+C134R(1,1)*P(203)+C134R(2,1)*P(204)
       FI(138) = (B013I+8*(C134I(5,3)-C134I(6,3)-C234I(5,3)+C234I(6,3)
     -   ))*P(8)+B023I*P(41)+C234I(2,1)*P(132)+C234I(1,1)*P(135)-M**2
     -   *P(9)*(D01234I*P(191)+(D1234I(2,1)-D1234I(3,1))*P(199))+C013
     -   4I*P(200)+C0234I*P(201)+C134I(1,1)*P(203)+C134I(2,1)*P(204)
       F(138)=DCMPLX(FR(138),FI(138))
       FR(139) = 4*((-4*C134R(3,2)+C234R(2,2)-C234R(3,2)+2*(C134R(1,2)
     -   +C134R(1,3)+C134R(2,2)-C134R(2,3)-3*(C134R(3,3)-C134R(4,3))+
     -   C234R(2,3)-C234R(4,3)))*P(8)+M**2*(-D1234R(2,1)+D1234R(3,1)-
     -   D1234R(3,2)-D1234R(4,2)+D1234R(5,2)+D1234R(6,2))*P(9)*P(25)*
     -   P(127))
       FI(139) = 4*((-4*C134I(3,2)+C234I(2,2)-C234I(3,2)+2*(C134I(1,2)
     -   +C134I(1,3)+C134I(2,2)-C134I(2,3)-3*(C134I(3,3)-C134I(4,3))+
     -   C234I(2,3)-C234I(4,3)))*P(8)+M**2*(-D1234I(2,1)+D1234I(3,1)-
     -   D1234I(3,2)-D1234I(4,2)+D1234I(5,2)+D1234I(6,2))*P(9)*P(25)*
     -   P(127))
       F(139)=DCMPLX(FR(139),FI(139))
       FR(140) = 4*((-4*C134R(3,2)+C234R(1,1)+C234R(1,2)-C234R(2,1)-C2
     -   34R(3,2)+2*(C134R(1,2)+C134R(1,3)+C134R(2,2)-C134R(2,3)-3*(C
     -   134R(3,3)-C134R(4,3))+C234R(2,3)+C234R(3,3)-2*C234R(4,3)))*P
     -   (8)+M**2*P(9)*(-((D1234R(2,2)+D1234R(3,2)-2*D1234R(6,2))*P(2
     -   )*P(10))+(D1234R(2,1)-D1234R(3,1))*P(205)))
       FI(140) = 4*((-4*C134I(3,2)+C234I(1,1)+C234I(1,2)-C234I(2,1)-C2
     -   34I(3,2)+2*(C134I(1,2)+C134I(1,3)+C134I(2,2)-C134I(2,3)-3*(C
     -   134I(3,3)-C134I(4,3))+C234I(2,3)+C234I(3,3)-2*C234I(4,3)))*P
     -   (8)+M**2*P(9)*(-((D1234I(2,2)+D1234I(3,2)-2*D1234I(6,2))*P(2
     -   )*P(10))+(D1234I(2,1)-D1234I(3,1))*P(205)))
       F(140)=DCMPLX(FR(140),FI(140))
       FR(141) = 4*((-4*C134R(3,2)-C234R(2,2)+C234R(3,2)+2*(C134R(1,2)
     -   +C134R(1,3)+C134R(2,2)-C134R(2,3)-3*(C134R(3,3)-C134R(4,3))+
     -   C234R(2,3)+C234R(3,3)-2*C234R(4,3)))*P(8)+M**2*P(9)*(-((D123
     -   4R(2,2)+D1234R(3,2)-2*D1234R(6,2))*P(1)*P(35))+(-D1234R(2,1)
     -   +D1234R(3,1))*P(206)))
       FI(141) = 4*((-4*C134I(3,2)-C234I(2,2)+C234I(3,2)+2*(C134I(1,2)
     -   +C134I(1,3)+C134I(2,2)-C134I(2,3)-3*(C134I(3,3)-C134I(4,3))+
     -   C234I(2,3)+C234I(3,3)-2*C234I(4,3)))*P(8)+M**2*P(9)*(-((D123
     -   4I(2,2)+D1234I(3,2)-2*D1234I(6,2))*P(1)*P(35))+(-D1234I(2,1)
     -   +D1234I(3,1))*P(206)))
       F(141)=DCMPLX(FR(141),FI(141))
       FR(142) = 8*(C134R(1,2)+C134R(1,3)+C134R(2,2)-C134R(2,3)-C234R(
     -   1,2)-C234R(1,3)-C234R(2,2)+C234R(2,3)-2*(C134R(3,2)-C234R(3,
     -   2))-3*(C134R(3,3)-C134R(4,3)-C234R(3,3)+C234R(4,3)))*P(8)
       FI(142) = 8*(C134I(1,2)+C134I(1,3)+C134I(2,2)-C134I(2,3)-C234I(
     -   1,2)-C234I(1,3)-C234I(2,2)+C234I(2,3)-2*(C134I(3,2)-C234I(3,
     -   2))-3*(C134I(3,3)-C134I(4,3)-C234I(3,3)+C234I(4,3)))*P(8)
       F(142)=DCMPLX(FR(142),FI(142))
       FR(143) = 4*(C134R(2,2)-C134R(3,2)-C234R(2,2)+C234R(3,2)-2*(C13
     -   4R(2,3)+C134R(3,3)-C234R(2,3)-C234R(3,3))+4*(C134R(4,3)-C234
     -   R(4,3)))*P(8)
       FI(143) = 4*(C134I(2,2)-C134I(3,2)-C234I(2,2)+C234I(3,2)-2*(C13
     -   4I(2,3)+C134I(3,3)-C234I(2,3)-C234I(3,3))+4*(C134I(4,3)-C234
     -   I(4,3)))*P(8)
       F(143)=DCMPLX(FR(143),FI(143))
       FR(144) = -4*(C134R(1,1)+C134R(1,2)-C134R(2,1)-C134R(3,2)-C234R
     -   (1,1)-C234R(1,2)+C234R(2,1)+C234R(3,2)+2*(C134R(2,3)+C134R(3
     -   ,3)-C234R(2,3)-C234R(3,3)-2*(C134R(4,3)-C234R(4,3))))*P(8)
       FI(144) = -4*(C134I(1,1)+C134I(1,2)-C134I(2,1)-C134I(3,2)-C234I
     -   (1,1)-C234I(1,2)+C234I(2,1)+C234I(3,2)+2*(C134I(2,3)+C134I(3
     -   ,3)-C234I(2,3)-C234I(3,3)-2*(C134I(4,3)-C234I(4,3))))*P(8)
       F(144)=DCMPLX(FR(144),FI(144))
       FR(145) = 4*(-C134R(2,2)+C134R(3,2)+C234R(2,2)-C234R(3,2)-2*(C1
     -   34R(2,3)-C134R(4,3)-C234R(2,3)+C234R(4,3)))*P(8)
       FI(145) = 4*(-C134I(2,2)+C134I(3,2)+C234I(2,2)-C234I(3,2)-2*(C1
     -   34I(2,3)-C134I(4,3)-C234I(2,3)+C234I(4,3)))*P(8)
       F(145)=DCMPLX(FR(145),FI(145))
       FR(146) = -4*(C134R(1,1)+C134R(1,2)-C134R(2,1)-C134R(2,2)-C234R
     -   (1,1)-C234R(1,2)+C234R(2,1)+C234R(2,2)+2*(C134R(2,3)+C134R(3
     -   ,3)-C234R(2,3)-C234R(3,3)-2*(C134R(4,3)-C234R(4,3))))*P(8)
       FI(146) = -4*(C134I(1,1)+C134I(1,2)-C134I(2,1)-C134I(2,2)-C234I
     -   (1,1)-C234I(1,2)+C234I(2,1)+C234I(2,2)+2*(C134I(2,3)+C134I(3
     -   ,3)-C234I(2,3)-C234I(3,3)-2*(C134I(4,3)-C234I(4,3))))*P(8)
       F(146)=DCMPLX(FR(146),FI(146))
       FR(147) = (B023R-2*(B13R(1,1)-B23R(1,1))+4*(C134R(4,2)-C234R(4,
     -   2))+8*(C134R(5,3)-C134R(6,3)-C234R(5,3)+C234R(6,3)))*P(8)+B0
     -   13R*P(41)+C0134R*P(193)+C0234R*P(194)+M**2*P(9)*(D01234R*P(1
     -   91)+(D1234R(2,1)-D1234R(3,1))*P(208))+C134R(1,1)*P(210)+C234
     -   R(1,1)*P(212)+C134R(2,1)*P(213)+C234R(2,1)*P(214)
       FI(147) = (B023I-2*(B13I(1,1)-B23I(1,1))+4*(C134I(4,2)-C234I(4,
     -   2))+8*(C134I(5,3)-C134I(6,3)-C234I(5,3)+C234I(6,3)))*P(8)+B0
     -   13I*P(41)+C0134I*P(193)+C0234I*P(194)+M**2*P(9)*(D01234I*P(1
     -   91)+(D1234I(2,1)-D1234I(3,1))*P(208))+C134I(1,1)*P(210)+C234
     -   I(1,1)*P(212)+C134I(2,1)*P(213)+C234I(2,1)*P(214)
       F(147)=DCMPLX(FR(147),FI(147))
       FR(148) = (B013R-2*B24R(1,1)+4*(C234R(4,2)+2*(C134R(5,3)-C134R(
     -   6,3)+C234R(6,3))))*P(8)+B024R*P(41)+C234R(2,1)*P(132)+C0134R
     -   *P(200)+C134R(1,1)*P(203)+C134R(2,1)*P(204)+M**2*P(9)*(C0123
     -   R*P(1)*P(35)+(C0124R-4*D1234R(7,2))*P(19)*P(64)+(-D1234R(1,1
     -   )+D1234R(3,1))*P(199)+D01234R*P(215))+C0234R*P(216)
       FI(148) = (B013I-2*B24I(1,1)+4*(C234I(4,2)+2*(C134I(5,3)-C134I(
     -   6,3)+C234I(6,3))))*P(8)+B024I*P(41)+C234I(2,1)*P(132)+C0134I
     -   *P(200)+C134I(1,1)*P(203)+C134I(2,1)*P(204)+M**2*P(9)*(C0123
     -   I*P(1)*P(35)+(C0124I-4*D1234I(7,2))*P(19)*P(64)+(-D1234I(1,1
     -   )+D1234I(3,1))*P(199)+D01234I*P(215))+C0234I*P(216)
       F(148)=DCMPLX(FR(148),FI(148))
       FR(149) = 4*(C134R(4,2)+C234R(4,2)+2*(C134R(5,3)-C134R(6,3)+C23
     -   4R(6,3)))*P(8)+(B013R+B024R)*P(41)+C0134R*P(193)+C134R(1,1)*
     -   P(195)+C134R(2,1)*P(197)+C234R(2,1)*P(198)-M**2*P(9)*(4*D123
     -   4R(7,2)*P(2)*P(10)+C0123R*P(1)*P(35)+C0124R*P(19)*P(64)+(D12
     -   34R(1,1)-D1234R(3,1))*P(192)+D01234R*P(215))+C0234R*P(216)
       FI(149) = 4*(C134I(4,2)+C234I(4,2)+2*(C134I(5,3)-C134I(6,3)+C23
     -   4I(6,3)))*P(8)+(B013I+B024I)*P(41)+C0134I*P(193)+C134I(1,1)*
     -   P(195)+C134I(2,1)*P(197)+C234I(2,1)*P(198)-M**2*P(9)*(4*D123
     -   4I(7,2)*P(2)*P(10)+C0123I*P(1)*P(35)+C0124I*P(19)*P(64)+(D12
     -   34I(1,1)-D1234I(3,1))*P(192)+D01234I*P(215))+C0234I*P(216)
       F(149)=DCMPLX(FR(149),FI(149))
       FR(150) = (B024R-2*B13R(1,1)+4*(C134R(4,2)+2*(C134R(5,3)-C134R(
     -   6,3)+C234R(6,3))))*P(8)+B013R*P(41)+C0134R*P(193)+C134R(1,1)
     -   *P(210)+C134R(2,1)*P(213)+C234R(2,1)*P(214)+M**2*P(9)*((C012
     -   3R-4*D1234R(7,2))*P(1)*P(35)+C0124R*P(19)*P(64)+(D1234R(1,1)
     -   -D1234R(3,1))*P(208)+D01234R*P(217))+C0234R*P(218)
       FI(150) = (B024I-2*B13I(1,1)+4*(C134I(4,2)+2*(C134I(5,3)-C134I(
     -   6,3)+C234I(6,3))))*P(8)+B013I*P(41)+C0134I*P(193)+C134I(1,1)
     -   *P(210)+C134I(2,1)*P(213)+C234I(2,1)*P(214)+M**2*P(9)*((C012
     -   3I-4*D1234I(7,2))*P(1)*P(35)+C0124I*P(19)*P(64)+(D1234I(1,1)
     -   -D1234I(3,1))*P(208)+D01234I*P(217))+C0234I*P(218)
       F(150)=DCMPLX(FR(150),FI(150))
       FR(151) = -2*(2*(-C134R(2,2)+C134R(3,2)-C234R(2,2)+C234R(3,2)+2
     -   *(C134R(2,3)+C134R(3,3)-2*C134R(4,3)-C234R(2,3)+C234R(4,3)))
     -   *P(8)+M**2*P(9)*(2*(D1234R(3,2)-D1234R(6,2))*P(2)*P(10)+D123
     -   4R(2,1)*P(1)*P(35)+D1234R(3,1)*P(33)*P(87)+D01234R*P(206)))
       FI(151) = -2*(2*(-C134I(2,2)+C134I(3,2)-C234I(2,2)+C234I(3,2)+2
     -   *(C134I(2,3)+C134I(3,3)-2*C134I(4,3)-C234I(2,3)+C234I(4,3)))
     -   *P(8)+M**2*P(9)*(2*(D1234I(3,2)-D1234I(6,2))*P(2)*P(10)+D123
     -   4I(2,1)*P(1)*P(35)+D1234I(3,1)*P(33)*P(87)+D01234I*P(206)))
       F(151)=DCMPLX(FR(151),FI(151))
       FR(152) = 2*(-2*(C134R(1,1)+C134R(1,2)-C134R(2,1)-C134R(3,2)-C2
     -   34R(2,2)+C234R(3,2)+2*(C134R(2,3)+C134R(3,3)-2*C134R(4,3)-C2
     -   34R(2,3)+C234R(4,3)))*P(8)+M**2*P(9)*((D1234R(2,1)-2*(D1234R
     -   (3,2)-D1234R(6,2)))*P(1)*P(35)+D1234R(3,1)*P(33)*P(87)+D0123
     -   4R*P(206)))
       FI(152) = 2*(-2*(C134I(1,1)+C134I(1,2)-C134I(2,1)-C134I(3,2)-C2
     -   34I(2,2)+C234I(3,2)+2*(C134I(2,3)+C134I(3,3)-2*C134I(4,3)-C2
     -   34I(2,3)+C234I(4,3)))*P(8)+M**2*P(9)*((D1234I(2,1)-2*(D1234I
     -   (3,2)-D1234I(6,2)))*P(1)*P(35)+D1234I(3,1)*P(33)*P(87)+D0123
     -   4I*P(206)))
       F(152)=DCMPLX(FR(152),FI(152))
       FR(153) = (B014R+2*(B14R(1,1)-B24R(1,1))-4*(C134R(4,2)-C234R(4,
     -   2))-8*(C134R(6,3)-C234R(6,3)))*P(8)+B024R*P(41)+C234R(2,1)*P
     -   (132)+C134R(2,1)*P(204)+C0234R*P(216)+M**2*P(9)*(D1234R(3,1)
     -   *P(199)+D01234R*P(219))+C0134R*P(220)
       FI(153) = (B014I+2*(B14I(1,1)-B24I(1,1))-4*(C134I(4,2)-C234I(4,
     -   2))-8*(C134I(6,3)-C234I(6,3)))*P(8)+B024I*P(41)+C234I(2,1)*P
     -   (132)+C134I(2,1)*P(204)+C0234I*P(216)+M**2*P(9)*(D1234I(3,1)
     -   *P(199)+D01234I*P(219))+C0134I*P(220)
       F(153)=DCMPLX(FR(153),FI(153))
       FR(154) = (B014R-4*(C134R(4,2)-C234R(4,2))-8*(C134R(6,3)-C234R(
     -   6,3)))*P(8)+B024R*P(41)+C134R(2,1)*P(197)+C234R(2,1)*P(198)+
     -   C0234R*P(216)+M**2*P(9)*(D1234R(3,1)*P(192)+D01234R*P(219))+
     -   C0134R*P(220)
       FI(154) = (B014I-4*(C134I(4,2)-C234I(4,2))-8*(C134I(6,3)-C234I(
     -   6,3)))*P(8)+B024I*P(41)+C134I(2,1)*P(197)+C234I(2,1)*P(198)+
     -   C0234I*P(216)+M**2*P(9)*(D1234I(3,1)*P(192)+D01234I*P(219))+
     -   C0134I*P(220)
       F(154)=DCMPLX(FR(154),FI(154))
       FR(155) = (B024R-8*(C134R(6,3)-C234R(6,3)))*P(8)+B014R*P(41)+C1
     -   34R(2,1)*P(213)+C234R(2,1)*P(214)+C0234R*P(218)-M**2*P(9)*(D
     -   1234R(3,1)*P(208)+D01234R*P(219))+C0134R*P(221)
       FI(155) = (B024I-8*(C134I(6,3)-C234I(6,3)))*P(8)+B014I*P(41)+C1
     -   34I(2,1)*P(213)+C234I(2,1)*P(214)+C0234I*P(218)-M**2*P(9)*(D
     -   1234I(3,1)*P(208)+D01234I*P(219))+C0134I*P(221)
       F(155)=DCMPLX(FR(155),FI(155))
       FR(156) = -4*((4*C134R(3,2)+6*(C134R(3,3)-C134R(4,3))+C234R(2,1
     -   )+C234R(3,2)-2*(C134R(1,2)+C134R(1,3)+C134R(2,2)-C134R(2,3)+
     -   C234R(2,3)-C234R(4,3)))*P(8)-M**2*P(9)*((D1234R(3,2)+D1234R(
     -   4,2))*P(55)-D1234R(2,1)*P(19)*P(64)+(D1234R(5,2)+D1234R(6,2)
     -   )*P(68)-D1234R(1,1)*P(206)+D1234R(3,1)*P(222)))
       FI(156) = -4*((4*C134I(3,2)+6*(C134I(3,3)-C134I(4,3))+C234I(2,1
     -   )+C234I(3,2)-2*(C134I(1,2)+C134I(1,3)+C134I(2,2)-C134I(2,3)+
     -   C234I(2,3)-C234I(4,3)))*P(8)-M**2*P(9)*((D1234I(3,2)+D1234I(
     -   4,2))*P(55)-D1234I(2,1)*P(19)*P(64)+(D1234I(5,2)+D1234I(6,2)
     -   )*P(68)-D1234I(1,1)*P(206)+D1234I(3,1)*P(222)))
       F(156)=DCMPLX(FR(156),FI(156))
       FR(157) = 4*((-C234R(2,1)+C234R(2,2)+2*(C134R(1,2)+C134R(1,3)+C
     -   134R(2,2)-C134R(2,3)-2*C134R(3,2)-3*(C134R(3,3)-C134R(4,3))+
     -   C234R(2,3)-C234R(3,2)-C234R(4,3)))*P(8)+M**2*P(9)*((-D1234R(
     -   1,1)+D1234R(3,1))*P(206)+(-D1234R(3,2)-D1234R(4,2)+D1234R(5,
     -   2)+D1234R(6,2))*P(60)*P(207)))
       FI(157) = 4*((-C234I(2,1)+C234I(2,2)+2*(C134I(1,2)+C134I(1,3)+C
     -   134I(2,2)-C134I(2,3)-2*C134I(3,2)-3*(C134I(3,3)-C134I(4,3))+
     -   C234I(2,3)-C234I(3,2)-C234I(4,3)))*P(8)+M**2*P(9)*((-D1234I(
     -   1,1)+D1234I(3,1))*P(206)+(-D1234I(3,2)-D1234I(4,2)+D1234I(5,
     -   2)+D1234I(6,2))*P(60)*P(207)))
       F(157)=DCMPLX(FR(157),FI(157))
       FR(158) = 2*(-2*(C134R(1,1)+C134R(1,2)-C134R(2,1)-C134R(2,2)+C2
     -   34R(2,1)-C234R(2,2)+2*(C134R(2,3)+C134R(3,3)-2*C134R(4,3)-C2
     -   34R(2,3)+C234R(3,2)+C234R(4,3)))*P(8)+M**2*P(9)*(2*(-D1234R(
     -   3,2)+D1234R(6,2))*P(2)*P(10)-D1234R(3,1)*P(33)*P(87)+D01234R
     -   *P(206)+D1234R(2,1)*P(60)*P(207)))
       FI(158) = 2*(-2*(C134I(1,1)+C134I(1,2)-C134I(2,1)-C134I(2,2)+C2
     -   34I(2,1)-C234I(2,2)+2*(C134I(2,3)+C134I(3,3)-2*C134I(4,3)-C2
     -   34I(2,3)+C234I(3,2)+C234I(4,3)))*P(8)+M**2*P(9)*(2*(-D1234I(
     -   3,2)+D1234I(6,2))*P(2)*P(10)-D1234I(3,1)*P(33)*P(87)+D01234I
     -   *P(206)+D1234I(2,1)*P(60)*P(207)))
       F(158)=DCMPLX(FR(158),FI(158))
       FR(159) = 2*(-2*(C134R(1,1)+C134R(1,2)-C134R(2,1)-C134R(2,2)+C2
     -   34R(2,1)+C234R(3,2)+2*(C134R(2,3)+C134R(3,3)-2*C134R(4,3)-C2
     -   34R(2,3)+C234R(4,3)))*P(8)+M**2*P(9)*(D1234R(3,1)*P(2)*P(10)
     -   +(D1234R(2,1)-2*(D1234R(3,2)-D1234R(6,2)))*P(1)*P(35)+D01234
     -   R*P(206)))
       FI(159) = 2*(-2*(C134I(1,1)+C134I(1,2)-C134I(2,1)-C134I(2,2)+C2
     -   34I(2,1)+C234I(3,2)+2*(C134I(2,3)+C134I(3,3)-2*C134I(4,3)-C2
     -   34I(2,3)+C234I(4,3)))*P(8)+M**2*P(9)*(D1234I(3,1)*P(2)*P(10)
     -   +(D1234I(2,1)-2*(D1234I(3,2)-D1234I(6,2)))*P(1)*P(35)+D01234
     -   I*P(206)))
       F(159)=DCMPLX(FR(159),FI(159))
       FR(160) = 2*(-2*(-C134R(2,2)+C134R(3,2)+C234R(2,1)+C234R(3,2)+2
     -   *(C134R(2,3)+C134R(3,3)-2*C134R(4,3)-C234R(2,3)+C234R(4,3)))
     -   *P(8)+M**2*P(9)*(D1234R(2,1)*P(1)*P(35)+2*(-D1234R(3,2)+D123
     -   4R(6,2))*P(19)*P(64)+D1234R(3,1)*P(33)*P(87)+D01234R*P(206))
     -   )
       FI(160) = 2*(-2*(-C134I(2,2)+C134I(3,2)+C234I(2,1)+C234I(3,2)+2
     -   *(C134I(2,3)+C134I(3,3)-2*C134I(4,3)-C234I(2,3)+C234I(4,3)))
     -   *P(8)+M**2*P(9)*(D1234I(2,1)*P(1)*P(35)+2*(-D1234I(3,2)+D123
     -   4I(6,2))*P(19)*P(64)+D1234I(3,1)*P(33)*P(87)+D01234I*P(206))
     -   )
       F(160)=DCMPLX(FR(160),FI(160))
       FR(161) = 2*(-2*(C134R(1,1)+C134R(1,2)-C134R(2,1)-C134R(3,2)+C2
     -   34R(2,1)-C234R(2,2)+2*(C134R(2,3)+C134R(3,3)-2*C134R(4,3)-C2
     -   34R(2,3)+C234R(3,2)+C234R(4,3)))*P(8)+M**2*P(9)*(-(D1234R(3,
     -   1)*P(2)*P(10))+2*(-D1234R(3,2)+D1234R(6,2))*P(19)*P(64)+D012
     -   34R*P(206)+D1234R(2,1)*P(60)*P(207)))
       FI(161) = 2*(-2*(C134I(1,1)+C134I(1,2)-C134I(2,1)-C134I(3,2)+C2
     -   34I(2,1)-C234I(2,2)+2*(C134I(2,3)+C134I(3,3)-2*C134I(4,3)-C2
     -   34I(2,3)+C234I(3,2)+C234I(4,3)))*P(8)+M**2*P(9)*(-(D1234I(3,
     -   1)*P(2)*P(10))+2*(-D1234I(3,2)+D1234I(6,2))*P(19)*P(64)+D012
     -   34I*P(206)+D1234I(2,1)*P(60)*P(207)))
       F(161)=DCMPLX(FR(161),FI(161))
       FR(162) = 4*(C134R(2,1)-C134R(2,2)-C234R(2,1)+C234R(2,2)+2*(-C1
     -   34R(2,3)+C134R(3,2)+C134R(4,3)+C234R(2,3)-C234R(3,2)-C234R(4
     -   ,3)))*P(8)
       FI(162) = 4*(C134I(2,1)-C134I(2,2)-C234I(2,1)+C234I(2,2)+2*(-C1
     -   34I(2,3)+C134I(3,2)+C134I(4,3)+C234I(2,3)-C234I(3,2)-C234I(4
     -   ,3)))*P(8)
       F(162)=DCMPLX(FR(162),FI(162))
       FR(163) = 4*(C134R(2,1)+C134R(3,2)-C234R(2,1)-C234R(3,2)-2*(C13
     -   4R(2,3)-C134R(4,3)-C234R(2,3)+C234R(4,3)))*P(8)
       FI(163) = 4*(C134I(2,1)+C134I(3,2)-C234I(2,1)-C234I(3,2)-2*(C13
     -   4I(2,3)-C134I(4,3)-C234I(2,3)+C234I(4,3)))*P(8)
       F(163)=DCMPLX(FR(163),FI(163))
       FR(164) = 4*((C134R(2,2)-C134R(3,2)+4*C134R(4,3)+2*(-C134R(2,3)
     -   -C134R(3,3)+C234R(2,2)+C234R(2,3)))*P(8)+M**2*(-D1234R(3,2)+
     -   D1234R(5,2))*P(9)*P(60)*P(207))
       FI(164) = 4*((C134I(2,2)-C134I(3,2)+4*C134I(4,3)+2*(-C134I(2,3)
     -   -C134I(3,3)+C234I(2,2)+C234I(2,3)))*P(8)+M**2*(-D1234I(3,2)+
     -   D1234I(5,2))*P(9)*P(60)*P(207))
       F(164)=DCMPLX(FR(164),FI(164))
       FR(165) = 8*((C134R(1,2)+C134R(1,3)+C134R(2,2)-C134R(2,3)-2*C13
     -   4R(3,2)-3*(C134R(3,3)-C134R(4,3))+C234R(2,2)+C234R(2,3))*P(8
     -   )-M**2*(D1234R(1,1)+D1234R(1,2)-D1234R(3,1)+D1234R(3,2)-2*D1
     -   234R(5,2))*P(9)*P(206))
       FI(165) = 8*((C134I(1,2)+C134I(1,3)+C134I(2,2)-C134I(2,3)-2*C13
     -   4I(3,2)-3*(C134I(3,3)-C134I(4,3))+C234I(2,2)+C234I(2,3))*P(8
     -   )-M**2*(D1234I(1,1)+D1234I(1,2)-D1234I(3,1)+D1234I(3,2)-2*D1
     -   234I(5,2))*P(9)*P(206))
       F(165)=DCMPLX(FR(165),FI(165))
       FR(166) = 4*((-C134R(1,1)-C134R(1,2)+C134R(2,1)+C134R(2,2)+2*(-
     -   C134R(2,3)-C134R(3,3)+2*C134R(4,3)+C234R(2,2)+C234R(2,3)))*P
     -   (8)+M**2*P(9)*((-D1234R(3,2)+D1234R(5,2))*P(25)*P(127)-D1234
     -   R(3,1)*P(205)+(D01234R+D1234R(1,1))*P(206)))
       FI(166) = 4*((-C134I(1,1)-C134I(1,2)+C134I(2,1)+C134I(2,2)+2*(-
     -   C134I(2,3)-C134I(3,3)+2*C134I(4,3)+C234I(2,2)+C234I(2,3)))*P
     -   (8)+M**2*P(9)*((-D1234I(3,2)+D1234I(5,2))*P(25)*P(127)-D1234
     -   I(3,1)*P(205)+(D01234I+D1234I(1,1))*P(206)))
       F(166)=DCMPLX(FR(166),FI(166))
       FR(167) = 4*((-C134R(1,1)-C134R(1,2)+C134R(2,1)+C134R(3,2)+4*C1
     -   34R(4,3)+2*(-C134R(2,3)-C134R(3,3)+C234R(2,2)+C234R(2,3)))*P
     -   (8)+M**2*P(9)*(D1234R(3,2)*P(55)+D1234R(5,2)*P(68)+D1234R(3,
     -   1)*P(205)+(D01234R+D1234R(1,1))*P(206)))
       FI(167) = 4*((-C134I(1,1)-C134I(1,2)+C134I(2,1)+C134I(3,2)+4*C1
     -   34I(4,3)+2*(-C134I(2,3)-C134I(3,3)+C234I(2,2)+C234I(2,3)))*P
     -   (8)+M**2*P(9)*(D1234I(3,2)*P(55)+D1234I(5,2)*P(68)+D1234I(3,
     -   1)*P(205)+(D01234I+D1234I(1,1))*P(206)))
       F(167)=DCMPLX(FR(167),FI(167))
       FR(168) = 4*((C134R(2,1)+C134R(3,2)+2*(-C134R(2,3)+C134R(4,3)+C
     -   234R(2,2)+C234R(2,3)))*P(8)-M**2*P(9)*(D1234R(3,2)*P(2)*P(10
     -   )+D1234R(3,1)*P(223)))
       FI(168) = 4*((C134I(2,1)+C134I(3,2)+2*(-C134I(2,3)+C134I(4,3)+C
     -   234I(2,2)+C234I(2,3)))*P(8)-M**2*P(9)*(D1234I(3,2)*P(2)*P(10
     -   )+D1234I(3,1)*P(223)))
       F(168)=DCMPLX(FR(168),FI(168))
       FR(169) = -4*(C134R(2,2)*P(8)-(C134R(3,2)+2*(-C134R(2,3)+C134R(
     -   4,3)+C234R(2,2)+C234R(2,3)))*P(8)+M**2*P(9)*(D1234R(3,2)*P(1
     -   9)*P(64)+D1234R(3,1)*P(206)))
       FI(169) = -4*(C134I(2,2)*P(8)-(C134I(3,2)+2*(-C134I(2,3)+C134I(
     -   4,3)+C234I(2,2)+C234I(2,3)))*P(8)+M**2*P(9)*(D1234I(3,2)*P(1
     -   9)*P(64)+D1234I(3,1)*P(206)))
       F(169)=DCMPLX(FR(169),FI(169))
       FR(170) = 4*((C134R(2,1)-C134R(2,2)+2*(-C134R(2,3)+C134R(3,2)+C
     -   134R(4,3)+C234R(2,2)+C234R(2,3)))*P(8)-M**2*(D1234R(3,1)+D12
     -   34R(3,2))*P(1)*P(9)*P(35))
       FI(170) = 4*((C134I(2,1)-C134I(2,2)+2*(-C134I(2,3)+C134I(3,2)+C
     -   134I(4,3)+C234I(2,2)+C234I(2,3)))*P(8)-M**2*(D1234I(3,1)+D12
     -   34I(3,2))*P(1)*P(9)*P(35))
       F(170)=DCMPLX(FR(170),FI(170))
       FR(171) = 8*(-C134R(2,2)-C134R(2,3)+C234R(2,2)+C234R(2,3))*P(8)
       FI(171) = 8*(-C134I(2,2)-C134I(2,3)+C234I(2,2)+C234I(2,3))*P(8)
       F(171)=DCMPLX(FR(171),FI(171))
       FR(172) = (B012R+B024R+8*(C124R(5,3)-C124R(6,3)-C134R(5,3)+C134
     -   R(6,3)))*P(8)+(B013R+B034R)*P(41)+M**2*P(1)*(D01234R*P(227)+
     -   (D1234R(1,1)-D1234R(3,1))*P(228))+C0124R*P(229)+C0134R*P(230
     -   )+C124R(1,1)*P(231)+C134R(1,1)*P(232)+C124R(2,1)*P(233)+C134
     -   R(2,1)*P(234)
       FI(172) = (B012I+B024I+8*(C124I(5,3)-C124I(6,3)-C134I(5,3)+C134
     -   I(6,3)))*P(8)+(B013I+B034I)*P(41)+M**2*P(1)*(D01234I*P(227)+
     -   (D1234I(1,1)-D1234I(3,1))*P(228))+C0124I*P(229)+C0134I*P(230
     -   )+C124I(1,1)*P(231)+C134I(1,1)*P(232)+C124I(2,1)*P(233)+C134
     -   I(2,1)*P(234)
       F(172)=DCMPLX(FR(172),FI(172))
       FR(173) = 8*(C124R(1,2)+C124R(1,3)+C124R(2,2)-C124R(2,3)-C134R(
     -   1,2)-C134R(1,3)-C134R(2,2)+C134R(2,3)-2*(C124R(3,2)-C134R(3,
     -   2))-3*(C124R(3,3)-C124R(4,3)-C134R(3,3)+C134R(4,3)))*P(8)
       FI(173) = 8*(C124I(1,2)+C124I(1,3)+C124I(2,2)-C124I(2,3)-C134I(
     -   1,2)-C134I(1,3)-C134I(2,2)+C134I(2,3)-2*(C124I(3,2)-C134I(3,
     -   2))-3*(C124I(3,3)-C124I(4,3)-C134I(3,3)+C134I(4,3)))*P(8)
       F(173)=DCMPLX(FR(173),FI(173))
       FR(174) = -8*(C134R(1,2)+C134R(1,3)+C134R(2,2)-C134R(2,3)-2*C13
     -   4R(3,2)-3*(C134R(3,3)-C134R(4,3)))*P(8)+4*M**2*P(1)*(-((D123
     -   4R(1,2)+D1234R(3,2)-2*D1234R(5,2))*P(9)*P(35))+(D1234R(1,1)-
     -   D1234R(3,1))*P(225))
       FI(174) = -8*(C134I(1,2)+C134I(1,3)+C134I(2,2)-C134I(2,3)-2*C13
     -   4I(3,2)-3*(C134I(3,3)-C134I(4,3)))*P(8)+4*M**2*P(1)*(-((D123
     -   4I(1,2)+D1234I(3,2)-2*D1234I(5,2))*P(9)*P(35))+(D1234I(1,1)-
     -   D1234I(3,1))*P(225))
       F(174)=DCMPLX(FR(174),FI(174))
       FR(175) = 4*(C124R(2,2)-C124R(3,2)-C134R(2,2)+C134R(3,2)-2*(C12
     -   4R(2,3)+C124R(3,3)-C134R(2,3)-C134R(3,3))+4*(C124R(4,3)-C134
     -   R(4,3)))*P(8)
       FI(175) = 4*(C124I(2,2)-C124I(3,2)-C134I(2,2)+C134I(3,2)-2*(C12
     -   4I(2,3)+C124I(3,3)-C134I(2,3)-C134I(3,3))+4*(C124I(4,3)-C134
     -   I(4,3)))*P(8)
       F(175)=DCMPLX(FR(175),FI(175))
       FR(176) = -8*(C134R(1,2)+C134R(1,3)+C134R(2,2)-C134R(2,3)-2*C13
     -   4R(3,2)-3*(C134R(3,3)-C134R(4,3)))*P(8)+4*M**2*P(1)*(-((D123
     -   4R(1,2)+D1234R(3,2)-2*D1234R(5,2))*P(2)*P(3))+(D1234R(1,1)-D
     -   1234R(3,1))*P(224))
       FI(176) = -8*(C134I(1,2)+C134I(1,3)+C134I(2,2)-C134I(2,3)-2*C13
     -   4I(3,2)-3*(C134I(3,3)-C134I(4,3)))*P(8)+4*M**2*P(1)*(-((D123
     -   4I(1,2)+D1234I(3,2)-2*D1234I(5,2))*P(2)*P(3))+(D1234I(1,1)-D
     -   1234I(3,1))*P(224))
       F(176)=DCMPLX(FR(176),FI(176))
       FR(177) = 4*(-2*(C134R(1,2)+C134R(1,3)+C134R(2,2)-C134R(2,3)-2*
     -   C134R(3,2)-3*(C134R(3,3)-C134R(4,3)))*P(8)+M**2*P(1)*((D1234
     -   R(2,1)-D1234R(3,1))*P(19)*P(49)+(D1234R(3,2)+D1234R(4,2))*P(
     -   154)+(D1234R(5,2)+D1234R(6,2))*P(161)))
       FI(177) = 4*(-2*(C134I(1,2)+C134I(1,3)+C134I(2,2)-C134I(2,3)-2*
     -   C134I(3,2)-3*(C134I(3,3)-C134I(4,3)))*P(8)+M**2*P(1)*((D1234
     -   I(2,1)-D1234I(3,1))*P(19)*P(49)+(D1234I(3,2)+D1234I(4,2))*P(
     -   154)+(D1234I(5,2)+D1234I(6,2))*P(161)))
       F(177)=DCMPLX(FR(177),FI(177))
       FR(178) = -4*(C124R(1,1)+C124R(1,2)-C124R(2,1)-C124R(3,2)-C134R
     -   (1,1)-C134R(1,2)+C134R(2,1)+C134R(3,2)+2*(C124R(2,3)+C124R(3
     -   ,3)-C134R(2,3)-C134R(3,3)-2*(C124R(4,3)-C134R(4,3))))*P(8)
       FI(178) = -4*(C124I(1,1)+C124I(1,2)-C124I(2,1)-C124I(3,2)-C134I
     -   (1,1)-C134I(1,2)+C134I(2,1)+C134I(3,2)+2*(C124I(2,3)+C124I(3
     -   ,3)-C134I(2,3)-C134I(3,3)-2*(C124I(4,3)-C134I(4,3))))*P(8)
       F(178)=DCMPLX(FR(178),FI(178))
       FR(179) = 4*(-C124R(2,2)+C124R(3,2)+C134R(2,2)-C134R(3,2)-2*(C1
     -   24R(2,3)-C124R(4,3)-C134R(2,3)+C134R(4,3)))*P(8)
       FI(179) = 4*(-C124I(2,2)+C124I(3,2)+C134I(2,2)-C134I(3,2)-2*(C1
     -   24I(2,3)-C124I(4,3)-C134I(2,3)+C134I(4,3)))*P(8)
       F(179)=DCMPLX(FR(179),FI(179))
       FR(180) = 8*(-C134R(5,3)+C134R(6,3))*P(8)+(B013R+B034R)*P(41)+C
     -   0134R*P(230)+C134R(1,1)*P(232)+C134R(2,1)*P(234)-M**2*P(1)*(
     -   C0234R*P(2)*P(3)+C0123R*P(9)*P(35)+4*D1234R(7,2)*P(19)*P(49)
     -   +(-D1234R(2,1)+D1234R(3,1))*P(228)+D01234R*P(235))
       FI(180) = 8*(-C134I(5,3)+C134I(6,3))*P(8)+(B013I+B034I)*P(41)+C
     -   0134I*P(230)+C134I(1,1)*P(232)+C134I(2,1)*P(234)-M**2*P(1)*(
     -   C0234I*P(2)*P(3)+C0123I*P(9)*P(35)+4*D1234I(7,2)*P(19)*P(49)
     -   +(-D1234I(2,1)+D1234I(3,1))*P(228)+D01234I*P(235))
       F(180)=DCMPLX(FR(180),FI(180))
       FR(181) = -8*(C134R(1,2)+C134R(1,3)+C134R(2,2)-C134R(2,3)-2*C13
     -   4R(3,2)-3*(C134R(3,3)-C134R(4,3)))*P(8)-4*M**2*(D1234R(1,1)+
     -   D1234R(1,2)-D1234R(3,1)+D1234R(3,2)-2*D1234R(5,2))*P(1)*P(19
     -   )*P(49)
       FI(181) = -8*(C134I(1,2)+C134I(1,3)+C134I(2,2)-C134I(2,3)-2*C13
     -   4I(3,2)-3*(C134I(3,3)-C134I(4,3)))*P(8)-4*M**2*(D1234I(1,1)+
     -   D1234I(1,2)-D1234I(3,1)+D1234I(3,2)-2*D1234I(5,2))*P(1)*P(19
     -   )*P(49)
       F(181)=DCMPLX(FR(181),FI(181))
       FR(182) = 4*(-2*(C134R(1,2)+C134R(1,3)+C134R(2,2)-C134R(2,3)-2*
     -   C134R(3,2)-3*(C134R(3,3)-C134R(4,3)))*P(8)+M**2*P(1)*((-D123
     -   4R(3,2)-D1234R(4,2)+D1234R(5,2)+D1234R(6,2))*P(33)*P(61)+(-D
     -   1234R(2,1)+D1234R(3,1))*P(224)))
       FI(182) = 4*(-2*(C134I(1,2)+C134I(1,3)+C134I(2,2)-C134I(2,3)-2*
     -   C134I(3,2)-3*(C134I(3,3)-C134I(4,3)))*P(8)+M**2*P(1)*((-D123
     -   4I(3,2)-D1234I(4,2)+D1234I(5,2)+D1234I(6,2))*P(33)*P(61)+(-D
     -   1234I(2,1)+D1234I(3,1))*P(224)))
       F(182)=DCMPLX(FR(182),FI(182))
       FR(183) = 4*(-2*(C134R(1,2)+C134R(1,3)+C134R(2,2)-C134R(2,3)-2*
     -   C134R(3,2)-3*(C134R(3,3)-C134R(4,3)))*P(8)+M**2*P(1)*((-D123
     -   4R(3,2)-D1234R(4,2)+D1234R(5,2)+D1234R(6,2))*P(86)*P(207)+(-
     -   D1234R(2,1)+D1234R(3,1))*P(225)))
       FI(183) = 4*(-2*(C134I(1,2)+C134I(1,3)+C134I(2,2)-C134I(2,3)-2*
     -   C134I(3,2)-3*(C134I(3,3)-C134I(4,3)))*P(8)+M**2*P(1)*((-D123
     -   4I(3,2)-D1234I(4,2)+D1234I(5,2)+D1234I(6,2))*P(86)*P(207)+(-
     -   D1234I(2,1)+D1234I(3,1))*P(225)))
       F(183)=DCMPLX(FR(183),FI(183))
       FR(184) = -8*((C134R(1,2)+C134R(1,3)+C134R(2,2)-C134R(2,3)-2*C1
     -   34R(3,2)-3*(C134R(3,3)-C134R(4,3)))*P(8)+M**2*(D1234R(2,2)+D
     -   1234R(3,2)-2*D1234R(6,2))*P(1)*P(236))
       FI(184) = -8*((C134I(1,2)+C134I(1,3)+C134I(2,2)-C134I(2,3)-2*C1
     -   34I(3,2)-3*(C134I(3,3)-C134I(4,3)))*P(8)+M**2*(D1234I(2,2)+D
     -   1234I(3,2)-2*D1234I(6,2))*P(1)*P(236))
       F(184)=DCMPLX(FR(184),FI(184))
       FR(185) = 4*((-C134R(2,2)+C134R(3,2)+2*(C134R(2,3)+C134R(3,3)-2
     -   *C134R(4,3)))*P(8)+M**2*(-D1234R(3,2)+D1234R(6,2))*P(1)*P(86
     -   )*P(207))
       FI(185) = 4*((-C134I(2,2)+C134I(3,2)+2*(C134I(2,3)+C134I(3,3)-2
     -   *C134I(4,3)))*P(8)+M**2*(-D1234I(3,2)+D1234I(6,2))*P(1)*P(86
     -   )*P(207))
       F(185)=DCMPLX(FR(185),FI(185))
       FR(186) = 4*((C134R(1,1)+C134R(1,2)-C134R(2,1)-C134R(3,2)+2*(C1
     -   34R(2,3)+C134R(3,3)-2*C134R(4,3)))*P(8)+M**2*P(1)*((-D1234R(
     -   3,2)+D1234R(6,2))*P(33)*P(61)+(D1234R(2,1)-D1234R(3,1))*P(23
     -   6)))
       FI(186) = 4*((C134I(1,1)+C134I(1,2)-C134I(2,1)-C134I(3,2)+2*(C1
     -   34I(2,3)+C134I(3,3)-2*C134I(4,3)))*P(8)+M**2*P(1)*((-D1234I(
     -   3,2)+D1234I(6,2))*P(33)*P(61)+(D1234I(2,1)-D1234I(3,1))*P(23
     -   6)))
       F(186)=DCMPLX(FR(186),FI(186))
       FR(187) = -4*(C124R(1,1)+C124R(1,2)-C124R(2,1)-C124R(2,2)-C134R
     -   (1,1)-C134R(1,2)+C134R(2,1)+C134R(2,2)+2*(C124R(2,3)+C124R(3
     -   ,3)-C134R(2,3)-C134R(3,3)-2*(C124R(4,3)-C134R(4,3))))*P(8)
       FI(187) = -4*(C124I(1,1)+C124I(1,2)-C124I(2,1)-C124I(2,2)-C134I
     -   (1,1)-C134I(1,2)+C134I(2,1)+C134I(2,2)+2*(C124I(2,3)+C124I(3
     -   ,3)-C134I(2,3)-C134I(3,3)-2*(C124I(4,3)-C134I(4,3))))*P(8)
       F(187)=DCMPLX(FR(187),FI(187))
       FR(188) = 4*((C134R(1,1)+C134R(1,2)-C134R(2,1)-C134R(2,2)+2*(C1
     -   34R(2,3)+C134R(3,3)-2*C134R(4,3)))*P(8)+M**2*P(1)*(D1234R(3,
     -   2)*P(154)+D1234R(6,2)*P(161)+(D1234R(2,1)-D1234R(3,1))*P(236
     -   )))
       FI(188) = 4*((C134I(1,1)+C134I(1,2)-C134I(2,1)-C134I(2,2)+2*(C1
     -   34I(2,3)+C134I(3,3)-2*C134I(4,3)))*P(8)+M**2*P(1)*(D1234I(3,
     -   2)*P(154)+D1234I(6,2)*P(161)+(D1234I(2,1)-D1234I(3,1))*P(236
     -   )))
       F(188)=DCMPLX(FR(188),FI(188))
       FR(189) = (B013R-2*(B34R(1,1)+2*C134R(4,2)+4*(C134R(5,3)-C134R(
     -   6,3))))*P(8)+B034R*P(41)+C0134R*P(229)+M**2*P(1)*((C0234R-4*
     -   D1234R(7,2))*P(2)*P(3)+C0123R*P(9)*P(35)+D01234R*P(235)+(D12
     -   34R(2,1)-D1234R(3,1))*P(237))+C134R(1,1)*P(239)+C134R(2,1)*P
     -   (240)
       FI(189) = (B013I-2*(B34I(1,1)+2*C134I(4,2)+4*(C134I(5,3)-C134I(
     -   6,3))))*P(8)+B034I*P(41)+C0134I*P(229)+M**2*P(1)*((C0234I-4*
     -   D1234I(7,2))*P(2)*P(3)+C0123I*P(9)*P(35)+D01234I*P(235)+(D12
     -   34I(2,1)-D1234I(3,1))*P(237))+C134I(1,1)*P(239)+C134I(2,1)*P
     -   (240)
       F(189)=DCMPLX(FR(189),FI(189))
       FR(190) = (B013R+B034R+2*B13R(1,1)-4*C134R(4,2)-8*(C134R(5,3)-C
     -   134R(6,3)))*P(8)+C134R(2,1)*P(210)+C134R(1,1)*P(213)+C0134R*
     -   P(229)+M**2*P(1)*(C0234R*P(2)*P(3)+(C0123R-4*D1234R(7,2))*P(
     -   9)*P(35)+D01234R*P(235)+(-D1234R(2,1)+D1234R(3,1))*P(241))
       FI(190) = (B013I+B034I+2*B13I(1,1)-4*C134I(4,2)-8*(C134I(5,3)-C
     -   134I(6,3)))*P(8)+C134I(2,1)*P(210)+C134I(1,1)*P(213)+C0134I*
     -   P(229)+M**2*P(1)*(C0234I*P(2)*P(3)+(C0123I-4*D1234I(7,2))*P(
     -   9)*P(35)+D01234I*P(235)+(-D1234I(2,1)+D1234I(3,1))*P(241))
       F(190)=DCMPLX(FR(190),FI(190))
       FR(191) = (B013R+B034R-2*(B12R(1,1)-B13R(1,1))+4*(C124R(4,2)-C1
     -   34R(4,2))+8*(C124R(5,3)-C124R(6,3)-C134R(5,3)+C134R(6,3)))*P
     -   (8)+(B012R+B024R)*P(41)+C134R(2,1)*P(210)+C134R(1,1)*P(213)+
     -   C0134R*P(229)+C0124R*P(230)+M**2*P(1)*((-D1234R(1,1)+D1234R(
     -   3,1))*P(241)+D01234R*P(242))+C124R(1,1)*P(244)+C124R(2,1)*P(
     -   245)
       FI(191) = (B013I+B034I-2*(B12I(1,1)-B13I(1,1))+4*(C124I(4,2)-C1
     -   34I(4,2))+8*(C124I(5,3)-C124I(6,3)-C134I(5,3)+C134I(6,3)))*P
     -   (8)+(B012I+B024I)*P(41)+C134I(2,1)*P(210)+C134I(1,1)*P(213)+
     -   C0134I*P(229)+C0124I*P(230)+M**2*P(1)*((-D1234I(1,1)+D1234I(
     -   3,1))*P(241)+D01234I*P(242))+C124I(1,1)*P(244)+C124I(2,1)*P(
     -   245)
       F(191)=DCMPLX(FR(191),FI(191))
       FR(192) = (B013R+B024R+2*(B24R(1,1)-B34R(1,1))+4*(C124R(4,2)-C1
     -   34R(4,2))+8*(C124R(5,3)-C124R(6,3)-C134R(5,3)+C134R(6,3)))*P
     -   (8)+(B012R+B034R)*P(41)+C124R(2,1)*P(108)+C124R(1,1)*P(109)+
     -   C0134R*P(229)+C0124R*P(230)+C134R(1,1)*P(239)+C134R(2,1)*P(2
     -   40)+M**2*P(1)*((D1234R(1,1)-D1234R(3,1))*P(237)+D01234R*P(24
     -   2))
       FI(192) = (B013I+B024I+2*(B24I(1,1)-B34I(1,1))+4*(C124I(4,2)-C1
     -   34I(4,2))+8*(C124I(5,3)-C124I(6,3)-C134I(5,3)+C134I(6,3)))*P
     -   (8)+(B012I+B034I)*P(41)+C124I(2,1)*P(108)+C124I(1,1)*P(109)+
     -   C0134I*P(229)+C0124I*P(230)+C134I(1,1)*P(239)+C134I(2,1)*P(2
     -   40)+M**2*P(1)*((D1234I(1,1)-D1234I(3,1))*P(237)+D01234I*P(24
     -   2))
       F(192)=DCMPLX(FR(192),FI(192))
       FR(193) = -2*(2*(C134R(2,2)-C134R(3,2)-2*(C134R(2,3)+C134R(3,3)
     -   -2*C134R(4,3)))*P(8)+M**2*P(1)*(2*(D1234R(3,2)-D1234R(5,2))*
     -   P(2)*P(3)-D1234R(1,1)*P(9)*P(35)+D1234R(3,1)*P(33)*P(61)+D01
     -   234R*P(225)))
       FI(193) = -2*(2*(C134I(2,2)-C134I(3,2)-2*(C134I(2,3)+C134I(3,3)
     -   -2*C134I(4,3)))*P(8)+M**2*P(1)*(2*(D1234I(3,2)-D1234I(5,2))*
     -   P(2)*P(3)-D1234I(1,1)*P(9)*P(35)+D1234I(3,1)*P(33)*P(61)+D01
     -   234I*P(225)))
       F(193)=DCMPLX(FR(193),FI(193))
       FR(194) = -2*(-2*(C134R(1,1)+C134R(1,2)-C134R(2,1)-C134R(3,2)+2
     -   *(C134R(2,3)+C134R(3,3)-2*C134R(4,3)))*P(8)+M**2*P(1)*(-((D1
     -   234R(1,1)-2*(D1234R(3,2)-D1234R(5,2)))*P(9)*P(35))+D1234R(3,
     -   1)*P(33)*P(61)+D01234R*P(225)))
       FI(194) = -2*(-2*(C134I(1,1)+C134I(1,2)-C134I(2,1)-C134I(3,2)+2
     -   *(C134I(2,3)+C134I(3,3)-2*C134I(4,3)))*P(8)+M**2*P(1)*(-((D1
     -   234I(1,1)-2*(D1234I(3,2)-D1234I(5,2)))*P(9)*P(35))+D1234I(3,
     -   1)*P(33)*P(61)+D01234I*P(225)))
       F(194)=DCMPLX(FR(194),FI(194))
       FR(195) = -2*(-2*(C134R(1,1)+C134R(1,2)-C134R(2,1)-C134R(2,2)+2
     -   *(C134R(2,3)+C134R(3,3)-2*C134R(4,3)))*P(8)+M**2*P(1)*(D1234
     -   R(3,1)*P(2)*P(3)-(D1234R(1,1)-2*(D1234R(3,2)-D1234R(5,2)))*P
     -   (9)*P(35)+D01234R*P(225)))
       FI(195) = -2*(-2*(C134I(1,1)+C134I(1,2)-C134I(2,1)-C134I(2,2)+2
     -   *(C134I(2,3)+C134I(3,3)-2*C134I(4,3)))*P(8)+M**2*P(1)*(D1234
     -   I(3,1)*P(2)*P(3)-(D1234I(1,1)-2*(D1234I(3,2)-D1234I(5,2)))*P
     -   (9)*P(35)+D01234I*P(225)))
       F(195)=DCMPLX(FR(195),FI(195))
       FR(196) = 2*(2*(C134R(1,1)+C134R(1,2)-C134R(2,1)-C134R(2,2)+2*(
     -   C134R(2,3)+C134R(3,3)-2*C134R(4,3)))*P(8)+M**2*P(1)*(2*(-D12
     -   34R(3,2)+D1234R(5,2))*P(2)*P(3)-D1234R(3,1)*P(33)*P(61)+D123
     -   4R(1,1)*P(86)*P(207)+D01234R*P(225)))
       FI(196) = 2*(2*(C134I(1,1)+C134I(1,2)-C134I(2,1)-C134I(2,2)+2*(
     -   C134I(2,3)+C134I(3,3)-2*C134I(4,3)))*P(8)+M**2*P(1)*(2*(-D12
     -   34I(3,2)+D1234I(5,2))*P(2)*P(3)-D1234I(3,1)*P(33)*P(61)+D123
     -   4I(1,1)*P(86)*P(207)+D01234I*P(225)))
       F(196)=DCMPLX(FR(196),FI(196))
       FR(197) = 2*(2*(-C134R(2,2)+C134R(3,2)+2*(C134R(2,3)+C134R(3,3)
     -   -2*C134R(4,3)))*P(8)+M**2*P(1)*(-(D1234R(1,1)*P(9)*P(35))+2*
     -   (-D1234R(3,2)+D1234R(5,2))*P(19)*P(49)+D1234R(3,1)*P(33)*P(6
     -   1)+D01234R*P(225)))
       FI(197) = 2*(2*(-C134I(2,2)+C134I(3,2)+2*(C134I(2,3)+C134I(3,3)
     -   -2*C134I(4,3)))*P(8)+M**2*P(1)*(-(D1234I(1,1)*P(9)*P(35))+2*
     -   (-D1234I(3,2)+D1234I(5,2))*P(19)*P(49)+D1234I(3,1)*P(33)*P(6
     -   1)+D01234I*P(225)))
       F(197)=DCMPLX(FR(197),FI(197))
       FR(198) = 2*(2*(C134R(1,1)+C134R(1,2)-C134R(2,1)-C134R(3,2)+2*(
     -   C134R(2,3)+C134R(3,3)-2*C134R(4,3)))*P(8)+M**2*P(1)*(-(D1234
     -   R(3,1)*P(2)*P(3))+2*(-D1234R(3,2)+D1234R(5,2))*P(19)*P(49)+D
     -   1234R(1,1)*P(86)*P(207)+D01234R*P(225)))
       FI(198) = 2*(2*(C134I(1,1)+C134I(1,2)-C134I(2,1)-C134I(3,2)+2*(
     -   C134I(2,3)+C134I(3,3)-2*C134I(4,3)))*P(8)+M**2*P(1)*(-(D1234
     -   I(3,1)*P(2)*P(3))+2*(-D1234I(3,2)+D1234I(5,2))*P(19)*P(49)+D
     -   1234I(1,1)*P(86)*P(207)+D01234I*P(225)))
       F(198)=DCMPLX(FR(198),FI(198))
       FR(199) = (B034R-8*(C124R(6,3)-C134R(6,3)))*P(8)+B024R*P(41)+C1
     -   34R(2,1)*P(210)+C124R(2,1)*P(245)+M**2*P(1)*(D1234R(3,1)*P(2
     -   41)+D01234R*P(246))+C0124R*P(247)+C0134R*P(248)
       FI(199) = (B034I-8*(C124I(6,3)-C134I(6,3)))*P(8)+B024I*P(41)+C1
     -   34I(2,1)*P(210)+C124I(2,1)*P(245)+M**2*P(1)*(D1234I(3,1)*P(2
     -   41)+D01234I*P(246))+C0124I*P(247)+C0134I*P(248)
       F(199)=DCMPLX(FR(199),FI(199))
       FR(200) = (B024R-4*(C124R(4,2)-C134R(4,2))-8*(C124R(6,3)-C134R(
     -   6,3)))*P(8)+B034R*P(41)+C124R(2,1)*P(233)+C134R(2,1)*P(234)-
     -   M**2*P(1)*(D1234R(3,1)*P(228)+D01234R*P(246))+C0124R*P(249)+
     -   C0134R*P(250)
       FI(200) = (B024I-4*(C124I(4,2)-C134I(4,2))-8*(C124I(6,3)-C134I(
     -   6,3)))*P(8)+B034I*P(41)+C124I(2,1)*P(233)+C134I(2,1)*P(234)-
     -   M**2*P(1)*(D1234I(3,1)*P(228)+D01234I*P(246))+C0124I*P(249)+
     -   C0134I*P(250)
       F(200)=DCMPLX(FR(200),FI(200))
       FR(201) = (B024R+2*(B24R(1,1)-B34R(1,1))-4*(C124R(4,2)-C134R(4,
     -   2))-8*(C124R(6,3)-C134R(6,3)))*P(8)+B034R*P(41)+C124R(2,1)*P
     -   (108)+C134R(2,1)*P(240)-M**2*P(1)*(D1234R(3,1)*P(237)+D01234
     -   R*P(246))+C0124R*P(249)+C0134R*P(250)
       FI(201) = (B024I+2*(B24I(1,1)-B34I(1,1))-4*(C124I(4,2)-C134I(4,
     -   2))-8*(C124I(6,3)-C134I(6,3)))*P(8)+B034I*P(41)+C124I(2,1)*P
     -   (108)+C134I(2,1)*P(240)-M**2*P(1)*(D1234I(3,1)*P(237)+D01234
     -   I*P(246))+C0124I*P(249)+C0134I*P(250)
       F(201)=DCMPLX(FR(201),FI(201))
       FR(202) = 4*(C124R(2,1)+C124R(3,2)-C134R(2,1)-C134R(3,2)-2*(C12
     -   4R(2,3)-C124R(4,3)-C134R(2,3)+C134R(4,3)))*P(8)
       FI(202) = 4*(C124I(2,1)+C124I(3,2)-C134I(2,1)-C134I(3,2)-2*(C12
     -   4I(2,3)-C124I(4,3)-C134I(2,3)+C134I(4,3)))*P(8)
       F(202)=DCMPLX(FR(202),FI(202))
       FR(203) = 4*(C124R(2,1)-C124R(2,2)-C134R(2,1)+C134R(2,2)+2*(-C1
     -   24R(2,3)+C124R(3,2)+C124R(4,3)+C134R(2,3)-C134R(3,2)-C134R(4
     -   ,3)))*P(8)
       FI(203) = 4*(C124I(2,1)-C124I(2,2)-C134I(2,1)+C134I(2,2)+2*(-C1
     -   24I(2,3)+C124I(3,2)+C124I(4,3)+C134I(2,3)-C134I(3,2)-C134I(4
     -   ,3)))*P(8)
       F(203)=DCMPLX(FR(203),FI(203))
       FR(204) = -4*((-C134R(2,2)+C134R(3,2)-2*(C134R(2,3)-C134R(4,3))
     -   )*P(8)+M**2*P(1)*(D1234R(3,2)*P(19)*P(49)+D1234R(3,1)*P(225)
     -   ))
       FI(204) = -4*((-C134I(2,2)+C134I(3,2)-2*(C134I(2,3)-C134I(4,3))
     -   )*P(8)+M**2*P(1)*(D1234I(3,2)*P(19)*P(49)+D1234I(3,1)*P(225)
     -   ))
       F(204)=DCMPLX(FR(204),FI(204))
       FR(205) = -4*((C134R(2,1)+C134R(3,2)-2*(C134R(2,3)-C134R(4,3)))
     -   *P(8)+M**2*P(1)*(D1234R(3,2)*P(2)*P(3)+D1234R(3,1)*P(236)))
       FI(205) = -4*((C134I(2,1)+C134I(3,2)-2*(C134I(2,3)-C134I(4,3)))
     -   *P(8)+M**2*P(1)*(D1234I(3,2)*P(2)*P(3)+D1234I(3,1)*P(236)))
       F(205)=DCMPLX(FR(205),FI(205))
       FR(206) = -4*((C134R(2,1)-C134R(2,2)+2*(-C134R(2,3)+C134R(3,2)+
     -   C134R(4,3)))*P(8)+M**2*(D1234R(3,1)+D1234R(3,2))*P(1)*P(9)*P
     -   (35))
       FI(206) = -4*((C134I(2,1)-C134I(2,2)+2*(-C134I(2,3)+C134I(3,2)+
     -   C134I(4,3)))*P(8)+M**2*(D1234I(3,1)+D1234I(3,2))*P(1)*P(9)*P
     -   (35))
       F(206)=DCMPLX(FR(206),FI(206))
       FR(207) = 8*(-C124R(2,2)-C124R(2,3)+C134R(2,2)+C134R(2,3))*P(8)
       FI(207) = 8*(-C124I(2,2)-C124I(2,3)+C134I(2,2)+C134I(2,3))*P(8)
       F(207)=DCMPLX(FR(207),FI(207))
       FR(208) = 4*(M**2*P(2)*((D1234R(1,1)+D1234R(1,2)-D1234R(3,1))*P
     -   (1)*P(3)+2*D1234R(3,2)*P(4)-D1234R(5,2)*P(7))+(C124R(1,1)+C1
     -   24R(1,2)-C124R(2,1)-C124R(3,2)+2*(C124R(2,3)+C124R(3,3)-2*C1
     -   24R(4,3)))*P(8))
       FI(208) = 4*(M**2*P(2)*((D1234I(1,1)+D1234I(1,2)-D1234I(3,1))*P
     -   (1)*P(3)+2*D1234I(3,2)*P(4)-D1234I(5,2)*P(7))+(C124I(1,1)+C1
     -   24I(1,2)-C124I(2,1)-C124I(3,2)+2*(C124I(2,3)+C124I(3,3)-2*C1
     -   24I(4,3)))*P(8))
       F(208)=DCMPLX(FR(208),FI(208))
       FR(209) = 4*M**2*P(2)*(2*D1234R(3,2)*P(4)+(D1234R(2,1)+D1234R(2
     -   ,2)-D1234R(3,1))*P(9)*P(10)-D1234R(6,2)*P(15))
       FI(209) = 4*M**2*P(2)*(2*D1234I(3,2)*P(4)+(D1234I(2,1)+D1234I(2
     -   ,2)-D1234I(3,1))*P(9)*P(10)-D1234I(6,2)*P(15))
       F(209)=DCMPLX(FR(209),FI(209))
       FR(210) = 4*M**2*P(2)*(D1234R(2,2)*P(1)*P(3)+2*D1234R(3,2)*P(4)
     -   -D1234R(6,2)*P(7)+(-D1234R(2,1)+D1234R(3,1))*P(16))
       FI(210) = 4*M**2*P(2)*(D1234I(2,2)*P(1)*P(3)+2*D1234I(3,2)*P(4)
     -   -D1234I(6,2)*P(7)+(-D1234I(2,1)+D1234I(3,1))*P(16))
       F(210)=DCMPLX(FR(210),FI(210))
       FR(211) = 4*M**2*P(2)*(2*D1234R(3,2)*P(4)+(D1234R(2,1)-D1234R(3
     -   ,1))*P(17)+D1234R(2,2)*P(18)*P(19)-D1234R(6,2)*P(22))
       FI(211) = 4*M**2*P(2)*(2*D1234I(3,2)*P(4)+(D1234I(2,1)-D1234I(3
     -   ,1))*P(17)+D1234I(2,2)*P(18)*P(19)-D1234I(6,2)*P(22))
       F(211)=DCMPLX(FR(211),FI(211))
       FR(212) = -2*(B13R(1,1)-2*(C123R(4,2)+2*(C123R(5,3)-C124R(5,3)+
     -   C124R(6,3))))*P(8)-M**2*P(2)*(C0234R*P(1)*P(3)+(C0134R-4*D12
     -   34R(7,2))*P(18)*P(19)+D01234R*P(23)+(-D1234R(1,1)+D1234R(3,1
     -   ))*P(26))+C0123R*P(32)+C0124R*P(34)+C123R(1,1)*P(37)+C124R(1
     -   ,1)*P(39)+C124R(2,1)*P(40)+(B013R+B024R)*P(41)
       FI(212) = -2*(B13I(1,1)-2*(C123I(4,2)+2*(C123I(5,3)-C124I(5,3)+
     -   C124I(6,3))))*P(8)-M**2*P(2)*(C0234I*P(1)*P(3)+(C0134I-4*D12
     -   34I(7,2))*P(18)*P(19)+D01234I*P(23)+(-D1234I(1,1)+D1234I(3,1
     -   ))*P(26))+C0123I*P(32)+C0124I*P(34)+C123I(1,1)*P(37)+C124I(1
     -   ,1)*P(39)+C124I(2,1)*P(40)+(B013I+B024I)*P(41)
       F(212)=DCMPLX(FR(212),FI(212))
       FR(213) = (B013R+B024R+4*(C123R(4,2)-C124R(4,2)+2*(C123R(5,3)-C
     -   124R(5,3)+C124R(6,3))))*P(8)+M**2*P(2)*(C0234R*P(1)*P(3)+4*D
     -   1234R(7,2)*P(9)*P(10)+C0134R*P(18)*P(19)+D01234R*P(23)+(D123
     -   4R(1,1)-D1234R(3,1))*P(42))+C0123R*P(47)+C0124R*P(48)+C123R(
     -   1,1)*P(51)+C124R(1,1)*P(52)+C124R(2,1)*P(53)
       FI(213) = (B013I+B024I+4*(C123I(4,2)-C124I(4,2)+2*(C123I(5,3)-C
     -   124I(5,3)+C124I(6,3))))*P(8)+M**2*P(2)*(C0234I*P(1)*P(3)+4*D
     -   1234I(7,2)*P(9)*P(10)+C0134I*P(18)*P(19)+D01234I*P(23)+(D123
     -   4I(1,1)-D1234I(3,1))*P(42))+C0123I*P(47)+C0124I*P(48)+C123I(
     -   1,1)*P(51)+C124I(1,1)*P(52)+C124I(2,1)*P(53)
       F(213)=DCMPLX(FR(213),FI(213))
       FR(214) = C123R(2,1)*P(37)+(B013R+B024R)*P(41)-M**2*P(2)*((C013
     -   4R-4*D1234R(7,2))*P(18)*P(19)+(-D1234R(2,1)+D1234R(3,1))*P(2
     -   6)+D01234R*P(54))+C0123R*P(56)+C0234R*P(58)+P(8)*(4*C123R(4,
     -   2)+8*C123R(6,3)-2*(B13R(1,1)+p3sq*C234R(2,1))+C234R(1,1)*P(5
     -   9))
       FI(214) = C123I(2,1)*P(37)+(B013I+B024I)*P(41)-M**2*P(2)*((C013
     -   4I-4*D1234I(7,2))*P(18)*P(19)+(-D1234I(2,1)+D1234I(3,1))*P(2
     -   6)+D01234I*P(54))+C0123I*P(56)+C0234I*P(58)+P(8)*(4*C123I(4,
     -   2)+8*C123I(6,3)-2*(B13I(1,1)+p3sq*C234I(2,1))+C234I(1,1)*P(5
     -   9))
       F(214)=DCMPLX(FR(214),FI(214))
       FR(215) = (B013R+B024R)*P(41)+C0123R*P(56)+C0234R*P(58)+M**2*P(
     -   2)*(4*D1234R(7,2)*P(1)*P(3)-C0134R*P(18)*P(19)-D01234R*P(54)
     -   +(-D1234R(2,1)+D1234R(3,1))*P(62))+C123R(2,1)*P(66)+P(8)*(8*
     -   C123R(6,3)+4*(C123R(4,2)+p3sq*C234R(2,2)+C234R(4,2))+C234R(1
     -   ,1)*P(63)+2*((C234R(1,2)-C234R(2,1))*P(59)-C234R(3,2)*P(67))
     -   )
       FI(215) = (B013I+B024I)*P(41)+C0123I*P(56)+C0234I*P(58)+M**2*P(
     -   2)*(4*D1234I(7,2)*P(1)*P(3)-C0134I*P(18)*P(19)-D01234I*P(54)
     -   +(-D1234I(2,1)+D1234I(3,1))*P(62))+C123I(2,1)*P(66)+P(8)*(8*
     -   C123I(6,3)+4*(C123I(4,2)+p3sq*C234I(2,2)+C234I(4,2))+C234I(1
     -   ,1)*P(63)+2*((C234I(1,2)-C234I(2,1))*P(59)-C234I(3,2)*P(67))
     -   )
       F(215)=DCMPLX(FR(215),FI(215))
       FR(216) = C123R(2,1)*P(51)+M**2*P(2)*(4*D1234R(7,2)*P(9)*P(10)+
     -   C0134R*P(18)*P(19)+(D1234R(2,1)-D1234R(3,1))*P(42)+D01234R*P
     -   (54))+P(8)*(B013R+B024R+8*C123R(6,3)+2*p3sq*C234R(2,1)-C234R
     -   (1,1)*P(59))+C0123R*P(69)+C0234R*P(70)
       FI(216) = C123I(2,1)*P(51)+M**2*P(2)*(4*D1234I(7,2)*P(9)*P(10)+
     -   C0134I*P(18)*P(19)+(D1234I(2,1)-D1234I(3,1))*P(42)+D01234I*P
     -   (54))+P(8)*(B013I+B024I+8*C123I(6,3)+2*p3sq*C234I(2,1)-C234I
     -   (1,1)*P(59))+C0123I*P(69)+C0234I*P(70)
       F(216)=DCMPLX(FR(216),FI(216))
       FR(217) = 4*((-C124R(2,2)+C124R(3,2)+2*(C124R(2,3)+C124R(3,3)-2
     -   *C124R(4,3)))*P(8)+M**2*P(2)*(2*D1234R(3,2)*P(4)+D1234R(1,2)
     -   *P(9)*P(10)-D1234R(5,2)*P(15)-D01234R*P(17)+D1234R(1,1)*P(71
     -   )-D1234R(3,1)*P(74)))
       FI(217) = 4*((-C124I(2,2)+C124I(3,2)+2*(C124I(2,3)+C124I(3,3)-2
     -   *C124I(4,3)))*P(8)+M**2*P(2)*(2*D1234I(3,2)*P(4)+D1234I(1,2)
     -   *P(9)*P(10)-D1234I(5,2)*P(15)-D01234I*P(17)+D1234I(1,1)*P(71
     -   )-D1234I(3,1)*P(74)))
       F(217)=DCMPLX(FR(217),FI(217))
       FR(218) = 4*((C124R(1,1)+C124R(1,2)-C124R(2,1)-C124R(2,2)+2*(C1
     -   24R(2,3)+C124R(3,3)-2*C124R(4,3)))*P(8)+M**2*P(2)*(2*D1234R(
     -   3,2)*P(4)+D01234R*P(17)+D1234R(1,2)*P(18)*P(19)-D1234R(5,2)*
     -   P(22)+D1234R(1,1)*P(76)-D1234R(3,1)*P(77)))
       FI(218) = 4*((C124I(1,1)+C124I(1,2)-C124I(2,1)-C124I(2,2)+2*(C1
     -   24I(2,3)+C124I(3,3)-2*C124I(4,3)))*P(8)+M**2*P(2)*(2*D1234I(
     -   3,2)*P(4)+D01234I*P(17)+D1234I(1,2)*P(18)*P(19)-D1234I(5,2)*
     -   P(22)+D1234I(1,1)*P(76)-D1234I(3,1)*P(77)))
       F(218)=DCMPLX(FR(218),FI(218))
       FR(219) = 2*M**2*P(2)*(4*D1234R(3,2)*P(4)+(D1234R(1,1)-2*(D1234
     -   R(3,1)-D1234R(4,2)))*P(9)*P(10)-D01234R*P(17)-2*(D1234R(5,2)
     -   *P(24)*P(25)+D1234R(6,2)*P(60)*P(61))+D1234R(2,1)*P(81))
       FI(219) = 2*M**2*P(2)*(4*D1234I(3,2)*P(4)+(D1234I(1,1)-2*(D1234
     -   I(3,1)-D1234I(4,2)))*P(9)*P(10)-D01234I*P(17)-2*(D1234I(5,2)
     -   *P(24)*P(25)+D1234I(6,2)*P(60)*P(61))+D1234I(2,1)*P(81))
       F(219)=DCMPLX(FR(219),FI(219))
       FR(220) = 2*M**2*P(2)*(4*D1234R(3,2)*P(4)+(D1234R(1,1)+2*D1234R
     -   (4,2))*P(9)*P(10)-D01234R*P(17)-2*(D1234R(6,2)*P(24)*P(25)+D
     -   1234R(5,2)*P(60)*P(61)+D1234R(3,1)*P(74))+D1234R(2,1)*P(85))
       FI(220) = 2*M**2*P(2)*(4*D1234I(3,2)*P(4)+(D1234I(1,1)+2*D1234I
     -   (4,2))*P(9)*P(10)-D01234I*P(17)-2*(D1234I(6,2)*P(24)*P(25)+D
     -   1234I(5,2)*P(60)*P(61)+D1234I(3,1)*P(74))+D1234I(2,1)*P(85))
       F(220)=DCMPLX(FR(220),FI(220))
       FR(221) = -2*M**2*P(2)*(-((D1234R(2,1)-2*(D1234R(3,1)-D1234R(4,
     -   2)))*P(1)*P(3))-4*D1234R(3,2)*P(4)+D01234R*P(17)+D1234R(1,1)
     -   *P(86)*P(87)+2*(D1234R(6,2)*P(24)*P(25)+D1234R(5,2)*P(91)))
       FI(221) = -2*M**2*P(2)*(-((D1234I(2,1)-2*(D1234I(3,1)-D1234I(4,
     -   2)))*P(1)*P(3))-4*D1234I(3,2)*P(4)+D01234I*P(17)+D1234I(1,1)
     -   *P(86)*P(87)+2*(D1234I(6,2)*P(24)*P(25)+D1234I(5,2)*P(91)))
       F(221)=DCMPLX(FR(221),FI(221))
       FR(222) = 2*M**2*P(2)*(4*D1234R(3,2)*P(4)+D01234R*P(17)+D1234R(
     -   1,1)*P(86)*P(87)+2*((-D1234R(3,1)+D1234R(4,2))*P(18)*P(19)-D
     -   1234R(6,2)*P(60)*P(61)-D1234R(5,2)*P(91))+D1234R(2,1)*P(94))
       FI(222) = 2*M**2*P(2)*(4*D1234I(3,2)*P(4)+D01234I*P(17)+D1234I(
     -   1,1)*P(86)*P(87)+2*((-D1234I(3,1)+D1234I(4,2))*P(18)*P(19)-D
     -   1234I(6,2)*P(60)*P(61)-D1234I(5,2)*P(91))+D1234I(2,1)*P(94))
       F(222)=DCMPLX(FR(222),FI(222))
       FR(223) = 2*M**2*P(2)*(4*D1234R(3,2)*P(4)-D1234R(1,1)*P(9)*P(10
     -   )+D01234R*P(17)-2*((D1234R(3,1)-D1234R(4,2))*P(18)*P(19)+D12
     -   34R(5,2)*P(60)*P(61)+D1234R(6,2)*P(91))+D1234R(2,1)*P(97))
       FI(223) = 2*M**2*P(2)*(4*D1234I(3,2)*P(4)-D1234I(1,1)*P(9)*P(10
     -   )+D01234I*P(17)-2*((D1234I(3,1)-D1234I(4,2))*P(18)*P(19)+D12
     -   34I(5,2)*P(60)*P(61)+D1234I(6,2)*P(91))+D1234I(2,1)*P(97))
       F(223)=DCMPLX(FR(223),FI(223))
       FR(224) = 2*M**2*P(2)*((D1234R(2,1)+2*D1234R(4,2))*P(1)*P(3)+4*
     -   D1234R(3,2)*P(4)-D1234R(1,1)*P(9)*P(10)+D01234R*P(17)+2*(D12
     -   34R(3,1)*P(16)-D1234R(5,2)*P(24)*P(25)-D1234R(6,2)*P(91)))
       FI(224) = 2*M**2*P(2)*((D1234I(2,1)+2*D1234I(4,2))*P(1)*P(3)+4*
     -   D1234I(3,2)*P(4)-D1234I(1,1)*P(9)*P(10)+D01234I*P(17)+2*(D12
     -   34I(3,1)*P(16)-D1234I(5,2)*P(24)*P(25)-D1234I(6,2)*P(91)))
       F(224)=DCMPLX(FR(224),FI(224))
       FR(225) = (B014R+B024R+8*C124R(6,3))*P(8)+C124R(2,1)*P(53)+C012
     -   4R*P(69)+M**2*P(2)*(C0234R*P(1)*P(3)+4*D1234R(7,2)*P(9)*P(10
     -   )+C0134R*P(18)*P(19)-D1234R(3,1)*P(42)+D01234R*P(98))
       FI(225) = (B014I+B024I+8*C124I(6,3))*P(8)+C124I(2,1)*P(53)+C012
     -   4I*P(69)+M**2*P(2)*(C0234I*P(1)*P(3)+4*D1234I(7,2)*P(9)*P(10
     -   )+C0134I*P(18)*P(19)-D1234I(3,1)*P(42)+D01234I*P(98))
       F(225)=DCMPLX(FR(225),FI(225))
       FR(226) = -2*(B14R(1,1)-2*(C124R(4,2)+2*C124R(6,3)))*P(8)+C124R
     -   (2,1)*P(40)+(B014R+B024R)*P(41)+C0124R*P(56)-M**2*P(2)*(C023
     -   4R*P(1)*P(3)+(C0134R-4*D1234R(7,2))*P(18)*P(19)+D1234R(3,1)*
     -   P(26)+D01234R*P(98))
       FI(226) = -2*(B14I(1,1)-2*(C124I(4,2)+2*C124I(6,3)))*P(8)+C124I
     -   (2,1)*P(40)+(B014I+B024I)*P(41)+C0124I*P(56)-M**2*P(2)*(C023
     -   4I*P(1)*P(3)+(C0134I-4*D1234I(7,2))*P(18)*P(19)+D1234I(3,1)*
     -   P(26)+D01234I*P(98))
       F(226)=DCMPLX(FR(226),FI(226))
       FR(227) = (B013R+B024R)*P(41)+C0124R*P(48)+2*P(8)*(4*(C123R(4,2
     -   )+C123R(5,3)-C124R(5,3)+C124R(6,3))+2*(-C124R(4,2)+p3sq*C234
     -   R(2,2)+C234R(4,2))-C234R(2,1)*P(57)-C234R(3,2)*P(59))+C123R(
     -   1,1)*P(66)-M**2*P(2)*((C0234R-4*D1234R(7,2))*P(1)*P(3)+C0134
     -   R*P(18)*P(19)+(D1234R(1,1)-D1234R(3,1))*P(62)+D01234R*P(99))
     -   +C0123R*P(106)+C124R(1,1)*P(108)+C124R(2,1)*P(109)
       FI(227) = (B013I+B024I)*P(41)+C0124I*P(48)+2*P(8)*(4*(C123I(4,2
     -   )+C123I(5,3)-C124I(5,3)+C124I(6,3))+2*(-C124I(4,2)+p3sq*C234
     -   I(2,2)+C234I(4,2))-C234I(2,1)*P(57)-C234I(3,2)*P(59))+C123I(
     -   1,1)*P(66)-M**2*P(2)*((C0234I-4*D1234I(7,2))*P(1)*P(3)+C0134
     -   I*P(18)*P(19)+(D1234I(1,1)-D1234I(3,1))*P(62)+D01234I*P(99))
     -   +C0123I*P(106)+C124I(1,1)*P(108)+C124I(2,1)*P(109)
       F(227)=DCMPLX(FR(227),FI(227))
       FR(228) = (B014R+B024R)*P(41)+C0124R*P(56)+2*P(8)*(4*C124R(6,3)
     -   +2*(C124R(4,2)+p3sq*C234R(2,2)+C234R(4,2))-C234R(2,1)*P(57)-
     -   C234R(3,2)*P(59))-M**2*P(2)*((C0234R-4*D1234R(7,2))*P(1)*P(3
     -   )+C0134R*P(18)*P(19)-D1234R(3,1)*P(62)+D01234R*P(98))+C124R(
     -   2,1)*P(109)
       FI(228) = (B014I+B024I)*P(41)+C0124I*P(56)+2*P(8)*(4*C124I(6,3)
     -   +2*(C124I(4,2)+p3sq*C234I(2,2)+C234I(4,2))-C234I(2,1)*P(57)-
     -   C234I(3,2)*P(59))-M**2*P(2)*((C0234I-4*D1234I(7,2))*P(1)*P(3
     -   )+C0134I*P(18)*P(19)-D1234I(3,1)*P(62)+D01234I*P(98))+C124I(
     -   2,1)*P(109)
       F(228)=DCMPLX(FR(228),FI(228))
       FR(229) = 4*M**2*P(2)*(2*D1234R(3,2)*P(4)+D1234R(3,1)*P(9)*P(10
     -   )-D1234R(6,2)*P(91))
       FI(229) = 4*M**2*P(2)*(2*D1234I(3,2)*P(4)+D1234I(3,1)*P(9)*P(10
     -   )-D1234I(6,2)*P(91))
       F(229)=DCMPLX(FR(229),FI(229))
       FR(230) = -4*((C124R(2,1)-C124R(2,2)+2*(-C124R(2,3)+C124R(3,2)+
     -   C124R(4,3)))*P(8)+M**2*P(2)*(-2*D1234R(3,2)*P(4)+D1234R(5,2)
     -   *P(91)+D1234R(3,1)*P(113)))
       FI(230) = -4*((C124I(2,1)-C124I(2,2)+2*(-C124I(2,3)+C124I(3,2)+
     -   C124I(4,3)))*P(8)+M**2*P(2)*(-2*D1234I(3,2)*P(4)+D1234I(5,2)
     -   *P(91)+D1234I(3,1)*P(113)))
       F(230)=DCMPLX(FR(230),FI(230))
       FR(231) = 4*M**2*P(2)*(2*D1234R(3,2)*P(4)+D1234R(3,1)*P(17)-D12
     -   34R(6,2)*P(24)*P(25))
       FI(231) = 4*M**2*P(2)*(2*D1234I(3,2)*P(4)+D1234I(3,1)*P(17)-D12
     -   34I(6,2)*P(24)*P(25))
       F(231)=DCMPLX(FR(231),FI(231))
       FR(232) = 4*((C124R(2,2)-C124R(3,2)+2*(C124R(2,3)-C124R(4,3)))*
     -   P(8)+M**2*P(2)*(2*D1234R(3,2)*P(4)+D1234R(3,1)*P(17)-D1234R(
     -   5,2)*P(24)*P(25)))
       FI(232) = 4*((C124I(2,2)-C124I(3,2)+2*(C124I(2,3)-C124I(4,3)))*
     -   P(8)+M**2*P(2)*(2*D1234I(3,2)*P(4)+D1234I(3,1)*P(17)-D1234I(
     -   5,2)*P(24)*P(25)))
       F(232)=DCMPLX(FR(232),FI(232))
       FR(233) = 4*M**2*P(2)*(2*D1234R(3,2)*P(4)-D1234R(3,1)*P(16)-D12
     -   34R(6,2)*P(60)*P(61))
       FI(233) = 4*M**2*P(2)*(2*D1234I(3,2)*P(4)-D1234I(3,1)*P(16)-D12
     -   34I(6,2)*P(60)*P(61))
       F(233)=DCMPLX(FR(233),FI(233))
       FR(234) = 4*(-((C124R(2,1)+C124R(3,2)-2*(C124R(2,3)-C124R(4,3))
     -   )*P(8))+M**2*P(2)*(2*D1234R(3,2)*P(4)-(D1234R(3,1)+D1234R(5,
     -   2))*P(60)*P(61)))
       FI(234) = 4*(-((C124I(2,1)+C124I(3,2)-2*(C124I(2,3)-C124I(4,3))
     -   )*P(8))+M**2*P(2)*(2*D1234I(3,2)*P(4)-(D1234I(3,1)+D1234I(5,
     -   2))*P(60)*P(61)))
       F(234)=DCMPLX(FR(234),FI(234))
       FR(235) = 8*(M**2*(D1234R(1,1)+D1234R(1,2)-D1234R(3,1)+D1234R(3
     -   ,2)-2*D1234R(5,2))*P(2)*P(4)+(C123R(1,1)+C123R(1,3)-C124R(1,
     -   2)-C124R(1,3)-C124R(2,2)+C124R(2,3)+2*(C123R(1,2)+C124R(3,2)
     -   )+3*(C124R(3,3)-C124R(4,3)))*P(8))
       FI(235) = 8*(M**2*(D1234I(1,1)+D1234I(1,2)-D1234I(3,1)+D1234I(3
     -   ,2)-2*D1234I(5,2))*P(2)*P(4)+(C123I(1,1)+C123I(1,3)-C124I(1,
     -   2)-C124I(1,3)-C124I(2,2)+C124I(2,3)+2*(C123I(1,2)+C124I(3,2)
     -   )+3*(C124I(3,3)-C124I(4,3)))*P(8))
       F(235)=DCMPLX(FR(235),FI(235))
       FR(236) = 4*((C123R(1,1)+C123R(1,2)+2*(C123R(3,2)+C123R(3,3)))*
     -   P(8)+M**2*P(2)*(D1234R(1,2)*P(1)*P(3)+(D1234R(1,1)-D1234R(3,
     -   1)+2*D1234R(3,2))*P(4)-D1234R(5,2)*P(7)+(D1234R(4,2)-D1234R(
     -   6,2))*P(60)*P(61)))
       FI(236) = 4*((C123I(1,1)+C123I(1,2)+2*(C123I(3,2)+C123I(3,3)))*
     -   P(8)+M**2*P(2)*(D1234I(1,2)*P(1)*P(3)+(D1234I(1,1)-D1234I(3,
     -   1)+2*D1234I(3,2))*P(4)-D1234I(5,2)*P(7)+(D1234I(4,2)-D1234I(
     -   6,2))*P(60)*P(61)))
       F(236)=DCMPLX(FR(236),FI(236))
       FR(237) = 4*((C123R(2,1)+3*C123R(3,2)+2*C123R(3,3))*P(8)+M**2*P
     -   (2)*(-((D1234R(3,1)-2*D1234R(3,2))*P(4))+D1234R(1,2)*P(9)*P(
     -   10)-D1234R(5,2)*P(15)-D1234R(1,1)*P(17)+D1234R(2,1)*P(18)*P(
     -   19)+(D1234R(4,2)-D1234R(6,2))*P(91)))
       FI(237) = 4*((C123I(2,1)+3*C123I(3,2)+2*C123I(3,3))*P(8)+M**2*P
     -   (2)*(-((D1234I(3,1)-2*D1234I(3,2))*P(4))+D1234I(1,2)*P(9)*P(
     -   10)-D1234I(5,2)*P(15)-D1234I(1,1)*P(17)+D1234I(2,1)*P(18)*P(
     -   19)+(D1234I(4,2)-D1234I(6,2))*P(91)))
       F(237)=DCMPLX(FR(237),FI(237))
       FR(238) = 4*((C123R(1,1)+C123R(1,2)+C123R(2,1)+3*C123R(3,2)+2*C
     -   123R(3,3))*P(8)+M**2*P(2)*(2*(-D1234R(3,1)+D1234R(3,2))*P(4)
     -   +(D1234R(1,1)+D1234R(1,2))*P(18)*P(19)-D1234R(5,2)*P(22)+(D1
     -   234R(2,1)+D1234R(4,2)-D1234R(6,2))*P(24)*P(25)))
       FI(238) = 4*((C123I(1,1)+C123I(1,2)+C123I(2,1)+3*C123I(3,2)+2*C
     -   123I(3,3))*P(8)+M**2*P(2)*(2*(-D1234I(3,1)+D1234I(3,2))*P(4)
     -   +(D1234I(1,1)+D1234I(1,2))*P(18)*P(19)-D1234I(5,2)*P(22)+(D1
     -   234I(2,1)+D1234I(4,2)-D1234I(6,2))*P(24)*P(25)))
       F(238)=DCMPLX(FR(238),FI(238))
       FR(239) = 4*((C123R(2,2)+C123R(3,2)+2*C123R(4,3))*P(8)+M**2*P(2
     -   )*(2*D1234R(3,2)*P(4)+D1234R(2,2)*P(18)*P(19)-D1234R(6,2)*P(
     -   22)+(D1234R(4,2)-D1234R(5,2))*P(24)*P(25)))
       FI(239) = 4*((C123I(2,2)+C123I(3,2)+2*C123I(4,3))*P(8)+M**2*P(2
     -   )*(2*D1234I(3,2)*P(4)+D1234I(2,2)*P(18)*P(19)-D1234I(6,2)*P(
     -   22)+(D1234I(4,2)-D1234I(5,2))*P(24)*P(25)))
       F(239)=DCMPLX(FR(239),FI(239))
       FR(240) = 4*((C123R(2,1)+C123R(2,2)+2*(C123R(3,2)+C123R(4,3)))*
     -   P(8)+M**2*P(2)*((D1234R(2,1)-D1234R(3,1)+2*D1234R(3,2))*P(4)
     -   +D1234R(2,2)*P(9)*P(10)-D1234R(6,2)*P(15)+(D1234R(4,2)-D1234
     -   R(5,2))*P(91)))
       FI(240) = 4*((C123I(2,1)+C123I(2,2)+2*(C123I(3,2)+C123I(4,3)))*
     -   P(8)+M**2*P(2)*((D1234I(2,1)-D1234I(3,1)+2*D1234I(3,2))*P(4)
     -   +D1234I(2,2)*P(9)*P(10)-D1234I(6,2)*P(15)+(D1234I(4,2)-D1234
     -   I(5,2))*P(91)))
       F(240)=DCMPLX(FR(240),FI(240))
       FR(241) = 4*((C123R(2,1)+C123R(3,2)+2*(C123R(2,2)+C123R(4,3)))*
     -   P(8)+M**2*P(2)*(D1234R(2,2)*P(1)*P(3)+(D1234R(2,1)-D1234R(3,
     -   1)+2*D1234R(3,2))*P(4)-D1234R(6,2)*P(7)+(D1234R(4,2)-D1234R(
     -   5,2))*P(60)*P(61)))
       FI(241) = 4*((C123I(2,1)+C123I(3,2)+2*(C123I(2,2)+C123I(4,3)))*
     -   P(8)+M**2*P(2)*(D1234I(2,2)*P(1)*P(3)+(D1234I(2,1)-D1234I(3,
     -   1)+2*D1234I(3,2))*P(4)-D1234I(6,2)*P(7)+(D1234I(4,2)-D1234I(
     -   5,2))*P(60)*P(61)))
       F(241)=DCMPLX(FR(241),FI(241))
       FR(242) = 8*(M**2*(D1234R(2,2)+D1234R(3,2)-2*D1234R(6,2))*P(2)*
     -   P(4)+(C123R(2,2)+C123R(2,3))*P(8))
       FI(242) = 8*(M**2*(D1234I(2,2)+D1234I(3,2)-2*D1234I(6,2))*P(2)*
     -   P(4)+(C123I(2,2)+C123I(2,3))*P(8))
       F(242)=DCMPLX(FR(242),FI(242))
       FR(243) = 8*(M**2*D1234R(3,2)*P(2)*P(4)+(C124R(2,2)+C124R(2,3))
     -   *P(8))
       FI(243) = 8*(M**2*D1234I(3,2)*P(2)*P(4)+(C124I(2,2)+C124I(2,3))
     -   *P(8))
       F(243)=DCMPLX(FR(243),FI(243))
       FR(244) = 2*(2*(-C234R(2,2)+C234R(3,2)-2*(C234R(2,3)-C234R(4,3)
     -   ))*P(8)+M**2*P(19)*(-((D1234R(2,1)+2*(D1234R(3,2)+D1234R(4,2
     -   )-D1234R(5,2)-D1234R(6,2)))*P(1)*P(49))+D1234R(1,1)*P(9)*P(6
     -   4)+D01234R*P(114)+D1234R(3,1)*P(118)))
       FI(244) = 2*(2*(-C234I(2,2)+C234I(3,2)-2*(C234I(2,3)-C234I(4,3)
     -   ))*P(8)+M**2*P(19)*(-((D1234I(2,1)+2*(D1234I(3,2)+D1234I(4,2
     -   )-D1234I(5,2)-D1234I(6,2)))*P(1)*P(49))+D1234I(1,1)*P(9)*P(6
     -   4)+D01234I*P(114)+D1234I(3,1)*P(118)))
       F(244)=DCMPLX(FR(244),FI(244))
       FR(245) = -4*((C234R(1,1)+C234R(1,2)-C234R(2,1)-C234R(2,2)+2*(C
     -   234R(2,3)+C234R(3,3)-2*C234R(4,3)))*P(8)+M**2*(D1234R(2,1)+D
     -   1234R(2,2)-D1234R(3,1)+D1234R(3,2)-2*D1234R(6,2))*P(9)*P(19)
     -   *P(64))
       FI(245) = -4*((C234I(1,1)+C234I(1,2)-C234I(2,1)-C234I(2,2)+2*(C
     -   234I(2,3)+C234I(3,3)-2*C234I(4,3)))*P(8)+M**2*(D1234I(2,1)+D
     -   1234I(2,2)-D1234I(3,1)+D1234I(3,2)-2*D1234I(6,2))*P(9)*P(19)
     -   *P(64))
       F(245)=DCMPLX(FR(245),FI(245))
       FR(246) = 4*((C234R(2,2)-C234R(3,2)-2*(C234R(2,3)+C234R(3,3)-2*
     -   C234R(4,3)))*P(8)+M**2*P(19)*(-((D1234R(2,2)+D1234R(3,2)-2*D
     -   1234R(6,2))*P(1)*P(49))+(-D1234R(2,1)+D1234R(3,1))*P(119)))
       FI(246) = 4*((C234I(2,2)-C234I(3,2)-2*(C234I(2,3)+C234I(3,3)-2*
     -   C234I(4,3)))*P(8)+M**2*P(19)*(-((D1234I(2,2)+D1234I(3,2)-2*D
     -   1234I(6,2))*P(1)*P(49))+(-D1234I(2,1)+D1234I(3,1))*P(119)))
       F(246)=DCMPLX(FR(246),FI(246))
       FR(247) = -2*(2*(C234R(2,2)-C234R(3,2)+2*(C234R(2,3)-C234R(4,3)
     -   ))*P(8)+M**2*P(19)*(2*(D1234R(3,2)+D1234R(4,2)-D1234R(5,2)-D
     -   1234R(6,2))*P(2)*P(18)+D1234R(1,1)*P(9)*P(64)+D01234R*P(114)
     -   +D1234R(2,1)*P(123)-D1234R(3,1)*P(126)))
       FI(247) = -2*(2*(C234I(2,2)-C234I(3,2)+2*(C234I(2,3)-C234I(4,3)
     -   ))*P(8)+M**2*P(19)*(2*(D1234I(3,2)+D1234I(4,2)-D1234I(5,2)-D
     -   1234I(6,2))*P(2)*P(18)+D1234I(1,1)*P(9)*P(64)+D01234I*P(114)
     -   +D1234I(2,1)*P(123)-D1234I(3,1)*P(126)))
       F(247)=DCMPLX(FR(247),FI(247))
       FR(248) = 4*((-C234R(1,1)-C234R(1,2)+C234R(2,1)+C234R(3,2)-2*(C
     -   234R(2,3)+C234R(3,3)-2*C234R(4,3)))*P(8)+M**2*P(19)*(-((D123
     -   4R(2,2)+D1234R(3,2)-2*D1234R(6,2))*P(2)*P(18))+(-D1234R(2,1)
     -   +D1234R(3,1))*P(114)))
       FI(248) = 4*((-C234I(1,1)-C234I(1,2)+C234I(2,1)+C234I(3,2)-2*(C
     -   234I(2,3)+C234I(3,3)-2*C234I(4,3)))*P(8)+M**2*P(19)*(-((D123
     -   4I(2,2)+D1234I(3,2)-2*D1234I(6,2))*P(2)*P(18))+(-D1234I(2,1)
     -   +D1234I(3,1))*P(114)))
       F(248)=DCMPLX(FR(248),FI(248))
       FR(249) = 4*((-C234R(2,2)+C234R(3,2)-2*(C234R(2,3)-C234R(4,3)))
     -   *P(8)+M**2*(-D1234R(3,2)+D1234R(6,2))*P(19)*P(86)*P(127))
       FI(249) = 4*((-C234I(2,2)+C234I(3,2)-2*(C234I(2,3)-C234I(4,3)))
     -   *P(8)+M**2*(-D1234I(3,2)+D1234I(6,2))*P(19)*P(86)*P(127))
       F(249)=DCMPLX(FR(249),FI(249))
       FR(250) = (B034R+8*(-C123R(6,3)+C234R(5,3)-C234R(6,3)))*P(8)+B0
     -   13R*P(41)+M**2*P(19)*(D01234R*P(128)+(D1234R(2,1)-D1234R(3,1
     -   ))*P(129))+C0123R*P(130)+C0234R*P(131)+C234R(1,1)*P(132)+C12
     -   3R(2,1)*P(134)+C234R(2,1)*P(135)
       FI(250) = (B034I+8*(-C123I(6,3)+C234I(5,3)-C234I(6,3)))*P(8)+B0
     -   13I*P(41)+M**2*P(19)*(D01234I*P(128)+(D1234I(2,1)-D1234I(3,1
     -   ))*P(129))+C0123I*P(130)+C0234I*P(131)+C234I(1,1)*P(132)+C12
     -   3I(2,1)*P(134)+C234I(2,1)*P(135)
       F(250)=DCMPLX(FR(250),FI(250))
       FR(251) = (B013R-4*(C123R(4,2)-C234R(4,2)+2*(C123R(6,3)-C234R(5
     -   ,3)+C234R(6,3))))*P(8)+B034R*P(41)-M**2*P(19)*(D01234R*P(128
     -   )+(D1234R(2,1)-D1234R(3,1))*P(136))+C0123R*P(137)+C0234R*P(1
     -   38)+C234R(1,1)*P(140)+C123R(2,1)*P(141)+C234R(2,1)*P(142)
       FI(251) = (B013I-4*(C123I(4,2)-C234I(4,2)+2*(C123I(6,3)-C234I(5
     -   ,3)+C234I(6,3))))*P(8)+B034I*P(41)-M**2*P(19)*(D01234I*P(128
     -   )+(D1234I(2,1)-D1234I(3,1))*P(136))+C0123I*P(137)+C0234I*P(1
     -   38)+C234I(1,1)*P(140)+C123I(2,1)*P(141)+C234I(2,1)*P(142)
       F(251)=DCMPLX(FR(251),FI(251))
       FR(252) = (B013R+B034R+2*(B13R(1,1)+B34R(1,1)-2*(C123R(4,2)-C23
     -   4R(4,2))+4*(-C123R(6,3)+C234R(5,3)-C234R(6,3))))*P(8)+C0123R
     -   *P(137)+C0234R*P(138)-M**2*P(19)*(D01234R*P(128)+(D1234R(2,1
     -   )-D1234R(3,1))*P(143))+C234R(1,1)*P(145)+C123R(2,1)*P(146)+C
     -   234R(2,1)*P(147)
       FI(252) = (B013I+B034I+2*(B13I(1,1)+B34I(1,1)-2*(C123I(4,2)-C23
     -   4I(4,2))+4*(-C123I(6,3)+C234I(5,3)-C234I(6,3))))*P(8)+C0123I
     -   *P(137)+C0234I*P(138)-M**2*P(19)*(D01234I*P(128)+(D1234I(2,1
     -   )-D1234I(3,1))*P(143))+C234I(1,1)*P(145)+C123I(2,1)*P(146)+C
     -   234I(2,1)*P(147)
       F(252)=DCMPLX(FR(252),FI(252))
       FR(253) = (B012R+B013R-8*(C123R(4,2)+C123R(5,3)+C234R(6,3)))*P(
     -   8)+(B024R+B034R)*P(41)+C123R(1,1)*P(141)+C234R(2,1)*P(142)-M
     -   **2*P(19)*((D1234R(1,1)-D1234R(3,1))*P(136)+D01234R*P(148))+
     -   C0123R*P(151)+C0234R*P(155)
       FI(253) = (B012I+B013I-8*(C123I(4,2)+C123I(5,3)+C234I(6,3)))*P(
     -   8)+(B024I+B034I)*P(41)+C123I(1,1)*P(141)+C234I(2,1)*P(142)-M
     -   **2*P(19)*((D1234I(1,1)-D1234I(3,1))*P(136)+D01234I*P(148))+
     -   C0123I*P(151)+C0234I*P(155)
       F(253)=DCMPLX(FR(253),FI(253))
       FR(254) = (B013R+B024R+B034R+2*(B13R(1,1)+B34R(1,1)-2*(C123R(4,
     -   2)+C234R(4,2)+2*(C123R(5,3)+C234R(6,3)))))*P(8)+B012R*P(41)+
     -   C123R(1,1)*P(146)+C234R(2,1)*P(147)-M**2*P(19)*((D1234R(1,1)
     -   -D1234R(3,1))*P(143)+D01234R*P(156))+C0123R*P(158)+C0234R*P(
     -   162)
       FI(254) = (B013I+B024I+B034I+2*(B13I(1,1)+B34I(1,1)-2*(C123I(4,
     -   2)+C234I(4,2)+2*(C123I(5,3)+C234I(6,3)))))*P(8)+B012I*P(41)+
     -   C123I(1,1)*P(146)+C234I(2,1)*P(147)-M**2*P(19)*((D1234I(1,1)
     -   -D1234I(3,1))*P(143)+D01234I*P(156))+C0123I*P(158)+C0234I*P(
     -   162)
       F(254)=DCMPLX(FR(254),FI(254))
       FR(255) = (B012R+B024R+B034R+2*(B12R(1,1)+B24R(1,1)-2*(C123R(4,
     -   2)+C234R(4,2)+2*(C123R(5,3)+C234R(6,3)))))*P(8)+B013R*P(41)+
     -   C123R(1,1)*P(134)+C234R(2,1)*P(135)+M**2*P(19)*((D1234R(1,1)
     -   -D1234R(3,1))*P(129)+D01234R*P(156))+C0234R*P(162)+C0123R*P(
     -   163)
       FI(255) = (B012I+B024I+B034I+2*(B12I(1,1)+B24I(1,1)-2*(C123I(4,
     -   2)+C234I(4,2)+2*(C123I(5,3)+C234I(6,3)))))*P(8)+B013I*P(41)+
     -   C123I(1,1)*P(134)+C234I(2,1)*P(135)+M**2*P(19)*((D1234I(1,1)
     -   -D1234I(3,1))*P(129)+D01234I*P(156))+C0234I*P(162)+C0123I*P(
     -   163)
       F(255)=DCMPLX(FR(255),FI(255))
       FR(256) = (B024R+B034R+2*B34R(1,1)-4*C234R(4,2)-8*C234R(6,3))*P
     -   (8)+C234R(2,1)*P(147)+C0234R*P(162)+M**2*P(19)*((C0134R-4*D1
     -   234R(7,2))*P(2)*P(18)+C0124R*P(9)*P(64)+D1234R(3,1)*P(143)+D
     -   01234R*P(164))
       FI(256) = (B024I+B034I+2*B34I(1,1)-4*C234I(4,2)-8*C234I(6,3))*P
     -   (8)+C234I(2,1)*P(147)+C0234I*P(162)+M**2*P(19)*((C0134I-4*D1
     -   234I(7,2))*P(2)*P(18)+C0124I*P(9)*P(64)+D1234I(3,1)*P(143)+D
     -   01234I*P(164))
       F(256)=DCMPLX(FR(256),FI(256))
       FR(257) = -8*C234R(6,3)*P(8)+(B024R+B034R)*P(41)+C234R(2,1)*P(1
     -   42)+C0234R*P(155)-M**2*P(19)*(C0134R*P(2)*P(18)+4*D1234R(7,2
     -   )*P(1)*P(49)+C0124R*P(9)*P(64)-D1234R(3,1)*P(136)+D01234R*P(
     -   164))
       FI(257) = -8*C234I(6,3)*P(8)+(B024I+B034I)*P(41)+C234I(2,1)*P(1
     -   42)+C0234I*P(155)-M**2*P(19)*(C0134I*P(2)*P(18)+4*D1234I(7,2
     -   )*P(1)*P(49)+C0124I*P(9)*P(64)-D1234I(3,1)*P(136)+D01234I*P(
     -   164))
       F(257)=DCMPLX(FR(257),FI(257))
       FR(258) = (B024R+B034R+2*B24R(1,1)-4*C234R(4,2)-8*C234R(6,3))*P
     -   (8)+C234R(2,1)*P(135)+C0234R*P(162)+M**2*P(19)*(C0134R*P(2)*
     -   P(18)+(C0124R-4*D1234R(7,2))*P(9)*P(64)-D1234R(3,1)*P(129)+D
     -   01234R*P(164))
       FI(258) = (B024I+B034I+2*B24I(1,1)-4*C234I(4,2)-8*C234I(6,3))*P
     -   (8)+C234I(2,1)*P(135)+C0234I*P(162)+M**2*P(19)*(C0134I*P(2)*
     -   P(18)+(C0124I-4*D1234I(7,2))*P(9)*P(64)-D1234I(3,1)*P(129)+D
     -   01234I*P(164))
       F(258)=DCMPLX(FR(258),FI(258))
       FR(259) = -2*(-2*(C234R(2,1)+C234R(3,2)-2*(C234R(2,3)-C234R(4,3
     -   )))*P(8)+M**2*P(19)*((D1234R(2,1)+2*(D1234R(3,2)+D1234R(4,2)
     -   -D1234R(5,2)-D1234R(6,2)))*P(1)*P(49)+D01234R*P(114)+D1234R(
     -   1,1)*P(86)*P(127)-D1234R(3,1)*P(167)))
       FI(259) = -2*(-2*(C234I(2,1)+C234I(3,2)-2*(C234I(2,3)-C234I(4,3
     -   )))*P(8)+M**2*P(19)*((D1234I(2,1)+2*(D1234I(3,2)+D1234I(4,2)
     -   -D1234I(5,2)-D1234I(6,2)))*P(1)*P(49)+D01234I*P(114)+D1234I(
     -   1,1)*P(86)*P(127)-D1234I(3,1)*P(167)))
       F(259)=DCMPLX(FR(259),FI(259))
       FR(260) = 2*(2*(C234R(2,1)-C234R(2,2)+2*(-C234R(2,3)+C234R(3,2)
     -   +C234R(4,3)))*P(8)+M**2*P(19)*(2*(-D1234R(3,2)-D1234R(4,2)+D
     -   1234R(5,2)+D1234R(6,2))*P(2)*P(18)-D01234R*P(114)+D1234R(3,1
     -   )*P(126)-D1234R(1,1)*P(86)*P(127)+D1234R(2,1)*P(171)))
       FI(260) = 2*(2*(C234I(2,1)-C234I(2,2)+2*(-C234I(2,3)+C234I(3,2)
     -   +C234I(4,3)))*P(8)+M**2*P(19)*(2*(-D1234I(3,2)-D1234I(4,2)+D
     -   1234I(5,2)+D1234I(6,2))*P(2)*P(18)-D01234I*P(114)+D1234I(3,1
     -   )*P(126)-D1234I(1,1)*P(86)*P(127)+D1234I(2,1)*P(171)))
       F(260)=DCMPLX(FR(260),FI(260))
       FR(261) = 2*(2*(C234R(2,1)-C234R(2,2)+2*(-C234R(2,3)+C234R(3,2)
     -   +C234R(4,3)))*P(8)+M**2*P(19)*(-((D1234R(1,1)+2*(D1234R(3,2)
     -   +D1234R(4,2)-D1234R(5,2)-D1234R(6,2)))*P(9)*P(64))-D01234R*P
     -   (114)+D1234R(2,1)*P(171)+D1234R(3,1)*P(174)))
       FI(261) = 2*(2*(C234I(2,1)-C234I(2,2)+2*(-C234I(2,3)+C234I(3,2)
     -   +C234I(4,3)))*P(8)+M**2*P(19)*(-((D1234I(1,1)+2*(D1234I(3,2)
     -   +D1234I(4,2)-D1234I(5,2)-D1234I(6,2)))*P(9)*P(64))-D01234I*P
     -   (114)+D1234I(2,1)*P(171)+D1234I(3,1)*P(174)))
       F(261)=DCMPLX(FR(261),FI(261))
       FR(262) = -2*(-2*(C234R(2,1)+C234R(3,2)-2*(C234R(2,3)-C234R(4,3
     -   )))*P(8)+M**2*P(19)*((D1234R(1,1)+2*(D1234R(3,2)+D1234R(4,2)
     -   -D1234R(5,2)-D1234R(6,2)))*P(9)*P(64)+D01234R*P(114)+D1234R(
     -   2,1)*P(179)-D1234R(3,1)*P(184)))
       FI(262) = -2*(-2*(C234I(2,1)+C234I(3,2)-2*(C234I(2,3)-C234I(4,3
     -   )))*P(8)+M**2*P(19)*((D1234I(1,1)+2*(D1234I(3,2)+D1234I(4,2)
     -   -D1234I(5,2)-D1234I(6,2)))*P(9)*P(64)+D01234I*P(114)+D1234I(
     -   2,1)*P(179)-D1234I(3,1)*P(184)))
       F(262)=DCMPLX(FR(262),FI(262))
       FR(263) = 4*((C234R(2,1)-C234R(2,2)+2*(-C234R(2,3)+C234R(3,2)+C
     -   234R(4,3)))*P(8)+M**2*P(19)*(D1234R(3,1)*P(114)+D1234R(3,2)*
     -   P(171)+D1234R(6,2)*P(187)))
       FI(263) = 4*((C234I(2,1)-C234I(2,2)+2*(-C234I(2,3)+C234I(3,2)+C
     -   234I(4,3)))*P(8)+M**2*P(19)*(D1234I(3,1)*P(114)+D1234I(3,2)*
     -   P(171)+D1234I(6,2)*P(187)))
       F(263)=DCMPLX(FR(263),FI(263))
       FR(264) = 4*((C234R(2,1)+C234R(3,2)-2*(C234R(2,3)-C234R(4,3)))*
     -   P(8)+M**2*P(19)*((-D1234R(3,2)+D1234R(6,2))*P(24)*P(33)+D123
     -   4R(3,1)*P(114)))
       FI(264) = 4*((C234I(2,1)+C234I(3,2)-2*(C234I(2,3)-C234I(4,3)))*
     -   P(8)+M**2*P(19)*((-D1234I(3,2)+D1234I(6,2))*P(24)*P(33)+D123
     -   4I(3,1)*P(114)))
       F(264)=DCMPLX(FR(264),FI(264))
       FR(265) = -4*(C123R(2,1)+3*C123R(3,2)+C234R(2,2)-C234R(3,2)+2*(
     -   C123R(3,3)+C234R(2,3)-C234R(4,3)))*P(8)
       FI(265) = -4*(C123I(2,1)+3*C123I(3,2)+C234I(2,2)-C234I(3,2)+2*(
     -   C123I(3,3)+C234I(2,3)-C234I(4,3)))*P(8)
       F(265)=DCMPLX(FR(265),FI(265))
       FR(266) = -4*(C123R(1,1)+C123R(1,2)+C123R(2,1)+3*C123R(3,2)-C23
     -   4R(2,1)-C234R(3,2)+2*(C123R(3,3)+C234R(2,3)-C234R(4,3)))*P(8
     -   )
       FI(266) = -4*(C123I(1,1)+C123I(1,2)+C123I(2,1)+3*C123I(3,2)-C23
     -   4I(2,1)-C234I(3,2)+2*(C123I(3,3)+C234I(2,3)-C234I(4,3)))*P(8
     -   )
       F(266)=DCMPLX(FR(266),FI(266))
       FR(267) = -4*(C123R(1,1)+C123R(1,2)-C234R(2,1)+C234R(2,2)+2*(C1
     -   23R(3,2)+C123R(3,3)+C234R(2,3)-C234R(3,2)-C234R(4,3)))*P(8)
       FI(267) = -4*(C123I(1,1)+C123I(1,2)-C234I(2,1)+C234I(2,2)+2*(C1
     -   23I(3,2)+C123I(3,3)+C234I(2,3)-C234I(3,2)-C234I(4,3)))*P(8)
       F(267)=DCMPLX(FR(267),FI(267))
       FR(268) = -4*(C123R(2,1)+C123R(3,2)-C234R(2,2)+C234R(3,2)+2*(C1
     -   23R(2,2)+C123R(4,3)+C234R(2,3)+C234R(3,3)-2*C234R(4,3)))*P(8
     -   )
       FI(268) = -4*(C123I(2,1)+C123I(3,2)-C234I(2,2)+C234I(3,2)+2*(C1
     -   23I(2,2)+C123I(4,3)+C234I(2,3)+C234I(3,3)-2*C234I(4,3)))*P(8
     -   )
       F(268)=DCMPLX(FR(268),FI(268))
       FR(269) = -4*(C123R(2,1)+C123R(2,2)+C234R(1,1)+C234R(1,2)-C234R
     -   (2,1)-C234R(2,2)+2*(C123R(3,2)+C123R(4,3)+C234R(2,3)+C234R(3
     -   ,3)-2*C234R(4,3)))*P(8)
       FI(269) = -4*(C123I(2,1)+C123I(2,2)+C234I(1,1)+C234I(1,2)-C234I
     -   (2,1)-C234I(2,2)+2*(C123I(3,2)+C123I(4,3)+C234I(2,3)+C234I(3
     -   ,3)-2*C234I(4,3)))*P(8)
       F(269)=DCMPLX(FR(269),FI(269))
       FR(270) = -4*(C123R(2,2)+C123R(3,2)+C234R(1,1)+C234R(1,2)-C234R
     -   (2,1)-C234R(3,2)+2*(C123R(4,3)+C234R(2,3)+C234R(3,3)-2*C234R
     -   (4,3)))*P(8)
       FI(270) = -4*(C123I(2,2)+C123I(3,2)+C234I(1,1)+C234I(1,2)-C234I
     -   (2,1)-C234I(3,2)+2*(C123I(4,3)+C234I(2,3)+C234I(3,3)-2*C234I
     -   (4,3)))*P(8)
       F(270)=DCMPLX(FR(270),FI(270))
       FR(271) = 8*(-C123R(2,2)-C123R(2,3)+C234R(1,2)+C234R(1,3)+C234R
     -   (2,2)-C234R(2,3)-2*C234R(3,2)-3*(C234R(3,3)-C234R(4,3)))*P(8
     -   )
       FI(271) = 8*(-C123I(2,2)-C123I(2,3)+C234I(1,2)+C234I(1,3)+C234I
     -   (2,2)-C234I(2,3)-2*C234I(3,2)-3*(C234I(3,3)-C234I(4,3)))*P(8
     -   )
       F(271)=DCMPLX(FR(271),FI(271))
       FR(272) = -8*(C234R(2,2)+C234R(2,3))*P(8)-4*M**2*(D1234R(1,1)+D
     -   1234R(1,2)-D1234R(3,1)+D1234R(3,2)-2*D1234R(5,2))*P(1)*P(19)
     -   *P(49)
       FI(272) = -8*(C234I(2,2)+C234I(2,3))*P(8)-4*M**2*(D1234I(1,1)+D
     -   1234I(1,2)-D1234I(3,1)+D1234I(3,2)-2*D1234I(5,2))*P(1)*P(19)
     -   *P(49)
       F(272)=DCMPLX(FR(272),FI(272))
       FR(273) = -4*(2*(C234R(2,2)+C234R(2,3))*P(8)+M**2*P(19)*((D1234
     -   R(1,2)+D1234R(3,2)-2*D1234R(5,2))*P(2)*P(18)+D01234R*P(114)+
     -   (D1234R(1,1)-D1234R(3,1))*P(189)))
       FI(273) = -4*(2*(C234I(2,2)+C234I(2,3))*P(8)+M**2*P(19)*((D1234
     -   I(1,2)+D1234I(3,2)-2*D1234I(5,2))*P(2)*P(18)+D01234I*P(114)+
     -   (D1234I(1,1)-D1234I(3,1))*P(189)))
       F(273)=DCMPLX(FR(273),FI(273))
       FR(274) = -4*(2*(C234R(2,2)+C234R(2,3))*P(8)+M**2*P(19)*((D1234
     -   R(1,2)+D1234R(3,2)-2*D1234R(5,2))*P(9)*P(64)+D01234R*P(114)+
     -   (D1234R(1,1)-D1234R(3,1))*P(190)))
       FI(274) = -4*(2*(C234I(2,2)+C234I(2,3))*P(8)+M**2*P(19)*((D1234
     -   I(1,2)+D1234I(3,2)-2*D1234I(5,2))*P(9)*P(64)+D01234I*P(114)+
     -   (D1234I(1,1)-D1234I(3,1))*P(190)))
       F(274)=DCMPLX(FR(274),FI(274))
       FR(275) = 4*(-2*(C234R(2,2)+C234R(2,3))*P(8)+M**2*P(19)*((-D123
     -   4R(3,2)+D1234R(5,2))*P(24)*P(33)+D1234R(3,1)*P(114)))
       FI(275) = 4*(-2*(C234I(2,2)+C234I(2,3))*P(8)+M**2*P(19)*((-D123
     -   4I(3,2)+D1234I(5,2))*P(24)*P(33)+D1234I(3,1)*P(114)))
       F(275)=DCMPLX(FR(275),FI(275))
       FR(276) = 4*(-2*(C234R(2,2)+C234R(2,3))*P(8)+M**2*P(19)*(D1234R
     -   (3,1)*P(114)+(-D1234R(3,2)+D1234R(5,2))*P(86)*P(127)))
       FI(276) = 4*(-2*(C234I(2,2)+C234I(2,3))*P(8)+M**2*P(19)*(D1234I
     -   (3,1)*P(114)+(-D1234I(3,2)+D1234I(5,2))*P(86)*P(127)))
       F(276)=DCMPLX(FR(276),FI(276))
       FR(277) = 4*(-2*(C234R(2,2)+C234R(2,3))*P(8)+M**2*P(19)*(2*D123
     -   4R(3,1)*P(114)+D1234R(3,2)*P(171)+D1234R(5,2)*P(187)))
       FI(277) = 4*(-2*(C234I(2,2)+C234I(2,3))*P(8)+M**2*P(19)*(2*D123
     -   4I(3,1)*P(114)+D1234I(3,2)*P(171)+D1234I(5,2)*P(187)))
       F(277)=DCMPLX(FR(277),FI(277))
       FR(278) = -8*((C234R(2,2)+C234R(2,3))*P(8)+M**2*D1234R(3,2)*P(1
     -   9)*P(114))
       FI(278) = -8*((C234I(2,2)+C234I(2,3))*P(8)+M**2*D1234I(3,2)*P(1
     -   9)*P(114))
       F(278)=DCMPLX(FR(278),FI(278))
       FR(279) = -8*(C123R(1,1)+2*C123R(1,2)+C123R(1,3)+C234R(2,2)+C23
     -   4R(2,3))*P(8)
       FI(279) = -8*(C123I(1,1)+2*C123I(1,2)+C123I(1,3)+C234I(2,2)+C23
     -   4I(2,3))*P(8)
       F(279)=DCMPLX(FR(279),FI(279))
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
c************************************************************************************
c************************************************************************************
c************************************************************************************
c       Amplitude                         
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************


       resultn(1) = mup1mup4*mup2mup3*F(1)+mup1mup3*mup2mup4*F(2)+mup1
     -   mup2*mup3mup4*F(3)+mup2mup4*p1mup1*p1mup3*F(4)+mup1mup4*p2mu
     -   p2*p2mup3*F(5)+mup2mup3*p1mup1*p1mup4*F(6)+mup2mup3*p1mup4*p
     -   2mup1*F(7)+mup1mup2*p1mup3*p1mup4*F(8)+p1mup1*p1mup2*p1mup3*
     -   p1mup4*F(9)+mup1mup2*p1mup4*p2mup3*F(10)+p1mup4*p4mup1*p4mup
     -   2*p4mup3*F(11)+mup1mup3*p2mup2*p2mup4*F(12)+mup1mup2*p2mup3*
     -   p2mup4*F(13)+p1mup1*p1mup3*p4mup2*(p1mup4*F(14)+p2mup4*F(15)
     -   +p4mup4*F(16))+p1mup1*p1mup2*p1mup3*(p2mup4*F(17)+p4mup4*F(1
     -   8))+p2mup1*p4mup2*p4mup3*(p1mup4*F(19)+p2mup4*F(20)+p4mup4*F
     -   (21))+p1mup1*p1mup3*p2mup2*(p1mup4*F(22)+p2mup4*F(23)+p4mup4
     -   *F(24))+p1mup2*p1mup3*p2mup1*(p1mup4*F(25)+p2mup4*F(26)+p4mu
     -   p4*F(27))+p1mup1*p1mup2*p2mup3*(p1mup4*F(28)+p2mup4*F(29)+p4
     -   mup4*F(30))+p2mup1*p2mup2*p4mup3*(p1mup4*F(31)+p2mup4*F(32)+
     -   p4mup4*F(33))+p2mup1*p2mup3*p4mup2*(p1mup4*F(34)+p2mup4*F(35
     -   )+p4mup4*F(36))+p2mup2*p2mup3*p4mup1*(p1mup4*F(37)+p2mup4*F(
     -   38)+p4mup4*F(39))+p1mup3*p2mup1*p2mup2*(p1mup4*F(40)+p2mup4*
     -   F(41)+p4mup4*F(42))+p1mup1*p2mup2*p2mup3*(p1mup4*F(43)+p2mup
     -   4*F(44)+p4mup4*F(45))+p1mup2*p2mup1*p2mup3*(p1mup4*F(23)+p2m
     -   up4*F(46)+p4mup4*F(47))+p2mup1*p2mup2*p2mup3*(p1mup4*F(48)+p
     -   2mup4*F(49)+p4mup4*F(50))+mup2mup4*p1mup1*p2mup3*F(51)+mup2m
     -   up3*p1mup1*p2mup4*F(52)+mup2mup4*p2mup1*p2mup3*F(53)+mup2mup
     -   3*p2mup1*p2mup4*F(54)+mup1mup3*p1mup4*p2mup2*F(55)+mup2mup4*
     -   p1mup3*p2mup1*F(56)+mup1mup3*p1mup2*p2mup4*F(57)+mup1mup3*p1
     -   mup2*p1mup4*F(58)+mup1mup2*p1mup3*p2mup4*F(59)+mup1mup4*p1mu
     -   p2*p2mup3*F(60)+mup1mup4*p1mup3*p2mup2*F(61)+mup3mup4*p2mup2
     -   *p4mup1*F(62)+mup2mup4*p2mup3*p4mup1*F(63)+mup1mup2*p2mup3*p
     -   4mup4*F(64)+mup1mup3*p2mup2*p4mup4*F(65)+mup2mup4*p2mup1*p4m
     -   up3*F(66)+mup1mup3*p2mup4*p4mup2*F(67)+mup2mup4*p1mup1*p4mup
     -   3*F(68)+mup2mup4*p1mup3*p4mup1*F(69)+mup1mup3*p1mup2*p4mup4*
     -   F(70)+mup1mup3*p1mup4*p4mup2*F(71)+mup1mup2*p1mup3*p4mup4*F(
     -   72)+mup3mup4*p1mup2*p4mup1*F(73)+mup1mup2*p2mup4*p4mup3*F(74
     -   )+mup3mup4*p2mup1*p4mup2*F(75)+mup1mup4*p2mup2*p4mup3*F(76)+
     -   mup1mup4*p2mup3*p4mup2*F(77)+mup2mup3*p2mup1*p4mup4*F(78)+mu
     -   p2mup3*p2mup4*p4mup1*F(79)+mup2mup3*p1mup1*p4mup4*F(80)+mup1
     -   mup2*p1mup4*p4mup3*F(81)+mup2mup3*p1mup4*p4mup1*F(82)+mup3mu
     -   p4*p1mup1*p4mup2*F(83)+mup3mup4*p1mup1*p1mup2*F(84)+mup3mup4
     -   *p1mup1*p2mup2*F(85)+mup3mup4*p1mup2*p2mup1*F(86)+mup3mup4*p
     -   2mup1*p2mup2*F(87)+mup1mup4*p1mup2*p1mup3*F(88)+mup1mup4*p1m
     -   up2*p4mup3*F(89)+mup1mup4*p1mup3*p4mup2*F(90)+mup1mup2*p4mup
     -   3*p4mup4*F(91)+mup3mup4*p4mup1*p4mup2*F(92)+mup2mup3*p4mup1*
     -   p4mup4*F(93)+mup1mup4*p4mup2*p4mup3*F(94)+mup2mup4*p4mup1*p4
     -   mup3*F(95)+mup1mup3*p4mup2*p4mup4*F(96)+p1mup3*p2mup1*p4mup2
     -   *(p1mup4*F(97)+p2mup4*F(98)+p4mup4*F(99))+p1mup1*p2mup2*p4mu
     -   p3*(p1mup4*F(100)+p2mup4*F(101)+p4mup4*F(102))+p1mup2*p2mup1
     -   *p4mup3*(p1mup4*F(103)+p2mup4*F(104)+p4mup4*F(105))+p1mup3*p
     -   2mup2*p4mup1*(p1mup4*F(106)+p2mup4*F(107)+p4mup4*F(108))+p1m
     -   up1*p2mup3*p4mup2*(p1mup4*F(109)+p2mup4*F(110)+p4mup4*F(111)
     -   )+p1mup2*p2mup3*p4mup1*(p1mup4*F(112)+p2mup4*F(113)+p4mup4*F
     -   (114))+p1mup1*p1mup2*p4mup3*(p1mup4*F(115)+p2mup4*F(116)+p4m
     -   up4*F(117))+p1mup2*p1mup3*p4mup1*(p1mup4*F(118)+p2mup4*F(119
     -   )+p4mup4*F(120))+p1mup2*p4mup1*p4mup3*(p1mup4*F(16)+p2mup4*F
     -   (121)+p4mup4*F(122))+p1mup1*p4mup2*p4mup3*(p1mup4*F(123)+p2m
     -   up4*F(124)+p4mup4*F(125))+p1mup3*p4mup1*p4mup2*(p1mup4*F(126
     -   )+p2mup4*F(127)+p4mup4*F(128))+p2mup2*p4mup1*p4mup3*(p2mup4*
     -   F(36)+p1mup4*F(129)+p4mup4*F(130))+p2mup3*p4mup1*p4mup2*(p1m
     -   up4*F(131)+p2mup4*F(132)+p4mup4*F(133))+p4mup1*p4mup2*p4mup3
     -   *(p2mup4*F(134)+p4mup4*F(135))

       if(printmy)then
       print*, "RESULT",resultn(1)
       endif


        if (comp.gt.0) then    

       if(printmy)then
       resultn(2) = abs(mup2mup4*p2mup3*F(136))
     &    +abs(p1mup4*p2mup2*p2mup3*F(137))
     -   +abs(mup2mup3*p2mup4*F(138))
     &   +abs(p1mup2*p1mup3*p2mup4*F(139))
     &   +abs(p1mup3*p2mup2*p2mup4*F(140))
     &   +abs(p1mup2*p2mup3*p2mup4*F(141))
     &   +abs(p2mup2*p2mup3*p2mup4*F(142))
     &   +abs(p2mup3*p2mup4*p4mup2*F(143))
     &   +abs(p2mup2*p2mup4*p4mup3*F(144))
     &   +abs(p2mup4*p4mup2*p4mup3*F(145))
     &   +abs(p2mup2*p2mup3*p4mup4*F(146))
     &   +abs(mup3mup4*p2mup2*F(147))
     &   +abs(mup2mup3*p1mup4*F(148))
     &   +abs(mup2mup4*p1mup3*F(149))
     &   +abs(mup3mup4*p1mup2*F(150))
     &   +abs(p1mup3*p2mup4*p4mup2*F(151))
     &   +abs(p1mup2*p2mup4*p4mup3*F(152))
     &   +abs(mup2mup3*p4mup4*F(153))
     &   +abs(mup2mup4*p4mup3*F(154))
     &   +abs(mup3mup4*p4mup2*F(155))
     &   +abs(p1mup2*p1mup4*p2mup3*F(156))
     &   +abs(p1mup3*p1mup4*p2mup2*F(157))
     &   +abs(p1mup3*p2mup2*p4mup4*F(158))
     &   +abs(p1mup2*p2mup3*p4mup4*F(159))
     &   +abs(p1mup4*p2mup3*p4mup2*F(160))
     &   +abs(p1mup4*p2mup2*p4mup3*F(161))
     &   +abs(p2mup2*p4mup3*p4mup4*F(162))
     &   +abs(p2mup3*p4mup2*p4mup4*F(163))
     &   +abs(p1mup3*p1mup4*p4mup2*F(164))
     &   +abs(p1mup2*p1mup3*p1mup4*F(165))
     &   +abs(p1mup2*p1mup3*p4mup4*F(166))
     &   +abs(p1mup2*p1mup4*p4mup3*F(167))
     &   +abs(p1mup3*p4mup2*p4mup4*F(168))
     &   +abs(p1mup4*p4mup2*p4mup3*F(169))
     &   +abs(p1mup2*p4mup3*p4mup4*F(170))
     &   +abs(p4mup2*p4mup3*p4mup4*F(171))

       print*, "ABS resultn(2)",resultn(2)
       endif


       resultn(2) = mup2mup4*p2mup3*F(136)+p1mup4*p2mup2*p2mup3*F(137)
     -   +mup2mup3*p2mup4*F(138)+p1mup2*p1mup3*p2mup4*F(139)+p1mup3*p
     -   2mup2*p2mup4*F(140)+p1mup2*p2mup3*p2mup4*F(141)+p2mup2*p2mup
     -   3*p2mup4*F(142)+p2mup3*p2mup4*p4mup2*F(143)+p2mup2*p2mup4*p4
     -   mup3*F(144)+p2mup4*p4mup2*p4mup3*F(145)+p2mup2*p2mup3*p4mup4
     -   *F(146)+mup3mup4*p2mup2*F(147)+mup2mup3*p1mup4*F(148)+mup2mu
     -   p4*p1mup3*F(149)+mup3mup4*p1mup2*F(150)+p1mup3*p2mup4*p4mup2
     -   *F(151)+p1mup2*p2mup4*p4mup3*F(152)+mup2mup3*p4mup4*F(153)+m
     -   up2mup4*p4mup3*F(154)+mup3mup4*p4mup2*F(155)+p1mup2*p1mup4*p
     -   2mup3*F(156)+p1mup3*p1mup4*p2mup2*F(157)+p1mup3*p2mup2*p4mup
     -   4*F(158)+p1mup2*p2mup3*p4mup4*F(159)+p1mup4*p2mup3*p4mup2*F(
     -   160)+p1mup4*p2mup2*p4mup3*F(161)+p2mup2*p4mup3*p4mup4*F(162)
     -   +p2mup3*p4mup2*p4mup4*F(163)+p1mup3*p1mup4*p4mup2*F(164)+p1m
     -   up2*p1mup3*p1mup4*F(165)+p1mup2*p1mup3*p4mup4*F(166)+p1mup2*
     -   p1mup4*p4mup3*F(167)+p1mup3*p4mup2*p4mup4*F(168)+p1mup4*p4mu
     -   p2*p4mup3*F(169)+p1mup2*p4mup3*p4mup4*F(170)+p4mup2*p4mup3*p
     -   4mup4*F(171)


       if(printmy) then
       print*, "ABS resultn(2)",resultn(2)
       endif

       resultn(3) = mup1mup3*p1mup4*F(172)+p1mup1*p1mup3*p1mup4*F(173)
     -   +p1mup3*p1mup4*p2mup1*F(174)+p1mup3*p1mup4*p4mup1*F(175)+p1m
     -   up1*p1mup4*p2mup3*F(176)+p1mup4*p2mup1*p2mup3*F(177)+p1mup1*
     -   p1mup4*p4mup3*F(178)+p1mup4*p4mup1*p4mup3*F(179)+mup1mup3*p2
     -   mup4*F(180)+p1mup1*p1mup3*p2mup4*F(181)+p1mup3*p2mup1*p2mup4
     -   *F(182)+p1mup1*p2mup3*p2mup4*F(183)+p2mup1*p2mup3*p2mup4*F(1
     -   84)+p2mup3*p2mup4*p4mup1*F(185)+p2mup1*p2mup4*p4mup3*F(186)+
     -   p1mup1*p1mup3*p4mup4*F(187)+p2mup1*p2mup3*p4mup4*F(188)+mup1
     -   mup4*p2mup3*F(189)+mup3mup4*p2mup1*F(190)+mup3mup4*p1mup1*F(
     -   191)+mup1mup4*p1mup3*F(192)+p1mup4*p2mup3*p4mup1*F(193)+p1mu
     -   p4*p2mup1*p4mup3*F(194)+p1mup3*p2mup1*p4mup4*F(195)+p1mup1*p
     -   2mup3*p4mup4*F(196)+p1mup3*p2mup4*p4mup1*F(197)+p1mup1*p2mup
     -   4*p4mup3*F(198)+mup3mup4*p4mup1*F(199)+mup1mup3*p4mup4*F(200
     -   )+mup1mup4*p4mup3*F(201)+p1mup3*p4mup1*p4mup4*F(202)+p1mup1*
     -   p4mup3*p4mup4*F(203)+p2mup4*p4mup1*p4mup3*F(204)+p2mup3*p4mu
     -   p1*p4mup4*F(205)+p2mup1*p4mup3*p4mup4*F(206)+p4mup1*p4mup3*p
     -   4mup4*F(207)
       resultn(4) = p1mup1*p1mup4*p4mup2*F(208)+p2mup2*p2mup4*p4mup1*F
     -   (209)+p2mup1*p2mup4*p4mup2*F(210)+p2mup1*p2mup2*p4mup4*F(211
     -   )+mup1mup2*p1mup4*F(212)+mup2mup4*p1mup1*F(213)+mup1mup2*p2m
     -   up4*F(214)+mup1mup4*p2mup2*F(215)+mup2mup4*p2mup1*F(216)+p1m
     -   up2*p1mup4*p4mup1*F(217)+p1mup1*p1mup2*p4mup4*F(218)+p1mup4*
     -   p2mup2*p4mup1*F(219)+p1mup2*p2mup4*p4mup1*F(220)+p1mup1*p2mu
     -   p4*p4mup2*F(221)+p1mup1*p2mup2*p4mup4*F(222)+p1mup2*p2mup1*p
     -   4mup4*F(223)+p1mup4*p2mup1*p4mup2*F(224)+mup2mup4*p4mup1*F(2
     -   25)+mup1mup2*p4mup4*F(226)+mup1mup4*p1mup2*F(227)+mup1mup4*p
     -   4mup2*F(228)+p2mup1*p4mup2*p4mup4*F(229)+p1mup1*p4mup2*p4mup
     -   4*F(230)+p2mup4*p4mup1*p4mup2*F(231)+p1mup4*p4mup1*p4mup2*F(
     -   232)+p2mup2*p4mup1*p4mup4*F(233)+p1mup2*p4mup1*p4mup4*F(234)
     -   +p1mup1*p1mup2*p1mup4*F(235)+p1mup1*p1mup4*p2mup2*F(236)+p1m
     -   up2*p1mup4*p2mup1*F(237)+p1mup1*p1mup2*p2mup4*F(238)+p1mup4*
     -   p2mup1*p2mup2*F(239)+p1mup1*p2mup2*p2mup4*F(240)+p1mup2*p2mu
     -   p1*p2mup4*F(241)+p2mup1*p2mup2*p2mup4*F(242)+p4mup1*p4mup2*p
     -   4mup4*F(243)
       resultn(5) = p1mup3*p2mup1*p4mup2*F(244)+p2mup2*p2mup3*p4mup1*F
     -   (245)+p2mup1*p2mup3*p4mup2*F(246)+p1mup2*p2mup1*p4mup3*F(247
     -   )+p2mup1*p2mup2*p4mup3*F(248)+p2mup1*p4mup2*p4mup3*F(249)+mu
     -   p2mup3*p2mup1*F(250)+mup1mup3*p2mup2*F(251)+mup1mup2*p2mup3*
     -   F(252)+mup1mup3*p1mup2*F(253)+mup1mup2*p1mup3*F(254)+mup2mup
     -   3*p1mup1*F(255)+mup1mup2*p4mup3*F(256)+mup1mup3*p4mup2*F(257
     -   )+mup2mup3*p4mup1*F(258)+p1mup1*p2mup3*p4mup2*F(259)+p1mup1*
     -   p2mup2*p4mup3*F(260)+p1mup3*p2mup2*p4mup1*F(261)+p1mup2*p2mu
     -   p3*p4mup1*F(262)+p2mup2*p4mup1*p4mup3*F(263)+p2mup3*p4mup1*p
     -   4mup2*F(264)+p1mup2*p1mup3*p2mup1*F(265)+p1mup1*p1mup2*p2mup
     -   3*F(266)+p1mup1*p1mup3*p2mup2*F(267)+p1mup2*p2mup1*p2mup3*F(
     -   268)+p1mup1*p2mup2*p2mup3*F(269)+p1mup3*p2mup1*p2mup2*F(270)
     -   +p2mup1*p2mup2*p2mup3*F(271)+p1mup1*p1mup3*p4mup2*F(272)+p1m
     -   up1*p1mup2*p4mup3*F(273)+p1mup2*p1mup3*p4mup1*F(274)+p1mup3*
     -   p4mup1*p4mup2*F(275)+p1mup1*p4mup2*p4mup3*F(276)+p1mup2*p4mu
     -   p1*p4mup3*F(277)+p4mup1*p4mup2*p4mup3*F(278)+p1mup1*p1mup2*p
     -   1mup3*F(279)
       resultn(2) =-resultn(2)
       resultn(3) =-resultn(3)
       resultn(4) =-resultn(4)
       resultn(5) =-resultn(5)


        endif    


c Following the Peskin-S convention this contribution should get a factor i from the Feyman rules
c (including the factor -1 for closed fermion loop)
c and a factor i coming from the scalar integrals.
c -> i*i = -1, this facor is changed here
       result=-resultn(1)
c************************************************************************************
c************************************************************************************
c************************************************************************************
       Return
       End
