c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
       subroutine XEgetPenOddG_R5_Div(M5,M1,M2,M3,M4,p1sq,p2sq,p3sq,p4sq,p5sq,p1p2
     & ,p1p3,p1p4,p2p3,p2p4,p3p4,p1p5,p4p5,musq,Div)
       Implicit none
       Real*8   p1sq, p1p2, p1p3, p1p4, p1p5 
       Real*8   p2sq, p2p3, p2p4
       Real*8   p3sq, p3p4
       Real*8   p4sq, p4p5 
       Real*8   p5sq 
       Real*8   s12, s13, s14, s15 
       Real*8   s23, s24
       Real*8   s34
       Real*8   s45 
       Real*8   s234 
       Real*8 M5,M1,M2,M3,M4,musq
       Complex*16 A01, A02,A03,A04,A05
       Real*8      A01R, A02R,A03R,A04R,A05R
       Real*8      A01I, A02I,A03I,A04I,A05I
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
       Real*8  B12R(3,3),B13R(3,3),B14R(3,3),B15R(3,3) 
       Real*8  B23R(3,3),B24R(3,3),B25R(3,3) 
       Real*8  B34R(3,3),B35R(3,3) 
       Real*8  B45R(3,3) 
       Real*8  B12I(3,3),B13I(3,3),B14I(3,3),B15I(3,3) 
       Real*8  B23I(3,3),B24I(3,3),B25I(3,3) 
       Real*8  B34I(3,3),B35I(3,3) 
       Real*8  B45I(3,3) 

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
       Real*8 D1234R(34,5),D1235R(34,5) 
       Real*8 D1245R(34,5) 
       Real*8 D1345R(34,5) 
       Real*8 D2345R(34,5) 
       Real*8 D1234I(34,5),D1235I(34,5) 
       Real*8 D1245I(34,5) 
       Real*8 D1345I(34,5) 
       Real*8 D2345I(34,5) 
       Complex*16 E012345 
       Real*8 E012345R 
       Real*8 E012345I 
       Real*8 E12345R(80,5) 
       Real*8 E12345I(80,5) 

       Logical PrintB,PrintC,PrintD,PrintE 
       parameter(printB=.false., printC=.false.,printD=.false.,printE=.false.)
       COMMON/EVALFunctionsPenOddG/ E012345R,E12345R,E012345I,E12345I
       COMMON/DVALFunctionsPenOddG/ D01234R,D1234R,D01234I,D1234I, D01235R,D1
     -   235R,D01235I,D1235I, D01245R,D1245R,D01245I,D1245I, D01345R,
     -   D1345R,D01345I,D1345I, D02345R,D2345R,D02345I,D2345I
       COMMON/CVALFunctionsPenOddG/ C0123R,C123R,C0123I,C123I, C0124R,C124R,C
     -   0124I,C124I, C0125R,C125R,C0125I,C125I, C0134R,C134R,C0134I,
     -   C134I, C0135R,C135R,C0135I,C135I, C0145R,C145R,C0145I,C145I,
     -    C0234R,C234R,C0234I,C234I, C0235R,C235R,C0235I,C235I, C0245
     -   R,C245R,C0245I,C245I, C0345R,C345R,C0345I,C345I 
       COMMON/BVALFunctionsPenOddG/ B012R,B012I, B013R,B013I, B014R,B014I, B0
     -   15R,B015I, B023R,B023I, B024R,B024I, B025R,B025I, B034R,B034
     -   I, B035R,B035I, B045R,B045I,B12R,B12I, B13R,B13I, B14R,
     -   B14I, B15R,B15I, B23R,B23I, B24R,B24I, B25R,B25I, B34R
     -   ,B34I, B35R,B35I, B45R,B45I 

       complex*16 A0finGDiv,B0finGDiv,C0finGDiv,D0finGDiv,E0finG
       external A0finGDiv,B0finGDiv,C0finGDiv,D0finGDiv,E0finG
       integer Div

       s12 = (p1sq +p2sq+ 2*p1p2) 
       s13 = (p1sq +p3sq+ 2*p1p3) 
       s14 = (p1sq +p4sq+ 2*p1p4) 
       s23 = (p2sq +p3sq+ 2*p2p3) 
       s24 = (p2sq +p4sq+ 2*p2p4) 
       s34 = (p3sq +p4sq+ 2*p3p4) 

       s15 = (p1sq +p5sq+ 2*p1p5) 
       s45 = (p4sq +p5sq+ 2*p4p5) 

       s234=s15

        A01=A0finGDiv(M5,musq,Div)  
        A01R=Dble(A01)   
        A01I=DImag(A01)   
c
        A02=A0finGDiv(M1,musq,Div)  
        A02R=Dble(A02)   
        A02I=DImag(A02)   
c
        A03=A0finGDiv(M2,musq,Div)  
        A03R=Dble(A03)   
        A03I=DImag(A03)   
c
        A04=A0finGDiv(M3,musq,Div)  
        A04R=Dble(A04)   
        A04I=DImag(A04)  
c
        A05=A0finGDiv(M4,musq,Div)  
        A05R=Dble(A05)   
        A05I=DImag(A05)  
c
c
       B012=B0finGDiv(M5,M1,p1sq,musq,Div)
       B012R=Dble(B012)
       B012I=DImag(B012)
       call tens_red2_R3_G_Div(M5,M1,p1sq,A02,A02R,A02I,A01,A01R,A01I,B012
     -   ,B012R,B012I,B12R,B12I) 
c
       B023=B0finGDiv(M1,M2,p2sq,musq,Div)
       B023R=Dble(B023)
       B023I=DImag(B023)
       call tens_red2_R3_G_Div(M1,M2,p2sq,A03,A03R,A03I,A02,A02R,A02I,B023
     -   ,B023R,B023I,B23R,B23I) 

       B034=B0finGDiv(M2,M3,p3sq,musq,Div)
       B034R=Dble(B034)
       B034I=DImag(B034)
       call tens_red2_R3_G_Div(M2,M3,p3sq,A04,A04R,A04I,A03,A03R,A03I,B034
     -   ,B034R,B034I,B34R,B34I) 

       B045=B0finGDiv(M3,M4,p4sq,musq,Div)
       B045R=Dble(B045)
       B045I=DImag(B045)
       call tens_red2_R3_G_Div(M3,M4,p4sq,A05,A05R,A05I,A04,A04R,A04I,B045
     -   ,B045R,B045I,B45R,B45I) 

       B013=B0finGDiv(M5,M2,s12,musq,Div)
       B013R=Dble(B013)
       B013I=DImag(B013)
       call tens_red2_R3_G_Div(M5,M2,s12,A03,A03R,A03I,A01,A01R,A01I,B013,
     -   B013R,B013I,B13R,B13I) 

       B014=B0finGDiv(M5,M3,s45,musq,Div)
       B014R=Dble(B014)
       B014I=DImag(B014)
       call tens_red2_R3_G_Div(M5,M3,s45,A04,A04R,A04I,A01,A01R,A01I,B014
     -   ,B014R,B014I,B14R,B14I) 

       B024=B0finGDiv(M1,M3,s23,musq,Div)
       B024R=Dble(B024)
       B024I=DImag(B024)
       call tens_red2_R3_G_Div(M1,M3,s23,A04,A04R,A04I,A02,A02R,A02I,B024,
     -   B024R,B024I,B24R,B24I) 

       B025=B0finGDiv(M1,M4,s234,musq,Div)
       B025R=Dble(B025)
       B025I=DImag(B025)
       call tens_red2_R3_G_Div(M1,M4,s234,A05,A05R,A05I,A02,A02R,A02I,B025,
     -   B025R,B025I,B25R,B25I) 

       B035=B0finGDiv(M2,M4,s34,musq,Div)
       B035R=Dble(B035)
       B035I=DImag(B035)
       call tens_red2_R3_G_Div(M2,M4,s34,A05,A05R,A05I,A03,A03R,A03I,B035,
     -   B035R,B035I,B35R,B35I) 

       B015=B0finGDiv(M5,M4,p5sq,musq,Div)
       B015R=Dble(B015)
       B015I=DImag(B015)
       call tens_red2_R3_G_Div(M5,M4,p5sq,A05,A05R,A05I,A01,A01R,A01I,B015,
     -   B015R,B015I,B15R,B15I) 
c
c
c
       C0123=C0finGDiv(M5,M1,M2,p1sq,p2sq,s12,musq,Div)
       call tens_red3_new_Re_Com_R4_G_Div(M5,M1,M2,p1sq,p2sq,s12,
     &  B023R,B013R,B012R,
     &  B23R,B13R,B12R,
     &  B023I,B013I,B012I,
     &  B23I,B13I,B12I,
     &  C0123,C0123R,C0123I,C123R,C123I)

       C0124=C0finGDiv(M5,M1,M3,p1sq,s23,s45,musq,Div)
       call tens_red3_new_Re_Com_R4_G_Div(M5,M1,M3,p1sq,s23,s45,
     & B024R,B014R,B012R,
     & B24R,B14R,B12R,
     & B024I,B014I,B012I,
     & B24I,B14I,B12I,
     & C0124,C0124R,C0124I,C124R,C124I)

       C0125=C0finGDiv(M5,M1,M4,p1sq,s234,p5sq,musq,Div)
       call tens_red3_new_Re_Com_R4_G_Div(M5,M1,M4,p1sq,s234,p5sq,
     & B025R,B015R,B012R,
     & B25R,B15R,B12R,
     & B025I,B015I,B012I,
     & B25I,B15I,B12I,
     & C0125,C0125R,C0125I,C125R,C125I)

       C0234=C0finGDiv(M1,M2,M3,p2sq,p3sq,s23,musq,Div)
       call tens_red3_new_Re_Com_R4_G_Div(M1,M2,M3,p2sq,p3sq,s23,
     & B034R,B024R,B023R,
     & B34R,B24R,B23R,
     & B034I,B024I,B023I,
     & B34I,B24I,B23I,
     & C0234,C0234R,C0234I,C234R,C234I)

       C0235=C0finGDiv(M1,M2,M4,p2sq,s34,s234,musq,Div)
       call tens_red3_new_Re_Com_R4_G_Div(M1,M2,M4,p2sq,s34,s234,
     & B035R,B025R,B023R,
     & B35R,B25R,B23R,
     & B035I,B025I,B023I,
     & B35I,B25I,B23I,
     & C0235,C0235R,C0235I,C235R,C235I)

       C0345=C0finGDiv(M2,M3,M4,p3sq,p4sq,s34,musq,Div)
       call tens_red3_new_Re_Com_R4_G_Div(M2,M3,M4,p3sq,p4sq,s34,
     & B045R,B035R,B034R,
     & B45R,B35R,B34R,
     & B045I,B035I,B034I,
     & B45I,B35I,B34I,
     & C0345,C0345R,C0345I,C345R,C345I)

       C0134=C0finGDiv(M5,M2,M3,s12,p3sq,s45,musq,Div)
       call tens_red3_new_Re_Com_R4_G_Div(M5,M2,M3,s12,p3sq,s45,
     & B034R,B014R,B013R,
     & B34R,B14R,B13R,
     & B034I,B014I,B013I,
     & B34I,B14I,B13I,
     & C0134,C0134R,C0134I,C134R,C134I)

       C0135=C0finGDiv(M5,M2,M4,s12,s34,p5sq,musq,Div)
       call tens_red3_new_Re_Com_R4_G_Div(M5,M2,M4,s12,s34,p5sq,
     & B035R,B015R,B013R,
     & B35R,B15R,B13R,
     & B035I,B015I,B013I,
     & B35I,B15I,B13I,
     & C0135,C0135R,C0135I,C135R,C135I)

       C0145=C0finGDiv(M5,M3,M4,s45,p4sq,p5sq,musq,Div)
       call tens_red3_new_Re_Com_R4_G_Div(M5,M3,M4,s45,p4sq,p5sq,
     & B045R,B015R,B014R,
     & B45R,B15R,B14R,
     & B045I,B015I,B014I,
     & B45I,B15I,B14I,
     & C0145,C0145R,C0145I,C145R,C145I)

       C0245=C0finGDiv(M1,M3,M4,s23,p4sq,s234,musq,Div)
       call tens_red3_new_Re_Com_R4_G_Div(M1,M3,M4,s23,p4sq,s234,
     & B045R,B025R,B024R,
     & B45R,B25R,B24R,
     & B045I,B025I,B024I,
     & B45I,B25I,B24I,
     & C0245,C0245R,C0245I,C245R,C245I)

c
c
c
       D01234=D0finGDiv(M5,M1,M2,M3,s12,s23,p1sq,p2sq,p3sq,s45,musq,Div)
       call tens_red4_new_Re_Com_R5_G_Div(M5,M1,M2,M3,p1sq,p2sq,p3sq,p1p2
     - ,p1p3,p2p3,C02
     -   34R,C0134R,C0124R,C0123R,C234R,C134R,C124R,C123R,C0234I,C013
     -   4I,C0124I,C0123I,C234I,C134I,C124I,C123I,D01234,D01234R,D012
     -   34I,D1234R,D1234I)

       D01235=D0finGDiv(M5,M1,M2,M4,s12,s234,p1sq,p2sq,s34,p5sq,musq,Div)
       call tens_red4_new_Re_Com_R5_G_Div(M5,M1,M2,M4,p1sq,p2sq,s34,p1p2
     - ,p1p3+p1p4,p2p3
     -   +p2p4,C0235R,C0135R,C0125R,C0123R,C235R,C135R,C125R,C123R,C0
     -   235I,C0135I,C0125I,C0123I,C235I,C135I,C125I,C123I,D01235,D01
     -   235R,D01235I,D1235R,D1235I)

       D01245=D0finGDiv(M5,M1,M3,M4,s45,s234,p1sq,s23,p4sq,p5sq,musq,Div)
       call tens_red4_new_Re_Com_R5_G_Div(M5,M1,M3,M4,p1sq,s23,p4sq,
     -   p1p2+p1p3,p1p4,p2p4
     -   +p3p4,C0245R,C0145R,C0125R,C0124R,C245R,C145R,C125R,C124R,C0
     -   245I,C0145I,C0125I,C0124I,C245I,C145I,C125I,C124I,D01245,D01
     -   245R,D01245I,D1245R,D1245I)

       D02345=D0finGDiv(M1,M2,M3,M4,s23,s34,p2sq,p3sq,p4sq,s234,musq,Div)
       call tens_red4_new_Re_Com_R5_G_Div(M1,M2,M3,M4,p2sq,p3sq,p4sq,p2p3
     -  ,p2p4,p3p4,C03
     -   45R,C0245R,C0235R,C0234R,C345R,C245R,C235R,C234R,C0345I,C024
     -   5I,C0235I,C0234I,C345I,C245I,C235I,C234I,D02345,D02345R,D023
     -   45I,D2345R,D2345I)

       D01345=D0finGDiv(M5,M2,M3,M4,s45,s34,s12,p3sq,p4sq,p5sq,musq,Div)
       call tens_red4_new_Re_Com_R5_G_Div(M5,M2,M3,M4,s12,p3sq,p4sq,
     - p1p3+p2p3,p1p4+p2p4
     -   ,p3p4,C0345R,C0145R,C0135R,C0134R,C345R,C145R,C135R,C134R,C0
     -   345I,C0145I,C0135I,C0134I,C345I,C145I,C135I,C134I,D01345,D01
     -   345R,D01345I,D1345R,D1345I)

       E012345=E0finG(M5,M1,M2,M3,M4,p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23
     -  ,s34,s45,s234,D02345,D01345,D01245,D01235,D01234)

       E012345R=Dble(E012345)
       E012345I=DImag(E012345)
       call tens_red5_new_Re_Com_R5_G_Div(M5,M1,M2,M3,M4,p1sq,p2sq,p3sq,
     -   p4sq,p1p2,p1p3,p1p4,p2
     -   p3,p2p4,p3p4,D02345R,D01345R,D01245R,D01235R,D01234R,D2345R,
     -   D1345R,D1245R,D1235R,D1234R,D02345I,D01345I,D01245I,D01235I,
     -   D01234I,D2345I,D1345I,D1245I,D1235I,D1234I,E12345R,E12345I)

c FC%
c FC%       If(PrintB) then
c FC%       write(*,*)"B12"
c FC%       write(*,*) "B12=", B012
c FC%       endif
c FC%       If(PrintB) then
c FC%       write(*,*)"B23"
c FC%       write(*,*) "B23=", B023
c FC%       endif
c FC%       If(PrintB) then
c FC%       write(*,*)"B34"
c FC%       write(*,*) "B34=", B034
c FC%       endif
c FC%       If(PrintB) then
c FC%       write(*,*)"B45"
c FC%       write(*,*) "B45=", B045
c FC%       endif
c FC%       If(PrintB) then
c FC%       write(*,*)"B13"
c FC%       write(*,*) "B13=", B013
c FC%       endif
c FC%       If(PrintB) then
c FC%       write(*,*)"B14"
c FC%       write(*,*) "B14=", B014
c FC%       endif
c FC%       If(PrintB) then
c FC%       write(*,*)"B24"
c FC%       write(*,*) "B24=", B024
c FC%       endif
c FC%       If(PrintB) then
c FC%       write(*,*)"B25"
c FC%       write(*,*) "B25=", B025
c FC%       endif
c FC%       If(PrintB) then
c FC%       write(*,*)"B35"
c FC%       write(*,*) "B35=", B035
c FC%       endif
c FC%       If(PrintB) then
c FC%       write(*,*)"B15"
c FC%       write(*,*) "B15=", B015
c FC%       endif
c FC%       If(PrintC) then
c FC%       write(*,*)"C123"
c FC%       call print_CijR_R4(C0123,C123R,C123I)
c FC%       Endif
c FC%       If(PrintC) then
c FC%       write(*,*)"C124"
c FC%       call print_CijR_R4(C0124,C124R,C124I)
c FC%       Endif
c FC%       If(PrintC) then
c FC%       write(*,*)"C125"
c FC%       call print_CijR_R4(C0125,C125R,C125I)
c FC%       Endif
c FC%       If(PrintC) then
c FC%       write(*,*)"C234"
c FC%       call print_CijR_R4(C0234,C234R,C234I)
c FC%       Endif
c FC%       If(PrintC) then
c FC%       write(*,*)"C235"
c FC%       call print_CijR_R4(C0235,C235R,C235I)
c FC%       Endif
c FC%       If(PrintC) then
c FC%       write(*,*)"C345"
c FC%       call print_CijR_R4(C0345,C345R,C345I)
c FC%       Endif
c FC%       If(PrintC) then
c FC%       write(*,*)"C134"
c FC%       call print_CijR_R4(C0134,C134R,C134I)
c FC%       Endif
c FC%       If(PrintC) then
c FC%       write(*,*)"C135"
c FC%       call print_CijR_R4(C0135,C135R,C135I)
c FC%       Endif
c FC%       If(PrintC) then
c FC%       write(*,*)"C145"
c FC%       call print_CijR_R4(C0145,C145R,C145I)
c FC%       Endif
c FC%       If(PrintC) then
c FC%       write(*,*)"C245"
c FC%       call print_CijR_R4(C0245,C245R,C245I)
c FC%       Endif
c FC%       If(PrintD) then
c FC%       write(*,*)"D1234"
c FC%       call print_DijR_R5(D01234,D1234R,D1234I)
c FC%       endif
c FC%       If(PrintD) then
c FC%       write(*,*)"D1235"
c FC%       call print_DijR_R5(D01235,D1235R,D1235I)
c FC%       endif
c FC%       If(PrintD) then
c FC%       write(*,*)"D1245"
c FC%       call print_DijR_R5(D01245,D1245R,D1245I)
c FC%       endif
c FC%       If(PrintD) then
c FC%       write(*,*)"D2345"
c FC%       call print_DijR_R5(D02345,D2345R,D2345I)
c FC%       endif
c FC%       If(PrintD) then
c FC%       write(*,*)"D1345"
c FC%       call print_DijR_R5(D01345,D1345R,D1345I)
c FC%       endif
c FC%       If(PrintE) then
c FC%       write(*,*)"E12345"
c FC%       call print_EijR(E012345,E12345R,E12345I)
c FC%       endif
c FC%


       return
       end
