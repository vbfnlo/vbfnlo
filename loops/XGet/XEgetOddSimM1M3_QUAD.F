       subroutine XEgetPenOddGSim_R5_QUAD(M5,M1,M2,M3,M4,p1sq,p2sq,p3sq,p4sq,p5sq,p1p2
     & ,p1p3,p1p4,p2p3,p2p4,p3p4,p1p5,p4p5,musq)
       Implicit none
       Real*16   p1sq, p1p2, p1p3, p1p4, p1p5 
       Real*16   p2sq, p2p3, p2p4
       Real*16   p3sq, p3p4
       Real*16   p4sq, p4p5 
       Real*16   p5sq 
       Real*16   s12, s13, s14, s15 
       Real*16   s23, s24
       Real*16   s34
       Real*16   s45 
       Real*16   s234 
       Real*16 M5,M1,M2,M3,M4,musq
       complex*32 A01, A02,A03,A04,A05
       Real*16      A01R, A02R,A03R,A04R,A05R
       Real*16      A01I, A02I,A03I,A04I,A05I
       complex*32  B012,B013,B014,B015 
       complex*32  B023,B024,B025 
       complex*32  B034,B035 
       complex*32  B045 
       Real*16  B012R,B013R,B014R,B015R 
       Real*16  B023R,B024R,B025R 
       Real*16  B034R,B035R 
       Real*16  B045R 
       Real*16  B012I,B013I,B014I,B015I 
       Real*16  B023I,B024I,B025I 
       Real*16  B034I,B035I 
       Real*16  B045I 
       Real*16  B12R(3,3),B13R(3,3),B14R(3,3),B15R(3,3) 
       Real*16  B23R(3,3),B24R(3,3),B25R(3,3) 
       Real*16  B34R(3,3),B35R(3,3) 
       Real*16  B45R(3,3) 
       Real*16  B12I(3,3),B13I(3,3),B14I(3,3),B15I(3,3) 
       Real*16  B23I(3,3),B24I(3,3),B25I(3,3) 
       Real*16  B34I(3,3),B35I(3,3) 
       Real*16  B45I(3,3) 

       complex*32 C0123,C0124,C0125 
       complex*32 C0134,C0135 
       complex*32 C0145 
       complex*32 C0234,C0235 
       complex*32 C0245 
       complex*32 C0345   
       Real*16 C0123R,C0124R,C0125R 
       Real*16 C0134R,C0135R 
       Real*16 C0145R 
       Real*16 C0234R,C0235R 
       Real*16 C0245R 
       Real*16 C0345R 
       Real*16 C0123I,C0124I,C0125I 
       Real*16 C0134I,C0135I 
       Real*16 C0145I 
       Real*16 C0234I,C0235I 
       Real*16 C0245I 
       Real*16 C0345I  
       Real*16 C123R(9,4),C124R(9,4),C125R(9,4) 
       Real*16 C134R(9,4),C135R(9,4) 
       Real*16 C145R(9,4) 
       Real*16 C234R(9,4),C235R(9,4) 
       Real*16 C245R(9,4) 
       Real*16 C345R(9,4) 
       Real*16 C123I(9,4),C124I(9,4),C125I(9,4) 
       Real*16 C134I(9,4),C135I(9,4) 
       Real*16 C145I(9,4) 
       Real*16 C234I(9,4),C235I(9,4) 
       Real*16 C245I(9,4) 
       Real*16 C345I(9,4) 
       complex*32  D01234,D01235 
       complex*32 D01245 
       complex*32 D01345 
       complex*32 D02345    
       Real*16 D01234R,D01235R 
       Real*16 D01245R 
       Real*16 D01345R 
       Real*16 D02345R  
       Real*16 D01234I,D01235I 
       Real*16 D01245I 
       Real*16 D01345I 
       Real*16 D02345I 
       Real*16 D1234R(34,5),D1235R(34,5) 
       Real*16 D1245R(34,5) 
       Real*16 D1345R(34,5) 
       Real*16 D2345R(34,5) 
       Real*16 D1234I(34,5),D1235I(34,5) 
       Real*16 D1245I(34,5) 
       Real*16 D1345I(34,5) 
       Real*16 D2345I(34,5) 
       complex*32 E012345 
       Real*16 E012345R 
       Real*16 E012345I 
       Real*16 E12345R(80,5) 
       Real*16 E12345I(80,5) 


cfc       COMMON/EVALFunctionsPenoddgSim/ E012345R,E12345R,E012345I,E12345I
cfc       COMMON/DVALFunctionsPenoddgSim/ D01234R,D1234R,D01234I,D1234I, D01235R,D1
cfc     -   235R,D01235I,D1235I, D01245R,D1245R,D01245I,D1245I, D01345R,
cfc     -   D1345R,D01345I,D1345I, D02345R,D2345R,D02345I,D2345I
cfc       COMMON/CVALFunctionsPenoddgSim/ C0123R,C123R,C0123I,C123I, C0124R,C124R,C
cfc     -   0124I,C124I, C0125R,C125R,C0125I,C125I, C0134R,C134R,C0134I,
cfc     -   C134I, C0135R,C135R,C0135I,C135I, C0145R,C145R,C0145I,C145I,
cfc     -    C0234R,C234R,C0234I,C234I, C0235R,C235R,C0235I,C235I, C0245
cfc     -   R,C245R,C0245I,C245I, C0345R,C345R,C0345I,C345I 
cfc       COMMON/BVALFunctionsPenoddgSim/ B012R,B012I, B013R,B013I, B014R,B014I, B0
cfc     -   15R,B015I, B023R,B023I, B024R,B024I, B025R,B025I, B034R,B034
cfc     -   I, B035R,B035I, B045R,B045I,B12R,B12I, B13R,B13I, B14R,
cfc     -   B14I, B15R,B15I, B23R,B23I, B24R,B24I, B25R,B25I, B34R
cfc     -   ,B34I, B35R,B35I, B45R,B45I 

       COMMON/EVALFunctionsPenoddgsim/ E012345R_d,E12345R_d,
     &                               E012345I_d,E12345I_d
       COMMON/DVALFunctionsPenoddgsim/ D01234R_d,D1234R_d,D01234I_d,
     &                D1234I_d, D01235R_d,D1235R_d,
     -                D01235I_d,D1235I_d, D01245R_d,D1245R_d,
     &                D01245I_d,D1245I_d, D01345R_d,
     -                D1345R_d,D01345I_d,D1345I_d, D02345R_d,
     &                D2345R_d,D02345I_d,D2345I_d
       COMMON/CVALFunctionsPenoddgsim/ C0123R_d,C123R_d,C0123I_d,C123I_d, 
     &             C0124R_d,C124R_d,C0124I_d,
     -             C124I_d, C0125R_d,C125R_d,C0125I_d,C125I_d,
     &             C0134R_d,C134R_d,C0134I_d,
     -             C134I_d, C0135R_d,C135R_d,C0135I_d,C135I_d,
     &             C0145R_d,C145R_d,C0145I_d,C145I_d,
     -             C0234R_d,C234R_d,C0234I_d,C234I_d, C0235R_d,C235R_d,
     &             C0235I_d,C235I_d, C0245R_d,
     -             C245R_d,C0245I_d,C245I_d, C0345R_d,C345R_d,C0345I_d,
     &             C345I_d 
       COMMON/BVALFunctionsPenoddgsim/ B012R_d,B012I_d, B013R_d,B013I_d, 
     &           B014R_d,B014I_d, B015R_d,
     -           B015I_d, B023R_d,B023I_d, B024R_d,B024I_d, 
     &           B025R_d,B025I_d, B034R_d,B034I_d,
     -           B035R_d,B035I_d, B045R_d,B045I_d,B12R_d,B12I_d, B13R_d,
     &           B13I_d, B14R_d,
     -           B14I_d, B15R_d,B15I_d, B23R_d,B23I_d, B24R_d,B24I_d, 
     &           B25R_d,B25I_d, B34R_d,
     -           B34I_d, B35R_d,B35I_d, B45R_d,B45I_d 


       Real*8  B012R_d,B013R_d,B014R_d,B015R_d 
       Real*8  B023R_d,B024R_d,B025R_d 
       Real*8  B034R_d,B035R_d 
       Real*8  B045R_d 
       Real*8  B012I_d,B013I_d,B014I_d,B015I_d 
       Real*8  B023I_d,B024I_d,B025I_d 
       Real*8  B034I_d,B035I_d 
       Real*8  B045I_d 
       Real*8  B12R_d(3,3),B13R_d(3,3),B14R_d(3,3),B15R_d(3,3) 
       Real*8  B23R_d(3,3),B24R_d(3,3),B25R_d(3,3) 
       Real*8  B34R_d(3,3),B35R_d(3,3) 
       Real*8  B45R_d(3,3) 
       Real*8  B12I_d(3,3),B13I_d(3,3),B14I_d(3,3),B15I_d(3,3) 
       Real*8  B23I_d(3,3),B24I_d(3,3),B25I_d(3,3) 
       Real*8  B34I_d(3,3),B35I_d(3,3) 
       Real*8  B45I_d(3,3) 


       Real*8 C0123R_d,C0124R_d,C0125R_d 
       Real*8 C0134R_d,C0135R_d 
       Real*8 C0145R_d 
       Real*8 C0234R_d,C0235R_d 
       Real*8 C0245R_d 
       Real*8 C0345R_d 
       Real*8 C0123I_d,C0124I_d,C0125I_d 
       Real*8 C0134I_d,C0135I_d 
       Real*8 C0145I_d 
       Real*8 C0234I_d,C0235I_d 
       Real*8 C0245I_d 
       Real*8 C0345I_d  
       Real*8 C123R_d(9,4),C124R_d(9,4),C125R_d(9,4) 
       Real*8 C134R_d(9,4),C135R_d(9,4) 
       Real*8 C145R_d(9,4) 
       Real*8 C234R_d(9,4),C235R_d(9,4) 
       Real*8 C245R_d(9,4) 
       Real*8 C345R_d(9,4) 
       Real*8 C123I_d(9,4),C124I_d(9,4),C125I_d(9,4) 
       Real*8 C134I_d(9,4),C135I_d(9,4) 
       Real*8 C145I_d(9,4) 
       Real*8 C234I_d(9,4),C235I_d(9,4) 
       Real*8 C245I_d(9,4) 
       Real*8 C345I_d(9,4) 

       Real*8 D01234R_d,D01235R_d 
       Real*8 D01245R_d 
       Real*8 D01345R_d 
       Real*8 D02345R_d  
       Real*8 D01234I_d,D01235I_d 
       Real*8 D01245I_d 
       Real*8 D01345I_d 
       Real*8 D02345I_d 
       Real*8 D1234R_d(34,5),D1235R_d(34,5) 
       Real*8 D1245R_d(34,5) 
       Real*8 D1345R_d(34,5) 
       Real*8 D2345R_d(34,5) 
       Real*8 D1234I_d(34,5),D1235I_d(34,5) 
       Real*8 D1245I_d(34,5) 
       Real*8 D1345I_d(34,5) 
       Real*8 D2345I_d(34,5) 

       Real*8 E012345R_d 
       Real*8 E012345I_d 
       Real*8 E12345R_d(80,5) 
       Real*8 E12345I_d(80,5) 

       complex*32 A0finG_QUAD,B0finG_QUAD,C0finG_QUAD,D0finG_QUAD
       complex*32 E0finG_QUAD
       external A0finG_QUAD,B0finG_QUAD,C0finG_QUAD,D0finG_QUAD
       external E0finG_QUAD

c       complex*32 A0finG,B0finG,C0finG,D0finG,E0finG
c       external A0finG,B0finG,C0finG,D0finG,E0finG

       s12 = (p1sq +p2sq+ 2*p1p2) 
       s13 = (p1sq +p3sq+ 2*p1p3) 
       s14 = (p1sq +p4sq+ 2*p1p4) 
       s23 = (p2sq +p3sq+ 2*p2p3) 
       s24 = (p2sq +p4sq+ 2*p2p4) 
       s34 = (p3sq +p4sq+ 2*p3p4) 

       s15 = (p1sq +p5sq+ 2*p1p5) 
       s45 = (p4sq +p5sq+ 2*p4p5) 

       s234=s15


        A01=A0finG_QUAD(M5,musq)  
        A01R=REAL(A01)   
        A01I=IMAG(A01)   
c
        A02=A0finG_QUAD(M1,musq)  
        A02R=REAL(A02)   
        A02I=IMAG(A02)   
c
        A03=A0finG_QUAD(M2,musq)  
        A03R=REAL(A03)   
        A03I=IMAG(A03)   
c
        A04=A0finG_QUAD(M3,musq)  
        A04R=REAL(A04)   
        A04I=IMAG(A04)  
c
        A05=A0finG_QUAD(M4,musq)  
        A05R=REAL(A05)   
        A05I=IMAG(A05)  
c
c
       B012=B0finG_QUAD(M5,M1,p1sq,musq)
       B012R=REAL(B012)
       B012I=IMAG(B012)
       call tens_red2_R3_G_QUAD(M5,M1,p1sq,A02,A02R,A02I,A01,A01R,A01I,B012
     -   ,B012R,B012I,B12R,B12I) 
c
       B023=B0finG_QUAD(M1,M2,p2sq,musq)
       B023R=REAL(B023)
       B023I=IMAG(B023)
       call tens_red2_R3_G_QUAD(M1,M2,p2sq,A03,A03R,A03I,A02,A02R,A02I,B023
     -   ,B023R,B023I,B23R,B23I) 

       B034=B0finG_QUAD(M2,M3,p3sq,musq)
       B034R=REAL(B034)
       B034I=IMAG(B034)
       call tens_red2_R3_G_QUAD(M2,M3,p3sq,A04,A04R,A04I,A03,A03R,A03I,B034
     -   ,B034R,B034I,B34R,B34I) 

       B045=B0finG_QUAD(M3,M4,p4sq,musq)
       B045R=REAL(B045)
       B045I=IMAG(B045)
       call tens_red2_R3_G_QUAD(M3,M4,p4sq,A05,A05R,A05I,A04,A04R,A04I,B045
     -   ,B045R,B045I,B45R,B45I) 

       B013=B0finG_QUAD(M5,M2,s12,musq)
       B013R=REAL(B013)
       B013I=IMAG(B013)
       call tens_red2_R3_G_QUAD(M5,M2,s12,A03,A03R,A03I,A01,A01R,A01I,B013,
     -   B013R,B013I,B13R,B13I) 

       B014=B0finG_QUAD(M5,M3,s45,musq)
       B014R=REAL(B014)
       B014I=IMAG(B014)
       call tens_red2_R3_G_QUAD(M5,M3,s45,A04,A04R,A04I,A01,A01R,A01I,B014
     -   ,B014R,B014I,B14R,B14I) 

       B024=B0finG_QUAD(M1,M3,s23,musq)
       B024R=REAL(B024)
       B024I=IMAG(B024)
       call tens_red2_R3_G_QUAD(M1,M3,s23,A04,A04R,A04I,A02,A02R,A02I,B024,
     -   B024R,B024I,B24R,B24I) 

       B025=B0finG_QUAD(M1,M4,s234,musq)
       B025R=REAL(B025)
       B025I=IMAG(B025)
       call tens_red2_R3_G_QUAD(M1,M4,s234,A05,A05R,A05I,A02,A02R,A02I,B025,
     -   B025R,B025I,B25R,B25I) 

       B035=B0finG_QUAD(M2,M4,s34,musq)
       B035R=REAL(B035)
       B035I=IMAG(B035)
       call tens_red2_R3_G_QUAD(M2,M4,s34,A05,A05R,A05I,A03,A03R,A03I,B035,
     -   B035R,B035I,B35R,B35I) 

       B015=B0finG_QUAD(M5,M4,p5sq,musq)
       B015R=REAL(B015)
       B015I=IMAG(B015)
       call tens_red2_R3_G_QUAD(M5,M4,p5sq,A05,A05R,A05I,A01,A01R,A01I,B015,
     -   B015R,B015I,B15R,B15I) 
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

       C0124=C0finG_QUAD(M5,M1,M3,p1sq,s23,s45,musq)
       call tens_red3_new_Re_Com_R4_G_QUAD(M5,M1,M3,p1sq,s23,s45,
     & B024R,B014R,B012R,
     & B24R,B14R,B12R,
     & B024I,B014I,B012I,
     & B24I,B14I,B12I,
     & C0124,C0124R,C0124I,C124R,C124I)

       C0125=C0finG_QUAD(M5,M1,M4,p1sq,s234,p5sq,musq)
       call tens_red3_new_Re_Com_R4_G_QUAD(M5,M1,M4,p1sq,s234,p5sq,
     & B025R,B015R,B012R,
     & B25R,B15R,B12R,
     & B025I,B015I,B012I,
     & B25I,B15I,B12I,
     & C0125,C0125R,C0125I,C125R,C125I)

       C0234=C0finG_QUAD(M1,M2,M3,p2sq,p3sq,s23,musq)
       call tens_red3_new_Re_Com_R4_G_QUAD(M1,M2,M3,p2sq,p3sq,s23,
     & B034R,B024R,B023R,
     & B34R,B24R,B23R,
     & B034I,B024I,B023I,
     & B34I,B24I,B23I,
     & C0234,C0234R,C0234I,C234R,C234I)

       C0235=C0finG_QUAD(M1,M2,M4,p2sq,s34,s234,musq)
       call tens_red3_new_Re_Com_R4_G_QUAD(M1,M2,M4,p2sq,s34,s234,
     & B035R,B025R,B023R,
     & B35R,B25R,B23R,
     & B035I,B025I,B023I,
     & B35I,B25I,B23I,
     & C0235,C0235R,C0235I,C235R,C235I)

       C0345=C0finG_QUAD(M2,M3,M4,p3sq,p4sq,s34,musq)
       call tens_red3_new_Re_Com_R4_G_QUAD(M2,M3,M4,p3sq,p4sq,s34,
     & B045R,B035R,B034R,
     & B45R,B35R,B34R,
     & B045I,B035I,B034I,
     & B45I,B35I,B34I,
     & C0345,C0345R,C0345I,C345R,C345I)

       C0134=C0finG_QUAD(M5,M2,M3,s12,p3sq,s45,musq)
       call tens_red3_new_Re_Com_R4_G_QUAD(M5,M2,M3,s12,p3sq,s45,
     & B034R,B014R,B013R,
     & B34R,B14R,B13R,
     & B034I,B014I,B013I,
     & B34I,B14I,B13I,
     & C0134,C0134R,C0134I,C134R,C134I)

       C0135=C0finG_QUAD(M5,M2,M4,s12,s34,p5sq,musq)
       call tens_red3_new_Re_Com_R4_G_QUAD(M5,M2,M4,s12,s34,p5sq,
     & B035R,B015R,B013R,
     & B35R,B15R,B13R,
     & B035I,B015I,B013I,
     & B35I,B15I,B13I,
     & C0135,C0135R,C0135I,C135R,C135I)

       C0145=C0finG_QUAD(M5,M3,M4,s45,p4sq,p5sq,musq)
       call tens_red3_new_Re_Com_R4_G_QUAD(M5,M3,M4,s45,p4sq,p5sq,
     & B045R,B015R,B014R,
     & B45R,B15R,B14R,
     & B045I,B015I,B014I,
     & B45I,B15I,B14I,
     & C0145,C0145R,C0145I,C145R,C145I)

       C0245=C0finG_QUAD(M1,M3,M4,s23,p4sq,s234,musq)
       call tens_red3_new_Re_Com_R4_G_QUAD(M1,M3,M4,s23,p4sq,s234,
     & B045R,B025R,B024R,
     & B45R,B25R,B24R,
     & B045I,B025I,B024I,
     & B45I,B25I,B24I,
     & C0245,C0245R,C0245I,C245R,C245I)

c
c       print*, "here_XGETEM1M3after C"
c
       D01234=D0finG_QUAD(M5,M1,M2,M3,s12,s23,p1sq,p2sq,p3sq,s45,musq)
       call tens_red4_new_Re_Com_R5_G_QUAD(M5,M1,M2,M3,p1sq,p2sq,p3sq,p1p2
     - ,p1p3,p2p3,C02
     -   34R,C0134R,C0124R,C0123R,C234R,C134R,C124R,C123R,C0234I,C013
     -   4I,C0124I,C0123I,C234I,C134I,C124I,C123I,D01234,D01234R,D012
     -   34I,D1234R,D1234I)

       D01235=D0finG_QUAD(M5,M1,M2,M4,s12,s234,p1sq,p2sq,s34,p5sq,musq)
       call tens_red4_new_Re_Com_R5_G_QUAD(M5,M1,M2,M4,p1sq,p2sq,s34,p1p2
     - ,p1p3+p1p4,p2p3
     -   +p2p4,C0235R,C0135R,C0125R,C0123R,C235R,C135R,C125R,C123R,C0
     -   235I,C0135I,C0125I,C0123I,C235I,C135I,C125I,C123I,D01235,D01
     -   235R,D01235I,D1235R,D1235I)

       D01245=D0finG_QUAD(M5,M1,M3,M4,s45,s234,p1sq,s23,p4sq,p5sq,musq)
       call tens_red4_new_Re_Com_R5_G_QUAD(M5,M1,M3,M4,p1sq,s23,p4sq,
     -   p1p2+p1p3,p1p4,p2p4
     -   +p3p4,C0245R,C0145R,C0125R,C0124R,C245R,C145R,C125R,C124R,C0
     -   245I,C0145I,C0125I,C0124I,C245I,C145I,C125I,C124I,D01245,D01
     -   245R,D01245I,D1245R,D1245I)

       D02345=D0finG_QUAD(M1,M2,M3,M4,s23,s34,p2sq,p3sq,p4sq,s234,musq)
       call tens_red4_new_Re_Com_R5_G_QUAD(M1,M2,M3,M4,p2sq,p3sq,p4sq,p2p3
     -  ,p2p4,p3p4,C03
     -   45R,C0245R,C0235R,C0234R,C345R,C245R,C235R,C234R,C0345I,C024
     -   5I,C0235I,C0234I,C345I,C245I,C235I,C234I,D02345,D02345R,D023
     -   45I,D2345R,D2345I)

       D01345=D0finG_QUAD(M5,M2,M3,M4,s45,s34,s12,p3sq,p4sq,p5sq,musq)
       call tens_red4_new_Re_Com_R5_G_QUAD(M5,M2,M3,M4,s12,p3sq,p4sq,
     - p1p3+p2p3,p1p4+p2p4
     -   ,p3p4,C0345R,C0145R,C0135R,C0134R,C345R,C145R,C135R,C134R,C0
     -   345I,C0145I,C0135I,C0134I,C345I,C145I,C135I,C134I,D01345,D01
     -   345R,D01345I,D1345R,D1345I)


c
c       print*, "here_XGETEM1M3after D"
c


       E012345=E0finG_QUAD(M5,M1,M2,M3,M4,p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23
     -  ,s34,s45,s234,D02345,D01345,D01245,D01235,D01234)

       E012345R=REAL(E012345)
       E012345I=IMAG(E012345)
       call tens_red5_new_Re_Com_R5_G_QUAD(M5,M1,M2,M3,M4,p1sq,p2sq,p3sq,
     -   p4sq,p1p2,p1p3,p1p4,p2
     -   p3,p2p4,p3p4,D02345R,D01345R,D01245R,D01235R,D01234R,D2345R,
     -   D1345R,D1245R,D1235R,D1234R,D02345I,D01345I,D01245I,D01235I,
     -   D01234I,D2345I,D1345I,D1245I,D1235I,D1234I,E12345R,E12345I)

c       print*, "here_XGETEM1M3after E"
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                    B0
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       call R3_From_BijR_QUAD_to_Dbl(B012R_d,B012I_d,B12R_d,B12I_d,
     &    B012R,B012I,B12R,B12I)
       call R3_From_BijR_QUAD_to_Dbl(B013R_d,B013I_d,B13R_d,B13I_d,
     &    B013R,B013I,B13R,B13I)
       call R3_From_BijR_QUAD_to_Dbl(B014R_d,B014I_d,B14R_d,B14I_d,
     &    B014R,B014I,B14R,B14I)
       call R3_From_BijR_QUAD_to_Dbl(B015R_d,B015I_d,B15R_d,B15I_d,
     &    B015R,B015I,B15R,B15I)
       call R3_From_BijR_QUAD_to_Dbl(B023R_d,B023I_d,B23R_d,B23I_d,
     &    B023R,B023I,B23R,B23I)
       call R3_From_BijR_QUAD_to_Dbl(B024R_d,B024I_d,B24R_d,B24I_d,
     &    B024R,B024I,B24R,B24I)
       call R3_From_BijR_QUAD_to_Dbl(B025R_d,B025I_d,B25R_d,B25I_d,
     &    B025R,B025I,B25R,B25I)
       call R3_From_BijR_QUAD_to_Dbl(B034R_d,B034I_d,B34R_d,B34I_d,
     &    B034R,B034I,B34R,B34I)
       call R3_From_BijR_QUAD_to_Dbl(B035R_d,B035I_d,B35R_d,B35I_d,
     &    B035R,B035I,B35R,B35I)
       call R3_From_BijR_QUAD_to_Dbl(B045R_d,B045I_d,B45R_d,B45I_d,
     &    B045R,B045I,B45R,B45I)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                    C0
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call  R4_From_CijR_QUAD_to_Dbl(C0123R_d,C0123I_d,C123R_d,C123I_d,
     -              C0123R,C0123I,C123R,C123I)
      call  R4_From_CijR_QUAD_to_Dbl(C0124R_d,C0124I_d,C124R_d,C124I_d,
     -              C0124R,C0124I,C124R,C124I)
      call  R4_From_CijR_QUAD_to_Dbl(C0125R_d,C0125I_d,C125R_d,C125I_d,
     -              C0125R,C0125I,C125R,C125I)
      call  R4_From_CijR_QUAD_to_Dbl(C0134R_d,C0134I_d,C134R_d,C134I_d,
     -              C0134R,C0134I,C134R,C134I)
      call  R4_From_CijR_QUAD_to_Dbl(C0135R_d,C0135I_d,C135R_d,C135I_d,
     -              C0135R,C0135I,C135R,C135I)
      call  R4_From_CijR_QUAD_to_Dbl(C0145R_d,C0145I_d,C145R_d,C145I_d,
     -              C0145R,C0145I,C145R,C145I)
      call  R4_From_CijR_QUAD_to_Dbl(C0234R_d,C0234I_d,C234R_d,C234I_d,
     -              C0234R,C0234I,C234R,C234I)
      call  R4_From_CijR_QUAD_to_Dbl(C0235R_d,C0235I_d,C235R_d,C235I_d,
     -              C0235R,C0235I,C235R,C235I)
      call  R4_From_CijR_QUAD_to_Dbl(C0245R_d,C0245I_d,C245R_d,C245I_d,
     -              C0245R,C0245I,C245R,C245I)
      call  R4_From_CijR_QUAD_to_Dbl(C0345R_d,C0345I_d,C345R_d,C345I_d,
     -              C0345R,C0345I,C345R,C345I)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                D0
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call  R5_From_DijR_QUAD_to_Dbl(D01234R_d,D01234I_d,D1234R_d,D1234I_d,
     -            D01234R,D01234I,D1234R,D1234I)
      call  R5_From_DijR_QUAD_to_Dbl(D01235R_d,D01235I_d,D1235R_d,D1235I_d,
     -            D01235R,D01235I,D1235R,D1235I)
      call  R5_From_DijR_QUAD_to_Dbl(D01245R_d,D01245I_d,D1245R_d,D1245I_d,
     -            D01245R,D01245I,D1245R,D1245I)
      call  R5_From_DijR_QUAD_to_Dbl(D01345R_d,D01345I_d,D1345R_d,D1345I_d,
     -            D01345R,D01345I,D1345R,D1345I)
      call  R5_From_DijR_QUAD_to_Dbl(D02345R_d,D02345I_d,D2345R_d,D2345I_d,
     -            D02345R,D02345I,D2345R,D2345I)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                E0
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call  R5_From_EijR_QUAD_to_Dbl(E012345R_d,E012345I_d,E12345R_d,E12345I_d,
     -        E012345R,E012345I,E12345R,E12345I)


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
cfc       If(PrintB) then
cfc       write(*,*)"B45"
cfc       write(*,*) "B45=", B045
cfc       endif
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
cfc       If(PrintB) then
cfc       write(*,*)"B25"
cfc       write(*,*) "B25=", B025
cfc       endif
cfc       If(PrintB) then
cfc       write(*,*)"B35"
cfc       write(*,*) "B35=", B035
cfc       endif
cfc       If(PrintB) then
cfc       write(*,*)"B15"
cfc       write(*,*) "B15=", B015
cfc       endif
cfc       If(PrintC) then
cfc       write(*,*)"C123",C0123R_d,C0123I_d
cfc       
cfc       call print_CijR(DCMPLX(C0123),C123R_d,C123I_d)
cfc       Endif
cfc       If(PrintC) then
cfc       write(*,*)"C124",C0124R_d,C0124I_d
cfc       call print_CijR(DCMPLX(C0124),C124R_d,C124I_d)
cfc       Endif
cfc       If(PrintC) then
cfc       write(*,*)"C125",C0125R_d,C0125I_d
cfc       call print_CijR(DCMPLX(C0125),C125R_d,C125I_d)
cfc       Endif
cfc       If(PrintC) then
cfc       write(*,*)"C234",C0234R_d,C0234I_d
cfc       call print_CijR(DCMPLX(C0234),C234R_d,C234I_d)
cfc       Endif
cfc       If(PrintC) then
cfc       write(*,*)"C235",C0235R_d,C0235I_d
cfc       call print_CijR(DCMPLX(C0235),C235R_d,C235I_d)
cfc       Endif
cfc       If(PrintC) then
cfc       write(*,*)"C345",C0345R_d,C0345I_d
cfc       call print_CijR(DCMPLX(C0345),C345R_d,C345I_d)
cfc       Endif
cfc       If(PrintC) then
cfc       write(*,*)"C134",C0134R_d,C0134I_d
cfc       call print_CijR(DCMPLX(C0134),C134R_d,C134I_d)
cfc       Endif
cfc       If(PrintC) then
cfc       write(*,*)"C135",C0135R_d,C0135I_d
cfc       call print_CijR(DCMPLX(C0135),C135R_d,C135I_d)
cfc       Endif
cfc       If(PrintC) then
cfc       write(*,*)"C145",C0145R_d,C0145I_d
cfc       call print_CijR(DCMPLX(C0145),C145R_d,C145I_d)
cfc       Endif
cfc       If(PrintC) then
cfc       write(*,*)"C245",C0245R_d,C0245I_d
cfc       call print_CijR(DCMPLX(C0245),C245R_d,C245I_d)
cfc       Endif
cfc       If(PrintD) then
cfc       write(*,*)"D1234",D01234R_d,D01234I_d
cfc       write(*,*)"D1234",D01234R,D01234I
cfc       call print_DijR(DCMPLX(D01234),D1234R_d,D1234I_d)
cfc       endif
cfc       If(PrintD) then
cfc       write(*,*)"D1235",D01235R_d,D01235I_d
cfc       write(*,*)"D1235",D01235R,D01235I
cfc       call print_DijR(DCMPLX(D01235),D1235R_d,D1235I_d)
cfc       endif
cfc       If(PrintD) then
cfc       write(*,*)"D1245",D01245R_d,D01245I_d
cfc       write(*,*)"D1245",D01245R,D01245I
cfc       call print_DijR(DCMPLX(D01245),D1245R_d,D1245I_d)
cfc       endif
cfc       If(PrintD) then
cfc       write(*,*)"D2345",D02345R_d,D02345I_d
cfc       write(*,*)"D2345",D02345R,D02345I
cfc       call print_DijR(DCMPLX(D02345),D2345R_d,D2345I_d)
cfc       endif
cfc       If(PrintD) then
cfc       write(*,*)"D1345",D01345R_d,D01345I_d
cfc       write(*,*)"D1345",D01345R,D01345I
cfc       call print_DijR(DCMPLX(D01345),D1345R_d,D1345I_d)
cfc       endif
cfc       If(PrintE) then
cfc       write(*,*)"E12345",E012345R_d,E012345I_d
cfc       write(*,*)"E12345",E012345R,E012345I
cfc       call print_EijR(DCMPLX(E012345),E12345R_d,E12345I_d)
cfc       endif


       return
       end
