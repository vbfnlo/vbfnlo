      REAL*8 FUNCTION MATRIXguuuuH3(PBAR,NHEL,fsign)
C     
C     Generated by MadGraph5_aMC@NLO v. 2.2.3, 2015-02-10
C     By the MadGraph5_aMC@NLO Development Team
C     Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
C     
C     Returns amplitude squared summed/avg over colors
C     for the point with external lines W(0:6,NEXTERNAL)
C     
C     Process: g u > u u u~ h QED=0
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=24)

#include "VBFNLO/utilities/coupl.inc"
#include "ggfloMG.inc"  

      INTEGER    NWAVEFUNCS,     NCOLOR
      PARAMETER (NWAVEFUNCS=14, NCOLOR=4)
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)

      INTEGER NAMPSO, NSQAMPSO
      PARAMETER (NAMPSO=1, NSQAMPSO=1)

      integer    nexternal
      parameter (nexternal=7)

C     
C     ARGUMENTS 
C     
      REAL*8 PBAR(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), fsign(NEXTERNAL)

C     
C     LOCAL VARIABLES 
C     
      INTEGER I,J,M,N
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR,NAMPSO)
      COMPLEX*16 W(18,NWAVEFUNCS)

C     Needed for v4 models
      COMPLEX*16 DUM0,DUM1
      DATA DUM0, DUM1/(0D0, 0D0), (1D0, 0D0)/

C     
C     COLOR DATA
C     
      DATA DENOM(1)/1/
      DATA (CF(I,  1),I=  1,  4) /   12,    4,    4,    0/
C     1 T(1,3,2) T(4,5)
      DATA DENOM(2)/1/
      DATA (CF(I,  2),I=  1,  4) /    4,   12,    0,    4/
C     1 T(1,3,5) T(4,2)
      DATA DENOM(3)/1/
      DATA (CF(I,  3),I=  1,  4) /    4,    0,   12,    4/
C     1 T(1,4,2) T(3,5)
      DATA DENOM(4)/1/
      DATA (CF(I,  4),I=  1,  4) /    0,    4,    4,   12/
C     1 T(1,4,5) T(3,2)
C     ----------
C     BEGIN CODE
C     ----------

      CALL VXXXXX(PBAR(0,1),ZERO,NHEL(1)*FSIGN(1),FSIGN(1),W(1,1))
      CALL IXXXXX(PBAR(0,2),ZERO,NHEL(2)*FSIGN(2),FSIGN(2),W(1,2))

      CALL OXXXXX(PBAR(0,3),ZERO,NHEL(3)*FSIGN(3),FSIGN(3),W(1,3))
      CALL OXXXXX(PBAR(0,4),ZERO,NHEL(4)*FSIGN(4),FSIGN(4),W(1,4))

      CALL IXXXXX(PBAR(0,5),ZERO,NHEL(5)*FSIGN(5),FSIGN(5),W(1,5))
      CALL SXXXXX(PBAR(0,6),FSIGN(6),W(1,6))

cfc      print*
cfc      Print*, "incoming fermion", W(1:6,2)
cfc      print*
cfc 
cfc      Print*, "outgoing fermion 1 " ,W(1:6,3)
cfc      print*
cfc      Print*, "outgoing fermion 2", W(1:6,4)
cfc      print*
cfc      Print*, "out going anti-fermion", W(1:6,5) 
cfc      print*
cfc      Print*, "Higgs", W(1:6,6)
cfc      print*


      CALL FVIXXX(W(1,2),W(1,1),GG,ZERO,ZERO,W(1,7))
      CALL JIOXXX(W(1,5),W(1,3),GG,ZERO,ZERO,W(1,8))
      CALL JIOXXX(W(1,7),W(1,4),GG,ZERO,ZERO,W(1,9))
C     Amplitude(s) for diagram number 1
      CALL VVSHXX(W(1,8),W(1,9),W(1,6),GH,AMP(1))
      CALL JIOXXX(W(1,5),W(1,4),GG,ZERO,ZERO,W(1,9))
      CALL JIOXXX(W(1,7),W(1,3),GG,ZERO,ZERO,W(1,10))
C     Amplitude(s) for diagram number 2
      CALL VVSHXX(W(1,9),W(1,10),W(1,6),GH,AMP(2))
      CALL FVOXXX(W(1,3),W(1,1),GG,ZERO,ZERO,W(1,10))
      CALL JIOXXX(W(1,2),W(1,4),GG,ZERO,ZERO,W(1,7))
      CALL JIOXXX(W(1,5),W(1,10),GG,ZERO,ZERO,W(1,11))
C     Amplitude(s) for diagram number 3
      CALL VVSHXX(W(1,7),W(1,11),W(1,6),GH,AMP(3))
      CALL JIOXXX(W(1,2),W(1,10),GG,ZERO,ZERO,W(1,11))
C     Amplitude(s) for diagram number 4
      CALL VVSHXX(W(1,9),W(1,11),W(1,6),GH,AMP(4))
      CALL FVOXXX(W(1,4),W(1,1),GG,ZERO,ZERO,W(1,11))
      CALL JIOXXX(W(1,2),W(1,3),GG,ZERO,ZERO,W(1,10))
      CALL JIOXXX(W(1,5),W(1,11),GG,ZERO,ZERO,W(1,12))
C     Amplitude(s) for diagram number 5
      CALL VVSHXX(W(1,10),W(1,12),W(1,6),GH,AMP(5))
      CALL JIOXXX(W(1,2),W(1,11),GG,ZERO,ZERO,W(1,12))
C     Amplitude(s) for diagram number 6
      CALL VVSHXX(W(1,8),W(1,12),W(1,6),GH,AMP(6))
      CALL FVIXXX(W(1,5),W(1,1),GG,ZERO,ZERO,W(1,12))
      CALL JIOXXX(W(1,12),W(1,4),GG,ZERO,ZERO,W(1,11))
C     Amplitude(s) for diagram number 7
      CALL VVSHXX(W(1,10),W(1,11),W(1,6),GH,AMP(7))
      CALL JIOXXX(W(1,12),W(1,3),GG,ZERO,ZERO,W(1,11))
C     Amplitude(s) for diagram number 8
      CALL VVSHXX(W(1,7),W(1,11),W(1,6),GH,AMP(8))
      CALL JVSHXX(W(1,1),W(1,6),GH,ZERO,ZERO,W(1,11))
      CALL FVOXXX(W(1,4),W(1,11),GG,ZERO,ZERO,W(1,12))
C     Amplitude(s) for diagram number 9
      CALL IOVXXX(W(1,5),W(1,12),W(1,10),GG,AMP(9))
      CALL FVIXXX(W(1,5),W(1,11),GG,ZERO,ZERO,W(1,13))
C     Amplitude(s) for diagram number 10
      CALL IOVXXX(W(1,13),W(1,4),W(1,10),GG,AMP(10))
C     Amplitude(s) for diagram number 11
      CALL VVVXXX(W(1,11),W(1,10),W(1,9),G,AMP(11))
      CALL FVOXXX(W(1,3),W(1,11),GG,ZERO,ZERO,W(1,14))
C     Amplitude(s) for diagram number 12
      CALL IOVXXX(W(1,5),W(1,14),W(1,7),GG,AMP(12))
C     Amplitude(s) for diagram number 13
      CALL IOVXXX(W(1,13),W(1,3),W(1,7),GG,AMP(13))
C     Amplitude(s) for diagram number 14
      CALL VVVXXX(W(1,11),W(1,7),W(1,8),G,AMP(14))
      CALL FVIXXX(W(1,2),W(1,11),GG,ZERO,ZERO,W(1,13))
C     Amplitude(s) for diagram number 15
      CALL IOVXXX(W(1,13),W(1,4),W(1,8),GG,AMP(15))
C     Amplitude(s) for diagram number 16
      CALL IOVXXX(W(1,2),W(1,12),W(1,8),GG,AMP(16))
C     Amplitude(s) for diagram number 17
      CALL IOVXXX(W(1,13),W(1,3),W(1,9),GG,AMP(17))
C     Amplitude(s) for diagram number 18
      CALL IOVXXX(W(1,2),W(1,14),W(1,9),GG,AMP(18))
C     Amplitude(s) for diagram number 19
      CALL VVVSXX(W(1,1),W(1,10),W(1,9),W(1,6),DUM1,GH4,AMP(19))
      CALL JVVXXX(W(1,1),W(1,10),G,ZERO,ZERO,W(1,14))
C     Amplitude(s) for diagram number 20
      CALL VVSHXX(W(1,9),W(1,14),W(1,6),GH,AMP(20))
      CALL JVVXXX(W(1,1),W(1,9),G,ZERO,ZERO,W(1,14))
C     Amplitude(s) for diagram number 21
      CALL VVSHXX(W(1,10),W(1,14),W(1,6),GH,AMP(21))
C     Amplitude(s) for diagram number 22
      CALL VVVSXX(W(1,1),W(1,7),W(1,8),W(1,6),DUM1,GH4,AMP(22))
      CALL JVVXXX(W(1,1),W(1,7),G,ZERO,ZERO,W(1,14))
C     Amplitude(s) for diagram number 23
      CALL VVSHXX(W(1,8),W(1,14),W(1,6),GH,AMP(23))
      CALL JVVXXX(W(1,1),W(1,8),G,ZERO,ZERO,W(1,14))
C     Amplitude(s) for diagram number 24
      CALL VVSHXX(W(1,7),W(1,14),W(1,6),GH,AMP(24))

c ----------------- GGFLO check -----------------
c      do i=1,6
c         print*,'pbar',i,'=',(pbar(j,i),j=0,3)
c      enddo
c      do i=1,8
c         mgJAMP(i)=JAMP(i)
c         print*,'JAMP(',i,')',JAMP(i)
cfc      enddo
      do i=1,12
         MGamp1(NHEL(1),NHEL(2),NHEL(5),i)=AMP(i)
c         print*,'AMP(',i,')   ',MGamp1(NHEL(1),NHEL(2),NHEL(5),i)
      enddo
      do i=13,24
         MGamp1(NHEL(1),NHEL(2),NHEL(5),i)=AMP(i)
c         print*,'AMP(',i,')   ',MGamp1(NHEL(1),NHEL(2),NHEL(5),i)
      enddo
c------------------------------------------------







C     JAMPs contributing to orders ALL_ORDERS=1

      JAMP(1,1)=+1D0/2D0*(+AMP(1)+1D0/3D0*AMP(2)+AMP(3)+1D0/3D0*AMP(4)
     $ +AMP(12)+AMP(14)+AMP(15)+1D0/3D0*AMP(17)+1D0/3D0*AMP(18)
     $ +AMP(22)+AMP(23)-AMP(24))
      JAMP(2,1)=+1D0/2D0*(-1D0/3D0*AMP(3)-AMP(4)-AMP(7)-1D0/3D0*AMP(8)
     $ -AMP(10)+AMP(11)-1D0/3D0*AMP(12)-1D0/3D0*AMP(13)-AMP(18)
     $ +AMP(19)+AMP(20)-AMP(21))
      JAMP(3,1)=+1D0/2D0*(-1D0/3D0*AMP(1)-AMP(2)-AMP(5)-1D0/3D0*AMP(6)
     $ -AMP(9)-AMP(11)-1D0/3D0*AMP(15)-1D0/3D0*AMP(16)-AMP(17)-AMP(19)
     $ -AMP(20)+AMP(21))
      JAMP(4,1)=+1D0/2D0*(+1D0/3D0*AMP(5)+AMP(6)+1D0/3D0*AMP(7)+AMP(8)
     $ +1D0/3D0*AMP(9)+1D0/3D0*AMP(10)+AMP(13)-AMP(14)+AMP(16)-AMP(22)
     $ -AMP(23)+AMP(24))

      MATRIXguuuuH3 = 0.D0
      M = 1
      DO I = 1, NCOLOR
         ZTEMP = (0.D0,0.D0)
         DO J = 1, NCOLOR
            ZTEMP = ZTEMP + CF(J,I)*JAMP(J,M)
         ENDDO
         N = 1
         MATRIXguuuuH3 = MATRIXguuuuH3 + ZTEMP*DCONJG(JAMP(I,N))/DENOM(I)
      ENDDO

c ----------------- GGFLO check -----------------
cfc      print*,'MATRIXUUUUGH3',MATRIXUUUUGH3
      MGampsq1=MATRIXguuuuH3 
c------------------------------------------------



      END
