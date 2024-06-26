      REAL*8 FUNCTION MATRIXubuHubug3(PBAR,NHEL,fsign)
C     
C     Generated by MadGraph5_aMC@NLO v. 2.2.3, 2015-02-10
C     By the MadGraph5_aMC@NLO Development Team
C     Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
C     
C     Returns amplitude squared summed/avg over colors
C     for the point with external lines W(0:6,NEXTERNAL)
C     
C     Process: u u~ > u u~ g h QED=0
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
      integer    nexternal
      parameter (nexternal=7)

      INTEGER NAMPSO, NSQAMPSO
      PARAMETER (NAMPSO=1, NSQAMPSO=1)

C     
C     ARGUMENTS 
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
C     1 T(2,1) T(5,3,4)
      DATA DENOM(2)/1/
      DATA (CF(I,  2),I=  1,  4) /    4,   12,    0,    4/
C     1 T(2,4) T(5,3,1)
      DATA DENOM(3)/1/
      DATA (CF(I,  3),I=  1,  4) /    4,    0,   12,    4/
C     1 T(3,1) T(5,2,4)
      DATA DENOM(4)/1/
      DATA (CF(I,  4),I=  1,  4) /    0,    4,    4,   12/
C     1 T(3,4) T(5,2,1)
C     ----------
C     BEGIN CODE
C     ----------
      CALL IXXXXX(PBAR(0,1),ZERO,NHEL(1)*FSIGN(1),FSIGN(1),W(1,1))
      CALL IXXXXX(PBAR(0,2),ZERO,NHEL(2)*FSIGN(2),FSIGN(2),W(1,4))

      CALL OXXXXX(PBAR(0,3),ZERO,NHEL(3)*FSIGN(3),FSIGN(3),W(1,2))
      CALL OXXXXX(PBAR(0,4),ZERO,NHEL(4)*FSIGN(4),FSIGN(4),W(1,3))

      CALL VXXXXX(PBAR(0,5),ZERO,NHEL(5)*FSIGN(5),FSIGN(5),W(1,5))
      CALL SXXXXX(PBAR(0,6),FSIGN(6),W(1,6))

cfc      Print*, "Incoming fermion 1 " ,W(1:6,1)
cfc      print*
cfc      Print*, "outgoing fermion", W(1:6,3)
cfc      print*
cfc 
cfc      Print*, "Incoming anti-fermion 1 " ,W(1:6,2)
cfc      print*
cfc      Print*, "outgoing anti-fermion", W(1:6,4)
cfc      print*
cfc      Print*, " gluon polarization", W(1:6,5) 
cfc      print*
cfc      Print*, "Higgs", W(1:6,6)
cfc      print*

      CALL JIOXXX(W(1,1),W(1,2),GG,ZERO,ZERO,W(1,7))
      CALL JIOXXX(W(1,4),W(1,3),GG,ZERO,ZERO,W(1,8))
      
cfc      Print*, "Current J12", W(1:6,7)
cfc      print*
cfc      Print*, "Current J34", W(1:6,8)
cfc      print*

C     Amplitude(s) for diagram number 1
      CALL VVVSXX(W(1,7),W(1,8),W(1,5),W(1,6),DUM1,GH4,AMP(1))
      CALL JVVXXX(W(1,7),W(1,5),G,ZERO,ZERO,W(1,9))
C     Amplitude(s) for diagram number 2
      CALL VVSHXX(W(1,8),W(1,9),W(1,6),GH,AMP(2))
      CALL JVSHXX(W(1,7),W(1,6),GH,ZERO,ZERO,W(1,9))
C     Amplitude(s) for diagram number 3
      CALL VVVXXX(W(1,8),W(1,5),W(1,9),G,AMP(3))
      CALL JVSHXX(W(1,5),W(1,6),GH,ZERO,ZERO,W(1,10))
C     Amplitude(s) for diagram number 4
      CALL VVVXXX(W(1,7),W(1,8),W(1,10),G,AMP(4))
      CALL FVOXXX(W(1,3),W(1,5),GG,ZERO,ZERO,W(1,11))
C     Amplitude(s) for diagram number 5
      CALL IOVXXX(W(1,4),W(1,11),W(1,9),GG,AMP(5))
      CALL FVIXXX(W(1,4),W(1,5),GG,ZERO,ZERO,W(1,12))
C     Amplitude(s) for diagram number 6
      CALL IOVXXX(W(1,12),W(1,3),W(1,9),GG,AMP(6))
      CALL FVOXXX(W(1,3),W(1,7),GG,ZERO,ZERO,W(1,9))
C     Amplitude(s) for diagram number 7
      CALL IOVXXX(W(1,4),W(1,9),W(1,10),GG,AMP(7))
      CALL FVIXXX(W(1,4),W(1,7),GG,ZERO,ZERO,W(1,9))
C     Amplitude(s) for diagram number 8
      CALL IOVXXX(W(1,9),W(1,3),W(1,10),GG,AMP(8))
      CALL JIOXXX(W(1,1),W(1,3),GG,ZERO,ZERO,W(1,9))
      CALL JIOXXX(W(1,4),W(1,2),GG,ZERO,ZERO,W(1,7))
C     Amplitude(s) for diagram number 9
      CALL VVVSXX(W(1,9),W(1,7),W(1,5),W(1,6),DUM1,GH4,AMP(9))
      CALL JVVXXX(W(1,9),W(1,5),G,ZERO,ZERO,W(1,13))
C     Amplitude(s) for diagram number 10
      CALL VVSHXX(W(1,7),W(1,13),W(1,6),GH,AMP(10))
      CALL JVSHXX(W(1,9),W(1,6),GH,ZERO,ZERO,W(1,13))
C     Amplitude(s) for diagram number 11
      CALL VVVXXX(W(1,7),W(1,5),W(1,13),G,AMP(11))
C     Amplitude(s) for diagram number 12
      CALL VVVXXX(W(1,9),W(1,7),W(1,10),G,AMP(12))
      CALL FVOXXX(W(1,2),W(1,5),GG,ZERO,ZERO,W(1,14))
C     Amplitude(s) for diagram number 13
      CALL IOVXXX(W(1,4),W(1,14),W(1,13),GG,AMP(13))
C     Amplitude(s) for diagram number 14
      CALL IOVXXX(W(1,12),W(1,2),W(1,13),GG,AMP(14))
      CALL FVOXXX(W(1,2),W(1,9),GG,ZERO,ZERO,W(1,13))
C     Amplitude(s) for diagram number 15
      CALL IOVXXX(W(1,4),W(1,13),W(1,10),GG,AMP(15))
      CALL FVIXXX(W(1,4),W(1,9),GG,ZERO,ZERO,W(1,13))
C     Amplitude(s) for diagram number 16
      CALL IOVXXX(W(1,13),W(1,2),W(1,10),GG,AMP(16))
      CALL FVIXXX(W(1,1),W(1,5),GG,ZERO,ZERO,W(1,13))
      CALL JIOXXX(W(1,13),W(1,3),GG,ZERO,ZERO,W(1,5))
C     Amplitude(s) for diagram number 17
      CALL VVSHXX(W(1,7),W(1,5),W(1,6),GH,AMP(17))
      CALL JIOXXX(W(1,13),W(1,2),GG,ZERO,ZERO,W(1,5))
C     Amplitude(s) for diagram number 18
      CALL VVSHXX(W(1,8),W(1,5),W(1,6),GH,AMP(18))
      CALL JIOXXX(W(1,1),W(1,11),GG,ZERO,ZERO,W(1,5))
C     Amplitude(s) for diagram number 19
      CALL VVSHXX(W(1,7),W(1,5),W(1,6),GH,AMP(19))
      CALL FVIXXX(W(1,1),W(1,7),GG,ZERO,ZERO,W(1,5))
C     Amplitude(s) for diagram number 20
      CALL IOVXXX(W(1,5),W(1,3),W(1,10),GG,AMP(20))
      CALL FVIXXX(W(1,1),W(1,10),GG,ZERO,ZERO,W(1,5))
C     Amplitude(s) for diagram number 21
      CALL IOVXXX(W(1,5),W(1,3),W(1,7),GG,AMP(21))
      CALL JIOXXX(W(1,1),W(1,14),GG,ZERO,ZERO,W(1,7))
C     Amplitude(s) for diagram number 22
      CALL VVSHXX(W(1,8),W(1,7),W(1,6),GH,AMP(22))
      CALL FVIXXX(W(1,1),W(1,8),GG,ZERO,ZERO,W(1,7))
C     Amplitude(s) for diagram number 23
      CALL IOVXXX(W(1,7),W(1,2),W(1,10),GG,AMP(23))
C     Amplitude(s) for diagram number 24
      CALL IOVXXX(W(1,5),W(1,2),W(1,8),GG,AMP(24))
C     JAMPs contributing to orders ALL_ORDERS=1


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


      JAMP(1,1)=+1D0/2D0*(+1D0/3D0*AMP(5)+1D0/3D0*AMP(6)+1D0/3D0
     $ *AMP(7)+1D0/3D0*AMP(8)-AMP(9)+AMP(10)-AMP(11)-AMP(12)+AMP(14)
     $ +AMP(15)+AMP(19)+AMP(20))
      JAMP(2,1)=+1D0/2D0*(-AMP(1)+AMP(2)-AMP(3)-AMP(4)-AMP(5)-AMP(8)
     $ -1D0/3D0*AMP(17)-AMP(18)-1D0/3D0*AMP(19)-1D0/3D0*AMP(20)
     $ -1D0/3D0*AMP(21)-AMP(24))
      JAMP(3,1)=+1D0/2D0*(+AMP(1)-AMP(2)+AMP(3)+AMP(4)-AMP(6)-AMP(7)
     $ -1D0/3D0*AMP(13)-1D0/3D0*AMP(14)-1D0/3D0*AMP(15)-1D0/3D0
     $ *AMP(16)-AMP(22)-AMP(23))
      JAMP(4,1)=+1D0/2D0*(+AMP(9)-AMP(10)+AMP(11)+AMP(12)+AMP(13)
     $ +AMP(16)+AMP(17)+1D0/3D0*AMP(18)+AMP(21)+1D0/3D0*AMP(22)
     $ +1D0/3D0*AMP(23)+1D0/3D0*AMP(24))

      MATRIXubuHubug3 = 0.D0
      M = 1
      DO I = 1, NCOLOR
         ZTEMP = (0.D0,0.D0)
         DO J = 1, NCOLOR
            ZTEMP = ZTEMP + CF(J,I)*JAMP(J,M)
         ENDDO
         N = 1
         MATRIXubuHubug3 = MATRIXubuHubug3 + ZTEMP*DCONJG(JAMP(I,N))/DENOM(I)
      ENDDO

C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
c ----------------- GGFLO check -----------------
c      print*,'MATRIXucucgH2',MATRIXucucgH2
      MGampsq1=MATRIXubuHubug3
c------------------------------------------------
c      STOP

      END
