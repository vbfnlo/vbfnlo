c This file contain the functions required to evaluate all the matrix elements 
c that appear in the penlinexxx subroutines and in the boxline subroutines.
c They are extension of the function SC1 in the file brakets.f
c They evaluate matrix elements in the quark line with up to 5 insertion of gamma
c matrices. The suffix added to SC are: "1","3" or "5" number of gamma matrices
c "r" the position for a slash(p), and "c" the position of 
c a gamma matrix that has to be contracted with a external current of polarization vector.
c Alpha is the helicity of the spinors.
c Author: Francisco Campanario
c Date: 13/2/2008

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION SC1c(CHII,A1,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC1c
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A1(0:3), AUX(0:3)
   
C
     
      DO I = 0,3
        AUX(I) = A1(I)
      ENDDO
      
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0) - AUX(3)
            ASLASH(1,2) = -AUX(1) + ZI*AUX(2)
            ASLASH(2,1) = -AUX(1) - ZI*AUX(2)
            ASLASH(2,2) = AUX(0) + AUX(3)
         ELSE
            ASLASH(1,1) = AUX(0) + AUX(3)
            ASLASH(1,2) = AUX(1) - ZI*AUX(2)
            ASLASH(2,1) = AUX(1) + ZI*AUX(2)
            ASLASH(2,2) = AUX(0) - AUX(3)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP

C
      SC1c  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION SC1r(CHII,A1,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC1r
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  AUX(0:3)
      REAL*8      A1(0:3)   
C
     
      DO I = 0,3
        AUX(I) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0) - AUX(3)
            ASLASH(1,2) = -AUX(1) + ZI*AUX(2)
            ASLASH(2,1) = -AUX(1) - ZI*AUX(2)
            ASLASH(2,2) = AUX(0) + AUX(3)
         ELSE
            ASLASH(1,1) = AUX(0) + AUX(3)
            ASLASH(1,2) = AUX(1) - ZI*AUX(2)
            ASLASH(2,1) = AUX(1) + ZI*AUX(2)
            ASLASH(2,2) = AUX(0) - AUX(3)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP

C
      SC1r  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
      FUNCTION SC3ccc(CHII,A1,A2,A3,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC3ccc
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A1(0:3), A2(0:3), A3(0:3), AUX(0:3,3)
     
C
      N = 3
      DO I = 0,3
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC3ccc  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION SC3ccr(CHII,A1,A2,A3,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC3ccr
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A2(0:3), A1(0:3), AUX(0:3,3)
      REAL*8  A3(0:3)
C
      N = 3
      DO I = 0,3
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC3ccr  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      FUNCTION SC3crc(CHII,A1,A2,A3,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC3crc
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A1(0:3), AUX(0:3,3),A3(0:3)
      REAL*8  A2(0:3)
C
      N = 3
      DO I = 0,3
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC3crc  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION SC3crr(CHII,A1,A2,A3,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC3crr
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A1(0:3), AUX(0:3,3)
      REAL*8  A2(0:3),A3(0:3)
C
      N = 3
      DO I = 0,3
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC3crr  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION SC3rcc(CHII,A1,A2,A3,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC3rcc
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A2(0:3), A3(0:3), AUX(0:3,3)
      REAL*8  A1(0:3)
C
      N = 3
      DO I = 0,3
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC3rcc  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION SC3rrc(CHII,A1,A2,A3,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC3rrc
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A3(0:3), AUX(0:3,3)
      REAL*8      A1(0:3), A2(0:3)
C
      N = 3
      DO I = 0,3
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC3rrc  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION SC3rrr(CHII,A1,A2,A3,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC3rrr
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  AUX(0:3,3)
      REAL*8      A1(0:3), A2(0:3), A3(0:3)

c      print*, CHII(1),A1(1),A2(1),A3(1),CHIF(1),ALPHA
C
      N = 3
      DO I = 0,3
         AUX(I,3) = DCMPLX(A3(I),0d0)
         AUX(I,2) = DCMPLX(A2(I),0d0)
         AUX(I,1) = DCMPLX(A1(I),0d0)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC3rrr  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION SC5ccccc(CHII,A1,A2,A3,A4,A5,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC5ccccc
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A1(0:3), A2(0:3), A3(0:3), A4(0:3),  A5(0:3), AUX(0:3,5)

     
C
      N = 5
      DO I = 0,3
         AUX(I,5) = A5(I)
         AUX(I,4) = A4(I)
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC5ccccc  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION SC5ccccr(CHII,A1,A2,A3,A4,A5,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC5ccccr
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A1(0:3), A2(0:3), A3(0:3), A4(0:3), AUX(0:3,5)
      REAL*8   A5(0:3)
     
C
      N = 5
      DO I = 0,3
         AUX(I,5) = A5(I)
         AUX(I,4) = A4(I)
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC5ccccr  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION SC5cccrc(CHII,A1,A2,A3,A4,A5,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC5cccrc
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A1(0:3), A2(0:3),A3(0:3), A5(0:3) ,AUX(0:3,5)
      REAL*8    A4(0:3)
     
C
      N = 5
      DO I = 0,3
         AUX(I,5) = A5(I)
         AUX(I,4) = A4(I)
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC5cccrc  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION SC5cccrr(CHII,A1,A2,A3,A4,A5,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC5cccrr
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A1(0:3), A2(0:3), A3(0:3),  AUX(0:3,5)
      REAL*8   A4(0:3),A5(0:3)
     
C
      N = 5
      DO I = 0,3
         AUX(I,5) = A5(I)
         AUX(I,4) = A4(I)
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC5cccrr  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION SC5ccrrc(CHII,A1,A2,A3,A4,A5,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC5ccrrc
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A1(0:3), A2(0:3), A5(0:3) ,AUX(0:3,5)
      REAL*8    A3(0:3),A4(0:3)
     
C
      N = 5
      DO I = 0,3
         AUX(I,5) = A5(I)
         AUX(I,4) = A4(I)
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC5ccrrc  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION SC5ccrrr(CHII,A1,A2,A3,A4,A5,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC5ccrrr
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A1(0:3), A2(0:3),  AUX(0:3,5)
      REAL*8    A3(0:3),A4(0:3),A5(0:3)
     
C
      N = 5
      DO I = 0,3
         AUX(I,5) = A5(I)
         AUX(I,4) = A4(I)
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC5ccrrr  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       FUNCTION SC5crccc(CHII,A1,A2,A3,A4,A5,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC5crccc
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A1(0:3), A3(0:3), A5(0:3), A4(0:3), AUX(0:3,5)
      REAL*8   A2(0:3)
     
C
      N = 5
      DO I = 0,3
         AUX(I,5) = A5(I)
         AUX(I,4) = A4(I)
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC5crccc  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       FUNCTION SC5crcrc(CHII,A1,A2,A3,A4,A5,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC5crcrc
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A1(0:3), A3(0:3), A5(0:3), AUX(0:3,5)
      REAL*8   A2(0:3), A4(0:3)
     
C
      N = 5
      DO I = 0,3
         AUX(I,5) = A5(I)
         AUX(I,4) = A4(I)
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC5crcrc  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION SC5crrrc(CHII,A1,A2,A3,A4,A5,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC5crrrc
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A1(0:3), A5(0:3) ,AUX(0:3,5)
      REAL*8    A2(0:3), A3(0:3),A4(0:3)
     
C
      N = 5
      DO I = 0,3
         AUX(I,5) = A5(I)
         AUX(I,4) = A4(I)
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC5crrrc  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION SC5rrccc(CHII,A1,A2,A3,A4,A5,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC5rrccc
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A3(0:3), A4(0:3), A5(0:3), AUX(0:3,5)
      REAL*8      A1(0:3), A2(0:3)
     
C
      N = 5
      DO I = 0,3
         AUX(I,5) = A5(I)
         AUX(I,4) = A4(I)
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC5rrccc  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      FUNCTION SC5rrcrc(CHII,A1,A2,A3,A4,A5,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC5rrcrc
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A3(0:3), A5(0:3), AUX(0:3,5)
      REAL*8      A1(0:3), A2(0:3), A4(0:3)
     
C
      N = 5
      DO I = 0,3
         AUX(I,5) = A5(I)
         AUX(I,4) = A4(I)
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC5rrcrc  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      FUNCTION SC7ccccrrr(CHII,A1,A2,A3,A4,A5,A6,A7,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC7ccccrrr
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A1(0:3), A2(0:3), A3(0:3), A4(0:3), AUX(0:3,7)
      REAL*8      A5(0:3), A6(0:3), A7(0:3)
     
C
      N = 7
      DO I = 0,3
         AUX(I,7) = A7(I)
         AUX(I,6) = A6(I)
         AUX(I,5) = A5(I)
         AUX(I,4) = A4(I)
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC7ccccrrr  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION SC7ccrrccc(CHII,A1,A2,A3,A4,A5,A6,A7,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC7ccrrccc
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A1(0:3), A2(0:3), A5(0:3), A6(0:3), A7(0:3), AUX(0:3,7)
      REAL*8      A3(0:3), A4(0:3)
     
C
      N = 7
      DO I = 0,3
         AUX(I,7) = A7(I)
         AUX(I,6) = A6(I)
         AUX(I,5) = A5(I)
         AUX(I,4) = A4(I)
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC7ccrrccc  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION SC7ccrrcrc(CHII,A1,A2,A3,A4,A5,A6,A7,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC7ccrrcrc
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A1(0:3), A2(0:3), A5(0:3), A7(0:3), AUX(0:3,7)
      REAL*8      A3(0:3), A4(0:3), A6(0:3)
     
C
      N = 7
      DO I = 0,3
         AUX(I,7) = A7(I)
         AUX(I,6) = A6(I)
         AUX(I,5) = A5(I)
         AUX(I,4) = A4(I)
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC7ccrrcrc  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END


