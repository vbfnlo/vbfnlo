

      SUBROUTINE SDU_WMZgL(PBAR,fsign,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : u d~ -> ve e+ mu- mu+ g  
C  
C Crossing   1 is u d~ -> ve e+ mu- mu+ g  
      IMPLICIT NONE
C  
C CONSTANTS
C  
#include "VBFNLO/utilities/global.inc"
#include "genpsr.inc"
#include "tensor.inc"
!       Include "genps.inc"
!       Include "nexternal.inc"
!       Include "maxamps.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB= 128, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 PBAR(0:3,NEXTERNAL),ANS(NCROSS)
      integer fsign(7)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T     !, P(0:3,NEXTERNAL)
      REAL*8 DUWMZgL
      INTEGER IHEL,IDEN(NCROSS)!,IC(NEXTERNAL,NCROSS)
      INTEGER IPROC, I !,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2(maxamps), jamp2(0:maxamps)
      common/to_amps26/  amp2,       jamp2

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrix26/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
!       DATA warned, isum_hel/.false.,0/
      DATA warned /.false./
!       DATA multi_channel/.true./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /   24/          
      DATA jamp2(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 7) /-1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 7) /-1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 7) /-1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 7) /-1,-1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 7) /-1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 7) /-1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 7) /-1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 7) /-1,-1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 7) /-1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 7) /-1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 7) /-1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 7) /-1,-1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 7) /-1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 7) /-1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 7) /-1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 7) /-1,-1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 7) /-1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 7) /-1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 7) /-1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 7) /-1,-1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 7) /-1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 7) /-1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 7) /-1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 7) /-1,-1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 7) /-1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 7) /-1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 7) /-1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 7) /-1,-1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 7) /-1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 7) /-1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 7) /-1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 7) /-1,-1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 7) /-1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 7) /-1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 7) /-1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 7) /-1, 1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 7) /-1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 7) /-1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 7) /-1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 7) /-1, 1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 7) /-1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 7) /-1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 7) /-1, 1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 7) /-1, 1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 7) /-1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 7) /-1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 7) /-1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 7) /-1, 1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  49),IHEL=1, 7) /-1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  50),IHEL=1, 7) /-1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  51),IHEL=1, 7) /-1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  52),IHEL=1, 7) /-1, 1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  53),IHEL=1, 7) /-1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  54),IHEL=1, 7) /-1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  55),IHEL=1, 7) /-1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  56),IHEL=1, 7) /-1, 1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  57),IHEL=1, 7) /-1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  58),IHEL=1, 7) /-1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  59),IHEL=1, 7) /-1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  60),IHEL=1, 7) /-1, 1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  61),IHEL=1, 7) /-1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  62),IHEL=1, 7) /-1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  63),IHEL=1, 7) /-1, 1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  64),IHEL=1, 7) /-1, 1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  65),IHEL=1, 7) / 1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  66),IHEL=1, 7) / 1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  67),IHEL=1, 7) / 1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  68),IHEL=1, 7) / 1,-1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  69),IHEL=1, 7) / 1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  70),IHEL=1, 7) / 1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  71),IHEL=1, 7) / 1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  72),IHEL=1, 7) / 1,-1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  73),IHEL=1, 7) / 1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  74),IHEL=1, 7) / 1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  75),IHEL=1, 7) / 1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  76),IHEL=1, 7) / 1,-1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  77),IHEL=1, 7) / 1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  78),IHEL=1, 7) / 1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  79),IHEL=1, 7) / 1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  80),IHEL=1, 7) / 1,-1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  81),IHEL=1, 7) / 1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  82),IHEL=1, 7) / 1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  83),IHEL=1, 7) / 1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  84),IHEL=1, 7) / 1,-1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  85),IHEL=1, 7) / 1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  86),IHEL=1, 7) / 1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  87),IHEL=1, 7) / 1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  88),IHEL=1, 7) / 1,-1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  89),IHEL=1, 7) / 1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  90),IHEL=1, 7) / 1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  91),IHEL=1, 7) / 1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  92),IHEL=1, 7) / 1,-1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  93),IHEL=1, 7) / 1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  94),IHEL=1, 7) / 1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  95),IHEL=1, 7) / 1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  96),IHEL=1, 7) / 1,-1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  97),IHEL=1, 7) / 1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  98),IHEL=1, 7) / 1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  99),IHEL=1, 7) / 1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 100),IHEL=1, 7) / 1, 1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 101),IHEL=1, 7) / 1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 102),IHEL=1, 7) / 1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 103),IHEL=1, 7) / 1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 104),IHEL=1, 7) / 1, 1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 105),IHEL=1, 7) / 1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 106),IHEL=1, 7) / 1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 107),IHEL=1, 7) / 1, 1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 108),IHEL=1, 7) / 1, 1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 109),IHEL=1, 7) / 1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 110),IHEL=1, 7) / 1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 111),IHEL=1, 7) / 1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 112),IHEL=1, 7) / 1, 1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 113),IHEL=1, 7) / 1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 114),IHEL=1, 7) / 1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 115),IHEL=1, 7) / 1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 116),IHEL=1, 7) / 1, 1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 117),IHEL=1, 7) / 1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 118),IHEL=1, 7) / 1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 119),IHEL=1, 7) / 1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 120),IHEL=1, 7) / 1, 1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 121),IHEL=1, 7) / 1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 122),IHEL=1, 7) / 1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 123),IHEL=1, 7) / 1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 124),IHEL=1, 7) / 1, 1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 125),IHEL=1, 7) / 1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 126),IHEL=1, 7) / 1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 127),IHEL=1, 7) / 1, 1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 128),IHEL=1, 7) / 1, 1, 1, 1, 1, 1, 1/
!       DATA (  IC(IHEL,  1),IHEL=1, 7) / 1, 2, 3, 4, 5, 6, 7/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  1/
C ----------
C BEGIN CODE
C ----------
c neu
      iproc=1
c -------
      NTRY=NTRY+1
!       DO IPROC=1,NCROSS
!       CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
!       DO IHEL=1,NEXTERNAL
!          JC(IHEL) = +1
!       ENDDO
!        
!       IF (multi_channel) THEN
!           DO IHEL=1,NGRAPHS
!               amp2(ihel)=0d0
!           ENDDO
!           DO IHEL=1,int(jamp2(0))
!               jamp2(ihel)=0d0
!           ENDDO
!       ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=DUWMZgL(PBAR ,NHEL(1,IHEL),FSIGN(1))            
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=DUWMZgL(PBAR ,NHEL(1,IHEL),fsign(1))            
            ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2(MAPCONFIG(IHEL))
          ENDDO
           ANS(IPROC)=0
          
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
!       ENDDO
      END
       
       
      REAL*8 FUNCTION DUWMZgL(PBAR,NHEL,FSIGN)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : u d~ -> ve e+ mu- mu+ g  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=  24,NEIGEN=  1) 
#include "VBFNLO/utilities/global.inc"
#include "genpsr.inc"
#include "tensor.inc"
!       include "genps.inc"
!       include "nexternal.inc"
!       include "maxamps.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  46, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 PBAR(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), fsign(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)

      DOUBLE COMPLEX subprocess(3)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2(maxamps), jamp2(0:maxamps)
      common/to_amps26/  amp2,       jamp2
!       include "coupl.inc"
#include "VBFNLO/utilities/coupl.inc"

      integer h
      common /hcount / h
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     4/                                  
C               T[ 2, 1, 7]                                                
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(PBAR(0,1),ZERO ,NHEL(1)*FSIGN(1),FSIGN(1),W(1,1))        
      CALL OXXXXX(PBAR(0,2),ZERO ,NHEL(2)*FSIGN(2),FSIGN(2),W(1,2))        
      CALL OXXXXX(PBAR(0,3),ZERO ,NHEL(3)*FSIGN(3),FSIGN(3),W(1,3))        
      CALL IXXXXX(PBAR(0,4),ZERO ,NHEL(4)*FSIGN(4),FSIGN(4),W(1,4))        
      CALL OXXXXX(PBAR(0,5),ZERO ,NHEL(5)*FSIGN(5),FSIGN(5),W(1,5))        
      CALL IXXXXX(PBAR(0,6),ZERO ,NHEL(6)*FSIGN(6),FSIGN(6),W(1,6))        
      CALL VXXXXX(PBAR(0,7),ZERO ,NHEL(7)*FSIGN(7),FSIGN(7),W(1,7))  
      CALL JIOXXX(W(1,4   ),W(1,3   ),GWF ,WMASS   ,WWIDTH  ,W(1,8   ))    
      CALL FVIXXX(W(1,1   ),W(1,7   ),GG ,ZERO    ,ZERO    ,W(1,9   ))     
      CALL FVOXXX(W(1,2   ),W(1,8   ),GWF ,ZERO    ,ZERO    ,W(1,10  ))    
      CALL JIOXXX(W(1,9   ),W(1,10  ),GAD ,ZERO    ,AWIDTH  ,W(1,11  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,11  ),GAL ,AMP(1   ))            
      CALL JIOXXX(W(1,9   ),W(1,10  ),GZD ,ZMASS   ,ZWIDTH  ,W(1,12  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,12  ),GZL ,AMP(2   ))            
      CALL FVOXXX(W(1,2   ),W(1,7   ),GG ,ZERO    ,ZERO    ,W(1,13  ))     
      CALL FVOXXX(W(1,13  ),W(1,8   ),GWF ,ZERO    ,ZERO    ,W(1,14  ))    
      CALL JIOXXX(W(1,1   ),W(1,14  ),GAD ,ZERO    ,AWIDTH  ,W(1,15  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,15  ),GAL ,AMP(3   ))            
      CALL JIOXXX(W(1,1   ),W(1,14  ),GZD ,ZMASS   ,ZWIDTH  ,W(1,16  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,16  ),GZL ,AMP(4   ))            
      CALL FVOXXX(W(1,10  ),W(1,7   ),GG ,ZERO    ,ZERO    ,W(1,17  ))     
      CALL JIOXXX(W(1,1   ),W(1,17  ),GAD ,ZERO    ,AWIDTH  ,W(1,18  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,18  ),GAL ,AMP(5   ))            
      CALL JIOXXX(W(1,1   ),W(1,17  ),GZD ,ZMASS   ,ZWIDTH  ,W(1,19  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,19  ),GZL ,AMP(6   ))            
      CALL FVIXXX(W(1,9   ),W(1,8   ),GWF ,ZERO    ,ZERO    ,W(1,20  ))    
      CALL JIOXXX(W(1,20  ),W(1,2   ),GAU ,ZERO    ,AWIDTH  ,W(1,21  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,21  ),GAL ,AMP(7   ))            
      CALL JIOXXX(W(1,20  ),W(1,2   ),GZU ,ZMASS   ,ZWIDTH  ,W(1,22  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,22  ),GZL ,AMP(8   ))            
      CALL FVIXXX(W(1,1   ),W(1,8   ),GWF ,ZERO    ,ZERO    ,W(1,23  ))    
      CALL JIOXXX(W(1,23  ),W(1,13  ),GAU ,ZERO    ,AWIDTH  ,W(1,24  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,24  ),GAL ,AMP(9   ))            
      CALL JIOXXX(W(1,23  ),W(1,13  ),GZU ,ZMASS   ,ZWIDTH  ,W(1,25  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,25  ),GZL ,AMP(10  ))            
      CALL FVIXXX(W(1,23  ),W(1,7   ),GG ,ZERO    ,ZERO    ,W(1,26  ))     
      CALL JIOXXX(W(1,26  ),W(1,2   ),GAU ,ZERO    ,AWIDTH  ,W(1,27  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,27  ),GAL ,AMP(11  ))            
      CALL JIOXXX(W(1,26  ),W(1,2   ),GZU ,ZMASS   ,ZWIDTH  ,W(1,28  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,28  ),GZL ,AMP(12  ))            
      CALL JIOXXX(W(1,9   ),W(1,2   ),GWF ,WMASS   ,WWIDTH  ,W(1,29  ))    
      CALL FVIXXX(W(1,4   ),W(1,29  ),GWF ,ZERO    ,ZERO    ,W(1,30  ))    
      CALL JIOXXX(W(1,30  ),W(1,3   ),GAL ,ZERO    ,AWIDTH  ,W(1,31  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,31  ),GAL ,AMP(13  ))            
      CALL JIOXXX(W(1,30  ),W(1,3   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,32  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,32  ),GZL ,AMP(14  ))            
      CALL JIOXXX(W(1,1   ),W(1,13  ),GWF ,WMASS   ,WWIDTH  ,W(1,33  ))    
      CALL FVIXXX(W(1,4   ),W(1,33  ),GWF ,ZERO    ,ZERO    ,W(1,34  ))    
      CALL JIOXXX(W(1,34  ),W(1,3   ),GAL ,ZERO    ,AWIDTH  ,W(1,35  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,35  ),GAL ,AMP(15  ))            
      CALL JIOXXX(W(1,34  ),W(1,3   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,36  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,36  ),GZL ,AMP(16  ))            
      CALL FVOXXX(W(1,3   ),W(1,29  ),GWF ,ZERO    ,ZERO    ,W(1,37  ))    
      CALL JIOXXX(W(1,4   ),W(1,37  ),GZN ,ZMASS   ,ZWIDTH  ,W(1,38  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,38  ),GZL ,AMP(17  ))            
      CALL FVOXXX(W(1,3   ),W(1,33  ),GWF ,ZERO    ,ZERO    ,W(1,39  ))    
      CALL JIOXXX(W(1,4   ),W(1,39  ),GZN ,ZMASS   ,ZWIDTH  ,W(1,40  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,40  ),GZL ,AMP(18  ))            
      CALL JVVXXX(W(1,8   ),W(1,29  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     41  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,41  ),GAL ,AMP(19  ))            
      CALL JVVXXX(W(1,8   ),W(1,29  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     42  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,42  ),GZL ,AMP(20  ))            
      CALL JVVXXX(W(1,8   ),W(1,33  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     43  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,43  ),GAL ,AMP(21  ))            
      CALL JVVXXX(W(1,8   ),W(1,33  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     44  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,44  ),GZL ,AMP(22  ))            
      CALL FVOXXX(W(1,5   ),W(1,29  ),GWF ,ZERO    ,ZERO    ,W(1,45  ))    
      CALL IOVXXX(W(1,6   ),W(1,45  ),W(1,8   ),GWF ,AMP(23  ))            
      CALL FVOXXX(W(1,5   ),W(1,33  ),GWF ,ZERO    ,ZERO    ,W(1,46  ))    
      CALL IOVXXX(W(1,6   ),W(1,46  ),W(1,8   ),GWF ,AMP(24  )) 
      subprocess(1) =   -AMP(   1)-AMP(   2)-AMP(   3)-AMP(   4)
     &                  -AMP(   5)-AMP(   6)-AMP(   7)-AMP(   8)
     &                  -AMP(   9)-AMP(  10)-AMP(  11)-AMP(  12)

!         print*,subprocess(1)

      subprocess(2) =   -AMP(  13)-AMP(  14)-AMP(  15)-AMP(  16) 
     &                  -AMP(  17)-AMP(  18)-AMP(  19)-AMP(  20)
     &                  -AMP(  21)-AMP(  22)-AMP(  23)-AMP(  24)

       JAMP(   1) =  subprocess(1) + subprocess(2) 
c------------------------------------------------------------------
c comparison with hand made code
      do i = 1,2
         if ((nhel(1).eq.nhel(2)).and.
     2        (nhel(3).eq.-1).and.(nhel(4).eq.-1)) then
           if((h.eq.1.and.nhel(5).eq.-1.and.nhel(6).eq.-1).or.
     1        (h.eq.2.and.nhel(5).eq.+1.and.nhel(6).eq.+1)) then
!              if((h.eq.1.and.nhel(7).eq.1).or.
!      1          (h.eq.2.and.nhel(7).eq.1).or.
!      2          (h.eq.3.and.nhel(7).eq.-1) .or.
!      3          (h.eq.4.and.nhel(7).eq.-1)) then
               tampr(i,nhel(1),nhel(7),1) 
     &            = subprocess(i)/g ! gluon pol. in helicity basis
!                print*,tampr(i,nhel(1),nhel(7),1)
!              endif 
           endif
         endif
      enddo ! divide out strong coupling
c------------------------------------------------------------------


      DUWMZgL = 0.D0 
  
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          DUWMZgL =DUWMZgL+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
!           print*, "DUWMZgLbla",DUWMZgL
      ENDDO
      Do I = 1, NGRAPHS
          amp2(i)=amp2(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2(i)=Jamp2(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END



