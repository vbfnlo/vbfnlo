! This subroutine returns the full hadronic matrix element squared for the process
! p p --> e- ve gamma jet at leading order. It basically uses the MADGRAPH-generated 
! code which is subsequently interfaced with vbfnlo's phasespace


      function m2s_qqwzglu_hel(xi,p,v,nlo,pdfchange)

      use globalvars, only: ldoblha

      implicit none


#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

      integer nlo,k
      logical pdfchange
      double precision xi(nx), p(0:3,max_p), v(0:3,max_v),m2s_qqwzglu_hel

      double precision helastestmom1(0:3,7),m2helastest(1),pdf(-6:6,2)
      double precision m2s_qqwzglu,m2s_wgaglu

      double precision helasmom(0:3,8) !blha momenta
      ! -- one entry more due to additional unused entry in MG
      double precision res(1) !blha amplitude

      double precision dotrr

      external m2s_qqwzglu,m2s_wgaglu,dotrr

      INTEGER init/0/
      SAVE init

      call getyourscalesready(xi,p(0,1),v(0,1),pdfchange,pdf,nlo,1)!.true.)

      if ( init .eq. 0 .and. nlo.eq.0 .and. .not. ldoblha) then
        if(.not.pdfchange) then

          write(6,*) ' '
          write(6,*) 'W-Zj amplitude square information:'
          write(6,*) '-----------------------------------------------'
          write(6,*) ' '
          write(6,*) ' '
          write(6,*) '-----------------------------------------------'

          call printnfl(.true.)
        endif
      endif
      init=1

      m2s_qqwzglu_hel=0d0

!
! Assign Momenta
!
      do k=0,3
         helastestmom1(k,1)=p(k,1)
         helastestmom1(k,2)=p(k,2)
         helastestmom1(k,3)=p(k,3)
         helastestmom1(k,4)=v(k,1)
         helastestmom1(k,5)=v(k,2)
         helastestmom1(k,6)=v(k,3)
         helastestmom1(k,7)=v(k,4)
      enddo

         !TODO implement blha_ranhelcomb, helicity argument to SUDBGEPVEMUPMUM 
      if (ldoblha) then
        helasmom(0:3,1:3) = p(0:3, 1:3)
        helasmom(0:3,4:7) = v(0:3, 1:4) 

        if (blha_fsign(3,blha_cursubproc,blha_curproc) == 1) then ! ubar d -> g
          if (blha_physdiagmap(1,blha_cursubproc,blha_curproc) == 1) then ! d ubar -> g
            helasmom(0:3,1) = p(0:3,2)
            helasmom(0:3,2) = p(0:3,1)
          endif
          Call SUBDGEMVEBMUMMUP(helasmom, res)
          res = res * 36d0 / 2d0
        else 
          if (blha_physdiagmap(1,blha_cursubproc,blha_curproc) == 3) then ! g d -> u
            helasmom(0:3,1) = p(0:3,2)
            helasmom(0:3,2) = p(0:3,1)
          endif
          if (blha_fsign(1,blha_cursubproc,blha_curproc) == 1) then ! d g -> u
            Call SDGUBEMVEBMUMMUP(helasmom, res)
            res = res * 96d0 / 2d0
          else ! dbar g -> ubar
             Call SUBGDBEMVEBMUMMUP(helasmom, res)
             res = res * 96d0 / 2d0
          endif
        endif

        if (blha_amptype(blha_curproc).eq.0) then ! tree
          blha_amp(1)=blha_amp(1)+res(1)
        else if (blha_amptype(blha_curproc).eq.2) then ! cctree
          ! 1-2
          call BLHA_cctree(1,2,
     &      (blha_CA-2*blha_CF)*res(1)
     &    )
          ! 1-3
          call BLHA_cctree(1,3,
     &      -blha_CA*res(1)
     &    )
          ! 2-3
          call BLHA_cctree(1,3,
     &      -blha_CA*res(1)
     &    )
        else if (blha_amptype(blha_curproc).eq.3) then ! sctree
          ! TODO
          call BLHA_amptypeerror(blha_amptype(blha_curproc),__FILE__,__LINE__)
        else 
          call BLHA_amptypeerror(blha_amptype(blha_curproc),__FILE__,__LINE__)
        endif

        return
      endif


!
! u~ d -> g a e- ve~ 
! 
      Call SUBDGEMVEBMUMMUP(helastestmom1,m2helastest)       !! PDF(5, 4, 3, 2, 1, 0, -1, ......, -5)
                                                      !! for(b, c, s, u, d, g, d_bar, ..., b_bar)
      m2s_qqwzglu_hel=m2s_qqwzglu_hel+m2helastest(1)  !! (switched u <-> d comp. to CTEQ)
     &  *(pdf(-2,1)*pdf(1,2)+pdf(-4,1)*pdf(3,2))


! u~ g -> d~ ga e- ve~ 
!
      Call SUBGDBEMVEBMUMMUP(helastestmom1,m2helastest)

      m2s_qqwzglu_hel=m2s_qqwzglu_hel+m2helastest(1)
     &  *(pdf(-2,1)*pdf(0,2)+pdf(-4,1)*pdf(0,2))
 

! d g -> u ga e- ve~ 

      Call SDGUBEMVEBMUMMUP(helastestmom1,m2helastest)
 
      m2s_qqwzglu_hel=m2s_qqwzglu_hel+m2helastest(1)
     &  *(pdf(1,1)*pdf(0,2)+pdf(3,1)*pdf(0,2))

!
! switch initial state
!
      do k=0,3
         helastestmom1(k,1)=p(k,2)
         helastestmom1(k,2)=p(k,1)
      enddo

!
! d u~ -> g ga e- ve~ 
!
      Call SUBDGEMVEBMUMMUP(helastestmom1,m2helastest)
 
      m2s_qqwzglu_hel=m2s_qqwzglu_hel+m2helastest(1)  !! switched u <-> d in cteq!
     &  *(pdf(-2,2)*pdf(1,1)+pdf(-4,2)*pdf(3,1))

! g u~ -> d~ ga e- ve~ 
!
      Call SUBGDBEMVEBMUMMUP(helastestmom1,m2helastest)

      m2s_qqwzglu_hel=m2s_qqwzglu_hel+m2helastest(1)
     &  *(pdf(-2,2)*pdf(0,1)+pdf(-4,2)*pdf(0,1))
!
! g d -> u ga e- ve~ 
! !
      Call SDGUBEMVEBMUMMUP(helastestmom1,m2helastest)

      m2s_qqwzglu_hel=m2s_qqwzglu_hel+m2helastest(1)
     &  *(pdf(1,2)*pdf(0,1)+pdf(3,2)*pdf(0,1))

      end


********************************************************************************
********************************************************************************

      SUBROUTINE SUBDGEMVEBMUMMUP(P1,ANS)
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : u~ d -> g e- ve~ mu- mu+  
C  
C Crossing   1 is u~ d -> g e- ve~ mu- mu+  
      IMPLICIT NONE
C  
C CONSTANTS
C  
#include "genps.inc"

      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB= 128, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      Double Precision P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      Double Precision T
      Double Precision UBDGEMVEBMUMMU
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      Double Precision hwgt, xtry, xrej
      INTEGER idum, ngood, igood(ncomb), jhel, j
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2(maxamps), jamp2(0:maxamps)
      common/to_ampsUBDGEMVEBMUMMUP/  amp2,       jamp2

      character*79         hel_buff
      common/to_helicityUBDGEMVEBMUMMUP/  hel_buff

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixUBDGEMVEBMUMMUP/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigsUBDGEMVEBMUMMUP/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.true./
      SAVE igood, IDUM, jhel
      DATA NGRAPHS /   24/          
      DATA jamp2(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1,7) /-1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1,7) /-1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1,7) /-1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   4),IHEL=1,7) /-1,-1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   5),IHEL=1,7) /-1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1,7) /-1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1,7) /-1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1,7) /-1,-1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1,7) /-1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  10),IHEL=1,7) /-1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  11),IHEL=1,7) /-1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1,7) /-1,-1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1,7) /-1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1,7) /-1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1,7) /-1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  16),IHEL=1,7) /-1,-1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  17),IHEL=1,7) /-1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1,7) /-1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1,7) /-1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1,7) /-1,-1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1,7) /-1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  22),IHEL=1,7) /-1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  23),IHEL=1,7) /-1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1,7) /-1,-1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1,7) /-1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1,7) /-1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1,7) /-1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  28),IHEL=1,7) /-1,-1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  29),IHEL=1,7) /-1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1,7) /-1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1,7) /-1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1,7) /-1,-1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1,7) /-1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  34),IHEL=1,7) /-1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  35),IHEL=1,7) /-1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1,7) /-1, 1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1,7) /-1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1,7) /-1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1,7) /-1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  40),IHEL=1,7) /-1, 1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  41),IHEL=1,7) /-1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1,7) /-1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1,7) /-1, 1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1,7) /-1, 1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1,7) /-1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  46),IHEL=1,7) /-1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  47),IHEL=1,7) /-1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1,7) /-1, 1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  49),IHEL=1,7) /-1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  50),IHEL=1,7) /-1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  51),IHEL=1,7) /-1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  52),IHEL=1,7) /-1, 1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  53),IHEL=1,7) /-1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  54),IHEL=1,7) /-1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  55),IHEL=1,7) /-1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  56),IHEL=1,7) /-1, 1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  57),IHEL=1,7) /-1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  58),IHEL=1,7) /-1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  59),IHEL=1,7) /-1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  60),IHEL=1,7) /-1, 1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  61),IHEL=1,7) /-1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  62),IHEL=1,7) /-1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  63),IHEL=1,7) /-1, 1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  64),IHEL=1,7) /-1, 1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  65),IHEL=1,7) / 1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  66),IHEL=1,7) / 1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  67),IHEL=1,7) / 1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  68),IHEL=1,7) / 1,-1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  69),IHEL=1,7) / 1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  70),IHEL=1,7) / 1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  71),IHEL=1,7) / 1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  72),IHEL=1,7) / 1,-1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  73),IHEL=1,7) / 1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  74),IHEL=1,7) / 1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  75),IHEL=1,7) / 1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  76),IHEL=1,7) / 1,-1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  77),IHEL=1,7) / 1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  78),IHEL=1,7) / 1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  79),IHEL=1,7) / 1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  80),IHEL=1,7) / 1,-1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  81),IHEL=1,7) / 1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  82),IHEL=1,7) / 1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  83),IHEL=1,7) / 1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  84),IHEL=1,7) / 1,-1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  85),IHEL=1,7) / 1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  86),IHEL=1,7) / 1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  87),IHEL=1,7) / 1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  88),IHEL=1,7) / 1,-1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  89),IHEL=1,7) / 1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  90),IHEL=1,7) / 1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  91),IHEL=1,7) / 1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  92),IHEL=1,7) / 1,-1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  93),IHEL=1,7) / 1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  94),IHEL=1,7) / 1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  95),IHEL=1,7) / 1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  96),IHEL=1,7) / 1,-1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  97),IHEL=1,7) / 1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  98),IHEL=1,7) / 1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  99),IHEL=1,7) / 1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 100),IHEL=1,7) / 1, 1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 101),IHEL=1,7) / 1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 102),IHEL=1,7) / 1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 103),IHEL=1,7) / 1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 104),IHEL=1,7) / 1, 1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 105),IHEL=1,7) / 1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 106),IHEL=1,7) / 1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 107),IHEL=1,7) / 1, 1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 108),IHEL=1,7) / 1, 1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 109),IHEL=1,7) / 1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 110),IHEL=1,7) / 1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 111),IHEL=1,7) / 1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 112),IHEL=1,7) / 1, 1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 113),IHEL=1,7) / 1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 114),IHEL=1,7) / 1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 115),IHEL=1,7) / 1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 116),IHEL=1,7) / 1, 1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 117),IHEL=1,7) / 1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 118),IHEL=1,7) / 1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 119),IHEL=1,7) / 1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 120),IHEL=1,7) / 1, 1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 121),IHEL=1,7) / 1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 122),IHEL=1,7) / 1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 123),IHEL=1,7) / 1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 124),IHEL=1,7) / 1, 1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 125),IHEL=1,7) / 1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 126),IHEL=1,7) / 1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 127),IHEL=1,7) / 1, 1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 128),IHEL=1,7) / 1, 1, 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1,7) / 1, 2, 3, 4, 5, 6, 7/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
!       CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
!       IF (multi_channel) THEN
!           DO IHEL=1,NGRAPHS
!               amp2(ihel)=0d0
!               jamp2(ihel)=0d0
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
                 T=UBDGEMVEBMUMMU(P1 ,NHEL(1,IHEL),JC(1))            
                 ANS(IPROC)=ANS(IPROC)+T
                  IF (T .GT. 0D0 .AND. .NOT. GOODHEL(IHEL,IPROC)) THEN
                      GOODHEL(IHEL,IPROC)=.TRUE.
                      NGOOD = NGOOD +1
                      IGOOD(NGOOD) = IHEL
C                WRITE(*,*) ngood,IHEL,T
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
              T=UBDGEMVEBMUMMU(P1 ,NHEL(1,IHEL),JC(1))            
           ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
!       IF (MULTI_CHANNEL) THEN
!           XTOT=0D0
!           DO IHEL=1,MAPCONFIG(0)
!               XTOT=XTOT+AMP2(MAPCONFIG(IHEL))
!           ENDDO
!           ANS(IPROC)=ANS(IPROC)*AMP2(MAPCONFIG(ICONFIG))/XTOT
!       ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      Double Precision FUNCTION UBDGEMVEBMUMMU(P,NHEL,IC)
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : u~ d -> g e- ve~ mu- mu+  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=  24,NEIGEN=  1) 
#include "genps.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  38, NCOLOR=   1) 
      Double Precision     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      Double Precision P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      Double Complex ZTEMP
      Double Precision DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      Double Complex AMP(NGRAPHS), JAMP(NCOLOR)
      Double Complex W(6,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2(maxamps), jamp2(0:maxamps)
      common/to_ampsUBDGEMVEBMUMMUP/  amp2,       jamp2
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/global.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     4/                                  
C               T[1,2,3]                                                   
C ----------
C BEGIN CODE
C ----------
      CALL OXXXXX(P(0,1   ),ZERO ,NHEL(1   ),-1*IC(1   ),W(1,1   ))        
      CALL IXXXXX(P(0,2   ),ZERO ,NHEL(2   ),+1*IC(2   ),W(1,2   ))        
      CALL VXXXXX(P(0,3   ),ZERO ,NHEL(3   ),+1*IC(3   ),W(1,3   ))        
      CALL OXXXXX(P(0,4   ),ZERO ,NHEL(4   ),+1*IC(4   ),W(1,4   ))        
      CALL IXXXXX(P(0,5   ),ZERO ,NHEL(5   ),-1*IC(5   ),W(1,5   ))        
      CALL OXXXXX(P(0,6   ),ZERO ,NHEL(6   ),+1*IC(6   ),W(1,6   ))        
      CALL IXXXXX(P(0,7   ),ZERO ,NHEL(7   ),-1*IC(7   ),W(1,7   ))        
      CALL FVIXXX(W(1,2   ),W(1,3   ),GG ,ZERO    ,ZERO    ,W(1,8   ))     
      CALL JIOXXX(W(1,5   ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,9   ))    
      CALL JIOXXX(W(1,7   ),W(1,6   ),GAL ,ZERO    ,ZERO    ,W(1,10  ))    
      CALL FVIXXX(W(1,8   ),W(1,9   ),GWF ,ZERO    ,ZERO    ,W(1,11  ))    
      CALL IOVXXX(W(1,11  ),W(1,1   ),W(1,10  ),GAU ,AMP(1   ))            
      CALL JIOXXX(W(1,7   ),W(1,6   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,12  ))    
      CALL IOVXXX(W(1,11  ),W(1,1   ),W(1,12  ),GZU ,AMP(2   ))            
      CALL FVOXXX(W(1,1   ),W(1,9   ),GWF ,ZERO    ,ZERO    ,W(1,13  ))    
      CALL FVIXXX(W(1,2   ),W(1,10  ),GAD ,ZERO    ,ZERO    ,W(1,14  ))    
      CALL IOVXXX(W(1,14  ),W(1,13  ),W(1,3   ),GG ,AMP(3   ))             
      CALL FVIXXX(W(1,2   ),W(1,12  ),GZD ,ZERO    ,ZERO    ,W(1,15  ))    
      CALL IOVXXX(W(1,15  ),W(1,13  ),W(1,3   ),GG ,AMP(4   ))             
      CALL JIOXXX(W(1,8   ),W(1,1   ),GWF ,WMASS   ,WWIDTH  ,W(1,16  ))    
      CALL FVIXXX(W(1,5   ),W(1,16  ),GWF ,ZERO    ,ZERO    ,W(1,17  ))    
      CALL JIOXXX(W(1,17  ),W(1,4   ),GAL ,ZERO    ,ZERO    ,W(1,18  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,18  ),GAL ,AMP(5   ))            
      CALL JIOXXX(W(1,17  ),W(1,4   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,19  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,19  ),GZL ,AMP(6   ))            
      CALL FVOXXX(W(1,4   ),W(1,16  ),GWF ,ZERO    ,ZERO    ,W(1,20  ))    
      CALL JIOXXX(W(1,5   ),W(1,20  ),GZN ,ZMASS   ,ZWIDTH  ,W(1,21  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,21  ),GZL ,AMP(7   ))            
      CALL IOVXXX(W(1,8   ),W(1,13  ),W(1,10  ),GAD ,AMP(8   ))            
      CALL IOVXXX(W(1,8   ),W(1,13  ),W(1,12  ),GZD ,AMP(9   ))  
!      
      if (.not.with_anom) then      
      CALL JVVXXX(W(1,9   ),W(1,16  ),GWWA ,ZERO    ,ZERO    ,W(1,         
     &     22  ))
      elseif (with_anom) then
      CALL JVVXXX_anomwaj0(W(1,9   ),W(1,16  ),GWWA ,ZERO    ,ZERO    ,W(1,         
     &     22  ))
      endif
!
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,22  ),GAL ,AMP(10  ))
!
      if (.not.with_anom) then             
      CALL JVVXXX(W(1,9   ),W(1,16  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     23  )) 
      elseif (with_anom) then
      CALL JVVXXX_anomwzj(W(1,9   ),W(1,16  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     23  ))
      endif
!      print*, with_anom
!                                                         
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,23  ),GZL ,AMP(11  ))            
      CALL FVOXXX(W(1,6   ),W(1,16  ),GWF ,ZERO    ,ZERO    ,W(1,24  ))    
      CALL IOVXXX(W(1,7   ),W(1,24  ),W(1,9   ),GWF ,AMP(12  ))            
      CALL FVIXXX(W(1,2   ),W(1,9   ),GWF ,ZERO    ,ZERO    ,W(1,25  ))    
      CALL FVOXXX(W(1,1   ),W(1,10  ),GAU ,ZERO    ,ZERO    ,W(1,26  ))    
      CALL IOVXXX(W(1,25  ),W(1,26  ),W(1,3   ),GG ,AMP(13  ))             
      CALL FVOXXX(W(1,1   ),W(1,12  ),GZU ,ZERO    ,ZERO    ,W(1,27  ))    
      CALL IOVXXX(W(1,25  ),W(1,27  ),W(1,3   ),GG ,AMP(14  ))             
      CALL FVOXXX(W(1,1   ),W(1,3   ),GG ,ZERO    ,ZERO    ,W(1,28  ))     
      CALL FVOXXX(W(1,28  ),W(1,9   ),GWF ,ZERO    ,ZERO    ,W(1,29  ))    
      CALL IOVXXX(W(1,2   ),W(1,29  ),W(1,10  ),GAD ,AMP(15  ))            
      CALL IOVXXX(W(1,2   ),W(1,29  ),W(1,12  ),GZD ,AMP(16  ))            
      CALL JIOXXX(W(1,2   ),W(1,28  ),GWF ,WMASS   ,WWIDTH  ,W(1,30  ))    
      CALL FVIXXX(W(1,5   ),W(1,30  ),GWF ,ZERO    ,ZERO    ,W(1,31  ))    
      CALL JIOXXX(W(1,31  ),W(1,4   ),GAL ,ZERO    ,ZERO    ,W(1,32  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,32  ),GAL ,AMP(17  ))            
      CALL JIOXXX(W(1,31  ),W(1,4   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,33  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,33  ),GZL ,AMP(18  ))            
      CALL FVOXXX(W(1,4   ),W(1,30  ),GWF ,ZERO    ,ZERO    ,W(1,34  ))    
      CALL JIOXXX(W(1,5   ),W(1,34  ),GZN ,ZMASS   ,ZWIDTH  ,W(1,35  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,35  ),GZL ,AMP(19  ))            
      CALL IOVXXX(W(1,25  ),W(1,28  ),W(1,10  ),GAU ,AMP(20  ))            
      CALL IOVXXX(W(1,25  ),W(1,28  ),W(1,12  ),GZU ,AMP(21  ))
!
      if (.not.with_anom) then                  
      CALL JVVXXX(W(1,9   ),W(1,30  ),GWWA ,ZERO    ,ZERO    ,W(1,         
     &     36  ))                                                          
      elseif (with_anom) then
      CALL JVVXXX_anomwaj0(W(1,9   ),W(1,30  ),GWWA ,ZERO    ,ZERO    ,W(1,         
     &     36  ))                                                          
      endif
!
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,36  ),GAL ,AMP(22  ))
      if (.not.with_anom) then                              
      CALL JVVXXX(W(1,9   ),W(1,30  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     37  ))
      elseif (with_anom) then
      CALL JVVXXX_anomwzj(W(1,9   ),W(1,30  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     37  ))
      endif
!
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,37  ),GZL ,AMP(23  ))            
      CALL FVOXXX(W(1,6   ),W(1,30  ),GWF ,ZERO    ,ZERO    ,W(1,38  ))    
      CALL IOVXXX(W(1,7   ),W(1,38  ),W(1,9   ),GWF ,AMP(24  ))            
      JAMP(   1) = -AMP(   1)-AMP(   2)-AMP(   3)-AMP(   4)-AMP(   5)
     &             -AMP(   6)-AMP(   7)-AMP(   8)-AMP(   9)-AMP(  10)
     &             -AMP(  11)-AMP(  12)-AMP(  13)-AMP(  14)-AMP(  15)
     &             -AMP(  16)-AMP(  17)-AMP(  18)-AMP(  19)-AMP(  20)
     &             -AMP(  21)-AMP(  22)-AMP(  23)-AMP(  24)
      UBDGEMVEBMUMMU = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          UBDGEMVEBMUMMU =UBDGEMVEBMUMMU+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2(i)=amp2(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2(i)=Jamp2(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END


      Double complex FUNCTION AMPUBDGEMVEBMUMMU(P,NHEL)
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : u~ d -> g e- ve~ mu- mu+  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=  24,NEIGEN=  1) 
#include "genps.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  38, NCOLOR=   1) 
      Double Precision     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      Double Precision P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I

      Double Precision DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      Double Complex AMP(NGRAPHS)
      Double Complex W(6,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2(maxamps), jamp2(0:maxamps)
      common/to_ampsUBDGEMVEBMUMMUP/  amp2,       jamp2
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/global.inc"

C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     4/                                  
C               T[1,2,3]                                                   
C ----------
C BEGIN CODE
C ----------
      CALL OXXXXX(P(0,1   ),ZERO ,NHEL(1   ),-1,W(1,1   ))        
      CALL IXXXXX(P(0,2   ),ZERO ,NHEL(2   ),+1,W(1,2   ))        
      CALL VXXXXX(P(0,3   ),ZERO ,NHEL(3   ),+1,W(1,3   ))        
      CALL OXXXXX(P(0,4   ),ZERO ,NHEL(4   ),+1,W(1,4   ))        
      CALL IXXXXX(P(0,5   ),ZERO ,NHEL(5   ),-1,W(1,5   ))        
      CALL OXXXXX(P(0,6   ),ZERO ,NHEL(6   ),+1,W(1,6   ))        
      CALL IXXXXX(P(0,7   ),ZERO ,NHEL(7   ),-1,W(1,7   ))        
      CALL FVIXXX(W(1,2   ),W(1,3   ),GG ,ZERO    ,ZERO    ,W(1,8   ))     
      CALL JIOXXX(W(1,5   ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,9   ))    
      CALL JIOXXX(W(1,7   ),W(1,6   ),GAL ,ZERO    ,ZERO    ,W(1,10  ))    
      CALL FVIXXX(W(1,8   ),W(1,9   ),GWF ,ZERO    ,ZERO    ,W(1,11  ))    
      CALL IOVXXX(W(1,11  ),W(1,1   ),W(1,10  ),GAU ,AMP(1   ))            
      CALL JIOXXX(W(1,7   ),W(1,6   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,12  ))    
      CALL IOVXXX(W(1,11  ),W(1,1   ),W(1,12  ),GZU ,AMP(2   ))            
      CALL FVOXXX(W(1,1   ),W(1,9   ),GWF ,ZERO    ,ZERO    ,W(1,13  ))    
      CALL FVIXXX(W(1,2   ),W(1,10  ),GAD ,ZERO    ,ZERO    ,W(1,14  ))    
      CALL IOVXXX(W(1,14  ),W(1,13  ),W(1,3   ),GG ,AMP(3   ))             
      CALL FVIXXX(W(1,2   ),W(1,12  ),GZD ,ZERO    ,ZERO    ,W(1,15  ))    
      CALL IOVXXX(W(1,15  ),W(1,13  ),W(1,3   ),GG ,AMP(4   ))             
      CALL JIOXXX(W(1,8   ),W(1,1   ),GWF ,WMASS   ,WWIDTH  ,W(1,16  ))    
      CALL FVIXXX(W(1,5   ),W(1,16  ),GWF ,ZERO    ,ZERO    ,W(1,17  ))    
      CALL JIOXXX(W(1,17  ),W(1,4   ),GAL ,ZERO    ,ZERO    ,W(1,18  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,18  ),GAL ,AMP(5   ))            
      CALL JIOXXX(W(1,17  ),W(1,4   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,19  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,19  ),GZL ,AMP(6   ))            
      CALL FVOXXX(W(1,4   ),W(1,16  ),GWF ,ZERO    ,ZERO    ,W(1,20  ))    
      CALL JIOXXX(W(1,5   ),W(1,20  ),GZN ,ZMASS   ,ZWIDTH  ,W(1,21  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,21  ),GZL ,AMP(7   ))            
      CALL IOVXXX(W(1,8   ),W(1,13  ),W(1,10  ),GAD ,AMP(8   ))            
      CALL IOVXXX(W(1,8   ),W(1,13  ),W(1,12  ),GZD ,AMP(9   ))  
!      
      if (.not.with_anom) then      
      CALL JVVXXX(W(1,9   ),W(1,16  ),GWWA ,ZERO    ,ZERO    ,W(1,         
     &     22  ))
      elseif (with_anom) then
      CALL JVVXXX_anomwaj0(W(1,9   ),W(1,16  ),GWWA ,ZERO    ,ZERO    ,W(1,         
     &     22  ))
      endif
!
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,22  ),GAL ,AMP(10  ))
!
      if (.not.with_anom) then             
      CALL JVVXXX(W(1,9   ),W(1,16  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     23  )) 
      elseif (with_anom) then
      CALL JVVXXX_anomwzj(W(1,9   ),W(1,16  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     23  ))
      endif
!                                                         
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,23  ),GZL ,AMP(11  ))            
      CALL FVOXXX(W(1,6   ),W(1,16  ),GWF ,ZERO    ,ZERO    ,W(1,24  ))    
      CALL IOVXXX(W(1,7   ),W(1,24  ),W(1,9   ),GWF ,AMP(12  ))            
      CALL FVIXXX(W(1,2   ),W(1,9   ),GWF ,ZERO    ,ZERO    ,W(1,25  ))    
      CALL FVOXXX(W(1,1   ),W(1,10  ),GAU ,ZERO    ,ZERO    ,W(1,26  ))    
      CALL IOVXXX(W(1,25  ),W(1,26  ),W(1,3   ),GG ,AMP(13  ))             
      CALL FVOXXX(W(1,1   ),W(1,12  ),GZU ,ZERO    ,ZERO    ,W(1,27  ))    
      CALL IOVXXX(W(1,25  ),W(1,27  ),W(1,3   ),GG ,AMP(14  ))             
      CALL FVOXXX(W(1,1   ),W(1,3   ),GG ,ZERO    ,ZERO    ,W(1,28  ))     
      CALL FVOXXX(W(1,28  ),W(1,9   ),GWF ,ZERO    ,ZERO    ,W(1,29  ))    
      CALL IOVXXX(W(1,2   ),W(1,29  ),W(1,10  ),GAD ,AMP(15  ))            
      CALL IOVXXX(W(1,2   ),W(1,29  ),W(1,12  ),GZD ,AMP(16  ))            
      CALL JIOXXX(W(1,2   ),W(1,28  ),GWF ,WMASS   ,WWIDTH  ,W(1,30  ))    
      CALL FVIXXX(W(1,5   ),W(1,30  ),GWF ,ZERO    ,ZERO    ,W(1,31  ))    
      CALL JIOXXX(W(1,31  ),W(1,4   ),GAL ,ZERO    ,ZERO    ,W(1,32  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,32  ),GAL ,AMP(17  ))            
      CALL JIOXXX(W(1,31  ),W(1,4   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,33  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,33  ),GZL ,AMP(18  ))            
      CALL FVOXXX(W(1,4   ),W(1,30  ),GWF ,ZERO    ,ZERO    ,W(1,34  ))    
      CALL JIOXXX(W(1,5   ),W(1,34  ),GZN ,ZMASS   ,ZWIDTH  ,W(1,35  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,35  ),GZL ,AMP(19  ))            
      CALL IOVXXX(W(1,25  ),W(1,28  ),W(1,10  ),GAU ,AMP(20  ))            
      CALL IOVXXX(W(1,25  ),W(1,28  ),W(1,12  ),GZU ,AMP(21  ))
!
      if (.not.with_anom) then                  
      CALL JVVXXX(W(1,9   ),W(1,30  ),GWWA ,ZERO    ,ZERO    ,W(1,         
     &     36  ))                                                          
      elseif (with_anom) then
      CALL JVVXXX_anomwaj0(W(1,9   ),W(1,30  ),GWWA ,ZERO    ,ZERO    ,W(1,         
     &     36  ))                                                          
      endif
!
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,36  ),GAL ,AMP(22  ))
      if (.not.with_anom) then                              
      CALL JVVXXX(W(1,9   ),W(1,30  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     37  ))
      elseif (with_anom) then
      CALL JVVXXX_anomwzj(W(1,9   ),W(1,30  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     37  ))
      endif
!
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,37  ),GZL ,AMP(23  ))            
      CALL FVOXXX(W(1,6   ),W(1,30  ),GWF ,ZERO    ,ZERO    ,W(1,38  ))    
      CALL IOVXXX(W(1,7   ),W(1,38  ),W(1,9   ),GWF ,AMP(24  ))

      AMPUBDGEMVEBMUMMU =-1d0*( -AMP(   1)-AMP(   2)-AMP(   3)-AMP(   4)-AMP(   5)
     &             -AMP(   6)-AMP(   7)-AMP(   8)-AMP(   9)-AMP(  10)
     &             -AMP(  11)-AMP(  12)-AMP(  13)-AMP(  14)-AMP(  15)
     &             -AMP(  16)-AMP(  17)-AMP(  18)-AMP(  19)-AMP(  20)
     &             -AMP(  21)-AMP(  22)-AMP(  23)-AMP(  24))
      END


        subroutine UBDGEMVEBMUMMUDP(P,UDGAENU)
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : u~ d -> g e- ve~ mu- mu+  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=  24,NEIGEN=  1) 
#include "genps.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  38, NCOLOR=   1) 
      Double Precision     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      Double Precision P(0:3,NEXTERNAL)
      INTEGER mu,nu
C  
C LOCAL VARIABLES 
C


      Double Complex AMP(NGRAPHS)
      Double Complex W(6,NWAVEFUNCS),UDGAENU(0:3,-1:1)
C  
C GLOBAL VARIABLES
C  
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/global.inc"

C ----------
C BEGIN CODE
C ----------
      do mu=0,3
      do nu=-1,1,2
      CALL OXXXXX(P(0,1   ),ZERO ,1,-1,W(1,1   ))        
      CALL IXXXXX(P(0,2   ),ZERO ,-1,+1,W(1,2   ))        
      CALL VCARTXDP(P(0,3   ),mu,+1,W(1,3   ))       
      CALL OXXXXX(P(0,4   ),ZERO ,-1,+1,W(1,4   ))        
      CALL IXXXXX(P(0,5   ),ZERO ,1,-1,W(1,5   ))
      if (nu.eq.-1) then        
      CALL OXXXXX(P(0,6   ),ZERO ,-1,+1,W(1,6   ))        
      CALL IXXXXX(P(0,7   ),ZERO , 1,-1,W(1,7   ))
      elseif (nu.eq.1) then
      CALL OXXXXX(P(0,6   ),ZERO ,1,+1,W(1,6   ))        
      CALL IXXXXX(P(0,7   ),ZERO ,-1,-1,W(1,7   ))
      endif
        
      CALL FVIXXX(W(1,2   ),W(1,3   ),GG ,ZERO    ,ZERO    ,W(1,8   ))     
      CALL JIOXXX(W(1,5   ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,9   ))    
      CALL JIOXXX(W(1,7   ),W(1,6   ),GAL ,ZERO    ,ZERO    ,W(1,10  ))    
      CALL FVIXXX(W(1,8   ),W(1,9   ),GWF ,ZERO    ,ZERO    ,W(1,11  ))    
      CALL IOVXXX(W(1,11  ),W(1,1   ),W(1,10  ),GAU ,AMP(1   ))            
      CALL JIOXXX(W(1,7   ),W(1,6   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,12  ))    
      CALL IOVXXX(W(1,11  ),W(1,1   ),W(1,12  ),GZU ,AMP(2   ))            
      CALL FVOXXX(W(1,1   ),W(1,9   ),GWF ,ZERO    ,ZERO    ,W(1,13  ))    
      CALL FVIXXX(W(1,2   ),W(1,10  ),GAD ,ZERO    ,ZERO    ,W(1,14  ))    
      CALL IOVXXX(W(1,14  ),W(1,13  ),W(1,3   ),GG ,AMP(3   ))             
      CALL FVIXXX(W(1,2   ),W(1,12  ),GZD ,ZERO    ,ZERO    ,W(1,15  ))    
      CALL IOVXXX(W(1,15  ),W(1,13  ),W(1,3   ),GG ,AMP(4   ))             
      CALL JIOXXX(W(1,8   ),W(1,1   ),GWF ,WMASS   ,WWIDTH  ,W(1,16  ))    
      CALL FVIXXX(W(1,5   ),W(1,16  ),GWF ,ZERO    ,ZERO    ,W(1,17  ))    
      CALL JIOXXX(W(1,17  ),W(1,4   ),GAL ,ZERO    ,ZERO    ,W(1,18  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,18  ),GAL ,AMP(5   ))            
      CALL JIOXXX(W(1,17  ),W(1,4   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,19  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,19  ),GZL ,AMP(6   ))            
      CALL FVOXXX(W(1,4   ),W(1,16  ),GWF ,ZERO    ,ZERO    ,W(1,20  ))    
      CALL JIOXXX(W(1,5   ),W(1,20  ),GZN ,ZMASS   ,ZWIDTH  ,W(1,21  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,21  ),GZL ,AMP(7   ))            
      CALL IOVXXX(W(1,8   ),W(1,13  ),W(1,10  ),GAD ,AMP(8   ))            
      CALL IOVXXX(W(1,8   ),W(1,13  ),W(1,12  ),GZD ,AMP(9   ))  
!      
      if (.not.with_anom) then      
      CALL JVVXXX(W(1,9   ),W(1,16  ),GWWA ,ZERO    ,ZERO    ,W(1,         
     &     22  ))
      elseif (with_anom) then
      CALL JVVXXX_anomwaj0(W(1,9   ),W(1,16  ),GWWA ,ZERO    ,ZERO    ,W(1,         
     &     22  ))
      endif
!
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,22  ),GAL ,AMP(10  ))
!
      if (.not.with_anom) then             
      CALL JVVXXX(W(1,9   ),W(1,16  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     23  )) 
      elseif (with_anom) then
      CALL JVVXXX_anomwzj(W(1,9   ),W(1,16  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     23  ))
      endif
!                                                         
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,23  ),GZL ,AMP(11  ))            
      CALL FVOXXX(W(1,6   ),W(1,16  ),GWF ,ZERO    ,ZERO    ,W(1,24  ))    
      CALL IOVXXX(W(1,7   ),W(1,24  ),W(1,9   ),GWF ,AMP(12  ))            
      CALL FVIXXX(W(1,2   ),W(1,9   ),GWF ,ZERO    ,ZERO    ,W(1,25  ))    
      CALL FVOXXX(W(1,1   ),W(1,10  ),GAU ,ZERO    ,ZERO    ,W(1,26  ))    
      CALL IOVXXX(W(1,25  ),W(1,26  ),W(1,3   ),GG ,AMP(13  ))             
      CALL FVOXXX(W(1,1   ),W(1,12  ),GZU ,ZERO    ,ZERO    ,W(1,27  ))    
      CALL IOVXXX(W(1,25  ),W(1,27  ),W(1,3   ),GG ,AMP(14  ))             
      CALL FVOXXX(W(1,1   ),W(1,3   ),GG ,ZERO    ,ZERO    ,W(1,28  ))     
      CALL FVOXXX(W(1,28  ),W(1,9   ),GWF ,ZERO    ,ZERO    ,W(1,29  ))    
      CALL IOVXXX(W(1,2   ),W(1,29  ),W(1,10  ),GAD ,AMP(15  ))            
      CALL IOVXXX(W(1,2   ),W(1,29  ),W(1,12  ),GZD ,AMP(16  ))            
      CALL JIOXXX(W(1,2   ),W(1,28  ),GWF ,WMASS   ,WWIDTH  ,W(1,30  ))    
      CALL FVIXXX(W(1,5   ),W(1,30  ),GWF ,ZERO    ,ZERO    ,W(1,31  ))    
      CALL JIOXXX(W(1,31  ),W(1,4   ),GAL ,ZERO    ,ZERO    ,W(1,32  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,32  ),GAL ,AMP(17  ))            
      CALL JIOXXX(W(1,31  ),W(1,4   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,33  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,33  ),GZL ,AMP(18  ))            
      CALL FVOXXX(W(1,4   ),W(1,30  ),GWF ,ZERO    ,ZERO    ,W(1,34  ))    
      CALL JIOXXX(W(1,5   ),W(1,34  ),GZN ,ZMASS   ,ZWIDTH  ,W(1,35  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,35  ),GZL ,AMP(19  ))            
      CALL IOVXXX(W(1,25  ),W(1,28  ),W(1,10  ),GAU ,AMP(20  ))            
      CALL IOVXXX(W(1,25  ),W(1,28  ),W(1,12  ),GZU ,AMP(21  ))
!
      if (.not.with_anom) then                  
      CALL JVVXXX(W(1,9   ),W(1,30  ),GWWA ,ZERO    ,ZERO    ,W(1,         
     &     36  ))                                                          
      elseif (with_anom) then
      CALL JVVXXX_anomwaj0(W(1,9   ),W(1,30  ),GWWA ,ZERO    ,ZERO    ,W(1,         
     &     36  ))                                                          
      endif
!
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,36  ),GAL ,AMP(22  ))
      if (.not.with_anom) then                              
      CALL JVVXXX(W(1,9   ),W(1,30  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     37  ))
      elseif (with_anom) then
      CALL JVVXXX_anomwzj(W(1,9   ),W(1,30  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     37  ))
      endif
!
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,37  ),GZL ,AMP(23  ))            
      CALL FVOXXX(W(1,6   ),W(1,30  ),GWF ,ZERO    ,ZERO    ,W(1,38  ))    
      CALL IOVXXX(W(1,7   ),W(1,38  ),W(1,9   ),GWF ,AMP(24  ))

      UDGAENU(mu,nu) = -AMP(   1)-AMP(   2)-AMP(   3)-AMP(   4)-AMP(   5)
     &             -AMP(   6)-AMP(   7)-AMP(   8)-AMP(   9)-AMP(  10)
     &             -AMP(  11)-AMP(  12)-AMP(  13)-AMP(  14)-AMP(  15)
     &             -AMP(  16)-AMP(  17)-AMP(  18)-AMP(  19)-AMP(  20)
     &             -AMP(  21)-AMP(  22)-AMP(  23)-AMP(  24)
      enddo
      enddo
      END


      SUBROUTINE SUBGDBEMVEBMUMMUP(P1,ANS)
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : u~ g -> d~ e- ve~ mu- mu+  
C  
C Crossing   1 is u~ g -> d~ e- ve~ mu- mu+  
      IMPLICIT NONE
C  
C CONSTANTS
C  
#include "genps.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB= 128, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      Double Precision P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      Double Precision T
      Double Precision UBGDBEMVEBMUMM
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      Double Precision hwgt, xtry, xrej
      INTEGER idum, ngood, igood(ncomb), jhel, j
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2(maxamps), jamp2(0:maxamps)
      common/to_ampsUBGDBEMVEBMUMMUP/  amp2,       jamp2

      character*79         hel_buff
      common/to_helicityUBGDBEMVEBMUMMUP/  hel_buff

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixUBGDBEMVEBMUMMUP/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigsUBGDBEMVEBMUMMUP/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.true./
      SAVE igood, IDUM, jhel
      DATA NGRAPHS /   24/          
      DATA jamp2(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1,7) /-1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1,7) /-1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1,7) /-1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   4),IHEL=1,7) /-1,-1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   5),IHEL=1,7) /-1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1,7) /-1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1,7) /-1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1,7) /-1,-1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1,7) /-1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  10),IHEL=1,7) /-1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  11),IHEL=1,7) /-1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1,7) /-1,-1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1,7) /-1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1,7) /-1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1,7) /-1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  16),IHEL=1,7) /-1,-1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  17),IHEL=1,7) /-1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1,7) /-1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1,7) /-1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1,7) /-1,-1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1,7) /-1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  22),IHEL=1,7) /-1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  23),IHEL=1,7) /-1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1,7) /-1,-1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1,7) /-1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1,7) /-1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1,7) /-1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  28),IHEL=1,7) /-1,-1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  29),IHEL=1,7) /-1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1,7) /-1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1,7) /-1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1,7) /-1,-1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1,7) /-1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  34),IHEL=1,7) /-1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  35),IHEL=1,7) /-1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1,7) /-1, 1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1,7) /-1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1,7) /-1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1,7) /-1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  40),IHEL=1,7) /-1, 1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  41),IHEL=1,7) /-1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1,7) /-1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1,7) /-1, 1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1,7) /-1, 1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1,7) /-1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  46),IHEL=1,7) /-1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  47),IHEL=1,7) /-1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1,7) /-1, 1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  49),IHEL=1,7) /-1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  50),IHEL=1,7) /-1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  51),IHEL=1,7) /-1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  52),IHEL=1,7) /-1, 1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  53),IHEL=1,7) /-1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  54),IHEL=1,7) /-1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  55),IHEL=1,7) /-1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  56),IHEL=1,7) /-1, 1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  57),IHEL=1,7) /-1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  58),IHEL=1,7) /-1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  59),IHEL=1,7) /-1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  60),IHEL=1,7) /-1, 1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  61),IHEL=1,7) /-1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  62),IHEL=1,7) /-1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  63),IHEL=1,7) /-1, 1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  64),IHEL=1,7) /-1, 1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  65),IHEL=1,7) / 1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  66),IHEL=1,7) / 1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  67),IHEL=1,7) / 1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  68),IHEL=1,7) / 1,-1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  69),IHEL=1,7) / 1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  70),IHEL=1,7) / 1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  71),IHEL=1,7) / 1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  72),IHEL=1,7) / 1,-1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  73),IHEL=1,7) / 1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  74),IHEL=1,7) / 1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  75),IHEL=1,7) / 1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  76),IHEL=1,7) / 1,-1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  77),IHEL=1,7) / 1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  78),IHEL=1,7) / 1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  79),IHEL=1,7) / 1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  80),IHEL=1,7) / 1,-1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  81),IHEL=1,7) / 1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  82),IHEL=1,7) / 1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  83),IHEL=1,7) / 1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  84),IHEL=1,7) / 1,-1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  85),IHEL=1,7) / 1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  86),IHEL=1,7) / 1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  87),IHEL=1,7) / 1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  88),IHEL=1,7) / 1,-1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  89),IHEL=1,7) / 1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  90),IHEL=1,7) / 1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  91),IHEL=1,7) / 1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  92),IHEL=1,7) / 1,-1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  93),IHEL=1,7) / 1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  94),IHEL=1,7) / 1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  95),IHEL=1,7) / 1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  96),IHEL=1,7) / 1,-1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  97),IHEL=1,7) / 1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  98),IHEL=1,7) / 1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  99),IHEL=1,7) / 1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 100),IHEL=1,7) / 1, 1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 101),IHEL=1,7) / 1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 102),IHEL=1,7) / 1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 103),IHEL=1,7) / 1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 104),IHEL=1,7) / 1, 1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 105),IHEL=1,7) / 1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 106),IHEL=1,7) / 1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 107),IHEL=1,7) / 1, 1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 108),IHEL=1,7) / 1, 1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 109),IHEL=1,7) / 1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 110),IHEL=1,7) / 1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 111),IHEL=1,7) / 1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 112),IHEL=1,7) / 1, 1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 113),IHEL=1,7) / 1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 114),IHEL=1,7) / 1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 115),IHEL=1,7) / 1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 116),IHEL=1,7) / 1, 1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 117),IHEL=1,7) / 1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 118),IHEL=1,7) / 1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 119),IHEL=1,7) / 1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 120),IHEL=1,7) / 1, 1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 121),IHEL=1,7) / 1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 122),IHEL=1,7) / 1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 123),IHEL=1,7) / 1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 124),IHEL=1,7) / 1, 1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 125),IHEL=1,7) / 1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 126),IHEL=1,7) / 1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 127),IHEL=1,7) / 1, 1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 128),IHEL=1,7) / 1, 1, 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1,7) / 1, 2, 3, 4, 5, 6, 7/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  96/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
!       CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
!       IF (multi_channel) THEN
!           DO IHEL=1,NGRAPHS
!               amp2(ihel)=0d0
!               jamp2(ihel)=0d0
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
                 T=UBGDBEMVEBMUMM(P1 ,NHEL(1,IHEL),JC(1))            
                 ANS(IPROC)=ANS(IPROC)+T
                  IF (T .GT. 0D0 .AND. .NOT. GOODHEL(IHEL,IPROC)) THEN
                      GOODHEL(IHEL,IPROC)=.TRUE.
                      NGOOD = NGOOD +1
                      IGOOD(NGOOD) = IHEL
C                WRITE(*,*) ngood,IHEL,T
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
              T=UBGDBEMVEBMUMM(P1 ,NHEL(1,IHEL),JC(1))            
           ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
!       IF (MULTI_CHANNEL) THEN
!           XTOT=0D0
!           DO IHEL=1,MAPCONFIG(0)
!               XTOT=XTOT+AMP2(MAPCONFIG(IHEL))
!           ENDDO
!           ANS(IPROC)=ANS(IPROC)*AMP2(MAPCONFIG(ICONFIG))/XTOT
!       ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      Double Precision FUNCTION UBGDBEMVEBMUMM(P,NHEL,IC)
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : u~ g -> d~ e- ve~ mu- mu+  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=  24,NEIGEN=  1) 
#include "genps.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  38, NCOLOR=   1) 
      Double Precision     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      Double Precision P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      Double Complex ZTEMP
      Double Precision DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      Double Complex AMP(NGRAPHS), JAMP(NCOLOR)
      Double Complex W(6,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2(maxamps), jamp2(0:maxamps)
      common/to_ampsUBGDBEMVEBMUMMUP/  amp2,       jamp2
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/global.inc"

C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     4/                                  
C               T[1,3,2]                                                   
C ----------
C BEGIN CODE
C ----------
      CALL OXXXXX(P(0,1   ),ZERO ,NHEL(1   ),-1*IC(1   ),W(1,1   ))        
      CALL VXXXXX(P(0,2   ),ZERO ,NHEL(2   ),-1*IC(2   ),W(1,2   ))        
      CALL IXXXXX(P(0,3   ),ZERO ,NHEL(3   ),-1*IC(3   ),W(1,3   ))        
      CALL OXXXXX(P(0,4   ),ZERO ,NHEL(4   ),+1*IC(4   ),W(1,4   ))        
      CALL IXXXXX(P(0,5   ),ZERO ,NHEL(5   ),-1*IC(5   ),W(1,5   ))        
      CALL OXXXXX(P(0,6   ),ZERO ,NHEL(6   ),+1*IC(6   ),W(1,6   ))        
      CALL IXXXXX(P(0,7   ),ZERO ,NHEL(7   ),-1*IC(7   ),W(1,7   ))        
      CALL FVIXXX(W(1,3   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,8   ))     
      CALL JIOXXX(W(1,5   ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,9   ))    
      CALL JIOXXX(W(1,7   ),W(1,6   ),GAL ,ZERO    ,ZERO    ,W(1,10  ))    
      CALL FVIXXX(W(1,8   ),W(1,9   ),GWF ,ZERO    ,ZERO    ,W(1,11  ))    
      CALL IOVXXX(W(1,11  ),W(1,1   ),W(1,10  ),GAU ,AMP(1   ))            
      CALL JIOXXX(W(1,7   ),W(1,6   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,12  ))    
      CALL IOVXXX(W(1,11  ),W(1,1   ),W(1,12  ),GZU ,AMP(2   ))            
      CALL FVOXXX(W(1,1   ),W(1,9   ),GWF ,ZERO    ,ZERO    ,W(1,13  ))    
      CALL FVIXXX(W(1,3   ),W(1,10  ),GAD ,ZERO    ,ZERO    ,W(1,14  ))    
      CALL IOVXXX(W(1,14  ),W(1,13  ),W(1,2   ),GG ,AMP(3   ))             
      CALL FVIXXX(W(1,3   ),W(1,12  ),GZD ,ZERO    ,ZERO    ,W(1,15  ))    
      CALL IOVXXX(W(1,15  ),W(1,13  ),W(1,2   ),GG ,AMP(4   ))             
      CALL JIOXXX(W(1,8   ),W(1,1   ),GWF ,WMASS   ,WWIDTH  ,W(1,16  ))    
      CALL FVIXXX(W(1,5   ),W(1,16  ),GWF ,ZERO    ,ZERO    ,W(1,17  ))    
      CALL JIOXXX(W(1,17  ),W(1,4   ),GAL ,ZERO    ,ZERO    ,W(1,18  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,18  ),GAL ,AMP(5   ))            
      CALL JIOXXX(W(1,17  ),W(1,4   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,19  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,19  ),GZL ,AMP(6   ))            
      CALL FVOXXX(W(1,4   ),W(1,16  ),GWF ,ZERO    ,ZERO    ,W(1,20  ))    
      CALL JIOXXX(W(1,5   ),W(1,20  ),GZN ,ZMASS   ,ZWIDTH  ,W(1,21  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,21  ),GZL ,AMP(7   ))            
      CALL IOVXXX(W(1,8   ),W(1,13  ),W(1,10  ),GAD ,AMP(8   ))            
      CALL IOVXXX(W(1,8   ),W(1,13  ),W(1,12  ),GZD ,AMP(9   ))
!      
      if (.not.with_anom) then      
      CALL JVVXXX(W(1,9   ),W(1,16  ),GWWA ,ZERO    ,ZERO    ,W(1,         
     &     22  ))  
      elseif (with_anom) then
      CALL JVVXXX_anomwaj0(W(1,9   ),W(1,16  ),GWWA ,ZERO    ,ZERO    ,W(1,         
     &     22  ))  
      endif
!
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,22  ),GAL ,AMP(10  ))
!
      if (.not.with_anom) then                  
      CALL JVVXXX(W(1,9   ),W(1,16  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     23  ))                                                          
      elseif (with_anom) then
      CALL JVVXXX_anomwzj(W(1,9   ),W(1,16  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     23  ))                                                          
      endif
!
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,23  ),GZL ,AMP(11  ))            
      CALL FVOXXX(W(1,6   ),W(1,16  ),GWF ,ZERO    ,ZERO    ,W(1,24  ))    
      CALL IOVXXX(W(1,7   ),W(1,24  ),W(1,9   ),GWF ,AMP(12  ))            
      CALL FVIXXX(W(1,3   ),W(1,9   ),GWF ,ZERO    ,ZERO    ,W(1,25  ))    
      CALL FVOXXX(W(1,1   ),W(1,10  ),GAU ,ZERO    ,ZERO    ,W(1,26  ))    
      CALL IOVXXX(W(1,25  ),W(1,26  ),W(1,2   ),GG ,AMP(13  ))             
      CALL FVOXXX(W(1,1   ),W(1,12  ),GZU ,ZERO    ,ZERO    ,W(1,27  ))    
      CALL IOVXXX(W(1,25  ),W(1,27  ),W(1,2   ),GG ,AMP(14  ))             
      CALL FVOXXX(W(1,1   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,28  ))     
      CALL FVOXXX(W(1,28  ),W(1,9   ),GWF ,ZERO    ,ZERO    ,W(1,29  ))    
      CALL IOVXXX(W(1,3   ),W(1,29  ),W(1,10  ),GAD ,AMP(15  ))            
      CALL IOVXXX(W(1,3   ),W(1,29  ),W(1,12  ),GZD ,AMP(16  ))            
      CALL JIOXXX(W(1,3   ),W(1,28  ),GWF ,WMASS   ,WWIDTH  ,W(1,30  ))    
      CALL FVIXXX(W(1,5   ),W(1,30  ),GWF ,ZERO    ,ZERO    ,W(1,31  ))    
      CALL JIOXXX(W(1,31  ),W(1,4   ),GAL ,ZERO    ,ZERO    ,W(1,32  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,32  ),GAL ,AMP(17  ))            
      CALL JIOXXX(W(1,31  ),W(1,4   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,33  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,33  ),GZL ,AMP(18  ))            
      CALL FVOXXX(W(1,4   ),W(1,30  ),GWF ,ZERO    ,ZERO    ,W(1,34  ))    
      CALL JIOXXX(W(1,5   ),W(1,34  ),GZN ,ZMASS   ,ZWIDTH  ,W(1,35  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,35  ),GZL ,AMP(19  ))            
      CALL IOVXXX(W(1,25  ),W(1,28  ),W(1,10  ),GAU ,AMP(20  ))            
      CALL IOVXXX(W(1,25  ),W(1,28  ),W(1,12  ),GZU ,AMP(21  ))
!
      if (.not.with_anom) then      
      CALL JVVXXX(W(1,9   ),W(1,30  ),GWWA ,ZERO    ,ZERO    ,W(1,         
     &     36  ))                                                          
      elseif (with_anom) then
      CALL JVVXXX_anomwaj0(W(1,9   ),W(1,30  ),GWWA ,ZERO    ,ZERO    ,W(1,         
     &     36  ))                                                          
      endif
!
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,36  ),GAL ,AMP(22  ))            
!
      if (.not.with_anom) then            
      CALL JVVXXX(W(1,9   ),W(1,30  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     37  ))  
      elseif (with_anom) then
      CALL JVVXXX_anomwzj(W(1,9   ),W(1,30  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     37  ))  
      endif
!
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,37  ),GZL ,AMP(23  ))            
      CALL FVOXXX(W(1,6   ),W(1,30  ),GWF ,ZERO    ,ZERO    ,W(1,38  ))    
      CALL IOVXXX(W(1,7   ),W(1,38  ),W(1,9   ),GWF ,AMP(24  ))            
      JAMP(   1) = -AMP(   1)-AMP(   2)-AMP(   3)-AMP(   4)-AMP(   5)
     &             -AMP(   6)-AMP(   7)-AMP(   8)-AMP(   9)-AMP(  10)
     &             -AMP(  11)-AMP(  12)-AMP(  13)-AMP(  14)-AMP(  15)
     &             -AMP(  16)-AMP(  17)-AMP(  18)-AMP(  19)-AMP(  20)
     &             -AMP(  21)-AMP(  22)-AMP(  23)-AMP(  24)
      UBGDBEMVEBMUMM = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          UBGDBEMVEBMUMM =UBGDBEMVEBMUMM+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2(i)=amp2(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2(i)=Jamp2(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END


      Double complex FUNCTION AMPUBGDBEMVEBMUMM(P,NHEL)
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : u~ g -> d~ e- ve~ mu- mu+  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=  24,NEIGEN=  1) 
#include "genps.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  38, NCOLOR=   1) 
      Double Precision     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      Double Precision P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I

      Double Precision DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      Double Complex AMP(NGRAPHS)
      Double Complex W(6,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2(maxamps), jamp2(0:maxamps)
      common/to_ampsUBGDBEMVEBMUMMUP/  amp2,       jamp2
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/global.inc"

C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     4/                                  
C               T[1,3,2]                                                   
C ----------
C BEGIN CODE
C ----------
      CALL OXXXXX(P(0,1   ),ZERO ,NHEL(1   ),-1,W(1,1   ))        
      CALL VXXXXX(P(0,2   ),ZERO ,NHEL(2   ),-1,W(1,2   ))        
      CALL IXXXXX(P(0,3   ),ZERO ,NHEL(3   ),-1,W(1,3   ))        
      CALL OXXXXX(P(0,4   ),ZERO ,NHEL(4   ),+1,W(1,4   ))        
      CALL IXXXXX(P(0,5   ),ZERO ,NHEL(5   ),-1,W(1,5   ))        
      CALL OXXXXX(P(0,6   ),ZERO ,NHEL(6   ),+1,W(1,6   ))        
      CALL IXXXXX(P(0,7   ),ZERO ,NHEL(7   ),-1,W(1,7   ))        
      CALL FVIXXX(W(1,3   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,8   ))     
      CALL JIOXXX(W(1,5   ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,9   ))    
      CALL JIOXXX(W(1,7   ),W(1,6   ),GAL ,ZERO    ,ZERO    ,W(1,10  ))    
      CALL FVIXXX(W(1,8   ),W(1,9   ),GWF ,ZERO    ,ZERO    ,W(1,11  ))    
      CALL IOVXXX(W(1,11  ),W(1,1   ),W(1,10  ),GAU ,AMP(1   ))            
      CALL JIOXXX(W(1,7   ),W(1,6   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,12  ))    
      CALL IOVXXX(W(1,11  ),W(1,1   ),W(1,12  ),GZU ,AMP(2   ))            
      CALL FVOXXX(W(1,1   ),W(1,9   ),GWF ,ZERO    ,ZERO    ,W(1,13  ))    
      CALL FVIXXX(W(1,3   ),W(1,10  ),GAD ,ZERO    ,ZERO    ,W(1,14  ))    
      CALL IOVXXX(W(1,14  ),W(1,13  ),W(1,2   ),GG ,AMP(3   ))             
      CALL FVIXXX(W(1,3   ),W(1,12  ),GZD ,ZERO    ,ZERO    ,W(1,15  ))    
      CALL IOVXXX(W(1,15  ),W(1,13  ),W(1,2   ),GG ,AMP(4   ))             
      CALL JIOXXX(W(1,8   ),W(1,1   ),GWF ,WMASS   ,WWIDTH  ,W(1,16  ))    
      CALL FVIXXX(W(1,5   ),W(1,16  ),GWF ,ZERO    ,ZERO    ,W(1,17  ))    
      CALL JIOXXX(W(1,17  ),W(1,4   ),GAL ,ZERO    ,ZERO    ,W(1,18  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,18  ),GAL ,AMP(5   ))            
      CALL JIOXXX(W(1,17  ),W(1,4   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,19  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,19  ),GZL ,AMP(6   ))            
      CALL FVOXXX(W(1,4   ),W(1,16  ),GWF ,ZERO    ,ZERO    ,W(1,20  ))    
      CALL JIOXXX(W(1,5   ),W(1,20  ),GZN ,ZMASS   ,ZWIDTH  ,W(1,21  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,21  ),GZL ,AMP(7   ))            
      CALL IOVXXX(W(1,8   ),W(1,13  ),W(1,10  ),GAD ,AMP(8   ))            
      CALL IOVXXX(W(1,8   ),W(1,13  ),W(1,12  ),GZD ,AMP(9   ))
!      
      if (.not.with_anom) then      
      CALL JVVXXX(W(1,9   ),W(1,16  ),GWWA ,ZERO    ,ZERO    ,W(1,         
     &     22  ))  
      elseif (with_anom) then
      CALL JVVXXX_anomwaj0(W(1,9   ),W(1,16  ),GWWA ,ZERO    ,ZERO    ,W(1,         
     &     22  ))  
      endif
!
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,22  ),GAL ,AMP(10  ))
!
      if (.not.with_anom) then                  
      CALL JVVXXX(W(1,9   ),W(1,16  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     23  ))                                                          
      elseif (with_anom) then
      CALL JVVXXX_anomwzj(W(1,9   ),W(1,16  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     23  ))                                                          
      endif
!
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,23  ),GZL ,AMP(11  ))            
      CALL FVOXXX(W(1,6   ),W(1,16  ),GWF ,ZERO    ,ZERO    ,W(1,24  ))    
      CALL IOVXXX(W(1,7   ),W(1,24  ),W(1,9   ),GWF ,AMP(12  ))            
      CALL FVIXXX(W(1,3   ),W(1,9   ),GWF ,ZERO    ,ZERO    ,W(1,25  ))    
      CALL FVOXXX(W(1,1   ),W(1,10  ),GAU ,ZERO    ,ZERO    ,W(1,26  ))    
      CALL IOVXXX(W(1,25  ),W(1,26  ),W(1,2   ),GG ,AMP(13  ))             
      CALL FVOXXX(W(1,1   ),W(1,12  ),GZU ,ZERO    ,ZERO    ,W(1,27  ))    
      CALL IOVXXX(W(1,25  ),W(1,27  ),W(1,2   ),GG ,AMP(14  ))             
      CALL FVOXXX(W(1,1   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,28  ))     
      CALL FVOXXX(W(1,28  ),W(1,9   ),GWF ,ZERO    ,ZERO    ,W(1,29  ))    
      CALL IOVXXX(W(1,3   ),W(1,29  ),W(1,10  ),GAD ,AMP(15  ))            
      CALL IOVXXX(W(1,3   ),W(1,29  ),W(1,12  ),GZD ,AMP(16  ))            
      CALL JIOXXX(W(1,3   ),W(1,28  ),GWF ,WMASS   ,WWIDTH  ,W(1,30  ))    
      CALL FVIXXX(W(1,5   ),W(1,30  ),GWF ,ZERO    ,ZERO    ,W(1,31  ))    
      CALL JIOXXX(W(1,31  ),W(1,4   ),GAL ,ZERO    ,ZERO    ,W(1,32  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,32  ),GAL ,AMP(17  ))            
      CALL JIOXXX(W(1,31  ),W(1,4   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,33  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,33  ),GZL ,AMP(18  ))            
      CALL FVOXXX(W(1,4   ),W(1,30  ),GWF ,ZERO    ,ZERO    ,W(1,34  ))    
      CALL JIOXXX(W(1,5   ),W(1,34  ),GZN ,ZMASS   ,ZWIDTH  ,W(1,35  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,35  ),GZL ,AMP(19  ))            
      CALL IOVXXX(W(1,25  ),W(1,28  ),W(1,10  ),GAU ,AMP(20  ))            
      CALL IOVXXX(W(1,25  ),W(1,28  ),W(1,12  ),GZU ,AMP(21  ))
!
      if (.not.with_anom) then      
      CALL JVVXXX(W(1,9   ),W(1,30  ),GWWA ,ZERO    ,ZERO    ,W(1,         
     &     36  ))                                                          
      elseif (with_anom) then
      CALL JVVXXX_anomwaj0(W(1,9   ),W(1,30  ),GWWA ,ZERO    ,ZERO    ,W(1,         
     &     36  ))                                                          
      endif
!
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,36  ),GAL ,AMP(22  ))            
!
      if (.not.with_anom) then            
      CALL JVVXXX(W(1,9   ),W(1,30  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     37  ))  
      elseif (with_anom) then
      CALL JVVXXX_anomwzj(W(1,9   ),W(1,30  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     37  ))  
      endif
!
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,37  ),GZL ,AMP(23  ))            
      CALL FVOXXX(W(1,6   ),W(1,30  ),GWF ,ZERO    ,ZERO    ,W(1,38  ))    
      CALL IOVXXX(W(1,7   ),W(1,38  ),W(1,9   ),GWF ,AMP(24  ))            
      AMPUBGDBEMVEBMUMM = (-AMP(   1)-AMP(   2)-AMP(   3)-AMP(   4)-AMP(   5)
     &             -AMP(   6)-AMP(   7)-AMP(   8)-AMP(   9)-AMP(  10)
     &             -AMP(  11)-AMP(  12)-AMP(  13)-AMP(  14)-AMP(  15)
     &             -AMP(  16)-AMP(  17)-AMP(  18)-AMP(  19)-AMP(  20)
     &             -AMP(  21)-AMP(  22)-AMP(  23)-AMP(  24))
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END


       subroutine UBGDBEMVEBMUMMDP(P,UDGAENU)
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : u~ g -> d~ e- ve~ mu- mu+  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=  24,NEIGEN=  1) 
#include "genps.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  38, NCOLOR=   1) 
      Double Precision     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      Double Precision P(0:3,NEXTERNAL)

C  
C LOCAL VARIABLES 
C  
      INTEGER mu,nu

      Double Complex AMP(NGRAPHS)
      Double Complex W(6,NWAVEFUNCS),UDGAENU(0:3,-1:1)
C  
C GLOBAL VARIABLES
C  
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/global.inc"


C ----------
C BEGIN CODE
C ----------
      do mu=0,3
      do nu=-1,1,2

      CALL OXXXXX(P(0,1   ),ZERO ,1,-1,W(1,1   ))        
      CALL VCARTXDP(P(0,2   ),mu,-1,W(1,2   )) 
      CALL IXXXXX(P(0,3   ),ZERO ,1,-1,W(1,3   ))
      CALL OXXXXX(P(0,4   ),ZERO ,-1,1,W(1,4   ))        
      CALL IXXXXX(P(0,5   ),ZERO ,1,-1,W(1,5   ))                
      if (nu.eq.-1) then        
      CALL OXXXXX(P(0,6   ),ZERO ,-1,+1,W(1,6   ))        
      CALL IXXXXX(P(0,7   ),ZERO , 1,-1,W(1,7   ))
      elseif (nu.eq.1) then
      CALL OXXXXX(P(0,6   ),ZERO ,1,+1,W(1,6   ))        
      CALL IXXXXX(P(0,7   ),ZERO ,-1,-1,W(1,7   ))
      endif

      CALL FVIXXX(W(1,3   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,8   ))     
      CALL JIOXXX(W(1,5   ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,9   ))    
      CALL JIOXXX(W(1,7   ),W(1,6   ),GAL ,ZERO    ,ZERO    ,W(1,10  ))    
      CALL FVIXXX(W(1,8   ),W(1,9   ),GWF ,ZERO    ,ZERO    ,W(1,11  ))    
      CALL IOVXXX(W(1,11  ),W(1,1   ),W(1,10  ),GAU ,AMP(1   ))            
      CALL JIOXXX(W(1,7   ),W(1,6   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,12  ))    
      CALL IOVXXX(W(1,11  ),W(1,1   ),W(1,12  ),GZU ,AMP(2   ))            
      CALL FVOXXX(W(1,1   ),W(1,9   ),GWF ,ZERO    ,ZERO    ,W(1,13  ))    
      CALL FVIXXX(W(1,3   ),W(1,10  ),GAD ,ZERO    ,ZERO    ,W(1,14  ))    
      CALL IOVXXX(W(1,14  ),W(1,13  ),W(1,2   ),GG ,AMP(3   ))             
      CALL FVIXXX(W(1,3   ),W(1,12  ),GZD ,ZERO    ,ZERO    ,W(1,15  ))    
      CALL IOVXXX(W(1,15  ),W(1,13  ),W(1,2   ),GG ,AMP(4   ))             
      CALL JIOXXX(W(1,8   ),W(1,1   ),GWF ,WMASS   ,WWIDTH  ,W(1,16  ))    
      CALL FVIXXX(W(1,5   ),W(1,16  ),GWF ,ZERO    ,ZERO    ,W(1,17  ))    
      CALL JIOXXX(W(1,17  ),W(1,4   ),GAL ,ZERO    ,ZERO    ,W(1,18  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,18  ),GAL ,AMP(5   ))            
      CALL JIOXXX(W(1,17  ),W(1,4   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,19  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,19  ),GZL ,AMP(6   ))            
      CALL FVOXXX(W(1,4   ),W(1,16  ),GWF ,ZERO    ,ZERO    ,W(1,20  ))    
      CALL JIOXXX(W(1,5   ),W(1,20  ),GZN ,ZMASS   ,ZWIDTH  ,W(1,21  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,21  ),GZL ,AMP(7   ))            
      CALL IOVXXX(W(1,8   ),W(1,13  ),W(1,10  ),GAD ,AMP(8   ))            
      CALL IOVXXX(W(1,8   ),W(1,13  ),W(1,12  ),GZD ,AMP(9   ))
!      
      if (.not.with_anom) then      
      CALL JVVXXX(W(1,9   ),W(1,16  ),GWWA ,ZERO    ,ZERO    ,W(1,         
     &     22  ))  
      elseif (with_anom) then
      CALL JVVXXX_anomwaj0(W(1,9   ),W(1,16  ),GWWA ,ZERO    ,ZERO    ,W(1,         
     &     22  ))  
      endif
!
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,22  ),GAL ,AMP(10  ))
!
      if (.not.with_anom) then                  
      CALL JVVXXX(W(1,9   ),W(1,16  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     23  ))                                                          
      elseif (with_anom) then
      CALL JVVXXX_anomwzj(W(1,9   ),W(1,16  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     23  ))                                                          
      endif
!
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,23  ),GZL ,AMP(11  ))            
      CALL FVOXXX(W(1,6   ),W(1,16  ),GWF ,ZERO    ,ZERO    ,W(1,24  ))    
      CALL IOVXXX(W(1,7   ),W(1,24  ),W(1,9   ),GWF ,AMP(12  ))            
      CALL FVIXXX(W(1,3   ),W(1,9   ),GWF ,ZERO    ,ZERO    ,W(1,25  ))    
      CALL FVOXXX(W(1,1   ),W(1,10  ),GAU ,ZERO    ,ZERO    ,W(1,26  ))    
      CALL IOVXXX(W(1,25  ),W(1,26  ),W(1,2   ),GG ,AMP(13  ))             
      CALL FVOXXX(W(1,1   ),W(1,12  ),GZU ,ZERO    ,ZERO    ,W(1,27  ))    
      CALL IOVXXX(W(1,25  ),W(1,27  ),W(1,2   ),GG ,AMP(14  ))             
      CALL FVOXXX(W(1,1   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,28  ))     
      CALL FVOXXX(W(1,28  ),W(1,9   ),GWF ,ZERO    ,ZERO    ,W(1,29  ))    
      CALL IOVXXX(W(1,3   ),W(1,29  ),W(1,10  ),GAD ,AMP(15  ))            
      CALL IOVXXX(W(1,3   ),W(1,29  ),W(1,12  ),GZD ,AMP(16  ))            
      CALL JIOXXX(W(1,3   ),W(1,28  ),GWF ,WMASS   ,WWIDTH  ,W(1,30  ))    
      CALL FVIXXX(W(1,5   ),W(1,30  ),GWF ,ZERO    ,ZERO    ,W(1,31  ))    
      CALL JIOXXX(W(1,31  ),W(1,4   ),GAL ,ZERO    ,ZERO    ,W(1,32  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,32  ),GAL ,AMP(17  ))            
      CALL JIOXXX(W(1,31  ),W(1,4   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,33  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,33  ),GZL ,AMP(18  ))            
      CALL FVOXXX(W(1,4   ),W(1,30  ),GWF ,ZERO    ,ZERO    ,W(1,34  ))    
      CALL JIOXXX(W(1,5   ),W(1,34  ),GZN ,ZMASS   ,ZWIDTH  ,W(1,35  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,35  ),GZL ,AMP(19  ))            
      CALL IOVXXX(W(1,25  ),W(1,28  ),W(1,10  ),GAU ,AMP(20  ))            
      CALL IOVXXX(W(1,25  ),W(1,28  ),W(1,12  ),GZU ,AMP(21  ))
!
      if (.not.with_anom) then      
      CALL JVVXXX(W(1,9   ),W(1,30  ),GWWA ,ZERO    ,ZERO    ,W(1,         
     &     36  ))                                                          
      elseif (with_anom) then
      CALL JVVXXX_anomwaj0(W(1,9   ),W(1,30  ),GWWA ,ZERO    ,ZERO    ,W(1,         
     &     36  ))                                                          
      endif
!
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,36  ),GAL ,AMP(22  ))            
!
      if (.not.with_anom) then            
      CALL JVVXXX(W(1,9   ),W(1,30  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     37  ))  
      elseif (with_anom) then
      CALL JVVXXX_anomwzj(W(1,9   ),W(1,30  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     37  ))  
      endif
!
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,37  ),GZL ,AMP(23  ))            
      CALL FVOXXX(W(1,6   ),W(1,30  ),GWF ,ZERO    ,ZERO    ,W(1,38  ))    
      CALL IOVXXX(W(1,7   ),W(1,38  ),W(1,9   ),GWF ,AMP(24  ))        
      UDGAENU(mu,nu) = -AMP(   1)-AMP(   2)-AMP(   3)-AMP(   4)-AMP(   5)
     &             -AMP(   6)-AMP(   7)-AMP(   8)-AMP(   9)-AMP(  10)
     &             -AMP(  11)-AMP(  12)-AMP(  13)-AMP(  14)-AMP(  15)
     &             -AMP(  16)-AMP(  17)-AMP(  18)-AMP(  19)-AMP(  20)
     &             -AMP(  21)-AMP(  22)-AMP(  23)-AMP(  24)
      enddo
      enddo      
      END


      SUBROUTINE SDGUBEMVEBMUMMUP(P1,ANS)
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : d g -> u e- ve~ mu- mu+  
C  
C Crossing   1 is d g -> u e- ve~ mu- mu+  
      IMPLICIT NONE
C  
C CONSTANTS
C  
#include "genps.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB= 128, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      Double Precision P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      Double Precision T
      Double Precision DGUBEMVEBMUMMU
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      Double Precision hwgt, xtry, xrej
      INTEGER idum, ngood, igood(ncomb), jhel, j
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2(maxamps), jamp2(0:maxamps)
      common/to_ampsDGUBEMVEBMUMMUP/  amp2,       jamp2

      character*79         hel_buff
      common/to_helicityDGUBEMVEBMUMMUP/  hel_buff

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixDGUBEMVEBMUMMUP/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigsDGUBEMVEBMUMMUP/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.true./
      SAVE igood, IDUM, jhel
      DATA NGRAPHS /   24/          
      DATA jamp2(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1,7) /-1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1,7) /-1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1,7) /-1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   4),IHEL=1,7) /-1,-1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   5),IHEL=1,7) /-1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1,7) /-1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1,7) /-1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1,7) /-1,-1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1,7) /-1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  10),IHEL=1,7) /-1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  11),IHEL=1,7) /-1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1,7) /-1,-1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1,7) /-1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1,7) /-1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1,7) /-1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  16),IHEL=1,7) /-1,-1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  17),IHEL=1,7) /-1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1,7) /-1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1,7) /-1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1,7) /-1,-1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1,7) /-1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  22),IHEL=1,7) /-1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  23),IHEL=1,7) /-1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1,7) /-1,-1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1,7) /-1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1,7) /-1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1,7) /-1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  28),IHEL=1,7) /-1,-1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  29),IHEL=1,7) /-1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1,7) /-1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1,7) /-1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1,7) /-1,-1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1,7) /-1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  34),IHEL=1,7) /-1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  35),IHEL=1,7) /-1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1,7) /-1, 1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1,7) /-1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1,7) /-1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1,7) /-1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  40),IHEL=1,7) /-1, 1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  41),IHEL=1,7) /-1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1,7) /-1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1,7) /-1, 1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1,7) /-1, 1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1,7) /-1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  46),IHEL=1,7) /-1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  47),IHEL=1,7) /-1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1,7) /-1, 1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  49),IHEL=1,7) /-1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  50),IHEL=1,7) /-1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  51),IHEL=1,7) /-1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  52),IHEL=1,7) /-1, 1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  53),IHEL=1,7) /-1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  54),IHEL=1,7) /-1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  55),IHEL=1,7) /-1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  56),IHEL=1,7) /-1, 1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  57),IHEL=1,7) /-1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  58),IHEL=1,7) /-1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  59),IHEL=1,7) /-1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  60),IHEL=1,7) /-1, 1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  61),IHEL=1,7) /-1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  62),IHEL=1,7) /-1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  63),IHEL=1,7) /-1, 1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  64),IHEL=1,7) /-1, 1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  65),IHEL=1,7) / 1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  66),IHEL=1,7) / 1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  67),IHEL=1,7) / 1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  68),IHEL=1,7) / 1,-1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  69),IHEL=1,7) / 1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  70),IHEL=1,7) / 1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  71),IHEL=1,7) / 1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  72),IHEL=1,7) / 1,-1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  73),IHEL=1,7) / 1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  74),IHEL=1,7) / 1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  75),IHEL=1,7) / 1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  76),IHEL=1,7) / 1,-1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  77),IHEL=1,7) / 1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  78),IHEL=1,7) / 1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  79),IHEL=1,7) / 1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  80),IHEL=1,7) / 1,-1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  81),IHEL=1,7) / 1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  82),IHEL=1,7) / 1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  83),IHEL=1,7) / 1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  84),IHEL=1,7) / 1,-1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  85),IHEL=1,7) / 1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  86),IHEL=1,7) / 1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  87),IHEL=1,7) / 1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  88),IHEL=1,7) / 1,-1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  89),IHEL=1,7) / 1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  90),IHEL=1,7) / 1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  91),IHEL=1,7) / 1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  92),IHEL=1,7) / 1,-1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  93),IHEL=1,7) / 1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  94),IHEL=1,7) / 1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  95),IHEL=1,7) / 1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  96),IHEL=1,7) / 1,-1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  97),IHEL=1,7) / 1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  98),IHEL=1,7) / 1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  99),IHEL=1,7) / 1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 100),IHEL=1,7) / 1, 1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 101),IHEL=1,7) / 1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 102),IHEL=1,7) / 1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 103),IHEL=1,7) / 1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 104),IHEL=1,7) / 1, 1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 105),IHEL=1,7) / 1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 106),IHEL=1,7) / 1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 107),IHEL=1,7) / 1, 1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 108),IHEL=1,7) / 1, 1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 109),IHEL=1,7) / 1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 110),IHEL=1,7) / 1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 111),IHEL=1,7) / 1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 112),IHEL=1,7) / 1, 1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 113),IHEL=1,7) / 1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 114),IHEL=1,7) / 1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 115),IHEL=1,7) / 1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 116),IHEL=1,7) / 1, 1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 117),IHEL=1,7) / 1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 118),IHEL=1,7) / 1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 119),IHEL=1,7) / 1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 120),IHEL=1,7) / 1, 1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 121),IHEL=1,7) / 1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 122),IHEL=1,7) / 1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 123),IHEL=1,7) / 1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 124),IHEL=1,7) / 1, 1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 125),IHEL=1,7) / 1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 126),IHEL=1,7) / 1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 127),IHEL=1,7) / 1, 1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 128),IHEL=1,7) / 1, 1, 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1,7) / 1, 2, 3, 4, 5, 6, 7/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  96/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
!       CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
!       IF (multi_channel) THEN
!           DO IHEL=1,NGRAPHS
!               amp2(ihel)=0d0
!               jamp2(ihel)=0d0
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
                 T=DGUBEMVEBMUMMU(P1 ,NHEL(1,IHEL),JC(1))            
                 ANS(IPROC)=ANS(IPROC)+T
                  IF (T .GT. 0D0 .AND. .NOT. GOODHEL(IHEL,IPROC)) THEN
                      GOODHEL(IHEL,IPROC)=.TRUE.
                      NGOOD = NGOOD +1
                      IGOOD(NGOOD) = IHEL
C                WRITE(*,*) ngood,IHEL,T
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
              T=DGUBEMVEBMUMMU(P1 ,NHEL(1,IHEL),JC(1))            
           ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
!       IF (MULTI_CHANNEL) THEN
!           XTOT=0D0
!           DO IHEL=1,MAPCONFIG(0)
!               XTOT=XTOT+AMP2(MAPCONFIG(IHEL))
!           ENDDO
!           ANS(IPROC)=ANS(IPROC)*AMP2(MAPCONFIG(ICONFIG))/XTOT
!       ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      Double Precision FUNCTION DGUBEMVEBMUMMU(P,NHEL,IC)
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : d g -> u e- ve~ mu- mu+  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=  24,NEIGEN=  1) 
#include "genps.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  38, NCOLOR=   1) 
      Double Precision     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      Double Precision P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      Double Complex ZTEMP
      Double Precision DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      Double Complex AMP(NGRAPHS), JAMP(NCOLOR)
      Double Complex W(6,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2(maxamps), jamp2(0:maxamps)
      common/to_ampsDGUBEMVEBMUMMUP/  amp2,       jamp2
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/global.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     4/                                  
C               T[3,1,2]                                                   
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),ZERO ,NHEL(1   ),+1*IC(1   ),W(1,1   ))        
      CALL VXXXXX(P(0,2   ),ZERO ,NHEL(2   ),-1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),ZERO ,NHEL(3   ),+1*IC(3   ),W(1,3   ))        
      CALL OXXXXX(P(0,4   ),ZERO ,NHEL(4   ),+1*IC(4   ),W(1,4   ))        
      CALL IXXXXX(P(0,5   ),ZERO ,NHEL(5   ),-1*IC(5   ),W(1,5   ))        
      CALL OXXXXX(P(0,6   ),ZERO ,NHEL(6   ),+1*IC(6   ),W(1,6   ))        
      CALL IXXXXX(P(0,7   ),ZERO ,NHEL(7   ),-1*IC(7   ),W(1,7   ))        
      CALL FVOXXX(W(1,3   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,8   ))     
      CALL JIOXXX(W(1,5   ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,9   ))    
      CALL JIOXXX(W(1,7   ),W(1,6   ),GAL ,ZERO    ,ZERO    ,W(1,10  ))    
      CALL FVOXXX(W(1,8   ),W(1,9   ),GWF ,ZERO    ,ZERO    ,W(1,11  ))    
      CALL IOVXXX(W(1,1   ),W(1,11  ),W(1,10  ),GAD ,AMP(1   ))            
      CALL JIOXXX(W(1,7   ),W(1,6   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,12  ))    
      CALL IOVXXX(W(1,1   ),W(1,11  ),W(1,12  ),GZD ,AMP(2   ))            
      CALL FVIXXX(W(1,1   ),W(1,9   ),GWF ,ZERO    ,ZERO    ,W(1,13  ))    
      CALL FVOXXX(W(1,3   ),W(1,10  ),GAU ,ZERO    ,ZERO    ,W(1,14  ))    
      CALL IOVXXX(W(1,13  ),W(1,14  ),W(1,2   ),GG ,AMP(3   ))             
      CALL FVOXXX(W(1,3   ),W(1,12  ),GZU ,ZERO    ,ZERO    ,W(1,15  ))    
      CALL IOVXXX(W(1,13  ),W(1,15  ),W(1,2   ),GG ,AMP(4   ))             
      CALL JIOXXX(W(1,1   ),W(1,8   ),GWF ,WMASS   ,WWIDTH  ,W(1,16  ))    
      CALL FVIXXX(W(1,5   ),W(1,16  ),GWF ,ZERO    ,ZERO    ,W(1,17  ))    
      CALL JIOXXX(W(1,17  ),W(1,4   ),GAL ,ZERO    ,ZERO    ,W(1,18  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,18  ),GAL ,AMP(5   ))            
      CALL JIOXXX(W(1,17  ),W(1,4   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,19  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,19  ),GZL ,AMP(6   ))            
      CALL FVOXXX(W(1,4   ),W(1,16  ),GWF ,ZERO    ,ZERO    ,W(1,20  ))    
      CALL JIOXXX(W(1,5   ),W(1,20  ),GZN ,ZMASS   ,ZWIDTH  ,W(1,21  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,21  ),GZL ,AMP(7   ))            
      CALL IOVXXX(W(1,13  ),W(1,8   ),W(1,10  ),GAU ,AMP(8   ))            
      CALL IOVXXX(W(1,13  ),W(1,8   ),W(1,12  ),GZU ,AMP(9   ))
!
      if (.not.with_anom) then
      CALL JVVXXX(W(1,9   ),W(1,16  ),GWWA ,ZERO    ,ZERO    ,W(1,         
     &     22  ))  
      elseif (with_anom) then
      CALL JVVXXX_anomwaj0(W(1,9   ),W(1,16  ),GWWA ,ZERO    ,ZERO    ,W(1,         
     &     22  ))  
      endif
!
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,22  ),GAL ,AMP(10  ))
!            
      if (.not.with_anom) then
      CALL JVVXXX(W(1,9   ),W(1,16  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     23  ))
      elseif (with_anom) then
      CALL JVVXXX_anomwzj(W(1,9   ),W(1,16  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     23  ))
      endif
!
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,23  ),GZL ,AMP(11  ))            
      CALL FVOXXX(W(1,6   ),W(1,16  ),GWF ,ZERO    ,ZERO    ,W(1,24  ))    
      CALL IOVXXX(W(1,7   ),W(1,24  ),W(1,9   ),GWF ,AMP(12  ))            
      CALL FVOXXX(W(1,3   ),W(1,9   ),GWF ,ZERO    ,ZERO    ,W(1,25  ))    
      CALL FVIXXX(W(1,1   ),W(1,10  ),GAD ,ZERO    ,ZERO    ,W(1,26  ))    
      CALL IOVXXX(W(1,26  ),W(1,25  ),W(1,2   ),GG ,AMP(13  ))             
      CALL FVIXXX(W(1,1   ),W(1,12  ),GZD ,ZERO    ,ZERO    ,W(1,27  ))    
      CALL IOVXXX(W(1,27  ),W(1,25  ),W(1,2   ),GG ,AMP(14  ))             
      CALL FVIXXX(W(1,1   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,28  ))     
      CALL FVIXXX(W(1,28  ),W(1,9   ),GWF ,ZERO    ,ZERO    ,W(1,29  ))    
      CALL IOVXXX(W(1,29  ),W(1,3   ),W(1,10  ),GAU ,AMP(15  ))            
      CALL IOVXXX(W(1,29  ),W(1,3   ),W(1,12  ),GZU ,AMP(16  ))            
      CALL JIOXXX(W(1,28  ),W(1,3   ),GWF ,WMASS   ,WWIDTH  ,W(1,30  ))    
      CALL FVIXXX(W(1,5   ),W(1,30  ),GWF ,ZERO    ,ZERO    ,W(1,31  ))    
      CALL JIOXXX(W(1,31  ),W(1,4   ),GAL ,ZERO    ,ZERO    ,W(1,32  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,32  ),GAL ,AMP(17  ))            
      CALL JIOXXX(W(1,31  ),W(1,4   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,33  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,33  ),GZL ,AMP(18  ))            
      CALL FVOXXX(W(1,4   ),W(1,30  ),GWF ,ZERO    ,ZERO    ,W(1,34  ))    
      CALL JIOXXX(W(1,5   ),W(1,34  ),GZN ,ZMASS   ,ZWIDTH  ,W(1,35  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,35  ),GZL ,AMP(19  ))            
      CALL IOVXXX(W(1,28  ),W(1,25  ),W(1,10  ),GAD ,AMP(20  ))            
      CALL IOVXXX(W(1,28  ),W(1,25  ),W(1,12  ),GZD ,AMP(21  ))
!
      if (.not.with_anom) then
      CALL JVVXXX(W(1,9   ),W(1,30  ),GWWA ,ZERO    ,ZERO    ,W(1,         
     &     36  ))                 
      elseif (with_anom) then
      CALL JVVXXX_anomwaj0(W(1,9   ),W(1,30  ),GWWA ,ZERO    ,ZERO    ,W(1,         
     &     36  ))                 
      endif
!
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,36  ),GAL ,AMP(22  ))            
!
      if (.not.with_anom) then
      CALL JVVXXX(W(1,9   ),W(1,30  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     37  ))   
      elseif (with_anom) then
      CALL JVVXXX_anomwzj(W(1,9   ),W(1,30  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     37  ))   
      endif
!
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,37  ),GZL ,AMP(23  ))            
      CALL FVOXXX(W(1,6   ),W(1,30  ),GWF ,ZERO    ,ZERO    ,W(1,38  ))    
      CALL IOVXXX(W(1,7   ),W(1,38  ),W(1,9   ),GWF ,AMP(24  ))            
      JAMP(   1) = +AMP(   1)+AMP(   2)+AMP(   3)+AMP(   4)+AMP(   5)
     &             +AMP(   6)+AMP(   7)+AMP(   8)+AMP(   9)+AMP(  10)
     &             +AMP(  11)+AMP(  12)+AMP(  13)+AMP(  14)+AMP(  15)
     &             +AMP(  16)+AMP(  17)+AMP(  18)+AMP(  19)+AMP(  20)
     &             +AMP(  21)+AMP(  22)+AMP(  23)+AMP(  24)
      DGUBEMVEBMUMMU = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          DGUBEMVEBMUMMU =DGUBEMVEBMUMMU+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2(i)=amp2(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2(i)=Jamp2(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END


      Double complex FUNCTION AMPDGUBEMVEBMUMMU(P,NHEL)
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : d g -> u e- ve~ mu- mu+  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=  24,NEIGEN=  1) 
#include "genps.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  38, NCOLOR=   1) 
      Double Precision     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      Double Precision P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I

      Double Precision DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      Double Complex AMP(NGRAPHS)
      Double Complex W(6,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2(maxamps), jamp2(0:maxamps)
      common/to_ampsDGUBEMVEBMUMMUP/  amp2,       jamp2
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/global.inc"

C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     4/                                  
C               T[3,1,2]                                                   
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),ZERO ,NHEL(1   ),+1,W(1,1   ))        
      CALL VXXXXX(P(0,2   ),ZERO ,NHEL(2   ),-1,W(1,2   ))        
      CALL OXXXXX(P(0,3   ),ZERO ,NHEL(3   ),+1,W(1,3   ))        
      CALL OXXXXX(P(0,4   ),ZERO ,NHEL(4   ),+1,W(1,4   ))        
      CALL IXXXXX(P(0,5   ),ZERO ,NHEL(5   ),-1,W(1,5   ))        
      CALL OXXXXX(P(0,6   ),ZERO ,NHEL(6   ),+1,W(1,6   ))        
      CALL IXXXXX(P(0,7   ),ZERO ,NHEL(7   ),-1,W(1,7   ))        
      CALL FVOXXX(W(1,3   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,8   ))     
      CALL JIOXXX(W(1,5   ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,9   ))    
      CALL JIOXXX(W(1,7   ),W(1,6   ),GAL ,ZERO    ,ZERO    ,W(1,10  ))    
      CALL FVOXXX(W(1,8   ),W(1,9   ),GWF ,ZERO    ,ZERO    ,W(1,11  ))    
      CALL IOVXXX(W(1,1   ),W(1,11  ),W(1,10  ),GAD ,AMP(1   ))            
      CALL JIOXXX(W(1,7   ),W(1,6   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,12  ))    
      CALL IOVXXX(W(1,1   ),W(1,11  ),W(1,12  ),GZD ,AMP(2   ))            
      CALL FVIXXX(W(1,1   ),W(1,9   ),GWF ,ZERO    ,ZERO    ,W(1,13  ))    
      CALL FVOXXX(W(1,3   ),W(1,10  ),GAU ,ZERO    ,ZERO    ,W(1,14  ))    
      CALL IOVXXX(W(1,13  ),W(1,14  ),W(1,2   ),GG ,AMP(3   ))             
      CALL FVOXXX(W(1,3   ),W(1,12  ),GZU ,ZERO    ,ZERO    ,W(1,15  ))    
      CALL IOVXXX(W(1,13  ),W(1,15  ),W(1,2   ),GG ,AMP(4   ))             
      CALL JIOXXX(W(1,1   ),W(1,8   ),GWF ,WMASS   ,WWIDTH  ,W(1,16  ))    
      CALL FVIXXX(W(1,5   ),W(1,16  ),GWF ,ZERO    ,ZERO    ,W(1,17  ))    
      CALL JIOXXX(W(1,17  ),W(1,4   ),GAL ,ZERO    ,ZERO    ,W(1,18  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,18  ),GAL ,AMP(5   ))            
      CALL JIOXXX(W(1,17  ),W(1,4   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,19  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,19  ),GZL ,AMP(6   ))            
      CALL FVOXXX(W(1,4   ),W(1,16  ),GWF ,ZERO    ,ZERO    ,W(1,20  ))    
      CALL JIOXXX(W(1,5   ),W(1,20  ),GZN ,ZMASS   ,ZWIDTH  ,W(1,21  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,21  ),GZL ,AMP(7   ))            
      CALL IOVXXX(W(1,13  ),W(1,8   ),W(1,10  ),GAU ,AMP(8   ))            
      CALL IOVXXX(W(1,13  ),W(1,8   ),W(1,12  ),GZU ,AMP(9   ))
!
      if (.not.with_anom) then
      CALL JVVXXX(W(1,9   ),W(1,16  ),GWWA ,ZERO    ,ZERO    ,W(1,         
     &     22  ))  
      elseif (with_anom) then
      CALL JVVXXX_anomwaj0(W(1,9   ),W(1,16  ),GWWA ,ZERO    ,ZERO    ,W(1,         
     &     22  ))  
      endif
!
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,22  ),GAL ,AMP(10  ))
!            
      if (.not.with_anom) then
      CALL JVVXXX(W(1,9   ),W(1,16  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     23  ))
      elseif (with_anom) then
      CALL JVVXXX_anomwzj(W(1,9   ),W(1,16  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     23  ))
      endif
!
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,23  ),GZL ,AMP(11  ))            
      CALL FVOXXX(W(1,6   ),W(1,16  ),GWF ,ZERO    ,ZERO    ,W(1,24  ))    
      CALL IOVXXX(W(1,7   ),W(1,24  ),W(1,9   ),GWF ,AMP(12  ))            
      CALL FVOXXX(W(1,3   ),W(1,9   ),GWF ,ZERO    ,ZERO    ,W(1,25  ))    
      CALL FVIXXX(W(1,1   ),W(1,10  ),GAD ,ZERO    ,ZERO    ,W(1,26  ))    
      CALL IOVXXX(W(1,26  ),W(1,25  ),W(1,2   ),GG ,AMP(13  ))             
      CALL FVIXXX(W(1,1   ),W(1,12  ),GZD ,ZERO    ,ZERO    ,W(1,27  ))    
      CALL IOVXXX(W(1,27  ),W(1,25  ),W(1,2   ),GG ,AMP(14  ))             
      CALL FVIXXX(W(1,1   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,28  ))     
      CALL FVIXXX(W(1,28  ),W(1,9   ),GWF ,ZERO    ,ZERO    ,W(1,29  ))    
      CALL IOVXXX(W(1,29  ),W(1,3   ),W(1,10  ),GAU ,AMP(15  ))            
      CALL IOVXXX(W(1,29  ),W(1,3   ),W(1,12  ),GZU ,AMP(16  ))            
      CALL JIOXXX(W(1,28  ),W(1,3   ),GWF ,WMASS   ,WWIDTH  ,W(1,30  ))    
      CALL FVIXXX(W(1,5   ),W(1,30  ),GWF ,ZERO    ,ZERO    ,W(1,31  ))    
      CALL JIOXXX(W(1,31  ),W(1,4   ),GAL ,ZERO    ,ZERO    ,W(1,32  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,32  ),GAL ,AMP(17  ))            
      CALL JIOXXX(W(1,31  ),W(1,4   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,33  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,33  ),GZL ,AMP(18  ))            
      CALL FVOXXX(W(1,4   ),W(1,30  ),GWF ,ZERO    ,ZERO    ,W(1,34  ))    
      CALL JIOXXX(W(1,5   ),W(1,34  ),GZN ,ZMASS   ,ZWIDTH  ,W(1,35  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,35  ),GZL ,AMP(19  ))            
      CALL IOVXXX(W(1,28  ),W(1,25  ),W(1,10  ),GAD ,AMP(20  ))            
      CALL IOVXXX(W(1,28  ),W(1,25  ),W(1,12  ),GZD ,AMP(21  ))
!
      if (.not.with_anom) then
      CALL JVVXXX(W(1,9   ),W(1,30  ),GWWA ,ZERO    ,ZERO    ,W(1,         
     &     36  ))                 
      elseif (with_anom) then
      CALL JVVXXX_anomwaj0(W(1,9   ),W(1,30  ),GWWA ,ZERO    ,ZERO    ,W(1,         
     &     36  ))                 
      endif
!
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,36  ),GAL ,AMP(22  ))            
!
      if (.not.with_anom) then
      CALL JVVXXX(W(1,9   ),W(1,30  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     37  ))   
      elseif (with_anom) then
      CALL JVVXXX_anomwzj(W(1,9   ),W(1,30  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     37  ))   
      endif
!
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,37  ),GZL ,AMP(23  ))            
      CALL FVOXXX(W(1,6   ),W(1,30  ),GWF ,ZERO    ,ZERO    ,W(1,38  ))    
      CALL IOVXXX(W(1,7   ),W(1,38  ),W(1,9   ),GWF ,AMP(24  ))            
      AMPDGUBEMVEBMUMMU = +AMP(   1)+AMP(   2)+AMP(   3)+AMP(   4)+AMP(   5)
     &             +AMP(   6)+AMP(   7)+AMP(   8)+AMP(   9)+AMP(  10)
     &             +AMP(  11)+AMP(  12)+AMP(  13)+AMP(  14)+AMP(  15)
     &             +AMP(  16)+AMP(  17)+AMP(  18)+AMP(  19)+AMP(  20)
     &             +AMP(  21)+AMP(  22)+AMP(  23)+AMP(  24)
      END


      SUBROUTINE DGUBEMVEBMUMMUDP(P,UDGAENU)
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : d g -> u e- ve~ mu- mu+  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=  24,NEIGEN=  1) 
#include "genps.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  38, NCOLOR=   1) 
      Double Precision     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      Double Precision P(0:3,NEXTERNAL)
      INTEGER mu,nu
C  
C LOCAL VARIABLES 
C


      Double Complex AMP(NGRAPHS)
      Double Complex W(6,NWAVEFUNCS),UDGAENU(0:3,-1:1)
C  
C GLOBAL VARIABLES
C  
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/global.inc"

C ----------
C BEGIN CODE
C ----------
      do mu=0,3
      do nu=-1,1,2

      CALL IXXXXX(P(0,1   ),ZERO ,-1,+1,W(1,1   ))        
      CALL VCARTXDP(P(0,2   ),mu,-1,W(1,2   ))       
      CALL OXXXXX(P(0,3   ),ZERO ,-1,+1,W(1,3   ))        
      CALL OXXXXX(P(0,4   ),ZERO ,-1,+1,W(1,4   ))        
      CALL IXXXXX(P(0,5   ),ZERO ,1,-1,W(1,5   ))
      if (nu.eq.-1) then        
      CALL OXXXXX(P(0,6   ),ZERO ,-1,+1,W(1,6   ))        
      CALL IXXXXX(P(0,7   ),ZERO , 1,-1,W(1,7   ))
      elseif (nu.eq.1) then
      CALL OXXXXX(P(0,6   ),ZERO ,1,+1,W(1,6   ))        
      CALL IXXXXX(P(0,7   ),ZERO ,-1,-1,W(1,7   ))
      endif

      CALL FVOXXX(W(1,3   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,8   ))     
      CALL JIOXXX(W(1,5   ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,9   ))    
      CALL JIOXXX(W(1,7   ),W(1,6   ),GAL ,ZERO    ,ZERO    ,W(1,10  ))    
      CALL FVOXXX(W(1,8   ),W(1,9   ),GWF ,ZERO    ,ZERO    ,W(1,11  ))    
      CALL IOVXXX(W(1,1   ),W(1,11  ),W(1,10  ),GAD ,AMP(1   ))            
      CALL JIOXXX(W(1,7   ),W(1,6   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,12  ))    
      CALL IOVXXX(W(1,1   ),W(1,11  ),W(1,12  ),GZD ,AMP(2   ))            
      CALL FVIXXX(W(1,1   ),W(1,9   ),GWF ,ZERO    ,ZERO    ,W(1,13  ))    
      CALL FVOXXX(W(1,3   ),W(1,10  ),GAU ,ZERO    ,ZERO    ,W(1,14  ))    
      CALL IOVXXX(W(1,13  ),W(1,14  ),W(1,2   ),GG ,AMP(3   ))             
      CALL FVOXXX(W(1,3   ),W(1,12  ),GZU ,ZERO    ,ZERO    ,W(1,15  ))    
      CALL IOVXXX(W(1,13  ),W(1,15  ),W(1,2   ),GG ,AMP(4   ))             
      CALL JIOXXX(W(1,1   ),W(1,8   ),GWF ,WMASS   ,WWIDTH  ,W(1,16  ))    
      CALL FVIXXX(W(1,5   ),W(1,16  ),GWF ,ZERO    ,ZERO    ,W(1,17  ))    
      CALL JIOXXX(W(1,17  ),W(1,4   ),GAL ,ZERO    ,ZERO    ,W(1,18  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,18  ),GAL ,AMP(5   ))            
      CALL JIOXXX(W(1,17  ),W(1,4   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,19  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,19  ),GZL ,AMP(6   ))            
      CALL FVOXXX(W(1,4   ),W(1,16  ),GWF ,ZERO    ,ZERO    ,W(1,20  ))    
      CALL JIOXXX(W(1,5   ),W(1,20  ),GZN ,ZMASS   ,ZWIDTH  ,W(1,21  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,21  ),GZL ,AMP(7   ))            
      CALL IOVXXX(W(1,13  ),W(1,8   ),W(1,10  ),GAU ,AMP(8   ))            
      CALL IOVXXX(W(1,13  ),W(1,8   ),W(1,12  ),GZU ,AMP(9   ))
!
      if (.not.with_anom) then
      CALL JVVXXX(W(1,9   ),W(1,16  ),GWWA ,ZERO    ,ZERO    ,W(1,         
     &     22  ))  
      elseif (with_anom) then
      CALL JVVXXX_anomwaj0(W(1,9   ),W(1,16  ),GWWA ,ZERO    ,ZERO    ,W(1,         
     &     22  ))  
      endif
!
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,22  ),GAL ,AMP(10  ))
!            
      if (.not.with_anom) then
      CALL JVVXXX(W(1,9   ),W(1,16  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     23  ))
      elseif (with_anom) then
      CALL JVVXXX_anomwzj(W(1,9   ),W(1,16  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     23  ))
      endif
!
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,23  ),GZL ,AMP(11  ))            
      CALL FVOXXX(W(1,6   ),W(1,16  ),GWF ,ZERO    ,ZERO    ,W(1,24  ))    
      CALL IOVXXX(W(1,7   ),W(1,24  ),W(1,9   ),GWF ,AMP(12  ))            
      CALL FVOXXX(W(1,3   ),W(1,9   ),GWF ,ZERO    ,ZERO    ,W(1,25  ))    
      CALL FVIXXX(W(1,1   ),W(1,10  ),GAD ,ZERO    ,ZERO    ,W(1,26  ))    
      CALL IOVXXX(W(1,26  ),W(1,25  ),W(1,2   ),GG ,AMP(13  ))             
      CALL FVIXXX(W(1,1   ),W(1,12  ),GZD ,ZERO    ,ZERO    ,W(1,27  ))    
      CALL IOVXXX(W(1,27  ),W(1,25  ),W(1,2   ),GG ,AMP(14  ))             
      CALL FVIXXX(W(1,1   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,28  ))     
      CALL FVIXXX(W(1,28  ),W(1,9   ),GWF ,ZERO    ,ZERO    ,W(1,29  ))    
      CALL IOVXXX(W(1,29  ),W(1,3   ),W(1,10  ),GAU ,AMP(15  ))            
      CALL IOVXXX(W(1,29  ),W(1,3   ),W(1,12  ),GZU ,AMP(16  ))            
      CALL JIOXXX(W(1,28  ),W(1,3   ),GWF ,WMASS   ,WWIDTH  ,W(1,30  ))    
      CALL FVIXXX(W(1,5   ),W(1,30  ),GWF ,ZERO    ,ZERO    ,W(1,31  ))    
      CALL JIOXXX(W(1,31  ),W(1,4   ),GAL ,ZERO    ,ZERO    ,W(1,32  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,32  ),GAL ,AMP(17  ))            
      CALL JIOXXX(W(1,31  ),W(1,4   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,33  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,33  ),GZL ,AMP(18  ))            
      CALL FVOXXX(W(1,4   ),W(1,30  ),GWF ,ZERO    ,ZERO    ,W(1,34  ))    
      CALL JIOXXX(W(1,5   ),W(1,34  ),GZN ,ZMASS   ,ZWIDTH  ,W(1,35  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,35  ),GZL ,AMP(19  ))            
      CALL IOVXXX(W(1,28  ),W(1,25  ),W(1,10  ),GAD ,AMP(20  ))            
      CALL IOVXXX(W(1,28  ),W(1,25  ),W(1,12  ),GZD ,AMP(21  ))
!
      if (.not.with_anom) then
      CALL JVVXXX(W(1,9   ),W(1,30  ),GWWA ,ZERO    ,ZERO    ,W(1,         
     &     36  ))                 
      elseif (with_anom) then
      CALL JVVXXX_anomwaj0(W(1,9   ),W(1,30  ),GWWA ,ZERO    ,ZERO    ,W(1,         
     &     36  ))                 
      endif
!
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,36  ),GAL ,AMP(22  ))            
!
      if (.not.with_anom) then
      CALL JVVXXX(W(1,9   ),W(1,30  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     37  ))   
      elseif (with_anom) then
      CALL JVVXXX_anomwzj(W(1,9   ),W(1,30  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     37  ))   
      endif
!
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,37  ),GZL ,AMP(23  ))            
      CALL FVOXXX(W(1,6   ),W(1,30  ),GWF ,ZERO    ,ZERO    ,W(1,38  ))    
      CALL IOVXXX(W(1,7   ),W(1,38  ),W(1,9   ),GWF ,AMP(24  ))            
      UDGAENU(mu,nu) = -AMP(   1)-AMP(   2)-AMP(   3)-AMP(   4)-AMP(   5)
     &             -AMP(   6)-AMP(   7)-AMP(   8)-AMP(   9)-AMP(  10)
     &             -AMP(  11)-AMP(  12)-AMP(  13)-AMP(  14)-AMP(  15)
     &             -AMP(  16)-AMP(  17)-AMP(  18)-AMP(  19)-AMP(  20)
     &             -AMP(  21)-AMP(  22)-AMP(  23)-AMP(  24)
      enddo
      enddo
      END
