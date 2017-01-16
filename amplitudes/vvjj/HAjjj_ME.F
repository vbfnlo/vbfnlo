
      subroutine HAjjj_ME(Isign,gsign,uucc,uuss,ddcc,ddss,udsc,ducs)
          use globalvars, only: ldoblha

      IMPLICIT NONE

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/coupl.inc"
*#include "VBFNLO/utilities/mssm.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

      integer Isign(4+n_v),gsign,iflavF1,iflavF3,isigF,isigA,isigG,iorder,
     $     iinsA,iinsG,k,i,mu,iisig(3,4)
      double precision uucc(3),uuss(3),ddcc(3),ddss(3),udsc(3),ducs(3),
     $     res(6,0:2),fac,dotrr

      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)


      double complex Higgs(6),Curr12(6,3:4,-1:1,2:3),Curr21(6,3:4,-1:1,2:3),
     $     Curr34(6,3:4,-1:1,2:3),Curr43(6,3:4,-1:1,2:3),
     $     Curr12A(6,3:4,-1:1,2,-1:1,2:3),Curr21A(6,3:4,-1:1,2,-1:1,2:3),
     $     Curr34A(6,3:4,-1:1,3:4,-1:1,2:3),Curr43A(6,3:4,-1:1,3:4,-1:1,2:3),
     $     Curr12GO(6,3:4,-1:1,4,-1:1,2:3),Curr21GO(6,3:4,-1:1,4,-1:1,2:3),
     $     Curr34GO(6,3:4,-1:1,4,-1:1,2:3),Curr43GO(6,3:4,-1:1,4,-1:1,2:3),
     $     Curr21GI(6,3:4,-1:1,4,-1:1,2:3),Curr43GI(6,3:4,-1:1,4,-1:1,2:3),
     $     Curr12GOt1Z(6,3:4,-1:1,-1:1,-1:1,2,2),
     $     Curr21GOt1Z(6,3:4,-1:1,-1:1,-1:1,2,2),
     $     Curr34GOt1Z(6,3:4,-1:1,-1:1,-1:1,3:4,3:4),
     $     Curr43GOt1Z(6,3:4,-1:1,-1:1,-1:1,3:4,3:4),
     $     Curr12GOt1W(6,3:4,-1:1,-1:1,-1:1,2,2),
     $     Curr21GOt1W(6,3:4,-1:1,-1:1,-1:1,2,2),
     $     Curr34GOt1W(6,3:4,-1:1,-1:1,-1:1,3:4,3:4),
     $     Curr43GOt1W(6,3:4,-1:1,-1:1,-1:1,3:4,3:4),
     $     Curr12GOt2Z(6,3:4,-1:1,-1:1,-1:1,2,2),
     $     Curr21GOt2Z(6,3:4,-1:1,-1:1,-1:1,2,2),
     $     Curr34GOt2Z(6,3:4,-1:1,-1:1,-1:1,3:4,2),
     $     Curr43GOt2Z(6,3:4,-1:1,-1:1,-1:1,3:4,2),
     $     Curr12GOt2W(6,3:4,-1:1,-1:1,-1:1,2,2),
     $     Curr21GOt2W(6,3:4,-1:1,-1:1,-1:1,2,2),
     $     Curr34GOt2W(6,3:4,-1:1,-1:1,-1:1,3:4,2),
     $     Curr43GOt2W(6,3:4,-1:1,-1:1,-1:1,3:4,2),
     $     Curr21GIt1Z(6,3:4,-1:1,-1:1,-1:1,2,2),
     $     Curr43GIt1Z(6,3:4,-1:1,-1:1,-1:1,3:4,3:4),
     $     Curr21GIt1W(6,3:4,-1:1,-1:1,-1:1,2,2),
     $     Curr43GIt1W(6,3:4,-1:1,-1:1,-1:1,3:4,3:4),
     $     Curr21GIt2Z(6,3:4,-1:1,-1:1,-1:1,2,2),
     $     Curr43GIt2Z(6,3:4,-1:1,-1:1,-1:1,3:4,2),
     $     Curr21GIt2W(6,3:4,-1:1,-1:1,-1:1,2,2),
     $     Curr43GIt2W(6,3:4,-1:1,-1:1,-1:1,3:4,2),
     $     Curr12tbv(6,3:4,-1:1,-1:1),Curr21tbv(6,3:4,-1:1,-1:1),
     $     Curr34tbv(6,3:4,-1:1,-1:1),Curr43tbv(6,3:4,-1:1,-1:1),
     $     Curr12GOtbv(6,3:4,-1:1,-1:1,2,-1:1),
     $     Curr21GOtbv(6,3:4,-1:1,-1:1,2,-1:1),
     $     Curr34GOtbv(6,3:4,-1:1,-1:1,3:4,-1:1),
     $     Curr43GOtbv(6,3:4,-1:1,-1:1,3:4,-1:1),
     $     Curr21GItbv(6,3:4,-1:1,-1:1,2,-1:1),
     $     Curr43GItbv(6,3:4,-1:1,-1:1,3:4,-1:1),
     $     SubAmpW(3:4,3:4,4,-1:1,-1:1,4,4),
     $     SubAmpZ(3:4,3:4,4,-1:1,-1:1,4,4),
     $     AmpSumZZ(3:4,3:4,4,-1:1,-1:1,2),
     $     AmpSumWW(3:4,3:4,4,-1:1,-1:1,2),
     $     SubAmpTBVu(3:4,3:4,4,-1:1,4,-1:1),
     $     SubAmpTBVl(3:4,3:4,4,-1:1,4,-1:1),
     $     Ampbuf(2),czero

      COMMON /HAJJJCURRENTS/ fac,Higgs,Curr12,Curr21,Curr34,Curr43,
     $     Curr12A,Curr21A,Curr34A,Curr43A,
     $     Curr12GO,Curr21GO,Curr34GO,Curr43GO,Curr21GI,Curr43GI,
     $     Curr12GOt1Z,Curr21GOt1Z,Curr34GOt1Z,Curr43GOt1Z,Curr21GIt1Z,Curr43GIt1Z,
     $     Curr12GOt2Z,Curr21GOt2Z,Curr34GOt2Z,Curr43GOt2Z,Curr21GIt2Z,Curr43GIt2Z,
     $     Curr12GOt1W,Curr21GOt1W,Curr34GOt1W,Curr43GOt1W,Curr21GIt1W,Curr43GIt1W,
     $     Curr12GOt2W,Curr21GOt2W,Curr34GOt2W,Curr43GOt2W,Curr21GIt2W,Curr43GIt2W,
     $     Curr12tbv,Curr21tbv,Curr34tbv,Curr43tbv,
     $     Curr12GOtbv,Curr21GOtbv,Curr34GOtbv,Curr43GOtbv,Curr21GItbv,Curr43GItbv

      data iisig /-1,0,-1,  -1,0,1,   1,0,-1,   1,0,1/
      parameter(czero=(0d0,0d0))
*      parameter ( xmtau = 1.777d0 )
*      parameter ( xmmu = 1.106d0 )
      external dotrr

c$$$      ph(0)=-dreal(Higgs(2))
c$$$      ph(1)=-dreal(Higgs(3))
c$$$      ph(2)=-dimag(Higgs(3))
c$$$      ph(3)=-dimag(Higgs(2))
c$$$      ph(4)=dotrr(ph,ph)


*      xmtau = ML
*      xmmu = MM


      do iflavF1=3,4
         do iflavF3=3,4
            do isigF=1,4
               do isigA=-1,1,2
                  do isigG=-1,1,2
                     do k=1,2
                        AmpSumZZ(iflavF1,iflavF3,isigF,isigA,isigG,k)=czero
                        AmpSumWW(iflavF1,iflavF3,isigF,isigA,isigG,k)=czero
                     enddo
                     do iinsA=1,4
                        do iinsG=1,4
                           SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,
     -                          iinsA,iinsG) = czero
                           SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,
     -                          iinsA,iinsG) = czero
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      do i=1,6
         do k=1,2
            res(i,k)=0d0
         enddo
      enddo
      if(gsign.eq.1)then
*******************************************************************************
******************************************************************************* 
********************************************************************************

**SOPHY NOTE FOR IFORT:  originally, if line ended at 'and' and the 3rd comment 
**                       line was uncommented
      if((Isign(1).eq.1).and.(Isign(4).eq.1))then
c     $     Isign(2).eq.1.and.
c     $     Isign(3).eq.1.and.
cS     $     Isign(4).eq.1)then
********************************************************************************

c     now combine two currents with the higgs emission:
c     Curr12(34)A/GO cab be combined with Curr34(12)GO/A
c     Curr12(34)GOt1_ can be combined with Curr34(12) and
c     Curr12(34)GOt2_ can be combined with Curr34(12)
c     The result is
c     SubAmp_(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,iinsG)
         do isigA=-1,1,2
            do isigG=-1,1,2
               do iflavF1=3,4
                  do iflavF3=3,4
                  !Z contributions
                     do isigF=1,4
                        CALL VVSXXX(Curr12GOt1Z(1,iflavF1,iisig(1,isigF),
     -                    isigA,isigG,1,2),Curr34(1,iflavF3,
     -                    iisig(3,isigF),2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,
     -                    isigG,1,2)) !iinsA=1 iinsG=2
                        CALL VVSXXX(Curr12A(1,iflavF1,iisig(1,isigF),1,
     -                    isigA,2),Curr34GO(1,iflavF3,iisig(3,isigF),3,
     -                    isigG,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,
     -                    isigG,1,3)) !iinsA=1 iinsG=3
                        CALL VVSXXX(Curr12A(1,iflavF1,iisig(1,isigF),1,
     -                    isigA,2),Curr34GO(1,iflavF3,iisig(3,isigF),4,
     -                    isigG,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,
     -                    isigG,1,4)) !iinsA=1 iinsG=4
                        CALL VVSXXX(Curr12GOt1Z(1,iflavF1,
     -                    iisig(1,isigF),isigA,isigG,2,1),Curr34(1,
     -                    iflavF3,iisig(3,isigF),2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,
     -                    isigG,2,1)) !iinsA=2 iinsG=1
                        CALL VVSXXX(Curr12A(1,iflavF1,iisig(1,isigF),2,
     -                    isigA,2),Curr34GO(1,iflavF3,iisig(3,isigF),3,
     -                    isigG,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,
     -                    isigG,2,3)) !iinsA=2 iinsG=3
                        CALL VVSXXX(Curr12A(1,iflavF1,iisig(1,isigF),2,
     -                    isigA,2),Curr34GO(1,iflavF3,iisig(3,isigF),4,
     -                    isigG,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,
     -                    isigG,2,4)) !iinsA=2 iinsG=4
                        
                        CALL VVSXXX(Curr12GO(1,iflavF1,iisig(1,isigF),1,
     -                    isigG,2),Curr34A(1,iflavF3,iisig(3,isigF),3,
     -                    isigA,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,
     -                    isigG,3,1)) !iinsA=3 iinsG=1
                        CALL VVSXXX(Curr12GO(1,iflavF1,iisig(1,isigF),2,
     -                    isigG,2),Curr34A(1,iflavF3,iisig(3,isigF),3,
     -                    isigA,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,
     -                    isigG,3,2)) !iinsA=3 iinsG=2
                        CALL VVSXXX(Curr12(1,iflavF1,iisig(1,isigF),2),
     -                    Curr34GOt1Z(1,iflavF3,iisig(3,isigF),isigA,
     -                    isigG,3,4),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,
     -                    isigG,3,4)) !iinsA=3 iinsG=4
                        CALL VVSXXX(Curr12GO(1,iflavF1,iisig(1,isigF),1,
     -                    isigG,2),Curr34A(1,iflavF3,iisig(3,isigF),4,
     -                    isigA,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,
     -                    isigG,4,1)) !iinsA=4 iinsG=1
                        CALL VVSXXX(Curr12GO(1,iflavF1,iisig(1,isigF),2,
     -                    isigG,2),Curr34A(1,iflavF3,iisig(3,isigF),4,
     -                    isigA,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,
     -                    isigG,4,2)) !iinsA=4 iinsG=2
                        CALL VVSXXX(Curr12(1,iflavF1,iisig(1,isigF),2),
     -                    Curr34GOt1Z(1,iflavF3,iisig(3,isigF),isigA,
     -                    isigG,4,3),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,
     -                    isigG,4,3)) !iinsA=4 iinsG=3
                        do iinsA=1,2 !contributions with iinsA=iinsG
                           do iorder=1,2
                              CALL VVSXXX(Curr12GOt2Z(1,iflavF1,iisig(1,
     -                    isigF),isigA,isigG,iinsA,iorder),
     $                             Curr34(1,iflavF3,iisig(3,isigF),2),
     -                    Higgs,GZZH,AmpBuf(iorder))
                           enddo
                           SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,
     -                    iinsA,iinsA)=AmpBuf(1)+AmpBuf(2)
                        enddo
                        do iinsA=3,4 !contributions with iinsA=iinsG
                           do iorder=1,2
                              CALL VVSXXX(Curr12(1,iflavF1,iisig(1,isigF),2),
     $                             Curr34GOt2Z(1,iflavF3,iisig(3,isigF),
     -                    isigA,isigG,iinsA,iorder),Higgs,GZZH,AmpBuf(iorder))
                              
                           enddo
                           SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,
     -                    iinsA,iinsA)=AmpBuf(1)+AmpBuf(2)
                        enddo
                        
                     enddo
                  enddo

!W contributions
                  iflavF3=7-iflavF1
                  isigF=1       ! = (-)(-)
                  CALL VVSXXX(Curr12GOt1W(1,iflavF1,iisig(1,isigF),isigA,isigG,1,2),Curr34(1,iflavF3,iisig(3,isigF),3),Higgs,GWWH,
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,1,2)) !iinsA=1 iinsG=2
                  CALL VVSXXX(Curr12A(1,iflavF1,iisig(1,isigF),1,isigA,3),Curr34GO(1,iflavF3,iisig(3,isigF),3,isigG,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,1,3)) !iinsA=1 iinsG=3
                  CALL VVSXXX(Curr12A(1,iflavF1,iisig(1,isigF),1,isigA,3),Curr34GO(1,iflavF3,iisig(3,isigF),4,isigG,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,1,4)) !iinsA=1 iinsG=4
                  CALL VVSXXX(Curr12GOt1W(1,iflavF1,iisig(1,isigF),isigA,isigG,2,1),Curr34(1,iflavF3,iisig(3,isigF),3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,2,1)) !iinsA=2 iinsG=1
                  CALL VVSXXX(Curr12A(1,iflavF1,iisig(1,isigF),2,isigA,3),Curr34GO(1,iflavF3,iisig(3,isigF),3,isigG,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,2,3)) !iinsA=2 iinsG=3
                  CALL VVSXXX(Curr12A(1,iflavF1,iisig(1,isigF),2,isigA,3),Curr34GO(1,iflavF3,iisig(3,isigF),4,isigG,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,2,4)) !iinsA=2 iinsG=4
                  
                  CALL VVSXXX(Curr12GO(1,iflavF1,iisig(1,isigF),1,isigG,3),Curr34A(1,iflavF3,iisig(3,isigF),3,isigA,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,3,1)) !iinsA=3 iinsG=1
                  CALL VVSXXX(Curr12GO(1,iflavF1,iisig(1,isigF),2,isigG,3),Curr34A(1,iflavF3,iisig(3,isigF),3,isigA,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,3,2)) !iinsA=3 iinsG=2
                  CALL VVSXXX(Curr12(1,iflavF1,iisig(1,isigF),3),Curr34GOt1W(1,iflavF3,iisig(3,isigF),isigA,isigG,3,4),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,3,4)) !iinsA=3 iinsG=4
                  CALL VVSXXX(Curr12GO(1,iflavF1,iisig(1,isigF),1,isigG,3),Curr34A(1,iflavF3,iisig(3,isigF),4,isigA,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,4,1)) !iinsA=4 iinsG=1
                  CALL VVSXXX(Curr12GO(1,iflavF1,iisig(1,isigF),2,isigG,3),Curr34A(1,iflavF3,iisig(3,isigF),4,isigA,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,4,2)) !iinsA=4 iinsG=2
                  CALL VVSXXX(Curr12(1,iflavF1,iisig(1,isigF),3),Curr34GOt1W(1,iflavF3,iisig(3,isigF),isigA,isigG,4,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,4,3)) !iinsA=4 iinsG=3
                  do iinsA=1,2   !contributions with iinsA=iinsG
                     do iorder=1,2
                        CALL VVSXXX(Curr12GOt2W(1,iflavF1,iisig(1,isigF),isigA,isigG,iinsA,iorder),
     $                       Curr34(1,iflavF3,iisig(3,isigF),3),Higgs,GWWH,AmpBuf(iorder))
                     enddo      !SubAmpW(3,4,isigF,isigA,isigG,1,1)=0?
                     SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,iinsA)=AmpBuf(1)+AmpBuf(2)
                  enddo
                  do iinsA=3,4   !contributions with iinsA=iinsG
                     do iorder=1,2
                        CALL VVSXXX(Curr12(1,iflavF1,iisig(1,isigF),3),
     $                       Curr34GOt2W(1,iflavF3,iisig(3,isigF),isigA,isigG,iinsA,iorder),Higgs,GWWH,AmpBuf(iorder))

                     enddo
                     SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,iinsA)=AmpBuf(1)+AmpBuf(2)
                  enddo
!TBV contributions
                  CALL VVSXXX(Curr12GO(1,iflavF1,iisig(1,isigF),1,isigG,3),Curr34tbv(1,iflavF3,iisig(3,isigF),isigA),Higgs,GWWH,
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,1,isigG))
                  CALL VVSXXX(Curr12GO(1,iflavF1,iisig(1,isigF),2,isigG,3),Curr34tbv(1,iflavF3,iisig(3,isigF),isigA),Higgs,GWWH,
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,2,isigG))
                  CALL VVSXXX(Curr12tbv(1,iflavF1,iisig(1,isigF),isigA),Curr34GO(1,iflavF3,iisig(3,isigF),3,isigG,3),Higgs,GWWH,
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,3,isigG))
                  CALL VVSXXX(Curr12tbv(1,iflavF1,iisig(1,isigF),isigA),Curr34GO(1,iflavF3,iisig(3,isigF),4,isigG,3),Higgs,GWWH,
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,4,isigG))

                  CALL VVSXXX(Curr12(1,iflavF1,iisig(1,isigF),3),Curr34GOtbv(1,iflavF3,iisig(3,isigF),isigA,3,isigG),Higgs,GWWH,
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,3,isigG))
                  CALL VVSXXX(Curr12(1,iflavF1,iisig(1,isigF),3),Curr34GOtbv(1,iflavF3,iisig(3,isigF),isigA,4,isigG),Higgs,GWWH,
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,4,isigG))
                  CALL VVSXXX(Curr12GOtbv(1,iflavF1,iisig(1,isigF),isigA,1,isigG),Curr34(1,iflavF3,iisig(3,isigF),3),Higgs,GWWH,
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,1,isigG))
                  CALL VVSXXX(Curr12GOtbv(1,iflavF1,iisig(1,isigF),isigA,2,isigG),Curr34(1,iflavF3,iisig(3,isigF),3),Higgs,GWWH,
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,2,isigG))

                  AmpSumWW(iflavF1,iflavF3,isigF,isigA,isigG,1)=            !modifications with gluons on the 12-line
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,1,isigG)+
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,2,isigG)+
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,1,isigG)+
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,2,isigG)
                  AmpSumWW(iflavF1,iflavF3,isigF,isigA,isigG,2)=            !modifications with gluons on the 34-line
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,3,isigG)+
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,4,isigG)+
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,3,isigG)+
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,4,isigG)

                  do iinsA=1,4
                     AmpSumWW(iflavF1,iflavF3,isigF,isigA,isigG,1) =
     -                    AmpSumWW(iflavF1,iflavF3,isigF,isigA,isigG,1)+
     $                    SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,1)+
     $                    SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,2)
                     AmpSumWW(iflavF1,iflavF3,isigF,isigA,isigG,2)=
     -                    AmpSumWW(iflavF1,iflavF3,isigF,isigA,isigG,2)+
     $                    SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,3)+
     $                    SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,4) 
                    
c$$$                     if(iflavF1.eq.3.and.iflavF3.eq.4)then
c$$$                        print*,"SubAmpW(",3,4,isigF,isigA,isigG,iinsA,1,")",SubAmpW(3,4,isigF,isigA,isigG,iinsA,1)
c$$$                        print*,"SubAmpW(",3,4,isigF,isigA,isigG,iinsA,2,")",SubAmpW(3,4,isigF,isigA,isigG,iinsA,2)
c$$$                        print*,"SubAmpW(",3,4,isigF,isigA,isigG,iinsA,3,")",SubAmpW(3,4,isigF,isigA,isigG,iinsA,3)
c$$$                        print*,"SubAmpW(",3,4,isigF,isigA,isigG,iinsA,4,")",SubAmpW(3,4,isigF,isigA,isigG,iinsA,4)
c$$$                     endif
                  enddo
               enddo
            enddo
         enddo
         
         do isigF=1,4
            do isigA=-1,1,2
               do isigG=-1,1,2
                  do k=1,2
                     res(5,k)=res(5,k)+
     -                    dreal(AmpSumWW(3,4,isigF,isigA,isigG,k))**2+
     -                    dimag(AmpSumWW(3,4,isigF,isigA,isigG,k))**2
                     res(6,k)=res(6,k)+
     -                    dreal(AmpSumWW(4,3,isigF,isigA,isigG,k))**2+
     -                    dimag(AmpSumWW(4,3,isigF,isigA,isigG,k))**2
                  enddo
               enddo
            enddo
         enddo
**********************************************************************************************************************************************************
      elseif(Isign(1).eq.1.and.
c     $     Isign(2).eq.1.and.
c     $     Isign(3).eq.-1.and.
     $     Isign(4).eq.-1)then
**********************************************************************************************************************************************************

c     now combine two currents with the higgs emission:
c     Curr12(43)A/GO cab be combined with Curr43(12)GO/A
c     Curr12(43)GOt1_ can be combined with Curr43(12) and
c     Curr12(43)GOt2_ can be combined with Curr43(12)
c     The result is
c     SubAmp_(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,iinsG)
         do isigA=-1,1,2
            do isigG=-1,1,2
               do iflavF1=3,4
                  do iflavF3=3,4
                  !Z contributions
                     do isigF=1,4
                        CALL VVSXXX(Curr12GOt1Z(1,iflavF1,
     -                    iisig(1,isigF),isigA,isigG,1,2),
     -                    Curr43(1,iflavF3,iisig(3,isigF),2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,
     -                    isigG,1,2)) !iinsA=1 iinsG=2
                        CALL VVSXXX(Curr12A(1,iflavF1,iisig(1,isigF),1,
     -                    isigA,2),Curr43GO(1,iflavF3,iisig(3,isigF),3,
     -                    isigG,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,
     -                    isigG,1,3)) !iinsA=1 iinsG=3
                        CALL VVSXXX(Curr12A(1,iflavF1,iisig(1,isigF),1,
     -                    isigA,2),Curr43GO(1,iflavF3,iisig(3,isigF),4,
     -                    isigG,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,
     -                    isigG,1,4)) !iinsA=1 iinsG=4
                        CALL VVSXXX(Curr12GOt1Z(1,iflavF1,iisig(1,isigF),
     -                    isigA,isigG,2,1),Curr43(1,iflavF3,
     -                    iisig(3,isigF),2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,
     -                    isigG,2,1)) !iinsA=2 iinsG=1
                        CALL VVSXXX(Curr12A(1,iflavF1,iisig(1,isigF),2,
     -                    isigA,2),Curr43GO(1,iflavF3,iisig(3,isigF),3,
     -                    isigG,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,
     -                    isigG,2,3)) !iinsA=2 iinsG=3
                        CALL VVSXXX(Curr12A(1,iflavF1,iisig(1,isigF),2,
     -                    isigA,2),Curr43GO(1,iflavF3,iisig(3,isigF),4,
     -                    isigG,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,
     -                    isigG,2,4)) !iinsA=2 iinsG=4
                        
                        CALL VVSXXX(Curr12GO(1,iflavF1,iisig(1,isigF),1,
     -                    isigG,2),Curr43A(1,iflavF3,iisig(3,isigF),3,
     -                    isigA,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,
     -                    isigG,3,1)) !iinsA=3 iinsG=1
                        CALL VVSXXX(Curr12GO(1,iflavF1,iisig(1,isigF),2,
     -                    isigG,2),Curr43A(1,iflavF3,iisig(3,isigF),3,
     -                    isigA,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,
     -                    isigG,3,2)) !iinsA=3 iinsG=2
                        CALL VVSXXX(Curr12(1,iflavF1,iisig(1,isigF),2),
     -                    Curr43GOt1Z(1,iflavF3,iisig(3,isigF),isigA,
     -                    isigG,3,4),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,
     -                    isigG,3,4)) !iinsA=3 iinsG=4
                        CALL VVSXXX(Curr12GO(1,iflavF1,iisig(1,isigF),1,
     -                    isigG,2),Curr43A(1,iflavF3,iisig(3,isigF),4,
     -                    isigA,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,
     -                    isigG,4,1)) !iinsA=4 iinsG=1
                        CALL VVSXXX(Curr12GO(1,iflavF1,iisig(1,isigF),2,
     -                    isigG,2),Curr43A(1,iflavF3,iisig(3,isigF),4,
     -                    isigA,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,
     -                    isigG,4,2)) !iinsA=4 iinsG=2
                        CALL VVSXXX(Curr12(1,iflavF1,iisig(1,isigF),2),
     -                    Curr43GOt1Z(1,iflavF3,iisig(3,isigF),isigA,
     -                    isigG,4,3),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,
     -                    isigG,4,3)) !iinsA=4 iinsG=3
                        do iinsA=1,2 !contributions with iinsA=iinsG
                           do iorder=1,2
                              CALL VVSXXX(Curr12GOt2Z(1,iflavF1,
     -                    iisig(1,isigF),isigA,isigG,iinsA,iorder),
     $                             Curr43(1,iflavF3,iisig(3,isigF),2),
     -                    Higgs,GZZH,AmpBuf(iorder))
                           enddo
                           SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,
     -                    iinsA,iinsA)=AmpBuf(1)+AmpBuf(2)
                        enddo
                        do iinsA=3,4 !contributions with iinsA=iinsG
                           do iorder=1,2
                              CALL VVSXXX(Curr12(1,iflavF1,iisig(1,isigF),2),
     $                             Curr43GOt2Z(1,iflavF3,iisig(3,isigF),
     -                    isigA,isigG,iinsA,iorder),Higgs,GZZH,AmpBuf(iorder))
                              
                           enddo
                           SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,
     -                    iinsA,iinsA)=AmpBuf(1)+AmpBuf(2)
                        enddo
                        
                     enddo
                  enddo

!W contributions
                  iflavF3=iflavF1
                  isigF=2       ! = (-)(-)
                  CALL VVSXXX(Curr12GOt1W(1,iflavF1,iisig(1,isigF),
     -                    isigA,isigG,1,2),Curr43(1,iflavF3,iisig(3,
     -                    isigF),3),Higgs,GWWH,
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,1,2)) !iinsA=1 iinsG=2
                  CALL VVSXXX(Curr12A(1,iflavF1,iisig(1,isigF),1,
     -                    isigA,3),Curr43GO(1,iflavF3,iisig(3,isigF),3,
     -                    isigG,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,1,3)) !iinsA=1 iinsG=3
                  CALL VVSXXX(Curr12A(1,iflavF1,iisig(1,isigF),1,
     -                    isigA,3),Curr43GO(1,iflavF3,iisig(3,isigF),4,
     -                    isigG,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,1,4)) !iinsA=1 iinsG=4
                  CALL VVSXXX(Curr12GOt1W(1,iflavF1,iisig(1,isigF),
     -                    isigA,isigG,2,1),Curr43(1,iflavF3,
     -                    iisig(3,isigF),3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,2,1)) !iinsA=2 iinsG=1
                  CALL VVSXXX(Curr12A(1,iflavF1,iisig(1,isigF),2,isigA,
     -                    3),Curr43GO(1,iflavF3,iisig(3,isigF),3,isigG,
     -                    3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,2,3)) !iinsA=2 iinsG=3
                  CALL VVSXXX(Curr12A(1,iflavF1,iisig(1,isigF),2,isigA,
     -                    3),Curr43GO(1,iflavF3,iisig(3,isigF),4,isigG,
     -                    3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,2,4)) !iinsA=2 iinsG=4
                  
                  CALL VVSXXX(Curr12GO(1,iflavF1,iisig(1,isigF),1,isigG,
     -                    3),Curr43A(1,iflavF3,iisig(3,isigF),3,isigA,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,3,1)) !iinsA=3 iinsG=1
                  CALL VVSXXX(Curr12GO(1,iflavF1,iisig(1,isigF),2,isigG,
     -                    3),Curr43A(1,iflavF3,iisig(3,isigF),3,isigA,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,3,2)) !iinsA=3 iinsG=2
                  CALL VVSXXX(Curr12(1,iflavF1,iisig(1,isigF),3),
     -                    Curr43GOt1W(1,iflavF3,iisig(3,isigF),isigA,
     -                    isigG,3,4),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,3,4)) !iinsA=3 iinsG=4
                  CALL VVSXXX(Curr12GO(1,iflavF1,iisig(1,isigF),1,isigG,
     -                    3),Curr43A(1,iflavF3,iisig(3,isigF),4,isigA,
     -                    3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,4,1)) !iinsA=4 iinsG=1
                  CALL VVSXXX(Curr12GO(1,iflavF1,iisig(1,isigF),2,isigG,
     -                    3),Curr43A(1,iflavF3,iisig(3,isigF),4,isigA,
     -                    3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,4,2)) !iinsA=4 iinsG=2
                  CALL VVSXXX(Curr12(1,iflavF1,iisig(1,isigF),3),
     -                    Curr43GOt1W(1,iflavF3,iisig(3,isigF),isigA,
     -                    isigG,4,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,4,3)) !iinsA=4 iinsG=3
                  do iinsA=1,2   !contributions with iinsA=iinsG
                     do iorder=1,2
                        CALL VVSXXX(Curr12GOt2W(1,iflavF1,
     -                    iisig(1,isigF),isigA,isigG,iinsA,iorder),
     $                       Curr43(1,iflavF3,iisig(3,isigF),3),Higgs,
     -                    GWWH,AmpBuf(iorder))
                     enddo      !SubAmpW(3,4,isigF,isigA,isigG,1,1)=0?
                     SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,
     -                    iinsA)=AmpBuf(1)+AmpBuf(2)
                  enddo
                  do iinsA=3,4   !contributions with iinsA=iinsG
                     do iorder=1,2
                        CALL VVSXXX(Curr12(1,iflavF1,iisig(1,isigF),3),
     $                       Curr43GOt2W(1,iflavF3,iisig(3,isigF),isigA,
     -                    isigG,iinsA,iorder),Higgs,GWWH,AmpBuf(iorder))

                     enddo
                     SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,
     -                    iinsA)=AmpBuf(1)+AmpBuf(2)
                  enddo
!TBV contributions
                  CALL VVSXXX(Curr12GO(1,iflavF1,iisig(1,isigF),1,isigG,
     -                    3),Curr43tbv(1,iflavF3,iisig(3,isigF),isigA),
     -                    Higgs,GWWH,
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,1,isigG))
                  CALL VVSXXX(Curr12GO(1,iflavF1,iisig(1,isigF),2,isigG,
     -                    3),Curr43tbv(1,iflavF3,iisig(3,isigF),isigA),
     -                    Higgs,GWWH,
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,2,isigG))
                  CALL VVSXXX(Curr12tbv(1,iflavF1,iisig(1,isigF),isigA),
     -                    Curr43GO(1,iflavF3,iisig(3,isigF),3,isigG,3),
     -                    Higgs,GWWH,
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,3,isigG))
                  CALL VVSXXX(Curr12tbv(1,iflavF1,iisig(1,isigF),isigA),
     -                    Curr43GO(1,iflavF3,iisig(3,isigF),4,isigG,3),
     -                    Higgs,GWWH,
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,4,isigG))

                  CALL VVSXXX(Curr12(1,iflavF1,iisig(1,isigF),3),
     -                    Curr43GOtbv(1,iflavF3,iisig(3,isigF),isigA,3,
     -                    isigG),Higgs,GWWH,
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,3,isigG))
                  CALL VVSXXX(Curr12(1,iflavF1,iisig(1,isigF),3),
     -                    Curr43GOtbv(1,iflavF3,iisig(3,isigF),isigA,4,
     -                    isigG),Higgs,GWWH,
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,4,isigG))
                  CALL VVSXXX(Curr12GOtbv(1,iflavF1,iisig(1,isigF),isigA,1,isigG),Curr43(1,iflavF3,iisig(3,isigF),3),Higgs,GWWH,
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,1,isigG))
                  CALL VVSXXX(Curr12GOtbv(1,iflavF1,iisig(1,isigF),isigA,2,isigG),Curr43(1,iflavF3,iisig(3,isigF),3),Higgs,GWWH,
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,2,isigG))

                  AmpSumWW(iflavF1,iflavF3,isigF,isigA,isigG,1)=            !modifications with gluons on the 12-line
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,1,isigG)+
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,2,isigG)+
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,1,isigG)+
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,2,isigG)
                  AmpSumWW(iflavF1,iflavF3,isigF,isigA,isigG,2)=            !modifications with gluons on the 34-line
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,3,isigG)+
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,4,isigG)+
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,3,isigG)+
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,4,isigG)

                  do iinsA=1,4
                     AmpSumWW(iflavF1,iflavF3,isigF,isigA,isigG,1)=
     -                    AmpSumWW(iflavF1,iflavF3,isigF,isigA,isigG,1)+
     $                    SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,1)+
     $                    SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,2)
                     AmpSumWW(iflavF1,iflavF3,isigF,isigA,isigG,2)=
     -                    AmpSumWW(iflavF1,iflavF3,isigF,isigA,isigG,2)+
     $                    SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,3)+
     $                    SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,4)                     
                  enddo
               enddo
            enddo
         enddo
         
         do isigF=1,4
            do isigA=-1,1,2
               do isigG=-1,1,2
                  do k=1,2
                     res(5,k)=res(5,k)+dreal(AmpSumWW(3,3,isigF,isigA,isigG,k))**2+dimag(AmpSumWW(3,3,isigF,isigA,isigG,k))**2
                     res(6,k)=res(6,k)+dreal(AmpSumWW(4,4,isigF,isigA,isigG,k))**2+dimag(AmpSumWW(4,4,isigF,isigA,isigG,k))**2
                  enddo
               enddo
            enddo
         enddo

**********************************************************************************************************************************************************
      elseif(Isign(1).eq.-1.and.
c     $     Isign(2).eq.-1.and.
c     $     Isign(3).eq.1.and.
     $     Isign(4).eq.1)then
**********************************************************************************************************************************************************

c     now combine two currents with the higgs emission:
c     Curr21(34)A/GO cab be combined with Curr34(21)GO/A
c     Curr21(34)GOt1_ can be combined with Curr34(21) and
c     Curr21(34)GOt2_ can be combined with Curr34(21)
c     The result is
c     SubAmp_(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,iinsG)
         do isigA=-1,1,2
            do isigG=-1,1,2
               do iflavF1=3,4
                  do iflavF3=3,4
                  !Z contributions
                     do isigF=1,4
                        CALL VVSXXX(Curr21GOt1Z(1,iflavF1,iisig(1,isigF),isigA,isigG,1,2),Curr34(1,iflavF3,iisig(3,isigF),2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,1,2)) !iinsA=1 iinsG=2
                        CALL VVSXXX(Curr21A(1,iflavF1,iisig(1,isigF),1,isigA,2),Curr34GO(1,iflavF3,iisig(3,isigF),3,isigG,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,1,3)) !iinsA=1 iinsG=3
                        CALL VVSXXX(Curr21A(1,iflavF1,iisig(1,isigF),1,isigA,2),Curr34GO(1,iflavF3,iisig(3,isigF),4,isigG,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,1,4)) !iinsA=1 iinsG=4
                        CALL VVSXXX(Curr21GOt1Z(1,iflavF1,iisig(1,isigF),isigA,isigG,2,1),Curr34(1,iflavF3,iisig(3,isigF),2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,2,1)) !iinsA=2 iinsG=1
                        CALL VVSXXX(Curr21A(1,iflavF1,iisig(1,isigF),2,isigA,2),Curr34GO(1,iflavF3,iisig(3,isigF),3,isigG,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,2,3)) !iinsA=2 iinsG=3
                        CALL VVSXXX(Curr21A(1,iflavF1,iisig(1,isigF),2,isigA,2),Curr34GO(1,iflavF3,iisig(3,isigF),4,isigG,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,2,4)) !iinsA=2 iinsG=4
                        
                        CALL VVSXXX(Curr21GO(1,iflavF1,iisig(1,isigF),1,isigG,2),Curr34A(1,iflavF3,iisig(3,isigF),3,isigA,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,3,1)) !iinsA=3 iinsG=1
                        CALL VVSXXX(Curr21GO(1,iflavF1,iisig(1,isigF),2,isigG,2),Curr34A(1,iflavF3,iisig(3,isigF),3,isigA,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,3,2)) !iinsA=3 iinsG=2
                        CALL VVSXXX(Curr21(1,iflavF1,iisig(1,isigF),2),Curr34GOt1Z(1,iflavF3,iisig(3,isigF),isigA,isigG,3,4),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,3,4)) !iinsA=3 iinsG=4
                        CALL VVSXXX(Curr21GO(1,iflavF1,iisig(1,isigF),1,isigG,2),Curr34A(1,iflavF3,iisig(3,isigF),4,isigA,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,4,1)) !iinsA=4 iinsG=1
                        CALL VVSXXX(Curr21GO(1,iflavF1,iisig(1,isigF),2,isigG,2),Curr34A(1,iflavF3,iisig(3,isigF),4,isigA,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,4,2)) !iinsA=4 iinsG=2
                        CALL VVSXXX(Curr21(1,iflavF1,iisig(1,isigF),2),Curr34GOt1Z(1,iflavF3,iisig(3,isigF),isigA,isigG,4,3),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,4,3)) !iinsA=4 iinsG=3
                        do iinsA=1,2 !contributions with iinsA=iinsG
                           do iorder=1,2
                              CALL VVSXXX(Curr21GOt2Z(1,iflavF1,iisig(1,isigF),isigA,isigG,iinsA,iorder),
     $                             Curr34(1,iflavF3,iisig(3,isigF),2),Higgs,GZZH,AmpBuf(iorder))
                           enddo
                           SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,iinsA)=AmpBuf(1)+AmpBuf(2)
                        enddo
                        do iinsA=3,4 !contributions with iinsA=iinsG
                           do iorder=1,2
                              CALL VVSXXX(Curr21(1,iflavF1,iisig(1,isigF),2),
     $                             Curr34GOt2Z(1,iflavF3,iisig(3,isigF),isigA,isigG,iinsA,iorder),Higgs,GZZH,AmpBuf(iorder))
                              
                           enddo
                           SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,iinsA)=AmpBuf(1)+AmpBuf(2)
                        enddo
                        
                     enddo
                  enddo

!W contributions
                  iflavF3=iflavF1
                  isigF=3       ! = (-)(-)
                  CALL VVSXXX(Curr21GOt1W(1,iflavF1,iisig(1,isigF),isigA,isigG,1,2),Curr34(1,iflavF3,iisig(3,isigF),3),Higgs,GWWH,
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,1,2)) !iinsA=1 iinsG=2
                  CALL VVSXXX(Curr21A(1,iflavF1,iisig(1,isigF),1,isigA,3),Curr34GO(1,iflavF3,iisig(3,isigF),3,isigG,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,1,3)) !iinsA=1 iinsG=3
                  CALL VVSXXX(Curr21A(1,iflavF1,iisig(1,isigF),1,isigA,3),Curr34GO(1,iflavF3,iisig(3,isigF),4,isigG,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,1,4)) !iinsA=1 iinsG=4
                  CALL VVSXXX(Curr21GOt1W(1,iflavF1,iisig(1,isigF),isigA,isigG,2,1),Curr34(1,iflavF3,iisig(3,isigF),3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,2,1)) !iinsA=2 iinsG=1
                  CALL VVSXXX(Curr21A(1,iflavF1,iisig(1,isigF),2,isigA,3),Curr34GO(1,iflavF3,iisig(3,isigF),3,isigG,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,2,3)) !iinsA=2 iinsG=3
                  CALL VVSXXX(Curr21A(1,iflavF1,iisig(1,isigF),2,isigA,3),Curr34GO(1,iflavF3,iisig(3,isigF),4,isigG,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,2,4)) !iinsA=2 iinsG=4
                  
                  CALL VVSXXX(Curr21GO(1,iflavF1,iisig(1,isigF),1,isigG,3),Curr34A(1,iflavF3,iisig(3,isigF),3,isigA,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,3,1)) !iinsA=3 iinsG=1
                  CALL VVSXXX(Curr21GO(1,iflavF1,iisig(1,isigF),2,isigG,3),Curr34A(1,iflavF3,iisig(3,isigF),3,isigA,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,3,2)) !iinsA=3 iinsG=2
                  CALL VVSXXX(Curr21(1,iflavF1,iisig(1,isigF),3),Curr34GOt1W(1,iflavF3,iisig(3,isigF),isigA,isigG,3,4),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,3,4)) !iinsA=3 iinsG=4
                  CALL VVSXXX(Curr21GO(1,iflavF1,iisig(1,isigF),1,isigG,3),Curr34A(1,iflavF3,iisig(3,isigF),4,isigA,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,4,1)) !iinsA=4 iinsG=1
                  CALL VVSXXX(Curr21GO(1,iflavF1,iisig(1,isigF),2,isigG,3),Curr34A(1,iflavF3,iisig(3,isigF),4,isigA,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,4,2)) !iinsA=4 iinsG=2
                  CALL VVSXXX(Curr21(1,iflavF1,iisig(1,isigF),3),Curr34GOt1W(1,iflavF3,iisig(3,isigF),isigA,isigG,4,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,4,3)) !iinsA=4 iinsG=3
                  do iinsA=1,2   !contributions with iinsA=iinsG
                     do iorder=1,2
                        CALL VVSXXX(Curr21GOt2W(1,iflavF1,iisig(1,isigF),isigA,isigG,iinsA,iorder),
     $                       Curr34(1,iflavF3,iisig(3,isigF),3),Higgs,GWWH,AmpBuf(iorder))
                     enddo      !SubAmpW(3,4,isigF,isigA,isigG,1,1)=0?
                     SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,iinsA)=AmpBuf(1)+AmpBuf(2)
                  enddo
                  do iinsA=3,4   !contributions with iinsA=iinsG
                     do iorder=1,2
                        CALL VVSXXX(Curr21(1,iflavF1,iisig(1,isigF),3),
     $                       Curr34GOt2W(1,iflavF3,iisig(3,isigF),isigA,isigG,iinsA,iorder),Higgs,GWWH,AmpBuf(iorder))

                     enddo
                     SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,iinsA)=AmpBuf(1)+AmpBuf(2)
                  enddo
!TBV contributions
                  CALL VVSXXX(Curr21GO(1,iflavF1,iisig(1,isigF),1,isigG,3),Curr34tbv(1,iflavF3,iisig(3,isigF),isigA),Higgs,GWWH,
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,1,isigG))
                  CALL VVSXXX(Curr21GO(1,iflavF1,iisig(1,isigF),2,isigG,3),Curr34tbv(1,iflavF3,iisig(3,isigF),isigA),Higgs,GWWH,
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,2,isigG))
                  CALL VVSXXX(Curr21tbv(1,iflavF1,iisig(1,isigF),isigA),Curr34GO(1,iflavF3,iisig(3,isigF),3,isigG,3),Higgs,GWWH,
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,3,isigG))
                  CALL VVSXXX(Curr21tbv(1,iflavF1,iisig(1,isigF),isigA),Curr34GO(1,iflavF3,iisig(3,isigF),4,isigG,3),Higgs,GWWH,
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,4,isigG))

                  CALL VVSXXX(Curr21(1,iflavF1,iisig(1,isigF),3),Curr34GOtbv(1,iflavF3,iisig(3,isigF),isigA,3,isigG),Higgs,GWWH,
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,3,isigG))
                  CALL VVSXXX(Curr21(1,iflavF1,iisig(1,isigF),3),Curr34GOtbv(1,iflavF3,iisig(3,isigF),isigA,4,isigG),Higgs,GWWH,
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,4,isigG))
                  CALL VVSXXX(Curr21GOtbv(1,iflavF1,iisig(1,isigF),isigA,1,isigG),Curr34(1,iflavF3,iisig(3,isigF),3),Higgs,GWWH,
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,1,isigG))
                  CALL VVSXXX(Curr21GOtbv(1,iflavF1,iisig(1,isigF),isigA,2,isigG),Curr34(1,iflavF3,iisig(3,isigF),3),Higgs,GWWH,
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,2,isigG))

                  AmpSumWW(iflavF1,iflavF3,isigF,isigA,isigG,1)=            !modifications with gluons on the 12-line
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,1,isigG)+
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,2,isigG)+
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,1,isigG)+
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,2,isigG)
                  AmpSumWW(iflavF1,iflavF3,isigF,isigA,isigG,2)=            !modifications with gluons on the 34-line
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,3,isigG)+
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,4,isigG)+
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,3,isigG)+
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,4,isigG)

                  do iinsA=1,4
                     AmpSumWW(iflavF1,iflavF3,isigF,isigA,isigG,1)=AmpSumWW(iflavF1,iflavF3,isigF,isigA,isigG,1)+
     $                    SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,1)+
     $                    SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,2)
                     AmpSumWW(iflavF1,iflavF3,isigF,isigA,isigG,2)=AmpSumWW(iflavF1,iflavF3,isigF,isigA,isigG,2)+
     $                    SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,3)+
     $                    SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,4)                     
                  enddo
               enddo
            enddo
         enddo
         
         do isigF=1,4
            do isigA=-1,1,2
               do isigG=-1,1,2
                  do k=1,2
                     res(5,k)=res(5,k)+dreal(AmpSumWW(4,4,isigF,isigA,isigG,k))**2+dimag(AmpSumWW(4,4,isigF,isigA,isigG,k))**2
                     res(6,k)=res(6,k)+dreal(AmpSumWW(3,3,isigF,isigA,isigG,k))**2+dimag(AmpSumWW(3,3,isigF,isigA,isigG,k))**2
                  enddo
               enddo
            enddo
         enddo
**********************************************************************************************************************************************************
      elseif(Isign(1).eq.-1.and.
c     $     Isign(2).eq.-1.and.
c     $     Isign(3).eq.-1.and.
     $     Isign(4).eq.-1)then
**********************************************************************************************************************************************************

c     now combine two currents with the higgs emission:
c     Curr21(43)A/GO cab be combined with Curr43(21)GO/A
c     Curr21(43)GOt1_ can be combined with Curr43(21) and
c     Curr21(43)GOt2_ can be combined with Curr43(21)
c     The result is
c     SubAmp_(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,iinsG)
         do isigA=-1,1,2
            do isigG=-1,1,2
               do iflavF1=3,4
                  do iflavF3=3,4
                  !Z contributions
                     do isigF=1,4
                        CALL VVSXXX(Curr21GOt1Z(1,iflavF1,iisig(1,isigF),isigA,isigG,1,2),Curr43(1,iflavF3,iisig(3,isigF),2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,1,2)) !iinsA=1 iinsG=2
                        CALL VVSXXX(Curr21A(1,iflavF1,iisig(1,isigF),1,isigA,2),Curr43GO(1,iflavF3,iisig(3,isigF),3,isigG,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,1,3)) !iinsA=1 iinsG=3
                        CALL VVSXXX(Curr21A(1,iflavF1,iisig(1,isigF),1,isigA,2),Curr43GO(1,iflavF3,iisig(3,isigF),4,isigG,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,1,4)) !iinsA=1 iinsG=4
                        CALL VVSXXX(Curr21GOt1Z(1,iflavF1,iisig(1,isigF),isigA,isigG,2,1),Curr43(1,iflavF3,iisig(3,isigF),2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,2,1)) !iinsA=2 iinsG=1
                        CALL VVSXXX(Curr21A(1,iflavF1,iisig(1,isigF),2,isigA,2),Curr43GO(1,iflavF3,iisig(3,isigF),3,isigG,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,2,3)) !iinsA=2 iinsG=3
                        CALL VVSXXX(Curr21A(1,iflavF1,iisig(1,isigF),2,isigA,2),Curr43GO(1,iflavF3,iisig(3,isigF),4,isigG,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,2,4)) !iinsA=2 iinsG=4
                        
                        CALL VVSXXX(Curr21GO(1,iflavF1,iisig(1,isigF),1,isigG,2),Curr43A(1,iflavF3,iisig(3,isigF),3,isigA,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,3,1)) !iinsA=3 iinsG=1
                        CALL VVSXXX(Curr21GO(1,iflavF1,iisig(1,isigF),2,isigG,2),Curr43A(1,iflavF3,iisig(3,isigF),3,isigA,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,3,2)) !iinsA=3 iinsG=2
                        CALL VVSXXX(Curr21(1,iflavF1,iisig(1,isigF),2),Curr43GOt1Z(1,iflavF3,iisig(3,isigF),isigA,isigG,3,4),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,3,4)) !iinsA=3 iinsG=4
                        CALL VVSXXX(Curr21GO(1,iflavF1,iisig(1,isigF),1,isigG,2),Curr43A(1,iflavF3,iisig(3,isigF),4,isigA,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,4,1)) !iinsA=4 iinsG=1
                        CALL VVSXXX(Curr21GO(1,iflavF1,iisig(1,isigF),2,isigG,2),Curr43A(1,iflavF3,iisig(3,isigF),4,isigA,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,4,2)) !iinsA=4 iinsG=2
                        CALL VVSXXX(Curr21(1,iflavF1,iisig(1,isigF),2),Curr43GOt1Z(1,iflavF3,iisig(3,isigF),isigA,isigG,4,3),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,4,3)) !iinsA=4 iinsG=3
                        do iinsA=1,2 !contributions with iinsA=iinsG
                           do iorder=1,2
                              CALL VVSXXX(Curr21GOt2Z(1,iflavF1,iisig(1,isigF),isigA,isigG,iinsA,iorder),
     $                             Curr43(1,iflavF3,iisig(3,isigF),2),Higgs,GZZH,AmpBuf(iorder))
                           enddo
                           SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,iinsA)=AmpBuf(1)+AmpBuf(2)
                        enddo
                        do iinsA=3,4 !contributions with iinsA=iinsG
                           do iorder=1,2
                              CALL VVSXXX(Curr21(1,iflavF1,iisig(1,isigF),2),
     $                             Curr43GOt2Z(1,iflavF3,iisig(3,isigF),isigA,isigG,iinsA,iorder),Higgs,GZZH,AmpBuf(iorder))
                              
                           enddo
                           SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,iinsA)=AmpBuf(1)+AmpBuf(2)
                        enddo
                        
                     enddo
                  enddo

!W contributions
                  iflavF3=7-iflavF1
                  isigF=4       ! = (-)(-)
                  CALL VVSXXX(Curr21GOt1W(1,iflavF1,iisig(1,isigF),isigA,isigG,1,2),Curr43(1,iflavF3,iisig(3,isigF),3),Higgs,GWWH,
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,1,2)) !iinsA=1 iinsG=2
                  CALL VVSXXX(Curr21A(1,iflavF1,iisig(1,isigF),1,isigA,3),Curr43GO(1,iflavF3,iisig(3,isigF),3,isigG,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,1,3)) !iinsA=1 iinsG=3
                  CALL VVSXXX(Curr21A(1,iflavF1,iisig(1,isigF),1,isigA,3),Curr43GO(1,iflavF3,iisig(3,isigF),4,isigG,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,1,4)) !iinsA=1 iinsG=4
                  CALL VVSXXX(Curr21GOt1W(1,iflavF1,iisig(1,isigF),isigA,isigG,2,1),Curr43(1,iflavF3,iisig(3,isigF),3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,2,1)) !iinsA=2 iinsG=1
                  CALL VVSXXX(Curr21A(1,iflavF1,iisig(1,isigF),2,isigA,3),Curr43GO(1,iflavF3,iisig(3,isigF),3,isigG,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,2,3)) !iinsA=2 iinsG=3
                  CALL VVSXXX(Curr21A(1,iflavF1,iisig(1,isigF),2,isigA,3),Curr43GO(1,iflavF3,iisig(3,isigF),4,isigG,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,2,4)) !iinsA=2 iinsG=4
                  
                  CALL VVSXXX(Curr21GO(1,iflavF1,iisig(1,isigF),1,isigG,3),Curr43A(1,iflavF3,iisig(3,isigF),3,isigA,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,3,1)) !iinsA=3 iinsG=1
                  CALL VVSXXX(Curr21GO(1,iflavF1,iisig(1,isigF),2,isigG,3),Curr43A(1,iflavF3,iisig(3,isigF),3,isigA,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,3,2)) !iinsA=3 iinsG=2
                  CALL VVSXXX(Curr21(1,iflavF1,iisig(1,isigF),3),Curr43GOt1W(1,iflavF3,iisig(3,isigF),isigA,isigG,3,4),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,3,4)) !iinsA=3 iinsG=4
                  CALL VVSXXX(Curr21GO(1,iflavF1,iisig(1,isigF),1,isigG,3),Curr43A(1,iflavF3,iisig(3,isigF),4,isigA,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,4,1)) !iinsA=4 iinsG=1
                  CALL VVSXXX(Curr21GO(1,iflavF1,iisig(1,isigF),2,isigG,3),Curr43A(1,iflavF3,iisig(3,isigF),4,isigA,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,4,2)) !iinsA=4 iinsG=2
                  CALL VVSXXX(Curr21(1,iflavF1,iisig(1,isigF),3),Curr43GOt1W(1,iflavF3,iisig(3,isigF),isigA,isigG,4,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,4,3)) !iinsA=4 iinsG=3
                  do iinsA=1,2   !contributions with iinsA=iinsG
                     do iorder=1,2
                        CALL VVSXXX(Curr21GOt2W(1,iflavF1,iisig(1,isigF),isigA,isigG,iinsA,iorder),
     $                       Curr43(1,iflavF3,iisig(3,isigF),3),Higgs,GWWH,AmpBuf(iorder))
                     enddo      !SubAmpW(3,4,isigF,isigA,isigG,1,1)=0?
                     SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,iinsA)=AmpBuf(1)+AmpBuf(2)
                  enddo
                  do iinsA=3,4   !contributions with iinsA=iinsG
                     do iorder=1,2
                        CALL VVSXXX(Curr21(1,iflavF1,iisig(1,isigF),3),
     $                       Curr43GOt2W(1,iflavF3,iisig(3,isigF),isigA,isigG,iinsA,iorder),Higgs,GWWH,AmpBuf(iorder))

                     enddo
                     SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,iinsA)=AmpBuf(1)+AmpBuf(2)
                  enddo
!TBV contributions
                  CALL VVSXXX(Curr21GO(1,iflavF1,iisig(1,isigF),1,isigG,3),Curr43tbv(1,iflavF3,iisig(3,isigF),isigA),Higgs,GWWH,
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,1,isigG))
                  CALL VVSXXX(Curr21GO(1,iflavF1,iisig(1,isigF),2,isigG,3),Curr43tbv(1,iflavF3,iisig(3,isigF),isigA),Higgs,GWWH,
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,2,isigG))
                  CALL VVSXXX(Curr21tbv(1,iflavF1,iisig(1,isigF),isigA),Curr43GO(1,iflavF3,iisig(3,isigF),3,isigG,3),Higgs,GWWH,
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,3,isigG))
                  CALL VVSXXX(Curr21tbv(1,iflavF1,iisig(1,isigF),isigA),Curr43GO(1,iflavF3,iisig(3,isigF),4,isigG,3),Higgs,GWWH,
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,4,isigG))

                  CALL VVSXXX(Curr21(1,iflavF1,iisig(1,isigF),3),Curr43GOtbv(1,iflavF3,iisig(3,isigF),isigA,3,isigG),Higgs,GWWH,
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,3,isigG))
                  CALL VVSXXX(Curr21(1,iflavF1,iisig(1,isigF),3),Curr43GOtbv(1,iflavF3,iisig(3,isigF),isigA,4,isigG),Higgs,GWWH,
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,4,isigG))
                  CALL VVSXXX(Curr21GOtbv(1,iflavF1,iisig(1,isigF),isigA,1,isigG),Curr43(1,iflavF3,iisig(3,isigF),3),Higgs,GWWH,
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,1,isigG))
                  CALL VVSXXX(Curr21GOtbv(1,iflavF1,iisig(1,isigF),isigA,2,isigG),Curr43(1,iflavF3,iisig(3,isigF),3),Higgs,GWWH,
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,2,isigG))

                  AmpSumWW(iflavF1,iflavF3,isigF,isigA,isigG,1)=            !modifications with gluons on the 12-line
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,1,isigG)+
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,2,isigG)+
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,1,isigG)+
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,2,isigG)
                  AmpSumWW(iflavF1,iflavF3,isigF,isigA,isigG,2)=            !modifications with gluons on the 34-line
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,3,isigG)+
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,4,isigG)+
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,3,isigG)+
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,4,isigG)

                  do iinsA=1,4
                     AmpSumWW(iflavF1,iflavF3,isigF,isigA,isigG,1)=AmpSumWW(iflavF1,iflavF3,isigF,isigA,isigG,1)+
     $                    SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,1)+
     $                    SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,2)
                     AmpSumWW(iflavF1,iflavF3,isigF,isigA,isigG,2)=AmpSumWW(iflavF1,iflavF3,isigF,isigA,isigG,2)+
     $                    SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,3)+
     $                    SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,4)                     
                  enddo
               enddo
            enddo
         enddo
         
         do isigF=1,4
            do isigA=-1,1,2
               do isigG=-1,1,2
                  do k=1,2
                     res(5,k)=res(5,k)+dreal(AmpSumWW(4,3,isigF,isigA,isigG,k))**2+dimag(AmpSumWW(4,3,isigF,isigA,isigG,k))**2
                     res(6,k)=res(6,k)+dreal(AmpSumWW(3,4,isigF,isigA,isigG,k))**2+dimag(AmpSumWW(3,4,isigF,isigA,isigG,k))**2
                  enddo
               enddo
            enddo
         enddo
         endif
         else !=if gsign=-1
**********************************************************************************************************************************************************
      if(
     $           Isign(1).eq.-1.and.
c     $           Isign(2).eq.1.and.
     $           Isign(3).eq.1!.and.
c     $           Isign(4).eq.1
     $           )then
**********************************************************************************************************************************************************

c     now combine two currents with the higgs emission:
c     Curr21(34)A/GI cab be combined with Curr34(21)GI/A
c     Curr21(34)GIt1_ can be combined with Curr34(21) and
c     Curr21(34)GIt2_ can be combined with Curr34(21)
c     The result is
c     SubAmp_(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,iinsG)
         do isigA=-1,1,2
            do isigG=-1,1,2
               do iflavF1=3,4
                  do iflavF3=3,4
                  !Z contributions
                     do isigF=1,4
                        CALL VVSXXX(Curr21GIt1Z(1,iflavF1,iisig(1,isigF),isigA,isigG,1,2),Curr34(1,iflavF3,iisig(3,isigF),2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,1,2)) !iinsA=1 iinsG=2
                        CALL VVSXXX(Curr21GIt1Z(1,iflavF1,iisig(1,isigF),isigA,isigG,2,1),Curr34(1,iflavF3,iisig(3,isigF),2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,2,1)) !iinsA=2 iinsG=1
                        
                        CALL VVSXXX(Curr21GI(1,iflavF1,iisig(1,isigF),1,isigG,2),Curr34A(1,iflavF3,iisig(3,isigF),3,isigA,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,3,1)) !iinsA=3 iinsG=1
                        CALL VVSXXX(Curr21GI(1,iflavF1,iisig(1,isigF),2,isigG,2),Curr34A(1,iflavF3,iisig(3,isigF),3,isigA,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,3,2)) !iinsA=3 iinsG=2
                        CALL VVSXXX(Curr21GI(1,iflavF1,iisig(1,isigF),1,isigG,2),Curr34A(1,iflavF3,iisig(3,isigF),4,isigA,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,4,1)) !iinsA=4 iinsG=1
                        CALL VVSXXX(Curr21GI(1,iflavF1,iisig(1,isigF),2,isigG,2),Curr34A(1,iflavF3,iisig(3,isigF),4,isigA,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,4,2)) !iinsA=4 iinsG=2
                        do iinsA=1,2 !contributions with iinsA=iinsG
                           do iorder=1,2
                              CALL VVSXXX(Curr21GIt2Z(1,iflavF1,iisig(1,isigF),isigA,isigG,iinsA,iorder),
     $                             Curr34(1,iflavF3,iisig(3,isigF),2),Higgs,GZZH,AmpBuf(iorder))
                           enddo
                           SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,iinsA)=AmpBuf(1)+AmpBuf(2)
                        enddo
                        
                     enddo
                  enddo

!W contributions
                  iflavF3=iflavF1
                  isigF=3       ! = (+)(-)
                  CALL VVSXXX(Curr21GIt1W(1,iflavF1,iisig(1,isigF),isigA,isigG,1,2),Curr34(1,iflavF3,iisig(3,isigF),3),Higgs,GWWH,
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,1,2)) !iinsA=1 iinsG=2
                  CALL VVSXXX(Curr21GIt1W(1,iflavF1,iisig(1,isigF),isigA,isigG,2,1),Curr34(1,iflavF3,iisig(3,isigF),3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,2,1)) !iinsA=2 iinsG=1
                  
                  CALL VVSXXX(Curr21GI(1,iflavF1,iisig(1,isigF),1,isigG,3),Curr34A(1,iflavF3,iisig(3,isigF),3,isigA,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,3,1)) !iinsA=3 iinsG=1
                  CALL VVSXXX(Curr21GI(1,iflavF1,iisig(1,isigF),2,isigG,3),Curr34A(1,iflavF3,iisig(3,isigF),3,isigA,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,3,2)) !iinsA=3 iinsG=2
                  CALL VVSXXX(Curr21GI(1,iflavF1,iisig(1,isigF),1,isigG,3),Curr34A(1,iflavF3,iisig(3,isigF),4,isigA,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,4,1)) !iinsA=4 iinsG=1
                  CALL VVSXXX(Curr21GI(1,iflavF1,iisig(1,isigF),2,isigG,3),Curr34A(1,iflavF3,iisig(3,isigF),4,isigA,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,4,2)) !iinsA=4 iinsG=2
                  do iinsA=1,2   !contributions with iinsA=iinsG
                     do iorder=1,2
                        CALL VVSXXX(Curr21GIt2W(1,iflavF1,iisig(1,isigF),isigA,isigG,iinsA,iorder),
     $                       Curr34(1,iflavF3,iisig(3,isigF),3),Higgs,GWWH,AmpBuf(iorder))
                     enddo
                     SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,iinsA)=AmpBuf(1)+AmpBuf(2)
                  enddo
!TBV contributions
                  CALL VVSXXX(Curr21GI(1,iflavF1,iisig(1,isigF),1,isigG,3),Curr34tbv(1,iflavF3,iisig(3,isigF),isigA),Higgs,GWWH,
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,1,isigG))
                  CALL VVSXXX(Curr21GI(1,iflavF1,iisig(1,isigF),2,isigG,3),Curr34tbv(1,iflavF3,iisig(3,isigF),isigA),Higgs,GWWH,
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,2,isigG))
                  CALL VVSXXX(Curr21GItbv(1,iflavF1,iisig(1,isigF),isigA,1,isigG),Curr34(1,iflavF3,iisig(3,isigF),3),Higgs,GWWH,
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,1,isigG))
                  CALL VVSXXX(Curr21GItbv(1,iflavF1,iisig(1,isigF),isigA,2,isigG),Curr34(1,iflavF3,iisig(3,isigF),3),Higgs,GWWH,
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,2,isigG))

                  AmpSumWW(iflavF1,iflavF3,isigF,isigA,isigG,1)=            !modifications with gluons on the 12-line
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,1,isigG)+
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,2,isigG)+
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,1,isigG)+
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,2,isigG)
                  do iinsA=1,4
                     AmpSumWW(iflavF1,iflavF3,isigF,isigA,isigG,1)=AmpSumWW(iflavF1,iflavF3,isigF,isigA,isigG,1)+
     $                    SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,1)+
     $                    SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,2)                 
                  enddo
               enddo
            enddo
         enddo
         
         do isigF=1,4
            do isigA=-1,1,2
               do isigG=-1,1,2
                  do k=1,2
                     res(5,k)=res(5,k)+dreal(AmpSumWW(4,4,isigF,isigA,isigG,k))**2+dimag(AmpSumWW(4,4,isigF,isigA,isigG,k))**2
                     res(6,k)=res(6,k)+dreal(AmpSumWW(3,3,isigF,isigA,isigG,k))**2+dimag(AmpSumWW(3,3,isigF,isigA,isigG,k))**2
                  enddo
               enddo
            enddo
         enddo

****************************************************************************
      elseif(
     $        Isign(1).eq.1.and.
c     $        Isign(2).eq.1.and.
     $        Isign(3).eq.-1!.and.
c     $        Isign(4).eq.1
     $        )then
****************************************************************************

c     now combine two currents with the higgs emission:
c     Curr12(43)A/GO cab be combined with Curr43(12)GO/A
c     Curr12(43)GOt1_ can be combined with Curr43(12) and
c     Curr12(43)GOt2_ can be combined with Curr43(12)
c     The result is
c     SubAmp_(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,iinsG)
         do isigA=-1,1,2
            do isigG=-1,1,2
               do iflavF1=3,4
                  do iflavF3=3,4
                  !Z contributions
                     do isigF=1,4
                        CALL VVSXXX(Curr12A(1,iflavF1,iisig(1,isigF),1,isigA,2),Curr43GI(1,iflavF3,iisig(3,isigF),3,isigG,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,1,3)) !iinsA=1 iinsG=3
                        CALL VVSXXX(Curr12A(1,iflavF1,iisig(1,isigF),1,isigA,2),Curr43GI(1,iflavF3,iisig(3,isigF),4,isigG,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,1,4)) !iinsA=1 iinsG=4
                        CALL VVSXXX(Curr12A(1,iflavF1,iisig(1,isigF),2,isigA,2),Curr43GI(1,iflavF3,iisig(3,isigF),3,isigG,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,2,3)) !iinsA=2 iinsG=3
                        CALL VVSXXX(Curr12A(1,iflavF1,iisig(1,isigF),2,isigA,2),Curr43GI(1,iflavF3,iisig(3,isigF),4,isigG,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,2,4)) !iinsA=2 iinsG=4
                        
                        CALL VVSXXX(Curr12(1,iflavF1,iisig(1,isigF),2),Curr43GIt1Z(1,iflavF3,iisig(3,isigF),isigA,isigG,3,4),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,3,4)) !iinsA=3 iinsG=4
                        CALL VVSXXX(Curr12(1,iflavF1,iisig(1,isigF),2),Curr43GIt1Z(1,iflavF3,iisig(3,isigF),isigA,isigG,4,3),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,4,3)) !iinsA=4 iinsG=3
                        do iinsA=3,4 !contributions with iinsA=iinsG
                           do iorder=1,2
                              CALL VVSXXX(Curr12(1,iflavF1,iisig(1,isigF),2),
     $                             Curr43GIt2Z(1,iflavF3,iisig(3,isigF),isigA,isigG,iinsA,iorder),Higgs,GZZH,AmpBuf(iorder))
                              
                           enddo
                           SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,iinsA)=AmpBuf(1)+AmpBuf(2)
                        enddo
                        
                     enddo
                  enddo

!W contributions
                  iflavF3=iflavF1
                  isigF=2       ! = (-)(+)
                  CALL VVSXXX(Curr12A(1,iflavF1,iisig(1,isigF),1,isigA,3),Curr43GI(1,iflavF3,iisig(3,isigF),3,isigG,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,1,3)) !iinsA=1 iinsG=3
                  CALL VVSXXX(Curr12A(1,iflavF1,iisig(1,isigF),1,isigA,3),Curr43GI(1,iflavF3,iisig(3,isigF),4,isigG,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,1,4)) !iinsA=1 iinsG=4
                  CALL VVSXXX(Curr12A(1,iflavF1,iisig(1,isigF),2,isigA,3),Curr43GI(1,iflavF3,iisig(3,isigF),3,isigG,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,2,3)) !iinsA=2 iinsG=3
                  CALL VVSXXX(Curr12A(1,iflavF1,iisig(1,isigF),2,isigA,3),Curr43GI(1,iflavF3,iisig(3,isigF),4,isigG,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,2,4)) !iinsA=2 iinsG=4
                  
                  CALL VVSXXX(Curr12(1,iflavF1,iisig(1,isigF),3),Curr43GIt1W(1,iflavF3,iisig(3,isigF),isigA,isigG,3,4),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,3,4)) !iinsA=3 iinsG=4
                  CALL VVSXXX(Curr12(1,iflavF1,iisig(1,isigF),3),Curr43GIt1W(1,iflavF3,iisig(3,isigF),isigA,isigG,4,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,4,3)) !iinsA=4 iinsG=3
                  do iinsA=3,4   !contributions with iinsA=iinsG
                     do iorder=1,2
                        CALL VVSXXX(Curr12(1,iflavF1,iisig(1,isigF),3),
     $                       Curr43GIt2W(1,iflavF3,iisig(3,isigF),isigA,isigG,iinsA,iorder),Higgs,GWWH,AmpBuf(iorder))

                     enddo
                     SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,iinsA)=AmpBuf(1)+AmpBuf(2)
                  enddo
!TBV contributions
                  CALL VVSXXX(Curr12tbv(1,iflavF1,iisig(1,isigF),isigA),Curr43GI(1,iflavF3,iisig(3,isigF),3,isigG,3),Higgs,GWWH,
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,3,isigG))
                  CALL VVSXXX(Curr12tbv(1,iflavF1,iisig(1,isigF),isigA),Curr43GI(1,iflavF3,iisig(3,isigF),4,isigG,3),Higgs,GWWH,
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,4,isigG))

                  CALL VVSXXX(Curr12(1,iflavF1,iisig(1,isigF),3),Curr43GItbv(1,iflavF3,iisig(3,isigF),isigA,3,isigG),Higgs,GWWH,
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,3,isigG))
                  CALL VVSXXX(Curr12(1,iflavF1,iisig(1,isigF),3),Curr43GItbv(1,iflavF3,iisig(3,isigF),isigA,4,isigG),Higgs,GWWH,
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,4,isigG))

                  AmpSumWW(iflavF1,iflavF3,isigF,isigA,isigG,2)=            !modifications with gluons on the 34-line
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,3,isigG)+
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,4,isigG)+
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,3,isigG)+
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,4,isigG)

                  do iinsA=1,4
                     AmpSumWW(iflavF1,iflavF3,isigF,isigA,isigG,2)=AmpSumWW(iflavF1,iflavF3,isigF,isigA,isigG,2)+
     $                    SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,3)+
     $                    SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,4)                     
                  enddo
               enddo
            enddo
         enddo
         
         do isigF=1,4
            do isigA=-1,1,2
               do isigG=-1,1,2
                  do k=1,2
                     res(5,k)=res(5,k)+dreal(AmpSumWW(3,3,isigF,isigA,isigG,k))**2+dimag(AmpSumWW(3,3,isigF,isigA,isigG,k))**2
                     res(6,k)=res(6,k)+dreal(AmpSumWW(4,4,isigF,isigA,isigG,k))**2+dimag(AmpSumWW(4,4,isigF,isigA,isigG,k))**2
                  enddo
               enddo
            enddo
         enddo
****************************************************************************
      elseif(
     $        Isign(1).eq.-1.and.
c     $        Isign(2).eq.1.and.
c     $        Isign(3).eq.-1.and.
     $        Isign(4).eq.-1
     $        )then
****************************************************************************

c     now combine two currents with the higgs emission:
c     Curr21(43)A/GI cab be combined with Curr43(21)GI/A
c     Curr21(43)GIt1_ can be combined with Curr43(21) and
c     Curr21(43)GIt2_ can be combined with Curr43(21)
c     The result is
c     SubAmp_(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,iinsG)
         do isigA=-1,1,2
            do isigG=-1,1,2
               do iflavF1=3,4
                  do iflavF3=3,4
                  !Z contributions
                     do isigF=1,4
                        CALL VVSXXX(Curr21GIt1Z(1,iflavF1,iisig(1,isigF),isigA,isigG,1,2),Curr43(1,iflavF3,iisig(3,isigF),2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,1,2)) !iinsA=1 iinsG=2
                        CALL VVSXXX(Curr21GIt1Z(1,iflavF1,iisig(1,isigF),isigA,isigG,2,1),Curr43(1,iflavF3,iisig(3,isigF),2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,2,1)) !iinsA=2 iinsG=1
                        
                        CALL VVSXXX(Curr21GI(1,iflavF1,iisig(1,isigF),1,isigG,2),Curr43A(1,iflavF3,iisig(3,isigF),3,isigA,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,3,1)) !iinsA=3 iinsG=1
                        CALL VVSXXX(Curr21GI(1,iflavF1,iisig(1,isigF),2,isigG,2),Curr43A(1,iflavF3,iisig(3,isigF),3,isigA,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,3,2)) !iinsA=3 iinsG=2
                        CALL VVSXXX(Curr21GI(1,iflavF1,iisig(1,isigF),1,isigG,2),Curr43A(1,iflavF3,iisig(3,isigF),4,isigA,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,4,1)) !iinsA=4 iinsG=1
                        CALL VVSXXX(Curr21GI(1,iflavF1,iisig(1,isigF),2,isigG,2),Curr43A(1,iflavF3,iisig(3,isigF),4,isigA,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,4,2)) !iinsA=4 iinsG=2
                        do iinsA=1,2 !contributions with iinsA=iinsG
                           do iorder=1,2
                              CALL VVSXXX(Curr21GIt2Z(1,iflavF1,iisig(1,isigF),isigA,isigG,iinsA,iorder),
     $                             Curr43(1,iflavF3,iisig(3,isigF),2),Higgs,GZZH,AmpBuf(iorder))
                           enddo
                           SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,iinsA)=AmpBuf(1)+AmpBuf(2)
                        enddo
                        
                     enddo
                  enddo

!W contributions
                  iflavF3=7-iflavF1
                  isigF=4       ! = (-)(-)
                  CALL VVSXXX(Curr21GIt1W(1,iflavF1,iisig(1,isigF),isigA,isigG,1,2),Curr43(1,iflavF3,iisig(3,isigF),3),Higgs,GWWH,
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,1,2)) !iinsA=1 iinsG=2
                  CALL VVSXXX(Curr21GIt1W(1,iflavF1,iisig(1,isigF),isigA,isigG,2,1),Curr43(1,iflavF3,iisig(3,isigF),3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,2,1)) !iinsA=2 iinsG=1
                  
                  CALL VVSXXX(Curr21GI(1,iflavF1,iisig(1,isigF),1,isigG,3),Curr43A(1,iflavF3,iisig(3,isigF),3,isigA,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,3,1)) !iinsA=3 iinsG=1
                  CALL VVSXXX(Curr21GI(1,iflavF1,iisig(1,isigF),2,isigG,3),Curr43A(1,iflavF3,iisig(3,isigF),3,isigA,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,3,2)) !iinsA=3 iinsG=2
                  CALL VVSXXX(Curr21GI(1,iflavF1,iisig(1,isigF),1,isigG,3),Curr43A(1,iflavF3,iisig(3,isigF),4,isigA,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,4,1)) !iinsA=4 iinsG=1
                  CALL VVSXXX(Curr21GI(1,iflavF1,iisig(1,isigF),2,isigG,3),Curr43A(1,iflavF3,iisig(3,isigF),4,isigA,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,4,2)) !iinsA=4 iinsG=2
                  do iinsA=1,2   !contributions with iinsA=iinsG
                     do iorder=1,2
                        CALL VVSXXX(Curr21GIt2W(1,iflavF1,iisig(1,isigF),isigA,isigG,iinsA,iorder),
     $                       Curr43(1,iflavF3,iisig(3,isigF),3),Higgs,GWWH,AmpBuf(iorder))
                     enddo 
                     SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,iinsA)=AmpBuf(1)+AmpBuf(2)
                  enddo
!TBV contributions
                  CALL VVSXXX(Curr21GI(1,iflavF1,iisig(1,isigF),1,isigG,3),Curr43tbv(1,iflavF3,iisig(3,isigF),isigA),Higgs,GWWH,
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,1,isigG))
                  CALL VVSXXX(Curr21GI(1,iflavF1,iisig(1,isigF),2,isigG,3),Curr43tbv(1,iflavF3,iisig(3,isigF),isigA),Higgs,GWWH,
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,2,isigG))

                  CALL VVSXXX(Curr21GItbv(1,iflavF1,iisig(1,isigF),isigA,1,isigG),Curr43(1,iflavF3,iisig(3,isigF),3),Higgs,GWWH,
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,1,isigG))
                  CALL VVSXXX(Curr21GItbv(1,iflavF1,iisig(1,isigF),isigA,2,isigG),Curr43(1,iflavF3,iisig(3,isigF),3),Higgs,GWWH,
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,2,isigG))

                  AmpSumWW(iflavF1,iflavF3,isigF,isigA,isigG,1)=            !modifications with gluons on the 12-line
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,1,isigG)+
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,2,isigG)+
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,1,isigG)+
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,2,isigG)

                  do iinsA=1,4
                     AmpSumWW(iflavF1,iflavF3,isigF,isigA,isigG,1)=AmpSumWW(iflavF1,iflavF3,isigF,isigA,isigG,1)+
     $                    SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,1)+
     $                    SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,2)                    
                  enddo
               enddo
            enddo
         enddo
         
         do isigF=1,4
            do isigA=-1,1,2
               do isigG=-1,1,2
                  do k=1,2
                     res(5,k)=res(5,k)+dreal(AmpSumWW(4,3,isigF,isigA,isigG,k))**2+dimag(AmpSumWW(4,3,isigF,isigA,isigG,k))**2
                     res(6,k)=res(6,k)+dreal(AmpSumWW(3,4,isigF,isigA,isigG,k))**2+dimag(AmpSumWW(3,4,isigF,isigA,isigG,k))**2
                  enddo
               enddo
            enddo
         enddo
****************************************************************************
      elseif(
c     $        Isign(1).eq.-1.and.
     $        Isign(2).eq.-1.and.
     $        Isign(3).eq.-1!.and.
c     $        Isign(4).eq.1
     $        )then
****************************************************************************

c     now combine two currents with the higgs emission:
c     Curr21(43)A/GI cab be combined with Curr43(21)GI/A
c     Curr21(43)GIt1_ can be combined with Curr43(21) and
c     Curr21(43)GIt2_ can be combined with Curr43(21)
c     The result is
c     SubAmp_(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,iinsG)
         do isigA=-1,1,2
            do isigG=-1,1,2
               do iflavF1=3,4
                  do iflavF3=3,4
                  !Z contributions
                     do isigF=1,4
                        CALL VVSXXX(Curr21A(1,iflavF1,iisig(1,isigF),1,isigA,2),Curr43GI(1,iflavF3,iisig(3,isigF),3,isigG,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,1,3)) !iinsA=1 iinsG=3
                        CALL VVSXXX(Curr21A(1,iflavF1,iisig(1,isigF),1,isigA,2),Curr43GI(1,iflavF3,iisig(3,isigF),4,isigG,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,1,4)) !iinsA=1 iinsG=4
                        CALL VVSXXX(Curr21A(1,iflavF1,iisig(1,isigF),2,isigA,2),Curr43GI(1,iflavF3,iisig(3,isigF),3,isigG,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,2,3)) !iinsA=2 iinsG=3
                        CALL VVSXXX(Curr21A(1,iflavF1,iisig(1,isigF),2,isigA,2),Curr43GI(1,iflavF3,iisig(3,isigF),4,isigG,2),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,2,4)) !iinsA=2 iinsG=4
                        
                        CALL VVSXXX(Curr21(1,iflavF1,iisig(1,isigF),2),Curr43GIt1Z(1,iflavF3,iisig(3,isigF),isigA,isigG,3,4),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,3,4)) !iinsA=3 iinsG=4
                        CALL VVSXXX(Curr21(1,iflavF1,iisig(1,isigF),2),Curr43GIt1Z(1,iflavF3,iisig(3,isigF),isigA,isigG,4,3),Higgs,
     $                       GZZH,SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,4,3)) !iinsA=4 iinsG=3

                        do iinsA=3,4 !contributions with iinsA=iinsG
                           do iorder=1,2
                              CALL VVSXXX(Curr21(1,iflavF1,iisig(1,isigF),2),
     $                             Curr43GIt2Z(1,iflavF3,iisig(3,isigF),isigA,isigG,iinsA,iorder),Higgs,GZZH,AmpBuf(iorder))
                              
                           enddo
                           SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,iinsA)=AmpBuf(1)+AmpBuf(2)
                        enddo
                        
                     enddo
                  enddo

!W contributions
                  iflavF3=7-iflavF1
                  isigF=4       ! = (+)(+)
                  CALL VVSXXX(Curr21A(1,iflavF1,iisig(1,isigF),1,isigA,3),Curr43GI(1,iflavF3,iisig(3,isigF),3,isigG,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,1,3)) !iinsA=1 iinsG=3
                  CALL VVSXXX(Curr21A(1,iflavF1,iisig(1,isigF),1,isigA,3),Curr43GI(1,iflavF3,iisig(3,isigF),4,isigG,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,1,4)) !iinsA=1 iinsG=4
                  CALL VVSXXX(Curr21A(1,iflavF1,iisig(1,isigF),2,isigA,3),Curr43GI(1,iflavF3,iisig(3,isigF),3,isigG,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,2,3)) !iinsA=2 iinsG=3
                  CALL VVSXXX(Curr21A(1,iflavF1,iisig(1,isigF),2,isigA,3),Curr43GI(1,iflavF3,iisig(3,isigF),4,isigG,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,2,4)) !iinsA=2 iinsG=4
                  
                  CALL VVSXXX(Curr21(1,iflavF1,iisig(1,isigF),3),Curr43GIt1W(1,iflavF3,iisig(3,isigF),isigA,isigG,3,4),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,3,4)) !iinsA=3 iinsG=4
                  CALL VVSXXX(Curr21(1,iflavF1,iisig(1,isigF),3),Curr43GIt1W(1,iflavF3,iisig(3,isigF),isigA,isigG,4,3),Higgs,GWWH, 
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,4,3)) !iinsA=4 iinsG=3

                  do iinsA=3,4   !contributions with iinsA=iinsG
                     do iorder=1,2
                        CALL VVSXXX(Curr21(1,iflavF1,iisig(1,isigF),3),
     $                       Curr43GIt2W(1,iflavF3,iisig(3,isigF),isigA,isigG,iinsA,iorder),Higgs,GWWH,AmpBuf(iorder))

                     enddo
                     SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,iinsA)=AmpBuf(1)+AmpBuf(2)
                  enddo
!TBV contributions
                  CALL VVSXXX(Curr21tbv(1,iflavF1,iisig(1,isigF),isigA),Curr43GI(1,iflavF3,iisig(3,isigF),3,isigG,3),Higgs,GWWH,
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,3,isigG))
                  CALL VVSXXX(Curr21tbv(1,iflavF1,iisig(1,isigF),isigA),Curr43GI(1,iflavF3,iisig(3,isigF),4,isigG,3),Higgs,GWWH,
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,4,isigG))

                  CALL VVSXXX(Curr21(1,iflavF1,iisig(1,isigF),3),Curr43GItbv(1,iflavF3,iisig(3,isigF),isigA,3,isigG),Higgs,GWWH,
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,3,isigG))
                  CALL VVSXXX(Curr21(1,iflavF1,iisig(1,isigF),3),Curr43GItbv(1,iflavF3,iisig(3,isigF),isigA,4,isigG),Higgs,GWWH,
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,4,isigG))

                  AmpSumWW(iflavF1,iflavF3,isigF,isigA,isigG,2)=            !modifications with gluons on the 34-line
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,3,isigG)+
     $                 SubAmpTBVu(iflavF1,iflavF3,isigF,isigA,4,isigG)+
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,3,isigG)+
     $                 SubAmpTBVl(iflavF1,iflavF3,isigF,isigA,4,isigG)

                  do iinsA=1,4
                     AmpSumWW(iflavF1,iflavF3,isigF,isigA,isigG,2)=AmpSumWW(iflavF1,iflavF3,isigF,isigA,isigG,2)+
     $                    SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,3)+
     $                    SubAmpW(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,4)                     
                  enddo
               enddo
            enddo
         enddo
         
         do isigF=1,4
            do isigA=-1,1,2
               do isigG=-1,1,2
                  do k=1,2
                     res(5,k)=res(5,k)+dreal(AmpSumWW(4,3,isigF,isigA,isigG,k))**2+dimag(AmpSumWW(4,3,isigF,isigA,isigG,k))**2
                     res(6,k)=res(6,k)+dreal(AmpSumWW(3,4,isigF,isigA,isigG,k))**2+dimag(AmpSumWW(3,4,isigF,isigA,isigG,k))**2
                  enddo
               enddo
            enddo
         enddo


         endif !=endif gsign
      endif

*****************************************************************************
*****************************************************************************
*****************************************************************************


c     Sum NC contributions which interfere to get
c     AmpSumZZ(iflavF1,iflavF3,isigF,isigA,isigG,k)
c     where k=1/2 stands for gluon attatched to upper/lower line
      do iflavF1=3,4
         do iflavF3=3,4
            do isigF=1,4
               do isigA=-1,1,2
                  do isigG=-1,1,2
                     do iinsA=1,4
                        AmpSumZZ(iflavF1,iflavF3,isigF,isigA,isigG,1)=AmpSumZZ(iflavF1,IflavF3,isigF,isigA,isigG,1)+
     $                       SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,1)+
     $                       SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,2)
                        AmpSumZZ(iflavF1,iflavF3,isigF,isigA,isigG,2)=AmpSumZZ(iflavF1,IflavF3,isigF,isigA,isigG,2)+
     $                       SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,3)+
     $                       SubAmpZ(iflavF1,iflavF3,isigF,isigA,isigG,iinsA,4)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      
c     sum over helicities
      
      i=1
      do iflavF1=3,4
         do iflavF3=3,4
            do isigF=1,4
               do isigA=-1,1,2
                  do isigG=-1,1,2
                     do k=1,2
                        res(i,k)=res(i,k)+dreal(AmpSumZZ(iflavF1,iflavF3,isigF,isigA,isigG,k))**2+
     $                       dimag(AmpSumZZ(iflavF1,iflavF3,isigF,isigA,isigG,k))**2
                     enddo
                  enddo
               enddo
            enddo
            
            i=i+1
         enddo
      enddo
      
      do mu=1,6
         do k=1,2
            res(mu,k)=res(mu,k)*fac
         enddo
         res(mu,0)=res(mu,1)+res(mu,2)
      enddo
      do mu=1,3
         uucc(mu) = res(1,mu-1) 
         uuss(mu) = res(2,mu-1)
         ddcc(mu) = res(3,mu-1)
         ddss(mu) = res(4,mu-1)
         udsc(mu) = res(5,mu-1)
         ducs(mu) = res(6,mu-1)
      enddo

      if (ldoblha) then
        if (blha_amptype(blha_curproc).eq.0) then ! tree
          if ( (blha_idsubproc(blha_cursubproc,blha_curproc)/10 .eq. 1) .or.
     &         (blha_idsubproc(blha_cursubproc,blha_curproc)/10 .eq. 2) ) then
            blha_amp(1)=blha_amp(1)+blha_CF/(4d0/3d0)*
     &       res(mod(blha_idsubproc(blha_cursubproc,blha_curproc),10),1)
          endif
          if ( (blha_idsubproc(blha_cursubproc,blha_curproc)/10 .eq. 1) .or.
     &         (blha_idsubproc(blha_cursubproc,blha_curproc)/10 .eq. 3) ) then
            blha_amp(1)=blha_amp(1)+blha_CF/(4d0/3d0)*
     &       res(mod(blha_idsubproc(blha_cursubproc,blha_curproc),10),2)
          endif

        else 
          call BLHA_amptypeerror(blha_amptype(blha_curproc),__FILE__,__LINE__)
        endif
      endif


      return
      end
