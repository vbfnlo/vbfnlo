      subroutine HAjj_ME(pbar,pAbar,sign,nlo,uucc,uuss,ddcc,ddss,udsc,ducs)
          use globalvars, only: ldoblha

      IMPLICIT NONE

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/coupl.inc"

      integer sign(4+n_v),nlo,iisig(3,4),i,k,iflavF1,iflavF3,isigA,
     $     isigF,iinsA,iflav,isigF1,isigF3,iflavF2,iflavF4,mu,ibos
      double precision uucc,uuss,ddcc,ddss,udsc,ducs, res(6),
     $     fac,pbar(0:3,4+n_v-1),pAbar(0:3),
     $     p(0:3,4+n_v-1),pA(0:4),scale(2),pIntBosu(0:4),pIntBosl(0:4),
     $     resvirt(6),fac1,fac2,alsUpperc2o4pi,alsLowerc2o4pi,alsULcvirtc2o4pi
      double precision tree(6)
      
      double complex GBF(2,3,3:4),SpinorI(6,4,-1:1),SpinorO(6,4,-1:1),
     $     Photon(6,-1:1),Higgs(6),Curr12(6,3:4,-1:1,2:3),
     $     Curr21(6,3:4,-1:1,2:3),Curr34(6,3:4,-1:1,2:3),
     $     Curr43(6,3:4,-1:1,2:3),Curr1A2(6,3:4,-1:1,2,-1:1,2:3),
     $     Curr2A1(6,3:4,-1:1,2,-1:1,2:3),Curr3A4(6,3:4,-1:1,3:4,-1:1,2:3),
     $     Curr4A3(6,3:4,-1:1,3:4,-1:1,2:3),Curr12tbv(6,3:4,-1:1,-1:1),Curr21tbv(6,3:4,-1:1,-1:1),
     $     Curr34tbv(6,3:4,-1:1,-1:1),Curr43tbv(6,3:4,-1:1,-1:1),
     $     Curr12H(6,3:4,-1:1,2:3),Curr21H(6,3:4,-1:1,2:3),Curr34H(6,3:4,-1:1,2:3),
     $     Curr43H(6,3:4,-1:1,2:3),SubAmpZ(3:4,3:4,4,-1:1,4),SubAmpW(3:4,3:4,4,-1:1,4),
     $     AmpSumZZ(3:4,3:4,4,-1:1),AmpSumWW(3:4,3:4,4,-1:1),SubAmpTBVu(3:4,3:4,4,-1:1),
     $     SubAmpTBVl(3:4,3:4,4,-1:1),BoxL(4,3:4,3:4,-1:1,-1:1,-1:1,2:3),
     $     VirtUpper(3:4,3:4,4,-1:1,2:3),
     $     VirtLower(3:4,3:4,4,-1:1,2:3),
     $     
     $     
     $     resultgauge(2,4),resultgaugeb(2,4),
     $     BornCheck(4,3:4,3:4,-1:1,-1:1,-1:1,2:3), cte

      logical passcheck(4,3:4,3:4,-1:1,-1:1,-1:1,2:3)


      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)
      
      COMMON /HAJJCOUPL/ GBF

      COMMON /HAJJCURRENTS/ Curr12,Curr21,Curr34,Curr43,Curr1A2,Curr2A1,Curr3A4,
     $     Curr4A3,Curr12tbv,Curr21tbv,Curr34tbv,Curr43tbv,Higgs,fac
!$OMP THREADPRIVATE(/HAJJCURRENTS/)

      COMMON /HAJJCURRENTSVIRT/ SpinorI,SpinorO,Photon,Curr12H,Curr21H,Curr34H,Curr43H,
     $     pIntBosu,pIntBosl,scale,alsUpperc2o4pi,alsLowerc2o4pi,alsULcvirtc2o4pi
!$OMP THREADPRIVATE(/HAJJCURRENTSVIRT/)

      double precision c2,c2o2pi,c2o4pi
      parameter (c2=4d0/3d0, c2o2pi=c2/2d0/pi, c2o4pi=c2/4d0/pi, cte=-pi**2/3d0+7d0)
      double precision p21(0:4), p43(0:4), lnmusqp1p2, lnmusqp3p4

      double precision dotrr
      EXTERNAL dotrr
      data iisig /-1,0,-1,  -1,0,1,   1,0,-1,   1,0,1/
      parameter (fac1=1d0, fac2=1d0)
      logical doGaugecheck,resetamps   !reset amps for pspoints that fail the gauge check?
      double precision accuracy, gaugecheck
      parameter ( accuracy=1d-3, doGaugecheck=.true.,resetamps=.true..and.doGaugecheck)
      logical pspfails
      integer wccounter
      common /cwardcheck/ wccounter, pspfails

c     reset amplitudes

      do iflavF1=3,4
         do iflavF3=3,4
            do isigF=1,4
               do isigA=-1,1,2
                  AmpSumZZ(iflavF1,iflavF3,isigF,isigA)=(0d0,0d0)
                  AmpSumWW(iflavF1,iflavF3,isigF,isigA)=(0d0,0d0)
               enddo
            enddo
         enddo
      enddo
      do k=1,6
         res(k)= 0d0
         tree(k) = 0d0
      enddo

*****************************************************************************
***************************************************************************** 
*****************************************************************************
      if(sign(1).eq.1.and.
c     $     sign(2).eq.1.and.
c     $     sign(3).eq.1.and.
     $     sign(4).eq.1)then
*****************************************************************************
      do isigA=-1,1,2
c     Z contributions
         do iflavF1=3,4
            do iflavF3=3,4
               do isigF=1,4
                  CALL VVSXXX(Curr1A2(1,iflavF1,iisig(1,isigF),1,isigA,2),Curr34(1,iflavF3,iisig(3,isigF),2),Higgs,GZZH,
     $                 SubAmpZ(iflavF1,iflavF3,isigF,isigA,1)) !iinsA=1
                  CALL VVSXXX(Curr1A2(1,iflavF1,iisig(1,isigF),2,isigA,2),Curr34(1,iflavF3,iisig(3,isigF),2),Higgs,GZZH,
     $                 SubAmpZ(iflavF1,iflavF3,isigF,isigA,2)) !iinsA=2
                  CALL VVSXXX(Curr12(1,iflavF1,iisig(1,isigF),2),Curr3A4(1,iflavF3,iisig(3,isigF),3,isigA,2),Higgs,GZZH,
     $                 SubAmpZ(iflavF1,iflavF3,isigF,isigA,3)) !iinsA=3
                  CALL VVSXXX(Curr12(1,iflavF1,iisig(1,isigF),2),Curr3A4(1,iflavF3,iisig(3,isigF),4,isigA,2),Higgs,GZZH,
     $                 SubAmpZ(iflavF1,iflavF3,isigF,isigA,4)) !iinsA=4
               enddo
            enddo
c     W contributions
            iflavF3=7-iflavF1
            isigF=1             ! = (-)(-)
            CALL VVSXXX(Curr1A2(1,iflavF1,iisig(1,isigF),1,isigA,3),Curr34(1,iflavF3,iisig(3,isigF),3),Higgs,GWWH,
     $           SubAmpW(iflavF1,iflavF3,isigF,isigA,1)) !iinsA=1
            CALL VVSXXX(Curr1A2(1,iflavF1,iisig(1,isigF),2,isigA,3),Curr34(1,iflavF3,iisig(3,isigF),3),Higgs,GWWH,
     $           SubAmpW(iflavF1,iflavF3,isigF,isigA,2)) !iinsA=2
            CALL VVSXXX(Curr12(1,iflavF1,iisig(1,isigF),3),Curr3A4(1,iflavF3,iisig(3,isigF),3,isigA,3),Higgs,GWWH,
     $           SubAmpW(iflavF1,iflavF3,isigF,isigA,3)) !iinsA=3
            CALL VVSXXX(Curr12(1,iflavF1,iisig(1,isigF),3),Curr3A4(1,iflavF3,iisig(3,isigF),4,isigA,3),Higgs,GWWH,
     $           SubAmpW(iflavF1,iflavF3,isigF,isigA,4)) !iinsA=4
c     TBV contributions
            CALL VVSXXX(Curr12tbv(1,iflavF1,iisig(1,isigF),isigA),Curr34(1,iflavF3,iisig(3,isigF),3),Higgs,GWWH,
     $           SubAmpTBVu(iflavF1,iflavF3,isigF,isigA))
            CALL VVSXXX(Curr12(1,iflavF1,iisig(1,isigF),3),Curr34tbv(1,iflavF3,iisig(3,isigF),isigA),Higgs,GWWH,
     $           SubAmpTBVl(iflavF1,iflavF3,isigF,isigA))
            AmpSumWW(iflavF1,iflavF3,isigF,isigA)=SubAmpTBVu(iflavF1,iflavF3,isigF,isigA)+SubAmpTBVl(iflavF1,iflavF3,isigF,isigA)
            do iinsA=1,4
               AmpSumWW(iflavF1,iflavF3,isigF,isigA)=
     -              AmpSumWW(iflavF1,iflavF3,isigF,isigA)+
     -              SubAmpW(iflavF1,iflavF3,isigF,isigA,iinsA)
            enddo
         enddo
         res(5)=res(5)+dreal(AmpSumWW(3,4,isigF,isigA))**2+dimag(AmpSumWW(3,4,isigF,isigA))**2
         res(6)=res(6)+dreal(AmpSumWW(4,3,isigF,isigA))**2+dimag(AmpSumWW(4,3,isigF,isigA))**2
      enddo   
*****************************************************************************
      elseif(sign(1).eq.1.and.
c     $     sign(2).eq.1.and.
c     $     sign(3).eq.-1.and.
     $     sign(4).eq.-1)then
*****************************************************************************

         do isigA=-1,1,2
c     Z contributions
            do iflavF1=3,4
               do iflavF3=3,4
                  do isigF=1,4
                     CALL VVSXXX(Curr1A2(1,iflavF1,iisig(1,isigF),1,isigA,2),Curr43(1,iflavF3,iisig(3,isigF),2),Higgs,GZZH,
     $                    SubAmpZ(iflavF1,iflavF3,isigF,isigA,1)) !iinsA=1
                     CALL VVSXXX(Curr1A2(1,iflavF1,iisig(1,isigF),2,isigA,2),Curr43(1,iflavF3,iisig(3,isigF),2),Higgs,GZZH,
     $                    SubAmpZ(iflavF1,iflavF3,isigF,isigA,2)) !iinsA=2
                     CALL VVSXXX(Curr12(1,iflavF1,iisig(1,isigF),2),Curr4A3(1,iflavF3,iisig(3,isigF),3,isigA,2),Higgs,GZZH,
     $                    SubAmpZ(iflavF1,iflavF3,isigF,isigA,3)) !iinsA=3
                     CALL VVSXXX(Curr12(1,iflavF1,iisig(1,isigF),2),Curr4A3(1,iflavF3,iisig(3,isigF),4,isigA,2),Higgs,GZZH,
     $                    SubAmpZ(iflavF1,iflavF3,isigF,isigA,4)) !iinsA=4
                  enddo
               enddo
c     W contributions
               iflavF3=iflavF1
               isigF=2          ! = (-)(+)
               CALL VVSXXX(Curr1A2(1,iflavF1,iisig(1,isigF),1,isigA,3),Curr43(1,iflavF3,iisig(3,isigF),3),Higgs,GWWH,
     $              SubAmpW(iflavF1,iflavF3,isigF,isigA,1)) !iinsA=1
               CALL VVSXXX(Curr1A2(1,iflavF1,iisig(1,isigF),2,isigA,3),Curr43(1,iflavF3,iisig(3,isigF),3),Higgs,GWWH,
     $              SubAmpW(iflavF1,iflavF3,isigF,isigA,2)) !iinsA=2
               CALL VVSXXX(Curr12(1,iflavF1,iisig(1,isigF),3),Curr4A3(1,iflavF3,iisig(3,isigF),3,isigA,3),Higgs,GWWH,
     $              SubAmpW(iflavF1,iflavF3,isigF,isigA,3)) !iinsA=3
               CALL VVSXXX(Curr12(1,iflavF1,iisig(1,isigF),3),Curr4A3(1,iflavF3,iisig(3,isigF),4,isigA,3),Higgs,GWWH,
     $              SubAmpW(iflavF1,iflavF3,isigF,isigA,4)) !iinsA=4

c     TBV contributions
               CALL VVSXXX(Curr12tbv(1,iflavF1,iisig(1,isigF),isigA),Curr43(1,iflavF3,iisig(3,isigF),3),Higgs,GWWH,
     $              SubAmpTBVu(iflavF1,iflavF3,isigF,isigA))
               CALL VVSXXX(Curr12(1,iflavF1,iisig(1,isigF),3),Curr43tbv(1,iflavF3,iisig(3,isigF),isigA),Higgs,GWWH,
     $              SubAmpTBVl(iflavF1,iflavF3,isigF,isigA))
               AmpSumWW(iflavF1,iflavF3,isigF,isigA)=SubAmpTBVu(iflavF1,iflavF3,isigF,isigA)+SubAmpTBVl(iflavF1,iflavF3,isigF,isigA)
               do iinsA=1,4
                  AmpSumWW(iflavF1,iflavF3,isigF,isigA)=AmpSumWW(iflavF1,iflavF3,isigF,isigA)+
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,iinsA)
               enddo
            enddo
            res(5)=res(5)+dreal(AmpSumWW(3,3,isigF,isigA))**2+dimag(AmpSumWW(3,3,isigF,isigA))**2
            res(6)=res(6)+dreal(AmpSumWW(4,4,isigF,isigA))**2+dimag(AmpSumWW(4,4,isigF,isigA))**2
         enddo
         
*****************************************************************************
      elseif(sign(1).eq.-1.and.
c     $     sign(2).eq.-1.and.
c     $     sign(3).eq.1.and.
     $        sign(4).eq.1)then
*****************************************************************************
         do isigA=-1,1,2
c     Z contributions
            do iflavF1=3,4
               do iflavF3=3,4
                  do isigF=1,4
                     CALL VVSXXX(Curr2A1(1,iflavF1,iisig(1,isigF),1,isigA,2),Curr34(1,iflavF3,iisig(3,isigF),2),Higgs,GZZH,
     $                    SubAmpZ(iflavF1,iflavF3,isigF,isigA,1)) !iinsA=1
                     CALL VVSXXX(Curr2A1(1,iflavF1,iisig(1,isigF),2,isigA,2),Curr34(1,iflavF3,iisig(3,isigF),2),Higgs,GZZH,
     $                    SubAmpZ(iflavF1,iflavF3,isigF,isigA,2)) !iinsA=2
                     CALL VVSXXX(Curr21(1,iflavF1,iisig(1,isigF),2),Curr3A4(1,iflavF3,iisig(3,isigF),3,isigA,2),Higgs,GZZH,
     $                    SubAmpZ(iflavF1,iflavF3,isigF,isigA,3)) !iinsA=3
                     CALL VVSXXX(Curr21(1,iflavF1,iisig(1,isigF),2),Curr3A4(1,iflavF3,iisig(3,isigF),4,isigA,2),Higgs,GZZH,
     $                    SubAmpZ(iflavF1,iflavF3,isigF,isigA,4)) !iinsA=4
                  enddo
               enddo
c     W contributions
               iflavF3=iflavF1
               isigF=3          ! = (+)(-)
               CALL VVSXXX(Curr2A1(1,iflavF1,iisig(1,isigF),1,isigA,3),Curr34(1,iflavF3,iisig(3,isigF),3),Higgs,GWWH,
     $              SubAmpW(iflavF1,iflavF3,isigF,isigA,1)) !iinsA=1
               CALL VVSXXX(Curr2A1(1,iflavF1,iisig(1,isigF),2,isigA,3),Curr34(1,iflavF3,iisig(3,isigF),3),Higgs,GWWH,
     $              SubAmpW(iflavF1,iflavF3,isigF,isigA,2)) !iinsA=2
               CALL VVSXXX(Curr21(1,iflavF1,iisig(1,isigF),3),Curr3A4(1,iflavF3,iisig(3,isigF),3,isigA,3),Higgs,GWWH,
     $              SubAmpW(iflavF1,iflavF3,isigF,isigA,3)) !iinsA=3
               CALL VVSXXX(Curr21(1,iflavF1,iisig(1,isigF),3),Curr3A4(1,iflavF3,iisig(3,isigF),4,isigA,3),Higgs,GWWH,
     $              SubAmpW(iflavF1,iflavF3,isigF,isigA,4)) !iinsA=4

c     TBV contributions
               CALL VVSXXX(Curr21tbv(1,iflavF1,iisig(1,isigF),isigA),Curr34(1,iflavF3,iisig(3,isigF),3),Higgs,GWWH,
     $              SubAmpTBVu(iflavF1,iflavF3,isigF,isigA))
               CALL VVSXXX(Curr21(1,iflavF1,iisig(1,isigF),3),Curr34tbv(1,iflavF3,iisig(3,isigF),isigA),Higgs,GWWH,
     $              SubAmpTBVl(iflavF1,iflavF3,isigF,isigA))

               AmpSumWW(iflavF1,iflavF3,isigF,isigA)=SubAmpTBVu(iflavF1,iflavF3,isigF,isigA)+SubAmpTBVl(iflavF1,iflavF3,isigF,isigA)

               do iinsA=1,4
                  AmpSumWW(iflavF1,iflavF3,isigF,isigA)=AmpSumWW(iflavF1,iflavF3,isigF,isigA)+
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,iinsA)
               enddo
            enddo
            res(5)=res(5)+dreal(AmpSumWW(4,4,isigF,isigA))**2+dimag(AmpSumWW(4,4,isigF,isigA))**2
            res(6)=res(6)+dreal(AmpSumWW(3,3,isigF,isigA))**2+dimag(AmpSumWW(3,3,isigF,isigA))**2
         enddo
         
*****************************************************************************
      elseif(sign(1).eq.-1.and.
c     $     sign(2).eq.-1.and.
c     $     sign(3).eq.-1.and.
     $        sign(4).eq.-1)then
*****************************************************************************

         do isigA=-1,1,2
c     Z contributions
            do iflavF1=3,4
               do iflavF3=3,4
                  do isigF=1,4
                     CALL VVSXXX(Curr2A1(1,iflavF1,iisig(1,isigF),1,isigA,2),Curr43(1,iflavF3,iisig(3,isigF),2),Higgs,GZZH,
     $                    SubAmpZ(iflavF1,iflavF3,isigF,isigA,1)) !iinsA=1
                     CALL VVSXXX(Curr2A1(1,iflavF1,iisig(1,isigF),2,isigA,2),Curr43(1,iflavF3,iisig(3,isigF),2),Higgs,GZZH,
     $                    SubAmpZ(iflavF1,iflavF3,isigF,isigA,2)) !iinsA=2
                     CALL VVSXXX(Curr21(1,iflavF1,iisig(1,isigF),2),Curr4A3(1,iflavF3,iisig(3,isigF),3,isigA,2),Higgs,GZZH,
     $                    SubAmpZ(iflavF1,iflavF3,isigF,isigA,3)) !iinsA=3
                     CALL VVSXXX(Curr21(1,iflavF1,iisig(1,isigF),2),Curr4A3(1,iflavF3,iisig(3,isigF),4,isigA,2),Higgs,GZZH,
     $                    SubAmpZ(iflavF1,iflavF3,isigF,isigA,4)) !iinsA=4
                  enddo
               enddo
c     W contributions
               iflavF3=7-iflavF1
               isigF=4          ! = (-)(+)
               CALL VVSXXX(Curr2A1(1,iflavF1,iisig(1,isigF),1,isigA,3),Curr43(1,iflavF3,iisig(3,isigF),3),Higgs,GWWH,
     $              SubAmpW(iflavF1,iflavF3,isigF,isigA,1)) !iinsA=1
               CALL VVSXXX(Curr2A1(1,iflavF1,iisig(1,isigF),2,isigA,3),Curr43(1,iflavF3,iisig(3,isigF),3),Higgs,GWWH,
     $              SubAmpW(iflavF1,iflavF3,isigF,isigA,2)) !iinsA=2
               CALL VVSXXX(Curr21(1,iflavF1,iisig(1,isigF),3),Curr4A3(1,iflavF3,iisig(3,isigF),3,isigA,3),Higgs,GWWH,
     $              SubAmpW(iflavF1,iflavF3,isigF,isigA,3)) !iinsA=3
               CALL VVSXXX(Curr21(1,iflavF1,iisig(1,isigF),3),Curr4A3(1,iflavF3,iisig(3,isigF),4,isigA,3),Higgs,GWWH,
     $              SubAmpW(iflavF1,iflavF3,isigF,isigA,4)) !iinsA=4

c     TBV contributions
               CALL VVSXXX(Curr21tbv(1,iflavF1,iisig(1,isigF),isigA),Curr43(1,iflavF3,iisig(3,isigF),3),Higgs,GWWH,
     $              SubAmpTBVu(iflavF1,iflavF3,isigF,isigA))
               CALL VVSXXX(Curr21(1,iflavF1,iisig(1,isigF),3),Curr43tbv(1,iflavF3,iisig(3,isigF),isigA),Higgs,GWWH,
     $              SubAmpTBVl(iflavF1,iflavF3,isigF,isigA))

               AmpSumWW(iflavF1,iflavF3,isigF,isigA)=SubAmpTBVu(iflavF1,iflavF3,isigF,isigA)+SubAmpTBVl(iflavF1,iflavF3,isigF,isigA)

               do iinsA=1,4
                  AmpSumWW(iflavF1,iflavF3,isigF,isigA)=AmpSumWW(iflavF1,iflavF3,isigF,isigA)+
     $                 SubAmpW(iflavF1,iflavF3,isigF,isigA,iinsA)
               enddo
            enddo
            res(5)=res(5)+dreal(AmpSumWW(4,3,isigF,isigA))**2+dimag(AmpSumWW(4,3,isigF,isigA))**2
            res(6)=res(6)+dreal(AmpSumWW(3,4,isigF,isigA))**2+dimag(AmpSumWW(3,4,isigF,isigA))**2
         enddo
      endif
*******************************************************************************
*******************************************************************************
      
c     sum the NC contributions which interfere
      do iflavF1=3,4
         do iflavF3=3,4
            do isigF=1,4
               do isigA=-1,1,2
                  do iinsA=1,4
                     AmpSumZZ(iflavF1,iflavF3,isigF,isigA)=AmpSumZZ(iflavF1,iflavF3,isigF,isigA)+
     $                 SubAmpZ(iflavF1,iflavF3,isigF,isigA,iinsA)
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
                  res(i)=res(i)+dreal(AmpSumZZ(iflavF1,iflavF3,isigF,isigA))**2+dimag(AmpSumZZ(iflavF1,iflavF3,isigF,isigA))**2
               enddo
            enddo
            i=i+1
         enddo
      enddo

      do i=1,6
        tree(i) = res(i)*fac
      enddo

#ifdef WITH_NLO
      if(nlo.eq.1)then
         do i=1,6
            resvirt(i)=0d0
         enddo

         do iflavF1=3,4
            do iflavF3=3,4
               do isigF=1,4
                  do isigA=-1,1,2
                     do ibos=2,3
                        VirtUpper(iflavF1,iflavF3,isigF,isigA,iBos)=(0d0,0d0)
                     enddo
                  enddo
               enddo
            enddo
         enddo

         do mu = 0,3
            do i = 2,4,2
               p(mu,i) = -pbar(mu,i)
            enddo
            do i = 1,3,2
               p(mu,i) = pbar(mu,i)
            enddo
            pA(mu)=-pAbar(mu)
         enddo
         

      if(doGaugecheck)then
         do isigA=-1,1,2  
            do isigF=1,4
               do iflavF1=3,4
                  do iflavF3=3,4
                     do iinsA=1,4
                        do iBos=2,3
                           passcheck(iinsA,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,iBos)=.true.
c                              print*,sign(1),sign(2),sign(3),sign(4),iinsA,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,iBos
c                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif

*******************************************************************************
*******************************************************************************
         if(sign(1).eq.1.and.
c     $        sign(2).eq.1.and.
c     $        sign(3).eq.1.and.
     $        sign(4).eq.1)then
*******************************************************************************
c     
         i=1
         do isigA=-1,1,2 
            do isigF1=-1,1,2
               do isigF3=-1,1,2 
                  do iflavF1=3,4
                     do iflavF3=3,4
                        CALL BOXLINEMM(p(0,1),pA(0),pIntBosu(0),p(0,2),SpinorO(2-isigF1,2,isigF1),SpinorI(2+isigF1,1,isigF1),
     $                       Photon(1,isigA),Curr34H(1,iflavF3,isigF3,2),isigF1,scale(1),3,i,resultgauge(1,1),
     $                       BoxL(1,iflavF1,iflavF3,isigF1,isigF3,isigA,2),resultgaugeb(1,1),
     $                       BornCheck(1,iflavF1,iflavF3,isigF1,isigF3,isigA,2)) !iinsA=1
                        BoxL(1,iflavF1,iflavF3,isigF1,isigF3,isigA,2)=
     $                       BoxL(1,iflavF1,iflavF3,isigF1,isigF3,isigA,2)*GBF((isigF1+3)/2,1,iflavF1)*GBF((isigF1+3)/2,2,iflavF1)
                        if(doGaugecheck)then
                           call box_gaugeMM(accuracy,p(0,1),pA(0),pIntBosu(0),p(0,2),scale(1),resultgauge(1,1),resultgaugeb(1,1),
     $                          gaugecheck,passcheck(1,iflavF1,iflavF3,isigF1,isigF3,isigA,2))
                        endif
c$$$  BornCheck(1,iflavF1,iflavF3,isigF1,isigF3,isigA,2)=
c$$$  $                          BornCheck(1,iflavF1,iflavF3,isigF1,isigF3,isigA,2)*GBF((isigF1+3)/2,1,iflavF1)*GBF((isigF1+3)/2,2,iflavF1)
                        i=-1
                     enddo
                  enddo
               enddo
            enddo
            do iflav=3,4
               iflavF1=iflav
               iflavF2=7-iflav
               iflavF3=7-iflav
               iflavF4=iflav
               isigF1=-1
               isigF3=-1
               CALL BOXLINEMM(p(0,1),pA(0),pIntBosu(0),p(0,2),
     -              SpinorO(2-isigF1,2,isigF1),SpinorI(2+isigF1,1,
     -              isigF1),Photon(1,isigA),
     $              Curr34H(1,iflavF3,isigF3,3),isigF1,scale(1),3,i,
     -              resultgauge(1,1),BoxL(1,iflavF1,iflavF3,isigF1,
     -              isigF3,isigA,3),resultgaugeb(1,1),
     $              BornCheck(1,iflavF1,iflavF3,isigF1,isigF3,isigA,3)) !iinsA=1               
               BoxL(1,iflavF1,iflavF3,isigF1,isigF3,isigA,3)=
     $              BoxL(1,iflavF1,iflavF3,isigF1,isigF3,isigA,3)*
     -              GBF((isigF1+3)/2,1,iflavF1)*GBF((isigF1+3)/2,3,iflavF2)
               if(doGaugecheck)then
                  call box_gaugeMM(accuracy,p(0,1),pA(0),pIntBosu(0),
     -              p(0,2),scale(1),resultgauge(1,1),resultgaugeb(1,1),
     -              gaugecheck,
     $                 passcheck(1,iflavF1,iflavF3,isigF1,isigF3,isigA,3))  
c                  if(.not.passcheck(1,iflavF1,iflavF3,isigF1,isigF3,isigA,3))                  print*,"gaugecheck=",gaugecheck
               endif
c$$$  BornCheck(1,iflavF1,iflavF3,isigF1,isigF3,isigA,3)=
c$$$  $                 BornCheck(1,iflavF1,iflavF3,isigF1,isigF3,isigA,3)*GBF((isigF1+3)/2,1,iflavF1)*GBF((isigF1+3)/2,3,iflavF2) 
            enddo   
         enddo

         i=1
         do isigA=-1,1,2 
            do isigF1=-1,1,2
               do isigF3=-1,1,2 
                  do iflavF1=3,4
                     do iflavF3=3,4
                        CALL BOXLINEMM(p(0,1),pIntBosu(0),pA(0),p(0,2),
     -                       SpinorO(2-isigF1,2,isigF1),SpinorI(2+isigF1,
     -                       1,isigF1),
     $                       Curr34H(1,iflavF3,isigF3,2),Photon(1,isigA),
     -                       isigF1,scale(1),3,i,resultgauge(1,2),
     $                       BoxL(2,iflavF1,iflavF3,isigF1,isigF3,isigA,
     -                       2),resultgaugeb(1,2),
     $                       BornCheck(2,iflavF1,iflavF3,isigF1,isigF3,
     -                       isigA,2)) !iinsA=2
                        BoxL(2,iflavF1,iflavF3,isigF1,isigF3,isigA,2)=
     $                       BoxL(2,iflavF1,iflavF3,isigF1,isigF3,isigA,
     -                       2)*GBF((isigF1+3)/2,1,iflavF1)*
     -                       GBF((isigF1+3)/2,2,iflavF1)
                        if(doGaugecheck)
     $                       call box_gaugeMM(accuracy,p(0,1),pIntBosu(0),
     -                       pA(0),p(0,2),scale(1),resultgauge(1,2),
     -                       resultgaugeb(1,2),
     $                       gaugecheck,passcheck(2,iflavF1,iflavF3,
     -                       isigF1,isigF3,isigA,2))      
c$$$  BornCheck(2,iflavF1,iflavF3,isigF1,isigF3,isigA,2)=
c$$$  $                          BornCheck(2,iflavF1,iflavF3,isigF1,isigF3,isigA,2)*GBF((isigF1+3)/2,1,iflavF1)*GBF((isigF1+3)/2,2,iflavF1)
                        i=-1
                     enddo
                  enddo
               enddo
            enddo
            do iflav=3,4
               iflavF1=iflav
               iflavF2=7-iflav
               iflavF3=7-iflav
               iflavF4=iflav
               isigF1=-1
               isigF3=-1
               CALL BOXLINEMM(p(0,1),pIntBosu(0),pA(0),p(0,2),
     -              SpinorO(2-isigF1,2,isigF1),SpinorI(2+isigF1,1,isigF1),
     -              Curr34H(1,iflavF3,isigF3,3),
     $              Photon(1,isigA),isigF1,scale(1),3,i,resultgauge(1,2),
     -              BoxL(2,iflavF1,iflavF3,isigF1,isigF3,isigA,3),
     -              resultgaugeb(1,2),
     $              BornCheck(2,iflavF1,iflavF3,isigF1,isigF3,isigA,3)) !iinsA=2 
               BoxL(2,iflavF1,iflavF3,isigF1,isigF3,isigA,3)=
     $              BoxL(2,iflavF1,iflavF3,isigF1,isigF3,isigA,3)*
     -              GBF((isigF1+3)/2,1,iflavF2)*GBF((isigF1+3)/2,3,iflavF1)
                        if(doGaugecheck)
     $                       call box_gaugeMM(accuracy,p(0,1),
     -              pIntBosu(0),pA(0),p(0,2),scale(1),resultgauge(1,2),
     -              resultgaugeb(1,2),gaugecheck,
     $                       passcheck(2,iflavF1,iflavF3,isigF1,isigF3,
     -              isigA,3))  
c$$$  BornCheck(2,iflavF1,iflavF3,isigF1,isigF3,isigA,3)=
c$$$  $                 BornCheck(2,iflavF1,iflavF3,isigF1,isigF3,isigA,3)*GBF((isigF1+3)/2,1,iflavF2)*GBF((isigF1+3)/2,3,iflavF1)                
            enddo
         enddo 

         i=1
         do isigA=-1,1,2
            do isigF1=-1,1,2
               do isigF3=-1,1,2  
                  do iflavF1=3,4
                     do iflavF3=3,4
                        CALL BOXLINEMM(p(0,3),pA(0),pIntBosl(0),p(0,4),
     -                       SpinorO(2-isigF3,4,isigF3),
     -                       SpinorI(2+isigF3,3,isigF3),
     $                       Photon(1,isigA),Curr12H(1,iflavF1,isigF1,2),
     -                       isigF3,scale(2),3,i,resultgauge(1,3),
     $                       BoxL(3,iflavF1,iflavF3,isigF1,isigF3,
     -                       isigA,2),resultgaugeb(1,3),
     $                       BornCheck(3,iflavF1,iflavF3,isigF1,isigF3,
     -                       isigA,2)) !iinsA=3
                        BoxL(3,iflavF1,iflavF3,isigF1,isigF3,isigA,2)=
     $                       BoxL(3,iflavF1,iflavF3,isigF1,isigF3,isigA,
     -                       2)*GBF((isigF3+3)/2,1,iflavF3)*GBF((isigF3+
     -                       3)/2,2,iflavF3)
                        if(doGaugecheck)
     $                       call box_gaugeMM(accuracy,p(0,3),pA(0),
     -                       pIntBosl(0),p(0,4),scale(2),resultgauge(1,3),
     -                       resultgaugeb(1,3),
     $                       gaugecheck,passcheck(3,iflavF1,iflavF3,
     -                       isigF1,isigF3,isigA,2))  
c$$$  BornCheck(3,iflavF1,iflavF3,isigF1,isigF3,isigA,2)=
c$$$  $                          BornCheck(3,iflavF1,iflavF3,isigF1,isigF3,isigA,2)*GBF((isigF3+3)/2,1,iflavF3)*GBF((isigF3+3)/2,2,iflavF3)
                        i=-1
                     enddo
                  enddo
               enddo
            enddo
            do iflav=3,4
               iflavF1=iflav
               iflavF2=7-iflav
               iflavF3=7-iflav
               iflavF4=iflav
               isigF1=-1
               isigF3=-1
               CALL BOXLINEMM(p(0,3),pA(0),pIntBosl(0),p(0,4),
     -                       SpinorO(2-isigF3,4,isigF3),
     -                       SpinorI(2+isigF3,3,isigF3),Photon(1,isigA),
     $              Curr12H(1,iflavF1,isigF1,3),isigF3,scale(2),3,i,
     -                       resultgauge(1,3),BoxL(3,iflavF1,iflavF3,
     -                       isigF1,isigF3,isigA,3),
     $              resultgaugeb(1,3),BornCheck(3,iflavF1,iflavF3,
     -                       isigF1,isigF3,isigA,3)) !iinsA=3
               BoxL(3,iflavF1,iflavF3,isigF1,isigF3,isigA,3)=
     $              BoxL(3,iflavF1,iflavF3,isigF1,isigF3,isigA,3)*
     -                       GBF((isigF3+3)/2,1,iflavF3)*GBF((isigF3+3)/
     -                       2,3,iflavF4)
               if(doGaugecheck)
     $              call box_gaugeMM(accuracy,p(0,3),pA(0),pIntBosl(0),
     -                       p(0,4),scale(2),resultgauge(1,3),
     -                       resultgaugeb(1,3),
     $              gaugecheck,passcheck(3,iflavF1,iflavF3,isigF1,isigF3,
     -                       isigA,3))  
c$$$  BornCheck(3,iflavF1,iflavF3,isigF1,isigF3,isigA,3)=
c$$$  $                 BornCheck(3,iflavF1,iflavF3,isigF1,isigF3,isigA,3)*GBF((isigF3+3)/2,1,iflavF3)*GBF((isigF3+3)/2,3,iflavF4)
            enddo
         enddo

         i=1
         do isigA=-1,1,2
            do isigF1=-1,1,2
               do isigF3=-1,1,2
                  do iflavF1=3,4
                     do iflavF3=3,4  
                        CALL BOXLINEMM(p(0,3),pIntBosl(0),pA(0),p(0,4),SpinorO(2-isigF3,4,isigF3),SpinorI(2+isigF3,3,isigF3),
     $                       Curr12H(1,iflavF1,isigF1,2),Photon(1,isigA),isigF3,scale(2),3,i,resultgauge(1,4),
     $                       BoxL(4,iflavF1,iflavF3,isigF1,isigF3,isigA,2),resultgaugeb(1,4),
     $                       BornCheck(4,iflavF1,iflavF3,isigF1,isigF3,isigA,2)) !iinsA=4
                        BoxL(4,iflavF1,iflavF3,isigF1,isigF3,isigA,2)=
     $                       BoxL(4,iflavF1,iflavF3,isigF1,isigF3,isigA,2)*GBF((isigF3+3)/2,1,iflavF3)*GBF((isigF3+3)/2,2,iflavF3)
                        if(doGaugecheck)
     $                       call box_gaugeMM(accuracy,p(0,3),pIntBosl(0),pA(0),p(0,4),scale(2),resultgauge(1,4),resultgaugeb(1,4),
     $                       gaugecheck,passcheck(4,iflavF1,iflavF3,isigF1,isigF3,isigA,2))  
c$$$  BornCheck(4,iflavF1,iflavF3,isigF1,isigF3,isigA,2)=
c$$$  $                          BornCheck(4,iflavF1,iflavF3,isigF1,isigF3,isigA,2)*GBF((isigF3+3)/2,1,iflavF3)*GBF((isigF3+3)/2,2,iflavF3)
                        i=-1
                     enddo
                  enddo
               enddo
            enddo
            do iflav=3,4
               iflavF1=iflav
               iflavF2=7-iflav
               iflavF3=7-iflav
               iflavF4=iflav
               isigF1=-1
               isigF3=-1
               CALL BOXLINEMM(p(0,3),pIntBosl(0),pA(0),p(0,4),
     $                       SpinorO(2-isigF3,4,isigF3),
     $                       SpinorI(2+isigF3,3,isigF3),
     $              Curr12H(1,iflavF1,isigF1,3),Photon(1,isigA),
     $                       isigF3,scale(2),3,i,resultgauge(1,4),
     $              BoxL(4,iflavF1,iflavF3,isigF1,isigF3,isigA,3),
     $                       resultgaugeb(1,4),
     $              BornCheck(4,iflavF1,iflavF3,isigF1,isigF3,
     $                       isigA,3)) !iinsA=4
               BoxL(4,iflavF1,iflavF3,isigF1,isigF3,isigA,3)=
     $              BoxL(4,iflavF1,iflavF3,isigF1,isigF3,isigA,3)*
     $                       GBF((isigF3+3)/2,1,iflavF4)*GBF((isigF3+3)/
     $                       2,3,iflavF3)
               if(doGaugecheck)
     $              call box_gaugeMM(accuracy,p(0,3),pIntBosl(0),pA(0),
     $                       p(0,4),scale(2),resultgauge(1,4),
     $                       resultgaugeb(1,4),gaugecheck,
     $              passcheck(4,iflavF1,iflavF3,isigF1,isigF3,isigA,3))  
c$$$  BornCheck(4,iflavF1,iflavF3,isigF1,isigF3,isigA,3)=
c$$$  $                 BornCheck(4,iflavF1,iflavF3,isigF1,isigF3,isigA,3)*GBF((isigF3+3)/2,1,iflavF4)*GBF((isigF3+3)/2,3,iflavF3)
            enddo
         enddo


         i=1
         do isigA=-1,1,2  
            do isigF=1,4
               do iflavF1=3,4
                  do iflavF3=3,4
!Z contributions

                     if(resetamps)then
                        if(passcheck(1,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,2).and.
     $                       passcheck(2,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,2))
     $                       VirtUpper(iflavF1,iflavF3,isigF,isigA,2)=
     $                       (BoxL(1,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,2)+BoxL(2,iflavF1,iflavF3,iisig(1,isigF),
     $                       iisig(3,isigF),isigA,2))-(SubAmpZ(iflavF1,iflavF3,isigF,isigA,1)+
     $                       SubAmpZ(iflavF1,iflavF3,isigF,isigA,2))*cte
                        if(passcheck(3,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,2).and.
     $                       passcheck(4,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,2))
     $                       VirtLower(iflavF1,iflavF3,isigF,isigA,2)=
     $                       (BoxL(3,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,2)+BoxL(4,iflavF1,iflavF3,iisig(1,isigF),
     $                       iisig(3,isigF),isigA,2))-(SubAmpZ(iflavF1,iflavF3,isigF,isigA,3)+
     $                       SubAmpZ(iflavF1,iflavF3,isigF,isigA,4))*cte
                     else
                        VirtUpper(iflavF1,iflavF3,isigF,isigA,2)=
     $                       (BoxL(1,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,2)+BoxL(2,iflavF1,iflavF3,iisig(1,isigF),
     $                       iisig(3,isigF),isigA,2))-(SubAmpZ(iflavF1,iflavF3,isigF,isigA,1)+
     $                       SubAmpZ(iflavF1,iflavF3,isigF,isigA,2))*cte
                        VirtLower(iflavF1,iflavF3,isigF,isigA,2)=
     $                       (BoxL(3,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,2)+BoxL(4,iflavF1,iflavF3,iisig(1,isigF),
     $                       iisig(3,isigF),isigA,2))-(SubAmpZ(iflavF1,iflavF3,isigF,isigA,3)+
     $                       SubAmpZ(iflavF1,iflavF3,isigF,isigA,4))*cte
                     endif
                  enddo
               enddo
            enddo
            
            do iflav=3,4
               iflavF1=iflav
               iflavF2=7-iflav
               iflavF3=7-iflav
               iflavF4=iflav
               isigF=1
!W contributions
               if(doGaugecheck)then
                  if(passcheck(1,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3).and.
     $                 passcheck(2,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3))
     $                 VirtUpper(iflavF1,iflavF3,isigF,isigA,3)=
     $                 (BoxL(1,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3)+BoxL(2,iflavF1,iflavF3,iisig(1,isigF),
     $                 iisig(3,isigF),isigA,3))-(SubAmpW(iflavF1,iflavF3,isigF,isigA,1)+SubAmpW(iflavF1,iflavF3,isigF,isigA,2))*cte
                  
                  if(passcheck(3,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3).and.
     $                 passcheck(4,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3))
     $                 VirtLower(iflavF1,iflavF3,isigF,isigA,3)=
     $                 (BoxL(3,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3)+BoxL(4,iflavF1,iflavF3,iisig(1,isigF),
     $                 iisig(3,isigF),isigA,3))-(SubAmpW(iflavF1,iflavF3,isigF,isigA,3)+SubAmpW(iflavF1,iflavF3,isigF,isigA,4))*cte
               else
                  VirtUpper(iflavF1,iflavF3,isigF,isigA,3)=
     $                 (BoxL(1,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3)+BoxL(2,iflavF1,iflavF3,iisig(1,isigF),
     $                 iisig(3,isigF),isigA,3))-(SubAmpW(iflavF1,iflavF3,isigF,isigA,1)+SubAmpW(iflavF1,iflavF3,isigF,isigA,2))*cte
                  
                  VirtLower(iflavF1,iflavF3,isigF,isigA,3)=
     $                 (BoxL(3,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3)+BoxL(4,iflavF1,iflavF3,iisig(1,isigF),
     $                 iisig(3,isigF),isigA,3))-(SubAmpW(iflavF1,iflavF3,isigF,isigA,3)+SubAmpW(iflavF1,iflavF3,isigF,isigA,4))*cte
               endif
            enddo
         enddo

c     isigF=1
         do isigA=-1,1,2
            resvirt(5)=resvirt(5)
     $           +2*dreal(AmpSumWW(3,4,isigF,isigA)*conjg(-alsUpperc2o4pi*VirtUpper(3,4,isigF,isigA,3)
     $           -alsLowerc2o4pi*VirtLower(3,4,isigF,isigA,3)+alsULcvirtc2o4pi*AmpSumWW(3,4,isigF,isigA)))
            resvirt(6)=resvirt(6)
     $           +2*dreal(AmpSumWW(4,3,isigF,isigA)*conjg(-alsUpperc2o4pi*VirtUpper(4,3,isigF,isigA,3)
     $           -alsLowerc2o4pi*VirtLower(4,3,isigF,isigA,3)+alsULcvirtc2o4pi*AmpSumWW(4,3,isigF,isigA)))
         enddo
         

*******************************************************************************
      elseif(sign(1).eq.1.and.
c     $        sign(2).eq.1.and.
c     $        sign(3).eq.-1.and.
     $        sign(4).eq.-1)then
*******************************************************************************
c     
         i=1
         do isigA=-1,1,2 
            do isigF1=-1,1,2
               do isigF3=-1,1,2 
                  do iflavF1=3,4
                     do iflavF3=3,4
                        CALL BOXLINEMM(p(0,1),pA(0),pIntBosu(0),p(0,2),SpinorO(2-isigF1,2,isigF1),SpinorI(2+isigF1,1,isigF1),
     $                       Photon(1,isigA),Curr43H(1,iflavF3,isigF3,2),isigF1,scale(1),3,i,resultgauge(1,1),
     $                       BoxL(1,iflavF1,iflavF3,isigF1,isigF3,isigA,2),resultgaugeb(1,1),
     $                       BornCheck(1,iflavF1,iflavF3,isigF1,isigF3,isigA,2)) !iinsA=1
                        BoxL(1,iflavF1,iflavF3,isigF1,isigF3,isigA,2)=
     $                       BoxL(1,iflavF1,iflavF3,isigF1,isigF3,isigA,2)*GBF((isigF1+3)/2,1,iflavF1)*GBF((isigF1+3)/2,2,iflavF1)
c$$$                        print*,"resultgauge(1,1)",resultgauge(1,1)
c$$$                        passcheck(1,iflavF1,iflavF3,isigF1,isigF3,isigA,2)=.true.
                        if(doGaugecheck)
     $                       call box_gaugeMM(accuracy,p(0,1),pA(0),pIntBosu(0),p(0,2),scale(1),resultgauge(1,1),resultgaugeb(1,1),
     $                       gaugecheck,passcheck(1,iflavF1,iflavF3,isigF1,isigF3,isigA,2))  
c$$$  BornCheck(1,iflavF1,iflavF3,isigF1,isigF3,isigA,2)=
c$$$  $                          BornCheck(1,iflavF1,iflavF3,isigF1,isigF3,isigA,2)*GBF((isigF1+3)/2,1,iflavF1)*GBF((isigF1+3)/2,2,iflavF1)
                        i=-1
                     enddo
                  enddo
               enddo
            enddo
            do iflav=3,4
               iflavF1=iflav
               iflavF2=7-iflav
               iflavF3=iflav
               iflavF4=7-iflav
               isigF1=-1
               isigF3=1
               CALL BOXLINEMM(p(0,1),pA(0),pIntBosu(0),p(0,2),
     -              SpinorO(2-isigF1,2,isigF1),SpinorI(2+isigF1,1,isigF1),
     -              Photon(1,isigA),
     $              Curr43H(1,iflavF3,isigF3,3),isigF1,scale(1),3,i,
     -              resultgauge(1,1),BoxL(1,iflavF1,iflavF3,isigF1,
     -              isigF3,isigA,3),
     $              resultgaugeb(1,1),BornCheck(1,iflavF1,iflavF3,
     -              isigF1,isigF3,isigA,3)) !iinsA=1
               
               BoxL(1,iflavF1,iflavF3,isigF1,isigF3,isigA,3)=
     $              BoxL(1,iflavF1,iflavF3,isigF1,isigF3,isigA,3)*
     -              GBF((isigF1+3)/2,1,iflavF1)*GBF((isigF1+3)/2,3,iflavF2)
               if(doGaugecheck)
     $              call box_gaugeMM(accuracy,p(0,1),pA(0),pIntBosu(0),
     -              p(0,2),scale(1),resultgauge(1,1),resultgaugeb(1,1),
     -              gaugecheck,
     $              passcheck(1,iflavF1,iflavF3,isigF1,isigF3,isigA,3))  
c$$$  BornCheck(1,iflavF1,iflavF3,isigF1,isigF3,isigA,3)=
c$$$  $                 BornCheck(1,iflavF1,iflavF3,isigF1,isigF3,isigA,3)*GBF((isigF1+3)/2,1,iflavF1)*GBF((isigF1+3)/2,3,iflavF2) 
            enddo   
         enddo

         i=1
         do isigA=-1,1,2 
            do isigF1=-1,1,2
               do isigF3=-1,1,2 
                  do iflavF1=3,4
                     do iflavF3=3,4
                        CALL BOXLINEMM(p(0,1),pIntBosu(0),pA(0),p(0,2),
     -                       SpinorO(2-isigF1,2,isigF1),
     -                       SpinorI(2+isigF1,1,isigF1),
     $                       Curr43H(1,iflavF3,isigF3,2),
     -                       Photon(1,isigA),i
     -                       sigF1,scale(1),3,i,resultgauge(1,2),
     $                       BoxL(2,iflavF1,iflavF3,isigF1,isigF3,
     -                       isigA,2),resultgaugeb(1,2),
     $                       BornCheck(2,iflavF1,iflavF3,isigF1,isigF3,
     -                       isigA,2)) !iinsA=2
                        BoxL(2,iflavF1,iflavF3,isigF1,isigF3,isigA,2)=
     $                       BoxL(2,iflavF1,iflavF3,isigF1,isigF3,
     -                       isigA,2)*GBF((isigF1+3)/2,1,iflavF1)*
     -                       GBF((isigF1+3)/2,2,iflavF1)
                        if(doGaugecheck)
     $                       call box_gaugeMM(accuracy,p(0,1),
     -                       pIntBosu(0),pA(0),p(0,2),scale(1),
     -                       resultgauge(1,2),resultgaugeb(1,2),
     $                       gaugecheck,passcheck(2,iflavF1,iflavF3,
     -                       isigF1,isigF3,isigA,2))  
c$$$  BornCheck(2,iflavF1,iflavF3,isigF1,isigF3,isigA,2)=
c$$$  $                          BornCheck(2,iflavF1,iflavF3,isigF1,isigF3,isigA,2)*GBF((isigF1+3)/2,1,iflavF1)*GBF((isigF1+3)/2,2,iflavF1)
                        i=-1
                     enddo
                  enddo
               enddo
            enddo
            do iflav=3,4
               iflavF1=iflav
               iflavF2=7-iflav
               iflavF3=iflav
               iflavF4=7-iflav
               isigF1=-1
               isigF3=1
               CALL BOXLINEMM(p(0,1),pIntBosu(0),pA(0),p(0,2),
     -              SpinorO(2-isigF1,2,isigF1),SpinorI(2+isigF1,1,isigF1),
     $              Curr43H(1,iflavF3,isigF3,3),Photon(1,isigA),isigF1,
     -              scale(1),3,i,resultgauge(1,2),
     $              BoxL(2,iflavF1,iflavF3,isigF1,isigF3,isigA,3),
     -              resultgaugeb(1,2),
     $              BornCheck(2,iflavF1,iflavF3,isigF1,isigF3,isigA,3)) !iinsA=2 
               BoxL(2,iflavF1,iflavF3,isigF1,isigF3,isigA,3)=
     $              BoxL(2,iflavF1,iflavF3,isigF1,isigF3,isigA,3)*
     -              GBF((isigF1+3)/2,1,iflavF2)*GBF((isigF1+3)/2,3,iflavF1)
               if(doGaugecheck)
     $              call box_gaugeMM(accuracy,p(0,1),pIntBosu(0),pA(0),
     -              p(0,2),scale(1),resultgauge(1,2),resultgaugeb(1,2),
     -              gaugecheck,
     $              passcheck(2,iflavF1,iflavF3,isigF1,isigF3,isigA,3))  
c$$$  BornCheck(2,iflavF1,iflavF3,isigF1,isigF3,isigA,3)=
c$$$  $                 BornCheck(2,iflavF1,iflavF3,isigF1,isigF3,isigA,3)*GBF((isigF1+3)/2,1,iflavF2)*GBF((isigF1+3)/2,3,iflavF1)                
            enddo
         enddo 

         i=1
         do isigA=-1,1,2
            do isigF1=-1,1,2
               do isigF3=-1,1,2  
                  do iflavF1=3,4
                     do iflavF3=3,4
                        CALL BOXLINEMM(p(0,4),pIntBosl(0),pA(0),p(0,3),
     -              SpinorO(2+isigF3,3,isigF3),SpinorI(2-isigF3,4,isigF3),
     $                       Curr12H(1,iflavF1,isigF1,2),
     -              Photon(1,isigA),-isigF3,scale(2),3,i,resultgauge(1,3),
     $                       BoxL(3,iflavF1,iflavF3,isigF1,isigF3,
     -              isigA,2),resultgaugeb(1,3),
     $                       BornCheck(3,iflavF1,iflavF3,isigF1,isigF3,
     -              isigA,2)) !iinsA=3
                        BoxL(3,iflavF1,iflavF3,isigF1,isigF3,isigA,2)=
     $                       BoxL(3,iflavF1,iflavF3,isigF1,isigF3,
     -              isigA,2)*GBF((-isigF3+3)/2,1,iflavF3)*
     -              GBF((-isigF3+3)/2,2,iflavF3)
                        if(doGaugecheck)
     $                       call box_gaugeMM(accuracy,p(0,4),
     -                       pIntBosl(0),pA(0),p(0,3),scale(2),
     -                       resultgauge(1,3),resultgaugeb(1,3),
     $                       gaugecheck,passcheck(3,iflavF1,iflavF3,
     -                       isigF1,isigF3,isigA,2))
c$$$  BornCheck(3,iflavF1,iflavF3,isigF1,isigF3,isigA,2)=
c$$$  $                          BornCheck(3,iflavF1,iflavF3,isigF1,isigF3,isigA,2)*GBF((-isigF3+3)/2,1,iflavF3)*GBF((-isigF3+3)/2,2,iflavF3)
                        i=-1
                     enddo
                  enddo
               enddo
            enddo
            do iflav=3,4
               iflavF1=iflav
               iflavF2=7-iflav
               iflavF3=iflav
               iflavF4=7-iflav
               isigF1=-1
               isigF3=1
               CALL BOXLINEMM(p(0,4),pIntBosl(0),pA(0),p(0,3),
     -              SpinorO(2+isigF3,3,isigF3),SpinorI(2-isigF3,4,isigF3),
     $              Curr12H(1,iflavF1,isigF1,3),Photon(1,isigA),-isigF3,
     -              scale(2),3,i,resultgauge(1,3),
     $              BoxL(3,iflavF1,iflavF3,isigF1,isigF3,isigA,3),
     -              resultgaugeb(1,3),
     $              BornCheck(3,iflavF1,iflavF3,isigF1,isigF3,isigA,3)) !iinsA=3
               BoxL(3,iflavF1,iflavF3,isigF1,isigF3,isigA,3)=
     $              BoxL(3,iflavF1,iflavF3,isigF1,isigF3,isigA,3)*
     -              GBF((-isigF3+3)/2,1,iflavF3)*GBF((-isigF3+3)/2,3,iflavF4)
               if(doGaugecheck)
     $              call box_gaugeMM(accuracy,p(0,4),pIntBosl(0),pA(0),
     -              p(0,3),scale(2),resultgauge(1,3),resultgaugeb(1,3),
     -              gaugecheck,
     $              passcheck(3,iflavF1,iflavF3,isigF1,isigF3,isigA,3))
c$$$  BornCheck(3,iflavF1,iflavF3,isigF1,isigF3,isigA,3)=
c$$$  $                        BornCheck(3,iflavF1,iflavF3,isigF1,isigF3,isigA,3)*GBF((-isigF3+3)/2,1,iflavF3)*GBF((-isigF3+3)/2,3,iflavF4)
            enddo
         enddo

         i=1
         do isigA=-1,1,2
            do isigF1=-1,1,2
               do isigF3=-1,1,2
                  do iflavF1=3,4
                     do iflavF3=3,4                            
                        CALL BOXLINEMM(p(0,4),pA(0),pIntBosl(0),p(0,3),
     -              SpinorO(2+isigF3,3,isigF3),SpinorI(2-isigF3,4,isigF3),
     $                       Photon(1,isigA),
     -                       Curr12H(1,iflavF1,isigF1,2),-isigF3,
     -                       scale(2),3,i,resultgauge(1,4),
     $                       BoxL(4,iflavF1,iflavF3,isigF1,isigF3,
     -                       isigA,2),resultgaugeb(1,4),
     $                       BornCheck(4,iflavF1,iflavF3,isigF1,isigF3,
     -                       isigA,2)) !iinsA=4
                        BoxL(4,iflavF1,iflavF3,isigF1,isigF3,isigA,2)=
     $                       BoxL(4,iflavF1,iflavF3,isigF1,isigF3,
     -                       isigA,2)*GBF((-isigF3+3)/2,1,iflavF3)*
     -                       GBF((-isigF3+3)/2,2,iflavF3)
                        if(doGaugecheck)
     $                       call box_gaugeMM(accuracy,p(0,4),pA(0),
     -                       pIntBosl(0),p(0,3),scale(2),
     -                       resultgauge(1,4),resultgaugeb(1,4),
     $                       gaugecheck,passcheck(4,iflavF1,
     -                       iflavF3,isigF1,isigF3,isigA,2))
c$$$  BornCheck(4,iflavF1,iflavF3,isigF1,isigF3,isigA,2)=
c$$$  $                          BornCheck(4,iflavF1,iflavF3,isigF1,isigF3,isigA,2)*GBF((-isigF3+3)/2,1,iflavF3)*GBF((-isigF3+3)/2,2,iflavF3)
                        i=-1
                     enddo
                  enddo
               enddo
            enddo
            do iflav=3,4
               iflavF1=iflav
               iflavF2=7-iflav
               iflavF3=iflav
               iflavF4=7-iflav
               isigF1=-1
               isigF3=1
               CALL BOXLINEMM(p(0,4),pA(0),pIntBosl(0),p(0,3),
     -              SpinorO(2+isigF3,3,isigF3),SpinorI(2-
     -              isigF3,4,isigF3),Photon(1,isigA),
     $              Curr12H(1,iflavF1,isigF1,3),sign(3)*isigF3,scale(2),
     -              3,i,resultgauge(1,4),
     $              BoxL(4,iflavF1,iflavF3,isigF1,isigF3,isigA,3),
     -              resultgaugeb(1,4),
     $              BornCheck(4,iflavF1,iflavF3,isigF1,isigF3,isigA,3)) !iinsA=4
               BoxL(4,iflavF1,iflavF3,isigF1,isigF3,isigA,3)=
     $              BoxL(4,iflavF1,iflavF3,isigF1,isigF3,isigA,3)*
     -              GBF((-isigF3+3)/2,1,iflavF4)*GBF((-isigF3+3)/2,3,iflavF3)
               if(doGaugecheck)
     $              call box_gaugeMM(accuracy,p(0,4),pA(0),
     -              pIntBosl(0),p(0,3),scale(2),resultgauge(1,4),
     -              resultgaugeb(1,4),gaugecheck,
     $              passcheck(4,iflavF1,iflavF3,isigF1,isigF3,isigA,3))
c$$$  BornCheck(4,iflavF1,iflavF3,isigF1,isigF3,isigA,3)=
c$$$  $                 BornCheck(4,iflavF1,iflavF3,isigF1,isigF3,isigA,3)*GBF((-isigF3+3)/2,1,iflavF4)*GBF((-isigF3+3)/2,3,iflavF3)
            enddo
         enddo 


         i=1
         do isigA=-1,1,2  
            do isigF=1,4
               do iflavF1=3,4
                  do iflavF3=3,4
!Z contributions
                     if(resetamps)then
                        if(passcheck(1,iflavF1,iflavF3,iisig(1,isigF),
     -                       iisig(3,isigF),isigA,2).and.
     $                       passcheck(2,iflavF1,iflavF3,iisig(1,isigF),
     -                       iisig(3,isigF),isigA,2))
     $                       VirtUpper(iflavF1,iflavF3,isigF,isigA,2)=
     $                       (BoxL(1,iflavF1,iflavF3,iisig(1,isigF),
     -                       iisig(3,isigF),isigA,2)
     $                       +BoxL(2,iflavF1,iflavF3,iisig(1,isigF),
     -                       iisig(3,isigF),isigA,2))
     $                       -(SubAmpZ(iflavF1,iflavF3,isigF,isigA,1)+
     -                       SubAmpZ(iflavF1,iflavF3,isigF,isigA,2))*cte
                        if(passcheck(3,iflavF1,iflavF3,iisig(1,isigF),
     -                       iisig(3,isigF),isigA,2).and.
     $                       passcheck(4,iflavF1,iflavF3,iisig(1,isigF),
     -                       iisig(3,isigF),isigA,2))
     $                       VirtLower(iflavF1,iflavF3,isigF,isigA,2)=
     $                       (BoxL(3,iflavF1,iflavF3,iisig(1,isigF),
     -                       iisig(3,isigF),isigA,2)
     $                       +BoxL(4,iflavF1,iflavF3,iisig(1,isigF),
     -                       iisig(3,isigF),isigA,2))
     $                       -(SubAmpZ(iflavF1,iflavF3,isigF,isigA,3)+
     -                       SubAmpZ(iflavF1,iflavF3,isigF,isigA,4))*cte
                     else
                        VirtUpper(iflavF1,iflavF3,isigF,isigA,2)=
     $                       (BoxL(1,iflavF1,iflavF3,iisig(1,isigF),
     -                       iisig(3,isigF),isigA,2)
     $                       +BoxL(2,iflavF1,iflavF3,iisig(1,isigF),
     -                       iisig(3,isigF),isigA,2))
     $                       -(SubAmpZ(iflavF1,iflavF3,isigF,isigA,1)+
     -                       SubAmpZ(iflavF1,iflavF3,isigF,isigA,2))*cte
                        VirtLower(iflavF1,iflavF3,isigF,isigA,2)=
     $                       (BoxL(3,iflavF1,iflavF3,iisig(1,isigF),
     -                       iisig(3,isigF),isigA,2)
     $                       +BoxL(4,iflavF1,iflavF3,iisig(1,isigF),
     -                       iisig(3,isigF),isigA,2))
     $                       -(SubAmpZ(iflavF1,iflavF3,isigF,isigA,3)+
     -                       SubAmpZ(iflavF1,iflavF3,isigF,isigA,4))*cte
                     endif
                  enddo
               enddo
            enddo
            
            do iflav=3,4
               iflavF1=iflav
               iflavF2=7-iflav
               iflavF3=iflav
               iflavF4=7-iflav
               isigF=2
!W contributions
               if(resetamps)then
                  if(passcheck(1,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3).and.
     $                 passcheck(2,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3))
     $                 VirtUpper(iflavF1,iflavF3,isigF,isigA,3)=
     $                 (BoxL(1,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3)
     $                 +BoxL(2,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3))
     $                 -(SubAmpW(iflavF1,iflavF3,isigF,isigA,1)+SubAmpW(iflavF1,iflavF3,isigF,isigA,2))*cte
                  if(passcheck(3,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3).and.
     $                 passcheck(4,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3))
     $                 VirtLower(iflavF1,iflavF3,isigF,isigA,3)=
     $                 (BoxL(3,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3)
     $                 +BoxL(4,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3))
     $                 -(SubAmpW(iflavF1,iflavF3,isigF,isigA,3)+SubAmpW(iflavF1,iflavF3,isigF,isigA,4))*cte
               else
                  VirtUpper(iflavF1,iflavF3,isigF,isigA,3)=
     $                 (BoxL(1,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3)
     $                 +BoxL(2,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3))
     $                 -(SubAmpW(iflavF1,iflavF3,isigF,isigA,1)+SubAmpW(iflavF1,iflavF3,isigF,isigA,2))*cte
                  VirtLower(iflavF1,iflavF3,isigF,isigA,3)=
     $                 (BoxL(3,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3)
     $                 +BoxL(4,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3))
     $                 -(SubAmpW(iflavF1,iflavF3,isigF,isigA,3)+SubAmpW(iflavF1,iflavF3,isigF,isigA,4))*cte
               endif
            enddo
         enddo

c     isigF=2
         do isigA=-1,1,2
            resvirt(5)=resvirt(5)
     $           +2*dreal(AmpSumWW(3,3,isigF,isigA)*conjg(-alsUpperc2o4pi*VirtUpper(3,3,isigF,isigA,3)
     $           -alsLowerc2o4pi*VirtLower(3,3,isigF,isigA,3)+alsULcvirtc2o4pi*AmpSumWW(3,3,isigF,isigA)))
            resvirt(6)=resvirt(6)
     $           +2*dreal(AmpSumWW(4,4,isigF,isigA)*conjg(-alsUpperc2o4pi*VirtUpper(4,4,isigF,isigA,3)
     $           -alsLowerc2o4pi*VirtLower(4,4,isigF,isigA,3)+alsULcvirtc2o4pi*AmpSumWW(4,4,isigF,isigA)))
         enddo
         

*******************************************************************************
      elseif(sign(1).eq.-1.and.
c     $        sign(2).eq.-1.and.
c     $        sign(3).eq.1.and.
     $        sign(4).eq.1)then
*******************************************************************************
c     
         i=1
         do isigA=-1,1,2 
            do isigF1=-1,1,2
               do isigF3=-1,1,2 
                  do iflavF1=3,4
                     do iflavF3=3,4                           
                        CALL BOXLINEMM(p(0,2),pIntBosu(0),pA(0),p(0,1),SpinorO(2+isigF1,1,isigF1),SpinorI(2-isigF1,2,isigF1),
     $                       Curr34H(1,iflavF3,isigF3,2),Photon(1,isigA),-isigF1,scale(1),3,i,resultgauge(1,1),
     $                       BoxL(1,iflavF1,iflavF3,isigF1,isigF3,isigA,2),resultgaugeb(1,1),
     $                       BornCheck(1,iflavF1,iflavF3,isigF1,isigF3,isigA,2)) !iinsA=1
                        BoxL(1,iflavF1,iflavF3,isigF1,isigF3,isigA,2)=
     $                       BoxL(1,iflavF1,iflavF3,isigF1,isigF3,isigA,2)*GBF((-isigF1+3)/2,1,iflavF1)*GBF((-isigF1+3)/2,2,iflavF1)
                        if(doGaugecheck)
     $                       call box_gaugeMM(accuracy,p(0,2),pIntBosu(0),pA(0),p(0,1),scale(1),resultgauge(1,1),resultgaugeb(1,1),
     $                       gaugecheck,passcheck(1,iflavF1,iflavF3,isigF1,isigF3,isigA,2))

c$$$  BornCheck(1,iflavF1,iflavF3,isigF1,isigF3,isigA,2)=
c$$$  $                          BornCheck(1,iflavF1,iflavF3,isigF1,isigF3,isigA,2)*GBF((-isigF1+3)/2,1,iflavF1)*GBF((-isigF1+3)/2,2,iflavF1)
                        i=-1
                     enddo
                  enddo
               enddo
            enddo
            do iflav=3,4
               iflavF1=iflav
               iflavF2=7-iflav
               iflavF3=iflav
               iflavF4=7-iflav
               isigF1=1
               isigF3=-1
               CALL BOXLINEMM(p(0,2),pIntBosu(0),pA(0),p(0,1),
     -              SpinorO(2+isigF1,1,isigF1),SpinorI(2-isigF1,2,isigF1),
     $              Curr34H(1,iflavF3,isigF3,3),Photon(1,isigA),-isigF1,
     -              scale(1),3,i,resultgauge(1,1),
     $              BoxL(1,iflavF1,iflavF3,isigF1,isigF3,isigA,3),
     -              resultgaugeb(1,1),
     $              BornCheck(1,iflavF1,iflavF3,isigF1,isigF3,isigA,3)) !iinsA=1
               BoxL(1,iflavF1,iflavF3,isigF1,isigF3,isigA,3)=
     $              BoxL(1,iflavF1,iflavF3,isigF1,isigF3,isigA,3)*
     -              GBF((-isigF1+3)/2,1,iflavF1)*GBF((-isigF1+3)/2,3,iflavF2)
               if(doGaugecheck)
     $              call box_gaugeMM(accuracy,p(0,2),pIntBosu(0),pA(0),
     -              p(0,1),scale(1),resultgauge(1,1),resultgaugeb(1,1),
     -              gaugecheck,
     $              passcheck(1,iflavF1,iflavF3,isigF1,isigF3,isigA,3))
c$$$  BornCheck(1,iflavF1,iflavF3,isigF1,isigF3,isigA,3)=
c$$$  $                        BornCheck(1,iflavF1,iflavF3,isigF1,isigF3,isigA,3)*GBF((-isigF1+3)/2,1,iflavF1)*GBF((-isigF1+3)/2,3,iflavF2) 
            enddo   
         enddo

         i=1
         do isigA=-1,1,2 
            do isigF1=-1,1,2
               do isigF3=-1,1,2 
                  do iflavF1=3,4
                     do iflavF3=3,4
                        CALL BOXLINEMM(p(0,2),pA(0),pIntBosu(0),p(0,1),
     -                       SpinorO(2+isigF1,1,isigF1),SpinorI(2-isigF1,
     -                       2,isigF1),
     $                       Photon(1,isigA),Curr34H(1,iflavF3,isigF3,2),
     -                       -isigF1,scale(1),3,i,resultgauge(1,2),
     $                       BoxL(2,iflavF1,iflavF3,isigF1,isigF3,isigA,2),
     -                       resultgaugeb(1,2),
     $                       BornCheck(2,iflavF1,iflavF3,isigF1,isigF3,
     -                       isigA,2)) !iinsA=2
                        BoxL(2,iflavF1,iflavF3,isigF1,isigF3,isigA,2)=
     $                       BoxL(2,iflavF1,iflavF3,isigF1,isigF3,
     -                       isigA,2)*GBF((-isigF1+3)/2,1,iflavF1)*
     -                       GBF((-isigF1+3)/2,2,iflavF1)
                        if(doGaugecheck)
     $                       call box_gaugeMM(accuracy,p(0,2),pA(0),
     -                       pIntBosu(0),p(0,1),scale(1),resultgauge(1,2),
     -                       resultgaugeb(1,2),
     $                       gaugecheck,passcheck(2,iflavF1,iflavF3,
     -                       isigF1,isigF3,isigA,2))
c$$$  BornCheck(2,iflavF1,iflavF3,isigF1,isigF3,isigA,2)=
c$$$  $                          BornCheck(2,iflavF1,iflavF3,isigF1,isigF3,isigA,2)*GBF((-isigF1+3)/2,1,iflavF1)*GBF((-isigF1+3)/2,2,iflavF1)
                        i=-1
                     enddo
                  enddo
               enddo
            enddo
            do iflav=3,4
               iflavF1=iflav
               iflavF2=7-iflav
               iflavF3=iflav
               iflavF4=7-iflav
               isigF1=1
               isigF3=-1
               CALL BOXLINEMM(p(0,2),pA(0),pIntBosu(0),p(0,1),
     -              SpinorO(2+isigF1,1,isigF1),SpinorI(2-isigF1,2,isigF1),
     -              Photon(1,isigA),
     $              Curr34H(1,iflavF3,isigF3,3),-isigF1,scale(1),3,i,
     -              resultgauge(1,2),
     $              BoxL(2,iflavF1,iflavF3,isigF1,isigF3,isigA,3),
     -              resultgaugeb(1,2),
     $              BornCheck(2,iflavF1,iflavF3,isigF1,isigF3,isigA,3)) !iinsA=2
               BoxL(2,iflavF1,iflavF3,isigF1,isigF3,isigA,3)=
     $              BoxL(2,iflavF1,iflavF3,isigF1,isigF3,isigA,3)*
     -              GBF((-isigF1+3)/2,1,iflavF2)*GBF((-isigF1+3)/2,3,
     -              iflavF1)
               if(doGaugecheck)
     $              call box_gaugeMM(accuracy,p(0,2),pA(0),pIntBosu(0),
     -              p(0,1),scale(1),resultgauge(1,2),resultgaugeb(1,2),
     -              gaugecheck,
     $              passcheck(2,iflavF1,iflavF3,isigF1,isigF3,isigA,3))
c$$$  BornCheck(2,iflavF1,iflavF3,isigF1,isigF3,isigA,3)=
c$$$  $                 BornCheck(2,iflavF1,iflavF3,isigF1,isigF3,isigA,3)*GBF((-isigF1+3)/2,1,iflavF2)*GBF((-isigF1+3)/2,3,iflavF1)              
            enddo
         enddo 

         i=1
         do isigA=-1,1,2
            do isigF1=-1,1,2
               do isigF3=-1,1,2  
                  do iflavF1=3,4
                     do iflavF3=3,4
                        CALL BOXLINEMM(p(0,3),pA(0),pIntBosl(0),p(0,4),
     -                       SpinorO(2-isigF3,4,isigF3),SpinorI(2+
     -                       isigF3,3,isigF3),
     $                       Photon(1,isigA),Curr21H(1,iflavF1,isigF1,2),
     -                       isigF3,scale(2),3,i,resultgauge(1,3),
     $                       BoxL(3,iflavF1,iflavF3,isigF1,isigF3,isigA,
     -                       2),resultgaugeb(1,3),
     $                       BornCheck(3,iflavF1,iflavF3,isigF1,isigF3,
     -                       isigA,2)) !iinsA=3
                        BoxL(3,iflavF1,iflavF3,isigF1,isigF3,isigA,2)=
     $                       BoxL(3,iflavF1,iflavF3,isigF1,isigF3,isigA,
     -                       2)*GBF((isigF3+3)/2,1,iflavF3)*GBF((isigF3+
     -                       3)/2,2,iflavF3)
                        if(doGaugecheck)
     $                       call box_gaugeMM(accuracy,p(0,3),pA(0),
     -                       pIntBosl(0),p(0,4),scale(2),resultgauge(1,3),
     -                       resultgaugeb(1,3),
     $                       gaugecheck,passcheck(3,iflavF1,iflavF3,
     -                       isigF1,isigF3,isigA,2))
c$$$  BornCheck(3,iflavF1,iflavF3,isigF1,isigF3,isigA,2)=
c$$$  $                          BornCheck(3,iflavF1,iflavF3,isigF1,isigF3,isigA,2)*GBF((isigF3+3)/2,1,iflavF3)*GBF((isigF3+3)/2,2,iflavF3)
                        i=-1
                     enddo
                  enddo
               enddo
            enddo
            do iflav=3,4
               iflavF1=iflav
               iflavF2=7-iflav
               iflavF3=iflav
               iflavF4=7-iflav
               isigF1=1
               isigF3=-1
               CALL BOXLINEMM(p(0,3),pA(0),pIntBosl(0),p(0,4),
     -              SpinorO(2-isigF3,4,isigF3),SpinorI(2+isigF3,3,isigF3),
     -              Photon(1,isigA),
     $              Curr21H(1,iflavF1,isigF1,3),isigF3,scale(2),3,i,
     -              resultgauge(1,3),BoxL(3,iflavF1,iflavF3,isigF1,
     -              isigF3,isigA,3),
     $              resultgaugeb(1,3),BornCheck(3,iflavF1,iflavF3,isigF1,
     -              isigF3,isigA,3)) !iinsA=3
               BoxL(3,iflavF1,iflavF3,isigF1,isigF3,isigA,3)=
     $              BoxL(3,iflavF1,iflavF3,isigF1,isigF3,isigA,3)*
     -              GBF((isigF3+3)/2,1,iflavF3)*GBF((isigF3+3)/2,3,iflavF4)
               if(doGaugecheck)
     $              call box_gaugeMM(accuracy,p(0,3),pA(0),pIntBosl(0),
     -              p(0,4),scale(2),resultgauge(1,3),resultgaugeb(1,3),
     -              gaugecheck,
     $              passcheck(3,iflavF1,iflavF3,isigF1,isigF3,isigA,3))
c$$$  BornCheck(3,iflavF1,iflavF3,isigF1,isigF3,isigA,3)=
c$$$  $                 BornCheck(3,iflavF1,iflavF3,isigF1,isigF3,isigA,3)*GBF((isigF3+3)/2,1,iflavF3)*GBF((isigF3+3)/2,3,iflavF4)
            enddo
         enddo

         i=1
         do isigA=-1,1,2
            do isigF1=-1,1,2
               do isigF3=-1,1,2
                  do iflavF1=3,4
                     do iflavF3=3,4  
                        CALL BOXLINEMM(p(0,3),pIntBosl(0),pA(0),p(0,4),
     -                       SpinorO(2-isigF3,4,isigF3),SpinorI(2+isigF3,
     -                       3,isigF3),
     $                       Curr21H(1,iflavF1,isigF1,2),Photon(1,isigA),
     -                       isigF3,scale(2),3,i,resultgauge(1,4),
     $                       BoxL(4,iflavF1,iflavF3,isigF1,isigF3,isigA,
     -                       2),resultgaugeb(1,4),
     $                       BornCheck(4,iflavF1,iflavF3,isigF1,isigF3,
     -                       isigA,2)) !iinsA=4
                        BoxL(4,iflavF1,iflavF3,isigF1,isigF3,isigA,2)=
     $                       BoxL(4,iflavF1,iflavF3,isigF1,isigF3,isigA,
     -                       2)*GBF((isigF3+3)/2,1,iflavF3)*GBF((isigF3+
     -                       3)/2,2,iflavF3)
                        if(doGaugecheck)
     $                       call box_gaugeMM(accuracy,p(0,3),pIntBosl(0),
     -                       pA(0),p(0,4),scale(2),resultgauge(1,4),
     -                       resultgaugeb(1,4),
     $                       gaugecheck,passcheck(4,iflavF1,iflavF3,
     -                       isigF1,isigF3,isigA,2))
c$$$  BornCheck(4,iflavF1,iflavF3,isigF1,isigF3,isigA,2)=
c$$$  $                          BornCheck(4,iflavF1,iflavF3,isigF1,isigF3,isigA,2)*GBF((isigF3+3)/2,1,iflavF3)*GBF((isigF3+3)/2,2,iflavF3)
                        i=-1
                     enddo
                  enddo
               enddo
            enddo
            do iflav=3,4
               iflavF1=iflav
               iflavF2=7-iflav
               iflavF3=iflav
               iflavF4=7-iflav
               isigF1=1
               isigF3=-1
               CALL BOXLINEMM(p(0,3),pIntBosl(0),pA(0),p(0,4),
     -              SpinorO(2-isigF3,4,isigF3),SpinorI(2+isigF3,3,isigF3),
     $              Curr21H(1,iflavF1,isigF1,3),Photon(1,isigA),isigF3,
     -              scale(2),3,i,resultgauge(1,4),
     $              BoxL(4,iflavF1,iflavF3,isigF1,isigF3,isigA,3),
     -              resultgaugeb(1,4),
     $              BornCheck(4,iflavF1,iflavF3,isigF1,isigF3,isigA,3)) !iinsA=4
               BoxL(4,iflavF1,iflavF3,isigF1,isigF3,isigA,3)=
     $              BoxL(4,iflavF1,iflavF3,isigF1,isigF3,isigA,3)*
     -              GBF((isigF3+3)/2,1,iflavF4)*GBF((isigF3+3)/2,3,iflavF3)
               if(doGaugecheck)
     $              call box_gaugeMM(accuracy,p(0,3),pIntBosl(0),pA(0),
     -              p(0,4),scale(2),resultgauge(1,4),resultgaugeb(1,4),
     -              gaugecheck,
     $              passcheck(4,iflavF1,iflavF3,isigF1,isigF3,isigA,3))
c$$$  BornCheck(4,iflavF1,iflavF3,isigF1,isigF3,isigA,3)=
c$$$  $                 BornCheck(4,iflavF1,iflavF3,isigF1,isigF3,isigA,3)*GBF((isigF3+3)/2,1,iflavF4)*GBF((isigF3+3)/2,3,iflavF3)
            enddo
         enddo 


         i=1
         do isigA=-1,1,2  
            do isigF=1,4
               do iflavF1=3,4
                  do iflavF3=3,4
!Z contributions
                     if(resetamps)then
                        if(passcheck(1,iflavF1,iflavF3,iisig(1,isigF),
     -                       iisig(3,isigF),isigA,2).and.
     $                       passcheck(2,iflavF1,iflavF3,iisig(1,isigF),
     -                       iisig(3,isigF),isigA,2))
     $                       VirtUpper(iflavF1,iflavF3,isigF,isigA,2)=
     $                       (BoxL(1,iflavF1,iflavF3,iisig(1,isigF),
     -                       iisig(3,isigF),isigA,2)
     $                       +BoxL(2,iflavF1,iflavF3,iisig(1,isigF),
     -                       iisig(3,isigF),isigA,2))
     $                       -(SubAmpZ(iflavF1,iflavF3,isigF,isigA,1)+
     -                       SubAmpZ(iflavF1,iflavF3,isigF,isigA,2))*cte
                        if(passcheck(3,iflavF1,iflavF3,iisig(1,isigF),
     -                       iisig(3,isigF),isigA,2).and.
     $                       passcheck(4,iflavF1,iflavF3,iisig(1,isigF),
     -                       iisig(3,isigF),isigA,2))
     $                       VirtLower(iflavF1,iflavF3,isigF,isigA,2)=
     $                       (BoxL(3,iflavF1,iflavF3,iisig(1,isigF),
     -                       iisig(3,isigF),isigA,2)
     $                       +BoxL(4,iflavF1,iflavF3,iisig(1,isigF),
     -                       iisig(3,isigF),isigA,2))
     $                       -(SubAmpZ(iflavF1,iflavF3,isigF,isigA,3)+
     -                       SubAmpZ(iflavF1,iflavF3,isigF,isigA,4))*cte
                     else
                        VirtUpper(iflavF1,iflavF3,isigF,isigA,2)=
     $                       (BoxL(1,iflavF1,iflavF3,iisig(1,isigF),
     -                       iisig(3,isigF),isigA,2)
     $                       +BoxL(2,iflavF1,iflavF3,iisig(1,isigF),
     -                       iisig(3,isigF),isigA,2))
     $                       -(SubAmpZ(iflavF1,iflavF3,isigF,isigA,1)+
     -                       SubAmpZ(iflavF1,iflavF3,isigF,isigA,2))*cte
                        VirtLower(iflavF1,iflavF3,isigF,isigA,2)=
     $                       (BoxL(3,iflavF1,iflavF3,iisig(1,isigF),
     -                       iisig(3,isigF),isigA,2)
     $                       +BoxL(4,iflavF1,iflavF3,iisig(1,isigF),
     -                       iisig(3,isigF),isigA,2))
     $                       -(SubAmpZ(iflavF1,iflavF3,isigF,isigA,3)+
     -                       SubAmpZ(iflavF1,iflavF3,isigF,isigA,4))*cte
                     endif
                  enddo
               enddo
            enddo
            
            do iflav=3,4
               iflavF1=iflav
               iflavF2=7-iflav
               iflavF3=iflav
               iflavF4=7-iflav
               isigF=3
!W contributions
               if(doGaugecheck)then
                  if(passcheck(1,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3).and.
     $                 passcheck(2,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3))
     $                 VirtUpper(iflavF1,iflavF3,isigF,isigA,3)=
     $                 (BoxL(1,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3)
     $                 +BoxL(2,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3))
     $                 -(SubAmpW(iflavF1,iflavF3,isigF,isigA,1)+SubAmpW(iflavF1,iflavF3,isigF,isigA,2))*cte
                  if(passcheck(3,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3).and.
     $                 passcheck(4,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3))
     $                 VirtLower(iflavF1,iflavF3,isigF,isigA,3)=
     $                 (BoxL(3,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3)
     $                 +BoxL(4,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3))
     $                 -(SubAmpW(iflavF1,iflavF3,isigF,isigA,3)+SubAmpW(iflavF1,iflavF3,isigF,isigA,4))*cte
               else
                  VirtUpper(iflavF1,iflavF3,isigF,isigA,3)=
     $                 (BoxL(1,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3)
     $                 +BoxL(2,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3))
     $                 -(SubAmpW(iflavF1,iflavF3,isigF,isigA,1)+SubAmpW(iflavF1,iflavF3,isigF,isigA,2))*cte
                  
                  VirtLower(iflavF1,iflavF3,isigF,isigA,3)=
     $                 (BoxL(3,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3)
     $                 +BoxL(4,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3))
     $                 -(SubAmpW(iflavF1,iflavF3,isigF,isigA,3)+SubAmpW(iflavF1,iflavF3,isigF,isigA,4))*cte
               endif
            enddo
         enddo

c     isigF=3
         do isigA=-1,1,2
            resvirt(5)=resvirt(5)
     $           +2*dreal(AmpSumWW(4,4,isigF,isigA)*conjg(-alsUpperc2o4pi*VirtUpper(4,4,isigF,isigA,3)
     $           -alsLowerc2o4pi*VirtLower(4,4,isigF,isigA,3)+alsULcvirtc2o4pi*AmpSumWW(4,4,isigF,isigA)))
            resvirt(6)=resvirt(6)
     $           +2*dreal(AmpSumWW(3,3,isigF,isigA)*conjg(-alsUpperc2o4pi*VirtUpper(3,3,isigF,isigA,3)
     $           -alsLowerc2o4pi*VirtLower(3,3,isigF,isigA,3)+alsULcvirtc2o4pi*AmpSumWW(3,3,isigF,isigA)))
         enddo
*******************************************************************************
      elseif(sign(1).eq.-1.and.
c     $        sign(2).eq.-1.and.
c     $        sign(3).eq.-1.and.
     $        sign(4).eq.-1)then
*******************************************************************************
c     
         i=1
         do isigA=-1,1,2 
            do isigF1=-1,1,2
               do isigF3=-1,1,2 
                  do iflavF1=3,4
                     do iflavF3=3,4                           
                        CALL BOXLINEMM(p(0,2),pIntBosu(0),pA(0),p(0,1),SpinorO(2+isigF1,1,isigF1),SpinorI(2-isigF1,2,isigF1),
     $                       Curr43H(1,iflavF3,isigF3,2),Photon(1,isigA),-isigF1,scale(1),3,i,resultgauge(1,1),
     $                       BoxL(1,iflavF1,iflavF3,isigF1,isigF3,isigA,2),resultgaugeb(1,1),
     $                       BornCheck(1,iflavF1,iflavF3,isigF1,isigF3,isigA,2)) !iinsA=1

                        BoxL(1,iflavF1,iflavF3,isigF1,isigF3,isigA,2)=
     $                       BoxL(1,iflavF1,iflavF3,isigF1,isigF3,isigA,2)*GBF((-isigF1+3)/2,1,iflavF1)*GBF((-isigF1+3)/2,2,iflavF1)
                        if(doGaugecheck)
     $                       call box_gaugeMM(accuracy,p(0,2),pIntBosu(0),pA(0),p(0,1),scale(1),resultgauge(1,1),resultgaugeb(1,1),
     $                       gaugecheck,passcheck(1,iflavF1,iflavF3,isigF1,isigF3,isigA,2))
c$$$  BornCheck(1,iflavF1,iflavF3,isigF1,isigF3,isigA,2)=
c$$$  $                          BornCheck(1,iflavF1,iflavF3,isigF1,isigF3,isigA,2)*GBF((-isigF1+3)/2,1,iflavF1)*GBF((-isigF1+3)/2,2,iflavF1)
                        i=-1
                     enddo
                  enddo
               enddo
            enddo
            do iflav=3,4
               iflavF1=iflav
               iflavF2=7-iflav
               iflavF3=7-iflav
               iflavF4=iflav
               isigF1=1
               isigF3=1
               CALL BOXLINEMM(p(0,2),pIntBosu(0),pA(0),p(0,1),
     -              SpinorO(2+isigF1,1,isigF1),SpinorI(2-isigF1,2,isigF1),
     $              Curr43H(1,iflavF3,isigF3,3),Photon(1,isigA),-isigF1,
     -              scale(1),3,i,resultgauge(1,1),
     $              BoxL(1,iflavF1,iflavF3,isigF1,isigF3,isigA,3),
     -              resultgaugeb(1,1),
     $              BornCheck(1,iflavF1,iflavF3,isigF1,isigF3,isigA,3)) !iinsA=1
               BoxL(1,iflavF1,iflavF3,isigF1,isigF3,isigA,3)=
     $              BoxL(1,iflavF1,iflavF3,isigF1,isigF3,isigA,3)*
     -              GBF((-isigF1+3)/2,1,iflavF1)*GBF((-isigF1+3)/2,3,iflavF2)
               if(doGaugecheck)
     $              call box_gaugeMM(accuracy,p(0,2),pIntBosu(0),pA(0),
     -              p(0,1),scale(1),resultgauge(1,1),resultgaugeb(1,1),
     -              gaugecheck,
     $              passcheck(1,iflavF1,iflavF3,isigF1,isigF3,isigA,3))
c$$$  BornCheck(1,iflavF1,iflavF3,isigF1,isigF3,isigA,3)=
c$$$  $                        BornCheck(1,iflavF1,iflavF3,isigF1,isigF3,isigA,3)*GBF((-isigF1+3)/2,1,iflavF1)*GBF((-isigF1+3)/2,3,iflavF2) 
            enddo   
         enddo

         i=1
         do isigA=-1,1,2 
            do isigF1=-1,1,2
               do isigF3=-1,1,2 
                  do iflavF1=3,4
                     do iflavF3=3,4
                        CALL BOXLINEMM(p(0,2),pA(0),pIntBosu(0),p(0,1),
     -                       SpinorO(2+isigF1,1,isigF1),SpinorI(2-
     -                       isigF1,2,isigF1),
     $                       Photon(1,isigA),Curr43H(1,iflavF3,isigF3,2),
     -                       -isigF1,scale(1),3,i,resultgauge(1,2),
     $                       BoxL(2,iflavF1,iflavF3,isigF1,isigF3,
     -                       isigA,2),resultgaugeb(1,2),
     $                       BornCheck(2,iflavF1,iflavF3,isigF1,isigF3,
     -                       isigA,2)) !iinsA=2
                        BoxL(2,iflavF1,iflavF3,isigF1,isigF3,isigA,2)=
     $                       BoxL(2,iflavF1,iflavF3,isigF1,isigF3,isigA,
     -                       2)*GBF((-isigF1+3)/2,1,iflavF1)*
     -                       GBF((-isigF1+3)/2,2,iflavF1)
                        if(doGaugecheck)
     $                       call box_gaugeMM(accuracy,p(0,2),pA(0),
     -                       pIntBosu(0),p(0,1),scale(1),resultgauge(1,2),
     -                       resultgaugeb(1,2),
     $                       gaugecheck,passcheck(2,iflavF1,iflavF3,
     -                       isigF1,isigF3,isigA,2))
c$$$  BornCheck(2,iflavF1,iflavF3,isigF1,isigF3,isigA,2)=
c$$$  $                          BornCheck(2,iflavF1,iflavF3,isigF1,isigF3,isigA,2)*GBF((-isigF1+3)/2,1,iflavF1)*GBF((-isigF1+3)/2,2,iflavF1)
                        i=-1
                     enddo
                  enddo
               enddo
            enddo
            do iflav=3,4
               iflavF1=iflav
               iflavF2=7-iflav
               iflavF3=7-iflav
               iflavF4=iflav
               isigF1=1
               isigF3=1
               CALL BOXLINEMM(p(0,2),pA(0),pIntBosu(0),p(0,1),
     -              SpinorO(2+isigF1,1,isigF1),SpinorI(2-isigF1,2,isigF1),
     -              Photon(1,isigA),
     $              Curr43H(1,iflavF3,isigF3,3),-isigF1,scale(1),3,i,
     -              resultgauge(1,2),BoxL(2,iflavF1,iflavF3,isigF1,
     -              isigF3,isigA,3),
     $              resultgaugeb(1,2),BornCheck(2,iflavF1,iflavF3,isigF1,
     -              isigF3,isigA,3)) !iinsA=2
               BoxL(2,iflavF1,iflavF3,isigF1,isigF3,isigA,3)=
     $              BoxL(2,iflavF1,iflavF3,isigF1,isigF3,isigA,3)*
     -              GBF((-isigF1+3)/2,1,iflavF2)*GBF((-isigF1+3)/2,3,
     -              iflavF1)
               if(doGaugecheck)
     $              call box_gaugeMM(accuracy,p(0,2),pA(0),pIntBosu(0),
     -              p(0,1),scale(1),resultgauge(1,2),resultgaugeb(1,2),
     -              gaugecheck,
     $              passcheck(2,iflavF1,iflavF3,isigF1,isigF3,isigA,3))
c$$$  BornCheck(2,iflavF1,iflavF3,isigF1,isigF3,isigA,3)=
c$$$  $              BornCheck(2,iflavF1,iflavF3,isigF1,isigF3,isigA,3)*GBF((-isigF1+3)/2,1,iflavF2)*GBF((-isigF1+3)/2,3,iflavF1)              
            enddo
         enddo 

         i=1
         do isigA=-1,1,2
            do isigF1=-1,1,2
               do isigF3=-1,1,2  
                  do iflavF1=3,4
                     do iflavF3=3,4
                        CALL BOXLINEMM(p(0,4),pIntBosl(0),pA(0),p(0,3),
     -                       SpinorO(2+isigF3,3,isigF3),SpinorI(2-
     -                       isigF3,4,isigF3),
     $                       Curr21H(1,iflavF1,isigF1,2),Photon(1,isigA),
     -                       -isigF3,scale(2),3,i,resultgauge(1,3),
     $                       BoxL(3,iflavF1,iflavF3,isigF1,isigF3,isigA,2),
     -                       resultgaugeb(1,3),
     $                       BornCheck(3,iflavF1,iflavF3,isigF1,isigF3,
     -                       isigA,2)) !iinsA=3
                        BoxL(3,iflavF1,iflavF3,isigF1,isigF3,isigA,2)=
     $                       BoxL(3,iflavF1,iflavF3,isigF1,isigF3,
     -                       isigA,2)*GBF((-isigF3+3)/2,1,iflavF3)*
     -                       GBF((-isigF3+3)/2,2,iflavF3)
                        if(doGaugecheck)
     $                       call box_gaugeMM(accuracy,p(0,4),pIntBosl(0),
     -                       pA(0),p(0,3),scale(2),resultgauge(1,3),
     -                       resultgaugeb(1,3),
     $                       gaugecheck,passcheck(3,iflavF1,iflavF3,
     -                       isigF1,isigF3,isigA,2))
c$$$  BornCheck(3,iflavF1,iflavF3,isigF1,isigF3,isigA,2)=
c$$$  $                          BornCheck(3,iflavF1,iflavF3,isigF1,isigF3,isigA,2)*GBF((-isigF3+3)/2,1,iflavF3)*GBF((-isigF3+3)/2,2,iflavF3)
                        i=-1
                     enddo
                  enddo
               enddo
            enddo
            do iflav=3,4
               iflavF1=iflav
               iflavF2=7-iflav
               iflavF3=7-iflav
               iflavF4=iflav
               isigF1=1
               isigF3=1
               CALL BOXLINEMM(p(0,4),pIntBosl(0),pA(0),p(0,3),
     -              SpinorO(2+isigF3,3,isigF3),SpinorI(2-isigF3,4,isigF3),
     $              Curr21H(1,iflavF1,isigF1,3),Photon(1,isigA),
     -              -isigF3,scale(2),3,i,resultgauge(1,3),
     $              BoxL(3,iflavF1,iflavF3,isigF1,isigF3,isigA,3),
     -              resultgaugeb(1,3),
     $              BornCheck(3,iflavF1,iflavF3,isigF1,isigF3,isigA,3)) !iinsA=3
               BoxL(3,iflavF1,iflavF3,isigF1,isigF3,isigA,3)=
     $              BoxL(3,iflavF1,iflavF3,isigF1,isigF3,isigA,3)*
     -              GBF((-isigF3+3)/2,1,iflavF3)*GBF((-isigF3+3)/2,3,
     -              iflavF4)
               if(doGaugecheck)
     $              call box_gaugeMM(accuracy,p(0,4),pIntBosl(0),pA(0),
     -              p(0,3),scale(2),resultgauge(1,3),resultgaugeb(1,3),
     -              gaugecheck,
     $              passcheck(3,iflavF1,iflavF3,isigF1,isigF3,isigA,3))
c$$$  BornCheck(3,iflavF1,iflavF3,isigF1,isigF3,isigA,3)=
c$$$  $                        BornCheck(3,iflavF1,iflavF3,isigF1,isigF3,isigA,3)*GBF((-isigF3+3)/2,1,iflavF3)*GBF((-isigF3+3)/2,3,iflavF4)
            enddo
         enddo

         i=1
         do isigA=-1,1,2
            do isigF1=-1,1,2
               do isigF3=-1,1,2
                  do iflavF1=3,4
                     do iflavF3=3,4                            
                        CALL BOXLINEMM(p(0,4),pA(0),pIntBosl(0),p(0,3),
     -                       SpinorO(2+isigF3,3,isigF3),SpinorI(2-
     -                       isigF3,4,isigF3),
     $                       Photon(1,isigA),Curr21H(1,iflavF1,isigF1,
     -                       2),-isigF3,scale(2),3,i,resultgauge(1,4),
     $                       BoxL(4,iflavF1,iflavF3,isigF1,isigF3,isigA,
     -                       2),resultgaugeb(1,4),
     $                       BornCheck(4,iflavF1,iflavF3,isigF1,isigF3,
     -                       isigA,2)) !iinsA=4
                        BoxL(4,iflavF1,iflavF3,isigF1,isigF3,isigA,2)=
     $                       BoxL(4,iflavF1,iflavF3,isigF1,isigF3,isigA,
     -                       2)*GBF((-isigF3+3)/2,1,iflavF3)*
     -                       GBF((-isigF3+3)/2,2,iflavF3)
                        if(doGaugecheck)
     $                       call box_gaugeMM(accuracy,p(0,4),pA(0),
     -                       pIntBosl(0),p(0,3),scale(2),
     -                       resultgauge(1,4),resultgaugeb(1,4),
     $                       gaugecheck,passcheck(4,iflavF1,iflavF3,
     -                       isigF1,isigF3,isigA,2))
c$$$  BornCheck(4,iflavF1,iflavF3,isigF1,isigF3,isigA,2)=
c$$$  $                          BornCheck(4,iflavF1,iflavF3,isigF1,isigF3,isigA,2)*GBF((-isigF3+3)/2,1,iflavF3)*GBF((-isigF3+3)/2,2,iflavF3)
                        i=-1
                     enddo
                  enddo
               enddo
            enddo
            do iflav=3,4
               iflavF1=iflav
               iflavF2=7-iflav
               iflavF3=7-iflav
               iflavF4=iflav
               isigF1=1
               isigF3=1
               CALL BOXLINEMM(p(0,4),pA(0),pIntBosl(0),p(0,3),
     -              SpinorO(2+isigF3,3,isigF3),SpinorI(2-isigF3,4,
     -              isigF3),Photon(1,isigA),
     $              Curr21H(1,iflavF1,isigF1,3),-isigF3,scale(2),3,i,
     -              resultgauge(1,4),BoxL(4,iflavF1,iflavF3,isigF1,
     -              isigF3,isigA,3),
     $              resultgaugeb(1,4),BornCheck(4,iflavF1,iflavF3,
     -              isigF1,isigF3,isigA,3)) !iinsA=4
               BoxL(4,iflavF1,iflavF3,isigF1,isigF3,isigA,3)=
     $              BoxL(4,iflavF1,iflavF3,isigF1,isigF3,isigA,3)*
     -              GBF((-isigF3+3)/2,1,iflavF4)*GBF((-isigF3+3)/2,3,
     -              iflavF3)
               if(doGaugecheck)
     $              call box_gaugeMM(accuracy,p(0,4),pA(0),pIntBosl(0),
     -              p(0,3),scale(2),resultgauge(1,4),resultgaugeb(1,4),
     -              gaugecheck,
     $              passcheck(4,iflavF1,iflavF3,isigF1,isigF3,isigA,3))
c$$$  BornCheck(4,iflavF1,iflavF3,isigF1,isigF3,isigA,3)=
c$$$  $                 BornCheck(4,iflavF1,iflavF3,isigF1,isigF3,isigA,3)*GBF((-isigF3+3)/2,1,iflavF4)*GBF((-isigF3+3)/2,3,iflavF3)
            enddo
         enddo 

         
         i=1
         do isigA=-1,1,2  
            do isigF=1,4
               do iflavF1=3,4
                  do iflavF3=3,4
!Z contributions
                     if(resetamps)then
                        if(passcheck(1,iflavF1,iflavF3,iisig(1,isigF),
     -                       iisig(3,isigF),isigA,2).and.
     $                       passcheck(2,iflavF1,iflavF3,iisig(1,isigF),
     -                       iisig(3,isigF),isigA,2))
     $                       VirtUpper(iflavF1,iflavF3,isigF,isigA,2)=
     $                       (BoxL(1,iflavF1,iflavF3,iisig(1,isigF),
     -                       iisig(3,isigF),isigA,2)
     $                       +BoxL(2,iflavF1,iflavF3,iisig(1,isigF),
     -                       iisig(3,isigF),isigA,2))
     $                       -(SubAmpZ(iflavF1,iflavF3,isigF,isigA,1)+
     -                       SubAmpZ(iflavF1,iflavF3,isigF,isigA,2))*cte
                        if(passcheck(3,iflavF1,iflavF3,iisig(1,isigF),
     -                       iisig(3,isigF),isigA,2).and.
     $                       passcheck(4,iflavF1,iflavF3,iisig(1,isigF),
     -                       iisig(3,isigF),isigA,2))
     $                       VirtLower(iflavF1,iflavF3,isigF,isigA,2)=
     $                       (BoxL(3,iflavF1,iflavF3,iisig(1,isigF),
     -                       iisig(3,isigF),isigA,2)
     $                       +BoxL(4,iflavF1,iflavF3,iisig(1,isigF),
     -                       iisig(3,isigF),isigA,2))
     $                       -(SubAmpZ(iflavF1,iflavF3,isigF,isigA,3)+
     -                       SubAmpZ(iflavF1,iflavF3,isigF,isigA,4))*cte
                     else
                        VirtUpper(iflavF1,iflavF3,isigF,isigA,2)=
     $                       (BoxL(1,iflavF1,iflavF3,iisig(1,isigF),
     -                       iisig(3,isigF),isigA,2)
     $                       +BoxL(2,iflavF1,iflavF3,iisig(1,isigF),
     -                       iisig(3,isigF),isigA,2))
     $                       -(SubAmpZ(iflavF1,iflavF3,isigF,isigA,1)+
     -                       SubAmpZ(iflavF1,iflavF3,isigF,isigA,2))*cte
                        VirtLower(iflavF1,iflavF3,isigF,isigA,2)=
     $                       (BoxL(3,iflavF1,iflavF3,iisig(1,isigF),
     -                       iisig(3,isigF),isigA,2)
     $                       +BoxL(4,iflavF1,iflavF3,iisig(1,isigF),
     -                       iisig(3,isigF),isigA,2))
     $                       -(SubAmpZ(iflavF1,iflavF3,isigF,isigA,3)+
     -                       SubAmpZ(iflavF1,iflavF3,isigF,isigA,4))*cte
                     endif
                  enddo
               enddo
            enddo
            
            do iflav=3,4
               iflavF1=iflav
               iflavF2=7-iflav
               iflavF3=7-iflav
               iflavF4=iflav
               isigF=4
!W contributions
               if(resetamps)then
                  if(passcheck(1,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3).and.
     $                 passcheck(2,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3))
     $                 VirtUpper(iflavF1,iflavF3,isigF,isigA,3)=
     $                 (BoxL(1,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3)
     $                 +BoxL(2,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3))
     $                 -(SubAmpW(iflavF1,iflavF3,isigF,isigA,1)+SubAmpW(iflavF1,iflavF3,isigF,isigA,2))*cte
                  if(passcheck(3,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3).and.
     $                 passcheck(4,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3))               
     $                 VirtLower(iflavF1,iflavF3,isigF,isigA,3)=
     $                 (BoxL(3,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3)
     $                 +BoxL(4,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3))
     $                 -(SubAmpW(iflavF1,iflavF3,isigF,isigA,3)+SubAmpW(iflavF1,iflavF3,isigF,isigA,4))*cte
               else
                  VirtUpper(iflavF1,iflavF3,isigF,isigA,3)=
     $                 (BoxL(1,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3)
     $                 +BoxL(2,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3))
     $                 -(SubAmpW(iflavF1,iflavF3,isigF,isigA,1)+SubAmpW(iflavF1,iflavF3,isigF,isigA,2))*cte
                  VirtLower(iflavF1,iflavF3,isigF,isigA,3)=
     $                 (BoxL(3,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3)
     $                 +BoxL(4,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,3))
     $                 -(SubAmpW(iflavF1,iflavF3,isigF,isigA,3)+SubAmpW(iflavF1,iflavF3,isigF,isigA,4))*cte
               endif
            enddo
         enddo


c     isigF=4
         do isigA=-1,1,2
            resvirt(5)=resvirt(5)
     $           +2*dreal(AmpSumWW(4,3,isigF,isigA)*conjg(-alsUpperc2o4pi*VirtUpper(4,3,isigF,isigA,3)
     $           -alsLowerc2o4pi*VirtLower(4,3,isigF,isigA,3)+alsULcvirtc2o4pi*AmpSumWW(4,3,isigF,isigA)))
            resvirt(6)=resvirt(6)
     $           +2*dreal(AmpSumWW(3,4,isigF,isigA)*conjg(-alsUpperc2o4pi*VirtUpper(3,4,isigF,isigA,3)
     $           -alsLowerc2o4pi*VirtLower(3,4,isigF,isigA,3)+alsULcvirtc2o4pi*AmpSumWW(3,4,isigF,isigA)))
         enddo
      endif


*******************************************************************************
*******************************************************************************

      if(doGaugecheck)then
         do isigA=-1,1,2  
            do isigF=1,4
               do iflavF1=3,4
                  do iflavF3=3,4
                     do iinsA=1,4
                        do iBos=2,3
                           if(.not.passcheck(iinsA,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,iBos)) then
                              pspfails=.true.
c                              print*,sign(1),sign(2),sign(3),sign(4),iinsA,iflavF1,iflavF3,iisig(1,isigF),iisig(3,isigF),isigA,iBos
                           endif
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif

c$$$      print*,"NBoxUpper=",NboxUpper
c$$$      print*,"NBoxLower=",Nboxlower
c$$$      if(pspfails)then
c$$$         print*,"xU=",xU
c$$$         print*,"xL=",xL
c$$$         if(abs(xu-73d0).lt.1d0)print*,"NBoxUpper=",NboxUpper,"   NFailUpper=",NFailUpper
c$$$         if(abs(xl-73d0).lt.1d0)print*,"NBoxLower=",NboxLower,"   NFailLower=",NFailLower
c$$$      endif

      i=1
      do iflavF1=3,4
         do iflavF3=3,4
            do isigF=1,4
               do isigA=-1,1,2                     
                  resvirt(i)=resvirt(i)
     $                 +2*dreal(AmpSumZZ(iflavF1,iflavF3,isigF,isigA)*conjg(-alsUpperc2o4pi*VirtUpper(iflavF1,iflavF3,isigF,isigA,2)
     $                 -alsLowerc2o4pi*VirtLower(iflavF1,iflavF3,isigF,isigA,2)
     $                 +alsULcvirtc2o4pi*AmpSumZZ(iflavF1,iflavF3,isigF,isigA)))
c$$$                  if(abs(xu-73d0).lt.1d0)print*,"VirtUpper(",iflavF1,iflavF3,isigF,isigA,2,")",VirtUpper(iflavF1,iflavF3,isigF,isigA,2)
c$$$                  if(abs(xl-73d0).lt.1d0)print*,"VirtLower(",iflavF1,iflavF3,isigF,isigA,2,")",VirtLower(iflavF1,iflavF3,isigF,isigA,2)
               enddo
            enddo
            i=i+1
         enddo
      enddo
      
      do i=1,6
         
c            print*,"res(",i,")=",res(i)*fac,"+",resvirt(i)*fac

c         if(nlo.eq.0)then
            res(i)=res(i)+resvirt(i)
c         else
c            res(i)=resvirt(i)
c         endif
      enddo

      endif                     !end if(nlo)
#endif


      do i=1,6
         res(i)=res(i)*fac  
c         resvirt(i)=resvirt(i)*fac 
      enddo

      uucc = res(1) 
      uuss = res(2)
      ddcc = res(3)
      ddss = res(4)
      udsc = res(5)
      ducs = res(6)

      if (ldoblha) then
        if (blha_amptype(blha_curproc).eq.0) then ! tree
          blha_amp(1)=blha_amp(1)+
     &      res(blha_idsubproc(blha_cursubproc,blha_curproc))

        else if (blha_amptype(blha_curproc).eq.1) then ! loop
          do mu=0,3
            p21(mu) = pbar(mu,2)*sign(2) - pbar(mu,1)*sign(1)
            p43(mu) = pbar(mu,4)*sign(4) - pbar(mu,3)*sign(3)
          enddo
          p21(4) = p21(0)**2 - p21(1)**2 - p21(2)**2 - p21(3)**2
          p43(4) = p43(0)**2 - p43(1)**2 - p43(2)**2 - p43(3)**2
          lnmusqp1p2 = log(blha_scale**2/(-p21(4)))
          lnmusqp3p4 = log(blha_scale**2/(-p43(4)))

          blha_amp(1)=blha_amp(1)+
     &      blha_CF/(2d0*pi)*2d0*
     &       (als(1,1)+als(2,1))*(-1d0)*
     &       tree(blha_idsubproc(blha_cursubproc,blha_curproc))
          blha_amp(2)=blha_amp(2)+
     &      1d0/(2d0*pi)*2d0*
     &       ( als(1,1)*(-blha_gammaQuark-blha_CF*lnmusqp1p2)
     &        +als(2,1)*(-blha_gammaQuark-blha_CF*lnmusqp3p4) )*
     &       tree(blha_idsubproc(blha_cursubproc,blha_curproc))
c this is more complicated than necessary, 
c but should work for other VBF processes as well
          blha_amp(3)= blha_amp(3)+
     &      (res(blha_idsubproc(blha_cursubproc,blha_curproc))     ! tree + loop + I
     &       -tree(blha_idsubproc(blha_cursubproc,blha_curproc)))  ! -tree
     &       * blha_CF/c2                                          ! correct for Nc ?= inf
     &      + (-(als(1,1)+als(2,1))/2d0/pi*                        ! -I
     &          (blha_CF*cvirt+
     &            2d0*(blha_KQuark+blha_gammaQuark+
     &                 blha_tgammaQuark+(pi**2/6d0-1)*blha_CF) )
     &         +als(1,1)/pi*                                       ! +(mu^2/2p1.p2)^epsilon
     &         (-blha_gammaQuark*lnmusqp1p2-blha_CF/2d0*lnmusqp1p2**2)
     &         +als(2,1)/pi*                                       ! +(mu^2/2p3.p4)^epsilon
     &         (-blha_gammaQuark*lnmusqp3p4-blha_CF/2d0*lnmusqp3p4**2)
     &        )*tree(blha_idsubproc(blha_cursubproc,blha_curproc))
          blha_amp(4)= blha_amp(4)+
     &       tree(blha_idsubproc(blha_cursubproc,blha_curproc))

        else if (blha_amptype(blha_curproc).eq.2) then ! cctree
          ! 1-2
          call BLHA_cctree(1,2,
     &      -blha_CF*
     &       tree(blha_idsubproc(blha_cursubproc,blha_curproc))
     &    )
          ! 3-4
          call BLHA_cctree(3,4,
     &      -blha_CF*
     &       tree(blha_idsubproc(blha_cursubproc,blha_curproc))
     &    )

        else 
          call BLHA_amptypeerror(blha_amptype(blha_curproc),__FILE__,__LINE__)
        endif
      endif

      return
      end

*******************************************************************************
*******************************************************************************
*******************************************************************************
*******************************************************************************
*******************************************************************************


*******************************************************************************
*******************************************************************************
*******************************************************************************
*******************************************************************************
*******************************************************************************
      subroutine HAjj_subtr_ME(sign,kin,uucc,uuss,ddcc,ddss,udsc,ducs)

      IMPLICIT NONE

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/coupl.inc"

      integer sign(4+n_v),iisig(3,4),i,iflavF1,iflavF3,isigA,
     $     isigF,iinsA,kin
      double precision uucc,uuss,ddcc,ddss,udsc,ducs, res(6),fac
      
      double complex Higgs(6),Curr12(6,2:3,3:4,-1:1,2:3),
     $     Curr21(6,2:3,3:4,-1:1,2:3),Curr34(6,2:3,3:4,-1:1,2:3),
     $     Curr43(6,2:3,3:4,-1:1,2:3),Curr12A(6,2:3,3:4,-1:1,2,-1:1,2:3),
     $     Curr21A(6,2:3,3:4,-1:1,2,-1:1,2:3),Curr34A(6,2:3,3:4,-1:1,3:4,-1:1,2:3),
     $     Curr43A(6,2:3,3:4,-1:1,3:4,-1:1,2:3),Curr12tbv(6,2:3,3:4,-1:1,-1:1),Curr21tbv(6,2:3,3:4,-1:1,-1:1),
     $     Curr34tbv(6,2:3,3:4,-1:1,-1:1),Curr43tbv(6,2:3,3:4,-1:1,-1:1),
     $     SubAmpZ(3:4,3:4,4,-1:1,4),SubAmpW(3:4,3:4,4,-1:1,4),
     $     AmpSumZZ(3:4,3:4,4,-1:1),AmpSumWW(3:4,3:4,4,-1:1),SubAmpTBVu(3:4,3:4,4,-1:1),
     $     SubAmpTBVl(3:4,3:4,4,-1:1)


      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)

      COMMON /HAJJSUBTRCURRENTS/ Curr12,Curr21,Curr34,Curr43,Curr12A,Curr21A,Curr34A,
     $     Curr43A,Curr12tbv,Curr21tbv,Curr34tbv,Curr43tbv,Higgs,fac
      double precision c2,c2o2pi,c2o4pi
      parameter (c2=4d0/3d0, c2o2pi=c2/2d0/pi, c2o4pi=c2/4d0/pi)
      double precision dotrr
      EXTERNAL dotrr
      data iisig /-1,0,-1,  -1,0,1,   1,0,-1,   1,0,1/
      

c     reset amplitudes

      do iflavF1=3,4
         do iflavF3=3,4
            do isigF=1,4
               do isigA=-1,1,2
                  AmpSumZZ(iflavF1,iflavF3,isigF,isigA)=(0d0,0d0)
                  AmpSumWW(iflavF1,iflavF3,isigF,isigA)=(0d0,0d0)
               enddo
            enddo
         enddo
      enddo
      res(1)=0d0
      res(2)=0d0
      res(3)=0d0
      res(4)=0d0
      res(5)=0d0
      res(6)=0d0
*****************************************************************************
***************************************************************************** 
*****************************************************************************
      if(sign(1).eq.1.and.
c     $     sign(2).eq.1.and.
c     $     sign(3).eq.1.and.
     $     sign(4).eq.1)then
*****************************************************************************
         i=0
         do isigA=-1,1,2
c     Z contributions
            do iflavF1=3,4
               do iflavF3=3,4
                  do isigF=1,4
                     CALL VVSXXX(Curr12A(1,kin,iflavF1,iisig(1,isigF),1,isigA,2),Curr34(1,kin,iflavF3,iisig(3,isigF),2),Higgs,GZZH,
     $                    SubAmpZ(iflavF1,iflavF3,isigF,isigA,1)) !iinsA=1
                     CALL VVSXXX(Curr12A(1,kin,iflavF1,iisig(1,isigF),2,isigA,2),Curr34(1,kin,iflavF3,iisig(3,isigF),2),Higgs,GZZH,
     $                    SubAmpZ(iflavF1,iflavF3,isigF,isigA,2)) !iinsA=2
                     CALL VVSXXX(Curr12(1,kin,iflavF1,iisig(1,isigF),2),Curr34A(1,kin,iflavF3,iisig(3,isigF),3,isigA,2),Higgs,GZZH,
     $                    SubAmpZ(iflavF1,iflavF3,isigF,isigA,3)) !iinsA=3
                     CALL VVSXXX(Curr12(1,kin,iflavF1,iisig(1,isigF),2),Curr34A(1,kin,iflavF3,iisig(3,isigF),4,isigA,2),Higgs,GZZH,
     $                    SubAmpZ(iflavF1,iflavF3,isigF,isigA,4)) !iinsA=4
                  enddo
               enddo
c     W contributions
               iflavF3=7-iflavF1
               isigF=1    ! = (-)(-)
               CALL VVSXXX(Curr12A(1,kin,iflavF1,iisig(1,isigF),1,isigA,3),Curr34(1,kin,iflavF3,iisig(3,isigF),3),Higgs,GWWH,
     $              SubAmpW(iflavF1,iflavF3,isigF,isigA,1)) !iinsA=1
               CALL VVSXXX(Curr12A(1,kin,iflavF1,iisig(1,isigF),2,isigA,3),Curr34(1,kin,iflavF3,iisig(3,isigF),3),Higgs,GWWH,
     $              SubAmpW(iflavF1,iflavF3,isigF,isigA,2)) !iinsA=2
               CALL VVSXXX(Curr12(1,kin,iflavF1,iisig(1,isigF),3),Curr34A(1,kin,iflavF3,iisig(3,isigF),3,isigA,3),Higgs,GWWH,
     $              SubAmpW(iflavF1,iflavF3,isigF,isigA,3)) !iinsA=3
               CALL VVSXXX(Curr12(1,kin,iflavF1,iisig(1,isigF),3),Curr34A(1,kin,iflavF3,iisig(3,isigF),4,isigA,3),Higgs,GWWH,
     $              SubAmpW(iflavF1,iflavF3,isigF,isigA,4)) !iinsA=4
c     TBV contributions
               CALL VVSXXX(Curr12tbv(1,kin,iflavF1,iisig(1,isigF),isigA),Curr34(1,kin,iflavF3,iisig(3,isigF),3),Higgs,GWWH,
     $              SubAmpTBVu(iflavF1,iflavF3,isigF,isigA))
               CALL VVSXXX(Curr12(1,kin,iflavF1,iisig(1,isigF),3),Curr34tbv(1,kin,iflavF3,iisig(3,isigF),isigA),Higgs,GWWH,
     $              SubAmpTBVl(iflavF1,iflavF3,isigF,isigA))
               i=i+1
c               print*,"i=",i
               AmpSumWW(iflavF1,iflavF3,isigF,isigA)=SubAmpTBVu(iflavF1,iflavF3,isigF,isigA)+SubAmpTBVl(iflavF1,iflavF3,isigF,isigA)
               do iinsA=1,4
                  AmpSumWW(iflavF1,iflavF3,isigF,isigA)=AmpSumWW(iflavF1,iflavF3,isigF,isigA)
     $                 +SubAmpW(iflavF1,iflavF3,isigF,isigA,iinsA)
               enddo
            enddo
            res(5)=res(5)+dreal(AmpSumWW(3,4,isigF,isigA))**2+dimag(AmpSumWW(3,4,isigF,isigA))**2
            res(6)=res(6)+dreal(AmpSumWW(4,3,isigF,isigA))**2+dimag(AmpSumWW(4,3,isigF,isigA))**2
         enddo   
*****************************************************************************
      elseif(sign(1).eq.1.and.
c     $     sign(2).eq.1.and.
c     $     sign(3).eq.-1.and.
     $     sign(4).eq.-1)then
*****************************************************************************

         do isigA=-1,1,2
c     Z contributions
            do iflavF1=3,4
               do iflavF3=3,4
                  do isigF=1,4
                     CALL VVSXXX(Curr12A(1,kin,iflavF1,iisig(1,isigF),1,isigA,2),Curr43(1,kin,iflavF3,iisig(3,isigF),2),Higgs,GZZH,
     $                    SubAmpZ(iflavF1,iflavF3,isigF,isigA,1)) !iinsA=1
                     CALL VVSXXX(Curr12A(1,kin,iflavF1,iisig(1,isigF),2,isigA,2),Curr43(1,kin,iflavF3,iisig(3,isigF),2),Higgs,GZZH,
     $                    SubAmpZ(iflavF1,iflavF3,isigF,isigA,2)) !iinsA=2
                     CALL VVSXXX(Curr12(1,kin,iflavF1,iisig(1,isigF),2),Curr43A(1,kin,iflavF3,iisig(3,isigF),3,isigA,2),Higgs,GZZH,
     $                    SubAmpZ(iflavF1,iflavF3,isigF,isigA,3)) !iinsA=3
                     CALL VVSXXX(Curr12(1,kin,iflavF1,iisig(1,isigF),2),Curr43A(1,kin,iflavF3,iisig(3,isigF),4,isigA,2),Higgs,GZZH,
     $                    SubAmpZ(iflavF1,iflavF3,isigF,isigA,4)) !iinsA=4
                  enddo
               enddo
c     W contributions
               iflavF3=iflavF1
               isigF=2    ! = (-)(+)
               CALL VVSXXX(Curr12A(1,kin,iflavF1,iisig(1,isigF),1,isigA,3),Curr43(1,kin,iflavF3,iisig(3,isigF),3),Higgs,GWWH,
     $              SubAmpW(iflavF1,iflavF3,isigF,isigA,1)) !iinsA=1
               CALL VVSXXX(Curr12A(1,kin,iflavF1,iisig(1,isigF),2,isigA,3),Curr43(1,kin,iflavF3,iisig(3,isigF),3),Higgs,GWWH,
     $              SubAmpW(iflavF1,iflavF3,isigF,isigA,2)) !iinsA=2
               CALL VVSXXX(Curr12(1,kin,iflavF1,iisig(1,isigF),3),Curr43A(1,kin,iflavF3,iisig(3,isigF),3,isigA,3),Higgs,GWWH,
     $              SubAmpW(iflavF1,iflavF3,isigF,isigA,3)) !iinsA=3
               CALL VVSXXX(Curr12(1,kin,iflavF1,iisig(1,isigF),3),Curr43A(1,kin,iflavF3,iisig(3,isigF),4,isigA,3),Higgs,GWWH,
     $              SubAmpW(iflavF1,iflavF3,isigF,isigA,4)) !iinsA=4

c     TBV contributions
               CALL VVSXXX(Curr12tbv(1,kin,iflavF1,iisig(1,isigF),isigA),Curr43(1,kin,iflavF3,iisig(3,isigF),3),Higgs,GWWH,
     $              SubAmpTBVu(iflavF1,iflavF3,isigF,isigA))
               CALL VVSXXX(Curr12(1,kin,iflavF1,iisig(1,isigF),3),Curr43tbv(1,kin,iflavF3,iisig(3,isigF),isigA),Higgs,GWWH,
     $              SubAmpTBVl(iflavF1,iflavF3,isigF,isigA))
               AmpSumWW(iflavF1,iflavF3,isigF,isigA)=SubAmpTBVu(iflavF1,iflavF3,isigF,isigA)+SubAmpTBVl(iflavF1,iflavF3,isigF,isigA)
               do iinsA=1,4
                  AmpSumWW(iflavF1,iflavF3,isigF,isigA)=AmpSumWW(iflavF1,iflavF3,isigF,isigA)
     $                 +SubAmpW(iflavF1,iflavF3,isigF,isigA,iinsA)
               enddo
            enddo
            res(5)=res(5)+dreal(AmpSumWW(3,3,isigF,isigA))**2+dimag(AmpSumWW(3,3,isigF,isigA))**2
            res(6)=res(6)+dreal(AmpSumWW(4,4,isigF,isigA))**2+dimag(AmpSumWW(4,4,isigF,isigA))**2
         enddo
   
*****************************************************************************
      elseif(sign(1).eq.-1.and.
c     $     sign(2).eq.-1.and.
c     $     sign(3).eq.1.and.
     $     sign(4).eq.1)then
*****************************************************************************
         do isigA=-1,1,2
c     Z contributions
            do iflavF1=3,4
               do iflavF3=3,4
                  do isigF=1,4
                     CALL VVSXXX(Curr21A(1,kin,iflavF1,iisig(1,isigF),1,isigA,2),Curr34(1,kin,iflavF3,iisig(3,isigF),2),Higgs,GZZH,
     $                    SubAmpZ(iflavF1,iflavF3,isigF,isigA,1)) !iinsA=1
                     CALL VVSXXX(Curr21A(1,kin,iflavF1,iisig(1,isigF),2,isigA,2),Curr34(1,kin,iflavF3,iisig(3,isigF),2),Higgs,GZZH,
     $                    SubAmpZ(iflavF1,iflavF3,isigF,isigA,2)) !iinsA=2
                     CALL VVSXXX(Curr21(1,kin,iflavF1,iisig(1,isigF),2),Curr34A(1,kin,iflavF3,iisig(3,isigF),3,isigA,2),Higgs,GZZH,
     $                    SubAmpZ(iflavF1,iflavF3,isigF,isigA,3)) !iinsA=3
                     CALL VVSXXX(Curr21(1,kin,iflavF1,iisig(1,isigF),2),Curr34A(1,kin,iflavF3,iisig(3,isigF),4,isigA,2),Higgs,GZZH,
     $                    SubAmpZ(iflavF1,iflavF3,isigF,isigA,4)) !iinsA=4
                  enddo
               enddo
c     W contributions
               iflavF3=iflavF1
               isigF=3    ! = (+)(-)
               CALL VVSXXX(Curr21A(1,kin,iflavF1,iisig(1,isigF),1,isigA,3),Curr34(1,kin,iflavF3,iisig(3,isigF),3),Higgs,GWWH,
     $              SubAmpW(iflavF1,iflavF3,isigF,isigA,1)) !iinsA=1
               CALL VVSXXX(Curr21A(1,kin,iflavF1,iisig(1,isigF),2,isigA,3),Curr34(1,kin,iflavF3,iisig(3,isigF),3),Higgs,GWWH,
     $              SubAmpW(iflavF1,iflavF3,isigF,isigA,2)) !iinsA=2
               CALL VVSXXX(Curr21(1,kin,iflavF1,iisig(1,isigF),3),Curr34A(1,kin,iflavF3,iisig(3,isigF),3,isigA,3),Higgs,GWWH,
     $              SubAmpW(iflavF1,iflavF3,isigF,isigA,3)) !iinsA=3
               CALL VVSXXX(Curr21(1,kin,iflavF1,iisig(1,isigF),3),Curr34A(1,kin,iflavF3,iisig(3,isigF),4,isigA,3),Higgs,GWWH,
     $              SubAmpW(iflavF1,iflavF3,isigF,isigA,4)) !iinsA=4

c     TBV contributions
               CALL VVSXXX(Curr21tbv(1,kin,iflavF1,iisig(1,isigF),isigA),Curr34(1,kin,iflavF3,iisig(3,isigF),3),Higgs,GWWH,
     $              SubAmpTBVu(iflavF1,iflavF3,isigF,isigA))
               CALL VVSXXX(Curr21(1,kin,iflavF1,iisig(1,isigF),3),Curr34tbv(1,kin,iflavF3,iisig(3,isigF),isigA),Higgs,GWWH,
     $              SubAmpTBVl(iflavF1,iflavF3,isigF,isigA))

               AmpSumWW(iflavF1,iflavF3,isigF,isigA)=SubAmpTBVu(iflavF1,iflavF3,isigF,isigA)+SubAmpTBVl(iflavF1,iflavF3,isigF,isigA)

               do iinsA=1,4
                  AmpSumWW(iflavF1,iflavF3,isigF,isigA)=AmpSumWW(iflavF1,iflavF3,isigF,isigA)
     $                 +SubAmpW(iflavF1,iflavF3,isigF,isigA,iinsA)
               enddo
            enddo
            res(5)=res(5)+dreal(AmpSumWW(4,4,isigF,isigA))**2+dimag(AmpSumWW(4,4,isigF,isigA))**2
            res(6)=res(6)+dreal(AmpSumWW(3,3,isigF,isigA))**2+dimag(AmpSumWW(3,3,isigF,isigA))**2
         enddo
   
*****************************************************************************
      elseif(sign(1).eq.-1.and.
c     $     sign(2).eq.-1.and.
c     $     sign(3).eq.-1.and.
     $     sign(4).eq.-1)then
*****************************************************************************

         do isigA=-1,1,2
c     Z contributions
            do iflavF1=3,4
               do iflavF3=3,4
                  do isigF=1,4
                     CALL VVSXXX(Curr21A(1,kin,iflavF1,iisig(1,isigF),1,isigA,2),Curr43(1,kin,iflavF3,iisig(3,isigF),2),Higgs,GZZH,
     $                    SubAmpZ(iflavF1,iflavF3,isigF,isigA,1)) !iinsA=1
                     CALL VVSXXX(Curr21A(1,kin,iflavF1,iisig(1,isigF),2,isigA,2),Curr43(1,kin,iflavF3,iisig(3,isigF),2),Higgs,GZZH,
     $                    SubAmpZ(iflavF1,iflavF3,isigF,isigA,2)) !iinsA=2
                     CALL VVSXXX(Curr21(1,kin,iflavF1,iisig(1,isigF),2),Curr43A(1,kin,iflavF3,iisig(3,isigF),3,isigA,2),Higgs,GZZH,
     $                    SubAmpZ(iflavF1,iflavF3,isigF,isigA,3)) !iinsA=3
                     CALL VVSXXX(Curr21(1,kin,iflavF1,iisig(1,isigF),2),Curr43A(1,kin,iflavF3,iisig(3,isigF),4,isigA,2),Higgs,GZZH,
     $                    SubAmpZ(iflavF1,iflavF3,isigF,isigA,4)) !iinsA=4
                  enddo
               enddo
c     W contributions
               iflavF3=7-iflavF1
               isigF=4    ! = (-)(+)
               CALL VVSXXX(Curr21A(1,kin,iflavF1,iisig(1,isigF),1,isigA,3),Curr43(1,kin,iflavF3,iisig(3,isigF),3),Higgs,GWWH,
     $              SubAmpW(iflavF1,iflavF3,isigF,isigA,1)) !iinsA=1
               CALL VVSXXX(Curr21A(1,kin,iflavF1,iisig(1,isigF),2,isigA,3),Curr43(1,kin,iflavF3,iisig(3,isigF),3),Higgs,GWWH,
     $              SubAmpW(iflavF1,iflavF3,isigF,isigA,2)) !iinsA=2
               CALL VVSXXX(Curr21(1,kin,iflavF1,iisig(1,isigF),3),Curr43A(1,kin,iflavF3,iisig(3,isigF),3,isigA,3),Higgs,GWWH,
     $              SubAmpW(iflavF1,iflavF3,isigF,isigA,3)) !iinsA=3
               CALL VVSXXX(Curr21(1,kin,iflavF1,iisig(1,isigF),3),Curr43A(1,kin,iflavF3,iisig(3,isigF),4,isigA,3),Higgs,GWWH,
     $              SubAmpW(iflavF1,iflavF3,isigF,isigA,4)) !iinsA=4

c     TBV contributions
               CALL VVSXXX(Curr21tbv(1,kin,iflavF1,iisig(1,isigF),isigA),Curr43(1,kin,iflavF3,iisig(3,isigF),3),Higgs,GWWH,
     $              SubAmpTBVu(iflavF1,iflavF3,isigF,isigA))
               CALL VVSXXX(Curr21(1,kin,iflavF1,iisig(1,isigF),3),Curr43tbv(1,kin,iflavF3,iisig(3,isigF),isigA),Higgs,GWWH,
     $              SubAmpTBVl(iflavF1,iflavF3,isigF,isigA))

               AmpSumWW(iflavF1,iflavF3,isigF,isigA)=SubAmpTBVu(iflavF1,iflavF3,isigF,isigA)+SubAmpTBVl(iflavF1,iflavF3,isigF,isigA)

               do iinsA=1,4
                  AmpSumWW(iflavF1,iflavF3,isigF,isigA)=AmpSumWW(iflavF1,iflavF3,isigF,isigA)
     $                 +SubAmpW(iflavF1,iflavF3,isigF,isigA,iinsA)
               enddo
            enddo
            res(5)=res(5)+dreal(AmpSumWW(4,3,isigF,isigA))**2+dimag(AmpSumWW(4,3,isigF,isigA))**2
            res(6)=res(6)+dreal(AmpSumWW(3,4,isigF,isigA))**2+dimag(AmpSumWW(3,4,isigF,isigA))**2
         enddo
      endif
*******************************************************************************
*******************************************************************************

c     sum the NC contributions which interfere
      do iflavF1=3,4
         do iflavF3=3,4
            do isigF=1,4
               do isigA=-1,1,2
                  do iinsA=1,4
                     AmpSumZZ(iflavF1,iflavF3,isigF,isigA)=AmpSumZZ(iflavF1,iflavF3,isigF,isigA)
     $                 +SubAmpZ(iflavF1,iflavF3,isigF,isigA,iinsA)
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
                  res(i)=res(i)+dreal(AmpSumZZ(iflavF1,iflavF3,isigF,isigA))**2+dimag(AmpSumZZ(iflavF1,iflavF3,isigF,isigA))**2
               enddo
            enddo
            i=i+1
         enddo
      enddo
      
      do i=1,6
         res(i)=res(i)*fac  
      enddo
      uucc = res(1) 
      uuss = res(2)
      ddcc = res(3)
      ddss = res(4)
      udsc = res(5)
      ducs = res(6)

      return
      end


c$$$      subroutine box_gaugeMM(Accuracy,p1,p2,p3,p4,musq,
c$$$     -  result,resultb,td,WardB)
c$$$c
c$$$c Author: Francisco Campanario
c$$$C Date: 16/01/2009
c$$$C Compute the gauge  for the boxes
c$$$c When the ward test is not better than Accuraccy then WardB=.false.
c$$$c result: result(2) ==resultgauge(1) and resultgauge(2) out of boxline 
c$$$cvwith gaugetest=3
c$$$C resultb:resultb(2) ==resultgaugeb(1) and resultgaugeb(2) out of
c$$$c boxline with gaugetest=3
c$$$c p1,...p4, momenta of the boxline
c$$$C musq renormalization invariance. Same as used in boxline
c$$$c td: result of the gaugetest. For normal points returns ~1d-13 it is
c$$$c what the determines how good the point is and drives the result of
c$$$c WardB
c$$$c Accuracy: Accuracy demanded for the gauge test. With 1d-3 only 1/1d6 
c$$$c should be instable.
c$$$
c$$$      IMPLICIT NONE
c$$$      Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3)
c$$$      Real*8 musq,test1,test2
c$$$      Complex*16 cte
c$$$      Integer  i
c$$$      Complex*16 result(2)
c$$$      Complex*16 resultb(2)
c$$$      real*8 pi,pi2,test(3),td,tQ,Qt
c$$$      parameter (pi =3.14159265358979323846264338327950288d0)
c$$$      parameter (pi2 =9.86960440108935861883449099987615114d0)
c$$$      external dotrr,theta
c$$$      real*8 dotrr,s14,lm,ls,ths,theta
c$$$      Logical WardB
c$$$      double precision Accuracy
c$$$
c$$$      WardB=.False.
c$$$      td=1.d08
c$$$
c$$$
c$$$c        ctedieter = (3*pi*(0d0,1d0)-pi**2)
c$$$c        cte=-(3*pi*(0d0,1d0)-4.0d0*pi**2/3.0d0+7.0d0)
c$$$
c$$$      s14=2d0*dotrr(p1,p4)+dotrr(p1,p1)+dotrr(p4,p4)
c$$$      lm=Log(musq)
c$$$      ls=Log(Abs(s14))
c$$$      ths=theta(s14)
c$$$      
c$$$c       cte = -cte_dieter but general dependency in musq.
c$$$c         The sign is later on corrected in his code in
c$$$c          -dotcr(j,p)
c$$$      cte= -3d0*(Lm - Ls + (0d0,1d0)*pi*ths)
c$$$     $     - (Lm - Ls + (0d0,1d0)*pi*ths)**2
c$$$c     substracted the c_virt
c$$$     $     -7d0 + 1d0/3d0*pi2
c$$$c variables for gauge test:test1,test2
c$$$
c$$$      if(abs(resultb(1)).lt.1d-13) then
c$$$      test1=abs(result(1)- resultb(1)*(cte))
c$$$      else
c$$$      test1=abs(result(1)/(resultb(1)*(cte))-1d0)
c$$$      endif
c$$$
c$$$      if(abs(resultb(2)).lt.1d-13) then
c$$$       test2=abs(result(2)- (resultb(2)*(cte)))
c$$$      else
c$$$      test2=abs(result(2)/(resultb(2)*(cte))-1d0)
c$$$      endif
c$$$c$$$      print*,"test1",test1
c$$$c$$$      print*,"test2",test2
c$$$      td=max(test1,test2)
c$$$
c$$$
c$$$      IF (Td.LT.Accuracy) then
c$$$      WardB=.True.
c$$$      ENDIF
c$$$
c$$$! MOD:
c$$$c$$$      if(td.gt.1d-3.and.td.lt.1d-1)then
c$$$c$$$         WardB=.true.
c$$$c$$$         print*,"td=",td
c$$$c$$$      else
c$$$c$$$         WardB=.false.
c$$$c$$$      endif
c$$$
c$$$      Return
c$$$      end
