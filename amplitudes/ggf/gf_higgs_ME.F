***********************************************************************
      subroutine HiggsQQ (pbar,vsign,phbar,hsign,res12,res14,resfi) 
*----------------------------------------------------------------------
c... modified by: Michael Kubocz
c... date       : 17.10.2007
c... modified   : 30.01.2011
*----------------------------------------------------------------------

      IMPLICIT NONE

#include "ggf_process.inc"

c... Vertex functions
      complex*16 H2G_Loop,A2G_Loop,H2G_effective,A2G_effective


      real*8  ALFAS,XMT,ALFA,XMZ,XMW,SIN2W,XMH,gf,alfas_lo
      real*8 clr, xm2, xmg, b, v, a
      COMMON /BKOPIN/ ALFAS,XMT,ALFA,XMZ,XMW,SIN2W,XMH,gf,alfas_lo
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)

      real*8 pbar(0:3,4),phbar(0:4),p(0:3,4),p21(0:4),p43(0:4),
     &       p41(0:4),p23(0:4)


      complex*16 j21(0:3,-1:1),j43(0:3,-1:1),j41(0:3,-1:1),
     &     j23(0:3,-1:1),psi(2,-1:1,4) 

      real*8 pi
      parameter (pi=3.14159 26535 89793238d0)

      complex*16 THjjTopA,THjjTopB,THjjBottomA,THjjBottomB,
     &           TAjjTopA,TAjjTopB,TAjjBottomA,TAjjBottomB,
     &           THeffA,THeffB,TAeffA,TAeffB,
     &           TsHjjTopA,TsHjjTopB,TsHjjBottomA,TsHjjBottomB,
     &           TbHjjTopA,TbHjjTopB,TbHjjBottomA,TbHjjBottomB

      complex*16 BsHjjStopA(2),BsHjjStopB(2),TsHjjStopA(2),
     &           TsHjjStopB(2),BsHjjSbottomA(2),BsHjjSbottomB(2),
     &           TsHjjSbottomA(2),TsHjjSbottomB(2),
     &           BbHjjStopA(2),BbHjjStopB(2),TbHjjStopA(2),
     &           TbHjjStopB(2),BbHjjSbottomA(2),BbHjjSbottomB(2),
     &           TbHjjSbottomA(2),TbHjjSbottomB(2)

c... other local variables 
      integer i,mu,isig1,isig3,icall
      data icall /0/

c... variables for the main part of the program
      integer vsign(4),hsign
      real*8  res12, res14, resfi  
     
      complex*16 mm, mmf
    
      real*8 fac, vev
      save icall,fac

      real*8 dotrr
c$$$      EXTERNAL dotrr,H2G_Loop,A2G_Loop,H2G_effective,A2G_effective,
c$$$     &         H1G_SqLoop,H2G_SqLoop
      EXTERNAL dotrr,H2G_Loop,A2G_Loop,H2G_effective,A2G_effective
      
      if (icall.eq.0) then
          vev = 2d0*xmw/(v(1,3)*dsqrt(8d0))
c... common overall factor
          fac=2d0*(4d0*alfas**2/vev)**2   
          icall = 1
      endif
     
c... define the internal momenta
      do mu = 0,3
         do i = 1,4
            p(mu,i) = pbar(mu,i)*vsign(i)
         enddo
         p21(mu) = p(mu,2)-p(mu,1) ! virtual gluon -q1
         p43(mu) = p(mu,4)-p(mu,3) ! virtual gluon -q2
         p41(mu) = p(mu,4)-p(mu,1) ! virtual gluon -q1 for res14
         p23(mu) = p(mu,2)-p(mu,3) ! virtual gluon -q2 for res14
      enddo
      p21(4) = dotrr(p21,p21)
      p43(4) = dotrr(p43,p43)
      p41(4) = dotrr(p41,p41)
      p23(4) = dotrr(p23,p23)

c... get the external spinors (including factor sqrt(2E) )
      call psi0m(4,pbar(0,1),vsign(1),psi)

c... get the f-fbar currents J21^mu, J43^mu
      call curr(1,psi(1,-1,2),psi(1,-1,1),j21)
      call curr(1,psi(1,-1,4),psi(1,-1,3),j43)
      call curr(1,psi(1,-1,4),psi(1,-1,1),j41)
      call curr(1,psi(1,-1,2),psi(1,-1,3),j23)

c... initialize amplitudes with zero values
      THjjTopA    = DCMPLX(0d0,0d0)
      THjjBottomA = DCMPLX(0d0,0d0)
      TAjjTopA    = DCMPLX(0d0,0d0)
      TAjjBottomA = DCMPLX(0d0,0d0)
      THeffA      = DCMPLX(0d0,0d0)
      TAeffA      = DCMPLX(0d0,0d0)

      TsHjjTopA   = DCMPLX(0d0,0d0)
      TsHjjBottomA= DCMPLX(0d0,0d0)

      TbHjjTopA   = DCMPLX(0d0,0d0)
      TbHjjBottomA= DCMPLX(0d0,0d0)

      THjjTopB    = DCMPLX(0d0,0d0)
      THjjBottomB = DCMPLX(0d0,0d0)
      TAjjTopB    = DCMPLX(0d0,0d0)
      TAjjBottomB = DCMPLX(0d0,0d0)
      THeffB      = DCMPLX(0d0,0d0)
      TAeffB      = DCMPLX(0d0,0d0)
      
      TsHjjTopB   = DCMPLX(0d0,0d0)
      TsHjjBottomB= DCMPLX(0d0,0d0)

      TbHjjTopB   = DCMPLX(0d0,0d0)
      TbHjjBottomB= DCMPLX(0d0,0d0)

      do i=1,2
         BsHjjStopA(i)    = DCMPLX(0d0,0d0)
         TsHjjStopA(i)    = DCMPLX(0d0,0d0)
         BsHjjSbottomA(i) = DCMPLX(0d0,0d0)
         TsHjjSbottomA(i) = DCMPLX(0d0,0d0)
         BsHjjStopB(i)    = DCMPLX(0d0,0d0)
         TsHjjStopB(i)    = DCMPLX(0d0,0d0)
         BsHjjSbottomB(i) = DCMPLX(0d0,0d0)
         TsHjjSbottomB(i) = DCMPLX(0d0,0d0)

         BbHjjStopA(i)    = DCMPLX(0d0,0d0)
         TbHjjStopA(i)    = DCMPLX(0d0,0d0)
         BbHjjSbottomA(i) = DCMPLX(0d0,0d0)
         TbHjjSbottomA(i) = DCMPLX(0d0,0d0)
         BbHjjStopB(i)    = DCMPLX(0d0,0d0)
         TbHjjStopB(i)    = DCMPLX(0d0,0d0)
         BbHjjSbottomB(i) = DCMPLX(0d0,0d0)
         TbHjjSbottomB(i) = DCMPLX(0d0,0d0)
      enddo

c-----------------------------------------------
c... start: loop over the quark polarization 
c-----------------------------------------------
      res12 = 0.0d0
      res14 = 0.0d0
      resfi = 0.0d0
      do isig1 = -1,1,2
         do isig3 = -1,1,2

c     mm:                     mmf:
c     p1 -->--*-->-- p2        p1 -->--*-->-- p4 
c             |                        |
c             O-- H / A0               O-- H / A0
c             |                        |
c     p3 -->--*-->-- p4        p3 -->--*-->-- p2

c... cp-even Higgs vertex in effective limes
            if (Heff) then
               THeffa=cplModHeff
     &         *(H2G_effective(p21,p43,J21(0,isig1),J43(0,isig3)))
               THeffb=cplModHeff
     &         *(H2G_effective(p23,p41,J23(0,isig1),J41(0,isig3)))
            endif
c... cp-odd Higgs vertex in effective limes
            if (Aeff) then
               TAeffa=cplModAeff
     &         *(A2G_effective(p21,p43,J21(0,isig1),J43(0,isig3)))
               TAeffb=cplModAeff
     &         *(A2G_effective(p23,p41,J23(0,isig1),J41(0,isig3)))
            endif
c... cp-even Higgs vertex with top-loop (SM)
            if (HloopTop) then
               THjjTopa=cplModHtop
     &         *(H2G_Loop(p21,p43,J21(0,isig1),J43(0,isig3),m_top))
               THjjTopb=cplModHtop
     &         *(H2G_Loop(p23,p41,J23(0,isig1),J41(0,isig3),m_top))
            endif
c... cp-even Higgs vertex with bottom-loop (SM)
            if (HloopBottom) then   
               THjjBottoma=cplModHbottom
     &         *(H2G_Loop(p21,p43,J21(0,isig1),J43(0,isig3),m_bottom))
               THjjBottomb=cplModHbottom
     &         *(H2G_Loop(p23,p41,J23(0,isig1),J41(0,isig3),m_bottom))
            endif
c... cp-odd Higgs vertex with top-loop (MSSM)
            if (AloopTop) then
               TAjjTopa=cplModAtop
     &         *(A2G_Loop(p21,p43,J21(0,isig1),J43(0,isig3),m_top))
               TAjjTopb=cplModAtop
     &         *(A2G_Loop(p23,p41,J23(0,isig1),J41(0,isig3),m_top))
            endif
c... cp-odd Higgs vertex with bottom-loop (MSSM)
            if (AloopBottom) then
               TAjjBottoma=cplModAbottom
     &         *(A2G_Loop(p21,p43,J21(0,isig1),J43(0,isig3),m_bottom))
               TAjjBottomb=cplModAbottom
     &         *(A2G_Loop(p23,p41,J23(0,isig1),J41(0,isig3),m_bottom))
            endif
c... cp-even Higgs h vertex with top-loop (MSSM)
            if (sHloopTop) then
               TsHjjTopa=cplModsHtop
     &         *(H2G_Loop(p21,p43,J21(0,isig1),J43(0,isig3),m_top))
               TsHjjTopb=cplModsHtop
     &         *(H2G_Loop(p23,p41,J23(0,isig1),J41(0,isig3),m_top))
            endif
c... cp-even Higgs h vertex with bottom-loop (MSSM)
            if (sHloopBottom) then   
               TsHjjBottoma=cplModsHbottom
     &         *(H2G_Loop(p21,p43,J21(0,isig1),J43(0,isig3),m_bottom))
               TsHjjBottomb=cplModsHbottom
     &         *(H2G_Loop(p23,p41,J23(0,isig1),J41(0,isig3),m_bottom))
            endif
c... cp-even Higgs H vertex with top-loop (MSSM)
            if (bHloopTop) then
               TbHjjTopa=cplModbHtop
     &         *(H2G_Loop(p21,p43,J21(0,isig1),J43(0,isig3),m_top))
               TbHjjTopb=cplModbHtop
     &         *(H2G_Loop(p23,p41,J23(0,isig1),J41(0,isig3),m_top))
            endif
c... cp-even Higgs H vertex with bottom-loop (MSSM)
            if (bHloopBottom) then   
               TbHjjBottoma=cplModbHbottom
     &         *(H2G_Loop(p21,p43,J21(0,isig1),J43(0,isig3),m_bottom))
               TbHjjBottomb=cplModbHbottom
     &         *(H2G_Loop(p23,p41,J23(0,isig1),J41(0,isig3),m_bottom))
            endif
c... cp-even h^0 vertex with stop1-loop (MSSM)
c$$$            if (sHloopStop1) then
c$$$               TsHjjStopA(1)=cplModsHstop(1)*(H2G_SqLoop(p21,p43,
c$$$     &              J21(0,isig1),J43(0,isig3),m_stop(1)))
c$$$               TsHjjStopB(1)=cplModsHstop(1)*(H2G_SqLoop(p23,p41,
c$$$     &              J23(0,isig1),J41(0,isig3),m_stop(1)))
c$$$               BsHjjStopA(1)=cplModsHstop(1)*(H1G_SqLoop(p21,p43,
c$$$     &              J21(0,isig1),J43(0,isig3),m_stop(1)))
c$$$               BsHjjStopB(1)=cplModsHstop(1)*(H1G_SqLoop(p23,p41,
c$$$     &              J23(0,isig1),J41(0,isig3),m_stop(1)))
c$$$            endif
c... cp-even h^0 vertex with stop2-loop (MSSM)
c$$$            if (sHloopStop2) then
c$$$               TsHjjStopA(2)=cplModsHstop(2)*(H2G_SqLoop(p21,p43,
c$$$     &              J21(0,isig1),J43(0,isig3),m_stop(2)))
c$$$               TsHjjStopB(2)=cplModsHstop(2)*(H2G_SqLoop(p23,p41,
c$$$     &              J23(0,isig1),J41(0,isig3),m_stop(2)))
c$$$               BsHjjStopA(2)=cplModsHstop(2)*(H1G_SqLoop(p21,p43,
c$$$     &              J21(0,isig1),J43(0,isig3),m_stop(2)))
c$$$               BsHjjStopB(2)=cplModsHstop(2)*(H1G_SqLoop(p23,p41,
c$$$     &              J23(0,isig1),J41(0,isig3),m_stop(2)))
c$$$            endif
c... cp-even h^0 vertex with sbottom1-loop (MSSM)
c$$$            if (sHloopSbottom1) then
c$$$               TsHjjSbottomA(1)=cplModsHsbottom(1)*(H2G_SqLoop(p21,  
c$$$     &              p43,J21(0,isig1),J43(0,isig3),m_sbottom(1)))
c$$$               TsHjjSbottomB(1)=cplModsHsbottom(1)*(H2G_SqLoop(p23,
c$$$     &              p41,J23(0,isig1),J41(0,isig3),m_sbottom(1)))
c$$$               BsHjjSbottomA(1)=cplModsHsbottom(1)*(H1G_SqLoop(p21,
c$$$     &              p43,J21(0,isig1),J43(0,isig3),m_sbottom(1)))
c$$$               BsHjjSbottomB(1)=cplModsHsbottom(1)*(H1G_SqLoop(p23,
c$$$     &              p41,J23(0,isig1),J41(0,isig3),m_sbottom(1)))
c$$$            endif
c... cp-even h^0 vertex with sbottom2-loop (MSSM)
c$$$            if (sHloopSbottom2) then
c$$$               TsHjjSbottomA(2)=cplModsHsbottom(2)*(H2G_SqLoop(p21,
c$$$     &              p43,J21(0,isig1),J43(0,isig3),m_sbottom(2)))
c$$$               TsHjjSbottomB(2)=cplModsHsbottom(2)*(H2G_SqLoop(p23,
c$$$     &              p41,J23(0,isig1),J41(0,isig3),m_sbottom(2)))
c$$$               BsHjjSbottomA(2)=cplModsHsbottom(2)*(H1G_SqLoop(p21,
c$$$     &              p43,J21(0,isig1),J43(0,isig3),m_sbottom(2)))
c$$$               BsHjjSbottomB(2)=cplModsHsbottom(2)*(H1G_SqLoop(p23,
c$$$     &              p41,J23(0,isig1),J41(0,isig3),m_sbottom(2)))
c$$$            endif
c... cp-even H^0 vertex with stop1-loop (MSSM)
c$$$            if (bHloopStop1) then
c$$$               TbHjjStopA(1)=cplModbHstop(1)*(H2G_SqLoop(p21,p43,
c$$$     &              J21(0,isig1),J43(0,isig3),m_stop(1)))
c$$$               TbHjjStopB(1)=cplModbHstop(1)*(H2G_SqLoop(p23,p41,
c$$$     &              J23(0,isig1),J41(0,isig3),m_stop(1)))
c$$$               BbHjjStopA(1)=cplModbHstop(1)*(H1G_SqLoop(p21,p43,
c$$$     &              J21(0,isig1),J43(0,isig3),m_stop(1)))
c$$$               BbHjjStopB(1)=cplModbHstop(1)*(H1G_SqLoop(p23,p41,
c$$$     &              J23(0,isig1),J41(0,isig3),m_stop(1)))
c$$$            endif
c$$$c... cp-even H^0 vertex with stop2-loop (MSSM)
c$$$            if (bHloopStop2) then
c$$$               TbHjjStopA(2)=cplModbHstop(2)*(H2G_SqLoop(p21,p43,
c$$$     &              J21(0,isig1),J43(0,isig3),m_stop(2)))
c$$$               TbHjjStopB(2)=cplModbHstop(2)*(H2G_SqLoop(p23,p41,
c$$$     &              J23(0,isig1),J41(0,isig3),m_stop(2)))
c$$$               BbHjjStopA(2)=cplModbHstop(2)*(H1G_SqLoop(p21,p43,
c$$$     &              J21(0,isig1),J43(0,isig3),m_stop(2)))
c$$$               BbHjjStopB(2)=cplModbHstop(2)*(H1G_SqLoop(p23,p41,
c$$$     &              J23(0,isig1),J41(0,isig3),m_stop(2)))
c$$$            endif
c... cp-even H^0 vertex with sbottom1-loop (MSSM)
c$$$            if (bHloopSbottom1) then
c$$$               TbHjjSbottomA(1)=cplModbHsbottom(1)*(H2G_SqLoop(p21,
c$$$     &              p43,J21(0,isig1),J43(0,isig3),m_sbottom(1)))
c$$$               TbHjjSbottomB(1)=cplModbHsbottom(1)*(H2G_SqLoop(p23,
c$$$     &              p41,J23(0,isig1),J41(0,isig3),m_sbottom(1)))
c$$$               BbHjjSbottomA(1)=cplModbHsbottom(1)*(H1G_SqLoop(p21,
c$$$     &              p43,J21(0,isig1),J43(0,isig3),m_sbottom(1)))
c$$$               BbHjjSbottomB(1)=cplModbHsbottom(1)*(H1G_SqLoop(p23,
c$$$     &              p41,J23(0,isig1),J41(0,isig3),m_sbottom(1)))
c$$$            endif
c... cp-even H^0 vertex with sbottom2-loop (MSSM)
c$$$            if (bHloopSbottom2) then
c$$$               TbHjjSbottomA(2)=cplModbHsbottom(2)*(H2G_SqLoop(p21,
c$$$     &              p43,J21(0,isig1),J43(0,isig3),m_sbottom(2)))
c$$$               TbHjjSbottomB(2)=cplModbHsbottom(2)*(H2G_SqLoop(p23,
c$$$     &              p41,J23(0,isig1),J41(0,isig3),m_sbottom(2)))
c$$$               BbHjjSbottomA(2)=cplModbHsbottom(2)*(H1G_SqLoop(p21,
c$$$     &              p43,J21(0,isig1),J43(0,isig3),m_sbottom(2)))
c$$$               BbHjjSbottomB(2)=cplModbHsbottom(2)*(H1G_SqLoop(p23,
c$$$     &              p41,J23(0,isig1),J41(0,isig3),m_sbottom(2)))
c$$$            endif

c... calculation of squared amplitude
            mm=(THjjTopA+THjjBottomA+TAjjTopA+TAjjBottomA+THeffA+TAeffA
     &         +TsHjjTopA+TsHjjBottomA+TbHjjTopA+TbHjjBottomA
     &         +BsHjjStopA(1)+TsHjjStopA(1)+BsHjjStopA(2)+TsHjjStopA(2)
     &         +BsHjjSbottomA(1)+TsHjjSbottomA(1)+BsHjjSbottomA(2)
     &         +TsHjjSbottomA(2)
     &         +BbHjjStopA(1)+TbHjjStopA(1)+BbHjjStopA(2)+TbHjjStopA(2)
     &         +BbHjjSbottomA(1)+TbHjjSbottomA(1)+BbHjjSbottomA(2)
     &         +TbHjjSbottomA(2))
     &         /(p21(4)*p43(4))

            res12 = res12 + dreal(mm)**2+dimag(mm)**2 

            mmf=(THjjTopB+THjjBottomB+TAjjTopB+TAjjBottomB+THeffB+TAeffB
     &          +TsHjjTopB+TsHjjBottomB+TbHjjTopB+TbHjjBottomB
     &          +BsHjjStopB(1)+TsHjjStopB(1)+BsHjjStopB(2)+TsHjjStopB(2)
     &          +BsHjjSbottomB(1)+TsHjjSbottomB(1)+BsHjjSbottomB(2)
     &          +TsHjjSbottomB(2)
     &          +BbHjjStopB(1)+TbHjjStopB(1)+BbHjjStopB(2)+TbHjjStopB(2)
     &          +BbHjjSbottomB(1)+TbHjjSbottomB(1)+BbHjjSbottomB(2)
     &          +TbHjjSbottomB(2))
     &          /(p41(4)*p23(4))

            res14 = res14 + dreal(mmf)**2+dimag(mmf)**2 

            if (isig1.eq.isig3) then
               resfi = resfi + 0.6666666666666666666666d0*
     &                (dreal(mm)*dreal(mmf)+dimag(mm)*dimag(mmf))
            endif
         enddo
      enddo
c-----------------------------------------------
c... end: loop over the quark polarization 
c-----------------------------------------------

c... multiply by the coupling factor, which includes factor 2 for color
      resfi = (res12+res14+resfi)*fac
      res12 = res12*fac
      res14 = res14*fac

      END

***********************************************************************
      SUBROUTINE HiggsQG(PBAR,FSIGN,FPOL,KBAR,GSIGN,PHBAR,HSIGN,RH)
*----------------------------------------------------------------------
c... modified by: Michael Kubocz
c... date       : 17.10.2007
c... modified   : 30.01.2011
*----------------------------------------------------------------------
c... HQQGG_mod calculates the color-summed amplitude square of the Hqqgg
c... amplitude for different Hgg-couplings and fixed fermion helicities. 
c... All external fermions are assumed to be massless. Gluon 
c... polarizations are summed over. 

c... pbar(0:3,2)    external (physical) quark momenta
c... fsign(2)       quark(+1) vs antiquark (-1)
c... fpol(2)        physical helicity of 2 quarks: +-1

c... kbar(0:3,2)    external momenta of 2 gluons (massless)
c... gsign(2)       sign factor for gluons: -1 = initial
c...                                        +1 = final state gluon

c... phbar          Higgs boson momentum (assumed on-shell)
c... hsign            +1 for final state Higgs, -1 for decaying Higgs

c... OUTPUT:
c... rh(0:2)        resulting polarisation and color summed amplitude
c...                square for the two large N color flows:
c...                  rh(1) = q1--g1--g2--q2
c...                  rh(2) = q1--g2--g1--q2
c...                  rh(0) = rh(1)+rh(2)
***********************************************************************

      IMPLICIT NONE

#include "ggf_process.inc"

c... vertex functions      
      complex*16 H2G_Loop,A2G_Loop,H3G_Loop,A3G_Loop

      complex*16 H2G_effective,A2G_effective,H3G_effective,A3G_effective,
     &           GGG_vertex

c... information on coupling constants from KOPPEL
      real*8 ALFAS,XMT,ALFA,XMZ,XMW,SIN2W,XMH, gf, alfas_lo
      real*8 clr, xm2, xmg, b, v, a
      COMMON /BKOPIN/ ALFAS,XMT,ALFA,XMZ,XMW,SIN2W,XMH, gf, alfas_lo
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)
c$$$      complex*16 amplit(8,3,2)
c$$$      COMMON  /AMPS/ amplit

      integer NF,NG,NGPL1,NFPL1
      real*8 CPLUS, CMINUS
      parameter (NF=2,NG=2,NGPL1=NG+2,NFPL1=NF)    
      parameter (CPLUS=2.3333333333333d0, CMINUS=3.0d0)

      real*8 PBAR(0:3,NF),KBAR(0:3,NG),PHBAR(0:3)
      real*8 K(0:4,NG),P(0:4,NFPL1),PH(0:3),Q(0:4,NGPL1)
      real*8 K21(0:3),K_H(0:3,2)
      real*8 x(0:3),y(0:3),q1q2,k1k2,k1q1,k2q1

      complex*16 psi(2,-1:1,NF),j_21(0:3,-1:1),j21(0:3),Jgg(0:3)
      complex*16 BRAKET(2,NF,NG)
      complex*16 epsg(0:3,2),epsh(0:3,2)


      real*8 vev,fac

      save fac

      complex*16 ampAtop(8),ampAbottom(8),ampHtop(8),ampHbottom(8),
     &           ampAeff(8),ampHeff(8),m12,m21,
     &           ampsHtop(8),ampsHbottom(8),ampbHtop(8),ampbHbottom(8),
     &           ampsHstop1(15),ampsHstop2(15), 
     &           ampsHsbottom1(15),ampsHsbottom2(15),
     &           ampbHstop1(15),ampbHstop2(15), 
     &           ampbHsbottom1(15),ampbHsbottom2(15)

      real*8 pi
      parameter (pi=3.14159 26535 89793d0)
      integer FSIGN(4), FPOL(NF), SIGMAF(NF), GSIGN(4), HSIGN
      real*8 rh(0:3) 
*segorig      integer FSIGN(NF), FPOL(NF), SIGMAF(NF), GSIGN(NG), HSIGN, CPL
*segorig      real*8 rh(0:2) 
c... further local variables
      integer mu, sigma, i1, i, i1pol, i2pol
      real*8 res12, res21, resfi
      real*8 detx, detxt, dummy(0:4),srtdet
      real*8 col11_22, col12_21

      logical ldebug
      integer npolg1, ienter
c... default is 2, npolg1=3 allows to check gauge cancellations 
      parameter (ldebug = .false., npolg1=2)  
      data ienter /0/
      save ienter

c... external functions
      complex*16  S1C
      real*8 dotrr,srt
c      EXTERNAL S1C,dotrr,srt,
c     &         H2G_Loop,A2G_Loop,H3G_Loop,A3G_Loop,
c     &         H1G_SqLoop,H2G_SqLoop,H3G_SqLoop,
c     &         H2G_effective,A2G_effective,H3G_effective,A3G_effective,
c     &         GGG_vertex
      EXTERNAL S1C,dotrr,srt,
     &         H2G_Loop,A2G_Loop,H3G_Loop,A3G_Loop,
     &         H2G_effective,A2G_effective,H3G_effective,A3G_effective,
     &         GGG_vertex

      col11_22=16d0/3d0
      col12_21=-2d0/3d0
      if (ienter.eq.0) then
         ienter = ienter+1
         vev = 2d0*xmw/(v(1,3)*dsqrt(8d0))
c... common overall factor
         fac=(4d0*alfas**2/vev)**2   
      endif

      do i = 0,2
         rh(i) = 0d0
      enddo
      res12=0d0
      res21=0d0
      resfi=0d0

c... assign internal momentum variables
      SIGMAF(1) = FSIGN(1)*FPOL(1)
      SIGMAF(2) = FSIGN(2)*FPOL(2)
      IF (SIGMAF(1).NE.SIGMAF(2)) THEN
        WRITE(*,*) ' Helicity nonconserving amplitude called in HQQGG '
        RETURN
      ENDIF
      SIGMA = SIGMAF(1)

c... the K1,K2 (gluon momenta) are outgoing, not ingoing !!!!
c... this is necessary for using the bra-ket subroutines
      DO MU = 0,3
         P(MU,1) = FSIGN(1)*PBAR(MU,1) 
         P(MU,2) = FSIGN(2)*PBAR(MU,2)
         K(MU,1) = GSIGN(1)*KBAR(MU,1)
         K(MU,2) = GSIGN(2)*KBAR(MU,2)
         PH(MU)  = hsign*phbar(mu)
         Q(MU,1) = P(MU,2)-P(MU,1)
         Q(MU,2) = K(MU,1)+K(MU,2)
         Q(MU,3) = PH(MU)+K(MU,2)
         Q(MU,4) = PH(MU)+K(MU,1)
      ENDDO

      if (ldebug) then
         print *,'P1 :',(P(mu,1),mu=0,3)
         print *,'P2 :',(P(mu,2),mu=0,3)
         print *,'K1 :',(K(mu,1),mu=0,3)
         print *,'K2 :',(K(mu,2),mu=0,3)
         print *,'PH :',(PH(mu),mu=0,3)
         print *,'Q1 :',(Q(mu,1),mu=0,3)
         print *,'Q2 :',(Q(mu,2),mu=0,3)
         print *,'Q3 :',(Q(mu,3),mu=0,3)
         print *,'Q4 :',(Q(mu,4),mu=0,3)
         print *,'mom. cons. (ph+K1+K2+P2-P1) :',
     &        ((Q(MU,1)+Q(MU,2)+PH(MU)),mu=0,3)
      endif
      P(4,1)=0D0
      P(4,2)=0D0
      K(4,1)=0D0
      K(4,2)=0D0
      DO i1=1,4
         Q(4,i1)=dotrr(Q(0,i1),Q(0,i1))
      ENDDO
c... get the external spinors (including factor sqrt(2E))
      call psi0m(2,pbar(0,1),fsign(1),psi)

c... get the f-fbar current J21^mu
      call curr(1,psi(1,-1,2),psi(1,-1,1),j_21)
      do mu=0,3
         j21(mu)=j_21(mu,SIGMA)/Q(4,1)
      enddo

c... calculate some scalar products
      q1q2= dotrr(Q(0,1),Q(0,2))
      k1k2= dotrr(k(0,1),k(0,2))
      k1q1= dotrr(k(0,1),q(0,1))
      k2q1= dotrr(q(0,1),k(0,2))
      if (ldebug) then
         print *,'q1q2 : ',q1q2
         print *,'k1k2 : ',k1k2
         print *,'k1q1 : ',k1q1
         print *,'k2q1 : ',k2q1
      endif

c... calculate basis vectors for gluon polarization
      do mu = 0,3
         x(mu) = K(mu,1)*k2q1 + K(mu,2)*k1q1 - Q(mu,1)*k1k2
      enddo
      call EPSRRR(y,k(0,1),K(0,2),q(0,1))
      detxt = Q(4,1)*k1k2-2*k1q1*k2q1
      detx = k1k2*detxt
      srtdet=srt(-detx)
      if (ldebug) then
         print*,' q1q2 = ',q1q2
         print*,' detx = ',detx,' detx/k1.k2 = ', detxt
         print*,' x.x  = ',dotrr(x,x)
         print*,' y.y  = ',dotrr(y,y)
         print*,' detxt/(|k1k2*q1^2|+2|k1q1*k2q1|) ',
     &        detxt/(abs(k1k2*q(4,1))+2*abs(k1q1*k2q1))
         print*,'detxt/k^6 = ', detxt/(k(0,1)**6)
      endif

c... initialize amplitudes with zero values
      do i=1,8
         ampAeff(i)    = DCMPLX(0d0,0d0)
         ampHeff(i)    = DCMPLX(0d0,0d0)  

         ampAtop(i)    = DCMPLX(0d0,0d0)
         ampAbottom(i) = DCMPLX(0d0,0d0)
         ampHtop(i)    = DCMPLX(0d0,0d0)
         ampHbottom(i) = DCMPLX(0d0,0d0)

         ampsHtop(i)    = DCMPLX(0d0,0d0)
         ampsHbottom(i) = DCMPLX(0d0,0d0)
         ampbHtop(i)    = DCMPLX(0d0,0d0)
         ampbHbottom(i) = DCMPLX(0d0,0d0)
      enddo
      do i=1,15
         ampsHstop1(i) = DCMPLX(0d0,0d0)
         ampsHstop2(i) = DCMPLX(0d0,0d0)
         ampsHsbottom1(i)= DCMPLX(0d0,0d0)
         ampsHsbottom2(i)= DCMPLX(0d0,0d0)

         ampbHstop1(i) = DCMPLX(0d0,0d0)
         ampbHstop2(i) = DCMPLX(0d0,0d0)
         ampbHsbottom1(i)= DCMPLX(0d0,0d0)
         ampbHsbottom2(i)= DCMPLX(0d0,0d0)
      enddo
      
c-----------------------------------------------
c... start: sum over gluon polarizations
c-----------------------------------------------
      do i1pol= 1,npolg1
         do i2pol = 1,2
c... calculate gluon-polarizations : there are 3 possibilities 
c... 1st one : the cartesian basis
c            call POLVEC(kbar(0,2),i2pol,ep(0,2))
c            if (i1pol.ne.3) call POLVEC(kbar(0,1),i1pol,ep(0,1))
c            do mu=0,3
c               epsg(mu,1)=ep(mu,1)
c               epsg(mu,2)=ep(mu,2)
c            enddo

c...   2nd one : the helicity basis
            call HELVEC(kbar(0,2),gsign(2),(i2pol*2-3),EPsg(0,2))
            if (i1pol.ne.3) call HELVEC(kbar(0,1),gsign(1),(i1pol*2-3),
     &                                  EPsg(0,1))
            if (i1pol.eq.3) then
               do mu = 0,3
                  epsg(mu,1) = K(mu,1)/k(0,1)
               enddo
            endif

c...  3rd one : the basis of eqn. 3.21 - 3.28
c...  (here the basis vectors from above are needed)
c            if (i1pol.eq.1 .and. i2pol.eq.1) then !the xx case
c               do mu = 0,3
c                  epsg(mu,1) = x(mu)/srtdet
c                  epsg(mu,2) = x(mu)/srtdet
c               enddo
c            elseif( i1pol.eq.1 .and. i2pol.eq.2) then !the xy case
c               do mu = 0,3
c                  epsg(mu,1) = x(mu)/srtdet
c                  epsg(mu,2) = y(mu)/srtdet
c               enddo
c            elseif( i1pol.eq.2 .and. i2pol.eq.1) then !the yx case
c               do mu = 0,3
c                  epsg(mu,1) = y(mu)/srtdet
c                  epsg(mu,2) = x(mu)/srtdet
c               enddo
c            elseif( i1pol.eq.2 .and. i2pol.eq.2) then !the yy case
c               do mu = 0,3
c                  epsg(mu,1) = y(mu)/srtdet
c                  epsg(mu,2) = y(mu)/srtdet
c               enddo
c            elseif( i1pol.eq.3 .and. i2pol.eq.1) then !the kx case
c               do mu = 0,3
c                  epsg(mu,1) = K(mu,1)/k(0,1)
c                  epsg(mu,2) = x(mu)/srtdet
c               enddo
c            elseif( i1pol.eq.3 .and. i2pol.eq.2) then !the ky case
c               do mu = 0,3
c                  epsg(mu,1) = K(mu,1)/K(0,1)
c                  epsg(mu,2) = y(mu)/srtdet
c               enddo
c            endif

c... calculate bra & kets
            CALL KET2C(psi(1,SIGMA,1),.false.,P(0,1),SIGMA,K(0,1),
     &                 EPSG(0,1),BRAKET(1,1,1),DUMMY)
            CALL KET2C(psi(1,SIGMA,1),.false.,P(0,1),SIGMA,K(0,2),
     &                 EPSG(0,2),BRAKET(1,1,2),DUMMY)
            CALL BRA2C(psi(1,SIGMA,2),.false.,P(0,2),SIGMA,K(0,1),
     &                 EPSG(0,1),BRAKET(1,2,1),DUMMY)
            CALL BRA2C(psi(1,SIGMA,2),.false.,P(0,2),SIGMA,K(0,2),
     &                 EPSG(0,2),BRAKET(1,2,2),DUMMY)

c... calculate the two-gluon current for amlitude 'amp1'
            CALL GGG_current(K(0,2),epsg(0,2),K(0,1),epsg(0,1),K21,Jgg)

c... calculate eps_iH and then the subamplitudes      
c... cp-even Higgs vertex in effective limes    
            if (Heff) then
               CALL EpsH_effective(K(0,1),.false.,epsg(0,1),ph,K_H(0,1),
     &                             epsh(0,1))
               CALL EpsH_effective(K(0,2),.false.,epsg(0,2),ph,K_H(0,1),
     &                             epsh(0,2))
               ampHeff(1)=cplModHeff
     &         *H2G_effective(q(0,1),q(0,2),J21,Jgg)
               ampHeff(2)=cplModHeff
     &         *GGG_vertex(K(0,1),Q(0,3),Q(0,1),epsg(0,1),epsh(0,2),J21)
               ampHeff(3)=cplModHeff
     &         *GGG_vertex(K(0,2),Q(0,1),Q(0,4),epsg(0,2),J21,epsh(0,1))
               ampHeff(4)=cplModHeff
     &         *S1C(psi(1,SIGMA,2),EpsH(0,2),.TRUE.,SIGMA,BRAKET(1,1,1))
               ampHeff(5)=cplModHeff
     &         *S1C(BRAKET(1,2,2),EpsH(0,1),.TRUE.,SIGMA,psi(1,sigma,1))
               ampHeff(6)=cplModHeff
     &         *S1C(BRAKET(1,2,1),EpsH(0,2),.TRUE.,SIGMA,psi(1,sigma,1))
               ampHeff(7)=cplModHeff
     &         *S1C(psi(1,SIGMA,2),EpsH(0,1),.TRUE.,SIGMA,BRAKET(1,1,2))
               ampHeff(8)=cplModHeff*H3G_effective(K(0,2),Q(0,1),K(0,1),
     &                               epsg(0,2),J21,epsg(0,1))
            endif
c... cp-odd Higgs vertex in effective limes
            if (Aeff) then
               CALL EpsA_effective(K(0,1),epsg(0,1),ph,K_H(0,1),
     &                             epsh(0,1))
               CALL EpsA_effective(K(0,2),epsg(0,2),ph,K_H(0,1),
     &                             epsh(0,2))
               ampAeff(1)=cplModAeff
     &         *A2G_effective(q(0,1),q(0,2),J21,Jgg)
               ampAeff(2)=cplModAeff
     &         *GGG_vertex(K(0,1),Q(0,3),Q(0,1),epsg(0,1),epsh(0,2),J21)
               ampAeff(3)=cplModAeff
     &         *GGG_vertex(K(0,2),Q(0,1),Q(0,4),epsg(0,2),J21,epsh(0,1))
               ampAeff(4)=cplModAeff
     &         *S1C(psi(1,SIGMA,2),EpsH(0,2),.TRUE.,SIGMA,BRAKET(1,1,1))
               ampAeff(5)=cplModAeff
     &         *S1C(BRAKET(1,2,2),EpsH(0,1),.TRUE.,SIGMA,psi(1,sigma,1))
               ampAeff(6)=cplModAeff
     &         *S1C(BRAKET(1,2,1),EpsH(0,2),.TRUE.,SIGMA,psi(1,sigma,1))
               ampAeff(7)=cplModAeff
     &         *S1C(psi(1,SIGMA,2),EpsH(0,1),.TRUE.,SIGMA,BRAKET(1,1,2))
               ampAeff(8)=cplModAeff
     &         *A3G_effective(PH,epsg(0,2),J21,epsg(0,1))
            endif
c... cp-even Higgs vertex with top-loop (SM)
            if (HloopTop) then
               CALL EpsH_Loop(K(0,1),.false.,epsg(0,1),ph,K_H(0,1),
     &                        epsh(0,1),m_top)
               CALL EpsH_Loop(K(0,2),.false.,epsg(0,2),ph,K_H(0,1),
     &                        epsh(0,2),m_top)
               ampHtop(1)=cplModHtop
     &         *H2G_Loop(q(0,1),q(0,2),J21,Jgg,m_top)
               ampHtop(2)=cplModHtop
     &         *GGG_vertex(K(0,1),Q(0,3),Q(0,1),epsg(0,1),epsh(0,2),J21)
               ampHtop(3)=cplModHtop
     &         *GGG_vertex(K(0,2),Q(0,1),Q(0,4),epsg(0,2),J21,epsh(0,1))
               ampHtop(4)=cplModHtop
     &         *S1C(psi(1,SIGMA,2),EpsH(0,2),.TRUE.,SIGMA,BRAKET(1,1,1))
               ampHtop(5)=cplModHtop
     &         *S1C(BRAKET(1,2,2),EpsH(0,1),.TRUE.,SIGMA,psi(1,sigma,1))
               ampHtop(6)=cplModHtop
     &         *S1C(BRAKET(1,2,1),EpsH(0,2),.TRUE.,SIGMA,psi(1,sigma,1))
               ampHtop(7)=cplModHtop
     &         *S1C(psi(1,SIGMA,2),EpsH(0,1),.TRUE.,SIGMA,BRAKET(1,1,2))
               ampHtop(8)=cplModHtop*H3G_Loop(K(0,2),Q(0,1),K(0,1),
     &              epsg(0,2),J21,epsg(0,1),m_top)
            endif
c... cp-even Higgs vertex with bottom-loop (SM)
            if (HloopBottom) then
               CALL EpsH_Loop(K(0,1),.false.,epsg(0,1),ph,K_H(0,1),
     &                        epsh(0,1),m_bottom)
               CALL EpsH_Loop(K(0,2),.false.,epsg(0,2),ph,K_H(0,1),
     &                        epsh(0,2),m_bottom)
               ampHbottom(1)=cplModHbottom
     &         *H2G_Loop(q(0,1),q(0,2),J21,Jgg,m_bottom)
               ampHbottom(2)=cplModHbottom
     &         *GGG_vertex(K(0,1),Q(0,3),Q(0,1),epsg(0,1),epsh(0,2),J21)
               ampHbottom(3)=cplModHbottom
     &         *GGG_vertex(K(0,2),Q(0,1),Q(0,4),epsg(0,2),J21,epsh(0,1))
               ampHbottom(4)=cplModHbottom
     &         *S1C(psi(1,SIGMA,2),EpsH(0,2),.TRUE.,SIGMA,BRAKET(1,1,1))
               ampHbottom(5)=cplModHbottom
     &         *S1C(BRAKET(1,2,2),EpsH(0,1),.TRUE.,SIGMA,psi(1,sigma,1))
               ampHbottom(6)=cplModHbottom
     &         *S1C(BRAKET(1,2,1),EpsH(0,2),.TRUE.,SIGMA,psi(1,sigma,1))
               ampHbottom(7)=cplModHbottom
     &         *S1C(psi(1,SIGMA,2),EpsH(0,1),.TRUE.,SIGMA,BRAKET(1,1,2))
               ampHbottom(8)=cplModHbottom*H3G_Loop(K(0,2),Q(0,1),
     &              K(0,1),epsg(0,2),J21,epsg(0,1),m_bottom)
            endif
c... cp-odd Higgs vertex with top-loop (MSSM)                     
            if (AloopTop) then
               CALL EpsA_Loop(K(0,1),epsg(0,1),ph,K_H(0,1),epsh(0,1),
     &                        m_top)
               CALL EpsA_Loop(K(0,2),epsg(0,2),ph,K_H(0,1),epsh(0,2),
     &                        m_top)
               ampAtop(1)=cplModAtop
     &         *A2G_Loop(Q(0,1),Q(0,2),J21,Jgg,m_top) 
               ampAtop(2)=cplModAtop
     &         *GGG_vertex(K(0,1),Q(0,3),Q(0,1),epsg(0,1),epsh(0,2),J21)
               ampAtop(3)=cplModAtop
     &         *GGG_vertex(K(0,2),Q(0,1),Q(0,4),epsg(0,2),J21,epsh(0,1))
               ampAtop(4)=cplModAtop
     &         *S1C(psi(1,SIGMA,2),EpsH(0,2),.TRUE.,SIGMA,BRAKET(1,1,1))
               ampAtop(5)=cplModAtop
     &         *S1C(BRAKET(1,2,2),EpsH(0,1),.TRUE.,SIGMA,psi(1,sigma,1))
               ampAtop(6)=cplModAtop
     &         *S1C(BRAKET(1,2,1),EpsH(0,2),.TRUE.,SIGMA,psi(1,sigma,1))
               ampAtop(7)=cplModAtop
     &         *S1C(psi(1,SIGMA,2),EpsH(0,1),.TRUE.,SIGMA,BRAKET(1,1,2))
               ampAtop(8)=cplModAtop*A3G_Loop(K(0,1),K(0,2),Q(0,1),
     &              epsg(0,1),epsg(0,2),J21(0),m_top)  
            endif
c... cp-odd Higgs vertex with bottom-loop (MSSM)
            if (AloopBottom) then
               CALL EpsA_Loop(K(0,1),epsg(0,1),ph,K_H(0,1),epsh(0,1),
     &                        m_bottom)
               CALL EpsA_Loop(K(0,2),epsg(0,2),ph,K_H(0,1),epsh(0,2),
     &                        m_bottom)
               ampAbottom(1)=cplModAbottom
     &         *A2G_Loop(Q(0,1),Q(0,2),J21,Jgg,m_bottom) !7
               ampAbottom(2)=cplModAbottom
     &         *GGG_vertex(K(0,1),Q(0,3),Q(0,1),epsg(0,1),epsh(0,2),J21)
               ampAbottom(3)=cplModAbottom
     &         *GGG_vertex(K(0,2),Q(0,1),Q(0,4),epsg(0,2),J21,epsh(0,1))
               ampAbottom(4)=cplModAbottom
     &         *S1C(psi(1,SIGMA,2),EpsH(0,2),.TRUE.,SIGMA,BRAKET(1,1,1))
               ampAbottom(5)=cplModAbottom
     &         *S1C(BRAKET(1,2,2),EpsH(0,1),.TRUE.,SIGMA,psi(1,sigma,1))
               ampAbottom(6)=cplModAbottom
     &         *S1C(BRAKET(1,2,1),EpsH(0,2),.TRUE.,SIGMA,psi(1,sigma,1))
               ampAbottom(7)=cplModAbottom
     &         *S1C(psi(1,SIGMA,2),EpsH(0,1),.TRUE.,SIGMA,BRAKET(1,1,2))
               ampAbottom(8)=cplModAbottom*A3G_Loop(K(0,1),K(0,2),
     &              Q(0,1),epsg(0,1),epsg(0,2),J21(0),m_bottom)  
            endif
c... cp-even Higgs h vertex with top-loop (MSSM)
            if (sHloopTop) then
               CALL EpsH_Loop(K(0,1),.false.,epsg(0,1),ph,K_H(0,1),
     &                        epsh(0,1),m_top)
               CALL EpsH_Loop(K(0,2),.false.,epsg(0,2),ph,K_H(0,1),
     &                        epsh(0,2),m_top)
               ampsHtop(1)=cplModsHtop
     &         *H2G_Loop(q(0,1),q(0,2),J21,Jgg,m_top)
               ampsHtop(2)=cplModsHtop
     &         *GGG_vertex(K(0,1),Q(0,3),Q(0,1),epsg(0,1),epsh(0,2),J21)
               ampsHtop(3)=cplModsHtop
     &         *GGG_vertex(K(0,2),Q(0,1),Q(0,4),epsg(0,2),J21,epsh(0,1))
               ampsHtop(4)=cplModsHtop
     &         *S1C(psi(1,SIGMA,2),EpsH(0,2),.TRUE.,SIGMA,BRAKET(1,1,1))
               ampsHtop(5)=cplModsHtop
     &         *S1C(BRAKET(1,2,2),EpsH(0,1),.TRUE.,SIGMA,psi(1,sigma,1))
               ampsHtop(6)=cplModsHtop
     &         *S1C(BRAKET(1,2,1),EpsH(0,2),.TRUE.,SIGMA,psi(1,sigma,1))
               ampsHtop(7)=cplModsHtop
     &         *S1C(psi(1,SIGMA,2),EpsH(0,1),.TRUE.,SIGMA,BRAKET(1,1,2))
               ampsHtop(8)=cplModsHtop*H3G_Loop(K(0,2),Q(0,1),K(0,1),
     &              epsg(0,2),J21,epsg(0,1),m_top)
             endif
c... cp-even Higgs vertex with bottom-loop (SM)
            if (sHloopBottom) then
               CALL EpsH_Loop(K(0,1),.false.,epsg(0,1),ph,K_H(0,1),
     &                        epsh(0,1),m_bottom)
               CALL EpsH_Loop(K(0,2),.false.,epsg(0,2),ph,K_H(0,1),
     &                        epsh(0,2),m_bottom)
               ampsHbottom(1)=cplModsHbottom
     &         *H2G_Loop(q(0,1),q(0,2),J21,Jgg,m_bottom)
               ampsHbottom(2)=cplModsHbottom
     &         *GGG_vertex(K(0,1),Q(0,3),Q(0,1),epsg(0,1),epsh(0,2),J21)
               ampsHbottom(3)=cplModsHbottom
     &         *GGG_vertex(K(0,2),Q(0,1),Q(0,4),epsg(0,2),J21,epsh(0,1))
               ampsHbottom(4)=cplModsHbottom
     &         *S1C(psi(1,SIGMA,2),EpsH(0,2),.TRUE.,SIGMA,BRAKET(1,1,1))
               ampsHbottom(5)=cplModsHbottom
     &         *S1C(BRAKET(1,2,2),EpsH(0,1),.TRUE.,SIGMA,psi(1,sigma,1))
               ampsHbottom(6)=cplModsHbottom
     &         *S1C(BRAKET(1,2,1),EpsH(0,2),.TRUE.,SIGMA,psi(1,sigma,1))
               ampsHbottom(7)=cplModsHbottom
     &         *S1C(psi(1,SIGMA,2),EpsH(0,1),.TRUE.,SIGMA,BRAKET(1,1,2))
               ampsHbottom(8)=cplModsHbottom*H3G_Loop(K(0,2),Q(0,1),
     &              K(0,1),epsg(0,2),J21,epsg(0,1),m_bottom)
            endif
c... cp-even Higgs H vertex with top-loop (MSSM)
            if (bHloopTop) then
               CALL EpsH_Loop(K(0,1),.false.,epsg(0,1),ph,K_H(0,1),
     &                        epsh(0,1),m_top)
               CALL EpsH_Loop(K(0,2),.false.,epsg(0,2),ph,K_H(0,1),
     &                        epsh(0,2),m_top)
               ampbHtop(1)=cplModbHtop
     &         *H2G_Loop(q(0,1),q(0,2),J21,Jgg,m_top)
               ampbHtop(2)=cplModbHtop
     &         *GGG_vertex(K(0,1),Q(0,3),Q(0,1),epsg(0,1),epsh(0,2),J21)
               ampbHtop(3)=cplModbHtop
     &         *GGG_vertex(K(0,2),Q(0,1),Q(0,4),epsg(0,2),J21,epsh(0,1))
               ampbHtop(4)=cplModbHtop
     &         *S1C(psi(1,SIGMA,2),EpsH(0,2),.TRUE.,SIGMA,BRAKET(1,1,1))
               ampbHtop(5)=cplModbHtop
     &         *S1C(BRAKET(1,2,2),EpsH(0,1),.TRUE.,SIGMA,psi(1,sigma,1))
               ampbHtop(6)=cplModbHtop
     &         *S1C(BRAKET(1,2,1),EpsH(0,2),.TRUE.,SIGMA,psi(1,sigma,1))
               ampbHtop(7)=cplModbHtop
     &         *S1C(psi(1,SIGMA,2),EpsH(0,1),.TRUE.,SIGMA,BRAKET(1,1,2))
               ampbHtop(8)=cplModbHtop*H3G_Loop(K(0,2),Q(0,1),K(0,1),
     &              epsg(0,2),J21,epsg(0,1),m_top)
             endif
c... cp-even Higgs H vertex with bottom-loop (MSSM)
            if (bHloopBottom) then
               CALL EpsH_Loop(K(0,1),.false.,epsg(0,1),ph,K_H(0,1),
     &                        epsh(0,1),m_bottom)
               CALL EpsH_Loop(K(0,2),.false.,epsg(0,2),ph,K_H(0,1),
     &                        epsh(0,2),m_bottom)
               ampbHbottom(1)=cplModbHbottom
     &         *H2G_Loop(q(0,1),q(0,2),J21,Jgg,m_bottom)
               ampbHbottom(2)=cplModbHbottom
     &         *GGG_vertex(K(0,1),Q(0,3),Q(0,1),epsg(0,1),epsh(0,2),J21)
               ampbHbottom(3)=cplModbHbottom
     &         *GGG_vertex(K(0,2),Q(0,1),Q(0,4),epsg(0,2),J21,epsh(0,1))
               ampbHbottom(4)=cplModbHbottom
     &         *S1C(psi(1,SIGMA,2),EpsH(0,2),.TRUE.,SIGMA,BRAKET(1,1,1))
               ampbHbottom(5)=cplModbHbottom
     &         *S1C(BRAKET(1,2,2),EpsH(0,1),.TRUE.,SIGMA,psi(1,sigma,1))
               ampbHbottom(6)=cplModbHbottom
     &         *S1C(BRAKET(1,2,1),EpsH(0,2),.TRUE.,SIGMA,psi(1,sigma,1))
               ampbHbottom(7)=cplModbHbottom
     &         *S1C(psi(1,SIGMA,2),EpsH(0,1),.TRUE.,SIGMA,BRAKET(1,1,2))
               ampbHbottom(8)=cplModbHbottom*H3G_Loop(K(0,2),Q(0,1),
     &              K(0,1),epsg(0,2),J21,epsg(0,1),m_bottom)
            endif
c... cp-even h^0 vertex with stop1-loop (MSSM)
c$$$            if (sHloopStop1) then
c$$$               CALL EpsH_SqTriLoop(K(0,1),epsg(0,1),ph,K_H(0,1),
c$$$     &              epsh(0,1),m_stop(1))
c$$$               CALL EpsH_SqTriLoop(K(0,2),epsg(0,2),ph,K_H(0,1),
c$$$     &              epsh(0,2),m_stop(1))
c$$$               ampsHstop1(1)=cplModsHstop(1)
c$$$     &         *H2G_SqLoop(Q(0,1),Q(0,2),J21,Jgg,m_stop(1))
c$$$               ampsHstop1(2)=cplModsHstop(1)*GGG_vertex(K(0,1),
c$$$     &              Q(0,3),Q(0,1),epsg(0,1),epsh(0,2),J21)
c$$$               ampsHstop1(3)=cplModsHstop(1)*GGG_vertex(K(0,2),
c$$$     &              Q(0,1),Q(0,4),epsg(0,2),J21,epsh(0,1))
c$$$               ampsHstop1(4)=cplModsHstop(1)*S1C(psi(1,SIGMA,2),
c$$$     &              EpsH(0,2),.TRUE.,SIGMA,BRAKET(1,1,1))
c$$$               ampsHstop1(5)=cplModsHstop(1)*S1C(BRAKET(1,2,2),
c$$$     &              EpsH(0,1),.TRUE.,SIGMA,psi(1,sigma,1))
c$$$               ampsHstop1(6)=cplModsHstop(1)*S1C(BRAKET(1,2,1),
c$$$     &              EpsH(0,2),.TRUE.,SIGMA,psi(1,sigma,1))
c$$$               ampsHstop1(7)=cplModsHstop(1)*S1C(psi(1,SIGMA,2),
c$$$     &              EpsH(0,1),.TRUE.,SIGMA,BRAKET(1,1,2))
c$$$               ampsHstop1(8)=cplModsHstop(1)*H3G_SqLoop(K(0,1),
c$$$     &              K(0,2),Q(0,1),epsg(0,1),epsg(0,2),J21(0),
c$$$     &              m_stop(1))
c$$$               CALL Eps_SqBlobLoop(K(0,1),epsg(0,1),ph,K_H(0,1),
c$$$     &              epsh(0,1),m_stop(1))
c$$$               CALL Eps_SqBlobLoop(K(0,2),epsg(0,2),ph,K_H(0,1),
c$$$     &              epsh(0,2),m_stop(1))
c$$$               ampsHstop1(9)=cplModsHstop(1)
c$$$     &         *H1G_SqLoop(Q(0,1),Q(0,2),J21,Jgg,m_stop(1))               
c$$$               ampsHstop1(10)=cplModsHstop(1)*GGG_vertex(K(0,1),
c$$$     &              Q(0,3),Q(0,1),epsg(0,1),epsh(0,2),J21)
c$$$               ampsHstop1(11)=cplModsHstop(1)*GGG_vertex(K(0,2),
c$$$     &              Q(0,1),Q(0,4),epsg(0,2),J21,epsh(0,1))
c$$$               ampsHstop1(12)=cplModsHstop(1)*S1C(psi(1,SIGMA,2),
c$$$     &              EpsH(0,2),.TRUE.,SIGMA,BRAKET(1,1,1)) 
c$$$               ampsHstop1(13)=cplModsHstop(1)*S1C(BRAKET(1,2,2),
c$$$     &              EpsH(0,1),.TRUE.,SIGMA,psi(1,sigma,1))
c$$$               ampsHstop1(14)=cplModsHstop(1)*S1C(BRAKET(1,2,1),
c$$$     &              EpsH(0,2),.TRUE.,SIGMA,psi(1,sigma,1))
c$$$               ampsHstop1(15)=cplModsHstop(1)*S1C(psi(1,SIGMA,2),
c$$$     &              EpsH(0,1),.TRUE.,SIGMA,BRAKET(1,1,2)) 
c$$$            endif
c... cp-even h^0 vertex with stop2-loop (MSSM)
c$$$           if (sHloopStop2) then
c$$$               CALL EpsH_SqTriLoop(K(0,1),epsg(0,1),ph,K_H(0,1),
c$$$     &             epsh(0,1),m_stop(2))
c$$$               CALL EpsH_SqTriLoop(K(0,2),epsg(0,2),ph,K_H(0,1),
c$$$     &              epsh(0,2),m_stop(2))
c$$$               ampsHstop2(1)=cplModsHstop(2)
c$$$     &         *H2G_SqLoop(Q(0,1),Q(0,2),J21,Jgg,m_stop(2))
c$$$               ampsHstop2(2)=cplModsHstop(2)*GGG_vertex(K(0,1),
c$$$     &              Q(0,3),Q(0,1),epsg(0,1),epsh(0,2),J21)
c$$$               ampsHstop2(3)=cplModsHstop(2)*GGG_vertex(K(0,2),
c$$$     &              Q(0,1),Q(0,4),epsg(0,2),J21,epsh(0,1))
c$$$               ampsHstop2(4)=cplModsHstop(2)*S1C(psi(1,SIGMA,2),
c$$$     &              EpsH(0,2),.TRUE.,SIGMA,BRAKET(1,1,1))
c$$$               ampsHstop2(5)=cplModsHstop(2)*S1C(BRAKET(1,2,2),
c$$$     &              EpsH(0,1),.TRUE.,SIGMA,psi(1,sigma,1))
c$$$               ampsHstop2(6)=cplModsHstop(2)*S1C(BRAKET(1,2,1),
c$$$     &              EpsH(0,2),.TRUE.,SIGMA,psi(1,sigma,1))
c$$$               ampsHstop2(7)=cplModsHstop(2)*S1C(psi(1,SIGMA,2),
c$$$     &              EpsH(0,1),.TRUE.,SIGMA,BRAKET(1,1,2))
c$$$               ampsHstop2(8)=cplModsHstop(2)*H3G_SqLoop(K(0,1),
c$$$     &              K(0,2),Q(0,1),epsg(0,1),epsg(0,2),J21(0),
c$$$     &              m_stop(2))
c$$$               CALL Eps_SqBlobLoop(K(0,1),epsg(0,1),ph,K_H(0,1),
c$$$     &              epsh(0,1),m_stop(2))
c$$$               CALL Eps_SqBlobLoop(K(0,2),epsg(0,2),ph,K_H(0,1),
c$$$     &              epsh(0,2),m_stop(2))
c$$$               ampsHstop2(9)=cplModsHstop(2)
c$$$     &         *H1G_SqLoop(Q(0,1),Q(0,2),J21,Jgg,m_stop(2))               
c$$$               ampsHstop2(10)=cplModsHstop(2)*GGG_vertex(K(0,1),
c$$$     &              Q(0,3),Q(0,1),epsg(0,1),epsh(0,2),J21)
c$$$               ampsHstop2(11)=cplModsHstop(2)*GGG_vertex(K(0,2),
c$$$     &              Q(0,1),Q(0,4),epsg(0,2),J21,epsh(0,1))
c$$$               ampsHstop2(12)=cplModsHstop(2)*S1C(psi(1,SIGMA,2),
c$$$     &              EpsH(0,2),.TRUE.,SIGMA,BRAKET(1,1,1)) 
c$$$               ampsHstop2(13)=cplModsHstop(2)*S1C(BRAKET(1,2,2),
c$$$     &              EpsH(0,1),.TRUE.,SIGMA,psi(1,sigma,1))
c$$$               ampsHstop2(14)=cplModsHstop(2)*S1C(BRAKET(1,2,1),
c$$$     &              EpsH(0,2),.TRUE.,SIGMA,psi(1,sigma,1))
c$$$               ampsHstop2(15)=cplModsHstop(2)*S1C(psi(1,SIGMA,2),
c$$$     &              EpsH(0,1),.TRUE.,SIGMA,BRAKET(1,1,2)) 
c$$$           endif
c... cp-even h^0 vertex with sbottom1-loop (MSSM)
c$$$            if (sHloopSbottom1) then
c$$$               CALL EpsH_SqTriLoop(K(0,1),epsg(0,1),ph,K_H(0,1),
c$$$     &              epsh(0,1),m_sbottom(1))
c$$$               CALL EpsH_SqTriLoop(K(0,2),epsg(0,2),ph,K_H(0,1),
c$$$     &              epsh(0,2),m_sbottom(1))
c$$$               ampsHsbottom1(1)=cplModsHsbottom(1)
c$$$     &         *H2G_SqLoop(Q(0,1),Q(0,2),J21,Jgg,m_sbottom(1))
c$$$               ampsHsbottom1(2)=cplModsHsbottom(1)
c$$$     &         *GGG_vertex(K(0,1),Q(0,3),Q(0,1),epsg(0,1),
c$$$     &                     epsh(0,2),J21)
c$$$               ampsHsbottom1(3)=cplModsHsbottom(1)
c$$$     &         *GGG_vertex(K(0,2),Q(0,1),Q(0,4),epsg(0,2),
c$$$     &                     J21,epsh(0,1))
c$$$               ampsHsbottom1(4)=cplModsHsbottom(1)
c$$$     &         *S1C(psi(1,SIGMA,2),EpsH(0,2),.TRUE.,SIGMA,
c$$$     &              BRAKET(1,1,1))
c$$$               ampsHsbottom1(5)=cplModsHsbottom(1)
c$$$     &         *S1C(BRAKET(1,2,2),EpsH(0,1),.TRUE.,SIGMA,
c$$$     &              psi(1,sigma,1))
c$$$               ampsHsbottom1(6)=cplModsHsbottom(1)
c$$$     &         *S1C(BRAKET(1,2,1),EpsH(0,2),.TRUE.,SIGMA,
c$$$     &              psi(1,sigma,1))
c$$$               ampsHsbottom1(7)=cplModsHsbottom(1)
c$$$     &         *S1C(psi(1,SIGMA,2),EpsH(0,1),.TRUE.,SIGMA,
c$$$     &              BRAKET(1,1,2))
c$$$               ampsHsbottom1(8)=cplModsHsbottom(1)
c$$$     &         *H3G_SqLoop(K(0,1),K(0,2),Q(0,1),epsg(0,1),
c$$$     &              epsg(0,2),J21(0),m_sbottom(1))
c$$$               CALL Eps_SqBlobLoop(K(0,1),epsg(0,1),ph,K_H(0,1),
c$$$     &              epsh(0,1),m_sbottom(1))
c$$$               CALL Eps_SqBlobLoop(K(0,2),epsg(0,2),ph,K_H(0,1),
c$$$     &              epsh(0,2),m_sbottom(1))
c$$$               ampsHsbottom1(9)=cplModsHsbottom(1)
c$$$     &         *H1G_SqLoop(Q(0,1),Q(0,2),J21,Jgg,m_sbottom(1))            
c$$$               ampsHsbottom1(10)=cplModsHsbottom(1)
c$$$     &         *GGG_vertex(K(0,1),Q(0,3),Q(0,1),epsg(0,1),
c$$$     &              epsh(0,2),J21)
c$$$               ampsHsbottom1(11)=cplModsHsbottom(1)
c$$$     &         *GGG_vertex(K(0,2),Q(0,1),Q(0,4),epsg(0,2),
c$$$     &              J21,epsh(0,1))
c$$$               ampsHsbottom1(12)=cplModsHsbottom(1)
c$$$     &         *S1C(psi(1,SIGMA,2),EpsH(0,2),.TRUE.,SIGMA,
c$$$     &              BRAKET(1,1,1)) 
c$$$               ampsHsbottom1(13)=cplModsHsbottom(1)
c$$$     &         *S1C(BRAKET(1,2,2),EpsH(0,1),.TRUE.,SIGMA,
c$$$     &              psi(1,sigma,1))
c$$$               ampsHsbottom1(14)=cplModsHsbottom(1)
c$$$     &         *S1C(BRAKET(1,2,1),EpsH(0,2),.TRUE.,SIGMA,
c$$$     &              psi(1,sigma,1))
c$$$               ampsHsbottom1(15)=cplModsHsbottom(1)
c$$$     &         *S1C(psi(1,SIGMA,2),EpsH(0,1),.TRUE.,SIGMA,
c$$$     &              BRAKET(1,1,2)) 
c$$$            endif
c... cp-even h^0 vertex with sbottom2-loop (MSSM)
c$$$            if (sHloopSbottom2) then
c$$$               CALL EpsH_SqTriLoop(K(0,1),epsg(0,1),ph,K_H(0,1),
c$$$     &              epsh(0,1),m_sbottom(2))
c$$$               CALL EpsH_SqTriLoop(K(0,2),epsg(0,2),ph,K_H(0,1),
c$$$     &              epsh(0,2),m_sbottom(2))
c$$$               ampsHsbottom2(1)=cplModsHsbottom(2)
c$$$     &         *H2G_SqLoop(Q(0,1),Q(0,2),J21,Jgg,m_sbottom(2))
c$$$               ampsHsbottom2(2)=cplModsHsbottom(2)
c$$$     &         *GGG_vertex(K(0,1),Q(0,3),Q(0,1),epsg(0,1),
c$$$     &              epsh(0,2),J21)
c$$$               ampsHsbottom2(3)=cplModsHsbottom(2)
c$$$     &         *GGG_vertex(K(0,2),Q(0,1),Q(0,4),epsg(0,2),
c$$$     &              J21,epsh(0,1))
c$$$               ampsHsbottom2(4)=cplModsHsbottom(2)
c$$$     &         *S1C(psi(1,SIGMA,2),EpsH(0,2),.TRUE.,SIGMA,
c$$$     &              BRAKET(1,1,1))
c$$$               ampsHsbottom2(5)=cplModsHsbottom(2)
c$$$     &         *S1C(BRAKET(1,2,2),EpsH(0,1),.TRUE.,SIGMA,
c$$$     &              psi(1,sigma,1))
c$$$               ampsHsbottom2(6)=cplModsHsbottom(2)
c$$$     &         *S1C(BRAKET(1,2,1),EpsH(0,2),.TRUE.,SIGMA,
c$$$     &              psi(1,sigma,1))
c$$$               ampsHsbottom2(7)=cplModsHsbottom(2)
c$$$     &         *S1C(psi(1,SIGMA,2),EpsH(0,1),.TRUE.,SIGMA,
c$$$     &              BRAKET(1,1,2))
c$$$               ampsHsbottom2(8)=cplModsHsbottom(2)
c$$$     &         *H3G_SqLoop(K(0,1),K(0,2),Q(0,1),epsg(0,1),
c$$$     &              epsg(0,2),J21(0),m_sbottom(2))
c$$$               CALL Eps_SqBlobLoop(K(0,1),epsg(0,1),ph,K_H(0,1),
c$$$     &              epsh(0,1),m_sbottom(2))
c$$$               CALL Eps_SqBlobLoop(K(0,2),epsg(0,2),ph,K_H(0,1),
c$$$     &              epsh(0,2),m_sbottom(2))
c$$$               ampsHsbottom2(9)=cplModsHsbottom(2)
c$$$     &         *H1G_SqLoop(Q(0,1),Q(0,2),J21,Jgg,m_sbottom(2))            
c$$$               ampsHsbottom2(10)=cplModsHsbottom(2)
c$$$     &         *GGG_vertex(K(0,1),Q(0,3),Q(0,1),epsg(0,1),
c$$$     &              epsh(0,2),J21)
c$$$               ampsHsbottom2(11)=cplModsHsbottom(2)
c$$$     &         *GGG_vertex(K(0,2),Q(0,1),Q(0,4),epsg(0,2),
c$$$     &              J21,epsh(0,1))
c$$$               ampsHsbottom2(12)=cplModsHsbottom(2)
c$$$     &         *S1C(psi(1,SIGMA,2),EpsH(0,2),.TRUE.,SIGMA,
c$$$     &              BRAKET(1,1,1)) 
c$$$               ampsHsbottom2(13)=cplModsHsbottom(2)
c$$$     &         *S1C(BRAKET(1,2,2),EpsH(0,1),.TRUE.,SIGMA,
c$$$     &              psi(1,sigma,1))
c$$$               ampsHsbottom2(14)=cplModsHsbottom(2)
c$$$     &         *S1C(BRAKET(1,2,1),EpsH(0,2),.TRUE.,SIGMA,
c$$$     &              psi(1,sigma,1))
c$$$               ampsHsbottom2(15)=cplModsHsbottom(2)
c$$$     &         *S1C(psi(1,SIGMA,2),EpsH(0,1),.TRUE.,SIGMA,
c$$$     &              BRAKET(1,1,2)) 
c$$$            endif
c... cp-even H^0 vertex with stop1-loop (MSSM)
c$$$            if (bHloopStop1) then
c$$$               CALL EpsH_SqTriLoop(K(0,1),epsg(0,1),ph,K_H(0,1),
c$$$     &              epsh(0,1),m_stop(1))
c$$$               CALL EpsH_SqTriLoop(K(0,2),epsg(0,2),ph,K_H(0,1),
c$$$     &              epsh(0,2),m_stop(1))
c$$$               ampbHstop1(1)=cplModbHstop(1)
c$$$     &         *H2G_SqLoop(Q(0,1),Q(0,2),J21,Jgg,m_stop(1))
c$$$               ampbHstop1(2)=cplModbHstop(1)
c$$$     &         *GGG_vertex(K(0,1),Q(0,3),Q(0,1),epsg(0,1),
c$$$     &              epsh(0,2),J21)
c$$$               ampbHstop1(3)=cplModbHstop(1)
c$$$     &         *GGG_vertex(K(0,2),Q(0,1),Q(0,4),epsg(0,2),
c$$$     &              J21,epsh(0,1))
c$$$               ampbHstop1(4)=cplModbHstop(1)
c$$$     &         *S1C(psi(1,SIGMA,2),EpsH(0,2),.TRUE.,SIGMA,
c$$$     &              BRAKET(1,1,1))
c$$$               ampbHstop1(5)=cplModbHstop(1)
c$$$     &         *S1C(BRAKET(1,2,2),EpsH(0,1),.TRUE.,SIGMA,
c$$$     &              psi(1,sigma,1))
c$$$               ampbHstop1(6)=cplModbHstop(1)
c$$$     &         *S1C(BRAKET(1,2,1),EpsH(0,2),.TRUE.,SIGMA,
c$$$     &              psi(1,sigma,1))
c$$$               ampbHstop1(7)=cplModbHstop(1)
c$$$     &         *S1C(psi(1,SIGMA,2),EpsH(0,1),.TRUE.,SIGMA,
c$$$     &              BRAKET(1,1,2))
c$$$               ampbHstop1(8)=cplModbHstop(1)
c$$$     &         *H3G_SqLoop(K(0,1),K(0,2),Q(0,1),epsg(0,1),
c$$$     &              epsg(0,2),J21(0),m_stop(1))
c$$$               CALL Eps_SqBlobLoop(K(0,1),epsg(0,1),ph,K_H(0,1),
c$$$     &              epsh(0,1),m_stop(1))
c$$$               CALL Eps_SqBlobLoop(K(0,2),epsg(0,2),ph,K_H(0,1),
c$$$     &              epsh(0,2),m_stop(1))
c$$$               ampbHstop1(9)=cplModbHstop(1)
c$$$     &         *H1G_SqLoop(Q(0,1),Q(0,2),J21,Jgg,m_stop(1))               
c$$$               ampbHstop1(10)=cplModbHstop(1)
c$$$     &         *GGG_vertex(K(0,1),Q(0,3),Q(0,1),epsg(0,1),
c$$$     &              epsh(0,2),J21)
c$$$               ampbHstop1(11)=cplModbHstop(1)
c$$$     &         *GGG_vertex(K(0,2),Q(0,1),Q(0,4),epsg(0,2),
c$$$     &              J21,epsh(0,1))
c$$$               ampbHstop1(12)=cplModbHstop(1)
c$$$     &         *S1C(psi(1,SIGMA,2),EpsH(0,2),.TRUE.,SIGMA,
c$$$     &              BRAKET(1,1,1)) 
c$$$               ampbHstop1(13)=cplModbHstop(1)
c$$$     &         *S1C(BRAKET(1,2,2),EpsH(0,1),.TRUE.,SIGMA,
c$$$     &              psi(1,sigma,1))
c$$$               ampbHstop1(14)=cplModbHstop(1)
c$$$     &         *S1C(BRAKET(1,2,1),EpsH(0,2),.TRUE.,SIGMA,
c$$$     &              psi(1,sigma,1))
c$$$               ampbHstop1(15)=cplModbHstop(1)
c$$$     &         *S1C(psi(1,SIGMA,2),EpsH(0,1),.TRUE.,SIGMA,
c$$$     &              BRAKET(1,1,2)) 
c$$$            endif
c... cp-even H^0 vertex with stop2-loop (MSSM)
c$$$           if (bHloopStop2) then
c$$$               CALL EpsH_SqTriLoop(K(0,1),epsg(0,1),ph,K_H(0,1),
c$$$     &             epsh(0,1),m_stop(2))
c$$$               CALL EpsH_SqTriLoop(K(0,2),epsg(0,2),ph,K_H(0,1),
c$$$     &              epsh(0,2),m_stop(2))
c$$$               ampbHstop2(1)=cplModbHstop(2)
c$$$     &         *H2G_SqLoop(Q(0,1),Q(0,2),J21,Jgg,m_stop(2))
c$$$               ampbHstop2(2)=cplModbHstop(2)
c$$$     &         *GGG_vertex(K(0,1),Q(0,3),Q(0,1),epsg(0,1),
c$$$     &              epsh(0,2),J21)
c$$$               ampbHstop2(3)=cplModbHstop(2)
c$$$     &         *GGG_vertex(K(0,2),Q(0,1),Q(0,4),epsg(0,2),
c$$$     &              J21,epsh(0,1))
c$$$               ampbHstop2(4)=cplModbHstop(2)
c$$$     &         *S1C(psi(1,SIGMA,2),EpsH(0,2),.TRUE.,SIGMA,
c$$$     &              BRAKET(1,1,1))
c$$$               ampbHstop2(5)=cplModbHstop(2)
c$$$     &         *S1C(BRAKET(1,2,2),EpsH(0,1),.TRUE.,SIGMA,
c$$$     &              psi(1,sigma,1))
c$$$               ampbHstop2(6)=cplModbHstop(2)
c$$$     &         *S1C(BRAKET(1,2,1),EpsH(0,2),.TRUE.,SIGMA,
c$$$     &              psi(1,sigma,1))
c$$$               ampbHstop2(7)=cplModbHstop(2)
c$$$     &         *S1C(psi(1,SIGMA,2),EpsH(0,1),.TRUE.,SIGMA,
c$$$     &              BRAKET(1,1,2))
c$$$               ampbHstop2(8)=cplModbHstop(2)
c$$$     &         *H3G_SqLoop(K(0,1),K(0,2),Q(0,1),epsg(0,1),
c$$$     &              epsg(0,2),J21(0),m_stop(2))
c$$$               CALL Eps_SqBlobLoop(K(0,1),epsg(0,1),ph,K_H(0,1),
c$$$     &              epsh(0,1),m_stop(2))
c$$$               CALL Eps_SqBlobLoop(K(0,2),epsg(0,2),ph,K_H(0,1),
c$$$     &              epsh(0,2),m_stop(2))
c$$$               ampbHstop2(9)=cplModbHstop(2)
c$$$     &         *H1G_SqLoop(Q(0,1),Q(0,2),J21,Jgg,m_stop(2))               
c$$$               ampbHstop2(10)=cplModbHstop(2)
c$$$     &         *GGG_vertex(K(0,1),Q(0,3),Q(0,1),epsg(0,1),
c$$$     &              epsh(0,2),J21)
c$$$               ampbHstop2(11)=cplModbHstop(2)
c$$$     &         *GGG_vertex(K(0,2),Q(0,1),Q(0,4),epsg(0,2),
c$$$     &              J21,epsh(0,1))
c$$$               ampbHstop2(12)=cplModbHstop(2)
c$$$     &         *S1C(psi(1,SIGMA,2),EpsH(0,2),.TRUE.,SIGMA,
c$$$     &              BRAKET(1,1,1)) 
c$$$               ampbHstop2(13)=cplModbHstop(2)
c$$$     &         *S1C(BRAKET(1,2,2),EpsH(0,1),.TRUE.,SIGMA,
c$$$     &              psi(1,sigma,1))
c$$$               ampbHstop2(14)=cplModbHstop(2)
c$$$     &         *S1C(BRAKET(1,2,1),EpsH(0,2),.TRUE.,SIGMA,
c$$$     &              psi(1,sigma,1))
c$$$               ampbHstop2(15)=cplModbHstop(2)
c$$$     &         *S1C(psi(1,SIGMA,2),EpsH(0,1),.TRUE.,SIGMA,
c$$$     &              BRAKET(1,1,2)) 
c$$$           endif
c... cp-even H^0 vertex with sbottom1-loop (MSSM)
c$$$            if (bHloopSbottom1) then
c$$$               CALL EpsH_SqTriLoop(K(0,1),epsg(0,1),ph,K_H(0,1),
c$$$     &              epsh(0,1),m_sbottom(1))
c$$$               CALL EpsH_SqTriLoop(K(0,2),epsg(0,2),ph,K_H(0,1),
c$$$     &              epsh(0,2),m_sbottom(1))
c$$$               ampbHsbottom1(1)=cplModbHsbottom(1)
c$$$     &         *H2G_SqLoop(Q(0,1),Q(0,2),J21,Jgg,m_sbottom(1))
c$$$               ampbHsbottom1(2)=cplModbHsbottom(1)
c$$$     &         *GGG_vertex(K(0,1),Q(0,3),Q(0,1),epsg(0,1),
c$$$     &              epsh(0,2),J21)
c$$$               ampbHsbottom1(3)=cplModbHsbottom(1)
c$$$     &         *GGG_vertex(K(0,2),Q(0,1),Q(0,4),epsg(0,2),
c$$$     &              J21,epsh(0,1))
c$$$               ampbHsbottom1(4)=cplModbHsbottom(1)
c$$$     &         *S1C(psi(1,SIGMA,2),EpsH(0,2),.TRUE.,SIGMA,
c$$$     &              BRAKET(1,1,1))
c$$$               ampbHsbottom1(5)=cplModbHsbottom(1)
c$$$     &         *S1C(BRAKET(1,2,2),EpsH(0,1),.TRUE.,SIGMA,
c$$$     &              psi(1,sigma,1))
c$$$               ampbHsbottom1(6)=cplModbHsbottom(1)
c$$$     &         *S1C(BRAKET(1,2,1),EpsH(0,2),.TRUE.,SIGMA,
c$$$     &              psi(1,sigma,1))
c$$$               ampbHsbottom1(7)=cplModbHsbottom(1)
c$$$     &         *S1C(psi(1,SIGMA,2),EpsH(0,1),.TRUE.,SIGMA,
c$$$     &              BRAKET(1,1,2))
c$$$               ampbHsbottom1(8)=cplModbHsbottom(1)
c$$$     &         *H3G_SqLoop(K(0,1),K(0,2),Q(0,1),epsg(0,1),
c$$$     &              epsg(0,2),J21(0),m_sbottom(1))
c$$$               CALL Eps_SqBlobLoop(K(0,1),epsg(0,1),ph,K_H(0,1),
c$$$     &              epsh(0,1),m_sbottom(1))
c$$$               CALL Eps_SqBlobLoop(K(0,2),epsg(0,2),ph,K_H(0,1),
c$$$     &              epsh(0,2),m_sbottom(1))
c$$$               ampbHsbottom1(9)=cplModbHsbottom(1)
c$$$     &         *H1G_SqLoop(Q(0,1),Q(0,2),J21,Jgg,m_sbottom(1))            
c$$$               ampbHsbottom1(10)=cplModbHsbottom(1)
c$$$     &         *GGG_vertex(K(0,1),Q(0,3),Q(0,1),epsg(0,1),
c$$$     &              epsh(0,2),J21)
c$$$               ampbHsbottom1(11)=cplModbHsbottom(1)
c$$$     &         *GGG_vertex(K(0,2),Q(0,1),Q(0,4),epsg(0,2),
c$$$     &              J21,epsh(0,1))
c$$$               ampbHsbottom1(12)=cplModbHsbottom(1)
c$$$     &         *S1C(psi(1,SIGMA,2),EpsH(0,2),.TRUE.,SIGMA,
c$$$     &              BRAKET(1,1,1)) 
c$$$               ampbHsbottom1(13)=cplModbHsbottom(1)
c$$$     &         *S1C(BRAKET(1,2,2),EpsH(0,1),.TRUE.,SIGMA,
c$$$     &              psi(1,sigma,1))
c$$$               ampbHsbottom1(14)=cplModbHsbottom(1)
c$$$     &         *S1C(BRAKET(1,2,1),EpsH(0,2),.TRUE.,SIGMA,
c$$$     &              psi(1,sigma,1))
c$$$               ampbHsbottom1(15)=cplModbHsbottom(1)
c$$$     &         *S1C(psi(1,SIGMA,2),EpsH(0,1),.TRUE.,SIGMA,
c$$$     &              BRAKET(1,1,2)) 
c$$$            endif
c... cp-even H^0 vertex with sbottom2-loop (MSSM)
c$$$            if (bHloopSbottom2) then
c$$$               CALL EpsH_SqTriLoop(K(0,1),epsg(0,1),ph,K_H(0,1),
c$$$     &              epsh(0,1),m_sbottom(2))
c$$$               CALL EpsH_SqTriLoop(K(0,2),epsg(0,2),ph,K_H(0,1),
c$$$     &              epsh(0,2),m_sbottom(2))
c$$$               ampbHsbottom2(1)=cplModbHsbottom(2)
c$$$     &         *H2G_SqLoop(Q(0,1),Q(0,2),J21,Jgg,m_sbottom(2))
c$$$               ampbHsbottom2(2)=cplModbHsbottom(2)
c$$$     &         *GGG_vertex(K(0,1),Q(0,3),Q(0,1),epsg(0,1),
c$$$     &              epsh(0,2),J21)
c$$$               ampbHsbottom2(3)=cplModbHsbottom(2)
c$$$     &         *GGG_vertex(K(0,2),Q(0,1),Q(0,4),epsg(0,2),
c$$$     &              J21,epsh(0,1))
c$$$               ampbHsbottom2(4)=cplModbHsbottom(2)
c$$$     &         *S1C(psi(1,SIGMA,2),EpsH(0,2),.TRUE.,SIGMA,
c$$$     &              BRAKET(1,1,1))
c$$$               ampbHsbottom2(5)=cplModbHsbottom(2)
c$$$     &         *S1C(BRAKET(1,2,2),EpsH(0,1),.TRUE.,SIGMA,
c$$$     &              psi(1,sigma,1))
c$$$               ampbHsbottom2(6)=cplModbHsbottom(2)
c$$$     &         *S1C(BRAKET(1,2,1),EpsH(0,2),.TRUE.,SIGMA,
c$$$     &              psi(1,sigma,1))
c$$$               ampbHsbottom2(7)=cplModbHsbottom(2)
c$$$     &         *S1C(psi(1,SIGMA,2),EpsH(0,1),.TRUE.,SIGMA,
c$$$     &              BRAKET(1,1,2))
c$$$               ampbHsbottom2(8)=cplModbHsbottom(2)
c$$$     &         *H3G_SqLoop(K(0,1),K(0,2),Q(0,1),epsg(0,1),
c$$$     &              epsg(0,2),J21(0),m_sbottom(2))
c$$$               CALL Eps_SqBlobLoop(K(0,1),epsg(0,1),ph,K_H(0,1),
c$$$     &              epsh(0,1),m_sbottom(2))
c$$$               CALL Eps_SqBlobLoop(K(0,2),epsg(0,2),ph,K_H(0,1),
c$$$     &              epsh(0,2),m_sbottom(2))
c$$$               ampbHsbottom2(9)=cplModbHsbottom(2)
c$$$     &         *H1G_SqLoop(Q(0,1),Q(0,2),J21,Jgg,m_sbottom(2))            
c$$$               ampbHsbottom2(10)=cplModbHsbottom(2)
c$$$     &         *GGG_vertex(K(0,1),Q(0,3),Q(0,1),epsg(0,1),
c$$$     &              epsh(0,2),J21)
c$$$               ampbHsbottom2(11)=cplModbHsbottom(2)
c$$$     &         *GGG_vertex(K(0,2),Q(0,1),Q(0,4),epsg(0,2),
c$$$     &              J21,epsh(0,1))
c$$$               ampbHsbottom2(12)=cplModbHsbottom(2)
c$$$     &         *S1C(psi(1,SIGMA,2),EpsH(0,2),.TRUE.,SIGMA,
c$$$     &              BRAKET(1,1,1)) 
c$$$               ampbHsbottom2(13)=cplModbHsbottom(2)
c$$$     &         *S1C(BRAKET(1,2,2),EpsH(0,1),.TRUE.,SIGMA,
c$$$     &              psi(1,sigma,1))
c$$$               ampbHsbottom2(14)=cplModbHsbottom(2)
c$$$     &         *S1C(BRAKET(1,2,1),EpsH(0,2),.TRUE.,SIGMA,
c$$$     &              psi(1,sigma,1))
c$$$               ampbHsbottom2(15)=cplModbHsbottom(2)
c$$$     &         *S1C(psi(1,SIGMA,2),EpsH(0,1),.TRUE.,SIGMA,
c$$$     &              BRAKET(1,1,2)) 
c$$$            endif

c...  amp(1), amp(9):         amp(8):
c         p1 -->--*-->-- p2        p1 -->--*-->-- p2 
c                 |                        |
c                 O-- H / A0       k1 --<--O-->-- k2      
c                 |                        |
c         k1 --<--*-->-- k2               H/A0   
c   
c...  amp(2), amp(10):         amp(3), amp(11):
c        p1 -->--*--->-- p2         p1 -->---*-->-- p2 
c                |                           |
c                |  H/A0                H/A0 |
c                |   |                   |   |
c        k1 --<--*->-O-->- k2     k1 -<--O-<-*-->-- k2        
c
c...  amp(4), amp(12):         amp(5), amp(13):
c      p1 -->-*->-*-->-- p2         p1 -->-*->-*-->-- p2 
c            /    |                        |    \
c           /     O-- H/A0          H/A0 --O     \
c          /      |                        |      \
c        k1       k2                       k1      k2        
c
c...  amp(7), amp(14):         amp(6), amp(15):
c      p1 -->-*->-*-->-- p2         p1 -->-*->-*-->-- p2 
c            /    |                        |    \
c           /     O-- H/A0          H/A0 --O     \
c          /      |                        |      \
c        k2       k1                       k2      k1   
            
c... sum up all subamplitudes            
            m12=ampHtop(1)+ampHbottom(1)+ampHeff(1)+ampAeff(1)
     &         +ampAtop(1)+ampAbottom(1)
     &         +ampsHtop(1)+ampsHbottom(1)+ampbHtop(1)+ampbHbottom(1)
     &         +ampsHstop1(1)+ampsHstop2(1)+ampbHstop1(1)+ampbHstop2(1)
     &         +ampsHsbottom1(1)+ampsHsbottom2(1)+ampbHsbottom1(1)
     &         +ampbHsbottom2(1)
     &         +ampHtop(2)+ampHbottom(2)+ampHeff(2)+ampAeff(2)
     &         +ampAtop(2)+ampAbottom(2)
     &         +ampsHtop(2)+ampsHbottom(2)+ampbHtop(2)+ampbHbottom(2)
     &         +ampsHstop1(2)+ampsHstop2(2)+ampbHstop1(2)+ampbHstop2(2)
     &         +ampsHsbottom1(2)+ampsHsbottom2(2)+ampbHsbottom1(2)
     &         +ampbHsbottom2(2)
     &         +ampHtop(3)+ampHbottom(3)+ampHeff(3)+ampAeff(3)
     &         +ampAtop(3)+ampAbottom(3)
     &         +ampsHtop(3)+ampsHbottom(3)+ampbHtop(3)+ampbHbottom(3)
     &         +ampsHstop1(3)+ampsHstop2(3)+ampbHstop1(3)+ampbHstop2(3)
     &         +ampsHsbottom1(3)+ampsHsbottom2(3)+ampbHsbottom1(3)
     &         +ampbHsbottom2(3)
     &         +ampHtop(6)+ampHbottom(6)+ampHeff(6)+ampAeff(6)
     &         +ampAtop(6)+ampAbottom(6)
     &         +ampsHtop(6)+ampsHbottom(6)+ampbHtop(6)+ampbHbottom(6)
     &         +ampsHstop1(6)+ampsHstop2(6)+ampbHstop1(6)+ampbHstop2(6)
     &         +ampsHsbottom1(6)+ampsHsbottom2(6)+ampbHsbottom1(6)
     &         +ampbHsbottom2(6)
     &         +ampHtop(7)+ampHbottom(7)+ampHeff(7)+ampAeff(7)
     &         +ampAtop(7)+ampAbottom(7)
     &         +ampsHtop(7)+ampsHbottom(7)+ampbHtop(7)+ampbHbottom(7)
     &         +ampsHstop1(7)+ampsHstop2(7)+ampbHstop1(7)+ampbHstop2(7)
     &         +ampsHsbottom1(7)+ampsHsbottom2(7)+ampbHsbottom1(7)
     &         +ampbHsbottom2(7)
     &         +ampHtop(8)+ampHbottom(8)+ampHeff(8)+ampAeff(8)
     &         +ampAtop(8)+ampAbottom(8)
     &         +ampsHtop(8)+ampsHbottom(8)+ampbHtop(8)+ampbHbottom(8)
     &         +ampsHstop1(8)+ampsHstop2(8)+ampbHstop1(8)+ampbHstop2(8)
     &         +ampsHsbottom1(8)+ampsHsbottom2(8)+ampbHsbottom1(8)
     &         +ampbHsbottom2(8)
     &         +ampsHstop1(9)+ampsHstop2(9)+ampsHsbottom1(9)
     &         +ampsHsbottom2(9)
     &         +ampbHstop1(9)+ampbHstop2(9)+ampbHsbottom1(9)
     &         +ampbHsbottom2(9)
     &         +ampsHstop1(10)+ampsHstop2(10)+ampsHsbottom1(10)
     &         +ampsHsbottom2(10)
     &         +ampbHstop1(10)+ampbHstop2(10)+ampbHsbottom1(10)
     &         +ampbHsbottom2(10)
     &         +ampsHstop1(11)+ampsHstop2(11)+ampsHsbottom1(11)
     &         +ampsHsbottom2(11)
     &         +ampbHstop1(11)+ampbHstop2(11)+ampbHsbottom1(11)
     &         +ampbHsbottom2(11)
     &         +ampsHstop1(14)+ampsHstop2(14)+ampsHsbottom1(14)
     &         +ampsHsbottom2(14)
     &         +ampbHstop1(14)+ampbHstop2(14)+ampbHsbottom1(14)
     &         +ampbHsbottom2(14)
     &         +ampsHstop1(15)+ampsHstop2(15)+ampsHsbottom1(15)
     &         +ampsHsbottom2(15)
     &         +ampbHstop1(15)+ampbHstop2(15)+ampbHsbottom1(15)
     &         +ampbHsbottom2(15)

            m21=-ampHtop(1)-ampHbottom(1)-ampHeff(1)-ampAeff(1)
     &          -ampAtop(1)-ampAbottom(1)
     &          -ampsHtop(1)-ampsHbottom(1)-ampbHtop(1)-ampbHbottom(1)
     &          -ampsHstop1(1)-ampsHstop2(1)-ampbHstop1(1)-ampbHstop2(1)
     &          -ampsHsbottom1(1)-ampsHsbottom2(1)-ampbHsbottom1(1)
     &          -ampbHsbottom2(1)
     &          -ampHtop(2)-ampHbottom(2)-ampHeff(2)-ampAeff(2)
     &          -ampAtop(2)-ampAbottom(2)
     &          -ampsHtop(2)-ampsHbottom(2)-ampbHtop(2)-ampbHbottom(2)
     &          -ampsHstop1(2)-ampsHstop2(2)-ampbHstop1(2)-ampbHstop2(2)
     &          -ampsHsbottom1(2)-ampsHsbottom2(2)-ampbHsbottom1(2)
     &          -ampbHsbottom2(2)
     &          -ampHtop(3)-ampHbottom(3)-ampHeff(3)-ampAeff(3)
     &          -ampAtop(3)-ampAbottom(3)
     &          -ampsHtop(3)-ampsHbottom(3)-ampbHtop(3)-ampbHbottom(3)
     &          -ampsHstop1(3)-ampsHstop2(3)-ampbHstop1(3)-ampbHstop2(3)
     &          -ampsHsbottom1(3)-ampsHsbottom2(3)-ampbHsbottom1(3)
     &          -ampbHsbottom2(3)
     &          +ampHtop(4)+ampHbottom(4)+ampHeff(4)+ampAeff(4)
     &          +ampAtop(4)+ampAbottom(4)
     &          +ampsHtop(4)+ampsHbottom(4)+ampbHtop(4)+ampbHbottom(4)
     &          +ampsHstop1(4)+ampsHstop2(4)+ampbHstop1(4)+ampbHstop2(4)
     &          +ampsHsbottom1(4)+ampsHsbottom2(4)+ampbHsbottom1(4)
     &          +ampbHsbottom2(4)
     &          +ampHtop(5)+ampHbottom(5)+ampHeff(5)+ampAeff(5)
     &          +ampAtop(5)+ampAbottom(5)
     &          +ampsHtop(5)+ampsHbottom(5)+ampbHtop(5)+ampbHbottom(5)
     &          +ampsHstop1(5)+ampsHstop2(5)+ampbHstop1(5)+ampbHstop2(5)
     &          +ampsHsbottom1(5)+ampsHsbottom2(5)+ampbHsbottom1(5)
     &          +ampbHsbottom2(5)
     &          -ampHtop(8)-ampHbottom(8)-ampHeff(8)-ampAeff(8)
     &          -ampAtop(8)-ampAbottom(8)
     &          -ampsHtop(8)-ampsHbottom(8)-ampbHtop(8)-ampbHbottom(8)
     &          -ampsHstop1(8)-ampsHstop2(8)-ampbHstop1(8)-ampbHstop2(8)
     &          -ampsHsbottom1(8)-ampsHsbottom2(8)-ampbHsbottom1(8)
     &          -ampbHsbottom2(8)
     &          -ampsHstop1(9)-ampsHstop2(9)-ampsHsbottom1(9)
     &          -ampsHsbottom2(9)
     &          -ampbHstop1(9)-ampbHstop2(9)-ampbHsbottom1(9)
     &          -ampbHsbottom2(9)
     &          -ampsHstop1(10)-ampsHstop2(10)-ampsHsbottom1(10)
     &          -ampsHsbottom2(10)
     &          -ampbHstop1(10)-ampbHstop2(10)-ampbHsbottom1(10)
     &          -ampbHsbottom2(10)
     &          -ampsHstop1(11)-ampsHstop2(11)-ampsHsbottom1(11)
     &          -ampsHsbottom2(11)
     &          -ampbHstop1(11)-ampbHstop2(11)-ampbHsbottom1(11)
     &          -ampbHsbottom2(11)
     &          +ampsHstop1(12)+ampsHstop2(12)+ampsHsbottom1(12)
     &          +ampsHsbottom2(12)
     &          +ampbHstop1(12)+ampbHstop2(12)+ampbHsbottom1(12)
     &          +ampbHsbottom2(12)
     &          +ampsHstop1(13)+ampsHstop2(13)+ampsHsbottom1(13)
     &          +ampsHsbottom2(13)
     &          +ampbHstop1(13)+ampbHstop2(13)+ampbHsbottom1(13)
     &          +ampbHsbottom2(13)

c... square the amplitudes     
            if (i1pol.ne.3) then
               res12=res12 + dreal(m12)**2+dimag(m12)**2 
               res21=res21 + dreal(m21)**2+dimag(m21)**2 
               resfi=resfi + 2d0*( dreal(m12)*dreal(m21)
     &              +dimag(m12)*dimag(m21) )
            endif
c... gauge test
            if (ldebug) then
               if (npolg1.eq.3) then
                  if (i1pol.ge.3) then
                     if ((abs(m12)/dsqrt(resfi).gt.1e-2).or.(abs(m21)
     &                    /dsqrt(resfi).gt.1e-2)) then
                        print *,'Warning in HiggsQG:  bad gauge check!'
                        print 13,' m12 ',m12
                        print 13,' m21 ',m21
                     endif
 13                  format(a,2g15.5,a,g12.2)
                  endif
               endif
            endif
c-----------------------------------------------
c... end: sum over gluon polarizations
c-----------------------------------------------            
         enddo
      enddo
      
c... multiply with color factors
      rh(1)=res12*col11_22 * fac 
      rh(2)=res21*col11_22 * fac 
      rh(0)=rh(1)+rh(2)+resfi*col12_21*fac 
   
      if (ldebug) then
         print *,rh(0),rh(1),rh(2)
      endif

      END

***********************************************************************
      SUBROUTINE HiggsGG(Qphyin,gpol,gsign,rh,amp,tri,box,pen)
*----------------------------------------------------------------------
c... modified by: Michael Kubocz     
c... date       : 17.10.2007
c... modified   : 30.01.2011
*----------------------------------------------------------------------
c... This routine calculates the color summed amplitude square of the
c... process ' g g -> g g h ' . 

c... Input variables
c...  Qphy(0:3,4) physical momenta of the 4 external gluons 
c...  gpol(4 )    polarization (1 or 2) of the gluons
c...  sign(4)     sign of the gluons : +1 for outgoing, -1 for ingoing

c... Output variables
c...  rh(0 )   the color summed amplitude square 
c...  rh(1-3)  the relative portions for the 3 different colors
c...  amp(3)   the summed amplitude of all graphs for each color factor
c...  tri(3)   same, but summing only the graphs with a ggh vertex
c...  box(3)   same for the graphs with a gggh vertex
c...  pen(3)   same for the graphs with a ggggh vertex
***********************************************************************

      IMPLICIT NONE
 
#include "ggf_process.inc"
#include "ggf_dotproducts.inc"

c... input
      real*8 Qphyin(0:3,4), Qphy(0:4,4)
      integer gpol(4),gsign(4)
c... output
      real*8 rh(0:3)
      complex*16 amp(3)
      complex*16 tri(3),box(3),pen(3)
      complex*16 triAmp(3),boxAmp(3),penAmp(3)

c... global factors (f2 is just for x-check with HGGGG)      
      real*8 C1,C2,vev,fac
      PARAMETER (C1=7.6666666666666666d0, C2=-1.3333333333333333d0)
      
c... common block for coupling constants etc.
      real*8 ALFAS,XMT,ALFA,XMZ,XMW,SIN2W,XMH, gf, alfas_lo
      real*8 clr, xm2, xmg, b, v, a
      COMMON /BKOPIN/ ALFAS,XMT,ALFA,XMZ,XMW,SIN2W,XMH, gf, alfas_lo
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)

c... internal momenta, polarization, current and dot-product variables
      complex*16 epvec(0:3,4)
      real*8 epvecR(0:3,4) ! same as epvec but real
      real*8  K(0:4,5),K_ij(0:4,4,4),K_H(0:3,4),Kpen(0:4,4)
      complex*16 J_ij(0:3,4,4)
      complex*16 eps_effH(0:3,4),eps_effA(0:3,4),
     &     eps_tH(0:3,4),eps_bH(0:3,4),eps_tA(0:3,4),eps_bA(0:3,4),
     &     eps_st1bH(0:3,4),eps_st1tH(0:3,4),
     &     eps_st2bH(0:3,4),eps_st2tH(0:3,4),
     &     eps_sb1bH(0:3,4),eps_sb1tH(0:3,4),
     &     eps_sb2bH(0:3,4),eps_sb2tH(0:3,4)
      COMMON /AMPINFO/epvec,J_ij,
     &       eps_effH,eps_effA,eps_tH,eps_bH,eps_tA,eps_bA,
     &       eps_st1bH,eps_st1tH,eps_st2bH,eps_st2tH,
     &       eps_sb1bH,eps_sb1tH,eps_sb2bH,eps_sb2tH,K,K_ij,K_H
*      COMMON /AMPINFO/epvec,K,J_ij,K_ij,K_H,
*     &       eps_effH,eps_effA,eps_tH,eps_bH,eps_tA,eps_bA,
*     &       eps_st1bH,eps_st1tH,eps_st2bH,eps_st2tH,
*     &       eps_sb1bH,eps_sb1tH,eps_sb2bH,eps_sb2tH

c... auxiliary variables
      integer i, mu,j
      real*8 ampsqr,amp1sqr,amp2sqr,amp3sqr
      real*8 epx(0:3,4)


c... external dot-product functions and the m_i subamlitudes
      real*8 dotrr
      COMPLEX*16 MMa,MMb,MMc,MMd,MMe
      complex*16 dotrc,dotcc
      EXTERNAL dotrr,dotrc,dotcc

c... debugging variables
      LOGICAL ldebug
      INTEGER gluonid
      PARAMETER (ldebug=.false.,gluonid=3)


      real*8 pb(0:3),Qphyb(0:4,4)
      logical lboost
      parameter (lboost=.false.)

      real*8 eps, qi
      common /bh4g_mom_eps/ qi(0:3,4),eps(0:3,4)


c... global factors
      vev = 2d0*xmw/(v(1,3)*dsqrt(8d0))
c... factor 2 --> c_i color factors
      fac = 2d0*(4d0*alfas**2/vev)

c... gluon vectors (K)
      qphy(0:3,:) = qphyin(0:3,:)
      qphy(4,:) = 0
      qphyb(4,:) = 0

      if (lboost) then
         pb(0) =  Qphy(0,1)+Qphy(0,2)
         pb(1) =-(Qphy(1,1)+Qphy(1,2))
         pb(2) =-(Qphy(2,1)+Qphy(2,2))
         pb(3) =-(Qphy(3,1)+Qphy(3,2))
         do i = 1,4
            call boostn(Qphy(0,i), pb, Qphyb(0,i))
            do mu=0,3
               K(mu,i)=-gsign(i)*Qphyb(mu,i)
            enddo
            call POLVEC(Qphyb(0,i),gpol(i),epx(0,i))
            do mu=0,3
               epvec(mu,i)  =epx(mu,i)
               epvecR(mu,i) =epx(mu,i) ! real Pol.vector 
            enddo
         enddo
      else
         do i=1,4
            do mu=0,3
               K(mu,i)   = -gsign(i)*Qphy(mu,i)
               Kpen(mu,i)= -gsign(i)*Qphy(mu,i)
            enddo
         enddo

c... fill polarization vectors 
c... though in the cartesian base all the epvec's are real,
c... i keep them complex to allow x-checks with the helicity base 
         do i = 1,4
c... cartesian base 
            call POLVEC(Qphy(0,i),gpol(i),epx(0,i))
            do mu=0,3
               epvec(mu,i)  =epx(mu,i)
               epvecR(mu,i) =epx(mu,i) ! real Pol.vector 
            enddo
c... alternative : the helicity base
c     call HELVEC(qphy(0,i),gsign(i),(gpol(i)*2-3),EPvec(0,i))
         enddo
      endif
      do mu = 0,3
         do i = 1,4
            qi(mu,i) = K(mu,i)
            eps(mu,i) = epvecR(mu,i)
         enddo
      enddo

c... set one eps_i=k_i when checking gauge invarianz
c      if (lgauge) then
c         do mu=0,3
c            epvec(mu,gluonid)= K(mu,gluonid)/K(0,gluonid)
c         enddo
c      endif

c... The higgs momentum = K(mu,5)
      do mu=0,3
         K(mu,5)=-(K(mu,1)+K(mu,2)+K(mu,3)+K(mu,4))
      enddo

c... Calculate the 2-Gluon-Currents J_ij (e.g. J_34)
      do i=1,4
         do j=1,4
            if (i.ne.j) then
               CALL GGG_current(K(0,i),epvec(0,i),K(0,j),epvec(0,j),
     &                          K_ij(0,i,j),J_ij(0,i,j))
            endif
         enddo
      enddo

      do i=1,4
         K(4,i)   =0D0
         Kpen(4,i)=0D0
      enddo
 
      do i=1,4
         do j=1,4
            if (i.ne.j) then
               K_ij(4,i,j)=dotrr(K_ij(0,i,j),K_ij(0,i,j))
            endif
         enddo
      enddo

c... calculates dot products and E_0, E_1i for different loop masses
      if (.not.(Heff.or.Aeff)) then
         CALL auxiliary_calculation(epvec,Kpen)
      endif

c... initialize amplitudes with zero values
      do i=1,3
         triAmp(i)=DCMPLX(0d0,0d0)
         boxAmp(i)=DCMPLX(0d0,0d0)
         penAmp(i)=DCMPLX(0d0,0d0)
      enddo

c... prepare effective polarizations vectors for different loop masses
c... cp-even Higgs vertex in effective limes
      if (Heff) then
         do i=1,4
            CALL EpsH_effective(K(0,i),.false.,epvec(0,i),K(0,5),
     &           K_H(0,i),eps_effH(0,i))
         enddo
      endif
c... cp-odd Higgs vertex in effective limes
      if (Aeff) then
         do i=1,4
            CALL EpsA_effective(K(0,i),epvec(0,i),K(0,5),K_H(0,i),
     &                          eps_effA(0,i))
         enddo
      endif
c... cp-even Higgs vertex with top-loop
      if (HloopTop.or.sHloopTop.or.bHloopTop) then
         do i=1,4
            CALL EpsH_Loop(K(0,i),.false.,epvec(0,i),K(0,5),K_H(0,i),
     &                     eps_tH(0,i),m_top)
         enddo
      endif
c... cp-even Higgs vertex with bottom-loop
      if (HloopBottom.or.sHloopBottom.or.bHloopBottom) then
         do i=1,4
            CALL EpsH_Loop(K(0,i),.false.,epvec(0,i),K(0,5),K_H(0,i),
     &                     eps_bH(0,i),m_bottom)  
         enddo   
      endif
c... cp-odd Higgs vertex with top-loop
      if (AloopTop) then
         do i=1,4
            CALL EpsA_Loop(K(0,i),epvec(0,i),K(0,5),K_H(0,i),
     &                     eps_tA(0,i),m_top)
         enddo  
      endif
c... cp-odd Higgs vertex with bottom-loop
      if (AloopBottom) then
         do i=1,4
            CALL EpsA_Loop(K(0,i),epvec(0,i),K(0,5),K_H(0,i),
     &                     eps_bA(0,i),m_bottom)
         enddo   
      endif
c... cp-even Higgs vertex with stop1-loop
c$$$      if (sHloopStop1.or.bHloopStop1) then
c$$$         do i=1,4
c$$$            CALL Eps_SqBlobLoop(K(0,i),epvec(0,i),K(0,5),K_H(0,i),
c$$$     &           eps_st1bH(0,i),m_stop(1))
c$$$            CALL EpsH_SqTriLoop(K(0,i),epvec(0,i),K(0,5),K_H(0,i),
c$$$     &           eps_st1tH(0,i),m_stop(1))
c$$$         enddo
c$$$      endif
c... cp-even Higgs vertex with stop2-loop
c$$$      if (sHloopStop2.or.bHloopStop2) then
c$$$         do i=1,4
c$$$            CALL Eps_SqBlobLoop(K(0,i),epvec(0,i),K(0,5),K_H(0,i),
c$$$     &           eps_st2bH(0,i),m_stop(2))
c$$$            CALL EpsH_SqTriLoop(K(0,i),epvec(0,i),K(0,5),K_H(0,i),
c$$$     &           eps_st2tH(0,i),m_stop(2))
c$$$         enddo
c$$$      endif
c... cp-even Higgs vertex with sbottom1-loop
c$$$      if (sHloopSbottom1.or.bHloopSbottom1) then
c$$$         do i=1,4
c$$$            CALL Eps_SqBlobLoop(K(0,i),epvec(0,i),K(0,5),K_H(0,i),
c$$$     &           eps_sb1bH(0,i),m_sbottom(1))
c$$$            CALL EpsH_SqTriLoop(K(0,i),epvec(0,i),K(0,5),K_H(0,i),
c$$$     &           eps_sb1tH(0,i),m_sbottom(1))
c$$$         enddo
c$$$      endif
c... cp-even Higgs vertex with sbottom2-loop
c$$$      if (sHloopSbottom2.or.bHloopSbottom2) then
c$$$         do i=1,4
c$$$            CALL Eps_SqBlobLoop(K(0,i),epvec(0,i),K(0,5),K_H(0,i),
c$$$     &           eps_sb2bH(0,i),m_sbottom(2))
c$$$            CALL EpsH_SqTriLoop(K(0,i),epvec(0,i),K(0,5),K_H(0,i),
c$$$     &           eps_sb2tH(0,i),m_sbottom(2))
c$$$         enddo
c$$$      endif 

c... sum up the amplitudes for 3 different color structures
      triAmp(1)=-MMa(1,2,4,3)-MMa(4,3,1,2)-MMa(2,1,3,4)-MMa(3,4,2,1)
     &          +MMa(1,3,2,4)+MMa(2,4,1,3)+MMa(3,1,4,2)+MMa(4,2,3,1)
     &          -MMb(1,2,4,3)+MMb(1,3,2,4)
     &          -MMd(1,3,2,4)+MMd(1,2,4,3)-MMd(3,1,4,2)-MMd(3,1,2,4)
     &          -MMd(4,2,3,1)+MMd(4,3,1,2)-MMd(2,4,1,3)-MMd(2,4,3,1)
      triAmp(2)=-MMa(1,2,3,4)-MMa(3,4,1,2)-MMa(2,1,4,3)-MMa(4,3,2,1)
     &          -MMa(1,3,2,4)-MMa(2,4,1,3)-MMa(3,1,4,2)-MMa(4,2,3,1)
     &          -MMb(1,2,3,4)-MMb(1,3,2,4)
     &          +MMd(1,2,3,4)+MMd(1,3,2,4)+MMd(3,1,4,2)+MMd(3,4,1,2)
     &          +MMd(4,3,2,1)+MMd(4,2,3,1)+MMd(2,4,1,3)+MMd(2,1,4,3)
      triAmp(3)=+MMa(1,2,3,4)+MMa(3,4,1,2)+MMa(2,1,4,3)+MMa(4,3,2,1)
     &          +MMa(1,2,4,3)+MMa(4,3,1,2)+MMa(2,1,3,4)+MMa(3,4,2,1)
     &          +MMb(1,2,3,4)+MMb(1,2,4,3)
     &          -MMd(1,2,3,4)-MMd(1,2,4,3)-MMd(3,4,1,2)+MMd(3,1,2,4)
     &          -MMd(4,3,2,1)-MMd(4,3,1,2)-MMd(2,1,4,3)+MMd(2,4,3,1)

      boxAmp(1)=(-MMc(1,2,4,3)+MMc(1,3,2,4)-MMc(2,1,3,4)+MMc(3,1,4,2))
      boxAmp(2)=(-MMc(1,2,3,4)-MMc(1,3,2,4)-MMc(2,1,4,3)-MMc(3,1,4,2))
      boxAmp(3)=(+MMc(1,2,3,4)+MMc(1,2,4,3)+MMc(2,1,4,3)+MMc(2,1,3,4))

      CALL calculate_E_functions()
      penAmp(1)=-(MMe(2,3,4,1)+MMe(4,3,2,1)+MMe(2,1,4,3)+MMe(3,2,1,4))
     &            /2d0
      penAmp(2)=-(MMe(2,1,3,4)+MMe(3,4,2,1)+MMe(4,2,1,3)+MMe(2,4,3,1))
     &            /2d0
      penAmp(3)=-(MMe(2,3,1,4)+MMe(3,1,4,2)+MMe(4,2,3,1)+MMe(3,2,4,1))
     &            /2d0

      amp(1)=(triAmp(1)+boxAmp(1)+penAmp(1))*fac
      amp(2)=(triAmp(2)+boxAmp(2)+penAmp(2))*fac
      amp(3)=(triAmp(3)+boxAmp(3)+penAmp(3))*fac

      amp1sqr = (dreal(amp(1))**2 + dimag(amp(1))**2)
      amp2sqr = (dreal(amp(2))**2 + dimag(amp(2))**2)
      amp3sqr = (dreal(amp(3))**2 + dimag(amp(3))**2)

c... derive colour summed amplitude square      
      ampsqr = C1*(amp1sqr+amp2sqr+amp3sqr)+
     &         C2*2*(
     &         dreal(amp(1))*dreal(amp(2)) +
     &         dimag(amp(1))*dimag(amp(2)) +
     &         dreal(amp(1))*dreal(amp(3)) +
     &         dimag(amp(1))*dimag(amp(3)) +
     &         dreal(amp(2))*dreal(amp(3)) +
     &         dimag(amp(2))*dimag(amp(3))) 

      rh(0) = ampsqr
      rh(1) = amp1sqr
      rh(2) = amp2sqr
      rh(3) = amp3sqr

      END

***********************************************************************
      FUNCTION MMa(i1,i2,i3,i4)
*----------------------------------------------------------------------
c... author  : Gunnar Klaemke
c... modified: Michael Kubocz
c... date    : 29.08.07 (MK)
*----------------------------------------------------------------------

c              H/A0 
c               |
c     p_i1 --<--O--*-->-- p_i3        
c                  |                    
c                  |             
c     p_i2 --<-----*-->-- p_i4
c
***********************************************************************

      IMPLICIT NONE

#include "ggf_process.inc"

      complex*16 MMa
      complex*16 aux(26)

      complex*16 epvec(0:3,4),J_ij(0:3,4,4)
      real*8  K(0:4,5),K_ij(0:4,4,4),K_H(0:3,4) 
      complex*16 eps_effH(0:3,4),eps_effA(0:3,4),
     &     eps_tH(0:3,4),eps_bH(0:3,4),eps_tA(0:3,4),eps_bA(0:3,4),
     &     eps_st1bH(0:3,4),eps_st1tH(0:3,4),
     &     eps_st2bH(0:3,4),eps_st2tH(0:3,4),
     &     eps_sb1bH(0:3,4),eps_sb1tH(0:3,4),
     &     eps_sb2bH(0:3,4),eps_sb2tH(0:3,4)
*      COMMON /AMPINFO/epvec,K,J_ij,K_ij,K_H,
*     &       eps_effH,eps_effA,eps_tH,eps_bH,eps_tA,eps_bA,
*     &       eps_st1bH,eps_st1tH,eps_st2bH,eps_st2tH,
*     &       eps_sb1bH,eps_sb1tH,eps_sb2bH,eps_sb2tH
      COMMON /AMPINFO/epvec,J_ij,
     &       eps_effH,eps_effA,eps_tH,eps_bH,eps_tA,eps_bA,
     &       eps_st1bH,eps_st1tH,eps_st2bH,eps_st2tH,
     &       eps_sb1bH,eps_sb1tH,eps_sb2bH,eps_sb2tH,K,K_ij,K_H

      integer i1,i2,i3,i4,j

      complex*16 GGG_vertex

      EXTERNAL GGG_vertex
            
      MMa=DCMPLX(0d0,0d0)
      do j=1,26
         aux(j)=DCMPLX(0d0,0d0)
      enddo
c... cp-even Higgs vertex for effective limes
      if (Heff) then
         aux(1)=cplModHeff
     &   *GGG_vertex(K_ij(0,i2,i4),K_H(0,i1),K(0,i3),J_ij(0,i2,i4),
     &               eps_effH(0,i1),epvec(0,i3))
c... cp-odd Higgs vertex for effective limes
      endif
      if (Aeff) then
         aux(2)=cplModAeff
     &   *GGG_vertex(K_ij(0,i2,i4),K_H(0,i1),K(0,i3),J_ij(0,i2,i4),
     &               eps_effA(0,i1),epvec(0,i3))
      endif
c... cp-even Higgs vertex with top-loop (SM)
      if (HloopTop) then
         aux(3)=cplModHtop
     &   *GGG_vertex(K_ij(0,i2,i4),K_H(0,i1),K(0,i3),J_ij(0,i2,i4),
     &               eps_tH(0,i1),epvec(0,i3))
      endif
c... cp-even Higgs vertex with bottom-loop (SM)
      if (HloopBottom) then
         aux(4)=cplModHbottom
     &   *GGG_vertex(K_ij(0,i2,i4),K_H(0,i1),K(0,i3),J_ij(0,i2,i4),
     &               eps_bH(0,i1),epvec(0,i3))
      endif
c... cp-odd Higgs vertex with top-loop (MSSM)
      if (AloopTop) then
         aux(5)=cplModAtop
     &   *GGG_vertex(K_ij(0,i2,i4),K_H(0,i1),K(0,i3),J_ij(0,i2,i4),
     &               eps_tA(0,i1),epvec(0,i3))
      endif
c... cp-odd Higgs vertex with bottom-loop (MSSM)
      if (AloopBottom) then
         aux(6)=cplModAbottom
     &   *GGG_vertex(K_ij(0,i2,i4),K_H(0,i1),K(0,i3),J_ij(0,i2,i4),
     &               eps_bA(0,i1),epvec(0,i3))
      endif
c... cp-even Higgs h vertex with top-loop (MSSM)
      if (sHloopTop) then
         aux(7)=cplModsHtop
     &   *GGG_vertex(K_ij(0,i2,i4),K_H(0,i1),K(0,i3),J_ij(0,i2,i4),
     &               eps_tH(0,i1),epvec(0,i3))
      endif
c... cp-even Higgs h vertex with bottom-loop (MSSM)
      if (sHloopBottom) then
         aux(8)=cplModsHbottom
     &   *GGG_vertex(K_ij(0,i2,i4),K_H(0,i1),K(0,i3),J_ij(0,i2,i4),
     &               eps_bH(0,i1),epvec(0,i3))
      endif
c... cp-even Higgs H vertex with top-loop (MSSM)
      if (bHloopTop) then
         aux(9)=cplModbHtop
     &   *GGG_vertex(K_ij(0,i2,i4),K_H(0,i1),K(0,i3),J_ij(0,i2,i4),
     &               eps_tH(0,i1),epvec(0,i3))
      endif
c... cp-even Higgs H vertex with bottom-loop (MSSM)
      if (bHloopBottom) then
         aux(10)=cplModbHbottom
     &   *GGG_vertex(K_ij(0,i2,i4),K_H(0,i1),K(0,i3),J_ij(0,i2,i4),
     &               eps_bH(0,i1),epvec(0,i3))
      endif
c... cp-even h^0 vertex with stop1-loop (MSSM)
c$$$      if (sHloopStop1) then
c$$$         aux(11)=cplModsHstop(1)
c$$$     &   *GGG_vertex(K_ij(0,i2,i4),K_H(0,i1),K(0,i3),J_ij(0,i2,i4),
c$$$     &               eps_st1bH(0,i1),epvec(0,i3))
c$$$         aux(12)=cplModsHstop(1)
c$$$     &   *GGG_vertex(K_ij(0,i2,i4),K_H(0,i1),K(0,i3),J_ij(0,i2,i4),
c$$$     &               eps_st1tH(0,i1),epvec(0,i3))
c$$$      endif
c... cp-even h^0 vertex with stop2-loop (MSSM)
c$$$      if (sHloopStop2) then
c$$$         aux(13)=cplModsHstop(2)
c$$$     &   *GGG_vertex(K_ij(0,i2,i4),K_H(0,i1),K(0,i3),J_ij(0,i2,i4),
c$$$     &               eps_st2bH(0,i1),epvec(0,i3))
c$$$         aux(14)=cplModsHstop(2)
c$$$     &   *GGG_vertex(K_ij(0,i2,i4),K_H(0,i1),K(0,i3),J_ij(0,i2,i4),
c$$$     &               eps_st2tH(0,i1),epvec(0,i3))
c$$$      endif
c... cp-even H^0 vertex with stop1-loop (MSSM)
c$$$      if (bHloopStop1) then
c$$$         aux(15)=cplModbHstop(1)
c$$$     &   *GGG_vertex(K_ij(0,i2,i4),K_H(0,i1),K(0,i3),J_ij(0,i2,i4),
c$$$     &               eps_st1bH(0,i1),epvec(0,i3))
c$$$         aux(16)=cplModbHstop(1)
c$$$     &   *GGG_vertex(K_ij(0,i2,i4),K_H(0,i1),K(0,i3),J_ij(0,i2,i4),
c$$$     &               eps_st1tH(0,i1),epvec(0,i3))
c$$$      endif
c... cp-even H^0 vertex with stop2-loop (MSSM)
c$$$      if (bHloopStop2) then
c$$$         aux(17)=cplModbHstop(2)
c$$$     &   *GGG_vertex(K_ij(0,i2,i4),K_H(0,i1),K(0,i3),J_ij(0,i2,i4),
c$$$     &               eps_st2bH(0,i1),epvec(0,i3))
c$$$         aux(18)=cplModbHstop(2)
c$$$     &   *GGG_vertex(K_ij(0,i2,i4),K_H(0,i1),K(0,i3),J_ij(0,i2,i4),
c$$$     &               eps_st2tH(0,i1),epvec(0,i3))
c$$$      endif
c... cp-even h^0 vertex with sbottom1-loop (MSSM)
c$$$      if (sHloopSbottom1) then
c$$$         aux(19)=cplModsHsbottom(1)
c$$$     &   *GGG_vertex(K_ij(0,i2,i4),K_H(0,i1),K(0,i3),J_ij(0,i2,i4),
c$$$     &               eps_sb1bH(0,i1),epvec(0,i3))
c$$$         aux(20)=cplModsHsbottom(1)
c$$$     &   *GGG_vertex(K_ij(0,i2,i4),K_H(0,i1),K(0,i3),J_ij(0,i2,i4),
c$$$     &               eps_sb1tH(0,i1),epvec(0,i3))
c$$$      endif
c... cp-even h^0 vertex with sbottom2-loop (MSSM)
c$$$      if (sHloopSbottom2) then
c$$$         aux(21)=cplModsHsbottom(2)
c$$$     &   *GGG_vertex(K_ij(0,i2,i4),K_H(0,i1),K(0,i3),J_ij(0,i2,i4),
c$$$     &               eps_sb2bH(0,i1),epvec(0,i3))
c$$$         aux(22)=cplModsHsbottom(2)
c$$$     &   *GGG_vertex(K_ij(0,i2,i4),K_H(0,i1),K(0,i3),J_ij(0,i2,i4),
c$$$     &               eps_sb2tH(0,i1),epvec(0,i3))
c$$$      endif
c... cp-even H^0 vertex with sbottom1-loop (MSSM)
c$$$      if (bHloopSbottom1) then
c$$$         aux(23)=cplModbHsbottom(1)
c$$$     &   *GGG_vertex(K_ij(0,i2,i4),K_H(0,i1),K(0,i3),J_ij(0,i2,i4),
c$$$     &               eps_sb1bH(0,i1),epvec(0,i3))
c$$$         aux(24)=cplModbHsbottom(1)
c$$$     &   *GGG_vertex(K_ij(0,i2,i4),K_H(0,i1),K(0,i3),J_ij(0,i2,i4),
c$$$     &               eps_sb1tH(0,i1),epvec(0,i3))
c$$$      endif
c... cp-even H^0 vertex with sbottom2-loop (MSSM)
c$$$      if (bHloopSbottom2) then
c$$$         aux(25)=cplModbHsbottom(2)
c$$$     &   *GGG_vertex(K_ij(0,i2,i4),K_H(0,i1),K(0,i3),J_ij(0,i2,i4),
c$$$     &               eps_sb2bH(0,i1),epvec(0,i3))
c$$$         aux(26)=cplModbHsbottom(2)
c$$$     &   *GGG_vertex(K_ij(0,i2,i4),K_H(0,i1),K(0,i3),J_ij(0,i2,i4),
c$$$     &               eps_sb2tH(0,i1),epvec(0,i3))
c$$$      endif

      do j=1,26
         MMa=MMa+aux(j)
      enddo
      
      RETURN
      END

***********************************************************************
      FUNCTION MMb(i1,i2,i3,i4)
*----------------------------------------------------------------------
c... author  : Gunnar Klaemke
c... modified: Michael Kubocz
c... date    : 29.08.07 (MK)
*----------------------------------------------------------------------

c     p_i1 --<--*-->-- p_i3  
c               |            
c               O-- H / A0    
c               |           
c     p_i2 --<--*-->-- p_i4

***********************************************************************

      IMPLICIT NONE

#include "ggf_process.inc"

      complex*16 MMb
      complex*16 aux(26)

      complex*16 epvec(0:3,4),J_ij(0:3,4,4)
      real*8  K(0:4,5),K_ij(0:4,4,4),K_H(0:3,4) 
      complex*16 eps_effH(0:3,4),eps_effA(0:3,4),
     &     eps_tH(0:3,4),eps_bH(0:3,4),eps_tA(0:3,4),eps_bA(0:3,4),
     &     eps_st1bH(0:3,4),eps_st1tH(0:3,4),
     &     eps_st2bH(0:3,4),eps_st2tH(0:3,4),
     &     eps_sb1bH(0:3,4),eps_sb1tH(0:3,4),
     &     eps_sb2bH(0:3,4),eps_sb2tH(0:3,4)
*      COMMON /AMPINFO/epvec,K,J_ij,K_ij,K_H,
*     &       eps_effH,eps_effA,eps_tH,eps_bH,eps_tA,eps_bA,
*     &       eps_st1bH,eps_st1tH,eps_st2bH,eps_st2tH,
*     &       eps_sb1bH,eps_sb1tH,eps_sb2bH,eps_sb2tH
      COMMON /AMPINFO/epvec,J_ij,
     &       eps_effH,eps_effA,eps_tH,eps_bH,eps_tA,eps_bA,
     &       eps_st1bH,eps_st1tH,eps_st2bH,eps_st2tH,
     &       eps_sb1bH,eps_sb1tH,eps_sb2bH,eps_sb2tH,K,K_ij,K_H

      integer i1,i2,i3,i4,j

      complex*16 H2G_Loop,A2G_Loop,H2G_effective,A2G_effective


c$$$      EXTERNAL H2G_Loop,A2G_Loop,H2G_effective,A2G_effective,
c$$$     &         H1G_SqLoop,H2G_SqLoop
      EXTERNAL H2G_Loop,A2G_Loop,H2G_effective,A2G_effective

      MMb=DCMPLX(0d0,0d0)
      do j=1,26
         aux(j)=DCMPLX(0d0,0d0)
      enddo
c... cp-even Higgs vertex for effective limes
      if (Heff) then
         aux(1)=cplModHeff*H2G_effective(K_ij(0,i2,i4),K_ij(0,i3,i1),
     &        J_ij(0,i2,i4),J_ij(0,i3,i1))
      endif
c... cp-odd Higgs vertex for effective limes
      if (Aeff) then
         aux(2)=cplModAeff*A2G_effective(K_ij(0,i2,i4),K_ij(0,i3,i1),
     &        J_ij(0,i2,i4),J_ij(0,i3,i1))
      endif
c... cp-even Higgs vertex with top-loop (SM)
      if (HloopTop) then
         aux(3)=cplModHtop*H2G_Loop(K_ij(0,i2,i4),K_ij(0,i3,i1),
     &        J_ij(0,i2,i4),J_ij(0,i3,i1),m_top) 
      endif
c... cp-even Higgs vertex with bottom-loop (SM)
      if (HloopBottom) then
         aux(4)=cplModHbottom*H2G_Loop(K_ij(0,i2,i4),K_ij(0,i3,i1),
     &        J_ij(0,i2,i4),J_ij(0,i3,i1),m_bottom)
      endif
c... cp-odd Higgs vertex with top-loop (MSSM)
      if (AloopTop) then
         aux(5)=cplModAtop*A2G_Loop(K_ij(0,i2,i4),K_ij(0,i3,i1),
     &        J_ij(0,i2,i4),J_ij(0,i3,i1),m_top)  
      endif
c... cp-odd Higgs vertex with bottom-loop (MSSM)
      if (AloopBottom) then
         aux(6)=cplModAbottom*A2G_Loop(K_ij(0,i2,i4),K_ij(0,i3,i1),
     &        J_ij(0,i2,i4),J_ij(0,i3,i1),m_bottom)
      endif
c... cp-even Higgs h vertex with top-loop (MSSM)
      if (sHloopTop) then
         aux(7)=cplModsHtop*H2G_Loop(K_ij(0,i2,i4),K_ij(0,i3,i1),
     &        J_ij(0,i2,i4),J_ij(0,i3,i1),m_top) 
      endif
c... cp-even Higgs H vertex with top-loop (MSSM)
      if (bHloopTop) then
         aux(8)=cplModbHtop*H2G_Loop(K_ij(0,i2,i4),K_ij(0,i3,i1),
     &        J_ij(0,i2,i4),J_ij(0,i3,i1),m_top) 
      endif
c... cp-even Higgs h vertex with bottom-loop (MSSM)
      if (sHloopBottom) then
         aux(9)=cplModsHbottom*H2G_Loop(K_ij(0,i2,i4),K_ij(0,i3,i1),
     &        J_ij(0,i2,i4),J_ij(0,i3,i1),m_bottom)
      endif
c... cp-even Higgs H vertex with bottom-loop (MSSM)
      if (bHloopBottom) then
         aux(10)=cplModbHbottom*H2G_Loop(K_ij(0,i2,i4),K_ij(0,i3,i1),
     &        J_ij(0,i2,i4),J_ij(0,i3,i1),m_bottom)
      endif
c... cp-even h^0 vertex with stop1-loop (MSSM)
c$$$      if (sHloopStop1) then
c$$$         aux(11)=cplModsHstop(1)*H1G_SqLoop(K_ij(0,i2,i4),
c$$$     &        K_ij(0,i3,i1),J_ij(0,i2,i4),J_ij(0,i3,i1),m_stop(1)) 
c$$$         aux(12)=cplModsHstop(1)*H2G_SqLoop(K_ij(0,i2,i4),
c$$$     &        K_ij(0,i3,i1),J_ij(0,i2,i4),J_ij(0,i3,i1),m_stop(1))
c$$$      endif
c... cp-even h^0 vertex with stop2-loop (MSSM)
c$$$      if (sHloopStop2) then
c$$$         aux(13)=cplModsHstop(2)*H1G_SqLoop(K_ij(0,i2,i4),
c$$$     &        K_ij(0,i3,i1),J_ij(0,i2,i4),J_ij(0,i3,i1),m_stop(2)) 
c$$$         aux(14)=cplModsHstop(2)*H2G_SqLoop(K_ij(0,i2,i4),
c$$$     &        K_ij(0,i3,i1),J_ij(0,i2,i4),J_ij(0,i3,i1),m_stop(2))
c$$$      endif
c... cp-even H^0 vertex with stop1-loop (MSSM)
c$$$      if (bHloopStop1) then
c$$$         aux(15)=cplModbHstop(1)*H1G_SqLoop(K_ij(0,i2,i4),
c$$$     &        K_ij(0,i3,i1),J_ij(0,i2,i4),J_ij(0,i3,i1),m_stop(1)) 
c$$$         aux(16)=cplModbHstop(1)*H2G_SqLoop(K_ij(0,i2,i4),
c$$$     &        K_ij(0,i3,i1),J_ij(0,i2,i4),J_ij(0,i3,i1),m_stop(1))
c$$$      endif
c... cp-even H^0 vertex with stop2-loop (MSSM)
c$$$      if (bHloopStop2) then
c$$$         aux(17)=cplModbHstop(2)*H1G_SqLoop(K_ij(0,i2,i4),
c$$$     &        K_ij(0,i3,i1),J_ij(0,i2,i4),J_ij(0,i3,i1),m_stop(2)) 
c$$$         aux(18)=cplModbHstop(2)*H2G_SqLoop(K_ij(0,i2,i4),
c$$$     &        K_ij(0,i3,i1),J_ij(0,i2,i4),J_ij(0,i3,i1),m_stop(2))
c$$$      endif
c...  cp-even h^0 vertex with sbottom1-loop (MSSM)
c$$$      if (sHloopSbottom1) then
c$$$         aux(19)=cplModsHsbottom(1)*H1G_SqLoop(K_ij(0,i2,i4),
c$$$     &     K_ij(0,i3,i1),J_ij(0,i2,i4),J_ij(0,i3,i1),m_sbottom(1)) 
c$$$         aux(20)=cplModsHsbottom(1)*H2G_SqLoop(K_ij(0,i2,i4),
c$$$     &     K_ij(0,i3,i1),J_ij(0,i2,i4),J_ij(0,i3,i1),m_sbottom(1))
c$$$      endif
c... cp-even h^0 vertex with sbottom2-loop (MSSM)
c$$$      if (sHloopSbottom2) then
c$$$         aux(21)=cplModsHsbottom(2)*H1G_SqLoop(K_ij(0,i2,i4),
c$$$     &     K_ij(0,i3,i1),J_ij(0,i2,i4),J_ij(0,i3,i1),m_sbottom(2)) 
c$$$         aux(22)=cplModsHsbottom(2)*H2G_SqLoop(K_ij(0,i2,i4),
c$$$     &     K_ij(0,i3,i1),J_ij(0,i2,i4),J_ij(0,i3,i1),m_sbottom(2))
c$$$      endif
c... cp-even H^0 vertex with sbottom1-loop (MSSM)
c$$$      if (bHloopSbottom1) then
c$$$         aux(23)=cplModbHsbottom(1)*H1G_SqLoop(K_ij(0,i2,i4),
c$$$     &     K_ij(0,i3,i1),J_ij(0,i2,i4),J_ij(0,i3,i1),m_sbottom(1)) 
c$$$         aux(24)=cplModbHsbottom(1)*H2G_SqLoop(K_ij(0,i2,i4),
c$$$     &     K_ij(0,i3,i1),J_ij(0,i2,i4),J_ij(0,i3,i1),m_sbottom(1))
c$$$      endif
c... cp-even H^0 vertex with sbottom2-loop (MSSM)
c$$$      if (bHloopSbottom2) then
c$$$         aux(25)=cplModbHsbottom(2)*H1G_SqLoop(K_ij(0,i2,i4),
c$$$     &     K_ij(0,i3,i1),J_ij(0,i2,i4),J_ij(0,i3,i1),m_sbottom(2)) 
c$$$         aux(26)=cplModbHsbottom(2)*H2G_SqLoop(K_ij(0,i2,i4),
c$$$     &     K_ij(0,i3,i1),J_ij(0,i2,i4),J_ij(0,i3,i1),m_sbottom(2))
c$$$      endif

c... sum up all subamplitudes
      do j=1,26
         MMb=MMb+aux(j)
      enddo 

      RETURN
      END

***********************************************************************
      FUNCTION MMc(i1,i2,i3,i4)
*----------------------------------------------------------------------
c... author  : Gunnar Klaemke
c... modified: Michael Kubocz
c... date    : 29.08.07 (MK)
*----------------------------------------------------------------------

c              H/A0 
c               |
c     p_i1 --<--O-->-- p_i3  
c               |            
c               |           
c     p_i2 --<--*-->-- p_i4

***********************************************************************

      IMPLICIT NONE

#include "ggf_process.inc"

      complex*16 MMc
      complex*16 aux(18)

      complex*16  epvec(0:3,4),J_ij(0:3,4,4)
      real*8  K(0:4,5),K_ij(0:4,4,4),K_H(0:3,4) 
      complex*16 eps_effH(0:3,4),eps_effA(0:3,4),
     &     eps_tH(0:3,4),eps_bH(0:3,4),eps_tA(0:3,4),eps_bA(0:3,4),
     &     eps_st1bH(0:3,4),eps_st1tH(0:3,4),
     &     eps_st2bH(0:3,4),eps_st2tH(0:3,4),
     &     eps_sb1bH(0:3,4),eps_sb1tH(0:3,4),
     &     eps_sb2bH(0:3,4),eps_sb2tH(0:3,4)
*      COMMON /AMPINFO/epvec,K,J_ij,K_ij,K_H,
*     &     eps_effH,eps_effA,eps_tH,eps_bH,eps_tA,eps_bA,
*     &     eps_st1bH,eps_st1tH,eps_st2bH,eps_st2tH,
*     &     eps_sb1bH,eps_sb1tH,eps_sb2bH,eps_sb2tH
      COMMON /AMPINFO/epvec,J_ij,
     &       eps_effH,eps_effA,eps_tH,eps_bH,eps_tA,eps_bA,
     &       eps_st1bH,eps_st1tH,eps_st2bH,eps_st2tH,
     &       eps_sb1bH,eps_sb1tH,eps_sb2bH,eps_sb2tH,K,K_ij,K_H
 
      integer i1,i2,i3,i4,j
      
      complex*16 H3G_Loop,A3G_Loop,H3G_effective,A3G_effective

c$$$  EXTERNAL H3G_Loop,A3G_Loop,H3G_effective,A3G_effective,H3G_SqLoop
      EXTERNAL H3G_Loop,A3G_Loop,H3G_effective,A3G_effective

      MMc=DCMPLX(0d0,0d0)
      do j=1,18
         aux(j)=DCMPLX(0d0,0d0)
      enddo
c... cp-even Higgs vertex in effective limes
      if (Heff) then
         aux(1)=cplModHeff
     &   *H3G_effective(K(0,i1),K(0,i3),K_ij(0,i2,i4),epvec(0,i1),
     &                  epvec(0,i3),J_ij(0,i2,i4))
      endif
c... cp-odd Higgs vertex in effective limes
      if (Aeff) then
         aux(2)=cplModAeff
     &   *A3G_effective(K(0,5),epvec(0,i1),epvec(0,i3),J_ij(0,i2,i4))
      endif
c... cp-even Higgs vertex with top-loop (SM)
      if (HloopTop) then
         aux(3)=cplModHtop*H3G_Loop(K(0,i3),K_ij(0,i2,i4),K(0,i1),
     &        epvec(0,i3),J_ij(0,i2,i4),epvec(0,i1),m_top) 
      endif
c... cp-even Higgs vertex with bottom-loop (SM)
      if (HloopBottom) then
         aux(4)=cplModHbottom*H3G_Loop(K(0,i3),K_ij(0,i2,i4),K(0,i1),
     &        epvec(0,i3),J_ij(0,i2,i4),epvec(0,i1),m_bottom) 
      endif
c... cp-odd Higgs vertex with top-loop (MSSM)
      if (AloopTop) then
         aux(5)=cplModAtop*A3G_Loop(K(0,i1),K(0,i3),K_ij(0,i2,i4),
     &        epvec(0,i1),epvec(0,i3),J_ij(0,i2,i4),m_top)
       endif
c... cp-odd Higgs vertex with bottom-loop (MSSM)
      if (AloopBottom) then
         aux(6)=cplModAbottom*A3G_Loop(K(0,i1),K(0,i3),K_ij(0,i2,i4),
     &        epvec(0,i1),epvec(0,i3),J_ij(0,i2,i4),m_bottom)
      endif
c... cp-even Higgs h vertex with top-loop (MSSM)
      if (sHloopTop) then
         aux(7)=cplModsHtop*H3G_Loop(K(0,i3),K_ij(0,i2,i4),K(0,i1),
     &        epvec(0,i3),J_ij(0,i2,i4),epvec(0,i1),m_top) 
      endif
c... cp-even Higgs H vertex with top-loop (MSSM)
      if (bHloopTop) then
         aux(8)=cplModbHtop*H3G_Loop(K(0,i3),K_ij(0,i2,i4),K(0,i1),
     &        epvec(0,i3),J_ij(0,i2,i4),epvec(0,i1),m_top) 
      endif
c... cp-even Higgs h vertex with bottom-loop (MSSM)
      if (sHloopBottom) then
         aux(9)=cplModsHbottom*H3G_Loop(K(0,i3),K_ij(0,i2,i4),K(0,i1),
     &        epvec(0,i3),J_ij(0,i2,i4),epvec(0,i1),m_bottom) 
      endif
c... cp-even Higgs H vertex with bottom-loop (MSSM)
      if (bHloopBottom) then
         aux(10)=cplModbHbottom*H3G_Loop(K(0,i3),K_ij(0,i2,i4),K(0,i1),
     &        epvec(0,i3),J_ij(0,i2,i4),epvec(0,i1),m_bottom) 
      endif
c... cp-even h^0 vertex with stop1-loop (MSSM)
c$$$      if (sHloopStop1) then
c$$$         aux(11)=cplModsHstop(1)*H3G_SqLoop(K(0,i1),K(0,i3),
c$$$     &        K_ij(0,i2,i4),epvec(0,i1),epvec(0,i3),
c$$$     &        J_ij(0,i2,i4),m_stop(1))
c$$$      endif
c... cp-even h^0 vertex with stop2-loop (MSSM)
c$$$      if (sHloopStop2) then
c$$$         aux(12)=cplModsHstop(2)*H3G_SqLoop(K(0,i1),K(0,i3),
c$$$     &        K_ij(0,i2,i4),epvec(0,i1),epvec(0,i3),
c$$$     &        J_ij(0,i2,i4),m_stop(2))
c$$$      endif
c... cp-even H^0 vertex with stop1-loop (MSSM)
c$$$      if (bHloopStop1) then
c$$$         aux(13)=cplModbHstop(1)*H3G_SqLoop(K(0,i1),K(0,i3),
c$$$     &        K_ij(0,i2,i4),epvec(0,i1),epvec(0,i3),
c$$$     &        J_ij(0,i2,i4),m_stop(1))
c$$$      endif
c... cp-even H^0 vertex with stop2-loop (MSSM)
c$$$      if (bHloopStop2) then
c$$$         aux(14)=cplModbHstop(2)*H3G_SqLoop(K(0,i1),K(0,i3),
c$$$     &        K_ij(0,i2,i4),epvec(0,i1),epvec(0,i3),
c$$$     &        J_ij(0,i2,i4),m_stop(2))
c$$$      endif
c... cp-even h^0 vertex with sbottom1-loop (MSSM)
c$$$      if (sHloopSbottom1) then
c$$$         aux(15)=cplModsHsbottom(1)*H3G_SqLoop(K(0,i1),K(0,i3),
c$$$     &        K_ij(0,i2,i4),epvec(0,i1),epvec(0,i3),
c$$$     &        J_ij(0,i2,i4),m_sbottom(1))
c$$$      endif
c... cp-even h^0 vertex with sbottom2-loop (MSSM)
c$$$      if (sHloopSbottom2) then
c$$$         aux(16)=cplModsHsbottom(2)*H3G_SqLoop(K(0,i1),K(0,i3),
c$$$     &        K_ij(0,i2,i4),epvec(0,i1),epvec(0,i3),
c$$$     &        J_ij(0,i2,i4),m_sbottom(2))
c$$$      endif
c... cp-even H^0 vertex with sbottom1-loop (MSSM)
c$$$      if (bHloopSbottom1) then
c$$$         aux(17)=cplModbHsbottom(1)*H3G_SqLoop(K(0,i1),K(0,i3),
c$$$     &        K_ij(0,i2,i4),epvec(0,i1),epvec(0,i3),
c$$$     &        J_ij(0,i2,i4),m_sbottom(1))
c$$$      endif
c... cp-even H^0 vertex with sbottom2-loop (MSSM)
c$$$      if (bHloopSbottom2) then
c$$$         aux(18)=cplModbHsbottom(2)*H3G_SqLoop(K(0,i1),K(0,i3),
c$$$     &        K_ij(0,i2,i4),epvec(0,i1),epvec(0,i3),
c$$$     &        J_ij(0,i2,i4),m_sbottom(2))
c$$$      endif

c... sum up all subamplitudes
      do j=1,18
         MMc=MMc+aux(j)
      enddo 

      RETURN
      END

***********************************************************************
      FUNCTION MMd(i1,i2,i3,i4)
*----------------------------------------------------------------------
c... author  : Gunnar Klaemke
c... modified: Michael Kubocz
c... date    : 29.08.07 (MK)
*----------------------------------------------------------------------

c                  p_i3
c                   |
c              H/A0 A
c               |   |
c     p_i1 --<--O---*-->-- p_i4  
c                   |
c                   V
c                   |           
c                  p_i2  

***********************************************************************

      IMPLICIT NONE

#include "ggf_process.inc"

      complex*16 MMd
      complex*16 aux(26)

      complex*16 epvec(0:3,4),J_ij(0:3,4,4)
      real*8 K(0:4,5),K_ij(0:4,4,4),K_H(0:3,4) 
      complex*16 eps_effH(0:3,4),eps_effA(0:3,4),
     &     eps_tH(0:3,4),eps_bH(0:3,4),eps_tA(0:3,4),eps_bA(0:3,4),
     &     eps_st1bH(0:3,4),eps_st1tH(0:3,4),
     &     eps_st2bH(0:3,4),eps_st2tH(0:3,4),
     &     eps_sb1bH(0:3,4),eps_sb1tH(0:3,4),
     &     eps_sb2bH(0:3,4),eps_sb2tH(0:3,4)
*      COMMON /AMPINFO/epvec,K,J_ij,K_ij,K_H,
*     &       eps_effH,eps_effA,eps_tH,eps_bH,eps_tA,eps_bA,
*     &       eps_st1bH,eps_st1tH,eps_st2bH,eps_st2tH,
*     &       eps_sb1bH,eps_sb1tH,eps_sb2bH,eps_sb2tH
      COMMON /AMPINFO/epvec,J_ij,
     &       eps_effH,eps_effA,eps_tH,eps_bH,eps_tA,eps_bA,
     &       eps_st1bH,eps_st1tH,eps_st2bH,eps_st2tH,
     &       eps_sb1bH,eps_sb1tH,eps_sb2bH,eps_sb2tH,K,K_ij,K_H

      integer i1,i2,i3,i4,j

      complex*16 dotrc,dotcc

      EXTERNAL dotrc,dotcc
      
      MMd=DCMPLX(0d0,0d0)
      do j=1,26
         aux(j)=DCMPLX(0d0,0d0)
      enddo
c... cp-even Higgs vertex in effective limes
      if (Heff) then
         aux(1)=cplModHeff*(dotcc(eps_effH(0,i1),epvec(0,i2))
     &        *dotcc(epvec(0,i3),epvec(0,i4))-dotcc(eps_effH(0,i1),
     &        epvec(0,i4))*dotcc(epvec(0,i3),epvec(0,i2)))
      endif
c... cp-odd Higgs vertex in effective limes
      if (Aeff) then
         aux(2)=cplModAeff*(dotcc(eps_effA(0,i1),epvec(0,i2))
     &        *dotcc(epvec(0,i3),epvec(0,i4))-dotcc(eps_effA(0,i1),
     &        epvec(0,i4))*dotcc(epvec(0,i3),epvec(0,i2)))
      endif
c... cp-even Higgs vertex with top-loop (SM)
      if (HloopTop) then
         aux(3)=cplModHtop*(dotcc(eps_tH(0,i1),epvec(0,i2))
     &        *dotcc(epvec(0,i3),epvec(0,i4))-dotcc(eps_tH(0,i1),
     &        epvec(0,i4))*dotcc(epvec(0,i3),epvec(0,i2)))
      endif
c... cp-even Higgs vertex with bottom-loop (SM)
      if (HloopBottom) then
         aux(4)=cplModHbottom*(dotcc(eps_bH(0,i1),epvec(0,i2))
     &        *dotcc(epvec(0,i3),epvec(0,i4))-dotcc(eps_bH(0,i1),
     &        epvec(0,i4))*dotcc(epvec(0,i3),epvec(0,i2)))
      endif
c... cp-odd Higgs vertex with top-loop (MSSM)
      if (AloopTop) then
         aux(5)=cplModAtop*(dotcc(eps_tA(0,i1),epvec(0,i2))
     &        *dotcc(epvec(0,i3),epvec(0,i4))-dotcc(eps_tA(0,i1),
     &        epvec(0,i4))*dotcc(epvec(0,i3),epvec(0,i2)))
       endif
c... cp-odd Higgs vertex with bottom-loop (MSSM)
      if (AloopBottom) then
         aux(6)=cplModAbottom*(dotcc(eps_bA(0,i1),epvec(0,i2))
     &        *dotcc(epvec(0,i3),epvec(0,i4))-dotcc(eps_bA(0,i1),
     &        epvec(0,i4))*dotcc(epvec(0,i3),epvec(0,i2)))
      endif
c... cp-even Higgs h vertex with top-loop (MSSM)
      if (sHloopTop) then
         aux(7)=cplModsHtop*(dotcc(eps_tH(0,i1),epvec(0,i2))
     &        *dotcc(epvec(0,i3),epvec(0,i4))-dotcc(eps_tH(0,i1),
     &        epvec(0,i4))*dotcc(epvec(0,i3),epvec(0,i2)))
       endif
c... cp-even Higgs H vertex with top-loop (MSSM)
      if (bHloopTop) then
         aux(8)=cplModbHtop*(dotcc(eps_tH(0,i1),epvec(0,i2))
     &        *dotcc(epvec(0,i3),epvec(0,i4))-dotcc(eps_tH(0,i1),
     &        epvec(0,i4))*dotcc(epvec(0,i3),epvec(0,i2)))
       endif
c... cp-even Higgs h vertex with bottom-loop (MSSM)
      if (sHloopBottom) then
         aux(9)=cplModsHbottom*(dotcc(eps_bH(0,i1),epvec(0,i2))
     &        *dotcc(epvec(0,i3),epvec(0,i4))-dotcc(eps_bH(0,i1),
     &        epvec(0,i4))*dotcc(epvec(0,i3),epvec(0,i2)))
      endif
c... cp-even Higgs H vertex with bottom-loop (MSSM)
      if (bHloopBottom) then
         aux(10)=cplModbHbottom*(dotcc(eps_bH(0,i1),epvec(0,i2))
     &        *dotcc(epvec(0,i3),epvec(0,i4))-dotcc(eps_bH(0,i1),
     &        epvec(0,i4))*dotcc(epvec(0,i3),epvec(0,i2)))
      endif
c... cp-even h^0 vertex with stop1-loop (MSSM)
c$$$      if (sHloopStop1) then
c$$$         aux(11)=cplModsHstop(1)*(dotcc(eps_st1bH(0,i1),
c$$$     &        epvec(0,i2))*dotcc(epvec(0,i3),epvec(0,i4))
c$$$     &        -dotcc(eps_st1bH(0,i1),epvec(0,i4))
c$$$     &        *dotcc(epvec(0,i3),epvec(0,i2)))
c$$$         aux(12)=cplModsHstop(1)*(dotcc(eps_st1tH(0,i1),
c$$$     &        epvec(0,i2))*dotcc(epvec(0,i3),epvec(0,i4))
c$$$     &        -dotcc(eps_st1tH(0,i1),epvec(0,i4))
c$$$     &        *dotcc(epvec(0,i3),epvec(0,i2)))
c$$$      endif
c... cp-even h^0 vertex with stop2-loop (MSSM)
c$$$      if (sHloopStop2) then
c$$$         aux(13)=cplModsHstop(2)*(dotcc(eps_st2bH(0,i1),
c$$$     &        epvec(0,i2))*dotcc(epvec(0,i3),epvec(0,i4))
c$$$     &        -dotcc(eps_st2bH(0,i1),epvec(0,i4))
c$$$     &        *dotcc(epvec(0,i3),epvec(0,i2)))
c$$$         aux(14)=cplModsHstop(2)*(dotcc(eps_st2tH(0,i1),
c$$$     &        epvec(0,i2))*dotcc(epvec(0,i3),epvec(0,i4))
c$$$     &        -dotcc(eps_st2tH(0,i1),epvec(0,i4))
c$$$     &        *dotcc(epvec(0,i3),epvec(0,i2)))
c$$$      endif
c... cp-even H^0 vertex with stop1-loop (MSSM)
c$$$      if (bHloopStop1) then
c$$$         aux(15)=cplModbHstop(1)*(dotcc(eps_st1bH(0,i1),
c$$$     &        epvec(0,i2))*dotcc(epvec(0,i3),epvec(0,i4))
c$$$     &        -dotcc(eps_st1bH(0,i1),epvec(0,i4))
c$$$     &        *dotcc(epvec(0,i3),epvec(0,i2)))
c$$$         aux(16)=cplModbHstop(1)*(dotcc(eps_st1tH(0,i1),
c$$$     &        epvec(0,i2))*dotcc(epvec(0,i3),epvec(0,i4))
c$$$     &        -dotcc(eps_st1tH(0,i1),epvec(0,i4))
c$$$     &        *dotcc(epvec(0,i3),epvec(0,i2)))
c$$$      endif
c... cp-even H^0 vertex with stop2-loop (MSSM)
c$$$      if (bHloopStop2) then
c$$$         aux(17)=cplModbHstop(2)*(dotcc(eps_st2bH(0,i1),
c$$$     &        epvec(0,i2))*dotcc(epvec(0,i3),epvec(0,i4))
c$$$     &        -dotcc(eps_st2bH(0,i1),epvec(0,i4))
c$$$     &        *dotcc(epvec(0,i3),epvec(0,i2)))
c$$$         aux(18)=cplModbHstop(2)*(dotcc(eps_st2tH(0,i1),
c$$$     &        epvec(0,i2))*dotcc(epvec(0,i3),epvec(0,i4))
c$$$     &        -dotcc(eps_st2tH(0,i1),epvec(0,i4))
c$$$     &        *dotcc(epvec(0,i3),epvec(0,i2)))
c$$$      endif
c... cp-even h^0 vertex with sbottom1-loop (MSSM)
c$$$      if (sHloopSbottom1) then
c$$$         aux(19)=cplModsHsbottom(1)*(dotcc(eps_sb1bH(0,i1),
c$$$     &        epvec(0,i2))*dotcc(epvec(0,i3),epvec(0,i4))
c$$$     &        -dotcc(eps_sb1bH(0,i1),epvec(0,i4))
c$$$     &        *dotcc(epvec(0,i3),epvec(0,i2)))
c$$$         aux(20)=cplModsHsbottom(1)*(dotcc(eps_sb1tH(0,i1),
c$$$     &        epvec(0,i2))*dotcc(epvec(0,i3),epvec(0,i4))
c$$$     &        -dotcc(eps_sb1tH(0,i1),epvec(0,i4))
c$$$     &        *dotcc(epvec(0,i3),epvec(0,i2)))
c$$$      endif
c... cp-even h^0 vertex with sbottom2-loop (MSSM)
c$$$      if (sHloopSbottom2) then
c$$$         aux(21)=cplModsHsbottom(2)*(dotcc(eps_sb2bH(0,i1),
c$$$     &        epvec(0,i2))*dotcc(epvec(0,i3),epvec(0,i4))
c$$$     &        -dotcc(eps_sb2bH(0,i1),epvec(0,i4))
c$$$     &        *dotcc(epvec(0,i3),epvec(0,i2)))
c$$$         aux(22)=cplModsHsbottom(2)*(dotcc(eps_sb2tH(0,i1),
c$$$     &        epvec(0,i2))*dotcc(epvec(0,i3),epvec(0,i4))
c$$$     &        -dotcc(eps_sb2tH(0,i1),epvec(0,i4))
c$$$     &        *dotcc(epvec(0,i3),epvec(0,i2)))
c$$$      endif
c... cp-even H^0 vertex with sbottom1-loop (MSSM)
c$$$      if (bHloopSbottom1) then
c$$$         aux(23)=cplModbHsbottom(1)*(dotcc(eps_sb1bH(0,i1),
c$$$     &        epvec(0,i2))*dotcc(epvec(0,i3),epvec(0,i4))
c$$$     &        -dotcc(eps_sb1bH(0,i1),epvec(0,i4))
c$$$     &        *dotcc(epvec(0,i3),epvec(0,i2)))
c$$$         aux(24)=cplModbHsbottom(1)*(dotcc(eps_sb1tH(0,i1),
c$$$     &        epvec(0,i2))*dotcc(epvec(0,i3),epvec(0,i4))
c$$$     &        -dotcc(eps_sb1tH(0,i1),epvec(0,i4))
c$$$     &        *dotcc(epvec(0,i3),epvec(0,i2)))
c$$$      endif
c... cp-even H^0 vertex with sbottom2-loop (MSSM)
c$$$      if (bHloopSbottom2) then
c$$$         aux(25)=cplModbHsbottom(2)*(dotcc(eps_sb2bH(0,i1),
c$$$     &        epvec(0,i2))*dotcc(epvec(0,i3),epvec(0,i4))
c$$$     &        -dotcc(eps_sb2bH(0,i1),epvec(0,i4))
c$$$     &        *dotcc(epvec(0,i3),epvec(0,i2)))
c$$$         aux(26)=cplModbHsbottom(2)*(dotcc(eps_sb2tH(0,i1),
c$$$     &        epvec(0,i2))*dotcc(epvec(0,i3),epvec(0,i4))
c$$$     &        -dotcc(eps_sb2tH(0,i1),epvec(0,i4))
c$$$     &        *dotcc(epvec(0,i3),epvec(0,i2)))
c$$$      endif

c... sum up all subamplitudes
      do j=1,26
         MMd=MMd+aux(j)
      enddo

      RETURN
      END

***********************************************************************
      FUNCTION MMe(i1,i2,i3,i4)
*----------------------------------------------------------------------
c... author  : Gunnar Klaemke  
c... modified: Michael Kubocz
c... date    : 29.08.07 (MK)
*----------------------------------------------------------------------

c           p_i3
c            |
c            |
c   p_i1 ----O---- p_i4  
c           / \
c          /   \
c        p_i2   H         

***********************************************************************         

      IMPLICIT NONE

#include "ggf_process.inc"

      complex*16 MMe
      complex*16 aux(34),pentH_DD_top,pentH_DD_bottom

      complex*16 epvec(0:3,4),J_ij(0:3,4,4) 
      real*8 K(0:4,5),K_ij(0:4,4,4),K_H(0:3,4)
      complex*16 eps_effH(0:3,4),eps_effA(0:3,4),
     &     eps_tH(0:3,4),eps_bH(0:3,4),eps_tA(0:3,4),eps_bA(0:3,4),
     &     eps_st1bH(0:3,4),eps_st1tH(0:3,4),
     &     eps_st2bH(0:3,4),eps_st2tH(0:3,4),
     &     eps_sb1bH(0:3,4),eps_sb1tH(0:3,4),
     &     eps_sb2bH(0:3,4),eps_sb2tH(0:3,4)
      COMMON /AMPINFO/epvec,J_ij,
     &       eps_effH,eps_effA,eps_tH,eps_bH,eps_tA,eps_bA,
     &       eps_st1bH,eps_st1tH,eps_st2bH,eps_st2tH,
     &       eps_sb1bH,eps_sb1tH,eps_sb2bH,eps_sb2tH,K,K_ij,K_H
*      COMMON /AMPINFO/epvec,K,J_ij,K_ij,K_H,
*     &       eps_effH,eps_effA,eps_tH,eps_bH,eps_tA,eps_bA,
*     &       eps_st1bH,eps_st1tH,eps_st2bH,eps_st2tH,
*     &       eps_sb1bH,eps_sb1tH,eps_sb2bH,eps_sb2tH
      
      integer i1,i2,i3,i4,j

      complex*16 dotcc
      complex*16 
     &     A4G_LoopDD

c$$$      EXTERNAL dotcc,H4G_SqTriLoop,H4G_SqBoxLoop,
c$$$     &     H4G_SqLoop,A4G_LoopDD
      EXTERNAL dotcc,A4G_LoopDD

      MMe=DCMPLX(0d0,0d0)
      do j=1,34
         aux(j)=DCMPLX(0d0,0d0)
      enddo
c... cp-even Higgs vertex for effective limes
      if (Heff) then
         aux(1)=cplModHeff*(dotcc(epvec(0,i1),epvec(0,i3))
     &        *dotcc(epvec(0,i2),epvec(0,i4))-dotcc(epvec(0,i1),
     &        epvec(0,i4))*dotcc(epvec(0,i3),epvec(0,i2)))
      endif
c... cp-odd Higgs vertex for effective limes
      if (Aeff) then
         aux(2)=DCMPLX(0.0d0,0.0d0)
      endif
c... cp-even Higgs vertex with top-loop (SM)
      if (HloopTop) then
         CALL H4G_LoopDD(i1,i2,i3,i4,m_top,pentH_DD_top)
         aux(3)=cplModHtop*pentH_DD_top 
      endif
c... cp-even Higgs vertex with bottom-loop (SM)
      if (HloopBottom) then
         CALL H4G_LoopDD(i1,i2,i3,i4,m_bottom,pentH_DD_bottom)
         aux(4)=cplModHbottom*pentH_DD_bottom
      endif
c... cp-odd Higgs vertex with top-loop (MSSM)
      if (AloopTop) then
         aux(5)=cplModAtop*A4G_LoopDD(i1,i2,i3,i4,m_top)
      endif
c... cp-odd Higgs vertex with bottom-loop (MSSM)
      if (AloopBottom) then
         aux(6)=cplModAbottom*A4G_LoopDD(i1,i2,i3,i4,m_bottom)
      endif
c... cp-even Higgs h vertex with top-loop (MSSM)
      if (sHloopTop) then
         CALL H4G_LoopDD(i1,i2,i3,i4,m_top,pentH_DD_top)
         aux(7)=cplModsHtop*pentH_DD_top
      endif
c... cp-even Higgs H vertex with top-loop (MSSM)
      if (bHloopTop) then
         CALL H4G_LoopDD(i1,i2,i3,i4,m_top,pentH_DD_top)
         aux(8)=cplModbHtop*pentH_DD_top
      endif
c... cp-even Higgs h vertex with bottom-loop (MSSM)
      if (sHloopBottom) then
         CALL H4G_LoopDD(i1,i2,i3,i4,m_bottom,pentH_DD_bottom)
         aux(9)=cplModsHbottom*pentH_DD_bottom
      endif
c... cp-even Higgs H vertex with bottom-loop (MSSM)
      if (bHloopBottom) then
         CALL H4G_LoopDD(i1,i2,i3,i4,m_bottom,pentH_DD_bottom)
         aux(10)=cplModbHbottom*pentH_DD_bottom
      endif
c... cp-even h^0 vertex with stop1-loop (MSSM)
c$$$      if (sHloopStop1) then
c$$$       aux(11)=cplModsHstop(1)*H4G_SqTriLoop(i1,i2,i3,i4,m_stop(1))
c$$$       aux(12)=cplModsHstop(1)*H4G_SqBoxLoop(i4,i1,i2,i3,m_stop(1)) 
c$$$       aux(13)=cplModsHstop(1)*H4G_SqLoop(i1,i2,i3,i4,m_stop(1))
c$$$      endif
c... cp-even h^0 vertex with stop2-loop (MSSM)
c$$$      if (sHloopStop2) then
c$$$       aux(14)=cplModsHstop(2)*H4G_SqTriLoop(i1,i2,i3,i4,m_stop(2))
c$$$       aux(15)=cplModsHstop(2)*H4G_SqBoxLoop(i4,i1,i2,i3,m_stop(2))
c$$$       aux(16)=cplModsHstop(2)*H4G_SqLoop(i1,i2,i3,i4,m_stop(2))
c$$$      endif
c... cp-even H^0 vertex with stop1-loop (MSSM)
c$$$      if (bHloopStop1) then
c$$$       aux(17)=cplModbHstop(1)*H4G_SqTriLoop(i1,i2,i3,i4,m_stop(1))
c$$$       aux(18)=cplModbHstop(1)*H4G_SqBoxLoop(i4,i1,i2,i3,m_stop(1)) 
c$$$       aux(19)=cplModbHstop(1)*H4G_SqLoop(i1,i2,i3,i4,m_stop(1))
c$$$      endif
c... cp-even H^0 vertex with stop2-loop (MSSM)
c$$$      if (bHloopStop2) then
c$$$       aux(20)=cplModbHstop(2)*H4G_SqTriLoop(i1,i2,i3,i4,m_stop(2))
c$$$       aux(21)=cplModbHstop(2)*H4G_SqBoxLoop(i4,i1,i2,i3,m_stop(2)) 
c$$$       aux(22)=cplModbHstop(2)*H4G_SqLoop(i1,i2,i3,i4,m_stop(2))
c$$$      endif
c... cp-even h^0 vertex with sbottom1-loop (MSSM)
c$$$      if (sHloopSbottom1) then
c$$$         aux(23)=cplModsHsbottom(1)*H4G_SqTriLoop(i1,i2,i3,i4,
c$$$     &        m_sbottom(1))
c$$$         aux(24)=cplModsHsbottom(1)*H4G_SqBoxLoop(i4,i1,i2,i3,
c$$$     &        m_sbottom(1)) 
c$$$         aux(25)=cplModsHsbottom(1)*H4G_SqLoop(i1,i2,i3,i4,
c$$$     &        m_sbottom(1))
c$$$      endif
c... cp-even h^0 vertex with sbottom2-loop (MSSM)
c$$$      if (sHloopSbottom2) then
c$$$         aux(26)=cplModsHsbottom(2)*H4G_SqTriLoop(i1,i2,i3,i4,
c$$$     &        m_sbottom(2))
c$$$         aux(27)=cplModsHsbottom(2)*H4G_SqBoxLoop(i4,i1,i2,i3,
c$$$     &        m_sbottom(2)) 
c$$$         aux(28)=cplModsHsbottom(2)*H4G_SqLoop(i1,i2,i3,i4,
c$$$     &        m_sbottom(2))
c$$$      endif
c... cp-even H^0 vertex with sbottom1-loop (MSSM)
c$$$      if (bHloopSbottom1) then
c$$$         aux(29)=cplModbHsbottom(1)*H4G_SqTriLoop(i1,i2,i3,i4,
c$$$     &        m_sbottom(1))
c$$$         aux(30)=cplModbHsbottom(1)*H4G_SqBoxLoop(i4,i1,i2,i3,
c$$$     &        m_sbottom(1))
c$$$         aux(31)=cplModbHsbottom(1)*H4G_SqLoop(i1,i2,i3,i4,
c$$$     &        m_sbottom(1))
c$$$      endif
c... cp-even H^0 vertex with sbottom2-loop (MSSM)
c$$$      if (bHloopSbottom2) then
c$$$         aux(32)=cplModbHsbottom(2)*H4G_SqTriLoop(i1,i2,i3,i4,
c$$$     &        m_sbottom(2))
c$$$         aux(33)=cplModbHsbottom(2)*H4G_SqBoxLoop(i4,i1,i2,i3,
c$$$     &        m_sbottom(2)) 
c$$$         aux(34)=cplModbHsbottom(2)*H4G_SqLoop(i1,i2,i3,i4,
c$$$     &        m_sbottom(2))
c$$$      endif

c... sum up all subamplitudes
      do j=1,34
         MMe=MMe+aux(j)
      enddo 
      
      RETURN
      END
