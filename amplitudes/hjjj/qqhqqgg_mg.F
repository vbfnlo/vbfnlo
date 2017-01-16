c     All q(pbar1)+ q(pbar3) -> q(pbar2)+ q(pbar4) 
c     +g(qbar1)+ g(qbar2) + h
c                           |
c                           ---> pbar5 pbar6 
c     Decay products of the higgs are stored in 
c     pbar(0:3,5) and pbar(0:3,6)
c     
c
      subroutine qqhqqgg_mg(pbar,sign,qbar,gsign,
     $     uucc,uuss,ddcc,ddss,udsc,ducs)
c
      implicit none                !        for H -> xx (based on bbar)
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/mssm.inc"
      double precision  fpi
      parameter ( fpi=4d0*pi)
c     switch for interference  terms
      logical lintOff,lcol
      common/interference/ lintOff
c
c     Arguments
      real*8 pbar(0:3,4+max_v),qbar(0:4,2)
      real*8 uucc,uuss,ddcc,ddss,udsc,ducs
      real*8 uucc_c(0:6),uuss_c(0:6),ddcc_c(0:6),
     $     ddss_c(0:6),udsc_c(0:6),ducs_c(0:6)
      integer sign(4+max_v),gsign(2)
c and which are calculated from output of KOPPLN
c
      double precision  clr, xm2, xmg, b, v, a
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     1                  V(4,5),A(4,5)
#include "VBFNLO/utilities/koppln.inc"
c CPS scheme
      double precision qgammaq
      common /VBFNLO_HIGGS_CPS/ qgammaq
c     
c     Local variables
      real*8 ph(0:4),lpbar(0:3,4),lqbar(0:3,2),resnc(4,0:6),rescc(0:6)
      real*8 res(6,0:6)
      integer lfsign(4),lgsign(2)
      real*8 p65(0:4),p87(0:4),p(0:3,4+max_v),fac

      integer mu,i,k


c
      lcol =.false.
      goto 20

      entry qqhqqgg_mg_c(pbar,sign,qbar,gsign,
     $     uucc_c,uuss_c,ddcc_c,ddss_c,udsc_c,ducs_c)
      lcol = .true.

 20   continue
c
c     switch of interference terms
c
      lintOff = .true.
c
c
c     initialize res(k)
      do i = 0,6
         do k=1,6
            res(k,i)=0d0
         enddo
      enddo
c     fill local momenta and sign factors
c     diagramatic momenta
      do mu = 0,3
         do i = 1,4+n_v
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo         
         p65(mu) = p(mu,6) - p(mu,5)
         if (n_v.eq.4) then
            p87(mu) = p(mu,8) - p(mu,7)
            ph(mu) = p65(mu) + p87(mu)
         elseif (n_v.eq.2) then
            ph(mu) = p65(mu)
         endif
      enddo
      p65(4) = p65(0)**2 - p65(1)**2 - p65(2)**2 - p65(3)**2
      if (n_v.eq.4) p87(4) = 
     1         p87(0)**2 - p87(1)**2 - p87(2)**2 - p87(3)**2
      ph(4) = ph(0)**2 - ph(1)**2 - ph(2)**2 - ph(3)**2
c
      do i=1,4
         lfsign(i)=sign(i)
         do mu=0,3
            lpbar(mu,i)=pbar(mu,i)
         enddo
      enddo
      do i=1,2
         lgsign(i)=gsign(i)
         do mu=0,3
            lqbar(mu,i)=qbar(mu,i)
         enddo
      enddo
c
c     call madgraph routines
c 
      
c      call SUC_UCGGH(lpbar,lfsign,lqbar,lgsign,ph,res(1))
c      call SUS_USGGH(lpbar,lfsign,lqbar,lgsign,ph,res(2))
c      call Sdc_dcGGH(lpbar,lfsign,lqbar,lgsign,ph,res(3))
c      call SdS_dSGGH(lpbar,lfsign,lqbar,lgsign,ph,res(4))
      call SQQ_QQGGH_NC(lpbar,lfsign,lqbar,lgsign,ph,resnc)
      call Sus_dcggH(lpbar,lfsign,lqbar,lgsign,ph,rescc)

c      call Sdc_usggH(lpbar,lfsign,lqbar,lgsign,ph,res(6))
c
c     include higgs decay 
      
C    SUM_pol |D_H(p_6-p_5) * psi-bar(6) psi(5)|^2 * color factor
C    * g_s**2 * production color factor

      fac = 1.d0 / ( (ph(4)-xm2(6))**2 + xmg(6)**2 ) !propagator^2 only
c
c multiply by 16 pi*M Gamma to compensate for the phase space integral in
c LIPSN
c
c   Int dq^2/2pi * dOmega/(32 pi^2) * 1/(q^2-M*2)^2+(MGamma)^2) 
c
c   = 1/(16 pi MGamma)
c
c (where the propagator factor is fac above)
      fac = fac *16*pi*qgammaq

      do i=0,6                  !loop over color flows
         do k=1,4               !loop over subprocesses
            res(k,i) = resnc(k,i)*fac
         enddo
         res(5,i) = rescc(i)*fac
         res(6,i) = rescc(i)*fac
      enddo

c     return matrix element squareds

c      print*,'qqhqqgg_mg: res=',res
      if(lcol) then
         do i = 0,6
            uucc_c(i) = res(1,i)
            uuss_c(i) = res(2,i)
            ddcc_c(i) = res(3,i)
            ddss_c(i) = res(4,i)
            udsc_c(i) = res(5,i)
            ducs_c(i) = res(6,i)
         enddo

      else
         uucc = res(1,0)
         uuss = res(2,0)
         ddcc = res(3,0)
         ddss = res(4,0)
         udsc = res(5,0)
         ducs = res(6,0)
      endif
      return
      end
