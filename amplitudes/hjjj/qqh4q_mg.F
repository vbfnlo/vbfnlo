c     All q(pbar1)+ q(pbar3) -> q(pbar2)+ q(pbar4) 
c     +q(qbar5)+ qbar(qbar6) + h
c                           |
c                           ---> pbar7 pbar8 
c     Decay products of the higgs are stored in 
c     pbar(0:3,7) and pbar(0:3,8)
c     
c     Neutral current weak boson fusion processes
c
      subroutine qqh4qnc_mg(pbar,sign,
     $     uuccuu,uussuu,ddccuu,ddssuu,
     $     uuccdd,uussdd,ddccdd,ddssdd)
c
      implicit none                
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/mssm.inc"
      double precision  fpi
      parameter ( fpi=4d0*pi)
c
c     Arguments
      real*8 pbar(0:3,6+max_v)
      real*8 uuccuu,uussuu,ddccuu,ddssuu
      real*8 uuccdd,uussdd,ddccdd,ddssdd
      real*8 uuccuu_c(0:2),uussuu_c(0:2),ddccuu_c(0:2),ddssuu_c(0:2)
      real*8 uuccdd_c(0:2),uussdd_c(0:2),ddccdd_c(0:2),ddssdd_c(0:2)
      integer sign(6+max_v)
c
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
      real*8 ph(0:4),lpbar(0:3,6),res(0:2,8)
      integer lfsign(6)
      integer iflav(6,8)
      real*8 p87(0:4),p(0:3,6+max_v),fac
      real*8 p109(0:4)
      integer mu,i,k
      logical lcol
c     flavors of different subprocesses
c                            1 2 3 4 5 6
      DATA(iflav(i,1),i=1,6)/1,1,1,1,1,1/
      DATA(iflav(i,2),i=1,6)/1,1,1,1,2,2/
      DATA(iflav(i,3),i=1,6)/1,1,2,2,1,1/
      DATA(iflav(i,4),i=1,6)/1,1,2,2,2,2/
      DATA(iflav(i,5),i=1,6)/2,2,1,1,1,1/
      DATA(iflav(i,6),i=1,6)/2,2,1,1,2,2/
      DATA(iflav(i,7),i=1,6)/2,2,2,2,1,1/
      DATA(iflav(i,8),i=1,6)/2,2,2,2,2,2/


c
      lcol = .false. 
      goto 20 
      entry qqh4qnc_mg_c(pbar,sign,
     $     uuccuu_c,uussuu_c,ddccuu_c,ddssuu_c,
     $     uuccdd_c,uussdd_c,ddccdd_c,ddssdd_c)
      lcol = .true.
      
 20   continue

c
c     initialize res(k)
      do i = 0,2
         do k=1,8
            res(i,k)=0.0d0
         enddo
      enddo
c     fill local momenta and sign factors
c     diagramatic momenta
      do mu=0,3
         do i=1,6+n_v
            p(mu,i)=sign(i)*pbar(mu,i)
         enddo
         p87(mu) = p(mu,8) - p(mu,7)
         if (n_v.eq.4) then
            p109(mu) = p(mu,10) - p(mu,9)
            ph(mu) = p87(mu) + p109(mu)
         elseif (n_v.eq.2) then
            ph(mu) = p87(mu)
         endif
      enddo
      p87(4)=p87(0)**2-p87(1)**2-p87(2)**2-p87(3)**2
      p109(4)=p109(0)**2-p109(1)**2-p109(2)**2-p109(3)**2
      ph(4) = ph(0)**2 - ph(1)**2 - ph(2)**2 - ph(3)**2 
c
c
      do i=1,6
         lfsign(i)=sign(i)
         do mu=0,3
            lpbar(mu,i)=pbar(mu,i)
         enddo
      enddo
c
c     call madgraph routines
c 
      
      do i=1,8   
c         print*,'iflav=',(iflav(k,i),k=1,6),',i=',i
c         call Suc_ucH21(lpbar,lfsign,iflav(1,i),ph,msq21)
c         print*,'msq21=',msq21
c         call Suc_ucH43(lpbar,lfsign,iflav(1,i),ph,msq43)
c         res(i) = msq21 + msq43
c         print*,'msq43=',msq43
c         print*,'msq21 + msq43=',msq21 + msq43,'='
c         print*,'-------------'
cmr         call sus_usbbh(lpbar,lfsign,iflav(1,i),ph,res(0,i))
         call sus_usuuh(lpbar,lfsign,iflav(1,i),ph,res(0,i))
c         print*,res(i)
c         print*,(msq21 + msq43)/res(i)
c         print*,'-------------'
      enddo
c
c     include higgs decay 
      
C    SUM_pol |D_H(p_8-p_7) * psi-bar(8) psi(7)|^2 * color factor
C    * g_s**2 * production color factor
c      
c
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
c
      do i=0,2                  !loop over color flow
         do k=1,8               !loop over subprocesses
            res(i,k) = res(i,k)*fac
         enddo
      enddo
c     return matrix element squareds (so far just neutral currents)
      if(lcol) then
         do i = 0,2
            uuccuu_c(i) = res(i,1)
            uuccdd_c(i) = res(i,2)
            uussuu_c(i) = res(i,3)
            uussdd_c(i) = res(i,4)
            ddccuu_c(i) = res(i,5)
            ddccdd_c(i) = res(i,6)
            ddssuu_c(i) = res(i,7)
            ddssdd_c(i) = res(i,8)
         enddo

      else
         uuccuu = res(0,1)
         uuccdd = res(0,2)
         uussuu = res(0,3)
         uussdd = res(0,4)
         ddccuu = res(0,5)
         ddccdd = res(0,6)
         ddssuu = res(0,7)
         ddssdd = res(0,8)
      endif
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     All q(pbar1)+ q(pbar3) -> q(pbar2)+ q(pbar4) 
c     +q(qbar5)+ qbar(qbar6) + h
c                           |
c                           ---> pbar7 pbar8 
c     Decay products of the higgs are stored in 
c     pbar(0:3,7) and pbar(0:3,8)
c     
c     charged current weak boson fusion processes
c
      subroutine qqh4qcc_mg(pbar,sign,
     $     udsc,ducs)
c
      implicit none                !        for H -> xx (based on bbar)
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/mssm.inc"
      double precision  fpi
      parameter ( fpi=4d0*pi)
c
c     Arguments
      real*8 pbar(0:3,6+max_v)
      real*8 udsc,ducs
      real*8 udsc_c(0:2),ducs_c(0:2)
      integer sign(6+max_v)
c
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
      real*8 ph(0:4),lpbar(0:3,6),res(0:2,2)
      integer lfsign(6)
      real*8 p87(0:4),p(0:3,6+max_v),fac
      real*8 p109(0:4)
      integer mu,i,k
      logical lcol
c

      lcol = .false.
      goto 20 
c
      entry qqh4qcc_mg_c(pbar,sign,
     $     udsc_c,ducs_c)
      lcol=.true.

 20   continue

c
c     initialize res(k)
      do i = 0,2
         do k=1,2
            res(i,k)=0.0d0
         enddo
      enddo
c     fill local momenta and sign factors
c     diagramatic momenta
      do mu=0,3
         do i=1,6+n_v
            p(mu,i)=sign(i)*pbar(mu,i)
         enddo
         p87(mu) = p(mu,8) - p(mu,7)
         if (n_v.eq.4) then
            p109(mu) = p(mu,10) - p(mu,9)
            ph(mu) = p87(mu) + p109(mu)
         elseif (n_v.eq.2) then
            ph(mu) = p87(mu)
         endif
      enddo
      p87(4)=p87(0)**2-p87(1)**2-p87(2)**2-p87(3)**2
      p109(4)=p109(0)**2-p109(1)**2-p109(2)**2-p109(3)**2
      ph(4) = ph(0)**2 - ph(1)**2 - ph(2)**2 - ph(3)**2 
c
      do i=1,6
         lfsign(i)=sign(i)
         do mu=0,3
            lpbar(mu,i)=pbar(mu,i)
         enddo
      enddo
     
c
c     call madgraph routines
c 
      
      do i=1,2   
         call sdc_usbbh(lpbar,lfsign,ph,res(0,i))
      enddo
c
c     include higgs decay 
      
C    SUM_pol |D_H(p_8-p_7) * psi-bar(8) psi(7)|^2 * color factor
C    * g_s**2 * production color factor
c    
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
c
      do i=0,2
         do k=1,2               !loop over subprocesses
            res(i,k) = res(i,k)*fac
         enddo
      enddo
c     return matrix element squareds (charged currents)

      if(lcol) then
         do i=0,2
            udsc_c(i) = res(i,1)
            ducs_c(i) = res(i,2)
         enddo
      else
         udsc = res(0,1)
         ducs = res(0,2)
      endif
c      
      return
      end
ccc
