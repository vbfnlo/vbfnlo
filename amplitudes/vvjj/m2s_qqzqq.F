c
c  Dieter Zeppenfeld, <dieter@pheno.physics.wisc.edu>
c  Initial version:  2003 February 12
c  Last modified by Barbara Jager:  2006 July 11
c  Further modificatin by Julien Baglio: 2012 June 10
c  to implement HH production
c  Last modified by Julien Baglio:  2014 Oct. 24
c
c  this function computes |ME|^2 for qq->qq e+ve mu-vm~ ("WW")
c	or qq->qq e+e-mu+mu- or qq->qq e+e- vm vm~ ("ZZ")
c       or qq->qq e+e- A or qq->qq ve ve~ A ("ZA")
c       or qq->qq HH ("HH")
c	at LO and virtual corrections via call of "qqwwqq.f" or "qqzzqq.f"
c       or "qqzaqq.f" or "qqhhqq.f"
c
c	for qq->qqZZ random helicity summation is employed for decay leptons
c
c ------------------------------------------------------------------------
c
      real*8 function m2s_wbfvv(
     &                   bos,   !in:  Boson identifier, 1=A, 2 = Z, 6=H, 34=W+W-
     &                   nlo,   !in:  Born only for nlo=0, +interference for 1
     &                   xi,	!in:  Feynman x parameters of incoming fermions
     &                   p,	!in:  momenta associated with external fermions
     &                   v,	!in:  momenta associated with Z/H dec. fermions
     &                   rn,     !in:  random number array
     &                   ps_number  )
          use globalvars, only: ldoblha
          use VBFNLO_HIGGSTO, only: VBFNLO_HIGGSTO_gridHt
      implicit none
c
c declare input/output variables
c
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "tensor.inc"
#include "tensorz.inc"
#include "tensorsp2aa.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/tau_hel.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/process.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

      integer bos, nlo, ps_number

      real*8 xi(nx), p(0:3,max_p), v(0:3,max_v), m2s
      real*8 rn !, intent(in)
      real*8 vtemp(0:3,max_v)

c  helicity selection
      integer h
      integer jsig, jsig1, jsig3, jsig5
      common /chelsum/ jsig,jsig1,jsig3,jsig5

      complex*16 zero
      parameter (zero=(0d0,0d0))
c
c declare external functions
c
      integer FL_ZH
      double complex dotrc
      external FL_ZH, dotrc
c
c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"

      real*8 q12(0:4,3), q34(0:4,3), qvv(0:3)!, qww(0:3)

c
c declare local variables
c
      real*8 q_sf, qq

      integer init/0/, I, J, mu
      save init
c
c declare variables for summation process
c
      INTEGER ires,lsign(3:4)

      integer FSIGN(4+max_v), physToDiag(4), nmin, nmax, nproc(4)
      save nmax, nproc

c variables for hadronic decays
      integer N_gen_W
      real*8 fac_W, mjj2
      external mjj2
      integer N_gen_up, N_gen_down
      real*8 fac_Z_up, fac_Z_down

c
c store contributions from subprocess ip in res(ip) where
      real*8 res(maxnumsubproc)
c
c declare parton distribution variables
c
      real*8 pdf(-6:6,2)
c
c define program switches
c
      logical ldebug
      parameter (ldebug=.false.)		!output debug information

      common /hcount / h
      real*8 weight,rnumb,RandomNumber

      integer a,b
      double complex aazz1(0:3,0:3,3), 
     1                azzz1(0:3,0:3,3), zazz1(0:3,0:3,3), 
     2                zzzz1(0:3,0:3,3), wwzz51(0:3,0:3,3), 
     3                wwzz61(0:3,0:3,3)

      double complex aazz2(0:3,0:3,3), 
     1                azzz2(0:3,0:3,3), zazz2(0:3,0:3,3), 
     2                zzzz2(0:3,0:3,3), wwzz52(0:3,0:3,3), 
     3                wwzz62(0:3,0:3,3)

      double complex aazz3(0:3,0:3,3), 
     1                azzz3(0:3,0:3,3), zazz3(0:3,0:3,3), 
     2                zzzz3(0:3,0:3,3), wwzz53(0:3,0:3,3), 
     3                wwzz63(0:3,0:3,3)

      double complex aazz4(0:3,0:3,3), 
     1                azzz4(0:3,0:3,3), zazz4(0:3,0:3,3), 
     2                zzzz4(0:3,0:3,3), wwzz54(0:3,0:3,3), 
     3                wwzz64(0:3,0:3,3)


*** start sophy qed: parton distribution variables
      real*8 Qpdf(-7:7,2,3),zz_qed,x1,x2
#include "VBFNLO/utilities/Apdf.inc"
*** end sophy qed

c CPS scheme
      DOUBLE PRECISION CLR,XM2,XMG,Bcoupl,Vcoupl,Acoupl
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),Bcoupl(6,6,6),Vcoupl(4,5),Acoupl(4,5)
      double precision qinv
      double precision qgammaq
      common /VBFNLO_HIGGS_CPS/ qgammaq

c -----------------------------------------------------------------------
      if ( (init .eq. 0) .and. (.not. ldoblha) ) then
         photon_hel=0
         if (bos.eq.2) then
            write(6,*) " "
            write(6,*) "Zjj amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
         elseif (bos.eq.1) then
            write(6,*) " "
            write(6,*) "Ajj amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
         elseif (bos.eq.6) then
            write(6,*) " "
            write(6,*) "Hjj amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
         elseif (bos.eq.34 .or. bos.eq.43) then
            write(6,*) " "
            write(6,*) "W+W- jj amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
            call vtoww_reset
         elseif (bos.eq.22) then
            write(6,*) " "
            write(6,*) "ZZ jj -> 4l amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
            call vtozz_reset
         elseif (bos.eq.21) then
            write(6,*) " "
            write(6,*) "ZZ jj -> 2l 2v amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
            call vtozz_reset
         elseif (bos.eq.211) then
            write(6,*) " "
            write(6,*) "ZA jj -> 2l A amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
            call vtozz_reset
         elseif (bos.eq.212) then
            write(6,*) " "
            write(6,*) "ZA jj -> 2v A amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
            call vtozz_reset
         elseif (bos.eq.61 .or. bos.eq.16) then
            write(6,*) " "
            write(6,*) "HAjj amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
         elseif (bos.eq.66) then
            write(6,*) " "
            write(6,*) "HHjj amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
            call vtohh_reset
         elseif (bos.eq.11) then
            write(6,*) " "
            write(6,*) "VV (-> Spin2) -> A A amplitude square info:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
            call vtosp2aa_reset
         elseif (bos.eq.734) then
           if (with_spin2) then  
            write(6,*) " "
            write(6,*) "VV (-> Spin2) -> W+W- amplitude square info:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
           else
            write(6,*) " "
            write(6,*) "VV -> Higgs (not Spin2!) -> W+W- amplitude square info:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
           endif
            call vtosp2ww_reset
         elseif (bos.eq.722) then
           if (with_spin2) then 
            write(6,*) " "
            write(6,*) "VV (-> Spin2) -> ZZ jj -> 4l"
            write(6,*) "amplitude square information:"
            write(6,*) " "
           else
            write(6,*) " "
            write(6,*) "VV -> Higgs (not Spin2!) -> ZZ jj -> 4l"
            write(6,*) "amplitude square information:"
            write(6,*) " "
           endif
            call vtosp2zz_reset
         elseif (bos.eq.721) then
           if (with_spin2) then 
            write(6,*) " "
            write(6,*) "VV (-> Spin2) -> ZZ jj ->  2l 2v"
            write(6,*) "amplitude square information:"
            write(6,*) " "
           else
            write(6,*) " "
            write(6,*) "VV -> Higgs (not Spin2!) -> ZZ jj ->  2l 2v"
            write(6,*) "amplitude square information:"
            write(6,*) " "
           endif
            call vtosp2zz_reset
         endif
         call printnfl(.true.)
         h = 1
         init = 1
         do i = 1,maxnumsubproc
            res(i) = 0
         enddo
      else if ( (init .eq. 0) .and. ldoblha ) then
        photon_hel=0
        h = 1
        init = 1
      else ! init .ne. 0
         do i = 1,nmax
            res(i) = 0
         enddo
      endif    ! end of init
c -----------------------------------------------------------------------
      call Calc_Momentum_Transfer(p, v, q12,q34,1)

c
c for WW->4l precalculate A->WW,Z->WW,AZ->WW etc leptonic tensors
c
c      if (.false.) then
      if (bos.eq.34 .or. bos.eq.43) then !  34 = W+W-
c reset lfs to .true. to force recalculation of virtual contributions in qqwwqq
         do i = 1,4
            lfs(i) = .true.
         enddo

c...Les Houches interface         
         if ((lha.or.hepmc).and..not.doNLO) then
            helicity(1)=-1
            helicity(2)= 1
            helicity(3)= 1
            helicity(4)=-1
         endif


c lepton spinors and W+- polarization vectors
         CALL IXXXXX(v(0,2),ZERO ,1,-1,wep) !W(1,5))          !e+       
         CALL OXXXXX(v(0,1),ZERO ,-1,1,wve) !W(1,6))          !ve 
         CALL OXXXXX(v(0,4),ZERO ,-1,1,wmu) !W(1,7))          !mu-      
         CALL IXXXXX(v(0,3),ZERO ,1,-1,wvm) !W(1,8))          !vm~
         CALL JIOXXX(wep,wve,GWF ,WMASS,WWIDTH,wp) !W(1,10))    !W+
         CALL JIOXXX(wvm,wmu,GWF ,WMASS,WWIDTH,wm) !W(1,11))    !W-

         do mu = 0,3
            qp(mu) = v(mu,1)+v(mu,2)
            qm(mu) = v(mu,3)+v(mu,4)
            qww(mu) = qp(mu) + qm(mu)
         enddo
         qp(4) = qp(0)**2-qp(1)**2-qp(2)**2-qp(3)**2
         qm(4) = qm(0)**2-qm(1)**2-qm(2)**2-qm(3)**2

         if (abs(nlo).eq.5 .or. nlo.eq.-4) then     ! determine W+- 
                                                    ! polarization vectors 
c wpp(mu) = wp(mu) - xp*qp(mu) and wmp(mu) = wm(mu) - xm*qm(mu)
            qq = qp(0)*qww(0)-qp(1)*qww(1)-qp(2)*qww(2)-qp(3)*qww(3)
            xp = dotrc(qww,wp)/qq
            qq = qm(0)*qww(0)-qm(1)*qww(1)-qm(2)*qww(2)-qm(3)*qww(3)
            xm = dotrc(qww,wm)/qq
            do mu = 0,3
               wpp(mu+1) = wp(mu+1) - xp*qp(mu)
               wmp(mu+1) = wm(mu+1) - xm*qm(mu)
            enddo
            do mu = 5,6
               wpp(mu) = wp(mu)
               wmp(mu) = wm(mu)
            enddo
         elseif (nlo.eq.1) then
            do mu = 1,6
               wpp(mu) = wp(mu)
               wmp(mu) = wm(mu)
            enddo
         endif   ! end of nlo = -4 or +-5

c leptonic tensors
C for W+W-
         j = 1

         SELECT CASE(procid)
         CASE(WPhadWMjj, WPWMhadjj)            ! semileptonic decay modes

            call anomal_formfactor(q12(0,j),q34(0,j),qp(0),qm(0))
            call atoww_had(v,aww)
            call ztoww_had(v,zww)
            call aatoww_had(q12(0,j),q34(0,j),v,aaww(0,0,j))
            call aztoww_had(q12(0,j),q34(0,j),v,azww(0,0,j))
            call aztoww_had(q34(0,j),q12(0,j),v,zaww(0,0,j))
            call zztoww_had(q12(0,j),q34(0,j),v,zzww(0,0,j))
            call wwtoww_had(q12(0,j),q34(0,j),v,wwww6(0,0,j)) ! q12 = W-
            call wwtoww_had(q34(0,j),q12(0,j),v,wwww5(0,0,j)) ! q12 = W+

            call WVtoWP_had(2,q34(0,j),v,NCwpa(0,0,1,j),NCwpz(0,0,1,j)) !emit W- on upper
            call WVtoWP_had(2,q12(0,j),v,NCwpa(0,0,2,j),NCwpz(0,0,2,j)) !emit W- on lower
            call WVtoWP_had(1,q34(0,j),v,CCwpa(0,0,1,j),CCwpz(0,0,1,j)) !emit W- on upper
            call WVtoWP_had(1,q12(0,j),v,CCwpa(0,0,2,j),CCwpz(0,0,2,j)) !emit W- on lower
C     for WV --> mu- nu_mu
            call WVtoWM_had(2,q34(0,j),v,NCwma(0,0,1,j),NCwmz(0,0,1,j)) !emit W+ on upper
            call WVtoWM_had(2,q12(0,j),v,NCwma(0,0,2,j),NCwmz(0,0,2,j)) !emit W+ on lower
            call WVtoWM_had(1,q34(0,j),v,CCwma(0,0,1,j),CCwmz(0,0,1,j)) !emit W+ on upper
            call WVtoWM_had(1,q12(0,j),v,CCwma(0,0,2,j),CCwmz(0,0,2,j)) !emit W+ on lower


         CASE DEFAULT         ! leptonic decay modes

         if (with_anom) then ! anomalous gauge boson couplings
c           using global form factor for all tensors of one phase space point
c           this ensures proper cancellations for anomalous contributions
c           energy scale is invariant WW mass
            call anomal_formfactor(q12(0,j),q34(0,j),qp(0),qm(0))
            call atoww_anomal2(v,aww)
            call ztoww_anomal2(v,zww)
            call aatoww_anomal2(q12(0,j),q34(0,j),v,aaww(0,0,j))
            call aztoww_anomal2(q12(0,j),q34(0,j),v,azww(0,0,j))
            call aztoww_anomal2(q34(0,j),q12(0,j),v,zaww(0,0,j))
            call zztoww_anomal2(q12(0,j),q34(0,j),v,zzww(0,0,j))
            call wwtoww_anomal2(q12(0,j),q34(0,j),v,wwww6(0,0,j)) ! q12 = W-
            call wwtoww_anomal2(q34(0,j),q12(0,j),v,wwww5(0,0,j)) ! q12 = W+

            call WVtoWP_anomal2(2,q34(0,j),v,NCwpa(0,0,1,j),NCwpz(0,0,1,j)) !emit W- on upper
            call WVtoWP_anomal2(2,q12(0,j),v,NCwpa(0,0,2,j),NCwpz(0,0,2,j)) !emit W- on lower
            call WVtoWP_anomal2(1,q34(0,j),v,CCwpa(0,0,1,j),CCwpz(0,0,1,j)) !emit W- on upper
            call WVtoWP_anomal2(1,q12(0,j),v,CCwpa(0,0,2,j),CCwpz(0,0,2,j)) !emit W- on lower
C     for WV --> mu- nu_mu
            call WVtoWM_anomal2(2,q34(0,j),v,NCwma(0,0,1,j),NCwmz(0,0,1,j)) !emit W+ on upper
            call WVtoWM_anomal2(2,q12(0,j),v,NCwma(0,0,2,j),NCwmz(0,0,2,j)) !emit W+ on lower
            call WVtoWM_anomal2(1,q34(0,j),v,CCwma(0,0,1,j),CCwmz(0,0,1,j)) !emit W+ on upper
            call WVtoWM_anomal2(1,q12(0,j),v,CCwma(0,0,2,j),CCwmz(0,0,2,j)) !emit W+ on lower

         elseif (with_kk) then  ! kaluza-klein scenario
#ifdef WITH_KK
            call atoww_kk(v,aww)
            call ztoww_kk(v,zww)
            call aatoww_kk(q12(0,j),q34(0,j),v,aaww(0,0,j))
            call aztoww_kk(q12(0,j),q34(0,j),v,azww(0,0,j))
            call aztoww_kk(q34(0,j),q12(0,j),v,zaww(0,0,j))
            call zztoww_kk(q12(0,j),q34(0,j),v,zzww(0,0,j))
            call wwtoww_kk(q12(0,j),q34(0,j),v,wwww6(0,0,j)) ! q12 = W-
            call wwtoww_kk(q34(0,j),q12(0,j),v,wwww5(0,0,j)) ! q12 = W+

            call WVtoWP_kk(2,q34(0,j),v,NCwpa(0,0,1,j),NCwpz(0,0,1,j)) !emit W- on upper
            call WVtoWP_kk(2,q12(0,j),v,NCwpa(0,0,2,j),NCwpz(0,0,2,j)) !emit W- on lower
            call WVtoWP_kk(1,q34(0,j),v,CCwpa(0,0,1,j),CCwpz(0,0,1,j)) !emit W- on upper
            call WVtoWP_kk(1,q12(0,j),v,CCwpa(0,0,2,j),CCwpz(0,0,2,j)) !emit W- on lower
C     for WV --> mu- nu_mu
            call WVtoWM_kk(2,q34(0,j),v,NCwma(0,0,1,j),NCwmz(0,0,1,j)) !emit W+ on upper
            call WVtoWM_kk(2,q12(0,j),v,NCwma(0,0,2,j),NCwmz(0,0,2,j)) !emit W+ on lower
            call WVtoWM_kk(1,q34(0,j),v,CCwma(0,0,1,j),CCwmz(0,0,1,j)) !emit W+ on upper
            call WVtoWM_kk(1,q12(0,j),v,CCwma(0,0,2,j),CCwmz(0,0,2,j)) !emit W+ on lower
#endif

         elseif (with_spin2) then  
            call atoww(v,aww)
            call ztoww(v,zww)
                  
            call aatoww_spin2(q12(0,j),q34(0,j),v,aaww(0,0,j))
            call aztoww_spin2(q12(0,j),q34(0,j),v,azww(0,0,j))
            call aztoww_spin2(q34(0,j),q12(0,j),v,zaww(0,0,j))
            call zztoww_spin2(q12(0,j),q34(0,j),v,zzww(0,0,j))
            call wwtoww_spin2(q12(0,j),q34(0,j),v,wwww6(0,0,j)) ! q12 = W-
            call wwtoww_spin2(q34(0,j),q12(0,j),v,wwww5(0,0,j)) ! q12 = W+
C for WV --> e+ nu_e
c NCw tensors for NC process (k=1...4), CCw for CC (k=5,6)
            call WVtoWP(2,q34(0,j),v,NCwpa(0,0,1,j),NCwpz(0,0,1,j)) !emit W- on upper
            call WVtoWP(2,q12(0,j),v,NCwpa(0,0,2,j),NCwpz(0,0,2,j)) !emit W- on lower
            call WVtoWP(1,q34(0,j),v,CCwpa(0,0,1,j),CCwpz(0,0,1,j)) !emit W- on upper
            call WVtoWP(1,q12(0,j),v,CCwpa(0,0,2,j),CCwpz(0,0,2,j)) !emit W- on lower
C     for WV --> mu- nu_mu
            call WVtoWM(2,q34(0,j),v,NCwma(0,0,1,j),NCwmz(0,0,1,j)) !emit W+ on upper
            call WVtoWM(2,q12(0,j),v,NCwma(0,0,2,j),NCwmz(0,0,2,j)) !emit W+ on lower
            call WVtoWM(1,q34(0,j),v,CCwma(0,0,1,j),CCwmz(0,0,1,j)) !emit W+ on upper
            call WVtoWM(1,q12(0,j),v,CCwma(0,0,2,j),CCwmz(0,0,2,j)) !emit W+ on lower
         else                   !(SM)
            call atoww(v,aww)
c         z1 = qww(0)*aww(0) - qww(1)*aww(1) - qww(2)*aww(2) - qww(3)*aww(3)
c         print*," gauge check A: ",
c     1   z1/(qww(0)*max(abs(aww(0)),abs(aww(1)),abs(aww(2)),abs(aww(3))))
            call ztoww(v,zww)
c         z1 = qww(0)*zww(0) - qww(1)*zww(1) - qww(2)*zww(2) - qww(3)*zww(3)
c         print*," gauge check Z: ",
c     1   z1/(qww(0)*max(abs(zww(0)),abs(zww(1)),abs(zww(2)),abs(zww(3))))
                  
            call aatoww(q12(0,j),q34(0,j),v,aaww(0,0,j))
            call aztoww(q12(0,j),q34(0,j),v,azww(0,0,j))
            call aztoww(q34(0,j),q12(0,j),v,zaww(0,0,j))
            call zztoww(q12(0,j),q34(0,j),v,zzww(0,0,j))
            call wwtoww(q12(0,j),q34(0,j),v,wwww6(0,0,j)) ! q12 = W-
            call wwtoww(q34(0,j),q12(0,j),v,wwww5(0,0,j)) ! q12 = W+
C for WV --> e+ nu_e
c NCw tensors for NC process (k=1...4), CCw for CC (k=5,6)
            call WVtoWP(2,q34(0,j),v,NCwpa(0,0,1,j),NCwpz(0,0,1,j)) !emit W- on upper
            call WVtoWP(2,q12(0,j),v,NCwpa(0,0,2,j),NCwpz(0,0,2,j)) !emit W- on lower
            call WVtoWP(1,q34(0,j),v,CCwpa(0,0,1,j),CCwpz(0,0,1,j)) !emit W- on upper
            call WVtoWP(1,q12(0,j),v,CCwpa(0,0,2,j),CCwpz(0,0,2,j)) !emit W- on lower
C     for WV --> mu- nu_mu
            call WVtoWM(2,q34(0,j),v,NCwma(0,0,1,j),NCwmz(0,0,1,j)) !emit W+ on upper
            call WVtoWM(2,q12(0,j),v,NCwma(0,0,2,j),NCwmz(0,0,2,j)) !emit W+ on lower
            call WVtoWM(1,q34(0,j),v,CCwma(0,0,1,j),CCwmz(0,0,1,j)) !emit W+ on upper
            call WVtoWM(1,q12(0,j),v,CCwma(0,0,2,j),CCwmz(0,0,2,j)) !emit W+ on lower
         endif   ! end of WW->4l leptonic tensors(for loop for kk,anom,spin2,sm)

         END SELECT      ! case for semileptonic decays

c     
c-----------------------------
c

c for ZZ->4l precalculate  leptonic tensors
      elseif (bos.eq.22) then
c reset lzs to .true. to force recalculation of virtual contributions in qqzzqq
         do i = 1,4
            lzs(i) = .true.
         enddo
c lepton spinors and Z polarization vectors ! 
c	lepton helicities not fixed 
c 	-> sum over all possible helicites in |M|**2
c
c select helicity: h ... random number for lepton helicity
c			 combination (h=1:4) .. needed only for ZZ

         if (ldoblha) then
           h = blha_ranhelcomb
         endif
         h = mod(h,4) + 1
         
         ie = sign(1,2-h)
         iu = (-1)**(h+1)

c...Les Houches interface
         if ((lha.or.hepmc).and..not.doNLO) then
            helicity(1)= ie
            helicity(2)=-ie
            helicity(3)= iu
            helicity(4)=-iu
         endif

         CALL IXXXXX(v(0,1),ZERO ,+ie,-1,lep)         !e+
         CALL OXXXXX(v(0,2),ZERO ,-ie, 1,lem)         !e- 
         CALL IXXXXX(v(0,3),ZERO ,+iu,-1,lup)         !mu+
         CALL OXXXXX(v(0,4),ZERO ,-iu, 1,lum)         !mu- 

         SELECT CASE(procid)
         CASE(ZZhadjj)
           ! for general output to all combinations: up-type first
           if ((finalquarks(1).eq.93 .or. finalquarks(1).eq.94)) then
             ! set couplings for Z hadronic decay into uubar
             call setZtouu
           endif
           CALL JIOXXX(lep,lem,GZ_ZF ,ZMASS,ZWIDTH,ze)    !Zl
           CALL JIOXXX(lep,lem,GZ_AF ,ZERO ,ZERO  ,ae)    !Al

         CASE DEFAULT
           CALL JIOXXX(lep,lem,GZL ,ZMASS,ZWIDTH,ze)    !Zl
           CALL JIOXXX(lep,lem,GAL ,ZERO ,ZERO  ,ae)    !Al

         END SELECT

         CALL JIOXXX(lup,lum,GZL ,ZMASS,ZWIDTH,zu)    !Zl
         CALL JIOXXX(lup,lum,GAL ,ZERO ,ZERO  ,au)    !Al

         do mu = 0,3
            qe(mu) = v(mu,1)+v(mu,2)
            qu(mu) = v(mu,3)+v(mu,4)
            qvv(mu) = qe(mu) + qu(mu)
         enddo

         qe(4) = qe(0)**2-qe(1)**2-qe(2)**2-qe(3)**2
         qu(4) = qu(0)**2-qu(1)**2-qu(2)**2-qu(3)**2

         if (abs(nlo).eq.5 .or. nlo.eq.-4) then  ! determine V/V 
                                                 ! polarization vectors 
c zet(mu) = ze(mu) - xxez*qe(mu) and zut(mu) = zu(mu) - xxmz*qu(mu)
            qq = qe(0)*qvv(0)-qe(1)*qvv(1)-qe(2)*qvv(2)-qe(3)*qvv(3)
            xxez = dotrc(qvv,ze)/qq
            xxea = dotrc(qvv,ae)/qq
            
	    qq = qu(0)*qvv(0)-qu(1)*qvv(1)-qu(2)*qvv(2)-qu(3)*qvv(3)
            xxuz = dotrc(qvv,zu)/qq            
            xxua = dotrc(qvv,au)/qq
            do mu = 0,3
               zet(mu+1) = ze(mu+1) - xxez*qe(mu)
               zut(mu+1) = zu(mu+1) - xxuz*qu(mu)
               
	       aet(mu+1) = ae(mu+1) - xxea*qe(mu)
               aut(mu+1) = au(mu+1) - xxua*qu(mu)
            enddo
            do mu = 5,6
               zet(mu) = ze(mu)
               zut(mu) = zu(mu)

               aet(mu) = ae(mu)
               aut(mu) = au(mu)
            enddo
         elseif (nlo.eq.1) then
            do mu = 1,6
               zet(mu) = ze(mu)
               zut(mu) = zu(mu)

               aet(mu) = ae(mu)
               aut(mu) = au(mu)
            enddo
         endif  ! end of nlo = -4 or +-5

c leptonic tensors
         j = 1


         SELECT CASE(procid)
         CASE(ZZhadjj)

           call anomal_formfactor(q12(0,j),q34(0,j),qe(0),qu(0))

           ! for AA,AZ,ZA,ZZ -> e+ e- mu+ mu- :
           call vvtozz_had(q12(0,j),q34(0,j),v,h,aazz(0,0,j),azzz(0,0,j),zazz(0,0,j),zzzz(0,0,j))
           ! for W+W-
           call wwtozz_had(q12(0,j),q34(0,j),v,h,wwzz5(0,0,j))  ! q12 = W+
           call wwtozz_had(q34(0,j),q12(0,j),v,h,wwzz6(0,0,j))  ! q12 = W-
           ! for V1V2 --> e+ e- or  V1V2 -->  mu+ mu- (Vi=A,Z)
           ! NCl tensors for NC process (k=1...4)
           call VVtoll_had(2,1,h,q34(0,j),v,aaee(0,0,1,j),azee(0,0,1,j),zaee(0,0,1,j),zzee(0,0,1,j))     !emit V on upper    
           call VVtoll_had(1,1,h,q12(0,j),v,aaee(0,0,2,j),azee(0,0,2,j),zaee(0,0,2,j),zzee(0,0,2,j))     !emit V on lower
           call VVtoll_had(2,2,h,q34(0,j),v,aauu(0,0,1,j),azuu(0,0,1,j),zauu(0,0,1,j),zzuu(0,0,1,j))     !emit V on upper     
           call VVtoll_had(1,2,h,q12(0,j),v,aauu(0,0,2,j),azuu(0,0,2,j),zauu(0,0,2,j),zzuu(0,0,2,j))     !emit V on lower
           ! for W+W- --> e+ e- or  W+W- -->  mu+ mu-
           ! CCl tensors for CC (k=5)
           call WWtoll_had(2,1,h,q34(0,j),v,CCee(0,0,1,j))     !emit V on upper
           call WWtoll_had(1,1,h,q12(0,j),v,CCee(0,0,2,j))     !emit V on lower
           call WWtoll_had(2,2,h,q34(0,j),v,CCuu(0,0,1,j))     !emit V on upper
           call WWtoll_had(1,2,h,q12(0,j),v,CCuu(0,0,2,j))     !emit V on lower
           ! for W-W+ --> e+ e- or mu+ mu-
           ! CCl for CC (k=6)
           call WWtoll_had(1,1,h,q34(0,j),v,CCee6(0,0,1,j))     !emit V on upper
           call WWtoll_had(2,1,h,q12(0,j),v,CCee6(0,0,2,j))     !emit V on lower
           call WWtoll_had(1,2,h,q34(0,j),v,CCuu6(0,0,1,j))     !emit V on upper
           call WWtoll_had(2,2,h,q12(0,j),v,CCuu6(0,0,2,j))     !emit V on lower
           ! for V-> llll
           call vto4l_had(v,h,azz,zzztens)


         CASE DEFAULT

C for AA,AZ,ZA,ZZ -> e+ e- mu+ mu- :
         if (with_kk) then      ! Kaluza-Klein scenario
#ifdef WITH_KK
            call vvtozz_KK(q12(0,j),q34(0,j),v,h,aazz(0,0,j),azzz(0,0,j),
     - 					   zazz(0,0,j),zzzz(0,0,j))
c
C for W+W-
            call wwtozz_KK(q12(0,j),q34(0,j),v,h,wwzz5(0,0,j)) ! q12 = W+
            call wwtozz_KK(q34(0,j),q12(0,j),v,h,wwzz6(0,0,j)) ! q12 = W-

C for V1V2 --> e+ e- or  V1V2 -->  mu+ mu- (Vi=A,Z)
c NCl tensors for NC process (k=1...4)

            call VVtoll_KK(2,1,h,q34(0,j),v,aaee(0,0,1,j),azee(0,0,1,j),
     -	 	 zaee(0,0,1,j),zzee(0,0,1,j))    	 !emit V on upper       
            call VVtoll_KK(1,1,h,q12(0,j),v,aaee(0,0,2,j),azee(0,0,2,j),
     -		 zaee(0,0,2,j),zzee(0,0,2,j))    	 !emit V on lower
            call VVtoll_KK(2,2,h,q34(0,j),v,aauu(0,0,1,j),azuu(0,0,1,j),
     -		 zauu(0,0,1,j),zzuu(0,0,1,j))    	 !emit V on upper       
            call VVtoll_KK(1,2,h,q12(0,j),v,aauu(0,0,2,j),azuu(0,0,2,j),
     -		zauu(0,0,2,j),zzuu(0,0,2,j))     	!emit V on lower

C for W+W- --> e+ e- or  W+W- -->  mu+ mu-
c  CCl tensors for CC (k=5)

            call WWtoll_KK(2,1,h,q34(0,j),v,CCee(0,0,1,j)) !emit V on upper
            call WWtoll_KK(1,1,h,q12(0,j),v,CCee(0,0,2,j)) !emit V on lower
            call WWtoll_KK(2,2,h,q34(0,j),v,CCuu(0,0,1,j)) !emit V on upper
            call WWtoll_KK(1,2,h,q12(0,j),v,CCuu(0,0,2,j)) !emit V on lower
	 
C for W-W+ --> e+ e- or mu+ mu-
c CCl for CC (k=6)
	
            call WWtoll_KK(1,1,h,q34(0,j),v,CCee6(0,0,1,j)) !emit V on upper
            call WWtoll_KK(2,1,h,q12(0,j),v,CCee6(0,0,2,j)) !emit V on lower
            call WWtoll_KK(1,2,h,q34(0,j),v,CCuu6(0,0,1,j)) !emit V on upper
            call WWtoll_KK(2,2,h,q12(0,j),v,CCuu6(0,0,2,j)) !emit V on lower
	 
	 
c for V-> llll
	
            call vto4l_KK(v,h,azz,zzztens)	
#endif	   
    
         elseif (with_spin2) then
            call vvtozz_spin2(q12(0,j),q34(0,j),v,h,aazz(0,0,j),azzz(0,0,j),
     - 					   zazz(0,0,j),zzzz(0,0,j))
c
C for W+W-
            call wwtozz_spin2(q12(0,j),q34(0,j),v,h,wwzz5(0,0,j)) ! q12 = W+
            call wwtozz_spin2(q34(0,j),q12(0,j),v,h,wwzz6(0,0,j)) ! q12 = W-

c other processes do not change

C for V1V2 --> e+ e- or  V1V2 -->  mu+ mu- (Vi=A,Z)
c NCl tensors for NC process (k=1...4)

            call VVtoll(2,1,h,q34(0,j),v,aaee(0,0,1,j),azee(0,0,1,j),
     #	 	 zaee(0,0,1,j),zzee(0,0,1,j))    	 !emit V on upper        
            call VVtoll(1,1,h,q12(0,j),v,aaee(0,0,2,j),azee(0,0,2,j),
     #		 zaee(0,0,2,j),zzee(0,0,2,j))    	 !emit V on lower
            call VVtoll(2,2,h,q34(0,j),v,aauu(0,0,1,j),azuu(0,0,1,j),
     #		 zauu(0,0,1,j),zzuu(0,0,1,j))    	 !emit V on upper        
            call VVtoll(1,2,h,q12(0,j),v,aauu(0,0,2,j),azuu(0,0,2,j),
     #		zauu(0,0,2,j),zzuu(0,0,2,j))     	!emit V on lower

C for W+W- --> e+ e- or  W+W- -->  mu+ mu-
c  CCl tensors for CC (k=5)

            call WWtoll(2,1,h,q34(0,j),v,CCee(0,0,1,j)) !emit V on upper
            call WWtoll(1,1,h,q12(0,j),v,CCee(0,0,2,j)) !emit V on lower
            call WWtoll(2,2,h,q34(0,j),v,CCuu(0,0,1,j)) !emit V on upper
            call WWtoll(1,2,h,q12(0,j),v,CCuu(0,0,2,j)) !emit V on lower
	 
C for W-W+ --> e+ e- or mu+ mu-
c CCl for CC (k=6)
	
            call WWtoll(1,1,h,q34(0,j),v,CCee6(0,0,1,j)) !emit V on upper
            call WWtoll(2,1,h,q12(0,j),v,CCee6(0,0,2,j)) !emit V on lower
            call WWtoll(1,2,h,q34(0,j),v,CCuu6(0,0,1,j)) !emit V on upper
            call WWtoll(2,2,h,q12(0,j),v,CCuu6(0,0,2,j)) !emit V on lower
	 
	 
c for V-> llll
	
            call vto4l(v,h,azz,zzztens)


         elseif (with_anom) then ! anomalous gauge boson couplings
c           using global form factor for all tensors of one phase space point
c           this ensures proper cancellations for anomalous contributions
c           energy scale is invariant ZZ mass        
           call anomal_formfactor(q12(0,j),q34(0,j),qe(0),qu(0))


           call vvtozz_anomal(q12(0,j),q34(0,j),v,h,aazz(0,0,j),azzz(0,0,j),
     # 					   zazz(0,0,j),zzzz(0,0,j))
c
C for W+W-
           call wwtozz_anomal(q12(0,j),q34(0,j),v,h,wwzz5(0,0,j))  ! q12 = W+
           call wwtozz_anomal(q34(0,j),q12(0,j),v,h,wwzz6(0,0,j))  ! q12 = W-

C for V1V2 --> e+ e- or  V1V2 -->  mu+ mu- (Vi=A,Z)
c NCl tensors for NC process (k=1...4)

           call VVtoll_anomal(2,1,h,q34(0,j),v,aaee(0,0,1,j),azee(0,0,1,j),
     #	 	 zaee(0,0,1,j),zzee(0,0,1,j))    	 !emit V on upper        
	   call VVtoll_anomal(1,1,h,q12(0,j),v,aaee(0,0,2,j),azee(0,0,2,j),
     #		 zaee(0,0,2,j),zzee(0,0,2,j))    	 !emit V on lower
           call VVtoll_anomal(2,2,h,q34(0,j),v,aauu(0,0,1,j),azuu(0,0,1,j),
     #		 zauu(0,0,1,j),zzuu(0,0,1,j))    	 !emit V on upper        
	   call VVtoll_anomal(1,2,h,q12(0,j),v,aauu(0,0,2,j),azuu(0,0,2,j),
     #		zauu(0,0,2,j),zzuu(0,0,2,j))     	!emit V on lower

C for W+W- --> e+ e- or  W+W- -->  mu+ mu-
c  CCl tensors for CC (k=5)

           call WWtoll_anomal(2,1,h,q34(0,j),v,CCee(0,0,1,j))     !emit V on upper
           call WWtoll_anomal(1,1,h,q12(0,j),v,CCee(0,0,2,j))     !emit V on lower
           call WWtoll_anomal(2,2,h,q34(0,j),v,CCuu(0,0,1,j))     !emit V on upper
           call WWtoll_anomal(1,2,h,q12(0,j),v,CCuu(0,0,2,j))     !emit V on lower
	 
C for W-W+ --> e+ e- or mu+ mu-
c CCl for CC (k=6)
	
           call WWtoll_anomal(1,1,h,q34(0,j),v,CCee6(0,0,1,j))     !emit V on upper
           call WWtoll_anomal(2,1,h,q12(0,j),v,CCee6(0,0,2,j))     !emit V on lower
           call WWtoll_anomal(1,2,h,q34(0,j),v,CCuu6(0,0,1,j))     !emit V on upper
           call WWtoll_anomal(2,2,h,q12(0,j),v,CCuu6(0,0,2,j))     !emit V on lower
	 
	 
c for V-> llll
	
           call vto4l_anomal(v,h,azz,zzztens)


c-----------------------------
c
         else                   ! SM
            call vvtozz(q12(0,j),q34(0,j),v,h,aazz(0,0,j),azzz(0,0,j),
     - 					   zazz(0,0,j),zzzz(0,0,j))
c   

C for W+W-
            call wwtozz(q12(0,j),q34(0,j),v,h,wwzz5(0,0,j)) ! q12 = W+
            call wwtozz(q34(0,j),q12(0,j),v,h,wwzz6(0,0,j)) ! q12 = W-

C for V1V2 --> e+ e- or  V1V2 -->  mu+ mu- (Vi=A,Z)
c NCl tensors for NC process (k=1...4)

            call VVtoll(2,1,h,q34(0,j),v,aaee(0,0,1,j),azee(0,0,1,j),
     -	 	 zaee(0,0,1,j),zzee(0,0,1,j))    	 !emit V on upper       
            call VVtoll(1,1,h,q12(0,j),v,aaee(0,0,2,j),azee(0,0,2,j),
     -		 zaee(0,0,2,j),zzee(0,0,2,j))    	 !emit V on lower
            call VVtoll(2,2,h,q34(0,j),v,aauu(0,0,1,j),azuu(0,0,1,j),
     -		 zauu(0,0,1,j),zzuu(0,0,1,j))    	 !emit V on upper       
            call VVtoll(1,2,h,q12(0,j),v,aauu(0,0,2,j),azuu(0,0,2,j),
     -		zauu(0,0,2,j),zzuu(0,0,2,j))     	!emit V on lower

C for W+W- --> e+ e- or  W+W- -->  mu+ mu-
c  CCl tensors for CC (k=5)

            call WWtoll(2,1,h,q34(0,j),v,CCee(0,0,1,j)) !emit V on upper
            call WWtoll(1,1,h,q12(0,j),v,CCee(0,0,2,j)) !emit V on lower
            call WWtoll(2,2,h,q34(0,j),v,CCuu(0,0,1,j)) !emit V on upper
            call WWtoll(1,2,h,q12(0,j),v,CCuu(0,0,2,j)) !emit V on lower
	 
C for W-W+ --> e+ e- or mu+ mu-
c CCl for CC (k=6)
	
            call WWtoll(1,1,h,q34(0,j),v,CCee6(0,0,1,j)) !emit V on upper
            call WWtoll(2,1,h,q12(0,j),v,CCee6(0,0,2,j)) !emit V on lower
            call WWtoll(1,2,h,q34(0,j),v,CCuu6(0,0,1,j)) !emit V on upper
            call WWtoll(2,2,h,q12(0,j),v,CCuu6(0,0,2,j)) !emit V on lower
	 
	 
c for V-> llll
	
            call vto4l(v,h,azz,zzztens)	
	     
c-----------------------------
c     
         endif                  !end of withkk/spin2/sm for bos = 22

         END SELECT


c for VV (-> spin-2) -> gamma gamma  precalculate leptonic tensors
      elseif (bos.eq.11) then

c     polarization vectors ! 
c
c select helicity: h ... random number for lepton helicity
c			 combination (h=1:4) .. needed only for ZZ

         if (ldoblha) then
           h = blha_ranhelcomb
         endif
         h = mod(h,4) + 1

         lsign(3) = sign(1,2-h)
         lsign(4) = (-1)**(h+1)

         if ((lha.or.hepmc).and..not.doNLO) then
            helicity(1)= lsign(3)
            helicity(2)= lsign(4)
         endif

         CALL VXXXXX(v(0,1),ZERO ,lsign(3),1,sp2au(1)) !photon 1
         CALL VXXXXX(v(0,2),ZERO ,lsign(4),1,sp2ae(1)) !photon 2

c leptonic tensors
         j=1
!          CALL wwsp2toaa(q34(0,j),q12(0,j),v,lsign,sp2wwaa5(0,0,j))
!          CALL wwsp2toaa(q12(0,j),q34(0,j),v,lsign,sp2wwaa6(0,0,j))
!          CALL zzsp2toaa(q34(0,j),q12(0,j),v,lsign,sp2zzaa(0,0,j))
!          CALL azsp2toaa(q34(0,j),q12(0,j),v,lsign,sp2azaa(0,0,j))
!          CALL azsp2toaa(q12(0,j),q34(0,j),v,lsign,sp2zaaa(0,0,j))
!          CALL aasp2toaa(q34(0,j),q12(0,j),v,lsign,sp2aaaa(0,0,j))

       call vvsp2tovv(1,q34(0,j),q12(0,j),v(0,1),v(0,2),sp2au(1),sp2ae(1), 
     1           sp2wwaa5(0,0,j),sp2wwaa6(0,0,j),sp2zzaa(0,0,j),
     2           sp2azaa(0,0,j),sp2zaaa(0,0,j),sp2aaaa(0,0,j))
         

c for VV -> spin-2 -> WW -> 4l precalculate leptonic tensors
       elseif (bos.eq.734) then !  734 = spin-2 -> W+W-
         do i = 1,4
            lfs(i) = .true.
         enddo

c...Les Houches interface         
         if ((lha.or.hepmc).and..not.doNLO) then
            helicity(1)=-1
            helicity(2)= 1
            helicity(3)= 1
            helicity(4)=-1
         endif


c lepton spinors and W+- polarization vectors
         CALL IXXXXX(v(0,2),ZERO ,1,-1,wep) !W(1,5))          !e+       
         CALL OXXXXX(v(0,1),ZERO ,-1,1,wve) !W(1,6))          !ve 
         CALL OXXXXX(v(0,4),ZERO ,-1,1,wmu) !W(1,7))          !mu-      
         CALL IXXXXX(v(0,3),ZERO ,1,-1,wvm) !W(1,8))          !vm~
         CALL JIOXXX(wep,wve,GWF ,WMASS,WWIDTH,wp) !W(1,10))    !W+
         CALL JIOXXX(wvm,wmu,GWF ,WMASS,WWIDTH,wm) !W(1,11))    !W-

         do mu = 0,3
            qp(mu) = v(mu,1)+v(mu,2)
            qm(mu) = v(mu,3)+v(mu,4)
            qww(mu) = qp(mu) + qm(mu)
         enddo
         qp(4) = qp(0)**2-qp(1)**2-qp(2)**2-qp(3)**2
         qm(4) = qm(0)**2-qm(1)**2-qm(2)**2-qm(3)**2

         if (abs(nlo).eq.5 .or. nlo.eq.-4) then     ! determine W+- 
                                                    ! polarization vectors 
c wpp(mu) = wp(mu) - xp*qp(mu) and wmp(mu) = wm(mu) - xm*qm(mu)
            qq = qp(0)*qww(0)-qp(1)*qww(1)-qp(2)*qww(2)-qp(3)*qww(3)
            xp = dotrc(qww,wp)/qq
            qq = qm(0)*qww(0)-qm(1)*qww(1)-qm(2)*qww(2)-qm(3)*qww(3)
            xm = dotrc(qww,wm)/qq
            do mu = 0,3
               wpp(mu+1) = wp(mu+1) - xp*qp(mu)
               wmp(mu+1) = wm(mu+1) - xm*qm(mu)
            enddo
            do mu = 5,6
               wpp(mu) = wp(mu)
               wmp(mu) = wm(mu)
            enddo
         elseif (nlo.eq.1) then
            do mu = 1,6
               wpp(mu) = wp(mu)
               wmp(mu) = wm(mu)
            enddo
         endif   ! end of nlo = -4 or +-5

c leptonic tensors
C for spin-2 -> W+W-
         j = 1
         if (with_spin2) then  

          call vvsp2tovv(2,q12(0,j),q34(0,j),qp,qm,wp,wm, 
     1           sp2wwww6(0,0,j),sp2wwww5(0,0,j),sp2zzww(0,0,j),
     2           sp2azww(0,0,j),sp2zaww(0,0,j),sp2aaww(0,0,j))

         else  ! Higgs instead of spin-2 resonance     
           
            call aatoww_Hres(q12(0,j),q34(0,j),v,sp2aaww(0,0,j))
            call aztoww_Hres(q12(0,j),q34(0,j),v,sp2azww(0,0,j))
            call aztoww_Hres(q34(0,j),q12(0,j),v,sp2zaww(0,0,j))
            call zztoww_Hres(q12(0,j),q34(0,j),v,sp2zzww(0,0,j))
            call wwtoww_Hres(q12(0,j),q34(0,j),v,sp2wwww6(0,0,j)) ! q12 = W-
            call wwtoww_Hres(q34(0,j),q12(0,j),v,sp2wwww5(0,0,j)) ! q12 = W+

         endif

c for  VV -> spin-2 -> ZZ -> 4l precalculate  leptonic tensors
      elseif (bos.eq.722) then
c reset lzs to .true. to force recalculation of virtual contributions in qqsp2zzqq
         do i = 1,4
            lzs(i) = .true.
         enddo
c lepton spinors and Z polarization vectors ! 
c	lepton helicities not fixed 
c 	-> sum over all possible helicites in |M|**2
c
c select helicity: h ... random number for lepton helicity
c			 combination (h=1:4) .. needed only for ZZ

         if (ldoblha) then
           h = blha_ranhelcomb
         endif
         h = mod(h,4) + 1
         
         ie = sign(1,2-h)
         iu = (-1)**(h+1)

c...Les Houches interface
         if ((lha.or.hepmc).and..not.doNLO) then
            helicity(1)= ie
            helicity(2)=-ie
            helicity(3)= iu
            helicity(4)=-iu
         endif
		
         CALL IXXXXX(v(0,1),ZERO ,+ie,-1,lep)	      !e+	
         CALL OXXXXX(v(0,2),ZERO ,-ie, 1,lem)	      !e- 
         CALL IXXXXX(v(0,3),ZERO ,+iu,-1,lup)	      !mu+	
         CALL OXXXXX(v(0,4),ZERO ,-iu, 1,lum)	      !mu- 
	
	 CALL JIOXXX(lep,lem,GZL ,ZMASS,ZWIDTH,ze)     !Zl
	 CALL JIOXXX(lep,lem,GAL ,ZERO ,ZERO  ,ae)      !Al
         CALL JIOXXX(lup,lum,GZL ,ZMASS,ZWIDTH,zu)     !Zl
	 CALL JIOXXX(lup,lum,GAL ,ZERO ,ZERO  ,au)      !Al

         do mu = 0,3
            qe(mu) = v(mu,1)+v(mu,2)
            qu(mu) = v(mu,3)+v(mu,4)
            qvv(mu) = qe(mu) + qu(mu)
         enddo
	 	 
         qe(4) = qe(0)**2-qe(1)**2-qe(2)**2-qe(3)**2
         qu(4) = qu(0)**2-qu(1)**2-qu(2)**2-qu(3)**2

         if (abs(nlo).eq.5 .or. nlo.eq.-4) then  ! determine V/V 
                                                 ! polarization vectors 
c zet(mu) = ze(mu) - xxez*qe(mu) and zut(mu) = zu(mu) - xxmz*qu(mu)
            qq = qe(0)*qvv(0)-qe(1)*qvv(1)-qe(2)*qvv(2)-qe(3)*qvv(3)
            xxez = dotrc(qvv,ze)/qq
            xxea = dotrc(qvv,ae)/qq
            
	    qq = qu(0)*qvv(0)-qu(1)*qvv(1)-qu(2)*qvv(2)-qu(3)*qvv(3)
            xxuz = dotrc(qvv,zu)/qq            
            xxua = dotrc(qvv,au)/qq
            do mu = 0,3
               zet(mu+1) = ze(mu+1) - xxez*qe(mu)
               zut(mu+1) = zu(mu+1) - xxuz*qu(mu)
               
	       aet(mu+1) = ae(mu+1) - xxea*qe(mu)
               aut(mu+1) = au(mu+1) - xxua*qu(mu)
            enddo
            do mu = 5,6
               zet(mu) = ze(mu)
               zut(mu) = zu(mu)

               aet(mu) = ae(mu)
               aut(mu) = au(mu)
            enddo
         elseif (nlo.eq.1) then
            do mu = 1,6
               zet(mu) = ze(mu)
               zut(mu) = zu(mu)

               aet(mu) = ae(mu)
               aut(mu) = au(mu)
            enddo
         endif  ! end of nlo = -4 or +-5

c leptonic tensors
       j = 1
       if (with_spin2) then  

         call vvsp2tovv(3,q12(0,j),q34(0,j),qe,qu,ze,zu, 
     1           wwzz51(0,0,j),wwzz61(0,0,j),zzzz1(0,0,j),
     2           azzz1(0,0,j),zazz1(0,0,j),aazz1(0,0,j))

         call vvsp2tovv(1,q12(0,j),q34(0,j),qe,qu,ae,au, 
     1           wwzz52(0,0,j),wwzz62(0,0,j),zzzz2(0,0,j),
     2           azzz2(0,0,j),zazz2(0,0,j),aazz2(0,0,j))

         call vvsp2tovv(4,q12(0,j),q34(0,j),qe,qu,ze,au, 
     1           wwzz53(0,0,j),wwzz63(0,0,j),zzzz3(0,0,j),
     2           azzz3(0,0,j),zazz3(0,0,j),aazz3(0,0,j))

         call vvsp2tovv(4,q12(0,j),q34(0,j),qe,qu,ae,zu, 
     1           wwzz54(0,0,j),wwzz64(0,0,j),zzzz4(0,0,j),
     2           azzz4(0,0,j),zazz4(0,0,j),aazz4(0,0,j))

         do a=0,3
          do b=0,3
      sp2wwzz6(a,b,j)=wwzz61(a,b,j)+wwzz62(a,b,j)+wwzz63(a,b,j)+wwzz64(a,b,j)
      sp2wwzz5(a,b,j)=wwzz51(a,b,j)+wwzz52(a,b,j)+wwzz53(a,b,j)+wwzz54(a,b,j)
      sp2zzzz(a,b,j)=zzzz1(a,b,j)+zzzz2(a,b,j)+zzzz3(a,b,j)+zzzz4(a,b,j)
      sp2azzz(a,b,j)=azzz1(a,b,j)+azzz2(a,b,j)+azzz3(a,b,j)+azzz4(a,b,j)
      sp2zazz(a,b,j)=zazz1(a,b,j)+zazz2(a,b,j)+zazz3(a,b,j)+zazz4(a,b,j)
      sp2aazz(a,b,j)=aazz1(a,b,j)+aazz2(a,b,j)+aazz3(a,b,j)+aazz4(a,b,j)
          enddo
         enddo

       else  ! Higgs instead of spin-2 resonance     

       call vvtozz_Hres(q12(0,j),q34(0,j),v,h,sp2aazz(0,0,j),sp2azzz(0,0,j),
     - 					   sp2zazz(0,0,j),sp2zzzz(0,0,j))

            call wwtozz_Hres(q12(0,j),q34(0,j),v,h,sp2wwzz5(0,0,j)) ! q12 = W+
            call wwtozz_Hres(q34(0,j),q12(0,j),v,h,sp2wwzz6(0,0,j)) ! q12 = W-

       endif


c for  VV -> spin-2 -> ZZ -> 2l 2v precalculate  leptonic tensors
      elseif (bos.eq.721) then
c reset lzs to .true. to force recalculation of virtual contributions in qqzzqq
         do i = 1,4
            lzs(i) = .true.
         enddo
c lepton spinors and Z polarization vectors ! 
c	lepton helicities not fixed 
c 		-> sum over all possible helicites in |M|**2
c
c select helicity: h ... random number for lepton helicity
c			 combination (h=1:4), but only h=1 or 3 
c	   	         contribute for ZZ -> 2l 2v (no righthanded neutrinos)

         if (ldoblha) then
           h = 2*blha_ranhelcomb-1
         endif
         h = mod(h+2,4)
         
         ie = sign(1,2-h)

c...Les Houches interface
         if ((lha.or.hepmc).and..not.doNLO) then
            helicity(1)= ie
            helicity(2)=-ie
            helicity(3)=-1
            helicity(4)= 1
         endif

		
         CALL IXXXXX(v(0,1),ZERO ,+ie,-1,lep) !e+   
         CALL OXXXXX(v(0,2),ZERO ,-ie, 1,lem) !e- 
         CALL OXXXXX(v(0,3),ZERO ,-1,+1,lum) !vm	  
         CALL IXXXXX(v(0,4),ZERO ,+1,-1,lup) !vm~	
         
         CALL JIOXXX(lep,lem,GZL ,ZMASS,ZWIDTH,ze) !Ze
         CALL JIOXXX(lep,lem,GAL ,ZERO ,ZERO  ,ae) !Ae
         CALL JIOXXX(lup,lum,GZN ,ZMASS,ZWIDTH,zu) !Zv
         
         do i = 1,4
            au(i) = 0d0
            zu(i) = -zu(i)
         enddo   			
         do i = 5,6    !add momentum info to A->vv current (au vanishes anyway) 
            au(i) = zu(i)
         enddo    	
	
         do mu = 0,3
            qe(mu) = v(mu,1)+v(mu,2)
            qu(mu) = v(mu,3)+v(mu,4)
            qvv(mu) = qe(mu) + qu(mu)
         enddo
         
         qe(4) = qe(0)**2-qe(1)**2-qe(2)**2-qe(3)**2
         qu(4) = qu(0)**2-qu(1)**2-qu(2)**2-qu(3)**2

         if (abs(nlo).eq.5 .or. nlo.eq.-4) then ! determine V/V 
                                                    ! polarization vectors 
c zet(mu) = ze(mu) - xxez*qe(mu) and zut(mu) = zu(mu) - xxmz*qu(mu)
            qq = qe(0)*qvv(0)-qe(1)*qvv(1)-qe(2)*qvv(2)-qe(3)*qvv(3)
            xxez = dotrc(qvv,ze)/qq
            xxea = dotrc(qvv,ae)/qq
            
	    qq = qu(0)*qvv(0)-qu(1)*qvv(1)-qu(2)*qvv(2)-qu(3)*qvv(3)
            xxuz = dotrc(qvv,zu)/qq            
            xxua = dotrc(qvv,au)/qq
            do mu = 0,3
               zet(mu+1) = ze(mu+1) - xxez*qe(mu)
               zut(mu+1) = zu(mu+1) - xxuz*qu(mu)
               
	       aet(mu+1) = ae(mu+1) - xxea*qe(mu)
               aut(mu+1) = au(mu+1) - xxua*qu(mu)
            enddo
            do mu = 5,6
               zet(mu) = ze(mu)
               zut(mu) = zu(mu)
               
               aet(mu) = ae(mu)
               aut(mu) = au(mu)
            enddo
         elseif (nlo.eq.1) then
            do mu = 1,6
               zet(mu) = ze(mu)
               zut(mu) = zu(mu)

               aet(mu) = ae(mu)
               aut(mu) = au(mu)
            enddo
         endif  ! end if nlo = -4 or +-5 for bos=721
         
c leptonic tensors
         j = 1

        if (with_spin2) then 

         call vvsp2tovv(3,q12(0,j),q34(0,j),qe,qu,ze,zu, 
     1           wwzz51(0,0,j),wwzz61(0,0,j),zzzz1(0,0,j),
     2           azzz1(0,0,j),zazz1(0,0,j),aazz1(0,0,j))

         call vvsp2tovv(4,q12(0,j),q34(0,j),qe,qu,ae,zu, 
     1           wwzz52(0,0,j),wwzz62(0,0,j),zzzz2(0,0,j),
     2           azzz2(0,0,j),zazz2(0,0,j),aazz2(0,0,j))

         do a=0,3
          do b=0,3
      sp2wwzz6(a,b,j)=wwzz61(a,b,j)+wwzz62(a,b,j)
      sp2wwzz5(a,b,j)=wwzz51(a,b,j)+wwzz52(a,b,j)
      sp2zzzz(a,b,j)=zzzz1(a,b,j)+zzzz2(a,b,j)
      sp2azzz(a,b,j)=azzz1(a,b,j)+azzz2(a,b,j)
      sp2zazz(a,b,j)=zazz1(a,b,j)+zazz2(a,b,j)
      sp2aazz(a,b,j)=aazz1(a,b,j)+aazz2(a,b,j)
          enddo
         enddo

        else  ! Higgs instead of spin-2 resonance      
       
         call vvtozzn_Hres(q12(0,j),q34(0,j),v,h,sp2aazz(0,0,j),sp2azzz(0,0,j),
     - 					    sp2zazz(0,0,j),sp2zzzz(0,0,j))
         call wwtozzn_Hres(q12(0,j),q34(0,j),v,h,sp2wwzz5(0,0,j)) ! q12 = W+
         call wwtozzn_Hres(q34(0,j),q12(0,j),v,h,sp2wwzz6(0,0,j)) ! q12 = W-	 

        endif


c for ZZ->2l 2v precalculate  leptonic tensors
      elseif (bos.eq.21) then
c reset lzs to .true. to force recalculation of virtual contributions in qqzzqq
         do i = 1,4
            lzs(i) = .true.
         enddo
c lepton spinors and Z polarization vectors ! 
c	lepton helicities not fixed 
c 		-> sum over all possible helicites in |M|**2
c
c select helicity: h ... random number for lepton helicity
c			 combination (h=1:4), but only h=1 or 3 
c	   	         contribute for ZZ -> 2l 2v (no righthanded neutrinos)

         if (ldoblha) then
           h = 2*blha_ranhelcomb-1
         endif
         h = mod(h+2,4)
         
         ie = sign(1,2-h)

c...Les Houches interface
         if ((lha.or.hepmc).and..not.doNLO) then
            helicity(1)= ie
            helicity(2)=-ie
            helicity(3)=-1
            helicity(4)= 1
         endif

		
         CALL IXXXXX(v(0,1),ZERO ,+ie,-1,lep) !e+   
         CALL OXXXXX(v(0,2),ZERO ,-ie, 1,lem) !e- 
         CALL OXXXXX(v(0,3),ZERO ,-1,+1,lum) !vm	  
         CALL IXXXXX(v(0,4),ZERO ,+1,-1,lup) !vm~	
         
         CALL JIOXXX(lep,lem,GZL ,ZMASS,ZWIDTH,ze) !Ze
         CALL JIOXXX(lep,lem,GAL ,ZERO ,ZERO  ,ae) !Ae
         CALL JIOXXX(lup,lum,GZN ,ZMASS,ZWIDTH,zu) !Zv
         
         do i = 1,4
            au(i) = 0d0
            zu(i) = -zu(i)
         enddo   			
         do i = 5,6    !add momentum info to A->vv current (au vanishes anyway) 
            au(i) = zu(i)
         enddo    	
	
         do mu = 0,3
            qe(mu) = v(mu,1)+v(mu,2)
            qu(mu) = v(mu,3)+v(mu,4)
            qvv(mu) = qe(mu) + qu(mu)
         enddo
         
         qe(4) = qe(0)**2-qe(1)**2-qe(2)**2-qe(3)**2
         qu(4) = qu(0)**2-qu(1)**2-qu(2)**2-qu(3)**2

         if (abs(nlo).eq.5 .or. nlo.eq.-4) then ! determine V/V 
                                                    ! polarization vectors 
c zet(mu) = ze(mu) - xxez*qe(mu) and zut(mu) = zu(mu) - xxmz*qu(mu)
            qq = qe(0)*qvv(0)-qe(1)*qvv(1)-qe(2)*qvv(2)-qe(3)*qvv(3)
            xxez = dotrc(qvv,ze)/qq
            xxea = dotrc(qvv,ae)/qq
            
	    qq = qu(0)*qvv(0)-qu(1)*qvv(1)-qu(2)*qvv(2)-qu(3)*qvv(3)
            xxuz = dotrc(qvv,zu)/qq            
            xxua = dotrc(qvv,au)/qq
            do mu = 0,3
               zet(mu+1) = ze(mu+1) - xxez*qe(mu)
               zut(mu+1) = zu(mu+1) - xxuz*qu(mu)
               
	       aet(mu+1) = ae(mu+1) - xxea*qe(mu)
               aut(mu+1) = au(mu+1) - xxua*qu(mu)
            enddo
            do mu = 5,6
               zet(mu) = ze(mu)
               zut(mu) = zu(mu)
               
               aet(mu) = ae(mu)
               aut(mu) = au(mu)
            enddo
         elseif (nlo.eq.1) then
            do mu = 1,6
               zet(mu) = ze(mu)
               zut(mu) = zu(mu)

               aet(mu) = ae(mu)
               aut(mu) = au(mu)
            enddo
         endif  ! end if nlo = -4 or +-5 for bos=21
         
c leptonic tensors
         j = 1
       
         if (with_kk) then      ! Kaluza-Klein scenario for bos=21
#ifdef WITH_KK
C for AA,AZ,ZA,ZZ -> e+ e- vm vm~ :
            call vvtozzn_kk(q12(0,j),q34(0,j),v,h,aazz(0,0,j),azzz(0,0,j),
     - 					    zazz(0,0,j),zzzz(0,0,j))
c
C for W+W-
            call wwtozzn_kk(q12(0,j),q34(0,j),v,h,wwzz5(0,0,j)) ! q12 = W+
            call wwtozzn_kk(q34(0,j),q12(0,j),v,h,wwzz6(0,0,j)) ! q12 = W-
	 
C for V1V2 --> e+ e- or  V1V2 -->  vm vm~ (Vi=A,Z)
c NCl tensors for NC process (k=1...4)

            call VVtoll_kk(2,1,h,q34(0,j),v,aaee(0,0,1,j),azee(0,0,1,j),
     -	 	 zaee(0,0,1,j),zzee(0,0,1,j))    	 !emit V on upper       
            call VVtoll_kk(1,1,h,q12(0,j),v,aaee(0,0,2,j),azee(0,0,2,j),
     -		 zaee(0,0,2,j),zzee(0,0,2,j))    	 !emit V on lower
            call VVtolln_kk(2,q34(0,j),v(0,3),aauu(0,0,1,j),azuu(0,0,1,j),
     -		 zauu(0,0,1,j),zzuu(0,0,1,j))    	 !emit V on upper       
            call VVtolln_kk(1,q12(0,j),v(0,3),aauu(0,0,2,j),azuu(0,0,2,j),
     -		zauu(0,0,2,j),zzuu(0,0,2,j))     	!emit V on lower


C for W+W- --> e+ e- or  W+W- -->  vm vm~
c  CCl tensors for CC (k=5)

            call WWtoll_kk(2,1,h,q34(0,j),v,CCee(0,0,1,j)) !emit V on upper
            call WWtoll_kk(1,1,h,q12(0,j),v,CCee(0,0,2,j)) !emit V on lower
            call WWtolln_kk(2,q34(0,j),v(0,3),CCuu(0,0,1,j)) !emit V on upper
            call WWtolln_kk(1,q12(0,j),v(0,3),CCuu(0,0,2,j)) !emit V on lower
	 
C for W-W+ --> e+ e- or vm vm~
c CCl for CC (k=6)
	
            call WWtoll_kk(1,1,h,q34(0,j),v,CCee6(0,0,1,j)) !emit V on upper
            call WWtoll_kk(2,1,h,q12(0,j),v,CCee6(0,0,2,j)) !emit V on lower
            call WWtolln_kk(1,q34(0,j),v(0,3),CCuu6(0,0,1,j)) !emit V on upper
            call WWtolln_kk(2,q12(0,j),v(0,3),CCuu6(0,0,2,j)) !emit V on lower
	 	 
c for V-> llll
	
            call vto4ln_kk(v,h,azz,zzztens)	
#endif	

         elseif (with_spin2) then     ! for bos = 21
C for AA,AZ,ZA,ZZ -> e+ e- vm vm~ :
            call vvtozzn_spin2(q12(0,j),q34(0,j),v,h,aazz(0,0,j),
     # 		azzz(0,0,j),zazz(0,0,j),zzzz(0,0,j))
c
C for W+W-
            call wwtozzn_spin2(q12(0,j),q34(0,j),v,h,wwzz5(0,0,j)) ! q12 = W+
            call wwtozzn_spin2(q34(0,j),q12(0,j),v,h,wwzz6(0,0,j)) ! q12 = W-

c other processes do not change
	 
C for V1V2 --> e+ e- or  V1V2 -->  vm vm~ (Vi=A,Z)
c NCl tensors for NC process (k=1...4)

            call VVtoll(2,1,h,q34(0,j),v,aaee(0,0,1,j),azee(0,0,1,j),
     -	 	 zaee(0,0,1,j),zzee(0,0,1,j))    	 !emit V on upper       
            call VVtoll(1,1,h,q12(0,j),v,aaee(0,0,2,j),azee(0,0,2,j),
     -		 zaee(0,0,2,j),zzee(0,0,2,j))    	 !emit V on lower
            call VVtolln(2,q34(0,j),v(0,3),aauu(0,0,1,j),azuu(0,0,1,j),
     -		 zauu(0,0,1,j),zzuu(0,0,1,j))    	 !emit V on upper       
            call VVtolln(1,q12(0,j),v(0,3),aauu(0,0,2,j),azuu(0,0,2,j),
     -		zauu(0,0,2,j),zzuu(0,0,2,j))     	!emit V on lower


C for W+W- --> e+ e- or  W+W- -->  vm vm~
c  CCl tensors for CC (k=5)

            call WWtoll(2,1,h,q34(0,j),v,CCee(0,0,1,j)) !emit V on upper
            call WWtoll(1,1,h,q12(0,j),v,CCee(0,0,2,j)) !emit V on lower
            call WWtolln(2,q34(0,j),v(0,3),CCuu(0,0,1,j)) !emit V on upper
            call WWtolln(1,q12(0,j),v(0,3),CCuu(0,0,2,j)) !emit V on lower
	 
C for W-W+ --> e+ e- or vm vm~
c CCl for CC (k=6)
	
            call WWtoll(1,1,h,q34(0,j),v,CCee6(0,0,1,j)) !emit V on upper
            call WWtoll(2,1,h,q12(0,j),v,CCee6(0,0,2,j)) !emit V on lower
            call WWtolln(1,q34(0,j),v(0,3),CCuu6(0,0,1,j)) !emit V on upper
            call WWtolln(2,q12(0,j),v(0,3),CCuu6(0,0,2,j)) !emit V on lower
	 	 
c for V-> llll
	
            call vto4ln(v,h,azz,zzztens)	


         elseif (with_anom) then ! anomalous gauge boson couplings
c        using global form factor for all tensors of one phase space point
c        this ensures proper cancellations for anomalous contributions
c        energy scale is invariant ZZ mass        
           call anomal_formfactor(q12(0,j),q34(0,j),qe(0),qu(0))

C for AA,AZ,ZA,ZZ -> e+ e- vm vm~ :
           call vvtozzn_anomal(q12(0,j),q34(0,j),v,h,aazz(0,0,j),azzz(0,0,j),
     # 					    zazz(0,0,j),zzzz(0,0,j))
c
C for W+W-
           call wwtozzn_anomal(q12(0,j),q34(0,j),v,h,wwzz5(0,0,j))  ! q12 = W+
           call wwtozzn_anomal(q34(0,j),q12(0,j),v,h,wwzz6(0,0,j))  ! q12 = W-
	 
C for V1V2 --> e+ e- or  V1V2 -->  vm vm~ (Vi=A,Z)
c NCl tensors for NC process (k=1...4)

           call VVtoll_anomal(2,1,h,q34(0,j),v,aaee(0,0,1,j),azee(0,0,1,j),
     #	 	 zaee(0,0,1,j),zzee(0,0,1,j))    	 !emit V on upper        
	   call VVtoll_anomal(1,1,h,q12(0,j),v,aaee(0,0,2,j),azee(0,0,2,j),
     #		 zaee(0,0,2,j),zzee(0,0,2,j))    	 !emit V on lower
           call VVtolln_anomal(2,q34(0,j),v(0,3),aauu(0,0,1,j),azuu(0,0,1,j),
     #		 zauu(0,0,1,j),zzuu(0,0,1,j))    	 !emit V on upper        
	   call VVtolln_anomal(1,q12(0,j),v(0,3),aauu(0,0,2,j),azuu(0,0,2,j),
     #		zauu(0,0,2,j),zzuu(0,0,2,j))     	!emit V on lower


C for W+W- --> e+ e- or  W+W- -->  vm vm~
c  CCl tensors for CC (k=5)

           call WWtoll_anomal(2,1,h,q34(0,j),v,CCee(0,0,1,j))     !emit V on upper
           call WWtoll_anomal(1,1,h,q12(0,j),v,CCee(0,0,2,j))     !emit V on lower
           call WWtolln_anomal(2,q34(0,j),v(0,3),CCuu(0,0,1,j))     !emit V on upper
           call WWtolln_anomal(1,q12(0,j),v(0,3),CCuu(0,0,2,j))     !emit V on lower
	 
C for W-W+ --> e+ e- or vm vm~
c CCl for CC (k=6)
	
           call WWtoll_anomal(1,1,h,q34(0,j),v,CCee6(0,0,1,j))     !emit V on upper
           call WWtoll_anomal(2,1,h,q12(0,j),v,CCee6(0,0,2,j))     !emit V on lower
           call WWtolln_anomal(1,q34(0,j),v(0,3),CCuu6(0,0,1,j))     !emit V on upper
           call WWtolln_anomal(2,q12(0,j),v(0,3),CCuu6(0,0,2,j))     !emit V on lower
	 	 
c for V-> llll
	
	  call vto4ln_anomal(v,h,azz,zzztens)


c --------------------------------------------
         else                   !SM  for bos = 21
C for AA,AZ,ZA,ZZ -> e+ e- vm vm~ :
            call vvtozzn(q12(0,j),q34(0,j),v,h,aazz(0,0,j),azzz(0,0,j),
     - 					    zazz(0,0,j),zzzz(0,0,j))
c
C     for W+W-
            call wwtozzn(q12(0,j),q34(0,j),v,h,wwzz5(0,0,j)) ! q12 = W+
            call wwtozzn(q34(0,j),q12(0,j),v,h,wwzz6(0,0,j)) ! q12 = W-
	 
C for V1V2 --> e+ e- or  V1V2 -->  vm vm~ (Vi=A,Z)
c NCl tensors for NC process (k=1...4)

            call VVtoll(2,1,h,q34(0,j),v,aaee(0,0,1,j),azee(0,0,1,j),
     -	 	 zaee(0,0,1,j),zzee(0,0,1,j))    	 !emit V on upper       
            call VVtoll(1,1,h,q12(0,j),v,aaee(0,0,2,j),azee(0,0,2,j),
     -		 zaee(0,0,2,j),zzee(0,0,2,j))    	 !emit V on lower
            call VVtolln(2,q34(0,j),v(0,3),aauu(0,0,1,j),azuu(0,0,1,j),
     -		 zauu(0,0,1,j),zzuu(0,0,1,j))    	 !emit V on upper       
            call VVtolln(1,q12(0,j),v(0,3),aauu(0,0,2,j),azuu(0,0,2,j),
     -		zauu(0,0,2,j),zzuu(0,0,2,j))     	!emit V on lower


C for W+W- --> e+ e- or  W+W- -->  vm vm~
c  CCl tensors for CC (k=5)

            call WWtoll(2,1,h,q34(0,j),v,CCee(0,0,1,j)) !emit V on upper
            call WWtoll(1,1,h,q12(0,j),v,CCee(0,0,2,j)) !emit V on lower
            call WWtolln(2,q34(0,j),v(0,3),CCuu(0,0,1,j)) !emit V on upper
            call WWtolln(1,q12(0,j),v(0,3),CCuu(0,0,2,j)) !emit V on lower
	 
C for W-W+ --> e+ e- or vm vm~
c CCl for CC (k=6)
	
            call WWtoll(1,1,h,q34(0,j),v,CCee6(0,0,1,j)) !emit V on upper
            call WWtoll(2,1,h,q12(0,j),v,CCee6(0,0,2,j)) !emit V on lower
            call WWtolln(1,q34(0,j),v(0,3),CCuu6(0,0,1,j)) !emit V on upper
            call WWtolln(2,q12(0,j),v(0,3),CCuu6(0,0,2,j)) !emit V on lower
            
c for V-> llll
	
            call vto4ln(v,h,azz,zzztens)	
		
         endif   !  end of kk/spin2/sm for bos=21

c --------------------------------------------

      elseif (bos.eq.211) then   !ZA -> l+l- A
cjp     j?????
         j=1 
c     reset lzs to .true. to force recalculation of virtual contributions in qqwwqq
         do i = 1,4
            lzs(i) = .true.
         enddo


c     select helicity: h ... random number for lepton helicity
c     combination (h=1:2) 
      
      if (ldoblha) then
        h = blha_ranhelcomb
      endif
      if(photon_hel.eq.0) then
         h = mod(h,4) + 1
         ie = sign(1,2-h)
         iu = (-1)**(h+1)
      else
         h = mod(h,2) + 1
         ie = (-1)**h
         iu = photon_hel
      endif
         
cfc         print*, "h", h
cfc         print*, "ie_e+", ie
cfc         print*, "ie_e-", -ie
cfc         print*, "iu", iu

c...  Les Houches interface         
         if ((lha.or.hepmc).and..not.doNLO) then
            helicity(1)= ie
            helicity(2)=-ie 
            helicity(3)= iu
         endif

         CALL IXXXXX(v(0,1),ZERO ,+ie,-1,lep) !e+   
         CALL OXXXXX(v(0,2),ZERO ,-ie, 1,lem) !e- 
         CALL VXXXXX(v(0,3),ZERO ,+iu, 1,ea) !a


         CALL JIOXXX(lep,lem,GZL ,ZMASS,ZWIDTH,ze) !Zl
         CALL JIOXXX(lep,lem,GAL ,ZERO,ZERO,ae) !al

cfc to recover old notation
         do i = 1,4
            zu(i) = 0d0
            au(i) = ea(i)
         enddo  
         do i = 5,6    !add momentum info to Z->vv current (zu vanishes anyway) 
            au(i) = ea(i)
            zu(i)=  0d0 !ea(i)
         enddo    
         
         do mu = 0,3
            qe(mu) = v(mu,1)+v(mu,2)
            qa(mu) = v(mu,3)
            qvv(mu) = qe(mu) + qa(mu)
         enddo
         
         qe(4) = qe(0)**2-qe(1)**2-qe(2)**2-qe(3)**2
         qa(4) = qa(0)**2-qa(1)**2-qa(2)**2-qa(3)**2

cjp to recover zz notation
         do mu=0,4
            qu(mu)=qa(mu)
         enddo
         
c leptonic tensors
c we keep the name of the output variables.
c to reuse most of the original qqzaqq.f 

	if(with_anom) then! anomalous gauge boson couplings
c        using global form factor for all tensors of one phase space point
c        this ensures proper cancellations for anomalous contributions
c        energy scale is invariant ZA mass        
           call anomal_formfactor(q12(0,j),q34(0,j),qe(0),qu(0))

            call vvtoza_anomal(q12(0,j),q34(0,j),v,h,aazz(0,0,j),azzz(0,0,j),
     - zazz(0,0,j),zzzz(0,0,j))

C for W+W-
            call wwtoza_anomal(q12(0,j),q34(0,j),v,h,wwzz5(0,0,j)) ! q12 = W+
            call wwtoza_anomal(q34(0,j),q12(0,j),v,h,wwzz6(0,0,j)) ! q12 = W-

C for V1V2 --> e+ e- 
c NCl tensors for NC process (k=1...4)

            call VVtoll_anomal(2,1,h,q34(0,j),v,aaee(0,0,1,j),azee(0,0,1,j),
     - zaee(0,0,1,j),zzee(0,0,1,j))  !emit V on upper       
            call VVtoll_anomal(1,1,h,q12(0,j),v,aaee(0,0,2,j),azee(0,0,2,j),
     - zaee(0,0,2,j),zzee(0,0,2,j))  !emit V on lower
	    aauu=0d0
	    azuu=0d0
	    zauu=0d0
	    zzuu=0d0
	    

C for W+W- --> e+ e-
c  CCl tensors for CC (k=5)

            call WWtoll_anomal(2,1,h,q34(0,j),v,CCee(0,0,1,j)) !emit V on upper
            call WWtoll_anomal(1,1,h,q12(0,j),v,CCee(0,0,2,j)) !emit V on lower
C     for W+W- -->  A

cfc ???? k=5 => ud -> sc entiendo es 1???? 
cfc to check first argument

            call calc_wwtoa_anomal(2,iu,q34(0,j),v,CCuu(0,0,1,j)) !emit Z -> e+ e- on upper
            call calc_wwtoa_anomal(1,iu,q12(0,j),v,CCuu(0,0,2,j)) !emit Z -> e+ e- on lower
 
C for W-W+ --> e+ e- 
c CCl for CC (k=6)

            call WWtoll_anomal(1,1,h,q34(0,j),v,CCee6(0,0,1,j)) !emit A on upper
            call WWtoll_anomal(2,1,h,q12(0,j),v,CCee6(0,0,2,j)) !emit A on lower

C     for W-W+ -->  A

cfc  k=6 => du -> cs entiendo es 2
cfc to check first argument

            call calc_wwtoa_anomal(1,iu,q34(0,j),v,CCuu6(0,0,1,j)) !emit Z -> e+ e- on upper
            call calc_wwtoa_anomal(2,iu,q12(0,j),v,CCuu6(0,0,2,j)) !emit Z -> e+ e- on lower

        else ! SM

            call vvtoza(q12(0,j),q34(0,j),v,h,aazz(0,0,j),azzz(0,0,j),
     - zazz(0,0,j),zzzz(0,0,j))

C for W+W-
            call wwtoza(q12(0,j),q34(0,j),v,h,wwzz5(0,0,j)) ! q12 = W+
            call wwtoza(q34(0,j),q12(0,j),v,h,wwzz6(0,0,j)) ! q12 = W-

C for V1V2 --> e+ e- 
c NCl tensors for NC process (k=1...4)

            call VVtoll(2,1,h,q34(0,j),v,aaee(0,0,1,j),azee(0,0,1,j),
     - zaee(0,0,1,j),zzee(0,0,1,j))  !emit V on upper       
            call VVtoll(1,1,h,q12(0,j),v,aaee(0,0,2,j),azee(0,0,2,j),
     - zaee(0,0,2,j),zzee(0,0,2,j))  !emit V on lower
	    aauu=0d0
	    azuu=0d0
	    zauu=0d0
	    zzuu=0d0
	    

C for W+W- --> e+ e-
c  CCl tensors for CC (k=5)

            call WWtoll(2,1,h,q34(0,j),v,CCee(0,0,1,j)) !emit V on upper
            call WWtoll(1,1,h,q12(0,j),v,CCee(0,0,2,j)) !emit V on lower
C     for W+W- -->  A

cfc ???? k=5 => ud -> sc entiendo es 1???? 
cfc to check first argument

            call calc_wwtoa(2,iu,q34(0,j),v,CCuu(0,0,1,j)) !emit Z -> e+ e- on upper
            call calc_wwtoa(1,iu,q12(0,j),v,CCuu(0,0,2,j)) !emit Z -> e+ e- on lower
 
C for W-W+ --> e+ e- 
c CCl for CC (k=6)

            call WWtoll(1,1,h,q34(0,j),v,CCee6(0,0,1,j)) !emit A on upper
            call WWtoll(2,1,h,q12(0,j),v,CCee6(0,0,2,j)) !emit A on lower

C     for W-W+ -->  A

cfc  k=6 => du -> cs entiendo es 2
cfc to check first argument

            call calc_wwtoa(1,iu,q34(0,j),v,CCuu6(0,0,1,j)) !emit Z -> e+ e- on upper
            call calc_wwtoa(2,iu,q12(0,j),v,CCuu6(0,0,2,j)) !emit Z -> e+ e- on lower
	endif

C     for V --> e+ e- a ( doesn't change with AC)
         call vtolla(v,h,azz,zzztens)

c----------------------------------------

      elseif (bos.eq.212) then   !ZA -> vv~ A
c     reset lzs to .true. to force recalculation of virtual contributions in qqwwqq
      j=1
         do i = 1,4
            lzs(i) = .true.
         enddo
c     select helicity: h ... random number for lepton helicity
c     combination (h=1:2) 

         if (ldoblha) then
           h = blha_ranhelcomb
         endif
         h = mod(h,2) + 1
         ie = 1
         iu = sign(1,1-h)
     

c...  Les Houches interface         
         if ((lha.or.hepmc).and..not.doNLO) then
            helicity(1)= -1
            helicity(2)= +1 
            helicity(3)= iu
         endif

         CALL OXXXXX(v(0,1),ZERO ,-1,+1,lum) !vm	  
         CALL IXXXXX(v(0,2),ZERO ,+1,-1,lup) !vm~	
         CALL VXXXXX(v(0,3),ZERO ,+iu, 1,ea) !a
         
         au=ea
         zu=0d0

         CALL JIOXXX(lup,lum,GZN,ZMASS, ZWIDTH,ze) !Zv
         ze(1:4)=-ze(1:4) ! different fermion order -> sign change in mg amplitudes
         ae(1:4)=0d0
         ae(5:6)=ze(5:6)
         
         
         do mu = 0,3
            qe(mu) = v(mu,1)+v(mu,2)
            qa(mu) = v(mu,3)
            qvv(mu) = qe(mu) + qa(mu)
         enddo
         
         qe(4) = qe(0)**2-qe(1)**2-qe(2)**2-qe(3)**2
         qa(4) = qa(0)**2-qa(1)**2-qa(2)**2-qa(3)**2
         qu=qa

     
c leptonic tensors
      if(with_anom) then! anomalous gauge boson couplings
c        using global form factor for all tensors of one phase space point
c        this ensures proper cancellations for anomalous contributions
c        energy scale is invariant ZA mass        
      
           call anomal_formfactor(q12(0,j),q34(0,j),qe(0),qu(0))

c vv->vv (not in SM)
          call vvtozan_anomal(q12(0,j),q34(0,j),v,h,aazz(0,0,j),azzz(0,0,j),
     - zazz(0,0,j),zzzz(0,0,j))
     
C for W+W- --> ve ve~ a
            call wwtozan_anomal(q12(0,j),q34(0,j),v,h,wwzz5(0,0,j)) ! q12 = W+
            call wwtozan_anomal(q34(0,j),q12(0,j),v,h,wwzz6(0,0,j)) ! q12 = W-

C V1V2 --> ve ve~  
c NCl tensors for NC process (k=1...4)
      qu=qe ! VVtolln, WWtolln assume Z_l Z_n
            call VVtolln_anomal(2,q34(0,j),v(0,1),aaee(0,0,1,j),azee(0,0,1,j),
     - zaee(0,0,1,j),zzee(0,0,1,j))     !emit V on upper       
            call VVtolln_anomal(1,q12(0,j),v(0,1),aaee(0,0,2,j),azee(0,0,2,j),
     - zaee(0,0,2,j),zzee(0,0,2,j))     !emit V on lower

C for W+W- --> ve ve~
c  CCl tensors for CC (k=5)

            call WWtolln_anomal(2,q34(0,j),v(0,1),CCee(0,0,1,j)) !emit V on upper
            call WWtolln_anomal(1,q12(0,j),v(0,1),CCee(0,0,2,j)) !emit V on lower
C     for W+W- -->  A

cfc  k=5 => ud -> sc entiendo es 1

            call calc_wwtoa_anomal(2,iu,q34(0,j),v,CCuu(0,0,1,j)) !emit Z -> ve ve~ on upper
            call calc_wwtoa_anomal(1,iu,q12(0,j),v,CCuu(0,0,2,j)) !emit Z -> ve ve~ on lower

C for W-W+ --> ve ve~ 
c CCl for CC (k=6)

            call WWtolln_anomal(1,q34(0,j),v,CCee6(0,0,1,j)) !emit A on upper
            call WWtolln_anomal(2,q12(0,j),v,CCee6(0,0,2,j)) !emit A on lower

C     for W-W+ -->  A

cfc  k=6 => du -> cs entiendo es 2

            call calc_wwtoa_anomal(1,iu,q34(0,j),v,CCuu6(0,0,1,j)) !emit Z -> ve ve~ on upper
            call calc_wwtoa_anomal(2,iu,q12(0,j),v,CCuu6(0,0,2,j)) !emit Z -> ve ve~ on lower
            
         else ! SM

C for W+W- --> ve ve~ a
            call wwtozan(q12(0,j),q34(0,j),v,h,wwzz5(0,0,j)) ! q12 = W+
            call wwtozan(q34(0,j),q12(0,j),v,h,wwzz6(0,0,j)) ! q12 = W-

C V1V2 --> ve ve~  
c NCl tensors for NC process (k=1...4)
      qu=qe ! VVtolln, WWtolln assume Z_l Z_n
            call VVtolln(2,q34(0,j),v(0,1),aaee(0,0,1,j),azee(0,0,1,j),
     - zaee(0,0,1,j),zzee(0,0,1,j))     !emit V on upper       
            call VVtolln(1,q12(0,j),v(0,1),aaee(0,0,2,j),azee(0,0,2,j),
     - zaee(0,0,2,j),zzee(0,0,2,j))     !emit V on lower

C for W+W- --> ve ve~
c  CCl tensors for CC (k=5)

            call WWtolln(2,q34(0,j),v(0,1),CCee(0,0,1,j)) !emit V on upper
            call WWtolln(1,q12(0,j),v(0,1),CCee(0,0,2,j)) !emit V on lower
C     for W+W- -->  A

cfc  k=5 => ud -> sc entiendo es 1

            call calc_wwtoa(2,iu,q34(0,j),v,CCuu(0,0,1,j)) !emit Z -> ve ve~ on upper
            call calc_wwtoa(1,iu,q12(0,j),v,CCuu(0,0,2,j)) !emit Z -> ve ve~ on lower

C for W-W+ --> ve ve~ 
c CCl for CC (k=6)

            call WWtolln(1,q34(0,j),v,CCee6(0,0,1,j)) !emit A on upper
            call WWtolln(2,q12(0,j),v,CCee6(0,0,2,j)) !emit A on lower

C     for W-W+ -->  A

cfc  k=6 => du -> cs entiendo es 2

            call calc_wwtoa(1,iu,q34(0,j),v,CCuu6(0,0,1,j)) !emit Z -> ve ve~ on upper
            call calc_wwtoa(2,iu,q12(0,j),v,CCuu6(0,0,2,j)) !emit Z -> ve ve~ on lower
     
      endif
      azz=0d0
      zzztens=0d0
         qu=qa
c-----------------------------

      elseif(bos.eq.16.or.bos.eq.61) then    ! bos=61=HAjj
         CALL Get_HAjj_currents(p,v,nlo)

#ifdef WITH_NLO 
      elseif (nlo.eq.1 .and. bos.eq.2) then
         call BCD_fill(p(0,1),p(0,2), !beam
     1        p(0,3),p(0,4),    !jet
     2        v(0,1),v(0,2))    !decay     
#endif

c -----------------------------------------------------------------------
c
      elseif(bos.eq.66) then ! bos=66=HHjj
c for HH precalculate  leptonic tensors
c reset lzs to .true. to force recalculation of virtual contributions in qqHHqq
         do i = 1,4
            lzs(i) = .true.
         enddo

c...Les Houches interface
         if ((lha.or.hepmc).and..not.doNLO) then
            helicity(1)= 0
            helicity(2)=0
         endif

         if(procid.eq.HHjj) then
            CALL SXXXXX(v(0,1   ),1,wp)                         
            CALL SXXXXX(v(0,2   ),1,wm)
         else ! including H decays
            do mu = 0,3
               vtemp(mu,1) = v(mu,1)+v(mu,2)
               vtemp(mu,2) = v(mu,3)+v(mu,4)
            enddo
            CALL SXXXXX(vtemp(0,1   ),1,wp)                         
            CALL SXXXXX(vtemp(0,2   ),1,wm)
         endif
         
         do mu = 0,3
            if(procid.eq.HHjj) then
               qe(mu) = v(mu,1)
               qu(mu) = v(mu,2)
            else
               qe(mu) = v(mu,1)+v(mu,2)
               qu(mu) = v(mu,3)+v(mu,4)
            endif
            qvv(mu) = qe(mu) + qu(mu)
         enddo
	 	 
         qe(4) = qe(0)**2-qe(1)**2-qe(2)**2-qe(3)**2
         qu(4) = qu(0)**2-qu(1)**2-qu(2)**2-qu(3)**2

c leptonic tensors
         j = 1

c-----------------------------
c
C for W+W-
         call wwtohh(q12(0,j),q34(0,j),v,wwhh5(0,0,j)) ! q12 = W+
         call wwtohh(q34(0,j),q12(0,j),v,wwhh6(0,0,j)) ! q12 = W-
C for ZZ
         call zztohh(q12(0,j),q34(0,j),v,zzhh(0,0,j)) ! q12 = Z

c -----------------------------------------------------------------------
c
      endif

      if(higgsscheme.eq.2) then    ! higgsscheme=2=CPS
        SELECT CASE(procid)
        CASE(Hjj)
          do mu=0,3
            qvv(mu) = v(mu,1)+v(mu,2)+v(mu,3)+v(mu,4)
          enddo
          qinv = sqrt(qvv(0)**2-qvv(1)**2-qvv(2)**2-qvv(3)**2)
          call VBFNLO_HIGGSTO_gridHt(qinv,qgammaq)
          qgammaq = qinv*qgammaq
        CASE DEFAULT
          qgammaq = xmg(6)
        END SELECT
      else
        qgammaq = xmg(6)
      endif
c
c
c ------------------------------------------------------------------------
c
c if desired, print output information
c
       if ( ldebug ) then

 10      format( " p(", i1, ") = ", 4(f10.3, 2x) )
 20      format( " v(", i1, ") = ", 4(f10.3, 2x) )
c
         write(6,*) " "
         print*," vector boson decay products momenta "
         do i = 1, n_v
            write(6,20) i, v(0,i), v(1,i), v(2,i), v(3,i)
         end do
         print*," parton momenta "
         do i = 1, n_p
            write(6,10) i, p(0,i), p(1,i), p(2,i), p(3,i)
         end do
c
         write(6,*) " "
         write(6,*) "xi(1), xi(2) =", xi
         write(6,*) "mu_f =",sqrt(mufsq(1,1)),sqrt(mufsq(2,1))
       end if   ! end ldebug

ccccc no pdf-information required

       if (ldoblha) then

         physToDiag(1)=blha_physdiagmap(1,blha_cursubproc,blha_curproc)
         physToDiag(2)=blha_physdiagmap(2,blha_cursubproc,blha_curproc)
         physToDiag(3)=blha_physdiagmap(3,blha_cursubproc,blha_curproc)
         physToDiag(4)=blha_physdiagmap(4,blha_cursubproc,blha_curproc)
         fsign(1) = blha_fsign(1,blha_cursubproc,blha_curproc)
         fsign(2) = blha_fsign(2,blha_cursubproc,blha_curproc)
         fsign(3) = blha_fsign(3,blha_cursubproc,blha_curproc)
         fsign(4) = blha_fsign(4,blha_cursubproc,blha_curproc)

         call wbf_zh(xi,p,v,rn,physToDiag,fsign,bos,nlo,ps_number,
     1            pdf,Qpdf,res,nmin,nmax)

         return
       endif

c
c call PDF subroutine in order to determine parton
c distributions in the incoming (anti)protons.
c
      q_sf = sqrt(mufsq(1,1))
      call pdfproton( xi(1), q_sf, pdf(-6,1) )        ! f_a(x1)
      Qpdf(-7,1,1) = GAMpdf
      Qpdf(7,1,1) = GAMpdf

      q_sf = sqrt(mufsq(2,1))
      call pdfproton( xi(2), q_sf, pdf(-6,2) )        ! f_a(x2)     
      Qpdf(-7,2,1) = GAMpdf
      Qpdf(7,2,1) = GAMpdf


*** start sophy qed : setting PDFs.  NB: send Qpdf through wbf_zh
      q_sf = sqrt(mufsq(1,1))
      zz_qed = (1.0d0 - xi(1)) * rn**2 + xi(1) !for z integration
      x1 = xi(1)/zz_qed
      call pdfproton(x1, q_sf, Qpdf(-6,1,2) ) ! f_a(x1/z)
      Qpdf(-7,1,2) = GAMpdf
      Qpdf(7,1,2) = GAMpdf

      q_sf = sqrt(mufsq(2,1))
      zz_qed = (1.0d0 - xi(2)) * rn**2 + xi(2) !for z integration
      x2 = xi(2)/zz_qed
      call pdfproton(x2, q_sf, Qpdf(-6,2,2) ) ! f_a(x2/z)
      Qpdf(-7,2,2) = GAMpdf
      Qpdf(7,2,2) = GAMpdf

      do i = 1, 2
         do mu = -6, 6
            Qpdf(mu,i,1) = pdf(mu,i)
         end do
      end do
*** end sophy qed

c select helicity
      if (bos.eq.2) then
         jsig = min(8*rn+1d0,8.01d0)
c         jsig = 0
      else
         jsig = 0
      endif

!       pdf(-6:0,1)=0d0
!       pdf(2:6,1)=0d0
!       pdf(-6:3,2)=0d0
!       pdf(5:6,2)=0d0
c
c      pdf(2,1) = 0d0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     Now sum over the subprocesses contributing to Z production

      nmax = FL_ZH(fsign,-1)            !reset counter for subprocesses to 0

C*******************  q1 q3 ---> q2 q4  V V  **********************
         
c   physToDiag(ext.momentum label) = Feynman diagram label

      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=3
      physToDiag(3)=2    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=4
C NOTE: for call of wbf_zh it is important that p(*,1,*) and p(*,3,*)
c correspond to 1-2 fermion line ALWAYS, i.e physToDiag(1/2)={1,3} and 
c similarly physToDiag(3/4)={2,4} for the 3-4 fermion line
      fsign(1) = 1
      fsign(2) = 1
      fsign(3) = 1
      fsign(4) = 1

      call wbf_zh(xi,p,v,rn,physToDiag,fsign,bos,nlo,ps_number,
     1            pdf,Qpdf,res,nmin,nmax)

*      write(*,*)'1 res(1-4) =', res(1),res(2),res(3),res(4)
     

c      if (ldebug) call resprint2(nmin,nmax,res)
      if (init.eq.1) nproc(1) = nmax

C*******************  q1 qb4 ---> q2 qb3 V V  **********************

      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=4
      physToDiag(3)=2    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=3

      fsign(3) = -1
      fsign(4) = -1

      call wbf_zh(xi,p,v,rn,physToDiag,fsign,bos,nlo,ps_number,
     1         pdf,Qpdf,res,nmin,nmax)

*      write(*,*)'2 res(1-4) =', res(1),res(2),res(3),res(4)

c      if (ldebug) call resprint2(nmin,nmax,res)
      if (init.eq.1) nproc(2) = nmax

C*******************  qbar2 q3 ---> qbar1 q4 V V   **********************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=3
      physToDiag(3)=1    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=4

      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = 1
      fsign(4) = 1
      
      call wbf_zh(xi,p,v,rn,physToDiag,fsign,bos,nlo,ps_number,
     1         pdf,Qpdf,res,nmin,nmax)

*      write(*,*)'3 res(1-4) =', res(1),res(2),res(3),res(4)

c      if (ldebug) call resprint2(nmin,nmax,res)
      if (init.eq.1) nproc(3) = nmax

C*******************  qbar2 qb4 ---> qbar1 qb3 V V  ******************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=4
      physToDiag(3)=1    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=3

      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = -1 
      fsign(4) = -1
         
      call wbf_zh(xi,p,v,rn,physToDiag,fsign,bos,nlo,ps_number,
     1         pdf,Qpdf,res,nmin,nmax)

*      write(*,*)'4 res(1-4) =', res(1),res(2),res(3),res(4)

c      if (ldebug) call resprint2(nmin,nmax,res)
      if (init.eq.1) nproc(4) = nmax


c**************  end of process evaluation part one ******************


! now run again for down-type decay products if 93 93 / 94 94 is 
! requested and the z decays hadronically.


      SELECT CASE(procid)
      CASE(ZZhadjj)
       if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
        ! set couplings for Z hadronic decay into ddbar
        call setZtodd

        ! reset lzs to .true. to force recalculation of virtual contributions in qqzzqq
        do i = 1,4
           lzs(i) = .true.
        enddo

        CALL JIOXXX(lep,lem,GZ_ZF ,ZMASS,ZWIDTH,ze)    !Zl
        CALL JIOXXX(lep,lem,GZ_AF ,ZERO ,ZERO  ,ae)    !Al

         do mu = 0,3
            qe(mu) = v(mu,1)+v(mu,2)
            qu(mu) = v(mu,3)+v(mu,4)
            qvv(mu) = qe(mu) + qu(mu)
         enddo

         qe(4) = qe(0)**2-qe(1)**2-qe(2)**2-qe(3)**2
         qu(4) = qu(0)**2-qu(1)**2-qu(2)**2-qu(3)**2

         if (abs(nlo).eq.5 .or. nlo.eq.-4) then  ! determine V/V 
                                                 ! polarization vectors 
            ! zet(mu) = ze(mu) - xxez*qe(mu) and zut(mu) = zu(mu) - xxmz*qu(mu)
            qq = qe(0)*qvv(0)-qe(1)*qvv(1)-qe(2)*qvv(2)-qe(3)*qvv(3)
            xxez = dotrc(qvv,ze)/qq
            xxea = dotrc(qvv,ae)/qq
            
            qq = qu(0)*qvv(0)-qu(1)*qvv(1)-qu(2)*qvv(2)-qu(3)*qvv(3)
            xxuz = dotrc(qvv,zu)/qq            
            xxua = dotrc(qvv,au)/qq
            do mu = 0,3
               zet(mu+1) = ze(mu+1) - xxez*qe(mu)
               zut(mu+1) = zu(mu+1) - xxuz*qu(mu)
               
               aet(mu+1) = ae(mu+1) - xxea*qe(mu)
               aut(mu+1) = au(mu+1) - xxua*qu(mu)
            enddo
            do mu = 5,6
               zet(mu) = ze(mu)
               zut(mu) = zu(mu)

               aet(mu) = ae(mu)
               aut(mu) = au(mu)
            enddo
         elseif (nlo.eq.1) then
            do mu = 1,6
               zet(mu) = ze(mu)
               zut(mu) = zu(mu)

               aet(mu) = ae(mu)
               aut(mu) = au(mu)
            enddo
         endif  ! end of nlo = -4 or +-5

         j = 1

         ! for AA,AZ,ZA,ZZ -> e+ e- mu+ mu- :
         call vvtozz_had(q12(0,j),q34(0,j),v,h,aazz(0,0,j),azzz(0,0,j),zazz(0,0,j),zzzz(0,0,j))
         ! for W+W-
         call wwtozz_had(q12(0,j),q34(0,j),v,h,wwzz5(0,0,j))  ! q12 = W+
         call wwtozz_had(q34(0,j),q12(0,j),v,h,wwzz6(0,0,j))  ! q12 = W-
         ! for V1V2 --> e+ e- or  V1V2 -->  mu+ mu- (Vi=A,Z)
         ! NCl tensors for NC process (k=1...4)
         call VVtoll_had(2,1,h,q34(0,j),v,aaee(0,0,1,j),azee(0,0,1,j),zaee(0,0,1,j),zzee(0,0,1,j))     !emit V on upper    
         call VVtoll_had(1,1,h,q12(0,j),v,aaee(0,0,2,j),azee(0,0,2,j),zaee(0,0,2,j),zzee(0,0,2,j))     !emit V on lower
         call VVtoll_had(2,2,h,q34(0,j),v,aauu(0,0,1,j),azuu(0,0,1,j),zauu(0,0,1,j),zzuu(0,0,1,j))     !emit V on upper     
         call VVtoll_had(1,2,h,q12(0,j),v,aauu(0,0,2,j),azuu(0,0,2,j),zauu(0,0,2,j),zzuu(0,0,2,j))     !emit V on lower
         ! for W+W- --> e+ e- or  W+W- -->  mu+ mu-
         ! CCl tensors for CC (k=5)
         call WWtoll_had(2,1,h,q34(0,j),v,CCee(0,0,1,j))     !emit V on upper
         call WWtoll_had(1,1,h,q12(0,j),v,CCee(0,0,2,j))     !emit V on lower
         call WWtoll_had(2,2,h,q34(0,j),v,CCuu(0,0,1,j))     !emit V on upper
         call WWtoll_had(1,2,h,q12(0,j),v,CCuu(0,0,2,j))     !emit V on lower
         ! for W-W+ --> e+ e- or mu+ mu-
         ! CCl for CC (k=6)
         call WWtoll_had(1,1,h,q34(0,j),v,CCee6(0,0,1,j))     !emit V on upper
         call WWtoll_had(2,1,h,q12(0,j),v,CCee6(0,0,2,j))     !emit V on lower
         call WWtoll_had(1,2,h,q34(0,j),v,CCuu6(0,0,1,j))     !emit V on upper
         call WWtoll_had(2,2,h,q12(0,j),v,CCuu6(0,0,2,j))     !emit V on lower
         ! for V-> llll
         call vto4l_had(v,h,azz,zzztens)


C*******************  q1 q3 ---> q2 q4  V V  **********************
         
c   physToDiag(ext.momentum label) = Feynman diagram label

      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=3
      physToDiag(3)=2    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=4
C NOTE: for call of wbf_zh it is important that p(*,1,*) and p(*,3,*)
c correspond to 1-2 fermion line ALWAYS, i.e physToDiag(1/2)={1,3} and 
c similarly physToDiag(3/4)={2,4} for the 3-4 fermion line
      fsign(1) = 1
      fsign(2) = 1
      fsign(3) = 1
      fsign(4) = 1

      call wbf_zh(xi,p,v,rn,physToDiag,fsign,bos,nlo,ps_number,
     1            pdf,Qpdf,res,nmin,nmax)

*      write(*,*)'1 res(1-4) =', res(1),res(2),res(3),res(4)
     

c      if (ldebug) call resprint2(nmin,nmax,res)
      if (init.eq.1) nproc(1) = nmax

C*******************  q1 qb4 ---> q2 qb3 V V  **********************

      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=4
      physToDiag(3)=2    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=3

      fsign(3) = -1
      fsign(4) = -1

      call wbf_zh(xi,p,v,rn,physToDiag,fsign,bos,nlo,ps_number,
     1         pdf,Qpdf,res,nmin,nmax)

*      write(*,*)'2 res(1-4) =', res(1),res(2),res(3),res(4)

c      if (ldebug) call resprint2(nmin,nmax,res)
      if (init.eq.1) nproc(2) = nmax

C*******************  qbar2 q3 ---> qbar1 q4 V V   **********************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=3
      physToDiag(3)=1    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=4

      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = 1
      fsign(4) = 1
      
      call wbf_zh(xi,p,v,rn,physToDiag,fsign,bos,nlo,ps_number,
     1         pdf,Qpdf,res,nmin,nmax)

*      write(*,*)'3 res(1-4) =', res(1),res(2),res(3),res(4)

c      if (ldebug) call resprint2(nmin,nmax,res)
      if (init.eq.1) nproc(3) = nmax

C*******************  qbar2 qb4 ---> qbar1 qb3 V V  ******************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=4
      physToDiag(3)=1    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=3

      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = -1 
      fsign(4) = -1
         
      call wbf_zh(xi,p,v,rn,physToDiag,fsign,bos,nlo,ps_number,
     1         pdf,Qpdf,res,nmin,nmax)

*      write(*,*)'4 res(1-4) =', res(1),res(2),res(3),res(4)

c      if (ldebug) call resprint2(nmin,nmax,res)
      if (init.eq.1) nproc(4) = nmax


c***  end of extra process evaluation for hadronic decay  ***

       endif
      end select


c*****************  end of process evaluation  *********************


      SELECT CASE(procid)
      CASE(WPhadWMjj)
        call hadDecayFactor_W(finalquarks(1), abs(mjj2(v(0,1),v(0,2))), N_gen_W, fac_W)
        DO IRES = 1,NMAX
          res(IRES) = res(IRES) * fac_W            ! factor fac_W for hadronic decay
        ENDDO
      CASE(WPWMhadjj)
        call hadDecayFactor_W(finalquarks(1), abs(mjj2(v(0,3),v(0,4))), N_gen_W, fac_W)
        DO IRES = 1,NMAX
          res(IRES) = res(IRES) * fac_W            ! factor fac_W for hadronic decay
        ENDDO
      CASE(ZZhadjj)
        call hadDecayFactor_Z(finalquarks(1), abs(mjj2(v(0,1),v(0,2))), N_gen_up, N_gen_down, fac_Z_up, fac_Z_down)
        DO IRES = 1,NMAX
          if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
             if (ires.le.(NMAX/2)) then
                res(IRES) = res(IRES) * fac_Z_up      ! factor fac_Z_up for up-type with all combinations
             elseif (ires.gt.(NMAX/2)) then
                res(IRES) = res(IRES) * fac_Z_down    ! factor fac_Z_down for down-type with all combinations
             endif
             if (mod(nmax,2).ne.0) then
               print*, "Something went wrong with the hadronic decays!"
               stop
             endif
          elseif (mod(abs(finalquarks(1)),2).eq.0) then
             res(IRES) = res(IRES) * fac_Z_up         ! up-type final states
          else
             res(IRES) = res(IRES) * fac_Z_down       ! down-type final states
          endif
        ENDDO
      END SELECT


      if (init.eq.1) then
         init = init+1
         if (lwarn) print 199," proc #s for Z/Hjj are ",nproc
 199     format(a,4i5)
      endif

      m2s = 0
      DO IRES = 1,NMAX
         m2s = m2s + RES(IRES)
      ENDDO

c...Les Houches interface - the most propable subprocess   
      if ((lha.or.hepmc).and..not.doNLO) then
         i=0
         weight=0.d0
         rnumb=RandomNumber()
         do while((i.le.nmax).and.(weight.le.rnumb*m2s))
            i=i+1
            weight=weight+res(i)
            iprocess=i
         enddo
         SELECT CASE(process)
         CASE(WPhadWMjj)
            if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
               rnumb=RandomNumber()
               finalquarks_psp(1) =  2 + 2* INT(rnumb*2)
               finalquarks_psp(2) = -1 - 2* INT(rnumb*2)
            endif
         CASE(WPWMhadjj)
            if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
               rnumb=RandomNumber()
               finalquarks_psp(1) =  1 + 2* INT(rnumb*2)
               finalquarks_psp(2) = -2 - 2* INT(rnumb*2)
            endif
         CASE(ZZhadjj)
            if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
               rnumb=RandomNumber()
               if (i.le.(nmax/2)) then    ! up-type
                 finalquarks_psp(1) =  2 + 2* INT(rnumb*N_gen_up)
                 finalquarks_psp(2) = -2 - 2* INT(rnumb*N_gen_up)
               else                       ! down-type
                 finalquarks_psp(1) =  1 + 2* INT(rnumb*N_gen_down)
                 finalquarks_psp(2) = -1 - 2* INT(rnumb*N_gen_down)
               endif
            endif
         END SELECT
      endif


c      DO IRES = nproc(1)+1,nproc(3)
c         m2s = m2s + RES(IRES)
c      ENDDO
    
      if (bos.eq.2) then
         if (jsig.eq.0) then
            m2s_wbfvv = m2s     !*3  ! factor 3 for electrons and muons and tau
         else
            m2s_wbfvv = m2s*8d0   ! factor 8 for random helicity summation

            tauchi(1)=jsig5
            tauchi(2)=jsig5
         endif
      elseif (bos.eq.6) then   
         if (rn.gt.0.5d0) then
            tauchi(1)= 1
         else
            tauchi(1)=-1
         endif
         tauchi(2)= -1*tauchi(1)

         m2s_wbfvv = m2s
      elseif (bos.eq.61 .or. bos.eq.16) then   
         if (rn.gt.0.5d0) then
            tauchi(1)= 1
         else
            tauchi(1)=-1
         endif
         tauchi(2)= -1*tauchi(1)

	 
         m2s_wbfvv = m2s

      elseif ((bos.eq.22).or.(bos.eq.722).or.bos.eq.211) then ! qq->ZZqq -> 4l
         m2s_wbfvv = m2s*4d0   ! factor 4 for random helicity summation
      elseif (bos.eq.11) then ! qq (-> spin-2)-> gamma gamma qq
         m2s_wbfvv = m2s*4d0   ! factor 4 for random helicity summation
      elseif ((bos.eq.21).or.(bos.eq.721).or.bos.eq.212) then ! qq->ZZqq -> 2l 2v
         m2s_wbfvv = m2s*2d0   ! factor 2 for random helicity summation
      else
         m2s_wbfvv = m2s
      endif
      if(photon_hel.ne.0) m2s_wbfvv = m2s_wbfvv/2d0

      if ((lha.or.hepmc).and..not.doNLO) then
         if (bos.eq.2) then
            if (jsig.ne.0) then                  
               helicity(2)= jsig5 !particle is at the 2nd place
               helicity(1)=-jsig5 
            else
               helicity(2)= 9
               helicity(1)= 9
            endif
         endif
         if (bos.eq.1) then
            helicity(1) = 9
         endif
      endif

      RETURN
      END

c******************************************************************
c
c   begin subroutine wbf_zh
c
c*****************************************************************
      subroutine wbf_zh(xi,p,v,rn,physToDiag,fsign,bos, nlo, ps_number,
     1                  pdf,Qpdf,res,nmin,nmax)
          use globalvars, only: ldoblha

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/process.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

      real*8 p(0:3,max_p), v(0:3,max_v), xi(nx)
      real*8 pdf(-6:6,2), res(maxnumsubproc)
      integer physToDiag(4), fsign(4+max_v), bos, h,nlo, nmin, nmax,ps_number
c      
      common /hcount / h
c
c wbf_zh calls the amplitude square routines 
c             qqzqq    for qq-->qqZ    for bos=2      and
c             qqhqq    for qq-->qqH    for bos=6
c             qqwwqq   for qq-->qqWW   for bos=34.or.43
c             qqzzqq   for qq-->qqZZ   for bos=22 or 21
c             qqzaqq   for qq-->qqZA   for bos=211 or 212
c             qqsp2aaqq   for qq(-->spin2)-->qqaa    for bos=11
c             qqhhqq   for qq-->qqHH   for bos=66
c
c  INPUT:  p(0:3,4)        external physical parton momenta
c          v(0:3,nv)       Z decay momenta, for Higgs production only
c                          the sum q(mu) = v(mu,1)+...+v(mu,nv) is needed
c          physToDiag(4)   physToDiag(ext.mom. label) = Feynman diagram label
c          fsign           sign factors for amplitude calls; see qqZqq
c          nlo             nlo = 0: calculate LO only
c                          nlo = 1: calculate full NLO (Born + interference)
c	   h		   specify lepton helicity combination
c  OUTPUT:
c          uucc(ID)   output in format as in qqZqq(j), but corresponding
c          etc.            to res(*) in m2s_qqZqq
c     or   res(k)   
c
c  In and Output
c    nmin, nmax            range of process numbers for this call
c
      real*8 N                              ! color factors
      parameter(N=3d0)

c Note: Factor 9 from color sum included in qqbqq
c
c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"

c local variables
      real*8 uucc, uuss, ddcc, ddss, udsc, ducs
      real*8 pbar(0:3,4+max_v), pA(0:3), dotrr, polcol
      double precision NCmatrixelt(0:1,0:1)
      double precision CCmatrixelt(0:1)
      integer iflav(4), diagToPhys(4), FL_ZH
      external dotrr, FL_ZH

      logical ChargedCurrent, sametype, oneAntiparticle
      logical ldebug
      integer i,if1,if2,icc1,k,mu
      parameter (ldebug=.false.)

      LOGICAL ldebug_ltens
      PARAMETER (ldebug_ltens=.false.)

*** start : sophy qed 
** d,u,s,c,b relative electric charges 
      double precision QuarkECharge(1:5),Qf(1:5),sigmaf(1:4) !
      data QuarkECharge /-1.0D0,2.0D0,-1.0D0,2.0D0,-1.0D0/
      logical run_ew

** photon mass squared (massG) and quark mass (massQ)
#include "VBFNLO/utilities/mssm.inc"
*      real*8 fmass2, phtnmass2

** tree level amplitudes
      double precision NCborn(0:1,0:1)
      double precision CCborn(0:1)
      double precision tree(6)

** PDF array for qed part
      real*8 Qpdf(-7:7,2,3)

** random number 
      real*8 rn !, intent(int)


c$$$      write(*,*)'p(1) =', p(0,1), p(1,1), p(2,1), p(3,1)
c$$$      write(*,*)'p(2) =', p(0,2), p(1,2), p(2,2), p(3,2)
c$$$      write(*,*)'p(3) =', p(0,3), p(1,3), p(2,3), p(3,3)
c$$$      write(*,*)'p(4) =', p(0,4), p(1,4), p(2,4), p(3,4)
c$$$
c$$$      write(*,*)'v(1) =', v(0,1), v(1,1), v(2,1), v(3,1)
c$$$      write(*,*)'v(2) =', v(0,2), v(1,2), v(2,2), v(3,2)
c$$$      write(*,*)'v(3) =', v(0,3), v(1,3), v(2,3), v(3,3)
c$$$      write(*,*)'v(4) =', v(0,4), v(1,4), v(2,4), v(3,4)
     
      uucc = 0d0
      uuss = 0d0
      ddcc = 0d0
      ddss = 0d0
      udsc = 0d0
      ducs = 0d0

      nmin = nmax+1
      do i = 1,4
         diagToPhys(physToDiag(i)) = i
      enddo
      polcol = 1d0/(4d0*N**2*xi(1)*xi(2))
         
c get the real emission amplitude squared, store it in uucs(-1,1) etc.
      do mu = 0,3
         do i = 1,4
            pbar(mu,physToDiag(i))=p(mu,i)
         enddo
         do i = 5,4+max_v
            pbar(mu,i)=0d0
         enddo
      enddo

      if (bos.eq.2) then
         do mu = 0,3             ! kinematics for Z-->l+l- decay
            pbar(mu,5) = v(mu,1) ! ebar
            pbar(mu,6) = v(mu,2) ! e
         enddo
         fsign(5) = -1
         fsign(6) = 1
         call qqzqq(pbar,fsign, nlo,
     1              uucc,uuss,ddcc,ddss,udsc,ducs)  
      elseif (bos.eq.1) then
         do mu = 0,3 
            pbar(mu,5) = v(mu,n_v) ! photon
            pbar(mu,6) = 0d0
         enddo
         fsign(5) = -1
         fsign(6) = 0
         call qqAqq(pbar,fsign, nlo,
     1              uucc,uuss,ddcc,ddss,udsc,ducs)         

      elseif (bos.eq.6) then
         do mu = 0,3
            pbar(mu,5) = 0      ! dummy momentum
            pbar(mu,6) = 0
            do i = 1,n_v
               pbar(mu,6) = pbar(mu,6) + v(mu,i) ! Higgs momentum
            enddo
         enddo
         fsign(5) = 0
         fsign(6) = 1       
         call qqhqq(pbar,fsign, nlo, 
     1              uucc,uuss,ddcc,ddss,udsc,ducs,tree)
      elseif (bos.eq.61 .or. bos.eq.16) then   ! HAjj
         do mu=0,3
            pbar(mu,1)=p(mu,1)
            pbar(mu,2)=p(mu,3)
            pbar(mu,3)=p(mu,2)
            pbar(mu,4)=p(mu,4)
         enddo
         do mu = 0,3
            pbar(mu,5) = 0      ! dummy momentum
            pbar(mu,6) = 0
            do i = 1,n_v-1
               pbar(mu,6) = pbar(mu,6) + v(mu,i) ! Higgs momentum
            enddo
         enddo
         fsign(5) = 0
         fsign(6) = 1
         do mu=0,3
            pA(mu)=v(mu,n_v)
         enddo
         CALL HAjj_ME(pbar,pA,fsign, nlo, 
     1              uucc,uuss,ddcc,ddss,udsc,ducs)
      elseif (bos.eq.34 .or. bos.eq.43) then         ! W+W- to 4 leptons
         do mu = 0,3             ! kinematics for H-->WW -->4 lepton decay
            pbar(mu,5) = v(mu,2) ! l+
            pbar(mu,6) = v(mu,1) ! nu
            pbar(mu,7) = v(mu,3) ! nubar
            pbar(mu,8) = v(mu,4) ! l-
         enddo
         fsign(5) = -1
         fsign(6) = 1
         fsign(7) = -1
         fsign(8) = 1
c         call qqhqq(pbar,fsign, nlo, 
c     1              uucc,uuss,ddcc,ddss,udsc,ducs,tree)
         if (lwarn) then
            call qqwwqq_comp(pbar,fsign, nlo,
     1                       uucc,uuss,ddcc,ddss,udsc,ducs)
         endif
         call qqwwqq(pbar,fsign, nlo, 1,
     1                    uucc,uuss,ddcc,ddss,udsc,ducs)
#ifdef WITH_MADGRAPH
        if (ldebug_ltens) call qqwwqq_mg2(p,pbar,v,fsign, nlo, 
     1                     uucc,uuss,ddcc,ddss,udsc,ducs)
#endif
c         call qqwwqq_comp(pbar,fsign, nlo, 
c     1                    uucc,uuss,ddcc,ddss,udsc,ducs)
      elseif (bos.eq.22) then         ! ZZ to 4 leptons
         do mu = 0,3             ! kinematics for H-->ZZ -->4 lepton decay
            pbar(mu,5) = v(mu,1) ! l+
            pbar(mu,6) = v(mu,2) ! l-
            pbar(mu,7) = v(mu,3) ! l"+
            pbar(mu,8) = v(mu,4) ! l"-
         enddo
         fsign(5) = -1
         fsign(6) = 1
         fsign(7) = -1
         fsign(8) = 1
         if (lwarn) then
c            call qqzzqq_comp(pbar,fsign, nlo,  bos,
c     1                       uucc,uuss,ddcc,ddss,udsc,ducs)
         endif
         call qqzzqq(pbar,fsign, nlo, 1,bos,
     1                    uucc,uuss,ddcc,ddss,udsc,ducs)

#ifdef WITH_MADGRAPH
        if (ldebug_ltens) call qqzzqq_mg2(p,pbar,v,fsign, nlo, bos,
     1                     uucc,uuss,ddcc,ddss,udsc,ducs)
#endif
c        call qqzzqq_comp(pbar,fsign, nlo, bos,
c     1                    uucc,uuss,ddcc,ddss,udsc,ducs)

      elseif (bos.eq.66) then         ! HH production
         do mu = 0,3             ! kinematics
            if(procid.eq.HHjj) then
               pbar(mu,5) = v(mu,1) ! H1
               pbar(mu,6) = v(mu,2) ! H2
            else ! including H decays
               pbar(mu,5) = v(mu,1)+v(mu,2) ! H1*
               pbar(mu,6) = v(mu,3)+v(mu,4) ! H2*
            endif
         enddo
         fsign(5) = 1
         fsign(6) = 1
         call qqhhqq(pbar,fsign, nlo, 1,bos,
     1                    uucc,uuss,ddcc,ddss,udsc,ducs)

      elseif (bos.eq.11) then         ! qq(-->spin-2)-->qq gamma gamma
         do mu = 0,3             
            pbar(mu,5) = v(mu,1) ! photon 1
            pbar(mu,6) = v(mu,2) ! photon 2
         enddo

         fsign(5) = 1
         fsign(6) = 1

         call qqsp2aaqq(pbar,fsign, nlo, 1,bos,
     1                    uucc,uuss,ddcc,ddss,udsc,ducs)

      elseif (bos.eq.734) then   ! spin-2 --> W+W- --> 4 leptons
         do mu = 0,3             ! kinematics for spin-2-->WW -->4 lepton decay
            pbar(mu,5) = v(mu,2) ! l+
            pbar(mu,6) = v(mu,1) ! nu
            pbar(mu,7) = v(mu,3) ! nubar
            pbar(mu,8) = v(mu,4) ! l-
         enddo
         fsign(5) = -1
         fsign(6) = 1
         fsign(7) = -1
         fsign(8) = 1

         call qqsp2wwqq(pbar,fsign, nlo, 1,
     1                    uucc,uuss,ddcc,ddss,udsc,ducs)

      elseif (bos.eq.722) then   ! spin-2 --> ZZ --> 4 leptons
         do mu = 0,3             !  kinematics for spin-2-->ZZ -->4 lepton decay
            pbar(mu,5) = v(mu,1) ! l+
            pbar(mu,6) = v(mu,2) ! l-
            pbar(mu,7) = v(mu,3) ! l"+
            pbar(mu,8) = v(mu,4) ! l"-
         enddo
         fsign(5) = -1
         fsign(6) = 1
         fsign(7) = -1
         fsign(8) = 1

         call qqsp2zzqq(pbar,fsign, nlo, 1,bos,
     1                    uucc,uuss,ddcc,ddss,udsc,ducs)

      elseif (bos.eq.721) then         ! spin-2 --> ZZ to 2 leptons + 2 neutrinos
         do mu = 0,3             ! kinematics for spin-2-->ZZ-->2l 2v decay
            pbar(mu,5) = v(mu,1) ! l+
            pbar(mu,6) = v(mu,2) ! l-
            pbar(mu,7) = v(mu,3) ! v
            pbar(mu,8) = v(mu,4) ! v~
         enddo
         fsign(5) = -1
         fsign(6) =  1
         fsign(7) =  1
         fsign(8) = -1
        
	call qqsp2zzqq(pbar,fsign, nlo, 1,bos,
     1                    uucc,uuss,ddcc,ddss,udsc,ducs) 

      elseif (bos.eq.21) then         ! ZZ to 2 leptons + 2 neutrinos
         do mu = 0,3             ! kinematics for H-->ZZ-->2l 2v decay
            pbar(mu,5) = v(mu,1) ! l+
            pbar(mu,6) = v(mu,2) ! l-
            pbar(mu,7) = v(mu,3) ! v
            pbar(mu,8) = v(mu,4) ! v~
         enddo
         fsign(5) = -1
         fsign(6) =  1
         fsign(7) =  1
         fsign(8) = -1
         
         if (lwarn) then
c            call qqzzqq_comp(pbar,fsign, nlo,  bos,
c     1                       uucc,uuss,ddcc,ddss,udsc,ducs)
         endif
        
	call qqzzqq(pbar,fsign, nlo, 1,bos,
     1                    uucc,uuss,ddcc,ddss,udsc,ducs) 

c        call qqzzqq_comp(pbar,fsign, nlo, bos,
c     1                    uucc,uuss,ddcc,ddss,udsc,ducs)
      elseif (bos.eq.211) then         ! ZA to 2 leptons + A
         do mu = 0,3             ! kinematics for qq-->ZA-->2l A decay
            pbar(mu,5) = v(mu,1) ! l+
            pbar(mu,6) = v(mu,2) ! l-
            pbar(mu,7) = v(mu,3) ! a
         enddo
         fsign(5) = -1
         fsign(6) =  1
         fsign(7) =  1
!          udsc=temp
	call qqzaqq(pbar,fsign, nlo, 1,bos,ps_number,
     1                    uucc,uuss,ddcc,ddss,udsc,ducs)

!        uucc=0d0
!        uuss=0d0
!        ddcc=0d0
!        ddss=0d0
!        ducs=0d0
!#ifdef WITH_MADGRAPH
cfc         if(.true.) then
cfc	  pbar(0:3,1:4) = p(0:3,1:4)
cfc	  call MG_ZA_UDSC(Pbar,temp)
cfc         endif
! 	print*,fsign,nlo,bos
cfc        print*,"udsc",temp,udsc,temp/udsc
cfc         if(.true.) then
cfc	  pbar(0:3,1:4) = p(0:3,1:4)
cfc          temp=0
cfc	  call MG_ZA_UUCC(Pbar,fsign,temp)
cfc         endif
cfc        print*, "uucc",temp,uucc,temp/uucc
cfc        pause
!#endif
cjp  
      if (.false.) then
        print*, "uucc", uucc
        print*, "uuss", uuss
        print*, "ddcc", ddcc
        print*, "ddss", ddss
        print*, "udsc", udsc
        print*, "ducs", ducs
c        pause
      endif
      
      elseif (bos.eq.212) then         ! ZA to 2 leptons + A
         do mu = 0,3             ! kinematics for qq-->ZA-->2l A decay
            pbar(mu,5) = v(mu,1) ! nu
            pbar(mu,6) = v(mu,2) ! nu~
            pbar(mu,7) = v(mu,3) ! a
         enddo
         fsign(5) = +1
         fsign(6) = -1
         fsign(7) =  1
!          udsc=temp
	call qqzaqq(pbar,fsign, nlo, 1,212,ps_number,
     1                    uucc,uuss,ddcc,ddss,udsc,ducs)
!       if (.true.) then
!         print*, "uucc", uucc
!         print*, "uuss", uuss
!         print*, "ddcc", ddcc
!         print*, "ddss", ddss
!         print*, "udsc", udsc
!         print*, "ducs", ducs
! c        pause
!       endif
      else
         print*," Invalid entry BOS = ",bos
         print*," in call of wbf_zh. Must be 2,6, 34/43, 22 or 21"
         stop
      endif

      if (ldoblha) then
        return
      endif

      NCmatrixelt(0,0)=uucc
      NCmatrixelt(0,1)=uuss
      NCmatrixelt(1,0)=ddcc
      NCmatrixelt(1,1)=ddss
      CCmatrixelt(0)=udsc
      CCmatrixelt(1)=ducs

*** start sophy qed : setting born amplitude
      if(nlo.ge.1) then
         NCborn(0,0)=tree(1)    !uucc
         NCborn(0,1)=tree(2)    !uuss
         NCborn(1,0)=tree(3)    !ddcc
         NCborn(1,1)=tree(4)    !ddss
         CCborn(0)=tree(5)      !udsc
         CCborn(1)=tree(6)      !ducs
      else
         NCborn(0,0)=0
         NCborn(0,1)=0
         NCborn(1,0)=0
         NCborn(1,1)=0
         CCborn(0)=0
         CCborn(1)=0
      endif
*** end 


      do if1=1,nflVBF             !(nfl/2)*2
         do if2=1,nflVBF          !(nfl/2)*2
            iflav(1)=if1*fsign(physToDiag(1))
            iflav(3)=if1*fsign(physToDiag(3))
            iflav(2)=if2*fsign(physToDiag(2))
            iflav(4)=if2*fsign(physToDiag(4))
            k=FL_ZH(iflav,1)
            if (ldebug) 
     &           print*,k,":",iflav(1),iflav(2),"-->",iflav(3),iflav(4)

            res(k)=pdf(sign1*iflav(1),1)
     &           *pdf(sign2*iflav(2),2)
     &           *NCmatrixelt(mod(if1,2),mod(if2,2))*polcol

            run_ew = .false.
            if (ewcor_switch) then
               if (nlo .eq. 2) then
                  run_ew = .true.
               end if
            end if
*** start sophy qed: For Higgs production, with full ew corrections ONLY!
#ifdef WITH_LT
            if((bos .eq. 6) .and. run_ew .and. (sector .ge. 3)) then
c     res(k) = 0       !LO is zero for checking 
               
** charge flow signs
               sigmaf(1) = dble(fsign(physToDiag(1))) 
               sigmaf(3) = -dble(fsign(physToDiag(3))) 
               sigmaf(2) = dble(fsign(physToDiag(2))) 
               sigmaf(4) = -dble(fsign(physToDiag(4))) 

** set up charges
               Qf(1) = QuarkECharge(abs(iflav(1)))/3.0d0
               Qf(2) = QuarkECharge(abs(iflav(2)))/3.0d0
               Qf(3) = QuarkECharge(abs(iflav(3)))/3.0d0
               Qf(4) = QuarkECharge(abs(iflav(4)))/3.0d0
               
               call finitecollinear(Qpdf,xi(1),rn,p,1,Qf,sigmaf,
     $              mufsq(1,1),massQ**2,massG)             
               call finitecollinear(Qpdf,xi(2),rn,p,2,Qf,sigmaf,
     $              mufsq(1,1),massQ**2,massG) 
               
c     include finite collinear pieces 
               res(k) = res(k) + (Qpdf(sign1*iflav(1),1,1)*
     &              Qpdf(sign2*iflav(2),2,3) + Qpdf(sign1*iflav(1),1,3)*
     &              Qpdf(sign2*iflav(2),2,1))*
     &              NCborn(mod(if1,2),mod(if2,2))*polcol
            endif
*** end sophy qed
#endif

C Now check if there is a CC contribution for this choice of initial state
C flavors (i,j). First: check whether initial uu or dd, not ud
            sametype=(mod(if1,2)).eq.(mod(if2,2))
            oneAntiparticle=
     &           ((fsign(physToDiag(1))*fsign(physToDiag(2))).eq.-1) 
c true if only one particle is an antiparticle
            ChargedCurrent=(oneAntiparticle.and.sametype) .or.
     &           ( (.not.oneAntiparticle).and.(.not.sametype) )
            if ( ChargedCurrent .and. if1.le.4 .and. if2.le.4 ) then   
! external bottom quarks are never allowed with charged currents
c     change 1<-->2 and 3<-->4 in outgoing quark flavors
               iflav(3)=(if1+2*mod(if1,2)-1)*fsign(physToDiag(3))
               iflav(4)=(if2+2*mod(if2,2)-1)*fsign(physToDiag(4))
               icc1 = abs(iflav(diagtophys(1)))
               k=FL_ZH(iflav,1)
               if (ldebug) 
     &              print*,k,":",iflav(1),iflav(2),"-->",iflav(3),iflav(4)
               res(k)=pdf(sign1*iflav(1),1)
     &              *pdf(sign2*iflav(2),2)
     &              *CCmatrixelt(mod(icc1,2))*polcol
#ifdef WITH_LT
***   start: sophy new qed: For Higgs production, with full ew corrections ONLY!
               if((bos .eq. 6) .and. run_ew .and. (sector .ge. 3)) then

** set up charges
                  Qf(1) = QuarkECharge(abs(iflav(1)))/3.0d0
                  Qf(2) = QuarkECharge(abs(iflav(2)))/3.0d0
                  Qf(3) = QuarkECharge(abs(iflav(3)))/3.0d0
                  Qf(4) = QuarkECharge(abs(iflav(4)))/3.0d0
                     
c     res(k) = 0    !LO is zero for checking 
        
                  call finitecollinear(Qpdf,xi(1),rn,p,1,Qf,sigmaf,
     $                 mufsq(1,1),massQ**2,massG) 
                  call finitecollinear(Qpdf,xi(2),rn,p,2,Qf,sigmaf,
     $                 mufsq(1,1),massQ**2,massG) 

c     include finite collinear piece's 
                  res(k) = res(k) + (Qpdf(sign1*iflav(1),1,3)
     &                 *Qpdf(sign2*iflav(2),2,1) +
     &                 Qpdf(sign1*iflav(1),1,1)*
     &                 Qpdf(sign2*iflav(2),2,3))*CCborn(mod(icc1,2))*
     &                 polcol
               endif
#endif
*** end: sophy new qed
            endif               ! endif(CC)
         enddo
      enddo 

      nmax = FL_ZH(iflav,-2)

      end
c******************************************************************
c
c   end subroutine wbf_zh
c
c*****************************************************************


c******************************************************************
c
c   begin function FL_ZH
c
c*****************************************************************
      integer function FL_ZH(iflav,colstruc)
      implicit none
      integer iflav(4),colstruc    ! input for color and flavor assignment
c  fill hepup color and flavor assignments for WBF processes with
c  two quark lines only

c  id1,2 are flavor identifiers for incoming quarks
c  id3,4 are flavor identifiers for outgoing quarks
c  colstruc = is irrelevant here
c
c  color structure of these interactions is trivial
c  and colstruc was included only to make this program 
c  more consistent with the others.

c  Note that colstruc is also a flag:
c  if colstruc=-1, then we are resetting.

c  The color structure is always the same:
c  color always flows from id1->id3 and id2->id4
c  no gluons are involved.

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
      
      integer lkup
      common /localkup/ lkup(numParticles,7,maxNumSubProc)
 
      integer listposition
      integer numresets
      save listposition
      save numresets

      data numresets /-1/
      data listposition /0/

      if (colstruc.eq.-1) then  ! we are resetting
c       write(*,*) "we are resetting"
        listposition=0
        numresets=numresets+1
      else if (colstruc.eq.-2) then
c       don"t do anything, just return the number of subprocesses.
      else
c increment the counter regardless of whether or not it"s the 1st time through.
        listposition=listposition+1

c       we fill info for this subprocess,if required
        if(numresets.eq.0) then
          call fillColoredPartons_ZH(iflav,
     1                               listposition+numdecay)
        endif       ! numresets.eq.0
      endif         
      FL_ZH=listposition
      end
c******************************************************************************
c
c   end function FL_ZH
c
c******************************************************************************


c*****************************************************************************
c
c    begin  subroutine fillColoredPartons_ZH
c
c*****************************************************************************
      subroutine fillColoredPartons_ZH(iflav,listposition)
c  assigns values to the variables in the common block localHEPUP
c  in particular, this subroutine assigns values to those variables that 
c  will be stored in the lookup tables generated by writeHEPUPtable.  
c  As the name suggests, this routine only stores the information for the 
c  colored partons.  Particles without color will be dealt with in the 
c  subroutine fillColorless.

      implicit none
c  include files
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/process.inc"

c  input parameters 
      integer iflav(4),id1,id2,id3,id4,listposition   

c  local variables
      integer i

c 4 quark flavors and subprocess number      
      select case(process)
      case (WMWMjj,WMhadWMjj)         ! WMWMjj is obtained from W+W+jj, flip IDs
        id1=-iflav(1)                 ! 4 quarks
        id2=-iflav(2)
        id3=-iflav(3)
        id4=-iflav(4)
      case default
        id1=iflav(1)                  ! 4 quarks
        id2=iflav(2)
        id3=iflav(3)
        id4=iflav(4)
      end select

   
      select case(process)
      case(Hjj_WW, Hjj_ZZ_ll, Hjj_ZZ_lnu, WPWMjj, ZZjj_ll, ZZjj_lnu,
     &      WPWPjj, WMWMjj,
     &      WPhadWMjj, WPWMhadjj, ZZhadjj,
     &      Hjj_WPhadWM, Hjj_WPWMhad, Hjj_ZZhad,
     &      WPhadWPjj, WMhadWMjj, Sp2jj_WW)
         lnup(listposition)=numParticles+2
      case(ZAjj,ZAjj_n) 
         lnup(listposition)=numParticles
      case(Hjj_AA, Hjj_mu, Hjj_tau, Hjj_bbar, Zjj_l, Zjj_nu, WPjj,
     &      WMjj, AAjj) 
         lnup(listposition)=numParticles-1
      case(Hjj,Ajj)             
         lnup(listposition)=numParticles-3
      case(HHjj) ! modif julien
         lnup(listposition)=numParticles-2
      case(HHjj_bbtau,HHjj_bbAA) ! modif julien
         lnup(listposition)=numParticles+2
      case(HAjj)
         lnup(listposition)=numParticles-2
      case(HAjj_AA, HAjj_mu, HAjj_tau, HAjj_bbar)
         lnup(listposition)=numParticles
      case(HAjj_WW, HAjj_ZZ_ll, HAjj_ZZ_lnu)
         lnup(listposition)=numParticles+3
      end select


      listup(1,listposition)=-1   !incoming quarks
      listup(2,listposition)=-1
      listup(3,listposition)=1    !outgoing quarks
      listup(4,listposition)=1
 
      lidup(1,listposition)=id1     ! the 4 quarks
      lidup(2,listposition)=id2
      lidup(3,listposition)=id3
      lidup(4,listposition)=id4

      do i=3,n_p
        lmothup(1,i,listposition)=1
        lmothup(2,i,listposition)=2
      enddo

      if (id1.gt.0) then
        licolup(1,1,listposition)=501
        licolup(1,3,listposition)=501
        licolup(2,1,listposition)=0
        licolup(2,3,listposition)=0
      else
        licolup(1,1,listposition)=0
        licolup(1,3,listposition)=0
        licolup(2,1,listposition)=501
        licolup(2,3,listposition)=501
      endif

      if (id2.gt.0) then
        licolup(1,2,listposition)=502
        licolup(1,4,listposition)=502
        licolup(2,2,listposition)=0
        licolup(2,4,listposition)=0
      else
        licolup(1,2,listposition)=0
        licolup(1,4,listposition)=0
        licolup(2,2,listposition)=502
        licolup(2,4,listposition)=502
      endif

      end  ! fillColoredPartons_ZH
c*****************************************************************************
c
c    end subroutine fillColoredPartons_ZH
c
c*****************************************************************************

      subroutine resprint2(nmin,nmax,res)
      implicit none
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
      integer nmin,nmax,k,ic
      save ic
      real*8 resold
      common /bresprint/ resold(16)
      real*8 xa,xb,res(maxnumsubproc),fac
      if (nmin.eq.1) ic=1
      xa = 0
      xb = 0
      if (nmin.le.190) then 
         fac = 4*9
      else
         fac = 4d0*24d0
      endif
      do k = nmin,nmax
         xa = xa + res(k)
      enddo
      print*," res([",nmin,",",nmax,"])",xa,
     1         xa/(resold(ic)+resold(ic+1))*fac
      ic = ic+2
      end

c --------------------------------------------------------------------------

      subroutine vtoww_reset
      implicit none
#include "tensor.inc"
      complex*16 zero
      parameter (zero=(0d0,0d0))
      integer j,jj,mu,nu
      do mu = 0,3
         aww(mu) = zero
         zww(mu) = zero
         do nu = 0,3
            do j = 1,3
               aaww(mu,nu,j) = zero
               azww(mu,nu,j) = zero
               zaww(mu,nu,j) = zero
               zzww(mu,nu,j) = zero
               wwww5(mu,nu,j) = zero
               wwww6(mu,nu,j) = zero
               do jj = 1,2
                  NCwpa(mu,nu,jj,j) = zero 
                  NCwpz(mu,nu,jj,j) = zero
                  CCwpa(mu,nu,jj,j) = zero 
                  CCwpz(mu,nu,jj,j) = zero
                  NCwma(mu,nu,jj,j) = zero 
                  NCwmz(mu,nu,jj,j) = zero
                  CCwma(mu,nu,jj,j) = zero 
                  CCwmz(mu,nu,jj,j) = zero 
               enddo
            enddo
         enddo
      enddo
      print*," v(v) to WW tensors initialized to ",zero 
      return
      end

c --------------------------------------------------------------------------

      subroutine vtozz_reset
      implicit none
#include "tensorz.inc"
      complex*16 zero
      parameter (zero=(0d0,0d0))
      integer j,jj,mu,nu
      do mu = 0,3
         azz(mu) = zero
         zzztens(mu) = zero
         do nu = 0,3
            do j = 1,3

               wwzz5(mu,nu,j) = zero
               wwzz6(mu,nu,j) = zero
	       aazz(mu,nu,j) = zero
	       azzz(mu,nu,j) = zero
	       zazz(mu,nu,j) = zero
	       zzzz(mu,nu,j) = zero
	       
	       do jj = 1,2
	       
                  aaee(mu,nu,jj,j) = zero
                  azee(mu,nu,jj,j) = zero
                  zaee(mu,nu,jj,j) = zero
                  zzee(mu,nu,jj,j) = zero
                  aauu(mu,nu,jj,j) = zero
                  azuu(mu,nu,jj,j) = zero
                  zauu(mu,nu,jj,j) = zero
                  zzuu(mu,nu,jj,j) = zero
                  CCee(mu,nu,jj,j) = zero
                  CCuu(mu,nu,jj,j) = zero 
                  CCee6(mu,nu,jj,j) = zero 
                  CCuu6(mu,nu,jj,j) = zero
	       
	       enddo            !jj

            enddo               !j
         enddo                  !nu
      enddo                     !mu
      
      print*," v(v) to ZZ tensors initialized to ",zero 
      return
      end

c --------------------------------------------------------------------------

      subroutine vtosp2ww_reset
      implicit none
#include "tensorsp2aa.inc"
#include "VBFNLO/utilities/global.inc"

      complex*16 zero
      parameter (zero=(0d0,0d0))
      integer mu,nu,j
      do mu = 0,3
         do nu = 0,3
            do j = 1,3

               sp2wwww5(mu,nu,j) = zero
               sp2wwww6(mu,nu,j) = zero
	       sp2aaww(mu,nu,j) = zero
	       sp2azww(mu,nu,j) = zero
	       sp2zaww(mu,nu,j) = zero
	       sp2zzww(mu,nu,j) = zero
	  
            enddo               !j      
         enddo                  !nu
       enddo                    !mu

       if (with_spin2) then        
      print*," VV jj (to spin-2) to W+W- tensors" 
      print*, "initialized to ",zero 
       else
      print*," VV jj (to Higgs) to W+W- tensors" 
      print*, "initialized to ",zero 
       endif
 
      return
      end

c --------------------------------------------------------------------------

      subroutine vtosp2zz_reset
      implicit none
#include "tensorsp2aa.inc"
#include "VBFNLO/utilities/global.inc"

      complex*16 zero
      parameter (zero=(0d0,0d0))
      integer mu,nu,j
      do mu = 0,3
         do nu = 0,3
            do j = 1,3

               sp2wwzz5(mu,nu,j) = zero
               sp2wwzz6(mu,nu,j) = zero
	       sp2aazz(mu,nu,j) = zero
	       sp2azzz(mu,nu,j) = zero
	       sp2zazz(mu,nu,j) = zero
	       sp2zzzz(mu,nu,j) = zero
	  
            enddo               !j      
         enddo                  !nu
       enddo                    !mu

       if (with_spin2) then        
      print*," VV jj (to spin-2) to ZZ tensors" 
      print*, "initialized to ",zero 
       else
      print*," VV jj (to Higgs) to ZZ tensors" 
      print*, "initialized to ",zero 
       endif
 
      return
      end

c --------------------------------------------------------------------------

      subroutine vtohh_reset
      implicit none
#include "tensorz.inc"
      complex*16 zero
      parameter (zero=(0d0,0d0))
      integer j,mu,nu
      do mu = 0,3
c         azz(mu) = zero
c         zzz(mu) = zero
         do nu = 0,3
            do j = 1,3

               wwhh5(mu,nu,j) = zero
               wwhh6(mu,nu,j) = zero
	       zzhh(mu,nu,j)  = zero

            enddo               !j
         enddo                  !nu
      enddo                     !mu
      
      print*," v(v) to HH tensors initialized to ",zero 
      return
      end

c --------------------------------------------------------------------------

      subroutine vtosp2aa_reset
      implicit none
#include "tensorsp2aa.inc"

      complex*16 zero
      parameter (zero=(0d0,0d0))
      integer mu,nu,j
      do mu = 0,3
         do nu = 0,3
            do j = 1,3

               sp2wwaa5(mu,nu,j) = zero
               sp2wwaa6(mu,nu,j) = zero
	       sp2aaaa(mu,nu,j) = zero
	       sp2azaa(mu,nu,j) = zero
	       sp2zaaa(mu,nu,j) = zero
	       sp2zzaa(mu,nu,j) = zero
	  
            enddo               !j      
         enddo                  !nu
       enddo                    !mu

      
      print*," VV jj (to spin-2) to gamma gamma tensors" 
      print*, "initialized to ",zero 
 
      return
      end

#ifdef WITH_LT
*** start sophy qed  - to end
*** new SUBROUTINES needed for the QED part

c*****************************************************************************
c
c convolutions of structure functions with \mathcal{G}_{ij}^{sub} 
c renormalization of structure functions is included
c
c****************************************************************************
      subroutine counterparts_qq(Mf2,MP2,xx,zz,p,ip,jp,gsubij,cff) 
! for initial state quarks 
      implicit none
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/scales.inc"      
! need mu_{f} i.e. factorization scale 
      integer ip,jp             !emitter and spectator particle id
      double precision p(0:3,n_p-1),xx,zz
      double precision Mf2, MP2
      double precision gsubij(3)! integrated dipoles factors 
      double precision cff(2)     !scheme dependent pieces
      double precision zmin
      parameter(zmin=0.999995d0) !cutoff on z for num.stability
      logical ldebug
      parameter(ldebug=.false.)
c
c     cff(1) = C_{ff}(z)
c     cff(2) = Int_{0}^{x} dz C_{ff}(z)
c

c     gsubij(1) = \mathcal{G}^{sub}_{ij} where ij = a1,1a,ab
c     gsubij(2) = Int_{0}^{x} dz \mathcal{G}^{sub}_{ij}(P^2,z)
c     gsubij(3) = g^{sub}_{ij} , this is the endpoint   
c
      double precision dotrr,scriL
      external dotrr,scriL
c
      complex*16 vli2
      external vli2
c     vli2 is called with a complex argument !
c     
c
c
c     local variables
      double precision Pffz,pij,muf2,lnomz,omz,lnz,lnpijomuf2
c
      muf2 = mufsq(1,1) !factorization scale (? make separate for the QCD scale)
c
      Pffz = (1.0d0 + zz**2)/(1.0d0 - zz) !ff splitting function for QED
c      
c     dot product of momenta ip and jp
c
      pij = 2.d0*dotrr(p(0,ip),p(0,jp))

      lnomz = log(1-zz)
      lnz = log(zz)
      lnpijomuf2 = DLOG(abs(pij)/muf2) !Log[Abs[pij]/mu_{f}^{2}]
c
c
c arctanh[z] = 1/2*(ln(1+z)-ln(1-z))
c
      cff(1) = Pffz * (lnomz - lnz - 3.0d0/4.0d0) + (9+5*zz)/4
      cff(2) = (xx*(7 + 2*xx) - 2*xx*(2 + xx)*
     $     0.5d0*(DLOG(2.0d0 - 2.0d0*xx)-DLOG(2.0d0*xx)) -
     $     2*DLOG(1 - xx)*(-3 + DLOG(1 - xx) - 
     $     2*DLOG(xx)) + 4*dble(vli2(dcmplx(xx,0.0d0))))/2.

      if(ldebug) print*,'cff(2)=',cff(2)

c
      if(ip.le.2 .and. jp.le.2) then ! initial-initial 
         if(zz.lt.zmin) then
            lnz = DLOG(zz)      !Log[z]
            omz = 1.0d0 - zz    !1 - z 
            lnomz = DLOG(omz)   !Log[1-z]

            gsubij(1) = Pffz * lnpijomuf2 - Pffz + omz + 
     $           (2.0d0 * lnomz + 1.0d0) * Pffz
         else
            gsubij(1) = 0
         endif
         gsubij(3) = scriL(abs(pij),mf2,mp2) - 
     $        pi**2/3.0d0 + 2.0d0
         
         gsubij(2) = 4*xx - (-1 + xx)*(3 + xx)*DLOG(1 - xx) - 
     $        2*DLOG(1 - xx)**2 -
     $        (DLOG(pij/muf2)*(xx*(2 + xx) + 4*DLOG(1 - xx)))/2.

         if(ldebug) then
         print*,'ip,jp=',ip,jp
         print*,'gsubij(2)=',gsubij(2),xx
         print*,'pij=',pij
         print*,'muf2=',muf2
         endif
         
      elseif(ip.gt.2 .and. jp.le.2) then ! final-initial
         if(zz.lt.zmin) then
            lnz = DLOG(zz)      !Log[z]
            omz = 1.0d0 - zz    !1 - z 
            lnomz = DLOG(omz)   !Log[1-z]
            
            gsubij(1) = (2.0d0 * log(1.0d0 + omz) - 
     $           2.0d0 * lnomz - 3.0d0/2.0d0)/omz
         else
            gsubij(1) = 0
         endif
         gsubij(3) = scriL(abs(pij),mf2,mp2) -
     $        pi**2/2.0d0 + 3.0d0/2.0d0
         
         gsubij(2) = Pi**2/6. + (3*DLOG(1 - xx))/2. + 
     $        DLOG(1 - xx)**2 + 
     $        2*dble(vli2(dcmplx(-1 + xx,0.0d0)))
         if(ldebug) then
         print*,'ip,jp=',ip,jp
         print*,'gsubij(2)=',gsubij(2),xx
         print*,'pij=',pij
         print*,'muf2=',muf2
         endif

      elseif(ip.le.2 .and. jp.gt.2) then ! initial-final

         if(zz.lt.zmin) then
            lnz = DLOG(zz)      !Log[z]
            omz = 1.0d0 - zz    !1 - z 
            lnomz = DLOG(omz)   !Log[1-z]
            gsubij(1) = Pffz * ( lnpijomuf2 - lnz - 1.0d0) -
     $        2.0d0 * log(1.0d0 + omz)/omz + 
     $           (1.0d0 + zz) * lnomz + omz +
     $           (2.0d0*lnomz + 1.0d0) * Pffz
         else
            gsubij(1) = 0
         endif
         gsubij(3) = scriL(abs(pij),mf2,mp2) +
     $        pi**2/6.0d0 - 1.0d0

         gsubij(2) = (Pi**2 + 15*xx + 15*DLOG(1 - xx) - 
     $        6*xx*DLOG(1 - xx) 
     $        - 9*xx**2*DLOG(1 - xx) - 12*DLOG(1 - xx)**2 + 
     $        6*xx*DLOG(xx) + 
     $        3*xx**2*DLOG(xx) - 6*xx*DLOG(Abs(pij)/muf2) - 
     $        3*xx**2*DLOG(Abs(pij)/muf2) - 
     $        12*DLOG(1 - xx)*DLOG(Abs(pij)/muf2) 
     $        - 12*dble(vli2(dcmplx(1 - xx,0.0d0))) - 
     $        12*dble(vli2(dcmplx(-1 + xx,0.0d0))))/6. 
         if(ldebug) then
         print*,'ip,jp=',ip,jp
         print*,'gsubij(2)=',gsubij(2),xx
         print*,'pij=',pij
         print*,'muf2=',muf2
         endif

      elseif(ip.gt.2 .and. jp. gt.2) then !fmass2 and phtnmass2 are regulators
         
         gsubij(1) = 0.0d0

         gsubij(3) = scriL(pij,mf2,mp2) - pi**2/3.0d0 + 
     $        3.0d0/2.0d0

         gsubij(2) = 0.0d0

      endif
    
      return
      end


********************************************************************************
********************************************************************************

** NEW FUNCTION      

      double precision function scriL(P2,Mf2,Mp2) 
      implicit none
      double precision P2       !mom sq
      double precision Mf2      !fermion mass sq 
      double precision Mp2       !photon mass sq
c
      scriL = DLOG(mf2/P2)*DLOG(mp2/P2) + DLOG(mp2/P2)-
     $     1.0d0/2.0d0 * (DLOG(mf2/mp2)**2 - DLOG(mf2/P2))
c 
      return
      end


********************************************************************************
********************************************************************************

      subroutine counterpart_phtn(zz,p,ip,jp,Hsubij) ! for intial state photons
c     mass factorization has been performed 
c     DIS scheme has been adopted
      implicit none
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/scales.inc"      
! need mu_{f} i.e. factorization scale 
      integer ip,jp             !emitter and spectator particle id
      double precision p(0:3,n_p-1),zz,hsubij
c
      double precision dotrr,scriL
      external dotrr,scriL
c
c     local variables
      double precision Pfgammaz,pij,muf2,Cfgammaz
c
      muf2 = mufsq(1,1) !factorization scale (? make separate for the QCD scale)
c      
c     dot product of momenta ip and jp
c
      pij = 2.d0*dotrr(p(0,ip),p(0,jp))

      Pfgammaz = (1.0d0 - zz)**2 + zz**2 !splitting gamma -> ff
      Cfgammaz = Pfgammaz * (DLOG(1.0d0-zz)-DLOG(zz))-8.0d0*zz**2 +
     $     8.0d0*zz - 1.0d0 ! C_{f \gamma}[z] in DIS factorization scheme

      Hsubij = Pfgammaz * (DLOG(pij/muf2)+DLOG((1.0d0-zz)/zz)) + 
     $     2.0d0*zz*(1.0d0-zz) - Cfgammaz

      return
      end

********************************************************************************
********************************************************************************

***  Subroutine to compute x * f_{q/p}^{collinear}(x,muf^2;{Qf,\sigmaf})
c
      subroutine finitecollinear(pdfx,xx,rn,p,ap,Qf,sigmaf,muf2,fmass2,
     $     phtnmass2)

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/mssm.inc"

*      double precision p(0:3,n_p-1)
      double precision p(0:3,max_p)
      integer ap                !initial state emitter 
      double precision pdfx(-7:7,2,3),xx,zz,Qf(4),sigmaf(4),muf2
      double precision rn ! intent(in) ! random number between 0 and 1 for z integration

c     pdfx(flav_a,part_a,1) = xx * f_{a/p}(xx,muf2)
c     pdfx(flav_a,part_a,2) = xx/zz * f_{a/p}(xx/zz,muf2)  
c     pdfx(flav_a,part_a,3) = xx * f_{a/p}^{collinear}(xx,muf2)

      integer k,jp,ip ! flavors
      double precision gsubij(3),cff(2),Hsub
      double precision gsub(1:4,1:4,3),jac
      double precision constii,constff,constif,constfi
      parameter(constii=0.0d0)  !for checking of dipole subtraction
      parameter(constff=0.0d0)  !for checking of dipole subtraction
      parameter(constif=0.0d0)  !for checking of dipole subtraction
      parameter(constfi=0.0d0)  !for checking of dipole subtraction


* start sophy: fermion and photon square masses
      real*8 fmass2, phtnmass2
* end sophy

c      print*,'alpha_qed=',AlfaQED,EL**2/(4.0d0*pi)
      zz = (1.0d0 - xx) * rn**2 + xx !for z integration
c      print*,'zz=',zz
      jac = (1.0d0 - xx) * 2.0d0 * rn ! jacobian factor 

      do ip=1,4
         do jp=1,4
            gsub(ip,jp,1) = 0
            gsub(ip,jp,2) = 0
            gsub(ip,jp,3) = 0
         enddo
      enddo

      do k=-7,7
         pdfx(k,ap,3) = 0.0d0   !initialize to zero
      enddo

      do jp = 3,4   
** initial -final
         call  counterparts_qq(fmass2,phtnmass2,
     $        xx,zz,p,ap,jp,gsubij,cff)

         gsub(ap,jp,1) = gsubij(1)
         gsub(ap,jp,2) = gsubij(2)
         gsub(ap,jp,3) = gsubij(3)
     
** final-initial
         call  counterparts_qq(fmass2,phtnmass2,
     $        xx,zz,p,jp,ap,gsubij,cff)

         gsub(jp,ap,1) = gsubij(1)
         gsub(jp,ap,2) = gsubij(2)
         gsub(jp,ap,3) = gsubij(3)
      enddo

** final-final
      call  counterparts_qq(fmass2,phtnmass2, 
     $     xx,zz,p,3,4,gsubij,cff)

      gsub(3,4,1) = gsubij(1)
      gsub(3,4,2) = gsubij(2)
      gsub(3,4,3) = constff/2.0d0       !gsubij(3)

** initial-initial
      call  counterparts_qq(fmass2,phtnmass2, 
     $     xx,zz,p,1,2,gsubij,cff)

      gsub(1,2,1) = gsubij(1)
      gsub(1,2,2) = gsubij(2)
      gsub(1,2,3) = gsubij(3)

** photon-induced
      call counterpart_phtn(zz,p,ap,ap+1,Hsub) 

c     you will need to worry about the z integration
c     put in jacobian factor is correct places 
c     this can be MC style
c     Check sign factors 

      do k = -5,5               !quark flavors

         pdfx(k,ap,3) = -AlfaQED/(2.0d0*pi)*
     $        (pdfx(k,ap,2) - pdfx(k,ap,1))*jac*
     $        (Qf(ap)*sigmaf(1) * Qf(3) * sigmaf(3)*
     $        (gsub(ap,3,1) + gsub(3,ap,1)) + 
     $        Qf(ap)*sigmaf(ap) * Qf(4) * sigmaf(4)*
     $        (gsub(ap,4,1) + gsub(4,ap,1)) +
     $        Qf(1)*sigmaf(1) * Qf(2)*sigmaf(2) * gsub(1,2,1) 
     $        + Qf(ap)**2 * Cff(1))

         pdfx(k,ap,3) = pdfx(k,ap,3) + AlfaQED/(2.0d0*pi)*
     $        pdfx(k,ap,1) * 
     $        (Qf(ap)*sigmaf(1) * Qf(3) * sigmaf(3)*
     $        (gsub(ap,3,2) + gsub(3,ap,2)) + 
     $        Qf(ap)*sigmaf(ap) * Qf(4) * sigmaf(4)*
     $        (gsub(ap,4,2) + gsub(4,ap,2)) +
     $        Qf(1)*sigmaf(1) * Qf(2)*sigmaf(2) * gsub(1,2,2)
     $        + Qf(ap)**2 * Cff(2)) 


c uncomment the following if you want to use the endpoints
c
c      pdfx(k,ap,3) = pdfx(k,ap,3) - AlfaQED/(2.0d0*pi)*
c     $     pdfx(k,ap,1)*
c     $     (Qf(ap)*sigmaf(1) * Qf(3) * sigmaf(3)*
c     $        (gsub(ap,3,3) + gsub(3,ap,3)) + 
c     $        Qf(ap)*sigmaf(ap) * Qf(4) * sigmaf(4)*
c     $        (gsub(ap,4,3) + gsub(4,ap,3)) +
c     $     (Qf(1)*sigmaf(1) * Qf(2)*sigmaf(2) * gsub(1,2,3)+ ! endpoint contributions (ii,if,fi cases)
c     $     (Qf(3)*sigmaf(3) * Qf(4)*sigmaf(4) * gsub(3,4,3)) ! endpoint contributions (ff case)
cccccccccccccccccc


         pdfx(k,ap,3) = pdfx(k,ap,3) + AlfaQED/(2.0d0*pi)*
     $        pdfx(7,ap,2)*Hsub * jac * !initial photons
     $        Qf(ap)**2 * 3.0d0 ! 3 colors 
c     check consistency of subtraction: initial-initial
         pdfx(k,ap,3) = pdfx(k,ap,3) - pdfx(k,ap,2)*jac*(1.0d0 - zz)*
     $        AlfaQED/(2.0d0*pi)*Qf(1)*sigmaf(1) * Qf(2)*sigmaf(2)*constii
c     initial-final
         pdfx(k,ap,3) = pdfx(k,ap,3) - pdfx(k,ap,2)*jac*
     $        AlfaQED/(2.0d0*pi)* ( Qf(ap)*sigmaf(1) * Qf(3) * sigmaf(3) +
     $        Qf(ap)*sigmaf(ap) * Qf(4) * sigmaf(4)) * constif
c     final-initial 
         pdfx(k,ap,3) = pdfx(k,ap,3) - pdfx(k,ap,2)*jac*
     $        AlfaQED/(2.0d0*pi)* ( Qf(ap)*sigmaf(1) * Qf(3) * sigmaf(3) +
     $        Qf(ap)*sigmaf(ap) * Qf(4) * sigmaf(4)) * constfi

      enddo

      end
#endif



