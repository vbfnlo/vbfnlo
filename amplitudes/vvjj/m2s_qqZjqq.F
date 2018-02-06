c
c	Dieter Zeppenfeld, <dieter@pheno.physics.wisc.edu>
c	Initial version:  2003 February 12
c	Modified by Barbara Jager:   2006 July 11
c       Modified by Julien Baglio:  2014 Oct. 14 (HHjj is included)
c       Last modified by Francisco Campanario:  2015 Jul 2 (ZAjj is included)

c
c  this routine computes |ME|^2 for qq->qq e+ve mu-vm~g ("WW")
c	or qq->qq e+e-mu+mu-g or qq->qq e+e- vm vm~ g ("ZZ")
c	(real emission corr. to qq->qqVV, VV->leptons) 
c	 via call of "qqwwqqj.f" or "qqzzqqj.f"
c
c	flags can be set such that collin./soft subtraction is tested
c
c	for qq->qqZZg random helicity summation is employed for decay leptons
c
c ------------------------------------------------------------------------
c
      subroutine m2s_wbfzhg(
     &                   bos,   !in:  Boson identifier,
c     				!     2=Z,6=H,34/43=WW,22=ZZ->4l,21=ZZ->2l 2v 
                                !     33=W+W+
     &                   nlo,   !in:  NLO=1: create subtraction term; LO = 0
     &                   lok,   !in:  need to calculate m2s(1:3)? T or F 
     &                   xi,	!in:  Feynman x parameters of incoming fermions
     &                   p,	!in:  momenta associated with external fermions
     &                   v,	!in:  momenta associated with WW decay fermions
     &                   rn,    !in:  random number array
     &                   xuz,   !in:  x, u, z for subtraction terms
     &                   m2s    !out: |M|^2*pdf1*pdf2 
     &                        )
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
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/process.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

      integer bos
      real*8 xi(nx), p(0:3,max_p,max_kin), v(0:3,max_v), 
     1       rn(1), xuz(2,2:3), m2s(0:3)
      real*8 vtemp(0:3,max_v)
      logical lok(3)

c  helicity selection
      integer h
      integer jsig, jsig1, jsig3, jsig5
      common /chelsum/ jsig,jsig1,jsig3,jsig5
c
      complex*16 zero
      parameter (zero=(0d0,0d0))
c
c declare external functions
c
      integer FL_ZHg
      external FL_ZHg
c
c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"
      real*8 x1,x2
      real*8 q12(0:4,3), q34(0:4,3), qvv(0:3), !qww(0:3), 
     1 	     lnQomu(2:3), omxi(2:3), 
     1       Ax(2:3), Bx(2:3), Cx(2:3), Dxi(2:3), tgs2oqsq(2:3), 
     2       ln1mxi, lnrat, z, lnz
c
c declare local variables
c
c
      real*8 q_sf
      real*8 dotrr
      external dotrr

      integer init/0/, I, J, mu
      save init
c
c declare variables for summation process
c
      INTEGER ires, nmaxold, FSIGN(4+max_v),
     &        gsign, lsign(3:4)
      save nmaxold

      integer physToDiag(5), nlo, nmin, nmax, nproc(8)
      save nmin, nmax, nproc
c
c store contributions from subprocess ip in res(ip,ID) where
c  ID = 1  : the real emission |M|^2 * pdf
c     = 2,3: sutraction terms for emision off 12 or 34 line
c     = 0  : subtracted result which drives integration, i.e
c res(*,0) = res(*,1)-res(*,2)-res(*,3)
      real*8 res(maxnumsubproc,0:3)

c
c declare parton distribution variables
c
      real*8 pdf(-6:6,2,3)

c variables for hadronic decays
      integer N_gen_W
      real*8 fac_W, mjj2
      external mjj2
      integer N_gen_up, N_gen_down
      real*8 fac_Z_up, fac_Z_down

c
c define program switches
c
      logical ldebug
c      common/ dbug / ldebug	!store "debugging value"
      parameter (ldebug=.false.)
c      data ldebug /.true./		!output debug information
      logical lcoll,lsoft	! parameters for subtraction tests
      parameter (lcoll=.false.,lsoft=.false.)
c      
      common /hcount / h
      real*8 weight,rnumb,RandomNumber

c CPS scheme
      DOUBLE PRECISION CLR,XM2,XMG,Bcoupl,Vcoupl,Acoupl
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),Bcoupl(6,6,6),Vcoupl(4,5),Acoupl(4,5)
      double precision qinv
      double precision qgammaq
      common /VBFNLO_HIGGS_CPS/ qgammaq

      logical upline, lowline, initial
      common/uplow/ upline, lowline, initial  

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
c
c --------------------------------------------------------------
c
c if first iteration, output acceptance cuts and fix constant input 
c parameters for call of qqZWqqji.f
c
      if ( init .ne. 0) then
         do i = 1,nmax
            do j=0,3
               res(i,j) = 0
            enddo
         enddo
      else if (ldoblha) then
        photon_hel=0
        h = 1
        init = 1
      else
        photon_hel=0
         if (bos.eq.2) then
            write(6,*) " "
            write(6,*) "Zjjj amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
         elseif (bos.eq.1) then
            write(6,*) " "
            write(6,*) "Ajjj amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
         elseif (bos.eq.6) then
            write(6,*) " "
            write(6,*) "Hjjj amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
         elseif (bos.eq.34 .or. bos.eq.43) then
            write(6,*) " "
            write(6,*) "W+W- jjj amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
            call vtoww_reset
         elseif (bos.eq.16 .or. bos.eq.61) then
            write(6,*) " "
            write(6,*) "HAjjj amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
         elseif (bos.eq.22) then
            write(6,*) " "
            write(6,*) "ZZ jjj -> 4l amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
            call vtozz_reset
         elseif (bos.eq.21) then
            write(6,*) " "
            write(6,*) "ZZ jjj -> 2l 2v amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
            call vtozz_reset
         elseif (bos.eq.211) then
            write(6,*) " "
            write(6,*) "ZA jjj -> 2l A amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
            call vtozz_reset
         elseif (bos.eq.212) then
            write(6,*) " "
            write(6,*) "ZA jjj -> 2v A amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
            call vtozz_reset
         elseif (bos.eq.66) then
            write(6,*) " "
            write(6,*) "HH jjj amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
            call vtohh_reset
         elseif (bos.eq.11) then
            write(6,*) " "
            write(6,*) "VV jjj (-> spin-2) -> gamma gamma  amplitude"
            write(6,*) "square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
            call vtosp2aa_reset
         elseif (bos.eq.734) then
           if (with_spin2) then  
            write(6,*) " "
            write(6,*) "VV jjj (-> Spin2) -> W+W- amplitude square info:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
           else
            write(6,*) " "
            write(6,*) "VV jjj -> Higgs (not Spin2!) -> W+W- amplitude square info:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
           endif
            call vtosp2ww_reset
         elseif (bos.eq.722) then
           if (with_spin2) then  
            write(6,*) " "
            write(6,*) "VV jjj (-> Spin2) -> ZZ amplitude square info:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
           else
            write(6,*) " "
            write(6,*) "VV jjj -> Higgs (not Spin2!) -> ZZ amplitude square info:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
           endif
            call vtosp2zz_reset
         elseif (bos.eq.721) then
           if (with_spin2) then  
            write(6,*) " "
            write(6,*) "VV jjj (-> Spin2) -> ZZ amplitude square info:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
           else
            write(6,*) " "
            write(6,*) "VV jjj -> Higgs (not Spin2!) -> ZZ amplitude square info:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
           endif
            call vtosp2zz_reset
         endif
         call printnfl(.true.)
         print*," "
         print*," creal = ",creal," cvirtual = ",cvirt
         print*," virt factor for alphas = 0.12 is ",
     1            1+0.12/pi*4./3.*cvirt
         h = 1
         init = 1
         do i = 1,maxnumsubproc
            do j = 0,3
               res(i,j) = 0
            enddo
         enddo
      endif   ! end of init
      
c
c -------------------------------------------------------------------
      call Calc_Momentum_Transfer(p, v, q12,q34,3)

c
c for WW->4l precalculate A->WW,Z->WW,AZ->WW etc leptonic tensors
c
      if (bos.eq.34 .or. bos.eq.43) then
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
c leptonic tensors
C for W+W-

         SELECT CASE(procid)
         CASE(WPhadWMjj, WPWMhadjj)

            do j = 2,3
               call anomal_formfactor(q12(0,j),q34(0,j),qp(0),qm(0))
               call atoww_had(v,aww)
               call ztoww_had(v,zww)
               call aatoww_had(q12(0,j),q34(0,j),v,aaww(0,0,j))
               call aztoww_had(q12(0,j),q34(0,j),v,azww(0,0,j))
               call aztoww_had(q34(0,j),q12(0,j),v,zaww(0,0,j))
               call zztoww_had(q12(0,j),q34(0,j),v,zzww(0,0,j))
               call wwtoww_had(q12(0,j),q34(0,j),v,wwww6(0,0,j)) ! q12 = W-
               call wwtoww_had(q34(0,j),q12(0,j),v,wwww5(0,0,j)) ! q12 = W+
C for WV --> e+ nu_e
c NCw tensors for NC process (k=1...4), CCw for CC (k=5,6)
               call WVtoWP_had(2,q34(0,j),v,NCwpa(0,0,1,j),NCwpz(0,0,1,j)) !emit W- on upper
               call WVtoWP_had(2,q12(0,j),v,NCwpa(0,0,2,j),NCwpz(0,0,2,j)) !emit W- on lower
               call WVtoWP_had(1,q34(0,j),v,CCwpa(0,0,1,j),CCwpz(0,0,1,j)) !emit W- on upper
               call WVtoWP_had(1,q12(0,j),v,CCwpa(0,0,2,j),CCwpz(0,0,2,j)) !emit W- on lower
C for WV --> mu- nu_mu
               call WVtoWM_had(2,q34(0,j),v,NCwma(0,0,1,j),NCwmz(0,0,1,j)) !emit W+ on upper
               call WVtoWM_had(2,q12(0,j),v,NCwma(0,0,2,j),NCwmz(0,0,2,j)) !emit W+ on lower
               call WVtoWM_had(1,q34(0,j),v,CCwma(0,0,1,j),CCwmz(0,0,1,j)) !emit W+ on upper
               call WVtoWM_had(1,q12(0,j),v,CCwma(0,0,2,j),CCwmz(0,0,2,j)) !emit W+ on lower
            enddo

         CASE DEFAULT
         
         if (with_anom) then    ! anomalous gauge boson couplings
c           using global form factor for all tensors of one phase space point
c           this ensures proper cancellations for anomalous contributions
c           energy scale is invariant WW mass

            do j = 2,3
              call anomal_formfactor(q12(0,j),q34(0,j),qp(0),qm(0))
              call atoww_anomal2(v,aww)
              call ztoww_anomal2(v,zww)

               call aatoww_anomal2(q12(0,j),q34(0,j),v,aaww(0,0,j))
               call aztoww_anomal2(q12(0,j),q34(0,j),v,azww(0,0,j))
               call aztoww_anomal2(q34(0,j),q12(0,j),v,zaww(0,0,j))
               call zztoww_anomal2(q12(0,j),q34(0,j),v,zzww(0,0,j))
               call wwtoww_anomal2(q12(0,j),q34(0,j),v,wwww6(0,0,j)) ! q12 = W-
               call wwtoww_anomal2(q34(0,j),q12(0,j),v,wwww5(0,0,j)) ! q12 = W+
	    	    
C for WV --> e+ nu_e
c NCw tensors for NC process (k=1...4), CCw for CC (k=5,6)
               call WVtoWP_anomal2(2,q34(0,j),v,NCwpa(0,0,1,j),
     -              NCwpz(0,0,1,j)) !emit W- on upper
               call WVtoWP_anomal2(2,q12(0,j),v,NCwpa(0,0,2,j),
     -              NCwpz(0,0,2,j)) !emit W- on lower
               call WVtoWP_anomal2(1,q34(0,j),v,CCwpa(0,0,1,j),
     -              CCwpz(0,0,1,j)) !emit W- on upper
               call WVtoWP_anomal2(1,q12(0,j),v,CCwpa(0,0,2,j),
     -              CCwpz(0,0,2,j)) !emit W- on lower
C for WV --> mu- nu_mu
               call WVtoWM_anomal2(2,q34(0,j),v,NCwma(0,0,1,j),
     -              NCwmz(0,0,1,j)) !emit W+ on upper
               call WVtoWM_anomal2(2,q12(0,j),v,NCwma(0,0,2,j),
     -              NCwmz(0,0,2,j)) !emit W+ on lower
               call WVtoWM_anomal2(1,q34(0,j),v,CCwma(0,0,1,j),
     -              CCwmz(0,0,1,j)) !emit W+ on upper
               call WVtoWM_anomal2(1,q12(0,j),v,CCwma(0,0,2,j),
     -              CCwmz(0,0,2,j)) !emit W+ on lower
            enddo

         elseif (with_kk) then  ! kaluza-klein scattering
#ifdef WITH_KK
            call atoww_kk(v,aww)
            call ztoww_kk(v,zww)
            do j = 2,3
               call aatoww_kk(q12(0,j),q34(0,j),v,aaww(0,0,j))
               call aztoww_kk(q12(0,j),q34(0,j),v,azww(0,0,j))
               call aztoww_kk(q34(0,j),q12(0,j),v,zaww(0,0,j))
               call zztoww_kk(q12(0,j),q34(0,j),v,zzww(0,0,j))
               call wwtoww_kk(q12(0,j),q34(0,j),v,wwww6(0,0,j)) ! q12 = W-
               call wwtoww_kk(q34(0,j),q12(0,j),v,wwww5(0,0,j)) ! q12 = W+
	    	    
C for WV --> e+ nu_e
c NCw tensors for NC process (k=1...4), CCw for CC (k=5,6)
               call WVtoWP_kk(2,q34(0,j),v,NCwpa(0,0,1,j),
     -              NCwpz(0,0,1,j)) !emit W- on upper
               call WVtoWP_kk(2,q12(0,j),v,NCwpa(0,0,2,j),
     -              NCwpz(0,0,2,j)) !emit W- on lower
               call WVtoWP_kk(1,q34(0,j),v,CCwpa(0,0,1,j),
     -              CCwpz(0,0,1,j)) !emit W- on upper
               call WVtoWP_kk(1,q12(0,j),v,CCwpa(0,0,2,j),
     -              CCwpz(0,0,2,j)) !emit W- on lower
C for WV --> mu- nu_mu
               call WVtoWM_kk(2,q34(0,j),v,NCwma(0,0,1,j),
     -              NCwmz(0,0,1,j)) !emit W+ on upper
               call WVtoWM_kk(2,q12(0,j),v,NCwma(0,0,2,j),
     -              NCwmz(0,0,2,j)) !emit W+ on lower
               call WVtoWM_kk(1,q34(0,j),v,CCwma(0,0,1,j),
     -              CCwmz(0,0,1,j)) !emit W+ on upper
               call WVtoWM_kk(1,q12(0,j),v,CCwma(0,0,2,j),
     -              CCwmz(0,0,2,j)) !emit W+ on lower
            enddo
#endif
 
         elseif (with_spin2) then 
            call atoww(v,aww)
            call ztoww(v,zww)

            do j = 2,3                  
               call aatoww_spin2(q12(0,j),q34(0,j),v,aaww(0,0,j))
               call aztoww_spin2(q12(0,j),q34(0,j),v,azww(0,0,j))
               call aztoww_spin2(q34(0,j),q12(0,j),v,zaww(0,0,j))
               call zztoww_spin2(q12(0,j),q34(0,j),v,zzww(0,0,j))
               call wwtoww_spin2(q12(0,j),q34(0,j),v,wwww6(0,0,j)) ! q12 = W-
               call wwtoww_spin2(q34(0,j),q12(0,j),v,wwww5(0,0,j)) ! q12 = W+
C for WV --> e+ nu_e
c NCw tensors for NC process (k=1...4), CCw for CC (k=5,6)
               call WVtoWP(2,q34(0,j),v,NCwpa(0,0,1,j),
     -              NCwpz(0,0,1,j)) !emit W- on upper
               call WVtoWP(2,q12(0,j),v,NCwpa(0,0,2,j),
     -              NCwpz(0,0,2,j)) !emit W- on lower
               call WVtoWP(1,q34(0,j),v,CCwpa(0,0,1,j),
     -              CCwpz(0,0,1,j)) !emit W- on upper
               call WVtoWP(1,q12(0,j),v,CCwpa(0,0,2,j),
     -              CCwpz(0,0,2,j)) !emit W- on lower
C     for WV --> mu- nu_mu
               call WVtoWM(2,q34(0,j),v,NCwma(0,0,1,j),
     -              NCwmz(0,0,1,j)) !emit W+ on upper
               call WVtoWM(2,q12(0,j),v,NCwma(0,0,2,j),
     -              NCwmz(0,0,2,j)) !emit W+ on lower
               call WVtoWM(1,q34(0,j),v,CCwma(0,0,1,j),
     -              CCwmz(0,0,1,j)) !emit W+ on upper
               call WVtoWM(1,q12(0,j),v,CCwma(0,0,2,j),
     -              CCwmz(0,0,2,j)) !emit W+ on lower
            enddo           
         else                   ! SM
            call atoww(v,aww)
            call ztoww(v,zww)
            do j = 2,3
               call aatoww(q12(0,j),q34(0,j),v,aaww(0,0,j))
               call aztoww(q12(0,j),q34(0,j),v,azww(0,0,j))
               call aztoww(q34(0,j),q12(0,j),v,zaww(0,0,j))
               call zztoww(q12(0,j),q34(0,j),v,zzww(0,0,j))
               call wwtoww(q12(0,j),q34(0,j),v,wwww6(0,0,j)) ! q12 = W-
               call wwtoww(q34(0,j),q12(0,j),v,wwww5(0,0,j)) ! q12 = W+
	    	    
C for WV --> e+ nu_e
c NCw tensors for NC process (k=1...4), CCw for CC (k=5,6)
               call WVtoWP(2,q34(0,j),v,NCwpa(0,0,1,j),
     -              NCwpz(0,0,1,j)) !emit W- on upper
               call WVtoWP(2,q12(0,j),v,NCwpa(0,0,2,j),
     -              NCwpz(0,0,2,j)) !emit W- on lower
               call WVtoWP(1,q34(0,j),v,CCwpa(0,0,1,j),
     -              CCwpz(0,0,1,j)) !emit W- on upper
               call WVtoWP(1,q12(0,j),v,CCwpa(0,0,2,j),
     -              CCwpz(0,0,2,j)) !emit W- on lower
C for WV --> mu- nu_mu
               call WVtoWM(2,q34(0,j),v,NCwma(0,0,1,j),
     -              NCwmz(0,0,1,j)) !emit W+ on upper
               call WVtoWM(2,q12(0,j),v,NCwma(0,0,2,j),
     -              NCwmz(0,0,2,j)) !emit W+ on lower
               call WVtoWM(1,q34(0,j),v,CCwma(0,0,1,j),
     -              CCwmz(0,0,1,j)) !emit W+ on upper
               call WVtoWM(1,q12(0,j),v,CCwma(0,0,2,j),
     -              CCwmz(0,0,2,j)) !emit W+ on lower
            enddo
         endif   ! end of kk/anom/sm

         end select
	 
c
c-----------------------------
c
c
c for ZZ->4l precalculate  leptonic tensors:

      elseif (bos.eq.22) then

c reset lzs to .true. to force recalculation of virtual contributions in qqzzqq
         do i = 1,4
            lzs(i) = .true.
         enddo
c lepton spinors and Z polarization vectors ! 
c	lepton helicities not fixed 
c 		-> sum over all possible helicites in |M|**2
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


         CALL IXXXXX(v(0,1),ZERO ,+ie,-1,lep) !e+
         CALL OXXXXX(v(0,2),ZERO ,-ie, 1,lem) !e- 
         CALL IXXXXX(v(0,3),ZERO ,+iu,-1,lup) !mu+
         CALL OXXXXX(v(0,4),ZERO ,-iu, 1,lum) !mu- 


         SELECT CASE(procid)
         CASE(ZZhadjj)
            ! for general output to all combinations: up-type first
            if ((finalquarks(1).eq.93 .or. finalquarks(1).eq.94)) then
               ! set couplings for Z hadronic decay into uubar
               call setZtouu
            endif
            CALL JIOXXX(lep,lem,GZ_ZF ,ZMASS,ZWIDTH,ze)     !Zl
            CALL JIOXXX(lep,lem,GZ_AF ,ZERO ,ZERO  ,ae)      !Al
            CALL JIOXXX(lup,lum,GZL ,ZMASS,ZWIDTH,zu)     !Zl
            CALL JIOXXX(lup,lum,GAL ,ZERO ,ZERO  ,au)      !Al

            do mu = 0,3
               qe(mu) = v(mu,1)+v(mu,2)
               qu(mu) = v(mu,3)+v(mu,4)
               qvv(mu) = qe(mu) + qu(mu)
            enddo

            qe(4) = qe(0)**2-qe(1)**2-qe(2)**2-qe(3)**2
            qu(4) = qu(0)**2-qu(1)**2-qu(2)**2-qu(3)**2

c leptonic tensors:

            do j = 2,3

               call anomal_formfactor(q12(0,j),q34(0,j),qe(0),qu(0))

               ! for V-> llll:	
               call vto4l_had(v,h,azz,zzztens)	
               ! for AA,AZ,ZA,ZZ -> e+ e- mu+ mu-:
               call vvtozz_had(q12(0,j),q34(0,j),v,h,aazz(0,0,j),azzz(0,0,j),zazz(0,0,j),zzzz(0,0,j))
               ! for W+W-
               call wwtozz_had(q12(0,j),q34(0,j),v,h,wwzz5(0,0,j)) ! q12 = W+
               call wwtozz_had(q34(0,j),q12(0,j),v,h,wwzz6(0,0,j)) ! q12 = W-
               ! for V1V2 --> e+ e- or  V1V2 -->  mu+ mu- (Vi=A,Z)
               ! NCl tensors for NC process (k=1...4)
               call VVtoll_had(2,1,h,q34(0,j),v,aaee(0,0,1,j),azee(0,0,1,j),zaee(0,0,1,j),zzee(0,0,1,j)) !emit V on upper
               call VVtoll_had(1,1,h,q12(0,j),v,aaee(0,0,2,j),azee(0,0,2,j),zaee(0,0,2,j),zzee(0,0,2,j)) !emit V on lower
               call VVtoll_had(2,2,h,q34(0,j),v,aauu(0,0,1,j),azuu(0,0,1,j),zauu(0,0,1,j),zzuu(0,0,1,j)) !emit V on upper
               call VVtoll_had(1,2,h,q12(0,j),v,aauu(0,0,2,j),azuu(0,0,2,j),zauu(0,0,2,j),zzuu(0,0,2,j)) !emit V on lower
               ! for W+W- --> e+ e- or  W+W- -->  mu+ mu-
               ! CCl tensors for CC (k=5)
               call WWtoll_had(2,1,h,q34(0,j),v,CCee(0,0,1,j)) !emit V on upper
               call WWtoll_had(1,1,h,q12(0,j),v,CCee(0,0,2,j)) !emit V on lower
               call WWtoll_had(2,2,h,q34(0,j),v,CCuu(0,0,1,j)) !emit V on upper
               call WWtoll_had(1,2,h,q12(0,j),v,CCuu(0,0,2,j)) !emit V on lower
               ! for W-W+ --> e+ e- or mu+ mu-
               ! CCl for CC (k=6)
               call WWtoll_had(1,1,h,q34(0,j),v,CCee6(0,0,1,j)) !emit V on upper
               call WWtoll_had(2,1,h,q12(0,j),v,CCee6(0,0,2,j)) !emit V on lower
               call WWtoll_had(1,2,h,q34(0,j),v,CCuu6(0,0,1,j)) !emit V on upper
               call WWtoll_had(2,2,h,q12(0,j),v,CCuu6(0,0,2,j)) !emit V on lower
               
            enddo               !j=2,3


         CASE DEFAULT

         
         CALL JIOXXX(lep,lem,GZL ,ZMASS,ZWIDTH,ze)     !Zl
         CALL JIOXXX(lep,lem,GAL ,ZERO ,ZERO  ,ae)      !Al
         CALL JIOXXX(lup,lum,GZL ,ZMASS,ZWIDTH,zu)     !Zl
         CALL JIOXXX(lup,lum,GAL ,ZERO ,ZERO  ,au)      !Al

         if (with_kk) then      ! Kaluza-Klein scenario
#ifdef WITH_KK
c for V-> llll:	
            call vto4l_KK(v,h,azz,zzztens)	
         

            do mu = 0,3
               qe(mu) = v(mu,1)+v(mu,2)
               qu(mu) = v(mu,3)+v(mu,4)
               qvv(mu) = qe(mu) + qu(mu)
            enddo
	 	 
            qe(4) = qe(0)**2-qe(1)**2-qe(2)**2-qe(3)**2
            qu(4) = qu(0)**2-qu(1)**2-qu(2)**2-qu(3)**2

c leptonic tensors:

            do j = 2,3
C for AA,AZ,ZA,ZZ -> e+ e- mu+ mu-:
               call vvtozz_KK(q12(0,j),q34(0,j),v,h,aazz(0,0,j),
     -              azzz(0,0,j),zazz(0,0,j),zzzz(0,0,j))
c
C for W+W-
               call wwtozz_KK(q12(0,j),q34(0,j),v,h,wwzz5(0,0,j)) ! q12 = W+
               call wwtozz_KK(q34(0,j),q12(0,j),v,h,wwzz6(0,0,j)) ! q12 = W-

C for V1V2 --> e+ e- or  V1V2 -->  mu+ mu- (Vi=A,Z)
c NCl tensors for NC process (k=1...4)

               call VVtoll_KK(2,1,h,q34(0,j),v,aaee(0,0,1,j),azee(0,0,1,j),
     -              zaee(0,0,1,j),zzee(0,0,1,j)) !emit V on upper        
               call VVtoll_KK(1,1,h,q12(0,j),v,aaee(0,0,2,j),azee(0,0,2,j),
     -              zaee(0,0,2,j),zzee(0,0,2,j)) !emit V on lower
               call VVtoll_KK(2,2,h,q34(0,j),v,aauu(0,0,1,j),azuu(0,0,1,j),
     -              zauu(0,0,1,j),zzuu(0,0,1,j)) !emit V on upper        
               call VVtoll_KK(1,2,h,q12(0,j),v,aauu(0,0,2,j),azuu(0,0,2,j),
     -              zauu(0,0,2,j),zzuu(0,0,2,j)) !emit V on lower

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
               
	 
            enddo               !j=2,3
#endif
 
         elseif (with_spin2) then   
c for V-> llll:	
            call vto4l(v,h,azz,zzztens)	
         

            do mu = 0,3
               qe(mu) = v(mu,1)+v(mu,2)
               qu(mu) = v(mu,3)+v(mu,4)
               qvv(mu) = qe(mu) + qu(mu)
            enddo
	 	 
          qe(4) = qe(0)**2-qe(1)**2-qe(2)**2-qe(3)**2
          qu(4) = qu(0)**2-qu(1)**2-qu(2)**2-qu(3)**2

c leptonic tensors:

            do j = 2,3
C for AA,AZ,ZA,ZZ -> e+ e- mu+ mu-:
               call vvtozz_spin2(q12(0,j),q34(0,j),v,h,aazz(0,0,j),
     -              azzz(0,0,j),zazz(0,0,j),zzzz(0,0,j))
c
C for W+W-
               call wwtozz_spin2(q12(0,j),q34(0,j),v,h,wwzz5(0,0,j))  ! q12 = W+
               call wwtozz_spin2(q34(0,j),q12(0,j),v,h,wwzz6(0,0,j))  ! q12 = W-

c other processes do not change

C for V1V2 --> e+ e- or  V1V2 -->  mu+ mu- (Vi=A,Z)
c NCl tensors for NC process (k=1...4)

               call VVtoll(2,1,h,q34(0,j),v,aaee(0,0,1,j),
     -              azee(0,0,1,j),zaee(0,0,1,j),zzee(0,0,1,j))  !emit V on upper        
	       call VVtoll(1,1,h,q12(0,j),v,aaee(0,0,2,j),
     -              azee(0,0,2,j),zaee(0,0,2,j),zzee(0,0,2,j))  !emit V on lower
               call VVtoll(2,2,h,q34(0,j),v,aauu(0,0,1,j),
     -              azuu(0,0,1,j),zauu(0,0,1,j),zzuu(0,0,1,j))  !emit V on upper        
	       call VVtoll(1,2,h,q12(0,j),v,aauu(0,0,2,j),
     -              azuu(0,0,2,j),zauu(0,0,2,j),zzuu(0,0,2,j))  !emit V on lower

C for W+W- --> e+ e- or  W+W- -->  mu+ mu-
c  CCl tensors for CC (k=5)

               call WWtoll(2,1,h,q34(0,j),v,CCee(0,0,1,j))     !emit V on upper
               call WWtoll(1,1,h,q12(0,j),v,CCee(0,0,2,j))     !emit V on lower
               call WWtoll(2,2,h,q34(0,j),v,CCuu(0,0,1,j))     !emit V on upper
               call WWtoll(1,2,h,q12(0,j),v,CCuu(0,0,2,j))     !emit V on lower
	 
C for W-W+ --> e+ e- or mu+ mu-
c CCl for CC (k=6)
	
               call WWtoll(1,1,h,q34(0,j),v,CCee6(0,0,1,j))     !emit V on upper
               call WWtoll(2,1,h,q12(0,j),v,CCee6(0,0,2,j))     !emit V on lower
               call WWtoll(1,2,h,q34(0,j),v,CCuu6(0,0,1,j))     !emit V on upper
               call WWtoll(2,2,h,q12(0,j),v,CCuu6(0,0,2,j))     !emit V on lower
	 
	 
        enddo !j=2,3

         elseif (with_anom) then ! anomalous gauge boson couplings
c           using global form factor for all tensors of one phase space point
c           this ensures proper cancellations for anomalous contributions
c           energy scale is invariant ZZ mass        

            do mu = 0,3
               qe(mu) = v(mu,1)+v(mu,2)
               qu(mu) = v(mu,3)+v(mu,4)
               qvv(mu) = qe(mu) + qu(mu)
            enddo
	 	 
            qe(4) = qe(0)**2-qe(1)**2-qe(2)**2-qe(3)**2
            qu(4) = qu(0)**2-qu(1)**2-qu(2)**2-qu(3)**2

c leptonic tensors:

            do j = 2,3

               call anomal_formfactor(q12(0,j),q34(0,j),qe(0),qu(0))
c for V-> llll:	
               call vto4l_anomal(v,h,azz,zzztens)	
         

C for AA,AZ,ZA,ZZ -> e+ e- mu+ mu-:
               call vvtozz_anomal(q12(0,j),q34(0,j),v,h,aazz(0,0,j),azzz(0,0,j),
     -              zazz(0,0,j),zzzz(0,0,j))
c
C for W+W-
               call wwtozz_anomal(q12(0,j),q34(0,j),v,h,wwzz5(0,0,j)) ! q12 = W+
               call wwtozz_anomal(q34(0,j),q12(0,j),v,h,wwzz6(0,0,j)) ! q12 = W-

C for V1V2 --> e+ e- or  V1V2 -->  mu+ mu- (Vi=A,Z)
c NCl tensors for NC process (k=1...4)

               call VVtoll_anomal(2,1,h,q34(0,j),v,aaee(0,0,1,j),azee(0,0,1,j),
     -              zaee(0,0,1,j),zzee(0,0,1,j)) !emit V on upper        
               call VVtoll_anomal(1,1,h,q12(0,j),v,aaee(0,0,2,j),azee(0,0,2,j),
     -              zaee(0,0,2,j),zzee(0,0,2,j)) !emit V on lower
               call VVtoll_anomal(2,2,h,q34(0,j),v,aauu(0,0,1,j),azuu(0,0,1,j),
     -              zauu(0,0,1,j),zzuu(0,0,1,j)) !emit V on upper        
               call VVtoll_anomal(1,2,h,q12(0,j),v,aauu(0,0,2,j),azuu(0,0,2,j),
     -              zauu(0,0,2,j),zzuu(0,0,2,j)) !emit V on lower

C for W+W- --> e+ e- or  W+W- -->  mu+ mu-
c  CCl tensors for CC (k=5)

               call WWtoll_anomal(2,1,h,q34(0,j),v,CCee(0,0,1,j)) !emit V on upper
               call WWtoll_anomal(1,1,h,q12(0,j),v,CCee(0,0,2,j)) !emit V on lower
               call WWtoll_anomal(2,2,h,q34(0,j),v,CCuu(0,0,1,j)) !emit V on upper
               call WWtoll_anomal(1,2,h,q12(0,j),v,CCuu(0,0,2,j)) !emit V on lower
	 
C for W-W+ --> e+ e- or mu+ mu-
c CCl for CC (k=6)
	
               call WWtoll_anomal(1,1,h,q34(0,j),v,CCee6(0,0,1,j)) !emit V on upper
               call WWtoll_anomal(2,1,h,q12(0,j),v,CCee6(0,0,2,j)) !emit V on lower
               call WWtoll_anomal(1,2,h,q34(0,j),v,CCuu6(0,0,1,j)) !emit V on upper
               call WWtoll_anomal(2,2,h,q12(0,j),v,CCuu6(0,0,2,j)) !emit V on lower
               
	 
            enddo               !j=2,3


         else                   !SM
c for V-> llll:	
            call vto4l(v,h,azz,zzztens)	
         

            do mu = 0,3
               qe(mu) = v(mu,1)+v(mu,2)
               qu(mu) = v(mu,3)+v(mu,4)
               qvv(mu) = qe(mu) + qu(mu)
            enddo
	 	 
            qe(4) = qe(0)**2-qe(1)**2-qe(2)**2-qe(3)**2
            qu(4) = qu(0)**2-qu(1)**2-qu(2)**2-qu(3)**2

c leptonic tensors:

            do j = 2,3
C for AA,AZ,ZA,ZZ -> e+ e- mu+ mu-:
               call vvtozz(q12(0,j),q34(0,j),v,h,aazz(0,0,j),azzz(0,0,j),
     -              zazz(0,0,j),zzzz(0,0,j))
c
C for W+W-
               call wwtozz(q12(0,j),q34(0,j),v,h,wwzz5(0,0,j)) ! q12 = W+
               call wwtozz(q34(0,j),q12(0,j),v,h,wwzz6(0,0,j)) ! q12 = W-

C for V1V2 --> e+ e- or  V1V2 -->  mu+ mu- (Vi=A,Z)
c NCl tensors for NC process (k=1...4)

               call VVtoll(2,1,h,q34(0,j),v,aaee(0,0,1,j),azee(0,0,1,j),
     -              zaee(0,0,1,j),zzee(0,0,1,j)) !emit V on upper        
               call VVtoll(1,1,h,q12(0,j),v,aaee(0,0,2,j),azee(0,0,2,j),
     -              zaee(0,0,2,j),zzee(0,0,2,j)) !emit V on lower
               call VVtoll(2,2,h,q34(0,j),v,aauu(0,0,1,j),azuu(0,0,1,j),
     -              zauu(0,0,1,j),zzuu(0,0,1,j)) !emit V on upper        
               call VVtoll(1,2,h,q12(0,j),v,aauu(0,0,2,j),azuu(0,0,2,j),
     -              zauu(0,0,2,j),zzuu(0,0,2,j)) !emit V on lower

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
               
	 
            enddo               !j=2,3
         endif   ! end of kk/sm for bos=22


         end select  ! semileptonic decays

c-----------------------------
c

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

      CALL VXXXXX(v(0,1),ZERO ,lsign(3),1,sp2au(1))          !A1
      CALL VXXXXX(v(0,2),ZERO ,lsign(4),1,sp2ae(1))          !A2

c leptonic tensors
         do j = 2,3
!       CALL wwsp2toaa(q34(0,j),q12(0,j),v,lsign,sp2wwaa5(0,0,j))
!       CALL wwsp2toaa(q12(0,j),q34(0,j),v,lsign,sp2wwaa6(0,0,j))
!       CALL zzsp2toaa(q34(0,j),q12(0,j),v,lsign,sp2zzaa(0,0,j))
!       CALL azsp2toaa(q34(0,j),q12(0,j),v,lsign,sp2azaa(0,0,j))
!       CALL azsp2toaa(q12(0,j),q34(0,j),v,lsign,sp2zaaa(0,0,j))
!       CALL aasp2toaa(q34(0,j),q12(0,j),v,lsign,sp2aaaa(0,0,j))

       call vvsp2tovv(1,q34(0,j),q12(0,j),v(0,1),v(0,2),sp2au(1),sp2ae(1), 
     1           sp2wwaa5(0,0,j),sp2wwaa6(0,0,j),sp2zzaa(0,0,j),
     2           sp2azaa(0,0,j),sp2zaaa(0,0,j),sp2aaaa(0,0,j))
         enddo


c for VV -> spin-2 -> WW -> 4l precalculate leptonic tensors
       elseif (bos.eq.734) then !  734 = spin-2 -> W+W-
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

c leptonic tensors
C for spin-2 -> W+W-
          do j = 2,3               
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
          enddo


c for VV -> spin-2 -> ZZ -> 4l precalculate leptonic tensors
      elseif (bos.eq.722) then   !  722 = spin-2 -> ZZ

c reset lzs to .true. to force recalculation of virtual contributions in qqzzqq
         do i = 1,4
            lzs(i) = .true.
         enddo
c lepton spinors and Z polarization vectors ! 
c	lepton helicities not fixed 
c 		-> sum over all possible helicites in |M|**2
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
		
		
         CALL IXXXXX(v(0,1),ZERO ,+ie,-1,lep) !e+	
         CALL OXXXXX(v(0,2),ZERO ,-ie, 1,lem) !e- 
         CALL IXXXXX(v(0,3),ZERO ,+iu,-1,lup) !mu+	
         CALL OXXXXX(v(0,4),ZERO ,-iu, 1,lum) !mu- 
         
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

c leptonic tensors:

            do j = 2,3

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

      else          ! Higgs instead of spin-2 resonance     
    
       call vvtozz_Hres(q12(0,j),q34(0,j),v,h,sp2aazz(0,0,j),sp2azzz(0,0,j),
     -              sp2zazz(0,0,j),sp2zzzz(0,0,j))

       call wwtozz_Hres(q12(0,j),q34(0,j),v,h,sp2wwzz5(0,0,j)) ! q12 = W+
       call wwtozz_Hres(q34(0,j),q12(0,j),v,h,sp2wwzz6(0,0,j)) ! q12 = W-

      endif      

	    enddo !j=2,3


c for VV -> spin-2 -> ZZ -> 2l 2v precalculate  leptonic tensors
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
c
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
         do i = 5,6   !add momentum info to A->vv current (au vanishes anyway) 
            au(i) = zu(i)
         enddo    	
	
        
            do mu = 0,3
               qe(mu) = v(mu,1)+v(mu,2)
               qu(mu) = v(mu,3)+v(mu,4)
               qvv(mu) = qe(mu) + qu(mu)
            enddo
	 	 
            qe(4) = qe(0)**2-qe(1)**2-qe(2)**2-qe(3)**2
            qu(4) = qu(0)**2-qu(1)**2-qu(2)**2-qu(3)**2

c leptonic tensors
          do j = 2,3
 
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

        else          ! Higgs instead of spin-2 resonance     

         call vvtozzn_Hres(q12(0,j),q34(0,j),v,h,sp2aazz(0,0,j),sp2azzz(0,0,j),
     -        sp2zazz(0,0,j),sp2zzzz(0,0,j))

         call wwtozzn_Hres(q12(0,j),q34(0,j),v,h,sp2wwzz5(0,0,j)) ! q12 = W+
         call wwtozzn_Hres(q34(0,j),q12(0,j),v,h,sp2wwzz6(0,0,j)) ! q12 = W-
	     
         endif             
           enddo               !j = 2,3	

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
c
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
         do i = 5,6   !add momentum info to A->vv current (au vanishes anyway) 
            au(i) = zu(i)
         enddo    	
	
c for V-> llvv
         if (with_kk) then      ! Kaluza-Klein Scenario 
#ifdef WITH_KK
            call vto4ln_KK(v,h,azz,zzztens)	
        
            do mu = 0,3
               qe(mu) = v(mu,1)+v(mu,2)
               qu(mu) = v(mu,3)+v(mu,4)
               qvv(mu) = qe(mu) + qu(mu)
            enddo
	 	 
            qe(4) = qe(0)**2-qe(1)**2-qe(2)**2-qe(3)**2
            qu(4) = qu(0)**2-qu(1)**2-qu(2)**2-qu(3)**2
            
c leptonic tensors
            do j = 2,3
C     for AA,AZ,ZA,ZZ -> e+ e- vm vm~:
               call vvtozzn_KK(q12(0,j),q34(0,j),v,h,aazz(0,0,j),
     -              azzz(0,0,j),zazz(0,0,j),zzzz(0,0,j))
c
C for W+W-
               call wwtozzn_KK(q12(0,j),q34(0,j),v,h,wwzz5(0,0,j)) ! q12 = W+
               call wwtozzn_KK(q34(0,j),q12(0,j),v,h,wwzz6(0,0,j)) ! q12 = W-
	 
C for V1V2 --> e+ e- or  V1V2 -->  vm vm~ (Vi=A,Z)
c NCl tensors for NC process (k=1...4)

               call VVtoll_KK(2,1,h,q34(0,j),v,aaee(0,0,1,j),azee(0,0,1,j),
     -              zaee(0,0,1,j),zzee(0,0,1,j)) !emit V on upper        
               call VVtoll_KK(1,1,h,q12(0,j),v,aaee(0,0,2,j),azee(0,0,2,j),
     -              zaee(0,0,2,j),zzee(0,0,2,j)) !emit V on lower
               call VVtolln_KK(2,q34(0,j),v(0,3),aauu(0,0,1,j),azuu(0,0,1,j),
     -              zauu(0,0,1,j),zzuu(0,0,1,j)) !emit V on upper        
               call VVtolln_KK(1,q12(0,j),v(0,3),aauu(0,0,2,j),azuu(0,0,2,j),
     -              zauu(0,0,2,j),zzuu(0,0,2,j)) !emit V on lower


C for W+W- --> e+ e- or  W+W- -->  vm vm~
c  CCl tensors for CC (k=5)

               call WWtoll_KK(2,1,h,q34(0,j),v,CCee(0,0,1,j)) !emit V on upper
               call WWtoll_KK(1,1,h,q12(0,j),v,CCee(0,0,2,j)) !emit V on lower
               call WWtolln_KK(2,q34(0,j),v,CCuu(0,0,1,j)) !emit V on upper
               call WWtolln_KK(1,q12(0,j),v,CCuu(0,0,2,j)) !emit V on lower
	 
C for W-W+ --> e+ e- or vm vm~
c CCl for CC (k=6)
	
               call WWtoll_KK(1,1,h,q34(0,j),v,CCee6(0,0,1,j)) !emit V on upper
               call WWtoll_KK(2,1,h,q12(0,j),v,CCee6(0,0,2,j)) !emit V on lower
               call WWtolln_KK(1,q34(0,j),v,CCuu6(0,0,1,j)) !emit V on upper
               call WWtolln_KK(2,q12(0,j),v,CCuu6(0,0,2,j)) !emit V on lower
	 
            enddo               !j	
#endif
 

         elseif (with_spin2) then   
            call vto4ln(v,h,azz,zzztens)	
        
            do mu = 0,3
               qe(mu) = v(mu,1)+v(mu,2)
               qu(mu) = v(mu,3)+v(mu,4)
               qvv(mu) = qe(mu) + qu(mu)
            enddo
	 	 
            qe(4) = qe(0)**2-qe(1)**2-qe(2)**2-qe(3)**2
            qu(4) = qu(0)**2-qu(1)**2-qu(2)**2-qu(3)**2

c leptonic tensors
            do j = 2,3
C for AA,AZ,ZA,ZZ -> e+ e- vm vm~:
               call vvtozzn_spin2(q12(0,j),q34(0,j),v,h,aazz(0,0,j),
     -              azzz(0,0,j),zazz(0,0,j),zzzz(0,0,j))
c
C for W+W-
               call wwtozzn_spin2(q12(0,j),q34(0,j),v,h,wwzz5(0,0,j))  ! q12 = W+
               call wwtozzn_spin2(q34(0,j),q12(0,j),v,h,wwzz6(0,0,j))  ! q12 = W-

c other processes do not change
	 
C for V1V2 --> e+ e- or  V1V2 -->  vm vm~ (Vi=A,Z)
c NCl tensors for NC process (k=1...4)

               call VVtoll(2,1,h,q34(0,j),v,aaee(0,0,1,j),
     -              azee(0,0,1,j),zaee(0,0,1,j),zzee(0,0,1,j))   !emit V on upper        
	       call VVtoll(1,1,h,q12(0,j),v,aaee(0,0,2,j),
     -              azee(0,0,2,j),zaee(0,0,2,j),zzee(0,0,2,j))   !emit V on lower
               call VVtolln(2,q34(0,j),v(0,3),aauu(0,0,1,j),
     -              azuu(0,0,1,j),zauu(0,0,1,j),zzuu(0,0,1,j))   !emit V on upper        
	       call VVtolln(1,q12(0,j),v(0,3),aauu(0,0,2,j),
     -              azuu(0,0,2,j),zauu(0,0,2,j),zzuu(0,0,2,j))   !emit V on lower


C for W+W- --> e+ e- or  W+W- -->  vm vm~
c  CCl tensors for CC (k=5)

               call WWtoll(2,1,h,q34(0,j),v,CCee(0,0,1,j))     !emit V on upper
               call WWtoll(1,1,h,q12(0,j),v,CCee(0,0,2,j))     !emit V on lower
               call WWtolln(2,q34(0,j),v,CCuu(0,0,1,j))     !emit V on upper
               call WWtolln(1,q12(0,j),v,CCuu(0,0,2,j))     !emit V on lower
	 
C for W-W+ --> e+ e- or vm vm~
c CCl for CC (k=6)
	
               call WWtoll(1,1,h,q34(0,j),v,CCee6(0,0,1,j))     !emit V on upper
               call WWtoll(2,1,h,q12(0,j),v,CCee6(0,0,2,j))     !emit V on lower
               call WWtolln(1,q34(0,j),v,CCuu6(0,0,1,j))     !emit V on upper
               call WWtolln(2,q12(0,j),v,CCuu6(0,0,2,j))     !emit V on lower
	 
        enddo !j	


         elseif (with_anom) then ! anomalous gauge boson couplings
c        using global form factor for all tensors of one phase space point
c        this ensures proper cancellations for anomalous contributions
c        energy scale is invariant ZZ mass        
  
            do mu = 0,3
               qe(mu) = v(mu,1)+v(mu,2)
               qu(mu) = v(mu,3)+v(mu,4)
               qvv(mu) = qe(mu) + qu(mu)
            enddo
	 	 
            qe(4) = qe(0)**2-qe(1)**2-qe(2)**2-qe(3)**2
            qu(4) = qu(0)**2-qu(1)**2-qu(2)**2-qu(3)**2

c leptonic tensors
            do j = 2,3

               call anomal_formfactor(q12(0,j),q34(0,j),qe(0),qu(0))
               call vto4ln_anomal(v,h,azz,zzztens)	

C for AA,AZ,ZA,ZZ -> e+ e- vm vm~:
               call vvtozzn_anomal(q12(0,j),q34(0,j),v,h,aazz(0,0,j),azzz(0,0,j),
     -              zazz(0,0,j),zzzz(0,0,j))
c
C for W+W-
               call wwtozzn_anomal(q12(0,j),q34(0,j),v,h,wwzz5(0,0,j)) ! q12 = W+
               call wwtozzn_anomal(q34(0,j),q12(0,j),v,h,wwzz6(0,0,j)) ! q12 = W-
	 
C for V1V2 --> e+ e- or  V1V2 -->  vm vm~ (Vi=A,Z)
c NCl tensors for NC process (k=1...4)
               
               call VVtoll_anomal(2,1,h,q34(0,j),v,aaee(0,0,1,j),azee(0,0,1,j),
     -              zaee(0,0,1,j),zzee(0,0,1,j)) !emit V on upper        
               call VVtoll_anomal(1,1,h,q12(0,j),v,aaee(0,0,2,j),azee(0,0,2,j),
     -              zaee(0,0,2,j),zzee(0,0,2,j)) !emit V on lower
               call VVtolln_anomal(2,q34(0,j),v(0,3),aauu(0,0,1,j),azuu(0,0,1,j),
     -              zauu(0,0,1,j),zzuu(0,0,1,j)) !emit V on upper        
               call VVtolln_anomal(1,q12(0,j),v(0,3),aauu(0,0,2,j),azuu(0,0,2,j),
     -              zauu(0,0,2,j),zzuu(0,0,2,j)) !emit V on lower


C for W+W- --> e+ e- or  W+W- -->  vm vm~
c  CCl tensors for CC (k=5)

               call WWtoll_anomal(2,1,h,q34(0,j),v,CCee(0,0,1,j)) !emit V on upper
               call WWtoll_anomal(1,1,h,q12(0,j),v,CCee(0,0,2,j)) !emit V on lower
               call WWtolln_anomal(2,q34(0,j),v,CCuu(0,0,1,j)) !emit V on upper
               call WWtolln_anomal(1,q12(0,j),v,CCuu(0,0,2,j)) !emit V on lower
	 
C for W-W+ --> e+ e- or vm vm~
c CCl for CC (k=6)
	
               call WWtoll_anomal(1,1,h,q34(0,j),v,CCee6(0,0,1,j)) !emit V on upper
               call WWtoll_anomal(2,1,h,q12(0,j),v,CCee6(0,0,2,j)) !emit V on lower
               call WWtolln_anomal(1,q34(0,j),v,CCuu6(0,0,1,j)) !emit V on upper
               call WWtolln_anomal(2,q12(0,j),v,CCuu6(0,0,2,j)) !emit V on lower
               
            enddo               !j	


         else                   !SM
            call vto4ln(v,h,azz,zzztens)	
        
            do mu = 0,3
               qe(mu) = v(mu,1)+v(mu,2)
               qu(mu) = v(mu,3)+v(mu,4)
               qvv(mu) = qe(mu) + qu(mu)

            enddo
	 	 
            qe(4) = qe(0)**2-qe(1)**2-qe(2)**2-qe(3)**2
            qu(4) = qu(0)**2-qu(1)**2-qu(2)**2-qu(3)**2

c leptonic tensors
            do j = 2,3
C for AA,AZ,ZA,ZZ -> e+ e- vm vm~:
               call vvtozzn(q12(0,j),q34(0,j),v,h,aazz(0,0,j),azzz(0,0,j),
     -              zazz(0,0,j),zzzz(0,0,j))
c
C for W+W-
               call wwtozzn(q12(0,j),q34(0,j),v,h,wwzz5(0,0,j)) ! q12 = W+
               call wwtozzn(q34(0,j),q12(0,j),v,h,wwzz6(0,0,j)) ! q12 = W-
	 
C for V1V2 --> e+ e- or  V1V2 -->  vm vm~ (Vi=A,Z)
c NCl tensors for NC process (k=1...4)
               
               call VVtoll(2,1,h,q34(0,j),v,aaee(0,0,1,j),azee(0,0,1,j),
     -              zaee(0,0,1,j),zzee(0,0,1,j)) !emit V on upper        
               call VVtoll(1,1,h,q12(0,j),v,aaee(0,0,2,j),azee(0,0,2,j),
     -              zaee(0,0,2,j),zzee(0,0,2,j)) !emit V on lower
               call VVtolln(2,q34(0,j),v(0,3),aauu(0,0,1,j),azuu(0,0,1,j),
     -              zauu(0,0,1,j),zzuu(0,0,1,j)) !emit V on upper        
               call VVtolln(1,q12(0,j),v(0,3),aauu(0,0,2,j),azuu(0,0,2,j),
     -              zauu(0,0,2,j),zzuu(0,0,2,j)) !emit V on lower


C for W+W- --> e+ e- or  W+W- -->  vm vm~
c  CCl tensors for CC (k=5)

               call WWtoll(2,1,h,q34(0,j),v,CCee(0,0,1,j)) !emit V on upper
               call WWtoll(1,1,h,q12(0,j),v,CCee(0,0,2,j)) !emit V on lower
               call WWtolln(2,q34(0,j),v,CCuu(0,0,1,j)) !emit V on upper
               call WWtolln(1,q12(0,j),v,CCuu(0,0,2,j)) !emit V on lower
	 
C for W-W+ --> e+ e- or vm vm~
c CCl for CC (k=6)

               call WWtoll(1,1,h,q34(0,j),v,CCee6(0,0,1,j)) !emit V on upper
               call WWtoll(2,1,h,q12(0,j),v,CCee6(0,0,2,j)) !emit V on lower
               call WWtolln(1,q34(0,j),v,CCuu6(0,0,1,j)) !emit V on upper
               call WWtolln(2,q12(0,j),v,CCuu6(0,0,2,j)) !emit V on lower
               
            enddo               !j	
         endif    ! end of kk/sm for bos=21
         
         
      elseif (bos.eq.211) then   !ZA -> l+l- A
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
	  ie = (-1)**(h)
	  iu = photon_hel
	endif


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
            zu(i)=  ea(i)
         enddo    
         
         do mu = 0,3
            qe(mu) = v(mu,1)+v(mu,2)
            qa(mu) = v(mu,3)
            qvv(mu) = qe(mu) + qa(mu)
         enddo
         
         qe(4) = qe(0)**2-qe(1)**2-qe(2)**2-qe(3)**2
         qa(4) = abs(qa(0)**2-qa(1)**2-qa(2)**2-qa(3)**2)+1d-20

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

         do j=2,3
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
         enddo ! j

        else ! SM
         do j=2,3
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
         enddo ! j
      endif ! SM
      
C     for V --> e+ e- a
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
	if(photon_hel.eq.0) then
         h = mod(h,2) + 1
         ie = 1
         iu = sign(1,1-h)
        else
	 ie=1
	 iu=photon_hel
	endif
     

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

       do j=2,3
        qu=qa
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
        enddo ! j
        
      else ! SM
      
       do j=2,3
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
        enddo ! j
      endif ! SM
      
      azz=0d0
      zzztens=0d0
         qu=qa
         
c-----------------------------

      elseif (bos.eq.16.or.bos.eq.61) then
         if(lok(1))CALL Get_HAjjj_currents(p,v)
         if(lok(2))CALL Get_HAjj_subtr_currents(p,v,2)
         if(lok(3))CALL Get_HAjj_subtr_currents(p,v,3)
c --------------------------------------------
	 
      elseif (bos.eq.66) then   ! HHjjj
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
               vtemp(mu,1)=v(mu,1)+v(mu,2)
               vtemp(mu,2)=v(mu,3)+v(mu,4)
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
         do j=2,3

c-----------------------------
c
C for W+W-
            call wwtohh(q12(0,j),q34(0,j),v,wwhh5(0,0,j)) ! q12 = W+
            call wwtohh(q34(0,j),q12(0,j),v,wwhh6(0,0,j)) ! q12 = W-
C for ZZ
            call zztohh(q12(0,j),q34(0,j),v,zzhh(0,0,j)) ! q12 = Z

         enddo
c     -----------------------------------------------------------------------

      endif   ! end of bos = ...

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
c  ---------------------------------------------------------------------
c
cc  debugging:
cc
cc collinear subtraction:
c	if (lcoll) then 
c      	    ldebug = lok(1).and.lok(2).and.lok(3).and.
c     #      	    dotrr(p(0,5,1),p(0,1,1)).lt.1
c	    if (ldebug) print*,"coll.dot.:",lok(1),lok(2),lok(3),
c     #      		       dotrr(p(0,5,1),p(0,1,1))
cc      
cc soft subtraction:	
c       elseif (lsoft) then		       
c	   ldebug = lok(1).and.lok(2).and.lok(3).and.p(0,5,1).lt.1	
c	   if (ldebug) print*,"soft.mom.:",lok(1),lok(2),lok(3),p(0,5,1)   
c       endif
c
c  ---------------------------------------------------------------------
c
c scales and als:

	if (.false.) then
         write(6,*) " "
         write(6,*) "xi(1), xi(2) =", xi
         write(6,*) "mu_f L=1:",sqrt(mufsq(1,1)),sqrt(mufsq(2,1))
         write(6,*) "mu_f L=2:",sqrt(mufsq(1,2)),sqrt(mufsq(2,2))
         write(6,*) "mu_f L=3:",sqrt(mufsq(1,3)),sqrt(mufsq(2,3))
         write(6,*) " "
         write(6,*) "xi(1), xi(2) =", xi
         write(6,*) "alphas L=1:",als(1,1), als(2,1)
         write(6,*) "alphas L=2:",als(1,2), als(2,2)
         write(6,*) "alphas L=3:",als(1,3), als(2,3)
        end if
      if (lwarn) then    ! check connection between scales
         if(abs(mufsq(1,2)/mufsq(1,1)-1).gt.1d-8) then
            print*," check muf 1,2 vs 1,1 ",mufsq(1,2),mufsq(1,1)
         endif
         if(abs(mufsq(2,3)/mufsq(2,1)-1).gt.1d-8) then
            print*," check muf 2,3 vs 2,1 ",mufsq(2,3),mufsq(2,1)
         endif
      endif

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
            write(6,10) i, p(0,i,1), p(1,i,1), p(2,i,1), p(3,i,1)
         end do
c
         write(6,*) " "
         write(6,*) "xi(1), xi(2) =", xi
         write(6,*) "mu_f =",sqrt(mufsq(1,1)),sqrt(mufsq(2,1))
         read(*,*)

       end if


c
c  end debugging
c

       if (ldoblha) then

         physToDiag(1)=blha_physdiagmap(1,blha_cursubproc,blha_curproc)
         physToDiag(2)=blha_physdiagmap(2,blha_cursubproc,blha_curproc)
         physToDiag(3)=blha_physdiagmap(3,blha_cursubproc,blha_curproc)
         physToDiag(4)=blha_physdiagmap(4,blha_cursubproc,blha_curproc)
         physToDiag(5)=blha_physdiagmap(5,blha_cursubproc,blha_curproc)
         fsign(1) = blha_fsign(1,blha_cursubproc,blha_curproc)
         fsign(2) = blha_fsign(2,blha_cursubproc,blha_curproc)
         fsign(3) = blha_fsign(3,blha_cursubproc,blha_curproc)
         fsign(4) = blha_fsign(4,blha_cursubproc,blha_curproc)
         gsign    = blha_fsign(5,blha_cursubproc,blha_curproc)

         call wbf_zh3j(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1            pdf,res,nmin,nmax)

         return
       endif


c  ---------------------------------------------------------------------
c
c call PDF subroutine in order to determine parton
c distributions in the incoming (anti)protons.
c for the NLO contributions x1 = x*y in my notes with x=z=xuz(1,j), y=xi(j-1)
c
      pdf=0d0
      x1 = xi(1)*xuz(1,2)
      q_sf = sqrt(mufsq(1,1))
*      write(*,*)'1 x1, q_sf =', x1, q_sf
      call pdfproton( xi(1), q_sf, pdf(-6,1,1) )      ! f_a(y)=f_a(x1/z)

      if (nlo.gt.0) then
      q_sf = sqrt(mufsq(1,2))                         ! f_a(x1) for upper line
*      write(*,*)'2 q_sf =', q_sf 
      call pdfproton( x1, q_sf, pdf(-6,1,2) )         !   NLO correction
      endif

      if (mufsq(1,3).ne.mufsq(1,1)) then              ! f_a(x1) for lower line
         q_sf = sqrt(mufsq(1,3))                      !   NLO correction
*         write(*,*)'3 q_sf =', q_sf 
         call pdfproton( xi(1), q_sf, pdf(-6,1,3) )
      else
         do i = -6,6
            pdf(i,1,3) = pdf(i,1,1)
         enddo
      endif

      x2 = xi(2)*xuz(1,3)
      q_sf = sqrt(mufsq(2,1))
*      write(*,*)'4 x2, q_sf =', x2, q_sf
      call pdfproton( xi(2), q_sf, pdf(-6,2,1) )      ! f_b(y)=f_a(x2/z)

      if (nlo.gt.0) then
      q_sf = sqrt(mufsq(2,3))                         ! f_b(x2) for lower line 
*      write(*,*)'5 q_sf =', q_sf 
      call pdfproton( x2, q_sf, pdf(-6,2,3) )         !   NLO correction
      endif

      if (mufsq(2,2).ne.mufsq(2,1)) then              ! f_b(x2) for upper line
         q_sf = sqrt(mufsq(2,2))                      !   NLO correction
*         write(*,*)'6 q_sf =', q_sf 
         call pdfproton( xi(2), q_sf, pdf(-6,2,2) )
      else
         do i = -6,6
            pdf(i,2,2) = pdf(i,2,1)
         enddo
      endif


c$$$      do i = -6, 6
c$$$         write(*,*)'i, pdf(i,1,1-3) =', i, pdf(i,1,1), pdf(i,1,2), pdf(i,1,3)
c$$$         write(*,*)'i, pdf(i,2,1-3) =', i, pdf(i,2,1), pdf(i,2,2), pdf(i,2,3)
c$$$      end do
c$$$      stop

c
c and fill the coefficient functions for the finite subtraction terms
C Note that the color factors TR and C2 are NOT!!! included here
      lnQomu(2) = log(q12(4,2)/mufsq(1,2))
      lnQomu(3) = log(q34(4,3)/mufsq(2,3))
      omxi(2) = 1d0-x1
      omxi(3) = 1d0-x2
      tgs2oqsq(2) = 8d0*pi*als(1,2)/q12(4,2)         !2g_s^2/Q^2
      tgs2oqsq(3) = 8d0*pi*als(2,3)/q34(4,3)         !2g_s^2/Q^2

      do j = 2,3
         z = xuz(1,j)
         ln1mxi = log(omxi(j))
c         tgs2oqsq(j) = tgs2oqsq(j) * 6*xuz(2,j)*(1d0-xuz(2,j))
C cut off Int_x1^1 dz log(1-z) at z<1-5E-6. This assures that the relative
C error of Int_x1^1 dz log(1-z) is less than 1E-4/(1-x1). see notes p 28.4
         if (z.lt.0.999995) then
            lnrat = lnQomu(j) + log((1d0-z)/z)
            lnz = log(z)
            Ax(j) = (z**2+(1-z)**2)*lnrat + 2*z*(1-z)
            Bx(j) = ( 2d0*(lnrat+lnz)-1.5d0 )/(1d0-z)
            Cx(j) = 1-z - 2d0*lnz/(1d0-z) - (1+z)*lnrat
         else
            Ax(j) = 0d0
            Bx(j) = 0d0
            Cx(j) = 0d0
         endif
         Dxi(j) = 1.5d0*lnQomu(j) + 2d0*ln1mxi*lnQomu(j)
     1            + ln1mxi**2 - 1.5d0*ln1mxi 
     2            + creal
         do i = 1,5
            pdf(-i,j-1,j)=(pdf(-i,j-1,j)*(Dxi(j)/omxi(j)-Bx(j)) +
     1                     pdf(-i,j-1,1)*(Bx(j)+Cx(j)))*tgs2oqsq(j)
            pdf(i,j-1,j)= (pdf(i,j-1,j)*(Dxi(j)/omxi(j)-Bx(j)) +
     1                     pdf(i,j-1,1)*(Bx(j)+Cx(j)))*tgs2oqsq(j)
         enddo
         pdf(0,j-1,j)=pdf(0,j-1,1)*Ax(j)*tgs2oqsq(j)
      enddo
c select helicity
      if (bos.eq.2) then

         jsig = min(8*rn(1)+1d0,8.01d0)
c         jsig = 0
      else
         jsig = 0
      endif

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     Now sum over the subprocesses contributing to H production

      nmax = FL_ZHg(fsign,-1)            !reset counter for subprocesses to 0

C*******************  q1 q3 ---> q2 q4 g V V   **********************
         
c   physToDiag(ext.momentum label) = Feynman diagram label

      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=3
      physToDiag(3)=2    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=4
      physToDiag(5)=5   ! gluon
C NOTE: for call of wbf_zh3j it is important that p(*,1,*) and p(*,3,*)
c correspond to 1-2 fermion line ALWAYS, i.e physToDiag(1/2)={1,3} and 
c similarly physToDiag(3/4)={2,4} for the 3-4 fermion line
      fsign(1) = 1
      fsign(2) = 1
      fsign(3) = 1
      fsign(4) = 1
      gsign    = 1

      call wbf_zh3j(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     1  	    pdf,res,nmin,nmax)

*      write(*,*)'1, res =', res(

c$$$      do sj = 1, maxnumsubproc
c$$$         write(*,*)'res(sj,0-3) =', res(sj,0),res(sj,1),res(sj,2),res(sj,3)
c$$$      end do
c$$$      write(*,*)'res(1,0-3) =', res(1,0),res(1,1),res(1,2),res(1,3)
c$$$      write(*,*)'res(2,0-3) =', res(2,0),res(2,1),res(2,2),res(2,3)

c      if (ldebug) call resprint(nmin,nmax,res)
      if (init.eq.1) nproc(1) = nmax

C*******************  q1 qb4 ---> q2 qb3 g V V   **********************

      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=4
      physToDiag(3)=2    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=3
      physToDiag(5)=5   ! gluon

      fsign(3) = -1
      fsign(4) = -1

      call wbf_zh3j(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     1  	   pdf,res,nmin,nmax)

c      if (ldebug) call resprint(nmin,nmax,res)
      if (init.eq.1) nproc(2) = nmax

C*******************  qbar2 q3 ---> qbar1 q4 g V V   **********************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=3
      physToDiag(3)=1    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=4
      physToDiag(5)=5   ! gluon

      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = 1
      fsign(4) = 1
      
      call wbf_zh3j(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     1  	   pdf,res,nmin,nmax)

c      if (ldebug) call resprint(nmin,nmax,res)
      if (init.eq.1) nproc(3) = nmax

C*******************  qbar2 qb4 ---> qbar1 qb3 g V V  ******************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=4
      physToDiag(3)=1    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=3
      physToDiag(5)=5   ! gluon

      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = -1 
      fsign(4) = -1
         
      call wbf_zh3j(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1  	   pdf,res,nmin,nmax)

c      if (ldebug) call resprint(nmin,nmax,res)
      if (init.eq.1) nproc(4) = nmax

C*******************  g q3 ---> qb1 q2 q4 V V  **********************

      physToDiag(1)=5
      physToDiag(2)=3
      physToDiag(3)=1
      physToDiag(4)=4
      physToDiag(5)=2

      fsign(1) = -1
      fsign(2) =  1
      fsign(3) =  1
      fsign(4) =  1
      gsign    = -1

      call wbf_zh3j(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1  	   pdf,res,nmin,nmax)

c      if (ldebug) call resprint(nmin,nmax,res)
      if (init.eq.1) nproc(5) = nmax

C*******************  q1 g ---> q2 qb3 q4 V V   **********************

      physToDiag(1)=1
      physToDiag(2)=5
      physToDiag(3)=2
      physToDiag(4)=3
      physToDiag(5)=4

      fsign(1) =  1
      fsign(2) =  1
      fsign(3) = -1
      fsign(4) =  1
      gsign    = -1

      call wbf_zh3j(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1  	   pdf,res,nmin,nmax)

c      if (ldebug) call resprint(nmin,nmax,res)
      if (init.eq.1) nproc(6) = nmax

C*******************  g qbar4 ---> qbar3 qb1 q2 V V   **********************

      physToDiag(1)=5
      physToDiag(2)=4
      physToDiag(3)=1
      physToDiag(4)=3
      physToDiag(5)=2

      fsign(1) = -1
      fsign(2) =  1
      fsign(3) = -1
      fsign(4) = -1

      call wbf_zh3j(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1  	   pdf,res,nmin,nmax)

c      if (ldebug) call resprint(nmin,nmax,res)
      if (init.eq.1) nproc(7) = nmax

C*******************  qbar2 g ---> qbar1 qb3 q4 V V  **********************

      physToDiag(1)=2
      physToDiag(2)=5
      physToDiag(3)=1
      physToDiag(4)=3
      physToDiag(5)=4
      
      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = -1
      fsign(4) =  1

      call wbf_zh3j(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     1  	   pdf,res,nmin,nmax)

c      if (ldebug) call resprint(nmin,nmax,res)
      if (init.eq.1) nproc(8) = nmax


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

            CALL JIOXXX(lep,lem,GZ_ZF ,ZMASS,ZWIDTH,ze)     !Zl
            CALL JIOXXX(lep,lem,GZ_AF ,ZERO ,ZERO  ,ae)      !Al

c leptonic tensors:

            do j = 2,3

               call anomal_formfactor(q12(0,j),q34(0,j),qe(0),qu(0))

               ! for V-> llll:	
               call vto4l_had(v,h,azz,zzztens)	
               ! for AA,AZ,ZA,ZZ -> e+ e- mu+ mu-:
               call vvtozz_had(q12(0,j),q34(0,j),v,h,aazz(0,0,j),azzz(0,0,j),zazz(0,0,j),zzzz(0,0,j))
               ! for W+W-
               call wwtozz_had(q12(0,j),q34(0,j),v,h,wwzz5(0,0,j)) ! q12 = W+
               call wwtozz_had(q34(0,j),q12(0,j),v,h,wwzz6(0,0,j)) ! q12 = W-
               ! for V1V2 --> e+ e- or  V1V2 -->  mu+ mu- (Vi=A,Z)
               ! NCl tensors for NC process (k=1...4)
               call VVtoll_had(2,1,h,q34(0,j),v,aaee(0,0,1,j),azee(0,0,1,j),zaee(0,0,1,j),zzee(0,0,1,j)) !emit V on upper
               call VVtoll_had(1,1,h,q12(0,j),v,aaee(0,0,2,j),azee(0,0,2,j),zaee(0,0,2,j),zzee(0,0,2,j)) !emit V on lower
               call VVtoll_had(2,2,h,q34(0,j),v,aauu(0,0,1,j),azuu(0,0,1,j),zauu(0,0,1,j),zzuu(0,0,1,j)) !emit V on upper
               call VVtoll_had(1,2,h,q12(0,j),v,aauu(0,0,2,j),azuu(0,0,2,j),zauu(0,0,2,j),zzuu(0,0,2,j)) !emit V on lower
               ! for W+W- --> e+ e- or  W+W- -->  mu+ mu-
               ! CCl tensors for CC (k=5)
               call WWtoll_had(2,1,h,q34(0,j),v,CCee(0,0,1,j)) !emit V on upper
               call WWtoll_had(1,1,h,q12(0,j),v,CCee(0,0,2,j)) !emit V on lower
               call WWtoll_had(2,2,h,q34(0,j),v,CCuu(0,0,1,j)) !emit V on upper
               call WWtoll_had(1,2,h,q12(0,j),v,CCuu(0,0,2,j)) !emit V on lower
               ! for W-W+ --> e+ e- or mu+ mu-
               ! CCl for CC (k=6)
               call WWtoll_had(1,1,h,q34(0,j),v,CCee6(0,0,1,j)) !emit V on upper
               call WWtoll_had(2,1,h,q12(0,j),v,CCee6(0,0,2,j)) !emit V on lower
               call WWtoll_had(1,2,h,q34(0,j),v,CCuu6(0,0,1,j)) !emit V on upper
               call WWtoll_had(2,2,h,q12(0,j),v,CCuu6(0,0,2,j)) !emit V on lower
               
            enddo               !j=2,3


C*******************  q1 q3 ---> q2 q4 g V V   **********************
         
c   physToDiag(ext.momentum label) = Feynman diagram label

      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=3
      physToDiag(3)=2    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=4
      physToDiag(5)=5   ! gluon
C NOTE: for call of wbf_zh3j it is important that p(*,1,*) and p(*,3,*)
c correspond to 1-2 fermion line ALWAYS, i.e physToDiag(1/2)={1,3} and 
c similarly physToDiag(3/4)={2,4} for the 3-4 fermion line
      fsign(1) = 1
      fsign(2) = 1
      fsign(3) = 1
      fsign(4) = 1
      gsign    = 1

      call wbf_zh3j(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1  	    pdf,res,nmin,nmax)

*      write(*,*)'1, res =', res(

c$$$      do sj = 1, maxnumsubproc
c$$$         write(*,*)'res(sj,0-3) =', res(sj,0),res(sj,1),res(sj,2),res(sj,3)
c$$$      end do
c$$$      write(*,*)'res(1,0-3) =', res(1,0),res(1,1),res(1,2),res(1,3)
c$$$      write(*,*)'res(2,0-3) =', res(2,0),res(2,1),res(2,2),res(2,3)

c      if (ldebug) call resprint(nmin,nmax,res)
      if (init.eq.1) nproc(1) = nmax

C*******************  q1 qb4 ---> q2 qb3 g V V   **********************

      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=4
      physToDiag(3)=2    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=3
      physToDiag(5)=5   ! gluon

      fsign(3) = -1
      fsign(4) = -1

      call wbf_zh3j(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1  	   pdf,res,nmin,nmax)

c      if (ldebug) call resprint(nmin,nmax,res)
      if (init.eq.1) nproc(2) = nmax

C*******************  qbar2 q3 ---> qbar1 q4 g V V   **********************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=3
      physToDiag(3)=1    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=4
      physToDiag(5)=5   ! gluon

      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = 1
      fsign(4) = 1
      
      call wbf_zh3j(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     1  	   pdf,res,nmin,nmax)

c      if (ldebug) call resprint(nmin,nmax,res)
      if (init.eq.1) nproc(3) = nmax

C*******************  qbar2 qb4 ---> qbar1 qb3 g V V  ******************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=4
      physToDiag(3)=1    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=3
      physToDiag(5)=5   ! gluon

      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = -1 
      fsign(4) = -1
         
      call wbf_zh3j(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1  	   pdf,res,nmin,nmax)

c      if (ldebug) call resprint(nmin,nmax,res)
      if (init.eq.1) nproc(4) = nmax

C*******************  g q3 ---> qb1 q2 q4 V V  **********************

      physToDiag(1)=5
      physToDiag(2)=3
      physToDiag(3)=1
      physToDiag(4)=4
      physToDiag(5)=2

      fsign(1) = -1
      fsign(2) =  1
      fsign(3) =  1
      fsign(4) =  1
      gsign    = -1

      call wbf_zh3j(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1  	   pdf,res,nmin,nmax)

c      if (ldebug) call resprint(nmin,nmax,res)
      if (init.eq.1) nproc(5) = nmax

C*******************  q1 g ---> q2 qb3 q4 V V   **********************

      physToDiag(1)=1
      physToDiag(2)=5
      physToDiag(3)=2
      physToDiag(4)=3
      physToDiag(5)=4

      fsign(1) =  1
      fsign(2) =  1
      fsign(3) = -1
      fsign(4) =  1
      gsign    = -1

      call wbf_zh3j(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1  	   pdf,res,nmin,nmax)

c      if (ldebug) call resprint(nmin,nmax,res)
      if (init.eq.1) nproc(6) = nmax

C*******************  g qbar4 ---> qbar3 qb1 q2 V V   **********************

      physToDiag(1)=5
      physToDiag(2)=4
      physToDiag(3)=1
      physToDiag(4)=3
      physToDiag(5)=2

      fsign(1) = -1
      fsign(2) =  1
      fsign(3) = -1
      fsign(4) = -1

      call wbf_zh3j(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1  	   pdf,res,nmin,nmax)

c      if (ldebug) call resprint(nmin,nmax,res)
      if (init.eq.1) nproc(7) = nmax

C*******************  qbar2 g ---> qbar1 qb3 q4 V V  **********************

      physToDiag(1)=2
      physToDiag(2)=5
      physToDiag(3)=1
      physToDiag(4)=3
      physToDiag(5)=4
      
      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = -1
      fsign(4) =  1

      call wbf_zh3j(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     1  	   pdf,res,nmin,nmax)

c      if (ldebug) call resprint(nmin,nmax,res)
      if (init.eq.1) nproc(8) = nmax


c***  end of extra process evaluation for hadronic decay  ***

       endif
      end select


c*****************  end of process evaluation  *********************


      SELECT CASE(procid)
      CASE(WPhadWMjj)
       call hadDecayFactor_W(finalquarks(1), abs(mjj2(v(0,1),v(0,2))), N_gen_W, fac_W)
       do j=0,3
        DO IRES = 1,NMAX
          res(IRES,j) = res(IRES,j) * fac_W            ! factor fac_W for hadronic decay
        ENDDO
       enddo
      CASE(WPWMhadjj)
       call hadDecayFactor_W(finalquarks(1), abs(mjj2(v(0,3),v(0,4))), N_gen_W, fac_W)
       do j=0,3
        DO IRES = 1,NMAX
          res(IRES,j) = res(IRES,j) * fac_W            ! factor fac_W for hadronic decay
        ENDDO
       enddo
      CASE(ZZhadjj)
       call hadDecayFactor_Z(finalquarks(1), abs(mjj2(v(0,1),v(0,2))), N_gen_up, N_gen_down, fac_Z_up, fac_Z_down)
       do j=0,3
        DO IRES = 1,NMAX
          if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
             if (ires.le.(NMAX/2)) then
                res(IRES,j) = res(IRES,j) * fac_Z_up      ! factor fac_Z_up for up-type with all combinations
             elseif (ires.gt.(NMAX/2)) then
                res(IRES,j) = res(IRES,j) * fac_Z_down    ! factor fac_Z_down for down-type with all combinations
             endif
             if (mod(nmax,2).ne.0) then
               print*, "Something went wrong with the hadronic decays!"
               stop
             endif
          elseif (mod(abs(finalquarks(1)),2).eq.0) then
             res(IRES,j) = res(IRES,j) * fac_Z_up         ! up-type final states
          else
             res(IRES,j) = res(IRES,j) * fac_Z_down       ! down-type final states
          endif
        ENDDO
       enddo
      END SELECT


      if (init.eq.1) then
         init = init+1
         if (lwarn) print 199," proc #s for Z/Hjjj are ",nproc
 199     format(a,8i5)
      endif

      do j=0,3
         m2s(j) = 0      
         DO IRES = 1,NMAX
            m2s(j) = m2s(j) + RES(IRES,j)
*            write(*,*)'j,ires,res =', j,ires,res(ires,j)
         ENDDO

         if(j.eq.0)then
c...Les Houches interface - the most propable subprocess 3jets at LO  
            if ((lha.or.hepmc).and..not.doNLO) then
               i=0
               weight=0.d0
               rnumb=RandomNumber()
               do while((i.le.nmax).and.(weight.le.rnumb*m2s(0)))
                  i=i+1
                  weight=weight+res(i,0)
                  iprocess=i
               enddo
               SELECT CASE(procid)
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
         endif

c         DO IRES = nproc(4)+1,nproc(8)
c            m2s(j) = m2s(j) + RES(IRES,j)
c         ENDDO


         if (bos.eq.2) then
            if (jsig.eq.0) then
               m2s(j) = m2s(j)  !*2     ! factor 2 for electrons and muons
            else
               m2s(j) = m2s(j)*8d0   ! factor 8 for random helicity summation
            endif
         elseif ((bos.eq.22).or.(bos.eq.722).or.(bos.eq.211)) then ! qq->ZZqq (both Z->ll~)
               m2s(j) = m2s(j)*4d0   ! factor 4 for random helicity summation
         elseif (bos.eq.11) then ! qq (-> spin-2)-> gamma gamma qq
               m2s(j) = m2s(j)*4d0   ! factor 4 for random helicity summation
         elseif ((bos.eq.21).or.(bos.eq.721).or.(bos.eq.212)) then ! qq->ZZqq (one Z->ll~)
               m2s(j) = m2s(j)*2d0   ! factor 4 for random helicity summation
         endif
	if(photon_hel.ne.0) m2s(j)=m2s(j)/2d0
      enddo  ! j loop


      if ((lha.or.hepmc).and..not.doNLO) then
         if (jsig.ne.0) then
            if (bos.eq.2) then
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


*      stop
c
c ---------------------------------------------------------------
c
c debugging:
c
      if (ldebug) then
         do j=1,3
            if (lok(j) .and. m2s(j).eq.0 ) then
               print*,j," lok(j) = ",lok(j)," m2s = ",m2s(j)
            endif
         enddo
         if (m2s(0).ne.0) then
c            if (abs((m2s(1)+m2s(2)+m2s(3))/m2s(0)-1).gt.1d-5) then
               print*," m2s(0) = ",m2s(0),(m2s(1)+m2s(2)+m2s(3))/m2s(0)
               print*," ratio subtr/real = ",(m2s(2)+m2s(3))/m2s(1)
               print*," m2s(0,1,2,3) = ",m2s
	       print*
c            endif
         else
            print*," m2s = ",m2s
         endif
         if (nmax.ne.nmaxold) print*," nmax: ",nmaxold,nmax
         nmaxold = nmax
      endif
c      if (ldebug) print*," lokt = ",lok
c      if (ldebug) read(*,113) text
c 113  format(a)

      RETURN
      END


c******************************************************************
c
c   begin subroutine wbf_zh3j
c
c*****************************************************************
      subroutine wbf_zh3j(xi,p,v,physToDiag,fsign,gsign,bos,
     1                    nlo,lok,xuz,
     1                    pdf,res,nmin,nmax)
          use globalvars, only: ldoblha
      implicit none
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/process.inc"

      real*8 p(0:3,max_p,max_kin), v(0:3,max_v), xi(nx), xuz(2,2:3)
      real*8 pdf(-6:6,2,3), res(maxnumsubproc,0:3)
      integer physToDiag(5), fsign(4+max_v), 
     1        gsign, h, bos, nlo, nmin, nmax
      logical lok(3)
c      
      common /hcount / h
c
c wbf_zh3j calls the amplitude square routines 
c             qqzqqj    for qq-->qqZ g   for bos=2      and
c             qqhqq     for qq-->qqH     for bos=6
c             qqwwqq     for qq-->qqWW     for bos=34.or.43
c             qqzzqq     for qq-->qqZZ     for bos=22.or.21
c for the subtraction terms for the NLO cross section calculation
c
c  INPUT:  p(0:3,5,3)      external physical parton momenta
c          v(0:3,nv)       Z decay momenta, for Higgs production only
c                          the sum q(mu) = v(mu,1)+...+v(mu,nv) is needed
c          physToDiag(5)   physToDiag(ext.mom. label) = Feynman diagram label
c          fsign,gsign     sign factors for amplitude calls; see qqZqq(j)
c	   h		   specify lepton helicity combination
c          nlo             nlo = 0: calculate LO only (i.e. no subtraction)
c                          nlo = 1: calculate full NLO subtraction
c          lok(3)          lok(ID)=T means momenta set ID passes acceptance
c                          cuts, i.e res(k,ID) needs to be calculated
c  OUTPUT:
c          uucc(ID)   output in format as in qqZqq(j), but corresponding
c          etc.            to res(*,ID) in m2s_qqZqq
c                          ID = 1  : the real emission |M|^2 * pdf
c                             = 2,3: sutraction terms for emision off 12 or 
c     or   res(k,ID)                 34 line
c                             = 0  : subtracted result which drives 
c                                    integration, i.e 
c                                    res(*,0) = res(*,1)+res(*,2)+res(*,3)
c
c  In and Output
c    nmin, nmax            range of process numbers for this call
c
      real*8 C2, TR, N                              ! color factors
      parameter(N=3d0,TR=0.5d0,C2=TR*(N**2-1d0)/N)

c Note: Factor 9 from color sum included in qqbqq. I am using color summed
c amplitudes here, while CS use color averages. Thus I later divide by
c 8*3 for an initial qg state and by 3*3 for an initial qq state. So
c
c  2*gs2*C2*9/(8*3) = 2 gs^2 * TR *8/3*9/(8*3) = 2 gs^2 * TR
c
c consistent with Catani 
c
c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"
c
c  helicity selection
      integer jsig, jsig1, jsig3, jsig5
      common /chelsum/ jsig,jsig1,jsig3,jsig5
c
      real*8 uuccb(-1:1,-1:1,2:3), ddccb(-1:1,-1:1,2:3),
     1       uussb(-1:1,-1:1,2:3), ddssb(-1:1,-1:1,2:3),
     2       udscb(-1:1,-1:1,2:3), ducsb(-1:1,-1:1,2:3),
     3       uucc(0:3), uuss(0:3), ddcc(0:3),
     4       ddss(0:3), udsc(0:3), ducs(0:3)

      real*8 pbar(0:3,4+max_v), qbar(0:4),q2,v2,q3,v3,sub(2:3), 
     1       e_in(2), dotrr, gs2(2:3), polcolq, polcolg
      real*8 xa
      double precision NCmatrixelt(0:1,0:1,3,2)
      double precision CCmatrixelt(0:1,3,2)
      integer iflav(5), diagToPhys(5), FL_ZHg
      external dotrr, FL_ZHg

      logical ChargedCurrent, sametype, oneAntiparticle
      logical loldp, ldebug, lres
      data e_in/0d0,0d0/
      save e_in,loldp, uuccb,uussb,ddccb,ddssb,udscb,ducsb,
     1     gs2,polcolq,polcolg
      integer i,if1,if2,icc1,j,k,mu, fc(6), fcb(6)
      double precision temp2, temp3, temp1
      

      parameter (ldebug=.false.)
      PARAMETER (lres=.false.)

      double precision tree(6)


      if (lres) then
        open(unit=31,file="irfinite_VBF_VV.chk",ACCESS='APPEND')
      endif


      do j = 0,3
         uucc(j) = 0d0
         uuss(j) = 0d0
         ddcc(j) = 0d0
         ddss(j) = 0d0
         udsc(j) = 0d0
         ducs(j) = 0d0
      enddo


      nmin = nmax+1
      do i = 1,5
         diagToPhys(physToDiag(i)) = i
      enddo

      loldp = e_in(1).eq.p(0,3,1) .and. e_in(2).eq.p(0,4,1)
c  reset the LO amplitude to 0 to avoid wrong subtraction for initial gluon
      if (.not.loldp) then            ! this is a new phase space point
         e_in(1) = p(0,3,1)
         e_in(2) = p(0,4,1)
         do i = -1,1
            do j = -1,1
               do k = 2,3
                  uuccb(i,j,k) = 0d0
                  uussb(i,j,k) = 0d0
                  ddccb(i,j,k) = 0d0
                  ddssb(i,j,k) = 0d0
                  udscb(i,j,k) = 0d0
                  ducsb(i,j,k) = 0d0
               enddo
            enddo
         enddo
c  determine strong coupling gs for the two quark lines and factor for 
c  polarization and spin average
         gs2(2) = 4d0*pi*als(1,1)
         gs2(3) = 4d0*pi*als(2,1)
         polcolq = 1d0/(4d0*N**2*xi(1)*xi(2))
         polcolg = 1d0/(4d0*N*(N**2-1)*xi(1)*xi(2))
      endif
         
c get the real emission amplitude squared, store it in uucc(-1,1) etc.
      do mu = 0,3
         do i = 1,5
            pbar(mu,physToDiag(i))=p(mu,i,1)
         enddo
         do i = 6,4+max_v
            pbar(mu,i)=0d0
         enddo
         qbar(mu) = pbar(mu,5)
      enddo
      qbar(4)=0d0

      if (bos.eq.2) then  ! Z
         do mu = 0,3             ! kinematics for Z-->l+l- decay
            pbar(mu,5) = v(mu,1) ! ebar
            pbar(mu,6) = v(mu,2) ! e
         enddo
         fsign(5) = -1
         fsign(6) = 1
         if (lok(1) .or. nlo.eq.0) 
     1      call qqzqqj_c(pbar,fsign,qbar,gsign, 
     1            uucc(1),uuss(1),ddcc(1),ddss(1),udsc(1),ducs(1))

      elseif (bos.eq.1) then  ! photon
         do mu = 0,3
            pbar(mu,5) = v(mu,n_v) ! photon
            pbar(mu,6) = 0d0 
         enddo
         fsign(5) = -1
         fsign(6) = 0
         if (lok(1) .or. nlo.eq.0) then
            call qqAqqj_c(pbar, fsign, qbar, gsign,
     1           uucc(1),uuss(1),ddcc(1),ddss(1),udsc(1),ducs(1),tree)
         endif
      elseif (bos.eq.6) then  ! higgs
         do mu = 0,3
            pbar(mu,5) = 0      ! dummy momentum
            pbar(mu,6) = 0
            do i = 1,n_v
               pbar(mu,6) = pbar(mu,6) + v(mu,i) ! Higgs momentum
            enddo
         enddo
         fsign(5) = 0
         fsign(6) = 1         
         if (lok(1) .or. nlo.eq.0) 
     1      call qqhqqj_c(pbar,fsign,qbar,gsign, 
     1            uucc(1),uuss(1),ddcc(1),ddss(1),udsc(1),ducs(1))
      elseif (bos.eq.16.or.bos.eq.61) then  ! Higgs plus photon
         fsign(5) = 0
         fsign(6) = 1
         if (lok(1) .or. nlo.eq.0) then
            call HAjjj_ME(fsign,gsign,
     $           uucc(1),uuss(1),ddcc(1),ddss(1),udsc(1),ducs(1))
         endif
      elseif (bos.eq.34 .or. bos.eq.43) then         ! W+W- to 4 leptons
         do mu = 0,3            ! kinematics for H-->WW -->4 lepton decay
            pbar(mu,5) = v(mu,2) ! l+
            pbar(mu,6) = v(mu,1) ! nu
            pbar(mu,7) = v(mu,3) ! nubar
            pbar(mu,8) = v(mu,4) ! l-
         enddo
         fsign(5) = -1
         fsign(6) = 1
         fsign(7) = -1
         fsign(8) = 1
c
         if (lok(1) .or. nlo.eq.0) then
c original version for H->WW only:		 
c           call qqhqqj_c(pbar,fsign,qbar,gsign, 
c     1            uucc(1),uuss(1),ddcc(1),ddss(1),udsc(1),ducs(1))
c	
c new routine "bqqwwqqj" with madgraph-call:
c            call bqqwwqqj_comp(pbar,fsign,qbar,gsign,
c     &           uucc(1),uuss(1),ddcc(1),ddss(1),udsc(1),ducs(1))
     
c new routine (with speed-up): uucc(1) stands for uucc(1:3) here!    
            call qqwwqqj(pbar,fsign,qbar,gsign,
     &            uucc(1),uuss(1),ddcc(1),ddss(1),udsc(1),ducs(1))  
         endif	

      elseif (bos.eq.22) then         ! ZZ to 4 leptons
         do mu = 0,3            ! kinematics for H-->ZZ -->4 lepton decay
            pbar(mu,5) = v(mu,1) ! l+
            pbar(mu,6) = v(mu,2) ! l-
            pbar(mu,7) = v(mu,3) ! l"+
            pbar(mu,8) = v(mu,4) ! l"-
         enddo
         fsign(5) = -1
         fsign(6) = 1
         fsign(7) = -1
         fsign(8) = 1
c
         if (lok(1) .or. nlo.eq.0) then
c routine "bqqzzqqj" with madgraph-call:
c            call bqqzzqqj_comp(pbar,fsign,qbar,gsign,bos,
c     &           uucc(1),uuss(1),ddcc(1),ddss(1),udsc(1),ducs(1))
c         
c new routine (with speed-up): uucc(1) stands for uucc(1:3) here!      
            call qqzzqqj(pbar,fsign,qbar,gsign,bos,
     &            uucc(1),uuss(1),ddcc(1),ddss(1),udsc(1),ducs(1)) 
	
c	compare to madgraph:
c            call bqqzzqqj_comp(pbar,fsign,qbar,gsign,bos,
c     &           uucc(1),uuss(1),ddcc(1),ddss(1),udsc(1),ducs(1))
     
        endif	!lok

      elseif (bos.eq.66) then ! HH production
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

         if (lok(1) .or. nlo.eq.0) then
            call qqhhqqj(pbar,fsign,qbar,gsign,bos,
     &           uucc(1),uuss(1),ddcc(1),ddss(1),udsc(1),ducs(1))
         endif !lok
         
      elseif (bos.eq.11) then         ! qq(-->spin-2)-->qq gamma gamma
         do mu = 0,3             
            pbar(mu,5) = v(mu,1) ! photon 1
            pbar(mu,6) = v(mu,2) ! photon 2
         enddo

         fsign(5) = 1
         fsign(6) = 1

            call qqsp2aaqqj(pbar,fsign,qbar,gsign,bos,
     &            uucc(1),uuss(1),ddcc(1),ddss(1),udsc(1),ducs(1))         
        

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
c
         if (lok(1) .or. nlo.eq.0) then

            call qqsp2wwqqj(pbar,fsign,qbar,gsign,
     &            uucc(1),uuss(1),ddcc(1),ddss(1),udsc(1),ducs(1))  
         endif	 


      elseif (bos.eq.722) then   ! spin-2 --> ZZ- --> 4 leptons
         do mu = 0,3             ! kinematics for spin-2 ->ZZ -->4 lepton decay
            pbar(mu,5) = v(mu,1) ! l+
            pbar(mu,6) = v(mu,2) ! l-
            pbar(mu,7) = v(mu,3) ! l"+
            pbar(mu,8) = v(mu,4) ! l"-
         enddo
         fsign(5) = -1
         fsign(6) = 1
         fsign(7) = -1
         fsign(8) = 1
c
         if (lok(1) .or. nlo.eq.0) then   
            call qqsp2zzqqj(pbar,fsign,qbar,gsign,bos,
     &            uucc(1),uuss(1),ddcc(1),ddss(1),udsc(1),ducs(1))      
        endif	!lok

      elseif (bos.eq.721) then         ! spin-2 --> ZZ to 2 leptons 2 neutrinos
         do mu = 0,3            ! kinematics for spin-2 -->ZZ -->2l 2v lepton decay
            pbar(mu,5) = v(mu,1) ! l+
            pbar(mu,6) = v(mu,2) ! l-
            pbar(mu,7) = v(mu,3) ! v
            pbar(mu,8) = v(mu,4) ! v~
         enddo
         fsign(5) = -1
         fsign(6) =  1
         fsign(7) =  1
         fsign(8) = -1
c
         if (lok(1) .or. nlo.eq.0) then
            call qqsp2zzqqj(pbar,fsign,qbar,gsign,bos,
     &            uucc(1),uuss(1),ddcc(1),ddss(1),udsc(1),ducs(1)) 
         endif  !lok	

      elseif (bos.eq.21) then         ! ZZ to 2 leptons 2 neutrinos
         do mu = 0,3            ! kinematics for H-->ZZ -->4 lepton decay
            pbar(mu,5) = v(mu,1) ! l+
            pbar(mu,6) = v(mu,2) ! l-
            pbar(mu,7) = v(mu,3) ! v
            pbar(mu,8) = v(mu,4) ! v~
         enddo
         fsign(5) = -1
         fsign(6) =  1
         fsign(7) =  1
         fsign(8) = -1
c
         if (lok(1) .or. nlo.eq.0) then
c routine "bqqzzqqj" with madgraph-call:
c            call bqqzzqqj_comp(pbar,fsign,qbar,gsign,bos,
c     &           uucc(1),uuss(1),ddcc(1),ddss(1),udsc(1),ducs(1))
c         
c new routine (with speed-up): uucc(1) stands for uucc(1:3) here!      
            call qqzzqqj(pbar,fsign,qbar,gsign,bos,
     &            uucc(1),uuss(1),ddcc(1),ddss(1),udsc(1),ducs(1)) 
	
c	compare to madgraph:
c            call bqqzzqqj_comp(pbar,fsign,qbar,gsign,bos,
c     &           uucc(1),uuss(1),ddcc(1),ddss(1),udsc(1),ducs(1))
     
        
         endif  !lok	

      elseif (bos.eq.211) then         ! ZA to 2 leptons + A
         do mu = 0,3             ! kinematics for qq-->ZA-->2l A decay
            pbar(mu,5) = v(mu,1) ! l+
            pbar(mu,6) = v(mu,2) ! l-
            pbar(mu,7) = v(mu,3) ! a
            pbar(mu,8) = 0d0
         enddo
         fsign(5) = -1
         fsign(6) =  1
         fsign(7) =  1
         fsign(7) =  1
!          udsc=temp
         if (lok(1) .or. nlo.eq.0) then
	    call qqzzqqj(pbar,fsign,qbar,gsign,bos,
     1                    uucc(1),uuss(1),ddcc(1),ddss(1),udsc(1),ducs(1))
         endif

       
!          if(.true.) then
! 	  pbar(0:3,1:4) = p(0:3,1:4,1)
! 	  pbar(0:3,8)=qbar(0:3)
! 	  call MG_ZA_UDSCG(Pbar,temp)
! 	  pbar(0:3,8) = 0d0
!          endif
! ! 	print*,fsign,nlo,bos
! !       print*,udsc
!       
!         print*,"udsc",temp,udsc(0),temp/udsc(0)

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
         if (lok(1) .or. nlo.eq.0) then
	    call qqzzqqj(pbar,fsign,qbar,gsign,bos,
     1                    uucc(1),uuss(1),ddcc(1),ddss(1),udsc(1),ducs(1))
	 endif

      else
         print*," Invalid entry BOS = ",bos,"in call of wbf_zh3j"
         print*," Must be 2, 6, 34 / 43, 22, 11 or 21 "
         stop
      endif   ! end of if bos = ..

      if (ldoblha) then
        return
      endif


      do j=2,3
         NCmatrixelt(0,0,j,1)=uucc(j)
         NCmatrixelt(0,1,j,1)=uuss(j)
         NCmatrixelt(1,0,j,1)=ddcc(j)
         NCmatrixelt(1,1,j,1)=ddss(j)
         CCmatrixelt(0,j,1)=udsc(j)
         CCmatrixelt(1,j,1)=ducs(j)
      enddo

c for the NLO case get the subtraction terms; 
c first the case with a final state gluon
      if (gsign.eq.1) then
         if (ldebug) then
            print*," final state gluon section in wbf_ZH3j "
            print*," jsig = ",jsig," jsig1,3 = ",jsig1,jsig3
            print 101," fsign = ",fsign
 101        format(a,6i5,a,i5)
         endif
         if (nlo.eq.1) then
            do j = 2,3               ! j=2: emission off 1,2 line
               do mu = 0,3           ! j=3: emission off 3,4 line
                  do i = 1,4
                     pbar(mu,physToDiag(i))=p(mu,i,j)
                  enddo
               enddo
               if (bos.eq.2 .and. lok(j)) then  ! Z
                  call qqZqq(pbar,fsign,0,
     1                       uucc(j),uuss(j),ddcc(j),
     2                       ddss(j),udsc(j),ducs(j))              
               elseif ((bos.eq.1) .and. lok(j)) then  ! photon
                  call qqAqq(pbar,fsign, 0,
     1                 uucc(j),uuss(j),ddcc(j),ddss(j),udsc(j),ducs(j))
 
               elseif (bos.eq.6 .and. lok(j)) then  ! Higgs
                  call qqHqq(pbar,fsign,0,
     1                       uucc(j),uuss(j),ddcc(j),
     2                       ddss(j),udsc(j),ducs(j),tree)
               elseif ((bos.eq.34 .or. bos.eq.43) .and. lok(j)) then  ! WW
c H->WW only:                
c		 call qqHqq(pbar,fsign,0,
c     1                       uucc(j),uuss(j),ddcc(j),
c     2                       ddss(j),udsc(j),ducs(j),tree)
c
c general qq->qqWW (madgraph):
c        	  call qqwwqq_comp(pbar,fsign, 0, 
c     1                    uucc(j),uuss(j),ddcc(j),ddss(j),
c     2			  udsc(j),ducs(j))
c
c  new routine for speed-up:
         	  call qqwwqq(pbar,fsign, 0,j,
     1                 uucc(j),uuss(j),ddcc(j),ddss(j),
     2                 udsc(j),ducs(j))

               elseif ((bos.eq.22) .and. lok(j)) then ! ZZ -> 4l
c general qq->qqZZ (madgraph):
c        	  call qqzzqq_comp(pbar,fsign, 0, bos,
c     1                    uucc(j),uuss(j),ddcc(j),ddss(j),
c     2			  udsc(j),ducs(j))

c     new routine for speed-up:
         	  call qqzzqq(pbar,fsign, 0,j,bos,
     1                    uucc(j),uuss(j),ddcc(j),ddss(j),
     2			  udsc(j),ducs(j))

               elseif (bos.eq.66) then ! HH production
                  call qqhhqq(pbar,fsign, 0,j,bos,
     1                    uucc(j),uuss(j),ddcc(j),ddss(j),
     2			  udsc(j),ducs(j))

               elseif ((bos.eq.11) .and. lok(j)) then

         	  call qqsp2aaqq(pbar,fsign, 0,j,bos,
     1                    uucc(j),uuss(j),ddcc(j),ddss(j),
     2			  udsc(j),ducs(j))


               elseif ((bos.eq.734) .and. lok(j)) then  ! spin-2 --> W+W- --> 4 leptons

         	  call qqsp2wwqq(pbar,fsign, 0,j,
     1                 uucc(j),uuss(j),ddcc(j),ddss(j),
     2                 udsc(j),ducs(j))

               elseif ((bos.eq.722) .and. lok(j)) then ! spin-2 --> ZZ -> 4l

         	  call qqsp2zzqq(pbar,fsign, 0,j,bos,
     1                    uucc(j),uuss(j),ddcc(j),ddss(j),
     2			  udsc(j),ducs(j))

               elseif ((bos.eq.721) .and. lok(j)) then  ! spin-2 --> ZZ -> 2l 2neutrinos

         	  call qqsp2zzqq(pbar,fsign, 0,j,bos,
     1                    uucc(j),uuss(j),ddcc(j),ddss(j),
     2			  udsc(j),ducs(j))

               elseif ((bos.eq.21) .and. lok(j)) then  ! ZZ -> 2l 2neutrinos
c general qq->qqZZ (madgraph):
c        	  call qqzzqq_comp(pbar,fsign, 0, bos,
c     1                    uucc(j),uuss(j),ddcc(j),ddss(j),
c     2			  udsc(j),ducs(j))

c     new routine for speed-up:
         	  call qqzzqq(pbar,fsign, 0,j,bos,
     1                    uucc(j),uuss(j),ddcc(j),ddss(j),
     2			  udsc(j),ducs(j))
     
               elseif ((bos.eq.211) .and. lok(j)) then ! ZA -> 2l A
c     new routine for speed-up:
         	  call qqzaqq(pbar,fsign, 0,j,bos,1,
     1                    uucc(j),uuss(j),ddcc(j),ddss(j),
     2			  udsc(j),ducs(j))
     
               elseif ((bos.eq.212) .and. lok(j)) then ! ZA -> 2nu A
c     new routine for speed-up:
         	  call qqzaqq(pbar,fsign, 0,j,bos,1,
     1                    uucc(j),uuss(j),ddcc(j),ddss(j),
     2			  udsc(j),ducs(j))
               elseif ((bos.eq.16.or.bos.eq.61) .and. lok(j)) then  ! Higgs+ptn
                  CALL HAjj_subtr_ME(fsign,j,uucc(j),uuss(j),ddcc(j),
     -                 ddss(j),udsc(j),ducs(j))
               else
                  uucc(j) = 0d0
                  uuss(j) = 0d0
                  ddcc(j) = 0d0
                  ddss(j) = 0d0
                  udsc(j) = 0d0
                  ducs(j) = 0d0 
               endif
	              
               NCmatrixelt(0,0,j,2)=uucc(j)
               NCmatrixelt(0,1,j,2)=uuss(j)
               NCmatrixelt(1,0,j,2)=ddcc(j)
               NCmatrixelt(1,1,j,2)=ddss(j)
               CCmatrixelt(0,j,2)=udsc(j)
               CCmatrixelt(1,j,2)=ducs(j)

c save matrix elements for later use with initial gluons
               uuccb(fsign(1),fsign(3),j) = uucc(j)
               uussb(fsign(1),fsign(3),j) = uuss(j)
               ddccb(fsign(1),fsign(3),j) = ddcc(j)
               ddssb(fsign(1),fsign(3),j) = ddss(j)
               udscb(fsign(1),fsign(3),j) = udsc(j)
               ducsb(fsign(1),fsign(3),j) = ducs(j)
c
               q2 = 2*dotrr(qbar,p(0,j-1,1))*xuz(1,j) !p(mu,j-1,1) is inc.quark
               v2 = 2d0*gs2(j)*                      !i.e  j-1 = 1,2
     $            ( 2d0/(1-xuz(1,j)+xuz(2,j))-(1+xuz(1,j)))
               q3 = 2*dotrr(qbar,p(0,j+1,1))*xuz(1,j) !p(mu,j+1,1) is out.quark
               v3 = 2d0*gs2(j)*                      !i.e. j+1 = 3,4
     $            ( 2d0/(1-xuz(1,j)+xuz(2,j))-(2-xuz(2,j)))
               sub(j) = v2/q2+v3/q3
            enddo  ! loop j
         endif   ! nlo=1
         iflav(5) = 0           ! final state gluon id


         do if1=1,nflVBF         !(nfl/2)*2
            do if2=1,nflVBF      !(nfl/2)*2
               iflav(1)=if1*fsign(physToDiag(1))
               iflav(3)=if1*fsign(physToDiag(3))
               iflav(2)=if2*fsign(physToDiag(2))
               iflav(4)=if2*fsign(physToDiag(4))
               do j = 2,3
                  k=FL_ZHg(iflav,j)
                  if (lok(1)) then
                     res(k,1)=pdf(sign1*iflav(1),1,2*j-3) !1 for j=2;
!3 for j=3
     &                    *pdf(sign2*iflav(2),2,4-j  ) !2 for j=2;1 for j=3
     &                    *NCmatrixelt(mod(if1,2),mod(if2,2),j,1)*polcolq
                  else
                     res(k,1) = 0
                  endif
                  if (nlo.eq.1 .and. lok(j)) then
                     res(k,j)=(pdf(sign1*iflav(1),1,j)
     &                    *pdf(sign2*iflav(2),2,j)
     &                    -pdf(sign1*iflav(1),1,2*j-3)
     &                    *pdf(sign2*iflav(2),2,4-j  )*sub(j))*C2
     &                    *NCmatrixelt(mod(if1,2),mod(if2,2),j,2)*polcolq
                     res(k,5-j)=0
                     res(k,0) = res(k,1)+res(k,j)
                  else
                     res(k,0) = res(k,1)
                     res(k,2) = 0
                     res(k,3) = 0
                  endif

c debugging for collinear and soft divergences -----------------------
                  if (lres .and. nlo.eq.1 .and. lok(j)) then
                    temp1 = sqrt(abs(dotrr(qbar,p(0,j-1,1)))) ! collinear 1
                    temp2 = sqrt(abs(dotrr(qbar,p(0,j+1,1)))) ! collinear 2
                    temp3 = abs(res(k,1)/res(k,j))               ! ratio RE / CT
                    if ( (temp3.ne.0.d0) .and. (qbar(0).le.0.2 .or. min(temp1,temp2).le.0.2) )
     &               write(31,*) "fs_nc", 
     &               qbar(0),
     &               temp1,
     &               temp2,
     &               temp3, j, res(k,1), res(k,j)
                  endif

               enddo
               if (ldebug .and. lok(1) .and. nlo.eq.1) then !debug
                  xa = res(k-1,1)+res(k,1) !debug
c           print*," same for res(k,*), k=",k," if1,2= ",if1,if2         !debug
                  print*,k," NC:up,low,tot",res(k-1,2)/xa,res(k,3)/xa, !debug
     &                 (res(k-1,2)+res(k,3))/xa !debug
               endif            !debug
C Now check if there is a CC contribution for this choice of initial state
C flavors (i,j). First: check whether initial uu or dd, not ud
               sametype=(mod(if1,2)).eq.(mod(if2,2))
               oneAntiparticle=
     &              ((fsign(physToDiag(1))*fsign(physToDiag(2))).eq.-1) 
c true if only one particle is an antiparticle
               ChargedCurrent=(oneAntiparticle.and.sametype) .or.
     &              ( (.not.oneAntiparticle).and.(.not.sametype) )
! never allow external b quarks for charged currents:
               if ( ChargedCurrent .and. if1.le.4 .and. if2.le.4 ) then   
c change 1<-->2 and 3<-->4 in outgoing quark flavors
                  iflav(3)=(if1+2*mod(if1,2)-1)*fsign(physToDiag(3))
                  iflav(4)=(if2+2*mod(if2,2)-1)*fsign(physToDiag(4))
                  icc1 = abs(iflav(diagtophys(1)))
                  do j = 2,3
                     k=FL_ZHg(iflav,j)
                     if (lok(1)) then
                        res(k,1)=pdf(sign1*iflav(1),1,2*j-3)
     &                       *pdf(sign2*iflav(2),2,4-j  )
     &                       *CCmatrixelt(mod(icc1,2),j,1)*polcolq
                     else
                        res(k,1) = 0
                     endif
                     if (nlo.eq.1 .and. lok(j)) then
                        res(k,j)=(pdf(sign1*iflav(1),1,j)
     &                       *pdf(sign2*iflav(2),2,j)
     &                       -pdf(sign1*iflav(1),1,2*j-3)
     &                       *pdf(sign2*iflav(2),2,4-j)*sub(j))*C2
     &                       *CCmatrixelt(mod(icc1,2),j,2)*polcolq
                        res(k,5-j) = 0
                        res(k,0) = res(k,1)+res(k,j)

                     else
                        res(k,0) = res(k,1)
                        res(k,2) = 0
                        res(k,3) = 0
                     endif

c debugging for collinear and soft divergences -----------------------
                  if (lres .and. nlo.eq.1 .and. lok(j)) then
                    temp1 = sqrt(abs(dotrr(qbar,p(0,j-1,1)))) ! collinear 1
                    temp2 = sqrt(abs(dotrr(qbar,p(0,j+1,1)))) ! collinear 2
                    temp3 = abs(res(k,1)/res(k,j))               ! ratio RE / CT
                    if ( (temp3.ne.0.d0) .and. (qbar(0).le.0.2 .or. min(temp1,temp2).le.0.2) )
     &               write(31,*) "fs_cc", 
     &               qbar(0),
     &               temp1,
     &               temp2,
     &               temp3, j, res(k,1), res(k,j)
                  endif

                  enddo         ! loop over upper vs. lower line
                  if (ldebug .and. (jsig.le.2.or.jsig.gt.8) !debug
     &                 .and. lok(1) .and. nlo.eq.1    ) then !debug
                     xa = res(k-1,1)+res(k,1) !debug
                     print*,k," CC:up,low,tot",res(k-1,2)/xa, !debug
     &                    res(k,3)/xa,(res(k-1,2)+res(k,3))/xa !debug
                  endif         !debug
               endif            ! endif(CC)
            enddo
         enddo 

      elseif (gsign.eq.-1) then           !initial gluon section
         j = 0
         if (nlo.eq.0) then
            if (fsign(1).eq.-fsign(2)) then 
               j=2              ! j=2: emission off 1,2 line
            elseif(fsign(3).eq.-fsign(4)) then 
               j=3              ! j=3: emission off 3,4 line
            endif
         elseif (nlo.eq.1) then
            do i=1,6
               fc(i)=fsign(i)
               fcb(i)=fsign(i)
            enddo
            if (fsign(1).eq.-fsign(2)) then 
               j=2              ! j=2: emission off 1,2 line
               fc(2)=-fsign(2)
               q2 = 2*dotrr(qbar,pbar(0,2))*xuz(1,j)
               fcb(1)=-fsign(1)
               q3 = 2*dotrr(qbar,pbar(0,1))*xuz(1,j)
            elseif(fsign(3).eq.-fsign(4)) then 
               j=3              ! j=3: emission off 3,4 line
               fc(4)=-fsign(4)
               q2 = 2*dotrr(qbar,pbar(0,4))*xuz(1,j)
               fcb(3)=-fsign(3)
               q3 = 2*dotrr(qbar,pbar(0,3))*xuz(1,j)
            endif
            v2 = 2d0*gs2(j) * ( (1-xuz(1,j))**2 + xuz(1,j)**2 )
            v3 = v2

            NCmatrixelt(0,0,2,2)=uuccb(fc(1),fc(3),j)
            NCmatrixelt(0,1,2,2)=uussb(fc(1),fc(3),j)
            NCmatrixelt(1,0,2,2)=ddccb(fc(1),fc(3),j)
            NCmatrixelt(1,1,2,2)=ddssb(fc(1),fc(3),j)
            CCmatrixelt(0,2,2)=udscb(fc(1),fc(3),j)
            CCmatrixelt(1,2,2)=ducsb(fc(1),fc(3),j)

            NCmatrixelt(0,0,3,2)=uuccb(fcb(1),fcb(3),j)
            NCmatrixelt(0,1,3,2)=uussb(fcb(1),fcb(3),j)
            NCmatrixelt(1,0,3,2)=ddccb(fcb(1),fcb(3),j)
            NCmatrixelt(1,1,3,2)=ddssb(fcb(1),fcb(3),j)
            CCmatrixelt(0,3,2)=udscb(fcb(1),fcb(3),j)
            CCmatrixelt(1,3,2)=ducsb(fcb(1),fcb(3),j)
         end if


         iflav(diagtophys(5)) = 0 ! Gluon is index 0 in pdf
         do if1=1,nflVBF             !(nfl/2)*2
            iflav(diagToPhys(1)) = if1 * fsign(1)
            iflav(diagToPhys(2)) = if1 * fsign(2)
            do if2=1,nflVBF          !(nfl/2)*2
               iflav(diagToPhys(3)) = if2 * fsign(3)
               iflav(diagToPhys(4)) = if2 * fsign(4)
               k=FL_ZHg(iflav,j)
               if (lok(1)) then
                  res(k,1)=pdf(sign1*iflav(1),1,2*j-3) !1 for j=2;3 for j=3
     &                 *pdf(sign2*iflav(2),2,4-j) !2 for j=2;1 for j=3
     &                 *NCmatrixelt(mod(if1,2),mod(if2,2),j,1)*polcolg
               else
                  res(k,1) = 0
               endif
               if ( nlo.eq.1 .and. lok(j) ) then
                  res(k,j) = ( pdf(sign1*iflav(1),1,j)
     &                 *pdf(sign2*iflav(2),2,j)
     &                 *(NCmatrixelt(mod(if1,2),mod(if2,2),2,2) +
     &                 NCmatrixelt(mod(if1,2),mod(if2,2),3,2))
     &                 - pdf(sign1*iflav(1),1,2*j-3)
     &                 * pdf(sign2*iflav(2),2,4-j) *
     &                 (NCmatrixelt(mod(if1,2),mod(if2,2),2,2)*v2/q2 +
     &                 NCmatrixelt(mod(if1,2),mod(if2,2),3,2)*v3/q3) )
     &                 * C2*polcolg
                  res(k,0) = res(k,1)+res(k,j)
                  res(k,5-j)=0
               else
                  res(k,0) = res(k,1)
                  res(k,2) = 0
                  res(k,3) = 0
               endif

c debugging for collinear and soft divergences -----------------------
                  if (lres .and. nlo.eq.1 .and. lok(j)) then
                    temp1 = sqrt(abs(dotrr(p(0,j-1,1),p(0,5,1)))) ! collinear 1
                    temp2 = sqrt(abs(dotrr(p(0,j-1,1),p(0,j+1,1)))) ! collinear 2
                    temp3 = abs(res(k,1)/res(k,j))               ! ratio RE / CT
                    if ( (temp3.ne.0.d0) .and. (qbar(0).le.0.2 .or. min(temp1,temp2).le.0.2) )
     &               write(31,*) "is_nc", 
     &               qbar(0),
     &               temp1,
     &               temp2,
     &               temp3, j, res(k,1), res(k,j)
                  endif

               if (ldebug .and. lok(1) .and. nlo.eq.1) then !debug
                  print*," k=",k," if1,2=",if1,if2 !debug
                  if (j.eq.2) then !debug
                     print*," upper line NC ",res(k,2)/res(k,1) !debug
                  elseif (j.eq.3) then !debug
                     print*," lower line NC ",res(k,3)/res(k,1) !debug
                  endif         !debug
               endif            !debug
            enddo
         enddo
C Next the CC contributions; no b quarks allowed here
         iflav(diagToPhys(5)) = 0
         do if1 = 1,(nfl/2)*2    ! 4 sum over all flavors for quark 1 and 2
            iflav(diagToPhys(1)) = fsign(1)*if1
            iflav(diagToPhys(2)) = fsign(2)*(if1-1+2*mod(if1,2))
            do if2 = mod(if1,2)+1,(nfl/2)*2,2 ! flavor of q3 set by q2 (mod 2)
               iflav(diagToPhys(3)) = fsign(3)*if2
               iflav(diagToPhys(4)) = fsign(4)*(if2-1+2*mod(if2,2))
               k=FL_ZHg(iflav,j)
               if (lok(1)) then
                  res(k,1)=pdf(sign1*iflav(1),1,2*j-3)
     &                 *pdf(sign2*iflav(2),2,4-j)
     &                 *CCmatrixelt(mod(if1,2),j,1)*polcolg
               else
                  res(k,1) = 0
               endif
               if ( nlo.eq.1 .and. lok(j) ) then
                  res(k,j)=(pdf(sign1*iflav(1),1,j)*
     &                 pdf(sign2*iflav(2),2,j)*
     &                 (CCmatrixelt(mod(if1,2),2,2) + 
     &                 CCmatrixelt(mod(if1,2),3,2)) 
     &                 - pdf(sign1*iflav(1),1,2*j-3)*
     &                 pdf(sign2*iflav(2),2,4-j)*
     &                 (CCmatrixelt(mod(if1,2),2,2)*v2/q2 +
     &                 CCmatrixelt(mod(if1,2),3,2)*v3/q3) )
     &                 * C2*polcolg
                  res(k,0) = res(k,1)+res(k,j)
                  res(k,5-j)=0
               else
                  res(k,0) = res(k,1)
                  res(k,2) = 0
                  res(k,3) = 0
               endif

c debugging for collinear and soft divergences -----------------------
                  if (lres .and. nlo.eq.1 .and. lok(j)) then
                    temp1 = sqrt(abs(dotrr(p(0,j-1,1),p(0,5,1)))) ! collinear 1
                    temp2 = sqrt(abs(dotrr(p(0,j-1,1),p(0,j+1,1)))) ! collinear 2
                    temp3 = abs(res(k,1)/res(k,j))               ! ratio RE / CT
                    if ( (temp3.ne.0.d0) .and. (qbar(0).le.0.2 .or. min(temp1,temp2).le.0.2) )
     &               write(31,*) "is_cc", 
     &               qbar(0),
     &               temp1,
     &               temp2,
     &               temp3, j, res(k,1), res(k,j)
                  endif

               if (ldebug .and. (jsig.le.2.or.jsig.gt.8) !debug
     &              .and. lok(1) .and. nlo.eq.1    ) then !debug
                  print*," CC case for k = ",k," if1,2=",if1,if2 !debug
                  print*," iflav = ",iflav !debug
                  if (j.eq.2)   !debug
     &                 print*," upper line CC ",res(k,2)/res(k,1) !debug
                  if (j.eq.3)   !debug
     &                 print*," lower line CC ",res(k,3)/res(k,1) !debug
               endif            !debug
            enddo
         enddo 

      endif

      nmax = FL_ZHg(iflav,-2)

      if (lres) then
        close(31)
      endif

      end
c******************************************************************
c
c   end subroutine wbf_zh3j
c
c*****************************************************************


c******************************************************************
c
c   begin function FL_ZHg
c
c*****************************************************************
      integer function FL_ZHg(iflav,colstruc)
      implicit none
      integer iflav(5),colstruc    ! input for color and flavor assignment
c  fill hepup color and flavor assignments for WBF processes with
c  one attached external gluon and count subprocesses. 
c  There are 2 color structures depending into which fermion line the 
c  gluon is inserted
c
c  colstruc = 2:upper
c             3:lower fermion line
c
c  cases can be distinguished according to which iflav(i)=0, i.e. corresponds
c  to the gluon

c  id1,2 are flavor identifiers for incoming quarks
c  id3,4 are flavor identifiers for outgoing quarks
c  gsign=+1: outgoing gluon, gsign=-1 for incoming gluon
c  colstruc =2,3 determines whether the gluon is coupled to the 
c  uppper or lower fermion line:

c  Note that colstruc is also a flag:
c  if colstruc=-1, then we are resetting.

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
          call fillColoredPartons_ZHg(iflav,
     1                                listposition+numdecay,colstruc)
        endif       ! numresets.eq.0
      endif         
      FL_ZHg=listposition

      end
c******************************************************************************
c
c   end function FL_ZHg
c
c******************************************************************************


c*****************************************************************************
c
c    begin  subroutine fillColoredPartons_ZHg
c
c*****************************************************************************
      subroutine fillColoredPartons_ZHg(iflav,listposition,colstruc)
c  assigns values to the variables in the common block localHEPUP
c  in particular, this subroutine assigns values to those variables that 
c  will be stored in the lookup tables generated by writeHEPUPtable.  
c  As the name suggests, this routine only stores the information for the 
c  colored partons.  Particles without color will be dealt with in the 
c  subroutine fillColorless.

      implicit none
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/process.inc"

      integer iflav(5),id1,id2,id3,id4,id5,listposition,i,i1
      integer iflavour(5),colstruc

      logical ldebug
      parameter (ldebug=.false.)


      do i1=1,5
         iflavour(i1)=iflav(i1) 
         if(iflavour(i1).eq.0) iflavour(i1)=21   
      enddo

c...flavours of 4 quarks and gluon

      id1=iflavour(1)
      id2=iflavour(2)
      id3=iflavour(3)
      id4=iflavour(4)
      id5=iflavour(5)
      select case(process)
      case (WMWMjjjLO, WMhadWMjjjLO)            ! WMWMjjj is obtained from W+W+jjj, flip IDs
        if (id1.ne.21) id1=-id1    
        if (id2.ne.21) id2=-id2
        if (id3.ne.21) id3=-id3
        if (id4.ne.21) id4=-id4
        if (id5.ne.21) id5=-id5
      end select

      select case (process)
      case(HjjjLO_WW, HjjjLO_ZZ_ll, HjjjLO_ZZ_lnu, WPWMjjjLO,
     &        ZZjjjLO_ll, ZZjjjLO_lnu,
     &        WPWpjjjLO, WMWMjjjLO,
     &        WPhadWMjjjLO, WPWMhadjjjLO, ZZhadjjjLO,
     &        HjjjLO_WPhadWM, HjjjLO_WPWMhad, HjjjLO_ZZhad,
     &        WPhadWPjjjLO, WMhadWMjjjLO, Sp2jjjLO_WW)
         lnup(listposition)=numParticles+3 
      case(ZAjjjLO,ZAjjjLO_n)
         lnup(listposition)=numParticles+1
      case(HjjjLO,AjjjLO)
         lnup(listposition)=numParticles-2
      case(HHjjjLO) ! modif julien
         lnup(listposition)=numParticles-1
      case(HHjjjLO_bbtau,HHjjjLO_bbAA) ! modif julien
         lnup(listposition)=numParticles+3
      case(HjjjLO_AA, HjjjLO_mu, HjjjLO_tau, HjjjLO_bbar, ZjjjLO_l,
     1     ZjjjLO_nu, WPjjjLO, WMjjjLO, AAjjjLO)
         lnup(listposition)=numParticles       
      case(HAjjjLO)
         lnup(listposition)=numParticles-1
      case(HAjjjLO_AA, HAjjjLO_mu, HAjjjLO_tau, HAjjjLO_bbar)
         lnup(listposition)=numParticles+1
      case(HAjjjLO_WW, HAjjjLO_ZZ_ll, HAjjjLO_ZZ_lnu)
         lnup(listposition)=numParticles+4          
      end select
      
      listup(1,listposition)=-1 !incoming quarks
      listup(2,listposition)=-1
      listup(3,listposition)=1  !outgoing quarks
      listup(4,listposition)=1
      listup(5,listposition)=1
 
c...four quarks and gluon
      lidup(1,listposition)=id1
      lidup(2,listposition)=id2
      lidup(3,listposition)=id3
      lidup(4,listposition)=id4
      lidup(5,listposition)=id5
            
      do i=3,n_p
         lmothup(1,i,listposition)=1
         lmothup(2,i,listposition)=2
         lspinup(i,listposition)=9
      enddo

c...final state gluon
      if(id5.eq.21)then
         if(colstruc.eq.2)then
            if(id1.gt.0) then
               licolup(1,1,listposition)=501
               licolup(1,3,listposition)=502
               licolup(1,5,listposition)=501
               licolup(2,1,listposition)=0
               licolup(2,3,listposition)=0
               licolup(2,5,listposition)=502
            else
               licolup(1,1,listposition)=0
               licolup(1,3,listposition)=0
               licolup(1,5,listposition)=502
               licolup(2,1,listposition)=501
               licolup(2,3,listposition)=502
               licolup(2,5,listposition)=501               
            endif            
            if(id2.gt.0) then
               licolup(1,2,listposition)=503
               licolup(1,4,listposition)=503
               licolup(2,2,listposition)=0
               licolup(2,4,listposition)=0
            else
               licolup(1,2,listposition)=0
               licolup(1,4,listposition)=0
               licolup(2,2,listposition)=503
               licolup(2,4,listposition)=503
            endif
         elseif(colstruc.eq.3)then
            if(id1.gt.0) then
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
            if(id2.gt.0) then
               licolup(1,2,listposition)=502
               licolup(1,4,listposition)=503
               licolup(1,5,listposition)=502
               licolup(2,2,listposition)=0
               licolup(2,4,listposition)=0
               licolup(2,5,listposition)=503
            else
               licolup(1,2,listposition)=0
               licolup(1,4,listposition)=0
               licolup(1,5,listposition)=503
               licolup(2,2,listposition)=502
               licolup(2,4,listposition)=503
               licolup(2,5,listposition)=502
            endif
         endif
      endif

c...initial particles gluon and quark/antiquark 
      if(id1.eq.21)then 
         if(id3.gt.0)then
            licolup(1,1,listposition)=502
            licolup(1,3,listposition)=502
            licolup(2,1,listposition)=501
            licolup(2,3,listposition)=0
         else
            licolup(1,1,listposition)=501
            licolup(1,3,listposition)=0
            licolup(2,1,listposition)=502
            licolup(2,3,listposition)=502
         endif 
         if(id5.gt.0)then
            licolup(1,5,listposition)=501
            licolup(2,5,listposition)=0
         else
            licolup(1,5,listposition)=0
            licolup(2,5,listposition)=501            
         endif
         if(id2.gt.0) then
            licolup(1,2,listposition)=503
            licolup(1,4,listposition)=503
            licolup(2,2,listposition)=0
            licolup(2,4,listposition)=0
         else
            licolup(1,2,listposition)=0
            licolup(1,4,listposition)=0
            licolup(2,2,listposition)=503
            licolup(2,4,listposition)=503
         endif         
      endif
      
c...initial particles quark/antiquark and gluon 
      if(id2.eq.21)then 
         if(id1.gt.0)then
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
         if(id5.gt.0)then
            licolup(1,5,listposition)=502
            licolup(2,5,listposition)=0
         else
            licolup(1,5,listposition)=0
            licolup(2,5,listposition)=502            
         endif
         if(id4.gt.0) then
            licolup(1,2,listposition)=503
            licolup(1,4,listposition)=503
            licolup(2,2,listposition)=502
            licolup(2,4,listposition)=0
         else
            licolup(1,2,listposition)=502
            licolup(1,4,listposition)=0
            licolup(2,2,listposition)=503
            licolup(2,4,listposition)=503
         endif
      endif

      if(ldebug)then
         print*,"====================================="
         print("(1P,5I6)"),id1,id2,id3,id4,id5
         print*,"-------------------------------------"
         print("(1P,5I6)"),(licolup(1,i1,listposition),i1=1,5)
         print("(1P,5I6)"),(licolup(2,i1,listposition),i1=1,5)
         print*,"====================================="
         print*,""
      endif
    
      end  ! fillColoredPartons_ZHg
c*****************************************************************************
c
c    end subroutine fillColoredPartons_ZHg
c
c*****************************************************************************

      subroutine resprint(nmin,nmax,res)
      implicit none
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
      integer nmin,nmax,k,i,ic
      save ic
      real*8 resold
      common /bresprint/ resold(16)
      real*8 xa,xb,xc,xd,res(maxnumsubproc,0:3),fac
      if (nmin.eq.1) ic=1
      xa = 0
      xb = 0
      xc = 0
      xd = 0
      if (nmin.le.190) then 
         fac = 4*9
      else
         fac = 4d0*24d0
      endif
      do k = nmin,nmax
         xa = xa + res(k,1)
         xb = xb + res(k,2)
         xc = xc + res(k,3)
         xd = xd + res(k,0)
      enddo
      if (nmax.lt.-50) then
         do k = nmin,nmax
            print*," res(",k,") ",(res(k,i),i=1,3)
         enddo
      endif
      write(6,10)" res([",nmin,",",nmax,"])",xa,xb/xa,xc/xa,xd/xa
c     1         ,xa/(resold(ic)+resold(ic+1))*fac
 10   format(a,i3,a,i3,a,g11.3,3f18.10)
      ic = ic+2
      end
