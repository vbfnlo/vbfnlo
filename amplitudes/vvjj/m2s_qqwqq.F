c
c	Dieter Zeppenfeld, <dieter@pheno.physics.wisc.edu>
c       Rob Szalapski, <robs@theory.kek.jp>
c	Initial version:  1995 March 24
c	Last modified by Barbara Jager:  May 2006
c
      real*8 function m2s_wbfw(
     &     bos,                 !in:  Boson identifier, 2 = Z, 6=H
     &     nlo,                 !in:  NLO=1: calculate virtual graphs
     &     xi,                  !in:  Feynman x parameters of incoming fermions
     &     p,                   !in:  momenta associated with external fermions
     &     v,                   !in:  momenta associated with WW decay fermions
     &     rn,                   !in:  random number array
     &     ps_number            !in:  Number of PS channel
     &                        )
          use globalvars, only: ldoblha
      implicit none
c
c declare input/output variables
c
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/lha.inc"
#include "tensorw.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/process.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

c variables for hadronic decays
      integer N_gen_W
      real*8 fac_W
      integer N_gen_up, N_gen_down
      real*8 fac_Z_up, fac_Z_down, mjj2
      external mjj2

      integer bos, nlo, ps_number
      real*8 xi(nx), p(0:3,max_p), v(0:3,max_v), rn(1), m2s

c  helicity selection
      integer h
      integer jsig, jsig1, jsig3, jsig5
      common /chelsum/ jsig,jsig1,jsig3,jsig5

      complex*16 zero
      parameter (zero=(0d0,0d0))
cc
c declare external functions
c
      integer FL_W
      external FL_W
      double complex dotrc
      external dotrc
cc
c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"
      real*8 q12(0:4,3), q34(0:4,3),qvv(0:3),qq
c
c declare local variables
c
c
      real*8 q_sf

      integer init/0/, I, nproc(0:4)
      save init, nproc
      integer j,mu
c
c declare variables for summation process
c
      INTEGER ires

      integer FSIGN(4+max_v), physToDiag(4), nmin, nmax
      save nmin, nmax
c
c store contributions from subprocess ip in res(ip)
      real*8 res(maxnumsubproc)
c
c declare parton distribution variables
c
      real*8 pdf(-6:6,2)
c
c define program switches
c
      logical ldebug
      parameter (ldebug=.false.)
c      data ldebug /.false./		!output debug information

      common /hcount / h

      real*8 weight,rnumb,RandomNumber


c -----------------------------------------------------------------------

      if ( init .gt. 1 ) then
         do i = 1,nmax
            res(i) = 0
         enddo
      else if (ldoblha) then
        h = 1
        init = 2
      else
      
	if (bos.eq.3.or.bos.eq.4) then    
	  
         write(6,*) " "
         write(6,*) "      Wjj amplitude square information:"
         write(6,*) "-----------------------------------------------"
         write(6,*) " "
	elseif (bos.eq.32) then      
         write(6,*) " "
         write(6,*) "      W+Zjj amplitude square information:"
         write(6,*) "-----------------------------------------------"
         write(6,*) " "
         call vtowz_reset
	elseif (bos.eq.42) then      
         write(6,*) " "
         write(6,*) "      W-Zjj amplitude square information:"
         write(6,*) "-----------------------------------------------"
         write(6,*) " "
         call vtowz_reset
      elseif (bos.eq.31) then      
         write(6,*) " "
         write(6,*) "      W+Ajj amplitude square information:"
         write(6,*) "-----------------------------------------------"
         write(6,*) " "
         call vtowz_reset
      elseif (bos.eq.41) then      
         write(6,*) " "
         write(6,*) "      W-Ajj amplitude square information:"
         write(6,*) "-----------------------------------------------"
         write(6,*) " "
         call vtowz_reset
	else
	 write(6,*) " this process is not implemented in m2sw"
	 stop
	endif !bos   

         call printnfl(.true.)
	 
         print*," creal = ",creal," cvirtual = ",cvirt
         print*," virt factor for alphas = 0.12 is ",
     1            1+0.12/pi*4./3.*cvirt
         h = 1
         init = 1
         do i = 1,maxnumsubproc
            res(i) = 0
         enddo
      endif


c -----------------------------------------------------------------------
      call Calc_Momentum_Transfer(p, v, q12,q34,1) 

c     setting helicities for W+jj, W-jj
      if (bos .eq. 3) then
         helicity(1) =  1       ! e+
         helicity(2) = -1       ! nu_e
      else if (bos .eq. 4) then
         helicity(1) =  1       ! nu_e_bar
         helicity(2) = -1       ! e-
      end if
      
      if (bos.eq.32) then       ! W+Z
c     reset lfs to .true. to force recalculation of virtual contributions in qqwwqq
         do i = 1,4
            lws(i) = .true.
         enddo
c     select helicity: h ... random number for lepton helicity
c     combination (h=1:2) 

         if (ldoblha) then
           h = blha_ranhelcomb
         endif
         h = mod(h,2) + 1
         
         ie = 1
         iu = sign(1,1-h)
         
c...Les Houches interface         
         if ((lha.or.hepmc).and..not.doNLO) then
            helicity(1)= ie
            helicity(2)=-ie 
            helicity(3)= iu
            helicity(4)=-iu
         endif
		
         CALL IXXXXX(v(0,1),ZERO ,+ie,-1,wl) !e+	
         CALL OXXXXX(v(0,2),ZERO ,-ie, 1,wn) !ve 
         CALL IXXXXX(v(0,3),ZERO ,+iu,-1,lp) !mu+	
         CALL OXXXXX(v(0,4),ZERO ,-iu, 1,lm) !mu- 
         
         CALL JIOXXX(wl,wn,GWF ,WMASS,WWIDTH,we) !Wl

         SELECT CASE(procid)
         CASE(WPZhadjj)
! for general output to all combinations: up-type first
         if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
! set couplings for Z hadronic decay into uubar
            call setZtouu
         endif
         CALL JIOXXX(lp,lm,GZ_ZF ,ZMASS,ZWIDTH,zu) !Zl
         CALL JIOXXX(lp,lm,GZ_AF ,ZERO ,ZERO  ,au) !Al
         case default
         CALL JIOXXX(lp,lm,GZL ,ZMASS,ZWIDTH,zu) !Zl
         CALL JIOXXX(lp,lm,GAL ,ZERO ,ZERO  ,au) !Al
         end select
         
         do mu = 0,3
            qe(mu) = v(mu,1)+v(mu,2)
            qu(mu) = v(mu,3)+v(mu,4)
            qvv(mu) = qe(mu) + qu(mu)
         enddo
	 	 
         qe(4) = qe(0)**2-qe(1)**2-qe(2)**2-qe(3)**2
         qu(4) = qu(0)**2-qu(1)**2-qu(2)**2-qu(3)**2
 
	        
	 if (abs(nlo).eq.5 .or. nlo.eq.-4) then ! determine V/V polarization vectors 
c wet(mu) = we(mu) - xxew*qe(mu) and zut(mu) = zu(mu) - xxmz*qu(mu)
            qq = qe(0)*qvv(0)-qe(1)*qvv(1)-qe(2)*qvv(2)-qe(3)*qvv(3)
            xxew = dotrc(qvv,we)/qq
            
	    qq = qu(0)*qvv(0)-qu(1)*qvv(1)-qu(2)*qvv(2)-qu(3)*qvv(3)
            xxuz = dotrc(qvv,zu)/qq            
            xxua = dotrc(qvv,au)/qq
            do mu = 0,3
               wet(mu+1) = we(mu+1) - xxew*qe(mu)
               zut(mu+1) = zu(mu+1) - xxuz*qu(mu)
               aut(mu+1) = au(mu+1) - xxua*qu(mu)
            enddo
            do mu = 5,6
               wet(mu) = we(mu)
               zut(mu) = zu(mu)
               aut(mu) = au(mu)
            enddo
         elseif (nlo.eq.1) then
            do mu = 1,6
               wet(mu) = we(mu)
               zut(mu) = zu(mu)
               aut(mu) = au(mu)
            enddo
         endif

c     leptonic tensors:

         SELECT CASE(procid)
      CASE(WPhadZjj, WPZhadjj)

         j = 1

c     using global form factor for all tensors of one phase space point
c     this ensures proper cancellations for anomalous contributions
c     energy scale is invariant WZ mass
         call anomal_formfactor(q12(0,j),q34(0,j),qe(0),qu(0))

c     for ZW+,AW+ -> e+ ve mu+ mu-:
         call vwtovwp_had(q12(0,j),q34(0,j),v,h,zwwv(0,0,j),awwv(0,0,j)) !q12 = V
         call vwtovwp_had(q34(0,j),q12(0,j),v,h,wzwv(0,0,j),wawv(0,0,j)) !q12 = W+
C     for V1V2 -->  mu+ mu- (Vi=A,Z)
         call VVtomm_had(1,h,q12(0,j),v,aauu(0,0,2,j),azuu(0,0,2,j),
     &        zauu(0,0,2,j),zzuu(0,0,2,j)) !emit W on lower
         call VVtomm_had(2,h,q34(0,j),v,aauu(0,0,1,j),azuu(0,0,1,j),
     &        zauu(0,0,1,j),zzuu(0,0,1,j)) !emit W on upper
C     for WW -->  mu+ mu-
         call WWtomm_had(1,h,q34(0,j),v,wwuu(0,0,1,j)) !emit V->e nu on upper
         call WWtomm_had(1,h,q12(0,j),v,wwuu(0,0,2,j)) !emit V->e nu on lower
C     for WV --> e+ ve
         call WPVtoln_had(1,q34(0,j),v,awen(0,0,1,j),zwen(0,0,1,j)) !emit V->mu+mu- on upper
         call WPVtoln_had(2,q12(0,j),v,awen(0,0,2,j),zwen(0,0,2,j)) !emit V->mu+mu- on lower
         call WPVtoln_had(1,q12(0,j),v,waen(0,0,2,j),wzen(0,0,2,j)) !emit V->mu+mu- on lower
         call WPVtoln_had(2,q34(0,j),v,waen(0,0,1,j),wzen(0,0,1,j)) !emit V->mu+mu- on upper  
c     for W+ --> e+ ve mu+ mu-
         call wpto4l_had(v,h,wwv)

      CASE DEFAULT

         if (with_kk) then      !Kaluza-Klein
#ifdef WITH_KK
            j = 1
c     for ZW+,AW+ -> e+ ve mu+ mu-:
            call vwtovwp_kk(q12(0,j),q34(0,j),v,h,zwwv(0,0,j),awwv(0,0,j)) !q12 = V
            call vwtovwp_kk(q34(0,j),q12(0,j),v,h,wzwv(0,0,j),wawv(0,0,j)) !q12 = W+
            
C     for V1V2 -->  mu+ mu- (Vi=A,Z)
            call VVtomm_kk(1,h,q12(0,j),v,aauu(0,0,2,j),azuu(0,0,2,j),
     $           zauu(0,0,2,j),zzuu(0,0,2,j)) !emit W on lower
            call VVtomm_kk(2,h,q34(0,j),v,aauu(0,0,1,j),azuu(0,0,1,j),
     $           zauu(0,0,1,j),zzuu(0,0,1,j)) !emit W on upper
            
C     for WW -->  mu+ mu-
            call WWtomm_kk(1,h,q34(0,j),v,wwuu(0,0,1,j)) !emit V->e nu on upper
            call WWtomm_kk(1,h,q12(0,j),v,wwuu(0,0,2,j)) !emit V->e nu on lower

C     for WV --> e+ ve
            call WPVtoln_kk(1,q34(0,j),v,awen(0,0,1,j),zwen(0,0,1,j)) !emit V->mu+mu- on upper
            call WPVtoln_kk(2,q12(0,j),v,awen(0,0,2,j),zwen(0,0,2,j)) !emit V->mu+mu- on lower
c        
            call WPVtoln_kk(1,q12(0,j),v,waen(0,0,2,j),wzen(0,0,2,j)) !emit V->mu+mu- on lower
            call WPVtoln_kk(2,q34(0,j),v,waen(0,0,1,j),wzen(0,0,1,j)) !emit V->mu+mu- on upper
            
c     for W+ --> e+ ve mu+ mu-
	
            call wpto4l_kk(v,h,wwv)	
#endif
            
         elseif (with_spin2) then    
            j = 1
c for ZW+,AW+ -> e+ ve mu+ mu-:
            call vwtovwp_spin2(q12(0,j),q34(0,j),v,h,zwwv(0,0,j),awwv(0,0,j)) !q12 = V
            call vwtovwp_spin2(q34(0,j),q12(0,j),v,h,wzwv(0,0,j),wawv(0,0,j)) !q12 = W+
            
c     other processes do not change for spin2 
C     for V1V2 -->  mu+ mu- (Vi=A,Z)
            call VVtomm(1,h,q12(0,j),v,aauu(0,0,2,j),azuu(0,0,2,j),
     $		 zauu(0,0,2,j),zzuu(0,0,2,j)) !emit W on lower
            call VVtomm(2,h,q34(0,j),v,aauu(0,0,1,j),azuu(0,0,1,j),
     $           zauu(0,0,1,j),zzuu(0,0,1,j)) !emit W on upper
            
C     for WW -->  mu+ mu-
            call WWtomm(1,h,q34(0,j),v,wwuu(0,0,1,j)) !emit V->e nu on upper
            call WWtomm(1,h,q12(0,j),v,wwuu(0,0,2,j)) !emit V->e nu on lower
         
C for WV --> e+ ve
            call WPVtoln(1,q34(0,j),v,awen(0,0,1,j),zwen(0,0,1,j)) !emit V->mu+mu- on upper
            call WPVtoln(2,q12(0,j),v,awen(0,0,2,j),zwen(0,0,2,j)) !emit V->mu+mu- on lower
c        
            call WPVtoln(1,q12(0,j),v,waen(0,0,2,j),wzen(0,0,2,j)) !emit V->mu+mu- on lower
            call WPVtoln(2,q34(0,j),v,waen(0,0,1,j),wzen(0,0,1,j)) !emit V->mu+mu- on upper  
	 
c     for W+ --> e+ ve mu+ mu-
	
            call wpto4l(v,h,wwv)

         elseif (with_anom) then ! anomalous gauge boson couplings  

            j = 1

c     using global form factor for all tensors of one phase space point
c     this ensures proper cancellations for anomalous contributions
c     energy scale is invariant WZ mass
            call anomal_formfactor(q12(0,j),q34(0,j),qe(0),qu(0))
        
c for ZW+,AW+ -> e+ ve mu+ mu-:
            call vwtovwp_anomal(q12(0,j),q34(0,j),v,h,zwwv(0,0,j),awwv(0,0,j)) !q12 = V
            call vwtovwp_anomal(q34(0,j),q12(0,j),v,h,wzwv(0,0,j),wawv(0,0,j)) !q12 = W+

C for V1V2 -->  mu+ mu- (Vi=A,Z)
            call VVtomm_anomal(1,h,q12(0,j),v,aauu(0,0,2,j),azuu(0,0,2,j),
     $           zauu(0,0,2,j),zzuu(0,0,2,j)) !emit W on lower
            call VVtomm_anomal(2,h,q34(0,j),v,aauu(0,0,1,j),azuu(0,0,1,j),
     $           zauu(0,0,1,j),zzuu(0,0,1,j)) !emit W on upper

C for WW -->  mu+ mu-
            call WWtomm_anomal(1,h,q34(0,j),v,wwuu(0,0,1,j)) !emit V->e nu on upper
            call WWtomm_anomal(1,h,q12(0,j),v,wwuu(0,0,2,j)) !emit V->e nu on lower

C for WV --> e+ ve
            call WPVtoln_anomal(1,q34(0,j),v,awen(0,0,1,j),zwen(0,0,1,j)) !emit V->mu+mu- on upper
            call WPVtoln_anomal(2,q12(0,j),v,awen(0,0,2,j),zwen(0,0,2,j)) !emit V->mu+mu- on lower
c        
            call WPVtoln_anomal(1,q12(0,j),v,waen(0,0,2,j),wzen(0,0,2,j)) !emit V->mu+mu- on lower
            call WPVtoln_anomal(2,q34(0,j),v,waen(0,0,1,j),wzen(0,0,1,j)) !emit V->mu+mu- on upper  
	 
c     for W+ --> e+ ve mu+ mu-
	
            call wpto4l_anomal(v,h,wwv)


         else                   ! SM

            j = 1
c     for ZW+,AW+ -> e+ ve mu+ mu-:
            call vwtovwp(q12(0,j),q34(0,j),v,h,zwwv(0,0,j),awwv(0,0,j)) !q12 = V
            call vwtovwp(q34(0,j),q12(0,j),v,h,wzwv(0,0,j),wawv(0,0,j)) !q12 = W+
           
C for V1V2 -->  mu+ mu- (Vi=A,Z)
            call VVtomm(1,h,q12(0,j),v,aauu(0,0,2,j),azuu(0,0,2,j),
     $           zauu(0,0,2,j),zzuu(0,0,2,j)) !emit W on lower
            call VVtomm(2,h,q34(0,j),v,aauu(0,0,1,j),azuu(0,0,1,j),
     $		 zauu(0,0,1,j),zzuu(0,0,1,j)) !emit W on upper

C     for WW -->  mu+ mu-
            call WWtomm(1,h,q34(0,j),v,wwuu(0,0,1,j)) !emit V->e nu on upper
            call WWtomm(1,h,q12(0,j),v,wwuu(0,0,2,j)) !emit V->e nu on lower

C     for WV --> e+ ve
            call WPVtoln(1,q34(0,j),v,awen(0,0,1,j),zwen(0,0,1,j)) !emit V->mu+mu- on upper
            call WPVtoln(2,q12(0,j),v,awen(0,0,2,j),zwen(0,0,2,j)) !emit V->mu+mu- on lower
c     
            call WPVtoln(1,q12(0,j),v,waen(0,0,2,j),wzen(0,0,2,j)) !emit V->mu+mu- on lower
            call WPVtoln(2,q34(0,j),v,waen(0,0,1,j),wzen(0,0,1,j)) !emit V->mu+mu- on upper
         
	 
c     for W+ --> e+ ve mu+ mu-
         
            call wpto4l(v,h,wwv)	
         
      
         endif
         
      end select                ! semileptonic decay processes
      
c     -----------------------------------------------------------------------
      
      elseif ( bos.eq.42) then  !W-Z
c     reset lws to .true. to force recalculation of virtual contributions in qqwwqq
         do i = 1,4
            lws(i) = .true.
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
            helicity(1)=-ie
            helicity(2)= ie 
            helicity(3)= iu
            helicity(4)=-iu
         endif

	 CALL OXXXXX(v(0,1),ZERO ,-ie, 1,wl) !e-       
	 CALL IXXXXX(v(0,2),ZERO ,+ie,-1,wn) !ve~ 
	 CALL IXXXXX(v(0,3),ZERO ,+iu,-1,lp) !mu+      
	 CALL OXXXXX(v(0,4),ZERO ,-iu, 1,lm) !mu- 
         
	 CALL JIOXXX(wn,wl,GWF ,WMASS,WWIDTH,we) !Wl

         SELECT CASE(procid)
         CASE(WMZhadjj)
!     for general output to all combinations: up-type first
            if ((finalquarks(1).eq.93 .or. finalquarks(1).eq.94)) then
!     set couplings for Z hadronic decay into uubar
               call setZtouu
            endif
            CALL JIOXXX(lp,lm,GZ_ZF ,ZMASS,ZWIDTH,zu) !Zl
            CALL JIOXXX(lp,lm,GZ_AF ,ZERO ,ZERO  ,au) !Al	
         case default
            CALL JIOXXX(lp,lm,GZL ,ZMASS,ZWIDTH,zu) !Zl
            CALL JIOXXX(lp,lm,GAL ,ZERO ,ZERO  ,au) !Al	
         end select
c     
         do mu = 0,3
            qe(mu) = v(mu,1)+v(mu,2)
            qu(mu) = v(mu,3)+v(mu,4)
            qvv(mu) = qe(mu) + qu(mu)
         enddo
	 	 
         qe(4) = qe(0)**2-qe(1)**2-qe(2)**2-qe(3)**2
         qu(4) = qu(0)**2-qu(1)**2-qu(2)**2-qu(3)**2
 
	        
         if (abs(nlo).eq.5 .or. nlo.eq.-4) then ! determine V/V polarization vectors 
c     wet(mu) = we(mu) - xxew*qe(mu) and zut(mu) = zu(mu) - xxmz*qu(mu)
            qq = qe(0)*qvv(0)-qe(1)*qvv(1)-qe(2)*qvv(2)-qe(3)*qvv(3)
            xxew = dotrc(qvv,we)/qq
            
	    qq = qu(0)*qvv(0)-qu(1)*qvv(1)-qu(2)*qvv(2)-qu(3)*qvv(3)
            xxuz = dotrc(qvv,zu)/qq            
            xxua = dotrc(qvv,au)/qq
            do mu = 0,3
               wet(mu+1) = we(mu+1) - xxew*qe(mu)
               zut(mu+1) = zu(mu+1) - xxuz*qu(mu)
               aut(mu+1) = au(mu+1) - xxua*qu(mu)
            enddo
            do mu = 5,6
               wet(mu) = we(mu)
               zut(mu) = zu(mu)
               aut(mu) = au(mu)
            enddo
         elseif (nlo.eq.1) then
            do mu = 1,6
               wet(mu) = we(mu)
               zut(mu) = zu(mu)
               aut(mu) = au(mu)
            enddo
         endif
         
c     leptonic tensors:


         SELECT CASE(procid)
         CASE(WMhadZjj, WMZhadjj)
         
            j = 1

c     using global form factor for all tensors of one phase space point
c     this ensures proper cancellations for anomalous contributions
c     energy scale is invariant WZ mass
            call anomal_formfactor(q12(0,j),q34(0,j),qe(0),qu(0))
            
c     for ZW-,AW- -> e- ve~ mu+ mu-:
            call vwtovwm_had(q12(0,j),q34(0,j),v,h,zwwv(0,0,j),awwv(0,0,j)) !q12 = V
            call vwtovwm_had(q34(0,j),q12(0,j),v,h,wzwv(0,0,j),wawv(0,0,j)) !q12 = W-
C     for V1V2 -->  mu+ mu- (Vi=A,Z)
            call VVtomm_had(1,h,q12(0,j),v,aauu(0,0,2,j),azuu(0,0,2,j),
     $		 zauu(0,0,2,j),zzuu(0,0,2,j)) !emit W on lower
            call VVtomm_had(2,h,q34(0,j),v,aauu(0,0,1,j),azuu(0,0,1,j),
     $           zauu(0,0,1,j),zzuu(0,0,1,j)) !emit W on upper
C     for W+W- -->  mu+ mu-
            call WWtomm_had(2,h,q34(0,j),v,wwuu(0,0,1,j)) !emit V->e nu on upper
            call WWtomm_had(2,h,q12(0,j),v,wwuu(0,0,2,j)) !emit V->e nu on lower
C     for WV --> e- ve~: 
            call WMVtoln_had(1,q34(0,j),v,awen(0,0,1,j),zwen(0,0,1,j)) !emit V->mu+mu- on upper
            call WMVtoln_had(2,q12(0,j),v,awen(0,0,2,j),zwen(0,0,2,j)) !emit V->mu+mu- on lower
            call WMVtoln_had(1,q12(0,j),v,waen(0,0,2,j),wzen(0,0,2,j)) !emit V->mu+mu- on lower
            call WMVtoln_had(2,q34(0,j),v,waen(0,0,1,j),wzen(0,0,1,j)) !emit V->mu+mu- on upper
c     for W- --> e- ve mu+ mu-
            call wmto4l_had(v,h,wwv)


         CASE DEFAULT

            if (with_kk) then   !Kaluza-Klein
#ifdef WITH_KK
               j = 1
	  
c     for ZW-,AW- -> e- ve~ mu+ mu-:
               call vwtovwm_kk(q12(0,j),q34(0,j),v,h,zwwv(0,0,j),awwv(0,0,j)) !q12 = V
               call vwtovwm_kk(q34(0,j),q12(0,j),v,h,wzwv(0,0,j),wawv(0,0,j)) !q12 = W-
           
           
C     for V1V2 -->  mu+ mu- (Vi=A,Z)
               call VVtomm_kk(1,h,q12(0,j),v,aauu(0,0,2,j),azuu(0,0,2,j),
     $              zauu(0,0,2,j),zzuu(0,0,2,j)) !emit W on lower
               call VVtomm_kk(2,h,q34(0,j),v,aauu(0,0,1,j),azuu(0,0,1,j),
     $              zauu(0,0,1,j),zzuu(0,0,1,j)) !emit W on upper

C     for W+W- -->  mu+ mu-
               call WWtomm_kk(2,h,q34(0,j),v,wwuu(0,0,1,j)) !emit V->e nu on upper
               call WWtomm_kk(2,h,q12(0,j),v,wwuu(0,0,2,j)) !emit V->e nu on lower

C     for WV --> e- ve~: 
               call WMVtoln_kk(1,q34(0,j),v,awen(0,0,1,j),zwen(0,0,1,j)) !emit V->mu+mu- on upper
               call WMVtoln_kk(2,q12(0,j),v,awen(0,0,2,j),zwen(0,0,2,j)) !emit V->mu+mu- on lower
c        
               call WMVtoln_kk(1,q12(0,j),v,waen(0,0,2,j),wzen(0,0,2,j)) !emit V->mu+mu- on lower
               call WMVtoln_kk(2,q34(0,j),v,waen(0,0,1,j),wzen(0,0,1,j)) !emit V->mu+mu- on upper
		 
c     for W- --> e- ve mu+ mu-
	
               call wmto4l_kk(v,h,wwv)	
#endif

            elseif (with_spin2) then    
               j = 1
c     for ZW-,AW- -> e- ve~ mu+ mu-:
               call vwtovwm_spin2(q12(0,j),q34(0,j),v,h,zwwv(0,0,j),awwv(0,0,j)) !q12 = V
               call vwtovwm_spin2(q34(0,j),q12(0,j),v,h,wzwv(0,0,j),wawv(0,0,j)) !q12 = W-

c     other processes do not change for spin2
          
C     for V1V2 -->  mu+ mu- (Vi=A,Z)
        call VVtomm(1,h,q12(0,j),v,aauu(0,0,2,j),azuu(0,0,2,j),
     $              zauu(0,0,2,j),zzuu(0,0,2,j)) !emit W on lower
        call VVtomm(2,h,q34(0,j),v,aauu(0,0,1,j),azuu(0,0,1,j),
     $       zauu(0,0,1,j),zzuu(0,0,1,j)) !emit W on upper
        
C     for W+W- -->  mu+ mu-
        call WWtomm(2,h,q34(0,j),v,wwuu(0,0,1,j)) !emit V->e nu on upper
        call WWtomm(2,h,q12(0,j),v,wwuu(0,0,2,j)) !emit V->e nu on lower
         
C     for WV --> e- ve~: 
        call WMVtoln(1,q34(0,j),v,awen(0,0,1,j),zwen(0,0,1,j)) !emit V->mu+mu- on upper
        call WMVtoln(2,q12(0,j),v,awen(0,0,2,j),zwen(0,0,2,j)) !emit V->mu+mu- on lower
c     
        call WMVtoln(1,q12(0,j),v,waen(0,0,2,j),wzen(0,0,2,j)) !emit V->mu+mu- on lower
        call WMVtoln(2,q34(0,j),v,waen(0,0,1,j),wzen(0,0,1,j)) !emit V->mu+mu- on upper
		 
c     for W- --> e- ve mu+ mu-
	
	call wmto4l(v,h,wwv)	

      elseif (with_anom) then   ! anomalous gauge boson couplings  

               j = 1

c     using global form factor for all tensors of one phase space point
c     this ensures proper cancellations for anomalous contributions
c     energy scale is invariant WZ mass
               call anomal_formfactor(q12(0,j),q34(0,j),qe(0),qu(0))

c     for ZW-,AW- -> e- ve~ mu+ mu-:
               call vwtovwm_anomal(q12(0,j),q34(0,j),v,h,zwwv(0,0,j),awwv(0,0,j)) !q12 = V
        call vwtovwm_anomal(q34(0,j),q12(0,j),v,h,wzwv(0,0,j),wawv(0,0,j)) !q12 = W-
        
C     for V1V2 -->  mu+ mu- (Vi=A,Z)
        call VVtomm_anomal(1,h,q12(0,j),v,aauu(0,0,2,j),azuu(0,0,2,j),
     $       zauu(0,0,2,j),zzuu(0,0,2,j)) !emit W on lower
        call VVtomm_anomal(2,h,q34(0,j),v,aauu(0,0,1,j),azuu(0,0,1,j),
     $       zauu(0,0,1,j),zzuu(0,0,1,j)) !emit W on upper

C     for W+W- -->  mu+ mu-
        call WWtomm_anomal(2,h,q34(0,j),v,wwuu(0,0,1,j)) !emit V->e nu on upper
        call WWtomm_anomal(2,h,q12(0,j),v,wwuu(0,0,2,j)) !emit V->e nu on lower
        
C for WV --> e- ve~: 
        call WMVtoln_anomal(1,q34(0,j),v,awen(0,0,1,j),zwen(0,0,1,j)) !emit V->mu+mu- on upper
        call WMVtoln_anomal(2,q12(0,j),v,awen(0,0,2,j),zwen(0,0,2,j)) !emit V->mu+mu- on lower
c     
        call WMVtoln_anomal(1,q12(0,j),v,waen(0,0,2,j),wzen(0,0,2,j)) !emit V->mu+mu- on lower
        call WMVtoln_anomal(2,q34(0,j),v,waen(0,0,1,j),wzen(0,0,1,j)) !emit V->mu+mu- on upper
        
c     for W- --> e- ve mu+ mu-
	
	call wmto4l_anomal(v,h,wwv)	
        
        
      else                      !SM
         j = 1
         
c for ZW-,AW- -> e- ve~ mu+ mu-:
         call vwtovwm(q12(0,j),q34(0,j),v,h,zwwv(0,0,j),awwv(0,0,j)) !q12 = V
         call vwtovwm(q34(0,j),q12(0,j),v,h,wzwv(0,0,j),wawv(0,0,j)) !q12 = W-
         
C     for V1V2 -->  mu+ mu- (Vi=A,Z)
         call VVtomm(1,h,q12(0,j),v,aauu(0,0,2,j),azuu(0,0,2,j),
     $        zauu(0,0,2,j),zzuu(0,0,2,j)) !emit W on lower
        call VVtomm(2,h,q34(0,j),v,aauu(0,0,1,j),azuu(0,0,1,j),
     $        zauu(0,0,1,j),zzuu(0,0,1,j)) !emit W on upper

C     for W+W- -->  mu+ mu-
        call WWtomm(2,h,q34(0,j),v,wwuu(0,0,1,j)) !emit V->e nu on upper
        call WWtomm(2,h,q12(0,j),v,wwuu(0,0,2,j)) !emit V->e nu on lower

C     for WV --> e- ve~: 
        call WMVtoln(1,q34(0,j),v,awen(0,0,1,j),zwen(0,0,1,j)) !emit V->mu+mu- on upper
        call WMVtoln(2,q12(0,j),v,awen(0,0,2,j),zwen(0,0,2,j)) !emit V->mu+mu- on lower
c     
        call WMVtoln(1,q12(0,j),v,waen(0,0,2,j),wzen(0,0,2,j)) !emit V->mu+mu- on lower
        call WMVtoln(2,q34(0,j),v,waen(0,0,1,j),wzen(0,0,1,j)) !emit V->mu+mu- on upper
		 
c     for W- --> e- ve mu+ mu-
	
	call wmto4l(v,h,wwv)	

      endif

      end select                ! semileptonic decay processes
    

      elseif (bos.eq.31) then   !W+A
c     reset lws to .true. to force recalculation of virtual contributions in qqwpaqq
         do i = 1,4
            lws(i) = .true.
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
            helicity(1)= ie
            helicity(2)=-ie 
            helicity(3)= iu
         endif

         CALL IXXXXX(v(0,1),ZERO ,+ie,-1,wl) !e+	
         CALL OXXXXX(v(0,2),ZERO ,-ie, 1,wn) !ve 
         CALL VXXXXX(v(0,3),ZERO ,+iu, 1,ea) !a

         CALL JIOXXX(wl,wn,GWF,WMASS,WWIDTH,we) !Wl

         do mu = 0,3
            qe(mu) = v(mu,1)+v(mu,2)
            qa(mu) = v(mu,3)
            qvv(mu) = qe(mu) + qa(mu)
         enddo
         
         qe(4) = qe(0)**2-qe(1)**2-qe(2)**2-qe(3)**2
         qa(4) = qa(0)**2-qa(1)**2-qa(2)**2-qa(3)**2

c     leptonic tensors:

         if (with_kk) then      !Kaluza-Klein
#ifdef WITH_KK
            j = 1
#endif
         elseif (with_spin2) then    
            j = 1
         elseif (with_anom) then ! anomalous gauge boson couplings  

            j = 1

c     using global form factor for all tensors of one phase space point
c     this ensures proper cancellations for anomalous contributions
c     energy scale is invariant WA mass
            call anomal_formfactor(q12(0,j),q34(0,j),qe(0),qa(0))

c     for VW+ -> e+ ve A:
            call vwptowpa_anomal(q12(0,j),q34(0,j),v,h,zwtowa(0,0,j),awtowa(0,0,j)) !q12 = V
            call vwptowpa_anomal(q34(0,j),q12(0,j),v,h,wztowa(0,0,j),watowa(0,0,j)) !q12 = W+        


C     for W+W- -->  A
            call calc_wwtoa_anomal(1,iu,q34(0,j),v,wwtoa(0,0,1,j)) !emit W+ -> e+ ve on upper
            call calc_wwtoa_anomal(1,iu,q12(0,j),v,wwtoa(0,0,2,j)) !emit W+ -> e+ ve on lower

C     for WV --> e+ ve
            call WPVtoln_anomal(1,q34(0,j),v,awen(0,0,1,j),zwen(0,0,1,j)) !emit V->mu+mu- on upper
            call WPVtoln_anomal(2,q12(0,j),v,awen(0,0,2,j),zwen(0,0,2,j)) !emit V->mu+mu- on lower
c     
            call WPVtoln_anomal(1,q12(0,j),v,waen(0,0,2,j),wzen(0,0,2,j)) !emit V->mu+mu- on lower
            call WPVtoln_anomal(2,q34(0,j),v,waen(0,0,1,j),wzen(0,0,1,j)) !emit V->mu+mu- on upper  

c     for W+ --> e+ ve~ A
            call wptowpa_anomal(qvv,v,h,wtowa)	 

         else                   !SM
            j = 1

c     for VW+ -> e+ ve A:
            call vwptowpa(q12(0,j),q34(0,j),v,h,zwtowa(0,0,j),awtowa(0,0,j)) !q12 = V
            call vwptowpa(q34(0,j),q12(0,j),v,h,wztowa(0,0,j),watowa(0,0,j)) !q12 = W+

C     for W+W- -->  A
            call calc_wwtoa(1,iu,q34(0,j),v,wwtoa(0,0,1,j)) !emit W+ -> e+ ve on upper
            call calc_wwtoa(1,iu,q12(0,j),v,wwtoa(0,0,2,j)) !emit W+ -> e+ ve on lower

C     for WV --> e+ ve
            call WPVtoln(1,q34(0,j),v,awen(0,0,1,j),zwen(0,0,1,j)) !emit A on upper
            call WPVtoln(2,q12(0,j),v,awen(0,0,2,j),zwen(0,0,2,j)) !emit A on lower
c     
            call WPVtoln(1,q12(0,j),v,waen(0,0,2,j),wzen(0,0,2,j)) !emit A on lower
            call WPVtoln(2,q34(0,j),v,waen(0,0,1,j),wzen(0,0,1,j)) !emit A on upper

c     for W+ --> e+ ve~ A
            call wptowpa(qvv,v,h,wtowa)

         endif                  ! koppln  
         
      elseif (bos.eq.41) then   !W-A
c     reset lws to .true. to force recalculation of virtual contributions in qqwwqq
         do i = 1,4
            lws(i) = .true.
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
            helicity(1)=-ie
            helicity(2)= ie 
            helicity(3)= iu
         endif

	 CALL OXXXXX(v(0,1),ZERO ,-ie, 1,wl) !e-       
	 CALL IXXXXX(v(0,2),ZERO ,+ie,-1,wn) !ve~ 
	 CALL VXXXXX(v(0,3),ZERO ,+iu, 1,ea) !a

	 CALL JIOXXX(wn,wl,GWF ,WMASS,WWIDTH,we) !Wl
         
         do mu = 0,3
            qe(mu) = v(mu,1)+v(mu,2)
            qa(mu) = v(mu,3)
            qvv(mu) = qe(mu) + qa(mu)
         enddo
         
         qe(4) = qe(0)**2-qe(1)**2-qe(2)**2-qe(3)**2
         qa(4) = qa(0)**2-qa(1)**2-qa(2)**2-qa(3)**2

c     leptonic tensors:

         if (with_kk) then      !Kaluza-Klein
#ifdef WITH_KK
            j = 1
            
#endif

         elseif (with_spin2) then    
            j = 1 	
         elseif (with_anom) then ! anomalous gauge boson couplings  

            j = 1

c     using global form factor for all tensors of one phase space point
c     this ensures proper cancellations for anomalous contributions
c     energy scale is invariant WA mass
            call anomal_formfactor(q12(0,j),q34(0,j),qe(0),qa(0))

c     for ZW-,AW- -> e- ve~ A:
            call vwmtowma_anomal(q12(0,j),q34(0,j),v,h,zwtowa(0,0,j),awtowa(0,0,j)) !q12 = V
            call vwmtowma_anomal(q34(0,j),q12(0,j),v,h,wztowa(0,0,j),watowa(0,0,j)) !q12 = W+        


C     for W+W- -->  A
            call calc_wwtoa_anomal(2,iu,q34(0,j),v,wwtoa(0,0,1,j)) !emit W+ -> e+ ve on upper
            call calc_wwtoa_anomal(2,iu,q12(0,j),v,wwtoa(0,0,2,j)) !emit W+ -> e+ ve on lower

C     for WV --> e- ve~: 
            call WMVtoln_anomal(1,q34(0,j),v,awen(0,0,1,j),zwen(0,0,1,j)) !emit V->mu+mu- on upper
            call WMVtoln_anomal(2,q12(0,j),v,awen(0,0,2,j),zwen(0,0,2,j)) !emit V->mu+mu- on lower
c     
            call WMVtoln_anomal(1,q12(0,j),v,waen(0,0,2,j),wzen(0,0,2,j)) !emit V->mu+mu- on lower
            call WMVtoln_anomal(2,q34(0,j),v,waen(0,0,1,j),wzen(0,0,1,j)) !emit V->mu+mu- on upper  

c     for W- --> e- ve~ A
            call wmtowma_anomal(qvv,v,h,wtowa)	 

         else                   !SM
            j = 1
            
c     for ZW-,AW- -> e- ve~ A:
            call vwmtowma(q12(0,j),q34(0,j),v,h,zwtowa(0,0,j),awtowa(0,0,j)) !q12 = V
            call vwmtowma(q34(0,j),q12(0,j),v,h,wztowa(0,0,j),watowa(0,0,j)) !q12 = W-
            
C     for W+W- -->  A
            call calc_wwtoa(2,iu,q34(0,j),v,wwtoa(0,0,1,j)) !emit W- -> e- nu~ on upper
            call calc_wwtoa(2,iu,q12(0,j),v,wwtoa(0,0,2,j)) !emit W- -> e- nu~ on lower

C     for WV --> e- ve~: 
            call WMVtoln(1,q34(0,j),v,awen(0,0,1,j),zwen(0,0,1,j)) !emit V->mu+mu- on upper
            call WMVtoln(2,q12(0,j),v,awen(0,0,2,j),zwen(0,0,2,j)) !emit V->mu+mu- on lower
c     
            call WMVtoln(1,q12(0,j),v,waen(0,0,2,j),wzen(0,0,2,j)) !emit V->mu+mu- on lower
            call WMVtoln(2,q34(0,j),v,waen(0,0,1,j),wzen(0,0,1,j)) !emit V->mu+mu- on upper
            
c     for W- --> e- ve~ A
            call wmtowma(qvv,v,h,wtowa)

         endif
c -----------------------------------------------------------------------
     
#ifdef WITH_NLO      
      elseif ( nlo.eq.1.and.(bos.eq.3.or.bos.eq.4) ) then
         call BCD_fill(p(0,1),p(0,2),                   !beam
     1                 p(0,3),p(0,4),                   !jet
     2                 v(0,1),v(0,2))                       !decay
#endif
      endif
      
c
c if desired, print output information
c
      if ( ldebug ) then
 10      format( " p(", i1, ") = ", 4(f10.3, 2x) )
 20      format( " v(", i1, ") = ", 4(f10.3, 2x) )
c
         write(6,*) " "
         do i = 1, n_v
            write(6,20) i, v(0,i), v(1,i), v(2,i), v(3,i)
         end do
         print*," parton momenta "
         do i = 1, n_p-1
            write(6,10) i, p(0,i), p(1,i), p(2,i), p(3,i)
         end do
c
         write(6,*) " "
         write(6,*) "xi(1), xi(2) =", xi
         write(6,*) "mu_f =",sqrt(mufsq(1,1)),sqrt(mufsq(2,1))
      end if

       if (ldoblha) then

         physToDiag(1)=blha_physdiagmap(1,blha_cursubproc,blha_curproc)
         physToDiag(2)=blha_physdiagmap(2,blha_cursubproc,blha_curproc)
         physToDiag(3)=blha_physdiagmap(3,blha_cursubproc,blha_curproc)
         physToDiag(4)=blha_physdiagmap(4,blha_cursubproc,blha_curproc)
         fsign(1) = blha_fsign(1,blha_cursubproc,blha_curproc)
         fsign(2) = blha_fsign(2,blha_cursubproc,blha_curproc)
         fsign(3) = blha_fsign(3,blha_cursubproc,blha_curproc)
         fsign(4) = blha_fsign(4,blha_cursubproc,blha_curproc)

         call wbf_w2j(xi,p,v,physToDiag,fsign,bos, nlo,
     1             pdf,res,nmin,nmax,ps_number)

         return
       endif
c
c call PDF subroutine in order to determine parton
c distributions in the incoming (anti)protons.
c
      q_sf = sqrt(mufsq(1,1))
      call pdfproton( xi(1), q_sf, pdf(-6,1) )        ! f_a(x1)

      q_sf = sqrt(mufsq(2,1))
      call pdfproton( xi(2), q_sf, pdf(-6,2) )        ! f_b(x2)


c select helicity
c      jsig = min(3*rn(1)+1d0,3.01d0)
      jsig = 0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     Now sum over the subprocesses contributing to W production

      nmax = FL_W(fsign,-1)            !reset counter for subprocesses to 0

C*******************  q1 q3 ---> q2 q4  W (Z)  **********************

c   physToDiag(ext.momentum label) = Feynman diagram label

      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=3
      physToDiag(3)=2    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=4
C NOTE: for call of wbf_w2j it is important that p(*,1,*) and p(*,3,*)
c correspond to 1-2 fermion line ALWAYS, i.e physToDiag(1/2)={1,3} and 
c similarly physToDiag(3/4)={2,4} for the 3-4 fermion line
      fsign(1) = 1
      fsign(2) = 1
      fsign(3) = 1
      fsign(4) = 1

 
      call wbf_w2j(xi,p,v,physToDiag,fsign,bos, nlo,
     1             pdf,res,nmin,nmax,ps_number)
!       print*, xi
!       print*, p
!       print*, v
!       print*, bos
!       print*, nmin
!       print*, nmax
!       STOP


c      if (ldebug) call resprint2(nmin,nmax,res)
      if (init.eq.1) nproc(1) = nmax

C******************* q1 qb4 ---> q2 qb3 W  (Z)  **********************
      
      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=4
      physToDiag(3)=2    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=3
c
      fsign(1) =  1
      fsign(2) =  1
      fsign(3) = -1
      fsign(4) = -1
      
      call wbf_w2j(xi,p,v,physToDiag,fsign,bos, nlo,
     1  	  pdf,res,nmin,nmax,ps_number)
c      if (ldebug) call resprint2(nmin,nmax,res)
      if (init.eq.1) nproc(2) = nmax
      
C******************* qbar2 q3 ---> qbar1 q4 W (Z)  **********************
      
      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=3
      physToDiag(3)=1    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=4
c
      fsign(1) = -1
      fsign(2) = -1
      fsign(3) =  1
      fsign(4) =  1
      
      call wbf_w2j(xi,p,v,physToDiag,fsign,bos,nlo,
     1  	  pdf,res,nmin,nmax,ps_number)
c      if (ldebug) call resprint2(nmin,nmax,res)
      if (init.eq.1) nproc(3) = nmax

C*******************  qbar2 qb4 ---> qbar1 qb3 W (Z)  **********************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=4
      physToDiag(3)=1    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=3
c
      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = -1
      fsign(4) = -1
      
      call wbf_w2j(xi,p,v,physToDiag,fsign,bos,nlo,
     1  	  pdf,res,nmin,nmax,ps_number)
c      if (ldebug) call resprint2(nmin,nmax,res)
      if (init.eq.1) nproc(4) = nmax
      
c**************  end of process evaluation part one ******************


! now run again for down-type decay products if 93 93 / 94 94 is 
! requested and the z decays hadronically.

      SELECT CASE(procid)
      CASE(WPZhadjj, WMZhadjj)

       if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
        ! set couplings for Z hadronic decay into ddbar
        call setZtodd

c reset lws to .true. to force recalculation of virtual contributions in qqwwqq
         do i = 1,4
            lws(i) = .true.
         enddo

         j = 1

         if (procid.eq.WPZhadjj) then

         CALL JIOXXX(lp,lm,GZ_ZF ,ZMASS,ZWIDTH,zu)      !Zl
         CALL JIOXXX(lp,lm,GZ_AF ,ZERO ,ZERO  ,au)      !Al
         
         do mu = 0,3
            qe(mu) = v(mu,1)+v(mu,2)
            qu(mu) = v(mu,3)+v(mu,4)
            qvv(mu) = qe(mu) + qu(mu)
         enddo
	 	 
         qe(4) = qe(0)**2-qe(1)**2-qe(2)**2-qe(3)**2
         qu(4) = qu(0)**2-qu(1)**2-qu(2)**2-qu(3)**2
 
	        
	 if (abs(nlo).eq.5 .or. nlo.eq.-4) then     ! determine V/V polarization vectors 
c wet(mu) = we(mu) - xxew*qe(mu) and zut(mu) = zu(mu) - xxmz*qu(mu)
            qq = qe(0)*qvv(0)-qe(1)*qvv(1)-qe(2)*qvv(2)-qe(3)*qvv(3)
            xxew = dotrc(qvv,we)/qq
            
	    qq = qu(0)*qvv(0)-qu(1)*qvv(1)-qu(2)*qvv(2)-qu(3)*qvv(3)
            xxuz = dotrc(qvv,zu)/qq            
            xxua = dotrc(qvv,au)/qq
            do mu = 0,3
               wet(mu+1) = we(mu+1) - xxew*qe(mu)
               zut(mu+1) = zu(mu+1) - xxuz*qu(mu)
               aut(mu+1) = au(mu+1) - xxua*qu(mu)
           enddo
           do mu = 5,6
               wet(mu) = we(mu)
               zut(mu) = zu(mu)
               aut(mu) = au(mu)
           enddo
         elseif (nlo.eq.1) then
            do mu = 1,6
               wet(mu) = we(mu)
               zut(mu) = zu(mu)
               aut(mu) = au(mu)
            enddo
         endif

c for ZW+,AW+ -> e+ ve mu+ mu-:
         call vwtovwp_had(q12(0,j),q34(0,j),v,h,zwwv(0,0,j),awwv(0,0,j)) !q12 = V
         call vwtovwp_had(q34(0,j),q12(0,j),v,h,wzwv(0,0,j),wawv(0,0,j)) !q12 = W+
C for V1V2 -->  mu+ mu- (Vi=A,Z)
         call VVtomm_had(1,h,q12(0,j),v,aauu(0,0,2,j),azuu(0,0,2,j),
     &             zauu(0,0,2,j),zzuu(0,0,2,j))    	 !emit W on lower
         call VVtomm_had(2,h,q34(0,j),v,aauu(0,0,1,j),azuu(0,0,1,j),
     &             zauu(0,0,1,j),zzuu(0,0,1,j))    	 !emit W on upper
C for WW -->  mu+ mu-
         call WWtomm_had(1,h,q34(0,j),v,wwuu(0,0,1,j))     !emit V->e nu on upper
         call WWtomm_had(1,h,q12(0,j),v,wwuu(0,0,2,j))     !emit V->e nu on lower
C for WV --> e+ ve
         call WPVtoln_had(1,q34(0,j),v,awen(0,0,1,j),zwen(0,0,1,j))  !emit V->mu+mu- on upper
         call WPVtoln_had(2,q12(0,j),v,awen(0,0,2,j),zwen(0,0,2,j))  !emit V->mu+mu- on lower
         call WPVtoln_had(1,q12(0,j),v,waen(0,0,2,j),wzen(0,0,2,j))  !emit V->mu+mu- on lower
         call WPVtoln_had(2,q34(0,j),v,waen(0,0,1,j),wzen(0,0,1,j))  !emit V->mu+mu- on upper  
c for W+ --> e+ ve mu+ mu-
         call wpto4l_had(v,h,wwv)

         else
         
         CALL JIOXXX(lp,lm,GZ_ZF ,ZMASS,ZWIDTH,zu)      !Zl
         CALL JIOXXX(lp,lm,GZ_AF ,ZERO ,ZERO  ,au)      !Al	

         do mu = 0,3
            qe(mu) = v(mu,1)+v(mu,2)
            qu(mu) = v(mu,3)+v(mu,4)
            qvv(mu) = qe(mu) + qu(mu)
         enddo
	 	 
         qe(4) = qe(0)**2-qe(1)**2-qe(2)**2-qe(3)**2
         qu(4) = qu(0)**2-qu(1)**2-qu(2)**2-qu(3)**2
 
	 if (abs(nlo).eq.5 .or. nlo.eq.-4) then     ! determine V/V polarization vectors 
c wet(mu) = we(mu) - xxew*qe(mu) and zut(mu) = zu(mu) - xxmz*qu(mu)
            qq = qe(0)*qvv(0)-qe(1)*qvv(1)-qe(2)*qvv(2)-qe(3)*qvv(3)
            xxew = dotrc(qvv,we)/qq
            
	    qq = qu(0)*qvv(0)-qu(1)*qvv(1)-qu(2)*qvv(2)-qu(3)*qvv(3)
            xxuz = dotrc(qvv,zu)/qq            
            xxua = dotrc(qvv,au)/qq
            do mu = 0,3
               wet(mu+1) = we(mu+1) - xxew*qe(mu)
               zut(mu+1) = zu(mu+1) - xxuz*qu(mu)
               aut(mu+1) = au(mu+1) - xxua*qu(mu)
           enddo
           do mu = 5,6
               wet(mu) = we(mu)
               zut(mu) = zu(mu)
               aut(mu) = au(mu)
           enddo
         elseif (nlo.eq.1) then
            do mu = 1,6
               wet(mu) = we(mu)
               zut(mu) = zu(mu)
               aut(mu) = au(mu)
            enddo
         endif

c for ZW-,AW- -> e- ve~ mu+ mu-:
         call vwtovwm_had(q12(0,j),q34(0,j),v,h,zwwv(0,0,j),awwv(0,0,j)) !q12 = V
         call vwtovwm_had(q34(0,j),q12(0,j),v,h,wzwv(0,0,j),wawv(0,0,j)) !q12 = W-
C for V1V2 -->  mu+ mu- (Vi=A,Z)
         call VVtomm_had(1,h,q12(0,j),v,aauu(0,0,2,j),azuu(0,0,2,j),
     #		 zauu(0,0,2,j),zzuu(0,0,2,j))    	 !emit W on lower
         call VVtomm_had(2,h,q34(0,j),v,aauu(0,0,1,j),azuu(0,0,1,j),
     #		 zauu(0,0,1,j),zzuu(0,0,1,j))    	 !emit W on upper
C for W+W- -->  mu+ mu-
         call WWtomm_had(2,h,q34(0,j),v,wwuu(0,0,1,j))     !emit V->e nu on upper
         call WWtomm_had(2,h,q12(0,j),v,wwuu(0,0,2,j))     !emit V->e nu on lower
C for WV --> e- ve~: 
         call WMVtoln_had(1,q34(0,j),v,awen(0,0,1,j),zwen(0,0,1,j))  !emit V->mu+mu- on upper
         call WMVtoln_had(2,q12(0,j),v,awen(0,0,2,j),zwen(0,0,2,j))  !emit V->mu+mu- on lower
         call WMVtoln_had(1,q12(0,j),v,waen(0,0,2,j),wzen(0,0,2,j))  !emit V->mu+mu- on lower
         call WMVtoln_had(2,q34(0,j),v,waen(0,0,1,j),wzen(0,0,1,j))  !emit V->mu+mu- on upper
c for W- --> e- ve mu+ mu-
         call wmto4l_had(v,h,wwv)	
         
         endif


C*******************  q1 q3 ---> q2 q4  W (Z)  **********************

c   physToDiag(ext.momentum label) = Feynman diagram label

      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=3
      physToDiag(3)=2    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=4
C NOTE: for call of wbf_w2j it is important that p(*,1,*) and p(*,3,*)
c correspond to 1-2 fermion line ALWAYS, i.e physToDiag(1/2)={1,3} and 
c similarly physToDiag(3/4)={2,4} for the 3-4 fermion line
      fsign(1) = 1
      fsign(2) = 1
      fsign(3) = 1
      fsign(4) = 1

 
      call wbf_w2j(xi,p,v,physToDiag,fsign,bos, nlo,
     1             pdf,res,nmin,nmax,ps_number)


c      if (ldebug) call resprint2(nmin,nmax,res)
      if (init.eq.1) nproc(1) = nmax

C******************* q1 qb4 ---> q2 qb3 W  (Z)  **********************
      
      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=4
      physToDiag(3)=2    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=3
c
      fsign(1) =  1
      fsign(2) =  1
      fsign(3) = -1
      fsign(4) = -1
      
      call wbf_w2j(xi,p,v,physToDiag,fsign,bos, nlo,
     1  	  pdf,res,nmin,nmax,ps_number)
c      if (ldebug) call resprint2(nmin,nmax,res)
      if (init.eq.1) nproc(2) = nmax
      
C******************* qbar2 q3 ---> qbar1 q4 W (Z)  **********************
      
      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=3
      physToDiag(3)=1    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=4
c
      fsign(1) = -1
      fsign(2) = -1
      fsign(3) =  1
      fsign(4) =  1
      
      call wbf_w2j(xi,p,v,physToDiag,fsign,bos,nlo,
     1  	  pdf,res,nmin,nmax,ps_number)
c      if (ldebug) call resprint2(nmin,nmax,res)
      if (init.eq.1) nproc(3) = nmax

C*******************  qbar2 qb4 ---> qbar1 qb3 W (Z)  **********************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=4
      physToDiag(3)=1    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=3
c
      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = -1
      fsign(4) = -1
      
      call wbf_w2j(xi,p,v,physToDiag,fsign,bos,nlo,
     1  	  pdf,res,nmin,nmax,ps_number)
c      if (ldebug) call resprint2(nmin,nmax,res)
      if (init.eq.1) nproc(4) = nmax


c***  end of extra process evaluation for hadronic decay  ***

       endif
      end select

      if (init.eq.1) then
         nproc(0) = 0
c         print*," process ranges are ",nproc
         init = init+1
      endif

      SELECT CASE(procid)
      CASE(WPhadZjj, WMhadZjj)
        call hadDecayFactor_W(finalquarks(1), abs(mjj2(v(0,1),v(0,2))), N_gen_W, fac_W)
        DO IRES = 1,NMAX
          res(IRES) = res(IRES) * fac_W            ! factor fac_W for hadronic decay
        ENDDO
      CASE(WPZhadjj, WMZhadjj)
        call hadDecayFactor_Z(finalquarks(1), abs(mjj2(v(0,3),v(0,4))), N_gen_up, N_gen_down, fac_Z_up, fac_Z_down)
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
         CASE(WPhadZjj)
            if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
               rnumb=RandomNumber()
               finalquarks_psp(1) =  2 + 2* INT(rnumb*2)
               finalquarks_psp(2) = -1 - 2* INT(rnumb*2)
            endif
         CASE(WMhadZjj)
            if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
               rnumb=RandomNumber()
               finalquarks_psp(1) =  1 + 2* INT(rnumb*2)
               finalquarks_psp(2) = -2 - 2* INT(rnumb*2)
            endif
         CASE(WPZhadjj, WMZhadjj)
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


      if (jsig.ne.0) m2s = m2s*3d0

      if (bos.eq.3.or.bos.eq.4) then
      	m2s_wbfw = m2s
      elseif (bos.eq.32.or.bos.eq.42.or.bos.eq.31.or.bos.eq.41) then
      	m2s_wbfw = m2s*2  	! factor 2 for random helicity summation
      endif
 
      return
      end


c******************************************************************
c
c   begin subroutine wbf_w2j
c
c*****************************************************************
      subroutine wbf_w2j(xi,p,v,physToDiag,fsign,bos, nlo,
     1                   pdf,res,nmin,nmax,ps_number)
          use globalvars, only: ldoblha
      implicit none
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
      real*8 p(0:3,max_p-1), v(0:3,max_v), xi(nx)
      real*8 pdf(-6:6,2), res(maxnumsubproc)
      integer physToDiag(4), fsign(4+max_v), bos, h, nlo, nmin, nmax,ps_number
c      
      common /hcount / h
c
c wbf_w2j calls the amplitude square routines 
c             qqwpqq    for bos = 3 (W+ production)
c             qqwmqq    for bos = 4 (W- production)
c             qqwpzqq   for bos = 32 (W+Z production)
c             qqwmzqq   for bos = 42 (W-Z production)
c             qqwpaqq   for bos = 31 (W+A production)
c             qqwmaqq   for bos = 41 (W-A production)
c
c  INPUT:  p(0:3,4)        external physical parton momenta
c          v(0:3,nv)       Z decay momenta, for Higgs production only
c                          the sum q(mu) = v(mu,1)+...+v(mu,nv) is needed
c          physToDiag(4)   physToDiag(ext.mom. label) = Feynman diagram label
c          fsign           sign factors for amplitude calls; see qqZqq
c          nlo             nlo = 0: calculate LO only
c                          nlo = 1: calculate full NLO (Born + interference)
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
c
      real*8 uucs, ddcs, udcc, udss
      real*8 uusc, ddsc, ducc, duss
      real*8 pbar(0:3,4+max_v), polcol
      double precision matrsq(0:1,0:1,2)
      integer iflav(4), FL_W
      external FL_W

      logical ldebug
      parameter (ldebug=.false.)

      LOGICAL ldebug_ltens
      PARAMETER (ldebug_ltens=.false.)

      integer i,i1,i2,if1,if2,k,lw,mu,qw
      integer bosi
c
      uucs = 0d0
      ddcs = 0d0
      udcc = 0d0
      udss = 0d0
c
      uusc = 0d0
      ddsc = 0d0
      ducc = 0d0
      duss = 0d0

      do if1 = 0,1
         do if2 = 0,1
            do i = 1,2
               matrsq(if1,if2,i) = 0d0
            enddo
         enddo
      enddo
      nmin = nmax+1
      if (bos.lt.10) then
         qw = 7-2*bos           ! W charge, +1 for W+, -1 for W-
      	 bosi = bos
      elseif(bos.eq.32.or.bos.eq.42) then
         qw   =  37/5-bos/5
	 bosi = -1/5+bos/10
      elseif(bos.eq.31.or.bos.eq.41) then
         qw   =  36/5-bos/5
	 bosi = -1/10+bos/10
      endif 
      	
      polcol = 1d0/(4d0*N**2*xi(1)*xi(2))
         
c get the real emission amplitude squared, store it in uucs(1) etc.
      do mu = 0,3
         do i = 1,4
            pbar(mu,physToDiag(i))=p(mu,i)
         enddo
      enddo	 
	 
      if (bos.eq.3.or.bos.eq.4) then
         do mu = 0,3
            pbar(mu,5) = v(mu,1)   ! lbar   for decay
            pbar(mu,6) = v(mu,2)   ! l
 	 enddo
         fsign(5) = -1
         fsign(6) = 1
c
        if (ldebug) then
           print*," ------------------------------------------"
           print*," wbf_w2j info: "
           print*," enter with bos = ",bos," fsign = ",fsign
        endif

        if (bos.eq.3) then                   !the W+ case
           call qqwpqq(pbar,fsign,nlo, uucs,ddcs,udcc,udss)
           if (ldebug) then
              print*," uucs = ",uucs
              print*," ddcs = ",ddcs
              print*," udcc = ",udcc
              print*," udss = ",udss
           endif                       
        elseif (bos.eq.4) then               !the W- case
           call qqwmqq(pbar,fsign,nlo, uusc,ddsc,ducc,duss)
           if (ldebug) then
              print*," uusc = ",uusc
              print*," ddsc = ",ddsc
              print*," ducc = ",ducc
              print*," duss = ",duss
           endif      
        endif !bos
  
     
        elseif(bos.eq.32) then         ! W+Z 
         do mu = 0,3             ! kinematics for H-->W+Z  decay
            pbar(mu,5) = v(mu,1) ! l+
            pbar(mu,6) = v(mu,2) ! nu
            pbar(mu,7) = v(mu,3) ! l"+
            pbar(mu,8) = v(mu,4) ! l"-
         enddo
         fsign(5) = -1
         fsign(6) = 1
         fsign(7) = -1
         fsign(8) = 1

          call qqwpzqq(pbar,fsign, nlo, 1,bos,
     1                    uucs,ddcs,udcc,udss) 

#ifdef WITH_MADGRAPH
         if (ldebug_ltens) call qqwzqq_mg2(p,pbar,v,fsign, nlo,bos,
     1                          uucs,ddcs,udcc,udss) 
#endif

          							 							 
        elseif(bos.eq.42) then         ! W-Z 
         do mu = 0,3             ! kinematics for H-->W-Z decay
            pbar(mu,5) = v(mu,1) ! l-
            pbar(mu,6) = v(mu,2) ! nubar
            pbar(mu,7) = v(mu,3) ! l"+
            pbar(mu,8) = v(mu,4) ! l"-
         enddo
         fsign(5) = 1
         fsign(6) = -1
         fsign(7) = -1
         fsign(8) = 1

         call qqwmzqq(pbar,fsign, nlo, 1,bos,
     1                    uusc,ddsc,ducc,duss) 
     						
#ifdef WITH_MADGRAPH
         if (ldebug_ltens) call qqwzqq_mg2(p,pbar,v,fsign, nlo,bos,
     1                          uusc,ddsc,ducc,duss) 
#endif

      elseif(bos.eq.31) then    !W+A
         do mu = 0,3            ! kinematics for H-->W+Z  decay
            pbar(mu,5) = v(mu,1) ! e+
            pbar(mu,6) = v(mu,2) ! ve
            pbar(mu,7) = v(mu,3) ! a
         enddo
         fsign(5) = -1
         fsign(6) = 1
         fsign(7) = 1
         
         call qqwpaqq(pbar,fsign, nlo, 1,bos,ps_number,
     1        uucs,ddcs,udcc,udss) 

      elseif(bos.eq.41) then    ! W-A 
         do mu = 0,3            ! kinematics for H-->W-Z decay
            pbar(mu,5) = v(mu,1) ! l-
            pbar(mu,6) = v(mu,2) ! nubar
            pbar(mu,7) = v(mu,3) ! a

         enddo
         fsign(5) = 1
         fsign(6) = -1
         fsign(7) = 1

         call qqwmaqq(pbar,fsign, nlo, 1,bos,ps_number,
     1        uusc,ddsc,ducc,duss)
						    
      else
         print*," undefined variable bos in wbf_w2j "
         stop
	
      endif !bos	

      if (ldoblha) then
        return
      endif

c matrsq(i1,i2,lw):  index i1 = flavor of beam 1
c                          i2 = flavor of beam 2
c                    index lw = line to which W is attached (1=12, 2=34 line)
	
	if(bos.eq.3.or.bos.eq.32.or.bos.eq.31) then
         matrsq(0,(1-fsign(4))/2, 2) = uucs ! (1-fs(4))/2 = flavor id
         matrsq(1,(1-fsign(4))/2, 2) = ddcs ! of beam 2 (= 0,1)
         matrsq((1-fsign(2))/2,0, 1) = udcc ! 1,j id"s W attached
         matrsq((1-fsign(2))/2,1, 1) = udss ! to 12 line, g to j line
	
	elseif(bos.eq.4.or.bos.eq.42.or.bos.eq.41) then
         matrsq(0,(1+fsign(4))/2, 2) = uusc ! (1-fs(3))/2 = flavor id
         matrsq(1,(1+fsign(4))/2, 2) = ddsc ! of beam 2 (= 0,1)
         matrsq((1+fsign(2))/2,0, 1) = ducc ! 1,j id"s W attached
         matrsq((1+fsign(2))/2,1, 1) = duss ! to 12 line, g to j line
	
	endif !bos 
      
         
      do lw = 1,2               ! identify line W is attached to
         do if1 = 1,(nfl/2)*2   ! identify flavor of NC vertex
            iflav(3-lw)=if1*fsign(physToDiag(3-lw))
            iflav(5-lw)=if1*fsign(physToDiag(5-lw))
            do if2=bosi/2+(1+fsign(2*lw))/2*qw,(nfl/2)*2,2
               iflav(lw)=if2*fsign(physToDiag(lw))
               iflav(lw+2)=if2*fsign(physToDiag(lw+2))-qw
               i1 = mod(abs(iflav(1)),2)
               i2 = mod(abs(iflav(2)),2)
               k=FL_W(iflav,1)
               res(k) = pdf(sign1*iflav(1),1)      !1for j=2;3for j=3
     &                 *pdf(sign2*iflav(2),2) !2for j=2;1for j=3
     &                 *matrsq(i1,i2,lw)*polcol
     		
               if (ldebug) then                                          !debug
                  print*," iflav = ",iflav(1),iflav(2)," --> ",          !debug
     &                               iflav(3),iflav(4)                   !debug
                  print*," called matrix element ",i1,i2,lw              !debug
               endif
            enddo               ! if2
         enddo                  ! if1
      enddo                     ! lw

      nmax = FL_W(iflav,-2)
            
      
      end
c******************************************************************
c
c   end subroutine wbf_w2j
c
c*****************************************************************


c******************************************************************
c
c   begin function FL_W
c
c*****************************************************************
      integer function FL_W(iflav,colstruc)
      implicit none
      integer iflav(4),colstruc    ! input for color and flavor assignment
c  fill hepup color and flavor assignments for WBF processes with
c  two quark lines only

c  iflav1,2 are flavor identifiers for incoming quarks
c  iflav3,4 are flavor identifiers for outgoing quarks
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

      if (colstruc.eq.-1) then
c       write(*,*) "we are resetting"
         listposition=0
         numresets=numresets+1
      else if (colstruc.ne.-2) then
c  don"t do anything for colstruc=-2, this call simply 
c  returns a value for numcolor. Otherwise: increment 
c  the counter regardless of whether or not it"s the 1st 
c  time through.
         listposition=listposition+1

c  and fill info for this subprocess during first pass 
         if(numresets.eq.0) then
            call fillColoredPartons_W(iflav,
     &                                listposition+numdecay)
         endif        ! numresets.eq.0
      endif
      FL_W=listposition
      end
c*****************************************************************************
c
c   end function FL_W
c
c*****************************************************************************

c*****************************************************************************
c
c    begin  subroutine fillColoredPartons_W
c
c*********************************************************************
      subroutine fillColoredPartons_W(iflav,listposition)
c  assigns values to the variables in the common block localHEPUP
c  in particular, this subroutine assigns values to those variables that 
c  will be stored in the lookup tables generated by writeHEPUPtable.  
c  As the name suggests, this routine only stores the information for the 
c  colored partons.  Particles without color will be dealt with in the 
c  subroutine fillColorless.
      IMPLICIT NONE
c  include files
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/process.inc"

c     include "qqh.inc"
c  input parameters: 4 quark flavors and subprocess number
      integer iflav(4),id1,id2,id3,id4,listposition
c  local variables
      integer i

      id1=iflav(1)
      id2=iflav(2)
      id3=iflav(3)
      id4=iflav(4)

      select case(procID)
      case(WPZjj,WMZjj,WPhadZjj,WPZhadjj,WMhadZjj,WMZhadjj) !W+Zjj, W-Zjj
         lnup(listposition)=numParticles+2
      case(WPAjj,WMAjj) !W+Ajj, W-Ajj
         lnup(listposition)=numParticles
       case(WPjj,WMjj) ! W+jj, W-jj
         lnup(listposition)=numParticles-1
      end select
            
c  lidprup is same for all subprocesses

      listup(1,listposition)=-1   !incoming quarks
      listup(2,listposition)=-1
      listup(3,listposition)=1    !outgoing quarks
      listup(4,listposition)=1
c  listup etc. for W and decay leptons is set in subroutine fillColorless
c  which is called in decay setup

      lidup(1,listposition)=id1     ! the 4 quarks
      lidup(2,listposition)=id2
      lidup(3,listposition)=id3
      lidup(4,listposition)=id4

      do i=3,n_p
        lmothup(1,i,listposition)=1
        lmothup(2,i,listposition)=2
      enddo

      if (id1.gt.0) then     ! 1 and 3 are quarks
        licolup(1,1,listposition)=501
        licolup(1,3,listposition)=501
        licolup(2,1,listposition)=0
        licolup(2,3,listposition)=0
      else                   ! 1 and 3 are anti-quarks
        licolup(1,1,listposition)=0
        licolup(1,3,listposition)=0
        licolup(2,1,listposition)=501
        licolup(2,3,listposition)=501
      endif

      if (id2.gt.0) then     ! 2 and 4 are quarks
        licolup(1,2,listposition)=502
        licolup(1,4,listposition)=502
        licolup(2,2,listposition)=0
        licolup(2,4,listposition)=0
      else                   ! 2 and 4 are anti-quarks
        licolup(1,2,listposition)=0
        licolup(1,4,listposition)=0
        licolup(2,2,listposition)=502
        licolup(2,4,listposition)=502
      endif

      end  ! fillColoredPartons_Wg
c*****************************************************************************
c
c    end subroutine fillColoredPartons_Wg
c
c*****************************************************************************

      subroutine vtowz_reset
      implicit none
#include "tensorw.inc"
      complex*16 zero
      parameter (zero=(0d0,0d0))
      integer j,jj,mu,nu

      do mu = 0,3
	 do nu = 0,3
            do j = 1,3
               zwwv(mu,nu,j) = zero
               awwv(mu,nu,j) = zero
	       wzwv(mu,nu,j) = zero
	       wawv(mu,nu,j) = zero	       
               wztowa(mu,nu,j) = zero
               watowa(mu,nu,j) = zero
               awtowa(mu,nu,j) = zero
               zwtowa(mu,nu,j) = zero
	       
	       do jj = 1,2	       
     	       aauu(mu,nu,jj,j) = zero
	       azuu(mu,nu,jj,j) = zero
	       zauu(mu,nu,jj,j) = zero
	       zzuu(mu,nu,jj,j) = zero
	       wwuu(mu,nu,jj,j) = zero
	       awen(mu,nu,jj,j) = zero
	       zwen(mu,nu,jj,j) = zero
	       waen(mu,nu,jj,j) = zero
	       wzen(mu,nu,jj,j) = zero	   
               wwtoa(mu,nu,jj,j) = zero     
	       enddo !jj

            enddo !j
         enddo !nu
      enddo !mu
      
      do mu = 0,5
         wwv(mu) = 0d0
         wtowa(mu) = 0d0
      enddo !mu
      
      print*," v(v) to WZ/WA tensors initialized to ",zero 
      return
      end

