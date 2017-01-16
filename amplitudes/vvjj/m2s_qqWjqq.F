c
c	Dieter Zeppenfeld, <dieter@pheno.physics.wisc.edu>
c       Rob Szalapski, <robs@theory.kek.jp>
c	Initial version:  1995 March 24
c	Last modified:  
c
      subroutine m2s_wbfw3j(
     &                   bos,   !in:  Boson identifier, 2 = Z, 6=H
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

      integer N_gen_W
      real*8 fac_W
      integer N_gen_up, N_gen_down
      real*8 fac_Z_up, fac_Z_down, mjj2
      external mjj2

      integer bos
      real*8 xi(nx), p(0:3,max_p,3), v(0:3,max_v), 
     1       rn(1), xuz(2,2:3), m2s(0:3)
      logical lok(3)

c  helicity selection
      integer jsig, jsig1, jsig3, jsig5
      common /chelsum/ jsig,jsig1,jsig3,jsig5
c
      complex*16 zero
      parameter (zero=(0d0,0d0))
c
c declare external functions
c
      integer FL_Wg
      external FL_Wg
c
c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"
      real*8 x1,x2
      real*8 q12(0:4,3), q34(0:4,3),qvv(0:3), lnQomu(2:3), omxi(2:3), 
     1       Ax(2:3), Bx(2:3), Cx(2:3), Dxi(2:3), tgs2oqsq(2:3), 
     2       ln1mxi, lnrat, z, lnz
c
c for check only:


c declare local variables
c
c
      real*8 q_sf
      real*8 qa3,qa5,q35,rjj, rsep, dot0p, qmina
      external dot0p,rsep

      integer init/0/, I, J, nproc(8),mu
      save init, nproc
c
c declare variables for summation process
c
      INTEGER ires, nmaxold
      save nmaxold

      integer FSIGN(4+max_v), gsign, physToDiag(5), nlo, nmin, nmax,h
      save nmin, nmax
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
c
c define program switches
c
      logical ldebug,ldebug_sub

      parameter (ldebug=.false.)
      data ldebug_sub /.false./		!output debug information

      common /hcount / h
      real*8 weight,rnumb,RandomNumber


c -----------------------------------------------------------------------
c
c if first iteration, output acceptance cuts and fix constant input 
c parameters for call of qqZWqqji.f
c
      if ( init .ne. 0 ) then
         do i = 1,nmax
            do j=0,3
               res(i,j) = 0
            enddo
         enddo
      else if (ldoblha) then
        h = 1
        init = 1
      else
      
	if (bos.eq.3.or.bos.eq.4) then    	  
         write(6,*) " "
         write(6,*) "      Wjjj amplitude square information:"
         write(6,*) "-----------------------------------------------"
         write(6,*) " "
	elseif (bos.eq.32) then      
         write(6,*) " "
         write(6,*) "      W+Zjjj amplitude square information:"
         write(6,*) "-----------------------------------------------"
         write(6,*) " "
         call vtowz_reset
	elseif (bos.eq.42) then      
         write(6,*) " "
         write(6,*) "      W-Zjjj amplitude square information:"
         write(6,*) "-----------------------------------------------"
         write(6,*) " "
         call vtowz_reset
      elseif (bos.eq.31) then      
         write(6,*) " "
         write(6,*) "      W+Ajjj amplitude square information:"
         write(6,*) "-----------------------------------------------"
         write(6,*) " "
         call vtowz_reset
      elseif (bos.eq.41) then      
         write(6,*) " "
         write(6,*) "      W-Ajjj amplitude square information:"
         write(6,*) "-----------------------------------------------"
         write(6,*) " "
         write(6,*) " "
         call vtowz_reset
	else
	 write(6,*) " this process is not implemented in m2sWj"
	 stop
	endif !bos   

         call printnfl(.true.)

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
      endif
      
c     -----------------------------------------------------------------------
      call Calc_Momentum_Transfer(p, v, q12,q34,3)


c     setting helicities for W+jj, W-jj
      if (bos .eq. 3) then
         helicity(1) = -1       ! nu_e
         helicity(2) =  1       ! e+
      else if (bos .eq. 4) then
         helicity(1) =  1       ! nu_e_bar
         helicity(2) = -1       ! e-
      end if
      
      if (bos.eq.32) then       ! W+Z
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
c     fix hel. for checks:
c     h = 2
c     print*,"h= ",h

         ie = 1
         iu = sign(1,1-h)


c...  Les Houches interface         
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
!     for general output to all combinations: up-type first
            if ((finalquarks(1).eq.93 .or. finalquarks(1).eq.94)) then
!     set couplings for Z hadronic decay into uubar
               call setZtouu
            endif
            CALL JIOXXX(lp,lm,GZ_ZF ,ZMASS,ZWIDTH,zu) !Zl
            CALL JIOXXX(lp,lm,GZ_AF ,ZERO ,ZERO  ,au) !Al
         CASE DEFAULT
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
         
c     leptonic tensors:

         SELECT CASE(procid)
         CASE(WPhadZjj,WPZhadjj)


            do j = 2,3
               call anomal_formfactor(q12(0,j),q34(0,j),qe(0),qu(0))
c     for ZW+,AW+ -> e+ ve mu+ mu-:
               call vwtovwp_had(q12(0,j),q34(0,j),v,h,zwwv(0,0,j),awwv(0,0,j)) !q12 = V
               call vwtovwp_had(q34(0,j),q12(0,j),v,h,wzwv(0,0,j),wawv(0,0,j)) !q12 = W+
C     for V1V2 -->  mu+ mu- (Vi=A,Z)
               call VVtomm_had(1,h,q12(0,j),v,aauu(0,0,2,j),azuu(0,0,2,j),
     $              zauu(0,0,2,j),zzuu(0,0,2,j)) !emit W on lower
               call VVtomm_had(2,h,q34(0,j),v,aauu(0,0,1,j),azuu(0,0,1,j),
     $              zauu(0,0,1,j),zzuu(0,0,1,j)) !emit W on upper
C     for WW -->  mu+ mu-
               call WWtomm_had(1,h,q34(0,j),v,wwuu(0,0,1,j)) !emit V->e nu on upper
               call WWtomm_had(1,h,q12(0,j),v,wwuu(0,0,2,j)) !emit V->e nu on lower
C     for WV --> e+ ve
               call WPVtoln_had(1,q34(0,j),v,awen(0,0,1,j),zwen(0,0,1,j)) !emit V->mu+mu- on upper
               call WPVtoln_had(2,q12(0,j),v,awen(0,0,2,j),zwen(0,0,2,j)) !emit V->mu+mu- on lower
               call WPVtoln_had(1,q12(0,j),v,waen(0,0,2,j),wzen(0,0,2,j)) !emit V->mu+mu- on lower
               call WPVtoln_had(2,q34(0,j),v,waen(0,0,1,j),wzen(0,0,1,j)) !emit V->mu+mu- on upper
            enddo               !j
c     for W+ --> e+ ve mu+ mu-
            call wpto4l_had(v,h,wwv)	


         CASE DEFAULT

            if (with_kk) then
#ifdef WITH_KK
               do j = 2,3
      
c     for ZW+,AW+ -> e+ ve mu+ mu-:
                  call vwtovwp_kk(q12(0,j),q34(0,j),v,h,zwwv(0,0,j),awwv(0,0,j)) !q12 = V
                  call vwtovwp_kk(q34(0,j),q12(0,j),v,h,wzwv(0,0,j),wawv(0,0,j)) !q12 = W+
           
C     for V1V2 -->  mu+ mu- (Vi=A,Z)
                  call VVtomm_kk(1,h,q12(0,j),v,aauu(0,0,2,j),azuu(0,0,2,j),
     $                 zauu(0,0,2,j),zzuu(0,0,2,j)) !emit W on lower
                  call VVtomm_kk(2,h,q34(0,j),v,aauu(0,0,1,j),azuu(0,0,1,j),
     $                 zauu(0,0,1,j),zzuu(0,0,1,j)) !emit W on upper

C     for WW -->  mu+ mu-
                  call WWtomm_kk(1,h,q34(0,j),v,wwuu(0,0,1,j)) !emit V->e nu on upper
                  call WWtomm_kk(1,h,q12(0,j),v,wwuu(0,0,2,j)) !emit V->e nu on lower

C     for WV --> e+ ve
                  call WPVtoln_kk(1,q34(0,j),v,awen(0,0,1,j),zwen(0,0,1,j)) !emit V->mu+mu- on upper
                  call WPVtoln_kk(2,q12(0,j),v,awen(0,0,2,j),zwen(0,0,2,j)) !emit V->mu+mu- on lower
c     
                  call WPVtoln_kk(1,q12(0,j),v,waen(0,0,2,j),wzen(0,0,2,j)) !emit V->mu+mu- on lower
                  call WPVtoln_kk(2,q34(0,j),v,waen(0,0,1,j),wzen(0,0,1,j)) !emit V->mu+mu- on upper
                  
               enddo            !j
c     for W+ --> e+ ve mu+ mu-
	
               call wpto4l_kk(v,h,wwv)	
#endif

            elseif (with_spin2) then    
               do j = 2,3
                  
c     for ZW+,AW+ -> e+ ve mu+ mu-:
                  call vwtovwp_spin2(q12(0,j),q34(0,j),v,h,zwwv(0,0,j),awwv(0,0,j)) !q12 = V
                  call vwtovwp_spin2(q34(0,j),q12(0,j),v,h,wzwv(0,0,j),wawv(0,0,j)) !q12 = W+

c     other processes do not change for spin2            
C     for V1V2 -->  mu+ mu- (Vi=A,Z)
                  call VVtomm(1,h,q12(0,j),v,aauu(0,0,2,j),azuu(0,0,2,j),
     $                 zauu(0,0,2,j),zzuu(0,0,2,j)) !emit W on lower
                  call VVtomm(2,h,q34(0,j),v,aauu(0,0,1,j),azuu(0,0,1,j),
     $                 zauu(0,0,1,j),zzuu(0,0,1,j)) !emit W on upper

C     for WW -->  mu+ mu-
                  call WWtomm(1,h,q34(0,j),v,wwuu(0,0,1,j)) !emit V->e nu on upper
                  call WWtomm(1,h,q12(0,j),v,wwuu(0,0,2,j)) !emit V->e nu on lower
                  
C     for WV --> e+ ve
                  call WPVtoln(1,q34(0,j),v,awen(0,0,1,j),zwen(0,0,1,j)) !emit V->mu+mu- on upper
                  call WPVtoln(2,q12(0,j),v,awen(0,0,2,j),zwen(0,0,2,j)) !emit V->mu+mu- on lower
c     
                  call WPVtoln(1,q12(0,j),v,waen(0,0,2,j),wzen(0,0,2,j)) !emit V->mu+mu- on lower
                  call WPVtoln(2,q34(0,j),v,waen(0,0,1,j),wzen(0,0,1,j)) !emit V->mu+mu- on upper
  
               enddo            !j
               
c     for W+ --> e+ ve mu+ mu-
	
               call wpto4l(v,h,wwv)	
               

            elseif (with_anom) then ! anomalous gauge boson couplings
c     using global form factor for all tensors of one phase space point
c     this ensures proper cancellations for anomalous contributions
c     energy scale is invariant WZ mass        
               
               do j = 2,3
                  
                  call anomal_formfactor(q12(0,j),q34(0,j),qe(0),qu(0))
                  
c     for ZW+,AW+ -> e+ ve mu+ mu-:
                  call vwtovwp_anomal(q12(0,j),q34(0,j),v,h,zwwv(0,0,j),awwv(0,0,j)) !q12 = V
                  call vwtovwp_anomal(q34(0,j),q12(0,j),v,h,wzwv(0,0,j),wawv(0,0,j)) !q12 = W+
           
C     for V1V2 -->  mu+ mu- (Vi=A,Z)
                  call VVtomm_anomal(1,h,q12(0,j),v,aauu(0,0,2,j),azuu(0,0,2,j),
     $                 zauu(0,0,2,j),zzuu(0,0,2,j)) !emit W on lower
                  call VVtomm_anomal(2,h,q34(0,j),v,aauu(0,0,1,j),azuu(0,0,1,j),
     $                 zauu(0,0,1,j),zzuu(0,0,1,j)) !emit W on upper
                  
C     for WW -->  mu+ mu-
                  call WWtomm_anomal(1,h,q34(0,j),v,wwuu(0,0,1,j)) !emit V->e nu on upper
                  call WWtomm_anomal(1,h,q12(0,j),v,wwuu(0,0,2,j)) !emit V->e nu on lower
                  
C     for WV --> e+ ve
                  call WPVtoln_anomal(1,q34(0,j),v,awen(0,0,1,j),zwen(0,0,1,j)) !emit V->mu+mu- on upper
                  call WPVtoln_anomal(2,q12(0,j),v,awen(0,0,2,j),zwen(0,0,2,j)) !emit V->mu+mu- on lower
c     
                  call WPVtoln_anomal(1,q12(0,j),v,waen(0,0,2,j),wzen(0,0,2,j)) !emit V->mu+mu- on lower
                  call WPVtoln_anomal(2,q34(0,j),v,waen(0,0,1,j),wzen(0,0,1,j)) !emit V->mu+mu- on upper
                  
               enddo            !j
	 
c     for W+ --> e+ ve mu+ mu-
               
               call wpto4l_anomal(v,h,wwv)	
               
               
            else                !SM
	     
               do j = 2,3
                  
c     for ZW+,AW+ -> e+ ve mu+ mu-:
                  call vwtovwp(q12(0,j),q34(0,j),v,h,zwwv(0,0,j),awwv(0,0,j)) !q12 = V
                  call vwtovwp(q34(0,j),q12(0,j),v,h,wzwv(0,0,j),wawv(0,0,j)) !q12 = W+
           
C     for V1V2 -->  mu+ mu- (Vi=A,Z)
                  call VVtomm(1,h,q12(0,j),v,aauu(0,0,2,j),azuu(0,0,2,j),
     $                 zauu(0,0,2,j),zzuu(0,0,2,j)) !emit W on lower
                  call VVtomm(2,h,q34(0,j),v,aauu(0,0,1,j),azuu(0,0,1,j),
     $                 zauu(0,0,1,j),zzuu(0,0,1,j)) !emit W on upper
        
C     for WW -->  mu+ mu-
                  call WWtomm(1,h,q34(0,j),v,wwuu(0,0,1,j)) !emit V->e nu on upper
                  call WWtomm(1,h,q12(0,j),v,wwuu(0,0,2,j)) !emit V->e nu on lower
                  
C     for WV --> e+ ve
                  call WPVtoln(1,q34(0,j),v,awen(0,0,1,j),zwen(0,0,1,j)) !emit V->mu+mu- on upper
                  call WPVtoln(2,q12(0,j),v,awen(0,0,2,j),zwen(0,0,2,j)) !emit V->mu+mu- on lower
c     
                  call WPVtoln(1,q12(0,j),v,waen(0,0,2,j),wzen(0,0,2,j)) !emit V->mu+mu- on lower
                  call WPVtoln(2,q34(0,j),v,waen(0,0,1,j),wzen(0,0,1,j)) !emit V->mu+mu- on upper
        
               enddo            !j
      
c     for W+ --> e+ ve mu+ mu-
      
               call wpto4l(v,h,wwv)	
            endif
        
         end select             !semilep
      
c -----------------------------------------------------------------------
      
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
! for general output to all combinations: up-type first
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
         
         
c     leptonic tensors:


         SELECT CASE(procid)
         CASE(WMhadZjj,WMZhadjj)


            do j = 2,3

c     using global form factor for all tensors of one phase space point
c     this ensures proper cancellations for anomalous contributions
c     energy scale is invariant WZ mass
               call anomal_formfactor(q12(0,j),q34(0,j),qe(0),qu(0))
               call vwtovwm_had(q12(0,j),q34(0,j),v,h,zwwv(0,0,j),awwv(0,0,j)) !q12 = V
               call vwtovwm_had(q34(0,j),q12(0,j),v,h,wzwv(0,0,j),wawv(0,0,j)) !q12 = W-
C     for V1V2 -->  mu+ mu- (Vi=A,Z)
               call VVtomm_had(1,h,q12(0,j),v,aauu(0,0,2,j),azuu(0,0,2,j),
     $              zauu(0,0,2,j),zzuu(0,0,2,j)) !emit W on lower
               call VVtomm_had(2,h,q34(0,j),v,aauu(0,0,1,j),azuu(0,0,1,j),
     $              zauu(0,0,1,j),zzuu(0,0,1,j)) !emit W on upper
C     for W+W- -->  mu+ mu-
               call WWtomm_had(2,h,q34(0,j),v,wwuu(0,0,1,j)) !emit V->e nu on upper
               call WWtomm_had(2,h,q12(0,j),v,wwuu(0,0,2,j)) !emit V->e nu on lower
C     for WV --> e- ve~: 
               call WMVtoln_had(1,q34(0,j),v,awen(0,0,1,j),zwen(0,0,1,j)) !emit V->mu+mu- on upper
               call WMVtoln_had(2,q12(0,j),v,awen(0,0,2,j),zwen(0,0,2,j)) !emit V->mu+mu- on lower
               call WMVtoln_had(1,q12(0,j),v,waen(0,0,2,j),wzen(0,0,2,j)) !emit V->mu+mu- on lower
               call WMVtoln_had(2,q34(0,j),v,waen(0,0,1,j),wzen(0,0,1,j)) !emit V->mu+mu- on upper
            enddo               !j		 
c     for W- --> e- ve mu+ mu-
            call wmto4l_had(v,h,wwv)	


         CASE DEFAULT
            
            if (with_kk) then
#ifdef WITH_KK
               do j = 2,3

c     for ZW-,AW- -> e- ve~ mu+ mu-:
                  call vwtovwm_kk(q12(0,j),q34(0,j),v,h,zwwv(0,0,j),awwv(0,0,j)) !q12 = V
                  call vwtovwm_kk(q34(0,j),q12(0,j),v,h,wzwv(0,0,j),wawv(0,0,j)) !q12 = W-
                  
                  
C     for V1V2 -->  mu+ mu- (Vi=A,Z)
                  call VVtomm_kk(1,h,q12(0,j),v,aauu(0,0,2,j),azuu(0,0,2,j),
     $                 zauu(0,0,2,j),zzuu(0,0,2,j)) !emit W on lower
                  call VVtomm_kk(2,h,q34(0,j),v,aauu(0,0,1,j),azuu(0,0,1,j),
     $                 zauu(0,0,1,j),zzuu(0,0,1,j)) !emit W on upper

C     for W+W- -->  mu+ mu-
                  call WWtomm_kk(2,h,q34(0,j),v,wwuu(0,0,1,j)) !emit V->e nu on upper
                  call WWtomm_kk(2,h,q12(0,j),v,wwuu(0,0,2,j)) !emit V->e nu on lower

C     for WV --> e- ve~: 
                  call WMVtoln_kk(1,q34(0,j),v,awen(0,0,1,j),zwen(0,0,1,j)) !emit V->mu+mu- on upper
                  call WMVtoln_kk(2,q12(0,j),v,awen(0,0,2,j),zwen(0,0,2,j)) !emit V->mu+mu- on lower
c     
                  call WMVtoln_kk(1,q12(0,j),v,waen(0,0,2,j),wzen(0,0,2,j)) !emit V->mu+mu- on lower
                  call WMVtoln_kk(2,q34(0,j),v,waen(0,0,1,j),wzen(0,0,1,j)) !emit V->mu+mu- on upper
                  
               enddo            !j		 
#endif		 
c     for W- --> e- ve mu+ mu-
               
               call wmto4l(v,h,wwv)	

            elseif (with_spin2) then    
               do j = 2,3

c     for ZW-,AW- -> e- ve~ mu+ mu-:
                  call vwtovwm_spin2(q12(0,j),q34(0,j),v,h,zwwv(0,0,j),awwv(0,0,j)) !q12 = V
                  call vwtovwm_spin2(q34(0,j),q12(0,j),v,h,wzwv(0,0,j),wawv(0,0,j)) !q12 = W-
                  
c     other processes do not change for spin2
                  
C     for V1V2 -->  mu+ mu- (Vi=A,Z)
                  call VVtomm(1,h,q12(0,j),v,aauu(0,0,2,j),azuu(0,0,2,j),
     $                 zauu(0,0,2,j),zzuu(0,0,2,j)) !emit W on lower
                  call VVtomm(2,h,q34(0,j),v,aauu(0,0,1,j),azuu(0,0,1,j),
     $                 zauu(0,0,1,j),zzuu(0,0,1,j)) !emit W on upper

C     for W+W- -->  mu+ mu-
                  call WWtomm(2,h,q34(0,j),v,wwuu(0,0,1,j)) !emit V->e nu on upper
                  call WWtomm(2,h,q12(0,j),v,wwuu(0,0,2,j)) !emit V->e nu on lower

C     for WV --> e- ve~: 
                  call WMVtoln(1,q34(0,j),v,awen(0,0,1,j),zwen(0,0,1,j)) !emit V->mu+mu- on upper
                  call WMVtoln(2,q12(0,j),v,awen(0,0,2,j),zwen(0,0,2,j)) !emit V->mu+mu- on lower
c     
                  call WMVtoln(1,q12(0,j),v,waen(0,0,2,j),wzen(0,0,2,j)) !emit V->mu+mu- on lower
                  call WMVtoln(2,q34(0,j),v,waen(0,0,1,j),wzen(0,0,1,j)) !emit V->mu+mu- on upper
                  
               enddo            !j		 
               
c     for W- --> e- ve mu+ mu-
               
               call wmto4l(v,h,wwv)	
               
               
            elseif (with_anom) then ! anomalous gauge boson couplings

               do j = 2,3

c     using global form factor for all tensors of one phase space point
c     this ensures proper cancellations for anomalous contributions
c     energy scale is invariant WZ mass
                  call anomal_formfactor(q12(0,j),q34(0,j),qe(0),qu(0))

                  call vwtovwm_anomal(q12(0,j),q34(0,j),v,h,zwwv(0,0,j),awwv(0,0,j)) !q12 = V
                  call vwtovwm_anomal(q34(0,j),q12(0,j),v,h,wzwv(0,0,j),wawv(0,0,j)) !q12 = W-
                  
                  
C     for V1V2 -->  mu+ mu- (Vi=A,Z)
                  call VVtomm_anomal(1,h,q12(0,j),v,aauu(0,0,2,j),azuu(0,0,2,j),
     $                 zauu(0,0,2,j),zzuu(0,0,2,j)) !emit W on lower
                  call VVtomm_anomal(2,h,q34(0,j),v,aauu(0,0,1,j),azuu(0,0,1,j),
     $                 zauu(0,0,1,j),zzuu(0,0,1,j)) !emit W on upper

C     for W+W- -->  mu+ mu-
                  call WWtomm_anomal(2,h,q34(0,j),v,wwuu(0,0,1,j)) !emit V->e nu on upper
                  call WWtomm_anomal(2,h,q12(0,j),v,wwuu(0,0,2,j)) !emit V->e nu on lower

C     for WV --> e- ve~: 
                  call WMVtoln_anomal(1,q34(0,j),v,awen(0,0,1,j),zwen(0,0,1,j)) !emit V->mu+mu- on upper
                  call WMVtoln_anomal(2,q12(0,j),v,awen(0,0,2,j),zwen(0,0,2,j)) !emit V->mu+mu- on lower
c     
                  call WMVtoln_anomal(1,q12(0,j),v,waen(0,0,2,j),wzen(0,0,2,j)) !emit V->mu+mu- on lower
                  call WMVtoln_anomal(2,q34(0,j),v,waen(0,0,1,j),wzen(0,0,1,j)) !emit V->mu+mu- on upper
                  
               enddo            !j		 
               
c     for W- --> e- ve mu+ mu-
               
               call wmto4l_anomal(v,h,wwv)	


            else                !SM

               do j = 2,3

                  call vwtovwm(q12(0,j),q34(0,j),v,h,zwwv(0,0,j),awwv(0,0,j)) !q12 = V
                  call vwtovwm(q34(0,j),q12(0,j),v,h,wzwv(0,0,j),wawv(0,0,j)) !q12 = W-
                  
                  
C     for V1V2 -->  mu+ mu- (Vi=A,Z)
                  call VVtomm(1,h,q12(0,j),v,aauu(0,0,2,j),azuu(0,0,2,j),
     $                 zauu(0,0,2,j),zzuu(0,0,2,j)) !emit W on lower
                  call VVtomm(2,h,q34(0,j),v,aauu(0,0,1,j),azuu(0,0,1,j),
     $                 zauu(0,0,1,j),zzuu(0,0,1,j)) !emit W on upper

C     for W+W- -->  mu+ mu-
                  call WWtomm(2,h,q34(0,j),v,wwuu(0,0,1,j)) !emit V->e nu on upper
                  call WWtomm(2,h,q12(0,j),v,wwuu(0,0,2,j)) !emit V->e nu on lower

C     for WV --> e- ve~: 
                  call WMVtoln(1,q34(0,j),v,awen(0,0,1,j),zwen(0,0,1,j)) !emit V->mu+mu- on upper
                  call WMVtoln(2,q12(0,j),v,awen(0,0,2,j),zwen(0,0,2,j)) !emit V->mu+mu- on lower
c     
                  call WMVtoln(1,q12(0,j),v,waen(0,0,2,j),wzen(0,0,2,j)) !emit V->mu+mu- on lower
                  call WMVtoln(2,q34(0,j),v,waen(0,0,1,j),wzen(0,0,1,j)) !emit V->mu+mu- on upper
                  
               enddo            !j		 
               
c     for W- --> e- ve mu+ mu-
               
               call wmto4l(v,h,wwv)	

            endif
            
         end select             !semilep

      elseif (bos.eq.31) then   ! W+A
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
            helicity(1)= ie
            helicity(2)=-ie
            helicity(3)= iu 
         endif

         CALL IXXXXX(v(0,1),ZERO ,+ie,-1,wl) !e+	
         CALL OXXXXX(v(0,2),ZERO ,-ie, 1,wn) !ve 
         CALL VXXXXX(v(0,3),ZERO ,+iu, 1,ea) !a

         CALL JIOXXX(wl,wn,GWF ,WMASS,WWIDTH,we) !Wl

         do mu = 0,3
            qe(mu) = v(mu,1)+v(mu,2)
            qa(mu) = v(mu,3)
            qvv(mu) = qe(mu) + qa(mu)
         enddo
         
         qe(4) = qe(0)**2-qe(1)**2-qe(2)**2-qe(3)**2
         qa(4) = qa(0)**2-qa(1)**2-qa(2)**2-qa(3)**2

c     leptonic tensors:


         if (with_anom) then ! anomalous gauge boson couplings
c     using global form factor for all tensors of one phase space point
c     this ensures proper cancellations for anomalous contributions
c     energy scale is invariant WA mass        
            do j = 2,3
               call anomal_formfactor(q12(0,j),q34(0,j),qe(0),qa(0))

c     for VW+ -> e+ ve A:
               call vwptowpa_anomal(q12(0,j),q34(0,j),v,h,zwtowa(0,0,j),awtowa(0,0,j)) !q12 = V
               call vwptowpa_anomal(q34(0,j),q12(0,j),v,h,wztowa(0,0,j),watowa(0,0,j)) !q12 = W+
               
C     for W+W- -->  A
               call calc_wwtoa_anomal(1,iu,q34(0,j),v,wwtoa(0,0,1,j)) !emit W+ -> e+ nu on upper
               call calc_wwtoa_anomal(1,iu,q12(0,j),v,wwtoa(0,0,2,j)) !emit W+ -> e+ nu on lower

C     for WV --> e+ ve
               call WPVtoln_anomal(1,q34(0,j),v,awen(0,0,1,j),zwen(0,0,1,j)) !emit A on upper
               call WPVtoln_anomal(2,q12(0,j),v,awen(0,0,2,j),zwen(0,0,2,j)) !emit A on lower
c     
               call WPVtoln_anomal(1,q12(0,j),v,waen(0,0,2,j),wzen(0,0,2,j)) !emit A on lower
               call WPVtoln_anomal(2,q34(0,j),v,waen(0,0,1,j),wzen(0,0,1,j)) !emit A on upper
               
            enddo               !j
c     for W+ --> e+ ve A
            call wptowpa_anomal(qvv,v,h,wtowa)

         else                   !SM
            
            do j = 2,3
c     for VW+ -> e+ ve A:
               call vwptowpa(q12(0,j),q34(0,j),v,h,zwtowa(0,0,j),awtowa(0,0,j)) !q12 = V
               call vwptowpa(q34(0,j),q12(0,j),v,h,wztowa(0,0,j),watowa(0,0,j)) !q12 = W+

C     for W+W- -->  A
               call calc_wwtoa(1,iu,q34(0,j),v,wwtoa(0,0,1,j)) !emit W+ -> e+ nu on upper
               call calc_wwtoa(1,iu,q12(0,j),v,wwtoa(0,0,2,j)) !emit W+ -> e+ nu on lower

C     for WV --> e+ ve
               call WPVtoln(1,q34(0,j),v,awen(0,0,1,j),zwen(0,0,1,j)) !emit A on upper
               call WPVtoln(2,q12(0,j),v,awen(0,0,2,j),zwen(0,0,2,j)) !emit A on lower
               
               call WPVtoln(1,q12(0,j),v,waen(0,0,2,j),wzen(0,0,2,j)) !emit A on lower
               call WPVtoln(2,q34(0,j),v,waen(0,0,1,j),wzen(0,0,1,j)) !emit A on upper
               
            enddo               !j
c     for W+ --> e+ ve A
            call wptowpa(qvv,v,h,wtowa)
            
         endif

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
            
#endif

         elseif (with_anom) then ! anomalous gauge boson couplings
c     using global form factor for all tensors of one phase space point
c     this ensures proper cancellations for anomalous contributions
c     energy scale is invariant WA mass        
            do j = 2,3
               call anomal_formfactor(q12(0,j),q34(0,j),qe(0),qa(0))

c     for VW+ -> e+ ve A:
               call vwmtowma_anomal(q12(0,j),q34(0,j),v,h,zwtowa(0,0,j),awtowa(0,0,j)) !q12 = V
               call vwmtowma_anomal(q34(0,j),q12(0,j),v,h,wztowa(0,0,j),watowa(0,0,j)) !q12 = W+
               
C     for W+W- -->  A
               call calc_wwtoa_anomal(2,iu,q34(0,j),v,wwtoa(0,0,1,j)) !emit W+ -> e+ nu on upper
               call calc_wwtoa_anomal(2,iu,q12(0,j),v,wwtoa(0,0,2,j)) !emit W+ -> e+ nu on lower

C     for WV --> e+ ve
               call WMVtoln_anomal(1,q34(0,j),v,awen(0,0,1,j),zwen(0,0,1,j)) !emit A on upper
               call WMVtoln_anomal(2,q12(0,j),v,awen(0,0,2,j),zwen(0,0,2,j)) !emit A on lower
c     
               call WMVtoln_anomal(1,q12(0,j),v,waen(0,0,2,j),wzen(0,0,2,j)) !emit A on lower
               call WMVtoln_anomal(2,q34(0,j),v,waen(0,0,1,j),wzen(0,0,1,j)) !emit A on upper
               
            enddo               !j
c     for W- --> e- ve A
            call wmtowma_anomal(qvv,v,h,wtowa)

         else                   !SM
            do j = 2,3
               
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
            enddo
            call wmtowma(qvv,v,h,wtowa)

         endif

      endif     

c -----------------------------------------------------------------------
c
c debugging:	
c
c subtraction for single W production:
c
c      ldebug = lok(1).and.lok(2).and.lok(3).and.p(0,5,1).gt.125
      if (bos.lt.10) then
      if (ldebug) then
         do j = 3,3
            if(j.eq.2) then
               qa3 = sqrt(2*dot0p(p(0,1,1),p(0,3,1)))
               qa5 = sqrt(2*dot0p(p(0,1,1),p(0,5,1)))
               q35 = sqrt(2*dot0p(p(0,3,1),p(0,5,1)))
               rjj = rsep(p(0,3,1),p(0,5,1))
            elseif(j.eq.3) then
               qa3 = sqrt(2*dot0p(p(0,2,1),p(0,4,1)))
               qa5 = sqrt(2*dot0p(p(0,2,1),p(0,5,1)))
               q35 = sqrt(2*dot0p(p(0,4,1),p(0,5,1)))
               rjj = rsep(p(0,4,1),p(0,5,1))
            endif
         enddo
      endif
c      ldebug = ldebug .and. rjj.lt.0.06 !q35/qa3.lt.0.01 !.and.q35/qa3.lt.0.1
c      ldebug = .false.
      if (ldebug) then
         print*,"----------------next event---------------------------"
         print 193," qa5/qa3 = ",qa5/qa3," q35/qa3 = ",q35/qa3,
     1             " rjj = ",rjj
 193     format(3(a,f10.6))
         print*," q12 = ",sqrt(q12(4,1))," q34 = ",sqrt(q34(4,1))
         q_sf = sqrt(mufsq(1,1))
         call pdfproton( xi(1), q_sf, pdf(-6,1,1) )        
         q_sf = sqrt(mufsq(2,1))
         call pdfproton( xi(2), q_sf, pdf(-6,2,1) )        
c         ldebug = ldebug .and. pdf(4,1,1).ne.0 .and. pdf(4,2,1).ne.0
      endif

      if ( ldebug ) then
 10      format( " p(", i1, ") = ", 4(f10.3, 2x) )
 20      format( " v(", i1, ") = ", 4(f10.3, 2x) )
c
         write(6,*) " "
         do i = 1, n_v
            write(6,20) i, v(0,i), v(1,i), v(2,i), v(3,i)
         end do
         qmina = 1d30
         do j = 1, 3
            if (j.eq.1) print*," (n+1) parton momenta "
            if (j.eq.2) print*," n-parton momenta for upper NLO corr."
            if (j.eq.3) print*," n-parton momenta for lower NLO corr."
            do i = 1, n_p-min(j-1,1)
               write(6,10) i, p(0,i,j), p(1,i,j), p(2,i,j), p(3,i,j)
            end do
            if(j.eq.2) then
               qa3 = sqrt(2*dot0p(p(0,1,1),p(0,3,1)))
               qa5 = sqrt(2*dot0p(p(0,1,1),p(0,5,1)))
               q35 = sqrt(2*dot0p(p(0,3,1),p(0,5,1)))
               print*," qa3 = ",qa3," qa5 = ",qa5," q35 = ",q35
               qmina = min(qmina,qa3,qa5)
            elseif(j.eq.3) then
               qa3 = sqrt(2*dot0p(p(0,2,1),p(0,4,1)))
               qa5 = sqrt(2*dot0p(p(0,2,1),p(0,5,1)))
               q35 = sqrt(2*dot0p(p(0,4,1),p(0,5,1)))
               print*," qa4 = ",qa3," qa5 = ",qa5," q45 = ",q35
               qmina = min(qmina,qa3,qa5)
            endif
         enddo
c
c debugging: scale info
c
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
      endif !bos<10
c
c  ----------- end debugging for single W/Z production -----------------
c  ---------------------------------------------------------------------
c
c  debugging for WZ:
c
c	if (bos.gt.10) then
c collinear subtraction:
c      nq = 4
c      ldebug_sub = lok(1).and.lok(2).and.lok(3).and.
c     #         dot0p(p(0,5,1),p(0,nq,1)).lt.1
c
c      if (ldebug_sub) print*,"coll.dot.:",lok(1),lok(2),lok(3),
c     #  			    dot0p(p(0,5,1),p(0,nq,1))
c      
c soft subtraction:
c
c      ldebug_sub = lok(1).and.lok(2).and.lok(3).and.p(0,5,1).lt.0.1
c      
c      if (ldebug_sub) print*,"soft.mom.:",lok(1),lok(2),lok(3),p(0,5,1)
c
c	endif !bos>10 
	
c  ---------------------------------------------------------------------

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

         call wbf_w3j(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1            pdf,res,nmin,nmax)

         return
       endif

c
c
c call PDF subroutine in order to determine parton
c distributions in the incoming (anti)protons.
c for the NLO contributions x1 = x*y in my notes with x=z=xuz(1,j), y=xi(j-1)
c
      x1 = xi(1)*xuz(1,2)
      q_sf = sqrt(mufsq(1,1))
      call pdfproton( xi(1), q_sf, pdf(-6,1,1) )      ! f_a(y)=f_a(x1/z)

      if (nlo.gt.0) then
      q_sf = sqrt(mufsq(1,2))                         ! f_a(x1) for upper line 
      call pdfproton( x1, q_sf, pdf(-6,1,2) )         !   NLO correction
      endif

      if (mufsq(1,3).ne.mufsq(1,1)) then              ! f_a(x1) for lower line
         q_sf = sqrt(mufsq(1,3))                      !   NLO correction
         call pdfproton( xi(1), q_sf, pdf(-6,1,3) )
      else
         do i = -6,6
            pdf(i,1,3) = pdf(i,1,1)
         enddo
      endif

      x2 = xi(2)*xuz(1,3)
      q_sf = sqrt(mufsq(2,1))
      call pdfproton( xi(2), q_sf, pdf(-6,2,1) )      ! f_b(y)=f_a(x2/z)

      if (nlo.gt.0) then
      q_sf = sqrt(mufsq(2,3))                         ! f_b(x2) for lower line 
      call pdfproton( x2, q_sf, pdf(-6,2,3) )         !   NLO correction
      endif

      if (mufsq(2,2).ne.mufsq(2,1)) then              ! f_b(x2) for upper line
         q_sf = sqrt(mufsq(2,2))                      !   NLO correction
         call pdfproton( xi(2), q_sf, pdf(-6,2,2) )
      else
         do i = -6,6
            pdf(i,2,2) = pdf(i,2,1)
         enddo
      endif
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
c      jsig = min(3*rn(1)+1d0,3.01d0)
      jsig = 0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     Now sum over the subprocesses contributing to W production

      nmax = FL_Wg(fsign,-1)            !reset counter for subprocesses to 0

C*******************  q1 q3 ---> q2 q4 g W (Z)   **********************

c   physToDiag(ext.momentum label) = Feynman diagram label

      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=3
      physToDiag(3)=2    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=4
      physToDiag(5)=5   ! gluon
C NOTE: for call of wbf_w3j it is important that p(*,1,*) and p(*,3,*)
c correspond to 1-2 fermion line ALWAYS, i.e physToDiag(1/2)={1,3} and 
c similarly physToDiag(3/4)={2,4} for the 3-4 fermion line
      fsign(1) = 1
      fsign(2) = 1
      fsign(3) = 1
      fsign(4) = 1
      gsign    = 1
      
      call wbf_w3j(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1             pdf,res,nmin,nmax)
c      if (ldebug) call resprintw(nmin,nmax,res)
      if (init.eq.1) nproc(1) = nmax
            
C******************* q1 qb4 ---> q2 qb3 g W (Z)   **********************
      
      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=4
      physToDiag(3)=2    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=3
      physToDiag(5)=5   ! gluon
c
      fsign(1) =  1
      fsign(2) =  1
      fsign(3) = -1
      fsign(4) = -1
      gsign    =  1
      
      call wbf_w3j(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1  	 pdf,res,nmin,nmax)
c      if (ldebug) call resprintw(nmin,nmax,res)
      if (init.eq.1) nproc(2) = nmax
      
C******************* qbar2 q3 ---> qbar1 q4 g W (Z)   **********************
      
      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=3
      physToDiag(3)=1    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=4
      physToDiag(5)=5   ! gluon
c
      fsign(1) = -1
      fsign(2) = -1
      fsign(3) =  1
      fsign(4) =  1
      gsign    =  1
      
      call wbf_w3j(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1  	 pdf,res,nmin,nmax)
c      if (ldebug) call resprintw(nmin,nmax,res)
      if (init.eq.1) nproc(3) = nmax

C*******************  qbar2 qb4 ---> qbar1 qb3 g W (Z)   **********************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=4
      physToDiag(3)=1    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=3
      physToDiag(5)=5   ! gluon
c
      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = -1
      fsign(4) = -1
      gsign    =  1
      
      call wbf_w3j(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1  	 pdf,res,nmin,nmax)
c      if (ldebug) call resprintw(nmin,nmax,res)
      if (init.eq.1) nproc(4) = nmax
         
c*******************  g q3 ---> qb1 q2 q4 W (Z)   **********************

c 899  continue
      physToDiag(1)=5             ! note: phystodiag(5) must be fsign=1
      physToDiag(2)=3             ! for convention in wbf_w3j
      physToDiag(3)=1             ! for ALL GLUON INITIATED PROCESSES
      physToDiag(4)=4
      physToDiag(5)=2
c
      fsign(1) = -1
      fsign(2) =  1
      fsign(3) =  1
      fsign(4) =  1
      gsign    = -1
      
      call wbf_w3j(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1  	 pdf,res,nmin,nmax)
c      if (ldebug) call resprintw(nmin,nmax,res)
      if (init.eq.1) nproc(5) = nmax
      
C*******************  q1 g ---> q2 qb3 q4 W (Z)   **********************
      
      physToDiag(1)=1
      physToDiag(2)=5
      physToDiag(3)=2
      physToDiag(4)=3
      physToDiag(5)=4
c
      fsign(1) =  1
      fsign(2) =  1
      fsign(3) = -1
      fsign(4) =  1
      gsign    = -1
      
      call wbf_w3j(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1  	 pdf,res,nmin,nmax)
c      if (ldebug) call resprintw(nmin,nmax,res)
      if (init.eq.1) nproc(6) = nmax
         
C*******************  g qbar4 ---> qb1 q2 qbar3 W (Z)  **********************

      physToDiag(1)=5
      physToDiag(2)=4
      physToDiag(3)=1
      physToDiag(4)=3
      physToDiag(5)=2
c
      fsign(1) = -1
      fsign(2) =  1
      fsign(3) = -1
      fsign(4) = -1
      gsign    = -1
         
      call wbf_w3j(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1  	 pdf,res,nmin,nmax)
c      if (ldebug) call resprintw(nmin,nmax,res)
      if (init.eq.1) nproc(7) = nmax
         
C*******************  qbar2 g ---> qbar1 qb3 q4 W (Z)   **********************

      physToDiag(1)=2
      physToDiag(2)=5
      physToDiag(3)=1
      physToDiag(4)=3
      physToDiag(5)=4
c
      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = -1
      fsign(4) =  1
      gsign    = -1
      
      call wbf_w3j(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1  	 pdf,res,nmin,nmax)
c      if (ldebug) call resprintw(nmin,nmax,res)
      if (init.eq.1) nproc(8) = nmax
         

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

        if (procID.eq.WPZhadjj) then


         CALL JIOXXX(lp,lm,GZ_ZF ,ZMASS,ZWIDTH,zu)      !Zl
         CALL JIOXXX(lp,lm,GZ_AF ,ZERO ,ZERO  ,au)      !Al
	 
         do mu = 0,3
            qe(mu) = v(mu,1)+v(mu,2)
            qu(mu) = v(mu,3)+v(mu,4)
            qvv(mu) = qe(mu) + qu(mu)
         enddo
	 	 
         qe(4) = qe(0)**2-qe(1)**2-qe(2)**2-qe(3)**2
         qu(4) = qu(0)**2-qu(1)**2-qu(2)**2-qu(3)**2
 
         do j = 2,3
         call anomal_formfactor(q12(0,j),q34(0,j),qe(0),qu(0))
c for ZW+,AW+ -> e+ ve mu+ mu-:
         call vwtovwp_had(q12(0,j),q34(0,j),v,h,zwwv(0,0,j),awwv(0,0,j)) !q12 = V
         call vwtovwp_had(q34(0,j),q12(0,j),v,h,wzwv(0,0,j),wawv(0,0,j)) !q12 = W+
C for V1V2 -->  mu+ mu- (Vi=A,Z)
         call VVtomm_had(1,h,q12(0,j),v,aauu(0,0,2,j),azuu(0,0,2,j),
     #		 zauu(0,0,2,j),zzuu(0,0,2,j))    	 !emit W on lower
         call VVtomm_had(2,h,q34(0,j),v,aauu(0,0,1,j),azuu(0,0,1,j),
     #		 zauu(0,0,1,j),zzuu(0,0,1,j))    	 !emit W on upper
C for WW -->  mu+ mu-
         call WWtomm_had(1,h,q34(0,j),v,wwuu(0,0,1,j))     !emit V->e nu on upper
         call WWtomm_had(1,h,q12(0,j),v,wwuu(0,0,2,j))     !emit V->e nu on lower
C for WV --> e+ ve
         call WPVtoln_had(1,q34(0,j),v,awen(0,0,1,j),zwen(0,0,1,j))  !emit V->mu+mu- on upper
         call WPVtoln_had(2,q12(0,j),v,awen(0,0,2,j),zwen(0,0,2,j))  !emit V->mu+mu- on lower
         call WPVtoln_had(1,q12(0,j),v,waen(0,0,2,j),wzen(0,0,2,j))  !emit V->mu+mu- on lower
         call WPVtoln_had(2,q34(0,j),v,waen(0,0,1,j),wzen(0,0,1,j))  !emit V->mu+mu- on upper
  	   enddo !j
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
	        
        do j = 2,3
         call anomal_formfactor(q12(0,j),q34(0,j),qe(0),qu(0))
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
	  enddo !j		 
c for W- --> e- ve mu+ mu-
	  call wmto4l_had(v,h,wwv)	


        endif


C*******************  q1 q3 ---> q2 q4 g W (Z)   **********************

c   physToDiag(ext.momentum label) = Feynman diagram label

      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=3
      physToDiag(3)=2    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=4
      physToDiag(5)=5   ! gluon
C NOTE: for call of wbf_w3j it is important that p(*,1,*) and p(*,3,*)
c correspond to 1-2 fermion line ALWAYS, i.e physToDiag(1/2)={1,3} and 
c similarly physToDiag(3/4)={2,4} for the 3-4 fermion line
      fsign(1) = 1
      fsign(2) = 1
      fsign(3) = 1
      fsign(4) = 1
      gsign    = 1
      
      call wbf_w3j(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1             pdf,res,nmin,nmax)
c      if (ldebug) call resprintw(nmin,nmax,res)
      if (init.eq.1) nproc(1) = nmax
            
C******************* q1 qb4 ---> q2 qb3 g W (Z)   **********************
      
      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=4
      physToDiag(3)=2    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=3
      physToDiag(5)=5   ! gluon
c
      fsign(1) =  1
      fsign(2) =  1
      fsign(3) = -1
      fsign(4) = -1
      gsign    =  1
      
      call wbf_w3j(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1  	 pdf,res,nmin,nmax)
c      if (ldebug) call resprintw(nmin,nmax,res)
      if (init.eq.1) nproc(2) = nmax
      
C******************* qbar2 q3 ---> qbar1 q4 g W (Z)   **********************
      
      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=3
      physToDiag(3)=1    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=4
      physToDiag(5)=5   ! gluon
c
      fsign(1) = -1
      fsign(2) = -1
      fsign(3) =  1
      fsign(4) =  1
      gsign    =  1
      
      call wbf_w3j(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1  	 pdf,res,nmin,nmax)
c      if (ldebug) call resprintw(nmin,nmax,res)
      if (init.eq.1) nproc(3) = nmax

C*******************  qbar2 qb4 ---> qbar1 qb3 g W (Z)   **********************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=4
      physToDiag(3)=1    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=3
      physToDiag(5)=5   ! gluon
c
      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = -1
      fsign(4) = -1
      gsign    =  1
      
      call wbf_w3j(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1  	 pdf,res,nmin,nmax)
c      if (ldebug) call resprintw(nmin,nmax,res)
      if (init.eq.1) nproc(4) = nmax
         
c*******************  g q3 ---> qb1 q2 q4 W (Z)   **********************

c 899  continue
      physToDiag(1)=5             ! note: phystodiag(5) must be fsign=1
      physToDiag(2)=3             ! for convention in wbf_w3j
      physToDiag(3)=1             ! for ALL GLUON INITIATED PROCESSES
      physToDiag(4)=4
      physToDiag(5)=2
c
      fsign(1) = -1
      fsign(2) =  1
      fsign(3) =  1
      fsign(4) =  1
      gsign    = -1
      
      call wbf_w3j(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1  	 pdf,res,nmin,nmax)
c      if (ldebug) call resprintw(nmin,nmax,res)
      if (init.eq.1) nproc(5) = nmax
      
C*******************  q1 g ---> q2 qb3 q4 W (Z)   **********************
      
      physToDiag(1)=1
      physToDiag(2)=5
      physToDiag(3)=2
      physToDiag(4)=3
      physToDiag(5)=4
c
      fsign(1) =  1
      fsign(2) =  1
      fsign(3) = -1
      fsign(4) =  1
      gsign    = -1
      
      call wbf_w3j(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1  	 pdf,res,nmin,nmax)
c      if (ldebug) call resprintw(nmin,nmax,res)
      if (init.eq.1) nproc(6) = nmax
         
C*******************  g qbar4 ---> qb1 q2 qbar3 W (Z)  **********************

      physToDiag(1)=5
      physToDiag(2)=4
      physToDiag(3)=1
      physToDiag(4)=3
      physToDiag(5)=2
c
      fsign(1) = -1
      fsign(2) =  1
      fsign(3) = -1
      fsign(4) = -1
      gsign    = -1
         
      call wbf_w3j(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1  	 pdf,res,nmin,nmax)
c      if (ldebug) call resprintw(nmin,nmax,res)
      if (init.eq.1) nproc(7) = nmax
         
C*******************  qbar2 g ---> qbar1 qb3 q4 W (Z)   **********************

      physToDiag(1)=2
      physToDiag(2)=5
      physToDiag(3)=1
      physToDiag(4)=3
      physToDiag(5)=4
c
      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = -1
      fsign(4) =  1
      gsign    = -1
      
      call wbf_w3j(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1  	 pdf,res,nmin,nmax)
c      if (ldebug) call resprintw(nmin,nmax,res)
      if (init.eq.1) nproc(8) = nmax

c***  end of extra process evaluation for hadronic decay  ***

       endif
      end select

      SELECT CASE(procid)
      CASE(WPhadZjj, WMhadZjj)
       call hadDecayFactor_W(finalquarks(1), abs(mjj2(v(0,1),v(0,2))), N_gen_W, fac_W)
       do j=0,3
        DO IRES = 1,NMAX
          res(IRES,j) = res(IRES,j) * fac_W            ! factor fac_W for hadronic decay
        ENDDO
       enddo
      CASE(WPZhadjj, WMZhadjj)
       call hadDecayFactor_Z(finalquarks(1), abs(mjj2(v(0,3),v(0,4))), N_gen_up, N_gen_down, fac_Z_up, fac_Z_down)
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

      init = init+1
      do j=0,3
         m2s(j) = 0
         DO IRES = 1,NMAX
            m2s(j) = m2s(j) + RES(IRES,j)
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
         endif


      	if (bos.eq.3.or.bos.eq.4) then
       	    if (jsig.ne.0) m2s(j) = m2s(j)*3d0  ! Insert weight factor 3 
	    			 		! for random helicity summation
      	elseif (bos.eq.32.or.bos.eq.42.or.bos.eq.31.or.bos.eq.41) then
      	    m2s(j) = m2s(j)*2d0  ! factor 2 for random helicity summation
      	endif
	 
      enddo


      if (lwarn .and. .false.) then
         do j=1,3
            if (lok(j) .and. m2s(j).eq.0 ) then
               print*,j," lok(j) = ",lok(j)," m2s = ",m2s(j)
            endif
         enddo
         if (m2s(0).ne.0) then
            if (abs((m2s(1)+m2s(2)+m2s(3))/m2s(0)-1).gt.1d-5) then
               print*," m2s(0) = ",m2s(0),(m2s(1)+m2s(2)+m2s(3))/m2s(0)
            endif
         else
            print*," m2s = ",m2s
         endif
         if (nmax.ne.nmaxold) print*," nmax: ",nmaxold,nmax
         nmaxold = nmax
      endif
      if (ldebug) then
         print*," lokt = ",lok
         if (qmina.lt.2d0) read(*,*)
      endif
c      if (ldebug) read(*,113) text
c 113  format(a)
c ---------------------------------------------------------------
c
c debugging (check subtraction):
c
      if (ldebug_sub) then
         do j=1,3
            if (lok(j) .and. m2s(j).eq.0 ) then
               print*,j," lok(j) = ",lok(j)," m2s = ",m2s(j)
            endif
         enddo
         if (m2s(0).ne.0) then
c            if (abs((m2s(1)+m2s(2)+m2s(3))/m2s(0)-1).gt.1d-5) then
               print*," m2s(0) = ",m2s(0),(m2s(1)+m2s(2)+m2s(3))/m2s(0)
               print*," m2s(0,1,2+3) = ",m2s(0),m2s(1),(m2s(2)+m2s(3))
               print*," ratio subtr/real = ",(m2s(2)+m2s(3))/m2s(1)
               print*," m2s(0,1,2,3) = ",m2s
	       if ( abs(1-abs((m2s(2)+m2s(3))/m2s(1))).gt.1e-2) read(*,*)
	       print*
c            endif
         else
            print*," m2s = ",m2s
         endif
      endif
c      if (ldebug) print*," lokt = ",lok
c      if (ldebug) read(*,113) text
c 113  format(a)

      return
      end


c******************************************************************
c
c   begin subroutine wbf_w3j
c
c*****************************************************************
      subroutine wbf_w3j(xi,p,v,physToDiag,fsign,gsign,bos, 
     1                    nlo,lok,xuz,
     1                    pdf,res,nmin,nmax)
          use globalvars, only: ldoblha
      implicit none
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
      real*8 p(0:3,max_p,3), v(0:3,max_v), xi(nx), xuz(2,2:3)
      real*8 pdf(-6:6,2,3), res(maxnumsubproc,0:3)
      integer physToDiag(5), fsign(4+max_v), gsign, 
     1	      bos, nlo, nmin, nmax,h
      logical lok(3)
c      
      common /hcount / h
c
c wbf_w3j calls the amplitude square routines 
c             qqwpqqj    for qq-->qqW+ g   and
c             qqwpqq     for qq-->qqW+ 
c 	or 
c             qqwpzqqj  for qq-->qqW+Z g   and
c             qqwpzqq   for qq-->qqW+Z 
c 	or 
c             qqwpaqqj  for qq-->qqW+A g   and
c             qqwpaqq   for qq-->qqW+A 
c
c     (and similarly for W- (Z/A) final states)
c
c     for the subtraction terms for the NLO cross section calculation
c
c  INPUT:  p(0:3,5,3)      external physical parton momenta
c          v(0:3,nv)        W decay momenta
c          physToDiag(5)   physToDiag(ext.mom. label) = Feynman diagram label
c          fsign,gsign     sign factors for amplitude calls; see qqwpqq(j)
c          nlo             nlo = 0: calculate LO only (i.e. no subtraction)
c                          nlo = 1: calculate full NLO subtraction
c          lok(3)          lok(ID)=T means momenta set ID passes acceptance
c                          cuts, i.e res(k,ID) needs to be calculated
c  OUTPUT:
c          res(k,ID)       output corresponding to res(*,ID) in m2s_qqWjqq
c                          ID = 1  : the real emission |M|^2 * pdf
c                             = 2,3: sutraction terms for emision off 12 or 
c                                    34 line
c                             = 0  : subtracted result which drives 
c                                    integration, i.e 
c                                    res(*,0) = res(*,1)+res(*,2)+res(*,3)
c
c  In and Output
c    nmin, nmax            range of process numbers for this call
c
      real*8 C2, TR, N                  ! color factors
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
c      integer jsig, jsig1, jsig3, jsig5
c      common /chelsum/ jsig,jsig1,jsig3,jsig5
c
      real*8 uucsb(-1:1,-1:1,2:3), ddcsb(-1:1,-1:1,2:3),
     1       udccb(-1:1,-1:1,2:3), udssb(-1:1,-1:1,2:3),
     1       uucs(0:3),ddcs(0:3),udcc(0:3),udss(0:3)
      real*8 uuscb(-1:1,-1:1,2:3), ddscb(-1:1,-1:1,2:3),
     1       duccb(-1:1,-1:1,2:3), dussb(-1:1,-1:1,2:3),
     1       uusc(0:3),ddsc(0:3),ducc(0:3),duss(0:3)
      real*8 pbar(0:3,4+max_v), qbar(0:4), q2,v2,q3,v3,sub(2:3), 
     1       e_in(2), dot0p, gs2(2:3), polcolq, polcolg
      real*8 xa
      double precision matr3sq(0:1,0:1,2,2:3)
      double precision matr2sq(0:1,0:1,2,2:3)
      integer iflav(5), FL_Wg
      external dot0p, FL_Wg

      logical loldp, ldebug
      parameter (ldebug=.false.)
      data e_in/0d0,0d0/
      save e_in,loldp, uucsb, ddcsb, udccb, udssb,
     1     uuscb,ddscb,duccb,dussb,
     1     gs2,polcolq,polcolg
      integer i,i1,i2,if1,if2,ifmin,j,k,lw,mu,qw,bosi, 
     1        fc(6), fcb(6)

      do j = 0,3
         uucs(j) = 0d0
         ddcs(j) = 0d0
         udcc(j) = 0d0
         udss(j) = 0d0
c
         uusc(j) = 0d0
         ddsc(j) = 0d0
         ducc(j) = 0d0
         duss(j) = 0d0
      enddo
      do if1 = 0,1
         do if2 = 0,1
            do i = 1,2
               do j = 2,3
                  matr3sq(if1,if2,i,j) = 0d0
                  matr2sq(if1,if2,i,j) = 0d0
               enddo
            enddo
         enddo
      enddo
      nmin = nmax+1
      if (bos.lt.10) then
         qw = 7-2*bos              ! W charge, +1 for W+, -1 for W-
      	 bosi = bos
      elseif(bos.eq.32.or.bos.eq.42) then
         qw   =  37/5-bos/5
	 bosi = -1/5+bos/10
      elseif(bos.eq.31.or.bos.eq.41) then
         qw   =  36/5-bos/5
	 bosi = -1/10+bos/10
      endif 

      loldp = e_in(1).eq.p(0,3,1) .and. e_in(2).eq.p(0,4,1)
c  reset the LO amplitude to 0 to avoid wrong subtraction for initial gluon
      if (.not.loldp) then            ! this is a new phase space point
         e_in(1) = p(0,3,1)
         e_in(2) = p(0,4,1)
         do i = -1,1
            do j = -1,1
               do k = 2,3
                  uucsb(i,j,k) = 0d0
                  ddcsb(i,j,k) = 0d0
                  udccb(i,j,k) = 0d0
                  udssb(i,j,k) = 0d0
c
                  uuscb(i,j,k) = 0d0
                  ddscb(i,j,k) = 0d0
                  duccb(i,j,k) = 0d0
                  dussb(i,j,k) = 0d0
               enddo
            enddo
         enddo
c  determine strong coupling gs for the two quark lines
         gs2(2) = 4d0*pi*als(1,1)
         gs2(3) = 4d0*pi*als(2,1)
         polcolq = 1d0/(4d0*N**2*xi(1)*xi(2))
         polcolg = 1d0/(4d0*N*(N**2-1)*xi(1)*xi(2))
      endif
         
c
c get the real emission amplitude squared, store it in uucs(1) etc.
      do mu = 0,3
         do i = 1,5
            pbar(mu,physToDiag(i))=p(mu,i,1)
         enddo
         qbar(mu) = pbar(mu,5)
      enddo	 
      qbar(4) = 0d0
c

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
         print*," wbf_w3j info: "
         print*," enter with bos = ",bos," fsign = ",fsign
      endif     
c
      if (bos.eq.3) then                   !the W+ case
         call qqwpqqj_c(pbar,fsign,qbar,gsign, 
     1                  uucs(1),ddcs(1),udcc(1),udss(1))
         if (ldebug) then
            print*," uucs = ",uucs
            print*," ddcs = ",ddcs
            print*," udcc = ",udcc
            print*," udss = ",udss
         endif
c                       
      elseif (bos.eq.4) then               !the W- case
         call qqwmqqj_c(pbar,fsign,qbar,gsign, 
     1                  uusc(1),ddsc(1),ducc(1),duss(1))
         if (ldebug) then
            print*," uusc = ",uusc
            print*," ddsc = ",ddsc
            print*," ducc = ",ducc
            print*," duss = ",duss
         endif
	 
      endif !bos = 4
c
      elseif(bos.eq.32) then         ! W+Z 
c
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
         
	 if (lok(1).or.nlo.eq.0) then	 
          call qqwpzqqj(pbar,fsign,qbar,gsign,bos,
     1                    uucs(1),ddcs(1),udcc(1),udss(1)) 	
	 endif
          						
           
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
	 if (lok(1).or.nlo.eq.0) then
           call qqwmzqqj(pbar,fsign,qbar,gsign,bos,
     1                    uusc(1),ddsc(1),ducc(1),duss(1)) 
 
         endif   
	 
      elseif(bos.eq.41) then         ! W-A 
         do mu = 0,3             
            pbar(mu,5) = v(mu,1) ! l-
            pbar(mu,6) = v(mu,2) ! nubar
            pbar(mu,7) = v(mu,3) ! a
         enddo
         fsign(5) = 1
         fsign(6) = -1
         fsign(7) = 1
	 if (lok(1).or.nlo.eq.0) then
           call qqwmaqqj(pbar,fsign,qbar,gsign,bos,
     1                    uusc(1),ddsc(1),ducc(1),duss(1))  
         endif 

      elseif(bos.eq.31) then         ! W+A 

         do mu = 0,3             
            pbar(mu,5) = v(mu,1) ! l+
            pbar(mu,6) = v(mu,2) ! nu
            pbar(mu,7) = v(mu,3) ! a
         enddo
         fsign(5) = -1
         fsign(6) = 1
         fsign(7) = 1


	 if (lok(1).or.nlo.eq.0) then	 
          call qqwpaqqj(pbar,fsign,qbar,gsign,bos,
     1                    uucs(1),ddcs(1),udcc(1),udss(1)) 
          endif

      else
         print*," undefined variable bos in wbf_w3j "
         stop
      endif !bos

      if (ldoblha) then
        return
      endif
      
      if (bos.eq.3.or.bos.eq.32.or.bos.eq.31) then
         do j = 2,3
            matr3sq(0,(1-fsign(4))/2, 2,j) = uucs(j) ! (1-fs(4))/2 = flavor id
            matr3sq(1,(1-fsign(4))/2, 2,j) = ddcs(j) ! of beam 2 (= 0,1)
            matr3sq((1-fsign(2))/2,0, 1,j) = udcc(j) ! 1,j id"s W attached
            matr3sq((1-fsign(2))/2,1, 1,j) = udss(j) ! to 12 line, g to j line
         enddo
c matr3sq(i1,i2,lw,j): index i1 = flavor of beam 1 for gluon not beam 1
c                            i2 = flavor of beam 2 for gluon not beam 2
c                         i1,i2 = flavor of first,third quark otherwise
c                      index lw = line to which W is attached (1=12, 2=34 line)
c                      index j  = line to which gluon is attached (2=12,3=34)
      elseif (bos.eq.4.or.bos.eq.42.or.bos.eq.41) then
         do j = 2,3
            matr3sq(0,(1+fsign(4))/2, 2,j) = uusc(j) ! (1-fs(3))/2 = flavor id
            matr3sq(1,(1+fsign(4))/2, 2,j) = ddsc(j) ! of beam 2 (= 0,1)
            matr3sq((1+fsign(2))/2,0, 1,j) = ducc(j) ! 1,j id"s W attached
            matr3sq((1+fsign(2))/2,1, 1,j) = duss(j) ! to 12 line, g to j line
         enddo
      endif !bos	 
   
      
      if (ldebug) then
         print*," real emission M^2 = "
         print 104," matr3sq = ",((matr3sq(0,i,2,j),i=0,1),j=2,3)
         print 104,"           ",((matr3sq(1,i,2,j),i=0,1),j=2,3)
         print 104,"           ",((matr3sq(i,0,1,j),i=0,1),j=2,3)
         print 104,"           ",((matr3sq(i,1,1,j),i=0,1),j=2,3)
 104     format(a,4g12.4)
 
c        read (*,113) text
c 113     format(a)
      endif
         

c for the NLO case get the subtraction terms; 
c first the case with a final state gluon
      if (gsign.eq.1) then
         if (ldebug) then
            print*," final state gluon section in wbf_W3j, lok = ",lok
c            print*," jsig = ",jsig," jsig1,3 = ",jsig1,jsig3
            print 101," fsign = ",fsign
            print*," x1,x2 = ",xi
 101        format(a,6i5,a,i5)
         endif
         if (nlo.eq.1) then
            do j = 2,3          ! j=2: emission off 1,2 line
               do mu = 0,3      ! j=3: emission off 3,4 line
                  do i = 1,4
                     pbar(mu,physToDiag(i))=p(mu,i,j)
                  enddo
               enddo
               if (bos.eq.3.and.lok(j)) then
                  call qqwpqq(pbar,fsign,0,
     1                 uucs(j),ddcs(j),udcc(j),udss(j))
               elseif (bos.eq.4.and.lok(j)) then
                  call qqwmqq(pbar,fsign,0,
     1                 uusc(j),ddsc(j),ducc(j),duss(j))
               elseif (bos.eq.32.and.lok(j)) then
                  call qqwpzqq(pbar,fsign,0,j,bos,
     1                 uucs(j),ddcs(j),udcc(j),udss(j))
               elseif (bos.eq.42.and.lok(j)) then
                  call qqwmzqq(pbar,fsign,0,j,bos,
     1                 uusc(j),ddsc(j),ducc(j),duss(j))
               elseif (bos.eq.31.and.lok(j)) then
                  call qqwpaqq(pbar,fsign,0,j,bos,1, !ps_number only relevant for virtuals
     1                 uucs(j),ddcs(j),udcc(j),udss(j))
               elseif (bos.eq.41.and.lok(j)) then
                  call qqwmaqq(pbar,fsign,0,j,bos,1,
     1                 uusc(j),ddsc(j),ducc(j),duss(j))
               endif
	       
	       
	       if (bos.eq.3.or.bos.eq.32.or.bos.eq.31) then
                  matr2sq(0,(1-fsign(4))/2, 2,j) = uucs(j)
                  matr2sq(1,(1-fsign(4))/2, 2,j) = ddcs(j)
                  matr2sq((1-fsign(2))/2,0, 1,j) = udcc(j)
                  matr2sq((1-fsign(2))/2,1, 1,j) = udss(j)
c save matrix elements for later use with initial gluons
                  uucsb(fsign(1),fsign(3),j) = uucs(j)
                  ddcsb(fsign(1),fsign(3),j) = ddcs(j)
                  udccb(fsign(1),fsign(3),j) = udcc(j)
                  udssb(fsign(1),fsign(3),j) = udss(j)
               elseif (bos.eq.4.or.bos.eq.42.or.bos.eq.41) then
                  matr2sq(0,(1+fsign(4))/2, 2,j) = uusc(j)
                  matr2sq(1,(1+fsign(4))/2, 2,j) = ddsc(j)
                  matr2sq((1+fsign(2))/2,0, 1,j) = ducc(j)
                  matr2sq((1+fsign(2))/2,1, 1,j) = duss(j)
c save matrix elements for later use with initial gluons
                  uuscb(fsign(1),fsign(3),j) = uusc(j)
                  ddscb(fsign(1),fsign(3),j) = ddsc(j)
                  duccb(fsign(1),fsign(3),j) = ducc(j)
                  dussb(fsign(1),fsign(3),j) = duss(j)
               endif          ! W+ (Z)  or  W-(Z)
c
               q2 = 2*dot0p(qbar,p(0,j-1,1))*xuz(1,j)!p(mu,j-1,1) is incom.quark
               v2 = 2d0*gs2(j)*                    !i.e  j-1 = 1,2
     $              ( 2d0/(1-xuz(1,j)+xuz(2,j))-(1+xuz(1,j)))
               q3 = 2*dot0p(qbar,p(0,j+1,1))*xuz(1,j)!p(mu,j+1,1) is outgo.quark
               v3 = 2d0*gs2(j)*                    !i.e. j+1 = 3,4
     $              ( 2d0/(1-xuz(1,j)+xuz(2,j))-(2-xuz(2,j)))
               sub(j) = v2/q2+v3/q3
            enddo     ! loop j
         endif        ! nlo=1
         iflav(5) = 0           ! final state gluon id
         do lw = 1,2            ! identify line W is attached to
            do if1 = 1,(nfl/2)*2  ! identify flavor of NC vertex
               iflav(3-lw)=if1*fsign(physToDiag(3-lw))
               iflav(5-lw)=if1*fsign(physToDiag(5-lw))
               do if2=bosi/2+(1+fsign(2*lw))/2*qw,(nfl/2)*2,2
                  iflav(lw)=if2*fsign(physToDiag(lw))
                  iflav(lw+2)=if2*fsign(physToDiag(lw+2))-qw
                  i1 = mod(abs(iflav(1)),2)
                  i2 = mod(abs(iflav(2)),2)
                  do j = 2,3
                     k=FL_Wg(iflav,j)
                     if (lok(1)) then
                        res(k,1)=pdf(sign1*iflav(1),1,2*j-3)!1for j=2;3for j=3
     &                          *pdf(sign2*iflav(2),2,4-j  )!2for j=2;1for j=3
     &                          *matr3sq(i1,i2,lw,j)*polcolq
                     else
                        res(k,1) = 0
                     endif
                     if (nlo.eq.1 .and. lok(j)) then
                        res(k,j)=(pdf(sign1*iflav(1),1,j)
     &                           *pdf(sign2*iflav(2),2,j)
     &                           -pdf(sign1*iflav(1),1,2*j-3)
     &                           *pdf(sign2*iflav(2),2,4-j  )*sub(j))*C2
     &                           *matr2sq(i1,i2,lw,j)*polcolq
                        res(k,5-j)=0
                        res(k,0) = res(k,1)+res(k,j)
                     else
                        res(k,0) = res(k,1)
                        res(k,2) = 0
                        res(k,3) = 0
                     endif
                  enddo    ! j = 2,3
                  if (ldebug) then                                       !debug
                     print*," iflav = ",iflav(1),iflav(3),               !debug
     &                                  iflav(2),iflav(4)                !debug
                     print*," called matrix element ",i1,i2,lw," j"      !debug
                     if (lok(1) .and. nlo.eq.1) then                     !debug
                        xa = res(k-1,1)+res(k,1)                         !debug
                        print*,k,xa," up,low,tot",                       !debug
     &                         res(k-1,2)/xa,res(k,3)/xa,                !debug
     &                         (res(k-1,2)+res(k,3))/xa                  !debug
                     endif                                               !debug
                  endif                                                  !debug
               enddo          ! if2
            enddo             ! if1
         enddo                ! lw
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
               q2 = 2*dot0p(qbar,pbar(0,2))*xuz(1,j)
               fcb(1)=-fsign(1)
               q3 = 2*dot0p(qbar,pbar(0,1))*xuz(1,j)
            elseif(fsign(3).eq.-fsign(4)) then 
               j=3              ! j=3: emission off 3,4 line
               fc(4)=-fsign(4)
               q2 = 2*dot0p(qbar,pbar(0,4))*xuz(1,j)
               fcb(3)=-fsign(3)
               q3 = 2*dot0p(qbar,pbar(0,3))*xuz(1,j)
            endif
            v2 = 2d0*gs2(j) * ( (1-xuz(1,j))**2 + xuz(1,j)**2 )
            v3 = v2

            if (bos.eq.3.or.bos.eq.32.or.bos.eq.31) then
               matr2sq(0,(1-fsign(4))/2, 2,2) = uucsb(fc(1),fc(3),j)
               matr2sq(1,(1-fsign(4))/2, 2,2) = ddcsb(fc(1),fc(3),j)
               matr2sq((1-fsign(2))/2,0, 1,2) = udccb(fc(1),fc(3),j)
               matr2sq((1-fsign(2))/2,1, 1,2) = udssb(fc(1),fc(3),j)

               matr2sq(0,(1-fsign(4))/2, 2,3) = uucsb(fcb(1),fcb(3),j)
               matr2sq(1,(1-fsign(4))/2, 2,3) = ddcsb(fcb(1),fcb(3),j)
               matr2sq((1-fsign(2))/2,0, 1,3) = udccb(fcb(1),fcb(3),j)
               matr2sq((1-fsign(2))/2,1, 1,3) = udssb(fcb(1),fcb(3),j)
            elseif (bos.eq.4.or.bos.eq.42.or.bos.eq.41) then
               matr2sq(0,(1+fsign(4))/2, 2,2) = uuscb(fc(1),fc(3),j)
               matr2sq(1,(1+fsign(4))/2, 2,2) = ddscb(fc(1),fc(3),j)
               matr2sq((1+fsign(2))/2,0, 1,2) = duccb(fc(1),fc(3),j)
               matr2sq((1+fsign(2))/2,1, 1,2) = dussb(fc(1),fc(3),j)

               matr2sq(0,(1+fsign(4))/2, 2,3) = uuscb(fcb(1),fcb(3),j)
               matr2sq(1,(1+fsign(4))/2, 2,3) = ddscb(fcb(1),fcb(3),j)
               matr2sq((1+fsign(2))/2,0, 1,3) = duccb(fcb(1),fcb(3),j)
               matr2sq((1+fsign(2))/2,1, 1,3) = dussb(fcb(1),fcb(3),j)
            endif               ! W+  or  W-
         endif                  ! nlo = 0 or 1

         if (lwarn .and. fsign(phystodiag(5)).ne.1) then
            print*," WARNING: p5 is not quark" 
            print*," fsign(p5) = ",fsign(phystodiag(5))
            print*," fsign = ",fsign
            print*," phystodiag = ",phystodiag
         endif

         iflav(j-1) = 0
         do lw = 1,2
            if (j.eq.lw+1) then       ! gluon attached to CC line
               ifmin = 5-bosi          ! 2 for W+, 1 for W-
            else
               ifmin = 1+(1+qw*fsign(phystodiag(lw)))/2
            endif
            if (ldebug)print*," j = ",j," lw = ",lw," ifmin = ",ifmin
            do if1 = 1,(nfl/2)*2
               do if2=ifmin,(nfl/2)*2,2
                  if (j.eq.lw+1) then       ! gluon attached to CC line
                     iflav(3-lw)=if1*fsign(physToDiag(3-lw))
                     iflav(5-lw)=if1*fsign(physToDiag(5-lw))
                     iflav(j+1)=if2*fsign(physToDiag(j+1))
                     iflav(5)  =if2-qw
                  else                      ! gluon attached to NC line
                     iflav(j+1)=if1*fsign(physToDiag(j+1))
                     iflav(5)=if1*fsign(physToDiag(5))
                     iflav(lw)=if2*fsign(physToDiag(lw))
                     iflav(lw+2)=if2*fsign(physToDiag(lw+2))-qw
                  endif
                  i1 = mod(abs(iflav(7-2*j)),2)  ! 7-2j= 3,1 for j=2,3
                  i2 = mod(abs(iflav(2*j-2)),2)  ! 2j-2= 2,4 for j=2,3
                  k=FL_Wg(iflav,j)
                  if (lok(1)) then
                     res(k,1)=pdf(sign1*iflav(1),1,2*j-3) !1 for j=2;3 for j=3
     &                       *pdf(sign2*iflav(2),2,4-j)   !2 for j=2;1 for j=3
     &                       *matr3sq(i1,i2,lw,j)*polcolg
                  else
                     res(k,1) = 0
                  endif
                  if ( nlo.eq.1 .and. lok(j) ) then
                     res(k,j) = ( pdf(sign1*iflav(1),1,j)
     &                           *pdf(sign2*iflav(2),2,j)
     &               * (matr2sq(i1,i2,lw,2) + matr2sq(i1,i2,lw,3))
     &               - pdf(sign1*iflav(1),1,2*j-3)
     &               * pdf(sign2*iflav(2),2,4-j) *
     &                (matr2sq(i1,i2,lw,2)*v2/q2 +
     &                                matr2sq(i1,i2,lw,3)*v3/q3) )
     &                           * C2*polcolg
                     res(k,0) = res(k,1)+res(k,j)
                     res(k,5-j)=0
                  else
                     res(k,0) = res(k,1)
                     res(k,2) = 0
                     res(k,3) = 0
                  endif
c 116              continue
                  if (ldebug) then                                       !debug
c                     if (j.eq.2) then
c                        print 117," iflav = ",iflav(1),iflav(3),
c     &                                IFLAV(5),iflav(2),iflav(4)
c                     else
c                        print 118," iflav = ",iflav(1),iflav(3),
c     &                                iflav(2),iflav(4),IFLAV(5)
c                     endif
 119                 format(a,3i3,a,3g10.3)
                     print 119," called matrix element ",i1,i2,lw," j:",
     &  matr2sq(i1,i2,lw,2),matr2sq(i1,i2,lw,3),matr3sq(i1,i2,lw,j)
                     if (lok(1) .and. nlo.eq.1) then                     !debug
                        xa = res(k-1,1)+res(k,1)                         !debug
                        print*,k,xa," up,low,tot",                       !debug
     &                         res(k-1,2)/xa,res(k,3)/xa,                !debug
     &                         (res(k-1,2)+res(k,3))/xa                  !debug
                     endif                                               !debug
                  endif                                                  !debug
               enddo            !if2
            enddo               !if1
         enddo                  !lw
      endif                     !initial vs. final gluon

      nmax = FL_Wg(iflav,-2)

      end

c******************************************************************
c
c   begin function FL_Wg
c
c*****************************************************************
      integer function FL_Wg(iflav,colstruc)
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

      if (colstruc.eq.-1) then
c       write(*,*) "we are resetting"
         listposition=0
         numresets=numresets+1
      else if (colstruc.ne.-2) then
c  don"t do anything for colstruc=-2, this call simply 
c  returns a value for numcolor. Otherwise: incremeilnt 
c  the counter regardless of whether or not it"s the 1st 
c  time through.
         listposition=listposition+1

c  and fill info for this subprocess during first pass 
         if(numresets.eq.0) then
            call fillColoredPartons_Wg(iflav,
     &                                listposition+numdecay,colstruc)
         endif        ! numresets.eq.0
      endif
      FL_Wg=listposition
      end
c*****************************************************************************
c
c   end function FL_Wg
c
c*****************************************************************************

c*****************************************************************************
c
c    begin  subroutine fillColoredPartons_Wg
c
c*********************************************************************
      subroutine fillColoredPartons_Wg(iflav,listposition,colstruc)
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
      integer iflav(5),id1,id2,id3,id4,id5,listposition
c  local variables
      integer i,i1
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
      case(WPZjjjLO,WMZjjjLO,WPhadZjjjLO,WPZhadjjjLO,WMhadZjjjLO,WMZhadjjjLO) !W+Zjj, W-Zjj
         lnup(listposition)=numParticles+3
      case(WPAjjjLO,WMAjjjLO) ! W+Ajj, W-Ajj
         lnup(listposition)=numParticles+1
      case(WPjjjLO,WMjjjLO) ! W+jj, W-jj
         lnup(listposition)=numParticles
      end select

c  lidprup is same for all subprocesses

      listup(1,listposition)=-1   !incoming quarks
      listup(2,listposition)=-1
      listup(3,listposition)=1    !outgoing quarks
      listup(4,listposition)=1
      listup(5,listposition)=1

c  listup etc. for W and decay leptons is set in subroutine fillColorless
c  which is called in decay setup

      lidup(1,listposition)=id1     ! the 4 quarks
      lidup(2,listposition)=id2
      lidup(3,listposition)=id3
      lidup(4,listposition)=id4
      lidup(5,listposition)=id5     ! and the gluon

      do i=3,n_p
        lmothup(1,i,listposition)=1
        lmothup(2,i,listposition)=2
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


      end  ! fillColoredPartons_Wg
c*****************************************************************************
c
c    end subroutine fillColoredPartons_Wg
c
c*****************************************************************************
 

      subroutine resprintw(nmin,nmax,res)

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
      if (nmin.le.125) then 
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
      if (nmax.lt.-129) then
         do k = nmin,nmax
            print*," res(",k,") ",(res(k,i),i=1,3)
         enddo
      endif
      write(6,10)" res([",nmin,",",nmax,"])",xa, !resold(ic),resold(ic+1),
     1         xb/xa,xc/xa,xd/xa
c     1         xa/(resold(ic)+resold(ic+1))*fac
 10   format(a,i3,a,i3,a,3g11.3,3f18.10)
      ic = ic+2
      end
