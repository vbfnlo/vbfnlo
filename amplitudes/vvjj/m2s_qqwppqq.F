c
c  Initial version:  Dieter Zeppenfeld, 2003 February 12
c  Last modified by sophy@particle.uni-karlsruhe.de 2012 February 24
c     From code by Barbara Jager
c
c  this function computes |ME|^2 for qq->qq e+ve mu+vm ("W+W+")
c	at LO and virtual corrections via call of "qqwwqq.f" 
c
c ------------------------------------------------------------------------
c
      real*8 function m2s_wbfwpp(
     &                   bos,   !in:  Boson identifier, 2 = Z, 6=H, 34=W+W-
     &                   nlo,   !in:  Born only for nlo=0, +interference for 1
     &                   xi,	!in:  Feynman x parameters of incoming fermions
     &                   p,	!in:  momenta associated with external fermions
     &                   v,	!in:  momenta associated with Z/H dec. fermions
     &                   rn     !in:  random number array
     &                        )
          use globalvars, only: ldoblha
      implicit none
c
c declare input/output variables
c
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/process.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"
#include "tensor.inc"
#include "tensorz.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/tau_hel.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"

c variables for hadronic decays
      integer N_gen_W
      real*8 fac_W, mjj2
      external mjj2

      integer bos, nlo
      real*8 xi(nx), p(0:3,max_p), v(0:3,max_v), rn(1), m2s

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

      real*8 q12(0:4,3), q34(0:4,3)
c
c declare local variables
c
      real*8 q_sf, qq

      integer init/0/, I, J, mu
      save init
c
c declare variables for summation process
c
      INTEGER ires
      integer FSIGN(4+max_v), physToDiag(4), nmin, nmax, nproc(4)
      save nmax, nproc
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

c -----------------------------------------------------------------------

      if ( (init .eq. 0) .and. (.not. ldoblha) ) then

         if (bos .eq. 33) then
            write(6,*) ' '
            if (procID .eq. WPWPjj .or. procID .eq. WPhadWPjj) then
               write(6,*) 'W+W+ jj amplitude square information:'
            else
               write(6,*) 'W-W- jj amplitude square information:'
            end if
            write(6,*) '-----------------------------------------------'
            write(6,*) ' '
            call vtowpm_reset
         endif
         call printnfl(.true.)
         h = 1
         init = 1
         do i = 1,maxnumsubproc
            res(i) = 0
         enddo
      else if ( (init .eq. 0) .and. ldoblha ) then
        h = 1
        init = 1
      else ! init .ne. 0
         do i = 1,nmax
            res(i) = 0
         enddo
      endif    ! end of init

c -----------------------------------------------------------------------

      call Calc_Momentum_Transfer(p, v, q12,q34,1)


* leptonic tensors for W+W+
      if (bos .eq. 33) then   

c reset lfs to .true. to force recalculation of virtual contributions in qqwwqq
         do i = 1,4
            lfp(i) = .true.
            lfm(i) = .true.
         enddo

         if ((lha.or.hepmc).and..not.doNLO) then
            if (procID.eq.WPWPjj) then
               helicity(1)= 1
               helicity(2)=-1
               helicity(3)= 1
               helicity(4)=-1
            else
               helicity(1)=-1
               helicity(2)= 1
               helicity(3)=-1
               helicity(4)= 1
            endif
         endif

         CALL IXXXXX(v(0,1),ZERO , 1,-1,wel) !e+ 
         CALL OXXXXX(v(0,2),ZERO ,-1, 1,wve) !ve 
         CALL IXXXXX(v(0,3),ZERO , 1,-1,wmu) !mu+ 
         CALL OXXXXX(v(0,4),ZERO ,-1, 1,wvu) !vm 
         CALL JIOXXX(wel,wve,GWF ,WMASS,WWIDTH,we)    !W+->e+ve
         CALL JIOXXX(wmu,wvu,GWF ,WMASS,WWIDTH,wu)    !W+->mu+vm
    
         do mu = 0,3
            qe(mu)  = v(mu,1)+v(mu,2)
            qu(mu)  = v(mu,3)+v(mu,4)
            qww(mu) = qe(mu) + qu(mu)
         enddo
         qe(4)  = qe(0)**2-qe(1)**2-qe(2)**2-qe(3)**2
         qu(4)  = qu(0)**2-qu(1)**2-qu(2)**2-qu(3)**2
         qww(4) = qww(0)**2-qww(1)**2-qww(2)**2-qww(3)**2
	 
	 
	 ! determine "reduced" polarization vectors:
	 if (abs(nlo).eq.5 .or. nlo.eq.-4) then   
            qq = qe(0)*qww(0)-qe(1)*qww(1)-qe(2)*qww(2)-qe(3)*qww(3)
            xxew = dotrc(qww,we)/qq
            
	    qq = qu(0)*qww(0)-qu(1)*qww(1)-qu(2)*qww(2)-qu(3)*qww(3)
            xxuw = dotrc(qww,wu)/qq  
	    
            do mu = 0,3
               wet(mu+1) = we(mu+1) - xxew*qe(mu)
               wut(mu+1) = wu(mu+1) - xxuw*qu(mu)
            enddo
            do mu = 5,6
               wet(mu) = we(mu)
               wut(mu) = wu(mu)
            enddo
         elseif (nlo.eq.1) then
            do mu = 1,6
               wet(mu) = we(mu)
               wut(mu) = wu(mu)
            enddo
            
         endif                  ! end nlo = -4,+-5	 for bos=33
	
         j = 1


         SELECT CASE(procid)
         CASE(WPhadWPjj, WMhadWMjj)
         
c           using global form factor for all tensors of one phase space point
c           this ensures proper cancellations for anomalous contributions
c           energy scale is invariant WW mass
            call anomal_formfactor(q12(0,j),q34(0,j),qe(0),qu(0))

c leptonic tensors
C for W+W+->e+ve mu+vm
            call wpptowpp_had(q12(0,j),q34(0,j),v,wppwpp(0,0,j))  
 
c t-channel processes:	 
C for V1 V2->e+ve
            call WPVtoWP_had(1,1,q34(0,j),v,awen(0,0,1,j),zwen(0,0,1,j)) !Wu on upper
            call WPVtoWP_had(1,2,q12(0,j),v,awen(0,0,2,j),zwen(0,0,2,j)) !Wu on lower
c        
            call WPVtoWP_had(1,1,q12(0,j),v,waen(0,0,2,j),wzen(0,0,2,j)) !Wu on lower
            call WPVtoWP_had(1,2,q34(0,j),v,waen(0,0,1,j),wzen(0,0,1,j)) !Wu on upper

c for V1 V2->mu+vm
            call WPVtoWP_had(2,1,q34(0,j),v,awmu(0,0,1,j),zwmu(0,0,1,j)) !We on upper
            call WPVtoWP_had(2,2,q12(0,j),v,awmu(0,0,2,j),zwmu(0,0,2,j)) !We on lower
c     
            call WPVtoWP_had(2,1,q12(0,j),v,wamu(0,0,2,j),wzmu(0,0,2,j)) !We on lower
            call WPVtoWP_had(2,2,q34(0,j),v,wamu(0,0,1,j),wzmu(0,0,1,j)) !We on upper

         CASE DEFAULT

         if (with_anom) then ! anomalous gauge boson couplings
c           using global form factor for all tensors of one phase space point
c           this ensures proper cancellations for anomalous contributions
c           energy scale is invariant WW mass
            call anomal_formfactor(q12(0,j),q34(0,j),qe(0),qu(0))

c leptonic tensors
C for W+W+->e+ve mu+vm
            call wpptowpp_anomal(q12(0,j),q34(0,j),v,wppwpp(0,0,j))  
 
c t-channel processes:	 
C for V1 V2->e+ve
            call WPVtoWP_anomal(1,1,q34(0,j),v,awen(0,0,1,j),zwen(0,0,1,j)) !Wu on upper
            call WPVtoWP_anomal(1,2,q12(0,j),v,awen(0,0,2,j),zwen(0,0,2,j)) !Wu on lower
c        
            call WPVtoWP_anomal(1,1,q12(0,j),v,waen(0,0,2,j),wzen(0,0,2,j)) !Wu on lower
            call WPVtoWP_anomal(1,2,q34(0,j),v,waen(0,0,1,j),wzen(0,0,1,j)) !Wu on upper

c for V1 V2->mu+vm
            call WPVtoWP_anomal(2,1,q34(0,j),v,awmu(0,0,1,j),zwmu(0,0,1,j)) !We on upper
            call WPVtoWP_anomal(2,2,q12(0,j),v,awmu(0,0,2,j),zwmu(0,0,2,j)) !We on lower
c     
            call WPVtoWP_anomal(2,1,q12(0,j),v,wamu(0,0,2,j),wzmu(0,0,2,j)) !We on lower
            call WPVtoWP_anomal(2,2,q34(0,j),v,wamu(0,0,1,j),wzmu(0,0,1,j)) !We on upper


         else !SM

c leptonic tensors
C for W+W+->e+ve mu+vm
            call wpptowpp(q12(0,j),q34(0,j),v,wppwpp(0,0,j))  
 
c t-channel processes:	 
C for V1 V2->e+ve
            call WPVtoWP(1,1,q34(0,j),v,awen(0,0,1,j),zwen(0,0,1,j)) !Wu on upper
            call WPVtoWP(1,2,q12(0,j),v,awen(0,0,2,j),zwen(0,0,2,j)) !Wu on lower
c        
            call WPVtoWP(1,1,q12(0,j),v,waen(0,0,2,j),wzen(0,0,2,j)) !Wu on lower
            call WPVtoWP(1,2,q34(0,j),v,waen(0,0,1,j),wzen(0,0,1,j)) !Wu on upper

c for V1 V2->mu+vm
            call WPVtoWP(2,1,q34(0,j),v,awmu(0,0,1,j),zwmu(0,0,1,j)) !We on upper
            call WPVtoWP(2,2,q12(0,j),v,awmu(0,0,2,j),zwmu(0,0,2,j)) !We on lower
c     
            call WPVtoWP(2,1,q12(0,j),v,wamu(0,0,2,j),wzmu(0,0,2,j)) !We on lower
            call WPVtoWP(2,2,q34(0,j),v,wamu(0,0,1,j),wzmu(0,0,1,j)) !We on upper
        
         endif !anom couplings or SM

         END SELECT


      endif


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
    
       if (ldoblha) then

         if (procid .eq. WmWmjj) then
           physToDiag(1)=blha_physdiagmap(3,blha_cursubproc,blha_curproc) 
           physToDiag(2)=blha_physdiagmap(4,blha_cursubproc,blha_curproc)
           physToDiag(3)=blha_physdiagmap(1,blha_cursubproc,blha_curproc)
           physToDiag(4)=blha_physdiagmap(2,blha_cursubproc,blha_curproc)
c first argument of blha_fsign is diagtophys(i)!
           fsign(1) = -blha_fsign(2,blha_cursubproc,blha_curproc)
           fsign(2) = -blha_fsign(1,blha_cursubproc,blha_curproc)
           fsign(3) = -blha_fsign(4,blha_cursubproc,blha_curproc)
           fsign(4) = -blha_fsign(3,blha_cursubproc,blha_curproc)
         else
           physToDiag(1)=blha_physdiagmap(1,blha_cursubproc,blha_curproc)
           physToDiag(2)=blha_physdiagmap(2,blha_cursubproc,blha_curproc)
           physToDiag(3)=blha_physdiagmap(3,blha_cursubproc,blha_curproc)
           physToDiag(4)=blha_physdiagmap(4,blha_cursubproc,blha_curproc)
           fsign(1) = blha_fsign(1,blha_cursubproc,blha_curproc)
           fsign(2) = blha_fsign(2,blha_cursubproc,blha_curproc)
           fsign(3) = blha_fsign(3,blha_cursubproc,blha_curproc)
           fsign(4) = blha_fsign(4,blha_cursubproc,blha_curproc)
         endif

         call wbf_wpp(xi,p,v,physToDiag,fsign,bos,nlo,
     1            pdf,res,nmin,nmax)

         return
       endif

c
c call PDF subroutine in order to determine parton
c distributions in the incoming (anti)protons.
c
      q_sf = sqrt(mufsq(1,1))
      call pdfproton( xi(1), q_sf, pdf(-6,1) )        ! f_a(x1)

      q_sf = sqrt(mufsq(2,1))
      call pdfproton( xi(2), q_sf, pdf(-6,2) )        ! f_a(x2)  


c select helicity
      jsig = 0


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     Now sum over the subprocesses contributing to W+W+ production

      nmax = FL_ZH(fsign,-1)            !reset counter for subprocesses to 0

C*******************  q1 q3 ---> q2 q4  V V  **********************
         
c   physToDiag(ext.momentum label) = Feynman diagram label

      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=3
      physToDiag(3)=2    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=4
C NOTE: for call of wbf_wpp it is important that p(*,1,*) and p(*,3,*)
c correspond to 1-2 fermion line ALWAYS, i.e physToDiag(1/2)={1,3} and 
c similarly physToDiag(3/4)={2,4} for the 3-4 fermion line
      fsign(1) = 1
      fsign(2) = 1
      fsign(3) = 1
      fsign(4) = 1

      call wbf_wpp(xi,p,v,physToDiag,fsign,bos,nlo,
     1            pdf,res,nmin,nmax)
     

c      if (ldebug) call resprint2(nmin,nmax,res)
      if (init.eq.1) nproc(1) = nmax

C*******************  q1 qb4 ---> q2 qb3 V V  **********************

      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=4
      physToDiag(3)=2    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=3

      fsign(3) = -1
      fsign(4) = -1

      call wbf_wpp(xi,p,v,physToDiag,fsign,bos, nlo,
     1         pdf,res,nmin,nmax)

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
      
      call wbf_wpp(xi,p,v,physToDiag,fsign,bos, nlo,
     1         pdf,res,nmin,nmax)

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
         
      call wbf_wpp(xi,p,v,physToDiag,fsign,bos, nlo,
     1         pdf,res,nmin,nmax)

c      if (ldebug) call resprint2(nmin,nmax,res)
      if (init.eq.1) nproc(4) = nmax
 
c*****************  end of process evaluation  *********************
      if (init.eq.1) then
         init = init+1
         if (lwarn) print 199," proc #s for W+W+jj are ",nproc
 199     format(a,4i5)
      endif

      SELECT CASE(process)
      CASE(WPhadWPjj,WMhadWMjj)
        call hadDecayFactor_W(finalquarks(1), abs(mjj2(v(0,1),v(0,2))), N_gen_W, fac_W)
        DO IRES = 1,NMAX
          res(IRES) = res(IRES) * fac_W            ! factor fac_W for hadronic decay
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
         if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
            rnumb=RandomNumber()
            SELECT CASE(process)
            CASE(WPhadWPjj)
               finalquarks_psp(1) =  2 + 2* INT(rnumb*2)
               finalquarks_psp(2) = -1 - 2* INT(rnumb*2)
            CASE(WMhadWMjj)
               finalquarks_psp(1) =  1 + 2* INT(rnumb*2)
               finalquarks_psp(2) = -2 - 2* INT(rnumb*2)
            END SELECT
         endif
      endif


c      DO IRES = nproc(1)+1,nproc(3)
c         m2s = m2s + RES(IRES)
c      ENDDO
    
      m2s_wbfwpp = m2s

      RETURN
      END


c******************************************************************
c
c   begin subroutine wbf_wpp
c
c*****************************************************************
      subroutine wbf_wpp(xi,p,v,physToDiag,fsign,bos, nlo,
     1                  pdf,res,nmin,nmax)
          use globalvars, only: ldoblha

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

      real*8 p(0:3,max_p), v(0:3,max_v), xi(nx)
      real*8 pdf(-6:6,2), res(maxnumsubproc)
      integer physToDiag(4), fsign(4+max_v), bos, h,nlo, nmin, nmax
c      
      common /hcount / h
c
c wbf_wpp calls the amplitude square routines 
c             qqwppqq  for qq-->qqW+W+ for bos=33
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
      real*8 pbar(0:3,4+max_v), dotrr, polcol
      integer iflav(4), diagToPhys(4), FL_ZH
      external dotrr, FL_ZH


      integer i,if1,if2,k,mu

      integer if1min,if2min
      double precision matrsq
      real*8 ucds ! for W+W+
      real*8 dsuc ! for W-W-

      LOGICAL ldebug_ltens
      PARAMETER (ldebug_ltens=.false.)

      ucds = 0d0  ! for W+W+
      dsuc = 0d0 ! for W+W+

      matrsq = 0d0

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


      if (bos .eq. 33) then  ! W+ W+
         do mu = 0,3             ! kinematics for 4-lepton final state
            pbar(mu,5) = v(mu,1) ! l+
            pbar(mu,6) = v(mu,2) ! nu
            pbar(mu,7) = v(mu,3) ! l'+
            pbar(mu,8) = v(mu,4) ! nu'
         enddo

         fsign(5) = -1
         fsign(6) = 1
         fsign(7) = -1
         fsign(8) = 1
!         if (lwarn) then
!           call qqwpmqq_mg(pbar,fsign,bos,ucds)
!         endif

         call qqwppqq(pbar,fsign, nlo,1,bos,ucds) 

#ifdef WITH_MADGRAPH
        if (ldebug_ltens) call qqwppqq_mg2(p,pbar,v,fsign, nlo, 
     1                     bos,ucds)
#endif

c           call qqwpmqq_mg(pbar,fsign,bos,ucds)

      else
         print*," Invalid entry BOS = ",bos
         print*," in call of wbf_wpp. Must be 33"
         stop
      endif


      if (ldoblha) then
        return
      endif

      
      if (bos.eq.33) then
      	matrsq = ucds !t-channel
      endif

c FLAVOR SUMMATION CHECKED FOR W+W+ CHANNEL ONLY!	
      IF (bos.eq.33) THEN
	
      ! ifmin = 2 for up & charm (anti-)quark, 1 for down & strange (anti-)quark
         if1min = (fsign(physToDiag(1))+3)/2
         if2min = (fsign(physToDiag(2))+3)/2
  
         do if1 = if1min,int(nfl/2)*2,2
            do if2 = if2min,(nfl/2)*2,2
               iflav(1)=if1*fsign(physToDiag(1)) !in
               iflav(2)=if2*fsign(physToDiag(2)) !in
               
               iflav(3)=iflav(1)-1 !out
               iflav(4)=iflav(2)-1 !out
               
               k=FL_ZH(iflav,1)
               
               res(k)=pdf(sign1*iflav(1),1)*pdf(sign2*iflav(2),2) 
     &              *matrsq*polcol  
	   
            enddo 
         enddo	

      end if  ! if proc. = 33 (W+W+) or not

      nmax = FL_ZH(iflav,-2)

      end

c******************************************************************
c
c   end subroutine wbf_wpp
c
c*****************************************************************

** Clearing the values of the leptonic tensors for W+W+

      subroutine vtowpm_reset

      implicit none

#include "tensor.inc"
#include "tensorz.inc"

      complex*16 zero
      parameter (zero=(0d0,0d0))
      integer j,jj,mu,nu

      do mu = 0,3
	 do nu = 0,3
            do j = 1,3
               wppwpp(mu,nu,j) = zero
               wmmwmm(mu,nu,j) = zero
	       
	       do jj = 1,2
		  
     		  awen(mu,nu,jj,j) = zero
		  zwen(mu,nu,jj,j) = zero
     		  waen(mu,nu,jj,j) = zero
		  wzen(mu,nu,jj,j) = zero
     		  awmu(mu,nu,jj,j) = zero
		  zwmu(mu,nu,jj,j) = zero
     		  wamu(mu,nu,jj,j) = zero
		  wzmu(mu,nu,jj,j) = zero
		  
		enddo !jj

            enddo !j
         enddo !nu
      enddo !mu
            
      print*,' VV -> WW tensors initialized to ',zero 
      return
      end
