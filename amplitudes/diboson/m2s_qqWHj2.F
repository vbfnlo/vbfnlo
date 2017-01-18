*****************************************************************************
c
c    begin function m2s_qqWHj2
c
c*****************************************************************************
c
c       Robin Roth <robin@particle.uni-karlsruhe.de>
c       Adapted from m2s_qqWAAj2 by Michael Rauch
c       Initial version: November 2012      
c       Last modified November 2013
c
c       this function computes |ME|^2 for qq-> ve(ve~) e+(-) H j ("W+-Hj")
c	at LO and virtual corrections via call of 'qqWHj.f'
c
c ------------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION m2s_qqWHj2(
     &      bos,  !in:  XXX Boson identifier,X=1 photon, X=2 Z, X=3,4 W+-
     &      nlo,  !in:  Born only for nlo=0, +interference for 1
     &      xi,   !in:  Feynman x parameters of incoming fermions
     &      p,    !in:  momenta associated with external fermions
     &      v,    !in:  momenta associated with gauge boson dec. fermions
     &      rn,    !in:  random number array
     &      NW )

      use globalvars, only: ldoblha

      IMPLICIT NONE

c declare global variables

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/lha.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

c declare input/output variables

      !TODO: remove nw, rn(1)?
      INTEGER bos, nlo,NW
      DOUBLE PRECISION xi(nx), p(0:3,max_p), v(0:3,n_v,3), rn(1), m2s

c declare external functions

      DOUBLE PRECISION dotrr
      DOUBLE COMPLEX dotrc
      INTEGER FL_VVg
      EXTERNAL dotrr, dotrc
      EXTERNAL FL_VVg

c declare local variables

c  helicity selection
 
c     INTEGER h
c     COMMON /hcount / h
c     COMMON /icount / ie,iu,it

      DOUBLE PRECISION zero
      PARAMETER (zero=0d0)

      DOUBLE PRECISION q_sf 
      INTEGER init/0/, I,  mu 
      SAVE init

c declare variables for summation process

      INTEGER ires

      INTEGER FSIGN(6), gsign, physToDiag(3), nmin, nmax
      SAVE nmax

c store contributions from subprocess ip in res(ip) where

      DOUBLE PRECISION res(maxnumsubproc)

c declare parton distribution variables

      DOUBLE PRECISION pdf(-6:6,2)
      DOUBLE PRECISION pdfckm(-6:6,2)

c variable to switch between whj and wj
      logical hashiggs

      ! Higgs propagator factor
      double precision fac
      ! constants needed for fac from common block
      DOUBLE PRECISION CLR,XM2,XMG,B
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6)

c define program switches

      LOGICAL ldebug
      PARAMETER (ldebug=.false.)

c variables for leshouches interface
      DOUBLE PRECISION weight,rnumb,RandomNumber

c ------------- BEGIN CODE -------------------

      if (ldebug) then      
         write(6,*) "------------from m2s_qqWHj2.F-----------------"
         write(6,*) " "
      endif

      if (bos == 36 .or. bos == 46) then
         hashiggs = .true.
      else
         hashiggs = .false.
      endif

      if ( init .eq. 0 .and. .not. ldoblha) then
            write(6,*) " "
         if (bos.eq.36) then
            write(6,*) "W+Hj amplitude square information:"
         elseif (bos.eq.46) then
            write(6,*) "W-Hj amplitude square information:"
         elseif (bos.eq.3) then
            write(6,*) "W+j amplitude square information:"
         elseif (bos.eq.4) then
            write(6,*) "W-j amplitude square information:"
         endif
            write(6,*) "--------------------------------------------"
            write(6,*) " "


         call printnfl(.true.)
         init = 1
      else
      endif

      res(:) = 0
      call wtowh_reset !implemented in m2s_qqWH

      if ( ldebug ) then
         write(6,*) " m2s_qqWHj2.F, part1 -------------------------- "
 10      format( " p(", i1, ") = ", 4(f10.3, 2x) )
 20      format( " v(", i1, ") = ", 4(f10.3, 2x) )

         write(6,*) " "
         do i = 1, n_v
            write(6,20) i, v(0,i,1), v(1,i,1), v(2,i,1), v(3,i,1)
         enddo
         print*," parton momenta "
         do i = 1, max_p-1
            write(6,10) i, p(0,i), p(1,i), p(2,i), p(3,i)
         enddo

         write(6,*) " "
         write(6,*) "xi(1), xi(2) =", xi
         write(6,*) "mu_f =",sqrt(mufsq(1,1))
      endif !ldebug

      do mu = 0,3
         qw0(mu,1) = v(mu,1,1)+v(mu,2,1)   !W+-
         if(hashiggs) then
            qh(mu,1) = 0d0              !H
            do i=3,n_v
               qh(mu,1) = qh(mu,1)+v(mu,i,1)
            enddo
            qwh(mu,1) = qh(mu,1)+qw0(mu,1)
         else !only W
            qwh(mu,1) = qw0(mu,1)
         endif
      enddo


      qw0(4,1) = dotrr(qw0(0,1),qw0(0,1))
      if (hashiggs) then
         qh(4,1) = dotrr(qh(0,1),qh(0,1))
      endif
      qwh(4,1) = dotrr(qwh(0,1),qwh(0,1))

c for WH->2lH precalculate leptonic tensors

c lepton spinors and W+, A polarisation vectors
 
!     if (bos.eq.3) then
!        CALL OXXXXX(v(0,1,1),ZERO ,-1,1,wve(1,1))           !ve 
!        CALL IXXXXX(v(0,2,1),ZERO ,1,-1,wep(1,1))           !e+     
!     elseif (bos.eq.4) then
!        CALL IXXXXX(v(0,2,1),ZERO ,1,-1,wve(1,1))          !ve~  
!        CALL OXXXXX(v(0,1,1),ZERO ,-1,1,wep(1,1))            !e-    
!     endif

      if(hashiggs) then
         CALL SXXXXX(qh(0,1),1,hig(1,1))                     !H
      endif

      if(hashiggs) then
         fac = 1.d0 / ( (qh(4,1)-xm2(6))**2 + xmg(6)**2 ) !propagator^2 only
         fac = fac *16*pi*xmg(6)
      else
         fac = 1d0
      endif
      
c  original      
      CALL OXXXXX(v(0,1,1),ZERO ,-1,1,wve(1,1))           !ve 
      CALL IXXXXX(v(0,2,1),ZERO ,1,-1,wep(1,1))           !e+     


      if ( hashiggs ) then
      if (with_anom) then ! anomalous gauge boson couplings
         if (bos.eq.36) then
            call jwpmtowh_an(v(0,1,1),1,wtowh(0,1),3)  !W+ -> ve e+ h
         elseif (bos.eq.46) then
            call jwpmtowh_an(v(0,1,1),1,wtowh(0,1),4)  !W- -> ve e- h
         endif
      else !SM
         if (bos.eq.36) then
            call jwpmtowh(v(0,1,1),1,wtowh(0,1),3)  !W+ -> ve e+ h
         elseif (bos.eq.46) then
            call jwpmtowh(v(0,1,1),1,wtowh(0,1),4)    
         endif
      endif !anom or SM
      CALL JIOXXX(wep(1,1),wve(1,1),GWF,WMASS,WWIDTH, wtoonlyw(0,1))!W+- 
      else
         CALL JIOXXX(wep(1,1),wve(1,1),GWF,WMASS,WWIDTH, wtowh(0,1))!W+- 
      endif

c if desired, print output information

      if ( ldebug ) then
c10      format( " p(", i1, ") = ", 4(f10.3, 2x) )
c20      format( " v(", i1, ") = ", 4(f10.3, 2x) )

         write(6,*) " "
         do i = 1, n_v
            write(6,20) i, v(0,i,1), v(1,i,1), v(2,i,1), v(3,i,1)
         end do
         print*," parton momenta "
         do i = 1, max_p-1
            write(6,10) i, p(0,i), p(1,i), p(2,i), p(3,i)
         enddo

         write(6,*) " "
         write(6,*) "xi(1), xi(2) =", xi
         write(6,*) "mu_f =",sqrt(mufsq(1,1))
         write(6,*) "sub_number ", sub_number
      endif !ldebug


      if (ldoblha) then
         physToDiag(1:3)=blha_physdiagmap(1:3,blha_cursubproc,blha_curproc)
         fsign(1:5) = blha_fsign(1:5,blha_cursubproc,blha_curproc)
         gsign = fsign(3)

         ! fsign(1) = 1          !incoming quark
         ! fsign(2) = -1         !incoming antiquark
         ! fsign(3) = 1          !outgoing ve
         ! fsign(4) = -1         !outgoing ep
         ! fsign(5) = 1          !outgoing H !TODO: not needed
         ! gsign    = 1          !outgoing gluon

         call qq_WHj2(xi,p,v,physToDiag,fsign,gsign,bos,nlo,
     &              pdf,pdfckm,res,nmin,nmax)
         return
      endif


c
c call PDF subroutine in order to determine parton
c distributions in the incoming (anti)protons.
c
      q_sf = sqrt(mufsq(1,1))

      call pdfproton( xi(1), q_sf, pdf(-6,1) )        ! f_a(x1)
      call pdfproton( xi(2), q_sf, pdf(-6,2) )        ! f_a(x2)

      call pdfckmtransform(pdf(-6,1),pdfckm(-6,1))

C********************************************************************

C     Now sum over the subprocesses contributing to W+ H production

      nmax = FL_VVg(fsign,-1)!reset counter for subprocesses to 0

      fsign(6) = 1 !TODO: reduce fsign size to 5 (or 4), check for array bounds

C*******************  q1 qbar2 ---> W H g   **********************

c   physToDiag(ext.momentum label) = Feynman diagram label

      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=2
      physToDiag(3)=3    ! gluon

      fsign(1) = 1          !incoming quark
      fsign(2) = -1         !incoming antiquark
      fsign(3) = 1          !outgoing ve
      fsign(4) = -1         !outgoing ep
      fsign(5) = 1          !outgoing H !TODO: not needed
      gsign    = 1          !outgoing gluon


c      if (sub_number.eq.1) then
         call qq_WHj2(xi,p,v,physToDiag,fsign,gsign,bos,nlo,
     &              pdf,pdfckm,res,nmin,nmax)
c      endif

C*******************  qbar2 q1 ---> W H g   **********************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=1
      physToDiag(3)=3    ! gluon

      fsign(1) = 1
      fsign(2) = -1
      gsign    = 1
c      if (sub_number.eq.2) then
         call qq_WHj2(xi,p,v,physToDiag,fsign,gsign,bos,nlo,
     &              pdf,pdfckm,res,nmin,nmax)
c      endif


C*******************  qbar2 g ---> W H qbar1   **********************

      physToDiag(1)=2
      physToDiag(2)=3
      physToDiag(3)=1

      fsign(1) = -1
      fsign(2) = -1
      gsign    = -1
c      if (sub_number.eq.3) then
         call qq_WHj2(xi,p,v,physToDiag,fsign,gsign,bos,nlo,
     1              pdf,pdfckm,res,nmin,nmax)
c      endif
C*******************  g q2bar ---> W H q1bar  **********************

      physToDiag(1)=3
      physToDiag(2)=2
      physToDiag(3)=1

      fsign(1) = -1
      fsign(2) = -1
      gsign    = -1
c      if (sub_number.eq.4) then
         call qq_WHj2(xi,p,v,physToDiag,fsign,gsign,bos,nlo,
     1              pdf,pdfckm,res,nmin,nmax)
c      endif

C*******************  g q1 ---> W H q2  ***********************************

      physToDiag(1)=3
      physToDiag(2)=1
      physToDiag(3)=2

      fsign(1) = 1
      fsign(2) = 1
      gsign    = -1
c      if (sub_number.eq.5) then
         call qq_WHj2(xi,p,v,physToDiag,fsign,gsign,bos,nlo,
     1              pdf,pdfckm,res,nmin,nmax)
c      endif

C*******************  q1 g ---> W H q2   **********************

      physToDiag(1)=1
      physToDiag(2)=3
      physToDiag(3)=2

      fsign(1) = 1
      fsign(2) = 1
      gsign    = -1
c      if (sub_number.eq.6) then
         call qq_WHj2(xi,p,v,physToDiag,fsign,gsign,bos,nlo,
     1              pdf,pdfckm,res,nmin,nmax)
c      endif
c*****************  end of process evaluation  *********************


      if (init.eq.1) then
         init = init+1
      endif

      m2s = 0
      DO IRES = 1,NMAX
         m2s = m2s + RES(IRES)
      ENDDO

      if (ldebug) then
         print*, 'nmax', nmax
         print*, 'm2s', m2s
      endif

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
      endif

      m2s_qqWHj2 = m2s*fac
      
      RETURN
      END

c*****************************************************************************
c
c    end function m2s_qqWHj2
c
c*****************************************************************************

c******************************************************************
c
c   begin subroutine qq_WHj2
c
c*****************************************************************

      subroutine qq_WHj2(xi,p,v,physToDiag,fsign,gsign,bos, nlo,
     1                  pdf,pdfckm,res,nmin,nmax)

      use globalvars, only: ldoblha

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"


      DOUBLE PRECISION p(0:3,n_p), v(0:3,n_v), xi(nx)
      DOUBLE PRECISION pdf(-6:6,2), pdfckm(-6:6,2), res(maxnumsubproc)
      INTEGER physToDiag(3), fsign(6), gsign, bos, nlo, nmin, nmax
c  INPUT:  p(0:3,4)        external physical parton momenta
c          v(0:3,nv)       leptons and photon momenta
c          physToDiag(4)   physToDiag(ext.mom. label) = Feynman diagram label
c          fsign           sign factors for amplitude calls; see qqZqq
c          nlo             nlo = 0: calculate LO only
c                          nlo = 1: calculate full NLO (Born + interference)
c  OUTPUT: res(k)
c
c  In and Output:
c          nmin, nmax      range of process numbers for this call
c
      DOUBLE PRECISION N                              ! color factors
      PARAMETER(N=3d0)
c
c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/coupl.inc"

c local variables
      DOUBLE PRECISION udwh, duwh, udwh2, duwh2
      DOUBLE PRECISION pbar(0:3,6), qbar(0:4), dotrr, polcolq, polcolg
      DOUBLE PRECISION NCmatrixelt(0:1)
      DOUBLE COMPLEX zero
      INTEGER iflav(4), diagToPhys(3), FL_VVg
      EXTERNAL dotrr, FL_VVg

      LOGICAL ldebug, ldebugm
      PARAMETER (ldebug=.false.)
      PARAMETER (ldebugm=.false.)
      PARAMETER (zero = (0d0,0d0) )
      INTEGER i,if1,j,k,mu

      ! variables used for pdf ckm calls
      integer f1, f2
      double precision getpdfckm
      
      ! switch between WHj and Wj
      logical hashiggs

      double precision vv

      !for madgraph comparison of gg-h-vertex
      double precision  clr, xm2, xmg, b, cpl_v, a
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     1                  cpl_V(4,5),A(4,5)
      udwh = 0d0
      duwh = 0d0
      udwh2 = 0d0
      duwh2 = 0d0 
    
      if(ldebug) then 
         print*, "--------- enter m2s_qqWHj2  qq_whj2----------------"
c         print*, "lok ", lok
         print*, "nlo ", nlo
c         print*, "xuz ", xuz
      endif

      if (bos.eq.36 .or. bos.eq.46) then
         hashiggs = .true.
      elseif (bos.eq.3 .or. bos.eq.4) then
         hashiggs = .false.
      else
         print*, 'unknown bosid ', bos
         stop
      endif

      nmin = nmax+1
      do i = 1,3
         diagToPhys(physToDiag(i)) = i
      enddo

      polcolq = 1d0/(4d0*N**2*xi(1)*xi(2))
      polcolg = 1d0/(4d0*N*(N**2-1)*xi(1)*xi(2))
      if (physToDiag(1).eq.3) then
         j=2           ! j=2: emission off 1,2 line
      elseif (physToDiag(2).eq.3) then
         j=3           ! j=3: emission off 3,4 line
      endif

      do mu = 0,3
         do i = 1,3
            pbar(mu,physToDiag(i))=p(mu,i)
         enddo
         qbar(mu) = pbar(mu,3)
      enddo
!      qbar(4)=0d0 ! massless jet ! not in here???

      do mu = 0,3             ! kinematics for WH
         pbar(mu,3) = v(mu,1) ! ve
         pbar(mu,4) = v(mu,2) ! e+
         !TODO: fix or remove
         if (hashiggs) then
            pbar(mu,5) = 0d0 ! higgs
            do i=3,n_v
               pbar(mu,5) = pbar(mu,5) + v(mu,i)
            enddo
         endif
         pbar(mu,6) = qbar(mu)
      enddo
      fsign(6) = gsign

#ifdef WITH_MADGRAPH
      if (ldebugm) then
         !set all required couplings
       vv = 2.d0*sqrt(xm2(3)*(1-xm2(3)/xm2(2)))/(-clr(2,1,-1)) ! 2*MW*SW/e !higgs vev = 246
      ! vv  = 2.d0*wmass*sw/ee
         print*, 'tmass', tmass
         print*, 'mh', sqrt(xm2(6))
         print*, 'vv', vv
         mad_tau = xm2(6)/(4d0*tmass**2)
         !mad_tau = hmass**2/(4d0*tmass**2)
         mad_series_t = 1d0 + mad_tau*7d0/30d0 + mad_tau**2*2d0/21d0 + mad_tau**3*26d0/525d0
         !mad_series_t = 1d0 
         !mad_series_p = 1d0 + mad_tau/3d0      + mad_tau**2*8d0/45d0 + mad_tau**3*4d0/35d0
!         mad_gh(1) = dcmplx( 1*g**2/4d0/PI/(3d0*PI*vv)*mad_series_t, 0) !scalarf=1
         mad_gh(1) = dcmplx( 1*als(1,1)/(3d0*PI*vv)*mad_series_t, 0) !scalarf=1
         mad_gh(2) = dcmplx( 0,0) !axialf = 0
         !mad_gh(2) = dcmplx( axialf *g**2/4d0/PI/(2d0*PI*V)*series_p, Zero)


         if (bos.eq.36) then
            call sud_wphgL(pbar,fsign,udwh2)
         else if (bos.eq.46) then
            call sdu_wmhgL(pbar,fsign,duwh2)
         endif
      endif
#endif

      if (bos.eq.3 .or. bos.eq.36) then
         call qqHWj(pbar,fsign,qbar,gsign,nlo,1,udwh,bos)
      elseif (bos.eq.4 .or. bos.eq.46) then
         call qqHWj(pbar,fsign,qbar,gsign,nlo,1,duwh,bos)
      endif

      if (ldoblha) then
         return
      endif

      if (ldebug) then
         if (bos.eq.36) then
            print*,"udwh =   ",udwh
         elseif (bos.eq.46) then
            print*,"duwh =   ",duwh
         endif
      endif

      NCmatrixelt(0)=udwh
      NCmatrixelt(1)=duwh

      if (ldebug) then 
         print*, "nfl ", nfl, (nfl/2)*2
         print*, "phystodiag", phystodiag
         print*, "fsign", fsign
      endif

      do if1=1,(nfl/2)*2,2
         if (gsign .eq. 1) then
            if (bos.eq.3 .or. bos.eq.36) then
               iflav(1)=(if1-physToDiag(1)+2)*fsign(physToDiag(1))
               iflav(2)=(if1+physToDiag(1)-1)*fsign(physToDiag(2))
               iflav(3)=0
               k=FL_VVg(iflav,1)
               f1 = sign1*iflav(1)
               f2 = sign2*iflav(2)
               res(k)=getpdfckm(pdf,pdfckm,f1,f2) 
     &               *NCmatrixelt(0)*polcolq 
!              res(k)=( pdf(sign1*iflav(1),1)*pdfckm(sign2*iflav(2),2)
!    &                 +pdfckm(sign1*iflav(1),1)*pdf(sign2*iflav(2),2))*0.5d0
!    &               *NCmatrixelt(0)*polcolq 
            elseif (bos.eq.4 .or. bos.eq.46) then
               iflav(1)=(if1+physToDiag(1)-1)*fsign(physToDiag(1))
               iflav(2)=(if1-physToDiag(1)+2)*fsign(physToDiag(2))
               iflav(3)=0
               k=FL_VVg(iflav,1)
               f1 = sign1*iflav(1)
               f2 = sign2*iflav(2)
!              res(k)=( pdf(sign1*iflav(1),1)*pdfckm(sign2*iflav(2),2)
!    &                 +pdfckm(sign1*iflav(1),1)*pdf(sign2*iflav(2),2))*0.5d0
               res(k)=getpdfckm(pdf,pdfckm,f1,f2) 
     &               *NCmatrixelt(1)*polcolq
            endif
         elseif (gsign.eq.-1) then !initial gluon section
            ! here getpdfckm gets factor 2 to compensate for 0.5 factor
            ! in definition
            if (bos.eq.3 .or. bos.eq.36) then
               iflav(4-j)=(if1+PhysToDiag(3)-1)*fsign(physToDiag(4-j))
               iflav(3)  =(if1-PhysToDiag(3)+2)*fsign(physToDiag(3))
               iflav(j-1)=0
               k=FL_VVg(iflav,j)
               f1 = sign1*iflav(1)
               f2 = sign2*iflav(2)
!              res(k)=( pdf(sign1*iflav(1),1)*pdfckm(sign2*iflav(2),2)
!    &                 +pdfckm(sign1*iflav(1),1)*pdf(sign2*iflav(2),2))*0.5d0
               res(k)=2d0*getpdfckm(pdf,pdfckm,f1,f2) 
     &               *NCmatrixelt(0)*polcolg
            elseif (bos.eq.4 .or. bos.eq.46) then
               iflav(4-j)=(if1-PhysToDiag(3)+2)*fsign(physToDiag(4-j))
               iflav(3)  =(if1+PhysToDiag(3)-1)*fsign(physToDiag(3))
               iflav(j-1)=0
               k=FL_VVg(iflav,j)
               f1 = sign1*iflav(1)
               f2 = sign2*iflav(2)
!              res(k)=( pdf(sign1*iflav(1),1)*pdfckm(sign2*iflav(2),2)
!    &                 +pdfckm(sign1*iflav(1),1)*pdf(sign2*iflav(2),2))*0.5d0
               res(k)=2d0*getpdfckm(pdf,pdfckm,f1,f2) 
     &               *NCmatrixelt(1)*polcolg
            endif
         endif
         if (ldebug) then
            print*, "iflav  ", iflav
            print*, 'res', res(1:8)
         endif
      enddo 

      if (ldebug) then
!         read(*,*)
      endif
      nmax = FL_VVg(iflav,-2)

      end
c******************************************************************
c
c   end subroutine qq_WH
c
c*****************************************************************

