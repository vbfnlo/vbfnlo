
c
c phase space generator for pp --> H + H + (np-nx) jets
c
c       Julien Baglio, <julien.baglio@kit.edu>
c       after modifications on a file created by
c	Dieter Zeppenfeld, <dieter@pheno.physics.wisc.edu>
c	Initial version:  1995 March 25 (original file)
c	Last modified:  2005 Dec. (original fine)
c       Last modified: 2014 Jul. (this routine)       
c
      subroutine def_ps_hhnj(bos0,bos1,bos2, iwidth0, iwidth1, iwidth2)
      implicit none 
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/cuts.inc"
      integer N1, bos0, bos1, bos2, bos(0:2), i
      integer iwidth0, iwidth1, iwidth2 
      parameter ( N1=max_p-nx )
c
c set up variables in common block BLIPSIvvnj for call of LIPSNHH
c 
      DOUBLE PRECISION RM2,RMG,RM2MIN,RM2MAX, S,M2MIN, YJMIN,YJMAX,PTJMIN,EJMIN
      INTEGER INFOJ, iwidth
      COMMON /BLIPSIvvnj/ RM2(0:2),RMG(0:2),RM2MIN(0:2),RM2MAX(0:2),
     1                S,M2MIN,YJMIN(N1),YJMAX(N1),PTJMIN(N1),EJMIN(N1),
     2                INFOJ(N1), iwidth(0:2)
      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)
      
      DOUBLE PRECISION xg ! parameter for simulation of VV distribution
      integer bos_var
      common /bosin/ xg,bos_var,bos

      s = ecm**2
C Choose cutoff on diparton invariant mass and Q^2 for initial-final 
C collinear cutoff: for numerical purposes only
      m2min = 0.1d0**2

      bos(0) = bos0
      bos(1) = bos1
      bos(2) = bos2
      
      bos_var = bos0      
      xg = 0d0	! for initialization only (overruled for bos0=-6)
      
c set iwidth(i) = 0 to generate distributions in zero width approcimation.
c the amplitude squared should have the normal Breit-Wigner factor included, 
c which then gets evaluted at q^2 = m^2, i.e. the program expects a factor 
c 1/(mGamma)^2 * |M_decay|^2 from the matrix element squared
      iwidth(0) = iwidth0
      iwidth(1) = iwidth1
      iwidth(2) = iwidth2
C set resonances: 
c       bos0 = 0    ==>  broad H1 + H2 mass distribution

      if(bos0.eq.0) then
         do i = 1,2
            rm2(i) = xm2(bos(i))
            rmg(i) = xmg(bos(i))
c
            if (bos(i).eq.2) then
               rm2min(i) = 1d-3 
               rm2max(i) = s/2D0
            elseif (bos(i).eq.3.or.bos(i).eq.4)then
               rm2min(i) = 1d-3 
               rm2max(i) = s/2D0
            else
               rm2min(i) = 1d-3
               rm2max(i) = s/2D0
            endif
         enddo
cc         rm2(0) = 500.**2
cc         rmg(0) = 500*500  !300.*500.
         rm2(0) = (sqrt(xm2(bos(1)))+sqrt(xm2(bos(2))))**2 ! intermediate particle, continuum
         rmg(0) = sqrt(rm2(0))*500d0
cc         rm2min(0) = (sqrt(rm2min(1))+sqrt(rm2min(2)))**2
         rm2min(0) = 1d-3    !min( (sqrt(xm2(6))+10d0)**2, 4d0*xm2(2) )	
c comment: check below if it changes between s and s/2
         rm2max(0) = s/2  !*0.95

      else
         print*," bos0 = ",bos0," not implemented in ps_hhnj.f "
         stop
      endif

      do i = 1,N1
         yjmin(i) = 0d0
         yjmax(i) = ymax_def
         ptjmin(i) = ptmin_def(i)
         ejmin(i) = 0d0
         infoj(i) = -1
      enddo
 60   format ( " machine energy = ", f10.1, "  GeV" )
 70   format ( " integration range for boson ",i2," = ",f10.3,
     1         "  to ",f12.3," GeV" )
 72   format ( " boson ",i2," has mass ",f8.2," GeV and width ",
     1         f10.4," GeV" )
 80   format ( " Minimal mom.transfer or mass for parton pair: ",f8.2,
     1         " MeV" )
      write(6,*)
      write(6,*) "       phase space generator information  "
      write(6,*) "-----------------------------------------------"
      write(6,*)
      do i = 0,2
         write(6,72) bos(i),sqrt(rm2(i)),rmg(i)/sqrt(rm2(i))
         write(6,70) bos(i),sqrt(rm2min(i)),sqrt(rm2max(i))
      enddo
      write(6,*)
c      write(6,10) ptjmin
c      write(6,20) ejmin
c      write(6,30) yjmin
c      write(6,40) yjmax
c      write(6,50) infoj
      write(6,60) sqrt(s)
      write(6,80) sqrt(m2min)*1000.d0
      write(6,*)
      write(6,*) "-----------------------------------------------"
      write(6,*)
      
      end

c
c*****************************************************************************
c
c    begin subroutine LIPSNHH
c
c*****************************************************************************
c
c       Michael Rauch, <rauch@particle.uni-karlsruhe.de>
c       Initial version: May 2012 
c       Last modified by Michael Rauch: May 2012
c       Last modified by Julien Baglio: July 2014
c       Based on LIPSVAJ by Michael Rauch
c
c phase space generator for the production of 2 Higgs + 2 jets (+ 1 jet)
c at hadron colliders
c
c  INPUT:   RD         random number array (between 0 and 1)
c           RN         additional random number which will not be driven by
c                      VEGAS
c           OS         integer flag to manage the decay of the H pair
c
c  OUTPUT:  K1,K2      4-vectors for 2 incoming partons
c           X1,X2      corresponding values of Feynman x
c           V(0:3,2)   4-vectors for the 2 final-state Higgs
c                      Q(*,1) = V(*,1) + V(*,2)
c           Pjets(0:3,N)  4-vector for the jets, assumed massless
c           W          the weight of the phase space point
c
c****************************************************************************

      SUBROUTINE LIPSNHH(N,RD,RN,
     1                   K1,K2,D,Pjets,
     2                   X1,X2,W,OS)
      IMPLICIT NONE
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/process.inc"

      INTEGER N,N1,OS
      REAL*8 TPI
c      PARAMETER (nrandl=12, TPI=2D0*PI )
c      PARAMETER (nrandl=16, TPI=2D0*PI )
      parameter (TPI=2.d0*PI)
      parameter ( N1=max_p-nx )
c      REAL*8 RD(nrandl),RN
      REAL*8 RD(*),RN
      REAL*8 K1(0:3), K2(0:3), D(0:3,4), X1, X2, W
      REAL*8 DTEMP(0:3,2)
      REAL*8 Q(0:4)
      REAL*8 Pjets(0:3,N)

c
c set up variables in common block BLIPSIvvnj for call of LIPSNHH
c 
      DOUBLE PRECISION RM2,RMG,RM2MIN,RM2MAX, S,M2MIN, YJMIN,YJMAX,PTJMIN,EJMIN
      INTEGER INFOJ, iwidth
      COMMON /BLIPSIvvnj/ RM2(0:2),RMG(0:2),RM2MIN(0:2),RM2MAX(0:2),
     1                S,M2MIN,YJMIN(N1),YJMAX(N1),PTJMIN(N1),EJMIN(N1),
     2                INFOJ(N1), iwidth(0:2)
      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)

C  local variables


      LOGICAL ldebug
      PARAMETER (ldebug=.false.)

      LOGICAL NONRESONANCE, RESONANCE, TWOBODYDECAY, TWOTOJETSPLUSX
      EXTERNAL NONRESONANCE, RESONANCE, TWOBODYDECAY, TWOTOJETSPLUSX

c--------------------------------------------------------------------
c---------------------------------------------------------------------

      if (ldebug) print*," enter LIPSNHH "

      W = 1d0             !init of weightfactor

c      write(*,*) "XXX rnd", (rd(i),i=1,12)
c generate q^2 of intermediate pseudo-state
      if (.not. Resonance(rm2(0), rmg(0), rm2min(0), rm2max(0),
     1                     rd(1), W, q(4))) return
c      write(*,*) "XXX w1", w,rm2, rmg

      if(OS.eq.1) then ! H pair on-shell
c generate 2 -> 1+N process for k1 k2 -> Pjets q_X from rd(4-9/12)
         if (.not.TwoToJetsPlusX(N, rd(4), 0d0, ecm, pTjmin, q(4),
     1        k1(0), k2(0), x1, x2, q(0), pjets(0,1), W)) return
c      write(*,*) "XXX w2", w, q

c-------------------------------------------------------------------

c generate phase space for decay of X -> H1 + H2 out of rd(2-3)
         if (.not. TwoBodyDecay(rd(2), rd(3), q(0), q(4), rm2(1),
     1         rm2(2), d(0,1),d(0,2),W)) return


      else ! H pair is decaying in Narrow-Width approximation
c generate 2 -> 1+N process for k1 k2 -> Pjets q_X from rd(8-13/16)
         if (.not.TwoToJetsPlusX(N, rd(8), 0d0, ecm, pTjmin, q(4),
     1        k1(0), k2(0), x1, x2, q(0), pjets(0,1), W)) return

c-------------------------------------------------------------------

c generate phase space for decay of X -> H1 + H2 out of rd(2-3)
         if (.not. TwoBodyDecay(rd(2), rd(3), q(0), q(4), rm2(1),
     1        rm2(2), dtemp(0,1),dtemp(0,2),W)) return
c generate phase space for decay of H1 -> l1 + l2 out of rd(4-5)
         call TwoBodyDecay0(rd(4), rd(5), dtemp(0,1), rm2(1),
     1        d(0,1),d(0,2),W)
c generate phase space for decay of H2 -> l3 + l4 out of rd(6-7)
         call TwoBodyDecay0(rd(6), rd(7), dtemp(0,2), rm2(2),
     1        d(0,3),d(0,4),W)
c- rescale the weight by the factor (8*pi^2) coming from the decay matrix elements
         W = W * (4.d0*TPI)**2


      endif ! end OS test
c      write(*,*) "XXX w3", w

c -----------------------------------------------------

c$$$      if (ldebug) then
c$$$         print*, "--------- new event -----------"
c$$$         print*, ""
c$$$         print*, "random numbers"
c$$$         print*, rn
c$$$         print*, rd(1), rd(2), rd(3)
c$$$         print*, rd(4), rd(5), rd(6)
c$$$         print*, rd(7), rd(8), rd(9)
c$$$         if (N.eq.1) print*, rd(10), rd(11), rd(12)
c$$$         print*, "y     = ", y
c$$$         print*, "tau   = ", tau
c$$$         print*, "x1    = ", x1
c$$$         print*, "x2    = ", x2
c$$$         print*, "wgt   = ", w
c$$$         print*, "---------------------------"
c$$$
c$$$         print*, "--- momenta conservation ? ---"
c$$$         do mu=0,3
c$$$            test = k1(mu)+k2(mu)-Pjet(mu,1)-q(mu,1)-q(mu,2)
c$$$            if (test.gt.1E-9) print*, test
c$$$            test = q(mu,1)-v(mu,1)-v(mu,2) 
c$$$            if (test.gt.1E-9) print*, test
c$$$            test = k1(mu)+k2(mu)-Pjet(mu,1)-v(mu,1)-v(mu,2)
c$$$     &             -v(mu,3)
c$$$             if (test.gt.1E-9) print*, test
c$$$         enddo
c$$$         print*,"tau*s = (k1+k2)^2"
c$$$         test = tau*s - dotrr(k1,k1)-dotrr(k2,k2)-2*dotrr(k1,k2)
c$$$         print*, test
c$$$         print*,"---tau*s = (q1 + q2 + q3 + pjet)^2 ?---"
c$$$         test = tau*s - dotrr(q(0,1),q(0,1)) - dotrr(q(0,2),q(0,2))
c$$$     &     - 2*dotrr(q(0,1),q(0,2))
c$$$     &     - dotrr(Pjet(0,1),Pjet(0,1)) - 2*dotrr(Pjet(0,1),q(0,1))
c$$$     &     - 2*dotrr(Pjet(0,1),q(0,2)) 
c$$$         print*, test
c$$$         print*, "--- mass equal to zero ? ---"
c$$$         print*,"p1^2 = ",dotrr(k1(0),k1(0))
c$$$         print*,"p2^2 = ",dotrr(k2(0),k2(0))
c$$$         print*,"Pjet^2 = ",dotrr(Pjet(0,1),Pjet(0,1))
c$$$         print*,"v1^2 = ",dotrr(v(0,1),v(0,1))
c$$$         print*,"v2^2 = ",dotrr(v(0,2),v(0,2))
c$$$         print*,"v3^2 = ",dotrr(v(0,3),v(0,3))
c$$$         print*, "-----------------------"
c$$$         read(*,*)
c$$$      endif

      return
      end

c*****************************************************************************
c
c    end subroutine LIPSNHH
c
c*****************************************************************************
