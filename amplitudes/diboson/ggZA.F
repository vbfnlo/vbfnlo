c****************subroutine GGZA *************************************
c
c	Michael Rauch, <rauch@particle.uni-karlsruhe.de>
c	Last modified: Jul 2011 (rauch)
c       Including virtual corrections
C
C  GGZA calculates the matrix elements**2 for 
C        g g   ---->   Z A   
C
C
C  and crossing related processes. Pauli interference terms for
c  identical fermions are neglected.
c
c          N = 1         Momentum configuration 1 (LO)
c
c*********************************************************************

      subroutine ggza(pbar,fsign, nlo, N, uuzz)

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/scales.inc"
c#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/koppln.inc"

c     electroweak couplings are taken from KOPPLN
      double precision  clr, xm2, xmg, b, v, a
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     1                  V(4,5),A(4,5)
      double precision  pbar(0:3,2+n_v)
      double precision  uuzz, res(3)
      double precision  p(0:3,2+n_v)


      double precision  eps(0:3,2,2) ! g in kartesian basis
      integer  fsign(2+n_v), sign(2+n_v), nlo, mu, i, j, k, isig, N
      integer  l1, l2


      double complex ceps(0:3,2,2) ! g in complex kartesian basis


      double complex mat(2,-1:1,2,2)
      double complex mm(2,-1:1), zero, im
      double complex m1, m2

      double complex dotcc, s1c, dotrc
      double precision dotrr, psumsq
      external dotcc, s1c, dotrr, dotrc, psumsq
      parameter (im = (0d0,1d0))
      parameter (zero = (0d0,0d0) )

c  helicity selection
      INTEGER h
      COMMON /hcount / h

c variables for virtual corrections
      double precision scale
      double precision c2,c2o4pi,oneo4pi,TR,tf
      parameter (c2=4d0/3d0, c2o4pi=c2/4d0/pi, oneo4pi=1d0/4d0/pi)
      parameter (TR=1d0/2d0,tf=1d0/2d0)
      logical lbox, ldebuglt, linit, lhiggs

c Variables for fermion loops
      integer cpl(4), Div
      double precision mf(4), vv
      double complex rfbox(6,2,4,2,2)
      double precision rgfbox(6,2,4,2,2) 
      double complex rfvert(2,3,4,2,2)


      integer*8 thirtytwo


      parameter(thirtytwo=32)

      Logical wardidtest1

      real*8 minus_qaa(0:3)
      real*8 minus_qa1z(0:3),minus_qa2z(0:3)

      data linit /.true./
      save linit
      parameter(ldebuglt = .false.)

c---------------------------- BEGIN CODE -----------------------

c whether to include diagrams with boxes or Higgs bosons
      lbox   = (floops.eq.1).or.(floops.eq.3)
      lhiggs = (floops.eq.2).or.(floops.eq.3)

      do k = 1,2
         do isig = -1,1,2
            do i=1,2
               do j=1,2
                  mat(k,isig,i,j) = (0.0d0,0.0d0)
               enddo
            enddo
         enddo
      enddo

c define the internal momenta
      do i = 1,2+n_v
         sign(i) = fsign(i)
      enddo
      do mu = 0,3
         do i = 1,2+n_v
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo
      enddo

c------------------------------------------------

      do l1 = 1,2                ! 2 gluon polarizations
         call polvec(pbar(0,1),l1,eps(0,l1,1)) ! get gluon pol.vectors
         call polvec(pbar(0,2),l1,eps(0,l1,2)) ! get gluon pol.vectors
         do mu = 0,3                   
c           eps(mu,l1,1) = pbar(mu,1)         !QCD gaugecheck
c           eps(mu,l1,2) = pbar(mu,2)         !QCD gaugecheck
           ceps(mu,l1,1)=eps(mu,l1,1)
           ceps(mu,l1,2)=eps(mu,l1,2)
         enddo
      enddo

*******************************************************
C Calculation of born matrix elements
*******************************************************

c g g -> Z A 
      scale=-2d0*dotrr(p(0,1),p(0,2))
      Div = 0

      mf(1) = 0d0 ! light up quarks
      mf(2) = 0d0 ! light down quarks
      mf(3) = xmt ! top
      mf(4) = bmass_h ! bottom
      vv = 2.d0*sqrt(xm2(3)*(1-xm2(3)/xm2(2)))/(-clr(2,1,-1)) ! 2*MW*SW/e

      cpl(1) = 3 ! light up quarks
      cpl(2) = 4 ! light down quarks
      cpl(3) = 3 ! top
      cpl(4) = 4 ! bottom

      do mu = 0,3
         minus_qa1z(mu) = -qp(mu,N)
         minus_qa2z(mu) = -qm(mu,N)
         minus_qaa(mu)  = -qaa(mu,N)
      enddo

      if (lbox) then 
c box loops
ccc first anchor point is always first gluon
ccc permutate second g and the two Z/photons over the other ones
ccc strong coupling is added later
cc g-1-2
       do j=1,4
       do l1=1,2
       do l2=1,2
         call fermionboxcouplingT(
     &     1d0,1d0,clr(cpl(j),2,+1),clr(cpl(j),1,+1),
     &     1d0,1d0,clr(cpl(j),2,-1),clr(cpl(j),1,-1),
     &     mf(j),
     &     pbar(0,1),pbar(0,2),minus_qa1z(0),minus_qa2z(0),
     &     ceps(0,l1,1),ceps(0,l2,2),wp(1,N),atau2(1,N),
     &     scale,+1,+1,rfbox(1,1,j,l1,l2),rgfbox(1,1,j,l1,l2),
     &     Div)
         call fermionboxcouplingT(
     &     1d0,1d0,clr(cpl(j),1,+1),clr(cpl(j),1,+1),
     &     1d0,1d0,clr(cpl(j),1,-1),clr(cpl(j),1,-1),
     &     mf(j),
     &     pbar(0,1),pbar(0,2),minus_qa1z(0),minus_qa2z(0),
     &     ceps(0,l1,1),ceps(0,l2,2),atau1(1,N),atau2(1,N),
     &     scale,+1,+1,rfbox(1,2,j,l1,l2),rgfbox(1,2,j,l1,l2),
     &     Div)
       enddo
       enddo
       enddo
cc g-2-1
       do j=1,4
       do l1=1,2
       do l2=1,2
         call fermionboxcouplingT(
     &     1d0,1d0,clr(cpl(j),1,+1),clr(cpl(j),2,+1),
     &     1d0,1d0,clr(cpl(j),1,-1),clr(cpl(j),2,-1),
     &     mf(j),
     &     pbar(0,1),pbar(0,2),minus_qa2z(0),minus_qa1z(0),
     &     ceps(0,l1,1),ceps(0,l2,2),atau2(1,N),wp(1,N),
     &     scale,+1,+1,rfbox(2,1,j,l1,l2),rgfbox(2,1,j,l1,l2),
     &     Div)
         call fermionboxcouplingT(
     &     1d0,1d0,clr(cpl(j),1,+1),clr(cpl(j),1,+1),
     &     1d0,1d0,clr(cpl(j),1,-1),clr(cpl(j),1,-1),
     &     mf(j),
     &     pbar(0,1),pbar(0,2),minus_qa2z(0),minus_qa1z(0),
     &     ceps(0,l1,1),ceps(0,l2,2),atau2(1,N),atau1(1,N),
     &     scale,+1,+1,rfbox(2,2,j,l1,l2),rgfbox(2,2,j,l1,l2),
     &     Div)
       enddo
       enddo
       enddo
cc 1-g-2
       do j=1,4
       do l1=1,2
       do l2=1,2
         call fermionboxcouplingT(
     &     1d0,clr(cpl(j),2,+1),1d0,clr(cpl(j),1,+1),
     &     1d0,clr(cpl(j),2,-1),1d0,clr(cpl(j),1,-1),
     &     mf(j),
     &     pbar(0,1),minus_qa1z(0),pbar(0,2),minus_qa2z(0),
     &     ceps(0,l1,1),wp(1,N),ceps(0,l2,2),atau2(1,N),
     &     scale,+1,+1,rfbox(3,1,j,l1,l2),rgfbox(3,1,j,l1,l2),
     &     Div)
         call fermionboxcouplingT(
     &     1d0,clr(cpl(j),1,+1),1d0,clr(cpl(j),1,+1),
     &     1d0,clr(cpl(j),1,-1),1d0,clr(cpl(j),1,-1),
     &     mf(j),
     &     pbar(0,1),minus_qa1z(0),pbar(0,2),minus_qa2z(0),
     &     ceps(0,l1,1),atau1(1,N),ceps(0,l2,2),atau2(1,N),
     &     scale,+1,+1,rfbox(3,2,j,l1,l2),rgfbox(3,2,j,l1,l2),
     &     Div)
       enddo
       enddo
       enddo
cc 2-g-1
       do j=1,4
       do l1=1,2
       do l2=1,2
         call fermionboxcouplingT(
     &     1d0,clr(cpl(j),1,+1),1d0,clr(cpl(j),2,+1),
     &     1d0,clr(cpl(j),1,-1),1d0,clr(cpl(j),2,-1),
     &     mf(j),
     &     pbar(0,1),minus_qa2z(0),pbar(0,2),minus_qa1z(0),
     &     ceps(0,l1,1),atau2(1,N),ceps(0,l2,2),wp(1,N),
     &     scale,+1,+1,rfbox(4,1,j,l1,l2),rgfbox(4,1,j,l1,l2),
     &     Div)
         call fermionboxcouplingT(
     &     1d0,clr(cpl(j),1,+1),1d0,clr(cpl(j),1,+1),
     &     1d0,clr(cpl(j),1,-1),1d0,clr(cpl(j),1,-1),
     &     mf(j),
     &     pbar(0,1),minus_qa2z(0),pbar(0,2),minus_qa1z(0),
     &     ceps(0,l1,1),atau2(1,N),ceps(0,l2,2),atau1(1,N),
     &     scale,+1,+1,rfbox(4,2,j,l1,l2),rgfbox(4,2,j,l1,l2),
     &     Div)
       enddo
       enddo
       enddo
cc 1-2-g
       do j=1,4
       do l1=1,2
       do l2=1,2
         call fermionboxcouplingT(
     &     1d0,clr(cpl(j),2,+1),clr(cpl(j),1,+1),1d0,
     &     1d0,clr(cpl(j),2,-1),clr(cpl(j),1,-1),1d0,
     &     mf(j),
     &     pbar(0,1),minus_qa1z(0),minus_qa2z(0),pbar(0,2),
     &     ceps(0,l1,1),wp(1,N),atau2(1,N),ceps(0,l2,2),
     &     scale,+1,+1,rfbox(5,1,j,l1,l2),rgfbox(5,1,j,l1,l2),
     &     Div)
         call fermionboxcouplingT(
     &     1d0,clr(cpl(j),1,+1),clr(cpl(j),1,+1),1d0,
     &     1d0,clr(cpl(j),1,-1),clr(cpl(j),1,-1),1d0,
     &     mf(j),
     &     pbar(0,1),minus_qa1z(0),minus_qa2z(0),pbar(0,2),
     &     ceps(0,l1,1),atau1(1,N),atau2(1,N),ceps(0,l2,2),
     &     scale,+1,+1,rfbox(5,2,j,l1,l2),rgfbox(5,2,j,l1,l2),
     &     Div)
       enddo
       enddo
       enddo
cc 2-1-g
       do j=1,4
       do l1=1,2
       do l2=1,2
         call fermionboxcouplingT(
     &     1d0,clr(cpl(j),1,+1),clr(cpl(j),2,+1),1d0,
     &     1d0,clr(cpl(j),1,-1),clr(cpl(j),2,-1),1d0,
     &     mf(j),
     &     pbar(0,1),minus_qa2z(0),minus_qa1z(0),pbar(0,2),
     &     ceps(0,l1,1),atau2(1,N),wp(1,N),ceps(0,l2,2),
     &     scale,+1,+1,rfbox(6,1,j,l1,l2),rgfbox(6,1,j,l1,l2),
     &     Div)
         call fermionboxcouplingT(
     &     1d0,clr(cpl(j),1,+1),clr(cpl(j),1,+1),1d0,
     &     1d0,clr(cpl(j),1,-1),clr(cpl(j),1,-1),1d0,
     &     mf(j),
     &     pbar(0,1),minus_qa2z(0),minus_qa1z(0),pbar(0,2),
     &     ceps(0,l1,1),atau2(1,N),atau1(1,N),ceps(0,l2,2),
     &     scale,+1,+1,rfbox(6,2,j,l1,l2),rgfbox(6,2,j,l1,l2),
     &     Div)
       enddo
       enddo
       enddo

       wardidtest1 = .true.
       do i=1,6
           if (rgfbox(i,1,1,1,1).gt.1d-1) then 
c            if (wardidtest1) then
c              write(*,*) "MARKER1: gauge test failed", j,i,rgfbox(i,j,1)
c              if (Div.ne.0) then
c              write(*,*) "MARKER4: gauge test failed", j,i,rfbox(i,j,1)
c              endif
c            else
c              write(*,*) "MARKER2: gauge test failed", j,i,rgfbox(i,j,1)
c              if (Div.ne.0) then
c              write(*,*) "MARKER5: gauge test failed", j,i,rfbox(i,j,1)
c              endif
c            endif
             wardidtest1=.false.
           endif
       enddo
c add everything together
       if (wardidtest1) then
cmr             box2count2m = box2count2m+1
cmr             box2count2ps(ps_number) = box2count2ps(ps_number)+1   
ccc         write(*,*) "MARKER3: gauge test ok"
         do l1=1,2
         do l2=1,2
           m1 = 0d0
           m2 = 0d0
           do i=1,6
           do j=1,2
             m1 = m1 + 2*rfbox(i,j,1,l1,l2) ! u,c
     &               + 1*rfbox(i,j,3,l1,l2) ! t
             m2 = m2 + 2*rfbox(i,j,2,l1,l2) ! d,s
     &               + 1*rfbox(i,j,4,l1,l2) ! b
           enddo
           enddo
c and multiply with color factor TR
           mat(1,-1,l1,l2) = (m1 + m2) * TR 
         enddo
         enddo
       endif
      endif

      if (lhiggs) then
c vertex diagrams
ccc strong coupling is added later
ccc triangles with Z/photon vanish, for Higgs only top relevant
ccc Both fermion-line directions identical - factor 2 added later
         do j=3,4 ! no contribution from ~ massless loops -- bottom interference important
         do l1=1,2 
         l2=l1 ! 2-1 and 1-2 vanish
           call gluonverHiggsEvenT(
     &       mf(j),
     &       pbar(0,1),pbar(0,2),minus_qaa(0),
     &       ceps(0,l1,1),ceps(0,l2,2),
     &       scale,1,rfvert(1,3,j,l1,l2),
     &       Div)
         enddo
         rfvert(1,3,j,1,2) = 0d0
         rfvert(1,3,j,2,1) = 0d0
         enddo

c add everything together
         do l1=1,2
         do l2=1,2
           mat(2,-1,l1,l2) = 0d0 
         do j=3,4 
           mat(2,-1,l1,l2) = mat(2,-1,l1,l2)
     &           - rfvert(1,3,j,l1,l2) ! *2*TR ! amplitude (factor 2 for both fermion directions) * colour factor
     &           * mf(j)/vv    ! Yukawa
     &           * w4l(h,1)/dcmplx(qaa(4,1)-xm2(6),xmg(6)) ! decay matrix element * Higgs propagator
         enddo
         enddo
         enddo
       endif

c------------------------------------------------------

c sum the graphs, square them and map them onto uuzz and ddzz

       res(1) = 0
       do l1=1,2
       do l2=1,2
         mm(1,-1) = 0d0
         do i=1,2
           mm(1,-1) = mm(1,-1) + mat(i,-1,l1,l2)
         enddo
         mm(1,-1) = mm(1,-1)*als(1,1)*oneo4pi
         res(1)= res(1) + dreal(mm(1,-1))**2
     &                  + dimag(mm(1,-1))**2
       enddo
       enddo

       res(1) = res(1)*8d0 ! 8 is the color sum factor
      
       uuzz = res(1)

      return
      end

