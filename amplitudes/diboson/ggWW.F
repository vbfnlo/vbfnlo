c****************subroutine GGWW *************************************
c
c	Michael Rauch, <rauch@particle.uni-karlsruhe.de>
c	Last modified: May 2012 (rauch)
c       Including virtual corrections
C
C  GGWW calculates the matrix elements**2 for 
C        g g   ---->   W W   
C
C
C  and crossing related processes. Pauli interference terms for
c  identical fermions are neglected.
c
c          N = 1         Momentum configuration 1 (LO)
c
c*********************************************************************

      subroutine ggww(pbar,fsign, nlo, N, uuww)

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
      double precision  uuww, res(3)
      double precision  p(0:3,2+n_v)


      double precision  eps(0:3,2,2) ! g in kartesian basis
      integer  fsign(2+n_v), sign(2+n_v), nlo, mu, i, j, k, isig, N
      integer  l1, l2


      double complex ceps(0:3,2,2) ! g in complex kartesian basis


      double complex mat(2,-1:1,2,2)
      double complex mm(2,-1:1), zero, im
      double complex m1

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
      double complex rfbox(6,1,4,2,2)
      double precision rgfbox(6,1,4,2,2) 
      double complex rfvert(2,3,4,2,2)


      integer*8 thirtytwo


      parameter(thirtytwo=32)

      Logical wardidtest1

      real*8 minus_qp(0:3),minus_qm(0:3),minus_qww(0:3)


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

c g g -> W W 
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
         minus_qp(mu) = -qp(mu,N)
         minus_qm(mu) = -qm(mu,N)
         minus_qww(mu) = -qww(mu,N)
      enddo

      if (lbox) then
c box loops
ccc first anchor point is always first gluon
ccc permutate second g and the two Z/photons over the other ones
ccc strong coupling is added later
cc g-1-2
       do j=1,3,2
       do l1=1,2
       do l2=1,2
         call fermionboxcouplingT1_GM(
     &     1d0,1d0,clr(cpl(j),3,+1),clr(cpl(j),3,+1),
     &     1d0,1d0,clr(cpl(j),3,-1),clr(cpl(j),3,-1),
     &     mf(j),mf(j),mf(j+1),mf(j),
     &     pbar(0,1),pbar(0,2),minus_qp(0),minus_qm(0),
     &     ceps(0,l1,1),ceps(0,l2,2),wp(1,N),wm(1,N),
     &     scale,+1,+1,rfbox(1,1,j,l1,l2),rgfbox(1,1,j,l1,l2),
     &     Div)
       enddo
       enddo
       enddo
cc g-2-1
       do j=1,3,2
       do l1=1,2
       do l2=1,2
         call fermionboxcouplingT1_GM(
     &     1d0,1d0,clr(cpl(j),3,+1),clr(cpl(j),3,+1),
     &     1d0,1d0,clr(cpl(j),3,-1),clr(cpl(j),3,-1),
     &     mf(j+1),mf(j+1),mf(j),mf(j+1),
     &     pbar(0,1),pbar(0,2),minus_qm(0),minus_qp(0),
     &     ceps(0,l1,1),ceps(0,l2,2),wm(1,N),wp(1,N),
     &     scale,+1,+1,rfbox(2,1,j,l1,l2),rgfbox(2,1,j,l1,l2),
     &     Div)
       enddo
       enddo
       enddo
cc 1-g-2
       do j=1,3,2
       do l1=1,2
       do l2=1,2
         call fermionboxcouplingT1_GM(
     &     1d0,clr(cpl(j),3,+1),1d0,clr(cpl(j),3,+1),
     &     1d0,clr(cpl(j),3,-1),1d0,clr(cpl(j),3,-1),
     &     mf(j),mf(j+1),mf(j+1),mf(j),
     &     pbar(0,1),minus_qp(0),pbar(0,2),minus_qm(0),
     &     ceps(0,l1,1),wp(1,N),ceps(0,l2,2),wm(1,N),
     &     scale,+1,+1,rfbox(3,1,j,l1,l2),rgfbox(3,1,j,l1,l2),
     &     Div)
       enddo
       enddo
       enddo
cc 2-g-1
       do j=1,3,2
       do l1=1,2
       do l2=1,2
         call fermionboxcouplingT1_GM(
     &     1d0,clr(cpl(j),3,+1),1d0,clr(cpl(j),3,+1),
     &     1d0,clr(cpl(j),3,-1),1d0,clr(cpl(j),3,-1),
     &     mf(j+1),mf(j),mf(j),mf(j+1),
     &     pbar(0,1),minus_qm(0),pbar(0,2),minus_qp(0),
     &     ceps(0,l1,1),wm(1,N),ceps(0,l2,2),wp(1,N),
     &     scale,+1,+1,rfbox(4,1,j,l1,l2),rgfbox(4,1,j,l1,l2),
     &     Div)
       enddo
       enddo
       enddo
cc 1-2-g
       do j=1,3,2
       do l1=1,2
       do l2=1,2
         call fermionboxcouplingT1_GM(
     &     1d0,clr(cpl(j),3,+1),clr(cpl(j),3,+1),1d0,
     &     1d0,clr(cpl(j),3,-1),clr(cpl(j),3,-1),1d0,
     &     mf(j),mf(j+1),mf(j),mf(j),
     &     pbar(0,1),minus_qp(0),minus_qm(0),pbar(0,2),
     &     ceps(0,l1,1),wp(1,N),wm(1,N),ceps(0,l2,2),
     &     scale,+1,+1,rfbox(5,1,j,l1,l2),rgfbox(5,1,j,l1,l2),
     &     Div)
       enddo
       enddo
       enddo
cc 2-1-g
       do j=1,3,2
       do l1=1,2
       do l2=1,2
         call fermionboxcouplingT1_GM(
     &     1d0,clr(cpl(j),3,+1),clr(cpl(j),3,+1),1d0,
     &     1d0,clr(cpl(j),3,-1),clr(cpl(j),3,-1),1d0,
     &     mf(j+1),mf(j),mf(j+1),mf(j+1),
     &     pbar(0,1),minus_qm(0),minus_qp(0),pbar(0,2),
     &     ceps(0,l1,1),wm(1,N),wp(1,N),ceps(0,l2,2),
     &     scale,+1,+1,rfbox(6,1,j,l1,l2),rgfbox(6,1,j,l1,l2),
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
           do i=1,6
           do j=1,1
             m1 = m1 + 2*rfbox(i,j,1,l1,l2) ! u,c
     &               + 1*rfbox(i,j,3,l1,l2) ! t
           enddo
           enddo
c and multiply with color factor TR
           mat(1,-1,l1,l2) = (m1) * TR 
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
     &       pbar(0,1),pbar(0,2),minus_qww(0),
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
     &           * w4l(h,1)/dcmplx(qww(4,1)-xm2(6),xmg(6)) ! decay matrix element * Higgs propagator
         enddo
         enddo
         enddo
       endif
c	write(*,*) "PPP", pbar
c       do i=1,6
c       do j=1,3,2
c       do l1=1,2
c       do l2=1,2
c	write(*,*) "GGG", i,j,l1,l2,rfbox(i,1,j,l1,l2) 
c       enddo                
c       enddo                
c       enddo                
c       enddo
c	read(*,*) 

c------------------------------------------------------

c sum the graphs, square them and map them onto uuww 

       res(1) = 0
       do l1=1,2
       do l2=1,2
         mm(1,-1) = 0d0
         do i=1,2
           mm(1,-1) = mm(1,-1) + mat(i,-1,l1,l2)
         enddo
         mm(1,-1) = mm(1,-1)*als(1,1)*oneo4pi
clt         mm(1,-1) = mm(1,-1) + fckamp(l1,l2)
         res(1)= res(1) + dreal(mm(1,-1))**2
     &                  + dimag(mm(1,-1))**2
       enddo
       enddo

       res(1) = res(1)*8d0 ! 8 is the color sum factor
clt       res(1) = res(1)*2d0 ! 8 is the color sum factor
      
       uuww = res(1)

      return
      end

