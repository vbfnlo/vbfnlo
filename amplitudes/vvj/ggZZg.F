c****************subroutine GGZZG *************************************
c
c	Michael Rauch, <rauch@particle.physik.uni-karlsruhe.de>
c	Last modified: June 2012
C
C  GGZZG calculates the matrix elements**2 for 
C        g g   ---->   (Z/A) (Z/A) g, 
C
C  and crossing related processes. Pauli interference terms for
c  identical fermions are neglected.
c
c  Gluon polarization vectors are given in the kartesian basis
c  (comparison with madgraph amplitudes tampr requires to express these in
c   kart. basis also, even though they are given in the helicity basis
c
c*********************************************************************

      subroutine ggZZg(pbar,sign, qbar, gsign, bos, nlo, N, ggZZj)

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/scales.inc"
c#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/koppln.inc"
#include "VBFNLO/utilities/process.inc"
#include "VBFNLO/loops/Utilities/GaugeTest.inc"

c electroweak couplings are taken from KOPPLN
      double precision  clr, xm2, xmg, b, v, a
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     1                  V(4,5),A(4,5)

c variables for the main part of the program
      double precision fpi
      parameter (fpi=4d0*pi)

      double precision  pbar(0:3,2+n_v), qbar(0:4), q(0:3)
      double precision  pglu(0:4,3), ptglu(0:4,3)
      double precision  minus_qp(0:3,1),minus_qm(0:3,1),minus_qaa(0:3,1)
      double precision  ggZZj, res(2)
      double precision  p(0:3,2+n_v)

      integer  sign(2+n_v), gsign, mu, i, j, k, bos, nlo, N
      integer  l(3), ll(3), lh(3), l1, l2, l3, k1, k2
      integer perm(6,3), signperm(6), redo, helmult, signglu(3)

      logical lvert, lbox, lpent,lbox2, lhiggs
      logical lhaszvv
      logical wardidtest1

      double precision eps(0:3,2,3) ! g in kartesian basis
      double complex ceps(0:5,2,3), cteps(0:5,2,2,3)
      double complex mat(6,2,2,2,2)
      double complex mm(2), zero
      double complex m1, m2, m3, m4

      double complex im
      double complex dotcc, s1c, dotrc
      double precision dotrr, psumsq
      external dotcc, s1c, dotrr, dotrc, psumsq
      parameter (im = (0d0,1d0))
      parameter (zero = (0d0,0d0) )

c Variables for fermion loops

      integer Div, cpl(4)
      double precision mf(4), vv, scale, musq
      double complex rfpenE(6,4,4,4,2,2,2)
      double complex rfpenO(6,4,4,4,2,2,2)
      double precision rgfpen(6,4,4,4,2,2,2) 
      double complex rfbox(6,4,4,3,2,2,2)
      double precision rgfbox(6,4,4,3,2,2,2) 
      double complex rfvert(6,4,4,4,2,2,2)
      double precision rgfvert(6,4,4,4,2,2,2) 
      double precision TR, oneo4pi
      parameter (TR=1d0/2d0)
      parameter (oneo4pi=1d0/4d0/pi)

      double precision RandomNumber
      external RandomNumber

c variables for ZZ,ZA,AA 
      integer typ, typl
      integer cpp(4), cpm(4)
      double complex cepsp(6,4), cepsm(6,4)

c  helicity selection
      INTEGER h
      COMMON /hcount / h

c Gram determinant
      logical singular
      common/DetCountML/singular

      logical ldebugm, linit, lhelrand
      data linit /.true./
      save linit
      parameter(ldebugm = .false.)
      parameter(lhelrand = .true.)


c---------------------------- BEGIN CODE -----------------------

c initialize & precompute stuff needed below:

      musq=mursq(1,1) ! renormalization scale

      singular = .false.

c whether to include diagrams with boxes or Higgs bosons
      lpent  = (floops.eq.1).or.(floops.eq.3) 
      lbox   = (floops.eq.1).or.(floops.eq.3)
      lbox2  = (floops.eq.1).or.(floops.eq.3)
      lvert  = (floops.eq.1).or.(floops.eq.3)
      lhiggs = (floops.eq.2).or.(floops.eq.3)

c map polarisation vectors and indices onto the arrays, so for the loops we can use do-loops
c makes debugging way easier
c typ: 1:AA, 2: ZA, 3: AZ, 4: ZZ (so can just stop at right index for AA and ZA)
      if (linit) then
        cpp(1) = 1
        cpp(2) = 2
        cpp(3) = 1
        cpp(4) = 2
        cpm(1) = 1
        cpm(2) = 1
        cpm(3) = 2
        cpm(4) = 2
        if (bos.eq.22) then
          typ = 4
          lhaszvv = .true.
        else if (bos.eq.21) then
          typ = 2
          lhaszvv = .true.
        else if (bos.eq.11) then
          typ = 1
          lhaszvv = .false.
        else
          write(*,*) "Error in ggZZg: wrong bos ",bos
          stop
        endif
        if (narrowwidth.eq.0) then
          typl = 1
        else
          typl = typ
        endif
        linit = .false.
      endif
      do mu=1,6
        cepsp(mu,1) = atau1(mu,N)
        cepsp(mu,2) = wp(mu,N)
        cepsp(mu,3) = atau1(mu,N)
        cepsp(mu,4) = wp(mu,N)
        cepsm(mu,1) = atau2(mu,N)
        cepsm(mu,2) = atau2(mu,N)
        cepsm(mu,3) = wm(mu,N)
        cepsm(mu,4) = wm(mu,N)
      enddo
      if (bos.eq.11) then
        do mu=0,3
          qp(mu,N) = qa1(mu,N)
          qm(mu,N) = qa2(mu,N)
        enddo
      endif

      do k = 1,6
      do j = 1,2
        do l1 = 1,2
        do l2 = 1,2
        do l3 = 1,2
              mat(k,j,l1,l2,l3) = (0.0d0,0.0d0)
        enddo
        enddo
        enddo
      enddo
      enddo

c define the internal momenta
      do mu = 0,3
         do i = 1,2+n_v
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo
         q(mu) = qbar(mu)*gsign
c pglu vector is all particles outgoing 
c -- initials are incoming (sign modifies which ones are really incoming)
         do i = 1,2 
           pglu(mu,i) = p(mu,i)
           signglu(i) = sign(i)
         enddo
         pglu(mu,3) = -q(mu)
         signglu(3) = gsign
         ptglu(mu,1) = pglu(mu,2)+pglu(mu,3)
         ptglu(mu,2) = pglu(mu,1)+pglu(mu,3)
         ptglu(mu,3) = pglu(mu,1)+pglu(mu,2)
         minus_qp(mu,1)  = -qp(mu,1)
         minus_qm(mu,1)  = -(pglu(mu,1)+pglu(mu,2)+ pglu(mu,3)+ minus_qp(mu,1)) !-qm(mu,1)
         minus_qaa(mu,1) = -(pglu(mu,1)+pglu(mu,2)+ pglu(mu,3))    !-qaa(mu,1)
      enddo
      do i = 1,3 
         pglu(4,i) = 0d0
         ptglu(4,i) = dotrr(ptglu(0,i),ptglu(0,i))
      enddo
c  Get the gluon polarization vector and the gluon emission spinors

      do k = 1,2        ! 2 gluon polarizations
      do i = 1,3   
        call helvec(pglu(0,i),+1,2*k-3,ceps(0,k,i))  ! get gluon pol.vectors
      enddo
      enddo
c combine two gluons 1=23, 2=31, 3=12
      do l1 = 1,2        
      do l2 = 1,2        
      do mu = 0,3    
        cteps(mu,l1,l2,1) = (
     &        2d0*dotrc(pglu(0,2),ceps(0,l2,3))*ceps(mu,l1,2)
     &       -2d0*dotrc(pglu(0,3),ceps(0,l1,2))*ceps(mu,l2,3)
     &       +dotcc(ceps(0,l1,2),ceps(0,l2,3))*(pglu(mu,3)-pglu(mu,2)))
     &       /(ptglu(4,1))
        cteps(mu,l1,l2,2) = (
     &        2d0*dotrc(pglu(0,3),ceps(0,l2,1))*ceps(mu,l1,3)
     &       -2d0*dotrc(pglu(0,1),ceps(0,l1,3))*ceps(mu,l2,1)
     &       +dotcc(ceps(0,l1,3),ceps(0,l2,1))*(pglu(mu,1)-pglu(mu,3)))
     &       /(ptglu(4,2))
        cteps(mu,l1,l2,3) = (
     &        2d0*dotrc(pglu(0,1),ceps(0,l2,2))*ceps(mu,l1,1)
     &       -2d0*dotrc(pglu(0,2),ceps(0,l1,1))*ceps(mu,l2,2)
     &       +dotcc(ceps(0,l1,1),ceps(0,l2,2))*(pglu(mu,2)-pglu(mu,1)))
     &       /(ptglu(4,3))
      enddo
      enddo
      enddo
c define all permutation of 1,2,3
      j=1
      do l1=1,3
      do l2=1,3
        if (l2.eq.l1) cycle
      do l3=1,3
        if ((l3.eq.l1).or.(l3.eq.l2)) cycle
        perm(j,1)=l1
        perm(j,2)=l2
        perm(j,3)=l3
        if ( (mod(l1,3)+1.eq.l2) .and. (mod(l2,3)+1.eq.l3) ) then
! cyclic permutation
          signperm(j) = +1
        else
          signperm(j) = -1
        endif
        j=j+1
      enddo
      enddo
      enddo

cc assertion check - must be 6 permutations (and +1 at the end)
c      if (j-1.ne.6) then
c        write(*,*) "Error in ggWWg: number of permutations wrong: ",j-1
c        do k=1,j-1
c          write(*,*) "Perm.",k,": ",perm(k,1),perm(k,2),perm(k,3), signperm(k)
c        enddo
c        stop
c      endif

c---- QCD gaugecheck --- eps_mu -> k_mu --- k_mu M^mu = 0 ------------

c         do mu = 0,3                   !QCD gaugecheck
c            eps(mu,l) = qbar(mu)
c         enddo

c--- End of QCD gaugecheck -------------------------------------------

c switch if random helicity summation for gluons
       helmult = 1
       if (lhelrand) then
         helmult = 2**3
         if (RandomNumber().lt.0.5) then
           ll(1) = 1
         else
           ll(1) = 2
         endif
         if (RandomNumber().lt.0.5) then
           ll(2) = 1
         else
           ll(2) = 2
         endif
         if (RandomNumber().lt.0.5) then
           ll(3) = 1
         else
           ll(3) = 2
         endif
ccc fixed values for debugging
c           ll(1) = 1
c           ll(2) = 1
c           ll(3) = 1
ccc---
         lh(1) = ll(1)
         lh(2) = ll(2)
         lh(3) = ll(3)
       else
         ll(1)=1
         ll(2)=1
         ll(3)=1
         lh(1)=2
         lh(2)=2
         lh(3)=2
       endif

c colour structure:
c All diagrams can be mapped onto two different structures:
c (gluon colours g1^a g2^b g3^c)
c 1: tr(t^a t^b t^c) 
c 2: tr(t^c t^b t^a) 

c--- start calculation ---
       scale=musq  !mursq(1,1)!-2d0*dotrr(p(0,1),p(0,2))
       Div=0

       mf(1) = 0d0 ! light quarks
       mf(2) = 0d0 ! light quarks
       mf(3) = xmt ! top
       mf(4) = bmass_h ! bottom
       vv = 2.d0*sqrt(xm2(3)*(1-xm2(3)/xm2(2)))/(-clr(2,1,-1)) ! 2*MW*SW/e

       cpl(1) = 3 ! light up quarks
       cpl(2) = 4 ! light down quarks
       cpl(3) = 3 ! top
       cpl(4) = 4 ! bottom

      if (lhiggs) then
c box diagrams
ccc strong coupling is added later
         do j=3,4 ! no contribution from ~ massless loops -- bottom interference important
         do k=1,6
         if(signperm(k).eq.1) then ! signperm(k) == -1 has only relative minus sign
         do l1=ll(1),lh(1) 
         l(1) = l1
         do l2=ll(2),lh(2) 
         l(2) = l2
         do l3=ll(3),lh(3) 
         l(3) = l3
           call gluonboxHiggsEven(
     &       mf(j),
     &       pglu(0,perm(k,1)),pglu(0,perm(k,2)),pglu(0,perm(k,3)),minus_qaa(0,1),
     &       ceps(0,l(perm(k,1)),perm(k,1)),ceps(0,l(perm(k,2)),perm(k,2)),ceps(0,l(perm(k,3)),perm(k,3)),
     &       rfbox(k,1,j,1,l(1),l(2),l(3)))
         enddo
         enddo
         enddo
         endif
         enddo
         enddo
c add everything together
         do l1=ll(1),lh(1)
         do l2=ll(2),lh(2)
         do l3=ll(3),lh(3)
           m1 = 0d0 
           m2 = 0d0 
           do j=3,4 
           do k=1,6 
c take out factor 3*m from Higgs-box routine
            if(signperm(k).eq.1) then
             m1 = m1 
     &           - rfbox(k,1,j,1,l1,l2,l3) 
     &           * 1d0/3d0/vv    ! Yukawa  (mf cancels with factor!)
     &           * w4l(h,1)/dcmplx(qaa(4,1)-xm2(6),xmg(6)) ! decay matrix element * Higgs propagator
            endif
           enddo
           enddo
c colour assignment depending on permutation
           mat(5,1,l1,l2,l3) =  m1 
           mat(5,2,l1,l2,l3) = -m1 
         enddo
         enddo
         enddo

c vertex diagrams
         do j=3,4 ! no contribution from ~ massless loops -- bottom interference important
         do k=1,3
           k1 = mod(k  ,3)+1
           k2 = mod(k+1,3)+1
         do l1=ll(k ),lh(k )
         do l2=ll(k1),lh(k1)
         do l3=ll(k2),lh(k2)
           l(k ) = l1
           l(k1) = l2
           l(k2) = l3
c g-(gg)-H
           call gluonverHiggsEvenT(
     &       mf(j),
     &       pglu(0,k),ptglu(0,k),minus_qaa(0,1),
     &       ceps(0,l1,k),cteps(0,l2,l3,k),
     &       scale,1,rfvert(1,1,j,k,l(1),l(2),l(3)),
     &       Div)
cc (gg)-g-H
c           call gluonverHiggsEvenT(
c     &       mf(j),
c     &       ptglu(0,k),pglu(0,k),minus_qaa(0,1),
c     &       cteps(0,l2,l3,k),ceps(0,l1,k),
c     &       scale,1,rfvert(2,1,j,k,l(1),l(2),l(3)),
c     &       Div)
         enddo
         enddo
         enddo
         enddo
         enddo

c add everything together
         do l1=ll(1),lh(1)
         do l2=ll(2),lh(2)
         do l3=ll(3),lh(3)
           m1 = 0d0
         do i=1,1 !2
         do j=3,4 
         do k=1,3 
           m1 = m1
     &           - rfvert(i,1,j,k,l1,l2,l3) * 2d0 ! both diagrams give the same value
     &           * mf(j)/vv    ! Yukawa
     &           * w4l(h,1)/dcmplx(qaa(4,1)-xm2(6),xmg(6)) ! decay matrix element * Higgs propagator
         enddo
         enddo
         enddo
           mat(6,1,l1,l2,l3) =  m1 
           mat(6,2,l1,l2,l3) = -m1 
         enddo
         enddo
         enddo
      endif

      if (lvert .and. lhaszvv) then
c vertex diagrams gg->g*-(triangle)>gZ,Z->VV
       do j=1,4
       do k=1,3
         k1 = mod(k  ,3)+1
         k2 = mod(k+1,3)+1
       do l1=ll(k ),lh(k )
       do l2=ll(k1),lh(k1)
       do l3=ll(k2),lh(k2)
         l(k ) = l1
         l(k1) = l2
         l(k2) = l3
c g-(gg)-Z
         call fermionvertcouplingT(
     &     1d0,1d0,clr(cpl(j),2,+1),
     &     1d0,1d0,clr(cpl(j),2,-1),
     &     mf(j),
     &     pglu(0,k),ptglu(0,k),minus_qaa(0,1),
     &     ceps(0,l1,k),cteps(0,l2,l3,k),wz(1,N),
     &     scale,+1,+1,rfvert(1,1,j,k,l(1),l(2),l(3)),
     &                rgfvert(1,1,j,k,l(1),l(2),l(3)),Div)
         call fermionvertcouplingT(
     &     1d0,1d0,clr(cpl(j),1,+1),
     &     1d0,1d0,clr(cpl(j),1,-1),
     &     mf(j),
     &     pglu(0,k),ptglu(0,k),minus_qaa(0,1),
     &     ceps(0,l1,k),cteps(0,l2,l3,k),wa(1,N),
     &     scale,+1,+1,rfvert(2,1,j,k,l(1),l(2),l(3)),
     &                rgfvert(2,1,j,k,l(1),l(2),l(3)),Div)
c (gg)-g-Z
         call fermionvertcouplingT(
     &     1d0,1d0,clr(cpl(j),2,+1),
     &     1d0,1d0,clr(cpl(j),2,-1),
     &     mf(j),
     &     ptglu(0,k),pglu(0,k),minus_qaa(0,1),
     &     cteps(0,l2,l3,k),ceps(0,l1,k),wz(1,N),
     &     scale,+1,+1,rfvert(3,1,j,k,l(1),l(2),l(3)),
     &                rgfvert(3,1,j,k,l(1),l(2),l(3)),Div)
         call fermionvertcouplingT(
     &     1d0,1d0,clr(cpl(j),1,+1),
     &     1d0,1d0,clr(cpl(j),1,-1),
     &     mf(j),
     &     ptglu(0,k),pglu(0,k),minus_qaa(0,1),
     &     cteps(0,l2,l3,k),ceps(0,l1,k),wa(1,N),
     &     scale,+1,+1,rfvert(4,1,j,k,l(1),l(2),l(3)),
     &                rgfvert(4,1,j,k,l(1),l(2),l(3)),Div)
       enddo
       enddo
       enddo
       enddo
       enddo

       wardidtest1 = .true.
       do i=1,4
         do j=1,3
           if (rgfvert(i,1,1,j,ll(1),ll(2),ll(3)).gt.GaugeAccuracyBoxL) then 
             wardidtest1=.false.
c             write(*,*) "Gauge Vert:", i,j,rgfvert(i,1,1,j,ll(1),ll(2),ll(3)), singular
             ggZZj = 0d0
             return
           endif
         enddo
       enddo
c add everything together
       do l1=ll(1),lh(1)
       do l2=ll(2),lh(2)
       do l3=ll(3),lh(3)
         m1 = 0d0
       do j=1,4
       do k=1,3
           m1 = m1 + 2*rfvert(j,1,1,k,l1,l2,l3) ! u,c
     &             + 1*rfvert(j,1,3,k,l1,l2,l3) ! t
     &             + 2*rfvert(j,1,2,k,l1,l2,l3) ! d,s
     &             + 1*rfvert(j,1,4,k,l1,l2,l3) ! b
       enddo
       enddo
c and multiply with color factor 
         mat(4,1,l1,l2,l3) =  m1 
         mat(4,2,l1,l2,l3) = -m1 
       enddo
       enddo
       enddo
      endif

      if (lbox) then
c box gauge loops
       do j=1,4 !Change for checking
       do k=1,3
         k1 = mod(k  ,3)+1
         k2 = mod(k+1,3)+1
       do i=typ,typl,-1
       do l1=ll(k ),lh(k )
       do l2=ll(k1),lh(k1)
       do l3=ll(k2),lh(k2)
         l(k ) = l1
         l(k1) = l2
         l(k2) = l3
cc g-(gg)-p-m
         call fermionboxcouplingT(
     &     1d0,1d0,clr(cpl(j),cpp(i),+1),clr(cpl(j),cpm(i),+1),
     &     1d0,1d0,clr(cpl(j),cpp(i),-1),clr(cpl(j),cpm(i),-1),
     &     mf(j),
     &     pglu(0,k),ptglu(0,k),minus_qp(0,1),minus_qm(0,1),
     &     ceps(0,l1,k),cteps(0,l2,l3,k),cepsp(1,i),cepsm(1,i),
     &     scale,+1,+1,rfbox(1,i,j,k,l(1),l(2),l(3)),
     &                rgfbox(1,i,j,k,l(1),l(2),l(3)),Div)
cc (gg)-g-p-m
         call fermionboxcouplingT(
     &     1d0,1d0,clr(cpl(j),cpp(i),+1),clr(cpl(j),cpm(i),+1),
     &     1d0,1d0,clr(cpl(j),cpp(i),-1),clr(cpl(j),cpm(i),-1),
     &     mf(j),
     &     ptglu(0,k),pglu(0,k),minus_qp(0,1),minus_qm(0,1),
     &     cteps(0,l2,l3,k),ceps(0,l1,k),cepsp(1,i),cepsm(1,i),
     &     scale,+1,+1,rfbox(2,i,j,k,l(1),l(2),l(3)),
     &                rgfbox(2,i,j,k,l(1),l(2),l(3)),Div)
cc g-(gg)-m-p
         call fermionboxcouplingT(
     &     1d0,1d0,clr(cpl(j),cpm(i),+1),clr(cpl(j),cpp(i),+1),
     &     1d0,1d0,clr(cpl(j),cpm(i),-1),clr(cpl(j),cpp(i),-1),
     &     mf(j),
     &     pglu(0,k),ptglu(0,k),minus_qm(0,1),minus_qp(0,1),
     &     ceps(0,l1,k),cteps(0,l2,l3,k),cepsm(1,i),cepsp(1,i),
     &     scale,+1,+1,rfbox(3,i,j,k,l(1),l(2),l(3)),
     &                rgfbox(3,i,j,k,l(1),l(2),l(3)),Div)
cc (gg)-g-m-p
         call fermionboxcouplingT(
     &     1d0,1d0,clr(cpl(j),cpm(i),+1),clr(cpl(j),cpp(i),+1),
     &     1d0,1d0,clr(cpl(j),cpm(i),-1),clr(cpl(j),cpp(i),-1),
     &     mf(j),
     &     ptglu(0,k),pglu(0,k),minus_qm(0,1),minus_qp(0,1),
     &     cteps(0,l2,l3,k),ceps(0,l1,k),cepsm(1,i),cepsp(1,i),
     &     scale,+1,+1,rfbox(4,i,j,k,l(1),l(2),l(3)),
     &                rgfbox(4,i,j,k,l(1),l(2),l(3)),Div)
cc g-p-(gg)-m
         call fermionboxcouplingT(
     &     1d0,clr(cpl(j),cpp(i),+1),1d0,clr(cpl(j),cpm(i),+1),
     &     1d0,clr(cpl(j),cpp(i),-1),1d0,clr(cpl(j),cpm(i),-1),
     &     mf(j),
     &     pglu(0,k),minus_qp(0,1),ptglu(0,k),minus_qm(0,1),
     &     ceps(0,l1,k),cepsp(1,i),cteps(0,l2,l3,k),cepsm(1,i),
     &     scale,+1,+1,rfbox(5,i,j,k,l(1),l(2),l(3)),
     &                rgfbox(5,i,j,k,l(1),l(2),l(3)),Div)
cc (gg)-p-g-m
         call fermionboxcouplingT(
     &     1d0,clr(cpl(j),cpp(i),+1),1d0,clr(cpl(j),cpm(i),+1),
     &     1d0,clr(cpl(j),cpp(i),-1),1d0,clr(cpl(j),cpm(i),-1),
     &     mf(j),
     &     ptglu(0,k),minus_qp(0,1),pglu(0,k),minus_qm(0,1),
     &     cteps(0,l2,l3,k),cepsp(1,i),ceps(0,l1,k),cepsm(1,i),
     &     scale,+1,+1,rfbox(6,i,j,k,l(1),l(2),l(3)),
     &                rgfbox(6,i,j,k,l(1),l(2),l(3)),Div)
       enddo
       enddo
       enddo
       enddo
       enddo
       enddo

       wardidtest1 = .true.
      do k=1,4
          do j=1,3
             do i=1,6
           if (rgfbox(i,typ,k,j,ll(1),ll(2),ll(3)).gt.GaugeAccuracyBoxL) then 
             wardidtest1=.false.
c             write(*,*) "Gauge Box1:",k,i,j,rgfbox(i,typ,k,j,ll(1),ll(2),ll(3)), singular
              ggZZj = 0d0
            return
           endif
             enddo
         enddo
        enddo
     

c add everything together
       do l1=ll(1),lh(1)
       do l2=ll(2),lh(2)
       do l3=ll(3),lh(3)
         m1 = 0d0
       do j=1,6
       do k=1,3
       do i=typl,typ
           m1 = m1 + 2*rfbox(j,i,1,k,l1,l2,l3) ! u,c
     &             + 1*rfbox(j,i,3,k,l1,l2,l3) ! t
     &             + 2*rfbox(j,i,2,k,l1,l2,l3) ! d,s
     &             + 1*rfbox(j,i,4,k,l1,l2,l3) ! b
       enddo
       enddo
       enddo
c gluon triangle maps onto both with relative sign
         mat(2,1,l1,l2,l3) =  m1 
         mat(2,2,l1,l2,l3) = -m1 
       enddo
       enddo
       enddo
      endif ! lbox

      if (lbox2 .and. lhaszvv) then
       do j=1,4
       do k=1,6
       do l1=ll(1),lh(1)
         l(1) = l1
       do l2=ll(2),lh(2)
         l(2) = l2
       do l3=ll(3),lh(3)
         l(3) = l3
         call fermionboxcouplingT(
     &     1d0,1d0,1d0,clr(cpl(j),2,+1),
     &     1d0,1d0,1d0,clr(cpl(j),2,-1),
     &     mf(j),
     &     pglu(0,perm(k,1)),pglu(0,perm(k,2)),pglu(0,perm(k,3)),minus_qaa(0,1),
     &     ceps(0,l(perm(k,1)),perm(k,1)),ceps(0,l(perm(k,2)),perm(k,2)),ceps(0,l(perm(k,3)),perm(k,3)),wz(1,1),
     &     scale,+1,+1,rfbox(k,2,j,1,l(1),l(2),l(3)),rgfbox(k,2,j,1,l(1),l(2),l(3)),
     &     Div)
         call fermionboxcouplingT(
     &     1d0,1d0,1d0,clr(cpl(j),1,+1),
     &     1d0,1d0,1d0,clr(cpl(j),1,-1),
     &     mf(j),
     &     pglu(0,perm(k,1)),pglu(0,perm(k,2)),pglu(0,perm(k,3)),minus_qaa(0,1),
     &     ceps(0,l(perm(k,1)),perm(k,1)),ceps(0,l(perm(k,2)),perm(k,2)),ceps(0,l(perm(k,3)),perm(k,3)),wa(1,1),
     &     scale,+1,+1,rfbox(k,2,j,2,l(1),l(2),l(3)),rgfbox(k,2,j,2,l(1),l(2),l(3)),
     &     Div)
       enddo
       enddo
       enddo
       enddo
       enddo

       wardidtest1 = .true.
      do k=1,4
        do i=1,6
         do j=1,2
             if (rgfbox(i,2,k,j,ll(1),ll(2),ll(3)).gt.GaugeAccuracyBoxL) then 
             wardidtest1=.false.
c             write(*,*) "Gauge Box2:", k,i,j,rgfbox(i,2,k,j,ll(1),ll(2),ll(3)), singular
             ggZZj = 0d0
             return
           endif
         enddo
         enddo
        enddo

c add everything together
       do l1=ll(1),lh(1)
       do l2=ll(2),lh(2)
       do l3=ll(3),lh(3)
         m1 = 0d0
         m2 = 0d0
         do k=1,6
         do i=1,2
          if(signperm(k).eq.1) then
           m1 = m1 + 2*rfbox(k,2,1,i,l1,l2,l3) ! u,c
     &             + 1*rfbox(k,2,3,i,l1,l2,l3) ! t
     &             + 2*rfbox(k,2,2,i,l1,l2,l3) ! d,s
     &             + 1*rfbox(k,2,4,i,l1,l2,l3) ! b
          else
           m2 = m2 + 2*rfbox(k,2,1,i,l1,l2,l3) ! u,c
     &             + 1*rfbox(k,2,3,i,l1,l2,l3) ! t
     &             + 2*rfbox(k,2,2,i,l1,l2,l3) ! d,s
     &             + 1*rfbox(k,2,4,i,l1,l2,l3) ! b
          endif
         enddo
         enddo
c same colour structure as defining diagram
         mat(3,1,l1,l2,l3) = m1
         mat(3,2,l1,l2,l3) = m2
       enddo
       enddo
       enddo

      endif !lbox2

      if (lpent) then
c non-Higgs diagrams
c pentagon loops
ccc p/m is Z_e/Z_mu
cc g-g-g-p-m
       redo=1
       do j=1,4
       do k=1,6
       do i=typ,typl,-1 ! loop reversed so gauge test done on last (no gamma* in narrowwidth)
       do l1=ll(1),lh(1)
         l(1) = l1
       do l2=ll(2),lh(2)
         l(2) = l2
       do l3=ll(3),lh(3)
         l(3) = l3
         call FermionPenCouplingT1_GM345(
     &     1d0,clr(cpl(j),cpp(i),+1),clr(cpl(j),cpm(i),+1),
     &     1d0,clr(cpl(j),cpp(i),-1),clr(cpl(j),cpm(i),-1),
     &     mf(j),mf(j),mf(j),
     &     pglu(0,perm(k,1)),pglu(0,perm(k,2)),pglu(0,perm(k,3)),minus_qp(0,1),minus_qm(0,1),
     &     ceps(0,l(perm(k,1)),perm(k,1)),ceps(0,l(perm(k,2)),perm(k,2)),ceps(0,l(perm(k,3)),perm(k,3)),cepsp(1,i),cepsm(1,i),
     &     scale,redo,redo,rfpenE(k,i,j,1,l(1),l(2),l(3)),rfpenO(k,i,j,1,l(1),l(2),l(3)),rgfpen(k,i,j,1,l(1),l(2),l(3)),
     &     Div)
           redo=-1
       enddo
       enddo
       enddo
       enddo
       redo=1
       enddo
       enddo
ccc g-g-g-m-p
c       redo=1
c       do j=1,4
c       do k=1,6
c       do i=typ,typl,-1
c       do l1=ll(1),lh(1)
c         l(1) = l1
c       do l2=ll(2),lh(2)
c         l(2) = l2
c       do l3=ll(3),lh(3)
c         l(3) = l3
c    
c          rfpenE(6,typ,j,2,ll(1),ll(2),ll(3))=-rfpenE(1,typ,j,1,ll(1),ll(2),ll(3))
c          rfpenE(4,typ,j,2,ll(1),ll(2),ll(3))=-rfpenE(2,typ,j,1,ll(1),ll(2),ll(3))
c          rfpenE(5,typ,j,2,ll(1),ll(2),ll(3))=-rfpenE(3,typ,j,1,ll(1),ll(2),ll(3))
c          rfpenE(2,typ,j,2,ll(1),ll(2),ll(3))=-rfpenE(4,typ,j,1,ll(1),ll(2),ll(3))
c          rfpenE(3,typ,j,2,ll(1),ll(2),ll(3))=-rfpenE(5,typ,j,1,ll(1),ll(2),ll(3))
c          rfpenE(1,typ,j,2,ll(1),ll(2),ll(3))=-rfpenE(6,typ,j,1,ll(1),ll(2),ll(3))
c          
c          rfpenO(6,typ,j,2,ll(1),ll(2),ll(3))=rfpenO(1,typ,j,1,ll(1),ll(2),ll(3))
c          rfpenO(4,typ,j,2,ll(1),ll(2),ll(3))=rfpenO(2,typ,j,1,ll(1),ll(2),ll(3))
c          rfpenO(5,typ,j,2,ll(1),ll(2),ll(3))=rfpenO(3,typ,j,1,ll(1),ll(2),ll(3))
c          rfpenO(2,typ,j,2,ll(1),ll(2),ll(3))=rfpenO(4,typ,j,1,ll(1),ll(2),ll(3))
c          rfpenO(3,typ,j,2,ll(1),ll(2),ll(3))=rfpenO(5,typ,j,1,ll(1),ll(2),ll(3))
c          rfpenO(1,typ,j,2,ll(1),ll(2),ll(3))=rfpenO(6,typ,j,1,ll(1),ll(2),ll(3))

cfc         call FermionPenCouplingT1_GM345(
cfc     &     1d0,clr(cpl(j),cpm(i),+1),clr(cpl(j),cpp(i),+1),
cfc     &     1d0,clr(cpl(j),cpm(i),-1),clr(cpl(j),cpp(i),-1),
cfc     &     mf(j),mf(j),mf(j),
cfc     &     pglu(0,perm(k,1)),pglu(0,perm(k,2)),pglu(0,perm(k,3)),minus_qm(0,1),minus_qp(0,1),
cfc     &     ceps(0,l(perm(k,1)),perm(k,1)),ceps(0,l(perm(k,2)),perm(k,2)),ceps(0,l(perm(k,3)),perm(k,3)),cepsm(1,i),cepsp(1,i),
cfc     &     scale,redo,redo,rfpenE(k,i,j,2,l(1),l(2),l(3)),rfpenO(k,i,j,2,l(1),l(2),l(3)),rgfpen(k,i,j,2,l(1),l(2),l(3)),
cfc     &     Div)
cfc           redo=-1
c       enddo
c       enddo
c       enddo
c       enddo
c       redo=1
c       enddo
c       enddo
cc g-g-p-g-m
       redo=1
       do j=1,4
       do k=1,6
       do i=typ,typl,-1
       do l1=ll(1),lh(1)
         l(1) = l1
       do l2=ll(2),lh(2)
         l(2) = l2
       do l3=ll(3),lh(3)
         l(3) = l3
         call FermionPenCouplingT1_GM345(
     &     clr(cpl(j),cpp(i),+1),1d0,clr(cpl(j),cpm(i),+1),
     &     clr(cpl(j),cpp(i),-1),1d0,clr(cpl(j),cpm(i),-1),
     &     mf(j),mf(j),mf(j),
     &     pglu(0,perm(k,1)),pglu(0,perm(k,2)),minus_qp(0,1),pglu(0,perm(k,3)),minus_qm(0,1),
     &     ceps(0,l(perm(k,1)),perm(k,1)),ceps(0,l(perm(k,2)),perm(k,2)),cepsp(1,i),ceps(0,l(perm(k,3)),perm(k,3)),cepsm(1,i),
     &     scale,redo,redo,rfpenE(k,i,j,3,l(1),l(2),l(3)),rfpenO(k,i,j,3,l(1),l(2),l(3)),rgfpen(k,i,j,3,l(1),l(2),l(3)),
     &     Div)
           redo=-1
       enddo
       enddo
       enddo
       enddo
       redo=1
       enddo
       enddo
ccc g-g-m-g-p
c       redo=1
c       do j=1,4
c       do k=1,6
c       do i=typ,typl,-1
c       do l1=ll(1),lh(1)
c         l(1) = l1
c       do l2=ll(2),lh(2)
c         l(2) = l2
c       do l3=ll(3),lh(3)
c         l(3) = l3
c
c          rfpenE(3,typ,j,4,ll(1),ll(2),ll(3))=-rfpenE(1,typ,j,3,ll(1),ll(2),ll(3))
c          rfpenE(5,typ,j,4,ll(1),ll(2),ll(3))=-rfpenE(2,typ,j,3,ll(1),ll(2),ll(3))
c          rfpenE(1,typ,j,4,ll(1),ll(2),ll(3))=-rfpenE(3,typ,j,3,ll(1),ll(2),ll(3))
c          rfpenE(6,typ,j,4,ll(1),ll(2),ll(3))=-rfpenE(4,typ,j,3,ll(1),ll(2),ll(3))
c          rfpenE(2,typ,j,4,ll(1),ll(2),ll(3))=-rfpenE(5,typ,j,3,ll(1),ll(2),ll(3))
c          rfpenE(4,typ,j,4,ll(1),ll(2),ll(3))=-rfpenE(6,typ,j,3,ll(1),ll(2),ll(3))
c          
c          rfpenO(3,typ,j,4,ll(1),ll(2),ll(3))=rfpenO(1,typ,j,3,ll(1),ll(2),ll(3))
c          rfpenO(5,typ,j,4,ll(1),ll(2),ll(3))=rfpenO(2,typ,j,3,ll(1),ll(2),ll(3))
c          rfpenO(1,typ,j,4,ll(1),ll(2),ll(3))=rfpenO(3,typ,j,3,ll(1),ll(2),ll(3))
c          rfpenO(6,typ,j,4,ll(1),ll(2),ll(3))=rfpenO(4,typ,j,3,ll(1),ll(2),ll(3))
c          rfpenO(2,typ,j,4,ll(1),ll(2),ll(3))=rfpenO(5,typ,j,3,ll(1),ll(2),ll(3))
c          rfpenO(4,typ,j,4,ll(1),ll(2),ll(3))=rfpenO(6,typ,j,3,ll(1),ll(2),ll(3))

cfc         call FermionPenCouplingT1_GM345(
cfc     &     clr(cpl(j),cpm(i),+1),1d0,clr(cpl(j),cpp(i),+1),
cfc     &     clr(cpl(j),cpm(i),-1),1d0,clr(cpl(j),cpp(i),-1),
cfc     &     mf(j),mf(j),mf(j),
cfc     &     pglu(0,perm(k,1)),pglu(0,perm(k,2)),minus_qm(0,1),pglu(0,perm(k,3)),minus_qp(0,1),
cfc     &     ceps(0,l(perm(k,1)),perm(k,1)),ceps(0,l(perm(k,2)),perm(k,2)),cepsm(1,i),ceps(0,l(perm(k,3)),perm(k,3)),cepsp(1,i),
cfc     &     scale,redo,redo,rfpenE(k,i,j,4,l(1),l(2),l(3)),rfpenO(k,i,j,4,l(1),l(2),l(3)),rgfpen(k,i,j,4,l(1),l(2),l(3)),
cfc     &     Div)
cfc           redo=-1
c      enddo
c       enddo
c       enddo
c       enddo
c       redo=1
c       enddo
c       enddo
  
cfc       do j=1,4
cfc
cfc          print*,
cfc          print*, "j",j
cfc          PRINT*, "1 ODD"      
cfc       write(*,*) "Gauge Pent:", 1,j,1,rfpenO(1,typ,j,1,ll(1),ll(2),ll(3))/rfpenO(6,typ,j,2,ll(1),ll(2),ll(3))-1d0
cfc       write(*,*) "Gauge Pent:", 2,j,1,rfpenO(2,typ,j,1,ll(1),ll(2),ll(3))/rfpenO(4,typ,j,2,ll(1),ll(2),ll(3))-1d0
cfc       write(*,*) "Gauge Pent:", 3,j,1,rfpenO(3,typ,j,1,ll(1),ll(2),ll(3))/rfpenO(5,typ,j,2,ll(1),ll(2),ll(3))-1d0
cfc       write(*,*) "Gauge Pent:", 4,j,1,rfpenO(4,typ,j,1,ll(1),ll(2),ll(3))/rfpenO(2,typ,j,2,ll(1),ll(2),ll(3))-1d0
cfc       write(*,*) "Gauge Pent:", 5,j,1,rfpenO(5,typ,j,1,ll(1),ll(2),ll(3))/rfpenO(3,typ,j,2,ll(1),ll(2),ll(3))-1d0
cfc       write(*,*) "Gauge Pent:", 6,j,1,rfpenO(6,typ,j,1,ll(1),ll(2),ll(3))/rfpenO(1,typ,j,2,ll(1),ll(2),ll(3))-1d0
cfc       PRINT*, "3"
cfc       write(*,*) "Gauge Pent:", 1,j,3,rfpenO(1,typ,j,3,ll(1),ll(2),ll(3))/rfpenO(3,typ,j,4,ll(1),ll(2),ll(3))-1d0
cfc       write(*,*) "Gauge Pent:", 2,j,3,rfpenO(2,typ,j,3,ll(1),ll(2),ll(3))/rfpenO(5,typ,j,4,ll(1),ll(2),ll(3))-1d0
cfc       write(*,*) "Gauge Pent:", 3,j,3,rfpenO(3,typ,j,3,ll(1),ll(2),ll(3))/rfpenO(1,typ,j,4,ll(1),ll(2),ll(3))-1d0
cfc       write(*,*) "Gauge Pent:", 4,j,3,rfpenO(4,typ,j,3,ll(1),ll(2),ll(3))/rfpenO(6,typ,j,4,ll(1),ll(2),ll(3))-1d0
cfc       write(*,*) "Gauge Pent:", 5,j,3,rfpenO(5,typ,j,3,ll(1),ll(2),ll(3))/rfpenO(2,typ,j,4,ll(1),ll(2),ll(3))-1d0
cfc       write(*,*) "Gauge Pent:", 6,j,3,rfpenO(6,typ,j,3,ll(1),ll(2),ll(3))/rfpenO(4,typ,j,4,ll(1),ll(2),ll(3))-1d0
cfc         PRINT*, "1 EVEN"      
cfc       write(*,*) "Gauge Pent:", 1,j,1,-rfpenE(1,typ,j,1,ll(1),ll(2),ll(3))/rfpenE(6,typ,j,2,ll(1),ll(2),ll(3))-1d0
cfc       write(*,*) "Gauge Pent:", 2,j,1,-rfpenE(2,typ,j,1,ll(1),ll(2),ll(3))/rfpenE(4,typ,j,2,ll(1),ll(2),ll(3))-1d0
cfc       write(*,*) "Gauge Pent:", 3,j,1,-rfpenE(3,typ,j,1,ll(1),ll(2),ll(3))/rfpenE(5,typ,j,2,ll(1),ll(2),ll(3))-1d0
cfc       write(*,*) "Gauge Pent:", 4,j,1,-rfpenE(4,typ,j,1,ll(1),ll(2),ll(3))/rfpenE(2,typ,j,2,ll(1),ll(2),ll(3))-1d0
cfc       write(*,*) "Gauge Pent:", 5,j,1,-rfpenE(5,typ,j,1,ll(1),ll(2),ll(3))/rfpenE(3,typ,j,2,ll(1),ll(2),ll(3))-1d0
cfc       write(*,*) "Gauge Pent:", 6,j,1,-rfpenE(6,typ,j,1,ll(1),ll(2),ll(3))/rfpenE(1,typ,j,2,ll(1),ll(2),ll(3))-1d0
cfc       PRINT*, "3"                                                              
cfc       write(*,*) "Gauge Pent:", 1,j,3,-rfpenE(1,typ,j,3,ll(1),ll(2),ll(3))/rfpenE(3,typ,j,4,ll(1),ll(2),ll(3))-1d0
cfc       write(*,*) "Gauge Pent:", 2,j,3,-rfpenE(2,typ,j,3,ll(1),ll(2),ll(3))/rfpenE(5,typ,j,4,ll(1),ll(2),ll(3))-1d0
cfc       write(*,*) "Gauge Pent:", 3,j,3,-rfpenE(3,typ,j,3,ll(1),ll(2),ll(3))/rfpenE(1,typ,j,4,ll(1),ll(2),ll(3))-1d0
cfc       write(*,*) "Gauge Pent:", 4,j,3,-rfpenE(4,typ,j,3,ll(1),ll(2),ll(3))/rfpenE(6,typ,j,4,ll(1),ll(2),ll(3))-1d0
cfc       write(*,*) "Gauge Pent:", 5,j,3,-rfpenE(5,typ,j,3,ll(1),ll(2),ll(3))/rfpenE(2,typ,j,4,ll(1),ll(2),ll(3))-1d0
cfc       write(*,*) "Gauge Pent:", 6,j,3,-rfpenE(6,typ,j,3,ll(1),ll(2),ll(3))/rfpenE(4,typ,j,4,ll(1),ll(2),ll(3))-1d0
cfc       enddo
cfc       pause

       wardidtest1 = .true.
      do j=1,4
       do i=1,4,2 ! Change here, only 1 and 3 are evaluated!!!
         do k=1,6           
           if (rgfpen(k,typ,j,i,ll(1),ll(2),ll(3)).gt.GaugeAccuracyPenL) then 
             wardidtest1=.false.
c             write(*,*) "Gauge Pent:", j,k,i,rgfpen(k,typ,j,i,ll(1),ll(2),ll(3)), singular
             ggZZj = 0d0
             return
           endif
         enddo
        enddo
       enddo
cc add everything together
c       do l1=ll(1),lh(1)
c       do l2=ll(2),lh(2)
c       do l3=ll(3),lh(3)
c         m1 = 0d0
c         m2 = 0d0
c         do i=1,4
c         do j=typl,typ
c         do k=1,6
c          if(signperm(k).eq.1) then
c           m1 = m1 + 2*(rfpenE(k,j,1,i,l1,l2,l3) + rfpenO(k,j,1,i,l1,l2,l3) )! u,c
c     &             + 1*(rfpenE(k,j,3,i,l1,l2,l3) + rfpenO(k,j,3,i,l1,l2,l3) )! t
c     &             + 2*(rfpenE(k,j,2,i,l1,l2,l3) + rfpenO(k,j,2,i,l1,l2,l3) )! d,s
c     &             + 1*(rfpenE(k,j,4,i,l1,l2,l3) + rfpenO(k,j,4,i,l1,l2,l3) )! b
c          else                                  
c           m2 = m2 + 2*(rfpenE(k,j,1,i,l1,l2,l3) + rfpenO(k,j,1,i,l1,l2,l3) )! u,c
c     &             + 1*(rfpenE(k,j,3,i,l1,l2,l3) + rfpenO(k,j,3,i,l1,l2,l3) )! t
c     &             + 2*(rfpenE(k,j,2,i,l1,l2,l3) + rfpenO(k,j,2,i,l1,l2,l3) )! d,s
c     &             + 1*(rfpenE(k,j,4,i,l1,l2,l3) + rfpenO(k,j,4,i,l1,l2,l3) )! b
c          endif
c         enddo
c         enddo
c         enddo
cc colour assignment depending on permutation
c         mat(1,1,l1,l2,l3) = m1 
c         mat(1,2,l1,l2,l3) = m2 
c       enddo
c       enddo
c       enddo
c add everything together
       do l1=ll(1),lh(1)
       do l2=ll(2),lh(2)
       do l3=ll(3),lh(3)
         m1 = 0d0
         m2 = 0d0
         do i=1,3,2
         do j=typl,typ
         do k=1,6
c first sum the even and odd parts individually
           m3 = + 2*rfpenE(k,j,1,i,l1,l2,l3) ! u,c
     &          + 1*rfpenE(k,j,3,i,l1,l2,l3) ! t
     &          + 2*rfpenE(k,j,2,i,l1,l2,l3) ! d,s
     &          + 1*rfpenE(k,j,4,i,l1,l2,l3) ! b
           m4 = + 2*rfpenO(k,j,1,i,l1,l2,l3) ! u,c
     &          + 1*rfpenO(k,j,3,i,l1,l2,l3) ! t
     &          + 2*rfpenO(k,j,2,i,l1,l2,l3) ! d,s
     &          + 1*rfpenO(k,j,4,i,l1,l2,l3) ! b
c then assign to colour structures
          if(signperm(k).eq.1) then
            m1 = m1 + m3 + m4
            m2 = m2 - m3 + m4
          else                                  
            m2 = m2 + m3 + m4
            m1 = m1 - m3 + m4
          endif
         enddo
         enddo
         enddo
c colour assignment depending on permutation
         mat(1,1,l1,l2,l3) = m1 
         mat(1,2,l1,l2,l3) = m2 
       enddo
       enddo
       enddo
      endif ! lpent

      if (singular) then 
c        write(*,*) "Gauge Sing:", singular
        ggZZj = 0d0
        return
      endif

c------------------------------------------------------
c sum the graphs, square them 
      res(1) = 0
      do l1 = ll(1),lh(1)
      do l2 = ll(2),lh(2)
      do l3 = ll(3),lh(3)
       do i=1,2
        mm(i) = 0
        do j = 1,6
           mm(i) = mm(i) +  mat(j,i,l1,l2,l3)
cFC
c         print*, "CCCMAT", i,j, mat(j,i,l1,l2,l3)
cFC           
        enddo
c multiply with (g_s)^3 / (16 pi^2)
        mm(i) = mm(i)*sqrt(als(1,1)**3d0*oneo4pi)
       enddo
       res(1)= res(1) 
     &          + 7/3d0 * ( dreal(mm(1))**2 + dimag(mm(1))**2 )
     &          + 7/3d0 * ( dreal(mm(2))**2 + dimag(mm(2))**2 )
     &          - 4/3d0 * ( dreal(mm(1))*dreal(mm(2)) + dimag(mm(1))*dimag(mm(2)) )
      enddo
      enddo
      enddo
      ggzzj = res(1) * helmult
      return
      end
