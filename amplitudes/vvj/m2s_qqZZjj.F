c*********************  subroutine m2s_qqZZjj *********************************
c
c INPUT     bos       Boson identifier, 2 = Z, 6=H
c           nlo       NLO=1: create subtraction term; LO = 0
c           lok       need to calculate m2s(1:n_kin)? T or F 
c           xi        Feynman x parameters of incoming fermions
c           p         momenta associated with external fermions
c           v         momenta associated with ZZ decay fermions
c           rn        random number array
c           xuz       x_ija for subtraction terms
c          
c
c OUTPUT    m2s       |M|^2*pdf1*pdf2 
c
c       Adapted from m2s_qqVVjj.F, which is adapted from m2s_wbfzhg by D. Zeppenfeld
c
c       Michael Rauch, <michael.rauch@kit.edu>
c	Initial version: October 2014
c	Last modified: October 2014
c
c***************************************************************************

      subroutine m2s_qqZZjj(bos,nlo,lok,xi,p,v,rn,xuz,m2s)

      implicit none

c declare input/output variables

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/lha.inc"

      integer bos
      real*8 xi(nx), p(0:3,max_p,max_kin), v(0:3,max_v,max_kin), 
     &       pch(0:3,max_p,max_kin),
     1       rn(1), xuz(2,2:11), m2s(0:n_kin)
      logical lok(max_kin)

      double precision N
      complex*16 zero
      parameter(N=3d0)
      parameter (zero=(0d0,0d0))

c  helicity selection
      INTEGER h
      COMMON /hcount / h

c declare external functions

      integer FL_VVg
      external FL_VVg

c alfas, scales etc


      double precision  clr, xm2, xmg, b
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6)

c declare local variables

      real*8 q_sf
      real*8 dotrr
      double complex help1c, zm2i(2:3), dotrc
      external dotrr, dotrc
      save zm2i

      integer init/0/, I, J, mu
      save init

c matrix elements
c indices denote 
c   mgg(i,j,k): i: up(1)/down(2) quark line
c               j: permutation initial-state
c               k: permutation in/out particles
c   mqq(i,j,k,l): 
c               i: diagram structures 
c               j: particle 2 of internal structure is up(1)/down(2) 
c               k: particle 4 of internal structure is up(1)/down(2) 
c               l: permutation in/out particles
      double precision mgg(2,2,4), mqq(3,2,2,4)

c declare variables for summation process

      INTEGER ires

      integer nlo, nmax
      save nmax

c store contributions from subprocess ip in res(ip,ID) where
c  ID = 1      : the real emission |M|^2 * pdf
c     = 2,n_kin: subtraction terms for emisions
c     = 0      : subtracted result which drives integration, i.e
c res(*,0) = sum_i=1^n_kin res(*,i)

      real*8 res(maxnumsubproc,0:n_kin)
      integer diag(59)
      double precision polcolqq, polcolgq, polcolgg
      double precision dipole(18)
      double precision fincol(21), fincolglobal, fincolfac
      double precision dipoleVVj_CS39, dipoleVVj_CS40, dipoleVVj_CS41
      double precision dipoleVVj_CS65, dipoleVVj_CS66, dipoleVVj_CS67, dipoleVVj_CS68
      double precision dipoleVVj_CS145, dipoleVVj_CS146, dipoleVVj_CS147, dipoleVVj_CS148
      double complex borncurr(0:3,-1:1,2,62)

c declare parton distribution variables

      real*8 pdf(-6:6,2,n_kin)

c define program switches

      logical ldebug, lfincolls
      double precision kininv(4,4), ratio(6)
      logical kintrue
      parameter (ldebug=.false., lfincolls=.true.)

      real*8 weight,rnumb,RandomNumber



c --------------------------------------------------------------

c if first iteration, output acceptance cuts and fix constant input 
c parameters for call of qqZZjj.f

      if ( init .ne. 0) then
         do i = 1,nmax
            do j=0,n_kin
               res(i,j) = 0
            enddo
         enddo
      else if (init.eq.0) then
         write(6,*) " "
         if (bos.eq.22) then
           write(6,*) "Z Z j j amplitude square information:"
         else if (bos.eq.21) then
           write(6,*) "Z A j j amplitude square information:"
         else if (bos.eq.11) then
           write(6,*) "A A j j amplitude square information:"
         else 
           write(*,*) "Called m2s_qqZZjj.F with invalid boson", bos
           stop
         endif
         write(6,*) "-----------------------------------------------"
         write(6,*) " "
         call dibostoww_reset
         call printnfl(.false.)
         if (nfl.ge.5) then
           lbquark = 1
         else
           lbquark = 0
         endif
         h=1
         init = 1
         do i = 1,maxnumsubproc
            do j = 0,n_kin
               res(i,j) = 0
            enddo
         enddo
         zm2i(2) = 1d0/dcmplx(xm2(2),-xmg(2))
         zm2i(3) = 1d0/dcmplx(xm2(3),-xmg(3))

         nmax = 59
         do i=1,nmax
           diag(i) = 1
         enddo

      endif

c ------------------------------------------------------------------- 

c polarisation vectors
c select helicity: h = random number for helicity (h=1:4)
      h = mod(h,4) + 1
      ie = (-1)**(h+1)
      iu = (-1)**((h-1)/2)

c precalculate leptonic tensors
c lepton spinors and W+- polarization vectors

c...Les Houches interface
      if ((lha.or.hepmc).and..not.doNLO) then
         helicity(1)=-ie
         helicity(2)= ie
         helicity(3)=-iu
         helicity(4)= iu
      endif

c leptonic tensors
c all kinematics which have unshifted lepton momenta (all but initial-initial)
      if (lok(1).or.lok(2).or.lok(3).or.lok(6).or.lok(7).or.lok(10).or.lok(11)) then
         J=1
         if ((bos.eq.22) .or. (bos.eq.21)) then
           CALL OXXXXX(v(0,1,J),ZERO ,-ie,+1,wve(1,J)) !W(1,3))          !e-
           CALL IXXXXX(v(0,2,J),ZERO ,+ie,-1,wep(1,J)) !W(1,4))          !e+
           CALL JIOXXX(wep(1,J),wve(1,J),GZL ,ZMASS,ZWIDTH,wp(1,J))    !Z_e
           CALL JIOXXX(wep(1,J),wve(1,J),GAL ,ZERO ,ZERO  ,atau1(1,J)) !A_e
           if (bos.eq.22) then
             CALL OXXXXX(v(0,3,J),ZERO ,-iu,+1,wmu(1,J)) !W(1,5))          !mu-
             CALL IXXXXX(v(0,4,J),ZERO ,+iu,-1,wvm(1,J)) !W(1,6))          !mu+
             CALL JIOXXX(wvm(1,J),wmu(1,J),GZL ,ZMASS,ZWIDTH,wm(1,J))    !Z_mu
             CALL JIOXXX(wvm(1,J),wmu(1,J),GAL ,ZERO ,ZERO  ,atau2(1,J)) !A_mu
             call zz_ztozz(v(0,1,J),J,wz(1,J),0) !Z -> e+ e- mu+ mu-
             call zz_atozz(v(0,1,J),J,wa(1,J),0) !A -> e+ e- mu+ mu-
           else if (bos.eq.21) then
             CALL VXXXXX(v(0,3,J),ZERO,iu,1,atau2(1,J))                  !A2
             call za_ztoza(v(0,1,J),J,wz(1,J),0) !Z -> e+ e- a
             call za_atoza(v(0,1,J),J,wa(1,J),0) !A -> e+ e- a
           endif

           do mu = 0,3
              qp(mu,J) = v(mu,1,J)+v(mu,2,J)   !Z/A_e
              qm(mu,J) = v(mu,3,J)+v(mu,4,J)   !Z/A_mu
              qaa(mu,J) = qp(mu,J)+qm(mu,J)    !Z_all
           enddo
           qp(4,J) = dotrr(qp(0,J),qp(0,J))
           qm(4,J) = dotrr(qm(0,J),qm(0,J))
           qaa(4,J) = dotrr(qaa(0,J),qaa(0,J))         

c define epsilon~ = espilon - k * k.epsilon/M**2 
           help1c = -zm2i(2) * dotrc(qaa(0,J),wz(1,J))
           do mu=0,3
              wz(mu+1,J) = wz(mu+1,J) + qaa(mu,J) * help1c
           enddo

         else if (bos.eq.11) then
           do mu = 0,3
              qa1(mu,J) = v(mu,1,J)! A1
              qa2(mu,J) = v(mu,2,J)! A2
              qaa(mu,J) = v(mu,1,J) + v(mu,2,J)! A1+A2
           enddo
           qa1(4,1) = dotrr(qa1(0,J),qa1(0,J))
           qa2(4,1) = dotrr(qa2(0,J),qa2(0,J))
           qaa(4,1) = dotrr(qaa(0,J),qaa(0,J))

           CALL VXXXXX(qa1(0,J),ZERO,ie,J,atau1(1,J)) !A1
           CALL VXXXXX(qa2(0,J),ZERO,iu,J,atau2(1,J)) !A2

         endif

         do i=1,6
c copy to others
           if (i.eq.1) j = 2
           if (i.eq.2) j = 3
           if (i.eq.3) j = 6
           if (i.eq.4) j = 7
           if (i.eq.5) j = 10
           if (i.eq.6) j = 11

           do mu=1,6
             wp(mu,j) = wp(mu,1)
             atau1(mu,j) = atau1(mu,1)
             wm(mu,j) = wm(mu,1)
             atau2(mu,j) = atau2(mu,1)
             wz(mu,j) = wz(mu,1)
             wa(mu,j) = wa(mu,1)
           enddo
           do mu = 0,4
              qp(mu,J)  = qp(mu,1)
              qm(mu,J)  = qm(mu,1)
              qaa(mu,J) = qaa(mu,1)
           enddo

         enddo
      endif

c then the four initial-initial kinematics with shifted lepton momenta
      do i=1,4
         if (i.eq.1) j = 4
         if (i.eq.2) j = 5
         if (i.eq.3) j = 8
         if (i.eq.4) j = 9

         if (lok(j)) then
           if ((bos.eq.22) .or. (bos.eq.21)) then
             CALL OXXXXX(v(0,1,J),ZERO ,-ie,+1,wve(1,J)) !W(1,3))          !e-
             CALL IXXXXX(v(0,2,J),ZERO ,+ie,-1,wep(1,J)) !W(1,4))          !e+
             CALL JIOXXX(wep(1,J),wve(1,J),GZL ,ZMASS,ZWIDTH,wp(1,J))    !Z_e
             CALL JIOXXX(wep(1,J),wve(1,J),GAL ,ZERO ,ZERO  ,atau1(1,J)) !A_e
             if (bos.eq.22) then
               CALL OXXXXX(v(0,3,J),ZERO ,-iu,+1,wmu(1,J)) !W(1,5))          !mu-
               CALL IXXXXX(v(0,4,J),ZERO ,+iu,-1,wvm(1,J)) !W(1,6))          !mu+
               CALL JIOXXX(wvm(1,J),wmu(1,J),GZL ,ZMASS,ZWIDTH,wm(1,J))    !Z_mu
               CALL JIOXXX(wvm(1,J),wmu(1,J),GAL ,ZERO ,ZERO  ,atau2(1,J)) !A_mu
               call zz_ztozz(v(0,1,J),J,wz(1,J),0) !Z -> e+ e- mu+ mu-
               call zz_atozz(v(0,1,J),J,wa(1,J),0) !A -> e+ e- mu+ mu-
             else if (bos.eq.21) then
               CALL VXXXXX(v(0,3,J),ZERO,iu,1,atau2(1,J))                  !A2
               call za_ztoza(v(0,1,J),J,wz(1,J),0) !Z -> e+ e- a
               call za_atoza(v(0,1,J),J,wa(1,J),0) !A -> e+ e- a
             endif

             do mu = 0,3
                qp(mu,J) = v(mu,1,J)+v(mu,2,J)   !Z/A_e
                qm(mu,J) = v(mu,3,J)+v(mu,4,J)   !Z/A_mu
                qaa(mu,J) = qp(mu,J)+qm(mu,J)    !Z_all
             enddo
             qp(4,J) = dotrr(qp(0,J),qp(0,J))
             qm(4,J) = dotrr(qm(0,J),qm(0,J))
             qaa(4,J) = dotrr(qaa(0,J),qaa(0,J))         

c define epsilon~ = espilon - k * k.epsilon/M**2 
             help1c = -zm2i(2) * dotrc(qaa(0,J),wz(1,J))
             do mu=0,3
                wz(mu+1,J) = wz(mu+1,J) + qaa(mu,J) * help1c
             enddo

           else if (bos.eq.11) then
             do mu = 0,3
                qa1(mu,J) = v(mu,1,J)! A1
                qa2(mu,J) = v(mu,2,J)! A2
                qaa(mu,J) = v(mu,1,J) + v(mu,2,J)! A1+A2
             enddo
             qa1(4,1) = dotrr(qa1(0,J),qa1(0,J))
             qa2(4,1) = dotrr(qa2(0,J),qa2(0,J))
             qaa(4,1) = dotrr(qaa(0,J),qaa(0,J))

             CALL VXXXXX(qa1(0,J),ZERO,ie,J,atau1(1,J)) !A1
             CALL VXXXXX(qa2(0,J),ZERO,iu,J,atau2(1,J)) !A2

           endif
         endif

      enddo

c call PDF subroutine in order to determine parton
c distributions in the incoming (anti)protons.
c only calculate for the active kinematics

      do i=1,n_kin
        if (lok(i)) then
          q_sf = sqrt(mufsq(1,i))
          call pdfproton( xi(1), q_sf, pdf(-6,1,i) )      ! f_a(y)=f_a(x1/z)

          q_sf = sqrt(mufsq(2,i))
          call pdfproton( xi(2), q_sf, pdf(-6,2,i) )      ! f_b(y)=f_a(x2/z)
        endif
      enddo

      polcolqq = 1d0/(4d0*N**2*xi(1)*xi(2))
      polcolgq = 1d0/(4d0*N*(N**2-1)*xi(1)*xi(2))
      polcolgg = 1d0/(4d0*(N**2-1)**2*xi(1)*xi(2))

C calculate matrix elements

      do mu=0,3
        pch(mu,1,1)=p(mu,2,1)
        pch(mu,2,1)=p(mu,1,1)
        pch(mu,3,1)=p(mu,3,1)
        pch(mu,4,1)=p(mu,4,1)
      enddo

C********************************************************************

C     Now sum over the subprocesses contributing to Z Z production

C********************************************************************

      if (lok(1)) then

         do i=1,4
           call qqZAZAqq(p,v,i,bos,mqq(1,1,1,i))
           call qqZAZAgg(p,v,i,bos,mgg(1,1,i))  
           if (i.ne.4) then
             call qqZAZAgg(pch,v,i,bos,mgg(1,2,i))  
           endif
         enddo

c (uu->uu) + (cc->cc)
         res(1,1) = diag(1) * mqq(3,1,1,1) * 
     &      ( pdf(2*sign1,1,1)*pdf(2*sign2,2,1) + 
     &        pdf(4*sign1,1,1)*pdf(4*sign2,2,1) )/2d0 * polcolqq
c (ud->du) + (cs->sc)
         res(2,1) = diag(2) * mqq(2,2,1,1) * 
     &      ( pdf(2*sign1,1,1)*pdf(1*sign2,2,1) + 
     &        pdf(4*sign1,1,1)*pdf(3*sign2,2,1) ) * polcolqq
c (uc->uc) + (cu->uc)
         res(3,1) = diag(3) * mqq(1,1,1,1) * 
     &      ( pdf(2*sign1,1,1)*pdf(4*sign2,2,1) + 
     &        pdf(4*sign1,1,1)*pdf(2*sign2,2,1) ) * polcolqq
c (us->us) + (ub->ub) + (cd->cd) + (cb->cb)
         res(4,1) = diag(4) * mqq(1,1,2,1) * 
     &      ( pdf(2*sign1,1,1)*pdf(3*sign2,2,1) + 
     &        lbquark*pdf(2*sign1,1,1)*pdf(5*sign2,2,1) + 
     &        pdf(4*sign1,1,1)*pdf(1*sign2,2,1) + 
     &        lbquark*pdf(4*sign1,1,1)*pdf(5*sign2,2,1) ) * polcolqq
c (uS->uS) + (uB->uB) + (cD->cD) + (cB->cB)
         res(6,1) = diag(6) * mqq(1,1,2,2) * 
     &      ( pdf(2*sign1,1,1)*pdf(-3*sign2,2,1) + 
     &        lbquark*pdf(2*sign1,1,1)*pdf(-5*sign2,2,1) + 
     &        pdf(4*sign1,1,1)*pdf(-1*sign2,2,1) + 
     &        lbquark*pdf(4*sign1,1,1)*pdf(-5*sign2,2,1) ) * polcolqq
c (uC->uC) + (cU->cU)
         res(7,1) = diag(7) * mqq(1,1,1,2) * 
     &      ( pdf(2*sign1,1,1)*pdf(-4*sign2,2,1) + 
     &        pdf(4*sign1,1,1)*pdf(-2*sign2,2,1) ) * polcolqq
c (uD->uD) + (cS->cS)
         res(8,1) = diag(8) * mqq(1,1,2,2) * 
     &      ( pdf(2*sign1,1,1)*pdf(-1*sign2,2,1) + 
     &        pdf(4*sign1,1,1)*pdf(-3*sign2,2,1) ) * polcolqq
c (uU->uU) + (cC->cC)
         res(9,1) = diag(9) * mqq(3,1,1,2) * 
     &      ( pdf(2*sign1,1,1)*pdf(-2*sign2,2,1) + 
     &        pdf(4*sign1,1,1)*pdf(-4*sign2,2,1) ) * polcolqq
c (uU->dD) + (cC->sS)
         res(10,1) = diag(10) * mqq(2,2,1,2) * 
     &      ( pdf(2*sign1,1,1)*pdf(-2*sign2,2,1) + 
     &        pdf(4*sign1,1,1)*pdf(-4*sign2,2,1) ) * polcolqq
c (uU->cC) + (cC->uU)
         res(11,1) = diag(11) * mqq(2,1,1,2) * 
     &      ( pdf(2*sign1,1,1)*pdf(-2*sign2,2,1) + 
     &        pdf(4*sign1,1,1)*pdf(-4*sign2,2,1) ) * polcolqq
c (uU->sS) + (uU->bB) + (cC->dD) + (cC->bB)
         res(12,1) = diag(12) * mqq(2,2,1,2) * 
     &      ( pdf(2*sign1,1,1)*pdf(-2*sign2,2,1) + 
     &        lbquark*pdf(2*sign1,1,1)*pdf(-2*sign2,2,1) + 
     &        pdf(4*sign1,1,1)*pdf(-4*sign2,2,1) + 
     &        lbquark*pdf(4*sign1,1,1)*pdf(-4*sign2,2,1) ) * polcolqq
c (du->du) + (sc->sc)
         res(14,1) = diag(14) * mqq(1,2,1,1) * 
     &      ( pdf(1*sign1,1,1)*pdf(2*sign2,2,1) + 
     &        pdf(3*sign1,1,1)*pdf(4*sign2,2,1) ) * polcolqq
c (dd->dd) + (ss->ss) + (bb->bb)
         res(15,1) = diag(15) * mqq(3,2,2,1) * 
     &      ( pdf(1*sign1,1,1)*pdf(1*sign2,2,1) + 
     &        pdf(3*sign1,1,1)*pdf(3*sign2,2,1) + 
     &        lbquark*pdf(5*sign1,1,1)*pdf(5*sign2,2,1) )/2d0 * polcolqq
c (dc->dc) + (su->su) + (bu->bu) + (bc->bc)
         res(16,1) = diag(16) * mqq(1,2,1,1) * 
     &      ( pdf(1*sign1,1,1)*pdf(4*sign2,2,1) + 
     &        pdf(3*sign1,1,1)*pdf(2*sign2,2,1) + 
     &        lbquark*pdf(5*sign1,1,1)*pdf(2*sign2,2,1) + 
     &        lbquark*pdf(5*sign1,1,1)*pdf(4*sign2,2,1) ) * polcolqq
c (ds->ds) + (db->db) + (sd->ds) + (sb->sb) + (bd->db) + (bs->sb)
         res(17,1) = diag(17) * mqq(1,2,2,1) * 
     &      ( pdf(1*sign1,1,1)*pdf(3*sign2,2,1) + 
     &        lbquark*pdf(1*sign1,1,1)*pdf(5*sign2,2,1) + 
     &        pdf(3*sign1,1,1)*pdf(1*sign2,2,1) + 
     &        lbquark*pdf(3*sign1,1,1)*pdf(5*sign2,2,1) + 
     &        lbquark*pdf(5*sign1,1,1)*pdf(1*sign2,2,1) + 
     &        lbquark*pdf(5*sign1,1,1)*pdf(3*sign2,2,1) ) * polcolqq
c (dS->dS) + (dB->dB) + (sD->sD) + (sB->sB) + (bD->bD) + (bS->bS)
         res(19,1) = diag(19) * mqq(1,2,2,2) * 
     &      ( pdf(1*sign1,1,1)*pdf(-3*sign2,2,1) + 
     &        lbquark*pdf(1*sign1,1,1)*pdf(-5*sign2,2,1) + 
     &        pdf(3*sign1,1,1)*pdf(-1*sign2,2,1) + 
     &        lbquark*pdf(3*sign1,1,1)*pdf(-5*sign2,2,1) + 
     &        lbquark*pdf(5*sign1,1,1)*pdf(-1*sign2,2,1) + 
     &        lbquark*pdf(5*sign1,1,1)*pdf(-3*sign2,2,1) ) * polcolqq
c (dC->dC) + (sU->sU) + (bU->bU) + (bC->bC)
         res(20,1) = diag(20) * mqq(1,2,1,2) * 
     &      ( pdf(1*sign1,1,1)*pdf(-4*sign2,2,1) + 
     &        pdf(3*sign1,1,1)*pdf(-2*sign2,2,1) + 
     &        lbquark*pdf(5*sign1,1,1)*pdf(-2*sign2,2,1) + 
     &        lbquark*pdf(5*sign1,1,1)*pdf(-4*sign2,2,1) ) * polcolqq
c (dD->uU) + (sS->cC)
         res(21,1) = diag(21) * mqq(2,1,2,2) * 
     &      ( pdf(1*sign1,1,1)*pdf(-1*sign2,2,1) + 
     &        pdf(3*sign1,1,1)*pdf(-3*sign2,2,1) ) * polcolqq
c (dD->dD) + (sS->sS) + (bB->bB)
         res(22,1) = diag(22) * mqq(3,2,2,2) * 
     &      ( pdf(1*sign1,1,1)*pdf(-1*sign2,2,1) + 
     &        pdf(3*sign1,1,1)*pdf(-3*sign2,2,1) + 
     &        lbquark*pdf(5*sign1,1,1)*pdf(-5*sign2,2,1) ) * polcolqq
c (dD->cC) + (sS->uU) + (bB->uU) + (bB->cC)
         res(23,1) = diag(23) * mqq(2,1,2,2) * 
     &      ( pdf(1*sign1,1,1)*pdf(-1*sign2,2,1) + 
     &        pdf(3*sign1,1,1)*pdf(-3*sign2,2,1) + 
     &        lbquark*pdf(5*sign1,1,1)*pdf(-5*sign2,2,1) + 
     &        lbquark*pdf(5*sign1,1,1)*pdf(-5*sign2,2,1) ) * polcolqq
c (dD->sS) + (dD->bB) + (sS->dD) + (sS->bB) + (bB->dD) + (bB->sS)
         res(24,1) = diag(24) * mqq(2,2,2,2) * 
     &      ( pdf(1*sign1,1,1)*pdf(-1*sign2,2,1) + 
     &        lbquark*pdf(1*sign1,1,1)*pdf(-1*sign2,2,1) + 
     &        pdf(3*sign1,1,1)*pdf(-3*sign2,2,1) + 
     &        lbquark*pdf(3*sign1,1,1)*pdf(-3*sign2,2,1) + 
     &        lbquark*pdf(5*sign1,1,1)*pdf(-5*sign2,2,1) + 
     &        lbquark*pdf(5*sign1,1,1)*pdf(-5*sign2,2,1) ) * polcolqq
c (dU->dU) + (sC->sC)
         res(26,1) = diag(26) * mqq(1,2,1,2) * 
     &      ( pdf(1*sign1,1,1)*pdf(-2*sign2,2,1) + 
     &        pdf(3*sign1,1,1)*pdf(-4*sign2,2,1) ) * polcolqq
c (Dc->cD) + (Su->uS) + (Bu->uB) + (Bc->cB)
         res(33,1) = diag(33) * mqq(2,2,1,3) * 
     &      ( pdf(-1*sign1,1,1)*pdf(4*sign2,2,1) + 
     &        pdf(-3*sign1,1,1)*pdf(2*sign2,2,1) + 
     &        lbquark*pdf(-5*sign1,1,1)*pdf(2*sign2,2,1) + 
     &        lbquark*pdf(-5*sign1,1,1)*pdf(4*sign2,2,1) ) * polcolqq
c (Ds->sD) + (Db->bD) + (Sd->dS) + (Sb->bS) + (Bd->dB) + (Bs->sB)
         res(34,1) = diag(34) * mqq(2,2,2,3) * 
     &      ( pdf(-1*sign1,1,1)*pdf(3*sign2,2,1) + 
     &        lbquark*pdf(-1*sign1,1,1)*pdf(5*sign2,2,1) + 
     &        pdf(-3*sign1,1,1)*pdf(1*sign2,2,1) + 
     &        lbquark*pdf(-3*sign1,1,1)*pdf(5*sign2,2,1) + 
     &        lbquark*pdf(-5*sign1,1,1)*pdf(1*sign2,2,1) + 
     &        lbquark*pdf(-5*sign1,1,1)*pdf(3*sign2,2,1) ) * polcolqq
c (Du->uD) + (Sc->cS)
         res(35,1) = diag(35) * mqq(2,2,1,3) * 
     &      ( pdf(-1*sign1,1,1)*pdf(2*sign2,2,1) + 
     &        pdf(-3*sign1,1,1)*pdf(4*sign2,2,1) ) * polcolqq
c (Dd->cC) + (Ss->uU) + (Bb->uU) + (Bb->cC)
         res(36,1) = diag(36) * mqq(1,2,1,3) * 
     &      ( pdf(-1*sign1,1,1)*pdf(1*sign2,2,1) + 
     &        pdf(-3*sign1,1,1)*pdf(3*sign2,2,1) + 
     &        lbquark*pdf(-5*sign1,1,1)*pdf(5*sign2,2,1) + 
     &        lbquark*pdf(-5*sign1,1,1)*pdf(5*sign2,2,1) ) * polcolqq
c (Dd->sS) + (Dd->bB) + (Ss->dD) + (Ss->bB) + (Bb->dD) + (Bb->sS)
         res(37,1) = diag(37) * mqq(1,2,2,3) * 
     &      ( pdf(-1*sign1,1,1)*pdf(1*sign2,2,1) + 
     &        lbquark*pdf(-1*sign1,1,1)*pdf(1*sign2,2,1) + 
     &        pdf(-3*sign1,1,1)*pdf(3*sign2,2,1) + 
     &        lbquark*pdf(-3*sign1,1,1)*pdf(3*sign2,2,1) + 
     &        lbquark*pdf(-5*sign1,1,1)*pdf(5*sign2,2,1) + 
     &        lbquark*pdf(-5*sign1,1,1)*pdf(5*sign2,2,1) ) * polcolqq
c (Dd->uU) + (Ss->cC)
         res(38,1) = diag(38) * mqq(1,2,1,3) * 
     &      ( pdf(-1*sign1,1,1)*pdf(1*sign2,2,1) + 
     &        pdf(-3*sign1,1,1)*pdf(3*sign2,2,1) ) * polcolqq
c (Dd->dD) + (Ss->sS) + (Bb->bB)
         res(39,1) = diag(39) * mqq(3,2,2,3) * 
     &      ( pdf(-1*sign1,1,1)*pdf(1*sign2,2,1) + 
     &        pdf(-3*sign1,1,1)*pdf(3*sign2,2,1) + 
     &        lbquark*pdf(-5*sign1,1,1)*pdf(5*sign2,2,1) ) * polcolqq
c (DD->DD) + (SS->SS) + (BB->BB)
         res(42,1) = diag(42) * mqq(3,2,2,4) * 
     &      ( pdf(-1*sign1,1,1)*pdf(-1*sign2,2,1) + 
     &        pdf(-3*sign1,1,1)*pdf(-3*sign2,2,1) + 
     &        lbquark*pdf(-5*sign1,1,1)*pdf(-5*sign2,2,1) )/2d0 * polcolqq
c (DU->UD) + (SC->CS)
         res(43,1) = diag(43) * mqq(1,1,2,4) * 
     &      ( pdf(-1*sign1,1,1)*pdf(-2*sign2,2,1) + 
     &        pdf(-3*sign1,1,1)*pdf(-4*sign2,2,1) ) * polcolqq
c (DS->SD) + (DB->BD) + (SD->SD) + (SB->BS) + (BD->BD) + (BS->BS)
         res(44,1) = diag(44) * mqq(1,2,2,4) * 
     &      ( pdf(-1*sign1,1,1)*pdf(-3*sign2,2,1) + 
     &        lbquark*pdf(-1*sign1,1,1)*pdf(-5*sign2,2,1) + 
     &        pdf(-3*sign1,1,1)*pdf(-1*sign2,2,1) + 
     &        lbquark*pdf(-3*sign1,1,1)*pdf(-5*sign2,2,1) + 
     &        lbquark*pdf(-5*sign1,1,1)*pdf(-1*sign2,2,1) + 
     &        lbquark*pdf(-5*sign1,1,1)*pdf(-3*sign2,2,1) ) * polcolqq
c (DC->CD) + (SU->US) + (BU->UB) + (BC->CB)
         res(45,1) = diag(45) * mqq(1,1,2,4) * 
     &      ( pdf(-1*sign1,1,1)*pdf(-4*sign2,2,1) + 
     &        pdf(-3*sign1,1,1)*pdf(-2*sign2,2,1) + 
     &        lbquark*pdf(-5*sign1,1,1)*pdf(-2*sign2,2,1) + 
     &        lbquark*pdf(-5*sign1,1,1)*pdf(-4*sign2,2,1) ) * polcolqq
c (Uc->cU) + (Cu->uC)
         res(46,1) = diag(46) * mqq(2,1,1,3) * 
     &      ( pdf(-2*sign1,1,1)*pdf(4*sign2,2,1) + 
     &        pdf(-4*sign1,1,1)*pdf(2*sign2,2,1) ) * polcolqq
c (Us->sU) + (Ub->bU) + (Cd->dC) + (Cb->bC)
         res(47,1) = diag(47) * mqq(2,1,2,3) * 
     &      ( pdf(-2*sign1,1,1)*pdf(3*sign2,2,1) + 
     &        lbquark*pdf(-2*sign1,1,1)*pdf(5*sign2,2,1) + 
     &        pdf(-4*sign1,1,1)*pdf(1*sign2,2,1) + 
     &        lbquark*pdf(-4*sign1,1,1)*pdf(5*sign2,2,1) ) * polcolqq
c (Uu->cC) + (Cc->uU)
         res(48,1) = diag(48) * mqq(1,1,1,3) * 
     &      ( pdf(-2*sign1,1,1)*pdf(2*sign2,2,1) + 
     &        pdf(-4*sign1,1,1)*pdf(4*sign2,2,1) ) * polcolqq
c (Uu->sS) + (Uu->bB) + (Cc->dD) + (Cc->bB)
         res(49,1) = diag(49) * mqq(1,1,2,3) * 
     &      ( pdf(-2*sign1,1,1)*pdf(2*sign2,2,1) + 
     &        lbquark*pdf(-2*sign1,1,1)*pdf(2*sign2,2,1) + 
     &        pdf(-4*sign1,1,1)*pdf(4*sign2,2,1) + 
     &        lbquark*pdf(-4*sign1,1,1)*pdf(4*sign2,2,1) ) * polcolqq
c (Uu->uU) + (Cc->cC)
         res(50,1) = diag(50) * mqq(3,1,1,3) * 
     &      ( pdf(-2*sign1,1,1)*pdf(2*sign2,2,1) + 
     &        pdf(-4*sign1,1,1)*pdf(4*sign2,2,1) ) * polcolqq
c (Uu->dD) + (Cc->sS)
         res(51,1) = diag(51) * mqq(1,1,2,3) * 
     &      ( pdf(-2*sign1,1,1)*pdf(2*sign2,2,1) + 
     &        pdf(-4*sign1,1,1)*pdf(4*sign2,2,1) ) * polcolqq
c (Ud->dU) + (Cs->sC)
         res(53,1) = diag(53) * mqq(2,1,2,3) * 
     &      ( pdf(-2*sign1,1,1)*pdf(1*sign2,2,1) + 
     &        pdf(-4*sign1,1,1)*pdf(3*sign2,2,1) ) * polcolqq
c (UD->UD) + (CS->CS)
         res(55,1) = diag(55) * mqq(2,2,1,4) * 
     &      ( pdf(-2*sign1,1,1)*pdf(-1*sign2,2,1) + 
     &        pdf(-4*sign1,1,1)*pdf(-3*sign2,2,1) ) * polcolqq
c (UU->UU) + (CC->CC)
         res(56,1) = diag(56) * mqq(3,1,1,4) * 
     &      ( pdf(-2*sign1,1,1)*pdf(-2*sign2,2,1) + 
     &        pdf(-4*sign1,1,1)*pdf(-4*sign2,2,1) )/2d0 * polcolqq
c (US->SU) + (UB->BU) + (CD->DC) + (CB->BC)
         res(57,1) = diag(57) * mqq(1,2,1,4) * 
     &      ( pdf(-2*sign1,1,1)*pdf(-3*sign2,2,1) + 
     &        lbquark*pdf(-2*sign1,1,1)*pdf(-5*sign2,2,1) + 
     &        pdf(-4*sign1,1,1)*pdf(-1*sign2,2,1) + 
     &        lbquark*pdf(-4*sign1,1,1)*pdf(-5*sign2,2,1) ) * polcolqq
c (UC->CU) + (CU->CU)
         res(58,1) = diag(58) * mqq(1,1,1,4) * 
     &      ( pdf(-2*sign1,1,1)*pdf(-4*sign2,2,1) + 
     &        pdf(-4*sign1,1,1)*pdf(-2*sign2,2,1) ) * polcolqq
c (uU->gg) + (cC->gg)
         res(13,1) = diag(13) * mgg(1,1,1) * 
     &      ( pdf(2*sign1,1,1)*pdf(-2*sign2,2,1) + 
     &        pdf(4*sign1,1,1)*pdf(-4*sign2,2,1) )/2d0 * polcolqq
c (dD->gg) + (sS->gg) + (bB->gg)
         res(25,1) = diag(25) * mgg(2,1,1) * 
     &      ( pdf(1*sign1,1,1)*pdf(-1*sign2,2,1) + 
     &        pdf(3*sign1,1,1)*pdf(-3*sign2,2,1) + 
     &        lbquark*pdf(5*sign1,1,1)*pdf(-5*sign2,2,1) )/2d0 * polcolqq
c (Dd->gg) + (Ss->gg) + (Bb->gg)
         res(40,1) = diag(40) * mgg(2,2,1) * 
     &      ( pdf(-1*sign1,1,1)*pdf(1*sign2,2,1) + 
     &        pdf(-3*sign1,1,1)*pdf(3*sign2,2,1) + 
     &        lbquark*pdf(-5*sign1,1,1)*pdf(5*sign2,2,1) )/2d0 * polcolqq
c (Uu->gg) + (Cc->gg)
         res(52,1) = diag(52) * mgg(1,2,1) * 
     &      ( pdf(-2*sign1,1,1)*pdf(2*sign2,2,1) + 
     &        pdf(-4*sign1,1,1)*pdf(4*sign2,2,1) )/2d0 * polcolqq

c (ug->ug) + (cg->cg)
         res(5,1) = diag(5) * mgg(1,1,2) * 
     &      ( pdf(2*sign1,1,1)*pdf(0*sign2,2,1) + 
     &        pdf(4*sign1,1,1)*pdf(0*sign2,2,1) ) * polcolgq
c (dg->dg) + (sg->sg) + (bg->bg)
         res(18,1) = diag(18) * mgg(2,1,2) * 
     &      ( pdf(1*sign1,1,1)*pdf(0*sign2,2,1) + 
     &        pdf(3*sign1,1,1)*pdf(0*sign2,2,1) + 
     &        lbquark*pdf(5*sign1,1,1)*pdf(0*sign2,2,1) ) * polcolgq
c (gu->ug) + (gc->cg)
         res(27,1) = diag(27) * mgg(1,2,2) * 
     &      ( pdf(0*sign1,1,1)*pdf(2*sign2,2,1) + 
     &        pdf(0*sign1,1,1)*pdf(4*sign2,2,1) ) * polcolgq
c (gd->dg) + (gs->sg) + (gb->bg)
         res(28,1) = diag(28) * mgg(2,2,2) * 
     &      ( pdf(0*sign1,1,1)*pdf(1*sign2,2,1) + 
     &        pdf(0*sign1,1,1)*pdf(3*sign2,2,1) + 
     &        lbquark*pdf(0*sign1,1,1)*pdf(5*sign2,2,1) ) * polcolgq
c (gD->Dg) + (gS->Sg) + (gB->Bg)
         res(31,1) = diag(31) * mgg(2,2,3) * 
     &      ( pdf(0*sign1,1,1)*pdf(-1*sign2,2,1) + 
     &        pdf(0*sign1,1,1)*pdf(-3*sign2,2,1) + 
     &        lbquark*pdf(0*sign1,1,1)*pdf(-5*sign2,2,1) ) * polcolgq
c (gU->Ug) + (gC->Cg)
         res(32,1) = diag(32) * mgg(1,2,3) * 
     &      ( pdf(0*sign1,1,1)*pdf(-2*sign2,2,1) + 
     &        pdf(0*sign1,1,1)*pdf(-4*sign2,2,1) ) * polcolgq
c (Dg->Dg) + (Sg->Sg) + (Bg->Bg)
         res(41,1) = diag(41) * mgg(2,1,3) * 
     &      ( pdf(-1*sign1,1,1)*pdf(0*sign2,2,1) + 
     &        pdf(-3*sign1,1,1)*pdf(0*sign2,2,1) + 
     &        lbquark*pdf(-5*sign1,1,1)*pdf(0*sign2,2,1) ) * polcolgq
c (Ug->Ug) + (Cg->Cg)
         res(54,1) = diag(54) * mgg(1,1,3) * 
     &      ( pdf(-2*sign1,1,1)*pdf(0*sign2,2,1) + 
     &        pdf(-4*sign1,1,1)*pdf(0*sign2,2,1) ) * polcolgq

c (gg->uU) + (gg->cC)
         res(29,1) = diag(29) * mgg(1,1,4) * 
     &      ( pdf(0*sign1,1,1)*pdf(0*sign2,2,1) + 
     &        pdf(0*sign1,1,1)*pdf(0*sign2,2,1) ) * polcolgg
c (gg->dD) + (gg->sS) + (gg->bB)
         res(30,1) = diag(30) * mgg(2,1,4) * 
     &      ( pdf(0*sign1,1,1)*pdf(0*sign2,2,1) + 
     &        pdf(0*sign1,1,1)*pdf(0*sign2,2,1) + 
     &        lbquark*pdf(0*sign1,1,1)*pdf(0*sign2,2,1) ) * polcolgg

       endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c dipoles
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       if (lok(2)) then 
         call qqgZAZA_curr(borncurr(0,-1,1,1),2,p(0,3,2),p(0,1,2),p(0,2,2),-1,-1,-1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,2),2,p(0,2,2),p(0,1,2),p(0,3,2),1,-1,1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,3),2,p(0,3,2),p(0,2,2),p(0,1,2),-1,-1,-1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,4),2,p(0,2,2),p(0,3,2),p(0,1,2),1,1,-1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,5),2,p(0,1,2),p(0,2,2),p(0,3,2),1,-1,1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,6),2,p(0,1,2),p(0,3,2),p(0,2,2),1,1,-1,bos)

         dipole(1) = - dipoleVVj_CS67(2,1,p(0,1,1),p(0,4,1),p(0,3,1),xuz(1,2),borncurr(0,-1,1,4))
         dipole(2) = - dipoleVVj_CS67(2,1,p(0,1,1),p(0,4,1),p(0,3,1),xuz(1,2),borncurr(0,-1,2,4))
         dipole(3) = - dipoleVVj_CS65(2,1,p(0,1,1),p(0,4,1),p(0,3,1),xuz(1,2),borncurr(0,0,1,6))
         dipole(4) = - dipoleVVj_CS65(2,0,p(0,1,1),p(0,4,1),p(0,3,1),xuz(1,2),borncurr(0,0,1,5))
         dipole(5) = - dipoleVVj_CS65(2,1,p(0,1,1),p(0,4,1),p(0,3,1),xuz(1,2),borncurr(0,0,2,6))
         dipole(6) = - dipoleVVj_CS65(2,0,p(0,1,1),p(0,4,1),p(0,3,1),xuz(1,2),borncurr(0,0,2,5))
         dipole(7) = - dipoleVVj_CS68(2,1,p(0,1,1),p(0,4,1),p(0,3,1),xuz(1,2),borncurr(0,-1,1,4))
         dipole(8) = - dipoleVVj_CS68(2,1,p(0,1,1),p(0,4,1),p(0,3,1),xuz(1,2),borncurr(0,-1,2,4))
         dipole(9) = - dipoleVVj_CS66(2,1,p(0,1,1),p(0,4,1),p(0,3,1),xuz(1,2),borncurr(0,0,1,6))
         dipole(10) = - dipoleVVj_CS66(2,1,p(0,1,1),p(0,4,1),p(0,3,1),xuz(1,2),borncurr(0,0,2,6))
         dipole(11) = - dipoleVVj_CS68(2,1,p(0,1,1),p(0,4,1),p(0,3,1),xuz(1,2),borncurr(0,-1,2,3))
         dipole(12) = - dipoleVVj_CS68(2,1,p(0,1,1),p(0,4,1),p(0,3,1),xuz(1,2),borncurr(0,-1,1,3))
         dipole(13) = - dipoleVVj_CS65(2,0,p(0,1,1),p(0,4,1),p(0,3,1),xuz(1,2),borncurr(0,0,2,2))
         dipole(14) = - dipoleVVj_CS65(2,1,p(0,1,1),p(0,4,1),p(0,3,1),xuz(1,2),borncurr(0,0,2,1))
         dipole(15) = - dipoleVVj_CS67(2,1,p(0,1,1),p(0,4,1),p(0,3,1),xuz(1,2),borncurr(0,-1,2,3))
         dipole(16) = - dipoleVVj_CS67(2,1,p(0,1,1),p(0,4,1),p(0,3,1),xuz(1,2),borncurr(0,-1,1,3))
         dipole(17) = - dipoleVVj_CS65(2,0,p(0,1,1),p(0,4,1),p(0,3,1),xuz(1,2),borncurr(0,0,1,2))
         dipole(18) = - dipoleVVj_CS65(2,1,p(0,1,1),p(0,4,1),p(0,3,1),xuz(1,2),borncurr(0,0,1,1))

         res(1,2) = diag(1) * dipole(1) * 
     &      ( pdf(2*sign1,1,2)*pdf(2*sign2,2,2) + 
     &        pdf(4*sign1,1,2)*pdf(4*sign2,2,2) )/2d0*polcolgq
         res(2,2) = diag(2) * dipole(2) * 
     &      ( pdf(2*sign1,1,2)*pdf(1*sign2,2,2) + 
     &        pdf(4*sign1,1,2)*pdf(3*sign2,2,2) )*polcolgq
         res(5,2) = diag(5) * dipole(3) * 
     &      ( pdf(2*sign1,1,2)*pdf(0*sign2,2,2) + 
     &        pdf(4*sign1,1,2)*pdf(0*sign2,2,2) )*polcolgq
         res(13,2) = diag(13) * dipole(4) * 
     &      ( pdf(2*sign1,1,2)*pdf(-2*sign2,2,2) + 
     &        pdf(4*sign1,1,2)*pdf(-4*sign2,2,2) )/2d0*polcolqq
         res(15,2) = diag(15) * dipole(2) * 
     &      ( pdf(1*sign1,1,2)*pdf(1*sign2,2,2) + 
     &        pdf(3*sign1,1,2)*pdf(3*sign2,2,2) + 
     &        lbquark*pdf(5*sign1,1,2)*pdf(5*sign2,2,2) )/2d0*polcolgq
         res(18,2) = diag(18) * dipole(5) * 
     &      ( pdf(1*sign1,1,2)*pdf(0*sign2,2,2) + 
     &        pdf(3*sign1,1,2)*pdf(0*sign2,2,2) + 
     &        lbquark*pdf(5*sign1,1,2)*pdf(0*sign2,2,2) )*polcolgq
         res(25,2) = diag(25) * dipole(6) * 
     &      ( pdf(1*sign1,1,2)*pdf(-1*sign2,2,2) + 
     &        pdf(3*sign1,1,2)*pdf(-3*sign2,2,2) + 
     &        lbquark*pdf(5*sign1,1,2)*pdf(-5*sign2,2,2) )/2d0*polcolqq
         res(27,2) = diag(27) * dipole(7) * 
     &      ( pdf(0*sign1,1,2)*pdf(2*sign2,2,2) + 
     &        pdf(0*sign1,1,2)*pdf(4*sign2,2,2) )*polcolgq
         res(28,2) = diag(28) * dipole(8) * 
     &      ( pdf(0*sign1,1,2)*pdf(1*sign2,2,2) + 
     &        pdf(0*sign1,1,2)*pdf(3*sign2,2,2) + 
     &        lbquark*pdf(0*sign1,1,2)*pdf(5*sign2,2,2) )*polcolgq
         res(29,2) = diag(29) * dipole(9) * 
     &      ( pdf(0*sign1,1,2)*pdf(0*sign2,2,2) + 
     &        pdf(0*sign1,1,2)*pdf(0*sign2,2,2) )*polcolgq
         res(30,2) = diag(30) * dipole(10) * 
     &      ( pdf(0*sign1,1,2)*pdf(0*sign2,2,2) + 
     &        pdf(0*sign1,1,2)*pdf(0*sign2,2,2) + 
     &        lbquark*pdf(0*sign1,1,2)*pdf(0*sign2,2,2) )*polcolgq
         res(31,2) = diag(31) * dipole(11) * 
     &      ( pdf(0*sign1,1,2)*pdf(-1*sign2,2,2) + 
     &        pdf(0*sign1,1,2)*pdf(-3*sign2,2,2) + 
     &        lbquark*pdf(0*sign1,1,2)*pdf(-5*sign2,2,2) )*polcolgq
         res(32,2) = diag(32) * dipole(12) * 
     &      ( pdf(0*sign1,1,2)*pdf(-2*sign2,2,2) + 
     &        pdf(0*sign1,1,2)*pdf(-4*sign2,2,2) )*polcolgq
         res(33,2) = diag(33) * dipole(1) * 
     &      ( pdf(-1*sign1,1,2)*pdf(4*sign2,2,2) + 
     &        pdf(-3*sign1,1,2)*pdf(2*sign2,2,2) + 
     &        lbquark*pdf(-5*sign1,1,2)*pdf(2*sign2,2,2) + 
     &        lbquark*pdf(-5*sign1,1,2)*pdf(4*sign2,2,2) )*polcolgq
         res(34,2) = diag(34) * dipole(2) * 
     &      ( pdf(-1*sign1,1,2)*pdf(3*sign2,2,2) + 
     &        lbquark*pdf(-1*sign1,1,2)*pdf(5*sign2,2,2) + 
     &        pdf(-3*sign1,1,2)*pdf(1*sign2,2,2) + 
     &        lbquark*pdf(-3*sign1,1,2)*pdf(5*sign2,2,2) + 
     &        lbquark*pdf(-5*sign1,1,2)*pdf(1*sign2,2,2) + 
     &        lbquark*pdf(-5*sign1,1,2)*pdf(3*sign2,2,2) )*polcolgq
         res(35,2) = diag(35) * dipole(1) * 
     &      ( pdf(-1*sign1,1,2)*pdf(2*sign2,2,2) + 
     &        pdf(-3*sign1,1,2)*pdf(4*sign2,2,2) )*polcolgq
         res(39,2) = diag(39) * dipole(2) * 
     &      ( pdf(-1*sign1,1,2)*pdf(1*sign2,2,2) + 
     &        pdf(-3*sign1,1,2)*pdf(3*sign2,2,2) + 
     &        lbquark*pdf(-5*sign1,1,2)*pdf(5*sign2,2,2) )*polcolgq
         res(40,2) = diag(40) * dipole(13) * 
     &      ( pdf(-1*sign1,1,2)*pdf(1*sign2,2,2) + 
     &        pdf(-3*sign1,1,2)*pdf(3*sign2,2,2) + 
     &        lbquark*pdf(-5*sign1,1,2)*pdf(5*sign2,2,2) )/2d0*polcolqq
         res(41,2) = diag(41) * dipole(14) * 
     &      ( pdf(-1*sign1,1,2)*pdf(0*sign2,2,2) + 
     &        pdf(-3*sign1,1,2)*pdf(0*sign2,2,2) + 
     &        lbquark*pdf(-5*sign1,1,2)*pdf(0*sign2,2,2) )*polcolgq
         res(42,2) = diag(42) * dipole(15) * 
     &      ( pdf(-1*sign1,1,2)*pdf(-1*sign2,2,2) + 
     &        pdf(-3*sign1,1,2)*pdf(-3*sign2,2,2) + 
     &        lbquark*pdf(-5*sign1,1,2)*pdf(-5*sign2,2,2) )/2d0*polcolgq
         res(43,2) = diag(43) * dipole(16) * 
     &      ( pdf(-1*sign1,1,2)*pdf(-2*sign2,2,2) + 
     &        pdf(-3*sign1,1,2)*pdf(-4*sign2,2,2) )*polcolgq
         res(44,2) = diag(44) * dipole(15) * 
     &      ( pdf(-1*sign1,1,2)*pdf(-3*sign2,2,2) + 
     &        lbquark*pdf(-1*sign1,1,2)*pdf(-5*sign2,2,2) + 
     &        pdf(-3*sign1,1,2)*pdf(-1*sign2,2,2) + 
     &        lbquark*pdf(-3*sign1,1,2)*pdf(-5*sign2,2,2) + 
     &        lbquark*pdf(-5*sign1,1,2)*pdf(-1*sign2,2,2) + 
     &        lbquark*pdf(-5*sign1,1,2)*pdf(-3*sign2,2,2) )*polcolgq
         res(45,2) = diag(45) * dipole(16) * 
     &      ( pdf(-1*sign1,1,2)*pdf(-4*sign2,2,2) + 
     &        pdf(-3*sign1,1,2)*pdf(-2*sign2,2,2) + 
     &        lbquark*pdf(-5*sign1,1,2)*pdf(-2*sign2,2,2) + 
     &        lbquark*pdf(-5*sign1,1,2)*pdf(-4*sign2,2,2) )*polcolgq
         res(46,2) = diag(46) * dipole(1) * 
     &      ( pdf(-2*sign1,1,2)*pdf(4*sign2,2,2) + 
     &        pdf(-4*sign1,1,2)*pdf(2*sign2,2,2) )*polcolgq
         res(47,2) = diag(47) * dipole(2) * 
     &      ( pdf(-2*sign1,1,2)*pdf(3*sign2,2,2) + 
     &        lbquark*pdf(-2*sign1,1,2)*pdf(5*sign2,2,2) + 
     &        pdf(-4*sign1,1,2)*pdf(1*sign2,2,2) + 
     &        lbquark*pdf(-4*sign1,1,2)*pdf(5*sign2,2,2) )*polcolgq
         res(50,2) = diag(50) * dipole(1) * 
     &      ( pdf(-2*sign1,1,2)*pdf(2*sign2,2,2) + 
     &        pdf(-4*sign1,1,2)*pdf(4*sign2,2,2) )*polcolgq
         res(52,2) = diag(52) * dipole(17) * 
     &      ( pdf(-2*sign1,1,2)*pdf(2*sign2,2,2) + 
     &        pdf(-4*sign1,1,2)*pdf(4*sign2,2,2) )/2d0*polcolqq
         res(53,2) = diag(53) * dipole(2) * 
     &      ( pdf(-2*sign1,1,2)*pdf(1*sign2,2,2) + 
     &        pdf(-4*sign1,1,2)*pdf(3*sign2,2,2) )*polcolgq
         res(54,2) = diag(54) * dipole(18) * 
     &      ( pdf(-2*sign1,1,2)*pdf(0*sign2,2,2) + 
     &        pdf(-4*sign1,1,2)*pdf(0*sign2,2,2) )*polcolgq
         res(56,2) = diag(56) * dipole(16) * 
     &      ( pdf(-2*sign1,1,2)*pdf(-2*sign2,2,2) + 
     &        pdf(-4*sign1,1,2)*pdf(-4*sign2,2,2) )/2d0*polcolgq
         res(57,2) = diag(57) * dipole(15) * 
     &      ( pdf(-2*sign1,1,2)*pdf(-3*sign2,2,2) + 
     &        lbquark*pdf(-2*sign1,1,2)*pdf(-5*sign2,2,2) + 
     &        pdf(-4*sign1,1,2)*pdf(-1*sign2,2,2) + 
     &        lbquark*pdf(-4*sign1,1,2)*pdf(-5*sign2,2,2) )*polcolgq
         res(58,2) = diag(58) * dipole(16) * 
     &      ( pdf(-2*sign1,1,2)*pdf(-4*sign2,2,2) + 
     &        pdf(-4*sign1,1,2)*pdf(-2*sign2,2,2) )*polcolgq
       endif 

       if (lok(3)) then 
         call qqgZAZA_curr(borncurr(0,-1,1,7),3,p(0,3,3),p(0,1,3),p(0,2,3),-1,-1,-1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,8),3,p(0,2,3),p(0,1,3),p(0,3,3),1,-1,1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,9),3,p(0,3,3),p(0,2,3),p(0,1,3),-1,-1,-1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,10),3,p(0,2,3),p(0,3,3),p(0,1,3),1,1,-1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,11),3,p(0,1,3),p(0,2,3),p(0,3,3),1,-1,1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,12),3,p(0,1,3),p(0,3,3),p(0,2,3),1,1,-1,bos)

         dipole(1) = - dipoleVVj_CS67(3,1,p(0,2,1),p(0,4,1),p(0,3,1),xuz(1,3),borncurr(0,-1,1,12))
         dipole(2) = - dipoleVVj_CS68(3,1,p(0,2,1),p(0,4,1),p(0,3,1),xuz(1,3),borncurr(0,-1,1,12))
         dipole(3) = - dipoleVVj_CS65(3,0,p(0,2,1),p(0,4,1),p(0,3,1),xuz(1,3),borncurr(0,0,1,11))
         dipole(4) = - dipoleVVj_CS67(3,1,p(0,2,1),p(0,4,1),p(0,3,1),xuz(1,3),borncurr(0,-1,2,12))
         dipole(5) = - dipoleVVj_CS68(3,1,p(0,2,1),p(0,4,1),p(0,3,1),xuz(1,3),borncurr(0,-1,2,12))
         dipole(6) = - dipoleVVj_CS65(3,0,p(0,2,1),p(0,4,1),p(0,3,1),xuz(1,3),borncurr(0,0,2,11))
         dipole(7) = - dipoleVVj_CS65(3,1,p(0,2,1),p(0,4,1),p(0,3,1),xuz(1,3),borncurr(0,0,1,10))
         dipole(8) = - dipoleVVj_CS65(3,1,p(0,2,1),p(0,4,1),p(0,3,1),xuz(1,3),borncurr(0,0,2,10))
         dipole(9) = - dipoleVVj_CS66(3,1,p(0,2,1),p(0,4,1),p(0,3,1),xuz(1,3),borncurr(0,0,1,10))
         dipole(10) = - dipoleVVj_CS66(3,1,p(0,2,1),p(0,4,1),p(0,3,1),xuz(1,3),borncurr(0,0,2,10))
         dipole(11) = - dipoleVVj_CS65(3,1,p(0,2,1),p(0,4,1),p(0,3,1),xuz(1,3),borncurr(0,0,2,9))
         dipole(12) = - dipoleVVj_CS65(3,1,p(0,2,1),p(0,4,1),p(0,3,1),xuz(1,3),borncurr(0,0,1,9))
         dipole(13) = - dipoleVVj_CS65(3,0,p(0,2,1),p(0,4,1),p(0,3,1),xuz(1,3),borncurr(0,0,2,8))
         dipole(14) = - dipoleVVj_CS68(3,1,p(0,2,1),p(0,4,1),p(0,3,1),xuz(1,3),borncurr(0,-1,2,7))
         dipole(15) = - dipoleVVj_CS67(3,1,p(0,2,1),p(0,4,1),p(0,3,1),xuz(1,3),borncurr(0,-1,2,7))
         dipole(16) = - dipoleVVj_CS65(3,0,p(0,2,1),p(0,4,1),p(0,3,1),xuz(1,3),borncurr(0,0,1,8))
         dipole(17) = - dipoleVVj_CS68(3,1,p(0,2,1),p(0,4,1),p(0,3,1),xuz(1,3),borncurr(0,-1,1,7))
         dipole(18) = - dipoleVVj_CS67(3,1,p(0,2,1),p(0,4,1),p(0,3,1),xuz(1,3),borncurr(0,-1,1,7))

         res(1,3) = diag(1) * dipole(1) * 
     &      ( pdf(2*sign1,1,3)*pdf(2*sign2,2,3) + 
     &        pdf(4*sign1,1,3)*pdf(4*sign2,2,3) )/2d0*polcolgq
         res(3,3) = diag(3) * dipole(1) * 
     &      ( pdf(2*sign1,1,3)*pdf(4*sign2,2,3) + 
     &        pdf(4*sign1,1,3)*pdf(2*sign2,2,3) )*polcolgq
         res(4,3) = diag(4) * dipole(1) * 
     &      ( pdf(2*sign1,1,3)*pdf(3*sign2,2,3) + 
     &        lbquark*pdf(2*sign1,1,3)*pdf(5*sign2,2,3) + 
     &        pdf(4*sign1,1,3)*pdf(1*sign2,2,3) + 
     &        lbquark*pdf(4*sign1,1,3)*pdf(5*sign2,2,3) )*polcolgq
         res(5,3) = diag(5) * dipole(2) * 
     &      ( pdf(2*sign1,1,3)*pdf(0*sign2,2,3) + 
     &        pdf(4*sign1,1,3)*pdf(0*sign2,2,3) )*polcolgq
         res(6,3) = diag(6) * dipole(1) * 
     &      ( pdf(2*sign1,1,3)*pdf(-3*sign2,2,3) + 
     &        lbquark*pdf(2*sign1,1,3)*pdf(-5*sign2,2,3) + 
     &        pdf(4*sign1,1,3)*pdf(-1*sign2,2,3) + 
     &        lbquark*pdf(4*sign1,1,3)*pdf(-5*sign2,2,3) )*polcolgq
         res(7,3) = diag(7) * dipole(1) * 
     &      ( pdf(2*sign1,1,3)*pdf(-4*sign2,2,3) + 
     &        pdf(4*sign1,1,3)*pdf(-2*sign2,2,3) )*polcolgq
         res(8,3) = diag(8) * dipole(1) * 
     &      ( pdf(2*sign1,1,3)*pdf(-1*sign2,2,3) + 
     &        pdf(4*sign1,1,3)*pdf(-3*sign2,2,3) )*polcolgq
         res(9,3) = diag(9) * dipole(1) * 
     &      ( pdf(2*sign1,1,3)*pdf(-2*sign2,2,3) + 
     &        pdf(4*sign1,1,3)*pdf(-4*sign2,2,3) )*polcolgq
         res(13,3) = diag(13) * dipole(3) * 
     &      ( pdf(2*sign1,1,3)*pdf(-2*sign2,2,3) + 
     &        pdf(4*sign1,1,3)*pdf(-4*sign2,2,3) )/2d0*polcolqq
         res(14,3) = diag(14) * dipole(4) * 
     &      ( pdf(1*sign1,1,3)*pdf(2*sign2,2,3) + 
     &        pdf(3*sign1,1,3)*pdf(4*sign2,2,3) )*polcolgq
         res(15,3) = diag(15) * dipole(4) * 
     &      ( pdf(1*sign1,1,3)*pdf(1*sign2,2,3) + 
     &        pdf(3*sign1,1,3)*pdf(3*sign2,2,3) + 
     &        lbquark*pdf(5*sign1,1,3)*pdf(5*sign2,2,3) )/2d0*polcolgq
         res(16,3) = diag(16) * dipole(4) * 
     &      ( pdf(1*sign1,1,3)*pdf(4*sign2,2,3) + 
     &        pdf(3*sign1,1,3)*pdf(2*sign2,2,3) + 
     &        lbquark*pdf(5*sign1,1,3)*pdf(2*sign2,2,3) + 
     &        lbquark*pdf(5*sign1,1,3)*pdf(4*sign2,2,3) )*polcolgq
         res(17,3) = diag(17) * dipole(4) * 
     &      ( pdf(1*sign1,1,3)*pdf(3*sign2,2,3) + 
     &        lbquark*pdf(1*sign1,1,3)*pdf(5*sign2,2,3) + 
     &        pdf(3*sign1,1,3)*pdf(1*sign2,2,3) + 
     &        lbquark*pdf(3*sign1,1,3)*pdf(5*sign2,2,3) + 
     &        lbquark*pdf(5*sign1,1,3)*pdf(1*sign2,2,3) + 
     &        lbquark*pdf(5*sign1,1,3)*pdf(3*sign2,2,3) )*polcolgq
         res(18,3) = diag(18) * dipole(5) * 
     &      ( pdf(1*sign1,1,3)*pdf(0*sign2,2,3) + 
     &        pdf(3*sign1,1,3)*pdf(0*sign2,2,3) + 
     &        lbquark*pdf(5*sign1,1,3)*pdf(0*sign2,2,3) )*polcolgq
         res(19,3) = diag(19) * dipole(4) * 
     &      ( pdf(1*sign1,1,3)*pdf(-3*sign2,2,3) + 
     &        lbquark*pdf(1*sign1,1,3)*pdf(-5*sign2,2,3) + 
     &        pdf(3*sign1,1,3)*pdf(-1*sign2,2,3) + 
     &        lbquark*pdf(3*sign1,1,3)*pdf(-5*sign2,2,3) + 
     &        lbquark*pdf(5*sign1,1,3)*pdf(-1*sign2,2,3) + 
     &        lbquark*pdf(5*sign1,1,3)*pdf(-3*sign2,2,3) )*polcolgq
         res(20,3) = diag(20) * dipole(4) * 
     &      ( pdf(1*sign1,1,3)*pdf(-4*sign2,2,3) + 
     &        pdf(3*sign1,1,3)*pdf(-2*sign2,2,3) + 
     &        lbquark*pdf(5*sign1,1,3)*pdf(-2*sign2,2,3) + 
     &        lbquark*pdf(5*sign1,1,3)*pdf(-4*sign2,2,3) )*polcolgq
         res(22,3) = diag(22) * dipole(4) * 
     &      ( pdf(1*sign1,1,3)*pdf(-1*sign2,2,3) + 
     &        pdf(3*sign1,1,3)*pdf(-3*sign2,2,3) + 
     &        lbquark*pdf(5*sign1,1,3)*pdf(-5*sign2,2,3) )*polcolgq
         res(25,3) = diag(25) * dipole(6) * 
     &      ( pdf(1*sign1,1,3)*pdf(-1*sign2,2,3) + 
     &        pdf(3*sign1,1,3)*pdf(-3*sign2,2,3) + 
     &        lbquark*pdf(5*sign1,1,3)*pdf(-5*sign2,2,3) )/2d0*polcolqq
         res(26,3) = diag(26) * dipole(4) * 
     &      ( pdf(1*sign1,1,3)*pdf(-2*sign2,2,3) + 
     &        pdf(3*sign1,1,3)*pdf(-4*sign2,2,3) )*polcolgq
         res(27,3) = diag(27) * dipole(7) * 
     &      ( pdf(0*sign1,1,3)*pdf(2*sign2,2,3) + 
     &        pdf(0*sign1,1,3)*pdf(4*sign2,2,3) )*polcolgq
         res(28,3) = diag(28) * dipole(8) * 
     &      ( pdf(0*sign1,1,3)*pdf(1*sign2,2,3) + 
     &        pdf(0*sign1,1,3)*pdf(3*sign2,2,3) + 
     &        lbquark*pdf(0*sign1,1,3)*pdf(5*sign2,2,3) )*polcolgq
         res(29,3) = diag(29) * dipole(9) * 
     &      ( pdf(0*sign1,1,3)*pdf(0*sign2,2,3) + 
     &        pdf(0*sign1,1,3)*pdf(0*sign2,2,3) )*polcolgq
         res(30,3) = diag(30) * dipole(10) * 
     &      ( pdf(0*sign1,1,3)*pdf(0*sign2,2,3) + 
     &        pdf(0*sign1,1,3)*pdf(0*sign2,2,3) + 
     &        lbquark*pdf(0*sign1,1,3)*pdf(0*sign2,2,3) )*polcolgq
         res(31,3) = diag(31) * dipole(11) * 
     &      ( pdf(0*sign1,1,3)*pdf(-1*sign2,2,3) + 
     &        pdf(0*sign1,1,3)*pdf(-3*sign2,2,3) + 
     &        lbquark*pdf(0*sign1,1,3)*pdf(-5*sign2,2,3) )*polcolgq
         res(32,3) = diag(32) * dipole(12) * 
     &      ( pdf(0*sign1,1,3)*pdf(-2*sign2,2,3) + 
     &        pdf(0*sign1,1,3)*pdf(-4*sign2,2,3) )*polcolgq
         res(40,3) = diag(40) * dipole(13) * 
     &      ( pdf(-1*sign1,1,3)*pdf(1*sign2,2,3) + 
     &        pdf(-3*sign1,1,3)*pdf(3*sign2,2,3) + 
     &        lbquark*pdf(-5*sign1,1,3)*pdf(5*sign2,2,3) )/2d0*polcolqq
         res(41,3) = diag(41) * dipole(14) * 
     &      ( pdf(-1*sign1,1,3)*pdf(0*sign2,2,3) + 
     &        pdf(-3*sign1,1,3)*pdf(0*sign2,2,3) + 
     &        lbquark*pdf(-5*sign1,1,3)*pdf(0*sign2,2,3) )*polcolgq
         res(42,3) = diag(42) * dipole(15) * 
     &      ( pdf(-1*sign1,1,3)*pdf(-1*sign2,2,3) + 
     &        pdf(-3*sign1,1,3)*pdf(-3*sign2,2,3) + 
     &        lbquark*pdf(-5*sign1,1,3)*pdf(-5*sign2,2,3) )/2d0*polcolgq
         res(52,3) = diag(52) * dipole(16) * 
     &      ( pdf(-2*sign1,1,3)*pdf(2*sign2,2,3) + 
     &        pdf(-4*sign1,1,3)*pdf(4*sign2,2,3) )/2d0*polcolqq
         res(54,3) = diag(54) * dipole(17) * 
     &      ( pdf(-2*sign1,1,3)*pdf(0*sign2,2,3) + 
     &        pdf(-4*sign1,1,3)*pdf(0*sign2,2,3) )*polcolgq
         res(55,3) = diag(55) * dipole(18) * 
     &      ( pdf(-2*sign1,1,3)*pdf(-1*sign2,2,3) + 
     &        pdf(-4*sign1,1,3)*pdf(-3*sign2,2,3) )*polcolgq
         res(56,3) = diag(56) * dipole(18) * 
     &      ( pdf(-2*sign1,1,3)*pdf(-2*sign2,2,3) + 
     &        pdf(-4*sign1,1,3)*pdf(-4*sign2,2,3) )/2d0*polcolgq
       endif 

       if (lok(4)) then 
         call qqgZAZA_curr(borncurr(0,-1,1,13),4,p(0,3,4),p(0,1,4),p(0,2,4),-1,-1,-1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,14),4,p(0,2,4),p(0,1,4),p(0,3,4),1,-1,1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,15),4,p(0,3,4),p(0,2,4),p(0,1,4),-1,-1,-1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,16),4,p(0,2,4),p(0,3,4),p(0,1,4),1,1,-1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,17),4,p(0,1,4),p(0,2,4),p(0,3,4),1,-1,1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,18),4,p(0,1,4),p(0,3,4),p(0,2,4),1,1,-1,bos)

         dipole(1) = - dipoleVVj_CS147(4,1,p(0,1,1),p(0,4,1),p(0,2,1),xuz(1,4),borncurr(0,-1,1,16))
         dipole(2) = - dipoleVVj_CS147(4,1,p(0,1,1),p(0,4,1),p(0,2,1),xuz(1,4),borncurr(0,-1,2,16))
         dipole(3) = - dipoleVVj_CS145(4,0,p(0,1,1),p(0,4,1),p(0,2,1),xuz(1,4),borncurr(0,0,1,18))
         dipole(4) = - dipoleVVj_CS145(4,1,p(0,1,1),p(0,4,1),p(0,2,1),xuz(1,4),borncurr(0,0,1,17))
         dipole(5) = - dipoleVVj_CS145(4,0,p(0,1,1),p(0,4,1),p(0,2,1),xuz(1,4),borncurr(0,0,2,18))
         dipole(6) = - dipoleVVj_CS145(4,1,p(0,1,1),p(0,4,1),p(0,2,1),xuz(1,4),borncurr(0,0,2,17))
         dipole(7) = - dipoleVVj_CS148(4,1,p(0,1,1),p(0,4,1),p(0,2,1),xuz(1,4),borncurr(0,-1,1,16))
         dipole(8) = - dipoleVVj_CS148(4,1,p(0,1,1),p(0,4,1),p(0,2,1),xuz(1,4),borncurr(0,-1,2,16))
         dipole(9) = - dipoleVVj_CS146(4,0,p(0,1,1),p(0,4,1),p(0,2,1),xuz(1,4),borncurr(0,0,1,18))
         dipole(10) = - dipoleVVj_CS146(4,0,p(0,1,1),p(0,4,1),p(0,2,1),xuz(1,4),borncurr(0,0,2,18))
         dipole(11) = - dipoleVVj_CS148(4,1,p(0,1,1),p(0,4,1),p(0,2,1),xuz(1,4),borncurr(0,-1,2,15))
         dipole(12) = - dipoleVVj_CS148(4,1,p(0,1,1),p(0,4,1),p(0,2,1),xuz(1,4),borncurr(0,-1,1,15))
         dipole(13) = - dipoleVVj_CS145(4,1,p(0,1,1),p(0,4,1),p(0,2,1),xuz(1,4),borncurr(0,0,2,14))
         dipole(14) = - dipoleVVj_CS145(4,0,p(0,1,1),p(0,4,1),p(0,2,1),xuz(1,4),borncurr(0,0,2,13))
         dipole(15) = - dipoleVVj_CS147(4,1,p(0,1,1),p(0,4,1),p(0,2,1),xuz(1,4),borncurr(0,-1,2,15))
         dipole(16) = - dipoleVVj_CS147(4,1,p(0,1,1),p(0,4,1),p(0,2,1),xuz(1,4),borncurr(0,-1,1,15))
         dipole(17) = - dipoleVVj_CS145(4,1,p(0,1,1),p(0,4,1),p(0,2,1),xuz(1,4),borncurr(0,0,1,14))
         dipole(18) = - dipoleVVj_CS145(4,0,p(0,1,1),p(0,4,1),p(0,2,1),xuz(1,4),borncurr(0,0,1,13))

         res(1,4) = diag(1) * dipole(1) * 
     &      ( pdf(2*sign1,1,4)*pdf(2*sign2,2,4) + 
     &        pdf(4*sign1,1,4)*pdf(4*sign2,2,4) )/2d0*polcolgq
         res(2,4) = diag(2) * dipole(2) * 
     &      ( pdf(2*sign1,1,4)*pdf(1*sign2,2,4) + 
     &        pdf(4*sign1,1,4)*pdf(3*sign2,2,4) )*polcolgq
         res(5,4) = diag(5) * dipole(3) * 
     &      ( pdf(2*sign1,1,4)*pdf(0*sign2,2,4) + 
     &        pdf(4*sign1,1,4)*pdf(0*sign2,2,4) )*polcolgq
         res(13,4) = diag(13) * dipole(4) * 
     &      ( pdf(2*sign1,1,4)*pdf(-2*sign2,2,4) + 
     &        pdf(4*sign1,1,4)*pdf(-4*sign2,2,4) )/2d0*polcolqq
         res(15,4) = diag(15) * dipole(2) * 
     &      ( pdf(1*sign1,1,4)*pdf(1*sign2,2,4) + 
     &        pdf(3*sign1,1,4)*pdf(3*sign2,2,4) + 
     &        lbquark*pdf(5*sign1,1,4)*pdf(5*sign2,2,4) )/2d0*polcolgq
         res(18,4) = diag(18) * dipole(5) * 
     &      ( pdf(1*sign1,1,4)*pdf(0*sign2,2,4) + 
     &        pdf(3*sign1,1,4)*pdf(0*sign2,2,4) + 
     &        lbquark*pdf(5*sign1,1,4)*pdf(0*sign2,2,4) )*polcolgq
         res(25,4) = diag(25) * dipole(6) * 
     &      ( pdf(1*sign1,1,4)*pdf(-1*sign2,2,4) + 
     &        pdf(3*sign1,1,4)*pdf(-3*sign2,2,4) + 
     &        lbquark*pdf(5*sign1,1,4)*pdf(-5*sign2,2,4) )/2d0*polcolqq
         res(27,4) = diag(27) * dipole(7) * 
     &      ( pdf(0*sign1,1,4)*pdf(2*sign2,2,4) + 
     &        pdf(0*sign1,1,4)*pdf(4*sign2,2,4) )*polcolgq
         res(28,4) = diag(28) * dipole(8) * 
     &      ( pdf(0*sign1,1,4)*pdf(1*sign2,2,4) + 
     &        pdf(0*sign1,1,4)*pdf(3*sign2,2,4) + 
     &        lbquark*pdf(0*sign1,1,4)*pdf(5*sign2,2,4) )*polcolgq
         res(29,4) = diag(29) * dipole(9) * 
     &      ( pdf(0*sign1,1,4)*pdf(0*sign2,2,4) + 
     &        pdf(0*sign1,1,4)*pdf(0*sign2,2,4) )*polcolgq
         res(30,4) = diag(30) * dipole(10) * 
     &      ( pdf(0*sign1,1,4)*pdf(0*sign2,2,4) + 
     &        pdf(0*sign1,1,4)*pdf(0*sign2,2,4) + 
     &        lbquark*pdf(0*sign1,1,4)*pdf(0*sign2,2,4) )*polcolgq
         res(31,4) = diag(31) * dipole(11) * 
     &      ( pdf(0*sign1,1,4)*pdf(-1*sign2,2,4) + 
     &        pdf(0*sign1,1,4)*pdf(-3*sign2,2,4) + 
     &        lbquark*pdf(0*sign1,1,4)*pdf(-5*sign2,2,4) )*polcolgq
         res(32,4) = diag(32) * dipole(12) * 
     &      ( pdf(0*sign1,1,4)*pdf(-2*sign2,2,4) + 
     &        pdf(0*sign1,1,4)*pdf(-4*sign2,2,4) )*polcolgq
         res(33,4) = diag(33) * dipole(1) * 
     &      ( pdf(-1*sign1,1,4)*pdf(4*sign2,2,4) + 
     &        pdf(-3*sign1,1,4)*pdf(2*sign2,2,4) + 
     &        lbquark*pdf(-5*sign1,1,4)*pdf(2*sign2,2,4) + 
     &        lbquark*pdf(-5*sign1,1,4)*pdf(4*sign2,2,4) )*polcolgq
         res(34,4) = diag(34) * dipole(2) * 
     &      ( pdf(-1*sign1,1,4)*pdf(3*sign2,2,4) + 
     &        lbquark*pdf(-1*sign1,1,4)*pdf(5*sign2,2,4) + 
     &        pdf(-3*sign1,1,4)*pdf(1*sign2,2,4) + 
     &        lbquark*pdf(-3*sign1,1,4)*pdf(5*sign2,2,4) + 
     &        lbquark*pdf(-5*sign1,1,4)*pdf(1*sign2,2,4) + 
     &        lbquark*pdf(-5*sign1,1,4)*pdf(3*sign2,2,4) )*polcolgq
         res(35,4) = diag(35) * dipole(1) * 
     &      ( pdf(-1*sign1,1,4)*pdf(2*sign2,2,4) + 
     &        pdf(-3*sign1,1,4)*pdf(4*sign2,2,4) )*polcolgq
         res(39,4) = diag(39) * dipole(2) * 
     &      ( pdf(-1*sign1,1,4)*pdf(1*sign2,2,4) + 
     &        pdf(-3*sign1,1,4)*pdf(3*sign2,2,4) + 
     &        lbquark*pdf(-5*sign1,1,4)*pdf(5*sign2,2,4) )*polcolgq
         res(40,4) = diag(40) * dipole(13) * 
     &      ( pdf(-1*sign1,1,4)*pdf(1*sign2,2,4) + 
     &        pdf(-3*sign1,1,4)*pdf(3*sign2,2,4) + 
     &        lbquark*pdf(-5*sign1,1,4)*pdf(5*sign2,2,4) )/2d0*polcolqq
         res(41,4) = diag(41) * dipole(14) * 
     &      ( pdf(-1*sign1,1,4)*pdf(0*sign2,2,4) + 
     &        pdf(-3*sign1,1,4)*pdf(0*sign2,2,4) + 
     &        lbquark*pdf(-5*sign1,1,4)*pdf(0*sign2,2,4) )*polcolgq
         res(42,4) = diag(42) * dipole(15) * 
     &      ( pdf(-1*sign1,1,4)*pdf(-1*sign2,2,4) + 
     &        pdf(-3*sign1,1,4)*pdf(-3*sign2,2,4) + 
     &        lbquark*pdf(-5*sign1,1,4)*pdf(-5*sign2,2,4) )/2d0*polcolgq
         res(43,4) = diag(43) * dipole(16) * 
     &      ( pdf(-1*sign1,1,4)*pdf(-2*sign2,2,4) + 
     &        pdf(-3*sign1,1,4)*pdf(-4*sign2,2,4) )*polcolgq
         res(44,4) = diag(44) * dipole(15) * 
     &      ( pdf(-1*sign1,1,4)*pdf(-3*sign2,2,4) + 
     &        lbquark*pdf(-1*sign1,1,4)*pdf(-5*sign2,2,4) + 
     &        pdf(-3*sign1,1,4)*pdf(-1*sign2,2,4) + 
     &        lbquark*pdf(-3*sign1,1,4)*pdf(-5*sign2,2,4) + 
     &        lbquark*pdf(-5*sign1,1,4)*pdf(-1*sign2,2,4) + 
     &        lbquark*pdf(-5*sign1,1,4)*pdf(-3*sign2,2,4) )*polcolgq
         res(45,4) = diag(45) * dipole(16) * 
     &      ( pdf(-1*sign1,1,4)*pdf(-4*sign2,2,4) + 
     &        pdf(-3*sign1,1,4)*pdf(-2*sign2,2,4) + 
     &        lbquark*pdf(-5*sign1,1,4)*pdf(-2*sign2,2,4) + 
     &        lbquark*pdf(-5*sign1,1,4)*pdf(-4*sign2,2,4) )*polcolgq
         res(46,4) = diag(46) * dipole(1) * 
     &      ( pdf(-2*sign1,1,4)*pdf(4*sign2,2,4) + 
     &        pdf(-4*sign1,1,4)*pdf(2*sign2,2,4) )*polcolgq
         res(47,4) = diag(47) * dipole(2) * 
     &      ( pdf(-2*sign1,1,4)*pdf(3*sign2,2,4) + 
     &        lbquark*pdf(-2*sign1,1,4)*pdf(5*sign2,2,4) + 
     &        pdf(-4*sign1,1,4)*pdf(1*sign2,2,4) + 
     &        lbquark*pdf(-4*sign1,1,4)*pdf(5*sign2,2,4) )*polcolgq
         res(50,4) = diag(50) * dipole(1) * 
     &      ( pdf(-2*sign1,1,4)*pdf(2*sign2,2,4) + 
     &        pdf(-4*sign1,1,4)*pdf(4*sign2,2,4) )*polcolgq
         res(52,4) = diag(52) * dipole(17) * 
     &      ( pdf(-2*sign1,1,4)*pdf(2*sign2,2,4) + 
     &        pdf(-4*sign1,1,4)*pdf(4*sign2,2,4) )/2d0*polcolqq
         res(53,4) = diag(53) * dipole(2) * 
     &      ( pdf(-2*sign1,1,4)*pdf(1*sign2,2,4) + 
     &        pdf(-4*sign1,1,4)*pdf(3*sign2,2,4) )*polcolgq
         res(54,4) = diag(54) * dipole(18) * 
     &      ( pdf(-2*sign1,1,4)*pdf(0*sign2,2,4) + 
     &        pdf(-4*sign1,1,4)*pdf(0*sign2,2,4) )*polcolgq
         res(56,4) = diag(56) * dipole(16) * 
     &      ( pdf(-2*sign1,1,4)*pdf(-2*sign2,2,4) + 
     &        pdf(-4*sign1,1,4)*pdf(-4*sign2,2,4) )/2d0*polcolgq
         res(57,4) = diag(57) * dipole(15) * 
     &      ( pdf(-2*sign1,1,4)*pdf(-3*sign2,2,4) + 
     &        lbquark*pdf(-2*sign1,1,4)*pdf(-5*sign2,2,4) + 
     &        pdf(-4*sign1,1,4)*pdf(-1*sign2,2,4) + 
     &        lbquark*pdf(-4*sign1,1,4)*pdf(-5*sign2,2,4) )*polcolgq
         res(58,4) = diag(58) * dipole(16) * 
     &      ( pdf(-2*sign1,1,4)*pdf(-4*sign2,2,4) + 
     &        pdf(-4*sign1,1,4)*pdf(-2*sign2,2,4) )*polcolgq
       endif 

       if (lok(5)) then 
         call qqgZAZA_curr(borncurr(0,-1,1,19),5,p(0,3,5),p(0,1,5),p(0,2,5),-1,-1,-1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,20),5,p(0,2,5),p(0,1,5),p(0,3,5),1,-1,1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,21),5,p(0,3,5),p(0,2,5),p(0,1,5),-1,-1,-1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,22),5,p(0,2,5),p(0,3,5),p(0,1,5),1,1,-1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,23),5,p(0,1,5),p(0,2,5),p(0,3,5),1,-1,1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,24),5,p(0,1,5),p(0,3,5),p(0,2,5),1,1,-1,bos)

         dipole(1) = - dipoleVVj_CS147(5,1,p(0,2,1),p(0,4,1),p(0,1,1),xuz(1,5),borncurr(0,-1,1,24))
         dipole(2) = - dipoleVVj_CS148(5,1,p(0,2,1),p(0,4,1),p(0,1,1),xuz(1,5),borncurr(0,-1,1,24))
         dipole(3) = - dipoleVVj_CS145(5,1,p(0,2,1),p(0,4,1),p(0,1,1),xuz(1,5),borncurr(0,0,1,23))
         dipole(4) = - dipoleVVj_CS147(5,1,p(0,2,1),p(0,4,1),p(0,1,1),xuz(1,5),borncurr(0,-1,2,24))
         dipole(5) = - dipoleVVj_CS148(5,1,p(0,2,1),p(0,4,1),p(0,1,1),xuz(1,5),borncurr(0,-1,2,24))
         dipole(6) = - dipoleVVj_CS145(5,1,p(0,2,1),p(0,4,1),p(0,1,1),xuz(1,5),borncurr(0,0,2,23))
         dipole(7) = - dipoleVVj_CS145(5,0,p(0,2,1),p(0,4,1),p(0,1,1),xuz(1,5),borncurr(0,0,1,22))
         dipole(8) = - dipoleVVj_CS145(5,0,p(0,2,1),p(0,4,1),p(0,1,1),xuz(1,5),borncurr(0,0,2,22))
         dipole(9) = - dipoleVVj_CS146(5,0,p(0,2,1),p(0,4,1),p(0,1,1),xuz(1,5),borncurr(0,0,1,22))
         dipole(10) = - dipoleVVj_CS146(5,0,p(0,2,1),p(0,4,1),p(0,1,1),xuz(1,5),borncurr(0,0,2,22))
         dipole(11) = - dipoleVVj_CS145(5,0,p(0,2,1),p(0,4,1),p(0,1,1),xuz(1,5),borncurr(0,0,2,21))
         dipole(12) = - dipoleVVj_CS145(5,0,p(0,2,1),p(0,4,1),p(0,1,1),xuz(1,5),borncurr(0,0,1,21))
         dipole(13) = - dipoleVVj_CS145(5,1,p(0,2,1),p(0,4,1),p(0,1,1),xuz(1,5),borncurr(0,0,2,20))
         dipole(14) = - dipoleVVj_CS148(5,1,p(0,2,1),p(0,4,1),p(0,1,1),xuz(1,5),borncurr(0,-1,2,19))
         dipole(15) = - dipoleVVj_CS147(5,1,p(0,2,1),p(0,4,1),p(0,1,1),xuz(1,5),borncurr(0,-1,2,19))
         dipole(16) = - dipoleVVj_CS145(5,1,p(0,2,1),p(0,4,1),p(0,1,1),xuz(1,5),borncurr(0,0,1,20))
         dipole(17) = - dipoleVVj_CS148(5,1,p(0,2,1),p(0,4,1),p(0,1,1),xuz(1,5),borncurr(0,-1,1,19))
         dipole(18) = - dipoleVVj_CS147(5,1,p(0,2,1),p(0,4,1),p(0,1,1),xuz(1,5),borncurr(0,-1,1,19))

         res(1,5) = diag(1) * dipole(1) * 
     &      ( pdf(2*sign1,1,5)*pdf(2*sign2,2,5) + 
     &        pdf(4*sign1,1,5)*pdf(4*sign2,2,5) )/2d0*polcolgq
         res(3,5) = diag(3) * dipole(1) * 
     &      ( pdf(2*sign1,1,5)*pdf(4*sign2,2,5) + 
     &        pdf(4*sign1,1,5)*pdf(2*sign2,2,5) )*polcolgq
         res(4,5) = diag(4) * dipole(1) * 
     &      ( pdf(2*sign1,1,5)*pdf(3*sign2,2,5) + 
     &        lbquark*pdf(2*sign1,1,5)*pdf(5*sign2,2,5) + 
     &        pdf(4*sign1,1,5)*pdf(1*sign2,2,5) + 
     &        lbquark*pdf(4*sign1,1,5)*pdf(5*sign2,2,5) )*polcolgq
         res(5,5) = diag(5) * dipole(2) * 
     &      ( pdf(2*sign1,1,5)*pdf(0*sign2,2,5) + 
     &        pdf(4*sign1,1,5)*pdf(0*sign2,2,5) )*polcolgq
         res(6,5) = diag(6) * dipole(1) * 
     &      ( pdf(2*sign1,1,5)*pdf(-3*sign2,2,5) + 
     &        lbquark*pdf(2*sign1,1,5)*pdf(-5*sign2,2,5) + 
     &        pdf(4*sign1,1,5)*pdf(-1*sign2,2,5) + 
     &        lbquark*pdf(4*sign1,1,5)*pdf(-5*sign2,2,5) )*polcolgq
         res(7,5) = diag(7) * dipole(1) * 
     &      ( pdf(2*sign1,1,5)*pdf(-4*sign2,2,5) + 
     &        pdf(4*sign1,1,5)*pdf(-2*sign2,2,5) )*polcolgq
         res(8,5) = diag(8) * dipole(1) * 
     &      ( pdf(2*sign1,1,5)*pdf(-1*sign2,2,5) + 
     &        pdf(4*sign1,1,5)*pdf(-3*sign2,2,5) )*polcolgq
         res(9,5) = diag(9) * dipole(1) * 
     &      ( pdf(2*sign1,1,5)*pdf(-2*sign2,2,5) + 
     &        pdf(4*sign1,1,5)*pdf(-4*sign2,2,5) )*polcolgq
         res(13,5) = diag(13) * dipole(3) * 
     &      ( pdf(2*sign1,1,5)*pdf(-2*sign2,2,5) + 
     &        pdf(4*sign1,1,5)*pdf(-4*sign2,2,5) )/2d0*polcolqq
         res(14,5) = diag(14) * dipole(4) * 
     &      ( pdf(1*sign1,1,5)*pdf(2*sign2,2,5) + 
     &        pdf(3*sign1,1,5)*pdf(4*sign2,2,5) )*polcolgq
         res(15,5) = diag(15) * dipole(4) * 
     &      ( pdf(1*sign1,1,5)*pdf(1*sign2,2,5) + 
     &        pdf(3*sign1,1,5)*pdf(3*sign2,2,5) + 
     &        lbquark*pdf(5*sign1,1,5)*pdf(5*sign2,2,5) )/2d0*polcolgq
         res(16,5) = diag(16) * dipole(4) * 
     &      ( pdf(1*sign1,1,5)*pdf(4*sign2,2,5) + 
     &        pdf(3*sign1,1,5)*pdf(2*sign2,2,5) + 
     &        lbquark*pdf(5*sign1,1,5)*pdf(2*sign2,2,5) + 
     &        lbquark*pdf(5*sign1,1,5)*pdf(4*sign2,2,5) )*polcolgq
         res(17,5) = diag(17) * dipole(4) * 
     &      ( pdf(1*sign1,1,5)*pdf(3*sign2,2,5) + 
     &        lbquark*pdf(1*sign1,1,5)*pdf(5*sign2,2,5) + 
     &        pdf(3*sign1,1,5)*pdf(1*sign2,2,5) + 
     &        lbquark*pdf(3*sign1,1,5)*pdf(5*sign2,2,5) + 
     &        lbquark*pdf(5*sign1,1,5)*pdf(1*sign2,2,5) + 
     &        lbquark*pdf(5*sign1,1,5)*pdf(3*sign2,2,5) )*polcolgq
         res(18,5) = diag(18) * dipole(5) * 
     &      ( pdf(1*sign1,1,5)*pdf(0*sign2,2,5) + 
     &        pdf(3*sign1,1,5)*pdf(0*sign2,2,5) + 
     &        lbquark*pdf(5*sign1,1,5)*pdf(0*sign2,2,5) )*polcolgq
         res(19,5) = diag(19) * dipole(4) * 
     &      ( pdf(1*sign1,1,5)*pdf(-3*sign2,2,5) + 
     &        lbquark*pdf(1*sign1,1,5)*pdf(-5*sign2,2,5) + 
     &        pdf(3*sign1,1,5)*pdf(-1*sign2,2,5) + 
     &        lbquark*pdf(3*sign1,1,5)*pdf(-5*sign2,2,5) + 
     &        lbquark*pdf(5*sign1,1,5)*pdf(-1*sign2,2,5) + 
     &        lbquark*pdf(5*sign1,1,5)*pdf(-3*sign2,2,5) )*polcolgq
         res(20,5) = diag(20) * dipole(4) * 
     &      ( pdf(1*sign1,1,5)*pdf(-4*sign2,2,5) + 
     &        pdf(3*sign1,1,5)*pdf(-2*sign2,2,5) + 
     &        lbquark*pdf(5*sign1,1,5)*pdf(-2*sign2,2,5) + 
     &        lbquark*pdf(5*sign1,1,5)*pdf(-4*sign2,2,5) )*polcolgq
         res(22,5) = diag(22) * dipole(4) * 
     &      ( pdf(1*sign1,1,5)*pdf(-1*sign2,2,5) + 
     &        pdf(3*sign1,1,5)*pdf(-3*sign2,2,5) + 
     &        lbquark*pdf(5*sign1,1,5)*pdf(-5*sign2,2,5) )*polcolgq
         res(25,5) = diag(25) * dipole(6) * 
     &      ( pdf(1*sign1,1,5)*pdf(-1*sign2,2,5) + 
     &        pdf(3*sign1,1,5)*pdf(-3*sign2,2,5) + 
     &        lbquark*pdf(5*sign1,1,5)*pdf(-5*sign2,2,5) )/2d0*polcolqq
         res(26,5) = diag(26) * dipole(4) * 
     &      ( pdf(1*sign1,1,5)*pdf(-2*sign2,2,5) + 
     &        pdf(3*sign1,1,5)*pdf(-4*sign2,2,5) )*polcolgq
         res(27,5) = diag(27) * dipole(7) * 
     &      ( pdf(0*sign1,1,5)*pdf(2*sign2,2,5) + 
     &        pdf(0*sign1,1,5)*pdf(4*sign2,2,5) )*polcolgq
         res(28,5) = diag(28) * dipole(8) * 
     &      ( pdf(0*sign1,1,5)*pdf(1*sign2,2,5) + 
     &        pdf(0*sign1,1,5)*pdf(3*sign2,2,5) + 
     &        lbquark*pdf(0*sign1,1,5)*pdf(5*sign2,2,5) )*polcolgq
         res(29,5) = diag(29) * dipole(9) * 
     &      ( pdf(0*sign1,1,5)*pdf(0*sign2,2,5) + 
     &        pdf(0*sign1,1,5)*pdf(0*sign2,2,5) )*polcolgq
         res(30,5) = diag(30) * dipole(10) * 
     &      ( pdf(0*sign1,1,5)*pdf(0*sign2,2,5) + 
     &        pdf(0*sign1,1,5)*pdf(0*sign2,2,5) + 
     &        lbquark*pdf(0*sign1,1,5)*pdf(0*sign2,2,5) )*polcolgq
         res(31,5) = diag(31) * dipole(11) * 
     &      ( pdf(0*sign1,1,5)*pdf(-1*sign2,2,5) + 
     &        pdf(0*sign1,1,5)*pdf(-3*sign2,2,5) + 
     &        lbquark*pdf(0*sign1,1,5)*pdf(-5*sign2,2,5) )*polcolgq
         res(32,5) = diag(32) * dipole(12) * 
     &      ( pdf(0*sign1,1,5)*pdf(-2*sign2,2,5) + 
     &        pdf(0*sign1,1,5)*pdf(-4*sign2,2,5) )*polcolgq
         res(40,5) = diag(40) * dipole(13) * 
     &      ( pdf(-1*sign1,1,5)*pdf(1*sign2,2,5) + 
     &        pdf(-3*sign1,1,5)*pdf(3*sign2,2,5) + 
     &        lbquark*pdf(-5*sign1,1,5)*pdf(5*sign2,2,5) )/2d0*polcolqq
         res(41,5) = diag(41) * dipole(14) * 
     &      ( pdf(-1*sign1,1,5)*pdf(0*sign2,2,5) + 
     &        pdf(-3*sign1,1,5)*pdf(0*sign2,2,5) + 
     &        lbquark*pdf(-5*sign1,1,5)*pdf(0*sign2,2,5) )*polcolgq
         res(42,5) = diag(42) * dipole(15) * 
     &      ( pdf(-1*sign1,1,5)*pdf(-1*sign2,2,5) + 
     &        pdf(-3*sign1,1,5)*pdf(-3*sign2,2,5) + 
     &        lbquark*pdf(-5*sign1,1,5)*pdf(-5*sign2,2,5) )/2d0*polcolgq
         res(52,5) = diag(52) * dipole(16) * 
     &      ( pdf(-2*sign1,1,5)*pdf(2*sign2,2,5) + 
     &        pdf(-4*sign1,1,5)*pdf(4*sign2,2,5) )/2d0*polcolqq
         res(54,5) = diag(54) * dipole(17) * 
     &      ( pdf(-2*sign1,1,5)*pdf(0*sign2,2,5) + 
     &        pdf(-4*sign1,1,5)*pdf(0*sign2,2,5) )*polcolgq
         res(55,5) = diag(55) * dipole(18) * 
     &      ( pdf(-2*sign1,1,5)*pdf(-1*sign2,2,5) + 
     &        pdf(-4*sign1,1,5)*pdf(-3*sign2,2,5) )*polcolgq
         res(56,5) = diag(56) * dipole(18) * 
     &      ( pdf(-2*sign1,1,5)*pdf(-2*sign2,2,5) + 
     &        pdf(-4*sign1,1,5)*pdf(-4*sign2,2,5) )/2d0*polcolgq
       endif 

       if (lok(6)) then 
         call qqgZAZA_curr(borncurr(0,-1,1,25),6,p(0,3,6),p(0,1,6),p(0,2,6),-1,-1,-1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,26),6,p(0,2,6),p(0,1,6),p(0,3,6),1,-1,1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,27),6,p(0,3,6),p(0,2,6),p(0,1,6),-1,-1,-1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,28),6,p(0,2,6),p(0,3,6),p(0,1,6),1,1,-1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,29),6,p(0,1,6),p(0,2,6),p(0,3,6),1,-1,1,bos)

         dipole(1) = - dipoleVVj_CS67(6,1,p(0,1,1),p(0,3,1),p(0,4,1),xuz(1,6),borncurr(0,-1,1,28))
         dipole(2) = - dipoleVVj_CS67(6,1,p(0,1,1),p(0,3,1),p(0,4,1),xuz(1,6),borncurr(0,-1,2,28))
         dipole(3) = - dipoleVVj_CS67(6,1,p(0,1,1),p(0,3,1),p(0,4,1),xuz(1,6),borncurr(0,-1,2,27))
         dipole(4) = - dipoleVVj_CS67(6,1,p(0,1,1),p(0,3,1),p(0,4,1),xuz(1,6),borncurr(0,-1,1,27))
         dipole(5) = - dipoleVVj_CS65(6,0,p(0,1,1),p(0,3,1),p(0,4,1),xuz(1,6),borncurr(0,0,1,29))
         dipole(6) = - dipoleVVj_CS65(6,0,p(0,1,1),p(0,3,1),p(0,4,1),xuz(1,6),borncurr(0,0,2,29))
         dipole(7) = - dipoleVVj_CS66(6,0,p(0,1,1),p(0,3,1),p(0,4,1),xuz(1,6),borncurr(0,0,1,26))
         dipole(8) = - dipoleVVj_CS66(6,0,p(0,1,1),p(0,3,1),p(0,4,1),xuz(1,6),borncurr(0,0,2,26))
         dipole(9) = - dipoleVVj_CS66(6,1,p(0,1,1),p(0,3,1),p(0,4,1),xuz(1,6),borncurr(0,0,1,25))
         dipole(10) = - dipoleVVj_CS66(6,1,p(0,1,1),p(0,3,1),p(0,4,1),xuz(1,6),borncurr(0,0,2,25))
         dipole(11) = - dipoleVVj_CS66(6,0,p(0,1,1),p(0,3,1),p(0,4,1),xuz(1,6),borncurr(0,0,2,29))
         dipole(12) = - dipoleVVj_CS66(6,0,p(0,1,1),p(0,3,1),p(0,4,1),xuz(1,6),borncurr(0,0,1,29))
         dipole(13) = - dipoleVVj_CS65(6,0,p(0,1,1),p(0,3,1),p(0,4,1),xuz(1,6),borncurr(0,0,2,26))
         dipole(14) = - dipoleVVj_CS65(6,0,p(0,1,1),p(0,3,1),p(0,4,1),xuz(1,6),borncurr(0,0,1,26))

         res(1,6) = diag(1) * dipole(1) * 
     &      ( pdf(2*sign1,1,6)*pdf(2*sign2,2,6) + 
     &        pdf(4*sign1,1,6)*pdf(4*sign2,2,6) )/2d0*polcolgq
         res(3,6) = diag(3) * dipole(1) * 
     &      ( pdf(2*sign1,1,6)*pdf(4*sign2,2,6) + 
     &        pdf(4*sign1,1,6)*pdf(2*sign2,2,6) )*polcolgq
         res(4,6) = diag(4) * dipole(2) * 
     &      ( pdf(2*sign1,1,6)*pdf(3*sign2,2,6) + 
     &        lbquark*pdf(2*sign1,1,6)*pdf(5*sign2,2,6) + 
     &        pdf(4*sign1,1,6)*pdf(1*sign2,2,6) + 
     &        lbquark*pdf(4*sign1,1,6)*pdf(5*sign2,2,6) )*polcolgq
         res(6,6) = diag(6) * dipole(3) * 
     &      ( pdf(2*sign1,1,6)*pdf(-3*sign2,2,6) + 
     &        lbquark*pdf(2*sign1,1,6)*pdf(-5*sign2,2,6) + 
     &        pdf(4*sign1,1,6)*pdf(-1*sign2,2,6) + 
     &        lbquark*pdf(4*sign1,1,6)*pdf(-5*sign2,2,6) )*polcolgq
         res(7,6) = diag(7) * dipole(4) * 
     &      ( pdf(2*sign1,1,6)*pdf(-4*sign2,2,6) + 
     &        pdf(4*sign1,1,6)*pdf(-2*sign2,2,6) )*polcolgq
         res(8,6) = diag(8) * dipole(3) * 
     &      ( pdf(2*sign1,1,6)*pdf(-1*sign2,2,6) + 
     &        pdf(4*sign1,1,6)*pdf(-3*sign2,2,6) )*polcolgq
         res(9,6) = diag(9) * dipole(4) * 
     &      ( pdf(2*sign1,1,6)*pdf(-2*sign2,2,6) + 
     &        pdf(4*sign1,1,6)*pdf(-4*sign2,2,6) )*polcolgq
         res(13,6) = diag(13) * dipole(5) * 
     &      ( pdf(2*sign1,1,6)*pdf(-2*sign2,2,6) + 
     &        pdf(4*sign1,1,6)*pdf(-4*sign2,2,6) )/2d0*polcolqq
         res(14,6) = diag(14) * dipole(1) * 
     &      ( pdf(1*sign1,1,6)*pdf(2*sign2,2,6) + 
     &        pdf(3*sign1,1,6)*pdf(4*sign2,2,6) )*polcolgq
         res(15,6) = diag(15) * dipole(2) * 
     &      ( pdf(1*sign1,1,6)*pdf(1*sign2,2,6) + 
     &        pdf(3*sign1,1,6)*pdf(3*sign2,2,6) + 
     &        lbquark*pdf(5*sign1,1,6)*pdf(5*sign2,2,6) )/2d0*polcolgq
         res(16,6) = diag(16) * dipole(1) * 
     &      ( pdf(1*sign1,1,6)*pdf(4*sign2,2,6) + 
     &        pdf(3*sign1,1,6)*pdf(2*sign2,2,6) + 
     &        lbquark*pdf(5*sign1,1,6)*pdf(2*sign2,2,6) + 
     &        lbquark*pdf(5*sign1,1,6)*pdf(4*sign2,2,6) )*polcolgq
         res(17,6) = diag(17) * dipole(2) * 
     &      ( pdf(1*sign1,1,6)*pdf(3*sign2,2,6) + 
     &        lbquark*pdf(1*sign1,1,6)*pdf(5*sign2,2,6) + 
     &        pdf(3*sign1,1,6)*pdf(1*sign2,2,6) + 
     &        lbquark*pdf(3*sign1,1,6)*pdf(5*sign2,2,6) + 
     &        lbquark*pdf(5*sign1,1,6)*pdf(1*sign2,2,6) + 
     &        lbquark*pdf(5*sign1,1,6)*pdf(3*sign2,2,6) )*polcolgq
         res(19,6) = diag(19) * dipole(3) * 
     &      ( pdf(1*sign1,1,6)*pdf(-3*sign2,2,6) + 
     &        lbquark*pdf(1*sign1,1,6)*pdf(-5*sign2,2,6) + 
     &        pdf(3*sign1,1,6)*pdf(-1*sign2,2,6) + 
     &        lbquark*pdf(3*sign1,1,6)*pdf(-5*sign2,2,6) + 
     &        lbquark*pdf(5*sign1,1,6)*pdf(-1*sign2,2,6) + 
     &        lbquark*pdf(5*sign1,1,6)*pdf(-3*sign2,2,6) )*polcolgq
         res(20,6) = diag(20) * dipole(4) * 
     &      ( pdf(1*sign1,1,6)*pdf(-4*sign2,2,6) + 
     &        pdf(3*sign1,1,6)*pdf(-2*sign2,2,6) + 
     &        lbquark*pdf(5*sign1,1,6)*pdf(-2*sign2,2,6) + 
     &        lbquark*pdf(5*sign1,1,6)*pdf(-4*sign2,2,6) )*polcolgq
         res(22,6) = diag(22) * dipole(3) * 
     &      ( pdf(1*sign1,1,6)*pdf(-1*sign2,2,6) + 
     &        pdf(3*sign1,1,6)*pdf(-3*sign2,2,6) + 
     &        lbquark*pdf(5*sign1,1,6)*pdf(-5*sign2,2,6) )*polcolgq
         res(25,6) = diag(25) * dipole(6) * 
     &      ( pdf(1*sign1,1,6)*pdf(-1*sign2,2,6) + 
     &        pdf(3*sign1,1,6)*pdf(-3*sign2,2,6) + 
     &        lbquark*pdf(5*sign1,1,6)*pdf(-5*sign2,2,6) )/2d0*polcolqq
         res(26,6) = diag(26) * dipole(4) * 
     &      ( pdf(1*sign1,1,6)*pdf(-2*sign2,2,6) + 
     &        pdf(3*sign1,1,6)*pdf(-4*sign2,2,6) )*polcolgq
         res(27,6) = diag(27) * dipole(7) * 
     &      ( pdf(0*sign1,1,6)*pdf(2*sign2,2,6) + 
     &        pdf(0*sign1,1,6)*pdf(4*sign2,2,6) )*polcolqq
         res(28,6) = diag(28) * dipole(8) * 
     &      ( pdf(0*sign1,1,6)*pdf(1*sign2,2,6) + 
     &        pdf(0*sign1,1,6)*pdf(3*sign2,2,6) + 
     &        lbquark*pdf(0*sign1,1,6)*pdf(5*sign2,2,6) )*polcolqq
         res(29,6) = diag(29) * dipole(9) * 
     &      ( pdf(0*sign1,1,6)*pdf(0*sign2,2,6) + 
     &        pdf(0*sign1,1,6)*pdf(0*sign2,2,6) )*polcolgq
         res(30,6) = diag(30) * dipole(10) * 
     &      ( pdf(0*sign1,1,6)*pdf(0*sign2,2,6) + 
     &        pdf(0*sign1,1,6)*pdf(0*sign2,2,6) + 
     &        lbquark*pdf(0*sign1,1,6)*pdf(0*sign2,2,6) )*polcolgq
         res(31,6) = diag(31) * dipole(11) * 
     &      ( pdf(0*sign1,1,6)*pdf(-1*sign2,2,6) + 
     &        pdf(0*sign1,1,6)*pdf(-3*sign2,2,6) + 
     &        lbquark*pdf(0*sign1,1,6)*pdf(-5*sign2,2,6) )*polcolqq
         res(32,6) = diag(32) * dipole(12) * 
     &      ( pdf(0*sign1,1,6)*pdf(-2*sign2,2,6) + 
     &        pdf(0*sign1,1,6)*pdf(-4*sign2,2,6) )*polcolqq
         res(40,6) = diag(40) * dipole(13) * 
     &      ( pdf(-1*sign1,1,6)*pdf(1*sign2,2,6) + 
     &        pdf(-3*sign1,1,6)*pdf(3*sign2,2,6) + 
     &        lbquark*pdf(-5*sign1,1,6)*pdf(5*sign2,2,6) )/2d0*polcolqq
         res(42,6) = diag(42) * dipole(3) * 
     &      ( pdf(-1*sign1,1,6)*pdf(-1*sign2,2,6) + 
     &        pdf(-3*sign1,1,6)*pdf(-3*sign2,2,6) + 
     &        lbquark*pdf(-5*sign1,1,6)*pdf(-5*sign2,2,6) )/2d0*polcolgq
         res(52,6) = diag(52) * dipole(14) * 
     &      ( pdf(-2*sign1,1,6)*pdf(2*sign2,2,6) + 
     &        pdf(-4*sign1,1,6)*pdf(4*sign2,2,6) )/2d0*polcolqq
         res(55,6) = diag(55) * dipole(3) * 
     &      ( pdf(-2*sign1,1,6)*pdf(-1*sign2,2,6) + 
     &        pdf(-4*sign1,1,6)*pdf(-3*sign2,2,6) )*polcolgq
         res(56,6) = diag(56) * dipole(4) * 
     &      ( pdf(-2*sign1,1,6)*pdf(-2*sign2,2,6) + 
     &        pdf(-4*sign1,1,6)*pdf(-4*sign2,2,6) )/2d0*polcolgq
       endif 

       if (lok(7)) then 
         call qqgZAZA_curr(borncurr(0,-1,1,30),7,p(0,3,7),p(0,1,7),p(0,2,7),-1,-1,-1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,31),7,p(0,2,7),p(0,1,7),p(0,3,7),1,-1,1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,32),7,p(0,3,7),p(0,2,7),p(0,1,7),-1,-1,-1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,33),7,p(0,1,7),p(0,2,7),p(0,3,7),1,-1,1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,34),7,p(0,1,7),p(0,3,7),p(0,2,7),1,1,-1,bos)

         dipole(1) = - dipoleVVj_CS67(7,1,p(0,2,1),p(0,3,1),p(0,4,1),xuz(1,7),borncurr(0,-1,1,34))
         dipole(2) = - dipoleVVj_CS66(7,0,p(0,2,1),p(0,3,1),p(0,4,1),xuz(1,7),borncurr(0,0,1,33))
         dipole(3) = - dipoleVVj_CS65(7,0,p(0,2,1),p(0,3,1),p(0,4,1),xuz(1,7),borncurr(0,0,1,33))
         dipole(4) = - dipoleVVj_CS67(7,1,p(0,2,1),p(0,3,1),p(0,4,1),xuz(1,7),borncurr(0,-1,2,34))
         dipole(5) = - dipoleVVj_CS66(7,0,p(0,2,1),p(0,3,1),p(0,4,1),xuz(1,7),borncurr(0,0,2,33))
         dipole(6) = - dipoleVVj_CS65(7,0,p(0,2,1),p(0,3,1),p(0,4,1),xuz(1,7),borncurr(0,0,2,33))
         dipole(7) = - dipoleVVj_CS66(7,1,p(0,2,1),p(0,3,1),p(0,4,1),xuz(1,7),borncurr(0,0,1,32))
         dipole(8) = - dipoleVVj_CS66(7,1,p(0,2,1),p(0,3,1),p(0,4,1),xuz(1,7),borncurr(0,0,2,32))
         dipole(9) = - dipoleVVj_CS67(7,1,p(0,2,1),p(0,3,1),p(0,4,1),xuz(1,7),borncurr(0,-1,2,30))
         dipole(10) = - dipoleVVj_CS65(7,0,p(0,2,1),p(0,3,1),p(0,4,1),xuz(1,7),borncurr(0,0,2,31))
         dipole(11) = - dipoleVVj_CS66(7,0,p(0,2,1),p(0,3,1),p(0,4,1),xuz(1,7),borncurr(0,0,2,31))
         dipole(12) = - dipoleVVj_CS67(7,1,p(0,2,1),p(0,3,1),p(0,4,1),xuz(1,7),borncurr(0,-1,1,30))
         dipole(13) = - dipoleVVj_CS65(7,0,p(0,2,1),p(0,3,1),p(0,4,1),xuz(1,7),borncurr(0,0,1,31))
         dipole(14) = - dipoleVVj_CS66(7,0,p(0,2,1),p(0,3,1),p(0,4,1),xuz(1,7),borncurr(0,0,1,31))

         res(1,7) = diag(1) * dipole(1) * 
     &      ( pdf(2*sign1,1,7)*pdf(2*sign2,2,7) + 
     &        pdf(4*sign1,1,7)*pdf(4*sign2,2,7) )/2d0*polcolgq
         res(2,7) = diag(2) * dipole(1) * 
     &      ( pdf(2*sign1,1,7)*pdf(1*sign2,2,7) + 
     &        pdf(4*sign1,1,7)*pdf(3*sign2,2,7) )*polcolgq
         res(5,7) = diag(5) * dipole(2) * 
     &      ( pdf(2*sign1,1,7)*pdf(0*sign2,2,7) + 
     &        pdf(4*sign1,1,7)*pdf(0*sign2,2,7) )*polcolqq
         res(13,7) = diag(13) * dipole(3) * 
     &      ( pdf(2*sign1,1,7)*pdf(-2*sign2,2,7) + 
     &        pdf(4*sign1,1,7)*pdf(-4*sign2,2,7) )/2d0*polcolqq
         res(15,7) = diag(15) * dipole(4) * 
     &      ( pdf(1*sign1,1,7)*pdf(1*sign2,2,7) + 
     &        pdf(3*sign1,1,7)*pdf(3*sign2,2,7) + 
     &        lbquark*pdf(5*sign1,1,7)*pdf(5*sign2,2,7) )/2d0*polcolgq
         res(18,7) = diag(18) * dipole(5) * 
     &      ( pdf(1*sign1,1,7)*pdf(0*sign2,2,7) + 
     &        pdf(3*sign1,1,7)*pdf(0*sign2,2,7) + 
     &        lbquark*pdf(5*sign1,1,7)*pdf(0*sign2,2,7) )*polcolqq
         res(25,7) = diag(25) * dipole(6) * 
     &      ( pdf(1*sign1,1,7)*pdf(-1*sign2,2,7) + 
     &        pdf(3*sign1,1,7)*pdf(-3*sign2,2,7) + 
     &        lbquark*pdf(5*sign1,1,7)*pdf(-5*sign2,2,7) )/2d0*polcolqq
         res(29,7) = diag(29) * dipole(7) * 
     &      ( pdf(0*sign1,1,7)*pdf(0*sign2,2,7) + 
     &        pdf(0*sign1,1,7)*pdf(0*sign2,2,7) )*polcolgq
         res(30,7) = diag(30) * dipole(8) * 
     &      ( pdf(0*sign1,1,7)*pdf(0*sign2,2,7) + 
     &        pdf(0*sign1,1,7)*pdf(0*sign2,2,7) + 
     &        lbquark*pdf(0*sign1,1,7)*pdf(0*sign2,2,7) )*polcolgq
         res(33,7) = diag(33) * dipole(9) * 
     &      ( pdf(-1*sign1,1,7)*pdf(4*sign2,2,7) + 
     &        pdf(-3*sign1,1,7)*pdf(2*sign2,2,7) + 
     &        lbquark*pdf(-5*sign1,1,7)*pdf(2*sign2,2,7) + 
     &        lbquark*pdf(-5*sign1,1,7)*pdf(4*sign2,2,7) )*polcolgq
         res(34,7) = diag(34) * dipole(9) * 
     &      ( pdf(-1*sign1,1,7)*pdf(3*sign2,2,7) + 
     &        lbquark*pdf(-1*sign1,1,7)*pdf(5*sign2,2,7) + 
     &        pdf(-3*sign1,1,7)*pdf(1*sign2,2,7) + 
     &        lbquark*pdf(-3*sign1,1,7)*pdf(5*sign2,2,7) + 
     &        lbquark*pdf(-5*sign1,1,7)*pdf(1*sign2,2,7) + 
     &        lbquark*pdf(-5*sign1,1,7)*pdf(3*sign2,2,7) )*polcolgq
         res(35,7) = diag(35) * dipole(9) * 
     &      ( pdf(-1*sign1,1,7)*pdf(2*sign2,2,7) + 
     &        pdf(-3*sign1,1,7)*pdf(4*sign2,2,7) )*polcolgq
         res(39,7) = diag(39) * dipole(9) * 
     &      ( pdf(-1*sign1,1,7)*pdf(1*sign2,2,7) + 
     &        pdf(-3*sign1,1,7)*pdf(3*sign2,2,7) + 
     &        lbquark*pdf(-5*sign1,1,7)*pdf(5*sign2,2,7) )*polcolgq
         res(40,7) = diag(40) * dipole(10) * 
     &      ( pdf(-1*sign1,1,7)*pdf(1*sign2,2,7) + 
     &        pdf(-3*sign1,1,7)*pdf(3*sign2,2,7) + 
     &        lbquark*pdf(-5*sign1,1,7)*pdf(5*sign2,2,7) )/2d0*polcolqq
         res(41,7) = diag(41) * dipole(11) * 
     &      ( pdf(-1*sign1,1,7)*pdf(0*sign2,2,7) + 
     &        pdf(-3*sign1,1,7)*pdf(0*sign2,2,7) + 
     &        lbquark*pdf(-5*sign1,1,7)*pdf(0*sign2,2,7) )*polcolqq
         res(42,7) = diag(42) * dipole(9) * 
     &      ( pdf(-1*sign1,1,7)*pdf(-1*sign2,2,7) + 
     &        pdf(-3*sign1,1,7)*pdf(-3*sign2,2,7) + 
     &        lbquark*pdf(-5*sign1,1,7)*pdf(-5*sign2,2,7) )/2d0*polcolgq
         res(43,7) = diag(43) * dipole(9) * 
     &      ( pdf(-1*sign1,1,7)*pdf(-2*sign2,2,7) + 
     &        pdf(-3*sign1,1,7)*pdf(-4*sign2,2,7) )*polcolgq
         res(44,7) = diag(44) * dipole(9) * 
     &      ( pdf(-1*sign1,1,7)*pdf(-3*sign2,2,7) + 
     &        lbquark*pdf(-1*sign1,1,7)*pdf(-5*sign2,2,7) + 
     &        pdf(-3*sign1,1,7)*pdf(-1*sign2,2,7) + 
     &        lbquark*pdf(-3*sign1,1,7)*pdf(-5*sign2,2,7) + 
     &        lbquark*pdf(-5*sign1,1,7)*pdf(-1*sign2,2,7) + 
     &        lbquark*pdf(-5*sign1,1,7)*pdf(-3*sign2,2,7) )*polcolgq
         res(45,7) = diag(45) * dipole(9) * 
     &      ( pdf(-1*sign1,1,7)*pdf(-4*sign2,2,7) + 
     &        pdf(-3*sign1,1,7)*pdf(-2*sign2,2,7) + 
     &        lbquark*pdf(-5*sign1,1,7)*pdf(-2*sign2,2,7) + 
     &        lbquark*pdf(-5*sign1,1,7)*pdf(-4*sign2,2,7) )*polcolgq
         res(46,7) = diag(46) * dipole(12) * 
     &      ( pdf(-2*sign1,1,7)*pdf(4*sign2,2,7) + 
     &        pdf(-4*sign1,1,7)*pdf(2*sign2,2,7) )*polcolgq
         res(47,7) = diag(47) * dipole(12) * 
     &      ( pdf(-2*sign1,1,7)*pdf(3*sign2,2,7) + 
     &        lbquark*pdf(-2*sign1,1,7)*pdf(5*sign2,2,7) + 
     &        pdf(-4*sign1,1,7)*pdf(1*sign2,2,7) + 
     &        lbquark*pdf(-4*sign1,1,7)*pdf(5*sign2,2,7) )*polcolgq
         res(50,7) = diag(50) * dipole(12) * 
     &      ( pdf(-2*sign1,1,7)*pdf(2*sign2,2,7) + 
     &        pdf(-4*sign1,1,7)*pdf(4*sign2,2,7) )*polcolgq
         res(52,7) = diag(52) * dipole(13) * 
     &      ( pdf(-2*sign1,1,7)*pdf(2*sign2,2,7) + 
     &        pdf(-4*sign1,1,7)*pdf(4*sign2,2,7) )/2d0*polcolqq
         res(53,7) = diag(53) * dipole(12) * 
     &      ( pdf(-2*sign1,1,7)*pdf(1*sign2,2,7) + 
     &        pdf(-4*sign1,1,7)*pdf(3*sign2,2,7) )*polcolgq
         res(54,7) = diag(54) * dipole(14) * 
     &      ( pdf(-2*sign1,1,7)*pdf(0*sign2,2,7) + 
     &        pdf(-4*sign1,1,7)*pdf(0*sign2,2,7) )*polcolqq
         res(56,7) = diag(56) * dipole(12) * 
     &      ( pdf(-2*sign1,1,7)*pdf(-2*sign2,2,7) + 
     &        pdf(-4*sign1,1,7)*pdf(-4*sign2,2,7) )/2d0*polcolgq
         res(57,7) = diag(57) * dipole(12) * 
     &      ( pdf(-2*sign1,1,7)*pdf(-3*sign2,2,7) + 
     &        lbquark*pdf(-2*sign1,1,7)*pdf(-5*sign2,2,7) + 
     &        pdf(-4*sign1,1,7)*pdf(-1*sign2,2,7) + 
     &        lbquark*pdf(-4*sign1,1,7)*pdf(-5*sign2,2,7) )*polcolgq
         res(58,7) = diag(58) * dipole(12) * 
     &      ( pdf(-2*sign1,1,7)*pdf(-4*sign2,2,7) + 
     &        pdf(-4*sign1,1,7)*pdf(-2*sign2,2,7) )*polcolgq
       endif 

       if (lok(8)) then 
         call qqgZAZA_curr(borncurr(0,-1,1,35),8,p(0,3,8),p(0,1,8),p(0,2,8),-1,-1,-1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,36),8,p(0,2,8),p(0,1,8),p(0,3,8),1,-1,1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,37),8,p(0,3,8),p(0,2,8),p(0,1,8),-1,-1,-1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,38),8,p(0,2,8),p(0,3,8),p(0,1,8),1,1,-1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,39),8,p(0,1,8),p(0,2,8),p(0,3,8),1,-1,1,bos)

         dipole(1) = - dipoleVVj_CS147(8,1,p(0,1,1),p(0,3,1),p(0,2,1),xuz(1,8),borncurr(0,-1,1,38))
         dipole(2) = - dipoleVVj_CS147(8,1,p(0,1,1),p(0,3,1),p(0,2,1),xuz(1,8),borncurr(0,-1,2,38))
         dipole(3) = - dipoleVVj_CS147(8,1,p(0,1,1),p(0,3,1),p(0,2,1),xuz(1,8),borncurr(0,-1,2,37))
         dipole(4) = - dipoleVVj_CS147(8,1,p(0,1,1),p(0,3,1),p(0,2,1),xuz(1,8),borncurr(0,-1,1,37))
         dipole(5) = - dipoleVVj_CS145(8,1,p(0,1,1),p(0,3,1),p(0,2,1),xuz(1,8),borncurr(0,0,1,39))
         dipole(6) = - dipoleVVj_CS145(8,1,p(0,1,1),p(0,3,1),p(0,2,1),xuz(1,8),borncurr(0,0,2,39))
         dipole(7) = - dipoleVVj_CS146(8,1,p(0,1,1),p(0,3,1),p(0,2,1),xuz(1,8),borncurr(0,0,1,36))
         dipole(8) = - dipoleVVj_CS146(8,1,p(0,1,1),p(0,3,1),p(0,2,1),xuz(1,8),borncurr(0,0,2,36))
         dipole(9) = - dipoleVVj_CS146(8,0,p(0,1,1),p(0,3,1),p(0,2,1),xuz(1,8),borncurr(0,0,1,35))
         dipole(10) = - dipoleVVj_CS146(8,0,p(0,1,1),p(0,3,1),p(0,2,1),xuz(1,8),borncurr(0,0,2,35))
         dipole(11) = - dipoleVVj_CS146(8,1,p(0,1,1),p(0,3,1),p(0,2,1),xuz(1,8),borncurr(0,0,2,39))
         dipole(12) = - dipoleVVj_CS146(8,1,p(0,1,1),p(0,3,1),p(0,2,1),xuz(1,8),borncurr(0,0,1,39))
         dipole(13) = - dipoleVVj_CS145(8,1,p(0,1,1),p(0,3,1),p(0,2,1),xuz(1,8),borncurr(0,0,2,36))
         dipole(14) = - dipoleVVj_CS145(8,1,p(0,1,1),p(0,3,1),p(0,2,1),xuz(1,8),borncurr(0,0,1,36))

         res(1,8) = diag(1) * dipole(1) * 
     &      ( pdf(2*sign1,1,8)*pdf(2*sign2,2,8) + 
     &        pdf(4*sign1,1,8)*pdf(4*sign2,2,8) )/2d0*polcolgq
         res(3,8) = diag(3) * dipole(1) * 
     &      ( pdf(2*sign1,1,8)*pdf(4*sign2,2,8) + 
     &        pdf(4*sign1,1,8)*pdf(2*sign2,2,8) )*polcolgq
         res(4,8) = diag(4) * dipole(2) * 
     &      ( pdf(2*sign1,1,8)*pdf(3*sign2,2,8) + 
     &        lbquark*pdf(2*sign1,1,8)*pdf(5*sign2,2,8) + 
     &        pdf(4*sign1,1,8)*pdf(1*sign2,2,8) + 
     &        lbquark*pdf(4*sign1,1,8)*pdf(5*sign2,2,8) )*polcolgq
         res(6,8) = diag(6) * dipole(3) * 
     &      ( pdf(2*sign1,1,8)*pdf(-3*sign2,2,8) + 
     &        lbquark*pdf(2*sign1,1,8)*pdf(-5*sign2,2,8) + 
     &        pdf(4*sign1,1,8)*pdf(-1*sign2,2,8) + 
     &        lbquark*pdf(4*sign1,1,8)*pdf(-5*sign2,2,8) )*polcolgq
         res(7,8) = diag(7) * dipole(4) * 
     &      ( pdf(2*sign1,1,8)*pdf(-4*sign2,2,8) + 
     &        pdf(4*sign1,1,8)*pdf(-2*sign2,2,8) )*polcolgq
         res(8,8) = diag(8) * dipole(3) * 
     &      ( pdf(2*sign1,1,8)*pdf(-1*sign2,2,8) + 
     &        pdf(4*sign1,1,8)*pdf(-3*sign2,2,8) )*polcolgq
         res(9,8) = diag(9) * dipole(4) * 
     &      ( pdf(2*sign1,1,8)*pdf(-2*sign2,2,8) + 
     &        pdf(4*sign1,1,8)*pdf(-4*sign2,2,8) )*polcolgq
         res(13,8) = diag(13) * dipole(5) * 
     &      ( pdf(2*sign1,1,8)*pdf(-2*sign2,2,8) + 
     &        pdf(4*sign1,1,8)*pdf(-4*sign2,2,8) )/2d0*polcolqq
         res(14,8) = diag(14) * dipole(1) * 
     &      ( pdf(1*sign1,1,8)*pdf(2*sign2,2,8) + 
     &        pdf(3*sign1,1,8)*pdf(4*sign2,2,8) )*polcolgq
         res(15,8) = diag(15) * dipole(2) * 
     &      ( pdf(1*sign1,1,8)*pdf(1*sign2,2,8) + 
     &        pdf(3*sign1,1,8)*pdf(3*sign2,2,8) + 
     &        lbquark*pdf(5*sign1,1,8)*pdf(5*sign2,2,8) )/2d0*polcolgq
         res(16,8) = diag(16) * dipole(1) * 
     &      ( pdf(1*sign1,1,8)*pdf(4*sign2,2,8) + 
     &        pdf(3*sign1,1,8)*pdf(2*sign2,2,8) + 
     &        lbquark*pdf(5*sign1,1,8)*pdf(2*sign2,2,8) + 
     &        lbquark*pdf(5*sign1,1,8)*pdf(4*sign2,2,8) )*polcolgq
         res(17,8) = diag(17) * dipole(2) * 
     &      ( pdf(1*sign1,1,8)*pdf(3*sign2,2,8) + 
     &        lbquark*pdf(1*sign1,1,8)*pdf(5*sign2,2,8) + 
     &        pdf(3*sign1,1,8)*pdf(1*sign2,2,8) + 
     &        lbquark*pdf(3*sign1,1,8)*pdf(5*sign2,2,8) + 
     &        lbquark*pdf(5*sign1,1,8)*pdf(1*sign2,2,8) + 
     &        lbquark*pdf(5*sign1,1,8)*pdf(3*sign2,2,8) )*polcolgq
         res(19,8) = diag(19) * dipole(3) * 
     &      ( pdf(1*sign1,1,8)*pdf(-3*sign2,2,8) + 
     &        lbquark*pdf(1*sign1,1,8)*pdf(-5*sign2,2,8) + 
     &        pdf(3*sign1,1,8)*pdf(-1*sign2,2,8) + 
     &        lbquark*pdf(3*sign1,1,8)*pdf(-5*sign2,2,8) + 
     &        lbquark*pdf(5*sign1,1,8)*pdf(-1*sign2,2,8) + 
     &        lbquark*pdf(5*sign1,1,8)*pdf(-3*sign2,2,8) )*polcolgq
         res(20,8) = diag(20) * dipole(4) * 
     &      ( pdf(1*sign1,1,8)*pdf(-4*sign2,2,8) + 
     &        pdf(3*sign1,1,8)*pdf(-2*sign2,2,8) + 
     &        lbquark*pdf(5*sign1,1,8)*pdf(-2*sign2,2,8) + 
     &        lbquark*pdf(5*sign1,1,8)*pdf(-4*sign2,2,8) )*polcolgq
         res(22,8) = diag(22) * dipole(3) * 
     &      ( pdf(1*sign1,1,8)*pdf(-1*sign2,2,8) + 
     &        pdf(3*sign1,1,8)*pdf(-3*sign2,2,8) + 
     &        lbquark*pdf(5*sign1,1,8)*pdf(-5*sign2,2,8) )*polcolgq
         res(25,8) = diag(25) * dipole(6) * 
     &      ( pdf(1*sign1,1,8)*pdf(-1*sign2,2,8) + 
     &        pdf(3*sign1,1,8)*pdf(-3*sign2,2,8) + 
     &        lbquark*pdf(5*sign1,1,8)*pdf(-5*sign2,2,8) )/2d0*polcolqq
         res(26,8) = diag(26) * dipole(4) * 
     &      ( pdf(1*sign1,1,8)*pdf(-2*sign2,2,8) + 
     &        pdf(3*sign1,1,8)*pdf(-4*sign2,2,8) )*polcolgq
         res(27,8) = diag(27) * dipole(7) * 
     &      ( pdf(0*sign1,1,8)*pdf(2*sign2,2,8) + 
     &        pdf(0*sign1,1,8)*pdf(4*sign2,2,8) )*polcolqq
         res(28,8) = diag(28) * dipole(8) * 
     &      ( pdf(0*sign1,1,8)*pdf(1*sign2,2,8) + 
     &        pdf(0*sign1,1,8)*pdf(3*sign2,2,8) + 
     &        lbquark*pdf(0*sign1,1,8)*pdf(5*sign2,2,8) )*polcolqq
         res(29,8) = diag(29) * dipole(9) * 
     &      ( pdf(0*sign1,1,8)*pdf(0*sign2,2,8) + 
     &        pdf(0*sign1,1,8)*pdf(0*sign2,2,8) )*polcolgq
         res(30,8) = diag(30) * dipole(10) * 
     &      ( pdf(0*sign1,1,8)*pdf(0*sign2,2,8) + 
     &        pdf(0*sign1,1,8)*pdf(0*sign2,2,8) + 
     &        lbquark*pdf(0*sign1,1,8)*pdf(0*sign2,2,8) )*polcolgq
         res(31,8) = diag(31) * dipole(11) * 
     &      ( pdf(0*sign1,1,8)*pdf(-1*sign2,2,8) + 
     &        pdf(0*sign1,1,8)*pdf(-3*sign2,2,8) + 
     &        lbquark*pdf(0*sign1,1,8)*pdf(-5*sign2,2,8) )*polcolqq
         res(32,8) = diag(32) * dipole(12) * 
     &      ( pdf(0*sign1,1,8)*pdf(-2*sign2,2,8) + 
     &        pdf(0*sign1,1,8)*pdf(-4*sign2,2,8) )*polcolqq
         res(40,8) = diag(40) * dipole(13) * 
     &      ( pdf(-1*sign1,1,8)*pdf(1*sign2,2,8) + 
     &        pdf(-3*sign1,1,8)*pdf(3*sign2,2,8) + 
     &        lbquark*pdf(-5*sign1,1,8)*pdf(5*sign2,2,8) )/2d0*polcolqq
         res(42,8) = diag(42) * dipole(3) * 
     &      ( pdf(-1*sign1,1,8)*pdf(-1*sign2,2,8) + 
     &        pdf(-3*sign1,1,8)*pdf(-3*sign2,2,8) + 
     &        lbquark*pdf(-5*sign1,1,8)*pdf(-5*sign2,2,8) )/2d0*polcolgq
         res(52,8) = diag(52) * dipole(14) * 
     &      ( pdf(-2*sign1,1,8)*pdf(2*sign2,2,8) + 
     &        pdf(-4*sign1,1,8)*pdf(4*sign2,2,8) )/2d0*polcolqq
         res(55,8) = diag(55) * dipole(3) * 
     &      ( pdf(-2*sign1,1,8)*pdf(-1*sign2,2,8) + 
     &        pdf(-4*sign1,1,8)*pdf(-3*sign2,2,8) )*polcolgq
         res(56,8) = diag(56) * dipole(4) * 
     &      ( pdf(-2*sign1,1,8)*pdf(-2*sign2,2,8) + 
     &        pdf(-4*sign1,1,8)*pdf(-4*sign2,2,8) )/2d0*polcolgq
       endif 

       if (lok(9)) then 
         call qqgZAZA_curr(borncurr(0,-1,1,40),9,p(0,3,9),p(0,1,9),p(0,2,9),-1,-1,-1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,41),9,p(0,2,9),p(0,1,9),p(0,3,9),1,-1,1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,42),9,p(0,3,9),p(0,2,9),p(0,1,9),-1,-1,-1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,43),9,p(0,1,9),p(0,2,9),p(0,3,9),1,-1,1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,44),9,p(0,1,9),p(0,3,9),p(0,2,9),1,1,-1,bos)

         dipole(1) = - dipoleVVj_CS147(9,1,p(0,2,1),p(0,3,1),p(0,1,1),xuz(1,9),borncurr(0,-1,1,44))
         dipole(2) = - dipoleVVj_CS146(9,1,p(0,2,1),p(0,3,1),p(0,1,1),xuz(1,9),borncurr(0,0,1,43))
         dipole(3) = - dipoleVVj_CS145(9,1,p(0,2,1),p(0,3,1),p(0,1,1),xuz(1,9),borncurr(0,0,1,43))
         dipole(4) = - dipoleVVj_CS147(9,1,p(0,2,1),p(0,3,1),p(0,1,1),xuz(1,9),borncurr(0,-1,2,44))
         dipole(5) = - dipoleVVj_CS146(9,1,p(0,2,1),p(0,3,1),p(0,1,1),xuz(1,9),borncurr(0,0,2,43))
         dipole(6) = - dipoleVVj_CS145(9,1,p(0,2,1),p(0,3,1),p(0,1,1),xuz(1,9),borncurr(0,0,2,43))
         dipole(7) = - dipoleVVj_CS146(9,0,p(0,2,1),p(0,3,1),p(0,1,1),xuz(1,9),borncurr(0,0,1,42))
         dipole(8) = - dipoleVVj_CS146(9,0,p(0,2,1),p(0,3,1),p(0,1,1),xuz(1,9),borncurr(0,0,2,42))
         dipole(9) = - dipoleVVj_CS147(9,1,p(0,2,1),p(0,3,1),p(0,1,1),xuz(1,9),borncurr(0,-1,2,40))
         dipole(10) = - dipoleVVj_CS145(9,1,p(0,2,1),p(0,3,1),p(0,1,1),xuz(1,9),borncurr(0,0,2,41))
         dipole(11) = - dipoleVVj_CS146(9,1,p(0,2,1),p(0,3,1),p(0,1,1),xuz(1,9),borncurr(0,0,2,41))
         dipole(12) = - dipoleVVj_CS147(9,1,p(0,2,1),p(0,3,1),p(0,1,1),xuz(1,9),borncurr(0,-1,1,40))
         dipole(13) = - dipoleVVj_CS145(9,1,p(0,2,1),p(0,3,1),p(0,1,1),xuz(1,9),borncurr(0,0,1,41))
         dipole(14) = - dipoleVVj_CS146(9,1,p(0,2,1),p(0,3,1),p(0,1,1),xuz(1,9),borncurr(0,0,1,41))

         res(1,9) = diag(1) * dipole(1) * 
     &      ( pdf(2*sign1,1,9)*pdf(2*sign2,2,9) + 
     &        pdf(4*sign1,1,9)*pdf(4*sign2,2,9) )/2d0*polcolgq
         res(2,9) = diag(2) * dipole(1) * 
     &      ( pdf(2*sign1,1,9)*pdf(1*sign2,2,9) + 
     &        pdf(4*sign1,1,9)*pdf(3*sign2,2,9) )*polcolgq
         res(5,9) = diag(5) * dipole(2) * 
     &      ( pdf(2*sign1,1,9)*pdf(0*sign2,2,9) + 
     &        pdf(4*sign1,1,9)*pdf(0*sign2,2,9) )*polcolqq
         res(13,9) = diag(13) * dipole(3) * 
     &      ( pdf(2*sign1,1,9)*pdf(-2*sign2,2,9) + 
     &        pdf(4*sign1,1,9)*pdf(-4*sign2,2,9) )/2d0*polcolqq
         res(15,9) = diag(15) * dipole(4) * 
     &      ( pdf(1*sign1,1,9)*pdf(1*sign2,2,9) + 
     &        pdf(3*sign1,1,9)*pdf(3*sign2,2,9) + 
     &        lbquark*pdf(5*sign1,1,9)*pdf(5*sign2,2,9) )/2d0*polcolgq
         res(18,9) = diag(18) * dipole(5) * 
     &      ( pdf(1*sign1,1,9)*pdf(0*sign2,2,9) + 
     &        pdf(3*sign1,1,9)*pdf(0*sign2,2,9) + 
     &        lbquark*pdf(5*sign1,1,9)*pdf(0*sign2,2,9) )*polcolqq
         res(25,9) = diag(25) * dipole(6) * 
     &      ( pdf(1*sign1,1,9)*pdf(-1*sign2,2,9) + 
     &        pdf(3*sign1,1,9)*pdf(-3*sign2,2,9) + 
     &        lbquark*pdf(5*sign1,1,9)*pdf(-5*sign2,2,9) )/2d0*polcolqq
         res(29,9) = diag(29) * dipole(7) * 
     &      ( pdf(0*sign1,1,9)*pdf(0*sign2,2,9) + 
     &        pdf(0*sign1,1,9)*pdf(0*sign2,2,9) )*polcolgq
         res(30,9) = diag(30) * dipole(8) * 
     &      ( pdf(0*sign1,1,9)*pdf(0*sign2,2,9) + 
     &        pdf(0*sign1,1,9)*pdf(0*sign2,2,9) + 
     &        lbquark*pdf(0*sign1,1,9)*pdf(0*sign2,2,9) )*polcolgq
         res(33,9) = diag(33) * dipole(9) * 
     &      ( pdf(-1*sign1,1,9)*pdf(4*sign2,2,9) + 
     &        pdf(-3*sign1,1,9)*pdf(2*sign2,2,9) + 
     &        lbquark*pdf(-5*sign1,1,9)*pdf(2*sign2,2,9) + 
     &        lbquark*pdf(-5*sign1,1,9)*pdf(4*sign2,2,9) )*polcolgq
         res(34,9) = diag(34) * dipole(9) * 
     &      ( pdf(-1*sign1,1,9)*pdf(3*sign2,2,9) + 
     &        lbquark*pdf(-1*sign1,1,9)*pdf(5*sign2,2,9) + 
     &        pdf(-3*sign1,1,9)*pdf(1*sign2,2,9) + 
     &        lbquark*pdf(-3*sign1,1,9)*pdf(5*sign2,2,9) + 
     &        lbquark*pdf(-5*sign1,1,9)*pdf(1*sign2,2,9) + 
     &        lbquark*pdf(-5*sign1,1,9)*pdf(3*sign2,2,9) )*polcolgq
         res(35,9) = diag(35) * dipole(9) * 
     &      ( pdf(-1*sign1,1,9)*pdf(2*sign2,2,9) + 
     &        pdf(-3*sign1,1,9)*pdf(4*sign2,2,9) )*polcolgq
         res(39,9) = diag(39) * dipole(9) * 
     &      ( pdf(-1*sign1,1,9)*pdf(1*sign2,2,9) + 
     &        pdf(-3*sign1,1,9)*pdf(3*sign2,2,9) + 
     &        lbquark*pdf(-5*sign1,1,9)*pdf(5*sign2,2,9) )*polcolgq
         res(40,9) = diag(40) * dipole(10) * 
     &      ( pdf(-1*sign1,1,9)*pdf(1*sign2,2,9) + 
     &        pdf(-3*sign1,1,9)*pdf(3*sign2,2,9) + 
     &        lbquark*pdf(-5*sign1,1,9)*pdf(5*sign2,2,9) )/2d0*polcolqq
         res(41,9) = diag(41) * dipole(11) * 
     &      ( pdf(-1*sign1,1,9)*pdf(0*sign2,2,9) + 
     &        pdf(-3*sign1,1,9)*pdf(0*sign2,2,9) + 
     &        lbquark*pdf(-5*sign1,1,9)*pdf(0*sign2,2,9) )*polcolqq
         res(42,9) = diag(42) * dipole(9) * 
     &      ( pdf(-1*sign1,1,9)*pdf(-1*sign2,2,9) + 
     &        pdf(-3*sign1,1,9)*pdf(-3*sign2,2,9) + 
     &        lbquark*pdf(-5*sign1,1,9)*pdf(-5*sign2,2,9) )/2d0*polcolgq
         res(43,9) = diag(43) * dipole(9) * 
     &      ( pdf(-1*sign1,1,9)*pdf(-2*sign2,2,9) + 
     &        pdf(-3*sign1,1,9)*pdf(-4*sign2,2,9) )*polcolgq
         res(44,9) = diag(44) * dipole(9) * 
     &      ( pdf(-1*sign1,1,9)*pdf(-3*sign2,2,9) + 
     &        lbquark*pdf(-1*sign1,1,9)*pdf(-5*sign2,2,9) + 
     &        pdf(-3*sign1,1,9)*pdf(-1*sign2,2,9) + 
     &        lbquark*pdf(-3*sign1,1,9)*pdf(-5*sign2,2,9) + 
     &        lbquark*pdf(-5*sign1,1,9)*pdf(-1*sign2,2,9) + 
     &        lbquark*pdf(-5*sign1,1,9)*pdf(-3*sign2,2,9) )*polcolgq
         res(45,9) = diag(45) * dipole(9) * 
     &      ( pdf(-1*sign1,1,9)*pdf(-4*sign2,2,9) + 
     &        pdf(-3*sign1,1,9)*pdf(-2*sign2,2,9) + 
     &        lbquark*pdf(-5*sign1,1,9)*pdf(-2*sign2,2,9) + 
     &        lbquark*pdf(-5*sign1,1,9)*pdf(-4*sign2,2,9) )*polcolgq
         res(46,9) = diag(46) * dipole(12) * 
     &      ( pdf(-2*sign1,1,9)*pdf(4*sign2,2,9) + 
     &        pdf(-4*sign1,1,9)*pdf(2*sign2,2,9) )*polcolgq
         res(47,9) = diag(47) * dipole(12) * 
     &      ( pdf(-2*sign1,1,9)*pdf(3*sign2,2,9) + 
     &        lbquark*pdf(-2*sign1,1,9)*pdf(5*sign2,2,9) + 
     &        pdf(-4*sign1,1,9)*pdf(1*sign2,2,9) + 
     &        lbquark*pdf(-4*sign1,1,9)*pdf(5*sign2,2,9) )*polcolgq
         res(50,9) = diag(50) * dipole(12) * 
     &      ( pdf(-2*sign1,1,9)*pdf(2*sign2,2,9) + 
     &        pdf(-4*sign1,1,9)*pdf(4*sign2,2,9) )*polcolgq
         res(52,9) = diag(52) * dipole(13) * 
     &      ( pdf(-2*sign1,1,9)*pdf(2*sign2,2,9) + 
     &        pdf(-4*sign1,1,9)*pdf(4*sign2,2,9) )/2d0*polcolqq
         res(53,9) = diag(53) * dipole(12) * 
     &      ( pdf(-2*sign1,1,9)*pdf(1*sign2,2,9) + 
     &        pdf(-4*sign1,1,9)*pdf(3*sign2,2,9) )*polcolgq
         res(54,9) = diag(54) * dipole(14) * 
     &      ( pdf(-2*sign1,1,9)*pdf(0*sign2,2,9) + 
     &        pdf(-4*sign1,1,9)*pdf(0*sign2,2,9) )*polcolqq
         res(56,9) = diag(56) * dipole(12) * 
     &      ( pdf(-2*sign1,1,9)*pdf(-2*sign2,2,9) + 
     &        pdf(-4*sign1,1,9)*pdf(-4*sign2,2,9) )/2d0*polcolgq
         res(57,9) = diag(57) * dipole(12) * 
     &      ( pdf(-2*sign1,1,9)*pdf(-3*sign2,2,9) + 
     &        lbquark*pdf(-2*sign1,1,9)*pdf(-5*sign2,2,9) + 
     &        pdf(-4*sign1,1,9)*pdf(-1*sign2,2,9) + 
     &        lbquark*pdf(-4*sign1,1,9)*pdf(-5*sign2,2,9) )*polcolgq
         res(58,9) = diag(58) * dipole(12) * 
     &      ( pdf(-2*sign1,1,9)*pdf(-4*sign2,2,9) + 
     &        pdf(-4*sign1,1,9)*pdf(-2*sign2,2,9) )*polcolgq
       endif 

       if (lok(10)) then 
         call qqgZAZA_curr(borncurr(0,-1,1,45),10,p(0,3,10),p(0,1,10),p(0,2,10),-1,-1,-1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,46),10,p(0,2,10),p(0,1,10),p(0,3,10),1,-1,1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,47),10,p(0,3,10),p(0,2,10),p(0,1,10),-1,-1,-1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,48),10,p(0,2,10),p(0,3,10),p(0,1,10),1,1,-1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,49),10,p(0,1,10),p(0,2,10),p(0,3,10),1,-1,1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,50),10,p(0,1,10),p(0,3,10),p(0,2,10),1,1,-1,bos)

         dipole(1) = - dipoleVVj_CS39(10,1,p(0,3,1),p(0,4,1),p(0,1,1),xuz(1,10),borncurr(0,0,1,50))
         dipole(2) = - dipoleVVj_CS41(10,1,p(0,3,1),p(0,4,1),p(0,1,1),xuz(1,10),borncurr(0,-1,1,49))
         dipole(3) = - dipoleVVj_CS40(10,1,p(0,3,1),p(0,4,1),p(0,1,1),xuz(1,10),borncurr(0,-1,1,49))
         dipole(4) = - dipoleVVj_CS39(10,1,p(0,3,1),p(0,4,1),p(0,1,1),xuz(1,10),borncurr(0,0,2,50))
         dipole(5) = - dipoleVVj_CS41(10,1,p(0,3,1),p(0,4,1),p(0,1,1),xuz(1,10),borncurr(0,-1,2,49))
         dipole(6) = - dipoleVVj_CS40(10,1,p(0,3,1),p(0,4,1),p(0,1,1),xuz(1,10),borncurr(0,-1,2,49))
         dipole(7) = - dipoleVVj_CS39(10,0,p(0,3,1),p(0,4,1),p(0,1,1),xuz(1,10),borncurr(0,0,1,48))
         dipole(8) = - dipoleVVj_CS39(10,0,p(0,3,1),p(0,4,1),p(0,1,1),xuz(1,10),borncurr(0,0,2,48))
         dipole(9) = - dipoleVVj_CS39(10,0,p(0,3,1),p(0,4,1),p(0,1,1),xuz(1,10),borncurr(0,0,2,47))
         dipole(10) = - dipoleVVj_CS39(10,0,p(0,3,1),p(0,4,1),p(0,1,1),xuz(1,10),borncurr(0,0,1,47))
         dipole(11) = - dipoleVVj_CS41(10,1,p(0,3,1),p(0,4,1),p(0,1,1),xuz(1,10),borncurr(0,-1,2,46))
         dipole(12) = - dipoleVVj_CS40(10,1,p(0,3,1),p(0,4,1),p(0,1,1),xuz(1,10),borncurr(0,-1,2,46))
         dipole(13) = - dipoleVVj_CS39(10,1,p(0,3,1),p(0,4,1),p(0,1,1),xuz(1,10),borncurr(0,0,2,45))
         dipole(14) = - dipoleVVj_CS41(10,1,p(0,3,1),p(0,4,1),p(0,1,1),xuz(1,10),borncurr(0,-1,1,46))
         dipole(15) = - dipoleVVj_CS40(10,1,p(0,3,1),p(0,4,1),p(0,1,1),xuz(1,10),borncurr(0,-1,1,46))
         dipole(16) = - dipoleVVj_CS39(10,1,p(0,3,1),p(0,4,1),p(0,1,1),xuz(1,10),borncurr(0,0,1,45))

         res(5,10) = diag(5) * dipole(1) * 
     &      ( pdf(2*sign1,1,10)*pdf(0*sign2,2,10) + 
     &        pdf(4*sign1,1,10)*pdf(0*sign2,2,10) )*polcolgq
         res(9,10) = diag(9) * dipole(2) * 
     &      ( pdf(2*sign1,1,10)*pdf(-2*sign2,2,10) + 
     &        pdf(4*sign1,1,10)*pdf(-4*sign2,2,10) )*polcolqq
         res(10,10) = diag(10) * dipole(2) * 
     &      ( pdf(2*sign1,1,10)*pdf(-2*sign2,2,10) + 
     &        pdf(4*sign1,1,10)*pdf(-4*sign2,2,10) )*polcolqq
         res(11,10) = diag(11) * dipole(2) * 
     &      ( pdf(2*sign1,1,10)*pdf(-2*sign2,2,10) + 
     &        pdf(4*sign1,1,10)*pdf(-4*sign2,2,10) )*polcolqq
         res(12,10) = diag(12) * dipole(2) * 
     &      ( pdf(2*sign1,1,10)*pdf(-2*sign2,2,10) + 
     &        lbquark*pdf(2*sign1,1,10)*pdf(-2*sign2,2,10) + 
     &        pdf(4*sign1,1,10)*pdf(-4*sign2,2,10) + 
     &        lbquark*pdf(4*sign1,1,10)*pdf(-4*sign2,2,10) )*polcolqq
         res(13,10) = diag(13) * dipole(3) * 
     &      ( pdf(2*sign1,1,10)*pdf(-2*sign2,2,10) + 
     &        pdf(4*sign1,1,10)*pdf(-4*sign2,2,10) )/2d0*polcolqq
         res(18,10) = diag(18) * dipole(4) * 
     &      ( pdf(1*sign1,1,10)*pdf(0*sign2,2,10) + 
     &        pdf(3*sign1,1,10)*pdf(0*sign2,2,10) + 
     &        lbquark*pdf(5*sign1,1,10)*pdf(0*sign2,2,10) )*polcolgq
         res(21,10) = diag(21) * dipole(5) * 
     &      ( pdf(1*sign1,1,10)*pdf(-1*sign2,2,10) + 
     &        pdf(3*sign1,1,10)*pdf(-3*sign2,2,10) )*polcolqq
         res(22,10) = diag(22) * dipole(5) * 
     &      ( pdf(1*sign1,1,10)*pdf(-1*sign2,2,10) + 
     &        pdf(3*sign1,1,10)*pdf(-3*sign2,2,10) + 
     &        lbquark*pdf(5*sign1,1,10)*pdf(-5*sign2,2,10) )*polcolqq
         res(23,10) = diag(23) * dipole(5) * 
     &      ( pdf(1*sign1,1,10)*pdf(-1*sign2,2,10) + 
     &        pdf(3*sign1,1,10)*pdf(-3*sign2,2,10) + 
     &        lbquark*pdf(5*sign1,1,10)*pdf(-5*sign2,2,10) + 
     &        lbquark*pdf(5*sign1,1,10)*pdf(-5*sign2,2,10) )*polcolqq
         res(24,10) = diag(24) * dipole(5) * 
     &      ( pdf(1*sign1,1,10)*pdf(-1*sign2,2,10) + 
     &        lbquark*pdf(1*sign1,1,10)*pdf(-1*sign2,2,10) + 
     &        pdf(3*sign1,1,10)*pdf(-3*sign2,2,10) + 
     &        lbquark*pdf(3*sign1,1,10)*pdf(-3*sign2,2,10) + 
     &        lbquark*pdf(5*sign1,1,10)*pdf(-5*sign2,2,10) + 
     &        lbquark*pdf(5*sign1,1,10)*pdf(-5*sign2,2,10) )*polcolqq
         res(25,10) = diag(25) * dipole(6) * 
     &      ( pdf(1*sign1,1,10)*pdf(-1*sign2,2,10) + 
     &        pdf(3*sign1,1,10)*pdf(-3*sign2,2,10) + 
     &        lbquark*pdf(5*sign1,1,10)*pdf(-5*sign2,2,10) )/2d0*polcolqq
         res(27,10) = diag(27) * dipole(7) * 
     &      ( pdf(0*sign1,1,10)*pdf(2*sign2,2,10) + 
     &        pdf(0*sign1,1,10)*pdf(4*sign2,2,10) )*polcolgq
         res(28,10) = diag(28) * dipole(8) * 
     &      ( pdf(0*sign1,1,10)*pdf(1*sign2,2,10) + 
     &        pdf(0*sign1,1,10)*pdf(3*sign2,2,10) + 
     &        lbquark*pdf(0*sign1,1,10)*pdf(5*sign2,2,10) )*polcolgq
         res(31,10) = diag(31) * dipole(9) * 
     &      ( pdf(0*sign1,1,10)*pdf(-1*sign2,2,10) + 
     &        pdf(0*sign1,1,10)*pdf(-3*sign2,2,10) + 
     &        lbquark*pdf(0*sign1,1,10)*pdf(-5*sign2,2,10) )*polcolgq
         res(32,10) = diag(32) * dipole(10) * 
     &      ( pdf(0*sign1,1,10)*pdf(-2*sign2,2,10) + 
     &        pdf(0*sign1,1,10)*pdf(-4*sign2,2,10) )*polcolgq
         res(36,10) = diag(36) * dipole(11) * 
     &      ( pdf(-1*sign1,1,10)*pdf(1*sign2,2,10) + 
     &        pdf(-3*sign1,1,10)*pdf(3*sign2,2,10) + 
     &        lbquark*pdf(-5*sign1,1,10)*pdf(5*sign2,2,10) + 
     &        lbquark*pdf(-5*sign1,1,10)*pdf(5*sign2,2,10) )*polcolqq
         res(37,10) = diag(37) * dipole(11) * 
     &      ( pdf(-1*sign1,1,10)*pdf(1*sign2,2,10) + 
     &        lbquark*pdf(-1*sign1,1,10)*pdf(1*sign2,2,10) + 
     &        pdf(-3*sign1,1,10)*pdf(3*sign2,2,10) + 
     &        lbquark*pdf(-3*sign1,1,10)*pdf(3*sign2,2,10) + 
     &        lbquark*pdf(-5*sign1,1,10)*pdf(5*sign2,2,10) + 
     &        lbquark*pdf(-5*sign1,1,10)*pdf(5*sign2,2,10) )*polcolqq
         res(38,10) = diag(38) * dipole(11) * 
     &      ( pdf(-1*sign1,1,10)*pdf(1*sign2,2,10) + 
     &        pdf(-3*sign1,1,10)*pdf(3*sign2,2,10) )*polcolqq
         res(39,10) = diag(39) * dipole(11) * 
     &      ( pdf(-1*sign1,1,10)*pdf(1*sign2,2,10) + 
     &        pdf(-3*sign1,1,10)*pdf(3*sign2,2,10) + 
     &        lbquark*pdf(-5*sign1,1,10)*pdf(5*sign2,2,10) )*polcolqq
         res(40,10) = diag(40) * dipole(12) * 
     &      ( pdf(-1*sign1,1,10)*pdf(1*sign2,2,10) + 
     &        pdf(-3*sign1,1,10)*pdf(3*sign2,2,10) + 
     &        lbquark*pdf(-5*sign1,1,10)*pdf(5*sign2,2,10) )/2d0*polcolqq
         res(41,10) = diag(41) * dipole(13) * 
     &      ( pdf(-1*sign1,1,10)*pdf(0*sign2,2,10) + 
     &        pdf(-3*sign1,1,10)*pdf(0*sign2,2,10) + 
     &        lbquark*pdf(-5*sign1,1,10)*pdf(0*sign2,2,10) )*polcolgq
         res(48,10) = diag(48) * dipole(14) * 
     &      ( pdf(-2*sign1,1,10)*pdf(2*sign2,2,10) + 
     &        pdf(-4*sign1,1,10)*pdf(4*sign2,2,10) )*polcolqq
         res(49,10) = diag(49) * dipole(14) * 
     &      ( pdf(-2*sign1,1,10)*pdf(2*sign2,2,10) + 
     &        lbquark*pdf(-2*sign1,1,10)*pdf(2*sign2,2,10) + 
     &        pdf(-4*sign1,1,10)*pdf(4*sign2,2,10) + 
     &        lbquark*pdf(-4*sign1,1,10)*pdf(4*sign2,2,10) )*polcolqq
         res(50,10) = diag(50) * dipole(14) * 
     &      ( pdf(-2*sign1,1,10)*pdf(2*sign2,2,10) + 
     &        pdf(-4*sign1,1,10)*pdf(4*sign2,2,10) )*polcolqq
         res(51,10) = diag(51) * dipole(14) * 
     &      ( pdf(-2*sign1,1,10)*pdf(2*sign2,2,10) + 
     &        pdf(-4*sign1,1,10)*pdf(4*sign2,2,10) )*polcolqq
         res(52,10) = diag(52) * dipole(15) * 
     &      ( pdf(-2*sign1,1,10)*pdf(2*sign2,2,10) + 
     &        pdf(-4*sign1,1,10)*pdf(4*sign2,2,10) )/2d0*polcolqq
         res(54,10) = diag(54) * dipole(16) * 
     &      ( pdf(-2*sign1,1,10)*pdf(0*sign2,2,10) + 
     &        pdf(-4*sign1,1,10)*pdf(0*sign2,2,10) )*polcolgq
       endif 

       if (lok(11)) then 
         call qqgZAZA_curr(borncurr(0,-1,1,51),11,p(0,3,11),p(0,1,11),p(0,2,11),-1,-1,-1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,52),11,p(0,2,11),p(0,1,11),p(0,3,11),1,-1,1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,53),11,p(0,3,11),p(0,2,11),p(0,1,11),-1,-1,-1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,54),11,p(0,2,11),p(0,3,11),p(0,1,11),1,1,-1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,55),11,p(0,1,11),p(0,2,11),p(0,3,11),1,-1,1,bos)
         call qqgZAZA_curr(borncurr(0,-1,1,56),11,p(0,1,11),p(0,3,11),p(0,2,11),1,1,-1,bos)

         dipole(1) = - dipoleVVj_CS39(11,0,p(0,3,1),p(0,4,1),p(0,2,1),xuz(1,11),borncurr(0,0,1,56))
         dipole(2) = - dipoleVVj_CS41(11,1,p(0,3,1),p(0,4,1),p(0,2,1),xuz(1,11),borncurr(0,-1,1,55))
         dipole(3) = - dipoleVVj_CS40(11,1,p(0,3,1),p(0,4,1),p(0,2,1),xuz(1,11),borncurr(0,-1,1,55))
         dipole(4) = - dipoleVVj_CS39(11,0,p(0,3,1),p(0,4,1),p(0,2,1),xuz(1,11),borncurr(0,0,2,56))
         dipole(5) = - dipoleVVj_CS41(11,1,p(0,3,1),p(0,4,1),p(0,2,1),xuz(1,11),borncurr(0,-1,2,55))
         dipole(6) = - dipoleVVj_CS40(11,1,p(0,3,1),p(0,4,1),p(0,2,1),xuz(1,11),borncurr(0,-1,2,55))
         dipole(7) = - dipoleVVj_CS39(11,1,p(0,3,1),p(0,4,1),p(0,2,1),xuz(1,11),borncurr(0,0,1,54))
         dipole(8) = - dipoleVVj_CS39(11,1,p(0,3,1),p(0,4,1),p(0,2,1),xuz(1,11),borncurr(0,0,2,54))
         dipole(9) = - dipoleVVj_CS39(11,1,p(0,3,1),p(0,4,1),p(0,2,1),xuz(1,11),borncurr(0,0,2,53))
         dipole(10) = - dipoleVVj_CS39(11,1,p(0,3,1),p(0,4,1),p(0,2,1),xuz(1,11),borncurr(0,0,1,53))
         dipole(11) = - dipoleVVj_CS41(11,1,p(0,3,1),p(0,4,1),p(0,2,1),xuz(1,11),borncurr(0,-1,2,52))
         dipole(12) = - dipoleVVj_CS40(11,1,p(0,3,1),p(0,4,1),p(0,2,1),xuz(1,11),borncurr(0,-1,2,52))
         dipole(13) = - dipoleVVj_CS39(11,0,p(0,3,1),p(0,4,1),p(0,2,1),xuz(1,11),borncurr(0,0,2,51))
         dipole(14) = - dipoleVVj_CS41(11,1,p(0,3,1),p(0,4,1),p(0,2,1),xuz(1,11),borncurr(0,-1,1,52))
         dipole(15) = - dipoleVVj_CS40(11,1,p(0,3,1),p(0,4,1),p(0,2,1),xuz(1,11),borncurr(0,-1,1,52))
         dipole(16) = - dipoleVVj_CS39(11,0,p(0,3,1),p(0,4,1),p(0,2,1),xuz(1,11),borncurr(0,0,1,51))

         res(5,11) = diag(5) * dipole(1) * 
     &      ( pdf(2*sign1,1,11)*pdf(0*sign2,2,11) + 
     &        pdf(4*sign1,1,11)*pdf(0*sign2,2,11) )*polcolgq
         res(9,11) = diag(9) * dipole(2) * 
     &      ( pdf(2*sign1,1,11)*pdf(-2*sign2,2,11) + 
     &        pdf(4*sign1,1,11)*pdf(-4*sign2,2,11) )*polcolqq
         res(10,11) = diag(10) * dipole(2) * 
     &      ( pdf(2*sign1,1,11)*pdf(-2*sign2,2,11) + 
     &        pdf(4*sign1,1,11)*pdf(-4*sign2,2,11) )*polcolqq
         res(11,11) = diag(11) * dipole(2) * 
     &      ( pdf(2*sign1,1,11)*pdf(-2*sign2,2,11) + 
     &        pdf(4*sign1,1,11)*pdf(-4*sign2,2,11) )*polcolqq
         res(12,11) = diag(12) * dipole(2) * 
     &      ( pdf(2*sign1,1,11)*pdf(-2*sign2,2,11) + 
     &        lbquark*pdf(2*sign1,1,11)*pdf(-2*sign2,2,11) + 
     &        pdf(4*sign1,1,11)*pdf(-4*sign2,2,11) + 
     &        lbquark*pdf(4*sign1,1,11)*pdf(-4*sign2,2,11) )*polcolqq
         res(13,11) = diag(13) * dipole(3) * 
     &      ( pdf(2*sign1,1,11)*pdf(-2*sign2,2,11) + 
     &        pdf(4*sign1,1,11)*pdf(-4*sign2,2,11) )/2d0*polcolqq
         res(18,11) = diag(18) * dipole(4) * 
     &      ( pdf(1*sign1,1,11)*pdf(0*sign2,2,11) + 
     &        pdf(3*sign1,1,11)*pdf(0*sign2,2,11) + 
     &        lbquark*pdf(5*sign1,1,11)*pdf(0*sign2,2,11) )*polcolgq
         res(21,11) = diag(21) * dipole(5) * 
     &      ( pdf(1*sign1,1,11)*pdf(-1*sign2,2,11) + 
     &        pdf(3*sign1,1,11)*pdf(-3*sign2,2,11) )*polcolqq
         res(22,11) = diag(22) * dipole(5) * 
     &      ( pdf(1*sign1,1,11)*pdf(-1*sign2,2,11) + 
     &        pdf(3*sign1,1,11)*pdf(-3*sign2,2,11) + 
     &        lbquark*pdf(5*sign1,1,11)*pdf(-5*sign2,2,11) )*polcolqq
         res(23,11) = diag(23) * dipole(5) * 
     &      ( pdf(1*sign1,1,11)*pdf(-1*sign2,2,11) + 
     &        pdf(3*sign1,1,11)*pdf(-3*sign2,2,11) + 
     &        lbquark*pdf(5*sign1,1,11)*pdf(-5*sign2,2,11) + 
     &        lbquark*pdf(5*sign1,1,11)*pdf(-5*sign2,2,11) )*polcolqq
         res(24,11) = diag(24) * dipole(5) * 
     &      ( pdf(1*sign1,1,11)*pdf(-1*sign2,2,11) + 
     &        lbquark*pdf(1*sign1,1,11)*pdf(-1*sign2,2,11) + 
     &        pdf(3*sign1,1,11)*pdf(-3*sign2,2,11) + 
     &        lbquark*pdf(3*sign1,1,11)*pdf(-3*sign2,2,11) + 
     &        lbquark*pdf(5*sign1,1,11)*pdf(-5*sign2,2,11) + 
     &        lbquark*pdf(5*sign1,1,11)*pdf(-5*sign2,2,11) )*polcolqq
         res(25,11) = diag(25) * dipole(6) * 
     &      ( pdf(1*sign1,1,11)*pdf(-1*sign2,2,11) + 
     &        pdf(3*sign1,1,11)*pdf(-3*sign2,2,11) + 
     &        lbquark*pdf(5*sign1,1,11)*pdf(-5*sign2,2,11) )/2d0*polcolqq
         res(27,11) = diag(27) * dipole(7) * 
     &      ( pdf(0*sign1,1,11)*pdf(2*sign2,2,11) + 
     &        pdf(0*sign1,1,11)*pdf(4*sign2,2,11) )*polcolgq
         res(28,11) = diag(28) * dipole(8) * 
     &      ( pdf(0*sign1,1,11)*pdf(1*sign2,2,11) + 
     &        pdf(0*sign1,1,11)*pdf(3*sign2,2,11) + 
     &        lbquark*pdf(0*sign1,1,11)*pdf(5*sign2,2,11) )*polcolgq
         res(31,11) = diag(31) * dipole(9) * 
     &      ( pdf(0*sign1,1,11)*pdf(-1*sign2,2,11) + 
     &        pdf(0*sign1,1,11)*pdf(-3*sign2,2,11) + 
     &        lbquark*pdf(0*sign1,1,11)*pdf(-5*sign2,2,11) )*polcolgq
         res(32,11) = diag(32) * dipole(10) * 
     &      ( pdf(0*sign1,1,11)*pdf(-2*sign2,2,11) + 
     &        pdf(0*sign1,1,11)*pdf(-4*sign2,2,11) )*polcolgq
         res(36,11) = diag(36) * dipole(11) * 
     &      ( pdf(-1*sign1,1,11)*pdf(1*sign2,2,11) + 
     &        pdf(-3*sign1,1,11)*pdf(3*sign2,2,11) + 
     &        lbquark*pdf(-5*sign1,1,11)*pdf(5*sign2,2,11) + 
     &        lbquark*pdf(-5*sign1,1,11)*pdf(5*sign2,2,11) )*polcolqq
         res(37,11) = diag(37) * dipole(11) * 
     &      ( pdf(-1*sign1,1,11)*pdf(1*sign2,2,11) + 
     &        lbquark*pdf(-1*sign1,1,11)*pdf(1*sign2,2,11) + 
     &        pdf(-3*sign1,1,11)*pdf(3*sign2,2,11) + 
     &        lbquark*pdf(-3*sign1,1,11)*pdf(3*sign2,2,11) + 
     &        lbquark*pdf(-5*sign1,1,11)*pdf(5*sign2,2,11) + 
     &        lbquark*pdf(-5*sign1,1,11)*pdf(5*sign2,2,11) )*polcolqq
         res(38,11) = diag(38) * dipole(11) * 
     &      ( pdf(-1*sign1,1,11)*pdf(1*sign2,2,11) + 
     &        pdf(-3*sign1,1,11)*pdf(3*sign2,2,11) )*polcolqq
         res(39,11) = diag(39) * dipole(11) * 
     &      ( pdf(-1*sign1,1,11)*pdf(1*sign2,2,11) + 
     &        pdf(-3*sign1,1,11)*pdf(3*sign2,2,11) + 
     &        lbquark*pdf(-5*sign1,1,11)*pdf(5*sign2,2,11) )*polcolqq
         res(40,11) = diag(40) * dipole(12) * 
     &      ( pdf(-1*sign1,1,11)*pdf(1*sign2,2,11) + 
     &        pdf(-3*sign1,1,11)*pdf(3*sign2,2,11) + 
     &        lbquark*pdf(-5*sign1,1,11)*pdf(5*sign2,2,11) )/2d0*polcolqq
         res(41,11) = diag(41) * dipole(13) * 
     &      ( pdf(-1*sign1,1,11)*pdf(0*sign2,2,11) + 
     &        pdf(-3*sign1,1,11)*pdf(0*sign2,2,11) + 
     &        lbquark*pdf(-5*sign1,1,11)*pdf(0*sign2,2,11) )*polcolgq
         res(48,11) = diag(48) * dipole(14) * 
     &      ( pdf(-2*sign1,1,11)*pdf(2*sign2,2,11) + 
     &        pdf(-4*sign1,1,11)*pdf(4*sign2,2,11) )*polcolqq
         res(49,11) = diag(49) * dipole(14) * 
     &      ( pdf(-2*sign1,1,11)*pdf(2*sign2,2,11) + 
     &        lbquark*pdf(-2*sign1,1,11)*pdf(2*sign2,2,11) + 
     &        pdf(-4*sign1,1,11)*pdf(4*sign2,2,11) + 
     &        lbquark*pdf(-4*sign1,1,11)*pdf(4*sign2,2,11) )*polcolqq
         res(50,11) = diag(50) * dipole(14) * 
     &      ( pdf(-2*sign1,1,11)*pdf(2*sign2,2,11) + 
     &        pdf(-4*sign1,1,11)*pdf(4*sign2,2,11) )*polcolqq
         res(51,11) = diag(51) * dipole(14) * 
     &      ( pdf(-2*sign1,1,11)*pdf(2*sign2,2,11) + 
     &        pdf(-4*sign1,1,11)*pdf(4*sign2,2,11) )*polcolqq
         res(52,11) = diag(52) * dipole(15) * 
     &      ( pdf(-2*sign1,1,11)*pdf(2*sign2,2,11) + 
     &        pdf(-4*sign1,1,11)*pdf(4*sign2,2,11) )/2d0*polcolqq
         res(54,11) = diag(54) * dipole(16) * 
     &      ( pdf(-2*sign1,1,11)*pdf(0*sign2,2,11) + 
     &        pdf(-4*sign1,1,11)*pdf(0*sign2,2,11) )*polcolgq
       endif 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c finite collinears
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c put everything in entry 59, split equally between 4/5 and 8/9

       fincolglobal = 2d0*pi/dotrr(p(0,1,1),p(0,2,1))

       if (lok(4).and.lfincolls) then 
         q_sf = sqrt(mufsq(1,4))
         fincolfac = fincolglobal*als(1,4)/xuz(1,4)*xi(1)

         call finitecollgqREM(xi(1),xuz(1,4),p(0,1,4),p(0,2,4),p(0,3,4),q_sf,fincol(1),sign1)
         call finitecollqgREM(xi(1),xuz(1,4),p(0,1,4),p(0,2,4),p(0,3,4),q_sf,2,fincol(2),sign1)
         call finitecollqgREM(xi(1),xuz(1,4),p(0,1,4),p(0,2,4),p(0,3,4),q_sf,4,fincol(3),sign1)
         call finitecollqqREM(xi(1),xuz(1,4),p(0,1,4),p(0,2,4),p(0,3,4),q_sf,2,fincol(4),sign1)
         call finitecollqqREM(xi(1),xuz(1,4),p(0,1,4),p(0,2,4),p(0,3,4),q_sf,4,fincol(5),sign1)
         call finitecollqgREM(xi(1),xuz(1,4),p(0,1,4),p(0,2,4),p(0,3,4),q_sf,1,fincol(6),sign1)
         call finitecollqgREM(xi(1),xuz(1,4),p(0,1,4),p(0,2,4),p(0,3,4),q_sf,3,fincol(7),sign1)
         call finitecollqqREM(xi(1),xuz(1,4),p(0,1,4),p(0,2,4),p(0,3,4),q_sf,1,fincol(8),sign1)
         call finitecollqqREM(xi(1),xuz(1,4),p(0,1,4),p(0,2,4),p(0,3,4),q_sf,3,fincol(9),sign1)
         call finitecollqqREM(xi(1),xuz(1,4),p(0,1,4),p(0,2,4),p(0,3,4),q_sf,-1,fincol(10),sign1)
         call finitecollqqREM(xi(1),xuz(1,4),p(0,1,4),p(0,2,4),p(0,3,4),q_sf,-3,fincol(11),sign1)
         call finitecollqgREM(xi(1),xuz(1,4),p(0,1,4),p(0,2,4),p(0,3,4),q_sf,-1,fincol(12),sign1)
         call finitecollqgREM(xi(1),xuz(1,4),p(0,1,4),p(0,2,4),p(0,3,4),q_sf,-3,fincol(13),sign1)
         call finitecollqqREM(xi(1),xuz(1,4),p(0,1,4),p(0,2,4),p(0,3,4),q_sf,-2,fincol(14),sign1)
         call finitecollqqREM(xi(1),xuz(1,4),p(0,1,4),p(0,2,4),p(0,3,4),q_sf,-4,fincol(15),sign1)
         call finitecollqgREM(xi(1),xuz(1,4),p(0,1,4),p(0,2,4),p(0,3,4),q_sf,-2,fincol(16),sign1)
         call finitecollqgREM(xi(1),xuz(1,4),p(0,1,4),p(0,2,4),p(0,3,4),q_sf,-4,fincol(17),sign1)
         if (lbquark .ne. 0) then
           call finitecollqgREM(xi(1),xuz(1,4),p(0,1,4),p(0,2,4),p(0,3,4),q_sf,5,fincol(18),sign1)
           call finitecollqqREM(xi(1),xuz(1,4),p(0,1,4),p(0,2,4),p(0,3,4),q_sf,5,fincol(19),sign1)
           call finitecollqgREM(xi(1),xuz(1,4),p(0,1,4),p(0,2,4),p(0,3,4),q_sf,-5,fincol(20),sign1)
           call finitecollqqREM(xi(1),xuz(1,4),p(0,1,4),p(0,2,4),p(0,3,4),q_sf,-5,fincol(21),sign1)
         endif

         res(59,4) = diag(59) * fincolfac * (
c g u/c: 16 1
     &     + borncurr(0,0,1,16) * fincol(1) *
     &       (pdf(2*sign2,2,4) + pdf(4*sign2,2,4))*polcolgq
c g d/s/b: 16 1
     &     + borncurr(0,0,2,16) * fincol(1) *
     &       (pdf(1*sign2,2,4) + pdf(3*sign2,2,4) + lbquark*pdf(5*sign2,2,4))*polcolgq
c g U/C: 15 1
     &     + borncurr(0,0,1,15) * fincol(1) *
     &       (pdf(-2*sign2,2,4) + pdf(-4*sign2,2,4))*polcolgq
c g D/S/B: 15 1
     &     + borncurr(0,0,2,15) * fincol(1) *
     &       (pdf(-1*sign2,2,4) + pdf(-3*sign2,2,4) + lbquark*pdf(-5*sign2,2,4))*polcolgq
c u g: 18 2
     &     + borncurr(0,0,1,18) * fincol(2) *
     &       pdf(0*sign2,2,4)*polcolgq
c c g: 18 3
     &     + borncurr(0,0,1,18) * fincol(3) *
     &       pdf(0*sign2,2,4)*polcolgq
c d g: 18 6
     &     + borncurr(0,0,2,18) * fincol(6) *
     &       pdf(0*sign2,2,4)*polcolgq
c s g: 18 7
     &     + borncurr(0,0,2,18) * fincol(7) *
     &       pdf(0*sign2,2,4)*polcolgq
c b g: 18 18
     &     + borncurr(0,0,2,18) * fincol(18) *
     &       lbquark*pdf(0*sign2,2,4)*polcolgq
c D g: 13 12
     &     + borncurr(0,0,2,13) * fincol(12) *
     &       pdf(0*sign2,2,4)*polcolgq
c S g: 13 13
     &     + borncurr(0,0,2,13) * fincol(13) *
     &       pdf(0*sign2,2,4)*polcolgq
c B g: 13 20
     &     + borncurr(0,0,2,13) * fincol(20) *
     &       lbquark*pdf(0*sign2,2,4)*polcolgq
c U g: 13 16
     &     + borncurr(0,0,1,13) * fincol(16) *
     &       pdf(0*sign2,2,4)*polcolgq
c C g: 13 17
     &     + borncurr(0,0,1,13) * fincol(17) *
     &       pdf(0*sign2,2,4)*polcolgq
c u U: 17 4
     &     + borncurr(0,0,1,17) * fincol(4) *
     &       pdf(-2*sign2,2,4)*polcolqq
c c C: 17 5
     &     + borncurr(0,0,1,17) * fincol(5) *
     &       pdf(-4*sign2,2,4)*polcolqq
c d D: 17 8
     &     + borncurr(0,0,2,17) * fincol(8) *
     &       pdf(-1*sign2,2,4)*polcolqq
c s S: 17 9
     &     + borncurr(0,0,2,17) * fincol(9) *
     &       pdf(-3*sign2,2,4)*polcolqq
c b B: 17 19
     &     + borncurr(0,0,2,17) * fincol(19) *
     &       lbquark*pdf(-5*sign2,2,4)*polcolqq
c D d: 14 10
     &     + borncurr(0,0,2,14) * fincol(10) *
     &       pdf(1*sign2,2,4)*polcolqq
c S s: 14 11
     &     + borncurr(0,0,2,14) * fincol(11) *
     &       pdf(3*sign2,2,4)*polcolqq
c B b: 14 21
     &     + borncurr(0,0,2,14) * fincol(21) *
     &       lbquark*pdf(5*sign2,2,4)*polcolqq
c U u: 14 14
     &     + borncurr(0,0,1,14) * fincol(14) *
     &       pdf(2*sign2,2,4)*polcolqq
c C c: 14 15
     &     + borncurr(0,0,1,14) * fincol(15) *
     &       pdf(4*sign2,2,4)*polcolqq
     &   )
       endif 

       if (lok(5).and.lfincolls) then 
         q_sf = sqrt(mufsq(2,5))
         fincolfac = fincolglobal*als(1,5)/xuz(1,5)*xi(2)

         call finitecollgqREM(xi(2),xuz(1,5),p(0,2,5),p(0,1,5),p(0,3,5),q_sf,fincol(1),sign2)
         call finitecollqgREM(xi(2),xuz(1,5),p(0,2,5),p(0,1,5),p(0,3,5),q_sf,2,fincol(2),sign2)
         call finitecollqgREM(xi(2),xuz(1,5),p(0,2,5),p(0,1,5),p(0,3,5),q_sf,4,fincol(3),sign2)
         call finitecollqqREM(xi(2),xuz(1,5),p(0,2,5),p(0,1,5),p(0,3,5),q_sf,2,fincol(4),sign2)
         call finitecollqqREM(xi(2),xuz(1,5),p(0,2,5),p(0,1,5),p(0,3,5),q_sf,4,fincol(5),sign2)
         call finitecollqgREM(xi(2),xuz(1,5),p(0,2,5),p(0,1,5),p(0,3,5),q_sf,1,fincol(6),sign2)
         call finitecollqgREM(xi(2),xuz(1,5),p(0,2,5),p(0,1,5),p(0,3,5),q_sf,3,fincol(7),sign2)
         call finitecollqqREM(xi(2),xuz(1,5),p(0,2,5),p(0,1,5),p(0,3,5),q_sf,1,fincol(8),sign2)
         call finitecollqqREM(xi(2),xuz(1,5),p(0,2,5),p(0,1,5),p(0,3,5),q_sf,3,fincol(9),sign2)
         call finitecollqqREM(xi(2),xuz(1,5),p(0,2,5),p(0,1,5),p(0,3,5),q_sf,-1,fincol(10),sign2)
         call finitecollqqREM(xi(2),xuz(1,5),p(0,2,5),p(0,1,5),p(0,3,5),q_sf,-3,fincol(11),sign2)
         call finitecollqgREM(xi(2),xuz(1,5),p(0,2,5),p(0,1,5),p(0,3,5),q_sf,-1,fincol(12),sign2)
         call finitecollqgREM(xi(2),xuz(1,5),p(0,2,5),p(0,1,5),p(0,3,5),q_sf,-3,fincol(13),sign2)
         call finitecollqqREM(xi(2),xuz(1,5),p(0,2,5),p(0,1,5),p(0,3,5),q_sf,-2,fincol(14),sign2)
         call finitecollqqREM(xi(2),xuz(1,5),p(0,2,5),p(0,1,5),p(0,3,5),q_sf,-4,fincol(15),sign2)
         call finitecollqgREM(xi(2),xuz(1,5),p(0,2,5),p(0,1,5),p(0,3,5),q_sf,-2,fincol(16),sign2)
         call finitecollqgREM(xi(2),xuz(1,5),p(0,2,5),p(0,1,5),p(0,3,5),q_sf,-4,fincol(17),sign2)
         if (lbquark .ne. 0) then
           call finitecollqgREM(xi(2),xuz(1,5),p(0,2,5),p(0,1,5),p(0,3,5),q_sf,5,fincol(18),sign2)
           call finitecollqqREM(xi(2),xuz(1,5),p(0,2,5),p(0,1,5),p(0,3,5),q_sf,5,fincol(19),sign2)
           call finitecollqgREM(xi(2),xuz(1,5),p(0,2,5),p(0,1,5),p(0,3,5),q_sf,-5,fincol(20),sign2)
           call finitecollqqREM(xi(2),xuz(1,5),p(0,2,5),p(0,1,5),p(0,3,5),q_sf,-5,fincol(21),sign2)
         endif

         res(59,5) = diag(59) * fincolfac * (
c u/c g: 24 1
     &     + borncurr(0,0,1,24) * fincol(1) *
     &       (pdf(2*sign1,1,5) + pdf(4*sign1,1,5))*polcolgq
c d/s/b g: 24 1
     &     + borncurr(0,0,2,24) * fincol(1) *
     &       (pdf(1*sign1,1,5) + pdf(3*sign1,1,5) + lbquark*pdf(5*sign1,1,5))*polcolgq
c U/C g: 19 1
     &     + borncurr(0,0,1,19) * fincol(1) *
     &       (pdf(-2*sign1,1,5) + pdf(-4*sign1,1,5))*polcolgq
c D/S/B g: 19 1
     &     + borncurr(0,0,2,19) * fincol(1) *
     &       (pdf(-1*sign1,1,5) + pdf(-3*sign1,1,5) + lbquark*pdf(-5*sign1,1,5))*polcolgq
c g u: 22 2
     &     + borncurr(0,0,1,22) * fincol(2) *
     &       pdf(0*sign1,1,5)*polcolgq
c g c: 22 3
     &     + borncurr(0,0,1,22) * fincol(3) *
     &       pdf(0*sign1,1,5)*polcolgq
c g d: 22 6
     &     + borncurr(0,0,2,22) * fincol(6) *
     &       pdf(0*sign1,1,5)*polcolgq
c g s: 22 7
     &     + borncurr(0,0,2,22) * fincol(7) *
     &       pdf(0*sign1,1,5)*polcolgq
c g b: 22 18
     &     + borncurr(0,0,2,22) * fincol(18) *
     &       lbquark*pdf(0*sign1,1,5)*polcolgq
c g D: 21 12
     &     + borncurr(0,0,2,21) * fincol(12) *
     &       pdf(0*sign1,1,5)*polcolgq
c g S: 21 13
     &     + borncurr(0,0,2,21) * fincol(13) *
     &       pdf(0*sign1,1,5)*polcolgq
c g B: 21 20
     &     + borncurr(0,0,2,21) * fincol(20) *
     &       lbquark*pdf(0*sign1,1,5)*polcolgq
c g U: 21 16
     &     + borncurr(0,0,1,21) * fincol(16) *
     &       pdf(0*sign1,1,5)*polcolgq
c g C: 21 17
     &     + borncurr(0,0,1,21) * fincol(17) *
     &       pdf(0*sign1,1,5)*polcolgq
c U u: 20 4
     &     + borncurr(0,0,1,20) * fincol(4) *
     &       pdf(-2*sign1,1,5)*polcolqq
c C c: 20 5
     &     + borncurr(0,0,1,20) * fincol(5) *
     &       pdf(-4*sign1,1,5)*polcolqq
c D d: 20 8
     &     + borncurr(0,0,2,20) * fincol(8) *
     &       pdf(-1*sign1,1,5)*polcolqq
c S s: 20 9
     &     + borncurr(0,0,2,20) * fincol(9) *
     &       pdf(-3*sign1,1,5)*polcolqq
c B b: 20 19
     &     + borncurr(0,0,2,20) * fincol(19) *
     &       lbquark*pdf(-5*sign1,1,5)*polcolqq
c d D: 23 10
     &     + borncurr(0,0,2,23) * fincol(10) *
     &       pdf(1*sign1,1,5)*polcolqq
c s S: 23 11
     &     + borncurr(0,0,2,23) * fincol(11) *
     &       pdf(3*sign1,1,5)*polcolqq
c b B: 23 21
     &     + borncurr(0,0,2,23) * fincol(21) *
     &       lbquark*pdf(5*sign1,1,5)*polcolqq
c u U: 23 14
     &     + borncurr(0,0,1,23) * fincol(14) *
     &       pdf(2*sign1,1,5)*polcolqq
c c C: 23 15
     &     + borncurr(0,0,1,23) * fincol(15) *
     &       pdf(4*sign1,1,5)*polcolqq
     &   )
       endif 

       if (lok(8).and.lfincolls) then 
         q_sf = sqrt(mufsq(1,8))
         fincolfac = fincolglobal*als(1,8)/xuz(1,8)*xi(1)
         call qqgZAZA_curr(borncurr(0,-1,1,61),8,p(0,1,8),p(0,3,8),p(0,2,8),1,1,-1,bos)

         call finitecollgqREM(xi(1),xuz(1,8),p(0,1,8),p(0,2,8),p(0,3,8),q_sf,fincol(1),sign1)
         call finitecollqgREM(xi(1),xuz(1,8),p(0,1,8),p(0,2,8),p(0,3,8),q_sf,2,fincol(2),sign1)
         call finitecollqgREM(xi(1),xuz(1,8),p(0,1,8),p(0,2,8),p(0,3,8),q_sf,4,fincol(3),sign1)
         call finitecollqqREM(xi(1),xuz(1,8),p(0,1,8),p(0,2,8),p(0,3,8),q_sf,2,fincol(4),sign1)
         call finitecollqqREM(xi(1),xuz(1,8),p(0,1,8),p(0,2,8),p(0,3,8),q_sf,4,fincol(5),sign1)
         call finitecollqgREM(xi(1),xuz(1,8),p(0,1,8),p(0,2,8),p(0,3,8),q_sf,1,fincol(6),sign1)
         call finitecollqgREM(xi(1),xuz(1,8),p(0,1,8),p(0,2,8),p(0,3,8),q_sf,3,fincol(7),sign1)
         call finitecollqqREM(xi(1),xuz(1,8),p(0,1,8),p(0,2,8),p(0,3,8),q_sf,1,fincol(8),sign1)
         call finitecollqqREM(xi(1),xuz(1,8),p(0,1,8),p(0,2,8),p(0,3,8),q_sf,3,fincol(9),sign1)
         call finitecollqqREM(xi(1),xuz(1,8),p(0,1,8),p(0,2,8),p(0,3,8),q_sf,-1,fincol(10),sign1)
         call finitecollqqREM(xi(1),xuz(1,8),p(0,1,8),p(0,2,8),p(0,3,8),q_sf,-3,fincol(11),sign1)
         call finitecollqgREM(xi(1),xuz(1,8),p(0,1,8),p(0,2,8),p(0,3,8),q_sf,-1,fincol(12),sign1)
         call finitecollqgREM(xi(1),xuz(1,8),p(0,1,8),p(0,2,8),p(0,3,8),q_sf,-3,fincol(13),sign1)
         call finitecollqqREM(xi(1),xuz(1,8),p(0,1,8),p(0,2,8),p(0,3,8),q_sf,-2,fincol(14),sign1)
         call finitecollqqREM(xi(1),xuz(1,8),p(0,1,8),p(0,2,8),p(0,3,8),q_sf,-4,fincol(15),sign1)
         call finitecollqgREM(xi(1),xuz(1,8),p(0,1,8),p(0,2,8),p(0,3,8),q_sf,-2,fincol(16),sign1)
         call finitecollqgREM(xi(1),xuz(1,8),p(0,1,8),p(0,2,8),p(0,3,8),q_sf,-4,fincol(17),sign1)
         if (lbquark .ne. 0) then
           call finitecollqgREM(xi(1),xuz(1,8),p(0,1,8),p(0,2,8),p(0,3,8),q_sf,5,fincol(18),sign1)
           call finitecollqqREM(xi(1),xuz(1,8),p(0,1,8),p(0,2,8),p(0,3,8),q_sf,5,fincol(19),sign1)
           call finitecollqgREM(xi(1),xuz(1,8),p(0,1,8),p(0,2,8),p(0,3,8),q_sf,-5,fincol(20),sign1)
           call finitecollqqREM(xi(1),xuz(1,8),p(0,1,8),p(0,2,8),p(0,3,8),q_sf,-5,fincol(21),sign1)
         endif

         res(59,8) = diag(59) * fincolfac * (
c g u/c: 38 1
     &     + borncurr(0,0,1,38) * fincol(1) *
     &       (pdf(2*sign2,2,8) + pdf(4*sign2,2,8))*polcolgq
c g d/s/b: 38 1
     &     + borncurr(0,0,2,38) * fincol(1) *
     &       (pdf(1*sign2,2,8) + pdf(3*sign2,2,8) + lbquark*pdf(5*sign2,2,8))*polcolgq
c g U/C: 37 1
     &     + borncurr(0,0,1,37) * fincol(1) *
     &       (pdf(-2*sign2,2,8) + pdf(-4*sign2,2,8))*polcolgq
c g D/S/B: 37 1
     &     + borncurr(0,0,2,37) * fincol(1) *
     &       (pdf(-1*sign2,2,8) + pdf(-3*sign2,2,8) + lbquark*pdf(-5*sign2,2,8))*polcolgq
c u g: 61 2
     &     + borncurr(0,0,1,61) * fincol(2) *
     &       pdf(0*sign2,2,8)*polcolgq
c c g: 61 3
     &     + borncurr(0,0,1,61) * fincol(3) *
     &       pdf(0*sign2,2,8)*polcolgq
c d g: 61 6
     &     + borncurr(0,0,2,61) * fincol(6) *
     &       pdf(0*sign2,2,8)*polcolgq
c s g: 61 7
     &     + borncurr(0,0,2,61) * fincol(7) *
     &       pdf(0*sign2,2,8)*polcolgq
c b g: 61 18
     &     + borncurr(0,0,2,61) * fincol(18) *
     &       lbquark*pdf(0*sign2,2,8)*polcolgq
c D g: 35 12
     &     + borncurr(0,0,2,35) * fincol(12) *
     &       pdf(0*sign2,2,8)*polcolgq
c S g: 35 13
     &     + borncurr(0,0,2,35) * fincol(13) *
     &       pdf(0*sign2,2,8)*polcolgq
c B g: 35 20
     &     + borncurr(0,0,2,35) * fincol(20) *
     &       lbquark*pdf(0*sign2,2,8)*polcolgq
c U g: 35 16
     &     + borncurr(0,0,1,35) * fincol(16) *
     &       pdf(0*sign2,2,8)*polcolgq
c C g: 35 17
     &     + borncurr(0,0,1,35) * fincol(17) *
     &       pdf(0*sign2,2,8)*polcolgq
c u U: 39 4
     &     + borncurr(0,0,1,39) * fincol(4) *
     &       pdf(-2*sign2,2,8)*polcolqq
c c C: 39 5
     &     + borncurr(0,0,1,39) * fincol(5) *
     &       pdf(-4*sign2,2,8)*polcolqq
c d D: 39 8
     &     + borncurr(0,0,2,39) * fincol(8) *
     &       pdf(-1*sign2,2,8)*polcolqq
c s S: 39 9
     &     + borncurr(0,0,2,39) * fincol(9) *
     &       pdf(-3*sign2,2,8)*polcolqq
c b B: 39 19
     &     + borncurr(0,0,2,39) * fincol(19) *
     &       lbquark*pdf(-5*sign2,2,8)*polcolqq
c D d: 36 10
     &     + borncurr(0,0,2,36) * fincol(10) *
     &       pdf(1*sign2,2,8)*polcolqq
c S s: 36 11
     &     + borncurr(0,0,2,36) * fincol(11) *
     &       pdf(3*sign2,2,8)*polcolqq
c B b: 36 21
     &     + borncurr(0,0,2,36) * fincol(21) *
     &       lbquark*pdf(5*sign2,2,8)*polcolqq
c U u: 36 14
     &     + borncurr(0,0,1,36) * fincol(14) *
     &       pdf(2*sign2,2,8)*polcolqq
c C c: 36 15
     &     + borncurr(0,0,1,36) * fincol(15) *
     &       pdf(4*sign2,2,8)*polcolqq
     &   )
       endif 

       if (lok(9).and.lfincolls) then 
         q_sf = sqrt(mufsq(2,9))
         fincolfac = fincolglobal*als(1,9)/xuz(1,9)*xi(2)
         call qqgZAZA_curr(borncurr(0,-1,1,62),9,p(0,2,9),p(0,3,9),p(0,1,9),1,1,-1,bos)

         call finitecollgqREM(xi(2),xuz(1,9),p(0,2,9),p(0,1,9),p(0,3,9),q_sf,fincol(1),sign2)
         call finitecollqgREM(xi(2),xuz(1,9),p(0,2,9),p(0,1,9),p(0,3,9),q_sf,2,fincol(2),sign2)
         call finitecollqgREM(xi(2),xuz(1,9),p(0,2,9),p(0,1,9),p(0,3,9),q_sf,4,fincol(3),sign2)
         call finitecollqqREM(xi(2),xuz(1,9),p(0,2,9),p(0,1,9),p(0,3,9),q_sf,2,fincol(4),sign2)
         call finitecollqqREM(xi(2),xuz(1,9),p(0,2,9),p(0,1,9),p(0,3,9),q_sf,4,fincol(5),sign2)
         call finitecollqgREM(xi(2),xuz(1,9),p(0,2,9),p(0,1,9),p(0,3,9),q_sf,1,fincol(6),sign2)
         call finitecollqgREM(xi(2),xuz(1,9),p(0,2,9),p(0,1,9),p(0,3,9),q_sf,3,fincol(7),sign2)
         call finitecollqqREM(xi(2),xuz(1,9),p(0,2,9),p(0,1,9),p(0,3,9),q_sf,1,fincol(8),sign2)
         call finitecollqqREM(xi(2),xuz(1,9),p(0,2,9),p(0,1,9),p(0,3,9),q_sf,3,fincol(9),sign2)
         call finitecollqqREM(xi(2),xuz(1,9),p(0,2,9),p(0,1,9),p(0,3,9),q_sf,-1,fincol(10),sign2)
         call finitecollqqREM(xi(2),xuz(1,9),p(0,2,9),p(0,1,9),p(0,3,9),q_sf,-3,fincol(11),sign2)
         call finitecollqgREM(xi(2),xuz(1,9),p(0,2,9),p(0,1,9),p(0,3,9),q_sf,-1,fincol(12),sign2)
         call finitecollqgREM(xi(2),xuz(1,9),p(0,2,9),p(0,1,9),p(0,3,9),q_sf,-3,fincol(13),sign2)
         call finitecollqqREM(xi(2),xuz(1,9),p(0,2,9),p(0,1,9),p(0,3,9),q_sf,-2,fincol(14),sign2)
         call finitecollqqREM(xi(2),xuz(1,9),p(0,2,9),p(0,1,9),p(0,3,9),q_sf,-4,fincol(15),sign2)
         call finitecollqgREM(xi(2),xuz(1,9),p(0,2,9),p(0,1,9),p(0,3,9),q_sf,-2,fincol(16),sign2)
         call finitecollqgREM(xi(2),xuz(1,9),p(0,2,9),p(0,1,9),p(0,3,9),q_sf,-4,fincol(17),sign2)
         if (lbquark .ne. 0) then
           call finitecollqgREM(xi(2),xuz(1,9),p(0,2,9),p(0,1,9),p(0,3,9),q_sf,5,fincol(18),sign2)
           call finitecollqqREM(xi(2),xuz(1,9),p(0,2,9),p(0,1,9),p(0,3,9),q_sf,5,fincol(19),sign2)
           call finitecollqgREM(xi(2),xuz(1,9),p(0,2,9),p(0,1,9),p(0,3,9),q_sf,-5,fincol(20),sign2)
           call finitecollqqREM(xi(2),xuz(1,9),p(0,2,9),p(0,1,9),p(0,3,9),q_sf,-5,fincol(21),sign2)
         endif

         res(59,9) = diag(59) * fincolfac * (
c u/c g: 44 1
     &     + borncurr(0,0,1,44) * fincol(1) *
     &       (pdf(2*sign1,1,9) + pdf(4*sign1,1,9))*polcolgq
c d/s/b g: 44 1
     &     + borncurr(0,0,2,44) * fincol(1) *
     &       (pdf(1*sign1,1,9) + pdf(3*sign1,1,9) + lbquark*pdf(5*sign1,1,9))*polcolgq
c U/C g: 40 1
     &     + borncurr(0,0,1,40) * fincol(1) *
     &       (pdf(-2*sign1,1,9) + pdf(-4*sign1,1,9))*polcolgq
c D/S/B g: 40 1
     &     + borncurr(0,0,2,40) * fincol(1) *
     &       (pdf(-1*sign1,1,9) + pdf(-3*sign1,1,9) + lbquark*pdf(-5*sign1,1,9))*polcolgq
c g u: 62 2
     &     + borncurr(0,0,1,62) * fincol(2) *
     &       pdf(0*sign1,1,9)*polcolgq
c g c: 62 3
     &     + borncurr(0,0,1,62) * fincol(3) *
     &       pdf(0*sign1,1,9)*polcolgq
c g d: 62 6
     &     + borncurr(0,0,2,62) * fincol(6) *
     &       pdf(0*sign1,1,9)*polcolgq
c g s: 62 7
     &     + borncurr(0,0,2,62) * fincol(7) *
     &       pdf(0*sign1,1,9)*polcolgq
c g b: 62 18
     &     + borncurr(0,0,2,62) * fincol(18) *
     &       lbquark*pdf(0*sign1,1,9)*polcolgq
c g D: 42 12
     &     + borncurr(0,0,2,42) * fincol(12) *
     &       pdf(0*sign1,1,9)*polcolgq
c g S: 42 13
     &     + borncurr(0,0,2,42) * fincol(13) *
     &       pdf(0*sign1,1,9)*polcolgq
c g B: 42 20
     &     + borncurr(0,0,2,42) * fincol(20) *
     &       lbquark*pdf(0*sign1,1,9)*polcolgq
c g U: 42 16
     &     + borncurr(0,0,1,42) * fincol(16) *
     &       pdf(0*sign1,1,9)*polcolgq
c g C: 42 17
     &     + borncurr(0,0,1,42) * fincol(17) *
     &       pdf(0*sign1,1,9)*polcolgq
c U u: 41 4
     &     + borncurr(0,0,1,41) * fincol(4) *
     &       pdf(-2*sign1,1,9)*polcolqq
c C c: 41 5
     &     + borncurr(0,0,1,41) * fincol(5) *
     &       pdf(-4*sign1,1,9)*polcolqq
c D d: 41 8
     &     + borncurr(0,0,2,41) * fincol(8) *
     &       pdf(-1*sign1,1,9)*polcolqq
c S s: 41 9
     &     + borncurr(0,0,2,41) * fincol(9) *
     &       pdf(-3*sign1,1,9)*polcolqq
c B b: 41 19
     &     + borncurr(0,0,2,41) * fincol(19) *
     &       lbquark*pdf(-5*sign1,1,9)*polcolqq
c d D: 43 10
     &     + borncurr(0,0,2,43) * fincol(10) *
     &       pdf(1*sign1,1,9)*polcolqq
c s S: 43 11
     &     + borncurr(0,0,2,43) * fincol(11) *
     &       pdf(3*sign1,1,9)*polcolqq
c b B: 43 21
     &     + borncurr(0,0,2,43) * fincol(21) *
     &       lbquark*pdf(5*sign1,1,9)*polcolqq
c u U: 43 14
     &     + borncurr(0,0,1,43) * fincol(14) *
     &       pdf(2*sign1,1,9)*polcolqq
c c C: 43 15
     &     + borncurr(0,0,1,43) * fincol(15) *
     &       pdf(4*sign1,1,9)*polcolqq
     &   )
       endif 

c*****************  end of process evaluation  *********************

c sum into column 0 which steers integration
      do ires = 1,nmax
        res(ires,0) = 0d0
        do j=1,n_kin
          res(ires,0) = res(ires,0) + res(ires,j) 
        enddo
      enddo

      if (init.eq.1) then
         init = init+1
      endif
      do j=0,n_kin
         m2s(j) = 0
         DO IRES = 1,NMAX
            m2s(j) = m2s(j) + RES(IRES,j)
         ENDDO

         if(j.eq.0)then
c...Les Houches interface - the most propable subprocess at LO  
            if ((lha.or.hepmc).and..not.doNLO) then
               i=0
               weight=0.d0
               rnumb=RandomNumber()
               do while((i.le.nmax).and.(weight.le.rnumb*m2s(0)))
                  i=i+1
                  weight=weight+res(i,0)
                  iprocess=i
               enddo
            endif
         endif
        
         m2s(j) = m2s(j) * 4d0     ! factor 4 for random helicity summation

      enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c checks for IR finiteness

      if (ldebug) then
      if (abs(dotrr(p(0,2,1),p(0,4,1))).le.1e-1.or.abs(dotrr(p(0,1,1),p(0,4,1))).le.1e-1.or.
     &    abs(dotrr(p(0,2,1),p(0,3,1))).le.1e-1.or.abs(dotrr(p(0,1,1),p(0,3,1))).le.1e-1.or.
     &    abs(dotrr(p(0,3,1),p(0,4,1))).le.1e-1) then
!      if (abs(m2s(1)/(m2s(10)+m2s(11))).ge.0.9.and.abs(m2s(1)/(m2s(3)+m2s(5))).le.1.1) then
!      if (abs(dotrr(p(0,2,1),p(0,4,1))).le.1e-1) then
      print*, "p1p2=",dotrr(p(0,1,1),p(0,2,1))
      print*, "p1p3=",dotrr(p(0,1,1),p(0,3,1))
      print*, "p1p4=",dotrr(p(0,1,1),p(0,4,1))
      print*, "p2p3=",dotrr(p(0,2,1),p(0,3,1))
      print*, "p2p4=",dotrr(p(0,2,1),p(0,4,1))
      print*, "p3p4=",dotrr(p(0,3,1),p(0,4,1))
      print*, "drive integration via ",  m2s(0), " Fract ", m2s(0)/m2s(1)
      print*, "real emission ME      ",  m2s(1) , lok( 1)
      print*, "Dipole 1  (14)        ", -m2s(2) , lok( 2)
      print*, "Dipole 2  (24)        ", -m2s(3) , lok( 3)
      print*, "Dipole 3  (14)        ", -m2s(4) , lok( 4)
      print*, "Dipole 4  (24)        ", -m2s(5) , lok( 5)
      print*, "Dipole 5  (13)        ", -m2s(6) , lok( 6)
      print*, "Dipole 6  (23)        ", -m2s(7) , lok( 7)
      print*, "Dipole 7  (13)        ", -m2s(8) , lok( 8)
      print*, "Dipole 8  (23)        ", -m2s(9) , lok( 9)
      print*, "Dipole 9  (34)        ", -m2s(10), lok(10)
      print*, "Dipole 10 (34)        ", -m2s(11), lok(11)
      print*, "REM/Dipole 1+3  (14)  ", -m2s(1)/(m2s( 2)+m2s( 4))
      print*, "REM/Dipole 2+4  (24)  ", -m2s(1)/(m2s( 3)+m2s( 5))
      print*, "REM/Dipole 5+7  (13)  ", -m2s(1)/(m2s( 6)+m2s( 8))
      print*, "REM/Dipole 6+8  (23)  ", -m2s(1)/(m2s( 7)+m2s( 9))
      print*, "REM/Dipole 9+10 (34)  ", -m2s(1)/(m2s(10)+m2s(11))
      print*, "REM/Dipole tot        ", -m2s(1)/(m2s(2)+m2s(3)+m2s(4)
     &    +m2s(5)+m2s(6)+m2s(7)+m2s(8)+m2s(9)+m2s(10)+m2s(11))
      print*, "---------------------------------"
      read(*,*)
c      do i=1,70
c        print *, "Dipole diag ", i, " : "
c        print*, "real emission ME      ",  res(i,1) , lok( 1)
c        print*, "Dipole 1  (14)        ", -res(i,2) , lok( 2)
c        print*, "Dipole 2  (24)        ", -res(i,3) , lok( 3)
c        print*, "Dipole 3  (14)        ", -res(i,4) , lok( 4)
c        print*, "Dipole 4  (24)        ", -res(i,5) , lok( 5)
c        print*, "Dipole 5  (13)        ", -res(i,6) , lok( 6)
c        print*, "Dipole 6  (23)        ", -res(i,7) , lok( 7)
c        print*, "Dipole 7  (13)        ", -res(i,8) , lok( 8)
c        print*, "Dipole 8  (23)        ", -res(i,9) , lok( 9)
c        print*, "Dipole 9  (34)        ", -res(i,10), lok(10)
c        print*, "Dipole 10 (34)        ", -res(i,11), lok(11)
c        print*, "REM/Dipole 1+3  (14)  ", -res(i,1)/(res(i, 2)+res(i, 4))
c        print*, "REM/Dipole 2+4  (24)  ", -res(i,1)/(res(i, 3)+res(i, 5))
c        print*, "REM/Dipole 5+7  (13)  ", -res(i,1)/(res(i, 6)+res(i, 8))
c        print*, "REM/Dipole 6+8  (23)  ", -res(i,1)/(res(i, 7)+res(i, 9))
c        print*, "REM/Dipole 9+10 (34)  ", -res(i,1)/(res(i,10)+res(i,11))
c        print*, "REM/Dipole tot        ", -res(i,1)/(res(i,2)+res(i,3)+res(i,4)
c     &      +res(i,5)+res(i,6)+res(i,7)+res(i,8)+res(i,9)+res(i,10)+res(i,11))
c        read(*,*) 
c      enddo
      endif
      endif

      if (.false.) then
        open(unit=31,file="irfinite.chk",ACCESS='APPEND')
        kintrue=.false.
        do i=1,4
          do j=i,4
            if (i.eq.j) then ! soft
              kininv(i,i) = p(0,i,1)
            else !collinear
              kininv(i,j) = sqrt(abs(dotrr(p(0,i,1),p(0,j,1))))
            endif
            if (kininv(i,j) .lt. 1d0) kintrue=.true.
          enddo
        enddo
        ratio(1) = - m2s(1)/(m2s( 2)+m2s( 4))
        ratio(2) = - m2s(1)/(m2s( 3)+m2s( 5))
        ratio(3) = - m2s(1)/(m2s( 6)+m2s( 8))
        ratio(4) = - m2s(1)/(m2s( 7)+m2s( 9))
        ratio(5) = - m2s(1)/(m2s(10)+m2s(11))
        ratio(6) =  -m2s(1)/(m2s(2)+m2s(3)+m2s(4)
     &    +m2s(5)+m2s(6)+m2s(7)+m2s(8)+m2s(9)+m2s(10)+m2s(11))
c        do i=1,5
c          if ( (ratio(i).gt.0) .and. (ratio(i).lt.5) ) kintrue = .true.
c        enddo
        if (kintrue) then ! write it out
          write(31,*) (kininv(i,i),i=3,4), 
     &                ((kininv(i,j),j=i+1,4),i=1,4), 
     &                ratio
        endif
        close(31)
      endif

      end

c*****************************************************************************
c
c    end subroutine m2s_qqVVjj
c
c*****************************************************************************

