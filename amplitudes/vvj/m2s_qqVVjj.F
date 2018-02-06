c*********************  subroutine m2s_qqVVjj *********************************
c
c INPUT     bos       Boson identifier, 2 = Z, 6=H
c           nlo       NLO=1: create subtraction term; LO = 0
c           lok       need to calculate m2s(1:n_kin)? T or F 
c           xi        Feynman x parameters of incoming fermions
c           p         momenta associated with external fermions
c           v         momenta associated with WW decay fermions
c           rn        random number array
c           xuz       x_ija for subtraction terms
c          
c
c OUTPUT    m2s       |M|^2*pdf1*pdf2 
c
c       Adapted from m2s_wbfzhg by D. Zeppenfeld
c
c       Michael Rauch, <rauch@particle.physik.uni-karlsruhe.de>
c	Initial version: July 2011
c	Last modified: December 2012
c
c***************************************************************************

      subroutine m2s_qqVVjj(bos,nlo,lok,xi,p,v,rn,xuz,m2s)

      implicit none

c declare input/output variables

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/process.inc"

      integer bos
      real*8 xi(nx), p(0:3,max_p,max_kin), v(0:3,max_v,max_kin), 
     &       pch(0:3,max_p,max_kin),
     1       rn(1), xuz(2,2:11), m2s(0:n_kin)
      logical lok(max_kin)

      double precision N
      parameter(N=3d0)
      DOUBLE PRECISION zero
      PARAMETER (zero=0d0)

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

c declare variables for hadronic decay 
      integer N_gen_W
      real*8 fac_W, mjj2
      external mjj2

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
      double precision mgg(2,2,4), mqq(7,2,2,4)

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
      integer diag(71)
      double precision polcolqq, polcolgq, polcolgg
      double precision dipole(18)
      double precision fincol(17), fincolglobal, fincolfac
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
c parameters for call of qqVVjj.f

      if ( init .ne. 0) then
         do i = 1,nmax
            do j=0,n_kin
               res(i,j) = 0
            enddo
         enddo
      else if (init.eq.0) then
         write(6,*) " "
         write(6,*) "W+W- j j amplitude square information:"
         write(6,*) "-----------------------------------------------"
         write(6,*) " "
         call dibostoww_reset
         call printnfl(.true.)
         lbquark = 0
         init = 1
         do i = 1,maxnumsubproc
            do j = 0,n_kin
               res(i,j) = 0
            enddo
         enddo
         zm2i(2) = 1d0/dcmplx(xm2(2),-xmg(2))
         zm2i(3) = 1d0/dcmplx(xm2(3),-xmg(3))

         nmax = 71
         do i=1,nmax
           diag(i) = 1
         enddo
cmr gg
c         diag(35) = 0
c         diag(36) = 0
cmr gq et al
c         diag( 6) = 0
c         diag(22) = 0 
c         diag(33) = 0 
c         diag(34) = 0 
c         diag(37) = 0 
c         diag(38) = 0 
c         diag(49) = 0 
c         diag(65) = 0 
cmr qq->gg 
c         diag(16) = 0
c         diag(30) = 0
c         diag(48) = 0
c         diag(62) = 0
cmr qq->qq
c         diag(13) = 1
c         diag(14) = 1
c dipoledebug
c         diag(1) = 1

      endif

c ------------------------------------------------------------------- 

c precalculate leptonic tensors
c lepton spinors and W+- polarization vectors

c...Les Houches interface
      if ((lha.or.hepmc).and..not.doNLO) then
         helicity(1)=-1
         helicity(2)= 1
         helicity(3)=-1
         helicity(4)= 1
      endif

c leptonic tensors
c all kinematics which have unshifted lepton momenta (all but initial-initial)
      if (lok(1).or.lok(2).or.lok(3).or.lok(6).or.lok(7).or.lok(10).or.lok(11)) then
         J=1
         CALL IXXXXX(v(0,2,J),ZERO ,1,-1,wep(1,J)) !W(1,5))          !e+       
         CALL OXXXXX(v(0,1,J),ZERO ,-1,1,wve(1,J)) !W(1,4))          !ve 
         CALL OXXXXX(v(0,3,J),ZERO ,-1,1,wmu(1,J)) !W(1,6))          !mu-      
         CALL IXXXXX(v(0,4,J),ZERO ,1,-1,wvm(1,J)) !W(1,7))          !vm~
         CALL JIOXXX(wep(1,J),wve(1,J),GWF ,WMASS,WWIDTH,wp(1,J)) !W(1,8))!W+
         CALL JIOXXX(wvm(1,J),wmu(1,J),GWF ,WMASS,WWIDTH,wm(1,J)) !W(1,9))!W-

         do mu = 0,3
            qp(mu,J) = v(mu,1,J)+v(mu,2,J)
            qm(mu,J) = v(mu,3,J)+v(mu,4,J)
            qww(mu,J) = qp(mu,J) + qm(mu,J) 
         enddo
         qp(4,J) = dotrr(qp(0,J),qp(0,J))
         qm(4,J) = dotrr(qm(0,J),qm(0,J))
         qww(4,J) = dotrr(qww(0,J),qww(0,J))         

         SELECT CASE(process)     
         CASE(WPHADWMJ, WPWMHADJ, WpHadWmJJLO, WpwMHadJJLO)
c hadronic tensors
           call DIBOSZtoQQll(v(0,1,J), J, zww(0,J))
           call DIBOSAtoQQll(v(0,1,J), J, aww(0,J))

         CASE DEFAULT

c leptonic tensors

           call DIBOSatoww(v(0,1,J),J,aww(0,J))
           call DIBOSztoww(v(0,1,J),J,zww(0,J))
         END SELECT
c define epsilon~ = espilon - k * k.epsilon/M**2 
         help1c = -zm2i(2) * dotrc(qww(0,J),zww(0,J))
         do mu=0,3
            zww(mu,J) = zww(mu,J) + qww(mu,J) * help1c
         enddo

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
             wm(mu,j) = wm(mu,1)
             aww(mu-1,j) = aww(mu-1,1)
             zww(mu-1,j) = zww(mu-1,1)
           enddo
           do mu = 0,4
              qp(mu,J)  = qp(mu,1)
              qm(mu,J)  = qm(mu,1)
              qww(mu,J) = qww(mu,1)
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
           CALL IXXXXX(v(0,2,J),ZERO ,1,-1,wep(1,J)) !W(1,5))          !e+       
           CALL OXXXXX(v(0,1,J),ZERO ,-1,1,wve(1,J)) !W(1,4))          !ve 
           CALL OXXXXX(v(0,3,J),ZERO ,-1,1,wmu(1,J)) !W(1,6))          !mu-      
           CALL IXXXXX(v(0,4,J),ZERO ,1,-1,wvm(1,J)) !W(1,7))          !vm~
           CALL JIOXXX(wep(1,J),wve(1,J),GWF ,WMASS,WWIDTH,wp(1,J)) !W(1,8))!W+
           CALL JIOXXX(wvm(1,J),wmu(1,J),GWF ,WMASS,WWIDTH,wm(1,J)) !W(1,9))!W-

           do mu = 0,3
              qp(mu,J) = v(mu,1,J)+v(mu,2,J)
              qm(mu,J) = v(mu,3,J)+v(mu,4,J)
              qww(mu,J) = qp(mu,J) + qm(mu,J) 
           enddo
           qp(4,J) = dotrr(qp(0,J),qp(0,J))
           qm(4,J) = dotrr(qm(0,J),qm(0,J))
           qww(4,J) = dotrr(qww(0,J),qww(0,J))    

           SELECT CASE(process)     
           CASE(WPHADWMJ, WPWMHADJ, WpHadWmJJLO, WpwMHadJJLO)
c hadronic tensors
             call DIBOSZtoQQll(v(0,1,J), J, zww(0,J))
             call DIBOSAtoQQll(v(0,1,J), J, aww(0,J))

           CASE DEFAULT     
c leptonic tensors
             call DIBOSatoww(v(0,1,J),J,aww(0,J))
             call DIBOSztoww(v(0,1,J),J,zww(0,J))
           END SELECT

c define epsilon~ = espilon - k * k.epsilon/M**2 
           help1c = -zm2i(2) * dotrc(qww(0,J),zww(0,J))
           do mu=0,3
              zww(mu,J) = zww(mu,J) + qww(mu,J) * help1c
           enddo
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

C     Now sum over the subprocesses contributing to W+ W- production

C********************************************************************

      if (lok(1)) then

         do i=1,4
           call qqVVqq(p,v,i,mqq(1,1,1,i))
           call qqVVgg(p,v,i,mgg(1,1,i))  
           if (i.ne.4) then
             call qqVVgg(pch,v,i,mgg(1,2,i))  
           endif
         enddo
cmr      do mu=0,3
cmr        pbar(mu,1)=p(mu,1,1)
cmr        pbar(mu,2)=p(mu,2,1)
cmr        pbar(mu,3)=p(mu,3,1)
cmr        pbar(mu,4)=p(mu,4,1)
cmr        pbar(mu,5)=v(mu,2,1)
cmr        pbar(mu,6)=v(mu,1,1)
cmr        pbar(mu,7)=v(mu,3,1)
cmr        pbar(mu,8)=v(mu,4,1)
cmr      enddo
cmr      call SUUWWUU(pbar,mg2res)
cmr      write(*,*) "MMM", mqq(3,1,1,1)/72d0, mg2res/g**4*(4*pi*als(1,1))**2, mg2res/g**4*(4*pi*als(1,1))**2/mqq(3,1,1,1)*72d0

c (uu->uu) + (cc->cc)
         res(1,1) = diag(1) * mqq(3,1,1,1) * 
     &      ( pdf(2*sign1,1,1)*pdf(2*sign2,2,1) + pdf(4*sign1,1,1)*pdf(4*sign2,2,1))/2d0 * polcolqq
c (ud->du) + (cs->sc)
         res(2,1) = diag(2) * mqq(7,2,1,1) * 
     &      ( pdf(2*sign1,1,1)*pdf(1*sign2,2,1) + pdf(4*sign1,1,1)*pdf(3*sign2,2,1)) * polcolqq
c (uc->uc) + (cu->cu)
         res(3,1) = diag(3) * mqq(1,1,1,1) * 
     &      ( pdf(2*sign1,1,1)*pdf(4*sign2,2,1) + pdf(4*sign1,1,1)*pdf(2*sign2,2,1)) * polcolqq
c (us->us) + (cd->cd)
         res(4,1) = diag(4) * mqq(1,1,2,1) * 
     &      ( pdf(2*sign1,1,1)*pdf(3*sign2,2,1) + pdf(4*sign1,1,1)*pdf(1*sign2,2,1)) * polcolqq
c (us->dc) + (cd->su)
         res(5,1) = diag(5) * mqq(4,2,1,1) * 
     &      ( pdf(2*sign1,1,1)*pdf(3*sign2,2,1) + pdf(4*sign1,1,1)*pdf(1*sign2,2,1)) * polcolqq
c (uS->uS) + (cD->cD)
         res(7,1) = diag(7) * mqq(1,1,2,2) * 
     &      ( pdf(2*sign1,1,1)*pdf(-3*sign2,2,1) + pdf(4*sign1,1,1)*pdf(-1*sign2,2,1)) * polcolqq
c (uC->uC) + (cU->cU)
         res(8,1) = diag(8) * mqq(1,1,1,2) * 
     &      ( pdf(2*sign1,1,1)*pdf(-4*sign2,2,1) + pdf(4*sign1,1,1)*pdf(-2*sign2,2,1)) * polcolqq
c (uC->dS) + (cU->sD)
         res(9,1) = diag(9) * mqq(4,2,1,2) * 
     &      ( pdf(2*sign1,1,1)*pdf(-4*sign2,2,1) + pdf(4*sign1,1,1)*pdf(-2*sign2,2,1)) * polcolqq
c (uD->uD) + (cS->cS)
         res(10,1) = diag(10) * mqq(6,1,2,2) * 
     &      ( pdf(2*sign1,1,1)*pdf(-1*sign2,2,1) + pdf(4*sign1,1,1)*pdf(-3*sign2,2,1)) * polcolqq
c (uD->cS) + (cS->uD)
         res(11,1) = diag(11) * mqq(5,1,2,2) * 
     &      ( pdf(2*sign1,1,1)*pdf(-1*sign2,2,1) + pdf(4*sign1,1,1)*pdf(-3*sign2,2,1)) * polcolqq
c (uU->uU) + (cC->cC)
         res(12,1) = diag(12) * mqq(3,1,1,2) * 
     &      ( pdf(2*sign1,1,1)*pdf(-2*sign2,2,1) + pdf(4*sign1,1,1)*pdf(-4*sign2,2,1)) * polcolqq
c (uU->dD) + (cC->sS)
         res(13,1) = diag(13) * mqq(7,2,1,2) * 
     &      ( pdf(2*sign1,1,1)*pdf(-2*sign2,2,1) + pdf(4*sign1,1,1)*pdf(-4*sign2,2,1)) * polcolqq
c (uU->cC) + (cC->uU)
         res(14,1) = diag(14) * mqq(2,1,1,2) * 
     &      ( pdf(2*sign1,1,1)*pdf(-2*sign2,2,1) + pdf(4*sign1,1,1)*pdf(-4*sign2,2,1)) * polcolqq
c (uU->sS) + (cC->dD)
         res(15,1) = diag(15) * mqq(2,2,1,2) * 
     &      ( pdf(2*sign1,1,1)*pdf(-2*sign2,2,1) + pdf(4*sign1,1,1)*pdf(-4*sign2,2,1)) * polcolqq
c (du->du) + (sc->sc)
         res(17,1) = diag(17) * mqq(6,2,1,1) * 
     &      ( pdf(1*sign1,1,1)*pdf(2*sign2,2,1) + pdf(3*sign1,1,1)*pdf(4*sign2,2,1)) * polcolqq
c (dd->dd) + (ss->ss)
         res(18,1) = diag(18) * mqq(3,2,2,1) * 
     &      ( pdf(1*sign1,1,1)*pdf(1*sign2,2,1) + pdf(3*sign1,1,1)*pdf(3*sign2,2,1))/2d0 * polcolqq
c (dc->us) + (su->cd)
         res(19,1) = diag(19) * mqq(4,1,2,1) * 
     &      ( pdf(1*sign1,1,1)*pdf(4*sign2,2,1) + pdf(3*sign1,1,1)*pdf(2*sign2,2,1)) * polcolqq
c (dc->dc) + (su->su)
         res(20,1) = diag(20) * mqq(1,2,1,1) * 
     &      ( pdf(1*sign1,1,1)*pdf(4*sign2,2,1) + pdf(3*sign1,1,1)*pdf(2*sign2,2,1)) * polcolqq
c (ds->ds) + (sd->sd)
         res(21,1) = diag(21) * mqq(1,2,2,1) * 
     &      ( pdf(1*sign1,1,1)*pdf(3*sign2,2,1) + pdf(3*sign1,1,1)*pdf(1*sign2,2,1)) * polcolqq
c (dS->uC) + (sD->cU)
         res(23,1) = diag(23) * mqq(4,1,2,2) * 
     &      ( pdf(1*sign1,1,1)*pdf(-3*sign2,2,1) + pdf(3*sign1,1,1)*pdf(-1*sign2,2,1)) * polcolqq
c (dS->dS) + (sD->sD)
         res(24,1) = diag(24) * mqq(1,2,2,2) * 
     &      ( pdf(1*sign1,1,1)*pdf(-3*sign2,2,1) + pdf(3*sign1,1,1)*pdf(-1*sign2,2,1)) * polcolqq
c (dC->dC) + (sU->sU)
         res(25,1) = diag(25) * mqq(1,2,1,2) * 
     &      ( pdf(1*sign1,1,1)*pdf(-4*sign2,2,1) + pdf(3*sign1,1,1)*pdf(-2*sign2,2,1)) * polcolqq
c (dD->uU) + (sS->cC)
         res(26,1) = diag(26) * mqq(7,1,2,2) * 
     &      ( pdf(1*sign1,1,1)*pdf(-1*sign2,2,1) + pdf(3*sign1,1,1)*pdf(-3*sign2,2,1)) * polcolqq
c (dD->dD) + (sS->sS)
         res(27,1) = diag(27) * mqq(3,2,2,2) * 
     &      ( pdf(1*sign1,1,1)*pdf(-1*sign2,2,1) + pdf(3*sign1,1,1)*pdf(-3*sign2,2,1)) * polcolqq
c (dD->cC) + (sS->uU)
         res(28,1) = diag(28) * mqq(2,1,2,2) * 
     &      ( pdf(1*sign1,1,1)*pdf(-1*sign2,2,1) + pdf(3*sign1,1,1)*pdf(-3*sign2,2,1)) * polcolqq
c (dD->sS) + (sS->dD)
         res(29,1) = diag(29) * mqq(2,2,2,2) * 
     &      ( pdf(1*sign1,1,1)*pdf(-1*sign2,2,1) + pdf(3*sign1,1,1)*pdf(-3*sign2,2,1)) * polcolqq
c (dU->dU) + (sC->sC)
         res(31,1) = diag(31) * mqq(6,2,1,2) * 
     &      ( pdf(1*sign1,1,1)*pdf(-2*sign2,2,1) + pdf(3*sign1,1,1)*pdf(-4*sign2,2,1)) * polcolqq
c (dU->sC) + (sC->dU)
         res(32,1) = diag(32) * mqq(5,2,1,2) * 
     &      ( pdf(1*sign1,1,1)*pdf(-2*sign2,2,1) + pdf(3*sign1,1,1)*pdf(-4*sign2,2,1)) * polcolqq
c (Dc->cD) + (Su->uS)
         res(39,1) = diag(39) * mqq(2,2,1,3) * 
     &      ( pdf(-1*sign1,1,1)*pdf(4*sign2,2,1) + pdf(-3*sign1,1,1)*pdf(2*sign2,2,1)) * polcolqq
c (Ds->cU) + (Sd->uC)
         res(40,1) = diag(40) * mqq(5,2,1,3) * 
     &      ( pdf(-1*sign1,1,1)*pdf(3*sign2,2,1) + pdf(-3*sign1,1,1)*pdf(1*sign2,2,1)) * polcolqq
c (Ds->sD) + (Sd->dS)
         res(41,1) = diag(41) * mqq(2,2,2,3) * 
     &      ( pdf(-1*sign1,1,1)*pdf(3*sign2,2,1) + pdf(-3*sign1,1,1)*pdf(1*sign2,2,1)) * polcolqq
c (Du->cS) + (Sc->uD)
         res(42,1) = diag(42) * mqq(4,2,1,3) * 
     &      ( pdf(-1*sign1,1,1)*pdf(2*sign2,2,1) + pdf(-3*sign1,1,1)*pdf(4*sign2,2,1)) * polcolqq
c (Du->uD) + (Sc->cS)
         res(43,1) = diag(43) * mqq(7,2,1,3) * 
     &      ( pdf(-1*sign1,1,1)*pdf(2*sign2,2,1) + pdf(-3*sign1,1,1)*pdf(4*sign2,2,1)) * polcolqq
c (Dd->cC) + (Ss->uU)
         res(44,1) = diag(44) * mqq(1,2,1,3) * 
     &      ( pdf(-1*sign1,1,1)*pdf(1*sign2,2,1) + pdf(-3*sign1,1,1)*pdf(3*sign2,2,1)) * polcolqq
c (Dd->sS) + (Ss->dD)
         res(45,1) = diag(45) * mqq(1,2,2,3) * 
     &      ( pdf(-1*sign1,1,1)*pdf(1*sign2,2,1) + pdf(-3*sign1,1,1)*pdf(3*sign2,2,1)) * polcolqq
c (Dd->uU) + (Ss->cC)
         res(46,1) = diag(46) * mqq(6,2,1,3) * 
     &      ( pdf(-1*sign1,1,1)*pdf(1*sign2,2,1) + pdf(-3*sign1,1,1)*pdf(3*sign2,2,1)) * polcolqq
c (Dd->dD) + (Ss->sS)
         res(47,1) = diag(47) * mqq(3,2,2,3) * 
     &      ( pdf(-1*sign1,1,1)*pdf(1*sign2,2,1) + pdf(-3*sign1,1,1)*pdf(3*sign2,2,1)) * polcolqq
c (DD->DD) + (SS->SS)
         res(50,1) = diag(50) * mqq(3,2,2,4) * 
     &      ( pdf(-1*sign1,1,1)*pdf(-1*sign2,2,1) + pdf(-3*sign1,1,1)*pdf(-3*sign2,2,1))/2d0 * polcolqq
c (DU->UD) + (SC->CS)
         res(51,1) = diag(51) * mqq(6,1,2,4) * 
     &      ( pdf(-1*sign1,1,1)*pdf(-2*sign2,2,1) + pdf(-3*sign1,1,1)*pdf(-4*sign2,2,1)) * polcolqq
c (DS->SD) + (SD->DS)
         res(52,1) = diag(52) * mqq(1,2,2,4) * 
     &      ( pdf(-1*sign1,1,1)*pdf(-3*sign2,2,1) + pdf(-3*sign1,1,1)*pdf(-1*sign2,2,1)) * polcolqq
c (DC->CD) + (SU->US)
         res(53,1) = diag(53) * mqq(1,1,2,4) * 
     &      ( pdf(-1*sign1,1,1)*pdf(-4*sign2,2,1) + pdf(-3*sign1,1,1)*pdf(-2*sign2,2,1)) * polcolqq
c (DC->SU) + (SU->DC)
         res(54,1) = diag(54) * mqq(4,1,2,4) * 
     &      ( pdf(-1*sign1,1,1)*pdf(-4*sign2,2,1) + pdf(-3*sign1,1,1)*pdf(-2*sign2,2,1)) * polcolqq
c (Uc->cU) + (Cu->uC)
         res(55,1) = diag(55) * mqq(2,1,1,3) * 
     &      ( pdf(-2*sign1,1,1)*pdf(4*sign2,2,1) + pdf(-4*sign1,1,1)*pdf(2*sign2,2,1)) * polcolqq
c (Uc->sD) + (Cu->dS)
         res(56,1) = diag(56) * mqq(5,1,2,3) * 
     &      ( pdf(-2*sign1,1,1)*pdf(4*sign2,2,1) + pdf(-4*sign1,1,1)*pdf(2*sign2,2,1)) * polcolqq
c (Us->sU) + (Cd->dC)
         res(57,1) = diag(57) * mqq(2,1,2,3) * 
     &      ( pdf(-2*sign1,1,1)*pdf(3*sign2,2,1) + pdf(-4*sign1,1,1)*pdf(1*sign2,2,1)) * polcolqq
c (Uu->cC) + (Cc->uU)
         res(58,1) = diag(58) * mqq(1,1,1,3) * 
     &      ( pdf(-2*sign1,1,1)*pdf(2*sign2,2,1) + pdf(-4*sign1,1,1)*pdf(4*sign2,2,1)) * polcolqq
c (Uu->sS) + (Cc->dD)
         res(59,1) = diag(59) * mqq(1,1,2,3) * 
     &      ( pdf(-2*sign1,1,1)*pdf(2*sign2,2,1) + pdf(-4*sign1,1,1)*pdf(4*sign2,2,1)) * polcolqq
c (Uu->uU) + (Cc->cC)
         res(60,1) = diag(60) * mqq(3,1,1,3) * 
     &      ( pdf(-2*sign1,1,1)*pdf(2*sign2,2,1) + pdf(-4*sign1,1,1)*pdf(4*sign2,2,1)) * polcolqq
c (Uu->dD) + (Cc->sS)
         res(61,1) = diag(61) * mqq(6,1,2,3) * 
     &      ( pdf(-2*sign1,1,1)*pdf(2*sign2,2,1) + pdf(-4*sign1,1,1)*pdf(4*sign2,2,1)) * polcolqq
c (Ud->sC) + (Cs->dU)
         res(63,1) = diag(63) * mqq(4,1,2,3) * 
     &      ( pdf(-2*sign1,1,1)*pdf(1*sign2,2,1) + pdf(-4*sign1,1,1)*pdf(3*sign2,2,1)) * polcolqq
c (Ud->dU) + (Cs->sC)
         res(64,1) = diag(64) * mqq(7,1,2,3) * 
     &      ( pdf(-2*sign1,1,1)*pdf(1*sign2,2,1) + pdf(-4*sign1,1,1)*pdf(3*sign2,2,1)) * polcolqq
c (UD->UD) + (CS->CS)
         res(66,1) = diag(66) * mqq(7,2,1,4) * 
     &      ( pdf(-2*sign1,1,1)*pdf(-1*sign2,2,1) + pdf(-4*sign1,1,1)*pdf(-3*sign2,2,1)) * polcolqq
c (UU->UU) + (CC->CC)
         res(67,1) = diag(67) * mqq(3,1,1,4) * 
     &      ( pdf(-2*sign1,1,1)*pdf(-2*sign2,2,1) + pdf(-4*sign1,1,1)*pdf(-4*sign2,2,1))/2d0 * polcolqq
c (US->CD) + (CD->US)
         res(68,1) = diag(68) * mqq(4,2,1,4) * 
     &      ( pdf(-2*sign1,1,1)*pdf(-3*sign2,2,1) + pdf(-4*sign1,1,1)*pdf(-1*sign2,2,1)) * polcolqq
c (US->SU) + (CD->DC)
         res(69,1) = diag(69) * mqq(1,2,1,4) * 
     &      ( pdf(-2*sign1,1,1)*pdf(-3*sign2,2,1) + pdf(-4*sign1,1,1)*pdf(-1*sign2,2,1)) * polcolqq
c (UC->CU) + (CU->UC)
         res(70,1) = diag(70) * mqq(1,1,1,4) * 
     &      ( pdf(-2*sign1,1,1)*pdf(-4*sign2,2,1) + pdf(-4*sign1,1,1)*pdf(-2*sign2,2,1)) * polcolqq
c (uU->gg) + (cC->gg)
         res(16,1) = diag(16) * mgg(1,1,1) * 
     &      ( pdf(2*sign1,1,1)*pdf(-2*sign2,2,1) + pdf(4*sign1,1,1)*pdf(-4*sign2,2,1))/2d0 * polcolqq
c (dD->gg) + (sS->gg)
         res(30,1) = diag(30) * mgg(2,1,1) * 
     &      ( pdf(1*sign1,1,1)*pdf(-1*sign2,2,1) + pdf(3*sign1,1,1)*pdf(-3*sign2,2,1))/2d0 * polcolqq
c (Dd->gg) + (Ss->gg)
         res(48,1) = diag(48) * mgg(2,2,1) * 
     &      ( pdf(-1*sign1,1,1)*pdf(1*sign2,2,1) + pdf(-3*sign1,1,1)*pdf(3*sign2,2,1))/2d0 * polcolqq
c (Uu->gg) + (Cc->gg)
         res(62,1) = diag(62) * mgg(1,2,1) * 
     &      ( pdf(-2*sign1,1,1)*pdf(2*sign2,2,1) + pdf(-4*sign1,1,1)*pdf(4*sign2,2,1))/2d0 * polcolqq

c (ug->ug) + (cg->cg)
         res(6,1) = diag(6) * mgg(1,1,2) * 
     &      ( pdf(2*sign1,1,1)*pdf(0*sign2,2,1) + pdf(4*sign1,1,1)*pdf(0*sign2,2,1)) * polcolgq
c (dg->dg) + (sg->sg)
         res(22,1) = diag(22) * mgg(2,1,2) * 
     &      ( pdf(1*sign1,1,1)*pdf(0*sign2,2,1) + pdf(3*sign1,1,1)*pdf(0*sign2,2,1)) * polcolgq
c (gu->ug) + (gc->cg)
         res(33,1) = diag(33) * mgg(1,2,2) * 
     &      ( pdf(0*sign1,1,1)*pdf(2*sign2,2,1) + pdf(0*sign1,1,1)*pdf(4*sign2,2,1)) * polcolgq
c (gd->dg) + (gs->sg)
         res(34,1) = diag(34) * mgg(2,2,2) * 
     &      ( pdf(0*sign1,1,1)*pdf(1*sign2,2,1) + pdf(0*sign1,1,1)*pdf(3*sign2,2,1)) * polcolgq
c (gD->Dg) + (gS->Sg)
         res(37,1) = diag(37) * mgg(2,2,3) * 
     &      ( pdf(0*sign1,1,1)*pdf(-1*sign2,2,1) + pdf(0*sign1,1,1)*pdf(-3*sign2,2,1)) * polcolgq
c (gU->Ug) + (gC->Cg)
         res(38,1) = diag(38) * mgg(1,2,3) * 
     &      ( pdf(0*sign1,1,1)*pdf(-2*sign2,2,1) + pdf(0*sign1,1,1)*pdf(-4*sign2,2,1)) * polcolgq
c (Dg->Dg) + (Sg->Sg)
         res(49,1) = diag(49) * mgg(2,1,3) * 
     &      ( pdf(-1*sign1,1,1)*pdf(0*sign2,2,1) + pdf(-3*sign1,1,1)*pdf(0*sign2,2,1)) * polcolgq
c (Ug->Ug) + (Cg->Cg)
         res(65,1) = diag(65) * mgg(1,1,3) * 
     &      ( pdf(-2*sign1,1,1)*pdf(0*sign2,2,1) + pdf(-4*sign1,1,1)*pdf(0*sign2,2,1)) * polcolgq

c (gg->uU) + (gg->cC)
         res(35,1) = diag(35) * mgg(1,1,4) * 
     &      ( pdf(0*sign1,1,1)*pdf(0*sign2,2,1) + pdf(0*sign1,1,1)*pdf(0*sign2,2,1)) * polcolgg
c (gg->dD) + (gg->sS)
         res(36,1) = diag(36) * mgg(2,1,4) * 
     &      ( pdf(0*sign1,1,1)*pdf(0*sign2,2,1) + pdf(0*sign1,1,1)*pdf(0*sign2,2,1)) * polcolgg

       endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c dipoles
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       if (lok(2)) then 
         call qqgWW_curr(borncurr(0,-1,1,1),2,p(0,3,2),p(0,1,2),p(0,2,2),-1,-1,-1)
         call qqgWW_curr(borncurr(0,-1,1,2),2,p(0,2,2),p(0,1,2),p(0,3,2),1,-1,1)
         call qqgWW_curr(borncurr(0,-1,1,3),2,p(0,3,2),p(0,2,2),p(0,1,2),-1,-1,-1)
         call qqgWW_curr(borncurr(0,-1,1,4),2,p(0,2,2),p(0,3,2),p(0,1,2),1,1,-1)
         call qqgWW_curr(borncurr(0,-1,1,5),2,p(0,1,2),p(0,2,2),p(0,3,2),1,-1,1)
         call qqgWW_curr(borncurr(0,-1,1,6),2,p(0,1,2),p(0,3,2),p(0,2,2),1,1,-1)

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
     &      ( pdf(2*sign1,1,2)*pdf(2*sign2,2,2) + pdf(4*sign1,1,2)*pdf(4*sign2,2,2))/2d0*polcolgq
         res(2,2) = diag(2) * dipole(2) * 
     &      ( pdf(2*sign1,1,2)*pdf(1*sign2,2,2) + pdf(4*sign1,1,2)*pdf(3*sign2,2,2))*polcolgq
         res(6,2) = diag(6) * dipole(3) * 
     &      ( pdf(2*sign1,1,2)*pdf(0*sign2,2,2) + pdf(4*sign1,1,2)*pdf(0*sign2,2,2))*polcolgq
         res(16,2) = diag(16) * dipole(4) * 
     &      ( pdf(2*sign1,1,2)*pdf(-2*sign2,2,2) + pdf(4*sign1,1,2)*pdf(-4*sign2,2,2))/2d0*polcolqq
         res(18,2) = diag(18) * dipole(2) * 
     &      ( pdf(1*sign1,1,2)*pdf(1*sign2,2,2) + pdf(3*sign1,1,2)*pdf(3*sign2,2,2))/2d0*polcolgq
         res(22,2) = diag(22) * dipole(5) * 
     &      ( pdf(1*sign1,1,2)*pdf(0*sign2,2,2) + pdf(3*sign1,1,2)*pdf(0*sign2,2,2))*polcolgq
         res(30,2) = diag(30) * dipole(6) * 
     &      ( pdf(1*sign1,1,2)*pdf(-1*sign2,2,2) + pdf(3*sign1,1,2)*pdf(-3*sign2,2,2))/2d0*polcolqq
         res(33,2) = diag(33) * dipole(7) * 
     &      ( pdf(0*sign1,1,2)*pdf(2*sign2,2,2) + pdf(0*sign1,1,2)*pdf(4*sign2,2,2))*polcolgq
         res(34,2) = diag(34) * dipole(8) * 
     &      ( pdf(0*sign1,1,2)*pdf(1*sign2,2,2) + pdf(0*sign1,1,2)*pdf(3*sign2,2,2))*polcolgq
         res(35,2) = diag(35) * dipole(9) * 
     &      ( pdf(0*sign1,1,2)*pdf(0*sign2,2,2) + pdf(0*sign1,1,2)*pdf(0*sign2,2,2))*polcolgq
         res(36,2) = diag(36) * dipole(10) * 
     &      ( pdf(0*sign1,1,2)*pdf(0*sign2,2,2) + pdf(0*sign1,1,2)*pdf(0*sign2,2,2))*polcolgq
         res(37,2) = diag(37) * dipole(11) * 
     &      ( pdf(0*sign1,1,2)*pdf(-1*sign2,2,2) + pdf(0*sign1,1,2)*pdf(-3*sign2,2,2))*polcolgq
         res(38,2) = diag(38) * dipole(12) * 
     &      ( pdf(0*sign1,1,2)*pdf(-2*sign2,2,2) + pdf(0*sign1,1,2)*pdf(-4*sign2,2,2))*polcolgq
         res(39,2) = diag(39) * dipole(1) * 
     &      ( pdf(-1*sign1,1,2)*pdf(4*sign2,2,2) + pdf(-3*sign1,1,2)*pdf(2*sign2,2,2))*polcolgq
         res(41,2) = diag(41) * dipole(2) * 
     &      ( pdf(-1*sign1,1,2)*pdf(3*sign2,2,2) + pdf(-3*sign1,1,2)*pdf(1*sign2,2,2))*polcolgq
         res(43,2) = diag(43) * dipole(1) * 
     &      ( pdf(-1*sign1,1,2)*pdf(2*sign2,2,2) + pdf(-3*sign1,1,2)*pdf(4*sign2,2,2))*polcolgq
         res(47,2) = diag(47) * dipole(2) * 
     &      ( pdf(-1*sign1,1,2)*pdf(1*sign2,2,2) + pdf(-3*sign1,1,2)*pdf(3*sign2,2,2))*polcolgq
         res(48,2) = diag(48) * dipole(13) * 
     &      ( pdf(-1*sign1,1,2)*pdf(1*sign2,2,2) + pdf(-3*sign1,1,2)*pdf(3*sign2,2,2))/2d0*polcolqq
         res(49,2) = diag(49) * dipole(14) * 
     &      ( pdf(-1*sign1,1,2)*pdf(0*sign2,2,2) + pdf(-3*sign1,1,2)*pdf(0*sign2,2,2))*polcolgq
         res(50,2) = diag(50) * dipole(15) * 
     &      ( pdf(-1*sign1,1,2)*pdf(-1*sign2,2,2) + pdf(-3*sign1,1,2)*pdf(-3*sign2,2,2))/2d0*polcolgq
         res(51,2) = diag(51) * dipole(16) * 
     &      ( pdf(-1*sign1,1,2)*pdf(-2*sign2,2,2) + pdf(-3*sign1,1,2)*pdf(-4*sign2,2,2))*polcolgq
         res(52,2) = diag(52) * dipole(15) * 
     &      ( pdf(-1*sign1,1,2)*pdf(-3*sign2,2,2) + pdf(-3*sign1,1,2)*pdf(-1*sign2,2,2))*polcolgq
         res(53,2) = diag(53) * dipole(16) * 
     &      ( pdf(-1*sign1,1,2)*pdf(-4*sign2,2,2) + pdf(-3*sign1,1,2)*pdf(-2*sign2,2,2))*polcolgq
         res(55,2) = diag(55) * dipole(1) * 
     &      ( pdf(-2*sign1,1,2)*pdf(4*sign2,2,2) + pdf(-4*sign1,1,2)*pdf(2*sign2,2,2))*polcolgq
         res(57,2) = diag(57) * dipole(2) * 
     &      ( pdf(-2*sign1,1,2)*pdf(3*sign2,2,2) + pdf(-4*sign1,1,2)*pdf(1*sign2,2,2))*polcolgq
         res(60,2) = diag(60) * dipole(1) * 
     &      ( pdf(-2*sign1,1,2)*pdf(2*sign2,2,2) + pdf(-4*sign1,1,2)*pdf(4*sign2,2,2))*polcolgq
         res(62,2) = diag(62) * dipole(17) * 
     &      ( pdf(-2*sign1,1,2)*pdf(2*sign2,2,2) + pdf(-4*sign1,1,2)*pdf(4*sign2,2,2))/2d0*polcolqq
         res(64,2) = diag(64) * dipole(2) * 
     &      ( pdf(-2*sign1,1,2)*pdf(1*sign2,2,2) + pdf(-4*sign1,1,2)*pdf(3*sign2,2,2))*polcolgq
         res(65,2) = diag(65) * dipole(18) * 
     &      ( pdf(-2*sign1,1,2)*pdf(0*sign2,2,2) + pdf(-4*sign1,1,2)*pdf(0*sign2,2,2))*polcolgq
         res(67,2) = diag(67) * dipole(16) * 
     &      ( pdf(-2*sign1,1,2)*pdf(-2*sign2,2,2) + pdf(-4*sign1,1,2)*pdf(-4*sign2,2,2))/2d0*polcolgq
         res(69,2) = diag(69) * dipole(15) * 
     &      ( pdf(-2*sign1,1,2)*pdf(-3*sign2,2,2) + pdf(-4*sign1,1,2)*pdf(-1*sign2,2,2))*polcolgq
         res(70,2) = diag(70) * dipole(16) * 
     &      ( pdf(-2*sign1,1,2)*pdf(-4*sign2,2,2) + pdf(-4*sign1,1,2)*pdf(-2*sign2,2,2))*polcolgq
       endif 

       if (lok(3)) then 
         call qqgWW_curr(borncurr(0,-1,1,7),3,p(0,3,3),p(0,1,3),p(0,2,3),-1,-1,-1)
         call qqgWW_curr(borncurr(0,-1,1,8),3,p(0,2,3),p(0,1,3),p(0,3,3),1,-1,1)
         call qqgWW_curr(borncurr(0,-1,1,9),3,p(0,3,3),p(0,2,3),p(0,1,3),-1,-1,-1)
         call qqgWW_curr(borncurr(0,-1,1,10),3,p(0,2,3),p(0,3,3),p(0,1,3),1,1,-1)
         call qqgWW_curr(borncurr(0,-1,1,11),3,p(0,1,3),p(0,2,3),p(0,3,3),1,-1,1)
         call qqgWW_curr(borncurr(0,-1,1,12),3,p(0,1,3),p(0,3,3),p(0,2,3),1,1,-1)

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
     &      ( pdf(2*sign1,1,3)*pdf(2*sign2,2,3) + pdf(4*sign1,1,3)*pdf(4*sign2,2,3))/2d0*polcolgq
         res(3,3) = diag(3) * dipole(1) * 
     &      ( pdf(2*sign1,1,3)*pdf(4*sign2,2,3) + pdf(4*sign1,1,3)*pdf(2*sign2,2,3))*polcolgq
         res(4,3) = diag(4) * dipole(1) * 
     &      ( pdf(2*sign1,1,3)*pdf(3*sign2,2,3) + pdf(4*sign1,1,3)*pdf(1*sign2,2,3))*polcolgq
         res(6,3) = diag(6) * dipole(2) * 
     &      ( pdf(2*sign1,1,3)*pdf(0*sign2,2,3) + pdf(4*sign1,1,3)*pdf(0*sign2,2,3))*polcolgq
         res(7,3) = diag(7) * dipole(1) * 
     &      ( pdf(2*sign1,1,3)*pdf(-3*sign2,2,3) + pdf(4*sign1,1,3)*pdf(-1*sign2,2,3))*polcolgq
         res(8,3) = diag(8) * dipole(1) * 
     &      ( pdf(2*sign1,1,3)*pdf(-4*sign2,2,3) + pdf(4*sign1,1,3)*pdf(-2*sign2,2,3))*polcolgq
         res(10,3) = diag(10) * dipole(1) * 
     &      ( pdf(2*sign1,1,3)*pdf(-1*sign2,2,3) + pdf(4*sign1,1,3)*pdf(-3*sign2,2,3))*polcolgq
         res(12,3) = diag(12) * dipole(1) * 
     &      ( pdf(2*sign1,1,3)*pdf(-2*sign2,2,3) + pdf(4*sign1,1,3)*pdf(-4*sign2,2,3))*polcolgq
         res(16,3) = diag(16) * dipole(3) * 
     &      ( pdf(2*sign1,1,3)*pdf(-2*sign2,2,3) + pdf(4*sign1,1,3)*pdf(-4*sign2,2,3))/2d0*polcolqq
         res(17,3) = diag(17) * dipole(4) * 
     &      ( pdf(1*sign1,1,3)*pdf(2*sign2,2,3) + pdf(3*sign1,1,3)*pdf(4*sign2,2,3))*polcolgq
         res(18,3) = diag(18) * dipole(4) * 
     &      ( pdf(1*sign1,1,3)*pdf(1*sign2,2,3) + pdf(3*sign1,1,3)*pdf(3*sign2,2,3))/2d0*polcolgq
         res(20,3) = diag(20) * dipole(4) * 
     &      ( pdf(1*sign1,1,3)*pdf(4*sign2,2,3) + pdf(3*sign1,1,3)*pdf(2*sign2,2,3))*polcolgq
         res(21,3) = diag(21) * dipole(4) * 
     &      ( pdf(1*sign1,1,3)*pdf(3*sign2,2,3) + pdf(3*sign1,1,3)*pdf(1*sign2,2,3))*polcolgq
         res(22,3) = diag(22) * dipole(5) * 
     &      ( pdf(1*sign1,1,3)*pdf(0*sign2,2,3) + pdf(3*sign1,1,3)*pdf(0*sign2,2,3))*polcolgq
         res(24,3) = diag(24) * dipole(4) * 
     &      ( pdf(1*sign1,1,3)*pdf(-3*sign2,2,3) + pdf(3*sign1,1,3)*pdf(-1*sign2,2,3))*polcolgq
         res(25,3) = diag(25) * dipole(4) * 
     &      ( pdf(1*sign1,1,3)*pdf(-4*sign2,2,3) + pdf(3*sign1,1,3)*pdf(-2*sign2,2,3))*polcolgq
         res(27,3) = diag(27) * dipole(4) * 
     &      ( pdf(1*sign1,1,3)*pdf(-1*sign2,2,3) + pdf(3*sign1,1,3)*pdf(-3*sign2,2,3))*polcolgq
         res(30,3) = diag(30) * dipole(6) * 
     &      ( pdf(1*sign1,1,3)*pdf(-1*sign2,2,3) + pdf(3*sign1,1,3)*pdf(-3*sign2,2,3))/2d0*polcolqq
         res(31,3) = diag(31) * dipole(4) * 
     &      ( pdf(1*sign1,1,3)*pdf(-2*sign2,2,3) + pdf(3*sign1,1,3)*pdf(-4*sign2,2,3))*polcolgq
         res(33,3) = diag(33) * dipole(7) * 
     &      ( pdf(0*sign1,1,3)*pdf(2*sign2,2,3) + pdf(0*sign1,1,3)*pdf(4*sign2,2,3))*polcolgq
         res(34,3) = diag(34) * dipole(8) * 
     &      ( pdf(0*sign1,1,3)*pdf(1*sign2,2,3) + pdf(0*sign1,1,3)*pdf(3*sign2,2,3))*polcolgq
         res(35,3) = diag(35) * dipole(9) * 
     &      ( pdf(0*sign1,1,3)*pdf(0*sign2,2,3) + pdf(0*sign1,1,3)*pdf(0*sign2,2,3))*polcolgq
         res(36,3) = diag(36) * dipole(10) * 
     &      ( pdf(0*sign1,1,3)*pdf(0*sign2,2,3) + pdf(0*sign1,1,3)*pdf(0*sign2,2,3))*polcolgq
         res(37,3) = diag(37) * dipole(11) * 
     &      ( pdf(0*sign1,1,3)*pdf(-1*sign2,2,3) + pdf(0*sign1,1,3)*pdf(-3*sign2,2,3))*polcolgq
         res(38,3) = diag(38) * dipole(12) * 
     &      ( pdf(0*sign1,1,3)*pdf(-2*sign2,2,3) + pdf(0*sign1,1,3)*pdf(-4*sign2,2,3))*polcolgq
         res(48,3) = diag(48) * dipole(13) * 
     &      ( pdf(-1*sign1,1,3)*pdf(1*sign2,2,3) + pdf(-3*sign1,1,3)*pdf(3*sign2,2,3))/2d0*polcolqq
         res(49,3) = diag(49) * dipole(14) * 
     &      ( pdf(-1*sign1,1,3)*pdf(0*sign2,2,3) + pdf(-3*sign1,1,3)*pdf(0*sign2,2,3))*polcolgq
         res(50,3) = diag(50) * dipole(15) * 
     &      ( pdf(-1*sign1,1,3)*pdf(-1*sign2,2,3) + pdf(-3*sign1,1,3)*pdf(-3*sign2,2,3))/2d0*polcolgq
         res(62,3) = diag(62) * dipole(16) * 
     &      ( pdf(-2*sign1,1,3)*pdf(2*sign2,2,3) + pdf(-4*sign1,1,3)*pdf(4*sign2,2,3))/2d0*polcolqq
         res(65,3) = diag(65) * dipole(17) * 
     &      ( pdf(-2*sign1,1,3)*pdf(0*sign2,2,3) + pdf(-4*sign1,1,3)*pdf(0*sign2,2,3))*polcolgq
         res(66,3) = diag(66) * dipole(18) * 
     &      ( pdf(-2*sign1,1,3)*pdf(-1*sign2,2,3) + pdf(-4*sign1,1,3)*pdf(-3*sign2,2,3))*polcolgq
         res(67,3) = diag(67) * dipole(18) * 
     &      ( pdf(-2*sign1,1,3)*pdf(-2*sign2,2,3) + pdf(-4*sign1,1,3)*pdf(-4*sign2,2,3))/2d0*polcolgq
       endif 

       if (lok(4)) then 
         call qqgWW_curr(borncurr(0,-1,1,13),4,p(0,3,4),p(0,1,4),p(0,2,4),-1,-1,-1)
         call qqgWW_curr(borncurr(0,-1,1,14),4,p(0,2,4),p(0,1,4),p(0,3,4),1,-1,1)
         call qqgWW_curr(borncurr(0,-1,1,15),4,p(0,3,4),p(0,2,4),p(0,1,4),-1,-1,-1)
         call qqgWW_curr(borncurr(0,-1,1,16),4,p(0,2,4),p(0,3,4),p(0,1,4),1,1,-1)
         call qqgWW_curr(borncurr(0,-1,1,17),4,p(0,1,4),p(0,2,4),p(0,3,4),1,-1,1)
         call qqgWW_curr(borncurr(0,-1,1,18),4,p(0,1,4),p(0,3,4),p(0,2,4),1,1,-1)

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
     &      ( pdf(2*sign1,1,4)*pdf(2*sign2,2,4) + pdf(4*sign1,1,4)*pdf(4*sign2,2,4))/2d0*polcolgq
         res(2,4) = diag(2) * dipole(2) * 
     &      ( pdf(2*sign1,1,4)*pdf(1*sign2,2,4) + pdf(4*sign1,1,4)*pdf(3*sign2,2,4))*polcolgq
         res(6,4) = diag(6) * dipole(3) * 
     &      ( pdf(2*sign1,1,4)*pdf(0*sign2,2,4) + pdf(4*sign1,1,4)*pdf(0*sign2,2,4))*polcolgq
         res(16,4) = diag(16) * dipole(4) * 
     &      ( pdf(2*sign1,1,4)*pdf(-2*sign2,2,4) + pdf(4*sign1,1,4)*pdf(-4*sign2,2,4))/2d0*polcolqq
         res(18,4) = diag(18) * dipole(2) * 
     &      ( pdf(1*sign1,1,4)*pdf(1*sign2,2,4) + pdf(3*sign1,1,4)*pdf(3*sign2,2,4))/2d0*polcolgq
         res(22,4) = diag(22) * dipole(5) * 
     &      ( pdf(1*sign1,1,4)*pdf(0*sign2,2,4) + pdf(3*sign1,1,4)*pdf(0*sign2,2,4))*polcolgq
         res(30,4) = diag(30) * dipole(6) * 
     &      ( pdf(1*sign1,1,4)*pdf(-1*sign2,2,4) + pdf(3*sign1,1,4)*pdf(-3*sign2,2,4))/2d0*polcolqq
         res(33,4) = diag(33) * dipole(7) * 
     &      ( pdf(0*sign1,1,4)*pdf(2*sign2,2,4) + pdf(0*sign1,1,4)*pdf(4*sign2,2,4))*polcolgq
         res(34,4) = diag(34) * dipole(8) * 
     &      ( pdf(0*sign1,1,4)*pdf(1*sign2,2,4) + pdf(0*sign1,1,4)*pdf(3*sign2,2,4))*polcolgq
         res(35,4) = diag(35) * dipole(9) * 
     &      ( pdf(0*sign1,1,4)*pdf(0*sign2,2,4) + pdf(0*sign1,1,4)*pdf(0*sign2,2,4))*polcolgq
         res(36,4) = diag(36) * dipole(10) * 
     &      ( pdf(0*sign1,1,4)*pdf(0*sign2,2,4) + pdf(0*sign1,1,4)*pdf(0*sign2,2,4))*polcolgq
         res(37,4) = diag(37) * dipole(11) * 
     &      ( pdf(0*sign1,1,4)*pdf(-1*sign2,2,4) + pdf(0*sign1,1,4)*pdf(-3*sign2,2,4))*polcolgq
         res(38,4) = diag(38) * dipole(12) * 
     &      ( pdf(0*sign1,1,4)*pdf(-2*sign2,2,4) + pdf(0*sign1,1,4)*pdf(-4*sign2,2,4))*polcolgq
         res(39,4) = diag(39) * dipole(1) * 
     &      ( pdf(-1*sign1,1,4)*pdf(4*sign2,2,4) + pdf(-3*sign1,1,4)*pdf(2*sign2,2,4))*polcolgq
         res(41,4) = diag(41) * dipole(2) * 
     &      ( pdf(-1*sign1,1,4)*pdf(3*sign2,2,4) + pdf(-3*sign1,1,4)*pdf(1*sign2,2,4))*polcolgq
         res(43,4) = diag(43) * dipole(1) * 
     &      ( pdf(-1*sign1,1,4)*pdf(2*sign2,2,4) + pdf(-3*sign1,1,4)*pdf(4*sign2,2,4))*polcolgq
         res(47,4) = diag(47) * dipole(2) * 
     &      ( pdf(-1*sign1,1,4)*pdf(1*sign2,2,4) + pdf(-3*sign1,1,4)*pdf(3*sign2,2,4))*polcolgq
         res(48,4) = diag(48) * dipole(13) * 
     &      ( pdf(-1*sign1,1,4)*pdf(1*sign2,2,4) + pdf(-3*sign1,1,4)*pdf(3*sign2,2,4))/2d0*polcolqq
         res(49,4) = diag(49) * dipole(14) * 
     &      ( pdf(-1*sign1,1,4)*pdf(0*sign2,2,4) + pdf(-3*sign1,1,4)*pdf(0*sign2,2,4))*polcolgq
         res(50,4) = diag(50) * dipole(15) * 
     &      ( pdf(-1*sign1,1,4)*pdf(-1*sign2,2,4) + pdf(-3*sign1,1,4)*pdf(-3*sign2,2,4))/2d0*polcolgq
         res(51,4) = diag(51) * dipole(16) * 
     &      ( pdf(-1*sign1,1,4)*pdf(-2*sign2,2,4) + pdf(-3*sign1,1,4)*pdf(-4*sign2,2,4))*polcolgq
         res(52,4) = diag(52) * dipole(15) * 
     &      ( pdf(-1*sign1,1,4)*pdf(-3*sign2,2,4) + pdf(-3*sign1,1,4)*pdf(-1*sign2,2,4))*polcolgq
         res(53,4) = diag(53) * dipole(16) * 
     &      ( pdf(-1*sign1,1,4)*pdf(-4*sign2,2,4) + pdf(-3*sign1,1,4)*pdf(-2*sign2,2,4))*polcolgq
         res(55,4) = diag(55) * dipole(1) * 
     &      ( pdf(-2*sign1,1,4)*pdf(4*sign2,2,4) + pdf(-4*sign1,1,4)*pdf(2*sign2,2,4))*polcolgq
         res(57,4) = diag(57) * dipole(2) * 
     &      ( pdf(-2*sign1,1,4)*pdf(3*sign2,2,4) + pdf(-4*sign1,1,4)*pdf(1*sign2,2,4))*polcolgq
         res(60,4) = diag(60) * dipole(1) * 
     &      ( pdf(-2*sign1,1,4)*pdf(2*sign2,2,4) + pdf(-4*sign1,1,4)*pdf(4*sign2,2,4))*polcolgq
         res(62,4) = diag(62) * dipole(17) * 
     &      ( pdf(-2*sign1,1,4)*pdf(2*sign2,2,4) + pdf(-4*sign1,1,4)*pdf(4*sign2,2,4))/2d0*polcolqq
         res(64,4) = diag(64) * dipole(2) * 
     &      ( pdf(-2*sign1,1,4)*pdf(1*sign2,2,4) + pdf(-4*sign1,1,4)*pdf(3*sign2,2,4))*polcolgq
         res(65,4) = diag(65) * dipole(18) * 
     &      ( pdf(-2*sign1,1,4)*pdf(0*sign2,2,4) + pdf(-4*sign1,1,4)*pdf(0*sign2,2,4))*polcolgq
         res(67,4) = diag(67) * dipole(16) * 
     &      ( pdf(-2*sign1,1,4)*pdf(-2*sign2,2,4) + pdf(-4*sign1,1,4)*pdf(-4*sign2,2,4))/2d0*polcolgq
         res(69,4) = diag(69) * dipole(15) * 
     &      ( pdf(-2*sign1,1,4)*pdf(-3*sign2,2,4) + pdf(-4*sign1,1,4)*pdf(-1*sign2,2,4))*polcolgq
         res(70,4) = diag(70) * dipole(16) * 
     &      ( pdf(-2*sign1,1,4)*pdf(-4*sign2,2,4) + pdf(-4*sign1,1,4)*pdf(-2*sign2,2,4))*polcolgq
       endif 

       if (lok(5)) then 
         call qqgWW_curr(borncurr(0,-1,1,19),5,p(0,3,5),p(0,1,5),p(0,2,5),-1,-1,-1)
         call qqgWW_curr(borncurr(0,-1,1,20),5,p(0,2,5),p(0,1,5),p(0,3,5),1,-1,1)
         call qqgWW_curr(borncurr(0,-1,1,21),5,p(0,3,5),p(0,2,5),p(0,1,5),-1,-1,-1)
         call qqgWW_curr(borncurr(0,-1,1,22),5,p(0,2,5),p(0,3,5),p(0,1,5),1,1,-1)
         call qqgWW_curr(borncurr(0,-1,1,23),5,p(0,1,5),p(0,2,5),p(0,3,5),1,-1,1)
         call qqgWW_curr(borncurr(0,-1,1,24),5,p(0,1,5),p(0,3,5),p(0,2,5),1,1,-1)

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
     &      ( pdf(2*sign1,1,5)*pdf(2*sign2,2,5) + pdf(4*sign1,1,5)*pdf(4*sign2,2,5))/2d0*polcolgq
         res(3,5) = diag(3) * dipole(1) * 
     &      ( pdf(2*sign1,1,5)*pdf(4*sign2,2,5) + pdf(4*sign1,1,5)*pdf(2*sign2,2,5))*polcolgq
         res(4,5) = diag(4) * dipole(1) * 
     &      ( pdf(2*sign1,1,5)*pdf(3*sign2,2,5) + pdf(4*sign1,1,5)*pdf(1*sign2,2,5))*polcolgq
         res(6,5) = diag(6) * dipole(2) * 
     &      ( pdf(2*sign1,1,5)*pdf(0*sign2,2,5) + pdf(4*sign1,1,5)*pdf(0*sign2,2,5))*polcolgq
         res(7,5) = diag(7) * dipole(1) * 
     &      ( pdf(2*sign1,1,5)*pdf(-3*sign2,2,5) + pdf(4*sign1,1,5)*pdf(-1*sign2,2,5))*polcolgq
         res(8,5) = diag(8) * dipole(1) * 
     &      ( pdf(2*sign1,1,5)*pdf(-4*sign2,2,5) + pdf(4*sign1,1,5)*pdf(-2*sign2,2,5))*polcolgq
         res(10,5) = diag(10) * dipole(1) * 
     &      ( pdf(2*sign1,1,5)*pdf(-1*sign2,2,5) + pdf(4*sign1,1,5)*pdf(-3*sign2,2,5))*polcolgq
         res(12,5) = diag(12) * dipole(1) * 
     &      ( pdf(2*sign1,1,5)*pdf(-2*sign2,2,5) + pdf(4*sign1,1,5)*pdf(-4*sign2,2,5))*polcolgq
         res(16,5) = diag(16) * dipole(3) * 
     &      ( pdf(2*sign1,1,5)*pdf(-2*sign2,2,5) + pdf(4*sign1,1,5)*pdf(-4*sign2,2,5))/2d0*polcolqq
         res(17,5) = diag(17) * dipole(4) * 
     &      ( pdf(1*sign1,1,5)*pdf(2*sign2,2,5) + pdf(3*sign1,1,5)*pdf(4*sign2,2,5))*polcolgq
         res(18,5) = diag(18) * dipole(4) * 
     &      ( pdf(1*sign1,1,5)*pdf(1*sign2,2,5) + pdf(3*sign1,1,5)*pdf(3*sign2,2,5))/2d0*polcolgq
         res(20,5) = diag(20) * dipole(4) * 
     &      ( pdf(1*sign1,1,5)*pdf(4*sign2,2,5) + pdf(3*sign1,1,5)*pdf(2*sign2,2,5))*polcolgq
         res(21,5) = diag(21) * dipole(4) * 
     &      ( pdf(1*sign1,1,5)*pdf(3*sign2,2,5) + pdf(3*sign1,1,5)*pdf(1*sign2,2,5))*polcolgq
         res(22,5) = diag(22) * dipole(5) * 
     &      ( pdf(1*sign1,1,5)*pdf(0*sign2,2,5) + pdf(3*sign1,1,5)*pdf(0*sign2,2,5))*polcolgq
         res(24,5) = diag(24) * dipole(4) * 
     &      ( pdf(1*sign1,1,5)*pdf(-3*sign2,2,5) + pdf(3*sign1,1,5)*pdf(-1*sign2,2,5))*polcolgq
         res(25,5) = diag(25) * dipole(4) * 
     &      ( pdf(1*sign1,1,5)*pdf(-4*sign2,2,5) + pdf(3*sign1,1,5)*pdf(-2*sign2,2,5))*polcolgq
         res(27,5) = diag(27) * dipole(4) * 
     &      ( pdf(1*sign1,1,5)*pdf(-1*sign2,2,5) + pdf(3*sign1,1,5)*pdf(-3*sign2,2,5))*polcolgq
         res(30,5) = diag(30) * dipole(6) * 
     &      ( pdf(1*sign1,1,5)*pdf(-1*sign2,2,5) + pdf(3*sign1,1,5)*pdf(-3*sign2,2,5))/2d0*polcolqq
         res(31,5) = diag(31) * dipole(4) * 
     &      ( pdf(1*sign1,1,5)*pdf(-2*sign2,2,5) + pdf(3*sign1,1,5)*pdf(-4*sign2,2,5))*polcolgq
         res(33,5) = diag(33) * dipole(7) * 
     &      ( pdf(0*sign1,1,5)*pdf(2*sign2,2,5) + pdf(0*sign1,1,5)*pdf(4*sign2,2,5))*polcolgq
         res(34,5) = diag(34) * dipole(8) * 
     &      ( pdf(0*sign1,1,5)*pdf(1*sign2,2,5) + pdf(0*sign1,1,5)*pdf(3*sign2,2,5))*polcolgq
         res(35,5) = diag(35) * dipole(9) * 
     &      ( pdf(0*sign1,1,5)*pdf(0*sign2,2,5) + pdf(0*sign1,1,5)*pdf(0*sign2,2,5))*polcolgq
         res(36,5) = diag(36) * dipole(10) * 
     &      ( pdf(0*sign1,1,5)*pdf(0*sign2,2,5) + pdf(0*sign1,1,5)*pdf(0*sign2,2,5))*polcolgq
         res(37,5) = diag(37) * dipole(11) * 
     &      ( pdf(0*sign1,1,5)*pdf(-1*sign2,2,5) + pdf(0*sign1,1,5)*pdf(-3*sign2,2,5))*polcolgq
         res(38,5) = diag(38) * dipole(12) * 
     &      ( pdf(0*sign1,1,5)*pdf(-2*sign2,2,5) + pdf(0*sign1,1,5)*pdf(-4*sign2,2,5))*polcolgq
         res(48,5) = diag(48) * dipole(13) * 
     &      ( pdf(-1*sign1,1,5)*pdf(1*sign2,2,5) + pdf(-3*sign1,1,5)*pdf(3*sign2,2,5))/2d0*polcolqq
         res(49,5) = diag(49) * dipole(14) * 
     &      ( pdf(-1*sign1,1,5)*pdf(0*sign2,2,5) + pdf(-3*sign1,1,5)*pdf(0*sign2,2,5))*polcolgq
         res(50,5) = diag(50) * dipole(15) * 
     &      ( pdf(-1*sign1,1,5)*pdf(-1*sign2,2,5) + pdf(-3*sign1,1,5)*pdf(-3*sign2,2,5))/2d0*polcolgq
         res(62,5) = diag(62) * dipole(16) * 
     &      ( pdf(-2*sign1,1,5)*pdf(2*sign2,2,5) + pdf(-4*sign1,1,5)*pdf(4*sign2,2,5))/2d0*polcolqq
         res(65,5) = diag(65) * dipole(17) * 
     &      ( pdf(-2*sign1,1,5)*pdf(0*sign2,2,5) + pdf(-4*sign1,1,5)*pdf(0*sign2,2,5))*polcolgq
         res(66,5) = diag(66) * dipole(18) * 
     &      ( pdf(-2*sign1,1,5)*pdf(-1*sign2,2,5) + pdf(-4*sign1,1,5)*pdf(-3*sign2,2,5))*polcolgq
         res(67,5) = diag(67) * dipole(18) * 
     &      ( pdf(-2*sign1,1,5)*pdf(-2*sign2,2,5) + pdf(-4*sign1,1,5)*pdf(-4*sign2,2,5))/2d0*polcolgq
       endif 

       if (lok(6)) then 
         call qqgWW_curr(borncurr(0,-1,1,25),6,p(0,3,6),p(0,1,6),p(0,2,6),-1,-1,-1)
         call qqgWW_curr(borncurr(0,-1,1,26),6,p(0,2,6),p(0,1,6),p(0,3,6),1,-1,1)
         call qqgWW_curr(borncurr(0,-1,1,27),6,p(0,3,6),p(0,2,6),p(0,1,6),-1,-1,-1)
         call qqgWW_curr(borncurr(0,-1,1,28),6,p(0,2,6),p(0,3,6),p(0,1,6),1,1,-1)
         call qqgWW_curr(borncurr(0,-1,1,29),6,p(0,1,6),p(0,2,6),p(0,3,6),1,-1,1)

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
     &      ( pdf(2*sign1,1,6)*pdf(2*sign2,2,6) + pdf(4*sign1,1,6)*pdf(4*sign2,2,6))/2d0*polcolgq
         res(3,6) = diag(3) * dipole(1) * 
     &      ( pdf(2*sign1,1,6)*pdf(4*sign2,2,6) + pdf(4*sign1,1,6)*pdf(2*sign2,2,6))*polcolgq
         res(4,6) = diag(4) * dipole(2) * 
     &      ( pdf(2*sign1,1,6)*pdf(3*sign2,2,6) + pdf(4*sign1,1,6)*pdf(1*sign2,2,6))*polcolgq
         res(7,6) = diag(7) * dipole(3) * 
     &      ( pdf(2*sign1,1,6)*pdf(-3*sign2,2,6) + pdf(4*sign1,1,6)*pdf(-1*sign2,2,6))*polcolgq
         res(8,6) = diag(8) * dipole(4) * 
     &      ( pdf(2*sign1,1,6)*pdf(-4*sign2,2,6) + pdf(4*sign1,1,6)*pdf(-2*sign2,2,6))*polcolgq
         res(10,6) = diag(10) * dipole(3) * 
     &      ( pdf(2*sign1,1,6)*pdf(-1*sign2,2,6) + pdf(4*sign1,1,6)*pdf(-3*sign2,2,6))*polcolgq
         res(12,6) = diag(12) * dipole(4) * 
     &      ( pdf(2*sign1,1,6)*pdf(-2*sign2,2,6) + pdf(4*sign1,1,6)*pdf(-4*sign2,2,6))*polcolgq
         res(16,6) = diag(16) * dipole(5) * 
     &      ( pdf(2*sign1,1,6)*pdf(-2*sign2,2,6) + pdf(4*sign1,1,6)*pdf(-4*sign2,2,6))/2d0*polcolqq
         res(17,6) = diag(17) * dipole(1) * 
     &      ( pdf(1*sign1,1,6)*pdf(2*sign2,2,6) + pdf(3*sign1,1,6)*pdf(4*sign2,2,6))*polcolgq
         res(18,6) = diag(18) * dipole(2) * 
     &      ( pdf(1*sign1,1,6)*pdf(1*sign2,2,6) + pdf(3*sign1,1,6)*pdf(3*sign2,2,6))/2d0*polcolgq
         res(20,6) = diag(20) * dipole(1) * 
     &      ( pdf(1*sign1,1,6)*pdf(4*sign2,2,6) + pdf(3*sign1,1,6)*pdf(2*sign2,2,6))*polcolgq
         res(21,6) = diag(21) * dipole(2) * 
     &      ( pdf(1*sign1,1,6)*pdf(3*sign2,2,6) + pdf(3*sign1,1,6)*pdf(1*sign2,2,6))*polcolgq
         res(24,6) = diag(24) * dipole(3) * 
     &      ( pdf(1*sign1,1,6)*pdf(-3*sign2,2,6) + pdf(3*sign1,1,6)*pdf(-1*sign2,2,6))*polcolgq
         res(25,6) = diag(25) * dipole(4) * 
     &      ( pdf(1*sign1,1,6)*pdf(-4*sign2,2,6) + pdf(3*sign1,1,6)*pdf(-2*sign2,2,6))*polcolgq
         res(27,6) = diag(27) * dipole(3) * 
     &      ( pdf(1*sign1,1,6)*pdf(-1*sign2,2,6) + pdf(3*sign1,1,6)*pdf(-3*sign2,2,6))*polcolgq
         res(30,6) = diag(30) * dipole(6) * 
     &      ( pdf(1*sign1,1,6)*pdf(-1*sign2,2,6) + pdf(3*sign1,1,6)*pdf(-3*sign2,2,6))/2d0*polcolqq
         res(31,6) = diag(31) * dipole(4) * 
     &      ( pdf(1*sign1,1,6)*pdf(-2*sign2,2,6) + pdf(3*sign1,1,6)*pdf(-4*sign2,2,6))*polcolgq
         res(33,6) = diag(33) * dipole(7) * 
     &      ( pdf(0*sign1,1,6)*pdf(2*sign2,2,6) + pdf(0*sign1,1,6)*pdf(4*sign2,2,6))*polcolqq
         res(34,6) = diag(34) * dipole(8) * 
     &      ( pdf(0*sign1,1,6)*pdf(1*sign2,2,6) + pdf(0*sign1,1,6)*pdf(3*sign2,2,6))*polcolqq
         res(35,6) = diag(35) * dipole(9) * 
     &      ( pdf(0*sign1,1,6)*pdf(0*sign2,2,6) + pdf(0*sign1,1,6)*pdf(0*sign2,2,6))*polcolgq
         res(36,6) = diag(36) * dipole(10) * 
     &      ( pdf(0*sign1,1,6)*pdf(0*sign2,2,6) + pdf(0*sign1,1,6)*pdf(0*sign2,2,6))*polcolgq
         res(37,6) = diag(37) * dipole(11) * 
     &      ( pdf(0*sign1,1,6)*pdf(-1*sign2,2,6) + pdf(0*sign1,1,6)*pdf(-3*sign2,2,6))*polcolqq
         res(38,6) = diag(38) * dipole(12) * 
     &      ( pdf(0*sign1,1,6)*pdf(-2*sign2,2,6) + pdf(0*sign1,1,6)*pdf(-4*sign2,2,6))*polcolqq
         res(48,6) = diag(48) * dipole(13) * 
     &      ( pdf(-1*sign1,1,6)*pdf(1*sign2,2,6) + pdf(-3*sign1,1,6)*pdf(3*sign2,2,6))/2d0*polcolqq
         res(50,6) = diag(50) * dipole(3) * 
     &      ( pdf(-1*sign1,1,6)*pdf(-1*sign2,2,6) + pdf(-3*sign1,1,6)*pdf(-3*sign2,2,6))/2d0*polcolgq
         res(62,6) = diag(62) * dipole(14) * 
     &      ( pdf(-2*sign1,1,6)*pdf(2*sign2,2,6) + pdf(-4*sign1,1,6)*pdf(4*sign2,2,6))/2d0*polcolqq
         res(66,6) = diag(66) * dipole(3) * 
     &      ( pdf(-2*sign1,1,6)*pdf(-1*sign2,2,6) + pdf(-4*sign1,1,6)*pdf(-3*sign2,2,6))*polcolgq
         res(67,6) = diag(67) * dipole(4) * 
     &      ( pdf(-2*sign1,1,6)*pdf(-2*sign2,2,6) + pdf(-4*sign1,1,6)*pdf(-4*sign2,2,6))/2d0*polcolgq
       endif 

       if (lok(7)) then 
         call qqgWW_curr(borncurr(0,-1,1,30),7,p(0,3,7),p(0,1,7),p(0,2,7),-1,-1,-1)
         call qqgWW_curr(borncurr(0,-1,1,31),7,p(0,2,7),p(0,1,7),p(0,3,7),1,-1,1)
         call qqgWW_curr(borncurr(0,-1,1,32),7,p(0,3,7),p(0,2,7),p(0,1,7),-1,-1,-1)
         call qqgWW_curr(borncurr(0,-1,1,33),7,p(0,1,7),p(0,2,7),p(0,3,7),1,-1,1)
         call qqgWW_curr(borncurr(0,-1,1,34),7,p(0,1,7),p(0,3,7),p(0,2,7),1,1,-1)

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
     &      ( pdf(2*sign1,1,7)*pdf(2*sign2,2,7) + pdf(4*sign1,1,7)*pdf(4*sign2,2,7))/2d0*polcolgq
         res(2,7) = diag(2) * dipole(1) * 
     &      ( pdf(2*sign1,1,7)*pdf(1*sign2,2,7) + pdf(4*sign1,1,7)*pdf(3*sign2,2,7))*polcolgq
         res(6,7) = diag(6) * dipole(2) * 
     &      ( pdf(2*sign1,1,7)*pdf(0*sign2,2,7) + pdf(4*sign1,1,7)*pdf(0*sign2,2,7))*polcolqq
         res(16,7) = diag(16) * dipole(3) * 
     &      ( pdf(2*sign1,1,7)*pdf(-2*sign2,2,7) + pdf(4*sign1,1,7)*pdf(-4*sign2,2,7))/2d0*polcolqq
         res(18,7) = diag(18) * dipole(4) * 
     &      ( pdf(1*sign1,1,7)*pdf(1*sign2,2,7) + pdf(3*sign1,1,7)*pdf(3*sign2,2,7))/2d0*polcolgq
         res(22,7) = diag(22) * dipole(5) * 
     &      ( pdf(1*sign1,1,7)*pdf(0*sign2,2,7) + pdf(3*sign1,1,7)*pdf(0*sign2,2,7))*polcolqq
         res(30,7) = diag(30) * dipole(6) * 
     &      ( pdf(1*sign1,1,7)*pdf(-1*sign2,2,7) + pdf(3*sign1,1,7)*pdf(-3*sign2,2,7))/2d0*polcolqq
         res(35,7) = diag(35) * dipole(7) * 
     &      ( pdf(0*sign1,1,7)*pdf(0*sign2,2,7) + pdf(0*sign1,1,7)*pdf(0*sign2,2,7))*polcolgq
         res(36,7) = diag(36) * dipole(8) * 
     &      ( pdf(0*sign1,1,7)*pdf(0*sign2,2,7) + pdf(0*sign1,1,7)*pdf(0*sign2,2,7))*polcolgq
         res(39,7) = diag(39) * dipole(9) * 
     &      ( pdf(-1*sign1,1,7)*pdf(4*sign2,2,7) + pdf(-3*sign1,1,7)*pdf(2*sign2,2,7))*polcolgq
         res(41,7) = diag(41) * dipole(9) * 
     &      ( pdf(-1*sign1,1,7)*pdf(3*sign2,2,7) + pdf(-3*sign1,1,7)*pdf(1*sign2,2,7))*polcolgq
         res(43,7) = diag(43) * dipole(9) * 
     &      ( pdf(-1*sign1,1,7)*pdf(2*sign2,2,7) + pdf(-3*sign1,1,7)*pdf(4*sign2,2,7))*polcolgq
         res(47,7) = diag(47) * dipole(9) * 
     &      ( pdf(-1*sign1,1,7)*pdf(1*sign2,2,7) + pdf(-3*sign1,1,7)*pdf(3*sign2,2,7))*polcolgq
         res(48,7) = diag(48) * dipole(10) * 
     &      ( pdf(-1*sign1,1,7)*pdf(1*sign2,2,7) + pdf(-3*sign1,1,7)*pdf(3*sign2,2,7))/2d0*polcolqq
         res(49,7) = diag(49) * dipole(11) * 
     &      ( pdf(-1*sign1,1,7)*pdf(0*sign2,2,7) + pdf(-3*sign1,1,7)*pdf(0*sign2,2,7))*polcolqq
         res(50,7) = diag(50) * dipole(9) * 
     &      ( pdf(-1*sign1,1,7)*pdf(-1*sign2,2,7) + pdf(-3*sign1,1,7)*pdf(-3*sign2,2,7))/2d0*polcolgq
         res(51,7) = diag(51) * dipole(9) * 
     &      ( pdf(-1*sign1,1,7)*pdf(-2*sign2,2,7) + pdf(-3*sign1,1,7)*pdf(-4*sign2,2,7))*polcolgq
         res(52,7) = diag(52) * dipole(9) * 
     &      ( pdf(-1*sign1,1,7)*pdf(-3*sign2,2,7) + pdf(-3*sign1,1,7)*pdf(-1*sign2,2,7))*polcolgq
         res(53,7) = diag(53) * dipole(9) * 
     &      ( pdf(-1*sign1,1,7)*pdf(-4*sign2,2,7) + pdf(-3*sign1,1,7)*pdf(-2*sign2,2,7))*polcolgq
         res(55,7) = diag(55) * dipole(12) * 
     &      ( pdf(-2*sign1,1,7)*pdf(4*sign2,2,7) + pdf(-4*sign1,1,7)*pdf(2*sign2,2,7))*polcolgq
         res(57,7) = diag(57) * dipole(12) * 
     &      ( pdf(-2*sign1,1,7)*pdf(3*sign2,2,7) + pdf(-4*sign1,1,7)*pdf(1*sign2,2,7))*polcolgq
         res(60,7) = diag(60) * dipole(12) * 
     &      ( pdf(-2*sign1,1,7)*pdf(2*sign2,2,7) + pdf(-4*sign1,1,7)*pdf(4*sign2,2,7))*polcolgq
         res(62,7) = diag(62) * dipole(13) * 
     &      ( pdf(-2*sign1,1,7)*pdf(2*sign2,2,7) + pdf(-4*sign1,1,7)*pdf(4*sign2,2,7))/2d0*polcolqq
         res(64,7) = diag(64) * dipole(12) * 
     &      ( pdf(-2*sign1,1,7)*pdf(1*sign2,2,7) + pdf(-4*sign1,1,7)*pdf(3*sign2,2,7))*polcolgq
         res(65,7) = diag(65) * dipole(14) * 
     &      ( pdf(-2*sign1,1,7)*pdf(0*sign2,2,7) + pdf(-4*sign1,1,7)*pdf(0*sign2,2,7))*polcolqq
         res(67,7) = diag(67) * dipole(12) * 
     &      ( pdf(-2*sign1,1,7)*pdf(-2*sign2,2,7) + pdf(-4*sign1,1,7)*pdf(-4*sign2,2,7))/2d0*polcolgq
         res(69,7) = diag(69) * dipole(12) * 
     &      ( pdf(-2*sign1,1,7)*pdf(-3*sign2,2,7) + pdf(-4*sign1,1,7)*pdf(-1*sign2,2,7))*polcolgq
         res(70,7) = diag(70) * dipole(12) * 
     &      ( pdf(-2*sign1,1,7)*pdf(-4*sign2,2,7) + pdf(-4*sign1,1,7)*pdf(-2*sign2,2,7))*polcolgq
       endif 

       if (lok(8)) then 
         call qqgWW_curr(borncurr(0,-1,1,35),8,p(0,3,8),p(0,1,8),p(0,2,8),-1,-1,-1)
         call qqgWW_curr(borncurr(0,-1,1,36),8,p(0,2,8),p(0,1,8),p(0,3,8),1,-1,1)
         call qqgWW_curr(borncurr(0,-1,1,37),8,p(0,3,8),p(0,2,8),p(0,1,8),-1,-1,-1)
         call qqgWW_curr(borncurr(0,-1,1,38),8,p(0,2,8),p(0,3,8),p(0,1,8),1,1,-1)
         call qqgWW_curr(borncurr(0,-1,1,39),8,p(0,1,8),p(0,2,8),p(0,3,8),1,-1,1)

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
     &      ( pdf(2*sign1,1,8)*pdf(2*sign2,2,8) + pdf(4*sign1,1,8)*pdf(4*sign2,2,8))/2d0*polcolgq
         res(3,8) = diag(3) * dipole(1) * 
     &      ( pdf(2*sign1,1,8)*pdf(4*sign2,2,8) + pdf(4*sign1,1,8)*pdf(2*sign2,2,8))*polcolgq
         res(4,8) = diag(4) * dipole(2) * 
     &      ( pdf(2*sign1,1,8)*pdf(3*sign2,2,8) + pdf(4*sign1,1,8)*pdf(1*sign2,2,8))*polcolgq
         res(7,8) = diag(7) * dipole(3) * 
     &      ( pdf(2*sign1,1,8)*pdf(-3*sign2,2,8) + pdf(4*sign1,1,8)*pdf(-1*sign2,2,8))*polcolgq
         res(8,8) = diag(8) * dipole(4) * 
     &      ( pdf(2*sign1,1,8)*pdf(-4*sign2,2,8) + pdf(4*sign1,1,8)*pdf(-2*sign2,2,8))*polcolgq
         res(10,8) = diag(10) * dipole(3) * 
     &      ( pdf(2*sign1,1,8)*pdf(-1*sign2,2,8) + pdf(4*sign1,1,8)*pdf(-3*sign2,2,8))*polcolgq
         res(12,8) = diag(12) * dipole(4) * 
     &      ( pdf(2*sign1,1,8)*pdf(-2*sign2,2,8) + pdf(4*sign1,1,8)*pdf(-4*sign2,2,8))*polcolgq
         res(16,8) = diag(16) * dipole(5) * 
     &      ( pdf(2*sign1,1,8)*pdf(-2*sign2,2,8) + pdf(4*sign1,1,8)*pdf(-4*sign2,2,8))/2d0*polcolqq
         res(17,8) = diag(17) * dipole(1) * 
     &      ( pdf(1*sign1,1,8)*pdf(2*sign2,2,8) + pdf(3*sign1,1,8)*pdf(4*sign2,2,8))*polcolgq
         res(18,8) = diag(18) * dipole(2) * 
     &      ( pdf(1*sign1,1,8)*pdf(1*sign2,2,8) + pdf(3*sign1,1,8)*pdf(3*sign2,2,8))/2d0*polcolgq
         res(20,8) = diag(20) * dipole(1) * 
     &      ( pdf(1*sign1,1,8)*pdf(4*sign2,2,8) + pdf(3*sign1,1,8)*pdf(2*sign2,2,8))*polcolgq
         res(21,8) = diag(21) * dipole(2) * 
     &      ( pdf(1*sign1,1,8)*pdf(3*sign2,2,8) + pdf(3*sign1,1,8)*pdf(1*sign2,2,8))*polcolgq
         res(24,8) = diag(24) * dipole(3) * 
     &      ( pdf(1*sign1,1,8)*pdf(-3*sign2,2,8) + pdf(3*sign1,1,8)*pdf(-1*sign2,2,8))*polcolgq
         res(25,8) = diag(25) * dipole(4) * 
     &      ( pdf(1*sign1,1,8)*pdf(-4*sign2,2,8) + pdf(3*sign1,1,8)*pdf(-2*sign2,2,8))*polcolgq
         res(27,8) = diag(27) * dipole(3) * 
     &      ( pdf(1*sign1,1,8)*pdf(-1*sign2,2,8) + pdf(3*sign1,1,8)*pdf(-3*sign2,2,8))*polcolgq
         res(30,8) = diag(30) * dipole(6) * 
     &      ( pdf(1*sign1,1,8)*pdf(-1*sign2,2,8) + pdf(3*sign1,1,8)*pdf(-3*sign2,2,8))/2d0*polcolqq
         res(31,8) = diag(31) * dipole(4) * 
     &      ( pdf(1*sign1,1,8)*pdf(-2*sign2,2,8) + pdf(3*sign1,1,8)*pdf(-4*sign2,2,8))*polcolgq
         res(33,8) = diag(33) * dipole(7) * 
     &      ( pdf(0*sign1,1,8)*pdf(2*sign2,2,8) + pdf(0*sign1,1,8)*pdf(4*sign2,2,8))*polcolqq
         res(34,8) = diag(34) * dipole(8) * 
     &      ( pdf(0*sign1,1,8)*pdf(1*sign2,2,8) + pdf(0*sign1,1,8)*pdf(3*sign2,2,8))*polcolqq
         res(35,8) = diag(35) * dipole(9) * 
     &      ( pdf(0*sign1,1,8)*pdf(0*sign2,2,8) + pdf(0*sign1,1,8)*pdf(0*sign2,2,8))*polcolgq
         res(36,8) = diag(36) * dipole(10) * 
     &      ( pdf(0*sign1,1,8)*pdf(0*sign2,2,8) + pdf(0*sign1,1,8)*pdf(0*sign2,2,8))*polcolgq
         res(37,8) = diag(37) * dipole(11) * 
     &      ( pdf(0*sign1,1,8)*pdf(-1*sign2,2,8) + pdf(0*sign1,1,8)*pdf(-3*sign2,2,8))*polcolqq
         res(38,8) = diag(38) * dipole(12) * 
     &      ( pdf(0*sign1,1,8)*pdf(-2*sign2,2,8) + pdf(0*sign1,1,8)*pdf(-4*sign2,2,8))*polcolqq
         res(48,8) = diag(48) * dipole(13) * 
     &      ( pdf(-1*sign1,1,8)*pdf(1*sign2,2,8) + pdf(-3*sign1,1,8)*pdf(3*sign2,2,8))/2d0*polcolqq
         res(50,8) = diag(50) * dipole(3) * 
     &      ( pdf(-1*sign1,1,8)*pdf(-1*sign2,2,8) + pdf(-3*sign1,1,8)*pdf(-3*sign2,2,8))/2d0*polcolgq
         res(62,8) = diag(62) * dipole(14) * 
     &      ( pdf(-2*sign1,1,8)*pdf(2*sign2,2,8) + pdf(-4*sign1,1,8)*pdf(4*sign2,2,8))/2d0*polcolqq
         res(66,8) = diag(66) * dipole(3) * 
     &      ( pdf(-2*sign1,1,8)*pdf(-1*sign2,2,8) + pdf(-4*sign1,1,8)*pdf(-3*sign2,2,8))*polcolgq
         res(67,8) = diag(67) * dipole(4) * 
     &      ( pdf(-2*sign1,1,8)*pdf(-2*sign2,2,8) + pdf(-4*sign1,1,8)*pdf(-4*sign2,2,8))/2d0*polcolgq
       endif 

       if (lok(9)) then 
         call qqgWW_curr(borncurr(0,-1,1,40),9,p(0,3,9),p(0,1,9),p(0,2,9),-1,-1,-1)
         call qqgWW_curr(borncurr(0,-1,1,41),9,p(0,2,9),p(0,1,9),p(0,3,9),1,-1,1)
         call qqgWW_curr(borncurr(0,-1,1,42),9,p(0,3,9),p(0,2,9),p(0,1,9),-1,-1,-1)
         call qqgWW_curr(borncurr(0,-1,1,43),9,p(0,1,9),p(0,2,9),p(0,3,9),1,-1,1)
         call qqgWW_curr(borncurr(0,-1,1,44),9,p(0,1,9),p(0,3,9),p(0,2,9),1,1,-1)

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
     &      ( pdf(2*sign1,1,9)*pdf(2*sign2,2,9) + pdf(4*sign1,1,9)*pdf(4*sign2,2,9))/2d0*polcolgq
         res(2,9) = diag(2) * dipole(1) * 
     &      ( pdf(2*sign1,1,9)*pdf(1*sign2,2,9) + pdf(4*sign1,1,9)*pdf(3*sign2,2,9))*polcolgq
         res(6,9) = diag(6) * dipole(2) * 
     &      ( pdf(2*sign1,1,9)*pdf(0*sign2,2,9) + pdf(4*sign1,1,9)*pdf(0*sign2,2,9))*polcolqq
         res(16,9) = diag(16) * dipole(3) * 
     &      ( pdf(2*sign1,1,9)*pdf(-2*sign2,2,9) + pdf(4*sign1,1,9)*pdf(-4*sign2,2,9))/2d0*polcolqq
         res(18,9) = diag(18) * dipole(4) * 
     &      ( pdf(1*sign1,1,9)*pdf(1*sign2,2,9) + pdf(3*sign1,1,9)*pdf(3*sign2,2,9))/2d0*polcolgq
         res(22,9) = diag(22) * dipole(5) * 
     &      ( pdf(1*sign1,1,9)*pdf(0*sign2,2,9) + pdf(3*sign1,1,9)*pdf(0*sign2,2,9))*polcolqq
         res(30,9) = diag(30) * dipole(6) * 
     &      ( pdf(1*sign1,1,9)*pdf(-1*sign2,2,9) + pdf(3*sign1,1,9)*pdf(-3*sign2,2,9))/2d0*polcolqq
         res(35,9) = diag(35) * dipole(7) * 
     &      ( pdf(0*sign1,1,9)*pdf(0*sign2,2,9) + pdf(0*sign1,1,9)*pdf(0*sign2,2,9))*polcolgq
         res(36,9) = diag(36) * dipole(8) * 
     &      ( pdf(0*sign1,1,9)*pdf(0*sign2,2,9) + pdf(0*sign1,1,9)*pdf(0*sign2,2,9))*polcolgq
         res(39,9) = diag(39) * dipole(9) * 
     &      ( pdf(-1*sign1,1,9)*pdf(4*sign2,2,9) + pdf(-3*sign1,1,9)*pdf(2*sign2,2,9))*polcolgq
         res(41,9) = diag(41) * dipole(9) * 
     &      ( pdf(-1*sign1,1,9)*pdf(3*sign2,2,9) + pdf(-3*sign1,1,9)*pdf(1*sign2,2,9))*polcolgq
         res(43,9) = diag(43) * dipole(9) * 
     &      ( pdf(-1*sign1,1,9)*pdf(2*sign2,2,9) + pdf(-3*sign1,1,9)*pdf(4*sign2,2,9))*polcolgq
         res(47,9) = diag(47) * dipole(9) * 
     &      ( pdf(-1*sign1,1,9)*pdf(1*sign2,2,9) + pdf(-3*sign1,1,9)*pdf(3*sign2,2,9))*polcolgq
         res(48,9) = diag(48) * dipole(10) * 
     &      ( pdf(-1*sign1,1,9)*pdf(1*sign2,2,9) + pdf(-3*sign1,1,9)*pdf(3*sign2,2,9))/2d0*polcolqq
         res(49,9) = diag(49) * dipole(11) * 
     &      ( pdf(-1*sign1,1,9)*pdf(0*sign2,2,9) + pdf(-3*sign1,1,9)*pdf(0*sign2,2,9))*polcolqq
         res(50,9) = diag(50) * dipole(9) * 
     &      ( pdf(-1*sign1,1,9)*pdf(-1*sign2,2,9) + pdf(-3*sign1,1,9)*pdf(-3*sign2,2,9))/2d0*polcolgq
         res(51,9) = diag(51) * dipole(9) * 
     &      ( pdf(-1*sign1,1,9)*pdf(-2*sign2,2,9) + pdf(-3*sign1,1,9)*pdf(-4*sign2,2,9))*polcolgq
         res(52,9) = diag(52) * dipole(9) * 
     &      ( pdf(-1*sign1,1,9)*pdf(-3*sign2,2,9) + pdf(-3*sign1,1,9)*pdf(-1*sign2,2,9))*polcolgq
         res(53,9) = diag(53) * dipole(9) * 
     &      ( pdf(-1*sign1,1,9)*pdf(-4*sign2,2,9) + pdf(-3*sign1,1,9)*pdf(-2*sign2,2,9))*polcolgq
         res(55,9) = diag(55) * dipole(12) * 
     &      ( pdf(-2*sign1,1,9)*pdf(4*sign2,2,9) + pdf(-4*sign1,1,9)*pdf(2*sign2,2,9))*polcolgq
         res(57,9) = diag(57) * dipole(12) * 
     &      ( pdf(-2*sign1,1,9)*pdf(3*sign2,2,9) + pdf(-4*sign1,1,9)*pdf(1*sign2,2,9))*polcolgq
         res(60,9) = diag(60) * dipole(12) * 
     &      ( pdf(-2*sign1,1,9)*pdf(2*sign2,2,9) + pdf(-4*sign1,1,9)*pdf(4*sign2,2,9))*polcolgq
         res(62,9) = diag(62) * dipole(13) * 
     &      ( pdf(-2*sign1,1,9)*pdf(2*sign2,2,9) + pdf(-4*sign1,1,9)*pdf(4*sign2,2,9))/2d0*polcolqq
         res(64,9) = diag(64) * dipole(12) * 
     &      ( pdf(-2*sign1,1,9)*pdf(1*sign2,2,9) + pdf(-4*sign1,1,9)*pdf(3*sign2,2,9))*polcolgq
         res(65,9) = diag(65) * dipole(14) * 
     &      ( pdf(-2*sign1,1,9)*pdf(0*sign2,2,9) + pdf(-4*sign1,1,9)*pdf(0*sign2,2,9))*polcolqq
         res(67,9) = diag(67) * dipole(12) * 
     &      ( pdf(-2*sign1,1,9)*pdf(-2*sign2,2,9) + pdf(-4*sign1,1,9)*pdf(-4*sign2,2,9))/2d0*polcolgq
         res(69,9) = diag(69) * dipole(12) * 
     &      ( pdf(-2*sign1,1,9)*pdf(-3*sign2,2,9) + pdf(-4*sign1,1,9)*pdf(-1*sign2,2,9))*polcolgq
         res(70,9) = diag(70) * dipole(12) * 
     &      ( pdf(-2*sign1,1,9)*pdf(-4*sign2,2,9) + pdf(-4*sign1,1,9)*pdf(-2*sign2,2,9))*polcolgq
       endif 

       if (lok(10)) then 
         call qqgWW_curr(borncurr(0,-1,1,45),10,p(0,3,10),p(0,1,10),p(0,2,10),-1,-1,-1)
         call qqgWW_curr(borncurr(0,-1,1,46),10,p(0,2,10),p(0,1,10),p(0,3,10),1,-1,1)
         call qqgWW_curr(borncurr(0,-1,1,47),10,p(0,3,10),p(0,2,10),p(0,1,10),-1,-1,-1)
         call qqgWW_curr(borncurr(0,-1,1,48),10,p(0,2,10),p(0,3,10),p(0,1,10),1,1,-1)
         call qqgWW_curr(borncurr(0,-1,1,49),10,p(0,1,10),p(0,2,10),p(0,3,10),1,-1,1)
         call qqgWW_curr(borncurr(0,-1,1,50),10,p(0,1,10),p(0,3,10),p(0,2,10),1,1,-1)

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

         res(6,10) = diag(6) * dipole(1) * 
     &      ( pdf(2*sign1,1,10)*pdf(0*sign2,2,10) + pdf(4*sign1,1,10)*pdf(0*sign2,2,10))*polcolgq
         res(12,10) = diag(12) * dipole(2) * 
     &      ( pdf(2*sign1,1,10)*pdf(-2*sign2,2,10) + pdf(4*sign1,1,10)*pdf(-4*sign2,2,10))*polcolqq
         res(13,10) = diag(13) * dipole(2) * 
     &      ( pdf(2*sign1,1,10)*pdf(-2*sign2,2,10) + pdf(4*sign1,1,10)*pdf(-4*sign2,2,10))*polcolqq
         res(14,10) = diag(14) * dipole(2) * 
     &      ( pdf(2*sign1,1,10)*pdf(-2*sign2,2,10) + pdf(4*sign1,1,10)*pdf(-4*sign2,2,10))*polcolqq
         res(15,10) = diag(15) * dipole(2) * 
     &      ( pdf(2*sign1,1,10)*pdf(-2*sign2,2,10) + pdf(4*sign1,1,10)*pdf(-4*sign2,2,10))*polcolqq
         res(16,10) = diag(16) * dipole(3) * 
     &      ( pdf(2*sign1,1,10)*pdf(-2*sign2,2,10) + pdf(4*sign1,1,10)*pdf(-4*sign2,2,10))/2d0*polcolqq
         res(22,10) = diag(22) * dipole(4) * 
     &      ( pdf(1*sign1,1,10)*pdf(0*sign2,2,10) + pdf(3*sign1,1,10)*pdf(0*sign2,2,10))*polcolgq
         res(26,10) = diag(26) * dipole(5) * 
     &      ( pdf(1*sign1,1,10)*pdf(-1*sign2,2,10) + pdf(3*sign1,1,10)*pdf(-3*sign2,2,10))*polcolqq
         res(27,10) = diag(27) * dipole(5) * 
     &      ( pdf(1*sign1,1,10)*pdf(-1*sign2,2,10) + pdf(3*sign1,1,10)*pdf(-3*sign2,2,10))*polcolqq
         res(28,10) = diag(28) * dipole(5) * 
     &      ( pdf(1*sign1,1,10)*pdf(-1*sign2,2,10) + pdf(3*sign1,1,10)*pdf(-3*sign2,2,10))*polcolqq
         res(29,10) = diag(29) * dipole(5) * 
     &      ( pdf(1*sign1,1,10)*pdf(-1*sign2,2,10) + pdf(3*sign1,1,10)*pdf(-3*sign2,2,10))*polcolqq
         res(30,10) = diag(30) * dipole(6) * 
     &      ( pdf(1*sign1,1,10)*pdf(-1*sign2,2,10) + pdf(3*sign1,1,10)*pdf(-3*sign2,2,10))/2d0*polcolqq
         res(33,10) = diag(33) * dipole(7) * 
     &      ( pdf(0*sign1,1,10)*pdf(2*sign2,2,10) + pdf(0*sign1,1,10)*pdf(4*sign2,2,10))*polcolgq
         res(34,10) = diag(34) * dipole(8) * 
     &      ( pdf(0*sign1,1,10)*pdf(1*sign2,2,10) + pdf(0*sign1,1,10)*pdf(3*sign2,2,10))*polcolgq
         res(37,10) = diag(37) * dipole(9) * 
     &      ( pdf(0*sign1,1,10)*pdf(-1*sign2,2,10) + pdf(0*sign1,1,10)*pdf(-3*sign2,2,10))*polcolgq
         res(38,10) = diag(38) * dipole(10) * 
     &      ( pdf(0*sign1,1,10)*pdf(-2*sign2,2,10) + pdf(0*sign1,1,10)*pdf(-4*sign2,2,10))*polcolgq
         res(44,10) = diag(44) * dipole(11) * 
     &      ( pdf(-1*sign1,1,10)*pdf(1*sign2,2,10) + pdf(-3*sign1,1,10)*pdf(3*sign2,2,10))*polcolqq
         res(45,10) = diag(45) * dipole(11) * 
     &      ( pdf(-1*sign1,1,10)*pdf(1*sign2,2,10) + pdf(-3*sign1,1,10)*pdf(3*sign2,2,10))*polcolqq
         res(46,10) = diag(46) * dipole(11) * 
     &      ( pdf(-1*sign1,1,10)*pdf(1*sign2,2,10) + pdf(-3*sign1,1,10)*pdf(3*sign2,2,10))*polcolqq
         res(47,10) = diag(47) * dipole(11) * 
     &      ( pdf(-1*sign1,1,10)*pdf(1*sign2,2,10) + pdf(-3*sign1,1,10)*pdf(3*sign2,2,10))*polcolqq
         res(48,10) = diag(48) * dipole(12) * 
     &      ( pdf(-1*sign1,1,10)*pdf(1*sign2,2,10) + pdf(-3*sign1,1,10)*pdf(3*sign2,2,10))/2d0*polcolqq
         res(49,10) = diag(49) * dipole(13) * 
     &      ( pdf(-1*sign1,1,10)*pdf(0*sign2,2,10) + pdf(-3*sign1,1,10)*pdf(0*sign2,2,10))*polcolgq
         res(58,10) = diag(58) * dipole(14) * 
     &      ( pdf(-2*sign1,1,10)*pdf(2*sign2,2,10) + pdf(-4*sign1,1,10)*pdf(4*sign2,2,10))*polcolqq
         res(59,10) = diag(59) * dipole(14) * 
     &      ( pdf(-2*sign1,1,10)*pdf(2*sign2,2,10) + pdf(-4*sign1,1,10)*pdf(4*sign2,2,10))*polcolqq
         res(60,10) = diag(60) * dipole(14) * 
     &      ( pdf(-2*sign1,1,10)*pdf(2*sign2,2,10) + pdf(-4*sign1,1,10)*pdf(4*sign2,2,10))*polcolqq
         res(61,10) = diag(61) * dipole(14) * 
     &      ( pdf(-2*sign1,1,10)*pdf(2*sign2,2,10) + pdf(-4*sign1,1,10)*pdf(4*sign2,2,10))*polcolqq
         res(62,10) = diag(62) * dipole(15) * 
     &      ( pdf(-2*sign1,1,10)*pdf(2*sign2,2,10) + pdf(-4*sign1,1,10)*pdf(4*sign2,2,10))/2d0*polcolqq
         res(65,10) = diag(65) * dipole(16) * 
     &      ( pdf(-2*sign1,1,10)*pdf(0*sign2,2,10) + pdf(-4*sign1,1,10)*pdf(0*sign2,2,10))*polcolgq
       endif 

       if (lok(11)) then 
         call qqgWW_curr(borncurr(0,-1,1,51),11,p(0,3,11),p(0,1,11),p(0,2,11),-1,-1,-1)
         call qqgWW_curr(borncurr(0,-1,1,52),11,p(0,2,11),p(0,1,11),p(0,3,11),1,-1,1)
         call qqgWW_curr(borncurr(0,-1,1,53),11,p(0,3,11),p(0,2,11),p(0,1,11),-1,-1,-1)
         call qqgWW_curr(borncurr(0,-1,1,54),11,p(0,2,11),p(0,3,11),p(0,1,11),1,1,-1)
         call qqgWW_curr(borncurr(0,-1,1,55),11,p(0,1,11),p(0,2,11),p(0,3,11),1,-1,1)
         call qqgWW_curr(borncurr(0,-1,1,56),11,p(0,1,11),p(0,3,11),p(0,2,11),1,1,-1)

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

         res(6,11) = diag(6) * dipole(1) * 
     &      ( pdf(2*sign1,1,11)*pdf(0*sign2,2,11) + pdf(4*sign1,1,11)*pdf(0*sign2,2,11))*polcolgq
         res(12,11) = diag(12) * dipole(2) * 
     &      ( pdf(2*sign1,1,11)*pdf(-2*sign2,2,11) + pdf(4*sign1,1,11)*pdf(-4*sign2,2,11))*polcolqq
         res(13,11) = diag(13) * dipole(2) * 
     &      ( pdf(2*sign1,1,11)*pdf(-2*sign2,2,11) + pdf(4*sign1,1,11)*pdf(-4*sign2,2,11))*polcolqq
         res(14,11) = diag(14) * dipole(2) * 
     &      ( pdf(2*sign1,1,11)*pdf(-2*sign2,2,11) + pdf(4*sign1,1,11)*pdf(-4*sign2,2,11))*polcolqq
         res(15,11) = diag(15) * dipole(2) * 
     &      ( pdf(2*sign1,1,11)*pdf(-2*sign2,2,11) + pdf(4*sign1,1,11)*pdf(-4*sign2,2,11))*polcolqq
         res(16,11) = diag(16) * dipole(3) * 
     &      ( pdf(2*sign1,1,11)*pdf(-2*sign2,2,11) + pdf(4*sign1,1,11)*pdf(-4*sign2,2,11))/2d0*polcolqq
         res(22,11) = diag(22) * dipole(4) * 
     &      ( pdf(1*sign1,1,11)*pdf(0*sign2,2,11) + pdf(3*sign1,1,11)*pdf(0*sign2,2,11))*polcolgq
         res(26,11) = diag(26) * dipole(5) * 
     &      ( pdf(1*sign1,1,11)*pdf(-1*sign2,2,11) + pdf(3*sign1,1,11)*pdf(-3*sign2,2,11))*polcolqq
         res(27,11) = diag(27) * dipole(5) * 
     &      ( pdf(1*sign1,1,11)*pdf(-1*sign2,2,11) + pdf(3*sign1,1,11)*pdf(-3*sign2,2,11))*polcolqq
         res(28,11) = diag(28) * dipole(5) * 
     &      ( pdf(1*sign1,1,11)*pdf(-1*sign2,2,11) + pdf(3*sign1,1,11)*pdf(-3*sign2,2,11))*polcolqq
         res(29,11) = diag(29) * dipole(5) * 
     &      ( pdf(1*sign1,1,11)*pdf(-1*sign2,2,11) + pdf(3*sign1,1,11)*pdf(-3*sign2,2,11))*polcolqq
         res(30,11) = diag(30) * dipole(6) * 
     &      ( pdf(1*sign1,1,11)*pdf(-1*sign2,2,11) + pdf(3*sign1,1,11)*pdf(-3*sign2,2,11))/2d0*polcolqq
         res(33,11) = diag(33) * dipole(7) * 
     &      ( pdf(0*sign1,1,11)*pdf(2*sign2,2,11) + pdf(0*sign1,1,11)*pdf(4*sign2,2,11))*polcolgq
         res(34,11) = diag(34) * dipole(8) * 
     &      ( pdf(0*sign1,1,11)*pdf(1*sign2,2,11) + pdf(0*sign1,1,11)*pdf(3*sign2,2,11))*polcolgq
         res(37,11) = diag(37) * dipole(9) * 
     &      ( pdf(0*sign1,1,11)*pdf(-1*sign2,2,11) + pdf(0*sign1,1,11)*pdf(-3*sign2,2,11))*polcolgq
         res(38,11) = diag(38) * dipole(10) * 
     &      ( pdf(0*sign1,1,11)*pdf(-2*sign2,2,11) + pdf(0*sign1,1,11)*pdf(-4*sign2,2,11))*polcolgq
         res(44,11) = diag(44) * dipole(11) * 
     &      ( pdf(-1*sign1,1,11)*pdf(1*sign2,2,11) + pdf(-3*sign1,1,11)*pdf(3*sign2,2,11))*polcolqq
         res(45,11) = diag(45) * dipole(11) * 
     &      ( pdf(-1*sign1,1,11)*pdf(1*sign2,2,11) + pdf(-3*sign1,1,11)*pdf(3*sign2,2,11))*polcolqq
         res(46,11) = diag(46) * dipole(11) * 
     &      ( pdf(-1*sign1,1,11)*pdf(1*sign2,2,11) + pdf(-3*sign1,1,11)*pdf(3*sign2,2,11))*polcolqq
         res(47,11) = diag(47) * dipole(11) * 
     &      ( pdf(-1*sign1,1,11)*pdf(1*sign2,2,11) + pdf(-3*sign1,1,11)*pdf(3*sign2,2,11))*polcolqq
         res(48,11) = diag(48) * dipole(12) * 
     &      ( pdf(-1*sign1,1,11)*pdf(1*sign2,2,11) + pdf(-3*sign1,1,11)*pdf(3*sign2,2,11))/2d0*polcolqq
         res(49,11) = diag(49) * dipole(13) * 
     &      ( pdf(-1*sign1,1,11)*pdf(0*sign2,2,11) + pdf(-3*sign1,1,11)*pdf(0*sign2,2,11))*polcolgq
         res(58,11) = diag(58) * dipole(14) * 
     &      ( pdf(-2*sign1,1,11)*pdf(2*sign2,2,11) + pdf(-4*sign1,1,11)*pdf(4*sign2,2,11))*polcolqq
         res(59,11) = diag(59) * dipole(14) * 
     &      ( pdf(-2*sign1,1,11)*pdf(2*sign2,2,11) + pdf(-4*sign1,1,11)*pdf(4*sign2,2,11))*polcolqq
         res(60,11) = diag(60) * dipole(14) * 
     &      ( pdf(-2*sign1,1,11)*pdf(2*sign2,2,11) + pdf(-4*sign1,1,11)*pdf(4*sign2,2,11))*polcolqq
         res(61,11) = diag(61) * dipole(14) * 
     &      ( pdf(-2*sign1,1,11)*pdf(2*sign2,2,11) + pdf(-4*sign1,1,11)*pdf(4*sign2,2,11))*polcolqq
         res(62,11) = diag(62) * dipole(15) * 
     &      ( pdf(-2*sign1,1,11)*pdf(2*sign2,2,11) + pdf(-4*sign1,1,11)*pdf(4*sign2,2,11))/2d0*polcolqq
         res(65,11) = diag(65) * dipole(16) * 
     &      ( pdf(-2*sign1,1,11)*pdf(0*sign2,2,11) + pdf(-4*sign1,1,11)*pdf(0*sign2,2,11))*polcolgq
       endif 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c finite collinears
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c put everything in entry 71, split equally between 4/5 and 8/9
c split into individual contributions to be done later

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

         res(71,4) = diag(71) * fincolfac * (
c g u/c: 16 1
     &     + borncurr(0,0,1,16) * fincol(1) *
     &       (pdf(2*sign2,2,4) + pdf(4*sign2,2,4))*polcolgq
c g d/s: 16 1
     &     + borncurr(0,0,2,16) * fincol(1) *
     &       (pdf(1*sign2,2,4) + pdf(3*sign2,2,4))*polcolgq
c g U/C: 15 1
     &     + borncurr(0,0,1,15) * fincol(1) *
     &       (pdf(-2*sign2,2,4) + pdf(-4*sign2,2,4))*polcolgq
c g D/S: 15 1
     &     + borncurr(0,0,2,15) * fincol(1) *
     &       (pdf(-1*sign2,2,4) + pdf(-3*sign2,2,4))*polcolgq
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
c D g: 13 12
     &     + borncurr(0,0,2,13) * fincol(12) *
     &       pdf(0*sign2,2,4)*polcolgq
c S g: 13 13
     &     + borncurr(0,0,2,13) * fincol(13) *
     &       pdf(0*sign2,2,4)*polcolgq
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
c D d: 14 10
     &     + borncurr(0,0,2,14) * fincol(10) *
     &       pdf(1*sign2,2,4)*polcolqq
c S s: 14 11
     &     + borncurr(0,0,2,14) * fincol(11) *
     &       pdf(3*sign2,2,4)*polcolqq
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

         res(71,5) = diag(71) * fincolfac * (
c u/c g: 24 1
     &     + borncurr(0,0,1,24) * fincol(1) *
     &       (pdf(2*sign1,1,5) + pdf(4*sign1,1,5))*polcolgq
c d/s g: 24 1
     &     + borncurr(0,0,2,24) * fincol(1) *
     &       (pdf(1*sign1,1,5) + pdf(3*sign1,1,5))*polcolgq
c U/C g: 19 1
     &     + borncurr(0,0,1,19) * fincol(1) *
     &       (pdf(-2*sign1,1,5) + pdf(-4*sign1,1,5))*polcolgq
c D/S g: 19 1
     &     + borncurr(0,0,2,19) * fincol(1) *
     &       (pdf(-1*sign1,1,5) + pdf(-3*sign1,1,5))*polcolgq
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
c g D: 21 12
     &     + borncurr(0,0,2,21) * fincol(12) *
     &       pdf(0*sign1,1,5)*polcolgq
c g S: 21 13
     &     + borncurr(0,0,2,21) * fincol(13) *
     &       pdf(0*sign1,1,5)*polcolgq
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
c d D: 23 10
     &     + borncurr(0,0,2,23) * fincol(10) *
     &       pdf(1*sign1,1,5)*polcolqq
c s S: 23 11
     &     + borncurr(0,0,2,23) * fincol(11) *
     &       pdf(3*sign1,1,5)*polcolqq
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
         call qqgWW_curr(borncurr(0,-1,1,61),8,p(0,1,8),p(0,3,8),p(0,2,8),1,1,-1)

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

         res(71,8) = diag(71) * fincolfac * (
c g u/c: 38 1
     &     + borncurr(0,0,1,38) * fincol(1) *
     &       (pdf(2*sign2,2,8) + pdf(4*sign2,2,8))*polcolgq
c g d/s: 38 1
     &     + borncurr(0,0,2,38) * fincol(1) *
     &       (pdf(1*sign2,2,8) + pdf(3*sign2,2,8))*polcolgq
c g U/C: 37 1
     &     + borncurr(0,0,1,37) * fincol(1) *
     &       (pdf(-2*sign2,2,8) + pdf(-4*sign2,2,8))*polcolgq
c g D/S: 37 1
     &     + borncurr(0,0,2,37) * fincol(1) *
     &       (pdf(-1*sign2,2,8) + pdf(-3*sign2,2,8))*polcolgq
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
c D g: 35 12
     &     + borncurr(0,0,2,35) * fincol(12) *
     &       pdf(0*sign2,2,8)*polcolgq
c S g: 35 13
     &     + borncurr(0,0,2,35) * fincol(13) *
     &       pdf(0*sign2,2,8)*polcolgq
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
c D d: 36 10
     &     + borncurr(0,0,2,36) * fincol(10) *
     &       pdf(1*sign2,2,8)*polcolqq
c S s: 36 11
     &     + borncurr(0,0,2,36) * fincol(11) *
     &       pdf(3*sign2,2,8)*polcolqq
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
         call qqgWW_curr(borncurr(0,-1,1,62),9,p(0,2,9),p(0,3,9),p(0,1,9),1,1,-1)

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

         res(71,9) = diag(71) * fincolfac * (
c u/c g: 44 1
     &     + borncurr(0,0,1,44) * fincol(1) *
     &       (pdf(2*sign1,1,9) + pdf(4*sign1,1,9))*polcolgq
c d/s g: 44 1
     &     + borncurr(0,0,2,44) * fincol(1) *
     &       (pdf(1*sign1,1,9) + pdf(3*sign1,1,9))*polcolgq
c U/C g: 40 1
     &     + borncurr(0,0,1,40) * fincol(1) *
     &       (pdf(-2*sign1,1,9) + pdf(-4*sign1,1,9))*polcolgq
c D/S g: 40 1
     &     + borncurr(0,0,2,40) * fincol(1) *
     &       (pdf(-1*sign1,1,9) + pdf(-3*sign1,1,9))*polcolgq
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
c g D: 42 12
     &     + borncurr(0,0,2,42) * fincol(12) *
     &       pdf(0*sign1,1,9)*polcolgq
c g S: 42 13
     &     + borncurr(0,0,2,42) * fincol(13) *
     &       pdf(0*sign1,1,9)*polcolgq
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
c d D: 43 10
     &     + borncurr(0,0,2,43) * fincol(10) *
     &       pdf(1*sign1,1,9)*polcolqq
c s S: 43 11
     &     + borncurr(0,0,2,43) * fincol(11) *
     &       pdf(3*sign1,1,9)*polcolqq
c u U: 43 14
     &     + borncurr(0,0,1,43) * fincol(14) *
     &       pdf(2*sign1,1,9)*polcolqq
c c C: 43 15
     &     + borncurr(0,0,1,43) * fincol(15) *
     &       pdf(4*sign1,1,9)*polcolqq
     &   )
       endif 

c*****************  end of process evaluation  *********************

c multiply results with hadronic factor
      SELECT CASE(procid)
      CASE(WPhadWMj)
       call hadDecayFactor_W(finalquarks(1), abs(mjj2(v(0,1,1),v(0,2,1))), N_gen_W, fac_W)
       do j=0,n_kin
        DO IRES = 1,NMAX
          res(IRES,j) = res(IRES,j) * fac_W            ! factor fac_W for hadronic decay
        ENDDO
       enddo
      CASE(WPWMhadj)
       call hadDecayFactor_W(finalquarks(1), abs(mjj2(v(0,3,1),v(0,4,1))), N_gen_W, fac_W)
       do j=0,n_kin
        DO IRES = 1,NMAX
          res(IRES,j) = res(IRES,j) * fac_W
        ENDDO
       enddo
      END SELECT

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


c*****************************************************************************
      double precision function dipoleVVj_CS39(kin,col,pmi,pmj,pma,xuz,borncurr)
c*****************************************************************************
      implicit none
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/scales.inc"

c INPUT
      integer kin, col
      double precision pmi(0:3), pmj(0:3), pma(0:3), xuz(2)
      double complex borncurr
c EXTERNAL & CONSTANTS
      double precision N, TR, CF, CA
      parameter(N=3d0, TR=1/2d0, CA=N, CF=(N**2-1)/(2d0*N))
      double precision dotrr
      external dotrr
c LOCAL

      double precision dipole, color
      double precision pmipmj, xija, zi

      pmipmj = dotrr(pmi,pmj)
      xija = xuz(1)
      zi   = xuz(2)

      dipole = -1/(2d0*pmipmj*xija)*
     &            8d0*pi*als(1,kin)*CF*
     &           (2/(1-zi+(1-xija))-(1+zi))
      if (col .eq. 0) then
        color = -CA/2d0/CF
      else
        color = (CA-2*CF)/2d0/CF
      endif
      
      dipoleVVj_CS39 = dipole*borncurr*color

      end
c*****************************************************************************

c*****************************************************************************
      double precision function dipoleVVj_CS40(kin,col,pmi,pmj,pma,xuz,borncurr)
c*****************************************************************************
      implicit none
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/scales.inc"

c INPUT
      integer kin, col
      double precision pmi(0:3), pmj(0:3), pma(0:3), xuz(2)
      double complex borncurr(0:3,-1:1)
c EXTERNAL & CONSTANTS
      double precision N, TR, CF, CA
      parameter(N=3d0, TR=1/2d0, CA=N, CF=(N**2-1)/(2d0*N))
      double precision dotrr, contract_Trjcj
      external dotrr, contract_Trjcj
c LOCAL
      integer isig, mu, nu
      double precision dipole(0:3,0:3), color
      double precision pmipmj, xija, zi, zj

      pmipmj = dotrr(pmi,pmj)
      xija = xuz(1)
      zi   = xuz(2)
      zj   = 1 - zi

      do mu=0,3
      do nu=0,mu
        dipole(mu,nu) = 
     &    -1/(2d0*pmipmj*xija)*16d0*pi*als(1,kin)*CA*
     &    ( 1/pmipmj*
     &      (zi*pmi(mu)-zj*pmj(mu))*(zi*pmi(nu)-zj*pmj(nu)))
        dipole(nu,mu) = dipole(mu,nu)
        if (mu.eq.nu) then
          dipole(mu,nu) = dipole(mu,nu) + sign(1d0,0.5d0-mu)*
     &      1/(2d0*pmipmj*xija)*16d0*pi*als(1,kin)*CA*
     &      (1/(1-zi+(1-xija)) + 1/(1-zj+(1-xija)) - 2)
        endif
      enddo
      enddo

      if (col .eq. 0) then
        write(*,*) "Error in dipole 40: 2 gluons in born ME"
        stop
      else
        color = -1/2d0 *4d0 ! -CA/2d0/CA * borncolour
      endif
      
      dipoleVVj_CS40 = 0d0
      do isig=-1,1,2
        dipoleVVj_CS40 = dipoleVVj_CS40 
     &    + contract_Trjcj(dipole,borncurr(0,isig))*color
      enddo

      end
c*****************************************************************************

c*****************************************************************************
      double precision function dipoleVVj_CS41(kin,col,pmi,pmj,pma,xuz,borncurr)
c*****************************************************************************
      implicit none
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/scales.inc"

c INPUT
      integer kin, col
      double precision pmi(0:3), pmj(0:3), pma(0:3), xuz(2)
      double complex borncurr(0:3,-1:1)
c EXTERNAL & CONSTANTS
      double precision N, TR, CF, CA
      parameter(N=3d0, TR=1/2d0, CA=N, CF=(N**2-1)/(2d0*N))
      double precision dotrr, contract_Trjcj
      external dotrr, contract_Trjcj
c LOCAL
      integer isig, mu, nu
      double precision dipole(0:3,0:3), color
      double precision pmipmj, xija, zi, zj

      pmipmj = dotrr(pmi,pmj)
      xija = xuz(1)
      zi   = xuz(2)
      zj   = 1 - zi

      do mu=0,3
      do nu=0,mu
        dipole(mu,nu) = 
     &    -1/(2d0*pmipmj*xija)*8d0*pi*als(1,kin)*TR*
     &    ( -2/pmipmj*
     &      (zi*pmi(mu)-zj*pmj(mu))*(zi*pmi(nu)-zj*pmj(nu)))
        dipole(nu,mu) = dipole(mu,nu)
        if (mu.eq.nu) then
          dipole(mu,nu) = dipole(mu,nu) + sign(1d0,0.5d0-mu)*
     &      1/(2d0*pmipmj*xija)*8d0*pi*als(1,kin)*TR
        endif
      enddo
      enddo

      if (col .eq. 0) then
        write(*,*) "Error in dipole 41: 2 gluons in born ME"
        stop
      else
        color = -1/2d0 *4d0 ! -CA/2d0/CA * borncolour
      endif
      
      dipoleVVj_CS41 = 0d0
      do isig=-1,1,2
        dipoleVVj_CS41 = dipoleVVj_CS41 
     &    + contract_Trjcj(dipole,borncurr(0,isig))*color
      enddo

      end
c*****************************************************************************

c*****************************************************************************
      double precision function dipoleVVj_CS65(kin,col,pma,pmi,pmk,xuz,borncurr)
c*****************************************************************************
      implicit none
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/scales.inc"

c INPUT
      integer kin, col
      double precision pma(0:3), pmi(0:3), pmk(0:3), xuz(2)
      double complex borncurr
c EXTERNAL & CONSTANTS
      double precision N, TR, CF, CA
      parameter(N=3d0, TR=1/2d0, CA=N, CF=(N**2-1)/(2d0*N))
      double precision dotrr
      external dotrr
c LOCAL

      double precision dipole, color
      double precision pmapmi, xika, ui

      pmapmi = dotrr(pma,pmi)
      xika = xuz(1)
      ui   = xuz(2)

      dipole = -1/(2d0*pmapmi*xika)*
     &            8d0*pi*als(1,kin)*CF*
     &           (2/(1-xika+ui)-(1+xika))
      if (col .eq. 0) then
        color = -CA/2d0/CF
      else
        color = (CA-2*CF)/2d0/CF
      endif
      
      dipoleVVj_CS65 = dipole*borncurr*color

      end
c*****************************************************************************

c*****************************************************************************
      double precision function dipoleVVj_CS66(kin,col,pma,pmi,pmk,xuz,borncurr)
c*****************************************************************************
      implicit none
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/scales.inc"

c INPUT
      integer kin, col
      double precision pma(0:3), pmi(0:3), pmk(0:3), xuz(2)
      double complex borncurr
c EXTERNAL & CONSTANTS
      double precision N, TR, CF, CA
      parameter(N=3d0, TR=1/2d0, CA=N, CF=(N**2-1)/(2d0*N))
      double precision dotrr
      external dotrr
c LOCAL

      double precision dipole, color
      double precision pmapmi, xika, ui

      pmapmi = dotrr(pma,pmi)
      xika = xuz(1)
      ui   = xuz(2)

      dipole = -1/(2d0*pmapmi*xika)*
     &            8d0*pi*als(1,kin)*TR*
     &           (1-2*xika*(1-xika))
      if (col .eq. 0) then
        color = -CA/2d0/CF
      else
        color = (CA-2*CF)/2d0/CF
      endif
      
      dipoleVVj_CS66 = dipole*borncurr*color

      end
c*****************************************************************************

c*****************************************************************************
      double precision function dipoleVVj_CS67(kin,col,pma,pmi,pmk,xuz,borncurr)
c*****************************************************************************
      implicit none
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/scales.inc"

c INPUT
      integer kin, col
      double precision pma(0:3), pmi(0:3), pmk(0:3), xuz(2)
      double complex borncurr(0:3,-1:1)
c EXTERNAL & CONSTANTS
      double precision N, TR, CF, CA
      parameter(N=3d0, TR=1/2d0, CA=N, CF=(N**2-1)/(2d0*N))
      double precision dotrr, contract_Trjcj
      external dotrr, contract_Trjcj
c LOCAL
      integer isig, mu, nu
      double precision dipole(0:3,0:3), color
      double precision pmapmi, pmipmk, xika, ui

      pmapmi = dotrr(pma,pmi)
      pmipmk = dotrr(pmi,pmk)
      xika = xuz(1)
      ui   = xuz(2)

      do mu=0,3
      do nu=0,mu
        dipole(mu,nu) = 
     &    -1/(2d0*pmapmi*xika)*8d0*pi*als(1,kin)*CF*
     &    ( (1d0-xika)/xika*2d0*ui*(1-ui)/pmipmk*
     &      (pmi(mu)/ui-pmk(mu)/(1-ui))*(pmi(nu)/ui-pmk(nu)/(1-ui)))
        dipole(nu,mu) = dipole(mu,nu)
        if (mu.eq.nu) then
          dipole(mu,nu) = dipole(mu,nu) + sign(1d0,0.5d0-mu)*
     &      1/(2d0*pmapmi*xika)*8d0*pi*als(1,kin)*CF*xika
        endif
      enddo
      enddo

      if (col .eq. 0) then
        write(*,*) "Error in dipole 67: 2 gluons in born ME"
        stop
      else
        color = -1/2d0 *4d0 ! -CA/2d0/CA * borncolour
      endif
      
      dipoleVVj_CS67 = 0d0
      do isig=-1,1,2
        dipoleVVj_CS67 = dipoleVVj_CS67 
     &    + contract_Trjcj(dipole,borncurr(0,isig))*color
      enddo

      end
c*****************************************************************************

c*****************************************************************************
      double precision function dipoleVVj_CS68(kin,col,pma,pmi,pmk,xuz,borncurr)
c*****************************************************************************
      implicit none
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/scales.inc"

c INPUT
      integer kin, col
      double precision pma(0:3), pmi(0:3), pmk(0:3), xuz(2)
      double complex borncurr(0:3,-1:1)
c EXTERNAL & CONSTANTS
      double precision N, TR, CF, CA
      parameter(N=3d0, TR=1/2d0, CA=N, CF=(N**2-1)/(2d0*N))
      double precision dotrr, contract_Trjcj
      external dotrr, contract_Trjcj
c LOCAL
      integer isig, mu, nu
      double precision dipole(0:3,0:3), color
      double precision pmapmi, pmipmk, xika, ui

      pmapmi = dotrr(pma,pmi)
      pmipmk = dotrr(pmi,pmk)
      xika = xuz(1)
      ui   = xuz(2)

      do mu=0,3
      do nu=0,mu
        dipole(mu,nu) = 
     &    -1/(2d0*pmapmi*xika)*16d0*pi*als(1,kin)*CA*
     &    ( (1d0-xika)/xika*ui*(1-ui)/pmipmk*
     &      (pmi(mu)/ui-pmk(mu)/(1-ui))*(pmi(nu)/ui-pmk(nu)/(1-ui)))
        dipole(nu,mu) = dipole(mu,nu)
        if (mu.eq.nu) then
          dipole(mu,nu) = dipole(mu,nu) + sign(1d0,0.5d0-mu)*
     &      1/(2d0*pmapmi*xika)*16d0*pi*als(1,kin)*CA*
     &      (1/(1-xika+ui) - 1 + xika*(1-xika))
        endif
      enddo
      enddo

      if (col .eq. 0) then
        write(*,*) "Error in dipole 68: 2 gluons in born ME"
        stop
      else
        color = -1/2d0 *4d0 ! -CA/2d0/CA * borncolour
      endif
      
      dipoleVVj_CS68 = 0d0
      do isig=-1,1,2
        dipoleVVj_CS68 = dipoleVVj_CS68 
     &    + contract_Trjcj(dipole,borncurr(0,isig))*color
      enddo

      end
c*****************************************************************************

c*****************************************************************************
      double precision function dipoleVVj_CS145(kin,col,pma,pmi,pmb,xuz,borncurr)
c*****************************************************************************
      implicit none
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/scales.inc"

c INPUT
      integer kin, col
      double precision pma(0:3), pmi(0:3), pmb(0:3), xuz(2)
      double complex borncurr
c EXTERNAL & CONSTANTS
      double precision N, TR, CF, CA
      parameter(N=3d0, TR=1/2d0, CA=N, CF=(N**2-1)/(2d0*N))
      double precision dotrr
      external dotrr
c LOCAL

      double precision dipole, color
      double precision pmapmb, pmipma, pmipmb, xiab

      pmapmb = dotrr(pma,pmb)
      pmipma = dotrr(pmi,pma)
      pmipmb = dotrr(pmi,pmb)
      xiab = xuz(1)

      dipole = -1/(2d0*pmipma*xiab)*
     &            8d0*pi*als(1,kin)*CF*
     &           (2/(1-xiab)-(1+xiab))
      if (col .eq. 0) then
        color = -CA/2d0/CF
      else
        color = (CA-2*CF)/2d0/CF
      endif
      
      dipoleVVj_CS145 = dipole*borncurr*color

      end
c*****************************************************************************

c*****************************************************************************
      double precision function dipoleVVj_CS146(kin,col,pma,pmi,pmb,xuz,borncurr)
c*****************************************************************************
      implicit none
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/scales.inc"

c INPUT
      integer kin, col
      double precision pma(0:3), pmi(0:3), pmb(0:3), xuz(2)
      double complex borncurr
c EXTERNAL & CONSTANTS
      double precision N, TR, CF, CA
      parameter(N=3d0, TR=1/2d0, CA=N, CF=(N**2-1)/(2d0*N))
      double precision dotrr
      external dotrr
c LOCAL

      double precision dipole, color
      double precision pmapmb, pmipma, pmipmb, xiab

      pmapmb = dotrr(pma,pmb)
      pmipma = dotrr(pmi,pma)
      pmipmb = dotrr(pmi,pmb)
      xiab = xuz(1)

      dipole = -1/(2d0*pmipma*xiab)*
     &            8d0*pi*als(1,kin)*TR*
     &           (1-2*xiab*(1-xiab))
      if (col .eq. 0) then
        color = -CA/2d0/CF
      else
        color = (CA-2*CF)/2d0/CF
      endif
      
      dipoleVVj_CS146 = dipole*borncurr*color

      end
c*****************************************************************************

c*****************************************************************************
      double precision function dipoleVVj_CS147(kin,col,pma,pmi,pmb,xuz,borncurr)
c*****************************************************************************
      implicit none
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/scales.inc"

c INPUT
      integer kin, col
      double precision pma(0:3), pmi(0:3), pmb(0:3), xuz(2)
      double complex borncurr(0:3,-1:1)
c EXTERNAL & CONSTANTS
      double precision N, TR, CF, CA
      parameter(N=3d0, TR=1/2d0, CA=N, CF=(N**2-1)/(2d0*N))
      double precision dotrr, contract_Trjcj
      external dotrr, contract_Trjcj
c LOCAL
      integer isig, mu, nu
      double precision dipole(0:3,0:3), color
      double precision pmapmb, pmipma, pmipmb, xiab

      pmapmb = dotrr(pma,pmb)
      pmipma = dotrr(pmi,pma)
      pmipmb = dotrr(pmi,pmb)
      xiab = xuz(1)

      do mu=0,3
      do nu=0,mu
        dipole(mu,nu) = 
     &    -1/(2d0*pmipma*xiab)*8d0*pi*als(1,kin)*CF*
     &    ( (1d0-xiab)/xiab*2d0*pmapmb/(pmipma*pmipmb)*
     &      (pmi(mu)-pmipma/pmapmb*pmb(mu))*
     &      (pmi(nu)-pmipma/pmapmb*pmb(nu)))
        dipole(nu,mu) = dipole(mu,nu)
        if (mu.eq.nu) then
          dipole(mu,nu) = dipole(mu,nu) + sign(1d0,0.5d0-mu)*
     &      1/(2d0*pmipma*xiab)*8d0*pi*als(1,kin)*CF*xiab
        endif
      enddo
      enddo

      if (col .eq. 0) then
        write(*,*) "Error in dipole 147: 2 gluons in born ME"
        stop
      else
        color = -1/2d0 *4d0 ! -CA/2d0/CA * borncolour
      endif
      
      dipoleVVj_CS147 = 0d0
      do isig=-1,1,2
        dipoleVVj_CS147 = dipoleVVj_CS147 
     &    + contract_Trjcj(dipole,borncurr(0,isig))*color
      enddo

      end
c*****************************************************************************

c*****************************************************************************
      double precision function dipoleVVj_CS148(kin,col,pma,pmi,pmb,xuz,borncurr)
c*****************************************************************************
      implicit none
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/scales.inc"

c INPUT
      integer kin, col
      double precision pma(0:3), pmi(0:3), pmb(0:3), xuz(2)
      double complex borncurr(0:3,-1:1)
c EXTERNAL & CONSTANTS
      double precision N, TR, CF, CA
      parameter(N=3d0, TR=1/2d0, CA=N, CF=(N**2-1)/(2d0*N))
      double precision dotrr, contract_Trjcj
      external dotrr, contract_Trjcj
c LOCAL
      integer isig, mu, nu
      double precision dipole(0:3,0:3), color
      double precision pmapmb, pmipma, pmipmb, xiab

      pmapmb = dotrr(pma,pmb)
      pmipma = dotrr(pmi,pma)
      pmipmb = dotrr(pmi,pmb)
      xiab = xuz(1)

      do mu=0,3
      do nu=0,mu
        dipole(mu,nu) = 
     &    -1/(2d0*pmipma*xiab)*16d0*pi*als(1,kin)*CA*
     &    ( (1d0-xiab)/xiab*pmapmb/(pmipma*pmipmb)*
     &      (pmi(mu)-pmipma/pmapmb*pmb(mu))*
     &      (pmi(nu)-pmipma/pmapmb*pmb(nu)))
        dipole(nu,mu) = dipole(mu,nu)
        if (mu.eq.nu) then
          dipole(mu,nu) = dipole(mu,nu) + sign(1d0,0.5d0-mu)*
     &      1/(2d0*pmipma*xiab)*16d0*pi*als(1,kin)*CA*
     &      (xiab/(1-xiab) + xiab*(1-xiab))
        endif
      enddo
      enddo

      if (col .eq. 0) then
        write(*,*) "Error in dipole 148: 2 gluons in born ME"
        stop
      else
        color = -1/2d0 *4d0 ! -CA/2d0/CA * borncolour
      endif
      
      dipoleVVj_CS148 = 0d0
      do isig=-1,1,2
        dipoleVVj_CS148 = dipoleVVj_CS148 
     &    + contract_Trjcj(dipole,borncurr(0,isig))*color
      enddo

      end
c*****************************************************************************

