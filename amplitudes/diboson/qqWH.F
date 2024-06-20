c****************subroutine QQWH *************************************
c
c  Adapted from qqWA by Johannes Bellm and qqVV by Vera Hankele
c
c  Robin Roth <robin@particle.physik.uni-karlsruhe.de>
c  Initial version: November 2012
c  Last modified: November 2012
C
C  QQWH calculates the matrix elements**2 for qq -> WH -> ve e+ H
C
C        q1 q2    ---->   W H,   W ---> f4-bar f3 
C
C  and crossing related processes. Pauli interference terms for
c  identical fermions are neglected.
c
C  This code is modified to allow for virtual corrections, more precisely
C  the interference of Born with the finite part of virtual diagrams
C  for 
c
c  INPUT:  NLO = 1       return udwa = |M_born|^2 + 2Re(M_Born^* M_virt)
c          NLO = 0       return udwa = |M_born|^2   etc.
c     
c          N = 1         Momentum configuration 1 (LO)
c          N = 2         Momentum configuration 2: emission of 1 line
c          N = 3         Momentum configuration 3: emission of 2 line
c
c*********************************************************************

      subroutine qqWH(pbar, sign, nlo, N, udwh, bos)

      use globalvars, only: ldoblha

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"


c     electroweak couplings are taken from KOPPLN
      double precision  clr, xm2, xmg, b, v, a
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     1                  V(4,5),A(4,5)
      double precision  pbar(0:3,2+n_v)
      double precision  udwh, res(2),resv(2)

      !BLHA interface
      double precision tree(2) ! save tree amplitude
      double precision lnmusqp1p2 ! running logs

      double precision  p(0:3,2+n_v) !n_v=5=2in + 3 out: l v h
      integer  sign(2+n_v), nlo, mu, i, k, isig 
      integer  N, bos
      double complex mat(2,-1:1,3), matv(2,-1:1,3) 
      double complex mm(2,-1:1), mv(2,-1:1), zero
      double complex mwp !, m12,  m22 
      double complex psi(2,-1:1,4), jqq(0:5,-1:1,2)
      double complex dotcc, s1c
      double precision dotrr
      external dotcc, s1c, dotrr
      parameter (zero=(0d0,0d0))

c variables for virtual corrections
      double precision c2,c2o4pi
      parameter (c2=4d0/3d0, c2o4pi=c2/4d0/pi)
      logical lnlo, lbox, ldebug, linit
      data linit /.true./
      save linit
      parameter(ldebug = .false.)

      lnlo = NLO.ne.0       ! include some virtual stuff if T
      lbox = NLO.eq.-4.or.NLO.eq.1      ! call boxline if T

c define flavors of external quarks
      if (linit) then
         linit = .false.
         if (.not. (bos.eq.3 .or. bos.eq.36 .or. bos.eq.4 .or. bos.eq.46)) then
           print *, "Wrong boson ids:", bos
           stop
         endif
      endif

      do k = 1,2
         do isig = -1,1,2
            do i = 1,3
               mat(k,isig,i) = 0
               matv(k,isig,i) = 0
            enddo
         enddo
      enddo

c define the internal momenta
      do mu = 0,3
         do i = 1,2+n_v
            p(mu,i) = pbar(mu,i) !*sign(i) ! signs are fixed before
         enddo
      enddo

      if (ldebug) then
         write(*,*) "in qqWH.F -----------------------------------"
         do i=1,2+n_v
            write(*,*) p(0,i), " ", p(1,i), " ",p(2,i), " ",p(3,i)
         enddo
      endif


c get the external quark spinors (including factor sqrt(2E) )
      call psi0m(2,pbar(0,1),sign(1),psi)

c       c       clr( ferm , bos , hel )


c       C  The fermion indices are:  1: neutrino of electron (muon,tau)
c       C                            2: electron (muon,tau)
c       C                            3: up quark (charm,top)
c       C                            4: down quark (strange,bottom)
c       C  The boson indices are     1: photon
c       C                            2: Z0 boson
c       C                            3: W+ boson
c       C                            4: W- boson
c       C                            5: gluon
c       C                            6: Higgs
c       C  The helicity indices are  -1: lefthanded
c       C                            +1: righthanded


c-------  S-CHANNEL --------------------------------------!

c get the f-fbar currents J21^mu=jqq(mu,*,1) 

c      if(ldebug) print*,"psi ",psi

      call curr6(1,psi(1,-1,2),p(0,2),psi(1,-1,1),p(0,1),jqq(0,-1,1))
c      if(ldebug) print*,"jqq ",jqq

c contract with wtowh tensor


c      if(ldebug) print*,"wtowh ",wtowh
      mwp = -dotcc(wtowh(0,N),jqq(0,-1,1))
c      if(ldebug) print*,"mwp ",mwp

      mat(1,-1,2) = mwp*clr(3,3,-1)
c      if(ldebug) print*,"clr ",clr
      matv(1,-1,1) = (0d0,0d0)

c sum the graphs, square them and map them onto udwa 

      resv(1) = 0
      mm(1,-1) = 0
      mv(1,-1) = (0d0, 0d0)
      res(1)= 0

      mm(1,-1) = mat(1,-1,2)  !+mat(1,-1,1)

      res(1)=   dreal(mm(1,-1))**2 + dimag(mm(1,-1))**2 


c  add Born type term and multiply by F_q = alphas*C_2/4pi
c  the additional factor for the born term is after adding the 
c  subtraction term and the counter term for the renormalization
c  of the pdfs (see also D(x) in m2s_qqVVj.f)

#ifdef WITH_NLO
      if (lnlo) then
         mv(1,-1) = +  matv(1,-1,2)  !virtual not prop. to born
      endif
#endif
      if (nlo.gt.0) then
#ifdef WITH_NLO
         mv(1,-1) = als(1,1)*c2o4pi*
     1                  (mv(1,-1) + mm(1,-1)* cvirtVV) 
#endif
      else if(nlo.lt.0) then
#ifdef WITH_NLO
         mv(1,-1) = als(1,1)*c2o4pi*mv(1,-1)
#endif
      else
         mv(1,-1) = (0d0, 0d0)
      endif

      resv(1) =   2d0*dreal( mm(1,-1)   *
     1              conjg( mv(1,-1) ) )

      if (nlo.ge.0) then
         tree(1) = res(1)*3d0
         res(1) = (res(1)+resv(1))*3d0      ! 3 is the color sum factor
      else
         res(1) = resv(1)*3d0               ! 3 is the color sum factor
      endif  


      udwh = res(1)

      if (ldoblha) then
         select case(blha_amptype(blha_curproc)) 
         case(0) ! tree
            ! blha_amp(1)=blha_amp(1)+
            ! &      res(blha_idsubproc(blha_cursubproc,blha_curproc))
            blha_amp(1)=blha_amp(1) + res(1)

         case(1) ! loop
            lnmusqp1p2 = log(blha_scale**2/(2d0*dotrr(p(0,1),p(0,2)) )) ! sign not used here

            ! eps**-2
            blha_amp(1)=blha_amp(1)+
     &      blha_CF/(2d0*pi)*2d0*
     &       als(1,1)*(-1d0)*
     &       tree(blha_idsubproc(blha_cursubproc,blha_curproc))
            ! eps**-1
            blha_amp(2)=blha_amp(2)+
     &      1d0/(2d0*pi)*2d0*
     &       als(1,1)*(-blha_gammaQuark-blha_CF*lnmusqp1p2)*
     &       tree(blha_idsubproc(blha_cursubproc,blha_curproc))
c this is more complicated than necessary, 
c but should work for other VBF processes as well
            ! eps**0
            ! idsubproc = 1 here
            blha_amp(3)= blha_amp(3)+
     &      resv(1)*blha_CA ! loop + I (no tree)
     &      + (-(als(1,1))/2d0/pi*                        ! -I
     &          (blha_CF*cvirtVV+  ! our cvirt contribution to I operator
     &            2d0*(blha_KQuark+blha_gammaQuark+ ! K, gamma from CS (C.28)
     &                 blha_tgammaQuark+ ! possible conversion term MSbar to DRbar
     &                  (-pi**2/3d0-1)*blha_CF) ) ! remainder to get correct analytic result (pi**2-8)
     &         +als(1,1)/pi*                                       ! +(mu^2/2p1.p2)^epsilon
     &         (-blha_gammaQuark*lnmusqp1p2-blha_CF/2d0*lnmusqp1p2**2)
     &        )*tree(blha_idsubproc(blha_cursubproc,blha_curproc))
            ! tree
            blha_amp(4)= blha_amp(4)+
     &       tree(blha_idsubproc(blha_cursubproc,blha_curproc))

         case(2) ! cctree
            ! 1-2
            call BLHA_cctree(1,2,
     &      -blha_CF*
     &       tree(blha_idsubproc(blha_cursubproc,blha_curproc))
     &    )

         case default 
            call BLHA_amptypeerror(blha_amptype(blha_curproc),__FILE__,__LINE__)
         end select
      endif



      return
      end
