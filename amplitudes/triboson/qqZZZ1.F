c****************subroutine QQZZZ1 *************************************
c
c	Vera Hankele, <vera@particle.physik.uni-karlsruhe.de>
c	Initial version:  September 2007
c	Last modified: June 2009
C
C  QQZZZ calculates the matrix elements**2 for the Higgs contributions in 
C
C        q1 q2    ---->   Z Z Z,   
C              Z ---> f1-bar f2, Z ---> f3-bar f4, Z ---> f5-bar f6,
C
C         i.e. output is |M12|^2 + |M23|^2 + |M13|^2 where
C         M12 = M(q qbar -> H Z3 -> Z1 Z2 Z3).
c
C  This code is modified to allow for virtual corrections, more precisely
C  the interference of Born with the finite part of virtual diagrams
C  for 
c
c  INPUT:  NLO = 1       return uuww = |M_born|^2 + 2Re(M_Born^* M_virt)
c          NLO = 0       return uuww = |M_born|^2   etc.
c  for steering the calculation of finite box and pentagon contributions
c  the following additional options are implemented
c          NLO = +4      set all finite box and pentagon contributions to 0
c                        but include cvirt*M_Born contribution
c          NLO = +5      true finite pentagons = 0, all other virtual
c                        contributions included, including leading q^mu/m_V
c                        terms of pentagons which can be related to boxes
c       remaining cases:   return uuww = 2Re(M_born^* M_virt). more precisely
c       the piece which is missing for the full NLO=1 when calling +|NLO|     
c          NLO = -4      finite Box + q^mu/m_V terms of Pentagon contributions
c          NLO = -5      true finite Pentagon contributions only
c
c          N = 1         Momentum configuration 1 (LO)
c          N = 2         Momentum configuration 2: emission of parton 1
c          N = 3         Momentum configuration 3: emission of parton 2
c
c*********************************************************************

      subroutine qqZZZ1(pbar,sign, nlo, N, uuzzz, ddzzz, NW)

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/scales.inc"


c     electroweak couplings are taken from KOPPLN
      double precision  clr, xm2, xmg, b, v, a
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     1                  V(4,5),A(4,5)

      double precision  pbar(0:3,2+n_v)
      double precision  uuzzz, ddzzz, res(2),resv(2)
      double precision  p(0:3,2+n_v)

      integer  sign(2+n_v), nlo, mu, i, k, isig, N
      integer  ifl(2)
      double complex temp, temp1
      double complex mat(2,-1:1,3), matv(2,-1:1,3)
      double complex madgr(3,-1:1,2)
      double complex mm(2,-1:1), zero
      double complex mz1, mz2, mz3
      double complex psi(2,-1:1,2), jqq(0:5,-1:1)
      double complex dotcc, s1c, dotrc
      double precision dotrr, psumsq
      external dotcc, s1c, dotrr, dotrc, psumsq
      save ifl
      parameter (zero = (0d0,0d0) )

c variable for q_mu q_nu terms in the propagators of the massive vector bosons

c  helicity selection
      INTEGER h
      COMMON /hcount / h


      double precision c2,c2o4pi
      parameter (c2=4d0/3d0, c2o4pi=c2/4d0/pi)
      logical lnlo, lbox, lpq, lpt, lpent, ldebug, ldebugm, linit
      data linit /.true./
      save linit
      parameter(ldebugm = .false., ldebug=.false.)
      integer NW !NW = 1 means Narrow width approximaion
                 !NW = 0 means full calculation


c---------------------------- BEGIN CODE -----------------------
      lnlo = NLO.ne.0       ! include some virtual stuff if T
      lbox = (NLO.eq.-4) .or. (NLO.eq.1) .or. (NLO.eq.5)    ! call boxline if T
      lpt = (NLO.eq.1) .or. (NLO.eq.-5)   ! true pentagon contribution
      lpq = .false.
      lpent = lpt .or. lpq

c define flavors of external quarks for the 2 subprocesses
      if (linit) then
         linit = .false.
         k = 1                  ! uuwwz: u = 3
         ifl(k) = 3
         k = 2                  ! ddwwz: d = 4
         ifl(k) = 4
      endif

      do k = 1,2
         do isig = -1,1,2
            do i = 1,3
                  mat(k,isig,i) = (0.0d0,0.0d0)
                  matv(k,isig,i) = (0.0d0,0.0d0)
            enddo
         enddo
      enddo

c define the internal momenta
      do mu = 0,3
         do i = 1,2+n_v
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo
      enddo

c------------------------------------------------

c get the external quark spinors (including factor sqrt(2E) )
      call psi0m(2,pbar(0,1),sign(1),psi)

c*******************************************************************
c
c          Tree-level contributions
c
c*******************************************************************

c  Calculation of processes qq -> A/Z -> Z Z Z ---------------------
      call curr6(1,psi(1,-1,2),p(0,2),psi(1,-1,1),p(0,1),jqq(0,-1))

c contract with atowwz/ztowwz tensor
      do isig = -1,1,2
         mz1 = dotcc(zzzz_higgs1(0,N),jqq(0,isig))
         mz2 = dotcc(zzzz_higgs2(0,N),jqq(0,isig))
         mz3 = dotcc(zzzz_higgs3(0,N),jqq(0,isig))
         do k = 1,2
            mat(k,isig,1)  = mz1*clr(ifl(k),2,isig)
            mat(k,isig,2)  = mz2*clr(ifl(k),2,isig)
            mat(k,isig,3)  = mz3*clr(ifl(k),2,isig)
         enddo !k
      enddo !isig

c-----------------------------------------------------
c------------------------------------------------------
c comparison with madgraph

#ifdef WITH_MADGRAPH
      if (ldebugm) then

         do isig = -1,1,2
            do i = 1,3
               do k = 1,2

                  madgr(i,isig,k) = 0d0

               enddo
            enddo
         enddo

         do isig = -1,1,2
            do k = 1,2
               do i =1,3

                  madgr(i,isig,k) = tamp(i,isig,k)

               enddo
            enddo
         enddo

c print results from checks:

         do i=1,3
            do k = 1,2
               do isig = -1,1,2

    	    print*,"ratio for k = ",k," isig = ",isig," i = ",i,":",
     &	    	    mat(k,isig,i)/madgr(i,isig,k)
            print*,"mat(",k,",",isig,",",i,") = ",mat(k,isig,i)
            temp = dreal(mat(k,isig,i))**2+dimag(mat(k,isig,i))**2
            print*,temp
            print*,"madgr(",i,",",isig,",",k,") = ",madgr(i,isig,k)
            temp1 = dreal(madgr(i,isig,k))**2+dimag(madgr(i,isig,k))**2
            print*,temp1
            print*,"temp/temp1 = ",temp/temp1
            print*," "

               enddo
            enddo
         enddo

      endif ! end of debugging block
#endif

c------------------------------------------------------


c sum the graphs, square them and map them onto uuww and ddww

      do k = 1,2
         res(k) = 0
         resv(k) = 0
         do isig = -1,1,2
            mm(k,isig) = 0
            do i = 1,3
               res(k)= res(k) + dreal(mat(k,isig,i))**2
     &                        + dimag(mat(k,isig,i))**2
            enddo
         enddo

c  add Born type term and multiply by F_q = alphas*C_2/4pi
c  the additional factor for the born term is after adding the 
c  subtraction term and the counter term for the renormalization
c  of the pdfs (see also D(x) in m2s_qqWWZj.f)

         if (nlo.ge.0) then
            res(k) = (res(k)+resv(k))*3d0      ! 3 is the color sum factor
         else
            res(k) = resv(k)*3d0               ! 3 is the color sum factor
         endif  
      enddo

      uuzzz = res(1)
      ddzzz = res(2)

      return
      end

c --------------------------------------------------------------------- 
c --------------------------------------------------------------------- 
