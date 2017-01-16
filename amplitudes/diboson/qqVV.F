c****************subroutine QQVV *************************************
c
c	Vera Hankele, <vera@particle.uni-karlsruhe.de>
c	Initial version: April 2006
c	Last modified: February 2009
C
C  QQVV calculates the matrix elements**2 for qq -> WW -> ve e+ mu- vm~
C
C        q1 q2    ---->   W+ W-,   W+ ---> f4-bar f3, W- ---> f6-bar f5
C
C  and crossing related processes. Pauli interference terms for
c  identical fermions are neglected.
c
C  This code is modified to allow for virtual corrections, more precisely
C  the interference of Born with the finite part of virtual diagrams
C  for 
c
c  INPUT:  NLO = 1       return uuww = |M_born|^2 + 2Re(M_Born^* M_virt)
c          NLO = 0       return uuww = |M_born|^2   etc.
c     
c          N = 1         Momentum configuration 1 (LO)
c          N = 2         Momentum configuration 2: emission of 1 line
c          N = 3         Momentum configuration 3: emission of 2 line
c
c*********************************************************************

      subroutine qqVV(pbar, sign, nlo, N, uuww, ddww)

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/coupl.inc"


c     electroweak couplings are taken from KOPPLN
      double precision  clr, xm2, xmg, b, v, a
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     1                  V(4,5),A(4,5)
      double precision  pbar(0:3,2+n_v)
      double precision  uuww, ddww, res(2),resv(2)
      double precision  p(0:3,2+n_v),
     1                  fqp(0:4,2), fqm(0:4,2)
      integer  sign(2+n_v), nlo, mu, i, k,
     1         isig1, isig, N
      integer  ifl(2,2)
      double complex psiwp(2,2), psiwm(2,2),
     1               jwp(0:5,-1:1), jwm(0:5,-1:1)
      double complex temp, temp1
      double complex mat(2,-1:1,3), matv(2,-1:1,3), mborn1, mborn2
      double complex madgr(3,-1:1,2)
      double complex mm(2,-1:1), mv(2,-1:1), zero
      double complex ma, mz, m1, m2

      double complex psi(2,-1:1,4), jqq(0:5,-1:1,2)
      double complex dotcc, s1c
      double precision dotrr
      external dotcc, s1c, dotrr
      save ifl
      parameter (zero = (0d0,0d0) )

c variables for virtual corrections
      double precision scale
      double precision c2,c2o4pi
      parameter (c2=4d0/3d0, c2o4pi=c2/4d0/pi)
      logical lnlo, lbox, ldebugm, linit
      double precision p1(0:3), p2(0:3) 
      double precision minus_qp(0:3),minus_qm(0:3)
      integer wardidtest1, wardidtest2
      double complex gaugebox1(2),gaugebox1_born(2)
      double complex gaugebox2(2), gaugebox2_born(2)
      double complex result(2),result_born(2)
      double complex cte
      data linit /.true./
      save linit
      parameter(ldebugm = .false.)

      lnlo = NLO.ne.0       ! include some virtual stuff if T
      lbox = NLO.eq.-4.or.NLO.eq.1      ! call boxline if T

c define flavors of external quarks
      if (linit) then
         linit = .false.
         k = 1                  ! uuww
         ifl(1,k) = 3
         ifl(2,k) = 3
         k = 2                  ! ddww
         ifl(1,k) = 4
         ifl(2,k) = 4
      endif

      do k = 1,2
         do isig1 = -1,1,2
            do i = 1,3
                  mat(k,isig1,i) = 0
                  matv(k,isig1,i) = 0
            enddo
         enddo
      enddo

c define the internal momenta
      do mu = 0,3
         do i = 1,2+n_v
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo
      enddo

c get the external quark spinors (including factor sqrt(2E) )
      call psi0m(2,pbar(0,1),sign(1),psi)

c------- T-CHANNEL ----------------------------------------------
   
      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,qp(0,N),wp(1,N),
     &              psiwp(1,2),fqp(0,2))
      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,qm(0,N),wm(1,N),
     &              psiwm(1,2),fqm(0,2))

      call curr6(-1,psiwp(1,2),fqp(0,2),
     &              psi(1,-1,1),p(0,1), jwp(0,-1) )
      call curr6(-1,psiwm(1,2),fqm(0,2),
     &              psi(1,-1,1),p(0,1), jwm(0,-1) )

      m1 = dotcc(jwm(0,-1),wp(1,N))   
      m2 = dotcc(jwp(0,-1),wm(1,N))
      mat(1,-1,1) = clr(3,3,-1)**2 * m1
      mat(2,-1,1) = clr(3,3,-1)**2 * m2

c-------  S-CHANNEL --------------------------------------

c get the f-fbar currents J21^mu=jqq(mu,*,1) 

      call curr6(1,psi(1,-1,2),p(0,2),psi(1,-1,1),p(0,1),jqq(0,-1,1))

c contract with atoww/ztoww tensor

      do isig1 = -1,1,2
         ma = -dotcc(aww(0,N),jqq(0,isig1,1))
         mz = -dotcc(zww(0,N),jqq(0,isig1,1))
         do k = 1,2
            mat(k,isig1,2) = ma*clr(ifl(1,k),1,isig1)
            mat(k,isig1,3) = mz*clr(ifl(1,k),2,isig1)
            matv(k,isig1,1) = (0d0,0d0)
         enddo
      enddo

c****************************************************************
c
c      BOX CONTRIBUTIONS
c
c****************************************************************

#ifdef WITH_NLO
      if (lbox) then
         scale = -2d0*dotrr(p(0,1),p(0,2))

c calculate born matrix element and store it in mborn for gauge check:
         mborn1=mat(1,-1,1)+mat(1,-1,2)+mat(1,-1,3)
         mborn2=mat(2,-1,1)+mat(2,-1,2)+mat(2,-1,3)

         do mu = 0,3
            p1(mu) = pbar(mu,1)
            p2(mu) = pbar(mu,2)
            minus_qp(mu) = -qp(mu,1)
            minus_qm(mu) = -p1(mu)-p2(mu)+qp(mu,1)
         enddo

         cte=(3*pi*(0d0,1d0)-4.0d0*pi**2/3.0d0+7.0d0)

         call boxlinemm(p1,minus_qp,minus_qm,p2,
     1    psi(1,-1,2),psi(1,-1,1),wp(1,1), wm(1,1),
     2    -1,scale,3,1,gaugebox1,result(1),gaugebox1_born,
     3    result_born(1))
         call boxlinemm(p1,minus_qm,minus_qp,p2,
     1    psi(1,-1,2),psi(1,-1,1),wm(1,1), wp(1,1),
     2    -1,scale,3,1,gaugebox2,result(2),gaugebox2_born,
     3    result_born(2))

         call Wardtest_Box(wp(1,1), qp(0,1), wm(1,1), qm(0,1),
     1            gaugebox1, gaugebox1_born, mborn1, wardidtest1, als(1,1))

         call Wardtest_Box(wm(1,1), qm(0,1), wp(1,1), qp(0,1),
     1            gaugebox2, gaugebox2_born, mborn2, wardidtest2, als(1,1))

         if (wardidtest1.eq.1) then
                  matv(1,-1,2) = (0.0d0, 0.0d0)
         else
            result(1)=result(1)+cte*result_born(1)
            matv(1,-1,2) = clr(3,3,-1)**2 * result(1) 
         endif

         if (wardidtest2.eq.1) then
                  matv(2,-1,2) = (0.0d0, 0.0d0)
         else
            result(2)=result(2)+cte*result_born(2)
            matv(2,-1,2) = clr(3,3,-1)**2 * result(2) 
         endif

      endif
#endif

#ifdef WITH_MADGRAPH
c------------------------------------------------------
c comparison with madgraph

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

         madgr(1,isig,k) = tamp(1,isig,k)
         madgr(2,isig,k) = tamp(2,isig,k) + tamp(5,isig,k)
     &                   + tamp(7,isig,k)
         madgr(3,isig,k) = tamp(3,isig,k) + tamp(4,isig,k)
     &                   + tamp(6,isig,k) + tamp(8,isig,k)
     &                   + tamp(9,isig,k)

            enddo
         enddo

c print results from checks:

         do isig = -1,1,2
            do k = 1,2
               do i=1,3

    	    print*,"ratio for k = ",k," isig = ",isig," i = ",i,":",
     &	    	mat(k,isig,i)/madgr(i,isig,k)
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

c------------------------------------------------------
#endif

c sum the graphs, square them and map them onto uuww and ddww
      do k = 1,2
         res(k) = 0
         resv(k) = 0
         do isig1 = -1,1,2
            mm(k,isig1) = 0
            do i = 1,3
               mm(k,isig1) = mm(k,isig1) +  mat(k,isig1,i)
            enddo
            res(k)= res(k) + dreal(mm(k,isig1))**2
     &                     + dimag(mm(k,isig1))**2

c  add Born type term and multiply by F_q = alphas*C_2/4pi
c  the additional factor for the born term is after adding the 
c  subtraction term and the counter term for the renormalization
c  of the pdfs (see also D(x) in m2s_qqVVj.f)

#ifdef WITH_NLO
            if (lnlo) then
               mv(k,isig1) = +  matv(k,isig1,2)  !virtual not prop. to born
            endif
#endif
            if (nlo.gt.0) then
#ifdef WITH_NLO
               mv(k,isig1) = als(1,1)*c2o4pi*
     1                  (mv(k,isig1) + mm(k,isig1)* cvirtVV) 
#endif
            else if(nlo.lt.0) then
#ifdef WITH_NLO
               mv(k,isig1) = als(1,1)*c2o4pi*mv(k,isig1)
#endif
            else
                mv(k,isig1) = (0d0, 0d0)
            endif

           resv(k) = resv(k) + 2d0*dreal( mm(k,isig1)   *
     1              conjg( mv(k,isig1) ) )

         enddo

         if (nlo.ge.0) then
            res(k) = (res(k)+resv(k))*3d0      ! 3 is the color sum factor
         else
            res(k) = resv(k)*3d0               ! 3 is the color sum factor
         endif  
      enddo

      uuww = res(1)
      ddww = res(2)

      return
      end


