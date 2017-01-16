c****************subroutine QQVVA *************************************
c
c	giuseppe bozzi, <giuseppe@particle.physik.uni-karlsruhe.de>
c	Initial version:  2008, January
c	Last modified: June 2009 by V. Hankele <vera@particle.uni-karlsruhe.de>
C
C  QQVVA calculates the matrix elements**2 for 
C     q1 q2    ---->   W+ W- A   
C     W+ ---> f1-bar f2, W- ---> f3-bar f4
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

      SUBROUTINE qqVVA(pbar,sign, nlo, N, uuwwa, ddwwa)

      IMPLICIT NONE

#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/scales.inc"

c     electroweak couplings are taken from KOPPLN

      DOUBLE PRECISION  clr, xm2, xmg, b, v, a
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     1                  V(4,5),A(4,5)
      DOUBLE PRECISION  pbar(0:3,2+n_v)
      DOUBLE PRECISION  uuwwa, ddwwa, res(2),resv(2)
      DOUBLE PRECISION  p(0:3,2+n_v), pwpz(0:4), pwmz(0:4),
     1                  fqp(0:4,2), fqm(0:4,2),
     2                  fqz(0:4,2)

      INTEGER sign(2+n_v), nlo, mu, i, k, isig, N
      INTEGER ifl(2)
      DOUBLE COMPLEX psiwp(2,-1:1,2), psiwm(2,-1:1,2),
     1               psia(2,-1:1,2), psiz(2,-1:1,2),
     2               jwp(0:5,-1:1,2), jwm(0:5,-1:1,2),
     3               ja(0:5,-1:1,2)


      DOUBLE COMPLEX j5pm(3),j5mp(3), jbpm(3), jbmp(3)

      DOUBLE COMPLEX temp, temp1
      DOUBLE COMPLEX mat(2,-1:1,10), matv(2,-1:1,10), mborn1, mborn2
      DOUBLE COMPLEX madgr(10,-1:1,2)
      DOUBLE COMPLEX mm(2,-1:1), mv(2,-1:1), zero
      DOUBLE COMPLEX ma, mz, m1, m2, m3, m4
      DOUBLE COMPLEX mv1, mv2, mv3, mv4
      DOUBLE COMPLEX psi(2,-1:1,2), jqq(0:5,-1:1)
      DOUBLE COMPLEX dotcc, s1c, dotrc
      DOUBLE PRECISION dotrr, psumsq
      EXTERNAL dotcc, s1c, dotrr, dotrc, psumsq
      SAVE ifl
      PARAMETER (zero = (0d0,0d0) )

c variable for q_mu q_nu terms in the propagators of the massive vector bosons

      DOUBLE COMPLEX zm2i(2:3)
      SAVE zm2i

c  helicity selection

      INTEGER h
      COMMON /hcount / h

c variables for virtual corrections
      DOUBLE PRECISION scale
      double precision c2,c2o4pi
      parameter (c2=4d0/3d0, c2o4pi=c2/4d0/pi)
      logical lnlo, lbox, lpq, lpt, lpent, ldebug, ldebugm, linit
      integer wardidtest1, wardidtest2, wardidtest3

      data linit /.true./
      save linit
      parameter(ldebugm = .false., ldebug=.false.)

c variables for implementation of the pentagon routine
      double complex borntest(3)
      double precision minus_qp(0:3),minus_qm(0:3),minus_qz(0:3)
      double precision minus_qww(0:3), minus_qwmz(0:3), minus_qwpz(0:3)
      double precision p1(0:3), p2(0:3) 

      double complex pentcontr1(3), pentcontr2(3), pentcontr3(3)

c variables for implementation of the box routine
      double complex gaugebox1(2),gaugebox1_born(2), boxtest(2), treetest(2)
      double complex gaugebox2(2), gaugebox2_born(2)
      double complex result(-1:1,2),result_born(-1:1,2)
      double complex cte

c variables for test purposes


c---------------------------- BEGIN CODE -----------------------

      lnlo = NLO.ne.0       ! include some virtual stuff if T
      lbox = (NLO.eq.-4) .or. (NLO.eq.1) .or. (NLO.eq.5)    ! call boxline if T
      lpt = (NLO.eq.1) .or. (NLO.eq.-5)   ! true pentagon contribution
      lpq = (NLO.eq.5) .or. (NLO.eq.-4)   ! q^mu/m_V terms of Pentagon
      lpent = lpt .or. lpq
c define flavors of external quarks for the 2 subprocesses

      if (linit) then
         linit = .false.
         k = 1                  ! uuwwa: u = 3
         ifl(k) = 3
         k = 2                  ! ddwwa: d = 4
         ifl(k) = 4
         zm2i(2) = 1d0/dcmplx(xm2(2),-xmg(2))
         zm2i(3) = 1d0/dcmplx(xm2(3),-xmg(3))
      endif

      do k = 1,2
         do isig = -1,1,2
            do i = 1,10
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
         pwpz(mu) = qp(mu,N) + qz(mu,N)
         pwmz(mu) = qm(mu,N) + qz(mu,N)
      enddo
      pwpz(4) = pwpz(0)**2-pwpz(1)**2-pwpz(2)**2-pwpz(3)**2
      pwmz(4) = pwmz(0)**2-pwmz(1)**2-pwmz(2)**2-pwmz(3)**2

c------------------------------------------------

c get the external quark spinors (including factor sqrt(2E) )

      call psi0m(2,pbar(0,1),sign(1),psi)

c get other bra and ket vectors

      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,qp(0,N),wp(1,N),
     &           psiwp(1,-1,2),fqp(0,2))
      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,qm(0,N),wm(1,N),
     &           psiwm(1,-1,2),fqm(0,2))
      call ket2c(psi(1,-1,1),.true.,p(0,1),-1,qp(0,N),wp(1,N),
     &           psiwp(1,-1,1),fqp(0,1))
      call ket2c(psi(1,-1,1),.true.,p(0,1),-1,qm(0,N),wm(1,N),
     &           psiwm(1,-1,1),fqm(0,1))

      do isig =-1,1,2

         call bra2c(psi(1,isig,2),.true.,p(0,2),isig,qz(0,N),
     &                 atau(1,N), psia(1,isig,2),fqz(0,2))
         call bra2c(psi(1,isig,2),.true.,p(0,2),isig,qz(0,N),
     &                 ztau(1,N), psiz(1,isig,2),fqz(0,2))
         call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qz(0,N),
     &                 atau(1,N), psia(1,isig,1),fqz(0,1))
         call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qz(0,N),
     &                 ztau(1,N), psiz(1,isig,1),fqz(0,1))
      enddo

c***************************************************************
c
c        Tree-level contributions
c
c***************************************************************

c Calculation of processes qq -> W W V ---------------------

c   u u -> W W A
      m3= -s1c(psiwm(1,-1,2),atau(1,N),.true.,-1,psiwp(1,-1,1))
     1          *clr(4,1,-1)
     2    -s1c(psiwm(1,-1,2),wp(1,N),.true.,-1,psia(1,-1,1))*clr(3,1,-1)
     3    -s1c(psia(1,-1,2),wm(1,N),.true.,-1,psiwp(1,-1,1))*clr(3,1,-1)

      mat(1,-1,1) = m3 * clr(3,3,-1)**2 ! uu -> WWA

c   d d -> W W A
      m4= -s1c(psiwp(1,-1,2),atau(1,N),.true.,-1,psiwm(1,-1,1))
     1         *clr(3,1,-1)
     2    -s1c(psiwp(1,-1,2),wm(1,N),.true.,-1,psia(1,-1,1))*clr(4,1,-1)
     3    -s1c(psia(1,-1,2),wp(1,N),.true.,-1,psiwm(1,-1,1))*clr(4,1,-1)

      mat(2,-1,1) = m4 * clr(3,3,-1)**2 ! dd -> WWA

c Calculation of processes qq -> V V -> W W Z -------------

c qq -> A A/Z(-> W W)
      call curr6(1,psi(1,-1,2),p(0,2),
     &           psia(1,-1,1),fqz(0,1), ja(0,-1,1) )

c contract with atoww/ztoww tensor
      do isig = -1,1,2
         m1 = -dotcc(aww(0,N),ja(0,isig,1))
         m2 = -dotcc(zww(0,N),ja(0,isig,1))

         do k = 1,2
            mat(k,isig,2) = (m1*clr(ifl(k),1,isig)
     1                     +  m2*clr(ifl(k),2,isig))
     2                       * clr(ifl(k),1,isig)  ! qq -> AWW
         enddo
      enddo

c qq -> A/Z(-> W W) A
      call curr6(1,psia(1,-1,2),fqz(0,2),
     &           psi(1,-1,1),p(0,1), ja(0,-1,2) )

c contract with atoww/ztoww tensor

      do isig = -1,1,2
         m1 = -dotcc(aww(0,N),ja(0,isig,2))
         m2 = -dotcc(zww(0,N),ja(0,isig,2))
         do k = 1,2
            mat(k,isig,3) = (m1*clr(ifl(k),1,isig)
     1                     +  m2*clr(ifl(k),2,isig))
     2                       * clr(ifl(k),1,isig)  ! qq -> WWA
         enddo
      enddo

c qq -> W+ W-(-> W A)
      call curr6(-1, psiwp(1,-1,2), fqp(0,2),
     &              psi(1,-1,1), p(0,1), jwp(0,-1,2) )

      call curr6(-1, psi(1,-1,2), p(0,2),
     &              psiwp(1,-1,1), fqp(0,1), jwp(0,-1,1) )
      
c contract with wmwa tensor
      m1 = -dotcc(wmwa(0,N),jwp(0,-1,1)) ! uu -> W+ W-(-> W A)
      m2 = -dotcc(wmwa(0,N),jwp(0,-1,2)) ! dd -> W-(-> W A) W+

c additional contributions due to q_mu q_nu/ m^2 terms in the propagator
      m3 = zm2i(3)*dotrc(pwmz(0),wmwa(0,N))*dotrc(pwmz(0),jwp(0,-1,1))
      m4 = zm2i(3)*dotrc(pwmz(0),wmwa(0,N))*dotrc(pwmz(0),jwp(0,-1,2))

      mat(1,-1,4) = m1 * clr(3,3,-1)**2 + m3 * clr(3,3,-1)**2
      mat(2,-1,4) = m2 * clr(3,3,-1)**2 + m4 * clr(3,3,-1)**2

c qq -> W- W+(-> W A)
      call curr6(-1, psiwm(1,-1,2), fqm(0,2),
     &              psi(1,-1,1), p(0,1), jwm(0,-1,2) )

      call curr6(-1, psi(1,-1,2), p(0,2),
     &              psiwm(1,-1,1), fqm(0,1), jwm(0,-1,1) )

c contract with wpwa tensor
      m1 = -dotcc(wpwa(0,N),jwm(0,-1,2)) ! uu -> W+(-> W A) W-
      m2 = -dotcc(wpwa(0,N),jwm(0,-1,1)) ! dd -> W- W+(-> W A)

c additional contributions due to q_mu q_nu/ m^2 terms in the propagator
      m3 = zm2i(3)*dotrc(pwpz(0),wpwa(0,N))*dotrc(pwpz(0),jwm(0,-1,2))
      m4 = zm2i(3)*dotrc(pwpz(0),wpwa(0,N))*dotrc(pwpz(0),jwm(0,-1,1))

      mat(1,-1,5) = m1 * clr(3,3,-1)**2 + m3 * clr(3,3,-1)**2
      mat(2,-1,5) = m2 * clr(3,3,-1)**2 + m4 * clr(3,3,-1)**2

c Calculation of processes qq -> A/Z -> W W Z ---------------------
      call curr6(1,psi(1,-1,2),p(0,2),psi(1,-1,1),p(0,1),jqq(0,-1))

c contract with atowwz/ztowwz tensor
      do isig = -1,1,2
         ma = -dotcc(awwa(0,N),jqq(0,isig))
         mz = -dotcc(zwwa(0,N),jqq(0,isig))
         do k = 1,2
            mat(k,isig,6)  = ma*clr(ifl(k),1,isig)
            mat(k,isig,7)  = mz*clr(ifl(k),2,isig)
         enddo
      enddo

c**************************************************************
c
c    Virtual contributions
c
c**************************************************************

#ifdef WITH_NLO
c Calculation of processes qq -> W W V ---------------------

c calculate pentagon contributions
      if (lpent) then
         scale = -2d0*dotrr(p(0,1),p(0,2))

c calculate born matrix element and store it in mborn for gauge check:
      mborn1=mat(1,-1,1)+mat(1,-1,2)+mat(1,-1,3)+mat(1,-1,4)+mat(1,-1,5)
     1     +mat(1,-1,6)+mat(1,-1,7)
      mborn2=mat(2,-1,1)+mat(2,-1,2)+mat(2,-1,3)+mat(2,-1,4)+mat(2,-1,5)
     1     +mat(2,-1,6)+mat(2,-1,7)

         do mu = 0,3
            p1(mu) = pbar(mu,1)
            p2(mu) = pbar(mu,2)
            minus_qp(mu) = -qp(mu,1)
            minus_qm(mu) = -qm(mu,1)
            minus_qz(mu) = -p1(mu)-p2(mu)+qp(mu,1)+qm(mu,1)
         enddo

         cte=(3*pi*(0d0,1d0)-4.0d0*pi**2/3.0d0+7.0d0)

         if (NLO.eq.1) then ! complete NLO
            call penlinemmm(p1(0),minus_qp(0),minus_qm(0),minus_qz(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wp(1,1),wm(1,1),atau(1,1),
     2         -1,scale,4,1,pentcontr1,j5pm(3),borntest,jbpm(3))
            call penlinemmm(p1(0),minus_qp(0),minus_qz(0),minus_qm(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wp(1,1),atau(1,1),wm(1,1),
     2         -1,scale,4,1,pentcontr2,j5pm(2),borntest,jbpm(2))
            call penlinemmm(p1(0),minus_qz(0),minus_qp(0),minus_qm(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau(1,1),wp(1,1),wm(1,1),
     2         -1,scale,4,1,pentcontr3,j5pm(1),borntest,jbpm(1))

            call Wardtest_Pent(psi, p1, p2, minus_qp, minus_qm, minus_qz,
     1                  wp(1,1), wm(1,1), atau(1,1), scale, als(1,1),
     2                  pentcontr1, mborn1, wardidtest1)
            call Wardtest_Pent(psi, p1, p2, minus_qp, minus_qz, minus_qm,
     1                  wp(1,1), atau(1,1), wm(1,1), scale, als(1,1),
     2                  pentcontr2, mborn1, wardidtest2)
            call Wardtest_Pent(psi, p1, p2, minus_qz, minus_qp, minus_qm,
     1                  atau(1,1), wp(1,1), wm(1,1), scale, als(1,1), 
     2                  pentcontr3, mborn1, wardidtest3)

            if ( (wardidtest1.eq.1) .or. (wardidtest2.eq.1) 
     1                     .or. (wardidtest3.eq.1) ) then
               do i = 1,3
                  j5pm(i)= (0.0d0, 0.0d0)
               enddo
            else
               do i = 1,3
                  j5pm(i)=j5pm(i)+cte*jbpm(i)
               enddo
            endif

            call penlinemmm(p1(0),minus_qm(0),minus_qp(0),minus_qz(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wm(1,1),wp(1,1),atau(1,1),
     2         -1,scale,4,1,pentcontr1,j5mp(3),borntest,jbmp(3))
            call penlinemmm(p1(0),minus_qm(0),minus_qz(0),minus_qp(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wm(1,1),atau(1,1),wp(1,1),
     2         -1,scale,4,1,pentcontr2,j5mp(2),borntest,jbmp(2))
            call penlinemmm(p1(0),minus_qz(0),minus_qm(0),minus_qp(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau(1,1),wm(1,1),wp(1,1),
     2        -1,scale,4,1,pentcontr3,j5mp(1),borntest,jbmp(1))

            call Wardtest_Pent(psi, p1, p2, minus_qm, minus_qp, minus_qz,
     1                  wm(1,1), wp(1,1), atau(1,1), scale, als(1,1),
     2                  pentcontr1, mborn2, wardidtest1)
            call Wardtest_Pent(psi, p1, p2, minus_qm, minus_qz, minus_qp,
     1                  wm(1,1), atau(1,1), wp(1,1), scale, als(1,1),
     2                  pentcontr2, mborn2, wardidtest2)
            call Wardtest_Pent(psi, p1, p2, minus_qz, minus_qm, minus_qp, 
     1                  atau(1,1), wm(1,1), wp(1,1), scale, als(1,1),
     2                  pentcontr3, mborn2, wardidtest3)

            if ( (wardidtest1.eq.1) .or. (wardidtest2.eq.1) 
     1                     .or. (wardidtest3.eq.1) ) then
               do i = 1,3
                  j5mp(i)= (0.0d0, 0.0d0)
               enddo
            else
                do i = 1,3
                   j5mp(i)=j5mp(i)+cte*jbmp(i)
                enddo
            endif

         elseif (NLO.eq.-5) then ! only true pentagon
            call penlinemmm(p1(0),minus_qp(0),minus_qm(0),minus_qz(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wpp(1,1),wmp(1,1),ataup(1,1),
     2         -1,scale,4,1,pentcontr1,j5pm(3),borntest,jbpm(3))
            call penlinemmm(p1(0),minus_qp(0),minus_qz(0),minus_qm(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wpp(1,1),ataup(1,1),wmp(1,1),
     2         -1,scale,4,1,pentcontr2,j5pm(2),borntest,jbpm(2))
            call penlinemmm(p1(0),minus_qz(0),minus_qp(0),minus_qm(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),ataup(1,1),wpp(1,1),wmp(1,1),
     2         -1,scale,4,1,pentcontr3,j5pm(1),borntest,jbpm(1))

            call Wardtest_Pent(psi, p1, p2, minus_qp, minus_qm, minus_qz,
     1                  wpp(1,1), wmp(1,1), ataup(1,1), scale, als(1,1),
     2                  pentcontr1, mborn1, wardidtest1)
            call Wardtest_Pent(psi, p1, p2, minus_qp, minus_qz, minus_qm,
     1                  wpp(1,1), ataup(1,1), wmp(1,1), scale, als(1,1),
     2                  pentcontr2, mborn1, wardidtest2)
            call Wardtest_Pent(psi, p1, p2, minus_qz, minus_qp, minus_qm,
     1                  ataup(1,1), wpp(1,1), wmp(1,1), scale, als(1,1), 
     2                  pentcontr3, mborn1, wardidtest3)

            if ( (wardidtest1.eq.1) .or. (wardidtest2.eq.1) 
     1                     .or. (wardidtest3.eq.1) ) then
               do i = 1,3
                  j5pm(i)= (0.0d0, 0.0d0)
               enddo
            else
               do i = 1,3
                  j5pm(i)=j5pm(i)+cte*jbpm(i)
               enddo
            endif

            call penlinemmm(p1(0),minus_qm(0),minus_qp(0),minus_qz(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wmp(1,1),wpp(1,1),ataup(1,1),
     2         -1,scale,4,1,pentcontr1,j5mp(3),borntest,jbmp(3))
            call penlinemmm(p1(0),minus_qm(0),minus_qz(0),minus_qp(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wmp(1,1),ataup(1,1),wpp(1,1),
     2         -1,scale,4,1,pentcontr2,j5mp(2),borntest,jbmp(2))
            call penlinemmm(p1(0),minus_qz(0),minus_qm(0),minus_qp(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),ataup(1,1),wmp(1,1),wpp(1,1),
     2        -1,scale,4,1,pentcontr3,j5mp(1),borntest,jbmp(1))

            call Wardtest_Pent(psi, p1, p2, minus_qm, minus_qp, minus_qz,
     1                  wmp(1,1), wpp(1,1), ataup(1,1), scale, als(1,1),
     2                  pentcontr1, mborn2, wardidtest1)
            call Wardtest_Pent(psi, p1, p2, minus_qm, minus_qz, minus_qp,
     1                  wmp(1,1), ataup(1,1), wpp(1,1), scale, als(1,1),
     2                  pentcontr2, mborn2, wardidtest2)
            call Wardtest_Pent(psi, p1, p2, minus_qz, minus_qm, minus_qp, 
     1                  ataup(1,1), wmp(1,1), wpp(1,1), scale, als(1,1),
     2                  pentcontr3, mborn2, wardidtest3)

            if ( (wardidtest1.eq.1) .or. (wardidtest2.eq.1) 
     1                     .or. (wardidtest3.eq.1) ) then
               do i = 1,3
                  j5mp(i)= (0.0d0, 0.0d0)
               enddo
            else
                do i = 1,3
                   j5mp(i)=j5mp(i)+cte*jbmp(i)
                enddo
            endif

         else ! pentagon parts ~ born
            call pent_to_box(p1, minus_qp, minus_qm, minus_qz, p2,
     1         psi(1,-1,2), psi(1,-1,1), wp(1,1), wm(1,1),atau(1,1),
     2         xp, xm, xa, -1, scale, j5pm(1), jbpm(1))
            call pent_to_box(p1, minus_qp, minus_qz, minus_qm, p2,
     1         psi(1,-1,2), psi(1,-1,1), wp(1,1), atau(1,1),wm(1,1),
     2         xp, xa, xm, -1, scale, j5pm(2), jbpm(2))
            call pent_to_box(p1, minus_qz, minus_qp, minus_qm, p2,
     1         psi(1,-1,2), psi(1,-1,1), atau(1,1), wp(1,1),wm(1,1),
     2         xa, xp, xm, -1, scale, j5pm(3), jbpm(3))

            call pent_to_box(p1, minus_qm, minus_qp, minus_qz, p2,
     1         psi(1,-1,2), psi(1,-1,1), wm(1,1), wp(1,1),atau(1,1),
     2         xm, xp, xa, -1, scale, j5mp(1), jbmp(1))
            call pent_to_box(p1, minus_qm, minus_qz, minus_qp, p2,
     1         psi(1,-1,2), psi(1,-1,1), wm(1,1), atau(1,1),wp(1,1),
     2         xm, xa, xp, -1, scale, j5mp(2), jbmp(2))
            call pent_to_box(p1, minus_qz, minus_qm, minus_qp, p2,
     1         psi(1,-1,2), psi(1,-1,1), atau(1,1), wm(1,1),wp(1,1),
     2         xa, xm, xp, -1, scale, j5mp(3), jbmp(3))

            do i = 1,3
               j5pm(i) = j5pm(i)+cte*jbpm(i)
               j5mp(i) = j5mp(i)+cte*jbmp(i)
            enddo

         endif !NLO

c   u u -> W W A
         mv3 = j5pm(1)*clr(3,1,-1)
     1       + j5pm(2)*clr(4,1,-1)
     2       + j5pm(3)*clr(3,1,-1)
         matv(1,-1,1) = mv3 * clr(3,3,-1)**2

c   d d -> W W A
         mv4 = j5mp(1)*clr(4,1,-1)
     1       + j5mp(2)*clr(3,1,-1)
     2       + j5mp(3)*clr(4,1,-1)

         matv(2,-1,1) = mv4 * clr(3,3,-1)**2

      endif !lpent

c Calculation of processes qq -> V V -> W W Z -------------

      if (lbox) then
         scale = -2d0*dotrr(p(0,1),p(0,2))

c calculate born matrix element and store it in mborn for gauge check:
      mborn1=mat(1,-1,1)+mat(1,-1,2)+mat(1,-1,3)+mat(1,-1,4)+mat(1,-1,5)
     1     +mat(1,-1,6)+mat(1,-1,7)
      mborn2=mat(2,-1,1)+mat(2,-1,2)+mat(2,-1,3)+mat(2,-1,4)+mat(2,-1,5)
     1     +mat(2,-1,6)+mat(2,-1,7)

         do mu = 0,3
            p1(mu) = pbar(mu,1)
            p2(mu) = pbar(mu,2)
            minus_qp(mu) = -qp(mu,1)
            minus_qm(mu) = -qm(mu,1)
            minus_qwmz(mu) = -p1(mu)-p2(mu)+qp(mu,1)
            minus_qwpz(mu) = -p1(mu)-p2(mu)+qm(mu,1)
            minus_qww(mu) = -qp(mu,1)-qm(mu,1)
            minus_qz(mu) = -p1(mu)-p2(mu)+qp(mu,1)+qm(mu,1)
         enddo

         cte=(3*pi*(0d0,1d0)-4.0d0*pi**2/3.0d0+7.0d0)

c qq -> A A/Z(-> W W)
         call boxlinemm(p1,minus_qz,minus_qww,p2,
     1    psi(1,-1,2),psi(1,-1,1),atau(1,1),zww(0,1),
     2    -1,scale,9,1,boxtest,result(-1,2),treetest,
     3    result_born(-1,2))
         call boxlinemm(p1,minus_qz,minus_qww,p2,
     1    psi(1,1,2),psi(1,1,1),atau(1,1),zww(0,1),
     2    1,scale,9,-1,boxtest,result(1,2),treetest,
     3    result_born(1,2))
         call boxlinemm(p1,minus_qz,minus_qww,p2,
     1    psi(1,-1,2),psi(1,-1,1),atau(1,1),aww(0,1),
     2    -1,scale,3,-1,gaugebox1,result(-1,1),gaugebox1_born,
     3    result_born(-1,1))
         call boxlinemm(p1,minus_qz,minus_qww,p2,
     1    psi(1,1,2),psi(1,1,1),atau(1,1),aww(0,1),
     2    1,scale,9,-1,boxtest,result(1,1),treetest,
     3    result_born(1,1))

         call Wardtest_Box(atau(1,1), qz(0,1), aww(0,1), qww(0,1),
     1            gaugebox1, gaugebox1_born, mborn1, wardidtest1, als(1,1))

         if (wardidtest1.eq.1) then
            do isig = -1,1,2
               do k = 1,2
                  matv(k,isig,2) = (0.0d0, 0.0d0)
               enddo
            enddo
         else
            do isig = -1,1,2
               mv1=-(result(isig,1)+cte*result_born(isig,1))
               mv2=-(result(isig,2)+cte*result_born(isig,2))
               do k = 1,2
                  matv(k,isig,2) = (mv1*clr(ifl(k),1,isig)
     1                     +  mv2*clr(ifl(k),2,isig))
     2                       * clr(ifl(k),1,isig)
               enddo
            enddo
         endif

c qq -> A/Z(-> W W) A
         call boxlinemm(p1,minus_qww,minus_qz,p2,
     1    psi(1,-1,2),psi(1,-1,1),aww(0,1), atau(1,1),
     2    -1,scale,3,1,gaugebox1,result(-1,1),gaugebox1_born,
     3    result_born(-1,1))
         call boxlinemm(p1,minus_qww,minus_qz,p2,
     1    psi(1,1,2),psi(1,1,1),aww(0,1), atau(1,1),
     2    1,scale,9,-1,boxtest,result(1,1),treetest,
     3    result_born(1,1))
         call boxlinemm(p1,minus_qww,minus_qz,p2,
     1    psi(1,-1,2),psi(1,-1,1),zww(0,1), atau(1,1),
     2    -1,scale,9,-1,boxtest,result(-1,2),treetest,
     3    result_born(-1,2))
         call boxlinemm(p1,minus_qww,minus_qz,p2,
     1    psi(1,1,2),psi(1,1,1),zww(0,1), atau(1,1),
     2    1,scale,9,-1,boxtest,result(1,2),treetest,
     3    result_born(1,2))

         call Wardtest_Box(aww(0,1), qww(0,1), atau(1,1), qz(0,1),
     1            gaugebox1, gaugebox1_born, mborn1, wardidtest1, als(1,1))

         if (wardidtest1.eq.1) then
            do isig = -1,1,2
               do k = 1,2
                  matv(k,isig,3) = (0.0d0, 0.0d0)
               enddo
            enddo
         else
            do isig = -1,1,2
               mv1=-(result(isig,1)+cte*result_born(isig,1))
               mv2=-(result(isig,2)+cte*result_born(isig,2))
               do k = 1,2
                  matv(k,isig,3) = (mv1*clr(ifl(k),1,isig)
     1                     +  mv2*clr(ifl(k),2,isig))
     2                       * clr(ifl(k),1,isig)
               enddo
            enddo
         endif

c qq -> W+ W-(-> W A)
         call boxlinemm(p1,minus_qp,minus_qwmz,p2,
     1    psi(1,-1,2),psi(1,-1,1),wp(1,1), wmwa(0,1),
     2    -1,scale,3,1,gaugebox1,result(-1,1),gaugebox1_born,
     3    result_born(-1,1))
         call boxlinemm(p1,minus_qwmz,minus_qp,p2,
     1    psi(1,-1,2),psi(1,-1,1),wmwa(0,1), wp(1,1),
     2    -1,scale,3,1,gaugebox2,result(-1,2),gaugebox2_born,
     3    result_born(-1,2))

         call Wardtest_Box(wp(1,1), qp(0,1), wmwa(0,1), pwmz(0),
     1            gaugebox1, gaugebox1_born, mborn1, wardidtest1, als(1,1))
         call Wardtest_Box(wmwa(0,1), pwmz(0), wp(1,1), qp(0,1),
     1            gaugebox2, gaugebox2_born, mborn2, wardidtest2, als(1,1))

         if (wardidtest1.eq.1) then
            matv(1,-1,4) = (0.0d0, 0.0d0)
         else
            mv1=-(result(-1,1)+cte*result_born(-1,1))
            matv(1,-1,4) = mv1 * clr(3,3,-1)**2
         endif
         if (wardidtest2.eq.1) then
            matv(2,-1,4) = (0.0d0, 0.0d0)
         else
            mv2=-(result(-1,2)+cte*result_born(-1,2))
            matv(2,-1,4) = mv2 * clr(3,3,-1)**2
         endif

c qq -> W- W+(-> W A)
          call boxlinemm(p1,minus_qm,minus_qwpz,p2,
     1    psi(1,-1,2),psi(1,-1,1),wm(1,1), wpwa(0,1),
     2    -1,scale,3,1,gaugebox1,result(-1,1),gaugebox1_born,
     3    result_born(-1,1))
         call boxlinemm(p1,minus_qwpz,minus_qm,p2,
     1    psi(1,-1,2),psi(1,-1,1),wpwa(0,1), wm(1,1),
     2    -1,scale,3,1,gaugebox2,result(-1,2),gaugebox2_born,
     3    result_born(-1,2))

         call Wardtest_Box(wm(1,1), qm(0,1), wpwa(0,1), pwpz(0),
     1            gaugebox1, gaugebox1_born, mborn1, wardidtest1, als(1,1))
         call Wardtest_Box(wpwa(0,1), pwpz(0), wm(1,1), qm(0,1),
     1            gaugebox2, gaugebox2_born, mborn2, wardidtest2, als(1,1))

         if (wardidtest1.eq.1) then
            matv(1,-1,5) = (0.0d0, 0.0d0)
         else
            mv1=-(result(-1,1)+cte*result_born(-1,1))
            matv(1,-1,5) = mv1 * clr(3,3,-1)**2
         endif
         if (wardidtest2.eq.1) then
            matv(2,-1,5) = (0.0d0, 0.0d0)
         else
            mv2=-(result(-1,2)+cte*result_born(-1,2))
            matv(2,-1,5) = mv2 * clr(3,3,-1)**2
         endif

      endif

c Calculation of processes qq -> A/Z -> W W Z ---------------------
      do isig = -1,1,2
         do k = 1,2
            matv(k,isig,6) = (0.0d0, 0.0d0)
            matv(k,isig,7) = (0.0d0,0.0d0)
         enddo
      enddo

#endif /* WITH_NLO */
      
c-----------------------------------------------------
c------------------------------------------------------
#ifdef WITH_MADGRAPH
c comparison with madgraph
      if (ldebugm) then

         do isig = -1,1,2
            do i = 1,7
               do k = 1,2

                  madgr(i,isig,k) = 0d0

               enddo
            enddo
         enddo

         do isig = -1,1,2
            do k = 1,2
               do i =1,7

                  madgr(i,isig,k) = tamp(i,isig,k)

               enddo
            enddo
         enddo

c print results from checks:

         do i=1,7
            do k = 1,2
               do isig = -1,1,2
                  print*,'mat(',k,',',isig,',',i,') = ',mat(k,isig,i)
                  temp = dreal(mat(k,isig,i))**2+dimag(mat(k,isig,i))**2
                  print*,'madgr(',i,',',isig,',',k,') = ',madgr(i,isig,k)
                  temp1 = dreal(madgr(i,isig,k))**2+dimag(madgr(i,isig,k))**2
                  print*,'temp/temp1 = ',temp/temp1
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
            do i = 1,10
               mm(k,isig) = mm(k,isig) +  mat(k,isig,i)
            enddo
            res(k)= res(k) + dreal(mm(k,isig))**2
     &                     + dimag(mm(k,isig))**2

c  add Born type term and multiply by F_q = alphas*C_2/4pi
c  the additional factor for the born term is after adding the 
c  subtraction term and the counter term for the renormalization
c  of the pdfs (see also D(x) in m2s_qqWWZj.f)

#ifdef WITH_NLO
            if (lnlo) then
               mv(k,isig)=matv(k,isig,1)+matv(k,isig,2)
     #              +matv(k,isig,3)
     #              +matv(k,isig,4)+matv(k,isig,5)+matv(k,isig,6)
     #              +matv(k,isig,7)
            endif
#endif
            if (nlo.gt.0) then
#ifdef WITH_NLO
               mv(k,isig) = als(1,1)*c2o4pi*
     1                  (mv(k,isig) + mm(k,isig)* cvirtVV) 
#endif
            else if(nlo.lt.0) then
#ifdef WITH_NLO
               mv(k,isig) = als(1,1)*c2o4pi*mv(k,isig)
#endif
            else
               mv(k,isig) = (0d0, 0d0)
            endif

           resv(k) = resv(k) + 2d0*dreal( mm(k,isig)   *
     1              conjg( mv(k,isig) ) )

         enddo

         if (nlo.ge.0) then
            res(k) = (res(k)+resv(k))*3d0      ! 3 is the color sum factor
         else
            res(k) = resv(k)*3d0               ! 3 is the color sum factor
         endif  
      enddo

      uuwwa = res(1)
      ddwwa = res(2)

      return
      end

