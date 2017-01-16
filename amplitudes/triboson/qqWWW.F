


c----------------------------------------------------------------------------
c----------------------------------------------------------------------------
c
c        Stefan Prestel : sprestel@particle.uni-karlsruhe.de
c        Initial version: February 2008
c
c        This file contains subroutines to calculate matrix elememts**2 for:
c
c     q1 q2~ --> W+W-W+, 
c       W+(1) ---> f1-bar f2, W-(1) ---> f3-bar f4, W+(2) ---> f5-bar f6
c     subroutine qqWWWp
c
c     q1 q2~ --> W-W+W-, 
c       W-(1) ---> f1-bar f2, W+(1) ---> f3-bar f4, W-(2) ---> f5-bar f6
c     subroutine qqWWWm
c
c----------------------------------------------------------------------------
c----------------------------------------------------------------------------


c****************subroutine QQWWWp********************************************
c
c  QQWWWp calculates the matrix elements**2 for 
C        q1 q2    ---->   W+ W- W+,
C         W+ ---> f1-bar f2, W- ---> f3-bar f4, W+ ---> f5-bar f6
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
c          NW = 1        calculation in narrow width approximation
c             = 0        full calculation
c******************************************************************************

      subroutine qqWWWp(pbar,sign, nlo, N, qqbarWWW,NW)

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/scales.inc"

c     electroweak couplings are taken from KOPPLN
      double precision  clr, xm2, xmg, b, v, a
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     1                  V(4,5),A(4,5)
      double precision  pbar(0:3,2+n_v)                              ! input momenta
      double precision  qqbarWWW, res,resv                           ! output M^2 of process
      double precision  p(0:3,2+n_v), pwp1wm1(0:4), pwp2wm1(0:4),    ! different momenta
     1                  fqp1(0:4,2), fqp2(0:4,2), fqm1(0:4,2)        ! and momentum combinations
      integer  sign(2+n_v), nlo, mu, i, N
      integer NW
      double complex psiwp1(2,-1:1,2), psiwp2(2,-1:1,2),             ! q+W spinors, e.g. <q~W+(1)|,
     1               psiwm1(2,-1:1,2),                               ! and currents
     2               jwp1(0:5,-1:1,2), jwp2(0:5,-1:1,2),
     3               jwm1(0:5,-1:1,2)
      double complex temp, temp1
      double complex mat(6), matv(6), mborn
      double complex madgr(6)
      double complex mm, mv
      double complex m3line1, m2linea1, m2linez1, m1line1
      double complex psi(2,-1:1,2), jqqbar(0:5,-1:1)                 ! spinors <q~| and |q>, current <q~|gamma_mu|q>
      double complex m3line2, m2linez2                     ! matrix variables for virtual corrections
      double complex dotcc, s1c, dotrc                               ! "dotproducts"
      double precision dotrr, psumsq
      external dotcc, s1c, dotrr, dotrc, psumsq

      double precision rzero,rone, epsilon
      double complex czero, cone
      parameter (rzero =  0d0 ,      rone =  1d0, epsilon = 1.0d-1)
      parameter (czero = (0d0,0d0) , cone = (1d0,0d0))

      double complex zm2i(2:3)                                       !variable for q_mu q_nu terms in the
      save zm2i                                                      !propagators of the massive vector bosons

c variables for virtual corrections
      double precision scale
      double precision c2,c2o4pi
      parameter (c2=4d0/3d0, c2o4pi=c2/4d0/pi)
      logical lnlo, lbox, lpq, lpt, lpent, linit
      data linit /.true./
      save linit
      double complex j5wpwm(1:6)
      double complex jbwpwm(1:6)
      double complex cte
      double complex resultgauge1(3), resultgauge3(3)
      double complex resultgaugeb1(3),resultgaugeb3(3)
      double precision p1(0:3),p2(0:3)
      double precision minus_qp1(0:3),minus_qp2(0:3),minus_q1m(0:3)
      double precision minus_qp1m(0:3), minus_qp2m(0:3)
      double complex result(8),result_born(8)
      double complex gaugebox1(2), gaugebox1_born(2), boxtest(2), treetest(2)
      double complex gaugebox2(2), gaugebox2_born(2)
      double complex gaugebox3(2), gaugebox3_born(2)
      double complex gaugebox4(2), gaugebox4_born(2)
      integer wardidtest1, wardidtest2, wardidtest3, wardidtest4

c variables for test purposes
      logical ldebugm
      parameter( ldebugm = .false.)

c---------------------------- BEGIN CODE --------------------------------------

      lnlo = NLO.ne.0       ! include some virtual stuff if T
      lbox = (NLO.eq.-4) .or. (NLO.eq.1) .or. (NLO.eq.5)    ! call boxline if T
      lpt = (NLO.eq.1) .or. (NLO.eq.-5)   ! true pentagon contribution
      lpq = (NLO.eq.5) .or. (NLO.eq.-4)   ! q^mu/m_V terms of Pentagon
      lpent = lpt .or. lpq

c define flavors of external quarks for the 2 subprocesses

      if (linit) then
         linit = .false.
         zm2i(2) = 1d0/dcmplx(xm2(2),-xmg(2))
         zm2i(3) = 1d0/dcmplx(xm2(3),-xmg(3))
      endif

      do i = 1,6
         mat(i)  = czero
         matv(i) = czero
      enddo

c define the internal momenta
      do mu = 0,3
         do i = 1,2+n_v
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo
         pwp1wm1(mu) = qp1(mu,N) + q1m(mu,N)
         pwp2wm1(mu) = qp2(mu,N) + q1m(mu,N)
      enddo
      pwp1wm1(4) = dotrr(pwp1wm1(0),pwp1wm1(0))
      pwp2wm1(4) = dotrr(pwp2wm1(0),pwp2wm1(0))

c------------------------------------------------------------------------------

c get the external quark spinors (including factor sqrt(2E) )
      call psi0m(2,pbar(0,1),sign(1),psi)

c get other bra and ket vectors
      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,qp1(0,N),wp1(1,N),psiwp1(1,-1,2),fqp1(0,2))     ! <q~W+(1)|
      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,q1m(0,N),w1m(1,N),psiwm1(1,-1,2),fqm1(0,2))     ! <q~W-(1)|
      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,qp2(0,N),wp2(1,N),psiwp2(1,-1,2),fqp2(0,2))     ! <q~W+(2)|

      call ket2c(psi(1,-1,1),.true.,p(0,1),-1,qp1(0,N),wp1(1,N),psiwp1(1,-1,1),fqp1(0,1))     ! |W+(1)q>
      call ket2c(psi(1,-1,1),.true.,p(0,1),-1,q1m(0,N),w1m(1,N),psiwm1(1,-1,1),fqm1(0,1))     ! |W-(1)q>
      call ket2c(psi(1,-1,1),.true.,p(0,1),-1,qp2(0,N),wp2(1,N),psiwp2(1,-1,1),fqp2(0,1))     ! |W+(2)q>

c get current for contraction with W->WWW leptonic tensor
      call curr6(-1,psi(1,-1,2),p(0,2),psi(1,-1,1),p(0,1),jqqbar(0,-1))                       !<q~|gamma_mu|q>

c get currents for contraction with A/Z->WW leptonic tensor
      call curr6(-1,psi(1,-1,2),p(0,2),psiwp1(1,-1,1),fqp1(0,1),jwp1(0,-1,1))                 !<q~|gamma_mu|W+(1)q>
      call curr6(-1,psi(1,-1,2),p(0,2),psiwp2(1,-1,1),fqp2(0,1),jwp2(0,-1,1))                 !<q~|gamma_mu|W+(2)q>
      call curr6(-1,psi(1,-1,2),p(0,2),psiwm1(1,-1,1),fqm1(0,1),jwm1(0,-1,1))                 !<q~|gamma_mu|W-(1)q>

      call curr6(-1,psiwp1(1,-1,2),fqp1(0,2),psi(1,-1,1),p(0,1),jwp1(0,-1,2))                 !<q~W+(1)|gamma_mu|q>
      call curr6(-1,psiwp2(1,-1,2),fqp2(0,2),psi(1,-1,1),p(0,1),jwp2(0,-1,2))                 !<q~W+(2)|gamma_mu|q>
      call curr6(-1,psiwm1(1,-1,2),fqm1(0,2),psi(1,-1,1),p(0,1),jwm1(0,-1,2))                 !<q~W-(1)|gamma_mu|q>


C------------------------------------------------------------------------------
C CALCULATION OF UD~ --> W+W-W+ --> LEPTONS PROCESSES

c    U D~ --> W+W-W+ --> LEPTONS
      m3line1 = -s1c(psiwp1(1,-1,2),w1m(1,N),.true.,-1,psiwp2(1,-1,1))  !W+(2)W-(1)W+(1)
     1          -s1c(psiwp2(1,-1,2),w1m(1,N),.true.,-1,psiwp1(1,-1,1))  !W+(1)W-(1)W+(2)
      m3line1 = m3line1*clr(3,3,-1)**3

      MAT(1)  = M3LINE1                                               !M(UD~->W+W-W+)

C------------------------------------------------------------------------------
C CALCULATION OF U D~ --> W+ V --> LEPTONS PROCESSES

c    U D~ --> ZW+ --> LEPTONS
      m2linez1 = -dotcc(WWW_zww(0,N,2),jwp1(0,-1,2))*clr(3,2,-1)  !ZW
     1           -dotcc(WWW_zww(0,N,1),jwp2(0,-1,2))*clr(3,2,-1)  !ZW
      MAT(2)   = M2LINEZ1*CLR(3,3,-1)                      !M(UD~->ZW+ ->L)
      m2linez1 = czero

c    additional contributions due to q_mu q_nu/ m^2 terms in the propagator
      m2linez1 = zm2i(2)*dotrc(pwp2wm1(0),WWW_zww(0,N,2))*dotrc(pwp2wm1(0),jwp1(0,-1,2))
     &         + zm2i(2)*dotrc(pwp1wm1(0),WWW_zww(0,N,1))*dotrc(pwp1wm1(0),jwp2(0,-1,2))

      MAT(2)   = MAT(2) + M2LINEZ1*CLR(3,2,-1)*clr(3,3,-1)

c    U D~ --> AW+ --> LEPTONS
      m2linea1 = -dotcc(WWW_aww(0,N,2),jwp1(0,-1,2))*clr(3,1,-1)  !AW
     1           -dotcc(WWW_aww(0,N,1),jwp2(0,-1,2))*clr(3,1,-1)  !AW

      MAT(3)   = M2LINEA1*CLR(3,3,-1)                      !M(UD~->AW+ ->L)

c    U D~ --> W+Z --> LEPTONS
      m2linez1 = -dotcc(WWW_zww(0,N,2),jwp1(0,-1,1))*clr(4,2,-1)  !WZ
     1           -dotcc(WWW_zww(0,N,1),jwp2(0,-1,1))*clr(4,2,-1)  !WZ

      MAT(4)   = M2LINEZ1*CLR(3,3,-1)                      !M(UD~->W+Z ->W+W-W+)

c    additional contributions due to q_mu q_nu/ m^2 terms in the propagator
      m2linez1 = zm2i(2)*dotrc(pwp2wm1(0),WWW_zww(0,N,2))*dotrc(pwp2wm1(0),jwp1(0,-1,1))
     &         + zm2i(2)*dotrc(pwp1wm1(0),WWW_zww(0,N,1))*dotrc(pwp1wm1(0),jwp2(0,-1,1))
      MAT(4)   = MAT(4) + M2LINEZ1*CLR(4,2,-1)*clr(3,3,-1)

c    U D~ --> W+A --> LEPTONS
      m2linea1 = -dotcc(WWW_aww(0,N,2),jwp1(0,-1,1))*clr(4,1,-1)  !WA
     1           -dotcc(WWW_aww(0,N,1),jwp2(0,-1,1))*clr(4,1,-1)  !WA

      MAT(5)   = M2LINEA1*CLR(3,3,-1)                      !M(UD~->W+A ->W+W-W+)

C------------------------------------------------------------------------------
C CALCULATION OF UD~ --> W+ --> LEPTONS PROCESSES

c    U D~ --> W+ --> LEPTONS
      m1line1 = -dotcc(wptowww(0,N),jqqbar(0,-1))
      m1line1 = m1line1*clr(3,3,-1)

      MAT(6)  = - M1LINE1                                     !M(UD~->W+->L)

c*****************************************************************************
c
c                        PENTAGON AND BOX CONTRIBUTIONS
c
c*****************************************************************************

#ifdef WITH_NLO
      if(lpent) then

         scale = -2.0d0*dotrr(p(0,1),p(0,2))

c calculate born matrix element and store it in mborn for gauge check:
         mborn=mat(1)+mat(2)+mat(3)+mat(4)+mat(5)+mat(6)

         do mu=0,3
            p1(mu) = pbar(mu,1)
            p2(mu) = pbar(mu,2)
            minus_qp1(mu) = -qp1(mu,1)
            minus_qp2(mu) = -qp2(mu,1)
            minus_q1m(mu) = -p1(mu) -p2(mu) -minus_qp1(mu) -minus_qp2(mu)
         enddo

         cte = (3*pi*(0d0,1d0)-4.0d0*pi**2/3.0d0+7.0d0)

         if(nlo.eq.1) then
            call penlinemmm(p1(0),minus_qp1(0),minus_q1m(0),minus_qp2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wp1(1,1),w1m(1,1),wp2(1,1),              ! W+(1) W- W+(2)
     2         -1,scale,4,1,resultgauge3(1),j5wpwm(3),resultgaugeb3(1),jbwpwm(3))
            call penlinemmm(p1(0),minus_qp2(0),minus_q1m(0),minus_qp1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wp2(1,1),w1m(1,1),wp1(1,1),              ! W+(2) W- W+(1)
     2         -1,scale,4,1,resultgauge1(1),j5wpwm(1),resultgaugeb1(1),jbwpwm(1))

            call Wardtest_Pent(psi, p1, p2, minus_qp1, minus_q1m, minus_qp2,
     1               wp1(1,1), w1m(1,1), wp2(1,1), scale, als(1,1),
     2               resultgauge3(1), mborn, wardidtest1)

            call Wardtest_Pent(psi, p1, p2, minus_qp2, minus_q1m, minus_qp1,
     1               wp2(1,1), w1m(1,1), wp1(1,1), scale, als(1,1),
     2               resultgauge1(1), mborn, wardidtest2)

            if (wardidtest1.eq.1 .or. wardidtest2.eq.1) then
               j5wpwm(1) = (0.0d0, 0.0d0)
               j5wpwm(3) = (0.0d0, 0.0d0)
            else
               j5wpwm(1)=  j5wpwm(1)+cte*jbwpwm(1)
               j5wpwm(3)=  j5wpwm(3)+cte*jbwpwm(3)
            endif

         elseif(nlo.eq.-5) then
            call penlinemmm(p1(0),minus_qp1(0),minus_q1m(0),minus_qp2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wpp1(1,1),w1mp(1,1),wpp2(1,1),              ! W+(1) W- W+(2)
     2         -1,scale,4,1,resultgauge3(1),j5wpwm(3),resultgaugeb3(1),jbwpwm(3))
            call penlinemmm(p1(0),minus_qp2(0),minus_q1m(0),minus_qp1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wpp2(1,1),w1mp(1,1),wpp1(1,1),              ! W+(2) W- W+(1)
     2         -1,scale,4,1,resultgauge1(1),j5wpwm(1),resultgaugeb1(1),jbwpwm(1))

            call Wardtest_Pent(psi, p1, p2, minus_qp1, minus_q1m, minus_qp2,
     1               wpp1(1,1), w1mp(1,1), wpp2(1,1), scale, als(1,1),
     2               resultgauge3(1), mborn, wardidtest1)

            call Wardtest_Pent(psi, p1, p2, minus_qp2, minus_q1m, minus_qp1,
     1               wpp2(1,1), w1mp(1,1), wpp1(1,1), scale, als(1,1),
     2               resultgauge1(1), mborn, wardidtest2)

            if (wardidtest1.eq.1 .or. wardidtest2.eq.1) then
               j5wpwm(1) = (0.0d0, 0.0d0)
               j5wpwm(3) = (0.0d0, 0.0d0)
            else
               j5wpwm(1)=  j5wpwm(1)+cte*jbwpwm(1)
               j5wpwm(3)=  j5wpwm(3)+cte*jbwpwm(3)
            endif
            
         else
c pentagon contributions shifted to boxes by shifting polarisation vectors
            call pent_to_box(p1, minus_qp1, minus_q1m, minus_qp2, p2,
     1         psi(1,-1,2), psi(1,-1,1), wp1(1,1), w1m(1,1),wp2(1,1),
     2         xp1, x1m, xp2, -1, scale, j5wpwm(3), jbwpwm(3))
            call pent_to_box(p1, minus_qp2, minus_q1m, minus_qp1, p2,
     1         psi(1,-1,2), psi(1,-1,1), wp2(1,1), w1m(1,1),wp1(1,1),
     2         xp2, x1m, xp1, -1, scale, j5wpwm(1), jbwpwm(1))

               j5wpwm(1)=  j5wpwm(1)+cte*jbwpwm(1)
               j5wpwm(3)=  j5wpwm(3)+cte*jbwpwm(3)
         endif

c    U D~ --> W+W-W+ --> LEPTONS
         m3line2 = j5wpwm(1) + j5wpwm(3)
         m3line2 = m3line2*clr(3,3,-1)**3

         MATV(1) = M3LINE2

      endif !lpent


      if(lbox) then
         scale = -2.0d0*dotrr(p(0,1),p(0,2))

c calculate born matrix element and store it in mborn for gauge check:
         mborn=mat(1)+mat(2)+mat(3)+mat(4)+mat(5)+mat(6)

         do mu=0,3
            p1(mu) = pbar(mu,1)
            p2(mu) = pbar(mu,2)
            minus_qp1(mu) = -qp1(mu,1)
            minus_qp2(mu) = -qp2(mu,1)
            minus_qp2m(mu) = -p1(mu) -p2(mu) +qp1(mu,1)
            minus_qp1m(mu) = -p1(mu) -p2(mu) +qp2(mu,1)
         enddo

         cte = (3*pi*(0d0,1d0)-4.0d0*pi**2/3.0d0+7.0d0)

         call boxlinemm(p1,minus_qp2m,minus_qp1,p2,
     1    psi(1,-1,2),psi(1,-1,1),WWW_zww(0,1,2),wp1(1,1),
     2    -1,scale,3,1,gaugebox1,result(1),gaugebox1_born,
     3    result_born(1))
         call boxlinemm(p1,minus_qp2m,minus_qp1,p2,
     1    psi(1,-1,2),psi(1,-1,1),WWW_aww(0,1,2),wp1(1,1),
     2    -1,scale,9,-1,boxtest,result(2),treetest,
     3    result_born(2))

         call boxlinemm(p1,minus_qp1m,minus_qp2,p2,
     1    psi(1,-1,2),psi(1,-1,1),WWW_zww(0,1,1),wp2(1,1),
     2    -1,scale,3,1,gaugebox2,result(3),gaugebox2_born,
     3    result_born(3))
         call boxlinemm(p1,minus_qp1m,minus_qp2,p2,
     1    psi(1,-1,2),psi(1,-1,1),WWW_aww(0,1,1),wp2(1,1),
     2    -1,scale,9,-1,boxtest,result(4),treetest,
     3    result_born(4))

         call boxlinemm(p1,minus_qp1,minus_qp2m,p2,
     1    psi(1,-1,2),psi(1,-1,1),wp1(1,1),WWW_zww(0,1,2),
     2    -1,scale,3,1,gaugebox3,result(5),gaugebox3_born,
     3    result_born(5))
         call boxlinemm(p1,minus_qp1,minus_qp2m,p2,
     1    psi(1,-1,2),psi(1,-1,1),wp1(1,1),WWW_aww(0,1,2),
     2    -1,scale,9,-1,boxtest,result(6),treetest,
     3    result_born(6))

         call boxlinemm(p1,minus_qp2,minus_qp1m,p2,
     1    psi(1,-1,2),psi(1,-1,1),wp2(1,1),WWW_zww(0,1,1),
     2    -1,scale,3,1,gaugebox4,result(7),gaugebox4_born,
     3    result_born(7))
         call boxlinemm(p1,minus_qp2,minus_qp1m,p2,
     1    psi(1,-1,2),psi(1,-1,1),wp2(1,1),WWW_aww(0,1,1),
     2    -1,scale,9,-1,boxtest,result(8),treetest,
     3    result_born(8))

         call Wardtest_Box(WWW_zww(0,1,2), pwp2wm1(0), wp1(1,1), qp1(0,1),
     1            gaugebox1, gaugebox1_born, mborn, wardidtest1, als(1,1))

         call Wardtest_Box(WWW_zww(0,1,1), pwp1wm1(0), wp2(1,1), qp2(0,1),
     1            gaugebox2, gaugebox2_born, mborn, wardidtest2, als(1,1))

         call Wardtest_Box(wp1(1,1), qp1(0,1), WWW_zww(0,1,2), pwp2wm1(0),
     1            gaugebox3, gaugebox3_born, mborn, wardidtest3, als(1,1))

         call Wardtest_Box(wp2(1,1), qp2(0,1), WWW_zww(0,1,1), pwp1wm1(0),
     1            gaugebox4, gaugebox4_born, mborn, wardidtest4, als(1,1))

          result(1) = -(result(1)+cte*result_born(1))
          result(2) = -(result(2)+cte*result_born(2))
          result(3) = -(result(3)+cte*result_born(3))
          result(4) = -(result(4)+cte*result_born(4))
          result(5) = -(result(5)+cte*result_born(5))
          result(6) = -(result(6)+cte*result_born(6))
          result(7) = -(result(7)+cte*result_born(7))
          result(8) = -(result(8)+cte*result_born(8))

         if (wardidtest1.eq.0) then
            m2linez2 = result(1)*clr(3,2,-1)
            MATV(2)  = M2LINEZ2*CLR(3,3,-1)
            m2linez2 = result(2)*clr(3,1,-1)
            MATV(3)  = M2LINEZ2*CLR(3,3,-1)
         endif

         if (wardidtest2.eq.0) then
            m2linez2 = result(3)*clr(3,2,-1)
            MATV(2)  = M2LINEZ2*CLR(3,3,-1)
            m2linez2 = result(4)*clr(3,1,-1)
            MATV(3)  = M2LINEZ2*CLR(3,3,-1)
         endif

         if (wardidtest3.eq.0) then
            m2linez2 = result(5)*clr(4,2,-1)
            MATV(4)  = M2LINEZ2*CLR(3,3,-1)
            m2linez2 = result(6)*clr(4,1,-1)
            MATV(5)  = M2LINEZ2*CLR(3,3,-1)
         endif

         if (wardidtest4.eq.0) then
            m2linez2 = result(7)*clr(4,2,-1)
            MATV(4)  = M2LINEZ2*CLR(3,3,-1)
            m2linez2 = result(8)*clr(4,1,-1)
            MATV(5)  = M2LINEZ2*CLR(3,3,-1)
         endif

      endif !lbox

      MATV(6) = CZERO
#endif

c------------------------------------------------------------------------

c comparison with madgraph

#ifdef WITH_MADGRAPH
      if (ldebugm) then
         do i = 1,6
            madgr(i) = czero
         enddo

         do i =1,6
            madgr(i) = WWWtamp(i,-1,1)
         enddo

c print results from checks:

         print*,'------------------------------------------------------'
         do i = 1,6
    	    print*,'ratio for i = ',i,':',
     &	       '(', dreal(mat(i))/dreal(madgr(i)) ,
     &              dimag(mat(i))/dimag(madgr(i)), ')'
            print*
            print*,'mat(',i,') = ',mat(i)
            temp = dreal(mat(i))**2+dimag(mat(i))**2
            print*,temp
            print*,'madgr(',i,') = ',madgr(i)
            temp1 = dreal(madgr(i))**2+dimag(madgr(i))**2
            print*,temp1
            print*,'temp/temp1 = ',temp/temp1
            print*,' '
         enddo
      endif ! end of debugging block
#endif

c------------------------------------------------------------------------------
c------------------------------------------------------------------------------

c sum the graphs, square them and map them onto qqbarWWW

      res = 0.0d0
      resv = 0.0d0
      mm = czero
      mv = czero
      do i = 1,6
         mm = mm +  mat(i)
      enddo
      res = res + dreal(mm)**2 + dimag(mm)**2

      if (lnlo) then
         do i=1,6
            mv = mv + matv(i)           !virtual not prop. to born
         enddo
      endif

      if (nlo.gt.0) then
         mv = als(1,1)*c2o4pi*(mv + mm*cvirtVV)
      else if(nlo.lt.0) then
         mv = als(1,1)*c2o4pi*mv
      else
         mv = czero
      endif

      resv = resv + 2d0*dreal( mm*dconjg(mv)) !+ dreal(mv)**2 + dimag(mv)**2  ! M_V**2 zero

      if (nlo.ge.0) then
         res = (res+resv)*3d0              ! 3 is the color sum factor
      else
         res = resv*3d0                    ! 3 is the color sum factor
      endif
      qqbarWWW = res

      return
      end

c*******************************************************************************
c  end of subroutine qqWWWp
c*******************************************************************************


c****************subroutine QQWWWm*********************************************
C
C  QQWWWm calculates the matrix elements**2 for
C        q1bar q2    ---->   W- W+ W-,
C         W- ---> f1-bar f2, W+ ---> f3-bar f4, W- ---> f5-bar f6
C
C  and crossing related processes. Pauli interference terms for
c  identical fermions are neglected.
c
C  This code is modified to allow for virtual corrections, more precisely
C  the interference of Born with the finite part of virtual diagrams
C  for 
c
c  INPUT:  see QQWWWp
c*******************************************************************************

      subroutine qqWWWm(pbar,sign, nlo, N, qqbarWWW,NW)

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/scales.inc"

c     electroweak couplings are taken from KOPPLN

      double precision  clr, xm2, xmg, b, v, a
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     1                  V(4,5),A(4,5)
      double precision  pbar(0:3,2+n_v)                         ! input momenta
      double precision  qqbarWWW, res,resv                      ! output M^2 of process
      double precision  p(0:3,2+n_v), pwp1wm1(0:4),pwp1wm2(0:4),! different momenta
     1                  fqp1(0:4,2), fqm1(0:4,2),fqm2(0:4,2)    ! and momentum combinations
      integer  sign(2+n_v), nlo, mu, i, N
      integer NW
      double complex psiwp1(2,-1:1,2),                          ! q+W spinors, e.g. <q~W+(1)|,
     1               psiwm1(2,-1:1,2), psiwm2(2,-1:1,2),        ! and currents
     2               jwp1(0:5,-1:1,2),
     3               jwm1(0:5,-1:1,2), jwm2(0:5,-1:1,2)
      double complex temp, temp1
      double complex mat(6), matv(6), mborn
      double complex madgr(6)
      double complex mm, mv
      double complex m3line1, m2linea1,m2linez1, m1line1
      double complex psi(2,-1:1,2), jqqbar(0:5,-1:1)            ! spinors <q~| and |q>, current <q~|gamma_mu|q>

      double complex m3line2, m2linez2                ! matrix variables for virtual corrections

      double complex dotcc, s1c, dotrc                          ! "dotproducts"
      double precision dotrr, psumsq
      external dotcc, s1c, dotrr, dotrc, psumsq

      double precision rzero,rone,epsilon
      double complex czero, cone
      parameter (rzero =  0d0 ,      rone =  1d0, epsilon = 1.0d-1)
      parameter (czero = (0d0,0d0) , cone = (1d0,0d0))

      double complex zm2i(2:3)                                  !variable for q_mu q_nu terms in the
      save zm2i                                                 !propagators of the massive vector bosons

c variables for virtual corrections

      double precision scale
      double precision c2,c2o4pi
      parameter (c2=4d0/3d0, c2o4pi=c2/4d0/pi)
      logical lnlo, lbox, lpq, lpt, lpent, linit
      data linit /.true./
      save linit
      double complex j5wmwp(1:6)
      double complex jbwmwp(1:6)

      double complex cte
      double complex resultgauge1(3), resultgauge3(3)
      double complex resultgaugeb1(3),resultgaugeb3(3)
      double precision p1(0:3),p2(0:3)
      double precision minus_qm1(0:3),minus_qm2(0:3),minus_q1p(0:3)
      double precision minus_qm1p(0:3), minus_qm2p(0:3)
      double complex result(8),result_born(8)
      double complex gaugebox1(2),gaugebox1_born(2), boxtest(2), treetest(2)
      double complex gaugebox2(2), gaugebox2_born(2)
      double complex gaugebox3(2), gaugebox3_born(2)
      double complex gaugebox4(2), gaugebox4_born(2)
      integer wardidtest1, wardidtest2, wardidtest3, wardidtest4

c variables for test purposes


      logical ldebugm, ldebugs,pentshift
      parameter( ldebugm = .false., ldebugs=.false., pentshift = .true.)

c---------------------------- BEGIN CODE --------------------------------------

      lnlo = NLO.ne.0       ! include some virtual stuff if T
      lbox = (NLO.eq.-4) .or. (NLO.eq.1) .or. (NLO.eq.5)   ! call boxline if T
      lpt = (NLO.eq.1) .or. (NLO.eq.-5)   ! true pentagon contribution
      lpq = (NLO.eq.5) .or. (NLO.eq.-4)   ! q^mu/m_V terms of Pentagon
      lpent = lpt .or. lpq

c define flavors of external quarks for the 2 subprocesses

      if (linit) then
         linit = .false.
         zm2i(2) = 1d0/dcmplx(xm2(2),-xmg(2))
         zm2i(3) = 1d0/dcmplx(xm2(3),-xmg(3))
      endif

      do i = 1,6
         mat(i)  = czero
         matv(i) = czero
      enddo

c define the internal momenta
      do mu = 0,3
         do i = 1,2+n_v
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo
         pwp1wm1(mu) = q1p(mu,N) + qm1(mu,N)
         pwp1wm2(mu) = q1p(mu,N) + qm2(mu,N)
      enddo
      pwp1wm1(4) = dotrr(pwp1wm1(0),pwp1wm1(0))
      pwp1wm2(4) = dotrr(pwp1wm2(0),pwp1wm2(0))

c------------------------------------------------------------------------------

c get the external quark spinors (including factor sqrt(2E) )

      call psi0m(2,pbar(0,1),sign(1),psi)

c get other bra and ket vectors

      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,q1p(0,N),w1p(1,N),psiwp1(1,-1,2),fqp1(0,2))     ! <q~W+(1)|
      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,qm1(0,N),wm1(1,N),psiwm1(1,-1,2),fqm1(0,2))     ! <q~W-(1)|
      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,qm2(0,N),wm2(1,N),psiwm2(1,-1,2),fqm2(0,2))     ! <q~W-(2)|

      call ket2c(psi(1,-1,1),.true.,p(0,1),-1,q1p(0,N),w1p(1,N),psiwp1(1,-1,1),fqp1(0,1))     ! |W+(1)q>
      call ket2c(psi(1,-1,1),.true.,p(0,1),-1,qm1(0,N),wm1(1,N),psiwm1(1,-1,1),fqm1(0,1))     ! |W-(1)q>
      call ket2c(psi(1,-1,1),.true.,p(0,1),-1,qm2(0,N),wm2(1,N),psiwm2(1,-1,1),fqm2(0,1))     ! |W-(2)q>

c get current for contraction with W->WWW leptonic tensor

      call curr6(-1,psi(1,-1,2),p(0,2),psi(1,-1,1),p(0,1), jqqbar(0,-1))                      !<q~|gamma_mu|q>

c get currents for contraction with A/Z->WW leptonic tensor

      call curr6(-1,psi(1,-1,2),p(0,2),psiwp1(1,-1,1),fqp1(0,1),jwp1(0,-1,1))                 !<q~|gamma_mu|W+(1)q>
      call curr6(-1,psi(1,-1,2),p(0,2),psiwm1(1,-1,1),fqm1(0,1),jwm1(0,-1,1))                 !<q~|gamma_mu|W-(1)q>
      call curr6(-1,psi(1,-1,2),p(0,2),psiwm2(1,-1,1),fqm2(0,1),jwm2(0,-1,1))                 !<q~|gamma_mu|W-(2)q>

      call curr6(-1,psiwp1(1,-1,2),fqp1(0,2),psi(1,-1,1),p(0,1),jwp1(0,-1,2))                 !<q~W+(1)|gamma_mu|q>
      call curr6(-1,psiwm1(1,-1,2),fqm1(0,2),psi(1,-1,1),p(0,1),jwm1(0,-1,2))                 !<q~W-(1)|gamma_mu|q>
      call curr6(-1,psiwm2(1,-1,2),fqm2(0,2),psi(1,-1,1),p(0,1),jwm2(0,-1,2))                 !<q~W-(2)|gamma_mu|q>

C------------------------------------------------------------------------------
C CALCULATION OF DU~ --> W-W+W- PROCESSES (TREE LEVEL)

C U~ D --> W-W+W- --> LEPTONS
      m3line1 = -s1c(psiwm1(1,-1,2),w1p(1,N),.true.,-1,psiwm2(1,-1,1))  !W-(2)W+(1)W-(1)
     1          -s1c(psiwm2(1,-1,2),w1p(1,N),.true.,-1,psiwm1(1,-1,1))  !W-(1)W+(1)W-(2)
      m3line1 = m3line1*clr(3,3,-1)**3
      MAT(1)  = M3LINE1                                               !M(DU~->W-W-W+)

C------------------------------------------------------------------------------
C CALCULATION OF DU~ --> W-V --> LEPTONS PROCESSES (TREE LEVEL)

c     DU~ --> ZW- --> LEPTONS
      m2linez1 = -dotcc(WWW_zww(0,N,3),jwm1(0,-1,2))*clr(4,2,-1)  !ZW
     1           -dotcc(WWW_zww(0,N,1),jwm2(0,-1,2))*clr(4,2,-1)  !ZW
      MAT(2)   = M2LINEZ1*CLR(3,3,-1)                      !M(DU~->ZW- ->L)

c    additional contributions due to q_mu q_nu/ m^2 terms in the propagator
      m2linez1 = zm2i(2)*dotrc(pwp1wm2(0),WWW_zww(0,N,3))*dotrc(pwp1wm2(0),jwm1(0,-1,2))
     &         + zm2i(2)*dotrc(pwp1wm1(0),WWW_zww(0,N,1))*dotrc(pwp1wm1(0),jwm2(0,-1,2))
      MAT(2)   = MAT(2) + M2LINEZ1*CLR(4,2,-1)*clr(3,3,-1)

c     DU~ --> AW- --> LEPTON
      m2linea1 = -dotcc(WWW_aww(0,N,3),jwm1(0,-1,2))*clr(4,1,-1)  !AW
     1           -dotcc(WWW_aww(0,N,1),jwm2(0,-1,2))*clr(4,1,-1)  !AW
      MAT(3)   = M2LINEA1*CLR(3,3,-1)                      !M(DU~->AW- ->L)

c     DU~ --> W-Z --> LEPTONS
      m2linez1 = -dotcc(WWW_zww(0,N,3),jwm1(0,-1,1))*clr(3,2,-1)  !WZ
     1           -dotcc(WWW_zww(0,N,1),jwm2(0,-1,1))*clr(3,2,-1)  !WZ
      MAT(4)   = M2LINEZ1*CLR(3,3,-1)                      !M(DU~->W-Z ->L)

c    additional contributions due to q_mu q_nu/ m^2 terms in the propagator
      m2linez1 = zm2i(2)*dotrc(pwp1wm2(0),WWW_zww(0,N,3))*dotrc(pwp1wm2(0),jwm1(0,-1,1))
     &         + zm2i(2)*dotrc(pwp1wm1(0),WWW_zww(0,N,1))*dotrc(pwp1wm1(0),jwm2(0,-1,1))
      MAT(4)   = MAT(4) + M2LINEZ1*CLR(3,2,-1)*clr(3,3,-1)

c     DU~ --> W-A --> LEPTONS
      m2linea1 = -dotcc(WWW_aww(0,N,3),jwm1(0,-1,1))*clr(3,1,-1)  !WA
     1           -dotcc(WWW_aww(0,N,1),jwm2(0,-1,1))*clr(3,1,-1)  !WA
      MAT(5)   = M2LINEA1*CLR(3,3,-1)                      !M(DU~->W-A ->L)

C------------------------------------------------------------------------------
C CALCULATION OF DU~ --> W- --> LEPTONS PROCESSES (TREE LEVEL)

c     DU~  --> W- --> LEPTONS
      m1line1 = -dotcc(wmtowww(0,N),jqqbar(0,-1))
      m1line1 = m1line1*clr(3,3,-1)
      MAT(6)  = -M1LINE1                                       !M(DU~->W- ->L)

c**************************************************************************
c
c                   PENTAGON AND BOX CONTRIBUTIONS
c
c**************************************************************************

#ifdef WITH_NLO
      if(lpent) then
         scale = -2.0d0*dotrr(p(0,1),p(0,2))

c calculate born matrix element and store it in mborn for gauge check:
         mborn=mat(1)+mat(2)+mat(3)+mat(4)+mat(5)+mat(6)

         do mu=0,3
            p1(mu) = pbar(mu,1)
            p2(mu) = pbar(mu,2)
            minus_qm1(mu) = -qm1(mu,1)
            minus_qm2(mu) = -qm2(mu,1)
            minus_q1p(mu) = -p1(mu) -p2(mu) -minus_qm1(mu) -minus_qm2(mu)
         enddo

         cte = (3*pi*(0d0,1d0)-4.0d0*pi**2/3.0d0+7.0d0)

         if(nlo.eq.1) then
               call penlinemmm(p1(0),minus_qm1(0),minus_q1p(0),minus_qm2(0),p2(0),
     1            psi(1,-1,2),psi(1,-1,1),wm1(1,1),w1p(1,1),wm2(1,1),              ! W+(1) W- W+(2)
     2            -1,scale,4,1,resultgauge3(1),j5wmwp(3),resultgaugeb3(1),jbwmwp(3))
               call penlinemmm(p1(0),minus_qm2(0),minus_q1p(0),minus_qm1(0),p2(0),
     1            psi(1,-1,2),psi(1,-1,1),wm2(1,1),w1p(1,1),wm1(1,1),              ! W+(2) W- W+(1)
     2            -1,scale,4,1,resultgauge1(1),j5wmwp(1),resultgaugeb1(1),jbwmwp(1))
            call Wardtest_Pent(psi, p1, p2, minus_qm1, minus_q1p, minus_qm2,
     1               wm1(1,1), w1p(1,1), wm2(1,1), scale, als(1,1),
     2               resultgauge3(1), mborn, wardidtest1)

            call Wardtest_Pent(psi, p1, p2, minus_qm2, minus_q1p, minus_qm1,
     1               wm2(1,1), w1p(1,1), wm1(1,1), scale, als(1,1),
     2               resultgauge1(1), mborn, wardidtest2)

            if (wardidtest1.eq.1 .or. wardidtest2.eq.1) then
               j5wmwp(1) = (0.0d0, 0.0d0)
               j5wmwp(3) = (0.0d0, 0.0d0)
            else
               j5wmwp(1)=  j5wmwp(1)+cte*jbwmwp(1)
               j5wmwp(3)=  j5wmwp(3)+cte*jbwmwp(3)
            endif

         elseif(nlo.eq.-5) then
            call penlinemmm(p1(0),minus_qm1(0),minus_q1p(0),minus_qm2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wmm1(1,1),w1pm(1,1),wmm2(1,1),              ! W+(1) W- W+(2)
     2         -1,scale,4,1,resultgauge3(1),j5wmwp(3),resultgaugeb3(1),jbwmwp(3))
            call penlinemmm(p1(0),minus_qm2(0),minus_q1p(0),minus_qm1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wmm2(1,1),w1pm(1,1),wmm1(1,1),              ! W+(2) W- W+(1)
     2         -1,scale,4,1,resultgauge1(1),j5wmwp(1),resultgaugeb1(1),jbwmwp(1))

            call Wardtest_Pent(psi, p1, p2, minus_qm1, minus_q1p, minus_qm2,
     1               wmm1(1,1), w1pm(1,1), wmm2(1,1), scale, als(1,1),
     2               resultgauge3(1), mborn, wardidtest1)

            call Wardtest_Pent(psi, p1, p2, minus_qm2, minus_q1p, minus_qm1,
     1               wmm2(1,1), w1pm(1,1), wmm1(1,1), scale, als(1,1),
     2               resultgauge1(1), mborn, wardidtest2)

            if (wardidtest1.eq.1 .or. wardidtest2.eq.1) then
               j5wmwp(1) = (0.0d0, 0.0d0)
               j5wmwp(3) = (0.0d0, 0.0d0)
            else
               j5wmwp(1)=  j5wmwp(1)+cte*jbwmwp(1)
               j5wmwp(3)=  j5wmwp(3)+cte*jbwmwp(3)
            endif

         else
c pentagon contributions shifted to boxes by shifting polarisation vectors
            call pent_to_box(p1, minus_qm1, minus_q1p, minus_qm2, p2,
     1         psi(1,-1,2), psi(1,-1,1), wm1(1,1), w1p(1,1),wm2(1,1),
     2         xm1, x1p, xmm2, -1, scale, j5wmwp(3), jbwmwp(3))
            call pent_to_box(p1, minus_qm2, minus_q1p, minus_qm1, p2,
     1         psi(1,-1,2), psi(1,-1,1), wm2(1,1), w1p(1,1),wm1(1,1),
     2         xmm2, x1p, xm1, -1, scale, j5wmwp(1), jbwmwp(1))

            j5wmwp(1)=  j5wmwp(1)+cte*jbwmwp(1)
            j5wmwp(3)=  j5wmwp(3)+cte*jbwmwp(3)

         endif       ! end nlo id selection


         m3line2 = j5wmwp(1) + j5wmwp(3)
         m3line2 = m3line2*clr(3,3,-1)**3
         MATV(1) = M3LINE2

      endif          ! end lpent


      if(lbox) then
         scale = -2.0d0*dotrr(p(0,1),p(0,2))

c calculate born matrix element and store it in mborn for gauge check:
         mborn=mat(1)+mat(2)+mat(3)+mat(4)+mat(5)+mat(6)

         do mu=0,3
            p1(mu) = pbar(mu,1)
            p2(mu) = pbar(mu,2)
            minus_qm1(mu) = -qm1(mu,1)
            minus_qm2(mu) = -qm2(mu,1)
            minus_qm2p(mu) = -p1(mu) -p2(mu) +qm1(mu,1)
            minus_qm1p(mu) = -p1(mu) -p2(mu) +qm2(mu,1)
         enddo

         cte = (3*pi*(0d0,1d0)-4.0d0*pi**2/3.0d0+7.0d0)

         call boxlinemm(p1,minus_qm2p,minus_qm1,p2,
     1    psi(1,-1,2),psi(1,-1,1),WWW_zww(0,1,3),wm1(1,1),
     2    -1,scale,3,1,gaugebox1,result(1),gaugebox1_born,
     3    result_born(1))
         call boxlinemm(p1,minus_qm2p,minus_qm1,p2,
     1    psi(1,-1,2),psi(1,-1,1),WWW_aww(0,1,3),wm1(1,1),
     2    -1,scale,9,-1,boxtest,result(2),treetest,
     3    result_born(2))

         call boxlinemm(p1,minus_qm1p,minus_qm2,p2,
     1    psi(1,-1,2),psi(1,-1,1),WWW_zww(0,1,1),wm2(1,1),
     2    -1,scale,3,1,gaugebox2,result(3),gaugebox2_born,
     3    result_born(3))
         call boxlinemm(p1,minus_qm1p,minus_qm2,p2,
     1    psi(1,-1,2),psi(1,-1,1),WWW_aww(0,1,1),wm2(1,1),
     2    -1,scale,9,-1,boxtest,result(4),treetest,
     3    result_born(4))

         call boxlinemm(p1,minus_qm1,minus_qm2p,p2,
     1    psi(1,-1,2),psi(1,-1,1),wm1(1,1),WWW_zww(0,1,3),
     2    -1,scale,3,1,gaugebox3,result(5),gaugebox3_born,
     3    result_born(5))
         call boxlinemm(p1,minus_qm1,minus_qm2p,p2,
     1    psi(1,-1,2),psi(1,-1,1),wm1(1,1),WWW_aww(0,1,3),
     2    -1,scale,9,-1,boxtest,result(6),treetest,
     3    result_born(6))

         call boxlinemm(p1,minus_qm2,minus_qm1p,p2,
     1    psi(1,-1,2),psi(1,-1,1),wm2(1,1),WWW_zww(0,1,1),
     2    -1,scale,3,1,gaugebox4,result(7),gaugebox4_born,
     3    result_born(7))
         call boxlinemm(p1,minus_qm2,minus_qm1p,p2,
     1    psi(1,-1,2),psi(1,-1,1),wm2(1,1),WWW_aww(0,1,1),
     2    -1,scale,9,-1,boxtest,result(8),treetest,
     3    result_born(8))

         call Wardtest_Box(WWW_zww(0,1,3), pwp1wm2(0), wm1(1,1), qm1(0,1),
     1            gaugebox1, gaugebox1_born, mborn, wardidtest1, als(1,1))

         call Wardtest_Box(WWW_zww(0,1,1), pwp1wm1(0), wm2(1,1), qm2(0,1),
     1            gaugebox2, gaugebox2_born, mborn, wardidtest2, als(1,1))

         call Wardtest_Box(wm1(1,1), qm1(0,1), WWW_zww(0,1,3), pwp1wm2(0),
     1            gaugebox3, gaugebox3_born, mborn, wardidtest3, als(1,1))

         call Wardtest_Box(wm2(1,1), qm2(0,1), WWW_zww(0,1,1), pwp1wm1(0),
     1            gaugebox4, gaugebox4_born, mborn, wardidtest4, als(1,1))

          result(1) = -(result(1)+cte*result_born(1))
          result(2) = -(result(2)+cte*result_born(2))
          result(3) = -(result(3)+cte*result_born(3))
          result(4) = -(result(4)+cte*result_born(4))
          result(5) = -(result(5)+cte*result_born(5))
          result(6) = -(result(6)+cte*result_born(6))
          result(7) = -(result(7)+cte*result_born(7))
          result(8) = -(result(8)+cte*result_born(8))

         if (wardidtest1.eq.0) then
            m2linez2 = result(1)*clr(3,2,-1)
            MATV(2)  = M2LINEZ2*CLR(3,3,-1)
            m2linez2 = result(2)*clr(3,1,-1)
            MATV(3)  = M2LINEZ2*CLR(3,3,-1)
         endif

         if (wardidtest2.eq.0) then
            m2linez2 = result(3)*clr(3,2,-1)
            MATV(2)  = M2LINEZ2*CLR(3,3,-1)
            m2linez2 = result(4)*clr(3,1,-1)
            MATV(3)  = M2LINEZ2*CLR(3,3,-1)
         endif

         if (wardidtest3.eq.0) then
            m2linez2 = result(5)*clr(4,2,-1)
            MATV(4)  = M2LINEZ2*CLR(3,3,-1)
            m2linez2 = result(6)*clr(4,1,-1)
            MATV(5)  = M2LINEZ2*CLR(3,3,-1)
         endif

         if (wardidtest4.eq.0) then
            m2linez2 = result(7)*clr(4,2,-1)
            MATV(4)  = M2LINEZ2*CLR(3,3,-1)
            m2linez2 = result(8)*clr(4,1,-1)
            MATV(5)  = M2LINEZ2*CLR(3,3,-1)
         endif

      endif

      MATV(6) = CZERO
#endif
c------------------------------------------------------------------------------

c comparison with madgraph

#ifdef WITH_MADGRAPH
      if (ldebugm) then
         do i = 1,6
            madgr(i) = czero
         enddo

         do i =1,6
            madgr(i) = WWWtamp(i,-1,2)
         enddo

c print results from checks:

         print*,'---------------------------------------------------'
         do i = 1,6
    	    print*,'ratio for i = ',i,':',
     &	       '(', dreal(mat(i))/dreal(madgr(i)) ,
     &              dimag(mat(i))/dimag(madgr(i)), ')'
            print*
            print*,'mat(',i,') = ',mat(i)
            temp = dreal(mat(i))**2+dimag(mat(i))**2
            print*,temp
            print*,'madgr(',i,') = ',madgr(i)
            temp1 = dreal(madgr(i))**2+dimag(madgr(i))**2
            print*,temp1
            print*,'temp/temp1 = ',temp/temp1
            print*,' '
         enddo
      endif ! end of debugging block
#endif

c------------------------------------------------------------------------------
c------------------------------------------------------------------------------

c sum the graphs, square them and map them onto qqbarWWW

      res = 0.0d0
      resv = 0.0d0
      mm = czero
      mv = czero
      do i = 1,6
         mm = mm +  mat(i)
      enddo
      res = res + dreal(mm)**2 + dimag(mm)**2

      if (lnlo) then
         do i=1,6
            mv = mv + matv(i)           !virtual not prop. to born
         enddo
      endif

      if (nlo.gt.0) then
         mv = als(1,1)*c2o4pi*(mv + mm*cvirtVV)
      else if(nlo.lt.0) then
         mv = als(1,1)*c2o4pi*mv
      else
         mv = czero
      endif
      resv = resv + 2d0*dreal( mm*dconjg(mv))

      if (nlo.ge.0) then
         res = (res+resv)*3d0              ! 3 is the color sum factor
      else
         res = resv*3d0                    ! 3 is the color sum factor
      endif

      qqbarWWW = res

      return
      end

c*******************************************************************************
c  end of subroutine qqWWWm
c*******************************************************************************

