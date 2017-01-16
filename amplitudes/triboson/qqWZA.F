c****************subroutine QQWZA *************************************
c
c       Michael Rauch, <rauch@particle.uni-karlsruhe.de>
c       Initial version:  2009, January 
c       Last modified: 2009, July
c       Adapted from qqVVA.f by G. Bozzi
C
C  QQWZA calculates the matrix elements**2 for 
C     q1 q2    ---->   W+ Z A   
C     W+ ---> f1-bar f2, Z ---> f3-bar f4
C
C  and crossing related processes. Pauli interference terms for
c  identical fermions are neglected.
c
C  This code is modified to allow for virtual corrections, more precisely
C  the interference of Born with the finite part of virtual diagrams
C  for 
c
c  INPUT:  NLO = 1       return udwza = |M_born|^2 + 2Re(M_Born^* M_virt)
c          NLO = 0       return udwza = |M_born|^2   etc.
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

      SUBROUTINE qqWZA(pbar,sign, nlo, N, udwza, bos, ps_number)

      IMPLICIT NONE

#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/scales.inc"


c     electroweak couplings are taken from KOPPLN

      DOUBLE PRECISION  clr, xm2, xmg, b, v, a
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     1                  V(4,5),A(4,5)
      DOUBLE PRECISION  pbar(0:3,2+n_v)
      DOUBLE PRECISION  udwza, res(2),resv(2)
      DOUBLE PRECISION  p(0:3,2+n_v), pwpz(0:4), pwpa(0:4), pzpa(0:4),
     1                  fqp(0:4,2), fqz(0:4,2),
     2                  fqa(0:4,2)

      INTEGER  sign(2+n_v), nlo, mu, i, j, isig, N, ps_number
      INTEGER  ifl(2), bos
      DOUBLE COMPLEX psiwp(2,-1:1,2), psia(2,-1:1,2),
     1               psiaz(2,-1:1,2), psiz(2,-1:1,2), 
     2               jwp(0:5,-1:1,2)

      DOUBLE COMPLEX j5pa(3,2),j5ap(3,2), jbpa(3,2), jbap(3,2)
      DOUBLE COMPLEX j4pa(2),jdpa(2), j4ap(2),jdap(2)
      DOUBLE COMPLEX j4pz(2),jdpz(2), j4zp(2),jdzp(2)
      DOUBLE COMPLEX j4paz(1),jdpaz(1), j4azp(1), jdazp(1)

      DOUBLE COMPLEX temp, temp1
      DOUBLE COMPLEX mat(1,-1:1,3), matv(1,-1:1,3), mborn
      DOUBLE COMPLEX madgr(3,-1:1,2)
      DOUBLE COMPLEX mm(2,-1:1), mv(2,-1:1), zero
      DOUBLE COMPLEX mw, m1, m2, m3, m4, mm1, mm2, mm3, mm4
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
      logical wardresult 
      data linit /.true./
      save linit
      parameter(ldebugm = .false., ldebug=.false.)

c variables for implementation of the pentagon routine

      double precision minus_qp(0:3),minus_qz(0:3),minus_qa(0:3)
      double precision minus_qza(0:3)
      double precision minus_qpz(0:3), minus_qpa(0:3)
      double precision p1(0:3), p2(0:3) 
      double complex pentcontr1(3), pentcontr2(3), pentcontr3(3)
      double complex pentcontr4(3), pentcontr5(3), pentcontr6(3)
      double complex pentcontr1b(3), pentcontr2b(3), pentcontr3b(3)
      double complex pentcontr4b(3), pentcontr5b(3), pentcontr6b(3)
      double complex dummy, dummyb

c variables for implementation of the box routine
      double complex boxcontr1(2)
      double complex boxcontr1b(2)

      double complex cte

c variables for test purposes


c---------------------------- BEGIN CODE -----------------------

      lnlo = NLO.ne.0       ! include some virtual stuff if T
      lbox = (NLO.eq.-4) .or. (NLO.eq.1) .or. (NLO.eq.5)    ! call boxline if T
      lpt = (NLO.eq.1) .or. (NLO.eq.-5)   ! true pentagon contribution
      lpq = (NLO.eq.5) .or. (NLO.eq.-4)   ! q^mu/m_V terms of Pentagon
      lpent = lpt .or. lpq
c define flavors of external quarks for the subprocess

      if (linit) then
         linit = .false.
         if (bos.eq.321) then
           ifl(1) = 3             ! udwza: u = 3; d = 4
           ifl(2) = 4
         else if (bos.eq.421) then
           ifl(1) = 4             
           ifl(2) = 3
         else 
           print *, "Wrong boson ids:", bos
           stop
         endif
         zm2i(2) = 1d0/dcmplx(xm2(2),-xmg(2))
         zm2i(3) = 1d0/dcmplx(xm2(3),-xmg(3))
      endif

c only left-handed quarks can appear due to W and m_q=0
      isig = -1

c       do isig = -1,-1,2
        do i = 1,3
              mat(1,isig,i) = (0.0d0,0.0d0)
              matv(1,isig,i) = (0.0d0,0.0d0)
        enddo
c      enddo

c define the internal momenta

      do mu = 0,3
         do i = 1,2+n_v
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo
         pwpz(mu) = qp(mu,N) + qz(mu,N)
         pwpa(mu) = qp(mu,N) + qa(mu,N)
         pzpa(mu) = qz(mu,N) + qa(mu,N)
      enddo
      pwpz(4) = pwpz(0)**2-pwpz(1)**2-pwpz(2)**2-pwpz(3)**2
      pwpa(4) = pwpa(0)**2-pwpa(1)**2-pwpa(2)**2-pwpa(3)**2
      pzpa(4) = pzpa(0)**2-pzpa(1)**2-pzpa(2)**2-pzpa(3)**2

c------------------------------------------------

c get the external quark spinors (including factor sqrt(2E) )

      call psi0m(2,pbar(0,1),sign(1),psi)

c get other bra and ket vectors

      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,qp(0,N),
     &           wp(1,N),psiwp(1,-1,2),fqp(0,2)) ! <d|W+
      call ket2c(psi(1,-1,1),.true.,p(0,1),-1,qp(0,N),
     &           wp(1,N),psiwp(1,-1,1),fqp(0,1)) ! W+|u>

c      do isig =-1,1,2

         call bra2c(psi(1,isig,2),.true.,p(0,2),isig,qa(0,N),
     &                 wa(1,N), psia(1,isig,2),fqa(0,2)) ! <d|A
         call bra2c(psi(1,isig,2),.true.,p(0,2),isig,qz(0,N),
     &                 ztau(1,N), psiz(1,isig,2),fqz(0,2)) ! <d|Z
         call bra2c(psi(1,isig,2),.true.,p(0,2),isig,qz(0,N),
     &                 atau(1,N), psiaz(1,isig,2),fqz(0,2)) ! <d|A(Z)
         call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qa(0,N),
     &                 wa(1,N), psia(1,isig,1),fqa(0,1)) ! A|u>
         call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qz(0,N),
     &                 ztau(1,N), psiz(1,isig,1),fqz(0,1)) ! Z|u>
         call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qz(0,N),
     &                 atau(1,N), psiaz(1,isig,1),fqz(0,1)) ! A(Z)|u>
c      enddo

c***************************************************************
c
c        Tree-level contributions
c
c***************************************************************

c Calculation of processes qq -> W Z A ---------------------

c   u d -> W Z A
      m1=(-s1c(psiwp(1,-1,2),ztau(1,N),.true.,-1,psia(1,-1,1) ) !AZW
     1    -s1c(psiwp(1,-1,2),wa(1,N)  ,.true.,-1,psiz(1,-1,1) ) !ZAW
     A   )*clr(ifl(1),2,-1)*clr(ifl(1),1,-1)+
     2   (-s1c(psiz(1,-1,2) ,wp(1,N)  ,.true.,-1,psia(1,-1,1) ) !AWZ
     B   )*clr(ifl(2),2,-1)*clr(ifl(1),1,-1)+
     3   (-s1c(psia(1,-1,2) ,wp(1,N)  ,.true.,-1,psiz(1,-1,1) ) !ZWA
     C   )*clr(ifl(1),2,-1)*clr(ifl(2),1,-1)+
     4   (-s1c(psiz(1,-1,2) ,wa(1,N)  ,.true.,-1,psiwp(1,-1,1)) !WAZ
     5    -s1c(psia(1,-1,2) ,ztau(1,N),.true.,-1,psiwp(1,-1,1)) !WZA
     D   )*clr(ifl(2),2,-1)*clr(ifl(2),1,-1)

      m2=(-s1c(psiwp(1,-1,2),atau(1,N),.true.,-1,psia(1,-1,1) ) !AzW
     1    -s1c(psiwp(1,-1,2),wa(1,N)  ,.true.,-1,psiaz(1,-1,1)) !zAW
     A   )*clr(ifl(1),1,-1)*clr(ifl(1),1,-1)+
     2   (-s1c(psiaz(1,-1,2),wp(1,N)  ,.true.,-1,psia(1,-1,1) ) !AWz
     B   )*clr(ifl(2),1,-1)*clr(ifl(1),1,-1)+
     3   (-s1c(psia(1,-1,2) ,wp(1,N)  ,.true.,-1,psiaz(1,-1,1)) !zWA
     C   )*clr(ifl(1),1,-1)*clr(ifl(2),1,-1)+
     4   (-s1c(psiaz(1,-1,2),wa(1,N)  ,.true.,-1,psiwp(1,-1,1)) !WAz
     5    -s1c(psia(1,-1,2) ,atau(1,N),.true.,-1,psiwp(1,-1,1)) !WzA
     D   )*clr(ifl(2),1,-1)*clr(ifl(2),1,-1)

      mat(1,-1,1) = (m1 + m2)*clr(3,3,-1)


c Calculation of processes qq -> V V -> W Z A -------------

c qq -> W A/Z (-> W Z/A)

c W A
      call curr6(-1,psia(1,-1,2),fqa(0,2),
     &           psi(1,-1,1),p(0,1), jwp(0,-1,1) ) ! W -> lt
      m1 = -dotcc(wpwz(0,N),jwp(0,-1,1))
c additional contributions due to q_mu q_nu/ m^2 terms in the propagator
      mm1 = zm2i(3)*dotrc(pwpz(0),wpwz(0,N))*dotrc(pwpz(0),jwp(0,-1,1))
      call curr6(-1,psi(1,-1,2),p(0,2),
     &           psiwp(1,-1,1),fqp(0,1), jwp(0,-1,1) ) ! A -> lt
      m1 = m1 -dotcc(aza(0,N),jwp(0,-1,1))
c W Z
ccc same current as before
c      call curr6(-1,psi(1,-1,2),p(0,2),
c     &           psiwp(1,-1,1),fqp(0,1), jwp(0,-1,1) ) ! Z -> lt
      m2 = -dotcc(zaz(0,N),jwp(0,-1,1))
      call curr6(-1,psiz(1,-1,2),fqz(0,2),
     &           psi(1,-1,1),p(0,1), jwp(0,-1,1) ) ! W -> lt
      m2 = m2 -dotcc(wpwa(0,N),jwp(0,-1,1))
      mm2 = zm2i(3)*dotrc(pwpa(0),wpwa(0,N))*dotrc(pwpa(0),jwp(0,-1,1))
c W z
      call curr6(-1,psiaz(1,-1,2),fqz(0,2),
     &           psi(1,-1,1),p(0,1), jwp(0,-1,1) ) ! W -> lt
      m1 = m1 -dotcc(wpwa(0,N),jwp(0,-1,1))
      mm1 = mm1 + zm2i(3)*dotrc(pwpa(0),wpwa(0,N))*dotrc(pwpa(0),jwp(0,-1,1))
c A W
      call curr6(-1,psi(1,-1,2),p(0,2),
     &           psia(1,-1,1),fqa(0,1), jwp(0,-1,1) ) ! W -> lt
      m3 = -dotcc(wpwz(0,N),jwp(0,-1,1))
      mm3 = zm2i(3)*dotrc(pwpz(0),wpwz(0,N))*dotrc(pwpz(0),jwp(0,-1,1))
      call curr6(-1,psiwp(1,-1,2),fqp(0,2),
     &           psi(1,-1,1),p(0,1), jwp(0,-1,1) ) ! A -> lt
      m3 = m3 -dotcc(aza(0,N),jwp(0,-1,1))
c Z W
ccc same current as before
c      call curr6(-1,psiwp(1,-1,2),fqp(0,2),
c     &           psi(1,-1,1),p(0,1), jwp(0,-1,1) ) ! Z -> lt
      m4 = -dotcc(zaz(0,N),jwp(0,-1,1))
      call curr6(-1,psi(1,-1,2),p(0,2),
     &           psiz(1,-1,1),fqz(0,1), jwp(0,-1,1) ) ! W -> lt
      m4 = m4 -dotcc(wpwa(0,N),jwp(0,-1,1))
      mm4 = zm2i(3)*dotrc(pwpa(0),wpwa(0,N))*dotrc(pwpa(0),jwp(0,-1,1))
c z W
      call curr6(-1,psi(1,-1,2),p(0,2),
     &           psiaz(1,-1,1),fqz(0,1), jwp(0,-1,1) ) ! W -> lt
      m3 = m3 -dotcc(wpwa(0,N),jwp(0,-1,1))
      mm3 = mm3 + zm2i(3)*dotrc(pwpa(0),wpwa(0,N))*dotrc(pwpa(0),jwp(0,-1,1))

         mat(1,isig,2) =( (m1+mm1)*clr(ifl(2),1,-1)
     1                   +(m2+mm2)*clr(ifl(2),2,-1)
     2                   +(m3+mm3)*clr(ifl(1),1,-1)
     3                   +(m4+mm4)*clr(ifl(1),2,-1)
     4                  )*clr(3,3,-1)

c Calculation of processes qq -> W -> W Z A ---------------------

      call curr6(-1,psi(1,-1,2),p(0,2),psi(1,-1,1),p(0,1),jqq(0,-1))

c contract with wtowza tensor

         mw = -dotcc(zwwa(0,N),jqq(0,isig))
         mat(1,isig,3)  = mw*clr(3,3,-1)
         matv(1,isig,3) = (0.0d0,0.0d0)
      
c**************************************************************
c
c    Virtual contributions
c
c**************************************************************

#ifdef WITH_NLO
c Calculation of processes qq -> W Z A ---------------------

c calculate pentagon contributions

      if (lpent) then
         scale = -2d0*dotrr(p(0,1),p(0,2))

         mborn=mat(1,-1,1)+mat(1,-1,2)+mat(1,-1,3)

         do mu = 0,3
            p1(mu) = pbar(mu,1)
            p2(mu) = pbar(mu,2)
            minus_qp(mu) = -qp(mu,1)
            minus_qz(mu) = -qz(mu,1)
            minus_qa(mu) = -p1(mu)-p2(mu)+qp(mu,1)+qz(mu,1)
         enddo

         cte=(3*pi*(0d0,1d0)-4.0d0*pi**2/3.0d0+7.0d0)

         if (NLO.eq.1) then ! complete NLO
c ZWA
           call penlinemmm(p1(0),minus_qz(0),minus_qp(0),minus_qa(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),ztau(1,1),wp(1,1),wa(1,1),
     2         -1,scale,4,1,pentcontr1,j5pa(1,1),pentcontr1b,jbpa(1,1))
           call pent_gauge(p1(0),minus_qz(0),minus_qp(0),minus_qa(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),ztau(1,1),wp(1,1),wa(1,1),
     2         -1,scale,pentcontr1,pentcontr1b,wardresult)
c          print *, "Wardresult 1", wardresult
           if (wardresult) then
c zWA
             call penlinemmm(p1(0),minus_qz(0),minus_qp(0),minus_qa(0),p2(0),
     1           psi(1,-1,2),psi(1,-1,1),atau(1,1),wp(1,1),wa(1,1),
     2           -1,scale,9,-1,dummy,j5pa(1,2),dummyb,jbpa(1,2))

c WZA
           call penlinemmm(p1(0),minus_qp(0),minus_qz(0),minus_qa(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wp(1,1),ztau(1,1),wa(1,1),
     2         -1,scale,4,1,pentcontr2,j5pa(2,1),pentcontr2b,jbpa(2,1))
           call pent_gauge(p1(0),minus_qp(0),minus_qz(0),minus_qa(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wp(1,1),ztau(1,1),wa(1,1),
     2         -1,scale,pentcontr2,pentcontr2b,wardresult)
c          print *, "Wardresult 2", wardresult
           if (wardresult) then
c WzA
             call penlinemmm(p1(0),minus_qp(0),minus_qz(0),minus_qa(0),p2(0),
     1           psi(1,-1,2),psi(1,-1,1),wp(1,1),atau(1,1),wa(1,1),
     2           -1,scale,9,-1,dummy,j5pa(2,2),dummyb,jbpa(2,2))

c WAZ
           call penlinemmm(p1(0),minus_qp(0),minus_qa(0),minus_qz(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wp(1,1),wa(1,1),ztau(1,1),
     2         -1,scale,4,1,pentcontr3,j5pa(3,1),pentcontr3b,jbpa(3,1))
           call pent_gauge(p1(0),minus_qp(0),minus_qa(0),minus_qz(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wp(1,1),wa(1,1),ztau(1,1),
     2         -1,scale,pentcontr3,pentcontr3b,wardresult)
c          print *, "Wardresult 3", wardresult
           if (wardresult) then
c WAz
             call penlinemmm(p1(0),minus_qp(0),minus_qa(0),minus_qz(0),p2(0),
     1           psi(1,-1,2),psi(1,-1,1),wp(1,1),wa(1,1),atau(1,1),
     2           -1,scale,9,-1,dummy,j5pa(3,2),dummyb,jbpa(3,2))

c ZAW
           call penlinemmm(p1(0),minus_qz(0),minus_qa(0),minus_qp(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),ztau(1,1),wa(1,1),wp(1,1),
     2         -1,scale,4,1,pentcontr4,j5ap(1,1),pentcontr4b,jbap(1,1))
           call pent_gauge(p1(0),minus_qz(0),minus_qa(0),minus_qp(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),ztau(1,1),wa(1,1),wp(1,1),
     2         -1,scale,pentcontr4,pentcontr4b,wardresult)
c          print *, "Wardresult 4", wardresult
           if (wardresult) then
c zAW
             call penlinemmm(p1(0),minus_qz(0),minus_qa(0),minus_qp(0),p2(0),
     1           psi(1,-1,2),psi(1,-1,1),atau(1,1),wa(1,1),wp(1,1),
     2           -1,scale,9,-1,dummy,j5ap(1,2),dummyb,jbap(1,2))
c AZW
           call penlinemmm(p1(0),minus_qa(0),minus_qz(0),minus_qp(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wa(1,1),ztau(1,1),wp(1,1),
     2         -1,scale,4,1,pentcontr5,j5ap(2,1),pentcontr5b,jbap(2,1))
           call pent_gauge(p1(0),minus_qa(0),minus_qz(0),minus_qp(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wa(1,1),ztau(1,1),wp(1,1),
     2         -1,scale,pentcontr5,pentcontr5b,wardresult)
c          print *, "Wardresult 5", wardresult
           if (wardresult) then
c AzW
             call penlinemmm(p1(0),minus_qa(0),minus_qz(0),minus_qp(0),p2(0),
     1           psi(1,-1,2),psi(1,-1,1),wa(1,1),atau(1,1),wp(1,1),
     2           -1,scale,9,-1,dummy,j5ap(2,2),dummyb,jbap(2,2))
c AWZ
           call penlinemmm(p1(0),minus_qa(0),minus_qp(0),minus_qz(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wa(1,1),wp(1,1),ztau(1,1),
     2         -1,scale,4,1,pentcontr6,j5ap(3,1),pentcontr6b,jbap(3,1))
           call pent_gauge(p1(0),minus_qa(0),minus_qp(0),minus_qz(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wa(1,1),wp(1,1),ztau(1,1),
     2         -1,scale,pentcontr6,pentcontr6b,wardresult)
c          print *, "Wardresult 6", wardresult
           if (wardresult) then
c AWz
             call penlinemmm(p1(0),minus_qa(0),minus_qp(0),minus_qz(0),p2(0),
     1           psi(1,-1,2),psi(1,-1,1),wa(1,1),wp(1,1),atau(1,1),
     2           -1,scale,9,-1,dummy,j5ap(3,2),dummyb,jbap(3,2))

           endif
           endif
           endif
           endif
           endif
           endif

           if (wardresult) then
             do i = 1,3
               do j = 1,2
                 j5pa(i,j)=j5pa(i,j)+cte*jbpa(i,j)
                 j5ap(i,j)=j5ap(i,j)+cte*jbap(i,j)
               enddo
             enddo
             pentcount2m = pentcount2m+1
             pentcount2ps(ps_number) = pentcount2ps(ps_number)+1
           else
             do i = 1,3
               do j = 1,2
                 j5pa(i,j)=(0.0d0, 0.0d0)
                 j5ap(i,j)=(0.0d0, 0.0d0)
               enddo
             enddo
             pentcountm = pentcountm+1
             pentcountps(ps_number) = pentcountps(ps_number)+1
           endif ! wardresult

         elseif (NLO.eq.-5) then ! only true pentagon
c ZWA
           call penlinemmm(p1(0),minus_qz(0),minus_qp(0),minus_qa(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),ztaup(1,1),wpp(1,1),wap(1,1),
     2         -1,scale,4,1,pentcontr1,j5pa(1,1),pentcontr1b,jbpa(1,1))
           call pent_gauge(p1(0),minus_qz(0),minus_qp(0),minus_qa(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),ztaup(1,1),wpp(1,1),wap(1,1),
     2         -1,scale,pentcontr1,pentcontr1b,wardresult)
c          print *, "Wardresult 1", wardresult
           if (wardresult) then
c zWA
             call penlinemmm(p1(0),minus_qz(0),minus_qp(0),minus_qa(0),p2(0),
     1           psi(1,-1,2),psi(1,-1,1),ataup(1,1),wpp(1,1),wap(1,1),
     2           -1,scale,9,-1,dummy,j5pa(1,2),dummyb,jbpa(1,2))

c WZA
           call penlinemmm(p1(0),minus_qp(0),minus_qz(0),minus_qa(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wpp(1,1),ztaup(1,1),wap(1,1),
     2         -1,scale,4,1,pentcontr2,j5pa(2,1),pentcontr2b,jbpa(2,1))
           call pent_gauge(p1(0),minus_qp(0),minus_qz(0),minus_qa(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wpp(1,1),ztaup(1,1),wap(1,1),
     2         -1,scale,pentcontr2,pentcontr2b,wardresult)
c          print *, "Wardresult 2", wardresult
           if (wardresult) then
c WzA
             call penlinemmm(p1(0),minus_qp(0),minus_qz(0),minus_qa(0),p2(0),
     1           psi(1,-1,2),psi(1,-1,1),wpp(1,1),ataup(1,1),wap(1,1),
     2           -1,scale,9,-1,dummy,j5pa(2,2),dummyb,jbpa(2,2))

c WAZ
           call penlinemmm(p1(0),minus_qp(0),minus_qa(0),minus_qz(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wpp(1,1),wap(1,1),ztaup(1,1),
     2         -1,scale,4,1,pentcontr3,j5pa(3,1),pentcontr3b,jbpa(3,1))
           call pent_gauge(p1(0),minus_qp(0),minus_qa(0),minus_qz(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wpp(1,1),wap(1,1),ztaup(1,1),
     2         -1,scale,pentcontr3,pentcontr3b,wardresult)
c          print *, "Wardresult 3", wardresult
           if (wardresult) then
c WAz
             call penlinemmm(p1(0),minus_qp(0),minus_qa(0),minus_qz(0),p2(0),
     1           psi(1,-1,2),psi(1,-1,1),wpp(1,1),wap(1,1),ataup(1,1),
     2           -1,scale,9,-1,dummy,j5pa(3,2),dummyb,jbpa(3,2))

c ZAW
           call penlinemmm(p1(0),minus_qz(0),minus_qa(0),minus_qp(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),ztaup(1,1),wap(1,1),wpp(1,1),
     2         -1,scale,4,1,pentcontr4,j5ap(1,1),pentcontr4b,jbap(1,1))
           call pent_gauge(p1(0),minus_qz(0),minus_qa(0),minus_qp(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),ztaup(1,1),wap(1,1),wpp(1,1),
     2         -1,scale,pentcontr4,pentcontr4b,wardresult)
c          print *, "Wardresult 4", wardresult
           if (wardresult) then
c zAW
             call penlinemmm(p1(0),minus_qz(0),minus_qa(0),minus_qp(0),p2(0),
     1           psi(1,-1,2),psi(1,-1,1),ataup(1,1),wap(1,1),wpp(1,1),
     2           -1,scale,9,-1,dummy,j5ap(1,2),dummyb,jbap(1,2))
c AZW
           call penlinemmm(p1(0),minus_qa(0),minus_qz(0),minus_qp(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wap(1,1),ztaup(1,1),wpp(1,1),
     2         -1,scale,4,1,pentcontr5,j5ap(2,1),pentcontr5b,jbap(2,1))
           call pent_gauge(p1(0),minus_qa(0),minus_qz(0),minus_qp(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wap(1,1),ztaup(1,1),wpp(1,1),
     2         -1,scale,pentcontr5,pentcontr5b,wardresult)
c          print *, "Wardresult 5", wardresult
           if (wardresult) then
c AzW
             call penlinemmm(p1(0),minus_qa(0),minus_qz(0),minus_qp(0),p2(0),
     1           psi(1,-1,2),psi(1,-1,1),wap(1,1),ataup(1,1),wpp(1,1),
     2           -1,scale,9,-1,dummy,j5ap(2,2),dummyb,jbap(2,2))
c AWZ
           call penlinemmm(p1(0),minus_qa(0),minus_qp(0),minus_qz(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wap(1,1),wpp(1,1),ztaup(1,1),
     2         -1,scale,4,1,pentcontr6,j5ap(3,1),pentcontr6b,jbap(3,1))
           call pent_gauge(p1(0),minus_qa(0),minus_qp(0),minus_qz(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wap(1,1),wpp(1,1),ztaup(1,1),
     2         -1,scale,pentcontr6,pentcontr6b,wardresult)
c          print *, "Wardresult 6", wardresult
           if (wardresult) then
c AWz
             call penlinemmm(p1(0),minus_qa(0),minus_qp(0),minus_qz(0),p2(0),
     1           psi(1,-1,2),psi(1,-1,1),wap(1,1),wpp(1,1),ataup(1,1),
     2           -1,scale,9,-1,dummy,j5ap(3,2),dummyb,jbap(3,2))

           endif
           endif
           endif
           endif
           endif
           endif

           if (wardresult) then
             do i = 1,3
               do j = 1,2
                 j5pa(i,j)=j5pa(i,j)+cte*jbpa(i,j)
                 j5ap(i,j)=j5ap(i,j)+cte*jbap(i,j)
               enddo
             enddo
             pentcount2m = pentcount2m+1
             pentcount2ps(ps_number) = pentcount2ps(ps_number)+1
           else
             do i = 1,3
               do j = 1,2
                 j5pa(i,j)=(0.0d0, 0.0d0)
                 j5ap(i,j)=(0.0d0, 0.0d0)
               enddo
             enddo
             pentcountm = pentcountm+1
             pentcountps(ps_number) = pentcountps(ps_number)+1
           endif ! wardresult

         else ! pentagon part ~ boxes
c ZWA
           call pent_to_box(p1(0),minus_qz(0),minus_qp(0),minus_qa(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),ztau(1,1),wp(1,1),wa(1,1),xz,xp,xa,
     2         -1,scale,j5pa(1,1),jbpa(1,1))
c zWA
           call pent_to_box(p1(0),minus_qz(0),minus_qp(0),minus_qa(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau(1,1),wp(1,1),wa(1,1),xaz,xp,xa,
     2         -1,scale,j5pa(1,2),jbpa(1,2))
c WZA
           call pent_to_box(p1(0),minus_qp(0),minus_qz(0),minus_qa(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wp(1,1),ztau(1,1),wa(1,1),xp,xz,xa,
     2         -1,scale,j5pa(2,1),jbpa(2,1))
c WzA
           call pent_to_box(p1(0),minus_qp(0),minus_qz(0),minus_qa(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wp(1,1),atau(1,1),wa(1,1),xp,xaz,xa,
     2         -1,scale,j5pa(2,2),jbpa(2,2))
c WAZ
           call pent_to_box(p1(0),minus_qp(0),minus_qa(0),minus_qz(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wp(1,1),wa(1,1),ztau(1,1),xp,xa,xz,
     2         -1,scale,j5pa(3,1),jbpa(3,1))
c WAz
           call pent_to_box(p1(0),minus_qp(0),minus_qa(0),minus_qz(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wp(1,1),wa(1,1),atau(1,1),xp,xa,xaz,
     2         -1,scale,j5pa(3,2),jbpa(3,2))
c ZAW
           call pent_to_box(p1(0),minus_qz(0),minus_qa(0),minus_qp(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),ztau(1,1),wa(1,1),wp(1,1),xz,xa,xp,
     2         -1,scale,j5ap(1,1),jbap(1,1))
c zAW
           call pent_to_box(p1(0),minus_qz(0),minus_qa(0),minus_qp(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau(1,1),wa(1,1),wp(1,1),xaz,xa,xp,
     2         -1,scale,j5ap(1,2),jbap(1,2))
c AZW
           call pent_to_box(p1(0),minus_qa(0),minus_qz(0),minus_qp(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wa(1,1),ztau(1,1),wp(1,1),xa,xz,xp,
     2         -1,scale,j5ap(2,1),jbap(2,1))
c AzW
           call pent_to_box(p1(0),minus_qa(0),minus_qz(0),minus_qp(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wa(1,1),atau(1,1),wp(1,1),xa,xaz,xp,
     2         -1,scale,j5ap(2,2),jbap(2,2))
c AWZ
           call pent_to_box(p1(0),minus_qa(0),minus_qp(0),minus_qz(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wa(1,1),wp(1,1),ztau(1,1),xa,xp,xz,
     2         -1,scale,j5ap(3,1),jbap(3,1))
c AWz
           call pent_to_box(p1(0),minus_qa(0),minus_qp(0),minus_qz(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wa(1,1),wp(1,1),atau(1,1),xa,xp,xaz,
     2         -1,scale,j5ap(3,2),jbap(3,2))

           do i = 1,3
             do j = 1,2
               j5pa(i,j)=j5pa(i,j)+cte*jbpa(i,j)
               j5ap(i,j)=j5ap(i,j)+cte*jbap(i,j)
             enddo
           enddo

         endif  ! nlo

c u d -> W Z A
         mv1 = j5pa(1,1) * clr(ifl(1),2,-1)*clr(ifl(2),1,-1) !ZWA
     1       + j5ap(3,1) * clr(ifl(2),2,-1)*clr(ifl(1),1,-1) !AWZ
     2       + (j5pa(2,1)+j5pa(3,1)) * clr(ifl(2),2,-1)*clr(ifl(2),1,-1) !WZA+WAZ
     3       + (j5ap(1,1)+j5ap(2,1)) * clr(ifl(1),2,-1)*clr(ifl(1),1,-1) !ZAW+AZW
         mv2 = j5pa(1,2) * clr(ifl(1),1,-1)*clr(ifl(2),1,-1) !zWA
     1       + j5ap(3,2) * clr(ifl(2),1,-1)*clr(ifl(1),1,-1) !AWz
     2       + (j5pa(2,2)+j5pa(3,2)) * clr(ifl(2),1,-1)*clr(ifl(2),1,-1) !WzA+WAz
     3       + (j5ap(1,2)+j5ap(2,2)) * clr(ifl(1),1,-1)*clr(ifl(1),1,-1) !zAW+AzW

         matv(1,-1,1) = (mv1 + mv2) * clr(3,3,-1)

      endif !lpent


c Calculation of processes qq -> V V -> W Z A -------------

      if (lbox) then
         scale = -2d0*dotrr(p(0,1),p(0,2))

c calculate born matrix element and store it in mborn for gauge check
         mborn = mat(1,-1,1)+mat(1,-1,2)+mat(1,-1,3)

         do mu = 0,3
            p1(mu) = pbar(mu,1)
            p2(mu) = pbar(mu,2)
            minus_qp(mu)  = -qp(mu,1)
            minus_qz(mu)  = -qz(mu,1)
            minus_qa(mu)  = -p1(mu)-p2(mu)+qp(mu,1)+qz(mu,1)
            minus_qpz(mu) = -qp(mu,1)-qz(mu,1)
            minus_qpa(mu) = -p1(mu)-p2(mu)+qz(mu,1)
            minus_qza(mu) = -p1(mu)-p2(mu)+qp(mu,1)
         enddo

         cte=(3*pi*(0d0,1d0)-4.0d0*pi**2/3.0d0+7.0d0)

c qq -> W A/Z (-> W Z/A)

c W> A
         call boxlinemm(p1,minus_qpz,minus_qa,p2,
     1    psi(1,-1,2),psi(1,-1,1),wpwz(0,N),wa(1,N),
     2    -1,scale,3,1,boxcontr1,j4pa(1),boxcontr1b,jdpa(1))
         call box_gauge(p1,minus_qpz,minus_qa,p2,
     1    scale,boxcontr1,boxcontr1b,wardresult)
c          print *, "Wardresult 1", wardresult
         if (wardresult) then
c W  A>
         call boxlinemm(p1,minus_qp,minus_qza,p2,
     1    psi(1,-1,2),psi(1,-1,1),wp(1,N),aza(0,N),
     2    -1,scale,3,1,boxcontr1,j4pa(2),boxcontr1b,jdpa(2))
         call box_gauge(p1,minus_qp,minus_qza,p2,
     1    scale,boxcontr1,boxcontr1b,wardresult)
c          print *, "Wardresult 2", wardresult
         if (wardresult) then
c W  Z>
         call boxlinemm(p1,minus_qp,minus_qza,p2,
     1    psi(1,-1,2),psi(1,-1,1),wp(1,N),zaz(0,N),
     2    -1,scale,8,-1,dummy,j4pz(2),dummyb,jdpz(2))
c W> Z
         call boxlinemm(p1,minus_qpa,minus_qz,p2,
     1    psi(1,-1,2),psi(1,-1,1),wpwa(0,N),ztau(1,N),
     2    -1,scale,3,1,boxcontr1,j4pz(1),boxcontr1b,jdpz(1))
         call box_gauge(p1,minus_qpa,minus_qz,p2,
     1    scale,boxcontr1,boxcontr1b,wardresult)
c          print *, "Wardresult 3", wardresult
         if (wardresult) then
c W> z
         call boxlinemm(p1,minus_qpa,minus_qz,p2,
     1    psi(1,-1,2),psi(1,-1,1),wpwa(0,N),atau(1,N),
     2    -1,scale,8,-1,boxcontr1,j4paz(1),boxcontr1b,jdpaz(1))
c A  W>
         call boxlinemm(p1,minus_qa,minus_qpz,p2,
     1    psi(1,-1,2),psi(1,-1,1),wa(1,N),wpwz(0,N),
     2    -1,scale,3,1,boxcontr1,j4ap(2),boxcontr1b,jdap(2))
         call box_gauge(p1,minus_qa,minus_qpz,p2,
     1    scale,boxcontr1,boxcontr1b,wardresult)
c          print *, "Wardresult 4", wardresult
         if (wardresult) then
c A> W
         call boxlinemm(p1,minus_qza,minus_qp,p2,
     1    psi(1,-1,2),psi(1,-1,1),aza(0,N),wp(1,N),
     2    -1,scale,3,1,boxcontr1,j4ap(1),boxcontr1b,jdap(1))
         call box_gauge(p1,minus_qza,minus_qp,p2,
     1    scale,boxcontr1,boxcontr1b,wardresult)
c          print *, "Wardresult 5", wardresult
         if (wardresult) then
c Z> W
         call boxlinemm(p1,minus_qza,minus_qp,p2,
     1    psi(1,-1,2),psi(1,-1,1),zaz(0,N),wp(1,N),
     2    -1,scale,8,-1,boxcontr1,j4zp(1),boxcontr1b,jdzp(1))
c Z  W>
         call boxlinemm(p1,minus_qz,minus_qpa,p2,
     1    psi(1,-1,2),psi(1,-1,1),ztau(1,N),wpwa(0,N),
     2    -1,scale,3,1,boxcontr1,j4zp(2),boxcontr1b,jdzp(2))
         call box_gauge(p1,minus_qz,minus_qpa,p2,
     1    scale,boxcontr1,boxcontr1b,wardresult)
c          print *, "Wardresult 6", wardresult
         if (wardresult) then
c z  W>
         call boxlinemm(p1,minus_qz,minus_qpa,p2,
     1    psi(1,-1,2),psi(1,-1,1),atau(1,N),wpwa(0,N),
     2    -1,scale,8,-1,boxcontr1,j4azp(1),boxcontr1b,jdazp(1))

         endif
         endif
         endif
         endif
         endif
         endif

         if (wardresult) then
           do i=1,2
             j4pa(i) = -(j4pa(i)+cte*jdpa(i))
             j4pz(i) = -(j4pz(i)+cte*jdpz(i))
             j4ap(i) = -(j4ap(i)+cte*jdap(i))
             j4zp(i) = -(j4zp(i)+cte*jdzp(i))
           enddo
           j4paz(1) = -(j4paz(1)+cte*jdpaz(1))
           j4azp(1) = -(j4azp(1)+cte*jdazp(1))

           mv1 = (j4pa(1)+j4pa(2)+j4paz(1))*clr(ifl(2),1,-1)
           mv2 = (j4pz(1)+j4pz(2))         *clr(ifl(2),2,-1)
           mv3 = (j4ap(1)+j4ap(2)+j4azp(1))*clr(ifl(1),1,-1)
           mv4 = (j4zp(1)+j4zp(2))         *clr(ifl(1),2,-1)
           matv(1,-1,2) = (mv1+mv2+mv3+mv4)*clr(3,3,-1)
         else 
           matv(1,-1,2) = (0.0d0, 0.0d0)
         endif

      endif

c Calculation of processes qq -> W -> W Z A ---------------------

      matv(1,isig,3) = (0.0d0,0.0d0)

#endif /* WITH_NLO */
c-----------------------------------------------------
c------------------------------------------------------
#ifdef WITH_MADGRAPH
c comparison with madgraph
      if (ldebugm) then

c         do isig = -1,1,2
            do i = 1,3

                  madgr(i,isig,1) = 0d0

            enddo
c         enddo

c         do isig = -1,1,2
            do i =1,3

                 madgr(i,isig,1) = tamp(i,isig,1)

            enddo
c         enddo

c print results from checks:

         do i=1,3
             print*,'mat(',1,',',isig,',',i,') = ',mat(1,isig,i)
             temp = dreal(mat(1,isig,i))**2+dimag(mat(1,isig,i))**2
             print*,'madgr(',i,',',isig,',',1,') = ',madgr(i,isig,1)
             temp1 = dreal(madgr(i,isig,1))**2+dimag(madgr(i,isig,1))**2
             print*,'temp/temp1 = ',temp/temp1
         enddo

      endif ! end of debugging block
#endif
c------------------------------------------------------


c sum the graphs, square them and map them onto uuww and ddww

       res(1) = 0
       resv(1) = 0
c       do isig = -1,1,2
          mm(1,isig) = 0
          do i = 1,3
             mm(1,isig) = mm(1,isig) +  mat(1,isig,i)
          enddo
          res(1)= res(1) + dreal(mm(1,isig))**2
     &                   + dimag(mm(1,isig))**2

c  add Born type term and multiply by F_q = alphas*C_2/4pi
c  the additional factor for the born term is after adding the 
c  subtraction term and the counter term for the renormalization
c  of the pdfs (see also D(x) in m2s_qqWWZj.f)

#ifdef WITH_NLO
            if (lnlo) then
               mv(1,isig)=matv(1,isig,1)+matv(1,isig,2)+matv(1,isig,3)
            endif
#endif
            if (nlo.gt.0) then
#ifdef WITH_NLO
               mv(1,isig) = als(1,1)*c2o4pi*
     1                  (mv(1,isig) + mm(1,isig)* cvirtVV) 
#endif
            else if(nlo.lt.0) then
#ifdef WITH_NLO
               mv(1,isig) = als(1,1)*c2o4pi*mv(1,isig)
#endif
            else
               mv(1,isig) = (0d0, 0d0)
            endif

           resv(1) = resv(1) + 2d0*dreal( mm(1,isig)   *
     1              conjg( mv(1,isig) ) )

c         enddo

      if (nlo.ge.0) then
         res(1) = (res(1)+resv(1))*3d0      ! 3 is the color sum factor
      else
         res(1) = resv(1)*3d0               ! 3 is the color sum factor
      endif  

      udwza = res(1)

      return
      end


