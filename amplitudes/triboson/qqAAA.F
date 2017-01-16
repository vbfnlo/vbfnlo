c****************subroutine QQAAA *************************************
c
c	Giuseppe Bozzi, <giuseppe@particle.uni-karlsruhe.de>
c	Last modified: Nov. 2010 (giuseppe)
c       Francisco Campanario, <francam@particle.uni-karlsruhe.de>
c	Last modified: Dec. 2010 (paco)
c       Including virtual corrections
C
C  QQAAA calculates the matrix elements**2 for 
C        q1 q1bar    ---->   A A A  
C
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

      subroutine qqAAA(pbar,fsign, nlo, N, uuaaa, ddaaa,ps_number)

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/scales.inc"

c     electroweak couplings are taken from KOPPLN
      double precision  clr, xm2, xmg, b, v, a
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     1                  V(4,5),A(4,5)
      double precision  pbar(0:3,2+n_v)
      double precision  uuaaa, ddaaa, res(2),resv(2)
      double precision  p(0:3,2+n_v)
      double precision  fqa1(0:4,2), fqa2(0:4,2), fqa3(0:4,2)
      integer  fsign(2+n_v), sign(2+n_v), nlo, mu, i, k, isig, N
      integer  ifl(2), ps_number
      double complex psia1(2,-1:1,2),psia2(2,-1:1,2),psia3(2,-1:1,2)


      double complex temp, temp1
      double complex mat(2,-1:1,1:5)
      double complex madgr(-1:1,2)
      double complex mm(2,-1:1), mv(2,-1:1), zero
      double complex ma, mz, m1
      double complex psi(2,-1:1,2), jqq(0:5,-1:1)
      double complex dotcc, s1c, dotrc
      double precision dotrr, psumsq
      external dotcc, s1c, dotrr, dotrc, psumsq
      save ifl
      parameter (zero = (0d0,0d0) )

c  helicity selection
      INTEGER h
      COMMON /hcount / h

c variables for virtual corrections
      double precision scale
      double precision c2,c2o4pi
      parameter (c2=4d0/3d0, c2o4pi=c2/4d0/pi)
      logical lnlo, lbox, lpq, lpt, lpent, ldebugm, linit


c   pentagons 
 
      integer i1,i2
      double complex cte
      double complex pentcontr1(3), pentcontr2(3), pentcontr3(3)
      double complex pentcontr4(3), pentcontr5(3), pentcontr6(3)
      double complex pentcontr1b(3), pentcontr2b(3), pentcontr3b(3)
      double complex pentcontr4b(3), pentcontr5b(3), pentcontr6b(3)
      double complex testg(3)
      Logical wardidtest1,wardidtest2,wardidtest3
      Logical wardidtest4,wardidtest5,wardidtest6


      real*8 p1(0:3),p2(0:3),minus_qa1(0:3),minus_qa2(0:3),minus_qz(0:3)


      double complex j5a1a2(3,-1:1,2),jba1a2(3,-1:1,2)
      double complex j5a2a1(3,-1:1,2),jba2a1(3,-1:1,2)


      double complex mv1, mv2
      double complex  mv1m,mv2m

      double complex  matv(2,-1:1,5),matm(2,-1:1,5)


      data linit /.true./
      save linit
      parameter(ldebugm = .false.)

c---------------------------- BEGIN CODE -----------------------

      lnlo = NLO.ne.0           ! include some virtual stuff if T
      lbox = (NLO.eq.-4) .or. (NLO.eq.1) .or. (NLO.eq.5) ! call boxline if T
      lpt = (NLO.eq.1) .or. (NLO.eq.-5) ! true pentagon contribution
      lpq = (NLO.eq.5) .or. (NLO.eq.-4) ! q^mu/m_V terms of Pentagon
      lpq=.false.
      lpent = lpt .or. lpq

c     define flavors of external quarks for the 2 subprocesses
      if (linit) then
         linit = .false.
         k = 1                  ! uuaaa: u = 3
         ifl(k) = 3
         k = 2                  ! ddaaa: d = 4
         ifl(k) = 4
      endif

      do k = 1,2
         do isig = -1,1,2
            do i2=1,5
            mat(k,isig,i2) = (0.0d0,0.0d0)
            matv(k,isig,i2) = (0.0d0,0.0d0)
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

c get the external quark spinors (including factor sqrt(2E) )
      call psi0m(2,pbar(0,1),sign(1),psi)

c get other bra and ket vectors
      do isig =-1,1,2

         call bra2c(psi(1,isig,2),.true.,p(0,2),isig,qa1(0,N),
     &              atau1(1,N), psia1(1,isig,2),fqa1(0,2))
         call bra2c(psi(1,isig,2),.true.,p(0,2),isig,qa2(0,N),
     &              atau2(1,N), psia2(1,isig,2),fqa2(0,2))
         call bra2c(psi(1,isig,2),.true.,p(0,2),isig,qa3(0,N),
     &              atau3(1,N), psia3(1,isig,2),fqa3(0,2))
         call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qa1(0,N),
     &              atau1(1,N), psia1(1,isig,1),fqa1(0,1))
         call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qa2(0,N),
     &              atau2(1,N), psia2(1,isig,1),fqa2(0,1))
         call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qa3(0,N),
     &              atau3(1,N), psia3(1,isig,1),fqa3(0,1))

      enddo

*******************************************************
C Calculation of born matrix elements
*******************************************************

c q qbar -> A A A
      do isig = -1,1,2
         m1=-s1c(psia3(1,isig,2),atau1(1,N),.true.,isig,psia2(1,isig,1))
     1      -s1c(psia3(1,isig,2),atau2(1,N),.true.,isig,psia1(1,isig,1))
     2      -s1c(psia1(1,isig,2),atau3(1,N),.true.,isig,psia2(1,isig,1))
     3      -s1c(psia1(1,isig,2),atau2(1,N),.true.,isig,psia3(1,isig,1))
     4      -s1c(psia2(1,isig,2),atau3(1,N),.true.,isig,psia1(1,isig,1))
     5      -s1c(psia2(1,isig,2),atau1(1,N),.true.,isig,psia3(1,isig,1))

         mat(1,isig,1) =  m1 * clr(3,1,isig)**3                 ! uu -> AAA
         mat(2,isig,1) =  m1 * clr(4,1,isig)**3                 ! dd -> AAA
      enddo

c qqbar -> A/Z -> Z A A 
      call curr6(1,psi(1,-1,2),p(0,2),psi(1,-1,1),p(0,1),jqq(0,-1))

c contract with atozaa/ztozaa tensor
      do isig = -1,1,2
         ma = -dotcc(aaaa(0,N),jqq(0,isig))
         mz = -dotcc(zaaa(0,N),jqq(0,isig))
         do k = 1,2
            mat(k,isig,4)  = ma*clr(ifl(k),1,isig)
            mat(k,isig,5)  = mz*clr(ifl(k),2,isig)
         enddo
      enddo

*******************************************************
* Calculation of virtual contribution matrix elements
*******************************************************
c#ifdef WITH_NLO
c PACO
#ifdef WITH_NLO
      if (lpent) then
c Calculation of processes qq -> Z/A A A ---------------------
      scale=-2d0*dotrr(p(0,1),p(0,2))
!        Print*, "scale1",scale

            do mu = 0,3
               p1(mu) = pbar(mu,1)
               p2(mu) = pbar(mu,2)
               minus_qa1(mu) = -qa1(mu,1)
               minus_qa2(mu) = -qa2(mu,1)
               minus_qz(mu) = -p1(mu)-p2(mu)+qa1(mu,1)+qa2(mu,1)
            enddo
      
      cte=(3*pi*(0d0,1d0)-4.0d0*pi*pi/3.0d0+7.0d0)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CHANGE Should be -5
        if (NLO.eq.-5) then ! only true pentagon 
c A1 A2 Z	 
            call penlinemmm(p1(0),minus_qa1(0),minus_qa2(0),minus_qz(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau1(1,1),atau2(1,1),atau3(1,1),
     2         -1,scale,4,1,pentcontr1,j5a1a2(3,-1,1),pentcontr1b,jba1a2(3,-1,1))
            call penlinemmm(p1(0),minus_qa1(0),minus_qa2(0),minus_qz(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),atau1(1,1),atau2(1,1),atau3(1,1),
     2         1,scale,9,-1,testg,j5a1a2(3,1,1),testg,jba1a2(3,1,1))
c A1 Z A2
            call penlinemmm(p1(0),minus_qa1(0),minus_qz(0),minus_qa2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau1(1,1),atau3(1,1),atau2(1,1),
     2         -1,scale,4,1,pentcontr2,j5a1a2(2,-1,1),pentcontr2b,jba1a2(2,-1,1))
            call penlinemmm(p1(0),minus_qa1(0),minus_qz(0),minus_qa2(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),atau1(1,1),atau3(1,1),atau2(1,1),
     2         1,scale,9,-1,testg,j5a1a2(2,1,1),testg,jba1a2(2,1,1))
c Z A1 A2
            call penlinemmm(p1(0),minus_qz(0),minus_qa1(0),minus_qa2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau3(1,1),atau1(1,1),atau2(1,1),
     2         -1,scale,4,1,pentcontr3,j5a1a2(1,-1,1),pentcontr3b,jba1a2(1,-1,1))
            call penlinemmm(p1(0),minus_qz(0),minus_qa1(0),minus_qa2(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),atau3(1,1),atau1(1,1),atau2(1,1),
     2         1,scale,9,-1,testg,j5a1a2(1,1,1),testg,jba1a2(1,1,1))
c A2 A1 Z
            call penlinemmm(p1(0),minus_qa2(0),minus_qa1(0),minus_qz(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau2(1,1),atau1(1,1),atau3(1,1),
     2         -1,scale,4,1,pentcontr4,j5a2a1(3,-1,1),pentcontr4b,jba2a1(3,-1,1))
            call penlinemmm(p1(0),minus_qa2(0),minus_qa1(0),minus_qz(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),atau2(1,1),atau1(1,1),atau3(1,1),
     2         1,scale,9,-1,testg,j5a2a1(3,1,1),testg,jba2a1(3,1,1))
c A2 Z A1
            call penlinemmm(p1(0),minus_qa2(0),minus_qz(0),minus_qa1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau2(1,1),atau3(1,1),atau1(1,1),
     2         -1,scale,4,1,pentcontr5,j5a2a1(2,-1,1),pentcontr5b,jba2a1(2,-1,1))
            call penlinemmm(p1(0),minus_qa2(0),minus_qz(0),minus_qa1(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),atau2(1,1),atau3(1,1),atau1(1,1),
     2         1,scale,9,-1,testg,j5a2a1(2,1,1),testg,jba2a1(2,1,1))
c Z A2 A1
            call penlinemmm(p1(0),minus_qz(0),minus_qa2(0),minus_qa1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau3(1,1),atau2(1,1),atau1(1,1),
     2         -1,scale,4,1,pentcontr6,j5a2a1(1,-1,1),pentcontr6b,jba2a1(1,-1,1))
            call penlinemmm(p1(0),minus_qz(0),minus_qa2(0),minus_qa1(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),atau3(1,1),atau2(1,1),atau1(1,1),
     2         1,scale,9,-1,testg,j5a2a1(1,1,1),testg,jba2a1(1,1,1))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c A1 A2 Z
           call pent_gauge(p1(0),minus_qa1(0),minus_qa2(0),minus_qz(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau1(1,1),atau2(1,1),atau3(1,1),
     2         -1,scale,pentcontr1,pentcontr1b,wardidtest1)
c A1 Z A2
            call pent_gauge(p1(0),minus_qa1(0),minus_qz(0),minus_qa2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau1(1,1),atau3(1,1),atau2(1,1),
     2         -1,scale,pentcontr2,pentcontr2b,wardidtest2)
c Z A1 A2 
            call pent_gauge(p1(0),minus_qz(0),minus_qa1(0),minus_qa2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau3(1,1),atau1(1,1),atau2(1,1),
     2         -1,scale,pentcontr3,pentcontr3b,wardidtest3)
c A2 A1 Z
            call pent_gauge(p1(0),minus_qa2(0),minus_qa1(0),minus_qz(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau2(1,1),atau1(1,1),atau3(1,1),
     2         -1,scale,pentcontr4,pentcontr4b,wardidtest4)
c A2 Z A1 
            call pent_gauge(p1(0),minus_qa2(0),minus_qz(0),minus_qa1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau2(1,1),atau3(1,1),atau1(1,1),
     2         -1,scale,pentcontr5,pentcontr5b,wardidtest5)
c Z A2 A1
            call pent_gauge(p1(0),minus_qz(0),minus_qa2(0),minus_qa1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau3(1,1),atau2(1,1),atau1(1,1),
     2         -1,scale,pentcontr6,pentcontr6b,wardidtest6)

  
       wardidtest1=wardidtest1.and.wardidtest2.and.wardidtest3.and.
     1             wardidtest4.and.wardidtest5.and.wardidtest6

 
            if (.not.wardidtest1) then
                do i = 1,3
                  do i1=-1,1,2
                     do i2=1,2
                       j5a1a2(i,i1,i2) = (0.0d0, 0.0d0)
                       j5a2a1(i,i1,i2) = (0.0d0, 0.0d0)
                     enddo
                  enddo
                 enddo
	     pentcountm = pentcountm+1
             pentcountps(ps_number) = pentcountps(ps_number)+1
c               pentcount = pentcount+1
            else
             do i = 1,3
                  do i1=-1,1,2
                     do i2=1,1
                        j5a1a2(i,i1,i2)=j5a1a2(i,i1,i2)+cte*jba1a2(i,i1,i2)
                        j5a2a1(i,i1,i2)=j5a2a1(i,i1,i2)+cte*jba2a1(i,i1,i2)
                     enddo
                 enddo
              enddo
	     pentcount2m = pentcount2m+1
             pentcount2ps(ps_number) = pentcount2ps(ps_number)+1
c               pentcount2 = pentcount2+1
            endif
c   q qbar -> A A W


            do isig=-1,1,2
            mv1=0d0
            mv2=0d0
            mv1m=0d0
            mv2m=0d0
           
           do i1=1,3 
             mv1=mv1+ j5a1a2(i1,isig,1)+ j5a2a1(i1,isig,1)
             mv1m=mv1m+ jba1a2(i1,isig,1)+ jba2a1(i1,isig,1)
           enddo

         matv(1,isig,1) =  mv1 * clr(3,1,isig)**3                 ! uu -> AAA
         matv(2,isig,1) =  mv1 * clr(4,1,isig)**3                 ! dd -> AAA
         matm(1,isig,1) =  mv1m * clr(3,1,isig)**3                 ! uu -> AAA
         matm(2,isig,1) =  mv1m * clr(4,1,isig)**3                 ! dd -> AAA


c$$$         print*, "mv1m",mv1m
c$$$         print*, "mv2m",mv2m
c$$$!         print*, "jba1a2",jba1a2
c$$$!         print*, "jba2a1",jba2a1
c$$$         print*, "mat(1,isig,1)",mat(1,isig,1)
c$$$         print*, "mat(2,isig,1)",mat(2,isig,1)
c$$$         print*, "matm(1,isig,1)",matm(1,isig,1)
c$$$         print*, "matm(2,isig,1)",matm(2,isig,1)
c$$$         print*, "ratio",matm(1,isig,1)/mat(1,isig,1)
c$$$         print*, "ratio",matm(2,isig,1)/mat(2,isig,1)
c$$$         print*, "matv(1,isig,1)",matv(1,isig,1)
c$$$         print*, "matv(2,isig,1)",matv(2,isig,1)
            enddo
c$$$            stop

!         Print*, 'Full Pen matv(k,-1,1)',matv(k,-1,1)
!         Print*, ''
         endif !NLO

      endif  !lpent

#endif /* WITH_NLO */
      
c Calculation of processes qq -> A V -> A A + l+neutriono ****************
#ifdef WITH_NLO
      

      if (lbox) then

          do k = 1,2
              do isig=-1,1,2
                 matv(k,isig,2) = (0.0d0, 0.0d0)
                 matv(k,isig,3) = (0.0d0, 0.0d0)
              enddo
           enddo

      endif !lbox
#endif
c Calculation of processes qq -> Z/A -> A A Z ********************
C Factorize to the born
        do k = 1,2
        do isig = -1,1,2
           matv(k,isig,4) = (0.0d0, 0.0d0)     
           matv(k,isig,5) = (0.0d0, 0.0d0) 
        enddo
        enddo
c#endif /* WITH_NLO */

c**************************************************************
c------------------------------------------------------------
c
c
c comparison with madgraph
c
c
c------------------------------------------------------------
c**************************************************************
c-----------------------------------------------------
c------------------------------------------------------
cc----------------------------------------------------------------------
      if (ldebugm) then
         do isig = -1,1,2
            do k = 1,2
               madgr(isig,k) = tamp(1,isig,k)
            enddo
         enddo

c print results from checks:

         do k = 1,2
            do isig = -1,1,2
              print*,'mat(',k,',',isig,',',1,') = ',mat(k,isig,1)
              temp = dreal(mat(k,isig,1))**2+dimag(mat(k,isig,1))**2
              print*,'madgr(',1,',',isig,',',k,') = ',madgr(isig,k)
              temp1 = dreal(madgr(isig,k))**2+dimag(madgr(isig,k))**2
              print*,'temp/temp1 = ',temp/temp1
            enddo
         enddo
         
      endif                     ! end of debugging block

c------------------------------------------------------

c sum the graphs, square them and map them onto uuaaa and ddaaa

      do k = 1,2
         res(k) = 0
         resv(k) = 0
         do isig = -1,1,2
            mm(k,isig) = 0
            mm(k,isig) = mm(k,isig) +  mat(k,isig,1) +  mat(k,isig,4) +  mat(k,isig,5)
            res(k)= res(k) + dreal(mm(k,isig))**2
     &                     + dimag(mm(k,isig))**2

c  add Born type term and multiply by F_q = alphas*C_2/4pi
c  the additional factor for the born term is after adding the 
c  subtraction term and the counter term for the renormalization
c  of the pdfs (see also D(x) in m2s_qqWWZj.f)

            if (lnlo) then
               mv(k,isig)=matv(k,isig,1)
            endif
            if (nlo.gt.0) then
               mv(k,isig) = als(1,1)*c2o4pi*
     1                  (mv(k,isig) + mm(k,isig)* cvirtVV)
            else if(nlo.lt.0) then
               mv(k,isig) = als(1,1)*c2o4pi*mv(k,isig)
            else
                mv(k,isig) = (0d0, 0d0)
            endif

           resv(k) = resv(k) + 2d0*dreal( mm(k,isig)   *
     1              conjg( mv(k,isig) ) )

         enddo

         if (nlo.ge.0) then
            res(k) = (res(k)+resv(k))*3d0 ! 3 is the color sum factor
         else
            res(k) = resv(k)*3d0 ! 3 is the color sum factor
         endif  
      enddo
      
      uuaaa = res(1)
      ddaaa = res(2)

      return
      end

