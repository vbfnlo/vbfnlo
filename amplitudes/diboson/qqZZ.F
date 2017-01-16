c****************subroutine QQZZ *************************************
c
c	Michael Rauch, <rauch@particle.uni-karlsruhe.de>
c	Last modified: Jun 2011 (rauch)
c       Including virtual corrections
C
C  QQZZ calculates the matrix elements**2 for 
C        q1 q1bar    ---->   Z Z   
C
C
C  and crossing related processes. Pauli interference terms for
c  identical fermions are neglected.
c
C  This code is modified to allow for virtual corrections, more precisely
C  the interference of Born with the finite part of virtual diagrams
C  for 
c
c  INPUT:  NLO = 1       return uuzz = |M_born|^2 + 2Re(M_Born^* M_virt)
c          NLO = 0       return uuzz = |M_born|^2   etc.
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

      subroutine qqzz(pbar,fsign, nlo, N, uuzz, ddzz)

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/scales.inc"

c     electroweak couplings are taken from KOPPLN
      double precision  clr, xm2, xmg, b, v, a
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     1                  V(4,5),A(4,5)
      double precision  pbar(0:3,2+n_v)
      double precision  uuzz, ddzz, res(3),resv(3)
      double precision  p(0:3,2+n_v)
      double precision  fqz1(0:4,2), fqz2(0:4,2)
      double precision  fqa1(0:4,2), fqa2(0:4,2)
      integer  fsign(2+n_v), sign(2+n_v), nlo, mu, i, k, isig, N
      integer  ifl(2)
      double complex psiz1(2,-1:1,2),psiz2(2,-1:1,2)
      double complex psia1(2,-1:1,2),psia2(2,-1:1,2)
      double complex temp, temp1
      double complex mat(2,-1:1,2)
      double complex mm(2,-1:1), mv(2,-1:1), zero
      double complex m1, m2, m3
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
      logical lnlo, lbox, ldebugm, linit


c   pentagons 
 
      integer i1,i2
      double complex cte
      double complex boxcontr1(3), boxcontr2(3)
      double complex boxcontr3(3), boxcontr4(3)
      double complex boxcontr5(3), boxcontr6(3)
      double complex boxcontr7(3), boxcontr8(3)
      double complex boxcontr1b(3), boxcontr2b(3)
      double complex boxcontr3b(3), boxcontr4b(3)
      double complex boxcontr5b(3), boxcontr6b(3)
      double complex boxcontr7b(3), boxcontr8b(3)

      Logical wardidtest1,wardidtest2


      real*8 p1(0:3),p2(0:3),minus_qa1(0:3),minus_qa2(0:3)
      real*8 minus_qa1z(0:3),minus_qa2z(0:3)


      double complex j5z1z2(3,-1:1,2),jbz1z2(3,-1:1,2)
      double complex j5z2z1(3,-1:1,2),jbz2z1(3,-1:1,2)
      double complex j5z1a2(3,-1:1,2),jbz1a2(3,-1:1,2)
      double complex j5z2a1(3,-1:1,2),jbz2a1(3,-1:1,2)
      double complex j5a1z2(3,-1:1,2),jba1z2(3,-1:1,2)
      double complex j5a2z1(3,-1:1,2),jba2z1(3,-1:1,2)
      double complex j5a1a2(3,-1:1,2),jba1a2(3,-1:1,2)
      double complex j5a2a1(3,-1:1,2),jba2a1(3,-1:1,2)


      double complex mv1, mv2, mv3
      double complex  mv1m,mv2m,mv3m

      double complex  matv(2,-1:1,5),matm(2,-1:1,5)


      data linit /.true./
      save linit
      parameter(ldebugm = .false.)

c---------------------------- BEGIN CODE -----------------------

      lnlo = NLO.ne.0       ! include some virtual stuff if T
      lbox = NLO.eq.-4.or.NLO.eq.1      ! call boxline if T

c     define flavors of external quarks for the 2 subprocesses
      if (linit) then
         linit = .false.
         k = 1                  ! uuzz: u = 3
         ifl(k) = 3
         k = 2                  ! ddzz: d = 4
         ifl(k) = 4
      endif

      do k = 1,2
         do isig = -1,1,2
            do i2=1,2
            mat(k,isig,i2) = (0.0d0,0.0d0)
            enddo
            do i2=1,5
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

         call bra2c(psi(1,isig,2),.true.,p(0,2),isig,qp(0,N),
     &              wp(1,N), psiz1(1,isig,2),fqz1(0,2))
         call bra2c(psi(1,isig,2),.true.,p(0,2),isig,qp(0,N),
     &              atau1(1,N), psia1(1,isig,2),fqa1(0,2))
         call bra2c(psi(1,isig,2),.true.,p(0,2),isig,qm(0,N),
     &              wm(1,N), psiz2(1,isig,2),fqz2(0,2))
         call bra2c(psi(1,isig,2),.true.,p(0,2),isig,qm(0,N),
     &              atau2(1,N), psia2(1,isig,2),fqa2(0,2))
         call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qp(0,N),
     &              wp(1,N), psiz1(1,isig,1),fqz1(0,1))
         call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qp(0,N),
     &              atau1(1,N), psia1(1,isig,1),fqa1(0,1))
         call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qm(0,N),
     &              wm(1,N), psiz2(1,isig,1),fqz2(0,1))
         call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qm(0,N),
     &              atau2(1,N), psia2(1,isig,1),fqa2(0,1))

      enddo

*******************************************************
C Calculation of born matrix elements
*******************************************************

c q qbar -> Z Z 
      do isig = -1,1,2
cc ZZ
         m1=+s1c(psi(1,isig,2),wp(1,N),.true.,isig,psiz2(1,isig,1))
     1      +s1c(psi(1,isig,2),wm(1,N),.true.,isig,psiz1(1,isig,1))
cc ZA+AZ 
         m2=+s1c(psi(1,isig,2),wp(1,N),.true.,isig,psia2(1,isig,1))
     1      +s1c(psi(1,isig,2),wm(1,N),.true.,isig,psia1(1,isig,1))
     1      +s1c(psi(1,isig,2),atau1(1,N),.true.,isig,psiz2(1,isig,1))
     1      +s1c(psi(1,isig,2),atau2(1,N),.true.,isig,psiz1(1,isig,1))
cc AA
         m3=+s1c(psi(1,isig,2),atau1(1,N),.true.,isig,psia2(1,isig,1))
     1      +s1c(psi(1,isig,2),atau2(1,N),.true.,isig,psia1(1,isig,1))

         mat(1,isig,1) =  m1 * clr(3,2,isig)**2                 ! uu -> ZZ
     1                  + m2 * clr(3,2,isig)*clr(3,1,isig)
     1                  + m3 * clr(3,1,isig)**2
         mat(2,isig,1) =  m1 * clr(4,2,isig)**2                 ! dd -> ZZ
     1                  + m2 * clr(4,2,isig)*clr(4,1,isig)
     1                  + m3 * clr(4,1,isig)**2
      enddo

c q qbar -> Z/A 

      call curr6(1,psi(1,-1,2),p(0,2),psi(1,-1,1),p(0,1),jqq(0,-1))

c contract with wz/wa tensor
      do isig = -1,1,2
        mat(1,isig,2) =  -dotcc(wz(1,N),jqq(0,isig))*clr(3,2,isig)
     &                   -dotcc(wa(1,N),jqq(0,isig))*clr(3,1,isig)
        mat(2,isig,2) =  -dotcc(wz(1,N),jqq(0,isig))*clr(4,2,isig)
     &                   -dotcc(wa(1,N),jqq(0,isig))*clr(4,1,isig)
        matv(1,isig,2) = (0.0d0,0.0d0)
        matv(2,isig,2) = (0.0d0,0.0d0)
      enddo

*******************************************************
* Calculation of virtual contribution matrix elements
*******************************************************
#ifdef WITH_NLO
      if (lbox) then
c Calculation of processes qq -> A A ---------------------
        scale=-2d0*dotrr(p(0,1),p(0,2))

        do mu = 0,3
           p1(mu) = pbar(mu,1)
           p2(mu) = pbar(mu,2)
           minus_qa1z(mu) = -qp(mu,1)
           minus_qa2z(mu) = -qm(mu,1)
        enddo
      
        cte=(3*pi*(0d0,1d0)-4.0d0*pi*pi/3.0d0+7.0d0)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        do isig = -1,1,2
c Z1 Z2        
          call boxlinemm(p1,minus_qa1z,minus_qa2z,p2,
     1       psi(1,isig,2),psi(1,isig,1),wp(1,N),wm(1,N),
     2       isig,scale,-3*isig,-isig,boxcontr1,j5z1z2(1,isig,1),
     3       boxcontr1b,jbz1z2(1,isig,1))
c Z1 A2        
          call boxlinemm(p1,minus_qa1z,minus_qa2z,p2,
     1       psi(1,isig,2),psi(1,isig,1),wp(1,N),atau2(1,N),
     2       isig,scale,-1,-1,boxcontr3,j5z1a2(1,isig,1),
     3       boxcontr3b,jbz1a2(1,isig,1))
c A1 Z2        
          call boxlinemm(p1,minus_qa1z,minus_qa2z,p2,
     1       psi(1,isig,2),psi(1,isig,1),atau1(1,N),wm(1,N),
     2       isig,scale,-1,-1,boxcontr5,j5a1z2(1,isig,1),
     3       boxcontr5b,jba1z2(1,isig,1))
c A1 A2        
          call boxlinemm(p1,minus_qa1z,minus_qa2z,p2,
     1       psi(1,isig,2),psi(1,isig,1),atau1(1,N),atau2(1,N),
     2       isig,scale,-1,-1,boxcontr7,j5a1a2(1,isig,1),
     3       boxcontr7b,jba1a2(1,isig,1))
        enddo
        do isig = -1,1,2
c Z2 Z1        
          call boxlinemm(p1,minus_qa2z,minus_qa1z,p2,
     1       psi(1,isig,2),psi(1,isig,1),wm(1,N),wp(1,N),
     2       isig,scale,-3*isig,-isig,boxcontr2,j5z2z1(1,isig,1),
     3       boxcontr2b,jbz2z1(1,isig,1))
c A2 Z1        
          call boxlinemm(p1,minus_qa2z,minus_qa1z,p2,
     1       psi(1,isig,2),psi(1,isig,1),atau2(1,N),wp(1,N),
     2       isig,scale,-1,-1,boxcontr4,j5a2z1(1,isig,1),
     3       boxcontr4b,jba2z1(1,isig,1))
c Z2 A1        
          call boxlinemm(p1,minus_qa2z,minus_qa1z,p2,
     1       psi(1,isig,2),psi(1,isig,1),wm(1,N),atau1(1,N),
     2       isig,scale,-1,-1,boxcontr6,j5z2a1(1,isig,1),
     3       boxcontr6b,jbz2a1(1,isig,1))
c A2 A1        
          call boxlinemm(p1,minus_qa2z,minus_qa1z,p2,
     1       psi(1,isig,2),psi(1,isig,1),atau2(1,N),atau1(1,N),
     2       isig,scale,-1,-1,boxcontr8,j5a2a1(1,isig,1),
     3       boxcontr8b,jba2a1(1,isig,1))
        enddo
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

         call box_gauge(p1,minus_qa1,minus_qa2,p2,
     2       scale,boxcontr1,boxcontr1b,wardidtest1)
         call box_gauge(p1,minus_qa2,minus_qa1,p2,
     2       scale,boxcontr2,boxcontr2b,wardidtest2)
  
         wardidtest1=wardidtest1.and.wardidtest2
 
         if (.not.wardidtest1) then
           do i1=-1,1,2
             do i2=1,2
               j5z1z2(1,i1,i2) = (0.0d0, 0.0d0)
               j5z2z1(1,i1,i2) = (0.0d0, 0.0d0)
               j5z1a2(1,i1,i2) = (0.0d0, 0.0d0)
               j5z2a1(1,i1,i2) = (0.0d0, 0.0d0)
               j5a1z2(1,i1,i2) = (0.0d0, 0.0d0)
               j5a2z1(1,i1,i2) = (0.0d0, 0.0d0)
               j5a1a2(1,i1,i2) = (0.0d0, 0.0d0)
               j5a2a1(1,i1,i2) = (0.0d0, 0.0d0)
             enddo
           enddo
         else
           do i1=-1,1,2
             do i2=1,1
                j5z1z2(1,i1,i2)=j5z1z2(1,i1,i2)+cte*jbz1z2(1,i1,i2)
                j5z2z1(1,i1,i2)=j5z2z1(1,i1,i2)+cte*jbz2z1(1,i1,i2)
                j5z1a2(1,i1,i2)=j5z1a2(1,i1,i2)+cte*jbz1a2(1,i1,i2)
                j5z2a1(1,i1,i2)=j5z2a1(1,i1,i2)+cte*jbz2a1(1,i1,i2)
                j5a1z2(1,i1,i2)=j5a1z2(1,i1,i2)+cte*jba1z2(1,i1,i2)
                j5a2z1(1,i1,i2)=j5a2z1(1,i1,i2)+cte*jba2z1(1,i1,i2)
                j5a1a2(1,i1,i2)=j5a1a2(1,i1,i2)+cte*jba1a2(1,i1,i2)
                j5a2a1(1,i1,i2)=j5a2a1(1,i1,i2)+cte*jba2a1(1,i1,i2)
             enddo
           enddo
         endif

c   q qbar -> Z Z 

         do isig=-1,1,2
           mv1=0d0
           mv2=0d0
           mv3=0d0
           mv1m=0d0
           mv2m=0d0
           mv3m=0d0
           
           mv1  = j5z1z2(1,isig,1) + j5z2z1(1,isig,1)
           mv1m = jbz1z2(1,isig,1) + jbz2z1(1,isig,1)
           mv2  = j5z1a2(1,isig,1) + j5z2a1(1,isig,1)
     &          + j5a1z2(1,isig,1) + j5a2z1(1,isig,1)
           mv2m = jbz1a2(1,isig,1) + jbz2a1(1,isig,1)
     &          + jba1z2(1,isig,1) + jba2z1(1,isig,1)
           mv3  = j5a1a2(1,isig,1) + j5a2a1(1,isig,1)
           mv3m = jba1a2(1,isig,1) + jba2a1(1,isig,1)

           matv(1,isig,1) =  mv1 * clr(3,2,isig)**2   ! uu -> ZZ
     &                     + mv2 * clr(3,2,isig)*clr(3,1,isig)
     &                     + mv3 * clr(3,1,isig)**2
           matv(2,isig,1) =  mv1 * clr(4,2,isig)**2   ! dd -> ZZ
     &                     + mv2 * clr(4,2,isig)*clr(4,1,isig)
     &                     + mv3 * clr(4,1,isig)**2
           matm(1,isig,1) =  mv1m * clr(3,2,isig)**2   ! uu -> ZZ
     &                     + mv2m * clr(3,2,isig)*clr(3,1,isig)
     &                     + mv3m * clr(3,1,isig)**2
           matm(2,isig,1) =  mv1m * clr(4,2,isig)**2   ! dd -> ZZ
     &                     + mv2m * clr(4,2,isig)*clr(4,1,isig)
     &                     + mv3m * clr(4,1,isig)**2
 
         enddo

      endif  !lbox

#endif /* WITH_NLO */
      
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
c print results from checks:

         do k = 1,2
           do isig = -1,1,2
             do i2 = 1,2
              print*,'mat(',k,',',isig,',',i2,') = ',mat(k,isig,i2)
              temp = dreal(mat(k,isig,i2))**2+dimag(mat(k,isig,i2))**2
              print*,'madgr(',i2,',',isig,',',k,') = ',tamp(i2,isig,k)
              temp1 = dreal(tamp(i2,isig,k))**2+dimag(tamp(i2,isig,k))**2
              print*,'temp/temp1 = ',temp/temp1
             enddo
           enddo
         enddo
         
      endif                     ! end of debugging block

c------------------------------------------------------

c sum the graphs, square them and map them onto uuzz and ddzz

      do k = 1,2
         res(k) = 0
         resv(k) = 0
         do isig = -1,1,2
            mm(k,isig) = 0
            mm(k,isig) = mm(k,isig) +  mat(k,isig,1) + mat(k,isig,2) 
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
      
      uuzz = res(1)
      ddzz = res(2)

      return
      end

