c****************subroutine QQWA *************************************
c
c	Johannes Bellm, <johannes@particle.uni-karlsruhe.de>
c	Initial version: Juli 2011
c	Last modified: 
C
c        Adapted from qqVV by Vera Hankele
c
C  QQWA calculates the matrix elements**2 for qq -> WA -> ve e+ A
C
C        q1 q2    ---->   W A,   W ---> f4-bar f3 
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

      subroutine qqWA(pbar, sign, nlo, N, udwa, bos)

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
      double precision  udwa, res(2),resv(2)
      double precision  p(0:3,2+n_v),
     1                  fqp(0:4,2),fqa(0:4,2)
      integer  sign(2+n_v), nlo, mu, i, k,
     1         isig1, N
      integer  ifl(2), bos
      double complex psiwp(2,2) , psian(2,2), 
     1               jwp(0:5,-1:1),jan(0:5,-1:1)

      double complex mat(2,-1:1,3), matv(2,-1:1,3), mborn1

      double complex mm(2,-1:1), mv(2,-1:1), zero
      double complex mwp, m12, m22 

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
      double precision minus_qp(0:3),minus_qa(0:3)
      integer wardidtest3, wardidtest4


      double complex gaugebox3(2), gaugebox3_born(2)
      double complex gaugebox4(2), gaugebox4_born(2)
      double complex result(4),result_born(4)
      double complex cte
      data linit /.true./
      save linit
      parameter(ldebugm = .false.)

      lnlo = NLO.ne.0       ! include some virtual stuff if T
      lbox = NLO.eq.-4.or.NLO.eq.1      ! call boxline if T

c define flavors of external quarks
      if (linit) then
         linit = .false.
         if (bos.eq.31) then
           ifl(1) = 3             
           ifl(2) = 4
         else if (bos.eq.41) then
           ifl(1) = 4             
           ifl(2) = 3
         else 
           print *, "Wrong boson ids:", bos
           stop
         endif
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

      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,qa(0,N),an(1,N),
     &              psian(1,2),fqa(0,2))
        
      call curr6(-1,psiwp(1,2),fqp(0,2),
     &              psi(1,-1,1),p(0,1), jwp(0,-1) )

      call curr6(-1,psian(1,2),fqa(0,2),
     &              psi(1,-1,1),p(0,1), jan(0,-1) )

      m12 = dotcc(jwp(0,-1),an(1,N))
      m22 = dotcc(jan(0,-1),wp(1,N))   


      mat(1,-1,1) = clr(ifl(2),3,-1)* clr(ifl(1),1,-1) * m12 
     1            + clr(ifl(2),1,-1)* clr(ifl(2),3,-1) * m22


c       clr( ferm , bos , hel )


C  The fermion indices are:  1: neutrino of electron (muon,tau)
C                            2: electron (muon,tau)
C                            3: up quark (charm,top)
C                            4: down quark (strange,bottom)
C  The boson indices are     1: photon
C                            2: Z0 boson
C                            3: W+ boson
C                            4: W- boson
C                            5: gluon
C                            6: Higgs
C  The helicity indices are  -1: lefthanded
C                            +1: righthanded


c-------  S-CHANNEL --------------------------------------!

c get the f-fbar currents J21^mu=jqq(mu,*,1) 

      call curr6(1,psi(1,-1,2),p(0,2),psi(1,-1,1),p(0,1),jqq(0,-1,1))

c contract with wtowz tensor

      
         mwp = -dotcc(wtowa(0,N),jqq(0,-1,1))

         
            mat(1,-1,2) = mwp*clr(3,3,-1)
            
            matv(1,-1,1) = (0d0,0d0)


c****************************************************************
c
c      BOX CONTRIBUTIONS
c
c****************************************************************

#ifdef WITH_NLO
      if (lbox) then
         scale = -2d0*dotrr(p(0,1),p(0,2))

c calculate born matrix element and store it in mborn for gauge check:
         mborn1=mat(1,-1,1)+mat(1,-1,2)
         

         do mu = 0,3
            p1(mu) = pbar(mu,1)
            p2(mu) = pbar(mu,2)
            minus_qp(mu) = -qp(mu,1)
!             minus_qz(mu) = -p1(mu)-p2(mu)+qp(mu,1)
            minus_qa(mu) = -p1(mu)-p2(mu)+qp(mu,1)
         enddo

         cte=(3*pi*(0d0,1d0)-4.0d0*pi**2/3.0d0+7.0d0)

         call boxlinemm(p1,minus_qp,minus_qa,p2,
     1    psi(1,-1,2),psi(1,-1,1),wp(1,1), an(1,1),
     2    -1,scale,3,1,gaugebox3,result(3),gaugebox3_born,
     3    result_born(3))

         call boxlinemm(p1,minus_qa,minus_qp,p2,
     1    psi(1,-1,2),psi(1,-1,1),an(1,1), wp(1,1),
     2    -1,scale,3,1,gaugebox4,result(4),gaugebox4_born,
     3    result_born(4))


         call Wardtest_Box(wp(1,1), qp(0,1), an(1,1), qa(0,1),
     1        gaugebox3, gaugebox3_born, mborn1, wardidtest3, als(1,1))

         call Wardtest_Box(an(1,1), qa(0,1), wp(1,1), qp(0,1),
     1        gaugebox4, gaugebox4_born, mborn1, wardidtest4, als(1,1))

         if (wardidtest3.eq.1 .or. wardidtest4.eq.1) then
                  matv(1,-1,2) = (0.0d0, 0.0d0)
                
         else           
            result(3)=result(3)+cte*result_born(3)
            result(4)=result(4)+cte*result_born(4)

            matv(1,-1,2) = clr(ifl(1),3,-1) * clr(ifl(2),1,-1) * result(3) !1u 2w+ 3a  4d~ 
     1                   + clr(ifl(1),1,-1) * clr(ifl(1),3,-1) * result(4) !1u 2a  3w+ 4d~ 

c
c                         1-------------~~~~~~~2
c                              (       |
c                              )glu    |                              
c                              (       |
c                         4-------------~~~~~~~3
c
c
         endif
      endif
#endif


c sum the graphs, square them and map them onto udwa 
      

           resv(1) = 0
            mm(1,-1) = 0
            mv(1,-1) = (0d0, 0d0)
           res(1)= 0

               mm(1,-1) = mat(1,-1,1)+mat(1,-1,2)
            
            res(1)=   dreal(mm(1,-1))**2
     &                     + dimag(mm(1,-1))**2 
                         

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
            res(1) = (res(1)+resv(1))*3d0      ! 3 is the color sum factor
         else
            res(1) = resv(1)*3d0               ! 3 is the color sum factor
         endif  
      

      udwa = res(1)


      return
      end



