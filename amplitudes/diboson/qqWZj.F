****************subroutine QQWZj *************************************
c
c	Johannes Bellm, <johannes@particle.physik.uni-karlsruhe.de>
c	Initial version:  2011 Juli
c	Last modified:    
c
cc       Adapted from qqVV by Vera Hankele
c
c  QQWZj calculates the matrix elements**2 for qq -> WZ g -> ve e+ mu- mu+ g
c
c        q1 q2 ----> W Z g,   W+ ---> f4-bar f3, Z ---> f6-bar f5
c
c  and crossing related processes. Pauli interference terms for
c  identical fermions are neglected.
c
c  Gluon polarization vectors are given in the kartesian basis 
c  (comparison with madgraph amplitudes tampr requires to express these in 
c   kart. basis also, even though they are given in the helicity basis 
c   a priori). 
c
c  This code includes only real emission contributions, i.e.
c
c      return udwz = |M_real|^2   etc.
c
c	fpials is attached only in the end of the code
c
c index j = 2:3 indicates, whether g is emitted from 
c		upper line (j=2) or lower line (j=3)
c	l is the gluon polarization in the kartesian basis (l=1,2)
c		l=0 stands for building blocks without gluon emission
c	k is the process ID (1:udwz)
c	isig1/isig3 are the helicities of parton1/2 
c
c*********************************************************************


      subroutine qqWZj(pbar,sign, qbar,gsign, udwz, bos)

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/coupl.inc"

c electroweak couplings are taken from KOPPLN

      double precision  clr, xm2, xmg, b, v, a, udwzg
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     1                  V(4,5),A(4,5)

c variables for the main part of the program      

      real*8 fpials(2:3), fpi
      parameter (fpi=4d0*pi)

      double precision  pbar(0:3,2+n_v), qbar(0:4), q(0:3)
      double precision  udwz, res(2)
      double precision  p(0:3,2+n_v), pq(0:4,2)
      double precision pgWP(0:4,2), pWP(0:4,2),
     1                 pgZb(0:4,2), pgAb(0:4,2), pZ(0:4,2), pA(0:4,2)
      integer  sign(2+n_v), gsign, mu, i, k, isig
      integer  ifl(2), bos
      integer  l   ! gluon polariz. (l=0:no g, l=1,2:g in kartesian basis)
      double complex mat(2,-1:1,0:2,9),temp, temp1
      double complex mata(2,-1:1,0:2,9),
     1 		     matz(2,-1:1,0:2,9)
      double complex mm(2,-1:1,2)
      double complex  m1(-1:1,2,3), m2(-1:1,2,3),m3(-1:1,2,3),m4(-1:1,2,3)
      double complex  raz! for checks only
      double precision eps(0:3,2) ! g in kartesian basis 
      double complex psi(2,-1:1,4), psiWP(2,-1:1,2),
     1               psiZb(2,-1:1,2), psiAb(2,-1:1,2)
      double complex psiWPg(2,-1:1, 2,2),
     1               psiZg(2,-1:1, 2,2),psiAg(2,-1:1, 2,2)
      double complex braketg(2,-1:1,2,2)
      double complex jqgaq(0:5,-1:1,2), jqagq(0:5,-1:1,2), 
     &               jqqa(0:5,-1:1,2)
      double complex jqqgWP(0:5,-1:1,2,0:3),
     1               jqqgZ(0:5,-1:1,2,0:3), jqqgA(0:5,-1:1,2,0:3)
      double complex braketgWP(2,-1:1,2,2),
     1               braketgZ(2,-1:1,2,2),braketgA(2,-1:1,2,2)   
      double precision  pwpg(0:4,4)
      double precision  pag(0:4,4),pzg(0:4,4)
      double complex mw(-1:1,2,0:2)

      double complex contract_Tjj,! contract_T1j,contract_T2j,
     1 		     dotcc, dotrc, dotqj, s1c
      external contract_Tjj, !contract_T1j,contract_T2j,
     1 	       dotcc, dotrc, dotqj, s1c
      integer lh
      double complex madgr(9,-1:1,-1:1,2),m1kb(2,-1:1,9),m2kb(2,-1:1,9) !checks

      double complex im
      parameter (im = (0d0,1d0))
      logical ldebug, linit, lerror, ldebugm
     
      save ifl, linit
	
      data linit /.true./, lerror /.true./, ldebug /.false./
      parameter(ldebugm = .false.) 

c  ---------------------------------------------------------------------
 
c initialize & precompute stuff needed below:

c fix strong coupling gs**2 for the 2 quark lines:

      fpials(2) = fpi*als(1,1)
      fpials(3) = fpi*als(2,1)
c       fpials(2) = fpi*alfas

c define flavors of external quarks for the 2 NC subprocesses

      if (linit) then
         linit = .false.
         if (bos.eq.32) then
           ifl(1) = 3             
           ifl(2) = 4
         else if (bos.eq.42) then
           ifl(1) = 4             
           ifl(2) = 3
         else 
           print *, "Wrong boson ids:", bos
           stop
         endif
      endif

      do k = 1,1
         do isig = -1,1,2
	    do l = 0,2
               do i = 1,9
                  mat(k,isig,l,i)  = 0
                  mata(k,isig,l,i) = 0
                  matz(k,isig,l,i) = 0
               enddo
            enddo
         enddo
      enddo

c define the internal momenta

      do mu = 0,3
         do i = 1,2+n_v
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo

	 q(mu) = qbar(mu)*gsign

      enddo

c---------calculation of matrix elements--------------------

c get the external quark spinors (including factor sqrt(2E) )

      call psi0m(2,pbar(0,1),sign(1),psi)

c     NOTE: psi(2,-1:1,2) is a two component spinor with helicity -1
c     or 1.  The last entry identifies the fermion.  If this entry is
c     odd psi is a ket, if even psi is a bra.
c     psi(1,isig1,1) = |1>_isig1
c     psi(1,isig1,2) = <2|_isig

c  Get the gluon polarization vector and the gluon emission spinors

      do l = 1,2	! 2 gluon polarizations

         call polvec(qbar,l,eps(0,l))  ! get gluon pol.vectors

c---- QCD gaugecheck --- eps_mu -> k_mu --- k_mu M^mu = 0 ------------

c         do mu = 0,3                   !QCD gaugecheck      
c            eps(mu,l) = qbar(mu)
c         enddo

c--- End of QCD gaugecheck -------------------------------------------    

         do isig = -1,-1         ! not needed in this process	! fermion helicity 

c     NOTES for bras and kets: .true. if psi is a 2-spinor of the chi
c     form as output by psi0m, .false. otherwise.  the last entry is
c     the sum of the two momenta (p plus q) and effectively the
c     momentum of the new spinor.
 
            call ket2r(psi(1,isig,1),.true.,p(0,1),isig,q,eps(0,l),
     $           braketg(1,isig,1,l),pq(0,1))      ! |q,1>_l,isig
            call bra2r(psi(1,isig,2),.true.,p(0,2),isig,q,eps(0,l),
     $           braketg(1,isig,2,l),pq(0,2))      ! <2,q|_l,isig2
         enddo
      enddo

c     braketg contains the free quark spinors multiplied by a fermion
c     propagator and a gluon eps_slash. 
c     NOTATION: braketg(2 component spinor, isig =-1 or 1 (fermion hel.),
c     fermion ID = 1:4, gluon polarization l=1:2)

c------- T-CHANNEL -----------------

c possibility 1 - emission from upper line

      do l=1,2
         do isig = -1,-1         ! not needed in this process
   
         call ket2c(braketg(1,isig,1,l),.false.,pq(0,1),isig,
     &        qp(0,1),wp(1,1),braketgWP(1,isig,1,l),pgWP(0,1)) !|W+,g,1>_l,isig

         call ket2c(braketg(1,isig,1,l),.false.,pq(0,1),isig,
     &        qz(0,1),zn(1,1),braketgZ(1,isig,1,l),pgZb(0,1)) !|Z,g,1>_l,isig
         call ket2c(braketg(1,isig,1,l),.false.,pq(0,1),isig,
     &        qa(0,1),an(1,1),braketgA(1,isig,1,l),pgAb(0,1)) !|A,g,1>_l,isig

         call curr6(1,psi(1,isig,2),p(0,2),
     &           braketgWP(1,isig,1,l),pgWP(0,1),jqqgWP(0,isig,l,1))
         call curr6(1,psi(1,isig,2),p(0,2),
     &           braketgZ(1,isig,1,l),pgZb(0,1),jqqgZ(0,isig,l,1))
         call curr6(1,psi(1,isig,2),p(0,2),
     &           braketgA(1,isig,1,l),pgAb(0,1),jqqgA(0,isig,l,1))


      m1(isig,l,1) = dotcc(jqqgWP(0,isig,l,1),zn(1,1))
      m2(isig,l,1) = dotcc(jqqgWP(0,isig,l,1),an(1,1))
      m3(isig,l,1) = dotcc(jqqgZ(0,isig,l,1),wp(1,1))
      m4(isig,l,1) = dotcc(jqqgA(0,isig,l,1),wp(1,1))      

      mat(1,isig,l,1) = clr(ifl(1),3,isig)*clr(ifl(2),2,isig) * m1(isig,l,1)
     1                + clr(ifl(1),3,isig)*clr(ifl(2),1,isig) * m2(isig,l,1)
     2                + clr(ifl(1),2,isig)*clr(ifl(1),3,isig) * m3(isig,l,1)
     3                + clr(ifl(1),1,isig)*clr(ifl(1),3,isig) * m4(isig,l,1)

c                   clr(ferm,bos,hel)
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

         enddo
      enddo

c             - emission from lower line

      do l=1,2
         do isig = -1,-1         ! not needed in this process

      call bra2c(braketg(1,isig,2,l),.false.,pq(0,2),isig,
     &       qp(0,1),wp(1,1),braketgWP(1,isig,2,l),pgWP(0,2)) !<2,g,W+|_l,isig

      call bra2c(braketg(1,isig,2,l),.false.,pq(0,2),isig,
     &       qz(0,1),zn(1,1),braketgZ(1,isig,2,l),pgZb(0,2)) !<2,g,Z|_l,isig

      call bra2c(braketg(1,isig,2,l),.false.,pq(0,2),isig,
     &       qa(0,1),an(1,1),braketgA(1,isig,2,l),pgAb(0,2)) !<2,g,A|_l,isig

      call curr6(1,braketgWP(1,isig,2,l),pgWP(0,2),
     &            psi(1,isig,1),p(0,1),jqqgWP(0,isig,l,2))

      call curr6(1,braketgZ(1,isig,2,l),pgWP(0,2),
     &            psi(1,isig,1),p(0,1),jqqgZ(0,isig,l,2))

      call curr6(1,braketgA(1,isig,2,l),pgWP(0,2),
     &            psi(1,isig,1),p(0,1),jqqgA(0,isig,l,2))

      m1(isig,l,2) = dotcc(jqqgWP(0,isig,l,2),zn(1,1))
      m2(isig,l,2) = dotcc(jqqgWP(0,isig,l,2),an(1,1))
      m3(isig,l,2) = dotcc(jqqgZ(0,isig,l,2),wp(1,1))
      m4(isig,l,2) = dotcc(jqqgA(0,isig,l,2),wp(1,1))

      mat(1,isig,l,2) = clr(ifl(2),3,isig)*clr(ifl(1),2,isig) * m1(isig,l,2)
     1                + clr(ifl(2),3,isig)*clr(ifl(1),1,isig) * m2(isig,l,2)
     2                + clr(ifl(2),2,isig)*clr(ifl(2),3,isig) * m3(isig,l,2)
     3                + clr(ifl(2),1,isig)*clr(ifl(2),3,isig) * m4(isig,l,2)

         enddo
      enddo

c             - emission from intermediate quark line

      do l=1,2
         do isig = -1,-1         ! not needed in this process

      call bra2c(psi(1,isig,2),.true.,p(0,2),isig,
     &       qp(0,1),wp(1,1),psiWP(1,isig,2),pWP(0,2)) !<2,W+|_isig

      call bra2c(psi(1,isig,2),.true.,p(0,2),isig,
     &       qz(0,1),zn(1,1),psiZb(1,isig,2),pZ(0,2)) !<2,Z|_isig1

      call bra2c(psi(1,isig,2),.true.,p(0,2),isig,
     &       qa(0,1),an(1,1),psiAb(1,isig,2),pA(0,2))  !<2,A|_isig

      call bra2r(psiWP(1,isig,2),.false.,pWP(0,2),isig,
     &            q,eps(0,l),psiWPg(1,isig,2,l),pWPg(0,2)) !<2,W+,g|_l,isig

      call bra2r(psiZb(1,isig,2),.false.,pZ(0,2),isig,
     &            q,eps(0,l),psiZg(1,isig,2,l),pZg(0,2)) !<2,Z,g|_l,isig

      call bra2r(psiAb(1,isig,2),.false.,pA(0,2),isig,
     &            q,eps(0,l),psiAg(1,isig,2,l),pAg(0,2)) !<2,A,g|_l,isig

      call curr6(1,psiWPg(1,isig,2,l),pWPg(0,2),
     &            psi(1,isig,1),p(0,1),jqqgWP(0,isig,l,3))

      call curr6(1,psiZg(1,isig,2,l),pZg(0,2),
     &            psi(1,isig,1),p(0,1),jqqgZ(0,isig,l,3))

      call curr6(1,psiAg(1,isig,2,l),pAg(0,2),
     &            psi(1,isig,1),p(0,1),jqqgA(0,isig,l,3))

      m1(isig,l,3) = dotcc(jqqgWP(0,isig,l,3),zn(1,1))
      m2(isig,l,3) = dotcc(jqqgWP(0,isig,l,3),an(1,1))
      m3(isig,l,3) = dotcc(jqqgZ(0,isig,l,3),wp(1,1))
      m4(isig,l,3) = dotcc(jqqgA(0,isig,l,3),wp(1,1))

      mat(1,isig,l,3) = clr(ifl(2),3,isig)*clr(ifl(1),2,isig) * m1(isig,l,3)
     1                + clr(ifl(2),3,isig)*clr(ifl(1),1,isig) * m2(isig,l,3)
     2                + clr(ifl(2),2,isig)*clr(ifl(2),3,isig) * m3(isig,l,3)
     3                + clr(ifl(2),1,isig)*clr(ifl(2),3,isig) * m4(isig,l,3)

         enddo
      enddo

      do l=1,2
         do isig = -1,-1         ! not needed in this process

            mat(1,isig,l,8) = mat(1,isig,l,1) + mat(1,isig,l,2)
     &                      + mat(1,isig,l,3)


         enddo
      enddo

c------- end of T-CHANNEL computation

c-------  S-CHANNEL ----------------

c     Get the f-fbar currents with one gluon radiated from the
c     current line.  There are two terms, one for gluon emission to
c     either side of the ffV vertex:

      do l=1,2
         do isig = -1,-1         ! not needed in this process

         call curr6(1,psi(1,-1,2),p(0,2),
     &              braketg(1,-1,1,l),pq(0,1),jqgaq(0,-1,l))

         call curr6(1,braketg(1,-1,2,l),pq(0,2),
     &              psi(1,-1,1),p(0,1),jqagq(0,-1,l))

            do mu = 0,5                                               
	       jqqa(mu,isig,l) = jqgaq(mu,isig,l) + jqagq(mu,isig,l)
c                            = (<2|gam.mu|q,1>+<2,q|gam.mu|1>)_l,isig)
            enddo
         enddo
      enddo

      do l=1,2
         do isig = -1,-1         ! not needed in this process
            mw(isig,l,1) = -dotcc(wtowz(0,1),jqqa(0,isig,l))

           mat(1,isig,l,4)= mw(isig,l,1)*clr(3,3,isig)
           mat(1,isig,l,9)= mat(1,isig,l,4)

         enddo
      enddo   !end gluon polarization l

c-------- end of S-CHANNEL computation ---------------------------

c------------------------------------------------------
c comparison with madgraph

#ifdef WITH_MADGRAPH
      if (ldebugm) then
         if (bos.eq.32) then
            call sud_wpzgL(pbar,sign,udwzg)          !madgraph
         else if (bos.eq.42) then
            call sdu_wmzgL(pbar,sign,udwzg)          !madgraph
         endif

            isig=-1 

            do i = 1,2 
                  m1kb(1,isig,i) = 0d0 ! kartesian basis l = 1
                  m2kb(1,isig,i) = 0d0 ! kartesian basis l = 2
                  do lh = -1,1,2
                     madgr(i,isig,lh,1) = 0d0
                  enddo        
            enddo

            do lh = -1,1,2
                  do i =1,2
                     
                     madgr(i,isig,lh,1) = tampr(i,isig,lh,1)
                     
                  enddo
            enddo

            do i=1,2     
                  m1kb(1,isig,i) = (madgr(i,isig,-1,1) - 
     &                 madgr(i,isig,1,1))/sqrt(2d0)       
                  m2kb(1,isig,i) = -im*(madgr(i,isig,-1,1) + 
     &                 madgr(i,isig,1,1))/sqrt(2d0)
            enddo

c print results from checks:

         do i=1,2
                  temp = dreal(mat(1,isig,1,7+i))**2+dimag(mat(1,isig,1,7+i))**2
                  temp1 = dreal(m1kb(1,isig,i))**2+dimag(m1kb(1,isig,i))**2
            if(abs(temp/temp1 -1).ge.1.0d-6) then
                  print*,"mat(",1,",",isig,",1,",7+i,") = ",mat(1,isig,1,7+i)
                  print*,"m1kb(",1,",",isig,",",i,") = ",m1kb(1,isig,i)
                  print*,"temp/temp1 = ",temp/temp1
            endif
                  temp = dreal(mat(1,isig,2,7+i))**2+dimag(mat(1,isig,2,7+i))**2
                  temp1 = dreal(m2kb(1,isig,i))**2+dimag(m2kb(1,isig,i))**2
            if(abs(temp/temp1 -1).ge.1.0d-6) then
                  print*,"mat(",1,",",isig,",2,",7+i,") = ",mat(1,isig,2,7+i)
                  print*,"m2kb(",1,",",isig,",",i,") = ",m2kb(1,isig,i)
                  print*,"temp/temp1 = ",temp/temp1
            endif
         enddo
      endif ! end of debugging block
#endif

c------------------------------------------------------

c-----------------END of calculation of matrix elements------------------

c sum the graphs, square them and map them onto udwz and ddww
         res(1) = 0
         do isig = -1,-1         ! not needed in this process
 	    do l = 1,2
               mm(1,isig,l) = 0
               do i = 1,7
                  mm(1,isig,l) = mm(1,isig,l) +  mat(1,isig,l,i)                
               enddo
       

               res(1)= res(1) + dreal(mm(1,isig,l))**2
     &                        + dimag(mm(1,isig,l))**2

            enddo
         enddo
         res(1) = res(1)*4d0*fpials(2)  ! C_2*3 is the color factor     


      udwz = res(1)


      return
      end


