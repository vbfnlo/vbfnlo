****************subroutine QQWHj *************************************
c
c  Adapted from qqWAj by Johannes Bellm and qqVVj by Vera Hankele
c
c  Robin Roth <robin@particle.physik.uni-karlsruhe.de>
c  Initial version: November 2012
c  Last modified: November 2012
c
c  QQWHj calculates the matrix elements**2 for qq -> WH g -> ve e+ H g
c
c        q1 q2 ----> W H g,   W ---> f4-bar f3
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
c      return udwa = |M_real|^2   etc.
c
c	fpials is attached only in the end of the code
c
c index j = 2:3 indicates, whether g is emitted from 
c		upper line (j=2) or lower line (j=3)
c	glupol is the gluon polarization in the kartesian basis (glupol=1,2)
c		glupol=0 stands for building blocks without gluon emission
c	k is the process ID (1:udwa)
c	isig1/isig3 are the helicities of parton1/2 
c
c*********************************************************************


      subroutine qqWHj(pbar,sign, qbar,gsign, udwh, bos)

c*********************************************************************         
c pbar(0:3,1) and pbar(0:3,2) are incoming partons
c pbar(0:3, 3) ve
c pbar(0:3, 4) e
c pbar(0:3, 5) h 
c qbar(0:3) gluon momentum
c*********************************************************************

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/coupl.inc"

c electroweak couplings are taken from KOPPLN

      double precision  clr, xm2, xmg, b, v, a
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     1                  V(4,5),A(4,5)

c variables for the main part of the program      

      real*8 fpials(2:3), fpi
      parameter (fpi=4d0*pi)

      double precision  pbar(0:3,2+n_v), qbar(0:4), q(0:3)

      double precision  udwh, res(2),udwhg
      double precision  p(0:3,2+n_v), pq(0:4,2)

      integer sign(2+n_v), gsign, mu, i, k, isig
      integer bos
      integer glupol   ! gluon polariz. (glupol=0:no g, glupol=1,2:g in kartesian basis)

      double complex mat(2,-1:1,0:2,9),temp, temp1
      double complex mm(2,-1:1,2)
      double precision eps(0:3,2) ! g in kartesian basis 
      double complex psi(2,-1:1,4)
      double complex braketg(2,-1:1,2,2)
      double complex jqghq(0:5,-1:1,2), jqhgq(0:5,-1:1,2), 
     &               jqqwh(0:5,-1:1,2)
      double complex mw(-1:1,2,0:2)
      double complex dotcc
      external   dotcc
      integer lh
      double complex madgr(9,-1:1,-1:1,2),m1kb(2,-1:1,9),m2kb(2,-1:1,9) !checks
      double complex im
      parameter (im = (0d0,1d0))
      logical ldebug, linit, ldebugm
     
      save linit

      data linit /.true./
      parameter(ldebug = .false.)
      parameter(ldebugm= .false.) 

c for debugging purposes
      double precision momsum(0:3)
      double precision dotrr
      external dotrr

c  ---------------------------------------------------------------------

c      if (ldebug) print*,"enter qqWHj ---------------------"

c initialize & precompute stuff needed below:

c fix strong coupling gs**2 for the 2 quark lines:
      fpials(2) = fpi*als(1,1)
      fpials(3) = fpi*als(2,1)
c       fpials(2) = fpi*alfas

c define flavors of external quarks for the 2 NC subprocesses

      if (linit) then
         linit = .false.
         if (.not. (bos.eq.3 .or. bos.eq.36 .or. bos.eq.4 .or. bos.eq.46)) then
           print *, "Wrong boson ids:", bos
           stop
         endif
      endif

      do k = 1,1
         do isig = -1,1,2
            do glupol = 0,2
               do i = 1,9
                  mat(k,isig,glupol,i)  = 0
               enddo
            enddo
         enddo
      enddo

c define the internal momenta

      do mu = 0,3
         do i = 1,2+n_v ! = 5 in WHj@LO
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo
         q(mu) = qbar(mu)*gsign
      enddo

      if (ldebug) then
         print*, gsign
         print*, pbar(0:3,1)
         print*, pbar(0:3,2)
         print*, pbar(0:3,3)
         print*, pbar(0:3,6)
         do mu=0,3
         ! momentum conservation test for the simple case 
         ! q qbar -> ... g only !
            momsum(mu) = pbar(mu,1)+pbar(mu,2)-pbar(mu,3)
     1                  -pbar(mu,4)-pbar(mu,5)-pbar(mu,6)
         enddo
         if (dotrr(momsum,momsum) > 1e-10) then
         print*, "momsum ", momsum(0) ,momsum(1), momsum(2), momsum(3)
         endif
      endif 

c---------calculation of matrix elements--------------------

c get the external quark spinors (including factor sqrt(2E) )

      call psi0m(2,pbar(0,1),sign(1),psi)
c     if (ldebug) print*, "psi1 ", psi(1:2,-1,1)
c     if (ldebug) print*, "psi2 ", psi(1:2,-1,2)

c     NOTE: psi(2,-1:1,2) is a two component spinor with helicity -1
c     or 1.  The last entry identifies the fermion.  If this entry is
c     odd psi is a ket, if even psi is a bra.
c     psi(1,isig1,1) = |1>_isig1
c     psi(1,isig1,2) = <2|_isig

c  Get the gluon polarization vector and the gluon emission spinors

      do glupol = 1,2	! 2 gluon polarizations

         call polvec(qbar,glupol,eps(0,glupol))  ! get gluon pol.vectors
c         if (ldebug) print*, "eps ",  eps

c---- QCD gaugecheck --- eps_mu -> k_mu --- k_mu M^mu = 0 ------------

c         do mu = 0,3                   !QCD gaugecheck      
c            eps(mu,glupol) = qbar(mu)
c         enddo

c--- End of QCD gaugecheck -------------------------------------------    

        isig = -1         ! not needed in this process	! fermion helicity 

c     NOTES for bras and kets: .true. if psi is a 2-spinor of the chi
c     form as output by psi0m, .false. otherwise.  the last entry is
c     the sum of the two momenta (p plus q) and effectively the
c     momentum of the new spinor.
 
            call ket2r(psi(1,isig,1),.true.,p(0,1),isig,q,eps(0,glupol),
     $           braketg(1,isig,1,glupol),pq(0,1))      ! |q,1>_l,isig

            call bra2r(psi(1,isig,2),.true.,p(0,2),isig,q,eps(0,glupol),
     $           braketg(1,isig,2,glupol),pq(0,2))      ! <2,q|_l,isig2

      enddo

c     braketg contains the free quark spinors multiplied by a fermion
c     propagator and a gluon eps_slash. 
c     NOTATION: braketg(2 component spinor, isig =-1 or 1 (fermion hel.),
c     fermion ID = 1:4, gluon polarization glupol=1:2)

c  
c  c                   clr(ferm,bos,hel)
c  C  The fermion indices are:  1: neutrino of electron (muon,tau)
c  C                            2: electron (muon,tau)
c  C                            3: up quark (charm,top)
c  C                            4: down quark (strange,bottom)
c  C  The boson indices are     1: photon
c  C                            2: Z0 boson
c  C                            3: W+ boson
c  C                            4: W- boson
c  C                            5: gluon
c  C                            6: Higgs
c  C  The helicity indices are  -1: lefthanded
c  C                            +1: righthanded
c  

c-------  S-CHANNEL ----------------

c     Get the f-fbar currents with one gluon radiated from the
c     current line.  There are two terms, one for gluon emission to
c     either side of the ffV vertex:

      do glupol=1,2
         isig = -1         ! not needed in this process

         call curr6(1,psi(1,-1,2),p(0,2),
     &              braketg(1,-1,1,glupol),pq(0,1),jqghq(0,-1,glupol))

         call curr6(1,braketg(1,-1,2,glupol),pq(0,2),
     &              psi(1,-1,1),p(0,1),jqhgq(0,-1,glupol))

         do mu = 0,5                                               
            jqqwh(mu,isig,glupol) = jqghq(mu,isig,glupol) + 
     1                             jqhgq(mu,isig,glupol)
c                            = (<2|gam.mu|q,1>+<2,q|gam.mu|1>)_l,isig)
         enddo
      enddo

      do glupol=1,2
         isig = -1         ! not needed in this process

c         if (ldebug)  print*, "wtowh", wtowh(0:5,1)
c         if (ldebug)  print*, "jqqwh", jqqwh(0:5,isig,glupol)
         mw(isig,glupol,1) = -dotcc(wtowh(0,1),jqqwh(0,isig,glupol))
         mat(1,isig,glupol,4)= mw(isig,glupol,1)*clr(3,3,isig)
         mat(1,isig,glupol,9)= mat(1,isig,glupol,4)
c         if (ldebug) print*, "jqqwh ", jqqwh(0,isig,glupol)
c         if (ldebug) print*, "mat ", mat(1,isig,glupol,9)
      enddo   !end gluon polarization l
      
      if (ldebug) then
         print*, "qqWHj"
         print*, psi(1:2,-1,1)
         print*, psi(1:2,-1,2)
c        print*, pq(0:4,2)
         print*, jqghq(0:5,-1,1)
         print*, jqghq(0:5,-1,2)
         print*, jqhgq(0:5,-1,1)
         print*, jqhgq(0:5,-1,2)
         print*, wtowh(0:3,1)
         print*, mat(1,-1,1:2,4)
      endif


c-------- end of S-CHANNEL computation ---------------------------

c------------------------------------------------------
c comparison with madgraph

      if (ldebugm) then

#ifdef WITH_MADGRAPH
         if (bos.eq.36) then
            call sud_wphgL(pbar,sign,udwhg)          !madgraph
c            print*, "mg ", udwhg
c            print*, "vb ", mat(1,-1,1,9)
c            print*, "vb ", mat(1,-1,2,9)
            ! does not yet exist:
         else if (bos.eq.46) then
            call sdu_wmhgL(pbar,sign,udwhg)          !madgraph
         endif

         if (ldebug) then
            print*, "psi 1 -1 1: ", psi(1:2,-1,1)
            print*, "psi 1 -1 2: ", psi(1:2,-1,2)

            do i=1,5 
               print*, p(0:3,i)
            enddo
            do i= 1,5 
               print*, pbar(0:3,i)
            enddo
            print*, qbar
         endif

c      endif
c      if (.False.) then ! comment second half of madgraph tests 


         isig=-1

         do i = 1,2                
            m1kb(1,isig,i) = 0d0 ! kartesian basis glupol = 1
            m2kb(1,isig,i) = 0d0 ! kartesian basis glupol = 2

            do lh = -1,1,2
               madgr(i,isig,lh,1) = 0d0
            enddo                
         enddo


         do glupol = -1,1,2
            do i =1,2                    
               ! nonzero for i=1,2; isig = -1, glupol = +1,-1
               madgr(i,isig,glupol,1) = tampr(i,isig,glupol,1)
               if (ldebug) print*, "i, l, tampr", i, glupol, tampr(i,isig,glupol,1)
            enddo
         enddo


         do i=1,2               
            m1kb(1,isig,i) = (madgr(i,isig,-1,1) - 
     &                 madgr(i,isig,1,1))/sqrt(2d0)       
            m2kb(1,isig,i) = -im*(madgr(i,isig,-1,1) + 
     &                 madgr(i,isig,1,1))/sqrt(2d0)                 
         enddo

c print results from checks:

         i=2 
         !do i=1,2                 
            temp = dreal(mat(1,isig,1,7+i))**2+dimag(mat(1,isig,1,7+i))**2
            temp1 = dreal(m1kb(1,isig,i))**2+dimag(m1kb(1,isig,i))**2
            if(abs(temp/temp1 -1).ge.1.0d-6) then
               print*,"mat(",1,",",isig,",1,",7+i,") = ",mat(1,isig,1,7+i)
               print*,"m1kb(",1,",",isig,",",i,") = ",m1kb(1,isig,i)
               print*,"temp/temp1 = ",temp/temp1
               print*, "matrix elements qqWHj are only within", abs(temp/temp1-1) 
c           else if (ldebug) then
c              print*, "matrix elements qqWHj are within", abs(temp/temp1-1) 
            endif      

            temp = dreal(mat(1,isig,2,7+i))**2+dimag(mat(1,isig,2,7+i))**2
            temp1 = dreal(m2kb(1,isig,i))**2+dimag(m2kb(1,isig,i))**2
            if(abs(temp/temp1 -1).ge.1.0d-6) then
               print*,"mat(",1,",",isig,",2,",7+i,") = ",mat(1,isig,2,7+i)
               print*,"m2kb(",1,",",isig,",",i,") = ",m2kb(1,isig,i)
               print*,"temp/temp1 = ",temp/temp1  
               print*, "matrix elements qqWHj are only within", abs(temp/temp1-1) 
c           else if (ldebug) then
c              print*, "matrix elements qqWHj are within", abs(temp/temp1-1) 
            endif 
         !enddo
#endif  
      endif ! end of debugging block for madgraph

c------------------------------------------------------


c-----------------END of calculation of matrix elements------------------


c sum the graphs, square them and map them onto udwa 

      res(1) = 0
      do isig = -1,-1         ! not needed in this process
         do glupol = 1,2
            mm(1,isig,glupol) = 0
            do i = 1,7
               mm(1,isig,glupol) = mm(1,isig,glupol) +  mat(1,isig,glupol,i)  

            enddo


            res(1)= res(1) + dreal(mm(1,isig,glupol))**2
     &                        + dimag(mm(1,isig,glupol))**2

         enddo
      enddo
      res(1) = res(1)*4d0*fpials(2)  ! C_2*3 is the color factor


      udwh = res(1)
c      if (ldebug) print*, "res ", res(1)

      return
      end


