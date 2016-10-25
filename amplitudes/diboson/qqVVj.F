****************subroutine QQVVj *************************************
c
c	Vera Hankele, <vera@particle.physik.uni-karlsruhe.de>
c	Initial version:  2006 April
c	Last modified: 
c
c  QQWWj calculates the matrix elements**2 for qq -> WW g -> ve e+ mu- vm~ g
c
c        q1 q2 ----> W+ W- g,   W+ ---> f4-bar f3, W- ---> f6-bar f5
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
c      return uucc = |M_real|^2   etc.
c
c	fpials is attached only in the end of the code
c
c index j = 2:3 indicates, whether g is emitted from 
c		upper line (j=2) or lower line (j=3)
c	l is the gluon polarization in the kartesian basis (l=1,2)
c		l=0 stands for building blocks without gluon emission
c	k is the process ID (1:uuww,2:ddww)
c	isig1/isig3 are the helicities of parton1/2 
c
c*********************************************************************

      subroutine qqVVj(pbar,sign, qbar,gsign, uuww, ddww)

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

      double precision  pbar(0:3,2+n_v), qbar(0:4), fac, musq, q(0:3)
      double precision  uuww, ddww, res(2),resv(2)
      double precision  p(0:3,2+n_v), pq(0:4,2),p21(0:4),
     1                  pww(0:4), fqp(0:4,2), fqm(0:4,2)
      double precision pgWP(0:4,2), pgWM(0:4,2), pWP(0:4,2), pWM(0:4,2)

      integer  sign(2+n_v), gsign, mu, i, jj, k, kk, id, isig
      integer  ifl(2,2), js1, js3, is
      integer  l   ! gluon polariz. (l=0:no g, l=1,2:g in kartesian basis)
      integer jmin, jmax
      logical jlog2,jlog3
      double complex prop21(4,2:3), prop43(4,2:3),
     1		     prop_pm(4,2:3),prop_mp(4,2:3)


      double complex mat(2,-1:1,0:2,9),temp, temp1
      double complex mata(2,-1:1,0:2,9),
     1 		     matz(2,-1:1,0:2,9)
      double complex mm(2,-1:1,2),m5(3,3:4,2:3,2)
      double complex ga,gb,gc,gd
      double complex maa, maz, mza, mzz, mww5, mww6
      double complex  mpm(4,2:3,2), mmp(4,2:3,2)
      double complex  m1u,m2u,m1l,m2l,z1u,z1l
      double complex  m1(-1:1,2,3), m2(-1:1,2,3)
      double complex  matot,ratot,raz,raa! for checks only
      double precision eps(0:3,2) ! g in kartesian basis 
      double complex psi(2,-1:1,4), psiWM(2,-1:1,2), psiWP(2,-1:1,2)
      double complex psiWMg(2,-1:1, 2,2), psiWPg(2,-1:1, 2,2)
      double complex jqq(0:5,-1:1,2,0:2), braketg(2,-1:1,2,2)
      double complex jqgaq(0:5,-1:1,2), jqagq(0:5,-1:1,2), 
     &               jqqa(0:5,-1:1,2), jh1(0:5,-1:1), jh2(0:5,-1:1)
      double complex jqqgWP(0:5,-1:1,2,0:3), jqqgWM(0:5,-1:1,2,0:3)
      double complex jwp(0:5,-1:1,-1:1,4,0:2), jwm(0:5,-1:1,-1:1,4,0:2)
      double precision fq(0:4,4),dummy(0:4),
     1 		       bq(0:4,4)
      double complex bkjqq(2,-1:1,-1:1,4,2:3,0:2),
     1		     bkjqqg(2,-1:1,-1:1,4,2:3,0:2),
     1 		     gbkjqq(2,-1:1,-1:1,4,2:3,0:2)
      double complex braketgWP(2,-1:1,2,2),braketgWM(2,-1:1,2,2)         
      double precision  pwpg(0:4,4),pwmg(0:4,4)
      double precision  pga(0:4,4),pgz(0:4,4),
     3			pag(0:4,4),pzg(0:4,4)
      double complex psia(2,-1:1,-1:1,4), psiz(2,-1:1,-1:1,4)
      double complex ja(0:5,-1:1,-1:1,4,0:2), jz(0:5,-1:1,-1:1,4,0:2),
     1		     jag(0:5,-1:1),jga(0:5,-1:1),jg0(0:5,-1:1),
     1		     jgz(0:5,-1:1),jzg(0:5,-1:1),jgm(0:5,-1:1)
      double complex jwgi(0:5,-1:1),jgwi(0:5,-1:1),jwg0i(0:5,-1:1),
     1		     jwgii(0:5,-1:1),jgwii(0:5,-1:1),jwg0ii(0:5,-1:1)
      double complex epsa0(0:3),epsz0(0:3),epsa(0:3),epsz(0:3),
     1               epsNCwp(0:5,-1:1,3:4,2,2:3,0:2),
     2		     epsNCwm(0:5,-1:1,3:4,2,2:3,0:2)
      double complex qepswp(3:4,2:3,0:2),qepswm(3:4,2:3,0:2)
      double complex jj21m(0:2),jj21p(0:2),jj43m(0:2),jj43p(0:2)
      double complex zp(4:5,2:3),zm(4:5,2:3)
      double complex zpm(2:3,2),zmp(2:3,2)
      double complex ma(-1:1,2,0:2), mz(-1:1,2,0:2)
      double complex propt(-1:1,-1:1,6,2:3,2)
      double complex  zm2i(2:4)
      double complex contract_Tjj,! contract_T1j,contract_T2j,
     1 		     dotcc, dotrc, dotqj, s1c
      external contract_Tjj, !contract_T1j,contract_T2j,
     1 	       dotcc, dotrc, dotqj, s1c
      integer lh
      double complex madgr(9,-1:1,-1:1,2),m1kb(2,-1:1,9),m2kb(2,-1:1,9) !checks
      double complex m1hj(-1:1,2:3),m1kj(2,2:3)
      double complex m1haj(-1:1,2:3),m1hzj(-1:1,2:3),
     1 		     m1hu(-1:1,2:3),m1hl(-1:1,2:3),m1hc(-1:1),
     1 		     m1hua(-1:1,2:3),m1hla(-1:1,2:3),
     1 		     m1huz(-1:1,2:3),m1hlz(-1:1,2:3),
     1 		     m1kaj(2,2:3),m1kzj(2,2:3),
     1 		     m1ku(2,2:3),m1kl(2,2:3),m1kc(2),
     1 		     m1kua(2,2:3),m1kla(2,2:3),m1kc2(2),m1kc3(2),
     1 		     m1kuz(2,2:3),m1klz(2,2:3),m1hc2(-1:1),m1hc3(-1:1)
      double complex im
      parameter (im = (0d0,1d0))
      logical ldebug, linit, lerror, ldebugm
      integer n,m
     
      save ifl, zm2i, linit
	
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

         k = 1                  ! uuww
         ifl(1,k) = 3
         ifl(2,k) = 3
         k = 2                  ! ddww
         ifl(1,k) = 4
         ifl(2,k) = 4

c         zm2i(2) = 1/dcmplx(xm2(2),-xmg(2))
c         zm2i(3) = 1/dcmplx(xm2(3),-xmg(3))
c         zm2i(4) = 1/dcmplx(xm2(4),-xmg(4))

      endif ! linit

      do k = 1,2
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

         do isig = -1,1,2	! fermion helicity 

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

      call ket2c(braketg(1,-1,1,l),.false.,pq(0,1),-1,
     &        qp(0,1),wp(1,1),braketgWP(1,-1,1,l),pgWP(0,1)) !|W+,g,1>_l,-1

      call ket2c(braketg(1,-1,1,l),.false.,pq(0,1),-1,
     &        qm(0,1),wm(1,1),braketgWM(1,-1,1,l),pgWM(0,1)) !|W-,g,1>_l,-1

      call curr6(-1,psi(1,-1,2),p(0,2),
     &           braketgWP(1,-1,1,l),pgWP(0,1),jqqgWP(0,-1,l,1))

      call curr6(-1,psi(1,-1,2),p(0,2),
     &           braketgWM(1,-1,1,l),pgWM(0,1),jqqgWM(0,-1,l,1))

      m1(-1,l,1) = dotcc(jqqgWP(0,-1,l,1),wm(1,1))
      m2(-1,l,1) = dotcc(jqqgWM(0,-1,l,1),wp(1,1))

      mat(1,-1,l,1) = clr(3,3,-1)**2 * m1(-1,l,1)
      mat(2,-1,l,1) = clr(3,3,-1)**2 * m2(-1,l,1)

      enddo

c             - emission from lower line

      do l=1,2

      call bra2c(braketg(1,-1,2,l),.false.,pq(0,2),-1,
     &       qm(0,1),wm(1,1),braketgWM(1,-1,2,l),pgWM(0,2)) !<2,g,W-|_l,-1

      call bra2c(braketg(1,-1,2,l),.false.,pq(0,2),-1,
     &       qp(0,1),wp(1,1),braketgWP(1,-1,2,l),pgWP(0,2)) !<2,g,W+|_l,-1

      call curr6(1,braketgWM(1,-1,2,l),pgWM(0,2),
     &            psi(1,-1,1),p(0,1),jqqgWM(0,-1,l,2))

      call curr6(1,braketgWP(1,-1,2,l),pgWP(0,2),
     &            psi(1,-1,1),p(0,1),jqqgWP(0,-1,l,2))

      m1(-1,l,2) = dotcc(jqqgWM(0,-1,l,2),wp(1,1))
      m2(-1,l,2) = dotcc(jqqgWP(0,-1,l,2),wm(1,1))

      mat(1,-1,l,2) = clr(3,3,-1)**2 * m1(-1,l,2)
      mat(2,-1,l,2) = clr(3,3,-1)**2 * m2(-1,l,2)

      enddo

c             - emission from intermediate quark line

      do l=1,2

      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,
     &       qm(0,1),wm(1,1),psiWM(1,-1,2),pWM(0,2)) !<2,W-,g|_-1

      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,
     &       qp(0,1),wp(1,1),psiWP(1,-1,2),pWP(0,2)) !<2,W+,g|_-1

      call bra2r(psiWM(1,-1,2),.false.,pWM(0,2),-1,
     &            q,eps(0,l),psiWMg(1,-1,2,l),pWMg(0,2)) !<2,W-,g|_l,-1

      call bra2r(psiWP(1,-1,2),.false.,pWP(0,2),-1,
     &            q,eps(0,l),psiWPg(1,-1,2,l),pWPg(0,2)) !<2,W+,g|_l,-1

      call curr6(1,psiWMg(1,-1,2,l),pWMg(0,2),
     &            psi(1,-1,1),p(0,1),jqqgWM(0,-1,l,3))

      call curr6(1,psiWPg(1,-1,2,l),pWPg(0,2),
     &            psi(1,-1,1),p(0,1),jqqgWP(0,-1,l,3))

      m1(-1,l,3) = dotcc(jqqgWM(0,-1,l,3),wp(1,1))
      m2(-1,l,3) = dotcc(jqqgWP(0,-1,l,3),wm(1,1))

      mat(1,-1,l,3) = clr(3,3,-1)**2 * m1(-1,l,3)
      mat(2,-1,l,3) = clr(3,3,-1)**2 * m2(-1,l,3)

      enddo

      do l=1,2

            mat(1,-1,l,8) = mat(1,-1,l,1) + mat(1,-1,l,2)
     &                      + mat(1,-1,l,3)

            mat(2,-1,l,8) = mat(2,-1,l,1) + mat(2,-1,l,2)
     &                      + mat(2,-1,l,3)

      enddo

c------- end of T-CHANNEL computation

c-------  S-CHANNEL ----------------

c     Get the f-fbar currents with one gluon radiated from the
c     current line.  There are two terms, one for gluon emission to
c     either side of the ffV vertex:

      do l=1,2
         do isig = -1,1,2

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
         do isig=-1,1,2

            ma(isig,l,1) = -dotcc(aww(0,1),jqgaq(0,isig,l))
            ma(isig,l,2) = -dotcc(aww(0,1),jqagq(0,isig,l))
            ma(isig,l,0) = -dotcc(aww(0,1),jqqa(0,isig,l))

            mz(isig,l,1) = -dotcc(zww(0,1),jqgaq(0,isig,l))
            mz(isig,l,2) = -dotcc(zww(0,1),jqagq(0,isig,l))
            mz(isig,l,0) = -dotcc(zww(0,1),jqqa(0,isig,l))

            do k = 1,2          !k=1: uuww, k=2: ddww
               mat(k,isig,l,4) = ma(isig,l,1)*clr(ifl(1,k),1,isig)
               mat(k,isig,l,5) = ma(isig,l,2)*clr(ifl(1,k),1,isig)
               mat(k,isig,l,6) = mz(isig,l,1)*clr(ifl(1,k),2,isig)
               mat(k,isig,l,7) = mz(isig,l,2)*clr(ifl(1,k),2,isig)
               mat(k,isig,l,9) = mat(k,isig,l,4) + mat(k,isig,l,5)
     &                         + mat(k,isig,l,6) + mat(k,isig,l,7)
            enddo
         enddo
      enddo   !end gluon polarization l

c-------- end of S-CHANNEL computation ---------------------------

#ifdef WITH_MADGRAPH
c------------------------------------------------------
c comparison with madgraph

      if (ldebugm) then

         do isig = -1,1,2
            do i = 1,9
               do k = 1,2

                  m1kb(k,isig,i) = 0d0 ! kartesian basis l = 1
                  m2kb(k,isig,i) = 0d0 ! kartesian basis l = 2

                  do lh = -1,1,2
                     madgr(i,isig,lh,k) = 0d0
                  enddo

               enddo
            enddo
         enddo

c compute m1/m2 in kartesian basis:

         do isig = -1,1,2
            do lh = -1,1,2
               do k = 1,2

         madgr(1,isig,lh,k) = tampr(1,isig,lh,k)
         madgr(2,isig,lh,k) = tampr(2,isig,lh,k)
         madgr(3,isig,lh,k) = tampr(3,isig,lh,k)
         madgr(4,isig,lh,k) = tampr(4,isig,lh,k)
     &                 + tampr(10,isig,lh,k) + tampr(14,isig,lh,k)
         madgr(5,isig,lh,k) = tampr(6,isig,lh,k)
     &                 + tampr(12,isig,lh,k) + tampr(16,isig,lh,k)
         madgr(6,isig,lh,k) = tampr(5,isig,lh,k)
     &                 + tampr(8,isig,lh,k) + tampr(11,isig,lh,k)
     &                 + tampr(15,isig,lh,k) + tampr(18,isig,lh,k)
         madgr(7,isig,lh,k) = tampr(7,isig,lh,k)
     &                 + tampr(9,isig,lh,k) + tampr(13,isig,lh,k)
     &                 + tampr(17,isig,lh,k) + tampr(19,isig,lh,k)
         madgr(8,isig,lh,k) = madgr(1,isig,lh,k)
     &                 + madgr(2,isig,lh,k) + madgr(3,isig,lh,k)
         madgr(9,isig,lh,k) = madgr(4,isig,lh,k) + madgr(5,isig,lh,k)
     &                 + madgr(6,isig,lh,k) + madgr(7,isig,lh,k)

               enddo
            enddo
         enddo


         do isig = -1,1,2
            do i=1,9
               do k=1,2
               
	 m1kb(k,isig,i) = (madgr(i,isig,-1,k) - 
     &            madgr(i,isig,1,k))/sqrt(2d0)       
	 m2kb(k,isig,i) = -im*(madgr(i,isig,-1,k) + 
     &            madgr(i,isig,1,k))/sqrt(2d0)

               enddo
            enddo
         enddo


c print results from checks:
         do isig = -1,1,2
            do k = 2,2
               do i = 4,7

    	    print*,"ratio for l=1, isig = ",isig,"i = ",i,":",
     &	    	mat(k,isig,1,i)/m1kb(k,isig,i)
            print*,"mat(",k,",",isig,",1,",i,") = ",mat(k,isig,1,i)
            temp = dreal(mat(k,isig,1,i))**2+dimag(mat(k,isig,1,i))**2
            print*,temp
            print*,"m1kb(",k,",",isig,",",i,") = ",m1kb(k,isig,i)
            temp1 = dreal(m1kb(k,isig,i))**2+dimag(m1kb(k,isig,i))**2
            print*,temp1
            print*,"temp/temp1 = ",temp/temp1
            print*," "

    	    print*,"ratio for l=2, isig = ",isig,"i = ",i,":",
     &	    	mat(k,isig,2,i)/m2kb(k,isig,i)
            print*,"mat(",k,",",isig,",2,",i,") = ",mat(k,isig,2,i)
            temp = dreal(mat(k,isig,2,i))**2+dimag(mat(k,isig,2,i))**2
            print*,temp
            print*,"m2kb(",k,",",isig,",",i,") = ",m2kb(k,isig,i)
            temp1 = dreal(m2kb(k,isig,i))**2+dimag(m2kb(k,isig,i))**2
            print*,temp1
            print*,"temp/temp1 = ",temp/temp1
            print*," "

               enddo
            enddo
         enddo

      endif ! end of debugging block
c------------------------------------------------------
#endif /* MADGRAPH */

c-----------------END of calculation of matrix elements------------------


c sum the graphs, square them and map them onto uuww and ddww
      do k = 1,2
         res(k) = 0
         do isig = -1,1,2
 	    do l = 1,2
               mm(k,isig,l) = 0
               do i = 1,7
                  mm(k,isig,l) = mm(k,isig,l) +  mat(k,isig,l,i)
               enddo

               res(k)= res(k) + dreal(mm(k,isig,l))**2
     &                        + dimag(mm(k,isig,l))**2

            enddo
         enddo
         res(k) = res(k)*4d0*fpials(2)  ! C_2*3 is the color factor     
      enddo

      uuww = res(1)
      ddww = res(2)

      return
      end


