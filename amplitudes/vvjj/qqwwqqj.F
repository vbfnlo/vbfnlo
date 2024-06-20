c---------------------------------------------------------------------------
c baj:  modify Dieter"s version from Apr. 12, 2005 such that 
c       virtual contributions are removed, but real emission terms included
c---------------------------------------------------------------------------

      subroutine qqwwqqj(pbar,sign,qbar,gsign,
     1                  uucc,uuss,ddcc,ddss,udsc,ducs)
          use globalvars, only: ldoblha
      implicit none
c
c	Barbara Jaeger, <baj@particle.physik.uni-karlsruhe.de>
c	Initial version:  2005 May 09
c	Last modified:    2005 
C
C  QQWWQQJ calculates the matrix elements**2 for electroweak
c  weak boson pair production by quark quark scattering
C
C        q1 q3    ---->   q2 q4 g W+W-,   W+ ---> f5-bar f6, W- ---> f7-bar f8
C
C  and crossing related processes. Pauli interference terms for
c  identical fermions are neglected. In particular, only the
c  t-channel exchange of elctroweak bosons is considered. s-channel
c  production of 3 weak bosons is NOT implemented.
c  Gluon polarization vectors are given in the kartesian basis 
c  (comparison with madgraph amplitudes tampr requires to express these in 
c   kart. basis also, even though they are given in the helicity basis a priori). 
c
C  This code includes only real emission contributions, i.e.
c
c      return uucc = |M_real|^2   etc.
c
c  debugging: included for NC (works for incoming g only) 
c	      and CC (arbitrary crossing)
c
c	fpials is attached only in the end of the code
c
c index j = 2:3 indicates, whether g is emitted from 
c		upper (12) line (j=2) or lower (34) line (j=3)
c	l is the gluon polarization in the kartesian basis (l=1,2)
c		l=0 stands for building blocks without gluon emission
c	k is the process ID (1:uucc,2:uuss,3:ddcc;4:ddss;5:udsc;6:ducs)
c	isig1/isig3 are the helicities of parton1/3 
c
c
c---------------------------------------------------------------------
c
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/an_couplings.inc"
c
c electroweak couplings are taken from KOPPLN
c
      double precision  clr, xm2, xmg, b, v, a
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     1                  V(4,5),A(4,5)
c
c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"
c
c variables for the main part of the program
c
      real*8 fpials(2:3), fpi
      parameter (fpi=4d0*pi)

      double precision  pbar(0:3,4+max_v),qbar(0:4), fac
      double precision uucc(3),uuss(3),ddcc(3),ddss(3),udsc(3),ducs(3),
     1 		       res(6,2:3),resv(6)
      double precision  p(0:3,4+max_v),q(0:4), p21(0:4,2:3), p43(0:4,2:3),
     1                  pq(0:4,4),pwp(0:4),pwm(0:4),pww(0:4),
     2			qpm(0:4,2:3),qmp(0:4,2:3)!,
c     2			qpml(0:4,2:3),qmpl(0:4,2:3)
      integer  sign(4+max_v),gsign, mu, i, j, k, kk,
     1         isig, isig1, isig3,is,ka,kb,kf
      integer  ifl(4,6) !, is1, is3 taken from "tensor.inc"
      integer  l   ! gluon polariz. (l=0:no g, l=1,2:g in kartesian basis)
      integer jmin, jmax
      logical jlog2,jlog3
      double complex prop21(4,2:3), prop43(4,2:3),
     1		     prop_pm(4,2:3),prop_mp(4,2:3)
      double complex mat(6,-1:1,-1:1,2:3,0:2,9),temp
      double complex mata(6,-1:1,-1:1,2:3,0:2,9),
     1 		     matz(6,-1:1,-1:1,2:3,0:2,9)
      double complex mm(6,-1:1,-1:1,2:3,2),m5(3,3:4,2:3,2)
      double complex ga,gb,gc,gd
      double complex maa, maz, mza, mzz, mww5, mww6
      double complex  mpm(4,2:3,2), mmp(4,2:3,2)
      double complex  m1u,m2u,m1l,m2l,z1u,z1l
c      double complex  m1,m2 ! for checks only
      double complex  matot,ratot,raz,raa! for checks only
      double precision eps(0:3,2) ! g in kartesian basis 
      double complex psi(2,-1:1,4),jqq(0:5,-1:1,2,-1:1,0:2), 
     1 		     braketg(2,-1:1,4,2), jh1(0:5,-1:1), jh2(0:5,-1:1)
      double complex psiwp(2,4), psiwm(2,4),
     1               jwp(0:5,-1:1,-1:1,4,0:2), jwm(0:5,-1:1,-1:1,4,0:2)
      double precision fqp(0:4,4), fqm(0:4,4), fq(0:4,4),dummy(0:4),
     1 		       bq(0:4,4)
      double complex bkjqq(2,-1:1,-1:1,4,2:3,0:2),
     1		     bkjqqg(2,-1:1,-1:1,4,2:3,0:2),
     1 		     gbkjqq(2,-1:1,-1:1,4,2:3,0:2)
      double complex braketgWP(2,-1:1,4,0:2),braketWPg(2,-1:1,4,0:2),
     3		     braketgWM(2,-1:1,4,0:2),braketWMg(2,-1:1,4,0:2)      
       double complex braketgA(2,-1:1,0:2),braketAg(2,-1:1,0:2),
     3		      braketgZ(2,-1:1,0:2),braketZg(2,-1:1,0:2)      
      double precision  pgwp(0:4,4),pgwm(0:4,4),
     3			pwpg(0:4,4),pwmg(0:4,4)
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
     2		     epsNCwm(0:5,-1:1,3:4,2,2:3,0:2),
     2               epsCCwp(0:5,-1:1,3:4,2,2:3,0:2), 
     2		     epsCCwm(0:5,-1:1,3:4,2,2:3,0:2)
      double complex qepsCCwp(2,2:3,0:2),qepsCCwm(2,2:3,0:2),
     1	 	     qepswp(3:4,2:3,0:2),qepswm(3:4,2:3,0:2)
      double complex jj21m(0:2),jj21p(0:2),jj43m(0:2),jj43p(0:2)
      double complex zp(4:5,2:3),zm(4:5,2:3)
      double complex zpm(2:3,2),zmp(2:3,2)
      double complex  ma(4,2:3,0:2), mz(4,2:3,0:2)
      double complex propt(-1:1,-1:1,6,2:3,2)
      double complex  zm2i(2:4)
      double complex contract_Tjj,! contract_T1j,contract_T2j,
     1 		     dotcc, dotrc, s1c !, dotqj
      external contract_Tjj, !contract_T1j,contract_T2j,
     1 	       dotcc, dotrc, s1c !, dotqj
      integer lh
      double complex m1hj(-1:1,2:3),m1kj(2,2:3)
      double complex m1haj(-1:1,2:3),m1hzj(-1:1,2:3),
     1 		     m1hu(-1:1,2:3),m1hl(-1:1,2:3),m1hc(-1:1),
     1 		     m1hua(-1:1,2:3),m1hla(-1:1,2:3),
     1 		     m1huz(-1:1,2:3),m1hlz(-1:1,2:3),
     1 		     m1kaj(2,2:3),m1kzj(2,2:3),
     1 		     m1ku(2,2:3),m1kl(2,2:3),m1kc(2),
     1 		     m1kua(2,2:3),m1kla(2,2:3),
     1 		     m1kuz(2,2:3),m1klz(2,2:3)
      double complex im
      parameter (im = (0d0,1d0))
      logical ldebug, linit,lerror,ldebugm
      integer n,m
     
      save ifl, zm2i, linit, lerror,ldebug,ldebugm
c stuff for comparison with madgraph :
#include "tampr-lists.inc" 
	
      data linit /.true./, lerror /.true./, ldebug /.false./, ldebugm /.false./

      real*8 widthfactor

      widthfactor=1.0d0
c      if (ofs.and.with_anom) then
c         widthfactor=0.d0
c      endif
c
c  ---------------------------------------------------------------------
c 
c initialize & precompute stuff needed below:
c

c fix strong coupling gs**2 for the 2 quark lines:

      fpials(2) = fpi*als(1,1)
      fpials(3) = fpi*als(2,1)
	  
c define flavors of external quarks for the 4 NC and 2 CC subprocesses
c
      if (linit) then
         linit = .false.
         k = 1                  ! uucc
         ifl(1,k) = 3
         ifl(2,k) = 3
         ifl(3,k) = 3
         ifl(4,k) = 3
         k = 2                  ! uuss
         ifl(1,k) = 3
         ifl(2,k) = 3
         ifl(3,k) = 4
         ifl(4,k) = 4
         k = 3                  ! ddcc
         ifl(1,k) = 4
         ifl(2,k) = 4
         ifl(3,k) = 3
         ifl(4,k) = 3
         k = 4                  ! ddss
         ifl(1,k) = 4
         ifl(2,k) = 4
         ifl(3,k) = 4
         ifl(4,k) = 4
         k = 5                  ! udsc
         ifl(1,k) = 3
         ifl(2,k) = 4
         ifl(3,k) = 4
         ifl(4,k) = 3
         k =6                   ! ducs
         ifl(1,k) = 4
         ifl(2,k) = 3
         ifl(3,k) = 3
         ifl(4,k) = 4
	 
         zm2i(2) = 1/dcmplx(xm2(2),-xmg(2)*widthfactor)
         zm2i(3) = 1/dcmplx(xm2(3),-xmg(3)*widthfactor)
         zm2i(4) = 1/dcmplx(xm2(4),-xmg(4)*widthfactor)
      endif

 701  if (ldebug) then
         print*," " 
         print*," enter qqwwqq "
      endif
      
      if (gsign.eq.1) then	!final state gluon
        jlog2 = .true.		! can couple to upper/lower line
	jlog3 = .true. 
      	jmin = 2
	jmax = 3
      else		       !initial state gluon -> only:
	 if (sign(1).ne.sign(2)) then  !gluon from upper line
	       jlog2 = .true.
	       jlog3 = .false. 
	       jmin = 2
	       jmax = 2
	 else			       !gluon from lower line
	       jlog2 = .false.
	       jlog3 = .true. 
	       jmin = 3
	       jmax = 3
       endif
      endif	
                  
      do k = 1,6
         do isig1 = -1,1,2
            do isig3 = -1,1,2
	       do j = 2,3
	          do l = 0,2
                     do i = 1,9
                  	mat(k,isig1,isig3,j,l,i)  = 0
                  	mata(k,isig1,isig3,j,l,i) = 0
                  	matz(k,isig1,isig3,j,l,i) = 0
     		     enddo
		  enddo	
               enddo
            enddo
         enddo
      enddo
c
c identify fermion line sign factors (for 1 3 -> 2 4 etc.)
c
c      is1 = sign(1)
c      is3 = sign(3)
c      js1 = (3+sign(1))/2       ! 1 for sign1=-1,2 for sign1=+1
c      js3 = (7+sign(3))/2       ! 3 for sign3=-1,4 for sign3=+1

c fix is1 such that is1 = +1 for   q1 ->  q2 g 
c		    is1 = -1 for  q1b -> q2b g 
c		    is1 =  0 for    g -> q1b q2
c
c      for is3:     is3 = +1 for   q3 ->  q4 g 
c		    is3 = -1 for  q3b -> q4b g 
c		    is3 =  0 for    g -> q3b q4

	is1 = (sign(1)+sign(2))/2  
	is3 = (sign(3)+sign(4))/2
		
c (is1,is3 are fixed here and don"t change throughout this run of the program)


c define the internal momenta
c
      do mu = 0,3
         do i = 1,4+n_v
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo
	 q(mu) = qbar(mu)*gsign
	 
         p21(mu,3) = p(mu,2) - p(mu,1)	! g from 34 line
         p21(mu,2) = p21(mu,3) + q(mu) 	! g from 12 line
         p43(mu,2) = p(mu,4) - p(mu,3)
         p43(mu,3) = p43(mu,2) + q(mu)	 
	 
         pwp(mu) = p(mu,6) - p(mu,5)
         pwm(mu) = p(mu,8) - p(mu,7)
         pww(mu) = pwp(mu) + pwm(mu)
	 
      enddo
	
      do j= 2,3
        p21(4,j) = p21(0,j)**2 - p21(1,j)**2 - p21(2,j)**2 - p21(3,j)**2
        p43(4,j) = p43(0,j)**2 - p43(1,j)**2 - p43(2,j)**2 - p43(3,j)**2     
      enddo
	
      q(4)   = 0d0
      pwp(4) = pwp(0)**2 - pwp(1)**2 - pwp(2)**2 - pwp(3)**2
      pwm(4) = pwm(0)**2 - pwm(1)**2 - pwm(2)**2 - pwm(3)**2
      pww(4) = pww(0)**2 - pww(1)**2 - pww(2)**2 - pww(3)**2

c  ---------------------------------------------------------------------
c
c get the vector boson propagator factors
c
c depending on value of j, gluon is attached to respective quark line or not;
c no V is attached here 

      do j = 2,3	
       prop21(1,j) = 1/p21(4,j)		!prop21(bos=1:4,j=2:3)
c        if (p21(4,j).le.0d0) then
c           prop21(2,j) = 1/(p21(4,j)-xm2(2))
c           prop21(3,j) = 1/(p21(4,j)-xm2(3))
c        else
           prop21(2,j) = 1/dcmplx(p21(4,j)-xm2(2),xmg(2)*widthfactor)
           prop21(3,j) = 1/dcmplx(p21(4,j)-xm2(3),xmg(3)*widthfactor)
c        endif
       prop21(4,j) = prop21(3,j)

       prop43(1,j) = 1/p43(4,j)
c        if (p43(4,j).le.0d0) then
c           prop43(2,j) = 1/(p43(4,j)-xm2(2))
c           prop43(3,j) = 1/(p43(4,j)-xm2(3))
c        else
           prop43(2,j) = 1/dcmplx(p43(4,j)-xm2(2),xmg(2)*widthfactor)
           prop43(3,j) = 1/dcmplx(p43(4,j)-xm2(3),xmg(3)*widthfactor)
c        endif
       prop43(4,j) = prop43(3,j)
      enddo
c
c  ---------------------------------------------------------------------

c for box-box and BV graphs we need the propagators for t-channel bosons between quark lines
c as seen from upper line these W momenta are INCOMING. They are OUTGOINg as seen from lower 
c line
      
      do j = 2,3 ! g from upper/lower line		
      
      do mu = 0,3
         qpm(mu,j) = pwp(mu)+p21(mu,j)  !W+ emitted on upper line
         qmp(mu,j) = pwm(mu)+p21(mu,j)  !W- emitted on upper line
      enddo
      
      qpm(4,j) = qpm(0,j)**2-qpm(1,j)**2-qpm(2,j)**2-qpm(3,j)**2
      qmp(4,j) = qmp(0,j)**2-qmp(1,j)**2-qmp(2,j)**2-qmp(3,j)**2

      prop_pm(1,j) = clr(3,3,-1)**2/qpm(4,j)
      prop_pm(2,j) = clr(3,3,-1)**2/dcmplx(qpm(4,j)-xm2(2),xmg(2)*widthfactor)
      prop_pm(3,j) = clr(3,3,-1)**4/dcmplx(qpm(4,j)-xm2(3),xmg(3)*widthfactor)
      prop_mp(1,j) = clr(3,3,-1)**2/qmp(4,j)
      prop_mp(2,j) = clr(3,3,-1)**2/dcmplx(qmp(4,j)-xm2(2),xmg(2)*widthfactor)
      prop_mp(3,j) = clr(3,3,-1)**4/dcmplx(qmp(4,j)-xm2(3),xmg(3)*widthfactor)
      
      prop_pm(4,j) = prop_pm(3,j)
      prop_mp(4,j) = prop_mp(3,j)
      
      enddo


c  ---------------------------------------------------------------------

c
c get the external quark spinors (including factor sqrt(2E) )
c
      call psi0m(4,pbar(0,1),sign(1),psi)

c     NOTE: psi(2,-1:1,4) is a two component spinor with helicity -1
c     or 1.  The last entry identifies the fermion.  If this entry is
c     odd psi is a ket, if even psi is a bra.
c     psi(1,isig1,1) = |1>_isig1
c     psi(1,isig3,3) = |3>_isig3
c     psi(1,isig1,2) = <2|_isig1
c     psi(1,isig3,4) = <4|_isig3
c
c get the f-fbar currents (with no gluon radiation) 
c	J21^mu=jqq(mu,*,1,is1,0), J43^mu=jqq(mu,*,2,is3,0) 
c
        call curr6(1,psi(1,-1,2),p(0,2),psi(1,-1,1),p(0,1),
     #						jqq(0,-1,1,is1,0))      
        call curr6(1,psi(1,-1,4),p(0,4),psi(1,-1,3),p(0,3),
     #						jqq(0,-1,2,is3,0))
c
c nomenclature: jqq(mu,		...	Lorentz index (0:3),(4:5) moment. info
c		    hel,	...	quark helicity (+-1)
c		    u/l,	...	upper(1)/lower(2) quark line
c		    is,		...	current for qq(g),qbqb(g),or gqbq line
c		     l)		...	g polarization (l=0:no g,l=1,2 with g)
c
c
c  Get the gluon polarization vector and the gluon emission spinors
      do l = 1,2	! 2 gluon polarizations
         call polvec(qbar,l,eps(0,l))  ! get gluon pol.vectors
	 	 
         do isig = -1,1,2	! fermion helicity 
 
c     NOTES for bras and kets: .true. if psi is a 2-spinor of the chi
c     form as output by psi0m, .false. otherwise.  the last entry is
c     the sum of the two momenta (p plus q) and effectively the
c     momentum of the new spinor.
            
            call ket2r(psi(1,isig,1),.true.,p(0,1),isig,q,eps(0,l),
     $           braketg(1,isig,1,l),pq(0,1))      ! |q,1>_l,isig1
            call bra2r(psi(1,isig,2),.true.,p(0,2),isig,q,eps(0,l),
     $           braketg(1,isig,2,l),pq(0,2))      ! <2,q|_l,isig2
            call ket2r(psi(1,isig,3),.true.,p(0,3),isig,q,eps(0,l),
     $           braketg(1,isig,3,l),pq(0,3))      ! |q,3>_l,isig3
            call bra2r(psi(1,isig,4),.true.,p(0,4),isig,q,eps(0,l),
     $           braketg(1,isig,4,l),pq(0,4))      ! <4,q|_l,isig4

         enddo
	          
c     braketg contains the free quark spinors multiplied by a fermion
c     propagator and a gluon eps_slash. 
c     NOTATION: braketg(2 component spinor, isig =-1 or 1 (fermion hel.),
c     fermion ID = 1:4, gluon polarization l=1:2)
 
      enddo
       
c     Get the f-fbar currents with one gluon radiated from the
c     current line.  There are two terms, one for gluon emission to
c     either side of the ffV vertex:
c
c	gluon from upper line:
      do l = 1, 2 ! gluon polarizations
         call curr6(1,psi(1,-1,2),p(0,2),braketg(1,-1,1,l),pq(0,1),jh1)	
c                                            =  <2|vertex|q,1>_l,isig1
         call curr6(1,braketg(1,-1,2,l),pq(0,2),psi(1,-1,1),p(0,1),jh2)	
c                                            =  <2,q|vertex|1>_l,isig1
         do isig = -1,1,2 ! fermion helicity
            do mu = 0,5
	       jqq(mu,isig,1,is1,l) = jh1(mu,isig) + jh2(mu,isig)
c                            = (<2|gam.mu|q,1>+<2,q|gam.mu|1>)_l,isig1
            enddo
         enddo
         
c	gluon from lower line:
         call curr6(1,psi(1,-1,4),p(0,4),braketg(1,-1,3,l),pq(0,3),jh1)	
c                                           =   <4|gam.mu|q,3>_l,isig3
         call curr6(1,braketg(1,-1,4,l),pq(0,4),psi(1,-1,3),p(0,3),jh2)	
c                                           =   <4,q|gam.mu|3>_l,isig3
         do isig = -1,1,2
            do mu = 0,5
               jqq(mu,isig,2,is3,l) = jh1(mu,isig) + jh2(mu,isig)
c                            = (<4|gam.mu|q,3>+<4,q|gam.mu|3>)_l,isig3
            enddo
         enddo
      enddo


c -----------------------------------------------------------------------
c -----------------------------------------------------------------------
c
c contract with vvtoww tensors to get Vertex-Vertex scattering diagrams
c
c neutral current first:
c
      do l = 1,2	! gluon polarization
        do isig1 = -1,1,2  ! fermion1 helicity
          do isig3 = -1,1,2  ! fermion3 helicity
	 
           if (jlog2) then
	   j = 2 ! g from upper line
	   
	   maa = contract_Tjj(aaww(0,0,j),jqq(0,isig1,1,is1,l),
     #	   					jqq(0,isig3,2,is3,0))
           maz = contract_Tjj(azww(0,0,j),jqq(0,isig1,1,is1,l),
     #	   					jqq(0,isig3,2,is3,0))
           mza = contract_Tjj(zaww(0,0,j),jqq(0,isig3,2,is3,0),
     #	   					jqq(0,isig1,1,is1,l))
           mzz = contract_Tjj(zzww(0,0,j),jqq(0,isig1,1,is1,l),
     #	   					jqq(0,isig3,2,is3,0))
 
           do k = 1,4
              mat(k,isig1,isig3,j,l,1) = 
     1    	   maa*clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3)
     2    	 + maz*clr(ifl(1,k),1,isig1)*clr(ifl(3,k),2,isig3)
     3    	 + mza*clr(ifl(1,k),2,isig1)*clr(ifl(3,k),1,isig3)
     4    	 + mzz*clr(ifl(1,k),2,isig1)*clr(ifl(3,k),2,isig3)
           enddo
	   endif

           if (jlog3) then
	   j = 3 ! g from lower line
	   
	   maa = contract_Tjj(aaww(0,0,j),jqq(0,isig1,1,is1,0),
     #	   					jqq(0,isig3,2,is3,l))
           maz = contract_Tjj(azww(0,0,j),jqq(0,isig1,1,is1,0),
     #	   					jqq(0,isig3,2,is3,l))
           mza = contract_Tjj(zaww(0,0,j),jqq(0,isig3,2,is3,l),
     #	   					jqq(0,isig1,1,is1,0))
           mzz = contract_Tjj(zzww(0,0,j),jqq(0,isig1,1,is1,0),
     #	   					jqq(0,isig3,2,is3,l))
 
           do k = 1,4
              mat(k,isig1,isig3,j,l,1) = 
     1    	   maa*clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3)
     2    	 + maz*clr(ifl(1,k),1,isig1)*clr(ifl(3,k),2,isig3)
     3    	 + mza*clr(ifl(1,k),2,isig1)*clr(ifl(3,k),1,isig3)
     4    	 + mzz*clr(ifl(1,k),2,isig1)*clr(ifl(3,k),2,isig3)
           enddo
	   endif

c
	   
          enddo
        enddo
      enddo
c 
c -----------------------
c
c charged current (k=5,6):
c      
      do l = 1,2	! gluon polarization (kart. basis)
      
      if (jlog2) then
      j =  2 ! g from upper line
            
      mww6 = contract_Tjj(wwww6(0,0,j),jqq(0,-1,1,is1,l),
     #      					jqq(0,-1,2,is3,0))
      mww5 = contract_Tjj(wwww5(0,0,j),jqq(0,-1,2,is3,0),
     #						jqq(0,-1,1,is1,l))
      mat(5,-1,-1,j,l,1) = mww5*clr(3,3,-1)**2!*sqrt(fpials(j))
      mat(6,-1,-1,j,l,1) = mww6*clr(3,3,-1)**2
      

      endif
  
      if (jlog3) then
      j =  3 ! g from lower line
            
      mww6 = contract_Tjj(wwww6(0,0,j),jqq(0,-1,1,is1,0),
     #						jqq(0,-1,2,is3,l))
      mww5 = contract_Tjj(wwww5(0,0,j),jqq(0,-1,2,is3,l),
     #						jqq(0,-1,1,is1,0))
      mat(5,-1,-1,j,l,1) = mww5*clr(3,3,-1)**2!*sqrt(fpials(j))
      mat(6,-1,-1,j,l,1) = mww6*clr(3,3,-1)**2
      
     
      endif
      enddo ! l-loop
     
c ----------------------------------------------------------------    
c
c check result for NC against madgraph output:

      if (ldebugm) then
        do k = 1,1
	  do lh = -1,1,2	! gluon pol. in hel. basis
 	    
	! need to distingiush between j=2,3 for gsign = -1, since 
 	! for incoming g: no g emission from other line
       	!    ... not accounted for here (take only sum of j=2 and j=3 pieces)
		
            m1hc(lh) = 0d0

c compute reference result from individual madgraph diagrams tampr(i,...):
 
  	    do  i = 475,630
		m1hc(lh) = m1hc(lh) + tampr(i,-1,-1,lh,k) 
	    enddo	
  	    do  i = 646,733
		m1hc(lh) = m1hc(lh) + tampr(i,-1,-1,lh,k) 
	    enddo	
  	    do  i = 749,836
		m1hc(lh) = m1hc(lh) + tampr(i,-1,-1,lh,k) 
	    enddo	
  
	  enddo ! end of lh loop 
	    
	    ! compute m1 in kartesian basis:
	    m1kc(1) =     (m1hc(-1)-m1hc(1))/sqrt(2d0)
	    m1kc(2) = -im*(m1hc(-1)+m1hc(1))/sqrt(2d0)
		    
		    
c print results from checks:

	    print*,"results for vertex-vertex contributions:"	
	    print*,"signs:",sign	
	    print*,"gign:",gsign	    
  	    print*," vbf for k=",k,"      "
	    
    	    print*,"ratio for l=1:",
     #	    	(mat(k,-1,-1,2,1,1)+mat(k,-1,-1,3,1,1))/
     #		(m1kc(1))
      	    print*,"ratio for l=2:",
     #	    	(mat(k,-1,-1,2,2,1)+mat(k,-1,-1,3,2,1))/
     #		(m1kc(2))
    	    print*,"res for l=1:",
     #	    	(mat(k,-1,-1,2,1,1)+mat(k,-1,-1,3,1,1)),
     #		(m1kc(1))
      	    print*,"res for l=2:",
     #	    	(mat(k,-1,-1,2,2,1)+mat(k,-1,-1,3,2,1)),
     #		(m1kc(2))
	  	
	   
        enddo ! end of k loop
      endif ! end of debugging block
      
c
c -----------------------------------------------------------------------
c ----------------------------------------------------------------    
c
c check result for CC against madgraph output:

      if (ldebugm) then
        do k = 5,6
	  do lh = -1,1,2	! gluon pol. in hel. basis
 	    
	! need to distingiush between j=2,3 for gsign = -1
	! for m1hj(lh,j):	
             m1hj(lh,2) = 0  ! g to upper line (j=2)
             m1hj(lh,3) = 0  ! g to lower line (j=3)	     
	     temp = m1hj(lh,2)

c compute reference result from individual madgraph diagrams tampr(i,...):
 
       	     do i = 303,374
	       m = 0
	  	do n = 1,36
	     	  if (i.eq.list1(n)) then
             	    temp = m1hj(lh,2) + tampr(i,-1,-1,lh,k)
	    	    m = m+1
	     	  endif    
	  	enddo
 	       m1hj(lh,2) = temp	
	       if(m.eq.0) then
	  	  m1hj(lh,3) = m1hj(lh,3) + tampr(i,-1,-1,lh,k) 
	       endif
       	     enddo	  
       	     do i = 385,400
	       m = 0
	  	do n = 1,8
	    	  if (i.eq.list2(n)) then
             	    temp = m1hj(lh,2) + tampr(i,-1,-1,lh,k)
	     	    m = m+1
	     	  endif    
	  	enddo
 	       m1hj(lh,2) = temp	
	       if(m.eq.0) then
	  	  m1hj(lh,3) = m1hj(lh,3) + tampr(i,-1,-1,lh,k) 
	       endif
       	     enddo	  
       	     do i = 421,444
	       m = 0
	  	do n = 1,12
	    	  if (i.eq.list3(n)) then
             	    temp = m1hj(lh,2) + tampr(i,-1,-1,lh,k)
	     	    m = m+1
	     	  endif    
	  	enddo
 	       m1hj(lh,2) = temp	
	       if(m.eq.0) then
	  	  m1hj(lh,3) = m1hj(lh,3) + tampr(i,-1,-1,lh,k) 
	       endif
       	     enddo	  
       
	     if (k.eq.5) then  ! for incoming g: no g emission from other line
       		if(.not.jlog2) m1hj(lh,2) = 0d0
       		if(.not.jlog3) m1hj(lh,3) = 0d0
	     else
       		if(.not.jlog2) m1hj(lh,3) = 0d0
       		if(.not.jlog3) m1hj(lh,2) = 0d0
             endif
  
	  enddo ! end of lh loop 
	    
	    ! compute m1 in kartesian basis:
	  do j = 2,3
	    m1kj(1,j) = (m1hj(-1,j)-m1hj(1,j))/sqrt(2d0)
	    m1kj(2,j) = -im*(m1hj(-1,j)+m1hj(1,j))/sqrt(2d0)
	  enddo    	  
		    
		    
c print results from checks:

	    print*,"signs:",sign	
	    print*,"gign:",gsign	    
  	    print*," vbf for k=",k,"      "
	    
    	    print*,"ratio for l=1:",
     #	    	(mat(k,-1,-1,2,1,1)+mat(k,-1,-1,3,1,1))/
     #		(m1kj(1,2)+m1kj(1,3))
      	    print*,"ratio for l=2:",
     #	    	(mat(k,-1,-1,2,2,1)+mat(k,-1,-1,3,2,1))/
     #		(m1kj(2,2)+m1kj(2,3))
    	    print*,"res for l=1:",
     #	    	(mat(k,-1,-1,2,1,1)+mat(k,-1,-1,3,1,1)),
     #		(m1kj(1,2)+m1kj(1,3))
      	    print*,"res for l=2:",
     #	    	(mat(k,-1,-1,2,2,1)+mat(k,-1,-1,3,2,1)),
     #		(m1kj(2,2)+m1kj(2,3))
c	  	
	   
        enddo ! end of k loop
      endif ! end of debugging block
      
c
c -----------------------------------------------------------------------
c -----------------------------------------------------------------------
c
c  prepare box diagrams: attach W+ and W- to external spinors
c 
c      isig = -1   : lefthanded spinors only coupling to W"s
c
c  Notation for virtual 2-component spinors and momenta
c
c  W+ attached to quark number i: psiwp(*,i) with momentum fqp(mu,i)
c  W- attached to quark number i: psiwm(*,i) with momentum fqm(mu,i)
c  
c  the fermion current corresponding to a quark line with the real emitted W+
c  attached next to quark number i is stored in jwp(mu,isig,is,i,l). Similarly
c  jwm(mu,isig,i,l) is the corresponding current for a W- attached next to 
c  quark i; for l = 0 jwp/jwm doesn"t contain a gluon; for l = 1,2 jwp/jwm 
c  includes sum of possible gluon couplings on quark leg i
c 
      l = 0    ! no gluon 
 
 	do i = 1,3,2  ! fermion ID
		
		if (i.eq.1) then
		   is = is1
		else
		   is = is3   
		endif   
 
	! only fermions with isig = -1 couple to W (-> no loop for isig)
	
	 call ket2c(psi(1,-1,i),.true.,p(0,i),-1,qp,wp,
     1  	    psiwp(1,i),fqp(0,i))
	 call ket2c(psi(1,-1,i),.true.,p(0,i),-1,qm,wm,
     1  	    psiwm(1,i),fqm(0,i))
	 call bra2c(psi(1,-1,i+1),.true.,p(0,i+1),-1,qp,wp,
     1  	    psiwp(1,i+1),fqp(0,i+1))
	 call bra2c(psi(1,-1,i+1),.true.,p(0,i+1),-1,qm,wm,
     1  	    psiwm(1,i+1),fqm(0,i+1))
     
	  
	 call curr6(-1,psi(1,-1,i+1),p(0,i+1),
     1  	       psiwp(1,i),fqp(0,i), jwp(0,-1,is,i,l)   )
	 call curr6(-1,psiwp(1,i+1),fqp(0,i+1),
     1  	       psi(1,-1,i),p(0,i), jwp(0,-1,is,i+1,l) )
	 call curr6(-1,psi(1,-1,i+1),p(0,i+1),
     1  	       psiwm(1,i),fqm(0,i), jwm(0,-1,is,i,l)   )
	 call curr6(-1,psiwm(1,i+1),fqm(0,i+1),
     1  	       psi(1,-1,i),p(0,i), jwm(0,-1,is,i+1,l) )
      

 	enddo ! i loop
 

c
c -----------------------------------------------------------------------
c
c keep structure of LO code, but replace jwp/jwm with 
c jwp/jwm(mu,isig,is,i,l=0:2) to take gluon radiation into account


c attach W+ to f1 or f2 / f3 or f4:

 	do i = 1,3,2  ! fermion ID (isigi=-1)  

		if (i.eq.1) then
		   is = is1
		else
		   is = is3   
		endif   

 
 ! gluon radiation from fermion i / i+1
 	 do l = 1,2
            	call ket2c(braketg(1,-1,i,l),.false.,pq(0,i),
     $                     -1,qp,wp,braketgWP(1,-1,i,l),pgwp(0,i))
            	call bra2c(braketg(1,-1,i+1,l),.false.,pq(0,i+1),
     $                 -1,qp,wp,braketgWP(1,-1,i+1,l),pgwp(0,i+1))
		
		call ket2r(psiwp(1,i),.false.,fqp(0,i),-1,
     $	    		    q,eps(0,l),braketWPg(1,-1,i,l),pwpg(0,i))      
            	call bra2r(psiwp(1,i+1),.false.,fqp(0,i+1),-1,
     $	    		q,eps(0,l),braketWPg(1,-1,i+1,l),pwpg(0,i+1))      


   	    call curr6(1,psi(1,-1,i+1),p(0,i+1),
     $	    			braketgWP(1,-1,i,l),pgwp(0,i),jwgi) 			
      	    call curr6(1,braketgWP(1,-1,i+1,l),pgwp(0,i+1),
     $	    			psi(1,-1,i),p(0,i),jwgii)		       
      	    call curr6(1,psi(1,-1,i+1),p(0,i+1),
     $			        braketWPg(1,-1,i,l),pwpg(0,i),jgwi)
      	    call curr6(1,braketWPg(1,-1,i+1,l),pwpg(0,i+1),
     $			        psi(1,-1,i),p(0,i),jgwii)
 
  ! gluon radiation from fermion i+1 / i
     	    call curr6(1,braketg(1,-1,i+1,l),pq(0,i+1),
     $	    			 psiwp(1,i),fqp(0,i),jwg0i)
     	    call curr6(1,psiwp(1,i+1),fqp(0,i+1),
     $	    			braketg(1,-1,i,l),pq(0,i),jwg0ii)


	    do mu = 0,5
	    	   
 		   jwp(mu,-1,is,i,l) = jwgi(mu,-1)+
     $			jgwi(mu,-1)+jwg0i(mu,-1)  ! W+ & g emission from i/i+1 line 
  		   
		   jwp(mu,-1,is,i+1,l) = jwgii(mu,-1)+
     $			jgwii(mu,-1)+jwg0ii(mu,-1)  ! W+ & g emission from i/i+1 line 
    	    	   
	    enddo  !mu       
       
       enddo ! l = 1,2
c
       enddo ! i loop
c


c
c ----------------------------------
c

c attach W- to f1 or f2 / f3 or f4:

 	do i = 1,3,2  ! fermion ID (isigi=-1)  
		
		if (i.eq.1) then
		   is = is1
		else
		   is = is3   
		endif   

 
 ! gluon radiation from fermion i / i+1
 	 do l = 1,2
            	call ket2c(braketg(1,-1,i,l),.false.,pq(0,i),
     $                 -1,qm,wm,braketgWM(1,-1,i,l),pgwm(0,i))
            	call bra2c(braketg(1,-1,i+1,l),.false.,pq(0,i+1),
     $                 -1,qm,wm,braketgWM(1,-1,i+1,l),pgwm(0,i+1))

            	call ket2r(psiwm(1,i),.false.,fqm(0,i),-1,
     $	    		q,eps(0,l),braketWMg(1,-1,i,l),pwmg(0,i))      
            	call bra2r(psiwm(1,i+1),.false.,fqm(0,i+1),-1,
     $	    		q,eps(0,l),braketWMg(1,-1,i+1,l),pwmg(0,i+1))      

   	    call curr6(1,psi(1,-1,i+1),p(0,i+1),
     $	    			braketgWM(1,-1,i,l),pgwm(0,i),jwgi) 			
      	    call curr6(1,braketgWM(1,-1,i+1,l),pgwm(0,i+1),
     $	    			psi(1,-1,i),p(0,i),jwgii)		       
      	    call curr6(1,psi(1,-1,i+1),p(0,i+1),
     $			        braketWMg(1,-1,i,l),pwmg(0,i),jgwi)
      	    call curr6(1,braketWMg(1,-1,i+1,l),pwmg(0,i+1),
     $			        psi(1,-1,i),p(0,i),jgwii)
 
  ! gluon radiation from fermion i+1 / i
     	    call curr6(1,braketg(1,-1,i+1,l),pq(0,i+1),
     $	    			psiwm(1,i),fqm(0,i),jwg0i)
     	    call curr6(1,psiwm(1,i+1),fqm(0,i+1),
     $	    			braketg(1,-1,i,l),pq(0,i),jwg0ii)


	    do mu = 0,5
	    	   
 		   jwm(mu,-1,is,i,l) = jwgi(mu,-1)+
     $			jgwi(mu,-1)+jwg0i(mu,-1)  ! W- & g emission from i/i+1 line 
  		   
		   jwm(mu,-1,is,i+1,l) = jwgii(mu,-1)+
     $			jgwii(mu,-1)+jwg0ii(mu,-1)  ! W- & g emission from i/i+1 line 
    	    	   
	    enddo  !mu       
       
       enddo ! l = 1,2
c
       enddo ! i loop
c

c 
c -----------------------------------------------------------------------
c -----------------------------------------------------------------------
c
c now calculate the Vertex-box diagrams; get t-channel W currents first

     	
      do ka = 1,2
         kb = 3-ka
	 
	 if (kb.eq.1) then
	 	is = is1
	 else
	 	is = is3
	 endif	
	 
c ka=1 and kb=2 is for "box correction" to upper line
c ka=2 and kb=1 is for "box correction" to lower line

      do j = 2,3 ! g from upper/lower line

         if (ka.eq.1) then
            zp(4,j) = -dcmplx(qmp(0,j),qmp(3,j))
            zp(5,j) = -dcmplx(qmp(1,j),qmp(2,j))
            zm(4,j) = -dcmplx(qpm(0,j),qpm(3,j))
            zm(5,j) = -dcmplx(qpm(1,j),qpm(2,j))
         else
            zm(4,j) = dcmplx(qmp(0,j),qmp(3,j))
            zm(5,j) = dcmplx(qmp(1,j),qmp(2,j))
            zp(4,j) = dcmplx(qpm(0,j),qpm(3,j))
            zp(5,j) = dcmplx(qpm(1,j),qpm(2,j))
         endif
	 	 
      enddo ! j
	 
         do isig = -1,1,2

c--------
	 
	 ! l = 0:
	   
	    call contract_T2j(NCwpa(0,0,ka,1+ka),
     #			    jqq(0,isig,kb,is,0),epsa0) !for W+ current
            call contract_T2j(NCwpz(0,0,ka,1+ka),
     #			    jqq(0,isig,kb,is,0),epsz0)
	   
	   do l = 1,2
            
	    call contract_T2j(NCwpa(0,0,ka,4-ka),
     #			    jqq(0,isig,kb,is,l),epsa) !for W+ current
            call contract_T2j(NCwpz(0,0,ka,4-ka),
     #			    jqq(0,isig,kb,is,l),epsz)
            
	    do mu = 0,3
               epsNCwp(mu,isig,3,ka,1+ka,0) =
     1              epsa0(mu)*clr(3,1,isig)+epsz0(mu)*clr(3,2,isig)
               epsNCwp(mu,isig,4,ka,1+ka,0) =
     1              epsa0(mu)*clr(4,1,isig)+epsz0(mu)*clr(4,2,isig)
 
               epsNCwp(mu,isig,3,ka,4-ka,l) =
     1              epsa(mu)*clr(3,1,isig)+epsz(mu)*clr(3,2,isig)
               epsNCwp(mu,isig,4,ka,4-ka,l) =
     1              epsa(mu)*clr(4,1,isig)+epsz(mu)*clr(4,2,isig)
            enddo ! mu
	    
            enddo !l

c--------
	    
	 ! l = 0:
	   
	    call contract_T2j(NCwma(0,0,ka,1+ka),
     #		    		jqq(0,isig,kb,is,0),epsa0) !for W- current
            call contract_T2j(NCwmz(0,0,ka,1+ka),
     #		    		jqq(0,isig,kb,is,0),epsz0)
	   
	   do l = 1,2
            
	    call contract_T2j(NCwma(0,0,ka,4-ka),
     #			    jqq(0,isig,kb,is,l),epsa) !for W+ current
            call contract_T2j(NCwmz(0,0,ka,4-ka),
     #			    jqq(0,isig,kb,is,l),epsz)
            

	    do mu = 0,3
               epsNCwm(mu,isig,3,ka,1+ka,0) =
     1              epsa0(mu)*clr(3,1,isig)+epsz0(mu)*clr(3,2,isig)
               epsNCwm(mu,isig,4,ka,1+ka,0) =
     1              epsa0(mu)*clr(4,1,isig)+epsz0(mu)*clr(4,2,isig)
 
               epsNCwm(mu,isig,3,ka,4-ka,l) =
     1              epsa(mu)*clr(3,1,isig)+epsz(mu)*clr(3,2,isig)
               epsNCwm(mu,isig,4,ka,4-ka,l) =
     1              epsa(mu)*clr(4,1,isig)+epsz(mu)*clr(4,2,isig)
            enddo ! mu
	    
           enddo !l
	    
	    
 	  do l = 0,2
	  do j = 2,3
		  
            do mu = 4,5           ! add momentum info to the currents
               epsNCwp(mu,isig,3,ka,j,l) = zp(mu,j)
               epsNCwp(mu,isig,4,ka,j,l) = zp(mu,j)
               epsNCwm(mu,isig,3,ka,j,l) = zm(mu,j)
               epsNCwm(mu,isig,4,ka,j,l) = zm(mu,j)
            enddo
	 
	   enddo !j
           enddo ! l
         enddo !isig


c -----------------------	 
	 
c and same for the CC processes (W attached to j43 or j21 current)

         isig = -1

	! l = 0:
	   
         call contract_T1j(CCwpa(0,0,ka,1+ka),
     #			 jqq(0,isig,kb,is,0),epsa0) !for W+ exchange
         call contract_T1j(CCwpz(0,0,ka,1+ka),
     #			 jqq(0,isig,kb,is,0),epsz0)
     
   	  do l = 1,2
	   
         call contract_T1j(CCwpa(0,0,ka,4-ka),
     #			 jqq(0,isig,kb,is,l),epsa) 
         call contract_T1j(CCwpz(0,0,ka,4-ka),
     #			 jqq(0,isig,kb,is,l),epsz)
   
 
         do mu = 0,3 
            epsCCwp(mu,isig,3,ka,1+ka,0) =
     1           epsa0(mu)*clr(3,1,isig)+epsz0(mu)*clr(3,2,isig)
            epsCCwp(mu,isig,4,ka,1+ka,0) =
     1           epsa0(mu)*clr(4,1,isig)+epsz0(mu)*clr(4,2,isig)
     
            epsCCwp(mu,isig,3,ka,4-ka,l) =
     1           epsa(mu)*clr(3,1,isig)+epsz(mu)*clr(3,2,isig)
            epsCCwp(mu,isig,4,ka,4-ka,l) =
     1           epsa(mu)*clr(4,1,isig)+epsz(mu)*clr(4,2,isig)
            enddo ! mu

	 if (ka.eq.1) then
            qepsCCwp(ka,1+ka,0) = -dotrc(qmp(0,1+ka),epsz0)*zm2i(2)
            qepsCCwp(ka,4-ka,l) = -dotrc(qmp(0,4-ka),epsz)*zm2i(2)
        else
            qepsCCwp(ka,1+ka,0) = dotrc(qpm(0,1+ka),epsz0)*zm2i(2)
            qepsCCwp(ka,4-ka,l) = dotrc(qpm(0,4-ka),epsz)*zm2i(2)
        endif !ka
	    
            enddo !l
         
	 	 
c ------------
     
	! l = 0:
	   
         call contract_T1j(CCwma(0,0,ka,1+ka),
     #			 jqq(0,isig,kb,is,0),epsa0) !for W- exchange
         call contract_T1j(CCwmz(0,0,ka,1+ka),
     #			 jqq(0,isig,kb,is,0),epsz0)
     
   	  do l = 1,2
	   
         call contract_T1j(CCwma(0,0,ka,4-ka),
     #			 jqq(0,isig,kb,is,l),epsa) 
         call contract_T1j(CCwmz(0,0,ka,4-ka),
     #			 jqq(0,isig,kb,is,l),epsz)
   
         do mu = 0,3 
            epsCCwm(mu,isig,3,ka,1+ka,0) =
     1           epsa0(mu)*clr(3,1,isig)+epsz0(mu)*clr(3,2,isig)
            epsCCwm(mu,isig,4,ka,1+ka,0) =
     1           epsa0(mu)*clr(4,1,isig)+epsz0(mu)*clr(4,2,isig)
     
            epsCCwm(mu,isig,3,ka,4-ka,l) =
     1           epsa(mu)*clr(3,1,isig)+epsz(mu)*clr(3,2,isig)
            epsCCwm(mu,isig,4,ka,4-ka,l) =
     1           epsa(mu)*clr(4,1,isig)+epsz(mu)*clr(4,2,isig)
            enddo ! mu

	
	if (ka.eq.1) then
            qepsCCwm(ka,1+ka,0) = -dotrc(qpm(0,1+ka),epsz0)*zm2i(2)
            qepsCCwm(ka,4-ka,l) = -dotrc(qpm(0,4-ka),epsz)*zm2i(2)
        else
            qepsCCwm(ka,1+ka,0) = dotrc(qmp(0,1+ka),epsz0)*zm2i(2)
            qepsCCwm(ka,4-ka,l) = dotrc(qmp(0,4-ka),epsz)*zm2i(2)
        endif !ka
	
	    
            enddo !l
	 
	 	 
 	do l = 0,2
	do j = 2,3
        
	 do mu = 4,5            ! add momentum info to the currents
            epsCCwp(mu,isig,3,ka,j,l) = zp(mu,j)
            epsCCwp(mu,isig,4,ka,j,l) = zp(mu,j)
            epsCCwm(mu,isig,3,ka,j,l) = zm(mu,j)
            epsCCwm(mu,isig,4,ka,j,l) = zm(mu,j)
         enddo 
 
 	enddo !j   	 
        enddo ! l		 
	 
	 
      enddo ! ka loop
       
    
      do l = 0,2
      jj21m(l) = dotcc(wm,jqq(0,-1,1,is1,l))
      jj43p(l) = dotcc(wp,jqq(0,-1,2,is3,l))
     
      jj21p(l) = dotcc(wp,jqq(0,-1,1,is1,l))
      jj43m(l) = dotcc(wm,jqq(0,-1,2,is3,l))
      enddo ! l		 
      	
	
c --------------------------------------------------------------------
c
c now construct the contribution to the amplitude by current contraction
c  
c neutral current first:  
      
      do isig = -1,1,2          ! 2 bosons attached to 12 line  
	
	 do kf=3,4 !for wp or wm
            qepswm(kf,2,0) = -dotrc(qpm(0,2),
     #		    epsNCwm(0,isig,kf,1,2,0))*zm2i(3)
            qepswp(kf,2,0) = -dotrc(qmp(0,2),
     #	    	    epsNCwp(0,isig,kf,1,2,0))*zm2i(3)

            
	    do l = 1,2
	    qepswm(kf,3,l) = -dotrc(qpm(0,3),
     #		    epsNCwm(0,isig,kf,1,3,l))*zm2i(3)
            qepswp(kf,3,l) = -dotrc(qmp(0,3),
     #	    	    epsNCwp(0,isig,kf,1,3,l))*zm2i(3)
            enddo !l
	    
         enddo !kf
	 
	 
        do l = 1,2

         do k = 1,2                !uucc, uuss, g from upper line
           m1u = dotcc(jwp(0,-1,is1,1,l),epsNCwm(0,isig,ifl(3,k),1,2,0))      !uW+(tW-)u   
     1         - qepswm(ifl(3,k),2,0)*jj21p(l)
           m2u = dotcc(jwm(0,-1,is1,2,l),epsNCwp(0,isig,ifl(3,k),1,2,0))      !u(tW+)W-u
     1         + qepswp(ifl(3,k),2,0)*jj21m(l)
            
         	                   !uucc, uuss, g from lower line
	   m1l = dotcc(jwp(0,-1,is1,1,0),epsNCwm(0,isig,ifl(3,k),1,3,l))      !uW+(tW-)u   
     1         - qepswm(ifl(3,k),3,l)*jj21p(0)
           m2l = dotcc(jwm(0,-1,is1,2,0),epsNCwp(0,isig,ifl(3,k),1,3,l))      !u(tW+)W-u
     1         + qepswp(ifl(3,k),3,l)*jj21m(0)
            
	    
	    mat(k,-1,isig,2,l,2) = clr(3,3,-1)**2 * (m1u+m2u)
	    mat(k,-1,isig,3,l,2) = clr(3,3,-1)**2 * (m1l+m2l)	      
         enddo !k

         do k = 3,4             !ddcc, ddss, g from upper line
           m1u = dotcc(jwp(0,-1,is1,2,l),epsNCwm(0,isig,ifl(3,k),1,2,0))     !d(tW-)W+d
     1         + qepswm(ifl(3,k),2,0)*jj21p(l)
           m2u = dotcc(jwm(0,-1,is1,1,l),epsNCwp(0,isig,ifl(3,k),1,2,0))     !dW-(tW+)d
     1         - qepswp(ifl(3,k),2,0)*jj21m(l)
            
        	                   !ddcc, ddss, g from lower line
	   m1l = dotcc(jwp(0,-1,is1,2,0),epsNCwm(0,isig,ifl(3,k),1,3,l))     !d(tW-)W+d
     1         + qepswm(ifl(3,k),3,l)*jj21p(0)
           m2l = dotcc(jwm(0,-1,is1,1,0),epsNCwp(0,isig,ifl(3,k),1,3,l))     !dW-(tW+)d
     1         - qepswp(ifl(3,k),3,l)*jj21m(0)


           mat(k,-1,isig,2,l,2) = clr(3,3,-1)**2 * (m1u+m2u)
           mat(k,-1,isig,3,l,2) = clr(3,3,-1)**2 * (m1l+m2l)
         enddo !k

      enddo	! l
      enddo ! isig
 
c ---------------
c
c charged current with 2 bosons attached to 12 line:
c
 
      do l = 1,2
  
      isig = -1
      k = 5                     !udsc, g from upper line
      m1u = dotcc(jwp(0,-1,is1,2,l),epsCCwm(0,isig,ifl(1,k),1,2,0))            !u(tAZ)W+d
      m2u = dotcc(jwp(0,-1,is1,1,l),epsCCwm(0,isig,ifl(2,k),1,2,0))            !uW+(tAZ)d
      z1u = qepsCCwm(1,2,0)*jj21p(l)*
     # 		(clr(ifl(1,k),2,isig)-clr(ifl(2,k),2,isig))
  
                               !udsc, g from lower line
      m1l = dotcc(jwp(0,-1,is1,2,0),epsCCwm(0,isig,ifl(1,k),1,3,l))            !u(tAZ)W+d
      m2l = dotcc(jwp(0,-1,is1,1,0),epsCCwm(0,isig,ifl(2,k),1,3,l))            !uW+(tAZ)d
      z1l = qepsCCwm(1,3,l)*jj21p(0)*
     # 		(clr(ifl(1,k),2,isig)-clr(ifl(2,k),2,isig))
     
      mat(k,-1,isig,2,l,2) = clr(3,3,-1)**2 * (m1u+m2u+z1u)
      mat(k,-1,isig,3,l,2) = clr(3,3,-1)**2 * (m1l+m2l+z1l)


c ------------
      
      k = 6                     !ducs, g from upper line
      m1u = dotcc(jwm(0,-1,is1,2,l),epsCCwp(0,isig,ifl(1,k),1,2,0))            !d(tAZ)W-u   
      m2u = dotcc(jwm(0,-1,is1,1,l),epsCCwp(0,isig,ifl(2,k),1,2,0))            !dW-(tAZ)u
      z1u =qepsCCwp(1,2,0)*jj21m(l)*
     # 		(clr(ifl(1,k),2,isig)-clr(ifl(2,k),2,isig))
    
     	                        !ducs, g from lower line
      m1l = dotcc(jwm(0,-1,is1,2,0),epsCCwp(0,isig,ifl(1,k),1,3,l))            !d(tAZ)W-u   
      m2l = dotcc(jwm(0,-1,is1,1,0),epsCCwp(0,isig,ifl(2,k),1,3,l))            !dW-(tAZ)u
      z1l =qepsCCwp(1,3,l)*jj21m(0)*
     # 		(clr(ifl(1,k),2,isig)-clr(ifl(2,k),2,isig))

      mat(k,-1,isig,2,l,2) = clr(3,3,-1)**2 * (m1u+m2u+z1u)
      mat(k,-1,isig,3,l,2) = clr(3,3,-1)**2 * (m1l+m2l+z1l)
    
      
      enddo	! l 	 
 
c

c -----------------------------------------------------------------
c

c repeat the same for 2 bosons attached to 34 line
c
c 	Neutral current:

      do isig = -1,1,2

        do kf=3,4
 	    do l = 1,2
            qepswm(kf,2,l) = dotrc(qmp(0,2),
     #		        epsNCwm(0,isig,kf,2,2,l))*zm2i(3)
            qepswp(kf,2,l) = dotrc(qpm(0,2),
     #	    		epsNCwp(0,isig,kf,2,2,l))*zm2i(3)
             enddo !l	    
       
	    qepswm(kf,3,0) = dotrc(qmp(0,3),
     #		        epsNCwm(0,isig,kf,2,3,0))*zm2i(3)
            qepswp(kf,3,0) = dotrc(qpm(0,3),
     #	    		epsNCwp(0,isig,kf,2,3,0))*zm2i(3)
	enddo !kf
       
	
      do l = 1,2

         do k = 1,3,2                  !uucc, ddcc, g from upper line
           m1u = dotcc(jwp(0,-1,is3,3,0),epsNCwm(0,isig,ifl(1,k),2,2,l)) !cW+(tW-)c   
     1         - qepswm(ifl(1,k),2,l)*jj43p(0)
           m2u = dotcc(jwm(0,-1,is3,4,0),epsNCwp(0,isig,ifl(1,k),2,2,l)) !c(tW+)W-c
     1         + qepswp(ifl(1,k),2,l)*jj43m(0)
             
	     
         	                   	!uucc, ddcc, g from lower line
	   m1l = dotcc(jwp(0,-1,is3,3,l),epsNCwm(0,isig,ifl(1,k),2,3,0)) !cW+(tW-)c   
     1         - qepswm(ifl(1,k),3,0)*jj43p(l)
           m2l = dotcc(jwm(0,-1,is3,4,l),epsNCwp(0,isig,ifl(1,k),2,3,0)) !c(tW+)W-c
     1         + qepswp(ifl(1,k),3,0)*jj43m(l)
            
	    mat(k,isig,-1,2,l,3) = clr(3,3,-1)**2 * (m1u+m2u)
            mat(k,isig,-1,3,l,3) = clr(3,3,-1)**2 * (m1l+m2l)
         enddo !k
	
	 
         do k = 2,4,2                  !uuss, ddss, g from upper line
           m1u = dotcc(jwp(0,-1,is3,4,0),epsNCwm(0,isig,ifl(1,k),2,2,l))      !s(tW-)W+s
     1         + qepswm(ifl(1,k),2,l)*jj43p(0)
           m2u = dotcc(jwm(0,-1,is3,3,0),epsNCwp(0,isig,ifl(1,k),2,2,l))      !sW-(tW+)s
     1         - qepswp(ifl(1,k),2,l)*jj43m(0)

          	                   	!uuss, ddss, g from lower line
          m1l = dotcc(jwp(0,-1,is3,4,l),epsNCwm(0,isig,ifl(1,k),2,3,0))      !s(tW-)W+s
     1         + qepswm(ifl(1,k),3,0)*jj43p(l)
          m2l = dotcc(jwm(0,-1,is3,3,l),epsNCwp(0,isig,ifl(1,k),2,3,0))      !sW-(tW+)s
     1         - qepswp(ifl(1,k),3,0)*jj43m(l)
     
            mat(k,isig,-1,2,l,3) = clr(3,3,-1)**2 * (m1u+m2u)
            mat(k,isig,-1,3,l,3) = clr(3,3,-1)**2 * (m1l+m2l)
         enddo !k
	 
      
       enddo ! l
       
      enddo !isig
      
c --------------------------      
    
c charged current, 2 bosons from lower line
      
      do l = 1,2
    
      isig = -1
      k = 5                     !udsc, g from upper line
     
      m1u = dotcc(jwm(0,-1,is3,4,0),epsCCwp(0,isig,ifl(3,k),2,2,l))            !s(tAZ)W-c
      m2u = dotcc(jwm(0,-1,is3,3,0),epsCCwp(0,isig,ifl(4,k),2,2,l))            !sW-(tAZ)c
      z1u = qepsCCwp(2,2,l)*jj43m(0)*
     #  	(clr(ifl(3,k),2,isig)-clr(ifl(4,k),2,isig))
     
      		                !udsc, g from lower line
      m1l = dotcc(jwm(0,-1,is3,4,l),epsCCwp(0,isig,ifl(3,k),2,3,0))            !s(tAZ)W-c
      m2l = dotcc(jwm(0,-1,is3,3,l),epsCCwp(0,isig,ifl(4,k),2,3,0))            !sW-(tAZ)c
      z1l = qepsCCwp(2,3,0)*jj43m(l)*
     #  	(clr(ifl(3,k),2,isig)-clr(ifl(4,k),2,isig))
     
      mat(k,isig,-1,2,l,3) = clr(3,3,-1)**2 * (m1u+m2u+z1u)
      mat(k,isig,-1,3,l,3) = clr(3,3,-1)**2 * (m1l+m2l+z1l)
      
     
c ------------------------      

     
      k = 6                     !ducs, g from upper line
     
      m1u = dotcc(jwp(0,-1,is3,4,0),epsCCwm(0,isig,ifl(3,k),2,2,l))            !c(tAZ)W+s   
      m2u = dotcc(jwp(0,-1,is3,3,0),epsCCwm(0,isig,ifl(4,k),2,2,l))            !cW+(tAZ)s
      z1u = qepsCCwm(2,2,l)*jj43p(0)*
     # 		(clr(ifl(3,k),2,isig)-clr(ifl(4,k),2,isig))
      
                           	!ducs, g from lower line
      m1l = dotcc(jwp(0,-1,is3,4,l),epsCCwm(0,isig,ifl(3,k),2,3,0))            !c(tAZ)W+s   
      m2l = dotcc(jwp(0,-1,is3,3,l),epsCCwm(0,isig,ifl(4,k),2,3,0))            !cW+(tAZ)s
      z1l = qepsCCwm(2,3,0)*jj43p(l)*
     # 		(clr(ifl(3,k),2,isig)-clr(ifl(4,k),2,isig))
     
      mat(k,isig,-1,2,l,3) = clr(3,3,-1)**2 * (m1u+m2u+z1u)
      mat(k,isig,-1,3,l,3) = clr(3,3,-1)**2 * (m1l+m2l+z1l)


      enddo	! l 	 

c ----------------------------------------------------------------    
c
c take special care of processes with incoming gluons:

	if (jmin.eq.3) then
	  do k = 1,6
	    do isig1 = -1,1,2
	      do isig3 = -1,1,2
	      do l = 1,2
	        do i = 2,3
		  mat(k,isig1,isig3,2,l,i) = 0d0
		enddo
	      enddo
	      enddo
	    enddo
	  enddo	  
	endif	
	
	if (jmax.eq.2) then
	  do k = 1,6
	    do isig1 = -1,1,2
	      do isig3 = -1,1,2
	      do l = 1,2
	        do i = 2,3
		  mat(k,isig1,isig3,3,l,i) = 0d0
		enddo
	      enddo
	      enddo
	    enddo
	  enddo	  
	endif	


c ----------------------------------------------------------------    
c
c check result for NC as obtained above against madgraph output 
c		(sum up relevant mg amplitudes in the following):

	! need to distingiush between j=2,3 to account for crossed processes
	! in the way computed here only uncrossed processes with outgoing 
	!	gluons can be cecked!

      if (ldebugm) then

       do k = 1,2
c       do k = 2,2
	
	  do lh = -1,1,2	! gluon pol. in hel. basis
 	    			
c compute reference result from individual madgraph diagrams tampr(i,...):
	
	    do j = 2,3
	    	m1haj(lh,j) = 0d0  ! graphs with A exchange
	    	m1hzj(lh,j) = 0d0  !		 Z exchange
           enddo 
	   m1hc(lh) = 0d0 ! don"t distingiush between j=2,3

					
	  if (k.eq.1) then  
	    do  i = 49,73
		m1hc(lh) = m1hc(lh) + tampr(i,-1,-1,lh,k) 
	    enddo	
	    do  i = 184,208
		m1hc(lh) = m1hc(lh) + tampr(i,-1,-1,lh,k) 
	    enddo
	    do  i = 455,474
		m1hc(lh) = m1hc(lh) + tampr(i,-1,-1,lh,k) 
	    enddo
	    do  i = 631,645
		m1hc(lh) = m1hc(lh) + tampr(i,-1,-1,lh,k) 
	    enddo	
	    do  i = 734,748
		m1hc(lh) = m1hc(lh) + tampr(i,-1,-1,lh,k) 
	    enddo
	    
	    do m = 1,20
		i = lncba2(m)
		m1haj(lh,2) = m1haj(lh,2) + tampr(i,-1,-1,lh,k) 
	    enddo
	    do m = 1,20
		i = lncba3(m)
		m1haj(lh,3) = m1haj(lh,3) + tampr(i,-1,-1,lh,k) 
	    enddo
	  
	  elseif (k.eq.2) then  
	    do  i = 49,73
		m1hc(lh) = m1hc(lh) + tampr(i,-1,-1,lh,k) 
	    enddo	
	    do  i = 368,392
		m1hc(lh) = m1hc(lh) + tampr(i,-1,-1,lh,k) 
	    enddo
	    do  i = 455,474
		m1hc(lh) = m1hc(lh) + tampr(i,-1,-1,lh,k) 
	    enddo
	    do  i = 631,645
		m1hc(lh) = m1hc(lh) + tampr(i,-1,-1,lh,k) 
	    enddo	
	    do  i = 734,748
		m1hc(lh) = m1hc(lh) + tampr(i,-1,-1,lh,k) 
	    enddo
	  
	    do m = 1,30
		i = lncbz2k(m)
		m1hzj(lh,2) = m1hzj(lh,2) + tampr(i,-1,-1,lh,k) 
	    enddo
	    do m = 1,30
		i = lncbz3k(m)
		m1hzj(lh,3) = m1hzj(lh,3) + tampr(i,-1,-1,lh,k) 
	    enddo
	    do m = 1,20
		i = lncba2k(m)
		m1haj(lh,2) = m1haj(lh,2) + tampr(i,-1,-1,lh,k) 
	    enddo
	    do m = 1,20
		i = lncba3k(m)
		m1haj(lh,3) = m1haj(lh,3) + tampr(i,-1,-1,lh,k) 
	    enddo
	  
	  endif  !k
	
	  
	     if (k.ne.3) then  ! for incoming g: no g emission from other line
       		if(.not.jlog2) m1haj(lh,2) = 0d0
       		if(.not.jlog3) m1haj(lh,3) = 0d0
       		if(.not.jlog2) m1hzj(lh,2) = 0d0
       		if(.not.jlog3) m1hzj(lh,3) = 0d0
 	     else
       		if(.not.jlog2) m1haj(lh,3) = 0d0
       		if(.not.jlog3) m1haj(lh,2) = 0d0
       		if(.not.jlog2) m1hzj(lh,3) = 0d0
       		if(.not.jlog3) m1hzj(lh,2) = 0d0
	     endif
	
	  enddo ! end of lh loop 	  
	    
		    
    ! compute m1 in kartesian basis:
	    
	    do l = 1,2
	      do j = 2,3
	    	m1kaj(l,j) = 0d0
	    	m1kzj(l,j) = 0d0
	      enddo
	      m1kc(l) = 0d0
	    enddo  	    
	    
	  do j = 2,3
	    m1kaj(1,j) =     (m1haj(-1,j)-m1haj(1,j))/sqrt(2d0)
	    m1kaj(2,j) = -im*(m1haj(-1,j)+m1haj(1,j))/sqrt(2d0)
	    m1kzj(1,j) =     (m1hzj(-1,j)-m1hzj(1,j))/sqrt(2d0)
	    m1kzj(2,j) = -im*(m1hzj(-1,j)+m1hzj(1,j))/sqrt(2d0)
	  enddo    	
	  m1kc(1) =	(m1hc(-1)-m1hc(1))/sqrt(2d0)
	  m1kc(2) = -im*(m1hc(-1)+m1hc(1))/sqrt(2d0)

c
c ---------------------------------		 
c
       if (ldebugm) then

        do l = 1,2 
	isig = -1
	  	  
	  matot = m1kc(l)
	  ratot = mat(k,-1,isig,2,l,3) + mat(k,-1,isig,2,l,2) +
     # 		  mat(k,-1,isig,3,l,3) + mat(k,-1,isig,3,l,2)
     
	  
          print*,"TEST A & Z contributions for vertex-box:"
	    print*,"signs:",sign	
	    print*,"gign:",gsign	    
          print*," k,l = ",k,l,", isig= ",isig
	  print*,"jmin,jmax:",jmin,jmax
          print*,"follows: coded mat,mg"
	  print*
	  print*,"rat: ",ratot/matot
    	  print*,"res: ",ratot,matot
	  print*


 	enddo !l
      endif ! "local" debugging block box/vertex (NC)

c ---------------------------------
    
      enddo ! k
          
      endif ! debugging


c --------------------------------------------------------------------------
c ----------------------------------------------------------------    
c
c check result for CC as obtained above against madgraph output 
c		(sum up relevant mg amplitudes in the following):


      if (ldebugm) then

       do k = 5,6
	
	  do lh = -1,1,2	! gluon pol. in hel. basis
 	    
	! need to distingiush between j=2,3
		
c compute reference result from individual madgraph diagrams tampr(i,...):
	
	    do j = 2,3
	    	m1haj(lh,j) = 0d0  ! graphs with A exchange
	    	m1hzj(lh,j) = 0d0  !		 Z exchange
           enddo


	    do m = 1,20
		i = listba2(m)
		m1haj(lh,2) = m1haj(lh,2) + tampr(i,-1,-1,lh,k) 
	    enddo
		
	    do m = 1,20
		i = listba3(m)
		m1haj(lh,3) = m1haj(lh,3) + tampr(i,-1,-1,lh,k) 
	    enddo
	
	    do m = 1,30
		i = listbz2(m)
		m1hzj(lh,2) = m1hzj(lh,2) + tampr(i,-1,-1,lh,k) 
	    enddo
	
	    do m = 1,30
		i = listbz3(m)
		m1hzj(lh,3) = m1hzj(lh,3) + tampr(i,-1,-1,lh,k) 
	    enddo
	     
			 
	     if (k.eq.5) then  ! for incoming g: no g emission from other line
       		if(.not.jlog2) m1haj(lh,2) = 0d0
       		if(.not.jlog3) m1haj(lh,3) = 0d0
       		if(.not.jlog2) m1hzj(lh,2) = 0d0
       		if(.not.jlog3) m1hzj(lh,3) = 0d0
	     else
       		if(.not.jlog2) m1haj(lh,3) = 0d0
       		if(.not.jlog3) m1haj(lh,2) = 0d0
       		if(.not.jlog2) m1hzj(lh,3) = 0d0
       		if(.not.jlog3) m1hzj(lh,2) = 0d0
             endif
	  
	  
	  enddo ! end of lh loop


    ! compute m1 in kartesian basis:
	    
	    do l = 1,2
	      do j = 2,3
	    	m1kaj(l,j) = 0d0
	    	m1kzj(l,j) = 0d0
	      enddo
	    enddo  	    
	    
	  do j = 2,3
	    m1kaj(1,j) =     (m1haj(-1,j)-m1haj(1,j))/sqrt(2d0)
	    m1kaj(2,j) = -im*(m1haj(-1,j)+m1haj(1,j))/sqrt(2d0)
	    m1kzj(1,j) =     (m1hzj(-1,j)-m1hzj(1,j))/sqrt(2d0)
	    m1kzj(2,j) = -im*(m1hzj(-1,j)+m1hzj(1,j))/sqrt(2d0)
	  enddo    	

c
c
c ---------------------------------		 

      if (ldebugm) then

        do l = 1,2 
	isig = -1

	  
	  matot = m1kzj(l,2)+m1kzj(l,3)+m1kaj(l,2)+m1kaj(l,3)
	  ratot = mat(k,-1,isig,2,l,3) + mat(k,-1,isig,2,l,2) +
     # 		  mat(k,-1,isig,3,l,3) + mat(k,-1,isig,3,l,2)
     
	  
          print*,"TEST A & Z contributions:"
	    print*,"signs:",sign	
	    print*,"gign:",gsign	    
          print*," k,l = ",k,l,", isig= ",isig
          print*,"follows: coded mat,mg"
	  print*,"rat: ",ratot/matot
    	  print*,"res: ",ratot,matot
	  print*


 	enddo !l
      endif ! "local" debugging block for box/vertex (CC)

c ---------------------------------
    
      enddo ! k
          
      endif ! debugging


c --------------------------------------------------------------------------
c -----------------------------------------------------------------------
c --------------------------------------------------------------------------
c
c next come the A/Z-->WW currents attached to the quark lines. 
c The most effective structure is the contraction of two 
c polarization vectors with one fermion line. First build these effective 
c polarization vectors from the currents aww(mu) and azz(mu)
c
c NOTE: the aww and azz currents are conserved (checked numerically). 
C	Hence there is no need to consider q^mu * q^nu/m_Z^2  terms 
C	in the Z boson propagator

c box correction to upper line: polarization vectors are 
c    jqq(mu,isig3,2,is3,l) = j43  with momentum    p43        and
c    aww/zww(mu)         	  with momentum    pww


c attach V=A/Z to f1:
	 
       
      i=1
      
         do isig = -1,1,2 ! fermion helicity
	 	 
            call ket2c(psi(1,isig,i),.true.,p(0,i),isig,pww,aww,
     1                 psia(1,isig,is1,i),fq(0,i))
            call ket2c(psi(1,isig,i),.true.,p(0,i),isig,pww,zww,
     1                 psiz(1,isig,is1,i),fq(0,i))
     
          
 ! gluon radiation from fermion i 
 	    do l = 1,2
            	call ket2c(braketg(1,isig,i,l),.false.,pq(0,i),
     $                 isig,pww,aww,braketgA(1,isig,l),pga(0,i))
                                            ! |aww,q,i>_l,isigi
           	call ket2c(braketg(1,isig,i,l),.false.,pq(0,i),
     $                 isig,pww,zww,braketgZ(1,isig,l),pgz(0,i))
                                            ! |zww,q,i>_l,isigi
            	call ket2r(psia(1,isig,is1,i),.false.,fq(0,i),isig,
     $	    		q,eps(0,l),braketAg(1,isig,l),pag(0,i))      
     					    ! |q,aww,i>_l,isigi
            	call ket2r(psiz(1,isig,is1,i),.false.,fq(0,i),isig,
     $	    		q,eps(0,l),braketZg(1,isig,l),pzg(0,i))      
     					    ! |q,zww,i>_l,isigi
	    enddo
	    	 
	 enddo ! isig loop   			    					    
			
							    
 	 do l = 1,2
   ! gluon radiation from fermion i 
     	    call curr6(1,psi(1,-1,i+1),p(0,i+1),
     $	    				braketgA(1,-1,l),pga(0,i),jag) 			
      	    call curr6(1,psi(1,-1,i+1),p(0,i+1),
     $			        	braketAg(1,-1,l),pag(0,i),jga)
  ! gluon radiation from fermion i+1 
     	    call curr6(1,braketg(1,-1,i+1,l),pq(0,i+1),
     $	    				psia(1,-1,is1,i),fq(0,i),jg0)
     
          		    
      	    call curr6(1,psi(1,-1,i+1),p(0,i+1),
     $	    				braketgZ(1,-1,l),pgz(0,i),jzg)  		       
       	    call curr6(1,psi(1,-1,i+1),p(0,i+1),
     $	    				braketZg(1,-1,l),pzg(0,i),jgz)  		       
  ! gluon radiation from fermion i+1
     	    call curr6(1,braketg(1,-1,i+1,l),pq(0,i+1),
     $	    				psiz(1,-1,is1,i),fq(0,i),jgm)

       	    do mu = 0,5
	    	do isig = -1,1,2
	    	   
 		   ja(mu,isig,is1,i,l) = jag(mu,isig)+
     $			jga(mu,isig)+jg0(mu,isig)  ! A+g emission from i/i+1 line 
     	    	   
		   jz(mu,isig,is1,i,l) = jzg(mu,isig)+
     $	    		jgz(mu,isig)+jgm(mu,isig)  ! Z+g emission from i/i+1 line 
		
		enddo !isig
	    enddo  !mu

   	 enddo ! l = 1,2


  ! gluon radiation from the non-i line:
  	 do l = 1,2
            call curr6(1,psi(1,-1,i+1),p(0,i+1),
     1              psia(1,-1,is1,i),fq(0,i), ja(0,-1,is1,i,0) )
     
            call curr6(1,psi(1,-1,i+1),p(0,i+1),
     1              psiz(1,-1,is1,i),fq(0,i), jz(0,-1,is1,i,0) )
     	 enddo
        
            
c ------------------------------------      
      
c attach V=A/Z to f3:

     
      i = 3
      
         do isig = -1,1,2 ! fermion helicity
	 
            call ket2c(psi(1,isig,i),.true.,p(0,i),isig,pww,aww,
     1                 psia(1,isig,is3,i),fq(0,i))
            call ket2c(psi(1,isig,i),.true.,p(0,i),isig,pww,zww,
     1                 psiz(1,isig,is3,i),fq(0,i))
     
 ! gluon radiation from fermion i 
 	    do l = 1,2
            	call ket2c(braketg(1,isig,i,l),.false.,pq(0,i),
     $                 isig,pww,aww,braketgA(1,isig,l),pga(0,i))
                                            ! |aww,q,i>_l,isigi
           	call ket2c(braketg(1,isig,i,l),.false.,pq(0,i),
     $                 isig,pww,zww,braketgZ(1,isig,l),pgz(0,i))
                                            ! |zww,q,i>_l,isigi
            	call ket2r(psia(1,isig,is3,i),.false.,fq(0,i),isig,
     $	    		q,eps(0,l),braketAg(1,isig,l),pag(0,i))      
     					    ! |q,aww,i>_l,isigi
            	call ket2r(psiz(1,isig,is3,i),.false.,fq(0,i),isig,
     $	    		q,eps(0,l),braketZg(1,isig,l),pzg(0,i))      
     					    ! |q,zww,i>_l,isigi
	    enddo			    					    

	 enddo ! isig loop   			    					    
			
							    
 	 do l = 1,2
      	    call curr6(1,psi(1,-1,i+1),p(0,i+1),
     $	    				braketgA(1,-1,l),pga(0,i),jag)			   
      	    call curr6(1,psi(1,-1,i+1),p(0,i+1),
     $	    				braketAg(1,-1,l),pag(0,i),jga)
  ! gluon radiation from fermion i+1 
     	    call curr6(1,braketg(1,-1,i+1,l),pq(0,i+1),
     $	    				psia(1,-1,is3,i),fq(0,i),jg0)


      	    call curr6(1,psi(1,-1,i+1),p(0,i+1),
     $	    				braketgZ(1,-1,l),pgz(0,i),jzg)  		       
       	    call curr6(1,psi(1,-1,i+1),p(0,i+1),
     $	    				braketZg(1,-1,l),pzg(0,i),jgz)  		       
  ! gluon radiation from fermion i+1
     	    call curr6(1,braketg(1,-1,i+1,l),pq(0,i+1),
     $	    				psiz(1,-1,is3,i),fq(0,i),jgm)

       	    do mu = 0,5
	    	do isig = -1,1,2
		
      	    	   ja(mu,isig,is3,i,l) = jag(mu,isig)+
     $			jga(mu,isig)+jg0(mu,isig)  ! A+g emission from i/i+1 line 
		
       	    	   jz(mu,isig,is3,i,l) = jzg(mu,isig)+
     $	    		jgz(mu,isig)+jgm(mu,isig)  ! Z+g emission from i/i+1 line 
	   	
		enddo ! isig 
	    enddo ! mu	

   	 enddo  ! l=1,2 loop


  ! gluon radiation from the non-i line:
  	 do l = 1,2
            call curr6(1,psi(1,-1,i+1),p(0,i+1),
     1              psia(1,-1,is3,i),fq(0,i), ja(0,-1,is3,i,0) )
     
            call curr6(1,psi(1,-1,i+1),p(0,i+1),
     1              psiz(1,-1,is3,i),fq(0,i), jz(0,-1,is3,i,0) )
     	 enddo
        
      
c      
c  
c ------------------------------------      
c        
c attach V=A/Z to f2:
  
     
      i=1
      
         do isig = -1,1,2 ! fermion helicity
	 
            call bra2c(psi(1,isig,i+1),.true.,p(0,i+1),isig,pww,aww,
     1                 psia(1,isig,is1,i+1),fq(0,i+1))
            call bra2c(psi(1,isig,i+1),.true.,p(0,i+1),isig,pww,zww,
     1                 psiz(1,isig,is1,i+1),fq(0,i+1))
     
 ! gluon radiation from fermion i+1 
 	    do l = 1,2
            	call bra2c(braketg(1,isig,i+1,l),.false.,pq(0,i+1),
     $                 isig,pww,aww,braketgA(1,isig,l),pga(0,i+1))
                                            ! <aww,q,i+1|_l,isig
           	call bra2c(braketg(1,isig,i+1,l),.false.,pq(0,i+1),
     $                 isig,pww,zww,braketgZ(1,isig,l),pgz(0,i+1))
                                            ! <zww,q,i+1|_l,isig
            	call bra2r(psia(1,isig,is1,i+1),.false.,fq(0,i+1),isig,
     $	    		q,eps(0,l),braketAg(1,isig,l),pag(0,i+1))      
     					    ! <q,aww,i+1|_l,isig
            	call bra2r(psiz(1,isig,is1,i+1),.false.,fq(0,i+1),isig,
     $	    		q,eps(0,l),braketZg(1,isig,l),pzg(0,i+1))      
     					    ! <q,zww,i+1|_l,isig
	    enddo			    					    

	 enddo ! isig loop   			    					    
			
							    
 	 do l = 1,2
      	    call curr6(1,braketgA(1,-1,l),pga(0,i+1),
     $	    			psi(1,-1,i),p(0,i),jag)		       
      	    call curr6(1,braketAg(1,-1,l),pag(0,i+1),
     $			        psi(1,-1,i),p(0,i),jga)
  ! gluon radiation from fermion i 
     	    call curr6(1,psia(1,-1,is1,i+1),fq(0,i+1),
     $	    			braketg(1,-1,i,l),pq(0,i),jg0)


      	    call curr6(1,braketgZ(1,-1,l),pgz(0,i+1),
     $	    				psi(1,-1,i),p(0,i),jzg)		       
       	    call curr6(1,braketZg(1,-1,l),pzg(0,i+1),
     $	    				psi(1,-1,i),p(0,i),jgz)		       
  ! gluon radiation from fermion i
     	    call curr6(1,psiz(1,-1,is1,i+1),fq(0,i+1),
     $	    			braketg(1,-1,i,l),pq(0,i),jgm)

       	    do mu = 0,5
	    	do isig = -1,1,2
      	    
	    	   ja(mu,isig,is1,i+1,l) = jag(mu,isig)+
     $			jga(mu,isig)+jg0(mu,isig)  ! A+g emission from i/i+1 line 
      	    	   jz(mu,isig,is1,i+1,l) = jzg(mu,isig)+
     $	    		jgz(mu,isig)+jgm(mu,isig)  ! Z+g emission from i/i+1 line 
	    
	    	enddo ! isig 
	    enddo ! mu					 

   	 enddo ! l=1,2 loop

  ! gluon radiation from the non-i line:
  	 do l = 1,2
            call curr6(1,psia(1,-1,is1,i+1),fq(0,i+1),
     1              psi(1,-1,i),p(0,i), ja(0,-1,is1,i+1,0) )
     
            call curr6(1,psiz(1,-1,is1,i+1),fq(0,i+1),
     1              psi(1,-1,i),p(0,i), jz(0,-1,is1,i+1,0) )
     	 enddo


c ------------------------------------      
c
c attach V=A/Z to f4:
  
      
      i=3
      
         do isig = -1,1,2 ! fermion helicity
	 
            call bra2c(psi(1,isig,i+1),.true.,p(0,i+1),isig,pww,aww,
     1                 psia(1,isig,is3,i+1),fq(0,i+1))
            call bra2c(psi(1,isig,i+1),.true.,p(0,i+1),isig,pww,zww,
     1                 psiz(1,isig,is3,i+1),fq(0,i+1))
     
 ! gluon radiation from fermion i+1 
 	    do l = 1,2
            	call bra2c(braketg(1,isig,i+1,l),.false.,pq(0,i+1),
     $                 isig,pww,aww,braketgA(1,isig,l),pga(0,i+1))
                                            ! <aww,q,i+1|_l,isig
           	call bra2c(braketg(1,isig,i+1,l),.false.,pq(0,i+1),
     $                 isig,pww,zww,braketgZ(1,isig,l),pgz(0,i+1))
                                            ! <zww,q,i+1|_l,isig
            	call bra2r(psia(1,isig,is3,i+1),.false.,fq(0,i+1),isig,
     $	    		q,eps(0,l),braketAg(1,isig,l),pag(0,i+1))      
     					    ! <q,aww,i+1|_l,isig
            	call bra2r(psiz(1,isig,is3,i+1),.false.,fq(0,i+1),isig,
     $	    		q,eps(0,l),braketZg(1,isig,l),pzg(0,i+1))      
     					    ! <q,zww,i+1|_l,isig
	    enddo			    					    
			
	 enddo ! isig loop   			    					    
							    
 	 do l = 1,2
      	    call curr6(1,braketgA(1,-1,l),pga(0,i+1),
     $	    				psi(1,-1,i),p(0,i),jag)		       
      	    call curr6(1,braketAg(1,-1,l),pag(0,i+1),
     $	    				psi(1,-1,i),p(0,i),jga)
  ! gluon radiation from fermion i 
     	    call curr6(1,psia(1,-1,is3,i+1),fq(0,i+1),
     $	    			braketg(1,-1,i,l),pq(0,i),jg0)

  
      	    call curr6(1,braketgZ(1,-1,l),pgz(0,i+1),
     $	    				psi(1,-1,i),p(0,i),jzg)		       
       	    call curr6(1,braketZg(1,-1,l),pzg(0,i+1),
     $	    				psi(1,-1,i),p(0,i),jgz)		       
  ! gluon radiation from fermion i
     	    call curr6(1,psiz(1,-1,is3,i+1),fq(0,i+1),
     $	    			braketg(1,-1,i,l),pq(0,i),jgm)

      	    do mu = 0,5
	    	do isig = -1,1,2
		
      	    	   ja(mu,isig,is3,i+1,l) = jag(mu,isig)+
     $			jga(mu,isig)+jg0(mu,isig)  ! A+g emission from i/i+1 line 
      	    	   jz(mu,isig,is3,i+1,l) = jzg(mu,isig)+
     $	    		jgz(mu,isig)+jgm(mu,isig)  ! Z+g emission from i/i+1 line 
	    	
		enddo	!isig
	    enddo !mu	

   	 enddo ! l=1,2 loop
	 					 

  ! gluon radiation from the non-i line:
  	 do l = 1,2
            call curr6(1,psia(1,-1,is3,i+1),fq(0,i+1),
     1              psi(1,-1,i),p(0,i), ja(0,-1,is3,i+1,0) )
     
            call curr6(1,psiz(1,-1,is3,i+1),fq(0,i+1),
     1              psi(1,-1,i),p(0,i), jz(0,-1,is3,i+1,0) )
     	 enddo
        
    
c -----------------------------------------------------------------------
           
      
      do isig1 = -1,1,2
        do isig3 = -1,1,2
     
       	  do l = 1,2
            do i = 1,2
	
c	A/Z from fermion i (=1,2) and g upper line:
c		ja and jz(0,isig1,is1,i,l)
		
c	A/Z from fermion i (=1,2 upper line), g lower line:
c		ja and jz(0,isig1,is1,i,0)
		
c	g upper line, A/Z from fermion i (= 3,4 lower line):
c	        ja and jz(0,isig3,is3,2+i,0)
	
c	g lower line and A/Z from fermion i (= 3,4 lower line):
c	        ja and jz(0,isig3,is3,2+i,l)
	
	
	! ma(i,j,l) for emission of A/Z from fermion i (i = 1,2,3,4), 
	!	and gluon with pol l=1,2 from line 12 (j=2) or 34 (j=3):
		
		ma(i,2,l) = dotcc(ja(0,isig1,is1,i,l),
     #					jqq(0,isig3,2,is3,0))
		ma(i,3,l) = dotcc(ja(0,isig1,is1,i,0),
     #					jqq(0,isig3,2,is3,l))	
		ma(i+2,2,l) = dotcc(jqq(0,isig1,1,is1,l),
     #					ja(0,isig3,is3,i+2,0))
		ma(i+2,3,l) = dotcc(jqq(0,isig1,1,is1,0),
     #					ja(0,isig3,is3,i+2,l))
		
		mz(i,2,l) = dotcc(jz(0,isig1,is1,i,l),
     #					jqq(0,isig3,2,is3,0))
		mz(i,3,l) = dotcc(jz(0,isig1,is1,i,0),
     #					jqq(0,isig3,2,is3,l))	
		mz(i+2,2,l) = dotcc(jqq(0,isig1,1,is1,l),
     #					jz(0,isig3,is3,i+2,0))
		mz(i+2,3,l) = dotcc(jqq(0,isig1,1,is1,0),
     #					jz(0,isig3,is3,i+2,l))

	
	    enddo  !i=1,2 loop
		

c --------------------------------------------------------------------
c
c  NEUTRAL CURRENT:
c
c	A/Z prop for aww piece contained in aww tensor,
c	prop43/21(1/2,j) takes care of weak boson exchange
	
	 do j = jmin,jmax
           do k = 1,4	
	       ! V from upper line:
              propt(isig1,isig3,k,j,2) = 
     1   	 clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3)*prop43(1,j) 
     2       +   clr(ifl(1,k),2,isig1)*clr(ifl(3,k),2,isig3)*prop43(2,j)
              
	      mata(k,isig1,isig3,j,l,4) = propt(isig1,isig3,k,j,2) *
     1   	 ( (ma(1,j,l)+ma(2,j,l))*clr(ifl(1,k),1,isig1) )
     
	      matz(k,isig1,isig3,j,l,4) = propt(isig1,isig3,k,j,2) *
     1   	 ( (mz(1,j,l)+mz(2,j,l))*clr(ifl(1,k),2,isig1)   )
     
	      mat(k,isig1,isig3,j,l,4) = mata(k,isig1,isig3,j,l,4)+
     1	      				 matz(k,isig1,isig3,j,l,4)
     

 	       ! V from lower line:
             propt(isig1,isig3,k,j,1) = 
     1   	 clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3)*prop21(1,j) 
     2       +   clr(ifl(1,k),2,isig1)*clr(ifl(3,k),2,isig3)*prop21(2,j)
              
	      mata(k,isig1,isig3,j,l,5) = propt(isig1,isig3,k,j,1) *
     1   	 ( (ma(3,j,l)+ma(4,j,l))*clr(ifl(3,k),1,isig3))
     
              matz(k,isig1,isig3,j,l,5) = propt(isig1,isig3,k,j,1) *
     1   	 ( (mz(3,j,l)+mz(4,j,l))*clr(ifl(3,k),2,isig3)   )
     
              mat(k,isig1,isig3,j,l,5) = mata(k,isig1,isig3,j,l,5) + 
     1	      				 matz(k,isig1,isig3,j,l,5)
	      
         enddo ! k loop
	 enddo  ! j loop	

c	enddo  ! l loop
	              
c -----------------------------------------------------------
c
c  CHARGED CURRENT:

	  if (isig1.eq.-1 .and. isig3.eq.-1) then
	    do j = jmin,jmax
             do k = 5,6
	     ! A/Z from upper line:
         	propt(-1,-1,k,j,2) = clr(ifl(1,k),3,-1)**2*prop43(3,j)
         	
	    mata(k,-1,-1,j,l,4) = propt(-1,-1,k,j,2) *
     1       (ma(2,j,l)*clr(ifl(2,k),1,-1)+ma(1,j,l)*clr(ifl(1,k),1,-1))		    
         	
	    matz(k,-1,-1,j,l,4) = propt(-1,-1,k,j,2) *
     1       (mz(2,j,l)*clr(ifl(2,k),2,-1)+mz(1,j,l)*clr(ifl(1,k),2,-1))		    
         	
	    mat(k,-1,-1,j,l,4) = mata(k,-1,-1,j,l,4)+matz(k,-1,-1,j,l,4)   
          
	     enddo ! k
	    enddo ! j
          else
	    do j = jmin,jmax
             do k = 5,6
         	propt(isig1,isig3,k,j,2) = 0d0
         	mat(k,isig1,isig3,j,l,4) = 0d0
             enddo ! k
	    enddo ! j
          endif


c -----------------------------------------------------------           
	    
	  if (isig1.eq.-1 .and. isig3.eq.-1) then
            do j = jmin,jmax
	     do k = 5,6
	     ! A/Z from lower line:
         	 propt(-1,-1,k,j,1) = clr(ifl(1,k),3,-1)**2*prop21(3,j)
            
	    mata(k,-1,-1,j,l,5) = propt(-1,-1,k,j,1) *
     1       (ma(4,j,l)*clr(ifl(4,k),1,-1)+ma(3,j,l)*clr(ifl(3,k),1,-1))		  
 
            matz(k,-1,-1,j,l,5) = propt(-1,-1,k,j,1) *
     1       (mz(4,j,l)*clr(ifl(4,k),2,-1)+mz(3,j,l)*clr(ifl(3,k),2,-1))		  
             
	    mat(k,-1,-1,j,l,5) = mata(k,-1,-1,j,l,5)+matz(k,-1,-1,j,l,5)
		
	     enddo ! k
	    enddo ! j
          else
            do j = jmin,jmax
             do k = 5,6
         	propt(isig1,isig3,k,j,1) = 0d0
         	mat(k,isig1,isig3,j,l,5) = 0d0
             enddo
	    enddo ! j
          endif
	      
	      enddo ! l loop

	  enddo ! isig3
      enddo ! isig1	
      	
c ----------------------------------------------------------------    
c
c check result for NC as obtained above against madgraph output 
c		(sum up relevant mg amplitudes in the following):

       if (ldebugm) then
c        do k = 1,1
        do k = 1,2
	  do lh = -1,1,2	! gluon pol. in hel. basis
 	    
	! need to distingiush between j=2,3 to account for crossed processes;
	! in the way computed here only uncrossed processes with outgoing 
	!	gluons can be cecked!
		
c compute reference result from individual madgraph diagrams tampr(i,...):
	
	  m1hc(lh) = 0d0

	  if(k.eq.1) then
	    
	    do  i = 19,48
		m1hc(lh) = m1hc(lh) + tampr(i,-1,-1,lh,k) 
	    enddo		
	    do  i = 86,135
		m1hc(lh) = m1hc(lh) + tampr(i,-1,-1,lh,k) 
	    enddo	
	    do  i = 154,183
		m1hc(lh) = m1hc(lh) + tampr(i,-1,-1,lh,k) 
	    enddo	
	    do  i = 221,270 
		m1hc(lh) = m1hc(lh) + tampr(i,-1,-1,lh,k) 
	    enddo	
	    do  i = 283,362 
		m1hc(lh) = m1hc(lh) + tampr(i,-1,-1,lh,k) 
	    enddo	
	    do  i = 375,454 
		m1hc(lh) = m1hc(lh) + tampr(i,-1,-1,lh,k) 
	    enddo
	    
	    
	  elseif (k.eq.2) then
	  
	    do  i = 19,48
		m1hc(lh) = m1hc(lh) + tampr(i,-1,-1,lh,k) 
	    enddo		
	    do  i = 86,135
		m1hc(lh) = m1hc(lh) + tampr(i,-1,-1,lh,k) 
	    enddo	
	    do  i = 148,227
		m1hc(lh) = m1hc(lh) + tampr(i,-1,-1,lh,k) 
	    enddo	
	    do  i = 240,319
		m1hc(lh) = m1hc(lh) + tampr(i,-1,-1,lh,k) 
	    enddo	
	    do  i = 338,367
		m1hc(lh) = m1hc(lh) + tampr(i,-1,-1,lh,k) 
	    enddo	
	    do  i =  405,454
		m1hc(lh) = m1hc(lh) + tampr(i,-1,-1,lh,k) 
	    enddo
	  
	  endif  !k
	
			 
	  enddo ! end of lh loop 	  
	    
		    
    ! compute m1 in kartesian basis:
	    
	    do l = 1,2
	      m1kc(l) = 0d0 
	    enddo  	    	    
	   
	    m1kc(1) =	  (m1hc(-1)-m1hc(1))/sqrt(2d0)
	    m1kc(2) = -im*(m1hc(-1)+m1hc(1))/sqrt(2d0)
c
c
c ---------------------------------		 
c		 		    
c perform checks and print results:

	    print*,"check for sum of j=2,3" 
 	    print*,"av/zv debugging:"	
	    print*,"signs:",sign	
	    print*,"gign:",gsign	    
  	    print*," vbf for k=",k,"      "
	    
	    
      	    print*,"A&Z from upper or lower line for A/Z->WW"
 	    print*,"ratio for l=1:",
     #	    	(mat(k,-1,-1,2,1,4)+mat(k,-1,-1,3,1,4)+
     #	    	 mat(k,-1,-1,2,1,5)+mat(k,-1,-1,3,1,5))/
     #		(m1kc(1))
      	    print*,"ratio for l=2:",
     #	    	(mat(k,-1,-1,2,2,4)+mat(k,-1,-1,3,2,4)+
     #	    	 mat(k,-1,-1,2,2,5)+mat(k,-1,-1,3,2,5))/
     #		(m1kc(2))
    	    print*,"res for l=1:",
c     #	    	(mat(k,-1,-1,2,1,4)+mat(k,-1,-1,3,1,4)),
c     #	    	(mat(k,-1,-1,2,1,5)+mat(k,-1,-1,3,1,5)),
     #	    	(mat(k,-1,-1,2,1,4)+mat(k,-1,-1,3,1,4)+
     #	    	 mat(k,-1,-1,2,1,5)+mat(k,-1,-1,3,1,5)),
     #		(m1kc(1))
      	    print*,"res for l=2:",
c     #	    	(mat(k,-1,-1,2,2,4)+mat(k,-1,-1,3,2,4)),
c     #	    	(mat(k,-1,-1,2,2,5)+mat(k,-1,-1,3,2,5)),
     #	    	(mat(k,-1,-1,2,2,4)+mat(k,-1,-1,3,2,4)+
     #	    	 mat(k,-1,-1,2,2,5)+mat(k,-1,-1,3,2,5)),
     #		(m1kc(2))
	    print*
	

        enddo ! end of k loop
      endif ! end of debugging block
      
c	fpials has to be attached below!      
      
c
c -----------------------------------------------------------------------
c ----------------------------------------------------------------    
c
c check result for CC as obtained above against madgraph output 
c		(sum up relevant mg amplitudes in the following):

      if (ldebugm) then
c        do k = 5,5
        do k = 5,6
	  do lh = -1,1,2	! gluon pol. in hel. basis
 	    
	! need to distingiush between j=2,3
		
c compute reference result from individual madgraph diagrams tampr(i,...):
	
	    do j = 2,3
	    	m1haj(lh,j) = 0d0  ! graphs with A (->WW) emission
	    	m1hzj(lh,j) = 0d0  !		 Z (->WW) emission
            enddo  

	    do m = 1,30
		i = lista2(m)
		m1haj(lh,2) = m1haj(lh,2) + tampr(i,-1,-1,lh,k) 
	    enddo
		
	    do m = 1,30
		i = lista3(m)
		m1haj(lh,3) = m1haj(lh,3) + tampr(i,-1,-1,lh,k) 
	    enddo
	
	    do m = 1,50
		i = listz2(m)
		m1hzj(lh,2) = m1hzj(lh,2) + tampr(i,-1,-1,lh,k) 
	    enddo
	
	    do m = 1,50
		i = listz3(m)
		m1hzj(lh,3) = m1hzj(lh,3) + tampr(i,-1,-1,lh,k) 
	    enddo
	
	     if (k.eq.5) then  ! for incoming g: no g emission from other line
       		if(.not.jlog2) m1haj(lh,2) = 0d0
       		if(.not.jlog3) m1haj(lh,3) = 0d0
       		if(.not.jlog2) m1hzj(lh,2) = 0d0
       		if(.not.jlog3) m1hzj(lh,3) = 0d0
	     else
       		if(.not.jlog2) m1haj(lh,3) = 0d0
       		if(.not.jlog3) m1haj(lh,2) = 0d0
       		if(.not.jlog2) m1hzj(lh,3) = 0d0
       		if(.not.jlog3) m1hzj(lh,2) = 0d0
             endif
			 
	  enddo ! end of lh loop 	  
	    
		    
    ! compute m1 in kartesian basis:
	    
	    do l = 1,2
	      do j = 2,3
	    	m1kaj(l,j) = 0d0
	    	m1kzj(l,j) = 0d0
	      enddo
	    enddo  	    
	    
	  do j = 2,3
	    m1kaj(1,j) =     (m1haj(-1,j)-m1haj(1,j))/sqrt(2d0)
	    m1kaj(2,j) = -im*(m1haj(-1,j)+m1haj(1,j))/sqrt(2d0)
	    m1kzj(1,j) =     (m1hzj(-1,j)-m1hzj(1,j))/sqrt(2d0)
	    m1kzj(2,j) = -im*(m1hzj(-1,j)+m1hzj(1,j))/sqrt(2d0)
	  enddo    	  
c
c
c ---------------------------------		 
c		 		    
c perform checks and print results:


 	if(.true.) then
	    print*,"check for sum of j=2,3" 
 	    print*,"av/zv debugging:"	
	    print*,"signs:",sign	
	    print*,"gign:",gsign	    
  	    print*," vbf for k=",k,"      "
	    
	    
      	    print*,"A/Z from upper or lower line"
 	    print*,"ratio for l=1:",
     #	    	(mat(k,-1,-1,2,1,4)+mat(k,-1,-1,3,1,4)+
     #	    	 mat(k,-1,-1,2,1,5)+mat(k,-1,-1,3,1,5))/
     #		(m1kaj(1,2)+m1kaj(1,3)+m1kzj(1,2)+m1kzj(1,3))
      	    print*,"ratio for l=2:",
     #	    	(mat(k,-1,-1,2,2,4)+mat(k,-1,-1,3,2,4)+
     #	    	 mat(k,-1,-1,2,2,5)+mat(k,-1,-1,3,2,5))/
     #		(m1kaj(2,2)+m1kaj(2,3)+m1kzj(2,2)+m1kzj(2,3))
    	    print*,"res for l=1:",
c     #	    	(mat(k,-1,-1,2,1,4)+mat(k,-1,-1,3,1,4)),
c     #	    	(mat(k,-1,-1,2,1,5)+mat(k,-1,-1,3,1,5)),
     #	    	(mat(k,-1,-1,2,1,4)+mat(k,-1,-1,3,1,4)+
     #	    	 mat(k,-1,-1,2,1,5)+mat(k,-1,-1,3,1,5)),
     #		(m1kaj(1,2)+m1kaj(1,3)+m1kzj(1,2)+m1kzj(1,3))
      	    print*,"res for l=2:",
c     #	    	(mat(k,-1,-1,2,2,4)+mat(k,-1,-1,3,2,4)),
c     #	    	(mat(k,-1,-1,2,2,5)+mat(k,-1,-1,3,2,5)),
     #	    	(mat(k,-1,-1,2,2,4)+mat(k,-1,-1,3,2,4)+
     #	    	 mat(k,-1,-1,2,2,5)+mat(k,-1,-1,3,2,5)),
     #		(m1kaj(2,2)+m1kaj(2,3)+m1kzj(2,2)+m1kzj(2,3))
	    print*
	
	if(.false.) then
     	    print*,"Z from upper or lower line"
  	    print*,"ratio for l=1:",
     #	    	(matz(k,-1,-1,2,1,4)+matz(k,-1,-1,3,1,4)+
     #	    	 matz(k,-1,-1,2,1,5)+matz(k,-1,-1,3,1,5))/
     #		(m1kzj(1,2)+m1kzj(1,3))
      	    print*,"ratio for l=2:",
     #	    	(matz(k,-1,-1,2,2,4)+matz(k,-1,-1,3,2,4)+
     #	    	 matz(k,-1,-1,2,2,5)+matz(k,-1,-1,3,2,5))/
     #		(m1kzj(2,2)+m1kzj(2,3))
    	    print*,"res for l=1:",
c     #	    	(matz(k,-1,-1,2,1,4)+matz(k,-1,-1,3,1,4)),
c     #	    	(matz(k,-1,-1,2,1,5)+matz(k,-1,-1,3,1,5)),
     #	    	(matz(k,-1,-1,2,1,4)+matz(k,-1,-1,3,1,4)+
     #	    	 matz(k,-1,-1,2,1,5)+matz(k,-1,-1,3,1,5)),
     #		(m1kzj(1,2)+m1kzj(1,3))
      	    print*,"res for l=2:",
c     #	    	(matz(k,-1,-1,2,2,4)+matz(k,-1,-1,3,2,4)),
c     #	    	(matz(k,-1,-1,2,2,5)+matz(k,-1,-1,3,2,5)),
     #	    	(matz(k,-1,-1,2,2,4)+matz(k,-1,-1,3,2,4)+
     #	    	 matz(k,-1,-1,2,2,5)+matz(k,-1,-1,3,2,5)),
     #		(m1kzj(2,2)+m1kzj(2,3))
	    print*
	endif    
	
	endif
c	   

        enddo ! end of k loop
      endif ! end of debugging block
      
c	fpials has to be attached below!      
      
c
c -----------------------------------------------------------------------
c -----------------------------------------------------------------------
c
c  next do the box-box graphs with one W emitted from the upper and the 
c  other from the lower line. These are only possible for lefthanded quarks
c  on both lines, i.e. isig1 = -1 and isig3 = -1
      
      do l = 1,2
         j = 2	!g from upper line 
	
c  upper line box: eps1 = wp, eps2 = jwm(mu,-1,3)
c  lower line box: eps1 = wm, eps2 = jwp(mu,-1,1)
      mpm(1,j,l) = -dotcc(jwp(0,-1,is1,1,l),jwm(0,-1,is3,3,0))  !- zpm
      mmp(1,j,l) = -dotcc(jwm(0,-1,is1,1,l),jwp(0,-1,is3,3,0))  ! -zmp
c  upper line box: eps1 = wp, eps2 = jwm(mu,-1,4)
c  lower line box: eps1 = jwp(mu,-1,1), eps2 = wm
      mpm(2,j,l) = -dotcc(jwp(0,-1,is1,1,l),jwm(0,-1,is3,4,0)) ! + zpm
      mmp(2,j,l) = -dotcc(jwm(0,-1,is1,1,l),jwp(0,-1,is3,4,0)) ! + zmp
c  upper line box: eps1 = jwm(mu,-1,3), eps2 = wp
c  lower line box: eps1 = wm, eps2 = jwp(mu,-1,2)
      mpm(3,j,l) = -dotcc(jwp(0,-1,is1,2,l),jwm(0,-1,is3,3,0)) ! + zpm 
      mmp(3,j,l) = -dotcc(jwm(0,-1,is1,2,l),jwp(0,-1,is3,3,0)) ! + zmp
c  upper line box: eps1 = jwm(mu,-1,4), eps2 = wp
c  lower line box: eps1 = jwp(mu,-1,2), eps2 = wm
      mpm(4,j,l) = -dotcc(jwp(0,-1,is1,2,l),jwm(0,-1,is3,4,0)) ! - zpm
      mmp(4,j,l) = -dotcc(jwm(0,-1,is1,2,l),jwp(0,-1,is3,4,0)) ! - zmp
 
         j = 3	!g from lower line 

      mpm(1,j,l) = -dotcc(jwp(0,-1,is1,1,0),jwm(0,-1,is3,3,l))  !- zpm
      mmp(1,j,l) = -dotcc(jwm(0,-1,is1,1,0),jwp(0,-1,is3,3,l))  ! -zmp
c  upper line box: eps1 = wp, eps2 = jwm(mu,-1,4)
c  lower line box: eps1 = jwp(mu,-1,1), eps2 = wm
      mpm(2,j,l) = -dotcc(jwp(0,-1,is1,1,0),jwm(0,-1,is3,4,l)) ! + zpm
      mmp(2,j,l) = -dotcc(jwm(0,-1,is1,1,0),jwp(0,-1,is3,4,l)) ! + zmp
c  upper line box: eps1 = jwm(mu,-1,3), eps2 = wp
c  lower line box: eps1 = wm, eps2 = jwp(mu,-1,2)
      mpm(3,j,l) = -dotcc(jwp(0,-1,is1,2,0),jwm(0,-1,is3,3,l)) ! + zpm 
      mmp(3,j,l) = -dotcc(jwm(0,-1,is1,2,0),jwp(0,-1,is3,3,l)) ! + zmp
c  upper line box: eps1 = jwm(mu,-1,4), eps2 = wp
c  lower line box: eps1 = jwp(mu,-1,2), eps2 = wm
      mpm(4,j,l) = -dotcc(jwp(0,-1,is1,2,0),jwm(0,-1,is3,4,l)) ! - zpm
      mmp(4,j,l) = -dotcc(jwm(0,-1,is1,2,0),jwp(0,-1,is3,4,l)) ! - zmp

      enddo !l

c ------------------------------

c  for the q^mu*q^nu/M_V^2 terms in the gauge boson propagators we need zpm/zmp
c
c
c neutral current:      

      do l = 1,2
      	
      zpm(2,l) = jj21p(l)*jj43m(0)*zm2i(3)
      zmp(2,l) = jj21m(l)*jj43p(0)*zm2i(3)
c
      zpm(3,l) = jj21p(0)*jj43m(l)*zm2i(3)
      zmp(3,l) = jj21m(0)*jj43p(l)*zm2i(3)


      do j = jmin,jmax

      k = 1
      mat(k,-1,-1,j,l,6) = (mpm(2,j,l)+zpm(j,l))*prop_pm(3,j) + 
     # 			   (mmp(3,j,l)+zmp(j,l))*prop_mp(3,j)

      k = 2
      mat(k,-1,-1,j,l,6) = (mpm(1,j,l)-zpm(j,l))*prop_pm(3,j) + 
     #  		   (mmp(4,j,l)-zmp(j,l))*prop_mp(3,j) 
      k = 3
      mat(k,-1,-1,j,l,6) = (mpm(4,j,l)-zpm(j,l))*prop_pm(3,j) + 
     # 			   (mmp(1,j,l)-zmp(j,l))*prop_mp(3,j) 
      k = 4
      mat(k,-1,-1,j,l,6) = (mpm(3,j,l)+zpm(j,l))*prop_pm(3,j) + 
     # 			   (mmp(2,j,l)+zmp(j,l))*prop_mp(3,j)
     
      enddo !j
      
      enddo !l 
c
c no debugging done yet for NC!
c      
c
c ---------
c
c charged current:      
      
      do l = 1,2
      	
      zpm(2,l) = jj21p(l)*jj43m(0)*zm2i(2)
      zmp(2,l) = jj21m(l)*jj43p(0)*zm2i(2)

      zpm(3,l) = jj21p(0)*jj43m(l)*zm2i(2)
      zmp(3,l) = jj21m(0)*jj43p(l)*zm2i(2)
      
     
      do j = jmin,jmax

      k = 5
      mat(k,-1,-1,j,l,6) = 
     1  mpm(1,j,l)*(prop_pm(1,j)*clr(ifl(2,k),1,-1)*clr(ifl(4,k),1,-1) + 
     2  	prop_pm(2,j)*clr(ifl(2,k),2,-1)*clr(ifl(4,k),2,-1) ) +
     1  mpm(2,j,l)*(prop_pm(1,j)*clr(ifl(2,k),1,-1)*clr(ifl(3,k),1,-1) + 
     2  	prop_pm(2,j)*clr(ifl(2,k),2,-1)*clr(ifl(3,k),2,-1) ) +
     1  mpm(3,j,l)*(prop_pm(1,j)*clr(ifl(1,k),1,-1)*clr(ifl(4,k),1,-1) + 
     2  	prop_pm(2,j)*clr(ifl(1,k),2,-1)*clr(ifl(4,k),2,-1) ) +
     1  mpm(4,j,l)*(prop_pm(1,j)*clr(ifl(1,k),1,-1)*clr(ifl(3,k),1,-1) + 
     2  	prop_pm(2,j)*clr(ifl(1,k),2,-1)*clr(ifl(3,k),2,-1) ) 
     3  + zpm(j,l)*prop_pm(2,j)*(clr(ifl(2,k),2,-1)-clr(ifl(1,k),2,-1))*
     4  		 (clr(ifl(3,k),2,-1)-clr(ifl(4,k),2,-1))

      mata(k,-1,-1,j,l,6) = 
     1  mpm(1,j,l)*(prop_pm(1,j)*clr(ifl(2,k),1,-1)*clr(ifl(4,k),1,-1))+
     1  mpm(2,j,l)*(prop_pm(1,j)*clr(ifl(2,k),1,-1)*clr(ifl(3,k),1,-1))+
     1  mpm(3,j,l)*(prop_pm(1,j)*clr(ifl(1,k),1,-1)*clr(ifl(4,k),1,-1))+
     1  mpm(4,j,l)*(prop_pm(1,j)*clr(ifl(1,k),1,-1)*clr(ifl(3,k),1,-1))
      
      matz(k,-1,-1,j,l,6) = 
     1  mpm(1,j,l)*(prop_pm(2,j)*clr(ifl(2,k),2,-1)*clr(ifl(4,k),2,-1))+
     1  mpm(2,j,l)*(prop_pm(2,j)*clr(ifl(2,k),2,-1)*clr(ifl(3,k),2,-1))+
     1  mpm(3,j,l)*(prop_pm(2,j)*clr(ifl(1,k),2,-1)*clr(ifl(4,k),2,-1))+
     1  mpm(4,j,l)*(prop_pm(2,j)*clr(ifl(1,k),2,-1)*clr(ifl(3,k),2,-1))
     3  + zpm(j,l)*prop_pm(2,j)*(clr(ifl(2,k),2,-1)-clr(ifl(1,k),2,-1))*
     4  		 (clr(ifl(3,k),2,-1)-clr(ifl(4,k),2,-1))     
     
      mat(k,-1,-1,j,l,6) = mata(k,-1,-1,j,l,6) + matz(k,-1,-1,j,l,6)


      k = 6
      mat(k,-1,-1,j,l,6) = 
     1  mmp(1,j,l)*(prop_mp(1,j)*clr(ifl(2,k),1,-1)*clr(ifl(4,k),1,-1) + 
     2  	prop_mp(2,j)*clr(ifl(2,k),2,-1)*clr(ifl(4,k),2,-1) ) +
     1  mmp(2,j,l)*(prop_mp(1,j)*clr(ifl(2,k),1,-1)*clr(ifl(3,k),1,-1) + 
     2  	prop_mp(2,j)*clr(ifl(2,k),2,-1)*clr(ifl(3,k),2,-1) ) +
     1  mmp(3,j,l)*(prop_mp(1,j)*clr(ifl(1,k),1,-1)*clr(ifl(4,k),1,-1) + 
     2  	prop_mp(2,j)*clr(ifl(1,k),2,-1)*clr(ifl(4,k),2,-1) ) +
     1  mmp(4,j,l)*(prop_mp(1,j)*clr(ifl(1,k),1,-1)*clr(ifl(3,k),1,-1) + 
     2  	prop_mp(2,j)*clr(ifl(1,k),2,-1)*clr(ifl(3,k),2,-1) ) 
     3  + zmp(j,l)*prop_mp(2,j)*(clr(ifl(2,k),2,-1)-clr(ifl(1,k),2,-1))*
     4  		 (clr(ifl(3,k),2,-1)-clr(ifl(4,k),2,-1))

      mata(k,-1,-1,j,l,6) = 
     1  mmp(1,j,l)*(prop_mp(1,j)*clr(ifl(2,k),1,-1)*clr(ifl(4,k),1,-1))+
     1  mmp(2,j,l)*(prop_mp(1,j)*clr(ifl(2,k),1,-1)*clr(ifl(3,k),1,-1))+
     1  mmp(3,j,l)*(prop_mp(1,j)*clr(ifl(1,k),1,-1)*clr(ifl(4,k),1,-1))+
     1  mmp(4,j,l)*(prop_mp(1,j)*clr(ifl(1,k),1,-1)*clr(ifl(3,k),1,-1))

      matz(k,-1,-1,j,l,6) = 
     1  mmp(1,j,l)*(prop_mp(2,j)*clr(ifl(2,k),2,-1)*clr(ifl(4,k),2,-1))+
     1  mmp(2,j,l)*(prop_mp(2,j)*clr(ifl(2,k),2,-1)*clr(ifl(3,k),2,-1))+
     1  mmp(3,j,l)*(prop_mp(2,j)*clr(ifl(1,k),2,-1)*clr(ifl(4,k),2,-1))+
     1  mmp(4,j,l)*(prop_mp(2,j)*clr(ifl(1,k),2,-1)*clr(ifl(3,k),2,-1))
     3  + zmp(j,l)*prop_mp(2,j)*(clr(ifl(2,k),2,-1)-clr(ifl(1,k),2,-1))*
     4  		 (clr(ifl(3,k),2,-1)-clr(ifl(4,k),2,-1))
      
       mat(k,-1,-1,j,l,6) = mata(k,-1,-1,j,l,6) + matz(k,-1,-1,j,l,6)


      enddo !j  
      enddo !l
c
c ----------------------------------------------------------------    
c
c debugging:
c
c check result for NC as obtained above against madgraph output 
c		(sum up relevant mg amplitudes in the following):

      if(ldebugm) then
c       do k = 1,1
       do k = 1,2
	
	  do lh = -1,1,2	! gluon pol. in hel. basis
 	    
	! need to distingiush between j=2,3
		
c compute reference result from individual madgraph diagrams tampr(i,...):
	
	    do j = 2,3
	    	m1hj(lh,j) = 0d0 
            enddo  
	    m1hc(lh) = 0d0
		
	    m1hj(lh,2) = tampr(13,-1,-1,lh,k)+tampr(15,-1,-1,lh,k)+
     #			 tampr(17,-1,-1,lh,k)+tampr(136,-1,-1,lh,k)+
     #			 tampr(138,-1,-1,lh,k)+tampr(141,-1,-1,lh,k)
	   
	    m1hj(lh,3) = tampr(14,-1,-1,lh,k)+tampr(16,-1,-1,lh,k)+
     #			 tampr(18,-1,-1,lh,k)+tampr(137,-1,-1,lh,k)+
     #			 tampr(139,-1,-1,lh,k)+tampr(140,-1,-1,lh,k)
     
     	    do  i = 1,6
		m1hc(lh) = m1hc(lh) + tampr(i,-1,-1,lh,k) 
	    enddo
	    do  i = 332,337
		m1hc(lh) = m1hc(lh) + tampr(i,-1,-1,lh,k) 
	    enddo


	    if (k.ne.3) then  ! for incoming g: no g emission from other line
		      if(.not.jlog2) m1hj(lh,2) = 0d0
		      if(.not.jlog3) m1hj(lh,3) = 0d0
	    else
		      if(.not.jlog2) m1hj(lh,3) = 0d0
		      if(.not.jlog3) m1hj(lh,2) = 0d0
	     endif


	  enddo ! end of lh loop 	  
	    
		    
    ! compute m1 in kartesian basis:
	    
	    do l = 1,2
	      do j = 2,3
	    	m1kj(l,j) = 0d0
	      enddo
	      m1kc(l) = 0d0
	    enddo  	    
	    
	  do j = 2,3
	    m1kj(1,j) =     (m1hj(-1,j)-m1hj(1,j))/sqrt(2d0)
	    m1kj(2,j) = -im*(m1hj(-1,j)+m1hj(1,j))/sqrt(2d0)
	  enddo    	
	  m1kc(1) =	(m1hc(-1)-m1hc(1))/sqrt(2d0)
	  m1kc(2) = -im*(m1hc(-1)+m1hc(1))/sqrt(2d0)
	  

c ---------------------------------		 

c       if (.false.) then
      if (.true.) then

c        do l = 1,1
        do l = 1,2 
	  
	  if (k.eq.1) then
	  
	  matot = m1kj(l,2)+m1kj(l,3)
	  ratot = mat(k,-1,-1,2,l,6) + mat(k,-1,-1,3,l,6)    
	 
	  elseif(k.eq.2) then
	  
	  matot = m1kc(l)
	  ratot = mat(k,-1,-1,2,l,6) + mat(k,-1,-1,3,l,6)    
	  
	  endif
	  
	  
          print*,"TEST A & Z contributions for box-box:"
	    print*,"signs:",sign	
	    print*,"gign:",gsign	    
          print*," k,l = ",k,l
	  print*,"jmin,jmax",jmin,jmax
          print*,"follows: coded mat,mg"
	  print*
	  print*,"rat: ",ratot/matot
    	  print*,"res: ",ratot,matot

 	  print*


 	enddo !l
      endif ! "local" debugging block for A&Z

c ---------------------------------
    
      enddo ! k
 
c ---------------------------------
      
c check result for CC as obtained above against madgraph output 
c		(sum up relevant mg amplitudes in the following):

      do k = 5,6
c       do k = 6,6
	
	  do lh = -1,1,2	! gluon pol. in hel. basis
 	    
	! need to distingiush between j=2,3
		
c compute reference result from individual madgraph diagrams tampr(i,...):
	
	    do j = 2,3
	    	m1haj(lh,j) = 0d0  ! graphs with A exchange
	    	m1hzj(lh,j) = 0d0  !		 Z exchange
            enddo  
		
		
	    do m = 1,12
		i = listab2(m)
		m1haj(lh,2) = m1haj(lh,2) + tampr(i,-1,-1,lh,k) 
	    enddo
		
	    do m = 1,12
		i = listab3(m)
		m1haj(lh,3) = m1haj(lh,3) + tampr(i,-1,-1,lh,k) 
	    enddo
	
	    do m = 1,12
		i = listzb2(m)
		m1hzj(lh,2) = m1hzj(lh,2) + tampr(i,-1,-1,lh,k) 
	    enddo
	
	    do m = 1,12
		i = listzb3(m)
		m1hzj(lh,3) = m1hzj(lh,3) + tampr(i,-1,-1,lh,k) 
	    enddo
	     

	    if (k.eq.5) then  ! for incoming g: no g emission from other line
		      if(.not.jlog2) m1haj(lh,2) = 0d0
		      if(.not.jlog3) m1haj(lh,3) = 0d0
		      if(.not.jlog2) m1hzj(lh,2) = 0d0
		      if(.not.jlog3) m1hzj(lh,3) = 0d0
	    else
		       if(.not.jlog2) m1haj(lh,3) = 0d0
		       if(.not.jlog3) m1haj(lh,2) = 0d0
		       if(.not.jlog2) m1hzj(lh,3) = 0d0
		       if(.not.jlog3) m1hzj(lh,2) = 0d0
	     endif


	  enddo ! end of lh loop 	  
	    
		    
    ! compute m1 in kartesian basis:
	    
	    do l = 1,2
	      do j = 2,3
	    	m1kaj(l,j) = 0d0
	    	m1kzj(l,j) = 0d0
	      enddo
	    enddo  	    
	    
	  do j = 2,3
	    m1kaj(1,j) =     (m1haj(-1,j)-m1haj(1,j))/sqrt(2d0)
	    m1kaj(2,j) = -im*(m1haj(-1,j)+m1haj(1,j))/sqrt(2d0)
	    m1kzj(1,j) =     (m1hzj(-1,j)-m1hzj(1,j))/sqrt(2d0)
	    m1kzj(2,j) = -im*(m1hzj(-1,j)+m1hzj(1,j))/sqrt(2d0)
	  enddo    	
	  

c ---------------------------------		 

       if (.false.) then
c      if (.true.) then

c        do l = 1,1
        do l = 1,2 
	  
	  matot = m1kzj(l,2)+m1kzj(l,3)+m1kaj(l,2)+m1kaj(l,3)
	  ratot = mat(k,-1,-1,2,l,6) + mat(k,-1,-1,3,l,6)
     
	  maz = m1kzj(l,2)+m1kzj(l,3)
	  raz = matz(k,-1,-1,2,l,6) + matz(k,-1,-1,3,l,6)
     
	  maa = m1kaj(l,2)+m1kaj(l,3)
	  raa = mata(k,-1,-1,2,l,6) + mata(k,-1,-1,3,l,6)
     
	  
          print*,"TEST A & Z contributions:"
	    print*,"signs:",sign	
	    print*,"gign:",gsign	    
          print*," k,l = ",k,l
	  print*,"jmin,jmax",jmin,jmax
          print*,"follows: coded mat,mg"
	  print*
 	  print*,"rata: ",raa/maa
    	  print*,"resa: ",raa,maa	  
	  print*
	  print*,"ratz: ",raz/maz
    	  print*,"resz: ",raz,maz
	  print*
	  print*,"rat: ",ratot/matot
    	  print*,"res: ",ratot,matot

 	  print*


 	enddo !l
      endif ! "local" debugging block for A&Z

c ---------------------------------		 

c ---------------------------------
    
      enddo ! k
          
      endif ! debugging


c --------------------------------------------------------------------------

c -----------------------------------------------------------------------
c
c
c and now, finally, the pentagon contributions, i.e. two W"s emitted from the
c  same quark line


      fac = clr(3,3,-1)**2
      
      
c upper line:
      do isig3 = 1,-1,-2
         
	 ! g from upper line, but not in bkjqq:
	 call ket2c(psi(1,-1,1),.true.,p(0,1),-1,p43(0,2),
     #	 	jqq(0,isig3,2,is3,0),bkjqq(1,-1,isig3,1,2,0),bq(0,1))
         call bra2c(psi(1,-1,2),.true.,p(0,2),-1,p43(0,2),
     #	 	jqq(0,isig3,2,is3,0),bkjqq(1,-1,isig3,2,2,0),bq(0,2))


       do l = 1,2
	
	 ! g from lower line:
	 call ket2c(psi(1,-1,1),.true.,p(0,1),-1,p43(0,3),
     #		 jqq(0,isig3,2,is3,l),bkjqq(1,-1,isig3,1,3,l),dummy)
         call bra2c(psi(1,-1,2),.true.,p(0,2),-1,p43(0,3),
     #		 jqq(0,isig3,2,is3,l),bkjqq(1,-1,isig3,2,3,l),dummy)
           
 	 ! g from upper line, in bkjqq:
	 call ket2c(braketg(1,-1,1,l),.false.,pq(0,1),-1,p43(0,2),
     #		 jqq(0,isig3,2,is3,0),gbkjqq(1,-1,isig3,1,2,l),dummy)
         call bra2c(braketg(1,-1,2,l),.false.,pq(0,2),-1,p43(0,2),
     #		 jqq(0,isig3,2,is3,0),gbkjqq(1,-1,isig3,2,2,l),dummy)
 	 
	 call ket2r(bkjqq(1,-1,isig3,1,2,0),.false.,bq(0,1),-1,q,
     #		 eps(0,l),bkjqqg(1,-1,isig3,1,2,l),dummy)
         call bra2r(bkjqq(1,-1,isig3,2,2,0),.false.,bq(0,2),-1,q,
     #		 eps(0,l),bkjqqg(1,-1,isig3,2,2,l),dummy)
          
    
c WW from different fermions (f1/f2):     
c  eps1=wp,eps2=j43,eps3=wm, g from lower line
         m5(2,3,3,l) = -s1c(psiwm(1,2),jqq(0,isig3,2,is3,l),.true.,
     #	 		    -1,psiwp(1,1))
c  eps1=wp,eps2=j43,eps3=wm, g from upper line
         ga = -s1c(braketgWM(1,-1,2,l),jqq(0,isig3,2,is3,0),.true.,
     #	 		    -1,psiwp(1,1))
         gb = -s1c(braketWMg(1,-1,2,l),jqq(0,isig3,2,is3,0),.true.,
     #	 		    -1,psiwp(1,1))
         gc = -s1c(psiwm(1,2),jqq(0,isig3,2,is3,0),.true.,
     #	 		    -1,braketWPg(1,-1,1,l))
         gd = -s1c(psiwm(1,2),jqq(0,isig3,2,is3,0),.true.,
     #	 		    -1,braketgWP(1,-1,1,l))
     
         m5(2,3,2,l) = ga+gb+gc+gd

     
c WW from f2:
c  eps1=j43,eps2=wp,eps3=wm, g from lower line
         m5(1,3,3,l) = -s1c(psiwm(1,2),wp,.true.,-1,
     #					 bkjqq(1,-1,isig3,1,3,l))
c  eps1=j43,eps2=wp,eps3=wm, g from upper line
         ga= -s1c(braketgWM(1,-1,2,l),wp,.true.,-1,
     #					 bkjqq(1,-1,isig3,1,2,0))
         gb= -s1c(braketWMg(1,-1,2,l),wp,.true.,-1,
     #				 	 bkjqq(1,-1,isig3,1,2,0))
	 gc = -s1c(psiwm(1,2),wp,.true.,-1,bkjqqg(1,-1,isig3,1,2,l))
	 gd = -s1c(psiwm(1,2),wp,.true.,-1,gbkjqq(1,-1,isig3,1,2,l))         
	 
	 m5(1,3,2,l) = ga+gb+gc+gd


c WW from f1:
c  eps1=wp,eps2=wm,eps3=j43, g from lower line
         m5(3,3,3,l) = -s1c(bkjqq(1,-1,isig3,2,3,l),wm,
     #					 .true.,-1,psiwp(1,1))
c  eps1=wp,eps2=wm,eps3=j43, g from upper line
         ga = -s1c(gbkjqq(1,-1,isig3,2,2,l),wm,.true.,-1,psiwp(1,1))
         gb = -s1c(bkjqqg(1,-1,isig3,2,2,l),wm,.true.,-1,psiwp(1,1))
         gc = -s1c(bkjqq(1,-1,isig3,2,2,0),wm,
     #					 .true.,-1,braketWPg(1,-1,1,l))
         gd = -s1c(bkjqq(1,-1,isig3,2,2,0),wm,
     #					 .true.,-1,braketgWP(1,-1,1,l))
	 
	 m5(3,3,2,l) = ga+gb+gc+gd
	

c -------------------

    
c WW from different fermions (f1/f2):     
c  eps1=wm,eps2=j43,eps3=wp, g from lower line
         m5(2,4,3,l) = -s1c(psiwp(1,2),jqq(0,isig3,2,is3,l),.true.,
     #	 		    -1,psiwm(1,1))
c  eps1=wm,eps2=j43,eps3=wp, g from upper line
         ga = -s1c(braketgWP(1,-1,2,l),jqq(0,isig3,2,is3,0),.true.,
     #	 		    -1,psiwm(1,1))
         gb = -s1c(braketWPg(1,-1,2,l),jqq(0,isig3,2,is3,0),.true.,
     #	 		    -1,psiwm(1,1))
         gc = -s1c(psiwp(1,2),jqq(0,isig3,2,is3,0),.true.,
     #	 		    -1,braketWMg(1,-1,1,l))
         gd = -s1c(psiwp(1,2),jqq(0,isig3,2,is3,0),.true.,
     #	 		    -1,braketgWM(1,-1,1,l))
     
         m5(2,4,2,l) = ga+gb+gc+gd
	 

c wp, wm from f2     
c  eps1=j43,eps2=wm,eps3=wp, g from lower line
         m5(1,4,3,l) = -s1c(psiwp(1,2),wm,.true.,-1,
     #				 bkjqq(1,-1,isig3,1,3,l))
c  eps1=j43,eps2=wm,eps3=wp, g from upper line
         ga= -s1c(braketgWP(1,-1,2,l),wm,.true.,-1,
     #				 bkjqq(1,-1,isig3,1,2,0))
         gb= -s1c(braketWPg(1,-1,2,l),wm,.true.,-1,
     #				 bkjqq(1,-1,isig3,1,2,0))
	 gc = -s1c(psiwp(1,2),wm,.true.,-1,bkjqqg(1,-1,isig3,1,2,l))
	 gd = -s1c(psiwp(1,2),wm,.true.,-1,gbkjqq(1,-1,isig3,1,2,l))         
	 
	 m5(1,4,2,l) = ga+gb+gc+gd


c wp, wm from f1     
c  eps1=wm,eps2=wp,eps3=j43, g from lower line
         m5(3,4,3,l) = -s1c(bkjqq(1,-1,isig3,2,3,l),wp,
     #					 .true.,-1,psiwm(1,1))
c  eps1=wm,eps2=wp,eps3=j43, g from upper line
         ga = -s1c(gbkjqq(1,-1,isig3,2,2,l),wp,.true.,-1,psiwm(1,1))
         gb = -s1c(bkjqqg(1,-1,isig3,2,2,l),wp,.true.,-1,psiwm(1,1))
         gc = -s1c(bkjqq(1,-1,isig3,2,2,0),wp,.true.,
     #					 -1,braketWMg(1,-1,1,l))
         gd = -s1c(bkjqq(1,-1,isig3,2,2,0),wp,.true.,
     #					 -1,braketgWM(1,-1,1,l))
	 
	 m5(3,4,2,l) = ga+gb+gc+gd

c ----------------	 

         do k = 1,4 ! W+,W- from same or different quark legs
         do j = jmin,jmax
	    kk = mod(k+2,4)
            if (k.eq.2) kk = 4 
            mat(k,-1,isig3,j,l,7) = 
     #	       fac*m5(2,ifl(1,k),j,l)*propt(-1,isig3,kk,j,2)+
     1         fac*(m5(1,ifl(2,k),j,l)+m5(3,ifl(1,k),j,l))*
     2						propt(-1,isig3,k,j,2)
         enddo !j
	 enddo !k
       
      
      do k = 5,6 ! W+,W- from same quark leg
       do j = jmin,jmax
        mat(k,-1,-1,j,l,7) = propt(-1,-1,k,j,2) * fac *
     1                    ( m5(1,ifl(2,k),j,l)+m5(3,ifl(1,k),j,l) )
       enddo !j
      enddo !k


       enddo !l
      enddo !isig3

c -----------
 
c lower line:
      do isig1 = 1,-1,-2
         
	 ! g from lower line, but not in bkjqq:
	 call ket2c(psi(1,-1,3),.true.,p(0,3),-1,p21(0,3),
     #		 jqq(0,isig1,1,is1,0),bkjqq(1,-1,isig1,3,3,0),bq(0,3))
         call bra2c(psi(1,-1,4),.true.,p(0,4),-1,p21(0,3),
     #		 jqq(0,isig1,1,is1,0),bkjqq(1,-1,isig1,4,3,0),bq(0,4))


       do l = 1,2
	
	 ! g from upper line:
	 call ket2c(psi(1,-1,3),.true.,p(0,3),-1,p21(0,2),
     #		 jqq(0,isig1,1,is1,l),bkjqq(1,-1,isig1,3,2,l),dummy)
         call bra2c(psi(1,-1,4),.true.,p(0,4),-1,p21(0,2),
     #		 jqq(0,isig1,1,is1,l),bkjqq(1,-1,isig1,4,2,l),dummy)
           
 	 ! g from lower line, in bkjqq:
	 call ket2c(braketg(1,-1,3,l),.false.,pq(0,3),-1,p21(0,3),
     #		 jqq(0,isig1,1,is1,0),gbkjqq(1,-1,isig1,3,3,l),dummy)
         call bra2c(braketg(1,-1,4,l),.false.,pq(0,4),-1,p21(0,3),
     #		 jqq(0,isig1,1,is1,0),gbkjqq(1,-1,isig1,4,3,l),dummy)
 	 
	 call ket2r(bkjqq(1,-1,isig1,3,3,0),.false.,bq(0,3),-1,q,
     #		 eps(0,l),bkjqqg(1,-1,isig1,3,3,l),dummy)
         call bra2r(bkjqq(1,-1,isig1,4,3,0),.false.,bq(0,4),-1,q,
     #		 eps(0,l),bkjqqg(1,-1,isig1,4,3,l),dummy)


c wp, wm from different fermion line (f3/f4)     
c  eps1=wp,eps2=j21,eps3=wm, g from upper line
         m5(2,3,2,l) = -s1c(psiwm(1,4),jqq(0,isig1,1,is1,l),.true.,
     #	 		    -1,psiwp(1,3))
c  eps1=wp,eps2=j21,eps3=wm, g from lower line
         ga = -s1c(braketgWM(1,-1,4,l),jqq(0,isig1,1,is1,0),.true.,
     #	 		    -1,psiwp(1,3))
         gb = -s1c(braketWMg(1,-1,4,l),jqq(0,isig1,1,is1,0),.true.,
     #	 		    -1,psiwp(1,3))
         gc = -s1c(psiwm(1,4),jqq(0,isig1,1,is1,0),.true.,
     #	 		    -1,braketWPg(1,-1,3,l))
         gd = -s1c(psiwm(1,4),jqq(0,isig1,1,is1,0),.true.,
     #	 		    -1,braketgWP(1,-1,3,l))
     
         m5(2,3,3,l) = ga+gb+gc+gd

     
c wp, wm from f4     
c  eps1=j21,eps2=wp,eps3=wm, g from upper line
         m5(1,3,2,l) = -s1c(psiwm(1,4),wp,.true.,-1,
     #					 bkjqq(1,-1,isig1,3,2,l))
c  eps1=j21,eps2=wp,eps3=wm, g from lower line
         ga= -s1c(braketgWM(1,-1,4,l),wp,.true.,-1,
     #					 bkjqq(1,-1,isig1,3,3,0))
         gb= -s1c(braketWMg(1,-1,4,l),wp,.true.,-1,
     #					 bkjqq(1,-1,isig1,3,3,0))
	 gc = -s1c(psiwm(1,4),wp,.true.,-1,bkjqqg(1,-1,isig1,3,3,l))
	 gd = -s1c(psiwm(1,4),wp,.true.,-1,gbkjqq(1,-1,isig1,3,3,l))         
	 
	 m5(1,3,3,l) = ga+gb+gc+gd


c wp, wm from f3     
c  eps1=wp,eps2=wm,eps3=j21, g from upper line
         m5(3,3,2,l) = -s1c(bkjqq(1,-1,isig1,4,2,l),wm,
     #					 .true.,-1,psiwp(1,3))
c  eps1=wp,eps2=wm,eps3=j21, g from lower line
         ga = -s1c(gbkjqq(1,-1,isig1,4,3,l),wm,.true.,-1,psiwp(1,3))
         gb = -s1c(bkjqqg(1,-1,isig1,4,3,l),wm,.true.,-1,psiwp(1,3))
         gc = -s1c(bkjqq(1,-1,isig1,4,3,0),wm,.true.,-1,
     #					 braketWPg(1,-1,3,l))
         gd = -s1c(bkjqq(1,-1,isig1,4,3,0),wm,.true.,-1,
     #			  		 braketgWP(1,-1,3,l))
	 
	 m5(3,3,3,l) = ga+gb+gc+gd
	 
c -------------------
     
     
c WW from different fermions (f3/f4):     
c  eps1=wm,eps2=j21,eps3=wp, g from upper line
         m5(2,4,2,l) = -s1c(psiwp(1,4),jqq(0,isig1,1,is1,l),.true.,
     #	 		    -1,psiwm(1,3))
c  eps1=wm,eps2=j21,eps3=wp, g from lower line
         ga = -s1c(braketgWP(1,-1,4,l),jqq(0,isig1,1,is1,0),.true.,
     #	 		    -1,psiwm(1,3))
         gb = -s1c(braketWPg(1,-1,4,l),jqq(0,isig1,1,is1,0),.true.,
     #	 		    -1,psiwm(1,3))
         gc = -s1c(psiwp(1,4),jqq(0,isig1,1,is1,0),.true.,
     #	 		    -1,braketWMg(1,-1,3,l))
         gd = -s1c(psiwp(1,4),jqq(0,isig1,1,is1,0),.true.,
     #	 		    -1,braketgWM(1,-1,3,l))
     
         m5(2,4,3,l) = ga+gb+gc+gd


c wp, wm from f4     
c  eps1=j21,eps2=wp,eps3=wm, g from upper line
         m5(1,4,2,l) = -s1c(psiwp(1,4),wm,.true.,-1,
     #					 bkjqq(1,-1,isig1,3,2,l))
c  eps1=j21,eps2=wp,eps3=wm, g from lower line
         ga= -s1c(braketgWP(1,-1,4,l),wm,.true.,-1,
     #					 bkjqq(1,-1,isig1,3,3,0))
         gb= -s1c(braketWPg(1,-1,4,l),wm,.true.,-1,
     #					 bkjqq(1,-1,isig1,3,3,0))
	 gc = -s1c(psiwp(1,4),wm,.true.,-1,bkjqqg(1,-1,isig1,3,3,l))
	 gd = -s1c(psiwp(1,4),wm,.true.,-1,gbkjqq(1,-1,isig1,3,3,l))         
	 
	 m5(1,4,3,l) = ga+gb+gc+gd


c wp, wm from f3     
c  eps1=wp,eps2=wm,eps3=j21, g from upper line
         m5(3,4,2,l) = -s1c(bkjqq(1,-1,isig1,4,2,l),wp,
     #					 .true.,-1,psiwm(1,3))
c  eps1=wp,eps2=wm,eps3=j21, g from lower line
         ga = -s1c(gbkjqq(1,-1,isig1,4,3,l),wp,.true.,-1,psiwm(1,3))
         gb = -s1c(bkjqqg(1,-1,isig1,4,3,l),wp,.true.,-1,psiwm(1,3))
         gc = -s1c(bkjqq(1,-1,isig1,4,3,0),wp,.true.,-1,
     #					 braketWMg(1,-1,3,l))
         gd = -s1c(bkjqq(1,-1,isig1,4,3,0),wp,.true.,-1,
     #					 braketgWM(1,-1,3,l))
	 
	 m5(3,4,3,l) = ga+gb+gc+gd
	

c ----------------	 
 
         do k = 1,4
         do j = jmin,jmax
            kk = k+1  - 2*mod(k+1,2)
            mat(k,isig1,-1,j,l,8) = 
     #	       fac*m5(2,ifl(3,k),j,l)*propt(isig1,-1,kk,j,1)+
     1         fac*( m5(1,ifl(4,k),j,l)+m5(3,ifl(3,k),j,l) )*
     2						propt(isig1,-1,k,j,1)
         enddo !j
         enddo
 
      do k = 5,6
         do j = jmin,jmax
         mat(k,-1,-1,j,l,8) = propt(-1,-1,k,j,1) * fac *
     1                    ( m5(1,ifl(4,k),j,l)+m5(3,ifl(3,k),j,l) )
         enddo !j
         
      enddo !k
       
	    
      enddo !l     
      
      enddo !isig

c 
c -----------------------------------------------------------------------

c ----------------------------------------------------------------    
c
c debugging:
c
c check result for NC as obtained above against madgraph output 
c		(sum up relevant mg amplitudes in the following):

      if (ldebugm) then

       do k = 1,2
c       do k = 1,1
	
	  do lh = -1,1,2	! gluon pol. in hel. basis
 	    
	! need to distingiush between j=2,3
		
c compute reference result from individual madgraph diagrams tampr(i,...):
	
	    do j = 2,3
	    	m1hua(lh,j) = 0d0  ! graphs with p from upper line
	    	m1hla(lh,j) = 0d0  !	-""- 	       lower
 	    	m1huz(lh,j) = 0d0  ! graphs with p from upper line
	    	m1hlz(lh,j) = 0d0  !	-""- 	       lower
		
		m1hc(lh) = 0d0
            enddo  
		
		
	    do m = 1,12
		i = lncapu2(m)
		m1hua(lh,2) = m1hua(lh,2) + tampr(i,-1,-1,lh,k) 
	    enddo		
	    do m = 1,12
		i = lnczpu2(m)
		m1huz(lh,2) = m1huz(lh,2) + tampr(i,-1,-1,lh,k) 
	    enddo	
	    do m = 1,6
		i = lncapu3(m)
		m1hua(lh,3) = m1hua(lh,3) + tampr(i,-1,-1,lh,k) 
	    enddo	
	    do m = 1,6
		i = lnczpu3(m)
		m1huz(lh,3) = m1huz(lh,3) + tampr(i,-1,-1,lh,k) 
	    enddo		
	    do m = 1,6
		i = lncapl2(m)
		m1hla(lh,2) = m1hla(lh,2) + tampr(i,-1,-1,lh,k) 
	    enddo		
	    do m = 1,6
		i = lnczpl2(m)
		m1hlz(lh,2) = m1hlz(lh,2) + tampr(i,-1,-1,lh,k) 
	    enddo	
	    do m = 1,12
		i = lncapl3(m)
		m1hla(lh,3) = m1hla(lh,3) + tampr(i,-1,-1,lh,k) 
	    enddo	
	    do m = 1,12
		i = lnczpl3(m)
		m1hlz(lh,3) = m1hlz(lh,3) + tampr(i,-1,-1,lh,k) 
	    enddo

	    do  i = 7,18
		m1hc(lh) = m1hc(lh) + tampr(i,-1,-1,lh,k) 
	    enddo
	    do  i =  74,85
		m1hc(lh) = m1hc(lh) + tampr(i,-1,-1,lh,k) 
	    enddo
	    do  i =  136,147
		m1hc(lh) = m1hc(lh) + tampr(i,-1,-1,lh,k) 
	    enddo
	    do  i =  228,239
		m1hc(lh) = m1hc(lh) + tampr(i,-1,-1,lh,k) 
	    enddo
	    do  i =  320,331
		m1hc(lh) = m1hc(lh) + tampr(i,-1,-1,lh,k) 
	    enddo
	    do  i =  393,404
		m1hc(lh) = m1hc(lh) + tampr(i,-1,-1,lh,k) 
	    enddo


	   if (k.ne.3) then  ! for incoming g: no g emission from other line
		     if(.not.jlog2) m1hua(lh,2) = 0d0
		     if(.not.jlog3) m1hua(lh,3) = 0d0
		     if(.not.jlog2) m1hla(lh,2) = 0d0
		     if(.not.jlog3) m1hla(lh,3) = 0d0
		     if(.not.jlog2) m1huz(lh,2) = 0d0
		     if(.not.jlog3) m1huz(lh,3) = 0d0
		     if(.not.jlog2) m1hlz(lh,2) = 0d0
		     if(.not.jlog3) m1hlz(lh,3) = 0d0
	   else
		      if(.not.jlog2) m1hua(lh,3) = 0d0
		      if(.not.jlog3) m1hua(lh,2) = 0d0
		      if(.not.jlog2) m1hla(lh,3) = 0d0
		      if(.not.jlog3) m1hla(lh,2) = 0d0
		      if(.not.jlog2) m1huz(lh,3) = 0d0
		      if(.not.jlog3) m1huz(lh,2) = 0d0
		      if(.not.jlog2) m1hlz(lh,3) = 0d0
		      if(.not.jlog3) m1hlz(lh,2) = 0d0
	    endif


	  enddo ! end of lh loop 	  
	    
		    
    ! compute m1 in kartesian basis:
	    
	    do l = 1,2
	      do j = 2,3
	    	m1kua(l,j) = 0d0
	    	m1kla(l,j) = 0d0
	    	m1kuz(l,j) = 0d0
	    	m1klz(l,j) = 0d0
	      enddo
	      m1kc(l) = 0d0	      
	    enddo  	    
	    
	  do j = 2,3
	    m1kua(1,j) =     (m1hua(-1,j)-m1hua(1,j))/sqrt(2d0)
	    m1kua(2,j) = -im*(m1hua(-1,j)+m1hua(1,j))/sqrt(2d0)
	    m1kla(1,j) =     (m1hla(-1,j)-m1hla(1,j))/sqrt(2d0)
	    m1kla(2,j) = -im*(m1hla(-1,j)+m1hla(1,j))/sqrt(2d0)	    
	    m1kuz(1,j) =     (m1huz(-1,j)-m1huz(1,j))/sqrt(2d0)
	    m1kuz(2,j) = -im*(m1huz(-1,j)+m1huz(1,j))/sqrt(2d0)
	    m1klz(1,j) =     (m1hlz(-1,j)-m1hlz(1,j))/sqrt(2d0)
	    m1klz(2,j) = -im*(m1hlz(-1,j)+m1hlz(1,j))/sqrt(2d0)	    
	  enddo    	

	  m1kc(1) =	(m1hc(-1)-m1hc(1))/sqrt(2d0)
	  m1kc(2) = -im*(m1hc(-1)+m1hc(1))/sqrt(2d0)
	  
c ---------------------------------		 

c       if (.false.) then
      if (.true.) then

c        do l = 1,1
        do l = 1,2 
	  
	 if(k.eq.1) then
	 
	  matot = m1kua(l,2)+m1kua(l,3)+m1kla(l,2)+m1kla(l,3)+
     # 		  m1kuz(l,2)+m1kuz(l,3)+m1klz(l,2)+m1klz(l,3)
	  ratot = mat(k,-1,-1,2,l,7) + mat(k,-1,-1,2,l,8) +
     #		  mat(k,-1,-1,3,l,7) + mat(k,-1,-1,3,l,8)
     
	 
	 elseif(k.eq.2) then
	 
	  matot = m1kc(l)
	  ratot = mat(k,-1,-1,2,l,7) + mat(k,-1,-1,2,l,8) +
     #		  mat(k,-1,-1,3,l,7) + mat(k,-1,-1,3,l,8)
     	  
	  endif
         
	  
          print*,"TEST WW contributions for pentagons:"
	    print*,"signs:",sign	
	    print*,"gign:",gsign	    
          print*," k,l = ",k,l
	  print*,"jmin,jmax",jmin,jmax
          print*,"follows: coded mat,mg"
	  print*
 	  print*,"rat: ",ratot/matot
	  print*,"res: ",ratot,matot
 	  print*


 	enddo !l
      endif ! "local" debugging block for pentagon

c ---------------------------------
    
      enddo ! k
      
c -------------------------------------------------------------------------
c
c CHARGED CURRENT:	

       do k = 5,6
c       do k = 5,5
	
	  do lh = -1,1,2	! gluon pol. in hel. basis
 	    
	! need to distingiush between j=2,3
		
c compute reference result from individual madgraph diagrams tampr(i,...):
	
	    do j = 2,3
	    	m1hu(lh,j) = 0d0  ! graphs with p from upper line
	    	m1hl(lh,j) = 0d0  !	-""- 	       lower
            enddo  
		
		
	    do m = 1,8
		i = listpu2(m)
		m1hu(lh,2) = m1hu(lh,2) + tampr(i,-1,-1,lh,k) 
	    enddo		
	    do m = 1,4
		i = listpu3(m)
		m1hu(lh,3) = m1hu(lh,3) + tampr(i,-1,-1,lh,k) 
	    enddo	
	    do m = 1,4
		i = listpl2(m)
		m1hl(lh,2) = m1hl(lh,2) + tampr(i,-1,-1,lh,k) 
	    enddo	
	    do m = 1,8
		i = listpl3(m)
		m1hl(lh,3) = m1hl(lh,3) + tampr(i,-1,-1,lh,k) 
	    enddo


	   if (k.eq.5) then  ! for incoming g: no g emission from other line
		     if(.not.jlog2) m1hu(lh,2) = 0d0
		     if(.not.jlog3) m1hu(lh,3) = 0d0
		     if(.not.jlog2) m1hl(lh,2) = 0d0
		     if(.not.jlog3) m1hl(lh,3) = 0d0
	   else
		      if(.not.jlog2) m1hu(lh,3) = 0d0
		      if(.not.jlog3) m1hu(lh,2) = 0d0
		      if(.not.jlog2) m1hl(lh,3) = 0d0
		      if(.not.jlog3) m1hl(lh,2) = 0d0
	    endif


	  enddo ! end of lh loop 	  
	    
		    
    ! compute m1 in kartesian basis:
	    
	    do l = 1,2
	      do j = 2,3
	    	m1ku(l,j) = 0d0
	    	m1kl(l,j) = 0d0
	      enddo
	    enddo  	    
	    
	  do j = 2,3
	    m1ku(1,j) =     (m1hu(-1,j)-m1hu(1,j))/sqrt(2d0)
	    m1ku(2,j) = -im*(m1hu(-1,j)+m1hu(1,j))/sqrt(2d0)
	    m1kl(1,j) =     (m1hl(-1,j)-m1hl(1,j))/sqrt(2d0)
	    m1kl(2,j) = -im*(m1hl(-1,j)+m1hl(1,j))/sqrt(2d0)	    
	  enddo    	
	  
c ---------------------------------		 

       if (.false.) then
c      if (.true.) then

        do l = 1,1
c        do l = 1,2 
	  
	  matot = m1ku(l,2)+m1ku(l,3)+m1kl(l,2)+m1kl(l,3)
	  ratot = mat(k,-1,-1,2,l,7) + mat(k,-1,-1,2,l,8) +
     #		  mat(k,-1,-1,3,l,7) + mat(k,-1,-1,3,l,8)
          
	  
          print*,"TEST WW contributions:"
	    print*,"signs:",sign	
	    print*,"gign:",gsign	    
          print*," k,l = ",k,l
	  print*,"jmin,jmax",jmin,jmax
          print*,"follows: coded mat,mg"
	  print*
  	  print*,"rat: ",ratot/matot
	  print*,"res: ",ratot,matot
 	  print*


 	enddo !l
      endif ! "local" debugging block for pentagon

c ---------------------------------
    
      enddo ! k
          
      endif ! debugging


c --------------------------------------------------------------------------
c --------------------------------------------------------------------------

c sum the graphs, square them and map them onto uucc, uuss etc.
      do k = 1,4
 	 do j = 2,3     
 	    res(k,j) = 0
	    do isig1 = -1,1,2
	       do isig3 = -1,1,2
 	    	  do l = 1,2
              	     mm(k,isig1,isig3,j,l) = 0
               	     do i = 1,8
                        mm(k,isig1,isig3,j,l) = 
     1                 		 mm(k,isig1,isig3,j,l) + 
     1		    	       (mat(k,isig1,isig3,j,l,i))
      		     enddo !i
              	     res(k,j) = res(k,j) 
     &		       	       + dreal(mm(k,isig1,isig3,j,l))**2
     &                         + dimag(mm(k,isig1,isig3,j,l))**2
	          enddo !l
	      enddo !isig3		     
           enddo !isig1
           res(k,j) = res(k,j)*12d0*fpials(j)   ! C_2*9 is the color factor
	 enddo !j
      enddo !k      
      
      do k = 5,6
 	 do j = 2,3     
 	    res(k,j) = 0
 	    do l = 1,2
               mm(k,-1,-1,j,l) = 0
               do i = 1,8
            	  mm(k,-1,-1,j,l) = 
     1      		   mm(k,-1,-1,j,l) + 
     1	    		 (mat(k,-1,-1,j,l,i))
     	       enddo !i
               res(k,j) = res(k,j) 
     &	    		 + dreal(mm(k,-1,-1,j,l))**2
     &      		 + dimag(mm(k,-1,-1,j,l))**2
	    enddo !l
            res(k,j) = res(k,j)*12d0*fpials(j)   ! C_2*9 is the color factor
 	 enddo !j
      enddo !k      
      
      if (jmin.eq.3) then
      	uucc(2) = 0d0
      	uuss(2) = 0d0
      	ddcc(2) = 0d0
      	ddss(2) = 0d0
      	udsc(2) = 0d0
      	ducs(2) = 0d0
      elseif (jmax.eq.2) then
      	uucc(3) = 0d0
      	uuss(3) = 0d0
      	ddcc(3) = 0d0
      	ddss(3) = 0d0
      	udsc(3) = 0d0
      	ducs(3) = 0d0
      endif	
		    
      if (ldebug) then
        do j = jmin,jmax
	if(udsc(j).ne.0d0) then
	 print*,"gsign:",gsign
	 print*,"jmin,jmax=",jmin,jmax
	 print*,"j=",j
         print*,"uucc",uucc(j),(res(1,j))/uucc(j)
         print*,"uuss",uuss(j),(res(2,j))/uuss(j)
         print*,"ddcc",ddcc(j),(res(3,j))/ddcc(j)
         print*,"ddss",ddss(j),(res(4,j))/ddss(j)
         print*,"udsc",udsc(j),(res(5,j))/udsc(j)
         print*,"ducs",ducs(j),(res(6,j))/ducs(j)
c         read(*,"(a)") text
	endif
	enddo 
 	 print*,"sum j = 2,3:"
         print*,"uucc",uucc(1),(res(1,2)+res(1,3))/uucc(1)
c         print*,"    ",res(1,2),res(1,3)
         print*,"uuss",uuss(1),(res(2,2)+res(2,3))/uuss(1)
         print*,"ddcc",ddcc(1),(res(3,2)+res(3,3))/ddcc(1)
         print*,"ddss",ddss(1),(res(4,2)+res(4,3))/ddss(1)
         print*,"udsc",udsc(1),(res(5,2)+res(5,3))/udsc(1)
         print*,"ducs",ducs(1),(res(6,2)+res(6,3))/ducs(1)
      endif

c      if (lwarn) then
      if (.false.) then
         resv(1) = uucc(1)
         resv(2) = uuss(1)
         resv(3) = ddcc(1)
         resv(4) = ddss(1)
         resv(5) = udsc(1)
         resv(6) = ducs(1)
         do k = 1,6
            if (((res(k,2)+res(k,3))/resv(k)-1).gt.1d-6) then
               if (.not.ldebug) then
                  print*," k =",k,":",(res(k,2)+res(k,3))/resv(k)-1
                  ldebug = .true.
                  lerror = .true.
                  goto 701
               endif
            endif
         enddo !k
         if (lerror) then
            print*," "
            lerror = .false.
            ldebug = .false.
	    print*,"set ldebug false"
         endif
      endif ! lwarn

c  set NC processes to zero if photon virtuality falls below cutoff
      if ( abs(p21(4,3)).lt.qsqAmin .or. abs(p43(4,3)).lt.qsqAmin ) then
         do k = 1,4
            res(k,3) = 0      ! no graphs with gluon on lower line
         enddo
      endif
      if ( abs(p43(4,2)).lt.qsqAmin .or. abs(p21(4,2)).lt.qsqAmin ) then
         do k = 1,4
            res(k,2) = 0      ! no graphs with gluon on upper line
         enddo
      endif

      if (ldoblha) then
        if (blha_amptype(blha_curproc).eq.0) then ! tree
          if ( (blha_idsubproc(blha_cursubproc,blha_curproc)/10 .eq. 1) .or.
     &         (blha_idsubproc(blha_cursubproc,blha_curproc)/10 .eq. 2) ) then
            blha_amp(1)=blha_amp(1)+blha_CF/(4d0/3d0)*
     &       res(mod(blha_idsubproc(blha_cursubproc,blha_curproc),10),2)
          endif
          if ( (blha_idsubproc(blha_cursubproc,blha_curproc)/10 .eq. 1) .or.
     &         (blha_idsubproc(blha_cursubproc,blha_curproc)/10 .eq. 3) ) then
            blha_amp(1)=blha_amp(1)+blha_CF/(4d0/3d0)*
     &       res(mod(blha_idsubproc(blha_cursubproc,blha_curproc),10),3)
          endif

        else 
          call BLHA_amptypeerror(blha_amptype(blha_curproc),__FILE__,__LINE__)
        endif
      endif

      do j = 2,3
         uucc(j) = res(1,j) 
         uuss(j) = res(2,j)
         ddcc(j) = res(3,j)
         ddss(j) = res(4,j)
         udsc(j) = res(5,j)
         ducs(j) = res(6,j)
      enddo
      uucc(1) = res(1,3) +  res(1,2)
      uuss(1) = res(2,3) +  res(2,2)
      ddcc(1) = res(3,3) +  res(3,2)
      ddss(1) = res(4,3) +  res(4,2)
      udsc(1) = res(5,3) +  res(5,2)
      ducs(1) = res(6,3) +  res(6,2)

      return
      end

