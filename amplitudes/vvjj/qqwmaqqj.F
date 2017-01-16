c---------------------------------------------------------------------------
c     lepton helicities have to be summed over externally 
c     ( this routine is for one fixed set of lep. helicities only )
c 
c     start to implement real emission in LO code here
c
c---------------------------------------------------------------------------
      
      subroutine qqwmaqqj(pbar,sign,qbar,gsign,bos,
     1                  ans1,ans2,ans3,ans4)
          use globalvars, only: ldoblha
      implicit none
c
c     Nicolas Kaiser <Nicolas.Kaiser@student.kit.edu>  
c     based on qqwmzqqj by Barbara Jaeger
c     Initial version:  July 2013
c     Last modified:     
C
C     qqwmaqqj calculates the matrix elements**2 for electroweak
c     weak boson pair production by quark quark scattering
C
C        q1 q3    ---->   q2 q4 g W- A
C
c     with W-A -> e+veA  (bos = 43) 	
C     Crossing related processes are also computed. Pauli interference terms for
c     identical fermions are neglected. In particular, only the
c     t-channel exchange of elctroweak bosons is considered. s-channel
c     production of 3 weak bosons is NOT implemented.
c     Gluon polarization vectors are given in the kartesian basis 
c     (comparison with madgraph amplitudes tamprw requires trafo to 
c     kart. basis also). 
c
C     This code includes only real emission contributions, i.e.
c     
c     return uucc = |M_real|^2   etc.
c
c     fpials is attached only in the end of the code
c
c     index j = 2:3 indicates, whether g is emitted from 
c     upper (12) line (j=2) or lower (34) line (j=3)
c
c     l is the gluon polarization in the kartesian basis (l=1,2)
c     l=0 stands for building blocks without gluon emission
c     k is the process ID (1:uucc,2:uuss,3:ddcc;4:ddss;5:udsc;6:ducs)
c     isig1/isig3 are the helicities of parton1/3 
c     
c
c---------------------------------------------------------------------
c
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"
#include "tensorw.inc"
c
c electroweak couplings are taken from KOPPLN
c
      double precision  clr, xm2, xmg, b, v, a
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     1                  V(4,5),A(4,5)
c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"
c
c variables for the main part of the program
c
c
      real*8 fpials(2:3), fpi
      parameter (fpi=4d0*pi)

      double precision  pbar(0:3,4+max_v),qbar(0:4)
      double precision ans1(3), ans2(3), ans3(3), ans4(3),res(4,2:3)

      double precision  p(0:3,4+max_v),q(0:4), p21(0:4,2:3), p43(0:4,2:3),
     1                  pq(0:4,4),pew(0:4),pa(0:4),pwa(0:4)
      integer  sign(4+max_v),gsign, mu, i, j, k, kk,is,
     1         isig, isig1, isig3,bos,il,arg,arg2,jout
      integer  ifl(4,6), is1, is3
      integer  l   ! gluon polariz. (l=0:no g, l=1,2:g in kartesian basis)
      integer jmin, jmax
      logical jlog2,jlog3
      double complex prop21(4,2:3), prop43(4,2:3),prop_ee(3,2:3),prop_uu(3,2:3)
      double complex mat(4,-1:1,-1:1,2:3,0:2,9)
      double complex mm(4,-1:1,-1:1,2:3,2)
      double complex maw,mzw,mwa,mwz

      double complex mw(2,2:3,2)
      double complex mbwa(2,3:4,-1:1,-1:1,2:3,2),mbaw(2,3:4,-1:1,-1:1,2:3,2),
     1		     zwa(-1:1,-1:1,3,2:3,2),zaw(-1:1,-1:1,3,2:3,2),
     1		     mbwv,mbvw,gbwv,gbvw
      double complex m5(3,3,3,2:3,2)
      double complex ga,gb,gc,gd
      double complex 
     1		     m1we,m2we,m3we,m4we,w1we,w2we,w3we,w4we,
     2		     m1e,m2e,m3e,m4e,
     3	             m1aau,m1azu,z1azu,
     4		     mw1aau,mw1azu,
     3	             m2aau,m2azu,
     4		     mw2aau,mw2azu,
     3	             m3aau,m3azu,z3azu,
     4		     mw3aau,mw3azu,
     3	             m4aau,m4azu,
     4		     mw4aau,mw4azu,
     2		     m1u,m2u,m3u,m4u,
     2		     m1awu,z1awu,m1wau,
     2		     m2awu,m2wau,
     2		     m3awu,z3awu,m3wau,
     2		     m4awu,m4wau

      double precision eps(0:3,2) ! g in kartesian basis 
      double complex psi(2,-1:1,4), jqq(0:5,-1:1,2,-1:1,0:2), 
     1 		     braketg(2,-1:1,4,2), jh1(0:5,-1:1), jh2(0:5,-1:1)
      double complex braketgw(2,-1:1,0:2,4),braketwg(2,-1:1,0:2,4)      
      double precision  pgw(0:4,4),pwg(0:4,4)
      double complex psiwe(2,-1:1,4),psia(2,-1:1,4),
     1               jew(0:5,-1:1,-1:1,4,0:2),jua(0:5,-1:1,-1:1,4,0:2)
      double complex bkjqq(2,-1:1,-1:1,4,2:3,0:2),
     1		     bkjqqg(2,-1:1,-1:1,4,2:3,0:2),
     1 		     gbkjqq(2,-1:1,-1:1,4,2:3,0:2)
      double complex jw(0:5,-1:1,-1:1,4,0:2),
     1		     jwg(0:5,-1:1),jgw(0:5,-1:1),jg0(0:5,-1:1)
      double complex psiw(2,-1:1,-1:1,4)
      double complex braketgwe(2,-1:1,4,0:2),braketweg(2,-1:1,4,0:2),
     3		     braketgau(2,-1:1,4,0:2),braketaug(2,-1:1,4,0:2)      
      double precision  pgwe(0:4,4),pweg(0:4,4),pgau(0:4,4),paug(0:4,4)
      double complex jwegi(0:5,-1:1),jgwei(0:5,-1:1),jweg0i(0:5,-1:1),
     1		     jwegii(0:5,-1:1),jgweii(0:5,-1:1),jweg0ii(0:5,-1:1),
     1		     jaugi(0:5,-1:1),jgaui(0:5,-1:1),jaug0i(0:5,-1:1),
     1		     jaugii(0:5,-1:1),jgauii(0:5,-1:1),jaug0ii(0:5,-1:1)

      double complex epsww(0:3),
     1 		     epsaw(0:3),epszw(0:3),epswa(0:3),epswz(0:3),
     1		     epsawe(0:3,-1:1,3:4,2,2:3,0:2), epszwe(0:3,-1:1,3:4,2,2:3,0:2),
     1		     epswae(0:3,-1:1,3:4,2,2:3,0:2), epswze(0:3,-1:1,3:4,2,2:3,0:2),
     1		     epswee(0:3,-1:1,3:4,2,2:3,0:2), epswuu(0:3,-1:1,3:4,2,2:3,0:2),
     1		     qepswee(3:4,-1:1,2:3,0:2),qepszwe(3:4,-1:1,2:3,0:2),
     1		     qepswze(3:4,-1:1,2:3,0:2),
     1		     qepswuu(3:4,-1:1,2:3,0:2)	 
      double precision fqwe(0:4,4),fqau(0:4,4)
      double precision fq(0:4,4),qee(0:4,2:3), qaa(0:4,2:3),dummy(0:4),bq(0:4,4)
      double complex zm2i(2:4)
      double complex jj21we(-1:1,0:2),jj43we(-1:1,0:2),
     1 		     jj21au(-1:1,0:2),jj43au(-1:1,0:2)

      double complex propt1(-1:1,-1:1,4,2,2:3), propt2(-1:1,-1:1,4,2,2:3),
     1		     propw(2,2:3)
c
       double complex psiwg,jwgg(2,-1:1,-1:1,3,2:3,2)
       integer kb,kp
c
      double complex im
      parameter (im = (0d0,1d0))

      integer idu,idl
      double complex contract_Tjj, dotcc, dotrc, dotqj, s1c
      external contract_Tjj, dotcc, dotrc, dotqj, s1c
      logical ldebug, linit,lerror,lmgdebug
      
      data linit /.true./, lerror /.false./, ldebug /.false./
      parameter (lmgdebug = .false.)
     
      integer*8 icb1, icount1, icb2, icount2
      data icb1/0/,icount1/0/,icb2/0/,icount2/0/

      save ifl, zm2i, linit, ldebug
      save jw,psiw,fq
      double complex  zero
      parameter (zero = (0d0,0d0) )
      integer ll

      if (bos.ne.41) then
         print*, "subroutine qqwmzqqj contains only",
     $        " ME for W-Ajjj production"
         print*, "but you called it with bos = ",bos
         stop
      endif	
      
c-------------------------
c     
c     
c     initialize & precompute stuff needed below:
c     
c     fix strong coupling gs**2 for the 2 quark lines:

      fpials(2) = fpi*als(1,1)
      fpials(3) = fpi*als(2,1)
c     
c     define flavors of external quarks for the 4 subprocesses
c     
      if (linit) then
         linit = .false.
         k = 1                  ! uusc
         ifl(1,k) = 3
         ifl(2,k) = 3
         ifl(3,k) = 4
         ifl(4,k) = 3
         k = 2                  ! ddsc
         ifl(1,k) = 4
         ifl(2,k) = 4
         ifl(3,k) = 4
         ifl(4,k) = 3
         k = 3                  ! ducc
         ifl(1,k) = 4
         ifl(2,k) = 3
         ifl(3,k) = 3
         ifl(4,k) = 3
         k = 4                  ! duss
         ifl(1,k) = 4
         ifl(2,k) = 3
         ifl(3,k) = 4
         ifl(4,k) = 4
c     
         zm2i(2) = 1/dcmplx(xm2(2),-xmg(2))
         zm2i(3) = 1/dcmplx(xm2(3),-xmg(3))
         zm2i(4) = 1/dcmplx(xm2(4),-xmg(4))

      endif

      if (ldebug) then
         print*," " 
         print*," enter qqwmaqqj"
      endif
      
      if (gsign.eq.1) then	!final state gluon
         jlog2 = .true.		! can couple to upper/lower line
         jlog3 = .true. 
         jmin = 2
         jmax = 3
      else                      !initial state gluon -> only:
	 if (sign(1).ne.sign(2)) then !gluon from upper line
            jlog2 = .true.
            jlog3 = .false. 
            jmin = 2
            jmax = 2
	 else                   !gluon from lower line
            jlog2 = .false.
            jlog3 = .true. 
            jmin = 3
            jmax = 3
         endif
      endif	

      do k = 1,4
         do isig1 = -1,1,2
            do isig3 = -1,1,2
	       do j = 2,3
	          do l = 0,2
               	     do i = 1,9
                  	mat(k,isig1,isig3,j,l,i) = 0
    		     enddo
		  enddo	
               enddo
            enddo
         enddo
      enddo
c     
c     identify fermion line sign factors (for 1 3 -> 2 4 etc.)
c     
c     fix is1 such that is1 = +1 for   q1 ->  q2 g 
c     is1 = -1 for  q1b -> q2b g 
c     is1 =  0 for    g -> q1b q2
c     
c     for is3:     is3 = +1 for   q3 ->  q4 g 
c     is3 = -1 for  q3b -> q4b g 
c     is3 =  0 for    g -> q3b q4

      is1 = (sign(1)+sign(2))/2  
      is3 = (sign(3)+sign(4))/2
      
c     (is1,is3 are fixed here and don"t change throughout this run of the program)

c     
c     define the internal momenta
c     
      do mu = 0,3
         do i = 1,4+n_v
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo
         q(mu) = qbar(mu)*gsign
         
         p21(mu,3) = p(mu,2) - p(mu,1) ! g from 34 line
         p21(mu,2) = p21(mu,3) + q(mu) ! g from 12 line
         p43(mu,2) = p(mu,4) - p(mu,3)
         p43(mu,3) = p43(mu,2) + q(mu)	 
         
         pew(mu) = -(p(mu,6) - p(mu,5))
         pa(mu) =   p(mu,7)
         pwa(mu) =  pew(mu) + pa(mu)  	       
         
      enddo                     !mu
      
      do j = 2,3
         p21(4,j) = p21(0,j)**2 - p21(1,j)**2 - p21(2,j)**2 - p21(3,j)**2
         p43(4,j) = p43(0,j)**2 - p43(1,j)**2 - p43(2,j)**2 - p43(3,j)**2
      enddo
      
      q(4) = 0d0
      pew(4) = pew(0)**2 - pew(1)**2 - pew(2)**2 - pew(3)**2
      pa(4)  = 0d0
      pwa(4) = pwa(0)**2 - pwa(1)**2 - pwa(2)**2 - pwa(3)**2
      

c     ---------------------------------------------------------------
c     
c     get the vector boson propagator factors
c     
c     depending on value of j, gluon is attached to respective quark line or not;
c     no V is attached here 
c     
      do j = 2,3	
         prop21(1,j) = 1/p21(4,j)
         prop21(2,j) = 1/dcmplx(p21(4,j)-xm2(2),xmg(2))
         prop21(3,j) = 1/dcmplx(p21(4,j)-xm2(3),xmg(3))
         prop21(4,j) = prop21(3,j)

         prop43(1,j) = 1/p43(4,j)
         prop43(2,j) = 1/dcmplx(p43(4,j)-xm2(2),xmg(2))
         prop43(3,j) = 1/dcmplx(p43(4,j)-xm2(3),xmg(3))
         prop43(4,j) = prop43(3,j)
      enddo                     !j

c     
c     ---------------------------------------------------------------------
c     
c     for box-box and BV graphs we need the propagators for t-channel bosons between quark lines
c     as seen from upper line these W momenta are INCOMING. They are OUTGOING as seen from lower 
c     line
      do j = 2,3                ! g from upper/lower line		
         
         do mu = 0,3
            qee(mu,j) = pew(mu)+p21(mu,j) ! e ve  emitted on upper line
            qaa(mu,j) = pa(mu)+p21(mu,j) ! A emitted on upper line
         enddo
         qee(4,j) = qee(0,j)**2-qee(1,j)**2-qee(2,j)**2-qee(3,j)**2
         qaa(4,j) = qaa(0,j)**2-qaa(1,j)**2-qaa(2,j)**2-qaa(3,j)**2

         prop_ee(1,j) = 1d0/qee(4,j)
         prop_ee(2,j) = 1d0/dcmplx(qee(4,j)-xm2(2),xmg(2))
         prop_ee(3,j) = 1d0/dcmplx(qee(4,j)-xm2(3),xmg(3))
         prop_uu(1,j) = 1d0/qaa(4,j)
         prop_uu(2,j) = 1d0/dcmplx(qaa(4,j)-xm2(2),xmg(2))
         prop_uu(3,j) = 1d0/dcmplx(qaa(4,j)-xm2(3),xmg(3))
      enddo


c  ---------------------------------------------------------------------
c
c     get the external quark spinors (including factor sqrt(2E) )
c     
      call psi0m(4,pbar(0,1),sign(1),psi)
c     
c     get the f-fbar currents (with no gluon radiation) 
c     J21^mu=jqq(mu,isig1,1,is1,0), J43^mu=jqq(mu,isig3,2,is3,0) 
c     
      call curr6(1,psi(1,-1,2),p(0,2),psi(1,-1,1),p(0,1),
     $     jqq(0,-1,1,is1,0))      
      call curr6(1,psi(1,-1,4),p(0,4),psi(1,-1,3),p(0,3),
     $     jqq(0,-1,2,is3,0))
c     
c     nomenclature: jqq(mu,		...	Lorentz index (0:3),(4:5) moment. info
c     hel,	...	quark helicity (+-1)
c     u/l,	...	upper(1)/lower(2) quark line
c     is,		...	current for qq(g),qbqb(g),or gqbq line
c     l)		...	g polarization (l=0:no g,l=1,2 with g)
c     
c     Get the gluon polarization vector and the gluon emission spinors
      do l = 1,2                ! 2 gluon polarizations
         call polvec(qbar,l,eps(0,l)) ! get gluon pol.vectors         
         
         do isig = -1,1,2	! fermion helicity 
            
            call ket2r(psi(1,isig,1),.true.,p(0,1),isig,q,eps(0,l),
     $           braketg(1,isig,1,l),pq(0,1)) ! |q,1>_l,isig1
            call bra2r(psi(1,isig,2),.true.,p(0,2),isig,q,eps(0,l),
     $           braketg(1,isig,2,l),pq(0,2)) ! <2,q|_l,isig2
            call ket2r(psi(1,isig,3),.true.,p(0,3),isig,q,eps(0,l),
     $           braketg(1,isig,3,l),pq(0,3)) ! |q,3>_l,isig3
            call bra2r(psi(1,isig,4),.true.,p(0,4),isig,q,eps(0,l),
     $           braketg(1,isig,4,l),pq(0,4)) ! <4,q|_l,isig4

         enddo

c     NOTATION: braketg(2 component spinor, fermion hel.,
c     fermion ID, gluon polarization)
         
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
         
c     gluon from lower line:
         call curr6(1,psi(1,-1,4),p(0,4),braketg(1,-1,3,l),pq(0,3),jh1)	
c     =   <4|gam.mu|q,3>_l,isig3
         call curr6(1,braketg(1,-1,4,l),pq(0,4),psi(1,-1,3),p(0,3),jh2)	
c     =   <4,q|gam.mu|3>_l,isig3
         do isig = -1,1,2
            do mu = 0,5
               jqq(mu,isig,2,is3,l) = jh1(mu,isig) + jh2(mu,isig)
c     = (<4|gam.mu|q,3>+<4,q|gam.mu|3>)_l,isig3
            enddo
         enddo
      enddo


c     -----------------------------------------------------------------------
c     
c     contract with vvtoww tensors to get Vertex-Vertex scattering diagrams
c     
      do l = 1,2                ! gluon polarization
         do isig = -1,1,2
            
            if (jlog2) then
               j = 2            ! g from upper line
               
               maw = contract_Tjj(awtowa(0,0,j),jqq(0,isig,1,is1,l),
     $              jqq(0,-1,2,is3,0))
               mzw = contract_Tjj(zwtowa(0,0,j),jqq(0,isig,1,is1,l),
     $              jqq(0,-1,2,is3,0))
               mwa = contract_Tjj(watowa(0,0,j),jqq(0,isig,2,is3,0),
     $              jqq(0,-1,1,is1,l))
               mwz = contract_Tjj(wztowa(0,0,j),jqq(0,isig,2,is3,0),
     $              jqq(0,-1,1,is1,l))
               
               do k = 1,2
                  mat(k,isig,-1,j,l,1) = 
     1                 maw*clr(ifl(1,k),1,isig)*clr(ifl(3,k),3,-1)
     2                 + mzw*clr(ifl(1,k),2,isig)*clr(ifl(3,k),3,-1)
               enddo
               do k = 3,4
                  mat(k,-1,isig,j,l,1) = 
     3                 mwa*clr(ifl(1,k),3,-1)*clr(ifl(3,k),1,isig)
     4                 + mwz*clr(ifl(1,k),3,-1)*clr(ifl(3,k),2,isig)
               enddo	    
            endif               !jlog2
            
            if (jlog3) then	
               j = 3            ! g from lower line

               maw = contract_Tjj(awtowa(0,0,j),jqq(0,isig,1,is1,0),
     $              jqq(0,-1,2,is3,l))
               mzw = contract_Tjj(zwtowa(0,0,j),jqq(0,isig,1,is1,0),
     $              jqq(0,-1,2,is3,l))
               mwa = contract_Tjj(watowa(0,0,j),jqq(0,isig,2,is3,l),
     $              jqq(0,-1,1,is1,0))
               mwz = contract_Tjj(wztowa(0,0,j),jqq(0,isig,2,is3,l),
     $              jqq(0,-1,1,is1,0))
               do k = 1,2
                  mat(k,isig,-1,j,l,1) = 
     1                 maw*clr(ifl(1,k),1,isig)*clr(ifl(3,k),3,-1)
     2                 + mzw*clr(ifl(1,k),2,isig)*clr(ifl(3,k),3,-1)
               enddo
               do k = 3,4
                  mat(k,-1,isig,j,l,1) = 
     3                 mwa*clr(ifl(1,k),3,-1)*clr(ifl(3,k),1,isig)
     4                 + mwz*clr(ifl(1,k),3,-1)*clr(ifl(3,k),2,isig)
               enddo
            endif               !jlog3
	    
         enddo                  !isig
      enddo                     ! l 
      
     
c     -------------------------------------------------------------------
c     -------------------------------------------------------------------
c     
c     prepare box diagrams: attach A to external spinors
c     
c     isig = +-1   : left- and righthanded spinors coupling to As
c     
c     Notation for virtual 2-component spinors and momenta
c     
c     W+ -> e+ve   attached to quark number i: psiwe(*,isig,i) with momentum fqwe(mu,i)
c     A            attached to quark number i: psia(*,isig,i) with momentum fqau(mu,i)
c     
c     the fermion current corresponding to a quark line with the real emitted W+
c     attached next to quark number i is stored in jlv(mu,isig,is,i,l); 
c     for l = 0 jwp/jwm doesn"t contain a gluon; for l = 1,2 jlv
c     includes sum of possible gluon couplings on quark leg i
c     
      l = 0                     ! no gluon 
      
      
      do i = 1,3,2

         if (i.eq.1) then
            is = is1
         else
            is = is3   
         endif   

c     W couples only to left-handed fermions:
	 isig = -1
         call ket2c(psi(1,isig,i),.true.,p(0,i),isig,qe,we,
     1        psiwe(1,isig,i),fqwe(0,i))
	 call bra2c(psi(1,isig,i+1),.true.,p(0,i+1),isig,qe,we,
     1        psiwe(1,isig,i+1),fqwe(0,i+1))
         
         call curr6(-1,psi(1,-1,i+1),p(0,i+1),
     1        psiwe(1,-1,i),fqwe(0,i), jew(0,-1,is,i,l)   )
         call curr6(-1,psiwe(1,-1,i+1),fqwe(0,i+1),
     1        psi(1,-1,i),p(0,i), jew(0,-1,is,i+1,l) )
         

c     A couples also to right-handed fermions:
         do isig = -1,1,2
            
            call ket2c(psi(1,isig,i),.true.,p(0,i),isig,qa,ea,
     1           psia(1,isig,i),fqau(0,i))
            call bra2c(psi(1,isig,i+1),.true.,p(0,i+1),isig,qa,ea,
     1           psia(1,isig,i+1),fqau(0,i+1))

         enddo                  !isig
         
         call curr6(1,psi(1,-1,i+1),p(0,i+1),
     1        psia(1,-1,i),fqau(0,i), jua(0,-1,is,i,l)   )
         call curr6(1,psia(1,-1,i+1),fqau(0,i+1),
     1        psi(1,-1,i),p(0,i), jua(0,-1,is,i+1,l) )
         
         
      enddo                     !i
c     
c -------------------------------------------------------------------
c 
c     now currents with gluon radiation (l=1:2):
      
      do i = 1,3,2              ! fermion ID (isigi=-1 or 1)  

         if (i.eq.1) then
            is = is1
         else
            is = is3   
         endif   

 	 do l = 1,2

c     start with W:
	    isig = -1 ! W
            call ket2c(braketg(1,isig,i,l),.false.,pq(0,i),
     $                 isig,qe,we,braketgwe(1,isig,i,l),pgwe(0,i))
            call bra2c(braketg(1,isig,i+1,l),.false.,pq(0,i+1),
     $             isig,qe,we,braketgwe(1,isig,i+1,l),pgwe(0,i+1))
	
	    call ket2r(psiwe(1,isig,i),.false.,fqwe(0,i),isig,
     $	    		q,eps(0,l),braketweg(1,isig,i,l),pweg(0,i))	 
            call bra2r(psiwe(1,isig,i+1),.false.,fqwe(0,i+1),isig,
     $	    		q,eps(0,l),braketweg(1,isig,i+1,l),pweg(0,i+1))      
  
! gluon radiation from fermion i / i+1
     	    call curr6(-1,psi(1,-1,i+1),p(0,i+1),
     $	    			braketgwe(1,-1,i,l),pgwe(0,i),jwegi) 			
      	    call curr6(-1,braketgwe(1,-1,i+1,l),pgwe(0,i+1),
     $	    			psi(1,-1,i),p(0,i),jwegii)		       
      	    call curr6(-1,psi(1,-1,i+1),p(0,i+1),
     $			        braketweg(1,-1,i,l),pweg(0,i),jgwei)
      	    call curr6(-1,braketweg(1,-1,i+1,l),pweg(0,i+1),	
     $			        psi(1,-1,i),p(0,i),jgweii)
    	    
! gluon radiation from fermion i+1 / i
     	    call curr6(-1,braketg(1,-1,i+1,l),pq(0,i+1),
     $	    			 psiwe(1,-1,i),fqwe(0,i),jweg0i)
     	    call curr6(-1,psiwe(1,-1,i+1),fqwe(0,i+1),
     $	    			braketg(1,-1,i,l),pq(0,i),jweg0ii)
     	 
	    do mu = 0,5
 		   jew(mu,isig,is,i,l) = jwegi(mu,isig)+
     $              jgwei(mu,isig)+jweg0i(mu,isig) ! We & g emission from i/i+1 line 
		   jew(mu,isig,is,i+1,l) = jwegii(mu,isig)+
     $			jgweii(mu,isig)+jweg0ii(mu,isig)
     
     	    enddo   !mu

c--------------
c
c     A (plus g )radiation:
c
  	  do isig = -1,1,2
	
            	call ket2c(braketg(1,isig,i,l),.false.,pq(0,i),
     $                     isig,qa,ea,braketgau(1,isig,i,l),pgau(0,i))
            	call bra2c(braketg(1,isig,i+1,l),.false.,pq(0,i+1),
     $                 isig,qa,ea,braketgau(1,isig,i+1,l),pgau(0,i+1))
		
		call ket2r(psia(1,isig,i),.false.,fqau(0,i),isig,
     $	    		    q,eps(0,l),braketaug(1,isig,i,l),paug(0,i))      
            	call bra2r(psia(1,isig,i+1),.false.,fqau(0,i+1),isig,
     $	    		q,eps(0,l),braketaug(1,isig,i+1,l),paug(0,i+1))  
     
             enddo              !isig    
    	    
! gluon radiation from fermion i / i+1
       	 
	  call curr6(1,psi(1,-1,i+1),p(0,i+1),
     $			     braketgau(1,-1,i,l),pgau(0,i),jaugi)		     
      	  call curr6(1,braketgau(1,-1,i+1,l),pgau(0,i+1),
     $			     psi(1,-1,i),p(0,i),jaugii) 		    
      	  call curr6(1,psi(1,-1,i+1),p(0,i+1),
     $	     		     braketaug(1,-1,i,l),paug(0,i),jgaui)
      	  call curr6(1,braketaug(1,-1,i+1,l),paug(0,i+1),     
     $			      psi(1,-1,i),p(0,i),jgauii)

 
! gluon radiation from fermion i+1 / i

    	  call curr6(1,braketg(1,-1,i+1,l),pq(0,i+1),
     $			      psia(1,-1,i),fqau(0,i),jaug0i)
     	  call curr6(1,psia(1,-1,i+1),fqau(0,i+1),
     $			     braketg(1,-1,i,l),pq(0,i),jaug0ii)
 
          do mu = 0,5
             do isig = -1,1,2
		   
                jua(mu,isig,is,i,l) = jaugi(mu,isig)+
     $               jgaui(mu,isig)+jaug0i(mu,isig)   
  		   
		   jua(mu,isig,is,i+1,l) = jaugii(mu,isig)+
     $			jgauii(mu,isig)+jaug0ii(mu,isig)   
	    	       	    	   
                enddo           !isig
             enddo              !mu       

          enddo                 ! l = 1,2
       enddo                    ! i loop


c ------------------------------------------------------------------------------
c -------------------------------------------------------------------
c
c     now calculate the Vertex-box diagrams; 
       do k = 1,2
          kk = 3-k
c     k=1 and kk=2 is for "box correction" to upper line
c     k=2 and kk=1 is for "box correction" to lower line
          
          if (kk.eq.1) then
             is = is1
          else
             is = is3
          endif	

          do isig = -1,1,2
             
!     VW-> e nu: (V = A or Z)
             
             do l = 0,2
	  	
                arg2 = 2+2*k+l*(2*k-3)*(l-3) ! arg = 1+k for l = 0
                arg = arg2/2    ! arg = 4-k for l = 1,2
		
                if (k.eq.1) then
                   call contract_T1j(awen(0,0,k,arg),jqq(0,isig,kk,is,l),epsaw) !for AW
                   call contract_T1j(zwen(0,0,k,arg),jqq(0,isig,kk,is,l),epszw) !for ZW
                   call contract_T2j(waen(0,0,k,arg),jqq(0,isig,kk,is,l),epswa) !for WA
                   call contract_T2j(wzen(0,0,k,arg),jqq(0,isig,kk,is,l),epswz) !for WZ
                else
                   call contract_T2j(awen(0,0,k,arg),jqq(0,isig,kk,is,l),epsaw) !for AW
                   call contract_T2j(zwen(0,0,k,arg),jqq(0,isig,kk,is,l),epszw) !for ZW  	 
                   call contract_T1j(waen(0,0,k,arg),jqq(0,isig,kk,is,l),epswa) !for WA
                   call contract_T1j(wzen(0,0,k,arg),jqq(0,isig,kk,is,l),epswz) !for WZ
                endif
                
! WW -> A:

                call contract_T2j(wwtoa(0,0,k,arg),jqq(0,isig,kk,is,l),epsww) !for WW 	  	
                do mu = 0,3     ! graphs with uncontracted V from VV"->mu+mu-
! (epsvuu) or from VV->e nu (epsvwe/epswve/epswee) 
! always: k=3/4 for quark type,arg = 2:3 (upper/lower),l=0:2

                   
                   epswuu(mu,isig,3,k,arg,l) = epsww(mu)*clr(3,3,isig) !k=3/4 for quark type
                   epswuu(mu,isig,4,k,arg,l) = epswuu(mu,isig,3,k,arg,l)
                   
                   epsawe(mu,isig,3,k,arg,l) = epsaw(mu)*clr(3,3,isig)  
                   epszwe(mu,isig,3,k,arg,l) = epszw(mu)*clr(3,3,isig)	       
                   epsawe(mu,isig,4,k,arg,l) = epsawe(mu,isig,3,k,arg,l)  
                   epszwe(mu,isig,4,k,arg,l) = epszwe(mu,isig,3,k,arg,l)
                   
                   epswae(mu,isig,3,k,arg,l) = epswa(mu)*clr(3,3,isig)  
                   epswze(mu,isig,3,k,arg,l) = epswz(mu)*clr(3,3,isig)	     
                   epswae(mu,isig,4,k,arg,l) = epswae(mu,isig,3,k,arg,l)  
                   epswze(mu,isig,4,k,arg,l) = epswze(mu,isig,3,k,arg,l)
                   
                   if (k.eq.2) then 
                      epswee(mu,isig,3,k,arg,l) =
     1                     epsaw(mu)*clr(3,1,isig)+epszw(mu)*clr(3,2,isig)  
                      epswee(mu,isig,4,k,arg,l) =
     1                     epsaw(mu)*clr(4,1,isig)+epszw(mu)*clr(4,2,isig)  
                   else         !k=1
                      epswee(mu,isig,3,k,arg,l) =
     1                     epswa(mu)*clr(3,1,isig)+epswz(mu)*clr(3,2,isig)  
                      epswee(mu,isig,4,k,arg,l) =
     1                     epswa(mu)*clr(4,1,isig)+epswz(mu)*clr(4,2,isig)  

                   endif

                enddo           !mu
             enddo              !l
          enddo                 !isig
       enddo                    !k
      
      do l = 0,2
         do isig = -1,1,2
            
            jj21au(isig,l) = dotcc(ea,jqq(0,isig,1,is1,l))
            jj21we(isig,l) = dotcc(we,jqq(0,isig,1,is1,l))
            
            jj43au(isig,l) = dotcc(ea,jqq(0,isig,2,is3,l))
            jj43we(isig,l) = dotcc(we,jqq(0,isig,2,is3,l))
            
         enddo                  !isig
      enddo                     ! l		 
c     
c --------------
c
c     now construct the contribution to the amplitude by current contraction 
c     virtual contributions are assembled in subroutine boxline
c     2 bosons attached to 12 line     
c

	   
c     gauge term for Z/W propagator:
      do isig3 = -1,1,2
         do l = 0,2
            do k=3,4            !up/down type
               do j = 2,3        
                  qepszwe(k,isig3,j,l) = -dotrc(qaa(0,j),
     $                 epszwe(0,isig3,k,1,j,l))*zm2i(2)	   
                  qepswee(k,isig3,j,l) = -dotrc(qaa(0,j),
     $                 epswee(0,isig3,k,1,j,l))*zm2i(3) 
                  qepswuu(k,isig3,j,l) = -dotrc(qee(0,j),
     $                 epswuu(0,isig3,k,1,j,l))*zm2i(3)	 		
               enddo            !j   	 
            enddo               !k
         enddo                  !l
      enddo                     !isig3
	   
      do isig1 = -1,1,2
         do isig3 = -1,1,2
            do l = 1,2
               
!     gluon from upper line:
               do k = 1,4       ! e+ ve from upper line 
                  
                  if (isig1.eq.-1) then
                     if (k.eq.1) then  
                        
                        m1we = dotcc(jew(0,isig1,is1,2,l),epswuu(0,isig3,ifl(3,k),1,2,0))
                        w1we = qepswuu(ifl(3,k),isig3,2,0)*jj21we(isig1,l)    
                        
                        m1e = (m1we+w1we)*clr(ifl(1,k),3,isig1)*clr(ifl(2,k),4,isig1)
                        m2e = 0d0

                     elseif (k.eq.2) then  
                        
                        m2we = dotcc(jew(0,isig1,is1,1,l),epswuu(0,isig3,ifl(3,k),1,2,0))
                        w2we = qepswuu(ifl(3,k),isig3,2,0)*jj21we(isig1,l)    
                        
                        m2e = (m2we-w2we)*clr(ifl(1,k),3,isig1)*clr(ifl(2,k),4,isig1)
                        m1e = 0d0
                        
                     else       ! k=3,4 			     
                        
                        m1e = 0d0 
                        
                        m2e = 0d0
                        
                     endif      !k
                     
                  else
                     m1e  =  0d0
                     m2e  =  0d0
                  endif         !isig1
                  
                  
!     A from upper line:
                  
                  if (k.lt.3.and.isig3.eq.-1) then !k=1,2
                     m1aau = dotcc(jua(0,isig1,is1,2,l),epsawe(0,isig3,ifl(3,k),1,2,0))
                     m1azu = dotcc(jua(0,isig1,is1,2,l),epszwe(0,isig3,ifl(3,k),1,2,0))
                     
                     z1azu = qepszwe(ifl(3,k),isig3,2,0)*jj21au(isig1,l)  
                     
                     mw1aau =  m1aau*clr(ifl(1,k),1,isig1)*clr(ifl(2,k),1,isig1)
                     mw1azu = (m1azu+z1azu)*clr(ifl(1,k),2,isig1)*clr(ifl(2,k),1,isig1)
                     
                     m1u = mw1aau+mw1azu
                     
                     m2aau = dotcc(jua(0,isig1,is1,1,l),epsawe(0,isig3,ifl(3,k),1,2,0))
                     m2azu = dotcc(jua(0,isig1,is1,1,l),epszwe(0,isig3,ifl(3,k),1,2,0))
                     
                     mw2aau =  m2aau*clr(ifl(2,k),1,isig1)*clr(ifl(1,k),1,isig1)
                     mw2azu = (m2azu-z1azu)*clr(ifl(2,k),2,isig1)*clr(ifl(1,k),1,isig1)
                     
                     m2u = mw2aau+mw2azu
                     
                  elseif (k.gt.2) then !k=3,4
                     m1awu = dotcc(jua(0,isig1,is1,2,l),epswee(0,isig3,ifl(3,k),1,2,0))
                     
                     z1awu = qepswee(ifl(3,k),isig3,2,0)*jj21au(isig1,l)  
                     
                     m1wau = (m1awu+z1awu)*clr(ifl(1,k),3,isig1)*clr(ifl(2,k),1,isig1)
                     
                     m1u = m1wau
                     
                     m2awu = dotcc(jua(0,isig1,is1,1,l),epswee(0,isig3,ifl(3,k),1,2,0))
                     
                     m2wau = (m2awu-z1awu)*clr(ifl(2,k),3,isig1)*clr(ifl(1,k),1,isig1)
                     
                     m2u = m2wau
                     
                  else
                     m1u = 0d0
                     m2u = 0d0	  
                  endif         !k
                  
                  mat(k,isig1,isig3,2,l,2) =  -(m1e-m1u+m2e-m2u)
                  
               enddo            !k
               
c-----------
!     gluon from lower line:
               do k = 1,4       ! e+ ve from upper line 
                  
                  if (isig1.eq.-1) then
                     if (k.eq.1) then  
                        
                        m1we = dotcc(jew(0,isig1,is1,2,0),epswuu(0,isig3,ifl(3,k),1,3,l))	   
                        w1we = qepswuu(ifl(3,k),isig3,3,l)*jj21we(isig1,0)    
                        
                        m1e = (m1we+w1we)*clr(ifl(1,k),3,isig1)*clr(ifl(2,k),4,isig1)
                        m2e = 0d0

                     elseif (k.eq.2) then  
                        
                        m2we = dotcc(jew(0,isig1,is1,1,0),epswuu(0,isig3,ifl(3,k),1,3,l))	   
                        w2we = qepswuu(ifl(3,k),isig3,3,l)*jj21we(isig1,0)    
                        
                        m2e = (m2we-w2we)*clr(ifl(1,k),3,isig1)*clr(ifl(2,k),4,isig1)
                        m1e = 0d0
                        
                     else       ! k=3,4 			     
                        
                        m1e = 0d0
                        
                        m2e = 0d0 
                        
                     endif      !k
                     
                  else
                     m1e  =  0d0
                     m2e  =  0d0
                  endif         !isig1
                  
                  
! A from upper line:
                  
                  if (k.lt.3.and.isig3.eq.-1) then !k=1,2
                     m1aau = dotcc(jua(0,isig1,is1,2,0),epsawe(0,isig3,ifl(3,k),1,3,l))
                     m1azu = dotcc(jua(0,isig1,is1,2,0),epszwe(0,isig3,ifl(3,k),1,3,l))

                     z1azu = qepszwe(ifl(3,k),isig3,3,l)*jj21au(isig1,0)  

                     
                     mw1aau =  m1aau*clr(ifl(1,k),1,isig1)*clr(ifl(2,k),1,isig1)
                     mw1azu = (m1azu+z1azu)*clr(ifl(1,k),2,isig1)*clr(ifl(2,k),1,isig1)
                     
                     m1u = mw1aau+mw1azu
                     
                     m2aau = dotcc(jua(0,isig1,is1,1,0),epsawe(0,isig3,ifl(3,k),1,3,l))
                     m2azu = dotcc(jua(0,isig1,is1,1,0),epszwe(0,isig3,ifl(3,k),1,3,l))

                     
                     mw2aau =  m2aau*clr(ifl(2,k),1,isig1)*clr(ifl(1,k),1,isig1)
                     mw2azu = (m2azu-z1azu)*clr(ifl(2,k),2,isig1)*clr(ifl(1,k),1,isig1)
                     
                     m2u = mw2aau+mw2azu
                     
                  elseif (k.gt.2) then !k=3,4
                     m1awu = dotcc(jua(0,isig1,is1,2,0),epswee(0,isig3,ifl(3,k),1,3,l))
                     
                     z1awu = qepswee(ifl(3,k),isig3,3,l)*jj21au(isig1,0)  
                     
                     m1wau = (m1awu+z1awu)*clr(ifl(1,k),3,isig1)*clr(ifl(2,k),1,isig1)
                     
                     m1u = m1wau
                     
                     m2awu = dotcc(jua(0,isig1,is1,1,0),epswee(0,isig3,ifl(3,k),1,3,l))
                     
                     m2wau = (m2awu-z1awu)*clr(ifl(2,k),3,isig1)*clr(ifl(1,k),1,isig1)
                     
                     m2u = m2wau
                     
                  else
                     m1u = 0d0
                     m2u = 0d0	  
                  endif         !k
                  
                  mat(k,isig1,isig3,3,l,2) =  -(m1e-m1u+m2e-m2u)	   
                  
               enddo            !k
               
            enddo               ! l 	   
         enddo                  !isig3
      enddo                     ! isig1  
      

c------------------------------------------------------------
      
c     repeat the same for 2 bosons attached to 34 line:
c     
c     gauge term for Z/W propagator:
      do isig1 = -1,1,2  
         do l = 0,2
            do j = 2,3
               do k=3,4         !up/down type
                  qepswze(k,isig1,j,l) = dotrc(qee(0,j),
     $                 epswze(0,isig1,k,2,j,l))*zm2i(2) 
                  qepswee(k,isig1,j,l) = dotrc(qee(0,j),
     $                 epswee(0,isig1,k,2,j,l))*zm2i(3) 
                  qepswuu(k,isig1,j,l) = dotrc(qaa(0,j),
     $                 epswuu(0,isig1,k,2,j,l))*zm2i(3)	   
               enddo            !k	   
            enddo               !j
         enddo                  !l
      enddo                     ! isig1  
      
      do isig1 = -1,1,2
         do isig3 = -1,1,2	   
            do l = 1,2
               
!     gluon from upper line:
               do k = 1,4       ! e+ve from lower line:
                  
                  if(isig3.eq.-1) then
                     if (k.lt.3) then !k=1,2
                        
                        m3e = 0d0
                        
                        m4e = 0d0
                        
                     elseif (k.eq.4) then
                        
                        m3e = 0d0

                        m4we = dotcc(jew(0,isig3,is3,3,0),epswuu(0,isig1,ifl(1,k),2,2,l))
                        w4we = qepswuu(ifl(1,k),isig1,2,l)*jj43we(isig3,0)    	   	   	   	   	   
                        m4e = (m4we-w4we)*clr(ifl(4,k),4,isig3)*clr(ifl(3,k),3,isig3)

                     elseif (k.eq.3) then
                        
                        m4e = 0d0

                        m3we = dotcc(jew(0,isig3,is3,4,0),epswuu(0,isig1,ifl(1,k),2,2,l))
                        w3we = qepswuu(ifl(1,k),isig1,2,l)*jj43we(isig3,0)    	   	   	   	   	   
                        m3e = (m3we+w3we)*clr(ifl(4,k),3,isig3)*clr(ifl(3,k),3,isig3)
                        
                     endif
                     
                  else
                     m3e = 0d0
                     m4e = 0d0
                  endif         !isig3  
                  
                  
! mu+mu- from lower line:
                  
                  if (k.lt.3.and.isig3.eq.-1) then ! k =1,2		
                     m3awu = dotcc(jua(0,isig3,is3,4,0),epswee(0,isig1,ifl(1,k),2,2,l))
                     
                     z3awu = qepswee(ifl(1,k),isig1,2,l)*jj43au(isig3,0)  
                     
                     m3wau = (m3awu+z3awu)*clr(ifl(3,k),3,isig3)*clr(ifl(4,k),1,isig3)
                     
                     m3u = m3wau
                     
                     m4awu = dotcc(jua(0,isig3,is3,3,0),epswee(0,isig1,ifl(1,k),2,2,l))
                     
                     m4wau = (m4awu-z3awu)*clr(ifl(4,k),3,isig3)*clr(ifl(3,k),1,isig3)
	     	   
                     m4u = m4wau
                     
                  elseif (k.gt.2) then !k=3,4
                     m3aau = dotcc(jua(0,isig3,is3,4,0),epswae(0,isig1,ifl(1,k),2,2,l))
                     m3azu = dotcc(jua(0,isig3,is3,4,0),epswze(0,isig1,ifl(1,k),2,2,l))

                     z3azu = qepswze(ifl(1,k),isig1,2,l)*jj43au(isig3,0)  
                     
                     mw3aau =  m3aau*clr(ifl(3,k),1,isig3)*clr(ifl(4,k),1,isig3)
                     mw3azu = (m3azu+z3azu)*clr(ifl(3,k),2,isig3)*clr(ifl(4,k),1,isig3)

                     m3u = mw3aau+mw3azu
                     
                     m4aau = dotcc(jua(0,isig3,is3,3,0),epswae(0,isig1,ifl(1,k),2,2,l))
                     m4azu = dotcc(jua(0,isig3,is3,3,0),epswze(0,isig1,ifl(1,k),2,2,l))
                     
                     mw4aau =  m4aau*clr(ifl(4,k),1,isig3)*clr(ifl(3,k),1,isig3)
                     mw4azu = (m4azu-z3azu)*clr(ifl(4,k),2,isig3)*clr(ifl(3,k),1,isig3)
                     
                     m4u = mw4aau+mw4azu
                     
                  else
                     m3u = 0d0
                     m4u = 0d0
                  endif         !k
                  
                  mat(k,isig1,isig3,2,l,3) = -(m3e-m3u+m4e-m4u)
                  
               enddo            !k
               
!     gluon from lower line:
               do k = 1,4       ! e+ve from lower line:
                  
                  if(isig3.eq.-1) then
                     if (k.lt.3) then !k=1,2
                        
                        m3e = 0d0
                        
                        m4e = 0d0
                        
                     elseif (k.eq.4) then
                        
                        m3e = 0d0

                        m4we = dotcc(jew(0,isig3,is3,3,l),epswuu(0,isig1,ifl(1,k),2,3,0))
                        w4we = qepswuu(ifl(1,k),isig1,3,0)*jj43we(isig3,l)    	   	   	   	   	   
                        m4e = (m4we-w4we)*clr(ifl(4,k),4,isig3)*clr(ifl(3,k),3,isig3)

                     elseif (k.eq.3) then
                        
                        m4e = 0d0

                        m3we = dotcc(jew(0,isig3,is3,4,l),epswuu(0,isig1,ifl(1,k),2,3,0))
                        w3we = qepswuu(ifl(1,k),isig1,3,0)*jj43we(isig3,l)    	   	   	   	   	   
                        m3e = (m3we+w3we)*clr(ifl(4,k),3,isig3)*clr(ifl(3,k),3,isig3)
                        
                     endif
                     
                  else
                     m3e = 0d0
                     m4e = 0d0
                  endif         !isig3  
                  
                  
! A from lower line:
                  
                  if (k.lt.3.and.isig3.eq.-1) then ! k =1,2		
                     m3awu = dotcc(jua(0,isig3,is3,4,l),epswee(0,isig1,ifl(1,k),2,3,0))
                     
                     z3awu = qepswee(ifl(1,k),isig1,3,0)*jj43au(isig3,l)  
                     
                     m3wau = (m3awu+z3awu)*clr(ifl(3,k),3,isig3)*clr(ifl(4,k),1,isig3)
                     
                     m3u = m3wau
                     
                     m4awu = dotcc(jua(0,isig3,is3,3,l),epswee(0,isig1,ifl(1,k),2,3,0))
                     
                     m4wau = (m4awu-z3awu)*clr(ifl(4,k),3,isig3)*clr(ifl(3,k),1,isig3)
                     
                     m4u = m4wau
                     
                  elseif (k.gt.2) then !k=3,4
                     m3aau = dotcc(jua(0,isig3,is3,4,l),epswae(0,isig1,ifl(1,k),2,3,0))
                     m3azu = dotcc(jua(0,isig3,is3,4,l),epswze(0,isig1,ifl(1,k),2,3,0))

                     z3azu = qepswze(ifl(1,k),isig1,3,0)*jj43au(isig3,l)  
                     
                     mw3aau =  m3aau*clr(ifl(3,k),1,isig3)*clr(ifl(4,k),1,isig3)
                     mw3azu = (m3azu+z3azu)*clr(ifl(3,k),2,isig3)*clr(ifl(4,k),1,isig3)
                     
                     m3u = mw3aau+mw3azu
                     
                     m4aau = dotcc(jua(0,isig3,is3,3,l),epswae(0,isig1,ifl(1,k),2,3,0))
                     m4azu = dotcc(jua(0,isig3,is3,3,l),epswze(0,isig1,ifl(1,k),2,3,0))
                     
                     mw4aau =  m4aau*clr(ifl(4,k),1,isig3)*clr(ifl(3,k),1,isig3)
                     mw4azu = (m4azu-z3azu)*clr(ifl(4,k),2,isig3)*clr(ifl(3,k),1,isig3)
                     
                     m4u = mw4aau+mw4azu
                     
                  else
                     m3u = 0d0
                     m4u = 0d0
                  endif         !k
                  
                  mat(k,isig1,isig3,3,l,3) = -(m3e-m3u+m4e-m4u)
               enddo            !k
               
            enddo               !l	   
         enddo                  !isig3
      enddo                     ! isig1  

c     ---------------------   
c     
c     take special care of processes with incoming gluons:

      if (jmin.eq.3.or.jmax.eq.2) then
         if (jmin.eq.3) jout = 2
         if (jmax.eq.2) jout = 3
         do k = 1,4
	    do isig1 = -1,1,2
               do isig3 = -1,1,2
                  do l = 1,2
                     do i = 2,3
                        mat(k,isig1,isig3,jout,l,i) = 0d0
                     enddo
                  enddo
               enddo
	    enddo
         enddo	  
      endif	
      

c     
c     -----------------------------------------------------------------------
c     -----------------------------------------------------------------------
c     
c     next come the W-->WA currents attached to the quark lines. For the virtual 
c     corrections the most effective structure is the contraction of two polarization 
c     vectors with one fermion line. First build these effective polarization vectors
c     from the current wtowa(mu)
c     
c     NOTE: the wwv currents are NOT conserved. Hence it is necessary 
c     to consider q^mu * q^nu/m_W^2  terms in the W boson propagator

      do i = 1,3,2
	 
	 if (i.eq.1) then
            is = is1
	 else
            is = is3
	 endif	       
	 isig = -1              ! only left-handed q-W coupling
         
         call ket2c(psi(1,isig,i),.true.,p(0,i),isig,pwa,wtowa,
     1        psiw(1,isig,is,i),fq(0,i))
         call bra2c(psi(1,isig,i+1),.true.,p(0,i+1),isig,pwa,wtowa,
     1        psiw(1,isig,is,i+1),fq(0,i+1))
         
!     gluon radiation from fermion i / i+1
         do l = 1,2
            call ket2c(braketg(1,isig,i,l),.false.,pq(0,i),
     $           isig,pwa,wtowa,braketgw(1,isig,l,i),pgw(0,i))
!     |wwz,q,i>_l,isigi
            call ket2r(psiw(1,isig,is,i),.false.,fq(0,i),isig,
     $           q,eps(0,l),braketwg(1,isig,l,i),pwg(0,i))      
!     |q,wwz,i>_l,isigi
            
            call bra2c(braketg(1,isig,i+1,l),.false.,pq(0,i+1),
     $           isig,pwa,wtowa,braketgw(1,isig,l,i+1),pgw(0,i+1))
!     <wwz,q,i+1|_l,isig
            call bra2r(psiw(1,isig,is,i+1),.false.,fq(0,i+1),isig,
     $           q,eps(0,l),braketwg(1,isig,l,i+1),pwg(0,i+1))      
!     <q,wwz,i+1|_l,isig
         enddo

         
c     W from fermion i first:			
         
 	 do l = 1,2
!     W from fermion i; gluon radiation from fermion i 
            call curr6(-1,psi(1,-1,i+1),p(0,i+1),
     $           braketgw(1,-1,l,i),pgw(0,i),jwg) 			
            call curr6(-1,psi(1,-1,i+1),p(0,i+1),
     $           braketwg(1,-1,l,i),pwg(0,i),jgw)
!     W from fermion i; gluon radiation from fermion i+1 
            call curr6(-1,braketg(1,-1,i+1,l),pq(0,i+1),
     $           psiw(1,-1,is,i),fq(0,i),jg0)

            do mu = 0,5
               do isig = -1,1,2 ! component for isig = +1 -> zero 
                  jw(mu,isig,is,i,l) = jwg(mu,isig)+
     $                 jgw(mu,isig)+jg0(mu,isig) ! W+g emission from i/i+1 line 
               enddo            !isig
            enddo               !mu

            
!     gluon radiation from the non-i line:
            call curr6(-1,psi(1,-1,i+1),p(0,i+1),
     1           psiw(1,-1,is,i),fq(0,i), jw(0,-1,is,i,0) )
            
c     now do W from fermion i+1:			

!     W from fermion i+1; gluon radiation from fermion i+1 
            call curr6(-1,braketgw(1,-1,l,i+1),pgw(0,i+1),
     $           psi(1,-1,i),p(0,i),jwg)		       
            call curr6(-1,braketwg(1,-1,l,i+1),pwg(0,i+1),
     $           psi(1,-1,i),p(0,i),jgw)

!     W from fermion i+1, gluon radiation from fermion i 
            call curr6(-1,psiw(1,-1,is,i+1),fq(0,i+1),
     $           braketg(1,-1,i,l),pq(0,i),jg0)

            do mu = 0,5
               do isig = -1,1,2
                  jw(mu,isig,is,i+1,l) = jwg(mu,isig)+
     $                 jgw(mu,isig)+jg0(mu,isig) ! A+g emission from i/i+1 line 
               enddo            ! isig 
            enddo               ! mu					 

!     gluon radiation from the non-i line:
            call curr6(-1,psiw(1,-1,is,i+1),fq(0,i+1),
     1           psi(1,-1,i),p(0,i), jw(0,-1,is,i+1,0) )
         enddo                  ! l=1,2 loop
      enddo                     !i
      
      
c---------------------------
      
      
c     extra terms for massive W propagator:      
      psiwg = -dotrc(pwa,wtowa(0))*zm2i(3)*clr(3,3,-1)
      
!     compute jwg(boson from upper/lower line,isig1,isig3,boson id,j,l)
      do isig1 = -1,1,2
         do isig3 = -1,1,2
            do kb = 1,3         ! boson id      
               do l = 1,2  
!     gluon from upper line:
                  jwgg(1,isig1,isig3,kb,2,l) = psiwg*dotcc(jqq(0,isig1,1,is1,l),
     $                 jqq(0,isig3,2,is3,0))*prop43(kb,2)
                  jwgg(2,isig1,isig3,kb,2,l) = psiwg*dotcc(jqq(0,isig1,1,is1,l),
     $                 jqq(0,isig3,2,is3,0))*prop21(kb,2)
                  
!     gluon from lower line:
                  jwgg(1,isig1,isig3,kb,3,l) = psiwg*dotcc(jqq(0,isig1,1,is1,0),
     $                 jqq(0,isig3,2,is3,l))*prop43(kb,3)
                  jwgg(2,isig1,isig3,kb,3,l) = psiwg*dotcc(jqq(0,isig1,1,is1,0),
     $                 jqq(0,isig3,2,is3,l))*prop21(kb,3)
               enddo	
            enddo
         enddo
      enddo
      
c---------------------------

      do j = 2,3
         propw(1,j) = clr(3,3,-1)**2*prop21(3,j)
         propw(2,j) = clr(3,3,-1)**2*prop43(3,j)
      enddo                     !j      
      
      do isig3 = -1,1,2
         isig1 = -1             ! only left-handed quarks couple to W
         

c     box correction to upper line: polarization vectors are 
c     jqq(mu,isig3,2,is3,l)=j43 with momentum    p43        and
c     wwv(mu)              with momentum    pwz
         
         do i = 1,2
            do l = 1,2
               
               mw(i,2,l) =  dotcc(jw(0,isig1,is1,3-i,l),jqq(0,isig3,2,is3,0))
               mw(i,3,l) =  dotcc(jw(0,isig1,is1,3-i,0),jqq(0,isig3,2,is3,l))
               
            enddo               !l
         enddo                  !i
         
         
         do k = 3,4
            do j = jmin,jmax
               propt1(isig1,isig3,k,2,j) = 
     1              clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3)*prop43(1,j) 
     2              + clr(ifl(1,k),2,isig1)*clr(ifl(3,k),2,isig3)*prop43(2,j)
               propt2(isig1,isig3,k,2,j) = 
     1              clr(ifl(2,k),1,isig1)*clr(ifl(3,k),1,isig3)*prop43(1,j) 
     2              + clr(ifl(2,k),2,isig1)*clr(ifl(3,k),2,isig3)*prop43(2,j)
            enddo               !j
         enddo                  !k
         
         
         do l = 1,2    
            
            if (isig3.eq.-1) then ! only left-handed quarks couple to W
               do j = jmin,jmax
                  do k = 1,2
                     kk = 3-k
                     mat(k,isig1,isig3,j,l,4) = propw(2,j) *
     1                    mw(k,j,l)*clr(3,3,isig1) 
     1                    - jwgg(1,isig1,isig3,3,j,l)*clr(3,3,isig1)*
     1                    clr(3,3,isig3)*(3-2*kk)	
                  enddo         !k
               enddo            !j
            endif               !isig3
            
            
            do j = jmin,jmax
               do k = 3,4 
                  mat(k,isig1,isig3,j,l,4) = propt2(isig1,isig3,k,2,j) *
     1                 mw(2,j,l)*clr(ifl(1,k),3,isig1)   
     1                 + propt1(isig1,isig3,k,2,j) *
     1                 mw(1,j,l)*clr(ifl(2,k),3,isig1)        
     1                 + jwgg(1,isig1,isig3,1,j,l)*(clr(ifl(1,k),1,isig1)
     1                 -clr(ifl(2,k),1,isig1))*clr(ifl(3,k),1,isig3)
     1                 + jwgg(1,isig1,isig3,2,j,l)*(clr(ifl(1,k),2,isig1)
     1                 -clr(ifl(2,k),2,isig1))*clr(ifl(3,k),2,isig3)

               enddo            !k
            enddo               !j
         enddo                  !l
      enddo                     !isig3	    
      
      
c     ------------------	    
      
      do isig1 = -1,1,2
         isig3 = -1             ! only left-handed quarks couple to W
         
c     box correction to lower line: polarization vectors are 
c     jqq(mu,isig1,1,is1,l)=j21 with momentum p21 and
c     wtowa(mu) with momentum pwa
         
         do i = 1,2
            do l = 1,2
               
               mw(i,2,l) =  dotcc(jw(0,isig3,is3,5-i,0),jqq(0,isig1,1,is1,l))
               mw(i,3,l) =  dotcc(jw(0,isig3,is3,5-i,l),jqq(0,isig1,1,is1,0))
               
            enddo               !l
         enddo                  !i

         do k = 1,2
            do j = jmin,jmax
               propt1(isig1,isig3,k,1,j) = 
     1              clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3)*prop21(1,j) 
     2              + clr(ifl(1,k),2,isig1)*clr(ifl(3,k),2,isig3)*prop21(2,j)
               propt2(isig1,isig3,k,1,j) = 
     1              clr(ifl(1,k),1,isig1)*clr(ifl(4,k),1,isig3)*prop21(1,j) 
     2              + clr(ifl(1,k),2,isig1)*clr(ifl(4,k),2,isig3)*prop21(2,j)
            enddo               !j
         enddo                  !k

         do l = 1,2

            do k = 1,2 
               do j = jmin,jmax
                  mat(k,isig1,isig3,j,l,5) = propt2(isig1,isig3,k,1,j) *
     1                 mw(2,j,l)*clr(ifl(3,k),3,isig3)   
     1                 + propt1(isig1,isig3,k,1,j) *
     1                 mw(1,j,l)*clr(ifl(4,k),3,isig3)   
     1                 + jwgg(2,isig1,isig3,1,j,l)*(clr(ifl(3,k),1,isig3)-clr(ifl(4,k),1,isig3))*
     1                 clr(ifl(1,k),1,isig1)
     1                 + jwgg(2,isig1,isig3,2,j,l)*(clr(ifl(3,k),2,isig3)-clr(ifl(4,k),2,isig3))*
     1                 clr(ifl(1,k),2,isig1)
               enddo            !j
            enddo
            

            if (isig1.eq.-1) then ! only left-handed quarks couple to W
               do j = jmin,jmax
                  do k = 3,4
                     kk = 5-k
                     kk = k-2
                     mat(k,isig1,isig3,j,l,5) = propw(1,j) *
     1                    mw(kk,j,l)*clr(3,3,isig1)
     1                    - jwgg(2,isig1,isig3,3,j,l)*
     1                    clr(3,3,isig1)*clr(3,3,isig3)*(-7+2*k)
                  enddo         !j
               enddo            !k
            endif               !isig1
         enddo                  !l
      enddo                     !isig1


c------------------------------------------------------------------------------
c------------------------------------------------------------------------------
c     
c     next do the box-box graphs with one boson (W/A) emitted from the 
c     upper and the other from the lower line. 
c     
c     
      do il = 1,2               ! two possible g pols.
         
         do j = jmin,jmax       ! g from upper or lower line
            
            if (j.eq.2) then
               l = il           ! loop over pol. of g from upper line
               ll = 0           ! no g from lower line
            else                !j=3
               ll = il          ! loop over pol. of g from lower line
               l = 0            ! no g from upper line
            endif
c     
            
c     W from upper line, V from lower line:
            do isig3 = -1,1,2
               isig1 = -1       !(only left-handed fermions couple to W)
               
               do idu = 1,2     ! fermion id for external W coupling
                  do idl = 3,4  ! fermion id for external V coupling
                     mbwa(idu,idl,isig1,isig3,j,il) = dotcc(jew(0,isig1,is1,idu,l),
     $                    jua(0,isig3,is3,idl,ll))*clr(3,3,isig1)
                  enddo         !idl
               enddo            !idu
            enddo               !isig3
            
c     W from lower line, V from upper line:
            
            do isig1 = -1,1,2
               isig3 = -1
c     
               do idu = 1,2     ! fermion id for external V coupling
                  do idl = 3,4  ! fermion id for external W coupling
                     mbaw(idu,idl,isig1,isig3,j,il) = dotcc(jew(0,isig3,is3,idl,ll),
     $                    jua(0,isig1,is1,idu,l))*clr(3,3,isig3)
                  enddo         !idl
               enddo            !idu
            enddo               !isig3

c     for the q^mu*q^nu/M_V^2 terms in the gauge boson propagators we need
            
            i = 3               ! type of exchanged boson (Z,W,A)  (Z/A contributions cancel)
            do isig3 = -1,1,2
               isig1 = -1
               
               zwa(isig1,isig3,i,j,il) = jj21we(isig1,l)*jj43au(isig3,ll)*
     $              zm2i(i)*clr(3,3,isig1)

            enddo               !isig3
            do isig1 = -1,1,2
               isig3 = -1
               
               zaw(isig1,isig3,i,j,il) = jj21au(isig1,l)*jj43we(isig3,ll)*
     $              zm2i(i)*clr(3,3,isig3)
               
            enddo               !isig1
         enddo                  !j
      enddo                     !il

c     -----------------------------
c     
c     now construct amplitudes:	
c     
      do j = jmin,jmax
         do l = 1,2
            
            
            do isig1 = -1,1,2
               isig3 = -1

               do k = 1,2
                  
                  if (isig1.eq.-1) then  
                     kk = 3-k
                     mbwv = (mbwa(kk,3,isig1,isig3,j,l)*clr(ifl(3,kk),1,isig3)+
     $                    mbwa(kk,4,isig1,isig3,j,l)*clr(ifl(4,kk),1,isig3)
     $                    )*clr(3,3,isig1)*clr(3,3,isig3)*prop_ee(3,j)
                     
                     gbwv = ( zwa(isig1,isig3,3,j,l)*
     $                    (clr(ifl(3,kk),1,isig3)-clr(ifl(4,kk),1,isig3))
     $                    )*clr(3,3,isig1)*clr(3,3,isig3)*prop_ee(3,j)
     $                    *(3-2*kk) ! +/- for k=2/1	
                     
                  else
                     mbwv = 0d0
                     gbwv = 0d0
                  endif

                  mbvw = 0d0
                  gbvw = 0d0
                  do kb = 1,2	
                     mbvw = mbvw +(
     $                    mbaw(1,3,isig1,isig3,j,l)*clr(ifl(1,k),1,isig1)*
     $                    clr(ifl(2,k),kb,isig1)*clr(ifl(4,k),kb,isig3)+
     $                    mbaw(2,3,isig1,isig3,j,l)*clr(ifl(2,k),1,isig1)*
     $                    clr(ifl(1,k),kb,isig1)*clr(ifl(4,k),kb,isig3)+
     $                    mbaw(1,4,isig1,isig3,j,l)*clr(ifl(1,k),1,isig1)*
     $                    clr(ifl(2,k),kb,isig1)*clr(ifl(3,k),kb,isig3)+
     $                    mbaw(2,4,isig1,isig3,j,l)*clr(ifl(2,k),1,isig1)*
     $                    clr(ifl(1,k),kb,isig1)*clr(ifl(3,k),kb,isig3) )
     $                    *prop_uu(kb,j)
                     
                  enddo         !kb
                  
                  mat(k,isig1,isig3,j,l,6) = -(mbwv+mbvw+gbvw+gbwv)
                  
                  
               enddo            !k
            enddo               !isig1


c     --------------
            
            
            do isig3 = -1,1,2
               isig1 = -1
               do k = 3,4
                  
                  if (isig3.eq.-1) then  
                     kk = 7-k
                     mbvw = (mbaw(1,kk,isig1,isig3,j,l)*clr(ifl(1,kk),1,isig1)+
     $                    mbaw(2,kk,isig1,isig3,j,l)*clr(ifl(2,kk),1,isig1)
     $                    )*clr(3,3,isig3)*clr(3,3,isig1)*prop_uu(3,j)
                     
                     gbvw = ( zaw(isig1,isig3,3,j,l)*
     $                    (clr(ifl(1,kk),1,isig1)-clr(ifl(2,kk),1,isig1))
     $                    )*clr(3,3,isig3)*clr(3,3,isig1)*prop_uu(3,j)
     $                    *(7-2*kk) ! +/- for k=4/3	
                     
                  else
                     mbvw = 0d0
                     gbvw = 0d0
                  endif	

                  mbwv = 0d0
                  gbwv = 0d0
                  do kb = 1,2	
                     mbwv = mbwv +(
     $                    mbwa(1,3,isig1,isig3,j,l)*clr(ifl(3,k),1,isig3)*
     $                    clr(ifl(2,k),kb,isig1)*clr(ifl(4,k),kb,isig3)+
     $                    mbwa(2,3,isig1,isig3,j,l)*clr(ifl(4,k),1,isig3)*
     $                    clr(ifl(1,k),kb,isig1)*clr(ifl(4,k),kb,isig3)+
     $                    mbwa(1,4,isig1,isig3,j,l)*clr(ifl(3,k),1,isig3)*
     $                    clr(ifl(2,k),kb,isig1)*clr(ifl(3,k),kb,isig3)+
     $                    mbwa(2,4,isig1,isig3,j,l)*clr(ifl(4,k),1,isig3)*
     $                    clr(ifl(1,k),kb,isig1)*clr(ifl(3,k),kb,isig3) )
     $                    *prop_ee(kb,j)
                     
                     gbwv = 0d0

                  enddo         !kb


                  mat(k,isig1,isig3,j,l,6) = -(mbwv+mbvw+gbvw+gbwv)
                  
               enddo            !k
            enddo               !isig
         enddo                  !j
      enddo                     !l


c------------------------------------------------------------------------------
c------------------------------------------------------------------------------
c     
c     and now, finally, the pentagon contributions, i.e. three bosons emitted 
c     from the same quark line
c     
c     compute matrix elements for pentagon on upper line:
c     
c     need m5(pos.j43=1:3,vext1=1:3,vext2=1:3), id of internal boson arbitrary 


      do isig3 = -1,1,2 
	 isig1 = -1

! g from upper line, but not in bkjqq:
         call ket2c(psi(1,isig1,1),.true.,p(0,1),isig1,p43(0,2),
     1        jqq(0,isig3,2,is3,0),bkjqq(1,isig1,isig3,1,2,0),bq(0,1))
         call bra2c(psi(1,isig1,2),.true.,p(0,2),isig1,p43(0,2),
     1        jqq(0,isig3,2,is3,0),bkjqq(1,isig1,isig3,2,2,0),bq(0,2))
         
         do l = 1,2
            
! g from lower line:
            call ket2c(psi(1,isig1,1),.true.,p(0,1),isig1,p43(0,3),
     1           jqq(0,isig3,2,is3,l),bkjqq(1,isig1,isig3,1,3,l),dummy)
            call bra2c(psi(1,isig1,2),.true.,p(0,2),isig1,p43(0,3),
     1           jqq(0,isig3,2,is3,l),bkjqq(1,isig1,isig3,2,3,l),dummy)
            
! g from upper line, in gbkjqq or bkjqqg:
            call ket2c(braketg(1,isig1,1,l),.false.,pq(0,1),isig1,p43(0,2),
     $		 jqq(0,isig3,2,is3,0),gbkjqq(1,isig1,isig3,1,2,l),dummy)
            call bra2c(braketg(1,isig1,2,l),.false.,pq(0,2),isig1,p43(0,2),
     $		 jqq(0,isig3,2,is3,0),gbkjqq(1,isig1,isig3,2,2,l),dummy)
            
            call ket2r(bkjqq(1,isig1,isig3,1,2,0),.false.,bq(0,1),isig1,q,
     $		 eps(0,l),bkjqqg(1,isig1,isig3,1,2,l),dummy)
            call bra2r(bkjqq(1,isig1,isig3,2,2,0),.false.,bq(0,2),isig1,q,
     $		 eps(0,l),bkjqqg(1,isig1,isig3,2,2,l),dummy)

c     ------------          
            
c     VV from different fermions (f1/f2):     
c     
c     eps1=v1,eps2=j43 with v2 arbitrary,eps3=v3, g from lower line: 
c     use m5(contr.=2,v1,v2,j=3,l=1:2)
            m5(2,3,1,3,l) = -s1c(psia(1,isig1,2),jqq(0,isig3,2,is3,l),
     $           .true.,isig1,psiwe(1,isig1,1))

            m5(2,1,3,3,l) = -s1c(psiwe(1,isig1,2),jqq(0,isig3,2,is3,l),
     $           .true.,isig1,psia(1,isig1,1))

c     eps1=v1,eps2=j43 with v2 arbitrary,eps3=v3, g from upper line: 
c     use m5(contr.=2,v1,v2,j=2,l=1:2)

c     WA
            ga = -s1c(braketgau(1,isig1,2,l),jqq(0,isig3,2,is3,0),.true.,
     $           isig1,psiwe(1,isig1,1))
            gb = -s1c(braketaug(1,isig1,2,l),jqq(0,isig3,2,is3,0),.true.,
     $           isig1,psiwe(1,isig1,1))
            gc = -s1c(psia(1,isig1,2),jqq(0,isig3,2,is3,0),.true.,
     $           isig1,braketweg(1,isig1,1,l))
            gd = -s1c(psia(1,isig1,2),jqq(0,isig3,2,is3,0),.true.,
     $           isig1,braketgwe(1,isig1,1,l))
            
            m5(2,3,1,2,l) = ga+gb+gc+gd


c     AW:
            ga = -s1c(braketgwe(1,isig1,2,l),jqq(0,isig3,2,is3,0),.true.,
     $           isig1,psia(1,isig1,1))
            gb = -s1c(braketweg(1,isig1,2,l),jqq(0,isig3,2,is3,0),.true.,
     $           isig1,psia(1,isig1,1))
            gc = -s1c(psiwe(1,isig1,2),jqq(0,isig3,2,is3,0),.true.,
     $           isig1,braketaug(1,isig1,1,l))
            gd = -s1c(psiwe(1,isig1,2),jqq(0,isig3,2,is3,0),.true.,
     $           isig1,braketgau(1,isig1,1,l))
            
            m5(2,1,3,2,l) = ga+gb+gc+gd


c---------------
            
c     VV from f2:
c     eps1=j43,eps2=v2,eps3=v3, g from lower line
            m5(1,3,1,3,l) = -s1c(psia(1,isig1,2),we,.true.,isig1,
     $           bkjqq(1,isig1,isig3,1,3,l))
            m5(1,1,3,3,l) = -s1c(psiwe(1,isig1,2),ea,.true.,isig1,
     $           bkjqq(1,isig1,isig3,1,3,l))
            
c     eps1=j43,eps2=v2,eps3=v3, g from upper line
c     
c     WA: 
            ga = -s1c(braketgau(1,isig1,2,l),we,.true.,isig1,
     $           bkjqq(1,isig1,isig3,1,2,0))
            gb = -s1c(braketaug(1,isig1,2,l),we,.true.,isig1,
     $           bkjqq(1,isig1,isig3,1,2,0))
            gc = -s1c(psia(1,isig1,2),we,.true.,isig1,bkjqqg(1,isig1,isig3,1,2,l))
            gd = -s1c(psia(1,isig1,2),we,.true.,isig1,gbkjqq(1,isig1,isig3,1,2,l))         
            
            m5(1,3,1,2,l) = ga+gb+gc+gd
            

c     AW: 
            ga = -s1c(braketgwe(1,isig1,2,l),ea,.true.,isig1,
     $           bkjqq(1,isig1,isig3,1,2,0))
            gb = -s1c(braketweg(1,isig1,2,l),ea,.true.,isig1,
     $           bkjqq(1,isig1,isig3,1,2,0))
            gc = -s1c(psiwe(1,isig1,2),ea,.true.,isig1,bkjqqg(1,isig1,isig3,1,2,l))
            gd = -s1c(psiwe(1,isig1,2),ea,.true.,isig1,gbkjqq(1,isig1,isig3,1,2,l))         
            
            m5(1,1,3,2,l) = ga+gb+gc+gd


c---------------------------
c     
c     VV from f1:
c     eps1=v1,eps2=v2,eps3=j43, g from lower line
            m5(3,3,1,3,l) = -s1c(bkjqq(1,isig1,isig3,2,3,l),ea,
     $           .true.,isig1,psiwe(1,isig1,1))

            m5(3,1,3,3,l) = -s1c(bkjqq(1,isig1,isig3,2,3,l),we,
     $           .true.,isig1,psia(1,isig1,1))
            
c     eps1=v1,eps2=v2,eps3=j43, g from upper line
c     
c     WA:
            ga = -s1c(gbkjqq(1,isig1,isig3,2,2,l),ea,.true.,isig1,psiwe(1,isig1,1))
            gb = -s1c(bkjqqg(1,isig1,isig3,2,2,l),ea,.true.,isig1,psiwe(1,isig1,1))
            gc = -s1c(bkjqq(1,isig1,isig3,2,2,0),ea,
     $           .true.,isig1,braketweg(1,isig1,1,l))
            gd = -s1c(bkjqq(1,isig1,isig3,2,2,0),ea,
     $           .true.,isig1,braketgwe(1,isig1,1,l))
            
            m5(3,3,1,2,l) = ga+gb+gc+gd
            
            
c     AW:
            ga = -s1c(gbkjqq(1,isig1,isig3,2,2,l),we,.true.,isig1,psia(1,isig1,1))
            gb = -s1c(bkjqqg(1,isig1,isig3,2,2,l),we,.true.,isig1,psia(1,isig1,1))
            gc = -s1c(bkjqq(1,isig1,isig3,2,2,0),we,
     $           .true.,isig1,braketaug(1,isig1,1,l))
            gd = -s1c(bkjqq(1,isig1,isig3,2,2,0),we,
     $           .true.,isig1,braketgau(1,isig1,1,l))
            
            m5(3,1,3,2,l) = ga+gb+gc+gd
            
            
c--------------------------------------------

            do j = jmin,jmax
               kb = 1
               
               if (isig3.eq.-1) then 
                  
                  k = 2         !1
                  
                  mat(k,isig1,isig3,j,l,7) = mat(k,isig1,isig3,j,l,7) 
     $                 + clr(3,3,-1)**3*prop43(3,j)*(
     $                 m5(2,3,kb,j,l)*clr(ifl(2,k),kb,isig1)
     $                 +m5(3,3,kb,j,l)*clr(7-ifl(1,k),kb,isig1)
     $                 +m5(3,kb,3,j,l)*clr(ifl(1,k),kb,isig1))
                  
                  k = 1         !2
                  
                  mat(k,isig1,isig3,j,l,7) = mat(k,isig1,isig3,j,l,7) 
     $                 + clr(3,3,-1)**3*prop43(3,j)*(
     $                 m5(1,3,kb,j,l)*clr(ifl(2,k),kb,isig1)
     $                 +m5(1,kb,3,j,l)*clr(7-ifl(2,k),kb,isig1)
     $                 +m5(2,kb,3,j,l)*clr(ifl(1,k),kb,isig1))
                  
               endif            !isig3
               
               
               do k= 3,4		      
                  do kp = 1,2
                     mat(k,isig1,isig3,j,l,7) = mat(k,isig1,isig3,j,l,7) + 
     $                    clr(3,3,-1)*prop43(kp,j)*clr(ifl(3,k),kp,isig3)*(
     $                    m5(1,3,kb,j,l)*clr(ifl(2,k),kb,isig1)*clr(ifl(1,k),kp,isig1)
     $                    +m5(1,kb,3,j,l)*clr(7-ifl(2,k),kb,isig1)*clr(ifl(1,k),kp,isig1)
     $                    +m5(2,3,kb,j,l)*clr(ifl(2,k),kb,isig1)*clr(ifl(2,k),kp,isig1)
     $                    +m5(2,kb,3,j,l)*clr(ifl(1,k),kb,isig1)*clr(7-ifl(2,k),kp,isig1)
     $                    +m5(3,3,kb,j,l)*clr(7-ifl(1,k),kb,isig1)*clr(ifl(2,k),kp,isig1)
     $                    +m5(3,kb,3,j,l)*clr(ifl(1,k),kb,isig1)*clr(ifl(2,k),kp,isig1))  
                  enddo         !kp
               enddo            !k
            enddo               !j
         enddo                  !l
      enddo                     !isig3
      
c     ------------------------------
c     
c     same for lower line:
c     
c     need m5(pos.j21=1:3,vext1=1:3,vext2=1:3), id of internal boson arbitrary 

      do isig1 = -1,1,2 
	 isig3 = -1

!     g from lower line, but not in bkjqq:
	 call ket2c(psi(1,isig3,3),.true.,p(0,3),isig3,p21(0,3),
     $        jqq(0,isig1,1,is1,0),bkjqq(1,isig3,isig1,3,3,0),bq(0,3))
         call bra2c(psi(1,isig3,4),.true.,p(0,4),isig3,p21(0,3),
     $        jqq(0,isig1,1,is1,0),bkjqq(1,isig3,isig1,4,3,0),bq(0,4))

         
         do l = 1,2
            
!     g from upper line:
            call ket2c(psi(1,isig3,3),.true.,p(0,3),isig3,p21(0,2),
     1           jqq(0,isig1,1,is1,l),bkjqq(1,isig3,isig1,3,2,l),dummy)
            call bra2c(psi(1,isig3,4),.true.,p(0,4),isig3,p21(0,2),
     1           jqq(0,isig1,1,is1,l),bkjqq(1,isig3,isig1,4,2,l),dummy)
            
!     g from lower line, in gbkjqq or bkjqqg:
            call ket2c(braketg(1,isig3,3,l),.false.,pq(0,3),isig3,p21(0,3),
     $		 jqq(0,isig1,1,is1,0),gbkjqq(1,isig3,isig1,3,3,l),dummy)
            call bra2c(braketg(1,isig3,4,l),.false.,pq(0,4),isig3,p21(0,3),
     $		 jqq(0,isig1,1,is1,0),gbkjqq(1,isig3,isig1,4,3,l),dummy)
            
            call ket2r(bkjqq(1,isig3,isig1,3,3,0),.false.,bq(0,3),isig3,q,
     $		 eps(0,l),bkjqqg(1,isig3,isig1,3,3,l),dummy)
            call bra2r(bkjqq(1,isig3,isig1,4,3,0),.false.,bq(0,4),isig3,q,
     $		 eps(0,l),bkjqqg(1,isig3,isig1,4,3,l),dummy)

c     ------------          
            
c     VV from different fermions (f3/f4):     
c     
c     
c     eps1=v1,eps2=j21 with v2 arbitrary,eps3=v3, g from upper line: 
c     use m5(contr.=2,v1,v2,j=3,l=1:2)
            m5(2,3,1,2,l) = -s1c(psia(1,isig3,4),jqq(0,isig1,1,is1,l),
     $           .true.,isig3,psiwe(1,isig3,3))

            m5(2,1,3,2,l) = -s1c(psiwe(1,isig3,4),jqq(0,isig1,1,is1,l),
     $           .true.,isig3,psia(1,isig3,3))

c     eps1=v1,eps2=j21 with v2 arbitrary,eps3=v3, g from lower line: 
c     use m5(contr.=2,v1,v2,j=2,l=1:2)
c     
c     WA
            ga = -s1c(braketgau(1,isig3,4,l),jqq(0,isig1,1,is1,0),.true.,
     $           isig3,psiwe(1,isig3,3))
            gb = -s1c(braketaug(1,isig3,4,l),jqq(0,isig1,1,is1,0),.true.,
     $           isig3,psiwe(1,isig3,3))
            gc = -s1c(psia(1,isig3,4),jqq(0,isig1,1,is1,0),.true.,
     $           isig3,braketweg(1,isig3,3,l))
            gd = -s1c(psia(1,isig3,4),jqq(0,isig1,1,is1,0),.true.,
     $           isig3,braketgwe(1,isig3,3,l))
            
            m5(2,3,1,3,l) = ga+gb+gc+gd


c     AW:
            ga = -s1c(braketgwe(1,isig3,4,l),jqq(0,isig1,1,is1,0),.true.,
     $           isig3,psia(1,isig3,3))
            gb = -s1c(braketweg(1,isig3,4,l),jqq(0,isig1,1,is1,0),.true.,
     $           isig3,psia(1,isig3,3))
            gc = -s1c(psiwe(1,isig3,4),jqq(0,isig1,1,is1,0),.true.,
     $           isig3,braketaug(1,isig3,3,l))
            gd = -s1c(psiwe(1,isig3,4),jqq(0,isig1,1,is1,0),.true.,
     $           isig3,braketgau(1,isig3,3,l))
            
            m5(2,1,3,3,l) = ga+gb+gc+gd


c---------------
c     
c     VV from f4:
c     eps1=j43,eps2=v2,eps3=v3, g from upper line
            m5(1,3,1,2,l) = -s1c(psia(1,isig3,4),we,.true.,isig3,
     $           bkjqq(1,isig3,isig1,3,2,l))

            m5(1,1,3,2,l) = -s1c(psiwe(1,isig3,4),ea,.true.,isig3,
     $           bkjqq(1,isig3,isig1,3,2,l))
            
c     eps1=j43,eps2=v2,eps3=v3, g from lower line
c     
c     WA: 
            ga = -s1c(braketgau(1,isig3,4,l),we,.true.,isig3,
     $           bkjqq(1,isig3,isig1,3,3,0))
            gb = -s1c(braketaug(1,isig3,4,l),we,.true.,isig3,
     $           bkjqq(1,isig3,isig1,3,3,0))
            gc = -s1c(psia(1,isig3,4),we,.true.,isig3,bkjqqg(1,isig3,isig1,3,3,l))
            gd = -s1c(psia(1,isig3,4),we,.true.,isig3,gbkjqq(1,isig3,isig1,3,3,l))       
            
            m5(1,3,1,3,l) = ga+gb+gc+gd
            

c     AW: 
            ga = -s1c(braketgwe(1,isig3,4,l),ea,.true.,isig3,
     $           bkjqq(1,isig3,isig1,3,3,0))
            gb = -s1c(braketweg(1,isig3,4,l),ea,.true.,isig3,
     $           bkjqq(1,isig3,isig1,3,3,0))
            gc = -s1c(psiwe(1,isig3,4),ea,.true.,isig3,bkjqqg(1,isig3,isig1,3,3,l))
            gd = -s1c(psiwe(1,isig3,4),ea,.true.,isig3,gbkjqq(1,isig3,isig1,3,3,l))     
            
            m5(1,1,3,3,l) = ga+gb+gc+gd


c---------------------------

c     VV from f3:
c     eps1=v1,eps2=v2,eps3=j21, g from upper line
            m5(3,3,1,2,l) = -s1c(bkjqq(1,isig3,isig1,4,2,l),ea,
     $           .true.,isig3,psiwe(1,isig3,3))

            m5(3,1,3,2,l) = -s1c(bkjqq(1,isig3,isig1,4,2,l),we,
     $           .true.,isig3,psia(1,isig3,3))
            
c     eps1=v1,eps2=v2,eps3=j21, g from lower line
c     
c     WA:
            ga = -s1c(gbkjqq(1,isig3,isig1,4,3,l),ea,.true.,isig3,psiwe(1,isig3,3))
            gb = -s1c(bkjqqg(1,isig3,isig1,4,3,l),ea,.true.,isig3,psiwe(1,isig3,3))
            gc = -s1c(bkjqq(1,isig3,isig1,4,3,0),ea,
     $           .true.,isig3,braketweg(1,isig3,3,l))
            gd = -s1c(bkjqq(1,isig3,isig1,4,3,0),ea,
     $           .true.,isig3,braketgwe(1,isig3,3,l))
            
            m5(3,3,1,3,l) = ga+gb+gc+gd
            
            
c     AW:
            ga = -s1c(gbkjqq(1,isig3,isig1,4,3,l),we,.true.,isig3,psia(1,isig3,3))
            gb = -s1c(bkjqqg(1,isig3,isig1,4,3,l),we,.true.,isig3,psia(1,isig3,3))
            gc = -s1c(bkjqq(1,isig3,isig1,4,3,0),we,
     $           .true.,isig3,braketaug(1,isig3,3,l))
            gd = -s1c(bkjqq(1,isig3,isig1,4,3,0),we,
     $           .true.,isig3,braketgau(1,isig3,3,l))
            
            m5(3,1,3,3,l) = ga+gb+gc+gd
            
            
c--------------------------------------------

            do j = jmin,jmax
               kb = 1	
               if (isig1.eq.-1) then 
                  
                  k = 4
                  
                  mat(k,isig1,isig3,j,l,8) = mat(k,isig1,isig3,j,l,8)
     $                 + clr(3,3,-1)**3*prop21(3,j)*(
     $                 m5(2,3,kb,j,l)*clr(ifl(4,k),kb,isig3)
     $                 +m5(3,3,kb,j,l)*clr(7-ifl(3,k),kb,isig3)
     $                 +m5(3,kb,3,j,l)*clr(ifl(3,k),kb,isig3))
                  
                  k = 3
                  
                  mat(k,isig1,isig3,j,l,8) = mat(k,isig1,isig3,j,l,8)
     $                 + clr(3,3,-1)**3*prop21(3,j)*(
     $                 m5(1,3,kb,j,l)*clr(ifl(4,k),kb,isig3)
     $                 +m5(1,kb,3,j,l)*clr(7-ifl(4,k),kb,isig3)
     $                 +m5(2,kb,3,j,l)*clr(ifl(3,k),kb,isig3))
                  
               endif            !isig1
               
               do k= 1,2	
                  do kp = 1,2
                     
                     mat(k,isig1,isig3,j,l,8) = mat(k,isig1,isig3,j,l,8) + 
     $                    clr(3,3,-1)*prop21(kp,j)*clr(ifl(1,k),kp,isig1)*(
     $                    m5(1,3,kb,j,l)*clr(ifl(4,k),kb,isig3)*clr(ifl(3,k),kp,isig3)
     $                    +m5(1,kb,3,j,l)*clr(7-ifl(4,k),kb,isig3)*clr(ifl(3,k),kp,isig3)
     $                    +m5(2,3,kb,j,l)*clr(ifl(4,k),kb,isig3)*clr(ifl(4,k),kp,isig3)
     $                    +m5(2,kb,3,j,l)*clr(ifl(3,k),kb,isig3)*clr(7-ifl(4,k),kp,isig3)
     $                    +m5(3,3,kb,j,l)*clr(7-ifl(3,k),kb,isig3)*clr(ifl(4,k),kp,isig3)
     $                    +m5(3,kb,3,j,l)*clr(ifl(3,k),kb,isig3)*clr(ifl(4,k),kp,isig3))

                  enddo         !kp
               enddo            !k
            enddo               !j
         enddo                  !l
      enddo                     !isig1


c ------------------------------------------------------------------------------
c     ------------------------------------------------------------------------------
      
c     sum the graphs, square them and map them onto ansi etc.

c     i = 1		VV
c     i = 2,3	BV
c     i = 4,5	WZZ
c     i = 6,9	BB
c     i = 7,8	P
      
      do k = 1,4
 	 do j = 2,3     
 	    res(k,j) = 0
	    do isig1 = -1,1,2
	       do isig3 = -1,1,2
 	    	  do l = 1,2
              	     mm(k,isig1,isig3,j,l) = 0
               	     do i = 1,8
                        mm(k,isig1,isig3,j,l) = 
     1                       mm(k,isig1,isig3,j,l) + 
     1                       (mat(k,isig1,isig3,j,l,i))
      		     enddo      !i
              	     res(k,j) = res(k,j) 
     &                    + dreal(mm(k,isig1,isig3,j,l))**2
     &                    + dimag(mm(k,isig1,isig3,j,l))**2
	          enddo         !l
               enddo            !isig3		     
            enddo               !isig1
            res(k,j) = res(k,j)*12d0*fpials(j) ! C_2*9 is the color factor
	 enddo                  !j
      enddo                     !k      

#ifdef WITH_MADGRAPH      
      if (lmgdebug) then
         if(gsign.eq.-1) then
            call mgr_uuscgi(pbar(0,1),sign(1),qbar,gsign,jmin,ans1(1))
            call mgr_ddscgi(pbar(0,1),sign(1),qbar,gsign,jmin,ans2(1))
            call mgr_duccgi(pbar(0,1),sign(1),qbar,gsign,jmin,ans3(1))
            call mgr_dussgi(pbar(0,1),sign(1),qbar,gsign,jmin,ans4(1))

 	    print*," gsign = ",gsign
            do j = jmin,jmax
               print*," mg, mat, mat/mg for j ",j

               print*,"ddcsg",ans1(1),(res(1,j)/fpials(j))/3/8*2,
     &              (res(1,j)/fpials(j))/ans1(1)/3/8*2, !  /3/8 for color
     &              abs(1-abs((res(1,j)/fpials(j))/3/8*2/ans1(1)))
               print*,"ddcsg",ans2(1),(res(2,j)/fpials(j))/3/8*2,
     &              (res(2,j)/fpials(j))/ans2(1)/3/8*2, !  /3/8 for color
     &              abs(1-abs((res(2,j)/fpials(j))/3/8*2/ans2(1)))
               print*,"duccg",ans3(1),(res(3,j)/fpials(j))/3/8*2,
     &              (res(3,j)/fpials(j))/ans3(1)/3/8*2, !  /3/8 for color
     &              abs(1-abs((res(3,j)/fpials(j))/3/8*2/ans3(1)))
               print*,"dussg",ans4(1),(res(4,j)/fpials(j))/3/8*2,
     &              (res(4,j)/fpials(j))/ans4(1)/3/8*2, !  /3/8 for color
     &              abs(1-abs((res(4,j)/fpials(j))/3/8*2/ans4(1)))
               print*
            enddo
         else
            call mgr_uuscg(pbar(0,1),sign(1),gsign,qbar,ans1(1))
            call mgr_ddscg(pbar(0,1),sign(1),gsign,qbar,ans2(1))
            call mgr_duccg(pbar(0,1),sign(1),gsign,qbar,ans3(1))
            call mgr_dussg(pbar(0,1),sign(1),gsign,qbar,ans4(1))

            print*," gsign = ",gsign
            print*," mg, mat, mat/mg for sum j = 2,3:"
            print*,"uuscg",ans1(1),(res(1,2)/fpials(2)+res(1,3)/fpials(3))/9*2,
     &           (res(1,2)/fpials(2)+res(1,3)/fpials(3))/ans1(1)/9*2, !  /9 for color *2 for spin
     &           abs(1-abs((res(1,2)/fpials(2)+res(1,3)/fpials(3))/9*2/ans1(1)))
            print*,"ddscg",ans2(1),(res(2,2)/fpials(2)+res(2,3)/fpials(3))/9*2,
     &           (res(2,2)/fpials(2)+res(2,3)/fpials(3))/ans2(1)/9*2,  
     &           abs(1-abs((res(2,2)/fpials(2)+res(2,3)/fpials(3))/9*2/ans2(1)))
            print*,"duccg",ans3(1),(res(3,2)/fpials(2)+res(3,3)/fpials(3))/9*2,
     &           (res(3,2)/fpials(2)+res(3,3)/fpials(3))/ans3(1)/9*2,  
     &           abs(1-abs((res(3,2)/fpials(2)+res(3,3)/fpials(3))/9*2/ans3(1)))
            print*,"dussg",ans4(1),(res(4,2)/fpials(2)+res(4,3)/fpials(3))/9*2,
     &           (res(4,2)/fpials(2)+res(4,3)/fpials(3))/ans4(1)/9*2,  
     &           abs(1-abs((res(4,2)/fpials(2)+res(4,3)/fpials(3))/9*2/ans4(1)))
            print*
         endif
         write(*,*)
         read(*,*)
      endif                     !ldebug
#endif      

      if (jmin.eq.3) then
         ans1(2) = 0d0
         ans2(2) = 0d0
         ans3(2) = 0d0
         ans4(2) = 0d0
      elseif (jmax.eq.2) then
         ans1(3) = 0d0
         ans2(3) = 0d0
         ans3(3) = 0d0
         ans4(3) = 0d0
      endif


c     set processes to zero if photon virtuality falls below cutoff
      if ( abs(p21(4,3)).lt.qsqAmin ) then
         do k = 1,4
            res(k,3) = 0        ! no graphs with gluon on lower line
         enddo
      endif
      if ( abs(p43(4,2)).lt.qsqAmin ) then
         do k = 1,4
            res(k,2) = 0        ! no graphs with gluon on upper line
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
         ans1(j) = res(1,j) 
         ans2(j) = res(2,j)
         ans3(j) = res(3,j)
         ans4(j) = res(4,j)
      enddo
      ans1(1) = res(1,3) +  res(1,2)
      ans2(1) = res(2,3) +  res(2,2)
      ans3(1) = res(3,3) +  res(3,2)
      ans4(1) = res(4,3) +  res(4,2)

      return
      end

