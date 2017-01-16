      subroutine qqsp2wwqqj(pbar,sign,qbar,gsign,
     1                  uucc,uuss,ddcc,ddss,udsc,ducs)
      implicit none
c
c  modified from qqwwqqj.F by Jessica Frank
C
C  QQSP2WWQQJ calculates the matrix elements**2 for electroweak
c  weak boson pair production via a spin-2 resonance by quark quark scattering
C
C        q1 q3    ---->   q2 q4 g W+W-,   W+ ---> f5-bar f6, W- ---> f7-bar f8
c
c Only resonant spin-2 contributions are considered!
c
C  This code includes only real emission contributions, i.e.
c
c      return uucc = |M_real|^2   etc.
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
#include "tensor.inc"
#include "tensorsp2aa.inc"
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

      double precision  pbar(0:3,4+max_v),qbar(0:4)
      double precision uucc(3),uuss(3),ddcc(3),ddss(3),udsc(3),ducs(3),
     1 		       res(6,2:3),resv(6)
      double precision  p(0:3,4+max_v),q(0:4), p21(0:4,2:3), p43(0:4,2:3),
     1                  pq(0:4,4),pwp(0:4),pwm(0:4),pww(0:4),
     2			qpm(0:4,2:3),qmp(0:4,2:3)!,
c     2			qpml(0:4,2:3),qmpl(0:4,2:3)
      integer  sign(4+max_v),gsign, mu, i, j, k,
     1         isig, isig1, isig3
      integer  ifl(4,6) !, is1, is3 taken from "tensor.inc"
      integer  l   ! gluon polariz. (l=0:no g, l=1,2:g in kartesian basis)
      integer jmin, jmax
      logical jlog2,jlog3
      double complex prop21(4,2:3), prop43(4,2:3),
     1		     prop_pm(4,2:3),prop_mp(4,2:3)
      double complex mat(6,-1:1,-1:1,2:3,0:2,9)
      double complex mata(6,-1:1,-1:1,2:3,0:2,9),
     1 		     matz(6,-1:1,-1:1,2:3,0:2,9)
      double complex mm(6,-1:1,-1:1,2:3,2)

      double complex maa, maz, mza, mzz, mww5, mww6


c      double complex  m1,m2 ! for checks only
      double complex  raz! for checks only
      double precision eps(0:3,2) ! g in kartesian basis 
      double complex psi(2,-1:1,4),jqq(0:5,-1:1,2,-1:1,0:2), 
     1 		     braketg(2,-1:1,4,2), jh1(0:5,-1:1), jh2(0:5,-1:1)

      double complex  zm2i(2:4)
      double complex contract_Tjj,! contract_T1j,contract_T2j,
     1 		     dotcc, dotrc, s1c !, dotqj
      external contract_Tjj, !contract_T1j,contract_T2j,
     1 	       dotcc, dotrc, s1c !, dotqj

      double complex im
      parameter (im = (0d0,1d0))
      logical ldebug, linit,lerror,ldebugm
     
      save ifl, zm2i, linit, lerror,ldebug,ldebugm
c stuff for comparison with madgraph :
#include "tampr-lists.inc" 
	
      data linit /.true./, lerror /.true./, ldebug /.false./, ldebugm /.false./

      real*8 widthfactor

      widthfactor=1.0d0
!       if (ofs.and.with_anom) then
!          widthfactor=0.d0
!       endif
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
         print*," enter qqsp2wwqqj "
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
	   
	   maa = contract_Tjj(sp2aaww(0,0,j),jqq(0,isig1,1,is1,l),
     #	   					jqq(0,isig3,2,is3,0))
           maz = contract_Tjj(sp2azww(0,0,j),jqq(0,isig1,1,is1,l),
     #	   					jqq(0,isig3,2,is3,0))
           mza = contract_Tjj(sp2zaww(0,0,j),jqq(0,isig3,2,is3,0),
     #	   					jqq(0,isig1,1,is1,l))
           mzz = contract_Tjj(sp2zzww(0,0,j),jqq(0,isig1,1,is1,l),
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
	   
	   maa = contract_Tjj(sp2aaww(0,0,j),jqq(0,isig1,1,is1,0),
     #	   					jqq(0,isig3,2,is3,l))
           maz = contract_Tjj(sp2azww(0,0,j),jqq(0,isig1,1,is1,0),
     #	   					jqq(0,isig3,2,is3,l))
           mza = contract_Tjj(sp2zaww(0,0,j),jqq(0,isig3,2,is3,l),
     #	   					jqq(0,isig1,1,is1,0))
           mzz = contract_Tjj(sp2zzww(0,0,j),jqq(0,isig1,1,is1,0),
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
            
      mww6 = contract_Tjj(sp2wwww6(0,0,j),jqq(0,-1,1,is1,l),
     #      					jqq(0,-1,2,is3,0))
      mww5 = contract_Tjj(sp2wwww5(0,0,j),jqq(0,-1,2,is3,0),
     #						jqq(0,-1,1,is1,l))
      mat(5,-1,-1,j,l,1) = mww5*clr(3,3,-1)**2!*sqrt(fpials(j))
      mat(6,-1,-1,j,l,1) = mww6*clr(3,3,-1)**2
      

      endif
  
      if (jlog3) then
      j =  3 ! g from lower line
            
      mww6 = contract_Tjj(sp2wwww6(0,0,j),jqq(0,-1,1,is1,0),
     #						jqq(0,-1,2,is3,l))
      mww5 = contract_Tjj(sp2wwww5(0,0,j),jqq(0,-1,2,is3,l),
     #						jqq(0,-1,1,is1,0))
      mat(5,-1,-1,j,l,1) = mww5*clr(3,3,-1)**2!*sqrt(fpials(j))
      mat(6,-1,-1,j,l,1) = mww6*clr(3,3,-1)**2
      
     
      endif
      enddo ! l-loop
     

c sum the graphs, square them and map them onto uucc, uuss etc.
      do k = 1,4
 	 do j = 2,3     
 	    res(k,j) = 0
	    do isig1 = -1,1,2
	       do isig3 = -1,1,2
 	    	  do l = 1,2
              	     mm(k,isig1,isig3,j,l) = 0
               	        i = 1
                        mm(k,isig1,isig3,j,l) = 
     1                 		 mm(k,isig1,isig3,j,l) + 
     1		    	       (mat(k,isig1,isig3,j,l,i))
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
      if ( abs(p21(4,3)).lt.qsqAmin ) then
         do k = 1,4
            res(k,3) = 0      ! no graphs with gluon on lower line
         enddo
      endif
      if ( abs(p43(4,2)).lt.qsqAmin ) then
         do k = 1,4
            res(k,2) = 0      ! no graphs with gluon on upper line
         enddo
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

