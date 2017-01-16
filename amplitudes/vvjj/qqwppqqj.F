c
      subroutine qqwppqqj(pbar,sign,qbar,gsign,bos,ucds)
          use globalvars, only: ldoblha

      implicit none
c
c	Barbara Jaeger, <bjaeger@physik.uni-wuerzburg.de>
c	Initial version:  Oct. 2008 
C
C  QQWPPQQJ calculates the matrix elements**2 for electroweak
c  weak boson pair production by quark quark scattering
C
C        q1 q3    ---->   q2 q4 W+ W+ g,   W ---> f5 f6, W ---> f7 f8
C
c  where the final state leptons are 2 charged leptons + 2 neutrinos. 
C  Crossing related processes can be computed as well. Pauli interference 
c  terms for identical fermions are neglected. In particular, only the
c  t-channel exchange of elctroweak bosons is considered. s-channel
c  production of 3 weak bosons is NOT implemented.
c
c
C  This code includes only real emission contributions, i.e.
c
c      return ucds = |M_real|^2   etc.
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
c	isig1/isig3 are the helicities of parton1/3 
c
c
c---------------------------------------------------------------------
c
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"
#include "tensor.inc"
#include "tensorz.inc"
c
c electroweak couplings are taken from KOPPLN
c
      double precision  clr, xm2, xmg, b, v, a
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     1     V(4,5),A(4,5)
c
c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"
c      
c for comparison with madgraph:      
#include "tampr_wpp_list.inc"
c
c variables for the main part of the program
c
      real*8 fpials(2:3), fpi
      parameter (fpi=4d0*pi)

      double precision  pbar(0:3,4+n_v),qbar(0:4)
      double precision ucds(3),res(2:3)
      double precision  p(0:3,4+n_v),q(0:4), p21(0:4,2:3), p43(0:4,2:3),
     1                  pq(0:4,4),pew(0:4),puw(0:4),pww(0:4),
     2			qee(0:4,2:3),quu(0:4,2:3)
 
      integer  sign(4+n_v),gsign, bos
      integer  mu, i, j, isig, is
      integer  ifl(4) 
      integer  l   ! gluon polariz. (l=0:no g, l=1,2:g in kartesian basis)
      integer jmin, jmax
      logical jlog2,jlog3
      integer ide,idu
      integer lk,kk
      integer lh
 
      double complex prop21(4,2:3), prop43(4,2:3),
     1		     prop_ee(4,2:3),prop_uu(4,2:3)

      double complex mat(2:3,0:2,9),temp
      double complex mm(2:3,2)
      double complex maeu(2:3,2),mzeu(2:3,2),maue(2:3,2),mzue(2:3,2)
      double complex mb(4,4,2:3,2)
      double complex m1we(2:3,2),m2we(2:3,2),m3we(2:3,2),m4we(2:3,2),  
     $  	     m1ze(2:3,2),m2ze(2:3,2),m3ze(2:3,2),m4ze(2:3,2),   
     $  	     m1ae(2:3,2),m2ae(2:3,2),m3ae(2:3,2),m4ae(2:3,2),  
     $  	     m1e(2:3,2),m2e(2:3,2),m3e(2:3,2),m4e(2:3,2),  
     $  	     w2we(2:3,2),w4we(2:3,2)
      double complex m1wu(2:3,2),m2wu(2:3,2),m3wu(2:3,2),m4wu(2:3,2),  
     $  	     m1zu(2:3,2),m2zu(2:3,2),m3zu(2:3,2),m4zu(2:3,2),   
     $  	     m1au(2:3,2),m2au(2:3,2),m3au(2:3,2),m4au(2:3,2),  
     $  	     m1u(2:3,2),m2u(2:3,2),m3u(2:3,2),m4u(2:3,2),  
     $  	     w2wu(2:3,2),w4wu(2:3,2)     
      double complex m5eu(2:3,2),m5ue(2:3,2)
      double complex mww
        
      double complex epsawe(0:3,2,2:3,0:2),epszwe(0:3,2,2:3,0:2),
     $ 		     epsawu(0:3,2,2:3,0:2),epszwu(0:3,2,2:3,0:2),
     $ 		     epswae(0:3,2,2:3,0:2),epswze(0:3,2,2:3,0:2),
     $ 		     epswau(0:3,2,2:3,0:2),epswzu(0:3,2,2:3,0:2)
      double complex qepszwe(2,2:3,0:2),qepszwu(2,2:3,0:2)


      double precision eps(0:3,2) ! g in kartesian basis 
      double complex psi(2,-1:1,4),jqq(0:5,-1:1,2,-1:1,0:2), 
     1 		     braketg(2,-1:1,4,2), jh1(0:5,-1:1), jh2(0:5,-1:1)
      double complex psiwe(2,4), psiwu(2,4),
     1               jwe(0:5,-1:1,-1:1,4,0:2), jwu(0:5,-1:1,-1:1,4,0:2)
      double precision fqe(0:4,4), fqu(0:4,4)
      double complex braketgwe(2,-1:1,4,0:2),braketweg(2,-1:1,4,0:2),
     3		     braketgwu(2,-1:1,4,0:2),braketwug(2,-1:1,4,0:2)
      double precision  pgwe(0:4,4),pgwu(0:4,4),
     3			pewg(0:4,4),puwg(0:4,4)
      double complex jwgi(0:5,-1:1),jgwi(0:5,-1:1),jwg0i(0:5,-1:1),
     1		     jwgii(0:5,-1:1),jgwii(0:5,-1:1),jwg0ii(0:5,-1:1)
      double complex jj21u(0:2),jj21e(0:2),jj43u(0:2),jj43e(0:2)

      double complex zeu(2:3,2),zue(2:3,2)
      double complex  zm2i(2:4)
      double complex contract_Tjj,dotcc, dotrc, s1c
      external contract_Tjj, dotcc, dotrc, s1c 

      double complex m1hj(-1:1,2:3),m1kj(2,2:3)
      double complex m1haj(-1:1,2:3),m1hzj(-1:1,2:3),
     1               m1kaj(2,2:3),m1kzj(2,2:3)
      double complex im
      parameter (im = (0d0,1d0))
      logical ldebug, linit,lerror
      integer n,m
     
      save ifl, zm2i, linit, lerror,ldebug
	
      data linit /.true./, lerror /.true./, ldebug /.false./
c      data linit /.true./, lerror /.true./, ldebug /.true./
      
      logical vvdebug,bbdebug,vbdebug,pdebug
      parameter (vvdebug = .false.,bbdebug = .false.,
     &  	 pdebug = .false.,vbdebug = .false.) 
c      parameter (vvdebug = .false.,vbdebug = .false.,
c     &  	 pdebug = .true.,bbdebug = .false.) 
c      parameter (vvdebug = .true.,vbdebug = .true.,
c     &  	 pdebug = .true.,bbdebug = .true.) 
c
c  ---------------------------------------------------------------------
c 
c initialize & precompute stuff needed below:
c
c fix strong coupling gs**2 for the 2 quark lines:

      fpials(2) = fpi*als(1,1)
      fpials(3) = fpi*als(2,1)

*      write(*,*)'fpials =', fpials

c$$$      write(*,*)'fpi, als(1,1) =', fpi, als(1,1)
c$$$      write(*,*)'fpi, als(2,1) =', fpi, als(2,1)
c$$$      stop
	  
c define flavors of external quarks for the 4 NC and 2 CC subprocesses
c
      if (linit) then           ! u(1)c(3)->d(2)s(4) etc.
         ifl(1) = 3
         ifl(2) = 4
         ifl(3) = 3
         ifl(4) = 4
c	 
         zm2i(2) = 1/dcmplx(xm2(2),-xmg(2))
         zm2i(3) = 1/dcmplx(xm2(3),-xmg(3))
         zm2i(4) = 1/dcmplx(xm2(4),-xmg(4))

         linit = .false.	 
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
      endif	! end gluon = final/initial
      
      do j = 2,3
         do l = 0,2
            do i = 1,9
               mat(j,l,i)  = 0
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
		
c (is1,is3 are fixed here and don't change throughout this run of the program)
c
c
c define the internal momenta
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
	 
         pew(mu) =   p(mu,6) - p(mu,5)
         puw(mu) =   p(mu,8) - p(mu,7)
         pww(mu) =   pew(mu) + puw(mu)
	 
      enddo
	
      do j= 2,3
         p21(4,j) = p21(0,j)**2 - p21(1,j)**2 - p21(2,j)**2 - p21(3,j)**2
         p43(4,j) = p43(0,j)**2 - p43(1,j)**2 - p43(2,j)**2 - p43(3,j)**2     
      enddo
	
      q(4)   = 0d0
      pew(4) = pew(0)**2 - pew(1)**2 - pew(2)**2 - pew(3)**2
      puw(4) = puw(0)**2 - puw(1)**2 - puw(2)**2 - puw(3)**2
      pww(4) = pww(0)**2 - pww(1)**2 - pww(2)**2 - pww(3)**2

c  ---------------------------------------------------------------------
c
c get the vector boson propagator factors
c
c depending on value of j, gluon is attached to respective quark line or not;
c no V is attached here 

      do j = 2,3	
         prop21(1,j) = 1/p21(4,j) !prop21(bos=1:4,j=2:3)
c         if (p21(4,j).le.0d0) then
c            prop21(2,j) = 1/(p21(4,j)-xm2(2))
c            prop21(3,j) = 1/(p21(4,j)-xm2(3))
c         else
            prop21(2,j) = 1/dcmplx(p21(4,j)-xm2(2),xmg(2))
            prop21(3,j) = 1/dcmplx(p21(4,j)-xm2(3),xmg(3))
c         endif
         prop21(4,j) = prop21(3,j)

         prop43(1,j) = 1/p43(4,j)
c         if (p43(4,j).le.0d0) then
c            prop43(2,j) = 1/(p43(4,j)-xm2(2))
c            prop43(3,j) = 1/(p43(4,j)-xm2(3))
c         else
            prop43(2,j) = 1/dcmplx(p43(4,j)-xm2(2),xmg(2))
            prop43(3,j) = 1/dcmplx(p43(4,j)-xm2(3),xmg(3))
c         endif
         prop43(4,j) = prop43(3,j)

      enddo
c
c  ---------------------------------------------------------------------

c for box-box and BV graphs we need the propagators for t-channel bosons between quark lines
c as seen from upper line these W momenta are INCOMING. They are OUTGOINg as seen from lower 
c line
      
      do j = 2,3                ! g from upper/lower line		
      
         do mu = 0,3
            qee(mu,j) = pew(mu)+p21(mu,j) !We emitted on upper line
            quu(mu,j) = puw(mu)+p21(mu,j) !Wu emitted on upper line
         enddo
      
         qee(4,j) = qee(0,j)**2-qee(1,j)**2-qee(2,j)**2-qee(3,j)**2
         quu(4,j) = quu(0,j)**2-quu(1,j)**2-quu(2,j)**2-quu(3,j)**2

         prop_ee(1,j) = clr(3,3,-1)**2/qee(4,j)
         prop_ee(2,j) = clr(3,3,-1)**2/dcmplx(qee(4,j)-xm2(2),xmg(2))
         prop_ee(3,j) = clr(3,3,-1)**4/dcmplx(qee(4,j)-xm2(3),xmg(3))
         prop_uu(1,j) = clr(3,3,-1)**2/quu(4,j)
         prop_uu(2,j) = clr(3,3,-1)**2/dcmplx(quu(4,j)-xm2(2),xmg(2))
         prop_uu(3,j) = clr(3,3,-1)**4/dcmplx(quu(4,j)-xm2(3),xmg(3))
         
         prop_ee(4,j) = prop_ee(3,j)
         prop_uu(4,j) = prop_uu(3,j)
         
      enddo  ! j loop


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
        call curr6(-1,psi(1,-1,2),p(0,2),psi(1,-1,1),p(0,1),
     -     jqq(0,-1,1,is1,0))      
        call curr6(-1,psi(1,-1,4),p(0,4),psi(1,-1,3),p(0,3),
     -       jqq(0,-1,2,is3,0))
c
c nomenclature: jqq(mu,		...	Lorentz index (0:3),(4:5) moment. info
c		    hel,	...	quark helicity (+-1)
c		    u/l,	...	upper(1)/lower(2) quark line
c		    is,		...	current for qq(g),qbqb(g),or gqbq line
c		     l)		...	g polarization (l=0:no g,l=1,2 with g)
c
c
c  Get the gluon polarization vector and the gluon emission spinors
        do l = 1,2              ! 2 gluon polarizations
           call polvec(qbar,l,eps(0,l)) ! get gluon pol.vectors
	 	 
           isig = -1            ! fermion helicity 
 
c     NOTES for bras and kets: .true. if psi is a 2-spinor of the chi
c     form as output by psi0m, .false. otherwise.  the last entry is
c     the sum of the two momenta (p plus q) and effectively the
c     momentum of the new spinor.
            
           call ket2r(psi(1,isig,1),.true.,p(0,1),isig,q,eps(0,l),
     $          braketg(1,isig,1,l),pq(0,1)) ! |q,1>_l,isig1
           call bra2r(psi(1,isig,2),.true.,p(0,2),isig,q,eps(0,l),
     $          braketg(1,isig,2,l),pq(0,2)) ! <2,q|_l,isig2
           call ket2r(psi(1,isig,3),.true.,p(0,3),isig,q,eps(0,l),
     $          braketg(1,isig,3,l),pq(0,3)) ! |q,3>_l,isig3
           call bra2r(psi(1,isig,4),.true.,p(0,4),isig,q,eps(0,l),
     $          braketg(1,isig,4,l),pq(0,4)) ! <4,q|_l,isig4

	          
c     braketg contains the free quark spinors multiplied by a fermion
c     propagator and a gluon eps_slash. 
c     NOTATION: braketg(2 component spinor, isig =-1 or 1 (fermion hel.),
c     fermion ID = 1:4, gluon polarization l=1:2)
 
        enddo  ! l loop
       
c     Get the f-fbar currents with one gluon radiated from the
c     current line.  There are two terms, one for gluon emission to
c     either side of the ffV vertex:
c
c	gluon from upper line:
        do l = 1, 2             ! gluon polarizations
           call curr6(-1,psi(1,-1,2),p(0,2),braketg(1,-1,1,l),pq(0,1),jh1)	
c                                            =  <2|vertex|q,1>_l,isig1
           call curr6(-1,braketg(1,-1,2,l),pq(0,2),psi(1,-1,1),p(0,1),jh2)	
c                                            =  <2,q|vertex|1>_l,isig1
           isig = -1            ! fermion helicity
           do mu = 0,5
              jqq(mu,isig,1,is1,l) = jh1(mu,isig) + jh2(mu,isig)
c        		  = (<2|gam.mu|q,1>+<2,q|gam.mu|1>)_l,isig1
           enddo
         
c	gluon from lower line:
           call curr6(-1,psi(1,-1,4),p(0,4),braketg(1,-1,3,l),pq(0,3),jh1)	
c                                           =   <4|gam.mu|q,3>_l,isig3
           call curr6(-1,braketg(1,-1,4,l),pq(0,4),psi(1,-1,3),p(0,3),jh2)	
c                                           =   <4,q|gam.mu|3>_l,isig3
           isig = -1
           do mu = 0,5
              jqq(mu,isig,2,is3,l) = jh1(mu,isig) + jh2(mu,isig)
c       		  = (<4|gam.mu|q,3>+<4,q|gam.mu|3>)_l,isig3
           enddo
        enddo  ! l loop


c -----------------------------------------------------------------------
c -----------------------------------------------------------------------
c
c contract with vvtoww tensors to get Vertex-Vertex scattering diagrams
c
        do l = 1,2              ! gluon polarization (kart. basis)
      
           if (jlog2) then
              j =  2            ! g from upper line
            
              mww = contract_Tjj(wppwpp(0,0,j),jqq(0,-1,1,is1,l),
     -             jqq(0,-1,2,is3,0))
              mat(j,l,1) = mww*clr(3,3,-1)**2   

           endif
  
           if (jlog3) then
              j =  3            ! g from lower line
              
              mww = contract_Tjj(wppwpp(0,0,j),jqq(0,-1,1,is1,0),
     -             jqq(0,-1,2,is3,l))
              mat(j,l,1) = mww*clr(3,3,-1)**2
              
           endif
        enddo                   ! l-loop
     
c ----------------------------------------------------------------    
c
c check result against madgraph output:

        if (vvdebug) then
           do lh = -1,1,2	! gluon pol. in hel. basis
              
	! need to distingiush between j=2,3 for gsign = -1
              m1hj(lh,2) = 0    ! g to upper line (j=2)
              m1hj(lh,3) = 0    ! g to lower line (j=3)	     
              temp = m1hj(lh,2)
              
c compute reference result from individual madgraph diagrams tamwpr(i,...):
 
              do i = 261,320
                 m = 0
                 do n = 1,30
                    if (i.eq.list1(n)) then ! g from upper line
                       temp = m1hj(lh,2) + tamwpr(i,-1,-1,lh)
                       m = m+1
                    endif    
                 enddo          !n
                 m1hj(lh,2) = temp	
                 if(m.eq.0) then ! g from lower line
                    m1hj(lh,3) = m1hj(lh,3) + tamwpr(i,-1,-1,lh) 
                 endif
              enddo             !i
              
              do i = 351,386
                 m = 0
                 do n = 1,18
                    if (i.eq.list2(n)) then
                       temp = m1hj(lh,2) + tamwpr(i,-1,-1,lh)
                       m = m+1
                    endif    
                 enddo          !n
                 m1hj(lh,2) = temp	
                 if(m.eq.0) then
                    m1hj(lh,3) = m1hj(lh,3) + tamwpr(i,-1,-1,lh) 
                 endif
              enddo             !i	  
              do i = 417,452
                 m = 0
                 do n = 1,18
                    if (i.eq.list3(n)) then
                       temp = m1hj(lh,2) + tamwpr(i,-1,-1,lh)
                       m = m+1
                    endif    
                 enddo          !n
                 m1hj(lh,2) = temp	
                 if(m.eq.0) then
                    m1hj(lh,3) = m1hj(lh,3) + tamwpr(i,-1,-1,lh) 
                 endif
              enddo             !i	  
              
              if(.not.jlog2) m1hj(lh,2) = 0d0
              if(.not.jlog3) m1hj(lh,3) = 0d0
              
           enddo                ! end of lh loop 
           
! compute m1 in kartesian basis:
           do j = 2,3
              m1kj(1,j) = (m1hj(-1,j)-m1hj(1,j))/sqrt(2d0)
              m1kj(2,j) = -im*(m1hj(-1,j)+m1hj(1,j))/sqrt(2d0)
           enddo    	  
		    
		    
c print results from checks:

           print*,'signs:',sign	
           print*,'gign:',gsign	    
           print*,' vv results  '
           
           do l = 1,2
              print*,'ratio for l=',l
              print*,
     -             (mat(2,l,1))/(m1kj(l,2)),
     -             (mat(3,l,1))/(m1kj(l,3))
c    	    print*,'qq results for l=',l
c	    print*,
c     -	    	(mat(2,l,1)),(mat(3,l,1))
c    	    print*,'mg results for l=',l
c	    print*,
c     -	    	(m1kj(l,2)),(m1kj(l,3))
           enddo                !l	  
           print*	
            
        endif                   ! end of debugging block vdebug
      
c -----------------------------------------------------------------------
c -----------------------------------------------------------------------
c
c  prepare box diagrams: attach Ws to external spinors
c 
c      isig = -1   : lefthanded spinors only coupling to W's
c
c  Notation for virtual 2-component spinors and momenta
c
c  W->ev   attached to quark number i: psiwe(*,i) with momentum fqe(mu,i)
c  W->muvm attached to quark number i: psiwu(*,i) with momentum fqu(mu,i)
c  
c  the fermion current corresponding to a quark line with the real emitted W->ev
c  attached next to quark number i is stored in jwe(mu,isig,is,i,l). Similarly
c  jwu(mu,isig,i,l) is the corresponding current for a W->muvm attached next to 
c  quark i; for l = 0 jwe/jwu doesn't contain a gluon; for l = 1,2 jwe/jwu 
c  includes sum of possible gluon couplings on quark leg i
c 
        l = 0                   ! no gluon 
 
 	do i = 1,3,2            ! fermion ID
           
           if (i.eq.1) then
              is = is1
           else
              is = is3   
           endif   
           
	! only fermions with isig = -1 couple to W (-> no loop for isig)
c$$$           write(*,*)'i, p(0,i) =', i, p(0,i)
c$$$           write(*,*)'i, qe,we =', i, qe,we
c$$$           write(*,*)'i, psi(1,-1,i) =', i, psi(1,-1,i)
c$$$           write(*,*)'i, psiwe(1,i) =', i, psiwe(1,i)
	
           call ket2c(psi(1,-1,i),.true.,p(0,i),-1,qe,we,
     1          psiwe(1,i),fqe(0,i))
           call ket2c(psi(1,-1,i),.true.,p(0,i),-1,qu,wu,
     1          psiwu(1,i),fqu(0,i))
           call bra2c(psi(1,-1,i+1),.true.,p(0,i+1),-1,qe,we,
     1          psiwe(1,i+1),fqe(0,i+1))
           call bra2c(psi(1,-1,i+1),.true.,p(0,i+1),-1,qu,wu,
     1          psiwu(1,i+1),fqu(0,i+1))

c$$$           write(*,*)'i, fqe(0,i) =', i, fqe(0,i)
c$$$           write(*,*)'i, psi(1,-1,i) =', i, psi(1,-1,i)
c$$$           write(*,*)'i, psiwe(1,i) =', i, psiwe(1,i)
           
	  
           call curr6(-1,psi(1,-1,i+1),p(0,i+1),
     1          psiwe(1,i),fqe(0,i), jwe(0,-1,is,i,l)   )
           call curr6(-1,psiwe(1,i+1),fqe(0,i+1),
     1          psi(1,-1,i),p(0,i), jwe(0,-1,is,i+1,l) )
           call curr6(-1,psi(1,-1,i+1),p(0,i+1),
     1          psiwu(1,i),fqu(0,i), jwu(0,-1,is,i,l)   )
           call curr6(-1,psiwu(1,i+1),fqu(0,i+1),
     1          psi(1,-1,i),p(0,i), jwu(0,-1,is,i+1,l) )

c$$$           write(*,*)'i, psi(1,-1,i+1) =', i, psi(1,-1,i+1)
c$$$           write(*,*)'i, p(0,i+1) =', i, p(0,i+1)
c$$$           write(*,*)'i, jwe(0,-1,is,i,l) =', i, jwe(0,-1,is,i,l)
      
      
 	enddo                   ! i loop
*        stop
 
c -----------------------------------------------------------------------
c
c keep structure of LO code, but replace jwe/jwu with 
c jwe/jwu(mu,isig,is,i,l=0:2) to take gluon radiation into account


c attach W->ev to f1 or f2 / f3 or f4:

 	do i = 1,3,2            ! fermion ID (isigi=-1)  
           
           if (i.eq.1) then
              is = is1
           else
              is = is3   
           endif   
           
 
 ! gluon radiation from fermion i / i+1
           do l = 1,2
              call ket2c(braketg(1,-1,i,l),.false.,pq(0,i),
     $             -1,qe,we,braketgwe(1,-1,i,l),pgwe(0,i))
              call bra2c(braketg(1,-1,i+1,l),.false.,pq(0,i+1),
     $             -1,qe,we,braketgwe(1,-1,i+1,l),pgwe(0,i+1))
              
              call ket2r(psiwe(1,i),.false.,fqe(0,i),-1,
     $             q,eps(0,l),braketweg(1,-1,i,l),pewg(0,i))      
              call bra2r(psiwe(1,i+1),.false.,fqe(0,i+1),-1,
     $             q,eps(0,l),braketweg(1,-1,i+1,l),pewg(0,i+1))      
              
              
              call curr6(-1,psi(1,-1,i+1),p(0,i+1),
     $             braketgwe(1,-1,i,l),pgwe(0,i),jwgi) 		
              call curr6(-1,braketgwe(1,-1,i+1,l),pgwe(0,i+1),
     $             psi(1,-1,i),p(0,i),jwgii)		       
              
              call curr6(-1,psi(1,-1,i+1),p(0,i+1),
     $             braketweg(1,-1,i,l),pewg(0,i),jgwi)
              call curr6(-1,braketweg(1,-1,i+1,l),pewg(0,i+1),
     $             psi(1,-1,i),p(0,i),jgwii)
 
  ! gluon radiation from fermion i+1 / i
              call curr6(-1,braketg(1,-1,i+1,l),pq(0,i+1),
     $             psiwe(1,i),fqe(0,i),jwg0i)
              call curr6(-1,psiwe(1,i+1),fqe(0,i+1),
     $             braketg(1,-1,i,l),pq(0,i),jwg0ii)
              

              do mu = 0,5
                 
                 jwe(mu,-1,is,i,l) = jwgi(mu,-1)+
     $                jgwi(mu,-1)+jwg0i(mu,-1) ! W+ & g emission from i line 
                 
                 jwe(mu,-1,is,i+1,l) = jwgii(mu,-1)+
     $                jgwii(mu,-1)+jwg0ii(mu,-1) ! W+ & g emission from i+1 line 
                 
              enddo             !mu       
              
           enddo                ! l = 1,2
c     
        enddo                   ! i loop
c
c ----------------------------------
c
c attach W->muvm to f1 or f2 / f3 or f4:

 	do i = 1,3,2            ! fermion ID (isigi=-1)  
		
           if (i.eq.1) then
              is = is1
           else
              is = is3   
           endif   
           
           
 ! gluon radiation from fermion i / i+1
           do l = 1,2
              call ket2c(braketg(1,-1,i,l),.false.,pq(0,i),
     $             -1,qu,wu,braketgwu(1,-1,i,l),pgwu(0,i))
              call bra2c(braketg(1,-1,i+1,l),.false.,pq(0,i+1),
     $             -1,qu,wu,braketgwu(1,-1,i+1,l),pgwu(0,i+1))
              
              call ket2r(psiwu(1,i),.false.,fqu(0,i),-1,
     $             q,eps(0,l),braketwug(1,-1,i,l),puwg(0,i))      
              call bra2r(psiwu(1,i+1),.false.,fqu(0,i+1),-1,
     $             q,eps(0,l),braketwug(1,-1,i+1,l),puwg(0,i+1))      
              
              call curr6(-1,psi(1,-1,i+1),p(0,i+1),
     $             braketgwu(1,-1,i,l),pgwu(0,i),jwgi) 		
              call curr6(-1,braketgwu(1,-1,i+1,l),pgwu(0,i+1),
     $             psi(1,-1,i),p(0,i),jwgii)		       
              
              call curr6(-1,psi(1,-1,i+1),p(0,i+1),
     $             braketwug(1,-1,i,l),puwg(0,i),jgwi)
              call curr6(-1,braketwug(1,-1,i+1,l),puwg(0,i+1),
     $             psi(1,-1,i),p(0,i),jgwii)
              
  ! gluon radiation from fermion i+1 / i
              call curr6(-1,braketg(1,-1,i+1,l),pq(0,i+1),
     $             psiwu(1,i),fqu(0,i),jwg0i)
              call curr6(-1,psiwu(1,i+1),fqu(0,i+1),
     $             braketg(1,-1,i,l),pq(0,i),jwg0ii)
              
              
              do mu = 0,5
                 
                 jwu(mu,-1,is,i,l) = jwgi(mu,-1)+jgwi(mu,-1)+
     $                jwg0i(mu,-1) ! W- & g emission from i line 
                 
                 jwu(mu,-1,is,i+1,l) = jwgii(mu,-1)+jgwii(mu,-1)+
     $                jwg0ii(mu,-1) ! W- & g emission from i+1 line 
                 
              enddo             !mu       
              
           enddo                ! l = 1,2
c     
        enddo                   ! i loop
c
c
c -----------------------------------------------------------------------
c -----------------------------------------------------------------------
c
c  next do the box-box graphs with one W emitted from the upper and the 
c  other from the lower line. These are only possible for lefthanded quarks
c  on both lines, i.e. isig1 = -1 and isig3 = -1
c
c  for the q^mu*q^nu/M_V^2 terms in the gauge boson propagators we need zeu/zue:
c
        l = 0
        
        jj21u(l) = dotcc(wu,jqq(0,-1,1,is1,l))
        jj43e(l) = dotcc(we,jqq(0,-1,2,is3,l))
        
        jj21e(l) = dotcc(we,jqq(0,-1,1,is1,l))
        jj43u(l) = dotcc(wu,jqq(0,-1,2,is3,l))
        
        do l = 1,2
           
           jj21u(l) = dotcc(wu,jqq(0,-1,1,is1,l))
           jj43e(l) = dotcc(we,jqq(0,-1,2,is3,l))
           
           jj21e(l) = dotcc(we,jqq(0,-1,1,is1,l))
           jj43u(l) = dotcc(wu,jqq(0,-1,2,is3,l))
           
           
           zeu(2,l) = jj21e(l)*jj43u(0)*zm2i(2)
           zue(2,l) = jj21u(l)*jj43e(0)*zm2i(2)
           
           zeu(3,l) = jj21e(0)*jj43u(l)*zm2i(2)
           zue(3,l) = jj21u(0)*jj43e(l)*zm2i(2)
           
        enddo                   !l
c
      
        do l = 1,2
c We from upper line, Wu from lower line:
           do ide = 1,2         ! fermion id for external We coupling
              do idu = 3,4      ! fermion id for external Wu coupling
                 j = 2
                 mb(ide,idu,j,l) = -dotcc(jwe(0,-1,is1,ide,l),
     -                jwu(0,-1,is3,idu,0))
                 j = 3
                 mb(ide,idu,j,l) = -dotcc(jwe(0,-1,is1,ide,0),
     -                jwu(0,-1,is3,idu,l))

c$$$             write(*,*)'ide,idu =', ide,idu
c$$$             write(*,*)'jwe(0,-1,is1,ide,l) =', jwe(0,-1,is1,ide,l)
c$$$             write(*,*)'jwu(0,-1,is3,idu,0) =', jwu(0,-1,is3,idu,0)
c$$$             write(*,*)'jwe(0,-1,is1,ide,0) =', jwe(0,-1,is1,ide,0)
c$$$             write(*,*)'jwu(0,-1,is3,idu,l) =', jwu(0,-1,is3,idu,l)

              enddo             !idu
           enddo                !ide
           
c Wu from upper line, We from lower line:
           do ide = 3,4         ! fermion id for external We coupling
              do idu = 1,2      ! fermion id for external Wu coupling
                 j = 2
                 mb(ide,idu,j,l) = -dotcc(jwu(0,-1,is1,idu,l),
     -                jwe(0,-1,is3,ide,0))
                 j = 3
                 mb(ide,idu,j,l) = -dotcc(jwu(0,-1,is1,idu,0),
     -                jwe(0,-1,is3,ide,l))

              enddo             !idu
           enddo                !ide
           
        enddo                   !l
      
      
        do j = jmin,jmax
           do l = 1,2
              maeu(j,l) = 0d0
              mzeu(j,l) = 0d0
              maue(j,l) = 0d0
              mzue(j,l) = 0d0
              
              do ide = 1,2      ! fermion id for external We coupling
                 do idu = 3,4   ! fermion id for external Wu coupling
                    
                    maeu(j,l)= maeu(j,l) + mb(ide,idu,j,l)*
     &                   clr(ifl(3-ide),1,-1)*clr(ifl(7-idu),1,-1)
                    
                    mzeu(j,l) = mzeu(j,l) + ( mb(ide,idu,j,l) -
     &                   zeu(j,l)*(-2*ide+3)*(-2*idu+7) )*
     &                   clr(ifl(3-ide),2,-1)*clr(ifl(7-idu),2,-1)
                    
                 enddo          !idu
              enddo             !ide
              
              do idu = 1,2 
                 do ide = 3,4 
                    
                    maue(j,l) = maue(j,l) + mb(ide,idu,j,l)*
     &                   clr(ifl(7-ide),1,-1)*clr(ifl(3-idu),1,-1)
                    mzue(j,l) = mzue(j,l) + ( mb(ide,idu,j,l) -
     &                   zue(j,l)*(-2*idu+3)*(-2*ide+7) )*
     &                   clr(ifl(7-ide),2,-1)*clr(ifl(3-idu),2,-1)
c$$$                    write(*,*)'j,l,maue =', j,l,maue(j,l)
c$$$                    write(*,*)'j,l,mzue =', j,l,mzue(j,l)
                    
                 enddo          !ide
              enddo             !idu
              
      
              mat(j,l,6) = -(maeu(j,l)*prop_ee(1,j)+mzeu(j,l)*prop_ee(2,j)+
     $             maue(j,l)*prop_uu(1,j)+mzue(j,l)*prop_uu(2,j))

*              write(*,*)'j,l,mat =', j,l,mat(j,l,6)
      
           enddo                !l      
        enddo                   !j
*        stop
      
c ----------------------------------------------------------------    
c
c debugging:
c
c check result as obtained above against madgraph output 
c		(sum up relevant mg amplitudes in the following):

        if(bbdebug) then
c		(sum up relevant mg amplitudes in the following):
	
           do lh = -1,1,2	! gluon pol. in hel. basis
 	    
	! need to distingiush between j=2,3
		
c compute reference result from individual madgraph diagrams tamwpr(i,...):
	
              do j = 2,3
                 m1haj(lh,j) = 0d0 ! graphs with A exchange
                 m1hzj(lh,j) = 0d0 !		 Z exchange
              enddo  
              
		
              do m = 1,24       ! g from upper line
                 i = listab2(m)
                 m1haj(lh,2) = m1haj(lh,2) + tamwpr(i,-1,-1,lh) 
              enddo
              
              do m = 1,24       ! g from lower line
                 i = listab3(m)
                 m1haj(lh,3) = m1haj(lh,3) + tamwpr(i,-1,-1,lh) 
              enddo
              
              do m = 1,24       ! g from upper line
                 i = listzb2(m)
                 m1hzj(lh,2) = m1hzj(lh,2) + tamwpr(i,-1,-1,lh) 
              enddo
              
              do m = 1,24       ! g from lower line
                 i = listzb3(m)
                 m1hzj(lh,3) = m1hzj(lh,3) + tamwpr(i,-1,-1,lh) 
              enddo
              
	    ! for incoming g: no g emission from other line
              if(.not.jlog2) m1haj(lh,2) = 0d0
              if(.not.jlog3) m1haj(lh,3) = 0d0
              if(.not.jlog2) m1hzj(lh,2) = 0d0
              if(.not.jlog3) m1hzj(lh,3) = 0d0
              
           enddo                ! end of lh loop 	  
           
	  		    
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
           
           
c     ---------------------------------		 

           do l = 1,2 
              
              print*,'box-box contributions:'
c     print*,'signs:',sign	
c     print*,'gign:',gsign	    
              print*,' l = ',l     
              
              do j = jmin,jmax	    
                 
                 print*,' j = ',j	  
                 print*,'mg:',m1kaj(l,j)+m1kzj(l,j)
                 print*,'mat = ',mat(j,l,6)
                 print*,'ratio = ',(m1kaj(l,j)+m1kzj(l,j))/mat(j,l,6)
                 print*
                 
              enddo             !j
              
           enddo                !l
           
          
        endif                   ! debugging  


c --------------------------------------------------------------------------
c --------------------------------------------------------------------------
c
c  next do the vertex-box graphs 
c 	with one boson (We/Wu) emitted from a quark line and
c	one lepton-neutrino pair from the t-channel:
c     	
c
c lk=1 and kk=2 is for external boson from upper line
c lk=2 and kk=1 is for external boson from lower line

      ! l = 0:
          
       lk = 1
       kk = 3-lk
       is = is3
        ! VW-> e ve:	
        call contract_T1j(awen(0,0,lk,1+lk),jqq(0,-1,kk,is,0),
     $       				    epsawe(0,lk,1+lk,0)) !for AW
        call contract_T1j(zwen(0,0,lk,1+lk),jqq(0,-1,kk,is,0),
     $       				    epszwe(0,lk,1+lk,0)) !for ZW
        ! VW-> mu vm:
        call contract_T1j(awmu(0,0,lk,1+lk),jqq(0,-1,kk,is,0),
     $       				    epsawu(0,lk,1+lk,0)) !for AW
        call contract_T1j(zwmu(0,0,lk,1+lk),jqq(0,-1,kk,is,0),
     $       				    epszwu(0,lk,1+lk,0)) !for ZW
     
        qepszwe(lk,1+lk,0) = -dotrc(quu(0,1+lk),epszwe(0,lk,1+lk,0))*zm2i(2)
        qepszwu(lk,1+lk,0) = -dotrc(qee(0,1+lk),epszwu(0,lk,1+lk,0))*zm2i(2)
     
        do l = 1,2
           call contract_T1j(awen(0,0,lk,4-lk),jqq(0,-1,kk,is,l),
     $  				     epsawe(0,lk,4-lk,l)) 
           call contract_T1j(zwen(0,0,lk,4-lk),jqq(0,-1,kk,is,l),
     $  				     epszwe(0,lk,4-lk,l)) 
           call contract_T1j(awmu(0,0,lk,4-lk),jqq(0,-1,kk,is,l),
     $  				     epsawu(0,lk,4-lk,l)) 
           call contract_T1j(zwmu(0,0,lk,4-lk),jqq(0,-1,kk,is,l),
     $  				     epszwu(0,lk,4-lk,l)) 
         	 
	   qepszwe(lk,4-lk,l) = -dotrc(quu(0,4-lk),epszwe(0,lk,4-lk,l))*zm2i(2)
	   qepszwu(lk,4-lk,l) = -dotrc(qee(0,4-lk),epszwu(0,lk,4-lk,l))*zm2i(2)
	   
		
	enddo !l

        	          	   
	lk = 2
        kk = 3-lk
        is = is1
	 ! VW-> e ve: 	
	 call contract_T1j(waen(0,0,lk,1+lk),jqq(0,-1,kk,is,0),
     $       				    epswae(0,lk,1+lk,0)) !for WA
	 call contract_T1j(wzen(0,0,lk,1+lk),jqq(0,-1,kk,is,0),
     $        				    epswze(0,lk,1+lk,0)) !for WZ
	 ! VW-> mu vm:
	 call contract_T1j(wamu(0,0,lk,1+lk),jqq(0,-1,kk,is,0),
     $       				    epswau(0,lk,1+lk,0)) !for WA
	 call contract_T1j(wzmu(0,0,lk,1+lk),jqq(0,-1,kk,is,0),
     $       				    epswzu(0,lk,1+lk,0)) !for WZ
	
         qepszwe(lk,1+lk,0) = -dotrc(qee(0,1+lk),epswze(0,lk,1+lk,0))*zm2i(2)
         qepszwu(lk,1+lk,0) = -dotrc(quu(0,1+lk),epswzu(0,lk,1+lk,0))*zm2i(2)
    
        do l = 1,2
           call contract_T1j(waen(0,0,lk,4-lk),jqq(0,-1,kk,is,l),
     $  		        	    epswae(0,lk,4-lk,l)) 
           call contract_T1j(wzen(0,0,lk,4-lk),jqq(0,-1,kk,is,l),
     $  		        	    epswze(0,lk,4-lk,l)) 
           call contract_T1j(wamu(0,0,lk,4-lk),jqq(0,-1,kk,is,l),
     $  		        	    epswau(0,lk,4-lk,l)) 
           call contract_T1j(wzmu(0,0,lk,4-lk),jqq(0,-1,kk,is,l),
     $  				    epswzu(0,lk,4-lk,l))
         	 
	   qepszwe(lk,4-lk,l) = -dotrc(qee(0,4-lk),epswze(0,lk,4-lk,l))*zm2i(2)
	   qepszwu(lk,4-lk,l) = -dotrc(quu(0,4-lk),epswzu(0,lk,4-lk,l))*zm2i(2)
	   
	enddo !l
	

ccccccccccccccc

	   
c gauge boson from upper quark line:	   
        lk = 1 
          
        do l = 1,2
	  
           if (jlog2) then
              j = 2
              m2we(j,l) = dotcc(jwe(0,-1,is1,1,l),epszwu(0,lk,j,0))    
              m1we(j,l) = dotcc(jwe(0,-1,is1,2,l),epszwu(0,lk,j,0))	    
              w2we(j,l) = qepszwu(lk,j,0)*jj21e(l)    
              
              m2ze(j,l) = (m2we(j,l)-w2we(j,l))*clr(ifl(1),3,-1)*
     -             clr(ifl(2),2,-1)	    
              m2ae(j,l) = dotcc(jwe(0,-1,is1,1,l),epsawu(0,lk,j,0))*
     &             clr(ifl(1),3,-1)*clr(ifl(2),1,-1)	  
              m2e(j,l)  = m2ze(j,l)+m2ae(j,l)
	      
              m1ze(j,l) = (m1we(j,l)+w2we(j,l))*clr(ifl(2),3,-1)*
     -             clr(ifl(1),2,-1)	  
              m1ae(j,l) = dotcc(jwe(0,-1,is1,2,l),epsawu(0,lk,j,0))*
     &             clr(ifl(2),3,-1)*clr(ifl(1),1,-1)	  
              m1e(j,l) = m1ze(j,l)+m1ae(j,l)
              
              
              m2wu(j,l) = dotcc(jwu(0,-1,is1,1,l),epszwe(0,lk,j,0))    
              m1wu(j,l) = dotcc(jwu(0,-1,is1,2,l),epszwe(0,lk,j,0))	    
              w2wu(j,l) = qepszwe(lk,j,0)*jj21u(l)    
              
              m2zu(j,l) = (m2wu(j,l)-w2wu(j,l))*clr(ifl(1),3,-1)*
     -             clr(ifl(2),2,-1)	    
              m2au(j,l) = dotcc(jwu(0,-1,is1,1,l),epsawe(0,lk,j,0))*
     &             clr(ifl(1),3,-1)*clr(ifl(2),1,-1)	  
              m2u(j,l)  = m2zu(j,l)+m2au(j,l)
	      
              m1zu(j,l) = (m1wu(j,l)+w2wu(j,l))*clr(ifl(2),3,-1)*
     -             clr(ifl(1),2,-1)	  
              m1au(j,l) = dotcc(jwu(0,-1,is1,2,l),epsawe(0,lk,j,0))*
     &             clr(ifl(2),3,-1)*clr(ifl(1),1,-1)	  
              m1u(j,l) = m1zu(j,l)+m1au(j,l)


              mat(j,l,2) =  -(m1e(j,l)+m1u(j,l)+m2e(j,l)+m2u(j,l))*clr(3,3,-1)
*              write(*,*)'j,l, mat(j,l,2) =', j,l, mat(j,l,2)
           endif                !jmin=2
	  
           if (jlog3) then
              j = 3
              m2we(j,l) = dotcc(jwe(0,-1,is1,1,0),epszwu(0,lk,j,l))    
              m1we(j,l) = dotcc(jwe(0,-1,is1,2,0),epszwu(0,lk,j,l))	    
              w2we(j,l) = qepszwu(lk,j,l)*jj21e(0)    
              
              m2ze(j,l) = (m2we(j,l)-w2we(j,l))*clr(ifl(1),3,-1)*
     -             clr(ifl(2),2,-1)	    
              m2ae(j,l) = dotcc(jwe(0,-1,is1,1,0),epsawu(0,lk,j,l))*
     &             clr(ifl(1),3,-1)*clr(ifl(2),1,-1)	  
              m2e(j,l)  = m2ze(j,l)+m2ae(j,l)
	      
              m1ze(j,l) = (m1we(j,l)+w2we(j,l))*clr(ifl(2),3,-1)*
     -             clr(ifl(1),2,-1)	  
              m1ae(j,l) = dotcc(jwe(0,-1,is1,2,0),epsawu(0,lk,j,l))*
     &             clr(ifl(2),3,-1)*clr(ifl(1),1,-1)	  
              m1e(j,l) = m1ze(j,l)+m1ae(j,l)
              
	   
              m2wu(j,l) = dotcc(jwu(0,-1,is1,1,0),epszwe(0,lk,j,l))    
              m1wu(j,l) = dotcc(jwu(0,-1,is1,2,0),epszwe(0,lk,j,l))	    
              w2wu(j,l) = qepszwe(lk,j,l)*jj21u(0)    
              
              m2zu(j,l) = (m2wu(j,l)-w2wu(j,l))*clr(ifl(1),3,-1)*
     -             clr(ifl(2),2,-1)	    
              m2au(j,l) = dotcc(jwu(0,-1,is1,1,0),epsawe(0,lk,j,l))*
     &             clr(ifl(1),3,-1)*clr(ifl(2),1,-1)	  
              m2u(j,l)  = m2zu(j,l)+m2au(j,l)
	      
              m1zu(j,l) = (m1wu(j,l)+w2wu(j,l))*clr(ifl(2),3,-1)*
     -             clr(ifl(1),2,-1)	  
              m1au(j,l) = dotcc(jwu(0,-1,is1,2,0),epsawe(0,lk,j,l))*
     &             clr(ifl(2),3,-1)*clr(ifl(1),1,-1)	  
              m1u(j,l) = m1zu(j,l)+m1au(j,l)
	  
              
              mat(j,l,2) =  -(m1e(j,l)+m1u(j,l)+m2e(j,l)+m2u(j,l))*
     -             clr(3,3,-1)
*              write(*,*)'j,l, mat(j,l,2) =', j,l, mat(j,l,2)
           endif                !jmax=3
           
        enddo                   !l 
 	   	 	 
ccccc
c	   
c gauge boson from lower quark line:	   
        lk = 2 
	   	   
        do l = 1,2
           
           if (jlog2) then
              j = 2
              m4we(j,l) = dotcc(jwe(0,-1,is3,3,0),epswzu(0,lk,j,l))    
              m3we(j,l) = dotcc(jwe(0,-1,is3,4,0),epswzu(0,lk,j,l))	    
              w4we(j,l) = -qepszwu(lk,j,l)*jj43e(0)    
              
              m4ze(j,l) = (m4we(j,l)-w4we(j,l))*clr(ifl(3),3,-1)*
     -             clr(ifl(4),2,-1)	    
              m4ae(j,l) = dotcc(jwe(0,-1,is3,3,0),epswau(0,lk,j,l))*
     &             clr(ifl(3),3,-1)*clr(ifl(4),1,-1)	  
              m4e(j,l)  = m4ze(j,l)+m4ae(j,l)
	      
              m3ze(j,l) = (m3we(j,l)+w4we(j,l))*clr(ifl(4),3,-1)*
     -             clr(ifl(3),2,-1)	  
              m3ae(j,l) = dotcc(jwe(0,-1,is3,4,0),epswau(0,lk,j,l))*
     &             clr(ifl(4),3,-1)*clr(ifl(3),1,-1)	  
              m3e(j,l) = m3ze(j,l)+m3ae(j,l)
              
              
              m4wu(j,l) = dotcc(jwu(0,-1,is3,3,0),epswze(0,lk,j,l))    
              m3wu(j,l) = dotcc(jwu(0,-1,is3,4,0),epswze(0,lk,j,l))	    
              w4wu(j,l) = -qepszwe(lk,j,l)*jj43u(0)    
              
              m4zu(j,l) = (m4wu(j,l)-w4wu(j,l))*clr(ifl(3),3,-1)*
     -             clr(ifl(4),2,-1)	    
              m4au(j,l) = dotcc(jwu(0,-1,is3,3,0),epswae(0,lk,j,l))*
     &             clr(ifl(3),3,-1)*clr(ifl(4),1,-1)	  
              m4u(j,l)  = m4zu(j,l)+m4au(j,l)
	      
              m3zu(j,l) = (m3wu(j,l)+w4wu(j,l))*clr(ifl(4),3,-1)*
     -             clr(ifl(3),2,-1)	  
              m3au(j,l) = dotcc(jwu(0,-1,is3,4,0),epswae(0,lk,j,l))*
     &             clr(ifl(4),3,-1)*clr(ifl(3),1,-1)	  
              m3u(j,l) = m3zu(j,l)+m3au(j,l)
              

              mat(j,l,3) =  -(m3e(j,l)+m3u(j,l)+m4e(j,l)+m4u(j,l))*
     -             clr(3,3,-1)
*              write(*,*)'j,l,mat(j,l,3) =', j,l,mat(j,l,3)
           endif                !jmin
 	  
           if (jlog3) then
              j = 3
              m4we(j,l) = dotcc(jwe(0,-1,is3,3,l),epswzu(0,lk,j,0))    
              m3we(j,l) = dotcc(jwe(0,-1,is3,4,l),epswzu(0,lk,j,0))	    
              w4we(j,l) = -qepszwu(lk,j,0)*jj43e(l)    
              
              m4ze(j,l) = (m4we(j,l)-w4we(j,l))*clr(ifl(3),3,-1)*
     -             clr(ifl(4),2,-1)	    
              m4ae(j,l) = dotcc(jwe(0,-1,is3,3,l),epswau(0,lk,j,0))*
     &             clr(ifl(3),3,-1)*clr(ifl(4),1,-1)	  
              m4e(j,l)  = m4ze(j,l)+m4ae(j,l)
	      
              m3ze(j,l) = (m3we(j,l)+w4we(j,l))*clr(ifl(4),3,-1)*
     -             clr(ifl(3),2,-1)	  
              m3ae(j,l) = dotcc(jwe(0,-1,is3,4,l),epswau(0,lk,j,0))*
     &             clr(ifl(4),3,-1)*clr(ifl(3),1,-1)	  
              m3e(j,l) = m3ze(j,l)+m3ae(j,l)
              
              
              m4wu(j,l) = dotcc(jwu(0,-1,is3,3,l),epswze(0,lk,j,0))    
              m3wu(j,l) = dotcc(jwu(0,-1,is3,4,l),epswze(0,lk,j,0))	    
              w4wu(j,l) = -qepszwe(lk,j,0)*jj43u(l)    
              
              m4zu(j,l) = (m4wu(j,l)-w4wu(j,l))*clr(ifl(3),3,-1)*
     -             clr(ifl(4),2,-1)	    
              m4au(j,l) = dotcc(jwu(0,-1,is3,3,l),epswae(0,lk,j,0))*
     &             clr(ifl(3),3,-1)*clr(ifl(4),1,-1)	  
              m4u(j,l)  = m4zu(j,l)+m4au(j,l)
	      
              m3zu(j,l) = (m3wu(j,l)+w4wu(j,l))*clr(ifl(4),3,-1)*
     -             clr(ifl(3),2,-1)	  
              m3au(j,l) = dotcc(jwu(0,-1,is3,4,l),epswae(0,lk,j,0))*
     &             clr(ifl(4),3,-1)*clr(ifl(3),1,-1)	  
              m3u(j,l) = m3zu(j,l)+m3au(j,l)	   
              
              mat(j,l,3) =  -(m3e(j,l)+m3u(j,l)+m4e(j,l)+m4u(j,l))*
     -             clr(3,3,-1)
*              write(*,*)'j,l,mat(j,l,3) =', j,l,mat(j,l,3)
           endif                !jmax
           
        enddo                   !l 
 	   

c ----------------------------------------------------------------    
c
c check result against madgraph output 
c		(sum up relevant mg amplitudes in the following):

        
        if (vbdebug) then

	
	  do lh = -1,1,2	! gluon pol. in hel. basis
 	    
	! need to distingiush between j=2,3
		
c compute reference result from individual madgraph diagrams tampr(i,...):
	
	    do j = 2,3
	    	m1haj(lh,j) = 0d0  ! graphs with A exchange
	    	m1hzj(lh,j) = 0d0  !		 Z exchange
           enddo  
	   
		
	  do m = 1,40
	      i = listba2(m)
	      m1haj(lh,2) = m1haj(lh,2) + tamwpr(i,-1,-1,lh) 
	  enddo
	  	     
	 do m = 1,40
	     i = listba3(m)
	     m1haj(lh,3) = m1haj(lh,3) + tamwpr(i,-1,-1,lh) 
	 enddo
	  
	   do m = 1,60
	       i = listbz2(m)
	       m1hzj(lh,2) = m1hzj(lh,2) + tamwpr(i,-1,-1,lh) 
	   enddo
c       
	   do m = 1,60
	       i = listbz3(m)
	       m1hzj(lh,3) = m1hzj(lh,3) + tamwpr(i,-1,-1,lh) 
	   enddo
     
 	! for incoming g: no g emission from other line
       	     if(.not.jlog2) m1haj(lh,2) = 0d0
       	     if(.not.jlog3) m1haj(lh,3) = 0d0
       	     if(.not.jlog2) m1hzj(lh,2) = 0d0
       	     if(.not.jlog3) m1hzj(lh,3) = 0d0
	  
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


	print*,'compare vertex-box graphs'
	do j = jmin,jmax
	   do l = 1,2	
	      print*,'j,l=',j,l
       	      print*,'mg  = ',m1kaj(l,j)+m1kzj(l,j)
       	      print*,'mat = ',mat(j,l,2)+mat(j,l,3)
       	      print*,'ratio = ',(m1kaj(l,j)+m1kzj(l,j))/
     % 	        		(mat(j,l,2)+mat(j,l,3))
     	      print*			      
	   enddo !l
	enddo !j

          
      endif ! debugging
c	
c --------------------------------------------------------------------
c --------------------------------------------------------------------
c
c compute pentagon contributions, i.e. three bosons emitted 
c from the same quark line
c
      do l = 1,2
      
         if(jlog2) then
c pentagon on upper line:
c  eps1=W,eps2=j_lower,eps3=W'
            j = 2
            m5ue(j,l) = s1c(braketgwu(1,-1,2,l),jqq(0,-1,2,is3,0),
     -           .true.,-1,psiwe(1,1))+
     $           s1c(braketwug(1,-1,2,l),jqq(0,-1,2,is3,0),.true.,-1,
     -           psiwe(1,1))+
     $           s1c(psiwu(1,2),jqq(0,-1,2,is3,0),.true.,-1,
     -           braketgwe(1,-1,1,l))+
     $           s1c(psiwu(1,2),jqq(0,-1,2,is3,0),.true.,-1,
     -           braketweg(1,-1,1,l))
            
            m5eu(j,l) = s1c(braketgwe(1,-1,2,l),jqq(0,-1,2,is3,0),
     -           .true.,-1,psiwu(1,1))+
     $           s1c(braketweg(1,-1,2,l),jqq(0,-1,2,is3,0),.true.,
     -           -1,psiwu(1,1))+
     $           s1c(psiwe(1,2),jqq(0,-1,2,is3,0),.true.,-1,
     -           braketgwu(1,-1,1,l))+
     $           s1c(psiwe(1,2),jqq(0,-1,2,is3,0),.true.,-1,
     -           braketwug(1,-1,1,l)) 

            mat(j,l,7) = clr(3,3,-1)**4*prop43(3,j)*(m5ue(j,l)+m5eu(j,l))
*            write(*,*)'j,l,mat(j,l,7) =', j,l,mat(j,l,7)

c  pentagon on lower line:      
c  eps1=W,eps2=j_upper,eps3=W'

            m5ue(j,l) = s1c(psiwu(1,4),jqq(0,-1,1,is1,l),.true.,-1,psiwe(1,3))
            m5eu(j,l) = s1c(psiwe(1,4),jqq(0,-1,1,is1,l),.true.,-1,psiwu(1,3))
            
            mat(j,l,8) = clr(3,3,-1)**4*prop21(3,j)*(m5ue(j,l)+m5eu(j,l))
*            write(*,*)'j,l,mat(j,l,8) =', j,l,mat(j,l,8)
      endif !jmin


      if(jlog3) then
c pentagon on upper line:
c  eps1=v,eps2=j_lower,eps3=v'
         j = 3
         m5ue(j,l) = s1c(psiwu(1,2),jqq(0,-1,2,is3,l),.true.,-1,psiwe(1,1))
         m5eu(j,l) = s1c(psiwe(1,2),jqq(0,-1,2,is3,l),.true.,-1,psiwu(1,1))

         mat(j,l,7) = clr(3,3,-1)**4*prop43(3,j)*(m5ue(j,l)+m5eu(j,l))
*            write(*,*)'j,l,mat(j,l,7) =', j,l,mat(j,l,7)
         
c     pentagon on lower line:
c  eps1=v,eps2=j_upper,eps3=v'
         m5ue(j,l) = s1c(braketgwu(1,-1,4,l),jqq(0,-1,1,is1,0),.true.,
     -           -1,psiwe(1,3))+
     $        s1c(braketwug(1,-1,4,l),jqq(0,-1,1,is1,0),.true.,-1,
     -           psiwe(1,3))+
     $        s1c(psiwu(1,4),jqq(0,-1,1,is1,0),.true.,-1,
     -           braketgwe(1,-1,3,l))+
     $        s1c(psiwu(1,4),jqq(0,-1,1,is1,0),.true.,-1,
     -           braketweg(1,-1,3,l))
         
         m5eu(j,l) = s1c(braketgwe(1,-1,4,l),jqq(0,-1,1,is1,0),.true.,
     -           -1,psiwu(1,3))+
     $        s1c(braketweg(1,-1,4,l),jqq(0,-1,1,is1,0),.true.,-1,
     -           psiwu(1,3))+
     $        s1c(psiwe(1,4),jqq(0,-1,1,is1,0),.true.,-1,
     -           braketgwu(1,-1,3,l))+
     $        s1c(psiwe(1,4),jqq(0,-1,1,is1,0),.true.,-1,
     -           braketwug(1,-1,3,l))

         mat(j,l,8) = clr(3,3,-1)**4*prop21(3,j)*(m5ue(j,l)+m5eu(j,l))
*            write(*,*)'j,l,mat(j,l,8) =', j,l,mat(j,l,8)
         endif                  !jmax
      
      
      enddo                     !l
c
c ----------------------------------------------------------------    
c
c check result against madgraph output 
c		(sum up relevant mg amplitudes in the following):


      if (pdebug) then
	
	  do lh = -1,1,2	! gluon pol. in hel. basis
 	    
c compute reference result from individual madgraph diagrams tampr(i,...):
	
	   do j = 2,3
	    	m1hj(lh,j) = 0d0  
           enddo  
	   
	   do m = 1,12
	      i = listp3(m)
	      m1hj(lh,3) = m1hj(lh,3) + tamwpr(i,-1,-1,lh) 
	   enddo
	  
	   do m = 1,12
	      i = listp2(m)
	      m1hj(lh,2) = m1hj(lh,2) + tamwpr(i,-1,-1,lh) 
	   enddo
 
     
 	! for incoming g: no g emission from other line
       	     if(.not.jlog2) m1hj(lh,2) = 0d0
       	     if(.not.jlog3) m1hj(lh,3) = 0d0
	  
	  enddo ! end of lh loop 	  	    
		    
    ! compute m1 in kartesian basis:
     		    
	    do l = 1,2
	      do j = 2,3
	    	m1kj(l,j) = 0d0
	      enddo
	    enddo  	    
	    
	  do j = 2,3
	    m1kj(1,j) =     (m1hj(-1,j)-m1hj(1,j))/sqrt(2d0)
	    m1kj(2,j) = -im*(m1hj(-1,j)+m1hj(1,j))/sqrt(2d0)
	  enddo    	


	print*,'compare pentagon graphs'
	do j = jmin,jmax
	   do l = 1,2	
	      print*,'j,l=',j,l
       	      print*,'mg  = ',m1kj(l,j)
       	      print*,'mat = ',mat(j,l,7)+mat(j,l,8)
       	      print*,'ratio = ',(m1kj(l,j))/
     % 	        		(mat(j,l,7)+mat(j,l,8))
     	      print*			      
	   enddo !l
	enddo !j

          
      endif ! debugging
c --------------------------------------------------------------------------
c --------------------------------------------------------------------------

c sum the graphs, square them and map them onto ucds, uuss etc.
      do j = 2,3     
         res(j) = 0
         do l = 1,2
            mm(j,l) = 0
            do i = 1,8
               mm(j,l) = mm(j,l) + mat(j,l,i)
*               write(*,*)'j,l,i, mat =', j, l, i, mat(j,l,i)
            enddo !i
            res(j) = res(j) 
     &  	      + dreal(mm(j,l))**2
     &  	      + dimag(mm(j,l))**2
         enddo !l
         res(j) = res(j)*12d0*fpials(j)   ! C_2*9 is the color factor
c$$$         write(*,*)'j, fpials(j)  =', j, fpials(j) 
c$$$         write(*,*)'j, res =', j, res(j)
      enddo !j

*      stop
		    
      if (ldebug) then
        do j = jmin,jmax
	 print*,'j=',j
         print*,'ucds(j):',ucds(j),res(j),(res(j))/ucds(j)
	enddo 
 	print*,'sum j = 2,3:'
        print*,'ucds(2+3):',ucds(1),(res(2)+res(3))/ucds(1)
      endif

      if (ldebug) then
            if (((res(2)+res(3))/ucds(1)-1).gt.1d-6) then
	    	  print*,'gsign=',gsign
                  print*,' qq,mg',':',(res(2)+res(3)),ucds(1)
                  print*,' res-ratio',':',(res(2)+res(3))/ucds(1)-1
            endif
      endif ! ldebug

      if (ldoblha) then
        if (blha_amptype(blha_curproc).eq.0) then ! tree
          if ( (blha_idsubproc(blha_cursubproc,blha_curproc)/10 .eq. 1) .or.
     &         (blha_idsubproc(blha_cursubproc,blha_curproc)/10 .eq. 2) ) then
            blha_amp(1)=blha_amp(1)+blha_CF/(4d0/3d0)*
     &       res(2)
          endif
          if ( (blha_idsubproc(blha_cursubproc,blha_curproc)/10 .eq. 1) .or.
     &         (blha_idsubproc(blha_cursubproc,blha_curproc)/10 .eq. 3) ) then
            blha_amp(1)=blha_amp(1)+blha_CF/(4d0/3d0)*
     &       res(3)
          endif

        else 
          call BLHA_amptypeerror(blha_amptype(blha_curproc),__FILE__,__LINE__)
        endif
      endif

      do j = 2,3
         ucds(j) = res(j) 
      enddo
      ucds(1) = res(3) +  res(2)

      return
      end

