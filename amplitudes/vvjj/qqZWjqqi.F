C********************  qqVqqj ***********************************

      subroutine qqZqqj(pbar, fsign, qbar, gsign, 
     1                  uucc, uuss, ddcc, ddss, udsc, ducs)
          use globalvars, only: ldoblha
      implicit none
      real*8 uucc,uuss,ddcc,ddss,
     1       udsc,ducs, pbar(0:3,6), qbar(0:4), 
     2       uucs,ddcs,udcc,udss,
     3       uusc,ddsc,ducc,duss
      real*8 uucc_c(3),uuss_c(3),ddcc_c(3),ddss_c(3),
     1       udsc_c(3),ducs_c(3),
     2       uucs_c(3),ddcs_c(3),udcc_c(3),udss_c(3),
     3       uusc_c(3),ddsc_c(3),ducc_c(3),duss_c(3)
      integer fsign(6), gsign
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

C     QQVQQj calculates of matrix elements for electroweak single weak 
c     boson production by quark quark scattering
c     i.e. weak boson fusion
c
C        q1 q3    ---->   q2 q4 V gluon ,   V ---> f5-bar f6
c
C     and crossing related processes. Pauli interference terms for
c     identical fermions are neglected. In particular, only the
c     t-channel exchange of elctroweak bosons is considered.
c     s-channel production of weak boson pairs is NOT implemented.
c
c  There are two options for selecting helicties: 1) do loop summation
c  to provide unpolarized cross sections, and 2) calculation of 1 chirality
c  combination only. See common block /chelsum/ below. uucc etc. are the 
c  1) polarization summed cross section or 2) the result for the selected
c  helicity combination
c
c	Dieter Zeppenfeld, <dieter@pheno.physics.wisc.edu>
c	Initial version:  2002 October 24
c	Last modified:    2003 January 5
c                         2010 January 13 (implemented Ajj process)
c
C     QQBQQI must be called first to initialize the desired
c     production process. The main task of QQBQQI is to set up 
c     the products of coupling constants needed in Feynman graphs 
c     a, ..., g which are stored in fcpl(4,4,8,6,7) and which are 
c     calculated from output of KOPPLN
      real*8 fcpl(4,4,8,6,7)
      integer vlmin(6,7), vlmax(6,7), vhmin(6,7), vhmax(6,7), 
     1     nprmax, sigmax, sig1(8,4), sig3(8,4), sig5(8,4)
      common /cqqbqq/ fcpl, nprmax, sigmax,
     1     vlmin, vlmax, vhmin, vhmax, sig1, sig3, sig5
      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)

c  helicity selection
      integer jsig, jsig1, jsig3, jsig5, jsigmin, jsigmax
      common /chelsum/ jsig,jsig1,jsig3,jsig5
C Select helicity combination by setting jsig externally. This program will
c then fill jsig1,jsig3,jsig5 to correspond to the chiralities of the 1-2,
c 3-4 and 5-6 fermion lines for external use, e.g. in calculating polarized
c cross sections. This is option 2 and it requires that jsig lies within
c the allowed sig range [1,sigmax] for the process. If jsig is outside this
c range the sum, NOT average!, over all polarization states will be calculated
c
c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"
      real*8 fpials(3), fpi
      PARAMETER (fpi=4d0*pi)
c   local variables
      real*8 p(0:3,6), p21(0:4), p43(0:4), p65(0:4), dummy(0:4),
     $     q(0:4), p21g(0:4), p43g(0:4), pq(0:4,4),
     $     eps(0:3,3), pk(0:4), res(6,2,-1:1)
      complex*16 psi(2,-1:1,6), j21(0:3,-1:1), j43(0:3,-1:1),
     $     j65(0:3,-1:1), prop65(4), propl(4,3,3), proph(4,3,3),
     $     mm3(6,8,3), matr(3), fac(6,8,7,3,3), mm2(6,8), s1c, matr1,
     $     mm4(6,8,3,3), matr2, matr3, matr4, matr7(2,3,3), corr56,
     $     prop21(4), prop43(4), prop21g(4), prop43g(4),
     $     braketg(2,-1:1,4,3), jh1(0:3,-1:1), jh2(0:3,-1:1),
     $     e21(0:3,-1:1,3), e43(0:3,-1:1,3), braketV(2), braketVg(2),
     $     braketgV(2), braket(2,-1:1,4)
      external s1c
      integer i, ibos, ic, j, jj, k, l, mu, vl, vh, 
     1        isig, isig1, isig3, isig5, imax, jtest, lmin, lmax

      logical ldebug, lprint, lcol
      parameter (ldebug=.false.)

C     Entry for Z production
      ibos = 2
      lcol = .false.
      goto 40
      
      entry qqZqqj_c(pbar, fsign, qbar, gsign, 
     1               uucc_c, uuss_c, ddcc_c, ddss_c, udsc_c, ducs_c)
      ibos = 2
      lcol = .true.
      goto 40

C  Entry for W+ production
      entry qqWPqqj(pbar,fsign,qbar,gsign, 
     1              uucs,ddcs,udcc,udss)
      ibos = 3
      lcol = .false.
      goto 40
      
      entry qqWPqqj_c(pbar,fsign,qbar,gsign, 
     1                uucs_c,ddcs_c,udcc_c,udss_c)
      ibos = 3
      lcol = .true.
      goto 40
      
C  Entry for W- production
      entry qqWMqqj(pbar,fsign,qbar,gsign, 
     1              uusc,ddsc,ducc,duss)
      ibos = 4
      lcol = .false.
      goto 40

      entry qqWMqqj_c(pbar,fsign,qbar,gsign, 
     1                uusc_c,ddsc_c,ducc_c,duss_c)
      ibos = 4
      lcol = .true.
      goto 40

C  Entry for A production
      entry qqAqqj(pbar,fsign,qbar,gsign, 
     1              uucc,uuss,ddcc,ddss,udsc,ducs)
      ibos = 1
      lcol = .false.
      goto 40

      entry qqAqqj_c(pbar,fsign,qbar,gsign, 
     1                uucc_c,uuss_c,ddcc_c,ddss_c,udsc_c,ducs_c)
      ibos = 1
      lcol = .true.    

 40   continue
      

c
c for initial gluon, i.e. gsign=-1, eliminate diagrams with V-->q qbar decay
      lmin = 1
      lmax = 3
      if (gsign.eq.-1) then
         if (fsign(1).eq.-fsign(2) .and. 
     1       fsign(3).eq.fsign(4) ) then
            lmax = 1            !  initial gluon attached to 1-2 line
         elseif (fsign(1).eq.fsign(2) .and. 
     1           fsign(3).eq.-fsign(4) ) then
            lmin = 3            !  initial gluon attached to 3-4 line
         endif
      endif
c fix strong coupling gs**2 for the 2 quark lines
      fpials(1) = fpi*als(1,1)
      fpials(3) = fpi*als(2,1)
c assign chirality summation range
      if (jsig.ge.1 .and. jsig.le.sigmax) then !fixed chirality selection
         jsigmin = jsig
         jsigmax = jsig
         jsig1 = sig1(jsig,3)
         jsig3 = sig3(jsig,3)
         jsig5 = sig5(jsig,3)
      else                                     !do loop over chiralities
         jsigmin = 1
         jsigmax = sigmax
         jsig1 = 0
         jsig3 = 0
         jsig5 = 0
      endif
     
c reset summation array fac to zero
      do j = 1,7                  ! all graphs
         do k = 1,nprmax          ! all flavors
            do isig = 1,sigmax    ! all helicity combinations
               do i = 1,3
                  do l = 1,3
                     fac(k,isig,j,i,l) = (0d0,0d0)
                  enddo
               enddo
            enddo
         enddo
      enddo
                
C  Define the internal momenta
      
      do mu = 0,3
         do i = 1,6                       ! the six fermions, 1, 3
            p(mu,i) = pbar(mu,i)*fsign(i) !    are incoming, 2, 4 are 
         enddo                            !    outgoing, V1 ---> f5, f6
         q(mu) = gsign*qbar(mu)      ! the gluon	
      enddo
      if (lwarn) then
         do mu=0,3
            dummy(mu) = p(mu,1)+p(mu,3)+p(mu,5)-
     1                 (p(mu,2)+p(mu,4)+p(mu,6)+q(mu))
         enddo
         if (abs(dummy(0)).gt.1d-8 .or.
     1       abs(dummy(1)).gt.1d-8 .or.
     2       abs(dummy(2)).gt.1d-8 .or.
     3       abs(dummy(3)).gt.1d-8 ) then
            print*," momentum mismatch in qqZWjqq "
            print 102," dummy = ",(dummy(mu),mu=0,3)
 101        format( " p(", i1, ") = ", 4(f10.3, 2x) )
 102        format(a,4f10.3)
            do i = 1,6
               write(6,101) i, p(0,i), p(1,i), p(2,i), p(3,i)
            end do
            write(6,101) 7, q(0), q(1), q(2), q(3)
         endif
      endif
      
      do mu = 0, 3
         p21(mu) = p(mu,2)-p(mu,1)
         p43(mu) = p(mu,4)-p(mu,3)
         p65(mu) = p(mu,6)-p(mu,5)
         p21g(mu) = p21(mu) + q(mu)   ! with gluon emission from 1-2 line
         p43g(mu) = p43(mu) + q(mu)   ! with gluon emission from 3-4 line
      enddo
      
      p21(4) = p21(0)**2-p21(1)**2-p21(2)**2-p21(3)**2 
      p43(4) = p43(0)**2-p43(1)**2-p43(2)**2-p43(3)**2
      p65(4) = p65(0)**2-p65(1)**2-p65(2)**2-p65(3)**2
      p21g(4) = p21g(0)**2-p21g(1)**2-p21g(2)**2-p21g(3)**2
      p43g(4) = p43g(0)**2-p43g(1)**2-p43g(2)**2-p43g(3)**2

      q(4) = 0.0d0
      if (lwarn) then
         lprint = (gsign.eq.1 .and. 
     1   (p21(4).gt.0.or.p21g(4).gt.0.or.p43(4).gt.0.or.p43g(4).gt.0) ) 
     2       .or. (gsign.eq.-1 .and.
     3   ((lmin.eq.1 .and. (p21g(4).gt.0 .or. p43(4).gt.0)) .or.
     4    (lmax.eq.3 .and. (p21(4).gt.0 .or. p43g(4).gt.0)) ) )
         if (lprint) then
            print*," timelike gauge boson in qqZWjqq: "
            print*," p21^2 = ",p21(4), p21g(4)
            print*," p43^2 = ",p43(4), p43g(4)
            print*," p65^2 = ",p65(4), p65(4)-xm2(3)
         endif
      endif
      
C  Get the vector boson propagator factors. For initial gluons, only
C  half the propagators will be used since the others correspond to 
C  V--> qbar q decay. To avoid zero divides when these propagators 
C  are accidentally on mass-shell, these propagators are set to zero
      
c igauge allows to select different ways to insure EW gauge invariance
c in the face of finite W and Z width effects for EW W and Z production:
c  
c   igauge = 1       overall factor scheme
c            2       complex mass scheme, m**2-->m**2-i*m*Gamma
c                    in all propagators
c            3       for ibos=2 only: naive Z Breit Wigner
c           -1       width = 0, for testing purposes only
c
c igauge is set in global.inc

c     prop21 assumes that neither V1 nor the gluon are attached to
c     the 1-2 line
      if (lmax.eq.3) then
         prop21(1) = 1/p21(4)   ! gamma
         if (abs(igauge).eq.1 .or. igauge.eq.3) then
            prop21(2) = 1/(p21(4)-xm2(2)) ! Z
            prop21(3) = 1/(p21(4)-xm2(3)) ! W+
         elseif (igauge.eq.2) then
            prop21(2) = 1/dcmplx( p21(4)-xm2(2), xmg(2) ) ! Z
            prop21(3) = 1/dcmplx(p21(4)-xm2(3), xmg(3) ) ! W+
         endif
         prop21(4) = prop21(3)  ! W-
      else
         do i = 1,4
            prop21(i) = 0
         enddo
      endif
      
c     prop43 assumes that neither V1 nor the gluon are attached to
c     the 3-4 line
      
      if (lmin.eq.1) then
         prop43(1) = 1/p43(4)                            ! gamma
         if (abs(igauge).eq.1 .or. igauge.eq.3) then
            prop43(2) = 1/(p43(4)-xm2(2))                ! Z
            prop43(3) = 1/(p43(4)-xm2(3))                ! W+
         elseif (igauge.eq.2) then
            prop43(2) = 1/dcmplx(p43(4)-xm2(2), xmg(2) ) ! Z
            prop43(3) = 1/dcmplx(p43(4)-xm2(3), xmg(3) ) ! W+
         endif
         prop43(4) = prop43(3)                           ! W-
      else
         do i = 1,4
            prop43(i) = 0
         enddo
      endif
      
c     prop21g assumes that the gluon but not V1 is attached to the
c     1-2 line
      if (lmin.eq.1) then
         prop21g(1) = 1/p21g(4)                          ! gamma
         if (abs(igauge).eq.1 .or. igauge.eq.3) then
            prop21g(2) = 1/(p21g(4)-xm2(2))              ! Z
            prop21g(3) = 1/(p21g(4)-xm2(3))              ! W+
         elseif (igauge.eq.2) then
            prop21g(2) = 1/dcmplx( p21g(4)-xm2(2), xmg(2) ) ! Z
            prop21g(3) = 1/dcmplx( p21g(4)-xm2(3), xmg(3) ) ! W+
         endif
         prop21g(4) = prop21g(3)                         ! W-
      else
         do i = 1,4
            prop21g(i) = 0
         enddo
      endif
      
c     prop43g assumes that the gluon but not V1 is attached
c     to the 3-4 line
      if (lmax.eq.3) then
         prop43g(1) = 1/p43g(4)                          ! gamma
         if (abs(igauge).eq.1 .or. igauge.eq.3) then
            prop43g(2) = 1/(p43g(4)-xm2(2))              ! Z
            prop43g(3) = 1/(p43g(4)-xm2(3))              ! W+
         elseif (igauge.eq.2) then
            prop43g(2) = 1/dcmplx(p43g(4)-xm2(2), xmg(2) ) ! Z
            prop43g(3) = 1/dcmplx(p43g(4)-xm2(3), xmg(3) ) ! W+
         endif
         prop43g(4) = prop43g(3)                         ! W-
      else
         do i = 1,4
            prop43g(i) = 0
         enddo
      endif

      if (ibos.eq.2) then
         prop65(1) = 1/p65(4)                            ! gamma
         if (igauge.ne.-1) then
            prop65(2) = 1/dcmplx(p65(4)-xm2(2),xmg(2))   ! Z
         else
            prop65(2) = 1/(p65(4)-xm2(2))               ! Z
         endif
         prop65(3) = 0                                   ! No W+
         prop65(4) = 0                                   ! No W-
C  correction factor for diagrams 5,6
         if (igauge.eq.1) then
            corr56 = (p65(4)-xm2(2))*prop65(2)
            prop65(1) = prop65(1) * corr56
         else
            corr56 = 1
         endif
      elseif (ibos.eq.1) then


         prop21(2) = 1/dcmplx(p21(4)-xm2(2),xmg(2))
         prop21(3) = 1/dcmplx(p21(4)-xm2(3),xmg(3))
         prop21(4) = prop21(3)
         
         prop43(2) = 1/dcmplx(p43(4)-xm2(2),xmg(2))
         prop43(3) = 1/dcmplx(p43(4)-xm2(3),xmg(3))
         prop43(4) = prop43(3)

         prop65(1) = 1d0        ! real gamma
         prop65(2) = 0d0        ! no Z
         prop65(3) = 0d0        ! No W+
         prop65(4) = 0d0        ! No W-

         corr56 = 0d0
      else
         if (igauge.ne.-1) then
            prop65(3) = 1/dcmplx(p65(4)-xm2(3),xmg(3))
         else
            prop65(3) = 1/(p65(4)-xm2(3))
         endif
         prop65(4) = prop65(3)
         prop65(1) = 0
         prop65(2) = 0
C  correction factor for diagrams 5,6
         if (igauge.eq.1) then
            corr56 = (p65(4)-xm2(3))*prop65(3)
         else
            corr56 = 1
         endif
      endif
      
      do vh = 1,4
         
c     The following notation is used.  propl(boson ID, (j+1)/2,
c     fermion line ID)Fermion line ID = 1 for the 1-2 line, 3 for the
c     3-4 line, identifies the line to which the gluon is attached.
         
c     For graphs 1,...,4 propl = prop65, proph ==> prop21 or prop43,
c     whichever corresponds to the fermion line NOT connected to V1.
         
c     V1 attached to 3-4 line, gluon radiation from the 1-2 line
         propl(vh,1,1) = prop65(vh)
         proph(vh,1,1) = prop21g(vh)
         
c     V1 attached to 3-4 line, gluon radiation from the 3-4 line
         propl(vh,1,3) = prop65(vh)
         proph(vh,1,3) = prop21(vh)
         
c     V1 attached to 1-2 line, gluon radiation from the 1-2 line
         propl(vh,2,1) = prop65(vh)
         proph(vh,2,1) = prop43(vh)
         
c     V1 attached to 1-2 line, gluon radiation from the 3-4 line
         propl(vh,2,3) = prop65(vh)
         proph(vh,2,3) = prop43g(vh)
         
c     For graphs 5, 6 there is no V1.  propl = prop21
c     and proph = prop43.
         
c     No V1, gluon radiation from the 1-2 line
         propl(vh,3,1) = prop21g(vh)
         proph(vh,3,1) = prop43(vh)
         
c     No V1, gluon radiation from the 3-4 line
         propl(vh,3,3) = prop21(vh)
         proph(vh,3,3) = prop43g(vh)
         
      enddo
      
c      print*,"get spinors"
C     Get the external spinors (including factor sqrt(2E) and sign
c     factor)
      if (ibos.eq.1) then 
        call psi0m(4,pbar,fsign,psi) 
      else
        call psi0m(6,pbar,fsign,psi) ! 6 external fermions with
                                             ! phsical momenta pbar
      endif
c      print*,"got spinors"
      
c     NOTE: psi(2,-1:1,6) is a two component spinor with helicity -1
c     or 1.  The last entry identifies the fermion.  If this entry is
c     odd psi is a ket, if even psi is a bra.
c     psi(1,isig1,1) = |1>_isig1
c     psi(1,isig3,3) = |3>_isig3
c     psi(1,isig5,5) = |5>_isig5
c     psi(1,isig1,2) = <2|_isig1
c     psi(1,isig3,4) = <4|_isig3
c     psi(1,isig5,6) = <6|_isig5
      
c      print*,"call curr"
C     get the f-fbar currents J21^mu, J43^mu, J65^mu with no gluon
c     radiation
c     NOTE: curr(sigmax,psibar,psi,j)
c     NOTE: j(0:3,-1:1) is a 4-vector with a helicity index
      call curr(1,psi(1,-1,2),psi(1,-1,1),j21) ! <2|gam.mu|1>_isig1
      call curr(1,psi(1,-1,4),psi(1,-1,3),j43) ! <4|gam.mu|3>_isig3
      if (ibos.eq.2) then ! if V_1 = Z both helicities possible
         call curr(1,psi(1,-1,6),psi(1,-1,5),j65)
                                               ! <6|gam.mu|5>_isig5
      elseif (ibos.eq.1) then
         call helvec(p65,1,-1,j65(0,-1)) ! helicity=-1
         j65(0,0)=1d0
         j65(1,0)=1d0
         j65(2,0)=1d0
         j65(3,0)=1d0
         call helvec(p65,1, 1,j65(0, 1)) ! helicity=+1

      else    ! only left-handed fermions couple to W
         call curr(-1,psi(1,-1,6),psi(1,-1,5),j65)
                                               ! <6|gam.mu|5>_isig5
      endif                                    ! gam.mu = dirac matrix
      
      jtest = 0
      
      if (jtest .eq. 1) then
         do isig = -1, 1, 2
            do mu = 0, 4
               j21(mu,isig) = p21(mu) * (1.0d0, 0.0d0)
            enddo
         enddo
      endif
      
      if (jtest .eq. 2) then
         do isig = -1, 1, 2
            do mu = 0, 4
               j43(mu,isig) = p43(mu) * (1.0d0, 0.0d0)
            enddo
         enddo
      endif
      
      imax = 2
c      print*,"call polvec"
      do i = 1, 2 ! 2 gluon polarizations
         call polvec(qbar,i,eps(0,i)) ! get gluon pol. vectors
      enddo
c      print*,"back from polvec"
      
c     In order to test gauge invariance.
      if (imax .eq. 3) then
         do mu = 0, 3
            eps(mu,3) = q(mu)/abs(q(0))
         enddo
      endif
      
      do i = 1, imax
         do isig = -1,1,2 ! fermion helicity = p/m

c     NOTES for bras and kets: .true. if psi is a 2-spinor of the chi
c     form as output by psi0m, .false. otherwise.  the last entry is
c     the sum of the two momenta ( p plus q) and effectively the
c     momentum of the new spinor.
            
c            print*,"call bras and kets"
            call ket2r(psi(1,isig,1),.true.,p(0,1),isig,q,eps(0,i),
     $           braketg(1,isig,1,i),pq(0,1))      ! |q,1>_i,isig1
            call bra2r(psi(1,isig,2),.true.,p(0,2),isig,q,eps(0,i),
     $           braketg(1,isig,2,i),pq(0,2))      ! <2,q|_i,isig2
            call ket2r(psi(1,isig,3),.true.,p(0,3),isig,q,eps(0,i),
     $           braketg(1,isig,3,i),pq(0,3))      ! |q,3>_i,isig3
            call bra2r(psi(1,isig,4),.true.,p(0,4),isig,q,eps(0,i),
     $           braketg(1,isig,4,i),pq(0,4))      ! <4,q|_i,isig4
         enddo
c         print*,"back from bras and kets"
         
c     Braketg contains the free spinors multiplied by a fermion
c     propagator and a gluon eps_slash.  We assume that f5 and f6 are
c     leptons, hence do not couple to gluons.
c     NOTATION: braketg(2 component spinor, sig =-1 or 1,
c     fermion ID, gluon polarization)
         
      enddo
      
      
c     Get the f-fbar currents with one gluon radiated from the
c     current line.  There are two terms, one for gluon emission to
c     either side of the ffV vertex.
      do i = 1, imax ! gluon polarizations
         call curr(1,psi(1,-1,2),braketg(1,-1,1,i),jh1)	
c                                            =  <2|vertex|q,1>_i,isig1
         call curr(1,braketg(1,-1,2,i),psi(1,-1,1),jh2)	
c                                            =  <2,q|vertex|1>_i,isig1
         do isig = -1,1,2
            do mu = 0,3
               e21(mu,isig,i) = jh1(mu,isig) + jh2(mu,isig)
c                            = (<2|gam.mu|q,1>+<2,q|gam.mu|1>)_i,isig1
            enddo
         enddo
         
         call curr(1,psi(1,-1,4),braketg(1,-1,3,i),jh1)	
c                                           =   <4|gam.mu|q,3>_i,isig3
         call curr(1,braketg(1,-1,4,i),psi(1,-1,3),jh2)	
c                                           =   <4,q|gam.mu|3>_i,isig3
         do mu = 0,3
            do isig = -1,1,2
               e43(mu,isig,i) = jh1(mu,isig) + jh2(mu,isig)
c                            = (<4|gam.mu|q,3>+<4,q|gam.mu|3>)_i,isig3
            enddo
         enddo
      enddo
      
      if (jtest .eq. 3) then
         do isig = -1, 1, 2
            do mu = 0, 4
               do i = 1, 2
                  e21(mu,isig,i) = p21g(mu) * (1.0d0, 0.0d0)
               enddo
            enddo
         enddo
      endif
      
      if (jtest .eq. 4) then
         do isig = -1, 1, 2
            do mu = 0, 4
               do i = 1, 2
                  e43(mu,isig,i) = p43g(mu) * (1.0d0, 0.0d0)
               enddo
            enddo
         enddo
      endif
      
c      print*,"get the amplitudes of the first 6 Feynman graphs"
C     get the amplitudes of the first 6 Feynman graphs for all
c     allowed helicity combinations.  matr(1) ==> gluon radiation
c     from 1-2 line, matr(3) ==> gluon radiation from 3-4 line.
      do j = 1,6 ! Feynman graphs w/o triple vertex
         if(ibos.eq.1.and.j.eq.5) goto 55
         jj = (j+1)/2                           ! pairs 1&2, 3&4, 5&6
         do isig = jsigmin,jsigmax              ! helicity combinations
            isig1 = sig1(isig,jj)
            if (isig1.eq.0) goto 50   ! some combinations not possible
            isig3 = sig3(isig,jj)
            isig5 = sig5(isig,jj)
            do i = 1, imax ! gluon polarizations
               
C     graphs a and b (1,2)
               if (j.eq.1) then ! attach V1 to f4
                  
c     gluon radiation from 1-2 fermion line
                  call bra2c(psi(1,isig3,4),.true.,p(0,4),isig3,p65,
     $                 j65(0,isig5),braketV,pk) ! <4,p65|_isig3,isig5
                  matr1 = -s1c(braketV,e21(0,isig1,i),.true.,isig3,
     $                 psi(1,isig3,3))
c    			=  <4,p65|e21|3>_i,isig1,isig3,isig5
c			=  <4,p65|e21|3>_i,isig
c     NOTE: isig = ( isig1,isig2,isig3 ), i is the gluon index
                  
c     gluon radiation from 3-4 fermion line
                  call bra2c(braketg(1,isig3,4,i),.false.,pq(0,4),
     $                 isig3,p65,j65(0,isig5),braketgV,dummy)
                                           ! <4,q,p65|_i,isig3,isig5
                  call bra2r(braketV,.false.,pk,isig3,q,eps(0,i),
     $                 braketVg,dummy)     ! <4,p65,q|_i,isig3,isig5
                  matr2 = -s1c(braketgV,j21(0,isig1),.true.,isig3,
     $                 psi(1,isig3,3))
c                                           = <4,q,p65|j21|3>_i,isig
                  matr3 = -s1c(braketVg,j21(0,isig1),.true.,isig3,
     $                 psi(1,isig3,3))
c                                           = <4,p65,q|j21|3>_i,isig
                  matr4 = -s1c(braketV,j21(0,isig1),.true.,isig3,
     $                 braketg(1,isig3,3,i))
c                                           = <4,p65|j21|q,3>_i,isig
                  
                  matr(1) = matr1 ! gluon attached to 1-2 line
                  matr(3) = matr2 + matr3 + matr4
                                        ! gluon attached to 3-4 line
                  
               elseif (j.eq.2) then ! attach V1 to f3
                  
c     gluon radiation from the 1-2 line
                  call ket2c(psi(1,isig3,3),.true.,p(0,3),isig3,p65,
     $                 j65(0,isig5),braketV,pk)
                                               ! |p65,3>_isig3,isig5
                  matr1 = -s1c(psi(1,isig3,4),e21(0,isig1,i),.true.,
     $                 isig3,braketV)
c                                              = <4|e21|p65,3>_i,isig
                  
                  
c     gluon radiation from the 3-4 line
                  call ket2c(braketg(1,isig3,3,i),.false.,pq(0,3),
     $                 isig3,p65,j65(0,isig5),braketgV,dummy)
                                           ! |3,q,p65>_i,isig3,isig5
                  call ket2r(braketV,.false.,pk,isig3,q,eps(0,i),
     $                 braketVg,dummy)
                                           ! |3,p65,q>_i,isig3,isig5
                  matr2 = -s1c(psi(1,isig3,4),j21(0,isig1),.true.,
     $                 isig3,braketgV)
c                                           = <4|j21|p65,q,3>_i,isig
                  matr3 = -s1c(psi(1,isig3,4),j21(0,isig1),.true.,
     $                 isig3,braketVg)
c                                           = <4|j21|q,p65,3>_i,isig
                  matr4 = -s1c(braketg(1,isig3,4,i),j21(0,isig1),
     $                 .true.,isig3,braketV)
c                                           = <4,q|j21|p65,3>_i,isig
                  
                  matr(1) = matr1       ! gluon attached to 1-2 line
                  matr(3) = matr2 + matr3 + matr4
                                        ! gluon attached to 3-4 line
                  
C     graphs c and d (3,4)
               elseif (j.eq.3) then             ! attach V1 to f2
                  
c     gluon radiation from 3-4 fermion line
                  call bra2c(psi(1,isig1,2),.true.,p(0,2),isig1,p65,
     $                 j65(0,isig5),braketV,pk) ! <2,p65|_isig1,isig5
                  matr1 = -s1c(braketV,e43(0,isig3,i),.true.,isig1,
     $                 psi(1,isig1,1))
c                                              = <2,p65|e43|1>_i,isig
                  
c     gluon radiation from 1-2 line
                  call bra2c(braketg(1,isig1,2,i),.false.,pq(0,2),
     $                 isig1,p65,j65(0,isig5),braketgV,dummy)
                                           ! <2,q,p65|_i,isig1,isig5
                  call bra2r(braketV,.false.,pk,isig1,q,eps(0,i),
     $                 braketVg,dummy)
                                           ! <2,p65,q|_i,isig1,isig5
                  matr2 = -s1c(braketgV,j43(0,isig3),.true.,isig1,
     $                 psi(1,isig1,1))
c                                           = <2,q,p65|j43|1>_i,isig
                  matr3 = -s1c(braketVg,j43(0,isig3),.true.,isig1,
     $                 psi(1,isig1,1))
c                                           = <2,p65,q|j43|1>_i,isig
                  matr4 = -s1c(braketV,j43(0,isig3),.true.,isig1,
     $                 braketg(1,isig1,1,i))
c                                           = <2,p65|j43|q,1>_i,isig
                  
                  matr(1) = matr2 + matr3 + matr4
                                        ! gluon attached to 1-2 line
                  matr(3) = matr1
                                        ! gluon attached to 3-4 line
                  
               elseif (j.eq.4) then          ! attach V1 to f1
                  
c     gluon radiation from the 3-4 line
                  call ket2c(psi(1,isig1,1),.true.,p(0,1),isig1,p65,
     $                 j65(0,isig5),braketV,pk) ! |1,p65>_isig1,isig5
                  matr1 = -s1c(psi(1,isig1,2),e43(0,isig3,i),.true.,
     $                 isig1,braketV)
c                                             = <2|e43|p65,1>_i,isig
                  
c     gluon radiation from the 1-2 line
                  call ket2c(braketg(1,isig1,1,i),.false.,pq(0,1),
     $                 isig1,p65,j65(0,isig5),braketgV,dummy)
                                            ! |p65,q,1>_i,isig1,isig5
                  call ket2r(braketV,.false.,pk,isig1,q,eps(0,i),
     $                 braketVg,dummy)      ! |q,p65,1>_i,isig1,isig5
                  matr2 = -s1c(psi(1,isig1,2),j43(0,isig3),.true.,
     $                 isig1,braketgV)
c                                           =  <2|j43|p65,q,1>_i,isig
                  matr3 = -s1c(psi(1,isig1,2),j43(0,isig3),.true.,
     $                 isig1,braketVg)
c                                            = <2|j43|q,p65,1>_i,isig
                  matr4 = -s1c(braketg(1,isig1,2,i),j43(0,isig3),
     $                 .true.,isig1,braketV)
c                                            = <2,q|j43|p65,1>_i,isig
                  
                  matr(1) = matr2 + matr3 + matr4
                                         ! gluon attached to 1-2 line
                  matr(3) = matr1        ! gluon attached to 3-4 line
                  
C     graphs e and f (5,6)
               elseif (j.eq.5) then                       ! No V1
                  
c     gluon radiation from line 1-2 line
                  call bra2c(psi(1,isig5,6),.true.,p(0,6),isig5,p43,
     $                 j43(0,isig3),braket,dummy)
                                               ! <6,p43|_isig3,isig5
                  matr(1) = -s1c(braket,e21(0,isig1,i),.true.,isig5,
     $                 psi(1,isig5,5))
                  matr(1) = matr(1)*corr56
c                                             = <6,p43|e21|5>_i,isig
                  
c     gluon radiation from 3-4 line
                  call bra2c(psi(1,isig5,6),.true.,p(0,6),isig5,p43g,
     $                 e43(0,isig3,i), braket,dummy)
                                            ! <6,p43g|_i,isig3,isig5
                  matr(3) = -s1c(braket,j21(0,isig1),.true.,isig5,
     $                 psi(1,isig5,5))
                  matr(3) = matr(3)*corr56
c                                            = <6,p43g|j21|5>_i,isig
                  
               elseif (j.eq.6) then                   ! No V1
                  
c     gluon radiation from 1-2 line
                  call ket2c(psi(1,isig5,5),.true.,p(0,5),isig5,p43,
     $                 j43(0,isig3),braket,dummy)
                                               ! |p43,5>_isig3,isig5
                  matr(1) = -s1c(psi(1,isig5,6),e21(0,isig1,i),
     $                 .true.,isig5,braket)
                  matr(1) = matr(1)*corr56
c                                             = <6|e21|p43,5>_i,isig
                  
                  
c     gluon radiation from 3-4 line
                  call ket2c(psi(1,isig5,5),.true.,p(0,5),isig5,p43g,
     $                 e43(0,isig3,i),braket,dummy)
                                              ! |p43g,5>_isig3,isig5
                  matr(3) = -s1c(psi(1,isig5,6),j21(0,isig1),.true.,
     $                 isig5,braket)
                  matr(3) = matr(3)*corr56
c                                             = <6|j21|p43g,5>_i,isig
                  
               endif

c               print*,"here is a big loop"
C     now get the coupling*propagator factors for subprocess k,
c     helicity combination isig, Feynman graph j
               do k = 1,nprmax ! flavors of quarks
                  do l = lmin, lmax, 2 ! gluon attached to 1-2 or 3-4
                     fac(k,isig,j,i,l) = 0
                     
c     NOTATION: fac(flavor combo, helicity combo, Feynman graph, gluon
c     pol., l = 1 if gluon attached to 1-2 line or 3 for 3-4 line)
                     
                     do vh = vhmin(k,j),vhmax(k,j)
                        do vl = vlmin(k,j),vlmax(k,j)
                           fac(k,isig,j,i,l) = fac(k,isig,j,i,l) +
     $                          fcpl(vl,vh,isig,k,j)
     $                          *propl(vl,jj,l)*proph(vh,jj,l)
                        enddo
                     enddo
                     fac(k,isig,j,i,l) = fac(k,isig,j,i,l)*matr(l)

                  enddo
               enddo
            enddo
 50         continue
         enddo
      enddo

 55   continue
      
C     and finally the three boson vertex. An extra factor of -1 is
c     included because OUTGOING momenta p21,p43,p65 are used in the
c     call of tbv
C     The correspondence of the boson identifiers V_i and the
c     currents J_ij is
C     V_1  =  J_65
C     V_2  =  J_21
C     V_3  =  J_43
      
      do isig = 1, 8                  ! all helicity combinations
         do k = 1, nprmax             ! all flavor combinations
            do i = 1, imax            ! all gluon polarizations
               do l = 1, 3, 2         ! gluon attached to 1-2 or 3-4
                  fac(k,isig,7,i,l) = 0
               enddo
            enddo
         enddo
      enddo
      
      if (ibos.eq.2.or.ibos.eq.1) then             ! Z production via WWZ vertex
                                      ! W+/W- attached to 1-2 or 3-4
         do isig = max(jsigmin,1),min(jsigmax,2),1
            isig1 = sig1(isig,4)
            isig3 = sig3(isig,4)
            isig5 = sig5(isig,4)
            
C     process k = 5 has V2 = W+, V3 = W-
                     
            do i = 1, imax                   ! gluon polarizations
               call tbv(e21(0,isig1,i),j43(0,isig3),j65(0,isig5),
     $              p21g,p43,p65,matr7(1,i,1))
                        
               call tbv(j21(0,isig1),e43(0,isig3,i),j65(0,isig5),
     $              p21,p43g,p65,matr7(1,i,3))
               
               do  jj = 1, 2
                  fac(5,isig,7,i,1) = fac(5,isig,7,i,1)
     $                 + matr7(jj,i,1)*prop65(jj)*fcpl(jj,1,isig,5,7)
     $                 * prop21g(3) * prop43(4)
                  fac(5,isig,7,i,3) = fac(5,isig,7,i,3)
     $                 + matr7(jj,i,3)*prop65(jj)*fcpl(jj,1,isig,5,7)
     $                 * prop21(3) * prop43g(4)
               enddo
            enddo
            
C     process 6 has V2 = W-, V3 = W+
                     
            do i = 1, imax                     ! gluon polarizations
               call tbv(j43(0,isig3),e21(0,isig1,i),j65(0,isig5),
     $              p43,p21g,p65,matr7(1,i,1))
                        
               call tbv(e43(0,isig3,i),j21(0,isig1),j65(0,isig5),
     $              p43g,p21,p65,matr7(1,i,3))
                        
               do jj = 1, 2 ! 1 = gamma, 2 = Z
                  fac(6,isig,7,i,1) = fac(6,isig,7,i,1)
     $                 + matr7(jj,i,1)*prop65(jj)*fcpl(jj,1,isig,6,7)
     $                 * prop21g(4) * prop43(3)
                  fac(6,isig,7,i,3) = fac(6,isig,7,i,3)
     $                 + matr7(jj,i,3)*prop65(jj)*fcpl(jj,1,isig,6,7)
     $                 * prop21(4) * prop43g(3)
               enddo
            enddo
         enddo
         
      elseif (ibos.eq.3) then                       ! W+ production
         do k = 1,4                                 ! 4 flavor combos.
            do isig = jsigmin,jsigmax
               isig1 = sig1(isig,4)
               isig3 = sig3(isig,4)
               isig5 = sig5(isig,4)
               
               do i = 1, imax                   ! gluon polarizations
                  if (k.le.2 .and. isig.ne.3) then
                                                 ! non-vanishing only
c     processes 1 and 2 have V1 = W-, V2 = Z/gamma, V3 = W+
                     
                     call tbv(e43(0,isig3,i),j65(0,isig5),
     $                    j21(0,isig1),p43g,p65,p21,matr7(1,i,3))
                     call tbv(j43(0,isig3),j65(0,isig5),
     $                    e21(0,isig1,i),p43,p65,p21g,matr7(1,i,1))
                     
                     do jj = 1, 2
                        fac(k,isig,7,i,1) = fac(k,isig,7,i,1)
     $                       + matr7(jj,i,1) * fcpl(jj,1,isig,k,7)
     $                       * prop65(3) * prop21g(jj) * prop43(4)
                        fac(k,isig,7,i,3) = fac(k,isig,7,i,3)
     $                       + matr7(jj,i,3) * fcpl(jj,1,isig,k,7)
     $                       * prop65(3) * prop21(jj) * prop43g(4)
                     enddo
                     
                  elseif (k.gt.2 .and. isig.ne.2) then
                                                  ! non-vanishing only
c     processes 3 and 4 have V1 = W-, V2 = W+, V3 = Z/gamma
                     
                     call tbv(e21(0,isig1,i),j65(0,isig5),
     $                    j43(0,isig3),p21g,p65,p43,matr7(1,i,1))
                     call tbv(j21(0,isig1),j65(0,isig5),
     $                    e43(0,isig3,i),p21,p65,p43g,matr7(1,i,3))
                     
                     do jj = 1, 2
                        fac(k,isig,7,i,1) = fac(k,isig,7,i,1)
     $                       + matr7(jj,i,1) * fcpl(jj,1,isig,k,7)
     $                       * prop65(3) * prop21g(4) * prop43(jj)
                        fac(k,isig,7,i,3) = fac(k,isig,7,i,3)
     $                       + matr7(jj,i,3) * fcpl(jj,1,isig,k,7)
     $                       * prop65(3) * prop21(4) * prop43g(jj)
                     enddo
                  endif
               enddo               
            enddo
         enddo
         
      elseif (ibos.eq.4) then                        ! W- production
         do k = 1,4                                 ! 4 flavor combos.
            do isig = jsigmin,jsigmax
               isig1 = sig1(isig,4)
               isig3 = sig3(isig,4)
               isig5 = sig5(isig,4)
               
               do i = 1, imax                   ! gluon polarizations
                  
                  if (k.le.2 .and. isig.ne.3) then
                                                 ! non-vanishing only
c     processes 1 and 2 have V1 = W+, V2 = Z/gamma, V3 = W-
                     
                     call tbv(j65(0,isig5),j43(0,isig3),
     $                    e21(0,isig1,i),p65,p43,p21g,matr7(1,i,1))
                     call tbv(j65(0,isig5),e43(0,isig3,i),
     $                    j21(0,isig1),p65,p43g,p21,matr7(1,i,3))
                     
                     do jj = 1, 2
                        fac(k,isig,7,i,1) = fac(k,isig,7,i,1)
     $                       + matr7(jj,i,1) * fcpl(jj,1,isig,k,7)
     $                       * prop65(4) * prop21g(jj) * prop43(3)
                        fac(k,isig,7,i,3) = fac(k,isig,7,i,3)
     $                       + matr7(jj,i,3) * fcpl(jj,1,isig,k,7)
     $                       * prop65(4) * prop21(jj) * prop43g(3)
                     enddo
                     
                  elseif (k.gt.2 .and. isig.ne.2) then
                                                 ! non-vanishing only
c     processes 3 and 4 have V1 = W+, V2 = W-, V3 = Z/gamma
                     
                     call tbv(j65(0,isig5),e21(0,isig1,i),
     $                    j43(0,isig3),p65,p21g,p43,matr7(1,i,1))
                     call tbv(j65(0,isig5),j21(0,isig1),
     $                    e43(0,isig3,i),p65,p21,p43g,matr7(1,i,3))
                     
                     do jj = 1, 2
                        fac(k,isig,7,i,1) = fac(k,isig,7,i,1)
     $                       + matr7(jj,i,1) * fcpl(jj,1,isig,k,7)
     $                       * prop65(4) * prop21g(3) * prop43(jj)
                        fac(k,isig,7,i,3) = fac(k,isig,7,i,3)
     $                       + matr7(jj,i,3) * fcpl(jj,1,isig,k,7)
     $                       * prop65(4) * prop21(3) * prop43g(jj)
                     enddo
                  endif
               enddo
            enddo
         enddo
         
      endif
      
      do k = 1, nprmax
         do isig = 1, sigmax
            do i = 1, 3
               do l = 1, 3
                  mm4(k,isig,i,l) = 0.0d0
               enddo
               mm3(k,isig,i) = 0.0d0
            enddo
            mm2(k,isig) = 0.0d0
         enddo
      enddo


c     Now sum the contributions from the seven Feynman graphs for
c     all subprocesses
       
      if (jtest .eq. 1) then
         print*,"answer answer answer answer answer answer answer "
         print*,(fac(1,1,1,1,3) + fac(1,1,2,1,3) + fac(1,1,5,1,3) +
     $        fac(1,1,6,1,3))/fac(1,1,1,1,3) 
         print*,"answer answer answer answer answer answer answer "
      endif
      
      if (jtest .eq. 2) then
         print*,"answer answer answer answer answer answer answer "
         print*,(fac(1,1,3,1,1) + fac(1,1,4,1,1) + fac(1,1,5,1,1)
     $        + fac(1,1,6,1,1))/fac(1,1,3,1,1) 
         print*,"answer answer answer answer answer answer answer "
      endif
      
      if (jtest .eq. 3) then
         print*,"answer answer answer answer answer answer answer "
         print*,(fac(1,1,1,1,1) + fac(1,1,2,1,1) + fac(1,1,5,1,1)
     $        + fac(1,1,6,1,1))/fac(1,1,1,1,1) 
         print*,"answer answer answer answer answer answer answer "
      endif
      
      if (jtest .eq. 4) then
         print*,"answer answer answer answer answer answer answer "
         print*,(fac(1,1,3,1,3) + fac(1,1,4,1,3) + fac(1,1,5,1,3)
     $        + fac(1,1,6,1,3))/fac(1,1,3,1,3) 
         print*,"answer answer answer answer answer answer answer "
      endif
 
      if (ldebug) print*," sigmax = ",sigmax
      do k = 1, nprmax
         do isig5 = -1,1
            res(k,1,isig5) = 0
            res(k,2,isig5) = 0
         enddo
         do isig = jsigmin,jsigmax
            isig5 = sig5(isig,3)
            if (ldebug) then
               if( (sig5(isig,1).ne.isig5 .and. sig5(isig,1).ne.0) .or.
     1             (sig5(isig,2).ne.isig5 .and. sig5(isig,2).ne.0) .or.
     2             (sig5(isig,4).ne.isig5 .and. sig5(isig,4).ne.0) ) 
     3           then
                  print*," check isig5 assignment for isig = ",isig
                  print*," sig5 = ",(sig5(isig,ic),ic=1,4)
               endif
            endif
            do i = 1, imax
               do l = lmin, lmax, 2
                  do j = 1, 7
                     mm4(k,isig,i,l) = mm4(k,isig,i,l) +
     &                                 fac(k,isig,j,i,l)
                  enddo
                  res(k,(l+1)/2,isig5) = res(k,(l+1)/2,isig5) + 
     &                 fpials(l) * (dreal(mm4(k,isig,i,l))**2 +
     &                              dimag(mm4(k,isig,i,l))**2)
               enddo
            enddo
         enddo
         do isig5 = -1,1,2            !   * color factor
            res(k,1,isig5) = res(k,1,isig5) * 12.0d0
            res(k,2,isig5) = res(k,2,isig5) * 12.0d0
         enddo
         res(k,1,0) = res(k,1,-1)+res(k,1,1)
         res(k,2,0) = res(k,2,-1)+res(k,2,1)
      enddo
c     NOTE:  It is neccesary to square before summing over the line
c     index (do l = 1, 3, 2) because there is no interference between
c     radiation from the upper line and radiation from the lower line.

      if (imax .eq. 3) then
         print*,"gauge invariance test"
         do k = 1, nprmax
            do isig = jsigmin,jsigmax
               do l = lmin, lmax, 2
                  do i = 1, 3
                     print*,"k,isig,l,i"
                     print*,k,isig,l,i
c                      print*,mm4(k,isig,i,l),mm4(k,isig,1,l)
                     print*,mm4(k,isig,i,l)
     &                   /( mm4(k,isig,1,l) + dcmplx(1.d-30,1.d-30) )
                  enddo
               enddo
            enddo
         enddo
      endif

      if (ldebug) then
         print*," gs = ",b(5,5,5)," alfas = ",b(5,5,5)**2/4/3.14159
         do k = 1,4
            print*,k," res = ",res(k,1,0),res(k,2,0)
         enddo
      endif
c  set processes to zero if photon virtuality falls below cutoff
      if (ibos .eq.2 .or.ibos.eq.1) then
         if (   abs(p21(4)) .lt.qsqAmin 
c     1     .or. abs(p43g(4)).lt.qsqAmin 
     2       ) then
            do k = 1,4
               res(k,2,0) = 0   ! no graphs with gluon on lower line
            enddo
         endif
         if (    abs(p43(4)) .lt.qsqAmin 
c     1      .or. abs(p21g(4)).lt.qsqAmin 
     2       ) then
            do k = 1,4
               res(k,1,0) = 0   ! no graphs with gluon on upper line
            enddo
         endif
      else
         if (abs(p21(4)) .lt.qsqAmin ) then
            res(1,2,0) = 0
            res(2,2,0) = 0
         endif
         if (abs(p43g(4)).lt.qsqAmin ) then
            res(3,2,0) = 0
            res(4,2,0) = 0
            if (lwarn) then
               print*," WARNING from qqZWjqqi.f: "
               print*," photon virtuality p43g = ",p43g(4),
     1             " below cutoff +- ",qsqAmin
            endif
         endif
         if (abs(p21g(4)).lt.qsqAmin ) then
            res(1,1,0) = 0
            res(2,1,0) = 0
            if (lwarn) then
               print*," WARNING from qqZWjqqi.f: "
               print*," photon virtuality p21g = ",p21g(4),
     1             " below cutoff +- ",qsqAmin
            endif
         endif
         if (abs(p43(4)) .lt.qsqAmin ) then
            res(3,1,0) = 0
            res(4,1,0) = 0
         endif
      endif

      if (ldoblha) then
        if (blha_amptype(blha_curproc).eq.0) then ! tree
          if ( (blha_idsubproc(blha_cursubproc,blha_curproc)/10 .eq. 1) .or.
     &         (blha_idsubproc(blha_cursubproc,blha_curproc)/10 .eq. 2) ) then
            blha_amp(1)=blha_amp(1)+blha_CF/(4d0/3d0)*
     &       res(mod(blha_idsubproc(blha_cursubproc,blha_curproc),10),1,0)
          endif
          if ( (blha_idsubproc(blha_cursubproc,blha_curproc)/10 .eq. 1) .or.
     &         (blha_idsubproc(blha_cursubproc,blha_curproc)/10 .eq. 3) ) then
            blha_amp(1)=blha_amp(1)+blha_CF/(4d0/3d0)*
     &       res(mod(blha_idsubproc(blha_cursubproc,blha_curproc),10),2,0)
          endif

        else 
          call BLHA_amptypeerror(blha_amptype(blha_curproc),__FILE__,__LINE__)
        endif
      endif

C  Assign the final results for the subprocesses
      if ((ibos.eq.2.or.ibos.eq.1) .and. .not.lcol) then
         uucc = res(1,1,0) +  res(1,2,0)
         uuss = res(2,1,0) +  res(2,2,0)
         ddcc = res(3,1,0) +  res(3,2,0)
         ddss = res(4,1,0) +  res(4,2,0)
         udsc = res(5,1,0) +  res(5,2,0)
         ducs = res(6,1,0) +  res(6,2,0)
      elseif ((ibos.eq.2.or.ibos.eq.1) .and. lcol) then
         do i = 1,2
            uucc_c(i+1) = res(1,i,0) 
            uuss_c(i+1) = res(2,i,0)
            ddcc_c(i+1) = res(3,i,0)
            ddss_c(i+1) = res(4,i,0)
            udsc_c(i+1) = res(5,i,0)
            ducs_c(i+1) = res(6,i,0)
         enddo
         uucc_c(1) = res(1,1,0) +  res(1,2,0)
         uuss_c(1) = res(2,1,0) +  res(2,2,0)
         ddcc_c(1) = res(3,1,0) +  res(3,2,0)
         ddss_c(1) = res(4,1,0) +  res(4,2,0)
         udsc_c(1) = res(5,1,0) +  res(5,2,0)
         ducs_c(1) = res(6,1,0) +  res(6,2,0)
      elseif (ibos.eq.3 .and. .not.lcol) then
         uucs = res(1,1,0)+res(1,2,0) 
         ddcs = res(2,1,0)+res(2,2,0)
         udcc = res(3,1,0)+res(3,2,0)
         udss = res(4,1,0)+res(4,2,0)
      elseif (ibos.eq.3 .and. lcol) then
         do i = 1,2
            uucs_c(i+1) = res(1,i,0) 
            ddcs_c(i+1) = res(2,i,0)
            udcc_c(i+1) = res(3,i,0)
            udss_c(i+1) = res(4,i,0)
         enddo
         uucs_c(1) = res(1,1,0)+res(1,2,0) 
         ddcs_c(1) = res(2,1,0)+res(2,2,0)
         udcc_c(1) = res(3,1,0)+res(3,2,0)
         udss_c(1) = res(4,1,0)+res(4,2,0)
      elseif (ibos.eq.4 .and. .not.lcol) then
         uusc = res(1,1,0)+res(1,2,0)
         ddsc = res(2,1,0)+res(2,2,0)
         ducc = res(3,1,0)+res(3,2,0)
         duss = res(4,1,0)+res(4,2,0)
      elseif (ibos.eq.4 .and. lcol) then
         do i = 1,2
            uusc_c(i+1) = res(1,i,0) 
            ddsc_c(i+1) = res(2,i,0)
            ducc_c(i+1) = res(3,i,0)
            duss_c(i+1) = res(4,i,0)
         enddo
         uusc_c(1) = res(1,1,0)+res(1,2,0)
         ddsc_c(1) = res(2,1,0)+res(2,2,0)
         ducc_c(1) = res(3,1,0)+res(3,2,0)
         duss_c(1) = res(4,1,0)+res(4,2,0)
      endif
                                                               
      end
      
