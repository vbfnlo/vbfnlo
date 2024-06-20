!****************subroutine QQHWJ *************************************
c
c        giuseppe bozzi, <giuseppe@particle.uni-karlsruhe.de>
c        Last modified: December 2009
C
c       Francisco Campanario, <francam@particle.uni-karlsruhe.de>
c       Last modified: Jul 2010
C       Mod: including NLO: hexagons,pentagons, boxes and different 
C      checks, including gauge and factorization of the IR 
C      divergences.....
c
c
c      
c    
c
c
C
C  QQVHj calculates the matrix elements**2 for 
C     q1 q2    ---->   W+- A A g,   
C     W+- ---> f1-bar f2
C
C  and crossing related processes. Pauli interference terms for
c  identical fermions are neglected.
c
c  Gluon polarization vectors are given in the kartesian basis
c  (comparison with madgraph amplitudes tampr requires to express these in
c   kart. basis also, even though they are given in the helicity basis
c   a priori).
c
c  This code includes only real emission contributions, i.e.
c
c      return udwwz = |M_real|^2   etc.
c
c       fpials is attached only in the end of the code
c
c       glupol is the gluon polarization in the kartesian basis (glupol=1,2)
c       k is the process ID (1:uuwwaj,2:ddwwaj)
c       isig is the helicity of partons 1 and 2
c
c*********************************************************************

      subroutine qqHWj(pbar, fsign, qbar, gsign, nlo, N, final, sbos)

      use globalvars, only: ldoblha

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

c     electroweak couplings are taken from KOPPLN

      double precision  clr, xm2, xmg, b, v, a
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     1                  V(4,5),A(4,5)

c variables for the main part of the program      

      DOUBLE PRECISION fpials(2:3), fpi
      PARAMETER (fpi=4d0*pi)

      ! pbar needs to hold 6 to allow most general calling case (WHjj)
      ! TODO: change 5 to either include or exclude all leptons
      double precision  pbar(0:3,6), qbar(0:4), q(0:3)
      double precision  res(2), resv(2), final
      double precision  p(0:3,6)

      integer fsign(6), gsign, mu, i, j, k, isig, nlo, N
      integer lh, glupol, sbos, bos

      double precision pq(0:4,2)

      double precision eps(0:3,2) ! external gluon

      ! quarks
      double complex psi(2,-1:1,2)
      double complex braketg(2,-1:1,2,2)
      double complex mat(2,-1:1,2,3), madgr(3,-1:1,-1:1,2),
     1               matb(2,-1:1,2,3)
      double complex matv(2,-1:1,2,4)
      double complex mm(2, -1:1,2), mv(2, -1:1,2)
      double complex mmb(2, -1:1,2)
      double complex m1kb(2,-1:1,3), m2kb(2,-1:1,3), temp, temp1
      double complex ma
      double complex jqgwhq(0:5,-1:1,2), jqwhgq(0:5,-1:1,2)
      double complex jqqwh(0:5,-1:1,2)
   
      real*8 dotrr,theta
      double complex im
      double complex dotcc, s1c, dotrc,cnf
      double complex zero 
      external dotcc, s1c, dotrc,dotrr,theta
      save k
      parameter (im = (0d0,1d0))
      parameter (zero = (0d0,0d0) )

c variable for q_mu q_nu terms in the propagators of the massive vector bosons

      double complex zm2i(2:3)
      save zm2i

      logical ldebug, ldebugm, linit
      logical ldebugtop, ldebugmtop
      logical lnlo, lbox
      data linit /.true./
      save linit
      parameter(ldebugm   = .false.)  !run madgraph checks
      parameter(ldebug    = .false.) !print debug output
      parameter(ldebugtop = .false.) !run toploop with debug info
      parameter(ldebugmtop= .false.)!do madgraph check with toploops
      real*8 tmp ! for madgraph result

c variables for virtual corrections
      DOUBLE PRECISION murbox
      double precision c2,c2o4pi,oneo4pi
      parameter (c2=4d0/3d0, c2o4pi=c2/4d0/pi, oneo4pi=1d0/4d0/pi)
c Variables for boxes
      complex*16 j5g(5,2,2),jbg(2,2),jg5n(2)
      complex*16 j5gg(5,2,2),jbgg(2,2),jg5ng(2)
      complex*16 j5gg2(5,2,2),jbgg2(2,2),jg5ng2(2)
      complex*16 nothing(2)
      complex*16 boxcontr1(2),boxcontr1b(2)
      complex*16 boxcontr2(2),boxcontr2b(2)
      complex*16 boxcontrNoAbe3(2)
      real*8 mqwh(0:3),td
      complex*16 wwh(0:3,3)
      complex*16 emqwh(0:3),empg(0:3)
c Variables for pentagons      
      real*8 p1(0:3),p2(0:3)
      real*8 mpg(0:3)

c Variables for hexagons
      Integer Div,PerformGauge,mutemp
      complex*16 mv1CA,mv1b
      complex*16 lns,lnt,lnu,cteEps1

      Logical wardidtest1,wardidtest2,wardidtest3
c  Colors Parameters       
      real*8 NF,inv2,inv3,CA,TR,tf,CF,eulergamma
!      parameter (NF=4d0) ! nf gets set using the global nfl
      parameter (inv2=1d0/2d0,inv3=1d0/3d0)
      parameter (CA=3d0)
      parameter (TR=1d0/2d0,tf=1d0/2d0)
      parameter (CF=4d0/3d0)
      parameter (eulergamma=5.772156649015328D-1)
      complex*16 cfc,cac,cfcac,pic2,cacgraph,convfact
      real*8 mans,mant,manu,musq
c     Variables for virtual
      complex*16 mv1,mv2
      complex*16 mmv1,mmv2,ratioFACT,mmv1b
      Complex*16 mmv1CA
      !complex*16 Myeps(0:3,2) ! g in kartesian basis 
      complex*16 Myeps(0:5,2) ! g in kartesian basis 
c  Compare born from penlines, hexline,boxline 
C against NLO_real
      real*8 ratioComp
      complex*16 K_cont,K_contcheck,K_contcheck1,K_contcheck2,K_contcheck3,K_contcheck4
      integer Reno_inv
C For adding Christoph factor from cancellation of 1/eps  divergences

      complex*16 factortoadd
      real*8 s,t,mmu,mursq1
      complex*16 check(2),check1(2)

      double precision RandomNumber
      external RandomNumber

      ! top Loop contributions
      logical doTopLoop
      parameter (doTopLoop = .true.)
      double complex epsgh(0:3,2) ! gluon coupling to quark line
      double complex epsext(0:3,2) ! gluon coupling to quark line
      ! those two are different if there is a top-loop radiating off a higgs


      double precision fq0(0:4,2)
      logical ltoploop
      double complex jqqw(0:5,-1:1)
      double precision pqqw(0:4)
      integer flgauge
      double complex psiw0(2,-1:1, 2)
      integer maxfloopmasses
      parameter (maxfloopmasses=1) !20)
      double precision mf(maxfloopmasses)
      double complex rfvert(2,maxfloopmasses,2)
      double precision vv
      double precision qin(0:3), qhin(0:3)
      integer DivTop, PerformGaugeTop

      integer glupolchoice

!      double complex jqqw1(0:5,-1:1)
!      double complex jqqw2(0:5,-1:1)
!      double precision pqqw1(0:4)
!      double precision pqqw2(0:4)


c for debugging purposes
      double precision momsum(0:3)
      real*8 mass2

      !BLHA interface
      double precision tree(2) ! save tree amplitude
      double precision lnmusqp1p2 ! running logs
              
c -------------- BEGIN CODE --------------------------------

c initialize & precompute stuff needed below:

c fix strong coupling gs**2 for the two quarks:

      fpials(2) = fpi*als(1,N)
      fpials(3) = fpi*als(2,N)
c      fpials(2) = fpi*alfas

      ! get number of flavors from global variable
      nf = dble(nfl)

c---------------------------- BEGIN CODE -----------------------

      lnlo = NLO.ne.0       ! include some virtual stuff if T
      lbox = (NLO.eq.-4) ! only needed if mapping wrong

!     lbox = (NLO.eq.-4).or.(NLO.eq.1).or.(NLO.eq.5).or.(NLO.eq.7).or.(NLO.eq.-7)  ! call boxline if T
!      lpt = (NLO.eq.1).or.(NLO.eq.-5).or.(NLO.eq.7).or.(NLO.eq.-7)    ! true pentagon contribution
!      !lpq = (NLO.eq.5) .or. (NLO.eq.-4)   ! q^mu/m_V terms of Pentagon
!      lpent = lpt !.or. lpq
!!      Hexagons
!      lhex= (NLO.eq.-6) .or. (NLO.eq.7).or. (NLO.eq.-7)      
!      lbox2=(NLO.eq.-8)

      ltoploop = (NLO.eq.-5)
!      ltoploop=.true.

      if (ldoblha .and. lnlo) then
        mursq(1,1) = blha_scale**2
        lbox =  .true. 
        ltoploop =  .true. 
      endif

!     sbos < 0 used to signal call from WHjj-part
      bos=abs(sbos)

c      if (linit) then
c         linit = .false.
      if (bos.eq.3 .or. bos.eq.36) then
         k = 1               ! udwhj
      elseif (bos.eq.4 .or. bos.eq.46) then
         k = 2               ! duwhj
      else 
         print *, "Wrong boson ids:", bos
         stop
      endif
      zm2i(3) = 1d0/dcmplx(xm2(3),-xmg(3))
c      endif

      do isig = -1,1,2
         do glupol = 1,2
            do i = 1,3
               mat(k,isig,glupol,i) = (0.0d0, 0.0d0)
               matv(k,isig,glupol,i) = (0.0d0, 0.0d0)
               matb(k,isig,glupol,i) = (0.0d0, 0.0d0)
            enddo
         enddo
      enddo

c define the internal momenta
 
      
      do mu = 0,3
         do i = 1,6 !qq v e h g
            p(mu,i) = pbar(mu,i)*fsign(i)
         enddo
         q(mu) = qbar(mu)*gsign
      enddo

      if (ldebug) then
         print*, "____ enter __FILE__ _____"
        !print*, gsign
         print*, 'q', pbar(0:3,1)
         print*, 'q', pbar(0:3,2)
!         print*, pbar(0:3,n_p+n_v) !6 in qq->WHqq case

         print*, 'w', pbar(0:3,4)+ pbar(0:3,3) !W
         print*, 'h', pbar(0:3,5)
         print*, 'g', pbar(0:3,6)


         !TODO: fix or remove with higgsdecay
         do mu=0,3
         ! momentum conservation test for the simple case 
         ! q qbar -> ... g only !
            momsum(mu) = pbar(mu,1)+pbar(mu,2)
            
            momsum(mu) = momsum(mu) - pbar(mu,3)- pbar(mu,4)- pbar(mu,6)
            if(bos.eq.36 .or. bos.eq.46) then
               momsum(mu) = momsum(mu) - pbar(mu,5)
            endif
         enddo
         if (dotrr(momsum,momsum) > 1e-10) then
         print*, "momsum ", momsum(0) ,momsum(1), momsum(2), momsum(3)
         endif
      endif 
c get the external quark spinors (including factor sqrt(2E) )

      call psi0m(2,pbar(0,1),fsign(1),psi)

c get other bra and ket vectors

c  Get the gluon polarization vector and the gluon emission spinors


      do glupol = 1,2                ! 2 gluon polarizations

         call polvec(qbar,glupol,eps(0,glupol)) ! get gluon pol.vectors

c---- QCD gaugecheck --- eps_mu -> k_mu --- k_mu M^mu = 0 ------------

         do mu = 0,3                   !QCD gaugecheck      
c           eps(mu,glupol) = qbar(mu)
            Myeps(mu,glupol)=eps(mu,glupol)
         enddo

c--- End of QCD gaugecheck -------------------------------------------    

c     NOTES for bras and kets: .true. if psi is a 2-spinor of the chi
c     form as output by psi0m, .false. otherwise.  the last entry is
c     the sum of the two momenta (p plus q) and effectively the
c     momentum of the new spinor.
 
         call ket2r(psi(1,-1,1),.true.,p(0,1),-1,q,eps(0,glupol),
     &        braketg(1,-1,1,glupol),pq(0,1)) ! |q,1>_l,isig
         
         call bra2r(psi(1,-1,2),.true.,p(0,2),-1,q,eps(0,glupol),
     &        braketg(1,-1,2,glupol),pq(0,2)) ! <2,q|_l,isig2
         
c     braketg contains the free quark spinors multiplied by a fermion
c     propagator and a gluon eps_slash. 
c     NOTATION: braketg(2 component spinor, isig =-1 or 1 (fermion hel.),
c     fermion ID = 1:4, gluon polarization glupol=1:2)
         
      enddo    ! 2 gluon polarizations

c Calculation of processes qq -> W g -> W H g ---------------------
      do glupol =1,2
         
         call curr6(-1,psi(1,-1,2),p(0,2),
     &        braketg(1,-1,1,glupol),pq(0,1),jqgwhq(0,-1,glupol))
         
         call curr6(-1,braketg(1,-1,2,glupol),pq(0,2),
     &        psi(1,-1,1),p(0,1),jqwhgq(0,-1,glupol))
         
         do mu = 0,5
            jqqwh(mu,-1,glupol) = jqgwhq(mu,-1,glupol) + jqwhgq(mu,-1,glupol)
         enddo
         
c     contract with wtowh tensor
         ma = -dotcc(wtowh(0,N),jqqwh(0,-1,glupol))
c     wtowh is universal (both signs)

         mat(k,-1,glupol,3) = ma*clr(3,3,-1)
         

      enddo                     ! gluon polarization
  
      if (ldebug .and. .false.) then
         print*, psi(1:2,-1,1)
         print*, psi(1:2,-1,2)
c        print*, pq(0:4,2)
         print*, jqgwhq(0:5,-1,1)
         print*, jqgwhq(0:5,-1,2)
         print*, jqwhgq(0:5,-1,1)
         print*, jqwhgq(0:5,-1,2)
         print*, wtowh(0:3,N)
         print*, mat(k,-1,1:2,3)
      endif

c**************************************************************
c
c    Virtual contributions
c
c**************************************************************


c**************************************************************
c
c    Gluon/Quark loop contributions
c
c**************************************************************
!#ifdef WITH_DIBOSONJET
#ifdef WITH_NLO 
      if (lnlo) then
         if(ldebug) print*, "in lnlo"

C Som variables to perform checks of the virtual
C contributions
C To checl Factorization of Divergences
C  Div=1 -> 1/Eps;   Div=2 ->1/Eps^2
         Div=0          ! different checks for values 1 or 2
         PerformGauge=0 ! different checks for values 1 or 2
         Reno_inv=0     ! checks for value 1 

C The only piece that is not checked with this variables 
C is the finite terms from the cancellation of the poles
C against the Real part. Ask Christoph for this

c Some factors
         cfc=dcmplx(cf)
         cnf=dcmplx(nf)
         cac=dcmplx(ca)
         cacgraph=cac
         cfcac=(cfc-cac/2d0)
         pic2=dcmplx(pi**2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            Finite Terms from Pacos expansion     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         do mu = 0,3
            p1(mu) = p(mu,1)
            p2(mu) = -p(mu,2)
            mpg(mu) = -q(mu)
         enddo

         mans=2d0*dotrr(p1,p2) !assuming massless (p1)**2,(p2)**2
         manu=2d0*dotrr(p2,mpg)!assuming massless (p2)**2,(mpg)**2
         mant=2d0*dotrr(p1,mpg)!assuming massless (p1)**2,(mpg)**2
         musq=mursq(1,1) ! renormalization scale


 100     lns=log(Abs(mans))-(0,1)*pi*theta(mans)-log(musq)
         lnt=log(Abs(mant))-(0,1)*pi*theta(mant)-log(musq)
         lnu=log(Abs(manu))-(0,1)*pi*theta(manu)-log(musq)

! This terms is a guess to make a substraction of each of the contributions to make
! them small, so that, we can decrease statistics in boxline,penlines, hexlines.
! So, it will give the M_v_born: K_cont*M_born
!       K_cont=(-2d0*pi*pi/3d0+50d0/9d0)*CAC-16d0/9d0*TF*cnf
!     1               +2d0*CFC*(-pi*pi+3d0)        
         K_cont=0d0
! This is the finite piece that has to be added to the finite contributions.
! It includes also the counterterms. We have factor out musq^-eps which goes to the 
! running of the alpha_s

         convfact=-2*EulerGamma*Inv3*cnf*TF + 2*Inv3*cnf*TF*Log(4*Pi) + CFC*(-1 + 3*EulerGamma - EulerGamma**2 - Inv2*Inv3*Pi**2 
     &   + Log(musq)*(2*EulerGamma - 2*Log(4*Pi)) - 
     -     Log(4*Pi)**2 - 2*EulerGamma*Log(Abs(mans)) + (0,2)*EulerGamma*Pi*theta(mans) + Log(4*Pi)
     & *(-3 + 2*EulerGamma + 2*Log(Abs(mans)) - (0,2)*Pi*theta(mans))) + 
     -     CAC*((11 - 3*EulerGamma)*EulerGamma*Inv2*Inv3 - Inv2**2*Inv3*Pi**2 + Log(musq)*(EulerGamma - Log(4*Pi)) - 
     &  Inv2*Log(4*Pi)**2 + EulerGamma*Log(Abs(mans)) - 
     -     EulerGamma*Log(Abs(mant)) - EulerGamma*Log(Abs(manu)) + (0,1)*EulerGamma*Pi*(-theta(mans)
     &    + theta(mant) + theta(manu)) + 
     -     Log(4*Pi)*(EulerGamma - 11*Inv2*Inv3 - Log(Abs(mans)) + Log(Abs(mant)) + Log(Abs(manu))
     &  - (0,1)*Pi*(-theta(mans) + theta(mant) + theta(manu))))


! The K_contcheck is used to check the renormalization invariance of the hexline
! hexlinenoAbe, penline, penlinenoAbe, boxline, boxlineNoABe
! Mv_(musq=s) + (convfact+ Kcontcheck= Kcontcheck4-K_cont)*mborn = Mv_(musq=s^prime) + (convfact+ Kcontcheck+ Kcontcheck4-K_cont)*mborn
! Basically, K_contcheck, include additional terms from the musq^-eps that is factor 
! out to define alpha_s(musq). Still, we have to subtract the counterterms to get
! the finite part of the hexlines,penline's independent of musq.
      If(Reno_inv.gt.0) then
         K_contcheck=-4*Inv2*Inv3*cnf*TF*Log(musq) + 
     1              CFC*Log(musq)*(3 - 2*EulerGamma + Log(musq) + 2*Log(4*Pi) - 
     2               2*Log(Abs(mans)) + (0,2)*Pi*theta(mans)) + 
     3              CAC*Log(musq)*(-EulerGamma + 11*Inv2*Inv3 + Inv2*Log(musq) + 
     4         Log((4*Pi*Abs(mans))/(Abs(mant)*Abs(manu))) - (0,1)*Pi*theta(mans) + 
     5            (0,1)*Pi*theta(mant) + (0,1)*Pi*theta(manu))
! This term cancel the musq dependence of the counterterms.
! It should be added
! to prove renormalization invariance of the finite contribution.
         K_contcheck4 =4*Inv2*Inv3*cnf*TF*(EulerGamma + Log(musq) - Log(4*Pi)) + 
     -  Inv2*Inv3*(-11*EulerGamma + Log(4194304d0) - 11*Log(musq) + 11*Log(Pi))*CAC
!    This terms provide the full finite part that have to be added to the result.
!  Nevertheless, it is not the term that we have to add since musq^-eps must be 
!  factorize. It is only to check that convfact is right: 
!   K_contcheck1=convfact + K_contcheck
         K_contcheck1=  -2*EulerGamma*Inv3*cnf*TF - 2*Inv3*cnf*TF*Log(musq) + 2*Inv3*cnf*TF*Log(4*Pi) + 
     -  CFC*(-1 + 3*EulerGamma - EulerGamma**2 - Inv2*Inv3*Pi**2 + Log(musq)**2 - 
     -     Log(4*Pi)**2 - 2*EulerGamma*Log(Abs(mans)) + (0,2)*EulerGamma*Pi*theta(mans) + 
     -     Log(4*Pi)*(-3 + 2*EulerGamma + 2*Log(Abs(mans)) - (0,2)*Pi*theta(mans)) + 
     -     Log(musq)*(3 - 2*Log(Abs(mans)) + (0,2)*Pi*theta(mans))) + 
     -  CAC*((11 - 3*EulerGamma)*EulerGamma*Inv2*Inv3 - Inv2**2*Inv3*Pi**2 + 
     -     Inv2*Log(musq)**2 - Inv2*Log(4*Pi)**2 + EulerGamma*Log(Abs(mans)) - 
     -     EulerGamma*Log(Abs(mant)) - EulerGamma*Log(Abs(manu)) + 
     -     (0,1)*EulerGamma*Pi*(-theta(mans) + theta(mant) + theta(manu)) + 
     -     Log(4*Pi)*(EulerGamma - 11*Inv2*Inv3 - Log(Abs(mans)) + Log(Abs(mant)) + Log(Abs(manu)) - 
     -        (0,1)*Pi*(-theta(mans) + theta(mant) + theta(manu))) + 
     -     Log(musq)*(11*Inv2*Inv3 + Log(Abs(mans)) - Log(Abs(mant)) - Log(Abs(manu)) + 
     -        (0,1)*Pi*(-theta(mans) + theta(mant) + theta(manu))))
         print*, ''
         print*, 'convfact',convfact
         Print*, 'K_contcheck',K_contcheck
         Print*, 'K_contcheck+convfact',K_contcheck+convfact
         Print*, 'K_contcheck1',K_contcheck1
         Print*, 'K_contcheck4',K_contcheck4
         Print*, 'K_contcheck+convfact-K_contcheck4',K_contcheck+convfact-K_contcheck4
!   This part proof renormalization invariance following Eq.2.38 of my paper 
!   
         K_contcheck2=CAC*(1/2d0*(lns*lns-lnt*lnt-lnu*lnu)+
     &  3D0*(-2*lns+lnt+lnu))
        K_contcheck3=-CFC*(lns-3d0)*lns
        Print*, 'K_contcheck2',K_contcheck2
      endif ! Finish the scale renormalization check

! This is to add the finite terms coming from the cancelation of the poles 
!  against the Real part. It comes from Christoph's code of WAj, WZj
! whether external particles are initial/final does not matter, so we use 
! the kinematics of the internal ME ordering
!      
! this has a difference of order 1e-12 compared to the
!  initial/final-dependend version used in WAAj@NLO (qqWAAj.F)

      s= 2d0*dotrr(pbar(0,1),pbar(0,2))
      t=-2d0*dotrr(pbar(0,1),qbar(0)) 
      mmu=-2d0*dotrr(pbar(0,2),qbar(0)) 
      mursq1=4.d0*pi*mursq(1,1)

      factortoadd=
     &  (25*als(1,1)*CA)/(9D0*Pi) + (5*als(1,1)*CF)/Pi - 
     &  (11*als(1,1)*CA*EulerGamma)/(12d0*Pi) - 
     &  (3*als(1,1)*CF*EulerGamma)/(2D0*Pi) + 
     &  (als(1,1)*CA*EulerGamma**2)/(4D0*Pi) + 
     &  (als(1,1)*CF*EulerGamma**2)/(2D0*Pi) - 
     &  (7*als(1,1)*CA*Pi)/24d0 - 
     &  (7*als(1,1)*CF*Pi)/12d0 - (8*als(1,1)*Nf*TR)/(9D0*Pi) + 
     &  (als(1,1)*EulerGamma*Nf*TR)/(3D0*Pi) - 
     &  (3*als(1,1)*CA*Log(mursq1/S))/(4D0*Pi) + 
     &  (3*als(1,1)*CF*Log(mursq1/S))/(2D0*Pi) + 
     &  (als(1,1)*CA*EulerGamma*Log(mursq1/S))/(2D0*Pi) - 
     &  (als(1,1)*CF*EulerGamma*Log(mursq1/S))/Pi - 
     &  (als(1,1)*CA*Log(mursq1/S)**2)/(4D0*Pi) + 
     &  (als(1,1)*CF*Log(mursq1/S)**2)/(2D0*Pi) + 
     &  (5*als(1,1)*CA*Log(-(mursq1/T)))/(6D0*Pi) - 
     &  (als(1,1)*CA*EulerGamma*Log(-(mursq1/T)))/(2D0*Pi) - 
     &  (als(1,1)*Nf*TR*Log(-(mursq1/T)))/(6D0*Pi) + 
     &  (als(1,1)*CA*Log(-(mursq1/T))**2)/(4D0*Pi) + 
     &  (5*als(1,1)*CA*Log(-(mursq1/mmU)))/(6D0*Pi) - 
     &  (als(1,1)*CA*EulerGamma*Log(-(mursq1/mmU)))/(2D0*Pi) - 
     &  (als(1,1)*Nf*TR*Log(-(mursq1/mmU)))/(6D0*Pi) + 
     &  (als(1,1)*CA*Log(-(mursq1/mmU))**2)/(4D0*Pi)

      factortoadd=factortoadd/2d0


!     print*, 'fac, conv'
!     print*, factortoadd
!     print*, convfact
!     
!     convfacttmp = (-2d0*(-540 + 16*NF + 68*Pi**2 
!    &      + 9*Log(2d0*mursq(1,1)/s) 
!    &      + 3*Log(2d0*mursq(1,1)/s)**2 
!    &      + 3*(-30 + NF)*Log(-2d0*mursq(1,1)/mmu) 
!    &      - 27*Log(-2d0*mursq(1,1)/mmu)**2 
!    &      + 3*(-30 + NF)*Log(-2d0*mursq(1,1)/t) 
!    &      - 27*Log(-2d0*mursq(1,1)/t)**2))/9d0

!     print*, convfacttmp
!     print*, convfacttmp*als(1,1)
!     print*, convfacttmp*als(1,1)/(2d0*pi)
!     print*, convfacttmp*als(1,1)/(4*pi)


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C                 Here start the virtual stuff
C
!         LBOX=.FALSE.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Start BOXES
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C               qq -> W g -> W H g -----
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if(lbox) then

!         murbox=musq  !mursq(1,1)!-2d0*dotrr(p(0,1),p(0,2))
         murbox=mursq(1,1)
         if(ldebug) then
            print*, "in lbox"
            print*, 'scale', murbox
            print*, 'scale2', mursq(1,1)
         endif

         do mu = 0,3
            p1(mu) = p(mu,1)
            p2(mu) = -p(mu,2)
            mpg(mu) = -q(mu)
            mqwh(mu) = -p1(mu)-p2(mu)+q(mu)
            ! only considers partons and therefore independent of weak part/decay
c for gauge test               
            emqwh(mu)=mqwh(mu)
            empg(mu)=mpg(mu)
         enddo


         do mu=0,3
            wwh(mu,1)=wtowh(mu,1)
         enddo
!!! Check Gauge invariance
         If(PerformGauge.eq.1) then
            do mu = 0,3
               Myeps(mu,1)=mpg(mu)
               Myeps(mu,2)=mpg(mu)
            enddo
         elseif(PerformGauge.eq.2) then
            do mu = 0,3
               mutemp=mu
               wwh(mutemp,N) = mqwh(mu)
            enddo
         endif    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C   Abelian Part
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C W g 
C
         !TODO glupol loop!

      call boxlineABETotal(p1,mqwh,mpg,p2,psi(1,-1,2),psi(1,-1,1),
     1    wwh(0,1),Myeps(0,1),-1,murbox,1,3,
     2    3,1,boxcontr1,j5g(1,1,1),boxcontr1b,jbg(1,1),Div)
c For non-abelian gauge and kk/M terms
      call boxlineABETotal(p1,mqwh,mpg,p2,psi(1,-1,2),psi(1,-1,1),
     1    emqwh,Myeps(0,1),-1,murbox,1,3,
     2    9,-1,nothing,j5gg(1,1,1),nothing,jbgg(1,1),Div)
      call boxlineABETotal(p1,mqwh,mpg,p2,psi(1,-1,2),psi(1,-1,1),
     1    wwh(0,1),empg,-1,murbox,1,3,
     2    9,-1,nothing,j5gg2(1,1,1),nothing,jbgg2(1,1),Div)
c Gluon Helicty
       call boxlineABETotal(p1,mqwh,mpg,p2,psi(1,-1,2),psi(1,-1,1),
     1    wwh(0,1),Myeps(0,2),-1,murbox,1,3,
     2    9,-1,boxcontr1,j5g(1,1,2),boxcontr1b,jbg(1,2),Div)
      call boxlineABETotal(p1,mqwh,mpg,p2,psi(1,-1,2),psi(1,-1,1),
     1    emqwh,Myeps(0,2),-1,murbox,1,3,
     2    9,-1,nothing,j5gg(1,1,2),nothing,jbgg(1,2),Div)

c Gluon Helicty
C
C g W
C
         call boxlineABETotal(p1,mpg,mqwh,p2,psi(1,-1,2),psi(1,-1,1),
     1    Myeps(0,1),wwh(0,1),-1,murbox,1,2,
     2    3,1,boxcontr2,j5g(1,2,1),boxcontr2b,jbg(2,1),Div)
c For non-abelian gauge and kk/M terms     
         call boxlineABETotal(p1,mpg,mqwh,p2,psi(1,-1,2),psi(1,-1,1),
     1    Myeps(0,1),emqwh,-1,murbox,1,2,
     2    9,-1,nothing,j5gg(1,2,1),nothing,jbgg(2,1),Div)
         call boxlineABETotal(p1,mpg,mqwh,p2,psi(1,-1,2),psi(1,-1,1),
     1    empg,wwh(0,1),-1,murbox,1,2,
     2    9,-1,nothing,j5gg2(1,2,1),nothing,jbgg2(2,1),Div)
c Gluon Helicty
         call boxlineABETotal(p1,mpg,mqwh,p2,psi(1,-1,2),psi(1,-1,1),
     1    Myeps(0,2),wwh(0,1),-1,murbox,1,2,
     2    9,-1,boxcontr2,j5g(1,2,2),boxcontr2b,jbg(2,2),Div)
         call boxlineABETotal(p1,mpg,mqwh,p2,psi(1,-1,2),psi(1,-1,1),
     1    Myeps(0,2),emqwh,-1,murbox,1,2,
     2    9,-1,nothing,j5gg(1,2,2),nothing,jbgg(2,2),Div)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Non-abelian Part
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C W
C
         call boxlineNoAbeTotal(p1,mqwh,p2,mpg,psi(1,-1,2),psi(1,-1,1),
     1    wwh(0,1),Myeps(0,1),-1,murbox,2,1,
     2    boxcontrNoAbe3,jg5n(1),Div)
c For non-abelian gauge and kk/M terms          
         call boxlineNoAbeTotal(p1,mqwh,p2,mpg,psi(1,-1,2),psi(1,-1,1),
     1    emqwh,Myeps(0,1),-1,murbox,9,-1,
     2    nothing,jg5ng(1),Div)
         call boxlineNoAbeTotal(p1,mqwh,p2,mpg,psi(1,-1,2),psi(1,-1,1),
     1    wwh(0,1),empg,-1,murbox,9,-1,
     2    nothing,jg5ng2(1),Div)
c Gluon Helicty
         call boxlineNoAbeTotal(p1,mqwh,p2,mpg,psi(1,-1,2),psi(1,-1,1),
     1    wwh(0,1),Myeps(0,2),-1,murbox,9,-1,
     2    boxcontrNoAbe3,jg5n(2),Div)
         call boxlineNoAbeTotal(p1,mqwh,p2,mpg,psi(1,-1,2),psi(1,-1,1),
     1    emqwh,Myeps(0,2),-1,murbox,9,-1,
     2    nothing,jg5ng(2),Div)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C   Gauge Boxes Abelian Part
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
C W g 
         call box_gauge(p1,mqwh,mpg,p2,murbox,
     2    boxcontr1,boxcontr1b,wardidtest1)
C g W
         call box_gauge(p1,mpg,mqwh,p2,murbox,
     2    boxcontr2,boxcontr2b,wardidtest2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Non-abelian Part
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C W
         wardidtest3=.false.
         if (abs(jg5ng(1)).lt.1d-3) then
            td=abs(-(j5gg(2,2,1)+j5gg(2,1,1))/2d0-jg5ng(1))
         else
            td=abs(-(j5gg(2,2,1)+j5gg(2,1,1))/2d0/jg5ng(1)+1d0)
         endif
         IF(td.le.1.5d-1) then
            if (abs(jg5ng2(1)).lt.1d-3) then
               td=abs(-(j5gg2(2,2,1)+j5gg2(2,1,1))/2d0-jg5ng2(1))
            else
               td=abs(-(j5gg2(2,2,1)+j5gg2(2,1,1))/2d0/jg5ng2(1)+1d0)
            endif
            IF(td.le.1.5d-1) wardidtest3=.true.
         endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!        Print*, 'wardidtest1',wardidtest1
!        Print*, 'wardidtest2',wardidtest2
!        Print*, 'wardidtest3',wardidtest3
         wardidtest1=wardidtest1.and.wardidtest2.and.wardidtest3
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
         If(Div.ne.0) then
! Check factorization of Divergences. Set the color equal.
! So that I can check separately the CF and CA factorization
!        Print*,'CFCFCF',cfc
            cfc=1d0
            cacgraph=1d0
!        Print*,'CFCFCF',cfc
            cfcac=cfc
            wardidtest1=.True.
         endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     

         if (.not.wardidtest1) then
            do glupol =1,2
               do i = 1,2
                  do mu=1,2
                     j5g(mu,i,glupol) = (0.0d0, 0.0d0)

                     j5gg(mu,i,glupol) = (0.0d0, 0.0d0)
                  enddo
               enddo        
               jg5n(glupol) =     (0.0d0, 0.0d0)
               jg5ng(glupol)=     (0.0d0, 0.0d0)
            enddo                   
         else
            do glupol=1,2
               do i = 1,2
                  j5g(1,i,glupol)=(cfc*j5g(1,i,glupol)+cfcac*j5g(2,i,glupol))
                  j5gg(1,i,glupol)=(cfc*j5gg(1,i,glupol)+cfcac*j5gg(2,i,glupol))
               enddo
               jg5n(glupol)=cacgraph*jg5n(glupol)
               jg5ng(glupol)=cacgraph*jg5ng(glupol)
            enddo 
         endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c   q qbar -> W g
c CHECK BORN 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
         do glupol=1,2

            mv1b =(jbg(2,glupol)+  !g W^*
     2       jbg(1,glupol))   !W^* g
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Additional KK terms
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
            mmv1b =((jbgg(2,glupol))*dotrc(mqwh,wwh(0,1))+  !g W^*
     2       (jbgg(1,glupol))*dotrc(mqwh,wwh(0,1)))*zm2i(3) !W^* g
C Change of sign to agree with Madgraph
            mv1b =-mv1b
            mmv1b =-mmv1b

            matb(k,-1,glupol,3)=(mv1b+mmv1b)* clr(3,3,-1)


            ratioComp=abs(mat(k,-1,glupol,3)/matb(k,-1,glupol,3)  -1d0)

            if(ratioComp.gt.1d-4) then        
               Print*, 'born amplitude comparison with box routines failed'
               Print*,"Born/b_box-1d0",ratioComp
               Print*, "Born",mat(k,-1,glupol,3)
               Print*, 'mv1',mv1*clr(3,3,-1)
               Print*, 'mmv1',mmv1
               !stop
            endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                        VIRTUAL
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Add here the Virtual part. Same as born, but replacing 
C  in the above eq. jb -> j5.
C  In addition the non-Abe terms have to be added.
C No way to check the part.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            mv1 =(j5g(1,2,glupol)+  !g W^*
     2       j5g(1,1,glupol))   !W^* g
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Additional KK terms
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
            mmv1 =((j5gg(1,2,glupol))*dotrc(mqwh,wwh(0,1))+  !g W^*
     2       (j5gg(1,1,glupol))*dotrc(mqwh,wwh(0,1)))*zm2i(3)  !W^* g
C Change of sign to agree with Madgraph
            mv1 =-mv1
            mmv1 =-mmv1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCcc
C  No abelian part
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
            mv2=jg5n(glupol)
            mmv2=(jg5ng(glupol)*dotrc(mqwh,wwh(0,1)))*zm2i(3) 
            mv2 =-mv2
            mmv2 =-mmv2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            IF (wardidtest1 ) THEN

               IF (PerformGauge.eq.1) then
! gauge test for the gluon. The CF part of mv1,mv2,mv3,
! is already zero. mv1,mv2,mv3, then, only contains CA part which
! should cancel agains NOABE, mv4,mv5,mv6
                  Print*, ''
                  Print*, 'Div',Div
                  Print*, 'PerformGauge',PerformGauge

                  Print*, ''
                  Print*, 'wardidtest1',wardidtest1 
                  Print*, 'cacgraph',cacgraph
                  Print*, 'cfc',cfc
                  Print*, 'cfca',cfcac
                  Print*, 'cfca-CF',cfcac-cfc
                  Print*, ''
                  Print*, 'mv1b',mv1b
                  Print*, 'mv1',mv1
                  Print*, 'mv2',mv2
                  Print*, 'Gauge',-mv2/mv1-1d0
                  Print*, 'Gaugetest',abs(-mv2/mv1-1d0)
                  Print*, 'mmv1',mmv1
                  Print*, 'mmv2',mmv2

                  !stop
               elseif(PerformGauge.gt.1) then
!   Here only the sum of all the contributions is gauge
!   invariant, but, with coupling cte equal=1. 
!    Otherwise, there are reamining terms that cancel against
!    penline's and so on due to the difference of clr(u(k),1,-1)
!    and clr(d(k),1,-1)

                  Print*, ''
                  Print*, 'wardidtest1',wardidtest1 
                  Print*, 'cacgraph',cacgraph
                  Print*, 'cfc',cfc
                  Print*, 'cfca',cfcac
                  Print*, 'cfca-CF',cfcac-cfc
                  Print*, ''

                  Print*, 'mv1b',mv1b
                  Print*, ''

                  Print*, 'mv1    ',mv1
                  Print*, ''
                  Print*, 'mv2    ',mv2
                  Print*, ''
                  Print*, 'mmv1    ',mmv1
                  Print*, ''
                  Print*, 'mmv2    ',mmv2
                  Print*, ''
               endif
            endif ! Warditest


            If(wardidtest1) then

               If(Div.ne.0) then
!    Check Factorization of CF and CA part separately        
!    Mind that the born part is already checked against Madgraph.
!    So this test is quite good!!
!    At this point mv1,mv2,mv3, only contains the CF part.
!    The color factor of result(2) has been set to CF above.
!    I have to duplicate the virtual part now but with the 
!    CA part: Copy and paste of above expression. and replace
!    in the first argument 1->2       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                  mv1CA =(j5g(2,2,glupol)+  !g W^*
     2       j5g(2,1,glupol))   !W^* g
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Additional KK terms
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
                  mmv1CA =(j5gg(2,2,glupol))*dotrc(mqwh,wwh(0,1))+  !g W^*
     2       (j5gg(2,1,glupol))*dotrc(mqwh,wwh(0,1))  !W^* g
C Change of sign to agree with Madgraph
                  mv1CA =-mv1CA
                  mmv1CA =-mmv1CA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                  If(Div.eq.2) then
! CF part
                     Print*, ""
                     Print*, "CF PART"
                     Print*, ""        
                     ratioFact=mv1/(-2d0*mv1b)
                     Print*, 'mv1',mv1
                     Print*, 'mv1b',mv1b
                     Print*, "CF Div/Eps^2", ratioFact
                     Print*, ""

                     ratioFact=mmv1/(-2d0*mmv1b)
                     Print*, 'mmv1',mmv1
                     Print*, 'mmv1b',mmv1b
!         Print*, "CF Div/Eps^2", ratioFact
                     Print*, "mmv* should be zero by wardidentities"

                     Print*, ""
                     !Print*, 'SUM of born', mmv1b+mmv2b         

! CA part
                     Print*, ""
                     Print*, ""
                     Print*, "CA PART"
                     Print*, ""
                     Print*, ""
                     Print*, 'mv1CA',mv1CA*(-1d0/2d0)
                     Print*, 'mvNoAbe',mv2
                     Print*, 'mv1b',mv1b
                     ratioFact=(mv1CA*(-1d0/2d0)+mv2)/(-1d0*mv1b)
                     Print*, "CA Div/Eps^2", ratioFact
                     Print*, ""

                     Print*, 'mmv1CA',mmv1CA*(-1d0/2d0)
                     Print*, 'mmvNoAbe',mmv2
                     Print*, 'mmv1b',mmv1b
                     Print*, "mmv* should be zero by wardidentities"
!         ratioFact=(mmv1CA*(-1d0/2d0)+mmv2)/(-1d0*mmv1b)
!         Print*, "CA Div/Eps^2", ratioFact
                     Print*, ""


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                  else   !  Div 1/Eps
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
                     lns=log(Abs(mans))-(0,1)*pi*theta(mans)-log(musq)
                     lnt=log(Abs(mant))-(0,1)*pi*theta(mant)-log(musq)
                     lnu=log(Abs(manu))-(0,1)*pi*theta(manu)-log(musq)

! CF PART
                     cteEps1=(-3d0+2d0*lns)
                     Print*, ""
                     Print*, "CF PART"
                     Print*, ""
                     Print*, 'CF cte',cteEps1
                     Print*, ""

                     ratioFact=mv1/(mv1b*cteEps1)
                     Print*, 'mv1',mv1
                     Print*, 'mv1b*cte',mv1b*cteEps1
                     Print*, "CF Div/Eps", ratioFact
                     Print*, ""

                     ratioFact=mmv1/(mmv1b*cteEps1)
                     Print*, 'mmv1',mmv1
                     Print*, 'mmv1b*cte',mmv1b*cteEps1
                     Print*, "mmv* should be zero by wardidentities"
!         Print*, "CF Div/Eps", ratioFact
                     Print*, ""
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CA part
                     cteEps1=(-lns+lnt+lnu)

                     Print*, "CA PART"
                     Print*, ""
                     Print*, 'cteEps1',cteEps1
                     Print*, ""
                     Print*, 'mv1CA',mv1CA*(-1d0/2d0)
                     Print*, 'mvNoAbe',mv2
                     Print*, 'mv1b*cte',mv1b*cteEps1
                     ratioFact=(mv1CA*(-1d0/2d0)+mv2)/(mv1b*cteEps1)
                     Print*, "CA Div/Eps", ratioFact
                     Print*, ""

                     Print*, 'mmv1CA',mmv1CA*(-1d0/2d0)
                     Print*, 'mmvNoAbe',mmv2
                     Print*, 'mmv1b*cte',mmv1b*cteEps1
                     Print*, "mmv* should be zero by wardidentities"
!         ratioFact=(mmv1CA*(-1d0/2d0)+mmv2)/(mmv1b*cteEps1)
!         Print*, "CA Div/Eps", ratioFact
                     Print*, ""
                  endif    !Div

               endif    ! Div, Compute CA part of ABe contributions
            endif   ! Warditest set to true for Div checks


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                   SUM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            matv(k,-1,glupol,3) = (mv1+mv2+mmv1+mmv2)*clr(3,3,-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

            IF (PerformGauge.ge.1) then
! then, matv(k,-1,glupol,1) should be zero
               Print*,''
               Print*,'matv(k,-1,glupol,3) should be zero'
               Print*,' matv(k,-1,glupol,3)', matv(k,-1,glupol,3)
               Print*, 'mv1*clr(3,3,-1)',mv1*clr(3,3,-1)
               Print*,''
!      elseif(PerformGauge.gt.1) then
! then, matv(k,-1,glupol,1) with couplings corrected should be zero
!      Print*,''
!      Print*,'then, matv(k,-1,glupol,3) with couplings corrected should be zero'
!      Print*,' matv(k,-1,glupol,1)', matv(k,-1,glupol,2)
!      cteEps1=clr(u(k),1,-1)/clr(d(k),1,-1) 
!      matv(k,-1,glupol,2) = (mv1+mv2+mmv1+mmv2)*clr(3,3,-1)
!      Print*,'matv(k,-1,glupol,2) with corrected couplings'
!      Print*,'matv(k,-1,glupol,2)', matv(k,-1,glupol,2)
!      Print*,''
            endif

            if (Div.eq.2) then
               Print*,''
               Print*,'The total contributions 1/eps^2'
               Print*, 'shoulb be -4*born since:'
               Print*, 'set CF=1 and CA=0 in ABE contribu'
               Print*,'matv(k,-1,glupol,3)',matv(k,-1,glupol,3)
               Print*, '-4*matvb(k,-1,glupol,3)',-4*matb(k,-1,glupol,3)
               Print*, 'ratio',matv(k,-1,glupol,3)/(-4*matb(k,-1,glupol,3))
            endif  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               Add and subtract cte's 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            matv(k,-1,glupol,3)=matv(k,-1,glupol,3)+matb(k,-1,glupol,3)*convfact

            ! K_cont is 0 here, so ignore this term
            matv(k,-1,glupol,3)=matv(k,-1,glupol,3)-matb(k,-1,glupol,3)*(K_cont)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         enddo  !  glupol helicity

!  Check renormalization invariance using K_contcheck      
         if(Reno_inv.eq.1) then
            check(1)=matv(k,-1,1,3)+(K_contcheck+K_contcheck4)*matb(k,-1,1,3)
            check(2)=matv(k,-1,2,3)+(K_contcheck+K_contcheck4)*matb(k,-1,2,3)
            check1(1)=matv(k,-1,1,3)
            check1(2)=matv(k,-1,2,3)

            Reno_inv=2
            musq=1d0
            goto 100
            ! goto makes code run again with new musq and reno_inv such
            ! that below check1 and check are comparing the
            ! renormalization invariance at musq=1 and musq given in
            ! config
         endif
         If (Reno_inv.gt.0) then
            Print*, ' '
            Print*, 'chek1 hel 1 ',check1(1)
            Print*, 'check2 hel -1',check1(2)
            Print*, 'chek1 hel 1 ',check(1)
            Print*, 'check2 hel -1',check(2)
            Print*, ''
            Print*, 'check1 hel 1 ',matv(k,-1,1,3)
            Print*, 'check2 hel -1',matv(k,-1,2,3)
            Print*, 'check1 hel 1 ',(matv(k,-1,1,3)+(K_contcheck+K_contcheck4)*matb(k,-1,1,3))
            Print*, 'check2 hel -1',(matv(k,-1,2,3)+(K_contcheck+K_contcheck4)*matb(k,-1,2,3))

            Print*, 'ratio1',check(1)/((matv(k,-1,1,3)+(K_contcheck+K_contcheck4)*matb(k,-1,1,3)))-1d0
            Print*, 'ratio2',check(2)/((matv(k,-1,2,3)+(K_contcheck+K_contcheck4)*matb(k,-1,2,3)))-1d0
       stop
         endif        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      endif ! end lbox
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C End Bosonic Virtual Corrections
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      endif ! lnlo 

      if ((lnlo .and. ltoploop) .or. ldebugmtop .or. ldebugtop) then
C Begin Fermionic Corrections
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                      / h                                   
!                g    /                                      
!       q ---*~~~~~~~/~~~~~~~ g                               
!           /      toploop                                   
!          /---- w                                                 
!       q /                                                  
!                                                            
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       musq=mursq(1,1) ! renormalization scale
       murbox=musq  !mursq(1,1)!-2d0*dotrr(p(0,1),p(0,2))
       DivTop=0 ! no poles in this case
       PerformGaugeTop=0 ! if 1: gauge check eps = p for glu, 2: gauge check for W
       vv = 2.d0*sqrt(xm2(3)*(1-xm2(3)/xm2(2)))/(-clr(2,1,-1)) ! 2*MW*SW/e !higgs vev = 246

       ! equal formulation
       ! 2d0*sqrt(xm2(3))/(v(1,3)*dsqrt(8d0))


!!! Check Gauge invariance
         If(PerformGaugeTop.eq.1) then
            do glupol = 1,2
               do mu = 0,3
                  Myeps(mu,glupol)=mpg(mu)
               enddo
            enddo
         elseif(PerformGaugeTop.eq.2) then
            do mu = 0,3
               !TODO: N.ne.1 even possible here?
               wtoonlyw(mu,N) = qw0(mu,N)
            enddo
         endif
 

!using routines from qqAAWj2.F, ggZZg.F
!============================================================         


       mf(1) = tmass ! top 
       ! can check mass dependence:
!       do i=2,maxfloopmasses
!          mf(i) = tmass*2**(i-1)
!       enddo
!      mf(2) = 1d2*tmass
!      mf(3) = 1d4*tmass


       ! switch between madgraph and vbfnlo polarization basis
!       do glupolchoice=1,1 !,2
       glupolchoice=1 
       
       ! reproducing madgraph gluon basis:
!       if (glupolchoice.eq.2) then
!          do glupol=1,2
!             CALL VXXXXX(qbar(0),0d0 ,(2*glupol-3)*gsign ,gsign,myeps(0,glupol))
!          enddo
!!       endif
!
!       CALL IXXXXX(PBAR(0,1),0d0 ,-fsign(1),fsign(1),fakemadg(1,1)) 
!       CALL OXXXXX(PBAR(0,2),0d0 ,-fsign(2),fsign(2),fakemadg(1,2))
!           
!       CALL FVOXXX(fakemadg(1,2   ),wtoonlyw(0,N),GWF ,0d0    ,0d0    ,fakemadg(1,8   ))
!       call jioxxx(fakemadg(1,1),fakemadg(1,8),GG, 0d0, 0d0 , fakemadg(1,9))
!       !CALL IOVXXX(fakemadg(1,1   ),fakemadg(1,9   ),myeps(0,glupol),GG ,fakemadgamp(glupol))        
!
!       if(ldebugtop .and. .false.) then
!          print*, "----------------- wavefunctions, quarks, vbfnlocode"
!          print*, fakemadg(1:6,1)
!          print*, fakemadg(1:6,2)
!          print*, 'glupol'
!          do glupol=1,2
!             print*, myeps(0:5,glupol)
!          enddo
!          print*, '--end glupol--'
!          print*, wtoonlyw(0:5,N)
!          print*, fakemadg(1:6,8)
!
!          print*, fakemadg(1:6,9)
!          print*, 'gluon prop comparison'
!          print*, jqqw1(0:5,-1) !written in glupolchoice=1
!          print*, jqqw2(0:5,-1) !written in glupolchoice=1
!       endif
!
!
!        else ! using vbfnlo polvec gluon basis


      ! only used for toploop contributions
      ! calculated here instead up with other bra/kets
      ! TODO: improve using chreal (.false.->.true.)
      call ket2c(psi(1,-1,1),.false.,p(0,1),-1,qw0(0,N),wtoonlyw(0,N),
     &              psiw0(1,-1,1),fq0(0,1))

      call bra2c(psi(1,-1,2),.false.,p(0,2),-1,qw0(0,N),wtoonlyw(0,N),
     &              psiw0(1,-1,2),fq0(0,2))


c prepare gluon current
cc <2W|1>
       call curr6(-1,psiw0(1,-1,2),fq0(0,2),
     &            psi(1,-1,1),p(0,1),jqqw(0,-1))
cc <2|W1>
       call curr6add(-1,psi(1,-1,2),p(0,2),
     &            psiw0(1,-1,1),fq0(0,1),jqqw(0,-1))

cc g propagator
       call propagate(-1,jqqw(0,-1),pqqw(0))


       ! to switch between gluon polarisation basis
!       endif !glupolchoice


       do mu=0,3
         qin(mu) = -q(mu)
         qhin(mu) = -qh(mu,1)
       enddo

!       if(ldebug.or. .true.) then
       if(ldebugtop) then
          print*, 'toploop contrib'
!          print*, 'fq0', fq0
!          print*, 'p',p
!          print*, 'jqqw',jqqw
!          print*, 'mf', mf)
!          print*, 'q', q
!          print*, 'qh', qh
!          print*, 'myeps', myep1
          print*, 'musq', musq
          print*, 'vv', vv
          if (sqrt(abs(mass2(pqqw))).lt. 30) print*, 'massglt'
!          if (sqrt(abs(mass2(pqqw))).lt. 200) print*, 'massglt200'
!          if (sqrt(abs(mass2(pqqw))).lt. 250) print*, 'massglt250'
!          if (sqrt(abs(mass2(pqqw))).lt. 300) print*, 'massglt300'
          print*, 'massg', sqrt(mass2(pqqw)), sqrt(-mass2(pqqw))
!          print*, 'mursq', mursq(1,1)

       !check if momenta sum to 0
         print*, 'momenta-sum check'
         print*, pqqw(0:3)+qin(0:3)+qhin(0:3)
         print*, 'pqqw, -qh, -q'
         print*, pqqw(0:3)
         print*, qhin(0:3)
         print*, qin(0:3)
       endif


! contract with fermion triangles
! structure from amplitudes/triboson/qqAAWj.F       

cc g-1-2
        !allow multiple masses       
!       do j=1,maxfloopmasses !different quarks (light, b, t) 
        j=1
          do glupol=1,2
             if (glupol.eq.1) then
                flgauge=1
             else
                flgauge=-1
             endif

! put momenta in myeps and then use:             
!      CALL VVSHXX(W(1,6   ),W(1,11  ),W(1,5   ),mad_GH ,AMP(2   ))     

             call gluonverHiggsEvenT(
     &       mf(j),
     &       pqqw(0),     qin(0),            qhin(0),
     &       jqqw(0,-1),  myeps(0,glupol),
     &       murbox,flgauge,rfvert(1,j,glupol),
     &       DivTop)

              ! numerically identically to above diagrams
              ! removed and added factor 2 below
              ! could be used for a numerical check
!             call gluonverHiggsEvenT(
!     &       mf(j),
!     &       qin(0),          pqqw(0),   qhin(0),
!     &       myeps(0,glupol), jqqw(0,-1),   
!     &       murbox,flgauge,rfvert(2,j,glupol),
!     &       DivTop)

          enddo
!       enddo !maxfloopmasses

!       if(ldebugtop) then
!       !if(ldebug .or. .true.) then
!!          print*, 'rfvert', rfvert
!          print*, 'tmass', tmass
!          print*, 'diag1-top  ', rfvert(1,1,1:2)*mf(1)/vv
!          print*, 'diag2-top  ', rfvert(2,1,1:2)*mf(1)/vv
!          print*, 'rel difference', maxval(abs((rfvert(1,1,1:2)-rfvert(2,1,1:2))/rfvert(1,1,1:2)))
!          do i=2,maxfloopmasses
!             print*, 'diag1-heavy',i, rfvert(1,i,1:2)*mf(i)/vv, 
!     &               maxval(abs((rfvert(1,i,1:2)-rfvert(2,i,1:2))/rfvert(1,i,1:2)))
!!             print*, 'diag2-heavy',i, rfvert(2,i,1:2)*mf(i)/vv
!          enddo
!       endif


       !reset matv to zero
       matv(k,-1,1:2,3) = 0d0

       do glupol=1,2
          i=1 ! possible loop directions
             matv(k,-1,glupol,3) = matv(k,-1,glupol,3) + (
             !from ggZZg.F        
     &       rfvert(i,1,glupol) * mf(1)/vv     !toploop
     &       *2  !both diagrams are the same 
!     &       rfvert(i,1,glupol) * mf(1)/vv     !bottomloop can be put here
     &                    ) * clr(3,3,-1) * TR ! W coupling to the quarkline

       enddo

       if (ldebugmtop) then
          print*, 'adding born and toploop result, born, matv, sum'
          print*, mat(k,-1,1:2,3)
          print*, matv(k,-1,1:2,3)*als(1,1)*oneo4pi

          ! if comparing the sum  of born and toploops:
!          mat(k,-1,1:2,3)=mat(k,-1,1:2,3)+matv(k,-1,1:2,3)*als(1,1)*oneo4pi
          ! if only comparing toploop
          mat(k,-1,1:2,3)=matv(k,-1,1:2,3)*als(1,1)*oneo4pi
          ! a change here has to correspond to a change in mg2_udwhg.F for madg


          print*, mat(k,-1,1:2,3)
       endif

       if(PerformGaugeTop.eq.1) then
          print*, 'Gaugetest for Toploops using Gluon, the following expressions should be zero'
          do glupol = 1,2
             print*, rfvert(1,1,glupol)
             print*, rfvert(2,1,glupol)
          enddo
       elseif(PerformGaugeTop.eq.2) then
          print*, 'Gaugetest for Toploops using W, the following expressions should be zero'
          do glupol = 1,2
             print*, rfvert(1,1,glupol)
             print*, rfvert(2,1,glupol)
          enddo
          print*, jqqw(0:5,-1)
       endif

!       enddo !glupolchoice

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      endif !ltoploop
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#endif  
!#endif dibosonjet
! end of with_nlo and dibosonjet

#ifdef WITH_MADGRAPH
c------------------------------------------------------
c comparison with madgraph

      !if (ldebugm.and..not.lnlo) then
      if (ldebugm .or. ldebugmtop) then

         if (ldebug) then 
         print*, "qqWHj"
!         print*, "psi 1 -1 1: ", psi(1:2,-1,1)
!         print*, "psi 1 -1 2: ", psi(1:2,-1,2)


c        print*, pq(0:4,2)
!        print*, jqghq(0:5,-1,1)
!        print*, jqghq(0:5,-1,2)
!        print*, jqhgq(0:5,-1,1)
!        print*, jqhgq(0:5,-1,2)
         print*, wtowh(0:3,1)
!         print*, mat(1,-1,1:2,4)
      endif

         do isig = -1,1,2
            do i = 1,3
                  
                  m1kb(1,isig,i) = 0d0 ! kartesian basis glupol = 1
                  m2kb(1,isig,i) = 0d0 ! kartesian basis glupol = 2
                  
                  do lh = -1,1,2
                     madgr(i,isig,lh,1) = 0d0
                  enddo
                  
            enddo
         enddo
         
         do isig = -1,1,2
            do lh = -1,1,2
               madgr(3,isig,lh,1) = tampr(2,isig,lh,1)
               madgr(3,isig,lh,1) = tampr(2,isig,lh,1)
!                 do i =1,3
!c                    madgr(i,isig,lh,1) = tampr(i,isig,lh,1)
!c                    madgr(i,isig,lh,1) = tampr(i,isig,lh,1)
!                    if (ldebug) print*, "i, isig, lh, tampr", i, isig, lh, tampr(i,isig,lh,1)
!                 enddo
            enddo
         enddo
         
         do isig = -1,1,2
            do i=1,3
               m1kb(1,isig,i) = (madgr(i,isig,-1,1) - 
     &                          madgr(i,isig,1,1))/sqrt(2d0)       
               m2kb(1,isig,i) = -im*(madgr(i,isig,-1,1) + 
     &                          madgr(i,isig,1,1))/sqrt(2d0)

            enddo
         enddo

         i=3
c        do i=1,3

         ! comparing with madgraphs gluon basis

!         if(ldebug) then
!            isig=-1
!            print*, 'madgraphs polarization basis'
!            do glupol=-1,1,2
!               print*, "madgr(",i,",",isig,",",glupol,") = ",madgr(i,isig,glupol,1)
!            enddo
!            print*, 'vbfnlo polarization basis'
!            print*,"m1kb(",1,",",isig,",",i,") = ",m1kb(1,isig,i)
!            print*,"m2kb(",1,",",isig,",",i,") = ",m2kb(1,isig,i)
!            do glupol=1,2
!               print*,"mat(",k,",",isig,",",glupol,",",i,") = ",mat(k,isig,glupol,i)
!            enddo
!         endif


c print results from checks:

         print*, 'start matrixelement comparison'

         do isig = -1,-1,2
            temp = dreal(mat(k,isig,1,i))**2+dimag(mat(k,isig,1,i))**2
            temp1 = dreal(m1kb(1,isig,i))**2+dimag(m1kb(1,isig,i))**2
            if(abs(temp/temp1 -1).ge.1d-12) then !was 1d-6
               print*,"mat(",k,",",isig,",1,",i,") = ",mat(k,isig,1,i)
               print*,"m1kb(",1,",",isig,",",i,")  = ",m1kb(1,isig,i)
               print*,'mat/m1kb = ', mat(k,isig,1,i)/m1kb(1,isig,i)
               print*,"abs(mat)/abs(m1kb) = ",temp/temp1
               print*, "matrix elements qqWHj2 are only within", abs(temp/temp1-1) 
               print*, "mat-m1kb", mat(k,isig,1,i)-m1kb(1,isig,i)
            else if (ldebug) then
               print*, "matrix elements qqWHj2 are within", abs(temp/temp1-1) 
            endif      

            temp = dreal(mat(k,isig,2,i))**2+dimag(mat(k,isig,2,i))**2
            temp1 = dreal(m2kb(1,isig,i))**2+dimag(m2kb(1,isig,i))**2
            if(abs(temp/temp1 -1).ge.1d-12) then
               print*,"mat(",k,",",isig,",2,",i,") = ",mat(k,isig,2,i)
               print*,"m2kb(",1,",",isig,",",i,")  = ",m2kb(1,isig,i)
               print*,'mat/m2kb = ', mat(k,isig,2,i)/m2kb(1,isig,i)
               print*,"abs(mat)/abs(m1kb) = ",temp/temp1
               print*, "matrix elements qqWHj2 are only within", abs(temp/temp1-1) 
               print*, "mat-m1kb", mat(k,isig,2,i)-m2kb(1,isig,i)
            else if (ldebug) then
               print*, "matrix elements qqWHj2 are within", abs(temp/temp1-1) 
            endif      

         enddo
c        enddo
         

      endif ! end of debugging block
c------------------------------------------------------
#endif


c sum the graphs, square them and map them onto uuwwaj and ddwwaj

      res(k) = 0
      resv(k) = 0
      do glupol = 1,2
         mm(k,-1,glupol) = 0
         mmb(k,-1,glupol) = 0
         do i = 1,3
            mm(k,-1,glupol) = mm(k,-1,glupol) +  mat(k,-1,glupol,i)
            mmb(k,-1,glupol) = mmb(k,-1,glupol) +  matb(k,-1,glupol,i)
         enddo
         
         res(k)= res(k) + dreal(mm(k,-1,glupol))**2
     &        + dimag(mm(k,-1,glupol))**2*(1d0)
         
      enddo
c  add Born type term and multiply by F_q = alphas*C_2/4pi
c  the additional factor for the born term is after adding the 
c  subtraction term and the counter term for the renormalization
c  of the pdfs (see also D(x) in m2s_qqWWZj.f)

#ifdef WITH_NLO
      if (lnlo) then
         do glupol = 1,2
            mv(k,-1,glupol) = 0
            do i = 1,3
               mv(k,-1,glupol) = mv(k,-1,glupol) +  matv(k,-1,glupol,i)
            enddo
         enddo
      endif
#endif
      do glupol = 1,2
         if (nlo.eq.4) then
#ifdef WITH_NLO
            mv(k,-1,glupol) = als(1,1)*oneo4pi*
     1         (mv(k,-1,glupol) + mm(k,-1,glupol)*K_cont) 
            !K_cont = 0 here
#endif
c         elseif(nlo.eq.7) then
c#ifdef WITH_NLO
c            mv(k,-1,glupol) = als(1,1)*oneo4pi*
c     1         (mv(k,-1,glupol) + mm(k,-1,glupol)*K_cont) +factortoadd*mmb(k,-1,glupol)
c#endif

            !Toploops (as fake pentagons, nlo.eq.-5)
         else if(nlo.eq.-5) then !TODO: do this or not? mmb=0 anyways
            mv(k,-1,glupol) = als(1,1)*oneo4pi*(mv(k,-1,glupol))
         !   print*, 'mmb'*mmb(k,-1,glupol)
         !   print*, 'factortoadd', factortoadd
         else if(nlo.lt.0 .or. nlo==1) then
#ifdef WITH_NLO
            mv(k,-1,glupol) = als(1,1)*oneo4pi*(mv(k,-1,glupol))+factortoadd*mmb(k,-1,glupol)
#endif
         else
            mv(k,-1,glupol) = (0d0, 0d0)
         endif


         ! take absolute value of matrix element to estimate maximum effect of toploops
         !
!         if (nlo.eq.-5) then
!            resv(k) = resv(k) + 2d0*abs(dreal(mm(k,-1,glupol)*conjg(mv(k,-1,glupol))))
!         else
         resv(k) = resv(k) + 2d0*dreal(mm(k,-1,glupol)*conjg(mv(k,-1,glupol)))
!         endif

      enddo

      if (nlo.ge.0) then
          res(k) = (res(k)+resv(k))*4d0*fpials(2)  ! C_2*3 is the color factor
          tree(k) = res(k)*4d0*fpials(2)
      else
         res(k) = resv(k)*4d0*fpials(2)           ! C_2*3 is the color factor
      endif  
      
      final = res(k)

      if(ldebug) then
         print*, 'mv', mv
         print*, 'res', res
         print*, 'final', final
      endif

      if (ldoblha) then
         select case(blha_amptype(blha_curproc)) 
         case(0) ! tree
            ! blha_amp(1)=blha_amp(1)+
            ! &      res(blha_idsubproc(blha_cursubproc,blha_curproc))
            blha_amp(1)=blha_amp(1) + res(k)

         case(1) ! loop
            lnmusqp1p2 = log(blha_scale**2/(-2d0*dotrr(p(0,1),p(0,2)) )) 

            ! eps**-2
            blha_amp(1)=blha_amp(1)+
     &      blha_CF/(2d0*pi)*2d0*
     &       als(1,1)*(-1d0)*
     &       tree(blha_idsubproc(blha_cursubproc,blha_curproc))
            ! eps**-1
            blha_amp(2)=blha_amp(2)+
     &      1d0/(2d0*pi)*2d0*
     &       als(1,1)*(-blha_gammaQuark-blha_CF*lnmusqp1p2)*
     &       tree(blha_idsubproc(blha_cursubproc,blha_curproc))
c this is more complicated than necessary, 
c but should work for other VBF processes as well
            ! eps**0
            ! idsubproc = 1 here
            blha_amp(3)= blha_amp(3)+
     &      resv(k)*blha_CA ! loop + I (no tree)
     &      + (-(als(1,1))/2d0/pi*                        ! -I
     &          ( ! cvirt not used here, =0 ? ! blha_CF*cvirtVV+  ! our cvirt contribution to I operator
     &            2d0*(blha_KQuark+blha_gammaQuark+ ! K, gamma from CS (C.28)
     &                 blha_tgammaQuark+ ! catani convention to I operator
     &                  (pi**2/6d0-1)*blha_CF) ) ! difference 1/Gamma(1-e) vs. Gamma(1+e)
     &         +als(1,1)/pi*                                       ! +(mu^2/2p1.p2)^epsilon
     &         (-blha_gammaQuark*lnmusqp1p2-blha_CF/2d0*lnmusqp1p2**2)
     &        )*tree(blha_idsubproc(blha_cursubproc,blha_curproc))
            ! tree
            blha_amp(4)= blha_amp(4)+
     &       tree(blha_idsubproc(blha_cursubproc,blha_curproc))

            if(ldebug) then
               print*, 'blha'
               print*, blha_amp(1:4)
               print*, tree, lnmusqp1p2
               print*, matv(k,-1,1:2,1:3)
            endif

         case(2) ! cctree
            ! 1-2
            call BLHA_cctree(1,2,
     &      -blha_CF*
     &       tree(blha_idsubproc(blha_cursubproc,blha_curproc))
     &    )

         case default 
            call BLHA_amptypeerror(blha_amptype(blha_curproc),__FILE__,__LINE__)
         end select
      endif
      
      return
      end


