      subroutine qqwwqq(pbar,sign, nlo, L,
     1                  uucc,uuss,ddcc,ddss,udsc,ducs)
          use globalvars, only: ldoblha
      implicit none
c
c	Dieter Zeppenfeld, <dieter@particle.physik.uni-karlsruhe.de>
c	Initial version:  2005 March 18
c	Last modified:    2005 
C
C  QQWWQQ calculates the matrix elements**2 for electroweak
c  weak boson pair production by quark quark scattering
C
C        q1 q3    ---->   q2 q4 W+W-,   W+ ---> f5-bar f6, W- ---> f7-bar f8
C
C  and crossing related processes. Pauli interference terms for
c  identical fermions are neglected. In particular, only the
c  t-channel exchange of elctroweak bosons is considered. s-channel
c  production of 3 weak bosons is NOT implemented.
c
C  This code is modified to allow for virtual corrections, more precisely
C  the interference of Born with the finite part of virtual diagrams
C  for 
c
c  INPUT:  NLO = 1       return uucc = |M_born|^2 + 2Re(M_Born^* M_virt)
c          NLO = 0       return uucc = |M_born|^2   etc.
c   for steering the calculation of finite box and pentagon contributions
c   the following additional options are implemented
c          NLO = +4      set all finite box and pentagon contributions to 0
c                        but include cvirt*M_Born contribution
c          NLO = +5      true finite pentagons = 0, all other virtual contributions
c                  #included, including leading q^mu/m_V terms of 
c                        pentagons which can be related to boxes
c        remaining cases:   return uucc = 2Re(M_Born^* M_virt), more precisely
c        the piece which is missing for the full NLO=1 when calling +|NLO|
c          NLO = -4      finite Born + q^mu/m_V terms of Pentagon contributions only, no cvirt
c          NLO = -5      true finite Pentagon contributions only
c
c  integer L (=1,2) is used for the subtraction code to select leptonic 
c  tensors for gluon emission off the upper Line (L=1) or lower quark 
c  line (L=2)
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
c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"
c
c variables for the main part of the program
c
      double precision  pbar(0:3,4+max_v), fac, musq
      double precision uucc,uuss,ddcc,ddss,udsc,ducs,res(6),resv(6)
      double precision tree(6)
      double precision  p(0:3,4+max_v), p21(0:4), p43(0:4), pwp(0:4),
     1                  pwm(0:4), pww(0:4)
      integer  sign(4+max_v), nlo, mu, i, j, jj, k, kk, id,
     1         isig, isig1, isig3
      integer  ifl(4,6), js1, js3, L, Lold  !, is1, is3
      double complex prop21(4), prop43(4)
      double complex mat(6,-1:1,-1:1,9), matv(6,-1:1,-1:1,9)
      double complex mm(6,-1:1,-1:1), 
     1               mv12(6,-1:1,-1:1), mv34(6,-1:1,-1:1)
      double complex maa, maz, mza, mzz, mww5, mww6
      double complex  m1, m2, ma(2), mz(2), mpm(4), mmp(4), m5(3,3:4)
      double complex mv1,mv2,mva(2),mvz(2), mv5(3,3:4)
      double complex psi(2,-1:1,4), jqq(0:5,-1:1,2)
      double complex psiwp(2,4), psiwm(2,4), bkjqq(2,-1:1,-1:1,4),
     1               jwp(0:5,-1:1,4), jwm(0:5,-1:1,4),
     1               jvwp(0:5,-1:1,-1:1,4), jvwm(0:5,-1:1,-1:1,4),
     2               qepswm(6), qepswp(6)
      double complex ja(0:5,-1:1,-1:1,4), jz(0:5,-1:1,-1:1,4), 
     1               jva(0:5,-1:1,-1:1,4), jvz(0:5,-1:1,-1:1,4),
     2               psia(2,-1:1,-1:1,4), psiz(2,-1:1,-1:1,4),
     3               mvpm(4,2), mvmp(4,2)
      double complex j5pm(0:3,3,-1:1,2), j5mp(0:3,3,-1:1,2)


      double complex epsa(0:3), epsz(0:3), 
     1               epsNCwp(0:5,-1:1,3:4,2), epsNCwm(0:5,-1:1,3:4,2),
     2               epsCCwp(0:5,-1:1,3:4,2), epsCCwm(0:5,-1:1,3:4,2)
      double precision fqp(0:4,4), fqm(0:4,4), fq(0:4,4),
     1                 qpm(0:4), qmp(0:4), dummy(0:4)
      double complex jj21m,jj21p,jj43m,jj43p, zpm, zmp, zm2i(2:4)
      double complex z1,zp(4:5),zm(4:5), qepsCCwp(2), qepsCCwm(2)
      double complex propt(-1:1,-1:1,6,2), prop_pm(3), prop_mp(3)
      double complex contract_Tjj, dotcc, dotrc, dotqj, s1c
      external contract_Tjj, dotcc, dotrc, dotqj, s1c
      logical ldebug, ldebugv, ldebugvf, linit,lerror,lgc(4)
      data linit /.true./, lerror /.false./, ldebug /.false./
      data lgc /4*.false./
      integer*8 icb1, icount1, icb2, icount2
      data icb1/0/,icount1/0/,icb2/0/,icount2/0/,Lold/1/
      double precision xv1, xv2, xgc1, xgc2
      parameter (ldebugv=.false., ldebugvf=.false.)
      save ifl, Lold, zm2i, linit, lgc, icb1, icb2, icount1, icount2
      save ja,jz,jva,jvz,psia,psiz,fq,j5pm,j5mp
      double complex zero
      parameter (zero = (0d0,0d0) )
      real*8 widthfactor
     
      integer icountmax
!      integer iterations1, iterations2, n2max1, n2max2
!      integer lflavr(5:4+nv)
!      character*50 gridname1, gridname2
!      common /bmainpr/ iterations1, iterations2, n2max1, n2max2,
!     2                 lflavr, gridname1, gridname2
c
c variables for virtual corrections
c
	logical bad_gauge,bad_gauge_sin
	common /gauge / bad_gauge,bad_gauge_sin
c
      double precision c2,c2o4pi     !,pi2o3, cvirtc
c      parameter (pi2o3=pi**2/3d0, cvirtc=pi2o3-7d0)
      parameter (c2=4d0/3d0, c2o4pi=c2/4d0/pi)
      double precision lnmusqp1p2, lnmusqp3p4

      logical lnlo, lbox, lpent, lpt, lpq
c      parameter (lnlo=.true.)
      lnlo = NLO.ne.0    ! include some virtual stuff if T
      lbox = NLO.eq.1 .or. NLO.eq.-4 .or. NLO.eq.5    ! call boxline if T
      lpt  = NLO.eq.1 .or. NLO.le.-5                  ! true pentagon contributions on
      lpq  = NLO.eq.5 .or. NLO.eq.-4                  ! q^mu/m_V terms of Pentagon contributions only
      lpent= lpt .or. lpq                             ! include pentagons if T
c
	bad_gauge = .false. ! set F at beginning of each run

        widthfactor=1.0d0
c        if (ofs.and.with_anom) then
c           widthfactor=0.d0
c        endif

c
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

 701  if (ldebug.or.ldebugv) then
         print*," " 
         print*," enter qqwwqq "
      endif
      do k = 1,6
         do isig1 = -1,1,2
            do isig3 = -1,1,2
               do i = 1,9
                  mat(k,isig1,isig3,i) = 0
                  matv(k,isig1,isig3,i) = 0
               enddo
            enddo
         enddo
      enddo
      do k = 0,5
         do isig1 = -1,1,2
            do isig3 = -1,1,2
               do i = 1,4
                  jvwp(k,isig1,isig3,i) = 0
                  jvwm(k,isig1,isig3,i) = 0
               enddo
            enddo
         enddo
      enddo
c
c identify fermion line sign factors
c
      is1 = sign(1)
      is3 = sign(3)
      js1 = (3+sign(1))/2       ! 1 for sign1=-1,2 for sign1=+1
      js3 = (7+sign(3))/2       ! 3 for sign3=-1,4 for sign3=+1
      if (ldebugv) then
         print*," fsign = ",(sign(i),i=1,4)
         print*," is1, is3 = ",is1,is3," js1, js3 = ",js1,js3
         print*," lfs = ",lfs," lfs1 = ",lfs(js1), " lfs3 = ",lfs(js3) 
         if (sign(1).ne.sign(2) .or. sign(3).ne.sign(4) ) then
            print*," quark sign factors ",(sign(i),i=1,4)
            print*," are not implemented yet in qqwwqq"
            stop
         endif
      endif

c
c define the internal momenta
c
      do mu = 0,3
         do i = 1,4+n_v
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo
         p21(mu) = p(mu,2) - p(mu,1)
         p43(mu) = p(mu,4) - p(mu,3)
         pwp(mu) = p(mu,6) - p(mu,5)
         pwm(mu) = p(mu,8) - p(mu,7)
         pww(mu) = pwp(mu) + pwm(mu)
      enddo
      p21(4) = p21(0)**2 - p21(1)**2 - p21(2)**2 - p21(3)**2
      p43(4) = p43(0)**2 - p43(1)**2 - p43(2)**2 - p43(3)**2
      pwp(4) = pwp(0)**2 - pwp(1)**2 - pwp(2)**2 - pwp(3)**2
      pwm(4) = pwm(0)**2 - pwm(1)**2 - pwm(2)**2 - pwm(3)**2
      pww(4) = pww(0)**2 - pww(1)**2 - pww(2)**2 - pww(3)**2

      if (lwarn) then
         do mu=0,3
            dummy(mu) = p(mu,1)+p(mu,3)+p(mu,5)+p(mu,7)-
     1                 (p(mu,2)+p(mu,4)+p(mu,6)+p(mu,8))
         enddo
         if (abs(dummy(0)).gt.1d-8 .or.
     1       abs(dummy(1)).gt.1d-8 .or.
     2       abs(dummy(2)).gt.1d-8 .or.
     3       abs(dummy(3)).gt.1d-8 ) then
            print*," momentum mismatch in qqwwqq "
            print 102," dummy = ",(dummy(mu),mu=0,3)
            print 103," sign  = ",sign
 101        format( " p(", i1, ") = ", 4(f10.3, 2x) )
 102        format(a,4f10.3)
 103        format(a,6i5)
            do i = 1,8
               write(6,101) i, p(0,i), p(1,i), p(2,i), p(3,i)
            end do
         endif 
      endif
c
c get the vector boson propagator factors
c
      prop21(1) = 1/p21(4)
      prop21(2) = 1/dcmplx(p21(4)-xm2(2),xmg(2)*widthfactor)
      prop21(3) = 1/dcmplx(p21(4)-xm2(3),xmg(3)*widthfactor)
      prop21(4) = prop21(3)

      prop43(1) = 1/p43(4)
      prop43(2) = 1/dcmplx(p43(4)-xm2(2),xmg(2)*widthfactor)
      prop43(3) = 1/dcmplx(p43(4)-xm2(3),xmg(3)*widthfactor)
      prop43(4) = prop43(3)
c
c for box-box and BV graphs we need the propagators for t-channel bosons between quark lines
c as seen from upper line these W momenta are INCOMING. They are OUTGOINg as seen from lower 
c line
      do mu = 0,3
         qpm(mu) = pwp(mu)+p(mu,2) - p(mu,1)    !W+ emitted on upper line
         qmp(mu) = pwm(mu)+p(mu,2) - p(mu,1)    !W- emitted on upper line
      enddo
      qpm(4) = qpm(0)**2-qpm(1)**2-qpm(2)**2-qpm(3)**2
      qmp(4) = qmp(0)**2-qmp(1)**2-qmp(2)**2-qmp(3)**2

      prop_pm(1) = clr(3,3,-1)**2/qpm(4)
      prop_pm(2) = clr(3,3,-1)**2/dcmplx(qpm(4)-xm2(2),xmg(2)*widthfactor)
      prop_pm(3) = clr(3,3,-1)**4/dcmplx(qpm(4)-xm2(3),xmg(3)*widthfactor)
      prop_mp(1) = clr(3,3,-1)**2/qmp(4)
      prop_mp(2) = clr(3,3,-1)**2/dcmplx(qmp(4)-xm2(2),xmg(2)*widthfactor)
      prop_mp(3) = clr(3,3,-1)**4/dcmplx(qmp(4)-xm2(3),xmg(3)*widthfactor)
c
c get the external quark spinors (including factor sqrt(2E) )
c
      call psi0m(4,pbar(0,1),sign(1),psi)
c
c get the f-fbar currents J21^mu=jqq(mu,*,1), J43^mu=jqq(mu,*,2) 
c
      call curr6(1,psi(1,-1,2),p(0,2),psi(1,-1,1),p(0,1),jqq(0,-1,1))
      call curr6(1,psi(1,-1,4),p(0,4),psi(1,-1,3),p(0,3),jqq(0,-1,2))
c
c contract with vvtoww tensors to get Vertex-Vertex scattering diagrams
c
      do isig1 = -1,1,2
         do isig3 = -1,1,2
            maa = contract_Tjj(aaww(0,0,L),jqq(0,isig1,1),jqq(0,isig3,2))
            maz = contract_Tjj(azww(0,0,L),jqq(0,isig1,1),jqq(0,isig3,2))
            mza = contract_Tjj(zaww(0,0,L),jqq(0,isig3,2),jqq(0,isig1,1))
            mzz = contract_Tjj(zzww(0,0,L),jqq(0,isig1,1),jqq(0,isig3,2))
            do k = 1,4
               mat(k,isig1,isig3,1) = 
     1              maa*clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3)
     2            + maz*clr(ifl(1,k),1,isig1)*clr(ifl(3,k),2,isig3)
     3            + mza*clr(ifl(1,k),2,isig1)*clr(ifl(3,k),1,isig3)
     4            + mzz*clr(ifl(1,k),2,isig1)*clr(ifl(3,k),2,isig3)
               matv(k,isig1,isig3,1) = (0d0,0d0)
            enddo
            if (ldebug) then
               do k = 1,4
                  m1 = 0
                  do i = 93,131
                     m1 = m1 + tamp(i,isig1,isig3,k)
                  enddo
                  do i = 135,156
                     m1 = m1 + tamp(i,isig1,isig3,k)
                  enddo
                  do i = 160,181
                     m1 = m1 + tamp(i,isig1,isig3,k)
                  enddo
                  print*," vbf for k=",k,"      ",mat(k,isig1,isig3,1)/m1-1
               enddo
            endif
         enddo
      enddo
      mww6 = contract_Tjj(wwww6(0,0,L),jqq(0,-1,1),jqq(0,-1,2))
      mww5 = contract_Tjj(wwww5(0,0,L),jqq(0,-1,2),jqq(0,-1,1))
      mat(5,-1,-1,1) = mww5*clr(3,3,-1)**2
      mat(6,-1,-1,1) = mww6*clr(3,3,-1)**2
      matv(5,-1,-1,1) = (0d0,0d0)
      matv(6,-1,-1,1) = (0d0,0d0)
      if (ldebug) then
         do k = 5,6
            m1 = 0
            do i = 59,76
               m1 = m1 + tamp(i,-1,-1,k)
            enddo
            do i = 79,82
               m1 = m1 + tamp(i,-1,-1,k)
            enddo
            do i = 87,92
               m1 = m1 + tamp(i,-1,-1,k)
            enddo
            print*," vbf for k=",k,"      ",mat(k,-1,-1,1)/m1-1
         enddo
      endif
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
c  attached next to quark number i is stored in jwp(mu,isig,i). Similarly
c  jwm(mu,isig,i) is the corresponding current for a W- attached next to quark i
c 
c  For the virtual amlitudes the notation, e.g.     jvwp(mu,isig,is,i)
c  is used for the boxline correction to a quark line with one W+ attached next to
c  quark #i and a free Lorentz index mu for the second attached EW boson. is=+-1
c  refers to the sign factor of this quark (vs.antiquark line). They are recalculated 
c  only if this quark line sign has not been calculated yet for this phase space point
c  (i.e. lfs = .true.)  Otherwise they are taken from saved previous calculation
 
      do i = 1,3,2
         call ket2c(psi(1,-1,i),.true.,p(0,i),-1,qp,wp,
     1              psiwp(1,i),fqp(0,i))
         call ket2c(psi(1,-1,i),.true.,p(0,i),-1,qm,wm,
     1              psiwm(1,i),fqm(0,i))
         call bra2c(psi(1,-1,i+1),.true.,p(0,i+1),-1,qp,wp,
     1              psiwp(1,i+1),fqp(0,i+1))
         call bra2c(psi(1,-1,i+1),.true.,p(0,i+1),-1,qm,wm,
     1              psiwm(1,i+1),fqm(0,i+1))
         call curr6(-1,psi(1,-1,i+1),p(0,i+1),
     1                 psiwp(1,i),fqp(0,i), jwp(0,-1,i)   )
         call curr6(-1,psiwp(1,i+1),fqp(0,i+1),
     1                 psi(1,-1,i),p(0,i), jwp(0,-1,i+1) )
         call curr6(-1,psi(1,-1,i+1),p(0,i+1),
     1                 psiwm(1,i),fqm(0,i), jwm(0,-1,i)   )
         call curr6(-1,psiwm(1,i+1),fqm(0,i+1),
     1                 psi(1,-1,i),p(0,i), jwm(0,-1,i+1) )
      enddo

#ifdef WITH_NLO 
      if (lbox .and. lfs(js1)) then
         do i = 1,2
            call boxlinec(-1,psi(1,-1,1),psi(1,-1,2),p(0,1),p(0,2),
     1                    .true., wp,3-i,   jwp(0,-1,i), 
     2                    jvwp(0,-1,is1,i) )
            call boxlinec(-1,psi(1,-1,1),psi(1,-1,2),p(0,1),p(0,2),
     1                    .true., wm,3-i,   jwm(0,-1,i), 
     2                    jvwm(0,-1,is1,i) )
         enddo
      endif
      if (lbox .and. lfs(js3)) then
         do i = 3,4
            call boxlinec(-1,psi(1,-1,3),psi(1,-1,4),p(0,3),p(0,4),
     1                    .true., wp,5-i,   jwp(0,-1,i), 
     2                    jvwp(0,-1,is3,i) )
            call boxlinec(-1,psi(1,-1,3),psi(1,-1,4),p(0,3),p(0,4),
     1                    .true., wm,5-i,   jwm(0,-1,i), 
     2                    jvwm(0,-1,is3,i) )
         enddo
      endif
#endif
c
c now calculate the Vertex-box diagrams; get t-channel W currents first
      do k = 1,2
         kk = 3-k
c k=1 and kk=2 is for "box correction" to upper line
c k=2 and kk=1 is for "box correction" to lower line
         if (k.eq.1) then
            zp(4) = -dcmplx(qmp(0),qmp(3))
            zp(5) = -dcmplx(qmp(1),qmp(2))
            zm(4) = -dcmplx(qpm(0),qpm(3))
            zm(5) = -dcmplx(qpm(1),qpm(2))
         else
            zm(4) = dcmplx(qmp(0),qmp(3))
            zm(5) = dcmplx(qmp(1),qmp(2))
            zp(4) = dcmplx(qpm(0),qpm(3))
            zp(5) = dcmplx(qpm(1),qpm(2))
         endif
         do isig = -1,1,2
            call contract_T2j(NCwpa(0,0,k,L),jqq(0,isig,kk),   epsa) !for W+ current
            call contract_T2j(NCwpz(0,0,k,L),jqq(0,isig,kk),   epsz)
            do mu = 0,3
               epsNCwp(mu,isig,3,k) =
     1              epsa(mu)*clr(3,1,isig)+epsz(mu)*clr(3,2,isig)
               epsNCwp(mu,isig,4,k) =
     1              epsa(mu)*clr(4,1,isig)+epsz(mu)*clr(4,2,isig)
            enddo
            call contract_T2j(NCwma(0,0,k,L),jqq(0,isig,kk),   epsa) !for W- current
            call contract_T2j(NCwmz(0,0,k,L),jqq(0,isig,kk),   epsz)
            do mu = 0,3 
               epsNCwm(mu,isig,3,k) =
     1              epsa(mu)*clr(3,1,isig)+epsz(mu)*clr(3,2,isig)
               epsNCwm(mu,isig,4,k) =
     1              epsa(mu)*clr(4,1,isig)+epsz(mu)*clr(4,2,isig)
            enddo
            do mu = 4,5           ! add momentum info to the currents
               epsNCwp(mu,isig,3,k) = zp(mu)
               epsNCwp(mu,isig,4,k) = zp(mu)
               epsNCwm(mu,isig,3,k) = zm(mu)
               epsNCwm(mu,isig,4,k) = zm(mu)
            enddo
         enddo

c and same for the CC processes (W attached to j43 or j21 current)

         isig = -1
         call contract_T1j(CCwpa(0,0,k,L),jqq(0,isig,kk),   epsa) !for W+ current
         call contract_T1j(CCwpz(0,0,k,L),jqq(0,isig,kk),   epsz)
         do mu = 0,3 
            epsCCwp(mu,isig,3,k) =
     1           epsa(mu)*clr(3,1,isig)+epsz(mu)*clr(3,2,isig)
            epsCCwp(mu,isig,4,k) =
     1           epsa(mu)*clr(4,1,isig)+epsz(mu)*clr(4,2,isig)
         enddo
         if (k.eq.1) then
            qepsCCwp(k) = -dotrc(qmp,epsz)*zm2i(2)
         else
            qepsCCwp(k) = dotrc(qpm,epsz)*zm2i(2)
         endif
         call contract_T1j(CCwma(0,0,k,L),jqq(0,isig,kk),   epsa) !for W- current
         call contract_T1j(CCwmz(0,0,k,L),jqq(0,isig,kk),   epsz)
         do mu = 0,3 
            epsCCwm(mu,isig,3,k) =
     1           epsa(mu)*clr(3,1,isig)+epsz(mu)*clr(3,2,isig)
            epsCCwm(mu,isig,4,k) =
     1           epsa(mu )*clr(4,1,isig)+epsz(mu)*clr(4,2,isig)
         enddo
         do mu = 4,5            ! add momentum info to the currents
            epsCCwp(mu,isig,3,k) = zp(mu)
            epsCCwp(mu,isig,4,k) = zp(mu)
            epsCCwm(mu,isig,3,k) = zm(mu)
            epsCCwm(mu,isig,4,k) = zm(mu)
         enddo 
         if (k.eq.1) then
            qepsCCwm(k) = -dotrc(qpm,epsz)*zm2i(2)
         else
            qepsCCwm(k) = dotrc(qmp,epsz)*zm2i(2)
         endif
      enddo
      
      jj21m = dotcc(wm,jqq(0,-1,1))
      jj21p = dotcc(wp,jqq(0,-1,1))
      jj43m = dotcc(wm,jqq(0,-1,2))
      jj43p = dotcc(wp,jqq(0,-1,2))

c now construct the contribution to the amplitude by current contraction
c virtual contributions are assembled in subroutine boxline
      do isig = -1,1,2            ! 2 bosons attached to 12 line
         do k=3,4
            qepswm(k) = -dotrc(qpm,epsNCwm(0,isig,k,1))*zm2i(3)
            qepswp(k) = -dotrc(qmp,epsNCwp(0,isig,k,1))*zm2i(3)
         enddo
         do k = 1,2                !uucc, uuss
            m1 = dotcc(jwp(0,-1,1),epsNCwm(0,isig,ifl(3,k),1))      !uW+(tW-)u   
     1         - qepswm(ifl(3,k))*jj21p
            m2 = dotcc(jwm(0,-1,2),epsNCwp(0,isig,ifl(3,k),1))      !u(tW+)W-u
     1         + qepswp(ifl(3,k))*jj21m
            mat(k,-1,isig,2) = clr(3,3,-1)**2 * (m1+m2)

#ifdef WITH_NLO 
            if (lbox) then
ctestb
               if (ldebugvf) then
            m1 = dotcc(jwp(0,-1,1),epsNCwm(0,isig,ifl(3,k),1))      !uW+(tW-)u   
            m2 = dotcc(jwm(0,-1,2),epsNCwp(0,isig,ifl(3,k),1))      !u(tW+)W-u
	    
                  call boxline6(id, 
     1              psi(1,-1,1),psi(1,-1,2),p(0,1),p(0,2),-1,
     2              wp,epsNCwm(0,isig,ifl(3,k),1),
     3              -m1,   mv1)
		print*,"mv1/m1=",mv1/m1
		stop
                  call boxline6(id, 
     1              psi(1,-1,1),psi(1,-1,2),p(0,1),p(0,2),-1,
     2              epsNCwp(0,isig,ifl(3,k),1),wm,
     3              -m2,   mv2)
                  m1 = -mv1
                  m2 = -mv2
               endif
cteste
               mv1 = dotcc(jvwp(0,-1,is1,1),epsNCwm(0,isig,ifl(3,k),1)) !uW+(tW-)u   
               mv2 = dotcc(jvwm(0,-1,is1,2),epsNCwp(0,isig,ifl(3,k),1)) !u(tW+)W-u
ctestb
               if (ldebugvf) then
                  print*," uucc, uuss: k = ",k
                  print*," m1 ratio ", mv1/m1
                  print*," m2 ratio ", mv2/m2
               endif
cteste
               matv(k,-1,isig,2) = clr(3,3,-1)**2 * (mv1+mv2)
            endif
#endif
         enddo
         do k = 3,4             !ddcc, ddss
            m1 = dotcc(jwp(0,-1,2),epsNCwm(0,isig,ifl(3,k),1))      !d(tW-)W+d
     1         + qepswm(ifl(3,k))*jj21p
            m2 = dotcc(jwm(0,-1,1),epsNCwp(0,isig,ifl(3,k),1))      !dW-(tW+)d
     1         - qepswp(ifl(3,k))*jj21m
            mat(k,-1,isig,2) = clr(3,3,-1)**2 * (m1+m2)
#ifdef WITH_NLO 
            if (lbox) then
               if (ldebugvf) then
            m1 = dotcc(jwp(0,-1,2),epsNCwm(0,isig,ifl(3,k),1)) !d(tW-)W+d
            m2 = dotcc(jwm(0,-1,1),epsNCwp(0,isig,ifl(3,k),1)) !dW-(tW+)d
                  call boxline6(id, 
     1              psi(1,-1,1),psi(1,-1,2),p(0,1),p(0,2),-1,
     2              epsNCwm(0,isig,ifl(3,k),1),wp,
     3              -m1,   mv1)
                  call boxline6(id, 
     1              psi(1,-1,1),psi(1,-1,2),p(0,1),p(0,2),-1,
     2              wm,epsNCwp(0,isig,ifl(3,k),1),
     3              -m2,   mv2)
                  m1 = -mv1
                  m2 = -mv2
               endif
               mv1 = dotcc(jvwp(0,-1,is1,2),epsNCwm(0,isig,ifl(3,k),1)) !d(tW-)W+d
               mv2 = dotcc(jvwm(0,-1,is1,1),epsNCwp(0,isig,ifl(3,k),1)) !dW-(tW+)d
               if (ldebugvf) then
                  print*," ddcc, ddss: k = ",k
                  print*," m1 ratio ", mv1/m1
                  print*," m2 ratio ", mv2/m2
               endif
               matv(k,-1,isig,2) = clr(3,3,-1)**2 * (mv1+mv2)
            endif
#endif
         enddo
         if (ldebug) then
            k = 1
            m1 = 0
            do i = 10,14
               m1 = m1 + tamp(i,-1,isig,k)
            enddo
            m1 = m1+tamp(89,-1,isig,k)+tamp(90,-1,isig,k)+tamp(132,-1,isig,k)+
     1           tamp(157,-1,isig,k)+tamp(158,-1,isig,k)
            print*," BVu for k=",k,"      ",mat(k,-1,isig,2)/m1
            k = 2
            m1 = 0
            do i = 10,14
               m1 = m1 + tamp(i,-1,isig,k)
            enddo
            m1 = m1+tamp(91,-1,isig,k)+tamp(92,-1,isig,k)+tamp(134,-1,isig,k)+
     1           tamp(158,-1,isig,k)+tamp(159,-1,isig,k)
            print*," BVu for k=",k,"      ",mat(k,-1,isig,2)/m1
            k = 3
            m1 = 0
            do i = 72,76
               m1 = m1 + tamp(i,-1,isig,k)
            enddo
            m1 = m1+tamp(89,-1,isig,k)+tamp(90,-1,isig,k)+tamp(132,-1,isig,k)+
     1           tamp(133,-1,isig,k)+tamp(157,-1,isig,k)
            print*," BVu for k=",k,"      ",mat(k,-1,isig,2)/m1
            k = 4
            m1 = 0
            do i = 46,50
               m1 = m1 + tamp(i,-1,isig,k)
            enddo
            m1 = m1+tamp(89,-1,isig,k)+tamp(90,-1,isig,k)+tamp(132,-1,isig,k)+
     1           tamp(158,-1,isig,k)+tamp(157,-1,isig,k)
            print*," BVu for k=",k,"      ",mat(k,-1,isig,2)/m1

         endif
      enddo
      isig = -1
      k = 5                     !udsc
      m1 = dotcc(jwp(0,-1,2),epsCCwm(0,isig,ifl(1,k),1))            !u(tAZ)W+d
      m2 = dotcc(jwp(0,-1,1),epsCCwm(0,isig,ifl(2,k),1))            !uW+(tAZ)d
      z1 =qepsCCwm(1)*jj21p*(clr(ifl(1,k),2,isig)-clr(ifl(2,k),2,isig))
      mat(k,-1,isig,2) = clr(3,3,-1)**2 * (m1+m2+z1)
#ifdef WITH_NLO 
      if (lbox) then
         if (ldebugvf) then
            call boxline6(id, 
     1        psi(1,-1,1),psi(1,-1,2),p(0,1),p(0,2),-1,
     2        epsCCwm(0,isig,ifl(1,k),1),wp,
     3        -m1,   mv1)
            call boxline6(id, 
     1        psi(1,-1,1),psi(1,-1,2),p(0,1),p(0,2),-1,
     2        wp,epsCCwm(0,isig,ifl(2,k),1),
     3        -m2,   mv2)
            m1 = - mv1
            m2 = - mv2
         endif
         mv1 = dotcc(jvwp(0,-1,is1,2),epsCCwm(0,isig,ifl(1,k),1)) !u(tAZ)W+d
         mv2 = dotcc(jvwp(0,-1,is1,1),epsCCwm(0,isig,ifl(2,k),1)) !uW+(tAZ)d
         if (ldebugvf) then
            print*," udsc: k = ",k
            print*," m1 ratio ", mv1/m1
            print*," m2 ratio ", mv2/m2
         endif
         matv(k,-1,isig,2) = clr(3,3,-1)**2 * (mv1+mv2)
      endif
#endif
      k = 6                     !ducs
      m1 = dotcc(jwm(0,-1,2),epsCCwp(0,isig,ifl(1,k),1))            !d(tAZ)W-u   
      m2 = dotcc(jwm(0,-1,1),epsCCwp(0,isig,ifl(2,k),1))            !dW-(tAZ)u
      z1 =qepsCCwp(1)*jj21m*(clr(ifl(1,k),2,isig)-clr(ifl(2,k),2,isig))
      mat(k,-1,isig,2) = clr(3,3,-1)**2 * (m1+m2+z1)
      if (lbox) then
         mv1 = dotcc(jvwm(0,-1,is1,2),epsCCwp(0,isig,ifl(1,k),1)) !d(tAZ)W-u   
         mv2 = dotcc(jvwm(0,-1,is1,1),epsCCwp(0,isig,ifl(2,k),1)) !dW-(tAZ)u
         matv(k,-1,isig,2) = clr(3,3,-1)**2 * (mv1+mv2)
      endif
c
c repeat the same for 2 bosons attached to 34 line
c
      do isig = -1,1,2
         do k=3,4
            qepswm(k) = dotrc(qmp,epsNCwm(0,isig,k,2))*zm2i(3)
            qepswp(k) = dotrc(qpm,epsNCwp(0,isig,k,2))*zm2i(3)
         enddo
         do k = 1,3,2                  !uucc, ddcc
            m1 = dotcc(jwp(0,-1,3),epsNCwm(0,isig,ifl(1,k),2)) !cW+(tW-)c   
     1         - qepswm(ifl(1,k))*jj43p
            m2 = dotcc(jwm(0,-1,4),epsNCwp(0,isig,ifl(1,k),2)) !c(tW+)W-c
     1         + qepswp(ifl(1,k))*jj43m
            mat(k,isig,-1,3) = clr(3,3,-1)**2 * (m1+m2)
#ifdef WITH_NLO 
            if (lbox) then
               if (ldebugvf) then
            m1 = dotcc(jwp(0,-1,3),epsNCwm(0,isig,ifl(1,k),2)) !cW+(tW-)c   
            m2 = dotcc(jwm(0,-1,4),epsNCwp(0,isig,ifl(1,k),2)) !c(tW+)W-c
                  call boxline6(id, 
     1              psi(1,-1,3),psi(1,-1,4),p(0,3),p(0,4),-1,
     2              wp,epsNCwm(0,isig,ifl(1,k),2),
     3              -m1,   mv1)
                  call boxline6(id, 
     1              psi(1,-1,3),psi(1,-1,4),p(0,3),p(0,4),-1,
     2              epsNCwp(0,isig,ifl(1,k),2),wm,
     3              -m2,   mv2)
                  m1 = - mv1
                  m2 = - mv2
               endif
               mv1 = dotcc(jvwp(0,-1,is3,3),epsNCwm(0,isig,ifl(1,k),2)) !cW+(tW-)c   
               mv2 = dotcc(jvwm(0,-1,is3,4),epsNCwp(0,isig,ifl(1,k),2)) !c(tW+)W-c
               if (ldebugvf) then
                  print*," uucc,ddcc lower line: k = ",k," isig1 = ",isig
                  print*," m1 ratio ", mv1/m1
                  print*," m2 ratio ", mv2/m2
               endif
               matv(k,isig,-1,3) = clr(3,3,-1)**2 * (mv1+mv2)
            endif
#endif
         enddo
         do k = 2,4,2                  !uuss, ddss
            m1 = dotcc(jwp(0,-1,4),epsNCwm(0,isig,ifl(1,k),2))      !s(tW-)W+s
     1         + qepswm(ifl(1,k))*jj43p
            m2 = dotcc(jwm(0,-1,3),epsNCwp(0,isig,ifl(1,k),2))      !sW-(tW+)s
     1         - qepswp(ifl(1,k))*jj43m
            mat(k,isig,-1,3) = clr(3,3,-1)**2 * (m1+m2)
            if (lbox) then
               mv1 = dotcc(jvwp(0,-1,is3,4),epsNCwm(0,isig,ifl(1,k),2)) !s(tW-)W+s
               mv2 = dotcc(jvwm(0,-1,is3,3),epsNCwp(0,isig,ifl(1,k),2)) !sW-(tW+)s
               matv(k,isig,-1,3) = clr(3,3,-1)**2 * (mv1+mv2)
            endif
         enddo
         if (ldebug) then
            k = 1
            m1 = 0
            do i = 36,40
               m1 = m1 + tamp(i,isig,-1,k)
            enddo
            m1 = m1+tamp(91,isig,-1,k)+tamp(92,isig,-1,k)+tamp(133,isig,-1,k)+
     1           tamp(134,isig,-1,k)+tamp(159,isig,-1,k)
            print*," BVl for k=",k,"      ",mat(k,isig,-1,3)/m1
            k = 2
            m1 = 0
            do i = 72,76
               m1 = m1 + tamp(i,isig,-1,k)
            enddo
            m1 = m1+tamp(89,isig,-1,k)+tamp(90,isig,-1,k)+tamp(132,isig,-1,k)+
     1           tamp(133,isig,-1,k)+tamp(157,isig,-1,k)
            print*," BVl for k=",k,"      ",mat(k,isig,-1,3)/m1
            k = 3
            m1 = 0
            do i = 10,14
               m1 = m1 + tamp(i,isig,-1,k)
            enddo
            m1 = m1+tamp(91,isig,-1,k)+tamp(92,isig,-1,k)+tamp(134,isig,-1,k)+
     1           tamp(158,isig,-1,k)+tamp(159,isig,-1,k)
            print*," BVl for k=",k,"      ",mat(k,isig,-1,3)/m1
            k = 4
            m1 = 0
            do i = 72,76
               m1 = m1 + tamp(i,isig,-1,k)
            enddo
            m1 = m1+tamp(91,isig,-1,k)+tamp(92,isig,-1,k)+tamp(134,isig,-1,k)+
     1           tamp(133,isig,-1,k)+tamp(159,isig,-1,k)
            print*," BVl for k=",k,"      ",mat(k,isig,-1,3)/m1

         endif
      enddo
      isig = -1
      k = 5                     !udsc
      m1 = dotcc(jwm(0,-1,4),epsCCwp(0,isig,ifl(3,k),2))            !s(tAZ)W-c
      m2 = dotcc(jwm(0,-1,3),epsCCwp(0,isig,ifl(4,k),2))            !sW-(tAZ)c
      z1 =qepsCCwp(2)*jj43m*(clr(ifl(3,k),2,isig)-clr(ifl(4,k),2,isig))
      mat(k,isig,-1,3) = clr(3,3,-1)**2 * (m1+m2+z1)
      if (lbox) then
         mv1 = dotcc(jvwm(0,-1,is3,4),epsCCwp(0,isig,ifl(3,k),2)) !s(tAZ)W-c
         mv2 = dotcc(jvwm(0,-1,is3,3),epsCCwp(0,isig,ifl(4,k),2)) !sW-(tAZ)c
         matv(k,isig,-1,3) = clr(3,3,-1)**2 * (mv1+mv2)
      endif
      k = 6                     !ducs
      m1 = dotcc(jwp(0,-1,4),epsCCwm(0,isig,ifl(3,k),2))            !c(tAZ)W+s   
      m2 = dotcc(jwp(0,-1,3),epsCCwm(0,isig,ifl(4,k),2))            !cW+(tAZ)s
      z1 =qepsCCwm(2)*jj43p*(clr(ifl(3,k),2,isig)-clr(ifl(4,k),2,isig))
      mat(k,isig,-1,3) = clr(3,3,-1)**2 * (m1+m2+z1)
      if (lbox) then
         mv1 = dotcc(jvwp(0,-1,is3,4),epsCCwm(0,isig,ifl(3,k),2)) !c(tAZ)W+s   
         mv2 = dotcc(jvwp(0,-1,is3,3),epsCCwm(0,isig,ifl(4,k),2)) !cW+(tAZ)s
         matv(k,isig,-1,3) = clr(3,3,-1)**2 * (mv1+mv2)
      endif
      if (ldebug) then
         do k = 5,6
            m1 = 0
            do i = 8,12
               m1 = m1 + tamp(i,-1,isig,k)
            enddo
            do i = 35,39
               m1 = m1 + tamp(i,-1,isig,k)
            enddo
            if (k.eq.5) then
               print*," BVu for k=",k,"      ",mat(k,-1,isig,2)/m1
            else
               print*," BVl for k=",k,"      ",mat(k,-1,isig,3)/m1
            endif
            z1 = m1
            m1 = 0
            do i = 55,58
               m1 = m1 + tamp(i,-1,isig,k)
            enddo
            m1 = m1 + tamp(77,-1,isig,k) + tamp(78,-1,isig,k)
            do i = 83,86
               m1 = m1 + tamp(i,-1,isig,k)
            enddo
            if (k.eq.5) then
               print*," BVl for k=",k,"      ",mat(k,-1,isig,3)/m1
            else
               print*," BVu for k=",k,"      ",mat(k,-1,isig,2)/m1
            endif
            print*," BV for k=",k,"       ",(mat(k,-1,isig,3)+mat(k,-1,isig,2))/(m1+z1)
         enddo
      endif

c next come the A/Z-->WW currents attached to the quark lines. For the virtual 
c corrections the most effective structure is the contraction of two polarization 
c vectors with one fermion line. First build these effective polarization vectors
c from the currents aww(mu) and azz(mu)
c
c NOTE: the aww and azz currents are conserved (checked numerically). Hence there is no 
c       need to consider q^mu * q^nu/m_Z^2  terms in the Z boson propagator
      if (lfs(js1) .or. L.ne.Lold) then
         i = 1
         do isig = -1,1,2
            call ket2c(psi(1,isig,i),.true.,p(0,i),isig,pww,aww,
     1                 psia(1,isig,is1,i),fq(0,i))
            call ket2c(psi(1,isig,i),.true.,p(0,i),isig,pww,zww,
     1                 psiz(1,isig,is1,i),fq(0,i))
            call bra2c(psi(1,isig,i+1),.true.,p(0,i+1),isig,pww,aww,
     1                 psia(1,isig,is1,i+1),fq(0,i+1))
            call bra2c(psi(1,isig,i+1),.true.,p(0,i+1),isig,pww,zww,
     1                 psiz(1,isig,is1,i+1),fq(0,i+1))
         enddo
         call curr6(1,psi(1,-1,i+1),p(0,i+1),
     1              psia(1,-1,is1,i),fq(0,i), ja(0,-1,is1,i)   )
         call curr6(1,psia(1,-1,is1,i+1),fq(0,i+1),
     1              psi(1,-1,i),p(0,i), ja(0,-1,is1,i+1) )
         call curr6(1,psi(1,-1,i+1),p(0,i+1),
     1              psiz(1,-1,is1,i),fq(0,i), jz(0,-1,is1,i)   )
         call curr6(1,psiz(1,-1,is1,i+1),fq(0,i+1),
     1              psi(1,-1,i),p(0,i), jz(0,-1,is1,i+1) )
#ifdef WITH_NLO 
         if (lbox) then
            do i = 1,2
               call boxlinec(1,psi(1,-1,1),psi(1,-1,2),p(0,1),p(0,2),
     1                       .true., aww,3-i,   ja(0,-1,is1,i), 
     2                       jva(0,-1,is1,i) )
               call boxlinec(1,psi(1,-1,1),psi(1,-1,2),p(0,1),p(0,2),
     1                       .false., zww,3-i,  jz(0,-1,is1,i), 
     2                       jvz(0,-1,is1,i) )
            enddo
         endif
#endif
      endif
      if (lfs(js3) .or. L.ne.Lold) then
         i = 3
         do isig = -1,1,2
            call ket2c(psi(1,isig,i),.true.,p(0,i),isig,pww,aww,
     1                 psia(1,isig,is3,i),fq(0,i))
            call ket2c(psi(1,isig,i),.true.,p(0,i),isig,pww,zww,
     1                 psiz(1,isig,is3,i),fq(0,i))
            call bra2c(psi(1,isig,i+1),.true.,p(0,i+1),isig,pww,aww,
     1                 psia(1,isig,is3,i+1),fq(0,i+1))
            call bra2c(psi(1,isig,i+1),.true.,p(0,i+1),isig,pww,zww,
     1                 psiz(1,isig,is3,i+1),fq(0,i+1))
         enddo
         call curr6(1,psi(1,-1,i+1),p(0,i+1),
     1              psia(1,-1,is3,i),fq(0,i), ja(0,-1,is3,i)   )
         call curr6(1,psia(1,-1,is3,i+1),fq(0,i+1),
     1              psi(1,-1,i),p(0,i), ja(0,-1,is3,i+1) )
         call curr6(1,psi(1,-1,i+1),p(0,i+1),
     1              psiz(1,-1,is3,i),fq(0,i), jz(0,-1,is3,i)   )
         call curr6(1,psiz(1,-1,is3,i+1),fq(0,i+1),
     1              psi(1,-1,i),p(0,i), jz(0,-1,is3,i+1) )
#ifdef WITH_NLO 
         if (lbox) then
            do i = 3,4
               call boxlinec(1,psi(1,-1,3),psi(1,-1,4),p(0,3),p(0,4),
     1                       .true., aww,5-i,   ja(0,-1,is3,i), 
     2                       jva(0,-1,is3,i) )
               call boxlinec(1,psi(1,-1,3),psi(1,-1,4),p(0,3),p(0,4),
     1                       .false., zww,5-i,   
     2                       jz(0,-1,is3,i), jvz(0,-1,is3,i) )
            enddo
         endif
#endif
      endif

      do isig1 = -1,1,2
         do isig3 = -1,1,2
c box correction to upper line: polarization vectors are 
c    jqq(mu,isig3,2)=j43 with momentum    p43        and
c    aww/zww(mu)         with momentum    pww

            ma(1) =      dotcc(ja(0,isig1,is1,2),jqq(0,isig3,2))
            ma(2) =      dotcc(ja(0,isig1,is1,1),jqq(0,isig3,2))
            mz(1) =      dotcc(jz(0,isig1,is1,2),jqq(0,isig3,2))
            mz(2) =      dotcc(jz(0,isig1,is1,1),jqq(0,isig3,2))
            if (lbox) then
               mva(1) = dotcc(jva(0,isig1,is1,2),jqq(0,isig3,2))
               mva(2) = dotcc(jva(0,isig1,is1,1),jqq(0,isig3,2))
               mvz(1) = dotcc(jvz(0,isig1,is1,2),jqq(0,isig3,2))
               mvz(2) = dotcc(jvz(0,isig1,is1,1),jqq(0,isig3,2))
            endif
            do k = 1,4
               propt(isig1,isig3,k,2) = 
     1            clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3)*prop43(1) 
     2          + clr(ifl(1,k),2,isig1)*clr(ifl(3,k),2,isig3)*prop43(2)
               mat(k,isig1,isig3,4) = propt(isig1,isig3,k,2) *
     1            ( (ma(1)+ma(2))*clr(ifl(1,k),1,isig1) +
     2              (mz(1)+mz(2))*clr(ifl(1,k),2,isig1)   )
               if (lbox) matv(k,isig1,isig3,4) = propt(isig1,isig3,k,2)*
     1            ( (mva(1)+mva(2))*clr(ifl(1,k),1,isig1) +
     2              (mvz(1)+mvz(2))*clr(ifl(1,k),2,isig1)   )
            enddo
c            if (lbox .and. ldebugv) then
c               id = -1
c               call boxline6(id, 
c     1              psi(1,isig1,1),psi(1,isig1,2),p(0,1),p(0,2),isig1,
c     2              jqq(0,isig3,2), aww, 
c     3              -ma(1),   mva(1))
c               m1 = -mva(1)
c               mva(1) = dotcc(jva(0,isig1,is1,2),jqq(0,isig3,2))
c               id = -1
c               call boxline6(id, 
c     1              psi(1,isig1,1),psi(1,isig1,2),p(0,1),p(0,2),isig1,
c     2              jqq(0,isig3,2), zww, 
c     3              -mz(1),   mvz(1))
c               m1 = -mvz(1)
c               mvz(1) = dotcc(jvz(0,isig1,is1,2),jqq(0,isig3,2))
c               id = -2
c               call boxline6(id, 
c     1              psi(1,isig1,1),psi(1,isig1,2),p(0,1),p(0,2),isig1,
c     2              aww, jqq(0,isig3,2),   
c     3              -ma(2),   mva(2))
c               m1 = -mva(2)
c               mva(2) = dotcc(jva(0,isig1,is1,1),jqq(0,isig3,2))
c               id = -2
c               call boxline6(id, 
c     1              psi(1,isig1,1),psi(1,isig1,2),p(0,1),p(0,2),isig1,
c     2              zww, jqq(0,isig3,2),  
c     3              -mz(2),   mvz(2))
c               mvz(2) = -mvz(2)
c               do k = 1,4
c                  m1 = matv(k,isig1,isig3,4)
c                  matv(k,isig1,isig3,4) = propt(isig1,isig3,k,2) *
c     1                ( (mva(1)+mva(2))*clr(ifl(1,k),1,isig1) +
c     2                  (mvz(1)+mvz(2))*clr(ifl(1,k),2,isig1)   )
c                  print*," aww/zww virtual upper line for k =",k,":",
c     1                   m1/matv(k,isig1,isig3,4)
c               enddo
c            endif
            if (ldebug) then
               k = 1
               m2 = 0
               do i = 4,9
                  m2 = m2 + tamp(i,isig1,isig3,k)
               enddo
               do i = 17,26
                  m2 = m2 + tamp(i,isig1,isig3,k)
               enddo
               do i = 55,70
                  m2 = m2 + tamp(i,isig1,isig3,k)
               enddo
               print*," BV4u for k=",k,"     ",mat(k,isig1,isig3,4)/m2
            endif
            if (isig1.eq.-1 .and. isig3.eq.-1) then
               do k = 5,6
                  propt(-1,-1,k,2) = clr(ifl(1,k),3,-1)**2*prop43(3)
                  mat(k,-1,-1,4) = propt(-1,-1,k,2) *
     1       ( ma(1)*clr(ifl(2,k),1,-1) + ma(2)*clr(ifl(1,k),1,-1) 
     2       + mz(1)*clr(ifl(2,k),2,-1) + mz(2)*clr(ifl(1,k),2,-1) )
                  if (lbox) matv(k,-1,-1,4) = propt(-1,-1,k,2) *
     1       ( mva(1)*clr(ifl(2,k),1,-1) + mva(2)*clr(ifl(1,k),1,-1) 
     2       + mvz(1)*clr(ifl(2,k),2,-1) + mvz(2)*clr(ifl(1,k),2,-1) )
               enddo
            else
               do k = 5,6
                  propt(isig1,isig3,k,2) = 0d0
                  mat(k,isig1,isig3,4) = 0d0
                  matv(k,isig1,isig3,4) = 0d0
               enddo
            endif
c box correction to lower line: polarization vectors are 
c    jqq(mu,isig1,1)=j21 with momentum    p21        and
c    aww/zww(mu)         with momentum    pww

            ma(1) = dotcc(ja(0,isig3,is3,4),jqq(0,isig1,1))
            ma(2) = dotcc(ja(0,isig3,is3,3),jqq(0,isig1,1))
            mz(1) = dotcc(jz(0,isig3,is3,4),jqq(0,isig1,1))
            mz(2) = dotcc(jz(0,isig3,is3,3),jqq(0,isig1,1))
            if (lbox) then
               mva(1) = dotcc(jva(0,isig3,is3,4),jqq(0,isig1,1))
               mva(2) = dotcc(jva(0,isig3,is3,3),jqq(0,isig1,1))
               mvz(1) = dotcc(jvz(0,isig3,is3,4),jqq(0,isig1,1))
               mvz(2) = dotcc(jvz(0,isig3,is3,3),jqq(0,isig1,1))
            endif

            do k = 1,4
               propt(isig1,isig3,k,1) = 
     1            clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3)*prop21(1)
     2          + clr(ifl(1,k),2,isig1)*clr(ifl(3,k),2,isig3)*prop21(2)
               mat(k,isig1,isig3,5) = propt(isig1,isig3,k,1) *
     1            ( (ma(1)+ma(2))*clr(ifl(3,k),1,isig3) +
     2              (mz(1)+mz(2))*clr(ifl(3,k),2,isig3)   )
               if (lbox) matv(k,isig1,isig3,5) = propt(isig1,isig3,k,1)*
     1            ( (mva(1)+mva(2))*clr(ifl(3,k),1,isig3) +
     2              (mvz(1)+mvz(2))*clr(ifl(3,k),2,isig3)   )
            enddo
c            if (lbox) then
c               id = -1
c               call boxline6(id, 
c     1              psi(1,isig3,3),psi(1,isig3,4),p(0,3),p(0,4),isig3,
c     2              jqq(0,isig1,1), aww, 
c     3              -ma(1),   mva(1))
c               id = -1
c               call boxline6(id, 
c     1              psi(1,isig3,3),psi(1,isig3,4),p(0,3),p(0,4),isig3,
c     2              jqq(0,isig1,1), zww, 
c     3              -mz(1),   mvz(1))
c               id = -2
c               call boxline6(id, 
c     1              psi(1,isig3,3),psi(1,isig3,4),p(0,3),p(0,4),isig3,
c     2              aww,  jqq(0,isig1,1),
c     3              -ma(2),   mva(2))
c               id = -2
c               call boxline6(id, 
c     1              psi(1,isig3,3),psi(1,isig3,4),p(0,3),p(0,4),isig3,
c     2              zww, jqq(0,isig1,1), 
c     3              -mz(2),   mvz(2))
c               do k = 1,4
c                  m1 = -propt(isig1,isig3,k,1) *
c     1                ( (mva(1)+mva(2))*clr(ifl(3,k),1,isig3) +
c     2                  (mvz(1)+mvz(2))*clr(ifl(3,k),2,isig3)   )
c                  print*," aww/zww virtual lower line for k =",k,":",
c     1                   m1/matv(k,isig1,isig3,5)
c               enddo
c            endif
            if (ldebug) then
               k = 1
               m1 = 0
               do i = 30,35
                  m1 = m1 + tamp(i,isig1,isig3,k)
               enddo
               do i = 43,52
                  m1 = m1 + tamp(i,isig1,isig3,k)
               enddo
               do i = 73,88
                  m1 = m1 + tamp(i,isig1,isig3,k)
               enddo
               print*," BV4l for k=",k,"     ",mat(k,isig1,isig3,5)/m1
               print*," BV4 for k =",k,"     ",(mat(k,isig1,isig3,4)+mat(k,isig1,isig3,5))/(m1+m2)
            endif
            if (isig1.eq.-1 .and. isig3.eq.-1) then
               do k = 5,6
                  propt(-1,-1,k,1) = clr(ifl(1,k),3,-1)**2*prop21(3)
                  mat(k,-1,-1,5) = propt(-1,-1,k,1) *
     1       ( ma(1)*clr(ifl(4,k),1,-1) + ma(2)*clr(ifl(3,k),1,-1) 
     2       + mz(1)*clr(ifl(4,k),2,-1) + mz(2)*clr(ifl(3,k),2,-1) )
                  if (lbox) matv(k,-1,-1,5) = propt(-1,-1,k,1) *
     1       ( mva(1)*clr(ifl(4,k),1,-1) + mva(2)*clr(ifl(3,k),1,-1) 
     2       + mvz(1)*clr(ifl(4,k),2,-1) + mvz(2)*clr(ifl(3,k),2,-1) )
               enddo
            else
               do k = 5,6
                  propt(isig1,isig3,k,1) = 0d0
                  mat(k,isig1,isig3,5) = 0d0
                  matv(k,isig1,isig3,5) = 0d0
               enddo
            endif
         enddo
      enddo
c
c  next do the box-box graphs with one W emitted from the upper and the 
c  other from the lower line. These are only possible for lefthanded quarks
c  on both lines, i.e. isig1 = -1 and isig3 = -1

c  upper line box: eps1 = wp, eps2 = jwm(mu,-1,3)
c  lower line box: eps1 = wm, eps2 = jwp(mu,-1,1)
      mpm(1) = -dotcc(jwp(0,-1,1),jwm(0,-1,3))  !- zpm
      mmp(1) = -dotcc(jwm(0,-1,1),jwp(0,-1,3))  ! -zmp
c  upper line box: eps1 = wp, eps2 = jwm(mu,-1,4)
c  lower line box: eps1 = jwp(mu,-1,1), eps2 = wm
      mpm(2) = -dotcc(jwp(0,-1,1),jwm(0,-1,4)) ! + zpm
      mmp(2) = -dotcc(jwm(0,-1,1),jwp(0,-1,4)) ! + zmp
c  upper line box: eps1 = jwm(mu,-1,3), eps2 = wp
c  lower line box: eps1 = wm, eps2 = jwp(mu,-1,2)
      mpm(3) = -dotcc(jwp(0,-1,2),jwm(0,-1,3)) ! + zpm 
      mmp(3) = -dotcc(jwm(0,-1,2),jwp(0,-1,3)) ! + zmp
c  upper line box: eps1 = jwm(mu,-1,4), eps2 = wp
c  lower line box: eps1 = jwp(mu,-1,2), eps2 = wm
      mpm(4) = -dotcc(jwp(0,-1,2),jwm(0,-1,4)) ! - zpm
      mmp(4) = -dotcc(jwm(0,-1,2),jwp(0,-1,4)) ! - zmp

#ifdef WITH_NLO 
      if (lbox) then
c 1 = upper line box: eps1 = wp, eps2 = jwm(mu,-1,3)
c 2 = lower line box: eps1 = wm, eps2 = jwp(mu,-1,1)
         mvpm(1,1) = -dotcc(jvwp(0,-1,is1,1),jwm(0,-1,3))
         mvmp(1,1) = -dotcc(jvwm(0,-1,is1,1),jwp(0,-1,3))
         mvpm(1,2) = -dotcc(jwp(0,-1,1),jvwm(0,-1,is3,3))
         mvmp(1,2) = -dotcc(jwm(0,-1,1),jvwp(0,-1,is3,3))
         if (ldebugvf) then
            id = -2
            call boxline6(id, 
     1           psi(1,-1,1),psi(1,-1,2),p(0,1),p(0,2),-1,
     2           wp,jwm(0,-1,3),
     3           mpm(1),   mv1)
            id = -2
            call boxline6(id, 
     1           psi(1,-1,1),psi(1,-1,2),p(0,1),p(0,2),-1,
     2           wm,jwp(0,-1,3),
     3           mmp(1),   mv2)
            print*," BB virtual 1,1: ",mv1/mvpm(1,1),mv2/mvmp(1,1)
            id = -2
            call boxline6(id, 
     1           psi(1,-1,3),psi(1,-1,4),p(0,3),p(0,4),-1,
     2           wm,jwp(0,-1,1),
     3           mpm(1),   mv1)
            id = -2
            call boxline6(id, 
     1           psi(1,-1,3),psi(1,-1,4),p(0,3),p(0,4),-1,
     2           wp,jwm(0,-1,1),
     3           mmp(1),   mv2)
            print*," BB virtual 1,2: ",mv1/mvpm(1,2),mv2/mvmp(1,2)
         endif

c  upper line box: eps1 = wp, eps2 = jwm(mu,-1,4)
c  lower line box: eps1 = jwp(mu,-1,1), eps2 = wm
         mvpm(2,1) = -dotcc(jvwp(0,-1,is1,1),jwm(0,-1,4))
         mvmp(2,1) = -dotcc(jvwm(0,-1,is1,1),jwp(0,-1,4))
         mvpm(2,2) = -dotcc(jwp(0,-1,1),jvwm(0,-1,is3,4))
         mvmp(2,2) = -dotcc(jwm(0,-1,1),jvwp(0,-1,is3,4))
         if (ldebugvf) then
            id = -2
            call boxline6(id, 
     1           psi(1,-1,1),psi(1,-1,2),p(0,1),p(0,2),-1,
     2           wp,jwm(0,-1,4),
     3           mpm(2),   mv1)
            id = -2
            call boxline6(id, 
     1           psi(1,-1,1),psi(1,-1,2),p(0,1),p(0,2),-1,
     2           wm,jwp(0,-1,4),
     3           mmp(2),   mv2)
            print*," BB virtual 2,1: ",mv1/mvpm(2,1),mv2/mvmp(2,1)
            id = -1
            call boxline6(id, 
     1           psi(1,-1,3),psi(1,-1,4),p(0,3),p(0,4),-1,
     2           jwp(0,-1,1), wm,
     3           mpm(2),   mv1)
            id = -1
            call boxline6(id, 
     1           psi(1,-1,3),psi(1,-1,4),p(0,3),p(0,4),-1,
     2           jwm(0,-1,1), wp,
     3           mmp(2),   mv2)
            print*," BB virtual 2,2: ",mv1/mvpm(2,2),mv2/mvmp(2,2)
         endif
c  upper line box: eps1 = jwm(mu,-1,3), eps2 = wp
c  lower line box: eps1 = wm, eps2 = jwp(mu,-1,2)
         mvpm(3,1) = -dotcc(jvwp(0,-1,is1,2),jwm(0,-1,3)) 
         mvmp(3,1) = -dotcc(jvwm(0,-1,is1,2),jwp(0,-1,3))
         mvpm(3,2) = -dotcc(jwp(0,-1,2),jvwm(0,-1,is3,3)) 
         mvmp(3,2) = -dotcc(jwm(0,-1,2),jvwp(0,-1,is3,3))
         if (ldebugvf) then
            id = -1
            call boxline6(id, 
     1           psi(1,-1,1),psi(1,-1,2),p(0,1),p(0,2),-1,
     2           jwm(0,-1,3), wp,
     3           mpm(3),   mv1)
            id = -1
            call boxline6(id, 
     1           psi(1,-1,1),psi(1,-1,2),p(0,1),p(0,2),-1,
     2           jwp(0,-1,3), wm,
     3           mmp(3),   mv2)
            print*," BB virtual 3,1: ",mv1/mvpm(3,1),mv2/mvmp(3,1)
            id = -2
            call boxline6(id, 
     1           psi(1,-1,3),psi(1,-1,4),p(0,3),p(0,4),-1,
     2           wm,jwp(0,-1,2),
     3           mpm(3),   mv1)
            id = -2
            call boxline6(id, 
     1           psi(1,-1,3),psi(1,-1,4),p(0,3),p(0,4),-1,
     2           wp,jwm(0,-1,2),
     3           mmp(3),   mv2)
            print*," BB virtual 3,2: ",mv1/mvpm(3,2),mv2/mvmp(3,2)
         endif
c  upper line box: eps1 = jwm(mu,-1,4), eps2 = wp
c  lower line box: eps1 = jwp(mu,-1,2), eps2 = wm
         mvpm(4,1) = -dotcc(jvwp(0,-1,is1,2),jwm(0,-1,4))
         mvmp(4,1) = -dotcc(jvwm(0,-1,is1,2),jwp(0,-1,4))
         mvpm(4,2) = -dotcc(jwp(0,-1,2),jvwm(0,-1,is3,4))
         mvmp(4,2) = -dotcc(jwm(0,-1,2),jvwp(0,-1,is3,4))
         if (ldebugvf) then
            id = -1
            call boxline6(id, 
     1           psi(1,-1,1),psi(1,-1,2),p(0,1),p(0,2),-1,
     2           jwm(0,-1,4), wp,
     3           mpm(4),   mv1)
            id = -1
            call boxline6(id, 
     1           psi(1,-1,1),psi(1,-1,2),p(0,1),p(0,2),-1,
     2           jwp(0,-1,4), wm,
     3           mmp(4),   mv2)
            print*," BB virtual 4,1: ",mv1/mvpm(4,1),mv2/mvmp(4,1)
            id = -1
            call boxline6(id, 
     1           psi(1,-1,3),psi(1,-1,4),p(0,3),p(0,4),-1,
     2           jwp(0,-1,2), wm,
     3           mpm(4),   mv1)
            id = -1
            call boxline6(id, 
     1           psi(1,-1,3),psi(1,-1,4),p(0,3),p(0,4),-1,
     2           jwm(0,-1,2), wp,
     3           mmp(4),   mv2)
            print*," BB virtual 4,2: ",mv1/mvpm(4,2),mv2/mvmp(4,2)
         endif
      endif

#endif
c  for the q^mu*q^nu/M_V^2 terms in the gauge boson propagators we need
      zpm = jj21p*jj43m*zm2i(3)
      zmp = jj21m*jj43p*zm2i(3)
c
      k = 1
      mat(k,-1,-1,6) = (mpm(2)+zpm)*prop_pm(3) + (mmp(3)+zmp)*prop_mp(3)
      if (ldebug) then
         m1 = tamp(3,-1,-1,k)+tamp(27,-1,-1,k)
         print*," 3+27,k=1          ",mat(k,-1,-1,6)/m1-1
      endif
      k = 2
      mat(k,-1,-1,6) = (mpm(1)-zpm)*prop_pm(3) + (mmp(4)-zmp)*prop_mp(3) 
      k = 3
      mat(k,-1,-1,6) = (mpm(4)-zpm)*prop_pm(3) + (mmp(1)-zmp)*prop_mp(3) 
      k = 4
      mat(k,-1,-1,6) = (mpm(3)+zpm)*prop_pm(3) + (mmp(2)+zmp)*prop_mp(3)
      k = 5
      zpm = jj21p*jj43m*zm2i(2)
      zmp = jj21m*jj43p*zm2i(2)
      mat(k,-1,-1,6) = 
     1    mpm(1)*(prop_pm(1)*clr(ifl(2,k),1,-1)*clr(ifl(4,k),1,-1) + 
     2            prop_pm(2)*clr(ifl(2,k),2,-1)*clr(ifl(4,k),2,-1) ) +
     1    mpm(2)*(prop_pm(1)*clr(ifl(2,k),1,-1)*clr(ifl(3,k),1,-1) + 
     2            prop_pm(2)*clr(ifl(2,k),2,-1)*clr(ifl(3,k),2,-1) ) +
     1    mpm(3)*(prop_pm(1)*clr(ifl(1,k),1,-1)*clr(ifl(4,k),1,-1) + 
     2            prop_pm(2)*clr(ifl(1,k),2,-1)*clr(ifl(4,k),2,-1) ) +
     1    mpm(4)*(prop_pm(1)*clr(ifl(1,k),1,-1)*clr(ifl(3,k),1,-1) + 
     2            prop_pm(2)*clr(ifl(1,k),2,-1)*clr(ifl(3,k),2,-1) ) 
     3    + zpm*prop_pm(2)*(clr(ifl(2,k),2,-1)-clr(ifl(1,k),2,-1))*
     4                     (clr(ifl(3,k),2,-1)-clr(ifl(4,k),2,-1))
      if (ldebug) then
         do k = 2,3
            m1 = tamp(1,-1,-1,k) + tamp(65,-1,-1,k)  
            print*," BB, k=",k,"          ",mat(k,-1,-1,6)/m1-1
         enddo
         k=4
         m1 = tamp(39,-1,-1,k) + tamp(63,-1,-1,k)  
         print*," BB, k=4           ",mat(k,-1,-1,6)/m1-1
         k=5
         m1 = tamp(1,-1,-1,k)+tamp(2,-1,-1,k)+tamp(28,-1,-1,k)+tamp(29,-1,-1,k)
     1      + tamp(30,-1,-1,k)+tamp(31,-1,-1,k)+tamp(3,-1,-1,k)+tamp(4,-1,-1,k)
         print*," BB, k=5           ",mat(k,-1,-1,6)/m1-1
      endif
      k = 6
      mat(k,-1,-1,6) = 
     1    mmp(1)*(prop_mp(1)*clr(ifl(2,k),1,-1)*clr(ifl(4,k),1,-1) + 
     2            prop_mp(2)*clr(ifl(2,k),2,-1)*clr(ifl(4,k),2,-1) ) +
     1    mmp(2)*(prop_mp(1)*clr(ifl(2,k),1,-1)*clr(ifl(3,k),1,-1) + 
     2            prop_mp(2)*clr(ifl(2,k),2,-1)*clr(ifl(3,k),2,-1) ) +
     1    mmp(3)*(prop_mp(1)*clr(ifl(1,k),1,-1)*clr(ifl(4,k),1,-1) + 
     2            prop_mp(2)*clr(ifl(1,k),2,-1)*clr(ifl(4,k),2,-1) ) +
     1    mmp(4)*(prop_mp(1)*clr(ifl(1,k),1,-1)*clr(ifl(3,k),1,-1) + 
     2            prop_mp(2)*clr(ifl(1,k),2,-1)*clr(ifl(3,k),2,-1) ) 
     3    + zmp*prop_mp(2)*(clr(ifl(2,k),2,-1)-clr(ifl(1,k),2,-1))*
     4                     (clr(ifl(3,k),2,-1)-clr(ifl(4,k),2,-1))
      if (ldebug) then
         m1 = tamp(1,-1,-1,k)+tamp(2,-1,-1,k)+tamp(28,-1,-1,k)+tamp(29,-1,-1,k)
     1      + tamp(30,-1,-1,k)+tamp(31,-1,-1,k)+tamp(3,-1,-1,k)+tamp(4,-1,-1,k)
         print*," BB, k=6           ",mat(k,-1,-1,6)/m1-1
      endif

      if (lbox) then
         do i = 1,2      ! 1 is for upper line, 2 for lower line QCD correction
            jj = 3 + 3*i ! stored in matv(...,6) and matv(...,9) respectively 
            k = 1
            matv(k,-1,-1,jj) = mvpm(2,i)*prop_pm(3)+mvmp(3,i)*prop_mp(3)
            k = 2
            matv(k,-1,-1,jj) = mvpm(1,i)*prop_pm(3)+mvmp(4,i)*prop_mp(3) 
            k = 3
            matv(k,-1,-1,jj) = mvpm(4,i)*prop_pm(3)+mvmp(1,i)*prop_mp(3) 
            k = 4
            matv(k,-1,-1,jj) = mvpm(3,i)*prop_pm(3)+mvmp(2,i)*prop_mp(3)
            k = 5
            matv(k,-1,-1,jj) = 
     1    mvpm(1,i)*(prop_pm(1)*clr(ifl(2,k),1,-1)*clr(ifl(4,k),1,-1) + 
     2            prop_pm(2)*clr(ifl(2,k),2,-1)*clr(ifl(4,k),2,-1) ) +
     1    mvpm(2,i)*(prop_pm(1)*clr(ifl(2,k),1,-1)*clr(ifl(3,k),1,-1) + 
     2            prop_pm(2)*clr(ifl(2,k),2,-1)*clr(ifl(3,k),2,-1) ) +
     1    mvpm(3,i)*(prop_pm(1)*clr(ifl(1,k),1,-1)*clr(ifl(4,k),1,-1) + 
     2            prop_pm(2)*clr(ifl(1,k),2,-1)*clr(ifl(4,k),2,-1) ) +
     1    mvpm(4,i)*(prop_pm(1)*clr(ifl(1,k),1,-1)*clr(ifl(3,k),1,-1) + 
     2            prop_pm(2)*clr(ifl(1,k),2,-1)*clr(ifl(3,k),2,-1) ) 
            k = 6
            matv(k,-1,-1,jj) = 
     1    mvmp(1,i)*(prop_mp(1)*clr(ifl(2,k),1,-1)*clr(ifl(4,k),1,-1) + 
     2            prop_mp(2)*clr(ifl(2,k),2,-1)*clr(ifl(4,k),2,-1) ) +
     1    mvmp(2,i)*(prop_mp(1)*clr(ifl(2,k),1,-1)*clr(ifl(3,k),1,-1) + 
     2            prop_mp(2)*clr(ifl(2,k),2,-1)*clr(ifl(3,k),2,-1) ) +
     1    mvmp(3,i)*(prop_mp(1)*clr(ifl(1,k),1,-1)*clr(ifl(4,k),1,-1) + 
     2            prop_mp(2)*clr(ifl(1,k),2,-1)*clr(ifl(4,k),2,-1) ) +
     1    mvmp(4,i)*(prop_mp(1)*clr(ifl(1,k),1,-1)*clr(ifl(3,k),1,-1) + 
     2            prop_mp(2)*clr(ifl(1,k),2,-1)*clr(ifl(3,k),2,-1) ) 
         enddo
      endif

c	if (lbox.and.id.ne.1) print*,"id=",id

c -------------------------------------------------------------------
c
c and now, finally, the pentagon contributions, i.e. two W"s emitted from the
c  same quark line

#ifdef WITH_NLO 	
      if (lpent .and. lfs(js1)) then     ! need new pentagon graphs for 12 line
         musq = -p21(4)
         if (NLO.eq.1) then
            do k = 1,3
               call J_virtual_tri_box_pent45(psi,p,1,2,sign, musq,
     1              wp, wm, zero, zero, k, -5, j5pm(0,k,is1,1))
               call J_virtual_tri_box_pent45(psi,p,1,2,sign, musq,
     1              wm, wp, zero, zero, k, -5, j5mp(0,k,is1,1))
            enddo
            icount1= icount1+1
         elseif (nlo.eq.-5) then
            do k = 1,3
               call J_virtual_tri_box_pent45(psi,p,1,2,sign, musq,
     1              wpp, wmp, xp, xm, k, -5, j5pm(0,k,is1,1) )
               call J_virtual_tri_box_pent45(psi,p,1,2,sign, musq,
     1              wmp, wpp, xm, xp, k, -5, j5mp(0,k,is1,1) )
            enddo
            icount1= icount1+1
         else
            do k = 1,3
               call J_virtual_tri_box_pent45(psi,p,1,2,sign, musq,
     1              wpp, wmp, xp, xm, k, 5, j5pm(0,k,is1,1) )
               call J_virtual_tri_box_pent45(psi,p,1,2,sign, musq,
     1              wmp, wpp, xm, xp, k, 5, j5mp(0,k,is1,1) )
            enddo
         endif

         lgc(js1) = .false.
         lfs(js1) = .false.
         if (lpt) then
	 
	 icountmax = 2**(N_points(0)-9)
	 
c  check numerical stability by gauge invariance
            mv1 = 0
            mv2 = 0
            xv1 = 0 
            xv2 = 0
            do k = 1,3
               mv1 = mv1+dotrc(p43,j5pm(0,k,is1,1))
               mv2 = mv2+dotrc(p43,j5mp(0,k,is1,1))
               xv1 = xv1+abs(dotrc(p43,j5pm(0,k,is1,1)))
               xv2 = xv2+abs(dotrc(p43,j5mp(0,k,is1,1)))
            enddo
c            if (sign(1).eq.2 .and. sign(3).eq.1) then
            if (.true.) then
               call gauge_check_5(j5pm(0,1,is1,1),p43,1,psi,p,1,2,wpp,wmp,0d0)
	         if(bad_gauge_sin) bad_gauge = .true.
               call gauge_check_5(j5mp(0,1,is1,1),p43,1,psi,p,1,2,wmp,wpp,0d0)
	         if(bad_gauge_sin) bad_gauge = .true.
            endif
c	    if (icb1.eq.0) then
c		print*, "criterion chosen: mv5/mvc.gt.1d-1"
c	    endif !icb1	
c            if (abs(mv1/xv1).gt.1d-3 .or. abs(mv2/xv2).gt.1d-3) then
	    if(bad_gauge) then
               icb1 = icb1+1
c               print*," pm upper: ",abs(mv1/xv1)
c               print*," mp upper: ",abs(mv2/xv2)
               lgc(js1) = .true.
            endif
            if (ldebug) then
               if (mod(icount1,icountmax).eq.0 ) then
                  print*," bad gauge check fraction ",
     1                             real(icb1)/real(icount1), js1
               endif
            endif
         endif   ! lpt
      endif

      if (lpent .and. lfs(js3)) then     ! need new pentagon graphs for 34 line
         musq = -p43(4)
         if (NLO.eq.1) then
            do k = 1,3
               call J_virtual_tri_box_pent45(psi,p,3,4,sign, musq,
     1              wp, wm, zero, zero, k, -5, j5pm(0,k,is3,2))
               call J_virtual_tri_box_pent45(psi,p,3,4,sign, musq,
     1              wm, wp, zero, zero, k, -5, j5mp(0,k,is3,2))
            enddo
            icount2= icount2+1
         elseif (nlo.eq.-5) then
            do k = 1,3
               call J_virtual_tri_box_pent45(psi,p,3,4,sign, musq,
     1              wpp, wmp, xp, xm, k, -5, j5pm(0,k,is3,2) )
                call J_virtual_tri_box_pent45(psi,p,3,4,sign, musq,
     1              wmp, wpp, xm, xp, k, -5, j5mp(0,k,is3,2) )
            enddo
            icount2= icount2+1
         else
            do k = 1,3
               call J_virtual_tri_box_pent45(psi,p,3,4,sign, musq,
     1              wpp, wmp, xp, xm, k, 5, j5pm(0,k,is3,2) )
              call J_virtual_tri_box_pent45(psi,p,3,4,sign, musq,
     1              wmp, wpp, xm, xp, k, 5, j5mp(0,k,is3,2) )
            enddo
         endif

         lgc(js3) = .false.
         lfs(js3) = .false.
         if (lpt) then
            mv1 = 0
            mv2 = 0
            xv1 = 0 
            xv2 = 0
            do k = 1,3
               mv1 = mv1+dotrc(p21,j5pm(0,k,is3,2))
               mv2 = mv2+dotrc(p21,j5mp(0,k,is3,2))
               xv1 = xv1+abs(dotrc(p21,j5pm(0,k,is3,2)))
               xv2 = xv2+abs(dotrc(p21,j5mp(0,k,is3,2)))
            enddo
            if (.true.) then
c            if (sign(1).eq.2 .and. sign(3).eq.1) then
               call gauge_check_5(j5pm(0,1,is3,2),p21,1,psi,p,3,4,wpp,wmp,0d0)
	         if(bad_gauge_sin) bad_gauge = .true.
               call gauge_check_5(j5mp(0,1,is3,2),p21,1,psi,p,3,4,wmp,wpp,0d0)
	         if(bad_gauge_sin) bad_gauge = .true.
            endif
c            if (abs(mv1/xv1).gt.1d-2 .or. abs(mv2/xv2).gt.1d-2) then
	    if(bad_gauge) then
               icb2 = icb2+1
c               print*," pm lower: ",abs(mv1/xv1)
c               print*," mp lower: ",abs(mv2/xv2)
               lgc(js3) = .true.
            endif
            if (ldebug) then
               if (mod(icount2,icountmax).eq.0 ) then
                  print*," bad gauge check fraction ",
     1                          real(icb2)/real(icount2), js3
               endif
            endif
         endif  !lpt
      endif
#endif

      fac = clr(3,3,-1)**2
c upper line:
      do isig3 = 1,-1,-2
         call ket2c(psi(1,-1,1),.true.,p(0,1),-1,
     1        p43,jqq(0,isig3,2),bkjqq(1,-1,isig3,1),dummy)
         call bra2c(psi(1,-1,2),.true.,p(0,2),-1,
     1        p43,jqq(0,isig3,2),bkjqq(1,-1,isig3,2),dummy)
c  eps1=j43,eps2=wp,eps3=wm
         m5(1,3) = -s1c(psiwm(1,2),wp,.true.,-1,bkjqq(1,-1,isig3,1))
c  eps1=wp,eps2=j43,eps3=wm
         m5(2,3) = -s1c(psiwm(1,2),jqq(0,isig3,2),.true.,-1,psiwp(1,1))
c  eps1=wp,eps2=wm,eps3=j43
         m5(3,3) = -s1c(bkjqq(1,-1,isig3,2),wm,.true.,-1,psiwp(1,1))
c  eps1=j43,eps2=wm,eps3=wp
         m5(1,4) = -s1c(psiwp(1,2),wm,.true.,-1,bkjqq(1,-1,isig3,1))
c  eps1=wm,eps2=j43,eps3=wp
         m5(2,4) = -s1c(psiwp(1,2),jqq(0,isig3,2),.true.,-1,psiwm(1,1))
c  eps1=wm,eps2=wp,eps3=j43
         m5(3,4) = -s1c(bkjqq(1,-1,isig3,2),wp,.true.,-1,psiwm(1,1))
         if (lpent) then
            do j = 1,3
               mv5(j,3) = -dotcc(j5pm(0,j,is1,1),jqq(0,isig3,2))
c     1                    -cvirtc*m5(j,3)
               mv5(j,4) = -dotcc(j5mp(0,j,is1,1),jqq(0,isig3,2))
c     1                    -cvirtc*m5(j,4)
            enddo
         endif

         do k = 1,4
            kk = mod(k+2,4)
            if (k.eq.2) kk = 4 
            mat(k,-1,isig3,7) = fac*m5(2,ifl(1,k))*propt(-1,isig3,kk,2)+
     1         fac*( m5(1,ifl(2,k))+m5(3,ifl(1,k)) )*propt(-1,isig3,k,2)
            if (lpent) matv(k,-1,isig3,7) = 
     1         fac*mv5(2,ifl(1,k))*propt(-1,isig3,kk,2)+
     2         fac*( mv5(1,ifl(2,k))+mv5(3,ifl(1,k)) )*propt(-1,isig3,k,2)
         enddo
      enddo
      do k = 5,6
         mat(k,-1,-1,7) = propt(-1,-1,k,2) * fac *
     1                    ( m5(1,ifl(2,k))+m5(3,ifl(1,k)) )
         if (lpent) matv(k,-1,-1,7) = propt(-1,-1,k,2) * fac *
     1                    ( mv5(1,ifl(2,k))+mv5(3,ifl(1,k)) )
         if (ldebug) then
            if (k.eq.5) then
               m1 = tamp(13,-1,-1,k)+tamp(40,-1,-1,k)
               print*," Pent, k=5         ",mat(k,-1,-1,7)/m1-1
            elseif (k.eq.6) then
               m1 = tamp(19,-1,-1,k)+tamp(46,-1,-1,k)
               print*," Pent, k=6         ",mat(k,-1,-1,7)/m1-1
            endif
         endif
      enddo
c lower line:
      do isig1 = 1,-1,-2
         call ket2c(psi(1,-1,3),.true.,p(0,3),-1,
     1        p21,jqq(0,isig1,1),bkjqq(1,isig1,-1,3),dummy)
         call bra2c(psi(1,-1,4),.true.,p(0,4),-1,
     1        p21,jqq(0,isig1,1),bkjqq(1,isig1,-1,4),dummy)
c  eps1=j21,eps2=wp,eps3=wm
         m5(1,3) = -s1c(psiwm(1,4),wp,.true.,-1,bkjqq(1,isig1,-1,3))
c  eps1=wp,eps2=j21,eps3=wm
         m5(2,3) = -s1c(psiwm(1,4),jqq(0,isig1,1),.true.,-1,psiwp(1,3))
c  eps1=wp,eps2=wm,eps3=j21
         m5(3,3) = -s1c(bkjqq(1,isig1,-1,4),wm,.true.,-1,psiwp(1,3))
c  eps1=j21,eps2=wm,eps3=wp
         m5(1,4) = -s1c(psiwp(1,4),wm,.true.,-1,bkjqq(1,isig1,-1,3))
c  eps1=wm,eps2=j21,eps3=wp
         m5(2,4) = -s1c(psiwp(1,4),jqq(0,isig1,1),.true.,-1,psiwm(1,3))
c  eps1=wm,eps2=wp,eps3=j21
         m5(3,4) = -s1c(bkjqq(1,isig1,-1,4),wp,.true.,-1,psiwm(1,3))
         if (lpent) then
            do j = 1,3
               mv5(j,3) = -dotcc(j5pm(0,j,is3,2),jqq(0,isig1,1))
c     1                    -cvirtc*m5(j,3)
               mv5(j,4) = -dotcc(j5mp(0,j,is3,2),jqq(0,isig1,1))
c     1                    -cvirtc*m5(j,4)
            enddo
         endif
         do k = 1,4
            kk = k+1  - 2*mod(k+1,2)
            mat(k,isig1,-1,8) = fac*m5(2,ifl(3,k))*propt(isig1,-1,kk,1)+
     1         fac*( m5(1,ifl(4,k))+m5(3,ifl(3,k)) )*propt(isig1,-1,k,1)
            if (lpent) matv(k,isig1,-1,8) = 
     1         fac*mv5(2,ifl(3,k))*propt(isig1,-1,kk,1)+
     2         fac*( mv5(1,ifl(4,k))+mv5(3,ifl(3,k)) )*propt(isig1,-1,k,1)
         enddo
      enddo
      do k = 5,6
         mat(k,-1,-1,8) = propt(-1,-1,k,1) * fac *
     1                    ( m5(1,ifl(4,k))+m5(3,ifl(3,k)) )
         if (lpent) matv(k,-1,-1,8) = propt(-1,-1,k,1) * fac *
     1                    ( mv5(1,ifl(4,k))+mv5(3,ifl(3,k)) )
         if (ldebug) then
            if (k.eq.5) then
               m1 = tamp(19,-1,-1,k)+tamp(46,-1,-1,k)
               print*," 19+46,k=5         ",mat(k,-1,-1,8)/m1-1
            elseif (k.eq.6) then
               m1 = tamp(13,-1,-1,k)+tamp(40,-1,-1,k)
               print*," 13+40,k=6         ",mat(k,-1,-1,8)/m1-1
            endif
         endif
      enddo

c
c ------------------------------------------------------------
c
c sum the graphs, square them and map them onto uucc, uuss etc.

      if (lpt) then
         xgc1 = real(icount1+1)/real(icount1-icb1+1)
         xgc2 = real(icount2+1)/real(icount2-icb2+1)
      endif

c
      do k = 1,6
         res(k) = 0
         resv(k) = 0
         tree(k) = 0
         do isig1 = -1,1,2
            do isig3 = -1,1,2
               mm(k,isig1,isig3) = 0
               do i = 1,8
                  mm(k,isig1,isig3) = mm(k,isig1,isig3) + 
     1                                mat(k,isig1,isig3,i)
               enddo
               res(k) = res(k) + dreal(mm(k,isig1,isig3))**2
     &                         + dimag(mm(k,isig1,isig3))**2
               if (lnlo) then
                  mv12(k,isig1,isig3) = 
     1         +  matv(k,isig1,isig3,2) + matv(k,isig1,isig3,4) + 
     2            matv(k,isig1,isig3,6) + matv(k,isig1,isig3,7) 
                  mv34(k,isig1,isig3) =
     1          + matv(k,isig1,isig3,3) + matv(k,isig1,isig3,5) +
     2            matv(k,isig1,isig3,9) + matv(k,isig1,isig3,8)

                  if (lpt .and. isig1.eq.-1) then   ! Pentagon contributes for q_L only
                     if (lgc(js1)) then
                        mv12(k,isig1,isig3)=0
                     else
                        mv12(k,isig1,isig3)=mv12(k,isig1,isig3)*xgc1
                     endif
                  endif
                  if (lpt .and. isig3.eq.-1) then   ! Pentagon contributes for q_L only
                     if (lgc(js3)) then
                        mv34(k,isig1,isig3)=0
                     else
                        mv34(k,isig1,isig3)=mv34(k,isig1,isig3)*xgc2
                     endif
                  endif

c  add Born type term and multiply by F_q = alphas*C_2/4pi
c  the factor pi^2/3+9/2 for the born term is after adding the subtraction term
c  and the counter term for the renormalization of the pdfs
                  if (nlo.gt.0) then
                     mv12(k,isig1,isig3) = als(1,1)*c2o4pi*
     1                ( mv12(k,isig1,isig3) + mm(k,isig1,isig3)*cvirt )
                     mv34(k,isig1,isig3) = als(2,1)*c2o4pi*
     1                ( mv34(k,isig1,isig3) + mm(k,isig1,isig3)*cvirt )
                  else
                     mv12(k,isig1,isig3) = 
     1                    als(1,1)*c2o4pi*mv12(k,isig1,isig3)
                     mv34(k,isig1,isig3) = 
     1                    als(2,1)*c2o4pi*mv34(k,isig1,isig3)
                  endif
                  resv(k) = resv(k) + 2*dreal(
     1                 mm(k,isig1,isig3)   *
     1              conjg( mv12(k,isig1,isig3)+mv34(k,isig1,isig3) )  )
               endif
            enddo
         enddo  
         tree(k) = res(k) * 9D0 
         if (nlo.ge.0) then
            res(k) = (res(k)+resv(k))*9d0      ! 9 is the color sum factor
         else
            res(k) = resv(k)*9d0               ! 9 is the color sum factor
         endif
      enddo
      if (ldebug) then
         print*,"uucc",uucc,res(1)/uucc
         print*,"uuss",uuss,res(2)/uuss
         print*,"ddcc",ddcc,res(3)/ddcc
         print*,"ddss",ddss,res(4)/ddss
         print*,"udsc",udsc,res(5)/udsc
         print*,"ducs",ducs,res(6)/ducs
      endif
c      if (ldebugv) 
c      read(*,"(a)") text
      if (lwarn) then
         resv(1) = uucc
         resv(2) = uuss
         resv(3) = ddcc
         resv(4) = ddss
         resv(5) = udsc
         resv(6) = ducs
         do k = 1,6
            if ((res(k)/resv(k)-1).gt.1d-7) then
               if (.not.ldebug) then
                  print*," k =",k,":",res(k)/resv(k)-1
                  ldebug = .true.
                  lerror = .true.
                  goto 701
               endif
            endif
         enddo
         if (lerror) then
            print*," "
            lerror = .false.
            ldebug = .false.
         endif
      endif
c eliminate processes with photon virtuality below cutoff
      if ( abs(p21(4)).lt.qsqAmin .or.  
     1     abs(p43(4)).lt.qsqAmin ) then
         do k = 1,4
            res(k) = 0
            tree(k) = 0
         enddo
      endif

      uucc = res(1)
      uuss = res(2)
      ddcc = res(3)
      ddss = res(4)
      udsc = res(5)
      ducs = res(6)

      if (ldoblha) then
        if (blha_amptype(blha_curproc).eq.0) then ! tree
          blha_amp(1)=blha_amp(1)+
     &      res(blha_idsubproc(blha_cursubproc,blha_curproc))

        else if (blha_amptype(blha_curproc).eq.1) then ! loop
          lnmusqp1p2 = log(blha_scale**2/(-p21(4)))
          lnmusqp3p4 = log(blha_scale**2/(-p43(4)))

          blha_amp(1)=blha_amp(1)+
     &      blha_CF/(2d0*pi)*2d0*
     &       (als(1,1)+als(2,1))*(-1d0)*
     &       tree(blha_idsubproc(blha_cursubproc,blha_curproc))
          blha_amp(2)=blha_amp(2)+
     &      1d0/(2d0*pi)*2d0*
     &       ( als(1,1)*(-blha_gammaQuark-blha_CF*lnmusqp1p2)
     &        +als(2,1)*(-blha_gammaQuark-blha_CF*lnmusqp3p4) )*
     &       tree(blha_idsubproc(blha_cursubproc,blha_curproc))
c this is more complicated than necessary, 
c but should work for other VBF processes as well
          blha_amp(3)= blha_amp(3)+
     &      (res(blha_idsubproc(blha_cursubproc,blha_curproc))     ! tree + loop + I
     &       -tree(blha_idsubproc(blha_cursubproc,blha_curproc)))  ! -tree
     &       * blha_CF/c2                                          ! correct for Nc ?= inf
     &      + (-(als(1,1)+als(2,1))/2d0/pi*                        ! -I
     &          (blha_CF*cvirt+
     &            2d0*(blha_KQuark+blha_gammaQuark+
     &                 blha_tgammaQuark+(pi**2/6d0-1)*blha_CF) )
     &         +als(1,1)/pi*                                       ! +(mu^2/2p1.p2)^epsilon
     &         (-blha_gammaQuark*lnmusqp1p2-blha_CF/2d0*lnmusqp1p2**2)
     &         +als(2,1)/pi*                                       ! +(mu^2/2p3.p4)^epsilon
     &         (-blha_gammaQuark*lnmusqp3p4-blha_CF/2d0*lnmusqp3p4**2)
     &        )*tree(blha_idsubproc(blha_cursubproc,blha_curproc))
          blha_amp(4)= blha_amp(4)+
     &       tree(blha_idsubproc(blha_cursubproc,blha_curproc))

        else if (blha_amptype(blha_curproc).eq.2) then ! cctree
          ! 1-2
          call BLHA_cctree(1,2,
     &        -blha_CF*
     &         tree(blha_idsubproc(blha_cursubproc,blha_curproc))
     &    )
          ! 3-4
          call BLHA_cctree(3,4,
     &        -blha_CF*
     &         tree(blha_idsubproc(blha_cursubproc,blha_curproc))
     &    )

        else 
          call BLHA_amptypeerror(blha_amptype(blha_curproc),__FILE__,__LINE__)
        endif
      endif

      Lold = L
      return
      end

c ------------------------------------------------------------


#ifdef WITH_NLO 
c
c----------  boxlinec(isigmax, ... ) ----------------------------------------
c
      subroutine boxlinec(isigmax, psi1,psi2,k1,k2,
     1                    ld0calc, jext, ivmu, 
     2                    jborn,    jvirt)
      implicit none
      integer isigmax, ivmu
      logical ld0calc
      double precision k1(0:3),k2(0:3)
      double complex psi1(2,-1:1),psi2(2,-1:1), jext(0:5), 
     1               jborn(0:5,-1:1),jvirt(0:5,-1:1)
c
c  variant of the boxline routines which returns the virtual corrections 
c  without contraction of one of the external gauge bosons, i.e. in form of a
c  current. Also, results for both fermion chiralities are returned for isgmax = +1
c 
c  determine the finite virtual corrections to the fundamental box-type
c  graph, i.e. the sum of all virtual corrections to the Born diagram
c
c    psi1(k1)    ---->---------->-------------->--  psi2(k2)
c                         S              S
c                         S              S
c                         S              S
c                     q1  V           q2 V
c                         S              S
c                         S              S
c
c                        mu           jext(q2)      for ivmu = 1 
c                     jext(q1)          mu          for ivmu = 2
c
c          ====       eps1(q1)        eps2(q2)
c
c   note: k1 = k2+q1+q2, k1.k1 = 0 = k2.k2, q1.q1,q2.q2 arbitrary
c
c  INPUT:      isigmax         -1: only Lefthanded fermions, +1: return L and R
c              psi1(2,-1:1)    the incoming Weyl spinor for both chiralities
c              psi2(2,-1:1)    the outgoing Weyl spinor, psi^+(k2)
c       k1(0:3),k2(0:3)        the corresponding momenta in direction of
c                              the fermion line
c       ld0calc                call BCD_fill_c if true
c       jext(0:3)              specified polarization vector/current" 
c                              corresponding to 1 external EW boson
c       jext(4:5)              outflowing momenta of this EW bosons
c                              q1(0:3) or q2(0:3) in DHELAS notation
c       jborn(0:5,-1:1)        the current of free index at ivmu for the 
c                              corresponding Born diagram and 2 quark chiralities
c
c  OUTPUT:   jvirt(0:5,-1:1)   the amplitude for the finite part of all 
c                              virtual 1-loop QCD corrections for free index at ivmu
c
      logical ldebug
      parameter (ldebug=.false.)
common block for D_ij etc. information via bcd_fill_c
      double complex D0c, Dijc(3,13), Teps1c, Teps2c, Tbornc,Tg1c,Tg2c
      double precision sc, tc, q1sqc, q2sqc
      common /bcd_qqv_c/ D0c,Dijc,Teps1c,Teps2c,Tbornc,Tg1c,Tg2c, 
     1                   sc,tc,q1sqc,q2sqc
C   local variables
      integer mu, nu, isig
      double precision s,t,mom(0:3),dotrr
      external dotrr
      double precision q1(0:4),q2(0:4)
      double complex eps1(0:3),eps2(0:3),mborn,mvirt
      double complex me1, me2, mq, d0t, s1c, s1r, sc3, z
      external d0t, s1c, s1r, sc3
      double complex ce1,ce2,cq,cborn
      double complex d0, d11,d12,d13, d21,d22,d23,d24,d25,d26,d27
      double complex d31,d32,d33,d34,d35,d36,d37,d38,d39,d310,
     1               d311,d312,d313
      double complex eps1k1, eps1k2, eps1q2, eps1eps2, 
     1               eps2k1, eps2q1, eps2k2, eps1q1, eps2q2,
     2               gce1, gce2, gcq

c variables for monitoring gauge invariance:
      double complex eps1g(0:3),eps2g(0:3),mborng,mvirtg
      double complex me1g, me2g, mqg
      double complex ce1g,ce2g,cqg
      double complex eps1gk1, eps1gk2, eps1gq2, eps1geps2g, 
     1               eps2gk1, eps2gq1, eps2gk2, eps1gq1, eps2gq2,
     2               gce1g, gce2g, gcqg
c
      integer icountmax


!      character*50 gridname1, gridname2
!      common /bmainpr/ iterations1, iterations2, n2max1, n2max2,
!     2                 lflavr, gridname1, gridname2
      logical breset
      parameter(breset = .true.)

#include "VBFNLO/utilities/global.inc"
	
      logical gimprove
      parameter (gimprove =.true.)
      integer bcount
      data bcount / 0/
      integer gcount
      data gcount / 0/

c transfer momentum information of EW bosons to q1,q2
      if (ivmu.eq.1) then     ! eps2 and q2 are specified, eps1 = mu = free index
         q2(0) = dreal(jext(4))
         q2(1) = dreal(jext(5))
         q2(2) = dimag(jext(5))
         q2(3) = dimag(jext(4))
         do mu = 0,3
            eps2(mu) = jext(mu)
            mom(mu) = k2(mu)+q2(mu)
            q1(mu) = k1(mu)-mom(mu)

         enddo
cfc         print*, "curr1"
cfc                print*, "k1", k1
cfc                print*, "k2", -k2
cfc                print*, "q1", -q1
cfc                print*, "q2", -q2

      elseif (ivmu.eq.2) then     ! eps1 and q1 are specified, eps2 = mu = free index
         q1(0) = dreal(jext(4))
         q1(1) = dreal(jext(5))
         q1(2) = dimag(jext(5))
         q1(3) = dimag(jext(4))
         do mu = 0,3
            eps1(mu) = jext(mu)
            mom(mu) = k1(mu)-q1(mu)
            q2(mu) = mom(mu)-k2(mu)
         enddo
cfc         print*, "curr2_old"
cfc                print*, "k1", k1
cfc                print*, "k2", -k2
cfc                print*, "q1", -q1
cfc                print*, "q2", -q2

      endif
      q1(4) = q1(0)**2-q1(1)**2-q1(2)**2-q1(3)**2
      q2(4) = q2(0)**2-q2(1)**2-q2(2)**2-q2(3)**2
c
c for gauge test:  replace eps(q) by q
      if (ivmu.eq.2) then
	do mu = 0,3
	   eps1g(mu) = q1(mu)	
	   eps2g(mu) = eps2(mu)
	enddo 
      elseif (ivmu.eq.1) then	      
	do mu = 0,3
	   eps2g(mu) = q2(mu)
	   eps1g(mu) = eps1(mu)
	enddo 
      else
        print*,"invalid value for ivmu"	       
      endif	

      if (ld0calc) call BCD_fill_c(k1,k2,q1,q2)
      d0 = d0c

      d11 = dijc(1,1)
      d12 = dijc(1,2)
      d13 = dijc(1,3)
      
      d21 = dijc(2,1)
      d22 = dijc(2,2)
      d23 = dijc(2,3)
      d24 = dijc(2,4)
      d25 = dijc(2,5)
      d26 = dijc(2,6)
      d27 = dijc(2,7)

      d31 = dijc(3,1)
      d32 = dijc(3,2)
      d33 = dijc(3,3)
      d34 = dijc(3,4)
      d35 = dijc(3,5)
      d36 = dijc(3,6)
      d37 = dijc(3,7)
      d38 = dijc(3,8)
      d39 = dijc(3,9)
      d310 = dijc(3,10)
      d311 = dijc(3,11)
      d312 = dijc(3,12)
      d313 = dijc(3,13)
c

Check consistency of kinematics
      if (ldebug) then
         print*," "
         print*," comparison of input to BOXLINE for ivmu = ",ivmu
         print("(a,3f16.4)")," q1sq = ",q1sqc,q1(4),q1sqc/q1(4)
         print("(a,3f16.4)")," q2sq = ",q2sqc,q2(4),q2sqc/q2(4)
         s = -2*dotrr(k1,k2)
         t =  2*dotrr(k2,q2)+q2(4)
         print("(a,3f16.4)")," s = ",sc,s,sc/s
         print("(a,3f16.4)")," t = ",tc,t,tc/t
c         do mu = 0,3
c            if (abs(k1(mu)-ak1(mu,id)).gt.(1d-10*k1(0)) .or.
c     1          abs(k2(mu)-ak2(mu,id)).gt.(1d-10*k2(0)) .or.
c     1          abs(q1(mu)-aq1(mu,id)).gt.(1d-10*k2(0)) .or.
c     1          abs(q2(mu)-aq2(mu,id)).gt.(1d-10*k2(0)) ) then
c               print*," WARNING: momentum mismatch in boxline, mu =",mu
c               print*," k1: ",k1(mu)/(ak1(mu,id)+1d-80)
c               print*," k2: ",k2(mu)/(ak2(mu,id)+1d-80)
c               print*," q1: ",q1(mu)/(aq1(mu,id)+1d-80)
c               print*," q2: ",q2(mu)/(aq2(mu,id)+1d-80)
c            endif
c         enddo
         print*," input psi1 ",psi1
         print*," input psi2 ",psi2
         print*," input jext ",jext
      endif

c  get matrix elements 
c
c  me1 = psibar(2) eps1slash psi(1)
c  me2 = psibar(2) eps2slash psi(1)
c  mq  = psibar(2) (q1-q2)slash psi(1)
      do mu = 0,3
         if (ivmu.eq.1) then
            do nu = 0,3
               eps1(nu)  = 0d0
               eps1g(nu) = 0d0
            enddo
            eps1(mu)  = -1d0
            eps1g(mu) = -1d0
            if (mu.eq.0) then
	    	eps1(mu)  = 1d0
	    	eps1g(mu) = 1d0
	    endif	
         elseif (ivmu.eq.2) then
            do nu = 0,3
               eps2(nu)  = 0d0
               eps2g(nu) = 0d0
            enddo
            eps2(mu) =  -1d0
            eps2g(mu) = -1d0
            if (mu.eq.0) then
	        eps2(mu)  = 1d0
	        eps2g(mu) = 1d0
	    endif	
         endif
	 
         do isig = -1,isigmax,2
            
            me1 = s1c(psi2(1,isig),eps1,.true.,isig,psi1(1,isig))
            me2 = s1c(psi2(1,isig),eps2,.true.,isig,psi1(1,isig))
            mq  = 2*s1r(psi2(1,isig),q1,.true.,isig,psi1(1,isig))
            mborn = jborn(mu,isig)
            	    
            if (ldebug) then
               print*," check amplitudes for mu =",mu," isig =",isig
               print*," check M_q: ",
     1         -mq/s1r(psi2(1,isig),q2,.true.,isig,psi1(1,isig)),
     1         " = 2?"
               t = mom(0)**2-mom(1)**2-mom(2)**2-mom(3)**2
               z = sc3(psi2(1,isig),eps2,mom,eps1,psi1(1,isig),isig)/t
               print*," M_born:    ",mborn/z," = 1? Abs(ratio) = ",
     &              dreal(mborn/z),abs(mborn/z)
            endif
c
c prepare the coefficients
            s = sc
            t = tc

      eps1k1 = eps1(0)*k1(0)-eps1(1)*k1(1)-eps1(2)*k1(2)-eps1(3)*k1(3)
      eps1k2 = eps1(0)*k2(0)-eps1(1)*k2(1)-eps1(2)*k2(2)-eps1(3)*k2(3)
      eps1q2 = eps1(0)*q2(0)-eps1(1)*q2(1)-eps1(2)*q2(2)-eps1(3)*q2(3)
      eps1eps2 = eps1(0)*eps2(0)-eps1(1)*eps2(1)-
     1           eps1(2)*eps2(2)-eps1(3)*eps2(3)
      eps2k1 = eps2(0)*k1(0)-eps2(1)*k1(1)-eps2(2)*k1(2)-eps2(3)*k1(3)
      eps2k2 = eps2(0)*k2(0)-eps2(1)*k2(1)-eps2(2)*k2(2)-eps2(3)*k2(3)
      eps2q1 = eps2(0)*q1(0)-eps2(1)*q1(1)-eps2(2)*q1(2)-eps2(3)*q1(3)

      eps1q1 = eps1(0)*q1(0)-eps1(1)*q1(1)-eps1(2)*q1(2)-eps1(3)*q1(3)
      eps2q2 = eps2(0)*q2(0)-eps2(1)*q2(1)-eps2(2)*q2(2)-eps2(3)*q2(3)
c      
c
      ce1 = -2*(d12+d24)*eps2k2*(q1sqc+q2sqc-3*s-4*t) - 
     1      2*(d12-d24)*eps2q1*(q2sqc-t) 
     2  +4*(- d11*eps2k2*s - d12*eps2k1*t + d13*eps2k2*
     3  (q2sqc-s-t) + d13*eps2q1*q2sqc - d21*eps2k2*s 
     4 - d22*eps2k2*t - d22*eps2q1*q2sqc + d23*eps2k2*q1sqc
     5 + d25*eps2k2*(q2sqc-s-2*t) 
     6 - d26*eps2k2*(q2sqc-s-t)+d26*eps2q1*t+2*d27*eps2q1
     7 - d32*eps2k2*q2sqc - d34*eps2k2*(q2sqc-t) 
     8 + d35*eps2k2*(q2sqc-s-t)+d36*eps2k2*(2*q2sqc-t) 
     9 + d37*eps2k2*q1sqc + d38*eps2k2*(q1sqc+q2sqc-s)
     & - d39*eps2k2*q1sqc 
     1 - d310*eps2k2*(q1sqc+2*q2sqc-2*s-t)
     2 - 4*d311*eps2k2 + 6*d312*eps2k2 + 2*d313*eps2q1)

      ce2= -2*(d12+d24)*(eps1k2*(q1sqc+q2sqc-s-2*t)
     1  +eps1q2*(q2sqc-s-3*t)) + 
     2  4*(d13*eps1k2*(q1sqc-2*s-t) - 
     3  d13*eps1q2*(2*s+t) + d22*eps1k2*q2sqc
     4  + d22*eps1q2*q2sqc - d23*eps1k2*t + 
     5    d23*eps1q2*(q1sqc-t) - d24*eps1k1*q2sqc 
     6  + d25*eps1k2*q1sqc + d25*eps1k1*(q2sqc-2*s-t) 
     7  + d26*eps1k2*t - d26*eps1k1*(q1sqc-s) - 
     8  2*d27*eps1q2 + d33*eps1k2*q1sqc + d33*eps1q2*q1sqc
     9  + d37*eps1k1*(q2sqc-s-t) + d38*eps1k1*q2sqc
     &  - d39*eps1k1*(q1sqc+q2sqc-s) - d310*eps1k1*(q2sqc-t)
     1  +2*d311*eps1k2+2*d312*eps1q2-6*d313*eps1k1)

C NOTE: cq = 1/2 * cq1 in my reduce code
      cq = (d12+d24)*eps1eps2*s + 2*(4*d12*eps2k2*eps1k2 + 
     1 3*d12*eps2k2*eps1q2+d12*eps2q1*eps1k2-4*d13*eps2k2*eps1k2-
     2 2*d13*eps2k2*eps1q2-2*d13*eps2q1*eps1k2-d13*eps1eps2*s+
     3 2*d22*eps2k2*eps1q2-d22*eps1eps2*t-2*d23*eps2q1*eps1k2-
     4 2*d23*eps2q1*eps1q2-d23*eps1eps2*t+6*d24*eps2k2*eps1k2+
     5 3*d24*eps2k2*eps1q2+d24*eps2q1*eps1k2-6*d25*eps2k2*eps1k2-
     6 2*d25*eps2k2*eps1q2-2*d25*eps2q1*eps1k2-d25*eps1eps2*s-
     7 4*d26*eps2k2*eps1q2+4*d26*eps2q1*eps1k2+2*d26*eps2q1*eps1q2+ 
     8 d26*eps1eps2*(s+2*t)-d32*eps1eps2*q2sqc+
     9 d33*eps1eps2*q1sqc+2*d34*eps2k2*eps1k2-2*d35*eps2k2*eps1k2+
     & 2*d36*eps2k2*eps1q2+d36*eps1eps2*(q2sqc-t)-
     1 2*d37*eps2q1*eps1k2+d37*eps1eps2*(q2sqc-s-t)+
     2 2*d38*eps2q1*eps1q2+d38*eps1eps2*(q1sqc+2*q2sqc-s) -
     3 2*d39*eps2q1*eps1q2-d39*eps1eps2*(2*q1sqc+q2sqc-s) -
     4 2*d310*eps2k2*eps1q2+2*d310*eps2q1*eps1k2-d310*eps1eps2*
     5 (2*q2sqc-s-2*t)+4*d312*eps1eps2-4*d313*eps1eps2)

      cborn = -2*((d36+d37-2*d310)*(q2sqc-t) + 
     1    d38*(q1sqc+2*q2sqc) - d39*(2*q1sqc+q2sqc) )
     2    -2*(d0+d11+d12-2*d13+d24-d25+d26-d37-d38+d39+d310)*s
     3    -2*(-(d22+d23-2*d26)*t + 2*d27 - d32*q2sqc
     4    +d33*q1sqc+ 6*(d312-d313))

      gce1 = -2*eps2q2*( d12*(2*q1sqc+q2sqc-s-3*t) - 
     1       2*(d13-d22)*(q1sqc+q2sqc-s-t)
     2     - d24*(q2sqc-s-t)-2*d26*q1sqc-4*d27-4*d312)

      gce2 = 2*eps1q1*(d12*(q2sqc+t)-2*d13*(q2sqc+s) + 
     1      2*d22*q2sqc + 2*d23*(q1sqc+q2sqc-s-2*t)
     2     +d24*(q2sqc+t)+2*d25*(-q2sqc+s+t)-2*d26*(q2sqc-t)
     3     -4*d27+2*d33*q1sqc + 4*d313)

      gcq = 4*eps1q1*eps2q2*(d12-d13+d22-d23+d38-d39) +
     1        2*eps1q1*eps2k2*(d12-4*d23+d24+2*d26-2*d37+2*d310) -
     2        4*eps1q1*eps2q1*(d23-d26+d33-d39) +
     3        2*eps2q2*eps1k2*(3*d12-4*d13+4*d22+d24-2*d25-
     4                                     2*d26+2*d36-2*d310) +
     5        4*eps2q2*eps1q2*(d12-d13+2*d22-2*d26+d32-d38)
     
            ce1 = ce1 + gce1
            ce2 = ce2 + gce2
            cq = cq + gcq
c
c  now add the terms from the vertex and propagator corrections
            ce1 = ce1 + 2*eps2k2*Teps2c + 2*eps2q2*Tg2c
            ce2 = ce2 + 2*eps1k1*Teps1c - 2*eps1q1*Tg1c
            cborn = cborn*t - Tbornc
c  and the final result
            mvirt = me1*ce1 + me2*ce2 + mq*cq + mborn*cborn
c
c -------------------------------------------
c
c "same" for gauge check: 
c
      if (gimprove) then      
     	
       t = mom(0)**2-mom(1)**2-mom(2)**2-mom(3)**2
     
       me1g = s1c(psi2(1,isig),eps1g,.true.,isig,psi1(1,isig))
       me2g = s1c(psi2(1,isig),eps2g,.true.,isig,psi1(1,isig))
       mqg = mq
       mborng = sc3(psi2(1,isig),eps2g,mom,eps1g,psi1(1,isig),isig)/t	       
      
       eps1gk1 = eps1g(0)*k1(0)-eps1g(1)*k1(1)-eps1g(2)*k1(2)-eps1g(3)*k1(3)
       eps1gk2 = eps1g(0)*k2(0)-eps1g(1)*k2(1)-eps1g(2)*k2(2)-eps1g(3)*k2(3)
       eps1gq2 = eps1g(0)*q2(0)-eps1g(1)*q2(1)-eps1g(2)*q2(2)-eps1g(3)*q2(3)
       eps1geps2g = eps1g(0)*eps2g(0)-eps1g(1)*eps2g(1)-
     1             eps1g(2)*eps2g(2)-eps1g(3)*eps2g(3)
       eps2gk1 = eps2g(0)*k1(0)-eps2g(1)*k1(1)-eps2g(2)*k1(2)-eps2g(3)*k1(3)
       eps2gk2 = eps2g(0)*k2(0)-eps2g(1)*k2(1)-eps2g(2)*k2(2)-eps2g(3)*k2(3)
       eps2gq1 = eps2g(0)*q1(0)-eps2g(1)*q1(1)-eps2g(2)*q1(2)-eps2g(3)*q1(3)

       eps1gq1 = eps1g(0)*q1(0)-eps1g(1)*q1(1)-eps1g(2)*q1(2)-eps1g(3)*q1(3)
       eps2gq2 = eps2g(0)*q2(0)-eps2g(1)*q2(1)-eps2g(2)*q2(2)-eps2g(3)*q2(3)
c
c
       ce1g = -2*(d12+d24)*eps2gk2*(q1sqc+q2sqc-3*s-4*t) - 
     1      2*(d12-d24)*eps2gq1*(q2sqc-t) 
     2  +4*(- d11*eps2gk2*s - d12*eps2gk1*t + d13*eps2gk2*
     3  (q2sqc-s-t) + d13*eps2gq1*q2sqc - d21*eps2gk2*s 
     4 - d22*eps2gk2*t - d22*eps2gq1*q2sqc + d23*eps2gk2*q1sqc
     5 + d25*eps2gk2*(q2sqc-s-2*t) 
     6 - d26*eps2gk2*(q2sqc-s-t)+d26*eps2gq1*t+2*d27*eps2gq1
     7 - d32*eps2gk2*q2sqc - d34*eps2gk2*(q2sqc-t) 
     8 + d35*eps2gk2*(q2sqc-s-t)+d36*eps2gk2*(2*q2sqc-t) 
     9 + d37*eps2gk2*q1sqc + d38*eps2gk2*(q1sqc+q2sqc-s)
     & - d39*eps2gk2*q1sqc 
     1 - d310*eps2gk2*(q1sqc+2*q2sqc-2*s-t)
     2 - 4*d311*eps2gk2 + 6*d312*eps2gk2 + 2*d313*eps2gq1)

       ce2g= -2*(d12+d24)*(eps1gk2*(q1sqc+q2sqc-s-2*t)
     1  +eps1gq2*(q2sqc-s-3*t)) + 
     2  4*(d13*eps1gk2*(q1sqc-2*s-t) - 
     3  d13*eps1gq2*(2*s+t) + d22*eps1gk2*q2sqc
     4  + d22*eps1gq2*q2sqc - d23*eps1gk2*t + 
     5    d23*eps1gq2*(q1sqc-t) - d24*eps1gk1*q2sqc 
     6  + d25*eps1gk2*q1sqc + d25*eps1gk1*(q2sqc-2*s-t) 
     7  + d26*eps1gk2*t - d26*eps1gk1*(q1sqc-s) - 
     8  2*d27*eps1gq2 + d33*eps1gk2*q1sqc + d33*eps1gq2*q1sqc
     9  + d37*eps1gk1*(q2sqc-s-t) + d38*eps1gk1*q2sqc
     &  - d39*eps1gk1*(q1sqc+q2sqc-s) - d310*eps1gk1*(q2sqc-t)
     1  +2*d311*eps1gk2+2*d312*eps1gq2-6*d313*eps1gk1)

C NOTE: cq = 1/2 * cq1 in my reduce code
       cqg = (d12+d24)*eps1geps2g*s + 2*(4*d12*eps2gk2*eps1gk2 + 
     1 3*d12*eps2gk2*eps1gq2+d12*eps2gq1*eps1gk2-4*d13*eps2gk2*eps1gk2-
     2 2*d13*eps2gk2*eps1gq2-2*d13*eps2gq1*eps1gk2-d13*eps1geps2g*s+
     3 2*d22*eps2gk2*eps1gq2-d22*eps1geps2g*t-2*d23*eps2gq1*eps1gk2-
     4 2*d23*eps2gq1*eps1gq2-d23*eps1geps2g*t+6*d24*eps2gk2*eps1gk2+
     5 3*d24*eps2gk2*eps1gq2+d24*eps2gq1*eps1gk2-6*d25*eps2gk2*eps1gk2-
     6 2*d25*eps2gk2*eps1gq2-2*d25*eps2gq1*eps1gk2-d25*eps1geps2g*s-
     7 4*d26*eps2gk2*eps1gq2+4*d26*eps2gq1*eps1gk2+2*d26*eps2gq1*eps1gq2+ 
     8 d26*eps1geps2g*(s+2*t)-d32*eps1geps2g*q2sqc+
     9 d33*eps1geps2g*q1sqc+2*d34*eps2gk2*eps1gk2-2*d35*eps2gk2*eps1gk2+
     & 2*d36*eps2gk2*eps1gq2+d36*eps1geps2g*(q2sqc-t)-
     1 2*d37*eps2gq1*eps1gk2+d37*eps1geps2g*(q2sqc-s-t)+
     2 2*d38*eps2gq1*eps1gq2+d38*eps1geps2g*(q1sqc+2*q2sqc-s) -
     3 2*d39*eps2gq1*eps1gq2-d39*eps1geps2g*(2*q1sqc+q2sqc-s) -
     4 2*d310*eps2gk2*eps1gq2+2*d310*eps2gq1*eps1gk2-d310*eps1geps2g*
     5 (2*q2sqc-s-2*t)+4*d312*eps1geps2g-4*d313*eps1geps2g)
      

       gce1g = -2*eps2gq2*( d12*(2*q1sqc+q2sqc-s-3*t) - 
     1       2*(d13-d22)*(q1sqc+q2sqc-s-t)
     2     - d24*(q2sqc-s-t)-2*d26*q1sqc-4*d27-4*d312)

       gce2g = 2*eps1gq1*(d12*(q2sqc+t)-2*d13*(q2sqc+s) + 
     1      2*d22*q2sqc + 2*d23*(q1sqc+q2sqc-s-2*t)
     2     +d24*(q2sqc+t)+2*d25*(-q2sqc+s+t)-2*d26*(q2sqc-t)
     3     -4*d27+2*d33*q1sqc + 4*d313)

       gcqg = 4*eps1gq1*eps2gq2*(d12-d13+d22-d23+d38-d39) +
     1        2*eps1gq1*eps2gk2*(d12-4*d23+d24+2*d26-2*d37+2*d310) -
     2        4*eps1gq1*eps2gq1*(d23-d26+d33-d39) +
     3        2*eps2gq2*eps1gk2*(3*d12-4*d13+4*d22+d24-2*d25-
     4                                     2*d26+2*d36-2*d310) +
     5        4*eps2gq2*eps1gq2*(d12-d13+2*d22-2*d26+d32-d38)
     
            ce1g = ce1g + gce1g
            ce2g = ce2g + gce2g
            cqg  = cqg + gcqg
c
c  now add the terms from the vertex and propagator corrections

            ce1g = ce1g + 2*eps2gk2*Teps2c  + 2*eps2gq2*Tg2c
            ce2g = ce2g + 2*eps1gk1*Teps1c  - 2*eps1gq1*Tg1c
c  and the final result
            mvirtg = me1g*ce1g + me2g*ce2g + mqg*cqg + mborng*cborn
	    
	endif !gimprove
	
c -------------------------------
	    
       if (gimprove) then
       bcount = bcount+1

	icountmax = 2**(N_points(0)+2)
	if(icountmax.eq.1) icountmax = 10
	
        if (abs(mvirtg)/abs(mborn).ge.1e-3) then		
          if ((abs((mvirtg))/(abs(me1g*ce1g)+abs(mqg*cqg)+
     #  		     abs(me2g*ce2g)+abs(mborng*cborn)).ge.1e-3)
     #		.or.(abs(mvirtg).ge.1e3))	 then
        
              gcount = gcount+1 		      
              mvirt = 0d0  !disregard this contribution
          
          endif       
        endif

        if (ldebug) then
           if (mod(bcount,icountmax).eq.0 ) then
              print*,"counting ratio:",dble(gcount)/dble(bcount)
           endif
	
           if(breset) then
	      if (bcount .eq. 1000000000 ) then!integer range exceeded -> reset
	         print*,"gcount,bcount = ", gcount,bcount," are reset now"
	         bcount = 0
	         gcount = 0
	      endif
           endif	  
        endif


       endif !gimprove

c ---------------------------------------------------------	   
	
            jvirt(mu,isig) = mvirt
         enddo !isig
      enddo !mu
            
      return
      end


c---------------------------------------------------------------------

      subroutine  J_virtual_tri_box_pent1(psi,p,nq1,nq2,sign, musq,
c     input currents (any ORDERED couple among [J1,J2,J3]) in DHELAS format
c     They can be  jmu1,jmu2  or  jmu1,jmu3   or   jmu2,jmu3
     #     jal,jbe,
c     uncontracted external current. It can be 1, 2 or 3
     #     num_ret_current,
c     the corresponding Born current
c           Jborn,
c     returned current
     #     Jret)

      implicit none
      complex * 16 psi(2,-1:1,4)
      real * 8 p(0:3,4),q(0:3),l(0:3),t(0:3)
      real * 8 musq
      complex * 16 jal(0:5),jbe(0:5)
      integer num_ret_current, nq1,nq2,nq3,nq4, sign(4)
      complex * 16 jret(0:3)
      logical lpr
      common /e0print/lpr

      integer mu
c map input momenta onto incoming standard and make use of 
c DHELAS information in currents
      if (nq1.eq.1 .and. nq2.eq.2) then
         nq3 = 3
         nq4 = 4
      elseif (nq1.eq.3 .and. nq2.eq.4) then
         nq3 = 1
         nq4 = 2
      else
         print*," check nq1 etc. in J_virtual_tri_box_pent1 call "
         print*,nq1,nq2,nq3,nq4
         stop
      endif
         
      if (num_ret_current.eq.1) then     ! p3=q+p4 for q incoming on 12 line
         do mu = 0,3
            q(mu) = p(mu,nq3)-p(mu,nq4)
         enddo
         l(0) = -dreal(jal(4))
         l(1) = -dreal(jal(5))
         l(2) = -dimag(jal(5))
         l(3) = -dimag(jal(4))
         t(0) = -dreal(jbe(4))
         t(1) = -dreal(jbe(5))
         t(2) = -dimag(jbe(5))
         t(3) = -dimag(jbe(4))
      elseif (num_ret_current.eq.2) then
         do mu = 0,3
            l(mu) = p(mu,nq3)-p(mu,nq4)
         enddo
         q(0) = -dreal(jal(4))
         q(1) = -dreal(jal(5))
         q(2) = -dimag(jal(5))
         q(3) = -dimag(jal(4))
         t(0) = -dreal(jbe(4))
         t(1) = -dreal(jbe(5))
         t(2) = -dimag(jbe(5))
         t(3) = -dimag(jbe(4))
      elseif (num_ret_current.eq.3) then
         do mu = 0,3
            t(mu) = p(mu,nq3)-p(mu,nq4)
         enddo
         q(0) = -dreal(jal(4))
         q(1) = -dreal(jal(5))
         q(2) = -dimag(jal(5))
         q(3) = -dimag(jal(4))
         l(0) = -dreal(jbe(4))
         l(1) = -dreal(jbe(5))
         l(2) = -dimag(jbe(5))
         l(3) = -dimag(jbe(4))
      endif

      do mu = 0,3
         jret(mu) = 0
      enddo
      call  J_virtual_tri_box_pent (
c     bra(p"), ket(p) and helicity 
     #     psi(1,-1,nq2),psi(1,-1,nq1),-1,
c     momenta and squared renormalization scale
     #     p(0,nq1),q,l,t,p(0,nq2),musq,
c     input currents (any ORDERED couple among [J1,J2,J3])
c     They can be  jmu1,jmu2  or  jmu1,jmu3   or   jmu2,jmu3
     #     jal,jbe,
c     uncontracted external current. It can be 1, 2 or 3
     #     num_ret_current,
c     returned current
     #     Jret)

c check momentum conservation and subtract Born contribution
      do mu = 0,3
         t(mu) = q(mu)+l(mu)+t(mu)+p(mu,nq1)-p(mu,nq2) 
         if (abs(t(mu)).gt.1d-8) then
            print*," no momentum conservation in j_virtual_tri_box_pent"
            print*," mom diff(",mu,") = ",t(mu)
         endif
      enddo
      return
      end

#endif


      subroutine jtest(jin,jout)
      implicit none
      double complex jin(0:5), jout(0:5), x, jq
      real*8 q(0:3), qq, mi, mo
      integer mu

      q(0) = dreal(jin(4))
      q(1) = dreal(jin(5))
      q(2) = dimag(jin(5))
      q(3) = dimag(jin(4))
      jq = jin(0)*q(0)+ jin(1)*q(1)+ jin(2)*q(2)+ jin(3)*q(3)
      qq = q(0)*q(0)+ q(1)*q(1)+ q(2)*q(2)+ q(3)*q(3)
      X = jq/qq
      mi = 0
      mo = 0
      do mu = 0,3
         jout(mu) = jin(mu) - x*q(mu)
         mi = max(mi,abs(jin(mu)))
         mo = max(mo,abs(jout(mu)))
         print*," mu = ",mu,": ",jout(mu),abs(jout(mu)/jin(mu))
      enddo
      print*," ratio = ",mo/mi
      jout(4) = jin(4)
      jout(5) = jin(5)

      end
