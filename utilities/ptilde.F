*************************************************************************
      SUBROUTINE ptilde(p, xuz, v)
c*************************************************************************
c     This is the manager of all ptilde subroutines and the (only) 
c     interface of those to the main program.  Depending on the process ID
c     the right ptilde routine is called.
c
c*************************************************************************
c     INPUT
c     p(mu,1:np,1) : np parton momenta for real emission graph
c     v(mu,1:nv,1) : nv momenta of other particles for real emission graph
c     
c     OUTPUT
c     p(mu,1:np,j) : np parton momenta in the ptilde kinematic
c     v(mu,1:nv,j) : nv momenta of other particles in the ptilde kinematic
c                    (j = 2, ..., n_kin)
c     xuz          : These are the z, u, x_ij,a of the Catani-Seymour algorithm
c*************************************************************************

      implicit none

#include "global.inc"
#include "process.inc"

      real*8 p(0:3,max_p,max_kin)
      real*8 v(0:3,max_v,max_kin)
      real*8 xuz(2,2:max_kin)

      
      SELECT CASE(procID)

c VBF: Hjj, Vjj, VVjj
      CASE(Hjj, Hjj_AA, Hjj_mu, Hjj_tau, Hjj_bbar, Hjj_WW, Hjj_ZZ_ll,
     1     Hjj_ZZ_lnu, HHjj, HHjj_bbtau, HHjj_bbAA, Ajj, Zjj_l, Zjj_nu, WPjj, WMjj, WPWMjj, ZZjj_ll,
     2     ZZjj_lnu, WPZjj, WMZjj, WPWPjj, WMWMjj, WPAjj, WMAjj,ZAJJ,ZAjj_n,
     3     HAjj, HAjj_AA, HAjj_mu, HAjj_tau, 
     4     HAjj_bbar, HAjj_WW, HAjj_ZZ_ll, HAjj_ZZ_lnu, AAjj,
     5     WPhadWMjj, WPWMhadjj, WPhadZjj, WPZhadjj, WMhadZjj, WMZhadjj, ZZhadjj,
     6     Hjj_WPhadWM, Hjj_WPWMhad, Hjj_ZZhad, WPhadWPjj, WMhadWMjj, 
     7     Sp2jj_WW, Sp2jj_ZZ_ll, Sp2jj_ZZ_lnu)
          call VBFptildemake(3,p(0,1,1),xuz,v(0,1,1))

c VBF: Hjjj
      CASE(Hjjj,Hjjj_AA, Hjjj_mu, Hjjj_tau, Hjjj_bbar, Hjjj_WW, 
     1         Hjjj_ZZ_ll,Hjjj_ZZ_lnu)
          call HJJJptildemake(n_kin,p(0,1,1),xuz(1,2),v(0,1,1))

c DIBOSON, TRIBOSON
      CASE( WW, WWZ, ZZZ, ZZWP, ZZWM, WWWP, WWWM , WWA, ZZA, ZZnA, WPZA, WMZA,
     1      WPAA, WMAA, ZAA, ZnAA, AAA, WPA, WMA, WPZ, WMZ,
     2      ZZ, ZA, AA,
     3      WPhadWMlep,WPlepWMhad,ZZhad,WPhadZ,WPZhad,WMhadZ,WMZhad,
     4      WPhadWMZ, WPWMhadZ, WWZhad, ZZhadWP, ZZWPhad, ZZhadWM, ZZWMhad, 
     5      WPhadWMWP, WPWMhadWP, WMhadWPWM, WMWPhadWM, ZZZhad,
     6      WPhadWMA, WPWMhadA, ZZhadA, WPhadZA, WPZhadA, WMhadZA, WMZhadA,
     6      WP_only, WM_only,
     7      WPH, WPH_AA, WPH_mu, WPH_tau, WPH_bbar, WPH_WW, WPH_ZZ_ll,  WPH_ZZ_lnu,
     7      WMH, WMH_AA, WMH_mu, WMH_tau, WMH_bbar, WMH_WW, WMH_ZZ_ll,  WMH_ZZ_lnu
     8      )

          call DIBOSptildemake(3,p(0,1,1),xuz(1,2),v(0,1,1))


c TRIBOSON plus jet
      CASE( WPAAj, WMAAj)
          call TriJetptildemake(3,p(0,1,1),xuz(1,2),v(0,1,1),4)


c DIBSON PLUS JET
      CASE( WMZJ, WPZJ, WMAJ, WPAJ, 
     -      WPhadZJ, WPZhadJ, WMhadZJ, WMZhadJ,
     7      WPHJ, WPHJ_AA, WPHJ_mu, WPHJ_tau, WPHJ_bbar, WPHJ_WW, WPHJ_ZZ_ll,  WPHJ_ZZ_lnu,
     7      WMHJ, WMHJ_AA, WMHJ_mu, WMHJ_tau, WMHJ_bbar, WMHJ_WW, WMHJ_ZZ_ll,  WMHJ_ZZ_lnu,
     -      WPJ, WMJ)
          call WGAJptildemake(p(0,1,1),xuz(1,2),v(0,1,1))

c qcdVVjj
      CASE( QCDWPZjj, QCDWMZjj, QCDWPWPjj,QCDWMWMjj,QCDWPjj,QCDWMjj,QCDZZjj_ll,QCDZZjj_lnu,QCDZjj_l,QCDZjj_nu 
     &     ,QCDWPAjj, QCDWMAjj,QCDZAjj_l,QCDZAjj_n,QCDAAjj)

	  call qcdVVJJptildemake(p(0,1,1),xuz(1,2),v(0,1,1))
      

c DIBOSON PLUS JET -- neutral ones
      CASE(WWJ, WPHADWMJ, WPWMHADJ, ZZJ)
          call VVJptildemake(p(0,1,1),xuz(1,2),v(0,1,1))

      CASE DEFAULT
         print *,"Error: ptilde routine for process "//
     &            "not implemented, ID = ",procID
         STOP
      END SELECT

      return
      end

c************************************************************************
      subroutine VBFptildemake(N,p,xuz,v)
c*************************************************************************
c
c assign ptilde momenta in notation of Catani&Seymour
c
c  input:  p(mu,1:np,1)      np parton momenta for real emission graph
c  output:                   np-1 parton momenta for Born amplitude
c          p(mu,1:np-1,2)        index=2: emission off 1-2 line
c          p(mu,1:np-1,3)        index=3: emission off 3-4 line
c
c          xuz(1,index)      x_ija etc
c          xuz(2,index)      u or z
c
c at the same time, fill the legoplot variables for the external momenta
c
c       Dieter Zeppenfeld, <dieter@particle.uni-karlsruhe.de>
c       Initial version:  November 2002
c       Last modified: September, 2008
c*************************************************************************

      implicit none

#include "global.inc"

      integer mu, a, i, j, N, ii, jj, k
      real*8 p(0:3,max_p,max_kin), xuz(2,2:3), v(0:3,max_v,max_kin)
      real*8 omx, pipj, pipa, pjpa
      real*8 tpi
      parameter (tpi=2d0*3.141592653589793d0)
      logical ldebug
      parameter (ldebug=.false.)

c  parton np is the potentially soft/collinear one
      do ii = 2,N
         if (ii.eq.2) then        ! set mom id for radiation off 1,2 line
            a = 1
            i = 5
            j = 3
         elseif (ii.eq.3) then
            a = 2
            i = 5
            j = 4
         else
            stop
         endif
         pipj = p(0,i,1)*p(0,j,1)-p(1,i,1)*p(1,j,1)-
     1          p(2,i,1)*p(2,j,1)-p(3,i,1)*p(3,j,1)
         pipa = p(0,i,1)*p(0,a,1)-p(1,i,1)*p(1,a,1)-
     1          p(2,i,1)*p(2,a,1)-p(3,i,1)*p(3,a,1)
         pjpa = p(0,j,1)*p(0,a,1)-p(1,j,1)*p(1,a,1)-
     1          p(2,j,1)*p(2,a,1)-p(3,j,1)*p(3,a,1)
         omx = pipj/(pipa+pjpa)
         xuz(1,ii) = 1d0 - omx
         xuz(2,ii) = pipa/(pipa+pjpa)
         if (ldebug) print*," assign momenta # ",a,j,3-a,7-j
         do mu = 0,3
            p(mu,5,ii) = 0
            p(mu,a,ii) = xuz(1,ii)*p(mu,a,1)
            p(mu,j,ii) = p(mu,j,1)+p(mu,i,1)-omx*p(mu,a,1)
            p(mu,3-a,ii)=p(mu,3-a,1)
            p(mu,7-j,ii)=p(mu,7-j,1)
         enddo

         if (ldebug) then
check momentum conservation
            do mu = 0,3
               omx = p(mu,n_p,1)
               do jj = 1,n_p-1
                  if(jj.le.2) omx = omx-p(mu,jj,1)+p(mu,jj,ii)
                  if(jj.gt.2) omx = omx+p(mu,jj,1)-p(mu,jj,ii)
               enddo
               print*," mom.cons. for mu=",mu,": sum= ",omx
            enddo
check masslessness of ptilde_a and ptilde_j
            print*," p_a^2 = ",p(0,a,ii)**2-p(1,a,ii)**2-
     1                         p(2,a,ii)**2-p(3,a,ii)**2
            print*," p_j^2 = ",p(0,j,ii)**2-p(1,j,ii)**2-
     1                         p(2,j,ii)**2-p(3,j,ii)**2
         endif
      enddo

      do i = 1, n_v
         do k = 1, n_kin
            do mu = 0,3
               v(mu,i,k) = v(mu,i,1)
            enddo
         enddo
      enddo

      if ( ldebug ) then
 10      format( " p(", i1, ") = ", 4(f10.3, 2x) )
c
         do j = 1,N
            if (j.eq.1) print*," (n+1) parton momenta "
            if (j.eq.2) print*," n-parton momenta for upper NLO corr."
            if (j.eq.3) print*," n-parton momenta for lower NLO corr."
            do i = 1, n_p-min(j-1,1)
               write(6,10) i, p(0,i,j), p(1,i,j), p(2,i,j), p(3,i,j)
            end do
         enddo
      end if
      end


c********************** subroutine wgajptildemake *************************
c
c    this subroutine assigns ptilde momenta to the partons and the final
c    state leptons in notation of Catani and Seymour.
c
c INPUT    p(mu,1:np,1)      np parton momenta for real emission graph
c          v(mu,1:nv,1)      nv lepton momenta for real emission graph
c
c OUTPUT   np-1 parton momenta for Born amplitude
c          p(mu,1:np-1,2...)    tilde momenta

c          lepton momenta for born amplitude
c          v(mu,1:nv,2.....)    tilde momenta

c
c          xuz               x_ija in the subtraction terms
c
c**************************************************************************
      subroutine wgajptildemake(p,xuz,v)

      implicit none

#include "global.inc"

      integer jj,ii,mu

      double precision pin(0:3),pk(0:3),pa(0:3),pj(0:3),v(0:3,max_v,max_kin),
     &    p(0:3,max_p,max_kin),pb(0:3),p4(0:3),p3(0:3),ptildea(0:3),ptildeij(0:3),
     &    pipb,papb,k(0:3),Kt(0:3),pipj,xija,zi,papj,papi,
     &    vstore(0:3,max_v)
c
      logical ldebug
      parameter (ldebug=.false.)

      double precision ptildeai(0:3),ptildek(0:3),
     &                 pipk,pkpa,pipa
      double precision xika,ui,xuz(2,2:max_kin),xiab,one
      parameter (one=1d0)


      double precision dotrr
      external dotrr

      if (ldebug) then
         print*, 'debug in ptilde.F'
c store elw parameters
         do mu=0,3
            do jj=1,n_v
               vstore(mu,jj)=v(mu,jj,1)
            enddo
         enddo
      endif

c  set the unused ptilde-momenta to zero
      do ii = 2,max_kin
         do jj=1,max_p
            do mu = 0,3
               p(mu,jj,ii) = 0d0
            enddo
         enddo
      enddo
      do ii = 2,max_kin
         do jj=1,n_v
            do mu = 0,3
               v(mu,jj,ii) = 0d0
            enddo
         enddo
      enddo


c output real configuration      
      if (ldebug) then
         print*, "real emission configuration"
         call printMomentumConservation(p,v,1)
      endif

c  store 3rd and 4th momenta
         do mu=0,3
         p4(mu)=p(mu,4,1)
         p3(mu)=p(mu,3,1)
         enddo

c

c   get ptilde kinematics for u~u~ > e- ve~ a u~ d~
c   this involves kinematics of Catani Seymour 5.67,5.147

c Dipoles 5.67

      do mu=0,3
      pa(mu) =p(mu,1,1)
      pin(mu)=p(mu,4,1)
      pk(mu) =p(mu,3,1)
      enddo

      pipk=dotrr(pin,pk)
      pkpa=dotrr(pk,pa)
      pipa=dotrr(pin,pa)

      xika=one - pipk/(pkpa+pipa)! (pkpa+pipa-pipk)/(pkpa+pipa)

      do mu=0,3
      ptildeai(mu)=xika*pa(mu)
      ptildek(mu) =pk(mu)+pin(mu)-(one-xika)*pa(mu)
      enddo
      ui=pipa/(pipa+pkpa)


c   generate 1. output
c
      do mu=0,3
      p(mu,1,2)=ptildeai(mu)       !this is the intermediate gluon, dipole (13,4)
      p(mu,2,2)=p(mu,2,1)
      p(mu,3,2)=ptildek(mu)
      do jj=1,n_v
         v(mu,jj,2)=v(mu,jj,1)
      enddo
      enddo

      xuz(1,2)=xika
      xuz(2,2)=ui

      if (ldebug) then
         print*, "New Dipole 1"
         call printMomentumConservation(p,v,2)
      endif

c
c start new dipole
c
      do mu=0,3
      pa(mu)=p(mu,2,1)
      pin(mu)=p(mu,4,1)
      pk(mu)=p(mu,3,1)
      enddo

      pipk=dotrr(pin,pk)
      pkpa=dotrr(pk,pa)
      pipa=dotrr(pin,pa)

      xika=one - pipk/(pkpa+pipa)!(pkpa+pipa-pipk)/(pkpa+pipa)

      do mu=0,3
      ptildeai(mu)=xika*pa(mu)
      ptildek(mu) =pk(mu)+pin(mu)-(one-xika)*pa(mu)
      enddo
      ui=pipa/(pipa+pkpa)

c   generate 2. output
      do mu=0,3
      p(mu,1,3)=ptildeai(mu)       !this is the intermediate gluon, dipole (23,4)
      p(mu,2,3)=p(mu,1,1)
      p(mu,3,3)=ptildek(mu)
      do jj=1,n_v
         v(mu,jj,3)=v(mu,jj,1)
      enddo
      enddo
      xuz(1,3)=xika
      xuz(2,3)=ui

c    test output
      if (ldebug) then
         print*, "New Dipole 2"
         call printMomentumConservation(p,v,3)
      endif


c Dipoles 5.147
c
c   start new dipole
c
      do mu=0,3
      pa(mu)=p(mu,1,1)
      pin(mu)=p(mu,4,1)
      pb(mu)=p(mu,2,1)
      enddo

      pipb=dotrr(pin,pb)
      papb=dotrr(pa,pb)
      pipa=dotrr(pin,pa)
      xiab=one- (pipa+pipb)/papb !(papb-pipa-pipb)/papb

      do mu=0,3
      ptildeai(mu)=xiab*pa(mu)
      K(mu)=pa(mu)+pb(mu)-pin(mu)
      Kt(mu)=ptildeai(mu)+pb(mu)
      enddo

      Call boostkinematics(k,kt,p(0,1,1),v(0,1,1),p(0,1,4),v(0,1,4))

c   generate 3. output
      do mu=0,3
      p(mu,1,4)=ptildeai(mu)       !this is the intermediate gluon, dipole (13,2)
      p(mu,2,4)=p(mu,2,1)
      enddo
      xuz(1,4)=xiab
      xuz(2,4)=0d0 !!!!!!!!!!!!!


c    test dipole
c
      if (ldebug) then
         print*, "New Dipole 3"
         call printMomentumConservation(p,v,4)
      endif

c   start new dipole
c
      do mu=0,3
      pa(mu)=p(mu,2,1)
      pin(mu)=p(mu,4,1)
      pb(mu)=p(mu,1,1)
      enddo

      pipb=dotrr(pin,pb)
      papb=dotrr(pa,pb)
      pipa=dotrr(pin,pa)
      xiab= one- (pipa+pipb)/papb!(papb-pipa-pipb)/(papb)

      do mu=0,3
      ptildeai(mu)=xiab*pa(mu)
      K(mu)=pa(mu)+pb(mu)-pin(mu)
      Kt(mu)=ptildeai(mu)+pb(mu)
      enddo

      Call boostkinematics(k,kt,p(0,1,1),v(0,1,1),p(0,1,5),v(0,1,5))

c   generate 4. output
      do mu=0,3
      p(mu,1,5)=ptildeai(mu)       !this is the intermediate gluon, dipole (23,1)
      p(mu,2,5)=p(mu,1,1)
      enddo
      xuz(1,5)=xiab
      xuz(2,5)=0d0 !!!!!!!!!!!!!!!!

c    test dipole
c
      if (ldebug) then
         print*, "New Dipole 4"
         call printMomentumConservation(p,v,5)
      endif


c   get ptilde kinematics for u d > e- ve~ a u u and inicial state crossed
c   this involves kinematics of Catani Seymour 5.67,5.147

c Dipoles 5.67

      do mu=0,3
      pa(mu) =p(mu,1,1)
      pin(mu)=p(mu,3,1)
      pk(mu) =p(mu,4,1)
      enddo

      pipk=dotrr(pin,pk)
      pkpa=dotrr(pk,pa)
      pipa=dotrr(pin,pa)

      xika=one - pipk/(pkpa+pipa)! (pkpa+pipa-pipk)/(pkpa+pipa)

      do mu=0,3
      ptildeai(mu)=xika*pa(mu)
      ptildek(mu) =pk(mu)+pin(mu)-(one-xika)*pa(mu)
      enddo
      ui=pipa/(pipa+pkpa)

c   generate 5. output
c
      do mu=0,3
      p(mu,1,6)=ptildeai(mu)       !this is the intermediate gluon, dipole (13,4)
      p(mu,2,6)=p(mu,2,1)
      p(mu,3,6)=ptildek(mu)
      do jj=1,n_v
         v(mu,jj,6)=v(mu,jj,1)
      enddo
      enddo
      xuz(1,6)=xika
      xuz(2,6)=ui

      if (ldebug) then
         print*, "New Dipole 5"
         call printMomentumConservation(p,v,6)
      endif

c
c start new dipole
c
      do mu=0,3
      pa(mu)=p(mu,2,1)
      pin(mu)=p(mu,3,1)
      pk(mu)=p(mu,4,1)
      enddo

      pipk=dotrr(pin,pk)
      pkpa=dotrr(pk,pa)
      pipa=dotrr(pin,pa)

      xika=one - pipk/(pkpa+pipa)! (pkpa+pipa-pipk)/(pkpa+pipa)

      do mu=0,3
      ptildeai(mu)=xika*pa(mu)
      ptildek(mu) =pk(mu)+pin(mu)-(one-xika)*pa(mu)
      enddo
      ui=pipa/(pipa+pkpa)

c   generate 6. output
      do mu=0,3
      p(mu,1,7)=ptildeai(mu)       !this is the intermediate gluon, dipole (23,4)
      p(mu,2,7)=p(mu,1,1)
      p(mu,3,7)=ptildek(mu)
      do jj=1,n_v
         v(mu,jj,7)=v(mu,jj,1)
      enddo
      enddo
      xuz(1,7)=xika
      xuz(2,7)=ui

c    test output
      if (ldebug) then
         print*, "New Dipole 6"
         call printMomentumConservation(p,v,7)
      endif


c Dipoles 5.147
c
c   start new dipole
c
      do mu=0,3
      pa(mu)=p(mu,1,1)
      pin(mu)=p(mu,3,1)
      pb(mu)=p(mu,2,1)
      p(mu,3,1)=p4(mu)
      enddo

      pipb=dotrr(pin,pb)
      papb=dotrr(pa,pb)
      pipa=dotrr(pin,pa)
      xiab= one- (pipa+pipb)/papb!(papb-pipa-pipb)/papb

      do mu=0,3
      ptildeai(mu)=xiab*pa(mu)
      K(mu)=pa(mu)+pb(mu)-pin(mu)
      Kt(mu)=ptildeai(mu)+pb(mu)
      enddo

      Call boostkinematics(k,kt,p(0,1,1),v(0,1,1),p(0,1,8),v(0,1,8))

      do mu=0,3
      p(mu,3,1)=p3(mu)
      enddo

c   generate 7. output
      do mu=0,3
      p(mu,1,8)=ptildeai(mu)       !this is the intermediate gluon, dipole (13,2)
      p(mu,2,8)=p(mu,2,1)
      enddo
      xuz(1,8)=xiab
      xuz(2,8)=0d0 !!!!!!!!!!!!!

c    test dipole
c
      if (ldebug) then
         print*, "New Dipole 7"
         call printMomentumConservation(p,v,8)
      endif

c   start new dipole
c
      do mu=0,3
      pa(mu)=p(mu,2,1)
      pin(mu)=p(mu,3,1)
      pb(mu)=p(mu,1,1)
      p(mu,3,1)=p4(mu)
      enddo

      pipb=dotrr(pin,pb)
      papb=dotrr(pa,pb)
      pipa=dotrr(pin,pa)
      xiab= one- (pipa+pipb)/papb!(papb-pipa-pipb)/(papb)

      do mu=0,3
      ptildeai(mu)=xiab*pa(mu)
      K(mu)=pa(mu)+pb(mu)-pin(mu)
      Kt(mu)=ptildeai(mu)+pb(mu)
      enddo

      Call boostkinematics(k,kt,p(0,1,1),v(0,1,1),p(0,1,9),v(0,1,9))

      do mu=0,3
      p(mu,3,1)=p3(mu)
      enddo

c   generate 8. output
      do mu=0,3
      p(mu,1,9)=ptildeai(mu)       !this is the intermediate gluon, dipole (23,1)
      p(mu,2,9)=p(mu,1,1)
      enddo
      xuz(1,9)=xiab
      xuz(2,9)=0d0 !!!!!!!!!!!!!!!!

c    test dipole
c
      if (ldebug) then
         print*, "New Dipole 8"
         call printMomentumConservation(p,v,9)
      endif

c Dipoles 5.40

      do mu=0,3
      pa(mu) =p(mu,1,1)
      pin(mu)=p(mu,3,1)
      pj(mu) =p(mu,4,1)
      enddo

      papj=dotrr(pa,pj)
      papi=dotrr(pa,pin)
      pipj=dotrr(pin,pj)

      xija=one - pipj/(papi+papj)! (papi+papj-pipj)/(papi+papj)

      do mu=0,3
      ptildea(mu)=xija*pa(mu)
      ptildeij(mu) =pin(mu)+pj(mu)-(one-xija)*pa(mu)
      enddo
      zi=papi/(papi+papj)

c   generate 9. output
c
      do mu=0,3
      p(mu,1,10)=ptildea(mu)       !this is the intermediate gluon, dipole (13,4)
      p(mu,2,10)=p(mu,2,1)
      p(mu,3,10)=ptildeij(mu)
      do jj=1,n_v
         v(mu,jj,10)=v(mu,jj,1)
      enddo
      enddo
      xuz(1,10)=xija
      xuz(2,10)=zi

      if (ldebug) then
         print*, "New Dipole 9"
         call printMomentumConservation(p,v,10)
      endif

      do mu=0,3
      pa(mu) =p(mu,2,1)
      pin(mu)=p(mu,3,1)
      pj(mu) =p(mu,4,1)
      enddo

      papj=dotrr(pa,pj)
      papi=dotrr(pa,pin)
      pipj=dotrr(pin,pj)

      xija=one - pipj/(papi+papj)! (papi+papj-pipj)/(papi+papj)

      do mu=0,3
      ptildea(mu)=xija*pa(mu)
      ptildeij(mu) =pin(mu)+pj(mu)-(one-xija)*pa(mu)
      enddo
      zi=papi/(papi+papj)

c   generate 10. output
c
      do mu=0,3
      p(mu,1,11)=p(mu,1,1)       !this is the intermediate gluon, dipole (13,4)
      p(mu,2,11)=ptildea(mu)
      p(mu,3,11)=ptildeij(mu)
      do jj=1,n_v
         v(mu,jj,11)=v(mu,jj,1)
      enddo
      enddo
      xuz(1,11)=xija
      xuz(2,11)=zi

      if (ldebug) then
         print*, "New Dipole 10"
         call printMomentumConservation(p,v,11)
         print*, "Elw input momenta invariant?"
         do ii=1,n_v
            print*, vstore(0:3,ii)-v(0:3,ii,1)
         enddo
         !read *
      endif


      end

      subroutine printMomentumConservation(p,v,i)
#include "global.inc"
         double precision p(0:3,max_p,max_kin), v(0:3,max_v,max_kin)
         double precision debugmom(0:3)
         integer i !which kinematic to print
         integer ii

         !print momenta for full debug
         print*, "Momenta for configuration ", i
         print*, 'n_p', n_p
         do ii=1,n_p
            print*, p(0:3,ii,i)
         enddo
         print*, 'n_v', n_v
         do ii=1,n_v
            print*, v(0:3,ii,i)
         enddo

         print*, "Momentum conservation of real em. kinematics (should be zero):"
         ! TODO: Fix output for higgs decay
         debugmom = 0d0
!         print*, p(0:3,1:4,1)
!         print*, v(0:3,1:n_v,1)
!         print*, debugmom
         do ii=1,2
            debugmom = debugmom + p(0:3,ii,1)
         enddo
         do ii=3,4
            debugmom = debugmom - p(0:3,ii,1)
         enddo
!         print*, debugmom
         do ii=1,n_v
            debugmom = debugmom - v(0:3,ii,1)
         enddo
         print*, debugmom
         print*, "Momentum conservation of Dipole kinematics (should be zero):"
         debugmom = 0d0
         do ii=1,2
            debugmom = debugmom + p(0:3,ii,i)
         enddo
            debugmom = debugmom - p(0:3,3,i)
         do ii=1,n_v
            debugmom = debugmom - v(0:3,ii,i)
         enddo
         print*, debugmom
         print*, "------------------------------------------------"
      end subroutine

c********************** subroutine qcdVVjjptildemake **********************
c
c    this subroutine assigns ptilde momenta to the partons and the final
c    state leptons in notation of Catani and Seymour.
c
c INPUT    p(mu,1:np,1)      np parton momenta for real emission graph
c          v(mu,1:nv,1)      nv lepton momenta for real emission graph
c
c OUTPUT   np-1 parton momenta for Born amplitude
c          p(mu,1:np-1,2...)    tilde momenta
c          lepton momenta for born amplitude
c          v(mu,1:nv,2.....)    tilde momenta
c
c          xyuz               x_ija in the subtraction terms
c          xyuz( ,kin+9)        1-xyuz for ff and id/fi kinematics
c
c
c    diplole type      ptilde configuration         D(ij,k) -> pt conf.
c
c     initial,initial            2-7                   3a+i-4
c     final,final                8-10                   k+5
c     i,f  / f,i                 11-16                 3a+i+j+1
c**************************************************************************
      subroutine qcdVVJJptildemake(p,xyuz,v)

      implicit none

#include "global.inc"

      integer jj,ii,mu,i,j,k,a,b, kin, u

      double precision p(0:3,max_p,max_kin), v(0:3,max_v,max_kin), 
     &    pipb,papb,pipj,pjpa,pipk,pipa,pjpk,
     &    kk(0:3),kt(0:3),kkt(0:3),k2,kkt2,
     &    xy,vstore(4),pstore(5), xyuz(2,2:max_kin),one
      real*8 alpha
      parameter (one=1d0,alpha=1d0)
c
      logical ldebug
      parameter (ldebug=.false.)
      
      double precision dotrr
      external dotrr


      if (ldebug) then
c store born kinematics energies
        do ii=1,4
          pstore(ii) = p(0,ii,1)
          vstore(ii) = v(0,ii,1)
        enddo
        pstore(5) = p(0,5,1)
      endif

c  set the unused ptilde-momenta to zero
         do ii = 2,max_kin
            do jj=1,max_p
            do mu = 0,3
               p(mu,jj,ii) = 0d0
             enddo
            enddo
         enddo
         do ii = 2,max_kin
            do jj=1,max_v
            do mu = 0,3
               v(mu,jj,ii) = 0d0
             enddo
            enddo
         enddo


c       initial initial D^ai,b
      do a=1,2
        b=3-a
        do i=3,5
          kin = 3*a+i-4
          
          papb = dotrr(p(0,a,1),p(0,b,1))
          pipa = dotrr(p(0,i,1),p(0,a,1))
          pipb = dotrr(p(0,i,1),p(0,b,1))
          xy = one - (pipa + pipb)/papb
          
          xyuz(1,kin) = xy
          xyuz(2,kin) = 0d0
          
          do mu=0,3
            p(mu,a,kin) = xy* p(mu,a,1)
            p(mu,b,kin) = p(mu,b,1)
            kk(mu) = p(mu,a,1) + p(mu,b,1) - p(mu,i,1)
            kt(mu) = p(mu,a,kin)+p(mu,b,1)
            kkt(mu)= kk(mu)+kt(mu)
          enddo
          
          k2 = dotrr(kk(0),kk(0))
          kkt2=dotrr(kkt(0),kkt(0))
          do ii=3,4
            jj= ii/i + ii           ! jj =  (ii < i) ? ii : ii+1
            do mu=0,3
              p(mu,ii,kin) = p(mu,jj,1) - 2d0*dotrr(p(0,jj,1),kkt(0))/kkt2*kkt(mu)
     &                       + 2d0* dotrr(p(0,jj,1),kk(0))/k2*kt(mu)
              v(mu,ii,kin) = v(mu,ii,1) - 2d0*dotrr(v(0,ii,1),kkt(0))/kkt2*kkt(mu)
     &                       + 2d0* dotrr(v(0,ii,1),kk(0))/k2*kt(mu)
              v(mu,ii-2,kin)=v(mu,ii-2,1) - 2d0*dotrr(v(0,ii-2,1),kkt(0))/kkt2*kkt(mu)
     &                       + 2d0* dotrr(v(0,ii-2,1),kk(0))/k2*kt(mu)
            enddo
          enddo
        enddo
      enddo


c       final final kinematics   D_ij,k  ^= D_ji,k
c       p(3) = pij,  p(4) = pk, z is calculated for i<j
      do kin=8,10
        if(kin.eq.8) then
            i = 4
            j = 5
            k = 3
        elseif(kin.eq.9) then
            i = 3
            j = 5
            k = 4
        else
            i = 3
            j = 4
            k = 5
        endif

        pipj=dotrr(p(0,i,1),p(0,j,1))
        pipk=dotrr(p(0,i,1),p(0,k,1))
        pjpk=dotrr(p(0,j,1),p(0,k,1))

        xy=pipj/(pipj+pipk+pjpk)
        
        do mu=0,3
          p(mu,1,kin) = p(mu,1,1)
          p(mu,2,kin) = p(mu,2,1)
          p(mu,3,kin) = xy/(xy-one)*p(mu,k,1)+p(mu,i,1)+p(mu,j,1)  ! pij
          p(mu,4,kin) = p(mu,k,1)/(one-xy)  ! pk
          do jj=1,4
            v(mu,jj,kin) = v(mu,jj,1)
          enddo
        enddo
        xyuz(1,kin)  = xy
        xyuz(1,kin+9)= one-xy
        xyuz(2,kin)  = pipk/(pjpk+pipk)
        xyuz(2,kin+9)= pjpk/(pjpk+pipk)
      enddo
        
      
c initail final  D^ai_k  AND  final initial D^a_ij   have identical kinematics
c and momenta are invariant under i <> j or  i <> k
c z,u will be calculated for i < j; i < k,
c 1-z , 1-u has to be used in the other case
c p(3) = ij or k
c p(4) = not participating parton
      do a = 1,2
        do ii=1,3
          if(ii.eq.1) then
              i = 4
              j = 5
              u = 3     ! unchanged particle
          elseif(ii.eq.2) then
              i = 3
              j = 5
              u = 4
          else
              i = 3
              j = 4
              u = 5
          endif
          kin = 3*a + i+j+1
          
          pipa = dotrr(p(0,i,1),p(0,a,1))
          pjpa = dotrr(p(0,j,1),p(0,a,1))
          pipj = dotrr(p(0,i,1),p(0,j,1))
          
          xyuz(1,kin+9) = pipj/(pipa+pjpa)
          xyuz(1,kin)  = one-xyuz(1,kin+9)
          xyuz(2,kin)  = pipa/(pipa+pjpa)
          xyuz(2,kin+9)= pjpa/(pipa+pjpa)
          
          xy = xyuz(1,kin)
          
          do mu=0,3
            p(mu, a, kin) = xy*p(mu, a, 1)
            p(mu,3-a,kin) =    p(mu,3-a,1)
            p(mu, 3, kin) = p(mu,i,1)+p(mu,j,1)-xyuz(1,kin+9)*p(mu,a,1)
            p(mu, 4, kin) = p(mu,u,1)
            do jj=1,4
              v(mu,jj,kin) = v(mu,jj,1)
            enddo
          enddo
        enddo
      enddo
          

      if (ldebug) then
      print*,'real emission kinematic'
 10      format( " p(", i1, ") = ", 4(f10.3, 2x) )
 20      format( " v(", i1, ") = ", 4(f10.3, 2x) )
c
         do i = 1, 4
            write(6,20) i, v(0,i,1), v(1,i,1), v(2,i,1), v(3,i,1)
         end do
         print*," parton momenta "
         do i = 1, 5
            write(6,10) i, p(0,i,1), p(1,i,1), p(2,i,1), p(3,i,1)
         end do
        do ii=1,16
        print*
          print*, "Dipole kinematic", ii
         do i = 1, 4
            write(6,20) i, v(0,i,ii), v(1,i,ii), v(2,i,ii), v(3,i,ii)
         end do
         print*," parton momenta "
         do i = 1, 4
            write(6,10) i, p(0,i,ii), p(1,i,ii), p(2,i,ii), p(3,i,ii)
         end do   
          print*, 'xyuz',(xyuz(i,ii),i=1,2)
          print*, "Massless particles:"
          do jj=1,4
            print*,dotrr(p(0,jj,ii),p(0,jj,ii)), " = 0 ?"
          enddo
          do jj=1,4
            print*,dotrr(v(0,jj,ii),v(0,jj,ii)), " = 0 ?"
          enddo
          print*,"Momentum conservation:"
          do mu=0,3
            kk(mu) = p(mu,1,ii)+p(mu,2,ii)
            if(ii.eq.1) kk(mu) = kk(mu)-p(mu,5,ii)
            do jj=3,4
              kk(mu)= kk(mu) - p(mu,jj,ii) - v(mu,jj,ii) - v(mu,jj-2,ii)
            enddo
            print*,kk(mu), " = 0 ?"
          enddo
          print*, "------------------------------------------------"
        enddo
        print*, "Did't change born kinematics? Test energies"
        print*, "pstore:  ", pstore
        print*, "changed? ", (sqrt(p(1,ii,1)**2+p(2,ii,1)**2+p(3,ii,1)**2) ,ii=1,5)
        print*, "vstore:  ", vstore
        print*, "changed? ", (sqrt(v(1,ii,1)**2+v(2,ii,1)**2+v(3,ii,1)**2) ,ii=1,4)
      endif

      end


c********************** subroutine vvjptildemake *************************
c
c    this subroutine assigns ptilde momenta to the partons and the final
c    state leptons in notation of Catani and Seymour.
c
c INPUT    p(mu,1:np,1)      np parton momenta for real emission graph
c          v(mu,1:nv,1)      nv lepton momenta for real emission graph
c
c OUTPUT   np-1 parton momenta for Born amplitude
c          p(mu,1:np-1,2...)    tilde momenta

c          lepton momenta for born amplitude
c          v(mu,1:nv,2.....)    tilde momenta

c
c          xuz               x_ija in the subtraction terms
c
c**************************************************************************
c
c based on wgajptildemake, but leaves position of momenta intact
c
c**************************************************************************
      subroutine vvjptildemake(p,xuz,v)

      implicit none

#include "global.inc"

      integer jj,ii,mu

      double precision pin(0:3),pk(0:3),pa(0:3),pj(0:3),v(0:3,max_v,max_kin),
     &    p(0:3,max_p,max_kin),pb(0:3),p4(0:3),p3(0:3),ptildea(0:3),ptildeij(0:3),
     &    pipb,papb,k(0:3),Kt(0:3),pipj,xija,zi,papj,papi
c
      logical ldebug
      parameter (ldebug=.false.)

      double precision ptildeai(0:3),ptildek(0:3)
     &    ,pipk,pkpa,pipa
      double precision xika,ui,xuz(2,2:max_kin),xiab,one
      parameter (one=1d0)

      double precision dotrr
      external dotrr

c  set the unused ptilde-momenta to zero
         do ii = 2,max_kin
            do jj=1,max_p
            do mu = 0,3
               p(mu,jj,ii) = 0d0
             enddo
            enddo
         enddo
         do ii = 2,max_kin
            do jj=1,max_v
            do mu = 0,3
               v(mu,jj,ii) = 0d0
             enddo
            enddo
         enddo

c  store 3rd and 4th momenta
         do mu=0,3
         p3(mu)=p(mu,3,1)
         p4(mu)=p(mu,4,1)
         enddo

c   get ptilde kinematics for p p -> j j + EW
c   this involves kinematics of Catani Seymour 5.67,5.147

c Dipoles 5.67

      do mu=0,3
      pa(mu) =p(mu,1,1)
      pin(mu)=p(mu,4,1)
      pk(mu) =p(mu,3,1)
      enddo

      pipk=dotrr(pin,pk)
      pkpa=dotrr(pk,pa)
      pipa=dotrr(pin,pa)

      xika=one - pipk/(pkpa+pipa)! (pkpa+pipa-pipk)/(pkpa+pipa)

      do mu=0,3
      ptildeai(mu)=xika*pa(mu)
      ptildek(mu) =pk(mu)+pin(mu)-(one-xika)*pa(mu)
      enddo
      ui=pipa/(pipa+pkpa)


c   generate 1. output
c
      do mu=0,3
      p(mu,1,2)=ptildeai(mu)       
      p(mu,2,2)=p(mu,2,1)
      p(mu,3,2)=ptildek(mu)
      do jj=1,n_v
         v(mu,jj,2)=v(mu,jj,1)
      enddo
      enddo

      xuz(1,2)=xika
      xuz(2,2)=ui

c
c start new dipole
c
      do mu=0,3
      pa(mu)=p(mu,2,1)
      pin(mu)=p(mu,4,1)
      pk(mu)=p(mu,3,1)
      enddo

      pipk=dotrr(pin,pk)
      pkpa=dotrr(pk,pa)
      pipa=dotrr(pin,pa)

      xika=one - pipk/(pkpa+pipa)!(pkpa+pipa-pipk)/(pkpa+pipa)

      do mu=0,3
      ptildeai(mu)=xika*pa(mu)
      ptildek(mu) =pk(mu)+pin(mu)-(one-xika)*pa(mu)
      enddo
      ui=pipa/(pipa+pkpa)

c   generate 2. output
      do mu=0,3
      p(mu,1,3)=p(mu,1,1)
      p(mu,2,3)=ptildeai(mu)       
      p(mu,3,3)=ptildek(mu)
      do jj=1,n_v
         v(mu,jj,3)=v(mu,jj,1)
      enddo
      enddo
      xuz(1,3)=xika
      xuz(2,3)=ui


c
c   start new dipole
c
      do mu=0,3
      pa(mu)=p(mu,1,1)
      pin(mu)=p(mu,4,1)
      pb(mu)=p(mu,2,1)
      enddo

      pipb=dotrr(pin,pb)
      papb=dotrr(pa,pb)
      pipa=dotrr(pin,pa)
      xiab=one- (pipa+pipb)/papb !(papb-pipa-pipb)/papb

      do mu=0,3
      ptildeai(mu)=xiab*pa(mu)
      K(mu)=pa(mu)+pb(mu)-pin(mu)
      Kt(mu)=ptildeai(mu)+pb(mu)
      enddo

      Call boostkinematics(k,kt,p(0,1,1),v(0,1,1),p(0,1,4),v(0,1,4))

c   generate 3. output
      do mu=0,3
      p(mu,1,4)=ptildeai(mu)  
      p(mu,2,4)=p(mu,2,1)
      enddo
      xuz(1,4)=xiab
      xuz(2,4)=0d0 ! unused

c   start new dipole
c
      do mu=0,3
      pa(mu)=p(mu,2,1)
      pin(mu)=p(mu,4,1)
      pb(mu)=p(mu,1,1)
      enddo

      pipb=dotrr(pin,pb)
      papb=dotrr(pa,pb)
      pipa=dotrr(pin,pa)
      xiab= one- (pipa+pipb)/papb!(papb-pipa-pipb)/(papb)

      do mu=0,3
      ptildeai(mu)=xiab*pa(mu)
      K(mu)=pa(mu)+pb(mu)-pin(mu)
      Kt(mu)=ptildeai(mu)+pb(mu)
      enddo

      Call boostkinematics(k,kt,p(0,1,1),v(0,1,1),p(0,1,5),v(0,1,5))

c   generate 4. output
      do mu=0,3
      p(mu,1,5)=p(mu,1,1)
      p(mu,2,5)=ptildeai(mu)  
      enddo
      xuz(1,5)=xiab
      xuz(2,5)=0d0 ! unused


c Dipoles 5.67

      do mu=0,3
      pa(mu) =p(mu,1,1)
      pin(mu)=p(mu,3,1)
      pk(mu) =p(mu,4,1)
      enddo

      pipk=dotrr(pin,pk)
      pkpa=dotrr(pk,pa)
      pipa=dotrr(pin,pa)

      xika=one - pipk/(pkpa+pipa)! (pkpa+pipa-pipk)/(pkpa+pipa)

      do mu=0,3
      ptildeai(mu)=xika*pa(mu)
      ptildek(mu) =pk(mu)+pin(mu)-(one-xika)*pa(mu)
      enddo
      ui=pipa/(pipa+pkpa)

c   generate 5. output
c
      do mu=0,3
      p(mu,1,6)=ptildeai(mu)       
      p(mu,2,6)=p(mu,2,1)
      p(mu,3,6)=ptildek(mu)
      do jj=1,n_v
         v(mu,jj,6)=v(mu,jj,1)
      enddo
      enddo
      xuz(1,6)=xika
      xuz(2,6)=ui

c
c start new dipole
c
      do mu=0,3
      pa(mu)=p(mu,2,1)
      pin(mu)=p(mu,3,1)
      pk(mu)=p(mu,4,1)
      enddo

      pipk=dotrr(pin,pk)
      pkpa=dotrr(pk,pa)
      pipa=dotrr(pin,pa)

      xika=one - pipk/(pkpa+pipa)! (pkpa+pipa-pipk)/(pkpa+pipa)

      do mu=0,3
      ptildeai(mu)=xika*pa(mu)
      ptildek(mu) =pk(mu)+pin(mu)-(one-xika)*pa(mu)
      enddo
      ui=pipa/(pipa+pkpa)

c   generate 6. output
      do mu=0,3
      p(mu,1,7)=p(mu,1,1)
      p(mu,2,7)=ptildeai(mu)       
      p(mu,3,7)=ptildek(mu)
      do jj=1,n_v
         v(mu,jj,7)=v(mu,jj,1)
      enddo
      enddo
      xuz(1,7)=xika
      xuz(2,7)=ui

c
c   start new dipole
c
      do mu=0,3
      pa(mu)=p(mu,1,1)
      pin(mu)=p(mu,3,1)
      pb(mu)=p(mu,2,1)
      p(mu,3,1)=p4(mu)
      enddo

      pipb=dotrr(pin,pb)
      papb=dotrr(pa,pb)
      pipa=dotrr(pin,pa)
      xiab= one- (pipa+pipb)/papb!(papb-pipa-pipb)/papb

      do mu=0,3
      ptildeai(mu)=xiab*pa(mu)
      K(mu)=pa(mu)+pb(mu)-pin(mu)
      Kt(mu)=ptildeai(mu)+pb(mu)
      enddo

      Call boostkinematics(k,kt,p(0,1,1),v(0,1,1),p(0,1,8),v(0,1,8))

      do mu=0,3
      p(mu,3,1)=p3(mu)
      enddo

c   generate 7. output
      do mu=0,3
      p(mu,1,8)=ptildeai(mu)       
      p(mu,2,8)=p(mu,2,1)
      enddo
      xuz(1,8)=xiab
      xuz(2,8)=0d0 ! unused

c   start new dipole
c
      do mu=0,3
      pa(mu)=p(mu,2,1)
      pin(mu)=p(mu,3,1)
      pb(mu)=p(mu,1,1)
      p(mu,3,1)=p4(mu)
      enddo

      pipb=dotrr(pin,pb)
      papb=dotrr(pa,pb)
      pipa=dotrr(pin,pa)
      xiab= one- (pipa+pipb)/papb!(papb-pipa-pipb)/(papb)

      do mu=0,3
      ptildeai(mu)=xiab*pa(mu)
      K(mu)=pa(mu)+pb(mu)-pin(mu)
      Kt(mu)=ptildeai(mu)+pb(mu)
      enddo

      Call boostkinematics(k,kt,p(0,1,1),v(0,1,1),p(0,1,9),v(0,1,9))

      do mu=0,3
      p(mu,3,1)=p3(mu)
      enddo

c   generate 8. output
      do mu=0,3
      p(mu,1,9)=p(mu,1,1)
      p(mu,2,9)=ptildeai(mu)       
      enddo
      xuz(1,9)=xiab
      xuz(2,9)=0d0 ! unused


c Dipoles 5.40

      do mu=0,3
      pa(mu) =p(mu,1,1)
      pin(mu)=p(mu,3,1)
      pj(mu) =p(mu,4,1)
      enddo

      papj=dotrr(pa,pj)
      papi=dotrr(pa,pin)
      pipj=dotrr(pin,pj)

      xija=one - pipj/(papi+papj)! (papi+papj-pipj)/(papi+papj)

      do mu=0,3
      ptildea(mu)=xija*pa(mu)
      ptildeij(mu) =pin(mu)+pj(mu)-(one-xija)*pa(mu)
      enddo
      zi=papi/(papi+papj)

c   generate 9. output
c
      do mu=0,3
      p(mu,1,10)=ptildea(mu)       !this is the intermediate gluon, dipole (13,4)
      p(mu,2,10)=p(mu,2,1)
      p(mu,3,10)=ptildeij(mu)
      do jj=1,n_v
         v(mu,jj,10)=v(mu,jj,1)
      enddo
      enddo
      xuz(1,10)=xija
      xuz(2,10)=zi


      do mu=0,3
      pa(mu) =p(mu,2,1)
      pin(mu)=p(mu,3,1)
      pj(mu) =p(mu,4,1)
      enddo

      papj=dotrr(pa,pj)
      papi=dotrr(pa,pin)
      pipj=dotrr(pin,pj)

      xija=one - pipj/(papi+papj)! (papi+papj-pipj)/(papi+papj)

      do mu=0,3
      ptildea(mu)=xija*pa(mu)
      ptildeij(mu) =pin(mu)+pj(mu)-(one-xija)*pa(mu)
      enddo
      zi=papi/(papi+papj)

c   generate 10. output
c
      do mu=0,3
      p(mu,1,11)=p(mu,1,1)       !this is the intermediate gluon, dipole (13,4)
      p(mu,2,11)=ptildea(mu)
      p(mu,3,11)=ptildeij(mu)
      do jj=1,n_v
         v(mu,jj,11)=v(mu,jj,1)
      enddo
      enddo
      xuz(1,11)=xija
      xuz(2,11)=zi

      end


      subroutine boostkinematics(k,kt,p,v,pt,vt)

      implicit none

#include "global.inc"

         integer mu,i,j
         ! todo: 4->n_p or max_p in whole function
         ! or merge with boostn functino in kinematics_aux.F
         double precision k(0:3),kt(0:3),p(0:3,4),pt(0:3,4),
     &         v(0:3,max_v),vt(0:3,max_v),dotrr,kpkt(0:3),kpkt2,k2

      logical ldebug
      Parameter (ldebug=.false.)

      external dotrr

         if(ldebug) print*, 'debug in ptilde.F boostkinematics'

         do i=1,4
            do mu=0,3
               pt(mu,i)=0d0
            enddo
         enddo

         do mu=0,3
            kpkt(mu)=k(mu)+kt(mu)
         enddo

         kpkt2=kpkt(0)**2-kpkt(1)**2-kpkt(2)**2-kpkt(3)**2
         k2=k(0)**2-k(1)**2-k(2)**2-k(3)**2

         do mu=0,3
            do i=1,4

               pt(mu,i)=p(mu,i)-2d0*dotrr(p(0,i),kpkt(0))/kpkt2*kpkt(mu)
     &      + 2d0* dotrr(p(0,i),k(0))/k2*kt(mu)
            enddo

            do i=1,n_v
               vt(mu,i)=v(mu,i)-2d0*dotrr(v(0,i),kpkt(0))/kpkt2*kpkt(mu)
     &      + 2d0* dotrr(v(0,i),k(0))/k2*kt(mu)
            enddo
         enddo

         if (ldebug) then
            do i=1,3
               do j=1,i
                  print*, i,j, "mom scps parton   ", dotrr(p(0,i),p(0,j)),dotrr(pt(0,i),pt(0,j)),
     &       dotrr(p(0,i),p(0,j))/dotrr(pt(0,i),pt(0,j))
               enddo
            enddo

            do i=1,n_v
               do j=1,i
                  print*, i,j, "mom scps leptons  ", dotrr(v(0,i),v(0,j)),dotrr(vt(0,i),vt(0,j)),
     &       dotrr(v(0,i),v(0,j))/dotrr(vt(0,i),vt(0,j))
               enddo
            enddo

            do i=1,3
               do j=1,n_v
                  print*, i,j, "mom scps parts leps", dotrr(p(0,i),v(0,j)), dotrr(pt(0,i),vt(0,j)),
     &       dotrr(p(0,i),v(0,j))/dotrr(pt(0,i),vt(0,j))
               enddo
            enddo
            print*, "------------------------------------------------------"
         endif

      end


c********************** subroutine DIBOSptildemake *************************
c
c    this subroutine assigns ptilde momenta to the partons and the final
c    state leptons in notation of Catani and Seymour.
c
c INPUT    p(mu,1:np,1)      np parton momenta for real emission graph
c          v(mu,1:nv,1)      nv lepton momenta for real emission graph
c
c OUTPUT   np-1 parton momenta for Born amplitude
c          p(mu,1:np-1,2)    index=2: emission off 1-2 line
c          p(mu,1:np-1,3)    index=3: emission off 3-4 line

c          lepton momenta for born amplitude
c          v(mu,1:nv,2)      index=2: emission off 1-2 line
c          v(mu,1:nv,3)      index=3: emission off 3-4 line
c
c          xuz               x_ija in the subtraction terms
c
c**************************************************************************

      subroutine DIBOSptildemake(N,p,xuz,v)

      implicit none

#include "global.inc"

      integer mu, a, i, j, N, ii, b
      real*8 p(0:3,max_p,max_kin), xuz, v(0:3,max_v,max_kin)
      real*8 omx, pipa, pipb, papb
      real*8 vK(1:n_v), vKtilde(1:n_v), KKtilde

      real*8 tpi
      real*8 K(0:4), Ktilde(0:4)
      double precision dotrr, temp(0:3), temp2(0:3)
      parameter (tpi=2d0*3.141592653589793)

      logical ldebug
      parameter (ldebug=.false.)

      external dotrr

      pipa = dotrr(p(0,3,1),p(0,1,1))
      pipb = dotrr(p(0,3,1),p(0,2,1))
      papb = dotrr(p(0,1,1),p(0,2,1))

      omx = (pipa + pipb)/papb
      xuz = 1d0 - omx

      do mu = 0,3
         K(mu) = p(mu,1,1) + p(mu,2,1) - p(mu,3,1)
      enddo
      K(4) = dotrr(K,K)

      do ii = 2,N
         if (ii.eq.2) then        ! set mom id for radiation off 1,2 line
            a = 1
            b = 2
            i = 3
         elseif (ii.eq.3) then
            a = 2
            b = 1
            i = 3
         else
            stop
         endif

         do mu = 0,3
            Ktilde(mu) = xuz*p(mu,a,1) + p(mu,b,1)
         enddo

         Ktilde(4) = dotrr(Ktilde,Ktilde)
         KKtilde = dotrr(K,Ktilde)

         if (ldebug) then
            print*," K^2/Ktilde^2 = 1 ? ", K(4)/Ktilde(4)
         endif

         if (ldebug) print*," assign tilde momenta for parton # ",a,i
         do mu = 0,3
            p(mu,i,ii) = 0.0d0
            p(mu,a,ii) = xuz*p(mu,a,1)
            p(mu,b,ii) = p(mu,b,1)
         enddo

         if (ldebug) print*," assign tilde momenta for all leptons: "

         do j = 1,n_v
            vK(j) = dotrr(v(0,j,1),K(0))
            vKtilde(j) = dotrr(v(0,j,1),Ktilde(0))

            do mu = 0,3
               v(mu,j,ii) = v(mu,j,1) - 2d0*(vK(j) + vKtilde(j))/
     &            (K(4) + Ktilde(4) + 2d0*KKtilde)*(K(mu) + Ktilde(mu))
     &            +2d0*vK(j)/K(4) * Ktilde(mu)
            enddo
         enddo

c check if ktilde_i*ktilde_j = k_i*k_j 
         if (ldebug) then
            do j = 1,n_v-1
               do mu = j+1,n_v
                  print*,"momentum configuration ",ii,": "
                  print*," vtilde",j,"*vtilde",mu,
     &                   "/v",j,"*v",mu," = 1 ?", 
     &            dotrr(v(0,j,ii),v(0,mu,ii))
     &            /dotrr( v(0,j,1),v(0,mu,1))
               enddo
            enddo
        endif

c check if v_i^2 = vtile_i^2
        if (ldebug) then
            do j =1,n_v
               if(dotrr(v(0,j,1),v(0,j,1)).ge.1d-7.or.
     &            dotrr(v(0,j,ii),v(0,j,ii)).ge.1d-7) then
               print*," v",j,"^2 = ", dotrr(v(0,j,1),v(0,j,1))
               print*,"momentum configuration ",ii,": "
               print*," vtilde",j,"^2 = ", dotrr(v(0,j,ii),v(0,j,ii))
               endif
            enddo
        endif

c check if (v1+v2)^2 = (vtilde1+vtilde2)^2
        if (ldebug) then
            print*,"momentum configuration ",ii,": "
            do mu = 0,3
               temp(mu) = v(mu,1,1)+v(mu,2,1)
               temp2(mu)= v(mu,1,ii)+v(mu,2,ii)
            enddo
            print*," (v1 + v2)^2/(vtilde1 + vtilde2)^2 = "
     &                  ,dotrr(temp,temp)/dotrr(temp2,temp2)
            do mu = 0,3
               temp(mu) = v(mu,3,1)+v(mu,4,1)
               temp2(mu)= v(mu,3,ii)+v(mu,4,ii)
            enddo
            print*," (v3 + v4)^2/(vtilde3 + vtilde4)^2 = "
     &                  ,dotrr(temp,temp)/dotrr(temp2,temp2)
        endif

c check momentum conservation
         if (ldebug) then
            do mu = 0,3
               temp(mu) = p(mu,a,1) + p(mu,b,1) - p(mu,i,1)
     &               - v(mu,1,1) - v(mu,2,1) - v(mu,3,1) - v(mu,4,1)
               print*," mom.cons. for original momentum configuration "
               print*," mu=",mu,": sum= ",temp(mu)
               temp(mu) = p(mu,a,ii) + p(mu,b,ii)
     &               - v(mu,1,ii) - v(mu,2,ii) - v(mu,3,ii) - v(mu,4,ii)
               print*," mom.cons. for new mom. configuration "
               print*,"mu=",mu,": sum= ",temp(mu)
            enddo

check masslessness of ptilde_a and ptilde_j
            print*," p_a^2 = ",p(0,a,ii)**2-p(1,a,ii)**2-
     1                         p(2,a,ii)**2-p(3,a,ii)**2
            print*," p_i^2 = ",p(0,i,ii)**2-p(1,i,ii)**2-
     1                         p(2,i,ii)**2-p(3,i,ii)**2
            do  j = 1 , n_v
               print*,"j = ",j, " v_j^2 = ",dotrr(v(0,j,ii),v(0,j,ii))
            enddo
         endif

      enddo

      if ( ldebug ) then
 10      format( " p(", i1, ") = ", 4(f10.3, 2x) )
c
         do j = 1,N
            if (j.eq.1) print*," (n+1) parton momenta "
            if (j.eq.2) print*," n-parton momenta for upper NLO corr."
            if (j.eq.3) print*," n-parton momenta for lower NLO corr."
            do i = 1, 3-min(j-1,1)
               write(6,10) i, p(0,i,j), p(1,i,j), p(2,i,j), p(3,i,j)
            enddo
         enddo

      end if
     
      return
      end


c*************************************************************************
      subroutine HJJJptildemake(N,p,xuz,v)
c*************************************************************************
c
c assign ptilde momenta in notation of Catani&Seymour
c
c  input:  p(mu,1:np,1)      np parton momenta for real emission graph
c          v(mu,nv,1)
c  output:                   np-1 parton momenta for Born amplitude
c          p(mu,1:np,nd)     nd is the dipole number
c
c          xuz(1,nd)      x_ija,y_ijk etc
c          xuz(2,nd)      u or z
c
c at the same time, fill the legoplot variables for the external momenta
c
c       Dieter Zeppenfeld, <dieter@pheno.physics.wisc.edu>
c       Initial version:  November 2002
c       Last modified: March 11, 2005
c*************************************************************************

      implicit none

#include "global.inc"

      integer mu,N,ii,j,i
      real*8 p(0:3,max_p,1:max_kin),xuz(2,2:max_kin),v(0:3,max_v,1:max_kin)

      integer N1,npart,ic
      integer kp(2:max_kin),jp(2:max_kin),ip(2:max_kin),nd
      real*8 tpi
      parameter (tpi=2d0*3.141592653589793)
      logical ldebug
      parameter (ldebug=.false.)
      real*8 ph(0:3)            ! 4-momentum of higgs
c     reconstructed from decay momentum

      N1=n_p-nx
      if (N.eq.1) then
         npart = N1-1
      else
         npart = N1
      endif
c
c  set superfluous variables to zero
      do ic = npart+1,N1
         do ii = 1,max_kin
            do mu = 0,3
               p(mu,nx+ic,ii) = 0
            enddo
         enddo
      enddo
c
c  put call to transform subroutine here
c     map ip,kp,jp to nd
c     loop over dipole configurations.
c     D^56,4 in the notation of campbell and ellis
c
      if(N.ge.2) then 
c     fill ip,kp,jp arrays
c
         call create_ikj(ip,kp,jp)   
c
c           
         do nd = 2,N
            call transform(nd,p,v,xuz,ip(nd),kp(nd),jp(nd))
         enddo

      endif


      if ( ldebug ) then
 10      format( ' p(', i1, ') = ', 4(f10.3, 2x) )
c
         do j = 1,N
            if (j.eq.1) print*,' (n+1) parton momenta '
            if (j.ge.2) print*,' n-parton momenta for',j,
     &   'momentum config'
            do i = 1, n_p-min(j-1,1)
               write(6,10) i, p(0,i,j), p(1,i,j), p(2,i,j), p(3,i,j)
            end do

         enddo         
      endif
      end

c********************** subroutine WGAJptildemake *************************
c**************************************************************************

      subroutine TriJetptildemake(N,p,xuz,v,nv)

      implicit none

      include 'global.inc'

      integer jj,ii,mu,N,nv

      double precision pin(0:3),pk(0:3),pa(0:3),pj(0:3),v(0:3,max_v,max_kin),
     &    p(0:3,max_p,max_kin),pb(0:3),p4(0:3),p3(0:3),ptildea(0:3),ptildeij(0:3),
     &    pipb,papb,k(0:3),Kt(0:3),pipj,xija,zi,papj,papi,
     &    vstore(0:3,max_v)
c
      logical ldebug
      parameter (ldebug=.false.)

      double precision ptildeai(0:3),ptildek(0:3)
     &    ,pipk,pkpa,pipa
      double precision xika,ui,xuz(2,2:max_kin),xiab,one
      parameter (one=1d0)

      double precision dotrr
      external dotrr

      if (ldebug) then
c store elw parameters
      do mu=0,3
      do ii=1,nv
       vstore(mu,ii)=v(mu,ii,1)
      enddo
      enddo
      endif

c  set the unused ptilde-momenta to zero
         do ii = 2,max_kin
            do jj=1,max_p
            do mu = 0,3
               p(mu,jj,ii) = 0d0
             enddo
            enddo
         enddo
         do ii = 2,max_kin
            do jj=1,max_v
            do mu = 0,3
               v(mu,jj,ii) = 0d0
             enddo
            enddo
         enddo


c  store 3rd and 4th momenta
         do mu=0,3
         p4(mu)=p(mu,4,1)
         p3(mu)=p(mu,3,1)
         enddo

c

c   get ptilde kinematics for u~u~ > e- ve~ a u~ d~
c   this involves kinematics of Catani Seymour 5.67,5.147

c Dipoles 5.67

      do mu=0,3
      pa(mu) =p(mu,1,1)
      pin(mu)=p(mu,4,1)
      pk(mu) =p(mu,3,1)
      enddo

      pipk=dotrr(pin,pk)
      pkpa=dotrr(pk,pa)
      pipa=dotrr(pin,pa)

      xika=one - pipk/(pkpa+pipa)! (pkpa+pipa-pipk)/(pkpa+pipa)

      do mu=0,3
      ptildeai(mu)=xika*pa(mu)
      ptildek(mu) =pk(mu)+pin(mu)-(one-xika)*pa(mu)
      enddo
      ui=pipa/(pipa+pkpa)


c   generate 1. output
c
      do mu=0,3
      p(mu,1,2)=ptildeai(mu)       !this is the intermediate gluon, dipole (13,4)
      p(mu,2,2)=p(mu,2,1)
      p(mu,3,2)=ptildek(mu)
      do ii=1,nv
        v(mu,ii,2)=v(mu,ii,1)
      enddo
      enddo

      xuz(1,2)=xika
      xuz(2,2)=ui

      if (ldebug) then
      print*, "New Dipole 1"
      print*, "Momentum conservation of real em. kinematics:"
      do mu=0,3
      print*, p(mu,1,1)+p(mu,2,1)-p(mu,3,1)-p(mu,4,1)-v(mu,1,1)-v(mu,2,1)-v(mu,3,1),"=0?"
      enddo
      print*, "Momentum conservation of Dipole kinematics:"
      do mu=0,3
      print*, p(mu,1,2)+p(mu,2,2)-p(mu,3,2)-v(mu,1,2)-v(mu,2,2)-v(mu,3,2),"=0?"
      enddo
      print*, "------------------------------------------------"
      endif

c
c start new dipole
c
      do mu=0,3
      pa(mu)=p(mu,2,1)
      pin(mu)=p(mu,4,1)
      pk(mu)=p(mu,3,1)
      enddo

      pipk=dotrr(pin,pk)
      pkpa=dotrr(pk,pa)
      pipa=dotrr(pin,pa)

      xika=one - pipk/(pkpa+pipa)!(pkpa+pipa-pipk)/(pkpa+pipa)

      do mu=0,3
      ptildeai(mu)=xika*pa(mu)
      ptildek(mu) =pk(mu)+pin(mu)-(one-xika)*pa(mu)
      enddo
      ui=pipa/(pipa+pkpa)

c   generate 2. output
      do mu=0,3
      p(mu,1,3)=ptildeai(mu)       !this is the intermediate gluon, dipole (23,4)
      p(mu,2,3)=p(mu,1,1)
      p(mu,3,3)=ptildek(mu)
      do ii=1,nv
        v(mu,ii,3)=v(mu,ii,1)
      enddo
      enddo
      xuz(1,3)=xika
      xuz(2,3)=ui

c    test output
      if (ldebug) then
      print*, "New Dipole 2"
      print*, "Momentum conservation of real em. kinematics:"
      do mu=0,3
      print*, p(mu,1,1)+p(mu,2,1)-p(mu,3,1)-p(mu,4,1)-v(mu,1,1)-v(mu,2,1)-v(mu,3,1),"=0?"
      enddo
      print*, "Momentum conservation of Dipole kinematics:"
      do mu=0,3
      print*, p(mu,1,3)+p(mu,2,3)-p(mu,3,3)-v(mu,1,3)-v(mu,2,3)-v(mu,3,3),"=0?"
      enddo
      print*, "------------------------------------------------"
      endif


c Dipoles 5.147
c
c   start new dipole
c
      do mu=0,3
      pa(mu)=p(mu,1,1)
      pin(mu)=p(mu,4,1)
      pb(mu)=p(mu,2,1)
      enddo

      pipb=dotrr(pin,pb)
      papb=dotrr(pa,pb)
      pipa=dotrr(pin,pa)
      xiab=one- (pipa+pipb)/papb !(papb-pipa-pipb)/papb

      do mu=0,3
      ptildeai(mu)=xiab*pa(mu)
      K(mu)=pa(mu)+pb(mu)-pin(mu)
      Kt(mu)=ptildeai(mu)+pb(mu)
      enddo

      Call boostkinematics(k,kt,p(0,1,1),v(0,1,1),p(0,1,4),v(0,1,4))

c   generate 3. output
      do mu=0,3
      p(mu,1,4)=ptildeai(mu)       !this is the intermediate gluon, dipole (13,2)
      p(mu,2,4)=p(mu,2,1)
      enddo
      xuz(1,4)=xiab
      xuz(2,4)=0d0 !!!!!!!!!!!!!


c    test dipole
c
      if (ldebug) then
      print*, "New Dipole 3"
      print*, "Momentum conservation of real em. kinematics:"
      do mu=0,3
      print*, p(mu,1,1)+p(mu,2,1)-p(mu,3,1)-p(mu,4,1)-v(mu,1,1)-v(mu,2,1)-v(mu,3,1),"=0?"
      enddo
      print*, "Momentum conservation of Dipole kinematics:"
      do mu=0,3
      print*, p(mu,1,4)+p(mu,2,4)-p(mu,3,4)-v(mu,1,4)-v(mu,2,4)-v(mu,3,4),"=0?"
      enddo
      print*, "------------------------------------------------"
      endif

c   start new dipole
c
      do mu=0,3
      pa(mu)=p(mu,2,1)
      pin(mu)=p(mu,4,1)
      pb(mu)=p(mu,1,1)
      enddo

      pipb=dotrr(pin,pb)
      papb=dotrr(pa,pb)
      pipa=dotrr(pin,pa)
      xiab= one- (pipa+pipb)/papb!(papb-pipa-pipb)/(papb)

      do mu=0,3
      ptildeai(mu)=xiab*pa(mu)
      K(mu)=pa(mu)+pb(mu)-pin(mu)
      Kt(mu)=ptildeai(mu)+pb(mu)
      enddo

      Call boostkinematics(k,kt,p(0,1,1),v(0,1,1),p(0,1,5),v(0,1,5))

c   generate 4. output
      do mu=0,3
      p(mu,1,5)=ptildeai(mu)       !this is the intermediate gluon, dipole (23,1)
      p(mu,2,5)=p(mu,1,1)
      enddo
      xuz(1,5)=xiab
      xuz(2,5)=0d0 !!!!!!!!!!!!!!!!

c    test dipole
c
      if (ldebug) then
      print*, "New Dipole 4"
      print*, "Momentum conservation of real em. kinematics:"
      do mu=0,3
      print*, p(mu,1,1)+p(mu,2,1)-p(mu,3,1)-p(mu,4,1)-v(mu,1,1)-v(mu,2,1)-v(mu,3,1),"=0?"
      enddo
      print*, "Momentum conservation of Dipole kinematics:"
      do mu=0,3
      print*, p(mu,1,5)+p(mu,2,5)-p(mu,3,5)-v(mu,1,5)-v(mu,2,5)-v(mu,3,5),"=0?"
      enddo
      print*, "------------------------------------------------"
      endif


c   get ptilde kinematics for u d > e- ve~ a u u and inicial state crossed
c   this involves kinematics of Catani Seymour 5.67,5.147

c Dipoles 5.67

      do mu=0,3
      pa(mu) =p(mu,1,1)
      pin(mu)=p(mu,3,1)
      pk(mu) =p(mu,4,1)
      enddo

      pipk=dotrr(pin,pk)
      pkpa=dotrr(pk,pa)
      pipa=dotrr(pin,pa)

      xika=one - pipk/(pkpa+pipa)! (pkpa+pipa-pipk)/(pkpa+pipa)

      do mu=0,3
      ptildeai(mu)=xika*pa(mu)
      ptildek(mu) =pk(mu)+pin(mu)-(one-xika)*pa(mu)
      enddo
      ui=pipa/(pipa+pkpa)

c   generate 5. output
c
      do mu=0,3
      p(mu,1,6)=ptildeai(mu)       !this is the intermediate gluon, dipole (13,4)
      p(mu,2,6)=p(mu,2,1)
      p(mu,3,6)=ptildek(mu)
      do ii=1,nv
        v(mu,ii,6)=v(mu,ii,1)
      enddo
      enddo
      xuz(1,6)=xika
      xuz(2,6)=ui

      if (ldebug) then
      print*, "New Dipole 5"
      print*, "Momentum conservation of real em. kinematics:"
      do mu=0,3
      print*, p(mu,1,1)+p(mu,2,1)-p(mu,3,1)-p(mu,4,1)-v(mu,1,1)-v(mu,2,1)-v(mu,3,1),"=0?"
      enddo
      print*, "Momentum conservation of Dipole kinematics:"
      do mu=0,3
      print*, p(mu,1,6)+p(mu,2,6)-p(mu,3,6)-v(mu,1,6)-v(mu,2,6)-v(mu,3,6),"=0?"
      enddo
      print*, "------------------------------------------------"
      endif

c
c start new dipole
c
      do mu=0,3
      pa(mu)=p(mu,2,1)
      pin(mu)=p(mu,3,1)
      pk(mu)=p(mu,4,1)
      enddo

      pipk=dotrr(pin,pk)
      pkpa=dotrr(pk,pa)
      pipa=dotrr(pin,pa)

      xika=one - pipk/(pkpa+pipa)! (pkpa+pipa-pipk)/(pkpa+pipa)

      do mu=0,3
      ptildeai(mu)=xika*pa(mu)
      ptildek(mu) =pk(mu)+pin(mu)-(one-xika)*pa(mu)
      enddo
      ui=pipa/(pipa+pkpa)

c   generate 6. output
      do mu=0,3
      p(mu,1,7)=ptildeai(mu)       !this is the intermediate gluon, dipole (23,4)
      p(mu,2,7)=p(mu,1,1)
      p(mu,3,7)=ptildek(mu)
      do ii=1,nv
        v(mu,ii,7)=v(mu,ii,1)
      enddo
      enddo
      xuz(1,7)=xika
      xuz(2,7)=ui

c    test output
      if (ldebug) then
      print*, "New Dipole 6"
      print*, "Momentum conservation of real em. kinematics:"
      do mu=0,3
      print*, p(mu,1,1)+p(mu,2,1)-p(mu,3,1)-p(mu,4,1)-v(mu,1,1)-v(mu,2,1)-v(mu,3,1),"=0?"
      enddo
      print*, "Momentum conservation of Dipole kinematics:"
      do mu=0,3
      print*, p(mu,1,7)+p(mu,2,7)-p(mu,3,7)-v(mu,1,7)-v(mu,2,7)-v(mu,3,7),"=0?"
      enddo
      print*, "------------------------------------------------"
      endif


c Dipoles 5.147
c
c   start new dipole
c
      do mu=0,3
      pa(mu)=p(mu,1,1)
      pin(mu)=p(mu,3,1)
      pb(mu)=p(mu,2,1)
      p(mu,3,1)=p4(mu)
      enddo

      pipb=dotrr(pin,pb)
      papb=dotrr(pa,pb)
      pipa=dotrr(pin,pa)
      xiab= one- (pipa+pipb)/papb!(papb-pipa-pipb)/papb

      do mu=0,3
      ptildeai(mu)=xiab*pa(mu)
      K(mu)=pa(mu)+pb(mu)-pin(mu)
      Kt(mu)=ptildeai(mu)+pb(mu)
      enddo

      Call boostkinematics(k,kt,p(0,1,1),v(0,1,1),p(0,1,8),v(0,1,8))

      do mu=0,3
      p(mu,3,1)=p3(mu)
      enddo

c   generate 7. output
      do mu=0,3
      p(mu,1,8)=ptildeai(mu)       !this is the intermediate gluon, dipole (13,2)
      p(mu,2,8)=p(mu,2,1)
      enddo
      xuz(1,8)=xiab
      xuz(2,8)=0d0 !!!!!!!!!!!!!

c    test dipole
c
      if (ldebug) then
      print*, "New Dipole 7"
      print*, "Momentum conservation of real em. kinematics:"
      do mu=0,3
      print*, p(mu,1,1)+p(mu,2,1)-p(mu,3,1)-p(mu,4,1)-v(mu,1,1)-v(mu,2,1)-v(mu,3,1),"=0?"
      enddo
      print*, "Momentum conservation of Dipole kinematics:"
      do mu=0,3
      print*, p(mu,1,8)+p(mu,2,8)-p(mu,3,8)-v(mu,1,8)-v(mu,2,8)-v(mu,3,8),"=0?"
      enddo
      print*, "------------------------------------------------"
      endif

c   start new dipole
c
      do mu=0,3
      pa(mu)=p(mu,2,1)
      pin(mu)=p(mu,3,1)
      pb(mu)=p(mu,1,1)
      p(mu,3,1)=p4(mu)
      enddo

      pipb=dotrr(pin,pb)
      papb=dotrr(pa,pb)
      pipa=dotrr(pin,pa)
      xiab= one- (pipa+pipb)/papb!(papb-pipa-pipb)/(papb)

      do mu=0,3
      ptildeai(mu)=xiab*pa(mu)
      K(mu)=pa(mu)+pb(mu)-pin(mu)
      Kt(mu)=ptildeai(mu)+pb(mu)
      enddo

      Call boostkinematics(k,kt,p(0,1,1),v(0,1,1),p(0,1,9),v(0,1,9))

      do mu=0,3
      p(mu,3,1)=p3(mu)
      enddo

c   generate 8. output
      do mu=0,3
      p(mu,1,9)=ptildeai(mu)       !this is the intermediate gluon, dipole (23,1)
      p(mu,2,9)=p(mu,1,1)
      enddo
      xuz(1,9)=xiab
      xuz(2,9)=0d0 !!!!!!!!!!!!!!!!

c    test dipole
c
      if (ldebug) then
      print*, "New Dipole 8"
      print*, "Momentum conservation of real em. kinematics:"
      do mu=0,3
      print*, p(mu,1,1)+p(mu,2,1)-p(mu,3,1)-p(mu,4,1)-v(mu,1,1)-v(mu,2,1)-v(mu,3,1),"=0?"
      enddo
      print*, "Momentum conservation of Dipole kinematics:"
      do mu=0,3
      print*, p(mu,1,9)+p(mu,2,9)-p(mu,3,9)-v(mu,1,9)-v(mu,2,9)-v(mu,3,9),"=0?"
      enddo
      print*, "------------------------------------------------"
      endif

c Dipoles 5.40

      do mu=0,3
      pa(mu) =p(mu,1,1)
      pin(mu)=p(mu,3,1)
      pj(mu) =p(mu,4,1)
      enddo

      papj=dotrr(pa,pj)
      papi=dotrr(pa,pin)
      pipj=dotrr(pin,pj)

      xija=one - pipj/(papi+papj)! (papi+papj-pipj)/(papi+papj)

      do mu=0,3
      ptildea(mu)=xija*pa(mu)
      ptildeij(mu) =pin(mu)+pj(mu)-(one-xija)*pa(mu)
      enddo
      zi=papi/(papi+papj)

c   generate 9. output
c
      do mu=0,3
      p(mu,1,10)=ptildea(mu)       !this is the intermediate gluon, dipole (13,4)
      p(mu,2,10)=p(mu,2,1)
      p(mu,3,10)=ptildeij(mu)
      do ii=1,nv
        v(mu,ii,10)=v(mu,ii,1)
      enddo
      enddo
      xuz(1,10)=xija
      xuz(2,10)=zi

      if (ldebug) then
      print*, "New Dipole 9"
      print*, "Momentum conservation of real em. kinematics:"
      do mu=0,3
      print*, p(mu,1,1)+p(mu,2,1)-p(mu,3,1)-p(mu,4,1)-v(mu,1,1)-v(mu,2,1)-v(mu,3,1),"=0?"
      enddo
      print*, "Momentum conservation of Dipole kinematics:"
      do mu=0,3
      print*, p(mu,1,10)+p(mu,2,10)-p(mu,3,10)-v(mu,1,10)-v(mu,2,10)-v(mu,3,10),"=0?"
      enddo
      print*, "------------------------------------------------"
      endif

      do mu=0,3
      pa(mu) =p(mu,2,1)
      pin(mu)=p(mu,3,1)
      pj(mu) =p(mu,4,1)
      enddo

      papj=dotrr(pa,pj)
      papi=dotrr(pa,pin)
      pipj=dotrr(pin,pj)

      xija=one - pipj/(papi+papj)! (papi+papj-pipj)/(papi+papj)

      do mu=0,3
      ptildea(mu)=xija*pa(mu)
      ptildeij(mu) =pin(mu)+pj(mu)-(one-xija)*pa(mu)
      enddo
      zi=papi/(papi+papj)

c   generate 10. output
c
      do mu=0,3
      p(mu,1,11)=p(mu,1,1)       !this is the intermediate gluon, dipole (13,4)
      p(mu,2,11)=ptildea(mu)
      p(mu,3,11)=ptildeij(mu)
      do ii=1,nv
        v(mu,ii,11)=v(mu,ii,1)
      enddo
      enddo
      xuz(1,11)=xija
      xuz(2,11)=zi

      if (ldebug) then
      print*, "New Dipole 10"
      print*, "Momentum conservation of real em. kinematics:"
      do mu=0,3
      print*, p(mu,1,1)+p(mu,2,1)-p(mu,3,1)-p(mu,4,1)-v(mu,1,1)-v(mu,2,1)-v(mu,3,1),"=0?"
      enddo
      print*, "Momentum conservation of Dipole kinematics:"
      do mu=0,3
      print*, p(mu,1,11)+p(mu,2,11)-p(mu,3,11)-v(mu,1,11)-v(mu,2,11)-v(mu,3,11),"=0?"
      enddo
      print*, "------------------------------------------------"
      print*, "Elw input momenta invariant?"
      do mu=0,3
      do ii=1,nv
      print*, vstore(mu,ii)-v(mu,1,ii)
      enddo
      enddo
      STOP
      endif

      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine create_ikj(i,k,j)
      implicit none
#include "VBFNLO/utilities/global.inc"
      integer i(2:max_kin),k(2:max_kin),j(2:max_kin)
      integer a,b,c(4)


c     map from CS notation to CE notation
c
      a = 1                     !initial state particle
      b = 2                     !initial state particle
      c(1) = 3
      c(2) = 4
      c(3) = 5
      c(4) = 6
c     maps to momentum configuration id
c     Please take note of this mapping! It's hard coded.
c     "upper-line dipoles"
      i(2) = c(1)
      k(2) = c(4)
      j(2) = c(3)
c
      i(3) = c(1)
      k(3) = c(3)
      j(3) = c(4)
c
      i(4) = c(3)
      k(4) = c(4)
      j(4) = c(1)
c
      i(5) = c(1)
      k(5) = c(4)
      j(5) = a
c
      i(6) = c(1)
      k(6) = c(3)
      j(6) = a
c      
      i(7) = c(3)
      k(7) = c(4)
      j(7) = a
c
      i(8) = a
      k(8) = c(4)
      j(8) = c(1)
c
      i(9) = a
      k(9) = c(3)
      j(9) = c(1)
c
      i(10) = a
      k(10) = c(4)
      j(10) = c(3)
c
      i(11) = a
      k(11) = c(3)
      j(11) = c(4)
c
c     "lower-line dipoles"
      i(12) = c(2)
      k(12) = c(4)
      j(12) = c(3)
c
      i(13) = c(2)
      k(13) = c(3)
      j(13) = c(4)
c
      i(14) = c(3)
      k(14) = c(4)
      j(14) = c(2)
c
      i(15) = c(2)
      k(15) = c(4)
      j(15) = b
c
      i(16) = c(2)
      k(16) = c(3)
      j(16) = b
c      
      i(17) = c(3)
      k(17) = c(4)
      j(17) = b
c
      i(18) = b
      k(18) = c(4)
      j(18) = c(2)
c
      i(19) = b
      k(19) = c(3)
      j(19) = c(2)
c
      i(20) = b
      k(20) = c(4)
      j(20) = c(3)
c
      i(21) = b
      k(21) = c(3)
      j(21) = c(4)
c     new dipoles for gg->q qb q qb
c     
      i(22) = b
      k(22) = c(2)
      j(22) = c(4)
c
      i(23) = a
      k(23) = c(1)
      j(23) = c(3)
c
      i(24) = a
      k(24) = c(1)
      j(24) = c(4)

      i(25) = b
      k(25) = c(2)
      j(25) = c(3)
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c****************************************************************************


c Terrance Figy <terrance@physics.wisc.edu>
c modified last on 17 April 2005
c Takes p1+p2 -> p3+p4 + ...px + .. qp(np+nv)
c to q1+q2 -> q2+q3 + ...qx + .. q(np+nv-1)
c This subroutine transforms p(0:3,np+nv) to q(0:3,np+nv)
c q(mu,np+nv) remains ordered after removal of emitted parton kp    
c
c 
      subroutine transform(nd,pmom,lmom,xuz,ip,kp,jp)
      implicit none
c     do interface with Dieter notation here
c     input: pmom(mu,np,1) parton momenta
c            lmom(mu,nv,1) decay momenta
c     pi is the emitter, pk is the emitted, and pj is the spectator
c     output:pmom(mu,np,nd)
c            lmom(mu,nv,nd)      
c            xuz(2,2:max_kin)
c     local variables:
c     p(mu,1:np+nv) This is filled by pmom and lmom.
c     q(mu,1:np+nv) Returned
c     xu(1) is x or y 
c     xu(2) is u or z
c     Declare variable used in subroutine
c
#include "VBFNLO/utilities/global.inc"
c     
c     you better put mxpt in global.inc
c     np is the number of partons and nv is the number of decay leptons
      integer ip,kp,jp,mu,ic,ipart,nd,j
      real*8 pmom(0:3,max_p,max_kin),lmom(0:3,max_v,max_kin),xuz(2,2:max_kin)
      real*8 p(0:3,max_p+max_v),q(0:3,max_p+max_v),xu(2)
      real*8 omy,y,x,omx,z,u
      real*8 sik,sjk,sij

      real*8 dotrr
      external dotrr

      real*8 k(0:3),kt(0:3),ks(0:3)
      real*8 kDk,ksDks,pDks,pDk
      real*8 one,two            !these can be put into global
      parameter(one=1.0d0,two=2.0d0)
      logical ldebug,ldebug1
      parameter(ldebug=.false.,ldebug1=.false.)
c
c
      if(ldebug1) print*,'i,k,j=',ip,kp,jp
c     fill p(mu,1:np+nv)
      do ic =1,n_p
         do mu=0,3
            p(mu,ic) = pmom(mu,ic,1)
         enddo
      enddo
      do ic =1,n_v
         do mu=0,3
            p(mu,ic+n_p) = lmom(mu,ic,1)
         enddo
      enddo
c
      if(ldebug1) then
         
         do j=1,n_p+n_v
            print*,'p',j,'=',(p(mu,j),mu=0,3)
         enddo
         if(ldebug) then
            do j=1,n_p
               print*,'pmom(',j,',1)=',(pmom(mu,j,1),mu=0,3)
            enddo
            do j=1,n_v
               print*,'lmom(',j,',1)=',(lmom(mu,j,1),mu=0,3)
            enddo
         endif
      endif
c
      do ic=1,n_p+n_v
         do mu=0,3
            q(mu,ic)=0d0
         enddo
      enddo
c
c      if(ldebug) print*,'before selection'
c      if(ldebug) print*,'i,k,j=',ip,kp,jp
      if(ldebug) call checkmomcon(p)
      
c     dot products maybe write a function dot for this
c      sik = p(0,ip)*p(0,kp) -p(1,ip)*p(1,kp)-
c     $     p(2,ip)*p(2,kp)-p(3,ip)*p(3,kp)
c      sjk = p(0,jp)*p(0,kp) -p(1,jp)*p(1,kp)-
c     $     p(2,jp)*p(2,kp)-p(3,jp)*p(3,kp)
c      sij = p(0,ip)*p(0,jp) -p(1,ip)*p(1,jp)-
c     $     p(2,ip)*p(2,jp)-p(3,ip)*p(3,jp)
c
      sik = dotrr(p(0,ip),p(0,kp))
      sjk = dotrr(p(0,jp),p(0,kp))
      sij = dotrr(p(0,ip),p(0,jp))
c
      if(ldebug) print*,'sij,sjk,sij=',sij,sjk,sij
c
c     Final-final
      
      if(ip.gt.2 .and. jp.gt.2) then
         if(ldebug1) print*,'final-final'
        y = sik/(sik+sjk+sij)
        omy = one - y
        z = sij/(sij+sjk)
c
        xu(1) = y
        xu(2) = z
c
        if(ldebug) then
           if(y.lt.0d0) print*,'y in ff is lt 0'
           if(z.lt.0d0) print*,'z in ff is lt 0'
        endif
c    
        if(ldebug) print*,'x=,',xu
        ipart=1
        do ic=1,n_p+n_v
           do mu =0,3
              if(ic.eq.jp) then
                 q(mu,ipart) = p(mu,jp)/omy
              elseif(ic.eq.ip) then
                 q(mu,ipart) = p(mu,ip) + p(mu,kp)-y*p(mu,jp)/omy
              elseif(ic.eq.kp) then
                 goto 20
              else
                 q(mu,ipart) = p(mu,ic)
              endif
           enddo
           ipart = ipart +1
 20        continue
        enddo
         if(ldebug) call checkmomcon(q)
         goto 60
c
c     Initial-final
      elseif(((ip.le.2) .and. (jp.gt.2))) then
         if(ldebug1) print*,'initial-final'
         omx = sjk/(sij+sik)
         x = one - omx
         u = sik/(sik+sij)
c         
         xu(1) = x           
         xu(2) = u              !1-z for fi case       
c
         if(ldebug) then
            if(x.lt.0d0) then
               print*,'x in if is lt 0'
               print*,'x=,',x
               read(*,*)
            endif
            if(u.lt.0d0) then
               print*,'u in if is lt 0'
               print*,'u=,',u
               read(*,*)
            endif
         endif
c
         if(ldebug) print*,'x=,',xu
         ipart = 1
         do ic=1,n_p+n_v
            do mu =0,3
               if(ic.eq.ip) then
                  q(mu,ipart) = x*p(mu,ip)
               elseif(ic.eq.jp) then
                  q(mu,ipart) = p(mu,jp) + p(mu,kp)-omx*p(mu,ip)
               elseif(ic.eq.kp) then
                  goto 30
               else
                  q(mu,ipart) = p(mu,ic)
               endif
            enddo
            ipart=ipart+1
 30         continue
         enddo
         if(ldebug) call checkmomcon(q)
         goto 60
c
c     Final - initial 
c     
      elseif(((ip.gt.2) .and. (jp.le.2))) then
         if(ldebug1) print*,'final-initial'
         omx = sik/(sij+sjk)
         x = one - omx
         z = sij/(sjk+sij)
         
         xu(1) = x           !=x for fi
         xu(2) = z           !=z for fi
c
         if(ldebug) then
            if(x.lt.0d0) then
               print*,'x in fi is lt 0'
               print*,'x=,',x
               read(*,*)
            endif
            if(z.lt.0d0) then
               print*,'z in fi is lt 0'
               print*,'z=,',z
               read(*,*)
            endif
         endif
c
         if(ldebug) print*,'x=,',xu
         ipart = 1
         do ic=1,n_p+n_v
            do mu =0,3
               if(ic.eq.jp) then
                  q(mu,ipart) = x*p(mu,jp)
               elseif(ic.eq.ip) then
                  q(mu,ipart) = p(mu,ip) + p(mu,kp)-omx*p(mu,jp)
               elseif(ic.eq.kp) then
                  goto 40
               else
                  q(mu,ipart) = p(mu,ic)
               endif
            enddo
            ipart=ipart+1
 40         continue
         enddo
         if(ldebug) call checkmomcon(q)
         goto 60       
c                  
c     Initial-initial 
      elseif((ip.le.2).and.(jp.le.2)) then
         if(ldebug1) print*,'initial-initial'
         omx = (sik+sjk)/sij
         x = one - omx
         xu(1) = x
         xu(2) = 0.0d0
         if(x.lt.0d0) print*,'x in ii is lt 0'
         if(ldebug) print*,'x=,',xu
c     
c     see 5.5 of Catani-seymour
         do mu = 0,3
            q(mu,ip) = x*p(mu,ip)
            q(mu,jp) = p(mu,jp)
            k(mu) = p(mu,ip)+p(mu,jp)-p(mu,kp)
            kt(mu) = q(mu,ip)+p(mu,jp)
            ks(mu) = k(mu)+kt(mu)
         enddo
         kDk=k(0)**2-k(1)**2-k(2)**2-k(3)**2
         ksDks=ks(0)**2-ks(1)**2-ks(2)**2-ks(3)**2
c     
         ipart = 3
         do ic = 3,n_p+n_v        !transform all final state particles
            if(ic.eq.kp) then 
               goto 50
            else
               pDks = p(0,ic)*ks(0)-p(1,ic)*ks(1) -
     $              p(2,ic)*ks(2)- p(3,ic)*ks(3)
               pDk =  p(0,ic)*k(0)-p(1,ic)*k(1) -
     $              p(2,ic)*k(2)- p(3,ic)*k(3)
               do mu =0,3
                  q(mu,ipart) = p(mu,ic) - 
     $                 two*pDks/ksDks*ks(mu)+two*pDk/kDk*kt(mu)
               enddo
               ipart = ipart + 1
            endif
 50         continue
         enddo   
         if(ldebug) call checkmomcon(q)
         goto 60
      else
         print*,'program will stop'
         stop       
      endif
c
 60   continue

      do mu=0,3
      if(q(mu,n_p+n_v).gt.0d0) then
         print*,'Last entry in q(mu,np+nv) not zero'
         print*,'program error'
         stop
      endif
      enddo
c     last entry in q is zero vector 
c     The higgs decay is the last nv entries in q. 
c     The first np entries are the partons.
c     fill pmom and lmom
      
      do j = 1,n_p+n_v-1
         do mu =0,3
            if(j.le.n_p-1) then
               pmom(mu,j,nd) = q(mu,j)
            else
               lmom(mu,j-n_p+1,nd) = q(mu,j)
            endif
         enddo
      enddo
      do mu=0,3
         pmom(mu,n_p,nd) = q(mu,n_p+n_v) !last entry is zero vector 
      enddo

      if(ldebug1) then
         
         do j=1,n_p+n_v           
            print*,'q(',j,')=',(q(mu,j),mu=0,3)
         enddo
c     print*,'q(',mu,np+nv,')=',q(mu,np+nv)
         if(ldebug) then
            do j=1,n_p
               print*,'pmom(',j,nd,')=',(pmom(mu,j,nd),mu=0,3)
            enddo
            do j=1,n_v
               print*,'lmom(',j,nd,')=',(lmom(mu,j,nd),mu=0,3)
            enddo
         endif
      endif
c     fill xuz
      xuz(1,nd) = xu(1)
      xuz(2,nd) = xu(2)
      if(ldebug1) then
         print*,'xuz(1,',nd,')',xuz(1,nd)
         print*,'xuz(2,',nd,')',xuz(2,nd)
         print*,'-----------------------'
c         read(*,*)
      endif
c
      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                
c     checkmomcon(p) this subr checks that momentum is conserved.
c                              and that partons are massless
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine checkmomcon(p)
      implicit none
#include "VBFNLO/utilities/global.inc"
      real*8 p(0:3,max_p+max_v),ptot(0:3)
      integer mu,j

c      First 2 in list are initial.

      do mu = 0,3
         ptot(mu) = p(mu,1)+p(mu,2)
      enddo
c     
      print*,'pa + pb =',ptot
      do j = 3,8
         do mu = 0,3
            ptot(mu) = ptot(mu) - p(mu,j)
         enddo
      enddo
      print*,'pa + pb - sum(p_i,i=3..8)=',(ptot(mu),mu=0,3)
      
c      print*,'tot mom=',(p(mu,1)+p(mu,2)-
c     $     p(mu,3)-p(mu,4)-p(mu,5)-p(mu,6)-
c     $     p(mu,7)-p(mu,8),mu=0,3)
      
c     check that partons are massless
      do j=1,n_p+n_v
         print*,'psq(',j,')',p(0,j)**2-p(1,j)**2-p(2,j)**2-p(3,j)**2
      enddo

      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

*************************************************************************
      SUBROUTINE ptilde_qed(p, xuz, v)
c*************************************************************************
c     This is the manager of all ptilde_qed subroutines and the (only) 
c     interface of those to the main program.  Depending on the process ID
c     the right ptilde_qed routine is called.
c
c*************************************************************************
c     INPUT
c     p(mu,1:np,1) : np parton momenta for real emission graph
c     v(mu,1:nv,1) : nv momenta of other particles for real emission graph
c     
c     OUTPUT
c     p(mu,1:np,j) : np parton momenta in the ptilde kinematic
c     v(mu,1:nv,j) : nv momenta of other particles in the ptilde kinematic
c                    (j = 2, ..., n_kin)
c     xuz          : These are the z, u, x_ij,a of the Catani-Seymour algorithm
c*************************************************************************

      implicit none

#include "global.inc"
#include "process.inc"

      real*8 p(0:3,max_p,max_kin)
      real*8 v(0:3,max_v,max_kin)
      real*8 xuz(2,2:max_kin)

      SELECT CASE(procID)

c VBF: Hjj, Vjj, VVjj
      CASE(Hjj, Hjj_AA, Hjj_mu, Hjj_tau, Hjj_bbar, Hjj_WW, Hjj_ZZ_ll,
     1     Hjj_ZZ_lnu, Ajj, Zjj_l, Zjj_nu, WPjj, WMjj, WPWMjj, ZZjj_ll,
     2     ZZjj_lnu, WPZjj, WMZjj, WPWPjj, WMWMjj, WPAjj, WMAjj,ZAJJ,ZAjj_n,
     3     HAjj, HAjj_AA, HAjj_mu, HAjj_tau, 
     4     HAjj_bbar, HAjj_WW, HAjj_ZZ_ll, HAjj_ZZ_lnu,
     5     Hjj_ZZhad, Hjj_WPhadWM, Hjj_WPWMhad)
          call VBFptildeqedmake(n_qed,p(0,1,1),xuz,v(0,1,1))

      CASE DEFAULT
         print *,"Error: ptilde_qed routine for process "//
     &            "not implemented, ID = ",procID
         STOP
      END SELECT

      return
      end


** The following subroutines are needed for real photon corrections.  
**        subroutines written by Terrance Figy.
c
c--------- ptilde_qed(p,xuz) ----------------------------------
c
      subroutine VBFptildeqedmake(N,p,xuz,v)
      implicit none
c
c assign ptilde momenta in notation of Catani&Seymour
c
c  input:  p(mu,1:np,1)      np parton momenta for real emission graph
c  output:                   np-1 parton momenta for Born amplitude
c          p(mu,1:np-1,2)        index=2: dipole subtraction
c          p(mu,1:np-1,3)        index=3: dipole subtraction
c                : 
c                :
c                :
c          p(mu,1:np-1,maxdips)
c
c          xuz(1,index)      x_ija etc
c          xuz(2,index)      u or z
c
c at the same time, fill the legoplot variables for the external momenta
c
c       Terrance Figy 
c       Initial version:  April 2008
c       Last modified: April 2008
c 
      integer mu,i,j,N,ii
#include "global.inc"
      real*8 p(0:3,max_p,max_kin),xuz(1:2,2:max_kin),v(0:3,max_v,max_kin)
c
c fill legoplot variables for final state particles and higgs
      integer N1,npart,ic
      real*8 tpi
*      parameter (N1=np-nx,tpi=2d0*3.141592653589793)
      parameter (N1=max_p-nx,tpi=2d0*Pi)
      REAL*8 Y,PT,PHI,PNJ
      COMMON /BLIPSO/ Y(N1), PT(N1), PHI(N1), PNJ(0:4)   ! LIPSN output 
      real*8 yp, pTp, phip, yv, pTv, phiv, yx, pTx, phix,
     $      yh,pTh,phih
      common /blegof/ yp(N1,max_kin), pTp(N1,max_kin), phip(N1,max_kin),
     1                yv(max_v,max_kin), pTv(max_v,max_kin), 
     2                phiv(max_v,max_kin),yh(max_kin),pTh(max_kin),
     3                phih(max_kin)
      logical ldebug
      parameter (ldebug=.false.)
      real*8 ph(0:3)            ! 4-momentum of higgs
      double precision mompart(0:3,1:8),momphtn(0:3)
c     this is the mapping to the dipole subtraction terms
c
      integer ip(2:13),jp(2:13)
      data ip /3,4,3,4,3,4,1,1,2,2,1,2/ ! emitter
      data jp /4,3,1,1,2,2,3,4,3,4,2,1/ ! spectator


*      write(*,*)'ptildeQED p(0,5,1) =', p(0,5,1)
*      write(*,*)'ptildeQED p(1,5,1) =', p(1,5,1)
*      write(*,*)'ptildeQED p(2,5,1) =', p(2,5,1)
*      write(*,*)'ptildeQED p(3,5,1) =', p(3,5,1)
c
c  fill legoplot variables for final state particles
      if (N.eq.1) then
         npart = N1-1
      else
         npart = N1
      endif
c  set superfluous variables to zero
      do ic = npart+1,N1
         do ii = 1,N
            yp(ic,ii) = 0
            pTp(ic,ii) = 0
            phip(ic,ii) = 0
            do mu = 0,3
               p(mu,nx+ic,ii) = 0
            enddo
         enddo
      enddo
      do ic = 1,npart
         yp(ic,1) = y(ic)
         pTp(ic,1) = PT(ic)
         phip(ic,1) = phi(ic)
         do ii = 2,N
            yp(ic,ii) = 0
            pTp(ic,ii) = 0
            phip(ic,ii) = 0
         enddo
         if (ldebug) then
            write(6,*) ' '
            call legoy(p(0,nx+ic,1),ptx,yx,phix)
c            print*,' parton ',ic,', 1:',y(ic)/yx,pt(ic)/ptx,phi(ic)/phix
         endif
      enddo
      do ic = 1,max_v
         call legoy(v(0,ic,1),ptv(ic,1),yv(ic,1),phiv(ic,1))
      enddo
      do mu = 0,3 ! compute momenta of higgs
         ph(mu) = 0.0d0
         do ic = 1,max_v
            ph(mu) = ph(mu) + v(mu,ic,1)
         enddo
      enddo
      call legoy(ph,pth(1),yh(1),phih(1)) ! fill lego plot for higgs
c            
c  fill legoplot variables for the ptilde
c  now that we have the ptilde we can loop over momentum configurations

      do ii=2,N                 ! loop over dipoles
c     compute transformed momenta
c     First fill momentum
         do mu = 0,3
            mompart(mu,1) = p(mu,1,1) ! quark
            mompart(mu,2) = p(mu,2,1) ! quark
            mompart(mu,3) = p(mu,3,1) ! quark
            mompart(mu,4) = p(mu,4,1) ! quark
            mompart(mu,5) = v(mu,1,1) ! lepton/b quark
            mompart(mu,6) = v(mu,2,1) ! lepton/b quark
            mompart(mu,7) = v(mu,3,1) ! lepton/b quark
            mompart(mu,8) = v(mu,4,1) ! lepton/b quark
            momphtn(mu)   = p(mu,5,1) ! photon
         enddo
         call get_qed_transform(mompart,momphtn,ip(ii),jp(ii),
     $        xuz(1,ii))
c     
c     Fill transformed momentum
         do mu = 0,3
            p(mu,1,ii) = mompart(mu,1) ! quark
            p(mu,2,ii) = mompart(mu,2) ! quark
            p(mu,3,ii) = mompart(mu,3) ! quark
            p(mu,4,ii) = mompart(mu,4) ! quark
            v(mu,1,ii) = mompart(mu,5) ! lepton/b quark
            v(mu,2,ii) = mompart(mu,6) ! lepton/b quark
            v(mu,3,ii) = mompart(mu,7) ! lepton/b quark
            v(mu,4,ii) = mompart(mu,8) ! lepton/b quark
         enddo
c     
c  check momentum conservation
c    
c            print*,'xuzmin=',xuzmin  
         do j=3,max_p-1         !loop over part. label j
            call lego(p(0,j,ii),ptp(j-2,ii),yp(j-2,ii),phip(j-2,ii)) 
         enddo
         do j =1,max_v
            call legoy(v(0,j,ii),ptv(j,ii),yv(j,ii),phiv(j,ii)) 
         enddo
         do mu = 0,3            ! compute momenta of higgs
            ph(mu) = 0.0d0
            do ic = 1,max_v
               ph(mu) = ph(mu) + v(mu,ic,ii)
            enddo
         enddo
         call legoy(ph,pth(ii),yh(ii),phih(ii)) ! fill lego plot for higgs
      enddo
c     Last entry is the soft/collinear photon 
      do ii=2,N
         ptp(max_p-2,ii)=0
         yp(max_p-2,ii)=0
         phip(max_p-2,ii)=0
      enddo


      if (ldebug) then
         do ii=2,N
            do ic=1,npart-1
               call lego(p(0,ic+2,ii),ptx,yx,phix)
               print*,' parton ',ic,ii,' :',yp(ic,ii)/yx,ptp(ic,ii)/ptx,
     1              phip(ic,ii)/phix
            enddo
         enddo
      endif
      if ( ldebug ) then
 10      format( ' p(', i1, ') = ', 4(f10.3, 2x) )
c     
         do j = 1,N
            if (j.eq.1) print*,' (n+1) parton momenta '
            if (j.ge.2) print*,' n-parton momenta for',
     -           j,'momentum config'
            do i = 1, max_p-min(j-1,1)
               write(6,10) i, p(0,i,j), p(1,i,j), p(2,i,j), p(3,i,j)
            end do
            print*,' pT = ',(ptp(i,j),i=1,4-min(j-1,1))
            print*,' y  = ',(yp(i,j),i=1,4-min(j-1,1))
            print*,' phi= ',(phip(i,j),i=1,4-min(j-1,1))
            print*, ' lepton pT = ',(pTv(i,j),i=1,2)
            print*, ' lepton y  = ',(yv(i,j),i=1,2)
            print*, ' lepton phi= ',(phiv(i,j),i=1,2)
         enddo         
      endif
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c  Real 
c  fermion(p(mu,1)) + fermion(p(mu,2)) -> fermion(p(mu,3)) + ... + photon(k(mu))
c  Born
c  fermion(p(mu,1)) + fermion(p(mu,2)) -> fermion(p(mu,3)) + ... 
c  i (emitter)
c  j (spectator) 
c  mompart(mu,n=1,6) parton + higgs decay products 4 vector momenta
c  momphtn(mu) photon momentum 4 vector
c  compute transformed momentum
c  dipole subtraction terms   
c 
c  Terrance Figy 
c  April 25, 2008
c  May 22, 2008
c  Reference: arXiv: hep-ph/9904440
c  "A General Approach to Photon Radiation off Fermions"
c

      subroutine get_qed_transform(mompart,momphtn,i,j,xuz)
      implicit none
      integer i,j ! emitter and spectator 
      integer ii,aa,bb,num,mu
      double precision mompart(0:3,1:8),momphtn(0:3)
      double precision k(0:3),p(0:3,1:8),ptilde(0:3,1:8)
      double precision xuz(1:2)
      double precision Pab(0:3),PTab(0:3),PplusPTab(0:3)
      double precision piDpj,pjDk,piDk,yy,zz
      double precision paDpi,paDk,xx
      double precision paDpb,pbDk
      double precision PabdotPab,PTabdotPTab,PabdotPTab 
      double precision zero,one,two
      parameter(zero=0.0d0,one=1.0d0,two=2.0d0)
      logical ldebug
      parameter (ldebug=.false.)

      double precision dotrr
      external dotrr
      
c     fill parton 4-vectors
      do num =1,8
         do mu=0,3
            p(mu,num) = mompart(mu,num) 
         enddo
      enddo
c     fill photon 4-vector
      do mu=0,3
         k(mu) = momphtn(mu)
      enddo
c

c     final-final case
      if(i.gt.2 .and. j.gt.2) then
c     dot products
         piDpj = dotrr(p(0,i),p(0,j))
         pjDk  = dotrr(p(0,j),k)
         piDk  = dotrr(p(0,i),k)  
c     x and y
         yy = piDk/(piDpj + piDk + pjDk)
         zz = piDpj/(piDpj + pjDk)
c     fill xuz
         xuz(1) = yy
         xuz(2) = zz
c
c     transformed momentum
c     eq (3.4)         

         do num=1,8
            do mu =0,3
               if(num.eq.i) then
                  ptilde(mu,num) = p(mu,i) + k(mu) - yy/(one-yy)*p(mu,j)
               elseif(num.eq.j) then 
                  ptilde(mu,num) = p(mu,j)/(one-yy)
               else
                  ptilde(mu,num) = p(mu,num)
               endif
            enddo
         enddo

         if(ldebug) then
            print*,'final-final'
            do mu=0,3
               print*,'p(',i,')+p(',j,')+k=',p(mu,i)+p(mu,j)+k(mu)
               print*,'pt(',i,')+pt(',j,')=',ptilde(mu,i)+ptilde(mu,j)
               
            enddo
*P            read(*,*)
         endif
         
c     replace partonic momenta
         do num =1,8
            do mu=0,3
               mompart(mu,num)=ptilde(mu,num)
               momphtn(mu) = zero
            enddo
         enddo
         
         return
      endif
c
c     final-initial      
c     i is a final state emitter p_i
c     j is an initial state spectator p_a
c
      if(i.gt.2 .and. j.le.2) then
c     set up labels (this should make formula checking easier)
         aa = j
         ii = i
c     dot products
         paDpi = dotrr(p(0,aa),p(0,ii))
         paDk  = dotrr(p(0,aa),k)
         piDk  = dotrr(p(0,ii),k)
c         
         xx = (paDpi + paDk - piDk)/(paDpi + paDk)
         zz = paDpi/(paDpi + paDk)
c     fill xuz
         xuz(1) = xx
         xuz(2) = zz
c
c     transformed momentum
c     eq (3.4)         
         do num=1,8
            do mu =0,3
               if(ii.eq.num) then
                  ptilde(mu,num) = p(mu,ii)+k(mu)-(one-xx)*p(mu,aa) ! pi
               elseif(aa.eq.num) then
                  ptilde(mu,num) = xx*p(mu,aa) ! pa
               else
                  ptilde(mu,num) = p(mu,num)
               endif
            enddo
         enddo
         if(ldebug) then
            print*,'final-initial'
            do mu=0,3
               print*,'p(',ii,')-p(',aa,')+k=',
     $              p(mu,ii)-p(mu,aa)+k(mu)
               print*,'pt(',ii,')-pt(',aa,')=',
     $              ptilde(mu,ii)-ptilde(mu,aa)
              
            enddo
*P            read(*,*)
         endif
c        replace partonic momenta
         do num =1,8
            do mu=0,3
               mompart(mu,num)=ptilde(mu,num)
               momphtn(mu) = zero
            enddo
         enddo
         
         return

      endif
c
c     initial-final
c     i is an initial state emitter p_a
c     j is a final state spectator  p_i
c
      if(i.le.2 .and. j.gt.2) then
c     set up labels (this should make formula checking easier)
         aa = i
         ii = j
c     dot products
         paDpi = dotrr(p(0,aa),p(0,ii))
         paDk  = dotrr(p(0,aa),k)
         piDk  = dotrr(p(0,ii),k)
c         
         xx = (paDpi + paDk - piDk)/(paDpi + paDk)
         zz = paDpi/(paDpi + paDk)
c     fill xuz
         xuz(1) = xx
         xuz(2) = zz
c
c     transformed momentum
c     eq (3.4)         
         do num=1,8
            do mu =0,3
               if(ii.eq.num) then
                  ptilde(mu,num) = p(mu,ii)+k(mu)-(one-xx)*p(mu,aa) ! pi
               elseif(aa.eq.num) then
                  ptilde(mu,num) = xx*p(mu,aa) ! pa
               else
                  ptilde(mu,num) = p(mu,num)
               endif
            enddo
         enddo
         if(ldebug) then
            print*,'initial-final'
            do mu=0,3
               print*,'p(',ii,')-p(',aa,')+k=',
     $              p(mu,ii)-p(mu,aa)+k(mu)
               print*,'pt(',ii,')-pt(',aa,')=',
     $              ptilde(mu,ii)-ptilde(mu,aa)
            enddo
*P            read(*,*)
         endif
c        replace partonic momenta
         do num =1,8
            do mu=0,3
               mompart(mu,num)=ptilde(mu,num)
               momphtn(mu) = zero
            enddo
         enddo
         
         return 
        
      endif

c     initial-initial

      if(i.le.2 .and. j.le.2) then
c     set up labels (this should make formula checking easier)
         aa = i
         bb = j
c     dot products
         paDpb = dotrr(p(0,aa),p(0,bb))
         paDk  = dotrr(p(0,aa),k)
         pbDk  = dotrr(p(0,bb),k)
c eq 3.21
         xx = (paDpb - paDk - pbDk)/paDpb
         yy = paDk/paDpb
c     fill xuz
         xuz(1) = xx
         xuz(2) = yy
c
c Pab = pa + pb - k eq 3.24
c PTab = x pa + pb eq 3.25
c
         do mu=0,3
            Pab(mu) = p(mu,aa) + p(mu,bb) - k(mu)
            PTab(mu) = xx*p(mu,aa) + p(mu,bb)
            PplusPTab(mu) = Pab(mu) + PTab(mu)
         enddo
         
         PabdotPab = dotrr(Pab,Pab)
         PTabdotPTab = dotrr(PTab,PTab)
         PabdotPTab = dotrr(Pab,PTab)
         
         do num=1,8
            do mu=0,3           !eq 3.26, 3.27
               if(num.eq.aa) then 
                  ptilde(mu,num)=xx*p(mu,aa) 
               elseif(num.eq.bb) then
                  ptilde(mu,num)=p(mu,bb) 
               else ! apply LT to all other particles
                  ptilde(mu,num) = p(mu,num) -PplusPTab(mu)*
     $                 dotrr(PplusPTab,p(0,num))/
     $                 (PabdotPab+PabdotPTab) +
     $                 two*PTab(mu)*dotrr(Pab,p(0,num))/PabdotPab
               endif
            enddo
         enddo
         if(ldebug) then
            print*,'initial-initial'
            do mu=0,3
               print*,'p(',aa,')+p(',bb,')-k-Sumjp=',
     $              p(mu,aa)+p(mu,bb)-k(mu)-p(mu,3)-
     $              p(mu,4)-p(mu,5)-p(mu,6)
               print*,'pt(',aa,')+pt(',bb,')-Sumjptilde=',
     $              ptilde(mu,aa)+ptilde(mu,bb)-ptilde(mu,3)-
     $              ptilde(mu,4)-ptilde(mu,5)-ptilde(mu,6)
            enddo
*P            read(*,*)
         endif
c        replace partonic momenta
         do num =1,8
            do mu=0,3
               mompart(mu,num)=ptilde(mu,num)
               momphtn(mu) = zero
            enddo
         enddo   
            
         return

      endif
      
      end 

c
c     computes g^{sub}_{ij,+}
c     g^{sub}_{ij,-} is zero for small fermion masses
c     We are working in the limit of small fermion masses.
c
      subroutine get_qed_dipole(mompart,momphtn,i,j,xuz,gsub)
      implicit none 
      integer i,j               ! emitter and spectator 
      integer ii,aa,bb,num,mu
      double precision mompart(0:3,1:8),momphtn(0:3)
      double precision k(0:3),p(0:3,1:8)
      double precision xuz(1:2)
      double precision gsub
      double precision piDk,yy,zz
      double precision paDk,xx
      double precision zero,one,two
      parameter(zero=0.0d0,one=1.0d0,two=2.0d0)
      double precision constii,constff,constif,constfi
      parameter(constii=0.0d0)  !for checking of dipole subtraction
      parameter(constff=0.0d0)  !for checking of dipole subtraction
      parameter(constif=0.0d0)  !for checking of dipole subtraction 
      parameter(constfi=0.0d0)  !for checking of dipole subtraction 
      logical ldebug
      parameter(ldebug=.false.)

      double precision dotrr
      external dotrr
c
c     
c      gsub1 = 0.0d0
c
c     fill parton 4-vectors
      do num =1,8
         do mu=0,3
            p(mu,num) = mompart(mu,num) 
         enddo
      enddo
c     fill photon 4-vector
      do mu=0,3
         k(mu) = momphtn(mu)
      enddo
c

c     final-final case
      if(i.gt.2 .and. j.gt.2) then
c     dot products
         piDk  = dotrr(p(0,i),k)  
c     x and y
         yy =  xuz(1)
         zz =  xuz(2)

c
c    compute g_{ij,+}(sub) gsub(1) 
c    compute g_{ij,-}(sub) gsub(2)
c    eq (3.1)     
         gsub = (two/(one-zz*(one-yy)) -one - zz)
         gsub = gsub/(one-yy)/piDk + yy*constff/piDk
c
         if(ldebug) then
            print*,'i=',i,'j=',j,'gsub=',gsub
            print*,'piDk=',piDk
         endif
         return
      endif
c
c     final-initial      
c     i is a final state emitter p_i
c     j is an initial state spectator p_a
c
      if(i.gt.2 .and. j.le.2) then
c     set up labels (this should make formula checking easier)
         aa = j
         ii = i
c     dot products
         
         piDk  = dotrr(p(0,ii),k)
c         
         xx = xuz(1)  
         zz = xuz(2)

c
c    compute g_{ij,+}(sub) gsub(1) 
c    compute g_{ij,-}(sub) gsub(2)
c    eq (3.1)     
         gsub = (two/(two-xx-zz) - one - zz)
         gsub = gsub/xx/piDk + constfi * (one - xx)/xx/piDk
c
c
         if(ldebug) then
            print*,'i=',i,'j=',j,'gsub=',gsub
            print*,'piDk=',piDk
         endif
         
         return

      endif
c
c     initial-final
c     i is an initial state emitter p_a
c     j is a final state spectator  p_i
c
      if(i.le.2 .and. j.gt.2) then
c     set up labels (this should make formula checking easier)
         aa = i
         ii = j
c     dot products
         
         paDk  = dotrr(p(0,aa),k)
        
c         
         xx = xuz(1)
         zz = xuz(2)

c
c    compute g_{ij,+}(sub) gsub(1) 
c    compute g_{ij,-}(sub) gsub(2)
c    eq (3.1)     
         gsub = (two/(two-xx-zz) - one - zz)
         gsub = gsub/xx/paDk + constif * (one - zz)/xx/paDk
c
        if(ldebug) then
            print*,'i=',i,'j=',j,'gsub=',gsub
            print*,'paDk=',paDk
         endif 
         return 
        
      endif

c     initial-initial

      if(i.le.2 .and. j.le.2) then
c     set up labels (this should make formula checking easier)
         aa = i
         bb = j
c     dot products
       
         paDk  = dotrr(p(0,aa),k)
        
c eq 3.21
         xx = xuz(1)
         yy = xuz(2) 

c eq 3.22
         gsub = (two/(one-xx) - one -xx)
         gsub = gsub/xx/paDk + constii*yy/xx/paDk
c
c
         if(ldebug) then
            print*,'i=',i,'j=',j,'gsub=',gsub
            print*,'paDk=',paDk
         endif

         return
      endif
      
      end
c
c
