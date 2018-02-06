c*************************************************************************
      SUBROUTINE InitPhaseSpace
c*************************************************************************
c     Setup phasespace generators.
c*************************************************************************

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/process.inc"
#include "VBFNLO/utilities/cuts.inc"
      integer WWJJ_interf
      common /cWWJJ_interf/ WWJJ_interf

      SELECT CASE(procID)
      CASE(Hjj, Hjj_AA, Hjj_mu, Hjj_tau, Hjj_bbar, Hjjj,
     1     Hjjj_AA, Hjjj_mu, Hjjj_tau, Hjjj_bbar, 
     2     HAjj, HAjj_AA, HAjj_mu, HAjj_tau, 
     3     HAjj_bbar)
         call def_ps_vnj(6)
      CASE(AAjj)
         call def_ps_vnj(39)
c         call def_ps_vnj(6)   !take this line instead of the previous one if 
                               !you want to run this process 
                               !with a Higgs instead of a spin-2 resonance 
                               !(second part of spin2toaa.F)
      CASE(Sp2jj_WW)
        if (with_spin2) then
         call def_ps_vvnj(39,3,4,1,1,1)
        else
         call def_ps_vvnj(6,3,4,1,1,1)  !Higgs instead of spin-2
        endif
      CASE(Sp2jj_ZZ_ll, Sp2jj_ZZ_lnu) 
        if (with_spin2) then
         call def_ps_vvnj(39,2,2,1,1,1) 
        else
         call def_ps_vvnj(6,2,2,1,1,1)  !Higgs instead of spin-2
        endif
      CASE(Hjj_WW,Hjjj_WW, HAjj_WW, Hjj_WPhadWM, Hjj_WPWMhad)
        if (higgsscheme .eq. 0 .or. higgsscheme .eq. -1) then
         call def_ps_vvnj(6,3,4,0,1,1)
        else
         call def_ps_vvnj(6,3,4,1,1,1)
        endif
      CASE(Hjj_ZZ_ll, Hjj_ZZ_lnu,Hjjj_ZZ_ll,Hjjj_ZZ_lnu,
     $        HAjj_ZZ_ll, HAjj_ZZ_lnu, Hjj_ZZhad)
        if (higgsscheme .eq. 0 .or. higgsscheme .eq. -1) then
         call def_ps_vvnj(6,2,2,0,1,1)
        else
         call def_ps_vvnj(6,2,2,1,1,1)
        endif
      CASE(Ajj)
         call def_ps_vnj(1)
      CASE(Zjj_l, Zjj_nu, QCDZjj_l, QCDZjj_nu)
         call def_ps_vnj(2)
      CASE(WPjj, QCDWPjj)
         call def_ps_vnj(3)
      CASE(WMjj, QCDWMjj)
         call def_ps_vnj(4)
      CASE(WPWMjj, WPhadWMjj, WPWMhadjj)
         call def_ps_vvnj(-6,3,4,1,1,1)
      CASE(WPWPjj,WMWMjj, WPhadWPjj, WMhadWMjj)
         call def_ps_vvnj(0,3,4,1,1,1)
      CASE(ZZjj_ll, ZZjj_lnu, ZZhadjj) 
         call def_ps_vvnj(-6,2,2,1,1,1)
      CASE(HHjj, HHjj_bbtau, HHjj_bbAA) 
         call def_ps_hhnj(0,6,6,1,0,0)
      CASE(WPZjj, WPhadZjj, WPZhadjj)
         call def_ps_vvnj(0,3,2,1,1,1)
      CASE(WMZjj, WMhadZjj, WMZhadjj)
         call def_ps_vvnj(0,4,2,1,1,1)
      CASE(WPAjj,QCDWPAjj)
         call def_ps_qqvaj(3,1)
      CASE(WMAjj,QCDWMAjj)
         call def_ps_qqvaj(4,1)
      CASE(ZAjj,ZAjj_n)
         call def_ps_qqvaj(2,1)
      CASE(GFHjj, GFHjj_AA, GFHjj_mu, GFHjj_tau, GFHjj_bbar)
         call def_ps_vnj(6) 
      CASE(GFHjj_WW)
         call def_ps_vvnj(6,3,4,0,1,1)
      CASE(GFHjj_ZZ_ll, GFHjj_ZZ_lnu)
         call def_ps_vvnj(6,2,2,0,1,1)

#ifdef WITH_DIBOSON
      CASE(WW,GFWW,WPhadWMlep,WPlepWMhad,GFWPhadWMlep,GFWPlepWMhad)
         call def_ps_qqVVj(3,4)
      CASE(WPZ,WPhadZ,WPZhad)
         call def_ps_qqVVj(3,2)
      CASE(WMZ,WMhadZ,WMZhad)
         call def_ps_qqVVj(4,2)
      CASE(ZZ,ZZhad,GFZZ,GFZZhad)
         call def_ps_qqVVj(2,2)
      CASE(WPA)
         call def_ps_qqVAj(3,1)
      CASE(WMA)
         call def_ps_qqVAj(4,1)

         ! could also use def_ps_vnj here
      CASE(WP_only)
         call def_ps_qqVj(3) 
      CASE(WM_only)
         call def_ps_qqVj(4)

      CASE( WPH, WPH_AA, WPH_mu, WPH_tau, WPH_bbar)
         call def_ps_qqVHj(3,6,-1,-1)  !-1 here means no H->VV
      CASE(WPH_WW)   
         call def_ps_qqVHj(3,6,3,4) ! H->WW
      CASE(WPH_ZZ_ll,  WPH_ZZ_lnu)
         call def_ps_qqVHj(3,6,2,2) ! W H->W ZZ

      CASE( WMH, WMH_AA, WMH_mu, WMH_tau, WMH_bbar)
         call def_ps_qqVHj(4,6,-1,-1)  !-1 here means no H->VV
      CASE(WMH_WW)   
         call def_ps_qqVHj(4,6,3,4) ! H->WW
      CASE(WMH_ZZ_ll,  WMH_ZZ_lnu)
         call def_ps_qqVHj(4,6,2,2) ! W H->W ZZ

      CASE(ZA,GFZA)
         call def_ps_qqVAj(2,1)
      CASE(AA,GFAA)
         call def_ps_qqAAj(1,1)
#endif

#ifdef WITH_DIBOSONJET
      CASE(WWJ,WPHADWMJ,WPWMHADJ,GFWWj,GFWPHADWMJ,GFWPWMHADJ)
         call def_ps_qqVVj(3,4)
      CASE(WMAJ,WPAJ)
         call def_ps_vnj(3)
      CASE(WMZJ, WMhadZJ, WMZhadJ)
         call def_ps_vvnj(0,4,2,1,1,1)
      CASE(WPZJ, WPhadZJ, WPZhadJ)
         call def_ps_vvnj(0,3,2,1,1,1)
!           call def_ps_qqVVj(3,2)         
      CASE(ZZJ,GFZZj)
         call def_ps_qqVVj(2,2)

         ! could also use def_ps_vnj here
      CASE(WPJ)
         call def_ps_qqVj(3) 
      CASE(WMJ)
         call def_ps_qqVj(4)

      CASE(WPHJ, WPHJ_AA, WPHJ_mu, WPHJ_tau, WPHJ_bbar )
         call def_ps_qqVHj(3,6,-1,-1)  !-1 here means no H->VV
      CASE(WPHJ_WW)   
         call def_ps_qqVHj(3,6,3,4) ! H->WW
      CASE(WPHJ_ZZ_ll,  WPHJ_ZZ_lnu)
         call def_ps_qqVHj(3,6,2,2) ! W H->W ZZ

      CASE(WMHJ, WMHJ_AA, WMHJ_mu, WMHJ_tau, WMHJ_bbar)
         call def_ps_qqVHj(4,6,-1,-1)  !-1 here means no H->VV
      CASE(WMHJ_WW)   
         call def_ps_qqVHj(4,6,3,4) ! H->WW
      CASE(WMHJ_ZZ_ll,  WMHJ_ZZ_lnu)
         call def_ps_qqVHj(4,6,2,2) ! W H->W ZZ
#endif

#ifdef WITH_QCDVV
      CASE(QCDWPZjj)
         call def_ps_qqVVj(3,2)
      CASE(QCDWMZjj)
         call def_ps_qqVVj(4,2)
      CASE(QCDAAjj)
         call def_ps_qqAAj(1,1)
      CASE(QCDZZjj_ll,QCDZZjj_lnu)
         call def_ps_qqVVj(2,2)
      CASE(QCDZAjj_l,QCDZAjj_n)
c         if(mza_min.ge.250d0) PS_loops=1
         call def_ps_qqvaj(2,1)
      CASE(QCDWPWPjj)
!       if(sub_number.eq.2.and.(WWJJ_interf.ge.2)) then ! include Higgs resonance
!          call def_ps_WWW(3,4,3)
!       else
         call def_ps_qqVVj(3,3)
!       endif
      CASE(QCDWMWMjj)
         call def_ps_qqvvj(4,4)
#endif

#ifdef WITH_TRIBOSON

        CASE(WWZ, WPhadWMZ, WPWMhadZ, WWZhad)
#ifdef WITH_KK
        if (with_kk) then
           call def_ps_WWZ_kk(3,2,4)
        else
#endif
           call def_ps_qqVVVj(3,4,2)
#ifdef WITH_KK
        endif
#endif
        CASE(ZZWP, ZZhadWP, ZZWPhad)
#ifdef WITH_KK
        if (with_kk) then
           call def_ps_ZZW_kk(2,3,2)
        else
#endif
           call def_ps_qqVVVj(2,2,3)
#ifdef WITH_KK
        endif   
#endif
      CASE(ZZWM, ZZhadWM, ZZWMhad)
#ifdef WITH_KK
        if (with_kk) then
         call def_ps_ZZW_kk(2,4,2)
        else
#endif
           call def_ps_qqVVVj(2,2,4)
#ifdef WITH_KK
        endif
#endif
      CASE(WWWP, WPhadWMWP, WPWMhadWP)
#ifdef WITH_KK
        if (with_kk) then
         call def_ps_WWW_kk(3,4,3)
        else
#endif
         call def_ps_WWW(3,4,3)
#ifdef WITH_KK
        endif
#endif
      CASE(WWWM, WMhadWPWM, WMWPhadWM)
#ifdef WITH_KK
        if (with_kk) then
         call def_ps_WWW_kk(4,3,4)
        else
#endif
         call def_ps_WWW(4,3,4)
#ifdef WITH_KK
        endif
#endif
      CASE(ZZZ, ZZZhad)
         call def_ps_ZZZj(2,2,2)
      CASE(WWA, WPhadWMA, WPWMhadA)
         call def_ps_qqVVAj(3,4,1)
      CASE(ZZA, ZZhadA, ZZnA)
         call def_ps_qqVVAj(2,2,1)
      CASE(WPZA, WPhadZA, WPZhadA)
         call def_ps_qqVVAj(3,2,1)
      CASE(WMZA, WMhadZA, WMZhadA)
         call def_ps_qqVVAj(4,2,1)
      CASE(WPAA)
         call def_ps_qqVAAj(3,1,1)
      CASE(WMAA)
         call def_ps_qqVAAj(4,1,1)
      CASE(ZAA,ZnAA)
         call def_ps_qqVAAj(2,1,1)
      CASE(AAA)
         call def_ps_qqAAAj(1,1,1)

#endif

#ifdef WITH_TRIBOSONJET
      CASE(WPAAj)
         call def_ps_qqVAAj(3,1,1)
      CASE(WMAAj)
         call def_ps_qqVAAj(4,1,1)
#endif

      CASE DEFAULT
         print*,"No initialization of phase space for this process"
         stop    
      END SELECT

      RETURN
      END


c*************************************************************************
      SUBROUTINE phasespace(rand,
     &                      p, x,
     &                      v, 
     $                      ps_number,  
     &                      dps  )
c*************************************************************************
c     This is the manager of all phasespace generators and the (only) 
c     interface of those to the main program.  Depending on the process ID
c     the right phasespace generator is called with the correct calling 
c     convention, which might well be individually different.
c
c*************************************************************************
c     INPUT
c     rand            : Array with random numbers driven by monaco. 
c     p               : 4-momenta of the partons involved in the basic
c                       QCD 2->2(+1)+X process.
c                       The first 2 entries are the incoming partons.
c     x               : the feynman x of the incoming partons 
c     v               : The 4-momenta of the other particles
c     ps_number       : Number specifying the desired phase space
c     
c     OUTPUT
c     dps             : The jacobi factor of the phasespace point.
c*************************************************************************
          use globalvars, only: ldoblha
      implicit none
      
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/process.inc"
#include "VBFNLO/utilities/cuts.inc"

      real*8 x(nx) 
      real*8 rand(max_PS_dim)
      real*8 p(0:3,max_p,max_kin)
      real*8 v(0:3,max_v,max_kin)
      real*8 dps

      real*8 rn
      real*8 Q(0:4)
      integer i,L,mu, ps_number
      common/MB1/ rn
      real*8 RandomNumber, mjj2
      external RandomNumber, mjj2
      double precision phirn, cosrn, sinrn, tmpx, tmpy
      real*8 TPI
      parameter (TPI=2D0*PI )

      integer nsing/0/

      integer WWJJ_interf
      common /cWWJJ_interf/ WWJJ_interf

      DOUBLE PRECISION CLR,XM2,XMG,CB,CV,CA
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),CB(6,6,6),CV(4,5),CA(4,5)

c fill momentum arrays with zero
      do L=1,max_kin
         do i = 1,max_p
            do mu=0,3
               p(mu,i,L) = 0d0
            enddo
         enddo
      enddo
      do L=1,max_kin
         do i=1,max_v
            do mu=0,3
               v(mu,i,L) = 0d0
            enddo
         enddo
      enddo
      p=0d0
      v=0d0

      rn = 0d0
c switch between different processes

      SELECT CASE(procID)
      CASE(Hjj, Hjj_AA, Hjj_mu, Hjj_tau, Hjj_bbar) ! pp-> Hjj(j) -> xxjj(j)
        if (higgsscheme .eq. 0 .or. higgsscheme .eq. -1) then
         call LIPSN(n_p-nx, rand, rn,
     1              p(0,1,1), p(0,2,1), 
     2              Q, v(0,1,1), v(0,2,1),
     3              p(0,3,1), X(1), X(2), dps, 1)
        else
         call LIPSN(n_p-nx, rand, rn,
     1              p(0,1,1), p(0,2,1), 
     2              Q, v(0,1,1), v(0,2,1),
     3              p(0,3,1), X(1), X(2), dps, 0)
        endif
        if (ldoblha .and. (procID.eq.Hjj)) then
          do i=0,3
            v(i,1,1) = v(i,1,1)+v(i,2,1)
            v(i,2,1) = 0d0
          enddo
          dps = dps / (xmg(6)/8.0d0/TPI)
        endif
      CASE(AAjj) ! pp-> spin-2 jj(j) -> AAjj(j)
         call LIPSN(n_p-nx, rand, rn,
     1              p(0,1,1), p(0,2,1), 
     2              Q, v(0,1,1), v(0,2,1),
     3              p(0,3,1), X(1), X(2), dps, 1)
      CASE(Sp2jj_WW, Sp2jj_ZZ_ll, Sp2jj_ZZ_lnu) ! pp -> spin-2 jj -> VVjj
         call LIPSNVV(n_p-nx,rand,rn,
     1                p(0,1,1),p(0,2,1),
     2                v(0,1,1),p(0,3,1),
     3                X(1),X(2),dps)
      CASE(Hjj_WW, Hjj_ZZ_ll, Hjj_ZZ_lnu,
     1     Hjj_WPhadWM, Hjj_WPWMhad, Hjj_ZZhad) ! pp -> Hjj(j) -> VVjj(j)
         call LIPSNVV(n_p-nx,rand,rn,
     1                p(0,1,1),p(0,2,1),
     2                v(0,1,1),p(0,3,1),
     3                X(1),X(2),dps)
      CASE(Ajj)
         call LIPSN0(n_p-nx, rand, rn,
     1               p(0,1,1), p(0,2,1), 
     2               Q, p(0,3,1), X(1), X(2), dps)
         v(0,1,1) = Q(0)
         v(1,1,1) = Q(1)
         v(2,1,1) = Q(2)
         v(3,1,1) = Q(3)   
      CASE(HAjj, HAjj_AA, HAjj_mu, HAjj_tau,  ! pp-> HAjj(j) -> Axxjj(j)
     3     HAjj_bbar)
         call LIPSN(n_p+1-nx, rand, rn,
     1              p(0,1,1), p(0,2,1), 
     2              Q, v(0,1,1), v(0,2,1),
     3              p(0,3,1), X(1), X(2), dps,1)
         v(0,3,1)=p(0,n_p+1,1)
         v(1,3,1)=p(1,n_p+1,1)
         v(2,3,1)=p(2,n_p+1,1)
         v(3,3,1)=p(3,n_p+1,1)
         p(0,n_p+1,1)=0d0
         p(1,n_p+1,1)=0d0
         p(2,n_p+1,1)=0d0
         p(3,n_p+1,1)=0d0
         if (ldoblha .and. (procID.eq.HAjj)) then
           do i=0,3
             v(i,1,1) = v(i,1,1)+v(i,2,1)
             v(i,2,1) = v(i,3,1)
             v(i,3,1) = 0d0
           enddo
           dps = dps / (xmg(6)/8.0d0/TPI)
         endif
      CASE(HAjj_WW, HAjj_ZZ_ll, HAjj_ZZ_lnu)
         call LIPSNVV(n_p+1-nx,rand,rn,
     1                p(0,1,1),p(0,2,1),
     2                v(0,1,1),p(0,3,1),
     3                X(1),X(2),dps)
         v(0,5,1)=p(0,n_p+1,1)
         v(1,5,1)=p(1,n_p+1,1)
         v(2,5,1)=p(2,n_p+1,1)
         v(3,5,1)=p(3,n_p+1,1)
         p(0,n_p+1,1)=0d0
         p(1,n_p+1,1)=0d0
         p(2,n_p+1,1)=0d0
         p(3,n_p+1,1)=0d0
      CASE(Zjj_l, Zjj_nu, WPjj, WMjj) ! pp -> Vjj(j) -> lljj(j)
         call LIPSN(n_p-nx, rand, rn,
     1              p(0,1,1), p(0,2,1), 
     2              Q, v(0,1,1), v(0,2,1),
     3              p(0,3,1), X(1), X(2), dps, 0)
      CASE(WPWMjj,WPWPjj,WMWMjj,ZZjj_ll,ZZjj_lnu,WPZjj,WMZjj,
     1     WPhadWMjj, WPWMhadjj, WPhadZjj, WPZhadjj, WMhadZjj, WMZhadjj, ZZhadjj,
     2     WPhadWPjj, WMhadWMjj) ! pp -> VVjj(j)
         call LIPSNVV(n_p-nx,rand,rn,
     1                p(0,1,1),p(0,2,1),
     2                v(0,1,1),p(0,3,1),
     3                X(1),X(2),dps)
      CASE(HHjj) ! pp -> HH
         call LIPSNHH(n_p-nx,rand,rn,
     1                p(0,1,1),p(0,2,1),
     2                v(0,1,1),p(0,3,1),
     3                X(1),X(2),dps,1)
      CASE(HHjj_bbtau, HHjj_bbAA) ! pp -> HH -> b bbar XX
         call LIPSNHH(n_p-nx,rand,rn,
     1                p(0,1,1),p(0,2,1),
     2                v(0,1,1),p(0,3,1),
     3                X(1),X(2),dps,0)
      CASE(WPAjj, WMAjj, ZAjj, QCDWPAjj,QCDWMAjj)        !pp -> VAjj(j)
         if (ps_number.eq.1.and.PS_loops.ne.1) then
            call LIPSVAJ1(n_p-nx,rand,rn,
     1           p(0,1,1), p(0,2,1),
     2           v(0,1,1),p(0,3,1),
     3           X(1),X(2),dps,narrowwidth)
         else 
            call LIPSVAJ2(n_p-nx,rand,rn,
     1           p(0,1,1), p(0,2,1),
     2           v(0,1,1),p(0,3,1),
     3           X(1),X(2),dps,narrowwidth)
         endif
      CASE(ZAjj_n)        !pp -> VAjj(j)
            call LIPSVAJ2(n_p-nx,rand,rn,
     1           p(0,1,1), p(0,2,1),
     2           v(0,1,1),p(0,3,1),
     3           X(1),X(2),dps,narrowwidth)
      CASE(HJJJ,Hjjj_AA, Hjjj_mu, Hjjj_tau, Hjjj_bbar) ! pp -> H + 3 jets
        if (higgsscheme .eq. 0 .or. higgsscheme .eq. -1) then
         call LIPSN(n_p-nx, rand, rn,
     1              p(0,1,1), p(0,2,1), 
     2              Q, v(0,1,1), v(0,2,1),
     3              p(0,3,1), X(1), X(2), dps, 1)
        else
         call LIPSN(n_p-nx, rand, rn,
     1              p(0,1,1), p(0,2,1), 
     2              Q, v(0,1,1), v(0,2,1),
     3              p(0,3,1), X(1), X(2), dps, 0)
        endif
         if (ldoblha .and. (procID.eq.Hjjj)) then
           do i=0,3
             v(i,1,1) = v(i,1,1)+v(i,2,1)
             v(i,2,1) = 0d0
           enddo
           dps = dps / (xmg(6)/8.0d0/TPI)
         endif
      CASE(Hjjj_WW, Hjjj_ZZ_ll, Hjjj_ZZ_lnu) ! pp -> Hjjj(j) -> VVjjj(j)
         call LIPSNVV(n_p-nx,rand,rn,
     1                p(0,1,1),p(0,2,1),
     2                v(0,1,1),p(0,3,1),
     3                X(1),X(2),dps)  
      CASE(GFHjj) ! pp -> Hjj in GF   
         call LIPSN0(n_p-nx, rand, rn,
     1               p(0,1,1), p(0,2,1), 
     2               Q, p(0,3,1), X(1), X(2), dps)
         v(0,1,1) = Q(0)
         v(1,1,1) = Q(1)
         v(2,1,1) = Q(2)
         v(3,1,1) = Q(3)   
      CASE(GFHjj_AA, GFHjj_mu, GFHjj_tau, GFHjj_bbar)     
         call LIPSN(n_p-nx, rand, rn,
     1              p(0,1,1), p(0,2,1), 
     2              Q, v(0,1,1), v(0,2,1),
     3              p(0,3,1), X(1), X(2), dps, 1)  
      CASE(GFHjj_WW, GFHjj_ZZ_ll, GFHjj_ZZ_lnu)           
         call LIPSNVV(n_p-nx,rand,rn,
     1                p(0,1,1),p(0,2,1),
     2                v(0,1,1),p(0,3,1),
     3                X(1),X(2),dps)
      CASE(QCDZjj_l,QCDZjj_nu,QCDWPjj,QCDWMjj)     
         call LIPSN(n_p-nx, rand, rn,
     1              p(0,1,1), p(0,2,1), 
     2              Q, v(0,1,1), v(0,2,1),
     3              p(0,3,1), X(1), X(2), dps, 0)  

#ifdef WITH_QCDVV
      CASE(QCDZAjj_l)        !pp -> ZAjj(j)
         if (ps_number.eq.1.and.PS_loops.ne.1) then
            call LIPSVAJ1(n_p-nx,rand,rn,
     1           p(0,1,1), p(0,2,1),
     2           v(0,1,1),p(0,3,1),
     3           X(1),X(2),dps,narrowwidth)
         else 
            call LIPSVAJ2(n_p-nx,rand,rn,
     1           p(0,1,1), p(0,2,1),
     2           v(0,1,1),p(0,3,1),
     3           X(1),X(2),dps,narrowwidth)
         endif
      CASE(QCDZAjj_n)        !pp -> ZAjj(j)
            call LIPSVAJ2(n_p-nx,rand,rn,
     1           p(0,1,1), p(0,2,1),
     2           v(0,1,1),p(0,3,1),
     3           X(1),X(2),dps,narrowwidth)
      CASE(QCDAAjj)
         call LIPSAAJ(n_p-nx,rand,rn,
     1                p(0,1,1), p(0,2,1),
     2                v(0,1,1),p(0,3,1),
     3                X(1),X(2),dps)
      CASE(QCDWPWPjj,QCDWPZjj,QCDWMZjj,
     1     QCDWMWMjj,QCDZZjj_ll,QCDZZjj_lnu)
!       if(isub.ne.sub_number) then
!         if(sub_number.eq.2.and.(WWJJ_interf.ge.2).and.mdijmin.lt.1d0) then
!          call def_ps_WWW(3,4,3)
! !          call def_ps_ZZZj(3,4,3)
!         else
!          call def_ps_qqVVj(3,3)
!         endif
!         isub=sub_number
!       endif
!       if(sub_number.eq.2.and.(WWJJ_interf.ge.2).and.mdijmin.lt.1d0) then ! include Higgs resonance
!       rn=RandomNumber()
!          call LIPSWWW(0,rand,rn,
!      1                p(0,1,1), p(0,2,1),
!      2                v(0,1,1),p(0,3,1),
!      3                X(1),X(2),dps,0)
!         p(:,3:4,1)=v(:,3:4,1)
!         v(:,3:4,1)=v(:,5:6,1)
! ! !       rn=RandomNumber()
! ! !          call LIPSZZZJ(0,rand,rn,
! ! !      1                p(0,1,1), p(0,2,1),
! ! !      2                v(0,1,1),p(0,3,1),
! ! !      3                X(1),X(2),dps,0, ps_number)
! ! !         p(:,3:4,1)=v(:,3:4,1)
! ! !         v(:,3:4,1)=v(:,5:6,1)
! !   
!       else
        call LIPSVVJ(n_p-nx,rand,RN,
     1                   p(0,1,1),p(0,2,1),V(0,1,1),P(0,3,1),
     2                   X(1),X(2),dps,0)
!       endif
#endif

#ifdef WITH_DIBOSON
      CASE(WW,GFWW,WPhadWMlep,WPlepWMhad,GFWPhadWMlep,GFWPlepWMhad)
         call LIPSVVJ(n_p-nx,rand,rn,
     1                p(0,1,1), p(0,2,1),
     2                v(0,1,1),p(0,3,1),
     3                X(1),X(2),dps,narrowwidth)

      CASE(WPZ,WPhadZ,WPZhad)
         call LIPSVVJ(n_p-nx,rand,rn,
     1        p(0,1,1), p(0,2,1),
     2        v(0,1,1),p(0,3,1),
     3        X(1),X(2),dps,narrowwidth)

      CASE(WMZ,WMhadZ,WMZhad)
         call LIPSVVJ(n_p-nx,rand,rn,
     1        p(0,1,1), p(0,2,1),
     2        v(0,1,1),p(0,3,1),
     3        X(1),X(2),dps,narrowwidth)
         
      CASE(ZZ,ZZhad,GFZZ,GFZZhad)
         call LIPSVVJ(n_p-nx,rand,rn,
     1                p(0,1,1), p(0,2,1),
     2                v(0,1,1),p(0,3,1),
     3                X(1),X(2),dps,narrowwidth)

      CASE(WPA)
         if (ps_number.eq.1) then
           call LIPSVAJ1(n_p-nx,rand,rn,
     1                   p(0,1,1), p(0,2,1),
     2                   v(0,1,1),p(0,3,1),
     3                   X(1),X(2),dps,narrowwidth)
         else if (ps_number.eq.2) then
           call LIPSVAJ2(n_p-nx,rand,rn,
     1                   p(0,1,1), p(0,2,1),
     2                   v(0,1,1),p(0,3,1),
     3                   X(1),X(2),dps,narrowwidth)
         endif

      CASE(WMA)
         if (ps_number.eq.1) then
           call LIPSVAJ1(n_p-nx,rand,rn,
     1                   p(0,1,1), p(0,2,1),
     2                   v(0,1,1),p(0,3,1),
     3                   X(1),X(2),dps,narrowwidth)
         else if (ps_number.eq.2) then
           call LIPSVAJ2(n_p-nx,rand,rn,
     1                   p(0,1,1), p(0,2,1),
     2                   v(0,1,1),p(0,3,1),
     3                   X(1),X(2),dps,narrowwidth)
         endif


c       n_p = number of partons  = 2
c       n_v = number of decay particles = 2 for W+-, 3 for W+-H
c       nx =  number of internal x paramaters = 2
      CASE(
     7      WPH, WPH_AA, WPH_mu, WPH_tau, WPH_bbar, WPH_WW, WPH_ZZ_ll,  WPH_ZZ_lnu,
     7      WMH, WMH_AA, WMH_mu, WMH_tau, WMH_bbar, WMH_WW, WMH_ZZ_ll,  WMH_ZZ_lnu
     &    )
         call LIPSVHJ(n_p-nx,rand,rn,
     1                   p(0,1,1), p(0,2,1),
     2                   v(0,1,1), p(0,3,1),
     3                   X(1),X(2),dps,narrowwidth,0) !last param: narrowwidthhiggs

      ! could also use LIPSVNJ here  
      CASE(WP_only,WM_only)
         call LIPSVJ(n_p-nx,rand,rn,
     1                   p(0,1,1), p(0,2,1),
     2                   v(0,1,1), p(0,3,1),
     3                   X(1),X(2),dps,narrowwidth)

      CASE(ZA,GFZA)
c flipped because for GFZA only PS1 relevant
         if (ps_number.eq.1) then
           call LIPSVAJ1(n_p-nx,rand,rn,
     1                   p(0,1,1), p(0,2,1),
     2                   v(0,1,1),p(0,3,1),
     3                   X(1),X(2),dps,narrowwidth)
         else if (ps_number.eq.2) then
           call LIPSVAJ2(n_p-nx,rand,rn,
     1                   p(0,1,1), p(0,2,1),
     2                   v(0,1,1),p(0,3,1),
     3                   X(1),X(2),dps,narrowwidth)
         endif

      CASE(AA,GFAA)
         call LIPSAAJ(n_p-nx,rand,rn,
     1                p(0,1,1), p(0,2,1),
     2                v(0,1,1),p(0,3,1),
     3                X(1),X(2),dps)
#endif

#ifdef WITH_DIBOSONJET
      CASE(WWJ,WPHADWMJ,WPWMHADJ,GFWWj,GFWPHADWMJ,GFWPWMHADJ)
         call LIPSVVJ(n_p-nx,rand,rn,
     1                p(0,1,1), p(0,2,1),
     2                v(0,1,1),p(0,3,1),
     3                X(1),X(2),dps,narrowwidth)
      CASE(WMZJ, WMhadZJ, WMZhadJ, WPZJ, WPhadZJ, WPZhadJ)

         call LIPSNVV(n_p-nx,rand,rn,
     1                p(0,1,1),p(0,2,1),
     2                v(0,1,1),p(0,3,1),
     3                X(1),X(2),dps)
!          call LIPSVVJ(n_p-nx,rand,rn,
!      1                p(0,1,1), p(0,2,1),
!      2                v(0,1,1),p(0,3,1),
!      3                X(1),X(2),dps,narrowwidth)
        
      CASE(WMAJ,WPAJ)
         if (ps_number.eq.1) then
            call LIPSNthreebod(n_p-nx, rand, rn,
     1           p(0,1,1), p(0,2,1), 
     2           Q, v(0,1,1), v(0,2,1),v(0,3,1),
     3           p(0,3,1), X(1), X(2), dps)
         else if (ps_number.eq.2) then
            call LIPSN(n_p-nx+1, rand, rn,
     1           p(0,1,1), p(0,2,1), 
     2           Q, v(0,1,1), v(0,2,1),
     3           p(0,3,1), X(1), X(2), dps, 0)
            do mu=0,3
               v(mu,3,1) = p(mu,n_p+1,1) !Generate photon momenta as jet momenta
               p(mu,n_p+1,1) = 0d0
            enddo        
         endif

      ! could also use LIPSVNJ here  
      CASE(WPJ,WMJ)
         call LIPSVJ(n_p-nx,rand,rn,
     1                   p(0,1,1), p(0,2,1),
     2                   v(0,1,1), p(0,3,1),
     3                   X(1),X(2),dps,narrowwidth)

      CASE(
     7      WPHJ, WPHJ_AA, WPHJ_mu, WPHJ_tau, WPHJ_bbar, WPHJ_WW, WPHJ_ZZ_ll,  WPHJ_ZZ_lnu,
     7      WMHJ, WMHJ_AA, WMHJ_mu, WMHJ_tau, WMHJ_bbar, WMHJ_WW, WMHJ_ZZ_ll,  WMHJ_ZZ_lnu
     &    )    
         call LIPSVHJ(n_p-nx,rand,rn,
     1                   p(0,1,1), p(0,2,1),
     2                   v(0,1,1), p(0,3,1),
     3                   X(1),X(2),dps,narrowwidth,0) !last param: narrowwidthhiggs

      CASE(ZZJ,GFZZj)
         call LIPSVVJ(n_p-nx,rand,rn,
     1                p(0,1,1), p(0,2,1),
     2                v(0,1,1),p(0,3,1),
     3                X(1),X(2),dps,narrowwidth)

#endif

#ifdef WITH_TRIBOSON

        CASE(WWZ, WPhadWMZ, WPWMhadZ, WWZhad)
#ifdef WITH_KK
        if (with_kk) then
         call LIPSWWZ_kk(n_p-nx,rand,rn,
     1                p(0,1,1), p(0,2,1),
     2                v(0,1,1),p(0,3,1),
     3                X(1),X(2),dps,narrowwidth)
        else
#endif
         call LIPSVVVJ(n_p-nx,rand,rn,
     1                p(0,1,1), p(0,2,1),
     2                v(0,1,1),p(0,3,1),
     3                X(1),X(2),dps,narrowwidth)
#ifdef WITH_KK
        endif
#endif
      CASE(ZZWP, ZZhadWP, ZZWPhad, ZZWM, ZZhadWM, ZZWMhad)
#ifdef WITH_KK
        if (with_kk) then
        
         call LIPSZZW_kk(n_p-nx,rand,rn,
     1                p(0,1,1), p(0,2,1),
     2                v(0,1,1),p(0,3,1),
     3                X(1),X(2),dps,narrowwidth)
        else
#endif
         call LIPSVVVJ(n_p-nx,rand,rn,
     1                p(0,1,1), p(0,2,1),
     2                v(0,1,1),p(0,3,1),
     3                X(1),X(2),dps,narrowwidth)
#ifdef WITH_KK
        endif
#endif

      CASE(WWWP, WPhadWMWP, WPWMhadWP, WWWM, WMhadWPWM, WMWPhadWM)
#ifdef WITH_KK
        if (with_kk) then
         call LIPSWWW_kk(n_p-nx,rand,rn,
     1                p(0,1,1), p(0,2,1),
     2                v(0,1,1),p(0,3,1),
     3                X(1),X(2),dps,narrowwidth)
        else
#endif
         call LIPSWWW(n_p-nx,rand,rn,
     1                p(0,1,1), p(0,2,1),
     2                v(0,1,1),p(0,3,1),
     3                X(1),X(2),dps,narrowwidth)
#ifdef WITH_KK
        endif
#endif

      CASE(ZZZ, ZZZhad)
         call LIPSZZZJ(n_p-nx,rand,rn,
     1                p(0,1,1), p(0,2,1),
     2                v(0,1,1),p(0,3,1),
     3                X(1),X(2),dps,narrowwidth, ps_number)


      CASE(ZZA, ZZhadA, ZZnA)
            if (ps_number.eq.1) then
               call LIPSVVAJ(n_p-nx,rand,rn,
     1              p(0,1,1), p(0,2,1),
     2              v(0,1,1),p(0,3,1),
     3              X(1),X(2),dps)
            elseif (ps_number.eq.2) then
               call LIPSVVAJ1(n_p-nx,rand,rn,
     1              p(0,1,1), p(0,2,1),
     2              v(0,1,1),p(0,3,1),
     3              X(1),X(2),dps)
            elseif (ps_number.eq.3) then
               call LIPSVVAJ2(n_p-nx,rand,rn,
     1              p(0,1,1), p(0,2,1),
     2              v(0,1,1),p(0,3,1),
     3              X(1),X(2),dps)
            endif


      CASE(WWA, WPhadWMA, WPWMhadA)
            if (ps_number.eq.1) then
               call LIPSVVAJ(n_p-nx,rand,rn,
     1              p(0,1,1), p(0,2,1),
     2              v(0,1,1),p(0,3,1),
     3              X(1),X(2),dps)
            elseif (ps_number.eq.2) then
               call LIPSVVAJ1(n_p-nx,rand,rn,
     1              p(0,1,1), p(0,2,1),
     2              v(0,1,1),p(0,3,1),
     3              X(1),X(2),dps)
            elseif (ps_number.eq.3) then
               call LIPSVVAJ2(n_p-nx,rand,rn,
     1              p(0,1,1), p(0,2,1),
     2              v(0,1,1),p(0,3,1),
     3              X(1),X(2),dps)
            endif


      CASE(WPZA, WMZA, WPhadZA, WPZhadA, WMhadZA, WMZhadA)
            if (ps_number.eq.1) then
               call LIPSVVAJ(n_p-nx,rand,rn,
     1              p(0,1,1), p(0,2,1),
     2              v(0,1,1),p(0,3,1),
     3              X(1),X(2),dps)
            elseif (ps_number.eq.2) then
               call LIPSVVAJ1(n_p-nx,rand,rn,
     1              p(0,1,1), p(0,2,1),
     2              v(0,1,1),p(0,3,1),
     3              X(1),X(2),dps)
            elseif (ps_number.eq.3) then
               call LIPSVVAJ2(n_p-nx,rand,rn,
     1              p(0,1,1), p(0,2,1),
     2              v(0,1,1),p(0,3,1),
     3              X(1),X(2),dps)
            endif

      CASE(WPAA, WMAA)
            if (ps_number.eq.1) then
               call LIPSVAAJ1(n_p-nx,rand,rn,
     1              p(0,1,1), p(0,2,1),
     2              v(0,1,1),p(0,3,1),
     3              X(1),X(2),dps,narrowwidth)
            elseif (ps_number.eq.2) then
               call LIPSVAAJ2(n_p-nx,rand,rn,
     1              p(0,1,1), p(0,2,1),
     2              v(0,1,1),p(0,3,1),
     3              X(1),X(2),dps,narrowwidth)
            elseif (ps_number.eq.3) then
               call LIPSVAAJ3(n_p-nx,rand,rn,
     1              p(0,1,1), p(0,2,1),
     2              v(0,1,1),p(0,3,1),
     3              X(1),X(2),dps,narrowwidth)
            elseif (ps_number.eq.4) then
               call LIPSVAAJ4(n_p-nx,rand,rn,
     1              p(0,1,1), p(0,2,1),
     2              v(0,1,1),p(0,3,1),
     3              X(1),X(2),dps,narrowwidth)
            endif
            if(with_anom .and. ps_number.eq.5) then   !one additional ps for H -> A A contribution
               call LIPSVAAJ5_anomal(n_p-nx,rand,rn,
     1              p(0,1,1), p(0,2,1),
     2              v(0,1,1),p(0,3,1),
     3              X(1),X(2),dps,narrowwidth)
            endif

      CASE(ZAA)
            if (ps_number.eq.1) then
               call LIPSVAAJ1(n_p-nx,rand,rn,
     1              p(0,1,1), p(0,2,1),
     2              v(0,1,1),p(0,3,1),
     3              X(1),X(2),dps,narrowwidth)
            elseif (ps_number.eq.2) then
               call LIPSVAAJ2(n_p-nx,rand,rn,
     1              p(0,1,1), p(0,2,1),
     2              v(0,1,1),p(0,3,1),
     3              X(1),X(2),dps,narrowwidth)
            elseif (ps_number.eq.3) then
               call LIPSVAAJ3(n_p-nx,rand,rn,
     1              p(0,1,1), p(0,2,1),
     2              v(0,1,1),p(0,3,1),
     3              X(1),X(2),dps,narrowwidth)
            elseif (ps_number.eq.4) then
               call LIPSVAAJ4(n_p-nx,rand,rn,
     1              p(0,1,1), p(0,2,1),
     2              v(0,1,1),p(0,3,1),
     3              X(1),X(2),dps,narrowwidth)
            endif
            if(with_anom .and. ps_number.eq.5) then   !one additional ps for H -> A A contribution
               call LIPSVAAJ5_anomal(n_p-nx,rand,rn,
     1              p(0,1,1), p(0,2,1),
     2              v(0,1,1),p(0,3,1),
     3              X(1),X(2),dps,narrowwidth)
            endif

      CASE(ZnAA)
         if(with_anom) then !multichannel due to new Higgs couplings
            if (ps_number.eq.1) then
               call LIPSVAAJ1(n_p-nx,rand,rn,
     1              p(0,1,1), p(0,2,1),
     2              v(0,1,1),p(0,3,1),
     3              X(1),X(2),dps,narrowwidth)
            elseif (ps_number.eq.2) then !for H -> Z A1
               call LIPSVAAJ2(n_p-nx,rand,rn,
     1              p(0,1,1), p(0,2,1),
     2              v(0,1,1),p(0,3,1),
     3              X(1),X(2),dps,narrowwidth)
            elseif (ps_number.eq.3) then !for H -> Z A2
               call LIPSVAAJ3(n_p-nx,rand,rn,
     1              p(0,1,1), p(0,2,1),
     2              v(0,1,1),p(0,3,1),
     3              X(1),X(2),dps,narrowwidth)
            elseif (ps_number.eq.4) then !for H -> A A
               call LIPSVAAJ5_anomal(n_p-nx,rand,rn,
     1              p(0,1,1), p(0,2,1),
     2              v(0,1,1),p(0,3,1),
     3              X(1),X(2),dps,narrowwidth)
            endif
         else
            call LIPSVAAJ1(n_p-nx,rand,rn,
     1           p(0,1,1), p(0,2,1),
     2           v(0,1,1),p(0,3,1),
     3           X(1),X(2),dps,narrowwidth)
         endif


      CASE(AAA)
         if(with_anom) then !multichannel due to new Higgs couplings
            if (ps_number.eq.1) then
               call LIPSAAAJ(n_p-nx,rand,rn,
     1              p(0,1,1), p(0,2,1),
     2              v(0,1,1),p(0,3,1),
     3              X(1),X(2),dps)
            elseif (ps_number.eq.2) then !for H -> A1 A2
               call LIPSAAAJ2_anomal(n_p-nx,rand,rn,
     1              p(0,1,1), p(0,2,1),
     2              v(0,1,1),p(0,3,1),
     3              X(1),X(2),dps)
            elseif (ps_number.eq.3) then !for H -> A1 A3
               call LIPSAAAJ3_anomal(n_p-nx,rand,rn,
     1              p(0,1,1), p(0,2,1),
     2              v(0,1,1),p(0,3,1),
     3              X(1),X(2),dps)
            elseif (ps_number.eq.4) then !for H -> A2 A3
               call LIPSAAAJ4_anomal(n_p-nx,rand,rn,
     1              p(0,1,1), p(0,2,1),
     2              v(0,1,1),p(0,3,1),
     3              X(1),X(2),dps)
            endif
         else
            call LIPSAAAJ(n_p-nx,rand,rn,
     1           p(0,1,1), p(0,2,1),
     2           v(0,1,1),p(0,3,1),
     3           X(1),X(2),dps)
         endif
#endif

#ifdef WITH_TRIBOSONJET
      CASE(WPAAj, WMAAj)
            if (ps_number.eq.1) then
               call LIPSVAAJ1(n_p-nx,rand,rn,
     1              p(0,1,1), p(0,2,1),
     2              v(0,1,1),p(0,3,1),
     3              X(1),X(2),dps,narrowwidth)
            elseif (ps_number.eq.2) then
               call LIPSVAAJ2(n_p-nx,rand,rn,
     1              p(0,1,1), p(0,2,1),
     2              v(0,1,1),p(0,3,1),
     3              X(1),X(2),dps,narrowwidth)
            elseif (ps_number.eq.3) then
               call LIPSVAAJ3(n_p-nx,rand,rn,
     1              p(0,1,1), p(0,2,1),
     2              v(0,1,1),p(0,3,1),
     3              X(1),X(2),dps,narrowwidth)
            elseif (ps_number.eq.4) then
               call LIPSVAAJ4(n_p-nx,rand,rn,
     1              p(0,1,1), p(0,2,1),
     2              v(0,1,1),p(0,3,1),
     3              X(1),X(2),dps,narrowwidth)
            endif
            if(with_anom .and. ps_number.eq.5) then   !one additional ps for H -> A A contribution
               call LIPSVAAJ5_anomal(n_p-nx,rand,rn,
     1              p(0,1,1), p(0,2,1),
     2              v(0,1,1),p(0,3,1),
     3              X(1),X(2),dps,narrowwidth)
            endif
#endif

      CASE DEFAULT
         print*," No phase space generator is called"
         stop   
      END SELECT       


c get additional random number not driven by VEGAS    
      if (.not.ldoblha) then
        rn = RandomNumber()
      else
        rn = rand(PS_dimension+1)
      endif
c globally rotate by angle(rn)
      phirn = TPI*rn
      cosrn = cos(phirn)
      sinrn = sin(phirn)
      do i=nx+1,n_p
        tmpx =  cosrn*p(1,i,1)+sinrn*p(2,i,1)
        tmpy = -sinrn*p(1,i,1)+cosrn*p(2,i,1)
        p(1,i,1) = tmpx
        p(2,i,1) = tmpy
      enddo
      do i=1,n_v
        tmpx =  cosrn*v(1,i,1)+sinrn*v(2,i,1)
        tmpy = -sinrn*v(1,i,1)+cosrn*v(2,i,1)
        v(1,i,1) = tmpx
        v(2,i,1) = tmpy
      enddo


      ! avoid singularity for small photon virtuality
      ! for semileptonic processes, when Z/gamma decays hadronically
      SELECT CASE(procid)
      CASE(ZZhadjj, ZZhad, GFZZhad, ZZhadWP, ZZhadWM, ZZZhad, Hjj_ZZhad, WPhadWMA, ZZhadA, WPhadZA, WMhadZA)
         if (mjj2(v(0,1,1),v(0,2,1)).lt.qsqAmin_had) dps=0.d0
      CASE(WPZhad, WMZhad, WPZhadjj, WMZhadjj, WPWMhadA, WPZhadA, WMZhadA)
         if (mjj2(v(0,3,1),v(0,4,1)).lt.qsqAmin_had) dps=0.d0
      CASE(WWZhad)
         if (mjj2(v(0,5,1),v(0,6,1)).lt.qsqAmin_had) dps=0.d0
      END SELECT
c      if (dps.gt.0.d0) print*, sqrt(mjj2(v(0,1,1),v(0,2,1)))
c      if ( sqrt((v(1,1,1)+v(1,2,1))**2 + (v(2,1,1)+v(2,2,1))**2) .lt. 20.0) dps=0.d0 


      if (.not.ldoblha) then
c Multiply with symmetry factor for identical leptons in the final state
        dps = dps*symmetryfactor
      endif

      RETURN
      END



