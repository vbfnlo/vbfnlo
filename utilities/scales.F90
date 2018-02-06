!***************************************************************************
!     This file contains routines for the calculation of the factorization
!     and renormalization scale, alpha_s and momentum transfers.
!***************************************************************************
!   LIST OF ALL FUNCTIONS AND SUBROUTINES IN THIS FILE:
!
!     SUBROUTINE InitScales
!     SUBROUTINE Scales
!     SUBROUTINE CalcVBFScales(p, v, Lmax, nlo, bos)
!     SUBROUTINE calcH3jscales(p, v, jets, njets, leptons, nl, photons, 
!                                             nphotons, Lmax, nlo)
!     SUBROUTINE calcGGFscales(p)
!     SUBROUTINE calcVXjscales(p, v, jets, njets, leptons, nl, photons, nphotons, nlo, Lmax)
!     SUBROUTINE calcQCDscales(p, v, jets, njets, leptons, nl, photons, 
!                                             nphotons, Lmax, nlo)
!     SUBROUTINE Calc_Momentum_Transfer(p, v, q12, q34, Lmax)
!     double precision function transverseEnergy(p)
!
!     FUNCTION alphas5(q2ren,iord1)
!     FUNCTION alphas5_hardwired(q2ren,iord1)
!     FUNCTION alphas(q2ren,xnf,iord1)
!     FUNCTION alphas_mstw(q2ren,iord1)  
!            ... followed by other downloaded mstw2008 alphas routines
!
!***************************************************************************
!   Last modified: 03.05.2011
!***************************************************************************

!***************************************************************************
      SUBROUTINE InitScales
!***************************************************************************
!     Initialization of renormalization and factorization scales
!     and printout for scales.
!***************************************************************************
          use globalvars, only: lglobalprint, ldoscales

      implicit none
#include "scales.inc"
#include "global.inc"
#include "koppln.inc"
#include "process.inc"

      real*8 qsq

      real*8 alphas5, alphas5_hardwired
      external alphas5, alphas5_hardwired

!* clrCT is the counterterm for qqV
      double complex clrCT
      common /bkopouCT/ clrCT(3:4,2:4,-1:1)
      DOUBLE PRECISION CLR,XM2,XMG,Bk,Vk,Ak
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),Bk(6,6,6), &
           Vk(4,5),Ak(4,5)

      if (ldoscales) then

         if (lglobalprint) then
            print*," "
            print*,"Information on scale setting in InitScales "
            print*,"-------------------------------------------"
            print*," "
         endif

!-------------------------------      
!   factorization scale
!-------------------------------
      SELECT CASE(id_muf)
! user defined constant scale 
      CASE(0)
         qsq = muf_user**2 * xif**2
         if (lglobalprint) print*,"Constant factorization scale chosen: mu_F =", &
                         sqrt(qsq)," GeV "
! dynamic scale = momentum transfer of W/Z
      CASE(1)
         SELECT CASE(procID)
         CASE(Hjj, Hjj_AA, Hjj_mu, Hjj_tau, Hjj_bbar, Hjj_WW, Hjj_ZZ_ll, &
              Hjj_ZZ_lnu, HHjj, HHjj_bbtau, HHjj_bbAA, Ajj, Zjj_l, Zjj_nu, WPjj, WMjj, WPWMjj, ZZjj_ll, &
              ZZjj_lnu, WPZjj, WMZjj, HAjj, HAjj_AA, HAjj_mu, HAjj_tau,  &
              HAjj_bbar, HAjj_WW, HAjj_ZZ_ll, HAjj_ZZ_lnu, AAjj, WPWPjj,  &
              WMWMjj, Sp2jj_WW, Sp2jj_ZZ_ll, Sp2jj_ZZ_lnu,  &
              WPhadWMjj, WPWMhadjj, WPhadZjj, WPZhadjj, WMhadZjj, WMZhadjj, ZZhadjj, &
              Hjj_WPhadWM, Hjj_WPWMhad, Hjj_ZZhad, WPhadWPjj, WMhadWMjj, &
              WPAjj, WMAjj, ZAjj, ZAjj_n)
            if (lglobalprint) print*,"Dynamical factorization scale chosen:"
            if (lglobalprint) print*,"        momentum transfer of W/Z  "
         CASE DEFAULT
            print *,"unreasonable scale choice for this process:  &
                  id_muf = ",id_muf
            STOP
         END SELECT

! dynamic scale = min(pt_j)
      CASE(2)
         SELECT CASE(procID)
         CASE(Hjj, Hjj_AA, Hjj_mu, Hjj_tau, Hjj_bbar, Hjj_WW, Hjj_ZZ_ll, &
              Hjj_ZZ_lnu, HHjj, HHjj_bbtau, HHjj_bbAA, Ajj, Zjj_l, Zjj_nu, WPjj, WMjj, WPWMjj, ZZjj_ll, &
              ZZjj_lnu, Hjjj, Hjjj_AA, Hjjj_mu, Hjjj_tau, Hjjj_bbar, Hjjj_WW,  &
              Hjjj_ZZ_ll,Hjjj_ZZ_lnu, WPZjj, WMZjj,WPAAj, WMAAj, HAjj, HAjj_AA, HAjj_mu, &
              HAjj_tau, HAjj_bbar, HAjj_WW, HAjj_ZZ_ll, HAjj_ZZ_lnu, &
              WPWPjj,WMWMjj,AAjj, &
              WPhadWMjj, WPWMhadjj, WPhadZjj, WPZhadjj, WMhadZjj, WMZhadjj, ZZhadjj, &
              Hjj_WPhadWM, Hjj_WPWMhad, Hjj_ZZhad, WPhadWPjj, WMhadWMjj,  &
              Sp2jj_WW, Sp2jj_ZZ_ll, Sp2jj_ZZ_lnu, WPAjj, WMAjj, ZAjj, ZAjj_n)


            if (lglobalprint) print*,"Dynamical factorization scale chosen:"
            if (lglobalprint) print*,"              mu_F = min(pt_j)"
         CASE DEFAULT
            print *,"unreasonable scale choice for this process:  &
                  id_muf = ",id_muf
            STOP
         END SELECT

! dynamic scale = invariant primary boson mass
      CASE(3,4)
         SELECT CASE(procID)
         CASE (WP_only,WM_only,WPJ,WMJ)
            if (lglobalprint) print*,"Dynamical factorization scale chosen:"
            if (lglobalprint) print*,"      invariant V mass "
         CASE ( &
           WPH, WPH_AA, WPH_mu, WPH_tau, WPH_bbar, WPH_WW, WPH_ZZ_ll,  WPH_ZZ_lnu, &
           WMH, WMH_AA, WMH_mu, WMH_tau, WMH_bbar, WMH_WW, WMH_ZZ_ll,  WMH_ZZ_lnu, &
           WPHJ, WPHJ_AA, WPHJ_mu, WPHJ_tau, WPHJ_bbar, WPHJ_WW, WPHJ_ZZ_ll,  WPHJ_ZZ_lnu, &
           WMHJ, WMHJ_AA, WMHJ_mu, WMHJ_tau, WMHJ_bbar, WMHJ_WW, WMHJ_ZZ_ll,  WMHJ_ZZ_lnu)
            if (lglobalprint) print*,"Dynamical factorization scale chosen:"
            if (lglobalprint) print*,"      invariant VH mass "
         CASE (WW, GFWW, WPZ, WMZ, WPA, WMA, ZZ, GFZZ, ZA, GFZA, AA, GFAA, &
              WPhadWMlep,WPlepWMhad,GFWPhadWMlep,GFWPlepWMhad, &
              ZZhad,GFZZhad, &
              WPhadZ,WPZhad,WMhadZ,WMZhad, &
              WPAJ, WMAJ, WPZJ, WMZJ, &
              WWJ, WPHADWMJ, WPWMHADJ, GFWWj, GFWPHADWMJ, GFWPWMHADJ, &
              ZZJ, GFZZj &
           )       
            if (lglobalprint) print*,"Dynamical factorization scale chosen:"
            if (lglobalprint) print*,"      invariant VV mass "
         CASE(ZZWP, ZZWM, WWZ, WWWP, WWWM, ZZZ, WWA, ZZA, ZZnA, WPZA, WMZA, &
              WPAA, WMAA, ZAA,ZnAA, AAA, WPAAj, WMAAj, &
              WPhadWMZ, WPWMhadZ, WWZhad, ZZhadWP, ZZWPhad, ZZhadWM,  &
              ZZWMhad, WPhadWMWP, WPWMhadWP, WMhadWPWM, WMWPhadWM, ZZZhad, &
              WPhadWMA, WPWMhadA, ZZhadA, WPhadZA, WPZhadA, WMhadZA, WMZhadA)
            if (lglobalprint) print*,"Dynamical factorization scale chosen:"
            if (lglobalprint) print*,"      invariant VVV mass "
         CASE DEFAULT
            print *,"unreasonable scale choice for this process: &
                  id_muf = ",id_muf
            STOP
         END SELECT

! dynamic scale = sqrt(pt_j1*pt_j2)
      CASE(5)
         SELECT CASE(procID)
         CASE(GFHjj, GFHjj_AA, GFHjj_mu, GFHjj_tau, GFHjj_bbar, &
              GFHjj_WW, GFHjj_ZZ_ll, GFHjj_ZZ_lnu)
            if (lglobalprint) print*,"Dynamical factorization scale chosen:"
            if (lglobalprint) print*,"              mu_F = sqrt(pt_j1*pt_j2)"
         CASE DEFAULT
            print *,"unreasonable scale choice for this process:  &
                  id_muf = ",id_muf
            STOP
         END SELECT

! constant scale = higgs mass
      CASE(6)
         qsq = xmh**2 * xif**2
         if (lglobalprint) print*,"Constant factorization scale chosen: mu_F = MH = ", &
                         sqrt(qsq)," GeV "

! dynamic scale = minimal ET
      CASE(7)
         SELECT CASE(procID)
         CASE(WW, GFWW, WWZ, ZZWP, ZZWM, WWWP, WWWM, ZZZ, WWA, ZZA, ZZnA, WPZA,  &
              WMZA, WPAA, WMAA, ZAA, ZnAA, AAA, WPAAj, WMAAj, &
              WPA,WMA,WPZ,WMZ,ZZ,GFZZ,ZA,GFZA,AA,GFAA, &
              WPhadWMlep,WPlepWMhad,GFWPhadWMlep,GFWPlepWMhad, &
              ZZhad,GFZZhad, &
              WPhadZ,WPZhad,WMhadZ,WMZhad, &
              WPhadWMZ, WPWMhadZ, WWZhad, ZZhadWP, ZZWPhad, ZZhadWM,  &
              ZZWMhad, WPhadWMWP, WPWMhadWP, WMhadWPWM, WMWPhadWM, ZZZhad, &
              WPhadWMA, WPWMhadA, ZZhadA, WPhadZA, WPZhadA, WMhadZA, WMZhadA, &
              WPAJ, WMAJ, WPZJ, WMZJ, &
              WWJ, WPHADWMJ, WPWMHADJ, GFWWj, GFWPHADWMJ, GFWPWMHADJ, &
              ZZJ, GFZZj, &
              WPJ,WMJ &
             )
            if (lglobalprint) print*,"Dynamical factorization scale chosen:"
            if (lglobalprint) print*,"      mu_F = min(ET_V) "
         CASE( &
           WPH, WPH_AA, WPH_mu, WPH_tau, WPH_bbar, WPH_WW, WPH_ZZ_ll,  WPH_ZZ_lnu, &
           WMH, WMH_AA, WMH_mu, WMH_tau, WMH_bbar, WMH_WW, WMH_ZZ_ll,  WMH_ZZ_lnu, &
           WPHJ, WPHJ_AA, WPHJ_mu, WPHJ_tau, WPHJ_bbar, WPHJ_WW, WPHJ_ZZ_ll,  WPHJ_ZZ_lnu, &
           WMHJ, WMHJ_AA, WMHJ_mu, WMHJ_tau, WMHJ_bbar, WMHJ_WW, WMHJ_ZZ_ll,  WMHJ_ZZ_lnu &
             )
            if (lglobalprint) print*,"Dynamical factorization scale chosen:"
            if (lglobalprint) print*,"      mu_F = min(ET_V, ET_H) "
         CASE DEFAULT
            print *,"unreasonable scale choice for this process:  &
                  id_muf = ",id_muf
            STOP
         END SELECT

! dynamical scale HT/2,  HT = sum(pT partons) + sum( ET V )
      CASE(8)
         SELECT CASE(procID)
         CASE(QCDWMWMjj,QCDWPWPjj,QCDWPjj,QCDWMjj,QCDWPZjj,QCDWMZjj,QCDAAjj, &
              QCDZjj_l, QCDZjj_nu,QCDWPAjj,QCDWMAjj,QCDZZjj_ll,QCDZZjj_lnu, &
              QCDZAjj_l,QCDZAjj_n, &
              WPZ,WMZ, WPZJ, WMZJ, &
              WPA,WPAJ,WMA,WMAJ, &
              WW,WWJ,GFWW,GFWWj, &
              ZZ,ZZJ,GFZZ,GFZZj, &
              ZA,GFZA, &
              AA,GFAA , &
              Hjj, Hjj_AA, Hjj_mu, Hjj_tau, Hjj_bbar, Hjj_WW, Hjj_ZZ_ll, &
              Hjj_ZZ_lnu, HHjj, HHjj_bbtau, HHjj_bbAA, Ajj, Zjj_l, Zjj_nu, WPjj, WMjj, WPWMjj, ZZjj_ll, &
              ZZjj_lnu, Hjjj, Hjjj_AA, Hjjj_mu, Hjjj_tau, Hjjj_bbar, Hjjj_WW,  &
              Hjjj_ZZ_ll,Hjjj_ZZ_lnu, WPZjj, WMZjj,WPAAj, WMAAj, HAjj, HAjj_AA, HAjj_mu, &
              HAjj_tau, HAjj_bbar, HAjj_WW, HAjj_ZZ_ll, HAjj_ZZ_lnu, &
              WPWPjj,WMWMjj,AAjj, &
              WPhadWMjj, WPWMhadjj, WPhadZjj, WPZhadjj, WMhadZjj, WMZhadjj, ZZhadjj, &
              Hjj_WPhadWM, Hjj_WPWMhad, Hjj_ZZhad, WPhadWPjj, WMhadWMjj,  &
              Sp2jj_WW, Sp2jj_ZZ_ll, Sp2jj_ZZ_lnu, WPAjj, WMAjj, ZAjj, ZAjj_n)
            if (lglobalprint) print*,"Dynamical factorization scale chosen:"
            if (lglobalprint) print*,"      mu_F = ( sum(pT partons) + sum( ET V ) )/2 "
       CASE DEFAULT
            print *,"unreasonable scale choice for this process:  &
                  id_muf = ",id_muf
            STOP
         END SELECT
         
! dynamical scale sum(pT_j*exp(|y_j-y12|)) + sum( ET V ) ( see arXiv 1311.6738 )
      CASE(9)
         SELECT CASE(procID)
         CASE(QCDWMWMjj,QCDWPWPjj,QCDWPjj,QCDWMjj,QCDWPZjj,QCDWMZjj,QCDAAjj, &
              QCDZjj_l, QCDZjj_nu,QCDWPAjj,QCDWMAjj,QCDZZjj_ll,QCDZZjj_lnu,QCDZAjj_l,QCDZAjj_n, &
              Hjj, Hjj_AA, Hjj_mu, Hjj_tau, Hjj_bbar, Hjj_WW, Hjj_ZZ_ll, &
              Hjj_ZZ_lnu, HHjj, HHjj_bbtau, HHjj_bbAA, Ajj, Zjj_l, Zjj_nu, WPjj, WMjj, WPWMjj, ZZjj_ll, &
              ZZjj_lnu, Hjjj, Hjjj_AA, Hjjj_mu, Hjjj_tau, Hjjj_bbar, Hjjj_WW,  &
              Hjjj_ZZ_ll,Hjjj_ZZ_lnu, WPZjj, WMZjj,WPAAj, WMAAj, HAjj, HAjj_AA, HAjj_mu, &
              HAjj_tau, HAjj_bbar, HAjj_WW, HAjj_ZZ_ll, HAjj_ZZ_lnu, &
              WPWPjj,WMWMjj,AAjj, &
              WPhadWMjj, WPWMhadjj, WPhadZjj, WPZhadjj, WMhadZjj, WMZhadjj, ZZhadjj, &
              Hjj_WPhadWM, Hjj_WPWMhad, Hjj_ZZhad, WPhadWPjj, WMhadWMjj,  &
              Sp2jj_WW, Sp2jj_ZZ_ll, Sp2jj_ZZ_lnu, WPAjj, WMAjj, ZAjj, ZAjj_n)
            if (lglobalprint) print*,"Dynamical factorization scale chosen:"
            if (lglobalprint) print*,"      mu_F = ( sum(pT_j*exp(|y_j-y12|)) + sum( ET V ) )/2"
            if (lglobalprint) print*,"( see arXiv 1311.6738 ) "
         CASE DEFAULT
            print *,"unreasonable scale choice for this process:  &
                  id_muf = ",id_muf
            STOP
         END SELECT
         
! dynamical scale HT/2,  HT = sum(pT partons) + sum( ET V )
      CASE(10)
         SELECT CASE(procID)
         CASE(QCDWMWMjj,QCDWPWPjj,QCDWPjj,QCDWMjj,QCDWPZjj,QCDWMZjj,QCDAAjj, &
              QCDZjj_l, QCDZjj_nu,QCDWPAjj,QCDWMAjj,QCDZZjj_ll,QCDZZjj_lnu,QCDZAjj_l,QCDZAjj_n, &
              Hjj, Hjj_AA, Hjj_mu, Hjj_tau, Hjj_bbar, Hjj_WW, Hjj_ZZ_ll, &
              Hjj_ZZ_lnu, HHjj, HHjj_bbtau, HHjj_bbAA, Ajj, Zjj_l, Zjj_nu, WPjj, WMjj, WPWMjj, ZZjj_ll, &
              ZZjj_lnu, Hjjj, Hjjj_AA, Hjjj_mu, Hjjj_tau, Hjjj_bbar, Hjjj_WW,  &
              Hjjj_ZZ_ll,Hjjj_ZZ_lnu, WPZjj, WMZjj,WPAAj, WMAAj, HAjj, HAjj_AA, HAjj_mu, &
              HAjj_tau, HAjj_bbar, HAjj_WW, HAjj_ZZ_ll, HAjj_ZZ_lnu, &
              WPWPjj,WMWMjj,AAjj, &
              WPhadWMjj, WPWMhadjj, WPhadZjj, WPZhadjj, WMhadZjj, WMZhadjj, ZZhadjj, &
              Hjj_WPhadWM, Hjj_WPWMhad, Hjj_ZZhad, WPhadWPjj, WMhadWMjj,  &
              Sp2jj_WW, Sp2jj_ZZ_ll, Sp2jj_ZZ_lnu, WPAjj, WMAjj, ZAjj, ZAjj_n)
            if (lglobalprint) print*,"Dynamical factorization scale chosen:"
            if (lglobalprint) print*,"      mu_F = ( ET(jj) + ET(VV) )/2"
         CASE DEFAULT
            print *,"unreasonable scale choice for this process:  &
                  id_muf = ",id_muf
            STOP
         END SELECT


      CASE DEFAULT
         print *,"unreasonable scale choice for this process:  &
               id_muf = ",id_muf
      END SELECT

!-------------------------------      
!   renormalization scale
!-------------------------------
      SELECT CASE(id_mur)

! user defined constant scale
      CASE(0)
         qsq = mur_user**2 * xir**2
         als(1,1) = alphas5(qsq,1)
         if (lglobalprint) print*,"Constant value for alphas chosen: alphas(NLO) =",als(1,1)
         call InitPDFs(0)
         als(1,1) = alphas5(qsq,0)
         if (lglobalprint) print*,"Constant value for alphas chosen: alphas(LO) =",als(1,1)
! dynamic scale = momentum transfer of W/Z
      CASE(1)
         SELECT CASE(procID)
         CASE(Hjj, Hjj_AA, Hjj_mu, Hjj_tau, Hjj_bbar, Hjj_WW, Hjj_ZZ_ll, &
              Hjj_ZZ_lnu, HHjj, HHjj_bbtau, HHjj_bbAA, Ajj, Zjj_l, Zjj_nu, WPjj, WMjj, WPWMjj, ZZjj_ll, &
              ZZjj_lnu, WPZjj, WMZjj, HAjj, HAjj_AA, HAjj_mu, HAjj_tau,  &
              HAjj_bbar, HAjj_WW, HAjj_ZZ_ll, HAjj_ZZ_lnu, AAjj, WPWPjj,  &
              WMWMjj, Sp2jj_WW, Sp2jj_ZZ_ll, Sp2jj_ZZ_lnu,  &
              WPhadWMjj, WPWMhadjj, WPhadZjj, WPZhadjj, WMhadZjj, WMZhadjj, ZZhadjj, &
              Hjj_WPhadWM, Hjj_WPWMhad, Hjj_ZZhad, WPhadWPjj, WMhadWMjj, &
              WPAjj, WMAjj, ZAjj, ZAjj_n)
            if (lglobalprint) print*,"Dynamical value for alphas chosen:"
            if (lglobalprint) print*,"alphas = alpha(momentum transfer of W/Z)"   
         CASE DEFAULT
            if (lglobalprint) print *,"unreasonable scale choice for this process:  &
                  id_mur = ",id_mur
            STOP
         END SELECT

! dynamic scale alpha(min(pt_j)) 
      CASE(2)
         SELECT CASE(procID)
         CASE(Hjj, Hjj_AA, Hjj_mu, Hjj_tau, Hjj_bbar, Hjj_WW, Hjj_ZZ_ll, &
              Hjj_ZZ_lnu, HHjj, HHjj_bbtau, HHjj_bbAA, Ajj, Zjj_l, Zjj_nu, WPjj, WMjj, WPWMjj, ZZjj_ll, &
              ZZjj_lnu, Hjjj, Hjjj_AA, Hjjj_mu, Hjjj_tau, Hjjj_bbar, Hjjj_WW,  &
              Hjjj_ZZ_ll,Hjjj_ZZ_lnu, WPZjj, WMZjj, WPAAj, WMAAj, HAjj, HAjj_AA, HAjj_mu,  &
              HAjj_tau, HAjj_bbar, HAjj_WW, HAjj_ZZ_ll, HAjj_ZZ_lnu, AAjj, &
              WPWPjj, WMWMjj, Sp2jj_WW, Sp2jj_ZZ_ll, Sp2jj_ZZ_lnu,  &
              WPhadWMjj, WPWMhadjj, WPhadZjj, WPZhadjj, WMhadZjj, WMZhadjj, ZZhadjj, &
              Hjj_WPhadWM, Hjj_WPWMhad, Hjj_ZZhad, WPhadWPjj, WMhadWMjj, &
              WPAjj, WMAjj, ZAjj,ZAjj_n)


            if (lglobalprint) print*,"Dynamical scale for alphas chosen:"
            if (lglobalprint) print*,"              alphas = alpha(min(pt_j))"
         CASE DEFAULT
            print *,"unreasonable scale choice for this process:  &
                  id_mur = ",id_mur
            STOP
         END SELECT

! dynamic scale = invariant VV mass
      CASE(3,4)
         SELECT CASE(procID)
         CASE(WP_only,WM_only,WPJ,WMJ)
            if (lglobalprint) print*,"Dynamical value for alphas chosen:"
            if (lglobalprint) print*,"alphas = alpha(invariant V mass)"
         CASE( &
              WPH, WPH_AA, WPH_mu, WPH_tau, WPH_bbar, WPH_WW, WPH_ZZ_ll,  WPH_ZZ_lnu, &
              WMH, WMH_AA, WMH_mu, WMH_tau, WMH_bbar, WMH_WW, WMH_ZZ_ll,  WMH_ZZ_lnu, &
              WPHJ, WPHJ_AA, WPHJ_mu, WPHJ_tau, WPHJ_bbar, WPHJ_WW, WPHJ_ZZ_ll,  WPHJ_ZZ_lnu, &
              WMHJ, WMHJ_AA, WMHJ_mu, WMHJ_tau, WMHJ_bbar, WMHJ_WW, WMHJ_ZZ_ll,  WMHJ_ZZ_lnu &
             )
            if (lglobalprint) print*,"Dynamical value for alphas chosen:"
            if (lglobalprint) print*,"alphas = alpha(invariant VH mass)"
         CASE(WW, GFWW, WPZ, WMZ, WPA, WMA, ZZ, GFZZ, ZA, GFZA, AA, GFAA, &
              WPhadWMlep,WPlepWMhad,GFWPhadWMlep,GFWPlepWMhad, &
              ZZhad,GFZZhad, &
              WPhadZ,WPZhad,WMhadZ,WMZhad, &
              WPAJ, WMAJ, WPZJ, WMZJ, &
              WWJ, WPHADWMJ, WPWMHADJ, GFWWj, GFWPHADWMJ, GFWPWMHADJ, &
              ZZJ, GFZZj)
            if (lglobalprint) print*,"Dynamical value for alphas chosen:"
            if (lglobalprint) print*,"alphas = alpha(invariant VV mass)"
         CASE(ZZWP, ZZWM, WWZ, WWWP, WWWM, ZZZ, WWA, ZZA, ZZnA, WPZA, WMZA, &
              WPAA, WMAA, ZAA, ZnAA,AAA, WPAAj, WMAAj, &
              WPhadWMZ, WPWMhadZ, WWZhad, ZZhadWP, ZZWPhad, ZZhadWM,  &
              ZZWMhad, WPhadWMWP, WPWMhadWP, WMhadWPWM, WMWPhadWM, ZZZhad, &
              WPhadWMA, WPWMhadA, ZZhadA, WPhadZA, WPZhadA, WMhadZA, WMZhadA)
            if (lglobalprint) print*,"Dynamical value for alphas chosen:"
            if (lglobalprint) print*,"alphas = alpha(invariant VVV mass)"
         CASE DEFAULT
            print *,"unreasonable scale choice for this process: &
                  id_mur = ",id_mur
            STOP
         END SELECT

! dynamic scale alpha(pt_j1) * alpha(pt_j2) * alpha(m_H)^2
      CASE(5)
         SELECT CASE(procID)
         CASE(GFHjj, GFHjj_AA, GFHjj_mu, GFHjj_tau, GFHjj_bbar, &
              GFHjj_WW, GFHjj_ZZ_ll, GFHjj_ZZ_lnu, WPZjj, WMZjj)
            if (lglobalprint) print*,"Dynamical scale for alphas chosen:"
            if (lglobalprint) print*,"   alphas = alpha(pt_j1)*alpha(pt_j2)* alpha(m_H)^2"
         CASE DEFAULT
            print *,"unreasonable scale choice for this process:  &
                  id_mur = ",id_mur
            STOP
         END SELECT

! constant scale = higgs mass
      CASE(6)
         qsq = xmh**2 * xir**2
         als(1,1) = alphas5(qsq,1)
         if (lglobalprint) print*,"Constant value for alphas chosen: alphas(NLO) =",als(1,1)
         call InitPDFs(0)
         als(1,1) = alphas5(qsq,0)
         if (lglobalprint) print*,"Constant value for alphas chosen: alphas(LO) =",als(1,1)

! dynamic scale = minimal ET
      CASE(7)
         SELECT CASE(procID)
         CASE(WW, GFWW, WWZ, ZZWP, ZZWM, WWWP, WWWM, ZZZ, WWA, ZZA, ZZnA, WPZA,  &
              WMZA, WPAA, WMAA, ZAA, ZnAA, AAA, WPAAj, WMAAj,  &
               WPA,WMA,WPZ,WMZ, ZZ, GFZZ, ZA, GFZA, AA, GFAA, &
              WPhadWMlep,WPlepWMhad,GFWPhadWMlep,GFWPlepWMhad, &
              ZZhad,GFZZhad, &
              WPhadZ,WPZhad,WMhadZ,WMZhad, &
              WPhadWMZ, WPWMhadZ, WWZhad, ZZhadWP, ZZWPhad, ZZhadWM,  &
              ZZWMhad, WPhadWMWP, WPWMhadWP, WMhadWPWM, WMWPhadWM, ZZZhad, &
              WPhadWMA, WPWMhadA, ZZhadA, WPhadZA, WPZhadA, WMhadZA, WMZhadA, &
              WPAJ, WMAJ, WPZJ, WMZJ, WPJ, WMJ, &
              WWJ, WPHADWMJ, WPWMHADJ, GFWWj, GFWPHADWMJ, GFWPWMHADJ, &
              ZZJ, GFZZj)
            if (lglobalprint) print*,"Dynamical value for alphas chosen:"
            if (lglobalprint) print*,"alphas = alpha(min(ET_V))"   
         CASE( &
           WPH, WPH_AA, WPH_mu, WPH_tau, WPH_bbar, WPH_WW, WPH_ZZ_ll,  WPH_ZZ_lnu, &
           WMH, WMH_AA, WMH_mu, WMH_tau, WMH_bbar, WMH_WW, WMH_ZZ_ll,  WMH_ZZ_lnu, &
           WPHJ, WPHJ_AA, WPHJ_mu, WPHJ_tau, WPHJ_bbar, WPHJ_WW, WPHJ_ZZ_ll,  WPHJ_ZZ_lnu, &
           WMHJ, WMHJ_AA, WMHJ_mu, WMHJ_tau, WMHJ_bbar, WMHJ_WW, WMHJ_ZZ_ll,  WMHJ_ZZ_lnu)
            if (lglobalprint) print*,"Dynamical value for alphas chosen:"
            if (lglobalprint) print*,"alphas = alpha(min(ET_V, ET_H))"   
         CASE DEFAULT
            print *,"unreasonable scale choice for this process:  &
                  id_mur = ",id_mur
            STOP
         END SELECT

! dynamical scale HT/2,  HT = sum(pT partons) + sum( ET V )
      CASE(8)
         SELECT CASE(procID)
         CASE(QCDWMWMjj,QCDWPWPjj,QCDWPjj,QCDWMjj,QCDWPZjj,QCDWMZjj,QCDAAjj, &
              QCDZjj_l, QCDZjj_nu,QCDWPAjj,QCDWMAjj,QCDZZjj_ll,QCDZZjj_lnu, &
              QCDZAjj_l,QCDZAjj_n, &
              WPZ,WMZ, WPZJ, WMZJ, &
              WPA,WPAJ,WMA,WMAJ, &
              WW,WWJ,GFWW,GFWWj, &
              ZZ,ZZJ,GFZZ,GFZZj, &
              ZA,GFZA, &
              AA,GFAA, &
              Hjj, Hjj_AA, Hjj_mu, Hjj_tau, Hjj_bbar, Hjj_WW, Hjj_ZZ_ll, &
              Hjj_ZZ_lnu, HHjj, HHjj_bbtau, HHjj_bbAA, Ajj, Zjj_l, Zjj_nu, WPjj, WMjj, WPWMjj, ZZjj_ll, &
              ZZjj_lnu, Hjjj, Hjjj_AA, Hjjj_mu, Hjjj_tau, Hjjj_bbar, Hjjj_WW,  &
              Hjjj_ZZ_ll,Hjjj_ZZ_lnu, WPZjj, WMZjj,WPAAj, WMAAj, HAjj, HAjj_AA, HAjj_mu, &
              HAjj_tau, HAjj_bbar, HAjj_WW, HAjj_ZZ_ll, HAjj_ZZ_lnu, &
              WPWPjj,WMWMjj,AAjj, &
              WPhadWMjj, WPWMhadjj, WPhadZjj, WPZhadjj, WMhadZjj, WMZhadjj, ZZhadjj, &
              Hjj_WPhadWM, Hjj_WPWMhad, Hjj_ZZhad, WPhadWPjj, WMhadWMjj,  &
              Sp2jj_WW, Sp2jj_ZZ_ll, Sp2jj_ZZ_lnu, WPAjj, WMAjj, ZAjj, ZAjj_n)
            if (lglobalprint) print*,"Dynamical renormalization scale chosen:"
            if (lglobalprint) print*,"      mu_R = ( sum(pT partons) + sum( ET V ) )/2 "
         CASE DEFAULT
            print *,"unreasonable scale choice for this process:  &
                  id_mur = ",id_mur
            STOP
         END SELECT

! dynamical scale sum(pT_j*exp(|y_j-y12|)) + sum( ET V ) ( see arXiv 1311.6738 )
      CASE(9)
         SELECT CASE(procID)
         CASE(QCDWMWMjj,QCDWPWPjj,QCDWPjj,QCDWMjj,QCDWPZjj,QCDWMZjj,QCDAAjj, &
              QCDZjj_l, QCDZjj_nu,QCDWPAjj,QCDWMAjj,QCDZZjj_ll,QCDZZjj_lnu,QCDZAjj_l,QCDZAjj_n, &
              Hjj, Hjj_AA, Hjj_mu, Hjj_tau, Hjj_bbar, Hjj_WW, Hjj_ZZ_ll, &
              Hjj_ZZ_lnu, HHjj, HHjj_bbtau, HHjj_bbAA, Ajj, Zjj_l, Zjj_nu, WPjj, WMjj, WPWMjj, ZZjj_ll, &
              ZZjj_lnu, Hjjj, Hjjj_AA, Hjjj_mu, Hjjj_tau, Hjjj_bbar, Hjjj_WW,  &
              Hjjj_ZZ_ll,Hjjj_ZZ_lnu, WPZjj, WMZjj,WPAAj, WMAAj, HAjj, HAjj_AA, HAjj_mu, &
              HAjj_tau, HAjj_bbar, HAjj_WW, HAjj_ZZ_ll, HAjj_ZZ_lnu, &
              WPWPjj,WMWMjj,AAjj, &
              WPhadWMjj, WPWMhadjj, WPhadZjj, WPZhadjj, WMhadZjj, WMZhadjj, ZZhadjj, &
              Hjj_WPhadWM, Hjj_WPWMhad, Hjj_ZZhad, WPhadWPjj, WMhadWMjj,  &
              Sp2jj_WW, Sp2jj_ZZ_ll, Sp2jj_ZZ_lnu, WPAjj, WMAjj, ZAjj, ZAjj_n)

            if (lglobalprint) print*,"Dynamical renormalization scale chosen:"
            if (lglobalprint) print*,"      mu_R = ( sum(pT_j*exp(|y_j-y12|)) + sum( ET V ) )/2"
            if (lglobalprint) print*,"( see arXiv 1311.6738 ) "
         CASE DEFAULT
            print *,"unreasonable scale choice for this process:  &
                  id_mur = ",id_mur
            STOP
         END SELECT

! dynamical scale HT/2,  HT = sum(pT partons) + sum( ET V )
      CASE(10)
         SELECT CASE(procID)
         CASE(QCDWMWMjj,QCDWPWPjj,QCDWPjj,QCDWMjj,QCDWPZjj,QCDWMZjj,QCDAAjj, &
              QCDZjj_l, QCDZjj_nu,QCDWPAjj,QCDWMAjj,QCDZZjj_ll,QCDZZjj_lnu,QCDZAjj_l,QCDZAjj_n, &
              Hjj, Hjj_AA, Hjj_mu, Hjj_tau, Hjj_bbar, Hjj_WW, Hjj_ZZ_ll, &
              Hjj_ZZ_lnu, HHjj, HHjj_bbtau, HHjj_bbAA, Ajj, Zjj_l, Zjj_nu, WPjj, WMjj, WPWMjj, ZZjj_ll, &
              ZZjj_lnu, Hjjj, Hjjj_AA, Hjjj_mu, Hjjj_tau, Hjjj_bbar, Hjjj_WW,  &
              Hjjj_ZZ_ll,Hjjj_ZZ_lnu, WPZjj, WMZjj,WPAAj, WMAAj, HAjj, HAjj_AA, HAjj_mu, &
              HAjj_tau, HAjj_bbar, HAjj_WW, HAjj_ZZ_ll, HAjj_ZZ_lnu, &
              WPWPjj,WMWMjj,AAjj, &
              WPhadWMjj, WPWMhadjj, WPhadZjj, WPZhadjj, WMhadZjj, WMZhadjj, ZZhadjj, &
              Hjj_WPhadWM, Hjj_WPWMhad, Hjj_ZZhad, WPhadWPjj, WMhadWMjj,  &
              Sp2jj_WW, Sp2jj_ZZ_ll, Sp2jj_ZZ_lnu, WPAjj, WMAjj, ZAjj, ZAjj_n)

            if (lglobalprint) print*,"Dynamical renormalization scale chosen:"
            if (lglobalprint) print*,"      mu_R = ( ET(jj) + ET(VV) )/2"
         CASE DEFAULT
            print *,"unreasonable scale choice for this process:  &
                  id_mur = ",id_mur
            STOP
         END SELECT


      CASE DEFAULT
         print *,"unreasonable scale choice for this process:  &
                 id_mur = ",id_mur
         STOP
      END SELECT


      SELECT CASE(id_mur)
      CASE(0,6)
         if (lglobalprint) print*," "
      CASE DEFAULT
         call InitPDFs(0)
         if (lglobalprint) print*," "
      END SELECT


!* Setting LO alfas(MZ)
      alfas_lo = alphas5(xmz**2D0,0)

      else 
! just some useful defaults
        alfas_lo = alphas5_hardwired(91.1876d0**2D0,0)
      endif

      RETURN
      END

!*************************************************************************
      Subroutine Scales(p, v, jets, njets, leptons, nl, photons,  &
                                              nphotons, Lmax, nlo)
!*************************************************************************
!     This is the general scales routine. It is intended to be the interface
!     to any basic, specialized or user-defined scales routines. 
!*************************************************************************
!     INPUT 
!     p       : 4-momenta of the partons
!     v       : 4-momenta of all other particles (leptons, photons, ...)
!     jets    : the array with the jet information as it is returned by the 
!               jetdefinition routine.
!     njets   : number of jets
!     leptons : the array with the charged lepton information.
!     nl      : number of charged leptons
!     photons : the array with the photon information.
!     nphotons: number of photons
!     Lmax    : number of different kinematics needed
!     nlo     : LO( = 0) or NLO(= 1)
!
!     OUTPUT
!     factorization scale and alfas(ren. scale) is transferred via 
!     common block in scales.inc
!*************************************************************************

      implicit none

#include "scales.inc"
#include "global.inc"
#include "koppln.inc"
#include "process.inc"

!     input variables
      real*8 p(0:3,max_p,max_kin)
      real*8 v(0:3,max_v,max_kin)
      real*8 jets(0:7,max_jets,max_kin) 
      real*8 leptons(0:8,max_v,max_kin)     
      real*8 photons(0:7,max_v,max_kin)     
      integer Lmax, njets(max_kin), nl(max_kin), nphotons(max_kin), nlo


      SELECT CASE(procID)

! VBF: Hjj, Vjj, VVjj
      CASE(Hjj, Hjj_AA, Hjj_mu, Hjj_tau, Hjj_bbar, Hjj_WW, Hjj_ZZ_ll, &
           Hjj_ZZ_lnu, Ajj, Zjj_l, Zjj_nu, WPjj, WMjj, WPWMjj, ZZjj_ll, &
           ZZjj_lnu, WPZjj, WMZjj, HAjj, HAjj_AA, HAjj_mu, HAjj_tau,  &
           HAjj_bbar, HAjj_WW, HAjj_ZZ_ll, HAjj_ZZ_lnu, AAjj,WPWPjj, &
           WMWMjj, HHjj, HHjj_bbtau, HHjj_bbAA, Sp2jj_WW, Sp2jj_ZZ_ll, Sp2jj_ZZ_lnu,  &
           WPhadWMjj, WPWMhadjj, WPhadZjj, WPZhadjj, WMhadZjj, WMZhadjj, ZZhadjj, &
           Hjj_WPhadWM, Hjj_WPWMhad, Hjj_ZZhad, WPhadWPjj, WMhadWMjj, &
           WPAjj, WMAjj, ZAjj, ZAjj_n)
          call calcVBFscales(p, v, jets, njets, leptons, nl, photons,  &
                                                  nphotons, 3, nlo)

! VBF: Hjjj
      CASE(Hjjj,Hjjj_AA, Hjjj_mu, Hjjj_tau, Hjjj_bbar, Hjjj_WW,  &
               Hjjj_ZZ_ll,Hjjj_ZZ_lnu)
          call calcH3Jscales(p, v, jets, njets, leptons, nl, photons,  &
                                                  nphotons, Lmax, nlo)

! DIBOSON and DIBOSON PLUS JET
! including WH, WHj and W, Wj          
      CASE(WW, GFWW, WPA, WMA, WPZ, WMZ,  &
           ZZ, GFZZ, ZA, GFZA, AA, GFAA, &
           WMAJ,WPAJ,WMZJ,WPZJ,WWJ,GFWWj,ZZJ,GFZZj, &
           WPhadWMlep,WPlepWMhad,GFWPhadWMlep,GFWPlepWMhad, &
           ZZhad,GFZZhad, &
           WPhadZ,WPZhad,WMhadZ,WMZhad, &
           WPHADWMJ, WPWMHADJ, GFWPHADWMJ, GFWPWMHADJ, &
           WPhadZJ, WPZhadJ, WMhadZJ, WMZhadJ, &
           WP_only, WM_only, WPJ, WMJ, &
           WPH, WPH_AA, WPH_mu, WPH_tau, WPH_bbar, WPH_WW, WPH_ZZ_ll,  WPH_ZZ_lnu, &
           WMH, WMH_AA, WMH_mu, WMH_tau, WMH_bbar, WMH_WW, WMH_ZZ_ll,  WMH_ZZ_lnu, &
           WPHJ, WPHJ_AA, WPHJ_mu, WPHJ_tau, WPHJ_bbar, WPHJ_WW, WPHJ_ZZ_ll,  WPHJ_ZZ_lnu, &
           WMHJ, WMHJ_AA, WMHJ_mu, WMHJ_tau, WMHJ_bbar, WMHJ_WW, WMHJ_ZZ_ll,  WMHJ_ZZ_lnu)
          call calcVXjScales(p, v, jets, njets, leptons, nl, photons, nphotons, nlo, Lmax)

! TRIBOSON
      CASE(WWZ, ZZZ, ZZWP, ZZWM, WWWP, WWWM, WWA, ZZA, ZZnA, WPZA, WMZA, &
              WPAA, WMAA, ZAA, ZnAA, AAA, &
              WPhadWMZ, WPWMhadZ, WWZhad, ZZhadWP, ZZWPhad, ZZhadWM,  &
              ZZWMhad, WPhadWMWP, WPWMhadWP, WMhadWPWM, WMWPhadWM, ZZZhad, &
              WPhadWMA, WPWMhadA, ZZhadA, WPhadZA, WPZhadA, WMhadZA, WMZhadA)
          call calcVXjScales(p, v, jets, njets, leptons, nl, photons, nphotons, nlo, Lmax)

! TRIBOSON plus jet
      CASE(WPAAj, WMAAj )
          call calcVXjScales(p, v, jets, njets, leptons, nl, photons, nphotons, nlo, Lmax)

! Gluon fusion
      CASE(GFHjj, GFHjj_AA, GFHjj_mu, GFHjj_tau, GFHjj_bbar, &
               GFHjj_WW, GFHjj_ZZ_ll, GFHjj_ZZ_lnu)
         call calcGGFScales(p)

! QCD induced processes
      CASE(QCDWMWMjj,QCDWPWPjj,QCDWPjj,QCDWMjj,QCDWPZjj,QCDWMZjj,QCDAAjj, &
       QCDZjj_l, QCDZjj_nu,QCDWPAjj,QCDWMAjj,QCDZZjj_ll,QCDZZjj_lnu,QCDZAjj_l,QCDZAjj_n)
         call calcQCDscales(p, v, jets, njets, leptons, nl, photons,  &
                                                  nphotons, Lmax, nlo)

      CASE DEFAULT
         print *,"Error: scales routine for process "// &
                  "not implemented, ID = ",procID
         STOP
      END SELECT

      RETURN
      END


!***************************************************************************
      SUBROUTINE calcVBFscales(p, v, jets, njets, leptons, nl, photons,  &
                                              nphotons, Lmax, nlo)
!***************************************************************************
!     INPUT 
!     p       : 4-momenta of the partons
!     v       : 4-momenta of all other particles (leptons, photons, ...)
!     jets    : the array with the jet information as it is returned by the 
!               jetdefinition routine.
!     njets   : number of jets
!     leptons : the array with the charged lepton information.
!     nl      : number of charged leptons
!     photons : the array with the photon information.
!     nphotons: number of photons
!     Lmax    : number of different kinematics needed
!     nlo     : LO( = 0) or NLO(= 1)
!
!     OUTPUT
!     factorization scale and alfas(ren. scale) is transferred via 
!     common block in scales.inc
!***************************************************************************
      implicit none

#include "scales.inc"
#include "global.inc"
#include "koppln.inc"

!     input variables
      real*8 p(0:3,max_p,max_kin)
      real*8 v(0:3,max_v,max_kin)
      real*8 jets(0:7,max_jets,max_kin) 
      real*8 leptons(0:8,max_v,max_kin)     
      real*8 photons(0:7,max_v,max_kin)     
      integer Lmax, njets(max_kin), nl(max_kin), nphotons(max_kin), nlo

      real*8 q12(0:4,max_kin), q34(0:4,max_kin), qsq, pt(max_jets), exp
      integer L, i
      real*8 alphas5,mass2,temp,transverseEnergy,tempvec(0:3),y12
      external alphas5,mass2,transverseEnergy

!----------------------------      
!   factorization scale
!----------------------------      
      SELECT CASE(id_muf)
! user defined constant scale 
      CASE(0)
         qsq = muf_user**2 * xif**2
         do L = 1,Lmax
            mufsq(1,L) = qsq
            mufsq(2,L) = qsq
         enddo

! dynamic scale = momentum transfer of W/Z
      CASE(1)
         call Calc_Momentum_Transfer(p, v, q12, q34, Lmax)
         do L = 1,Lmax
            mufsq(1,L) = max(xif**2*q12(4,L),2d0)
            mufsq(2,L) = max(xif**2*q34(4,L),2d0)
         enddo

! dynamic scale = min(pt_j1, pt_j2, ...)
      CASE(2)
         do L = 1,Lmax
         qsq = 1.0d12
         do i = 1,njets(L)
            pt(i) = jets(5,i,L)**2
            pt(i) = MAX(pt(i),4D0)
            qsq = min(qsq,pt(i))
         enddo
            qsq = qsq * xif**2         
            mufsq(1,L) = qsq
            mufsq(2,L) = qsq
         enddo

! constant scale = higgs mass
      CASE(6)
         qsq = xmh**2 * xif**2
         do L = 1,Lmax
            mufsq(1,L) = qsq
            mufsq(2,L) = qsq
         enddo

         
! HT/2 HT = sum(pT partons) + sum( ET V )
      CASE(8)
          do L=1,Lmax
            temp=0d0
            do i = 1,n_p ! sum over parton pt
              temp = temp + sqrt(p(1,i,L)*p(1,i,L) + p(2,i,L)*p(2,i,L))
            enddo
            do i = 1,nphotons(L) ! sum over photon pt
              temp = temp + photons(5,i,L)
            enddo
            do i = 1,n_v-nphotons(L),2 ! sum over lepton pairs
              tempvec(:)= v(:,i,L) + v(:,i+1,L)
              temp = temp + transverseEnergy(tempvec)
            enddo
            
            temp=temp/2d0*xif
            temp = temp**2
            mufsq(1,L) = temp
            mufsq(2,L) = temp
         enddo

! sum(pT_j*exp(|y_j-y12|)) + sum( ET V ) ( see arXiv 1311.6738 )
      CASE(9)
          do L=1,Lmax
            temp=0d0
            y12 = 0.5d0*(jets(6,1,L)+jets(6,2,L)) ! avg rapidity of two hardest jets
            do i = 1,njets(L) ! sum over jets
              temp = temp + jets(5,i,L)*exp(abs(jets(6,i,L)-y12))
            enddo
            do i = 1,nphotons(L) ! sum over photon pt
              temp = temp + photons(5,i,L)
            enddo
            do i = 1,n_v-nphotons(L),2 ! sum over lepton pairs
              tempvec(:)= v(:,i,L) + v(:,i+1,L)
              temp = temp + transverseEnergy(tempvec)
            enddo
            
            temp=temp/2d0*xif
            temp = temp**2
            mufsq(1,L) = temp
            mufsq(2,L) = temp
         enddo


! E_T(jj) + E_T(VV)
      CASE(10)
          do L=1,Lmax
            tempvec = jets(0:3,1,L) + jets(0:3,2,L)
            temp = transverseEnergy(tempvec) 
            
            tempvec=0d0
            do i=1,n_v 
              tempvec = tempvec + v(0:3,i,L)
            enddo
            temp = temp + transverseEnergy(tempvec)
            
            temp=temp/2d0*xif
            temp = temp**2
            mufsq(1,L) = temp
            mufsq(2,L) = temp
         enddo


      CASE DEFAULT
         print *,"unreasonable scale choice : id_muf = ",id_muf
         STOP
      END SELECT

!----------------------------
!   renormalization scale
!----------------------------
      SELECT CASE(id_mur)

! user defined constant scale
      CASE(0)
         qsq = mur_user**2 * xir**2
         if (nlo.eq.0) then
            als(1,1) = alphas5(qsq,0)
         else
            als(1,1) = alphas5(qsq,1)
         endif
         do L = 1,Lmax
            als(1,L) = als(1,1)
            als(2,L) = als(1,1)
         enddo

! dynamic scale = momentum transfer of W/Z
      CASE(1)
         call Calc_Momentum_Transfer(p, v, q12, q34, Lmax)
         do L = 1,Lmax
            qsq = max( xir**2*q12(4,L), 2d0 )
            if (nlo.eq.0) then
               als(1,L) = alphas5( qsq, 0 )
            else
               als(1,L) = alphas5( qsq, 1 )
            endif
            qsq = max( xir**2*q34(4,L), 2d0 )
            if (nlo.eq.0) then
               als(2,L) = alphas5( qsq, 0 )
            else
               als(2,L) = alphas5( qsq, 1 )
            endif
         enddo

! dynamic scale alpha(min(pt_j)) 
      CASE(2)
         do L = 1,Lmax
         qsq = 1.0d12
         do i = 1,njets(L)
            pt(i) = jets(5,i,L)**2
            pt(i) = MAX(pt(i),4D0)
            qsq = min(qsq,pt(i))
         enddo
         qsq = qsq * xir**2         
         if (nlo.eq.0) then
            als(1,L) = alphas5(qsq,0)
            als(2,L) = alphas5(qsq,0)
         else
            als(1,L) = alphas5(qsq,1)
            als(2,L) = alphas5(qsq,1)
         endif
         enddo

! constant scale = higgs mass
      CASE(6)
         qsq = xmh**2 * xir**2
         if (nlo.eq.0) then
            als(1,1) = alphas5(qsq,0)
         else
            als(1,1) = alphas5(qsq,1)
         endif
         do L = 1,Lmax
            als(1,L) = als(1,1)
            als(2,L) = als(1,1)
         enddo
         
! HT/2 HT = sum(pT partons) + sum( ET V )
      CASE(8)
          do L=1,Lmax
           if(Njets(L).ge.2) then
            temp=0d0
            do i = 1,n_p ! sum over parton pt
              temp = temp + sqrt(p(1,i,L)*p(1,i,L) + p(2,i,L)*p(2,i,L))
            enddo
            do i = 1,nphotons(L) ! sum over photon pt
              temp = temp + photons(5,i,L)
            enddo
            do i = 1,n_v-nphotons(L),2 ! sum over lepton pairs
              tempvec(:)= v(:,i,L) + v(:,i+1,L)
              temp = temp + transverseEnergy(tempvec)
            enddo
            
            temp=temp/2d0*xir
            temp = temp**2
            mursq(1,L) = temp
            if(nlo.eq.0) then
              als(1,L) = alphas5(mursq(1,L),0)
            else
              als(1,L) = alphas5(mursq(1,L),1)
            endif
           else
            als(1,L) = 0d0
           endif
           als(2,L) = als(1,L)
         enddo

! sum(pT_j*exp(|y_j-y12|)) + sum( ET V ) ( see arXiv 1311.6738 )
      CASE(9)
          do L=1,Lmax
           if(Njets(L).ge.2) then
            temp=0d0
            y12 = 0.5d0*(jets(6,1,L)+jets(6,2,L)) ! avg rapidity of two hardest jets
            do i = 1,njets(L) ! sum over jets
              temp = temp + jets(5,i,L)*exp(abs(jets(6,i,L)-y12))
            enddo
            do i = 1,nphotons(L) ! sum over photon pt
              temp = temp + photons(5,i,L)
            enddo
            do i = 1,n_v-nphotons(L),2 ! sum over lepton pairs
              tempvec(:)= v(:,i,L) + v(:,i+1,L)
              temp = temp + transverseEnergy(tempvec)
            enddo
            
            temp=temp/2d0*xir
            temp = temp**2
            mursq(1,L) = temp
            if(nlo.eq.0) then
              als(1,L) = alphas5(mursq(1,L),0)
            else
              als(1,L) = alphas5(mursq(1,L),1)
            endif
           else
            als(1,L) = 0d0
           endif
           als(2,L) = als(1,L)
         enddo


! E_T(jj) + E_T(VV)
      CASE(10)
          do L=1,Lmax
           if(Njets(L).ge.2) then
            tempvec = jets(0:3,1,L) + jets(0:3,2,L)
            temp = transverseEnergy(tempvec) 
            
            tempvec=0d0
            do i=1,n_v 
              tempvec = tempvec + v(0:3,i,L)
            enddo
            temp = temp + transverseEnergy(tempvec)
            
            temp=temp/2d0*xir
            temp = temp**2
            mursq(1,L) = temp
            if(nlo.eq.0) then
              als(1,L) = alphas5(mursq(1,L),0)
            else
              als(1,L) = alphas5(mursq(1,L),1)
            endif
           else
            als(1,L) = 0d0
           endif
           als(2,L) = als(1,L)
         enddo


      CASE DEFAULT
         print *,"unreasonable scale choice : id_mur = ",id_mur
         STOP
      END SELECT

      do L = 1, Lmax
         mursq(1,L) = qsq
         mursq(2,L) = qsq
      enddo

      RETURN
      END


!***************************************************************************
      SUBROUTINE calcH3jscales(p, v, jets, njets, leptons, nl, photons,  &
                                              nphotons, Lmax, nlo)
!***************************************************************************
!     INPUT 
!     p       : 4-momenta of the partons
!     v       : 4-momenta of all other particles (leptons, photons, ...)
!     jets    : the array with the jet information as it is returned by the 
!               jetdefinition routine.
!     njets   : number of jets
!     leptons : the array with the charged lepton information.
!     nl      : number of charged leptons
!     photons : the array with the photon information.
!     nphotons: number of photons
!     Lmax    : number of different kinematics needed
!     nlo     : LO( = 0) or NLO(= 1)
!
!     OUTPUT
!     factorization scale and alfas(ren. scale) is transferred via 
!     common block in scales.inc
!***************************************************************************
      implicit none

#include "scales.inc"
#include "global.inc"
#include "koppln.inc"

!     input variables
      real*8 p(0:3,max_p,max_kin)
      real*8 v(0:3,max_v,max_kin)
      real*8 jets(0:7,max_jets,max_kin) 
      real*8 leptons(0:8,max_v,max_kin)     
      real*8 photons(0:7,max_v,max_kin)     
      integer Lmax, njets(max_kin), nl(max_kin), nphotons(max_kin), nlo

      real*8 qsq, pt(max_jets)
      integer L, i
      real*8 alphas5
      external alphas5

!----------------------------      
!   factorization scale
!----------------------------      
      SELECT CASE(id_muf)

! user defined constant scale 
      CASE(0)
         qsq = muf_user**2 * xif**2
         do L = 1,Lmax
            mufsq(1,L) = qsq
            mufsq(2,L) = qsq
         enddo

! dynamic scale = min(pt_j1, pt_j2, ...)
      CASE(2)
         do L = 1,Lmax
           qsq = 1.0d12
            do i = 1,njets(L)
               pt(i) = jets(5,i,L)**2
               pt(i) = MAX(pt(i),4D0)
               qsq = min(qsq,pt(i))
            enddo
            qsq = qsq * xif**2         
            mufsq(1,L) = qsq
            mufsq(2,L) = qsq
         enddo

! constant scale = higgs mass
      CASE(6)
         qsq = xmh**2 * xif**2
         do L = 1,Lmax
            mufsq(1,L) = qsq
            mufsq(2,L) = qsq
         enddo
      CASE DEFAULT
         print *,"unreasonable scale choice : id_muf = ",id_muf
         STOP
      END SELECT

!----------------------------
!   renormalization scale
!----------------------------
      SELECT CASE(id_mur)

! user defined constant scale
      CASE(0)
         qsq = mur_user**2 * xir**2
         if (nlo.eq.0) then
            als(1,1) = alphas5(qsq,0)
         else
            als(1,1) = alphas5(qsq,1)
         endif
         do L = 1,Lmax
            als(1,L) = als(1,1)
            als(2,L) = als(1,1)
            mursq(1,L) = qsq
            mursq(2,L) = qsq
         enddo


! dynamic scale alpha(min(pt_j)) 
      CASE(2)
         do L = 1,Lmax
            qsq = 1.0d12
            do i = 1,njets(L)
               pt(i) = jets(5,i,L)**2
               pt(i) = MAX(pt(i),4D0)
               qsq = min(qsq,pt(i))
            enddo
            qsq = qsq * xir**2
            mursq(1,L) = qsq
            mursq(2,L) = qsq
            if (nlo.eq.0) then
               als(1,L) = alphas5(qsq,0)
               als(2,L) = alphas5(qsq,0)
            else
               als(1,L) = alphas5(qsq,1)
               als(2,L) = alphas5(qsq,1)
            endif
         enddo

! constant scale = higgs mass
      CASE(6)
         qsq = xmh**2 * xir**2
         if (nlo.eq.0) then
            als(1,1) = alphas5(qsq,0)
         else
            als(1,1) = alphas5(qsq,1)
         endif
         do L = 1,Lmax
            als(1,L) = als(1,1)
            als(2,L) = als(1,1)
            mursq(1,L) = qsq
            mursq(2,L) = qsq
         enddo
      CASE DEFAULT
         print *,"unreasonable scale choice : id_mur = ",id_mur
         STOP
      END SELECT


      RETURN
      END


!***************************************************************************
      SUBROUTINE calcGGFScales(p)
!***************************************************************************
!     Calculates factorization and renormalization scales for 
!     single top processes.
!***************************************************************************
      implicit none

#include "scales.inc"
#include "global.inc"
#include "koppln.inc"

      real*8 p(0:3,max_p)
      integer l
      real*8 alphas5, qsq, qsq2, pt1, pt2
      external alphas5

!* clrCT is the counterterm for qqV
      double complex clrCT
      common /bkopouCT/ clrCT(3:4,2:4,-1:1)
      DOUBLE PRECISION CLR,XM2,XMG,Bk,Vk,Ak
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),Bk(6,6,6), &
           Vk(4,5),Ak(4,5)

     
!----------------------------      
!   factorization scale
!----------------------------      
      SELECT CASE(id_muf)

! user defined constant scale 
      CASE(0)
         qsq = muf_user**2 * xif**2
         do L = 1,3
            mufsq(1,L) = qsq
            mufsq(2,L) = qsq
         enddo

! dynamic scale = sqrt(pt_j1*pt_j2)
      CASE(5)
         qsq = p(1,3)**2+p(2,3)**2
         qsq2 = p(1,4)**2+p(2,4)**2
         qsq = MAX(qsq,4D0)
         qsq2 = MAX(qsq2,4D0)
         qsq = sqrt(qsq*qsq2)
         qsq = MAX(qsq,5D0)* xif**2
          do L = 1,3
            mufsq(1,L) = qsq
            mufsq(2,L) = qsq
         enddo

! constant scale = higgs mass
      CASE(6)
         qsq = xmh**2 * xif**2
         do L = 1,3
            mufsq(1,L) = qsq
            mufsq(2,L) = qsq
         enddo

      CASE DEFAULT
         print *,"unreasonable scale choice : id_muf = ",id_muf
         STOP
      END SELECT

!----------------------------
!   renormalization scale
!----------------------------
      SELECT CASE(id_mur)

! user defined constant scale
      CASE(0)
         qsq = mur_user**2 * xir**2
         als(1,1) = alphas5(qsq,0)
         do L = 1,3
            als(1,L) = als(1,1)
            als(2,L) = als(1,1)
            als(3,L) = als(1,1)
         enddo

! dynamic scale alpha(pt_j1) * alpha(pt_j2) * alpha(m_H)^2
      CASE(5)
         PT1=p(1,3)**2+p(2,3)**2
         PT2=p(1,4)**2+p(2,4)**2
         PT1=MAX(PT1,4D0)* xir**2
         PT2=MAX(PT2,4D0)* xir**2
         qsq =  xm2(6) * xir**2
         do L = 1,3
            als(1,L) = alphas5(pt1,0)
            als(2,L) = alphas5(pt2,0)
            als(3,L) = alphas5(qsq,0)
         enddo

! constant scale = higgs mass
      CASE(6)
         qsq = xmh**2 * xir**2
         als(1,1) = alphas5(qsq,0)
         do L = 1,3
            als(1,L) = als(1,1)
            als(2,L) = als(1,1)
            als(3,L) = als(1,1)
         enddo
       
      CASE DEFAULT
         print *,"unreasonable scale choice : id_mur = ",id_mur
         STOP
      END SELECT


      RETURN
      END


!***************************************************************************
      SUBROUTINE CalcVXjScales(p, v, jets, njets, leptons, nl, photons, nphotons, nlo, Lmax)
!***************************************************************************
!  Calculates factorization and renormalization scales for diboson processes 
!  including W, Wj, WH, WHj.
!
!  INPUT
!     p    : 4-momenta of the partons involved in the basic process.
!            The first 2 entries are the incoming partons.
!     v    : The 4-momenta of the other particles
!     jets    : the array with the jet information as it is returned by the 
!               jetdefinition routine.
!     njets   : number of jets
!     leptons : the array with the charged lepton information.
!     nl      : number of charged leptons
!     photons : the array with the photon information.
!     nphotons: number of photons
!     Lmax    : number of different kinematics needed
!     nlo     : LO( = 0) or NLO(= 1)
!     Lmax : number of momentum configurations, L=1 for direct term
!            L=2,...,Lmax for ptilde momenta for subtraction terms
!     
!  OUTPUT
!     mufsq, mursq, als via common block (defined in scales.inc)
!***************************************************************************
      implicit none

#include "scales.inc"
#include "global.inc"
#include "coupl.inc"
#include "process.inc"

      integer nlo, i
      real*8 p(0:3,max_p,max_kin), v(0:3,max_v,max_kin)
      real*8 jets(0:7,max_jets,max_kin) 
      real*8 leptons(0:8,max_v,max_kin)     
      real*8 photons(0:7,max_v,max_kin)     
      integer njets(max_kin), nl(max_kin), nphotons(max_kin)
      real*8 qsq(1:lmax)

      integer Lmax

      real*8 alphas5
      external alphas5


!     CalcVXjScales supports only a common scale for both protons and all occurring alpha_s!


!----------------------------      
!   factorization scale
!----------------------------      

      call CalcVXjScale_helper(p, v, jets, njets, leptons, nl, photons, nphotons, nlo, lmax, xif, id_muf, muf_user, .false., qsq)

      mufsq(1,1:lmax) = qsq(1:lmax)
      mufsq(2,1:lmax) = qsq(1:lmax)


!----------------------------      
!   renormalization scale
!----------------------------      
      call CalcVXjScale_helper(p, v, jets, njets, leptons, nl, photons, nphotons, nlo, lmax, xir, id_mur, mur_user, .true., qsq)

      mursq(1,1:lmax) = qsq(1:lmax)
      mursq(2,1:lmax) = qsq(1:lmax)
      if (nlo.eq.0) then
         do i = 1, lmax
            als(1,i) = alphas5(qsq(i),0)
            als(2,i) = alphas5(qsq(i),0)
         enddo
      else
         do i = 1, lmax
            als(1,i) = alphas5(qsq(i),1)
            als(2,i) = alphas5(qsq(i),1)
         enddo
      endif        

      ! if the chosen scale is too small, alphas5 breaks and returns 'nan'
      ! this happens only in extreme cases so neglecting those
      if(isnan(als(1,1))) then
         print*, 'als nan in scales.F calcvxjscales'
         print*, 'setting als=0'
         als(1,1) = 0d0
      endif
        
      return

      end


      subroutine CalcVXjScale_helper(p, v, jets, njets, leptons, nl, photons, nphotons, &
                                     nlo, Lmax, xi, scaleid, mu_user, isRenScale, qsq)

      implicit none

#include "scales.inc"
#include "global.inc"
#include "coupl.inc"
#include "process.inc"

      real*8, intent(in) :: p(0:3,max_p,max_kin), v(0:3,max_v,max_kin)
      real*8, intent(in) :: jets(0:7,max_jets,max_kin) 
      real*8, intent(in) :: leptons(0:8,max_v,max_kin)     
      real*8, intent(in) :: photons(0:7,max_v,max_kin)     
      integer, intent(in) :: njets(max_kin), nl(max_kin), nphotons(max_kin)
      integer, intent(in) :: nlo, scaleid, lmax
      logical, intent(in) :: isRenScale         ! this can be used for scale choices which have different definitions
                                                ! for renormalization and factorization scale
      real*8, intent(in) :: mu_user, xi

      real*8, intent(out) :: qsq(max_kin)

      real*8 qvv(0:4,max_kin)
      real*8 qv(0:3)
      real*8 ptsq
      double precision etmin, transverseEnergy
      integer L, mu, i
      real*8 alphas5, mjj2, mass2
      external alphas5, mjj2, mass2, transverseEnergy

      logical ldebug
      parameter (ldebug=.false.)

      integer n_Vmassive, n_photons
      logical higgs_process, atLeastOneJet

      ! to get the masses
      real*8 clr(4,5,-1:1), xm2(6), xmg(6), Bk(6,6,6), Vk(4,5), Ak(4,5)
      COMMON /BKOPOU/ CLR, XM2, XMG, Bk, Vk, Ak


! Get number of massive vector bosons (W,Z) and massless bosons (photons)
! from n_lepton+n_decayquark and particle_IDs(n_lepton+n_decayquark+i).
! In case of the VH(j)-processes all particles from 3..n_v are counted towards the higgs.
      SELECT CASE (procID)
      CASE ( &
           WPH, WPH_AA, WPH_mu, WPH_tau, WPH_bbar, WPH_WW, WPH_ZZ_ll,  WPH_ZZ_lnu, &
           WMH, WMH_AA, WMH_mu, WMH_tau, WMH_bbar, WMH_WW, WMH_ZZ_ll,  WMH_ZZ_lnu, &
           WPHJ, WPHJ_AA, WPHJ_mu, WPHJ_tau, WPHJ_bbar, WPHJ_WW, WPHJ_ZZ_ll,  WPHJ_ZZ_lnu, &
           WMHJ, WMHJ_AA, WMHJ_mu, WMHJ_tau, WMHJ_bbar, WMHJ_WW, WMHJ_ZZ_ll,  WMHJ_ZZ_lnu &
           )
         n_Vmassive = 1
         n_photons = 0
         higgs_process = .true.
      case default
         n_Vmassive = (N_lepton + N_decayquarks) / 2
         higgs_process = .false.
         n_photons = 0
         do i=2*n_Vmassive+1, n_v        ! assumption: leptons/decay quarks first in array v
            if (particle_IDs(i).eq.22) n_photons = n_photons + 1
         enddo
      end select


!----------------------------      
!   list of available scale choices
!----------------------------      

      SELECT CASE(scaleid)

! user defined constant scale 
      CASE(0)
         qsq(1:Lmax) = mu_user**2 * xi**2


! dynamical scale = ptj_min
      CASE(2)                 ! this scale makes only sense for the processes with a jet at LO,
                              ! but this is checked in InitScales and again here...
         atLeastOneJet = .false.
         do L=1, Lmax
            atLeastOneJet = atLeastOneJet .or. (njets(L).gt.0)
            qsq(l) = 1.0d12
            do i = 1, njets(L)                        ! njets(L) = 0 may occur for some kinematics, but for these kinematics the
                                                      ! amplitude is not evaluated. For at least 1 kinematic njets(L) > 0...
               ptsq = MAX(jets(5,i,L)**2,4D0)         ! lower safety bound
               qsq(l) = min(qsq(l),ptsq)
            enddo
            qsq(l) = qsq(l) * xi**2
         enddo
         if (.not.atLeastOneJet) then
            print*, "The scale choice ID = 2 is not valid for this process!"
            print*, "There are no jets defined in any kinematic!"
            stop
         endif


! dynamical scale = invariant mass of final state from boson decay products
! allow for ID 3 and 4 to be compatible with old input-files
      CASE(3,4)
         do mu = 0,3
            qvv(mu,1:lmax) = 0d0
            do i = 1,n_v
               qvv(mu,1:lmax) = qvv(mu,1:lmax) + v(mu,i,1:lmax)
            enddo
         enddo
         qvv(4,1:lmax) = qvv(0,1:lmax)**2 - qvv(1,1:lmax)**2 - qvv(2,1:lmax)**2 - qvv(3,1:lmax)**2
         qsq(1:lmax) = xi**2 * qvv(4,1:lmax)


! constant scale = higgs mass
      CASE(6)
         qsq(1:Lmax) = hmass**2 * xi**2


! dynamical scale = minimal ET of "primary" bosons
      CASE(7)
      do L = 1,Lmax
         if (higgs_process) then
            qv(0:3) = 0d0
            do i=3,n_v
               qv(0:3) = qv(0:3) + v(0:3,i,l)
            enddo
            etmin = transverseEnergy(qv)
         else
            etmin = 1d20
         endif
         do i=1,n_Vmassive
            qv(0:3) = v(0:3,2*i-1,l) + v(0:3,2*i,l)
            etmin = min(etmin, transverseEnergy(qv))
         enddo
         do i=1,n_photons
            etmin = min(etmin, transverseEnergy(v(0,2*n_Vmassive+i,l)))
         enddo
         qsq(l) = etmin**2 * xi**2
         qsq(l) = max(qsq(l), 5D0 * xi**2)      ! safety lower bound
      enddo


! default
      CASE DEFAULT
         print*,"in calcvxjscales undefined scale choice "
         stop
      END SELECT
      

      RETURN
      END


!***************************************************************************
      SUBROUTINE calcQCDscales(p, v, jets, njets, leptons, nl, photons,  &
                                              nphotons, Lmax, nlo)
!***************************************************************************
!     INPUT 
!     p       : 4-momenta of the partons
!     v       : 4-momenta of all other particles (leptons, photons, ...)
!     jets    : the array with the jet information as it is returned by the 
!               jetdefinition routine.
!     njets   : number of jets
!     leptons : the array with the charged lepton information.
!     nl      : number of charged leptons
!     photons : the array with the photon information.
!     nphotons: number of photons
!     Lmax    : number of different kinematics needed
!     nlo     : LO( = 0) or NLO(= 1)
!
!     OUTPUT
!     factorization scale and alfas(ren. scale) is transferred via 
!     common block in scales.inc
!***************************************************************************
      implicit none
#include "scales.inc"
#include "global.inc"
#include "coupl.inc"

!     input variables
      real*8 p(0:3,max_p,max_kin)
      real*8 v(0:3,max_v,max_kin)
      real*8 jets(0:7,max_jets,max_kin) 
      real*8 leptons(0:8,max_v,max_kin)     
      real*8 photons(0:7,max_v,max_kin)     
      integer Lmax, njets(max_kin), nl(max_kin), nphotons(max_kin), nlo
      
      integer L,i
      real*8 tempvec(0:3),y12
      real*8 transverseEnergy
      external transverseEnergy
      real*8 temp,alphas5,mass2,mjj2
      external alphas5,mass2,mjj2
      
!----------------------------      
!   factorization scale
!----------------------------      
      SELECT CASE(id_muf)

! user defined constant scale 
      CASE(0)
        mufsq(1,1) = muf_user**2 * xif**2
        do L=1,Lmax
          mufsq(1,L) = mufsq(1,1)
          mufsq(2,L) = mufsq(1,1)
        enddo

! constant scale = higgs mass
      CASE(6)
         mufsq(1,1) = hmass**2 * xif**2
         do L = 1,Lmax
            mufsq(1,L) = mufsq(1,1)
            mufsq(2,L) = mufsq(1,1)
         enddo
         
! HT/2 HT = sum(pT partons) + sum( ET V )
      CASE(8)
          do L=1,Lmax
            temp=0d0
            do i = 1,n_p ! sum over parton pt
              temp = temp + sqrt(p(1,i,L)*p(1,i,L) + p(2,i,L)*p(2,i,L))
            enddo
            do i = 1,nphotons(L) ! sum over photon pt
              temp = temp + photons(5,i,L)
            enddo
            do i = 1,n_v-nphotons(L),2 ! sum over lepton pairs
              tempvec(:)= v(:,i,L) + v(:,i+1,L)
              temp = temp + transverseEnergy(tempvec)
            enddo
            
            temp=temp/2d0*xif
            temp = temp**2
            mufsq(1,L) = temp
            mufsq(2,L) = temp
         enddo

! sum(pT_j*exp(|y_j-y12|)) + sum( ET V ) ( see arXiv 1311.6738 )
      CASE(9)
          do L=1,Lmax
            temp=0d0
            y12 = 0.5d0*(jets(6,1,L)+jets(6,2,L)) ! avg rapidity of two hardest jets
            do i = 1,njets(L) ! sum over jets
              temp = temp + jets(5,i,L)*exp(abs(jets(6,i,L)-y12))
            enddo
            do i = 1,nphotons(L) ! sum over photon pt
              temp = temp + photons(5,i,L)
            enddo
            do i = 1,n_v-nphotons(L),2 ! sum over lepton pairs
              tempvec(:)= v(:,i,L) + v(:,i+1,L)
              temp = temp + transverseEnergy(tempvec)
            enddo
            
            temp=temp/2d0*xif
            temp = temp**2
            mufsq(1,L) = temp
            mufsq(2,L) = temp
         enddo


! E_T(jj) + E_T(VV)
      CASE(10)
          do L=1,Lmax
            tempvec = jets(0:3,1,L) + jets(0:3,2,L)
            temp = transverseEnergy(tempvec) 
            
            tempvec=0d0
            do i=1,n_v 
              tempvec = tempvec + v(0:3,i,L)
            enddo
            temp = temp + transverseEnergy(tempvec)
            
            temp=temp/2d0*xif
            temp = temp**2
            mufsq(1,L) = temp
            mufsq(2,L) = temp
         enddo


      end select
      
      
!----------------------------      
!   renormalization scale
!----------------------------      
      SELECT CASE(id_mur)

! user defined constant scale 
      CASE(0)
        mursq(1,1) = mur_user**2 * xir**2
        if(nlo.eq.0) then
          als(1,1) = alphas5(mursq(1,1),0)
        else
          als(1,1) = alphas5(mursq(1,1),1)
        endif
        do L=2,Lmax
          mursq(1,L) = mursq(1,1)
          als(1,L) = als(1,1)
        enddo

! constant scale = higgs mass
      CASE(6)
        mursq(1,1) = hmass**2 * xir**2
        if(nlo.eq.0) then
          als(1,1) = alphas5(mursq(1,1),0)
        else
          als(1,1) = alphas5(mursq(1,1),1)
        endif
        do L=2,Lmax
          mursq(1,L) = mursq(1,1)
          als(1,L) = als(1,1)
        enddo
         
! HT/2 HT = sum(pT partons) + sum( ET V )
      CASE(8)
          do L=1,Lmax
           if(Njets(L).ge.2) then
            temp=0d0
            do i = 1,n_p ! sum over parton pt
              temp = temp + sqrt(p(1,i,L)*p(1,i,L) + p(2,i,L)*p(2,i,L))
            enddo
            do i = 1,nphotons(L) ! sum over photon pt
              temp = temp + photons(5,i,L)
            enddo
            do i = 1,n_v-nphotons(L),2 ! sum over lepton pairs
              tempvec(:)= v(:,i,L) + v(:,i+1,L)
              temp = temp + transverseEnergy(tempvec)
            enddo
            
            temp=temp/2d0*xir
            temp = temp**2
            mursq(1,L) = temp
            if(nlo.eq.0) then
              als(1,L) = alphas5(mursq(1,L),0)
            else
              als(1,L) = alphas5(mursq(1,L),1)
            endif
           else
            als(1,L) = 0d0
           endif
         enddo

! sum(pT_j*exp(|y_j-y12|)) + sum( ET V ) ( see arXiv 1311.6738 )
      CASE(9)
          do L=1,Lmax
           if(Njets(L).ge.2) then
            temp=0d0
            y12 = 0.5d0*(jets(6,1,L)+jets(6,2,L)) ! avg rapidity of two hardest jets
            do i = 1,njets(L) ! sum over jets
              temp = temp + jets(5,i,L)*exp(abs(jets(6,i,L)-y12))
            enddo
            do i = 1,nphotons(L) ! sum over photon pt
              temp = temp + photons(5,i,L)
            enddo
            do i = 1,n_v-nphotons(L),2 ! sum over lepton pairs
              tempvec(:)= v(:,i,L) + v(:,i+1,L)
              temp = temp + transverseEnergy(tempvec)
            enddo
            
            temp=temp/2d0*xir
            temp = temp**2
            mursq(1,L) = temp
            if(nlo.eq.0) then
              als(1,L) = alphas5(mursq(1,L),0)
            else
              als(1,L) = alphas5(mursq(1,L),1)
            endif
           else
            als(1,L) = 0d0
           endif
         enddo


! E_T(jj) + E_T(VV)
      CASE(10)
          do L=1,Lmax
           if(Njets(L).ge.2) then
            tempvec = jets(0:3,1,L) + jets(0:3,2,L)
            temp = transverseEnergy(tempvec) 
            
            tempvec=0d0
            do i=1,n_v 
              tempvec = tempvec + v(0:3,i,L)
            enddo
            temp = temp + transverseEnergy(tempvec)
            
            temp=temp/2d0*xir
            temp = temp**2
            mursq(1,L) = temp
            if(nlo.eq.0) then
              als(1,L) = alphas5(mursq(1,L),0)
            else
              als(1,L) = alphas5(mursq(1,L),1)
            endif
           else
            als(1,L) = 0d0
           endif
         enddo

      end select

      end


!*****************************************************************
      SUBROUTINE Calc_Momentum_Transfer(p, v, q12, q34, Lmax)
!*****************************************************************

      implicit none

#include "global.inc"

      real*8 p(0:3,max_p,max_kin), v(0:3,max_v,max_kin)
      real*8 q12(0:4,max_kin), q34(0:4,max_kin)
      integer Lmax, L, mu

      
      do mu = 0,3
         q12(mu,1) = p(mu,1,1)-p(mu,5,1)-p(mu,3,1)
         q34(mu,1) = p(mu,2,1)-p(mu,5,1)-p(mu,4,1)         
         do L = 2,Lmax
            q12(mu,L) = p(mu,1,L)-p(mu,3,L)
            q34(mu,L) = p(mu,2,L)-p(mu,4,L)
         enddo
      enddo
      do L = 1,Lmax
         q12(4,L) = abs(q12(0,L)**2-q12(1,L)**2-q12(2,L)**2-q12(3,L)**2)
         q34(4,L) = abs(q34(0,L)**2-q34(1,L)**2-q34(2,L)**2-q34(3,L)**2)
      enddo

      RETURN
      END


!*******************************************************************************
!*** calculate transverse energy of four vector p
!*** E_T = sqrt( p_T^2 + m^2 )
!*******************************************************************************
      double precision function transverseEnergy(p)

      implicit none
      double precision p(0:3)
      
      transverseEnergy = sqrt(p(0)*p(0)-p(3)*p(3))
      return
      end


!*******************************************************************************
!*******************************************************************************
!**                                                                          ***
!** Now we move to the alpha_s routines.                                     ***
!**                                                                          ***
!*******************************************************************************
!*******************************************************************************

!****************************************************************************
!     alphas routine, now a wrapper for LHAPDF alphas and hard wired alphas
!     the iord parameter is ignored, when LHAPDF alphas is used.
!****************************************************************************

      DOUBLE PRECISION FUNCTION alphas5(q2ren,iord1)
      IMPLICIT NONE     
      REAL*8 q2ren
      INTEGER iord1
      double precision alphasPDF, alphas5_hardwired, alphas_mstw
      external alphasPDF, alphas5_hardwired, alphas_mstw
      
      integer pdflib
      COMMON/PDFparameters/pdflib


      SELECT CASE(pdflib)
      CASE(0)
         alphas5 = alphas5_hardwired(q2ren, iord1)
#ifdef WITH_LHA
      CASE(1)
         alphas5 = alphasPDF(sqrt(q2ren))
#endif
      CASE(3)
         alphas5 = alphas_mstw(q2ren, iord1)
      CASE default
         alphas5 = alphas5_hardwired(q2ren, iord1)
      END SELECT

      RETURN                                                                  
      END       

!=========================================================================
!===========  alphas5.for ================================================
!=========================================================================

      DOUBLE PRECISION FUNCTION alphas5_hardwired(q2ren,iord1)
! same as alphas but for fixed 5 flavors   
      IMPLICIT NONE                                                           
      INTEGER iord1
      REAL*8 q2ren, xnf, lnsc
      REAL*8 dlambda5
      REAL*8 pi
      parameter (pi=3.141592653589793238d0,xnf=5d0)

      LOGICAL*4 ldebug
      parameter (ldebug=.false.)	        !output debug information

      if (iord1.eq.1) then
!         dlambda5=0.2020            !for CTEQ4
         dlambda5=0.2260            !for CTEQ6 and CT10
      else
!         dlambda5=0.1810            !for CTEQ4
         dlambda5=0.1652            !for CTEQ6 and CT10
      endif
      if (ldebug) then
         print*," lambda 5 = ", DLAMBDA5
      endif

      IF ( ldebug ) THEN
         WRITE(6,*) "qren",sqrt(q2ren)
         WRITE(6,*) "iord1=",iord1
      ENDIF
       
      lnsc = LOG(q2ren/dlambda5**2)
      IF (iord1 .EQ. 0 ) THEN
         alphas5_hardwired = 12.d0*pi/( (33.d0-2.d0*xnf)*lnsc )   
      ELSEIF (iord1 .EQ. 1) THEN
         alphas5_hardwired = 12.d0*pi/( (33.d0-2.d0*xnf)*lnsc ) &
               *( 1.d0 - 6.d0*(153.d0-19.d0*xnf)/(33.d0-2.d0*xnf)**2 &
                        *LOG(lnsc)/lnsc )
         ! this breaks if lnsc<0 
         ! which happens for really small scales
! fixme: write workaround/warning
      ELSE
         alphas5_hardwired = 0.1185
      ENDIF                                                      

!ASFIX
!asfix      alphas5 = 0.1185
!ASFIX
      if (ldebug) print*," alphas = ",alphas5_hardwired
      RETURN                                                                  
      END                                                                     
!


!============================================================================
!===========  alphas.for ====================================================
!============================================================================

      DOUBLE PRECISION FUNCTION alphas(q2ren,xnf,iord1)   
      IMPLICIT NONE                                                           
      INTEGER iord1
      REAL*8 q2ren,xnf
      REAL*8 dlambda3,dlambda5,dlambda6,dlambda1
      REAL*8 pi
      parameter (pi=3.141592653589793238d0)
      REAL*8 amc,amb,amt
      parameter (amc=1.35d0,amb=4.5d0,amt=175.0d0)       

      LOGICAL*4 ldebug
!---------------------------
      REAL*8 DLAMBDA4
      COMMON/QCDLAM/DLAMBDA4
!---------------------------
      parameter (ldebug=.true.)	        !output debug information
      if (iord1.eq.1) then
         dlambda4=0.325
      else
         dlambda4=0.215
      endif
      if (ldebug) then
         print*," lambda 4 = ", DLAMBDA4," nf = ",xnf
      endif

      IF (iord1 .eq. 0) THEN                                                  
         dlambda3=dlambda4*(amc/dlambda4)**(2./27.)                           
         dlambda5=dlambda4*(dlambda4/amb)**(2./23.)                           
         dlambda6=dlambda5*(dlambda5/amt)**(2./21.)                           
      ELSEIF (iord1 .eq. 1) THEN
         dlambda3=dlambda4*(amc/dlambda4)**(2./27.)                            &
                 *(dlog(amc**2/dlambda4**2))**(107./2025)                     
         dlambda5=dlambda4*(dlambda4/amb)**(2./23.)                            &
                 *(dlog(amb**2/dlambda4**2))**(-963./13225.)                  
         dlambda6=dlambda5*(dlambda5/amt)**(2./21.)                            &
                 *(dlog(amt**2/dlambda5**2))**(-321./3381.)                   
      ENDIF                                                                   

      IF (xnf .EQ. 4) THEN                                                   
         dlambda1=dlambda4                                                    
      ELSEIF (xnf .EQ. 5) THEN                                               
         dlambda1=dlambda5                                                    
      ELSEIF (xnf .EQ. 6) THEN   
         dlambda1=dlambda6                                   
      ELSE
         WRITE(6,*)"warning: xnf not defined in alphas"
         STOP
      ENDIF                                                                   

      IF ( ldebug ) THEN
         WRITE(6,*) "lambda_3=",dlambda3                                      
         WRITE(6,*) "lambda_4=",dlambda4                                      
         WRITE(6,*) "lambda_5=",dlambda5                                      
         WRITE(6,*) "lambda_6=",dlambda6                                      
         WRITE(6,*) "dlambda1",dlambda1
         WRITE(6,*) "q2ren",q2ren
         WRITE(6,*) "xnf",xnf
         WRITE(6,*) "iord1=",iord1
       ENDIF
       
       IF (iord1 .EQ. 0 ) THEN
          alphas = 12.d0*pi/( (33.d0-2.d0*xnf)*LOG(q2ren/dlambda1**2) )   
       ELSEIF (iord1 .EQ. 1) THEN
          alphas = 12.d0*pi/( (33.d0-2.d0*xnf)*LOG(q2ren/dlambda1**2) ) &
               *( 1.d0 - 6.d0*(153.d0-19.d0*xnf)/(33.d0-2.d0*xnf)**2 &
                        *LOG(LOG(q2ren/dlambda1**2)) &
                        /LOG(q2ren/dlambda1**2)               )
      ELSE
         alphas = 0.1185
      ENDIF                                                      

!ASFIX
!asfix      alphas = 0.1185
!ASFIX
      RETURN                                                                  
      END


!*******************************************************************************
!*******************************************************************************

!*******************************************************************************
!** The next section of this file contains the evolution of alpha_s for      ***
!** MSTW2008.                                                                ***
!** The first function is the VBFNLO interface to the alphas routines from   ***
!** the downloaded mstw package:                                             ***
!**      http://projects.hepforge.org/mstwpdf/code/code.html                 ***
!** Following the function alphas_mstw, all from the mstw alphas routines.   ***
!** Only alterations are names of functions/subroutines, to avoid clashes    ***
!** and a few purely cosmetic things.                                        ***
!*******************************************************************************

!*******************************************************************************
!*******************************************************************************

      double precision function alphas_mstw(qsq,iord1)

      implicit none

!* iord1 = lo / nlo; qsq = q^2
      double precision qsq
      integer iord1

!* fr2 is muf / mur, but I've set it to 1 here to agree with LHAPDF
      double precision fr2

      double precision mstwALPHAS
      external mstwALPHAS


!** ALPHAS routines initialised with PDFs
      fr2 = 1D0     
      if (iord1 .eq. 0) then
         call mstwINITALPHAS(iord1,fr2,1D0,0.68183D0,1.40D0,4.75D0,1.D10)
      else
         call mstwINITALPHAS(iord1,fr2,1D0,0.49128D0,1.40D0,4.75D0,1.D10)
      end if

!* A CHECK:  
!      write(*,*)'order, alphasMZ =', iord1, mstwALPHAS(91.1876D0)
!     write(*,*)'LO alphas(MZ) should be: 0.13939'
!     write(*,*)'NLO alphas(MZ) should be: 0.12018'

      alphas_mstw = mstwALPHAS(sqrt(qsq))


      end 


!*******************************************************************************
!*******************************************************************************
!----------------------------------------------------------------------
!--   Stand-alone code for alpha_s cannibalised (with permission)
!--   from Andreas Vogt's QCD-PEGASUS package (hep-ph/0408244).
!--   The running coupling alpha_s is obtained at N^mLO (m = 0,1,2,3)
!--   by solving the renormalisation group equation in the MSbar scheme
!--   by a fourth-order Runge-Kutta integration.  Transitions from
!--   n_f to n_f+1 flavours are made when the factorisation scale
!--   mu_f equals the pole masses m_h (h = c,b,t).  At exactly
!--   the thresholds m_{c,b,t}, the number of flavours n_f = {3,4,5}.
!--   The top quark mass should be set to be very large to evolve with
!--   a maximum of five flavours.  The factorisation scale mu_f may be
!--   a constant multiple of the renormalisation scale mu_r.  The input
!--   factorisation scale mu_(f,0) should be less than or equal to
!--   the charm quark mass.  However, if it is greater than the
!--   charm quark mass, the value of alpha_s at mu_(f,0) = 1 GeV will
!--   be found using a root-finding algorithm.
!--
!--   Example of usage.
!--   First call the initialisation routine (only needed once):
!--
!--    IORD = 2                  ! perturbative order (N^mLO,m=0,1,2,3)
!--    FR2 = 1.D0                ! ratio of mu_f^2 to mu_r^2
!--    MUR = 1.D0                ! input mu_r in GeV
!--    ASMUR = 0.5D0             ! input value of alpha_s at mu_r
!--    MC = 1.4D0                ! charm quark mass
!--    MB = 4.75D0               ! bottom quark mass
!--    MT = 1.D10                ! top quark mass
!--    CALL mstwINITALPHAS(IORD, FR2, MUR, ASMUR, MC, MB, MT)
!--
!--   Then get alpha_s at a renormalisation scale mu_r with:
!--
!--    MUR = 100.D0              ! renormalisation scale in GeV
!--    ALFAS = mstwALPHAS(MUR)
!--
!----------------------------------------------------------------------
!--   Comments to Graeme Watt <Graeme.Watt(at)cern.ch>
!----------------------------------------------------------------------

      subroutine mstwINITALPHAS(IORD, FR2, MUR, ASMUR, MC, MB, MT)

!--   IORD = 0 (LO), 1 (NLO), 2 (NNLO), 3 (NNNLO).
!--   FR2 = ratio of mu_f^2 to mu_r^2 (must be a fixed value).
!--   MUR = input renormalisation scale (in GeV) for alpha_s.
!--   ASMUR = input value of alpha_s at the renormalisation scale MUR.
!--   MC,MB,MT = heavy quark masses in GeV.

      IMPLICIT NONE

      INTEGER IORD,IORDc,MAXF,MODE
      DOUBLE PRECISION FR2,MUR,ASMUR,MC,MB,MT,EPS,A,B,MSTW_DZEROX, &
           R0c,FR2c,MURc,ASMURc,MCc,MBc,MTc,mstwFINDALPHASR0,R0,ASI
      COMMON / MSTW_DZEROXcommon / FR2c,MURc,ASMURc,MCc,MBc,MTc,R0c,IORDc
      PARAMETER(EPS=1.D-10,MAXF=10000,MODE=1)
      EXTERNAL mstwFINDALPHASR0

      IF (MUR*sqrt(FR2).LE.MC) THEN ! Check that MUF <= MC.
         R0 = MUR
         ASI = ASMUR
      ELSE                      ! Solve for alpha_s at R0 = 1 GeV.
!--   Copy variables to common block.
         R0c = 1.D0/sqrt(FR2)
         IORDc = IORD
         FR2c = FR2
         MURc = MUR
         ASMURc = ASMUR
         MCc = MC
         MBc = MB
         MTc = MT
!--   Now get alpha_s(R0) corresponding to alpha_s(MUR).
         A = 0.02D0              ! lower bound for alpha_s(R0)
         B = 2.00D0              ! upper bound for alpha_s(R0)
         R0 = R0c
         ASI = MSTW_DZEROX(A,B,EPS,MAXF,mstwFINDALPHASR0,MODE)
      END IF

      CALL mstwINITALPHASR0(IORD, FR2, R0, ASI, MC, MB, MT)

      RETURN
      END


!*******************************************************************************
!*******************************************************************************

!--   Find the zero of this function using MSTW_DZEROX.

      DOUBLE PRECISION FUNCTION mstwFINDALPHASR0(ASI)

      IMPLICIT NONE

      INTEGER IORD
      DOUBLE PRECISION FR2, R0, ASI, MC, MB, MT, MUR, ASMUR, mstwALPHAS
      COMMON / MSTW_DZEROXcommon / FR2, MUR, ASMUR, MC, MB, MT, R0, IORD

      CALL mstwINITALPHASR0(IORD, FR2, R0, ASI, MC, MB, MT)
      mstwFINDALPHASR0 = mstwALPHAS(MUR) - ASMUR ! solve equal to zero

      RETURN
      END


!*******************************************************************************
!*******************************************************************************

      SUBROUTINE mstwINITALPHASR0(IORD, FR2, R0, ASI, MC, MB, MT)

!--   IORD = 0 (LO), 1 (NLO), 2 (NNLO), 3 (NNNLO).
!--   FR2 = ratio of mu_f^2 to mu_r^2 (must be a fixed value).
!--   R0 = input renormalisation scale (in GeV) for alphas_s.
!--   ASI = input value of alpha_s at the renormalisation scale R0.
!--   MC,MB,MT = heavy quark masses in GeV.
!--   Must have R0*sqrt(FR2) <= MC to call this subroutine.

      IMPLICIT NONE

      INTEGER IORD,NAORD,NASTPS,IVFNS,NFF
      DOUBLE PRECISION FR2,R0,ASI,MC,MB,MT,LOGFR,R20, &
           PI,ZETA,CF,CA,TR,AS0,M20,MC2,MB2,MT2
      PARAMETER(PI = 3.14159265358979D0)

      COMMON / RZETA  / ZETA(6)
      COMMON / COLOUR / CF, CA, TR
      COMMON / ASINP  / AS0, M20
      COMMON / ASPAR  / NAORD, NASTPS
      COMMON / VARFLV / IVFNS
      COMMON / NFFIX  / NFF
      COMMON / FRRAT  / LOGFR

!
! ..QCD colour factors
!
      CA = 3.D0
      CF = 4./3.D0
      TR = 0.5D0
!
! ..The lowest integer values of the Zeta function
!
      ZETA(1) = 0.57721566490153D0
      ZETA(2) = 1.644934066848226D0
      ZETA(3) = 1.202056903159594D0
      ZETA(4) = 1.082323233711138D0
      ZETA(5) = 1.036927755143370D0
      ZETA(6) = 1.017343061984449D0

      IVFNS = 1                 ! variable flavour-number scheme (VFNS)
!      IVFNS = 0                 ! fixed flavour-number scheme (FFNS)
      NFF = 4                   ! number of flavours for FFNS
      NAORD = IORD              ! perturbative order of alpha_s
      NASTPS = 20               ! num. steps in Runge-Kutta integration
      R20 = R0**2               ! input renormalisation scale
      MC2 = MC**2               ! mu_f^2 for charm threshold
      MB2 = MB**2               ! mu_f^2 for bottom threshold
      MT2 = MT**2               ! mu_f^2 for top threshold
      LOGFR = LOG(FR2)          ! log of ratio of mu_f^2 to mu_r^2
      M20 = R20 * FR2           ! input factorisation scale

!
! ..Stop some nonsense
!
      IF ( (IVFNS .EQ. 0) .AND. (NFF .LT. 3) ) THEN
         WRITE (6,*) 'Wrong flavour number for FFNS evolution. STOP'
         STOP
      END IF
      IF ( (IVFNS .EQ. 0) .AND. (NFF .GT. 5) ) THEN
         WRITE (6,*) 'Wrong flavour number for FFNS evolution. STOP'
         STOP
      END IF
!     
      IF ( NAORD .GT. 3 ) THEN
         WRITE (6,*) 'Specified order in a_s too high. STOP' 
         STOP
      END IF
!
      IF ( (IVFNS .NE. 0) .AND. (FR2 .GT. 4.001D0) ) THEN
         WRITE (6,*) 'Too low mu_r for VFNS evolution. STOP'
         STOP
      END IF
!
      IF ( (IVFNS .EQ. 1) .AND. (M20 .GT. MC2) ) THEN
         WRITE (6,*) 'Too high mu_0 for VFNS evolution. STOP'
         STOP
      END IF
!     
      IF ( (ASI .GT. 2.D0) .OR. (ASI .LT. 2.D-2) ) THEN
         WRITE (6,*) 'alpha_s out of range. STOP'
         STOP
      END IF
!     
      IF ( (IVFNS .EQ. 1) .AND. (MC2 .GT. MB2) ) THEN
         WRITE (6,*) 'Wrong charm-bottom mass hierarchy. STOP'
         STOP
      END IF
      IF ( (IVFNS .EQ. 1) .AND. (MB2 .GT. MT2) ) THEN
         WRITE (6,*) 'Wrong bottom-top mass hierarchy. STOP'
         STOP
      END IF
!

!--   Store the beta function coefficients in a COMMON block.
      CALL MSTWBETAFCT

!--   Store a_s = alpha_s(mu_r^2)/(4 pi) at the input scale R0.
      AS0 = ASI / (4.D0* PI)

!--   Store alpha_s at the heavy flavour thresholds in a COMMON block.
       IF (IVFNS .NE. 0) THEN
          CALL MSTWEVNFTHR (MC2, MB2, MT2)
       END IF

      RETURN
      END


!*******************************************************************************
!*******************************************************************************

      DOUBLE PRECISION FUNCTION mstwALPHAS(MUR)

      IMPLICIT NONE

      INTEGER NFF,IVFNS,NF
      DOUBLE PRECISION PI,LOGFR,AS0,M20,ASC,M2C,ASB,M2B,AST,M2T,M2,MUR, &
           R2,ASI,ASF,R20,R2T,R2B,R2C,MSTW_AS
      PARAMETER ( PI = 3.14159265358979D0 )
!
! ..Input common blocks 
! 
       COMMON / NFFIX  / NFF
       COMMON / VARFLV / IVFNS 
       COMMON / FRRAT  / LOGFR
       COMMON / ASINP  / AS0, M20
       COMMON / ASFTHR / ASC, M2C, ASB, M2B, AST, M2T

       R2 = MUR**2
       M2 = R2 * EXP(+LOGFR)
       IF (IVFNS .EQ. 0) THEN
!
!   Fixed number of flavours
!
          NF  = NFF
          R20 = M20 * R2/M2
          ASI = AS0
          ASF = MSTW_AS (R2, R20, AS0, NF)
!
       ELSE
!
! ..Variable number of flavours
!
          IF (M2 .GT. M2T) THEN
             NF = 6
             R2T = M2T * R2/M2
             ASI = AST
             ASF = MSTW_AS (R2, R2T, AST, NF)
!
          ELSE IF (M2 .GT. M2B) THEN
             NF = 5
             R2B = M2B * R2/M2
             ASI = ASB
             ASF = MSTW_AS (R2, R2B, ASB, NF)
!     
          ELSE IF (M2 .GT. M2C) THEN
             NF = 4
             R2C = M2C * R2/M2
             ASI = ASC
             ASF = MSTW_AS (R2, R2C, ASC, NF)
!     
          ELSE
             NF = 3
             R20 = M20 * R2/M2
             ASI = AS0
             ASF = MSTW_AS (R2, R20, AS0, NF)
!       
          END IF
!
       END IF
!
! ..Final value of alpha_s
!
       mstwALPHAS = 4.D0*PI*ASF

       RETURN
       END


!*******************************************************************************
!*******************************************************************************
! =====================================================================
!
! ..The threshold matching of the QCD coupling in the MS(bar) scheme,  
!    a_s = alpha_s(mu_r^2)/(4 pi),  for  NF -> NF + 1  active flavours 
!    up to order a_s^4 (NNNLO).
!
! ..The value  ASNF  of a_s for NF flavours at the matching scale, the 
!    logarithm  LOGRH = ln (mu_r^2/m_H^2) -- where m_H is the pole mass
!    of the heavy quark -- and  NF  are passed as arguments to the 
!    function  mstw_ASNF1.  The order of the expansion  NAORD  (defined as 
!    the 'n' in N^nLO) is provided by the common-block  ASPAR.
!
! ..The matching coefficients are inverted from Chetyrkin, Kniehl and
!    Steinhauser, Phys. Rev. Lett. 79 (1997) 2184. The QCD colour
!    factors have been hard-wired in these results. The lowest integer 
!    values of the Zeta function are given by the common-block  RZETA.
!
! =====================================================================
!
      DOUBLE PRECISION FUNCTION mstw_ASNF1 (ASNF, LOGRH, NF)
!
      IMPLICIT NONE
      INTEGER NF, NAORD, NASTPS, PRVCLL, K1, K2
      DOUBLE PRECISION ASNF,LOGRH,ZETA,CMC,CMCI30,CMCF30,CMCF31, &
           CMCI31,ASP,LRHP

      DIMENSION CMC(3,0:3)
!
! ---------------------------------------------------------------------
!
! ..Input common-blocks 
!
      COMMON / ASPAR  / NAORD, NASTPS
      COMMON / RZETA  / ZETA(6)
!
! ..Variables to be saved for the next call
!
      SAVE CMC, CMCI30, CMCF30, CMCF31, CMCI31, PRVCLL
!
! ---------------------------------------------------------------------
!
! ..The coupling-constant matching coefficients (CMC's) up to NNNLO 
!   (calculated and saved in the first call of this routine)
!
       IF (PRVCLL .NE. 1) THEN
!
         CMC(1,0) =  0.D0
         CMC(1,1) =  2./3.D0
!
         CMC(2,0) = 14./3.D0
         CMC(2,1) = 38./3.D0
         CMC(2,2) =  4./9.D0  
!
         CMCI30 = + 80507./432.D0 * ZETA(3) + 58933./1944.D0  &
                  + 128./3.D0 * ZETA(2) * (1.+ DLOG(2.D0)/3.D0)
         CMCF30 = - 64./9.D0 * (ZETA(2) + 2479./3456.D0)
         CMCI31 =   8941./27.D0
         CMCF31 = - 409./27.D0
         CMC(3,2) = 511./9.D0
         CMC(3,3) = 8./27.D0
!
         PRVCLL = 1
!
       END IF
!
! ---------------------------------------------------------------------
!
! ..The N_f dependent CMC's, and the alpha_s matching at order NAORD 
!
       CMC(3,0) = CMCI30 + NF * CMCF30
       CMC(3,1) = CMCI31 + NF * CMCF31
!
       mstw_ASNF1 = ASNF
       IF (NAORD .EQ. 0) GO TO 1
       ASP   = ASNF
!
       DO 11 K1 = 1, NAORD 
         ASP = ASP * ASNF
         LRHP = 1.D0
!
       DO 12 K2 = 0, K1
         mstw_ASNF1 = mstw_ASNF1 + ASP * CMC(K1,K2) * LRHP
         LRHP = LRHP * LOGRH
!
  12   CONTINUE
  11   CONTINUE
!
! ---------------------------------------------------------------------
!
   1   RETURN
       END


!*******************************************************************************
!*******************************************************************************
! =================================================================av==
!
! ..The subroutine  MSTWEVNFTHR  for the evolution of  a_s = alpha_s/(4 pi)
!    from a three-flavour initial scale to the four- to six-flavour
!    thresholds (identified with the squares of the corresponding quark
!    masses).  The results are written to the common-block  ASFTHR.
!
! ..The input scale  M20 = mu_(f,0)^2  and the corresponding value 
!    AS0  of a_s  are provided by  ASINP.  The fixed scale logarithm
!    LOGFR = ln (mu_f^2/mu_r^2) is specified in  FRRAT.  The alpha_s
!    matching is done by the function mstw_ASNF1.
!
! =====================================================================
!
       SUBROUTINE MSTWEVNFTHR (MC2, MB2, MT2)
!
       IMPLICIT NONE
       DOUBLE PRECISION MC2, MB2, MT2, M20, M2C, M2B, M2T, R20, R2C,  &
                        R2B, R2T, MSTW_AS, mstw_ASNF1, AS0, ASC, ASB, AST, &
                        ASC3, ASB4, AST5, LOGFR, SC, SB, ST
!
! ---------------------------------------------------------------------
! 
! ..Input common blocks
!  
       COMMON / ASINP  / AS0, M20
       COMMON / FRRAT  / LOGFR
!
! ..Output common blocks
!
       COMMON / ASFTHR / ASC, M2C, ASB, M2B, AST, M2T

! ---------------------------------------------------------------------
!
! ..Coupling constants at and evolution distances to/between thresholds
! 
       R20 = M20 * EXP(-LOGFR)
!
! ..Charm
!
       M2C  = MC2
       R2C  = M2C * R20/M20
       ASC3 = MSTW_AS (R2C, R20, AS0, 3)
       SC   = LOG (AS0 / ASC3)
       ASC  = mstw_ASNF1 (ASC3, -LOGFR, 3)
!
! ..Bottom 
!
       M2B  = MB2
       R2B  = M2B * R20/M20
       ASB4 = MSTW_AS (R2B, R2C, ASC, 4)
       SB   = LOG (ASC / ASB4)
       ASB  = mstw_ASNF1 (ASB4, -LOGFR, 4)
!
! ..Top
!
       M2T  = MT2
       R2T  = M2T * R20/M20
       AST5 = MSTW_AS (R2T, R2B, ASB, 5)
       ST   = LOG (ASB / AST5)
       AST  = mstw_ASNF1 (AST5, -LOGFR, 5)

       RETURN
       END


!*******************************************************************************
!*******************************************************************************
! =================================================================av==
!
! ..The running coupling of QCD,  
!
!         MSTW_AS  =  a_s  =  alpha_s(mu_r^2)/(4 pi),
!
!    obtained by integrating the evolution equation for a fixed number
!    of massless flavours  NF.  Except at leading order (LO),  MSTW_AS  is 
!    obtained using a fourth-order Runge-Kutta integration.
!
! ..The initial and final scales  R20  and  R2,  the value  AS0  at
!    R20, and  NF  are passed as function arguments.  The coefficients 
!    of the beta function up to  a_s^5 (N^3LO)  are provided by the 
!    common-block  BETACOM.  The order of the expansion  NAORD (defined
!    as the 'n' in N^nLO) and the number of steps  NASTPS  for the 
!    integration beyond LO are given by the common-block  ASPAR.
!
! =====================================================================
!
      DOUBLE PRECISION FUNCTION MSTW_AS (R2, R20, AS0, NF)
!
      IMPLICIT NONE

      INTEGER NFMIN, NFMAX, NF, NAORD, NASTPS, K1
      DOUBLE PRECISION R2, R20, AS0, SXTH, BETA0, BETA1, BETA2, BETA3, &
           FBETA1,FBETA2,FBETA3,A,LRRAT,DLR,XK0,XK1,XK2,XK3
      PARAMETER (NFMIN = 3, NFMAX = 6)
!
! ---------------------------------------------------------------------
!
! ..Input common-blocks 
!
       COMMON / ASPAR  / NAORD, NASTPS
       COMMON / BETACOM   / BETA0 (NFMIN:NFMAX), BETA1 (NFMIN:NFMAX), &
                         BETA2 (NFMIN:NFMAX), BETA3 (NFMIN:NFMAX)
!
! ..The beta functions FBETAn at N^nLO for n = 1, 2, and 3
!
       FBETA1(A) = - A**2 * ( BETA0(NF) + A *   BETA1(NF) )
       FBETA2(A) = - A**2 * ( BETA0(NF) + A * ( BETA1(NF) &
                              + A * BETA2(NF) ) )
       FBETA3(A) = - A**2 * ( BETA0(NF) + A * ( BETA1(NF) &
                              + A * (BETA2(NF) + A * BETA3(NF)) ) )
!
! ---------------------------------------------------------------------
!
! ..Initial value, evolution distance and step size
!
       MSTW_AS = AS0
       LRRAT = LOG (R2/R20)
       DLR = LRRAT / NASTPS
!
! ..Solution of the evolution equation depending on  NAORD
!   (fourth-order Runge-Kutta beyond the leading order)
!
       IF (NAORD .EQ. 0) THEN
!
         MSTW_AS = AS0 / (1.+ BETA0(NF) * AS0 * LRRAT)
!
       ELSE IF (NAORD .EQ. 1) THEN
!
       DO 2 K1 = 1, NASTPS
         XK0 = DLR * FBETA1 (MSTW_AS)
         XK1 = DLR * FBETA1 (MSTW_AS + 0.5 * XK0)
         XK2 = DLR * FBETA1 (MSTW_AS + 0.5 * XK1)
         XK3 = DLR * FBETA1 (MSTW_AS + XK2)
         MSTW_AS = MSTW_AS + 1d0/6d0 * (XK0 + 2.* XK1 + 2.* XK2 + XK3)
  2    CONTINUE
!
       ELSE IF (NAORD .EQ. 2) THEN
!
       DO 3 K1 = 1, NASTPS
         XK0 = DLR * FBETA2 (MSTW_AS)
         XK1 = DLR * FBETA2 (MSTW_AS + 0.5 * XK0)
         XK2 = DLR * FBETA2 (MSTW_AS + 0.5 * XK1)
         XK3 = DLR * FBETA2 (MSTW_AS + XK2)
         MSTW_AS = MSTW_AS + SXTH * (XK0 + 2.* XK1 + 2.* XK2 + XK3)
  3    CONTINUE
!  
       ELSE IF (NAORD .EQ. 3) THEN
!
       DO 4 K1 = 1, NASTPS
         XK0 = DLR * FBETA3 (MSTW_AS)
         XK1 = DLR * FBETA3 (MSTW_AS + 0.5 * XK0)
         XK2 = DLR * FBETA3 (MSTW_AS + 0.5 * XK1)
         XK3 = DLR * FBETA3 (MSTW_AS + XK2)
         MSTW_AS = MSTW_AS + SXTH * (XK0 + 2.* XK1 + 2.* XK2 + XK3)
  4    CONTINUE
       END IF
!
! ---------------------------------------------------------------------
!
       RETURN
       END


!*******************************************************************************
!*******************************************************************************
! =================================================================av==
!
! ..The subroutine MSTWBETAFCT for the coefficients  BETA0...BETA3  of the 
!    beta function of QCD up to order alpha_s^5 (N^3LO), normalized by 
!
!        d a_s / d ln mu_r^2  =  - BETA0 a_s^2 - BETA1 a_s^3 - ... 
!
!    with  a_s = alpha_s/(4*pi). 
!
! ..The MSbar coefficients are written to the common-block BETACOM for 
!   NF = 3...6  (parameters NFMIN, NFMAX) quark flavours.
!
! ..The factors CF, CA and TF  are taken from the common-block  COLOUR.
!    Beyond NLO the QCD colour factors are hard-wired in this routine,
!    and the numerical coefficients are truncated to six digits.
!
! =====================================================================
!
       SUBROUTINE MSTWBETAFCT

       implicit none

       double precision B00, B01, B10, B11
       double precision CF, CA, TR
       double precision BETA0, BETA1, BETA2, BETA3
!
       INTEGER NFMIN, NFMAX, NF
       PARAMETER (NFMIN = 3, NFMAX = 6)
!
! ---------------------------------------------------------------------
!
! ..Input common-block
!
       COMMON / COLOUR / CF, CA, TR
!
! ..Output common-block
!
       COMMON / BETACOM   / BETA0 (NFMIN:NFMAX), BETA1 (NFMIN:NFMAX), &
                         BETA2 (NFMIN:NFMAX), BETA3 (NFMIN:NFMAX)
!
! ..The full LO and NLO coefficients 
!
       B00 =  11./3.D0 * CA
       B01 =  -4./3.D0 * TR
       B10 =  34./3.D0 * CA**2
       B11 = -20./3.D0 * CA*TR - 4.* CF*TR
!
! ..Flavour-number loop and output to the array
!
       DO 1 NF = NFMIN, NFMAX
!
       BETA0(NF) = B00 + B01 * NF
       BETA1(NF) = B10 + B11 * NF
!
       BETA2(NF) = 1428.50 - 279.611 * NF + 6.01852 * NF**2
       BETA3(NF) = 29243.0 - 6946.30 * NF + 405.089 * NF**2  &
                   + 1.49931 * NF**3
!
  1    CONTINUE
!
       RETURN
       END


!*******************************************************************************
!*******************************************************************************

!--   G.W. MSTW_DZEROX taken from CERNLIB to find the zero of a function.
      DOUBLE PRECISION FUNCTION MSTW_DZEROX(A0,B0,EPS,MAXF,F,MODE)

      implicit none

      double precision Z1, A0, B0, EPS, F, FA, FB, FC, FD, FDA, FDB
      double precision HALF, A, B, C, D, ATL, H, P, Q, TOL, W, HB
      integer maxf, mode, ie, im1, im2, MF
!     Based on
!
!        J.C.P. Bus and T.J. Dekker, Two Efficient Algorithms with
!        Guaranteed Convergence for Finding a Zero of a Function,
!        ACM Trans. Math. Software 1 (1975) 330-345.
!
!        (MODE = 1: Algorithm M;    MODE = 2: Algorithm R)
      CHARACTER*80 ERRTXT
      LOGICAL LMT
      DIMENSION IM1(2),IM2(2),LMT(2)
      PARAMETER (Z1 = 1, HALF = Z1/2)
      DATA IM1 /2,3/, IM2 /-1,3/
      MSTW_DZEROX = 0.D0             ! G.W. to prevent compiler warning
      IF(MODE .NE. 1 .AND. MODE .NE. 2) THEN
       C=0
       WRITE(ERRTXT,101) MODE
       WRITE(6,*) ERRTXT
       GO TO 99
      ENDIF
      FA=F(B0)
      FB=F(A0)
      IF(FA*FB .GT. 0) THEN
       C=0
       WRITE(ERRTXT,102) A0,B0
       WRITE(6,*) ERRTXT
       GO TO 99
      ENDIF
      ATL=ABS(EPS)
      B=A0
      A=B0
      LMT(2)=.TRUE.
      MF=2
    1 C=A
      FC=FA
    2 IE=0
    3 IF(ABS(FC) .LT. ABS(FB)) THEN
       IF(C .NE. A) THEN
        D=A
        FD=FA
       END IF
       A=B
       B=C
       C=A
       FA=FB
       FB=FC
       FC=FA
      END IF
      TOL=ATL*(1+ABS(C))
      H=HALF*(C+B)
      HB=H-B
      IF(ABS(HB) .GT. TOL) THEN
       IF(IE .GT. IM1(MODE)) THEN
        W=HB
       ELSE
        TOL=TOL*SIGN(Z1,HB)
        P=(B-A)*FB
        LMT(1)=IE .LE. 1
        IF(LMT(MODE)) THEN
         Q=FA-FB
         LMT(2)=.FALSE.
        ELSE
         FDB=(FD-FB)/(D-B)
         FDA=(FD-FA)/(D-A)
         P=FDA*P
         Q=FDB*FA-FDA*FB
        END IF
        IF(P .LT. 0) THEN
         P=-P
         Q=-Q
        END IF
        IF(IE .EQ. IM2(MODE)) P=P+P
        IF(P .EQ. 0 .OR. P .LE. Q*TOL) THEN
         W=TOL
        ELSEIF(P .LT. HB*Q) THEN
         W=P/Q
        ELSE
         W=HB
        END IF
       END IF
       D=A
       A=B
       FD=FA
       FA=FB
       B=B+W
       MF=MF+1
       IF(MF .GT. MAXF) THEN
        WRITE(6,*) "Error in MSTW_DZEROX: TOO MANY FUNCTION CALLS"
        GO TO 99
       ENDIF
       FB=F(B)
       IF(FB .EQ. 0 .OR. SIGN(Z1,FC) .EQ. SIGN(Z1,FB)) GO TO 1
       IF(W .EQ. HB) GO TO 2
       IE=IE+1
       GO TO 3
      END IF
      MSTW_DZEROX=C
   99 CONTINUE
      RETURN
  101 FORMAT('Error in MSTW_DZEROX: MODE = ',I3,' ILLEGAL')
  102 FORMAT('Error in MSTW_DZEROX: F(A)&F(B) HAVE THE SAME SIGN, A = ', &
           1P,D15.8,', B = ',D15.8)
      END


!*******************************************************************************
