!******************************************************************************************
!
! Implementation of CPS Higgs scheme 
!
! Code by Giampiero Passarino
!   (see original header below)
!
! Adaptation for VBFNLO:
!   Michael Rauch <michael.rauch@kit.edu>
!   Last changed: January 2013
!
!******************************************************************************************

!-----------------------------------------------------------------------------------------*
!                                                                                         *
! cpHTO v 1.1 (May 2012) by Giampiero                                                     * 
!                                                                                         *
! History:                                                                                * 
! Version   Date            Comment                                                       *
! -------   ----            -------                                                       *
! 1.0       December 2011   Original code. Giampiero                                      *
!           March    2012   bug fixed in ghb                                              *
! 1.1       May      2012   Expanded solution for low muh                                 *
!                                                                                         *
! Given the parametrization s_H = muh^2 - i*muh*gh for the Higgs-boson complex pole       *
! computes gh at a given value of muh in the Standard Model                               *
! It also returns the OS width and mass and width in the bar-scheme                       *
!                                                                                         * 
!-----------------------------------------------------------------------*                 *
!     Comments to Giampiero Passarino <giampiero(at)to.infn.it>         *                 *
!-----------------------------------------------------------------------*                 *
!                                                                                         *
! %\cite{Goria:2011wa}                                                                    *   
! \bibitem{Goria:2011wa}                                                                  *
!   S.~Goria, G.~Passarino and D.~Rosco,                                                  *
!   %``The Higgs Boson Lineshape,''                                                       *
!   arXiv:1112.5517 [hep-ph].                                                             *
!   %%CITATION = ARXIV:1112.5517;%%                                                       *
!                                                                                         *
! %\cite{Passarino:2010qk}                                                                *
! \bibitem{Passarino:2010qk}                                                              *
!   G.~Passarino, C.~Sturm and S.~Uccirati,                                               *
!   %``Higgs Pseudo-Observables, Second Riemann Sheet and All That,''                     *
!   Nucl.\ Phys.\ B {\bf 834} (2010) 77                                                   *
!   [arXiv:1001.3360 [hep-ph]].                                                           *
!   %%CITATION = ARXIV:1001.3360;%%                                                       *
!                                                                                         *
!-----------------------------------------------------------------------------------------*

      MODULE VBFNLO_HIGGSTO
          implicit none
       INTEGER qcdc,pcnt,gtop
       REAL*8 cxe,cxmu,cxtau,cxu,cxd,cxc,cxs,cxt,cxb, &
     &        cxes,cxmus,cxtaus,cxus,cxds,cxcs,cxss,cxts,cxbs, &
     &        clxe,clxmu,clxtau,clxu,clxd,clxc,clxs,clxt,clxb,muhcp, &
     &        scalec,cxqs,cxq,cxls,cxl,cxtmb,cxtmbs,cxbi,cxtmbi,xq, &
     &        cxtc,cxbc,yimt
       REAL*8, dimension(2) :: cxwi,clcts,clwtb,cpw,cpz
      !END MODULE VBFNLO_HIGGSTO_aux_Hcp
!
!-----------------------------------------------------------------------------------------
!
      !MODULE VBFNLO_HIGGSTO_masses
       REAL*8 mt
       REAL*8, parameter :: m0= 0.d0
!mr
!       REAL*8, parameter :: mw= 80.398d0
!       REAL*8, parameter :: mz= 91.1876d0
       REAL*8, parameter :: me= 0.51099907d-3
       REAL*8, parameter :: mm= 0.105658389d0
       REAL*8, parameter :: mtl= 1.77684d0
       REAL*8, parameter :: muq= 0.190d0
       REAL*8, parameter :: mdq= 0.190d0
       REAL*8, parameter :: mcq= 1.55d0
       REAL*8, parameter :: msq= 0.190d0
!       REAL*8, parameter :: mbq= 4.69d0
!       REAL*8, parameter :: mb= 4.69d0
       REAL*8, parameter :: mtiny= 1.d-10
!       REAL*8, parameter :: imw= 2.08872d0
!       REAL*8, parameter :: imz= 2.4952d0
!       REAL*8, parameter :: swr= mw*mw-imw*imw
!       REAL*8, parameter :: swi= -mw*imw*(1.d0-0.5d0*(imw*imw)/(mw*mw))
!       REAL*8, parameter :: szr= mz*mz-imz*imz
!       REAL*8, parameter :: szi= -mz*imz*(1.d0-0.5d0*(imz*imz)/(mz*mz))
       REAL*8 mw, mz, mbq, mb, imw, imz, swr, swi, szr, szi
      !END MODULE VBFNLO_HIGGSTO_masses
!
!-----------------------------------------------------------------------------------------
!
      !MODULE VBFNLO_HIGGSTO_units
       INTEGER, parameter :: izer= 0
       INTEGER, parameter :: ione= 1
       INTEGER, parameter :: itwo= 2
       REAL*8, parameter :: eps= -1.d0
       REAL*8, parameter :: one= 1.d0
       REAL*8, parameter :: zero= 0.d0
       REAL*8, parameter :: qeps= 1.d-25
       REAL*8, parameter :: ez= 1.d-15
       REAL*8 :: co(1:2)=(/1.d0,0.d0/)
       REAL*8 :: ci(1:2)=(/0.d0,1.d0/)
       REAL*8 :: cz(1:2)=(/0.d0,0.d0/)
      !END MODULE VBFNLO_HIGGSTO_units
!
!-----------------------------------------------------------------------------------------
!
      !MODULE VBFNLO_HIGGSTO_riemann
       REAL*8, parameter :: pi= 3.141592653589793238462643d0
       REAL*8, parameter :: pis= pi*pi
       REAL*8, parameter :: piq= pis*pis
       REAL*8, parameter :: rz2= 1.64493406684823d0
       REAL*8, parameter :: rz3= 1.20205690315959d0
       REAL*8, parameter :: rz4= 1.08232323371114d0
       REAL*8, parameter :: rz5= 1.03692775514337d0
       REAL*8, parameter :: eg= 0.5772156649d0
       REAL*8, parameter :: ln_pi= 0.114472988584940016d1
      !END MODULE VBFNLO_HIGGSTO_riemann
!
!-----------------------------------------------------------------------------------------
!
      !MODULE VBFNLO_HIGGSTO_bernoulli
       REAL*8 :: b_num(0:18)=(/1.d0,-5.d-1,1.66666666666666d-1, &
     &  0.d0,-3.33333333333333d-2,0.d0,2.38095238095238d-2, &
     &  0.d0,-3.33333333333333d-2,0.d0,7.57575757575757d-2, &
     &  0.d0,-2.53113553113553d-1,0.d0,1.16666666666666d0, &
     &  0.d0,-7.09215686274509d0,0.d0,5.49711779448621d1/)
      !END MODULE VBFNLO_HIGGSTO_bernoulli
!
!-----------------------------------------------------------------------------------------
!
      !MODULE VBFNLO_HIGGSTO_ferbos
       INTEGER ifb
      !END MODULE VBFNLO_HIGGSTO_ferbos
!
!-----------------------------------------------------------------------------------------
!
      !MODULE VBFNLO_HIGGSTO_aux_hbb
       REAL*8 xb,xtop
       REAL*8, dimension(2) :: cxw,clw,ccts,csts,cctq,cswt,cdwt, &
     &                          cdzt,cxws,cxwc,cxz,rcmw,cctvi
      !END MODULE VBFNLO_HIGGSTO_aux_hbb
!
!-----------------------------------------------------------------------------------------
!
      !MODULE VBFNLO_HIGGSTO_optcp
       CHARACTER (len=3) ocp
      !END MODULE VBFNLO_HIGGSTO_optcp
!
!-----------------------------------------------------------------------------------------
!
      !MODULE VBFNLO_HIGGSTO_set_phys_const
       REAL*8 alpha_0,g_f,g_w,als,alpha_s_mh,alpha_s_kmh,imt
      !END MODULE VBFNLO_HIGGSTO_set_phys_const
!
!-----------------------------------------------------------------------------------------
!
      !MODULE VBFNLO_HIGGSTO_rootW
       INTEGER inc
      !END MODULE VBFNLO_HIGGSTO_rootW
!  
!-----------------------------------------------------------------------------------------
!
      !MODULE VBFNLO_HIGGSTO_transfmh
      REAL*8 tmuh
      !END MODULE VBFNLO_HIGGSTO_transfmh
!
!-----------------------------------------------------------------------------------------
!

!-----------------------------------------------------------------------
!--   Comments to Graeme Watt <watt(at)hep.ucl.ac.uk>
!----------------------------------------------------------------------
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
!--    FR2 = 1.d0                ! ratio of mu_f^2 to mu_r^2
!--    MUR = 1.d0                ! input mu_r in GeV
!--    ASMUR = 0.5d0             ! input value of alpha_s at mu_r
!--    MC = 1.4d0                ! charm quark mass
!--    MB = 4.75d0               ! bottom quark mass
!--    MT = 1.D10                ! top quark mass
!--    CALL VBFNLO_HIGGSTO_INITALPHAS(IORD, FR2, MUR, ASMUR, MC, MB, MT)
!--
!--   Then get alpha_s at a renormalisation scale mu_r with:
!--
!--    MUR = 100.d0              ! renormalisation scale in GeV
!--    ALFAS = VBFNLO_HIGGSTO_ALPHAS(MUR)
!--
!
      !MODULE VBFNLO_HIGGSTO_DZpar
      INTEGER IORDc
      REAL*8 FR2c,MURc,ASMURc,MCc,MBc,MTc,R0c
      !END MODULE VBFNLO_HIGGSTO_DZpar
!
!      MODULE VBFNLO_HIGGSTO_DZEROX
!      INTEGER IORD
!      REAL*8 FR2,MUR,ASMUR,MC,MB,MT,R0
!      END MODULE VBFNLO_HIGGSTO_DZEROX
!
      !MODULE VBFNLO_HIGGSTO_RZETA
      REAL*8, dimension(6) ::  ZETA
      !END MODULE VBFNLO_HIGGSTO_RZETA
      !MODULE VBFNLO_HIGGSTO_COLOUR 
      REAL*8 CF, CA, TR
      !END MODULE VBFNLO_HIGGSTO_COLOUR 
      !MODULE VBFNLO_HIGGSTO_ASINP  
      REAL*8 AS0, M20
      !END MODULE VBFNLO_HIGGSTO_ASINP  
      !MODULE VBFNLO_HIGGSTO_ASPAR  
      INTEGER NAORD, NASTPS
      !END MODULE VBFNLO_HIGGSTO_ASPAR  
      !MODULE VBFNLO_HIGGSTO_VARFLV 
      INTEGER IVFNS
      !END MODULE VBFNLO_HIGGSTO_VARFLV 
      !MODULE VBFNLO_HIGGSTO_NFFIX  
      INTEGER NFF
      !END MODULE VBFNLO_HIGGSTO_NFFIX  
      !MODULE VBFNLO_HIGGSTO_FRRAT  
      REAL*8 LOGFR
      !END MODULE VBFNLO_HIGGSTO_FRRAT  
      !MODULE VBFNLO_HIGGSTO_ASFTHR 
      REAL*8 ASC,M2C,ASB,M2B,AST,M2T
      !END MODULE VBFNLO_HIGGSTO_ASFTHR 
      !MODULE VBFNLO_HIGGSTO_BETACOM   
      REAL*8, dimension(3:6) :: BETA0,BETA1,BETA2,BETA3 
      !END MODULE VBFNLO_HIGGSTO_BETACOM   
!
!      MODULE VBFNLO_HIGGSTO_nmlo_alphas
!      CONTAINS
!
!-----------------------------------------------------------------------

      !MODULE VBFNLO_HIGGSTO_common_niels     
       REAL*8, dimension(3,0:15) :: plr
       REAL*8, dimension(3,0:15) :: plr_4
      !END MODULE VBFNLO_HIGGSTO_common_niels
!
!-----------------------------------------------------------------------------------------
!
      !MODULE VBFNLO_HIGGSTO_kountAC
       INTEGER kp,km
      !END MODULE VBFNLO_HIGGSTO_kountAC
!
!-----------------------------------------------------------------------------------------
!
      
      !MODULE VBFNLO_HIGGSTO_Solve_NonLin

! Corrections to FUNCTION Enorm - 28 November 2003

      ! IMPLICIT NONE
      PUBLIC  :: VBFNLO_HIGGSTO_hbrd,VBFNLO_HIGGSTO_hybrd


      !MODULE VBFNLO_HIGGSTO_acmplx_pro
      INTERFACE OPERATOR ( .CP. )
      MODULE PROCEDURE cp
      END INTERFACE

      INTERFACE OPERATOR ( .CQ. )
      MODULE PROCEDURE cq
      END INTERFACE

      INTERFACE OPERATOR ( .LCC. )
      MODULE PROCEDURE lcc
      END INTERFACE

      INTERFACE OPERATOR ( .RLNZ. )
      MODULE PROCEDURE rlnz
      END INTERFACE

      INTERFACE OPERATOR ( .ILNZ. )
      MODULE PROCEDURE ilnz
      END INTERFACE

      INTERFACE OPERATOR ( .RLN. )
      MODULE PROCEDURE rln
      END INTERFACE

      INTERFACE OPERATOR ( .CLN. )
      MODULE PROCEDURE cln
      END INTERFACE

      INTERFACE OPERATOR ( .FLN. )
      MODULE PROCEDURE fln
      END INTERFACE

      INTERFACE OPERATOR ( .LNSRS. )
      MODULE PROCEDURE lnsrs
      END INTERFACE

      INTERFACE OPERATOR ( .CRSRS. )
      MODULE PROCEDURE crsrs
      END INTERFACE

      INTERFACE OPERATOR ( .CR. )
      MODULE PROCEDURE cr
      END INTERFACE

      INTERFACE OPERATOR ( .CRZ. )
      MODULE PROCEDURE crz
      END INTERFACE

      INTERFACE OPERATOR ( .LC. )
      MODULE PROCEDURE lc
      END INTERFACE


      CONTAINS
          
       FUNCTION cp(x,y)
       !USE VBFNLO_HIGGSTO_riemann
       !USE VBFNLO_HIGGSTO_units
       IMPLICIT NONE
       REAL*8, dimension(2) :: cp,z1,z2
       REAL*8, INTENT(IN), dimension(2) :: x,y
       IF(abs(x(2)).eq.1.d0.and.abs(y(2)).eq.1.d0) THEN
        cp(1)= x(1)*y(1)
        cp(2)= 0.d0
       ELSE
        IF(abs(x(2)).eq.1.d0) THEN
         z1(1)= x(1)
         z1(2)= 0.d0
        ELSE
         z1= x
        ENDIF
        IF(abs(y(2)).eq.1.d0) THEN
         z2(1)= y(1)
         z2(2)= 0.d0
        ELSE
         z2= y        
        ENDIF
        cp(1)= z1(1)*z2(1)-z1(2)*z2(2)
        cp(2)= z1(1)*z2(2)+z1(2)*z2(1)
       ENDIF
       END FUNCTION cp


       FUNCTION cq( x,y )
       IMPLICIT NONE
       REAL*8, dimension(2) :: cq,z1,z2
       REAL*8, INTENT(IN), dimension(2) :: x,y
       REAL*8 theta
       IF(x(2).eq.0.d0.and.y(2).eq.0.d0) THEN
        cq(1)= x(1)/y(1)
        cq(2)= 0.d0
        RETURN
       ENDIF
       IF(abs(x(2)).eq.1.d0.and.abs(y(2)).eq.1.d0) THEN
        cq(1)= x(1)/y(1)
        cq(2)= 0.d0
       ELSE
        IF(abs(x(2)).eq.1.d0) then
         z1(1)= x(1)
         z1(2)= 0.d0
         z2(1)= y(1)
         z2(2)= y(2)         
        ELSEIF(abs(y(2)).eq.1.d0) then
         z1(1)= x(1)
         z1(2)= x(2)
         z2(1)= y(1)
         z2(2)= 0.d0         
        ELSE
         z1= x
         z2= y
        ENDIF
        IF(abs(z2(1)) > abs(z2(2))) THEN
         theta= z2(2)/z2(1)
         cq(1)= (z1(1)+theta*z1(2))/(theta*z2(2)+z2(1))
         cq(2)= (z1(2)-theta*z1(1))/(theta*z2(2)+z2(1))
        ELSE
         theta= z2(1)/z2(2)
         cq(1)= (theta*z1(1)+z1(2))/(theta*z2(1)+z2(2))
         cq(2)= (theta*z1(2)-z1(1))/(theta*z2(1)+z2(2)) 
        ENDIF     
       ENDIF
       END FUNCTION cq

      !END MODULE VBFNLO_HIGGSTO_acmplx_pro
!
!------------------------------------------------------------------
!
      !MODULE VBFNLO_HIGGSTO_acmplx_rat
      !END MODULE VBFNLO_HIGGSTO_acmplx_rat
!
!------------------------------------------------------------------
!
      !MODULE VBFNLO_HIGGSTO_linear_comb_c
       FUNCTION lcc(c,x)
       !USE VBFNLO_HIGGSTO_acmplx_pro
       !USE VBFNLO_HIGGSTO_units
       IMPLICIT NONE
       REAL*8, dimension(:) :: c
       REAL*8, dimension(2) :: x,lcc,sum
       INTENT(IN) c,x
       INTEGER n,i
       n= size(c)
       sum= c(n)*co
       DO i=n-1,1,-1
        sum= (sum.cp.x)+c(i)*co
       ENDDO
       lcc= sum      
       END FUNCTION lcc
      !END MODULE VBFNLO_HIGGSTO_linear_comb_c
!
!-----------------------------------------------------------------------------------------
!
      !MODULE VBFNLO_HIGGSTO_real_lnz
      FUNCTION rlnz(x,y)
      IMPLICIT NONE
      REAL*8 rlnz
      REAL*8, INTENT(IN) :: x,y
      rlnz= 0.5d0*log(x*x+y*y)
      END FUNCTION rlnz
      !END MODULE VBFNLO_HIGGSTO_real_lnz
!
!-----------------------------------------------------------------------------------------
!
      !MODULE VBFNLO_HIGGSTO_imag_lnz
      FUNCTION ilnz(x,y)
      !USE VBFNLO_HIGGSTO_riemann
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: x,y
      REAL*8 teta,tnteta,sr,si,ax,ay,ilnz
      ax= abs(x)
      ay= abs(y)
      IF(x.eq.0.d0) THEN
       IF(y > 0.d0) THEN
          teta= pi/2.d0
       ELSE
          teta= -pi/2.d0
       ENDIF
       ilnz= teta
       RETURN
      ELSE IF(y.eq.0.d0) THEN 
       IF(x > 0.d0) THEN
          teta= 0.d0
       ELSE
          teta= pi
       ENDIF
       ilnz= teta
       RETURN
      ELSE
       tnteta= ay/ax
       teta= atan(tnteta)
       sr= x/ax
       si= y/ay
       IF(sr > 0.d0) THEN
        ilnz= si*teta
       ELSE
        ilnz= si*(pi-teta)
       ENDIF
       RETURN
      ENDIF
      END FUNCTION ilnz
      !END MODULE VBFNLO_HIGGSTO_imag_lnz
!
!-----------------------------------------------------------------------------------------
!
      !MODULE VBFNLO_HIGGSTO_real_ln
      FUNCTION rln(x,ep)
      IMPLICIT NONE
      REAL*8 rln
      REAL*8, INTENT(IN) :: x,ep
      rln= log(abs(x))
      END FUNCTION rln
      !END MODULE VBFNLO_HIGGSTO_real_ln
!
!-----------------------------------------------------------------------------------------
!
      !MODULE VBFNLO_HIGGSTO_cmplx_ln
       FUNCTION cln(x,ep)
       !USE VBFNLO_HIGGSTO_riemann
       IMPLICIT NONE
       REAL*8, dimension(2) :: cln
       REAL*8, INTENT(IN) :: x,ep
       IF(abs(ep).ne.1.d0) THEN
        print*,' Wrong argument for CLN '
        STOP
       ENDIF
       cln(1)= log(abs(x))
       IF(x > 0.d0) THEN
        cln(2)= 0.d0
       ELSE
        cln(2)= ep*pi
       ENDIF 
       END FUNCTION cln
      !END MODULE VBFNLO_HIGGSTO_cmplx_ln
!
!-----------------------------------------------------------------------------------------
!
      !MODULE VBFNLO_HIGGSTO_full_ln
       FUNCTION fln(x,y)
       !USE VBFNLO_HIGGSTO_riemann
       IMPLICIT NONE
       REAL*8 ax,ay,teta,tnteta,sr,si
       REAL*8, INTENT(IN) :: x,y
       REAL*8, dimension(2) :: fln
       IF(abs(y).eq.1.d0) THEN
        fln(1)= log(abs(x))
        IF(x > 0.d0) THEN
         fln(2)= 0.d0
        ELSE
         fln(2)= y*pi
        ENDIF 
       ELSE
        fln(1)= 0.5d0*log(x*x+y*y)
        ax= abs(x)
        ay= abs(y)
        IF(x.eq.0.d0) THEN
         IF(y > 0.d0) THEN
          teta= pi/2.d0
         ELSE
          teta= -pi/2.d0
         ENDIF
         fln(2)= teta
        ELSE IF(y.eq.0.d0) THEN 
         IF(x > 0.d0) THEN
          teta= 0.d0
         ELSE
          teta= pi
         ENDIF
         fln(2)= teta
        ELSE
         tnteta= ay/ax
         teta= atan(tnteta)
         sr= x/ax
         si= y/ay
         IF(sr > 0.d0) THEN
          fln(2)= si*teta
         ELSE
          fln(2)= si*(pi-teta)
         ENDIF
        ENDIF
       ENDIF
       END FUNCTION fln
      !END MODULE VBFNLO_HIGGSTO_full_ln
!
!------------------------------------------------------------------------
!
      !MODULE VBFNLO_HIGGSTO_ln_2_riemann
       FUNCTION lnsrs(x,y)
       !USE VBFNLO_HIGGSTO_riemann
       !USE VBFNLO_HIGGSTO_full_ln
       !USE VBFNLO_HIGGSTO_units
       IMPLICIT NONE
       REAL*8, dimension(2) :: lnsrs,x,y,olnsrs
       REAL*8 ax,ay,teta,tnteta,ilnx,sr,si,xms
       INTENT(IN) x,y
!
       IF(abs(y(2)).ne.1.d0) THEN
!
        xms= x(1)*x(1)+x(2)*x(2)
        IF(abs(1.d0-sqrt(xms)).lt.1.d-12) THEN
         teta= atan(abs(x(2)/x(1)))
         sr= x(1)/abs(x(1))
         si= x(2)/abs(x(2))
         lnsrs(1)= 0.d0 
         IF(sr > 0.d0) THEN
          lnsrs(2)= si*teta
         ELSE
          lnsrs(2)= si*(pi-teta)
         ENDIF
         RETURN
        ELSE 
         lnsrs= x(1).fln.x(2)
         olnsrs= y(1).fln.y(2)
         RETURN
        ENDIF
       ENDIF 
!
       lnsrs(1)= 0.5d0*log(x(1)*x(1)+x(2)*x(2))
       ax= abs(x(1))
       ay= abs(x(2))
       IF(x(1).eq.0.d0) THEN
        IF(x(2) > 0.d0) THEN
           teta= pi/2.d0
        ELSE
           teta= -pi/2.d0
        ENDIF
        ilnx= teta
       ELSE IF(x(2).eq.0.d0) THEN 
        IF(x(1) > 0.d0) THEN
         teta= 0.d0
        ELSE
         teta= pi
        ENDIF
        ilnx= teta
       ELSE
        tnteta= ay/ax
        teta= atan(tnteta)
        sr= x(1)/ax
        si= x(2)/ay
        IF(sr > 0.d0) THEN
         ilnx= si*teta
        ELSE
         ilnx= si*(pi-teta)
        ENDIF
       ENDIF
!
       IF(x(1) > 0.d0) THEN
        IF(y(1) < 0.d0) THEN
         IF((x(2) > 0.d0.and.y(2) > 0).or. &
     &      (x(2) < 0.d0.and.y(2) < 0)) THEN
          lnsrs(2)= ilnx
         ELSEIF(x(2) > 0.d0.and.y(2) < 0) THEN
          lnsrs(2)= ilnx-pi
         ELSE
          print*,'+++++++++++++++++++++++++++++++++++'
          print*,' anomaly ln '
          print*,x
          print*,y
          print*,'+++++++++++++++++++++++++++++++++++'
         ENDIF
        ELSE
         lnsrs(2)= ilnx
        ENDIF
       ELSE
        IF(abs(y(2)).ne.1.d0) THEN
         print*,' Wrong argument for LNSRS '
         STOP
        ELSE
!
         IF((x(2)*y(2)) < 0.d0) THEN
          lnsrs(2)= ilnx+2.d0*y(2)*pi
         ELSE
          lnsrs(2)= ilnx
         ENDIF
!
        ENDIF
       ENDIF       
       RETURN
       END FUNCTION lnsrs
      !END MODULE VBFNLO_HIGGSTO_ln_2_riemann
!
!-----------------------------------------------------------------------------------------
!
      !MODULE VBFNLO_HIGGSTO_cmplx_srs_root
       FUNCTION crsrs(x,y)
       !USE VBFNLO_HIGGSTO_units
       !USE VBFNLO_HIGGSTO_riemann
       IMPLICIT NONE
       REAL*8 r,rs,u,v 
       REAL*8, dimension(2) :: crsrs
       REAL*8, INTENT(IN), dimension(2):: x,y
!
       IF(x(1) > 0.d0) THEN
        rs= x(1)*x(1)+x(2)*x(2)
        r= sqrt(rs)
        u= 0.5d0*(r+x(1))
        crsrs(1)= sqrt(u)
        r= sqrt(rs)
        v= 0.5d0*(r-x(1))
        crsrs(2)= sign(one,x(2))*sqrt(v)
       ELSE
        rs= x(1)*x(1)+x(2)*x(2)
        r= sqrt(rs)
        u= 0.5d0*(r+x(1))
        v= 0.5d0*(r-x(1))
!
        IF((x(2)*y(2)) > 0.d0) THEN
         crsrs(1)= sqrt(u)
         crsrs(2)= sign(one,x(2))*sqrt(v)
        ELSEIF((x(2)*y(2)) < 0.d0) THEN
         crsrs(1)= -sqrt(u)
         crsrs(2)= sign(one,x(2))*sqrt(v)
        ENDIF
!
       ENDIF
       END FUNCTION crsrs
      !END MODULE VBFNLO_HIGGSTO_cmplx_srs_root
!
!-----------------------------------------------------------------------------------------
!
      !MODULE VBFNLO_HIGGSTO_cmplx_root
       FUNCTION cr(x,ep)
       IMPLICIT NONE
       REAL*8, dimension(2) :: cr
       REAL*8, INTENT(IN) :: x,ep
       IF(x > 0.d0) THEN 
        cr(1)= sqrt(x)
        cr(2)= 0.d0       
       ELSE
        cr(1)= 0.d0
        cr(2)= ep*sqrt(abs(x))
       ENDIF
       END FUNCTION cr
      !END MODULE VBFNLO_HIGGSTO_cmplx_root  
!
!-----------------------------------------------------------------------------------------
!
      !MODULE VBFNLO_HIGGSTO_cmplx_rootz
       FUNCTION crz(x,y)
       !USE VBFNLO_HIGGSTO_units
       IMPLICIT NONE
       REAL*8, dimension(2) :: crz
       REAL*8, INTENT(IN) :: x,y
       REAL*8 rs,r,u,v,sax
       IF(abs(y/x).lt.1.d-8) THEN
        IF(x.gt.0.d0) THEN
         crz(1)= sqrt(x)*(1.d0+0.5d0*y*y/(x*x))
         crz(2)= 0.5d0*y/sqrt(x)
        ELSEIF(x.lt.0.d0) THEN
         sax= sqrt(abs(x))
         crz(1)= 0.5d0*abs(y)/sax
         crz(2)= sign(one,y)*sax*(1.d0+0.5d0*y*y/(x*x))
        ENDIF
       ELSE
        rs= x*x+y*y
        r= sqrt(rs)
        u= 0.5d0*(r+x)
        crz(1)= sqrt(u)
        IF(x > 0.d0.and.y.eq.0.d0) THEN
         crz(2)= 0.d0
        ELSE
         v= 0.5d0*(r-x)
         crz(2)= sign(one,y)*sqrt(v)
        ENDIF
       ENDIF
       END FUNCTION crz
      !END MODULE VBFNLO_HIGGSTO_cmplx_rootz
!
!-----------------------------------------------------------------------------------------
!
      !MODULE VBFNLO_HIGGSTO_linear_comb
       FUNCTION lc(c,x)
       IMPLICIT NONE
       REAL*8, dimension(:) :: c
       REAL*8 x,lc
       INTENT(IN) c,x
       INTEGER n,i
       REAL*8 sum      
       n= size(c)
       sum= c(n)
       DO i=n-1,1,-1
        sum= sum*x+c(i)
       ENDDO
       lc= sum      
       END FUNCTION lc
      !END MODULE VBFNLO_HIGGSTO_linear_comb
!
!-----------------------------------------------------------------------------------------
!
      !MODULE VBFNLO_HIGGSTO_sp_fun
!
! extension to n+p = 4
!
! 1 Li_2
! 2 Li_3
! 3 S_12
! 4 Li_4
! 5 S_22
! 6 S_13
!
      RECURSIVE FUNCTION VBFNLO_HIGGSTO_s_niels_up4(x) RESULT(res)
      !USE VBFNLO_HIGGSTO_riemann
      !USE VBFNLO_HIGGSTO_acmplx_pro
      !USE VBFNLO_HIGGSTO_cmplx_ln
      !USE VBFNLO_HIGGSTO_real_ln
      !USE VBFNLO_HIGGSTO_common_niels
      !USE VBFNLO_HIGGSTO_linear_comb
      !USE VBFNLO_HIGGSTO_units
!
      IMPLICIT NONE
!
      REAL*8, dimension(6,2) :: res
      REAL*8, dimension(2) :: x
      INTENT(IN) x
!
      INTEGER i,j,l
      REAL*8 ym,ym2,ymod,l2,xmo,z_exp,sum,s13m1,s22m1
      REAL*8, dimension(6) :: sign
      REAL*8, dimension(0:15) :: sumr
      REAL*8, dimension(6,2) :: aniels,bniels,add
      REAL*8, dimension(2) :: ln_omx,ln2_omx,ln3_omx,ln_omy,ln2_omy, &
     &  ln3_omy,ln_y,ln2_y,ln3_y,ln_yi,ln_myi,ln2_myi,ln3_myi,ln4_y, &
     &  y,omy,yi,yt,omx,myi,prod,li2,li3,s12,prod1,prod2,prodl, &
     &  prodl1,prodl2,add1,add2,add3,add4,add5,add6,add7,add8,add9, &
     &  add10,ln_my,ln2_my,ln3_my,ln4_my,add11,add12,add13,add14, &
     &  add15,add16,add17,add18,add19,add20,add21,ln4_omx,ln_ymo, &
     &  ln2_ymo,ln3_ymo
      REAL*8, parameter :: li4h=6.72510831970049127d-2
!
      IF(abs(x(2)).ne.1.d0) THEN
       print 1,x
       stop
      ENDIF
    1 format(' wrong input for niels ',2e20.5)
!
      l2= log(2.d0)
      s13m1= li4h+l2*(7.d0/8.d0*rz3-rz2*l2/4.d0+l2*l2*l2/24.d0)-rz4
      s22m1= 2.d0*li4h+l2*(7.d0/4.d0*rz3-0.5d0*rz2*l2+l2*l2*l2/12.d0)- &
     &       15.d0/8.d0*rz4
!
      IF(x(1).eq.1.d0) THEN
       res(1,1)= rz2
       res(2,1)= rz3
       res(3,1)= rz3
       res(4,1)= rz4
       res(5,1)= -1/2*(rz2**2-3.d0*rz4)
       res(6,1)= rz4 
       FORALL (i=1:6) res(i,2)= 0.d0
       RETURN
      ENDIF
      IF(x(1).eq.0.5d0) THEN
       res(1,1)= 0.5d0*(rz2-l2*l2)
       res(2,1)= 7.d0/8.d0*rz3-l2*(0.5d0*rz2-l2*l2/6.d0)
       res(3,1)= rz3/8.d0-l2*l2*l2/6.d0
       res(4,1)= li4h
       res(5,1)= -l2*res(3,1)-l2**4/8.d0-rz2*rz2/4.d0+3.d0/4.d0*rz4
       res(6,1)= -res(4,1)+l2*res(3,1)+l2**2/4.d0*(1.d0/3.d0+rz2)- &
     &           l2*rz3-rz4
       FORALL (i=1:6) res(i,2)= 0.d0
       RETURN
      ENDIF
!
! gets rid of arguments with negative real part
!
      IF(x(1) < 0.d0) THEN
       xmo= x(1)-1.d0
       ymod= xmo*xmo
       y(1)= x(1)*xmo/ymod
       y(2)= -x(2)
       omx= co-x
       ln_omx= omx(1).cln.omx(2)
       ln2_omx= ln_omx.cp.ln_omx
       ln3_omx= ln2_omx.cp.ln_omx
       ln4_omx= ln3_omx.cp.ln_omx
       aniels= VBFNLO_HIGGSTO_s_niels_up4(y) 
       li2(1:2)= aniels(1,1:2)
       li3(1:2)= aniels(2,1:2)
       s12(1:2)= aniels(3,1:2)
       add1= ln2_omx.cp.li2
       add2= ln_omx.cp.li3
       add3= ln_omx.cp.s12
       add(1,1:2)= -0.5d0*ln2_omx(1:2)
       sign(1)= -1.d0
       prod= ln_omx.cp.li2
       add(2,1:2)= aniels(3,1:2)-ln3_omx(1:2)/6.d0-prod(1:2)
       sign(2)= -1.d0
       add(3,1:2)= ln3_omx(1:2)/6.d0
       sign(3)= 1.d0
       add(4,1:2)= -0.5d0*add1(1:2)-add2(1:2)+add3(1:2) &
     &      -1.d0/24.d0*ln4_omx(1:2)-aniels(6,1:2)+aniels(5,1:2) &
     &      +co(1:2)*(2.d0*s13m1-s22m1+1.d0/8.d0*rz4)
       sign(4)= -1.d0
       add(5,1:2)= add3(1:2)+1.d0/24.d0*ln4_omx(1:2) &
     &      -2.d0*aniels(6,1:2)+aniels(5,1:2)+co(1:2)*(4.d0*s13m1 &
     &      -2.d0*s22m1+0.25d0*rz4)
       sign(5)= 1.d0
       add(6,1:2)= -1.d0/24.d0*ln4_omx(1:2)-aniels(6,1:2)+ &
     &       co(1:2)*(4.d0*s13m1-2.d0*s22m1+0.25d0*rz4)
       sign(6)= -1.d0
      ELSE
       y= x
       add= 0.d0
       sign= 1.d0
      ENDIF
!
      ym2= y(1)*y(1)
      ym= abs(y(1))
!
! |y| < 1 & Re(y) < 1/2
!
      IF(ym < 1.d0.and.y(1) < 0.5d0) THEN
       omy= co-y
       z_exp= -(omy(1).rln.omy(2))
       DO l=1,3
        sumr(0:15)= plr(l,0:15)
        sum= sumr.lc.z_exp
        res(l,1)= sum
        res(l,2)= 0.d0
       ENDDO
       DO l=4,6
        sumr(0:15)= plr_4(l-3,0:15)
        sum= sumr.lc.z_exp
        res(l,1)= sum
        res(l,2)= 0.d0
       ENDDO
       FORALL (i=1:6,j=1:2) res(i,j)= sign(i)*res(i,j)+add(i,j)
       RETURN
!
! |y| < 1 & Re(y) > 1/2
!
      ELSEIF(ym < 1.d0.and.y(1) > 0.5d0) THEN
       omy= co-y
       ln_omy= omy(1).cln.omy(2)
       ln_y= y(1).cln.y(2)
       ln2_omy= ln_omy.cp.ln_omy
       ln2_y= ln_y.cp.ln_y
       ln3_y= ln2_y.cp.ln_y
       aniels= VBFNLO_HIGGSTO_s_niels_up4(omy)
       s12(1:2)= aniels(3,1:2)
       li3(1:2)= aniels(2,1:2)
       li2(1:2)= aniels(1,1:2)
       prodl= ln_y.cp.ln_omy
       res(1,1:2)= -li2(1:2)-prodl(1:2)+rz2*co(1:2)
       prodl= ln_omy.cp.ln2_y
       prod= ln_y.cp.li2
       res(2,1:2)= -s12(1:2)-prod(1:2)-0.5d0*prodl(1:2)+ &
     &             rz2*ln_y(1:2)+rz3*co(1:2)
       prodl= ln_y.cp.ln2_omy
       prod= ln_omy.cp.li2
       res(3,1:2)= -li3(1:2)+prod(1:2)+0.5d0*prodl(1:2)+rz3*co(1:2)
       add1= ln2_y.cp.li2       
       add2= ln_y.cp.s12
       add3= ln_omy.cp.ln3_y
       add4= ln_omy.cp.s12
       add5= ln_y.cp.ln_omy
       add5= add5.cp.li2
       add6= ln_y.cp.li3
       add7= ln2_y.cp.ln2_omy
       add8= ln_omy.cp.li3
       add9= ln_y.cp.ln3_omy
       add10= ln2_omy.cp.li2
       res(4,1:2)= -aniels(6,1:2)-0.5d0*add1(1:2)-add2(1:2)- &
     &     add3(1:2)/6.d0+rz3*ln_y(1:2)+0.5d0*rz2*ln2_y+rz4*co(1:2)
       res(5,1:2)= aniels(5,1:2)-add4(1:2)-add5(1:2)+add6(1:2) &
     &    -add2(1:2)-0.25d0*add7(1:2)+0.5d0*co(1:2)*(rz2*rz2-3.d0*rz4)
       res(6,1:2)= aniels(4,1:2)-add8(1:2)+add9(1:2)/6.d0+ &
     &     0.5d0*add10(1:2)-co(1:2)*s13m1
       FORALL (i=1:6,j=1:2) res(i,j)= sign(i)*res(i,j)+add(i,j)
       RETURN
!
! |y| > 1 & Re(y^-1) < 1/2
!
      ELSEIF(ym > 1.d0.and.y(1) >= 2.d0) THEN
       yi(1)= y(1)/ym2
       yi(2)= -y(2)
       ln_yi= yi(1).cln.yi(2)
       ln_my= (-y(1)).cln.(-y(2))
       ln2_my= ln_my.cp.ln_my
       ln3_my= ln3_my.cp.ln_my
       ln4_my= ln3_my.cp.ln_my
       myi(1)= -y(1)/ym2
       myi(2)= y(2)
       ln_myi= myi(1).cln.myi(2)
       ln2_myi= ln_myi.cp.ln_myi
       ln3_myi= ln2_myi.cp.ln_myi
       aniels= VBFNLO_HIGGSTO_s_niels_up4(yi)
       li2(1:2)= aniels(1,1:2)
       li3(1:2)= aniels(2,1:2)
       s12(1:2)= aniels(3,1:2)
       res(1,1:2)= -aniels(1,1:2)-0.5d0*ln2_myi(1:2)-rz2*co(1:2)
       res(2,1:2)= aniels(2,1:2)+ln3_myi(1:2)/6.d0+rz2*ln_myi(1:2)
       li2(1:2)= aniels(1,1:2)
       prod= ln_myi.cp.li2
       res(3,1:2)= aniels(2,1:2)-aniels(3,1:2)-ln3_myi(1:2)/6.d0- &
     &             prod(1:2)+rz3*co(1:2)
       res(4,1:2)= aniels(4,1:2)+0.5d0*rz2*ln2_my(1:2)+ &
     &     1.d0/24.d0*ln4_my+7.d0/4.d0*rz4*co(1:2)
       add1= ln_my.cp.li3
       add1= ln_my.cp.s12
       add3= ln2_my.cp.li2
       res(5,1:2)= -aniels(5,1:2)+2.d0*aniels(4,1:2)+add1(1:2)- &
     &     rz3*ln_my-ln4_my/24.d0+7.d0/4.d0*rz4*co(1:2)
       res(6,1:2)= aniels(4,1:2)+aniels(6,1:2)-aniels(5,1:2)+ &
     &     add1(1:2)-add2(1:2)+0.5d0*add3(1:2)+ln4_my/24.d0+ &
     &     co(1:2)*(-2.d0*s13m1+s22m1+7.d0/8.d0*rz4)
       FORALL (i=1:6,j=1:2) res(i,j)= sign(i)*res(i,j)+add(i,j)
       RETURN
!
! |y| > 1 & Re(y^-1) > 1/2
!
      ELSEIF(ym > 1.d0.and.y(1) < 2.d0) THEN
       yt(1)= (y(1)-1.d0)/y(1)
       yt(2)= y(2)
       omy= co-y
       ln_y= y(1).cln.y(2)
       ln_omy= omy(1).cln.omy(2)
       ln2_y= ln_y.cp.ln_y
       ln3_y= ln2_y.cp.ln_y
       ln4_y= ln3_y.cp.ln_y
       ln2_omy= ln_omy.cp.ln_omy
       ln_my= (-y(1)).cln.(-y(2))
       ln2_my= ln_my.cp.ln_my
       ln3_my= ln2_my.cp.ln_my
       ln4_my= ln3_my.cp.ln_my
       aniels= VBFNLO_HIGGSTO_s_niels_up4(yt)
       bniels= 0.d0
       prod= ln_y.cp.ln_omy
       bniels(1,1:2)= aniels(1,1:2)-prod(1:2)+0.5d0*ln2_y(1:2)+ &
     &                rz2*co(1:2)
       res(1,1:2)= aniels(1,1:2)-prod(1:2)+0.5d0*ln2_y(1:2)+ &
     &             rz2*co(1:2)
       prodl= ln2_y.cp.ln_omy
       li2(1:2)= bniels(1,1:2)
       prod= ln_y.cp.li2  
       res(2,1:2)= -aniels(3,1:2)+prod(1:2)+0.5d0*prodl(1:2)- &
     &             ln3_y(1:2)/6.d0+rz3*co(1:2)
       li2(1:2)= aniels(1,1:2)
       prodl1= ln_omy.cp.ln2_y
       prodl2= ln_y.cp.ln2_omy
       prod1= ln_omy.cp.li2
       prod2= ln_y.cp.li2
       res(3,1:2)= rz3*co(1:2)-0.5d0*prodl1(1:2)-prod1(1:2)+ &
     &             0.5d0*prodl2(1:2)+ln3_y(1:2)/6.d0+prod2(1:2)+ &
     &             aniels(2,1:2)-aniels(3,1:2)
       add1= ln_ymo.cp.ln3_y
       add2= ln_y.cp.s12
       add3= ln2_y.cp.li2
       add4= ln_ymo.cp.s12
       add5= ln_my.cp.ln_ymo
       add5= add5.cp.ln2_y
       add6= ln_my.cp.s12
       add7= ln_y.cp.ln_ymo
       add7= add7.cp.li2
       add8= ln_my.cp.ln3_y
       add9= ln_y.cp.ln_my
       add9= add9.cp.li2
       add10= ln_y.cp.ln_my
       add11= ln_y.cp.li3
       add12= ln2_y.cp.ln2_ymo
       add13= ln_ymo.cp.li3
       add14= ln_my.cp.ln_ymo
       add14= add14.cp.li2
       add15= ln_my.cp.li3
       add16= ln_y.cp.ln_ymo
       add16= add16.cp.ln2_my
       add17= ln_y.cp.ln_my
       add17= add17.cp.ln2_ymo    
       add18= ln_y.cp.ln3_ymo
       add19= ln2_ymo.cp.li2
       add20= ln2_my.cp.li2
       add21= ln2_y.cp.ln2_my
       res(4,1:2)= -1.d0/6.d0*add1(1:2)-add2(1:2)+ln_y(1:2)*rz3 &
     &  -0.5d0*ln2_my*rz2+0.5d0*add3(1:2)-0.5d0*ln2_y(1:2)*rz2 &
     &  -1.d0/24.d0*ln4_my(1:2)+1.d0/6.d0*ln4_y(1:2)+aniels(6,1:2) &
     &  -11.d0/4.d0*rz4*co(1:2)
       res(5,1:2)= -5.d0/6.d0*add1(1:2)+add4(1:2)+0.5d0*add5(1:2) &
     &  -0.5d0*add8(1:2)+add6(1:2)-add7(1:2)-add9(1:2) &
     &  +add10(1:2)*rz2+add11(1:2)-3.d0*add2(1:2)+ln_y(1:2)*rz3 &
     &  +1.d0/4.d0*add12(1:2)+2.d0*add3(1:2)-ln2_y(1:2)*rz2 &
     &  +1.d0/24.d0*ln4_my(1:2)+7.d0/12.d0*ln4_y(1:2) &
     &  +2.d0*aniels(6,1:2)-aniels(5,1:2)-7.d0/2.d0*rz4*co(1:2)
       res(6,1:2)= -7.d0/6.d0*add1(1:2)-add13(1:2)+add4(1:2) &
     &  +3.d0/2.d0*add5(1:2)+add14(1:2)-add8(1:2)-add15(1:2) &
     &  +add6(1:2)-0.5d0*add16(1:2)-2.d0*add7(1:2)-0.5d0*add17(1:2) &
     &  -2.d0*add9(1:2)+add10(1:2)*rz2-1.d0/6.d0*add18(1:2) &
     &  +2.d0*add11(1:2)-2.d0*add2(1:2)+0.5d0*add19(1:2) &
     &  +0.5d0*add20(1:2)-0.5d0*ln2_my*rz2+3.d0/4.d0*add12(1:2) &
     &  +0.5d0*add21(1:2)+2.d0*add3(1:2)-0.5d0*ln2_y(1:2)*rz2 &
     &  -1.d0/24.d0*ln4_my(1:2)+ 7.d0/12.d0*ln4_y(1:2)+aniels(4,1:2) &
     &  +aniels(6,1:2)-aniels(5,1:2)+co(1:2)*(2.d0*s13m1-s22m1 &
     &  -21.d0/8.d0*rz4)
       FORALL (i=1:6,j=1:2) res(i,j)= sign(i)*res(i,j)+add(i,j)
       RETURN
      ENDIF
!
      END FUNCTION VBFNLO_HIGGSTO_s_niels_up4
!
!------------------------------------------------------------------
!
      FUNCTION VBFNLO_HIGGSTO_li2_srsz(x,y,unit) RESULT(value)
      !USE VBFNLO_HIGGSTO_acmplx_pro
      !USE VBFNLO_HIGGSTO_full_ln
      !USE VBFNLO_HIGGSTO_riemann
      !USE VBFNLO_HIGGSTO_units
      !USE VBFNLO_HIGGSTO_kountAC
!
      IMPLICIT NONE
!
      INTEGER unit
      REAL*8 xms,atheta,theta,sr,si
      REAL*8, dimension (:) :: x,y
      REAL*8, dimension(3,2) :: aux
      REAL*8, dimension(2) :: lnc,add,value
      INTENT(IN) x,y
!
      IF(abs(y(2)).ne.1.d0) THEN
       IF(unit.eq.1) THEN
        xms= x(1)*x(1)+x(2)*x(2)
        IF(abs(1.d0-sqrt(xms)).gt.1.d-12) THEN
         print*,' apparent inconsistency '
        ENDIF
        atheta= atan(abs(x(2)/x(1)))
        sr= x(1)/abs(x(1))
        si= x(2)/abs(x(2))
        IF(sr > 0.d0.and.si > 0.d0) THEN
         theta= atheta
        ELSEIF(sr > 0.d0.and.si < 0.d0) THEN
         theta= 2.d0*pi-atheta
        ELSEIF(sr < 0.d0.and.si > 0.d0) THEN
         theta= pi-atheta
        ELSEIF(sr < 0.d0.and.si < 0.d0) THEN
         theta= pi+atheta
        ENDIF
         aux= VBFNLO_HIGGSTO_poly_unit(theta)
         value(1:2)= aux(1,1:2)
       ELSE 
        value= VBFNLO_HIGGSTO_li2(x)
       ENDIF 
       RETURN 
      ENDIF
!
      IF(x(1) > 1.d0) THEN
       IF(y(1) < 1.d0) THEN
        value= VBFNLO_HIGGSTO_li2(x)
        print*,'+++++++++++++++++++++++++++++++++++'
        print*,' anomaly Li2 '
        print*,x
        print*,y
        print*,'+++++++++++++++++++++++++++++++++++'
       ENDIF
       lnc= x(1).fln.x(2)
       add(1)= -lnc(2)
       add(2)= lnc(1)
       IF(y(2) < 0.d0.and.x(2) > 0.d0) THEN
        value= VBFNLO_HIGGSTO_li2(x)-2.d0*pi*add
        km= km+1
       ELSEIF(y(2) > 0.d0.and.x(2) < 0.d0) THEN
        value= VBFNLO_HIGGSTO_li2(x)+2.d0*pi*add
        kp= kp+1
       ELSE
        value= VBFNLO_HIGGSTO_li2(x)
       ENDIF
      ELSE
       value= VBFNLO_HIGGSTO_li2(x)
      ENDIF
!       
      RETURN
!
      END FUNCTION VBFNLO_HIGGSTO_li2_srsz
!
!------------------------------------------------------------------
!
      FUNCTION VBFNLO_HIGGSTO_li3_srsz(x,y,unit) RESULT(value)
      !USE VBFNLO_HIGGSTO_acmplx_pro
      !USE VBFNLO_HIGGSTO_full_ln
      !USE VBFNLO_HIGGSTO_riemann
      !USE VBFNLO_HIGGSTO_units
!
      IMPLICIT NONE
!
      INTEGER unit
      REAL*8 xms,atheta,theta,sr,si
      REAL*8, dimension (:) :: x,y
      REAL*8, dimension(3,2) :: aux
      REAL*8, dimension(2) :: lnc,lncs,add,value
      INTENT(IN) x,y
!
      IF(abs(y(2)).ne.1.d0) THEN
       IF(unit.eq.1) THEN
        xms= x(1)*x(1)+x(2)*x(2)
        IF(abs(1.d0-sqrt(xms)).gt.1.d-12) THEN
         print*,' apparent inconsistency '
        ENDIF
        atheta= atan(abs(x(2)/x(1)))
        sr= x(1)/abs(x(1))
        si= x(2)/abs(x(2))
        IF(sr > 0.d0.and.si > 0.d0) THEN
         theta= atheta
        ELSEIF(sr > 0.d0.and.si < 0.d0) THEN
         theta= 2.d0*pi-atheta
        ELSEIF(sr < 0.d0.and.si > 0.d0) THEN
         theta= pi-atheta
        ELSEIF(sr < 0.d0.and.si < 0.d0) THEN
         theta= pi+atheta
        ENDIF
         aux= VBFNLO_HIGGSTO_poly_unit(theta)
         value(1:2)= aux(2,1:2)
       ELSE 
        value= VBFNLO_HIGGSTO_li3(x)
       ENDIF 
       RETURN 
      ENDIF
!
      IF(x(1) > 1.d0) THEN
       IF(y(1) < 1.d0) THEN
        print*,' anomaly'
        print*,x
        print*,y
        STOP
       ENDIF
       lnc= x(1).fln.x(2)
       lncs= lnc.cp.lnc
       add(1)= -lncs(2)
       add(2)= lncs(1)
       IF(y(2) < 0.d0.and.x(2) > 0.d0) THEN
        value= VBFNLO_HIGGSTO_li3(x)-pi*add
       ELSEIF(y(2) > 0.d0.and.x(2) < 0.d0) THEN
        value= VBFNLO_HIGGSTO_li3(x)+pi*add
       ELSE
        value= VBFNLO_HIGGSTO_li3(x)
       ENDIF
      ELSE
       value= VBFNLO_HIGGSTO_li3(x)
      ENDIF
!       
      RETURN
!
      END FUNCTION VBFNLO_HIGGSTO_li3_srsz
!
!
!------------------------------------------------------------------
!
      FUNCTION VBFNLO_HIGGSTO_poly_unit(theta) RESULT(value)
      !USE VBFNLO_HIGGSTO_riemann
      !USE VBFNLO_HIGGSTO_full_ln
      !USE VBFNLO_HIGGSTO_real_lnz
      !USE VBFNLO_HIGGSTO_imag_lnz
      !USE VBFNLO_HIGGSTO_linear_comb_c
      !USE VBFNLO_HIGGSTO_acmplx_pro
      !USE VBFNLO_HIGGSTO_bernoulli
      !USE VBFNLO_HIGGSTO_units
!
      IMPLICIT NONE
!
      INTEGER n
      REAL*8 theta
      REAL*8, dimension(2) :: arg,ln,z,sum2,sum3,sum4,pw2,pw3,pw4, &
     &        cb2,cb3,cb4
      REAL*8, dimension(4) :: psi
      REAL*8, dimension(0:15) :: rzf
      REAL*8, dimension(3,2) :: value
      REAL*8, dimension(0:17) :: c2,c3,c4,fac
!
      data rzf(0:15) /-0.5d0,0.d0,1.644934066848226d0, &
     &     1.202056903159594d0,1.082323233711138d0, &
     &     1.036927755143370d0,1.017343061984449d0, &
     &     1.008349277381923d0,1.004077356197944d0, &
     &     1.002008392826082d0,1.000994575127818d0, &
     &     1.000494188604119d0,1.000246086553308d0, &
     &     1.000122713346578d0,1.000061248135059d0, &
     &     1.000030588236307d0/
!
      psi(1)= eg
      psi(2)= psi(1)+1.d0
      psi(3)= psi(2)+0.5d0
      psi(4)= psi(3)+1.d0/3.d0
!
      fac(0)= 1.d0
      DO n=1,17
       fac(n)= n*fac(n-1)
      ENDDO
!
      arg(1)= 0.d0
      arg(2)= -theta
!
      ln= arg(1).fln.arg(2)
!
      DO n=0,17
       IF((2-n) >= 0) THEN
        c2(n)= rzf(2-n)/fac(n)
       ELSE
        c2(n)= -b_num(n-1)/(n-1)/fac(n)
       ENDIF
       IF((3-n) >= 0) THEN
        c3(n)= rzf(3-n)/fac(n)
       ELSE
        c3(n)= -b_num(n-2)/(n-2)/fac(n)
       ENDIF
       IF((4-n) >= 0) THEN
        c4(n)= rzf(4-n)/fac(n)
       ELSE
        c4(n)= -b_num(n-3)/(n-3)/fac(n)
       ENDIF
      ENDDO 
!
      z(1)= 0.d0
      z(2)= theta
!
      sum2= c2.lcc.z
      sum3= c3.lcc.z
      sum4= c4.lcc.z
!
      pw2= z
      pw3= z.cp.pw2
      pw4= z.cp.pw3
!
      cb2(1:2)= co(1:2)*(psi(2)-psi(1))-ln(1:2)
      cb3(1:2)= co(1:2)*(psi(3)-psi(1))-ln(1:2)
      cb4(1:2)= co(1:2)*(psi(4)-psi(1))-ln(1:2)
!
      cb2= pw2.cp.cb2
      cb3= pw3.cp.cb3
      cb4= pw4.cp.cb4
!
      value(1,1:2)= sum2(1:2)+cb2(1:2)
!
      value(2,1:2)= sum3(1:2)+0.5d0*cb3(1:2)
!
      value(3,1:2)= 1.d0/6.d0*(sum4(1:2)+1.d0/6.d0*cb4(1:2))
!
      RETURN
!
      END FUNCTION VBFNLO_HIGGSTO_poly_unit
! 
!-----------------------------------------------------------------
!
      FUNCTION VBFNLO_HIGGSTO_li2(x) RESULT(value)
      !USE VBFNLO_HIGGSTO_acmplx_pro
      !USE VBFNLO_HIGGSTO_acmplx_rat
      !USE VBFNLO_HIGGSTO_full_ln
      !USE VBFNLO_HIGGSTO_bernoulli
      !USE VBFNLO_HIGGSTO_riemann
      !USE VBFNLO_HIGGSTO_units
!
      IMPLICIT NONE
!
      REAL*8, dimension(2) ::x,value
      INTENT(IN) x
!
      INTEGER n
      REAL*8 ym,ym2,sign1,sign2,sign3,fact, parr,pari,parm,parm2,zr
      REAL*8, dimension(2) :: clnx,clnomx,clnoy,clnz,clnomz, &
     &        add1,add2,add3,par,res
      REAL*8, dimension(0:14) :: bf
      REAL*8, dimension(15) :: ct,sn
      REAL*8, dimension(2) :: omx,y,oy,omy,z,omz,t,omt
!
      omx= co-x
      IF(x(1) < 0.d0) THEN
       y= omx
       sign1= -1.d0
       clnx= x(1).fln.x(2)
       clnomx= omx(1).fln.omx(2)
       add1= pis/6.d0*co-(clnx.cp.clnomx)
      ELSE
       y= x
       sign1= 1.d0
       add1= 0.d0
      ENDIF
      omy= co-y
      ym2= y(1)*y(1)+y(2)*y(2)
      ym= sqrt(ym2)
      IF(ym > 1.d0) THEN
       z(1)= y(1)/ym2
       z(2)= -y(2)/ym2
       sign2= -1.d0
       oy= -y
       clnoy= oy(1).fln.oy(2)
       add2= -pis/6.d0*co-0.5d0*(clnoy.cp.clnoy)
      ELSE
       z= y
       sign2= 1.d0
       add2= 0.d0
      ENDIF
      omz= co-z
      zr= z(1)
      IF(zr > 0.5d0) THEN
       t= co-z
       omt= co-t
       sign3= -1.d0
       clnz= z(1).fln.z(2)
       clnomz= omz(1).fln.omz(2)
       add3= pis/6.d0*co-(clnz.cp.clnomz)
      ELSE
       t= z
       omt= co-t
       sign3= 1.d0
       add3= 0.d0
      ENDIF
      par= omt(1).fln.omt(2)
      fact= 1.d0
      DO n=0,14
       bf(n)= b_num(n)/fact
       fact= fact*(n+2.d0)
      ENDDO
      parr= par(1)
      pari= par(2)
      parm2= parr*parr+pari*pari
      parm= sqrt(parm2)
      ct(1)= parr/parm
      sn(1)= pari/parm
      DO n=2,15
       ct(n)= ct(1)*ct(n-1)-sn(1)*sn(n-1)
       sn(n)= sn(1)*ct(n-1)+ct(1)*sn(n-1)
      ENDDO
!      
      res(1)= -((((((((bf(14)*ct(15)*parm2+bf(12)*ct(13))*parm2+ &
     &                 bf(10)*ct(11))*parm2+bf(8)*ct(9))*parm2+ &
     &                 bf(6)*ct(7))*parm2+bf(4)*ct(5))*parm2+ &
     &                 bf(2)*ct(3))*(-parm)+bf(1)*ct(2))*(-parm)+ &
     &                 bf(0)*ct(1))*parm 
      res(2)= -((((((((bf(14)*sn(15)*parm2+bf(12)*sn(13))*parm2+ &
     &                 bf(10)*sn(11))*parm2+bf(8)*sn(9))*parm2+ &
     &                 bf(6)*sn(7))*parm2+bf(4)*sn(5))*parm2+ &
     &                 bf(2)*sn(3))*(-parm)+bf(1)*sn(2))*(-parm)+ &
     &                 bf(0)*sn(1))*parm 
!
      value= sign1*(sign2*(sign3*res+add3)+add2)+add1
!
      RETURN
!
      END FUNCTION VBFNLO_HIGGSTO_li2
!
!------------------------------------------------------------------
!
      FUNCTION VBFNLO_HIGGSTO_li3(x)
      !USE VBFNLO_HIGGSTO_riemann
      !USE VBFNLO_HIGGSTO_units
!
      IMPLICIT NONE
!
      REAL*8, dimension (:) :: x
      REAL*8 VBFNLO_HIGGSTO_li3(size(x))
      INTENT(IN) x
      INTEGER n
      REAL*8 xm,xm2,xtst,pr,pj,p2,pm,pr1,pj1,pm1,p12,pr2,pj2,p22,pm2, &
     &       rln2_x,iln2_x,tm2,y2r,ym2,ytst
      REAL*8, dimension(0:14) :: bf
      REAL*8, dimension(15) :: ct,sn,ct1,sn1,ct2,sn2
      REAL*8, dimension(2) :: y,addx,ox,clnx,par,res,u1,u2,ln_y,omy, &
     &  ln_omy,addy,par1,par2,res1,res2,t,resa,resb,ln_t,res3,res4, &
     &  ln_omt,addt,addt2,omt,omu1,omu2
      REAL*8 :: b(0:14)=(/1.d0,-0.75d0, &
     &       0.236111111111111111111111111111111d0, &
     &       -3.472222222222222222222222222222222d-2, &
     &       6.481481481481481481481481481481482d-4, &
     &       4.861111111111111111111111111111111d-4, &
     &       -2.393550012597631645250692869740488d-5, &
     &       -1.062925170068027210884353741496599d-5, &
     &       7.794784580498866213151927437641723d-7, &
     &       2.526087595532039976484420928865373d-7, &
     &       -2.359163915200471237027273583310139d-8, &
     &       -6.168132746415574698402981231264060d-9, &
     &       6.824456748981078267312315451125495d-10, &
     &       1.524285616929084572552216019859487d-10, &
     &       -1.916909414174054295837274763110831d-11/)
!
      FORALL (n=0:14) bf(n)= b(n)/(n+1.d0)
!
      xm2= x(1)*x(1)+x(2)*x(2)
      xm= sqrt(xm2)
!
!-----the modulus of x is checked
!
      xtst= xm-1.d0
      IF(xtst <= 1.d-20) THEN
       y= x
       addx= 0.d0        
      ELSE IF(xm > 1.d-20) THEN
       y(1)= x(1)/xm2
       y(2)= -x(2)/xm2
       ox= -x
       clnx= VBFNLO_HIGGSTO_cqlnx(ox)
       rln2_x= clnx(1)*clnx(1)
       iln2_x= clnx(2)*clnx(2)
       addx(1)= -clnx(1)*(rz2+1.d0/6.d0*(rln2_x-3.d0*iln2_x))
       addx(2)= -clnx(2)*(rz2+1.d0/6.d0*(3.d0*rln2_x-iln2_x))
      ENDIF
!
!-----once x --> y, |y|<1 the sign of re(y) is checked
!     if re(y)>0 a transformation is required for re(y)>1/2
!
      y2r= y(1)*y(1)-y(2)*y(2) 
      IF(y(1) >= 0.d0.or.y2r < 0.d0) THEN
       ytst= y(1)-0.5d0
       IF(ytst <= 1.d-20) THEN
!
!-----li_3(y) is computed
!
        omy= co-y
        par= VBFNLO_HIGGSTO_cqlnomx(y,omy)
        pr= -par(1)
        pj= -par(2)
        p2= pr*pr+pj*pj
        pm= sqrt(p2)
        ct(1)= pr/pm
        sn(1)= pj/pm
        DO n=2,15
         ct(n)= ct(1)*ct(n-1)-sn(1)*sn(n-1)
         sn(n)= sn(1)*ct(n-1)+ct(1)*sn(n-1)
        ENDDO
        res(1)= pm*(bf(0)*ct(1)+pm*(bf(1)*ct(2)+pm* &
     &              (bf(2)*ct(3)+pm*(bf(3)*ct(4)+pm* &
     &              (bf(4)*ct(5)+pm*(bf(5)*ct(6)+pm* &
     &              (bf(6)*ct(7)+pm*(bf(7)*ct(8)+pm* &
     &              (bf(8)*ct(9)+pm*(bf(9)*ct(10)+pm* &
     &              (bf(10)*ct(11)+pm*(bf(11)*ct(12)+pm* &
     &              (bf(12)*ct(13)+pm*(bf(13)*ct(14)+pm* &
     &              (bf(14)*ct(15))))))))))))))))
        res(2)= pm*(bf(0)*sn(1)+pm*(bf(1)*sn(2)+pm* &
     &              (bf(2)*sn(3)+pm*(bf(3)*sn(4)+pm* &
     &              (bf(4)*sn(5)+pm*(bf(5)*sn(6)+pm* &
     &              (bf(6)*sn(7)+pm*(bf(7)*sn(8)+pm* &
     &              (bf(8)*sn(9)+pm*(bf(9)*sn(10)+pm* &
     &              (bf(10)*sn(11)+pm*(bf(11)*sn(12)+pm* &
     &              (bf(12)*sn(13)+pm*(bf(13)*sn(14)+pm* &
     &              (bf(14)*sn(15))))))))))))))))
        VBFNLO_HIGGSTO_li3= res+addx
        RETURN
       ELSE IF(ytst > 1.d-20) THEN
        ym2= y(1)*y(1)+y(2)*y(2)
        u1(1)= 1.d0-y(1)/ym2
        u1(2)= y(2)/ym2
        u2= co-y
        ln_y= VBFNLO_HIGGSTO_cqlnx(y)
        omy= co-y
        ln_omy= VBFNLO_HIGGSTO_cqlnomx(y,omy)
        addy(1)= rz3+rz2*ln_y(1)+1.d0/6.d0*ln_y(1)* &
     &               (ln_y(1)*ln_y(1)-3.d0*ln_y(2)*ln_y(2))- &
     &               0.5d0*ln_omy(1)*(ln_y(1)*ln_y(1)-ln_y(2)* &
     &               ln_y(2))+ln_y(1)*ln_y(2)*ln_omy(2)
        addy(2)= rz2*ln_y(2)+1.d0/6.d0*ln_y(2)*(3.d0* &
     &               ln_y(1)*ln_y(1)-ln_y(2)*ln_y(2))-0.5d0* &
     &               ln_omy(2)*(ln_y(1)*ln_y(1)-ln_y(2)*ln_y(2))- &
     &               ln_y(1)*ln_omy(1)*ln_y(2)
!
!-----li_3(1-1/y) is computed
!
        omu1= co-u1
        par1= VBFNLO_HIGGSTO_cqlnomx(u1,omu1)
        pr1= -par1(1)
        pj1= -par1(2)
        p12= pr1*pr1+pj1*pj1
        pm1= sqrt(p12)
        ct1(1)= pr1/pm1
        sn1(1)= pj1/pm1
        DO n=2,15
         ct1(n)= ct1(1)*ct1(n-1)-sn1(1)*sn1(n-1)
         sn1(n)= sn1(1)*ct1(n-1)+ct1(1)*sn1(n-1)
        ENDDO
        res1(1)= pm1*(bf(0)*ct1(1)+pm1*(bf(1)*ct1(2)+pm1* &
     &               (bf(2)*ct1(3)+pm1*(bf(3)*ct1(4)+pm1* &
     &               (bf(4)*ct1(5)+pm1*(bf(5)*ct1(6)+pm1* &
     &               (bf(6)*ct1(7)+pm1*(bf(7)*ct1(8)+pm1* &
     &               (bf(8)*ct1(9)+pm1*(bf(9)*ct1(10)+pm1* &
     &               (bf(10)*ct1(11)+pm1*(bf(11)*ct1(12)+pm1* &
     &               (bf(12)*ct1(13)+pm1*(bf(13)*ct1(14)+pm1* &
     &               (bf(14)*ct1(15))))))))))))))))
        res1(2)= pm1*(bf(0)*sn1(1)+pm1*(bf(1)*sn1(2)+pm1* &
     &               (bf(2)*sn1(3)+pm1*(bf(3)*sn1(4)+pm1* &
     &               (bf(4)*sn1(5)+pm1*(bf(5)*sn1(6)+pm1* &
     &               (bf(6)*sn1(7)+pm1*(bf(7)*sn1(8)+pm1* &
     &               (bf(8)*sn1(9)+pm1*(bf(9)*sn1(10)+pm1* &
     &               (bf(10)*sn1(11)+pm1*(bf(11)*sn1(12)+pm1* &
     &               (bf(12)*sn1(13)+pm1*(bf(13)*sn1(14)+pm1* &
     &               (bf(14)*sn1(15))))))))))))))))
!
!-----li_3(1-y) is computed
!
        omu2= co-u2
        par2= VBFNLO_HIGGSTO_cqlnomx(u2,omu2)
        pr2= -par2(1)
        pj2= -par2(2)
        p22= pr2*pr2+pj2*pj2
        pm2= sqrt(p22)
        ct2(1)= pr2/pm2
        sn2(1)= pj2/pm2
        DO n=2,15
         ct2(n)= ct2(1)*ct2(n-1)-sn2(1)*sn2(n-1)
         sn2(n)= sn2(1)*ct2(n-1)+ct2(1)*sn2(n-1)
        ENDDO
        res2(1)= pm2*(bf(0)*ct2(1)+pm2*(bf(1)*ct2(2)+pm2* &
     &               (bf(2)*ct2(3)+pm2*(bf(3)*ct2(4)+pm2* &
     &               (bf(4)*ct2(5)+pm2*(bf(5)*ct2(6)+pm2* &
     &               (bf(6)*ct2(7)+pm2*(bf(7)*ct2(8)+pm2* &
     &               (bf(8)*ct2(9)+pm2*(bf(9)*ct2(10)+pm2* &
     &               (bf(10)*ct2(11)+pm2*(bf(11)*ct2(12)+pm2* &
     &               (bf(12)*ct2(13)+pm2*(bf(13)*ct2(14)+pm2* &
     &               (bf(14)*ct2(15))))))))))))))))
        res2(2)= pm2*(bf(0)*sn2(1)+pm2*(bf(1)*sn2(2)+pm2* &
     &               (bf(2)*sn2(3)+pm2*(bf(3)*sn2(4)+pm2* &
     &               (bf(4)*sn2(5)+pm2*(bf(5)*sn2(6)+pm2* &
     &               (bf(6)*sn2(7)+pm2*(bf(7)*sn2(8)+pm2* &
     &               (bf(8)*sn2(9)+pm2*(bf(9)*sn2(10)+pm2* &
     &               (bf(10)*sn2(11)+pm2*(bf(11)*sn2(12)+pm2* &
     &               (bf(12)*sn2(13)+pm2*(bf(13)*sn2(14)+pm2* &
     &               (bf(14)*sn2(15))))))))))))))))
        VBFNLO_HIGGSTO_li3= -res1-res2+addx+addy
        RETURN
       ENDIF   
!
!-----if re(y)<0 a transformation is required in terms of t = -y 
!     and of t^2
!
      ELSE IF(y(1) < 0.d0) THEN
!
!-----first t
!
       t= -y
       IF(t(1) <= 0.5d0) THEN
!
!-----li_3(t) is computed
!
        omt= co-t
        par= VBFNLO_HIGGSTO_cqlnomx(t,omt)
        pr= -par(1)
        pj= -par(2)
        p2= pr*pr+pj*pj
        pm= sqrt(p2)
        ct(1)= pr/pm
        sn(1)= pj/pm
        DO n=2,15
         ct(n)= ct(1)*ct(n-1)-sn(1)*sn(n-1)
         sn(n)= sn(1)*ct(n-1)+ct(1)*sn(n-1)
        ENDDO
        resa(1)= pm*(bf(0)*ct(1)+pm*(bf(1)*ct(2)+pm* &
     &              (bf(2)*ct(3)+pm*(bf(3)*ct(4)+pm* &
     &              (bf(4)*ct(5)+pm*(bf(5)*ct(6)+pm* &
     &              (bf(6)*ct(7)+pm*(bf(7)*ct(8)+pm* &
     &              (bf(8)*ct(9)+pm*(bf(9)*ct(10)+pm* &
     &              (bf(10)*ct(11)+pm*(bf(11)*ct(12)+pm* &
     &              (bf(12)*ct(13)+pm*(bf(13)*ct(14)+pm* &
     &              (bf(14)*ct(15))))))))))))))))
        resa(2)= pm*(bf(0)*sn(1)+pm*(bf(1)*sn(2)+pm* &
     &              (bf(2)*sn(3)+pm*(bf(3)*sn(4)+pm* &
     &              (bf(4)*sn(5)+pm*(bf(5)*sn(6)+pm* &
     &              (bf(6)*sn(7)+pm*(bf(7)*sn(8)+pm* &
     &              (bf(8)*sn(9)+pm*(bf(9)*sn(10)+pm* &
     &              (bf(10)*sn(11)+pm*(bf(11)*sn(12)+pm* &
     &              (bf(12)*sn(13)+pm*(bf(13)*sn(14)+pm* &
     &              (bf(14)*sn(15))))))))))))))))
       ELSE IF(t(1) > 0.5d0) THEN
        tm2= t(1)*t(1)+t(2)*t(2)
        u1(1)= 1.d0-t(1)/tm2
        u1(2)= t(2)/tm2
        u2= co-t
        ln_t= VBFNLO_HIGGSTO_cqlnx(t)
        omt= co-t
        ln_omt= VBFNLO_HIGGSTO_cqlnomx(t,omt)
        addt(1)= rz3+rz2*ln_t(1)+1.d0/6.d0*ln_t(1)* &
     &              (ln_t(1)*ln_t(1)-3.d0*ln_t(2)*ln_t(2))- &
     &              0.5d0*ln_omt(1)*(ln_t(1)*ln_t(1)-ln_t(2)* &
     &              ln_t(2))+ln_t(1)*ln_t(2)*ln_omt(2)
        addt(2)= rz2*ln_t(2)+1.d0/6.d0*ln_t(2)*(3.d0* &
     &               ln_t(1)*ln_t(1)-ln_t(2)*ln_t(2))-0.5d0* &
     &               ln_omt(2)*(ln_t(1)*ln_t(1)-ln_t(2)*ln_t(2))- &
     &               ln_t(1)*ln_omt(1)*ln_t(2)
!
!-----li3(1-1/t) is computed
!
        omu1= co-u1
        par1= VBFNLO_HIGGSTO_cqlnomx(u1,omu1)
        pr1= -par1(1)
        pj1= -par1(2)
        p12= pr1*pr1+pj1*pj1
        pm1= sqrt(p12)
        ct1(1)= pr1/pm1
        sn1(1)= pj1/pm1
        DO n=2,15
         ct1(n)= ct1(1)*ct1(n-1)-sn1(1)*sn1(n-1)
         sn1(n)= sn1(1)*ct1(n-1)+ct1(1)*sn1(n-1)
        ENDDO
        res1(1)= pm1*(bf(0)*ct1(1)+pm1*(bf(1)*ct1(2)+pm1* &
     &               (bf(2)*ct1(3)+pm1*(bf(3)*ct1(4)+pm1* &
     &               (bf(4)*ct1(5)+pm1*(bf(5)*ct1(6)+pm1* &
     &               (bf(6)*ct1(7)+pm1*(bf(7)*ct1(8)+pm1* &
     &               (bf(8)*ct1(9)+pm1*(bf(9)*ct1(10)+pm1* &
     &               (bf(10)*ct1(11)+pm1*(bf(11)*ct1(12)+pm1* &
     &               (bf(12)*ct1(13)+pm1*(bf(13)*ct1(14)+pm1* &
     &               (bf(14)*ct1(15))))))))))))))))
        res1(2)= pm1*(bf(0)*sn1(1)+pm1*(bf(1)*sn1(2)+pm1* &
     &               (bf(2)*sn1(3)+pm1*(bf(3)*sn1(4)+pm1* &
     &               (bf(4)*sn1(5)+pm1*(bf(5)*sn1(6)+pm1* &
     &               (bf(6)*sn1(7)+pm1*(bf(7)*sn1(8)+pm1* &
     &               (bf(8)*sn1(9)+pm1*(bf(9)*sn1(10)+pm1* &
     &               (bf(10)*sn1(11)+pm1*(bf(11)*sn1(12)+pm1* &
     &               (bf(12)*sn1(13)+pm1*(bf(13)*sn1(14)+pm1* &
     &               (bf(14)*sn1(15))))))))))))))))
!
!-----li3(1-t) is computed
!
        omu2= co-u2
        par2= VBFNLO_HIGGSTO_cqlnomx(u2,omu2)
        pr2= -par2(1)
        pj2= -par2(2)
        p22= pr2*pr2+pj2*pj2
        pm2= sqrt(p22)
        ct2(1)= pr2/pm2
        sn2(1)= pj2/pm2
        DO n=2,15
         ct2(n)= ct2(1)*ct2(n-1)-sn2(1)*sn2(n-1)
         sn2(n)= sn2(1)*ct2(n-1)+ct2(1)*sn2(n-1)
        ENDDO
        res2(1)= pm2*(bf(0)*ct2(1)+pm2*(bf(1)*ct2(2)+pm2* &
     &               (bf(2)*ct2(3)+pm2*(bf(3)*ct2(4)+pm2* &
     &               (bf(4)*ct2(5)+pm2*(bf(5)*ct2(6)+pm2* &
     &               (bf(6)*ct2(7)+pm2*(bf(7)*ct2(8)+pm2* &
     &               (bf(8)*ct2(9)+pm2*(bf(9)*ct2(10)+pm2* &
     &               (bf(10)*ct2(11)+pm2*(bf(11)*ct2(12)+pm2* &
     &               (bf(12)*ct2(13)+pm2*(bf(13)*ct2(14)+pm2* &
     &               (bf(14)*ct2(15))))))))))))))))
        res2(2)= pm2*(bf(0)*sn2(1)+pm2*(bf(1)*sn2(2)+pm2* &
     &               (bf(2)*sn2(3)+pm2*(bf(3)*sn2(4)+pm2* &
     &               (bf(4)*sn2(5)+pm2*(bf(5)*sn2(6)+pm2* &
     &               (bf(6)*sn2(7)+pm2*(bf(7)*sn2(8)+pm2* &
     &               (bf(8)*sn2(9)+pm2*(bf(9)*sn2(10)+pm2* &
     &               (bf(10)*sn2(11)+pm2*(bf(11)*sn2(12)+pm2* &
     &               (bf(12)*sn2(13)+pm2*(bf(13)*sn2(14)+pm2* &
     &               (bf(14)*sn2(15))))))))))))))))
        resa= -res1-res2+addt
       ENDIF
!
!-----THEN t^2
!
       t(1)= y(1)*y(1)-y(2)*y(2)
       t(2)= 2.d0*y(1)*y(2)
       IF(t(1) <= 0.5d0) THEN
!
!-----li_3(t^2) is computed
!
        omt= co-t
        par= VBFNLO_HIGGSTO_cqlnomx(t,omt)
        pr= -par(1)
        pj= -par(2)
        p2= pr*pr+pj*pj
        pm= sqrt(p2)
        ct(1)= pr/pm
        sn(1)= pj/pm
        DO n=2,15
         ct(n)= ct(1)*ct(n-1)-sn(1)*sn(n-1)
         sn(n)= sn(1)*ct(n-1)+ct(1)*sn(n-1)
        ENDDO
        resb(1)= pm*(bf(0)*ct(1)+pm*(bf(1)*ct(2)+pm* &
     &              (bf(2)*ct(3)+pm*(bf(3)*ct(4)+pm* &
     &              (bf(4)*ct(5)+pm*(bf(5)*ct(6)+pm* &
     &              (bf(6)*ct(7)+pm*(bf(7)*ct(8)+pm* &
     &              (bf(8)*ct(9)+pm*(bf(9)*ct(10)+pm* &
     &              (bf(10)*ct(11)+pm*(bf(11)*ct(12)+pm* &
     &              (bf(12)*ct(13)+pm*(bf(13)*ct(14)+pm* &
     &              (bf(14)*ct(15))))))))))))))))
        resb(2)= pm*(bf(0)*sn(1)+pm*(bf(1)*sn(2)+pm* &
     &              (bf(2)*sn(3)+pm*(bf(3)*sn(4)+pm* &
     &              (bf(4)*sn(5)+pm*(bf(5)*sn(6)+pm* &
     &              (bf(6)*sn(7)+pm*(bf(7)*sn(8)+pm* &
     &              (bf(8)*sn(9)+pm*(bf(9)*sn(10)+pm* &
     &              (bf(10)*sn(11)+pm*(bf(11)*sn(12)+pm* &
     &              (bf(12)*sn(13)+pm*(bf(13)*sn(14)+pm* &
     &              (bf(14)*sn(15))))))))))))))))
       ELSE IF(t(1) > 0.5d0) THEN
        tm2= t(1)*t(1)+t(2)*t(2)
        u1(1)= 1.d0-t(1)/tm2
        u1(2)= t(2)/tm2
        u2= co-t
        ln_t= VBFNLO_HIGGSTO_cqlnx(t)
        omt= co-t
        ln_omt= VBFNLO_HIGGSTO_cqlnomx(t,omt)
        addt2(1)= rz3+rz2*ln_t(1)+1.d0/6.d0*ln_t(1)* &
     &                (ln_t(1)*ln_t(1)-3.d0*ln_t(2)*ln_t(2))- &
     &                0.5d0*ln_omt(1)*(ln_t(1)*ln_t(1)-ln_t(2)* &
     &                ln_t(2))+ln_t(1)*ln_t(2)*ln_omt(2)
        addt2(2)= rz2*ln_t(2)+1.d0/6.d0*ln_t(2)*(3.d0* &
     &                ln_t(1)*ln_t(1)-ln_t(2)*ln_t(2))-0.5d0* &
     &                ln_omt(2)*(ln_t(1)*ln_t(1)-ln_t(2)*ln_t(2))- &
     &                ln_t(1)*ln_omt(1)*ln_t(2)
!
!-----li_3(1-1/t^2) is computed
!
        omu1= co-u1
        par1= VBFNLO_HIGGSTO_cqlnomx(u1,omu1)
        pr1= -par1(1)
        pj1= -par1(2)
        p12= pr1*pr1+pj1*pj1
        pm1= sqrt(p12)
        ct1(1)= pr1/pm1
        sn1(1)= pj1/pm1
        DO n=2,15
         ct1(n)= ct1(1)*ct1(n-1)-sn1(1)*sn1(n-1)
         sn1(n)= sn1(1)*ct1(n-1)+ct1(1)*sn1(n-1)
        ENDDO
        res3(1)= pm1*(bf(0)*ct1(1)+pm1*(bf(1)*ct1(2)+pm1* &
     &               (bf(2)*ct1(3)+pm1*(bf(3)*ct1(4)+pm1* &
     &               (bf(4)*ct1(5)+pm1*(bf(5)*ct1(6)+pm1* &
     &               (bf(6)*ct1(7)+pm1*(bf(7)*ct1(8)+pm1* &
     &               (bf(8)*ct1(9)+pm1*(bf(9)*ct1(10)+pm1* &
     &               (bf(10)*ct1(11)+pm1*(bf(11)*ct1(12)+pm1* &
     &               (bf(12)*ct1(13)+pm1*(bf(13)*ct1(14)+pm1* &
     &               (bf(14)*ct1(15))))))))))))))))
        res3(2)= pm1*(bf(0)*sn1(1)+pm1*(bf(1)*sn1(2)+pm1* &
     &               (bf(2)*sn1(3)+pm1*(bf(3)*sn1(4)+pm1* &
     &               (bf(4)*sn1(5)+pm1*(bf(5)*sn1(6)+pm1* &
     &               (bf(6)*sn1(7)+pm1*(bf(7)*sn1(8)+pm1* &
     &               (bf(8)*sn1(9)+pm1*(bf(9)*sn1(10)+pm1* &
     &               (bf(10)*sn1(11)+pm1*(bf(11)*sn1(12)+pm1* &
     &               (bf(12)*sn1(13)+pm1*(bf(13)*sn1(14)+pm1* &
     &               (bf(14)*sn1(15))))))))))))))))
!
!-----li_3(1-t^2) is computed
!
        omu2= co-u2
        par2= VBFNLO_HIGGSTO_cqlnomx(u2,omu2)
        pr2= -par2(1)
        pj2= -par2(2)
        p22= pr2*pr2+pj2*pj2
        pm2= sqrt(p22)
        ct2(1)= pr2/pm2
        sn2(1)= pj2/pm2
        DO n=2,15
         ct2(n)= ct2(1)*ct2(n-1)-sn2(1)*sn2(n-1)
         sn2(n)= sn2(1)*ct2(n-1)+ct2(1)*sn2(n-1)
        ENDDO
        res4(1)= pm2*(bf(0)*ct2(1)+pm2*(bf(1)*ct2(2)+pm2* &
     &               (bf(2)*ct2(3)+pm2*(bf(3)*ct2(4)+pm2* &
     &               (bf(4)*ct2(5)+pm2*(bf(5)*ct2(6)+pm2* &
     &               (bf(6)*ct2(7)+pm2*(bf(7)*ct2(8)+pm2* &
     &               (bf(8)*ct2(9)+pm2*(bf(9)*ct2(10)+pm2* &
     &               (bf(10)*ct2(11)+pm2*(bf(11)*ct2(12)+pm2* &
     &               (bf(12)*ct2(13)+pm2*(bf(13)*ct2(14)+pm2* &
     &               (bf(14)*ct2(15))))))))))))))))
        res4(2)= pm2*(bf(0)*sn2(1)+pm2*(bf(1)*sn2(2)+pm2* &
     &               (bf(2)*sn2(3)+pm2*(bf(3)*sn2(4)+pm2* &
     &               (bf(4)*sn2(5)+pm2*(bf(5)*sn2(6)+pm2* &
     &               (bf(6)*sn2(7)+pm2*(bf(7)*sn2(8)+pm2* &
     &               (bf(8)*sn2(9)+pm2*(bf(9)*sn2(10)+pm2* &
     &               (bf(10)*sn2(11)+pm2*(bf(11)*sn2(12)+pm2* &
     &               (bf(12)*sn2(13)+pm2*(bf(13)*sn2(14)+pm2* &
     &               (bf(14)*sn2(15))))))))))))))))
        resb= -res3-res4+addt2
       ENDIF
       VBFNLO_HIGGSTO_li3= -resa+0.25d0*resb+addx
       RETURN
      ENDIF
!
      END FUNCTION VBFNLO_HIGGSTO_li3
!
!-----------------------------------------------------------------------------------------
!
      SUBROUTINE VBFNLO_HIGGSTO_init_niels
      !USE VBFNLO_HIGGSTO_common_niels
!
      plr= 0.d0
!
      plr(1,1)= 1.d0
      plr(1,2)= -2.5d-1
      plr(1,3)= 2.77777777777777777d-2
      plr(1,5)= -2.77777777777777777d-4
      plr(1,7)= 4.72411186696900982d-6
      plr(1,9)= -9.18577307466196355d-8
      plr(1,11)= 1.89788699889709990d-9
      plr(1,13)= -4.06476164514422552d-11
      plr(1,15)= 8.92169102045645255d-13
!
      plr(2,1)= 1.d0
      plr(2,2)= -3.75d-1
      plr(2,3)= 7.87037037037037037d-2
      plr(2,4)= -8.68055555555555555d-3
      plr(2,5)= 1.29629629629629629d-4
      plr(2,6)= 8.10185185185185185d-5
      plr(2,7)= -3.41935716085375949d-6
      plr(2,8)= -1.32865646258503401d-6
      plr(2,9)= 8.66087175610985134d-8
      plr(2,10)= 2.52608759553203997d-8
      plr(2,11)= -2.14469446836406476d-9
      plr(2,12)= -5.14011062201297891d-10
      plr(2,13)= 5.24958211460082943d-11
      plr(2,14)= 1.08877544066363183d-11
      plr(2,15)= -1.27793960944936953d-12
!
      plr(3,2)= 2.5d-1
      plr(3,3)= -8.33333333333333333d-2
      plr(3,4)= 1.04166666666666666d-2
      plr(3,6)= -1.15740740740740740d-4
      plr(3,8)= 2.06679894179894179d-6
      plr(3,10)= -4.13359788359788359d-8
      plr(3,12)= 8.69864874494504124d-10
      plr(3,14)= -1.88721076381696185d-11
!
      plr_4= 0.d0
!
      plr_4(1,1)= 1.d0
      plr_4(1,2)= -4.375d-1
      plr_4(1,3)= 1.16512345679012345d-1
      plr_4(1,4)= -1.98206018518518518d-2
      plr_4(1,5)= 1.92793209876543209d-3
      plr_4(1,6)= -3.10570987654320987d-5
      plr_4(1,7)= -1.56240091148578352d-5
      plr_4(1,8)= 8.48512354677320663d-7
      plr_4(1,9)= 2.29096166031897114d-7
      plr_4(1,10)= -2.18326142185269169d-8
      plr_4(1,11)= -3.88282487917201557d-9
      plr_4(1,12)= 5.44629210322033211d-10
      plr_4(1,13)= 6.96080521068272540d-11
      plr_4(1,14)= -1.33757376864452151d-11
      plr_4(1,15)= -1.27848526852665716d-12
!
      plr_4(2,2)= 2.5d-1
      plr_4(2,3)= -1.25d-1
      plr_4(2,4)= 2.95138888888888888d-2
      plr_4(2,5)= -3.47222222222222222d-3
      plr_4(2,6)= 5.40123456790123456d-5
      plr_4(2,7)= 3.47222222222222222d-5
      plr_4(2,8)= -1.49596875787351977d-6
      plr_4(2,9)= -5.90513983371126228d-7
      plr_4(2,10)= 3.89739229024943310d-8
      plr_4(2,11)= 1.14822163433274544d-8
      plr_4(2,12)= -9.82984964666863015d-10
      plr_4(2,13)= -2.37235874862137488d-10
      plr_4(2,14)= 2.43730598177895652d-11
      plr_4(2,15)= 5.08095205643028190d-12
!
      plr_4(3,2)= 5.55555555555555555d-2
      plr_4(3,3)= -2.08333333333333333d-2
      plr_4(3,4)= 2.77777777777777777d-3
      plr_4(3,6)= -3.30687830687830687d-5
      plr_4(3,8)= 6.12384871644130903d-7
      plr_4(3,10)= -1.25260541927208593d-8
      plr_4(3,12)= 2.67650730613693576d-10
      plr_4(3,14)= -5.87132237631943687d-12
!
      RETURN
!
      END SUBROUTINE VBFNLO_HIGGSTO_init_niels
!
!-----CQLNX---------------------------------------------
!
!--- Computes  ln(z)                              
!
      FUNCTION VBFNLO_HIGGSTO_cqlnx(arg) RESULT(res)
      !USE VBFNLO_HIGGSTO_riemann
      !USE VBFNLO_HIGGSTO_units
      IMPLICIT NONE
!
      REAL*8 teta,zm,zm2,tnteta,sr,si
      REAL*8, dimension(2) :: arg,aarg,res
      INTENT(IN) arg
!
      IF((abs(arg(2))-eps).eq.0.d0) THEN
       res(1)= log(abs(arg(1)))
       IF(arg(1) > 0.d0) THEN
        res(2)= 0.d0
       ELSE
        res(2)= pi*sign(one,arg(2))
       ENDIF
       RETURN
      ENDIF
!
      aarg= abs(arg)
      zm2= (arg(1))**2+(arg(2))**2
      zm= sqrt(zm2)
      res(1)= log(zm)
      IF(arg(1).eq.0.d0) THEN
       IF(arg(2) > 0.d0) THEN
        teta= pi/2.d0
       ELSE
        teta= -pi/2.d0
       ENDIF
       res(2)= teta
       RETURN
      ELSE IF(arg(2).eq.0.d0) THEN 
       IF(arg(1) > 0.d0) THEN
        teta= 0.d0
       ELSE
        teta= pi
       ENDIF
       res(2)= teta
       RETURN
      ELSE
       tnteta= aarg(2)/aarg(1)
       teta= atan(tnteta)
       sr= arg(1)/aarg(1)
       si= arg(2)/aarg(2)
       IF(sr > 0.d0) THEN
        res(2)= si*teta
       ELSE
        res(2)= si*(pi-teta)
       ENDIF
       RETURN
      ENDIF
!
      END FUNCTION VBFNLO_HIGGSTO_cqlnx
!
!-----VBFNLO_HIGGSTO_cqlnomx---------------------------------------
!
!--- Computes ln(1-x), usually |x| << 1                 
!
      FUNCTION VBFNLO_HIGGSTO_cqlnomx(arg,omarg) RESULT(res)
!
      IMPLICIT NONE
!
      INTEGER n,k
      REAL*8 zr,zi,zm2,zm
      REAL*8, dimension(2) :: arg,omarg,res,ares
      REAL*8, dimension(10) :: ct,sn
      INTENT(IN) arg,omarg
!
      zr= arg(1)
      zi= arg(2)
      zm2= zr*zr+zi*zi
      zm= sqrt(zm2)
      IF(zm < 1.d-7) THEN
       ct(1)= zr/zm
       sn(1)= zi/zm
       DO n=2,10
        ct(n)= ct(1)*ct(n-1)-sn(1)*sn(n-1)
        sn(n)= sn(1)*ct(n-1)+ct(1)*sn(n-1)
       ENDDO
       ares(1)= ct(10)/10.d0
       ares(2)= sn(10)/10.d0
       DO k=9,1,-1
        ares(1)= ares(1)*zm+ct(k)/k
        ares(2)= ares(2)*zm+sn(k)/k
       ENDDO
       ares(1)= -ares(1)*zm
       ares(2)= -ares(2)*zm
      ELSE
       ares= VBFNLO_HIGGSTO_cqlnx(omarg)
      ENDIF
      res= ares
!
      RETURN
!
      END FUNCTION VBFNLO_HIGGSTO_cqlnomx
!
      !END MODULE VBFNLO_HIGGSTO_sp_fun
!
!
      SUBROUTINE VBFNLO_HIGGSTO_INITALPHAS(IORD, FR2, MUR, ASMUR, MC, MB, MT)
      !USE VBFNLO_HIGGSTO_DZpar 
!

!--   IORD = 0 (LO), 1 (NLO), 2 (NNLO), 3 (NNNLO).
!--   FR2 = ratio of mu_f^2 to mu_r^2 (must be a fixed value).
!--   MUR = input renormalisation scale (in GeV) for alpha_s.
!--   ASMUR = input value of alpha_s at the renormalisation scale MUR.
!--   MC,MB,MT = heavy quark masses in GeV.
!
      IMPLICIT NONE
!
      INTEGER IORD
      REAL*8 FR2,MUR,ASMUR,MC,MB,MT,A,B, & ! ,VBFNLO_HIGGSTO_DZERO, &
     &       R0,ASI
     ! &       VBFNLO_HIGGSTO_FINDALPHASR0,R0,ASI
      REAL*8, parameter :: EPS=1.d-10
      INTEGER, parameter :: MAXF=10000
      INTEGER, parameter :: MODE=1
      ! EXTERNAL VBFNLO_HIGGSTO_FINDALPHASR0

      IF(MUR*sqrt(FR2).LE.MC) THEN ! Check that MUF <= MC.
       R0 = MUR
       ASI = ASMUR
      ELSE                      ! Solve for alpha_s at R0 = 1 GeV.
!
!--   Copy variables to common block.
!
       R0c = 1.d0/sqrt(FR2)
       IORDc = IORD
       FR2c = FR2
       MURc = MUR
       ASMURc = ASMUR
       MCc = MC
       MBc = MB
       MTc = MT
!
!--   Now get alpha_s(R0) corresponding to alpha_s(MUR).
!
       A = 0.02d0              ! lower bound for alpha_s(R0)
       B = 2.00d0              ! upper bound for alpha_s(R0)
       R0 = R0c
       ASI = VBFNLO_HIGGSTO_DZERO(A,B,EPS,MAXF,VBFNLO_HIGGSTO_FINDALPHASR0,MODE)
      ENDIF

      CALL VBFNLO_HIGGSTO_INITALPHASR0(IORD,FR2,R0,ASI,MC,MB,MT)

      RETURN
      END SUBROUTINE VBFNLO_HIGGSTO_INITALPHAS
!
!----------------------------------------------------------------------
!
!--   Find the zero of this function using VBFNLO_HIGGSTO_DZEROX.
!
      FUNCTION VBFNLO_HIGGSTO_FINDALPHASR0(ASI)
      !USE VBFNLO_HIGGSTO_DZpar
!
      IMPLICIT NONE
!
      REAL*8 VBFNLO_HIGGSTO_FINDALPHASR0,ASI
!
      CALL VBFNLO_HIGGSTO_INITALPHASR0(IORDc,FR2c,R0c,ASI,MCc,MBc,MTc)
!
      VBFNLO_HIGGSTO_FINDALPHASR0 = VBFNLO_HIGGSTO_ALPHAS(MURc) - ASMURc ! solve equal to zero

 
      RETURN
      END FUNCTION VBFNLO_HIGGSTO_FINDALPHASR0
!
!----------------------------------------------------------------------
!
      SUBROUTINE VBFNLO_HIGGSTO_INITALPHASR0(IORD,FR2,R0,ASI,MC,MB,MT)
      !USE VBFNLO_HIGGSTO_RZETA  
      !USE VBFNLO_HIGGSTO_COLOUR 
      !USE VBFNLO_HIGGSTO_ASINP  
      !USE VBFNLO_HIGGSTO_ASPAR  
      !USE VBFNLO_HIGGSTO_VARFLV 
      !USE VBFNLO_HIGGSTO_NFFIX  
      !USE VBFNLO_HIGGSTO_FRRAT  
!
!--   IORD = 0 (LO), 1 (NLO), 2 (NNLO), 3 (NNNLO).
!--   FR2 = ratio of mu_f^2 to mu_r^2 (must be a fixed value).
!--   R0 = input renormalisation scale (in GeV) for alphas_s.
!--   ASI = input value of alpha_s at the renormalisation scale R0.
!--   MC,MB,MT = heavy quark masses in GeV.
!--   Must have R0*sqrt(FR2) <= MC to call this subroutine.
!
      IMPLICIT NONE
!
      INTEGER IORD
      REAL*8 FR2,R0,ASI,MC,MB,MT,R20,MC2,MB2,MT2
      REAL*8, parameter :: PI= 3.14159265358979d0
!
! ..QCD colour factors
!
      CA = 3.d0
      CF = 4.d0/3.d0
      TR = 0.5d0
!
! ..The lowest integer values of the Zeta function
!
      ZETA(1) = 0.57721566490153d0
      ZETA(2) = 1.644934066848226d0
      ZETA(3) = 1.202056903159594d0
      ZETA(4) = 1.082323233711138d0
      ZETA(5) = 1.036927755143370d0
      ZETA(6) = 1.017343061984449d0

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
      IF( (IVFNS .EQ. 0) .AND. (NFF .LT. 3) ) THEN
       print*, 'Wrong flavour number for FFNS evolution. STOP'
       STOP
      ENDIF
      IF( (IVFNS .EQ. 0) .AND. (NFF .GT. 5) ) THEN
       print*, 'Wrong flavour number for FFNS evolution. STOP'
       STOP
      ENDIF
!     
      IF( NAORD .GT. 3 ) THEN
       print*, 'Specified order in a_s too high. STOP' 
       STOP
      ENDIF
!
      IF( (IVFNS .NE. 0) .AND. (FR2 .GT. 4.001d0) ) THEN
       print*, 'Too low mu_r for VFNS evolution. STOP'
       STOP
      ENDIF
!
      IF( (IVFNS .EQ. 1) .AND. (M20 .GT. MC2) ) THEN
       print*, 'Too high mu_0 for VFNS evolution. STOP'
       STOP
      ENDIF
!     
      IF( (ASI .GT. 2.d0) .OR. (ASI .LT. 2.D-2) ) THEN
       print*, 'alpha_s out of range. STOP'
       STOP
      ENDIF
!     
      IF( (IVFNS .EQ. 1) .AND. (MC2 .GT. MB2) ) THEN
       print*, 'Wrong charm-bottom mass hierarchy. STOP'
       STOP
      ENDIF
      IF( (IVFNS .EQ. 1) .AND. (MB2 .GT. MT2) ) THEN
       print*, 'Wrong bottom-top mass hierarchy. STOP'
       STOP
      ENDIF
!
!
!--   Store the beta function coefficients in a COMMON block.
!
      CALL VBFNLO_HIGGSTO_BETAFCT
!
!--   Store a_s = alpha_s(mu_r^2)/(4 pi) at the input scale R0.
!
      AS0 = ASI / (4.d0* PI)
!
!--   Store alpha_s at the heavy flavour thresholds in a COMMON block.
!
      IF(IVFNS .NE. 0) THEN
       CALL VBFNLO_HIGGSTO_EVNFTHR(MC2,MB2,MT2)
      ENDIF

      RETURN
      END SUBROUTINE VBFNLO_HIGGSTO_INITALPHASR0
!
!----------------------------------------------------------------------
!
      FUNCTION VBFNLO_HIGGSTO_ALPHAS(MUR)
      !USE VBFNLO_HIGGSTO_NFFIX  
      !USE VBFNLO_HIGGSTO_VARFLV 
      !USE VBFNLO_HIGGSTO_FRRAT  
      !USE VBFNLO_HIGGSTO_ASINP  
      !USE VBFNLO_HIGGSTO_ASFTHR 
!
      IMPLICIT NONE
      INTEGER NF
!      REAL*8 M2,MUR,R2,ASI,ASF,R20,R2T,R2B,R2C,AS,VBFNLO_HIGGSTO_ALPHAS
      REAL*8 M2,MUR,R2,ASI,ASF,R20,R2T,R2B,R2C,VBFNLO_HIGGSTO_ALPHAS
      REAL*8, parameter :: PI= 3.14159265358979d0
!
      ! INTERFACE
      !  FUNCTION VBFNLO_HIGGSTO_AS(R2,R20,AS0,NF)
      !  !USE VBFNLO_HIGGSTO_ASPAR  
      !  !USE VBFNLO_HIGGSTO_BETACOM   
      !  IMPLICIT NONE
      !  INTEGER NF
      !  REAL*8 R2,R20,AS0,VBFNLO_HIGGSTO_AS
      !  END FUNCTION VBFNLO_HIGGSTO_AS
      ! END INTERFACE 
!
! ..Input common blocks 
! 

      R2 = MUR**2
      M2 = R2 * EXP(+LOGFR)
      IF(IVFNS .EQ. 0) THEN
!
!   Fixed number of flavours
!
       NF  = NFF
       R20 = M20 * R2/M2
       ASI = AS0
       ASF = VBFNLO_HIGGSTO_AS(R2,R20,AS0,NF)
!
      ELSE
!
! ..Variable number of flavours
!
       IF(M2 .GT. M2T) THEN
        NF = 6
        R2T = M2T * R2/M2
        ASI = AST
        ASF = VBFNLO_HIGGSTO_AS(R2,R2T,AST,NF)
!
       ELSE IF(M2 .GT. M2B) THEN
        NF = 5
        R2B = M2B * R2/M2
        ASI = ASB
        ASF = VBFNLO_HIGGSTO_AS(R2,R2B,ASB,NF)
!     
       ELSE IF(M2 .GT. M2C) THEN
        NF = 4
        R2C = M2C * R2/M2
        ASI = ASC
        ASF = VBFNLO_HIGGSTO_AS(R2,R2C,ASC,NF)
!     
       ELSE
        NF = 3
        R20 = M20 * R2/M2
        ASI = AS0
        ASF = VBFNLO_HIGGSTO_AS(R2,R20,AS0,NF)
!       
       ENDIF
!
      ENDIF
!
! ..Final value of alpha_s
!
      VBFNLO_HIGGSTO_ALPHAS = 4.d0*PI*ASF
!
      RETURN
      END FUNCTION VBFNLO_HIGGSTO_ALPHAS
!
! =================================================================av==


! =====================================================================
!
! ..The threshold matching of the QCD coupling in the MS(bar) scheme,  
!    a_s = alpha_s(mu_r^2)/(4 pi),  for  NF -> NF + 1  active flavours 
!    up to order a_s^4 (NNNLO).
!
! ..The value  ASNF  of a_s for NF flavours at the matching scale, the 
!    logarithm  LOGRH = ln (mu_r^2/m_H^2) -- where m_H is the pole mass
!    of the heavy quark -- and  NF  are passed as arguments to the 
!    function  VBFNLO_HIGGSTO_ASNF1.  The order of the expansion  NAORD  (defined as 
!    the 'n' in N^nLO) is provided by the common-block  ASPAR.
!
! ..The matching coefficients are inverted from Chetyrkin, Kniehl and
!    Steinhauser, Phys. Rev. Lett. 79 (1997) 2184. The QCD colour
!    factors have been hard-wired in these results. The lowest integer 
!    values of the Zeta function are given by the common-block  RZETA.
!
! =====================================================================
!
!
      FUNCTION VBFNLO_HIGGSTO_ASNF1(ASNF,LOGRH,NF)
      !USE VBFNLO_HIGGSTO_ASPAR 
      !USE VBFNLO_HIGGSTO_RZETA 
!
      IMPLICIT NONE
      INTEGER NF,PRVCLL,K1,K2
      REAL*8 ASNF,LOGRH,CMCI30,CMCF30,CMCF31, &
     &       CMCI31,ASP,LRHP,VBFNLO_HIGGSTO_ASNF1
      REAL*8, dimension(3,0:3) :: CMC
!*
! ---------------------------------------------------------------------
!
! ..Input common-blocks 
!
! ..Variables to be saved for the next call
!
      SAVE CMC,CMCI30,CMCF30,CMCF31,CMCI31,PRVCLL
!
! ---------------------------------------------------------------------
!
! ..The coupling-constant matching coefficients (CMC's) up to NNNLO 
!   (calculated and saved in the first call of this routine)
!
      IF(PRVCLL .NE. 1) THEN
!
       CMC(1,0) =  0.d0
       CMC(1,1) =  2.d0/3.d0
!
       CMC(2,0) = 14.d0/3.d0
       CMC(2,1) = 38.d0/3.d0
       CMC(2,2) =  4.d0/9.d0  
!
       CMCI30 = + 80507.d0/432.d0 * ZETA(3) + 58933.d0/1944.d0  &
     &          + 128.d0/3.d0 * ZETA(2) * (1.+ DLOG(2.d0)/3.d0)
       CMCF30 = - 64.d0/9.d0 * (ZETA(2) + 2479.d0/3456.d0)
       CMCI31 =   8941.d0/27.d0
       CMCF31 = - 409.d0/27.d0
       CMC(3,2) = 511.d0/9.d0
       CMC(3,3) = 8.d0/27.d0
!
       PRVCLL = 1
!
      ENDIF
!
! ---------------------------------------------------------------------
!
! ..The N_f dependent CMC's, and the alpha_s matching at order NAORD 
!
      CMC(3,0) = CMCI30 + NF * CMCF30
      CMC(3,1) = CMCI31 + NF * CMCF31
!
      VBFNLO_HIGGSTO_ASNF1 = ASNF
      IF(NAORD .EQ. 0) GO TO 1
       ASP   = ASNF
!
       DO K1 = 1,NAORD 
        ASP = ASP * ASNF
        LRHP = 1.d0
        DO K2 = 0,K1
         VBFNLO_HIGGSTO_ASNF1 = VBFNLO_HIGGSTO_ASNF1 + ASP * CMC(K1,K2) * LRHP
         LRHP = LRHP * LOGRH
        ENDDO
       ENDDO
!
   1  RETURN
!
      END FUNCTION VBFNLO_HIGGSTO_ASNF1
!
! =================================================================av==
!
! ..The subroutine  EVNFTHR  for the evolution of  a_s = alpha_s/(4 pi)
!    from a three-flavour initial scale to the four- to six-flavour
!    thresholds (identified with the squares of the corresponding quark
!    masses).  The results are written to the common-block  ASFTHR.
!
! ..The input scale  M20 = mu_(f,0)^2  and the corresponding value 
!    AS0  of a_s  are provided by  ASINP.  The fixed scale logarithm
!    LOGFR = ln (mu_f^2/mu_r^2) is specified in  FRRAT.  The alpha_s
!    matching is done by the function VBFNLO_HIGGSTO_ASNF1.
!
! =====================================================================
!
!
       SUBROUTINE VBFNLO_HIGGSTO_EVNFTHR(MC2,MB2,MT2)
       !USE VBFNLO_HIGGSTO_ASINP  
       !USE VBFNLO_HIGGSTO_FRRAT  
       !USE VBFNLO_HIGGSTO_ASFTHR 
!
       IMPLICIT NONE
       REAL*8 MC2,MB2,MT2,R20,R2C,R2B,R2T, & !VBFNLO_HIGGSTO_AS,VBFNLO_HIGGSTO_ASNF1, &
     &        ASC3,ASB4,AST5,SC,SB,ST
!
! ---------------------------------------------------------------------
! 
! ..Input common blocks
!  
!
! ..Output common blocks
!

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
       ASC3 = VBFNLO_HIGGSTO_AS(R2C,R20,AS0,3)
       SC   = log(AS0 / ASC3)
       ASC  = VBFNLO_HIGGSTO_ASNF1(ASC3,-LOGFR,3)
!
! ..Bottom 
!
       M2B  = MB2
       R2B  = M2B * R20/M20
       ASB4 = VBFNLO_HIGGSTO_AS(R2B,R2C,ASC,4)
       SB   = log(ASC / ASB4)
       ASB  = VBFNLO_HIGGSTO_ASNF1(ASB4,-LOGFR,4)
!
! ..Top
!
       M2T  = MT2
       R2T  = M2T * R20/M20
       AST5 = VBFNLO_HIGGSTO_AS(R2T,R2B,ASB,5)
       ST   = log(ASB / AST5)
       AST  = VBFNLO_HIGGSTO_ASNF1(AST5,-LOGFR,5)

       RETURN
       END SUBROUTINE VBFNLO_HIGGSTO_EVNFTHR 
!
! =================================================================av==
!
! ..The running coupling of QCD,  
!
!         AS  =  a_s  =  alpha_s(mu_r^2)/(4 pi),
!
!    obtained by integrating the evolution equation for a fixed number
!    of massless flavours  NF.  Except at leading order (LO),  AS  is 
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
!
      FUNCTION VBFNLO_HIGGSTO_AS(R2, R20,AS0,NF)
      !USE VBFNLO_HIGGSTO_ASPAR  
      !USE VBFNLO_HIGGSTO_BETACOM   
!
      IMPLICIT NONE
      INTEGER NF,K1
      REAL*8 R2,R20,AS0,AS1,AS2,AS3,AS4,FBETA1,FBETA2,FBETA3,FBETA4, &
     &       LRRAT,DLR,XK0,XK1,XK2,XK3,VBFNLO_HIGGSTO_AS
      INTEGER, parameter :: NFMIN =3
      INTEGER, parameter :: NFMAX =6
      REAL*8, parameter :: SXTH =0.166666666666666d0
!
! ---------------------------------------------------------------------
!
! ..Input common-blocks 
!
!
! ..The beta functions FBETAn at N^nLO for n = 1, 2, and 3
!
!
! ---------------------------------------------------------------------
!
! ..Initial value, evolution distance and step size
!
      VBFNLO_HIGGSTO_AS = AS0
      LRRAT = log(R2/R20)
      DLR = LRRAT / NASTPS
!
! ..Solution of the evolution equation depending on  NAORD
!   (fourth-order Runge-Kutta beyond the leading order)
!
       IF(NAORD .EQ. 0) THEN
!
         VBFNLO_HIGGSTO_AS = AS0 / (1.+ BETA0(NF) * AS0 * LRRAT)
!
       ELSE IF(NAORD .EQ. 1) THEN
!
        DO K1 = 1,NASTPS
         AS1= VBFNLO_HIGGSTO_AS
         FBETA1= - AS1**2 * ( BETA0(NF) + AS1 *   BETA1(NF) )
         XK0 = DLR * FBETA1 
         AS2= VBFNLO_HIGGSTO_AS + 0.5d0 * XK0     
         FBETA2= - AS2**2 * ( BETA0(NF) + AS2 *   BETA1(NF) )
         XK1 = DLR * FBETA2
         AS3= VBFNLO_HIGGSTO_AS + 0.5d0 * XK1
         FBETA3= - AS3**2 * ( BETA0(NF) + AS3 *   BETA1(NF) )
         XK2 = DLR * FBETA3
         AS4= VBFNLO_HIGGSTO_AS + XK2
         FBETA4= - AS4**2 * ( BETA0(NF) + AS4 *   BETA1(NF) )
         XK3 = DLR * FBETA4
         VBFNLO_HIGGSTO_AS = VBFNLO_HIGGSTO_AS + SXTH * (XK0 + 2.d0* XK1 + 2.d0* XK2 + XK3)
        ENDDO
!
       ELSE IF(NAORD .EQ. 2) THEN
!
       DO K1 = 1,NASTPS
         AS1= VBFNLO_HIGGSTO_AS
         FBETA1= - AS1**2 * ( BETA0(NF) + AS1 * ( BETA1(NF) &
     &           + AS1 * BETA2(NF) ) )
         XK0 = DLR * FBETA1 
         AS2= VBFNLO_HIGGSTO_AS + 0.5d0 * XK0     
         FBETA2= - AS2**2 * ( BETA0(NF) + AS2 * ( BETA1(NF) &
     &           + AS2 * BETA2(NF) ) )
         XK1 = DLR * FBETA2
         AS3= VBFNLO_HIGGSTO_AS + 0.5d0 * XK1
         FBETA3= - AS3**2 * ( BETA0(NF) + AS3 * ( BETA1(NF) &
     &           + AS3 * BETA2(NF) ) )
         XK2 = DLR * FBETA3
         AS4= VBFNLO_HIGGSTO_AS + XK2
         FBETA4= - AS4**2 * ( BETA0(NF) + AS4 * ( BETA1(NF) &
     &           + AS4 * BETA2(NF) ) )
         XK3 = DLR * FBETA4
         VBFNLO_HIGGSTO_AS = VBFNLO_HIGGSTO_AS + SXTH * (XK0 + 2.d0* XK1 + 2.d0* XK2 + XK3)
        ENDDO
!  
       ELSE IF(NAORD .EQ. 3) THEN
!
        DO K1 = 1,NASTPS
         AS1= VBFNLO_HIGGSTO_AS
         FBETA1= - AS1**2 * ( BETA0(NF) + AS1 * ( BETA1(NF) &
     &           + AS1 * (BETA2(NF) + AS1 * BETA3(NF)) ) )
         XK0 = DLR * FBETA1 
         AS2= VBFNLO_HIGGSTO_AS + 0.5d0 * XK0     
         FBETA2= - AS2**2 * ( BETA0(NF) + AS2 * ( BETA1(NF) &
     &           + AS2 * (BETA2(NF) + AS2 * BETA3(NF)) ) )
         XK1 = DLR * FBETA2
         AS3= VBFNLO_HIGGSTO_AS + 0.5d0 * XK1
         FBETA3= - AS3**2 * ( BETA0(NF) + AS3 * ( BETA1(NF) &
     &           + AS3 * (BETA2(NF) + AS3 * BETA3(NF)) ) )
         XK2 = DLR * FBETA3
         AS4= VBFNLO_HIGGSTO_AS + XK2
         FBETA4= - AS4**2 * ( BETA0(NF) + AS4 * ( BETA1(NF) &
     &           + AS4 * (BETA2(NF) + AS4 * BETA3(NF)) ) )
         XK3 = DLR * FBETA4
         VBFNLO_HIGGSTO_AS = VBFNLO_HIGGSTO_AS + SXTH * (XK0 + 2.d0* XK1 + 2.d0* XK2 + XK3)
        ENDDO
!
       ENDIF
!
! ---------------------------------------------------------------------
!
       RETURN
       END FUNCTION VBFNLO_HIGGSTO_AS
!
! =================================================================av==
!
! ..The subroutine BETAFCT for the coefficients  BETA0...BETA3  of the 
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
!
       SUBROUTINE VBFNLO_HIGGSTO_BETAFCT
       !USE VBFNLO_HIGGSTO_COLOUR 
       !USE VBFNLO_HIGGSTO_BETACOM   
!
       IMPLICIT NONE
       INTEGER NF
       REAL*8 B00,B01,B10,B11
       INTEGER, parameter :: NFMIN=3
       INTEGER, parameter :: NFMAX=6
!
! ---------------------------------------------------------------------
!
! ..The full LO and NLO coefficients 
!
       B00 =  11.d0/3.d0 * CA
       B01 =  -4.d0/3.d0 * TR
       B10 =  34.d0/3.d0 * CA**2
       B11 = -20.d0/3.d0 * CA*TR - 4.* CF*TR
!
! ..Flavour-number loop and output to the array
!
       DO NF = NFMIN,NFMAX
!
        BETA0(NF) = B00 + B01 * NF
        BETA1(NF) = B10 + B11 * NF
!
        BETA2(NF) = 1428.50d0 - 279.611d0 * NF + 6.01852d0 * NF**2
        BETA3(NF) = 29243.0d0 - 6946.30d0 * NF + 405.089d0 * NF**2  &
     &              + 1.49931d0 * NF**3
       ENDDO

!
       RETURN
       END SUBROUTINE VBFNLO_HIGGSTO_BETAFCT
!
! =================================================================av==
!
!--   G.W. VBFNLO_HIGGSTO_DZEROX taken from CERNLIB to find the zero of a function.
!
      FUNCTION VBFNLO_HIGGSTO_DZERO(A0,B0,EPS,MAXF,F,MODE)
      IMPLICIT REAL*8 (A-H,O-Z)
!     Based on
!
!        J.C.P. Bus and T.J. Dekker, Two Efficient Algorithms with
!        Guaranteed Convergence for Finding a Zero of a Function,
!        ACM Trans. Math. Software 1 (1975) 330-345.
!
!        (MODE = 1: Algorithm M;    MODE = 2: Algorithm R)
!
!      CHARACTER*80 ERRTXT
      REAL*8 VBFNLO_HIGGSTO_DZERO
      LOGICAL LMT
      DIMENSION IM1(2),IM2(2),LMT(2)
      integer mode, maxf, mf, im1, im2, ie
      PARAMETER (Z1 = 1, HALF = Z1/2)
      DATA IM1 /2,3/, IM2 /-1,3/
!
      VBFNLO_HIGGSTO_DZERO = 0.d0             ! G.W. to prevent compiler warning
      IF(MODE .NE. 1 .AND. MODE .NE. 2) THEN
       C=0
!       WRITE(ERRTXT,101) MODE
       print 101,MODE
!       WRITE(6,*) ERRTXT
       GO TO 99
      ENDIF
      FA=F(B0)
      FB=F(A0)
      IF(FA*FB .GT. 0) THEN
       C=0
!       WRITE(ERRTXT,102) A0,B0
       print 102,A0,B0
!       WRITE(6,*) ERRTXT
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
       ENDIF
       A=B
       B=C
       C=A
       FA=FB
       FB=FC
       FC=FA
      ENDIF
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
        ENDIF
        IF(P .LT. 0) THEN
         P=-P
         Q=-Q
        ENDIF
        IF(IE .EQ. IM2(MODE)) P=P+P
        IF(P .EQ. 0 .OR. P .LE. Q*TOL) THEN
         W=TOL
        ELSEIF(P .LT. HB*Q) THEN
         W=P/Q
        ELSE
         W=HB
        ENDIF
       ENDIF
       D=A
       A=B
       FD=FA
       FA=FB
       B=B+W
       MF=MF+1
       IF(MF .GT. MAXF) THEN
!        WRITE(6,*) "Error in VBFNLO_HIGGSTO_DZERO: TOO MANY FUNCTION CALLS"
        print*,"Error in VBFNLO_HIGGSTO_DZERO: TOO MANY FUNCTION CALLS"
        GO TO 99
       ENDIF
       FB=F(B)
       IF(FB .EQ. 0 .OR. SIGN(Z1,FC) .EQ. SIGN(Z1,FB)) GO TO 1
       IF(W .EQ. HB) GO TO 2
       IE=IE+1
       GO TO 3
      ENDIF
      VBFNLO_HIGGSTO_DZERO=C
   99 CONTINUE
      RETURN
  101 FORMAT('Error in DZERO: MODE = ',I3,' ILLEGAL')
  102 FORMAT('Error in DZERO: F(A) AND F(B) HAVE THE SAME SIGN, A = ', &
     &     1P,D15.8,', B = ',D15.8)
!
      END FUNCTION VBFNLO_HIGGSTO_DZERO
!
!------------------------------------------------------------------
!
      !MODULE VBFNLO_HIGGSTO_OLAS
      ! CONTAINS
!
      FUNCTION VBFNLO_HIGGSTO_b0af_em(scal,psi,ps0i,xmsi,xm0i) RESULT(value)
      !USE VBFNLO_HIGGSTO_acmplx_pro
      !USE VBFNLO_HIGGSTO_acmplx_rat
      !USE VBFNLO_HIGGSTO_cmplx_root
      !USE VBFNLO_HIGGSTO_cmplx_rootz
      !USE VBFNLO_HIGGSTO_ln_2_riemann
      !USE VBFNLO_HIGGSTO_full_ln
      !USE VBFNLO_HIGGSTO_units
!
      IMPLICIT NONE
!
      REAL*8 scal,ps0i,xm0i,scals,ps0,xm0
      REAL*8, dimension(2) :: value,psi,ps,xmsi,xms,betasc,betas, &
     &        betac,beta,argc,arg,lbet
!
      scals= scal*scal
      ps= psi/scals
      ps0= ps0i/scals
      xms= xmsi/scals
      xm0= xm0i/scal
!
      betasc= co-4.d0*(xms.cq.ps)
      IF(betasc(2).eq.0.d0) THEN
       betasc(2)= -eps
       betac= (betasc(1)).cr.(betasc(2))
      ELSE
       betac= (betasc(1)).crz.(betasc(2))
      ENDIF 
      argc= (betac+co).cq.(betac-co)
      IF(argc(2).eq.0.d0) THEN
       argc(2)= eps
       lbet= argc(1).fln.argc(2)
      ELSE
       betas(1)= 1.d0-4.d0*(xm0*xm0)/ps0
       betas(2)= -eps
       beta= (betas(1)).cr.(betas(2))
       arg= (beta+co).cq.(beta-co)
       IF(arg(2).eq.0.d0) arg(2)= eps
       lbet= argc.lnsrs.arg
      ENDIF
!
      value= 2.d0*co-(betac.cp.lbet)
!
      RETURN
!
      END FUNCTION VBFNLO_HIGGSTO_b0af_em
!
      !END MODULE VBFNLO_HIGGSTO_OLAS
!
!--------------------------------------------------------------------------
!
      FUNCTION VBFNLO_HIGGSTO_quarkQCD(scal,psi,ps0i,xmsi,xm0i,type) RESULT(value)
      !USE VBFNLO_HIGGSTO_riemann
      !USE VBFNLO_HIGGSTO_acmplx_pro
      !USE VBFNLO_HIGGSTO_acmplx_rat
      !USE VBFNLO_HIGGSTO_cmplx_root
      !USE VBFNLO_HIGGSTO_cmplx_rootz
      !USE VBFNLO_HIGGSTO_cmplx_srs_root
      !USE VBFNLO_HIGGSTO_ln_2_riemann
      !USE VBFNLO_HIGGSTO_full_ln
      !USE VBFNLO_HIGGSTO_sp_fun
      !USE VBFNLO_HIGGSTO_units
!
      IMPLICIT NONE
!
      INTEGER it,type,unit
      REAL*8 scal,ps0i,xm0i,scals,ps0,xm0,sgn
      REAL*8, dimension(2) :: value,psi,ps,xmsi,xms,betasc,betas, &
     &        betac,beta,argc,arg,lq,lm,cx,comx,cxs,x,omx,xs,clx,clomx, &
     &        clxs,li2cx,li3cx,li2cxs,li3cxs,copx,clopx,lqs,qcd,lms, &
     &        opx,tau,taus,clxx,clxxs
      REAL*8, dimension(6,2) :: aux,auxs
!
      scals= scal*scal
      ps= psi/scals
      ps0= ps0i/scals
      xms= xmsi/scals
      xm0= xm0i/scal
!
      IF(psi(2).eq.0.d0.and.xmsi(2).eq.0.d0 &
     &   .and.psi(1).le.4.d0*xmsi(1)) THEN
       unit= 1
      ELSE
       unit= 0
      ENDIF
!     
      IF(abs(ps(2)/ps(1)).lt.1.d-10.and.xms(2).eq.0.d0) THEN
       betasc(1)= 1.d0-4.d0*xms(1)/ps(1)
       betasc(2)= 4.d0/(ps(1)*ps(1))*xms(1)*ps(2)
      ELSE
       betasc= co-4.d0*(xms.cq.ps)
      ENDIF
      IF(betasc(2).eq.0.d0) THEN
       betasc(2)= -eps
       betac= (betasc(1)).cr.(betasc(2))
      ELSE
       betac= (betasc(1)).crz.(betasc(2))
      ENDIF
      argc= (betac+co).cq.(betac-co)
!
      betas(1)= 1.d0-4.d0*xm0*xm0/ps0
      betas(2)= -eps
      beta= (betas(1)).cr.(betas(2))
      arg= (beta+co).cq.(beta-co)
!
      IF(arg(2).eq.0.d0) THEN
       x(1)= 1.d0/arg(1) 
       x(2)= -eps      
       sgn= sign(one,x(1))
       xs(1)= x(1)*x(1)
       xs(2)= -sgn*eps
      ELSE
       x= (beta-co).cq.(beta+co)
       xs= x.cp.x  
      ENDIF 
      omx= co-x
      opx= co+x
!
      IF(arg(2).eq.0.d0) arg(2)= eps
      IF(argc(2).eq.0.d0) THEN
       it= 0
       argc(2)= eps
       lq= argc(1).fln.argc(2)
      ELSE
       it= 1
       lq= argc.lnsrs.arg
      ENDIF
      lqs= lq.cp.lq
!
      IF(it.eq.0.d0) THEN
       cx(1)= 1.d0/argc(1) 
       cx(2)= -eps      
       comx= co-cx
       copx= co+cx
       sgn= sign(one,cx(1))
       cxs(1)= cx(1)*cx(1)
       cxs(2)= -sgn*eps
       clx= cx(1).fln.cx(2)
       clxs= clx.cp.clx
       clxx= cxs(1).fln.cxs(2)
       clxxs= clxx.cp.clxx
       clomx= comx(1).fln.comx(2)
       clopx= copx(1).fln.copx(2)
       aux= VBFNLO_HIGGSTO_s_niels_up4(cx) 
       li2cx(1:2)= aux(1,1:2) 
       li3cx(1:2)= aux(2,1:2) 
       auxs= VBFNLO_HIGGSTO_s_niels_up4(cxs) 
       li2cxs(1:2)= auxs(1,1:2) 
       li3cxs(1:2)= auxs(2,1:2) 
      ELSEIF(it.eq.1.d0) THEN
       cx= (betac-co).cq.(betac+co)
       comx= co-cx
       copx= co+cx
       cxs= cx.cp.cx
       clx= cx.lnsrs.x
       clxs= clx.cp.clx
       clxx= cxs.lnsrs.xs
       clxxs= clxx.cp.clxx
       clomx= comx.lnsrs.omx
       clopx= copx.lnsrs.opx
       li2cx= VBFNLO_HIGGSTO_li2_srsz(cx,x,unit)
       li3cx= VBFNLO_HIGGSTO_li3_srsz(cx,x,unit)
       li2cxs= VBFNLO_HIGGSTO_li2_srsz(cxs,xs,unit)
       li3cxs= VBFNLO_HIGGSTO_li3_srsz(cxs,xs,unit)
      ENDIF
!
      IF(xms(2).eq.0.d0) THEN
       lm(1)= log(xms(1))
       lm(2)= 0.d0
      ELSE  
       lm= xms(1).fln.xms(2)
      ENDIF
      lms= lm.cp.lm
!
      tau= xms.cq.ps
      taus= tau.cp.tau
!
      IF(type.eq.0) THEN
!
       qcd= -3.d0/4.d0*(co-12.d0*tau)*rz2 &
     & +1.d0/16.d0*(3.d0*co+344.d0*tau) &
     & -3.d0/2.d0*((co-6.d0*tau).cp.lms) &
     & +3.d0/4.d0*((3.d0*co-14.d0*tau).cp.(betac.cp.clx)) &
     & -3.d0*((li3cxs-2.d0*li3cx+4.d0/3.d0*(li2cx.cp.clx) &
     & +1.d0/3.d0*(clxs.cp.clomx)+rz3*co &
     & -2.d0/3.d0*(clxx.cp.li2cxs)-1.d0/6.d0*(clxxs.cp.clomx) &
     & -1.d0/6.d0* &
     &  (clxxs.cp.clopx)).cp.((co-4.d0*tau).cp.(co-2.d0*tau))) &
     & -1.d0/2.d0*((4.d0*li2cxs-4.d0*li2cx-4.d0*(clomx.cp.clx) &
     & -2.d0*(cx.cp.(clxs.cq.comx))+4.d0*(clxx.cp.clomx) &
     & +4.d0*(clxx.cp.clopx)+(clxxs.cp.(cxs.cq.comx)) &
     & +(clxxs.cp.(cxs.cq.copx))).cp.(betac.cp.(co-4.d0*tau))) &
     & +1.d0/4.d0*(lm.cp.(11.d0*co-108.d0*tau)) &
     & +1.d0/4.d0*(clxs.cp.(3.d0*co+58.d0*taus-28.d0*tau))
!  
      ELSEIF(type.eq.1) THEN
!
      qcd= -3.d0/4.d0*(co-12.d0*tau)*rz2 &
     & +1.d0/16.d0*(67.d0*co-40.d0*tau) &
     & +3.d0/4.d0*((3.d0*co-14.d0*tau).cp.(betac.cp.clx)) &
     & -3.d0*((li3cxs-2.d0*li3cx+4.d0/3.d0*(li2cx.cp.clx) &
     & +1.d0/3.d0*(clxs.cp.clomx)+rz3*co &
     & -2.d0/3.d0*(clxx.cp.li2cxs)-1.d0/6.d0*(clxxs.cp.clomx) &
     & -1.d0/6.d0* &
     &  (clxxs.cp.clopx)).cp.((co-4.d0*tau).cp.(co-2.d0*tau))) &
     & -1.d0/2.d0*((4.d0*li2cxs-4.d0*li2cx-4.d0*(clomx.cp.clx) &
     & -2.d0*(cx.cp.(clxs.cq.comx))+4.d0*(clxx.cp.clomx) &
     & +4.d0*(clxx.cp.clopx)+(clxxs.cp.(cxs.cq.comx)) &
     & +(clxxs.cp.(cxs.cq.copx))).cp.(betac.cp.(co-4.d0*tau))) &
     & -2.d0* &
     &  (((co-8.d0*tau).cp.(betac.cp.lq)).cp.(co-3.d0/4.d0*lm)) &
     & -3.d0*((betac.cp.lq).cp.(lm.cp.tau)) &
     & +4.d0*((co.cp.betac).cp.(lq.cp.tau)) &
     & -9.d0/4.d0*(lm.cp.(co+4.d0*tau)) &
     & +9.d0*(lms.cp.tau) &
     & +1.d0/4.d0*(clxs.cp.(3.d0*co+58.d0*taus-28.d0*tau))
!
      ENDIF
!
      value= qcd
!
      RETURN
!
      END FUNCTION VBFNLO_HIGGSTO_quarkQCD
!
!------------------------------------------------------------------
!
!-----------------------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------------------
!
!
      !MODULE VBFNLO_HIGGSTO_qcd
      ! CONTAINS
!
       FUNCTION VBFNLO_HIGGSTO_ralphas(rs0,rs,als)
       !USE VBFNLO_HIGGSTO_riemann
       !USE VBFNLO_HIGGSTO_masses
       IMPLICIT NONE
       INTEGER i,nfe
       REAL*8 VBFNLO_HIGGSTO_ralphas,rs0,rs,als,x1,x2,xacc,qcdl,rat,rl,rll,rb,fac, &
     &        facp,rhs,ratl2,pqcdl4,qcdb0,qcdb1,qcdb2,qcda,pqcdl3, &
     &        fac2,pqcdl5
       REAL*8, dimension(5) :: b0,b1,b2
       INTEGER, parameter :: nf=5 
!
!-----limits for lambda_5 are 1 mev < lambda_5 < 10 gev
!
       IF(als.eq.0.d0) THEN
        VBFNLO_HIGGSTO_ralphas= 0.d0
       ELSE
        x1= 0.001d0
        x2= 10.0d0
        xacc= 1.d-12
        qcdl= VBFNLO_HIGGSTO_qcdlam(nf,als,rs0,x1,x2,xacc)
        pqcdl5= qcdl
        DO i=1,5
         b0(i)= (11.d0-2.d0/3.d0*i)/4.d0
         b1(i)= (102.d0-38.d0/3.d0*i)/16.d0
         b2(i)= 0.5d0*(2857.d0-i*(5033.d0/9.d0-325.d0/27.d0*i))/64.d0
        ENDDO
!
        IF(rs < mbq) THEN
         rat= mbq/qcdl
         rl= 2.d0*log(rat)
         rll= log(rl)
         rb= log(b0(5)/b0(4))
         fac= b1(5)/b0(5)-b1(4)/b0(4)
         facp= b2(5)/b0(5)-b2(4)/b0(4)
         fac2= b1(5)*b1(5)/b0(5)/b0(5)-b1(4)*b1(4)/b0(4)/b0(4)
         rhs= (b0(5)-b0(4))*rl+fac*rll-b1(4)/b0(4)*rb+ &
     &         b1(5)/b0(5)/b0(5)*fac*rll/rl+1.d0/b0(5)/rl*( &
     &         fac2-facp-7.d0/72.d0)
         rhs= rhs/b0(4)
         ratl2= exp(rhs)
         qcdl= qcdl/sqrt(ratl2)      
         pqcdl4= qcdl
         nfe= nf-1
         IF(rs < mcq) THEN
          rat= mcq/qcdl
          rl= 2.d0*log(rat)
          rll= log(rl)
          rb= log(b0(4)/b0(3))
          fac= b1(4)/b0(4)-b1(3)/b0(3)
          facp= b2(4)/b0(4)-b2(3)/b0(3)
          fac2= b1(4)*b1(4)/b0(4)/b0(4)-b1(3)*b1(3)/b0(3)/b0(3)
          rhs= (b0(4)-b0(3))*rl+fac*rll-b1(3)/b0(3)*rb+ &
     &          b1(4)/b0(4)/b0(4)*fac*rll/rl+1.d0/b0(4)/rl*( &
     &          fac2-facp-7.d0/72.d0)
          rhs= rhs/b0(3)
          ratl2= exp(rhs)
          qcdl= qcdl/sqrt(ratl2)      
          pqcdl3= qcdl
          nfe= nf-2
         ENDIF
        ELSE
         nfe= nf
        ENDIF
!
        qcdb0= 11.d0-2.d0/3.d0*nfe
        qcdb1= 102.d0-38.d0/3.d0*nfe
        qcdb2= 0.5d0*(2857.d0-5033.d0/9.d0*nfe+325.d0/27.d0*nfe*nfe)
        qcda= 2.d0*log(rs/qcdl)
!
        VBFNLO_HIGGSTO_ralphas= 4.d0*pi/qcdb0/qcda*(1.d0-qcdb1/qcdb0**2/qcda* &
     &           log(qcda)+(qcdb1/qcdb0**2/qcda)**2*((log(qcda)- &
     &           0.5d0)**2+qcdb2*qcdb0/qcdb1**2-5.d0/4.d0))
!
       ENDIF
       RETURN
       END FUNCTION VBFNLO_HIGGSTO_ralphas
!
!---------------------------------------------------------------------------
!
       FUNCTION VBFNLO_HIGGSTO_qcdlam(nf,als,rs,x1,x2,xacc)
       IMPLICIT NONE
       INTEGER nf,j
       REAL*8 VBFNLO_HIGGSTO_qcdlam,als,rs,x1,x2,xacc,fmid,f,dx,xmid
       INTEGER, parameter :: jmax=50
!
       fmid= VBFNLO_HIGGSTO_qcdscale(nf,als,rs,x2)
       f= VBFNLO_HIGGSTO_qcdscale(nf,als,rs,x1)
       IF(f*fmid >= 0.d0) STOP
       IF(f < 0.d0) THEN
        VBFNLO_HIGGSTO_qcdlam= x1
        dx= x2-x1
       ELSE
        VBFNLO_HIGGSTO_qcdlam= x2
        dx= x1-x2
       ENDIF
       DO j=1,jmax
        dx= dx*0.5d0
        xmid= VBFNLO_HIGGSTO_qcdlam+dx
        fmid= VBFNLO_HIGGSTO_qcdscale(nf,als,rs,xmid)
        IF(fmid <= 0.d0) VBFNLO_HIGGSTO_qcdlam= xmid
        IF(abs(dx) < xacc.or.fmid.eq.0.d0) RETURN
       ENDDO
       END FUNCTION VBFNLO_HIGGSTO_qcdlam  
!
       FUNCTION VBFNLO_HIGGSTO_qcdscale(nf,als,rs,x)
       !USE VBFNLO_HIGGSTO_riemann
       IMPLICIT NONE
       INTEGER nf
       REAL*8 als,rs,x,qcdb0,qcdb1,qcdb2,qcda,VBFNLO_HIGGSTO_qcdscale
!
       qcdb0= 11.d0-2.d0/3.d0*nf
       qcdb1= 102.d0-38.d0/3.d0*nf
       qcdb2= 0.5d0*(2857.d0-5033.d0/9.d0*nf+325.d0/27.d0*nf*nf)
       qcda= 2.d0*log(rs/x)
       VBFNLO_HIGGSTO_qcdscale= als-(4.d0*pi/qcdb0/qcda*(1.d0-qcdb1/qcdb0**2/qcda* &
     &             log(qcda)+(qcdb1/qcdb0**2/qcda)**2*((log(qcda)- &
     &             0.5d0)**2+qcdb2*qcdb0/qcdb1**2-5.d0/4.d0)))
       RETURN
       END FUNCTION VBFNLO_HIGGSTO_qcdscale
!
!--------------------------------------------------------------------
!
      FUNCTION VBFNLO_HIGGSTO_run_bc(scal)   RESULT(value)
      !USE VBFNLO_HIGGSTO_riemann
      !USE VBFNLO_HIGGSTO_masses
      !USE VBFNLO_HIGGSTO_set_phys_const
!
      IMPLICIT NONE
!
      REAL*8 scal,scal2,alsr,aexps,als1g,alsc4,alsc42, &
     &       als1g2,alsb,alsb2,alsz,alsz2,fn,b03,b13,b23,g03,g13, &
     &       g23,rex3,b03s,b03c,b13s,sm1g,cfm13,cfm23, &
     &       b04,b04c,b14,b24,g04,g14,g24,rex4,b04s,b14s,asldf,cfm14, &
     &       cfm24,rsmb,zero,x1,x2,xacc,cmm4,cmb,b05,b15,b25,g05,g15, &
     &       g25,rex5,b05s,b05c,b15s,cfm15,cfm25,rcqm,bmm5,rbqm,scal1g, &
     &       rsmc
      REAL*8, dimension(2) :: value
!
      scal2= scal*scal
      alsr= VBFNLO_HIGGSTO_ralphas(mz,scal,als)
      aexps= alsr/pi
!
!-----COMPUTES THE RUNNING CHARM MASS
!
      scal1g= 1.d0
      als1g= VBFNLO_HIGGSTO_ralphas(mz,scal1g,als)/pi
      alsc4= VBFNLO_HIGGSTO_ralphas(mz,mcq,als)/pi
      alsc42= alsc4*alsc4
      als1g2= als1g*als1g
      alsb= VBFNLO_HIGGSTO_ralphas(mz,mbq,als)/pi
      alsb2= alsb*alsb
      alsz= VBFNLO_HIGGSTO_ralphas(mz,scal,als)/pi
      alsz2= alsz*alsz
!
!-----FIRST THE RUNNING OF THE S-QUARK MASS
!     UP TO C/B-THRESHOLD
!
      fn= 3.d0
      b03= (11.d0-2.d0/3.d0*fn)/4.d0
      b13= (102.d0-38.d0/3.d0*fn)/16.d0
      b23= (2857.d0/2.d0-5033.d0/18.d0*fn+325.d0/54.d0*fn*fn)/64.d0
      g03= 1.d0
      g13= (202.d0/3.d0-20.d0/9.d0*fn)/16.d0
      g23= (1249.d0-(2216.d0/27.d0+160.d0/3.d0*rz3)*fn- &
     &     140.d0/81.d0*fn*fn)/64.d0
      rex3= g03/b03
      b03s= b03*b03
      b03c= b03s*b03
      b13s= b13*b13
      sm1g= 0.189d0
      asldf= alsc4-als1g
      cfm13= g13/b03-b13*g03/b03s
      cfm23= g23/b03-b13*g13/b03s-b23*g03/b03s+b13s*g03/b03c
      rsmc= sm1g*(alsc4/als1g)**rex3*(1.d0+cfm13* &
     &      asldf+0.5d0*cfm13*cfm13*asldf*asldf+0.5d0*cfm23* &
     &      (alsc42-als1g2))
!
      fn= 4.d0
      b04= (11.d0-2.d0/3.d0*fn)/4.d0
      b14= (102.d0-38.d0/3.d0*fn)/16.d0
      b24= (2857.d0/2.d0-5033.d0/18.d0*fn+325.d0/54.d0*fn*fn)/64.d0
      g04= 1.d0
      g14= (202.d0/3.d0-20.d0/9.d0*fn)/16.d0
      g24= (1249.d0-(2216.d0/27.d0+160.d0/3.d0*rz3)*fn- &
     &      140.d0/81.d0*fn*fn)/64.d0
      rex4= g04/b04
      b04s= b04*b04
      b04c= b04s*b04
      b14s= b14*b14
      asldf= alsb-alsc4
      cfm14= g14/b04-b14*g04/b04s
      cfm24= g24/b04-b14*g14/b04s-b24*g04/b04s+ &
     &       b14s*g04/b04c
      rsmb= rsmc*(alsb/alsc4)**rex4*(1.d0+cfm14* &
     &      asldf+0.5d0*cfm14*cfm14*asldf*asldf+0.5d0*cfm24* &
     &      (alsb2-alsc42))
!
!-----C QUARK mASS AT C-THRESHOLD
!
      zero= 0.d0
      x1= 0.5d0
      x2= 2.0d0
      xacc= 1.d-12
      fn= 4.d0
      cmm4= VBFNLO_HIGGSTO_rrunm(x1,x2,xacc,mcq,alsc4,rsmc,zero,fn)
!
!-----C QUARK MASS AT B-THRESHOLD
!
      cmb= cmm4*(alsb/alsc4)**rex4*(1.d0+cfm14* &
     &     (alsb-alsc4)+0.5d0*cfm14*cfm14* &
     &     (alsb-alsc4)**2+0.5d0*cfm24*(alsb2-alsc42))
!
!-----RUNNING CHARM MASS
!
      fn= 5.d0
      b05= (11.d0-2.d0/3.d0*fn)/4.d0
      b15= (102.d0-38.d0/3.d0*fn)/16.d0
      b25= (2857.d0/2.d0-5033.d0/18.d0*fn+325.d0/54.d0*fn*fn)/64.d0
      g05= 1.d0
      g15= (202.d0/3.d0-20.d0/9.d0*fn)/16.d0
      g25= (1249.d0-(2216.d0/27.d0+160.d0/3.d0*rz3)*fn- &
     &      140.d0/81.d0*fn*fn)/64.d0
      rex5= g05/b05
      b05s= b05*b05
      b05c= b05s*b05
      b15s= b15*b15
      cfm15= g15/b05-b15*g05/b05s
      cfm25= g25/b05-b15*g15/b05s-b25*g05/b05s+b15s*g05/b05c
      rcqm= cmb*(alsz/alsb)**rex5*(1.d0+cfm15* &
     &     (alsz-alsb)+0.5d0*cfm15*cfm15* &
     &     (alsz-alsb)**2+0.5d0*cfm25*(alsz2-alsb2))
!
!-----B QUARK MASS AT B-THRESHOLD
!
      x1= 0.5d0
      x2= 6.0d0
      xacc= 1.d-12
      fn= 5.d0
      bmm5= VBFNLO_HIGGSTO_rrunm(x1,x2,xacc,mbq,alsb,cmb,rsmb,fn)
!
!-----RUNNINg B MASS
!
      asldf= alsz-alsb
      rbqm= bmm5*(alsz/alsb)**rex5*(1.d0+cfm15* &
     &      asldf+0.5d0*cfm15*cfm15*asldf*asldf+0.5d0*cfm25* &
     &      (alsz2-alsb2))
!
      value(1)= rcqm
      value(2)= rbqm
!
      RETURN
      END FUNCTION VBFNLO_HIGGSTO_run_bc
!
!
!-----RRUNM--------------------------------------------------------
!     COMPUTES THE RUNNING QUARK MASS AT THE POLE MASS
!
      FUNCTION VBFNLO_HIGGSTO_rrunm(x1,x2,xacc,qm,als,rm1,rm2,fn)
      IMPLICIT NONE
!
      INTEGER, parameter :: jmax=50
      INTEGER j
      REAL*8 VBFNLO_HIGGSTO_rrunm,x1,x2,xacc,qm,als,rm1,rm2,fn,fmid,f,dx,xmid
!
      fmid= VBFNLO_HIGGSTO_qcdmass(qm,als,rm1,rm2,fn,x2)
      f= VBFNLO_HIGGSTO_qcdmass(qm,als,rm1,rm2,fn,x1)
      IF(f*fmid.ge.0.d0) then
       print*,'root must be bracketed for bisection'
       print 1,qm
    1  format(/' error detected by VBFNLO_HIGGSTO_rrunm ',/ &
     &         ' current value of quark mass  = ',e20.5)
       STOP
      ENDIF
      IF(f < 0.d0) then
       VBFNLO_HIGGSTO_rrunm= x1
       dx= x2-x1
      ELSE
       VBFNLO_HIGGSTO_rrunm= x2
       dx= x1-x2
      ENDIF
      DO j=1,jmax
       dx= dx*0.5d0
       xmid= VBFNLO_HIGGSTO_rrunm+dx
       fmid= VBFNLO_HIGGSTO_qcdmass(qm,als,rm1,rm2,fn,xmid)
       IF(fmid.le.0.d0) VBFNLO_HIGGSTO_rrunm= xmid
       IF(abs(dx) < xacc.or.fmid.eq.0.d0) RETURN
      ENDDO
!
      END FUNCTION VBFNLO_HIGGSTO_rrunm
!
!---------------------------------------------------------------------
!
      FUNCTION VBFNLO_HIGGSTO_qcdmass(qm,als,rm1,rm2,fn,x)
      !USE VBFNLO_HIGGSTO_riemann
!
      IMPLICIT NONE
!
      REAL*8 VBFNLO_HIGGSTO_qcdmass,qm,als,rm1,rm2,fn,x,rln,rlns,r1,r2,delta0, &
     &       delta1,delta2,rhs
!
      rln= 2.d0*log(qm/x)
      rlns= rln*rln
      r1= rm1/x
      r2= rm2/x
      delta0= 3.d0/4.d0*rz2-3.d0/8.d0
      delta1= r1*(pis/8.d0+r1*(-0.597d0+0.230d0*r1))
      delta2= r2*(pis/8.d0+r2*(-0.597d0+0.230d0*r2))
!
      rhs= 1.d0+als*(4.d0/3.d0+rln+als*(3817.d0/288.d0-8.d0/3.d0+ &
     &     2.d0/3.d0*(2.d0+log(2.d0))*rz2-rz3/6.d0-fn/3.d0*(rz2+ &
     &     71.d0/48.d0)+4.d0/3.d0*(delta0+delta1+delta2)+(173.d0/ &
     &     24.d0-13.d0/36.d0*fn)*rln+(15.d0/8.d0-fn/12.d0)*rlns))
      VBFNLO_HIGGSTO_qcdmass= qm-x*rhs
!
      RETURN
!
      END FUNCTION VBFNLO_HIGGSTO_qcdmass
!
      !END MODULE VBFNLO_HIGGSTO_qcd
!
!-----------------------------------------------------------------------------
!
      !MODULE VBFNLO_HIGGSTO_a_cmplx
      ! CONTAINS
!
      FUNCTION VBFNLO_HIGGSTO_b021_dm_cp(scal,psi,ps0i,xmsi,xm0i) RESULT(value)
      !USE VBFNLO_HIGGSTO_cmplx_root
      !USE VBFNLO_HIGGSTO_cmplx_rootz
      !USE VBFNLO_HIGGSTO_cmplx_srs_root
      !USE VBFNLO_HIGGSTO_ln_2_riemann
      !USE VBFNLO_HIGGSTO_acmplx_pro
      !USE VBFNLO_HIGGSTO_acmplx_rat
      !USE VBFNLO_HIGGSTO_full_ln
      !USE VBFNLO_HIGGSTO_units
!
      IMPLICIT NONE
!

      REAL*8 scal,ps0i,ps0,xm1,xm2
      REAL*8, intent(in), dimension(2) :: xm0i
      REAL*8, intent(in), dimension(2,2) :: xmsi
      REAL*8, dimension(2) :: value,xm1c,xm2c
      REAL*8, dimension(2,2) :: xms
      REAL*8, dimension(2) :: xm0
      REAL*8, dimension(2) :: psi,ps,aroot,root,lambdasc,lambdas, &
     &        lambdac,lambda,argc,arg,llam,l1,l2 
!
      ps= psi
      ps0= ps0i
      xms= xmsi
      xm0= xm0i
!
      xm1c(1:2)= xms(1,1:2)
      xm2c(1:2)= xms(2,1:2)
      xm1c= xm1c.cq.ps
      xm2c= xm2c.cq.ps
      xm1= xm0(1)*xm0(1)/ps0
      xm2= xm0(2)*xm0(2)/ps0
!
      aroot= xm1c.cp.xm2c
      root= (aroot(1).crz.aroot(2))
!
      lambdasc= co+(xm1c.cp.xm1c)+(xm2c.cp.xm2c)-2.d0*( &
     &          xm1c+xm2c+(xm1c.cp.xm2c))
      lambdas(1)= 1.d0+xm1*xm1+xm2*xm2-2.d0*( &
     &            xm1+xm2+xm1*xm2)
      lambdas(2)= -eps
      lambdac= (lambdasc(1)).crz.(lambdasc(2))
      lambda= (lambdas(1)).cr.(lambdas(2))
      IF(lambda(2).eq.0.d0) lambda(2)= -eps
!
      argc= 0.5d0*((-co+xm1c+xm2c-lambdac).cq.root)
!
      arg(1)= 0.5d0*(-1.d0+xm1+xm2-lambda(1))/sqrt(xm1*xm2)
      arg(2)= eps
!
      llam= argc.lnsrs.arg
!
      l1= xm1c(1).fln.xm1c(2)
      l2= xm2c(1).fln.xm2c(2)
!
      value= (((xm1c-xm2c-co).cq.lambdac).cp.llam)+0.5d0*(l1-l2) 
!
      RETURN
      END FUNCTION VBFNLO_HIGGSTO_b021_dm_cp
!
      !END MODULE VBFNLO_HIGGSTO_a_cmplx
!
!-----------------------------------------------------------------------------------------
!
      !MODULE VBFNLO_HIGGSTO_root_find2
      ! CONTAINS
!
      FUNCTION VBFNLO_HIGGSTO_zeroin(f, ax, bx, aerr, rerr) RESULT(fn_val)
!
      !USE VBFNLO_HIGGSTO_rootW
! 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-07-14  Time: 12:32:54
 
!-----------------------------------------------------------------------

!         FINDING A ZERO OF THE FUNCTION F(X) IN THE INTERVAL (AX,BX)

!                       ------------------------

!  INPUT...

!  F      FUNCTION SUBPROGRAM WHICH EVALUATES F(X) FOR ANY X IN THE
!         CLOSED INTERVAL (AX,BX).  IT IS ASSUMED THAT F IS CONTINUOUS,
!         AND THAT F(AX) AND F(BX) HAVE DIFFERENT SIGNS.
!  AX     LEFT ENDPOINT OF THE INTERVAL
!  BX     RIGHT ENDPOINT OF THE INTERVAL
!  AERR   THE ABSOLUTE ERROR TOLERANCE TO BE SATISFIED
!  RERR   THE RELATIVE ERROR TOLERANCE TO BE SATISFIED

!  OUTPUT...

!         ABCISSA APPROXIMATING A ZERO OF F IN THE INTERVAL (AX,BX)

!-----------------------------------------------------------------------
!  ZEROIN IS A SLIGHTLY MODIFIED TRANSLATION OF THE ALGOL PROCEDURE
!  ZERO GIVEN BY RICHARD BRENT IN ALGORITHMS FOR MINIMIZATION WITHOUT
!  DERIVATIVES, PRENTICE-HALL, INC. (1973).
!-----------------------------------------------------------------------
      
      IMPLICIT NONE
      REAL*8, INTENT(IN)  :: ax
      REAL*8, INTENT(IN)  :: bx
      REAL*8, INTENT(IN)  :: aerr
      REAL*8, INTENT(IN)  :: rerr
      REAL*8              :: fn_val
!
! EXTERNAL f
!
      INTERFACE
       FUNCTION f(x) RESULT(fn_val)
       IMPLICIT NONE
       REAL*8, INTENT(IN)  :: x
       REAL*8              :: fn_val
       END FUNCTION f
      END INTERFACE
!
      REAL*8 :: a,b,c,d,e,eps,fa,fb,fc,tol,xm,p,q,r,s,atol,rtol

!  COMPUTE EPS, THE RELATIVE MACHINE PRECISION

      eps= EPSILON(0.0d0)

! INITIALIZATION
      
      a= ax
      b= bx
      fa= f(a)
      fb= f(b)
      IF(fa*fb.gt.0.d0) THEN
       inc= 1
      ENDIF
      atol= 0.5d0*aerr
      rtol= MAX(0.5d0*rerr, 2.0d0*eps)

! BEGIN STEP

   10 c= a
      fc= fa
      d= b-a
      e= d
   20 IF(ABS(fc) < ABS(fb)) THEN
         a= b
         b= c
         c= a
         fa= fb
         fb= fc
         fc= fa
      END IF

! CONVERGENCE TEST

      tol= rtol*MAX(ABS(b),ABS(c))+atol
      xm= 0.5d0*(c-b)
      IF(ABS(xm) > tol) THEN
       IF(fb /= 0.0d0) THEN
            
! IS BISECTION NECESSARY
    
        IF(ABS(e) >= tol) THEN
         IF(ABS(fa) > ABS(fb)) THEN
        
! IS QUADRATIC INTERPOLATION POSSIBLE
        
          IF(a == c) THEN
          
! LINEAR INTERPOLATION
          
           s= fb/fc
           p= (c-b)*s
           q= 1.0d0-s
          ELSE
          
! INVERSE QUADRATIC INTERPOLATION
          
           q= fa/fc
           r= fb/fc
           s= fb/fa
           p= s*((c-b)*q*(q-r)-(b-a)*(r-1.0d0))
           q= (q-1.0d0)*(r-1.0d0)*(s-1.0d0)
          ENDIF
        
! ADJUST SIGNS
        
          IF(p > 0.0d0) q= -q
          p= ABS(p)
        
! IS INTERPOLATION ACCEPTABLE
        
          IF(2.0*p < (3.0*xm*q-ABS(tol*q))) THEN
           IF(p < ABS(0.5d0*e*q)) THEN
            e= d
            d= p/q
            GO TO 30
           ENDIF
          ENDIF
         ENDIF
        ENDIF
    
! BISECTION
    
        d= xm
        e= d
    
! COMPLETE STEP
    
  30    a= b
        fa= fb
        IF(ABS(d) > tol) b= b+d
        IF(ABS(d) <= tol) b= b+SIGN(tol,xm)
        fb= f(b)
        IF(fb*(fc/ABS(fc)) > 0.0d0) GO TO 10
        GO TO 20
       ENDIF
      ENDIF

! DONE

      fn_val= b
      RETURN
      END FUNCTION VBFNLO_HIGGSTO_zeroin
!
      !END MODULE VBFNLO_HIGGSTO_root_find2
!
!--------------------------------------------------------------------------
!
      !MODULE VBFNLO_HIGGSTO_HBB_cp
      ! CONTAINS
!   
      FUNCTION VBFNLO_HIGGSTO_SSHH(x)
      !USE VBFNLO_HIGGSTO_aux_Hcp
!
      IMPLICIT NONE
!
      REAL*8, INTENT(IN)  :: x
      REAL*8 VBFNLO_HIGGSTO_SSHH,rgh,muc,scalc
!
      muc= muhcp
      scalc= scalec
      rgh= x*muc
!
      VBFNLO_HIGGSTO_SSHH= VBFNLO_HIGGSTO_SHH(muc,scalc,rgh)
!
      RETURN
      END FUNCTION VBFNLO_HIGGSTO_SSHH
!
!---------------------------------------------------------------------------------------
!
       FUNCTION VBFNLO_HIGGSTO_SHH(muhr,scal,rgh)   
       !USE VBFNLO_HIGGSTO_ferbos
       !USE VBFNLO_HIGGSTO_a_cmplx
       !USE VBFNLO_HIGGSTO_aux_Hcp
       !USE VBFNLO_HIGGSTO_aux_Hbb
       !USE VBFNLO_HIGGSTO_units
       !USE VBFNLO_HIGGSTO_acmplx_pro
       !USE VBFNLO_HIGGSTO_acmplx_rat
       !USE VBFNLO_HIGGSTO_full_ln
       !USE VBFNLO_HIGGSTO_ln_2_riemann
       !USE VBFNLO_HIGGSTO_masses
       !USE VBFNLO_HIGGSTO_set_phys_const
       !USE VBFNLO_HIGGSTO_riemann
       !USE VBFNLO_HIGGSTO_optcp
       !USE VBFNLO_HIGGSTO_OLAS
       !USE VBFNLO_HIGGSTO_common_niels
       !USE VBFNLO_HIGGSTO_QCD
       !USE VBFNLO_HIGGSTO_cmplx_rootz
!
       IMPLICIT NONE
!
       INTEGER iz
       REAL*8 VBFNLO_HIGGSTO_SHH,muhr,rgh,muhs,scal,scals,p2,xm0,str,sti,EWC, &
     &        sconv,asmur,emc,emb,emt,as_NLO, &
     &        crmbs,crmcs,lcxb,lcxc,lcxbs,lcxcs,lclxb,lclxc,qcdtop,neg
       REAL*8, dimension(2) :: axm0
       REAL*8, dimension(2,2) :: axms


       REAL*8, dimension(2) :: sh,shs,clh,b0sumb,b0sumf,cxp,ksumb, &
     &         ksumf,coefB1,coefB2,coefB3,coefB4,coefB5,coefB6,coefB7, &
     &         coefB8,coefB9,coefB10,coefB11,coefB12, &
     &         totalf,totalt,totalb,total,b0sumt,ksumt,xms, &
     &         b0part,cpt,ccxt,deltag,sww0,coefW1,coefW2, &
     &         ksumw,ccxts,cctsi,shi,clt,cxhw,cstsi,csmcts,cltmw, &
     &         b0sumw,sww,DW,b0sumw1,b0sumw2,cmxw,clmw,cxtw, &
     &         nloqcd,runbc,ttqcd,wccxt,wclt,wccxts, &
     &         wcxtw,wcltmw
!     
       ! INTERFACE
       !  SUBROUTINE VBFNLO_HIGGSTO_INITALPHAS(asord,FR2,MUR,asmur,emc,emb,emt)
       !  !USE VBFNLO_HIGGSTO_DZpar 
       !  IMPLICIT NONE
       !  INTEGER asord
       !  REAL*8 FR2,MUR,asmur,emc,emb,emt,VBFNLO_HIGGSTO_FINDALPHASR0
       !  EXTERNAL VBFNLO_HIGGSTO_FINDALPHASR0
       !  END SUBROUTINE VBFNLO_HIGGSTO_INITALPHAS
       ! END INTERFACE
! !
       ! INTERFACE
       !  FUNCTION VBFNLO_HIGGSTO_ALPHAS(MUR)
       !  !USE VBFNLO_HIGGSTO_NFFIX  
       !  !USE VBFNLO_HIGGSTO_VARFLV 
       !  !USE VBFNLO_HIGGSTO_FRRAT  
       !  !USE VBFNLO_HIGGSTO_ASINP  
       !  !USE VBFNLO_HIGGSTO_ASFTHR 
       !  IMPLICIT NONE
       !  REAL*8 MUR,VBFNLO_HIGGSTO_ALPHAS
       !  END FUNCTION VBFNLO_HIGGSTO_ALPHAS
       ! END INTERFACE
! !
       ! INTERFACE
       !  FUNCTION VBFNLO_HIGGSTO_quarkQCD(scal,psi,ps0i,xmsi,xm0i,type)  &
     ! &                        RESULT(value)
       !  !USE VBFNLO_HIGGSTO_riemann
       !  !USE VBFNLO_HIGGSTO_acmplx_pro
       !  !USE VBFNLO_HIGGSTO_acmplx_rat
       !  !USE VBFNLO_HIGGSTO_cmplx_root
       !  !USE VBFNLO_HIGGSTO_cmplx_rootz
       !  !USE VBFNLO_HIGGSTO_cmplx_srs_root
       !  !USE VBFNLO_HIGGSTO_ln_2_riemann
       !  !USE VBFNLO_HIGGSTO_full_ln
       !  !USE VBFNLO_HIGGSTO_sp_fun
       !  !USE VBFNLO_HIGGSTO_units
       !  IMPLICIT NONE
       !  INTEGER type
       !  REAL*8 scal,ps0i,xm0i
       !  REAL*8, dimension(2) :: value,psi,xmsi
       !  END FUNCTION VBFNLO_HIGGSTO_quarkQCD
       ! END INTERFACE
!
       muhs= muhr*muhr
       scals= scal*scal
!
       asmur= 0.12018d0
       emc= 1.4d0
       emb= 4.75d0
       emt= mt
       iz= 1
       CALL VBFNLO_HIGGSTO_INITALPHAS(iz,one,mz,asmur,emc,emb,emt)
       als= VBFNLO_HIGGSTO_ALPHAS(scal)
       as_NLO= als/pi
!
       runbc= VBFNLO_HIGGSTO_run_bc(scal)
       crmbs= runbc(2)*runbc(2)
       crmcs= runbc(1)*runbc(1)         
!
       lcxb= crmbs/scals
       lcxc= crmcs/scals
!
       lcxbs= lcxb*lcxb
       lcxcs= lcxc*lcxc
       lclxb= log(lcxb)
       lclxc= log(lcxc)
!
       IF(gtop == 1) THEN
        imt= g_f/(8.d0*sqrt(2.d0)*pi)*(mt*mt-mw*mw)**2* &
     &      (mt*mt+2.d0*mw*mw)/mt**3
       ELSEIF(gtop == 0) THEN 
        imt= 0.d0
       ELSE
        imt= yimt
       ENDIF
!
       str= mt*mt-0.25d0*imt*imt
       sti= -mt*imt
       cpt(1)= str
       cpt(2)= sti
       wccxt= cpt/scals
       wccxts= wccxt.cp.wccxt
       ccxt= cpt/scals
       ccxts= ccxt.cp.ccxt
!
       cctsi= co.cq.ccts
       cstsi= co.cq.csts
       csmcts= csts-ccts
       cmxw= -cxw
!
       clmw= cmxw(1).fln.cmxw(2)
       clmw(2)= clmw(2)-2.d0*pi
!
       sh(1)= muhs/scals
       sh(2)= -muhr*rgh/scals
!
       cxhw= sh-cxw 
       shs= sh.cp.sh
       shi= co.cq.sh
!
       clh= sh(1).fln.sh(2)
       clt= ccxt(1).fln.ccxt(2)
       wclt= wccxt(1).fln.wccxt(2)
       cxtw= ccxt-cxw
       wcxtw= wccxt-cxw
       cltmw= cxtw(1).fln.cxtw(2)
       wcltmw= wcxtw(1).fln.wcxtw(2)
!
! W
       coefB1= (12.d0*cxw-4.d0*sh+(shs.cq.cxw))/64.d0
!
! Z
       coefB2= (-4.d0*(sh.cq.ccts)+(shs.cq.cxw)+ &
     &         12.d0*(cxw.cq.cctq))/128.d0
!  
! H
       coefB3= 9.d0/128.d0*(shs.cq.cxw)
!
! top
       coefB4= -3.d0/32.d0*((4.d0*ccxt-sh).cp.(ccxt.cq.cxw))
!
! light fermions
!
       coefB5= -3.d0/32.d0*((4.d0*co*lcxb-sh).cq.cxw)*lcxb
!
       coefB6= -1.d0/32.d0*((4.d0*co*cxtau-sh).cq.cxw)*cxtau
!
       coefB7= -3.d0/32.d0*((4.d0*co*lcxc-sh).cq.cxw)*lcxc
!
       coefB8= -1.d0/32.d0*((4.d0*co*cxmu-sh).cq.cxw)*cxmu
!
       coefB9= -3.d0/32.d0*((4.d0*co*cxs-sh).cq.cxw)*cxs
!
       coefB10= -3.d0/32.d0*((4.d0*co*cxd-sh).cq.cxw)*cxd
!
       coefB11= -3.d0/32.d0*((4.d0*co*cxu-sh).cq.cxw)*cxu
!
       coefB12= -1.d0/32.d0*((4.d0*co*cxe-sh).cq.cxw)*cxe
!
       cxp(1)= muhs
       cxp(2)= -muhr*rgh
       p2= muhs
!
       IF(ifb.eq.0) THEN
!
        xms(1)= swr
        xms(2)= swi
        xm0= mw
        b0part= VBFNLO_HIGGSTO_b0af_em(scal,cxp,p2,xms,xm0)
        b0sumb= (coefB1.cp.b0part)
        xms(1)= szr
        xms(2)= szi
        xm0= mz
        b0part= VBFNLO_HIGGSTO_b0af_em(scal,cxp,p2,xms,xm0)
        b0sumb= b0sumb+(coefB2.cp.b0part)
        xms(1)= muhs
        xms(2)= -muhr*rgh
        xm0= muhr
        b0part= VBFNLO_HIGGSTO_b0af_em(scal,cxp,p2,xms,xm0)
        b0sumb= b0sumb+(coefB3.cp.b0part)
!
        xms(1:2)= crmbs*co(1:2)
        xm0= sqrt(crmbs)
        b0part= VBFNLO_HIGGSTO_b0af_em(scal,cxp,p2,xms,xm0)
        b0sumf= (coefB5.cp.b0part)
        xms(1:2)= mtl*mtl*co(1:2)
        xm0= mtl
        b0part= VBFNLO_HIGGSTO_b0af_em(scal,cxp,p2,xms,xm0)
        b0sumf= b0sumf+(coefB6.cp.b0part)
        xms(1:2)= crmcs*co(1:2)
        xm0= sqrt(crmcs)
        b0part= VBFNLO_HIGGSTO_b0af_em(scal,cxp,p2,xms,xm0)
        b0sumf= b0sumf+(coefB7.cp.b0part)
        xms(1:2)= mm*mm*co(1:2)
        xm0= mm
        b0part= VBFNLO_HIGGSTO_b0af_em(scal,cxp,p2,xms,xm0)
        b0sumf= b0sumf+(coefB8.cp.b0part)
        xms(1:2)= msq*msq*co(1:2)
        xm0= msq
        b0part= VBFNLO_HIGGSTO_b0af_em(scal,cxp,p2,xms,xm0)
        b0sumf= b0sumf+(coefB9.cp.b0part)
        xms(1:2)= mdq*mdq*co(1:2)
        xm0= mdq
        b0part= VBFNLO_HIGGSTO_b0af_em(scal,cxp,p2,xms,xm0)
        b0sumf= b0sumf+(coefB10.cp.b0part)
        xms(1:2)= muq*muq*co(1:2)
        xm0= muq
        b0part= VBFNLO_HIGGSTO_b0af_em(scal,cxp,p2,xms,xm0)
        b0sumf= b0sumf+(coefB11.cp.b0part)
        xms(1:2)= me*me*co(1:2)
        xm0= me
        b0part= VBFNLO_HIGGSTO_b0af_em(scal,cxp,p2,xms,xm0)
        b0sumf= b0sumf+(coefB12.cp.b0part)
!
        xms(1:2)= cpt(1:2)
        xm0= mt
        b0part= VBFNLO_HIGGSTO_b0af_em(scal,cxp,p2,xms,xm0)
        b0sumt= (coefB4.cp.b0part)
!
        iz= 0
        xms(1:2)= crmbs*co(1:2)
        xm0= sqrt(crmbs)
        nloqcd= crmbs*VBFNLO_HIGGSTO_quarkQCD(scal,cxp,p2,xms,xm0,iz)
        xms(1:2)= crmcs*co(1:2)
        xm0= sqrt(crmcs)
        nloqcd= nloqcd+crmcs*VBFNLO_HIGGSTO_quarkQCD(scal,cxp,p2,xms,xm0,iz)
        IF(qcdc==0) nloqcd= 0.d0
!
        xms(1:2)= cpt(1:2)
        xm0= mt
        iz= 1
        ttqcd= VBFNLO_HIGGSTO_quarkQCD(scal,cxp,p2,xms,xm0,iz)
        ttqcd= cpt.cp.ttqcd
        IF(qcdc==0) ttqcd= 0.d0
!
       ELSEIF(ifb.eq.1) THEN
!
        xms(1:2)= crmbs*co(1:2)
        xm0= sqrt(crmbs)
        b0part= VBFNLO_HIGGSTO_b0af_em(scal,cxp,p2,xms,xm0)
        b0sumf= (coefB5.cp.b0part)
        xms(1:2)= mtl*mtl*co(1:2)
        xm0= mtl
        b0part= VBFNLO_HIGGSTO_b0af_em(scal,cxp,p2,xms,xm0)
        b0sumf= b0sumf+(coefB6.cp.b0part)
        xms(1:2)= crmcs*co(1:2)
        xm0= sqrt(crmcs)
        b0part= VBFNLO_HIGGSTO_b0af_em(scal,cxp,p2,xms,xm0)
        b0sumf= b0sumf+(coefB7.cp.b0part)
        xms(1:2)= mm*mm*co(1:2)
        xm0= mm
        b0part= VBFNLO_HIGGSTO_b0af_em(scal,cxp,p2,xms,xm0)
        b0sumf= b0sumf+(coefB8.cp.b0part)
        xms(1:2)= msq*msq*co(1:2)
        xm0= msq
        b0part= VBFNLO_HIGGSTO_b0af_em(scal,cxp,p2,xms,xm0)
        b0sumf= b0sumf+(coefB9.cp.b0part)
        xms(1:2)= mdq*mdq*co(1:2)
        xm0= mdq
        b0part= VBFNLO_HIGGSTO_b0af_em(scal,cxp,p2,xms,xm0)
        b0sumf= b0sumf+(coefB10.cp.b0part)
        xms(1:2)= muq*muq*co(1:2)
        xm0= muq
        b0part=VBFNLO_HIGGSTO_b0af_em(scal,cxp,p2,xms,xm0)
        b0sumf= b0sumf+(coefB11.cp.b0part)
        xms(1:2)= me*me*co(1:2)
        xm0= me
        b0part= VBFNLO_HIGGSTO_b0af_em(scal,cxp,p2,xms,xm0)
        b0sumf= b0sumf+(coefB12.cp.b0part)
!
        iz= 0
        xms(1:2)= crmbs*co(1:2)
        xm0= sqrt(crmbs)
        nloqcd= crmbs*VBFNLO_HIGGSTO_quarkQCD(scal,cxp,p2,xms,xm0,iz)
        xms(1:2)= crmcs*co(1:2)
        xm0= sqrt(crmcs)
        nloqcd= nloqcd+crmcs*VBFNLO_HIGGSTO_quarkQCD(scal,cxp,p2,xms,xm0,iz)
        IF(qcdc==0) nloqcd= 0.d0
!
       ELSEIF(ifb.eq.2) THEN
!
        xms(1:2)= cpt(1:2)
        xm0= mt
        b0part= VBFNLO_HIGGSTO_b0af_em(scal,cxp,p2,xms,xm0)
        b0sumt= (coefB4.cp.b0part)
!    
        xms(1:2)= cpt(1:2)
        xm0= mt
        iz= 1
        ttqcd= VBFNLO_HIGGSTO_quarkQCD(scal,cxp,p2,xms,xm0,iz)
        ttqcd= cpt.cp.ttqcd
        IF(qcdc==0) ttqcd= 0.d0
!
       ELSEIF(ifb.eq.3) THEN
!
        xms(1)= swr
        xms(2)= swi
        xm0= mw
        b0part= VBFNLO_HIGGSTO_b0af_em(scal,cxp,p2,xms,xm0)
        b0sumb= (coefB1.cp.b0part)
        xms(1)= szr
        xms(2)= szi
        xm0= mz
        b0part= VBFNLO_HIGGSTO_b0af_em(scal,cxp,p2,xms,xm0)
        b0sumb= b0sumb+(coefB2.cp.b0part)
        xms(1)= muhs
        xms(2)= -muhr*rgh
        xm0= muhr
        b0part= VBFNLO_HIGGSTO_b0af_em(scal,cxp,p2,xms,xm0)
        b0sumb= b0sumb+(coefB3.cp.b0part)
!
       ENDIF
!
       IF(ifb.eq.0.or.ifb.eq.1) THEN
        ksumf= -(sh.cq.cxw)*(cxe*clxe+cxmu*clxmu+ &
     &        cxtau*clxtau+3.d0*(cxd*clxd+cxs*clxs+lcxb*lclxb+ &
     &        cxu*clxu+lcxc*lclxc))/32.d0+ &
     &        cxwi*(cxes+cxmus+cxtaus+3.d0*( &
     &        cxds+cxss+lcxbs+cxus+lcxcs))/8.d0
       ENDIF
       IF(ifb.eq.0.or.ifb.eq.2) THEN
        ksumt= 3.d0/32.d0*(4.d0*((ccxt.cq.cxw).cp.ccxt) &
     &        -(clt.cp.((sh.cq.cxw).cp.ccxt)))
       ENDIF
!
       IF(ifb.eq.0.or.ifb.eq.3) THEN
!
        ksumb= (-2.d0*((2.d0*co+cctsi+3.d0*(sh.cq.cxw)).cp.sh) &
     &  -(clcts.cp.((6.d0*cctsi-(sh.cq.cxw)).cp.sh)) &
     &  +3.d0*(clw.cp.((4.d0*co+2.d0*cctsi-(sh.cq.cxw)).cp.sh)) &
     &  -3.d0*(clh.cp.(shs.cq.cxw)))/128.d0
!
       ENDIF
!
       EWC= 4.d0*sqrt(2.d0)*g_f/pis
       sconv= muhs/scals
!
       IF(ifb.eq.0.or.ifb.eq.1) THEN
        totalf= b0sumf+ksumf
       ENDIF
!
       IF(ifb.eq.0.or.ifb.eq.2) THEN
        totalt= b0sumt+ksumt
!
        neg= EWC*(swr*totalt(2)+swi*totalt(1))/sconv &
     &      -as_NLO*g_f/(sqrt(2.d0)*pis)* &
     &      (sh(1)*ttqcd(2)+sh(2)*ttqcd(1))/sconv
!
        IF(neg < 0.d0) THEN
         totalt= 0.d0
         qcdtop= 0.d0
         pcnt= -1
        ELSE
         qcdtop= as_NLO*g_f/(sqrt(2.d0)*pis)* &
     &          (sh(1)*ttqcd(2)+sh(2)*ttqcd(1))/sconv
         pcnt= 1
        ENDIF       
!
       ENDIF
!
       IF(ifb.eq.0.or.ifb.eq.3) THEN
        totalb= b0sumb+ksumb
        IF((swr*totalb(2)+swi*totalb(1)) < 0.d0) THEN
         totalb= 0.d0
         pcnt= -1
        ELSE
         pcnt= 1
        ENDIF
       ENDIF 
!
       total= totalf+totalt+totalb
!
!--- W self energies
!
       deltag= 6.d0*co+0.5d0*(((7.d0*co-4.d0*csts).cq.csts).cp.clcts)
!
       sww0= -(38.d0*cxw+6.d0*wccxt+7.d0*sh &
     & -48.d0*(((wccxt.cq.sh).cq.cxw).cp.wccxt)+8.d0*(cxw.cq.sh))/128.d0 &
     & -3.d0/64.d0*((cxw-sh+(cxws.cq.cxhw)).cp.clh) &
     & +3.d0/32.d0*(((co-4.d0*((wccxt.cq.sh).cq.cxw).cp.wccxt)).cp.wclt) &
     & +((((8.d0*co-17.d0*cstsi+3.d0*cctsi).cp.cxw) &
     & -6.d0*((cxw.cq.sh).cq.cctq)).cp.clcts)/64.d0 &
     & -((cxw.cq.sh).cq.cctq)/32.d0+5.d0/128.d0*(cxw.cq.ccts)
!
       coefW1= -(((8.d0*co-(sh.cq.cxw)).cp.sh)*sh &
     &       -4.d0*((-12.d0*cxw+7.d0*sh).cp.cxw))/192.d0
!
       coefW2= -((cxws.cq.csmcts).cp.(416.d0*co-192.d0*csts  &
     &       -((132.d0*co-((12.d0*co+cctsi).cq.ccts)).cq.ccts)))/192.d0
!
       cxp(1)= swr
       cxp(2)= swi
       p2= mw*mw
!
       axms(1,1)= swr
       axms(1,2)= swi
       axm0(1)= mw
       axms(2,1)= muhs
       axms(2,2)= -muhr*rgh
       axm0(2)= muhr
       b0part= VBFNLO_HIGGSTO_b021_dm_cp(scal,cxp,p2,axms,axm0)
       b0sumw1= (coefW1.cp.b0part)
       b0sumw= (coefW1.cp.b0part)
!
       axms(1,1)= szr
       axms(1,2)= szi
       axm0(1)= mz
       axms(2,1)= swr
       axms(2,2)= swi
       axm0(2)= mw
       b0part= VBFNLO_HIGGSTO_b021_dm_cp(scal,cxp,p2,axms,axm0)
       b0sumw2= (coefW2.cp.b0part)
       b0sumw= b0sumw+(coefW2.cp.b0part)
!
       ksumw= -12.d0*((cxw &
     &    -0.5d0*((3.d0*co-(wccxts.cq.cxws)).cp.wccxt)).cp.wcltmw) &
     &    -((24.d0*cxw-((14.d0*co-(sh.cq.cxw)).cp.sh)).cp.clh) &
     &    +((36.d0*cxw-14.d0*sh-18.d0*((co-4.d0*(wccxt.cq.sh)).cp.wccxt) &
     &    +(shs.cq.cxw)).cp.clw) &
     &    -6.d0*(((2.d0*co+((cxwi-12.d0*shi).cp.wccxt)).cp.wccxt) &
     &    +1.d0/6.d0*((15.d0*co-(sh.cq.cxw)).cp.sh) &
     &    +2.d0/9.d0*((97.d0*co+9.d0*(cxw.cq.sh)).cp.cxw)) &
     &    +(((cxw.cq.ccts).cp.(co-6.d0*(cxw.cq.sh))).cq.ccts) &
     &    -2.d0*(((cxw.cq.csmcts).cp.clcts).cp.(62.d0*co &
     &    -48.d0*csts-5.d0*cctsi)) &
     &    -18.d0*(((cxws.cq.sh).cq.cctq).cp.clcts) &
     &    -72.d0*((wclt.cp.wccxts).cp.(shi-1.d0/12.d0*(wccxt.cq.cxws))) &
     &    +23.d0*(cxw.cq.ccts)
!
       ksumw= ksumw/192.d0+3.d0/16.d0*(cxw.cp.(clw-clmw))
!
       sww= b0sumw+ksumw
!
       DW= (-sww+sww0)/sconv+deltag/16.d0
!       DW= 0.d0
!
       IF(ifb.eq.0) THEN
        VBFNLO_HIGGSTO_SHH= rgh/muhr* &
     &       (1.d0+EWC*(swr*DW(1)-swi*DW(2))) &
     &       -EWC*(swr*total(2)+swi*total(1))/sconv &
     &       +as_NLO*g_f/(sqrt(2.d0)*pis)* &
     &       (sh(1)*nloqcd(2)+sh(2)*nloqcd(1))/sconv &
     &       +qcdtop
!
       ELSEIF(ifb.eq.1) THEN
        VBFNLO_HIGGSTO_SHH= rgh/muhr* &
     &       (1.d0+EWC*(swr*DW(1)-swi*DW(2))) &
     &       -EWC*(swr*totalf(2)+swi*totalf(1))/sconv &
     &       +as_NLO*g_f/(sqrt(2.d0)*pis)* &
     &       (sh(1)*nloqcd(2)+sh(2)*nloqcd(1))/sconv
       ELSEIF(ifb.eq.2) THEN
        VBFNLO_HIGGSTO_SHH= rgh/muhr* &
     &       (1.d0+EWC*(swr*DW(1)-swi*DW(2))) &
     &       -EWC*(swr*totalt(2)+swi*totalt(1))/sconv &
     &       +qcdtop
       ELSEIF(ifb.eq.3) THEN
        VBFNLO_HIGGSTO_SHH= rgh/muhr* &
     &       (1.d0+EWC*(swr*DW(1)-swi*DW(2))) &
     &       -EWC*(swr*totalb(2)+swi*totalb(1))/sconv
       ENDIF
!
       RETURN
!
       END FUNCTION VBFNLO_HIGGSTO_SHH
!
      !END MODULE VBFNLO_HIGGSTO_HBB_cp
!
!--------------------------------------------------------------------------
!
!----------------------------------------------------------------
! H total grids ``Handbook of LHC Higgs Cross Sections: 1. Inclusive Observables,''
!----------------------------------------------------------------
!
      SUBROUTINE VBFNLO_HIGGSTO_gridHt(mass,evalue)
! 
      IMPLICIT NONE
!
      INTEGER top,gdim
      REAL*8 u,value,evalue,mass
      REAL*8, dimension(321) :: bc,cc,dc
! 
! u value of M_H at which the spline is to be evaluated
!
      gdim= 321
!
      CALL VBFNLO_HIGGSTO_FMMsplineSingleHt(bc,cc,dc,top,gdim)
!
      u= mass
      CALL VBFNLO_HIGGSTO_Seval3SingleHt(u,bc,cc,dc,top,gdim,value)
!
      evalue= value
!
      RETURN
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
      SUBROUTINE VBFNLO_HIGGSTO_FMMsplineSingleHt(b,c,d,top,gdim)
!
!---------------------------------------------------------------------------
!
      INTEGER k,n,i,top,gdim,l
!
      REAL*8, dimension(321) :: xc,yc
      REAL*8, dimension(321) :: x,y
!
      REAL*8, DIMENSION(gdim) :: b 
! linear coeff
!
      REAL*8, DIMENSION(gdim) :: c 
! quadratic coeff.
!
      REAL*8, DIMENSION(gdim) :: d 
! cubic coeff.
!
      REAL*8 :: t
      REAL*8,PARAMETER:: ZERO=0.0, TWO=2.0, THREE=3.0
!
! The grid
!
!
      DATA (xc(i),i=1,321)/ &
     &   80.0d0, 81.0d0, 82.0d0, 83.0d0, 84.0d0, 85.0d0, &
     &   86.0d0, 87.0d0, 88.0d0, 89.0d0, 90.0d0, 91.0d0, 92.0d0, &
     &   93.0d0, 94.0d0, 95.0d0, 96.0d0, 97.0d0, 98.0d0, 99.0d0, &
     &   100.0d0, 101.0d0, 102.0d0, 103.0d0, 104.0d0, 105.0d0, &
     &   106.0d0, 107.0d0, 108.0d0, 109.0d0, 110.0d0, 110.5d0, &
     &   111.0d0, 111.5d0, 112.0d0, 112.5d0, 113.0d0, 113.5d0, &
     &   114.0d0, 114.5d0, 115.0d0, 115.5d0, 116.0d0, 116.5d0, &
     &   117.0d0, 117.5d0, 118.0d0, 118.5d0, 119.0d0, 119.5d0, &
     &   120.0d0, 120.1d0, 120.2d0, 120.3d0, 120.4d0, 120.5d0, &
     &   120.6d0, 120.7d0, 120.8d0, 120.9d0, 121.0d0, 121.1d0, &
     &   121.2d0, 121.3d0, 121.4d0, 121.5d0, 121.6d0, 121.7d0, &
     &   121.8d0, 121.9d0, 122.0d0, 122.1d0, 122.2d0, 122.3d0, &
     &   122.4d0, 122.5d0, 122.6d0, 122.7d0, 122.8d0, 122.9d0, &
     &   123.0d0, 123.1d0, 123.2d0, 123.3d0, 123.4d0, 123.5d0, &
     &   123.6d0, 123.7d0, 123.8d0, 123.9d0, 124.0d0, 124.1d0, &
     &   124.2d0, 124.3d0, 124.4d0, 124.5d0, 124.6d0, 124.7d0, &
     &   124.8d0, 124.9d0, 125.0d0, 125.1d0, 125.2d0, 125.3d0, &
     &   125.4d0, 125.5d0, 125.6d0, 125.7d0, 125.8d0, 125.9d0, &
     &   126.0d0, 126.1d0, 126.2d0, 126.3d0, 126.4d0, 126.5d0, &
     &   126.6d0, 126.7d0, 126.8d0, 126.9d0, 127.0d0, 127.1d0, &
     &   127.2d0, 127.3d0, 127.4d0, 127.5d0, 127.6d0, 127.7d0, &
     &   127.8d0, 127.9d0, 128.0d0, 128.1d0, 128.2d0, 128.3d0, &
     &   128.4d0, 128.5d0, 128.6d0, 128.7d0, 128.8d0, 128.9d0, &
     &   129.0d0, 129.1d0, 129.2d0, 129.3d0, 129.4d0, 129.5d0, &
     &   129.6d0, 129.7d0, 129.8d0, 129.9d0, 130.0d0, 130.5d0, &
     &   131.0d0, 131.5d0, 132.0d0, 132.5d0, 133.0d0, 133.5d0, &
     &   134.0d0, 134.5d0, 135.0d0, 135.5d0, 136.0d0, 136.5d0, &
     &   137.0d0, 137.5d0, 138.0d0, 138.5d0, 139.0d0, 139.5d0, &
     &   140.0d0, 140.5d0, 141.0d0, 141.5d0, 142.0d0, 142.5d0, &
     &   143.0d0, 143.5d0, 144.0d0, 144.5d0, 145.0d0, 145.5d0, &
     &   146.0d0, 146.5d0, 147.0d0, 147.5d0, 148.0d0, 148.5d0, &
     &   149.0d0, 149.5d0, 150.0d0, 152.0d0, 154.0d0, 156.0d0, &
     &   158.0d0, 160.0d0, 162.0d0, 164.0d0, 165.0d0, 166.0d0, &
     &   168.0d0, 170.0d0, 172.0d0, 174.0d0, 175.0d0, 176.0d0, &
     &   178.0d0, 180.0d0, 182.0d0, 184.0d0, 185.0d0, 186.0d0, &
     &   188.0d0, 190.0d0, 192.0d0, 194.0d0, 195.0d0, 196.0d0, &
     &   198.0d0, 200.0d0, 202.0d0, 204.0d0, 206.0d0, 208.0d0, &
     &   210.0d0, 212.0d0, 214.0d0, 216.0d0, 218.0d0, 220.0d0, &
     &   222.0d0, 224.0d0, 226.0d0, 228.0d0, 230.0d0, 232.0d0, &
     &   234.0d0, 236.0d0, 238.0d0, 240.0d0, 242.0d0, 244.0d0, &
     &   246.0d0, 248.0d0, 250.0d0, 252.0d0, 254.0d0, 256.0d0, &
     &   258.0d0, 260.0d0, 262.0d0, 264.0d0, 266.0d0, 268.0d0, &
     &   270.0d0, 272.0d0, 274.0d0, 276.0d0, 278.0d0, 280.0d0, &
     &   282.0d0, 284.0d0, 286.0d0, 288.0d0, 290.0d0, 292.0d0, &
     &   294.0d0, 296.0d0, 298.0d0, 300.0d0, 305.0d0, 310.0d0, &
     &   315.0d0, 320.0d0, 325.0d0, 330.0d0, 335.0d0, 340.0d0, &
     &   345.0d0, 350.0d0, 360.0d0, 370.0d0, 380.0d0, 390.0d0, &
     &   400.0d0, 420.0d0, 440.0d0, 450.0d0, 460.0d0, 480.0d0, &
     &   500.0d0, 520.0d0, 540.0d0, 550.0d0, 560.0d0, 580.0d0, &
     &   600.0d0, 620.0d0, 640.0d0, 650.0d0, 660.0d0, 680.0d0, &
     &   700.0d0, 720.0d0, 740.0d0, 750.0d0, 760.0d0, 780.0d0, &
     &   800.0d0, 820.0d0, 840.0d0, 850.0d0, 860.0d0, 880.0d0, &
     &   900.0d0, 920.0d0, 940.0d0, 950.0d0, 960.0d0, 980.0d0, &
     &   1000.0d0/

!
      DATA (yc(i),i=1,321)/ &
     & 1.99d-03, 2.01d-03, 2.04d-03, 2.06d-03, 2.08d-03, 2.11d-03, &
     & 2.13d-03, 2.15d-03, 2.18d-03, 2.20d-03, 2.22d-03, 2.25d-03, &
     & 2.27d-03, 2.30d-03, 2.32d-03, 2.35d-03, 2.37d-03, 2.40d-03, &
     & 2.43d-03, 2.46d-03, 2.48d-03, 2.51d-03, 2.54d-03, 2.58d-03, &
     & 2.61d-03, 2.64d-03, 2.68d-03, 2.72d-03, 2.76d-03, 2.80d-03, &
     & 2.85d-03, 2.87d-03, 2.89d-03, 2.92d-03, 2.95d-03, 2.97d-03, &
     & 3.00d-03, 3.03d-03, 3.06d-03, 3.09d-03, 3.12d-03, 3.15d-03, &
     & 3.19d-03, 3.22d-03, 3.26d-03, 3.30d-03, 3.33d-03, 3.38d-03, &
     & 3.42d-03, 3.46d-03, 3.51d-03, 3.52d-03, 3.52d-03, 3.53d-03, &
     & 3.54d-03, 3.55d-03, 3.56d-03, 3.57d-03, 3.58d-03, 3.59d-03, &
     & 3.60d-03, 3.61d-03, 3.62d-03, 3.63d-03, 3.64d-03, 3.65d-03, &
     & 3.66d-03, 3.67d-03, 3.68d-03, 3.69d-03, 3.71d-03, 3.72d-03, &
     & 3.73d-03, 3.74d-03, 3.75d-03, 3.76d-03, 3.77d-03, 3.78d-03, &
     & 3.79d-03, 3.81d-03, 3.82d-03, 3.83d-03, 3.84d-03, 3.85d-03, &
     & 3.86d-03, 3.88d-03, 3.89d-03, 3.90d-03, 3.91d-03, 3.93d-03, &
     & 3.94d-03, 3.95d-03, 3.96d-03, 3.98d-03, 3.99d-03, 4.00d-03, &
     & 4.02d-03, 4.03d-03, 4.04d-03, 4.06d-03, 4.07d-03, 4.08d-03, &
     & 4.10d-03, 4.11d-03, 4.12d-03, 4.14d-03, 4.15d-03, 4.17d-03, &
     & 4.18d-03, 4.20d-03, 4.21d-03, 4.23d-03, 4.24d-03, 4.26d-03, &
     & 4.27d-03, 4.29d-03, 4.30d-03, 4.32d-03, 4.33d-03, 4.35d-03, &
     & 4.36d-03, 4.38d-03, 4.40d-03, 4.41d-03, 4.43d-03, 4.45d-03, &
     & 4.46d-03, 4.48d-03, 4.50d-03, 4.51d-03, 4.53d-03, 4.55d-03, &
     & 4.57d-03, 4.58d-03, 4.60d-03, 4.62d-03, 4.64d-03, 4.65d-03, &
     & 4.67d-03, 4.69d-03, 4.71d-03, 4.73d-03, 4.75d-03, 4.77d-03, &
     & 4.79d-03, 4.81d-03, 4.83d-03, 4.85d-03, 4.87d-03, 4.89d-03, &
     & 4.91d-03, 5.01d-03, 5.12d-03, 5.23d-03, 5.35d-03, 5.48d-03, &
     & 5.60d-03, 5.74d-03, 5.88d-03, 6.03d-03, 6.18d-03, 6.34d-03, &
     & 6.51d-03, 6.69d-03, 6.87d-03, 7.06d-03, 7.26d-03, 7.47d-03, &
     & 7.70d-03, 7.93d-03, 8.17d-03, 8.43d-03, 8.70d-03, 8.98d-03, &
     & 9.28d-03, 9.59d-03, 9.92d-03, 1.03d-02, 1.06d-02, 1.10d-02, &
     & 1.14d-02, 1.19d-02, 1.23d-02, 1.28d-02, 1.33d-02, 1.39d-02, &
     & 1.45d-02, 1.51d-02, 1.58d-02, 1.65d-02, 1.73d-02, 2.11d-02, &
     & 2.66d-02, 3.51d-02, 5.02d-02, 8.31d-02, 1.47d-01, 2.15d-01, &
     & 2.46d-01, 2.76d-01, 3.30d-01, 3.80d-01, 4.29d-01, 4.77d-01, &
     & 5.01d-01, 5.25d-01, 5.75d-01, 6.31d-01, 7.00d-01, 7.88d-01, &
     & 8.32d-01, 8.76d-01, 9.60d-01, 1.04d+00, 1.12d+00, 1.20d+00, &
     & 1.24d+00, 1.28d+00, 1.35d+00, 1.43d+00, 1.51d+00, 1.59d+00, &
     & 1.68d+00, 1.76d+00, 1.85d+00, 1.93d+00, 2.02d+00, 2.12d+00, &
     & 2.21d+00, 2.31d+00, 2.40d+00, 2.50d+00, 2.61d+00, 2.71d+00, &
     & 2.82d+00, 2.93d+00, 3.04d+00, 3.16d+00, 3.27d+00, 3.40d+00, &
     & 3.52d+00, 3.64d+00, 3.77d+00, 3.91d+00, 4.04d+00, 4.18d+00, &
     & 4.32d+00, 4.46d+00, 4.61d+00, 4.76d+00, 4.91d+00, 5.07d+00, &
     & 5.23d+00, 5.39d+00, 5.55d+00, 5.72d+00, 5.89d+00, 6.07d+00, &
     & 6.25d+00, 6.43d+00, 6.61d+00, 6.80d+00, 6.99d+00, 7.19d+00, &
     & 7.39d+00, 7.59d+00, 7.79d+00, 8.00d+00, 8.22d+00, 8.43d+00, &
     & 8.99d+00, 9.57d+00, 1.02d+01, 1.08d+01, 1.14d+01, 1.21d+01, &
     & 1.28d+01, 1.35d+01, 1.42d+01, 1.52d+01, 1.76d+01, 2.02d+01, &
     & 2.31d+01, 2.61d+01, 2.92d+01, 3.59d+01, 4.30d+01, 4.68d+01, &
     & 5.08d+01, 5.91d+01, 6.80d+01, 7.75d+01, 8.77d+01, 9.30d+01, &
     & 9.86d+01, 1.10d+02, 1.23d+02, 1.36d+02, 1.50d+02, 1.58d+02, &
     & 1.65d+02, 1.82d+02, 1.99d+02, 2.17d+02, 2.37d+02, 2.47d+02, &
     & 2.58d+02, 2.80d+02, 3.04d+02, 3.30d+02, 3.57d+02, 3.71d+02, &
     & 3.86d+02, 4.16d+02, 4.49d+02, 4.84d+02, 5.21d+02, 5.40d+02, &
     & 5.60d+02, 6.02d+02, 6.47d+02/

!
      n= 321
      FORALL(l=1:321)
       x(l)= xc(l)
       y(l)= yc(l)
      ENDFORALL

!.....Set up tridiagonal system.........................................
!     b=diagonal, d=offdiagonal, c=right-hand side
!
      d(1)= x(2)-x(1)
      c(2)= (y(2)-y(1))/d(1)
      DO k= 2,n-1
       d(k)= x(k+1)-x(k)
       b(k)= TWO*(d(k-1)+d(k))
       c(k+1)= (y(k+1)-y(k))/d(k)
       c(k)= c(k+1)-c(k)
      END DO
!
!.....End conditions.  third derivatives at x(1) and x(n) obtained
!     from divided differences.......................................
!
      b(1)= -d(1)
      b(n)= -d(n-1)
      c(1)= ZERO
      c(n)= ZERO
      IF (n > 3) THEN
       c(1)= c(3)/(x(4)-x(2))-c(2)/(x(3)-x(1))
       c(n)= c(n-1)/(x(n)-x(n-2))-c(n-2)/(x(n-1)-x(n-3))
       c(1)= c(1)*d(1)*d(1)/(x(4)-x(1))
       c(n)= -c(n)*d(n-1)*d(n-1)/(x(n)-x(n-3))
      END IF
!
      DO k=2,n    ! forward elimination
       t= d(k-1)/b(k-1)
       b(k)= b(k)-t*d(k-1)
       c(k)= c(k)-t*c(k-1)
      END DO
!
      c(n)= c(n)/b(n)   
!
! back substitution ( makes c the sigma of text)
!
      DO k=n-1,1,-1
       c(k)= (c(k)-d(k)*c(k+1))/b(k)
      END DO
!
!.....Compute polynomial coefficients...................................
!
      b(n)= (y(n)-y(n-1))/d(n-1)+d(n-1)*(c(n-1)+c(n)+c(n))
      DO k=1,n-1
       b(k)= (y(k+1)-y(k))/d(k)-d(k)*(c(k+1)+c(k)+c(k))
       d(k)= (c(k+1)-c(k))/d(k)
       c(k)= THREE*c(k)
      END DO
      c(n)= THREE*c(n)
      d(n)= d(n-1)
!
      RETURN
!
      END SUBROUTINE VBFNLO_HIGGSTO_FMMsplineSingleHt   
!
!------------------------------------------------------------------------
!
      SUBROUTINE VBFNLO_HIGGSTO_Seval3SingleHt(u,b,c,d,top,gdim,f,fp,fpp,fppp)
!
! ---------------------------------------------------------------------------
!
      REAL*8,INTENT(IN) :: u 
! abscissa at which the spline is to be evaluated
!
      INTEGER j,k,n,l,top,gdim
!
      REAL*8, dimension(321) :: xc,yc
      REAL*8, dimension(321) :: x,y
      REAL*8, DIMENSION(gdim) :: b,c,d 
! linear,quadratic,cubic coeff
!
      REAL*8,INTENT(OUT),OPTIONAL:: f,fp,fpp,fppp 
! function, 1st,2nd,3rd deriv
!
      INTEGER, SAVE :: i=1
      REAL*8    :: dx
      REAL*8,PARAMETER:: TWO=2.0, THREE=3.0, SIX=6.0
!
! The grid
!
      DATA (xc(i),i=1,321)/ &
     &   80.0d0, 81.0d0, 82.0d0, 83.0d0, 84.0d0, 85.0d0, &
     &   86.0d0, 87.0d0, 88.0d0, 89.0d0, 90.0d0, 91.0d0, 92.0d0, &
     &   93.0d0, 94.0d0, 95.0d0, 96.0d0, 97.0d0, 98.0d0, 99.0d0, &
     &   100.0d0, 101.0d0, 102.0d0, 103.0d0, 104.0d0, 105.0d0, &
     &   106.0d0, 107.0d0, 108.0d0, 109.0d0, 110.0d0, 110.5d0, &
     &   111.0d0, 111.5d0, 112.0d0, 112.5d0, 113.0d0, 113.5d0, &
     &   114.0d0, 114.5d0, 115.0d0, 115.5d0, 116.0d0, 116.5d0, &
     &   117.0d0, 117.5d0, 118.0d0, 118.5d0, 119.0d0, 119.5d0, &
     &   120.0d0, 120.1d0, 120.2d0, 120.3d0, 120.4d0, 120.5d0, &
     &   120.6d0, 120.7d0, 120.8d0, 120.9d0, 121.0d0, 121.1d0, &
     &   121.2d0, 121.3d0, 121.4d0, 121.5d0, 121.6d0, 121.7d0, &
     &   121.8d0, 121.9d0, 122.0d0, 122.1d0, 122.2d0, 122.3d0, &
     &   122.4d0, 122.5d0, 122.6d0, 122.7d0, 122.8d0, 122.9d0, &
     &   123.0d0, 123.1d0, 123.2d0, 123.3d0, 123.4d0, 123.5d0, &
     &   123.6d0, 123.7d0, 123.8d0, 123.9d0, 124.0d0, 124.1d0, &
     &   124.2d0, 124.3d0, 124.4d0, 124.5d0, 124.6d0, 124.7d0, &
     &   124.8d0, 124.9d0, 125.0d0, 125.1d0, 125.2d0, 125.3d0, &
     &   125.4d0, 125.5d0, 125.6d0, 125.7d0, 125.8d0, 125.9d0, &
     &   126.0d0, 126.1d0, 126.2d0, 126.3d0, 126.4d0, 126.5d0, &
     &   126.6d0, 126.7d0, 126.8d0, 126.9d0, 127.0d0, 127.1d0, &
     &   127.2d0, 127.3d0, 127.4d0, 127.5d0, 127.6d0, 127.7d0, &
     &   127.8d0, 127.9d0, 128.0d0, 128.1d0, 128.2d0, 128.3d0, &
     &   128.4d0, 128.5d0, 128.6d0, 128.7d0, 128.8d0, 128.9d0, &
     &   129.0d0, 129.1d0, 129.2d0, 129.3d0, 129.4d0, 129.5d0, &
     &   129.6d0, 129.7d0, 129.8d0, 129.9d0, 130.0d0, 130.5d0, &
     &   131.0d0, 131.5d0, 132.0d0, 132.5d0, 133.0d0, 133.5d0, &
     &   134.0d0, 134.5d0, 135.0d0, 135.5d0, 136.0d0, 136.5d0, &
     &   137.0d0, 137.5d0, 138.0d0, 138.5d0, 139.0d0, 139.5d0, &
     &   140.0d0, 140.5d0, 141.0d0, 141.5d0, 142.0d0, 142.5d0, &
     &   143.0d0, 143.5d0, 144.0d0, 144.5d0, 145.0d0, 145.5d0, &
     &   146.0d0, 146.5d0, 147.0d0, 147.5d0, 148.0d0, 148.5d0, &
     &   149.0d0, 149.5d0, 150.0d0, 152.0d0, 154.0d0, 156.0d0, &
     &   158.0d0, 160.0d0, 162.0d0, 164.0d0, 165.0d0, 166.0d0, &
     &   168.0d0, 170.0d0, 172.0d0, 174.0d0, 175.0d0, 176.0d0, &
     &   178.0d0, 180.0d0, 182.0d0, 184.0d0, 185.0d0, 186.0d0, &
     &   188.0d0, 190.0d0, 192.0d0, 194.0d0, 195.0d0, 196.0d0, &
     &   198.0d0, 200.0d0, 202.0d0, 204.0d0, 206.0d0, 208.0d0, &
     &   210.0d0, 212.0d0, 214.0d0, 216.0d0, 218.0d0, 220.0d0, &
     &   222.0d0, 224.0d0, 226.0d0, 228.0d0, 230.0d0, 232.0d0, &
     &   234.0d0, 236.0d0, 238.0d0, 240.0d0, 242.0d0, 244.0d0, &
     &   246.0d0, 248.0d0, 250.0d0, 252.0d0, 254.0d0, 256.0d0, &
     &   258.0d0, 260.0d0, 262.0d0, 264.0d0, 266.0d0, 268.0d0, &
     &   270.0d0, 272.0d0, 274.0d0, 276.0d0, 278.0d0, 280.0d0, &
     &   282.0d0, 284.0d0, 286.0d0, 288.0d0, 290.0d0, 292.0d0, &
     &   294.0d0, 296.0d0, 298.0d0, 300.0d0, 305.0d0, 310.0d0, &
     &   315.0d0, 320.0d0, 325.0d0, 330.0d0, 335.0d0, 340.0d0, &
     &   345.0d0, 350.0d0, 360.0d0, 370.0d0, 380.0d0, 390.0d0, &
     &   400.0d0, 420.0d0, 440.0d0, 450.0d0, 460.0d0, 480.0d0, &
     &   500.0d0, 520.0d0, 540.0d0, 550.0d0, 560.0d0, 580.0d0, &
     &   600.0d0, 620.0d0, 640.0d0, 650.0d0, 660.0d0, 680.0d0, &
     &   700.0d0, 720.0d0, 740.0d0, 750.0d0, 760.0d0, 780.0d0, &
     &   800.0d0, 820.0d0, 840.0d0, 850.0d0, 860.0d0, 880.0d0, &
     &   900.0d0, 920.0d0, 940.0d0, 950.0d0, 960.0d0, 980.0d0, &
     &   1000.0d0/

!
      DATA (yc(i),i=1,321)/ &
     & 1.99d-03, 2.01d-03, 2.04d-03, 2.06d-03, 2.08d-03, 2.11d-03, &
     & 2.13d-03, 2.15d-03, 2.18d-03, 2.20d-03, 2.22d-03, 2.25d-03, &
     & 2.27d-03, 2.30d-03, 2.32d-03, 2.35d-03, 2.37d-03, 2.40d-03, &
     & 2.43d-03, 2.46d-03, 2.48d-03, 2.51d-03, 2.54d-03, 2.58d-03, &
     & 2.61d-03, 2.64d-03, 2.68d-03, 2.72d-03, 2.76d-03, 2.80d-03, &
     & 2.85d-03, 2.87d-03, 2.89d-03, 2.92d-03, 2.95d-03, 2.97d-03, &
     & 3.00d-03, 3.03d-03, 3.06d-03, 3.09d-03, 3.12d-03, 3.15d-03, &
     & 3.19d-03, 3.22d-03, 3.26d-03, 3.30d-03, 3.33d-03, 3.38d-03, &
     & 3.42d-03, 3.46d-03, 3.51d-03, 3.52d-03, 3.52d-03, 3.53d-03, &
     & 3.54d-03, 3.55d-03, 3.56d-03, 3.57d-03, 3.58d-03, 3.59d-03, &
     & 3.60d-03, 3.61d-03, 3.62d-03, 3.63d-03, 3.64d-03, 3.65d-03, &
     & 3.66d-03, 3.67d-03, 3.68d-03, 3.69d-03, 3.71d-03, 3.72d-03, &
     & 3.73d-03, 3.74d-03, 3.75d-03, 3.76d-03, 3.77d-03, 3.78d-03, &
     & 3.79d-03, 3.81d-03, 3.82d-03, 3.83d-03, 3.84d-03, 3.85d-03, &
     & 3.86d-03, 3.88d-03, 3.89d-03, 3.90d-03, 3.91d-03, 3.93d-03, &
     & 3.94d-03, 3.95d-03, 3.96d-03, 3.98d-03, 3.99d-03, 4.00d-03, &
     & 4.02d-03, 4.03d-03, 4.04d-03, 4.06d-03, 4.07d-03, 4.08d-03, &
     & 4.10d-03, 4.11d-03, 4.12d-03, 4.14d-03, 4.15d-03, 4.17d-03, &
     & 4.18d-03, 4.20d-03, 4.21d-03, 4.23d-03, 4.24d-03, 4.26d-03, &
     & 4.27d-03, 4.29d-03, 4.30d-03, 4.32d-03, 4.33d-03, 4.35d-03, &
     & 4.36d-03, 4.38d-03, 4.40d-03, 4.41d-03, 4.43d-03, 4.45d-03, &
     & 4.46d-03, 4.48d-03, 4.50d-03, 4.51d-03, 4.53d-03, 4.55d-03, &
     & 4.57d-03, 4.58d-03, 4.60d-03, 4.62d-03, 4.64d-03, 4.65d-03, &
     & 4.67d-03, 4.69d-03, 4.71d-03, 4.73d-03, 4.75d-03, 4.77d-03, &
     & 4.79d-03, 4.81d-03, 4.83d-03, 4.85d-03, 4.87d-03, 4.89d-03, &
     & 4.91d-03, 5.01d-03, 5.12d-03, 5.23d-03, 5.35d-03, 5.48d-03, &
     & 5.60d-03, 5.74d-03, 5.88d-03, 6.03d-03, 6.18d-03, 6.34d-03, &
     & 6.51d-03, 6.69d-03, 6.87d-03, 7.06d-03, 7.26d-03, 7.47d-03, &
     & 7.70d-03, 7.93d-03, 8.17d-03, 8.43d-03, 8.70d-03, 8.98d-03, &
     & 9.28d-03, 9.59d-03, 9.92d-03, 1.03d-02, 1.06d-02, 1.10d-02, &
     & 1.14d-02, 1.19d-02, 1.23d-02, 1.28d-02, 1.33d-02, 1.39d-02, &
     & 1.45d-02, 1.51d-02, 1.58d-02, 1.65d-02, 1.73d-02, 2.11d-02, &
     & 2.66d-02, 3.51d-02, 5.02d-02, 8.31d-02, 1.47d-01, 2.15d-01, &
     & 2.46d-01, 2.76d-01, 3.30d-01, 3.80d-01, 4.29d-01, 4.77d-01, &
     & 5.01d-01, 5.25d-01, 5.75d-01, 6.31d-01, 7.00d-01, 7.88d-01, &
     & 8.32d-01, 8.76d-01, 9.60d-01, 1.04d+00, 1.12d+00, 1.20d+00, &
     & 1.24d+00, 1.28d+00, 1.35d+00, 1.43d+00, 1.51d+00, 1.59d+00, &
     & 1.68d+00, 1.76d+00, 1.85d+00, 1.93d+00, 2.02d+00, 2.12d+00, &
     & 2.21d+00, 2.31d+00, 2.40d+00, 2.50d+00, 2.61d+00, 2.71d+00, &
     & 2.82d+00, 2.93d+00, 3.04d+00, 3.16d+00, 3.27d+00, 3.40d+00, &
     & 3.52d+00, 3.64d+00, 3.77d+00, 3.91d+00, 4.04d+00, 4.18d+00, &
     & 4.32d+00, 4.46d+00, 4.61d+00, 4.76d+00, 4.91d+00, 5.07d+00, &
     & 5.23d+00, 5.39d+00, 5.55d+00, 5.72d+00, 5.89d+00, 6.07d+00, &
     & 6.25d+00, 6.43d+00, 6.61d+00, 6.80d+00, 6.99d+00, 7.19d+00, &
     & 7.39d+00, 7.59d+00, 7.79d+00, 8.00d+00, 8.22d+00, 8.43d+00, &
     & 8.99d+00, 9.57d+00, 1.02d+01, 1.08d+01, 1.14d+01, 1.21d+01, &
     & 1.28d+01, 1.35d+01, 1.42d+01, 1.52d+01, 1.76d+01, 2.02d+01, &
     & 2.31d+01, 2.61d+01, 2.92d+01, 3.59d+01, 4.30d+01, 4.68d+01, &
     & 5.08d+01, 5.91d+01, 6.80d+01, 7.75d+01, 8.77d+01, 9.30d+01, &
     & 9.86d+01, 1.10d+02, 1.23d+02, 1.36d+02, 1.50d+02, 1.58d+02, &
     & 1.65d+02, 1.82d+02, 1.99d+02, 2.17d+02, 2.37d+02, 2.47d+02, &
     & 2.58d+02, 2.80d+02, 3.04d+02, 3.30d+02, 3.57d+02, 3.71d+02, &
     & 3.86d+02, 4.16d+02, 4.49d+02, 4.84d+02, 5.21d+02, 5.40d+02, &
     & 5.60d+02, 6.02d+02, 6.47d+02/
!
      n= 321
      FORALL(l=1:321)
       x(l)= xc(l)
       y(l)= yc(l)
      ENDFORALL
!
!.....First check if u is in the same interval found on the
!     last call to Seval.............................................
!
      IF (  (i<1) .OR. (i >= n) ) i=1
      IF ( (u < x(i))  .OR.  (u >= x(i+1)) ) THEN
       i=1   
!
! binary search
!
       j= n+1
       DO
        k= (i+j)/2
        IF (u < x(k)) THEN
         j= k
        ELSE
         i= k
        ENDIF
        IF (j <= i+1) EXIT
       ENDDO
      ENDIF
!
      dx= u-x(i)   
!
! evaluate the spline
!
      IF (Present(f))    f= y(i)+dx*(b(i)+dx*(c(i)+dx*d(i)))
      IF (Present(fp))   fp= b(i)+dx*(TWO*c(i) + dx*THREE*d(i))
      IF (Present(fpp))  fpp= TWO*c(i) + dx*SIX*d(i)
      IF (Present(fppp)) fppp= SIX*d(i)
!
      RETURN
!
      END SUBROUTINE VBFNLO_HIGGSTO_Seval3SingleHt  
!
      END SUBROUTINE VBFNLO_HIGGSTO_gridHt
!
!----------------------------------------------------------------
! gH grids
!----------------------------------------------------------------
!
      SUBROUTINE VBFNLO_HIGGSTO_gridlow(mass,evalue)
! 
      IMPLICIT NONE
!
      INTEGER top,gdim
      REAL*8 u,value,evalue,mass
      REAL*8, dimension(22) :: bc,cc,dc
! 
! u value of M_H at which the spline is to be evaluated
! top= -1,0,1 lower, central, upper value for m_top
!
      gdim= 22
!
      CALL VBFNLO_HIGGSTO_FMMsplineSingleL(bc,cc,dc,top,gdim)
!
      u= mass
      CALL VBFNLO_HIGGSTO_Seval3SingleL(u,bc,cc,dc,top,gdim,value)
!
      evalue= value
!
      RETURN
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
      SUBROUTINE VBFNLO_HIGGSTO_FMMsplineSingleL(b,c,d,top,gdim)
!
!---------------------------------------------------------------------------
!
      INTEGER k,n,i,top,gdim,l
!
      REAL*8, dimension(22) :: xc,yc
      REAL*8, dimension(22) :: x,y
!
      REAL*8, DIMENSION(gdim) :: b 
! linear coeff
!
      REAL*8, DIMENSION(gdim) :: c 
! quadratic coeff.
!
      REAL*8, DIMENSION(gdim) :: d 
! cubic coeff.
!
      REAL*8 :: t
      REAL*8,PARAMETER:: ZERO=0.0, TWO=2.0, THREE=3.0
!
! The grid
!
!
      DATA (xc(i),i=1,22)/ &
     &   190.d0,191.d0,192.d0,193.d0,194.d0,195.d0,196.d0,197.d0,198.d0, &
     &   199.d0,200.d0,240.d0,241.d0,242.d0,243.d0,244.d0,245.d0,246.d0, &
     &   247.d0,248.d0,249.d0,250.d0/
!
      DATA (yc(i),i=1,22)/ &
     &   0.10346d1,0.10750d1,0.11151d1,0.11549d1,0.11942d1,0.12329d1, &
     &   0.12708d1,0.13083d1,0.13456d1,0.13832d1,0.14212d1, &
     &   0.32401d1,0.32994d1,0.33593d1,0.34200d1,0.34813d1,0.35433d1, &
     &   0.36061d1,0.36695d1,0.37336d1,0.37984d1,0.38640d1/
!
      n= 22
      FORALL(l=1:22)
       x(l)= xc(l)
       y(l)= yc(l)
      ENDFORALL

!.....Set up tridiagonal system.........................................
!     b=diagonal, d=offdiagonal, c=right-hand side
!
      d(1)= x(2)-x(1)
      c(2)= (y(2)-y(1))/d(1)
      DO k= 2,n-1
       d(k)= x(k+1)-x(k)
       b(k)= TWO*(d(k-1)+d(k))
       c(k+1)= (y(k+1)-y(k))/d(k)
       c(k)= c(k+1)-c(k)
      END DO
!
!.....End conditions.  third derivatives at x(1) and x(n) obtained
!     from divided differences.......................................
!
      b(1)= -d(1)
      b(n)= -d(n-1)
      c(1)= ZERO
      c(n)= ZERO
      IF (n > 3) THEN
       c(1)= c(3)/(x(4)-x(2))-c(2)/(x(3)-x(1))
       c(n)= c(n-1)/(x(n)-x(n-2))-c(n-2)/(x(n-1)-x(n-3))
       c(1)= c(1)*d(1)*d(1)/(x(4)-x(1))
       c(n)= -c(n)*d(n-1)*d(n-1)/(x(n)-x(n-3))
      END IF
!
      DO k=2,n    ! forward elimination
       t= d(k-1)/b(k-1)
       b(k)= b(k)-t*d(k-1)
       c(k)= c(k)-t*c(k-1)
      END DO
!
      c(n)= c(n)/b(n)   
!
! back substitution ( makes c the sigma of text)
!
      DO k=n-1,1,-1
       c(k)= (c(k)-d(k)*c(k+1))/b(k)
      END DO
!
!.....Compute polynnomial coefficients...................................
!
      b(n)= (y(n)-y(n-1))/d(n-1)+d(n-1)*(c(n-1)+c(n)+c(n))
      DO k=1,n-1
       b(k)= (y(k+1)-y(k))/d(k)-d(k)*(c(k+1)+c(k)+c(k))
       d(k)= (c(k+1)-c(k))/d(k)
       c(k)= THREE*c(k)
      END DO
      c(n)= THREE*c(n)
      d(n)= d(n-1)
!
      RETURN
!
      END SUBROUTINE VBFNLO_HIGGSTO_FMMsplineSingleL   
!
!------------------------------------------------------------------------
!
      SUBROUTINE VBFNLO_HIGGSTO_Seval3SingleL(u,b,c,d,top,gdim,f,fp,fpp,fppp)
!
! ---------------------------------------------------------------------------
!
      REAL*8,INTENT(IN) :: u 
! abscissa at which the spline is to be evaluated
!
      INTEGER j,k,n,l,top,gdim
!
      REAL*8, dimension(22) :: xc,yc
      REAL*8, dimension(22) :: x,y
      REAL*8, DIMENSION(gdim) :: b,c,d 
! linear,quadratic,cubic coeff
!
      REAL*8,INTENT(OUT),OPTIONAL:: f,fp,fpp,fppp 
! function, 1st,2nd,3rd deriv
!
      INTEGER, SAVE :: i=1
      REAL*8    :: dx
      REAL*8,PARAMETER:: TWO=2.0, THREE=3.0, SIX=6.0
!
! The grid
!
      DATA (xc(l),l=1,22)/ &
     &   190.d0,191.d0,192.d0,193.d0,194.d0,195.d0,196.d0,197.d0,198.d0, &
     &   199.d0,200.d0,240.d0,241.d0,242.d0,243.d0,244.d0,245.d0,246.d0, &
     &   247.d0,248.d0,249.d0,250.d0/
!
      DATA (yc(l),l=1,22)/ &
     &   0.10346d1,0.10750d1,0.11151d1,0.11549d1,0.11942d1,0.12329d1, &
     &   0.12708d1,0.13083d1,0.13456d1,0.13832d1,0.14212d1, &
     &   0.32401d1,0.32994d1,0.33593d1,0.34200d1,0.34813d1,0.35433d1, &
     &   0.36061d1,0.36695d1,0.37336d1,0.37984d1,0.38640d1/
!
      n= 22
      FORALL(l=1:22)
       x(l)= xc(l)
       y(l)= yc(l)
      ENDFORALL
!
!.....First check if u is in the same interval found on the
!     last call to Seval.............................................
!
      IF (  (i<1) .OR. (i >= n) ) i=1
      IF ( (u < x(i))  .OR.  (u >= x(i+1)) ) THEN
       i=1   
!
! binary search
!
       j= n+1
       DO
        k= (i+j)/2
        IF (u < x(k)) THEN
         j= k
        ELSE
         i= k
        ENDIF
        IF (j <= i+1) EXIT
       ENDDO
      ENDIF
!
      dx= u-x(i)   
!
! evaluate the spline
!
      IF (Present(f))    f= y(i)+dx*(b(i)+dx*(c(i)+dx*d(i)))
      IF (Present(fp))   fp= b(i)+dx*(TWO*c(i) + dx*THREE*d(i))
      IF (Present(fpp))  fpp= TWO*c(i) + dx*SIX*d(i)
      IF (Present(fppp)) fppp= SIX*d(i)
!
      RETURN
!
      END SUBROUTINE VBFNLO_HIGGSTO_Seval3SingleL  
!
      END SUBROUTINE VBFNLO_HIGGSTO_gridlow
!
!---------------------------------------------------------------------------------
!
      FUNCTION VBFNLO_HIGGSTO_deriv(scal,rhm) RESULT(value)
      !USE VBFNLO_HIGGSTO_masses
      !USE VBFNLO_HIGGSTO_acmplx_pro
      !USE VBFNLO_HIGGSTO_acmplx_rat
      !USE VBFNLO_HIGGSTO_cmplx_root
      !USE VBFNLO_HIGGSTO_full_ln
      !USE VBFNLO_HIGGSTO_qcd
      !USE VBFNLO_HIGGSTO_set_phys_const
      !USE VBFNLO_HIGGSTO_units
!
      IMPLICIT NONE
!
      INTEGER iz
      REAL*8 scal,scals,rhm,ps,xw,xz,lw,lz,ls,xe,lxe,xmu,lxmu,xtau, &
     &     lxtau,xu,lxu,xd,lxd,xc,lxc,xs,lxs,xt,lxt,xb,lxb,rmbs,rmcs, &
     &     xes,xmus,xtaus,xus,xds,xcs,xss,xts,xbs, &
     &     xec,xmuc,xtauc,xuc,xdc,xcc,xsc,xtc,xbc, &
     &     bwis,bzis,fw,fz, &
     &     bxeis,bxmuis,bxtauis,bxuis,bxdis,bxcis,bxsis,bxtis,bxbis, &
     &     bxes,bxmus,bxtaus,bxus,bxds,bxcs,bxss,bxts,bxbs,   &
     &     asmur,emc,emb,emt
      REAL*8, dimension(2) :: bs,bw,bz,arg,lbw,lbz,dmz,bzi,bh,lbh,bwi, &
     &     dmw,dmh,sh,bxe,lbxe,bxmu,lbxmu,bxtau,lbxtau,bxu,lbxu,bxd, &
     &     lbxd,bxc,lbxc,bxs,lbxs,bxt,lbxt,bxb,lbxb,runbc,bwic,bzic, &
     &     d2mh,dmwh,dmzh,dmzw,dmww,dmzz, &
     &     bxei,bxmui,bxtaui,bxui,bxdi,bxci,bxsi,bxti,bxbi
      REAL*8, dimension(10,2) :: value 
!
      ! INTERFACE
      !  SUBROUTINE VBFNLO_HIGGSTO_INITALPHAS(asord,FR2,MUR,asmur,emc,emb,emt)
      !  !USE VBFNLO_HIGGSTO_DZpar 
      !  IMPLICIT NONE
      !  INTEGER asord
      !  REAL*8 FR2,MUR,asmur,emc,emb,emt,VBFNLO_HIGGSTO_FINDALPHASR0
      !  EXTERNAL VBFNLO_HIGGSTO_FINDALPHASR0
      !  END SUBROUTINE VBFNLO_HIGGSTO_INITALPHAS
      ! END INTERFACE
! !
      ! INTERFACE
      !  FUNCTION VBFNLO_HIGGSTO_ALPHAS(MUR)
      !  !USE VBFNLO_HIGGSTO_NFFIX  
      !  !USE VBFNLO_HIGGSTO_VARFLV 
      !  !USE VBFNLO_HIGGSTO_FRRAT  
      !  !USE VBFNLO_HIGGSTO_ASINP  
      !  !USE VBFNLO_HIGGSTO_ASFTHR 
      !  IMPLICIT NONE
      !  REAL*8 MUR,VBFNLO_HIGGSTO_ALPHAS
      !  END FUNCTION VBFNLO_HIGGSTO_ALPHAS
      ! END INTERFACE
!
      scals= scal*scal
      ps= rhm*rhm
!
      xw= (mw*mw)/ps
      xz= (mz*mz)/ps
      lw= log(xw)
      lz= log(xz)
      ls= log(ps/scals)
      fw= 1.d0-4.d0*xw*(1.d0-3.d0*xw)
      fz= 1.d0-4.d0*xz*(1.d0-3.d0*xz)
!
      asmur= 0.12018d0
      emc= 1.4d0
      emb= 4.75d0
      emt= mt
      iz= 1
      CALL VBFNLO_HIGGSTO_INITALPHAS(iz,one,mz,asmur,emc,emb,emt)
      als= VBFNLO_HIGGSTO_ALPHAS(scal)
      runbc= VBFNLO_HIGGSTO_run_bc(scal)
      rmbs= runbc(2)*runbc(2)
      rmcs= runbc(1)*runbc(1)         
!
      xe= (me*me)/ps 
      lxe= log(xe) 
      xmu= (mm*mm)/ps 
      lxmu= log(xmu) 
      xtau= (mtl*mtl)/ps 
      lxtau= log(xtau) 
      xu= (muq*muq)/ps 
      lxu= log(xu) 
      xd= (mdq*mdq)/ps 
      lxd= log(xd) 
      xc= rmcs/ps 
      lxc= log(xc) 
      xs= (msq*msq)/ps 
      lxs= log(xs) 
      xt= (mt*mt)/ps 
      lxt= log(xt) 
      xb= rmbs/ps 
      lxb= log(xb) 
      xes= xe*xe
      xmus= xmu*xmu
      xtaus= xtau*xtau
      xus= xu*xu
      xds= xd*xd
      xcs= xc*xc
      xss= xs*xs
      xts= xt*xt
      xbs= xb*xb
      xec= xes*xe
      xmuc= xmus*xmu
      xtauc= xtaus*xtau
      xuc= xus*xu
      xdc= xds*xd
      xcc= xcs*xc
      xsc= xss*xs
      xtc= xts*xt
      xbc= xbs*xb
!
! H
!
      bs(1)= -3.d0
      bs(2)= -eps
      bh= (bs(1)).cr.(bs(2))
      arg= (bh+co).cq.(bh-co)
      IF(arg(2).eq.0.d0) arg(2)= eps
      lbh= arg(1).fln.arg(2)
!
!
! W
!
      bs(1)= 1.d0-4.d0*(mw*mw)/ps
      bs(2)= -eps
      bw= (bs(1)).cr.(bs(2))
      IF(bw(2).eq.0.d0) bw(2)= -eps 
      bwi= co.cq.bw
      IF(bwi(2).eq.0.d0) bwi(2)= eps
      bwis= 1.d0/(1.d0-4.d0*(mw*mw)/ps)
      bwic= bwis*bwi
      arg= (bw+co).cq.(bw-co)
      IF(arg(2).eq.0.d0) arg(2)= eps
      lbw= arg(1).fln.arg(2)
!
! Z
!
      bs(1)= 1.d0-4.d0*(mz*mz)/ps
      bs(2)= -eps
      bz= (bs(1)).cr.(bs(2))
      IF(bz(2).eq.0.d0) bz(2)= -eps 
      bzi= co.cq.bz
      IF(bzi(2).eq.0.d0) bzi(2)= eps
      bzis= 1.d0/(1.d0-4.d0*(mz*mz)/ps)
      bzic= bzis*bzi
      arg= (bz+co).cq.(bz-co)
      IF(arg(2).eq.0.d0) arg(2)= eps
      lbz= arg(1).fln.arg(2)
!
! f
!
      bs(1)= 1.d0-4.d0*(me*me)/ps
      bs(2)= -eps
      bxe= (bs(1)).cr.(bs(2))
      arg= (bxe+co).cq.(bxe-co)
      IF(arg(2).eq.0.d0) arg(2)= eps
      lbxe= arg(1).fln.arg(2)
!
      bs(1)= 1.d0-4.d0*(mm*mm)/ps
      bs(2)= -eps
      bxmu= (bs(1)).cr.(bs(2))
      arg= (bxmu+co).cq.(bxmu-co)
      IF(arg(2).eq.0.d0) arg(2)= eps
      lbxmu= arg(1).fln.arg(2)
!
      bs(1)= 1.d0-4.d0*(mtl*mtl)/ps
      bs(2)= -eps
      bxtau= (bs(1)).cr.(bs(2))
      arg= (bxtau+co).cq.(bxtau-co)
      IF(arg(2).eq.0.d0) arg(2)= eps
      lbxtau= arg(1).fln.arg(2)
!
      bs(1)= 1.d0-4.d0*(muq*muq)/ps
      bs(2)= -eps
      bxu= (bs(1)).cr.(bs(2))
      arg= (bxu+co).cq.(bxu-co)
      IF(arg(2).eq.0.d0) arg(2)= eps
      lbxu= arg(1).fln.arg(2)
!
      bs(1)= 1.d0-4.d0*(mdq*mdq)/ps
      bs(2)= -eps
      bxd= (bs(1)).cr.(bs(2))
      arg= (bxd+co).cq.(bxd-co)
      IF(arg(2).eq.0.d0) arg(2)= eps
      lbxd= arg(1).fln.arg(2)
!
      bs(1)= 1.d0-4.d0*rmcs/ps
      bs(2)= -eps
      bxc= (bs(1)).cr.(bs(2))
      arg= (bxc+co).cq.(bxc-co)
      IF(arg(2).eq.0.d0) arg(2)= eps
      lbxc= arg(1).fln.arg(2)
!
      bs(1)= 1.d0-4.d0*(msq*msq)/ps
      bs(2)= -eps
      bxs= (bs(1)).cr.(bs(2))
      arg= (bxs+co).cq.(bxs-co)
      IF(arg(2).eq.0.d0) arg(2)= eps
      lbxs= arg(1).fln.arg(2)
!
      bs(1)= 1.d0-4.d0*(mt*mt)/ps
      bs(2)= -eps
      bxt= (bs(1)).cr.(bs(2))
      arg= (bxt+co).cq.(bxt-co)
      IF(arg(2).eq.0.d0) arg(2)= eps
      lbxt= arg(1).fln.arg(2)
!
      bs(1)= 1.d0-4.d0*rmbs/ps
      bs(2)= -eps
      bxb= (bs(1)).cr.(bs(2))
      arg= (bxb+co).cq.(bxb-co)
      IF(arg(2).eq.0.d0) arg(2)= eps
      lbxb= arg(1).fln.arg(2)
!
      bxei= co.cq.bxe
      bxmui= co.cq.bxmu
      bxtaui= co.cq.bxtau
      bxui= co.cq.bxu
      bxdi= co.cq.bxd
      bxci= co.cq.bxc
      bxsi= co.cq.bxs
      bxti= co.cq.bxt
      bxbi= co.cq.bxb
      bxeis= 1.d0/(1.d0-4.d0*(me*me)/ps)
      bxmuis= 1.d0/(1.d0-4.d0*(mm*mm)/ps)
      bxtauis= 1.d0/(1.d0-4.d0*(mtl*mtl)/ps)
      bxuis= 1.d0/(1.d0-4.d0*(muq*muq)/ps)
      bxdis= 1.d0/(1.d0-4.d0*(mdq*mdq)/ps)
      bxcis= 1.d0/(1.d0-4.d0*rmcs/ps)
      bxsis= 1.d0/(1.d0-4.d0*(msq*msq)/ps)
      bxtis= 1.d0/(1.d0-4.d0*(mt*mt)/ps)
      bxbis= 1.d0/(1.d0-4.d0*rmbs/ps)
!
      bxes= 1.d0-4.d0*(me*me)/ps
      bxmus= 1.d0-4.d0*(mm*mm)/ps
      bxtaus= 1.d0-4.d0*(mtl*mtl)/ps
      bxus= 1.d0-4.d0*(muq*muq)/ps
      bxds= 1.d0-4.d0*(mdq*mdq)/ps
      bxcs= 1.d0-4.d0*rmcs/ps
      bxss= 1.d0-4.d0*(msq*msq)/ps
      bxts= 1.d0-4.d0*(mt*mt)/ps
      bxbs= 1.d0-4.d0*rmbs/ps
!
      sh= - 9.d0*(lbh.cp.bh) &
     &  - fz*(lbz.cp.bz) &
     &  - 2.d0*fw*(lbw.cp.bw) &
     & +co*( &
     &  - (1.d0-6.d0*xz)*lz &
     &  - 2.d0*(1.d0-6.d0*xw)*lw &
     &  - 6.d0*(1.d0-2.d0*xw)*ls &
     & +2.d0*(9.d0+xw*(-10.d0+24.d0*xw)) &
     &  - 2.d0*xz*(5.d0-12.d0*xz) &
     & +6.d0*xz*ls)
!
      sh= sh/(128.d0*xw)*ps+ps/(32.d0*xw)*( &
     &   - 3.d0*xt*bxts*(bxt.cp.lbxt)-3.d0*xt*lxt*co &
     &   - 3.d0*xb*bxbs*(bxb.cp.lbxb)-3.d0*xb*lxb*co &
     &   - 3.d0*xs*bxss*(bxs.cp.lbxs)-3.d0*xs*lxs*co &
     &   - 3.d0*xc*bxcs*(bxc.cp.lbxc)-3.d0*xc*lxc*co &
     &   - 3.d0*xd*bxds*(bxd.cp.lbxd)-3.d0*xd*lxd*co &
     &   - 3.d0*xu*bxus*(bxu.cp.lbxu)-3.d0*xd*lxu*co &
     &   - xtau*bxtaus*(bxtau.cp.lbxtau)-xtau*lxtau*co &
     &   - xmu*bxmus*(bxmu.cp.lbxmu)-xmu*lxmu*co &
     &   - xe*bxes*(bxe.cp.lbxe)-xe*lxe*co+co*( &
     &   - ls*(3.d0*(xt+xb+xs+xc+xd+xu)+xtau+xmu+xe) &
     &  +6.d0*xt*(1.d0-2.d0*xt) &
     &  +6.d0*xb*(1.d0-2.d0*xb) &
     &  +6.d0*xs*(1.d0-2.d0*xs) &
     &  +6.d0*xc*(1.d0-2.d0*xc) &
     &  +6.d0*xd*(1.d0-2.d0*xd) &
     &  +6.d0*xu*(1.d0-2.d0*xu) &
     &  +2.d0*xtau*(1.d0-2.d0*xtau) &
     &  +2.d0*xmu*(1.d0-2.d0*xmu) &
     &  +2.d0*xe*(1.d0-2.d0*xe)))
!
      dmh= - 18.d0*(bh.cp.lbh) &
     &  - 2.d0*xz*fz*(lbz.cq.bz) &
     &  - 2.d0*(1.d0-2.d0*xz)*(lbz.cp.bz) &
     &  - 4.d0*xw*fw*(lbw.cq.bw) &
     &  - 4.d0*(1.d0-2.d0*xw)*(lbw.cp.bw) &
     & +co*( &
     &  - 2.d0*(1.d0-3.d0*xz)*lz &
     &  - 4.d0*(1.d0-3.d0*xw)*lw &
     &  - 12.d0*(1.d0-xw)*ls &
     & +6.d0*(5.d0-2.d0*xw*(1.d0+2.d0*xw)) &
     &  - 6.d0*xz*(1.d0+2.d0*xz) &
     & +6.d0*xz*ls)
!
      dmh= dmh/(128.d0*xw) &
     & +1.d0/(32.d0*xw)*( &
     &  - 2.d0*xes*bxes*(lbxe.cq.bxe) &
     &  - 2.d0*xmus*bxmus*(lbxmu.cq.bxmu) &
     &  - 2.d0*xtaus*bxtaus*(lbxtau.cq.bxtau) &
     &  - 6.d0*xus*bxus*(lbxu.cq.bxu) &
     &  - 6.d0*xds*bxds*(lbxd.cq.bxd) &
     &  - 6.d0*xcs*bxcs*(lbxc.cq.bxc) &
     &  - 6.d0*xss*bxss*(lbxs.cq.bxs) &
     &  - 6.d0*xts*bxts*(lbxt.cq.bxt) &
     &  - 6.d0*xbs*bxbs*(lbxb.cq.bxb) &
     &  - 3.d0*xt*((lbxt.cp.bxt)+lxt*co) &
     &  - 3.d0*xb*((lbxb.cp.bxb)+lxb*co) &
     &  - 3.d0*xc*((lbxc.cp.bxc)+lxc*co) &
     &  - 3.d0*xs*((lbxs.cp.bxs)+lxs*co) &
     &  - 3.d0*xu*((lbxu.cp.bxu)+lxu*co) &
     &  - 3.d0*xd*((lbxd.cp.bxd)+lxd*co) &
     &  - xtau*((lbxtau.cp.bxtau)+lxtau*co) &
     &  - xmu*((lbxmu.cp.bxmu)+lxmu*co) &
     &  - xe*((lbxe.cp.bxe)+lxe*co)+co*( &
     &  - ls*(3.d0*(xt+xb+xs+xc+xd+xu)+xtau+xmu+xe) &
     & +3.d0*xt*(2.d0+xts*bxts) &
     & +3.d0*xb*(2.d0+xbs*bxbs) &
     & +3.d0*xc*(2.d0+xcs*bxcs) &
     & +3.d0*xs*(2.d0+xss*bxss) &
     & +3.d0*xu*(2.d0+xus*bxus) &
     & +3.d0*xd*(2.d0+xds*bxds) &
     & +xtau*(2.d0+xtaus*bxtaus) &
     & +xmu*(2.d0+xmus*bxmus) &
     & +xe*(2.d0+xes*bxes)))
!
      dmw= (9.d0*(bh.cp.lbh) &
     &    +fz*(bz.cp.lbz) &
     &     - xw*fw*(lbw.cq.bw) &
     &    +2.d0*(1.d0-12.d0*xw*xw)*(bw.cp.lbw) &
     &     +co*( &
     &     2.d0*lw+(1.d0-6.d0*xz)*lz+6.d0*(1.d0-xz)*ls &
     &     - 2.d0*(9.d0-xw*(2.d0+36.d0*xw)) &
     &    +2.d0*xz*(5.d0-12.d0*xz)))/(128.d0*xw*xw)
      dmw= dmw-1.d0/(32.d0*xw*xw)*( &
     &   - 3.d0*xt*bxts*(bxt.cp.lbxt)-3.d0*xt*lxt*co &
     &   - 3.d0*xb*bxbs*(bxb.cp.lbxb)-3.d0*xb*lxb*co &
     &   - 3.d0*xs*bxss*(bxs.cp.lbxs)-3.d0*xs*lxs*co &
     &   - 3.d0*xc*bxcs*(bxc.cp.lbxc)-3.d0*xc*lxc*co &
     &   - 3.d0*xd*bxds*(bxd.cp.lbxd)-3.d0*xd*lxd*co &
     &   - 3.d0*xu*bxus*(bxu.cp.lbxu)-3.d0*xd*lxu*co &
     &   - xtau*bxtaus*(bxtau.cp.lbxtau)-xtau*lxtau*co &
     &   - xmu*bxmus*(bxmu.cp.lbxmu)-xmu*lxmu*co &
     &   - xe*bxes*(bxe.cp.lbxe)-xe*lxe*co &
     &  +co*( &
     &   - ls*(3.d0*(xt+xb+xs+xc+xd+xu)+xtau+xmu+xe) &
     &  +6.d0*xt*(1.d0-2.d0*xt) &
     &  +6.d0*xb*(1.d0-2.d0*xb) &
     &  +6.d0*xs*(1.d0-2.d0*xs) &
     &  +6.d0*xc*(1.d0-2.d0*xc) &
     &  +6.d0*xd*(1.d0-2.d0*xd) &
     &  +6.d0*xu*(1.d0-2.d0*xu) &
     &  +2.d0*xtau*(1.d0-2.d0*xtau) &
     &  +2.d0*xmu*(1.d0-2.d0*xmu) &
     &  +2.d0*xe*(1.d0-2.d0*xe)))
!
      dmz= (-0.25d0+(xz-3.d0*xz))*(lbz.cq.bz) &
     &   +2.d0*(1-6.d0*xz)*(bz.cp.lbz)+ &
     &    co*(-4.d0+30.d0*xz+3.d0*(lz+ls))
!
      dmz= dmz/(64.d0*xw)
!
      d2mh= (2.d0*fw*bwis*co &
     &     +fz*bzis*co &
     &     -9.d0*(lbh.cp.bh) &
     &     +2.d0*xz*xz*fz*(lbz.cp.bzic) &
     &     -2.d0*xz*(1.d0-12.d0*xz)*(lbz.cp.bzi) &
     &     -(lbz.cp.bz) &
     &     +4.d0*xw*xw*fz*(lbw.cp.bwic) &
     &     -4.d0*xw*(1.d0-12.d0*xw)*(lbw.cp.bwi) &
     &     -2.d0*(lbw.cp.bw) &
     &     +co*( &
     &      -6.d0*ls-lz-2.d0*lw+9.d0+4.d0*xw*(1+3.d0*xw)+ &
     &       2.d0*xz*(1.d0+3.d0*xz)))/(64.d0*xw*ps)
!
      d2mh= d2mh+1.d0/(32.d0*xw*ps)*( &
     &      4.d0*xec*(lbxe*bxeis) &
     &     +2.d0*bxei &
     &     -16.d0*xec*(lbxe.cp.bxei) &
     &     +4.d0*xmuc*(lbxmu*bxmuis) &
     &     +2.d0*bxmui &
     &     -16.d0*xmuc*(lbxmu.cp.bxmui) &
     &     +4.d0*xtauc*(lbxtau*bxtauis) &
     &     +2.d0*bxtaui &
     &     -16.d0*xtauc*(lbxtau.cp.bxtaui) &
     &     +12.d0*xuc*(lbxu*bxuis) &
     &     +6.d0*bxui &
     &     -48.d0*xuc*(lbxu.cp.bxui) &
     &     +12.d0*xdc*(lbxd*bxdis) &
     &     +6.d0*bxdi &
     &     -48.d0*xdc*(lbxd.cp.bxdi) &
     &     +12.d0*xcc*(lbxc*bxcis) &
     &     +6.d0*bxci &
     &     -48.d0*xcc*(lbxc.cp.bxci) &
     &     +12.d0*xsc*(lbxs*bxsis) &
     &     +6.d0*bxsi &
     &     -48.d0*xsc*(lbxs.cp.bxsi) &
     &     +12.d0*xtc*(lbxt*bxtis) &
     &     +6.d0*bxti &
     &     -48.d0*xtc*(lbxu.cp.bxti) &
     &     +12.d0*xbc*(lbxb*bxbis) &
     &     +6.d0*bxbi &
     &     -48.d0*xbc*(lbxb.cp.bxbi)+ &
     &     co*( &
     &     xec*bxes &
     &     +xmuc*bxmus &
     &     +xtauc*bxtaus &
     &     +xuc*bxus &
     &     +xdc*bxds &
     &     +xcc*bxcs &
     &     +xsc*bxss &
     &     +xtc*bxts &
     &     +xbc*bxbs))
!
      dmwh= (2.d0*xw*fw*bwis*co &
     &   +9.d0*(lbh.cp.bh) &
     &   +xz*fz*(lbz.cq.bz) &
     &   +(1.d0-2.d0*xz)*(lbz.cp.bz) &
     &    - 4.d0*xw*xw*fw*(lbw.cp.bwic) &
     &    - xw*(1.d0+xw*xw*(-10.d0+48.d0*xw))*(lbw.cq.bw) &
     &   +2.d0*(lbw.cp.bw) &
     &   +co*( &
     &   +2.d0*lw &
     &   +(1.d0-3.d0*xz)*lz &
     &   +3.d0*(2.d0-xz)*ls &
     &    - (15.d0+xw*(-2.d0+12.d0*xw)) &
     &   +3.d0*xz*(1.d0+2.d0*xz)))/(64.d0*xw*xw*ps)
!
      dmzh= (2.d0*fz*bzis*co &
     &    - 4.d0*xs*fz*(lbz.cp.bzic) &
     &    - 3.d0*(1.d0+xz*(-6.d0+24.d0*xz))*(lbz.cq.bz) &
     &   +4.d0*(lbz.cp.bz) &
     &   +co*( &
     &   +6.d0*ls &
     &   +6.d0*lz &
     &    - 4.d0*(1.d0+6.d0*xz)))/(128.d0*xw*ps)
!
      dmzw= (fz*(lbz.cq.bz) &
     &    - 8.d0*(1.d0-6.d0*xz)*(lbz.cp.bz) &
     &   +co*( &
     &    - 12.d0*ls &
     &    - 12.d0*lz &
     &   +8.d0*(2.d0-15.d0*xz)))/(256.d0*xw*xw*ps)
!
      dmww= (xw*(1.d0-4.d0*xw*(1.d0-2.d0*xw))*bwis*co &
     &    - 18.d0*(lbh.cp.bh) &
     &    - 2.d0*fz*(lbz.cp.bz) &
     &    - 2.d0*xw*xw*fw*(lbw.cp.bwic) &
     &   +2.d0*xw*(1.d0-12.d0*xw*xw)*(lbw.cq.bw) &
     &    - 4.d0*(lbw.cp.bw) &
     &   +co*(  &
     &    - 4.d0*lw &
     &    - 2.d0*(1.d0-6.d0*xz)*lz &
     &    - 12.d0*(1.d0-xz)*ls &
     &   +4.d0*(9.d0-xw*(1.d0-6.d0*xw)) &
     &    - 4.d0*xz*(5.d0-12.d0*xz)))/(128.d0*xw*xw*xw*ps)
      dmww= dmww+1.d0/(162.d0*xw*xw*xw*ps)*( &
     &   - 3.d0*xt*bxts*(bxt.cp.lbxt)-3.d0*xt*lxt*co &
     &   - 3.d0*xb*bxbs*(bxb.cp.lbxb)-3.d0*xb*lxb*co &
     &   - 3.d0*xs*bxss*(bxs.cp.lbxs)-3.d0*xs*lxs*co &
     &   - 3.d0*xc*bxcs*(bxc.cp.lbxc)-3.d0*xc*lxc*co &
     &   - 3.d0*xd*bxds*(bxd.cp.lbxd)-3.d0*xd*lxd*co &
     &   - 3.d0*xu*bxus*(bxu.cp.lbxu)-3.d0*xd*lxu*co &
     &   - xtau*bxtaus*(bxtau.cp.lbxtau)-xtau*lxtau*co &
     &   - xmu*bxmus*(bxmu.cp.lbxmu)-xmu*lxmu*co &
     &   - xe*bxes*(bxe.cp.lbxe)-xe*lxe*co+co*( &
     &   - ls*(3.d0*(xt+xb+xs+xc+xd+xu)+xtau+xmu+xe) &
     &  +6.d0*xt*(1.d0-2.d0*xt) &
     &  +6.d0*xb*(1.d0-2.d0*xb) &
     &  +6.d0*xs*(1.d0-2.d0*xs) &
     &  +6.d0*xc*(1.d0-2.d0*xc) &
     &  +6.d0*xd*(1.d0-2.d0*xd) &
     &  +6.d0*xu*(1.d0-2.d0*xu) &
     &  +2.d0*xtau*(1.d0-2.d0*xtau) &
     &  +2.d0*xmu*(1.d0-2.d0*xmu) &
     &  +2.d0*xe*(1.d0-2.d0*xe)))
!
      dmzz= (fz*bzis*co &
     &    - 2.d0*xz*fz*(lbz.cp.bzic) &
     &   +8.d0*xz*(1.d0-6.d0*xz)*(lbz.cq.bz)  &
     &    - 48.d0*xz*(lbz.cp.bz) &
     &   +4.d0*(1.d0+42.d0*xz)*co)/(256.d0*xw*xz*ps)
!
      value(1,1:2)= sh(1:2)
      value(2,1:2)= dmh(1:2)
      value(3,1:2)= dmw(1:2)
      value(4,1:2)= dmz(1:2)
      value(5,1:2)= d2mh(1:2)
      value(6,1:2)= dmwh(1:2)
      value(7,1:2)= dmzh(1:2)
      value(8,1:2)= dmzw(1:2)
      value(9,1:2)= dmww(1:2)
      value(10,1:2)= dmzz(1:2)
!
      RETURN
!
      END FUNCTION VBFNLO_HIGGSTO_deriv
!
!--------------------------------------------------------------------------
!----- Main ---------------------------------------------------------------
!--------------------------------------------------------------------------
!
      SUBROUTINE VBFNLO_HIGGSTO_pole(m,cpgh)
      !USE VBFNLO_HIGGSTO
      ! USE VBFNLO_HIGGSTO_masses
      ! USE VBFNLO_HIGGSTO_riemann
      ! USE VBFNLO_HIGGSTO_sp_fun
      ! USE VBFNLO_HIGGSTO_aux_Hcp
!
      IMPLICIT NONE
!

      REAL*8 muh,cpgh,gos,mhb,ghb,m,expgHi,EWC,ghi
      REAL*8, dimension(10,2) :: expc
!
!mr
#include "VBFNLO/utilities/koppln.inc"
      DOUBLE PRECISION CLR,XM2,XMG,B,Vcoupl,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),Vcoupl(4,5),A(4,5)
!mr
      ! INTERFACE
      !  FUNCTION VBFNLO_HIGGSTO_deriv(scal,rhm) RESULT(value)
      !  !USE VBFNLO_HIGGSTO
      !  ! USE VBFNLO_HIGGSTO_masses
      !  ! USE VBFNLO_HIGGSTO_acmplx_pro
      !  ! USE VBFNLO_HIGGSTO_acmplx_rat
      !  ! USE VBFNLO_HIGGSTO_cmplx_root
      !  ! USE VBFNLO_HIGGSTO_full_ln
      !  ! USE VBFNLO_HIGGSTO_qcd
      !  ! USE VBFNLO_HIGGSTO_units
      !  IMPLICIT NONE
      !  REAL*8 scal,rhm
      !  REAL*8, dimension(10,2) :: value 
      !  END FUNCTION VBFNLO_HIGGSTO_deriv
      ! END INTERFACE
!
      CALL VBFNLO_HIGGSTO_init_niels
!
      muh= m 
!
!mr
      mt = xmt
      mw= sqrt(xm2(3))
      mz= sqrt(xm2(2))
      mbq= xmb
      mb= xmb
      imw= xmg(3)/mw
      imz= xmg(2)/mz
      swr= mw*mw-imw*imw
      swi= -mw*imw*(1.d0-0.5d0*(imw*imw)/(mw*mw))
      szr= mz*mz-imz*imz
      szi= -mz*imz*(1.d0-0.5d0*(imz*imz)/(mz*mz))
      qcdc = 1
      gtop = 1
      yimt = 0d0
!mr
!mr      print 1,mt
!mr 1    format(' top mass [GeV] ',e20.5)
!mr      print 2,muh
!mr 2    format(' muh      [GeV] ',e20.5)
!
      CALL VBFNLO_HIGGSTO_gridHt(muh,gos)
!
!mr      print 3,gos
!mr 3    format(' GH OS    [GeV] =',e20.7)
!
!mr      print*,'   '
!mr      print*,'*****************************************************'
!mr      print*,'   '
!mr      print*,'   '
!
      IF(muh.le.200.d0) THEN 
       expc= VBFNLO_HIGGSTO_deriv(muh,muh)
       EWC= 4.d0*sqrt(2.d0)*1.16637d-5*(mw*mw)/pis
       expgHi= gos &
     &      -EWC*EWC*expc(2,2)*expc(2,2)*gos &
     &      +0.5d0*EWC*(expc(2,2)/muh-expc(5,2)*muh)*gos*gos
       cpgh= expgHi
!mr       print 6,cpgH
      ELSEIF((muh > 200.d0).and.(muh < 240.d0)) THEN      
       CALL VBFNLO_HIGGSTO_gridlow(muh,ghi)
       cpgh= ghi
!mr       print 6,cpgH
      ELSE 
       CALL VBFNLO_HIGGSTO_GH(muh,cpgh)
      ENDIF
!
      mhb= sqrt(muh*muh+cpgh*cpgh)
      ghb= mhb/muh*cpgh
!
!mr      print*,'   '
!mr      print*,'*** Bar scheme **************************************'
!mr      print*,'   '
!mr      print 4,mhb
!mr 4    format(' MHB      [GeV] = ',e20.5)
!mr      print 5,ghb
!mr 5    format(' GHB      [GeV] = ',e20.5)
!mr 6    format(' gH Total [GeV] = ',e20.7)
!
!mr      print*,'   '
!mr      print*,'*****************************************************'
!mr      print*,'   '
!mr      print*,'   '
!
      RETURN
!
      END SUBROUTINE VBFNLO_HIGGSTO_pole
!
!---------------------------------------------------------------------------------------
!
      SUBROUTINE VBFNLO_HIGGSTO_GH(muh,cpgh)
      !USE VBFNLO_HIGGSTO_masses
      !USE VBFNLO_HIGGSTO_aux_Hcp
      !USE VBFNLO_HIGGSTO_aux_hbb
      !USE VBFNLO_HIGGSTO_cmplx_rootz
      !USE VBFNLO_HIGGSTO_acmplx_pro
      !USE VBFNLO_HIGGSTO_acmplx_rat
      !USE VBFNLO_HIGGSTO_full_ln
      !USE VBFNLO_HIGGSTO_units
      !USE VBFNLO_HIGGSTO_ferbos
      !USE VBFNLO_HIGGSTO_rootW
      !USE VBFNLO_HIGGSTO_HBB_cp
      !USE VBFNLO_HIGGSTO_set_phys_const
      !USE VBFNLO_HIGGSTO_root_find2
!
      IMPLICIT NONE
!
      INTEGER i
      REAL*8 muh,scals,x1,x2,xacc,tgH,itgH,gHf,gHb,gHt,cpgh
      REAL*8, dimension(2) :: cmw,cmz
!
!mr
#include "VBFNLO/utilities/koppln.inc"
      DOUBLE PRECISION CLR,XM2,XMG,B,Vcoupl,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),Vcoupl(4,5),A(4,5)
!
!mr
      mt = xmt
      mw= sqrt(xm2(3))
      mz= sqrt(xm2(2))
      mbq= xmb
      mb= xmb
      imw= xmg(3)/mw
      imz= xmg(2)/mz
      swr= mw*mw-imw*imw
      swi= -mw*imw*(1.d0-0.5d0*(imw*imw)/(mw*mw))
      szr= mz*mz-imz*imz
      szi= -mz*imz*(1.d0-0.5d0*(imz*imz)/(mz*mz))
      qcdc = 1
      gtop = 1
      yimt = 0d0
!mr
      scalec= muh
      scals= scalec*scalec
!
      xtop= mt*mt/(muh*muh)
      xb= mb*mb/(muh*muh)
!
      cxe= (me*me)/scals
      cxmu= (mm*mm)/scals
      cxtau= (mtl*mtl)/scals
      cxu= (muq*muq)/scals
      cxd= (mdq*mdq)/scals
      cxc= (mcq*mcq)/scals
      cxs= (msq*msq)/scals
      cxt= (mt*mt)/scals
      cxb= (mb*mb)/scals
!
      cxes= cxe*cxe
      cxmus= cxmu*cxmu
      cxtaus= cxtau*cxtau
      cxus= cxu*cxu
      cxds= cxd*cxd
      cxcs= cxc*cxc
      cxss= cxs*cxs
      cxts= cxt*cxt
      cxbs= cxb*cxb
      cxtmb= cxt-cxb
      cxtmbs= cxtmb*cxtmb
      cxtmbi= 1.d0/cxtmb
!
      cxtc= cxts*cxt
      cxbc= cxbs*cxb
!
      clxe= log(cxe)
      clxmu= log(cxmu)
      clxtau= log(cxtau)
      clxu= log(cxu)
      clxd= log(cxd)
      clxc= log(cxc)
      clxs= log(cxs)
      clxt= log(cxt)
      clxb= log(cxb)
!
      cmw(1)= swr
      cmw(2)= swi
      cmz(1)= szr
      cmz(2)= szi
!
      rcmw= cmw(1).crz.cmw(2)
!
      cxw(1)= swr/scals
      cxw(2)= swi/scals
      cxz(1)= szr/scals
      cxz(2)= szi/scals
      cxws= cxw.cp.cxw
      cxwc= cxws.cp.cxw
      cxwi= co.cq.cxw
!
      clw= cxw(1).fln.cxw(2)
!
      ccts= cmw.cq.cmz
      csts= co-ccts
      cctq= ccts.cp.ccts
      cctvi= cctq.cp.ccts
      clcts= ccts(1).fln.ccts(2)
!
      xq= cxt+cxb
      cxq= cxb+cxc+cxs+cxu+cxd
      cxqs= cxbs+cxcs+cxss+cxus+cxds
      cxl= cxtau+cxmu+cxe
      cxls= cxtaus+cxmus+cxes
!
      clwtb= cxtmbs*co+cxws-2.d0*xq*cxw
!
      cpw= swr*co+swi*ci
      cpz= szr*co+szi*ci
!
      g_f= 1.16637d-5
!
      muhcp= muh
!
      IF(muhcp.ge.900.d0) THEN
       x2= 1.d3/muh
      ELSEIF((muhcp.gt.600.d0).and.(muhcp.lt.900.d0)) THEN
       x2= 4.d2/muh
      ELSEIF(muhcp.lt.160.d0) THEN
       x2= 1.d0/muh
      ELSE
       x2= 2.d2/muh
      ENDIF
      xacc= 1.d-10
      DO i=3,0,-1
       ifb= i
       x1= -1.d0/muh        
!
       inc= 0
       tgH= muh*VBFNLO_HIGGSTO_zeroin(VBFNLO_HIGGSTO_SSHH,x1,x2,zero,xacc)
!
! T
!
       IF(i==0) THEN
!
! LF
!
       ELSEIF(i==1) THEN
        gHf= tgH
!
! top
!
       ELSEIF(i==2) THEN
        gHt= tgH
!
! bos
!
       ELSEIF(i==3) THEN
        gHb= tgH
!
       ENDIF

       IF(inc==1) THEN
        print*,'---- Warning ------------------------ ' 
        IF(i==0) THEN
         print*,' inconsistent interval for tot '
        ELSEIF(i==1) THEN
         print*,' inconsistent interval for lf '
        ELSEIF(i==2) THEN
         print*,' inconsistent interval for top '
        ELSEIF(i==3) THEN
         print*,' inconsistent interval for bos '
        ENDIF
        print*,'------------------------------------ '
       ENDIF
      ENDDO
!
!      print*,'+++++ gh details  +++++++++++++++++++++++++++++++++++'
!      print*,'   '
      xacc= 1.d-20
      DO i=3,0,-1
       ifb= i
       IF(i==0) THEN
        x1= 0.9d0*tgH/muh
        x2= 1.1d0*tgH/muh
       ELSE IF(i==1) THEN
        x1= 0.9d0*gHf/muh
        x2= 1.1d0*gHf/muh
       ELSE IF(i==2) THEN
        x1= 0.9d0*gHt/muh
        x2= 1.1d0*gHt/muh
       ELSE IF(i==3) THEN
        x1= 0.9d0*gHb/muh
        x2= 1.1d0*gHb/muh
       ENDIF
       inc= 0
       itgH= muh*VBFNLO_HIGGSTO_zeroin(VBFNLO_HIGGSTO_SSHH,x1,x2,zero,xacc)
       IF(i==0) THEN
!        print 116,itgH
        cpgh= itgh
       ELSEIF(i==1) THEN
!        print 117,itgH
       ELSEIF(i==2) THEN
!        print 1118,itgH
       ELSEIF(i==3) THEN
!        print 118,itgH
       ENDIF
       IF(inc==1) THEN
        print*,'---- Warning ------------------------ ' 
        IF(i==0) THEN
         print*,' inconsistent interval for tot '
        ELSEIF(i==1) THEN
         print*,' inconsistent interval for lf '
        ELSEIF(i==2) THEN
         print*,' inconsistent interval for top '
        ELSEIF(i==3) THEN
         print*,' inconsistent interval for bos '
        ENDIF
        print*,'------------------------------------ '
       ENDIF
      ENDDO
!      print*,'   '
!      print*,'+++++++++++++++++++++++++++++++++++++++++++++++++++++'
!
!
      RETURN
!
      END SUBROUTINE VBFNLO_HIGGSTO_GH
!
!---------------------------------------------------------------------------------------
!
      ! CONTAINS
      SUBROUTINE VBFNLO_HIGGSTO_hbrd(VBFNLO_HIGGSTO_FCN,n,x,fvec,epsfcn,tol,info,diag)
! 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-07-15  Time: 13:27:42
!
      INTEGER, INTENT(IN)     :: n
      REAL*8, INTENT(IN OUT)  :: x(n)
      REAL*8, INTENT(IN OUT)  :: fvec(n)
      REAL*8, INTENT(IN)      :: epsfcn
      REAL*8, INTENT(IN)      :: tol
      INTEGER, INTENT(OUT)    :: info
      REAL*8, INTENT(OUT)     :: diag(n)
!
! EXTERNAL fcn
!
      INTERFACE
       SUBROUTINE VBFNLO_HIGGSTO_FCN(N,X,FVEC,IFLAG)
       IMPLICIT NONE
       INTEGER, INTENT(IN)      :: n
       REAL*8, INTENT(IN)    :: x(n)
       REAL*8, INTENT(OUT)   :: fvec(n)
       INTEGER, INTENT(IN OUT)  :: iflag
       END SUBROUTINE VBFNLO_HIGGSTO_FCN
      END INTERFACE

!   **********

!   SUBROUTINE VBFNLO_HIGGSTO_HBRD

!   THE PURPOSE OF HBRD IS TO FIND A ZERO OF A SYSTEM OF N NONLINEAR
!   FUNCTIONS IN N VARIABLES BY A MODIFICATION OF THE POWELL HYBRID METHOD.
!   THIS IS DONE BY USING THE MORE GENERAL NONLINEAR EQUATION SOLVER HYBRD.
!   THE USER MUST PROVIDE A SUBROUTINE VBFNLO_HIGGSTO_WHICH CALCULATES THE FUNCTIONS.
!   THE JACOBIAN IS THEN CALCULATED BY A FORWARD-DIFFERENCE APPROXIMATION.

!   THE SUBROUTINE VBFNLO_HIGGSTO_STATEMENT IS

!     SUBROUTINE VBFNLO_HIGGSTO_HBRD(N,X,FVEC,EPSFCN,TOL,INFO,WA,LWA)

!   WHERE

!     FCN IS THE NAME OF THE USER-SUPPLIED SUBROUTINE VBFNLO_HIGGSTO_WHICH CALCULATES
!       THE FUNCTIONS.  FCN MUST BE DECLARED IN AN EXTERNAL STATEMENT
!       IN THE USER CALLING PROGRAM, AND SHOULD BE WRITTEN AS FOLLOWS.

!       SUBROUTINE VBFNLO_HIGGSTO_FCN(N,X,FVEC,IFLAG)
!       INTEGER N,IFLAG
!       REAL X(N),FVEC(N)
!       ----------
!       CALCULATE THE FUNCTIONS AT X AND RETURN THIS VECTOR IN FVEC.
!       ---------
!       RETURN
!       END

!       THE VALUE OF IFLAG NOT BE CHANGED BY FCN UNLESS
!       THE USER WANTS TO TERMINATE THE EXECUTION OF HBRD.
!       IN THIS CASE SET IFLAG TO A NEGATIVE INTEGER.

!     N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!       OF FUNCTIONS AND VARIABLES.

!     X IS AN ARRAY OF LENGTH N. ON INPUT X MUST CONTAIN AN INITIAL
!       ESTIMATE OF THE SOLUTION VECTOR.  ON OUTPUT X CONTAINS THE
!       FINAL ESTIMATE OF THE SOLUTION VECTOR.

!     FVEC IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS
!       THE FUNCTIONS EVALUATED AT THE OUTPUT X.

!     EPSFCN IS AN INPUT VARIABLE USED IN DETERMINING A SUITABLE STEP LENGTH
!       FOR THE FORWARD-DIFFERENCE APPROXIMATION.  THIS APPROXIMATION ASSUMES
!       THAT THE RELATIVE ERRORS IN THE FUNCTIONS ARE OF THE ORDER OF EPSFCN.
!       IF EPSFCN IS LESS THAN THE MACHINE PRECISION, IT IS ASSUMED THAT THE
!       RELATIVE ERRORS IN THE FUNCTIONS ARE OF THE ORDER OF THE MACHINE
!       PRECISION.

!     TOL IS A NONNEGATIVE INPUT VARIABLE.  TERMINATION OCCURS WHEN THE
!       ALGORITHM ESTIMATES THAT THE RELATIVE ERROR BETWEEN X AND THE SOLUTION
!       IS AT MOST TOL.

!     INFO IS AN INTEGER OUTPUT VARIABLE.  IF THE USER HAS TERMINATED
!       EXECUTION, INFO IS SET TO THE (NEGATIVE) VALUE OF IFLAG.
!       SEE DESCRIPTION OF FCN.  OTHERWISE, INFO IS SET AS FOLLOWS.

!       INFO= 0   IMPROPER INPUT PARAMETERS.

!       INFO= 1   ALGORITHM ESTIMATES THAT THE RELATIVE ERROR
!                  BETWEEN X AND THE SOLUTION IS AT MOST TOL.

!       INFO= 2   NUMBER OF CALLS TO FCN HAS REACHED OR EXCEEDED 200*(N+1).

!       INFO= 3   TOL IS TOO SMALL. NO FURTHER IMPROVEMENT IN
!                  THE APPROXIMATE SOLUTION X IS POSSIBLE.

!       INFO= 4   ITERATION IS NOT MAKING GOOD PROGRESS.

!   SUBPROGRAMS CALLED

!     USER-SUPPLIED ...... FCN

!     MINPACK-SUPPLIED ... HYBRD

!   ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!   BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE

! Reference:
! Powell, M.J.D. 'A hybrid method for nonlinear equations' in Numerical Methods
!      for Nonlinear Algebraic Equations', P.Rabinowitz (editor), Gordon and
!      Breach, London 1970.
!   **********
      INTEGER :: maxfev, ml,mode,mu,nfev,nprint
      REAL*8  :: xtol
      REAL*8, PARAMETER  :: factor= 100.0d0,zero= 0.0d0
      
      info= 0

!     CHECK THE INPUT PARAMETERS FOR ERRORS.

      IF(n <= 0.or.epsfcn < zero.or.tol < zero) GO TO 20

!     CALL VBFNLO_HIGGSTO_HYBRD.

      maxfev= 200*(n+1)
      xtol= tol
      ml= n-1
      mu= n-1
      mode= 2
      nprint= 0
      CALL VBFNLO_HIGGSTO_hybrd(VBFNLO_HIGGSTO_FCN,n,x,fvec,xtol,maxfev,ml,mu,epsfcn,diag, &
     &               mode,factor,nprint,info,nfev)
      IF(info == 5) info= 4
   20 RETURN
       
!     LAST CARD OF SUBROUTINE VBFNLO_HIGGSTO_HBRD.

      END SUBROUTINE VBFNLO_HIGGSTO_hbrd
!
!------------------------------------------------------------------------
!
      SUBROUTINE VBFNLO_HIGGSTO_hybrd(VBFNLO_HIGGSTO_FCN,n,x,fvec,xtol,maxfev,ml,mu,epsfcn, &
     &                     diag,mode,factor,nprint,info,nfev)

      INTEGER, INTENT(IN)        :: n
      REAL*8, INTENT(IN OUT)  :: x(n)
      REAL*8, INTENT(IN OUT)  :: fvec(n)
      REAL*8, INTENT(IN)      :: xtol
      INTEGER, INTENT(IN OUT)    :: maxfev
      INTEGER, INTENT(IN OUT)    :: ml
      INTEGER, INTENT(IN)        :: mu
      REAL*8, INTENT(IN)      :: epsfcn
      REAL*8, INTENT(OUT)     :: diag(n)
      INTEGER, INTENT(IN)        :: mode
      REAL*8, INTENT(IN)      :: factor
      INTEGER, INTENT(IN OUT)    :: nprint
      INTEGER, INTENT(OUT)       :: info
      INTEGER, INTENT(OUT)       :: nfev

! EXTERNAL fcn

      INTERFACE
       SUBROUTINE VBFNLO_HIGGSTO_FCN(N,X,FVEC,IFLAG)
       IMPLICIT NONE
       INTEGER, INTENT(IN)      :: n
       REAL*8, INTENT(IN)    :: x(n)
       REAL*8, INTENT(OUT)   :: fvec(n)
       INTEGER, INTENT(IN OUT)  :: iflag
       END SUBROUTINE VBFNLO_HIGGSTO_FCN
      END INTERFACE

!   **********

!   SUBROUTINE VBFNLO_HIGGSTO_HYBRD

!   THE PURPOSE OF HYBRD IS TO FIND A ZERO OF A SYSTEM OF N NONLINEAR
!   FUNCTIONS IN N VARIABLES BY A MODIFICATION OF THE POWELL HYBRID METHOD.
!   THE USER MUST PROVIDE A SUBROUTINE VBFNLO_HIGGSTO_WHICH CALCULATES THE FUNCTIONS.
!   THE JACOBIAN IS THEN CALCULATED BY A FORWARD-DIFFERENCE APPROXIMATION.

!   THE SUBROUTINE VBFNLO_HIGGSTO_STATEMENT IS

!     SUBROUTINE VBFNLO_HIGGSTO_HYBRD(FCN,N,X,FVEC,XTOL,MAXFEV,ML,MU,EPSFCN,
!                      DIAG,MODE,FACTOR,NPRINT,INFO,NFEV,FJAC,
!                      LDFJAC,R,LR,QTF,WA1,WA2,WA3,WA4)

!   WHERE

!     FCN IS THE NAME OF THE USER-SUPPLIED SUBROUTINE VBFNLO_HIGGSTO_WHICH CALCULATES
!       THE FUNCTIONS.  FCN MUST BE DECLARED IN AN EXTERNAL STATEMENT IN
!       THE USER CALLING PROGRAM, AND SHOULD BE WRITTEN AS FOLLOWS.

!       SUBROUTINE VBFNLO_HIGGSTO_FCN(N, X, FVEC, IFLAG)
!       INTEGER N, IFLAG
!       REAL X(N), FVEC(N)
!       ----------
!       CALCULATE THE FUNCTIONS AT X AND
!       RETURN THIS VECTOR IN FVEC.
!       ---------
!       RETURN
!       END

!       THE VALUE OF IFLAG SHOULD NOT BE CHANGED BY FCN UNLESS
!       THE USER WANTS TO TERMINATE EXECUTION OF HYBRD.
!       IN THIS CASE SET IFLAG TO A NEGATIVE INTEGER.

!     N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!       OF FUNCTIONS AND VARIABLES.

!     X IS AN ARRAY OF LENGTH N.  ON INPUT X MUST CONTAIN AN INITIAL
!       ESTIMATE OF THE SOLUTION VECTOR.  ON OUTPUT X CONTAINS THE FINAL
!       ESTIMATE OF THE SOLUTION VECTOR.

!     FVEC IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS
!       THE FUNCTIONS EVALUATED AT THE OUTPUT X.

!     XTOL IS A NONNEGATIVE INPUT VARIABLE.  TERMINATION OCCURS WHEN THE
!       RELATIVE ERROR BETWEEN TWO CONSECUTIVE ITERATES IS AT MOST XTOL.

!     MAXFEV IS A POSITIVE INTEGER INPUT VARIABLE.  TERMINATION OCCURS WHEN
!       THE NUMBER OF CALLS TO FCN IS AT LEAST MAXFEV BY THE END OF AN
!       ITERATION.

!     ML IS A NONNEGATIVE INTEGER INPUT VARIABLE WHICH SPECIFIES THE
!       NUMBER OF SUBDIAGONALS WITHIN THE BAND OF THE JACOBIAN MATRIX.
!       IF THE JACOBIAN IS NOT BANDED, SET ML TO AT LEAST N - 1.

!     MU IS A NONNEGATIVE INTEGER INPUT VARIABLE WHICH SPECIFIES THE NUMBER
!       OF SUPERDIAGONALS WITHIN THE BAND OF THE JACOBIAN MATRIX.
!       IF THE JACOBIAN IS NOT BANDED, SET MU TO AT LEAST N - 1.

!     EPSFCN IS AN INPUT VARIABLE USED IN DETERMINING A SUITABLE STEP LENGTH
!       FOR THE FORWARD-DIFFERENCE APPROXIMATION.  THIS APPROXIMATION
!       ASSUMES THAT THE RELATIVE ERRORS IN THE FUNCTIONS ARE OF THE ORDER
!       OF EPSFCN. IF EPSFCN IS LESS THAN THE MACHINE PRECISION,
!       IT IS ASSUMED THAT THE RELATIVE ERRORS IN THE FUNCTIONS ARE OF THE
!       ORDER OF THE MACHINE PRECISION.

!     DIAG IS AN ARRAY OF LENGTH N. IF MODE= 1 (SEE BELOW),
!       DIAG IS INTERNALLY SET.  IF MODE= 2, DIAG MUST CONTAIN POSITIVE
!       ENTRIES THAT SERVE AS MULTIPLICATIVE SCALE FACTORS FOR THE
!       VARIABLES.

!     MODE IS AN INTEGER INPUT VARIABLE. IF MODE= 1, THE VARIABLES WILL BE
!       SCALED INTERNALLY.  IF MODE= 2, THE SCALING IS SPECIFIED BY THE
!       INPUT DIAG.  OTHER VALUES OF MODE ARE EQUIVALENT TO MODE= 1.

!     FACTOR IS A POSITIVE INPUT VARIABLE USED IN DETERMINING THE
!       INITIAL STEP BOUND. THIS BOUND IS SET TO THE PRODUCT OF
!       FACTOR AND THE EUCLIDEAN NORM OF DIAG*X IF NONZERO, OR ELSE
!       TO FACTOR ITSELF. IN MOST CASES FACTOR SHOULD LIE IN THE
!       INTERVAL (.1,100.).100. IS A GENERALLY RECOMMENDED VALUE.

!     NPRINT IS AN INTEGER INPUT VARIABLE THAT ENABLES CONTROLLED
!       PRINTING OF ITERATES IF IT IS POSITIVE. IN THIS CASE,
!       FCN IS CALLED WITH IFLAG= 0 AT THE BEGINNING OF THE FIRST
!       ITERATION AND EVERY NPRINT ITERATIONS THEREAFTER AND
!       IMMEDIATELY PRIOR TO RETURN, WITH X AND FVEC AVAILABLE
!       FOR PRINTING. IF NPRINT IS NOT POSITIVE, NO SPECIAL CALLS
!       OF FCN WITH IFLAG= 0 ARE MADE.

!     INFO IS AN INTEGER OUTPUT VARIABLE. IF THE USER HAS
!       TERMINATED EXECUTION, INFO IS SET TO THE (NEGATIVE)
!       VALUE OF IFLAG. SEE DESCRIPTION OF FCN. OTHERWISE,
!       INFO IS SET AS FOLLOWS.

!       INFO= 0   IMPROPER INPUT PARAMETERS.

!       INFO= 1   RELATIVE ERROR BETWEEN TWO CONSECUTIVE ITERATES
!                  IS AT MOST XTOL.

!       INFO= 2   NUMBER OF CALLS TO FCN HAS REACHED OR EXCEEDED MAXFEV.

!       INFO= 3   XTOL IS TOO SMALL. NO FURTHER IMPROVEMENT IN
!                  THE APPROXIMATE SOLUTION X IS POSSIBLE.

!       INFO= 4   ITERATION IS NOT MAKING GOOD PROGRESS, AS
!                  MEASURED BY THE IMPROVEMENT FROM THE LAST
!                  FIVE JACOBIAN EVALUATIONS.

!       INFO= 5   ITERATION IS NOT MAKING GOOD PROGRESS, AS MEASURED BY
!                  THE IMPROVEMENT FROM THE LAST TEN ITERATIONS.

!     NFEV IS AN INTEGER OUTPUT VARIABLE SET TO THE NUMBER OF CALLS TO FCN.

!     FJAC IS AN OUTPUT N BY N ARRAY WHICH CONTAINS THE ORTHOGONAL MATRIX Q
!       PRODUCED BY THE QR FACTORIZATION OF THE FINAL APPROXIMATE JACOBIAN.

!     LDFJAC IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN N
!       WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY FJAC.

!     R IS AN OUTPUT ARRAY OF LENGTH LR WHICH CONTAINS THE
!       UPPER TRIANGULAR MATRIX PRODUCED BY THE QR FACTORIZATION
!       OF THE FINAL APPROXIMATE JACOBIAN, STORED ROWWISE.

!     LR IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN (N*(N+1))/2.

!     QTF IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS
!       THE VECTOR (Q TRANSPOSE)*FVEC.

!     WA1, WA2, WA3, AND WA4 ARE WORK ARRAYS OF LENGTH N.

!   SUBPROGRAMS CALLED

!     USER-SUPPLIED ...... FCN

!     MINPACK-SUPPLIED ... DOGLEG,SPMPAR,ENORM,FDJAC1,
!                          QFORM,QRFAC,R1MPYQ,R1UPDT

!     FORTRAN-SUPPLIED ... ABS,MAX,MIN,MIN,MOD

!   ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!   BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE

!   **********
      
      INTEGER    :: i,iflag,iter,j,jm1,l,lr,msum,ncfail,ncsuc, &
     &              nslow1,nslow2
      INTEGER    :: iwa(1)
      LOGICAL    :: jeval,sing
      REAL*8  :: actred,delta,epsmch,fnorm,fnorm1,pnorm,prered,   &
     &           ratio,sum,temp,xnorm
      REAL*8, PARAMETER  :: one= 1.0d0,p1= 0.1d0,p5= 0.5d0,   &
     &                      p001= 0.001d0,p0001= 0.0001d0,zero= 0.0d0
      REAL*8  :: fjac(n,n),r(n*(n+1)/2),qtf(n),wa1(n),wa2(n), &
     &           wa3(n),wa4(n)
!
      epsmch= EPSILON(1.0d0)
      info= 0
      iflag= 0
      nfev= 0
      lr= n*(n+1)/2
      IF(n > 0.and.xtol >= zero.and.maxfev > 0.and.ml >= 0 &
     &        .and.mu >= 0.and.factor > zero ) THEN
      IF(mode == 2) THEN
       diag(1:n)= one
      ENDIF
      iflag= 1
      CALL VBFNLO_HIGGSTO_fcn(n,x,fvec,iflag)
      nfev= 1
      IF(iflag >= 0) THEN
      fnorm= VBFNLO_HIGGSTO_enorm(n,fvec)
      msum= MIN(ml+mu+1,n)
      iter= 1
      ncsuc= 0
      ncfail= 0
      nslow1= 0
      nslow2= 0
  20  jeval= .true.
      iflag= 2
      CALL VBFNLO_HIGGSTO_fdjac1(VBFNLO_HIGGSTO_FCN,n,x,fvec,fjac,n,iflag,ml,mu,epsfcn, &
     &                wa1,wa2)
      nfev= nfev+msum
      IF(iflag >= 0) THEN
      CALL VBFNLO_HIGGSTO_qrfac(n,n,fjac,n,.false.,iwa,1,wa1,wa2,wa3)
      IF(iter == 1) THEN
      IF(mode /= 2) THEN
       DO j= 1,n
        diag(j)= wa2(j)
        IF(wa2(j) == zero) diag(j)= one
       ENDDO
      ENDIF
      wa3(1:n)= diag(1:n)*x(1:n)
      xnorm= VBFNLO_HIGGSTO_enorm(n,wa3)
      delta= factor*xnorm
      IF(delta == zero) delta= factor
      ENDIF
      qtf(1:n)= fvec(1:n)
      DO j= 1,n
       IF(fjac(j,j) /= zero) THEN
        sum= zero
        DO i= j,n
          sum= sum+fjac(i,j)*qtf(i)
        ENDDO
        temp= -sum/fjac(j,j)
        DO i= j,n
          qtf(i)= qtf(i)+fjac(i,j)*temp
        ENDDO
       ENDIF
      ENDDO
      sing= .false.
      DO j= 1,n
       l= j
       jm1= j-1
       IF(jm1 >= 1) THEN
        DO i= 1,jm1
         r(l)= fjac(i,j)
         l= l+n-i
        ENDDO
       ENDIF
       r(l)= wa1(j)
       IF(wa1(j) == zero) sing= .true.
      ENDDO
      CALL VBFNLO_HIGGSTO_qform(n,n,fjac,n,wa1)
      IF(mode /= 2) THEN
       DO j= 1,n
        diag(j)= MAX(diag(j),wa2(j))
       ENDDO
      ENDIF
  120 IF(nprint > 0) THEN
       iflag= 0
       IF(MOD(iter-1,nprint) == 0) CALL VBFNLO_HIGGSTO_fcn(n,x,fvec,iflag)
       IF(iflag < 0) GO TO 190
      ENDIF
      CALL VBFNLO_HIGGSTO_dogleg(n,r,lr,diag,qtf,delta,wa1,wa2,wa3)
      DO j= 1,n
       wa1(j)= -wa1(j)
       wa2(j)= x(j)+wa1(j)
       wa3(j)= diag(j)*wa1(j)
      ENDDO
      pnorm= VBFNLO_HIGGSTO_enorm(n,wa3)
      IF(iter == 1) delta= MIN(delta,pnorm)
      iflag= 1
      CALL VBFNLO_HIGGSTO_fcn(n,wa2,wa4,iflag)
      nfev= nfev+1
      IF(iflag >= 0) THEN
      fnorm1= VBFNLO_HIGGSTO_enorm(n,wa4)
      actred= -one
      IF(fnorm1 < fnorm) actred= one-(fnorm1/fnorm)**2
      l= 1
      DO i= 1,n
       sum= zero
       DO j= i,n
        sum= sum+r(l)*wa1(j)
        l= l+1
       ENDDO
       wa3(i)= qtf(i)+sum
      ENDDO
      temp= VBFNLO_HIGGSTO_enorm(n,wa3)
      prered= zero
      IF(temp < fnorm) prered= one-(temp/fnorm)**2
      ratio= zero
      IF(prered > zero) ratio= actred/prered
      IF(ratio < p1) THEN
       ncsuc= 0
       ncfail= ncfail+1
       delta= p5*delta
      ELSE
       ncfail= 0
       ncsuc= ncsuc+1
       IF(ratio >= p5.or.ncsuc > 1) delta= MAX(delta,pnorm/p5)
       IF(ABS(ratio-one) <= p1) delta= pnorm/p5
      ENDIF
      IF(ratio >= p0001) THEN
       DO j= 1,n
        x(j)= wa2(j)
        wa2(j)= diag(j)*x(j)
        fvec(j)= wa4(j)
       ENDDO
       xnorm= VBFNLO_HIGGSTO_enorm(n,wa2)
       fnorm= fnorm1
       iter= iter+1
      ENDIF
      nslow1= nslow1+1
      IF(actred >= p001) nslow1= 0
      IF(jeval) nslow2= nslow2+1
      IF(actred >= p1) nslow2= 0
      IF(delta <= xtol*xnorm.or.fnorm == zero) info= 1
      IF(info == 0) THEN
      IF(nfev >= maxfev) info= 2
      IF(p1*MAX(p1*delta,pnorm) <= epsmch*xnorm) info= 3
      IF(nslow2 == 5) info= 4
      IF(nslow1 == 10) info= 5
      IF(info == 0) THEN
      IF(ncfail /= 2) THEN
      DO j= 1,n
       sum= zero
       DO i= 1,n
        sum= sum+fjac(i,j)*wa4(i)
       ENDDO
       wa2(j)= (sum-wa3(j))/pnorm
       wa1(j)= diag(j)*((diag(j)*wa1(j))/pnorm)
       IF(ratio >= p0001) qtf(j)= sum
      ENDDO
      CALL VBFNLO_HIGGSTO_r1updt(n,n,r,lr,wa1,wa2,wa3,sing)
      CALL VBFNLO_HIGGSTO_r1mpyq(n,n,fjac,n,wa2,wa3)
      CALL VBFNLO_HIGGSTO_r1mpyq(1,n,qtf,1,wa2,wa3)
      jeval= .false.
      GO TO 120
      ENDIF
      GO TO 20
      ENDIF
      ENDIF
      ENDIF
      ENDIF
      ENDIF
      ENDIF
!
 190  IF(iflag < 0) info= iflag
      iflag= 0
      IF(nprint > 0) CALL VBFNLO_HIGGSTO_fcn(n,x,fvec,iflag)
      RETURN
!
      END SUBROUTINE VBFNLO_HIGGSTO_hybrd
!
!--------------------------------------------------------------------
!
      SUBROUTINE VBFNLO_HIGGSTO_dogleg(n,r,lr,diag,qtb,delta,x,wa1,wa2)
!
      INTEGER, INTENT(IN)        :: n
      INTEGER, INTENT(IN)        :: lr
      REAL*8, INTENT(IN)      :: r(lr)
      REAL*8, INTENT(IN)      :: diag(n)
      REAL*8, INTENT(IN)      :: qtb(n)
      REAL*8, INTENT(IN)      :: delta
      REAL*8, INTENT(IN OUT)  :: x(n)
      REAL*8, INTENT(OUT)     :: wa1(n)
      REAL*8, INTENT(OUT)     :: wa2(n)

!     **********

!     SUBROUTINE VBFNLO_HIGGSTO_DOGLEG

!     GIVEN AN M BY N MATRIX A, AN N BY N NONSINGULAR DIAGONAL
!     MATRIX D, AN M-VECTOR B, AND A POSITIVE NUMBER DELTA, THE
!     PROBLEM IS TO DETERMINE THE CONVEX COMBINATION X OF THE
!     GAUSS-NEWTON AND SCALED GRADIENT DIRECTIONS THAT MINIMIZES
!     (A*X - B) IN THE LEAST SQUARES SENSE, SUBJECT TO THE
!     RESTRICTION THAT THE EUCLIDEAN NORM OF D*X BE AT MOST DELTA.

!     THIS SUBROUTINE VBFNLO_HIGGSTO_COMPLETES THE SOLUTION OF THE PROBLEM
!     IF IT IS PROVIDED WITH THE NECESSARY INFORMATION FROM THE
!     QR FACTORIZATION OF A. THAT IS, IF A= Q*R, WHERE Q HAS
!     ORTHOGONAL COLUMNS AND R IS AN UPPER TRIANGULAR MATRIX,
!     THEN DOGLEG EXPECTS THE FULL UPPER TRIANGLE OF R AND
!     THE FIRST N COMPONENTS OF (Q TRANSPOSE)*B.

!     THE SUBROUTINE VBFNLO_HIGGSTO_STATEMENT IS

!       SUBROUTINE VBFNLO_HIGGSTO_DOGLEG(N,R,LR,DIAG,QTB,DELTA,X,WA1,WA2)

!     WHERE

!       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE ORDER OF R.

!       R IS AN INPUT ARRAY OF LENGTH LR WHICH MUST CONTAIN THE UPPER
!         TRIANGULAR MATRIX R STORED BY ROWS.

!       LR IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN
!         (N*(N+1))/2.

!       DIAG IS AN INPUT ARRAY OF LENGTH N WHICH MUST CONTAIN THE
!         DIAGONAL ELEMENTS OF THE MATRIX D.

!       QTB IS AN INPUT ARRAY OF LENGTH N WHICH MUST CONTAIN THE FIRST
!         N ELEMENTS OF THE VECTOR (Q TRANSPOSE)*B.

!       DELTA IS A POSITIVE INPUT VARIABLE WHICH SPECIFIES AN UPPER
!         BOUND ON THE EUCLIDEAN NORM OF D*X.

!       X IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE DESIRED
!         CONVEX COMBINATION OF THE GAUSS-NEWTON DIRECTION AND THE
!         SCALED GRADIENT DIRECTION.

!       WA1 AND WA2 ARE WORK ARRAYS OF LENGTH N.

!     SUBPROGRAMS CALLED

!       MINPACK-SUPPLIED ... SPMPAR,ENORM

!       FORTRAN-SUPPLIED ... ABS,MAX,MIN,SQRT

!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE

!     **********

      INTEGER    :: i,j,jj,jp1,k,l
      REAL*8  :: alpha,bnorm,epsmch,gnorm,qnorm,sgnorm,sum,temp
!
      epsmch= EPSILON(1.0d0)
      jj= (n*(n+1))/2+1
      DO k= 1,n
       j= n-k+1
       jp1= j+1
       jj= jj-k
       l= jj+1
       sum= 0.0
       IF(n >= jp1) THEN
        DO i= jp1,n
         sum= sum+r(l)*x(i)
         l= l+1
        ENDDO
       ENDIF
       temp= r(jj)
       IF(temp == 0.0) THEN
        l= j
        DO i= 1,j
         temp= MAX(temp,ABS(r(l)))
         l= l+n-i
        ENDDO
        temp= epsmch*temp
        IF(temp == 0.0) temp= epsmch
       ENDIF
       x(j)= (qtb(j)-sum)/temp
      ENDDO
!
      DO j= 1,n
       wa1(j)= 0.0
       wa2(j)= diag(j)*x(j)
      ENDDO
      qnorm= VBFNLO_HIGGSTO_enorm(n,wa2)
!
      IF(qnorm > delta) THEN
       l= 1
       DO j= 1,n
        temp= qtb(j)
        DO i= j,n
         wa1(i)= wa1(i)+r(l)*temp
         l= l+1
        ENDDO
        wa1(j)= wa1(j)/diag(j)
       ENDDO
       gnorm= VBFNLO_HIGGSTO_enorm(n,wa1)
       sgnorm= 0.0
       alpha= delta/qnorm
       IF(gnorm /= 0.0) THEN
        DO j= 1,n
         wa1(j)= (wa1(j)/gnorm)/diag(j)
        ENDDO
        l= 1
        DO j= 1,n
         sum= 0.0
         DO i= j,n
          sum= sum+r(l)*wa1(i)
          l= l+1
         ENDDO
         wa2(j)= sum
        ENDDO
        temp= VBFNLO_HIGGSTO_enorm(n,wa2)
        sgnorm= (gnorm/temp)/temp
        alpha= 0.0
        IF(sgnorm < delta) THEN
         bnorm= VBFNLO_HIGGSTO_enorm(n,qtb)
         temp= (bnorm/gnorm)*(bnorm/qnorm)*(sgnorm/delta)
         temp= temp-(delta/qnorm)*(sgnorm/delta)**2+ &
     &         SQRT((temp-(delta/qnorm))**2+(1.0d0-(delta/qnorm)**2)* &
     &         (1.0d0-( sgnorm/delta)**2))
         alpha= ((delta/qnorm)*(1.0d0-(sgnorm/delta)**2))/temp
        ENDIF
       ENDIF
       temp= (1.0d0-alpha)*MIN(sgnorm,delta)
       DO j= 1,n
        x(j)= temp*wa1(j)+alpha*x(j)
       ENDDO
      ENDIF
!
      RETURN
!
      END SUBROUTINE VBFNLO_HIGGSTO_dogleg
!      
!-------------------------------------------------------------------
!
      SUBROUTINE VBFNLO_HIGGSTO_fdjac1(VBFNLO_HIGGSTO_FCN,n,x,fvec,fjac,ldfjac,iflag,ml,mu, &
     &                      epsfcn,wa1,wa2)
!
      INTEGER, INTENT(IN)        :: n
      REAL*8, INTENT(IN OUT)  :: x(n)
      REAL*8, INTENT(IN)      :: fvec(n)
      INTEGER, INTENT(IN)        :: ldfjac
      REAL*8, INTENT(OUT)     :: fjac(ldfjac,n)
      INTEGER, INTENT(IN OUT)    :: iflag
      INTEGER, INTENT(IN)        :: ml
      INTEGER, INTENT(IN)        :: mu
      REAL*8, INTENT(IN)      :: epsfcn
      REAL*8, INTENT(IN OUT)  :: wa1(n)
      REAL*8, INTENT(OUT)     :: wa2(n)
       
! EXTERNAL fcn
      
      INTERFACE
       SUBROUTINE VBFNLO_HIGGSTO_FCN(N,X,FVEC,IFLAG)
       IMPLICIT NONE
       INTEGER, PARAMETER  :: dp= SELECTED_REAL_KIND(14,60)
       INTEGER,INTENT(IN)      :: n
       REAL*8,INTENT(IN)    :: x(n)
       REAL*8, INTENT(OUT)   :: fvec(n)
       INTEGER, INTENT(IN OUT)  :: iflag
       END SUBROUTINE VBFNLO_HIGGSTO_FCN
      END INTERFACE

!   **********

!   SUBROUTINE VBFNLO_HIGGSTO_FDJAC1

!   THIS SUBROUTINE VBFNLO_HIGGSTO_COMPUTES A FORWARD-DIFFERENCE APPROXIMATION TO THE N BY N
!   JACOBIAN MATRIX ASSOCIATED WITH A SPECIFIED PROBLEM OF N FUNCTIONS IN N
!   VARIABLES.  IF THE JACOBIAN HAS A BANDED FORM, THEN FUNCTION EVALUATIONS
!   ARE SAVED BY ONLY APPROXIMATING THE NONZERO TERMS.

!   THE SUBROUTINE VBFNLO_HIGGSTO_STATEMENT IS

!     SUBROUTINE VBFNLO_HIGGSTO_FDJAC1(FCN,N,X,FVEC,FJAC,LDFJAC,IFLAG,ML,MU,EPSFCN,
!                       WA1,WA2)

!   WHERE

!     FCN IS THE NAME OF THE USER-SUPPLIED SUBROUTINE VBFNLO_HIGGSTO_WHICH CALCULATES
!       THE FUNCTIONS.  FCN MUST BE DECLARED IN AN EXTERNAL STATEMENT IN
!       THE USER CALLING PROGRAM, AND SHOULD BE WRITTEN AS FOLLOWS.

!       SUBROUTINE VBFNLO_HIGGSTO_FCN(N,X,FVEC,IFLAG)
!       INTEGER N,IFLAG
!       REAL X(N),FVEC(N)
!       ----------
!       CALCULATE THE FUNCTIONS AT X AND
!       RETURN THIS VECTOR IN FVEC.
!       ----------
!       RETURN
!       END

!       THE VALUE OF IFLAG SHOULD NOT BE CHANGED BY FCN UNLESS
!       THE USER WANTS TO TERMINATE EXECUTION OF FDJAC1.
!       IN THIS CASE SET IFLAG TO A NEGATIVE INTEGER.

!     N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!       OF FUNCTIONS AND VARIABLES.

!     X IS AN INPUT ARRAY OF LENGTH N.

!     FVEC IS AN INPUT ARRAY OF LENGTH N WHICH MUST CONTAIN THE
!       FUNCTIONS EVALUATED AT X.

!     FJAC IS AN OUTPUT N BY N ARRAY WHICH CONTAINS THE
!       APPROXIMATION TO THE JACOBIAN MATRIX EVALUATED AT X.

!     LDFJAC IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN N
!       WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY FJAC.

!     IFLAG IS AN INTEGER VARIABLE WHICH CAN BE USED TO TERMINATE
!       THE EXECUTION OF FDJAC1.  SEE DESCRIPTION OF FCN.

!     ML IS A NONNEGATIVE INTEGER INPUT VARIABLE WHICH SPECIFIES
!       THE NUMBER OF SUBDIAGONALS WITHIN THE BAND OF THE
!       JACOBIAN MATRIX. IF THE JACOBIAN IS NOT BANDED, SET
!       ML TO AT LEAST N - 1.

!     EPSFCN IS AN INPUT VARIABLE USED IN DETERMINING A SUITABLE
!       STEP LENGTH FOR THE FORWARD-DIFFERENCE APPROXIMATION. THIS
!       APPROXIMATION ASSUMES THAT THE RELATIVE ERRORS IN THE
!       FUNCTIONS ARE OF THE ORDER OF EPSFCN. IF EPSFCN IS LESS
!       THAN THE MACHINE PRECISION, IT IS ASSUMED THAT THE RELATIVE
!       ERRORS IN THE FUNCTIONS ARE OF THE ORDER OF THE MACHINE PRECISION.

!     MU IS A NONNEGATIVE INTEGER INPUT VARIABLE WHICH SPECIFIES
!       THE NUMBER OF SUPERDIAGONALS WITHIN THE BAND OF THE
!       JACOBIAN MATRIX. IF THE JACOBIAN IS NOT BANDED, SET
!       MU TO AT LEAST N - 1.

!     WA1 AND WA2 ARE WORK ARRAYS OF LENGTH N.  IF ML + MU + 1 IS AT
!       LEAST N, THEN THE JACOBIAN IS CONSIDERED DENSE, AND WA2 IS
!       NOT REFERENCED.

!   SUBPROGRAMS CALLED

!     MINPACK-SUPPLIED ... SPMPAR

!     FORTRAN-SUPPLIED ... ABS,MAX,SQRT

!   ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!   BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE

!   **********
!
      
      INTEGER    :: i,j,k,msum
      REAL*8  :: eps,epsmch,h,temp
      REAL*8, PARAMETER  :: zero= 0.0d0
!
      epsmch= EPSILON(1.0d0)
      eps= SQRT(MAX(epsfcn,epsmch))
      msum= ml+mu+1
      IF(msum >= n) THEN
       DO j= 1,n
        temp= x(j)
        h= eps*ABS(temp)
        IF(h == zero) h= eps
        x(j)= temp+h
        CALL VBFNLO_HIGGSTO_fcn(n,x,wa1,iflag)
        IF(iflag < 0) EXIT
        x(j)= temp
        DO i= 1,n
         fjac(i,j)= (wa1(i)-fvec(i))/h
        ENDDO
       ENDDO
      ELSE
       DO k= 1,msum
        DO j= k,n,msum
         wa2(j)= x(j)
         h= eps*ABS(wa2(j))
         IF(h == zero) h= eps
         x(j)= wa2(j)+h
        ENDDO
        CALL VBFNLO_HIGGSTO_fcn(n,x,wa1,iflag)
        IF(iflag < 0) EXIT
        DO j= k,n,msum
         x(j)= wa2(j)
         h= eps*ABS(wa2(j))
         IF(h == zero) h= eps
         DO i= 1,n
          fjac(i,j)= zero
          IF(i >= j-mu.and.i <= j+ml) fjac(i,j)= (wa1(i)-fvec(i))/h
         ENDDO
        ENDDO
       ENDDO
      ENDIF
! 
      RETURN
!
      END SUBROUTINE VBFNLO_HIGGSTO_fdjac1
!
!-------------------------------------------------------------------
!
      SUBROUTINE VBFNLO_HIGGSTO_qform(m,n,q,ldq,wa)
!
      INTEGER, INTENT(IN)     :: m
      INTEGER, INTENT(IN)     :: n
      INTEGER, INTENT(IN)     :: ldq
      REAL*8, INTENT(OUT)  :: q(ldq,m)
      REAL*8, INTENT(OUT)  :: wa(m)

!   **********

!   SUBROUTINE VBFNLO_HIGGSTO_QFORM

!   THIS SUBROUTINE VBFNLO_HIGGSTO_PROCEEDS FROM THE COMPUTED QR FACTORIZATION OF AN M BY N
!   MATRIX A TO ACCUMULATE THE M BY M ORTHOGONAL MATRIX Q FROM ITS FACTORED FORM.

!   THE SUBROUTINE VBFNLO_HIGGSTO_STATEMENT IS

!     SUBROUTINE VBFNLO_HIGGSTO_QFORM(M,N,Q,LDQ,WA)

!   WHERE

!     M IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!       OF ROWS OF A AND THE ORDER OF Q.

!     N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER OF COLUMNS OF A.

!     Q IS AN M BY M ARRAY. ON INPUT THE FULL LOWER TRAPEZOID IN
!       THE FIRST MIN(M,N) COLUMNS OF Q CONTAINS THE FACTORED FORM.
!       ON OUTPUT Q HAS BEEN ACCUMULATED INTO A SQUARE MATRIX.

!     LDQ IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN M
!       WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY Q.

!     WA IS A WORK ARRAY OF LENGTH M.

!   SUBPROGRAMS CALLED

!     FORTRAN-SUPPLIED ... MIN

!   ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!   BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE

!   **********
!
      INTEGER    :: i,j,jm1,k,l,minmn,np1
      REAL*8  :: sum,temp
      REAL*8, PARAMETER  :: one= 1.0d0,zero= 0.0d0
!
      minmn= MIN(m,n)
      IF(minmn >= 2) THEN
       DO j= 2,minmn
        jm1= j-1
        DO i= 1,jm1
         q(i,j)= zero
        ENDDO
       ENDDO
      ENDIF
      np1= n+1
      IF(m >= np1) THEN
       DO j= np1,m
        DO i= 1,m
         q(i,j)= zero
        ENDDO
        q(j,j)= one
       ENDDO
      ENDIF
      DO l= 1,minmn
       k= minmn-l+1
       DO i= k,m
        wa(i)= q(i,k)
        q(i,k)= zero
       ENDDO
       q(k,k)= one
       IF(wa(k) /= zero) THEN
        DO j= k,m
         sum= zero
         DO i= k,m
          sum= sum+q(i,j)*wa(i)
         ENDDO
         temp= sum/wa(k)
         DO i= k,m
          q(i,j)= q(i,j)-temp*wa(i)
         ENDDO
        ENDDO
       ENDIF
      ENDDO
! 
      RETURN
!
      END SUBROUTINE VBFNLO_HIGGSTO_qform
!
!----------------------------------------------------------------
!
      SUBROUTINE VBFNLO_HIGGSTO_qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
!
      INTEGER, INTENT(IN)        :: m
      INTEGER, INTENT(IN)        :: n
      INTEGER, INTENT(IN)        :: lda
      REAL*8, INTENT(IN OUT)  :: a(lda,n)
      LOGICAL, INTENT(IN)        :: pivot
      INTEGER, INTENT(IN)        :: lipvt
      INTEGER, INTENT(OUT)       :: ipvt(lipvt)
      REAL*8, INTENT(OUT)     :: rdiag(n)
      REAL*8, INTENT(OUT)     :: acnorm(n)
      REAL*8, INTENT(OUT)     :: wa(n)

!   **********

!   SUBROUTINE VBFNLO_HIGGSTO_QRFAC

!   THIS SUBROUTINE VBFNLO_HIGGSTO_USES HOUSEHOLDER TRANSFORMATIONS WITH COLUMN PIVOTING
!   (OPTIONAL) TO COMPUTE A QR FACTORIZATION OF THE M BY N MATRIX A.
!   THAT IS, QRFAC DETERMINES AN ORTHOGONAL MATRIX Q, A PERMUTATION MATRIX P,
!   AND AN UPPER TRAPEZOIDAL MATRIX R WITH DIAGONAL ELEMENTS OF NONINCREASING
!   MAGNITUDE, SUCH THAT A*P= Q*R.  THE HOUSEHOLDER TRANSFORMATION FOR
!   COLUMN K, K= 1,2,...,MIN(M,N), IS OF THE FORM

!                         T
!         I - (1/U(K))*U*U

!   WHERE U HAS ZEROS IN THE FIRST K-1 POSITIONS.  THE FORM OF THIS
!   TRANSFORMATION AND THE METHOD OF PIVOTING FIRST APPEARED IN THE
!   CORRESPONDING LINPACK SUBROUTINE.

!   THE SUBROUTINE VBFNLO_HIGGSTO_STATEMENT IS

!     SUBROUTINE VBFNLO_HIGGSTO_QRFAC(M,N,A,LDA,PIVOT,IPVT,LIPVT,RDIAG,ACNORM,WA)

!   WHERE

!     M IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER OF ROWS OF A.

!     N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!       OF COLUMNS OF A.

!     A IS AN M BY N ARRAY.  ON INPUT A CONTAINS THE MATRIX FOR WHICH THE
!       QR FACTORIZATION IS TO BE COMPUTED.  ON OUTPUT THE STRICT UPPER
!       TRAPEZOIDAL PART OF A CONTAINS THE STRICT UPPER TRAPEZOIDAL PART OF R,
!       AND THE LOWER TRAPEZOIDAL PART OF A CONTAINS A FACTORED FORM OF Q
!       (THE NON-TRIVIAL ELEMENTS OF THE U VECTORS DESCRIBED ABOVE).

!     LDA IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN M
!       WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY A.

!     PIVOT IS A LOGICAL INPUT VARIABLE.  IF PIVOT IS SET TRUE,
!       THEN COLUMN PIVOTING IS ENFORCED.  IF PIVOT IS SET FALSE,
!       THEN NO COLUMN PIVOTING IS DONE.

!     IPVT IS AN INTEGER OUTPUT ARRAY OF LENGTH LIPVT.  IPVT DEFINES THE
!       PERMUTATION MATRIX P SUCH THAT A*P= Q*R.
!       COLUMN J OF P IS COLUMN IPVT(J) OF THE IDENTITY MATRIX.
!       IF PIVOT IS FALSE, IPVT IS NOT REFERENCED.

!     LIPVT IS A POSITIVE INTEGER INPUT VARIABLE.  IF PIVOT IS FALSE,
!       THEN LIPVT MAY BE AS SMALL AS 1.  IF PIVOT IS TRUE, THEN
!       LIPVT MUST BE AT LEAST N.

!     RDIAG IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE
!       DIAGONAL ELEMENTS OF R.

!     ACNORM IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE NORMS OF
!       THE CORRESPONDING COLUMNS OF THE INPUT MATRIX A.
!       IF THIS INFORMATION IS NOT NEEDED, THEN ACNORM CAN COINCIDE WITH RDIAG.

!     WA IS A WORK ARRAY OF LENGTH N. IF PIVOT IS FALSE, THEN WA
!       CAN COINCIDE WITH RDIAG.

!   SUBPROGRAMS CALLED

!     MINPACK-SUPPLIED ... SPMPAR,ENORM

!     FORTRAN-SUPPLIED ... MAX,SQRT,MIN

!   ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!   BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE

!   **********
      INTEGER    :: i,j,jp1,k,kmax,minmn
      REAL*8  :: ajnorm,epsmch,sum,temp
      REAL*8, PARAMETER  :: one= 1.0d0,p05= 0.05d0,zero= 0.0d0
!
      epsmch= EPSILON(1.0d0)
      DO j= 1,n
       acnorm(j)= VBFNLO_HIGGSTO_enorm(m,a(1:,j))
       rdiag(j)= acnorm(j)
       wa(j)= rdiag(j)
       IF(pivot) ipvt(j)= j
      ENDDO
      minmn= MIN(m,n)
      DO j= 1,minmn
       IF(pivot) THEN
        kmax= j
        DO k= j,n
         IF(rdiag(k) > rdiag(kmax)) kmax= k
        ENDDO
        IF(kmax /= j) THEN
         DO i= 1,m
          temp= a(i,j)
          a(i,j)= a(i,kmax)
          a(i,kmax)= temp
         ENDDO
         rdiag(kmax)= rdiag(j)
         wa(kmax)= wa(j)
         k= ipvt(j)
         ipvt(j)= ipvt(kmax)
         ipvt(kmax)= k
        ENDIF
       ENDIF
       ajnorm= VBFNLO_HIGGSTO_enorm(m-j+1,a(j:,j))
       IF(ajnorm /= zero) THEN
        IF(a(j,j) < zero) ajnorm= -ajnorm
        DO i= j,m
         a(i,j)= a(i,j)/ajnorm
        ENDDO
        a(j,j)= a(j,j)+one
        jp1= j+1
        IF(n >= jp1) THEN
         DO k= jp1,n
          sum= zero
          DO i= j,m
           sum= sum+a(i,j)*a(i,k)
          ENDDO
          temp= sum/a(j,j)
          DO i= j,m
           a(i,k)= a(i,k)-temp*a(i,j)
          ENDDO
          IF(.NOT.(.NOT.pivot.OR.rdiag(k) == zero)) THEN
           temp= a(j,k)/rdiag(k)
           rdiag(k)= rdiag(k)*SQRT(MAX(zero,one-temp**2))
           IF(p05*(rdiag(k)/wa(k))**2 <= epsmch) THEN
            rdiag(k)= VBFNLO_HIGGSTO_enorm(m-j,a(jp1:,k))
            wa(k)= rdiag(k)
           ENDIF
          ENDIF
         ENDDO
        ENDIF
       ENDIF
       rdiag(j)= -ajnorm
      ENDDO
!
      RETURN
!
      END SUBROUTINE VBFNLO_HIGGSTO_qrfac
!
!-------------------------------------------------------------
!
      SUBROUTINE VBFNLO_HIGGSTO_r1mpyq(m,n,a,lda,v,w)
!
      INTEGER, INTENT(IN)        :: m
      INTEGER, INTENT(IN)        :: n
      INTEGER, INTENT(IN)        :: lda
      REAL*8, INTENT(IN OUT)  :: a(lda,n)
      REAL*8, INTENT(IN)      :: v(n)
      REAL*8, INTENT(IN)      :: w(n)

!   **********

!   SUBROUTINE VBFNLO_HIGGSTO_R1MPYQ

!   GIVEN AN M BY N MATRIX A, THIS SUBROUTINE VBFNLO_HIGGSTO_COMPUTES A*Q WHERE
!   Q IS THE PRODUCT OF 2*(N - 1) TRANSFORMATIONS

!         GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)

!   AND GV(I), GW(I) ARE GIVENS ROTATIONS IN THE (I,N) PLANE WHICH
!   ELIMINATE ELEMENTS IN THE I-TH AND N-TH PLANES, RESPECTIVELY.
!   Q ITSELF IS NOT GIVEN, RATHER THE INFORMATION TO RECOVER THE
!   GV, GW ROTATIONS IS SUPPLIED.

!   THE SUBROUTINE VBFNLO_HIGGSTO_STATEMENT IS

!     SUBROUTINE VBFNLO_HIGGSTO_R1MPYQ(M, N, A, LDA, V, W)

!   WHERE

!     M IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER OF ROWS OF A.

!     N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER OF COLUMNS OF A.

!     A IS AN M BY N ARRAY.  ON INPUT A MUST CONTAIN THE MATRIX TO BE
!       POSTMULTIPLIED BY THE ORTHOGONAL MATRIX Q DESCRIBED ABOVE.
!       ON OUTPUT A*Q HAS REPLACED A.

!     LDA IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN M
!       WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY A.

!     V IS AN INPUT ARRAY OF LENGTH N. V(I) MUST CONTAIN THE INFORMATION
!       NECESSARY TO RECOVER THE GIVENS ROTATION GV(I) DESCRIBED ABOVE.

!     W IS AN INPUT ARRAY OF LENGTH N. W(I) MUST CONTAIN THE INFORMATION
!       NECESSARY TO RECOVER THE GIVENS ROTATION GW(I) DESCRIBED ABOVE.

!   SUBROUTINES CALLED

!     FORTRAN-SUPPLIED ... ABS, SQRT

!   ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!   BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE

!   **********
      
      INTEGER    :: i,j,nmj,nm1
      REAL*8  :: COS,SIN,temp
      REAL*8, PARAMETER  :: one= 1.0d0
!
      nm1= n-1
      IF(nm1 >= 1) THEN
       DO nmj= 1,nm1
        j= n-nmj
        IF(ABS(v(j)) > one) COS= one/v(j)
        IF(ABS(v(j)) > one) SIN= SQRT(one-COS**2)
        IF(ABS(v(j)) <= one) SIN= v(j)
        IF(ABS(v(j)) <= one) COS= SQRT(one-SIN**2)
        DO i= 1,m
         temp= COS*a(i,j)-SIN*a(i,n)
         a(i,n)= SIN*a(i,j)+COS*a(i,n)
         a(i,j)= temp
        ENDDO
       ENDDO
       DO j= 1,nm1
        IF(ABS(w(j)) > one) COS= one/w(j)
        IF(ABS(w(j)) > one) SIN= SQRT(one-COS**2)
        IF(ABS(w(j)) <= one) SIN= w(j)
        IF(ABS(w(j)) <= one) COS= SQRT(one-SIN**2)
        DO i= 1,m
         temp= COS*a(i,j)+SIN*a(i,n)
         a(i,n)= -SIN*a(i,j)+COS*a(i,n)
         a(i,j)= temp
        ENDDO
       ENDDO
      ENDIF
!
      RETURN
!
      END SUBROUTINE VBFNLO_HIGGSTO_r1mpyq
!
!--------------------------------------------------------------
!
      SUBROUTINE VBFNLO_HIGGSTO_r1updt(m,n,s,ls,u,v,w,sing)
!
      INTEGER, INTENT(IN)        :: m
      INTEGER, INTENT(IN)        :: n
      INTEGER, INTENT(IN)        :: ls
      REAL*8, INTENT(IN OUT)  :: s(ls)
      REAL*8, INTENT(IN)      :: u(m)
      REAL*8, INTENT(IN OUT)  :: v(n)
      REAL*8, INTENT(OUT)     :: w(m)
      LOGICAL, INTENT(OUT)       :: sing

!   **********

!   SUBROUTINE VBFNLO_HIGGSTO_R1UPDT

!   GIVEN AN M BY N LOWER TRAPEZOIDAL MATRIX S, AN M-VECTOR U,
!   AND AN N-VECTOR V, THE PROBLEM IS TO DETERMINE AN
!   ORTHOGONAL MATRIX Q SUCH THAT

!                 T
!         (S + U*V )*Q

!   IS AGAIN LOWER TRAPEZOIDAL.

!   THIS SUBROUTINE VBFNLO_HIGGSTO_DETERMINES Q AS THE PRODUCT OF 2*(N - 1) TRANSFORMATIONS

!         GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)

!   WHERE GV(I), GW(I) ARE GIVENS ROTATIONS IN THE (I,N) PLANE
!   WHICH ELIMINATE ELEMENTS IN THE I-TH AND N-TH PLANES, RESPECTIVELY.
!   Q ITSELF IS NOT ACCUMULATED, RATHER THE INFORMATION TO RECOVER THE GV,
!   GW ROTATIONS IS RETURNED.

!   THE SUBROUTINE VBFNLO_HIGGSTO_STATEMENT IS

!     SUBROUTINE VBFNLO_HIGGSTO_R1UPDT(M,N,S,LS,U,V,W,SING)

!   WHERE

!     M IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER OF ROWS OF S.

!     N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!       OF COLUMNS OF S.  N MUST NOT EXCEED M.

!     S IS AN ARRAY OF LENGTH LS. ON INPUT S MUST CONTAIN THE LOWER
!       TRAPEZOIDAL MATRIX S STORED BY COLUMNS. ON OUTPUT S CONTAINS
!       THE LOWER TRAPEZOIDAL MATRIX PRODUCED AS DESCRIBED ABOVE.

!     LS IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN
!       (N*(2*M-N+1))/2.

!     U IS AN INPUT ARRAY OF LENGTH M WHICH MUST CONTAIN THE VECTOR U.

!     V IS AN ARRAY OF LENGTH N. ON INPUT V MUST CONTAIN THE VECTOR V.
!       ON OUTPUT V(I) CONTAINS THE INFORMATION NECESSARY TO
!       RECOVER THE GIVENS ROTATION GV(I) DESCRIBED ABOVE.

!     W IS AN OUTPUT ARRAY OF LENGTH M. W(I) CONTAINS INFORMATION
!       NECESSARY TO RECOVER THE GIVENS ROTATION GW(I) DESCRIBED ABOVE.

!     SING IS A LOGICAL OUTPUT VARIABLE.  SING IS SET TRUE IF ANY OF THE
!       DIAGONAL ELEMENTS OF THE OUTPUT S ARE ZERO.  OTHERWISE SING IS
!       SET FALSE.

!   SUBPROGRAMS CALLED

!     MINPACK-SUPPLIED ... SPMPAR

!     FORTRAN-SUPPLIED ... ABS,SQRT

!   ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!   BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE, JOHN L. NAZARETH

!   **********

      INTEGER    :: i,j,jj,l,nmj,nm1
      REAL*8  :: COS,cotan,giant,SIN,TAN,tau,temp
      REAL*8, PARAMETER  :: one= 1.0d0,p5= 0.5d0,p25= 0.25d0,zero= 0.0d0
!
      giant= HUGE(1.0d0)
      jj= (n*(2*m-n+1))/2-(m-n)
      l= jj
      DO i= n,m
       w(i)= s(l)
       l= l+1
      ENDDO
      nm1= n-1
      IF(nm1 >= 1) THEN
       DO nmj= 1,nm1
        j= n-nmj
        jj= jj-(m-j+1)
        w(j)= zero
        IF(v(j) /= zero) THEN
         IF(ABS(v(n)) < ABS(v(j))) THEN
          cotan= v(n)/v(j)
          SIN= p5/SQRT(p25+p25*cotan**2)
          COS= SIN*cotan
          tau= one
          IF(ABS(COS)*giant > one) tau= one/COS
         ELSE
          TAN= v(j)/v(n)
          COS= p5/SQRT(p25+p25*TAN**2)
          SIN= COS*TAN
          tau= SIN
         ENDIF
         v(n)= SIN*v(j)+COS*v(n)
         v(j)= tau
         l= jj
         DO i= j,m
          temp= COS*s(l)-SIN*w(i)
          w(i)= SIN*s(l)+COS*w(i)
          s(l)= temp
          l= l+1
         ENDDO
        ENDIF
       ENDDO
      ENDIF
      DO i= 1,m
       w(i)= w(i)+v(n)*u(i)
      ENDDO
      sing= .false.
      IF(nm1 >= 1) THEN
       DO j= 1,nm1
        IF(w(j) /= zero) THEN
         IF(ABS(s(jj)) < ABS(w(j))) THEN
          cotan= s(jj)/w(j)
          SIN= p5/SQRT(p25+p25*cotan**2)
          COS= SIN*cotan
          tau= one
          IF(ABS(COS)*giant > one) tau= one/COS
         ELSE
          TAN= w(j)/s(jj)
          COS= p5/SQRT(p25+p25*TAN**2)
          SIN= COS*TAN
          tau= SIN
         ENDIF
         l= jj
         DO i= j,m
          temp= COS*s(l)+SIN*w(i)
          w(i)= -SIN*s(l)+COS*w(i)
          s(l)= temp
          l= l+1
         ENDDO
         w(j)= tau
        ENDIF
        IF(s(jj) == zero) sing= .true.
        jj= jj+(m-j+1)
       ENDDO
      ENDIF
      l= jj
      DO i= n,m
       s(l)= w(i)
       l= l+1
      ENDDO
      IF(s(jj) == zero) sing= .true.
!
      RETURN
!
      END SUBROUTINE VBFNLO_HIGGSTO_r1updt
!
!-----------------------------------------------------------------
!
      FUNCTION VBFNLO_HIGGSTO_enorm(n,x) RESULT(fn_val)
! 
      INTEGER, INTENT(IN)    :: n
      REAL*8, INTENT(IN)  :: x(n)
      REAL*8              :: fn_val

!   **********

!   FUNCTION ENORM

!   GIVEN AN N-VECTOR X, THIS FUNCTION CALCULATES THE EUCLIDEAN NORM OF X.

!   THE EUCLIDEAN NORM IS COMPUTED BY ACCUMULATING THE SUM OF SQUARES IN THREE
!   DIFFERENT SUMS.  THE SUMS OF SQUARES FOR THE SMALL AND LARGE COMPONENTS
!   ARE SCALED SO THAT NO OVERFLOWS OCCUR.  NON-DESTRUCTIVE UNDERFLOWS ARE
!   PERMITTED.  UNDERFLOWS AND OVERFLOWS DO NOT OCCUR IN THE COMPUTATION OF THE UNSCALED
!   SUM OF SQUARES FOR THE INTERMEDIATE COMPONENTS.
!   THE DEFINITIONS OF SMALL, INTERMEDIATE AND LARGE COMPONENTS DEPEND ON
!   TWO CONSTANTS, RDWARF AND RGIANT.  THE MAIN RESTRICTIONS ON THESE CONSTANTS
!   ARE THAT RDWARF**2 NOT UNDERFLOW AND RGIANT**2 NOT OVERFLOW.
!   THE CONSTANTS GIVEN HERE ARE SUITABLE FOR EVERY KNOWN COMPUTER.

!   THE FUNCTION STATEMENT IS

!     REAL FUNCTION ENORM(N, X)

!   WHERE

!     N IS A POSITIVE INTEGER INPUT VARIABLE.

!     X IS AN INPUT ARRAY OF LENGTH N.

!   SUBPROGRAMS CALLED

!     FORTRAN-SUPPLIED ... ABS,SQRT

!   ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!   BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE

!   **********

      INTEGER    :: i
      REAL*8  :: agiant,floatn,s1,s2,s3,xabs,x1max,x3max
      REAL*8, PARAMETER  :: rdwarf= 1.0D-100,rgiant= 1.0D+100
!
      s1= 0.0d0
      s2= 0.0d0
      s3= 0.0d0
      x1max= 0.0d0
      x3max= 0.0d0
      floatn= n
      agiant= rgiant/floatn
      DO i= 1,n
       xabs= ABS(x(i))
       IF(xabs <= rdwarf.or.xabs >= agiant) THEN
        IF(xabs > rdwarf) THEN
         IF(xabs > x1max) THEN
          s1= 1.0d0+s1*(x1max/xabs)**2
          x1max= xabs
         ELSE
          s1= s1+(xabs/x1max)**2
         ENDIF
        ELSE
         IF(xabs > x3max) THEN
          s3= 1.0d0+s3*(x3max/xabs)**2
          x3max= xabs
         ELSE
          IF(xabs /= 0.0d0) s3= s3+(xabs/x3max)**2
         ENDIF
        ENDIF
       ELSE
        s2= s2+xabs**2
       ENDIF
      ENDDO
      IF(s1 /= 0.0d0) THEN
       fn_val= x1max*SQRT(s1+(s2/x1max)/x1max)
      ELSE
       IF(s2 /= 0.0d0) THEN
        IF(s2 >= x3max) fn_val= SQRT(s2*(1.0d0+(x3max/s2)*(x3max*s3)))
        IF(s2 < x3max) fn_val= SQRT(x3max*((s2/x3max)+(x3max*s3)))
       ELSE
        fn_val= x3max*SQRT(s3)
       ENDIF
      ENDIF
!
      RETURN
!
      END FUNCTION VBFNLO_HIGGSTO_enorm
!
      !END MODULE VBFNLO_HIGGSTO_Solve_NonLin
!
!------------------------------------------------------------------------------
!
       SUBROUTINE VBFNLO_HIGGSTO_poles(m,nv)   
       !USE VBFNLO_HIGGSTO_transfmh
       !USE VBFNLO_HIGGSTO_masses
       !USE VBFNLO_HIGGSTO_set_phys_const
       !USE VBFNLO_HIGGSTO_Solve_NonLin
!
       IMPLICIT NONE
!
       INTEGER nv,n,info
       REAL*8 m,tol  
       REAL*8, dimension(3) :: xcp,fvcp,diag
       ! EXTERNAL VBFNLO_HIGGSTO_Cpoles
!
       n= nv
       tol= 1.d-8
       tmuh= m
       xcp(1)= sqrt(20.d0/mw)
       xcp(2)= 1.d0
       xcp(3)= 1.d0
       CALL VBFNLO_HIGGSTO_HBRD(VBFNLO_HIGGSTO_Cpoles,n,xcp,fvcp,tol,tol,info,diag)
       print 2003,info
       print*,'   '
       IF(info == 1) THEN
        print*,'ALGORITHM ESTIMATES THAT THE RELATIVE ERROR'
        print*,'BETWEEN X AND THE SOLUTION IS AT MOST TOL'
       ELSEIF(info == 2) THEN
        print*,'NUMBER OF CALLS TO FCN HAS REACHED'
        print*,'OR EXCEEDED 200*(N+1)'
       ELSEIF(info == 3) THEN
        print*,'TOL IS TOO SMALL. NO FURTHER IMPROVEMENT IN'
        print*,'THE APPROXIMATE SOLUTION X IS POSSIBLE'
       ELSEIF(info == 4) THEN 
        print*,'ITERATION IS NOT MAKING GOOD PROGRESS'
       ENDIF
       print*,'   '
       print 2005,1.d-2*xcp(2)*xcp(2)*m
       print 2004,1.d-1*xcp(1)*xcp(1)*mw, &
     &            imw*(1.d0-0.5d0*(imw*imw)/(mw*mw))
       IF(n == 3) THEN
        print 2006,1.d-2*xcp(3)*xcp(3)*mt,imt
       ENDIF
!
 2003  format(' info =',i2)
!
 2004   format(' gammaW = ',2e20.5)
 2005   format(' gammaH = ',e20.5)
 2006   format(' gammat = ',2e20.5)
!
        RETURN
!
       END SUBROUTINE VBFNLO_HIGGSTO_poles   
!
!------------------------------------------------------------------------------
!
       SUBROUTINE VBFNLO_HIGGSTO_Cpoles(n,xcp,fvcp,iflag)   
       !USE VBFNLO_HIGGSTO_ferbos
       !USE VBFNLO_HIGGSTO_aux_Hcp
       !USE VBFNLO_HIGGSTO_units
       !USE VBFNLO_HIGGSTO_acmplx_pro
       !USE VBFNLO_HIGGSTO_acmplx_rat
       !USE VBFNLO_HIGGSTO_full_ln
       !USE VBFNLO_HIGGSTO_ln_2_riemann
       !USE VBFNLO_HIGGSTO_masses
       !USE VBFNLO_HIGGSTO_set_phys_const
       !USE VBFNLO_HIGGSTO_riemann
       !USE VBFNLO_HIGGSTO_optcp
       !USE VBFNLO_HIGGSTO_transfmh
       !USE VBFNLO_HIGGSTO_QCD
!
       IMPLICIT NONE
!

       INTEGER, INTENT(IN)      :: n
       REAL*8, INTENT(IN)    :: xcp(n)
       REAL*8, INTENT(OUT)   :: fvcp(n)
       INTEGER, INTENT(IN OUT)  :: iflag

       INTEGER iz
       REAL*8 muh,rgh,muhs,scal,scals,p2,xm0,rgw, &
     &        lswr,lswi,emc,emb,emt,asmur,rgt,as_NLO,EWC, &
     &        crmbs,crmcs,lcxb,lcxc,lcxbs,lcxcs,lclxb,lclxc,as_LO
       REAL*8, dimension(2,2) :: axms
       REAL*8, dimension(2) :: axm0,runbc
       REAL*8, dimension(2) :: sh,shs,clh,b0sumb,b0sumf,cxp,ksumb, &
     &         ksumf,coefB1,coefB2,coefB3,coefB4,coefB5,coefB6,coefB7, &
     &         coefB8,coefB9,coefB10,coefB11,coefB12, &
     &         shhf,shht,shhb,shh,b0sumt,ksumt,xms, &
     &         b0part,cpt,cpts,ccxt,deltag,sww0,coefW1, &
     &         coefW2,ksumw,ccxts,cctsi,shi,clt,cxhw,cstsi,csmcts,cltmw, &
     &         b0sumw,sww,DW,b0sumw1,b0sumw2,cxw,cxz,ccts,clw,csts, &
     &         cctq,cxws,cmxw,clmw,cxtw, &
     &         coeft1,coeft2,coeft3,coeft4,b0sumtop,ksumtop,stt, &
     &         coeft1s,b0sumtops,ksumtops,stts,cctqi,ucpt,nloqcd,ttqcd, &
     &         DWW,sww0W,swwW,ksumwW,lcxwi,lclcts
!     
       ! INTERFACE
       !  SUBROUTINE VBFNLO_HIGGSTO_INITALPHAS(asord,FR2,MUR,asmur,emc,emb,emt)
       !  !USE VBFNLO_HIGGSTO_DZpar 
       !  IMPLICIT NONE
       !  INTEGER asord
       !  REAL*8 FR2,MUR,asmur,emc,emb,emt,VBFNLO_HIGGSTO_FINDALPHASR0
       !  EXTERNAL VBFNLO_HIGGSTO_FINDALPHASR0
       !  END SUBROUTINE VBFNLO_HIGGSTO_INITALPHAS
       ! END INTERFACE
! !
       ! INTERFACE
       !  FUNCTION VBFNLO_HIGGSTO_ALPHAS(MUR)
       !  !USE VBFNLO_HIGGSTO_NFFIX  
       !  !USE VBFNLO_HIGGSTO_VARFLV 
       !  !USE VBFNLO_HIGGSTO_FRRAT  
       !  !USE VBFNLO_HIGGSTO_ASINP  
       !  !USE VBFNLO_HIGGSTO_ASFTHR 
       !  IMPLICIT NONE
       !  REAL*8 MUR,VBFNLO_HIGGSTO_ALPHAS
       !  END FUNCTION VBFNLO_HIGGSTO_ALPHAS
       ! END INTERFACE
! !
       ! INTERFACE
       !  FUNCTION VBFNLO_HIGGSTO_lquarkQCD(scal,psi,ps0i,xmsi,xm0i,type)  &
     ! &                         RESULT(value)
       !  !USE VBFNLO_HIGGSTO_riemann
       !  !USE VBFNLO_HIGGSTO_acmplx_pro
       !  !USE VBFNLO_HIGGSTO_acmplx_rat
       !  !USE VBFNLO_HIGGSTO_cmplx_root
       !  !USE VBFNLO_HIGGSTO_cmplx_rootz
       !  !USE VBFNLO_HIGGSTO_cmplx_srs_root
       !  !USE VBFNLO_HIGGSTO_ln_2_riemann
       !  !USE VBFNLO_HIGGSTO_full_ln
       !  !USE VBFNLO_HIGGSTO_sp_fun
       !  !USE VBFNLO_HIGGSTO_units
       !  IMPLICIT NONE
       !  INTEGER type
       !  REAL*8 scal,ps0i,xm0i
       !  REAL*8, dimension(2) :: value,psi,xmsi
       !  END FUNCTION VBFNLO_HIGGSTO_lquarkQCD
       ! END INTERFACE
! !
       ! INTERFACE
       !  FUNCTION VBFNLO_HIGGSTO_lb0af_em(scal,psi,ps0i,xmsi,xm0i) RESULT(value)
       !  !USE VBFNLO_HIGGSTO_acmplx_pro
       !  !USE VBFNLO_HIGGSTO_acmplx_rat
       !  !USE VBFNLO_HIGGSTO_cmplx_root
       !  !USE VBFNLO_HIGGSTO_cmplx_rootz
       !  !USE VBFNLO_HIGGSTO_cmplx_srs_root
       !  !USE VBFNLO_HIGGSTO_ln_2_riemann
       !  !USE VBFNLO_HIGGSTO_full_ln
       !  !USE VBFNLO_HIGGSTO_units
       !  IMPLICIT NONE
       !  REAL*8 scal,ps0i,xm0i
       !  REAL*8, dimension(2) :: value,psi,xmsi
       !  END FUNCTION VBFNLO_HIGGSTO_lb0af_em
       ! END INTERFACE
! !
       ! INTERFACE
       !  FUNCTION VBFNLO_HIGGSTO_lb0af_dm(scal,psi,ps0i,xmsi,xm0i) RESULT(value)
       !  !USE VBFNLO_HIGGSTO_acmplx_pro
       !  !USE VBFNLO_HIGGSTO_acmplx_rat
       !  !USE VBFNLO_HIGGSTO_cmplx_root
       !  !USE VBFNLO_HIGGSTO_cmplx_rootz
       !  !USE VBFNLO_HIGGSTO_cmplx_srs_root
       !  !USE VBFNLO_HIGGSTO_ln_2_riemann
       !  !USE VBFNLO_HIGGSTO_full_ln
       !  !USE VBFNLO_HIGGSTO_units
       !  IMPLICIT NONE
       !  REAL*8 scal,ps0i
       !  REAL*8, dimension(2,2) :: xmsi
       !  REAL*8, dimension(2) :: value,psi,xm0i
       !  END FUNCTION VBFNLO_HIGGSTO_lb0af_dm
       ! END INTERFACE
! !
       ! INTERFACE
       !  FUNCTION VBFNLO_HIGGSTO_lb021_dm_cp(scal,psi,ps0i,xmsi,xm0i) RESULT(value)
       !  !USE VBFNLO_HIGGSTO_cmplx_root
       !  !USE VBFNLO_HIGGSTO_cmplx_rootz
       !  !USE VBFNLO_HIGGSTO_cmplx_srs_root
       !  !USE VBFNLO_HIGGSTO_ln_2_riemann
       !  !USE VBFNLO_HIGGSTO_acmplx_pro
       !  !USE VBFNLO_HIGGSTO_acmplx_rat
       !  !USE VBFNLO_HIGGSTO_full_ln
       !  !USE VBFNLO_HIGGSTO_units
       !  IMPLICIT NONE
       !  REAL*8 scal,ps0i
       !  REAL*8, intent(in), dimension(2) :: xm0i
       !  REAL*8, intent(in), dimension(2,2) :: xmsi
       !  REAL*8, dimension(2) :: value
       !  REAL*8, dimension(2) :: psi
       !  END FUNCTION VBFNLO_HIGGSTO_lb021_dm_cp
       ! END INTERFACE
!
       muh= tmuh
       scal= muh
       muhs= muh*muh
       scals= scal*scal
!
       cxw(1)= mw*mw/scals*(1.d0-0.5d-2*(xcp(1)*xcp(1))**2)
       cxw(2)= -1.d-1*mw*mw/scals*xcp(1)*xcp(1)
!
       cmxw= -cxw
       cxws= cxw.cp.cxw
       lcxwi= co.cq.cxw
!
       rgw= 1.d-1*xcp(1)*xcp(1)*mw
       rgh= 1.d-2*xcp(2)*xcp(2)*muh
       IF(n == 3) THEN
        rgt= 1.d-2*xcp(3)*xcp(3)*mt
       ELSEIF(n == 2) THEN
        rgt= imt
       ENDIF
!      
       lswr= mw*mw*(1.d0-0.5d-2*(xcp(1)*xcp(1))**2)
       lswi= -mw*rgw
!
       clw= cxw(1).fln.cxw(2)
       clmw= cmxw(1).fln.cmxw(2)
       clmw(2)= clmw(2)-2.d0*pi
!
       cxz(1)= szr/scals
       cxz(2)= szi/scals
       ccts= cxw.cq.cxz
       csts= co-ccts
       cctq= ccts.cp.ccts
       cctqi= co.cq.cctq
       lclcts= ccts(1).fln.ccts(2)
!
       asmur= 0.12018d0
       emc= 1.4d0
       emb= 4.75d0
       emt= mt
       iz= 1
       CALL VBFNLO_HIGGSTO_INITALPHAS(iz,one,mz,asmur,emc,emb,emt)
       as_NLO= VBFNLO_HIGGSTO_ALPHAS(scal)/pi
       runbc= VBFNLO_HIGGSTO_run_bc(scal)
       crmbs= runbc(2)*runbc(2)
       crmcs= runbc(1)*runbc(1)         
!
       lcxb= crmbs/scals
       lcxc= crmcs/scals
!
       lcxbs= lcxb*lcxb
       lcxcs= lcxc*lcxc
       lclxb= log(lcxb)
       lclxc= log(lcxc)
!
       ucpt(1)= mt*mt
       ucpt(2)= 0.d0
       cpt(1)= mt*mt/scals
       cpt(2)= 0.d0
       cpts= cpt.cp.cpt
!       
       ccxt= cpt
       ccxts= ccxt.cp.ccxt
       cctsi= co.cq.ccts
       cstsi= co.cq.csts
       csmcts= csts-ccts
!
       sh(1)= muhs/scals
       sh(2)= -muh*rgh/scals
!
       cxhw= sh-cxw 
       shs= sh.cp.sh
       shi= co.cq.sh
!
       cxtw= ccxt-cxw
       clh= sh(1).fln.sh(2)
       clt= ccxt(1).fln.ccxt(2)
       cltmw= cxtw(1).fln.cxtw(2)
!
! W
       coefB1= (12.d0*cxw-4.d0*sh+(shs.cq.cxw))/64.d0
!
! Z
       coefB2= (-4.d0*(sh.cq.ccts)+(shs.cq.cxw)+ &
     &         12.d0*(cxw.cq.cctq))/128.d0
!  
! H
       coefB3= 9.d0/128.d0*(shs.cq.cxw)
!
! top
       coefB4= -3.d0/32.d0*((4.d0*ccxt-sh).cp.(ccxt.cq.cxw))
!
! light fermions
!
       coefB5= -3.d0/32.d0*((4.d0*co*lcxb-sh).cq.cxw)*lcxb
!
       coefB6= -1.d0/32.d0*((4.d0*co*cxtau-sh).cq.cxw)*cxtau
!
       coefB7= -3.d0/32.d0*((4.d0*co*lcxc-sh).cq.cxw)*lcxc
!
       coefB8= -1.d0/32.d0*((4.d0*co*cxmu-sh).cq.cxw)*cxmu
!
       coefB9= -3.d0/32.d0*((4.d0*co*cxs-sh).cq.cxw)*cxs
!
       coefB10= -3.d0/32.d0*((4.d0*co*cxd-sh).cq.cxw)*cxd
!
       coefB11= -3.d0/32.d0*((4.d0*co*cxu-sh).cq.cxw)*cxu
!
       coefB12= -1.d0/32.d0*((4.d0*co*cxe-sh).cq.cxw)*cxe
!
       cxp(1)= muhs
       cxp(2)= -muh*rgh
       p2= muhs
!
       xms(1)= lswr
       xms(2)= lswi
       xm0= mw
       b0part= VBFNLO_HIGGSTO_lb0af_em(scal,cxp,p2,xms,xm0)
       b0sumb= coefB1.cp.b0part
!
       xms(1)= szr
       xms(2)= szi
       xm0= mz
       b0part= VBFNLO_HIGGSTO_lb0af_em(scal,cxp,p2,xms,xm0)
       b0sumb= b0sumb+(coefB2.cp.b0part)
!
       xms(1)= muhs
       xms(2)= -muh*rgh
       xm0= muh
       b0part= VBFNLO_HIGGSTO_lb0af_em(scal,cxp,p2,xms,xm0)
       b0sumb= b0sumb+(coefB3.cp.b0part)
!
       xms(1:2)= crmbs*co(1:2)
       xm0= sqrt(crmbs)
       b0part= VBFNLO_HIGGSTO_lb0af_em(scal,cxp,p2,xms,xm0)
       b0sumf= coefB5.cp.b0part
!
       xms(1:2)= mtl*mtl*co(1:2)
       xm0= mtl
       b0part= VBFNLO_HIGGSTO_lb0af_em(scal,cxp,p2,xms,xm0)
       b0sumf= b0sumf+(coefB6.cp.b0part)
!
       xms(1:2)= crmcs*co(1:2)
       xm0= sqrt(crmcs)
       b0part= VBFNLO_HIGGSTO_lb0af_em(scal,cxp,p2,xms,xm0)
       b0sumf= b0sumf+(coefB7.cp.b0part)
!
       xms(1:2)= mm*mm*co(1:2)
       xm0= mm
       b0part= VBFNLO_HIGGSTO_lb0af_em(scal,cxp,p2,xms,xm0)
       b0sumf= b0sumf+(coefB8.cp.b0part)
!
       xms(1:2)= msq*msq*co(1:2)
       xm0= msq
       b0part= VBFNLO_HIGGSTO_lb0af_em(scal,cxp,p2,xms,xm0)
       b0sumf= b0sumf+(coefB9.cp.b0part)
!
       xms(1:2)= mdq*mdq*co(1:2)
       xm0= mdq
       b0part= VBFNLO_HIGGSTO_lb0af_em(scal,cxp,p2,xms,xm0)
       b0sumf= b0sumf+(coefB10.cp.b0part)
!
       xms(1:2)= muq*muq*co(1:2)
       xm0= muq
       b0part=VBFNLO_HIGGSTO_lb0af_em(scal,cxp,p2,xms,xm0)
       b0sumf= b0sumf+(coefB11.cp.b0part)
!
       xms(1:2)= me*me*co(1:2)
       xm0= me
       b0part= VBFNLO_HIGGSTO_lb0af_em(scal,cxp,p2,xms,xm0)
       b0sumf= b0sumf+(coefB12.cp.b0part)
!
       xms(1:2)= ucpt(1:2)
       xm0= mt
       b0part= VBFNLO_HIGGSTO_lb0af_em(scal,cxp,p2,xms,xm0)
       b0sumt= coefB4.cp.b0part
!
       iz= 0
       xms(1:2)= crmbs*co(1:2)
       xm0= sqrt(crmbs)
       nloqcd= crmbs*VBFNLO_HIGGSTO_lquarkQCD(scal,cxp,p2,xms,xm0,iz)
       xms(1:2)= crmcs*co(1:2)
       xm0= sqrt(crmcs)
       nloqcd= nloqcd+crmcs*VBFNLO_HIGGSTO_lquarkQCD(scal,cxp,p2,xms,xm0,iz)
       IF(qcdc==0) nloqcd= 0.d0
!
       xms= ucpt
       xm0= mt
       iz= 1
       ttqcd= VBFNLO_HIGGSTO_lquarkQCD(scal,cxp,p2,xms,xm0,iz)
       ttqcd= cpt.cp.ttqcd
       IF(qcdc==0) ttqcd= 0.d0
!
       ksumf= -1.d0/32.d0*(sh.cq.cxw)*(cxe*clxe+cxmu*clxmu+ &
     &        cxtau*clxtau+3.d0*(cxd*clxd+cxs*clxs+lcxb*lclxb+ &
     &        cxu*clxu+lcxc*lclxc))+ &
     &        1.d0/8.d0*lcxwi*(cxes+cxmus+cxtaus+3.d0*( &
     &        cxds+cxss+lcxbs+cxus+lcxcs))
       ksumt= 3.d0/8.d0*((ccxt.cq.cxw).cp.ccxt) &
     &       -3.d0/32.d0*(clt.cp.((sh.cq.cxw).cp.ccxt))
!
       ksumb= -1.d0/64.d0*((2.d0*co+cctsi+3.d0*(sh.cq.cxw)).cp.sh) &
     &  -1.d0/128.d0*(lclcts.cp.((6.d0*cctsi-(sh.cq.cxw)).cp.sh)) &
     &  +3.d0/128.d0*(clw.cp.((4.d0*co+2.d0*cctsi-(sh.cq.cxw)).cp.sh)) &
     &  -3.d0/128.d0*(clh.cp.(shs.cq.cxw))
!
       shhf= b0sumf+ksumf
       shht= b0sumt+ksumt
       IF((lswr*shht(2)+lswi*shht(1)) < 0.d0) shht= 0.d0
       shhf= shhf+shht
       shhb= b0sumb+ksumb
       IF((lswr*shhb(2)+lswi*shhb(1)) < 0.d0) shhb= 0.d0
       shh= shhf+shhb
!
!--- W self energies
!
       deltag= 6.d0*co+0.5d0*(((7.d0*co-4.d0*csts).cq.csts).cp.lclcts)
!
       sww0= -(38.d0*cxw+6.d0*ccxt+7.d0*sh &
     &   -48.d0*(((ccxt.cq.sh).cq.cxw).cp.ccxt)+8.d0*(cxw.cq.sh))/128.d0 &
     &   -3.d0/64.d0*((cxw-sh+(cxws.cq.cxhw)).cp.clh) &
     &   +3.d0/32.d0*(((co-4.d0*((ccxt.cq.sh).cq.cxw).cp.ccxt)).cp.clt) &
     &   +((((8.d0*co-17.d0*cstsi+3.d0*cctsi).cp.cxw) &
     &   -6.d0*((cxw.cq.sh).cq.cctq)).cp.lclcts)/64.d0 &
     &    -((cxw.cq.sh).cq.cctq)/32.d0+5.d0/128.d0*(cxw.cq.ccts)
!
       sww0W= -(38.d0*cxw+6.d0*ccxt+7.d0*sh &
     &   -48.d0*(((ccxt.cq.sh).cq.cxw).cp.ccxt)+8.d0*(cxw.cq.sh))/128.d0 &
     &   -3.d0/64.d0*((cxw-sh+(cxws.cq.cxhw)).cp.(clh-clw)) &
     &   +3.d0/32.d0* &
     &    (((co-4.d0*((ccxt.cq.sh).cq.cxw).cp.ccxt)).cp.(clt-clw)) &
     &   +((((8.d0*co-17.d0*cstsi+3.d0*cctsi).cp.cxw) &
     &   -6.d0*((cxw.cq.sh).cq.cctq)).cp.lclcts)/64.d0 &
     &    -((cxw.cq.sh).cq.cctq)/32.d0+5.d0/128.d0*(cxw.cq.ccts)
!
       coefW1= -(((8.d0*co-(sh.cq.cxw)).cp.sh)*sh &
     &     -4.d0*((-12.d0*cxw+7.d0*sh).cp.cxw))/192.d0
!
       coefW2= -((cxws.cq.csmcts).cp.(416.d0*co-192.d0*csts  &
     &     -((132.d0*co-((12.d0*co+cctsi).cq.ccts)).cq.ccts)))/192.d0
!
       cxp(1)= lswr
       cxp(2)= lswi
       p2= mw*mw
!
       axms(1,1)= lswr
       axms(1,2)= lswi
       axm0(1)= mw
       axms(2,1)= muhs
       axms(2,2)= -muh*rgh
       axm0(2)= muh
       b0part= VBFNLO_HIGGSTO_lb021_dm_cp(scal,cxp,p2,axms,axm0)
       b0sumw1= (coefW1.cp.b0part)
       b0sumw= (coefW1.cp.b0part)
!
       axms(1,1)= szr
       axms(1,2)= szi
       axm0(1)= mz
       axms(2,1)= lswr
       axms(2,2)= lswi
       axm0(2)= mw
       b0part= VBFNLO_HIGGSTO_lb021_dm_cp(scal,cxp,p2,axms,axm0)
       b0sumw2= (coefW2.cp.b0part)
       b0sumw= b0sumw+(coefW2.cp.b0part)
!
       ksumw= -12.d0*((cxw &
     &    -0.5d0*((3.d0*co-(ccxts.cq.cxws)).cp.ccxt)).cp.cltmw) &
     &    -((24.d0*cxw-((14.d0*co-(sh.cq.cxw)).cp.sh)).cp.clh) &
     &    +((36.d0*cxw-14.d0*sh-18.d0*((co-4.d0*(ccxt.cq.sh)).cp.ccxt) &
     &    +(shs.cq.cxw)).cp.clw) &
     &    -6.d0*(((2.d0*co+((lcxwi-12.d0*shi).cp.ccxt)).cp.ccxt) &
     &    +1.d0/6.d0*((15.d0*co-(sh.cq.cxw)).cp.sh) &
     &    +2.d0/9.d0*((97.d0*co+9.d0*(cxw.cq.sh)).cp.cxw)) &
     &    +(((cxw.cq.ccts).cp.(co-6.d0*(cxw.cq.sh))).cq.ccts) &
     &    -2.d0*(((cxw.cq.csmcts).cp.lclcts).cp.(62.d0*co &
     &    -48.d0*csts-5.d0*cctsi)) &
     &    -18.d0*(((cxws.cq.sh).cq.cctq).cp.lclcts) &
     &    -72.d0*((clt.cp.ccxts).cp.(shi-1.d0/12.d0*(ccxt.cq.cxws))) &
     &    +23.d0*(cxw.cq.ccts)
       ksumw= ksumw/192.d0+3.d0/16.d0*(cxw.cp.(clw-clmw))
!
       ksumwW= -12.d0*((cxw &
     &    -0.5d0*((3.d0*co-(ccxts.cq.cxws)).cp.ccxt)).cp.(cltmw-clw)) &
     &    -((24.d0*cxw-((14.d0*co-(sh.cq.cxw)).cp.sh)).cp.(clh-clw)) &
     &    -6.d0*(((2.d0*co+((lcxwi-12.d0*shi).cp.ccxt)).cp.ccxt) &
     &    +1.d0/6.d0*((15.d0*co-(sh.cq.cxw)).cp.sh) &
     &    +2.d0/9.d0*((97.d0*co+9.d0*(cxw.cq.sh)).cp.cxw)) &
     &    +(((cxw.cq.ccts).cp.(co-6.d0*(cxw.cq.sh))).cq.ccts) &
     &    -2.d0*(((cxw.cq.csmcts).cp.lclcts).cp.(62.d0*co &
     &    -48.d0*csts-5.d0*cctsi)) &
     &    -18.d0*(((cxws.cq.sh).cq.cctq).cp.lclcts) &
     &    -72.d0* &
     &     (((clt-clw).cp.ccxts).cp.(shi-1.d0/12.d0*(ccxt.cq.cxws))) &
     &    +23.d0*(cxw.cq.ccts)
       ksumwW= ksumwW/192.d0+3.d0/16.d0*(cxw.cp.(clw-clmw))
!
       sww= b0sumw+ksumw
       swwW= b0sumw+ksumwW
!
       DW= -sww+sww0+deltag/16.d0
       DWW= -swwW+sww0W+deltag/16.d0
!       DW= 0.d0
!
       ksumtop= 1.d0/128.d0*(-(2.d0*cxw+3.d0*cxb*co+cxbs*lcxwi)+48.d0* &
     & ((cpt.cq.sh).cq.cxw)*cxbs+(cpts.cq.cxw)+(cpt.cp.(co &
     & -4.d0*cxb*lcxwi)))*clxb &
     & +1.d0/576.d0*(-(50.d0*cxw-9.d0*cxb*co-9.d0*cxbs*lcxwi &
     & +17.d0*(cxw.cq.cctq)-40.d0*(cxw.cq.ccts))+18.d0*(cpts.cq.cxw) &
     & +(cpt.cp.(co+17.d0*cctsi+9.d0*cxb*lcxwi-216.d0*(( &
     & cpts.cq.cxw).cq.sh)+54.d0*((cxw.cq.sh).cq.cctq) &
     & +18.d0*(sh.cq.cxw)-216.d0*(shi.cp.cxw)*cxbs+108.d0*(cxw.cq.sh)))) &
     & +1.d0/1152.d0*(clt.cp.(432.d0*((cpts.cq.cxw).cp.(cpt.cq.sh)) &
     & -(cxw.cp.(32.d0*co-40.d0*cctsi+17.d0*cctqi)) &
     & +(cpt.cp.(32.d0-64.d0*csts-41.d0*cctsi-9.d0*(sh.cq.cxw))))) &
     & -1.d0/128.d0*(clh.cp.((-4.d0*cpt+5.d0*sh).cp.(cpt.cq.cxw))) &
     & +1.d0/1152.d0*(clw.cp.((50.d0*cxw+27.d0*cxb*co+9.d0*cxbs* &
     & lcxwi+17.d0*(cxw.cq.cctq)-40.d0*(cxw.cq.ccts)) &
     & +9.d0*(cpts.cq.cxw)+(cpt.cp.(7.d0*co+64.d0*csts-7.d0*cctsi &
     & -18.d0*cxb*lcxwi-108.d0*((cxw.cq.sh).cq.cctq) &
     & -216.d0*(cxw.cq.sh))))) &
     & -1.d0/1152.d0*(lclcts.cp.((cxw.cp.(32.d0*co-40.d0* &
     & cctsi+17.d0*cctqi))+(cpt.cp.(16.d0*co+64.d0*csts-7.d0* &
     & cctsi-108.d0*((cxw.cq.sh).cq.cctq)))))
!
       coeft1= -1.d0/576.d0*(-(cxw.cp.(32.d0*co-40.d0*cctsi+17.d0* &
     &    cctqi))+(cpt.cp.(16.d0*co+64.d0*csts-7.d0*cctsi)))
!
       coeft2= 1.d0/64.d0*((-4.d0*cpt+sh).cp.(cpt.cq.cxw))
!
       coeft3= 1.d0/9.d0*(csts.cp.cpt)
!
       coeft4= 1.d0/64.d0*((-2.d0*cxw+cxb*co+cxbs*lcxwi) &
     & +(cpts.cq.cxw)+(cpt.cp.(co-2.d0*cxb*lcxwi)))
!
       cxp(1)= mt*mt
       cxp(2)= -mt*rgt
       p2= mt*mt
!
       axms(1,1)= mt*mt
       axms(1,2)= -mt*rgt
       axm0(1)= mt
       axms(2,1)= szr
       axms(2,2)= szi
       axm0(2)= mz
       b0part= VBFNLO_HIGGSTO_lb0af_dm(scal,cxp,p2,axms,axm0)
       b0sumtop= (coeft1.cp.b0part)
!
       axms(1,1)= mt*mt
       axms(1,2)= -mt*rgt
       axm0(1)= mt
       axms(2,1)= muhs
       axms(2,2)= -muh*rgh
       axm0(2)= muh
       b0part= VBFNLO_HIGGSTO_lb0af_dm(scal,cxp,p2,axms,axm0)
       b0sumtop= b0sumtop+(coeft2.cp.b0part)
!
       b0part= 2.d0*co
       b0sumtop= b0sumtop+(coeft3.cp.b0part)
!
       axms(1,1)= mb*mb
       axms(1,2)= 0.d0
       axm0(1)= mb
       axms(2,1)= lswr
       axms(2,2)= lswi
       axm0(2)= mw
       b0part= VBFNLO_HIGGSTO_lb0af_dm(scal,cxp,p2,axms,axm0)
       b0sumtop= b0sumtop+(coeft4.cp.b0part)
!
       ksumtops= cpt/6.d0-(clt.cp.cpt)/2.d0
!
       coeft1s= 1.d0/3.d0*cpt    
!
       b0part= 2.d0*co
       b0sumtops= coeft1s.cp.b0part
!
       stt= b0sumtop+ksumtop
       stts= b0sumtops+ksumtops
!
       EWC= 4.d0*sqrt(2.d0)*g_f/pis
!
       asmur= 0.13939d0
       emc= 1.4d0
       emb= 4.75d0
       emt= mt
       iz= 0
       CALL VBFNLO_HIGGSTO_INITALPHAS(iz,one,mz,asmur,emc,emb,emt)
       as_LO= VBFNLO_HIGGSTO_ALPHAS(scal)/pi
!
       fvcp(1)= rgw/mw*(1.d0+EWC*(lswr*DWW(1)-lswi*DWW(2))) &
     &       -EWC*(scal/mw)**2*(lswr*swwW(2)+lswi*swwW(1))
!
       fvcp(2)= rgh/muh*(1.d0+EWC*(lswr*DW(1)-lswi*DW(2))) &
     &       -EWC*(lswr*shh(2)+lswi*shh(1))*scals/muhs &
     &       +as_NLO*g_f/(sqrt(2.d0)*pis)* &
     &       (sh(1)*nloqcd(2)+sh(2)*nloqcd(1))*scals/muhs &
     &       +as_NLO*g_f/(sqrt(2.d0)*pis)* &
     &       (sh(1)*ttqcd(2)+sh(2)*ttqcd(1))*scals/muhs
!
       IF(n == 3) THEN
        fvcp(3)= rgt/mt*(1.d0+2.d0*EWC*(lswr*DW(1)-lswi*DW(2))) &
     &        -scals/(mt*mt)*(EWC*(lswr*stt(2)+lswi*stt(1))+ &
     &        4.d0*(scals/(mt*mt))*as_LO*stts(2))
       ENDIF
!
       RETURN
!
       END SUBROUTINE VBFNLO_HIGGSTO_Cpoles
!
!------------------------------------------------------------------
!
      FUNCTION VBFNLO_HIGGSTO_lquarkQCD(scal,psi,ps0i,xmsi,xm0i,type) RESULT(value)
      !USE VBFNLO_HIGGSTO_riemann
      !USE VBFNLO_HIGGSTO_acmplx_pro
      !USE VBFNLO_HIGGSTO_acmplx_rat
      !USE VBFNLO_HIGGSTO_cmplx_root
      !USE VBFNLO_HIGGSTO_cmplx_rootz
      !USE VBFNLO_HIGGSTO_cmplx_srs_root
      !USE VBFNLO_HIGGSTO_ln_2_riemann
      !USE VBFNLO_HIGGSTO_full_ln
      !USE VBFNLO_HIGGSTO_sp_fun
      !USE VBFNLO_HIGGSTO_units
!
      IMPLICIT NONE
!
      INTEGER it,type,unit
      REAL*8 scal,ps0i,xm0i,scals,ps0,xm0,sgn
      REAL*8, dimension(2) :: value,psi,ps,xmsi,xms,betasc,betas, &
     &        betac,beta,argc,arg,lq,lm,cx,comx,cxs,x,omx,xs,clx,clomx, &
     &        clxs,li2cx,li3cx,li2cxs,li3cxs,copx,clopx,lqs,qcd,lms, &
     &        opx,tau,taus,clxx,clxxs
      REAL*8, dimension(6,2) :: aux,auxs
!
      scals= scal*scal
      ps= psi/scals
      ps0= ps0i/scals
      xms= xmsi/scals
      xm0= xm0i/scal
!
      IF(psi(2).eq.0.d0.and.xmsi(2).eq.0.d0 &
     &   .and.psi(1).le.4.d0*xmsi(1)) THEN
       unit= 1
      ELSE
       unit= 0
      ENDIF
!     
      IF(abs(ps(2)/ps(1)).lt.1.d-10.and.xms(2).eq.0.d0) THEN
       betasc(1)= 1.d0-4.d0*xms(1)/ps(1)
       betasc(2)= 4.d0/(ps(1)*ps(1))*xms(1)*ps(2)
      ELSE
       betasc= co-4.d0*(xms.cq.ps)
      ENDIF
      IF(betasc(2).eq.0.d0) THEN
       betasc(2)= -eps
       betac= (betasc(1)).cr.(betasc(2))
      ELSE
       betac= (betasc(1)).crz.(betasc(2))
      ENDIF
      argc= (betac+co).cq.(betac-co)
!
      betas(1)= 1.d0-4.d0*xm0*xm0/ps0
      betas(2)= -eps
      beta= (betas(1)).cr.(betas(2))
      arg= (beta+co).cq.(beta-co)
!
      IF(arg(2).eq.0.d0) THEN
       x(1)= 1.d0/arg(1) 
       x(2)= -eps      
       sgn= sign(one,x(1))
       xs(1)= x(1)*x(1)
       xs(2)= -sgn*eps
      ELSE
       x= (beta-co).cq.(beta+co)
       xs= x.cp.x  
      ENDIF 
      omx= co-x
      opx= co+x
!
      IF(arg(2).eq.0.d0) arg(2)= eps
      IF(argc(2).eq.0.d0) THEN
       it= 0
       argc(2)= eps
       lq= argc(1).fln.argc(2)
      ELSE
       it= 1
       lq= argc.lnsrs.arg
      ENDIF
      lqs= lq.cp.lq
!
      IF(it.eq.0.d0) THEN
       cx(1)= 1.d0/argc(1) 
       cx(2)= -eps      
       comx= co-cx
       copx= co+cx
       sgn= sign(one,cx(1))
       cxs(1)= cx(1)*cx(1)
       cxs(2)= -sgn*eps
       clx= cx(1).fln.cx(2)
       clxs= clx.cp.clx
       clxx= cxs(1).fln.cxs(2)
       clxxs= clxx.cp.clxx
       clomx= comx(1).fln.comx(2)
       clopx= copx(1).fln.copx(2)
       aux= VBFNLO_HIGGSTO_s_niels_up4(cx) 
       li2cx(1:2)= aux(1,1:2) 
       li3cx(1:2)= aux(2,1:2) 
       auxs= VBFNLO_HIGGSTO_s_niels_up4(cxs) 
       li2cxs(1:2)= auxs(1,1:2) 
       li3cxs(1:2)= auxs(2,1:2) 
      ELSEIF(it.eq.1.d0) THEN
       cx= (betac-co).cq.(betac+co)
       comx= co-cx
       copx= co+cx
       cxs= cx.cp.cx
       clx= cx.lnsrs.x
       clxs= clx.cp.clx
       clxx= cxs.lnsrs.xs
       clxxs= clxx.cp.clxx
       clomx= comx.lnsrs.omx
       clopx= copx.lnsrs.opx
       li2cx= VBFNLO_HIGGSTO_li2_srsz(cx,x,unit)
       li3cx= VBFNLO_HIGGSTO_li3_srsz(cx,x,unit)
       li2cxs= VBFNLO_HIGGSTO_li2_srsz(cxs,xs,unit)
       li3cxs= VBFNLO_HIGGSTO_li3_srsz(cxs,xs,unit)
      ENDIF
!
      IF(xms(2).eq.0.d0) THEN
       lm(1)= log(xms(1))
       lm(2)= 0.d0
      ELSE  
       lm= xms(1).fln.xms(2)
      ENDIF
      lms= lm.cp.lm
!
      tau= xms.cq.ps
      taus= tau.cp.tau
!
      IF(type.eq.0) THEN
!
       qcd= -3.d0/4.d0*(co-12.d0*tau)*rz2 &
     & +1.d0/16.d0*(3.d0*co+344.d0*tau) &
     & -3.d0/2.d0*((co-6.d0*tau).cp.lms) &
     & +3.d0/4.d0*((3.d0*co-14.d0*tau).cp.(betac.cp.clx)) &
     & -3.d0*((li3cxs-2.d0*li3cx+4.d0/3.d0*(li2cx.cp.clx) &
     & +1.d0/3.d0*(clxs.cp.clomx)+rz3*co &
     & -2.d0/3.d0*(clxx.cp.li2cxs)-1.d0/6.d0*(clxxs.cp.clomx) &
     & -1.d0/6.d0* &
     &  (clxxs.cp.clopx)).cp.((co-4.d0*tau).cp.(co-2.d0*tau))) &
     & -1.d0/2.d0*((4.d0*li2cxs-4.d0*li2cx-4.d0*(clomx.cp.clx) &
     & -2.d0*(cx.cp.(clxs.cq.comx))+4.d0*(clxx.cp.clomx) &
     & +4.d0*(clxx.cp.clopx)+(clxxs.cp.(cxs.cq.comx)) &
     & +(clxxs.cp.(cxs.cq.copx))).cp.(betac.cp.(co-4.d0*tau))) &
     & +1.d0/4.d0*(lm.cp.(11.d0*co-108.d0*tau)) &
     & +1.d0/4.d0*(clxs.cp.(3.d0*co+58.d0*taus-28.d0*tau))
!  
      ELSEIF(type.eq.1) THEN
!
      qcd= -3.d0/4.d0*(co-12.d0*tau)*rz2 &
     & +1.d0/16.d0*(67.d0*co-40.d0*tau) &
     & +3.d0/4.d0*((3.d0*co-14.d0*tau).cp.(betac.cp.clx)) &
     & -3.d0*((li3cxs-2.d0*li3cx+4.d0/3.d0*(li2cx.cp.clx) &
     & +1.d0/3.d0*(clxs.cp.clomx)+rz3*co &
     & -2.d0/3.d0*(clxx.cp.li2cxs)-1.d0/6.d0*(clxxs.cp.clomx) &
     & -1.d0/6.d0* &
     &  (clxxs.cp.clopx)).cp.((co-4.d0*tau).cp.(co-2.d0*tau))) &
     & -1.d0/2.d0*((4.d0*li2cxs-4.d0*li2cx-4.d0*(clomx.cp.clx) &
     & -2.d0*(cx.cp.(clxs.cq.comx))+4.d0*(clxx.cp.clomx) &
     & +4.d0*(clxx.cp.clopx)+(clxxs.cp.(cxs.cq.comx)) &
     & +(clxxs.cp.(cxs.cq.copx))).cp.(betac.cp.(co-4.d0*tau))) &
     & -2.d0* &
     &  (((co-8.d0*tau).cp.(betac.cp.lq)).cp.(co-3.d0/4.d0*lm)) &
     & -3.d0*((betac.cp.lq).cp.(lm.cp.tau)) &
     & +4.d0*((co.cp.betac).cp.(lq.cp.tau)) &
     & -9.d0/4.d0*(lm.cp.(co+4.d0*tau)) &
     & +9.d0*(lms.cp.tau) &
     & +1.d0/4.d0*(clxs.cp.(3.d0*co+58.d0*taus-28.d0*tau))
!
      ENDIF
!
      value= qcd
!
      RETURN
!
      END FUNCTION VBFNLO_HIGGSTO_lquarkQCD
!
!------------------------------------------------------------------
!
      FUNCTION VBFNLO_HIGGSTO_lb0af_dm(scal,psi,ps0i,xmsi,xm0i) RESULT(value)
      !USE VBFNLO_HIGGSTO_acmplx_pro
      !USE VBFNLO_HIGGSTO_acmplx_rat
      !USE VBFNLO_HIGGSTO_cmplx_root
      !USE VBFNLO_HIGGSTO_cmplx_rootz
      !USE VBFNLO_HIGGSTO_cmplx_srs_root
      !USE VBFNLO_HIGGSTO_ln_2_riemann
      !USE VBFNLO_HIGGSTO_full_ln
      !USE VBFNLO_HIGGSTO_units
!
      IMPLICIT NONE
!
      REAL*8 scal,ps0i,scals,ps0,xm1,xm2
      REAL*8, dimension(2,2) :: xmsi,xms
      REAL*8, dimension(2) :: value,psi,ps,lambdasc,lambdas, &
     &        lambdac,lambda,argc,arg,llam,xm0,xm0i,aroot, &
     &        root,rat,lnr,xm1c,xm2c
!
      scals= scal*scal
      ps= psi/scals
      ps0= ps0i/scals
      xms= xmsi/scals
      xm0= xm0i/scal
      xm1c(1:2)= xms(1,1:2)
      xm2c(1:2)= xms(2,1:2)
      xm1= xm0(1)*xm0(1)
      xm2= xm0(2)*xm0(2)
      aroot= xm1c.cp.xm2c
      root= (aroot(1).crz.aroot(2))
!
      lambdasc= (ps.cp.ps)+(xm1c.cp.xm1c)+(xm2c.cp.xm2c)-2.d0*( &
     &          (ps.cp.xm1c)+(ps.cp.xm2c)+(xm1c.cp.xm2c))
      lambdas(1)= ps0*ps0+xm1*xm1+xm2*xm2-2.d0*( &
     &            ps0*xm1+ps0*xm2+xm1*xm2)
      lambdas(2)= -eps
      lambdac= (lambdasc(1)).crz.(lambdasc(2))
      lambda= (lambdas(1)).cr.(lambdas(2))
      IF(lambda(2).eq.0.d0) lambda(2)= -eps
!
      argc= 0.5d0*((-ps+xm1c+xm2c-lambdac).cq.root)
!
      arg(1)= 0.5d0*(-ps0+xm1+xm2-lambda(1))/sqrt(xm1*xm2)
      arg(2)= eps
!
      llam= argc.lnsrs.arg
!
      rat= xm1c.cq.xm2c 
      lnr= rat(1).fln.rat(2)
!
      value= 2.d0*co-0.5d0*(((xm1c-xm2c).cq.ps).cp.lnr)- &
     &       ((lambdac.cq.ps).cp.llam)
!
      RETURN
!
      END FUNCTION VBFNLO_HIGGSTO_lb0af_dm
!
!-------------------------------------------------------------------
!
      FUNCTION VBFNLO_HIGGSTO_lb0af_em(scal,psi,ps0i,xmsi,xm0i) RESULT(value)
      !USE VBFNLO_HIGGSTO_acmplx_pro
      !USE VBFNLO_HIGGSTO_acmplx_rat
      !USE VBFNLO_HIGGSTO_cmplx_root
      !USE VBFNLO_HIGGSTO_cmplx_rootz
      !USE VBFNLO_HIGGSTO_cmplx_srs_root
      !USE VBFNLO_HIGGSTO_ln_2_riemann
      !USE VBFNLO_HIGGSTO_full_ln
      !USE VBFNLO_HIGGSTO_units
!
      IMPLICIT NONE
!
      REAL*8 scal,ps0i,xm0i,scals,ps0,xm0
      REAL*8, dimension(2) :: value,psi,ps,xmsi,xms,betasc,betas, &
     &        betac,beta,argc,arg,lbet
!
      scals= scal*scal
      ps= psi/scals
      ps0= ps0i/scals
      xms= xmsi/scals
      xm0= xm0i/scal
!
      betasc= co-4.d0*(xms.cq.ps)
      IF(betasc(2).eq.0.d0) THEN
       betasc(2)= -eps
       betac= (betasc(1)).cr.(betasc(2))
      ELSE
       betac= (betasc(1)).crz.(betasc(2))
      ENDIF 
      argc= (betac+co).cq.(betac-co)
      IF(argc(2).eq.0.d0) THEN
       argc(2)= eps
       lbet= argc(1).fln.argc(2)
      ELSE
       betas(1)= 1.d0-4.d0*(xm0*xm0)/ps0
       betas(2)= -eps
       beta= (betas(1)).cr.(betas(2))
       arg= (beta+co).cq.(beta-co)
       IF(arg(2).eq.0.d0) arg(2)= eps
       lbet= argc.lnsrs.arg
      ENDIF
!
      value= 2.d0*co-(betac.cp.lbet)
!
      RETURN
!
      END FUNCTION VBFNLO_HIGGSTO_lb0af_em
!
!-----------------------------------------------------------------------
!
      FUNCTION VBFNLO_HIGGSTO_lb021_dm_cp(scal,psi,ps0i,xmsi,xm0i) RESULT(value)
      !USE VBFNLO_HIGGSTO_cmplx_root
      !USE VBFNLO_HIGGSTO_cmplx_rootz
      !USE VBFNLO_HIGGSTO_cmplx_srs_root
      !USE VBFNLO_HIGGSTO_ln_2_riemann
      !USE VBFNLO_HIGGSTO_acmplx_pro
      !USE VBFNLO_HIGGSTO_acmplx_rat
      !USE VBFNLO_HIGGSTO_full_ln
      !USE VBFNLO_HIGGSTO_units
!
      IMPLICIT NONE
!

      REAL*8 scal,ps0i,ps0,xm1,xm2
      REAL*8, intent(in), dimension(2) :: xm0i
      REAL*8, intent(in), dimension(2,2) :: xmsi
      REAL*8, dimension(2) :: value,xm1c,xm2c
      REAL*8, dimension(2,2) :: xms
      REAL*8, dimension(2) :: xm0
      REAL*8, dimension(2) :: psi,ps,aroot,root,lambdasc,lambdas, &
     &        lambdac,lambda,argc,arg,llam,l1,l2 
!
      ps= psi
      ps0= ps0i
      xms= xmsi
      xm0= xm0i
!
      xm1c(1:2)= xms(1,1:2)
      xm2c(1:2)= xms(2,1:2)
      xm1c= xm1c.cq.ps
      xm2c= xm2c.cq.ps
      xm1= xm0(1)*xm0(1)/ps0
      xm2= xm0(2)*xm0(2)/ps0
!
      aroot= xm1c.cp.xm2c
      root= (aroot(1).crz.aroot(2))
!
      lambdasc= co+(xm1c.cp.xm1c)+(xm2c.cp.xm2c)-2.d0*( &
     &          xm1c+xm2c+(xm1c.cp.xm2c))
      lambdas(1)= 1.d0+xm1*xm1+xm2*xm2-2.d0*( &
     &            xm1+xm2+xm1*xm2)
      lambdas(2)= -eps
      lambdac= (lambdasc(1)).crz.(lambdasc(2))
      lambda= (lambdas(1)).cr.(lambdas(2))
      IF(lambda(2).eq.0.d0) lambda(2)= -eps
!
      argc= 0.5d0*((-co+xm1c+xm2c-lambdac).cq.root)
!
      arg(1)= 0.5d0*(-1.d0+xm1+xm2-lambda(1))/sqrt(xm1*xm2)
      arg(2)= eps
!
      llam= argc.lnsrs.arg
!
      l1= xm1c(1).fln.xm1c(2)
      l2= xm2c(1).fln.xm2c(2)
!
      value= (((xm1c-xm2c-co).cq.lambdac).cp.llam)+0.5d0*(l1-l2) 
!
      RETURN
!
      END FUNCTION VBFNLO_HIGGSTO_lb021_dm_cp
!
!----------------------------------------------------------------------
!

      end module
