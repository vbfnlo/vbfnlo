      subroutine qqwwqq_mg2(p,pbar,v,sign, nlo, 
     1                  uucc,uuss,ddcc,ddss,udsc,ducs)

      implicit none                

#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"

#include "VBFNLO/utilities/process.inc"

      double precision pbar(0:3,4+max_v)
      double precision uucc,uuss,ddcc,ddss,udsc,ducs
      double precision uucc2,uuss2,ddcc2,ddss2,udsc2,ducs2
      double precision p(0:3,max_p),v(0:3,max_v)
      integer  sign(4+max_v), nlo, i, j,mu,nu
      real*8 q12(0:4,3), q34(0:4,3)
      complex*16 mg2_tens_1d(0:5)
      complex*16 mg2_tens_2d(0:3,0:3)
      complex*16 mg2_tens_2d_comb(0:3,0:3,2)
      double precision sumtest

      call Calc_Momentum_Transfer(p, v, q12,q34,1)

      j = 1
         
#ifdef WITH_MADGRAPH

         SELECT CASE(procid)
         CASE(WPhadWMjj, WPWMhadjj)            ! semileptonic decay modes

            call anomal_formfactor(q12(0,j),q34(0,j),qp(0),qm(0))
            call atoww_mg2_had(v,aww)
            call ztoww_mg2_had(v,zww)
            call aatoww_mg2_had(q12(0,j),q34(0,j),v,aaww(0,0,j))
            call aztoww_mg2_had(q12(0,j),q34(0,j),v,azww(0,0,j))
            call aztoww_mg2_had(q34(0,j),q12(0,j),v,zaww(0,0,j))
            call zztoww_mg2_had(q12(0,j),q34(0,j),v,zzww(0,0,j))
            call wwtoww_mg2_had(q12(0,j),q34(0,j),v,wwww6(0,0,j)) ! q12 = W-
            call wwtoww_mg2_had(q34(0,j),q12(0,j),v,wwww5(0,0,j)) ! q12 = W+

            call WVtoWP_mg2_had(2,q34(0,j),v,NCwpa(0,0,1,j),NCwpz(0,0,1,j)) !emit W- on upper
            call WVtoWP_mg2_had(2,q12(0,j),v,NCwpa(0,0,2,j),NCwpz(0,0,2,j)) !emit W- on lower
            call WVtoWP_mg2_had(1,q34(0,j),v,CCwpa(0,0,1,j),CCwpz(0,0,1,j)) !emit W- on upper
            call WVtoWP_mg2_had(1,q12(0,j),v,CCwpa(0,0,2,j),CCwpz(0,0,2,j)) !emit W- on lower
C     for WV --> mu- nu_mu
            call WVtoWM_mg2_had(2,q34(0,j),v,NCwma(0,0,1,j),NCwmz(0,0,1,j)) !emit W+ on upper
            call WVtoWM_mg2_had(2,q12(0,j),v,NCwma(0,0,2,j),NCwmz(0,0,2,j)) !emit W+ on lower
            call WVtoWM_mg2_had(1,q34(0,j),v,CCwma(0,0,1,j),CCwmz(0,0,1,j)) !emit W+ on upper
            call WVtoWM_mg2_had(1,q12(0,j),v,CCwma(0,0,2,j),CCwmz(0,0,2,j)) !emit W+ on lower

         call qqwwqq(pbar,sign, nlo, 1,
     1                    uucc2,uuss2,ddcc2,ddss2,udsc2,ducs2)

            call anomal_formfactor(q12(0,j),q34(0,j),qp(0),qm(0))
            call atoww_had(v,aww)
            call ztoww_had(v,zww)
            call aatoww_had(q12(0,j),q34(0,j),v,aaww(0,0,j))
            call aztoww_had(q12(0,j),q34(0,j),v,azww(0,0,j))
            call aztoww_had(q34(0,j),q12(0,j),v,zaww(0,0,j))
            call zztoww_had(q12(0,j),q34(0,j),v,zzww(0,0,j))
            call wwtoww_had(q12(0,j),q34(0,j),v,wwww6(0,0,j)) ! q12 = W-
            call wwtoww_had(q34(0,j),q12(0,j),v,wwww5(0,0,j)) ! q12 = W+

            call WVtoWP_had(2,q34(0,j),v,NCwpa(0,0,1,j),NCwpz(0,0,1,j)) !emit W- on upper
            call WVtoWP_had(2,q12(0,j),v,NCwpa(0,0,2,j),NCwpz(0,0,2,j)) !emit W- on lower
            call WVtoWP_had(1,q34(0,j),v,CCwpa(0,0,1,j),CCwpz(0,0,1,j)) !emit W- on upper
            call WVtoWP_had(1,q12(0,j),v,CCwpa(0,0,2,j),CCwpz(0,0,2,j)) !emit W- on lower
C     for WV --> mu- nu_mu
            call WVtoWM_had(2,q34(0,j),v,NCwma(0,0,1,j),NCwmz(0,0,1,j)) !emit W+ on upper
            call WVtoWM_had(2,q12(0,j),v,NCwma(0,0,2,j),NCwmz(0,0,2,j)) !emit W+ on lower
            call WVtoWM_had(1,q34(0,j),v,CCwma(0,0,1,j),CCwmz(0,0,1,j)) !emit W+ on upper
            call WVtoWM_had(1,q12(0,j),v,CCwma(0,0,2,j),CCwmz(0,0,2,j)) !emit W+ on lower


         do mu = 0,3
          do nu = 0,3
           mg2_tens_2d(mu,nu) = dcmplx(0.,0.)
           mg2_tens_2d_comb(mu,nu,1) = dcmplx(0.,0.)
           mg2_tens_2d_comb(mu,nu,2) = dcmplx(0.,0.)
          enddo
         enddo
         do i = 0,5
           mg2_tens_1d(i) = dcmplx(0.,0.)
         enddo

         print*, ""
         call atoww_mg2_had(v,mg2_tens_1d)
         sumtest = 0
         do i = 0,5
           sumtest = sumtest + abs( aww(i) / mg2_tens_1d(i) - 1)
!           print*, i, aww(i), mg2_tens_1d(i), abs( aww(i) / mg2_tens_1d(i) - 1)
           mg2_tens_1d(i) = dcmplx(0.,0.)
         enddo
         print*, "tensor a to ww:", sumtest

!         print*, ""
         call ztoww_mg2_had(v,mg2_tens_1d)
         sumtest = 0
         do i = 0,5
           sumtest = sumtest + abs( zww(i) / mg2_tens_1d(i) - 1)
!           print*, i, zww(i), mg2_tens_1d(i), abs( zww(i) / mg2_tens_1d(i) - 1)
           mg2_tens_1d(i) = dcmplx(0.,0.)
         enddo
         print*, "tensor z to ww:", sumtest

!         print*, ""
         call aatoww_mg2_had(q12(0,j),q34(0,j),v,mg2_tens_2d(0,0))
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( aaww(mu,nu,1) / mg2_tens_2d(mu,nu) - 1)
!           print*, mu,nu, aaww(mu,nu,1), mg2_tens_2d(mu,nu), abs( aaww(mu,nu,1) / mg2_tens_2d(mu,nu) - 1)
           mg2_tens_2d(mu,nu) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor aa to ww:", sumtest

!         print*, ""
         call aztoww_mg2_had(q12(0,j),q34(0,j),v,mg2_tens_2d(0,0))
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( azww(mu,nu,1) / mg2_tens_2d(mu,nu) - 1)
!           print*, mu,nu, azww(mu,nu,1), mg2_tens_2d(mu,nu), abs( azww(mu,nu,1) / mg2_tens_2d(mu,nu) - 1)
           mg2_tens_2d(mu,nu) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor az to ww:", sumtest

!         print*, ""
         call aztoww_mg2_had(q34(0,j),q12(0,j),v,mg2_tens_2d(0,0))
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( zaww(mu,nu,1) / mg2_tens_2d(mu,nu) - 1)
!           print*, mu,nu, zaww(mu,nu,1), mg2_tens_2d(mu,nu), abs( zaww(mu,nu,1) / mg2_tens_2d(mu,nu) - 1)
           mg2_tens_2d(mu,nu) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor za to ww:", sumtest

!         print*, ""
         call zztoww_mg2_had(q12(0,j),q34(0,j),v,mg2_tens_2d(0,0))
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( zzww(mu,nu,1) / mg2_tens_2d(mu,nu) - 1)
!           print*, mu,nu, zzww(mu,nu,1), mg2_tens_2d(mu,nu), abs( zzww(mu,nu,1) / mg2_tens_2d(mu,nu) - 1)
           mg2_tens_2d(mu,nu) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor zz to ww:", sumtest

!         print*, ""
         call wwtoww_mg2_had(q12(0,j),q34(0,j),v,mg2_tens_2d(0,0)) ! q12 = W-
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( wwww6(mu,nu,1) / mg2_tens_2d(mu,nu) - 1)
!           print*, mu,nu, wwww6(mu,nu,1), mg2_tens_2d(mu,nu), abs( wwww6(mu,nu,1) / mg2_tens_2d(mu,nu) - 1)
           mg2_tens_2d(mu,nu) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor ww to ww 6:", sumtest

!         print*, ""
         call wwtoww_mg2_had(q34(0,j),q12(0,j),v,mg2_tens_2d(0,0)) ! q12 = W+
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( wwww5(mu,nu,1) / mg2_tens_2d(mu,nu) - 1)
!           print*, mu,nu, wwww5(mu,nu,1), mg2_tens_2d(mu,nu), abs( wwww5(mu,nu,1) / mg2_tens_2d(mu,nu) - 1)
           mg2_tens_2d(mu,nu) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor ww to ww 5:", sumtest

!         print*, ""
         call WVtoWP_mg2_had(2,q34(0,j),v,mg2_tens_2d_comb(0,0,1),mg2_tens_2d_comb(0,0,2)) !emit W- on upper
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( NCwpa(mu,nu,1,1) / mg2_tens_2d_comb(mu,nu,1) - 1)
!           print*, mu,nu, NCwpa(mu,nu,1,1), mg2_tens_2d_comb(mu,nu,1), abs( NCwpa(mu,nu,1,1) / mg2_tens_2d_comb(mu,nu,1) - 1)
           mg2_tens_2d_comb(mu,nu,1) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor wa to w+ 1:", sumtest
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( NCwpz(mu,nu,1,1) / mg2_tens_2d_comb(mu,nu,2) - 1)
!           print*, mu,nu, NCwpz(mu,nu,1,1), mg2_tens_2d_comb(mu,nu,2), abs( NCwpz(mu,nu,1,1) / mg2_tens_2d_comb(mu,nu,2) - 1)
           mg2_tens_2d_comb(mu,nu,2) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor wz to w+ 1:", sumtest

!         print*, ""
         call WVtoWP_mg2_had(2,q12(0,j),v,mg2_tens_2d_comb(0,0,1),mg2_tens_2d_comb(0,0,2)) !emit W- on upper
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( NCwpa(mu,nu,2,1) / mg2_tens_2d_comb(mu,nu,1) - 1)
!           print*, mu,nu, NCwpa(mu,nu,2,1), mg2_tens_2d_comb(mu,nu,1), abs( NCwpa(mu,nu,2,1) / mg2_tens_2d_comb(mu,nu,1) - 1)
           mg2_tens_2d_comb(mu,nu,1) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor wa to w+ 2:", sumtest
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( NCwpz(mu,nu,2,1) / mg2_tens_2d_comb(mu,nu,2) - 1)
!           print*, mu,nu, NCwpz(mu,nu,2,1), mg2_tens_2d_comb(mu,nu,2), abs( NCwpz(mu,nu,2,1) / mg2_tens_2d_comb(mu,nu,2) - 1)
           mg2_tens_2d_comb(mu,nu,2) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor wz to w+ 2:", sumtest

!         print*, ""
         call WVtoWP_mg2_had(1,q34(0,j),v,mg2_tens_2d_comb(0,0,1),mg2_tens_2d_comb(0,0,2)) !emit W- on upper
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( CCwpa(mu,nu,1,1) / mg2_tens_2d_comb(mu,nu,1) - 1)
!           print*, mu,nu, CCwpa(mu,nu,1,1), mg2_tens_2d_comb(mu,nu,1), abs( CCwpa(mu,nu,1,1) / mg2_tens_2d_comb(mu,nu,1) - 1)
           mg2_tens_2d_comb(mu,nu,1) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor wa to w+ 3:", sumtest
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( CCwpz(mu,nu,1,1) / mg2_tens_2d_comb(mu,nu,2) - 1)
!           print*, mu,nu, CCwpz(mu,nu,1,1), mg2_tens_2d_comb(mu,nu,2), abs( CCwpz(mu,nu,1,1) / mg2_tens_2d_comb(mu,nu,2) - 1)
           mg2_tens_2d_comb(mu,nu,2) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor wz to w+ 3:", sumtest

!         print*, ""
         call WVtoWP_mg2_had(1,q12(0,j),v,mg2_tens_2d_comb(0,0,1),mg2_tens_2d_comb(0,0,2)) !emit W- on upper
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( CCwpa(mu,nu,2,1) / mg2_tens_2d_comb(mu,nu,1) - 1)
!           print*, mu,nu, CCwpa(mu,nu,2,1), mg2_tens_2d_comb(mu,nu,1), abs( CCwpa(mu,nu,2,1) / mg2_tens_2d_comb(mu,nu,1) - 1)
           mg2_tens_2d_comb(mu,nu,1) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor wa to w+ 4:", sumtest
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( CCwpz(mu,nu,2,1) / mg2_tens_2d_comb(mu,nu,2) - 1)
!           print*, mu,nu, CCwpz(mu,nu,2,1), mg2_tens_2d_comb(mu,nu,2), abs( CCwpz(mu,nu,2,1) / mg2_tens_2d_comb(mu,nu,2) - 1)
           mg2_tens_2d_comb(mu,nu,2) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor wz to w+ 4:", sumtest

!         print*, ""
         call WVtoWM_mg2_had(2,q34(0,j),v,mg2_tens_2d_comb(0,0,1),mg2_tens_2d_comb(0,0,2)) !emit W- on upper
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( NCwma(mu,nu,1,1) / mg2_tens_2d_comb(mu,nu,1) - 1)
!           print*, mu,nu, NCwma(mu,nu,1,1), mg2_tens_2d_comb(mu,nu,1), abs( NCwma(mu,nu,1,1) / mg2_tens_2d_comb(mu,nu,1) - 1)
           mg2_tens_2d_comb(mu,nu,1) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor wa to w- 1:", sumtest
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( NCwmz(mu,nu,1,1) / mg2_tens_2d_comb(mu,nu,2) - 1)
!           print*, mu,nu, NCwmz(mu,nu,1,1), mg2_tens_2d_comb(mu,nu,2), abs( NCwmz(mu,nu,1,1) / mg2_tens_2d_comb(mu,nu,2) - 1)
           mg2_tens_2d_comb(mu,nu,2) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor wz to w- 1:", sumtest

!         print*, ""
         call WVtoWM_mg2_had(2,q12(0,j),v,mg2_tens_2d_comb(0,0,1),mg2_tens_2d_comb(0,0,2)) !emit W- on upper
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( NCwma(mu,nu,2,1) / mg2_tens_2d_comb(mu,nu,1) - 1)
!           print*, mu,nu, NCwma(mu,nu,2,1), mg2_tens_2d_comb(mu,nu,1), abs( NCwma(mu,nu,2,1) / mg2_tens_2d_comb(mu,nu,1) - 1)
           mg2_tens_2d_comb(mu,nu,1) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor wa to w- 2:", sumtest
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( NCwmz(mu,nu,2,1) / mg2_tens_2d_comb(mu,nu,2) - 1)
!           print*, mu,nu, NCwmz(mu,nu,2,1), mg2_tens_2d_comb(mu,nu,2), abs( NCwmz(mu,nu,2,1) / mg2_tens_2d_comb(mu,nu,2) - 1)
           mg2_tens_2d_comb(mu,nu,2) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor wz to w- 2:", sumtest

!         print*, ""
         call WVtoWM_mg2_had(1,q34(0,j),v,mg2_tens_2d_comb(0,0,1),mg2_tens_2d_comb(0,0,2)) !emit W- on upper
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( CCwma(mu,nu,1,1) / mg2_tens_2d_comb(mu,nu,1) - 1)
!           print*, mu,nu, CCwma(mu,nu,1,1), mg2_tens_2d_comb(mu,nu,1), abs( CCwma(mu,nu,1,1) / mg2_tens_2d_comb(mu,nu,1) - 1)
           mg2_tens_2d_comb(mu,nu,1) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor wa to w- 3:", sumtest
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( CCwmz(mu,nu,1,1) / mg2_tens_2d_comb(mu,nu,2) - 1)
!           print*, mu,nu, CCwmz(mu,nu,1,1), mg2_tens_2d_comb(mu,nu,2), abs( CCwmz(mu,nu,1,1) / mg2_tens_2d_comb(mu,nu,2) - 1)
           mg2_tens_2d_comb(mu,nu,2) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor wz to w- 3:", sumtest

!         print*, ""
         call WVtoWM_mg2_had(1,q12(0,j),v,mg2_tens_2d_comb(0,0,1),mg2_tens_2d_comb(0,0,2)) !emit W- on upper
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( CCwma(mu,nu,2,1) / mg2_tens_2d_comb(mu,nu,1) - 1)
!           print*, mu,nu, CCwma(mu,nu,2,1), mg2_tens_2d_comb(mu,nu,1), abs( CCwma(mu,nu,2,1) / mg2_tens_2d_comb(mu,nu,1) - 1)
           mg2_tens_2d_comb(mu,nu,1) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor wa to w- 4:", sumtest
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( CCwmz(mu,nu,2,1) / mg2_tens_2d_comb(mu,nu,2) - 1)
!           print*, mu,nu, CCwmz(mu,nu,2,1), mg2_tens_2d_comb(mu,nu,2), abs( CCwmz(mu,nu,2,1) / mg2_tens_2d_comb(mu,nu,2) - 1)
           mg2_tens_2d_comb(mu,nu,2) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor wz to w- 4:", sumtest

         print*, ""
         print*,"1 - uucc/uucc2 =   ",abs(1-uucc/uucc2)
         print*,"1 - uuss/uuss2 =   ",abs(1-uuss/uuss2)
         print*,"1 - ddcc/ddcc2 =   ",abs(1-ddcc/ddcc2)
         print*,"1 - ddss/ddss2 =   ",abs(1-ddss/ddss2)
         print*,"1 - udsc/udsc2 =   ",abs(1-udsc/udsc2)
         print*,"1 - ducs/ducs2 =   ",abs(1-ducs/ducs2)
         read(*,*)

         end select

#endif


      end


      subroutine AAtoWW_mg2_had(Q1,Q2,P,aatens)
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : a a -> e+ ve mu- vm~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=  21) 
      INTEGER    NWAVEFUNCS
      PARAMETER (NWAVEFUNCS=  38) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 Q1(0:3), Q2(0:3), P(0:3,3:6)
      complex*16 aatens(0:3,0:3)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J, MU, NU
      COMPLEX*16 AMP(NGRAPHS)
      COMPLEX*16 W(6,NWAVEFUNCS),WX1(6),WX2(6)
      INTEGER NHEL(3:6), ihel
      DATA (NHEL(IHEL),IHEL=3,6) / +1,-1,-1,+1/
      Save NHEL
c  lepton spinors and W+- polarization vectors
#include "tensor.inc"

      logical lswitch
      parameter (lswitch=.true.)
C  
C GLOBAL VARIABLES
C  
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/process.inc"
C                                                                          
C ----------
C BEGIN CODE
C ----------
      if (lswitch) then      ! read in lepton an W states from common block
         do i=1,6
            w(i,3) = wep(i)
            w(i,4) = wve(i)
            w(i,5) = wmu(i)
            w(i,6) = wvm(i)
         enddo
      endif

      DO MU = 0,3
         DO NU = 0,3

      CALL VCARTX(Q1(0),ZERO ,ZERO, MU,-1,W(1,1   ))        
      CALL VCARTX(Q2(0),ZERO ,ZERO, NU,-1,W(1,2   )) 

      if (procid.eq.WPhadWMjj) then

      CALL FVIXXX(W(1,3   ),W(1,2   ),GAD ,ZERO    ,ZERO    ,W(1,7   ))    
      CALL JIOXXX(W(1,7   ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,8   ))    
      CALL JVVXXX(W(1,8   ),W(1,1   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     9   ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,9   ),GWF ,AMP(1   ))            
      CALL FVOXXX(W(1,5   ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,10  ))    
      CALL IOVXXX(W(1,6   ),W(1,10  ),W(1,8   ),GWF ,AMP(2   ))            
      CALL FVOXXX(W(1,4   ),W(1,1   ),GAU ,ZERO    ,ZERO    ,W(1,11  ))    
      CALL JIOXXX(W(1,3   ),W(1,11  ),GWF ,WMASS   ,WWIDTH  ,W(1,12  ))    
      CALL JVVXXX(W(1,12  ),W(1,2   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     13  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,13  ),GWF ,AMP(3   ))            
      CALL FVOXXX(W(1,5   ),W(1,2   ),GAL ,ZERO    ,ZERO    ,W(1,14  ))    
      CALL IOVXXX(W(1,6   ),W(1,14  ),W(1,12  ),GWF ,AMP(4   ))            
      CALL FVOXXX(W(1,11  ),W(1,2   ),GAU ,ZERO    ,ZERO    ,W(1,15  ))    
      CALL JIOXXX(W(1,3   ),W(1,15  ),GWF ,WMASS   ,WWIDTH  ,W(1,16  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,16  ),GWF ,AMP(5   ))            
      CALL FVIXXX(W(1,7   ),W(1,1   ),GAD ,ZERO    ,ZERO    ,W(1,17  ))    
      CALL JIOXXX(W(1,17  ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,18  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,18  ),GWF ,AMP(6   ))            
      CALL JIOXXX(W(1,7   ),W(1,11  ),GWF ,WMASS   ,WWIDTH  ,W(1,19  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,19  ),GWF ,AMP(7   ))            
      CALL FVOXXX(W(1,4   ),W(1,2   ),GAU ,ZERO    ,ZERO    ,W(1,20  ))    
      CALL JIOXXX(W(1,3   ),W(1,20  ),GWF ,WMASS   ,WWIDTH  ,W(1,21  ))    
      CALL JVVXXX(W(1,21  ),W(1,1   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     22  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,22  ),GWF ,AMP(8   ))            
      CALL IOVXXX(W(1,6   ),W(1,10  ),W(1,21  ),GWF ,AMP(9   ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GAD ,ZERO    ,ZERO    ,W(1,23  ))    
      CALL JIOXXX(W(1,23  ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,24  ))    
      CALL JVVXXX(W(1,24  ),W(1,2   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     25  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,25  ),GWF ,AMP(10  ))            
      CALL IOVXXX(W(1,6   ),W(1,14  ),W(1,24  ),GWF ,AMP(11  ))            
      CALL FVOXXX(W(1,20  ),W(1,1   ),GAU ,ZERO    ,ZERO    ,W(1,26  ))    
      CALL JIOXXX(W(1,3   ),W(1,26  ),GWF ,WMASS   ,WWIDTH  ,W(1,27  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,27  ),GWF ,AMP(12  ))            
      CALL FVIXXX(W(1,23  ),W(1,2   ),GAD ,ZERO    ,ZERO    ,W(1,28  ))    
      CALL JIOXXX(W(1,28  ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,29  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,29  ),GWF ,AMP(13  ))            
      CALL JIOXXX(W(1,23  ),W(1,20  ),GWF ,WMASS   ,WWIDTH  ,W(1,30  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,30  ),GWF ,AMP(14  ))            
      CALL JIOXXX(W(1,3   ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,31  ))    
      CALL JVVXXX(W(1,31  ),W(1,2   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     32  ))                                                          
      CALL JVVXXX(W(1,32  ),W(1,1   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     33  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,33  ),GWF ,AMP(15  ))            
      CALL FVOXXX(W(1,10  ),W(1,2   ),GAL ,ZERO    ,ZERO    ,W(1,34  ))    
      CALL IOVXXX(W(1,6   ),W(1,34  ),W(1,31  ),GWF ,AMP(16  ))            
      CALL IOVXXX(W(1,6   ),W(1,10  ),W(1,32  ),GWF ,AMP(17  ))            
      CALL JVVXXX(W(1,31  ),W(1,1   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     35  ))                                                          
      CALL JVVXXX(W(1,35  ),W(1,2   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     36  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,36  ),GWF ,AMP(18  ))            
      CALL FVOXXX(W(1,14  ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,37  ))    
      CALL IOVXXX(W(1,6   ),W(1,37  ),W(1,31  ),GWF ,AMP(19  ))            
      CALL IOVXXX(W(1,6   ),W(1,14  ),W(1,35  ),GWF ,AMP(20  ))            
      CALL JW3WXX(W(1,1   ),W(1,31  ),W(1,2   ),GWWA ,GWWA ,WMASS   ,      
     &     WWIDTH  ,W(1,38  ))                                             
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,38  ),GWF ,AMP(21  ))            
      aatens(MU,NU) = -AMP(   1)-AMP(   2)-AMP(   3)-AMP(   4)-AMP(   5)
     &             -AMP(   6)-AMP(   7)-AMP(   8)-AMP(   9)-AMP(  10)
     &             -AMP(  11)-AMP(  12)-AMP(  13)-AMP(  14)-AMP(  15)
     &             -AMP(  16)-AMP(  17)-AMP(  18)-AMP(  19)-AMP(  20)
     &             -AMP(  21)

      elseif (procid.eq.WPWMhadjj) then

      CALL FVIXXX(W(1,3   ),W(1,2   ),GAL ,ZERO    ,ZERO    ,W(1,7   ))    
      CALL JIOXXX(W(1,7   ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,8   ))    
      CALL FVOXXX(W(1,5   ),W(1,8   ),GWF ,ZERO    ,ZERO    ,W(1,9   ))    
      CALL IOVXXX(W(1,6   ),W(1,9   ),W(1,1   ),GAU ,AMP(1   ))            
      CALL JVVXXX(W(1,8   ),W(1,1   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     10  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,10  ),GWF ,AMP(2   ))            
      CALL FVOXXX(W(1,5   ),W(1,1   ),GAD ,ZERO    ,ZERO    ,W(1,11  ))    
      CALL IOVXXX(W(1,6   ),W(1,11  ),W(1,8   ),GWF ,AMP(3   ))            
      CALL FVIXXX(W(1,7   ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,12  ))    
      CALL JIOXXX(W(1,12  ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,13  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,13  ),GWF ,AMP(4   ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,14  ))    
      CALL JIOXXX(W(1,14  ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,15  ))    
      CALL FVOXXX(W(1,5   ),W(1,15  ),GWF ,ZERO    ,ZERO    ,W(1,16  ))    
      CALL IOVXXX(W(1,6   ),W(1,16  ),W(1,2   ),GAU ,AMP(5   ))            
      CALL JVVXXX(W(1,15  ),W(1,2   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     17  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,17  ),GWF ,AMP(6   ))            
      CALL FVOXXX(W(1,5   ),W(1,2   ),GAD ,ZERO    ,ZERO    ,W(1,18  ))    
      CALL IOVXXX(W(1,6   ),W(1,18  ),W(1,15  ),GWF ,AMP(7   ))            
      CALL FVIXXX(W(1,14  ),W(1,2   ),GAL ,ZERO    ,ZERO    ,W(1,19  ))    
      CALL JIOXXX(W(1,19  ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,20  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,20  ),GWF ,AMP(8   ))            
      CALL JIOXXX(W(1,3   ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,21  ))    
      CALL JVVXXX(W(1,21  ),W(1,2   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     22  ))                                                          
      CALL FVOXXX(W(1,5   ),W(1,22  ),GWF ,ZERO    ,ZERO    ,W(1,23  ))    
      CALL IOVXXX(W(1,6   ),W(1,23  ),W(1,1   ),GAU ,AMP(9   ))            
      CALL FVOXXX(W(1,11  ),W(1,21  ),GWF ,ZERO    ,ZERO    ,W(1,24  ))    
      CALL IOVXXX(W(1,6   ),W(1,24  ),W(1,2   ),GAU ,AMP(10  ))            
      CALL JVVXXX(W(1,22  ),W(1,1   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     25  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,25  ),GWF ,AMP(11  ))            
      CALL FVOXXX(W(1,11  ),W(1,2   ),GAD ,ZERO    ,ZERO    ,W(1,26  ))    
      CALL IOVXXX(W(1,6   ),W(1,26  ),W(1,21  ),GWF ,AMP(12  ))            
      CALL IOVXXX(W(1,6   ),W(1,11  ),W(1,22  ),GWF ,AMP(13  ))            
      CALL FVOXXX(W(1,18  ),W(1,21  ),GWF ,ZERO    ,ZERO    ,W(1,27  ))    
      CALL IOVXXX(W(1,6   ),W(1,27  ),W(1,1   ),GAU ,AMP(14  ))            
      CALL JVVXXX(W(1,21  ),W(1,1   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     28  ))                                                          
      CALL FVOXXX(W(1,5   ),W(1,28  ),GWF ,ZERO    ,ZERO    ,W(1,29  ))    
      CALL IOVXXX(W(1,6   ),W(1,29  ),W(1,2   ),GAU ,AMP(15  ))            
      CALL JVVXXX(W(1,28  ),W(1,2   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     30  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,30  ),GWF ,AMP(16  ))            
      CALL FVOXXX(W(1,18  ),W(1,1   ),GAD ,ZERO    ,ZERO    ,W(1,31  ))    
      CALL IOVXXX(W(1,6   ),W(1,31  ),W(1,21  ),GWF ,AMP(17  ))            
      CALL IOVXXX(W(1,6   ),W(1,18  ),W(1,28  ),GWF ,AMP(18  ))            
      CALL JW3WXX(W(1,1   ),W(1,21  ),W(1,2   ),GWWA ,GWWA ,WMASS   ,      
     &     WWIDTH  ,W(1,32  ))                                             
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,32  ),GWF ,AMP(19  ))            
      CALL FVOXXX(W(1,5   ),W(1,21  ),GWF ,ZERO    ,ZERO    ,W(1,33  ))    
      CALL FVOXXX(W(1,33  ),W(1,2   ),GAU ,ZERO    ,ZERO    ,W(1,34  ))    
      CALL IOVXXX(W(1,6   ),W(1,34  ),W(1,1   ),GAU ,AMP(20  ))            
      CALL FVOXXX(W(1,33  ),W(1,1   ),GAU ,ZERO    ,ZERO    ,W(1,35  ))    
      CALL IOVXXX(W(1,6   ),W(1,35  ),W(1,2   ),GAU ,AMP(21  ))            
      aatens(MU,NU) = -AMP(   1)-AMP(   2)-AMP(   3)-AMP(   4)-AMP(   5)
     &             -AMP(   6)-AMP(   7)-AMP(   8)-AMP(   9)-AMP(  10)
     &             -AMP(  11)-AMP(  12)-AMP(  13)-AMP(  14)-AMP(  15)
     &             -AMP(  16)-AMP(  17)-AMP(  18)-AMP(  19)-AMP(  20)
     &             -AMP(  21)

      else

       aatens(MU,NU) = dcmplx(0.,0.)

      endif

         enddo
      enddo
      END


      subroutine AZtoWW_mg2_had(Q1,Q2,P,az)
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : a z -> e+ ve mu- vm~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=  25) 
      INTEGER    NWAVEFUNCS
      PARAMETER (NWAVEFUNCS=  42) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 Q1(0:3), Q2(0:3), P(0:3,3:6)
      complex*16 az(0:3,0:3)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J, MU, NU
      COMPLEX*16 AMP(NGRAPHS)
      COMPLEX*16 W(6,NWAVEFUNCS),WX1(6),WX2(6)
      INTEGER NHEL(3:6), ihel
      DATA (NHEL(IHEL),IHEL=3,6) / +1,-1,-1,+1/
      Save NHEL
c  lepton spinors and W+- polarization vectors
#include "tensor.inc"

      logical lswitch
      parameter (lswitch=.true.)
C  
C GLOBAL VARIABLES
C  
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/process.inc"
C                                                                          
C ----------
C BEGIN CODE
C ----------
      if (lswitch) then      ! read in lepton an W states from common block
         do i=1,6
            w(i,3) = wep(i)
            w(i,4) = wve(i)
            w(i,5) = wmu(i)
            w(i,6) = wvm(i)
         enddo
      endif

      DO MU = 0,3
         DO NU = 0,3
      CALL VCARTX(Q1(0),ZERO ,ZERO, MU,-1,W(1,1   ))        
      CALL VCARTX(Q2(0),ZMASS ,Zwidth, NU,-1,W(1,2   )) 

      if (procid.eq.WPhadWMjj) then

      CALL FVIXXX(W(1,3   ),W(1,2   ),GZD ,ZERO    ,ZERO    ,W(1,7   ))    
      CALL JIOXXX(W(1,7   ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,8   ))    
      CALL JVVXXX(W(1,8   ),W(1,1   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     9   ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,9   ),GWF ,AMP(1   ))            
      CALL FVOXXX(W(1,5   ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,10  ))    
      CALL IOVXXX(W(1,6   ),W(1,10  ),W(1,8   ),GWF ,AMP(2   ))            
      CALL FVOXXX(W(1,4   ),W(1,1   ),GAU ,ZERO    ,ZERO    ,W(1,11  ))    
      CALL JIOXXX(W(1,3   ),W(1,11  ),GWF ,WMASS   ,WWIDTH  ,W(1,12  ))    
      CALL FVOXXX(W(1,5   ),W(1,12  ),GWF ,ZERO    ,ZERO    ,W(1,13  ))    
      CALL IOVXXX(W(1,6   ),W(1,13  ),W(1,2   ),GZN ,AMP(3   ))            
      CALL JVVXXX(W(1,12  ),W(1,2   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     14  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,14  ),GWF ,AMP(4   ))            
      CALL FVOXXX(W(1,5   ),W(1,2   ),GZL ,ZERO    ,ZERO    ,W(1,15  ))    
      CALL IOVXXX(W(1,6   ),W(1,15  ),W(1,12  ),GWF ,AMP(5   ))            
      CALL FVOXXX(W(1,11  ),W(1,2   ),GZU ,ZERO    ,ZERO    ,W(1,16  ))    
      CALL JIOXXX(W(1,3   ),W(1,16  ),GWF ,WMASS   ,WWIDTH  ,W(1,17  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,17  ),GWF ,AMP(6   ))            
      CALL FVIXXX(W(1,7   ),W(1,1   ),GAD ,ZERO    ,ZERO    ,W(1,18  ))    
      CALL JIOXXX(W(1,18  ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,19  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,19  ),GWF ,AMP(7   ))            
      CALL JIOXXX(W(1,7   ),W(1,11  ),GWF ,WMASS   ,WWIDTH  ,W(1,20  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,20  ),GWF ,AMP(8   ))            
      CALL FVOXXX(W(1,4   ),W(1,2   ),GZU ,ZERO    ,ZERO    ,W(1,21  ))    
      CALL JIOXXX(W(1,3   ),W(1,21  ),GWF ,WMASS   ,WWIDTH  ,W(1,22  ))    
      CALL JVVXXX(W(1,22  ),W(1,1   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     23  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,23  ),GWF ,AMP(9   ))            
      CALL IOVXXX(W(1,6   ),W(1,10  ),W(1,22  ),GWF ,AMP(10  ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GAD ,ZERO    ,ZERO    ,W(1,24  ))    
      CALL JIOXXX(W(1,24  ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,25  ))    
      CALL FVOXXX(W(1,5   ),W(1,25  ),GWF ,ZERO    ,ZERO    ,W(1,26  ))    
      CALL IOVXXX(W(1,6   ),W(1,26  ),W(1,2   ),GZN ,AMP(11  ))            
      CALL JVVXXX(W(1,25  ),W(1,2   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     27  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,27  ),GWF ,AMP(12  ))            
      CALL IOVXXX(W(1,6   ),W(1,15  ),W(1,25  ),GWF ,AMP(13  ))            
      CALL FVOXXX(W(1,21  ),W(1,1   ),GAU ,ZERO    ,ZERO    ,W(1,28  ))    
      CALL JIOXXX(W(1,3   ),W(1,28  ),GWF ,WMASS   ,WWIDTH  ,W(1,29  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,29  ),GWF ,AMP(14  ))            
      CALL FVIXXX(W(1,24  ),W(1,2   ),GZD ,ZERO    ,ZERO    ,W(1,30  ))    
      CALL JIOXXX(W(1,30  ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,31  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,31  ),GWF ,AMP(15  ))            
      CALL JIOXXX(W(1,24  ),W(1,21  ),GWF ,WMASS   ,WWIDTH  ,W(1,32  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,32  ),GWF ,AMP(16  ))            
      CALL JIOXXX(W(1,3   ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,33  ))    
      CALL FVOXXX(W(1,10  ),W(1,33  ),GWF ,ZERO    ,ZERO    ,W(1,34  ))    
      CALL IOVXXX(W(1,6   ),W(1,34  ),W(1,2   ),GZN ,AMP(17  ))            
      CALL JVVXXX(W(1,33  ),W(1,2   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     35  ))                                                          
      CALL JVVXXX(W(1,35  ),W(1,1   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     36  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,36  ),GWF ,AMP(18  ))            
      CALL FVOXXX(W(1,10  ),W(1,2   ),GZL ,ZERO    ,ZERO    ,W(1,37  ))    
      CALL IOVXXX(W(1,6   ),W(1,37  ),W(1,33  ),GWF ,AMP(19  ))            
      CALL IOVXXX(W(1,6   ),W(1,10  ),W(1,35  ),GWF ,AMP(20  ))            
      CALL JVVXXX(W(1,33  ),W(1,1   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     38  ))                                                          
      CALL FVOXXX(W(1,5   ),W(1,38  ),GWF ,ZERO    ,ZERO    ,W(1,39  ))    
      CALL IOVXXX(W(1,6   ),W(1,39  ),W(1,2   ),GZN ,AMP(21  ))            
      CALL JVVXXX(W(1,38  ),W(1,2   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     40  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,40  ),GWF ,AMP(22  ))            
      CALL FVOXXX(W(1,15  ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,41  ))    
      CALL IOVXXX(W(1,6   ),W(1,41  ),W(1,33  ),GWF ,AMP(23  ))            
      CALL IOVXXX(W(1,6   ),W(1,15  ),W(1,38  ),GWF ,AMP(24  ))            
      CALL JW3WXX(W(1,2   ),W(1,33  ),W(1,1   ),GWWZ ,GWWA ,WMASS   ,      
     &     WWIDTH  ,W(1,42  ))                                             
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,42  ),GWF ,AMP(25  ))            
      AZ(MU,NU) = -AMP(   1)-AMP(   2)-AMP(   3)-AMP(   4)-AMP(   5)
     &             -AMP(   6)-AMP(   7)-AMP(   8)-AMP(   9)-AMP(  10)
     &             -AMP(  11)-AMP(  12)-AMP(  13)-AMP(  14)-AMP(  15)
     &             -AMP(  16)-AMP(  17)-AMP(  18)-AMP(  19)-AMP(  20)
     &             -AMP(  21)-AMP(  22)-AMP(  23)-AMP(  24)-AMP(  25)

      elseif (procid.eq.WPWMhadjj) then

      CALL FVIXXX(W(1,3   ),W(1,2   ),GZL ,ZERO    ,ZERO    ,W(1,7   ))    
      CALL JIOXXX(W(1,7   ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,8   ))    
      CALL FVOXXX(W(1,5   ),W(1,8   ),GWF ,ZERO    ,ZERO    ,W(1,9   ))    
      CALL IOVXXX(W(1,6   ),W(1,9   ),W(1,1   ),GAU ,AMP(1   ))            
      CALL JVVXXX(W(1,8   ),W(1,1   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     10  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,10  ),GWF ,AMP(2   ))            
      CALL FVOXXX(W(1,5   ),W(1,1   ),GAD ,ZERO    ,ZERO    ,W(1,11  ))    
      CALL IOVXXX(W(1,6   ),W(1,11  ),W(1,8   ),GWF ,AMP(3   ))            
      CALL FVIXXX(W(1,7   ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,12  ))    
      CALL JIOXXX(W(1,12  ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,13  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,13  ),GWF ,AMP(4   ))            
      CALL FVOXXX(W(1,4   ),W(1,2   ),GZN ,ZERO    ,ZERO    ,W(1,14  ))    
      CALL JIOXXX(W(1,3   ),W(1,14  ),GWF ,WMASS   ,WWIDTH  ,W(1,15  ))    
      CALL FVOXXX(W(1,5   ),W(1,15  ),GWF ,ZERO    ,ZERO    ,W(1,16  ))    
      CALL IOVXXX(W(1,6   ),W(1,16  ),W(1,1   ),GAU ,AMP(5   ))            
      CALL JVVXXX(W(1,15  ),W(1,1   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     17  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,17  ),GWF ,AMP(6   ))            
      CALL IOVXXX(W(1,6   ),W(1,11  ),W(1,15  ),GWF ,AMP(7   ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,18  ))    
      CALL JIOXXX(W(1,18  ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,19  ))    
      CALL FVOXXX(W(1,5   ),W(1,19  ),GWF ,ZERO    ,ZERO    ,W(1,20  ))    
      CALL IOVXXX(W(1,6   ),W(1,20  ),W(1,2   ),GZU ,AMP(8   ))            
      CALL JVVXXX(W(1,19  ),W(1,2   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     21  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,21  ),GWF ,AMP(9   ))            
      CALL FVOXXX(W(1,5   ),W(1,2   ),GZD ,ZERO    ,ZERO    ,W(1,22  ))    
      CALL IOVXXX(W(1,6   ),W(1,22  ),W(1,19  ),GWF ,AMP(10  ))            
      CALL FVIXXX(W(1,18  ),W(1,2   ),GZL ,ZERO    ,ZERO    ,W(1,23  ))    
      CALL JIOXXX(W(1,23  ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,24  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,24  ),GWF ,AMP(11  ))            
      CALL JIOXXX(W(1,18  ),W(1,14  ),GWF ,WMASS   ,WWIDTH  ,W(1,25  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,25  ),GWF ,AMP(12  ))            
      CALL JIOXXX(W(1,3   ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,26  ))    
      CALL JVVXXX(W(1,26  ),W(1,2   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     27  ))                                                          
      CALL FVOXXX(W(1,5   ),W(1,27  ),GWF ,ZERO    ,ZERO    ,W(1,28  ))    
      CALL IOVXXX(W(1,6   ),W(1,28  ),W(1,1   ),GAU ,AMP(13  ))            
      CALL FVOXXX(W(1,11  ),W(1,26  ),GWF ,ZERO    ,ZERO    ,W(1,29  ))    
      CALL IOVXXX(W(1,6   ),W(1,29  ),W(1,2   ),GZU ,AMP(14  ))            
      CALL JVVXXX(W(1,27  ),W(1,1   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     30  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,30  ),GWF ,AMP(15  ))            
      CALL FVOXXX(W(1,11  ),W(1,2   ),GZD ,ZERO    ,ZERO    ,W(1,31  ))    
      CALL IOVXXX(W(1,6   ),W(1,31  ),W(1,26  ),GWF ,AMP(16  ))            
      CALL IOVXXX(W(1,6   ),W(1,11  ),W(1,27  ),GWF ,AMP(17  ))            
      CALL FVOXXX(W(1,22  ),W(1,26  ),GWF ,ZERO    ,ZERO    ,W(1,32  ))    
      CALL IOVXXX(W(1,6   ),W(1,32  ),W(1,1   ),GAU ,AMP(18  ))            
      CALL JVVXXX(W(1,26  ),W(1,1   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     33  ))                                                          
      CALL FVOXXX(W(1,5   ),W(1,33  ),GWF ,ZERO    ,ZERO    ,W(1,34  ))    
      CALL IOVXXX(W(1,6   ),W(1,34  ),W(1,2   ),GZU ,AMP(19  ))            
      CALL JVVXXX(W(1,33  ),W(1,2   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     35  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,35  ),GWF ,AMP(20  ))            
      CALL FVOXXX(W(1,22  ),W(1,1   ),GAD ,ZERO    ,ZERO    ,W(1,36  ))    
      CALL IOVXXX(W(1,6   ),W(1,36  ),W(1,26  ),GWF ,AMP(21  ))            
      CALL IOVXXX(W(1,6   ),W(1,22  ),W(1,33  ),GWF ,AMP(22  ))            
      CALL JW3WXX(W(1,2   ),W(1,26  ),W(1,1   ),GWWZ ,GWWA ,WMASS   ,      
     &     WWIDTH  ,W(1,37  ))                                             
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,37  ),GWF ,AMP(23  ))            
      CALL FVOXXX(W(1,5   ),W(1,26  ),GWF ,ZERO    ,ZERO    ,W(1,38  ))    
      CALL FVOXXX(W(1,38  ),W(1,2   ),GZU ,ZERO    ,ZERO    ,W(1,39  ))    
      CALL IOVXXX(W(1,6   ),W(1,39  ),W(1,1   ),GAU ,AMP(24  ))            
      CALL FVOXXX(W(1,38  ),W(1,1   ),GAU ,ZERO    ,ZERO    ,W(1,40  ))    
      CALL IOVXXX(W(1,6   ),W(1,40  ),W(1,2   ),GZU ,AMP(25  ))            
      AZ(MU,NU) = -AMP(   1)-AMP(   2)-AMP(   3)-AMP(   4)-AMP(   5)
     &             -AMP(   6)-AMP(   7)-AMP(   8)-AMP(   9)-AMP(  10)
     &             -AMP(  11)-AMP(  12)-AMP(  13)-AMP(  14)-AMP(  15)
     &             -AMP(  16)-AMP(  17)-AMP(  18)-AMP(  19)-AMP(  20)
     &             -AMP(  21)-AMP(  22)-AMP(  23)-AMP(  24)-AMP(  25)

      else

       AZ(MU,NU) = dcmplx(0.,0.)

      endif

         ENDDO
      ENDDO

      END


      subroutine ZZtoWW_mg2_had(Q1,Q2,P,zztens)
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : z z -> e+ ve mu- vm~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=  32) 
      INTEGER    NWAVEFUNCS
      PARAMETER (NWAVEFUNCS=  51) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 Q1(0:3), Q2(0:3), P(0:3,3:6)
      complex*16 zztens(0:3,0:3)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J, MU, NU
      COMPLEX*16 AMP(NGRAPHS)
      COMPLEX*16 W(6,NWAVEFUNCS),WX1(6),WX2(6)
      INTEGER NHEL(3:6), ihel
      DATA (NHEL(IHEL),IHEL=3,6) / +1,-1,-1,+1/
      Save NHEL
c  lepton spinors and W+- polarization vectors
#include "tensor.inc"

      logical lswitch
      parameter (lswitch=.true.)
C  
C GLOBAL VARIABLES
C  
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/process.inc"
C                                                                          
C ----------
C BEGIN CODE
C ----------
      if (lswitch) then      ! read in lepton an W states from common block
         do i=1,6
            w(i,3) = wep(i)
            w(i,4) = wve(i)
            w(i,5) = wmu(i)
            w(i,6) = wvm(i)
         enddo
      endif

      DO MU = 0,3
c         print*," mu = ",mu,":",(q1(mu)+q2(mu))/(p(mu,3)+p(mu,4)+p(mu,5)+p(mu,6))
         DO NU = 0,3
      CALL VCARTX(Q1(0),ZMASS ,Zwidth, MU,-1,W(1,1   ))        
      CALL VCARTX(Q2(0),ZMASS ,Zwidth, NU,-1,W(1,2   ))

      if (procid.eq.WPhadWMjj) then

      CALL FVIXXX(W(1,3   ),W(1,2   ),GZD ,ZERO    ,ZERO    ,W(1,7   ))    
      CALL JIOXXX(W(1,7   ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,8   ))    
      CALL FVOXXX(W(1,5   ),W(1,8   ),GWF ,ZERO    ,ZERO    ,W(1,9   ))    
      CALL IOVXXX(W(1,6   ),W(1,9   ),W(1,1   ),GZN ,AMP(1   ))            
      CALL JVVXXX(W(1,8   ),W(1,1   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     10  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,10  ),GWF ,AMP(2   ))            
      CALL FVOXXX(W(1,5   ),W(1,1   ),GZL ,ZERO    ,ZERO    ,W(1,11  ))    
      CALL IOVXXX(W(1,6   ),W(1,11  ),W(1,8   ),GWF ,AMP(3   ))            
      CALL FVOXXX(W(1,4   ),W(1,1   ),GZU ,ZERO    ,ZERO    ,W(1,12  ))    
      CALL JIOXXX(W(1,3   ),W(1,12  ),GWF ,WMASS   ,WWIDTH  ,W(1,13  ))    
      CALL FVOXXX(W(1,5   ),W(1,13  ),GWF ,ZERO    ,ZERO    ,W(1,14  ))    
      CALL IOVXXX(W(1,6   ),W(1,14  ),W(1,2   ),GZN ,AMP(4   ))            
      CALL JVVXXX(W(1,13  ),W(1,2   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     15  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,15  ),GWF ,AMP(5   ))            
      CALL FVOXXX(W(1,5   ),W(1,2   ),GZL ,ZERO    ,ZERO    ,W(1,16  ))    
      CALL IOVXXX(W(1,6   ),W(1,16  ),W(1,13  ),GWF ,AMP(6   ))            
      CALL FVOXXX(W(1,12  ),W(1,2   ),GZU ,ZERO    ,ZERO    ,W(1,17  ))    
      CALL JIOXXX(W(1,3   ),W(1,17  ),GWF ,WMASS   ,WWIDTH  ,W(1,18  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,18  ),GWF ,AMP(7   ))            
      CALL FVIXXX(W(1,7   ),W(1,1   ),GZD ,ZERO    ,ZERO    ,W(1,19  ))    
      CALL JIOXXX(W(1,19  ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,20  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,20  ),GWF ,AMP(8   ))            
      CALL JIOXXX(W(1,7   ),W(1,12  ),GWF ,WMASS   ,WWIDTH  ,W(1,21  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,21  ),GWF ,AMP(9   ))            
      CALL FVOXXX(W(1,4   ),W(1,2   ),GZU ,ZERO    ,ZERO    ,W(1,22  ))    
      CALL JIOXXX(W(1,3   ),W(1,22  ),GWF ,WMASS   ,WWIDTH  ,W(1,23  ))    
      CALL FVOXXX(W(1,5   ),W(1,23  ),GWF ,ZERO    ,ZERO    ,W(1,24  ))    
      CALL IOVXXX(W(1,6   ),W(1,24  ),W(1,1   ),GZN ,AMP(10  ))            
      CALL JVVXXX(W(1,23  ),W(1,1   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     25  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,25  ),GWF ,AMP(11  ))            
      CALL IOVXXX(W(1,6   ),W(1,11  ),W(1,23  ),GWF ,AMP(12  ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GZD ,ZERO    ,ZERO    ,W(1,26  ))    
      CALL JIOXXX(W(1,26  ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,27  ))    
      CALL FVOXXX(W(1,5   ),W(1,27  ),GWF ,ZERO    ,ZERO    ,W(1,28  ))    
      CALL IOVXXX(W(1,6   ),W(1,28  ),W(1,2   ),GZN ,AMP(13  ))            
      CALL JVVXXX(W(1,27  ),W(1,2   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     29  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,29  ),GWF ,AMP(14  ))            
      CALL IOVXXX(W(1,6   ),W(1,16  ),W(1,27  ),GWF ,AMP(15  ))            
      CALL FVOXXX(W(1,22  ),W(1,1   ),GZU ,ZERO    ,ZERO    ,W(1,30  ))    
      CALL JIOXXX(W(1,3   ),W(1,30  ),GWF ,WMASS   ,WWIDTH  ,W(1,31  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,31  ),GWF ,AMP(16  ))            
      CALL FVIXXX(W(1,26  ),W(1,2   ),GZD ,ZERO    ,ZERO    ,W(1,32  ))    
      CALL JIOXXX(W(1,32  ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,33  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,33  ),GWF ,AMP(17  ))            
      CALL JIOXXX(W(1,26  ),W(1,22  ),GWF ,WMASS   ,WWIDTH  ,W(1,34  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,34  ),GWF ,AMP(18  ))            
      CALL JIOXXX(W(1,3   ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,35  ))    
      CALL JVVXXX(W(1,35  ),W(1,2   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     36  ))                                                          
      CALL FVOXXX(W(1,5   ),W(1,36  ),GWF ,ZERO    ,ZERO    ,W(1,37  ))    
      CALL IOVXXX(W(1,6   ),W(1,37  ),W(1,1   ),GZN ,AMP(19  ))            
      CALL FVOXXX(W(1,11  ),W(1,35  ),GWF ,ZERO    ,ZERO    ,W(1,38  ))    
      CALL IOVXXX(W(1,6   ),W(1,38  ),W(1,2   ),GZN ,AMP(20  ))            
      CALL JVVXXX(W(1,36  ),W(1,1   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     39  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,39  ),GWF ,AMP(21  ))            
      CALL FVOXXX(W(1,11  ),W(1,2   ),GZL ,ZERO    ,ZERO    ,W(1,40  ))    
      CALL IOVXXX(W(1,6   ),W(1,40  ),W(1,35  ),GWF ,AMP(22  ))            
      CALL IOVXXX(W(1,6   ),W(1,11  ),W(1,36  ),GWF ,AMP(23  ))            
      CALL FVOXXX(W(1,16  ),W(1,35  ),GWF ,ZERO    ,ZERO    ,W(1,41  ))    
      CALL IOVXXX(W(1,6   ),W(1,41  ),W(1,1   ),GZN ,AMP(24  ))            
      CALL JVVXXX(W(1,35  ),W(1,1   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     42  ))                                                          
      CALL FVOXXX(W(1,5   ),W(1,42  ),GWF ,ZERO    ,ZERO    ,W(1,43  ))    
      CALL IOVXXX(W(1,6   ),W(1,43  ),W(1,2   ),GZN ,AMP(25  ))            
      CALL JVVXXX(W(1,42  ),W(1,2   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     44  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,44  ),GWF ,AMP(26  ))            
      CALL FVOXXX(W(1,16  ),W(1,1   ),GZL ,ZERO    ,ZERO    ,W(1,45  ))    
      CALL IOVXXX(W(1,6   ),W(1,45  ),W(1,35  ),GWF ,AMP(27  ))            
      CALL IOVXXX(W(1,6   ),W(1,16  ),W(1,42  ),GWF ,AMP(28  ))            
      CALL JW3WXX(W(1,1   ),W(1,35  ),W(1,2   ),GWWZ ,GWWZ ,WMASS   ,      
     &     WWIDTH  ,W(1,46  ))                                             
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,46  ),GWF ,AMP(29  ))            
      CALL FVOXXX(W(1,5   ),W(1,35  ),GWF ,ZERO    ,ZERO    ,W(1,47  ))    
      CALL FVOXXX(W(1,47  ),W(1,2   ),GZN ,ZERO    ,ZERO    ,W(1,48  ))    
      CALL IOVXXX(W(1,6   ),W(1,48  ),W(1,1   ),GZN ,AMP(30  ))            
      CALL FVOXXX(W(1,47  ),W(1,1   ),GZN ,ZERO    ,ZERO    ,W(1,49  ))    
      CALL IOVXXX(W(1,6   ),W(1,49  ),W(1,2   ),GZN ,AMP(31  ))            
      CALL HVVXXX(W(1,1   ),W(1,2   ),GZZH ,HMASS   ,HWIDTH  ,W(1,         
     &     50  ))                                                          
      CALL JVSXXX(W(1,35  ),W(1,50  ),GWWH ,WMASS   ,WWIDTH  ,W(1,         
     &     51  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,51  ),GWF ,AMP(32  ))            
      zztens(MU,NU) = -AMP(   1)-AMP(   2)-AMP(   3)-AMP(   4)-AMP(   5)
     &             -AMP(   6)-AMP(   7)-AMP(   8)-AMP(   9)-AMP(  10)
     &             -AMP(  11)-AMP(  12)-AMP(  13)-AMP(  14)-AMP(  15)
     &             -AMP(  16)-AMP(  17)-AMP(  18)-AMP(  19)-AMP(  20)
     &             -AMP(  21)-AMP(  22)-AMP(  23)-AMP(  24)-AMP(  25)
     &             -AMP(  26)-AMP(  27)-AMP(  28)-AMP(  29)-AMP(  30)
     &             -AMP(  31)-AMP(  32)

      elseif (procid.eq.WPWMhadjj) then

      CALL FVIXXX(W(1,3   ),W(1,2   ),GZL ,ZERO    ,ZERO    ,W(1,7   ))    
      CALL JIOXXX(W(1,7   ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,8   ))    
      CALL FVOXXX(W(1,5   ),W(1,8   ),GWF ,ZERO    ,ZERO    ,W(1,9   ))    
      CALL IOVXXX(W(1,6   ),W(1,9   ),W(1,1   ),GZU ,AMP(1   ))            
      CALL JVVXXX(W(1,8   ),W(1,1   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     10  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,10  ),GWF ,AMP(2   ))            
      CALL FVOXXX(W(1,5   ),W(1,1   ),GZD ,ZERO    ,ZERO    ,W(1,11  ))    
      CALL IOVXXX(W(1,6   ),W(1,11  ),W(1,8   ),GWF ,AMP(3   ))            
      CALL FVOXXX(W(1,4   ),W(1,1   ),GZN ,ZERO    ,ZERO    ,W(1,12  ))    
      CALL JIOXXX(W(1,3   ),W(1,12  ),GWF ,WMASS   ,WWIDTH  ,W(1,13  ))    
      CALL FVOXXX(W(1,5   ),W(1,13  ),GWF ,ZERO    ,ZERO    ,W(1,14  ))    
      CALL IOVXXX(W(1,6   ),W(1,14  ),W(1,2   ),GZU ,AMP(4   ))            
      CALL JVVXXX(W(1,13  ),W(1,2   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     15  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,15  ),GWF ,AMP(5   ))            
      CALL FVOXXX(W(1,5   ),W(1,2   ),GZD ,ZERO    ,ZERO    ,W(1,16  ))    
      CALL IOVXXX(W(1,6   ),W(1,16  ),W(1,13  ),GWF ,AMP(6   ))            
      CALL FVOXXX(W(1,12  ),W(1,2   ),GZN ,ZERO    ,ZERO    ,W(1,17  ))    
      CALL JIOXXX(W(1,3   ),W(1,17  ),GWF ,WMASS   ,WWIDTH  ,W(1,18  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,18  ),GWF ,AMP(7   ))            
      CALL FVIXXX(W(1,7   ),W(1,1   ),GZL ,ZERO    ,ZERO    ,W(1,19  ))    
      CALL JIOXXX(W(1,19  ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,20  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,20  ),GWF ,AMP(8   ))            
      CALL JIOXXX(W(1,7   ),W(1,12  ),GWF ,WMASS   ,WWIDTH  ,W(1,21  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,21  ),GWF ,AMP(9   ))            
      CALL FVOXXX(W(1,4   ),W(1,2   ),GZN ,ZERO    ,ZERO    ,W(1,22  ))    
      CALL JIOXXX(W(1,3   ),W(1,22  ),GWF ,WMASS   ,WWIDTH  ,W(1,23  ))    
      CALL FVOXXX(W(1,5   ),W(1,23  ),GWF ,ZERO    ,ZERO    ,W(1,24  ))    
      CALL IOVXXX(W(1,6   ),W(1,24  ),W(1,1   ),GZU ,AMP(10  ))            
      CALL JVVXXX(W(1,23  ),W(1,1   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     25  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,25  ),GWF ,AMP(11  ))            
      CALL IOVXXX(W(1,6   ),W(1,11  ),W(1,23  ),GWF ,AMP(12  ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GZL ,ZERO    ,ZERO    ,W(1,26  ))    
      CALL JIOXXX(W(1,26  ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,27  ))    
      CALL FVOXXX(W(1,5   ),W(1,27  ),GWF ,ZERO    ,ZERO    ,W(1,28  ))    
      CALL IOVXXX(W(1,6   ),W(1,28  ),W(1,2   ),GZU ,AMP(13  ))            
      CALL JVVXXX(W(1,27  ),W(1,2   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     29  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,29  ),GWF ,AMP(14  ))            
      CALL IOVXXX(W(1,6   ),W(1,16  ),W(1,27  ),GWF ,AMP(15  ))            
      CALL FVOXXX(W(1,22  ),W(1,1   ),GZN ,ZERO    ,ZERO    ,W(1,30  ))    
      CALL JIOXXX(W(1,3   ),W(1,30  ),GWF ,WMASS   ,WWIDTH  ,W(1,31  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,31  ),GWF ,AMP(16  ))            
      CALL FVIXXX(W(1,26  ),W(1,2   ),GZL ,ZERO    ,ZERO    ,W(1,32  ))    
      CALL JIOXXX(W(1,32  ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,33  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,33  ),GWF ,AMP(17  ))            
      CALL JIOXXX(W(1,26  ),W(1,22  ),GWF ,WMASS   ,WWIDTH  ,W(1,34  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,34  ),GWF ,AMP(18  ))            
      CALL JIOXXX(W(1,3   ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,35  ))    
      CALL JVVXXX(W(1,35  ),W(1,2   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     36  ))                                                          
      CALL FVOXXX(W(1,5   ),W(1,36  ),GWF ,ZERO    ,ZERO    ,W(1,37  ))    
      CALL IOVXXX(W(1,6   ),W(1,37  ),W(1,1   ),GZU ,AMP(19  ))            
      CALL FVOXXX(W(1,11  ),W(1,35  ),GWF ,ZERO    ,ZERO    ,W(1,38  ))    
      CALL IOVXXX(W(1,6   ),W(1,38  ),W(1,2   ),GZU ,AMP(20  ))            
      CALL JVVXXX(W(1,36  ),W(1,1   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     39  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,39  ),GWF ,AMP(21  ))            
      CALL FVOXXX(W(1,11  ),W(1,2   ),GZD ,ZERO    ,ZERO    ,W(1,40  ))    
      CALL IOVXXX(W(1,6   ),W(1,40  ),W(1,35  ),GWF ,AMP(22  ))            
      CALL IOVXXX(W(1,6   ),W(1,11  ),W(1,36  ),GWF ,AMP(23  ))            
      CALL FVOXXX(W(1,16  ),W(1,35  ),GWF ,ZERO    ,ZERO    ,W(1,41  ))    
      CALL IOVXXX(W(1,6   ),W(1,41  ),W(1,1   ),GZU ,AMP(24  ))            
      CALL JVVXXX(W(1,35  ),W(1,1   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     42  ))                                                          
      CALL FVOXXX(W(1,5   ),W(1,42  ),GWF ,ZERO    ,ZERO    ,W(1,43  ))    
      CALL IOVXXX(W(1,6   ),W(1,43  ),W(1,2   ),GZU ,AMP(25  ))            
      CALL JVVXXX(W(1,42  ),W(1,2   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     44  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,44  ),GWF ,AMP(26  ))            
      CALL FVOXXX(W(1,16  ),W(1,1   ),GZD ,ZERO    ,ZERO    ,W(1,45  ))    
      CALL IOVXXX(W(1,6   ),W(1,45  ),W(1,35  ),GWF ,AMP(27  ))            
      CALL IOVXXX(W(1,6   ),W(1,16  ),W(1,42  ),GWF ,AMP(28  ))            
      CALL JW3WXX(W(1,1   ),W(1,35  ),W(1,2   ),GWWZ ,GWWZ ,WMASS   ,      
     &     WWIDTH  ,W(1,46  ))                                             
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,46  ),GWF ,AMP(29  ))            
      CALL FVOXXX(W(1,5   ),W(1,35  ),GWF ,ZERO    ,ZERO    ,W(1,47  ))    
      CALL FVOXXX(W(1,47  ),W(1,2   ),GZU ,ZERO    ,ZERO    ,W(1,48  ))    
      CALL IOVXXX(W(1,6   ),W(1,48  ),W(1,1   ),GZU ,AMP(30  ))            
      CALL FVOXXX(W(1,47  ),W(1,1   ),GZU ,ZERO    ,ZERO    ,W(1,49  ))    
      CALL IOVXXX(W(1,6   ),W(1,49  ),W(1,2   ),GZU ,AMP(31  ))            
      CALL HVVXXX(W(1,1   ),W(1,2   ),GZZH ,HMASS   ,HWIDTH  ,W(1,         
     &     50  ))                                                          
      CALL JVSXXX(W(1,35  ),W(1,50  ),GWWH ,WMASS   ,WWIDTH  ,W(1,         
     &     51  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,51  ),GWF ,AMP(32  ))            
      zztens(MU,NU) = -AMP(   1)-AMP(   2)-AMP(   3)-AMP(   4)-AMP(   5)
     &             -AMP(   6)-AMP(   7)-AMP(   8)-AMP(   9)-AMP(  10)
     &             -AMP(  11)-AMP(  12)-AMP(  13)-AMP(  14)-AMP(  15)
     &             -AMP(  16)-AMP(  17)-AMP(  18)-AMP(  19)-AMP(  20)
     &             -AMP(  21)-AMP(  22)-AMP(  23)-AMP(  24)-AMP(  25)
     &             -AMP(  26)-AMP(  27)-AMP(  28)-AMP(  29)-AMP(  30)
     &             -AMP(  31)-AMP(  32)

      else

       zztens(MU,NU) = dcmplx(0.,0.)

      endif

          ENDDO
      ENDDO

      END


      subroutine WWtoWW_mg2_had(Q1,Q2,P,wwtens)
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : w- w+ -> e+ ve mu- vm~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=  31) 
      INTEGER    NWAVEFUNCS
      PARAMETER (NWAVEFUNCS=  54)
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 Q1(0:3), Q2(0:3), P(0:3,3:6)
      complex*16 wwtens(0:3,0:3)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J, MU, NU
      COMPLEX*16 AMP(NGRAPHS)
      COMPLEX*16 W(6,NWAVEFUNCS)
      INTEGER NHEL(3:6), ihel
      DATA (NHEL(IHEL),IHEL=3,6) / +1,-1,-1,+1/
      Save NHEL
c  lepton spinors and W+- polarization vectors
#include "tensor.inc"

      logical lswitch
      parameter (lswitch=.true.)
C  
C GLOBAL VARIABLES
C  
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/process.inc"
C                                                                          
C ----------
C BEGIN CODE
C ----------
      if (lswitch) then      ! read in lepton an W states from common block
         do i=1,6
            w(i,3) = wep(i)
            w(i,4) = wve(i)
            w(i,5) = wmu(i)
            w(i,6) = wvm(i)
         enddo
      endif

      DO MU = 0,3
         DO NU = 0,3
      CALL VCARTX(Q1(0),wMASS ,wwidth, MU,-1,W(1,1   ))        
      CALL VCARTX(Q2(0),wMASS ,wwidth, NU,-1,W(1,2   )) 

      if (procid.eq.WPhadWMjj) then

      CALL FVIXXX(W(1,3   ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,7   ))    
      CALL JIOXXX(W(1,7   ),W(1,4   ),GAU ,ZERO    ,AWIDTH  ,W(1,8   ))    
      CALL FVOXXX(W(1,5   ),W(1,8   ),GAL ,ZERO    ,ZERO    ,W(1,9   ))    
      CALL IOVXXX(W(1,6   ),W(1,9   ),W(1,1   ),GWF ,AMP(1   ))            
      CALL JIOXXX(W(1,7   ),W(1,4   ),GZU ,ZMASS   ,ZWIDTH  ,W(1,10  ))    
      CALL FVOXXX(W(1,5   ),W(1,10  ),GZL ,ZERO    ,ZERO    ,W(1,11  ))    
      CALL IOVXXX(W(1,6   ),W(1,11  ),W(1,1   ),GWF ,AMP(2   ))            
      CALL JVVXXX(W(1,1   ),W(1,8   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     12  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,12  ),GWF ,AMP(3   ))            
      CALL JVVXXX(W(1,1   ),W(1,10  ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     13  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,13  ),GWF ,AMP(4   ))            
      CALL FVOXXX(W(1,5   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,14  ))    
      CALL IOVXXX(W(1,6   ),W(1,14  ),W(1,10  ),GZN ,AMP(5   ))            
      CALL FVIXXX(W(1,7   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,15  ))    
      CALL JIOXXX(W(1,15  ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,16  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,16  ),GWF ,AMP(6   ))            
      CALL FVOXXX(W(1,4   ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,17  ))    
      CALL JIOXXX(W(1,3   ),W(1,17  ),GAD ,ZERO    ,AWIDTH  ,W(1,18  ))    
      CALL FVOXXX(W(1,5   ),W(1,18  ),GAL ,ZERO    ,ZERO    ,W(1,19  ))    
      CALL IOVXXX(W(1,6   ),W(1,19  ),W(1,1   ),GWF ,AMP(7   ))            
      CALL JIOXXX(W(1,3   ),W(1,17  ),GZD ,ZMASS   ,ZWIDTH  ,W(1,20  ))    
      CALL FVOXXX(W(1,5   ),W(1,20  ),GZL ,ZERO    ,ZERO    ,W(1,21  ))    
      CALL IOVXXX(W(1,6   ),W(1,21  ),W(1,1   ),GWF ,AMP(8   ))            
      CALL JVVXXX(W(1,1   ),W(1,18  ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     22  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,22  ),GWF ,AMP(9   ))            
      CALL JVVXXX(W(1,1   ),W(1,20  ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     23  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,23  ),GWF ,AMP(10  ))            
      CALL IOVXXX(W(1,6   ),W(1,14  ),W(1,20  ),GZN ,AMP(11  ))            
      CALL FVOXXX(W(1,17  ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,24  ))    
      CALL JIOXXX(W(1,3   ),W(1,24  ),GWF ,WMASS   ,WWIDTH  ,W(1,25  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,25  ),GWF ,AMP(12  ))            
      CALL JIOXXX(W(1,3   ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,26  ))    
      CALL JVVXXX(W(1,2   ),W(1,26  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     27  ))                                                          
      CALL FVOXXX(W(1,5   ),W(1,27  ),GAL ,ZERO    ,ZERO    ,W(1,28  ))    
      CALL IOVXXX(W(1,6   ),W(1,28  ),W(1,1   ),GWF ,AMP(13  ))            
      CALL JVVXXX(W(1,2   ),W(1,26  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     29  ))                                                          
      CALL FVOXXX(W(1,5   ),W(1,29  ),GZL ,ZERO    ,ZERO    ,W(1,30  ))    
      CALL IOVXXX(W(1,6   ),W(1,30  ),W(1,1   ),GWF ,AMP(14  ))            
      CALL JVVXXX(W(1,1   ),W(1,27  ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     31  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,31  ),GWF ,AMP(15  ))            
      CALL JVVXXX(W(1,1   ),W(1,29  ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     32  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,32  ),GWF ,AMP(16  ))            
      CALL HVVXXX(W(1,2   ),W(1,26  ),GWWH ,HMASS   ,HWIDTH  ,W(1,         
     &     33  ))                                                          
      CALL JVSXXX(W(1,1   ),W(1,33  ),GWWH ,WMASS   ,WWIDTH  ,W(1,         
     &     34  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,34  ),GWF ,AMP(17  ))            
      CALL FVOXXX(W(1,14  ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,35  ))    
      CALL IOVXXX(W(1,6   ),W(1,35  ),W(1,26  ),GWF ,AMP(18  ))            
      CALL IOVXXX(W(1,6   ),W(1,14  ),W(1,29  ),GZN ,AMP(19  ))            
      CALL JVVXXX(W(1,2   ),W(1,1   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     36  ))                                                          
      CALL FVOXXX(W(1,4   ),W(1,36  ),GAU ,ZERO    ,ZERO    ,W(1,37  ))    
      CALL JIOXXX(W(1,3   ),W(1,37  ),GWF ,WMASS   ,WWIDTH  ,W(1,38  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,38  ),GWF ,AMP(20  ))            
      CALL JVVXXX(W(1,2   ),W(1,1   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     39  ))                                                          
      CALL FVOXXX(W(1,4   ),W(1,39  ),GZU ,ZERO    ,ZERO    ,W(1,40  ))    
      CALL JIOXXX(W(1,3   ),W(1,40  ),GWF ,WMASS   ,WWIDTH  ,W(1,41  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,41  ),GWF ,AMP(21  ))            
      CALL FVIXXX(W(1,3   ),W(1,36  ),GAD ,ZERO    ,ZERO    ,W(1,42  ))    
      CALL JIOXXX(W(1,42  ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,43  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,43  ),GWF ,AMP(22  ))            
      CALL FVIXXX(W(1,3   ),W(1,39  ),GZD ,ZERO    ,ZERO    ,W(1,44  ))    
      CALL JIOXXX(W(1,44  ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,45  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,45  ),GWF ,AMP(23  ))            
      CALL JWWWXX(W(1,1   ),W(1,2   ),W(1,26  ),GWWA ,GWWZ ,WMASS   ,      
     &     WWIDTH  ,W(1,46  ))                                             
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,46  ),GWF ,AMP(24  ))            
      CALL FVOXXX(W(1,5   ),W(1,26  ),GWF ,ZERO    ,ZERO    ,W(1,47  ))    
      CALL FVOXXX(W(1,47  ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,48  ))    
      CALL IOVXXX(W(1,6   ),W(1,48  ),W(1,1   ),GWF ,AMP(25  ))            
      CALL JVVXXX(W(1,26  ),W(1,36  ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     49  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,49  ),GWF ,AMP(26  ))            
      CALL JVVXXX(W(1,26  ),W(1,39  ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     50  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,50  ),GWF ,AMP(27  ))            
      CALL HVVXXX(W(1,2   ),W(1,1   ),GWWH ,HMASS   ,HWIDTH  ,W(1,         
     &     51  ))                                                          
      CALL JVSXXX(W(1,26  ),W(1,51  ),GWWH ,WMASS   ,WWIDTH  ,W(1,         
     &     52  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,52  ),GWF ,AMP(28  ))            
      CALL IOVXXX(W(1,6   ),W(1,47  ),W(1,39  ),GZN ,AMP(29  ))            
      CALL FVOXXX(W(1,5   ),W(1,36  ),GAL ,ZERO    ,ZERO    ,W(1,53  ))    
      CALL IOVXXX(W(1,6   ),W(1,53  ),W(1,26  ),GWF ,AMP(30  ))            
      CALL FVOXXX(W(1,5   ),W(1,39  ),GZL ,ZERO    ,ZERO    ,W(1,54  ))    
      CALL IOVXXX(W(1,6   ),W(1,54  ),W(1,26  ),GWF ,AMP(31  ))            
      wwtens(MU,NU) = -AMP(   1)-AMP(   2)-AMP(   3)-AMP(   4)-AMP(   5)
     &             -AMP(   6)-AMP(   7)-AMP(   8)-AMP(   9)-AMP(  10)
     &             -AMP(  11)-AMP(  12)-AMP(  13)-AMP(  14)-AMP(  15)
     &             -AMP(  16)-AMP(  17)-AMP(  18)-AMP(  19)-AMP(  20)
     &             -AMP(  21)-AMP(  22)-AMP(  23)-AMP(  24)-AMP(  25)
     &             -AMP(  26)-AMP(  27)-AMP(  28)-AMP(  29)-AMP(  30)
     &             -AMP(  31)

      elseif (procid.eq.WPWMhadjj) then

      CALL FVIXXX(W(1,3   ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,7   ))    
      CALL JIOXXX(W(1,7   ),W(1,4   ),GZN ,ZMASS   ,ZWIDTH  ,W(1,8   ))    
      CALL FVOXXX(W(1,5   ),W(1,8   ),GZD ,ZERO    ,ZERO    ,W(1,9   ))    
      CALL IOVXXX(W(1,6   ),W(1,9   ),W(1,1   ),GWF ,AMP(1   ))            
      CALL JVVXXX(W(1,1   ),W(1,8   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     10  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,10  ),GWF ,AMP(2   ))            
      CALL FVOXXX(W(1,5   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,11  ))    
      CALL IOVXXX(W(1,6   ),W(1,11  ),W(1,8   ),GZU ,AMP(3   ))            
      CALL FVIXXX(W(1,7   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,12  ))    
      CALL JIOXXX(W(1,12  ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,13  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,13  ),GWF ,AMP(4   ))            
      CALL FVOXXX(W(1,4   ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,14  ))    
      CALL JIOXXX(W(1,3   ),W(1,14  ),GAL ,ZERO    ,AWIDTH  ,W(1,15  ))    
      CALL FVOXXX(W(1,5   ),W(1,15  ),GAD ,ZERO    ,ZERO    ,W(1,16  ))    
      CALL IOVXXX(W(1,6   ),W(1,16  ),W(1,1   ),GWF ,AMP(5   ))            
      CALL JIOXXX(W(1,3   ),W(1,14  ),GZL ,ZMASS   ,ZWIDTH  ,W(1,17  ))    
      CALL FVOXXX(W(1,5   ),W(1,17  ),GZD ,ZERO    ,ZERO    ,W(1,18  ))    
      CALL IOVXXX(W(1,6   ),W(1,18  ),W(1,1   ),GWF ,AMP(6   ))            
      CALL JVVXXX(W(1,1   ),W(1,15  ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     19  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,19  ),GWF ,AMP(7   ))            
      CALL JVVXXX(W(1,1   ),W(1,17  ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     20  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,20  ),GWF ,AMP(8   ))            
      CALL IOVXXX(W(1,6   ),W(1,11  ),W(1,15  ),GAU ,AMP(9   ))            
      CALL IOVXXX(W(1,6   ),W(1,11  ),W(1,17  ),GZU ,AMP(10  ))            
      CALL FVOXXX(W(1,14  ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,21  ))    
      CALL JIOXXX(W(1,3   ),W(1,21  ),GWF ,WMASS   ,WWIDTH  ,W(1,22  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,22  ),GWF ,AMP(11  ))            
      CALL JIOXXX(W(1,3   ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,23  ))    
      CALL JVVXXX(W(1,2   ),W(1,23  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     24  ))                                                          
      CALL FVOXXX(W(1,5   ),W(1,24  ),GAD ,ZERO    ,ZERO    ,W(1,25  ))    
      CALL IOVXXX(W(1,6   ),W(1,25  ),W(1,1   ),GWF ,AMP(12  ))            
      CALL JVVXXX(W(1,2   ),W(1,23  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     26  ))                                                          
      CALL FVOXXX(W(1,5   ),W(1,26  ),GZD ,ZERO    ,ZERO    ,W(1,27  ))    
      CALL IOVXXX(W(1,6   ),W(1,27  ),W(1,1   ),GWF ,AMP(13  ))            
      CALL JVVXXX(W(1,1   ),W(1,24  ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     28  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,28  ),GWF ,AMP(14  ))            
      CALL JVVXXX(W(1,1   ),W(1,26  ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     29  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,29  ),GWF ,AMP(15  ))            
      CALL HVVXXX(W(1,2   ),W(1,23  ),GWWH ,HMASS   ,HWIDTH  ,W(1,         
     &     30  ))                                                          
      CALL JVSXXX(W(1,1   ),W(1,30  ),GWWH ,WMASS   ,WWIDTH  ,W(1,         
     &     31  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,31  ),GWF ,AMP(16  ))            
      CALL FVOXXX(W(1,11  ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,32  ))    
      CALL IOVXXX(W(1,6   ),W(1,32  ),W(1,23  ),GWF ,AMP(17  ))            
      CALL IOVXXX(W(1,6   ),W(1,11  ),W(1,24  ),GAU ,AMP(18  ))            
      CALL IOVXXX(W(1,6   ),W(1,11  ),W(1,26  ),GZU ,AMP(19  ))            
      CALL JVVXXX(W(1,2   ),W(1,1   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     33  ))                                                          
      CALL FVOXXX(W(1,4   ),W(1,33  ),GZN ,ZERO    ,ZERO    ,W(1,34  ))    
      CALL JIOXXX(W(1,3   ),W(1,34  ),GWF ,WMASS   ,WWIDTH  ,W(1,35  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,35  ),GWF ,AMP(20  ))            
      CALL JVVXXX(W(1,2   ),W(1,1   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     36  ))                                                          
      CALL FVIXXX(W(1,3   ),W(1,36  ),GAL ,ZERO    ,ZERO    ,W(1,37  ))    
      CALL JIOXXX(W(1,37  ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,38  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,38  ),GWF ,AMP(21  ))            
      CALL FVIXXX(W(1,3   ),W(1,33  ),GZL ,ZERO    ,ZERO    ,W(1,39  ))    
      CALL JIOXXX(W(1,39  ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,40  ))    
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,40  ),GWF ,AMP(22  ))            
      CALL JWWWXX(W(1,1   ),W(1,2   ),W(1,23  ),GWWA ,GWWZ ,WMASS   ,      
     &     WWIDTH  ,W(1,41  ))                                             
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,41  ),GWF ,AMP(23  ))            
      CALL FVOXXX(W(1,5   ),W(1,23  ),GWF ,ZERO    ,ZERO    ,W(1,42  ))    
      CALL FVOXXX(W(1,42  ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,43  ))    
      CALL IOVXXX(W(1,6   ),W(1,43  ),W(1,1   ),GWF ,AMP(24  ))            
      CALL JVVXXX(W(1,23  ),W(1,36  ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     44  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,44  ),GWF ,AMP(25  ))            
      CALL JVVXXX(W(1,23  ),W(1,33  ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     45  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,45  ),GWF ,AMP(26  ))            
      CALL HVVXXX(W(1,2   ),W(1,1   ),GWWH ,HMASS   ,HWIDTH  ,W(1,         
     &     46  ))                                                          
      CALL JVSXXX(W(1,23  ),W(1,46  ),GWWH ,WMASS   ,WWIDTH  ,W(1,         
     &     47  ))                                                          
      CALL IOVXXX(W(1,6   ),W(1,5   ),W(1,47  ),GWF ,AMP(27  ))            
      CALL IOVXXX(W(1,6   ),W(1,42  ),W(1,36  ),GAU ,AMP(28  ))            
      CALL IOVXXX(W(1,6   ),W(1,42  ),W(1,33  ),GZU ,AMP(29  ))            
      CALL FVOXXX(W(1,5   ),W(1,36  ),GAD ,ZERO    ,ZERO    ,W(1,48  ))    
      CALL IOVXXX(W(1,6   ),W(1,48  ),W(1,23  ),GWF ,AMP(30  ))            
      CALL FVOXXX(W(1,5   ),W(1,33  ),GZD ,ZERO    ,ZERO    ,W(1,49  ))    
      CALL IOVXXX(W(1,6   ),W(1,49  ),W(1,23  ),GWF ,AMP(31  ))            
      wwtens(MU,NU) = -AMP(   1)-AMP(   2)-AMP(   3)-AMP(   4)-AMP(   5)
     &             -AMP(   6)-AMP(   7)-AMP(   8)-AMP(   9)-AMP(  10)
     &             -AMP(  11)-AMP(  12)-AMP(  13)-AMP(  14)-AMP(  15)
     &             -AMP(  16)-AMP(  17)-AMP(  18)-AMP(  19)-AMP(  20)
     &             -AMP(  21)-AMP(  22)-AMP(  23)-AMP(  24)-AMP(  25)
     &             -AMP(  26)-AMP(  27)-AMP(  28)-AMP(  29)-AMP(  30)
     &             -AMP(  31)

      else

       wwtens(MU,NU) = dcmplx(0.,0.)

      endif

          ENDDO
      ENDDO

      END


      subroutine ZtoWW_mg2_had(P,z)
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : z -> e+ ve mu- vm~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=   5) 
      INTEGER    NWAVEFUNCS
      PARAMETER (NWAVEFUNCS=  13) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 Q(0:3), P(0:3,2:5)
      complex*16 z(0:5)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J, mu
      COMPLEX*16 AMP(NGRAPHS)
      COMPLEX*16 W(6,NWAVEFUNCS)
      INTEGER NHEL(2:5), ihel
      DATA (NHEL(IHEL),IHEL=2,5) / +1,-1,-1,+1/
      Save NHEL
c  lepton spinors and W+- polarization vectors
#include "tensor.inc"

      logical lswitch
      parameter (lswitch=.true.)
C  
C GLOBAL VARIABLES
C  
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/process.inc"
C                                                                          
C ----------
C BEGIN CODE
C ----------
      if (lswitch) then      ! read in lepton an W states from common block
         do i=1,6
            w(i,2) = wep(i)
            w(i,3) = wve(i)
            w(i,4) = wmu(i)
            w(i,5) = wvm(i)
         enddo
      endif

      do mu = 0,3
         q(mu) = p(mu,2)+p(mu,3)+p(mu,4)+p(mu,5)
      enddo
c
      do mu = 0,3
      CALL VCARTX(Q(0     ),ZMASS ,Zwidth,MU,-1,W(1,1   )) 

      if (procid.eq.WPhadWMjj) then

      CALL JIOXXX(W(1,2   ),W(1,3   ),GWF ,WMASS   ,WWIDTH  ,W(1,6   ))    
      CALL FVOXXX(W(1,4   ),W(1,6   ),GWF ,ZERO    ,ZERO    ,W(1,7   ))    
      CALL IOVXXX(W(1,5   ),W(1,7   ),W(1,1   ),GZN ,AMP(1   ))            
      CALL JVVXXX(W(1,6   ),W(1,1   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     8   ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,8   ),GWF ,AMP(2   ))            
      CALL FVOXXX(W(1,4   ),W(1,1   ),GZL ,ZERO    ,ZERO    ,W(1,9   ))    
      CALL IOVXXX(W(1,5   ),W(1,9   ),W(1,6   ),GWF ,AMP(3   ))            
      CALL FVOXXX(W(1,3   ),W(1,1   ),GZU ,ZERO    ,ZERO    ,W(1,10  ))    
      CALL JIOXXX(W(1,2   ),W(1,10  ),GWF ,WMASS   ,WWIDTH  ,W(1,11  ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,11  ),GWF ,AMP(4   ))            
      CALL FVIXXX(W(1,2   ),W(1,1   ),GZD ,ZERO    ,ZERO    ,W(1,12  ))    
      CALL JIOXXX(W(1,12  ),W(1,3   ),GWF ,WMASS   ,WWIDTH  ,W(1,13  ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,13  ),GWF ,AMP(5   ))            
      Z(MU) = -AMP(   1)-AMP(   2)-AMP(   3)-AMP(   4)-AMP(   5)

      elseif (procid.eq.WPWMhadjj) then

      CALL JIOXXX(W(1,2   ),W(1,3   ),GWF ,WMASS   ,WWIDTH  ,W(1,6   ))    
      CALL FVOXXX(W(1,4   ),W(1,6   ),GWF ,ZERO    ,ZERO    ,W(1,7   ))    
      CALL IOVXXX(W(1,5   ),W(1,7   ),W(1,1   ),GZU ,AMP(1   ))            
      CALL JVVXXX(W(1,6   ),W(1,1   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     8   ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,8   ),GWF ,AMP(2   ))            
      CALL FVOXXX(W(1,4   ),W(1,1   ),GZD ,ZERO    ,ZERO    ,W(1,9   ))    
      CALL IOVXXX(W(1,5   ),W(1,9   ),W(1,6   ),GWF ,AMP(3   ))            
      CALL FVOXXX(W(1,3   ),W(1,1   ),GZN ,ZERO    ,ZERO    ,W(1,10  ))    
      CALL JIOXXX(W(1,2   ),W(1,10  ),GWF ,WMASS   ,WWIDTH  ,W(1,11  ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,11  ),GWF ,AMP(4   ))            
      CALL FVIXXX(W(1,2   ),W(1,1   ),GZL ,ZERO    ,ZERO    ,W(1,12  ))    
      CALL JIOXXX(W(1,12  ),W(1,3   ),GWF ,WMASS   ,WWIDTH  ,W(1,13  ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,13  ),GWF ,AMP(5   ))            
      Z(MU) = -AMP(   1)-AMP(   2)-AMP(   3)-AMP(   4)-AMP(   5)

      else

       Z(MU) = dcmplx(0.,0.)

      endif

      enddo
      z(4) = dcmplx(q(0),q(3))
      z(5) = dcmplx(q(1),q(2))
c
      return
      END


      subroutine AtoWW_mg2_had(P,a)
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : a -> e+ ve mu- vm~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=   4) 
      INTEGER    NWAVEFUNCS
      PARAMETER (NWAVEFUNCS=  12) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 Q(0:3), P(0:3,2:5)
      complex*16 a(0:5)
C  
C LOCAL VARIABLES 
C
      INTEGER I,J, mu
      COMPLEX*16 AMP(NGRAPHS)
      COMPLEX*16 W(6,NWAVEFUNCS)
      INTEGER NHEL(2:5), ihel
      DATA (NHEL(IHEL),IHEL=2,5) /+1,-1,-1,+1/
      Save NHEL
c  lepton spinors and W+- polarization vectors
#include "tensor.inc"

      logical lswitch
      parameter (lswitch=.true.)
C  
C GLOBAL VARIABLES
C  
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/process.inc"
C                                                                          
C ----------
C BEGIN CODE
C ----------
      if (lswitch) then      ! read in lepton an W states from common block
         do i=1,6
            w(i,2) = wep(i)
            w(i,3) = wve(i)
            w(i,4) = wmu(i)
            w(i,5) = wvm(i)
         enddo
      endif

      do mu = 0,3
         q(mu) = p(mu,2)+p(mu,3)+p(mu,4)+p(mu,5)
      enddo
c
      do mu = 0,3
      CALL VCARTX(Q(0),ZERO ,ZERO, MU,-1,W(1,1   ))

      if (procid.eq.WPhadWMjj) then

      CALL JIOXXX(W(1,2   ),W(1,3   ),GWF ,WMASS   ,WWIDTH  ,W(1,6   ))    
      CALL JVVXXX(W(1,6   ),W(1,1   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     7   ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,7   ),GWF ,AMP(1   ))            
      CALL FVOXXX(W(1,4   ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,8   ))    
      CALL IOVXXX(W(1,5   ),W(1,8   ),W(1,6   ),GWF ,AMP(2   ))            
      CALL FVOXXX(W(1,3   ),W(1,1   ),GAU ,ZERO    ,ZERO    ,W(1,9   ))    
      CALL JIOXXX(W(1,2   ),W(1,9   ),GWF ,WMASS   ,WWIDTH  ,W(1,10  ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,10  ),GWF ,AMP(3   ))            
      CALL FVIXXX(W(1,2   ),W(1,1   ),GAD ,ZERO    ,ZERO    ,W(1,11  ))    
      CALL JIOXXX(W(1,11  ),W(1,3   ),GWF ,WMASS   ,WWIDTH  ,W(1,12  ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,12  ),GWF ,AMP(4   ))            
      A(MU) = -AMP(   1)-AMP(   2)-AMP(   3)-AMP(   4)

      elseif (procid.eq.WPWMhadjj) then

      CALL JIOXXX(W(1,2   ),W(1,3   ),GWF ,WMASS   ,WWIDTH  ,W(1,6   ))    
      CALL FVOXXX(W(1,4   ),W(1,6   ),GWF ,ZERO    ,ZERO    ,W(1,7   ))    
      CALL IOVXXX(W(1,5   ),W(1,7   ),W(1,1   ),GAU ,AMP(1   ))            
      CALL JVVXXX(W(1,6   ),W(1,1   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     8   ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,8   ),GWF ,AMP(2   ))            
      CALL FVOXXX(W(1,4   ),W(1,1   ),GAD ,ZERO    ,ZERO    ,W(1,9   ))    
      CALL IOVXXX(W(1,5   ),W(1,9   ),W(1,6   ),GWF ,AMP(3   ))            
      CALL FVIXXX(W(1,2   ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,10  ))    
      CALL JIOXXX(W(1,10  ),W(1,3   ),GWF ,WMASS   ,WWIDTH  ,W(1,11  ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,11  ),GWF ,AMP(4   ))            
      A(MU) = -AMP(   1)-AMP(   2)-AMP(   3)-AMP(   4)

      else

       A(MU) = dcmplx(0.,0.)

      endif

      enddo
      a(4) = dcmplx(q(0),q(3))
      a(5) = dcmplx(q(1),q(2))
c
      END


      subroutine WVtoWM_mg2_had(id,Q,v,wmatens,wmztens)
C  Input:   Q(*)        one of the two incoming boson momenta 
c           id = 1      Q = q1 = W- momentum
c           id = 2      Q = q2 = gamma/Z momentum
c  OUPTPUT: wmatens(mu,nu)  lepton tensor for w gamma--> mu- nubar
c           wmztens(mu,nu)  lepton tensor for w Z    --> mu- nubar
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : w- z -> mu- vm~  
C and           w- a -> mu- vm~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=   3) 
      INTEGER    NWAVEFUNCS
      PARAMETER (NWAVEFUNCS=   7)
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      integer id
      real*8 q(0:3), v(0:3,4)
      complex*16 wmatens(0:3,0:3), wmztens(0:3,0:3)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J, MU, NU
      real*8 q2sq, Q1(0:3), Q2(0:3)
      COMPLEX*16 AMP(NGRAPHS)
      COMPLEX*16 W(6,NWAVEFUNCS), denz
c  lepton spinors and W+- polarization vectors
#include "tensor.inc"

      logical ldebug
      parameter (ldebug=.false.)
C  
C GLOBAL VARIABLES
C  
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/process.inc"
C                                                                          
C ----------
C BEGIN CODE
C ----------
      if (id.eq.1) then
         q2(0) = dreal(wm(5))
         q2(1) = dreal(wm(6))
         q2(2) = dimag(wm(6))
         q2(3) = dimag(wm(5))
         do mu = 0,3
            q1(mu) = q(mu)
            q2(mu) = q2(mu)-q(mu)
         enddo
      elseif (id.eq.2) then
         q1(0) = dreal(wm(5))
         q1(1) = dreal(wm(6))
         q1(2) = dimag(wm(6))
         q1(3) = dimag(wm(5))
         do mu = 0,3
            q2(mu) = q(mu)
            q1(mu) = q1(mu)-q(mu)
         enddo
      endif

      do i=1,6                  ! read in lepton states from common block
         w(i,3) = wmu(i)
         w(i,4) = wvm(i)
      enddo

      q2sq = q2(0)**2-q2(1)**2-q2(2)**2-q2(3)**2
      denz = q2sq/dcmplx(q2sq-zmass**2,zmass*zwidth)
      DO MU = 0,3
      CALL VCARTX(Q1(0),WMASS ,Wwidth, MU,-1,W(1,1))        
      DO NU = 0,3

      if (procid.eq.WPhadWMjj) then

      CALL VCARTX(Q2(0),ZMASS  ,ZWIDTH  , NU,-1,W(1,2))
      CALL FVOXXX(W(1,3   ),W(1,2   ),GZL ,ZERO    ,ZERO    ,W(1,5   ))    
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,1   ),GWF ,AMP(1   ))            
      CALL FVOXXX(W(1,3   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,6   ))    
      CALL IOVXXX(W(1,4   ),W(1,6   ),W(1,2   ),GZN ,AMP(2   ))            
      CALL JVVXXX(W(1,1   ),W(1,2   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     7   ))                                                          
      CALL IOVXXX(W(1,4   ),W(1,3   ),W(1,7   ),GWF ,AMP(3   ))            
      wmztens(mu,nu) = -AMP(   1)-AMP(   2)-AMP(   3)

      CALL VCARTX(Q2(0),ZERO  ,ZERO  , NU,-1,W(1,2))
      CALL FVOXXX(W(1,3   ),W(1,2   ),GAL ,ZERO    ,ZERO    ,W(1,5   ))    
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,1   ),GWF ,AMP(1   ))            
      CALL JVVXXX(W(1,1   ),W(1,2   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     6   ))                                                          
      CALL IOVXXX(W(1,4   ),W(1,3   ),W(1,6   ),GWF ,AMP(2   ))            
      wmatens(mu,nu) = -AMP(   1)-AMP(   2)

      elseif (procid.eq.WPWMhadjj) then

      CALL VCARTX(Q2(0),ZMASS  ,ZWIDTH  , NU,-1,W(1,2))
      CALL FVOXXX(W(1,3   ),W(1,2   ),GZD ,ZERO    ,ZERO    ,W(1,5   ))    
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,1   ),GWF ,AMP(1   ))            
      CALL FVOXXX(W(1,3   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,6   ))    
      CALL IOVXXX(W(1,4   ),W(1,6   ),W(1,2   ),GZU ,AMP(2   ))            
      CALL JVVXXX(W(1,1   ),W(1,2   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     7   ))                                                          
      CALL IOVXXX(W(1,4   ),W(1,3   ),W(1,7   ),GWF ,AMP(3   ))            
      wmztens(mu,nu) = -AMP(   1)-AMP(   2)-AMP(   3)

      CALL VCARTX(Q2(0),ZERO  ,ZERO  , NU,-1,W(1,2))
      CALL FVOXXX(W(1,3   ),W(1,2   ),GAD ,ZERO    ,ZERO    ,W(1,5   ))    
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,1   ),GWF ,AMP(1   ))            
      CALL FVOXXX(W(1,3   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,6   ))    
      CALL IOVXXX(W(1,4   ),W(1,6   ),W(1,2   ),GAU ,AMP(2   ))            
      CALL JVVXXX(W(1,1   ),W(1,2   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     7   ))                                                          
      CALL IOVXXX(W(1,4   ),W(1,3   ),W(1,7   ),GWF ,AMP(3   ))            
      wmatens(mu,nu) = -AMP(   1)-AMP(   2)-AMP(   3)

      else

       wmatens(mu,nu) = dcmplx(0.,0.)
       wmztens(mu,nu) = dcmplx(0.,0.)

      endif

         enddo
      enddo

      end


      subroutine WVtoWP_mg2_had(id,Q,v,wpatens,wpztens)
C  Input:   Q(*)        one of the two incoming boson momenta 
c           id = 1      Q = q1 = W+ momentum
c           id = 2      Q = q2 = gamma/Z momentum
c  OUPTPUT: wpatens(mu,nu)  lepton tensor for w gamma--> e+nu
c           wpztens(mu,nu)  lepton tensor for w Z    --> e+nu
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : w+(q1) z(q2) -> e+ ve  
C FOR PROCESS : w+(q1) a(q2) -> e+ ve  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=   3) 
      INTEGER    NWAVEFUNCS
      PARAMETER (NWAVEFUNCS=   7) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS
C
      integer id
      real*8 q(0:3), v(0:3,4)
      complex*16 wpatens(0:3,0:3), wpztens(0:3,0:3)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J, MU, NU
      real*8 q2sq, Q1(0:3), Q2(0:3)
      COMPLEX*16 AMP(NGRAPHS)
      COMPLEX*16 W(6,NWAVEFUNCS), denz
c  lepton spinors and W+- polarization vectors
#include "tensor.inc"

      logical ldebug
      parameter (ldebug=.false.)
C  
C GLOBAL VARIABLES
C  
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/process.inc"
C                                                                          
C ----------
C BEGIN CODE
C ----------
      
      if (id.eq.1) then
         q2(0) = dreal(wp(5))
         q2(1) = dreal(wp(6))
         q2(2) = dimag(wp(6))
         q2(3) = dimag(wp(5))
         do mu = 0,3
            q1(mu) = q(mu)
            q2(mu) = q2(mu)-q(mu)
         enddo
      elseif (id.eq.2) then
         q1(0) = dreal(wp(5))
         q1(1) = dreal(wp(6))
         q1(2) = dimag(wp(6))
         q1(3) = dimag(wp(5))
         do mu = 0,3
            q2(mu) = q(mu)
            q1(mu) = q1(mu)-q(mu)
         enddo
      endif

      do i=1,6                  ! read in lepton states from common block
         w(i,3) = wep(i)
         w(i,4) = wve(i)
      enddo

      q2sq = q2(0)**2-q2(1)**2-q2(2)**2-q2(3)**2
      denz = q2sq/dcmplx(q2sq-zmass**2,zmass*zwidth)
      DO MU = 0,3
      DO NU = 0,3
      CALL VCARTX(Q1(0),WMASS ,Wwidth, MU,-1,W(1,1))        

      if (procid.eq.WPhadWMjj) then

      CALL VCARTX(Q2(0),ZMASS  ,ZWIDTH  , NU,-1,W(1,2))
      CALL FVIXXX(W(1,3   ),W(1,2   ),GZD ,ZERO    ,ZERO    ,W(1,5   ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,1   ),GWF ,AMP(1   ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,6   ))    
      CALL IOVXXX(W(1,6   ),W(1,4   ),W(1,2   ),GZU ,AMP(2   ))            
      CALL JVVXXX(W(1,2   ),W(1,1   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     7   ))                                                          
      CALL IOVXXX(W(1,3   ),W(1,4   ),W(1,7   ),GWF ,AMP(3   ))            
      wpztens(mu,nu) = -AMP(   1)-AMP(   2)-AMP(   3)

      CALL VCARTX(Q2(0),ZERO  ,ZERO  , NU,-1,W(1,2))
      CALL FVIXXX(W(1,3   ),W(1,2   ),GAD ,ZERO    ,ZERO    ,W(1,5   ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,1   ),GWF ,AMP(1   ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,6   ))    
      CALL IOVXXX(W(1,6   ),W(1,4   ),W(1,2   ),GAU ,AMP(2   ))            
      CALL JVVXXX(W(1,2   ),W(1,1   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     7   ))                                                          
      CALL IOVXXX(W(1,3   ),W(1,4   ),W(1,7   ),GWF ,AMP(3   ))            
      wpatens(mu,nu) = -amp(1)-amp(2)-amp(3)

      elseif (procid.eq.WPWMhadjj) then

      CALL VCARTX(Q2(0),ZMASS  ,ZWIDTH  , NU,-1,W(1,2))
      CALL FVIXXX(W(1,3   ),W(1,2   ),GZL ,ZERO    ,ZERO    ,W(1,5   ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,1   ),GWF ,AMP(1   ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,6   ))    
      CALL IOVXXX(W(1,6   ),W(1,4   ),W(1,2   ),GZN ,AMP(2   ))            
      CALL JVVXXX(W(1,2   ),W(1,1   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     7   ))                                                          
      CALL IOVXXX(W(1,3   ),W(1,4   ),W(1,7   ),GWF ,AMP(3   ))            
      wpztens(mu,nu) = -AMP(   1)-AMP(   2)-AMP(   3)

      CALL VCARTX(Q2(0),ZERO  ,ZERO  , NU,-1,W(1,2))
      CALL FVIXXX(W(1,3   ),W(1,2   ),GAL ,ZERO    ,ZERO    ,W(1,5   ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,1   ),GWF ,AMP(1   ))            
      CALL JVVXXX(W(1,2   ),W(1,1   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     6   ))                                                          
      CALL IOVXXX(W(1,3   ),W(1,4   ),W(1,6   ),GWF ,AMP(2   ))            
      wpatens(mu,nu) = -amp(1)-amp(2)

      else

       wpatens(mu,nu) = dcmplx(0.,0.)
       wpztens(mu,nu) = dcmplx(0.,0.)

      endif

      enddo
      enddo

      end

