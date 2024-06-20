      subroutine qqwzqq_mg2(p,pbar,v,sign, nlo,bos,
     1                          uucs,ddcs,udcc,udss) 
      
      implicit none                

#include "VBFNLO/utilities/global.inc"
#include "tensorw.inc"
#include "VBFNLO/utilities/process.inc"

      double precision pbar(0:3,4+max_v)
      double precision uucs,ddcs,udcc,udss
      double precision uucs2,ddcs2,udcc2,udss2
      double precision p(0:3,max_p),v(0:3,max_v)
      integer  sign(4+max_v), nlo, i, j,mu,nu,bos
      real*8 q12(0:4,3), q34(0:4,3)
      complex*16 mg2_tens_1d(0:5)
      complex*16 mg2_tens_2d(0:3,0:3)
      complex*16 mg2_tens_2d_comb(0:3,0:3,4)
      double precision sumtest

      integer h
      common /hcount / h

      call Calc_Momentum_Transfer(p, v, q12,q34,1)

      j = 1
         
#ifdef WITH_MADGRAPH

         SELECT CASE(procid)
         CASE(WPhadZjj, WPZhadjj)            ! semileptonic decay modes

            call anomal_formfactor(q12(0,j),q34(0,j),qe(0),qu(0))
            call vwtovwp_mg2_had(q12(0,j),q34(0,j),v,h,zwwv(0,0,j),awwv(0,0,j)) !q12 = V
            call vwtovwp_mg2_had(q34(0,j),q12(0,j),v,h,wzwv(0,0,j),wawv(0,0,j)) !q12 = W+
            call VVtomm_mg2_had(1,h,q12(0,j),v,aauu(0,0,2,j),azuu(0,0,2,j),zauu(0,0,2,j),zzuu(0,0,2,j))    	 !emit W on lower
            call VVtomm_mg2_had(2,h,q34(0,j),v,aauu(0,0,1,j),azuu(0,0,1,j),zauu(0,0,1,j),zzuu(0,0,1,j))    	 !emit W on upper
            call WWtomm_mg2_had(1,h,q34(0,j),v,wwuu(0,0,1,j))     !emit V->e nu on upper
            call WWtomm_mg2_had(1,h,q12(0,j),v,wwuu(0,0,2,j))     !emit V->e nu on lower
            call WPVtoln_mg2_had(1,q34(0,j),v,awen(0,0,1,j),zwen(0,0,1,j))  !emit V->mu+mu- on upper
            call WPVtoln_mg2_had(2,q12(0,j),v,awen(0,0,2,j),zwen(0,0,2,j))  !emit V->mu+mu- on lower
            call WPVtoln_mg2_had(1,q12(0,j),v,waen(0,0,2,j),wzen(0,0,2,j))  !emit V->mu+mu- on lower
            call WPVtoln_mg2_had(2,q34(0,j),v,waen(0,0,1,j),wzen(0,0,1,j))  !emit V->mu+mu- on upper  
            call wpto4l_mg2_had(v,h,wwv)

            call qqwpzqq(pbar,sign, nlo,1,bos,uucs2,ddcs2,udcc2,udss2)

            print*, ""
            print*,"1 - uucs/uucs2 =   ",abs(1-uucs/uucs2)
            print*,"1 - ddcs/ddcs2 =   ",abs(1-ddcs/ddcs2)
            print*,"1 - udcc/udcc2 =   ",abs(1-udcc/udcc2)
            print*,"1 - udss/udss2 =   ",abs(1-udss/udss2)

            call anomal_formfactor(q12(0,j),q34(0,j),qe(0),qu(0))
            call vwtovwp_had(q12(0,j),q34(0,j),v,h,zwwv(0,0,j),awwv(0,0,j)) !q12 = V
            call vwtovwp_had(q34(0,j),q12(0,j),v,h,wzwv(0,0,j),wawv(0,0,j)) !q12 = W+
            call VVtomm_had(1,h,q12(0,j),v,aauu(0,0,2,j),azuu(0,0,2,j),zauu(0,0,2,j),zzuu(0,0,2,j))    	 !emit W on lower
            call VVtomm_had(2,h,q34(0,j),v,aauu(0,0,1,j),azuu(0,0,1,j),zauu(0,0,1,j),zzuu(0,0,1,j))    	 !emit W on upper
            call WWtomm_had(1,h,q34(0,j),v,wwuu(0,0,1,j))     !emit V->e nu on upper
            call WWtomm_had(1,h,q12(0,j),v,wwuu(0,0,2,j))     !emit V->e nu on lower
            call WPVtoln_had(1,q34(0,j),v,awen(0,0,1,j),zwen(0,0,1,j))  !emit V->mu+mu- on upper
            call WPVtoln_had(2,q12(0,j),v,awen(0,0,2,j),zwen(0,0,2,j))  !emit V->mu+mu- on lower
            call WPVtoln_had(1,q12(0,j),v,waen(0,0,2,j),wzen(0,0,2,j))  !emit V->mu+mu- on lower
            call WPVtoln_had(2,q34(0,j),v,waen(0,0,1,j),wzen(0,0,1,j))  !emit V->mu+mu- on upper  
            call wpto4l_had(v,h,wwv)


         do mu = 0,3
          do nu = 0,3
           mg2_tens_2d(mu,nu) = dcmplx(0.,0.)
           mg2_tens_2d_comb(mu,nu,1) = dcmplx(0.,0.)
           mg2_tens_2d_comb(mu,nu,2) = dcmplx(0.,0.)
           mg2_tens_2d_comb(mu,nu,3) = dcmplx(0.,0.)
           mg2_tens_2d_comb(mu,nu,4) = dcmplx(0.,0.)
          enddo
         enddo
         do i = 0,5
           mg2_tens_1d(i) = dcmplx(0.,0.)
         enddo


         print*, ""
         call wpto4l_mg2_had(v,h,mg2_tens_1d)
         sumtest = 0
         do i = 0,5
           sumtest = sumtest + abs( wwv(i) / mg2_tens_1d(i) - 1)
!           print*, i, wwv(i), mg2_tens_1d(i), abs( wwv(i) / mg2_tens_1d(i) - 1)
           mg2_tens_1d(i) = dcmplx(0.,0.)
         enddo
         print*, "tensor W+ to WZ:", sumtest


!         print*, ""
         call WPVtoln_mg2_had(1,q34(0,j),v,mg2_tens_2d_comb(0,0,1),mg2_tens_2d_comb(0,0,2))  !emit V->mu+mu- on upper
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( awen(mu,nu,1,j) / mg2_tens_2d_comb(mu,nu,1) - 1)
           ! print*, mu,nu, awen(mu,nu,1,j), mg2_tens_2d_comb(mu,nu,1), abs( awen(mu,nu,1,1) / mg2_tens_2d_comb(mu,nu,1) - 1)
           mg2_tens_2d_comb(mu,nu,1) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor W+Z to W+ 1 had:", sumtest
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( zwen(mu,nu,1,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           ! print*, mu,nu, zwen(mu,nu,1,j), mg2_tens_2d_comb(mu,nu,2), abs( zwen(mu,nu,1,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           mg2_tens_2d_comb(mu,nu,2) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor W+Z to W+ 1 had:", sumtest

!         print*, ""
         call WPVtoln_mg2_had(2,q12(0,j),v,mg2_tens_2d_comb(0,0,1),mg2_tens_2d_comb(0,0,2))  !emit V->mu+mu- on lower
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( awen(mu,nu,2,j) / mg2_tens_2d_comb(mu,nu,1) - 1)
           ! print*, mu,nu, awen(mu,nu,2,j), mg2_tens_2d_comb(mu,nu,1), abs( awen(mu,nu,2,1) / mg2_tens_2d_comb(mu,nu,1) - 1)
           mg2_tens_2d_comb(mu,nu,1) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor W+A to W+ 2 had:", sumtest
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( zwen(mu,nu,2,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           ! print*, mu,nu, zwen(mu,nu,2,j), mg2_tens_2d_comb(mu,nu,2), abs( zwen(mu,nu,2,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           mg2_tens_2d_comb(mu,nu,2) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor W+Z to W+ 2 had:", sumtest

!         print*, ""
         call WPVtoln_mg2_had(2,q34(0,j),v,mg2_tens_2d_comb(0,0,1),mg2_tens_2d_comb(0,0,2))  !emit V->mu+mu- on upper  
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( waen(mu,nu,1,j) / mg2_tens_2d_comb(mu,nu,1) - 1)
           ! print*, mu,nu, waen(mu,nu,1,j), mg2_tens_2d_comb(mu,nu,1), abs( waen(mu,nu,1,1) / mg2_tens_2d_comb(mu,nu,1) - 1)
           mg2_tens_2d_comb(mu,nu,1) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor W+A to W+ 4 had:", sumtest
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( wzen(mu,nu,1,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           ! print*, mu,nu, wzen(mu,nu,1,j), mg2_tens_2d_comb(mu,nu,2), abs( wzen(mu,nu,1,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           mg2_tens_2d_comb(mu,nu,2) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor W+Z to W+ 4 had:", sumtest

!         print*, ""
         call WPVtoln_mg2_had(1,q12(0,j),v,mg2_tens_2d_comb(0,0,1),mg2_tens_2d_comb(0,0,2))  !emit V->mu+mu- on lower
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( waen(mu,nu,2,j) / mg2_tens_2d_comb(mu,nu,1) - 1)
           ! print*, mu,nu, waen(mu,nu,2,j), mg2_tens_2d_comb(mu,nu,1), abs( waen(mu,nu,2,1) / mg2_tens_2d_comb(mu,nu,1) - 1)
           mg2_tens_2d_comb(mu,nu,1) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor W+A to W+ 3 had:", sumtest
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( wzen(mu,nu,2,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           ! print*, mu,nu, wzen(mu,nu,2,j), mg2_tens_2d_comb(mu,nu,2), abs( wzen(mu,nu,2,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           mg2_tens_2d_comb(mu,nu,2) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor W+Z to W+ 3 had:", sumtest


!         print*, ""
         call VVtomm_mg2_had(1,h,q12(0,j),v,mg2_tens_2d_comb(0,0,1),mg2_tens_2d_comb(0,0,2),
     &                                      mg2_tens_2d_comb(0,0,3),mg2_tens_2d_comb(0,0,4))    	 !emit W on lower
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( aauu(mu,nu,2,j) / mg2_tens_2d_comb(mu,nu,1) - 1)
           ! print*, mu,nu, aauu(mu,nu,2,j), mg2_tens_2d_comb(mu,nu,1), abs( aauu(mu,nu,2,j) / mg2_tens_2d_comb(mu,nu,1) - 1)
           mg2_tens_2d_comb(mu,nu,1) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor AA to VV 1:", sumtest
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( azuu(mu,nu,2,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           ! print*, mu,nu, azuu(mu,nu,2,j), mg2_tens_2d_comb(mu,nu,2), abs( azuu(mu,nu,2,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           mg2_tens_2d_comb(mu,nu,2) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor AZ to VV 1:", sumtest
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( zauu(mu,nu,2,j) / mg2_tens_2d_comb(mu,nu,3) - 1)
           ! print*, mu,nu, zauu(mu,nu,2,j), mg2_tens_2d_comb(mu,nu,3), abs( zauu(mu,nu,2,j) / mg2_tens_2d_comb(mu,nu,3) - 1)
           mg2_tens_2d_comb(mu,nu,3) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor ZA to VV 1:", sumtest
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( zzuu(mu,nu,2,j) / mg2_tens_2d_comb(mu,nu,4) - 1)
           ! print*, mu,nu, zzuu(mu,nu,2,j), mg2_tens_2d_comb(mu,nu,4), abs( zzuu(mu,nu,2,j) / mg2_tens_2d_comb(mu,nu,4) - 1)
           mg2_tens_2d_comb(mu,nu,4) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor ZZ to VV 1:", sumtest

!         print*, ""
         call VVtomm_mg2_had(2,h,q34(0,j),v,mg2_tens_2d_comb(0,0,1),mg2_tens_2d_comb(0,0,2),
     &                                      mg2_tens_2d_comb(0,0,3),mg2_tens_2d_comb(0,0,4))    	 !emit W on upper
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( aauu(mu,nu,1,j) / mg2_tens_2d_comb(mu,nu,1) - 1)
           ! print*, mu,nu, aauu(mu,nu,1,j), mg2_tens_2d_comb(mu,nu,1), abs( aauu(mu,nu,1,j) / mg2_tens_2d_comb(mu,nu,1) - 1)
           mg2_tens_2d_comb(mu,nu,1) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor AA to VV 2:", sumtest
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( azuu(mu,nu,1,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           ! print*, mu,nu, azuu(mu,nu,1,j), mg2_tens_2d_comb(mu,nu,2), abs( azuu(mu,nu,1,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           mg2_tens_2d_comb(mu,nu,2) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor AZ to VV 2:", sumtest
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( zauu(mu,nu,1,j) / mg2_tens_2d_comb(mu,nu,3) - 1)
           ! print*, mu,nu, zauu(mu,nu,1,j), mg2_tens_2d_comb(mu,nu,3), abs( zauu(mu,nu,1,j) / mg2_tens_2d_comb(mu,nu,3) - 1)
           mg2_tens_2d_comb(mu,nu,3) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor ZA to VV 2:", sumtest
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( zzuu(mu,nu,1,j) / mg2_tens_2d_comb(mu,nu,4) - 1)
           ! print*, mu,nu, zzuu(mu,nu,1,j), mg2_tens_2d_comb(mu,nu,4), abs( zzuu(mu,nu,1,j) / mg2_tens_2d_comb(mu,nu,4) - 1)
           mg2_tens_2d_comb(mu,nu,4) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor ZZ to VV 2:", sumtest


!         print*, ""
         call WWtomm_mg2_had(1,h,q34(0,j),v,mg2_tens_2d(0,0))     !emit V->e nu on upper
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( wwuu(mu,nu,1,j) / mg2_tens_2d(mu,nu) - 1)
           ! print*, mu,nu, wwuu(mu,nu,1,j), mg2_tens_2d(mu,nu), abs( wwuu(mu,nu,1,j) / mg2_tens_2d(mu,nu) - 1)
           mg2_tens_2d(mu,nu) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor WW to Z 1:", sumtest

!         print*, ""
         call WWtomm_mg2_had(1,h,q12(0,j),v,mg2_tens_2d(0,0))     !emit V->e nu on lower
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( wwuu(mu,nu,2,j) / mg2_tens_2d(mu,nu) - 1)
           ! print*, mu,nu, wwuu(mu,nu,2,j), mg2_tens_2d(mu,nu), abs( wwuu(mu,nu,2,j) / mg2_tens_2d(mu,nu) - 1)
           mg2_tens_2d(mu,nu) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor WW to Z 2:", sumtest


!         print*, ""
         call vwtovwp_mg2_had(q12(0,j),q34(0,j),v,h,mg2_tens_2d_comb(0,0,1),mg2_tens_2d_comb(0,0,2)) !q12 = V
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( zwwv(mu,nu,j) / mg2_tens_2d_comb(mu,nu,1) - 1)
           ! print*, mu,nu, zwwv(mu,nu,j), mg2_tens_2d_comb(mu,nu,1), abs( zwwv(mu,nu,j) / mg2_tens_2d_comb(mu,nu,1) - 1)
           mg2_tens_2d_comb(mu,nu,1) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor W+Z to W+V 1:", sumtest
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( awwv(mu,nu,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           ! print*, mu,nu, awwv(mu,nu,j), mg2_tens_2d_comb(mu,nu,2), abs( awwv(mu,nu,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           mg2_tens_2d_comb(mu,nu,2) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor W+A to W+V 1:", sumtest

!         print*, ""
         call vwtovwp_mg2_had(q34(0,j),q12(0,j),v,h,mg2_tens_2d_comb(0,0,1),mg2_tens_2d_comb(0,0,2)) !q12 = W+
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( wzwv(mu,nu,j) / mg2_tens_2d_comb(mu,nu,1) - 1)
           ! print*, mu,nu, wzwv(mu,nu,j), mg2_tens_2d_comb(mu,nu,1), abs( wzwv(mu,nu,j) / mg2_tens_2d_comb(mu,nu,1) - 1)
           mg2_tens_2d_comb(mu,nu,1) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor W+Z to W+V 2:", sumtest
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( wawv(mu,nu,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           ! print*, mu,nu, wawv(mu,nu,j), mg2_tens_2d_comb(mu,nu,2), abs( wawv(mu,nu,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           mg2_tens_2d_comb(mu,nu,2) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor W+A to W+V 2:", sumtest


         print*, ""
         print*, ""
         read(*,*)


         CASE(WMhadZjj, WMZhadjj)            ! semileptonic decay modes

            call anomal_formfactor(q12(0,j),q34(0,j),qe(0),qu(0))
            call vwtovwm_mg2_had(q12(0,j),q34(0,j),v,h,zwwv(0,0,j),awwv(0,0,j)) !q12 = V
            call vwtovwm_mg2_had(q34(0,j),q12(0,j),v,h,wzwv(0,0,j),wawv(0,0,j)) !q12 = W-
            call VVtomm_mg2_had(1,h,q12(0,j),v,aauu(0,0,2,j),azuu(0,0,2,j),zauu(0,0,2,j),zzuu(0,0,2,j))    	 !emit W on lower
            call VVtomm_mg2_had(2,h,q34(0,j),v,aauu(0,0,1,j),azuu(0,0,1,j),zauu(0,0,1,j),zzuu(0,0,1,j))    	 !emit W on upper
            call WWtomm_mg2_had(2,h,q34(0,j),v,wwuu(0,0,1,j))     !emit V->e nu on upper
            call WWtomm_mg2_had(2,h,q12(0,j),v,wwuu(0,0,2,j))     !emit V->e nu on lower
            call WMVtoln_mg2_had(1,q34(0,j),v,awen(0,0,1,j),zwen(0,0,1,j))  !emit V->mu+mu- on upper
            call WMVtoln_mg2_had(2,q12(0,j),v,awen(0,0,2,j),zwen(0,0,2,j))  !emit V->mu+mu- on lower
            call WMVtoln_mg2_had(1,q12(0,j),v,waen(0,0,2,j),wzen(0,0,2,j))  !emit V->mu+mu- on lower
            call WMVtoln_mg2_had(2,q34(0,j),v,waen(0,0,1,j),wzen(0,0,1,j))  !emit V->mu+mu- on upper
            call wmto4l_mg2_had(v,h,wwv)	

            call qqwmzqq(pbar,sign, nlo, 1,bos,uucs2,ddcs2,udcc2,udss2)

            print*, ""
            print*,"1 - uusc/uusc2 =   ",abs(1-uucs/uucs2)
            print*,"1 - ddsc/ddsc2 =   ",abs(1-ddcs/ddcs2)
            print*,"1 - ducc/ducc2 =   ",abs(1-udcc/udcc2)
            print*,"1 - duss/duss2 =   ",abs(1-udss/udss2)

            call anomal_formfactor(q12(0,j),q34(0,j),qe(0),qu(0))
            call vwtovwm_had(q12(0,j),q34(0,j),v,h,zwwv(0,0,j),awwv(0,0,j)) !q12 = V
            call vwtovwm_had(q34(0,j),q12(0,j),v,h,wzwv(0,0,j),wawv(0,0,j)) !q12 = W-
            call VVtomm_had(1,h,q12(0,j),v,aauu(0,0,2,j),azuu(0,0,2,j),zauu(0,0,2,j),zzuu(0,0,2,j))    	 !emit W on lower
            call VVtomm_had(2,h,q34(0,j),v,aauu(0,0,1,j),azuu(0,0,1,j),zauu(0,0,1,j),zzuu(0,0,1,j))    	 !emit W on upper
            call WWtomm_had(2,h,q34(0,j),v,wwuu(0,0,1,j))     !emit V->e nu on upper
            call WWtomm_had(2,h,q12(0,j),v,wwuu(0,0,2,j))     !emit V->e nu on lower
            call WMVtoln_had(1,q34(0,j),v,awen(0,0,1,j),zwen(0,0,1,j))  !emit V->mu+mu- on upper
            call WMVtoln_had(2,q12(0,j),v,awen(0,0,2,j),zwen(0,0,2,j))  !emit V->mu+mu- on lower
            call WMVtoln_had(1,q12(0,j),v,waen(0,0,2,j),wzen(0,0,2,j))  !emit V->mu+mu- on lower
            call WMVtoln_had(2,q34(0,j),v,waen(0,0,1,j),wzen(0,0,1,j))  !emit V->mu+mu- on upper
            call wmto4l_had(v,h,wwv)


         do mu = 0,3
          do nu = 0,3
           mg2_tens_2d(mu,nu) = dcmplx(0.,0.)
           mg2_tens_2d_comb(mu,nu,1) = dcmplx(0.,0.)
           mg2_tens_2d_comb(mu,nu,2) = dcmplx(0.,0.)
           mg2_tens_2d_comb(mu,nu,3) = dcmplx(0.,0.)
           mg2_tens_2d_comb(mu,nu,4) = dcmplx(0.,0.)
          enddo
         enddo
         do i = 0,5
           mg2_tens_1d(i) = dcmplx(0.,0.)
         enddo


         print*, ""
         call WMto4l_mg2_had(v,h,mg2_tens_1d)
         sumtest = 0
         do i = 0,5
           sumtest = sumtest + abs( wwv(i) / mg2_tens_1d(i) - 1)
!           print*, i, wwv(i), mg2_tens_1d(i), abs( wwv(i) / mg2_tens_1d(i) - 1)
           mg2_tens_1d(i) = dcmplx(0.,0.)
         enddo
         print*, "tensor W- to WZ:", sumtest


!         print*, ""
         call WMVtoln_mg2_had(1,q34(0,j),v,mg2_tens_2d_comb(0,0,1),mg2_tens_2d_comb(0,0,2))  !emit V->mu+mu- on upper
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( awen(mu,nu,1,j) / mg2_tens_2d_comb(mu,nu,1) - 1)
           ! print*, mu,nu, awen(mu,nu,1,j), mg2_tens_2d_comb(mu,nu,1), abs( awen(mu,nu,1,1) / mg2_tens_2d_comb(mu,nu,1) - 1)
           mg2_tens_2d_comb(mu,nu,1) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor W-Z to W- 1 had:", sumtest
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( zwen(mu,nu,1,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           ! print*, mu,nu, zwen(mu,nu,1,j), mg2_tens_2d_comb(mu,nu,2), abs( zwen(mu,nu,1,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           mg2_tens_2d_comb(mu,nu,2) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor W-Z to W- 1 had:", sumtest

!         print*, ""
         call WMVtoln_mg2_had(2,q12(0,j),v,mg2_tens_2d_comb(0,0,1),mg2_tens_2d_comb(0,0,2))  !emit V->mu+mu- on lower
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( awen(mu,nu,2,j) / mg2_tens_2d_comb(mu,nu,1) - 1)
           ! print*, mu,nu, awen(mu,nu,2,j), mg2_tens_2d_comb(mu,nu,1), abs( awen(mu,nu,2,1) / mg2_tens_2d_comb(mu,nu,1) - 1)
           mg2_tens_2d_comb(mu,nu,1) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor W-A to W- 2 had:", sumtest
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( zwen(mu,nu,2,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           ! print*, mu,nu, zwen(mu,nu,2,j), mg2_tens_2d_comb(mu,nu,2), abs( zwen(mu,nu,2,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           mg2_tens_2d_comb(mu,nu,2) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor W-Z to W- 2 had:", sumtest

!         print*, ""
         call WMVtoln_mg2_had(2,q34(0,j),v,mg2_tens_2d_comb(0,0,1),mg2_tens_2d_comb(0,0,2))  !emit V->mu+mu- on upper  
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( waen(mu,nu,1,j) / mg2_tens_2d_comb(mu,nu,1) - 1)
           ! print*, mu,nu, waen(mu,nu,1,j), mg2_tens_2d_comb(mu,nu,1), abs( waen(mu,nu,1,1) / mg2_tens_2d_comb(mu,nu,1) - 1)
           mg2_tens_2d_comb(mu,nu,1) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor W-A to W- 4 had:", sumtest
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( wzen(mu,nu,1,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           ! print*, mu,nu, wzen(mu,nu,1,j), mg2_tens_2d_comb(mu,nu,2), abs( wzen(mu,nu,1,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           mg2_tens_2d_comb(mu,nu,2) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor W-Z to W- 4 had:", sumtest

!         print*, ""
         call WMVtoln_mg2_had(1,q12(0,j),v,mg2_tens_2d_comb(0,0,1),mg2_tens_2d_comb(0,0,2))  !emit V->mu+mu- on lower
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( waen(mu,nu,2,j) / mg2_tens_2d_comb(mu,nu,1) - 1)
           ! print*, mu,nu, waen(mu,nu,2,j), mg2_tens_2d_comb(mu,nu,1), abs( waen(mu,nu,2,1) / mg2_tens_2d_comb(mu,nu,1) - 1)
           mg2_tens_2d_comb(mu,nu,1) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor W-A to W- 3 had:", sumtest
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( wzen(mu,nu,2,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           ! print*, mu,nu, wzen(mu,nu,2,j), mg2_tens_2d_comb(mu,nu,2), abs( wzen(mu,nu,2,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           mg2_tens_2d_comb(mu,nu,2) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor W-Z to W- 3 had:", sumtest


!         print*, ""
         call VVtomm_mg2_had(1,h,q12(0,j),v,mg2_tens_2d_comb(0,0,1),mg2_tens_2d_comb(0,0,2),
     &                                      mg2_tens_2d_comb(0,0,3),mg2_tens_2d_comb(0,0,4))    	 !emit W on lower
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( aauu(mu,nu,2,j) / mg2_tens_2d_comb(mu,nu,1) - 1)
           ! print*, mu,nu, aauu(mu,nu,2,j), mg2_tens_2d_comb(mu,nu,1), abs( aauu(mu,nu,2,j) / mg2_tens_2d_comb(mu,nu,1) - 1)
           mg2_tens_2d_comb(mu,nu,1) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor AA to VV 1:", sumtest
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( azuu(mu,nu,2,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           ! print*, mu,nu, azuu(mu,nu,2,j), mg2_tens_2d_comb(mu,nu,2), abs( azuu(mu,nu,2,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           mg2_tens_2d_comb(mu,nu,2) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor AZ to VV 1:", sumtest
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( zauu(mu,nu,2,j) / mg2_tens_2d_comb(mu,nu,3) - 1)
           ! print*, mu,nu, zauu(mu,nu,2,j), mg2_tens_2d_comb(mu,nu,3), abs( zauu(mu,nu,2,j) / mg2_tens_2d_comb(mu,nu,3) - 1)
           mg2_tens_2d_comb(mu,nu,3) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor ZA to VV 1:", sumtest
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( zzuu(mu,nu,2,j) / mg2_tens_2d_comb(mu,nu,4) - 1)
           ! print*, mu,nu, zzuu(mu,nu,2,j), mg2_tens_2d_comb(mu,nu,4), abs( zzuu(mu,nu,2,j) / mg2_tens_2d_comb(mu,nu,4) - 1)
           mg2_tens_2d_comb(mu,nu,4) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor ZZ to VV 1:", sumtest

!         print*, ""
         call VVtomm_mg2_had(2,h,q34(0,j),v,mg2_tens_2d_comb(0,0,1),mg2_tens_2d_comb(0,0,2),
     &                                      mg2_tens_2d_comb(0,0,3),mg2_tens_2d_comb(0,0,4))    	 !emit W on upper
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( aauu(mu,nu,1,j) / mg2_tens_2d_comb(mu,nu,1) - 1)
           ! print*, mu,nu, aauu(mu,nu,1,j), mg2_tens_2d_comb(mu,nu,1), abs( aauu(mu,nu,1,j) / mg2_tens_2d_comb(mu,nu,1) - 1)
           mg2_tens_2d_comb(mu,nu,1) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor AA to VV 2:", sumtest
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( azuu(mu,nu,1,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           ! print*, mu,nu, azuu(mu,nu,1,j), mg2_tens_2d_comb(mu,nu,2), abs( azuu(mu,nu,1,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           mg2_tens_2d_comb(mu,nu,2) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor AZ to VV 2:", sumtest
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( zauu(mu,nu,1,j) / mg2_tens_2d_comb(mu,nu,3) - 1)
           ! print*, mu,nu, zauu(mu,nu,1,j), mg2_tens_2d_comb(mu,nu,3), abs( zauu(mu,nu,1,j) / mg2_tens_2d_comb(mu,nu,3) - 1)
           mg2_tens_2d_comb(mu,nu,3) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor ZA to VV 2:", sumtest
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( zzuu(mu,nu,1,j) / mg2_tens_2d_comb(mu,nu,4) - 1)
           ! print*, mu,nu, zzuu(mu,nu,1,j), mg2_tens_2d_comb(mu,nu,4), abs( zzuu(mu,nu,1,j) / mg2_tens_2d_comb(mu,nu,4) - 1)
           mg2_tens_2d_comb(mu,nu,4) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor ZZ to VV 2:", sumtest


!         print*, ""
         call WWtomm_mg2_had(2,h,q34(0,j),v,mg2_tens_2d(0,0))     !emit V->e nu on upper
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( wwuu(mu,nu,1,j) / mg2_tens_2d(mu,nu) - 1)
           ! print*, mu,nu, wwuu(mu,nu,1,j), mg2_tens_2d(mu,nu), abs( wwuu(mu,nu,1,j) / mg2_tens_2d(mu,nu) - 1)
           mg2_tens_2d(mu,nu) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor WW to Z 1:", sumtest

!         print*, ""
         call WWtomm_mg2_had(2,h,q12(0,j),v,mg2_tens_2d(0,0))     !emit V->e nu on lower
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( wwuu(mu,nu,2,j) / mg2_tens_2d(mu,nu) - 1)
           ! print*, mu,nu, wwuu(mu,nu,2,j), mg2_tens_2d(mu,nu), abs( wwuu(mu,nu,2,j) / mg2_tens_2d(mu,nu) - 1)
           mg2_tens_2d(mu,nu) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor WW to Z 2:", sumtest


!         print*, ""
         call vwtovWM_mg2_had(q12(0,j),q34(0,j),v,h,mg2_tens_2d_comb(0,0,1),mg2_tens_2d_comb(0,0,2)) !q12 = V
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( zwwv(mu,nu,j) / mg2_tens_2d_comb(mu,nu,1) - 1)
           ! print*, mu,nu, zwwv(mu,nu,j), mg2_tens_2d_comb(mu,nu,1), abs( zwwv(mu,nu,j) / mg2_tens_2d_comb(mu,nu,1) - 1)
           mg2_tens_2d_comb(mu,nu,1) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor W-Z to W-V 1:", sumtest
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( awwv(mu,nu,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           ! print*, mu,nu, awwv(mu,nu,j), mg2_tens_2d_comb(mu,nu,2), abs( awwv(mu,nu,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           mg2_tens_2d_comb(mu,nu,2) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor W-A to W-V 1:", sumtest

!         print*, ""
         call vwtovWM_mg2_had(q34(0,j),q12(0,j),v,h,mg2_tens_2d_comb(0,0,1),mg2_tens_2d_comb(0,0,2)) !q12 = W-
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( wzwv(mu,nu,j) / mg2_tens_2d_comb(mu,nu,1) - 1)
           ! print*, mu,nu, wzwv(mu,nu,j), mg2_tens_2d_comb(mu,nu,1), abs( wzwv(mu,nu,j) / mg2_tens_2d_comb(mu,nu,1) - 1)
           mg2_tens_2d_comb(mu,nu,1) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor W-Z to W-V 2:", sumtest
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( wawv(mu,nu,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           ! print*, mu,nu, wawv(mu,nu,j), mg2_tens_2d_comb(mu,nu,2), abs( wawv(mu,nu,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           mg2_tens_2d_comb(mu,nu,2) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor W-A to W-V 2:", sumtest


         print*, ""
         print*, ""
         read(*,*)


         end select

#endif


      end


c ------------------------------------------------------------	

      subroutine wpto4l_mg2_had(p,h,wp4l)
	
      IMPLICIT NONE

      REAL*8 P(0:3,2:5)
      complex*16 wp4l(0:5)
      integer lsign(2:5,4),h,ll,mu   

	! possible lepton helicities:
	data (lsign(ll,1),ll=2,5) / -1, -1, -1, -1/  !h=1
	data (lsign(ll,2),ll=2,5) / -1, -1,  1,  1/  !h=2

	call wpto4ln_mg2_had(P,lsign(2,h),wp4l)
	end
		
c ------------------------------------------------------------	
c ------------------------------------------------------------	
   
      subroutine wpto4LN_mg2_had(P,lsign,wp)

C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : w+ -> e+ ve mu+ mu-  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=   7) 
      INTEGER    NWAVEFUNCS  
      PARAMETER (NWAVEFUNCS=  15) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,2:5),q(0:3)
      INTEGER lsign(2:5)
      complex*16 wp(0:5)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J,mu
      COMPLEX*16 AMP(NGRAPHS)
      COMPLEX*16 W(6,NWAVEFUNCS)
C  
c  lepton spinors and V polarization vectors
#include "tensorw.inc"
C  
C GLOBAL VARIABLES
c
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/process.inc"
   
      logical lswitch
      parameter (lswitch=.true.)
C                                                                          
C ----------
C BEGIN CODE
C ----------


      if (lswitch) then
	do i = 1,6
	   w(i,2) = wl(i)
	   w(i,3) = wn(i)
	   w(i,4) = lp(i)
	   w(i,5) = lm(i)
	   
	   w(i,6) = we(i)
	   w(i,10) = au(i)
	   w(i,11) = zu(i)
	enddo
      else
      	CALL IXXXXX(P(0,2   ),ZERO ,-lsign(2),-1,W(1,2   ))       
      	CALL OXXXXX(P(0,3   ),ZERO ,+lsign(3),+1,W(1,3   ))       
      	CALL IXXXXX(P(0,4   ),ZERO ,-lsign(4),-1,W(1,4   ))       
      	CALL OXXXXX(P(0,5   ),ZERO ,+lsign(5),+1,W(1,5   )) 

      	CALL JIOXXX(W(1,2   ),W(1,3   ),GWF ,WMASS   ,WWIDTH  ,W(1,6   ))    
      endif
             
      do mu = 0,3
         q(mu) = p(mu,2)+p(mu,3)+p(mu,4)+p(mu,5)
      enddo

      do mu = 0,3
      CALL VCARTX(Q(0),wmass ,wwidth, MU,-1,W(1,1   ))


      if (procid.eq.wphadzjj) then
      

      CALL JIOXXX(W(1,2   ),W(1,3   ),GWF ,WMASS   ,WWIDTH  ,W(1,6   ))    
      CALL JVVXXX(W(1,1   ),W(1,6   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     7   ))                                                          
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,7   ),GAL ,AMP(1   ))            
      CALL JVVXXX(W(1,1   ),W(1,6   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     8   ))                                                          
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,8   ),GZL ,AMP(2   ))            
      CALL FVIXXX(W(1,4   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,9   ))    
      CALL IOVXXX(W(1,9   ),W(1,5   ),W(1,6   ),GWF ,AMP(3   ))            
      CALL FVOXXX(W(1,3   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,10  ))    
      CALL JIOXXX(W(1,2   ),W(1,10  ),GAD ,ZERO    ,AWIDTH  ,W(1,11  ))    
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,11  ),GAL ,AMP(4   ))            
      CALL JIOXXX(W(1,2   ),W(1,10  ),GZD ,ZMASS   ,ZWIDTH  ,W(1,12  ))    
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,12  ),GZL ,AMP(5   ))            
      CALL FVIXXX(W(1,2   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,13  ))    
      CALL JIOXXX(W(1,13  ),W(1,3   ),GAU ,ZERO    ,AWIDTH  ,W(1,14  ))    
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,14  ),GAL ,AMP(6   ))            
      CALL JIOXXX(W(1,13  ),W(1,3   ),GZU ,ZMASS   ,ZWIDTH  ,W(1,15  ))    
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,15  ),GZL ,AMP(7   ))            
      wp(   mu)  = +AMP(   1)+AMP(   2)+AMP(   3)+AMP(   4)+AMP(   5)
     &             +AMP(   6)+AMP(   7)
      
      
      elseif (procid.eq.wpzhadjj .and. ztouu) then


      CALL JIOXXX(W(1,2   ),W(1,3   ),GWF ,WMASS   ,WWIDTH  ,W(1,6   ))    
      CALL FVIXXX(W(1,4   ),W(1,6   ),GWF ,ZERO    ,ZERO    ,W(1,7   ))    
      CALL IOVXXX(W(1,7   ),W(1,5   ),W(1,1   ),GWF ,AMP(1   ))            
      CALL JVVXXX(W(1,1   ),W(1,6   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     8   ))                                                          
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,8   ),GAU ,AMP(2   ))            
      CALL JVVXXX(W(1,1   ),W(1,6   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     9   ))                                                          
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,9   ),GZU ,AMP(3   ))            
      CALL FVOXXX(W(1,3   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,10  ))    
      CALL JIOXXX(W(1,2   ),W(1,10  ),GAL ,ZERO    ,AWIDTH  ,W(1,11  ))    
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,11  ),GAU ,AMP(4   ))            
      CALL JIOXXX(W(1,2   ),W(1,10  ),GZL ,ZMASS   ,ZWIDTH  ,W(1,12  ))    
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,12  ),GZU ,AMP(5   ))            
      CALL FVIXXX(W(1,2   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,13  ))    
      CALL JIOXXX(W(1,13  ),W(1,3   ),GZN ,ZMASS   ,ZWIDTH  ,W(1,14  ))    
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,14  ),GZU ,AMP(6   ))            
      wp(   mu)  = +AMP(   1)+AMP(   2)+AMP(   3)+AMP(   4)+AMP(   5)
     &             +AMP(   6)


      elseif (procid.eq.wpzhadjj .and. ztodd) then


      CALL JIOXXX(W(1,2   ),W(1,3   ),GWF ,WMASS   ,WWIDTH  ,W(1,6   ))    
      CALL JVVXXX(W(1,1   ),W(1,6   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     7   ))                                                          
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,7   ),GAD ,AMP(1   ))            
      CALL JVVXXX(W(1,1   ),W(1,6   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     8   ))                                                          
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,8   ),GZD ,AMP(2   ))            
      CALL FVIXXX(W(1,4   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,9   ))    
      CALL IOVXXX(W(1,9   ),W(1,5   ),W(1,6   ),GWF ,AMP(3   ))            
      CALL FVOXXX(W(1,3   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,10  ))    
      CALL JIOXXX(W(1,2   ),W(1,10  ),GAL ,ZERO    ,AWIDTH  ,W(1,11  ))    
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,11  ),GAD ,AMP(4   ))            
      CALL JIOXXX(W(1,2   ),W(1,10  ),GZL ,ZMASS   ,ZWIDTH  ,W(1,12  ))    
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,12  ),GZD ,AMP(5   ))            
      CALL FVIXXX(W(1,2   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,13  ))    
      CALL JIOXXX(W(1,13  ),W(1,3   ),GZN ,ZMASS   ,ZWIDTH  ,W(1,14  ))    
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,14  ),GZD ,AMP(6   ))            
      wp(   mu)  = +AMP(   1)+AMP(   2)+AMP(   3)+AMP(   4)+AMP(   5)
     &             +AMP(   6)


      endif


      enddo
             
      wp(4) = dcmplx(q(0),q(3))
      wp(5) = dcmplx(q(1),q(2))

      return
      END


c ------------------------------------------------------------	

      subroutine wmto4l_mg2_had(p,h,wm4l)
	
      IMPLICIT NONE

      REAL*8 P(0:3,2:5)
      complex*16 wm4l(0:5)
      integer lsign(2:5,4),h,ll,mu   

	! possible lepton helicities:
	data (lsign(ll,1),ll=2,5) / -1, -1, -1, -1/  !h=1
	data (lsign(ll,2),ll=2,5) / -1, -1,  1,  1/  !h=2

	call wmto4ln_mg2_had(P,lsign(2,h),wm4l)
	end
c ------------------------------------------------------------	
       
      subroutine WMTO4LN_mg2_had(P,lsign,wm)
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : w- -> e- ve~ mu+ mu-  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   7) 
      INTEGER    NWAVEFUNCS  
      PARAMETER (NWAVEFUNCS=  15) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,2:5),q(0:3)
      INTEGER lsign(2:5)
      complex*16 wm(0:5)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J,mu
      COMPLEX*16 ZTEMP
      COMPLEX*16 AMP(NGRAPHS)
      COMPLEX*16 W(6,NWAVEFUNCS)
C   
c  lepton spinors and V polarization vectors
#include "tensorw.inc"
c 
C GLOBAL VARIABLES
C  
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/process.inc"
   
      logical lswitch
      parameter (lswitch=.true.)
C  
c ----------
C BEGIN CODE
C ----------

	do i = 1,6
	   w(i,2) = wl(i)
	   w(i,3) = wn(i)
	   w(i,4) = lp(i)
	   w(i,5) = lm(i)
	   
	   w(i,6) = we(i)	   
 	   w(i,10) = au(i)	   
	   w(i,11) = zu(i)	   
	 enddo  

             
      do mu = 0,3
         q(mu) = p(mu,2)+p(mu,3)+p(mu,4)+p(mu,5)
      enddo


      do mu = 0,3
      CALL VCARTX(Q(0),wmass ,wwidth, MU,-1,W(1,1   ))


      if (procid.eq.wmhadzjj) then
      

      CALL JIOXXX(W(1,3   ),W(1,2   ),GWF ,WMASS   ,WWIDTH  ,W(1,6   ))    
      CALL FVIXXX(W(1,4   ),W(1,6   ),GWF ,ZERO    ,ZERO    ,W(1,7   ))    
      CALL IOVXXX(W(1,7   ),W(1,5   ),W(1,1   ),GWF ,AMP(1   ))            
      CALL JVVXXX(W(1,6   ),W(1,1   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     8   ))                                                          
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,8   ),GAL ,AMP(2   ))            
      CALL JVVXXX(W(1,6   ),W(1,1   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     9   ))                                                          
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,9   ),GZL ,AMP(3   ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,10  ))    
      CALL JIOXXX(W(1,10  ),W(1,2   ),GAD ,ZERO    ,AWIDTH  ,W(1,11  ))    
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,11  ),GAL ,AMP(4   ))            
      CALL JIOXXX(W(1,10  ),W(1,2   ),GZD ,ZMASS   ,ZWIDTH  ,W(1,12  ))    
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,12  ),GZL ,AMP(5   ))            
      CALL FVOXXX(W(1,2   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,13  ))    
      CALL JIOXXX(W(1,3   ),W(1,13  ),GAU ,ZERO    ,AWIDTH  ,W(1,14  ))    
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,14  ),GAL ,AMP(6   ))            
      CALL JIOXXX(W(1,3   ),W(1,13  ),GZU ,ZMASS   ,ZWIDTH  ,W(1,15  ))    
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,15  ),GZL ,AMP(7   ))            
      wm(   mu)  = -AMP(   1)-AMP(   2)-AMP(   3)-AMP(   4)-AMP(   5)
     &             -AMP(   6)-AMP(   7)
      
      
      elseif (procid.eq.wmzhadjj .and. ztouu) then


      CALL JIOXXX(W(1,3   ),W(1,2   ),GWF ,WMASS   ,WWIDTH  ,W(1,6   ))    
      CALL JVVXXX(W(1,6   ),W(1,1   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     7   ))                                                          
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,7   ),GAU ,AMP(1   ))            
      CALL JVVXXX(W(1,6   ),W(1,1   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     8   ))                                                          
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,8   ),GZU ,AMP(2   ))            
      CALL FVIXXX(W(1,4   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,9   ))    
      CALL IOVXXX(W(1,9   ),W(1,5   ),W(1,6   ),GWF ,AMP(3   ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,10  ))    
      CALL JIOXXX(W(1,10  ),W(1,2   ),GAL ,ZERO    ,AWIDTH  ,W(1,11  ))    
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,11  ),GAU ,AMP(4   ))            
      CALL JIOXXX(W(1,10  ),W(1,2   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,12  ))    
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,12  ),GZU ,AMP(5   ))            
      CALL FVOXXX(W(1,2   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,13  ))    
      CALL JIOXXX(W(1,3   ),W(1,13  ),GZN ,ZMASS   ,ZWIDTH  ,W(1,14  ))    
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,14  ),GZU ,AMP(6   ))            
      wm(   mu)  = -AMP(   1)-AMP(   2)-AMP(   3)-AMP(   4)-AMP(   5)
     &             -AMP(   6)


      elseif (procid.eq.wmzhadjj .and. ztodd) then


      CALL JIOXXX(W(1,3   ),W(1,2   ),GWF ,WMASS   ,WWIDTH  ,W(1,6   ))    
      CALL FVIXXX(W(1,4   ),W(1,6   ),GWF ,ZERO    ,ZERO    ,W(1,7   ))    
      CALL IOVXXX(W(1,7   ),W(1,5   ),W(1,1   ),GWF ,AMP(1   ))            
      CALL JVVXXX(W(1,6   ),W(1,1   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     8   ))                                                          
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,8   ),GAD ,AMP(2   ))            
      CALL JVVXXX(W(1,6   ),W(1,1   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     9   ))                                                          
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,9   ),GZD ,AMP(3   ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,10  ))    
      CALL JIOXXX(W(1,10  ),W(1,2   ),GAL ,ZERO    ,AWIDTH  ,W(1,11  ))    
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,11  ),GAD ,AMP(4   ))            
      CALL JIOXXX(W(1,10  ),W(1,2   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,12  ))    
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,12  ),GZD ,AMP(5   ))            
      CALL FVOXXX(W(1,2   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,13  ))    
      CALL JIOXXX(W(1,3   ),W(1,13  ),GZN ,ZMASS   ,ZWIDTH  ,W(1,14  ))    
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,14  ),GZD ,AMP(6   ))            
      wm(   mu)  = -AMP(   1)-AMP(   2)-AMP(   3)-AMP(   4)-AMP(   5)
     &             -AMP(   6)


      endif


      enddo
             
      wm(4) = dcmplx(q(0),q(3))
      wm(5) = dcmplx(q(1),q(2))

      END
		

c ------------------------------------------------------------	

      subroutine WMVtoln_mg2_had(id,Q,v,wmatens,wmztens)
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
#include "tensorw.inc"

      logical ldebug
      parameter (ldebug=.false.)
C  
C GLOBAL VARIABLES
C  
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/process.inc"
C                                                                          
C ----------
C BEGIN CODE
C ----------
      if (id.eq.1) then
         q2(0) = dreal(we(5))
         q2(1) = dreal(we(6))
         q2(2) = dimag(we(6))
         q2(3) = dimag(we(5))
         do mu = 0,3
            q1(mu) = q(mu)
            q2(mu) = q2(mu)-q(mu)
         enddo
      elseif (id.eq.2) then
         q1(0) = dreal(we(5))
         q1(1) = dreal(we(6))
         q1(2) = dimag(we(6))
         q1(3) = dimag(we(5))
         do mu = 0,3
            q2(mu) = q(mu)
            q1(mu) = q1(mu)-q(mu)
         enddo
      endif

      do i=1,6                  ! read in lepton states from common block
         w(i,3) = wl(i)
         w(i,4) = wn(i)
	 w(i,7) = we(i)
      enddo

      q2sq = q2(0)**2-q2(1)**2-q2(2)**2-q2(3)**2
      denz = q2sq/dcmplx(q2sq-zmass**2,zmass*zwidth)
      DO MU = 0,3
         DO NU = 0,3
            CALL VCARTX(Q1(0),WMASS ,Wwidth, MU,-1,W(1,1))        
            CALL VCARTX(Q2(0),ZERO  ,ZERO  , NU,-1,W(1,2))
            CALL FVOXXX(W(1,3),W(1,2),GZL ,ZERO,ZERO,W(1,5))    
            CALL FVOXXX(W(1,3),W(1,1),GWF ,ZERO,ZERO,W(1,6))    
            CALL IOVXXX(W(1,4),W(1,5),W(1,1),GWF ,AMP(1))   
            CALL IOVXXX(W(1,4),W(1,6),W(1,2),GZN ,AMP(2))            
            CALL VVVXXX (W(1,7),W(1,1),W(1,2),GWWZ,AMP(3))
            wmztens(mu,nu) = (-AMP(   1)-AMP(   2)-AMP(   3))*denz
            wmatens(mu,nu) = -(amp(1)*gal(1)/gzl(1)+amp(3)*gwwa/gwwz)
            

            if (procid.eq.wmhadzjj) then


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
      wmatens(mu,nu)  = -AMP(   1)-AMP(   2)-AMP(   3)

            endif


         enddo
      enddo

      end

c----------------------------------------------------------------

      subroutine WPVtoln_mg2_had(id,Q,v,wpatens,wpztens)
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
#include "tensorw.inc"

      logical ldebug
      parameter (ldebug=.false.)
C  
C GLOBAL VARIABLES
C  
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/process.inc"
C                                                                          
C ----------
C BEGIN CODE
C ----------
      
      if (id.eq.1) then
         q2(0) = dreal(we(5))
         q2(1) = dreal(we(6))
         q2(2) = dimag(we(6))
         q2(3) = dimag(we(5))
         do mu = 0,3
            q1(mu) = q(mu)
            q2(mu) = q2(mu)-q(mu)
         enddo
      elseif (id.eq.2) then
         q1(0) = dreal(we(5))
         q1(1) = dreal(we(6))
         q1(2) = dimag(we(6))
         q1(3) = dimag(we(5))
         do mu = 0,3
            q2(mu) = q(mu)
            q1(mu) = q1(mu)-q(mu)
         enddo
      endif
c
      do i=1,6                  ! read in lepton states from common block
         w(i,3) = wl(i)
         w(i,4) = wn(i)
	 w(i,7) = we(i)
      enddo

      q2sq = q2(0)**2-q2(1)**2-q2(2)**2-q2(3)**2
      denz = q2sq/dcmplx(q2sq-zmass**2,zmass*zwidth)
      DO MU = 0,3
         DO NU = 0,3
            CALL VCARTX(Q1(0),WMASS ,Wwidth, MU,-1,W(1,1))        
            CALL VCARTX(Q2(0),ZERO  ,ZERO  , NU,-1,W(1,2))
            CALL FVIXXX(W(1,3),W(1,2),GZL ,ZERO,ZERO,W(1,5))    
            CALL IOVXXX(W(1,5),W(1,4),W(1,1),GWF ,AMP(1))            
            CALL FVIXXX(W(1,3),W(1,1),GWF ,ZERO,ZERO,W(1,6))    
            CALL IOVXXX(W(1,6),W(1,4),W(1,2),GZN ,AMP(2))	                
            CALL VVVXXX(W(1,1),W(1,7),W(1,2),GWWZ ,AMP(3))            
            wpztens(mu,nu) = (-AMP(   1)-AMP(   2)-AMP(   3))*denz
            wpatens(mu,nu) = -(amp(1)*gal(1)/gzl(1)+amp(3)*gwwa/gwwz)


            if (procid.eq.wphadzjj) then


            CALL VCARTX(Q2(0),ZMASS  ,ZWIDTH  , NU,-1,W(1,2))
      CALL FVIXXX(W(1,3   ),W(1,2   ),GZD ,ZERO    ,ZERO    ,W(1,5   ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,1   ),GWF ,AMP(1   ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,6   ))    
      CALL IOVXXX(W(1,6   ),W(1,4   ),W(1,2   ),GZU ,AMP(2   ))            
      CALL JVVXXX(W(1,2   ),W(1,1   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     7   ))                                                          
      CALL IOVXXX(W(1,3   ),W(1,4   ),W(1,7   ),GWF ,AMP(3   ))            
      wpztens(mu,nu) = -(+AMP(   1)+AMP(   2)+AMP(   3))


            CALL VCARTX(Q2(0),ZERO  ,ZERO  , NU,-1,W(1,2))
      CALL FVIXXX(W(1,3   ),W(1,2   ),GAD ,ZERO    ,ZERO    ,W(1,5   ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,1   ),GWF ,AMP(1   ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,6   ))    
      CALL IOVXXX(W(1,6   ),W(1,4   ),W(1,2   ),GAU ,AMP(2   ))            
      CALL JVVXXX(W(1,2   ),W(1,1   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     7   ))                                                          
      CALL IOVXXX(W(1,3   ),W(1,4   ),W(1,7   ),GWF ,AMP(3   ))            
      wpatens(mu,nu) = -(+AMP(   1)+AMP(   2)+AMP(   3))


            endif


         enddo
      enddo

      end


c ------------------------------------------------------------	

      subroutine vvtomm_mg2_had(idb,h,Q,v,aa,az,za,zz)
      
c 	tensor for VV->mu+mu- with V=A,Z      
	
      IMPLICIT NONE

        integer idb,idbb,h
        real*8 q(0:3), v(0:3,4)
        complex*16 aa(0:3,0:3)
        complex*16 az(0:3,0:3)
        complex*16 za(0:3,0:3)
        complex*16 zz(0:3,0:3)
	integer lsign(2,2),ll,mu,nu
      	
	data (lsign(ll,1),ll=1,2) / -1, -1/
	data (lsign(ll,2),ll=1,2) / +1, +1/
	
	idbb = -idb+3

	call AAtoMM1_mg2_had(idb ,Q,v,lsign(1,h),aa(0,0))
	call AZtoMM1_mg2_had(idb ,Q,v,lsign(1,h),az(0,0))
	call AZtoMM1_mg2_had(idbb,Q,v,lsign(1,h),za(0,0)) 
	call ZZtoMM1_mg2_had(idb ,Q,v,lsign(1,h),zz(0,0))
				
	end
		
c ------------------------------------------------------------	
       
      subroutine AAtoMM1_mg2_had(idb,Q,v,lsign,aatens)
C  Input:   Q(*)        one of the two incoming boson momenta 
c           idb = 1      Q = q1 = A(1) momentum
c           idb = 2      Q = q2 = A(2) momentum
c  OUPTPUT: aatens(mu,nu)  lepton tensor for a a --> mu+mu-  
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : a a -> mu+ mu-  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS 
      PARAMETER (NGRAPHS=   2) 
      INTEGER    NWAVEFUNCS   
      PARAMETER (NWAVEFUNCS=   6) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      integer idb
      real*8 q(0:3), v(0:3,4), p(0:3,3:4)
      complex*16 aatens(0:3,0:3)
      integer lsign(3:4)
C    
C LOCAL VARIABLES 
C  
      INTEGER I,J,mu,nu
      real*8 Q1(0:3), Q2(0:3)
      COMPLEX*16 AMP(NGRAPHS)
      COMPLEX*16 W(6,NWAVEFUNCS)
C  
c  lepton spinors 
#include "tensorw.inc"
C  
C GLOBAL VARIABLES
C  
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/process.inc"
   
      logical lswitch
      parameter (lswitch=.true.)
C                                                                          
C ----------
C BEGIN CODE
C ----------
      
      if (idb.eq.1) then
         do mu = 0,3
            q1(mu) = q(mu)
	    q2(mu) = qu(mu)-q(mu)
	    p(mu,3) = v(mu,3)
	    p(mu,4) = v(mu,4)
         enddo
      elseif (idb.eq.2) then
         do mu = 0,3
            q2(mu) = q(mu)
            q1(mu) = qu(mu)-q(mu)
	    p(mu,3) = v(mu,3)
	    p(mu,4) = v(mu,4)
         enddo
      endif
	 
      if (lswitch) then
	do i = 1,6
	   w(i,3) = lp(i)
	   w(i,4) = lm(i)
	enddo
      else
       	CALL IXXXXX(P(0,3   ),ZERO ,-lsign(3),-1,W(1,3   ))        
       	CALL OXXXXX(P(0,4   ),ZERO ,+lsign(4),+1,W(1,4   )) 
      endif    
	     
      DO MU = 0,3
         DO NU = 0,3
      
      CALL VCARTX(Q1(0),zero ,zero, MU,-1,W(1,1	  ))	    
      CALL VCARTX(Q2(0),zero ,zero, NU,-1,W(1,2	  )) 
 
      CALL FVIXXX(W(1,3   ),W(1,2   ),GAL ,ZERO    ,ZERO    ,W(1,5   ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,1   ),GAL ,AMP(1   ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,6   ))    
      CALL IOVXXX(W(1,6   ),W(1,4   ),W(1,2   ),GAL ,AMP(2   ))            
      
      aatens(mu,nu) = +AMP(   1)+AMP(   2)


      if ( (procid.eq.wpzhadjj .or. procid.eq.wmzhadjj) .and. ztouu ) then
 
      CALL FVIXXX(W(1,3   ),W(1,2   ),GAU ,ZERO    ,ZERO    ,W(1,5   ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,1   ),GAU ,AMP(1   ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GAU ,ZERO    ,ZERO    ,W(1,6   ))    
      CALL IOVXXX(W(1,6   ),W(1,4   ),W(1,2   ),GAU ,AMP(2   ))            
      aatens(mu,nu) = +AMP(   1)+AMP(   2)

 
      elseif ( (procid.eq.wpzhadjj .or. procid.eq.wmzhadjj) .and. ztodd ) then

      CALL FVIXXX(W(1,3   ),W(1,2   ),GAD ,ZERO    ,ZERO    ,W(1,5   ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,1   ),GAD ,AMP(1   ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GAD ,ZERO    ,ZERO    ,W(1,6   ))    
      CALL IOVXXX(W(1,6   ),W(1,4   ),W(1,2   ),GAD ,AMP(2   ))            
      aatens(mu,nu) = +AMP(   1)+AMP(   2)

      
      endif
        
	 enddo !nu
      enddo !mu

      end

c ------------------------------------------------------------	
       
      subroutine AZtoMM1_mg2_had(idb,Q,v,lsign,aztens)
C  Input:   Q(*)        one of the two incoming boson momenta 
c           idb = 1      Q = q1 = A(1) momentum
c           idb = 2      Q = q2 = Z(2) momentum
c  OUPTPUT: aztens(mu,nu)  lepton tensor for a z --> l+l-  
C  
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : a z -> mu+ mu-  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS 
      PARAMETER (NGRAPHS=   2) 
      INTEGER    NWAVEFUNCS  
      PARAMETER (NWAVEFUNCS=   6) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      integer idb
      real*8 q(0:3), v(0:3,4), p(0:3,3:4)
      complex*16 aztens(0:3,0:3)
      integer lsign(3:4)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J,mu,nu
      real*8 Q1(0:3), Q2(0:3)
      COMPLEX*16 AMP(NGRAPHS)
      COMPLEX*16 W(6,NWAVEFUNCS)
C  
c  lepton spinors 
#include "tensorw.inc"
C  
C GLOBAL VARIABLES
C  
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/process.inc"
   
      logical lswitch
      parameter (lswitch=.true.)
C  
C ----------
C BEGIN CODE
C ----------
      
      if (idb.eq.1) then
         do mu = 0,3
            q1(mu) = q(mu)
	    q2(mu) = qu(mu)-q(mu)
	    p(mu,3) = v(mu,3)
	    p(mu,4) = v(mu,4)
         enddo
      elseif (idb.eq.2) then
         do mu = 0,3
            q2(mu) = q(mu)
            q1(mu) = qu(mu)-q(mu)
	    p(mu,3) = v(mu,3)
	    p(mu,4) = v(mu,4)
         enddo
      endif
	 

      if (lswitch) then
	do i = 1,6
	   w(i,3) = lp(i)
	   w(i,4) = lm(i)
	enddo
      else
       	CALL IXXXXX(P(0,3   ),ZERO ,-lsign(3),-1,W(1,3   ))        
       	CALL OXXXXX(P(0,4   ),ZERO ,+lsign(4),+1,W(1,4   )) 
      endif 
             
      DO MU = 0,3
         DO NU = 0,3
      
      CALL VCARTX(Q1(0),zero  ,zero, MU,-1,W(1,1   ))        
      CALL VCARTX(Q2(0),zMASS ,zwidth, NU,-1,W(1,2   )) 
 
      CALL FVIXXX(W(1,3   ),W(1,2   ),GZL ,ZERO    ,ZERO    ,W(1,5   ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,1   ),GAL ,AMP(1   ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,6   ))    
      CALL IOVXXX(W(1,6   ),W(1,4   ),W(1,2   ),GZL ,AMP(2   ))            
      
      aztens(mu,nu) = +AMP(   1)+AMP(   2)


      if ( (procid.eq.wpzhadjj .or. procid.eq.wmzhadjj) .and. ztouu ) then
 
      CALL FVIXXX(W(1,3   ),W(1,2   ),GZU ,ZERO    ,ZERO    ,W(1,5   ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,1   ),GAU ,AMP(1   ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GAU ,ZERO    ,ZERO    ,W(1,6   ))    
      CALL IOVXXX(W(1,6   ),W(1,4   ),W(1,2   ),GZU ,AMP(2   ))            
      aztens(mu,nu) = +AMP(   1)+AMP(   2)

 
      elseif ( (procid.eq.wpzhadjj .or. procid.eq.wmzhadjj) .and. ztodd ) then

      CALL FVIXXX(W(1,3   ),W(1,2   ),GZD ,ZERO    ,ZERO    ,W(1,5   ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,1   ),GAD ,AMP(1   ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GAD ,ZERO    ,ZERO    ,W(1,6   ))    
      CALL IOVXXX(W(1,6   ),W(1,4   ),W(1,2   ),GZD ,AMP(2   ))            
      aztens(mu,nu) = +AMP(   1)+AMP(   2)

      
      endif


	 enddo !nu
      enddo !mu

      end

c ------------------------------------------------------------	
       
      subroutine ZZtoMM1_mg2_had(idb,Q,v,lsign,zztens)
C  Input:   Q(*)        one of the two incoming boson momenta 
c           idb = 1      Q = q1 = Z(1) momentum
c           idb = 2      Q = q2 = Z(2) momentum
c  OUPTPUT: zztens(mu,nu)  lepton tensor for z z --> l+l-  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : z z -> mu+ mu-  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS 
      PARAMETER (NGRAPHS=   2) 
      INTEGER    NWAVEFUNCS   
      PARAMETER (NWAVEFUNCS=   6) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      integer idb
      real*8 q(0:3), v(0:3,4), p(0:3,3:4)
      complex*16 zztens(0:3,0:3)
      integer lsign(3:4)
C  
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J,mu,nu
      real*8 Q1(0:3), Q2(0:3)
      COMPLEX*16 AMP(NGRAPHS)
      COMPLEX*16 W(6,NWAVEFUNCS)
C  
c  lepton spinors 
#include "tensorw.inc"
C  
C GLOBAL VARIABLES
C  
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/process.inc"
   
      logical lswitch
      parameter (lswitch=.true.)
C  
C ----------
C BEGIN CODE
C ----------
      
      if (idb.eq.1) then
         do mu = 0,3
            q1(mu) = q(mu)
	    q2(mu) = qu(mu)-q(mu)
	    p(mu,3) = v(mu,3)
	    p(mu,4) = v(mu,4)
         enddo
      elseif (idb.eq.2) then
         do mu = 0,3
            q2(mu) = q(mu)
            q1(mu) = qu(mu)-q(mu)
	    p(mu,3) = v(mu,3)
	    p(mu,4) = v(mu,4)
         enddo
      endif 
 
      if (lswitch) then
	do i = 1,6
	   w(i,3) = lp(i)
	   w(i,4) = lm(i)
	enddo
      else
       	CALL IXXXXX(P(0,3   ),ZERO ,-lsign(3),-1,W(1,3   ))        
       	CALL OXXXXX(P(0,4   ),ZERO ,+lsign(4),+1,W(1,4   )) 
      endif 
             
      DO MU = 0,3
         DO NU = 0,3
      
      CALL VCARTX(Q1(0),zMASS ,zwidth, MU,-1,W(1,1   ))        
      CALL VCARTX(Q2(0),zMASS ,zwidth, NU,-1,W(1,2   )) 
 
      CALL FVIXXX(W(1,3   ),W(1,2   ),GZL ,ZERO    ,ZERO    ,W(1,5   ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,1   ),GZL ,AMP(1   ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GZL ,ZERO    ,ZERO    ,W(1,6   ))    
      CALL IOVXXX(W(1,6   ),W(1,4   ),W(1,2   ),GZL ,AMP(2   ))            
      
      zztens(mu,nu) = +AMP(   1)+AMP(   2)


      if ( (procid.eq.wpzhadjj .or. procid.eq.wmzhadjj) .and. ztouu ) then
 
      CALL FVIXXX(W(1,3   ),W(1,2   ),GZU ,ZERO    ,ZERO    ,W(1,5   ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,1   ),GZU ,AMP(1   ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GZU ,ZERO    ,ZERO    ,W(1,6   ))    
      CALL IOVXXX(W(1,6   ),W(1,4   ),W(1,2   ),GZU ,AMP(2   ))            
      zztens(mu,nu) = +AMP(   1)+AMP(   2)

 
      elseif ( (procid.eq.wpzhadjj .or. procid.eq.wmzhadjj) .and. ztodd ) then

      CALL FVIXXX(W(1,3   ),W(1,2   ),GZD ,ZERO    ,ZERO    ,W(1,5   ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,1   ),GZD ,AMP(1   ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GZD ,ZERO    ,ZERO    ,W(1,6   ))    
      CALL IOVXXX(W(1,6   ),W(1,4   ),W(1,2   ),GZD ,AMP(2   ))            
      zztens(mu,nu) = +AMP(   1)+AMP(   2)

      
      endif

        
	 enddo !nu
      enddo !mu

      end

  
c ------------------------------------------------------------


      subroutine wwtomm_mg2_had(idb,h,Q,v,wll)

c 	tensor for WW->mu+mu-       
	
      IMPLICIT NONE

        integer idb,h
        real*8 q(0:3), v(0:3,4)
        complex*16 wll(0:3,0:3)
	integer lsign(2,2),i,ll,mu,nu     

	data (lsign(ll,1),ll=1,2) / -1, -1/
	data (lsign(ll,2),ll=1,2) / +1, +1/
		
	call WWtoMM1_mg2_had(idb,Q,v,lsign(1,h),wll(0,0))
		
	end
		
c ------------------------------------------------------------	
       
      subroutine WWtoMM1_mg2_had(idb,Q,v,lsign,wll)
C  Input:   Q(*)        one of the two incoming boson momenta 
c           idb = 1      Q = q1 = W+ momentum
c           idb = 2      Q = q2 = W- momentum
c  OUPTPUT: wll(mu,nu)  lepton tensor for w w--> l+l-  
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : w+ w- -> l+ l-  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=   3) 
      INTEGER    NWAVEFUNCS  
      PARAMETER (NWAVEFUNCS=   9) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      integer idb
      real*8 q(0:3), v(0:3,4), p(0:3,3:4)
      complex*16 wll(0:3,0:3)
      integer lsign(3:4)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J, MU, NU
      real*8 Q1(0:3), Q2(0:3)
      COMPLEX*16 AMP(NGRAPHS)
      COMPLEX*16 W(6,NWAVEFUNCS)
c  lepton spinors 
#include "tensorw.inc"
C  
C GLOBAL VARIABLES
C  
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/process.inc"
   
      logical lswitch
      parameter (lswitch=.true.)
C  
C ----------
C BEGIN CODE
C ----------
      
      if (idb.eq.1) then
         do mu = 0,3
            q1(mu)  = q(mu)
	    q2(mu)  = qu(mu)-q(mu)
	    p(mu,3) = v(mu,3)
	    p(mu,4) = v(mu,4)
         enddo
      elseif (idb.eq.2) then
         do mu = 0,3
            q2(mu)  = q(mu)
            q1(mu)  = qu(mu)-q(mu)
	    p(mu,3) = v(mu,3)
	    p(mu,4) = v(mu,4)
         enddo
      endif
	  
      if (lswitch) then
      	do i = 1,6
         w(i,3) = lp(i)
         w(i,4) = lm(i)
   
       	 w(i,8) = au(i)
	 w(i,9) = zu(i)
	enddo
      else	
       CALL IXXXXX(P(0,3   ),ZERO ,-lsign(3),-1,W(1,3   ))        
       CALL OXXXXX(P(0,4   ),ZERO ,+lsign(4),+1,W(1,4   )) 
      endif       
	     
      DO MU = 0,3
         DO NU = 0,3
      
      CALL VCARTX(Q1(0),WMASS ,wwidth, MU,-1,W(1,1   ))        
      CALL VCARTX(Q2(0),WMASS ,wwidth, NU,-1,W(1,2   )) 
      CALL FVIXXX(W(1,3   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,5   ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,2   ),GWF ,AMP(1   ))       
      CALL VVVXXX(W(1,1   ),W(1,2   ),W(1,8   ),GWWA ,AMP(2   ))           
      CALL VVVXXX(W(1,1   ),W(1,2   ),W(1,9   ),GWWZ ,AMP(3   ))                 
      
      wll(mu,nu) = +AMP(   1)+AMP(   2)+AMP(   3)


      if ( (procid.eq.wpzhadjj .or. procid.eq.wmzhadjj) .and. ztouu ) then
 
      CALL FVIXXX(W(1,3   ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,5   ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,1   ),GWF ,AMP(1   ))            
      CALL JVVXXX(W(1,1   ),W(1,2   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     6   ))                                                          
      CALL IOVXXX(W(1,3   ),W(1,4   ),W(1,6   ),GAU ,AMP(2   ))            
      CALL JVVXXX(W(1,1   ),W(1,2   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     7   ))                                                          
      CALL IOVXXX(W(1,3   ),W(1,4   ),W(1,7   ),GZU ,AMP(3   ))            
      wll(mu,nu) = +AMP(   1)+AMP(   2)+AMP(   3)

 
      elseif ( (procid.eq.wpzhadjj .or. procid.eq.wmzhadjj) .and. ztodd ) then

      CALL FVIXXX(W(1,3   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,5   ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,2   ),GWF ,AMP(1   ))            
      CALL JVVXXX(W(1,1   ),W(1,2   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     6   ))                                                          
      CALL IOVXXX(W(1,3   ),W(1,4   ),W(1,6   ),GAD ,AMP(2   ))            
      CALL JVVXXX(W(1,1   ),W(1,2   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     7   ))                                                          
      CALL IOVXXX(W(1,3   ),W(1,4   ),W(1,7   ),GZD ,AMP(3   ))            
      wll(mu,nu) = +AMP(   1)+AMP(   2)+AMP(   3)

      
      endif


	 enddo !nu
      enddo !mu

      end
      
c ------------------------------------------------------------	

c ------------------------------------------------------------	

      subroutine vwtovwp_mg2_had(q1,q2,P,h,zwp,awp)
	
      IMPLICIT NONE

      REAL*8 Q1(0:3), Q2(0:3), P(0:3,3:6)
      complex*16 zwp(0:3,0:3), awp(0:3,0:3)
      integer lsign(3:6,2),ll,mu,nu,h
      

	! possible lepton helicities:
	data (lsign(ll,1),ll=3,6) / -1, -1, -1, -1/  !h=1
	data (lsign(ll,2),ll=3,6) / -1, -1,  1,  1/  !h=2
	
	call zwtovwp1_mg2_had(q1,q2,P,lsign(3,h),zwp(0,0))
	call awtovwp1_mg2_had(q1,q2,P,lsign(3,h),awp(0,0))
	
	! helicity combi nr.h is selected
	
	end
		
c ------------------------------------------------------------	

      subroutine zwtovwp1_mg2_had(q1,q2,P,lsign,zwp)
	
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : z w+ -> e+ ve mu+ mu-  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS 
      PARAMETER (NGRAPHS=  46) 
      INTEGER    NWAVEFUNCS  
      PARAMETER (NWAVEFUNCS=  67) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8  Q1(0:3), Q2(0:3), P(0:3,3:6)
      complex*16 zwp(0:3,0:3)
      integer lsign(3:6)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J,MU, NU
      COMPLEX*16 AMP(NGRAPHS)
      COMPLEX*16 W(6,NWAVEFUNCS)
C  
c  lepton spinors and Z polarization vectors
#include "tensorw.inc"
C  
C GLOBAL VARIABLES
C  
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/process.inc"
   
      logical lswitch
      parameter (lswitch=.true.)
C  
C                                                                          
C ----------
C BEGIN CODE
C ----------

      if (lswitch) then
      	do i = 1,6
         w(i,3) = wl(i)
         w(i,4) = wn(i)
         w(i,5) = lp(i)
         w(i,6) = lm(i)

	 w(i,12) = au(i)
	 w(i,14) = zu(i)
	 w(i,36) = we(i)
      	enddo
      else 
        CALL IXXXXX(P(0,3   ),ZERO ,-lsign(3   ),-1,W(1,3   ))	    
        CALL OXXXXX(P(0,4   ),ZERO ,+lsign(4   ),+1,W(1,4   ))	    
        CALL IXXXXX(P(0,5   ),ZERO ,-lsign(5   ),-1,W(1,5   ))	    
        CALL OXXXXX(P(0,6   ),ZERO ,+lsign(6   ),+1,W(1,6   ))	    
     
        CALL JIOXXX(W(1,5   ),W(1,6   ),GAL ,ZERO    ,ZERO    ,W(1,12  ))    
        CALL JIOXXX(W(1,5   ),W(1,6   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,14  ))    
        CALL JIOXXX(W(1,3   ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,36  ))    
      endif

      DO MU = 0,3
         DO NU = 0,3
       
      CALL VCARTX(Q1(0),ZMASS ,zwidth, MU,-1,W(1,1   ))        
      CALL VCARTX(Q2(0),WMASS ,wwidth, NU,-1,W(1,2   )) 


      if (procid.eq.wphadzjj) then


      CALL FVIXXX(W(1,3   ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,7   ))    
      CALL JIOXXX(W(1,7   ),W(1,4   ),GAU ,ZERO    ,AWIDTH  ,W(1,8   ))    
      CALL FVIXXX(W(1,5   ),W(1,8   ),GAL ,ZERO    ,ZERO    ,W(1,9   ))    
      CALL IOVXXX(W(1,9   ),W(1,6   ),W(1,1   ),GZL ,AMP(1   ))            
      CALL JIOXXX(W(1,7   ),W(1,4   ),GZU ,ZMASS   ,ZWIDTH  ,W(1,10  ))    
      CALL FVIXXX(W(1,5   ),W(1,10  ),GZL ,ZERO    ,ZERO    ,W(1,11  ))    
      CALL IOVXXX(W(1,11  ),W(1,6   ),W(1,1   ),GZL ,AMP(2   ))            
      CALL FVIXXX(W(1,5   ),W(1,1   ),GZL ,ZERO    ,ZERO    ,W(1,12  ))    
      CALL IOVXXX(W(1,12  ),W(1,6   ),W(1,8   ),GAL ,AMP(3   ))            
      CALL IOVXXX(W(1,12  ),W(1,6   ),W(1,10  ),GZL ,AMP(4   ))            
      CALL FVOXXX(W(1,4   ),W(1,1   ),GZU ,ZERO    ,ZERO    ,W(1,13  ))    
      CALL JIOXXX(W(1,3   ),W(1,13  ),GWF ,WMASS   ,WWIDTH  ,W(1,14  ))    
      CALL JVVXXX(W(1,2   ),W(1,14  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     15  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,15  ),GAL ,AMP(5   ))            
      CALL JVVXXX(W(1,2   ),W(1,14  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     16  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,16  ),GZL ,AMP(6   ))            
      CALL FVIXXX(W(1,5   ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,17  ))    
      CALL IOVXXX(W(1,17  ),W(1,6   ),W(1,14  ),GWF ,AMP(7   ))            
      CALL FVOXXX(W(1,13  ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,18  ))    
      CALL JIOXXX(W(1,3   ),W(1,18  ),GAD ,ZERO    ,AWIDTH  ,W(1,19  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,19  ),GAL ,AMP(8   ))            
      CALL JIOXXX(W(1,3   ),W(1,18  ),GZD ,ZMASS   ,ZWIDTH  ,W(1,20  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,20  ),GZL ,AMP(9   ))            
      CALL FVIXXX(W(1,7   ),W(1,1   ),GZU ,ZERO    ,ZERO    ,W(1,21  ))    
      CALL JIOXXX(W(1,21  ),W(1,4   ),GAU ,ZERO    ,AWIDTH  ,W(1,22  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,22  ),GAL ,AMP(10  ))            
      CALL JIOXXX(W(1,21  ),W(1,4   ),GZU ,ZMASS   ,ZWIDTH  ,W(1,23  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,23  ),GZL ,AMP(11  ))            
      CALL JIOXXX(W(1,7   ),W(1,13  ),GAU ,ZERO    ,AWIDTH  ,W(1,24  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,24  ),GAL ,AMP(12  ))            
      CALL JIOXXX(W(1,7   ),W(1,13  ),GZU ,ZMASS   ,ZWIDTH  ,W(1,25  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,25  ),GZL ,AMP(13  ))            
      CALL FVOXXX(W(1,4   ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,26  ))    
      CALL JIOXXX(W(1,3   ),W(1,26  ),GAD ,ZERO    ,AWIDTH  ,W(1,27  ))    
      CALL FVIXXX(W(1,5   ),W(1,27  ),GAL ,ZERO    ,ZERO    ,W(1,28  ))    
      CALL IOVXXX(W(1,28  ),W(1,6   ),W(1,1   ),GZL ,AMP(14  ))            
      CALL JIOXXX(W(1,3   ),W(1,26  ),GZD ,ZMASS   ,ZWIDTH  ,W(1,29  ))    
      CALL FVIXXX(W(1,5   ),W(1,29  ),GZL ,ZERO    ,ZERO    ,W(1,30  ))    
      CALL IOVXXX(W(1,30  ),W(1,6   ),W(1,1   ),GZL ,AMP(15  ))            
      CALL IOVXXX(W(1,12  ),W(1,6   ),W(1,27  ),GAL ,AMP(16  ))            
      CALL IOVXXX(W(1,12  ),W(1,6   ),W(1,29  ),GZL ,AMP(17  ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GZD ,ZERO    ,ZERO    ,W(1,31  ))    
      CALL JIOXXX(W(1,31  ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,32  ))    
      CALL JVVXXX(W(1,2   ),W(1,32  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     33  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,33  ),GAL ,AMP(18  ))            
      CALL JVVXXX(W(1,2   ),W(1,32  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     34  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,34  ),GZL ,AMP(19  ))            
      CALL IOVXXX(W(1,17  ),W(1,6   ),W(1,32  ),GWF ,AMP(20  ))            
      CALL FVOXXX(W(1,26  ),W(1,1   ),GZD ,ZERO    ,ZERO    ,W(1,35  ))    
      CALL JIOXXX(W(1,3   ),W(1,35  ),GAD ,ZERO    ,AWIDTH  ,W(1,36  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,36  ),GAL ,AMP(21  ))            
      CALL JIOXXX(W(1,3   ),W(1,35  ),GZD ,ZMASS   ,ZWIDTH  ,W(1,37  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,37  ),GZL ,AMP(22  ))            
      CALL FVIXXX(W(1,31  ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,38  ))    
      CALL JIOXXX(W(1,38  ),W(1,4   ),GAU ,ZERO    ,AWIDTH  ,W(1,39  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,39  ),GAL ,AMP(23  ))            
      CALL JIOXXX(W(1,38  ),W(1,4   ),GZU ,ZMASS   ,ZWIDTH  ,W(1,40  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,40  ),GZL ,AMP(24  ))            
      CALL JIOXXX(W(1,31  ),W(1,26  ),GAD ,ZERO    ,AWIDTH  ,W(1,41  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,41  ),GAL ,AMP(25  ))            
      CALL JIOXXX(W(1,31  ),W(1,26  ),GZD ,ZMASS   ,ZWIDTH  ,W(1,42  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,42  ),GZL ,AMP(26  ))            
      CALL JIOXXX(W(1,3   ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,43  ))    
      CALL JVVXXX(W(1,2   ),W(1,43  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     44  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,44  ),GAL ,ZERO    ,ZERO    ,W(1,45  ))    
      CALL IOVXXX(W(1,45  ),W(1,6   ),W(1,1   ),GZL ,AMP(27  ))            
      CALL JVVXXX(W(1,2   ),W(1,43  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     46  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,46  ),GZL ,ZERO    ,ZERO    ,W(1,47  ))    
      CALL IOVXXX(W(1,47  ),W(1,6   ),W(1,1   ),GZL ,AMP(28  ))            
      CALL HVVXXX(W(1,2   ),W(1,43  ),GWWH ,HMASS   ,HWIDTH  ,W(1,         
     &     48  ))                                                          
      CALL JVSXXX(W(1,1   ),W(1,48  ),GZZH ,ZMASS   ,ZWIDTH  ,W(1,         
     &     49  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,49  ),GZL ,AMP(29  ))            
      CALL FVIXXX(W(1,12  ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,50  ))    
      CALL IOVXXX(W(1,50  ),W(1,6   ),W(1,43  ),GWF ,AMP(30  ))            
      CALL IOVXXX(W(1,12  ),W(1,6   ),W(1,44  ),GAL ,AMP(31  ))            
      CALL IOVXXX(W(1,12  ),W(1,6   ),W(1,46  ),GZL ,AMP(32  ))            
      CALL FVIXXX(W(1,17  ),W(1,43  ),GWF ,ZERO    ,ZERO    ,W(1,51  ))    
      CALL IOVXXX(W(1,51  ),W(1,6   ),W(1,1   ),GZL ,AMP(33  ))            
      CALL JVVXXX(W(1,43  ),W(1,1   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     52  ))                                                          
      CALL JVVXXX(W(1,2   ),W(1,52  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     53  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,53  ),GAL ,AMP(34  ))            
      CALL JVVXXX(W(1,2   ),W(1,52  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     54  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,54  ),GZL ,AMP(35  ))            
      CALL FVIXXX(W(1,17  ),W(1,1   ),GZN ,ZERO    ,ZERO    ,W(1,55  ))    
      CALL IOVXXX(W(1,55  ),W(1,6   ),W(1,43  ),GWF ,AMP(36  ))            
      CALL IOVXXX(W(1,17  ),W(1,6   ),W(1,52  ),GWF ,AMP(37  ))            
      CALL JVVXXX(W(1,1   ),W(1,2   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     56  ))                                                          
      CALL FVOXXX(W(1,4   ),W(1,56  ),GWF ,ZERO    ,ZERO    ,W(1,57  ))    
      CALL JIOXXX(W(1,3   ),W(1,57  ),GAD ,ZERO    ,AWIDTH  ,W(1,58  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,58  ),GAL ,AMP(38  ))            
      CALL JIOXXX(W(1,3   ),W(1,57  ),GZD ,ZMASS   ,ZWIDTH  ,W(1,59  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,59  ),GZL ,AMP(39  ))            
      CALL FVIXXX(W(1,3   ),W(1,56  ),GWF ,ZERO    ,ZERO    ,W(1,60  ))    
      CALL JIOXXX(W(1,60  ),W(1,4   ),GAU ,ZERO    ,AWIDTH  ,W(1,61  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,61  ),GAL ,AMP(40  ))            
      CALL JIOXXX(W(1,60  ),W(1,4   ),GZU ,ZMASS   ,ZWIDTH  ,W(1,62  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,62  ),GZL ,AMP(41  ))            
      CALL JW3WXX(W(1,2   ),W(1,1   ),W(1,43  ),GWWZ ,GWWA ,ZERO    ,      
     &     AWIDTH  ,W(1,63  ))                                             
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,63  ),GAL ,AMP(42  ))            
      CALL JW3WXX(W(1,2   ),W(1,1   ),W(1,43  ),GWWZ ,GWWZ ,ZMASS   ,      
     &     ZWIDTH  ,W(1,64  ))                                             
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,64  ),GZL ,AMP(43  ))            
      CALL JVVXXX(W(1,56  ),W(1,43  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     65  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,65  ),GAL ,AMP(44  ))            
      CALL JVVXXX(W(1,56  ),W(1,43  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     66  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,66  ),GZL ,AMP(45  ))            
      CALL FVIXXX(W(1,5   ),W(1,56  ),GWF ,ZERO    ,ZERO    ,W(1,67  ))    
      CALL IOVXXX(W(1,67  ),W(1,6   ),W(1,43  ),GWF ,AMP(46  ))            
      zwp(mu,nu) = +AMP(   1)+AMP(   2)+AMP(   3)+AMP(   4)+AMP(   5)
     &             +AMP(   6)+AMP(   7)+AMP(   8)+AMP(   9)+AMP(  10)
     &             +AMP(  11)+AMP(  12)+AMP(  13)+AMP(  14)+AMP(  15)
     &             +AMP(  16)+AMP(  17)+AMP(  18)+AMP(  19)+AMP(  20)
     &             +AMP(  21)+AMP(  22)+AMP(  23)+AMP(  24)+AMP(  25)
     &             +AMP(  26)+AMP(  27)+AMP(  28)+AMP(  29)+AMP(  30)
     &             +AMP(  31)+AMP(  32)+AMP(  33)+AMP(  34)+AMP(  35)
     &             +AMP(  36)+AMP(  37)+AMP(  38)+AMP(  39)+AMP(  40)
     &             +AMP(  41)+AMP(  42)+AMP(  43)+AMP(  44)+AMP(  45)
     &             +AMP(  46)


      elseif (procid.eq.wpzhadjj .and. ztouu) then


      CALL FVIXXX(W(1,3   ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,7   ))    
      CALL JIOXXX(W(1,7   ),W(1,4   ),GZN ,ZMASS   ,ZWIDTH  ,W(1,8   ))    
      CALL FVIXXX(W(1,5   ),W(1,8   ),GZU ,ZERO    ,ZERO    ,W(1,9   ))    
      CALL IOVXXX(W(1,9   ),W(1,6   ),W(1,1   ),GZU ,AMP(1   ))            
      CALL FVIXXX(W(1,5   ),W(1,1   ),GZU ,ZERO    ,ZERO    ,W(1,10  ))    
      CALL IOVXXX(W(1,10  ),W(1,6   ),W(1,8   ),GZU ,AMP(2   ))            
      CALL FVOXXX(W(1,4   ),W(1,1   ),GZN ,ZERO    ,ZERO    ,W(1,11  ))    
      CALL JIOXXX(W(1,3   ),W(1,11  ),GWF ,WMASS   ,WWIDTH  ,W(1,12  ))    
      CALL FVIXXX(W(1,5   ),W(1,12  ),GWF ,ZERO    ,ZERO    ,W(1,13  ))    
      CALL IOVXXX(W(1,13  ),W(1,6   ),W(1,2   ),GWF ,AMP(3   ))            
      CALL JVVXXX(W(1,2   ),W(1,12  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     14  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,14  ),GAU ,AMP(4   ))            
      CALL JVVXXX(W(1,2   ),W(1,12  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     15  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,15  ),GZU ,AMP(5   ))            
      CALL FVOXXX(W(1,11  ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,16  ))    
      CALL JIOXXX(W(1,3   ),W(1,16  ),GAL ,ZERO    ,AWIDTH  ,W(1,17  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,17  ),GAU ,AMP(6   ))            
      CALL JIOXXX(W(1,3   ),W(1,16  ),GZL ,ZMASS   ,ZWIDTH  ,W(1,18  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,18  ),GZU ,AMP(7   ))            
      CALL FVIXXX(W(1,7   ),W(1,1   ),GZN ,ZERO    ,ZERO    ,W(1,19  ))    
      CALL JIOXXX(W(1,19  ),W(1,4   ),GZN ,ZMASS   ,ZWIDTH  ,W(1,20  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,20  ),GZU ,AMP(8   ))            
      CALL JIOXXX(W(1,7   ),W(1,11  ),GZN ,ZMASS   ,ZWIDTH  ,W(1,21  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,21  ),GZU ,AMP(9   ))            
      CALL FVOXXX(W(1,4   ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,22  ))    
      CALL JIOXXX(W(1,3   ),W(1,22  ),GAL ,ZERO    ,AWIDTH  ,W(1,23  ))    
      CALL FVIXXX(W(1,5   ),W(1,23  ),GAU ,ZERO    ,ZERO    ,W(1,24  ))    
      CALL IOVXXX(W(1,24  ),W(1,6   ),W(1,1   ),GZU ,AMP(10  ))            
      CALL JIOXXX(W(1,3   ),W(1,22  ),GZL ,ZMASS   ,ZWIDTH  ,W(1,25  ))    
      CALL FVIXXX(W(1,5   ),W(1,25  ),GZU ,ZERO    ,ZERO    ,W(1,26  ))    
      CALL IOVXXX(W(1,26  ),W(1,6   ),W(1,1   ),GZU ,AMP(11  ))            
      CALL IOVXXX(W(1,10  ),W(1,6   ),W(1,23  ),GAU ,AMP(12  ))            
      CALL IOVXXX(W(1,10  ),W(1,6   ),W(1,25  ),GZU ,AMP(13  ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GZL ,ZERO    ,ZERO    ,W(1,27  ))    
      CALL JIOXXX(W(1,27  ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,28  ))    
      CALL FVIXXX(W(1,5   ),W(1,28  ),GWF ,ZERO    ,ZERO    ,W(1,29  ))    
      CALL IOVXXX(W(1,29  ),W(1,6   ),W(1,2   ),GWF ,AMP(14  ))            
      CALL JVVXXX(W(1,2   ),W(1,28  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     30  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,30  ),GAU ,AMP(15  ))            
      CALL JVVXXX(W(1,2   ),W(1,28  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     31  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,31  ),GZU ,AMP(16  ))            
      CALL FVOXXX(W(1,22  ),W(1,1   ),GZL ,ZERO    ,ZERO    ,W(1,32  ))    
      CALL JIOXXX(W(1,3   ),W(1,32  ),GAL ,ZERO    ,AWIDTH  ,W(1,33  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,33  ),GAU ,AMP(17  ))            
      CALL JIOXXX(W(1,3   ),W(1,32  ),GZL ,ZMASS   ,ZWIDTH  ,W(1,34  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,34  ),GZU ,AMP(18  ))            
      CALL FVIXXX(W(1,27  ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,35  ))    
      CALL JIOXXX(W(1,35  ),W(1,4   ),GZN ,ZMASS   ,ZWIDTH  ,W(1,36  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,36  ),GZU ,AMP(19  ))            
      CALL JIOXXX(W(1,27  ),W(1,22  ),GAL ,ZERO    ,AWIDTH  ,W(1,37  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,37  ),GAU ,AMP(20  ))            
      CALL JIOXXX(W(1,27  ),W(1,22  ),GZL ,ZMASS   ,ZWIDTH  ,W(1,38  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,38  ),GZU ,AMP(21  ))            
      CALL JIOXXX(W(1,3   ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,39  ))    
      CALL JVVXXX(W(1,2   ),W(1,39  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     40  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,40  ),GAU ,ZERO    ,ZERO    ,W(1,41  ))    
      CALL IOVXXX(W(1,41  ),W(1,6   ),W(1,1   ),GZU ,AMP(22  ))            
      CALL JVVXXX(W(1,2   ),W(1,39  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     42  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,42  ),GZU ,ZERO    ,ZERO    ,W(1,43  ))    
      CALL IOVXXX(W(1,43  ),W(1,6   ),W(1,1   ),GZU ,AMP(23  ))            
      CALL FVIXXX(W(1,10  ),W(1,39  ),GWF ,ZERO    ,ZERO    ,W(1,44  ))    
      CALL IOVXXX(W(1,44  ),W(1,6   ),W(1,2   ),GWF ,AMP(24  ))            
      CALL HVVXXX(W(1,2   ),W(1,39  ),GWWH ,HMASS   ,HWIDTH  ,W(1,         
     &     45  ))                                                          
      CALL JVSXXX(W(1,1   ),W(1,45  ),GZZH ,ZMASS   ,ZWIDTH  ,W(1,         
     &     46  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,46  ),GZU ,AMP(25  ))            
      CALL IOVXXX(W(1,10  ),W(1,6   ),W(1,40  ),GAU ,AMP(26  ))            
      CALL IOVXXX(W(1,10  ),W(1,6   ),W(1,42  ),GZU ,AMP(27  ))            
      CALL JVVXXX(W(1,39  ),W(1,1   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     47  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,47  ),GWF ,ZERO    ,ZERO    ,W(1,48  ))    
      CALL IOVXXX(W(1,48  ),W(1,6   ),W(1,2   ),GWF ,AMP(28  ))            
      CALL JVVXXX(W(1,2   ),W(1,47  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     49  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,49  ),GAU ,AMP(29  ))            
      CALL JVVXXX(W(1,2   ),W(1,47  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     50  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,50  ),GZU ,AMP(30  ))            
      CALL JVVXXX(W(1,1   ),W(1,2   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     51  ))                                                          
      CALL FVOXXX(W(1,4   ),W(1,51  ),GWF ,ZERO    ,ZERO    ,W(1,52  ))    
      CALL JIOXXX(W(1,3   ),W(1,52  ),GAL ,ZERO    ,AWIDTH  ,W(1,53  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,53  ),GAU ,AMP(31  ))            
      CALL JIOXXX(W(1,3   ),W(1,52  ),GZL ,ZMASS   ,ZWIDTH  ,W(1,54  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,54  ),GZU ,AMP(32  ))            
      CALL FVIXXX(W(1,3   ),W(1,51  ),GWF ,ZERO    ,ZERO    ,W(1,55  ))    
      CALL JIOXXX(W(1,55  ),W(1,4   ),GZN ,ZMASS   ,ZWIDTH  ,W(1,56  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,56  ),GZU ,AMP(33  ))            
      CALL JW3WXX(W(1,2   ),W(1,1   ),W(1,39  ),GWWZ ,GWWA ,ZERO    ,      
     &     AWIDTH  ,W(1,57  ))                                             
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,57  ),GAU ,AMP(34  ))            
      CALL JW3WXX(W(1,2   ),W(1,1   ),W(1,39  ),GWWZ ,GWWZ ,ZMASS   ,      
     &     ZWIDTH  ,W(1,58  ))                                             
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,58  ),GZU ,AMP(35  ))            
      CALL FVIXXX(W(1,5   ),W(1,39  ),GWF ,ZERO    ,ZERO    ,W(1,59  ))    
      CALL FVIXXX(W(1,59  ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,60  ))    
      CALL IOVXXX(W(1,60  ),W(1,6   ),W(1,1   ),GZU ,AMP(36  ))            
      CALL FVIXXX(W(1,59  ),W(1,1   ),GZD ,ZERO    ,ZERO    ,W(1,61  ))    
      CALL IOVXXX(W(1,61  ),W(1,6   ),W(1,2   ),GWF ,AMP(37  ))            
      CALL JVVXXX(W(1,51  ),W(1,39  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     62  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,62  ),GAU ,AMP(38  ))            
      CALL JVVXXX(W(1,51  ),W(1,39  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     63  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,63  ),GZU ,AMP(39  ))            
      CALL IOVXXX(W(1,59  ),W(1,6   ),W(1,51  ),GWF ,AMP(40  ))            
      zwp(mu,nu) = +AMP(   1)+AMP(   2)+AMP(   3)+AMP(   4)+AMP(   5)
     &             +AMP(   6)+AMP(   7)+AMP(   8)+AMP(   9)+AMP(  10)
     &             +AMP(  11)+AMP(  12)+AMP(  13)+AMP(  14)+AMP(  15)
     &             +AMP(  16)+AMP(  17)+AMP(  18)+AMP(  19)+AMP(  20)
     &             +AMP(  21)+AMP(  22)+AMP(  23)+AMP(  24)+AMP(  25)
     &             +AMP(  26)+AMP(  27)+AMP(  28)+AMP(  29)+AMP(  30)
     &             +AMP(  31)+AMP(  32)+AMP(  33)+AMP(  34)+AMP(  35)
     &             +AMP(  36)+AMP(  37)+AMP(  38)+AMP(  39)+AMP(  40)


      elseif (procid.eq.wpzhadjj .and. ztodd) then


      CALL FVIXXX(W(1,3   ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,7   ))    
      CALL JIOXXX(W(1,7   ),W(1,4   ),GZN ,ZMASS   ,ZWIDTH  ,W(1,8   ))    
      CALL FVIXXX(W(1,5   ),W(1,8   ),GZD ,ZERO    ,ZERO    ,W(1,9   ))    
      CALL IOVXXX(W(1,9   ),W(1,6   ),W(1,1   ),GZD ,AMP(1   ))            
      CALL FVIXXX(W(1,5   ),W(1,1   ),GZD ,ZERO    ,ZERO    ,W(1,10  ))    
      CALL IOVXXX(W(1,10  ),W(1,6   ),W(1,8   ),GZD ,AMP(2   ))            
      CALL FVOXXX(W(1,4   ),W(1,1   ),GZN ,ZERO    ,ZERO    ,W(1,11  ))    
      CALL JIOXXX(W(1,3   ),W(1,11  ),GWF ,WMASS   ,WWIDTH  ,W(1,12  ))    
      CALL JVVXXX(W(1,2   ),W(1,12  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     13  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,13  ),GAD ,AMP(3   ))            
      CALL JVVXXX(W(1,2   ),W(1,12  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     14  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,14  ),GZD ,AMP(4   ))            
      CALL FVIXXX(W(1,5   ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,15  ))    
      CALL IOVXXX(W(1,15  ),W(1,6   ),W(1,12  ),GWF ,AMP(5   ))            
      CALL FVOXXX(W(1,11  ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,16  ))    
      CALL JIOXXX(W(1,3   ),W(1,16  ),GAL ,ZERO    ,AWIDTH  ,W(1,17  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,17  ),GAD ,AMP(6   ))            
      CALL JIOXXX(W(1,3   ),W(1,16  ),GZL ,ZMASS   ,ZWIDTH  ,W(1,18  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,18  ),GZD ,AMP(7   ))            
      CALL FVIXXX(W(1,7   ),W(1,1   ),GZN ,ZERO    ,ZERO    ,W(1,19  ))    
      CALL JIOXXX(W(1,19  ),W(1,4   ),GZN ,ZMASS   ,ZWIDTH  ,W(1,20  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,20  ),GZD ,AMP(8   ))            
      CALL JIOXXX(W(1,7   ),W(1,11  ),GZN ,ZMASS   ,ZWIDTH  ,W(1,21  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,21  ),GZD ,AMP(9   ))            
      CALL FVOXXX(W(1,4   ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,22  ))    
      CALL JIOXXX(W(1,3   ),W(1,22  ),GAL ,ZERO    ,AWIDTH  ,W(1,23  ))    
      CALL FVIXXX(W(1,5   ),W(1,23  ),GAD ,ZERO    ,ZERO    ,W(1,24  ))    
      CALL IOVXXX(W(1,24  ),W(1,6   ),W(1,1   ),GZD ,AMP(10  ))            
      CALL JIOXXX(W(1,3   ),W(1,22  ),GZL ,ZMASS   ,ZWIDTH  ,W(1,25  ))    
      CALL FVIXXX(W(1,5   ),W(1,25  ),GZD ,ZERO    ,ZERO    ,W(1,26  ))    
      CALL IOVXXX(W(1,26  ),W(1,6   ),W(1,1   ),GZD ,AMP(11  ))            
      CALL IOVXXX(W(1,10  ),W(1,6   ),W(1,23  ),GAD ,AMP(12  ))            
      CALL IOVXXX(W(1,10  ),W(1,6   ),W(1,25  ),GZD ,AMP(13  ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GZL ,ZERO    ,ZERO    ,W(1,27  ))    
      CALL JIOXXX(W(1,27  ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,28  ))    
      CALL JVVXXX(W(1,2   ),W(1,28  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     29  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,29  ),GAD ,AMP(14  ))            
      CALL JVVXXX(W(1,2   ),W(1,28  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     30  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,30  ),GZD ,AMP(15  ))            
      CALL IOVXXX(W(1,15  ),W(1,6   ),W(1,28  ),GWF ,AMP(16  ))            
      CALL FVOXXX(W(1,22  ),W(1,1   ),GZL ,ZERO    ,ZERO    ,W(1,31  ))    
      CALL JIOXXX(W(1,3   ),W(1,31  ),GAL ,ZERO    ,AWIDTH  ,W(1,32  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,32  ),GAD ,AMP(17  ))            
      CALL JIOXXX(W(1,3   ),W(1,31  ),GZL ,ZMASS   ,ZWIDTH  ,W(1,33  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,33  ),GZD ,AMP(18  ))            
      CALL FVIXXX(W(1,27  ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,34  ))    
      CALL JIOXXX(W(1,34  ),W(1,4   ),GZN ,ZMASS   ,ZWIDTH  ,W(1,35  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,35  ),GZD ,AMP(19  ))            
      CALL JIOXXX(W(1,27  ),W(1,22  ),GAL ,ZERO    ,AWIDTH  ,W(1,36  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,36  ),GAD ,AMP(20  ))            
      CALL JIOXXX(W(1,27  ),W(1,22  ),GZL ,ZMASS   ,ZWIDTH  ,W(1,37  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,37  ),GZD ,AMP(21  ))            
      CALL JIOXXX(W(1,3   ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,38  ))    
      CALL JVVXXX(W(1,2   ),W(1,38  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     39  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,39  ),GAD ,ZERO    ,ZERO    ,W(1,40  ))    
      CALL IOVXXX(W(1,40  ),W(1,6   ),W(1,1   ),GZD ,AMP(22  ))            
      CALL JVVXXX(W(1,2   ),W(1,38  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     41  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,41  ),GZD ,ZERO    ,ZERO    ,W(1,42  ))    
      CALL IOVXXX(W(1,42  ),W(1,6   ),W(1,1   ),GZD ,AMP(23  ))            
      CALL HVVXXX(W(1,2   ),W(1,38  ),GWWH ,HMASS   ,HWIDTH  ,W(1,         
     &     43  ))                                                          
      CALL JVSXXX(W(1,1   ),W(1,43  ),GZZH ,ZMASS   ,ZWIDTH  ,W(1,         
     &     44  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,44  ),GZD ,AMP(24  ))            
      CALL FVIXXX(W(1,10  ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,45  ))    
      CALL IOVXXX(W(1,45  ),W(1,6   ),W(1,38  ),GWF ,AMP(25  ))            
      CALL IOVXXX(W(1,10  ),W(1,6   ),W(1,39  ),GAD ,AMP(26  ))            
      CALL IOVXXX(W(1,10  ),W(1,6   ),W(1,41  ),GZD ,AMP(27  ))            
      CALL FVIXXX(W(1,15  ),W(1,38  ),GWF ,ZERO    ,ZERO    ,W(1,46  ))    
      CALL IOVXXX(W(1,46  ),W(1,6   ),W(1,1   ),GZD ,AMP(28  ))            
      CALL JVVXXX(W(1,38  ),W(1,1   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     47  ))                                                          
      CALL JVVXXX(W(1,2   ),W(1,47  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     48  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,48  ),GAD ,AMP(29  ))            
      CALL JVVXXX(W(1,2   ),W(1,47  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     49  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,49  ),GZD ,AMP(30  ))            
      CALL FVIXXX(W(1,15  ),W(1,1   ),GZU ,ZERO    ,ZERO    ,W(1,50  ))    
      CALL IOVXXX(W(1,50  ),W(1,6   ),W(1,38  ),GWF ,AMP(31  ))            
      CALL IOVXXX(W(1,15  ),W(1,6   ),W(1,47  ),GWF ,AMP(32  ))            
      CALL JVVXXX(W(1,1   ),W(1,2   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     51  ))                                                          
      CALL FVOXXX(W(1,4   ),W(1,51  ),GWF ,ZERO    ,ZERO    ,W(1,52  ))    
      CALL JIOXXX(W(1,3   ),W(1,52  ),GAL ,ZERO    ,AWIDTH  ,W(1,53  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,53  ),GAD ,AMP(33  ))            
      CALL JIOXXX(W(1,3   ),W(1,52  ),GZL ,ZMASS   ,ZWIDTH  ,W(1,54  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,54  ),GZD ,AMP(34  ))            
      CALL FVIXXX(W(1,3   ),W(1,51  ),GWF ,ZERO    ,ZERO    ,W(1,55  ))    
      CALL JIOXXX(W(1,55  ),W(1,4   ),GZN ,ZMASS   ,ZWIDTH  ,W(1,56  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,56  ),GZD ,AMP(35  ))            
      CALL JW3WXX(W(1,2   ),W(1,1   ),W(1,38  ),GWWZ ,GWWA ,ZERO    ,      
     &     AWIDTH  ,W(1,57  ))                                             
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,57  ),GAD ,AMP(36  ))            
      CALL JW3WXX(W(1,2   ),W(1,1   ),W(1,38  ),GWWZ ,GWWZ ,ZMASS   ,      
     &     ZWIDTH  ,W(1,58  ))                                             
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,58  ),GZD ,AMP(37  ))            
      CALL JVVXXX(W(1,51  ),W(1,38  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     59  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,59  ),GAD ,AMP(38  ))            
      CALL JVVXXX(W(1,51  ),W(1,38  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     60  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,60  ),GZD ,AMP(39  ))            
      CALL FVIXXX(W(1,5   ),W(1,51  ),GWF ,ZERO    ,ZERO    ,W(1,61  ))    
      CALL IOVXXX(W(1,61  ),W(1,6   ),W(1,38  ),GWF ,AMP(40  ))            
      zwp(mu,nu) = +AMP(   1)+AMP(   2)+AMP(   3)+AMP(   4)+AMP(   5)
     &             +AMP(   6)+AMP(   7)+AMP(   8)+AMP(   9)+AMP(  10)
     &             +AMP(  11)+AMP(  12)+AMP(  13)+AMP(  14)+AMP(  15)
     &             +AMP(  16)+AMP(  17)+AMP(  18)+AMP(  19)+AMP(  20)
     &             +AMP(  21)+AMP(  22)+AMP(  23)+AMP(  24)+AMP(  25)
     &             +AMP(  26)+AMP(  27)+AMP(  28)+AMP(  29)+AMP(  30)
     &             +AMP(  31)+AMP(  32)+AMP(  33)+AMP(  34)+AMP(  35)
     &             +AMP(  36)+AMP(  37)+AMP(  38)+AMP(  39)+AMP(  40)


      endif
          
          ENDDO
      ENDDO
 
      END

c ------------------------------------------------------------	

      subroutine awtovwp1_mg2_had(q1,q2,P,lsign,awp)
	  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : a w+ -> e+ ve mu+ mu-  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=  44) 
c#include "genps.inc"
      INTEGER    NWAVEFUNCS    
      PARAMETER (NWAVEFUNCS=  64) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 Q1(0:3), Q2(0:3), P(0:3,3:6)
      complex*16 awp(0:3,0:3)
      integer lsign(3:6)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J,MU, NU
      COMPLEX*16 AMP(NGRAPHS)
      COMPLEX*16 W(6,NWAVEFUNCS)
C  
c  lepton spinors and Z polarization vectors
#include "tensorw.inc"
C  
C GLOBAL VARIABLES
c
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/process.inc"
   
      logical lswitch
      parameter (lswitch=.true.)
C  
C                                                                          
C ----------
C BEGIN CODE
C ----------

      if (lswitch) then
      	do i = 1,6
         w(i,3) = wl(i)
         w(i,4) = wn(i)
         w(i,5) = lp(i)
         w(i,6) = lm(i)

	 w(i,22) = au(i)
	 w(i,24) = zu(i)
	 w(i,29) = we(i)
      	enddo
      else 
	CALL IXXXXX(P(0,3   ),ZERO ,-lsign(3   ),-1,W(1,3   ))	    
      	CALL OXXXXX(P(0,4   ),ZERO ,+lsign(4   ),+1,W(1,4   ))	    
      	CALL IXXXXX(P(0,5   ),ZERO ,-lsign(5   ),-1,W(1,5   ))	    
      	CALL OXXXXX(P(0,6   ),ZERO ,+lsign(6   ),+1,W(1,6   ))
     
      	CALL JIOXXX(W(1,5   ),W(1,6   ),GAL ,ZERO    ,ZERO    ,W(1,22  ))    
      	CALL JIOXXX(W(1,5   ),W(1,6   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,24  ))    
      	CALL JIOXXX(W(1,3   ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,29  ))         
      endif	    

      DO MU = 0,3
         DO NU = 0,3
       
      CALL VCARTX(Q1(0),ZERO ,zero, MU,-1,W(1,1   ))        
      CALL VCARTX(Q2(0),WMASS ,wwidth, NU,-1,W(1,2   )) 

      if (procid.eq.wphadzjj) then

      CALL FVIXXX(W(1,3   ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,7   ))    
      CALL JIOXXX(W(1,7   ),W(1,4   ),GAU ,ZERO    ,AWIDTH  ,W(1,8   ))    
      CALL FVIXXX(W(1,5   ),W(1,8   ),GAL ,ZERO    ,ZERO    ,W(1,9   ))    
      CALL IOVXXX(W(1,9   ),W(1,6   ),W(1,1   ),GAL ,AMP(1   ))            
      CALL JIOXXX(W(1,7   ),W(1,4   ),GZU ,ZMASS   ,ZWIDTH  ,W(1,10  ))    
      CALL FVIXXX(W(1,5   ),W(1,10  ),GZL ,ZERO    ,ZERO    ,W(1,11  ))    
      CALL IOVXXX(W(1,11  ),W(1,6   ),W(1,1   ),GAL ,AMP(2   ))            
      CALL FVIXXX(W(1,5   ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,12  ))    
      CALL IOVXXX(W(1,12  ),W(1,6   ),W(1,8   ),GAL ,AMP(3   ))            
      CALL IOVXXX(W(1,12  ),W(1,6   ),W(1,10  ),GZL ,AMP(4   ))            
      CALL FVOXXX(W(1,4   ),W(1,1   ),GAU ,ZERO    ,ZERO    ,W(1,13  ))    
      CALL JIOXXX(W(1,3   ),W(1,13  ),GWF ,WMASS   ,WWIDTH  ,W(1,14  ))    
      CALL JVVXXX(W(1,2   ),W(1,14  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     15  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,15  ),GAL ,AMP(5   ))            
      CALL JVVXXX(W(1,2   ),W(1,14  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     16  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,16  ),GZL ,AMP(6   ))            
      CALL FVIXXX(W(1,5   ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,17  ))    
      CALL IOVXXX(W(1,17  ),W(1,6   ),W(1,14  ),GWF ,AMP(7   ))            
      CALL FVOXXX(W(1,13  ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,18  ))    
      CALL JIOXXX(W(1,3   ),W(1,18  ),GAD ,ZERO    ,AWIDTH  ,W(1,19  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,19  ),GAL ,AMP(8   ))            
      CALL JIOXXX(W(1,3   ),W(1,18  ),GZD ,ZMASS   ,ZWIDTH  ,W(1,20  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,20  ),GZL ,AMP(9   ))            
      CALL FVIXXX(W(1,7   ),W(1,1   ),GAU ,ZERO    ,ZERO    ,W(1,21  ))    
      CALL JIOXXX(W(1,21  ),W(1,4   ),GAU ,ZERO    ,AWIDTH  ,W(1,22  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,22  ),GAL ,AMP(10  ))            
      CALL JIOXXX(W(1,21  ),W(1,4   ),GZU ,ZMASS   ,ZWIDTH  ,W(1,23  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,23  ),GZL ,AMP(11  ))            
      CALL JIOXXX(W(1,7   ),W(1,13  ),GAU ,ZERO    ,AWIDTH  ,W(1,24  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,24  ),GAL ,AMP(12  ))            
      CALL JIOXXX(W(1,7   ),W(1,13  ),GZU ,ZMASS   ,ZWIDTH  ,W(1,25  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,25  ),GZL ,AMP(13  ))            
      CALL FVOXXX(W(1,4   ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,26  ))    
      CALL JIOXXX(W(1,3   ),W(1,26  ),GAD ,ZERO    ,AWIDTH  ,W(1,27  ))    
      CALL FVIXXX(W(1,5   ),W(1,27  ),GAL ,ZERO    ,ZERO    ,W(1,28  ))    
      CALL IOVXXX(W(1,28  ),W(1,6   ),W(1,1   ),GAL ,AMP(14  ))            
      CALL JIOXXX(W(1,3   ),W(1,26  ),GZD ,ZMASS   ,ZWIDTH  ,W(1,29  ))    
      CALL FVIXXX(W(1,5   ),W(1,29  ),GZL ,ZERO    ,ZERO    ,W(1,30  ))    
      CALL IOVXXX(W(1,30  ),W(1,6   ),W(1,1   ),GAL ,AMP(15  ))            
      CALL IOVXXX(W(1,12  ),W(1,6   ),W(1,27  ),GAL ,AMP(16  ))            
      CALL IOVXXX(W(1,12  ),W(1,6   ),W(1,29  ),GZL ,AMP(17  ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GAD ,ZERO    ,ZERO    ,W(1,31  ))    
      CALL JIOXXX(W(1,31  ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,32  ))    
      CALL JVVXXX(W(1,2   ),W(1,32  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     33  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,33  ),GAL ,AMP(18  ))            
      CALL JVVXXX(W(1,2   ),W(1,32  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     34  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,34  ),GZL ,AMP(19  ))            
      CALL IOVXXX(W(1,17  ),W(1,6   ),W(1,32  ),GWF ,AMP(20  ))            
      CALL FVOXXX(W(1,26  ),W(1,1   ),GAD ,ZERO    ,ZERO    ,W(1,35  ))    
      CALL JIOXXX(W(1,3   ),W(1,35  ),GAD ,ZERO    ,AWIDTH  ,W(1,36  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,36  ),GAL ,AMP(21  ))            
      CALL JIOXXX(W(1,3   ),W(1,35  ),GZD ,ZMASS   ,ZWIDTH  ,W(1,37  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,37  ),GZL ,AMP(22  ))            
      CALL FVIXXX(W(1,31  ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,38  ))    
      CALL JIOXXX(W(1,38  ),W(1,4   ),GAU ,ZERO    ,AWIDTH  ,W(1,39  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,39  ),GAL ,AMP(23  ))            
      CALL JIOXXX(W(1,38  ),W(1,4   ),GZU ,ZMASS   ,ZWIDTH  ,W(1,40  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,40  ),GZL ,AMP(24  ))            
      CALL JIOXXX(W(1,31  ),W(1,26  ),GAD ,ZERO    ,AWIDTH  ,W(1,41  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,41  ),GAL ,AMP(25  ))            
      CALL JIOXXX(W(1,31  ),W(1,26  ),GZD ,ZMASS   ,ZWIDTH  ,W(1,42  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,42  ),GZL ,AMP(26  ))            
      CALL JIOXXX(W(1,3   ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,43  ))    
      CALL JVVXXX(W(1,2   ),W(1,43  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     44  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,44  ),GAL ,ZERO    ,ZERO    ,W(1,45  ))    
      CALL IOVXXX(W(1,45  ),W(1,6   ),W(1,1   ),GAL ,AMP(27  ))            
      CALL JVVXXX(W(1,2   ),W(1,43  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     46  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,46  ),GZL ,ZERO    ,ZERO    ,W(1,47  ))    
      CALL IOVXXX(W(1,47  ),W(1,6   ),W(1,1   ),GAL ,AMP(28  ))            
      CALL FVIXXX(W(1,12  ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,48  ))    
      CALL IOVXXX(W(1,48  ),W(1,6   ),W(1,43  ),GWF ,AMP(29  ))            
      CALL IOVXXX(W(1,12  ),W(1,6   ),W(1,44  ),GAL ,AMP(30  ))            
      CALL IOVXXX(W(1,12  ),W(1,6   ),W(1,46  ),GZL ,AMP(31  ))            
      CALL FVIXXX(W(1,17  ),W(1,43  ),GWF ,ZERO    ,ZERO    ,W(1,49  ))    
      CALL IOVXXX(W(1,49  ),W(1,6   ),W(1,1   ),GAL ,AMP(32  ))            
      CALL JVVXXX(W(1,43  ),W(1,1   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     50  ))                                                          
      CALL JVVXXX(W(1,2   ),W(1,50  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     51  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,51  ),GAL ,AMP(33  ))            
      CALL JVVXXX(W(1,2   ),W(1,50  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     52  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,52  ),GZL ,AMP(34  ))            
      CALL IOVXXX(W(1,17  ),W(1,6   ),W(1,50  ),GWF ,AMP(35  ))            
      CALL JVVXXX(W(1,1   ),W(1,2   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     53  ))                                                          
      CALL FVOXXX(W(1,4   ),W(1,53  ),GWF ,ZERO    ,ZERO    ,W(1,54  ))    
      CALL JIOXXX(W(1,3   ),W(1,54  ),GAD ,ZERO    ,AWIDTH  ,W(1,55  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,55  ),GAL ,AMP(36  ))            
      CALL JIOXXX(W(1,3   ),W(1,54  ),GZD ,ZMASS   ,ZWIDTH  ,W(1,56  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,56  ),GZL ,AMP(37  ))            
      CALL FVIXXX(W(1,3   ),W(1,53  ),GWF ,ZERO    ,ZERO    ,W(1,57  ))    
      CALL JIOXXX(W(1,57  ),W(1,4   ),GAU ,ZERO    ,AWIDTH  ,W(1,58  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,58  ),GAL ,AMP(38  ))            
      CALL JIOXXX(W(1,57  ),W(1,4   ),GZU ,ZMASS   ,ZWIDTH  ,W(1,59  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,59  ),GZL ,AMP(39  ))            
      CALL JW3WXX(W(1,2   ),W(1,1   ),W(1,43  ),GWWA ,GWWA ,ZERO    ,      
     &     AWIDTH  ,W(1,60  ))                                             
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,60  ),GAL ,AMP(40  ))            
      CALL JW3WXX(W(1,2   ),W(1,1   ),W(1,43  ),GWWA ,GWWZ ,ZMASS   ,      
     &     ZWIDTH  ,W(1,61  ))                                             
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,61  ),GZL ,AMP(41  ))            
      CALL JVVXXX(W(1,53  ),W(1,43  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     62  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,62  ),GAL ,AMP(42  ))            
      CALL JVVXXX(W(1,53  ),W(1,43  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     63  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,63  ),GZL ,AMP(43  ))            
      CALL FVIXXX(W(1,5   ),W(1,53  ),GWF ,ZERO    ,ZERO    ,W(1,64  ))    
      CALL IOVXXX(W(1,64  ),W(1,6   ),W(1,43  ),GWF ,AMP(44  ))            
      awp(mu,nu) = +AMP(   1)+AMP(   2)+AMP(   3)+AMP(   4)+AMP(   5)
     &             +AMP(   6)+AMP(   7)+AMP(   8)+AMP(   9)+AMP(  10)
     &             +AMP(  11)+AMP(  12)+AMP(  13)+AMP(  14)+AMP(  15)
     &             +AMP(  16)+AMP(  17)+AMP(  18)+AMP(  19)+AMP(  20)
     &             +AMP(  21)+AMP(  22)+AMP(  23)+AMP(  24)+AMP(  25)
     &             +AMP(  26)+AMP(  27)+AMP(  28)+AMP(  29)+AMP(  30)
     &             +AMP(  31)+AMP(  32)+AMP(  33)+AMP(  34)+AMP(  35)
     &             +AMP(  36)+AMP(  37)+AMP(  38)+AMP(  39)+AMP(  40)
     &             +AMP(  41)+AMP(  42)+AMP(  43)+AMP(  44)


      elseif (procid.eq.wpzhadjj .and. ztouu) then

      CALL FVIXXX(W(1,3   ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,7   ))    
      CALL JIOXXX(W(1,7   ),W(1,4   ),GZN ,ZMASS   ,ZWIDTH  ,W(1,8   ))    
      CALL FVIXXX(W(1,5   ),W(1,8   ),GZU ,ZERO    ,ZERO    ,W(1,9   ))    
      CALL IOVXXX(W(1,9   ),W(1,6   ),W(1,1   ),GAU ,AMP(1   ))            
      CALL FVIXXX(W(1,5   ),W(1,1   ),GAU ,ZERO    ,ZERO    ,W(1,10  ))    
      CALL IOVXXX(W(1,10  ),W(1,6   ),W(1,8   ),GZU ,AMP(2   ))            
      CALL FVOXXX(W(1,4   ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,11  ))    
      CALL JIOXXX(W(1,3   ),W(1,11  ),GAL ,ZERO    ,AWIDTH  ,W(1,12  ))    
      CALL FVIXXX(W(1,5   ),W(1,12  ),GAU ,ZERO    ,ZERO    ,W(1,13  ))    
      CALL IOVXXX(W(1,13  ),W(1,6   ),W(1,1   ),GAU ,AMP(3   ))            
      CALL JIOXXX(W(1,3   ),W(1,11  ),GZL ,ZMASS   ,ZWIDTH  ,W(1,14  ))    
      CALL FVIXXX(W(1,5   ),W(1,14  ),GZU ,ZERO    ,ZERO    ,W(1,15  ))    
      CALL IOVXXX(W(1,15  ),W(1,6   ),W(1,1   ),GAU ,AMP(4   ))            
      CALL IOVXXX(W(1,10  ),W(1,6   ),W(1,12  ),GAU ,AMP(5   ))            
      CALL IOVXXX(W(1,10  ),W(1,6   ),W(1,14  ),GZU ,AMP(6   ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,16  ))    
      CALL JIOXXX(W(1,16  ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,17  ))    
      CALL FVIXXX(W(1,5   ),W(1,17  ),GWF ,ZERO    ,ZERO    ,W(1,18  ))    
      CALL IOVXXX(W(1,18  ),W(1,6   ),W(1,2   ),GWF ,AMP(7   ))            
      CALL JVVXXX(W(1,2   ),W(1,17  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     19  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,19  ),GAU ,AMP(8   ))            
      CALL JVVXXX(W(1,2   ),W(1,17  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     20  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,20  ),GZU ,AMP(9   ))            
      CALL FVOXXX(W(1,11  ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,21  ))    
      CALL JIOXXX(W(1,3   ),W(1,21  ),GAL ,ZERO    ,AWIDTH  ,W(1,22  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,22  ),GAU ,AMP(10  ))            
      CALL JIOXXX(W(1,3   ),W(1,21  ),GZL ,ZMASS   ,ZWIDTH  ,W(1,23  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,23  ),GZU ,AMP(11  ))            
      CALL FVIXXX(W(1,16  ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,24  ))    
      CALL JIOXXX(W(1,24  ),W(1,4   ),GZN ,ZMASS   ,ZWIDTH  ,W(1,25  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,25  ),GZU ,AMP(12  ))            
      CALL JIOXXX(W(1,16  ),W(1,11  ),GAL ,ZERO    ,AWIDTH  ,W(1,26  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,26  ),GAU ,AMP(13  ))            
      CALL JIOXXX(W(1,16  ),W(1,11  ),GZL ,ZMASS   ,ZWIDTH  ,W(1,27  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,27  ),GZU ,AMP(14  ))            
      CALL JIOXXX(W(1,3   ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,28  ))    
      CALL JVVXXX(W(1,2   ),W(1,28  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     29  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,29  ),GAU ,ZERO    ,ZERO    ,W(1,30  ))    
      CALL IOVXXX(W(1,30  ),W(1,6   ),W(1,1   ),GAU ,AMP(15  ))            
      CALL JVVXXX(W(1,2   ),W(1,28  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     31  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,31  ),GZU ,ZERO    ,ZERO    ,W(1,32  ))    
      CALL IOVXXX(W(1,32  ),W(1,6   ),W(1,1   ),GAU ,AMP(16  ))            
      CALL FVIXXX(W(1,10  ),W(1,28  ),GWF ,ZERO    ,ZERO    ,W(1,33  ))    
      CALL IOVXXX(W(1,33  ),W(1,6   ),W(1,2   ),GWF ,AMP(17  ))            
      CALL IOVXXX(W(1,10  ),W(1,6   ),W(1,29  ),GAU ,AMP(18  ))            
      CALL IOVXXX(W(1,10  ),W(1,6   ),W(1,31  ),GZU ,AMP(19  ))            
      CALL JVVXXX(W(1,28  ),W(1,1   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     34  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,34  ),GWF ,ZERO    ,ZERO    ,W(1,35  ))    
      CALL IOVXXX(W(1,35  ),W(1,6   ),W(1,2   ),GWF ,AMP(20  ))            
      CALL JVVXXX(W(1,2   ),W(1,34  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     36  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,36  ),GAU ,AMP(21  ))            
      CALL JVVXXX(W(1,2   ),W(1,34  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     37  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,37  ),GZU ,AMP(22  ))            
      CALL JVVXXX(W(1,1   ),W(1,2   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     38  ))                                                          
      CALL FVOXXX(W(1,4   ),W(1,38  ),GWF ,ZERO    ,ZERO    ,W(1,39  ))    
      CALL JIOXXX(W(1,3   ),W(1,39  ),GAL ,ZERO    ,AWIDTH  ,W(1,40  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,40  ),GAU ,AMP(23  ))            
      CALL JIOXXX(W(1,3   ),W(1,39  ),GZL ,ZMASS   ,ZWIDTH  ,W(1,41  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,41  ),GZU ,AMP(24  ))            
      CALL FVIXXX(W(1,3   ),W(1,38  ),GWF ,ZERO    ,ZERO    ,W(1,42  ))    
      CALL JIOXXX(W(1,42  ),W(1,4   ),GZN ,ZMASS   ,ZWIDTH  ,W(1,43  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,43  ),GZU ,AMP(25  ))            
      CALL JW3WXX(W(1,2   ),W(1,1   ),W(1,28  ),GWWA ,GWWA ,ZERO    ,      
     &     AWIDTH  ,W(1,44  ))                                             
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,44  ),GAU ,AMP(26  ))            
      CALL JW3WXX(W(1,2   ),W(1,1   ),W(1,28  ),GWWA ,GWWZ ,ZMASS   ,      
     &     ZWIDTH  ,W(1,45  ))                                             
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,45  ),GZU ,AMP(27  ))            
      CALL FVIXXX(W(1,5   ),W(1,28  ),GWF ,ZERO    ,ZERO    ,W(1,46  ))    
      CALL FVIXXX(W(1,46  ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,47  ))    
      CALL IOVXXX(W(1,47  ),W(1,6   ),W(1,1   ),GAU ,AMP(28  ))            
      CALL FVIXXX(W(1,46  ),W(1,1   ),GAD ,ZERO    ,ZERO    ,W(1,48  ))    
      CALL IOVXXX(W(1,48  ),W(1,6   ),W(1,2   ),GWF ,AMP(29  ))            
      CALL JVVXXX(W(1,38  ),W(1,28  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     49  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,49  ),GAU ,AMP(30  ))            
      CALL JVVXXX(W(1,38  ),W(1,28  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     50  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,50  ),GZU ,AMP(31  ))            
      CALL IOVXXX(W(1,46  ),W(1,6   ),W(1,38  ),GWF ,AMP(32  ))            
      awp(mu,nu) = +AMP(   1)+AMP(   2)+AMP(   3)+AMP(   4)+AMP(   5)
     &             +AMP(   6)+AMP(   7)+AMP(   8)+AMP(   9)+AMP(  10)
     &             +AMP(  11)+AMP(  12)+AMP(  13)+AMP(  14)+AMP(  15)
     &             +AMP(  16)+AMP(  17)+AMP(  18)+AMP(  19)+AMP(  20)
     &             +AMP(  21)+AMP(  22)+AMP(  23)+AMP(  24)+AMP(  25)
     &             +AMP(  26)+AMP(  27)+AMP(  28)+AMP(  29)+AMP(  30)
     &             +AMP(  31)+AMP(  32)


      elseif (procid.eq.wpzhadjj .and. ztodd) then

      CALL FVIXXX(W(1,3   ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,7   ))    
      CALL JIOXXX(W(1,7   ),W(1,4   ),GZN ,ZMASS   ,ZWIDTH  ,W(1,8   ))    
      CALL FVIXXX(W(1,5   ),W(1,8   ),GZD ,ZERO    ,ZERO    ,W(1,9   ))    
      CALL IOVXXX(W(1,9   ),W(1,6   ),W(1,1   ),GAD ,AMP(1   ))            
      CALL FVIXXX(W(1,5   ),W(1,1   ),GAD ,ZERO    ,ZERO    ,W(1,10  ))    
      CALL IOVXXX(W(1,10  ),W(1,6   ),W(1,8   ),GZD ,AMP(2   ))            
      CALL FVOXXX(W(1,4   ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,11  ))    
      CALL JIOXXX(W(1,3   ),W(1,11  ),GAL ,ZERO    ,AWIDTH  ,W(1,12  ))    
      CALL FVIXXX(W(1,5   ),W(1,12  ),GAD ,ZERO    ,ZERO    ,W(1,13  ))    
      CALL IOVXXX(W(1,13  ),W(1,6   ),W(1,1   ),GAD ,AMP(3   ))            
      CALL JIOXXX(W(1,3   ),W(1,11  ),GZL ,ZMASS   ,ZWIDTH  ,W(1,14  ))    
      CALL FVIXXX(W(1,5   ),W(1,14  ),GZD ,ZERO    ,ZERO    ,W(1,15  ))    
      CALL IOVXXX(W(1,15  ),W(1,6   ),W(1,1   ),GAD ,AMP(4   ))            
      CALL IOVXXX(W(1,10  ),W(1,6   ),W(1,12  ),GAD ,AMP(5   ))            
      CALL IOVXXX(W(1,10  ),W(1,6   ),W(1,14  ),GZD ,AMP(6   ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,16  ))    
      CALL JIOXXX(W(1,16  ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,17  ))    
      CALL JVVXXX(W(1,2   ),W(1,17  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     18  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,18  ),GAD ,AMP(7   ))            
      CALL JVVXXX(W(1,2   ),W(1,17  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     19  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,19  ),GZD ,AMP(8   ))            
      CALL FVIXXX(W(1,5   ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,20  ))    
      CALL IOVXXX(W(1,20  ),W(1,6   ),W(1,17  ),GWF ,AMP(9   ))            
      CALL FVOXXX(W(1,11  ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,21  ))    
      CALL JIOXXX(W(1,3   ),W(1,21  ),GAL ,ZERO    ,AWIDTH  ,W(1,22  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,22  ),GAD ,AMP(10  ))            
      CALL JIOXXX(W(1,3   ),W(1,21  ),GZL ,ZMASS   ,ZWIDTH  ,W(1,23  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,23  ),GZD ,AMP(11  ))            
      CALL FVIXXX(W(1,16  ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,24  ))    
      CALL JIOXXX(W(1,24  ),W(1,4   ),GZN ,ZMASS   ,ZWIDTH  ,W(1,25  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,25  ),GZD ,AMP(12  ))            
      CALL JIOXXX(W(1,16  ),W(1,11  ),GAL ,ZERO    ,AWIDTH  ,W(1,26  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,26  ),GAD ,AMP(13  ))            
      CALL JIOXXX(W(1,16  ),W(1,11  ),GZL ,ZMASS   ,ZWIDTH  ,W(1,27  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,27  ),GZD ,AMP(14  ))            
      CALL JIOXXX(W(1,3   ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,28  ))    
      CALL JVVXXX(W(1,2   ),W(1,28  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     29  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,29  ),GAD ,ZERO    ,ZERO    ,W(1,30  ))    
      CALL IOVXXX(W(1,30  ),W(1,6   ),W(1,1   ),GAD ,AMP(15  ))            
      CALL JVVXXX(W(1,2   ),W(1,28  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     31  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,31  ),GZD ,ZERO    ,ZERO    ,W(1,32  ))    
      CALL IOVXXX(W(1,32  ),W(1,6   ),W(1,1   ),GAD ,AMP(16  ))            
      CALL FVIXXX(W(1,10  ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,33  ))    
      CALL IOVXXX(W(1,33  ),W(1,6   ),W(1,28  ),GWF ,AMP(17  ))            
      CALL IOVXXX(W(1,10  ),W(1,6   ),W(1,29  ),GAD ,AMP(18  ))            
      CALL IOVXXX(W(1,10  ),W(1,6   ),W(1,31  ),GZD ,AMP(19  ))            
      CALL FVIXXX(W(1,20  ),W(1,28  ),GWF ,ZERO    ,ZERO    ,W(1,34  ))    
      CALL IOVXXX(W(1,34  ),W(1,6   ),W(1,1   ),GAD ,AMP(20  ))            
      CALL JVVXXX(W(1,28  ),W(1,1   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     35  ))                                                          
      CALL JVVXXX(W(1,2   ),W(1,35  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     36  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,36  ),GAD ,AMP(21  ))            
      CALL JVVXXX(W(1,2   ),W(1,35  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     37  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,37  ),GZD ,AMP(22  ))            
      CALL FVIXXX(W(1,20  ),W(1,1   ),GAU ,ZERO    ,ZERO    ,W(1,38  ))    
      CALL IOVXXX(W(1,38  ),W(1,6   ),W(1,28  ),GWF ,AMP(23  ))            
      CALL IOVXXX(W(1,20  ),W(1,6   ),W(1,35  ),GWF ,AMP(24  ))            
      CALL JVVXXX(W(1,1   ),W(1,2   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     39  ))                                                          
      CALL FVOXXX(W(1,4   ),W(1,39  ),GWF ,ZERO    ,ZERO    ,W(1,40  ))    
      CALL JIOXXX(W(1,3   ),W(1,40  ),GAL ,ZERO    ,AWIDTH  ,W(1,41  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,41  ),GAD ,AMP(25  ))            
      CALL JIOXXX(W(1,3   ),W(1,40  ),GZL ,ZMASS   ,ZWIDTH  ,W(1,42  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,42  ),GZD ,AMP(26  ))            
      CALL FVIXXX(W(1,3   ),W(1,39  ),GWF ,ZERO    ,ZERO    ,W(1,43  ))    
      CALL JIOXXX(W(1,43  ),W(1,4   ),GZN ,ZMASS   ,ZWIDTH  ,W(1,44  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,44  ),GZD ,AMP(27  ))            
      CALL JW3WXX(W(1,2   ),W(1,1   ),W(1,28  ),GWWA ,GWWA ,ZERO    ,      
     &     AWIDTH  ,W(1,45  ))                                             
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,45  ),GAD ,AMP(28  ))            
      CALL JW3WXX(W(1,2   ),W(1,1   ),W(1,28  ),GWWA ,GWWZ ,ZMASS   ,      
     &     ZWIDTH  ,W(1,46  ))                                             
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,46  ),GZD ,AMP(29  ))            
      CALL JVVXXX(W(1,39  ),W(1,28  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     47  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,47  ),GAD ,AMP(30  ))            
      CALL JVVXXX(W(1,39  ),W(1,28  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     48  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,48  ),GZD ,AMP(31  ))            
      CALL FVIXXX(W(1,5   ),W(1,39  ),GWF ,ZERO    ,ZERO    ,W(1,49  ))    
      CALL IOVXXX(W(1,49  ),W(1,6   ),W(1,28  ),GWF ,AMP(32  ))            
      awp(mu,nu) = +AMP(   1)+AMP(   2)+AMP(   3)+AMP(   4)+AMP(   5)
     &             +AMP(   6)+AMP(   7)+AMP(   8)+AMP(   9)+AMP(  10)
     &             +AMP(  11)+AMP(  12)+AMP(  13)+AMP(  14)+AMP(  15)
     &             +AMP(  16)+AMP(  17)+AMP(  18)+AMP(  19)+AMP(  20)
     &             +AMP(  21)+AMP(  22)+AMP(  23)+AMP(  24)+AMP(  25)
     &             +AMP(  26)+AMP(  27)+AMP(  28)+AMP(  29)+AMP(  30)
     &             +AMP(  31)+AMP(  32)


      endif
                                            
     
          ENDDO
      ENDDO
c
      END
c ------------------------------------------------------------	

      subroutine vwtovwm_mg2_had(q1,q2,P,h,zwm,awm)
	
      IMPLICIT NONE

      REAL*8 Q1(0:3), Q2(0:3), P(0:3,3:6)
      complex*16 zwm(0:3,0:3), awm(0:3,0:3)
      integer lsign(3:6,2),ll,mu,nu,h
      

	! possible lepton helicities:
	data (lsign(ll,1),ll=3,6) /  -1,  -1, -1, -1/  !h=1
	data (lsign(ll,2),ll=3,6) /  -1,  -1,  1,  1/  !h=2

	! helicity combi nr.h is selected
	
	call zwtovwm1_mg2_had(q1,q2,P,lsign(3,h),zwm(0,0))
	call awtovwm1_mg2_had(q1,q2,P,lsign(3,h),awm(0,0))
		
	
	end
		
c ------------------------------------------------------------	

      subroutine zwtovwm1_mg2_had(q1,q2,P,lsign,zwm)
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : z w- -> e- ve~ mu+ mu-  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=  46) 
      INTEGER    NWAVEFUNCS     
      PARAMETER (NWAVEFUNCS=  69) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 Q1(0:3), Q2(0:3), P(0:3,3:6)
      complex*16 zwm(0:3,0:3)
      integer lsign(3:6)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J,MU, NU
      COMPLEX*16 AMP(NGRAPHS)
      COMPLEX*16 W(6,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
c  lepton spinors and Z polarization vectors
#include "tensorw.inc"
      
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/process.inc"
   
      logical lswitch
      parameter (lswitch=.true.)
C  
C ----------
C BEGIN CODE
C ----------

      if (lswitch) then
      	do i = 1,6
         w(i,3) = wl(i)
         w(i,4) = wn(i)
         w(i,5) = lp(i)
         w(i,6) = lm(i)

	 w(i,14) = au(i)
	 w(i,15) = zu(i)
	 w(i,36) = we(i)
      	enddo
      else 
      	CALL OXXXXX(P(0,3   ),ZERO ,+lsign(3   ),+1,W(1,3   ))	   
      	CALL IXXXXX(P(0,4   ),ZERO ,-lsign(4   ),-1,W(1,4   ))	   
      	CALL IXXXXX(P(0,5   ),ZERO ,-lsign(5   ),-1,W(1,5   ))	   
      	CALL OXXXXX(P(0,6   ),ZERO ,+lsign(6   ),+1,W(1,6   ))	   
      	
	CALL JIOXXX(W(1,5   ),W(1,6   ),GAL ,ZERO    ,ZERO    ,W(1,14  ))    
        CALL JIOXXX(W(1,5   ),W(1,6   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,15  ))    
        CALL JIOXXX(W(1,4   ),W(1,3   ),GWF ,WMASS   ,WWIDTH  ,W(1,36  ))    
      endif
	

      DO MU = 0,3
         DO NU = 0,3
      
      CALL VCARTX(Q1(0),ZMASS ,zwidth, MU,-1,W(1,1   ))        
      CALL VCARTX(Q2(0),WMASS ,wwidth, NU,-1,W(1,2   )) 
 
      if (procid.eq.wmhadzjj) then

      CALL FVOXXX(W(1,3   ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,7   ))    
      CALL JIOXXX(W(1,4   ),W(1,7   ),GAU ,ZERO    ,AWIDTH  ,W(1,8   ))    
      CALL FVIXXX(W(1,5   ),W(1,8   ),GAL ,ZERO    ,ZERO    ,W(1,9   ))    
      CALL IOVXXX(W(1,9   ),W(1,6   ),W(1,1   ),GZL ,AMP(1   ))            
      CALL JIOXXX(W(1,4   ),W(1,7   ),GZU ,ZMASS   ,ZWIDTH  ,W(1,10  ))    
      CALL FVIXXX(W(1,5   ),W(1,10  ),GZL ,ZERO    ,ZERO    ,W(1,11  ))    
      CALL IOVXXX(W(1,11  ),W(1,6   ),W(1,1   ),GZL ,AMP(2   ))            
      CALL FVIXXX(W(1,5   ),W(1,1   ),GZL ,ZERO    ,ZERO    ,W(1,12  ))    
      CALL IOVXXX(W(1,12  ),W(1,6   ),W(1,8   ),GAL ,AMP(3   ))            
      CALL IOVXXX(W(1,12  ),W(1,6   ),W(1,10  ),GZL ,AMP(4   ))            
      CALL FVIXXX(W(1,4   ),W(1,1   ),GZU ,ZERO    ,ZERO    ,W(1,13  ))    
      CALL JIOXXX(W(1,13  ),W(1,3   ),GWF ,WMASS   ,WWIDTH  ,W(1,14  ))    
      CALL FVIXXX(W(1,5   ),W(1,14  ),GWF ,ZERO    ,ZERO    ,W(1,15  ))    
      CALL IOVXXX(W(1,15  ),W(1,6   ),W(1,2   ),GWF ,AMP(5   ))            
      CALL JVVXXX(W(1,14  ),W(1,2   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     16  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,16  ),GAL ,AMP(6   ))            
      CALL JVVXXX(W(1,14  ),W(1,2   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     17  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,17  ),GZL ,AMP(7   ))            
      CALL FVIXXX(W(1,13  ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,18  ))    
      CALL JIOXXX(W(1,18  ),W(1,3   ),GAD ,ZERO    ,AWIDTH  ,W(1,19  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,19  ),GAL ,AMP(8   ))            
      CALL JIOXXX(W(1,18  ),W(1,3   ),GZD ,ZMASS   ,ZWIDTH  ,W(1,20  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,20  ),GZL ,AMP(9   ))            
      CALL FVOXXX(W(1,7   ),W(1,1   ),GZU ,ZERO    ,ZERO    ,W(1,21  ))    
      CALL JIOXXX(W(1,4   ),W(1,21  ),GAU ,ZERO    ,AWIDTH  ,W(1,22  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,22  ),GAL ,AMP(10  ))            
      CALL JIOXXX(W(1,4   ),W(1,21  ),GZU ,ZMASS   ,ZWIDTH  ,W(1,23  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,23  ),GZL ,AMP(11  ))            
      CALL JIOXXX(W(1,13  ),W(1,7   ),GAU ,ZERO    ,AWIDTH  ,W(1,24  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,24  ),GAL ,AMP(12  ))            
      CALL JIOXXX(W(1,13  ),W(1,7   ),GZU ,ZMASS   ,ZWIDTH  ,W(1,25  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,25  ),GZL ,AMP(13  ))            
      CALL FVIXXX(W(1,4   ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,26  ))    
      CALL JIOXXX(W(1,26  ),W(1,3   ),GAD ,ZERO    ,AWIDTH  ,W(1,27  ))    
      CALL FVIXXX(W(1,5   ),W(1,27  ),GAL ,ZERO    ,ZERO    ,W(1,28  ))    
      CALL IOVXXX(W(1,28  ),W(1,6   ),W(1,1   ),GZL ,AMP(14  ))            
      CALL JIOXXX(W(1,26  ),W(1,3   ),GZD ,ZMASS   ,ZWIDTH  ,W(1,29  ))    
      CALL FVIXXX(W(1,5   ),W(1,29  ),GZL ,ZERO    ,ZERO    ,W(1,30  ))    
      CALL IOVXXX(W(1,30  ),W(1,6   ),W(1,1   ),GZL ,AMP(15  ))            
      CALL IOVXXX(W(1,12  ),W(1,6   ),W(1,27  ),GAL ,AMP(16  ))            
      CALL IOVXXX(W(1,12  ),W(1,6   ),W(1,29  ),GZL ,AMP(17  ))            
      CALL FVOXXX(W(1,3   ),W(1,1   ),GZD ,ZERO    ,ZERO    ,W(1,31  ))    
      CALL JIOXXX(W(1,4   ),W(1,31  ),GWF ,WMASS   ,WWIDTH  ,W(1,32  ))    
      CALL FVIXXX(W(1,5   ),W(1,32  ),GWF ,ZERO    ,ZERO    ,W(1,33  ))    
      CALL IOVXXX(W(1,33  ),W(1,6   ),W(1,2   ),GWF ,AMP(18  ))            
      CALL JVVXXX(W(1,32  ),W(1,2   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     34  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,34  ),GAL ,AMP(19  ))            
      CALL JVVXXX(W(1,32  ),W(1,2   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     35  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,35  ),GZL ,AMP(20  ))            
      CALL FVIXXX(W(1,26  ),W(1,1   ),GZD ,ZERO    ,ZERO    ,W(1,36  ))    
      CALL JIOXXX(W(1,36  ),W(1,3   ),GAD ,ZERO    ,AWIDTH  ,W(1,37  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,37  ),GAL ,AMP(21  ))            
      CALL JIOXXX(W(1,36  ),W(1,3   ),GZD ,ZMASS   ,ZWIDTH  ,W(1,38  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,38  ),GZL ,AMP(22  ))            
      CALL FVOXXX(W(1,31  ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,39  ))    
      CALL JIOXXX(W(1,4   ),W(1,39  ),GAU ,ZERO    ,AWIDTH  ,W(1,40  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,40  ),GAL ,AMP(23  ))            
      CALL JIOXXX(W(1,4   ),W(1,39  ),GZU ,ZMASS   ,ZWIDTH  ,W(1,41  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,41  ),GZL ,AMP(24  ))            
      CALL JIOXXX(W(1,26  ),W(1,31  ),GAD ,ZERO    ,AWIDTH  ,W(1,42  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,42  ),GAL ,AMP(25  ))            
      CALL JIOXXX(W(1,26  ),W(1,31  ),GZD ,ZMASS   ,ZWIDTH  ,W(1,43  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,43  ),GZL ,AMP(26  ))            
      CALL JIOXXX(W(1,4   ),W(1,3   ),GWF ,WMASS   ,WWIDTH  ,W(1,44  ))    
      CALL JVVXXX(W(1,44  ),W(1,2   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     45  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,45  ),GAL ,ZERO    ,ZERO    ,W(1,46  ))    
      CALL IOVXXX(W(1,46  ),W(1,6   ),W(1,1   ),GZL ,AMP(27  ))            
      CALL JVVXXX(W(1,44  ),W(1,2   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     47  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,47  ),GZL ,ZERO    ,ZERO    ,W(1,48  ))    
      CALL IOVXXX(W(1,48  ),W(1,6   ),W(1,1   ),GZL ,AMP(28  ))            
      CALL FVIXXX(W(1,12  ),W(1,44  ),GWF ,ZERO    ,ZERO    ,W(1,49  ))    
      CALL IOVXXX(W(1,49  ),W(1,6   ),W(1,2   ),GWF ,AMP(29  ))            
      CALL HVVXXX(W(1,44  ),W(1,2   ),GWWH ,HMASS   ,HWIDTH  ,W(1,         
     &     50  ))                                                          
      CALL JVSXXX(W(1,1   ),W(1,50  ),GZZH ,ZMASS   ,ZWIDTH  ,W(1,         
     &     51  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,51  ),GZL ,AMP(30  ))            
      CALL IOVXXX(W(1,12  ),W(1,6   ),W(1,45  ),GAL ,AMP(31  ))            
      CALL IOVXXX(W(1,12  ),W(1,6   ),W(1,47  ),GZL ,AMP(32  ))            
      CALL JVVXXX(W(1,1   ),W(1,44  ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     52  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,52  ),GWF ,ZERO    ,ZERO    ,W(1,53  ))    
      CALL IOVXXX(W(1,53  ),W(1,6   ),W(1,2   ),GWF ,AMP(33  ))            
      CALL JVVXXX(W(1,52  ),W(1,2   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     54  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,54  ),GAL ,AMP(34  ))            
      CALL JVVXXX(W(1,52  ),W(1,2   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     55  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,55  ),GZL ,AMP(35  ))            
      CALL JVVXXX(W(1,2   ),W(1,1   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     56  ))                                                          
      CALL FVIXXX(W(1,4   ),W(1,56  ),GWF ,ZERO    ,ZERO    ,W(1,57  ))    
      CALL JIOXXX(W(1,57  ),W(1,3   ),GAD ,ZERO    ,AWIDTH  ,W(1,58  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,58  ),GAL ,AMP(36  ))            
      CALL JIOXXX(W(1,57  ),W(1,3   ),GZD ,ZMASS   ,ZWIDTH  ,W(1,59  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,59  ),GZL ,AMP(37  ))            
      CALL FVOXXX(W(1,3   ),W(1,56  ),GWF ,ZERO    ,ZERO    ,W(1,60  ))    
      CALL JIOXXX(W(1,4   ),W(1,60  ),GAU ,ZERO    ,AWIDTH  ,W(1,61  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,61  ),GAL ,AMP(38  ))            
      CALL JIOXXX(W(1,4   ),W(1,60  ),GZU ,ZMASS   ,ZWIDTH  ,W(1,62  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,62  ),GZL ,AMP(39  ))            
      CALL JW3WXX(W(1,44  ),W(1,1   ),W(1,2   ),GWWZ ,GWWA ,ZERO    ,      
     &     AWIDTH  ,W(1,63  ))                                             
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,63  ),GAL ,AMP(40  ))            
      CALL JW3WXX(W(1,44  ),W(1,1   ),W(1,2   ),GWWZ ,GWWZ ,ZMASS   ,      
     &     ZWIDTH  ,W(1,64  ))                                             
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,64  ),GZL ,AMP(41  ))            
      CALL FVIXXX(W(1,5   ),W(1,44  ),GWF ,ZERO    ,ZERO    ,W(1,65  ))    
      CALL FVIXXX(W(1,65  ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,66  ))    
      CALL IOVXXX(W(1,66  ),W(1,6   ),W(1,1   ),GZL ,AMP(42  ))            
      CALL FVIXXX(W(1,65  ),W(1,1   ),GZN ,ZERO    ,ZERO    ,W(1,67  ))    
      CALL IOVXXX(W(1,67  ),W(1,6   ),W(1,2   ),GWF ,AMP(43  ))            
      CALL JVVXXX(W(1,44  ),W(1,56  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     68  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,68  ),GAL ,AMP(44  ))            
      CALL JVVXXX(W(1,44  ),W(1,56  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     69  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,69  ),GZL ,AMP(45  ))            
      CALL IOVXXX(W(1,65  ),W(1,6   ),W(1,56  ),GWF ,AMP(46  ))            
      zwm(mu,nu) = -AMP(   1)-AMP(   2)-AMP(   3)-AMP(   4)-AMP(   5)
     &             -AMP(   6)-AMP(   7)-AMP(   8)-AMP(   9)-AMP(  10)
     &             -AMP(  11)-AMP(  12)-AMP(  13)-AMP(  14)-AMP(  15)
     &             -AMP(  16)-AMP(  17)-AMP(  18)-AMP(  19)-AMP(  20)
     &             -AMP(  21)-AMP(  22)-AMP(  23)-AMP(  24)-AMP(  25)
     &             -AMP(  26)-AMP(  27)-AMP(  28)-AMP(  29)-AMP(  30)
     &             -AMP(  31)-AMP(  32)-AMP(  33)-AMP(  34)-AMP(  35)
     &             -AMP(  36)-AMP(  37)-AMP(  38)-AMP(  39)-AMP(  40)
     &             -AMP(  41)-AMP(  42)-AMP(  43)-AMP(  44)-AMP(  45)
     &             -AMP(  46)


      elseif (procid.eq.wmzhadjj .and. ztouu) then

      CALL FVOXXX(W(1,3   ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,7   ))    
      CALL JIOXXX(W(1,4   ),W(1,7   ),GZN ,ZMASS   ,ZWIDTH  ,W(1,8   ))    
      CALL FVIXXX(W(1,5   ),W(1,8   ),GZU ,ZERO    ,ZERO    ,W(1,9   ))    
      CALL IOVXXX(W(1,9   ),W(1,6   ),W(1,1   ),GZU ,AMP(1   ))            
      CALL FVIXXX(W(1,5   ),W(1,1   ),GZU ,ZERO    ,ZERO    ,W(1,10  ))    
      CALL IOVXXX(W(1,10  ),W(1,6   ),W(1,8   ),GZU ,AMP(2   ))            
      CALL FVIXXX(W(1,4   ),W(1,1   ),GZN ,ZERO    ,ZERO    ,W(1,11  ))    
      CALL JIOXXX(W(1,11  ),W(1,3   ),GWF ,WMASS   ,WWIDTH  ,W(1,12  ))    
      CALL JVVXXX(W(1,12  ),W(1,2   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     13  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,13  ),GAU ,AMP(3   ))            
      CALL JVVXXX(W(1,12  ),W(1,2   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     14  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,14  ),GZU ,AMP(4   ))            
      CALL FVIXXX(W(1,5   ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,15  ))    
      CALL IOVXXX(W(1,15  ),W(1,6   ),W(1,12  ),GWF ,AMP(5   ))            
      CALL FVIXXX(W(1,11  ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,16  ))    
      CALL JIOXXX(W(1,16  ),W(1,3   ),GAL ,ZERO    ,AWIDTH  ,W(1,17  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,17  ),GAU ,AMP(6   ))            
      CALL JIOXXX(W(1,16  ),W(1,3   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,18  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,18  ),GZU ,AMP(7   ))            
      CALL FVOXXX(W(1,7   ),W(1,1   ),GZN ,ZERO    ,ZERO    ,W(1,19  ))    
      CALL JIOXXX(W(1,4   ),W(1,19  ),GZN ,ZMASS   ,ZWIDTH  ,W(1,20  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,20  ),GZU ,AMP(8   ))            
      CALL JIOXXX(W(1,11  ),W(1,7   ),GZN ,ZMASS   ,ZWIDTH  ,W(1,21  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,21  ),GZU ,AMP(9   ))            
      CALL FVIXXX(W(1,4   ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,22  ))    
      CALL JIOXXX(W(1,22  ),W(1,3   ),GAL ,ZERO    ,AWIDTH  ,W(1,23  ))    
      CALL FVIXXX(W(1,5   ),W(1,23  ),GAU ,ZERO    ,ZERO    ,W(1,24  ))    
      CALL IOVXXX(W(1,24  ),W(1,6   ),W(1,1   ),GZU ,AMP(10  ))            
      CALL JIOXXX(W(1,22  ),W(1,3   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,25  ))    
      CALL FVIXXX(W(1,5   ),W(1,25  ),GZU ,ZERO    ,ZERO    ,W(1,26  ))    
      CALL IOVXXX(W(1,26  ),W(1,6   ),W(1,1   ),GZU ,AMP(11  ))            
      CALL IOVXXX(W(1,10  ),W(1,6   ),W(1,23  ),GAU ,AMP(12  ))            
      CALL IOVXXX(W(1,10  ),W(1,6   ),W(1,25  ),GZU ,AMP(13  ))            
      CALL FVOXXX(W(1,3   ),W(1,1   ),GZL ,ZERO    ,ZERO    ,W(1,27  ))    
      CALL JIOXXX(W(1,4   ),W(1,27  ),GWF ,WMASS   ,WWIDTH  ,W(1,28  ))    
      CALL JVVXXX(W(1,28  ),W(1,2   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     29  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,29  ),GAU ,AMP(14  ))            
      CALL JVVXXX(W(1,28  ),W(1,2   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     30  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,30  ),GZU ,AMP(15  ))            
      CALL IOVXXX(W(1,15  ),W(1,6   ),W(1,28  ),GWF ,AMP(16  ))            
      CALL FVIXXX(W(1,22  ),W(1,1   ),GZL ,ZERO    ,ZERO    ,W(1,31  ))    
      CALL JIOXXX(W(1,31  ),W(1,3   ),GAL ,ZERO    ,AWIDTH  ,W(1,32  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,32  ),GAU ,AMP(17  ))            
      CALL JIOXXX(W(1,31  ),W(1,3   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,33  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,33  ),GZU ,AMP(18  ))            
      CALL FVOXXX(W(1,27  ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,34  ))    
      CALL JIOXXX(W(1,4   ),W(1,34  ),GZN ,ZMASS   ,ZWIDTH  ,W(1,35  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,35  ),GZU ,AMP(19  ))            
      CALL JIOXXX(W(1,22  ),W(1,27  ),GAL ,ZERO    ,AWIDTH  ,W(1,36  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,36  ),GAU ,AMP(20  ))            
      CALL JIOXXX(W(1,22  ),W(1,27  ),GZL ,ZMASS   ,ZWIDTH  ,W(1,37  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,37  ),GZU ,AMP(21  ))            
      CALL JIOXXX(W(1,4   ),W(1,3   ),GWF ,WMASS   ,WWIDTH  ,W(1,38  ))    
      CALL JVVXXX(W(1,38  ),W(1,2   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     39  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,39  ),GAU ,ZERO    ,ZERO    ,W(1,40  ))    
      CALL IOVXXX(W(1,40  ),W(1,6   ),W(1,1   ),GZU ,AMP(22  ))            
      CALL JVVXXX(W(1,38  ),W(1,2   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     41  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,41  ),GZU ,ZERO    ,ZERO    ,W(1,42  ))    
      CALL IOVXXX(W(1,42  ),W(1,6   ),W(1,1   ),GZU ,AMP(23  ))            
      CALL HVVXXX(W(1,38  ),W(1,2   ),GWWH ,HMASS   ,HWIDTH  ,W(1,         
     &     43  ))                                                          
      CALL JVSXXX(W(1,1   ),W(1,43  ),GZZH ,ZMASS   ,ZWIDTH  ,W(1,         
     &     44  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,44  ),GZU ,AMP(24  ))            
      CALL FVIXXX(W(1,10  ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,45  ))    
      CALL IOVXXX(W(1,45  ),W(1,6   ),W(1,38  ),GWF ,AMP(25  ))            
      CALL IOVXXX(W(1,10  ),W(1,6   ),W(1,39  ),GAU ,AMP(26  ))            
      CALL IOVXXX(W(1,10  ),W(1,6   ),W(1,41  ),GZU ,AMP(27  ))            
      CALL FVIXXX(W(1,15  ),W(1,38  ),GWF ,ZERO    ,ZERO    ,W(1,46  ))    
      CALL IOVXXX(W(1,46  ),W(1,6   ),W(1,1   ),GZU ,AMP(28  ))            
      CALL JVVXXX(W(1,1   ),W(1,38  ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     47  ))                                                          
      CALL JVVXXX(W(1,47  ),W(1,2   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     48  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,48  ),GAU ,AMP(29  ))            
      CALL JVVXXX(W(1,47  ),W(1,2   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     49  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,49  ),GZU ,AMP(30  ))            
      CALL FVIXXX(W(1,15  ),W(1,1   ),GZD ,ZERO    ,ZERO    ,W(1,50  ))    
      CALL IOVXXX(W(1,50  ),W(1,6   ),W(1,38  ),GWF ,AMP(31  ))            
      CALL IOVXXX(W(1,15  ),W(1,6   ),W(1,47  ),GWF ,AMP(32  ))            
      CALL JVVXXX(W(1,2   ),W(1,1   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     51  ))                                                          
      CALL FVIXXX(W(1,4   ),W(1,51  ),GWF ,ZERO    ,ZERO    ,W(1,52  ))    
      CALL JIOXXX(W(1,52  ),W(1,3   ),GAL ,ZERO    ,AWIDTH  ,W(1,53  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,53  ),GAU ,AMP(33  ))            
      CALL JIOXXX(W(1,52  ),W(1,3   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,54  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,54  ),GZU ,AMP(34  ))            
      CALL FVOXXX(W(1,3   ),W(1,51  ),GWF ,ZERO    ,ZERO    ,W(1,55  ))    
      CALL JIOXXX(W(1,4   ),W(1,55  ),GZN ,ZMASS   ,ZWIDTH  ,W(1,56  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,56  ),GZU ,AMP(35  ))            
      CALL JW3WXX(W(1,38  ),W(1,1   ),W(1,2   ),GWWZ ,GWWA ,ZERO    ,      
     &     AWIDTH  ,W(1,57  ))                                             
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,57  ),GAU ,AMP(36  ))            
      CALL JW3WXX(W(1,38  ),W(1,1   ),W(1,2   ),GWWZ ,GWWZ ,ZMASS   ,      
     &     ZWIDTH  ,W(1,58  ))                                             
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,58  ),GZU ,AMP(37  ))            
      CALL JVVXXX(W(1,38  ),W(1,51  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     59  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,59  ),GAU ,AMP(38  ))            
      CALL JVVXXX(W(1,38  ),W(1,51  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     60  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,60  ),GZU ,AMP(39  ))            
      CALL FVIXXX(W(1,5   ),W(1,51  ),GWF ,ZERO    ,ZERO    ,W(1,61  ))    
      CALL IOVXXX(W(1,61  ),W(1,6   ),W(1,38  ),GWF ,AMP(40  ))            
      zwm(mu,nu) = -AMP(   1)-AMP(   2)-AMP(   3)-AMP(   4)-AMP(   5)
     &             -AMP(   6)-AMP(   7)-AMP(   8)-AMP(   9)-AMP(  10)
     &             -AMP(  11)-AMP(  12)-AMP(  13)-AMP(  14)-AMP(  15)
     &             -AMP(  16)-AMP(  17)-AMP(  18)-AMP(  19)-AMP(  20)
     &             -AMP(  21)-AMP(  22)-AMP(  23)-AMP(  24)-AMP(  25)
     &             -AMP(  26)-AMP(  27)-AMP(  28)-AMP(  29)-AMP(  30)
     &             -AMP(  31)-AMP(  32)-AMP(  33)-AMP(  34)-AMP(  35)
     &             -AMP(  36)-AMP(  37)-AMP(  38)-AMP(  39)-AMP(  40)


      elseif (procid.eq.wmzhadjj .and. ztodd) then

      CALL FVOXXX(W(1,3   ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,7   ))    
      CALL JIOXXX(W(1,4   ),W(1,7   ),GZN ,ZMASS   ,ZWIDTH  ,W(1,8   ))    
      CALL FVIXXX(W(1,5   ),W(1,8   ),GZD ,ZERO    ,ZERO    ,W(1,9   ))    
      CALL IOVXXX(W(1,9   ),W(1,6   ),W(1,1   ),GZD ,AMP(1   ))            
      CALL FVIXXX(W(1,5   ),W(1,1   ),GZD ,ZERO    ,ZERO    ,W(1,10  ))    
      CALL IOVXXX(W(1,10  ),W(1,6   ),W(1,8   ),GZD ,AMP(2   ))            
      CALL FVIXXX(W(1,4   ),W(1,1   ),GZN ,ZERO    ,ZERO    ,W(1,11  ))    
      CALL JIOXXX(W(1,11  ),W(1,3   ),GWF ,WMASS   ,WWIDTH  ,W(1,12  ))    
      CALL FVIXXX(W(1,5   ),W(1,12  ),GWF ,ZERO    ,ZERO    ,W(1,13  ))    
      CALL IOVXXX(W(1,13  ),W(1,6   ),W(1,2   ),GWF ,AMP(3   ))            
      CALL JVVXXX(W(1,12  ),W(1,2   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     14  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,14  ),GAD ,AMP(4   ))            
      CALL JVVXXX(W(1,12  ),W(1,2   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     15  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,15  ),GZD ,AMP(5   ))            
      CALL FVIXXX(W(1,11  ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,16  ))    
      CALL JIOXXX(W(1,16  ),W(1,3   ),GAL ,ZERO    ,AWIDTH  ,W(1,17  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,17  ),GAD ,AMP(6   ))            
      CALL JIOXXX(W(1,16  ),W(1,3   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,18  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,18  ),GZD ,AMP(7   ))            
      CALL FVOXXX(W(1,7   ),W(1,1   ),GZN ,ZERO    ,ZERO    ,W(1,19  ))    
      CALL JIOXXX(W(1,4   ),W(1,19  ),GZN ,ZMASS   ,ZWIDTH  ,W(1,20  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,20  ),GZD ,AMP(8   ))            
      CALL JIOXXX(W(1,11  ),W(1,7   ),GZN ,ZMASS   ,ZWIDTH  ,W(1,21  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,21  ),GZD ,AMP(9   ))            
      CALL FVIXXX(W(1,4   ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,22  ))    
      CALL JIOXXX(W(1,22  ),W(1,3   ),GAL ,ZERO    ,AWIDTH  ,W(1,23  ))    
      CALL FVIXXX(W(1,5   ),W(1,23  ),GAD ,ZERO    ,ZERO    ,W(1,24  ))    
      CALL IOVXXX(W(1,24  ),W(1,6   ),W(1,1   ),GZD ,AMP(10  ))            
      CALL JIOXXX(W(1,22  ),W(1,3   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,25  ))    
      CALL FVIXXX(W(1,5   ),W(1,25  ),GZD ,ZERO    ,ZERO    ,W(1,26  ))    
      CALL IOVXXX(W(1,26  ),W(1,6   ),W(1,1   ),GZD ,AMP(11  ))            
      CALL IOVXXX(W(1,10  ),W(1,6   ),W(1,23  ),GAD ,AMP(12  ))            
      CALL IOVXXX(W(1,10  ),W(1,6   ),W(1,25  ),GZD ,AMP(13  ))            
      CALL FVOXXX(W(1,3   ),W(1,1   ),GZL ,ZERO    ,ZERO    ,W(1,27  ))    
      CALL JIOXXX(W(1,4   ),W(1,27  ),GWF ,WMASS   ,WWIDTH  ,W(1,28  ))    
      CALL FVIXXX(W(1,5   ),W(1,28  ),GWF ,ZERO    ,ZERO    ,W(1,29  ))    
      CALL IOVXXX(W(1,29  ),W(1,6   ),W(1,2   ),GWF ,AMP(14  ))            
      CALL JVVXXX(W(1,28  ),W(1,2   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     30  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,30  ),GAD ,AMP(15  ))            
      CALL JVVXXX(W(1,28  ),W(1,2   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     31  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,31  ),GZD ,AMP(16  ))            
      CALL FVIXXX(W(1,22  ),W(1,1   ),GZL ,ZERO    ,ZERO    ,W(1,32  ))    
      CALL JIOXXX(W(1,32  ),W(1,3   ),GAL ,ZERO    ,AWIDTH  ,W(1,33  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,33  ),GAD ,AMP(17  ))            
      CALL JIOXXX(W(1,32  ),W(1,3   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,34  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,34  ),GZD ,AMP(18  ))            
      CALL FVOXXX(W(1,27  ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,35  ))    
      CALL JIOXXX(W(1,4   ),W(1,35  ),GZN ,ZMASS   ,ZWIDTH  ,W(1,36  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,36  ),GZD ,AMP(19  ))            
      CALL JIOXXX(W(1,22  ),W(1,27  ),GAL ,ZERO    ,AWIDTH  ,W(1,37  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,37  ),GAD ,AMP(20  ))            
      CALL JIOXXX(W(1,22  ),W(1,27  ),GZL ,ZMASS   ,ZWIDTH  ,W(1,38  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,38  ),GZD ,AMP(21  ))            
      CALL JIOXXX(W(1,4   ),W(1,3   ),GWF ,WMASS   ,WWIDTH  ,W(1,39  ))    
      CALL JVVXXX(W(1,39  ),W(1,2   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     40  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,40  ),GAD ,ZERO    ,ZERO    ,W(1,41  ))    
      CALL IOVXXX(W(1,41  ),W(1,6   ),W(1,1   ),GZD ,AMP(22  ))            
      CALL JVVXXX(W(1,39  ),W(1,2   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     42  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,42  ),GZD ,ZERO    ,ZERO    ,W(1,43  ))    
      CALL IOVXXX(W(1,43  ),W(1,6   ),W(1,1   ),GZD ,AMP(23  ))            
      CALL FVIXXX(W(1,10  ),W(1,39  ),GWF ,ZERO    ,ZERO    ,W(1,44  ))    
      CALL IOVXXX(W(1,44  ),W(1,6   ),W(1,2   ),GWF ,AMP(24  ))            
      CALL HVVXXX(W(1,39  ),W(1,2   ),GWWH ,HMASS   ,HWIDTH  ,W(1,         
     &     45  ))                                                          
      CALL JVSXXX(W(1,1   ),W(1,45  ),GZZH ,ZMASS   ,ZWIDTH  ,W(1,         
     &     46  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,46  ),GZD ,AMP(25  ))            
      CALL IOVXXX(W(1,10  ),W(1,6   ),W(1,40  ),GAD ,AMP(26  ))            
      CALL IOVXXX(W(1,10  ),W(1,6   ),W(1,42  ),GZD ,AMP(27  ))            
      CALL JVVXXX(W(1,1   ),W(1,39  ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     47  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,47  ),GWF ,ZERO    ,ZERO    ,W(1,48  ))    
      CALL IOVXXX(W(1,48  ),W(1,6   ),W(1,2   ),GWF ,AMP(28  ))            
      CALL JVVXXX(W(1,47  ),W(1,2   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     49  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,49  ),GAD ,AMP(29  ))            
      CALL JVVXXX(W(1,47  ),W(1,2   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     50  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,50  ),GZD ,AMP(30  ))            
      CALL JVVXXX(W(1,2   ),W(1,1   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     51  ))                                                          
      CALL FVIXXX(W(1,4   ),W(1,51  ),GWF ,ZERO    ,ZERO    ,W(1,52  ))    
      CALL JIOXXX(W(1,52  ),W(1,3   ),GAL ,ZERO    ,AWIDTH  ,W(1,53  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,53  ),GAD ,AMP(31  ))            
      CALL JIOXXX(W(1,52  ),W(1,3   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,54  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,54  ),GZD ,AMP(32  ))            
      CALL FVOXXX(W(1,3   ),W(1,51  ),GWF ,ZERO    ,ZERO    ,W(1,55  ))    
      CALL JIOXXX(W(1,4   ),W(1,55  ),GZN ,ZMASS   ,ZWIDTH  ,W(1,56  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,56  ),GZD ,AMP(33  ))            
      CALL JW3WXX(W(1,39  ),W(1,1   ),W(1,2   ),GWWZ ,GWWA ,ZERO    ,      
     &     AWIDTH  ,W(1,57  ))                                             
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,57  ),GAD ,AMP(34  ))            
      CALL JW3WXX(W(1,39  ),W(1,1   ),W(1,2   ),GWWZ ,GWWZ ,ZMASS   ,      
     &     ZWIDTH  ,W(1,58  ))                                             
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,58  ),GZD ,AMP(35  ))            
      CALL FVIXXX(W(1,5   ),W(1,39  ),GWF ,ZERO    ,ZERO    ,W(1,59  ))    
      CALL FVIXXX(W(1,59  ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,60  ))    
      CALL IOVXXX(W(1,60  ),W(1,6   ),W(1,1   ),GZD ,AMP(36  ))            
      CALL FVIXXX(W(1,59  ),W(1,1   ),GZU ,ZERO    ,ZERO    ,W(1,61  ))    
      CALL IOVXXX(W(1,61  ),W(1,6   ),W(1,2   ),GWF ,AMP(37  ))            
      CALL JVVXXX(W(1,39  ),W(1,51  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     62  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,62  ),GAD ,AMP(38  ))            
      CALL JVVXXX(W(1,39  ),W(1,51  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     63  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,63  ),GZD ,AMP(39  ))            
      CALL IOVXXX(W(1,59  ),W(1,6   ),W(1,51  ),GWF ,AMP(40  ))            
      zwm(mu,nu) = -AMP(   1)-AMP(   2)-AMP(   3)-AMP(   4)-AMP(   5)
     &             -AMP(   6)-AMP(   7)-AMP(   8)-AMP(   9)-AMP(  10)
     &             -AMP(  11)-AMP(  12)-AMP(  13)-AMP(  14)-AMP(  15)
     &             -AMP(  16)-AMP(  17)-AMP(  18)-AMP(  19)-AMP(  20)
     &             -AMP(  21)-AMP(  22)-AMP(  23)-AMP(  24)-AMP(  25)
     &             -AMP(  26)-AMP(  27)-AMP(  28)-AMP(  29)-AMP(  30)
     &             -AMP(  31)-AMP(  32)-AMP(  33)-AMP(  34)-AMP(  35)
     &             -AMP(  36)-AMP(  37)-AMP(  38)-AMP(  39)-AMP(  40)


      endif
             
	   ENDDO
       ENDDO

     
      END
      
c ------------------------------------------------------------	

      subroutine awtovwm1_mg2_had(q1,q2,P,lsign,awm)
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : a w- -> e- ve~ mu+ mu-  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=  44) 
      INTEGER    NWAVEFUNCS     
      PARAMETER (NWAVEFUNCS=  66) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 Q1(0:3), Q2(0:3), P(0:3,3:6)
      complex*16 awm(0:3,0:3)
      integer lsign(3:6)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J,MU, NU
      COMPLEX*16 AMP(NGRAPHS)
      COMPLEX*16 W(6,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/process.inc"
   
      logical lswitch
      parameter (lswitch=.true.)
c
C ----------
C BEGIN CODE
C ----------

      CALL OXXXXX(P(0,3   ),ZERO ,+lsign(3   ),+1,W(1,3   ))        
      CALL IXXXXX(P(0,4   ),ZERO ,-lsign(4   ),-1,W(1,4   ))        
      CALL IXXXXX(P(0,5   ),ZERO ,-lsign(5   ),-1,W(1,5   ))        
      CALL OXXXXX(P(0,6   ),ZERO ,+lsign(6   ),+1,W(1,6   ))        

      DO MU = 0,3
         DO NU = 0,3
      
      CALL VCARTX(Q1(0),ZERO,zero, MU,-1,W(1,1   ))        
      CALL VCARTX(Q2(0),WMASS ,wwidth, NU,-1,W(1,2   )) 

      if (procid.eq.wmhadzjj) then

      CALL FVOXXX(W(1,3   ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,7   ))    
      CALL JIOXXX(W(1,4   ),W(1,7   ),GAU ,ZERO    ,AWIDTH  ,W(1,8   ))    
      CALL FVIXXX(W(1,5   ),W(1,8   ),GAL ,ZERO    ,ZERO    ,W(1,9   ))    
      CALL IOVXXX(W(1,9   ),W(1,6   ),W(1,1   ),GAL ,AMP(1   ))            
      CALL JIOXXX(W(1,4   ),W(1,7   ),GZU ,ZMASS   ,ZWIDTH  ,W(1,10  ))    
      CALL FVIXXX(W(1,5   ),W(1,10  ),GZL ,ZERO    ,ZERO    ,W(1,11  ))    
      CALL IOVXXX(W(1,11  ),W(1,6   ),W(1,1   ),GAL ,AMP(2   ))            
      CALL FVIXXX(W(1,5   ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,12  ))    
      CALL IOVXXX(W(1,12  ),W(1,6   ),W(1,8   ),GAL ,AMP(3   ))            
      CALL IOVXXX(W(1,12  ),W(1,6   ),W(1,10  ),GZL ,AMP(4   ))            
      CALL FVIXXX(W(1,4   ),W(1,1   ),GAU ,ZERO    ,ZERO    ,W(1,13  ))    
      CALL JIOXXX(W(1,13  ),W(1,3   ),GWF ,WMASS   ,WWIDTH  ,W(1,14  ))    
      CALL FVIXXX(W(1,5   ),W(1,14  ),GWF ,ZERO    ,ZERO    ,W(1,15  ))    
      CALL IOVXXX(W(1,15  ),W(1,6   ),W(1,2   ),GWF ,AMP(5   ))            
      CALL JVVXXX(W(1,14  ),W(1,2   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     16  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,16  ),GAL ,AMP(6   ))            
      CALL JVVXXX(W(1,14  ),W(1,2   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     17  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,17  ),GZL ,AMP(7   ))            
      CALL FVIXXX(W(1,13  ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,18  ))    
      CALL JIOXXX(W(1,18  ),W(1,3   ),GAD ,ZERO    ,AWIDTH  ,W(1,19  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,19  ),GAL ,AMP(8   ))            
      CALL JIOXXX(W(1,18  ),W(1,3   ),GZD ,ZMASS   ,ZWIDTH  ,W(1,20  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,20  ),GZL ,AMP(9   ))            
      CALL FVOXXX(W(1,7   ),W(1,1   ),GAU ,ZERO    ,ZERO    ,W(1,21  ))    
      CALL JIOXXX(W(1,4   ),W(1,21  ),GAU ,ZERO    ,AWIDTH  ,W(1,22  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,22  ),GAL ,AMP(10  ))            
      CALL JIOXXX(W(1,4   ),W(1,21  ),GZU ,ZMASS   ,ZWIDTH  ,W(1,23  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,23  ),GZL ,AMP(11  ))            
      CALL JIOXXX(W(1,13  ),W(1,7   ),GAU ,ZERO    ,AWIDTH  ,W(1,24  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,24  ),GAL ,AMP(12  ))            
      CALL JIOXXX(W(1,13  ),W(1,7   ),GZU ,ZMASS   ,ZWIDTH  ,W(1,25  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,25  ),GZL ,AMP(13  ))            
      CALL FVIXXX(W(1,4   ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,26  ))    
      CALL JIOXXX(W(1,26  ),W(1,3   ),GAD ,ZERO    ,AWIDTH  ,W(1,27  ))    
      CALL FVIXXX(W(1,5   ),W(1,27  ),GAL ,ZERO    ,ZERO    ,W(1,28  ))    
      CALL IOVXXX(W(1,28  ),W(1,6   ),W(1,1   ),GAL ,AMP(14  ))            
      CALL JIOXXX(W(1,26  ),W(1,3   ),GZD ,ZMASS   ,ZWIDTH  ,W(1,29  ))    
      CALL FVIXXX(W(1,5   ),W(1,29  ),GZL ,ZERO    ,ZERO    ,W(1,30  ))    
      CALL IOVXXX(W(1,30  ),W(1,6   ),W(1,1   ),GAL ,AMP(15  ))            
      CALL IOVXXX(W(1,12  ),W(1,6   ),W(1,27  ),GAL ,AMP(16  ))            
      CALL IOVXXX(W(1,12  ),W(1,6   ),W(1,29  ),GZL ,AMP(17  ))            
      CALL FVOXXX(W(1,3   ),W(1,1   ),GAD ,ZERO    ,ZERO    ,W(1,31  ))    
      CALL JIOXXX(W(1,4   ),W(1,31  ),GWF ,WMASS   ,WWIDTH  ,W(1,32  ))    
      CALL FVIXXX(W(1,5   ),W(1,32  ),GWF ,ZERO    ,ZERO    ,W(1,33  ))    
      CALL IOVXXX(W(1,33  ),W(1,6   ),W(1,2   ),GWF ,AMP(18  ))            
      CALL JVVXXX(W(1,32  ),W(1,2   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     34  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,34  ),GAL ,AMP(19  ))            
      CALL JVVXXX(W(1,32  ),W(1,2   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     35  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,35  ),GZL ,AMP(20  ))            
      CALL FVIXXX(W(1,26  ),W(1,1   ),GAD ,ZERO    ,ZERO    ,W(1,36  ))    
      CALL JIOXXX(W(1,36  ),W(1,3   ),GAD ,ZERO    ,AWIDTH  ,W(1,37  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,37  ),GAL ,AMP(21  ))            
      CALL JIOXXX(W(1,36  ),W(1,3   ),GZD ,ZMASS   ,ZWIDTH  ,W(1,38  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,38  ),GZL ,AMP(22  ))            
      CALL FVOXXX(W(1,31  ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,39  ))    
      CALL JIOXXX(W(1,4   ),W(1,39  ),GAU ,ZERO    ,AWIDTH  ,W(1,40  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,40  ),GAL ,AMP(23  ))            
      CALL JIOXXX(W(1,4   ),W(1,39  ),GZU ,ZMASS   ,ZWIDTH  ,W(1,41  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,41  ),GZL ,AMP(24  ))            
      CALL JIOXXX(W(1,26  ),W(1,31  ),GAD ,ZERO    ,AWIDTH  ,W(1,42  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,42  ),GAL ,AMP(25  ))            
      CALL JIOXXX(W(1,26  ),W(1,31  ),GZD ,ZMASS   ,ZWIDTH  ,W(1,43  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,43  ),GZL ,AMP(26  ))            
      CALL JIOXXX(W(1,4   ),W(1,3   ),GWF ,WMASS   ,WWIDTH  ,W(1,44  ))    
      CALL JVVXXX(W(1,44  ),W(1,2   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     45  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,45  ),GAL ,ZERO    ,ZERO    ,W(1,46  ))    
      CALL IOVXXX(W(1,46  ),W(1,6   ),W(1,1   ),GAL ,AMP(27  ))            
      CALL JVVXXX(W(1,44  ),W(1,2   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     47  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,47  ),GZL ,ZERO    ,ZERO    ,W(1,48  ))    
      CALL IOVXXX(W(1,48  ),W(1,6   ),W(1,1   ),GAL ,AMP(28  ))            
      CALL FVIXXX(W(1,12  ),W(1,44  ),GWF ,ZERO    ,ZERO    ,W(1,49  ))    
      CALL IOVXXX(W(1,49  ),W(1,6   ),W(1,2   ),GWF ,AMP(29  ))            
      CALL IOVXXX(W(1,12  ),W(1,6   ),W(1,45  ),GAL ,AMP(30  ))            
      CALL IOVXXX(W(1,12  ),W(1,6   ),W(1,47  ),GZL ,AMP(31  ))            
      CALL JVVXXX(W(1,1   ),W(1,44  ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     50  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,50  ),GWF ,ZERO    ,ZERO    ,W(1,51  ))    
      CALL IOVXXX(W(1,51  ),W(1,6   ),W(1,2   ),GWF ,AMP(32  ))            
      CALL JVVXXX(W(1,50  ),W(1,2   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     52  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,52  ),GAL ,AMP(33  ))            
      CALL JVVXXX(W(1,50  ),W(1,2   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     53  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,53  ),GZL ,AMP(34  ))            
      CALL JVVXXX(W(1,2   ),W(1,1   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     54  ))                                                          
      CALL FVIXXX(W(1,4   ),W(1,54  ),GWF ,ZERO    ,ZERO    ,W(1,55  ))    
      CALL JIOXXX(W(1,55  ),W(1,3   ),GAD ,ZERO    ,AWIDTH  ,W(1,56  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,56  ),GAL ,AMP(35  ))            
      CALL JIOXXX(W(1,55  ),W(1,3   ),GZD ,ZMASS   ,ZWIDTH  ,W(1,57  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,57  ),GZL ,AMP(36  ))            
      CALL FVOXXX(W(1,3   ),W(1,54  ),GWF ,ZERO    ,ZERO    ,W(1,58  ))    
      CALL JIOXXX(W(1,4   ),W(1,58  ),GAU ,ZERO    ,AWIDTH  ,W(1,59  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,59  ),GAL ,AMP(37  ))            
      CALL JIOXXX(W(1,4   ),W(1,58  ),GZU ,ZMASS   ,ZWIDTH  ,W(1,60  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,60  ),GZL ,AMP(38  ))            
      CALL JW3WXX(W(1,44  ),W(1,1   ),W(1,2   ),GWWA ,GWWA ,ZERO    ,      
     &     AWIDTH  ,W(1,61  ))                                             
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,61  ),GAL ,AMP(39  ))            
      CALL JW3WXX(W(1,44  ),W(1,1   ),W(1,2   ),GWWA ,GWWZ ,ZMASS   ,      
     &     ZWIDTH  ,W(1,62  ))                                             
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,62  ),GZL ,AMP(40  ))            
      CALL FVIXXX(W(1,5   ),W(1,44  ),GWF ,ZERO    ,ZERO    ,W(1,63  ))    
      CALL FVIXXX(W(1,63  ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,64  ))    
      CALL IOVXXX(W(1,64  ),W(1,6   ),W(1,1   ),GAL ,AMP(41  ))            
      CALL JVVXXX(W(1,44  ),W(1,54  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     65  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,65  ),GAL ,AMP(42  ))            
      CALL JVVXXX(W(1,44  ),W(1,54  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     66  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,66  ),GZL ,AMP(43  ))            
      CALL IOVXXX(W(1,63  ),W(1,6   ),W(1,54  ),GWF ,AMP(44  ))            
      awm(mu,nu) = -AMP(   1)-AMP(   2)-AMP(   3)-AMP(   4)-AMP(   5)
     &             -AMP(   6)-AMP(   7)-AMP(   8)-AMP(   9)-AMP(  10)
     &             -AMP(  11)-AMP(  12)-AMP(  13)-AMP(  14)-AMP(  15)
     &             -AMP(  16)-AMP(  17)-AMP(  18)-AMP(  19)-AMP(  20)
     &             -AMP(  21)-AMP(  22)-AMP(  23)-AMP(  24)-AMP(  25)
     &             -AMP(  26)-AMP(  27)-AMP(  28)-AMP(  29)-AMP(  30)
     &             -AMP(  31)-AMP(  32)-AMP(  33)-AMP(  34)-AMP(  35)
     &             -AMP(  36)-AMP(  37)-AMP(  38)-AMP(  39)-AMP(  40)
     &             -AMP(  41)-AMP(  42)-AMP(  43)-AMP(  44)


      elseif (procid.eq.wmzhadjj .and. ztouu) then

      CALL FVOXXX(W(1,3   ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,7   ))    
      CALL JIOXXX(W(1,4   ),W(1,7   ),GZN ,ZMASS   ,ZWIDTH  ,W(1,8   ))    
      CALL FVIXXX(W(1,5   ),W(1,8   ),GZU ,ZERO    ,ZERO    ,W(1,9   ))    
      CALL IOVXXX(W(1,9   ),W(1,6   ),W(1,1   ),GAU ,AMP(1   ))            
      CALL FVIXXX(W(1,5   ),W(1,1   ),GAU ,ZERO    ,ZERO    ,W(1,10  ))    
      CALL IOVXXX(W(1,10  ),W(1,6   ),W(1,8   ),GZU ,AMP(2   ))            
      CALL FVIXXX(W(1,4   ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,11  ))    
      CALL JIOXXX(W(1,11  ),W(1,3   ),GAL ,ZERO    ,AWIDTH  ,W(1,12  ))    
      CALL FVIXXX(W(1,5   ),W(1,12  ),GAU ,ZERO    ,ZERO    ,W(1,13  ))    
      CALL IOVXXX(W(1,13  ),W(1,6   ),W(1,1   ),GAU ,AMP(3   ))            
      CALL JIOXXX(W(1,11  ),W(1,3   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,14  ))    
      CALL FVIXXX(W(1,5   ),W(1,14  ),GZU ,ZERO    ,ZERO    ,W(1,15  ))    
      CALL IOVXXX(W(1,15  ),W(1,6   ),W(1,1   ),GAU ,AMP(4   ))            
      CALL IOVXXX(W(1,10  ),W(1,6   ),W(1,12  ),GAU ,AMP(5   ))            
      CALL IOVXXX(W(1,10  ),W(1,6   ),W(1,14  ),GZU ,AMP(6   ))            
      CALL FVOXXX(W(1,3   ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,16  ))    
      CALL JIOXXX(W(1,4   ),W(1,16  ),GWF ,WMASS   ,WWIDTH  ,W(1,17  ))    
      CALL JVVXXX(W(1,17  ),W(1,2   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     18  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,18  ),GAU ,AMP(7   ))            
      CALL JVVXXX(W(1,17  ),W(1,2   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     19  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,19  ),GZU ,AMP(8   ))            
      CALL FVIXXX(W(1,5   ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,20  ))    
      CALL IOVXXX(W(1,20  ),W(1,6   ),W(1,17  ),GWF ,AMP(9   ))            
      CALL FVIXXX(W(1,11  ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,21  ))    
      CALL JIOXXX(W(1,21  ),W(1,3   ),GAL ,ZERO    ,AWIDTH  ,W(1,22  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,22  ),GAU ,AMP(10  ))            
      CALL JIOXXX(W(1,21  ),W(1,3   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,23  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,23  ),GZU ,AMP(11  ))            
      CALL FVOXXX(W(1,16  ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,24  ))    
      CALL JIOXXX(W(1,4   ),W(1,24  ),GZN ,ZMASS   ,ZWIDTH  ,W(1,25  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,25  ),GZU ,AMP(12  ))            
      CALL JIOXXX(W(1,11  ),W(1,16  ),GAL ,ZERO    ,AWIDTH  ,W(1,26  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,26  ),GAU ,AMP(13  ))            
      CALL JIOXXX(W(1,11  ),W(1,16  ),GZL ,ZMASS   ,ZWIDTH  ,W(1,27  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,27  ),GZU ,AMP(14  ))            
      CALL JIOXXX(W(1,4   ),W(1,3   ),GWF ,WMASS   ,WWIDTH  ,W(1,28  ))    
      CALL JVVXXX(W(1,28  ),W(1,2   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     29  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,29  ),GAU ,ZERO    ,ZERO    ,W(1,30  ))    
      CALL IOVXXX(W(1,30  ),W(1,6   ),W(1,1   ),GAU ,AMP(15  ))            
      CALL JVVXXX(W(1,28  ),W(1,2   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     31  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,31  ),GZU ,ZERO    ,ZERO    ,W(1,32  ))    
      CALL IOVXXX(W(1,32  ),W(1,6   ),W(1,1   ),GAU ,AMP(16  ))            
      CALL FVIXXX(W(1,10  ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,33  ))    
      CALL IOVXXX(W(1,33  ),W(1,6   ),W(1,28  ),GWF ,AMP(17  ))            
      CALL IOVXXX(W(1,10  ),W(1,6   ),W(1,29  ),GAU ,AMP(18  ))            
      CALL IOVXXX(W(1,10  ),W(1,6   ),W(1,31  ),GZU ,AMP(19  ))            
      CALL FVIXXX(W(1,20  ),W(1,28  ),GWF ,ZERO    ,ZERO    ,W(1,34  ))    
      CALL IOVXXX(W(1,34  ),W(1,6   ),W(1,1   ),GAU ,AMP(20  ))            
      CALL JVVXXX(W(1,1   ),W(1,28  ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     35  ))                                                          
      CALL JVVXXX(W(1,35  ),W(1,2   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     36  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,36  ),GAU ,AMP(21  ))            
      CALL JVVXXX(W(1,35  ),W(1,2   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     37  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,37  ),GZU ,AMP(22  ))            
      CALL FVIXXX(W(1,20  ),W(1,1   ),GAD ,ZERO    ,ZERO    ,W(1,38  ))    
      CALL IOVXXX(W(1,38  ),W(1,6   ),W(1,28  ),GWF ,AMP(23  ))            
      CALL IOVXXX(W(1,20  ),W(1,6   ),W(1,35  ),GWF ,AMP(24  ))            
      CALL JVVXXX(W(1,2   ),W(1,1   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     39  ))                                                          
      CALL FVIXXX(W(1,4   ),W(1,39  ),GWF ,ZERO    ,ZERO    ,W(1,40  ))    
      CALL JIOXXX(W(1,40  ),W(1,3   ),GAL ,ZERO    ,AWIDTH  ,W(1,41  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,41  ),GAU ,AMP(25  ))            
      CALL JIOXXX(W(1,40  ),W(1,3   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,42  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,42  ),GZU ,AMP(26  ))            
      CALL FVOXXX(W(1,3   ),W(1,39  ),GWF ,ZERO    ,ZERO    ,W(1,43  ))    
      CALL JIOXXX(W(1,4   ),W(1,43  ),GZN ,ZMASS   ,ZWIDTH  ,W(1,44  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,44  ),GZU ,AMP(27  ))            
      CALL JW3WXX(W(1,28  ),W(1,1   ),W(1,2   ),GWWA ,GWWA ,ZERO    ,      
     &     AWIDTH  ,W(1,45  ))                                             
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,45  ),GAU ,AMP(28  ))            
      CALL JW3WXX(W(1,28  ),W(1,1   ),W(1,2   ),GWWA ,GWWZ ,ZMASS   ,      
     &     ZWIDTH  ,W(1,46  ))                                             
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,46  ),GZU ,AMP(29  ))            
      CALL JVVXXX(W(1,28  ),W(1,39  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     47  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,47  ),GAU ,AMP(30  ))            
      CALL JVVXXX(W(1,28  ),W(1,39  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     48  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,48  ),GZU ,AMP(31  ))            
      CALL FVIXXX(W(1,5   ),W(1,39  ),GWF ,ZERO    ,ZERO    ,W(1,49  ))    
      CALL IOVXXX(W(1,49  ),W(1,6   ),W(1,28  ),GWF ,AMP(32  ))            
      awm(mu,nu) = -AMP(   1)-AMP(   2)-AMP(   3)-AMP(   4)-AMP(   5)
     &             -AMP(   6)-AMP(   7)-AMP(   8)-AMP(   9)-AMP(  10)
     &             -AMP(  11)-AMP(  12)-AMP(  13)-AMP(  14)-AMP(  15)
     &             -AMP(  16)-AMP(  17)-AMP(  18)-AMP(  19)-AMP(  20)
     &             -AMP(  21)-AMP(  22)-AMP(  23)-AMP(  24)-AMP(  25)
     &             -AMP(  26)-AMP(  27)-AMP(  28)-AMP(  29)-AMP(  30)
     &             -AMP(  31)-AMP(  32)


      elseif (procid.eq.wmzhadjj .and. ztodd) then

      CALL FVOXXX(W(1,3   ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,7   ))    
      CALL JIOXXX(W(1,4   ),W(1,7   ),GZN ,ZMASS   ,ZWIDTH  ,W(1,8   ))    
      CALL FVIXXX(W(1,5   ),W(1,8   ),GZD ,ZERO    ,ZERO    ,W(1,9   ))    
      CALL IOVXXX(W(1,9   ),W(1,6   ),W(1,1   ),GAD ,AMP(1   ))            
      CALL FVIXXX(W(1,5   ),W(1,1   ),GAD ,ZERO    ,ZERO    ,W(1,10  ))    
      CALL IOVXXX(W(1,10  ),W(1,6   ),W(1,8   ),GZD ,AMP(2   ))            
      CALL FVIXXX(W(1,4   ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,11  ))    
      CALL JIOXXX(W(1,11  ),W(1,3   ),GAL ,ZERO    ,AWIDTH  ,W(1,12  ))    
      CALL FVIXXX(W(1,5   ),W(1,12  ),GAD ,ZERO    ,ZERO    ,W(1,13  ))    
      CALL IOVXXX(W(1,13  ),W(1,6   ),W(1,1   ),GAD ,AMP(3   ))            
      CALL JIOXXX(W(1,11  ),W(1,3   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,14  ))    
      CALL FVIXXX(W(1,5   ),W(1,14  ),GZD ,ZERO    ,ZERO    ,W(1,15  ))    
      CALL IOVXXX(W(1,15  ),W(1,6   ),W(1,1   ),GAD ,AMP(4   ))            
      CALL IOVXXX(W(1,10  ),W(1,6   ),W(1,12  ),GAD ,AMP(5   ))            
      CALL IOVXXX(W(1,10  ),W(1,6   ),W(1,14  ),GZD ,AMP(6   ))            
      CALL FVOXXX(W(1,3   ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,16  ))    
      CALL JIOXXX(W(1,4   ),W(1,16  ),GWF ,WMASS   ,WWIDTH  ,W(1,17  ))    
      CALL FVIXXX(W(1,5   ),W(1,17  ),GWF ,ZERO    ,ZERO    ,W(1,18  ))    
      CALL IOVXXX(W(1,18  ),W(1,6   ),W(1,2   ),GWF ,AMP(7   ))            
      CALL JVVXXX(W(1,17  ),W(1,2   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     19  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,19  ),GAD ,AMP(8   ))            
      CALL JVVXXX(W(1,17  ),W(1,2   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     20  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,20  ),GZD ,AMP(9   ))            
      CALL FVIXXX(W(1,11  ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,21  ))    
      CALL JIOXXX(W(1,21  ),W(1,3   ),GAL ,ZERO    ,AWIDTH  ,W(1,22  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,22  ),GAD ,AMP(10  ))            
      CALL JIOXXX(W(1,21  ),W(1,3   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,23  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,23  ),GZD ,AMP(11  ))            
      CALL FVOXXX(W(1,16  ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,24  ))    
      CALL JIOXXX(W(1,4   ),W(1,24  ),GZN ,ZMASS   ,ZWIDTH  ,W(1,25  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,25  ),GZD ,AMP(12  ))            
      CALL JIOXXX(W(1,11  ),W(1,16  ),GAL ,ZERO    ,AWIDTH  ,W(1,26  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,26  ),GAD ,AMP(13  ))            
      CALL JIOXXX(W(1,11  ),W(1,16  ),GZL ,ZMASS   ,ZWIDTH  ,W(1,27  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,27  ),GZD ,AMP(14  ))            
      CALL JIOXXX(W(1,4   ),W(1,3   ),GWF ,WMASS   ,WWIDTH  ,W(1,28  ))    
      CALL JVVXXX(W(1,28  ),W(1,2   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     29  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,29  ),GAD ,ZERO    ,ZERO    ,W(1,30  ))    
      CALL IOVXXX(W(1,30  ),W(1,6   ),W(1,1   ),GAD ,AMP(15  ))            
      CALL JVVXXX(W(1,28  ),W(1,2   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     31  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,31  ),GZD ,ZERO    ,ZERO    ,W(1,32  ))    
      CALL IOVXXX(W(1,32  ),W(1,6   ),W(1,1   ),GAD ,AMP(16  ))            
      CALL FVIXXX(W(1,10  ),W(1,28  ),GWF ,ZERO    ,ZERO    ,W(1,33  ))    
      CALL IOVXXX(W(1,33  ),W(1,6   ),W(1,2   ),GWF ,AMP(17  ))            
      CALL IOVXXX(W(1,10  ),W(1,6   ),W(1,29  ),GAD ,AMP(18  ))            
      CALL IOVXXX(W(1,10  ),W(1,6   ),W(1,31  ),GZD ,AMP(19  ))            
      CALL JVVXXX(W(1,1   ),W(1,28  ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     34  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,34  ),GWF ,ZERO    ,ZERO    ,W(1,35  ))    
      CALL IOVXXX(W(1,35  ),W(1,6   ),W(1,2   ),GWF ,AMP(20  ))            
      CALL JVVXXX(W(1,34  ),W(1,2   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     36  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,36  ),GAD ,AMP(21  ))            
      CALL JVVXXX(W(1,34  ),W(1,2   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     37  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,37  ),GZD ,AMP(22  ))            
      CALL JVVXXX(W(1,2   ),W(1,1   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     38  ))                                                          
      CALL FVIXXX(W(1,4   ),W(1,38  ),GWF ,ZERO    ,ZERO    ,W(1,39  ))    
      CALL JIOXXX(W(1,39  ),W(1,3   ),GAL ,ZERO    ,AWIDTH  ,W(1,40  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,40  ),GAD ,AMP(23  ))            
      CALL JIOXXX(W(1,39  ),W(1,3   ),GZL ,ZMASS   ,ZWIDTH  ,W(1,41  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,41  ),GZD ,AMP(24  ))            
      CALL FVOXXX(W(1,3   ),W(1,38  ),GWF ,ZERO    ,ZERO    ,W(1,42  ))    
      CALL JIOXXX(W(1,4   ),W(1,42  ),GZN ,ZMASS   ,ZWIDTH  ,W(1,43  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,43  ),GZD ,AMP(25  ))            
      CALL JW3WXX(W(1,28  ),W(1,1   ),W(1,2   ),GWWA ,GWWA ,ZERO    ,      
     &     AWIDTH  ,W(1,44  ))                                             
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,44  ),GAD ,AMP(26  ))            
      CALL JW3WXX(W(1,28  ),W(1,1   ),W(1,2   ),GWWA ,GWWZ ,ZMASS   ,      
     &     ZWIDTH  ,W(1,45  ))                                             
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,45  ),GZD ,AMP(27  ))            
      CALL FVIXXX(W(1,5   ),W(1,28  ),GWF ,ZERO    ,ZERO    ,W(1,46  ))    
      CALL FVIXXX(W(1,46  ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,47  ))    
      CALL IOVXXX(W(1,47  ),W(1,6   ),W(1,1   ),GAD ,AMP(28  ))            
      CALL FVIXXX(W(1,46  ),W(1,1   ),GAU ,ZERO    ,ZERO    ,W(1,48  ))    
      CALL IOVXXX(W(1,48  ),W(1,6   ),W(1,2   ),GWF ,AMP(29  ))            
      CALL JVVXXX(W(1,28  ),W(1,38  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     49  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,49  ),GAD ,AMP(30  ))            
      CALL JVVXXX(W(1,28  ),W(1,38  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     50  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,50  ),GZD ,AMP(31  ))            
      CALL IOVXXX(W(1,46  ),W(1,6   ),W(1,38  ),GWF ,AMP(32  ))            
      awm(mu,nu) = -AMP(   1)-AMP(   2)-AMP(   3)-AMP(   4)-AMP(   5)
     &             -AMP(   6)-AMP(   7)-AMP(   8)-AMP(   9)-AMP(  10)
     &             -AMP(  11)-AMP(  12)-AMP(  13)-AMP(  14)-AMP(  15)
     &             -AMP(  16)-AMP(  17)-AMP(  18)-AMP(  19)-AMP(  20)
     &             -AMP(  21)-AMP(  22)-AMP(  23)-AMP(  24)-AMP(  25)
     &             -AMP(  26)-AMP(  27)-AMP(  28)-AMP(  29)-AMP(  30)
     &             -AMP(  31)-AMP(  32)


      endif
      
          ENDDO
      ENDDO
      END
c ------------------------------------------------------------	

