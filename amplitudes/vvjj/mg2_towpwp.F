      subroutine qqwppqq_mg2(p,pbar,v,sign, nlo, 
     1                     bos,ucds)
      
      implicit none                

#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "tensorz.inc"

#include "VBFNLO/utilities/process.inc"

      double precision pbar(0:3,4+max_v)
      double precision ucds
      double precision ucds2
      double precision p(0:3,max_p),v(0:3,max_v)
      integer  sign(4+max_v), nlo, i, j,mu,nu,bos
      real*8 q12(0:4,3), q34(0:4,3)
      complex*16 mg2_tens_1d(0:5)
      complex*16 mg2_tens_2d(0:3,0:3)
      complex*16 mg2_tens_2d_comb(0:3,0:3,2)
      double precision sumtest

      call Calc_Momentum_Transfer(p, v, q12,q34,1)

      j = 1
         
#ifdef WITH_MADGRAPH

         SELECT CASE(procid)
         CASE(WPhadWPjj, WMhadWMjj)            ! semileptonic decay modes

            call anomal_formfactor(q12(0,j),q34(0,j),qe(0),qu(0))

            call wpptowpp_mg2_had(q12(0,j),q34(0,j),v,wppwpp(0,0,j))  

            call WPVtoWP_mg2_had(1,1,q34(0,j),v,awen(0,0,1,j),zwen(0,0,1,j)) !Wu on upper
            call WPVtoWP_mg2_had(1,2,q12(0,j),v,awen(0,0,2,j),zwen(0,0,2,j)) !Wu on lower
            call WPVtoWP_mg2_had(1,1,q12(0,j),v,waen(0,0,2,j),wzen(0,0,2,j)) !Wu on lower
            call WPVtoWP_mg2_had(1,2,q34(0,j),v,waen(0,0,1,j),wzen(0,0,1,j)) !Wu on upper

            call WPVtoWP_mg2_had(2,1,q34(0,j),v,awmu(0,0,1,j),zwmu(0,0,1,j)) !We on upper
            call WPVtoWP_mg2_had(2,2,q12(0,j),v,awmu(0,0,2,j),zwmu(0,0,2,j)) !We on lower
            call WPVtoWP_mg2_had(2,1,q12(0,j),v,wamu(0,0,2,j),wzmu(0,0,2,j)) !We on lower
            call WPVtoWP_mg2_had(2,2,q34(0,j),v,wamu(0,0,1,j),wzmu(0,0,1,j)) !We on upper

         call qqwppqq(pbar,sign, nlo,1,bos,ucds2)

            call anomal_formfactor(q12(0,j),q34(0,j),qe(0),qu(0))

            call wpptowpp_had(q12(0,j),q34(0,j),v,wppwpp(0,0,j))  

            call WPVtoWP_had(1,1,q34(0,j),v,awen(0,0,1,j),zwen(0,0,1,j)) !Wu on upper
            call WPVtoWP_had(1,2,q12(0,j),v,awen(0,0,2,j),zwen(0,0,2,j)) !Wu on lower
            call WPVtoWP_had(1,1,q12(0,j),v,waen(0,0,2,j),wzen(0,0,2,j)) !Wu on lower
            call WPVtoWP_had(1,2,q34(0,j),v,waen(0,0,1,j),wzen(0,0,1,j)) !Wu on upper

            call WPVtoWP_had(2,1,q34(0,j),v,awmu(0,0,1,j),zwmu(0,0,1,j)) !We on upper
            call WPVtoWP_had(2,2,q12(0,j),v,awmu(0,0,2,j),zwmu(0,0,2,j)) !We on lower
            call WPVtoWP_had(2,1,q12(0,j),v,wamu(0,0,2,j),wzmu(0,0,2,j)) !We on lower
            call WPVtoWP_had(2,2,q34(0,j),v,wamu(0,0,1,j),wzmu(0,0,1,j)) !We on upper


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
         call wpptowpp_mg2_had(q12(0,j),q34(0,j),v,mg2_tens_2d(0,0))  
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( wppwpp(mu,nu,1) / mg2_tens_2d(mu,nu) - 1)
           ! print*, mu,nu, wppwpp(mu,nu,1), mg2_tens_2d(mu,nu), abs( wppwpp(mu,nu,1) / mg2_tens_2d(mu,nu) - 1)
           mg2_tens_2d(mu,nu) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor w+w+ to w+w+ :", sumtest

!         print*, ""
         call WPVtoWP_mg2_had(1,1,q34(0,j),v,mg2_tens_2d_comb(0,0,1),mg2_tens_2d_comb(0,0,2)) !Wu on upper
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( awen(mu,nu,1,j) / mg2_tens_2d_comb(mu,nu,1) - 1)
           ! print*, mu,nu, awen(mu,nu,1,j), mg2_tens_2d_comb(mu,nu,1), abs( awen(mu,nu,1,1) / mg2_tens_2d_comb(mu,nu,1) - 1)
           mg2_tens_2d_comb(mu,nu,1) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor aw to w+ 1 had:", sumtest
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( zwen(mu,nu,1,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           ! print*, mu,nu, zwen(mu,nu,1,j), mg2_tens_2d_comb(mu,nu,2), abs( zwen(mu,nu,1,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           mg2_tens_2d_comb(mu,nu,2) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor zw to w+ 1 had:", sumtest

!         print*, ""
         call WPVtoWP_mg2_had(1,2,q12(0,j),v,mg2_tens_2d_comb(0,0,1),mg2_tens_2d_comb(0,0,2)) !Wu on lower
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( awen(mu,nu,2,j) / mg2_tens_2d_comb(mu,nu,1) - 1)
           ! print*, mu,nu, awen(mu,nu,2,j), mg2_tens_2d_comb(mu,nu,1), abs( awen(mu,nu,2,1) / mg2_tens_2d_comb(mu,nu,1) - 1)
           mg2_tens_2d_comb(mu,nu,1) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor aw to w+ 2 had:", sumtest
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( zwen(mu,nu,2,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           ! print*, mu,nu, zwen(mu,nu,2,j), mg2_tens_2d_comb(mu,nu,2), abs( zwen(mu,nu,2,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           mg2_tens_2d_comb(mu,nu,2) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor zw to w+ 2 had:", sumtest

!         print*, ""
         call WPVtoWP_mg2_had(1,2,q34(0,j),v,mg2_tens_2d_comb(0,0,1),mg2_tens_2d_comb(0,0,2)) !Wu on upper
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( waen(mu,nu,1,j) / mg2_tens_2d_comb(mu,nu,1) - 1)
           ! print*, mu,nu, waen(mu,nu,1,j), mg2_tens_2d_comb(mu,nu,1), abs( waen(mu,nu,1,1) / mg2_tens_2d_comb(mu,nu,1) - 1)
           mg2_tens_2d_comb(mu,nu,1) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor wa to w+ 1 had:", sumtest
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( wzen(mu,nu,1,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           ! print*, mu,nu, wzen(mu,nu,1,j), mg2_tens_2d_comb(mu,nu,2), abs( wzen(mu,nu,1,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           mg2_tens_2d_comb(mu,nu,2) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor wz to w+ 1 had:", sumtest

!         print*, ""
         call WPVtoWP_mg2_had(1,1,q12(0,j),v,mg2_tens_2d_comb(0,0,1),mg2_tens_2d_comb(0,0,2)) !Wu on lower
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( waen(mu,nu,2,j) / mg2_tens_2d_comb(mu,nu,1) - 1)
           ! print*, mu,nu, waen(mu,nu,2,j), mg2_tens_2d_comb(mu,nu,1), abs( waen(mu,nu,2,1) / mg2_tens_2d_comb(mu,nu,1) - 1)
           mg2_tens_2d_comb(mu,nu,1) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor wa to w+ 2 had:", sumtest
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( wzen(mu,nu,2,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           ! print*, mu,nu, wzen(mu,nu,2,j), mg2_tens_2d_comb(mu,nu,2), abs( wzen(mu,nu,2,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           mg2_tens_2d_comb(mu,nu,2) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor wz to w+ 2 had:", sumtest

!         print*, ""
         call WPVtoWP_mg2_had(2,1,q34(0,j),v,mg2_tens_2d_comb(0,0,1),mg2_tens_2d_comb(0,0,2)) !Wu on upper
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( awmu(mu,nu,1,j) / mg2_tens_2d_comb(mu,nu,1) - 1)
           ! print*, mu,nu, awmu(mu,nu,1,j), mg2_tens_2d_comb(mu,nu,1), abs( awmu(mu,nu,1,1) / mg2_tens_2d_comb(mu,nu,1) - 1)
           mg2_tens_2d_comb(mu,nu,1) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor aw to w+ 1 lep:", sumtest
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( zwmu(mu,nu,1,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           ! print*, mu,nu, zwmu(mu,nu,1,j), mg2_tens_2d_comb(mu,nu,2), abs( zwmu(mu,nu,1,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           mg2_tens_2d_comb(mu,nu,2) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor zw to w+ 1 lep:", sumtest

!         print*, ""
         call WPVtoWP_mg2_had(2,2,q12(0,j),v,mg2_tens_2d_comb(0,0,1),mg2_tens_2d_comb(0,0,2)) !Wu on lower
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( awmu(mu,nu,2,j) / mg2_tens_2d_comb(mu,nu,1) - 1)
           ! print*, mu,nu, awmu(mu,nu,2,j), mg2_tens_2d_comb(mu,nu,1), abs( awmu(mu,nu,2,1) / mg2_tens_2d_comb(mu,nu,1) - 1)
           mg2_tens_2d_comb(mu,nu,1) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor aw to w+ 2 lep:", sumtest
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( zwmu(mu,nu,2,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           ! print*, mu,nu, zwmu(mu,nu,2,j), mg2_tens_2d_comb(mu,nu,2), abs( zwmu(mu,nu,2,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           mg2_tens_2d_comb(mu,nu,2) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor zw to w+ 2 lep:", sumtest

!         print*, ""
         call WPVtoWP_mg2_had(2,2,q34(0,j),v,mg2_tens_2d_comb(0,0,1),mg2_tens_2d_comb(0,0,2)) !Wu on upper
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( wamu(mu,nu,1,j) / mg2_tens_2d_comb(mu,nu,1) - 1)
           ! print*, mu,nu, wamu(mu,nu,1,j), mg2_tens_2d_comb(mu,nu,1), abs( wamu(mu,nu,1,1) / mg2_tens_2d_comb(mu,nu,1) - 1)
           mg2_tens_2d_comb(mu,nu,1) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor wa to w+ 1 lep:", sumtest
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( wzmu(mu,nu,1,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           ! print*, mu,nu, wzmu(mu,nu,1,j), mg2_tens_2d_comb(mu,nu,2), abs( wzmu(mu,nu,1,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           mg2_tens_2d_comb(mu,nu,2) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor wz to w+ 1 lep:", sumtest

!         print*, ""
         call WPVtoWP_mg2_had(2,1,q12(0,j),v,mg2_tens_2d_comb(0,0,1),mg2_tens_2d_comb(0,0,2)) !Wu on lower
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( wamu(mu,nu,2,j) / mg2_tens_2d_comb(mu,nu,1) - 1)
           ! print*, mu,nu, wamu(mu,nu,2,j), mg2_tens_2d_comb(mu,nu,1), abs( wamu(mu,nu,2,1) / mg2_tens_2d_comb(mu,nu,1) - 1)
           mg2_tens_2d_comb(mu,nu,1) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor wa to w+ 2 lep:", sumtest
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( wzmu(mu,nu,2,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           ! print*, mu,nu, wzmu(mu,nu,2,j), mg2_tens_2d_comb(mu,nu,2), abs( wzmu(mu,nu,2,j) / mg2_tens_2d_comb(mu,nu,2) - 1)
           mg2_tens_2d_comb(mu,nu,2) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor wz to w+ 2 lep:", sumtest


         print*, ""
         print*,"1 - ucds/ucds2 =   ",abs(1-ucds/ucds2)
         read(*,*)

         end select

#endif


      end


      subroutine WPPTOWPP_mg2_had(q1,q2,p,ww)
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : w+ w+ -> e+ ve mu+ vm  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=  37) 
      INTEGER    NWAVEFUNCS 
      PARAMETER (NWAVEFUNCS=  58) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 Q1(0:3), Q2(0:3), P(0:3,3:6)
      complex*16 ww(0:3,0:3)
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
c
C ----------
C BEGIN CODE
C ----------
c read in lepton an W states from common block:
         do i=1,6
            w(i,3) = wel(i)
            w(i,4) = wve(i)
            w(i,5) = wmu(i)
            w(i,6) = wvu(i)
         enddo

c	print*,'q1,q2 in ww = ',q1(0),q2(0)

      DO MU = 0,3
         DO NU = 0,3
      
      CALL VCARTX(Q1(0),wMASS ,wwidth, MU,-1,W(1,1   ))        
      CALL VCARTX(Q2(0),wMASS ,wwidth, NU,-1,W(1,2   )) 
               
      CALL FVIXXX(W(1,3   ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,7   ))    
      CALL JIOXXX(W(1,7   ),W(1,4   ),GAU ,ZERO    ,AWIDTH  ,W(1,8   ))    
      CALL FVIXXX(W(1,5   ),W(1,8   ),GAL ,ZERO    ,ZERO    ,W(1,9   ))    
      CALL IOVXXX(W(1,9   ),W(1,6   ),W(1,1   ),GWF ,AMP(1   ))            
      CALL JIOXXX(W(1,7   ),W(1,4   ),GZU ,ZMASS   ,ZWIDTH  ,W(1,10  ))    
      CALL FVIXXX(W(1,5   ),W(1,10  ),GZL ,ZERO    ,ZERO    ,W(1,11  ))    
      CALL IOVXXX(W(1,11  ),W(1,6   ),W(1,1   ),GWF ,AMP(2   ))            
      CALL JVVXXX(W(1,8   ),W(1,1   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     12  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,12  ),GWF ,AMP(3   ))            
      CALL JVVXXX(W(1,10  ),W(1,1   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     13  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,13  ),GWF ,AMP(4   ))            
      CALL FVIXXX(W(1,5   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,14  ))    
      CALL IOVXXX(W(1,14  ),W(1,6   ),W(1,10  ),GZN ,AMP(5   ))            
      CALL FVOXXX(W(1,4   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,15  ))    
      CALL JIOXXX(W(1,3   ),W(1,15  ),GAD ,ZERO    ,AWIDTH  ,W(1,16  ))    
      CALL FVIXXX(W(1,5   ),W(1,16  ),GAL ,ZERO    ,ZERO    ,W(1,17  ))    
      CALL IOVXXX(W(1,17  ),W(1,6   ),W(1,2   ),GWF ,AMP(6   ))            
      CALL JIOXXX(W(1,3   ),W(1,15  ),GZD ,ZMASS   ,ZWIDTH  ,W(1,18  ))    
      CALL FVIXXX(W(1,5   ),W(1,18  ),GZL ,ZERO    ,ZERO    ,W(1,19  ))    
      CALL IOVXXX(W(1,19  ),W(1,6   ),W(1,2   ),GWF ,AMP(7   ))            
      CALL JVVXXX(W(1,16  ),W(1,2   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     20  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,20  ),GWF ,AMP(8   ))            
      CALL JVVXXX(W(1,18  ),W(1,2   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     21  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,21  ),GWF ,AMP(9   ))            
      CALL FVIXXX(W(1,5   ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,22  ))    
      CALL IOVXXX(W(1,22  ),W(1,6   ),W(1,18  ),GZN ,AMP(10  ))            
      CALL JIOXXX(W(1,7   ),W(1,15  ),GWF ,WMASS   ,WWIDTH  ,W(1,23  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,23  ),GWF ,AMP(11  ))            
      CALL FVOXXX(W(1,4   ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,24  ))    
      CALL JIOXXX(W(1,3   ),W(1,24  ),GAD ,ZERO    ,AWIDTH  ,W(1,25  ))    
      CALL FVIXXX(W(1,5   ),W(1,25  ),GAL ,ZERO    ,ZERO    ,W(1,26  ))    
      CALL IOVXXX(W(1,26  ),W(1,6   ),W(1,1   ),GWF ,AMP(12  ))            
      CALL JIOXXX(W(1,3   ),W(1,24  ),GZD ,ZMASS   ,ZWIDTH  ,W(1,27  ))    
      CALL FVIXXX(W(1,5   ),W(1,27  ),GZL ,ZERO    ,ZERO    ,W(1,28  ))    
      CALL IOVXXX(W(1,28  ),W(1,6   ),W(1,1   ),GWF ,AMP(13  ))            
      CALL JVVXXX(W(1,25  ),W(1,1   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     29  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,29  ),GWF ,AMP(14  ))            
      CALL JVVXXX(W(1,27  ),W(1,1   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     30  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,30  ),GWF ,AMP(15  ))            
      CALL IOVXXX(W(1,14  ),W(1,6   ),W(1,27  ),GZN ,AMP(16  ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,31  ))    
      CALL JIOXXX(W(1,31  ),W(1,4   ),GAU ,ZERO    ,AWIDTH  ,W(1,32  ))    
      CALL FVIXXX(W(1,5   ),W(1,32  ),GAL ,ZERO    ,ZERO    ,W(1,33  ))    
      CALL IOVXXX(W(1,33  ),W(1,6   ),W(1,2   ),GWF ,AMP(17  ))            
      CALL JIOXXX(W(1,31  ),W(1,4   ),GZU ,ZMASS   ,ZWIDTH  ,W(1,34  ))    
      CALL FVIXXX(W(1,5   ),W(1,34  ),GZL ,ZERO    ,ZERO    ,W(1,35  ))    
      CALL IOVXXX(W(1,35  ),W(1,6   ),W(1,2   ),GWF ,AMP(18  ))            
      CALL JVVXXX(W(1,32  ),W(1,2   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     36  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,36  ),GWF ,AMP(19  ))            
      CALL JVVXXX(W(1,34  ),W(1,2   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     37  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,37  ),GWF ,AMP(20  ))            
      CALL IOVXXX(W(1,22  ),W(1,6   ),W(1,34  ),GZN ,AMP(21  ))            
      CALL JIOXXX(W(1,31  ),W(1,24  ),GWF ,WMASS   ,WWIDTH  ,W(1,38  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,38  ),GWF ,AMP(22  ))            
      CALL JIOXXX(W(1,3   ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,39  ))    
      CALL JVVXXX(W(1,2   ),W(1,39  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     40  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,40  ),GAL ,ZERO    ,ZERO    ,W(1,41  ))    
      CALL IOVXXX(W(1,41  ),W(1,6   ),W(1,1   ),GWF ,AMP(23  ))            
      CALL JVVXXX(W(1,2   ),W(1,39  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     42  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,42  ),GZL ,ZERO    ,ZERO    ,W(1,43  ))    
      CALL IOVXXX(W(1,43  ),W(1,6   ),W(1,1   ),GWF ,AMP(24  ))            
      CALL FVIXXX(W(1,14  ),W(1,39  ),GWF ,ZERO    ,ZERO    ,W(1,44  ))    
      CALL IOVXXX(W(1,44  ),W(1,6   ),W(1,2   ),GWF ,AMP(25  ))            
      CALL JVVXXX(W(1,40  ),W(1,1   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     45  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,45  ),GWF ,AMP(26  ))            
      CALL JVVXXX(W(1,42  ),W(1,1   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     46  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,46  ),GWF ,AMP(27  ))            
      CALL HVVXXX(W(1,2   ),W(1,39  ),GWWH ,HMASS   ,HWIDTH  ,W(1,         
     &     47  ))                                                          
      CALL JVSXXX(W(1,1   ),W(1,47  ),GWWH ,WMASS   ,WWIDTH  ,W(1,         
     &     48  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,48  ),GWF ,AMP(28  ))            
      CALL IOVXXX(W(1,14  ),W(1,6   ),W(1,42  ),GZN ,AMP(29  ))            
      CALL FVIXXX(W(1,22  ),W(1,39  ),GWF ,ZERO    ,ZERO    ,W(1,49  ))    
      CALL IOVXXX(W(1,49  ),W(1,6   ),W(1,1   ),GWF ,AMP(30  ))            
      CALL JVVXXX(W(1,1   ),W(1,39  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     50  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,50  ),GAL ,ZERO    ,ZERO    ,W(1,51  ))    
      CALL IOVXXX(W(1,51  ),W(1,6   ),W(1,2   ),GWF ,AMP(31  ))            
      CALL JVVXXX(W(1,1   ),W(1,39  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     52  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,52  ),GZL ,ZERO    ,ZERO    ,W(1,53  ))    
      CALL IOVXXX(W(1,53  ),W(1,6   ),W(1,2   ),GWF ,AMP(32  ))            
      CALL JVVXXX(W(1,50  ),W(1,2   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     54  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,54  ),GWF ,AMP(33  ))            
      CALL JVVXXX(W(1,52  ),W(1,2   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     55  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,55  ),GWF ,AMP(34  ))            
      CALL HVVXXX(W(1,1   ),W(1,39  ),GWWH ,HMASS   ,HWIDTH  ,W(1,         
     &     56  ))                                                          
      CALL JVSXXX(W(1,2   ),W(1,56  ),GWWH ,WMASS   ,WWIDTH  ,W(1,         
     &     57  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,57  ),GWF ,AMP(35  ))            
      CALL IOVXXX(W(1,22  ),W(1,6   ),W(1,52  ),GZN ,AMP(36  ))            
      CALL JWWWXX(W(1,1   ),W(1,39  ),W(1,2   ),GWWZ ,GWWA ,WMASS   ,      
     &     WWIDTH  ,W(1,58  ))                                             
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,58  ),GWF ,AMP(37  ))            
      ww(mu,nu)  = +AMP(   1)+AMP(   2)+AMP(   3)+AMP(   4)+AMP(   5)
     &             +AMP(   6)+AMP(   7)+AMP(   8)+AMP(   9)+AMP(  10)
     &             +AMP(  11)+AMP(  12)+AMP(  13)+AMP(  14)+AMP(  15)
     &             +AMP(  16)+AMP(  17)+AMP(  18)+AMP(  19)+AMP(  20)
     &             +AMP(  21)+AMP(  22)+AMP(  23)+AMP(  24)+AMP(  25)
     &             +AMP(  26)+AMP(  27)+AMP(  28)+AMP(  29)+AMP(  30)
     &             +AMP(  31)+AMP(  32)+AMP(  33)+AMP(  34)+AMP(  35)
     &             +AMP(  36)+AMP(  37)
     
c test only (Higgs graphs):     
c      ww(mu,nu)  = +AMP(24)+AMP(31)
     
      ENDDO
      enddo
      
      END
      

c----------------------------------------------------------------

      subroutine WPVtoWP_mg2_had(idl,id,Q,v,wpa,wpz)
C  Input:   Q(*)        one of the two incoming boson momenta (t-channel)
c           id = 1      Q = q1 = W+ momentum 
c           id = 2      Q = q2 = gamma/Z momentum
c	    idl = 1	l,vu = e+,ve
c	    idl = 2	l,vu = mu+,vmu	
c	
c  OUPTPUT: wpa(mu,nu)  lepton tensor for w gamma--> l+nu
c           wpz(mu,nu)  lepton tensor for w Z    --> l+nu
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : w+(q1) z(q2) -> l+ vl  
C FOR PROCESS : w+(q1) a(q2) -> l+ vl  
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
      integer id,idl
      real*8 q(0:3), v(0:3,4)
      complex*16 wpa(0:3,0:3), wpz(0:3,0:3)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J, MU, NU
      real*8 q2sq, Q1(0:3), Q2(0:3)
      COMPLEX*16 AMP(NGRAPHS)
      COMPLEX*16 W(6,NWAVEFUNCS), denz
      complex*16 wwp(6),wl(6),wn(6)
c  lepton spinors and W+- polarization vectors
#include "tensor.inc"
#include "tensorz.inc"

      logical ldebug
      parameter (ldebug=.false.)
C  
C GLOBAL VARIABLES
C  
#include "VBFNLO/utilities/coupl.inc"
C                                                                          
C ----------
C BEGIN CODE
C ----------

      if (idl.eq.1) then
         do i = 1,6
	    wwp(i) = we(i)
	    wl(i) = wel(i)
	    wn(i) = wve(i)
	 enddo      
      else
         do i = 1,6
	    wwp(i) = wu(i)
	    wl(i) = wmu(i)
	    wn(i) = wvu(i)
	 enddo            
      endif 
      
      if (id.eq.1) then
         q2(0) = dreal(wwp(5))
         q2(1) = dreal(wwp(6))
         q2(2) = dimag(wwp(6))
         q2(3) = dimag(wwp(5))
         do mu = 0,3
            q1(mu) = q(mu)
            q2(mu) = q2(mu)-q(mu)
         enddo
      elseif (id.eq.2) then
         q1(0) = dreal(wwp(5))
         q1(1) = dreal(wwp(6))
         q1(2) = dimag(wwp(6))
         q1(3) = dimag(wwp(5))
         do mu = 0,3
            q2(mu) = q(mu)
            q1(mu) = q1(mu)-q(mu)
         enddo
      endif
c
      do i=1,6                  ! read in lepton states from common block
         w(i,3) = wl(i)
         w(i,4) = wn(i)
	 w(i,7) = wwp(i)
      enddo

      DO MU = 0,3
         CALL VCARTX(Q1(0),WMASS ,Wwidth, MU,-1,W(1,1))        
         DO NU = 0,3

         if (idl.eq.1) then

            CALL VCARTX(Q2(0),ZMASS  ,Zwidth  , NU,-1,W(1,2))
            CALL FVIXXX(W(1,3   ),W(1,2   ),GZD ,ZERO    ,ZERO    ,W(1,5   ))    
            CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,1   ),GWF ,AMP(1   ))            
            CALL FVIXXX(W(1,3   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,6   ))    
            CALL IOVXXX(W(1,6   ),W(1,4   ),W(1,2   ),GZU ,AMP(2   ))            
            CALL JVVXXX(W(1,2   ),W(1,1   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,7   ))                                                          
            CALL IOVXXX(W(1,3   ),W(1,4   ),W(1,7   ),GWF ,AMP(3   ))            
            wpz(mu,nu) = (-AMP(   1)-AMP(   2)-AMP(   3))

            CALL VCARTX(Q2(0),ZERO  ,ZERO  , NU,-1,W(1,2))
            CALL FVIXXX(W(1,3   ),W(1,2   ),GAD ,ZERO    ,ZERO    ,W(1,5   ))    
            CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,1   ),GWF ,AMP(1   ))            
            CALL FVIXXX(W(1,3   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,6   ))    
            CALL IOVXXX(W(1,6   ),W(1,4   ),W(1,2   ),GAU ,AMP(2   ))            
            CALL JVVXXX(W(1,2   ),W(1,1   ),GWWA ,WMASS   ,WWIDTH  ,W(1,7   ))                                                          
            CALL IOVXXX(W(1,3   ),W(1,4   ),W(1,7   ),GWF ,AMP(3   ))            
            wpa(mu,nu) = -AMP(   1)-AMP(   2)-AMP(   3)

         else

            CALL VCARTX(Q2(0),ZMASS  ,Zwidth  , NU,-1,W(1,2))
            CALL FVIXXX(W(1,3   ),W(1,2   ),GZL ,ZERO    ,ZERO    ,W(1,5   ))    
            CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,1   ),GWF ,AMP(1   ))            
            CALL FVIXXX(W(1,3   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,6   ))    
            CALL IOVXXX(W(1,6   ),W(1,4   ),W(1,2   ),GZN ,AMP(2   ))            
            CALL JVVXXX(W(1,2   ),W(1,1   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,7   ))                                                          
            CALL IOVXXX(W(1,3   ),W(1,4   ),W(1,7   ),GWF ,AMP(3   ))            
            wpz(mu,nu) = (-AMP(   1)-AMP(   2)-AMP(   3))

            CALL VCARTX(Q2(0),ZERO  ,ZERO  , NU,-1,W(1,2))
            CALL FVIXXX(W(1,3   ),W(1,2   ),GAL ,ZERO    ,ZERO    ,W(1,5   ))    
            CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,1   ),GWF ,AMP(1   ))            
            CALL JVVXXX(W(1,2   ),W(1,1   ),GWWA ,WMASS   ,WWIDTH  ,W(1,6   ))                                                          
            CALL IOVXXX(W(1,3   ),W(1,4   ),W(1,6   ),GWF ,AMP(2   ))            
            wpa(mu,nu) = -AMP(   1)-AMP(   2)

         endif

	    
         enddo
      enddo
 
      end


