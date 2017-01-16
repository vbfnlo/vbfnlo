


      subroutine qqzzqq_mg2(p,pbar,v,sign, nlo, bos,
     1                  uucc,uuss,ddcc,ddss,udsc,ducs)

      implicit none                

#include "VBFNLO/utilities/global.inc"
#include "tensorz.inc"

#include "VBFNLO/utilities/process.inc"

      double precision pbar(0:3,4+max_v)
      double precision uucc,uuss,ddcc,ddss,udsc,ducs
      double precision uucc2,uuss2,ddcc2,ddss2,udsc2,ducs2
      double precision p(0:3,max_p),v(0:3,max_v)
      integer  sign(4+max_v), nlo, i, j,mu,nu, bos
      real*8 q12(0:4,3), q34(0:4,3)
      complex*16 mg2_tens_1d(0:5), mg2_tens_1d_2(0:5)
      complex*16 mg2_tens_2d_1(0:3,0:3)
      complex*16 mg2_tens_2d_2(0:3,0:3)
      complex*16 mg2_tens_2d_3(0:3,0:3)
      complex*16 mg2_tens_2d_4(0:3,0:3)
      complex*16 mg2_tens_2d_comb(0:3,0:3,2)
      double precision sumtest

      integer h
      common /hcount / h

      call Calc_Momentum_Transfer(p, v, q12,q34,1)

      j = 1
         
#ifdef WITH_MADGRAPH

         SELECT CASE(procid)
         CASE(ZZhadjj)

c leptonic tensors
         j = 1

           ! for AA,AZ,ZA,ZZ -> e+ e- mu+ mu- :
           call vvtozz_mg2_had(q12(0,j),q34(0,j),v,h,aazz(0,0,j),azzz(0,0,j),zazz(0,0,j),zzzz(0,0,j))
           ! for W+W-
           call wwtozz_mg2_had(q12(0,j),q34(0,j),v,h,wwzz5(0,0,j))  ! q12 = W+
           call wwtozz_mg2_had(q34(0,j),q12(0,j),v,h,wwzz6(0,0,j))  ! q12 = W-
           ! for V1V2 --> e+ e- or  V1V2 -->  mu+ mu- (Vi=A,Z)
           ! NCl tensors for NC process (k=1...4)
           call VVtoll_mg2_had(2,1,h,q34(0,j),v,aaee(0,0,1,j),azee(0,0,1,j),zaee(0,0,1,j),zzee(0,0,1,j))     !emit V on upper    
           call VVtoll_mg2_had(1,1,h,q12(0,j),v,aaee(0,0,2,j),azee(0,0,2,j),zaee(0,0,2,j),zzee(0,0,2,j))     !emit V on lower
           call VVtoll_mg2_had(2,2,h,q34(0,j),v,aauu(0,0,1,j),azuu(0,0,1,j),zauu(0,0,1,j),zzuu(0,0,1,j))     !emit V on upper     
           call VVtoll_mg2_had(1,2,h,q12(0,j),v,aauu(0,0,2,j),azuu(0,0,2,j),zauu(0,0,2,j),zzuu(0,0,2,j))     !emit V on lower
           ! for W+W- --> e+ e- or  W+W- -->  mu+ mu-
           ! CCl tensors for CC (k=5)
           call WWtoll_mg2_had(2,1,h,q34(0,j),v,CCee(0,0,1,j))     !emit V on upper
           call WWtoll_mg2_had(1,1,h,q12(0,j),v,CCee(0,0,2,j))     !emit V on lower
           call WWtoll_mg2_had(2,2,h,q34(0,j),v,CCuu(0,0,1,j))     !emit V on upper
           call WWtoll_mg2_had(1,2,h,q12(0,j),v,CCuu(0,0,2,j))     !emit V on lower
           ! for W-W+ --> e+ e- or mu+ mu-
           ! CCl for CC (k=6)
           call WWtoll_mg2_had(1,1,h,q34(0,j),v,CCee6(0,0,1,j))     !emit V on upper
           call WWtoll_mg2_had(2,1,h,q12(0,j),v,CCee6(0,0,2,j))     !emit V on lower
           call WWtoll_mg2_had(1,2,h,q34(0,j),v,CCuu6(0,0,1,j))     !emit V on upper
           call WWtoll_mg2_had(2,2,h,q12(0,j),v,CCuu6(0,0,2,j))     !emit V on lower
           ! for V-> llll
           call vto4l_mg2_had(v,h,azz,zzztens)


         call qqzzqq(pbar,sign, nlo, 1,bos,
     1                    uucc2,uuss2,ddcc2,ddss2,udsc2,ducs2)


           ! for AA,AZ,ZA,ZZ -> e+ e- mu+ mu- :
           call vvtozz_had(q12(0,j),q34(0,j),v,h,aazz(0,0,j),azzz(0,0,j),zazz(0,0,j),zzzz(0,0,j))
           ! for W+W-
           call wwtozz_had(q12(0,j),q34(0,j),v,h,wwzz5(0,0,j))  ! q12 = W+
           call wwtozz_had(q34(0,j),q12(0,j),v,h,wwzz6(0,0,j))  ! q12 = W-
           ! for V1V2 --> e+ e- or  V1V2 -->  mu+ mu- (Vi=A,Z)
           ! NCl tensors for NC process (k=1...4)
           call VVtoll_had(2,1,h,q34(0,j),v,aaee(0,0,1,j),azee(0,0,1,j),zaee(0,0,1,j),zzee(0,0,1,j))     !emit V on upper    
           call VVtoll_had(1,1,h,q12(0,j),v,aaee(0,0,2,j),azee(0,0,2,j),zaee(0,0,2,j),zzee(0,0,2,j))     !emit V on lower
           call VVtoll_had(2,2,h,q34(0,j),v,aauu(0,0,1,j),azuu(0,0,1,j),zauu(0,0,1,j),zzuu(0,0,1,j))     !emit V on upper     
           call VVtoll_had(1,2,h,q12(0,j),v,aauu(0,0,2,j),azuu(0,0,2,j),zauu(0,0,2,j),zzuu(0,0,2,j))     !emit V on lower
           ! for W+W- --> e+ e- or  W+W- -->  mu+ mu-
           ! CCl tensors for CC (k=5)
           call WWtoll_had(2,1,h,q34(0,j),v,CCee(0,0,1,j))     !emit V on upper
           call WWtoll_had(1,1,h,q12(0,j),v,CCee(0,0,2,j))     !emit V on lower
           call WWtoll_had(2,2,h,q34(0,j),v,CCuu(0,0,1,j))     !emit V on upper
           call WWtoll_had(1,2,h,q12(0,j),v,CCuu(0,0,2,j))     !emit V on lower
           ! for W-W+ --> e+ e- or mu+ mu-
           ! CCl for CC (k=6)
           call WWtoll_had(1,1,h,q34(0,j),v,CCee6(0,0,1,j))     !emit V on upper
           call WWtoll_had(2,1,h,q12(0,j),v,CCee6(0,0,2,j))     !emit V on lower
           call WWtoll_had(1,2,h,q34(0,j),v,CCuu6(0,0,1,j))     !emit V on upper
           call WWtoll_had(2,2,h,q12(0,j),v,CCuu6(0,0,2,j))     !emit V on lower
           ! for V-> llll
           call vto4l_had(v,h,azz,zzztens)


         do mu = 0,3
          do nu = 0,3
           mg2_tens_2d_1(mu,nu) = dcmplx(0.,0.)
           mg2_tens_2d_2(mu,nu) = dcmplx(0.,0.)
           mg2_tens_2d_3(mu,nu) = dcmplx(0.,0.)
           mg2_tens_2d_4(mu,nu) = dcmplx(0.,0.)
           mg2_tens_2d_comb(mu,nu,1) = dcmplx(0.,0.)
           mg2_tens_2d_comb(mu,nu,2) = dcmplx(0.,0.)
          enddo
         enddo
         do i = 0,5
           mg2_tens_1d(i) = dcmplx(0.,0.)
           mg2_tens_1d_2(i) = dcmplx(0.,0.)
         enddo


         print*, ""
         call vto4l_mg2_had(v,h,mg2_tens_1d,mg2_tens_1d_2)
         sumtest = 0
         do i = 0,5
           sumtest = sumtest + abs( azz(i) / mg2_tens_1d(i) - 1)
!           print*, i, azz(i), mg2_tens_1d(i), abs( azz(i) / mg2_tens_1d(i) - 1)
           mg2_tens_1d(i) = dcmplx(0.,0.)
         enddo
         print*, "tensor a to zz:", sumtest

!         print*, ""
         sumtest = 0
         do i = 0,5
           sumtest = sumtest + abs( zzztens(i) / mg2_tens_1d_2(i) - 1)
!           print*, i, zzztens(i), mg2_tens_1d_2(i), abs( zzztens(i) / mg2_tens_1d_2(i) - 1)
           mg2_tens_1d_2(i) = dcmplx(0.,0.)
         enddo
         print*, "tensor z to zz:", sumtest


!         print*, ""
         call wwtozz_mg2_had(q12(0,j),q34(0,j),v,h,mg2_tens_2d_1(0,0))  ! q12 = W+
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( wwzz5(mu,nu,1) / mg2_tens_2d_1(mu,nu) - 1)
!           print*, mu,nu, wwzz5(mu,nu,1), mg2_tens_2d_1(mu,nu), abs( wwzz5(mu,nu,1) / mg2_tens_2d_1(mu,nu) - 1)
           mg2_tens_2d_1(mu,nu) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor ww to zz 5:", sumtest


!         print*, ""
         call vvtozz_mg2_had(q12(0,j),q34(0,j),v,h,mg2_tens_2d_1(0,0),mg2_tens_2d_2(0,0),mg2_tens_2d_3(0,0),mg2_tens_2d_4(0,0))
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( aazz(mu,nu,1) / mg2_tens_2d_1(mu,nu) - 1)
!           print*, mu,nu, aazz(mu,nu,1), mg2_tens_2d_1(mu,nu), abs( aazz(mu,nu,1) / mg2_tens_2d_1(mu,nu) - 1)
           mg2_tens_2d_1(mu,nu) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor aa to zz:", sumtest

!         print*, ""
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( azzz(mu,nu,1) / mg2_tens_2d_2(mu,nu) - 1)
!           print*, mu,nu, azzz(mu,nu,1), mg2_tens_2d_2(mu,nu), abs( azzz(mu,nu,1) / mg2_tens_2d_2(mu,nu) - 1)
           mg2_tens_2d_2(mu,nu) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor az to zz:", sumtest

!         print*, ""
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( zazz(mu,nu,1) / mg2_tens_2d_3(mu,nu) - 1)
!           print*, mu,nu, zazz(mu,nu,1), mg2_tens_2d_3(mu,nu), abs( zazz(mu,nu,1) / mg2_tens_2d_3(mu,nu) - 1)
           mg2_tens_2d_3(mu,nu) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor za to zz:", sumtest

!         print*, ""
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( zzzz(mu,nu,1) / mg2_tens_2d_4(mu,nu) - 1)
!           print*, mu,nu, zzzz(mu,nu,1), mg2_tens_2d_4(mu,nu), abs( zzzz(mu,nu,1) / mg2_tens_2d_4(mu,nu) - 1)
           mg2_tens_2d_4(mu,nu) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor zz to zz:", sumtest


!         print*, ""
         call VVtoll_mg2_had(1,1,h,q12(0,j),v,mg2_tens_2d_1(0,0),mg2_tens_2d_2(0,0),mg2_tens_2d_3(0,0),mg2_tens_2d_4(0,0))     !emit V on lower
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( aaee(mu,nu,2,j) / mg2_tens_2d_1(mu,nu) - 1)
!           print*, mu,nu, aaee(mu,nu,2,j), mg2_tens_2d_1(mu,nu), abs( aaee(mu,nu,2,j) / mg2_tens_2d_1(mu,nu) - 1)
           mg2_tens_2d_1(mu,nu) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor aa to z:", sumtest

!         print*, ""
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( azee(mu,nu,2,j) / mg2_tens_2d_2(mu,nu) - 1)
!           print*, mu,nu, azee(mu,nu,2,j), mg2_tens_2d_2(mu,nu), abs( azee(mu,nu,2,j) / mg2_tens_2d_2(mu,nu) - 1)
           mg2_tens_2d_2(mu,nu) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor az to z:", sumtest

!         print*, ""
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( zaee(mu,nu,2,j) / mg2_tens_2d_3(mu,nu) - 1)
!           print*, mu,nu, zaee(mu,nu,2,j), mg2_tens_2d_3(mu,nu), abs( zaee(mu,nu,2,j) / mg2_tens_2d_3(mu,nu) - 1)
           mg2_tens_2d_3(mu,nu) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor za to z:", sumtest

!         print*, ""
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( zzee(mu,nu,2,j) / mg2_tens_2d_4(mu,nu) - 1)
!           print*, mu,nu, zzee(mu,nu,2,j), mg2_tens_2d_4(mu,nu), abs( zzee(mu,nu,2,j) / mg2_tens_2d_4(mu,nu) - 1)
           mg2_tens_2d_4(mu,nu) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor zz to z:", sumtest


!         print*, ""
         call WWtoll_mg2_had(1,1,h,q12(0,j),v,mg2_tens_2d_1(0,0))     !emit V on lower
         sumtest = 0
         do mu = 0,3
          do nu = 0,3
           sumtest = sumtest + abs( CCee(mu,nu,2,j) / mg2_tens_2d_1(mu,nu) - 1)
!           print*, mu,nu, CCee(mu,nu,2,j), mg2_tens_2d_1(mu,nu), abs( CCee(mu,nu,2,j) / mg2_tens_2d_1(mu,nu) - 1)
           mg2_tens_2d_1(mu,nu) = dcmplx(0.,0.)
          enddo
         enddo
         print*, "tensor ww to z:", sumtest


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


      subroutine wwtozz_mg2_had(q1,q2,P,h,ww)
	
      IMPLICIT NONE

      REAL*8 Q1(0:3), Q2(0:3), P(0:3,3:6)
      complex*16 ww(0:3,0:3)
      integer lsign(3:6,4),ll,mu,nu,h
      

	! possible lepton helicities:
	data (lsign(ll,1),ll=3,6) / -1, -1, -1, -1/  !h=1
	data (lsign(ll,2),ll=3,6) / -1, -1,  1,  1/  !h=2
	data (lsign(ll,3),ll=3,6) /  1,  1, -1, -1/  !h=3
	data (lsign(ll,4),ll=3,6) /  1,  1,  1,  1/  !h=4
	
	call WWtoZZ1_mg2_had(q1,q2,P,lsign(3,h),ww(0,0))
	
	! helicity combi nr.h is selected
	
	end
		
c ------------------------------------------------------------	

      subroutine WWtoZZ1_mg2_had(q1,q2,P,lsign,ww)
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : w+ w- -> e+ e- mu+ mu-  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=  50) 
      INTEGER    NWAVEFUNCS 
      PARAMETER (NWAVEFUNCS=  68) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 Q1(0:3), Q2(0:3), P(0:3,3:6)
      complex*16 ww(0:3,0:3)
      integer lsign(3:6)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J,MU, NU
      COMPLEX*16 AMP(NGRAPHS)
      COMPLEX*16 W(6,NWAVEFUNCS)
c
c  lepton spinors and Z polarization vectors
#include "tensorz.inc"
C  
C  
C GLOBAL VARIABLES
C  
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
   
      logical lswitch
      parameter (lswitch=.true.)
C                                                                          
C ----------
C BEGIN CODE
C ----------
  
      	do i = 1,6
	 w(i,3) = lep(i)
	 w(i,4) = lem(i)
	 w(i,5) = lup(i)
	 w(i,6) = lum(i)
      	enddo

      DO MU = 0,3
         DO NU = 0,3
      
      CALL VCARTX(Q1(0),WMASS ,wwidth, MU,-1,W(1,1   ))        
      CALL VCARTX(Q2(0),WMASS ,wwidth, NU,-1,W(1,2   )) 
 
      if (Ztouu) then


      CALL FVIXXX(W(1,3   ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,7   ))    
      CALL JIOXXX(W(1,7   ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,8   ))    
      CALL JVVXXX(W(1,1   ),W(1,8   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     9   ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,9   ),GAL ,AMP(1   ))            
      CALL JVVXXX(W(1,1   ),W(1,8   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     10  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,10  ),GZL ,AMP(2   ))            
      CALL FVIXXX(W(1,5   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,11  ))    
      CALL IOVXXX(W(1,11  ),W(1,6   ),W(1,8   ),GWF ,AMP(3   ))            
      CALL FVOXXX(W(1,4   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,12  ))    
      CALL JIOXXX(W(1,3   ),W(1,12  ),GWF ,WMASS   ,WWIDTH  ,W(1,13  ))    
      CALL FVIXXX(W(1,5   ),W(1,13  ),GWF ,ZERO    ,ZERO    ,W(1,14  ))    
      CALL IOVXXX(W(1,14  ),W(1,6   ),W(1,2   ),GWF ,AMP(4   ))            
      CALL JVVXXX(W(1,13  ),W(1,2   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     15  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,15  ),GAL ,AMP(5   ))            
      CALL JVVXXX(W(1,13  ),W(1,2   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     16  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,16  ),GZL ,AMP(6   ))            
      CALL FVOXXX(W(1,12  ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,17  ))    
      CALL JIOXXX(W(1,3   ),W(1,17  ),GAU ,ZERO    ,AWIDTH  ,W(1,18  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,18  ),GAL ,AMP(7   ))            
      CALL JIOXXX(W(1,3   ),W(1,17  ),GZU ,ZMASS   ,ZWIDTH  ,W(1,19  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,19  ),GZL ,AMP(8   ))            
      CALL FVIXXX(W(1,7   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,20  ))    
      CALL JIOXXX(W(1,20  ),W(1,4   ),GAU ,ZERO    ,AWIDTH  ,W(1,21  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,21  ),GAL ,AMP(9   ))            
      CALL JIOXXX(W(1,20  ),W(1,4   ),GZU ,ZMASS   ,ZWIDTH  ,W(1,22  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,22  ),GZL ,AMP(10  ))            
      CALL JIOXXX(W(1,7   ),W(1,12  ),GAD ,ZERO    ,AWIDTH  ,W(1,23  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,23  ),GAL ,AMP(11  ))            
      CALL JIOXXX(W(1,7   ),W(1,12  ),GZD ,ZMASS   ,ZWIDTH  ,W(1,24  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,24  ),GZL ,AMP(12  ))            
      CALL JIOXXX(W(1,3   ),W(1,4   ),GZU ,ZMASS   ,ZWIDTH  ,W(1,25  ))    
      CALL FVIXXX(W(1,11  ),W(1,25  ),GZN ,ZERO    ,ZERO    ,W(1,26  ))    
      CALL IOVXXX(W(1,26  ),W(1,6   ),W(1,2   ),GWF ,AMP(13  ))            
      CALL JIOXXX(W(1,3   ),W(1,4   ),GAU ,ZERO    ,AWIDTH  ,W(1,27  ))    
      CALL JVVXXX(W(1,2   ),W(1,27  ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     28  ))                                                          
      CALL JVVXXX(W(1,1   ),W(1,28  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     29  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,29  ),GAL ,AMP(14  ))            
      CALL JVVXXX(W(1,1   ),W(1,28  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     30  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,30  ),GZL ,AMP(15  ))            
      CALL JVVXXX(W(1,2   ),W(1,25  ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     31  ))                                                          
      CALL JVVXXX(W(1,1   ),W(1,31  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     32  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,32  ),GAL ,AMP(16  ))            
      CALL JVVXXX(W(1,1   ),W(1,31  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     33  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,33  ),GZL ,AMP(17  ))            
      CALL FVIXXX(W(1,11  ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,34  ))    
      CALL IOVXXX(W(1,34  ),W(1,6   ),W(1,27  ),GAL ,AMP(18  ))            
      CALL IOVXXX(W(1,34  ),W(1,6   ),W(1,25  ),GZL ,AMP(19  ))            
      CALL IOVXXX(W(1,11  ),W(1,6   ),W(1,28  ),GWF ,AMP(20  ))            
      CALL IOVXXX(W(1,11  ),W(1,6   ),W(1,31  ),GWF ,AMP(21  ))            
      CALL JVVXXX(W(1,27  ),W(1,1   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     35  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,35  ),GWF ,ZERO    ,ZERO    ,W(1,36  ))    
      CALL IOVXXX(W(1,36  ),W(1,6   ),W(1,2   ),GWF ,AMP(22  ))            
      CALL JVVXXX(W(1,25  ),W(1,1   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     37  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,37  ),GWF ,ZERO    ,ZERO    ,W(1,38  ))    
      CALL IOVXXX(W(1,38  ),W(1,6   ),W(1,2   ),GWF ,AMP(23  ))            
      CALL JVVXXX(W(1,35  ),W(1,2   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     39  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,39  ),GAL ,AMP(24  ))            
      CALL JVVXXX(W(1,35  ),W(1,2   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     40  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,40  ),GZL ,AMP(25  ))            
      CALL JVVXXX(W(1,37  ),W(1,2   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     41  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,41  ),GAL ,AMP(26  ))            
      CALL JVVXXX(W(1,37  ),W(1,2   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     42  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,42  ),GZL ,AMP(27  ))            
      CALL JVVXXX(W(1,1   ),W(1,2   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     43  ))                                                          
      CALL FVOXXX(W(1,4   ),W(1,43  ),GAU ,ZERO    ,ZERO    ,W(1,44  ))    
      CALL JIOXXX(W(1,3   ),W(1,44  ),GAU ,ZERO    ,AWIDTH  ,W(1,45  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,45  ),GAL ,AMP(28  ))            
      CALL JIOXXX(W(1,3   ),W(1,44  ),GZU ,ZMASS   ,ZWIDTH  ,W(1,46  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,46  ),GZL ,AMP(29  ))            
      CALL JVVXXX(W(1,1   ),W(1,2   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     47  ))                                                          
      CALL FVOXXX(W(1,4   ),W(1,47  ),GZU ,ZERO    ,ZERO    ,W(1,48  ))    
      CALL JIOXXX(W(1,3   ),W(1,48  ),GAU ,ZERO    ,AWIDTH  ,W(1,49  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,49  ),GAL ,AMP(30  ))            
      CALL JIOXXX(W(1,3   ),W(1,48  ),GZU ,ZMASS   ,ZWIDTH  ,W(1,50  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,50  ),GZL ,AMP(31  ))            
      CALL FVIXXX(W(1,3   ),W(1,43  ),GAU ,ZERO    ,ZERO    ,W(1,51  ))    
      CALL JIOXXX(W(1,51  ),W(1,4   ),GAU ,ZERO    ,AWIDTH  ,W(1,52  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,52  ),GAL ,AMP(32  ))            
      CALL JIOXXX(W(1,51  ),W(1,4   ),GZU ,ZMASS   ,ZWIDTH  ,W(1,53  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,53  ),GZL ,AMP(33  ))            
      CALL FVIXXX(W(1,3   ),W(1,47  ),GZU ,ZERO    ,ZERO    ,W(1,54  ))    
      CALL JIOXXX(W(1,54  ),W(1,4   ),GAU ,ZERO    ,AWIDTH  ,W(1,55  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,55  ),GAL ,AMP(34  ))            
      CALL JIOXXX(W(1,54  ),W(1,4   ),GZU ,ZMASS   ,ZWIDTH  ,W(1,56  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,56  ),GZL ,AMP(35  ))            
      CALL JW3WXX(W(1,1   ),W(1,27  ),W(1,2   ),GWWA ,GWWA ,ZERO    ,      
     &     AWIDTH  ,W(1,57  ))                                             
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,57  ),GAL ,AMP(36  ))            
      CALL JW3WXX(W(1,1   ),W(1,27  ),W(1,2   ),GWWA ,GWWZ ,ZMASS   ,      
     &     ZWIDTH  ,W(1,58  ))                                             
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,58  ),GZL ,AMP(37  ))            
      CALL JW3WXX(W(1,1   ),W(1,25  ),W(1,2   ),GWWZ ,GWWA ,ZERO    ,      
     &     AWIDTH  ,W(1,59  ))                                             
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,59  ),GAL ,AMP(38  ))            
      CALL JW3WXX(W(1,1   ),W(1,25  ),W(1,2   ),GWWZ ,GWWZ ,ZMASS   ,      
     &     ZWIDTH  ,W(1,60  ))                                             
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,60  ),GZL ,AMP(39  ))            
      CALL FVIXXX(W(1,5   ),W(1,27  ),GAL ,ZERO    ,ZERO    ,W(1,61  ))    
      CALL FVIXXX(W(1,61  ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,62  ))    
      CALL IOVXXX(W(1,62  ),W(1,6   ),W(1,2   ),GWF ,AMP(40  ))            
      CALL FVIXXX(W(1,5   ),W(1,25  ),GZL ,ZERO    ,ZERO    ,W(1,63  ))    
      CALL FVIXXX(W(1,63  ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,64  ))    
      CALL IOVXXX(W(1,64  ),W(1,6   ),W(1,2   ),GWF ,AMP(41  ))            
      CALL HVVXXX(W(1,1   ),W(1,2   ),GWWH ,HMASS   ,HWIDTH  ,W(1,         
     &     65  ))                                                          
      CALL JVSXXX(W(1,25  ),W(1,65  ),GZZH ,ZMASS   ,ZWIDTH  ,W(1,         
     &     66  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,66  ),GZL ,AMP(42  ))            
      CALL IOVXXX(W(1,61  ),W(1,6   ),W(1,43  ),GAL ,AMP(43  ))            
      CALL IOVXXX(W(1,63  ),W(1,6   ),W(1,43  ),GAL ,AMP(44  ))            
      CALL IOVXXX(W(1,61  ),W(1,6   ),W(1,47  ),GZL ,AMP(45  ))            
      CALL IOVXXX(W(1,63  ),W(1,6   ),W(1,47  ),GZL ,AMP(46  ))            
      CALL FVIXXX(W(1,5   ),W(1,43  ),GAL ,ZERO    ,ZERO    ,W(1,67  ))    
      CALL IOVXXX(W(1,67  ),W(1,6   ),W(1,27  ),GAL ,AMP(47  ))            
      CALL IOVXXX(W(1,67  ),W(1,6   ),W(1,25  ),GZL ,AMP(48  ))            
      CALL FVIXXX(W(1,5   ),W(1,47  ),GZL ,ZERO    ,ZERO    ,W(1,68  ))    
      CALL IOVXXX(W(1,68  ),W(1,6   ),W(1,27  ),GAL ,AMP(49  ))            
      CALL IOVXXX(W(1,68  ),W(1,6   ),W(1,25  ),GZL ,AMP(50  ))            
      WW(MU,NU) = +AMP(   1)+AMP(   2)+AMP(   3)+AMP(   4)+AMP(   5)
     &             +AMP(   6)+AMP(   7)+AMP(   8)+AMP(   9)+AMP(  10)
     &             +AMP(  11)+AMP(  12)+AMP(  13)+AMP(  14)+AMP(  15)
     &             +AMP(  16)+AMP(  17)+AMP(  18)+AMP(  19)+AMP(  20)
     &             +AMP(  21)+AMP(  22)+AMP(  23)+AMP(  24)+AMP(  25)
     &             +AMP(  26)+AMP(  27)+AMP(  28)+AMP(  29)+AMP(  30)
     &             +AMP(  31)+AMP(  32)+AMP(  33)+AMP(  34)+AMP(  35)
     &             +AMP(  36)+AMP(  37)+AMP(  38)+AMP(  39)+AMP(  40)
     &             +AMP(  41)+AMP(  42)+AMP(  43)+AMP(  44)+AMP(  45)
     &             +AMP(  46)+AMP(  47)+AMP(  48)+AMP(  49)+AMP(  50)

      elseif (Ztodd) then


      CALL FVOXXX(W(1,4   ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,7   ))    
      CALL JIOXXX(W(1,3   ),W(1,7   ),GWF ,WMASS   ,WWIDTH  ,W(1,8   ))    
      CALL JVVXXX(W(1,1   ),W(1,8   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     9   ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,9   ),GAL ,AMP(1   ))            
      CALL JVVXXX(W(1,1   ),W(1,8   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     10  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,10  ),GZL ,AMP(2   ))            
      CALL FVIXXX(W(1,5   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,11  ))    
      CALL IOVXXX(W(1,11  ),W(1,6   ),W(1,8   ),GWF ,AMP(3   ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,12  ))    
      CALL JIOXXX(W(1,12  ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,13  ))    
      CALL FVIXXX(W(1,5   ),W(1,13  ),GWF ,ZERO    ,ZERO    ,W(1,14  ))    
      CALL IOVXXX(W(1,14  ),W(1,6   ),W(1,2   ),GWF ,AMP(4   ))            
      CALL JVVXXX(W(1,13  ),W(1,2   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     15  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,15  ),GAL ,AMP(5   ))            
      CALL JVVXXX(W(1,13  ),W(1,2   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     16  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,16  ),GZL ,AMP(6   ))            
      CALL FVOXXX(W(1,7   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,17  ))    
      CALL JIOXXX(W(1,3   ),W(1,17  ),GAD ,ZERO    ,AWIDTH  ,W(1,18  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,18  ),GAL ,AMP(7   ))            
      CALL JIOXXX(W(1,3   ),W(1,17  ),GZD ,ZMASS   ,ZWIDTH  ,W(1,19  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,19  ),GZL ,AMP(8   ))            
      CALL FVIXXX(W(1,12  ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,20  ))    
      CALL JIOXXX(W(1,20  ),W(1,4   ),GAD ,ZERO    ,AWIDTH  ,W(1,21  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,21  ),GAL ,AMP(9   ))            
      CALL JIOXXX(W(1,20  ),W(1,4   ),GZD ,ZMASS   ,ZWIDTH  ,W(1,22  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,22  ),GZL ,AMP(10  ))            
      CALL JIOXXX(W(1,12  ),W(1,7   ),GAU ,ZERO    ,AWIDTH  ,W(1,23  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,23  ),GAL ,AMP(11  ))            
      CALL JIOXXX(W(1,12  ),W(1,7   ),GZU ,ZMASS   ,ZWIDTH  ,W(1,24  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,24  ),GZL ,AMP(12  ))            
      CALL JIOXXX(W(1,3   ),W(1,4   ),GZD ,ZMASS   ,ZWIDTH  ,W(1,25  ))    
      CALL FVIXXX(W(1,11  ),W(1,25  ),GZN ,ZERO    ,ZERO    ,W(1,26  ))    
      CALL IOVXXX(W(1,26  ),W(1,6   ),W(1,2   ),GWF ,AMP(13  ))            
      CALL JIOXXX(W(1,3   ),W(1,4   ),GAD ,ZERO    ,AWIDTH  ,W(1,27  ))    
      CALL JVVXXX(W(1,2   ),W(1,27  ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     28  ))                                                          
      CALL JVVXXX(W(1,1   ),W(1,28  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     29  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,29  ),GAL ,AMP(14  ))            
      CALL JVVXXX(W(1,1   ),W(1,28  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     30  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,30  ),GZL ,AMP(15  ))            
      CALL JVVXXX(W(1,2   ),W(1,25  ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     31  ))                                                          
      CALL JVVXXX(W(1,1   ),W(1,31  ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     32  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,32  ),GAL ,AMP(16  ))            
      CALL JVVXXX(W(1,1   ),W(1,31  ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     33  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,33  ),GZL ,AMP(17  ))            
      CALL FVIXXX(W(1,11  ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,34  ))    
      CALL IOVXXX(W(1,34  ),W(1,6   ),W(1,27  ),GAL ,AMP(18  ))            
      CALL IOVXXX(W(1,34  ),W(1,6   ),W(1,25  ),GZL ,AMP(19  ))            
      CALL IOVXXX(W(1,11  ),W(1,6   ),W(1,28  ),GWF ,AMP(20  ))            
      CALL IOVXXX(W(1,11  ),W(1,6   ),W(1,31  ),GWF ,AMP(21  ))            
      CALL JVVXXX(W(1,27  ),W(1,1   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     35  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,35  ),GWF ,ZERO    ,ZERO    ,W(1,36  ))    
      CALL IOVXXX(W(1,36  ),W(1,6   ),W(1,2   ),GWF ,AMP(22  ))            
      CALL JVVXXX(W(1,25  ),W(1,1   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     37  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,37  ),GWF ,ZERO    ,ZERO    ,W(1,38  ))    
      CALL IOVXXX(W(1,38  ),W(1,6   ),W(1,2   ),GWF ,AMP(23  ))            
      CALL JVVXXX(W(1,35  ),W(1,2   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     39  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,39  ),GAL ,AMP(24  ))            
      CALL JVVXXX(W(1,35  ),W(1,2   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     40  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,40  ),GZL ,AMP(25  ))            
      CALL JVVXXX(W(1,37  ),W(1,2   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     41  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,41  ),GAL ,AMP(26  ))            
      CALL JVVXXX(W(1,37  ),W(1,2   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     42  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,42  ),GZL ,AMP(27  ))            
      CALL JVVXXX(W(1,1   ),W(1,2   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     43  ))                                                          
      CALL FVOXXX(W(1,4   ),W(1,43  ),GAD ,ZERO    ,ZERO    ,W(1,44  ))    
      CALL JIOXXX(W(1,3   ),W(1,44  ),GAD ,ZERO    ,AWIDTH  ,W(1,45  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,45  ),GAL ,AMP(28  ))            
      CALL JIOXXX(W(1,3   ),W(1,44  ),GZD ,ZMASS   ,ZWIDTH  ,W(1,46  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,46  ),GZL ,AMP(29  ))            
      CALL JVVXXX(W(1,1   ),W(1,2   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     47  ))                                                          
      CALL FVOXXX(W(1,4   ),W(1,47  ),GZD ,ZERO    ,ZERO    ,W(1,48  ))    
      CALL JIOXXX(W(1,3   ),W(1,48  ),GAD ,ZERO    ,AWIDTH  ,W(1,49  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,49  ),GAL ,AMP(30  ))            
      CALL JIOXXX(W(1,3   ),W(1,48  ),GZD ,ZMASS   ,ZWIDTH  ,W(1,50  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,50  ),GZL ,AMP(31  ))            
      CALL FVIXXX(W(1,3   ),W(1,43  ),GAD ,ZERO    ,ZERO    ,W(1,51  ))    
      CALL JIOXXX(W(1,51  ),W(1,4   ),GAD ,ZERO    ,AWIDTH  ,W(1,52  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,52  ),GAL ,AMP(32  ))            
      CALL JIOXXX(W(1,51  ),W(1,4   ),GZD ,ZMASS   ,ZWIDTH  ,W(1,53  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,53  ),GZL ,AMP(33  ))            
      CALL FVIXXX(W(1,3   ),W(1,47  ),GZD ,ZERO    ,ZERO    ,W(1,54  ))    
      CALL JIOXXX(W(1,54  ),W(1,4   ),GAD ,ZERO    ,AWIDTH  ,W(1,55  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,55  ),GAL ,AMP(34  ))            
      CALL JIOXXX(W(1,54  ),W(1,4   ),GZD ,ZMASS   ,ZWIDTH  ,W(1,56  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,56  ),GZL ,AMP(35  ))            
      CALL JW3WXX(W(1,1   ),W(1,27  ),W(1,2   ),GWWA ,GWWA ,ZERO    ,      
     &     AWIDTH  ,W(1,57  ))                                             
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,57  ),GAL ,AMP(36  ))            
      CALL JW3WXX(W(1,1   ),W(1,27  ),W(1,2   ),GWWA ,GWWZ ,ZMASS   ,      
     &     ZWIDTH  ,W(1,58  ))                                             
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,58  ),GZL ,AMP(37  ))            
      CALL JW3WXX(W(1,1   ),W(1,25  ),W(1,2   ),GWWZ ,GWWA ,ZERO    ,      
     &     AWIDTH  ,W(1,59  ))                                             
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,59  ),GAL ,AMP(38  ))            
      CALL JW3WXX(W(1,1   ),W(1,25  ),W(1,2   ),GWWZ ,GWWZ ,ZMASS   ,      
     &     ZWIDTH  ,W(1,60  ))                                             
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,60  ),GZL ,AMP(39  ))            
      CALL FVIXXX(W(1,5   ),W(1,27  ),GAL ,ZERO    ,ZERO    ,W(1,61  ))    
      CALL FVIXXX(W(1,61  ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,62  ))    
      CALL IOVXXX(W(1,62  ),W(1,6   ),W(1,2   ),GWF ,AMP(40  ))            
      CALL FVIXXX(W(1,5   ),W(1,25  ),GZL ,ZERO    ,ZERO    ,W(1,63  ))    
      CALL FVIXXX(W(1,63  ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,64  ))    
      CALL IOVXXX(W(1,64  ),W(1,6   ),W(1,2   ),GWF ,AMP(41  ))            
      CALL HVVXXX(W(1,1   ),W(1,2   ),GWWH ,HMASS   ,HWIDTH  ,W(1,         
     &     65  ))                                                          
      CALL JVSXXX(W(1,25  ),W(1,65  ),GZZH ,ZMASS   ,ZWIDTH  ,W(1,         
     &     66  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,66  ),GZL ,AMP(42  ))            
      CALL IOVXXX(W(1,61  ),W(1,6   ),W(1,43  ),GAL ,AMP(43  ))            
      CALL IOVXXX(W(1,63  ),W(1,6   ),W(1,43  ),GAL ,AMP(44  ))            
      CALL IOVXXX(W(1,61  ),W(1,6   ),W(1,47  ),GZL ,AMP(45  ))            
      CALL IOVXXX(W(1,63  ),W(1,6   ),W(1,47  ),GZL ,AMP(46  ))            
      CALL FVIXXX(W(1,5   ),W(1,43  ),GAL ,ZERO    ,ZERO    ,W(1,67  ))    
      CALL IOVXXX(W(1,67  ),W(1,6   ),W(1,27  ),GAL ,AMP(47  ))            
      CALL IOVXXX(W(1,67  ),W(1,6   ),W(1,25  ),GZL ,AMP(48  ))            
      CALL FVIXXX(W(1,5   ),W(1,47  ),GZL ,ZERO    ,ZERO    ,W(1,68  ))    
      CALL IOVXXX(W(1,68  ),W(1,6   ),W(1,27  ),GAL ,AMP(49  ))            
      CALL IOVXXX(W(1,68  ),W(1,6   ),W(1,25  ),GZL ,AMP(50  ))            
       WW(MU,NU) = +AMP(   1)+AMP(   2)+AMP(   3)+AMP(   4)+AMP(   5)
     &             +AMP(   6)+AMP(   7)+AMP(   8)+AMP(   9)+AMP(  10)
     &             +AMP(  11)+AMP(  12)+AMP(  13)+AMP(  14)+AMP(  15)
     &             +AMP(  16)+AMP(  17)+AMP(  18)+AMP(  19)+AMP(  20)
     &             +AMP(  21)+AMP(  22)+AMP(  23)+AMP(  24)+AMP(  25)
     &             +AMP(  26)+AMP(  27)+AMP(  28)+AMP(  29)+AMP(  30)
     &             +AMP(  31)+AMP(  32)+AMP(  33)+AMP(  34)+AMP(  35)
     &             +AMP(  36)+AMP(  37)+AMP(  38)+AMP(  39)+AMP(  40)
     &             +AMP(  41)+AMP(  42)+AMP(  43)+AMP(  44)+AMP(  45)
     &             +AMP(  46)+AMP(  47)+AMP(  48)+AMP(  49)+AMP(  50)


      endif

            ENDDO
      ENDDO

      END
c ------------------------------------------------------------	

c ------------------------------------------------------------	

      subroutine wwtoll_mg2_had(idb,idl,h,Q,v,wll)
	
      IMPLICIT NONE

        integer idb,idl,h
        real*8 q(0:3), v(0:3,4)
        complex*16 wll(0:3,0:3)
	integer lsign(2,2),i,ll,mu,nu
      

	data (lsign(ll,1),ll=1,2) / -1, -1/
	data (lsign(ll,2),ll=1,2) / +1, +1/
	
	if (idl.eq.1) then ! e+e-
	   if (h.eq.1.or.h.eq.2)  i = 1
	   if (h.eq.3.or.h.eq.4)  i = 2
	elseif (idl.eq.2) then ! mu+mu-
	   if (h.eq.1.or.h.eq.3)  i = 1
	   if (h.eq.2.or.h.eq.4)  i = 2
	endif   

	
	call WWtoLL1_mg2_had(idb,idl,Q,v,lsign(1,i),wll(0,0))


	end
		
c ------------------------------------------------------------	
       
      subroutine WWtoLL1_mg2_had(idb,idl,Q,v,lsign,wll)
C  Input:   Q(*)        one of the two incoming boson momenta 
c           idb = 1      Q = q1 = W+ momentum
c           idb = 2      Q = q2 = W- momentum
c	    idl = 1	 l is e
c	    idl = 2	 l is mu
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
      PARAMETER (NWAVEFUNCS=   7) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      integer idb,idl
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
#include "tensorz.inc"
C  
C GLOBAL VARIABLES
C  
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
   
      logical lswitch
      parameter (lswitch=.true.)
C  
C ----------
C BEGIN CODE
C ----------
      
      if (idb.eq.1) then
         do mu = 0,3
            q1(mu) = q(mu)
	    if (idl.eq.1) then
                q2(mu) = qe(mu)-q(mu)
		p(mu,3) = v(mu,1)
		p(mu,4) = v(mu,2)
	    else
	        q2(mu) = qu(mu)-q(mu)
		p(mu,3) = v(mu,3)
		p(mu,4) = v(mu,4)
	    endif
         enddo
      elseif (idb.eq.2) then
         do mu = 0,3
            q2(mu) = q(mu)
	    if (idl.eq.1) then
            	q1(mu) = qe(mu)-q(mu)
		p(mu,3) = v(mu,1)
		p(mu,4) = v(mu,2)
	    else
            	q1(mu) = qu(mu)-q(mu)
		p(mu,3) = v(mu,3)
		p(mu,4) = v(mu,4)
	    endif
         enddo
      endif

	  if (idl.eq.1) then
	    do i = 1,6
	    w(i,3) = lep(i)
	    w(i,4) = lem(i)
	    enddo
	  else
	    do i = 1,6
	    w(i,3) = lup(i)
	    w(i,4) = lum(i)
	    enddo
	  endif
	     
      DO MU = 0,3
         DO NU = 0,3
      
      CALL VCARTX(Q1(0),WMASS ,wwidth, MU,-1,W(1,1   ))        
      CALL VCARTX(Q2(0),WMASS ,wwidth, NU,-1,W(1,2   )) 

      CALL FVIXXX(W(1,3   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,5   ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,2   ),GWF ,AMP(1   ))            
      CALL JVVXXX(W(1,1   ),W(1,2   ),GWWA ,ZERO    ,ZERO    ,W(1,         
     &     6   ))                                                          
      CALL IOVXXX(W(1,3   ),W(1,4   ),W(1,6   ),GAL ,AMP(2   ))            
      CALL JVVXXX(W(1,1   ),W(1,2   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     7   ))                                                          
      CALL IOVXXX(W(1,3   ),W(1,4   ),W(1,7   ),GZL ,AMP(3   ))            
      
      
      wll(mu,nu) = +AMP(   1)+AMP(   2)+AMP(   3)


      if (Ztouu .and. idl.eq.1) then


      CALL FVIXXX(W(1,3   ),W(1,2   ),GWF ,ZERO    ,ZERO    ,W(1,5   ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,1   ),GWF ,AMP(1   ))            
      CALL JVVXXX(W(1,1   ),W(1,2   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     6   ))                                                          
      CALL IOVXXX(W(1,3   ),W(1,4   ),W(1,6   ),GAU ,AMP(2   ))            
      CALL JVVXXX(W(1,1   ),W(1,2   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     7   ))                                                          
      CALL IOVXXX(W(1,3   ),W(1,4   ),W(1,7   ),GZU ,AMP(3   ))            
      wll(MU,NU) = +AMP(   1)+AMP(   2)+AMP(   3)
      

      elseif (Ztodd .and. idl.eq.1) then


      CALL FVIXXX(W(1,3   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,5   ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,2   ),GWF ,AMP(1   ))            
      CALL JVVXXX(W(1,1   ),W(1,2   ),GWWA ,ZERO    ,AWIDTH  ,W(1,         
     &     6   ))                                                          
      CALL IOVXXX(W(1,3   ),W(1,4   ),W(1,6   ),GAD ,AMP(2   ))            
      CALL JVVXXX(W(1,1   ),W(1,2   ),GWWZ ,ZMASS   ,ZWIDTH  ,W(1,         
     &     7   ))                                                          
      CALL IOVXXX(W(1,3   ),W(1,4   ),W(1,7   ),GZD ,AMP(3   ))            
      wll(MU,NU) = +AMP(   1)+AMP(   2)+AMP(   3)


      endif
         
	 enddo !nu
      enddo !mu

      end
      
c ------------------------------------------------------------


c ------------------------------------------------------------	

      subroutine vvtozz_mg2_had(q1,q2,P,h,aa,az,za,zz)
	
      IMPLICIT NONE

      REAL*8 Q1(0:3), Q2(0:3), P(0:3,3:6)
      complex*16 aa(0:3,0:3)
      complex*16 az(0:3,0:3)
      complex*16 za(0:3,0:3)
      complex*16 zz(0:3,0:3)
      integer lsign(3:6,4),h,ll,mu,nu
      

	! possible lepton helicities:
	data (lsign(ll,1),ll=3,6) / -1, -1, -1, -1/
	data (lsign(ll,2),ll=3,6) / -1, -1,  1,  1/
	data (lsign(ll,3),ll=3,6) /  1,  1, -1, -1/
	data (lsign(ll,4),ll=3,6) /  1,  1,  1,  1/
	
	call AAtoZZ1_mg2_had(q1,q2,P,lsign(3,h),aa(0,0))
	call AZtoZZ1_mg2_had(q1,q2,P,lsign(3,h),az(0,0))
	call AZtoZZ1_mg2_had(q2,q1,P,lsign(3,h),za(0,0))
	call ZZtoZZ1_mg2_had(q1,q2,P,lsign(3,h),zz(0,0))
		
	end
		
c ------------------------------------------------------------	
       
      subroutine AAtoZZ1_mg2_had(q1,q2,P,lsign,aa)
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : a a -> e+ e- mu+ mu-  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=  40) 
      INTEGER    NWAVEFUNCS
      PARAMETER (NWAVEFUNCS=  58) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 Q1(0:3), Q2(0:3), P(0:3,3:6)
      complex*16 aa(0:3,0:3)
      integer lsign(3:6)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J, mu,nu
      COMPLEX*16 AMP(NGRAPHS)
      COMPLEX*16 W(6,NWAVEFUNCS)
C  
c  lepton spinors and Z polarization vectors
#include "tensorz.inc"
c
C GLOBAL VARIABLES
C  
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/coupl.inc"
   
      logical lswitch
      parameter (lswitch=.true.)
C                                                                          
C ----------
C BEGIN CODE
C ----------
     
      	do i = 1,6
	 w(i,3) = lep(i)
	 w(i,4) = lem(i)
	 w(i,5) = lup(i)
	 w(i,6) = lum(i)
      	enddo

      DO MU = 0,3
         DO NU = 0,3
      
      CALL VCARTX(Q1(0),zero ,zero, MU,-1,W(1,1	  ))	    
      CALL VCARTX(Q2(0),zero ,zero, NU,-1,W(1,2	  ))      
      
      if (Ztouu) then


      CALL FVIXXX(W(1,3   ),W(1,2   ),GAU ,ZERO    ,ZERO    ,W(1,7   ))    
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
      CALL JIOXXX(W(1,3   ),W(1,13  ),GAU ,ZERO    ,AWIDTH  ,W(1,14  ))    
      CALL FVIXXX(W(1,5   ),W(1,14  ),GAL ,ZERO    ,ZERO    ,W(1,15  ))    
      CALL IOVXXX(W(1,15  ),W(1,6   ),W(1,2   ),GAL ,AMP(5   ))            
      CALL JIOXXX(W(1,3   ),W(1,13  ),GZU ,ZMASS   ,ZWIDTH  ,W(1,16  ))    
      CALL FVIXXX(W(1,5   ),W(1,16  ),GZL ,ZERO    ,ZERO    ,W(1,17  ))    
      CALL IOVXXX(W(1,17  ),W(1,6   ),W(1,2   ),GAL ,AMP(6   ))            
      CALL FVIXXX(W(1,5   ),W(1,2   ),GAL ,ZERO    ,ZERO    ,W(1,18  ))    
      CALL IOVXXX(W(1,18  ),W(1,6   ),W(1,14  ),GAL ,AMP(7   ))            
      CALL IOVXXX(W(1,18  ),W(1,6   ),W(1,16  ),GZL ,AMP(8   ))            
      CALL FVOXXX(W(1,13  ),W(1,2   ),GAU ,ZERO    ,ZERO    ,W(1,19  ))    
      CALL JIOXXX(W(1,3   ),W(1,19  ),GAU ,ZERO    ,AWIDTH  ,W(1,20  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,20  ),GAL ,AMP(9   ))            
      CALL JIOXXX(W(1,3   ),W(1,19  ),GZU ,ZMASS   ,ZWIDTH  ,W(1,21  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,21  ),GZL ,AMP(10  ))            
      CALL FVIXXX(W(1,7   ),W(1,1   ),GAU ,ZERO    ,ZERO    ,W(1,22  ))    
      CALL JIOXXX(W(1,22  ),W(1,4   ),GAU ,ZERO    ,AWIDTH  ,W(1,23  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,23  ),GAL ,AMP(11  ))            
      CALL JIOXXX(W(1,22  ),W(1,4   ),GZU ,ZMASS   ,ZWIDTH  ,W(1,24  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,24  ),GZL ,AMP(12  ))            
      CALL JIOXXX(W(1,7   ),W(1,13  ),GAU ,ZERO    ,AWIDTH  ,W(1,25  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,25  ),GAL ,AMP(13  ))            
      CALL JIOXXX(W(1,7   ),W(1,13  ),GZU ,ZMASS   ,ZWIDTH  ,W(1,26  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,26  ),GZL ,AMP(14  ))            
      CALL FVOXXX(W(1,4   ),W(1,2   ),GAU ,ZERO    ,ZERO    ,W(1,27  ))    
      CALL JIOXXX(W(1,3   ),W(1,27  ),GAU ,ZERO    ,AWIDTH  ,W(1,28  ))    
      CALL FVIXXX(W(1,5   ),W(1,28  ),GAL ,ZERO    ,ZERO    ,W(1,29  ))    
      CALL IOVXXX(W(1,29  ),W(1,6   ),W(1,1   ),GAL ,AMP(15  ))            
      CALL JIOXXX(W(1,3   ),W(1,27  ),GZU ,ZMASS   ,ZWIDTH  ,W(1,30  ))    
      CALL FVIXXX(W(1,5   ),W(1,30  ),GZL ,ZERO    ,ZERO    ,W(1,31  ))    
      CALL IOVXXX(W(1,31  ),W(1,6   ),W(1,1   ),GAL ,AMP(16  ))            
      CALL IOVXXX(W(1,12  ),W(1,6   ),W(1,28  ),GAL ,AMP(17  ))            
      CALL IOVXXX(W(1,12  ),W(1,6   ),W(1,30  ),GZL ,AMP(18  ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GAU ,ZERO    ,ZERO    ,W(1,32  ))    
      CALL JIOXXX(W(1,32  ),W(1,4   ),GAU ,ZERO    ,AWIDTH  ,W(1,33  ))    
      CALL FVIXXX(W(1,5   ),W(1,33  ),GAL ,ZERO    ,ZERO    ,W(1,34  ))    
      CALL IOVXXX(W(1,34  ),W(1,6   ),W(1,2   ),GAL ,AMP(19  ))            
      CALL JIOXXX(W(1,32  ),W(1,4   ),GZU ,ZMASS   ,ZWIDTH  ,W(1,35  ))    
      CALL FVIXXX(W(1,5   ),W(1,35  ),GZL ,ZERO    ,ZERO    ,W(1,36  ))    
      CALL IOVXXX(W(1,36  ),W(1,6   ),W(1,2   ),GAL ,AMP(20  ))            
      CALL IOVXXX(W(1,18  ),W(1,6   ),W(1,33  ),GAL ,AMP(21  ))            
      CALL IOVXXX(W(1,18  ),W(1,6   ),W(1,35  ),GZL ,AMP(22  ))            
      CALL FVOXXX(W(1,27  ),W(1,1   ),GAU ,ZERO    ,ZERO    ,W(1,37  ))    
      CALL JIOXXX(W(1,3   ),W(1,37  ),GAU ,ZERO    ,AWIDTH  ,W(1,38  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,38  ),GAL ,AMP(23  ))            
      CALL JIOXXX(W(1,3   ),W(1,37  ),GZU ,ZMASS   ,ZWIDTH  ,W(1,39  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,39  ),GZL ,AMP(24  ))            
      CALL FVIXXX(W(1,32  ),W(1,2   ),GAU ,ZERO    ,ZERO    ,W(1,40  ))    
      CALL JIOXXX(W(1,40  ),W(1,4   ),GAU ,ZERO    ,AWIDTH  ,W(1,41  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,41  ),GAL ,AMP(25  ))            
      CALL JIOXXX(W(1,40  ),W(1,4   ),GZU ,ZMASS   ,ZWIDTH  ,W(1,42  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,42  ),GZL ,AMP(26  ))            
      CALL JIOXXX(W(1,32  ),W(1,27  ),GAU ,ZERO    ,AWIDTH  ,W(1,43  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,43  ),GAL ,AMP(27  ))            
      CALL JIOXXX(W(1,32  ),W(1,27  ),GZU ,ZMASS   ,ZWIDTH  ,W(1,44  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,44  ),GZL ,AMP(28  ))            
      CALL JIOXXX(W(1,3   ),W(1,4   ),GAU ,ZERO    ,AWIDTH  ,W(1,45  ))    
      CALL FVIXXX(W(1,12  ),W(1,45  ),GAL ,ZERO    ,ZERO    ,W(1,46  ))    
      CALL IOVXXX(W(1,46  ),W(1,6   ),W(1,2   ),GAL ,AMP(29  ))            
      CALL JIOXXX(W(1,3   ),W(1,4   ),GZU ,ZMASS   ,ZWIDTH  ,W(1,47  ))    
      CALL FVIXXX(W(1,12  ),W(1,47  ),GZL ,ZERO    ,ZERO    ,W(1,48  ))    
      CALL IOVXXX(W(1,48  ),W(1,6   ),W(1,2   ),GAL ,AMP(30  ))            
      CALL FVIXXX(W(1,12  ),W(1,2   ),GAL ,ZERO    ,ZERO    ,W(1,49  ))    
      CALL IOVXXX(W(1,49  ),W(1,6   ),W(1,45  ),GAL ,AMP(31  ))            
      CALL IOVXXX(W(1,49  ),W(1,6   ),W(1,47  ),GZL ,AMP(32  ))            
      CALL FVIXXX(W(1,18  ),W(1,45  ),GAL ,ZERO    ,ZERO    ,W(1,50  ))    
      CALL IOVXXX(W(1,50  ),W(1,6   ),W(1,1   ),GAL ,AMP(33  ))            
      CALL FVIXXX(W(1,18  ),W(1,47  ),GZL ,ZERO    ,ZERO    ,W(1,51  ))    
      CALL IOVXXX(W(1,51  ),W(1,6   ),W(1,1   ),GAL ,AMP(34  ))            
      CALL FVIXXX(W(1,18  ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,52  ))    
      CALL IOVXXX(W(1,52  ),W(1,6   ),W(1,45  ),GAL ,AMP(35  ))            
      CALL IOVXXX(W(1,52  ),W(1,6   ),W(1,47  ),GZL ,AMP(36  ))            
      CALL FVIXXX(W(1,5   ),W(1,45  ),GAL ,ZERO    ,ZERO    ,W(1,53  ))    
      CALL FVIXXX(W(1,53  ),W(1,2   ),GAL ,ZERO    ,ZERO    ,W(1,54  ))    
      CALL IOVXXX(W(1,54  ),W(1,6   ),W(1,1   ),GAL ,AMP(37  ))            
      CALL FVIXXX(W(1,5   ),W(1,47  ),GZL ,ZERO    ,ZERO    ,W(1,55  ))    
      CALL FVIXXX(W(1,55  ),W(1,2   ),GAL ,ZERO    ,ZERO    ,W(1,56  ))    
      CALL IOVXXX(W(1,56  ),W(1,6   ),W(1,1   ),GAL ,AMP(38  ))            
      CALL FVIXXX(W(1,53  ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,57  ))    
      CALL IOVXXX(W(1,57  ),W(1,6   ),W(1,2   ),GAL ,AMP(39  ))            
      CALL FVIXXX(W(1,55  ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,58  ))    
      CALL IOVXXX(W(1,58  ),W(1,6   ),W(1,2   ),GAL ,AMP(40  ))            
      AA(MU,NU) = +AMP(   1)+AMP(   2)+AMP(   3)+AMP(   4)+AMP(   5)
     &             +AMP(   6)+AMP(   7)+AMP(   8)+AMP(   9)+AMP(  10)
     &             +AMP(  11)+AMP(  12)+AMP(  13)+AMP(  14)+AMP(  15)
     &             +AMP(  16)+AMP(  17)+AMP(  18)+AMP(  19)+AMP(  20)
     &             +AMP(  21)+AMP(  22)+AMP(  23)+AMP(  24)+AMP(  25)
     &             +AMP(  26)+AMP(  27)+AMP(  28)+AMP(  29)+AMP(  30)
     &             +AMP(  31)+AMP(  32)+AMP(  33)+AMP(  34)+AMP(  35)
     &             +AMP(  36)+AMP(  37)+AMP(  38)+AMP(  39)+AMP(  40)


      elseif (Ztodd) then


      CALL FVIXXX(W(1,3   ),W(1,2   ),GAD ,ZERO    ,ZERO    ,W(1,7   ))    
      CALL JIOXXX(W(1,7   ),W(1,4   ),GAD ,ZERO    ,AWIDTH  ,W(1,8   ))    
      CALL FVIXXX(W(1,5   ),W(1,8   ),GAL ,ZERO    ,ZERO    ,W(1,9   ))    
      CALL IOVXXX(W(1,9   ),W(1,6   ),W(1,1   ),GAL ,AMP(1   ))            
      CALL JIOXXX(W(1,7   ),W(1,4   ),GZD ,ZMASS   ,ZWIDTH  ,W(1,10  ))    
      CALL FVIXXX(W(1,5   ),W(1,10  ),GZL ,ZERO    ,ZERO    ,W(1,11  ))    
      CALL IOVXXX(W(1,11  ),W(1,6   ),W(1,1   ),GAL ,AMP(2   ))            
      CALL FVIXXX(W(1,5   ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,12  ))    
      CALL IOVXXX(W(1,12  ),W(1,6   ),W(1,8   ),GAL ,AMP(3   ))            
      CALL IOVXXX(W(1,12  ),W(1,6   ),W(1,10  ),GZL ,AMP(4   ))            
      CALL FVOXXX(W(1,4   ),W(1,1   ),GAD ,ZERO    ,ZERO    ,W(1,13  ))    
      CALL JIOXXX(W(1,3   ),W(1,13  ),GAD ,ZERO    ,AWIDTH  ,W(1,14  ))    
      CALL FVIXXX(W(1,5   ),W(1,14  ),GAL ,ZERO    ,ZERO    ,W(1,15  ))    
      CALL IOVXXX(W(1,15  ),W(1,6   ),W(1,2   ),GAL ,AMP(5   ))            
      CALL JIOXXX(W(1,3   ),W(1,13  ),GZD ,ZMASS   ,ZWIDTH  ,W(1,16  ))    
      CALL FVIXXX(W(1,5   ),W(1,16  ),GZL ,ZERO    ,ZERO    ,W(1,17  ))    
      CALL IOVXXX(W(1,17  ),W(1,6   ),W(1,2   ),GAL ,AMP(6   ))            
      CALL FVIXXX(W(1,5   ),W(1,2   ),GAL ,ZERO    ,ZERO    ,W(1,18  ))    
      CALL IOVXXX(W(1,18  ),W(1,6   ),W(1,14  ),GAL ,AMP(7   ))            
      CALL IOVXXX(W(1,18  ),W(1,6   ),W(1,16  ),GZL ,AMP(8   ))            
      CALL FVOXXX(W(1,13  ),W(1,2   ),GAD ,ZERO    ,ZERO    ,W(1,19  ))    
      CALL JIOXXX(W(1,3   ),W(1,19  ),GAD ,ZERO    ,AWIDTH  ,W(1,20  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,20  ),GAL ,AMP(9   ))            
      CALL JIOXXX(W(1,3   ),W(1,19  ),GZD ,ZMASS   ,ZWIDTH  ,W(1,21  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,21  ),GZL ,AMP(10  ))            
      CALL FVIXXX(W(1,7   ),W(1,1   ),GAD ,ZERO    ,ZERO    ,W(1,22  ))    
      CALL JIOXXX(W(1,22  ),W(1,4   ),GAD ,ZERO    ,AWIDTH  ,W(1,23  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,23  ),GAL ,AMP(11  ))            
      CALL JIOXXX(W(1,22  ),W(1,4   ),GZD ,ZMASS   ,ZWIDTH  ,W(1,24  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,24  ),GZL ,AMP(12  ))            
      CALL JIOXXX(W(1,7   ),W(1,13  ),GAD ,ZERO    ,AWIDTH  ,W(1,25  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,25  ),GAL ,AMP(13  ))            
      CALL JIOXXX(W(1,7   ),W(1,13  ),GZD ,ZMASS   ,ZWIDTH  ,W(1,26  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,26  ),GZL ,AMP(14  ))            
      CALL FVOXXX(W(1,4   ),W(1,2   ),GAD ,ZERO    ,ZERO    ,W(1,27  ))    
      CALL JIOXXX(W(1,3   ),W(1,27  ),GAD ,ZERO    ,AWIDTH  ,W(1,28  ))    
      CALL FVIXXX(W(1,5   ),W(1,28  ),GAL ,ZERO    ,ZERO    ,W(1,29  ))    
      CALL IOVXXX(W(1,29  ),W(1,6   ),W(1,1   ),GAL ,AMP(15  ))            
      CALL JIOXXX(W(1,3   ),W(1,27  ),GZD ,ZMASS   ,ZWIDTH  ,W(1,30  ))    
      CALL FVIXXX(W(1,5   ),W(1,30  ),GZL ,ZERO    ,ZERO    ,W(1,31  ))    
      CALL IOVXXX(W(1,31  ),W(1,6   ),W(1,1   ),GAL ,AMP(16  ))            
      CALL IOVXXX(W(1,12  ),W(1,6   ),W(1,28  ),GAL ,AMP(17  ))            
      CALL IOVXXX(W(1,12  ),W(1,6   ),W(1,30  ),GZL ,AMP(18  ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GAD ,ZERO    ,ZERO    ,W(1,32  ))    
      CALL JIOXXX(W(1,32  ),W(1,4   ),GAD ,ZERO    ,AWIDTH  ,W(1,33  ))    
      CALL FVIXXX(W(1,5   ),W(1,33  ),GAL ,ZERO    ,ZERO    ,W(1,34  ))    
      CALL IOVXXX(W(1,34  ),W(1,6   ),W(1,2   ),GAL ,AMP(19  ))            
      CALL JIOXXX(W(1,32  ),W(1,4   ),GZD ,ZMASS   ,ZWIDTH  ,W(1,35  ))    
      CALL FVIXXX(W(1,5   ),W(1,35  ),GZL ,ZERO    ,ZERO    ,W(1,36  ))    
      CALL IOVXXX(W(1,36  ),W(1,6   ),W(1,2   ),GAL ,AMP(20  ))            
      CALL IOVXXX(W(1,18  ),W(1,6   ),W(1,33  ),GAL ,AMP(21  ))            
      CALL IOVXXX(W(1,18  ),W(1,6   ),W(1,35  ),GZL ,AMP(22  ))            
      CALL FVOXXX(W(1,27  ),W(1,1   ),GAD ,ZERO    ,ZERO    ,W(1,37  ))    
      CALL JIOXXX(W(1,3   ),W(1,37  ),GAD ,ZERO    ,AWIDTH  ,W(1,38  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,38  ),GAL ,AMP(23  ))            
      CALL JIOXXX(W(1,3   ),W(1,37  ),GZD ,ZMASS   ,ZWIDTH  ,W(1,39  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,39  ),GZL ,AMP(24  ))            
      CALL FVIXXX(W(1,32  ),W(1,2   ),GAD ,ZERO    ,ZERO    ,W(1,40  ))    
      CALL JIOXXX(W(1,40  ),W(1,4   ),GAD ,ZERO    ,AWIDTH  ,W(1,41  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,41  ),GAL ,AMP(25  ))            
      CALL JIOXXX(W(1,40  ),W(1,4   ),GZD ,ZMASS   ,ZWIDTH  ,W(1,42  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,42  ),GZL ,AMP(26  ))            
      CALL JIOXXX(W(1,32  ),W(1,27  ),GAD ,ZERO    ,AWIDTH  ,W(1,43  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,43  ),GAL ,AMP(27  ))            
      CALL JIOXXX(W(1,32  ),W(1,27  ),GZD ,ZMASS   ,ZWIDTH  ,W(1,44  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,44  ),GZL ,AMP(28  ))            
      CALL JIOXXX(W(1,3   ),W(1,4   ),GAD ,ZERO    ,AWIDTH  ,W(1,45  ))    
      CALL FVIXXX(W(1,12  ),W(1,45  ),GAL ,ZERO    ,ZERO    ,W(1,46  ))    
      CALL IOVXXX(W(1,46  ),W(1,6   ),W(1,2   ),GAL ,AMP(29  ))            
      CALL JIOXXX(W(1,3   ),W(1,4   ),GZD ,ZMASS   ,ZWIDTH  ,W(1,47  ))    
      CALL FVIXXX(W(1,12  ),W(1,47  ),GZL ,ZERO    ,ZERO    ,W(1,48  ))    
      CALL IOVXXX(W(1,48  ),W(1,6   ),W(1,2   ),GAL ,AMP(30  ))            
      CALL FVIXXX(W(1,12  ),W(1,2   ),GAL ,ZERO    ,ZERO    ,W(1,49  ))    
      CALL IOVXXX(W(1,49  ),W(1,6   ),W(1,45  ),GAL ,AMP(31  ))            
      CALL IOVXXX(W(1,49  ),W(1,6   ),W(1,47  ),GZL ,AMP(32  ))            
      CALL FVIXXX(W(1,18  ),W(1,45  ),GAL ,ZERO    ,ZERO    ,W(1,50  ))    
      CALL IOVXXX(W(1,50  ),W(1,6   ),W(1,1   ),GAL ,AMP(33  ))            
      CALL FVIXXX(W(1,18  ),W(1,47  ),GZL ,ZERO    ,ZERO    ,W(1,51  ))    
      CALL IOVXXX(W(1,51  ),W(1,6   ),W(1,1   ),GAL ,AMP(34  ))            
      CALL FVIXXX(W(1,18  ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,52  ))    
      CALL IOVXXX(W(1,52  ),W(1,6   ),W(1,45  ),GAL ,AMP(35  ))            
      CALL IOVXXX(W(1,52  ),W(1,6   ),W(1,47  ),GZL ,AMP(36  ))            
      CALL FVIXXX(W(1,5   ),W(1,45  ),GAL ,ZERO    ,ZERO    ,W(1,53  ))    
      CALL FVIXXX(W(1,53  ),W(1,2   ),GAL ,ZERO    ,ZERO    ,W(1,54  ))    
      CALL IOVXXX(W(1,54  ),W(1,6   ),W(1,1   ),GAL ,AMP(37  ))            
      CALL FVIXXX(W(1,5   ),W(1,47  ),GZL ,ZERO    ,ZERO    ,W(1,55  ))    
      CALL FVIXXX(W(1,55  ),W(1,2   ),GAL ,ZERO    ,ZERO    ,W(1,56  ))    
      CALL IOVXXX(W(1,56  ),W(1,6   ),W(1,1   ),GAL ,AMP(38  ))            
      CALL FVIXXX(W(1,53  ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,57  ))    
      CALL IOVXXX(W(1,57  ),W(1,6   ),W(1,2   ),GAL ,AMP(39  ))            
      CALL FVIXXX(W(1,55  ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,58  ))    
      CALL IOVXXX(W(1,58  ),W(1,6   ),W(1,2   ),GAL ,AMP(40  ))            
      AA(MU,NU) = +AMP(   1)+AMP(   2)+AMP(   3)+AMP(   4)+AMP(   5)
     &             +AMP(   6)+AMP(   7)+AMP(   8)+AMP(   9)+AMP(  10)
     &             +AMP(  11)+AMP(  12)+AMP(  13)+AMP(  14)+AMP(  15)
     &             +AMP(  16)+AMP(  17)+AMP(  18)+AMP(  19)+AMP(  20)
     &             +AMP(  21)+AMP(  22)+AMP(  23)+AMP(  24)+AMP(  25)
     &             +AMP(  26)+AMP(  27)+AMP(  28)+AMP(  29)+AMP(  30)
     &             +AMP(  31)+AMP(  32)+AMP(  33)+AMP(  34)+AMP(  35)
     &             +AMP(  36)+AMP(  37)+AMP(  38)+AMP(  39)+AMP(  40)


      endif

	   enddo
        enddo
	
	end
	
				
c ------------------------------------------------------------
	
       
      subroutine AZtoZZ1_mg2_had(q1,q2,P,lsign,az)
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : a z -> e+ e- mu+ mu-  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=  40) 
      INTEGER    NWAVEFUNCS
      PARAMETER (NWAVEFUNCS=  58) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 Q1(0:3), Q2(0:3), P(0:3,3:6)
      complex*16 az(0:3,0:3)
      integer lsign(3:6)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J,mu,nu
      COMPLEX*16 AMP(NGRAPHS)
      COMPLEX*16 W(6,NWAVEFUNCS)
C  
c  lepton spinors and Z polarization vectors
#include "tensorz.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
c
C GLOBAL VARIABLES
C  
#include "VBFNLO/utilities/coupl.inc"
   
      logical lswitch
      parameter (lswitch=.true.)
C                                                                          
C ----------
C BEGIN CODE
C ----------

      	do i = 1,6
	 w(i,3) = lep(i)
	 w(i,4) = lem(i)
	 w(i,5) = lup(i)
	 w(i,6) = lum(i)
      	enddo

      DO MU = 0,3
         DO NU = 0,3
      
      CALL VCARTX(Q1(0),zero  ,zero, MU,-1,W(1,1   ))        
      CALL VCARTX(Q2(0),zMASS ,zwidth, NU,-1,W(1,2   )) 
 
      if (Ztouu) then


      CALL FVIXXX(W(1,3   ),W(1,2   ),GZU ,ZERO    ,ZERO    ,W(1,7   ))    
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
      CALL JIOXXX(W(1,3   ),W(1,13  ),GAU ,ZERO    ,AWIDTH  ,W(1,14  ))    
      CALL FVIXXX(W(1,5   ),W(1,14  ),GAL ,ZERO    ,ZERO    ,W(1,15  ))    
      CALL IOVXXX(W(1,15  ),W(1,6   ),W(1,2   ),GZL ,AMP(5   ))            
      CALL JIOXXX(W(1,3   ),W(1,13  ),GZU ,ZMASS   ,ZWIDTH  ,W(1,16  ))    
      CALL FVIXXX(W(1,5   ),W(1,16  ),GZL ,ZERO    ,ZERO    ,W(1,17  ))    
      CALL IOVXXX(W(1,17  ),W(1,6   ),W(1,2   ),GZL ,AMP(6   ))            
      CALL FVIXXX(W(1,5   ),W(1,2   ),GZL ,ZERO    ,ZERO    ,W(1,18  ))    
      CALL IOVXXX(W(1,18  ),W(1,6   ),W(1,14  ),GAL ,AMP(7   ))            
      CALL IOVXXX(W(1,18  ),W(1,6   ),W(1,16  ),GZL ,AMP(8   ))            
      CALL FVOXXX(W(1,13  ),W(1,2   ),GZU ,ZERO    ,ZERO    ,W(1,19  ))    
      CALL JIOXXX(W(1,3   ),W(1,19  ),GAU ,ZERO    ,AWIDTH  ,W(1,20  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,20  ),GAL ,AMP(9   ))            
      CALL JIOXXX(W(1,3   ),W(1,19  ),GZU ,ZMASS   ,ZWIDTH  ,W(1,21  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,21  ),GZL ,AMP(10  ))            
      CALL FVIXXX(W(1,7   ),W(1,1   ),GAU ,ZERO    ,ZERO    ,W(1,22  ))    
      CALL JIOXXX(W(1,22  ),W(1,4   ),GAU ,ZERO    ,AWIDTH  ,W(1,23  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,23  ),GAL ,AMP(11  ))            
      CALL JIOXXX(W(1,22  ),W(1,4   ),GZU ,ZMASS   ,ZWIDTH  ,W(1,24  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,24  ),GZL ,AMP(12  ))            
      CALL JIOXXX(W(1,7   ),W(1,13  ),GAU ,ZERO    ,AWIDTH  ,W(1,25  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,25  ),GAL ,AMP(13  ))            
      CALL JIOXXX(W(1,7   ),W(1,13  ),GZU ,ZMASS   ,ZWIDTH  ,W(1,26  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,26  ),GZL ,AMP(14  ))            
      CALL FVOXXX(W(1,4   ),W(1,2   ),GZU ,ZERO    ,ZERO    ,W(1,27  ))    
      CALL JIOXXX(W(1,3   ),W(1,27  ),GAU ,ZERO    ,AWIDTH  ,W(1,28  ))    
      CALL FVIXXX(W(1,5   ),W(1,28  ),GAL ,ZERO    ,ZERO    ,W(1,29  ))    
      CALL IOVXXX(W(1,29  ),W(1,6   ),W(1,1   ),GAL ,AMP(15  ))            
      CALL JIOXXX(W(1,3   ),W(1,27  ),GZU ,ZMASS   ,ZWIDTH  ,W(1,30  ))    
      CALL FVIXXX(W(1,5   ),W(1,30  ),GZL ,ZERO    ,ZERO    ,W(1,31  ))    
      CALL IOVXXX(W(1,31  ),W(1,6   ),W(1,1   ),GAL ,AMP(16  ))            
      CALL IOVXXX(W(1,12  ),W(1,6   ),W(1,28  ),GAL ,AMP(17  ))            
      CALL IOVXXX(W(1,12  ),W(1,6   ),W(1,30  ),GZL ,AMP(18  ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GAU ,ZERO    ,ZERO    ,W(1,32  ))    
      CALL JIOXXX(W(1,32  ),W(1,4   ),GAU ,ZERO    ,AWIDTH  ,W(1,33  ))    
      CALL FVIXXX(W(1,5   ),W(1,33  ),GAL ,ZERO    ,ZERO    ,W(1,34  ))    
      CALL IOVXXX(W(1,34  ),W(1,6   ),W(1,2   ),GZL ,AMP(19  ))            
      CALL JIOXXX(W(1,32  ),W(1,4   ),GZU ,ZMASS   ,ZWIDTH  ,W(1,35  ))    
      CALL FVIXXX(W(1,5   ),W(1,35  ),GZL ,ZERO    ,ZERO    ,W(1,36  ))    
      CALL IOVXXX(W(1,36  ),W(1,6   ),W(1,2   ),GZL ,AMP(20  ))            
      CALL IOVXXX(W(1,18  ),W(1,6   ),W(1,33  ),GAL ,AMP(21  ))            
      CALL IOVXXX(W(1,18  ),W(1,6   ),W(1,35  ),GZL ,AMP(22  ))            
      CALL FVOXXX(W(1,27  ),W(1,1   ),GAU ,ZERO    ,ZERO    ,W(1,37  ))    
      CALL JIOXXX(W(1,3   ),W(1,37  ),GAU ,ZERO    ,AWIDTH  ,W(1,38  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,38  ),GAL ,AMP(23  ))            
      CALL JIOXXX(W(1,3   ),W(1,37  ),GZU ,ZMASS   ,ZWIDTH  ,W(1,39  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,39  ),GZL ,AMP(24  ))            
      CALL FVIXXX(W(1,32  ),W(1,2   ),GZU ,ZERO    ,ZERO    ,W(1,40  ))    
      CALL JIOXXX(W(1,40  ),W(1,4   ),GAU ,ZERO    ,AWIDTH  ,W(1,41  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,41  ),GAL ,AMP(25  ))            
      CALL JIOXXX(W(1,40  ),W(1,4   ),GZU ,ZMASS   ,ZWIDTH  ,W(1,42  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,42  ),GZL ,AMP(26  ))            
      CALL JIOXXX(W(1,32  ),W(1,27  ),GAU ,ZERO    ,AWIDTH  ,W(1,43  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,43  ),GAL ,AMP(27  ))            
      CALL JIOXXX(W(1,32  ),W(1,27  ),GZU ,ZMASS   ,ZWIDTH  ,W(1,44  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,44  ),GZL ,AMP(28  ))            
      CALL JIOXXX(W(1,3   ),W(1,4   ),GAU ,ZERO    ,AWIDTH  ,W(1,45  ))    
      CALL FVIXXX(W(1,12  ),W(1,45  ),GAL ,ZERO    ,ZERO    ,W(1,46  ))    
      CALL IOVXXX(W(1,46  ),W(1,6   ),W(1,2   ),GZL ,AMP(29  ))            
      CALL JIOXXX(W(1,3   ),W(1,4   ),GZU ,ZMASS   ,ZWIDTH  ,W(1,47  ))    
      CALL FVIXXX(W(1,12  ),W(1,47  ),GZL ,ZERO    ,ZERO    ,W(1,48  ))    
      CALL IOVXXX(W(1,48  ),W(1,6   ),W(1,2   ),GZL ,AMP(30  ))            
      CALL FVIXXX(W(1,12  ),W(1,2   ),GZL ,ZERO    ,ZERO    ,W(1,49  ))    
      CALL IOVXXX(W(1,49  ),W(1,6   ),W(1,45  ),GAL ,AMP(31  ))            
      CALL IOVXXX(W(1,49  ),W(1,6   ),W(1,47  ),GZL ,AMP(32  ))            
      CALL FVIXXX(W(1,18  ),W(1,45  ),GAL ,ZERO    ,ZERO    ,W(1,50  ))    
      CALL IOVXXX(W(1,50  ),W(1,6   ),W(1,1   ),GAL ,AMP(33  ))            
      CALL FVIXXX(W(1,18  ),W(1,47  ),GZL ,ZERO    ,ZERO    ,W(1,51  ))    
      CALL IOVXXX(W(1,51  ),W(1,6   ),W(1,1   ),GAL ,AMP(34  ))            
      CALL FVIXXX(W(1,18  ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,52  ))    
      CALL IOVXXX(W(1,52  ),W(1,6   ),W(1,45  ),GAL ,AMP(35  ))            
      CALL IOVXXX(W(1,52  ),W(1,6   ),W(1,47  ),GZL ,AMP(36  ))            
      CALL FVIXXX(W(1,5   ),W(1,45  ),GAL ,ZERO    ,ZERO    ,W(1,53  ))    
      CALL FVIXXX(W(1,53  ),W(1,2   ),GZL ,ZERO    ,ZERO    ,W(1,54  ))    
      CALL IOVXXX(W(1,54  ),W(1,6   ),W(1,1   ),GAL ,AMP(37  ))            
      CALL FVIXXX(W(1,5   ),W(1,47  ),GZL ,ZERO    ,ZERO    ,W(1,55  ))    
      CALL FVIXXX(W(1,55  ),W(1,2   ),GZL ,ZERO    ,ZERO    ,W(1,56  ))    
      CALL IOVXXX(W(1,56  ),W(1,6   ),W(1,1   ),GAL ,AMP(38  ))            
      CALL FVIXXX(W(1,53  ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,57  ))    
      CALL IOVXXX(W(1,57  ),W(1,6   ),W(1,2   ),GZL ,AMP(39  ))            
      CALL FVIXXX(W(1,55  ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,58  ))    
      CALL IOVXXX(W(1,58  ),W(1,6   ),W(1,2   ),GZL ,AMP(40  ))            
      AZ(MU,NU) = +AMP(   1)+AMP(   2)+AMP(   3)+AMP(   4)+AMP(   5)
     &             +AMP(   6)+AMP(   7)+AMP(   8)+AMP(   9)+AMP(  10)
     &             +AMP(  11)+AMP(  12)+AMP(  13)+AMP(  14)+AMP(  15)
     &             +AMP(  16)+AMP(  17)+AMP(  18)+AMP(  19)+AMP(  20)
     &             +AMP(  21)+AMP(  22)+AMP(  23)+AMP(  24)+AMP(  25)
     &             +AMP(  26)+AMP(  27)+AMP(  28)+AMP(  29)+AMP(  30)
     &             +AMP(  31)+AMP(  32)+AMP(  33)+AMP(  34)+AMP(  35)
     &             +AMP(  36)+AMP(  37)+AMP(  38)+AMP(  39)+AMP(  40)


      elseif (Ztodd) then


      CALL FVIXXX(W(1,3   ),W(1,2   ),GZD ,ZERO    ,ZERO    ,W(1,7   ))    
      CALL JIOXXX(W(1,7   ),W(1,4   ),GAD ,ZERO    ,AWIDTH  ,W(1,8   ))    
      CALL FVIXXX(W(1,5   ),W(1,8   ),GAL ,ZERO    ,ZERO    ,W(1,9   ))    
      CALL IOVXXX(W(1,9   ),W(1,6   ),W(1,1   ),GAL ,AMP(1   ))            
      CALL JIOXXX(W(1,7   ),W(1,4   ),GZD ,ZMASS   ,ZWIDTH  ,W(1,10  ))    
      CALL FVIXXX(W(1,5   ),W(1,10  ),GZL ,ZERO    ,ZERO    ,W(1,11  ))    
      CALL IOVXXX(W(1,11  ),W(1,6   ),W(1,1   ),GAL ,AMP(2   ))            
      CALL FVIXXX(W(1,5   ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,12  ))    
      CALL IOVXXX(W(1,12  ),W(1,6   ),W(1,8   ),GAL ,AMP(3   ))            
      CALL IOVXXX(W(1,12  ),W(1,6   ),W(1,10  ),GZL ,AMP(4   ))            
      CALL FVOXXX(W(1,4   ),W(1,1   ),GAD ,ZERO    ,ZERO    ,W(1,13  ))    
      CALL JIOXXX(W(1,3   ),W(1,13  ),GAD ,ZERO    ,AWIDTH  ,W(1,14  ))    
      CALL FVIXXX(W(1,5   ),W(1,14  ),GAL ,ZERO    ,ZERO    ,W(1,15  ))    
      CALL IOVXXX(W(1,15  ),W(1,6   ),W(1,2   ),GZL ,AMP(5   ))            
      CALL JIOXXX(W(1,3   ),W(1,13  ),GZD ,ZMASS   ,ZWIDTH  ,W(1,16  ))    
      CALL FVIXXX(W(1,5   ),W(1,16  ),GZL ,ZERO    ,ZERO    ,W(1,17  ))    
      CALL IOVXXX(W(1,17  ),W(1,6   ),W(1,2   ),GZL ,AMP(6   ))            
      CALL FVIXXX(W(1,5   ),W(1,2   ),GZL ,ZERO    ,ZERO    ,W(1,18  ))    
      CALL IOVXXX(W(1,18  ),W(1,6   ),W(1,14  ),GAL ,AMP(7   ))            
      CALL IOVXXX(W(1,18  ),W(1,6   ),W(1,16  ),GZL ,AMP(8   ))            
      CALL FVOXXX(W(1,13  ),W(1,2   ),GZD ,ZERO    ,ZERO    ,W(1,19  ))    
      CALL JIOXXX(W(1,3   ),W(1,19  ),GAD ,ZERO    ,AWIDTH  ,W(1,20  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,20  ),GAL ,AMP(9   ))            
      CALL JIOXXX(W(1,3   ),W(1,19  ),GZD ,ZMASS   ,ZWIDTH  ,W(1,21  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,21  ),GZL ,AMP(10  ))            
      CALL FVIXXX(W(1,7   ),W(1,1   ),GAD ,ZERO    ,ZERO    ,W(1,22  ))    
      CALL JIOXXX(W(1,22  ),W(1,4   ),GAD ,ZERO    ,AWIDTH  ,W(1,23  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,23  ),GAL ,AMP(11  ))            
      CALL JIOXXX(W(1,22  ),W(1,4   ),GZD ,ZMASS   ,ZWIDTH  ,W(1,24  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,24  ),GZL ,AMP(12  ))            
      CALL JIOXXX(W(1,7   ),W(1,13  ),GAD ,ZERO    ,AWIDTH  ,W(1,25  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,25  ),GAL ,AMP(13  ))            
      CALL JIOXXX(W(1,7   ),W(1,13  ),GZD ,ZMASS   ,ZWIDTH  ,W(1,26  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,26  ),GZL ,AMP(14  ))            
      CALL FVOXXX(W(1,4   ),W(1,2   ),GZD ,ZERO    ,ZERO    ,W(1,27  ))    
      CALL JIOXXX(W(1,3   ),W(1,27  ),GAD ,ZERO    ,AWIDTH  ,W(1,28  ))    
      CALL FVIXXX(W(1,5   ),W(1,28  ),GAL ,ZERO    ,ZERO    ,W(1,29  ))    
      CALL IOVXXX(W(1,29  ),W(1,6   ),W(1,1   ),GAL ,AMP(15  ))            
      CALL JIOXXX(W(1,3   ),W(1,27  ),GZD ,ZMASS   ,ZWIDTH  ,W(1,30  ))    
      CALL FVIXXX(W(1,5   ),W(1,30  ),GZL ,ZERO    ,ZERO    ,W(1,31  ))    
      CALL IOVXXX(W(1,31  ),W(1,6   ),W(1,1   ),GAL ,AMP(16  ))            
      CALL IOVXXX(W(1,12  ),W(1,6   ),W(1,28  ),GAL ,AMP(17  ))            
      CALL IOVXXX(W(1,12  ),W(1,6   ),W(1,30  ),GZL ,AMP(18  ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GAD ,ZERO    ,ZERO    ,W(1,32  ))    
      CALL JIOXXX(W(1,32  ),W(1,4   ),GAD ,ZERO    ,AWIDTH  ,W(1,33  ))    
      CALL FVIXXX(W(1,5   ),W(1,33  ),GAL ,ZERO    ,ZERO    ,W(1,34  ))    
      CALL IOVXXX(W(1,34  ),W(1,6   ),W(1,2   ),GZL ,AMP(19  ))            
      CALL JIOXXX(W(1,32  ),W(1,4   ),GZD ,ZMASS   ,ZWIDTH  ,W(1,35  ))    
      CALL FVIXXX(W(1,5   ),W(1,35  ),GZL ,ZERO    ,ZERO    ,W(1,36  ))    
      CALL IOVXXX(W(1,36  ),W(1,6   ),W(1,2   ),GZL ,AMP(20  ))            
      CALL IOVXXX(W(1,18  ),W(1,6   ),W(1,33  ),GAL ,AMP(21  ))            
      CALL IOVXXX(W(1,18  ),W(1,6   ),W(1,35  ),GZL ,AMP(22  ))            
      CALL FVOXXX(W(1,27  ),W(1,1   ),GAD ,ZERO    ,ZERO    ,W(1,37  ))    
      CALL JIOXXX(W(1,3   ),W(1,37  ),GAD ,ZERO    ,AWIDTH  ,W(1,38  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,38  ),GAL ,AMP(23  ))            
      CALL JIOXXX(W(1,3   ),W(1,37  ),GZD ,ZMASS   ,ZWIDTH  ,W(1,39  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,39  ),GZL ,AMP(24  ))            
      CALL FVIXXX(W(1,32  ),W(1,2   ),GZD ,ZERO    ,ZERO    ,W(1,40  ))    
      CALL JIOXXX(W(1,40  ),W(1,4   ),GAD ,ZERO    ,AWIDTH  ,W(1,41  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,41  ),GAL ,AMP(25  ))            
      CALL JIOXXX(W(1,40  ),W(1,4   ),GZD ,ZMASS   ,ZWIDTH  ,W(1,42  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,42  ),GZL ,AMP(26  ))            
      CALL JIOXXX(W(1,32  ),W(1,27  ),GAD ,ZERO    ,AWIDTH  ,W(1,43  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,43  ),GAL ,AMP(27  ))            
      CALL JIOXXX(W(1,32  ),W(1,27  ),GZD ,ZMASS   ,ZWIDTH  ,W(1,44  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,44  ),GZL ,AMP(28  ))            
      CALL JIOXXX(W(1,3   ),W(1,4   ),GAD ,ZERO    ,AWIDTH  ,W(1,45  ))    
      CALL FVIXXX(W(1,12  ),W(1,45  ),GAL ,ZERO    ,ZERO    ,W(1,46  ))    
      CALL IOVXXX(W(1,46  ),W(1,6   ),W(1,2   ),GZL ,AMP(29  ))            
      CALL JIOXXX(W(1,3   ),W(1,4   ),GZD ,ZMASS   ,ZWIDTH  ,W(1,47  ))    
      CALL FVIXXX(W(1,12  ),W(1,47  ),GZL ,ZERO    ,ZERO    ,W(1,48  ))    
      CALL IOVXXX(W(1,48  ),W(1,6   ),W(1,2   ),GZL ,AMP(30  ))            
      CALL FVIXXX(W(1,12  ),W(1,2   ),GZL ,ZERO    ,ZERO    ,W(1,49  ))    
      CALL IOVXXX(W(1,49  ),W(1,6   ),W(1,45  ),GAL ,AMP(31  ))            
      CALL IOVXXX(W(1,49  ),W(1,6   ),W(1,47  ),GZL ,AMP(32  ))            
      CALL FVIXXX(W(1,18  ),W(1,45  ),GAL ,ZERO    ,ZERO    ,W(1,50  ))    
      CALL IOVXXX(W(1,50  ),W(1,6   ),W(1,1   ),GAL ,AMP(33  ))            
      CALL FVIXXX(W(1,18  ),W(1,47  ),GZL ,ZERO    ,ZERO    ,W(1,51  ))    
      CALL IOVXXX(W(1,51  ),W(1,6   ),W(1,1   ),GAL ,AMP(34  ))            
      CALL FVIXXX(W(1,18  ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,52  ))    
      CALL IOVXXX(W(1,52  ),W(1,6   ),W(1,45  ),GAL ,AMP(35  ))            
      CALL IOVXXX(W(1,52  ),W(1,6   ),W(1,47  ),GZL ,AMP(36  ))            
      CALL FVIXXX(W(1,5   ),W(1,45  ),GAL ,ZERO    ,ZERO    ,W(1,53  ))    
      CALL FVIXXX(W(1,53  ),W(1,2   ),GZL ,ZERO    ,ZERO    ,W(1,54  ))    
      CALL IOVXXX(W(1,54  ),W(1,6   ),W(1,1   ),GAL ,AMP(37  ))            
      CALL FVIXXX(W(1,5   ),W(1,47  ),GZL ,ZERO    ,ZERO    ,W(1,55  ))    
      CALL FVIXXX(W(1,55  ),W(1,2   ),GZL ,ZERO    ,ZERO    ,W(1,56  ))    
      CALL IOVXXX(W(1,56  ),W(1,6   ),W(1,1   ),GAL ,AMP(38  ))            
      CALL FVIXXX(W(1,53  ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,57  ))    
      CALL IOVXXX(W(1,57  ),W(1,6   ),W(1,2   ),GZL ,AMP(39  ))            
      CALL FVIXXX(W(1,55  ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,58  ))    
      CALL IOVXXX(W(1,58  ),W(1,6   ),W(1,2   ),GZL ,AMP(40  ))            
      AZ(MU,NU) = +AMP(   1)+AMP(   2)+AMP(   3)+AMP(   4)+AMP(   5)
     &             +AMP(   6)+AMP(   7)+AMP(   8)+AMP(   9)+AMP(  10)
     &             +AMP(  11)+AMP(  12)+AMP(  13)+AMP(  14)+AMP(  15)
     &             +AMP(  16)+AMP(  17)+AMP(  18)+AMP(  19)+AMP(  20)
     &             +AMP(  21)+AMP(  22)+AMP(  23)+AMP(  24)+AMP(  25)
     &             +AMP(  26)+AMP(  27)+AMP(  28)+AMP(  29)+AMP(  30)
     &             +AMP(  31)+AMP(  32)+AMP(  33)+AMP(  34)+AMP(  35)
     &             +AMP(  36)+AMP(  37)+AMP(  38)+AMP(  39)+AMP(  40)


      endif
	
	  enddo
       enddo	

      END

c ------------------------------------------------------------	
c       
      subroutine ZZtoZZ1_mg2_had(q1,q2,P,lsign,zz)
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : z z -> e+ e- mu+ mu-  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS 
      PARAMETER (NGRAPHS=  43) 
      INTEGER    NWAVEFUNCS   
      PARAMETER (NWAVEFUNCS=  64) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 Q1(0:3), Q2(0:3), P(0:3,3:6)
      complex*16 zz(0:3,0:3)
      integer lsign(3:6)
C  
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J,mu,nu
      COMPLEX*16 AMP(NGRAPHS)
      COMPLEX*16 W(6,NWAVEFUNCS)
C  
c  lepton spinors and Z polarization vectors
#include "tensorz.inc"
C  
C  
C GLOBAL VARIABLES
C  
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
C   
      logical lswitch
      parameter (lswitch=.true.)
C                                                                          
C ----------
C BEGIN CODE
C ----------

      	do i = 1,6
	 w(i,3) = lep(i)
	 w(i,4) = lem(i)
	 w(i,5) = lup(i)
	 w(i,6) = lum(i)
      	enddo

      DO MU = 0,3
         DO NU = 0,3
      
      CALL VCARTX(Q1(0),ZMASS ,zwidth, MU,-1,W(1,1   ))        
      CALL VCARTX(Q2(0),ZMASS ,zwidth, NU,-1,W(1,2   )) 
 
      if (Ztouu) then


      CALL FVIXXX(W(1,3   ),W(1,2   ),GZU ,ZERO    ,ZERO    ,W(1,7   ))    
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
      CALL JIOXXX(W(1,3   ),W(1,13  ),GAU ,ZERO    ,AWIDTH  ,W(1,14  ))    
      CALL FVIXXX(W(1,5   ),W(1,14  ),GAL ,ZERO    ,ZERO    ,W(1,15  ))    
      CALL IOVXXX(W(1,15  ),W(1,6   ),W(1,2   ),GZL ,AMP(5   ))            
      CALL JIOXXX(W(1,3   ),W(1,13  ),GZU ,ZMASS   ,ZWIDTH  ,W(1,16  ))    
      CALL FVIXXX(W(1,5   ),W(1,16  ),GZL ,ZERO    ,ZERO    ,W(1,17  ))    
      CALL IOVXXX(W(1,17  ),W(1,6   ),W(1,2   ),GZL ,AMP(6   ))            
      CALL FVIXXX(W(1,5   ),W(1,2   ),GZL ,ZERO    ,ZERO    ,W(1,18  ))    
      CALL IOVXXX(W(1,18  ),W(1,6   ),W(1,14  ),GAL ,AMP(7   ))            
      CALL IOVXXX(W(1,18  ),W(1,6   ),W(1,16  ),GZL ,AMP(8   ))            
      CALL FVOXXX(W(1,13  ),W(1,2   ),GZU ,ZERO    ,ZERO    ,W(1,19  ))    
      CALL JIOXXX(W(1,3   ),W(1,19  ),GAU ,ZERO    ,AWIDTH  ,W(1,20  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,20  ),GAL ,AMP(9   ))            
      CALL JIOXXX(W(1,3   ),W(1,19  ),GZU ,ZMASS   ,ZWIDTH  ,W(1,21  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,21  ),GZL ,AMP(10  ))            
      CALL FVIXXX(W(1,7   ),W(1,1   ),GZU ,ZERO    ,ZERO    ,W(1,22  ))    
      CALL JIOXXX(W(1,22  ),W(1,4   ),GAU ,ZERO    ,AWIDTH  ,W(1,23  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,23  ),GAL ,AMP(11  ))            
      CALL JIOXXX(W(1,22  ),W(1,4   ),GZU ,ZMASS   ,ZWIDTH  ,W(1,24  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,24  ),GZL ,AMP(12  ))            
      CALL JIOXXX(W(1,7   ),W(1,13  ),GAU ,ZERO    ,AWIDTH  ,W(1,25  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,25  ),GAL ,AMP(13  ))            
      CALL JIOXXX(W(1,7   ),W(1,13  ),GZU ,ZMASS   ,ZWIDTH  ,W(1,26  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,26  ),GZL ,AMP(14  ))            
      CALL FVOXXX(W(1,4   ),W(1,2   ),GZU ,ZERO    ,ZERO    ,W(1,27  ))    
      CALL JIOXXX(W(1,3   ),W(1,27  ),GAU ,ZERO    ,AWIDTH  ,W(1,28  ))    
      CALL FVIXXX(W(1,5   ),W(1,28  ),GAL ,ZERO    ,ZERO    ,W(1,29  ))    
      CALL IOVXXX(W(1,29  ),W(1,6   ),W(1,1   ),GZL ,AMP(15  ))            
      CALL JIOXXX(W(1,3   ),W(1,27  ),GZU ,ZMASS   ,ZWIDTH  ,W(1,30  ))    
      CALL FVIXXX(W(1,5   ),W(1,30  ),GZL ,ZERO    ,ZERO    ,W(1,31  ))    
      CALL IOVXXX(W(1,31  ),W(1,6   ),W(1,1   ),GZL ,AMP(16  ))            
      CALL IOVXXX(W(1,12  ),W(1,6   ),W(1,28  ),GAL ,AMP(17  ))            
      CALL IOVXXX(W(1,12  ),W(1,6   ),W(1,30  ),GZL ,AMP(18  ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GZU ,ZERO    ,ZERO    ,W(1,32  ))    
      CALL JIOXXX(W(1,32  ),W(1,4   ),GAU ,ZERO    ,AWIDTH  ,W(1,33  ))    
      CALL FVIXXX(W(1,5   ),W(1,33  ),GAL ,ZERO    ,ZERO    ,W(1,34  ))    
      CALL IOVXXX(W(1,34  ),W(1,6   ),W(1,2   ),GZL ,AMP(19  ))            
      CALL JIOXXX(W(1,32  ),W(1,4   ),GZU ,ZMASS   ,ZWIDTH  ,W(1,35  ))    
      CALL FVIXXX(W(1,5   ),W(1,35  ),GZL ,ZERO    ,ZERO    ,W(1,36  ))    
      CALL IOVXXX(W(1,36  ),W(1,6   ),W(1,2   ),GZL ,AMP(20  ))            
      CALL IOVXXX(W(1,18  ),W(1,6   ),W(1,33  ),GAL ,AMP(21  ))            
      CALL IOVXXX(W(1,18  ),W(1,6   ),W(1,35  ),GZL ,AMP(22  ))            
      CALL FVOXXX(W(1,27  ),W(1,1   ),GZU ,ZERO    ,ZERO    ,W(1,37  ))    
      CALL JIOXXX(W(1,3   ),W(1,37  ),GAU ,ZERO    ,AWIDTH  ,W(1,38  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,38  ),GAL ,AMP(23  ))            
      CALL JIOXXX(W(1,3   ),W(1,37  ),GZU ,ZMASS   ,ZWIDTH  ,W(1,39  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,39  ),GZL ,AMP(24  ))            
      CALL FVIXXX(W(1,32  ),W(1,2   ),GZU ,ZERO    ,ZERO    ,W(1,40  ))    
      CALL JIOXXX(W(1,40  ),W(1,4   ),GAU ,ZERO    ,AWIDTH  ,W(1,41  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,41  ),GAL ,AMP(25  ))            
      CALL JIOXXX(W(1,40  ),W(1,4   ),GZU ,ZMASS   ,ZWIDTH  ,W(1,42  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,42  ),GZL ,AMP(26  ))            
      CALL JIOXXX(W(1,32  ),W(1,27  ),GAU ,ZERO    ,AWIDTH  ,W(1,43  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,43  ),GAL ,AMP(27  ))            
      CALL JIOXXX(W(1,32  ),W(1,27  ),GZU ,ZMASS   ,ZWIDTH  ,W(1,44  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,44  ),GZL ,AMP(28  ))            
      CALL JIOXXX(W(1,3   ),W(1,4   ),GAU ,ZERO    ,AWIDTH  ,W(1,45  ))    
      CALL FVIXXX(W(1,12  ),W(1,45  ),GAL ,ZERO    ,ZERO    ,W(1,46  ))    
      CALL IOVXXX(W(1,46  ),W(1,6   ),W(1,2   ),GZL ,AMP(29  ))            
      CALL JIOXXX(W(1,3   ),W(1,4   ),GZU ,ZMASS   ,ZWIDTH  ,W(1,47  ))    
      CALL FVIXXX(W(1,12  ),W(1,47  ),GZL ,ZERO    ,ZERO    ,W(1,48  ))    
      CALL IOVXXX(W(1,48  ),W(1,6   ),W(1,2   ),GZL ,AMP(30  ))            
      CALL HVVXXX(W(1,47  ),W(1,2   ),GZZH ,HMASS   ,HWIDTH  ,W(1,         
     &     49  ))                                                          
      CALL JVSXXX(W(1,1   ),W(1,49  ),GZZH ,ZMASS   ,ZWIDTH  ,W(1,         
     &     50  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,50  ),GZL ,AMP(31  ))            
      CALL FVIXXX(W(1,12  ),W(1,2   ),GZL ,ZERO    ,ZERO    ,W(1,51  ))    
      CALL IOVXXX(W(1,51  ),W(1,6   ),W(1,45  ),GAL ,AMP(32  ))            
      CALL IOVXXX(W(1,51  ),W(1,6   ),W(1,47  ),GZL ,AMP(33  ))            
      CALL FVIXXX(W(1,18  ),W(1,45  ),GAL ,ZERO    ,ZERO    ,W(1,52  ))    
      CALL IOVXXX(W(1,52  ),W(1,6   ),W(1,1   ),GZL ,AMP(34  ))            
      CALL FVIXXX(W(1,18  ),W(1,47  ),GZL ,ZERO    ,ZERO    ,W(1,53  ))    
      CALL IOVXXX(W(1,53  ),W(1,6   ),W(1,1   ),GZL ,AMP(35  ))            
      CALL HVVXXX(W(1,1   ),W(1,47  ),GZZH ,HMASS   ,HWIDTH  ,W(1,         
     &     54  ))                                                          
      CALL JVSXXX(W(1,2   ),W(1,54  ),GZZH ,ZMASS   ,ZWIDTH  ,W(1,         
     &     55  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,55  ),GZL ,AMP(36  ))            
      CALL FVIXXX(W(1,18  ),W(1,1   ),GZL ,ZERO    ,ZERO    ,W(1,56  ))    
      CALL IOVXXX(W(1,56  ),W(1,6   ),W(1,45  ),GAL ,AMP(37  ))            
      CALL IOVXXX(W(1,56  ),W(1,6   ),W(1,47  ),GZL ,AMP(38  ))            
      CALL FVIXXX(W(1,5   ),W(1,45  ),GAL ,ZERO    ,ZERO    ,W(1,57  ))    
      CALL FVIXXX(W(1,57  ),W(1,2   ),GZL ,ZERO    ,ZERO    ,W(1,58  ))    
      CALL IOVXXX(W(1,58  ),W(1,6   ),W(1,1   ),GZL ,AMP(39  ))            
      CALL FVIXXX(W(1,5   ),W(1,47  ),GZL ,ZERO    ,ZERO    ,W(1,59  ))    
      CALL FVIXXX(W(1,59  ),W(1,2   ),GZL ,ZERO    ,ZERO    ,W(1,60  ))    
      CALL IOVXXX(W(1,60  ),W(1,6   ),W(1,1   ),GZL ,AMP(40  ))            
      CALL FVIXXX(W(1,57  ),W(1,1   ),GZL ,ZERO    ,ZERO    ,W(1,61  ))    
      CALL IOVXXX(W(1,61  ),W(1,6   ),W(1,2   ),GZL ,AMP(41  ))            
      CALL FVIXXX(W(1,59  ),W(1,1   ),GZL ,ZERO    ,ZERO    ,W(1,62  ))    
      CALL IOVXXX(W(1,62  ),W(1,6   ),W(1,2   ),GZL ,AMP(42  ))            
      CALL HVVXXX(W(1,1   ),W(1,2   ),GZZH ,HMASS   ,HWIDTH  ,W(1,         
     &     63  ))                                                          
      CALL JVSXXX(W(1,47  ),W(1,63  ),GZZH ,ZMASS   ,ZWIDTH  ,W(1,         
     &     64  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,64  ),GZL ,AMP(43  ))            
      ZZ(MU,NU) = +AMP(   1)+AMP(   2)+AMP(   3)+AMP(   4)+AMP(   5)
     &             +AMP(   6)+AMP(   7)+AMP(   8)+AMP(   9)+AMP(  10)
     &             +AMP(  11)+AMP(  12)+AMP(  13)+AMP(  14)+AMP(  15)
     &             +AMP(  16)+AMP(  17)+AMP(  18)+AMP(  19)+AMP(  20)
     &             +AMP(  21)+AMP(  22)+AMP(  23)+AMP(  24)+AMP(  25)
     &             +AMP(  26)+AMP(  27)+AMP(  28)+AMP(  29)+AMP(  30)
     &             +AMP(  31)+AMP(  32)+AMP(  33)+AMP(  34)+AMP(  35)
     &             +AMP(  36)+AMP(  37)+AMP(  38)+AMP(  39)+AMP(  40)
     &             +AMP(  41)+AMP(  42)+AMP(  43)


      elseif (Ztodd) then


      CALL FVIXXX(W(1,3   ),W(1,2   ),GZD ,ZERO    ,ZERO    ,W(1,7   ))    
      CALL JIOXXX(W(1,7   ),W(1,4   ),GAD ,ZERO    ,AWIDTH  ,W(1,8   ))    
      CALL FVIXXX(W(1,5   ),W(1,8   ),GAL ,ZERO    ,ZERO    ,W(1,9   ))    
      CALL IOVXXX(W(1,9   ),W(1,6   ),W(1,1   ),GZL ,AMP(1   ))            
      CALL JIOXXX(W(1,7   ),W(1,4   ),GZD ,ZMASS   ,ZWIDTH  ,W(1,10  ))    
      CALL FVIXXX(W(1,5   ),W(1,10  ),GZL ,ZERO    ,ZERO    ,W(1,11  ))    
      CALL IOVXXX(W(1,11  ),W(1,6   ),W(1,1   ),GZL ,AMP(2   ))            
      CALL FVIXXX(W(1,5   ),W(1,1   ),GZL ,ZERO    ,ZERO    ,W(1,12  ))    
      CALL IOVXXX(W(1,12  ),W(1,6   ),W(1,8   ),GAL ,AMP(3   ))            
      CALL IOVXXX(W(1,12  ),W(1,6   ),W(1,10  ),GZL ,AMP(4   ))            
      CALL FVOXXX(W(1,4   ),W(1,1   ),GZD ,ZERO    ,ZERO    ,W(1,13  ))    
      CALL JIOXXX(W(1,3   ),W(1,13  ),GAD ,ZERO    ,AWIDTH  ,W(1,14  ))    
      CALL FVIXXX(W(1,5   ),W(1,14  ),GAL ,ZERO    ,ZERO    ,W(1,15  ))    
      CALL IOVXXX(W(1,15  ),W(1,6   ),W(1,2   ),GZL ,AMP(5   ))            
      CALL JIOXXX(W(1,3   ),W(1,13  ),GZD ,ZMASS   ,ZWIDTH  ,W(1,16  ))    
      CALL FVIXXX(W(1,5   ),W(1,16  ),GZL ,ZERO    ,ZERO    ,W(1,17  ))    
      CALL IOVXXX(W(1,17  ),W(1,6   ),W(1,2   ),GZL ,AMP(6   ))            
      CALL FVIXXX(W(1,5   ),W(1,2   ),GZL ,ZERO    ,ZERO    ,W(1,18  ))    
      CALL IOVXXX(W(1,18  ),W(1,6   ),W(1,14  ),GAL ,AMP(7   ))            
      CALL IOVXXX(W(1,18  ),W(1,6   ),W(1,16  ),GZL ,AMP(8   ))            
      CALL FVOXXX(W(1,13  ),W(1,2   ),GZD ,ZERO    ,ZERO    ,W(1,19  ))    
      CALL JIOXXX(W(1,3   ),W(1,19  ),GAD ,ZERO    ,AWIDTH  ,W(1,20  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,20  ),GAL ,AMP(9   ))            
      CALL JIOXXX(W(1,3   ),W(1,19  ),GZD ,ZMASS   ,ZWIDTH  ,W(1,21  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,21  ),GZL ,AMP(10  ))            
      CALL FVIXXX(W(1,7   ),W(1,1   ),GZD ,ZERO    ,ZERO    ,W(1,22  ))    
      CALL JIOXXX(W(1,22  ),W(1,4   ),GAD ,ZERO    ,AWIDTH  ,W(1,23  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,23  ),GAL ,AMP(11  ))            
      CALL JIOXXX(W(1,22  ),W(1,4   ),GZD ,ZMASS   ,ZWIDTH  ,W(1,24  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,24  ),GZL ,AMP(12  ))            
      CALL JIOXXX(W(1,7   ),W(1,13  ),GAD ,ZERO    ,AWIDTH  ,W(1,25  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,25  ),GAL ,AMP(13  ))            
      CALL JIOXXX(W(1,7   ),W(1,13  ),GZD ,ZMASS   ,ZWIDTH  ,W(1,26  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,26  ),GZL ,AMP(14  ))            
      CALL FVOXXX(W(1,4   ),W(1,2   ),GZD ,ZERO    ,ZERO    ,W(1,27  ))    
      CALL JIOXXX(W(1,3   ),W(1,27  ),GAD ,ZERO    ,AWIDTH  ,W(1,28  ))    
      CALL FVIXXX(W(1,5   ),W(1,28  ),GAL ,ZERO    ,ZERO    ,W(1,29  ))    
      CALL IOVXXX(W(1,29  ),W(1,6   ),W(1,1   ),GZL ,AMP(15  ))            
      CALL JIOXXX(W(1,3   ),W(1,27  ),GZD ,ZMASS   ,ZWIDTH  ,W(1,30  ))    
      CALL FVIXXX(W(1,5   ),W(1,30  ),GZL ,ZERO    ,ZERO    ,W(1,31  ))    
      CALL IOVXXX(W(1,31  ),W(1,6   ),W(1,1   ),GZL ,AMP(16  ))            
      CALL IOVXXX(W(1,12  ),W(1,6   ),W(1,28  ),GAL ,AMP(17  ))            
      CALL IOVXXX(W(1,12  ),W(1,6   ),W(1,30  ),GZL ,AMP(18  ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GZD ,ZERO    ,ZERO    ,W(1,32  ))    
      CALL JIOXXX(W(1,32  ),W(1,4   ),GAD ,ZERO    ,AWIDTH  ,W(1,33  ))    
      CALL FVIXXX(W(1,5   ),W(1,33  ),GAL ,ZERO    ,ZERO    ,W(1,34  ))    
      CALL IOVXXX(W(1,34  ),W(1,6   ),W(1,2   ),GZL ,AMP(19  ))            
      CALL JIOXXX(W(1,32  ),W(1,4   ),GZD ,ZMASS   ,ZWIDTH  ,W(1,35  ))    
      CALL FVIXXX(W(1,5   ),W(1,35  ),GZL ,ZERO    ,ZERO    ,W(1,36  ))    
      CALL IOVXXX(W(1,36  ),W(1,6   ),W(1,2   ),GZL ,AMP(20  ))            
      CALL IOVXXX(W(1,18  ),W(1,6   ),W(1,33  ),GAL ,AMP(21  ))            
      CALL IOVXXX(W(1,18  ),W(1,6   ),W(1,35  ),GZL ,AMP(22  ))            
      CALL FVOXXX(W(1,27  ),W(1,1   ),GZD ,ZERO    ,ZERO    ,W(1,37  ))    
      CALL JIOXXX(W(1,3   ),W(1,37  ),GAD ,ZERO    ,AWIDTH  ,W(1,38  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,38  ),GAL ,AMP(23  ))            
      CALL JIOXXX(W(1,3   ),W(1,37  ),GZD ,ZMASS   ,ZWIDTH  ,W(1,39  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,39  ),GZL ,AMP(24  ))            
      CALL FVIXXX(W(1,32  ),W(1,2   ),GZD ,ZERO    ,ZERO    ,W(1,40  ))    
      CALL JIOXXX(W(1,40  ),W(1,4   ),GAD ,ZERO    ,AWIDTH  ,W(1,41  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,41  ),GAL ,AMP(25  ))            
      CALL JIOXXX(W(1,40  ),W(1,4   ),GZD ,ZMASS   ,ZWIDTH  ,W(1,42  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,42  ),GZL ,AMP(26  ))            
      CALL JIOXXX(W(1,32  ),W(1,27  ),GAD ,ZERO    ,AWIDTH  ,W(1,43  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,43  ),GAL ,AMP(27  ))            
      CALL JIOXXX(W(1,32  ),W(1,27  ),GZD ,ZMASS   ,ZWIDTH  ,W(1,44  ))    
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,44  ),GZL ,AMP(28  ))            
      CALL JIOXXX(W(1,3   ),W(1,4   ),GAD ,ZERO    ,AWIDTH  ,W(1,45  ))    
      CALL FVIXXX(W(1,12  ),W(1,45  ),GAL ,ZERO    ,ZERO    ,W(1,46  ))    
      CALL IOVXXX(W(1,46  ),W(1,6   ),W(1,2   ),GZL ,AMP(29  ))            
      CALL JIOXXX(W(1,3   ),W(1,4   ),GZD ,ZMASS   ,ZWIDTH  ,W(1,47  ))    
      CALL FVIXXX(W(1,12  ),W(1,47  ),GZL ,ZERO    ,ZERO    ,W(1,48  ))    
      CALL IOVXXX(W(1,48  ),W(1,6   ),W(1,2   ),GZL ,AMP(30  ))            
      CALL HVVXXX(W(1,47  ),W(1,2   ),GZZH ,HMASS   ,HWIDTH  ,W(1,         
     &     49  ))                                                          
      CALL JVSXXX(W(1,1   ),W(1,49  ),GZZH ,ZMASS   ,ZWIDTH  ,W(1,         
     &     50  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,50  ),GZL ,AMP(31  ))            
      CALL FVIXXX(W(1,12  ),W(1,2   ),GZL ,ZERO    ,ZERO    ,W(1,51  ))    
      CALL IOVXXX(W(1,51  ),W(1,6   ),W(1,45  ),GAL ,AMP(32  ))            
      CALL IOVXXX(W(1,51  ),W(1,6   ),W(1,47  ),GZL ,AMP(33  ))            
      CALL FVIXXX(W(1,18  ),W(1,45  ),GAL ,ZERO    ,ZERO    ,W(1,52  ))    
      CALL IOVXXX(W(1,52  ),W(1,6   ),W(1,1   ),GZL ,AMP(34  ))            
      CALL FVIXXX(W(1,18  ),W(1,47  ),GZL ,ZERO    ,ZERO    ,W(1,53  ))    
      CALL IOVXXX(W(1,53  ),W(1,6   ),W(1,1   ),GZL ,AMP(35  ))            
      CALL HVVXXX(W(1,1   ),W(1,47  ),GZZH ,HMASS   ,HWIDTH  ,W(1,         
     &     54  ))                                                          
      CALL JVSXXX(W(1,2   ),W(1,54  ),GZZH ,ZMASS   ,ZWIDTH  ,W(1,         
     &     55  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,55  ),GZL ,AMP(36  ))            
      CALL FVIXXX(W(1,18  ),W(1,1   ),GZL ,ZERO    ,ZERO    ,W(1,56  ))    
      CALL IOVXXX(W(1,56  ),W(1,6   ),W(1,45  ),GAL ,AMP(37  ))            
      CALL IOVXXX(W(1,56  ),W(1,6   ),W(1,47  ),GZL ,AMP(38  ))            
      CALL FVIXXX(W(1,5   ),W(1,45  ),GAL ,ZERO    ,ZERO    ,W(1,57  ))    
      CALL FVIXXX(W(1,57  ),W(1,2   ),GZL ,ZERO    ,ZERO    ,W(1,58  ))    
      CALL IOVXXX(W(1,58  ),W(1,6   ),W(1,1   ),GZL ,AMP(39  ))            
      CALL FVIXXX(W(1,5   ),W(1,47  ),GZL ,ZERO    ,ZERO    ,W(1,59  ))    
      CALL FVIXXX(W(1,59  ),W(1,2   ),GZL ,ZERO    ,ZERO    ,W(1,60  ))    
      CALL IOVXXX(W(1,60  ),W(1,6   ),W(1,1   ),GZL ,AMP(40  ))            
      CALL FVIXXX(W(1,57  ),W(1,1   ),GZL ,ZERO    ,ZERO    ,W(1,61  ))    
      CALL IOVXXX(W(1,61  ),W(1,6   ),W(1,2   ),GZL ,AMP(41  ))            
      CALL FVIXXX(W(1,59  ),W(1,1   ),GZL ,ZERO    ,ZERO    ,W(1,62  ))    
      CALL IOVXXX(W(1,62  ),W(1,6   ),W(1,2   ),GZL ,AMP(42  ))            
      CALL HVVXXX(W(1,1   ),W(1,2   ),GZZH ,HMASS   ,HWIDTH  ,W(1,         
     &     63  ))                                                          
      CALL JVSXXX(W(1,47  ),W(1,63  ),GZZH ,ZMASS   ,ZWIDTH  ,W(1,         
     &     64  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,6   ),W(1,64  ),GZL ,AMP(43  ))            
      ZZ(MU,NU) = +AMP(   1)+AMP(   2)+AMP(   3)+AMP(   4)+AMP(   5)
     &             +AMP(   6)+AMP(   7)+AMP(   8)+AMP(   9)+AMP(  10)
     &             +AMP(  11)+AMP(  12)+AMP(  13)+AMP(  14)+AMP(  15)
     &             +AMP(  16)+AMP(  17)+AMP(  18)+AMP(  19)+AMP(  20)
     &             +AMP(  21)+AMP(  22)+AMP(  23)+AMP(  24)+AMP(  25)
     &             +AMP(  26)+AMP(  27)+AMP(  28)+AMP(  29)+AMP(  30)
     &             +AMP(  31)+AMP(  32)+AMP(  33)+AMP(  34)+AMP(  35)
     &             +AMP(  36)+AMP(  37)+AMP(  38)+AMP(  39)+AMP(  40)
     &             +AMP(  41)+AMP(  42)+AMP(  43)


      endif
 
	  enddo
       enddo

      END
c ------------------------------------------------------------	

      subroutine vvtoll_mg2_had(idb,idl,h,Q,v,aa,az,za,zz)
	
      IMPLICIT NONE

        integer idb,idbb,idl,h
        real*8 q(0:3), v(0:3,4)
        complex*16 aa(0:3,0:3)
        complex*16 az(0:3,0:3)
        complex*16 za(0:3,0:3)
c        complex*16 zan(0:3,0:3)
        complex*16 zz(0:3,0:3)
	integer lsign(2,2),i,ll,mu,nu
      
c	data (lsign(ll,1),ll=1,2) /  1,  1/
c	data (lsign(ll,2),ll=1,2) / -1, -1/
	
	data (lsign(ll,1),ll=1,2) / -1, -1/
	data (lsign(ll,2),ll=1,2) / +1, +1/
	
	idbb = -idb+3
	
	if (idl.eq.1) then ! e+e-
	   if (h.eq.1.or.h.eq.2)  i = 1
	   if (h.eq.3.or.h.eq.4)  i = 2
	elseif (idl.eq.2) then ! mu+mu-
	   if (h.eq.1.or.h.eq.3)  i = 1
	   if (h.eq.2.or.h.eq.4)  i = 2
	endif   

	        call AAtoLL1_mg2_had(idb ,idl,Q,v,lsign(1,i),aa(0,0))
	        call AZtoLL1_mg2_had(idb ,idl,Q,v,lsign(1,i),az(0,0))
	        call ZAtoLL1_mg2_had(idb,idl,Q,v,lsign(1,i),za(0,0)) 
	        call ZZtoLL1_mg2_had(idb ,idl,Q,v,lsign(1,i),zz(0,0))
			
	end
		
c ------------------------------------------------------------	
       
      subroutine AAtoLL1_mg2_had(idb,idl,Q,v,lsign,aa)
C  Input:   Q(*)        one of the two incoming boson momenta 
c           idb = 1      Q = q1 = A(1) momentum
c           idb = 2      Q = q2 = A(2) momentum
c	    idl = 1	 l is e
c	    idl = 2	 l is mu
c  OUPTPUT: aa(mu,nu)  lepton tensor for a a --> l+l-  
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : a a -> e+ e-  
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
      integer idb,idl
      real*8 q(0:3), v(0:3,4), p(0:3,3:4)
      complex*16 aa(0:3,0:3)
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
#include "tensorz.inc"
C  
C GLOBAL VARIABLES
C  
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/coupl.inc"
   
      logical lswitch
      parameter (lswitch=.true.)
C                                                                          
C ----------
C BEGIN CODE
C ----------
      
      if (idb.eq.1) then
         do mu = 0,3
            q1(mu) = q(mu)
	    if (idl.eq.1) then
                q2(mu) = qe(mu)-q(mu)
		p(mu,3) = v(mu,1)
		p(mu,4) = v(mu,2)
	    else
	        q2(mu) = qu(mu)-q(mu)
		p(mu,3) = v(mu,3)
		p(mu,4) = v(mu,4)
	    endif
         enddo
      elseif (idb.eq.2) then
         do mu = 0,3
            q2(mu) = q(mu)
	    if (idl.eq.1) then
            	q1(mu) = qe(mu)-q(mu)
		p(mu,3) = v(mu,1)
		p(mu,4) = v(mu,2)
	    else
            	q1(mu) = qu(mu)-q(mu)
		p(mu,3) = v(mu,3)
		p(mu,4) = v(mu,4)
	    endif
         enddo
      endif
	
	  if (idl.eq.1) then
	    do i = 1,6
	    w(i,3) = lep(i)
	    w(i,4) = lem(i)
	    enddo
	  else
	    do i = 1,6
	    w(i,3) = lup(i)
	    w(i,4) = lum(i)
	    enddo
	  endif
             
      DO MU = 0,3
         DO NU = 0,3
      
      CALL VCARTX(Q1(0),zero ,zero, MU,-1,W(1,1	  ))	    
      CALL VCARTX(Q2(0),zero ,zero, NU,-1,W(1,2	  )) 
 
 
      CALL FVIXXX(W(1,3   ),W(1,2   ),GAL ,ZERO    ,ZERO    ,W(1,5   ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,1   ),GAL ,AMP(1   ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,6   ))    
      CALL IOVXXX(W(1,6   ),W(1,4   ),W(1,2   ),GAL ,AMP(2   ))            
      
      aa(mu,nu) = +AMP(   1)+AMP(   2)


      if (Ztouu .and. idl.eq.1) then


      CALL FVIXXX(W(1,3   ),W(1,2   ),GAU ,ZERO    ,ZERO    ,W(1,5   ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,1   ),GAU ,AMP(1   ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GAU ,ZERO    ,ZERO    ,W(1,6   ))    
      CALL IOVXXX(W(1,6   ),W(1,4   ),W(1,2   ),GAU ,AMP(2   ))            
      aa(mu,nu) = +AMP(   1)+AMP(   2)

      elseif (Ztodd .and. idl.eq.1) then

      CALL FVIXXX(W(1,3   ),W(1,2   ),GAD ,ZERO    ,ZERO    ,W(1,5   ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,1   ),GAD ,AMP(1   ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GAD ,ZERO    ,ZERO    ,W(1,6   ))    
      CALL IOVXXX(W(1,6   ),W(1,4   ),W(1,2   ),GAD ,AMP(2   ))            
      aa(mu,nu) = +AMP(   1)+AMP(   2)


      endif


	 enddo !nu
      enddo !mu

      end

c ------------------------------------------------------------	
       
      subroutine AZtoLL1_mg2_had(idb,idl,Q,v,lsign,az)
C  Input:   Q(*)        one of the two incoming boson momenta 
c           idb = 1      Q = q1 = A(1) momentum
c           idb = 2      Q = q2 = A(2) momentum
c	    idl = 1	 l is e
c	    idl = 2	 l is mu
c  OUPTPUT: aa(mu,nu)  lepton tensor for a z --> l+l-  
C  
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : a z -> e+ e-  
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
      integer idb,idl
      real*8 q(0:3), v(0:3,4), p(0:3,3:4)
      complex*16 az(0:3,0:3)
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
#include "tensorz.inc"
C  
C GLOBAL VARIABLES
C  
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
   
      logical lswitch
      parameter (lswitch=.true.)
C  
C ----------
C BEGIN CODE
C ----------
      
      if (idb.eq.1) then
         do mu = 0,3
            q1(mu) = q(mu)
	    if (idl.eq.1) then
                q2(mu) = qe(mu)-q(mu)
		p(mu,3) = v(mu,1)
		p(mu,4) = v(mu,2)
	    else
	        q2(mu) = qu(mu)-q(mu)
		p(mu,3) = v(mu,3)
		p(mu,4) = v(mu,4)
	    endif
         enddo
      elseif (idb.eq.2) then
         do mu = 0,3
            q2(mu) = q(mu)
	    if (idl.eq.1) then
            	q1(mu) = qe(mu)-q(mu)
		p(mu,3) = v(mu,1)
		p(mu,4) = v(mu,2)
	    else
            	q1(mu) = qu(mu)-q(mu)
		p(mu,3) = v(mu,3)
		p(mu,4) = v(mu,4)
	    endif
         enddo
      endif
	 
	  if (idl.eq.1) then
	    do i = 1,6
	    w(i,3) = lep(i)
	    w(i,4) = lem(i)
	    enddo
	  else
	    do i = 1,6
	    w(i,3) = lup(i)
	    w(i,4) = lum(i)
	    enddo
	  endif
             
      DO MU = 0,3
         DO NU = 0,3
      
      CALL VCARTX(Q1(0),zero  ,zero, MU,-1,W(1,1   ))        
      CALL VCARTX(Q2(0),zMASS ,zwidth, NU,-1,W(1,2   )) 
 
      CALL FVIXXX(W(1,3   ),W(1,2   ),GZL ,ZERO    ,ZERO    ,W(1,5   ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,1   ),GAL ,AMP(1   ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,6   ))    
      CALL IOVXXX(W(1,6   ),W(1,4   ),W(1,2   ),GZL ,AMP(2   ))            
      
      az(mu,nu) = +AMP(   1)+AMP(   2)


      if (Ztouu .and. idl.eq.1) then


      CALL FVIXXX(W(1,3   ),W(1,2   ),GZU ,ZERO    ,ZERO    ,W(1,5   ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,1   ),GAU ,AMP(1   ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GAU ,ZERO    ,ZERO    ,W(1,6   ))    
      CALL IOVXXX(W(1,6   ),W(1,4   ),W(1,2   ),GZU ,AMP(2   ))            

      
      az(mu,nu) = +AMP(   1)+AMP(   2)

      elseif (Ztodd .and. idl.eq.1) then


      CALL FVIXXX(W(1,3   ),W(1,2   ),GZD ,ZERO    ,ZERO    ,W(1,5   ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,1   ),GAD ,AMP(1   ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GAD ,ZERO    ,ZERO    ,W(1,6   ))    
      CALL IOVXXX(W(1,6   ),W(1,4   ),W(1,2   ),GZD ,AMP(2   ))            


      az(mu,nu) = +AMP(   1)+AMP(   2)

      endif


	 enddo !nu
      enddo !mu

      end

c ------------------------------------------------------------	
       
      subroutine ZAtoLL1_mg2_had(idb,idl,Q,v,lsign,za)
C  Input:   Q(*)        one of the two incoming boson momenta 
c           idb = 1      Q = q1 = A(1) momentum
c           idb = 2      Q = q2 = A(2) momentum
c	    idl = 1	 l is e
c	    idl = 2	 l is mu
c  OUPTPUT: aa(mu,nu)  lepton tensor for a z --> l+l-  
C  
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : z a -> e+ e-  
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
      integer idb,idl
      real*8 q(0:3), v(0:3,4), p(0:3,3:4)
      complex*16 za(0:3,0:3)
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
#include "tensorz.inc"
C  
C GLOBAL VARIABLES
C  
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/coupl.inc"
   
      logical lswitch
      parameter (lswitch=.true.)
C                                                                          
C ----------
C BEGIN CODE
C ----------
      
      if (idb.eq.1) then
         do mu = 0,3
            q1(mu) = q(mu)
	    if (idl.eq.1) then
                q2(mu) = qe(mu)-q(mu)
		p(mu,3) = v(mu,1)
		p(mu,4) = v(mu,2)
	    else
	        q2(mu) = qu(mu)-q(mu)
		p(mu,3) = v(mu,3)
		p(mu,4) = v(mu,4)
	    endif
         enddo
      elseif (idb.eq.2) then
         do mu = 0,3
            q2(mu) = q(mu)
	    if (idl.eq.1) then
            	q1(mu) = qe(mu)-q(mu)
		p(mu,3) = v(mu,1)
		p(mu,4) = v(mu,2)
	    else
            	q1(mu) = qu(mu)-q(mu)
		p(mu,3) = v(mu,3)
		p(mu,4) = v(mu,4)
	    endif
         enddo
      endif

	  if (idl.eq.1) then
	    do i = 1,6
	    w(i,3) = lep(i)
	    w(i,4) = lem(i)
	    enddo
	  else
	    do i = 1,6
	    w(i,3) = lup(i)
	    w(i,4) = lum(i)
	    enddo
	  endif
	     
      DO MU = 0,3
         DO NU = 0,3
      
      CALL VCARTX(Q1(0),zmass ,zwidth, MU,-1,W(1,1   ))        
      CALL VCARTX(Q2(0),zero  ,zero , NU,-1,W(1,2   )) 
 
      CALL FVIXXX(W(1,3   ),W(1,2   ),GAL ,ZERO    ,ZERO    ,W(1,5   ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,1   ),GZL ,AMP(1   ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GZL ,ZERO    ,ZERO    ,W(1,6   ))    
      CALL IOVXXX(W(1,6   ),W(1,4   ),W(1,2   ),GAL ,AMP(2   ))            
      
      za(mu,nu) = +AMP(   1)+AMP(   2)


      if (Ztouu .and. idl.eq.1) then


      CALL FVIXXX(W(1,3   ),W(1,2   ),GAU ,ZERO    ,ZERO    ,W(1,5   ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,1   ),GZU ,AMP(1   ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GZU ,ZERO    ,ZERO    ,W(1,6   ))    
      CALL IOVXXX(W(1,6   ),W(1,4   ),W(1,2   ),GAU ,AMP(2   ))            


      za(mu,nu) = +AMP(   1)+AMP(   2)

      elseif (Ztodd .and. idl.eq.1) then


      CALL FVIXXX(W(1,3   ),W(1,2   ),GAD ,ZERO    ,ZERO    ,W(1,5   ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,1   ),GZD ,AMP(1   ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GZD ,ZERO    ,ZERO    ,W(1,6   ))    
      CALL IOVXXX(W(1,6   ),W(1,4   ),W(1,2   ),GAD ,AMP(2   ))            


      za(mu,nu) = +AMP(   1)+AMP(   2)

      endif


 	 enddo !nu
      enddo !mu

      end


c ------------------------------------------------------------	
       
      subroutine ZZtoLL1_mg2_had(idb,idl,Q,v,lsign,zz)
C  Input:   Q(*)        one of the two incoming boson momenta 
c           idb = 1      Q = q1 = A(1) momentum
c           idb = 2      Q = q2 = A(2) momentum
c	    idl = 1	 l is e
c	    idl = 2	 l is mu
c  OUPTPUT: aa(mu,nu)  lepton tensor for z z --> l+l-  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : z z -> e+ e-  
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
      integer idb,idl
      real*8 q(0:3), v(0:3,4), p(0:3,3:4)
      complex*16 zz(0:3,0:3)
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
#include "tensorz.inc"
C  
C GLOBAL VARIABLES
C  
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/coupl.inc"
   
      logical lswitch
      parameter (lswitch=.true.)
C  
C ----------
C BEGIN CODE
C ----------
      
      if (idb.eq.1) then
         do mu = 0,3
            q1(mu) = q(mu)
	    if (idl.eq.1) then
                q2(mu) = qe(mu)-q(mu)
		p(mu,3) = v(mu,1)
		p(mu,4) = v(mu,2)
	    else
	        q2(mu) = qu(mu)-q(mu)
		p(mu,3) = v(mu,3)
		p(mu,4) = v(mu,4)
	    endif
         enddo
      elseif (idb.eq.2) then
         do mu = 0,3
            q2(mu) = q(mu)
	    if (idl.eq.1) then
            	q1(mu) = qe(mu)-q(mu)
		p(mu,3) = v(mu,1)
		p(mu,4) = v(mu,2)
	    else
            	q1(mu) = qu(mu)-q(mu)
		p(mu,3) = v(mu,3)
		p(mu,4) = v(mu,4)
	    endif
         enddo
      endif

	  if (idl.eq.1) then
	    do i = 1,6
	    w(i,3) = lep(i)
	    w(i,4) = lem(i)
	    enddo
	  else
	    do i = 1,6
	    w(i,3) = lup(i)
	    w(i,4) = lum(i)
	    enddo
	  endif
	  
	     
      DO MU = 0,3
         DO NU = 0,3
      
      CALL VCARTX(Q1(0),zMASS ,zwidth, MU,-1,W(1,1   ))        
      CALL VCARTX(Q2(0),zMASS ,zwidth, NU,-1,W(1,2   )) 
 
      CALL FVIXXX(W(1,3   ),W(1,2   ),GZL ,ZERO    ,ZERO    ,W(1,5   ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,1   ),GZL ,AMP(1   ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GZL ,ZERO    ,ZERO    ,W(1,6   ))    
      CALL IOVXXX(W(1,6   ),W(1,4   ),W(1,2   ),GZL ,AMP(2   ))            
      
      zz(mu,nu) = +AMP(   1)+AMP(   2)


      if (Ztouu .and. idl.eq.1) then


      CALL FVIXXX(W(1,3   ),W(1,2   ),GZU ,ZERO    ,ZERO    ,W(1,5   ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,1   ),GZU ,AMP(1   ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GZU ,ZERO    ,ZERO    ,W(1,6   ))    
      CALL IOVXXX(W(1,6   ),W(1,4   ),W(1,2   ),GZU ,AMP(2   ))            


      zz(mu,nu) = +AMP(   1)+AMP(   2)

      elseif (Ztodd .and. idl.eq.1) then

      CALL FVIXXX(W(1,3   ),W(1,2   ),GZD ,ZERO    ,ZERO    ,W(1,5   ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,1   ),GZD ,AMP(1   ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GZD ,ZERO    ,ZERO    ,W(1,6   ))    
      CALL IOVXXX(W(1,6   ),W(1,4   ),W(1,2   ),GZD ,AMP(2   ))


      zz(mu,nu) = +AMP(   1)+AMP(   2)

      endif


	 enddo !nu
      enddo !mu

      end
c ------------------------------------------------------------	

c ------------------------------------------------------------	

      subroutine vto4l_mg2_had(p,h,a4l,z4l)
	
      IMPLICIT NONE

      REAL*8 P(0:3,2:5)
      complex*16 a4l(0:5),z4l(0:5)
      integer lsign(2:5,4),h,ll,mu   

	! possible lepton helicities:
	data (lsign(ll,1),ll=2,5) / -1, -1, -1, -1/  !h=1
	data (lsign(ll,2),ll=2,5) / -1, -1,  1,  1/  !h=2
	data (lsign(ll,3),ll=2,5) /  1,  1, -1, -1/  !h=3
	data (lsign(ll,4),ll=2,5) /  1,  1,  1,  1/  !h=4

	call Zto4l_mg2_had(P,lsign(2,h),z4l)
	call Ato4l_mg2_had(P,lsign(2,h),a4l)
		
			
	end
		
c ------------------------------------------------------------	
c ------------------------------------------------------------	
   
      subroutine Zto4L_mg2_had(P,lsign,z)
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : z -> e+ e- mu+ mu-  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=   8) 
      INTEGER    NWAVEFUNCS 
      PARAMETER (NWAVEFUNCS=  16) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,2:5),q(0:3)
      INTEGER lsign(2:5)
      complex*16 z(0:5)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J,mu
      COMPLEX*16 AMP(NGRAPHS)
      COMPLEX*16 W(6,NWAVEFUNCS)
C  
c  lepton spinors and Z polarization vectors
#include "tensorz.inc"

C GLOBAL VARIABLES
C  
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
   
      logical lswitch
      parameter (lswitch=.true.)
C  
C ----------
C BEGIN CODE
C ----------

      	do i = 1,6
	 w(i,2) = lep(i)
	 w(i,3) = lem(i)
	 w(i,4) = lup(i)
	 w(i,5) = lum(i)
      	enddo
             
      do mu = 0,3
         q(mu) = p(mu,2)+p(mu,3)+p(mu,4)+p(mu,5)
      enddo
c
      do mu = 0,3
      CALL VCARTX(Q(0),zmass ,zwidth, MU,-1,W(1,1   ))

      if (Ztouu) then


      CALL JIOXXX(W(1,2   ),W(1,3   ),GAU ,ZERO    ,AWIDTH  ,W(1,6   ))    
      CALL FVIXXX(W(1,4   ),W(1,6   ),GAL ,ZERO    ,ZERO    ,W(1,7   ))    
      CALL IOVXXX(W(1,7   ),W(1,5   ),W(1,1   ),GZL ,AMP(1   ))            
      CALL JIOXXX(W(1,2   ),W(1,3   ),GZU ,ZMASS   ,ZWIDTH  ,W(1,8   ))    
      CALL FVIXXX(W(1,4   ),W(1,8   ),GZL ,ZERO    ,ZERO    ,W(1,9   ))    
      CALL IOVXXX(W(1,9   ),W(1,5   ),W(1,1   ),GZL ,AMP(2   ))            
      CALL FVIXXX(W(1,4   ),W(1,1   ),GZL ,ZERO    ,ZERO    ,W(1,10  ))    
      CALL IOVXXX(W(1,10  ),W(1,5   ),W(1,6   ),GAL ,AMP(3   ))            
      CALL IOVXXX(W(1,10  ),W(1,5   ),W(1,8   ),GZL ,AMP(4   ))            
      CALL FVOXXX(W(1,3   ),W(1,1   ),GZU ,ZERO    ,ZERO    ,W(1,11  ))    
      CALL JIOXXX(W(1,2   ),W(1,11  ),GAU ,ZERO    ,AWIDTH  ,W(1,12  ))    
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,12  ),GAL ,AMP(5   ))            
      CALL JIOXXX(W(1,2   ),W(1,11  ),GZU ,ZMASS   ,ZWIDTH  ,W(1,13  ))    
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,13  ),GZL ,AMP(6   ))            
      CALL FVIXXX(W(1,2   ),W(1,1   ),GZU ,ZERO    ,ZERO    ,W(1,14  ))    
      CALL JIOXXX(W(1,14  ),W(1,3   ),GAU ,ZERO    ,AWIDTH  ,W(1,15  ))    
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,15  ),GAL ,AMP(7   ))            
      CALL JIOXXX(W(1,14  ),W(1,3   ),GZU ,ZMASS   ,ZWIDTH  ,W(1,16  ))    
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,16  ),GZL ,AMP(8   ))            
      Z(MU) = +AMP(   1)+AMP(   2)+AMP(   3)+AMP(   4)+AMP(   5)
     &             +AMP(   6)+AMP(   7)+AMP(   8)


      elseif (Ztodd) then


      CALL JIOXXX(W(1,2   ),W(1,3   ),GAD ,ZERO    ,AWIDTH  ,W(1,6   ))    
      CALL FVIXXX(W(1,4   ),W(1,6   ),GAL ,ZERO    ,ZERO    ,W(1,7   ))    
      CALL IOVXXX(W(1,7   ),W(1,5   ),W(1,1   ),GZL ,AMP(1   ))            
      CALL JIOXXX(W(1,2   ),W(1,3   ),GZD ,ZMASS   ,ZWIDTH  ,W(1,8   ))    
      CALL FVIXXX(W(1,4   ),W(1,8   ),GZL ,ZERO    ,ZERO    ,W(1,9   ))    
      CALL IOVXXX(W(1,9   ),W(1,5   ),W(1,1   ),GZL ,AMP(2   ))            
      CALL FVIXXX(W(1,4   ),W(1,1   ),GZL ,ZERO    ,ZERO    ,W(1,10  ))    
      CALL IOVXXX(W(1,10  ),W(1,5   ),W(1,6   ),GAL ,AMP(3   ))            
      CALL IOVXXX(W(1,10  ),W(1,5   ),W(1,8   ),GZL ,AMP(4   ))            
      CALL FVOXXX(W(1,3   ),W(1,1   ),GZD ,ZERO    ,ZERO    ,W(1,11  ))    
      CALL JIOXXX(W(1,2   ),W(1,11  ),GAD ,ZERO    ,AWIDTH  ,W(1,12  ))    
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,12  ),GAL ,AMP(5   ))            
      CALL JIOXXX(W(1,2   ),W(1,11  ),GZD ,ZMASS   ,ZWIDTH  ,W(1,13  ))    
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,13  ),GZL ,AMP(6   ))            
      CALL FVIXXX(W(1,2   ),W(1,1   ),GZD ,ZERO    ,ZERO    ,W(1,14  ))    
      CALL JIOXXX(W(1,14  ),W(1,3   ),GAD ,ZERO    ,AWIDTH  ,W(1,15  ))    
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,15  ),GAL ,AMP(7   ))            
      CALL JIOXXX(W(1,14  ),W(1,3   ),GZD ,ZMASS   ,ZWIDTH  ,W(1,16  ))    
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,16  ),GZL ,AMP(8   ))            
      Z(MU) = +AMP(   1)+AMP(   2)+AMP(   3)+AMP(   4)+AMP(   5)
     &             +AMP(   6)+AMP(   7)+AMP(   8)


      endif

      enddo
             
      z(4) = dcmplx(q(0),q(3))
      z(5) = dcmplx(q(1),q(2))

      return
      END
      
c------------------------------------------------------------------------      
          
      subroutine Ato4L_mg2_had(P,lsign,a)
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : a -> e+ e- mu+ mu-  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=   8) 
      INTEGER    NWAVEFUNCS    
      PARAMETER (NWAVEFUNCS=  16) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,2:5),q(0:3)
      INTEGER lsign(2:5)
      complex*16 a(0:5)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J,mu
      COMPLEX*16 AMP(NGRAPHS)
      COMPLEX*16 W(6,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
C  
c  lepton spinors and Z polarization vectors
#include "tensorz.inc"
   
      logical lswitch
      parameter (lswitch=.true.)
C                                                                          
C ----------
C BEGIN CODE
C ----------

      	do i = 1,6
	 w(i,2) = lep(i)
	 w(i,3) = lem(i)
	 w(i,4) = lup(i)
	 w(i,5) = lum(i)
      	enddo

      do mu = 0,3
         q(mu) = p(mu,2)+p(mu,3)+p(mu,4)+p(mu,5)
      enddo

c
      do mu = 0,3
      CALL VCARTX(Q(0),zero ,zero, MU,-1,W(1,1   ))

      if (Ztouu) then


      CALL JIOXXX(W(1,2   ),W(1,3   ),GAU ,ZERO    ,AWIDTH  ,W(1,6   ))    
      CALL FVIXXX(W(1,4   ),W(1,6   ),GAL ,ZERO    ,ZERO    ,W(1,7   ))    
      CALL IOVXXX(W(1,7   ),W(1,5   ),W(1,1   ),GAL ,AMP(1   ))            
      CALL JIOXXX(W(1,2   ),W(1,3   ),GZU ,ZMASS   ,ZWIDTH  ,W(1,8   ))    
      CALL FVIXXX(W(1,4   ),W(1,8   ),GZL ,ZERO    ,ZERO    ,W(1,9   ))    
      CALL IOVXXX(W(1,9   ),W(1,5   ),W(1,1   ),GAL ,AMP(2   ))            
      CALL FVIXXX(W(1,4   ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,10  ))    
      CALL IOVXXX(W(1,10  ),W(1,5   ),W(1,6   ),GAL ,AMP(3   ))            
      CALL IOVXXX(W(1,10  ),W(1,5   ),W(1,8   ),GZL ,AMP(4   ))            
      CALL FVOXXX(W(1,3   ),W(1,1   ),GAU ,ZERO    ,ZERO    ,W(1,11  ))    
      CALL JIOXXX(W(1,2   ),W(1,11  ),GAU ,ZERO    ,AWIDTH  ,W(1,12  ))    
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,12  ),GAL ,AMP(5   ))            
      CALL JIOXXX(W(1,2   ),W(1,11  ),GZU ,ZMASS   ,ZWIDTH  ,W(1,13  ))    
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,13  ),GZL ,AMP(6   ))            
      CALL FVIXXX(W(1,2   ),W(1,1   ),GAU ,ZERO    ,ZERO    ,W(1,14  ))    
      CALL JIOXXX(W(1,14  ),W(1,3   ),GAU ,ZERO    ,AWIDTH  ,W(1,15  ))    
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,15  ),GAL ,AMP(7   ))            
      CALL JIOXXX(W(1,14  ),W(1,3   ),GZU ,ZMASS   ,ZWIDTH  ,W(1,16  ))    
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,16  ),GZL ,AMP(8   ))            
      A(MU) = +AMP(   1)+AMP(   2)+AMP(   3)+AMP(   4)+AMP(   5)
     &             +AMP(   6)+AMP(   7)+AMP(   8)


      elseif (Ztodd) then


      CALL JIOXXX(W(1,2   ),W(1,3   ),GAD ,ZERO    ,AWIDTH  ,W(1,6   ))    
      CALL FVIXXX(W(1,4   ),W(1,6   ),GAL ,ZERO    ,ZERO    ,W(1,7   ))    
      CALL IOVXXX(W(1,7   ),W(1,5   ),W(1,1   ),GAL ,AMP(1   ))            
      CALL JIOXXX(W(1,2   ),W(1,3   ),GZD ,ZMASS   ,ZWIDTH  ,W(1,8   ))    
      CALL FVIXXX(W(1,4   ),W(1,8   ),GZL ,ZERO    ,ZERO    ,W(1,9   ))    
      CALL IOVXXX(W(1,9   ),W(1,5   ),W(1,1   ),GAL ,AMP(2   ))            
      CALL FVIXXX(W(1,4   ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,10  ))    
      CALL IOVXXX(W(1,10  ),W(1,5   ),W(1,6   ),GAL ,AMP(3   ))            
      CALL IOVXXX(W(1,10  ),W(1,5   ),W(1,8   ),GZL ,AMP(4   ))            
      CALL FVOXXX(W(1,3   ),W(1,1   ),GAD ,ZERO    ,ZERO    ,W(1,11  ))    
      CALL JIOXXX(W(1,2   ),W(1,11  ),GAD ,ZERO    ,AWIDTH  ,W(1,12  ))    
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,12  ),GAL ,AMP(5   ))            
      CALL JIOXXX(W(1,2   ),W(1,11  ),GZD ,ZMASS   ,ZWIDTH  ,W(1,13  ))    
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,13  ),GZL ,AMP(6   ))            
      CALL FVIXXX(W(1,2   ),W(1,1   ),GAD ,ZERO    ,ZERO    ,W(1,14  ))    
      CALL JIOXXX(W(1,14  ),W(1,3   ),GAD ,ZERO    ,AWIDTH  ,W(1,15  ))    
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,15  ),GAL ,AMP(7   ))            
      CALL JIOXXX(W(1,14  ),W(1,3   ),GZD ,ZMASS   ,ZWIDTH  ,W(1,16  ))    
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,16  ),GZL ,AMP(8   ))            
      A(MU) = +AMP(   1)+AMP(   2)+AMP(   3)+AMP(   4)+AMP(   5)
     &             +AMP(   6)+AMP(   7)+AMP(   8)


      endif
     
      enddo   
     
      a(4) = dcmplx(q(0),q(3))
      a(5) = dcmplx(q(1),q(2))

      return
      END



