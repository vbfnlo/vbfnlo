cccc Partial wave decomposition subroutine      

      SUBROUTINE vvtovv_Wigner(qp1,qp2,procnum,jNLO,aaww_anomal,azww_anomal,zaww_anomal,zzww_anomal,wwww5_anomal,wwww6_anomal,aazz_anomal,azzz_anomal,zazz_anomal,zzzz_anomal,wwzz5_anomal,wwzz6_anomal)
c This routine calculates vv->(vv)4l
c using polarization vectors and Wigner decomposition.
c
c IN:
c     qp1,qp2           : momenta of the incoming vector bosons in the LAB frame
c     procmum           : number of process to call the right outgoing momenta
c OUT:    
c     xxyy_anomal       : xx->yy tensor, including leptonic decay of yy
c
      IMPLICIT NONE

      INTEGER NGRAPHS
      PARAMETER( NGRAPHS = 20)
      DOUBLE PRECISION rZERO, rONE, rPi
      PARAMETER(rZERO = 0.0D0, rONE = 1.0D0, rPi = 3.14159265358979323846d0)
      DOUBLE COMPLEX cIMAG, cZERO, cONE
      PARAMETER( cIMAG=(0.0D0,1.0D0), cZERO=(0.0D0,0.0D0), cONE=(1.0D0,0.0D0))
      INTEGER procnum

      DOUBLE PRECISION QP1(0:3),QP2(0:3),QP3(0:3),QP4(0:3)
      DOUBLE COMPLEX aaww_anomal(0:3,0:3),azww_anomal(0:3,0:3),zaww_anomal(0:3,0:3),zzww_anomal(0:3,0:3),wwww5_anomal(0:3,0:3),wwww6_anomal(0:3,0:3)
      DOUBLE COMPLEX aazz_anomal(0:3,0:3),azzz_anomal(0:3,0:3),zazz_anomal(0:3,0:3),zzzz_anomal(0:3,0:3),wwzz5_anomal(0:3,0:3),wwzz6_anomal(0:3,0:3)    
      INTEGER jNLO

      DOUBLE PRECISION Q1(0:3), Q2(0:3), Q3(0:3), Q4(0:3) 
      DOUBLE PRECISION Q1ON(0:3), Q2ON(0:3)
      DOUBLE PRECISION smalld(0:2,-2:2,-2:2)
      DOUBLE PRECISION sameprod(1:11)

      DOUBLE PRECISION Q1T, Q2T, Q3T, Q4T        
      DOUBLE PRECISION Sq1, Sq2, Sq3, Sq4
      DOUBLE PRECISION Ecm_VV
      DOUBLE PRECISION COST
      
      DOUBLE COMPLEX D1W, D1A, D1Z,D2W, D2A, D2Z, D3, D4
      DOUBLE PRECISION Nqw, Nqz
      DOUBLE PRECISION Npw(4),Npz(4)
      
      DOUBLE COMPLEX WPB(6), AUB(6), ZUB(6), WMB(6), AEB(6), ZEB(6), EAB(6)



      DOUBLE COMPLEX OFFepsw(6,-1:1,2), OFFepspw(6,-1:1,4)   
      DOUBLE COMPLEX OFFeps2w(6,-1:1,4), ONONepsw(6,-1:1,2)  
      DOUBLE COMPLEX OFFepsCw(6,-1:1,2), OFFepsC2w(6,-1:1,2) 
      DOUBLE COMPLEX OFFepsz(6,-1:1,2), OFFepspz(6,-1:1,4)   
      DOUBLE COMPLEX OFFeps2z(6,-1:1,4), ONONepsz(6,-1:1,2)  
      DOUBLE COMPLEX OFFepsCz(6,-1:1,2), OFFepsC2z(6,-1:1,2) 
      DOUBLE COMPLEX U3a(-1:1), U4a(-1:1), U4z(-1:1),U3z(-1:1),U3wp(-1:1),U4wm(-1:1)
      DOUBLE COMPLEX M(-1:1,-1:1,-1:1,-1:1,1:11,1:11),M1(-1:1,-1:1,-1:1,1:11,1:11),Mtest(-1:1,-1:1,-1:1,-1:1,1:11,1:11)
      DOUBLE COMPLEX MM(-1:1,0:3,1:11,1:11),M2(-1:1,-1:1,1:11,1:11) 
      DOUBLE COMPLEX Aoff(-1:1,-1:1,-1:1,-1:1,1:11,1:11,0:2)
      DOUBLE COMPLEX Aofftest(-1:1,-1:1,-1:1,-1:1,1:11,1:11,0:2)
      DOUBLE COMPLEX Aoo(-1:1,-1:1,-1:1,-1:1,1:11,1:11,0:2)
      DOUBLE COMPLEX U(-1:1,-1:1,-1:1,-1:1,1:11,1:11,0:2)
      DOUBLE COMPLEX vvt(0:3,0:3,1:11,1:11)


      INTEGER I,J,K,N,MU,NU
      INTEGER LAMBDA, LAMBDA1, LAMBDA2, LAMBDA3, LAMBDA4, BLOCK1, BLOCK2
      logical tunitarize3
      DOUBLE COMPLEX DOTCC
      EXTERNAL DOTCC
           
#include "VBFNLO/utilities/global.inc"


#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"
#include "gmunu.inc"
#include "transm.inc"

#include "VBFNLO/utilities/lha.inc"

#include "tensor.inc"
#include "tensorz.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/process.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"



      tunitarize3=.true.
      Nqw = rZero
      Nqz = rZero      
      Npw(1:4) = rZero
      Npz(1:4) = rZero
      smalld(0:2,-2:2,-2:2) = rZero
      M(-1:1,-1:1,-1:1,-1:1,1:11,1:11) = cZero
      Mtest(-1:1,-1:1,-1:1,-1:1,1:11,1:11) = cZero
      Aoff(-1:1,-1:1,-1:1,-1:1,1:11,1:11,0:2) = cZero
      Aofftest(-1:1,-1:1,-1:1,-1:1,1:11,1:11,0:2) = cZero
      Aoo(-1:1,-1:1,-1:1,-1:1,1:11,1:11,0:2) = cZero
      U(-1:1,-1:1,-1:1,-1:1,1:11,1:11,0:2) = cZero
      M1(-1:1,-1:1,-1:1,1:11,1:11) = cZero
      M2(-1:1,-1:1,1:11,1:11) = cZero
      MM(-1:1,0:3,1:11,1:11) = cZero
      vvt(0:3,0:3,1:11,1:11) = cZero
      aaww_anomal(0:3,0:3)=cZero
      azww_anomal(0:3,0:3)=cZero
      zaww_anomal(0:3,0:3)=cZero
      zzww_anomal(0:3,0:3)=cZero
      wwww5_anomal(0:3,0:3)=cZero
      wwww6_anomal(0:3,0:3)=cZero
      aazz_anomal(0:3,0:3)=cZero
      azzz_anomal(0:3,0:3)=cZero
      zazz_anomal(0:3,0:3)=cZero
      zzzz_anomal(0:3,0:3)=cZero
      wwzz5_anomal(0:3,0:3)=cZero
      wwzz6_anomal(0:3,0:3)=cZero 
      ecm_real = rZERO
      Ecm_VV = rZERO
      COST=rZERO
      sameprod(1:11)=2d0
      sameprod(1)=2d0
      sameprod(5)=2d0
      sameprod(9)=2d0

ccccccccccccccccccccccccc   
c     call momenta of the outgoing particles according to process number in vbfnlo.dat
      select case (procnum)
      case (200)
        

          qp3(0) = dble(wp(5))
          qp3(1) = dble(wp(6))
          qp3(2) = dimag(wp(6))
          qp3(3) = dimag(wp(5))
         
          qp4(0) = dble(wm(5))
          qp4(1) = dble(wm(6))
          qp4(2) = dimag(wm(6))
          qp4(3) = dimag(wm(5))

      case (210)
        
          qp3(0) = dble(ze(5))
          qp3(1) = dble(ze(6))
          qp3(2) = dimag(ze(6))
          qp3(3) = dimag(ze(5))
         
          qp4(0) = dble(zu(5))
          qp4(1) = dble(zu(6))
          qp4(2) = dimag(zu(6))
          qp4(3) = dimag(zu(5))
        
      case (211)
        
          qp3(0) = dble(ze(5))
          qp3(1) = dble(ze(6))
          qp3(2) = dimag(ze(6))
          qp3(3) = dimag(ze(5))
         
          qp4(0) = dble(zu(5))
          qp4(1) = dble(zu(6))
          qp4(2) = dimag(zu(6))
          qp4(3) = dimag(zu(5))
         
      case (290)
        
          qp3(0) = dble(ze(5))
          qp3(1) = dble(ze(6))
          qp3(2) = dimag(ze(6))
          qp3(3) = dimag(ze(5))
         
          qp4(0) = dble(ea(5))
          qp4(1) = dble(ea(6))
          qp4(2) = dimag(ea(6))
          qp4(3) = dimag(ea(5))
          
      case (291)
        
          qp3(0) = dble(ze(5))
          qp3(1) = dble(ze(6))
          qp3(2) = dimag(ze(6))
          qp3(3) = dimag(ze(5))
         
          qp4(0) = dble(ea(5))
          qp4(1) = dble(ea(6))
          qp4(2) = dimag(ea(6))
          qp4(3) = dimag(ea(5))
          
      end select      

      call BOOSTCM(Qp1(0),Qp2(0),Qp3(0),Qp4(0), TT(0,0,jNLO),TTinv(0,0,jNLO))       
      ecm_real = dsqrt((qp1(0)+qp2(0))**2-((qp1(1)+qp2(1))**2+(qp1(2)+qp2(2))**2+(qp1(3)+qp2(3))**2))
      
**    boost to CM Frame
      do i = 0,3
         q1(i) = TT(i,0,jNLO)*qp1(0)-TT(i,1,jNLO)*qp1(1)-TT(i,2,jNLO)*qp1(2)-TT(i,3,jNLO)*qp1(3)
         if( abs(q1(i)).le.abs(1.0d-10*ecm_real) ) then
            q1(i) = 0.0d0
         endif         
         q2(i) = TT(i,0,jNLO)*qp2(0)-TT(i,1,jNLO)*qp2(1)-TT(i,2,jNLO)*qp2(2)-TT(i,3,jNLO)*qp2(3)
         if( abs(q2(i)).le.abs(1.0d-10*ecm_real) ) then
            q2(i) = 0.0d0
         endif         
         q3(i) = TT(i,0,jNLO)*qp3(0)-TT(i,1,jNLO)*qp3(1)-TT(i,2,jNLO)*qp3(2)-TT(i,3,jNLO)*qp3(3)
         if( abs(q3(i)).le.abs(1.0d-10*ecm_real) ) then
            q3(i) = 0.0d0
         endif         
         q4(i) = TT(i,0,jNLO)*qp4(0)-TT(i,1,jNLO)*qp4(1)-TT(i,2,jNLO)*qp4(2)-TT(i,3,jNLO)*qp4(3)         
         if( abs(q4(i)).le.abs(1.0d-10*ecm_real) ) then
            q4(i) = 0.0d0
         endif
        
      enddo


c     boost currents according to process number
      do i = 0,3   
        select case (procnum)
        case (200)   
                  
            wpb(i+1) = TT(i,0,jNLO)*wp(1)-TT(i,1,jNLO)*wp(2)-TT(i,2,jNLO)*wp(3)-TT(i,3,jNLO)*wp(4)
            wmb(i+1) = TT(i,0,jNLO)*wm(1)-TT(i,1,jNLO)*wm(2)-TT(i,2,jNLO)*wm(3)-TT(i,3,jNLO)*wm(4)
          
        case (210)
          
            zub(i+1) = TT(i,0,jNLO)*zu(1)-TT(i,1,jNLO)*zu(2)-TT(i,2,jNLO)*zu(3)-TT(i,3,jNLO)*zu(4)             
            zeb(i+1) = TT(i,0,jNLO)*ze(1)-TT(i,1,jNLO)*ze(2)-TT(i,2,jNLO)*ze(3)-TT(i,3,jNLO)*ze(4)
            aub(i+1) = TT(i,0,jNLO)*au(1)-TT(i,1,jNLO)*au(2)-TT(i,2,jNLO)*au(3)-TT(i,3,jNLO)*au(4)
            aeb(i+1) = TT(i,0,jNLO)*ae(1)-TT(i,1,jNLO)*ae(2)-TT(i,2,jNLO)*ae(3)-TT(i,3,jNLO)*ae(4)
          
        case (211)
          
            zub(i+1) = TT(i,0,jNLO)*zu(1)-TT(i,1,jNLO)*zu(2)-TT(i,2,jNLO)*zu(3)-TT(i,3,jNLO)*zu(4)             
            zeb(i+1) = TT(i,0,jNLO)*ze(1)-TT(i,1,jNLO)*ze(2)-TT(i,2,jNLO)*ze(3)-TT(i,3,jNLO)*ze(4)
            aub(i+1) = TT(i,0,jNLO)*au(1)-TT(i,1,jNLO)*au(2)-TT(i,2,jNLO)*au(3)-TT(i,3,jNLO)*au(4)
            aeb(i+1) = TT(i,0,jNLO)*ae(1)-TT(i,1,jNLO)*ae(2)-TT(i,2,jNLO)*ae(3)-TT(i,3,jNLO)*ae(4)
          
        case (290)
          
            zeb(i+1) = TT(i,0,jNLO)*ze(1)-TT(i,1,jNLO)*ze(2)-TT(i,2,jNLO)*ze(3)-TT(i,3,jNLO)*ze(4)
            aeb(i+1) = TT(i,0,jNLO)*ae(1)-TT(i,1,jNLO)*ae(2)-TT(i,2,jNLO)*ae(3)-TT(i,3,jNLO)*ae(4)
            eab(i+1) = TT(i,0,jNLO)*ea(1)-TT(i,1,jNLO)*ea(2)-TT(i,2,jNLO)*ea(3)-TT(i,3,jNLO)*ea(4)  
          
        case (291)
          
            zeb(i+1) = TT(i,0,jNLO)*ze(1)-TT(i,1,jNLO)*ze(2)-TT(i,2,jNLO)*ze(3)-TT(i,3,jNLO)*ze(4)
            aeb(i+1) = TT(i,0,jNLO)*ae(1)-TT(i,1,jNLO)*ae(2)-TT(i,2,jNLO)*ae(3)-TT(i,3,jNLO)*ae(4)
            eab(i+1) = TT(i,0,jNLO)*ea(1)-TT(i,1,jNLO)*ea(2)-TT(i,2,jNLO)*ea(3)-TT(i,3,jNLO)*ea(4)  
          
        end select                      
      enddo
 
                

      q1T = q1(1)**2+q1(2)**2+q1(3)**2
      q2T = q2(1)**2+q2(2)**2+q2(3)**2
      q3T = q3(1)**2+q3(2)**2+q3(3)**2
      q4T = q4(1)**2+q4(2)**2+q4(3)**2


      Sq1 = q1(0)**2- q1T
      Sq2 = q2(0)**2- q2T
      Sq3 = q3(0)**2- q3T
      Sq4 = q4(0)**2- q4T

      Ecm_VV = (q1(0)+q2(0))**2-(q1(1)+q2(1))**2-(q1(2)+q2(2))**2-(q1(3)+q2(3))**2
            
      COST = Q3(3) / dsqrt(q3T)
c     normalization for longitudinal polarization vectors for the xx->yy matrix element, off-shell photon with 1/MZ
       Nqz   = rONE/ZMASS
       Nqw   = rONE/WMASS

       call voffxx(q1,WMASS,-1,Nqw,OFFepsw(1,-1,1)) 
       call voffxx(q2,WMASS,-1,Nqw,OFFepsw(1,-1,2))
       call voffxx(q1,ZMASS,-1,Nqz,OFFepsz(1,-1,1)) 
       call voffxx(q2,ZMASS,-1,Nqz,OFFepsz(1,-1,2))
      
    
c     normalization for longitudinal EXTERNAL polarization vectors , off-shell photon with MZ/Q
       Npw(1) = WMASS/(dabs(Sq1))
       Npw(2) = WMASS/(dabs(Sq2))
       Npw(3) = WMASS/(dabs(Sq3))
       Npw(4) = WMASS/(dabs(Sq4))
       Npz(1) = ZMASS/(dabs(Sq1))
       Npz(2) = ZMASS/(dabs(Sq2))
       Npz(3) = ZMASS/(dabs(Sq3))
       Npz(4) = ZMASS/(dabs(Sq4))       

       call voffxx(q1,WMASS,-1,Npw(1),OFFepspw(1,-1,1)) 
       call voffxx(q2,WMASS,-1,Npw(2),OFFepspw(1,-1,2))
       call voffxx(q3,WMASS,-1,Npw(3),OFFepspw(1,-1,3)) 
       call voffxx(q4,WMASS,-1,Npw(4),OFFepspw(1,-1,4))
       call voffxx(q1,ZMASS,-1,Npz(1),OFFepspz(1,-1,1)) 
       call voffxx(q2,ZMASS,-1,Npz(2),OFFepspz(1,-1,2))
       call voffxx(q3,ZMASS,-1,Npz(3),OFFepspz(1,-1,3)) 
       call voffxx(q4,ZMASS,-1,Npz(4),OFFepspz(1,-1,4))
c      complex conjugated pol. vectors
       call voffxx(q1,WMASS,+1,Npw(1),OFFepsCw(1,-1,1)) 
       call voffxx(q2,WMASS,+1,Npw(2),OFFepsCw(1,-1,2)) 
       call voffxx(q1,ZMASS,+1,Npz(1),OFFepsCz(1,-1,1)) 
       call voffxx(q2,ZMASS,+1,Npz(2),OFFepsCz(1,-1,2)) 

**     Boosted back to Lab Frame  
       do i = 1,4
        do lambda = -1,1
           OFFeps2w(1,lambda,i) = TTinv(0,0,jNLO)*OFFepspw(1,lambda,i)+TTinv(0,1,jNLO)*OFFepspw(2,lambda,i)
     &                         + TTinv(0,2,jNLO)*OFFepspw(3,lambda,i)+TTinv(0,3,jNLO)*OFFepspw(4,lambda,i)  
           OFFeps2z(1,lambda,i) = TTinv(0,0,jNLO)*OFFepspz(1,lambda,i)+TTinv(0,1,jNLO)*OFFepspz(2,lambda,i)
     &                         + TTinv(0,2,jNLO)*OFFepspz(3,lambda,i)+TTinv(0,3,jNLO)*OFFepspz(4,lambda,i) 
   
           do k = 1,3
             OFFeps2w(k+1,lambda,i) = -TTinv(k,0,jNLO)* OFFepspw(1,lambda,i)-TTinv(k,1,jNLO)*OFFepspw(2,lambda,i)
     &                               - TTinv(k,2,jNLO)* OFFepspw(3,lambda,i)-TTinv(k,3,jNLO)*OFFepspw(4,lambda,i)  
             OFFeps2z(k+1,lambda,i) = -TTinv(k,0,jNLO)* OFFepspz(1,lambda,i)-TTinv(k,1,jNLO)*OFFepspz(2,lambda,i)
     &                               - TTinv(k,2,jNLO)* OFFepspz(3,lambda,i)-TTinv(k,3,jNLO)*OFFepspz(4,lambda,i)     
      
           enddo
        enddo
       enddo 
       do i = 1,2
        do lambda = -1,1
           OFFepsC2w(1,lambda,i) = TTinv(0,0,jNLO)*OFFepsCw(1,lambda,i)+TTinv(0,1,jNLO)*OFFepsCw(2,lambda,i)
     &                         + TTinv(0,2,jNLO)*OFFepsCw(3,lambda,i)+TTinv(0,3,jNLO)*OFFepsCw(4,lambda,i) 
           OFFepsC2z(1,lambda,i) = TTinv(0,0,jNLO)*OFFepsCz(1,lambda,i)+TTinv(0,1,jNLO)*OFFepsCz(2,lambda,i)
     &                         + TTinv(0,2,jNLO)*OFFepsCz(3,lambda,i)+TTinv(0,3,jNLO)*OFFepsCz(4,lambda,i)  
   
           do k = 1,3
             OFFepsC2w(k+1,lambda,i) = -TTinv(k,0,jNLO)*OFFepsCw(1,lambda,i)-TTinv(k,1,jNLO)*OFFepsCw(2,lambda,i)
     &                               - TTinv(k,2,jNLO)* OFFepsCw(3,lambda,i)-TTinv(k,3,jNLO)*OFFepsCw(4,lambda,i)   
             OFFepsC2z(k+1,lambda,i) = -TTinv(k,0,jNLO)*OFFepsCz(1,lambda,i)-TTinv(k,1,jNLO)*OFFepsCz(2,lambda,i)
     &                               - TTinv(k,2,jNLO)* OFFepsCz(3,lambda,i)-TTinv(k,3,jNLO)*OFFepsCz(4,lambda,i)    
      
           enddo
        enddo
       enddo      
      
     
      
**     Main Calculation: M(mu,nu), mu,nu=0,3      

      if(.not.isnan(abs(COST)))then
      CALL WIGNERdvv(COST,smalld)   
      else
        smalld(0,0,0)=1
      endif    

      if(Tunitarize3) then  ! T matrix unitarization      

        CALL vvvv_wigner_final(OFFepsz(1,-1,1),OFFepsw(1,-1,1),q3,q4,Aoff)   ! Voff(1,2) -> Von(3,4)
        
        ! On-shell to On-shell      
        q1On(0) =q3(0)
        q1On(1) =rZero
        q1On(2) =rZero
        q1On(3) =sqrt(q3T)*sign(1.0d0,Q1(3))

        q2On(0) =q4(0)
        q2On(1) =-q1On(1)
        q2On(2) =-q1On(2)
        q2On(3) =-q1On(3)
          
        call voffxx(q1On,WMASS,-1,Nqw,ONONepsw(1,-1,1))
        call voffxx(q2On,WMASS,-1,Nqw,ONONepsw(1,-1,2))
        call voffxx(q1On,ZMASS,-1,Nqz,ONONepsz(1,-1,1))
        call voffxx(q2On,ZMASS,-1,Nqz,ONONepsz(1,-1,2))

        CALL vvvv_wigner_final(ONONepsz(1,-1,1),ONONepsw(1,-1,1),q3,q4,Aoo)    ! Von(3,4)  -> Von(3,4)   
        
        ! T-matriz Unitarization calculation
        CALL tmatrixvv(Aoff,Aoo, Ecm_VV,Sq1,Sq2,Sq3,Sq4,67,11, U)   

!         M(lambda1,lambda2,lambda3,lambda4, block1, block2), lambdai=-1,1  , block accounts for the different 2->2 processes with the same total charge                   
        do lambda1 = -1,1
         do lambda2 = -1,1
          do lambda3 = -1,1
           do lambda4 = -1,1
            do block1 = 1,11
             do block2 = 1,11   !adjust coefficient for same particle
              M(lambda1,lambda2,lambda3,lambda4,block1,block2) = sameprod(block2)*8.0d0*rPi*(U(lambda1,lambda2,lambda3,lambda4,block1,block2,0)*smalld(0,lambda1-lambda2,lambda3-lambda4) 
     &          + 3.0d0*U(lambda1,lambda2,lambda3,lambda4,block1,block2,1)*smalld(1,lambda1-lambda2,lambda3-lambda4) 
     &          + 5.0d0*U(lambda1,lambda2,lambda3,lambda4,block1,block2,2)*smalld(2,lambda1-lambda2,lambda3-lambda4) )  
             enddo
            enddo
           enddo
          enddo
         enddo
        enddo
!!        
      else   ! non unitarization used   
      
        CALL vvvv_wigner_final(OFFepsz(1,-1,1), OFFepsw(1,-1,1),q3(0),q4(0),Aoff)   !Voff(1,2) -> Von(3,4)  


        do block1 = 1,11
         do block2 = 1,11
          do lambda1 = -1,1
           do lambda2 = -1,1
            do lambda3 = -1,1
             do lambda4 = -1,1
               M(lambda1,lambda2,lambda3,lambda4,block1,block2) = sameprod(block2)*8.0d0*rPi*( Aoff(lambda1,lambda2,lambda3,lambda4,block1,block2,0)*smalld(0,lambda1-lambda2,lambda3-lambda4) 
     &          + 3.0d0*Aoff(lambda1,lambda2,lambda3,lambda4,block1,block2,1)*smalld(1,lambda1-lambda2,lambda3-lambda4) 
     &          + 5.0d0*Aoff(lambda1,lambda2,lambda3,lambda4,block1,block2,2)*smalld(2,lambda1-lambda2,lambda3-lambda4) )

               
             enddo
            enddo
           enddo
          enddo
         enddo
        enddo


        endif
 
            
c     construction of the leptonic tensor from the individual parts according to process number
        
      select case (procnum)
      case (200)
        
          do lambda = -1,1
           u3wp(lambda) = dotcc(OFFepspw(1,lambda,3),Wpb(1))
           u4wm(lambda) = dotcc(OFFepspw(1,lambda,4),Wmb(1))
          enddo     

          do block1 =1,11
            do lambda1 = -1,1
              do lambda2 = -1,1
                do lambda3 = -1,1
                  do lambda = -1,1
                    M1(lambda1,lambda2,lambda3,block1,11) = M1(lambda1,lambda2,lambda3,block1,11) +
     &                                                      M(lambda1,lambda2,lambda3,lambda,block1,11)*u4wm(lambda)  
        
                  enddo            
               enddo
              enddo
            enddo
          enddo

          do block1 = 1,11
            do lambda1 = -1,1
              do lambda2 = -1,1
                do lambda = -1,1
                  M2(lambda1,lambda2,block1,11) = M2(lambda1,lambda2,block1,11) +
     &                                            M1(lambda1,lambda2,lambda,block1,11)*u3wp(lambda)      
 
                enddo            
              enddo
            enddo
          enddo
        
          D1Z = rOne/dcmplx(Sq1-ZMASS**2, ZMASS*ZWIDTH) 
          D2Z = rOne/dcmplx(Sq2-ZMASS**2, ZMASS*ZWIDTH)
          D1W = rOne/dcmplx(Sq1-WMASS**2, WMASS*WWIDTH) 
          D2W = rOne/dcmplx(Sq2-WMASS**2, WMASS*WWIDTH)
          D1A = rOne/dcmplx(Sq1, rZero) 
          D2A = rOne/dcmplx(Sq2, rZero)     
        
          do nu = 0,3
            
              do lambda1 = -1,1
                do lambda = -1,1
                  MM(lambda1,nu,1,11) = MM(lambda1,nu,1,11) + (-1)**(lambda+1)
     &                         *d2a*OFFepsC2z(nu+1,lambda,2)*M2(lambda1,lambda,1,11) 
                  MM(lambda1,nu,3,11) = MM(lambda1,nu,3,11) + (-1)**(lambda+1)
     &                         *d2z*OFFepsC2z(nu+1,lambda,2)*M2(lambda1,lambda,3,11)
                  MM(lambda1,nu,7,11) = MM(lambda1,nu,7,11) + (-1)**(lambda+1)
     &                         *d2a*OFFepsC2z(nu+1,lambda,2)*M2(lambda1,lambda,7,11) 
                  MM(lambda1,nu,9,11) = MM(lambda1,nu,9,11) + (-1)**(lambda+1)
     &                         *d2z*OFFepsC2z(nu+1,lambda,2)*M2(lambda1,lambda,9,11)
                  MM(lambda1,nu,10,11) = MM(lambda1,nu,10,11) + (-1)**(lambda+1)
     &                         *d2w*OFFepsC2w(nu+1,lambda,2)*M2(lambda1,lambda,10,11) 
                  MM(lambda1,nu,11,11) = MM(lambda1,nu,11,11) + (-1)**(lambda+1)
     &                         *d2w*OFFepsC2w(nu+1,lambda,2)*M2(lambda1,lambda,11,11)
                enddo
              enddo
            
          enddo
         
          do mu = 0,3
            do nu = 0,3
 
                do lambda1 = -1,1
                  vvt(mu,nu,1,11) = vvt(mu,nu,1,11) + (-1)**(lambda1+1)  
     &                     *d1a*OFFepsC2z(mu+1,lambda1,1)*MM(lambda1,nu,1,11)
                  vvt(mu,nu,3,11) = vvt(mu,nu,3,11) + (-1)**(lambda1+1)
     &                     *d1a*OFFepsC2z(mu+1,lambda1,1)*MM(lambda1,nu,3,11) 
                  vvt(mu,nu,7,11) = vvt(mu,nu,7,11) + (-1)**(lambda1+1)  
     &                     *d1z*OFFepsC2z(mu+1,lambda1,1)*MM(lambda1,nu,7,11)
                  vvt(mu,nu,9,11) = vvt(mu,nu,9,11) + (-1)**(lambda1+1)
     &                     *d1z*OFFepsC2z(mu+1,lambda1,1)*MM(lambda1,nu,9,11) 
                  vvt(mu,nu,10,11) = vvt(mu,nu,10,11) + (-1)**(lambda1+1)  
     &                     *d1w*OFFepsC2w(mu+1,lambda1,1)*MM(lambda1,nu,10,11)
                  vvt(mu,nu,11,11) = vvt(mu,nu,11,11) + (-1)**(lambda1+1)
     &                     *d1w*OFFepsC2w(mu+1,lambda1,1)*MM(lambda1,nu,11,11) 
                enddo

            enddo
          enddo         
          do mu = 0,3
            do nu = 0,3
              aaww_anomal(mu,nu)=vvt(mu,nu,1,11)
              azww_anomal(mu,nu)=vvt(mu,nu,3,11)
              zaww_anomal(mu,nu)=vvt(nu,mu,7,11)
              zzww_anomal(mu,nu)=vvt(mu,nu,9,11)
              
              wwww5_anomal(mu,nu)=vvt(nu,mu,11,11)
              wwww6_anomal(mu,nu)=vvt(mu,nu,10,11)
            enddo
          enddo
            
        

      case (210)
        
          do lambda = -1,1
           u3z(lambda) = dotcc(OFFepspz(1,lambda,3),Zeb(1))
           u4z(lambda) = dotcc(OFFepspz(1,lambda,4),Zub(1))           
           u3a(lambda) = dotcc(OFFepspz(1,lambda,3),Aeb(1))
           u4a(lambda) = dotcc(OFFepspz(1,lambda,4),Aub(1))
          enddo     

          do block1 =1,11

              do lambda1 = -1,1
                do lambda2 = -1,1
                  do lambda3 = -1,1
                    do lambda = -1,1
                      M1(lambda1,lambda2,lambda3,block1,1) = M1(lambda1,lambda2,lambda3,block1,1) +
     &                                                            M(lambda1,lambda2,lambda3,lambda,block1,1)*u4a(lambda) 
                      M1(lambda1,lambda2,lambda3,block1,3) = M1(lambda1,lambda2,lambda3,block1,3) +
     &                                                            M(lambda1,lambda2,lambda3,lambda,block1,3)*u4z(lambda) 
                      M1(lambda1,lambda2,lambda3,block1,7) = M1(lambda1,lambda2,lambda3,block1,7) +
     &                                                            M(lambda1,lambda2,lambda3,lambda,block1,7)*u4a(lambda) 
                      M1(lambda1,lambda2,lambda3,block1,9) = M1(lambda1,lambda2,lambda3,block1,9) +
     &                                                            M(lambda1,lambda2,lambda3,lambda,block1,9)*u4z(lambda)  
        
                    enddo            
                  enddo
                enddo
              enddo

          enddo

          do block1 = 1,11

              do lambda1 = -1,1
                do lambda2 = -1,1
                  do lambda = -1,1
                    M2(lambda1,lambda2,block1,1) = M2(lambda1,lambda2,block1,1) +
     &                                                  M1(lambda1,lambda2,lambda,block1,1)*u3a(lambda)      
                    M2(lambda1,lambda2,block1,3) = M2(lambda1,lambda2,block1,3) +
     &                                                  M1(lambda1,lambda2,lambda,block1,3)*u3a(lambda)  
                    M2(lambda1,lambda2,block1,7) = M2(lambda1,lambda2,block1,7) +
     &                                                  M1(lambda1,lambda2,lambda,block1,7)*u3z(lambda)  
                    M2(lambda1,lambda2,block1,9) = M2(lambda1,lambda2,block1,9) +
     &                                                  M1(lambda1,lambda2,lambda,block1,9)*u3z(lambda)  
  
                  enddo            
                enddo
              enddo

          enddo
        
          D1Z = rOne/dcmplx(Sq1-ZMASS**2, ZMASS*ZWIDTH) 
          D2Z = rOne/dcmplx(Sq2-ZMASS**2, ZMASS*ZWIDTH)
          D1W = rOne/dcmplx(Sq1-WMASS**2, WMASS*WWIDTH) 
          D2W = rOne/dcmplx(Sq2-WMASS**2, WMASS*WWIDTH)
          D1A = rOne/dcmplx(Sq1, rZero) 
          D2A = rOne/dcmplx(Sq2, rZero)     
        
          do nu = 0,3
            do block2 = 1,9
              do lambda1 = -1,1
                do lambda = -1,1
                  MM(lambda1,nu,1,block2) = MM(lambda1,nu,1,block2) + (-1)**(lambda+1)
     &                         *d2a*OFFepsC2z(nu+1,lambda,2)*M2(lambda1,lambda,1,block2) 
                  MM(lambda1,nu,3,block2) = MM(lambda1,nu,3,block2) + (-1)**(lambda+1)
     &                         *d2z*OFFepsC2z(nu+1,lambda,2)*M2(lambda1,lambda,3,block2)
                  MM(lambda1,nu,7,block2) = MM(lambda1,nu,7,block2) + (-1)**(lambda+1)
     &                         *d2a*OFFepsC2z(nu+1,lambda,2)*M2(lambda1,lambda,7,block2) 
                  MM(lambda1,nu,9,block2) = MM(lambda1,nu,9,block2) + (-1)**(lambda+1)
     &                         *d2z*OFFepsC2z(nu+1,lambda,2)*M2(lambda1,lambda,9,block2)
                  MM(lambda1,nu,10,block2) = MM(lambda1,nu,10,block2) + (-1)**(lambda+1)
     &                         *d2w*OFFepsC2w(nu+1,lambda,2)*M2(lambda1,lambda,10,block2) 
                  MM(lambda1,nu,11,block2) = MM(lambda1,nu,11,block2) + (-1)**(lambda+1)
     &                         *d2w*OFFepsC2w(nu+1,lambda,2)*M2(lambda1,lambda,11,block2)
                enddo
              enddo
            enddo
          enddo
         
          do mu = 0,3
            do nu = 0,3
              do block2 = 1,9
                do lambda1 = -1,1
                  vvt(mu,nu,1,block2) = vvt(mu,nu,1,block2) + (-1)**(lambda1+1)  
     &                     *d1a*OFFepsC2z(mu+1,lambda1,1)*MM(lambda1,nu,1,block2)
                  vvt(mu,nu,3,block2) = vvt(mu,nu,3,block2) + (-1)**(lambda1+1)
     &                     *d1a*OFFepsC2z(mu+1,lambda1,1)*MM(lambda1,nu,3,block2) 
                  vvt(mu,nu,7,block2) = vvt(mu,nu,7,block2) + (-1)**(lambda1+1)  
     &                     *d1z*OFFepsC2z(mu+1,lambda1,1)*MM(lambda1,nu,7,block2)
                  vvt(mu,nu,9,block2) = vvt(mu,nu,9,block2) + (-1)**(lambda1+1)
     &                     *d1z*OFFepsC2z(mu+1,lambda1,1)*MM(lambda1,nu,9,block2) 
                  vvt(mu,nu,10,block2) = vvt(mu,nu,10,block2) + (-1)**(lambda1+1)  
     &                     *d1w*OFFepsC2w(mu+1,lambda1,1)*MM(lambda1,nu,10,block2)
                  vvt(mu,nu,11,block2) = vvt(mu,nu,11,block2) + (-1)**(lambda1+1)
     &                     *d1w*OFFepsC2w(mu+1,lambda1,1)*MM(lambda1,nu,11,block2) 
                enddo
              enddo
            enddo
          enddo         
          do mu = 0,3
            do nu = 0,3
              aazz_anomal(mu,nu)=vvt(mu,nu,1,1)+vvt(mu,nu,1,3)+ vvt(mu,nu,1,7)+vvt(mu,nu,1,9)
              azzz_anomal(mu,nu)=vvt(mu,nu,3,1)+vvt(mu,nu,3,3)+ vvt(mu,nu,3,7)+vvt(mu,nu,3,9)
              zazz_anomal(mu,nu)=vvt(nu,mu,7,1)+vvt(nu,mu,7,3)+ vvt(nu,mu,7,7)+vvt(nu,mu,7,9)
              zzzz_anomal(mu,nu)=vvt(mu,nu,9,1)+vvt(mu,nu,9,3)+ vvt(mu,nu,9,7)+vvt(mu,nu,9,9)
              wwzz5_anomal(mu,nu)=vvt(mu,nu,11,1)+vvt(mu,nu,11,3)+ vvt(mu,nu,11,7)+
     &                  vvt(mu,nu,11,9)
              wwzz6_anomal(mu,nu)=vvt(nu,mu,10,1)+vvt(nu,mu,10,3)+ vvt(nu,mu,10,7)+
     &                  vvt(nu,mu,10,9)
            enddo
          enddo  



      case (211)
        
          do lambda = -1,1
           u3z(lambda) = dotcc(OFFepspz(1,lambda,3),Zeb(1))
           u3a(lambda) = dotcc(OFFepspz(1,lambda,3),Aeb(1))           
           u4z(lambda) = dotcc(OFFepspz(1,lambda,4),Zub(1))
          enddo     

          do block1 =1,11

              do lambda1 = -1,1
                do lambda2 = -1,1
                  do lambda3 = -1,1
                    do lambda = -1,1
                      M1(lambda1,lambda2,lambda3,block1,3) = M1(lambda1,lambda2,lambda3,block1,3) +
     &                                                            M(lambda1,lambda2,lambda3,lambda,block1,3)*u4z(lambda) 
                      M1(lambda1,lambda2,lambda3,block1,9) = M1(lambda1,lambda2,lambda3,block1,9) +
     &                                                            M(lambda1,lambda2,lambda3,lambda,block1,9)*u4z(lambda)  
        
                    enddo            
                  enddo
                enddo
              enddo

          enddo

          do block1 = 1,11

              do lambda1 = -1,1
                do lambda2 = -1,1
                  do lambda = -1,1     
                    M2(lambda1,lambda2,block1,3) = M2(lambda1,lambda2,block1,3) +
     &                                                  M1(lambda1,lambda2,lambda,block1,3)*u3a(lambda)   
                    M2(lambda1,lambda2,block1,9) = M2(lambda1,lambda2,block1,9) +
     &                                                  M1(lambda1,lambda2,lambda,block1,9)*u3z(lambda)     
  
                  enddo            
                enddo
              enddo

          enddo
        
          D1Z = rOne/dcmplx(Sq1-ZMASS**2, ZMASS*ZWIDTH) 
          D2Z = rOne/dcmplx(Sq2-ZMASS**2, ZMASS*ZWIDTH)
          D1W = rOne/dcmplx(Sq1-WMASS**2, WMASS*WWIDTH) 
          D2W = rOne/dcmplx(Sq2-WMASS**2, WMASS*WWIDTH)
          D1A = rOne/dcmplx(Sq1, rZero) 
          D2A = rOne/dcmplx(Sq2, rZero)     
        
          do nu = 0,3
            do block2 = 1,9
              do lambda1 = -1,1
                do lambda = -1,1
                  MM(lambda1,nu,1,block2) = MM(lambda1,nu,1,block2) + (-1)**(lambda+1)
     &                         *d2a*OFFepsC2z(nu+1,lambda,2)*M2(lambda1,lambda,1,block2) 
                  MM(lambda1,nu,3,block2) = MM(lambda1,nu,3,block2) + (-1)**(lambda+1)
     &                         *d2z*OFFepsC2z(nu+1,lambda,2)*M2(lambda1,lambda,3,block2)
                  MM(lambda1,nu,7,block2) = MM(lambda1,nu,7,block2) + (-1)**(lambda+1)
     &                         *d2a*OFFepsC2z(nu+1,lambda,2)*M2(lambda1,lambda,7,block2) 
                  MM(lambda1,nu,9,block2) = MM(lambda1,nu,9,block2) + (-1)**(lambda+1)
     &                         *d2z*OFFepsC2z(nu+1,lambda,2)*M2(lambda1,lambda,9,block2)
                  MM(lambda1,nu,10,block2) = MM(lambda1,nu,10,block2) + (-1)**(lambda+1)
     &                         *d2w*OFFepsC2w(nu+1,lambda,2)*M2(lambda1,lambda,10,block2) 
                  MM(lambda1,nu,11,block2) = MM(lambda1,nu,11,block2) + (-1)**(lambda+1)
     &                         *d2w*OFFepsC2w(nu+1,lambda,2)*M2(lambda1,lambda,11,block2)
                enddo
              enddo
            enddo
          enddo
         
          do mu = 0,3
            do nu = 0,3
              do block2 = 1,9
                do lambda1 = -1,1
                  vvt(mu,nu,1,block2) = vvt(mu,nu,1,block2) + (-1)**(lambda1+1)  
     &                     *d1a*OFFepsC2z(mu+1,lambda1,1)*MM(lambda1,nu,1,block2)
                  vvt(mu,nu,3,block2) = vvt(mu,nu,3,block2) + (-1)**(lambda1+1)
     &                     *d1a*OFFepsC2z(mu+1,lambda1,1)*MM(lambda1,nu,3,block2) 
                  vvt(mu,nu,7,block2) = vvt(mu,nu,7,block2) + (-1)**(lambda1+1)  
     &                     *d1z*OFFepsC2z(mu+1,lambda1,1)*MM(lambda1,nu,7,block2)
                  vvt(mu,nu,9,block2) = vvt(mu,nu,9,block2) + (-1)**(lambda1+1)
     &                     *d1z*OFFepsC2z(mu+1,lambda1,1)*MM(lambda1,nu,9,block2) 
                  vvt(mu,nu,10,block2) = vvt(mu,nu,10,block2) + (-1)**(lambda1+1)  
     &                     *d1w*OFFepsC2w(mu+1,lambda1,1)*MM(lambda1,nu,10,block2)
                  vvt(mu,nu,11,block2) = vvt(mu,nu,11,block2) + (-1)**(lambda1+1)
     &                     *d1w*OFFepsC2w(mu+1,lambda1,1)*MM(lambda1,nu,11,block2) 
                enddo
              enddo
            enddo
          enddo         
          do mu = 0,3
            do nu = 0,3
              aazz_anomal(mu,nu)=vvt(mu,nu,1,3)+vvt(mu,nu,1,9)
              azzz_anomal(mu,nu)=vvt(mu,nu,3,3)+vvt(mu,nu,3,9)
              zazz_anomal(mu,nu)=vvt(nu,mu,7,3)+vvt(nu,mu,7,9)
              zzzz_anomal(mu,nu)=vvt(mu,nu,9,3)+vvt(mu,nu,9,9)
              wwzz5_anomal(mu,nu)=vvt(mu,nu,11,3)+vvt(mu,nu,11,9)
              wwzz6_anomal(mu,nu)=vvt(nu,mu,10,3)+vvt(nu,mu,10,9)
            enddo
          enddo  
        

      case (290)
        
          do lambda = -1,1
           u3z(lambda) = dotcc(OFFepspz(1,lambda,3),Zeb(1))
           u3a(lambda) = dotcc(OFFepspz(1,lambda,3),Aeb(1))
           u4a(lambda) = dotcc(OFFepspz(1,lambda,4),eab(1))
          enddo     

          do block1 =1,11

              do lambda1 = -1,1
                do lambda2 = -1,1
                  do lambda3 = -1,1
                    do lambda = -1,1
                      M1(lambda1,lambda2,lambda3,block1,1) = M1(lambda1,lambda2,lambda3,block1,1) +
     &                                                            M(lambda1,lambda2,lambda3,lambda,block1,1)*u4a(lambda)  
                      M1(lambda1,lambda2,lambda3,block1,7) = M1(lambda1,lambda2,lambda3,block1,7) +
     &                                                            M(lambda1,lambda2,lambda3,lambda,block1,7)*u4a(lambda)  
        
                    enddo            
                  enddo
                enddo
              enddo

          enddo

          do block1 = 1,11

              do lambda1 = -1,1
                do lambda2 = -1,1
                  do lambda = -1,1
                    M2(lambda1,lambda2,block1,1) = M2(lambda1,lambda2,block1,1) +
     &                                                  M1(lambda1,lambda2,lambda,block1,1)*u3a(lambda)      
                    M2(lambda1,lambda2,block1,7) = M2(lambda1,lambda2,block1,7) +
     &                                                  M1(lambda1,lambda2,lambda,block1,7)*u3z(lambda)    
  
                  enddo            
                enddo
              enddo

          enddo
        
          D1Z = rOne/dcmplx(Sq1-ZMASS**2, ZMASS*ZWIDTH) 
          D2Z = rOne/dcmplx(Sq2-ZMASS**2, ZMASS*ZWIDTH)
          D1W = rOne/dcmplx(Sq1-WMASS**2, WMASS*WWIDTH) 
          D2W = rOne/dcmplx(Sq2-WMASS**2, WMASS*WWIDTH)
          D1A = rOne/dcmplx(Sq1, rZero) 
          D2A = rOne/dcmplx(Sq2, rZero)     
        
          do nu = 0,3
            do block2 = 1,7
              do lambda1 = -1,1
                do lambda = -1,1
                  MM(lambda1,nu,1,block2) = MM(lambda1,nu,1,block2) + (-1)**(lambda+1)
     &                         *d2a*OFFepsC2z(nu+1,lambda,2)*M2(lambda1,lambda,1,block2) 
                  MM(lambda1,nu,3,block2) = MM(lambda1,nu,3,block2) + (-1)**(lambda+1)
     &                         *d2z*OFFepsC2z(nu+1,lambda,2)*M2(lambda1,lambda,3,block2)
                  MM(lambda1,nu,7,block2) = MM(lambda1,nu,7,block2) + (-1)**(lambda+1)
     &                         *d2a*OFFepsC2z(nu+1,lambda,2)*M2(lambda1,lambda,7,block2) 
                  MM(lambda1,nu,9,block2) = MM(lambda1,nu,9,block2) + (-1)**(lambda+1)
     &                         *d2z*OFFepsC2z(nu+1,lambda,2)*M2(lambda1,lambda,9,block2)
                  MM(lambda1,nu,10,block2) = MM(lambda1,nu,10,block2) + (-1)**(lambda+1)
     &                         *d2w*OFFepsC2w(nu+1,lambda,2)*M2(lambda1,lambda,10,block2) 
                  MM(lambda1,nu,11,block2) = MM(lambda1,nu,11,block2) + (-1)**(lambda+1)
     &                         *d2w*OFFepsC2w(nu+1,lambda,2)*M2(lambda1,lambda,11,block2)
                enddo
              enddo
            enddo
          enddo
         
          do mu = 0,3
            do nu = 0,3
              do block2 = 1,7
                do lambda1 = -1,1
                  vvt(mu,nu,1,block2) = vvt(mu,nu,1,block2) + (-1)**(lambda1+1)  
     &                     *d1a*OFFepsC2z(mu+1,lambda1,1)*MM(lambda1,nu,1,block2)
                  vvt(mu,nu,3,block2) = vvt(mu,nu,3,block2) + (-1)**(lambda1+1)
     &                     *d1a*OFFepsC2z(mu+1,lambda1,1)*MM(lambda1,nu,3,block2) 
                  vvt(mu,nu,7,block2) = vvt(mu,nu,7,block2) + (-1)**(lambda1+1)  
     &                     *d1z*OFFepsC2z(mu+1,lambda1,1)*MM(lambda1,nu,7,block2)
                  vvt(mu,nu,9,block2) = vvt(mu,nu,9,block2) + (-1)**(lambda1+1)
     &                     *d1z*OFFepsC2z(mu+1,lambda1,1)*MM(lambda1,nu,9,block2) 
                  vvt(mu,nu,10,block2) = vvt(mu,nu,10,block2) + (-1)**(lambda1+1)  
     &                     *d1w*OFFepsC2w(mu+1,lambda1,1)*MM(lambda1,nu,10,block2)
                  vvt(mu,nu,11,block2) = vvt(mu,nu,11,block2) + (-1)**(lambda1+1)
     &                     *d1w*OFFepsC2w(mu+1,lambda1,1)*MM(lambda1,nu,11,block2) 
                enddo
              enddo
            enddo
          enddo         
          do mu = 0,3
            do nu = 0,3
              aazz_anomal(mu,nu)=vvt(mu,nu,1,1)+vvt(mu,nu,1,7)
              azzz_anomal(mu,nu)=vvt(mu,nu,3,1)+vvt(mu,nu,3,7)
              zazz_anomal(mu,nu)=vvt(nu,mu,7,1)+vvt(nu,mu,7,7)
              zzzz_anomal(mu,nu)=vvt(mu,nu,9,1)+vvt(mu,nu,9,7)
              wwzz5_anomal(mu,nu)=vvt(mu,nu,11,1)+vvt(mu,nu,11,7)
              wwzz6_anomal(mu,nu)=vvt(nu,mu,10,1)+vvt(nu,mu,10,7)
            enddo
          enddo  
        


      case (291)
        
          do lambda = -1,1
           u3z(lambda) = dotcc(OFFepspz(1,lambda,3),Zeb(1))
           u4a(lambda) = dotcc(OFFepspz(1,lambda,4),eab(1))
          enddo     

          do block1 =1,11

              do lambda1 = -1,1
                do lambda2 = -1,1
                  do lambda3 = -1,1
                    do lambda = -1,1
                      M1(lambda1,lambda2,lambda3,block1,7) = M1(lambda1,lambda2,lambda3,block1,7) +
     &                                                            M(lambda1,lambda2,lambda3,lambda,block1,7)*u4a(lambda)  
        
                    enddo            
                  enddo
                enddo
              enddo

          enddo

          do block1 = 1,11

              do lambda1 = -1,1
                do lambda2 = -1,1
                  do lambda = -1,1
                    M2(lambda1,lambda2,block1,7) = M2(lambda1,lambda2,block1,7) +
     &                                                  M1(lambda1,lambda2,lambda,block1,7)*u3z(lambda)      
  
                  enddo            
                enddo
              enddo

          enddo
        
          D1Z = rOne/dcmplx(Sq1-ZMASS**2, ZMASS*ZWIDTH) 
          D2Z = rOne/dcmplx(Sq2-ZMASS**2, ZMASS*ZWIDTH)
          D1W = rOne/dcmplx(Sq1-WMASS**2, WMASS*WWIDTH) 
          D2W = rOne/dcmplx(Sq2-WMASS**2, WMASS*WWIDTH)
          D1A = rOne/dcmplx(Sq1, rZero) 
          D2A = rOne/dcmplx(Sq2, rZero)     
        
          do nu = 0,3
            do block2 = 1,7
              do lambda1 = -1,1
                do lambda = -1,1
                  MM(lambda1,nu,1,block2) = MM(lambda1,nu,1,block2) + (-1)**(lambda+1)
     &                         *d2a*OFFepsC2z(nu+1,lambda,2)*M2(lambda1,lambda,1,block2) 
                  MM(lambda1,nu,3,block2) = MM(lambda1,nu,3,block2) + (-1)**(lambda+1)
     &                         *d2z*OFFepsC2z(nu+1,lambda,2)*M2(lambda1,lambda,3,block2)
                  MM(lambda1,nu,7,block2) = MM(lambda1,nu,7,block2) + (-1)**(lambda+1)
     &                         *d2a*OFFepsC2z(nu+1,lambda,2)*M2(lambda1,lambda,7,block2) 
                  MM(lambda1,nu,9,block2) = MM(lambda1,nu,9,block2) + (-1)**(lambda+1)
     &                         *d2z*OFFepsC2z(nu+1,lambda,2)*M2(lambda1,lambda,9,block2)
                  MM(lambda1,nu,10,block2) = MM(lambda1,nu,10,block2) + (-1)**(lambda+1)
     &                         *d2w*OFFepsC2w(nu+1,lambda,2)*M2(lambda1,lambda,10,block2) 
                  MM(lambda1,nu,11,block2) = MM(lambda1,nu,11,block2) + (-1)**(lambda+1)
     &                         *d2w*OFFepsC2w(nu+1,lambda,2)*M2(lambda1,lambda,11,block2)
                enddo
              enddo
            enddo
          enddo
         
          do mu = 0,3
            do nu = 0,3
              do block2 = 1,7
                do lambda1 = -1,1
                  vvt(mu,nu,1,block2) = vvt(mu,nu,1,block2) + (-1)**(lambda1+1)  
     &                     *d1a*OFFepsC2z(mu+1,lambda1,1)*MM(lambda1,nu,1,block2)
                  vvt(mu,nu,3,block2) = vvt(mu,nu,3,block2) + (-1)**(lambda1+1)
     &                     *d1a*OFFepsC2z(mu+1,lambda1,1)*MM(lambda1,nu,3,block2) 
                  vvt(mu,nu,7,block2) = vvt(mu,nu,7,block2) + (-1)**(lambda1+1)  
     &                     *d1z*OFFepsC2z(mu+1,lambda1,1)*MM(lambda1,nu,7,block2)
                  vvt(mu,nu,9,block2) = vvt(mu,nu,9,block2) + (-1)**(lambda1+1)
     &                     *d1z*OFFepsC2z(mu+1,lambda1,1)*MM(lambda1,nu,9,block2) 
                  vvt(mu,nu,10,block2) = vvt(mu,nu,10,block2) + (-1)**(lambda1+1)  
     &                     *d1w*OFFepsC2w(mu+1,lambda1,1)*MM(lambda1,nu,10,block2)
                  vvt(mu,nu,11,block2) = vvt(mu,nu,11,block2) + (-1)**(lambda1+1)
     &                     *d1w*OFFepsC2w(mu+1,lambda1,1)*MM(lambda1,nu,11,block2) 
                enddo
              enddo
            enddo
          enddo         
          do mu = 0,3
            do nu = 0,3
              aazz_anomal(mu,nu)=vvt(mu,nu,1,7)
              azzz_anomal(mu,nu)=vvt(mu,nu,3,7)
              zazz_anomal(mu,nu)=vvt(nu,mu,7,7)
              zzzz_anomal(mu,nu)=vvt(mu,nu,9,7)
              wwzz5_anomal(mu,nu)=vvt(mu,nu,11,7)
              wwzz6_anomal(mu,nu)=vvt(nu,mu,10,7)
            enddo
          enddo  
        

      end select   
      END      





