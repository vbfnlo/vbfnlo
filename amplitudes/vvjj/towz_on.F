cccc Partial wave decomposition subroutine      

      SUBROUTINE wpvtowpv_Wigner(Qp1,Qp2,Qp3,Qp4,jNLO, avT,zvT)
c This routine calculates WW->(WW)4l
c using polarization vectors and Wigner decomposition.
c
c IN:
c     Qp1, Qp2          : momenta of the incoming vector bosons in the LAB frame
c     Qp3, Qp4          : momenta of the outgoing vector bosons in the LAB frame
c OUT:    
c     avT          : AW+ -> e+ ve mu+ mu- tensor
c     zvT          : ZW+ -> e+ ve mu+ mu- tensor
c
      IMPLICIT NONE

      INTEGER NGRAPHS
      PARAMETER( NGRAPHS = 20)
      DOUBLE PRECISION rZERO, rONE, rPi
      PARAMETER(rZERO = 0.0D0, rONE = 1.0D0, rPi = 3.14159265358979323846d0)
      DOUBLE COMPLEX cIMAG, cZERO, cONE
      PARAMETER( cIMAG=(0.0D0,1.0D0), cZERO=(0.0D0,0.0D0), cONE=(1.0D0,0.0D0))

      DOUBLE PRECISION QP1(0:3),QP2(0:3),QP3(0:3),QP4(0:3)
      DOUBLE COMPLEX avT(0:3,0:3),zvT(0:3,0:3)    
      INTEGER jNLO
      
      DOUBLE PRECISION Q1(0:3), Q2(0:3), Q3(0:3), Q4(0:3) 
      DOUBLE PRECISION Q1ON(0:3), Q2ON(0:3)
      DOUBLE PRECISION smalld(0:2,-2:2,-2:2)

      DOUBLE PRECISION Q1T, Q2T, Q3T, Q4T        
      DOUBLE PRECISION Sq1, Sq2, Sq3, Sq4
      DOUBLE PRECISION Ecm_VV
      DOUBLE PRECISION COST
      
      DOUBLE COMPLEX D1, D2A, D2Z, D3, D4
      DOUBLE PRECISION Nqw,Nqz, Npw(4), Npz(4)
      
      DOUBLE COMPLEX WEB(6), AUB(6), ZUB(6)
      DOUBLE COMPLEX OFFepsw(6,-1:1,4),OFFepsz(6,-1:1,4), OFFepspz(6,-1:1,4) , OFFepspw(6,-1:1,4)  !ToDo
      DOUBLE COMPLEX ONONepsw(6,-1:1,NGRAPHS), ONONepsz(6,-1:1,NGRAPHS)   !ToDo
      DOUBLE COMPLEX OFFepsC(6,-1:1,2), OFFepsC2(6,-1:1,2) !ToDo
      DOUBLE COMPLEX U3(-1:1), U4a(-1:1), U4z(-1:1)
      DOUBLE COMPLEX M(-1:1,-1:1,-1:1,-1:1,1:3,1:3),M1(-1:1,-1:1,-1:1,1:3,1:3)
      DOUBLE COMPLEX MM(-1:1,0:3,1:3,1:3),M2(-1:1,-1:1,1:3,1:3) 
      DOUBLE COMPLEX Aoff(-1:1,-1:1,-1:1,-1:1,1:3,1:3,0:2)
      DOUBLE COMPLEX Aoo(-1:1,-1:1,-1:1,-1:1,1:3,1:3,0:2)
      DOUBLE COMPLEX U(-1:1,-1:1,-1:1,-1:1,1:3,1:3,0:2)
      DOUBLE COMPLEX vvt(0:3,0:3,1:3,1:3)


      INTEGER I,J,K,N,MU,NU
      INTEGER LAMBDA, LAMBDA1, LAMBDA2, LAMBDA3, LAMBDA4, BLOCK1, BLOCK2
      logical tunitarize3
      DOUBLE COMPLEX DOTCC
      EXTERNAL DOTCC
           
#include "VBFNLO/utilities/global.inc"

#include "gmunu.inc"
#include "transm.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"


#include "VBFNLO/utilities/lha.inc"

#include "tensorw.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/process.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"


      tunitarize3=.true.
      Nqw=rZERO
      Nqz=rZERO
      Npw(1:4) = rZero
      Npz(1:4) = rZero
      smalld(0:2,-2:2,-2:2) = rZero
      M(-1:1,-1:1,-1:1,-1:1,1:3,1:3) = cZero
      Aoff(-1:1,-1:1,-1:1,-1:1,1:3,1:3,0:2) = cZero
      Aoo(-1:1,-1:1,-1:1,-1:1,1:3,1:3,0:2) = cZero
      U(-1:1,-1:1,-1:1,-1:1,1:3,1:3,0:2) = cZero
      M1(-1:1,-1:1,-1:1,1:3,1:3) = cZero
      M2(-1:1,-1:1,1:3,1:3) = cZero
      MM(-1:1,0:3,1:3,1:3) = cZero
      vvt(0:3,0:3,1:3,1:3) = cZero 
      avt(0:3,0:3)=cZero 
      zvt(0:3,0:3)=cZero
      ecm_real = rZERO
      Ecm_VV = rZERO
      COST=rZERO

ccccccccccccccccccccccccc      
      call BOOSTCM(Qp1(0),Qp2(0),Qp3(0),Qp4(0), TT(0,0,jNLO),TTinv(0,0,jNLO))       
      ecm_real = dsqrt(qp1(0)+qp2(0))**2-((qp1(1)+qp2(1))**2+(qp1(2)+qp2(2))**2+(qp1(3)+qp2(3))**2)
      
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
      
      do i = 0,3               
         web(i+1) = TT(i,0,jNLO)*we(1)-TT(i,1,jNLO)*we(2)-TT(i,2,jNLO)*we(3)-TT(i,3,jNLO)*we(4)
         aub(i+1) = TT(i,0,jNLO)*au(1)-TT(i,1,jNLO)*au(2)-TT(i,2,jNLO)*au(3)-TT(i,3,jNLO)*au(4)
         zub(i+1) = TT(i,0,jNLO)*zu(1)-TT(i,1,jNLO)*zu(2)-TT(i,2,jNLO)*zu(3)-TT(i,3,jNLO)*zu(4)         
      enddo
      web(5) = dcmplx(q3(0),q3(3))
      web(6) = dcmplx(q3(1),q3(2))      
      aub(5) = dcmplx(q4(0),q4(3))
      aub(6) = dcmplx(q4(1),q4(2)) 
      zub(5) = dcmplx(q4(0),q4(3))
      zub(6) = dcmplx(q4(1),q4(2))                  

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



      
    

       call voffxx(q1,WMASS,+1,Npw(1),OFFepsC(1,-1,1)) 
       call voffxx(q2,ZMASS,+1,Npz(2),OFFepsC(1,-1,2)) 

**     Boosted back to Lab Frame Eps1 & Eps2 (incoming vector bosons)

       do i = 1,2
        do lambda = -1,1
           OFFepsC2(1,lambda,i) = TTinv(0,0,jNLO)*OFFepsC(1,lambda,i)+TTinv(0,1,jNLO)*OFFepsC(2,lambda,i)
     &                         + TTinv(0,2,jNLO)*OFFepsC(3,lambda,i)+TTinv(0,3,jNLO)*OFFepsC(4,lambda,i)  
   
           do k = 1,3
             OFFepsC2(k+1,lambda,i) = -TTinv(k,0,jNLO)*OFFepsC(1,lambda,i)-TTinv(k,1,jNLO)*OFFepsC(2,lambda,i)
     &                               - TTinv(k,2,jNLO)*OFFepsC(3,lambda,i)-TTinv(k,3,jNLO)*OFFepsC(4,lambda,i)     
      
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

        CALL wpvwpv_wigner_final(OFFepsz(1,-1,1),OFFepsw(1,-1,1),q3,q4,Aoff)   ! Voff(1,2) -> Von(3,4)

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



        CALL wpvwpv_wigner_final(ONONepsz(1,-1,1), ONONepsw(1,-1,1),q3,q4,Aoo)    ! Von(3,4)  -> Von(3,4)   
        
        ! T-matriz Unitarization calculation
        CALL tmatrixvv(Aoff,Aoo, Ecm_VV,Sq1,Sq2,Sq3,Sq4,21,3, U)   


          



!         M(lambda1,lambda2,lambda3,lambda4, block1, block2), lambdai=-1,1  , block accounts for the different 2->2 processes with the same total charge                       
        do lambda1 = -1,1
         do lambda2 = -1,1
          do lambda3 = -1,1
           do lambda4 = -1,1
            do block1 = 1,3
             do block2 = 1,3
             
              M(lambda1,lambda2,lambda3,lambda4,block1,block2) = 8.0d0*rPi*(U(lambda1,lambda2,lambda3,lambda4,block1,block2,0)*smalld(0,lambda1-lambda2,lambda3-lambda4) 
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
        CALL wpvwpv_wigner_final(OFFepsz(1,-1,1),OFFepsw(1,-1,1),q3,q4,Aoff)   !Voff(1,2) -> Von(3,4)  

 
        do block1 = 1,3
         do block2 = 1,3
          do lambda1 = -1,1
           do lambda2 = -1,1
            do lambda3 = -1,1
             do lambda4 = -1,1
               M(lambda1,lambda2,lambda3,lambda4,block1,block2) = 8.0d0*rPi*( Aoff(lambda1,lambda2,lambda3,lambda4,block1,block2,0)*smalld(0,lambda1-lambda2,lambda3-lambda4) 
     &          + 3.0d0*Aoff(lambda1,lambda2,lambda3,lambda4,block1,block2,1)*smalld(1,lambda1-lambda2,lambda3-lambda4) 
     &          + 5.0d0*Aoff(lambda1,lambda2,lambda3,lambda4,block1,block2,2)*smalld(2,lambda1-lambda2,lambda3-lambda4) )
              
             enddo
            enddo
           enddo
          enddo
         enddo
        enddo
      endif
        


c     construction of the leptonic tensor from the individual parts 
        

          do lambda = -1,1
           u3(lambda) = dotcc(OFFepspw(1,lambda,3),Web(1))
           u4a(lambda) = dotcc(OFFepspz(1,lambda,4),aub(1))
           u4z(lambda) = dotcc(OFFepspz(1,lambda,4),zub(1)) 

          enddo     

        
      do block1 =1,3
        do lambda1 = -1,1
          do lambda2 = -1,1
            do lambda3 = -1,1
              do lambda = -1,1
                M1(lambda1,lambda2,lambda3,block1,3) = M1(lambda1,lambda2,lambda3,block1,3) +
     &                                        M(lambda1,lambda2,lambda3,lambda,block1,3)*u4z(lambda)  
                M1(lambda1,lambda2,lambda3,block1,1) = M1(lambda1,lambda2,lambda3,block1,1) +
     &                                        M(lambda1,lambda2,lambda3,lambda,block1,1)*u4a(lambda)          
              enddo            
            enddo
          enddo
        enddo
      enddo

      do block1 = 1,3
        do lambda1 = -1,1
          do lambda2 = -1,1
            do lambda = -1,1
              M2(lambda1,lambda2,block1,1) = M2(lambda1,lambda2,block1,1) +
     &                              M1(lambda1,lambda2,lambda,block1,1)*u3(lambda)      
              M2(lambda1,lambda2,block1,3) = M2(lambda1,lambda2,block1,3) +
     &                              M1(lambda1,lambda2,lambda,block1,3)*u3(lambda)  
            enddo            
          enddo
        enddo
      enddo
        
      D1 = rOne/dcmplx(Sq1-WMASS**2, WMASS*WWIDTH) 
      D2Z = rOne/dcmplx(Sq2-ZMASS**2, ZMASS*ZWIDTH)   
      D2A = rOne/dcmplx(Sq2, rZero)     

      do nu = 0,3
        
          do lambda1 = -1,1
            do lambda = -1,1
               MM(lambda1,nu,3,1) = MM(lambda1,nu,3,1) + (-1)**(lambda+1)
     &                         *d2z*OFFepsC2(nu+1,lambda,2)*M2(lambda1,lambda,3,1) !
               MM(lambda1,nu,1,1) = MM(lambda1,nu,1,1) + (-1)**(lambda+1)
     &                         *d2a*OFFepsC2(nu+1,lambda,2)*M2(lambda1,lambda,1,1)   
               MM(lambda1,nu,3,3) = MM(lambda1,nu,3,3) + (-1)**(lambda+1)
     &                         *d2z*OFFepsC2(nu+1,lambda,2)*M2(lambda1,lambda,3,3) !
               MM(lambda1,nu,1,3) = MM(lambda1,nu,1,3) + (-1)**(lambda+1)
     &                         *d2a*OFFepsC2(nu+1,lambda,2)*M2(lambda1,lambda,1,3)  
            enddo
          enddo
        
      enddo
         
      do mu = 0,3
        do nu = 0,3
          
            do lambda1 = -1,1
                vvt(mu,nu,1,1) = vvt(mu,nu,1,1) + (-1)**(lambda1+1)  
     &                     *d1*OFFepsC2(mu+1,lambda1,1)*MM(lambda1,nu,1,1)
                vvt(mu,nu,3,1) = vvt(mu,nu,3,1) + (-1)**(lambda1+1)
     &                     *d1*OFFepsC2(mu+1,lambda1,1)*MM(lambda1,nu,3,1) 
                vvt(mu,nu,1,3) = vvt(mu,nu,1,3) + (-1)**(lambda1+1)  
     &                     *d1*OFFepsC2(mu+1,lambda1,1)*MM(lambda1,nu,1,3)
                vvt(mu,nu,3,3) = vvt(mu,nu,3,3) + (-1)**(lambda1+1)
     &                     *d1*OFFepsC2(mu+1,lambda1,1)*MM(lambda1,nu,3,3)      
            enddo
          
        enddo
      enddo  
     
      do mu = 0,3
        do nu = 0,3
            avt(mu,nu)=vvt(nu,mu,1,1)+ vvt(nu,mu,1,3)
            zvt(mu,nu)=vvt(nu,mu,3,1)+ vvt(nu,mu,3,3)

        enddo
      enddo  

      END      





cccc Partial wave decomposition subroutine      

      SUBROUTINE wmvtowmv_Wigner_final(Qp1,Qp2,Qp3,Qp4,jNLO, avT,zvT)
c This routine calculates WW->(WW)4l
c using polarization vectors and Wigner decomposition.
c
c IN:
c     Qp1, Qp2          : momenta of the incoming vector bosons in the LAB frame
c     Qp3, Qp4          : momenta of the outgoing vector bosons in the LAB frame
c OUT:    
c     avT          : AW+ -> e+ ve mu+ mu- tensor
c     zvT          : ZW+ -> e+ ve mu+ mu- tensor
c
      IMPLICIT NONE

      INTEGER NGRAPHS
      PARAMETER( NGRAPHS = 20)
      DOUBLE PRECISION rZERO, rONE, rPi
      PARAMETER(rZERO = 0.0D0, rONE = 1.0D0, rPi = 3.14159265358979323846d0)
      DOUBLE COMPLEX cIMAG, cZERO, cONE
      PARAMETER( cIMAG=(0.0D0,1.0D0), cZERO=(0.0D0,0.0D0), cONE=(1.0D0,0.0D0))

      DOUBLE PRECISION QP1(0:3),QP2(0:3),QP3(0:3),QP4(0:3)
      DOUBLE COMPLEX avT(0:3,0:3),zvT(0:3,0:3)    
      INTEGER jNLO
      
      DOUBLE PRECISION Q1(0:3), Q2(0:3), Q3(0:3), Q4(0:3) 
      DOUBLE PRECISION Q1ON(0:3), Q2ON(0:3)
      DOUBLE PRECISION smalld(0:2,-2:2,-2:2)

      DOUBLE PRECISION Q1T, Q2T, Q3T, Q4T        
      DOUBLE PRECISION Sq1, Sq2, Sq3, Sq4
      DOUBLE PRECISION Ecm_VV
      DOUBLE PRECISION COST
      
      DOUBLE COMPLEX D1, D2A, D2Z, D3, D4
      DOUBLE PRECISION Nqw,Nqz, Npw(4), Npz(4)
      
      DOUBLE COMPLEX WEB(6), AUB(6), ZUB(6)
      DOUBLE COMPLEX OFFepsw(6,-1:1,4),OFFepsz(6,-1:1,4), OFFepspz(6,-1:1,4) , OFFepspw(6,-1:1,4)  !ToDo
      DOUBLE COMPLEX ONONepsw(6,-1:1,NGRAPHS), ONONepsz(6,-1:1,NGRAPHS)   !ToDo
      DOUBLE COMPLEX OFFepsC(6,-1:1,2), OFFepsC2(6,-1:1,2) !ToDo
      DOUBLE COMPLEX U3(-1:1), U4a(-1:1), U4z(-1:1)
      DOUBLE COMPLEX M(-1:1,-1:1,-1:1,-1:1,1:3,1:3),M1(-1:1,-1:1,-1:1,1:3,1:3)
      DOUBLE COMPLEX MM(-1:1,0:3,1:3,1:3),M2(-1:1,-1:1,1:3,1:3) 
      DOUBLE COMPLEX Aoff(-1:1,-1:1,-1:1,-1:1,1:3,1:3,0:2)
      DOUBLE COMPLEX Aoo(-1:1,-1:1,-1:1,-1:1,1:3,1:3,0:2)
      DOUBLE COMPLEX U(-1:1,-1:1,-1:1,-1:1,1:3,1:3,0:2)
      DOUBLE COMPLEX vvt(0:3,0:3,1:3,1:3)


      INTEGER I,J,K,N,MU,NU
      INTEGER LAMBDA, LAMBDA1, LAMBDA2, LAMBDA3, LAMBDA4, BLOCK1, BLOCK2
      logical tunitarize3
      DOUBLE COMPLEX DOTCC
      EXTERNAL DOTCC
           
#include "VBFNLO/utilities/global.inc"

#include "gmunu.inc"
#include "transm.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"


#include "VBFNLO/utilities/lha.inc"

#include "tensorw.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/process.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"



      tunitarize3=.true.
      Nqw=rZERO
      Nqz=rZERO
      Npw(1:4) = rZero
      Npz(1:4) = rZero
      smalld(0:2,-2:2,-2:2) = rZero
      M(-1:1,-1:1,-1:1,-1:1,1:3,1:3) = cZero
      Aoff(-1:1,-1:1,-1:1,-1:1,1:3,1:3,0:2) = cZero
      Aoo(-1:1,-1:1,-1:1,-1:1,1:3,1:3,0:2) = cZero
      U(-1:1,-1:1,-1:1,-1:1,1:3,1:3,0:2) = cZero
      M1(-1:1,-1:1,-1:1,1:3,1:3) = cZero
      M2(-1:1,-1:1,1:3,1:3) = cZero
      MM(-1:1,0:3,1:3,1:3) = cZero
      vvt(0:3,0:3,1:3,1:3) = cZero 
      avt(0:3,0:3)=cZero 
      zvt(0:3,0:3)=cZero
      ecm_real = rZERO
      Ecm_VV = rZERO
      COST=rZERO

ccccccccccccccccccccccccc      
      call BOOSTCM(Qp1(0),Qp2(0),Qp3(0),Qp4(0), TT(0,0,jNLO),TTinv(0,0,jNLO))       
      ecm_real = dsqrt(qp1(0)+qp2(0))**2-((qp1(1)+qp2(1))**2+(qp1(2)+qp2(2))**2+(qp1(3)+qp2(3))**2)
      
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

      
      do i = 0,3               
         web(i+1) = TT(i,0,jNLO)*we(1)-TT(i,1,jNLO)*we(2)-TT(i,2,jNLO)*we(3)-TT(i,3,jNLO)*we(4)
         aub(i+1) = TT(i,0,jNLO)*au(1)-TT(i,1,jNLO)*au(2)-TT(i,2,jNLO)*au(3)-TT(i,3,jNLO)*au(4)
         zub(i+1) = TT(i,0,jNLO)*zu(1)-TT(i,1,jNLO)*zu(2)-TT(i,2,jNLO)*zu(3)-TT(i,3,jNLO)*zu(4)         
      enddo
      web(5) = dcmplx(q3(0),q3(3))
      web(6) = dcmplx(q3(1),q3(2))      
      aub(5) = dcmplx(q4(0),q4(3))
      aub(6) = dcmplx(q4(1),q4(2)) 
      zub(5) = dcmplx(q4(0),q4(3))
      zub(6) = dcmplx(q4(1),q4(2))                  

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



      
    

       call voffxx(q1,WMASS,+1,Npw(1),OFFepsC(1,-1,1)) 
       call voffxx(q2,ZMASS,+1,Npz(2),OFFepsC(1,-1,2)) 

**     Boosted back to Lab Frame Eps1 & Eps2 (incoming vector bosons)

       do i = 1,2
        do lambda = -1,1
           OFFepsC2(1,lambda,i) = TTinv(0,0,jNLO)*OFFepsC(1,lambda,i)+TTinv(0,1,jNLO)*OFFepsC(2,lambda,i)
     &                         + TTinv(0,2,jNLO)*OFFepsC(3,lambda,i)+TTinv(0,3,jNLO)*OFFepsC(4,lambda,i)  
   
           do k = 1,3
             OFFepsC2(k+1,lambda,i) = -TTinv(k,0,jNLO)*OFFepsC(1,lambda,i)-TTinv(k,1,jNLO)*OFFepsC(2,lambda,i)
     &                               - TTinv(k,2,jNLO)*OFFepsC(3,lambda,i)-TTinv(k,3,jNLO)*OFFepsC(4,lambda,i)     
      
           enddo
        enddo
       enddo      

   
      
**     Main Calculation: M(mu,nu), mu,nu=0,3      

      
      if(.not.isnan(abs(COST)))then
        CALL WIGNERdvv(COST,smalld)   
        else
          smalld(0,0,0)=1
        endif    
 
      if(Tunitarize3) then  ! T matrix unitarization     ToDo: change  
        CALL wmvwmv_wigner_final(OFFepsz(1,-1,1),OFFepsw(1,-1,1),q3,q4,Aoff)   ! Voff(1,2) -> Von(3,4)
        
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


        CALL wmvwmv_wigner_final(ONONepsz(1,-1,1), ONONepsw(1,-1,1),q3,q4,Aoo)    ! Von(3,4)  -> Von(3,4)   
        
        ! T-matriz Unitarization calculation
        CALL tmatrixvv(Aoff,Aoo, Ecm_VV,Sq1,Sq2,Sq3,Sq4,21,3, U)  
        
        

!         M(lambda1,lambda2,lambda3,lambda4, block1, block2), lambdai=-1,1  , block accounts for the different 2->2 processes with the same total charge                           
        do lambda1 = -1,1
         do lambda2 = -1,1
          do lambda3 = -1,1
           do lambda4 = -1,1
            do block1 = 1,3
             do block2 = 1,3
              M(lambda1,lambda2,lambda3,lambda4,block1,block2) = 8.0d0*rPi*(U(lambda1,lambda2,lambda3,lambda4,block1,block2,0)*smalld(0,lambda1-lambda2,lambda3-lambda4) 
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
        CALL wmvwmv_wigner_final(OFFepsz(1,-1,1),OFFepsw(1,-1,1),q3,q4,Aoff)   !Voff(1,2) -> Von(3,4)  


        do block1 = 1,3
         do block2 = 1,3
          do lambda1 = -1,1
           do lambda2 = -1,1
            do lambda3 = -1,1
             do lambda4 = -1,1
               M(lambda1,lambda2,lambda3,lambda4,block1,block2) = 8.0d0*rPi*( Aoff(lambda1,lambda2,lambda3,lambda4,block1,block2,0)*smalld(0,lambda1-lambda2,lambda3-lambda4) 
     &          + 3.0d0*Aoff(lambda1,lambda2,lambda3,lambda4,block1,block2,1)*smalld(1,lambda1-lambda2,lambda3-lambda4) 
     &          + 5.0d0*Aoff(lambda1,lambda2,lambda3,lambda4,block1,block2,2)*smalld(2,lambda1-lambda2,lambda3-lambda4) )
             enddo
            enddo
           enddo
          enddo
         enddo
        enddo
        endif

c     construction of the leptonic tensor from the individual parts 
        

          do lambda = -1,1
           u3(lambda) = dotcc(OFFepspw(1,lambda,3),Web(1))
           u4a(lambda) = dotcc(OFFepspz(1,lambda,4),aub(1))
           u4z(lambda) = dotcc(OFFepspz(1,lambda,4),zub(1)) 
          enddo     

        
      do block1 =1,3
        do lambda1 = -1,1
          do lambda2 = -1,1
            do lambda3 = -1,1
              do lambda = -1,1
                M1(lambda1,lambda2,lambda3,block1,3) = M1(lambda1,lambda2,lambda3,block1,3) +
     &                                        M(lambda1,lambda2,lambda3,lambda,block1,3)*u4z(lambda)  
                M1(lambda1,lambda2,lambda3,block1,1) = M1(lambda1,lambda2,lambda3,block1,1) +
     &                                        M(lambda1,lambda2,lambda3,lambda,block1,1)*u4a(lambda)          
              enddo            
            enddo
          enddo
        enddo
      enddo

      do block1 = 1,3
        do lambda1 = -1,1
          do lambda2 = -1,1
            do lambda = -1,1
              M2(lambda1,lambda2,block1,1) = M2(lambda1,lambda2,block1,1) +
     &                              M1(lambda1,lambda2,lambda,block1,1)*u3(lambda)      
              M2(lambda1,lambda2,block1,3) = M2(lambda1,lambda2,block1,3) +
     &                              M1(lambda1,lambda2,lambda,block1,3)*u3(lambda)  
            enddo            
          enddo
        enddo
      enddo
        
      D1 = rOne/dcmplx(Sq1-WMASS**2, WMASS*WWIDTH) 
      D2Z = rOne/dcmplx(Sq2-ZMASS**2, ZMASS*ZWIDTH)   
      D2A = rOne/dcmplx(Sq2, rZero)     
        
      do nu = 0,3
        
          do lambda1 = -1,1
            do lambda = -1,1
               MM(lambda1,nu,3,1) = MM(lambda1,nu,3,1) + (-1)**(lambda+1)
     &                         *d2z*OFFepsC2(nu+1,lambda,2)*M2(lambda1,lambda,3,1) !
               MM(lambda1,nu,1,1) = MM(lambda1,nu,1,1) + (-1)**(lambda+1)
     &                         *d2a*OFFepsC2(nu+1,lambda,2)*M2(lambda1,lambda,1,1)   
               MM(lambda1,nu,3,3) = MM(lambda1,nu,3,3) + (-1)**(lambda+1)
     &                         *d2z*OFFepsC2(nu+1,lambda,2)*M2(lambda1,lambda,3,3) !
               MM(lambda1,nu,1,3) = MM(lambda1,nu,1,3) + (-1)**(lambda+1)
     &                         *d2a*OFFepsC2(nu+1,lambda,2)*M2(lambda1,lambda,1,3)  
            enddo
          enddo
        
      enddo
         
      do mu = 0,3
        do nu = 0,3
          
            do lambda1 = -1,1
                vvt(mu,nu,1,1) = vvt(mu,nu,1,1) + (-1)**(lambda1+1)  
     &                     *d1*OFFepsC2(mu+1,lambda1,1)*MM(lambda1,nu,1,1)
                vvt(mu,nu,3,1) = vvt(mu,nu,3,1) + (-1)**(lambda1+1)
     &                     *d1*OFFepsC2(mu+1,lambda1,1)*MM(lambda1,nu,3,1) 
                vvt(mu,nu,1,3) = vvt(mu,nu,1,3) + (-1)**(lambda1+1)  
     &                     *d1*OFFepsC2(mu+1,lambda1,1)*MM(lambda1,nu,1,3)
                vvt(mu,nu,3,3) = vvt(mu,nu,3,3) + (-1)**(lambda1+1)
     &                     *d1*OFFepsC2(mu+1,lambda1,1)*MM(lambda1,nu,3,3)      
            enddo
          
        enddo
      enddo         
      do mu = 0,3
        do nu = 0,3
            avt(mu,nu)=vvt(nu,mu,1,1)+ vvt(nu,mu,1,3)
            zvt(mu,nu)=vvt(nu,mu,3,1)+ vvt(nu,mu,3,3)
        enddo
      enddo  
      
  
           
      END    