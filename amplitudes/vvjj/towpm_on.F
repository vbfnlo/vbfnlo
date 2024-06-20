cccc Partial wave decomposition subroutine      

      SUBROUTINE WPPTOWPP_Wigner(Qp1,Qp2,Qp3,Qp4,jNLO, WWT)
c This routine calculates WW->(WW)4l
c using polarization vectors and Wigner decomposition.
c
c IN:
c     Qp1, Qp2     : momenta of the incoming vector bosons in the LAB frame
c     Qp3, Qp4     : momenta of the outgoing vector bosons in the LAB frame
c OUT:    
c     WWT          : w+w+ -> w+(e+ ve) w+(mu+ nmu) tensor
c
      IMPLICIT NONE

      INTEGER NGRAPHS
      PARAMETER( NGRAPHS = 20)
      DOUBLE PRECISION rZERO, rONE, rPi
      PARAMETER(rZERO = 0.0D0, rONE = 1.0D0, rPi = 3.14159265358979323846d0)
      DOUBLE COMPLEX cIMAG, cZERO, cONE
      PARAMETER( cIMAG=(0.0D0,1.0D0), cZERO=(0.0D0,0.0D0), cONE=(1.0D0,0.0D0))

      DOUBLE PRECISION QP1(0:3),QP2(0:3),QP3(0:3),QP4(0:3)
      DOUBLE COMPLEX WWT(0:3,0:3)  
      INTEGER jNLO
      
      DOUBLE PRECISION Q1(0:3), Q2(0:3), Q3(0:3), Q4(0:3)
      DOUBLE PRECISION Q1ON(0:3), Q2ON(0:3)
      DOUBLE PRECISION smalld(0:2,-2:2,-2:2)

      DOUBLE PRECISION Q1T, Q2T, Q3T, Q4T      
      DOUBLE PRECISION Sq1, Sq2, Sq3, Sq4
      DOUBLE PRECISION Ecm_VV
      DOUBLE PRECISION COST
      
      DOUBLE COMPLEX D1, D2, D3, D4
      DOUBLE PRECISION Nq(4), Np(4)
      
      DOUBLE COMPLEX WEB(6), WUB(6)
      DOUBLE COMPLEX OFFeps(6,-1:1,NGRAPHS), OFFepsp(6,-1:1,NGRAPHS)
      DOUBLE COMPLEX OFFeps2(6,-1:1,NGRAPHS), ONONeps(6,-1:1,NGRAPHS)
      DOUBLE COMPLEX OFFepsC(6,-1:1,NGRAPHS), OFFepsC2(6,-1:1,NGRAPHS)
      DOUBLE COMPLEX U3(-1:1), U4(-1:1)
      DOUBLE COMPLEX M(-1:1,-1:1,-1:1,-1:1),M1(-1:1,-1:1,-1:1)
      DOUBLE COMPLEX MM(-1:1,0:3),M2(-1:1,-1:1) 
      DOUBLE COMPLEX Aoff(-1:1,-1:1,-1:1,-1:1,1:1,1:1,0:2)
      DOUBLE COMPLEX Aoo(-1:1,-1:1,-1:1,-1:1,1:1,1:1,0:2)
      DOUBLE COMPLEX U(-1:1,-1:1,-1:1,-1:1,1:1,1:1,0:2)

      INTEGER I,J,K,N,MU,NU
      INTEGER LAMBDA, LAMBDA1, LAMBDA2, LAMBDA3, LAMBDA4
      logical tunitarize3
      DOUBLE COMPLEX DOTCC
      EXTERNAL DOTCC
           
#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "gmunu.inc"
#include "transm.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"

      tunitarize3=.true.
      Nq(1:4) = rZero
      Np(1:4) = rZero
      smalld(0:2,-2:2,-2:2) = rZero
      M(-1:1,-1:1,-1:1,-1:1) = cZero
      Aoff(-1:1,-1:1,-1:1,-1:1,1:1,1:1,0:2) = cZero
      Aoo(-1:1,-1:1,-1:1,-1:1,1:1,1:1,0:2) = cZero
      U(-1:1,-1:1,-1:1,-1:1,1:1,1:1,0:2) = cZero
      M1(-1:1,-1:1,-1:1) = cZero
      M2(-1:1,-1:1) = cZero
      MM(-1:1,0:3) = cZero
      wwt(0:3,0:3) = cZero 
      
      ecm_real = rZERO
      Ecm_VV = rZERO

ccccccccccccccccccccccccc      
      
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

      do i = 0,3               
         web(i+1) = TT(i,0,jNLO)*we(1)-TT(i,1,jNLO)*we(2)-TT(i,2,jNLO)*we(3)-TT(i,3,jNLO)*we(4)
         wub(i+1) = TT(i,0,jNLO)*wu(1)-TT(i,1,jNLO)*wu(2)-TT(i,2,jNLO)*wu(3)-TT(i,3,jNLO)*wu(4)
      enddo
      web(5) = dcmplx(q3(0),q3(3))
      web(6) = dcmplx(q3(1),q3(2))      
      wub(5) = dcmplx(q4(0),q4(3))
      wub(6) = dcmplx(q4(1),q4(2))            

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
       Nq(1) = rONE/WMASS
       Nq(2) = rONE/WMASS
       Nq(3) = rONE/WMASS
       Nq(4) = rONE/WMASS
       call voffxx(q1,WMASS,-1,Nq(1),OFFeps(1,-1,1)) 
       call voffxx(q2,WMASS,-1,Nq(2),OFFeps(1,-1,2)) 

c     normalization for longitudinal EXTERNAL polarization vectors , off-shell photon with MZ/Q      
       Np(1) = WMASS/(dabs(Sq1))
       Np(2) = WMASS/(dabs(Sq2))
       Np(3) = WMASS/(dabs(Sq3))
       Np(4) = WMASS/(dabs(Sq4))
       call voffxx(q1,WMASS,-1,Np(1),OFFepsp(1,-1,1)) 
       call voffxx(q2,WMASS,-1,Np(2),OFFepsp(1,-1,2)) 
       call voffxx(q3,WMASS,-1,Np(3),OFFepsp(1,-1,3)) 
       call voffxx(q4,WMASS,-1,Np(4),OFFepsp(1,-1,4))   
**     Boosted back to Lab Frame Eps1 & Eps2 (incoming vector bosons)
       do i = 1,2
        do lambda = -1,1
           OFFeps2(1,lambda,i) = TTinv(0,0,jNLO)*OFFepsp(1,lambda,i)+TTinv(0,1,jNLO)*OFFepsp(2,lambda,i)
     &                         + TTinv(0,2,jNLO)*OFFepsp(3,lambda,i)+TTinv(0,3,jNLO)*OFFepsp(4,lambda,i)    
           do k = 1,3
             OFFeps2(k+1,lambda,i) = -TTinv(k,0,jNLO)*OFFepsp(1,lambda,i)-TTinv(k,1,jNLO)*OFFepsp(2,lambda,i)
     &                               - TTinv(k,2,jNLO)*OFFepsp(3,lambda,i)-TTinv(k,3,jNLO)*OFFepsp(4,lambda,i)     
           enddo
        enddo
       enddo      

      
      do k = 1,2
       do i = 1,4
         OFFepsC2(i,-1,k) = -(OFFeps2(i,1,k))
         OFFepsC2(i,0,k) = (OFFeps2(i,0,k))
         OFFepsC2(i,1,k) = -(OFFeps2(i,-1,k))
       enddo
      enddo      
      
**     Main Calculation: M(mu,nu), mu,nu=0,3      
          
      if(.not.isnan(abs(COST)))then
         CALL WIGNERdvv(COST,smalld)   
         else
           smalld(0,0,0)=1
         endif        
      if(Tunitarize3) then  ! T matrix unitarization       
       
        CALL wwww_wigner_final(OFFeps(1,-1,1),OFFeps(1,-1,2),q3,q4,Aoff)   ! Voff(1,2) -> Von(3,4)
        
        ! On-shell to On-shell      
        q1On(0) =q3(0)
        q1On(1) =rZero
        q1On(2) =rZero
        q1On(3) =sqrt(q3T)*sign(1.0d0,Q1(3))

        q2On(0) =q4(0)
        q2On(1) =-q1On(1)
        q2On(2) =-q1On(2)
        q2On(3) =-q1On(3)
          
        call voffxx(q1On,WMASS,-1,Nq(1),ONONeps(1,-1,1))
        call voffxx(q2On,WMASS,-1,Nq(2),ONONeps(1,-1,2))

        CALL wwww_wigner_final(ONONeps(1,-1,1),ONONeps(1,-1,2),q3,q4,Aoo)    ! Von(3,4)  -> Von(3,4)   
        
        ! T-matriz Unitarization calculation
        CALL tmatrixvv(Aoff,Aoo, Ecm_VV,Sq1,Sq2,Sq3,Sq4,9,1, U)   

!         M(lambda1,lambda2,lambda3,lambda4, block1, block2), lambdai=-1,1  , block accounts for the different 2->2 processes with the same total charge                      
        do lambda1 = -1,1
         do lambda2 = -1,1
          do lambda3 = -1,1
           do lambda4 = -1,1
             M(lambda1,lambda2,lambda3,lambda4) = 16.0d0*rPi*(U(lambda1,lambda2,lambda3,lambda4,1,1,0)*smalld(0,lambda1-lambda2,lambda3-lambda4) 
     &        + 3.0d0*U(lambda1,lambda2,lambda3,lambda4,1,1,1)*smalld(1,lambda1-lambda2,lambda3-lambda4) 
     &        + 5.0d0*U(lambda1,lambda2,lambda3,lambda4,1,1,2)*smalld(2,lambda1-lambda2,lambda3-lambda4) )  
           enddo
          enddo
         enddo
        enddo
!!        
      else   ! non unitarization used             
        CALL wwww_wigner_final(OFFeps(1,-1,1),OFFeps(1,-1,2),q3,q4,Aoff)   !Voff(1,2) -> Von(3,4)  
        
        do lambda1 = -1,1
         do lambda2 = -1,1
          do lambda3 = -1,1
           do lambda4 = -1,1
             M(lambda1,lambda2,lambda3,lambda4) = 16.0d0*rPi*( Aoff(lambda1,lambda2,lambda3,lambda4,1,1,0)*smalld(0,lambda1-lambda2,lambda3-lambda4)
     &        + 3.0d0*Aoff(lambda1,lambda2,lambda3,lambda4,1,1,1)*smalld(1,lambda1-lambda2,lambda3-lambda4) 
     &        + 5.0d0*Aoff(lambda1,lambda2,lambda3,lambda4,1,1,2)*smalld(2,lambda1-lambda2,lambda3-lambda4) )
           enddo
          enddo
         enddo
        enddo
        endif
        
c     construction of the leptonic tensor from the individual parts 
          do lambda = -1,1
           u3(lambda) = dotcc(OFFepsp(1,lambda,3),Web(1))
           u4(lambda) = dotcc(OFFepsp(1,lambda,4),Wub(1))
          enddo     

        
      do lambda1 = -1,1
          do lambda2 = -1,1
            do lambda3 = -1,1
              do lambda = -1,1
                M1(lambda1,lambda2,lambda3) = M1(lambda1,lambda2,lambda3) +
     &                                        M(lambda1,lambda2,lambda3,lambda)*u4(lambda)      
              enddo            
            enddo
          enddo
      enddo
        
      do lambda1 = -1,1
          do lambda2 = -1,1
            do lambda = -1,1
              M2(lambda1,lambda2) = M2(lambda1,lambda2) +
     &                              M1(lambda1,lambda2,lambda)*u3(lambda)      
            enddo            
          enddo
      enddo
        
      D1 = rOne/dcmplx(Sq1-WMASS**2, WMASS*WWIDTH)
      D2 = rOne/dcmplx(Sq2-WMASS**2, WMASS*WWIDTH)        
        
      do nu = 0,3
         do lambda1 = -1,1
          do lambda = -1,1
               MM(lambda1,nu) = MM(lambda1,nu) + (-1)**(lambda+1)
     &                         *d2*OFFepsC2(nu+1,lambda,2)*M2(lambda1,lambda)
          enddo
         enddo
      enddo
         
      do mu = 0,3
         do nu = 0,3
           do lambda1 = -1,1
                wwt(mu,nu) = wwt(mu,nu) + (-1)**(lambda1+1)
     &                     *d1*OFFepsC2(mu+1,lambda1,1)*MM(lambda1,nu)
           enddo
         enddo
      enddo         
              
      END      
      
      
 

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE VOFFXX(p,vmass,nsv,nx,vx)
c This subroutine computes a vector wavefunction with off-shell considerations.
c (Original version on-shell momenta, vxxxxx.F as calculated by helas)
c
c       p(0:3)              : four-momentum of the vector boson (off-shell)
c       vmass               : mass of the vector boson
c       nsv = -1,1          : +1 for final, -1 for initial
c       nx  = Nq, Nk or NN  : normalization factor for the off-shell polvec
c
c       vx(6,3)    : vector wavefunction epsilon^mu(v,lambda), with lambda=-1,0,1 (polarization of the vector boson)   

      IMPLICIT NONE

      DOUBLE COMPLEX VX(6,-1:1)
      DOUBLE PRECISION P(0:3), VMASS, NX
      DOUBLE PRECISION PP2, PT2, PT, PP, PV, PV2
      DOUBLE PRECISION HEL, HEL0, PMIN, PPS, PZPT, EMP
      INTEGER NSV, NSVAHL, I

      DOUBLE PRECISION rZERO, rHALF, rONE, rTWO
      PARAMETER( rZERO = 0.0D0, rHALF = 0.5D0 )
      PARAMETER( rONE = 1.0D0, rTWO = 2.0D0 )
      DOUBLE COMPLEX cIMAG, cZERO, cONE
      PARAMETER( cIMAG=(0.0D0,1.0D0), cZERO=(0.0D0,0.0D0) )
      PARAMETER( cONE=(1.0D0,0.0D0) )
                
      pt2 = p(1)**2+p(2)**2+p(3)**2
      pt = dsqrt(pt2)
      pv2 = p(1)**2+p(2)**2
      pv = dsqrt(pv2)

      pmin = min(dabs(p(0)),dabs(pt))
      pps = min(pmin,dabs(pv))

      do i = -1,1
         vx(5,i) = dcmplx(p(0),p(3))*nsv
         vx(6,i) = dcmplx(p(1),p(2))*nsv
      enddo

      if( vmass.ne.rZero ) then
         do i = -1,1
            hel = dble(i)
            nsvahl = nsv*iabs(i)
            hel0 = rOne-dabs(hel)
            if( pmin.ne.rZero ) then
               if( pps.eq.rZero ) then
                  vx(1,i) = dcmplx(hel0*dabs(pt)*Nx)
                  vx(2,i) = dcmplx(-hel*dsqrt(rHalf)) 
                  vx(3,i) = dcmplx(rZero,nsvahl*sign(dsqrt(rHalf),p(3)))      
                  vx(4,i) = dcmplx(hel0*p(0)*(p(3)/pt)*Nx)                  
               else
                  emp = (p(0)/pmin)*Nx
                  vx(1,i) = dcmplx(hel0*pmin*Nx)
                  pzpt = p(3)/(pmin*pps)*dsqrt(rHalf)*hel
                  vx(2,i) = dcmplx(hel0*p(1)*emp-p(1)*pzpt,
     &                 -nsvahl*p(2)/pps*dsqrt(rHalf))
                  vx(3,i) = dcmplx( hel0*p(2)*emp-p(2)*pzpt,
     &                 nsvahl*p(1)/pps*dsqrt(rHalf))
                  vx(4,i) = dcmplx(hel0*p(3)*emp+hel*pps/pmin*dsqrt(rHalf))
               endif
            else
               vx(1,i) = cZero
               vx(2,i) = dcmplx(-hel*dsqrt(rHalf))
               vx(3,i) = dsqrt(rHalf)*nsvahl*cImag
               vx(4,i) = dcmplx(hel0)
            endif
         enddo
      else
         write(*,*)
     &     'Error: mv = 0 is not a valid value for voffxx subroutine.'
         write(*,*)
     &     'It has not been yet implemented.'
      endif

      RETURN
      END

      
