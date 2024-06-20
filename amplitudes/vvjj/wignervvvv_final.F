c This file contains the partial wave decomposition
c for the quartic-vector boson couplings (and the Wigner subroutines)

      SUBROUTINE vvvv_wigner_final(polva,polvw,q3,q4,A)
ccccccccccccccc
c Partial wave decomposition for VVVV-vertex
c
c This subroutine calculates the vvvv-vertex with anomalous couplings
c and the on-shell polarization vectors for vector bosons
c
c IN:
c     polva,polvw,      : vector bosons polarization vectors (CM-frame)
c     q3, q4            : momenta of the outgoing vector bosons (CM-frame)

c OUT:
c     A(lambda1,lambda2,lambda3,lambda4) : Wigner decomposition
c
      IMPLICIT NONE

      INTEGER I,J,K
      INTEGER LAMBDA, LAMBDA1, LAMBDA2, LAMBDA3, LAMBDA4 , BLOCK1, BLOCK2     

      INTEGER NFUNCS
      PARAMETER( NFUNCS = 10)      
      DOUBLE PRECISION rZERO, rONE, rPI, rHALF
      PARAMETER(rZERO = 0.0D0,rONE=1.0D0,rPI=3.14159265358979323846d0,rHALF= 1.0D0/2.0D0)
      DOUBLE COMPLEX cZERO, cONE
      PARAMETER(cZERO=(0.0D0,0.0D0), cONE=(1.0D0,0.0D0))
      
      DOUBLE PRECISION Q3(0:3), Q4(0:3), P3(0:3,NFUNCS), P4(0:3,NFUNCS),sameprod(1:11)
      DOUBLE PRECISION Q3T, Q4T, QQ3, QQ4
      DOUBLE PRECISION COSP(NFUNCS), SINP(NFUNCS), COST
      DOUBLE PRECISION NQ3, NQ4,NQW, NQA
      DOUBLE PRECISION Dinv(-2:2,-2:2,0:2,1:3)
      DOUBLE PRECISION smalld(0:2,-2:2,-2:2)      
      DOUBLE COMPLEX M(-1:1,-1:1,-1:1,-1:1,1:11,1:11,1:3)
      DOUBLE COMPLEX Msorted(-1:1,-1:1,-1:1,-1:1,1:11,1:11,1:3)
      DOUBLE COMPLEX A(-1:1,-1:1,-1:1,-1:1,1:11,1:11,0:2)
   
      DOUBLE COMPLEX POLVA(6,-1:1,4), POLVA3(6,-1:1,3),POLVA4(6,-1:1,3),POLVW(6,-1:1,4), POLVW3(6,-1:1,3), POLVW4(6,-1:1,3)


#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"

      p3(0:3,1:6) = rZero
      p4(0:3,1:6) = rZero

      POLVA3(6,-1:1,3)=cZERO
      POLVA4(6,-1:1,3)=cZERO

      POLVW3(6,-1:1,3)=cZERO
      POLVW4(6,-1:1,3)=cZERO


      a(-1:1,-1:1,-1:1,-1:1,1:11,1:11,0:2) = cZero

      M(-1:1,-1:1,-1:1,-1:1,1:11,1:11,1:3) = cZero
      sameprod(1:11)=rHALF
      sameprod(1)=rHALF
      sameprod(5)=rHALF
      sameprod(9)=rHALF
!!!!!!!!!!!!!!!!!!!!!!!!



      
      q3T = q3(1)**2+q3(2)**2+q3(3)**2
      q4T = q4(1)**2+q4(2)**2+q4(3)**2

      
      qq3 = q3(0)**2- q3T
      qq4 = q4(0)**2- q4T

!Normalization for longitudinal polarization, Photon with 1/ZMASS      
      Nqa = rOne/ZMASS
      Nqw = rONE/WMASS

      
!theta2 = pi/3  (60)      
      cosP(1) = rHALF
      sinP(1) = sqrt(3.0d0)*rHalf
!theta3 = pi/2  (90)
      cosP(2) = rZERO
      sinP(2) = rONE
! theta6 = 2pi/3 (120)     
      cosP(3) = -rHALF
      sinP(3) = sqrt(3.0d0)*rHalf


c     momenta for three set angles
      do i = 1,3
         p3(0,i) = q3(0)
         p3(1,i) = sqrt(q3T)*sinP(i)
         p3(2,i) = rZero               
         p3(3,i) = sqrt(q3T)*cosP(i)
         
         p4(0,i) = q4(0)
         p4(1,i) = - p3(1,i)
         p4(2,i) = rZero        
         p4(3,i) = - p3(3,i)         
      enddo

c     outgoing polarization vectors
      do i = 1,3
         call voffxx(p3(0,i),ZMASS,+1,Nqa,polva3(1,-1,i))  
         call voffxx(p3(0,i),WMASS,+1,Nqw,polvw3(1,-1,i))
         call voffxx(p4(0,i),ZMASS,+1,Nqa,polva4(1,-1,i))
         call voffxx(p4(0,i),WMASS,+1,Nqw,polvw4(1,-1,i))              
      enddo    
            
      
! M(lambda1,lambda2,lambda3,lambda4,block1,block2,i), lambdai=-1,1, 
! block1 & block2 represent the different 2 particle states, i is for the different angles
!    

      do i=1,3
        CALL aaaa_ONanom4(polva(1,-1,1),polva(1,-1,2),polva3(1,-1,i),polva4(1,-1,i),M(-1,-1,-1,-1,1,1,i))      !aaaa   
                                                                                                               !aaah   
        CALL aaaz_ONanom4(polva(1,-1,1),polva(1,-1,2),polva3(1,-1,i),polva4(1,-1,i),M(-1,-1,-1,-1,1,3,i))      !aaaz
                                                                                                               !aaha
        CALL aahh_ONanom4(polva(1,-1,1),polva(1,-1,2),polva3(1,-1,i),polva4(1,-1,i),M(-1,-1,-1,-1,1,5,i))      !aahh
                                                                                                               !aahz
        CALL aaaz_ONanom4(polva(1,-1,1),polva(1,-1,2),polva4(1,-1,i),polva3(1,-1,i),M(-1,-1,-1,-1,1,7,i))      !aaza
                                                                                                               !aazh
        CALL aazz_ONanom4(polva(1,-1,1),polva(1,-1,2),polva3(1,-1,i),polva4(1,-1,i),M(-1,-1,-1,-1,1,9,i))      !aazz
        CALL aawmwp_ONanom4(polva(1,-1,1),polva(1,-1,2),polvw3(1,-1,i),polvw4(1,-1,i),M(-1,-1,-1,-1,1,10,i))   !aawmwp
        CALL aawmwp_ONanom4(polva(1,-1,1),polva(1,-1,2),polvw4(1,-1,i),polvw3(1,-1,i),M(-1,-1,-1,-1,1,11,i))   !aawpwm  11
                                                                                                               !ahaa
        CALL aahh_ONanom4(polva(1,-1,1),polva3(1,-1,i),polva(1,-1,2),polva4(1,-1,i),M(-1,-1,-1,-1,2,2,i))      !ahah
                                                                                                               !ahaz
        CALL aahh_ONanom4(polva(1,-1,1),polva4(1,-1,i),polva(1,-1,2),polva3(1,-1,i),M(-1,-1,-1,-1,2,4,i))      !ahha
                                                                                                               !ahhh
        CALL ahhz_ONanom4(polva(1,-1,1),polva(1,-1,2),polva3(1,-1,i),polva4(1,-1,i),M(-1,-1,-1,-1,2,6,i))      !ahhz
                                                                                                               !ahza
        CALL ahhz_ONanom4(polva(1,-1,1),polva(1,-1,2),polva4(1,-1,i),polva3(1,-1,i),M(-1,-1,-1,-1,2,8,i))      !ahzh
        CALL ahzz_ONanom4(polva(1,-1,1),polva(1,-1,2),polva3(1,-1,i),polva4(1,-1,i),M(-1,-1,-1,-1,2,9,i))      !ahzz       
        CALL ahwmwp_ONanom4(polva(1,-1,1),polva(1,-1,2),polvw3(1,-1,i),polvw4(1,-1,i),M(-1,-1,-1,-1,2,10,i))   !ahwmwp      
        CALL ahwmwp_ONanom4(polva(1,-1,1),polva(1,-1,2),polvw4(1,-1,i),polvw3(1,-1,i),M(-1,-1,-1,-1,2,11,i))   !ahwpwm  22
        CALL aaaz_ONanom4(polva(1,-1,1),polva3(1,-1,i),polva4(1,-1,i),polva(1,-1,2),M(-1,-1,-1,-1,3,1,i))      !azaa
                                                                                                               !azah
        CALL aazz_ONanom4(polva(1,-1,1),polva3(1,-1,i),polva(1,-1,2),polva4(1,-1,i),M(-1,-1,-1,-1,3,3,i))      !azaz
                                                                                                               !azha
        CALL ahhz_ONanom4(polva(1,-1,1),polva3(1,-1,i),polva4(1,-1,i),polva(1,-1,2),M(-1,-1,-1,-1,3,5,i))      !azhh
        CALL ahzz_ONanom4(polva(1,-1,1),polva3(1,-1,i),polva(1,-1,2),polva4(1,-1,i),M(-1,-1,-1,-1,3,6,i))      !azhz
        CALL aazz_ONanom4(polva(1,-1,1),polva4(1,-1,i),polva(1,-1,2),polva3(1,-1,i),M(-1,-1,-1,-1,3,7,i))      !azza
        CALL ahzz_ONanom4(polva(1,-1,1),polva4(1,-1,i),polva(1,-1,2),polva3(1,-1,i),M(-1,-1,-1,-1,3,8,i))      !azzh
        CALL azzz_ONanom4(polva(1,-1,1),polva(1,-1,2),polva3(1,-1,i),polva4(1,-1,i),M(-1,-1,-1,-1,3,9,i))      !azzz  
        CALL azwmwp_ONanom4(polva(1,-1,1),polva(1,-1,2),polvw3(1,-1,i),polvw4(1,-1,i),M(-1,-1,-1,-1,3,10,i))   !azwmwp  
        CALL azwmwp_ONanom4(polva(1,-1,1),polva(1,-1,2),polvw4(1,-1,i),polvw3(1,-1,i),M(-1,-1,-1,-1,3,11,i))   !azwpwm 33
                                                                                                               !haaa
        CALL aahh_ONanom4(polva(1,-1,2),polva3(1,-1,i),polva(1,-1,1),polva4(1,-1,i),M(-1,-1,-1,-1,4,2,i))      !haah
                                                                                                               !haaz
        CALL aahh_ONanom4(polva(1,-1,2),polva4(1,-1,i),polva(1,-1,1),polva3(1,-1,i),M(-1,-1,-1,-1,4,4,i))      !haha
                                                                                                               !hahh
        CALL ahhz_ONanom4(polva(1,-1,2),polva(1,-1,1),polva3(1,-1,i),polva4(1,-1,i),M(-1,-1,-1,-1,4,6,i))      !hahz
                                                                                                               !haza
        CALL ahhz_ONanom4(polva(1,-1,2),polva(1,-1,1),polva4(1,-1,i),polva3(1,-1,i),M(-1,-1,-1,-1,4,8,i))      !hazh
        CALL ahzz_ONanom4(polva(1,-1,2),polva(1,-1,1),polva3(1,-1,i),polva4(1,-1,i),M(-1,-1,-1,-1,4,9,i))      !hazz       
        CALL ahwmwp_ONanom4(polva(1,-1,2),polva(1,-1,1),polvw3(1,-1,i),polvw4(1,-1,i),M(-1,-1,-1,-1,4,10,i))   !hawmwp      
        CALL ahwmwp_ONanom4(polva(1,-1,2),polva(1,-1,1),polvw4(1,-1,i),polvw3(1,-1,i),M(-1,-1,-1,-1,4,11,i))   !hawpwm 44
        CALL aahh_ONanom4(polva3(1,-1,i),polva4(1,-1,i),polva(1,-1,1),polva(1,-1,2),M(-1,-1,-1,-1,5,1,i))      !hhaa
                                                                                                               !hhah
        CALL ahhz_ONanom4(polva3(1,-1,i),polva(1,-1,1),polva(1,-1,2),polva4(1,-1,i),M(-1,-1,-1,-1,5,3,i))      !hhaz
                                                                                                               !hhha
        CALL hhhh_ONanom4(polva(1,-1,1),polva(1,-1,2),polva3(1,-1,i),polva4(1,-1,i),M(-1,-1,-1,-1,5,5,i))      !hhhh
                                                                                                               !hhhz
        CALL ahhz_ONanom4(polva4(1,-1,i),polva(1,-1,1),polva(1,-1,2),polva3(1,-1,i),M(-1,-1,-1,-1,5,7,i))      !hhza
                                                                                                               !hhzh
        CALL hhzz_ONanom4(polva(1,-1,1),polva(1,-1,2),polva3(1,-1,i),polva4(1,-1,i),M(-1,-1,-1,-1,5,9,i))      !hhzz  
        CALL hhwmwp_ONanom4(polva(1,-1,1),polva(1,-1,2),polvw3(1,-1,i),polvw4(1,-1,i),M(-1,-1,-1,-1,5,10,i))   !hhwmwp 
        CALL hhwmwp_ONanom4(polva(1,-1,1),polva(1,-1,2),polvw4(1,-1,i),polvw3(1,-1,i),M(-1,-1,-1,-1,5,11,i))   !hhwpwm
                                                                                                               !hzaa
        CALL ahhz_ONanom4(polva3(1,-1,i),polva(1,-1,1),polva4(1,-1,i),polva(1,-1,2),M(-1,-1,-1,-1,6,2,i))      !hzah
        CALL ahzz_ONanom4(polva3(1,-1,i),polva(1,-1,1),polva(1,-1,2),polva4(1,-1,i),M(-1,-1,-1,-1,6,3,i))      !hzaz
        CALL ahhz_ONanom4(polva4(1,-1,i),polva(1,-1,1),polva3(1,-1,i),polva(1,-1,2),M(-1,-1,-1,-1,6,4,i))      !hzha 
                                                                                                               !hzhh
        CALL hhzz_ONanom4(polva(1,-1,1),polva3(1,-1,i),polva(1,-1,2),polva4(1,-1,i),M(-1,-1,-1,-1,6,6,i))      !hzhz 
        CALL ahzz_ONanom4(polva4(1,-1,i),polva(1,-1,1),polva(1,-1,2),polva3(1,-1,i),M(-1,-1,-1,-1,6,7,i))      !hzza
        CALL hhzz_ONanom4(polva(1,-1,1),polva4(1,-1,i),polva(1,-1,2),polva3(1,-1,i),M(-1,-1,-1,-1,6,8,i))      !hzzh 
                                                                                                               !hzzz
        CALL hzwmwp_ONanom4(polva(1,-1,1),polva(1,-1,2),polvw3(1,-1,i),polvw4(1,-1,i),M(-1,-1,-1,-1,6,10,i))   !hzwmwp  
        CALL hzwmwp_ONanom4(polva(1,-1,1),polva(1,-1,2),polvw4(1,-1,i),polvw3(1,-1,i),M(-1,-1,-1,-1,6,11,i))   !hzwpwm
        CALL aaaz_ONanom4(polva(1,-1,2),polva3(1,-1,i),polva4(1,-1,i),polva(1,-1,1),M(-1,-1,-1,-1,7,1,i))      !zaaa
                                                                                                               !zaah
        CALL aazz_ONanom4(polva(1,-1,2),polva3(1,-1,i),polva(1,-1,1),polva4(1,-1,i),M(-1,-1,-1,-1,7,3,i))      !zaaz
                                                                                                               !zaha
        CALL ahhz_ONanom4(polva(1,-1,2),polva3(1,-1,i),polva4(1,-1,i),polva(1,-1,1),M(-1,-1,-1,-1,7,5,i))      !zahh
        CALL ahzz_ONanom4(polva(1,-1,2),polva3(1,-1,i),polva(1,-1,1),polva4(1,-1,i),M(-1,-1,-1,-1,7,6,i))      !zahz
        CALL aazz_ONanom4(polva(1,-1,2),polva4(1,-1,i),polva(1,-1,1),polva3(1,-1,i),M(-1,-1,-1,-1,7,7,i))      !zaza
        CALL ahzz_ONanom4(polva(1,-1,2),polva4(1,-1,i),polva(1,-1,1),polva3(1,-1,i),M(-1,-1,-1,-1,7,8,i))      !zazh
        CALL azzz_ONanom4(polva(1,-1,2),polva(1,-1,1),polva3(1,-1,i),polva4(1,-1,i),M(-1,-1,-1,-1,7,9,i))      !zazz  
        CALL azwmwp_ONanom4(polva(1,-1,2),polva(1,-1,1),polvw3(1,-1,i),polvw4(1,-1,i),M(-1,-1,-1,-1,7,10,i))   !zawmwp  
        CALL azwmwp_ONanom4(polva(1,-1,2),polva(1,-1,1),polvw4(1,-1,i),polvw3(1,-1,i),M(-1,-1,-1,-1,7,11,i))   !zawpwm 
                                                                                                               !zhaa
        CALL ahhz_ONanom4(polva3(1,-1,i),polva(1,-1,2),polva4(1,-1,i),polva(1,-1,1),M(-1,-1,-1,-1,8,2,i))      !zhah
        CALL ahzz_ONanom4(polva3(1,-1,i),polva(1,-1,2),polva(1,-1,1),polva4(1,-1,i),M(-1,-1,-1,-1,8,3,i))      !zhaz
        CALL ahhz_ONanom4(polva4(1,-1,i),polva(1,-1,2),polva3(1,-1,i),polva(1,-1,1),M(-1,-1,-1,-1,8,4,i))      !zhha 
                                                                                                               !zhhh
        CALL hhzz_ONanom4(polva(1,-1,2),polva3(1,-1,i),polva(1,-1,1),polva4(1,-1,i),M(-1,-1,-1,-1,8,6,i))      !zhhz 
        CALL ahzz_ONanom4(polva4(1,-1,i),polva(1,-1,2),polva(1,-1,1),polva3(1,-1,i),M(-1,-1,-1,-1,8,7,i))      !zhza
        CALL hhzz_ONanom4(polva(1,-1,2),polva4(1,-1,i),polva(1,-1,1),polva3(1,-1,i),M(-1,-1,-1,-1,8,8,i))      !zhzh 
                                                                                                               !zhzz
        CALL hzwmwp_ONanom4(polva(1,-1,2),polva(1,-1,1),polvw3(1,-1,i),polvw4(1,-1,i),M(-1,-1,-1,-1,8,10,i))   !zhwmwp  
        CALL hzwmwp_ONanom4(polva(1,-1,2),polva(1,-1,1),polvw4(1,-1,i),polvw3(1,-1,i),M(-1,-1,-1,-1,8,11,i))   !zhwpwm 
        CALL aazz_ONanom4(polva3(1,-1,i),polva4(1,-1,i),polva(1,-1,1),polva(1,-1,2),M(-1,-1,-1,-1,9,1,i))      !zzaa
        CALL ahzz_ONanom4(polva3(1,-1,i),polva4(1,-1,i),polva(1,-1,1),polva(1,-1,2),M(-1,-1,-1,-1,9,2,i))      !zzah
        CALL azzz_ONanom4(polva3(1,-1,i),polva(1,-1,1),polva(1,-1,2),polva4(1,-1,i),M(-1,-1,-1,-1,9,3,i))      !zzaz
        CALL ahzz_ONanom4(polva3(1,-1,i),polva4(1,-1,i),polva(1,-1,1),polva(1,-1,2),M(-1,-1,-1,-1,9,4,i))      !zzha 
        CALL hhzz_ONanom4(polva3(1,-1,i),polva4(1,-1,i),polva(1,-1,1),polva(1,-1,2),M(-1,-1,-1,-1,9,5,i))      !zzhh  
                                                                                                               !zzhz
        CALL azzz_ONanom4(polva4(1,-1,i),polva(1,-1,1),polva(1,-1,2),polva3(1,-1,i),M(-1,-1,-1,-1,9,7,i))      !zzza 
                                                                                                               !zzzh
        CALL zzzz_ONanom4(polva(1,-1,1),polva(1,-1,2),polva3(1,-1,i),polva4(1,-1,i),M(-1,-1,-1,-1,9,9,i))      !zzzz  
        CALL zzwmwp_ONanom4(polva(1,-1,1),polva(1,-1,2),polvw3(1,-1,i),polvw4(1,-1,i),M(-1,-1,-1,-1,9,10,i))   !zzwmwp
        CALL zzwmwp_ONanom4(polva(1,-1,1),polva(1,-1,2),polvw4(1,-1,i),polvw3(1,-1,i),M(-1,-1,-1,-1,9,11,i))   !zzwpwm
        CALL aawmwp_ONanom4(polva3(1,-1,i),polva4(1,-1,i),polvw(1,-1,1),polvw(1,-1,2),M(-1,-1,-1,-1,10,1,i))   !wmwpaa
        CALL aawmwp_ONanom4(polva3(1,-1,i),polva4(1,-1,i),polvw(1,-1,2),polvw(1,-1,1),M(-1,-1,-1,-1,11,1,i))   !wpwmaa     
        CALL ahwmwp_ONanom4(polva3(1,-1,i),polva4(1,-1,i),polvw(1,-1,1),polvw(1,-1,2),M(-1,-1,-1,-1,10,2,i))   !wmwpah      
        CALL ahwmwp_ONanom4(polva3(1,-1,i),polva4(1,-1,i),polvw(1,-1,2),polvw(1,-1,1),M(-1,-1,-1,-1,11,2,i))   !wpwmah
        CALL azwmwp_ONanom4(polva3(1,-1,i),polva4(1,-1,i),polvw(1,-1,1),polvw(1,-1,2),M(-1,-1,-1,-1,10,3,i))   !wmwpaz  
        CALL azwmwp_ONanom4(polva3(1,-1,i),polva4(1,-1,i),polvw(1,-1,2),polvw(1,-1,1),M(-1,-1,-1,-1,11,3,i))   !wpwmaz    
        CALL ahwmwp_ONanom4(polva4(1,-1,i),polva3(1,-1,i),polvw(1,-1,1),polvw(1,-1,2),M(-1,-1,-1,-1,10,4,i))   !wmwpha      
        CALL ahwmwp_ONanom4(polva4(1,-1,i),polva3(1,-1,i),polvw(1,-1,2),polvw(1,-1,1),M(-1,-1,-1,-1,11,4,i))   !wpwmha
        CALL hhwmwp_ONanom4(polva3(1,-1,i),polva4(1,-1,i),polvw(1,-1,1),polvw(1,-1,2),M(-1,-1,-1,-1,10,5,i))   !wmwphh 
        CALL hhwmwp_ONanom4(polva3(1,-1,i),polva4(1,-1,i),polvw(1,-1,2),polvw(1,-1,1),M(-1,-1,-1,-1,11,5,i))   !wpwmhh
        CALL hzwmwp_ONanom4(polva3(1,-1,i),polva4(1,-1,i),polvw(1,-1,1),polvw(1,-1,2),M(-1,-1,-1,-1,10,6,i))   !wmwphz  
        CALL hzwmwp_ONanom4(polva3(1,-1,i),polva4(1,-1,i),polvw(1,-1,2),polvw(1,-1,1),M(-1,-1,-1,-1,11,6,i))   !wpwmhz
        CALL azwmwp_ONanom4(polva4(1,-1,i),polva3(1,-1,i),polvw(1,-1,1),polvw(1,-1,2),M(-1,-1,-1,-1,10,7,i))   !wmwpza  
        CALL azwmwp_ONanom4(polva4(1,-1,i),polva3(1,-1,i),polvw(1,-1,2),polvw(1,-1,1),M(-1,-1,-1,-1,11,7,i))   !wpwmza 
        CALL hzwmwp_ONanom4(polva4(1,-1,i),polva3(1,-1,i),polvw(1,-1,1),polvw(1,-1,2),M(-1,-1,-1,-1,10,8,i))   !wmwpzh  
        CALL hzwmwp_ONanom4(polva4(1,-1,i),polva3(1,-1,i),polvw(1,-1,2),polvw(1,-1,1),M(-1,-1,-1,-1,11,8,i))   !wpwmzh 
        CALL zzwmwp_ONanom4(polva3(1,-1,i),polva4(1,-1,i),polvw(1,-1,1),polvw(1,-1,2),M(-1,-1,-1,-1,10,9,i))   !wmwpzz
        CALL zzwmwp_ONanom4(polva3(1,-1,i),polva4(1,-1,i),polvw(1,-1,2),polvw(1,-1,1),M(-1,-1,-1,-1,11,9,i))   !wpwmzz
        CALL wpwmwmwp_ONanom4(polvw(1,-1,1),polvw(1,-1,2) ,polvw3(1,-1,i),polvw4(1,-1,i),M(-1,-1,-1,-1,11,10,i))!wpwmwmwp
        CALL wpwmwmwp_ONanom4(polvw(1,-1,2),polvw(1,-1,1) ,polvw3(1,-1,i),polvw4(1,-1,i),M(-1,-1,-1,-1,10,10,i))!wmwpwmwp
        CALL wpwmwmwp_ONanom4(polvw(1,-1,1),polvw(1,-1,2) ,polvw4(1,-1,i),polvw3(1,-1,i),M(-1,-1,-1,-1,11,11,i))!wpwmwpwm
        CALL wpwmwmwp_ONanom4(polvw(1,-1,2),polvw(1,-1,1) ,polvw4(1,-1,i),polvw3(1,-1,i),M(-1,-1,-1,-1,10,11,i))!wmwpwpwm
        
      enddo
!      from abcd_ONanom4: M(lambda(a),lambda(b),lambda(c),lambda(d)
c      swapping helicities so that Msorted= M(lambda(particle1),lambda(particle2),lambda(particle3),lambda(particle4)    
       call sortmatrixpol(0, 11, M(-1,-1,-1,-1,1,1,1), Msorted(-1,-1,-1,-1,1,1,1))

      CALL inv_Wignervv(Dinv)   



c     A (i) = Dinv(lambda1-lambda2, lambda3 - lambda4, i , j)*M(j)  j=theta(j), sameprod accounts for same particles in final sates as well as xy and yx being the same state after integration over theta                  
      do lambda1 = -1,1
       do lambda2 = -1,1
        do lambda3 = -1,1
         do lambda4 = -1,1
          do block1 = 1,11
           do block2 = 1,11
            do i= 0,2
             a(lambda1,lambda2,lambda3,lambda4,block1,block2,i) = sameprod(block2)*((Dinv(lambda1-lambda2,lambda3-lambda4,i,1))*Msorted(lambda1,lambda2,lambda3,lambda4,block1,block2,1) 
     &                                                            +(Dinv(lambda1-lambda2,lambda3-lambda4,i,2))*Msorted(lambda1,lambda2,lambda3,lambda4,block1,block2,2)
     &                                                            +(Dinv(lambda1-lambda2,lambda3-lambda4,i,3))*Msorted(lambda1,lambda2,lambda3,lambda4,block1,block2,3))

            enddo
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo       

      END








      SUBROUTINE wpvwpv_wigner_final(polva,polvw,q3,q4,A)
ccccccccccccccc
c Partial wave decomposition for WvWv-vertex
c
c This subroutine calculates the wwww-vertex with anomalous couplings
c and the on-shell polarization vectors for vector bosons
c
c IN:
c     polVa,polVw      : incoming vector bosons polarization vectors (CM-frame)
c     q3,q4            : outgoing vector bosons momenta in the CM-frame
c OUT:
c     A(lambda1,lambda2,lambda3,lambda4) : Wigner decomposition
c
      IMPLICIT NONE

      INTEGER I,J,K
      INTEGER LAMBDA, LAMBDA1, LAMBDA2, LAMBDA3, LAMBDA4 , BLOCK1, BLOCK2     

      INTEGER NFUNCS
      PARAMETER( NFUNCS = 10)      
      DOUBLE PRECISION rZERO, rONE, rPI, rHALF
      PARAMETER(rZERO = 0.0D0, rONE=1.0D0, rPI=3.14159265358979323846d0, rHALF=1.0D0/2.0D0)
      DOUBLE COMPLEX cZERO, cONE
      PARAMETER(cZERO=(0.0D0,0.0D0), cONE=(1.0D0,0.0D0))
      
      DOUBLE PRECISION Q3(0:3), Q4(0:3), P3(0:3,NFUNCS), P4(0:3,NFUNCS), sameprod(6)
      DOUBLE PRECISION Q3T, Q4T, QQ3, QQ4
      DOUBLE PRECISION COSP(NFUNCS), SINP(NFUNCS), COST
      DOUBLE PRECISION NQ3, NQ4, Nqa, Nqw
      DOUBLE PRECISION Dinv(-2:2,-2:2,0:2,1:3)
      DOUBLE PRECISION smalld(0:2,-2:2,-2:2)      
      DOUBLE COMPLEX M(-1:1,-1:1,-1:1,-1:1,1:3,1:3,1:3)
      DOUBLE COMPLEX Msorted(-1:1,-1:1,-1:1,-1:1,1:3,1:3,1:3)      
      DOUBLE COMPLEX A(-1:1,-1:1,-1:1,-1:1,1:3,1:3,0:2)  
      DOUBLE COMPLEX POLVA(6,-1:1,4), POLVA3(6,-1:1,3),POLVA4(6,-1:1,3),POLVW(6,-1:1,4), POLVW3(6,-1:1,3), POLVW4(6,-1:1,3)




#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"
      POLVA3(6,-1:1,3)=cZERO
      POLVA4(6,-1:1,3)=cZERO

      POLVW3(6,-1:1,3)=cZERO
      POLVW4(6,-1:1,3)=cZERO


      a(-1:1,-1:1,-1:1,-1:1,1:3,1:3,0:2) = cZero
      M(-1:1,-1:1,-1:1,-1:1,1:3,1:3,1:3) = cZero
      sameprod(1:3)=rONE


!!!!!!!!!!!!!!!!!!!!!!!!
! Structure identical to vvvv_wigner_final, look there for explanation
      q3T = q3(1)**2+q3(2)**2+q3(3)**2
      q4T = q4(1)**2+q4(2)**2+q4(3)**2

      
      qq3 = q3(0)**2- q3T
      qq4 = q4(0)**2- q4T

      
      Nqa = rOne/ZMASS
      Nqw = rONE/WMASS
      
!theta2 = pi/3  (60)      
      cosP(1) = rHALF
      sinP(1) = sqrt(3.0d0)*rHalf
!theta3 = pi/2  (90)
      cosP(2) = rZERO
      sinP(2) = rONE
! theta6 = 2pi/3 (120)     
      cosP(3) = -rHALF
      sinP(3) = sqrt(3.0d0)*rHalf



      do i = 1,3
         p3(0,i) = q3(0)
         p3(1,i) = sqrt(q3T)*sinP(i)
         p3(2,i) = rZero               
         p3(3,i) = sqrt(q3T)*cosP(i)
         
         p4(0,i) = q4(0)
         p4(1,i) = - p3(1,i)
         p4(2,i) = rZero        
         p4(3,i) = - p3(3,i)         
      enddo        
      
      do i = 1,3
         call voffxx(p3(0,i),ZMASS,+1,Nqa,polva3(1,-1,i))  
         call voffxx(p3(0,i),WMASS,+1,Nqw,polvw3(1,-1,i))
         call voffxx(p4(0,i),ZMASS,+1,Nqa,polva4(1,-1,i))
         call voffxx(p4(0,i),WMASS,+1,Nqw,polvw4(1,-1,i))              
      enddo              
      
      do i=1,3

        CALL aawmwp_ONanom4(polva(1,-1,2),polva4(1,-1,i),polvw3(1,-1,i),polvw(1,-1,1),M(-1,-1,-1,-1,1,1,i))   !wpawpa
        CALL ahwmwp_ONanom4(polva(1,-1,2),polva4(1,-1,i),polvw3(1,-1,i),polvw(1,-1,1),M(-1,-1,-1,-1,1,2,i))   !wpawph
        CALL azwmwp_ONanom4(polva(1,-1,2),polva4(1,-1,i),polvw3(1,-1,i),polvw(1,-1,1),M(-1,-1,-1,-1,1,3,i))   !wpawpz

        CALL ahwmwp_ONanom4(polva4(1,-1,i),polva(1,-1,2),polvw3(1,-1,i),polvw(1,-1,1),M(-1,-1,-1,-1,2,1,i))   !wphwpa
        CALL hhwmwp_ONanom4(polva(1,-1,2),polva4(1,-1,i),polvw3(1,-1,i),polvw(1,-1,1),M(-1,-1,-1,-1,2,2,i))   !wphwph
        CALL hzwmwp_ONanom4(polva(1,-1,2),polva4(1,-1,i),polvw3(1,-1,i),polvw(1,-1,1),M(-1,-1,-1,-1,2,3,i))   !wphwpz

        CALL azwmwp_ONanom4(polva4(1,-1,i),polva(1,-1,2),polvw3(1,-1,i),polvw(1,-1,1),M(-1,-1,-1,-1,3,1,i))   !wpzwpa
        CALL hzwmwp_ONanom4(polva4(1,-1,i),polva(1,-1,2),polvw3(1,-1,i),polvw(1,-1,1),M(-1,-1,-1,-1,3,2,i))   !wpzwph
        CALL zzwmwp_ONanom4(polva(1,-1,2),polva4(1,-1,i),polvw3(1,-1,i),polvw(1,-1,1),M(-1,-1,-1,-1,3,3,i))   !wpzwpz


      enddo


       call sortmatrixpol(1, 3, M(-1,-1,-1,-1,1,1,1), Msorted(-1,-1,-1,-1,1,1,1))

 
      CALL inv_Wignervv(Dinv)   

      
                 
      do lambda1 = -1,1
       do lambda2 = -1,1
        do lambda3 = -1,1
         do lambda4 = -1,1
          do block1 = 1,3
           do block2 = 1,3
            do i= 0,2
             a(lambda1,lambda2,lambda3, lambda4,block1,block2,i) =sameprod(block2)*((Dinv(lambda1-lambda2,lambda3-lambda4,i,1))*Msorted(lambda1,lambda2,lambda3,lambda4,block1,block2,1) 
     &                                                            +(Dinv(lambda1-lambda2,lambda3-lambda4,i,2))*Msorted(lambda1,lambda2,lambda3,lambda4,block1,block2,2)
     &                                                            +(Dinv(lambda1-lambda2,lambda3-lambda4,i,3))*Msorted(lambda1,lambda2,lambda3,lambda4,block1,block2,3))

            enddo
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo       


    
      END



      SUBROUTINE wmvwmv_wigner_final(polva,polvw,q3,q4,A)
ccccccccccccccc
c Partial wave decomposition for WvWv-vertex
c
c This subroutine calculates the wwww-vertex with anomalous couplings
c and the on-shell polarization vectors for vector bosons
c
c IN:
c     polVa,polVw      : incoming vector bosons polarization vectors (CM-frame)
c     q3,q4            : outgoing vector bosons momenta in the CM-frame
c OUT:
c     A(lambda1,lambda2,lambda3,lambda4) : Wigner decomposition
c
      IMPLICIT NONE

      INTEGER I,J,K
      INTEGER LAMBDA, LAMBDA1, LAMBDA2, LAMBDA3, LAMBDA4 , BLOCK1, BLOCK2     

      INTEGER NFUNCS
      PARAMETER( NFUNCS = 10)      
      DOUBLE PRECISION rZERO, rONE, rPI, rHALF
      PARAMETER(rZERO = 0.0D0, rONE=1.0D0, rPI=3.14159265358979323846d0, rHALF=1.0D0/2.0D0)
      DOUBLE COMPLEX cZERO, cONE
      PARAMETER(cZERO=(0.0D0,0.0D0), cONE=(1.0D0,0.0D0))
      
      DOUBLE PRECISION Q3(0:3), Q4(0:3), P3(0:3,NFUNCS), P4(0:3,NFUNCS), sameprod(6)
      DOUBLE PRECISION Q3T, Q4T, QQ3, QQ4
      DOUBLE PRECISION COSP(NFUNCS), SINP(NFUNCS), COST
      DOUBLE PRECISION NQ3, NQ4, Nqa, Nqw
      DOUBLE PRECISION Dinv(-2:2,-2:2,0:2,1:3)
      DOUBLE PRECISION smalld(0:2,-2:2,-2:2)      
      DOUBLE COMPLEX M(-1:1,-1:1,-1:1,-1:1,1:3,1:3,1:3)
      DOUBLE COMPLEX Msorted(-1:1,-1:1,-1:1,-1:1,1:3,1:3,1:3)      
      DOUBLE COMPLEX A(-1:1,-1:1,-1:1,-1:1,1:3,1:3,0:2)    
      DOUBLE COMPLEX POLVA(6,-1:1,4), POLVA3(6,-1:1,3),POLVA4(6,-1:1,3),POLVW(6,-1:1,4), POLVW3(6,-1:1,3), POLVW4(6,-1:1,3)



#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"
      POLVA3(6,-1:1,3)=cZERO
      POLVA4(6,-1:1,3)=cZERO

      POLVW3(6,-1:1,3)=cZERO
      POLVW4(6,-1:1,3)=cZERO


      a(-1:1,-1:1,-1:1,-1:1,1:3,1:3,0:2) = cZero
      M(-1:1,-1:1,-1:1,-1:1,1:3,1:3,1:3) = cZero
      sameprod(1:3)=rONE



!!!!!!!!!!!!!!!!!!!!!!!!
! Structure identical to vvvv_wigner_final, look there for explanation

      q3T = q3(1)**2+q3(2)**2+q3(3)**2
      q4T = q4(1)**2+q4(2)**2+q4(3)**2

      
      qq3 = q3(0)**2- q3T
      qq4 = q4(0)**2- q4T

      
      Nqa = rOne/ZMASS
      Nqw = rONE/WMASS
      
!theta2 = pi/3  (60)      
      cosP(1) = rHALF
      sinP(1) = sqrt(3.0d0)*rHalf
!theta3 = pi/2  (90)
      cosP(2) = rZERO
      sinP(2) = rONE
! theta6 = 2pi/3 (120)     
      cosP(3) = -rHALF
      sinP(3) = sqrt(3.0d0)*rHalf



      do i = 1,3
         p3(0,i) = q3(0)
         p3(1,i) = sqrt(q3T)*sinP(i)
         p3(2,i) = rZero               
         p3(3,i) = sqrt(q3T)*cosP(i)
         
         p4(0,i) = q4(0)
         p4(1,i) = - p3(1,i)
         p4(2,i) = rZero        
         p4(3,i) = - p3(3,i)         
      enddo        
      
      do i = 1,3
         call voffxx(p3(0,i),ZMASS,+1,Nqa,polva3(1,-1,i))  
         call voffxx(p3(0,i),WMASS,+1,Nqw,polvw3(1,-1,i))
         call voffxx(p4(0,i),ZMASS,+1,Nqa,polva4(1,-1,i))
         call voffxx(p4(0,i),WMASS,+1,Nqw,polvw4(1,-1,i))              
      enddo              
      

      do i=1,3
        CALL aawmwp_ONanom4(polva(1,-1,2),polva4(1,-1,i),polvw(1,-1,1),polvw3(1,-1,i),M(-1,-1,-1,-1,1,1,i))   !wmawma
        CALL ahwmwp_ONanom4(polva(1,-1,2),polva4(1,-1,i),polvw(1,-1,1),polvw3(1,-1,i),M(-1,-1,-1,-1,1,2,i))   !wmawmh
        CALL azwmwp_ONanom4(polva(1,-1,2),polva4(1,-1,i),polvw(1,-1,1),polvw3(1,-1,i),M(-1,-1,-1,-1,1,3,i))   !wmawmz

        CALL ahwmwp_ONanom4(polva4(1,-1,i),polva(1,-1,2),polvw(1,-1,1),polvw3(1,-1,i),M(-1,-1,-1,-1,2,1,i))   !wmhwma
        CALL hhwmwp_ONanom4(polva(1,-1,2),polva4(1,-1,i),polvw(1,-1,1),polvw3(1,-1,i),M(-1,-1,-1,-1,2,2,i))   !wmhwmh
        CALL hzwmwp_ONanom4(polva(1,-1,2),polva4(1,-1,i),polvw(1,-1,1),polvw3(1,-1,i),M(-1,-1,-1,-1,2,3,i))   !wmhwmz


        CALL azwmwp_ONanom4(polva4(1,-1,i),polva(1,-1,2),polvw(1,-1,1),polvw3(1,-1,i),M(-1,-1,-1,-1,3,1,i))   !wmzwma
        CALL hzwmwp_ONanom4(polva4(1,-1,i),polva(1,-1,2),polvw(1,-1,1),polvw3(1,-1,i),M(-1,-1,-1,-1,3,2,i))   !wmzwmh
        CALL zzwmwp_ONanom4(polva(1,-1,2),polva4(1,-1,i),polvw(1,-1,1),polvw3(1,-1,i),M(-1,-1,-1,-1,3,3,i))   !wmzwmz


      enddo
      

       call sortmatrixpol(-1, 3, M(-1,-1,-1,-1,1,1,1), Msorted(-1,-1,-1,-1,1,1,1))

      CALL inv_Wignervv(Dinv)   

     
  
                     
      do lambda1 = -1,1
       do lambda2 = -1,1
        do lambda3 = -1,1
         do lambda4 = -1,1
          do block1 = 1,3
           do block2 = 1,3
            do i= 0,2
             a(lambda1,lambda2,lambda3, lambda4,block1,block2,i) = sameprod(block2)*((Dinv(lambda1-lambda2,lambda3-lambda4,i,1))*Msorted(lambda1,lambda2,lambda3,lambda4,block1,block2,1) 
     &                                                            +(Dinv(lambda1-lambda2,lambda3-lambda4,i,2))*Msorted(lambda1,lambda2,lambda3,lambda4,block1,block2,2)
     &                                                            +(Dinv(lambda1-lambda2,lambda3-lambda4,i,3))*Msorted(lambda1,lambda2,lambda3,lambda4,block1,block2,3))

            enddo
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo       


    
      END      
      


      SUBROUTINE wwww_wigner_final(polv1,polv2,q3,q4,A)
ccccccccccccccc
c Partial wave decomposition for WWWW-vertex
c
c This subroutine calculates the wwww-vertex with anomalous couplings
c and the on-shell polarization vectors for vector bosons
c
c IN:
c     polV1,polV2      : incoming vector bosons polarization vectors (CM-frame)
c     k3,k4            : outgoing vector bosons momenta in the CM-frame
c OUT:
c     A(lambda1,lambda2) : Wigner decomposition
c
      IMPLICIT NONE

      INTEGER I,J,K
      INTEGER LAMBDA, LAMBDA1, LAMBDA2, LAMBDA3, LAMBDA4      

      INTEGER NFUNCS
      PARAMETER( NFUNCS = 10)      
      DOUBLE PRECISION rZERO, rONE, rPI, rHALF
      PARAMETER(rZERO = 0.0D0, rONE=1.0D0, rPI=3.14159265358979323846d0, rHALF=1.0D0/2.0D0)
      DOUBLE COMPLEX cZERO, cONE
      PARAMETER(cZERO=(0.0D0,0.0D0), cONE=(1.0D0,0.0D0))
      
      DOUBLE PRECISION Q3(0:3), Q4(0:3), P3(0:3,NFUNCS), P4(0:3,NFUNCS)
      DOUBLE PRECISION Q3T, Q4T, QQ3, QQ4
      DOUBLE PRECISION COSP(NFUNCS), SINP(NFUNCS), COST
      DOUBLE PRECISION NQ3, NQ4
      DOUBLE PRECISION Dinv(-2:2,-2:2,0:2,1:3)
      
      DOUBLE COMPLEX M(-1:1,-1:1,-1:1,-1:1,1:1,1:1,1:3)
      DOUBLE COMPLEX Msorted(-1:1,-1:1,-1:1,-1:1,1:1,1:1,1:3)
      DOUBLE COMPLEX A(-1:1,-1:1,-1:1,-1:1,1:1,1:1,0:2)
      DOUBLE COMPLEX POLV1(6,-1:1), POLV2(6,-1:1)
      DOUBLE COMPLEX POLV3(6,-1:1,NFUNCS), POLV4(6,-1:1,NFUNCS)


#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"

      p3(0:3,1:6) = rZero
      p4(0:3,1:6) = rZero
      polv3(1:6,-1:1,1:6) = cZero
      polv4(1:6,-1:1,1:6) = cZero
      a(-1:1,-1:1,-1:1,-1:1,1:1,1:1,0:2) = cZero
      M(-1:1,-1:1,-1:1,-1:1,1:1,1:1,1:3) = cZero
      

!!!!!!!!!!!!!!!!!!!!!!!!
! Structure identical to vvvv_wigner_final, look there for explanation
      q3T = q3(1)**2+q3(2)**2+q3(3)**2
      q4T = q4(1)**2+q4(2)**2+q4(3)**2
      
      qq3 = q3(0)**2- q3T
      qq4 = q4(0)**2- q4T
      

      Nq3 = rONE/WMASS
      Nq4 = rONE/WMASS

      
!theta2 = pi/3  (60)      
      cosP(1) = rHALF
      sinP(1) = sqrt(3.0d0)*rHalf
!theta3 = pi/2  (90)
      cosP(2) = rZERO
      sinP(2) = rONE
! theta6 = 2pi/3 (120)     
      cosP(3) = -rHALF
      sinP(3) = sqrt(3.0d0)*rHalf
                  
      do i = 1,3
         p3(0,i) = q3(0)
         p3(1,i) = sqrt(q3T)*sinP(i)
         p3(2,i) = rZero               
         p3(3,i) = sqrt(q3T)*cosP(i)
         
         p4(0,i) = q4(0)
         p4(1,i) = - p3(1,i)
         p4(2,i) = rZero        
         p4(3,i) = - p3(3,i)         
      enddo        
       
      do i = 1,3
         call voffxx(p3(0,i),WMASS,+1,Nq3,polv3(1,-1,i))
         call voffxx(p4(0,i),WMASS,+1,Nq4,polv4(1,-1,i))            
      enddo            
      
      do i=1,3
      CALL wpwmwmwp_ONanom4(polv1(1,-1),polv3(1,-1,i),polv2(1,-1),polv4(1,-1,i),M(-1,-1,-1,-1,1,1,i))
      enddo 
      do lambda1= -1,1
        do lambda2= -1,1
          do lambda3= -1,1
            do lambda4= -1,1
              do i=1,3
                Msorted(lambda1,lambda3,lambda2,lambda4,1,1, i)=M(lambda1,lambda2,lambda3,lambda4,1,1,i)
              enddo
            enddo
          enddo
        enddo
      enddo      
      
      CALL inv_Wignervv(Dinv)   !  For WW->WW Dinv->Dinv/2.0d0 (unitarity circle)
                        
      do lambda1 = -1,1
       do lambda2 = -1,1
        do lambda3 = -1,1
         do lambda4 = -1,1
          do i= 0,2
            a(lambda1,lambda2,lambda3, lambda4,1,1,i) = 1/2*((Dinv(lambda1-lambda2,lambda3-lambda4,i,1))*Msorted(lambda1,lambda2,lambda3,lambda4,1,1,1) 
     &                                                            +(Dinv(lambda1-lambda2,lambda3-lambda4,i,2))*Msorted(lambda1,lambda2,lambda3,lambda4,1,1,2)
     &                                                            +(Dinv(lambda1-lambda2,lambda3-lambda4,i,3))*Msorted(lambda1,lambda2,lambda3,lambda4,1,1,3))

          enddo
         enddo
        enddo
       enddo
      enddo       

      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Wigners d-matrix

      SUBROUTINE WIGNERdvv(costheta,smalld)
      
      IMPLICIT NONE
      
      DOUBLE PRECISION costheta  ! Cos of the angle
      DOUBLE PRECISION ct, st 
      DOUBLE PRECISION smalld(0:2,-2:2,-2:2) ! d-functions d(l,alpha,beta)
      
      INTEGER j,l1,l2
      
      DOUBLE PRECISION rZERO, rONE, rHalf
      PARAMETER( rZERO=0.0D0, rONE=1.0D0, rHALF = 0.5D0 )
      
      ct = costheta
      st = sqrt(rONE - costheta**2)
      
      do j = 0,2
       do l1 = -2,2
        do l2 = -2,2
         smalld(j,l1,l2) = rZero
        enddo
       enddo
      enddo
      
      smalld(0,0,0) = rONE
      
      smalld(1,0,0) = CT
      
      smalld(1,1,0) = - sqrt(rHalf)*ST
      smalld(1,0,1) = (-rONE)* smalld(1,1,0)
      smalld(1,-1,0) = (-rONE)*smalld(1,1,0) 
      smalld(1,0,-1) = smalld(1,1,0)
      
      smalld(1,1,1) = rHalf*(rONE + CT)
      smalld(1,1,-1) = rHalf*(rONE - CT)      
      smalld(1,-1,-1) = smalld(1,1,1)         
      smalld(1,-1,1) = smalld(1,1,-1)         
      
      smalld(2,0,0) = rHalf*(3.0d0*CT**2 - rONE)
      
      smalld(2,1,0) = - sqrt(3.0d0*rHalf)*ST*CT
      smalld(2,0,1) = (-rONE)*smalld(2,1,0)
      smalld(2,0,-1) = smalld(2,1,0)
      smalld(2,-1,0) = (-rONE)*smalld(2,1,0)
      
      smalld(2,1,1) = rHalf * (rONE + CT) * (2.0d0*CT - rONE)
      smalld(2,1,-1) = rHalf * (rONE - CT) * (2.0d0*CT + rONE)
      smalld(2,-1,-1) = smalld(2,1,1)            
      smalld(2,-1,1) = smalld(2,1,-1)            
      
      smalld(2,2,1) = - rHalf*(rONE + CT)* ST
      smalld(2,1,2) = - smalld(2,2,1)
      smalld(2,-2,-1) = -smalld(2,2,1)
      smalld(2,-1,-2) = smalld(2,2,1)      
      
      smalld(2,2,-1) = - rHalf*(rONE - CT)* ST
      smalld(2,-2,1) = - smalld(2,2,-1)
      smalld(2,1,-2) = smalld(2,2,-1)
      smalld(2,-1,2) = - smalld(2,2,-1)
      
      smalld(2,2,2) = (rHALF*(rONE + CT))**2
      smalld(2,2,-2) = (rHALF*(rONE - CT))**2
      smalld(2,-2,-2) = smalld(2,2,2)
      smalld(2,-2,2) = smalld(2,2,-2) 
      
      smalld(2,2,0) = (sqrt(6.0d0)/4.0d0) * ST**2
      smalld(2,-2,0) = smalld(2,2,0)
      smalld(2,0,2) = smalld(2,2,0)
      smalld(2,0,-2) = smalld(2,0,2)      
           
      END     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE inv_Wignervv(Dinv)

c  SUBROUTINE inverse Wigner d-Matrix for the wigner subroutine
c  for theta 1 = 0, 2 = PI/3, 3 = PI/2             
c  Dinv = (D(lambda1-lambda2, lambda3-lambda4, j, theta))^{-1}
c  A (i) = Dinv(lambda1-lambda2, lambda3 - lambda4, i , j)*M(j)  j=theta(j) 
      IMPLICIT NONE
      
      DOUBLE PRECISION Dinv(-2:2,-2:2,0:2,3)
      
      DOUBLE PRECISION COST1, COST2, COST3
      DOUBLE PRECISION SINT1, SINT2, SINT3
      DOUBLE PRECISION SECT
      DOUBLE PRECISION rOne, rZero, rPi, rHalf, cinv, sqrt2
      PARAMETER( rOne = 1.0d0, rZero = 0.0d0, rPi = 3.14159265358979323846d0 , rHalf = 0.5d0)
      PARAMETER( cinv = rOne/(8.0d0*rPi), sqrt2 = sqrt(2.0d0) ) 
      
      Dinv(-2:2, -2:2, 0:2, 1:3) = rZero
! 60 (pi/3)      
      COST1 = rHalf
      SINT1 = sqrt(3.0d0)/2.0d0
! 90 (pi/2)      
      COST2 = rZero
      SINT2 = rOne
! 120 (2pi/3)      
      COST3 = -rHalf
      SINT3 = sqrt(3.0d0)/2.0d0
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 0,0
      Dinv(0,0,0,1) = cinv*2.0d0/3.0d0
      Dinv(0,0,0,2) = -cinv*rONE/3.0d0
      Dinv(0,0,0,3) = cinv*2.0d0/3.0d0
      Dinv(0,0,1,1) = cinv*rOne/3.0d0
      Dinv(0,0,1,2) = rZero
      Dinv(0,0,1,3) = -cinv*rONE/3.0d0
      Dinv(0,0,2,1) = cinv*4.0d0/15.0d0
      Dinv(0,0,2,2) = -cinv*8.0d0/15.0d0
      Dinv(0,0,2,3) = cinv*4.0d0/15.0d0
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1,1
      Dinv(1,1,0,1) = rZERO
      Dinv(1,1,0,2) = rZero
      Dinv(1,1,0,3) = rZero
      Dinv(1,1,1,1) = cinv*4.0d0/9.0d0
      Dinv(1,1,1,2) = rZero
      Dinv(1,1,1,3) = rZero
      Dinv(1,1,2,1) = cinv*4.0d0/15.0d0
      Dinv(1,1,2,2) = -cinv*2.0d0/5.0d0
      Dinv(1,1,2,3) = rZero

      Dinv(-1,-1,0:2,1:3) = Dinv(1,1,0:2,1:3)    

! 1,-1      
      Dinv(1,-1,0,1) = rZERO
      Dinv(1,-1,0,2) = rZero
      Dinv(1,-1,0,3) = rZero
      Dinv(1,-1,1,1) = -cinv*4.0d0/3.0d0
      Dinv(1,-1,1,2) = cinv*4.0d0/3.0d0
      Dinv(1,-1,1,3) = rZero
      Dinv(1,-1,2,1) = cinv*4.0d0/5.0d0
      Dinv(1,-1,2,2) = -cinv*2.0d0/5.0d0
      Dinv(1,-1,2,3) = rZero
      
      Dinv(-1,1,1:2,1:3) = Dinv(1,-1,1:2,1:3)        

! 1,0
      Dinv(1,0,0,1) = rZERO
      Dinv(1,0,0,2) = rZero
      Dinv(1,0,0,3) = rZero
      Dinv(1,0,1,1) = rZero
      Dinv(1,0,1,2) = -cinv*sqrt(2.0d0)/3.0d0
      Dinv(1,0,1,3) = rZero
      Dinv(1,0,2,1) = -cinv*4.0d0*sqrt(2.0d0)/15.0d0
      Dinv(1,0,2,2) = cinv*2.0d0*sqrt(2.0d0/3.0d0)/5.0d0
      Dinv(1,0,2,3) = rZero
      
      Dinv(-1,0,0:2,1:3) = - Dinv(1,0,0:2,1:3)
      
      Dinv(0,1,0:2,1:3) = - Dinv(1,0,0:2,1:3)
      
      Dinv(0,-1,0:2,1:3) = Dinv(1,0,0:2,1:3)
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
! 2,2
      Dinv(2,2,2,1)   = cinv*(rOne/5.0d0)*(((rONE + COST1)*rHalf)**2)**(-1)

      Dinv(-2,-2,2,1) = Dinv(2,2,2,1)
! 2,1      
      Dinv(2,1,2,1)   = - cinv*(rOne/5.0d0)* (rHALF*(rONE+COST1)*SINT1)**(-1)      

      Dinv(-1,-2,2,1) = Dinv(2,1,2,1)
      
      Dinv(1,2,2,1)   = - Dinv(2,1,2,1)
      
      Dinv(-2,-1,2,1) = - Dinv(2,1,2,1)
! 2,0      
      Dinv(2,0,2,1)   = cinv*(rOne/5.0d0)*(sqrt(6.0D0)*SINT1**2 /4.0D0)**(-1)
      
      Dinv(0,-2,2,1)  = Dinv(2,0,2,1)
      
      Dinv(-2,0,2,1)  = Dinv(2,0,2,1)
      
      Dinv(0,2,2,1)   = Dinv(2,0,2,1)
! 2,-1      
      Dinv(2,-1,2,1)  = -cinv*(rOne/5.0d0)*((rONE-COST1)*rHALF*SINT1)**(-1)
      
      Dinv(-2,1,2,1)  = - Dinv(2,-1,2,1)
      
      Dinv(-1,2,2,1)  = -Dinv(2,-1,2,1)
      
      Dinv(1,-2,2,1)  = Dinv(2,-1,2,1)
! 2,-2
      Dinv(2,-2,2,1) = cinv*(rOne/5.0d0)*(((rONE-COST1)*rHALF)**2)**(-1)
      
      Dinv(-2,2,2,1) = Dinv(2,-2,2,1)
      
      END      


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE sortmatrixpol(charge,rank,M,Msorted)

      IMPLICIT NONE
      INTEGER RANK, lambda1, lambda2, lambda3, lambda4, block1, block2, i
      DOUBLE COMPLEX M(-1:1,-1:1,-1:1,-1:1,1:rank,1:rank,1:3)
      DOUBLE COMPLEX Msorted(-1:1,-1:1,-1:1,-1:1,1:rank,1:rank,1:3)      

      INTEGER casenumber,charge

      DOUBLE COMPLEX cZERO, cONE
      PARAMETER(cZERO=(0.0D0,0.0D0), cONE=(1.0D0,0.0D0))

      Msorted(-1:1,-1:1,-1:1,-1:1,1:rank,1:rank,1:3) = cZERO



      select case(charge) !0: WPWM, ZZ, HH, AA etc., 1: WPZ etc. , -1 WMZ etc., 2: WPWP, WMWM
      case(0)
      do lambda1= -1,1
        do lambda2= -1,1
          do lambda3= -1,1
            do lambda4= -1,1
              do i=1,3
                do block1= 1,rank
                  do block2= 1,rank
                    casenumber =block2 + rank*(block1-1)
                    select case(casenumber)
                    case(1:3,5:6,9:10,16:17,20:21,31,32,49,50,53,54,64,65,97,98,120)            
                    Msorted(lambda1,lambda2,lambda3,lambda4,block1,block2,i)=M(lambda1,lambda2,lambda3,lambda4,block1,block2,i)
                    case(4,7:8,11,19,22,33,52,55,66,99,121)            
                    Msorted(lambda1,lambda2,lambda4,lambda3,block1,block2,i)=M(lambda1,lambda2,lambda3,lambda4,block1,block2,i)
                    case(13:14,25,28,61)            
                    Msorted(lambda1,lambda3,lambda2,lambda4,block1,block2,i)=M(lambda1,lambda2,lambda3,lambda4,block1,block2,i)
                    case(15,18,29,30,63)            
                    Msorted(lambda1,lambda3,lambda4,lambda2,block1,block2,i)=M(lambda1,lambda2,lambda3,lambda4,block1,block2,i)
                    case(12,23:24,27,60)            
                    Msorted(lambda1,lambda4,lambda2,lambda3,block1,block2,i)=M(lambda1,lambda2,lambda3,lambda4,block1,block2,i)
                    case(26)            
                    Msorted(lambda1,lambda4,lambda3,lambda2,block1,block2,i)=M(lambda1,lambda2,lambda3,lambda4,block1,block2,i)
                    case(38,39,42,43,75,76,86,87,109)           
                    Msorted(lambda2,lambda1,lambda3,lambda4,block1,block2,i)=M(lambda1,lambda2,lambda3,lambda4,block1,block2,i)
                    case(41,44,77,88,110 )           
                    Msorted(lambda2,lambda1,lambda4,lambda3,block1,block2,i)=M(lambda1,lambda2,lambda3,lambda4,block1,block2,i)
                    case(46,47,58,91,94)            
                    Msorted(lambda2,lambda3,lambda1,lambda4,block1,block2,i)=M(lambda1,lambda2,lambda3,lambda4,block1,block2,i)
                    case(48,51,62,95,96)            
                    Msorted(lambda2,lambda3,lambda4,lambda1,block1,block2,i)=M(lambda1,lambda2,lambda3,lambda4,block1,block2,i)
                    case(57)            
                    Msorted(lambda2,lambda4,lambda1,lambda3,block1,block2,i)=M(lambda1,lambda2,lambda3,lambda4,block1,block2,i)
                    case(59)            
                    Msorted(lambda2,lambda4,lambda3,lambda1,block1,block2,i)=M(lambda1,lambda2,lambda3,lambda4,block1,block2,i)
                    case(35,36,69,72,83 )           
                    Msorted(lambda3,lambda1,lambda2,lambda4,block1,block2,i)=M(lambda1,lambda2,lambda3,lambda4,block1,block2,i)
                    case(37,40,73,74,85)            
                    Msorted(lambda3,lambda1,lambda4,lambda2,block1,block2,i)=M(lambda1,lambda2,lambda3,lambda4,block1,block2,i)
                    case(80 )           
                    Msorted(lambda3,lambda2,lambda1,lambda4,block1,block2,i)=M(lambda1,lambda2,lambda3,lambda4,block1,block2,i)
                    case(84 )           
                    Msorted(lambda3,lambda2,lambda4,lambda1,block1,block2,i)=M(lambda1,lambda2,lambda3,lambda4,block1,block2,i)
                    case(45,56,89,90,92,93,100,101,102,104,105,108)            
                    Msorted(lambda3,lambda4,lambda1,lambda2,block1,block2,i)=M(lambda1,lambda2,lambda3,lambda4,block1,block2,i)
                    case(103,106,107)           
                    Msorted(lambda3,lambda4,lambda2,lambda1,block1,block2,i)=M(lambda1,lambda2,lambda3,lambda4,block1,block2,i)
                    case(34,67,68,71,82 )           
                    Msorted(lambda4,lambda1,lambda2,lambda3,block1,block2,i)=M(lambda1,lambda2,lambda3,lambda4,block1,block2,i)
                    case(70 )           
                    Msorted(lambda4,lambda1,lambda3,lambda2,block1,block2,i)=M(lambda1,lambda2,lambda3,lambda4,block1,block2,i)
                    case(79)            
                    Msorted(lambda4,lambda2,lambda1,lambda3,block1,block2,i)=M(lambda1,lambda2,lambda3,lambda4,block1,block2,i)
                    case(81 )           
                    Msorted(lambda4,lambda2,lambda3,lambda1,block1,block2,i)=M(lambda1,lambda2,lambda3,lambda4,block1,block2,i)
                    case(78,111,112,113,115,116,119 )           
                    Msorted(lambda4,lambda3,lambda1,lambda2,block1,block2,i)=M(lambda1,lambda2,lambda3,lambda4,block1,block2,i)
                    case(114,117,118 )           
                    Msorted(lambda4,lambda3,lambda2,lambda1,block1,block2,i)=M(lambda1,lambda2,lambda3,lambda4,block1,block2,i)
                    end select
                  enddo
                enddo              
              enddo
            enddo
          enddo
        enddo
      enddo
      case(1)
      do lambda1= -1,1
        do lambda2= -1,1
          do lambda3= -1,1
            do lambda4= -1,1
              do i=1,3
                do block1= 1,rank
                  do block2= 1,rank
                    casenumber =block2 + rank*(block1-1)
                    select case(casenumber)
                    case(1,2,3,5,6,9)
                    Msorted(lambda4,lambda1,lambda3,lambda2,block1,block2,i)=M(lambda1,lambda2,lambda3,lambda4,block1,block2,i)
                    case(4,7,8)
                    Msorted(lambda4,lambda2,lambda3,lambda1,block1,block2,i)=M(lambda1,lambda2,lambda3,lambda4,block1,block2,i)                    
                    end select
                  enddo
                enddo              
              enddo
            enddo
          enddo
        enddo
      enddo
      case(-1)
        do lambda1= -1,1
          do lambda2= -1,1
            do lambda3= -1,1
              do lambda4= -1,1
                do i=1,3
                  do block1= 1,rank
                    do block2= 1,rank
                    casenumber =block2 + rank*(block1-1)
                    select case(casenumber)
                    case(1,2,3,5,6,9)
                    Msorted(lambda3,lambda1,lambda4,lambda2,block1,block2,i)=M(lambda1,lambda2,lambda3,lambda4,block1,block2,i)
                    case(4,7,8)
                    Msorted(lambda3,lambda2,lambda4,lambda1,block1,block2,i)=M(lambda1,lambda2,lambda3,lambda4,block1,block2,i)                    
                    end select
                    enddo
                  enddo              
                enddo
              enddo
            enddo
          enddo
        enddo      
      case(2)
        do lambda1= -1,1
          do lambda2= -1,1
            do lambda3= -1,1
              do lambda4= -1,1
                do i=1,3
                  do block1= 1,rank
                    do block2= 1,rank
                        Msorted(lambda1,lambda4,lambda2,lambda3,block1,block2,i)=M(lambda1,lambda2,lambda3,lambda4,block1,block2,i)   
                      

                    enddo
                  enddo              
                enddo
              enddo
            enddo
          enddo
        enddo             
      end select
      
      
      END

 