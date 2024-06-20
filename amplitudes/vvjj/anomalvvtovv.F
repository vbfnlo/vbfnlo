c
c This file contains routines for quartic-vector boson couplings
c with anomalous couplings for on-shell polarization vectors
c in the CM-Frame, as needed by ????_on.F.
c   subroutine aaaa_ONanom4(a1,a2,a3,a4,M)
c   subroutine aaah_ONanom4(a1,a2,a3,h,M)
c   subroutine aaaz_ONanom4(a1,a2,a3,z,M)
c   subroutine aahh_ONanom4(a1,a2,h1,h2,M)
c   subroutine aahz_ONanom4(a1,a2,h,z,M)
c   subroutine aazz_ONanom4(a1,a2,z1,z1,M)
c   subroutine aawmwp_ONanom4(a1,a2,wm1,wp1,M)
c   subroutine ahhh_ONanom4(a1,h1,h2,h3,M)
c   subroutine ahhz_ONanom4(a1,h1,h2,z,M)
c   subroutine ahzz_ONanom4(a1,h,z1,z2,M)
c   subroutine ahwmwp_ONanom4(a1,h,wm1,wp1,M)
c   subroutine azzz_ONanom4(a1,z1,z2,z3,M)
c   subroutine azwmwp_ONanom4(a1,z,wm1,wp1,M)
c   subroutine hhhh_ONanom4(h1,h2,h3,h4,M)
c   subroutine hhhz_ONanom4(h1,h2,h3,z,M)
c   subroutine hhzz_ONanom4(h1,h2,z1,z2,M)
c   subroutine hhwmwp_ONanom4(h1,h2,wm1,wz,M)
c   subroutine hzzz_ONanom4(h,z1,z2,z3,M)
c   subroutine hzwmwp_ONanom4(h,z,wm1,wp1,M)
c   subroutine zzzz_ONanom4(z1,z2,z3,z4,M)
c   subroutine zzwmwp_ONanom4(z1,z2,wm1,wp1,M)
c   subroutine wpwmwmwp_ONanom4(wp1,wm1,wm2,wp2,M)
c
c   
c
ccccccccccccccc
c aaaa-vertex c
ccccccccccccccc

      subroutine aaaa_ONanom4(a1,a2,a3,a4,M)
c
c This subroutine calculates the vvvv-vertex with anomalous couplings
c and the on-shell polarization vectors for vector bosons
c
c IN:
c     polv1, polv2      : the incoming polvector wz(bosons in the CM frame)
c     polv3, polv4      : the outgoing polvector wz(bosons in the CM frame)
c OUT:
c     M(lambda1,lambda2) 
c
      IMPLICIT NONE

      INTEGER NGRAPHS
      PARAMETER( NGRAPHS=20 )
      DOUBLE PRECISION rZERO, rONE, rHalf
      PARAMETER( rZERO=0.0D0, rONE=1.0D0, rHALF = 0.5D0 )
      DOUBLE COMPLEX cIMAG, cZERO, cONE
      PARAMETER( cIMAG=(0.0D0,1.0D0), cZERO=(0.0D0,0.0D0), cONE=(1.0D0,0.0D0))
      DOUBLE COMPLEX a1(6,-1:1), a2(6,-1:1), a3(6,-1:1), a4(6,-1:1)      
      DOUBLE COMPLEX M(-1:1,-1:1,-1:1,-1:1)
      
      DOUBLE COMPLEX OFFeps(6,-1:1,NGRAPHS), OFFepsC(6,-1:1,NGRAPHS)
      DOUBLE PRECISION K(0:3,4)      
      DOUBLE COMPLEX O(-1:1,4,-1:1,4), OP(4,-1:1,4)
      DOUBLE COMPLEXPP(1:4,1:4)
    
      INTEGER I,J,N,MU,NU
      INTEGER LAMBDA1, LAMBDA2, LAMBDA3, LAMBDA4, LAMBDA
      INTEGER L1,L2

      DOUBLE COMPLEX DOTCC, DOTRC,contract_Tjj
      DOUBLE PRECISION DOTRR
      EXTERNAL DOTCC, DOTRC, DOTRR, contract_Tjj
      double precision sw2, sw
#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "gmunu.inc"
#include "transm.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"
      sw2 = sin2w
      sw = sqrt(sin2w)
      OFFeps(6,-1:1,4) = cZero
      K(0:3,1:4) = rZero
      M(-1:1,-1:1,-1:1,-1:1) = cZero
      PP(1:4,1:4) = rZero
      O(-1:1,1:4,-1:1,1:4) = cZero
      OP(1:4,-1:1,1:4) = cZero  
      
      OFFeps(1:6,-1:1,1) = a1(1:6,-1:1)
      OFFeps(1:6,-1:1,2) = a2(1:6,-1:1)
      OFFeps(1:6,-1:1,3) = a3(1:6,-1:1)
      OFFeps(1:6,-1:1,4) = a4(1:6,-1:1) 


      do i = 1,4
        K(0,i) = dble(OFFeps(5,-1,i))*(-rONE) !K1(0:3)
        K(1,i) = dble(OFFeps(6,-1,i))*(-rONE) !K2(0:3)
        K(2,i) = dimag(OFFeps(6,-1,i))*(-rONE)!K3(0:3)
        K(3,i) = dimag(OFFeps(5,-1,i))*(-rONE)!K4(0:3)
      enddo      



      do i = 1,4
        do j = 1,4
          pp(i,j) = dotrr(k(0,i),k(0,j))
        enddo
      enddo

      do L1 = -1,1
        do L2 = -1,1 
          do i = 1,4
            do j = 1,4
              O(L1,i,L2,j) = dotcc(OFFeps(1,L1,i),OFFeps(1,L2,j))
            enddo
          enddo
        enddo
      enddo

      do LAMBDA = -1,1
        do i = 1,4 
          do j = 1,4
            OP(j,lambda,i) = dotrc(k(0,j),OFFeps(1,lambda,i))
          enddo
        enddo
      enddo


!!!! ===============================================================      
!!!! COUPLINGS CALCULATION
!!!! ===============================================================


                 

         if( ft0.ne.0d0 ) then  !ft0 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) + 8*ft0*(gw**4)*(sw**4)*(OP(1,lambda2,2)*OP(2,lambda1,1)*OP(3,lambda4,4)*OP(4,lambda3,3)-OP(3,lambda4,4)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda2,2)-OP(2,lambda4,4)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda3,3)-OP(2,lambda3,3)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*OP(4,lambda1,1)*(OP(2,lambda3,3)*OP(3,lambda2,2)-pp(2,3)*O(lambda2,2,lambda3,3))+pp(1,3)*pp(2,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*OP(3,lambda1,1)*(OP(2,lambda4,4)*OP(4,lambda2,2)-pp(2,4)*O(lambda2,2,lambda4,4))-OP(1,lambda2,2)*OP(2,lambda1,1)*pp(3,4)*O(lambda3,3,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4))
               enddo                      
              enddo
             enddo 
            enddo
         endif


         if( ft1.ne.0d0 ) then  !ft1 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +8*ft1*(gw**4)*(sw**4)*(OP(1,lambda2,2)*OP(2,lambda1,1)*OP(3,lambda4,4)*OP(4,lambda3,3)-OP(3,lambda4,4)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda2,2)-OP(2,lambda4,4)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda3,3)-OP(2,lambda3,3)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*OP(4,lambda1,1)*(OP(2,lambda3,3)*OP(3,lambda2,2)-pp(2,3)*O(lambda2,2,lambda3,3))+pp(1,3)*pp(2,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*OP(3,lambda1,1)*(OP(2,lambda4,4)*OP(4,lambda2,2)-pp(2,4)*O(lambda2,2,lambda4,4))-OP(1,lambda2,2)*OP(2,lambda1,1)*pp(3,4)*O(lambda3,3,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4))

               enddo                      
              enddo
             enddo 
            enddo
         endif


         if( ft2.ne.0d0 ) then  !ft2 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +2*ft2*(gw**4)*(sw**4)*(OP(1,lambda4,4)*OP(2,lambda3,3)*OP(3,lambda1,1)*OP(4,lambda2,2)+OP(1,lambda4,4)*OP(2,lambda1,1)*OP(3,lambda2,2)*OP(4,lambda3,3)-OP(2,lambda4,4)*OP(4,lambda3,3)*pp(1,3)*O(lambda1,1,lambda2,2)-OP(2,lambda3,3)*OP(3,lambda4,4)*pp(1,4)*O(lambda1,1,lambda2,2)-OP(1,lambda4,4)*OP(4,lambda3,3)*pp(2,3)*O(lambda1,1,lambda2,2)+OP(1,lambda4,4)*OP(2,lambda3,3)*pp(3,4)*O(lambda1,1,lambda2,2)-OP(3,lambda4,4)*OP(4,lambda2,2)*pp(1,2)*O(lambda1,1,lambda3,3)-OP(2,lambda4,4)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda3,3)-OP(1,lambda4,4)*OP(4,lambda2,2)*pp(2,3)*O(lambda1,1,lambda3,3)+OP(1,lambda4,4)*OP(3,lambda2,2)*pp(2,4)*O(lambda1,1,lambda3,3)-OP(3,lambda2,2)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda4,4)-OP(2,lambda3,3)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda4,4)-OP(3,lambda4,4)*OP(4,lambda1,1)*pp(1,2)*O(lambda2,2,lambda3,3)-OP(2,lambda4,4)*OP(4,lambda1,1)*pp(1,3)*O(lambda2,2,lambda3,3)+OP(2,lambda4,4)*OP(3,lambda1,1)*pp(1,4)*O(lambda2,2,lambda3,3)+OP(2,lambda1,1)*OP(3,lambda4,4)*pp(1,4)*O(lambda2,2,lambda3,3)-OP(1,lambda4,4)*OP(3,lambda1,1)*pp(2,4)*O(lambda2,2,lambda3,3)-OP(1,lambda4,4)*OP(2,lambda1,1)*pp(3,4)*O(lambda2,2,lambda3,3)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-OP(3,lambda1,1)*OP(4,lambda3,3)*pp(1,2)*O(lambda2,2,lambda4,4)+OP(2,lambda3,3)*OP(4,lambda1,1)*pp(1,3)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*OP(4,lambda3,3)*pp(1,3)*O(lambda2,2,lambda4,4)-OP(2,lambda3,3)*OP(3,lambda1,1)*pp(1,4)*O(lambda2,2,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(-(OP(3,lambda4,4)*pp(2,4)*O(lambda1,1,lambda2,2))+OP(2,lambda4,4)*(OP(3,lambda2,2)*OP(4,lambda1,1)+pp(3,4)*O(lambda1,1,lambda2,2))+OP(4,lambda2,2)*pp(2,3)*O(lambda1,1,lambda4,4)-OP(3,lambda2,2)*pp(2,4)*O(lambda1,1,lambda4,4)-OP(4,lambda1,1)*pp(2,3)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*(OP(3,lambda4,4)*OP(4,lambda2,2)-pp(3,4)*O(lambda2,2,lambda4,4)))+OP(3,lambda2,2)*OP(4,lambda1,1)*pp(1,2)*O(lambda3,3,lambda4,4)+OP(3,lambda1,1)*OP(4,lambda2,2)*pp(1,2)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(4,lambda2,2)*pp(1,3)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(3,lambda2,2)*pp(1,4)*O(lambda3,3,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+OP(1,lambda2,2)*(OP(3,lambda4,4)*pp(2,4)*O(lambda1,1,lambda3,3)+OP(2,lambda4,4)*(OP(3,lambda1,1)*OP(4,lambda3,3)-pp(3,4)*O(lambda1,1,lambda3,3))+OP(4,lambda3,3)*pp(2,3)*O(lambda1,1,lambda4,4)+OP(2,lambda3,3)*(OP(3,lambda4,4)*OP(4,lambda1,1)-pp(3,4)*O(lambda1,1,lambda4,4))-OP(4,lambda1,1)*pp(2,3)*O(lambda3,3,lambda4,4)-OP(3,lambda1,1)*pp(2,4)*O(lambda3,3,lambda4,4)))

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( ft3.ne.0d0 ) then  !ft3 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +2*ft3*gw**4*sw**4*(-pp(1,2)*OP(4,lambda2,2)*OP(3,lambda4,4)*O(lambda1,1,lambda3,3)-pp(2,3)*OP(4,lambda2,2)*OP(1,lambda4,4)*O(lambda1,1,lambda3,3)-pp(1,3)*OP(4,lambda2,2)*OP(2,lambda3,3)*O(lambda1,1,lambda4,4)+pp(1,2)*OP(3,lambda1,1)*OP(4,lambda2,2)*O(lambda3,3,lambda4,4)-pp(1,3)*OP(2,lambda1,1)*OP(4,lambda2,2)*O(lambda3,3,lambda4,4)-pp(1,3)*OP(4,lambda3,3)*OP(2,lambda4,4)*O(lambda1,1,lambda2,2)-pp(1,4)*OP(2,lambda3,3)*OP(3,lambda4,4)*O(lambda1,1,lambda2,2)-pp(2,3)*OP(4,lambda3,3)*OP(1,lambda4,4)*O(lambda1,1,lambda2,2)+pp(3,4)*OP(2,lambda3,3)*OP(1,lambda4,4)*O(lambda1,1,lambda2,2)-pp(1,4)*OP(3,lambda2,2)*OP(2,lambda4,4)*O(lambda1,1,lambda3,3)+pp(2,4)*OP(3,lambda2,2)*OP(1,lambda4,4)*O(lambda1,1,lambda3,3)-pp(1,2)*OP(3,lambda2,2)*OP(4,lambda3,3)*O(lambda1,1,lambda4,4)-pp(1,2)*OP(4,lambda1,1)*OP(3,lambda4,4)*O(lambda2,2,lambda3,3)-pp(1,3)*OP(4,lambda1,1)*OP(2,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,4)*OP(3,lambda1,1)*OP(2,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,4)*OP(2,lambda1,1)*OP(3,lambda4,4)*O(lambda2,2,lambda3,3)-pp(2,4)*OP(3,lambda1,1)*OP(1,lambda4,4)*O(lambda2,2,lambda3,3)-pp(3,4)*OP(2,lambda1,1)*OP(1,lambda4,4)*O(lambda2,2,lambda3,3)-pp(1,2)*OP(3,lambda1,1)*OP(4,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,3)*OP(4,lambda1,1)*OP(2,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,3)*OP(2,lambda1,1)*OP(4,lambda3,3)*O(lambda2,2,lambda4,4)-pp(1,4)*OP(3,lambda1,1)*OP(2,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(pp(2,4)*(-OP(3,lambda4,4))*O(lambda1,1,lambda2,2)+OP(2,lambda4,4)*(OP(4,lambda1,1)*OP(3,lambda2,2)+pp(3,4)*O(lambda1,1,lambda2,2))+pp(2,3)*OP(4,lambda2,2)*O(lambda1,1,lambda4,4)-pp(2,4)*OP(3,lambda2,2)*O(lambda1,1,lambda4,4)-pp(2,3)*OP(4,lambda1,1)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*(OP(4,lambda2,2)*OP(3,lambda4,4)-pp(3,4)*O(lambda2,2,lambda4,4)))+pp(1,2)*OP(4,lambda1,1)*OP(3,lambda2,2)*O(lambda3,3,lambda4,4)-pp(1,4)*OP(2,lambda1,1)*OP(3,lambda2,2)*O(lambda3,3,lambda4,4)+OP(1,lambda2,2)*(pp(2,4)*OP(3,lambda4,4)*O(lambda1,1,lambda3,3)+OP(2,lambda4,4)*(OP(3,lambda1,1)*OP(4,lambda3,3)-pp(3,4)*O(lambda1,1,lambda3,3))+pp(2,3)*OP(4,lambda3,3)*O(lambda1,1,lambda4,4)+OP(2,lambda3,3)*(OP(4,lambda1,1)*OP(3,lambda4,4)-pp(3,4)*O(lambda1,1,lambda4,4))-pp(2,3)*OP(4,lambda1,1)*O(lambda3,3,lambda4,4)-pp(2,4)*OP(3,lambda1,1)*O(lambda3,3,lambda4,4))+OP(3,lambda1,1)*OP(4,lambda2,2)*OP(2,lambda3,3)*OP(1,lambda4,4)+OP(2,lambda1,1)*OP(3,lambda2,2)*OP(4,lambda3,3)*OP(1,lambda4,4)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4))
               enddo                      
              enddo
             enddo 
            enddo
         endif


         if( ft4.ne.0d0 ) then  !ft4 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4)+gw**4*ft4*sw**4*(OP(1, lambda4, 4)*OP(2, lambda3, 3)*OP(3, lambda1, 1)*OP(4, lambda2, 2) + OP(1, lambda4, 4)*OP(2, lambda1, 1)*OP(3, lambda2, 2)*OP(4, lambda3, 3) + pp(3, 4)*OP(1, lambda4, 4)*OP(2, lambda3, 3)*O(lambda1, 1, lambda2, 2) - pp(1, 4)*OP(2, lambda3, 3)*OP(3, lambda4, 4)*O(lambda1, 1, lambda2, 2) - pp(2, 3)*OP(1, lambda4, 4)*OP(4, lambda3, 3)*O(lambda1, 1, lambda2, 2) - pp(1, 3)*OP(2, lambda4, 4)*OP(4, lambda3, 3)*O(lambda1, 1, lambda2, 2) + pp(2, 4)*OP(1, lambda4, 4)*OP(3, lambda2, 2)*O(lambda1, 1, lambda3, 3) - pp(1, 4)*OP(2, lambda4, 4)*OP(3, lambda2, 2)*O(lambda1, 1, lambda3, 3) - pp(2, 3)*OP(1, lambda4, 4)*OP(4, lambda2, 2)*O(lambda1, 1, lambda3, 3) - pp(1, 2)*OP(3, lambda4, 4)*OP(4, lambda2, 2)*O(lambda1, 1, lambda3, 3) - pp(1, 3)*OP(2, lambda3, 3)*OP(4, lambda2, 2)*O(lambda1, 1, lambda4, 4) - pp(1, 2)*OP(3, lambda2, 2)*OP(4, lambda3, 3)*O(lambda1, 1, lambda4, 4) - pp(3, 4)*OP(1, lambda4, 4)*OP(2, lambda1, 1)*O(lambda2, 2, lambda3, 3) - pp(2, 4)*OP(1, lambda4, 4)*OP(3, lambda1, 1)*O(lambda2, 2, lambda3, 3) + pp(1, 4)*OP(2, lambda4, 4)*OP(3, lambda1, 1)*O(lambda2, 2, lambda3, 3) + pp(1, 4)*OP(2, lambda1, 1)*OP(3, lambda4, 4)*O(lambda2, 2, lambda3, 3) - pp(1, 3)*OP(2, lambda4, 4)*OP(4, lambda1, 1)*O(lambda2, 2, lambda3, 3) - pp(1, 2)*OP(3, lambda4, 4)*OP(4, lambda1, 1)*O(lambda2, 2, lambda3, 3) + pp(1, 3)*pp(2, 4)*O(lambda1, 1, lambda4, 4)*O(lambda2, 2, lambda3, 3) + pp(1, 2)*pp(3, 4)*O(lambda1, 1, lambda4, 4)*O(lambda2, 2, lambda3, 3) - pp(1, 4)*OP(2, lambda3, 3)*OP(3, lambda1, 1)*O(lambda2, 2, lambda4, 4) + pp(1, 3)*OP(2, lambda3, 3)*OP(4, lambda1, 1)*O(lambda2, 2, lambda4, 4) + pp(1, 3)*OP(2, lambda1, 1)*OP(4, lambda3, 3)*O(lambda2, 2, lambda4, 4) - pp(1, 2)*OP(3, lambda1, 1)*OP(4, lambda3, 3)*O(lambda2, 2, lambda4, 4) + pp(1, 4)*pp(2, 3)*O(lambda1, 1, lambda3, 3)*O(lambda2, 2, lambda4, 4) + pp(1, 2)*pp(3, 4)*O(lambda1, 1, lambda3, 3)*O(lambda2, 2, lambda4, 4) + OP(1, lambda3, 3)*(-(pp(2, 4)*OP(3, lambda4, 4)*O(lambda1, 1, lambda2, 2)) + OP(2, lambda4, 4)*(OP(3, lambda2, 2)*OP(4, lambda1, 1) + pp(3, 4)*O(lambda1, 1, lambda2, 2)) -   pp(2, 4)*OP(3, lambda2, 2)*O(lambda1, 1, lambda4, 4) + pp(2, 3)*OP(4, lambda2, 2)*O(lambda1, 1, lambda4, 4) - pp(2, 3)*OP(4, lambda1, 1)*O(lambda2, 2, lambda4, 4) +   OP(2, lambda1, 1)*(OP(3, lambda4, 4)*OP(4, lambda2, 2) - pp(3, 4)*O(lambda2, 2, lambda4, 4))) - pp(1, 4)*OP(2, lambda1, 1)*OP(3, lambda2, 2)*O(lambda3, 3, lambda4, 4) + pp(1, 2)*OP(3, lambda2, 2)*OP(4, lambda1, 1)*O(lambda3, 3, lambda4, 4) - pp(1, 3)*OP(2, lambda1, 1)*OP(4, lambda2, 2)*O(lambda3, 3, lambda4, 4) + pp(1, 2)*OP(3, lambda1, 1)*OP(4, lambda2, 2)*O(lambda3, 3, lambda4, 4) + pp(1, 4)*pp(2, 3)*O(lambda1, 1, lambda2, 2)*O(lambda3, 3, lambda4, 4) + pp(1, 3)*pp(2, 4)*O(lambda1, 1, lambda2, 2)*O(lambda3, 3, lambda4, 4) + OP(1, lambda2, 2)*(pp(2, 4)*OP(3, lambda4, 4)*O(lambda1, 1, lambda3, 3) +   OP(2, lambda4, 4)*(OP(3, lambda1, 1)*OP(4, lambda3, 3) - pp(3, 4)*O(lambda1, 1, lambda3, 3)) + pp(2, 3)*OP(4, lambda3, 3)*O(lambda1, 1, lambda4, 4) +   OP(2, lambda3, 3)*(OP(3, lambda4, 4)*OP(4, lambda1, 1) - pp(3, 4)*O(lambda1, 1, lambda4, 4)) - pp(2, 4)*OP(3, lambda1, 1)*O(lambda3, 3, lambda4, 4) -   pp(2, 3)*OP(4, lambda1, 1)*O(lambda3, 3, lambda4, 4)))
 
               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( ft5.ne.0d0 ) then  !ft5 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(-4*ft5*(gw**4)*(sw**4)*(-1+SW2)*(ZMASS**2)*(OP(1,lambda2,2)*OP(2,lambda1,1)*OP(3,lambda4,4)*OP(4,lambda3,3)-OP(3,lambda4,4)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda2,2)-OP(2,lambda4,4)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda3,3)-OP(2,lambda3,3)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*OP(4,lambda1,1)*(OP(2,lambda3,3)*OP(3,lambda2,2)-pp(2,3)*O(lambda2,2,lambda3,3))+pp(1,3)*pp(2,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*OP(3,lambda1,1)*(OP(2,lambda4,4)*OP(4,lambda2,2)-pp(2,4)*O(lambda2,2,lambda4,4))-OP(1,lambda2,2)*OP(2,lambda1,1)*pp(3,4)*O(lambda3,3,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)))/(WMASS**2)

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( ft6.ne.0d0 ) then  !ft6 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(-4*ft6*(gw**4)*(sw**4)*(-1+SW2)*(ZMASS**2)*(OP(1,lambda2,2)*OP(2,lambda1,1)*OP(3,lambda4,4)*OP(4,lambda3,3)-OP(3,lambda4,4)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda2,2)-OP(2,lambda4,4)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda3,3)-OP(2,lambda3,3)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*OP(4,lambda1,1)*(OP(2,lambda3,3)*OP(3,lambda2,2)-pp(2,3)*O(lambda2,2,lambda3,3))+pp(1,3)*pp(2,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*OP(3,lambda1,1)*(OP(2,lambda4,4)*OP(4,lambda2,2)-pp(2,4)*O(lambda2,2,lambda4,4))-OP(1,lambda2,2)*OP(2,lambda1,1)*pp(3,4)*O(lambda3,3,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)))/(WMASS**2)

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( ft7.ne.0d0 ) then  !ft7 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) -((ft7*(gw**4)*(sw**4)*(-1+SW2)*(ZMASS**2)*(OP(1,lambda4,4)*OP(2,lambda3,3)*OP(3,lambda1,1)*OP(4,lambda2,2)+OP(1,lambda4,4)*OP(2,lambda1,1)*OP(3,lambda2,2)*OP(4,lambda3,3)-OP(2,lambda4,4)*OP(4,lambda3,3)*pp(1,3)*O(lambda1,1,lambda2,2)-OP(2,lambda3,3)*OP(3,lambda4,4)*pp(1,4)*O(lambda1,1,lambda2,2)-OP(1,lambda4,4)*OP(4,lambda3,3)*pp(2,3)*O(lambda1,1,lambda2,2)+OP(1,lambda4,4)*OP(2,lambda3,3)*pp(3,4)*O(lambda1,1,lambda2,2)-OP(3,lambda4,4)*OP(4,lambda2,2)*pp(1,2)*O(lambda1,1,lambda3,3)-OP(2,lambda4,4)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda3,3)-OP(1,lambda4,4)*OP(4,lambda2,2)*pp(2,3)*O(lambda1,1,lambda3,3)+OP(1,lambda4,4)*OP(3,lambda2,2)*pp(2,4)*O(lambda1,1,lambda3,3)-OP(3,lambda2,2)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda4,4)-OP(2,lambda3,3)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda4,4)-OP(3,lambda4,4)*OP(4,lambda1,1)*pp(1,2)*O(lambda2,2,lambda3,3)-OP(2,lambda4,4)*OP(4,lambda1,1)*pp(1,3)*O(lambda2,2,lambda3,3)+OP(2,lambda4,4)*OP(3,lambda1,1)*pp(1,4)*O(lambda2,2,lambda3,3)+OP(2,lambda1,1)*OP(3,lambda4,4)*pp(1,4)*O(lambda2,2,lambda3,3)-OP(1,lambda4,4)*OP(3,lambda1,1)*pp(2,4)*O(lambda2,2,lambda3,3)-OP(1,lambda4,4)*OP(2,lambda1,1)*pp(3,4)*O(lambda2,2,lambda3,3)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-OP(3,lambda1,1)*OP(4,lambda3,3)*pp(1,2)*O(lambda2,2,lambda4,4)+OP(2,lambda3,3)*OP(4,lambda1,1)*pp(1,3)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*OP(4,lambda3,3)*pp(1,3)*O(lambda2,2,lambda4,4)-OP(2,lambda3,3)*OP(3,lambda1,1)*pp(1,4)*O(lambda2,2,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(-(OP(3,lambda4,4)*pp(2,4)*O(lambda1,1,lambda2,2))+OP(2,lambda4,4)*(OP(3,lambda2,2)*OP(4,lambda1,1)+pp(3,4)*O(lambda1,1,lambda2,2))+OP(4,lambda2,2)*pp(2,3)*O(lambda1,1,lambda4,4)-OP(3,lambda2,2)*pp(2,4)*O(lambda1,1,lambda4,4)-OP(4,lambda1,1)*pp(2,3)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*(OP(3,lambda4,4)*OP(4,lambda2,2)-pp(3,4)*O(lambda2,2,lambda4,4)))+OP(3,lambda2,2)*OP(4,lambda1,1)*pp(1,2)*O(lambda3,3,lambda4,4)+OP(3,lambda1,1)*OP(4,lambda2,2)*pp(1,2)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(4,lambda2,2)*pp(1,3)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(3,lambda2,2)*pp(1,4)*O(lambda3,3,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+OP(1,lambda2,2)*(OP(3,lambda4,4)*pp(2,4)*O(lambda1,1,lambda3,3)+OP(2,lambda4,4)*(OP(3,lambda1,1)*OP(4,lambda3,3)-pp(3,4)*O(lambda1,1,lambda3,3))+OP(4,lambda3,3)*pp(2,3)*O(lambda1,1,lambda4,4)+OP(2,lambda3,3)*(OP(3,lambda4,4)*OP(4,lambda1,1)-pp(3,4)*O(lambda1,1,lambda4,4))-OP(4,lambda1,1)*pp(2,3)*O(lambda3,3,lambda4,4)-OP(3,lambda1,1)*pp(2,4)*O(lambda3,3,lambda4,4))))/(WMASS**2))

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( ft8.ne.0d0 ) then  !ft8 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(2*ft8*(gw**4)*(sw**4)*(-1+SW2)**2*(ZMASS**4)*(OP(1,lambda2,2)*OP(2,lambda1,1)*OP(3,lambda4,4)*OP(4,lambda3,3)-OP(3,lambda4,4)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda2,2)-OP(2,lambda4,4)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda3,3)-OP(2,lambda3,3)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*OP(4,lambda1,1)*(OP(2,lambda3,3)*OP(3,lambda2,2)-pp(2,3)*O(lambda2,2,lambda3,3))+pp(1,3)*pp(2,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*OP(3,lambda1,1)*(OP(2,lambda4,4)*OP(4,lambda2,2)-pp(2,4)*O(lambda2,2,lambda4,4))-OP(1,lambda2,2)*OP(2,lambda1,1)*pp(3,4)*O(lambda3,3,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)))/(WMASS**4)

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( ft9.ne.0d0 ) then  !ft9 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(ft9*(gw**4)*(sw**4)*(-1+SW2)**2*(ZMASS**4)*(OP(1,lambda4,4)*OP(2,lambda3,3)*OP(3,lambda1,1)*OP(4,lambda2,2)+OP(1,lambda4,4)*OP(2,lambda1,1)*OP(3,lambda2,2)*OP(4,lambda3,3)-OP(2,lambda4,4)*OP(4,lambda3,3)*pp(1,3)*O(lambda1,1,lambda2,2)-OP(2,lambda3,3)*OP(3,lambda4,4)*pp(1,4)*O(lambda1,1,lambda2,2)-OP(1,lambda4,4)*OP(4,lambda3,3)*pp(2,3)*O(lambda1,1,lambda2,2)+OP(1,lambda4,4)*OP(2,lambda3,3)*pp(3,4)*O(lambda1,1,lambda2,2)-OP(3,lambda4,4)*OP(4,lambda2,2)*pp(1,2)*O(lambda1,1,lambda3,3)-OP(2,lambda4,4)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda3,3)-OP(1,lambda4,4)*OP(4,lambda2,2)*pp(2,3)*O(lambda1,1,lambda3,3)+OP(1,lambda4,4)*OP(3,lambda2,2)*pp(2,4)*O(lambda1,1,lambda3,3)-OP(3,lambda2,2)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda4,4)-OP(2,lambda3,3)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda4,4)-OP(3,lambda4,4)*OP(4,lambda1,1)*pp(1,2)*O(lambda2,2,lambda3,3)-OP(2,lambda4,4)*OP(4,lambda1,1)*pp(1,3)*O(lambda2,2,lambda3,3)+OP(2,lambda4,4)*OP(3,lambda1,1)*pp(1,4)*O(lambda2,2,lambda3,3)+OP(2,lambda1,1)*OP(3,lambda4,4)*pp(1,4)*O(lambda2,2,lambda3,3)-OP(1,lambda4,4)*OP(3,lambda1,1)*pp(2,4)*O(lambda2,2,lambda3,3)-OP(1,lambda4,4)*OP(2,lambda1,1)*pp(3,4)*O(lambda2,2,lambda3,3)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-OP(3,lambda1,1)*OP(4,lambda3,3)*pp(1,2)*O(lambda2,2,lambda4,4)+OP(2,lambda3,3)*OP(4,lambda1,1)*pp(1,3)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*OP(4,lambda3,3)*pp(1,3)*O(lambda2,2,lambda4,4)-OP(2,lambda3,3)*OP(3,lambda1,1)*pp(1,4)*O(lambda2,2,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(-(OP(3,lambda4,4)*pp(2,4)*O(lambda1,1,lambda2,2))+OP(2,lambda4,4)*(OP(3,lambda2,2)*OP(4,lambda1,1)+pp(3,4)*O(lambda1,1,lambda2,2))+OP(4,lambda2,2)*pp(2,3)*O(lambda1,1,lambda4,4)-OP(3,lambda2,2)*pp(2,4)*O(lambda1,1,lambda4,4)-OP(4,lambda1,1)*pp(2,3)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*(OP(3,lambda4,4)*OP(4,lambda2,2)-pp(3,4)*O(lambda2,2,lambda4,4)))+OP(3,lambda2,2)*OP(4,lambda1,1)*pp(1,2)*O(lambda3,3,lambda4,4)+OP(3,lambda1,1)*OP(4,lambda2,2)*pp(1,2)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(4,lambda2,2)*pp(1,3)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(3,lambda2,2)*pp(1,4)*O(lambda3,3,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+OP(1,lambda2,2)*(OP(3,lambda4,4)*pp(2,4)*O(lambda1,1,lambda3,3)+OP(2,lambda4,4)*(OP(3,lambda1,1)*OP(4,lambda3,3)-pp(3,4)*O(lambda1,1,lambda3,3))+OP(4,lambda3,3)*pp(2,3)*O(lambda1,1,lambda4,4)+OP(2,lambda3,3)*(OP(3,lambda4,4)*OP(4,lambda1,1)-pp(3,4)*O(lambda1,1,lambda4,4))-OP(4,lambda1,1)*pp(2,3)*O(lambda3,3,lambda4,4)-OP(3,lambda1,1)*pp(2,4)*O(lambda3,3,lambda4,4))))/(2*(WMASS**4))

               enddo                      
              enddo
             enddo 
            enddo
         endif

      

      end



ccccccccccccccc
c aaah-vertex c
ccccccccccccccc

      subroutine aaah_ONanom4(a1,a2,a3,h,M)
c
c This subroutine calculates the vvvv-vertex with anomalous couplings
c and the on-shell polarization vectors for vector bosons
c
c IN:
c     polv1, polv2      : the incoming polvector wz(bosons in the CM frame)
c     polv3, polv4      : the outgoing polvector wz(bosons in the CM frame)
c OUT:
c     M(lambda1,lambda2) 
c
      IMPLICIT NONE

      INTEGER NGRAPHS
      PARAMETER( NGRAPHS=20 )
      DOUBLE PRECISION rZERO, rONE, rHalf
      PARAMETER( rZERO=0.0D0, rONE=1.0D0, rHALF = 0.5D0 )
      DOUBLE COMPLEX cIMAG, cZERO, cONE
      PARAMETER( cIMAG=(0.0D0,1.0D0), cZERO=(0.0D0,0.0D0), cONE=(1.0D0,0.0D0))
      DOUBLE COMPLEX a1(6,-1:1), a2(6,-1:1), a3(6,-1:1), h(6,-1:1)      
      DOUBLE COMPLEX M(-1:1,-1:1,-1:1,-1:1)
      
      DOUBLE COMPLEX OFFeps(6,-1:1,NGRAPHS), OFFepsC(6,-1:1,NGRAPHS)
      DOUBLE PRECISION K(0:3,4)      
      DOUBLE COMPLEX O(-1:1,4,-1:1,4), OP(4,-1:1,4)
      DOUBLE COMPLEXPP(1:4,1:4)
    
      INTEGER I,J,N,MU,NU
      INTEGER LAMBDA1, LAMBDA2, LAMBDA3, LAMBDA4, LAMBDA
      INTEGER L1,L2

      DOUBLE COMPLEX DOTCC, DOTRC,contract_Tjj
      DOUBLE PRECISION DOTRR
      EXTERNAL DOTCC, DOTRC, DOTRR, contract_Tjj
      double precision sw2, sw
#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "gmunu.inc"
#include "transm.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"
      sw2 = sin2w
      sw = sqrt(sin2w)
      OFFeps(6,-1:1,4) = cZero
      K(0:3,1:4) = rZero
      M(-1:1,-1:1,-1:1,-1:1) = cZero
      PP(1:4,1:4) = rZero
      O(-1:1,1:4,-1:1,1:4) = cZero
      OP(1:4,-1:1,1:4) = cZero  
      
      OFFeps(1:6,-1:1,1) = a1(1:6,-1:1)
      OFFeps(1:6,-1:1,2) = a2(1:6,-1:1)
      OFFeps(1:6,-1:1,3) = a3(1:6,-1:1)
      OFFeps(1:6,-1:1,4) = h(1:6,-1:1) 

      do i = 1,4
        K(0,i) = dble(OFFeps(5,-1,i))*(-rONE) !K1(0:3)
        K(1,i) = dble(OFFeps(6,-1,i))*(-rONE) !K2(0:3)
        K(2,i) = dimag(OFFeps(6,-1,i))*(-rONE)!K3(0:3)
        K(3,i) = dimag(OFFeps(5,-1,i))*(-rONE)!K4(0:3)
      enddo
      



      do i = 1,4
        do j = 1,4
          pp(i,j) = dotrr(k(0,i),k(0,j))
        enddo
      enddo

      do L1 = -1,1
        do L2 = -1,1 
          do i = 1,4
            do j = 1,4
              O(L1,i,L2,j) = dotcc(OFFeps(1,L1,i),OFFeps(1,L2,j))
            enddo
          enddo
        enddo
      enddo

      do LAMBDA = -1,1
        do i = 1,4 
          do j = 1,4
            OP(j,lambda,i) = dotrc(k(0,j),OFFeps(1,lambda,i))
          enddo
        enddo
      enddo


!!!! ===============================================================      
!!!! COUPLINGS CALCULATION
!!!! ===============================================================



      end



ccccccccccccccc
c aaaz-vertex c
ccccccccccccccc

      subroutine aaaz_ONanom4(a1,a2,a3,z,M)
c
c This subroutine calculates the vvvv-vertex with anomalous couplings
c and the on-shell polarization vectors for vector bosons
c
c IN:
c     polv1, polv2      : the incoming polvector wz(bosons in the CM frame)
c     polv3, polv4      : the outgoing polvector wz(bosons in the CM frame)
c OUT:
c     M(lambda1,lambda2) 
c
      IMPLICIT NONE

      INTEGER NGRAPHS
      PARAMETER( NGRAPHS=20 )
      DOUBLE PRECISION rZERO, rONE, rHalf
      PARAMETER( rZERO=0.0D0, rONE=1.0D0, rHALF = 0.5D0 )
      DOUBLE COMPLEX cIMAG, cZERO, cONE
      PARAMETER( cIMAG=(0.0D0,1.0D0), cZERO=(0.0D0,0.0D0), cONE=(1.0D0,0.0D0))
      DOUBLE COMPLEX a1(6,-1:1), a2(6,-1:1), a3(6,-1:1), z(6,-1:1)      
      DOUBLE COMPLEX M(-1:1,-1:1,-1:1,-1:1)
      
      DOUBLE COMPLEX OFFeps(6,-1:1,NGRAPHS), OFFepsC(6,-1:1,NGRAPHS)
      DOUBLE PRECISION K(0:3,4)      
      DOUBLE COMPLEX O(-1:1,4,-1:1,4), OP(4,-1:1,4)
      DOUBLE COMPLEXPP(1:4,1:4)
    
      INTEGER I,J,N,MU,NU
      INTEGER LAMBDA1, LAMBDA2, LAMBDA3, LAMBDA4, LAMBDA
      INTEGER L1,L2

      DOUBLE COMPLEX DOTCC, DOTRC,contract_Tjj
      DOUBLE PRECISION DOTRR
      EXTERNAL DOTCC, DOTRC, DOTRR, contract_Tjj
      double precision sw2, sw
#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "gmunu.inc"
#include "transm.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"
      sw2 = sin2w
      sw = sqrt(sin2w)
      OFFeps(6,-1:1,4) = cZero
      K(0:3,1:4) = rZero
      M(-1:1,-1:1,-1:1,-1:1) = cZero
      PP(1:4,1:4) = rZero
      O(-1:1,1:4,-1:1,1:4) = cZero
      OP(1:4,-1:1,1:4) = cZero  
      
      OFFeps(1:6,-1:1,1) = a1(1:6,-1:1)
      OFFeps(1:6,-1:1,2) = a2(1:6,-1:1)
      OFFeps(1:6,-1:1,3) = a3(1:6,-1:1)
      OFFeps(1:6,-1:1,4) = z(1:6,-1:1) 

      do i = 1,4
        K(0,i) = dble(OFFeps(5,-1,i))*(-rONE) !K1(0:3)
        K(1,i) = dble(OFFeps(6,-1,i))*(-rONE) !K2(0:3)
        K(2,i) = dimag(OFFeps(6,-1,i))*(-rONE)!K3(0:3)
        K(3,i) = dimag(OFFeps(5,-1,i))*(-rONE)!K4(0:3)
      enddo
      



      do i = 1,4
        do j = 1,4
          pp(i,j) = dotrr(k(0,i),k(0,j))
        enddo
      enddo

      do L1 = -1,1
        do L2 = -1,1 
          do i = 1,4
            do j = 1,4
              O(L1,i,L2,j) = dotcc(OFFeps(1,L1,i),OFFeps(1,L2,j))
            enddo
          enddo
        enddo
      enddo

      do LAMBDA = -1,1
        do i = 1,4 
          do j = 1,4
            OP(j,lambda,i) = dotrc(k(0,j),OFFeps(1,lambda,i))
          enddo
        enddo
      enddo


!!!! ===============================================================      
!!!! COUPLINGS CALCULATION
!!!! ===============================================================
                

         if( ft0.ne.0d0 ) then  !ft0 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(8*ft0*(gw**4)*(sw**3)*WMASS*(OP(1,lambda2,2)*OP(2,lambda1,1)*OP(3,lambda4,4)*OP(4,lambda3,3)-OP(3,lambda4,4)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda2,2)-OP(2,lambda4,4)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda3,3)-OP(2,lambda3,3)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*OP(4,lambda1,1)*(OP(2,lambda3,3)*OP(3,lambda2,2)-pp(2,3)*O(lambda2,2,lambda3,3))+pp(1,3)*pp(2,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*OP(3,lambda1,1)*(OP(2,lambda4,4)*OP(4,lambda2,2)-pp(2,4)*O(lambda2,2,lambda4,4))-OP(1,lambda2,2)*OP(2,lambda1,1)*pp(3,4)*O(lambda3,3,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)))/ZMASS

               enddo                      
              enddo
             enddo 
            enddo
         endif


         if( ft1.ne.0d0 ) then  !ft1 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(8*ft1*(gw**4)*(sw**3)*WMASS*(OP(1,lambda2,2)*OP(2,lambda1,1)*OP(3,lambda4,4)*OP(4,lambda3,3)-OP(3,lambda4,4)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda2,2)-OP(2,lambda4,4)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda3,3)-OP(2,lambda3,3)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*OP(4,lambda1,1)*(OP(2,lambda3,3)*OP(3,lambda2,2)-pp(2,3)*O(lambda2,2,lambda3,3))+pp(1,3)*pp(2,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*OP(3,lambda1,1)*(OP(2,lambda4,4)*OP(4,lambda2,2)-pp(2,4)*O(lambda2,2,lambda4,4))-OP(1,lambda2,2)*OP(2,lambda1,1)*pp(3,4)*O(lambda3,3,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)))/ZMASS

               enddo                      
              enddo
             enddo 
            enddo
         endif


         if( ft2.ne.0d0 ) then  !ft2 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(2*ft2*(gw**4)*(sw**3)*WMASS*(OP(1,lambda4,4)*OP(2,lambda3,3)*OP(3,lambda1,1)*OP(4,lambda2,2)+OP(1,lambda4,4)*OP(2,lambda1,1)*OP(3,lambda2,2)*OP(4,lambda3,3)-OP(2,lambda4,4)*OP(4,lambda3,3)*pp(1,3)*O(lambda1,1,lambda2,2)-OP(2,lambda3,3)*OP(3,lambda4,4)*pp(1,4)*O(lambda1,1,lambda2,2)-OP(1,lambda4,4)*OP(4,lambda3,3)*pp(2,3)*O(lambda1,1,lambda2,2)+OP(1,lambda4,4)*OP(2,lambda3,3)*pp(3,4)*O(lambda1,1,lambda2,2)-OP(3,lambda4,4)*OP(4,lambda2,2)*pp(1,2)*O(lambda1,1,lambda3,3)-OP(2,lambda4,4)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda3,3)-OP(1,lambda4,4)*OP(4,lambda2,2)*pp(2,3)*O(lambda1,1,lambda3,3)+OP(1,lambda4,4)*OP(3,lambda2,2)*pp(2,4)*O(lambda1,1,lambda3,3)-OP(3,lambda2,2)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda4,4)-OP(2,lambda3,3)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda4,4)-OP(3,lambda4,4)*OP(4,lambda1,1)*pp(1,2)*O(lambda2,2,lambda3,3)-OP(2,lambda4,4)*OP(4,lambda1,1)*pp(1,3)*O(lambda2,2,lambda3,3)+OP(2,lambda4,4)*OP(3,lambda1,1)*pp(1,4)*O(lambda2,2,lambda3,3)+OP(2,lambda1,1)*OP(3,lambda4,4)*pp(1,4)*O(lambda2,2,lambda3,3)-OP(1,lambda4,4)*OP(3,lambda1,1)*pp(2,4)*O(lambda2,2,lambda3,3)-OP(1,lambda4,4)*OP(2,lambda1,1)*pp(3,4)*O(lambda2,2,lambda3,3)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-OP(3,lambda1,1)*OP(4,lambda3,3)*pp(1,2)*O(lambda2,2,lambda4,4)+OP(2,lambda3,3)*OP(4,lambda1,1)*pp(1,3)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*OP(4,lambda3,3)*pp(1,3)*O(lambda2,2,lambda4,4)-OP(2,lambda3,3)*OP(3,lambda1,1)*pp(1,4)*O(lambda2,2,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(-(OP(3,lambda4,4)*pp(2,4)*O(lambda1,1,lambda2,2))+OP(2,lambda4,4)*(OP(3,lambda2,2)*OP(4,lambda1,1)+pp(3,4)*O(lambda1,1,lambda2,2))+OP(4,lambda2,2)*pp(2,3)*O(lambda1,1,lambda4,4)-OP(3,lambda2,2)*pp(2,4)*O(lambda1,1,lambda4,4)-OP(4,lambda1,1)*pp(2,3)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*(OP(3,lambda4,4)*OP(4,lambda2,2)-pp(3,4)*O(lambda2,2,lambda4,4)))+OP(3,lambda2,2)*OP(4,lambda1,1)*pp(1,2)*O(lambda3,3,lambda4,4)+OP(3,lambda1,1)*OP(4,lambda2,2)*pp(1,2)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(4,lambda2,2)*pp(1,3)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(3,lambda2,2)*pp(1,4)*O(lambda3,3,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+OP(1,lambda2,2)*(OP(3,lambda4,4)*pp(2,4)*O(lambda1,1,lambda3,3)+OP(2,lambda4,4)*(OP(3,lambda1,1)*OP(4,lambda3,3)-pp(3,4)*O(lambda1,1,lambda3,3))+OP(4,lambda3,3)*pp(2,3)*O(lambda1,1,lambda4,4)+OP(2,lambda3,3)*(OP(3,lambda4,4)*OP(4,lambda1,1)-pp(3,4)*O(lambda1,1,lambda4,4))-OP(4,lambda1,1)*pp(2,3)*O(lambda3,3,lambda4,4)-OP(3,lambda1,1)*pp(2,4)*O(lambda3,3,lambda4,4))))/ZMASS

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( ft3.ne.0d0 ) then  !ft3 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(1/ZMASS)*2*ft3*WMASS*gw**4*sw**3*(-pp(1,2)*OP(4,lambda2,2)*OP(3,lambda4,4)*O(lambda1,1,lambda3,3)-pp(2,3)*OP(4,lambda2,2)*OP(1,lambda4,4)*O(lambda1,1,lambda3,3)-pp(1,3)*OP(4,lambda2,2)*OP(2,lambda3,3)*O(lambda1,1,lambda4,4)+pp(1,2)*OP(3,lambda1,1)*OP(4,lambda2,2)*O(lambda3,3,lambda4,4)-pp(1,3)*OP(2,lambda1,1)*OP(4,lambda2,2)*O(lambda3,3,lambda4,4)-pp(1,3)*OP(4,lambda3,3)*OP(2,lambda4,4)*O(lambda1,1,lambda2,2)-pp(1,4)*OP(2,lambda3,3)*OP(3,lambda4,4)*O(lambda1,1,lambda2,2)-pp(2,3)*OP(4,lambda3,3)*OP(1,lambda4,4)*O(lambda1,1,lambda2,2)+pp(3,4)*OP(2,lambda3,3)*OP(1,lambda4,4)*O(lambda1,1,lambda2,2)-pp(1,4)*OP(3,lambda2,2)*OP(2,lambda4,4)*O(lambda1,1,lambda3,3)+pp(2,4)*OP(3,lambda2,2)*OP(1,lambda4,4)*O(lambda1,1,lambda3,3)-pp(1,2)*OP(3,lambda2,2)*OP(4,lambda3,3)*O(lambda1,1,lambda4,4)-pp(1,2)*OP(4,lambda1,1)*OP(3,lambda4,4)*O(lambda2,2,lambda3,3)-pp(1,3)*OP(4,lambda1,1)*OP(2,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,4)*OP(3,lambda1,1)*OP(2,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,4)*OP(2,lambda1,1)*OP(3,lambda4,4)*O(lambda2,2,lambda3,3)-pp(2,4)*OP(3,lambda1,1)*OP(1,lambda4,4)*O(lambda2,2,lambda3,3)-pp(3,4)*OP(2,lambda1,1)*OP(1,lambda4,4)*O(lambda2,2,lambda3,3)-pp(1,2)*OP(3,lambda1,1)*OP(4,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,3)*OP(4,lambda1,1)*OP(2,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,3)*OP(2,lambda1,1)*OP(4,lambda3,3)*O(lambda2,2,lambda4,4)-pp(1,4)*OP(3,lambda1,1)*OP(2,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(pp(2,4)*(-OP(3,lambda4,4))*O(lambda1,1,lambda2,2)+OP(2,lambda4,4)*(OP(4,lambda1,1)*OP(3,lambda2,2)+pp(3,4)*O(lambda1,1,lambda2,2))+pp(2,3)*OP(4,lambda2,2)*O(lambda1,1,lambda4,4)-pp(2,4)*OP(3,lambda2,2)*O(lambda1,1,lambda4,4)-pp(2,3)*OP(4,lambda1,1)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*(OP(4,lambda2,2)*OP(3,lambda4,4)-pp(3,4)*O(lambda2,2,lambda4,4)))+pp(1,2)*OP(4,lambda1,1)*OP(3,lambda2,2)*O(lambda3,3,lambda4,4)-pp(1,4)*OP(2,lambda1,1)*OP(3,lambda2,2)*O(lambda3,3,lambda4,4)+OP(1,lambda2,2)*(pp(2,4)*OP(3,lambda4,4)*O(lambda1,1,lambda3,3)+OP(2,lambda4,4)*(OP(3,lambda1,1)*OP(4,lambda3,3)-pp(3,4)*O(lambda1,1,lambda3,3))+pp(2,3)*OP(4,lambda3,3)*O(lambda1,1,lambda4,4)+OP(2,lambda3,3)*(OP(4,lambda1,1)*OP(3,lambda4,4)-pp(3,4)*O(lambda1,1,lambda4,4))-pp(2,3)*OP(4,lambda1,1)*O(lambda3,3,lambda4,4)-pp(2,4)*OP(3,lambda1,1)*O(lambda3,3,lambda4,4))+OP(3,lambda1,1)*OP(4,lambda2,2)*OP(2,lambda3,3)*OP(1,lambda4,4)+OP(2,lambda1,1)*OP(3,lambda2,2)*OP(4,lambda3,3)*OP(1,lambda4,4)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4))
               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( ft4.ne.0d0 ) then  !ft4 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) -(gw**4*ft4*sw**3*(-1 + 2*SW2)*ZMASS*(OP(1, lambda4, 4)*OP(2, lambda3, 3)*OP(3, lambda1, 1)*OP(4, lambda2, 2) + OP(1, lambda4, 4)*OP(2, lambda1, 1)*OP(3, lambda2, 2)*OP(4, lambda3, 3) + pp(3, 4)*OP(1, lambda4, 4)*OP(2, lambda3, 3)*O(lambda1, 1, lambda2, 2) - pp(1, 4)*OP(2, lambda3, 3)*OP(3, lambda4, 4)*O(lambda1, 1, lambda2, 2) - pp(2, 3)*OP(1, lambda4, 4)*OP(4, lambda3, 3)*O(lambda1, 1, lambda2, 2) - pp(1, 3)*OP(2, lambda4, 4)*OP(4, lambda3, 3)*O(lambda1, 1, lambda2, 2) + pp(2, 4)*OP(1, lambda4, 4)*OP(3, lambda2, 2)*O(lambda1, 1, lambda3, 3) - pp(1, 4)*OP(2, lambda4, 4)*OP(3, lambda2, 2)*O(lambda1, 1, lambda3, 3) - pp(2, 3)*OP(1, lambda4, 4)*OP(4, lambda2, 2)*O(lambda1, 1, lambda3, 3) - pp(1, 2)*OP(3, lambda4, 4)*OP(4, lambda2, 2)*O(lambda1, 1, lambda3, 3) - pp(1, 3)*OP(2, lambda3, 3)*OP(4, lambda2, 2)*O(lambda1, 1, lambda4, 4) - pp(1, 2)*OP(3, lambda2, 2)*OP(4, lambda3, 3)*O(lambda1, 1, lambda4, 4) - pp(3, 4)*OP(1, lambda4, 4)*OP(2, lambda1, 1)*O(lambda2, 2, lambda3, 3) - pp(2, 4)*OP(1, lambda4, 4)*OP(3, lambda1, 1)*O(lambda2, 2, lambda3, 3) + pp(1, 4)*OP(2, lambda4, 4)*OP(3, lambda1, 1)*O(lambda2, 2, lambda3, 3) + pp(1, 4)*OP(2, lambda1, 1)*OP(3, lambda4, 4)*O(lambda2, 2, lambda3, 3) - pp(1, 3)*OP(2, lambda4, 4)*OP(4, lambda1, 1)*O(lambda2, 2, lambda3, 3) - pp(1, 2)*OP(3, lambda4, 4)*OP(4, lambda1, 1)*O(lambda2, 2, lambda3, 3) + pp(1, 3)*pp(2, 4)*O(lambda1, 1, lambda4, 4)*O(lambda2, 2, lambda3, 3) + pp(1, 2)*pp(3, 4)*O(lambda1, 1, lambda4, 4)*O(lambda2, 2, lambda3, 3) - pp(1, 4)*OP(2, lambda3, 3)*OP(3, lambda1, 1)*O(lambda2, 2, lambda4, 4) + pp(1, 3)*OP(2, lambda3, 3)*OP(4, lambda1, 1)*O(lambda2, 2, lambda4, 4) + pp(1, 3)*OP(2, lambda1, 1)*OP(4, lambda3, 3)*O(lambda2, 2, lambda4, 4) - pp(1, 2)*OP(3, lambda1, 1)*OP(4, lambda3, 3)*O(lambda2, 2, lambda4, 4) + pp(1, 4)*pp(2, 3)*O(lambda1, 1, lambda3, 3)*O(lambda2, 2, lambda4, 4) + pp(1, 2)*pp(3, 4)*O(lambda1, 1, lambda3, 3)*O(lambda2, 2, lambda4, 4) + OP(1, lambda3, 3)*(-(pp(2, 4)*OP(3, lambda4, 4)*O(lambda1, 1, lambda2, 2)) + OP(2, lambda4, 4)*(OP(3, lambda2, 2)*OP(4, lambda1, 1) + pp(3, 4)*O(lambda1, 1, lambda2, 2)) -   pp(2, 4)*OP(3, lambda2, 2)*O(lambda1, 1, lambda4, 4) + pp(2, 3)*OP(4, lambda2, 2)*O(lambda1, 1, lambda4, 4) - pp(2, 3)*OP(4, lambda1, 1)*O(lambda2, 2, lambda4, 4) +   OP(2, lambda1, 1)*(OP(3, lambda4, 4)*OP(4, lambda2, 2) - pp(3, 4)*O(lambda2, 2, lambda4, 4))) - pp(1, 4)*OP(2, lambda1, 1)*OP(3, lambda2, 2)*O(lambda3, 3, lambda4, 4) + pp(1, 2)*OP(3, lambda2, 2)*OP(4, lambda1, 1)*O(lambda3, 3, lambda4, 4) - pp(1, 3)*OP(2, lambda1, 1)*OP(4, lambda2, 2)*O(lambda3, 3, lambda4, 4) + pp(1, 2)*OP(3, lambda1, 1)*OP(4, lambda2, 2)*O(lambda3, 3, lambda4, 4) + pp(1, 4)*pp(2, 3)*O(lambda1, 1, lambda2, 2)*O(lambda3, 3, lambda4, 4) + pp(1, 3)*pp(2, 4)*O(lambda1, 1, lambda2, 2)*O(lambda3, 3, lambda4, 4) + OP(1, lambda2, 2)*(pp(2, 4)*OP(3, lambda4, 4)*O(lambda1, 1, lambda3, 3) +   OP(2, lambda4, 4)*(OP(3, lambda1, 1)*OP(4, lambda3, 3) - pp(3, 4)*O(lambda1, 1, lambda3, 3)) + pp(2, 3)*OP(4, lambda3, 3)*O(lambda1, 1, lambda4, 4) +   OP(2, lambda3, 3)*(OP(3, lambda4, 4)*OP(4, lambda1, 1) - pp(3, 4)*O(lambda1, 1, lambda4, 4)) - pp(2, 4)*OP(3, lambda1, 1)*O(lambda3, 3, lambda4, 4) -   pp(2, 3)*OP(4, lambda1, 1)*O(lambda3, 3, lambda4, 4))))/(2*WMASS)

               enddo                      
              enddo
             enddo 
            enddo
         endif
         if( ft5.ne.0d0 ) then  !ft5 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(-2*ft5*(gw**4)*(sw**3)*(-1+2*SW2)*ZMASS*(OP(1,lambda2,2)*OP(2,lambda1,1)*OP(3,lambda4,4)*OP(4,lambda3,3)-OP(3,lambda4,4)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda2,2)-OP(2,lambda4,4)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda3,3)-OP(2,lambda3,3)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*OP(4,lambda1,1)*(OP(2,lambda3,3)*OP(3,lambda2,2)-pp(2,3)*O(lambda2,2,lambda3,3))+pp(1,3)*pp(2,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*OP(3,lambda1,1)*(OP(2,lambda4,4)*OP(4,lambda2,2)-pp(2,4)*O(lambda2,2,lambda4,4))-OP(1,lambda2,2)*OP(2,lambda1,1)*pp(3,4)*O(lambda3,3,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)))/WMASS

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( ft6.ne.0d0 ) then  !ft6 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(-2*ft6*(gw**4)*(sw**3)*(-1+2*SW2)*ZMASS*(OP(1,lambda2,2)*OP(2,lambda1,1)*OP(3,lambda4,4)*OP(4,lambda3,3)-OP(3,lambda4,4)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda2,2)-OP(2,lambda4,4)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda3,3)-OP(2,lambda3,3)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*OP(4,lambda1,1)*(OP(2,lambda3,3)*OP(3,lambda2,2)-pp(2,3)*O(lambda2,2,lambda3,3))+pp(1,3)*pp(2,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*OP(3,lambda1,1)*(OP(2,lambda4,4)*OP(4,lambda2,2)-pp(2,4)*O(lambda2,2,lambda4,4))-OP(1,lambda2,2)*OP(2,lambda1,1)*pp(3,4)*O(lambda3,3,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)))/WMASS

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( ft7.ne.0d0 ) then  !ft7 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) -(ft7*(gw**4)*(sw**3)*(-1+2*SW2)*ZMASS*(OP(1,lambda4,4)*OP(2,lambda3,3)*OP(3,lambda1,1)*OP(4,lambda2,2)+OP(1,lambda4,4)*OP(2,lambda1,1)*OP(3,lambda2,2)*OP(4,lambda3,3)-OP(2,lambda4,4)*OP(4,lambda3,3)*pp(1,3)*O(lambda1,1,lambda2,2)-OP(2,lambda3,3)*OP(3,lambda4,4)*pp(1,4)*O(lambda1,1,lambda2,2)-OP(1,lambda4,4)*OP(4,lambda3,3)*pp(2,3)*O(lambda1,1,lambda2,2)+OP(1,lambda4,4)*OP(2,lambda3,3)*pp(3,4)*O(lambda1,1,lambda2,2)-OP(3,lambda4,4)*OP(4,lambda2,2)*pp(1,2)*O(lambda1,1,lambda3,3)-OP(2,lambda4,4)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda3,3)-OP(1,lambda4,4)*OP(4,lambda2,2)*pp(2,3)*O(lambda1,1,lambda3,3)+OP(1,lambda4,4)*OP(3,lambda2,2)*pp(2,4)*O(lambda1,1,lambda3,3)-OP(3,lambda2,2)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda4,4)-OP(2,lambda3,3)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda4,4)-OP(3,lambda4,4)*OP(4,lambda1,1)*pp(1,2)*O(lambda2,2,lambda3,3)-OP(2,lambda4,4)*OP(4,lambda1,1)*pp(1,3)*O(lambda2,2,lambda3,3)+OP(2,lambda4,4)*OP(3,lambda1,1)*pp(1,4)*O(lambda2,2,lambda3,3)+OP(2,lambda1,1)*OP(3,lambda4,4)*pp(1,4)*O(lambda2,2,lambda3,3)-OP(1,lambda4,4)*OP(3,lambda1,1)*pp(2,4)*O(lambda2,2,lambda3,3)-OP(1,lambda4,4)*OP(2,lambda1,1)*pp(3,4)*O(lambda2,2,lambda3,3)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-OP(3,lambda1,1)*OP(4,lambda3,3)*pp(1,2)*O(lambda2,2,lambda4,4)+OP(2,lambda3,3)*OP(4,lambda1,1)*pp(1,3)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*OP(4,lambda3,3)*pp(1,3)*O(lambda2,2,lambda4,4)-OP(2,lambda3,3)*OP(3,lambda1,1)*pp(1,4)*O(lambda2,2,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(-(OP(3,lambda4,4)*pp(2,4)*O(lambda1,1,lambda2,2))+OP(2,lambda4,4)*(OP(3,lambda2,2)*OP(4,lambda1,1)+pp(3,4)*O(lambda1,1,lambda2,2))+OP(4,lambda2,2)*pp(2,3)*O(lambda1,1,lambda4,4)-OP(3,lambda2,2)*pp(2,4)*O(lambda1,1,lambda4,4)-OP(4,lambda1,1)*pp(2,3)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*(OP(3,lambda4,4)*OP(4,lambda2,2)-pp(3,4)*O(lambda2,2,lambda4,4)))+OP(3,lambda2,2)*OP(4,lambda1,1)*pp(1,2)*O(lambda3,3,lambda4,4)+OP(3,lambda1,1)*OP(4,lambda2,2)*pp(1,2)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(4,lambda2,2)*pp(1,3)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(3,lambda2,2)*pp(1,4)*O(lambda3,3,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+OP(1,lambda2,2)*(OP(3,lambda4,4)*pp(2,4)*O(lambda1,1,lambda3,3)+OP(2,lambda4,4)*(OP(3,lambda1,1)*OP(4,lambda3,3)-pp(3,4)*O(lambda1,1,lambda3,3))+OP(4,lambda3,3)*pp(2,3)*O(lambda1,1,lambda4,4)+OP(2,lambda3,3)*(OP(3,lambda4,4)*OP(4,lambda1,1)-pp(3,4)*O(lambda1,1,lambda4,4))-OP(4,lambda1,1)*pp(2,3)*O(lambda3,3,lambda4,4)-OP(3,lambda1,1)*pp(2,4)*O(lambda3,3,lambda4,4))))/(2*WMASS)

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( ft8.ne.0d0 ) then  !ft8 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(2*ft8*(gw**4)*(sw**5)*(-1+SW2)*(ZMASS**3)*(OP(1,lambda2,2)*OP(2,lambda1,1)*OP(3,lambda4,4)*OP(4,lambda3,3)-OP(3,lambda4,4)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda2,2)-OP(2,lambda4,4)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda3,3)-OP(2,lambda3,3)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*OP(4,lambda1,1)*(OP(2,lambda3,3)*OP(3,lambda2,2)-pp(2,3)*O(lambda2,2,lambda3,3))+pp(1,3)*pp(2,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*OP(3,lambda1,1)*(OP(2,lambda4,4)*OP(4,lambda2,2)-pp(2,4)*O(lambda2,2,lambda4,4))-OP(1,lambda2,2)*OP(2,lambda1,1)*pp(3,4)*O(lambda3,3,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)))/(WMASS**3)

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( ft9.ne.0d0 ) then  !ft9 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(ft9*(gw**4)*(sw**5)*(-1+SW2)*(ZMASS**3)*(OP(1,lambda4,4)*OP(2,lambda3,3)*OP(3,lambda1,1)*OP(4,lambda2,2)+OP(1,lambda4,4)*OP(2,lambda1,1)*OP(3,lambda2,2)*OP(4,lambda3,3)-OP(2,lambda4,4)*OP(4,lambda3,3)*pp(1,3)*O(lambda1,1,lambda2,2)-OP(2,lambda3,3)*OP(3,lambda4,4)*pp(1,4)*O(lambda1,1,lambda2,2)-OP(1,lambda4,4)*OP(4,lambda3,3)*pp(2,3)*O(lambda1,1,lambda2,2)+OP(1,lambda4,4)*OP(2,lambda3,3)*pp(3,4)*O(lambda1,1,lambda2,2)-OP(3,lambda4,4)*OP(4,lambda2,2)*pp(1,2)*O(lambda1,1,lambda3,3)-OP(2,lambda4,4)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda3,3)-OP(1,lambda4,4)*OP(4,lambda2,2)*pp(2,3)*O(lambda1,1,lambda3,3)+OP(1,lambda4,4)*OP(3,lambda2,2)*pp(2,4)*O(lambda1,1,lambda3,3)-OP(3,lambda2,2)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda4,4)-OP(2,lambda3,3)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda4,4)-OP(3,lambda4,4)*OP(4,lambda1,1)*pp(1,2)*O(lambda2,2,lambda3,3)-OP(2,lambda4,4)*OP(4,lambda1,1)*pp(1,3)*O(lambda2,2,lambda3,3)+OP(2,lambda4,4)*OP(3,lambda1,1)*pp(1,4)*O(lambda2,2,lambda3,3)+OP(2,lambda1,1)*OP(3,lambda4,4)*pp(1,4)*O(lambda2,2,lambda3,3)-OP(1,lambda4,4)*OP(3,lambda1,1)*pp(2,4)*O(lambda2,2,lambda3,3)-OP(1,lambda4,4)*OP(2,lambda1,1)*pp(3,4)*O(lambda2,2,lambda3,3)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-OP(3,lambda1,1)*OP(4,lambda3,3)*pp(1,2)*O(lambda2,2,lambda4,4)+OP(2,lambda3,3)*OP(4,lambda1,1)*pp(1,3)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*OP(4,lambda3,3)*pp(1,3)*O(lambda2,2,lambda4,4)-OP(2,lambda3,3)*OP(3,lambda1,1)*pp(1,4)*O(lambda2,2,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(-(OP(3,lambda4,4)*pp(2,4)*O(lambda1,1,lambda2,2))+OP(2,lambda4,4)*(OP(3,lambda2,2)*OP(4,lambda1,1)+pp(3,4)*O(lambda1,1,lambda2,2))+OP(4,lambda2,2)*pp(2,3)*O(lambda1,1,lambda4,4)-OP(3,lambda2,2)*pp(2,4)*O(lambda1,1,lambda4,4)-OP(4,lambda1,1)*pp(2,3)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*(OP(3,lambda4,4)*OP(4,lambda2,2)-pp(3,4)*O(lambda2,2,lambda4,4)))+OP(3,lambda2,2)*OP(4,lambda1,1)*pp(1,2)*O(lambda3,3,lambda4,4)+OP(3,lambda1,1)*OP(4,lambda2,2)*pp(1,2)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(4,lambda2,2)*pp(1,3)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(3,lambda2,2)*pp(1,4)*O(lambda3,3,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+OP(1,lambda2,2)*(OP(3,lambda4,4)*pp(2,4)*O(lambda1,1,lambda3,3)+OP(2,lambda4,4)*(OP(3,lambda1,1)*OP(4,lambda3,3)-pp(3,4)*O(lambda1,1,lambda3,3))+OP(4,lambda3,3)*pp(2,3)*O(lambda1,1,lambda4,4)+OP(2,lambda3,3)*(OP(3,lambda4,4)*OP(4,lambda1,1)-pp(3,4)*O(lambda1,1,lambda4,4))-OP(4,lambda1,1)*pp(2,3)*O(lambda3,3,lambda4,4)-OP(3,lambda1,1)*pp(2,4)*O(lambda3,3,lambda4,4))))/(2*(WMASS**3))

               enddo                      
              enddo
             enddo 
            enddo
         endif

      

      end



ccccccccccccccc
c aahh-vertex c
ccccccccccccccc

      subroutine aahh_ONanom4(a1,a2,h1,h2,M)
c
c This subroutine calculates the vvvv-vertex with anomalous couplings
c and the on-shell polarization vectors for vector bosons
c
c IN:
c     polv1, polv2      : the incoming polvector wz(bosons in the CM frame)
c     polv3, polv4      : the outgoing polvector wz(bosons in the CM frame)
c OUT:
c     M(lambda1,lambda2) 
c
      IMPLICIT NONE

      INTEGER NGRAPHS
      PARAMETER( NGRAPHS=20 )
      DOUBLE PRECISION rZERO, rONE, rHalf
      PARAMETER( rZERO=0.0D0, rONE=1.0D0, rHALF = 0.5D0 )
      DOUBLE COMPLEX cIMAG, cZERO, cONE
      PARAMETER( cIMAG=(0.0D0,1.0D0), cZERO=(0.0D0,0.0D0), cONE=(1.0D0,0.0D0))
      DOUBLE COMPLEX a1(6,-1:1), a2(6,-1:1), h1(6,-1:1), h2(6,-1:1)      
      DOUBLE COMPLEX M(-1:1,-1:1,-1:1,-1:1)
      
      DOUBLE COMPLEX OFFeps(6,-1:1,NGRAPHS), OFFepsC(6,-1:1,NGRAPHS)
      DOUBLE PRECISION K(0:3,4)      
      DOUBLE COMPLEX O(-1:1,4,-1:1,4), OP(4,-1:1,4)
      DOUBLE COMPLEXPP(1:4,1:4)
    
      INTEGER I,J,N,MU,NU
      INTEGER LAMBDA1, LAMBDA2, LAMBDA3, LAMBDA4, LAMBDA
      INTEGER L1,L2

      DOUBLE COMPLEX DOTCC, DOTRC,contract_Tjj
      DOUBLE PRECISION DOTRR
      EXTERNAL DOTCC, DOTRC, DOTRR, contract_Tjj
      double precision sw2, sw
#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "gmunu.inc"
#include "transm.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"
      sw2 = sin2w
      sw = sqrt(sin2w)
      OFFeps(6,-1:1,4) = cZero
      K(0:3,1:4) = rZero
      M(-1:1,-1:1,-1:1,-1:1) = cZero
      PP(1:4,1:4) = rZero
      O(-1:1,1:4,-1:1,1:4) = cZero
      OP(1:4,-1:1,1:4) = cZero  
      
      OFFeps(1:6,-1:1,1) = a1(1:6,-1:1)
      OFFeps(1:6,-1:1,2) = a2(1:6,-1:1)
      OFFeps(1:6,-1:1,3) = h1(1:6,-1:1)
      OFFeps(1:6,-1:1,4) = h2(1:6,-1:1) 

      do i = 1,4
        K(0,i) = dble(OFFeps(5,-1,i))*(-rONE) !K1(0:3)
        K(1,i) = dble(OFFeps(6,-1,i))*(-rONE) !K2(0:3)
        K(2,i) = dimag(OFFeps(6,-1,i))*(-rONE)!K3(0:3)
        K(3,i) = dimag(OFFeps(5,-1,i))*(-rONE)!K4(0:3)
      enddo
      



      do i = 1,4
        do j = 1,4
          pp(i,j) = dotrr(k(0,i),k(0,j))
        enddo
      enddo

      do L1 = -1,1
        do L2 = -1,1 
          do i = 1,4
            do j = 1,4
              O(L1,i,L2,j) = dotcc(OFFeps(1,L1,i),OFFeps(1,L2,j))
            enddo
          enddo
        enddo
      enddo

      do LAMBDA = -1,1
        do i = 1,4 
          do j = 1,4
            OP(j,lambda,i) = dotrc(k(0,j),OFFeps(1,lambda,i))
          enddo
        enddo
      enddo


!!!! ===============================================================      
!!!! COUPLINGS CALCULATION
!!!! ===============================================================





         if( fm0.ne.0d0 ) then  !fm0 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
!              do lambda3 = -1,1
!               do lambda4 = -1,1
                  M(lambda1,lambda2,0,0) = M(lambda1,lambda2,0,0) +2*fm0*(gw**2)*(sw**2)*1.0*1.0*pp(3,4)*(OP(1,lambda2,2)*OP(2,lambda1,1)-pp(1,2)*O(lambda1,1,lambda2,2))

!               enddo                      
!              enddo
             enddo 
            enddo
         endif


         if( fm1.ne.0d0 ) then  !fm1 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
!              do lambda3 = -1,1
!               do lambda4 = -1,1
                  M(lambda1,lambda2,0,0) = M(lambda1,lambda2,0,0) +(fm1*(gw**2)*(sw**2)*1.0*1.0*(-(OP(2,lambda1,1)*OP(4,lambda2,2)*pp(1,3))+OP(3,lambda2,2)*(OP(4,lambda1,1)*pp(1,2)-OP(2,lambda1,1)*pp(1,4))-OP(1,lambda2,2)*OP(4,lambda1,1)*pp(2,3)+OP(3,lambda1,1)*(OP(4,lambda2,2)*pp(1,2)-OP(1,lambda2,2)*pp(2,4))+pp(1,4)*pp(2,3)*O(lambda1,1,lambda2,2)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda2,2)))/2

!               enddo                      
!              enddo
             enddo 
            enddo
         endif


         if( fm2.ne.0d0 ) then  !fm2 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
!              do lambda3 = -1,1
!               do lambda4 = -1,1
                  M(lambda1,lambda2,0,0) = M(lambda1,lambda2,0,0) +fm2*(gw**2)*(sw**2)*1.0*1.0*pp(3,4)*(OP(1,lambda2,2)*OP(2,lambda1,1)-pp(1,2)*O(lambda1,1,lambda2,2))

!               enddo                      
!              enddo
             enddo 
            enddo
         endif

         if( fm3.ne.0d0 ) then  !fm3 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
!              do lambda3 = -1,1
!               do lambda4 = -1,1
                  M(lambda1,lambda2,0,0) = M(lambda1,lambda2,0,0) +(fm3*(gw**2)*(sw**2)*(-1+SW2)*(ZMASS**2)*1.0*1.0*(OP(2,lambda1,1)*OP(4,lambda2,2)*pp(1,3)+OP(3,lambda2,2)*(-(OP(4,lambda1,1)*pp(1,2))+OP(2,lambda1,1)*pp(1,4))+OP(1,lambda2,2)*OP(4,lambda1,1)*pp(2,3)+OP(3,lambda1,1)*(-(OP(4,lambda2,2)*pp(1,2))+OP(1,lambda2,2)*pp(2,4))-pp(1,4)*pp(2,3)*O(lambda1,1,lambda2,2)-pp(1,3)*pp(2,4)*O(lambda1,1,lambda2,2)))/(4*(WMASS**2))

!               enddo                      
!              enddo
             enddo 
            enddo
         endif

         if( fm4.ne.0d0 ) then  !fm4 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
!              do lambda3 = -1,1
!               do lambda4 = -1,1
                  M(lambda1,lambda2,0,0) = M(lambda1,lambda2,0,0) -(fm4*(gw**2)*(sw**2)*1.0*1.0*pp(3,4)*(OP(1,lambda2,2)*OP(2,lambda1,1)-pp(1,2)*O(lambda1,1,lambda2,2)))

!               enddo                      
!              enddo
             enddo 
            enddo
         endif

         if( fm5.ne.0d0 ) then  !fm5 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
!              do lambda3 = -1,1
!               do lambda4 = -1,1
                  M(lambda1,lambda2,0,0) = M(lambda1,lambda2,0,0) +(fm5*(gw**2)*(sw**2)*1.0*1.0*(-(OP(2,lambda1,1)*OP(4,lambda2,2)*pp(1,3))+OP(3,lambda2,2)*(OP(4,lambda1,1)*pp(1,2)-OP(2,lambda1,1)*pp(1,4))-OP(1,lambda2,2)*OP(4,lambda1,1)*pp(2,3)+OP(3,lambda1,1)*(OP(4,lambda2,2)*pp(1,2)-OP(1,lambda2,2)*pp(2,4))+pp(1,4)*pp(2,3)*O(lambda1,1,lambda2,2)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda2,2)))/4

!               enddo                      
!              enddo
             enddo 
            enddo
         endif

         if( fm5hc.ne.0d0 ) then  !fm5hc coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
!              do lambda3 = -1,1
!               do lambda4 = -1,1
                  M(lambda1,lambda2,0,0) = M(lambda1,lambda2,0,0) +(fm5hc*(gw**2)*(sw**2)*1.0*1.0*(-(OP(2,lambda1,1)*OP(4,lambda2,2)*pp(1,3))+OP(3,lambda2,2)*(OP(4,lambda1,1)*pp(1,2)-OP(2,lambda1,1)*pp(1,4))-OP(1,lambda2,2)*OP(4,lambda1,1)*pp(2,3)+OP(3,lambda1,1)*(OP(4,lambda2,2)*pp(1,2)-OP(1,lambda2,2)*pp(2,4))+pp(1,4)*pp(2,3)*O(lambda1,1,lambda2,2)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda2,2)))/4

!               enddo                      
!              enddo
             enddo 
            enddo
         endif

         if( fm6.ne.0d0 ) then  !fm6 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
!              do lambda3 = -1,1
!               do lambda4 = -1,1
                  M(lambda1,lambda2,0,0) = M(lambda1,lambda2,0,0) +fm6*(gw**2)*(sw**2)*1.0*1.0*pp(3,4)*(OP(1,lambda2,2)*OP(2,lambda1,1)-pp(1,2)*O(lambda1,1,lambda2,2))

!               enddo                      
!              enddo
             enddo 
            enddo
         endif

         if( fm7.ne.0d0 ) then  !fm7 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
!              do lambda3 = -1,1
!               do lambda4 = -1,1
                  M(lambda1,lambda2,0,0) = M(lambda1,lambda2,0,0) +(fm7*(gw**2)*(sw**2)*1.0*1.0*(OP(2,lambda1,1)*OP(4,lambda2,2)*pp(1,3)+OP(3,lambda2,2)*(-(OP(4,lambda1,1)*pp(1,2))+OP(2,lambda1,1)*pp(1,4))+OP(1,lambda2,2)*OP(4,lambda1,1)*pp(2,3)+OP(3,lambda1,1)*(-(OP(4,lambda2,2)*pp(1,2))+OP(1,lambda2,2)*pp(2,4))-pp(1,4)*pp(2,3)*O(lambda1,1,lambda2,2)-pp(1,3)*pp(2,4)*O(lambda1,1,lambda2,2)))/4

!               enddo                      
!              enddo
             enddo 
            enddo
         endif         

      end



ccccccccccccccc
c aahz-vertex c
ccccccccccccccc

      subroutine aahz_ONanom4(a1,a2,h,z,M)
c
c This subroutine calculates the vvvv-vertex with anomalous couplings
c and the on-shell polarization vectors for vector bosons
c
c IN:
c     polv1, polv2      : the incoming polvector wz(bosons in the CM frame)
c     polv3, polv4      : the outgoing polvector wz(bosons in the CM frame)
c OUT:
c     M(lambda1,lambda2) 
c
      IMPLICIT NONE

      INTEGER NGRAPHS
      PARAMETER( NGRAPHS=20 )
      DOUBLE PRECISION rZERO, rONE, rHalf
      PARAMETER( rZERO=0.0D0, rONE=1.0D0, rHALF = 0.5D0 )
      DOUBLE COMPLEX cIMAG, cZERO, cONE
      PARAMETER( cIMAG=(0.0D0,1.0D0), cZERO=(0.0D0,0.0D0), cONE=(1.0D0,0.0D0))
      DOUBLE COMPLEX a1(6,-1:1), a2(6,-1:1), h(6,-1:1), z(6,-1:1)      
      DOUBLE COMPLEX M(-1:1,-1:1,-1:1,-1:1)
      
      DOUBLE COMPLEX OFFeps(6,-1:1,NGRAPHS), OFFepsC(6,-1:1,NGRAPHS)
      DOUBLE PRECISION K(0:3,4)      
      DOUBLE COMPLEX O(-1:1,4,-1:1,4), OP(4,-1:1,4)
      DOUBLE COMPLEXPP(1:4,1:4)
    
      INTEGER I,J,N,MU,NU
      INTEGER LAMBDA1, LAMBDA2, LAMBDA3, LAMBDA4, LAMBDA
      INTEGER L1,L2

      DOUBLE COMPLEX DOTCC, DOTRC,contract_Tjj
      DOUBLE PRECISION DOTRR
      EXTERNAL DOTCC, DOTRC, DOTRR, contract_Tjj
      double precision sw2, sw
#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "gmunu.inc"
#include "transm.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"
      sw2 = sin2w
      sw = sqrt(sin2w)
      OFFeps(6,-1:1,4) = cZero
      K(0:3,1:4) = rZero
      M(-1:1,-1:1,-1:1,-1:1) = cZero
      PP(1:4,1:4) = rZero
      O(-1:1,1:4,-1:1,1:4) = cZero
      OP(1:4,-1:1,1:4) = cZero  
      
      OFFeps(1:6,-1:1,1) = a1(1:6,-1:1)
      OFFeps(1:6,-1:1,2) = a2(1:6,-1:1)
      OFFeps(1:6,-1:1,3) = h(1:6,-1:1)
      OFFeps(1:6,-1:1,4) = z(1:6,-1:1) 

      do i = 1,4
        K(0,i) = dble(OFFeps(5,-1,i))*(-rONE) !K1(0:3)
        K(1,i) = dble(OFFeps(6,-1,i))*(-rONE) !K2(0:3)
        K(2,i) = dimag(OFFeps(6,-1,i))*(-rONE)!K3(0:3)
        K(3,i) = dimag(OFFeps(5,-1,i))*(-rONE)!K4(0:3)
      enddo
      



      do i = 1,4
        do j = 1,4
          pp(i,j) = dotrr(k(0,i),k(0,j))
        enddo
      enddo

      do L1 = -1,1
        do L2 = -1,1 
          do i = 1,4
            do j = 1,4
              O(L1,i,L2,j) = dotcc(OFFeps(1,L1,i),OFFeps(1,L2,j))
            enddo
          enddo
        enddo
      enddo

      do LAMBDA = -1,1
        do i = 1,4 
          do j = 1,4
            OP(j,lambda,i) = dotrc(k(0,j),OFFeps(1,lambda,i))
          enddo
        enddo
      enddo


!!!! ===============================================================      
!!!! COUPLINGS CALCULATION
!!!! ===============================================================



        

      end


ccccccccccccccc
c aazz-vertex c
ccccccccccccccc

      subroutine aazz_ONanom4(a1,a2,z1,z2,M)
c
c This subroutine calculates the vvvv-vertex with anomalous couplings
c and the on-shell polarization vectors for vector bosons
c
c IN:
c     polv1, polv2      : the incoming polvector wz(bosons in the CM frame)
c     polv3, polv4      : the outgoing polvector wz(bosons in the CM frame)
c OUT:
c     M(lambda1,lambda2) 
c
      IMPLICIT NONE

      INTEGER NGRAPHS
      PARAMETER( NGRAPHS=20 )
      DOUBLE PRECISION rZERO, rONE, rHalf
      PARAMETER( rZERO=0.0D0, rONE=1.0D0, rHALF = 0.5D0 )
      DOUBLE COMPLEX cIMAG, cZERO, cONE
      PARAMETER( cIMAG=(0.0D0,1.0D0), cZERO=(0.0D0,0.0D0), cONE=(1.0D0,0.0D0))
      DOUBLE COMPLEX a1(6,-1:1), a2(6,-1:1), z1(6,-1:1), z2(6,-1:1)      
      DOUBLE COMPLEX M(-1:1,-1:1,-1:1,-1:1)
      
      DOUBLE COMPLEX OFFeps(6,-1:1,NGRAPHS), OFFepsC(6,-1:1,NGRAPHS)
      DOUBLE PRECISION K(0:3,4)      
      DOUBLE COMPLEX O(-1:1,4,-1:1,4), OP(4,-1:1,4)
      DOUBLE COMPLEXPP(1:4,1:4)
    
      INTEGER I,J,N,MU,NU
      INTEGER LAMBDA1, LAMBDA2, LAMBDA3, LAMBDA4, LAMBDA
      INTEGER L1,L2

      DOUBLE COMPLEX DOTCC, DOTRC,contract_Tjj
      DOUBLE PRECISION DOTRR
      EXTERNAL DOTCC, DOTRC, DOTRR, contract_Tjj
      double precision sw2, sw
#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "gmunu.inc"
#include "transm.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"
      sw2 = sin2w
      sw = sqrt(sin2w)
      OFFeps(6,-1:1,4) = cZero
      K(0:3,1:4) = rZero
      M(-1:1,-1:1,-1:1,-1:1) = cZero
      PP(1:4,1:4) = rZero
      O(-1:1,1:4,-1:1,1:4) = cZero
      OP(1:4,-1:1,1:4) = cZero  
      
      OFFeps(1:6,-1:1,1) = a1(1:6,-1:1)
      OFFeps(1:6,-1:1,2) = a2(1:6,-1:1)
      OFFeps(1:6,-1:1,3) = z1(1:6,-1:1)
      OFFeps(1:6,-1:1,4) = z2(1:6,-1:1) 

      do i = 1,4
        K(0,i) = dble(OFFeps(5,-1,i))*(-rONE) !K1(0:3)
        K(1,i) = dble(OFFeps(6,-1,i))*(-rONE) !K2(0:3)
        K(2,i) = dimag(OFFeps(6,-1,i))*(-rONE)!K3(0:3)
        K(3,i) = dimag(OFFeps(5,-1,i))*(-rONE)!K4(0:3)
      enddo
      



      do i = 1,4
        do j = 1,4
          pp(i,j) = dotrr(k(0,i),k(0,j))
        enddo
      enddo

      do L1 = -1,1
        do L2 = -1,1 
          do i = 1,4
            do j = 1,4
              O(L1,i,L2,j) = dotcc(OFFeps(1,L1,i),OFFeps(1,L2,j))
            enddo
          enddo
        enddo
      enddo

      do LAMBDA = -1,1
        do i = 1,4 
          do j = 1,4
            OP(j,lambda,i) = dotrc(k(0,j),OFFeps(1,lambda,i))
          enddo
        enddo
      enddo


!!!! ===============================================================      
!!!! COUPLINGS CALCULATION
!!!! ===============================================================



                 

         if( ft0.ne.0d0 ) then  !ft0 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) -8*ft0*(gw**4)*SW2*(-1+SW2)*(OP(1,lambda2,2)*OP(2,lambda1,1)*OP(3,lambda4,4)*OP(4,lambda3,3)-OP(3,lambda4,4)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda2,2)-OP(2,lambda4,4)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda3,3)-OP(2,lambda3,3)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*OP(4,lambda1,1)*(OP(2,lambda3,3)*OP(3,lambda2,2)-pp(2,3)*O(lambda2,2,lambda3,3))+pp(1,3)*pp(2,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*OP(3,lambda1,1)*(OP(2,lambda4,4)*OP(4,lambda2,2)-pp(2,4)*O(lambda2,2,lambda4,4))-OP(1,lambda2,2)*OP(2,lambda1,1)*pp(3,4)*O(lambda3,3,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4))

               enddo                      
              enddo
             enddo 
            enddo
         endif


         if( ft1.ne.0d0 ) then  !ft1 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) -8*ft1*(gw**4)*SW2*(-1+SW2)*(OP(1,lambda2,2)*OP(2,lambda1,1)*OP(3,lambda4,4)*OP(4,lambda3,3)-OP(3,lambda4,4)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda2,2)-OP(2,lambda4,4)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda3,3)-OP(2,lambda3,3)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*OP(4,lambda1,1)*(OP(2,lambda3,3)*OP(3,lambda2,2)-pp(2,3)*O(lambda2,2,lambda3,3))+pp(1,3)*pp(2,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*OP(3,lambda1,1)*(OP(2,lambda4,4)*OP(4,lambda2,2)-pp(2,4)*O(lambda2,2,lambda4,4))-OP(1,lambda2,2)*OP(2,lambda1,1)*pp(3,4)*O(lambda3,3,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4))

               enddo                      
              enddo
             enddo 
            enddo
         endif


         if( ft2.ne.0d0 ) then  !ft2 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) -2*ft2*(gw**4)*SW2*(-1+SW2)*(OP(1,lambda4,4)*OP(2,lambda3,3)*OP(3,lambda1,1)*OP(4,lambda2,2)+OP(1,lambda4,4)*OP(2,lambda1,1)*OP(3,lambda2,2)*OP(4,lambda3,3)-OP(2,lambda4,4)*OP(4,lambda3,3)*pp(1,3)*O(lambda1,1,lambda2,2)-OP(2,lambda3,3)*OP(3,lambda4,4)*pp(1,4)*O(lambda1,1,lambda2,2)-OP(1,lambda4,4)*OP(4,lambda3,3)*pp(2,3)*O(lambda1,1,lambda2,2)+OP(1,lambda4,4)*OP(2,lambda3,3)*pp(3,4)*O(lambda1,1,lambda2,2)-OP(3,lambda4,4)*OP(4,lambda2,2)*pp(1,2)*O(lambda1,1,lambda3,3)-OP(2,lambda4,4)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda3,3)-OP(1,lambda4,4)*OP(4,lambda2,2)*pp(2,3)*O(lambda1,1,lambda3,3)+OP(1,lambda4,4)*OP(3,lambda2,2)*pp(2,4)*O(lambda1,1,lambda3,3)-OP(3,lambda2,2)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda4,4)-OP(2,lambda3,3)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda4,4)-OP(3,lambda4,4)*OP(4,lambda1,1)*pp(1,2)*O(lambda2,2,lambda3,3)-OP(2,lambda4,4)*OP(4,lambda1,1)*pp(1,3)*O(lambda2,2,lambda3,3)+OP(2,lambda4,4)*OP(3,lambda1,1)*pp(1,4)*O(lambda2,2,lambda3,3)+OP(2,lambda1,1)*OP(3,lambda4,4)*pp(1,4)*O(lambda2,2,lambda3,3)-OP(1,lambda4,4)*OP(3,lambda1,1)*pp(2,4)*O(lambda2,2,lambda3,3)-OP(1,lambda4,4)*OP(2,lambda1,1)*pp(3,4)*O(lambda2,2,lambda3,3)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-OP(3,lambda1,1)*OP(4,lambda3,3)*pp(1,2)*O(lambda2,2,lambda4,4)+OP(2,lambda3,3)*OP(4,lambda1,1)*pp(1,3)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*OP(4,lambda3,3)*pp(1,3)*O(lambda2,2,lambda4,4)-OP(2,lambda3,3)*OP(3,lambda1,1)*pp(1,4)*O(lambda2,2,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(-(OP(3,lambda4,4)*pp(2,4)*O(lambda1,1,lambda2,2))+OP(2,lambda4,4)*(OP(3,lambda2,2)*OP(4,lambda1,1)+pp(3,4)*O(lambda1,1,lambda2,2))+OP(4,lambda2,2)*pp(2,3)*O(lambda1,1,lambda4,4)-OP(3,lambda2,2)*pp(2,4)*O(lambda1,1,lambda4,4)-OP(4,lambda1,1)*pp(2,3)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*(OP(3,lambda4,4)*OP(4,lambda2,2)-pp(3,4)*O(lambda2,2,lambda4,4)))+OP(3,lambda2,2)*OP(4,lambda1,1)*pp(1,2)*O(lambda3,3,lambda4,4)+OP(3,lambda1,1)*OP(4,lambda2,2)*pp(1,2)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(4,lambda2,2)*pp(1,3)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(3,lambda2,2)*pp(1,4)*O(lambda3,3,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+OP(1,lambda2,2)*(OP(3,lambda4,4)*pp(2,4)*O(lambda1,1,lambda3,3)+OP(2,lambda4,4)*(OP(3,lambda1,1)*OP(4,lambda3,3)-pp(3,4)*O(lambda1,1,lambda3,3))+OP(4,lambda3,3)*pp(2,3)*O(lambda1,1,lambda4,4)+OP(2,lambda3,3)*(OP(3,lambda4,4)*OP(4,lambda1,1)-pp(3,4)*O(lambda1,1,lambda4,4))-OP(4,lambda1,1)*pp(2,3)*O(lambda3,3,lambda4,4)-OP(3,lambda1,1)*pp(2,4)*O(lambda3,3,lambda4,4)))

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( ft3.ne.0d0 ) then  !ft3 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) -2*ft3*(SW2-1)*gw**4*sw**2*(-pp(1,2)*OP(4,lambda2,2)*OP(3,lambda4,4)*O(lambda1,1,lambda3,3)-pp(2,3)*OP(4,lambda2,2)*OP(1,lambda4,4)*O(lambda1,1,lambda3,3)-pp(1,3)*OP(4,lambda2,2)*OP(2,lambda3,3)*O(lambda1,1,lambda4,4)+pp(1,2)*OP(3,lambda1,1)*OP(4,lambda2,2)*O(lambda3,3,lambda4,4)-pp(1,3)*OP(2,lambda1,1)*OP(4,lambda2,2)*O(lambda3,3,lambda4,4)-pp(1,3)*OP(4,lambda3,3)*OP(2,lambda4,4)*O(lambda1,1,lambda2,2)-pp(1,4)*OP(2,lambda3,3)*OP(3,lambda4,4)*O(lambda1,1,lambda2,2)-pp(2,3)*OP(4,lambda3,3)*OP(1,lambda4,4)*O(lambda1,1,lambda2,2)+pp(3,4)*OP(2,lambda3,3)*OP(1,lambda4,4)*O(lambda1,1,lambda2,2)-pp(1,4)*OP(3,lambda2,2)*OP(2,lambda4,4)*O(lambda1,1,lambda3,3)+pp(2,4)*OP(3,lambda2,2)*OP(1,lambda4,4)*O(lambda1,1,lambda3,3)-pp(1,2)*OP(3,lambda2,2)*OP(4,lambda3,3)*O(lambda1,1,lambda4,4)-pp(1,2)*OP(4,lambda1,1)*OP(3,lambda4,4)*O(lambda2,2,lambda3,3)-pp(1,3)*OP(4,lambda1,1)*OP(2,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,4)*OP(3,lambda1,1)*OP(2,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,4)*OP(2,lambda1,1)*OP(3,lambda4,4)*O(lambda2,2,lambda3,3)-pp(2,4)*OP(3,lambda1,1)*OP(1,lambda4,4)*O(lambda2,2,lambda3,3)-pp(3,4)*OP(2,lambda1,1)*OP(1,lambda4,4)*O(lambda2,2,lambda3,3)-pp(1,2)*OP(3,lambda1,1)*OP(4,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,3)*OP(4,lambda1,1)*OP(2,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,3)*OP(2,lambda1,1)*OP(4,lambda3,3)*O(lambda2,2,lambda4,4)-pp(1,4)*OP(3,lambda1,1)*OP(2,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(pp(2,4)*(-OP(3,lambda4,4))*O(lambda1,1,lambda2,2)+OP(2,lambda4,4)*(OP(4,lambda1,1)*OP(3,lambda2,2)+pp(3,4)*O(lambda1,1,lambda2,2))+pp(2,3)*OP(4,lambda2,2)*O(lambda1,1,lambda4,4)-pp(2,4)*OP(3,lambda2,2)*O(lambda1,1,lambda4,4)-pp(2,3)*OP(4,lambda1,1)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*(OP(4,lambda2,2)*OP(3,lambda4,4)-pp(3,4)*O(lambda2,2,lambda4,4)))+pp(1,2)*OP(4,lambda1,1)*OP(3,lambda2,2)*O(lambda3,3,lambda4,4)-pp(1,4)*OP(2,lambda1,1)*OP(3,lambda2,2)*O(lambda3,3,lambda4,4)+OP(1,lambda2,2)*(pp(2,4)*OP(3,lambda4,4)*O(lambda1,1,lambda3,3)+OP(2,lambda4,4)*(OP(3,lambda1,1)*OP(4,lambda3,3)-pp(3,4)*O(lambda1,1,lambda3,3))+pp(2,3)*OP(4,lambda3,3)*O(lambda1,1,lambda4,4)+OP(2,lambda3,3)*(OP(4,lambda1,1)*OP(3,lambda4,4)-pp(3,4)*O(lambda1,1,lambda4,4))-pp(2,3)*OP(4,lambda1,1)*O(lambda3,3,lambda4,4)-pp(2,4)*OP(3,lambda1,1)*O(lambda3,3,lambda4,4))+OP(3,lambda1,1)*OP(4,lambda2,2)*OP(2,lambda3,3)*OP(1,lambda4,4)+OP(2,lambda1,1)*OP(3,lambda2,2)*OP(4,lambda3,3)*OP(1,lambda4,4)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4))

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( ft4.ne.0d0 ) then  !ft4 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(gw**4*ft4*sw**2*ZMASS**2*(sw**4*OP(1, lambda4, 4)*OP(2, lambda3, 3)*OP(3, lambda1, 1)*OP(4, lambda2, 2) + (WMASS**4*OP(1, lambda4, 4)*OP(2, lambda3, 3)*OP(3, lambda1, 1)*OP(4, lambda2, 2))/ZMASS**4 - (2*SW2*WMASS**2*OP(1, lambda4, 4)*OP(2, lambda1, 1)*OP(3, lambda2, 2)*OP(4, lambda3, 3))/ZMASS**2 - (2*SW2*WMASS**2*pp(3, 4)*OP(1, lambda4, 4)*OP(2, lambda3, 3)*O(lambda1, 1, lambda2, 2))/ZMASS**2 + (2*SW2*WMASS**2*pp(1, 4)*OP(2, lambda3, 3)*OP(3, lambda4, 4)*O(lambda1, 1, lambda2, 2))/ZMASS**2 + (2*SW2*WMASS**2*pp(2, 3)*OP(1, lambda4, 4)*OP(4, lambda3, 3)*O(lambda1, 1, lambda2, 2))/ZMASS**2 + (2*SW2*WMASS**2*pp(1, 3)*OP(2, lambda4, 4)*OP(4, lambda3, 3)*O(lambda1, 1, lambda2, 2))/ZMASS**2 + sw**4*pp(2, 4)*OP(1, lambda4, 4)*OP(3, lambda2, 2)*O(lambda1, 1, lambda3, 3) + (WMASS**4*pp(2, 4)*OP(1, lambda4, 4)*OP(3, lambda2, 2)*O(lambda1, 1, lambda3, 3))/ZMASS**4 - sw**4*pp(1, 4)*OP(2, lambda4, 4)*OP(3, lambda2, 2)*O(lambda1, 1, lambda3, 3) - (WMASS**4*pp(1, 4)*OP(2, lambda4, 4)*OP(3, lambda2, 2)*O(lambda1, 1, lambda3, 3))/ZMASS**4 - sw**4*pp(2, 3)*OP(1, lambda4, 4)*OP(4, lambda2, 2)*O(lambda1, 1, lambda3, 3) - (WMASS**4*pp(2, 3)*OP(1, lambda4, 4)*OP(4, lambda2, 2)*O(lambda1, 1, lambda3, 3))/ZMASS**4 + (2*SW2*WMASS**2*pp(1, 2)*OP(3, lambda4, 4)*OP(4, lambda2, 2)*O(lambda1, 1, lambda3, 3))/ZMASS**2 - sw**4*pp(1, 3)*OP(2, lambda3, 3)*OP(4, lambda2, 2)*O(lambda1, 1, lambda4, 4) - (WMASS**4*pp(1, 3)*OP(2, lambda3, 3)*OP(4, lambda2, 2)*O(lambda1, 1, lambda4, 4))/ZMASS**4 + (2*SW2*WMASS**2*pp(1, 2)*OP(3, lambda2, 2)*OP(4, lambda3, 3)*O(lambda1, 1, lambda4, 4))/ZMASS**2 + (2*SW2*WMASS**2*pp(3, 4)*OP(1, lambda4, 4)*OP(2, lambda1, 1)*O(lambda2, 2, lambda3, 3))/ZMASS**2 - sw**4*pp(2, 4)*OP(1, lambda4, 4)*OP(3, lambda1, 1)*O(lambda2, 2, lambda3, 3) - (WMASS**4*pp(2, 4)*OP(1, lambda4, 4)*OP(3, lambda1, 1)*O(lambda2, 2, lambda3, 3))/ZMASS**4 + sw**4*pp(1, 4)*OP(2, lambda4, 4)*OP(3, lambda1, 1)*O(lambda2, 2, lambda3, 3) + (WMASS**4*pp(1, 4)*OP(2, lambda4, 4)*OP(3, lambda1, 1)*O(lambda2, 2, lambda3, 3))/ZMASS**4 - (2*SW2*WMASS**2*pp(1, 4)*OP(2, lambda1, 1)*OP(3, lambda4, 4)*O(lambda2, 2, lambda3, 3))/ZMASS**2 - sw**4*pp(1, 3)*OP(2, lambda4, 4)*OP(4, lambda1, 1)*O(lambda2, 2, lambda3, 3) - (WMASS**4*pp(1, 3)*OP(2, lambda4, 4)*OP(4, lambda1, 1)*O(lambda2, 2, lambda3, 3))/ZMASS**4 + (2*SW2*WMASS**2*pp(1, 2)*OP(3, lambda4, 4)*OP(4, lambda1, 1)*O(lambda2, 2, lambda3, 3))/ZMASS**2 + sw**4*pp(1, 3)*pp(2, 4)*O(lambda1, 1, lambda4, 4)*O(lambda2, 2, lambda3, 3) + (WMASS**4*pp(1, 3)*pp(2, 4)*O(lambda1, 1, lambda4, 4)*O(lambda2, 2, lambda3, 3))/ZMASS**4 - (2*SW2*WMASS**2*pp(1, 2)*pp(3, 4)*O(lambda1, 1, lambda4, 4)*O(lambda2, 2, lambda3, 3))/ZMASS**2 - sw**4*pp(1, 4)*OP(2, lambda3, 3)*OP(3, lambda1, 1)*O(lambda2, 2, lambda4, 4) - (WMASS**4*pp(1, 4)*OP(2, lambda3, 3)*OP(3, lambda1, 1)*O(lambda2, 2, lambda4, 4))/ZMASS**4 + sw**4*pp(1, 3)*OP(2, lambda3, 3)*OP(4, lambda1, 1)*O(lambda2, 2, lambda4, 4) + (WMASS**4*pp(1, 3)*OP(2, lambda3, 3)*OP(4, lambda1, 1)*O(lambda2, 2, lambda4, 4))/ZMASS**4 - (2*SW2*WMASS**2*pp(1, 3)*OP(2, lambda1, 1)*OP(4, lambda3, 3)*O(lambda2, 2, lambda4, 4))/ZMASS**2 + (2*SW2*WMASS**2*pp(1, 2)*OP(3, lambda1, 1)*OP(4, lambda3, 3)*O(lambda2, 2, lambda4, 4))/ZMASS**2 + sw**4*pp(1, 4)*pp(2, 3)*O(lambda1, 1, lambda3, 3)*O(lambda2, 2, lambda4, 4) + (WMASS**4*pp(1, 4)*pp(2, 3)*O(lambda1, 1, lambda3, 3)*O(lambda2, 2, lambda4, 4))/ZMASS**4 - (2*SW2*WMASS**2*pp(1, 2)*pp(3, 4)*O(lambda1, 1, lambda3, 3)*O(lambda2, 2, lambda4, 4))/ZMASS**2 + OP(1, lambda3, 3)*((2*SW2*WMASS**2*pp(2, 4)*OP(3, lambda4, 4)*O(lambda1, 1, lambda2, 2))/ZMASS**2 + OP(2, lambda4, 4)*((1 + 2*sw**4 - 2*SW2)*OP(3, lambda2, 2)*OP(4, lambda1, 1) -     (2*SW2*WMASS**2*pp(3, 4)*O(lambda1, 1, lambda2, 2))/ZMASS**2) - sw**4*pp(2, 4)*OP(3, lambda2, 2)*O(lambda1, 1, lambda4, 4) - (WMASS**4*pp(2, 4)*OP(3, lambda2, 2)*O(lambda1, 1, lambda4, 4))/ZMASS**4 +   sw**4*pp(2, 3)*OP(4, lambda2, 2)*O(lambda1, 1, lambda4, 4) + (WMASS**4*pp(2, 3)*OP(4, lambda2, 2)*O(lambda1, 1, lambda4, 4))/ZMASS**4 - sw**4*pp(2, 3)*OP(4, lambda1, 1)*O(lambda2, 2, lambda4, 4) -   (WMASS**4*pp(2, 3)*OP(4, lambda1, 1)*O(lambda2, 2, lambda4, 4))/ZMASS**4 + 2*(-1 + SW2)*SW2*OP(2, lambda1, 1)*(OP(3, lambda4, 4)*OP(4, lambda2, 2) - pp(3, 4)*O(lambda2, 2, lambda4, 4))) + (2*SW2*WMASS**2*pp(1, 4)*OP(2, lambda1, 1)*OP(3, lambda2, 2)*O(lambda3, 3, lambda4, 4))/ZMASS**2 - (2*SW2*WMASS**2*pp(1, 2)*OP(3, lambda2, 2)*OP(4, lambda1, 1)*O(lambda3, 3, lambda4, 4))/ZMASS**2 + (2*SW2*WMASS**2*pp(1, 3)*OP(2, lambda1, 1)*OP(4, lambda2, 2)*O(lambda3, 3, lambda4, 4))/ZMASS**2 - (2*SW2*WMASS**2*pp(1, 2)*OP(3, lambda1, 1)*OP(4, lambda2, 2)*O(lambda3, 3, lambda4, 4))/ZMASS**2 - (2*SW2*WMASS**2*pp(1, 4)*pp(2, 3)*O(lambda1, 1, lambda2, 2)*O(lambda3, 3, lambda4, 4))/ZMASS**2 - (2*SW2*WMASS**2*pp(1, 3)*pp(2, 4)*O(lambda1, 1, lambda2, 2)*O(lambda3, 3, lambda4, 4))/ZMASS**2 + 2*(-1 + SW2)*SW2*OP(1, lambda2, 2)*(pp(2, 4)*OP(3, lambda4, 4)*O(lambda1, 1, lambda3, 3) + OP(2, lambda4, 4)*(OP(3, lambda1, 1)*OP(4, lambda3, 3) - pp(3, 4)*O(lambda1, 1, lambda3, 3)) +   pp(2, 3)*OP(4, lambda3, 3)*O(lambda1, 1, lambda4, 4) + OP(2, lambda3, 3)*(OP(3, lambda4, 4)*OP(4, lambda1, 1) - pp(3, 4)*O(lambda1, 1, lambda4, 4)) -   pp(2, 4)*OP(3, lambda1, 1)*O(lambda3, 3, lambda4, 4) - pp(2, 3)*OP(4, lambda1, 1)*O(lambda3, 3, lambda4, 4))))/(2*WMASS**2)

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( ft5.ne.0d0 ) then  !ft5 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(2*ft5*(gw**4)*SW2*(ZMASS**2)*((sw**4)*OP(1,lambda2,2)*OP(2,lambda1,1)*OP(3,lambda4,4)*OP(4,lambda3,3)+((WMASS**4)*OP(1,lambda2,2)*OP(2,lambda1,1)*OP(3,lambda4,4)*OP(4,lambda3,3))/(ZMASS**4)-(sw**4)*OP(3,lambda4,4)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda2,2)-((WMASS**4)*OP(3,lambda4,4)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda2,2))/(ZMASS**4)-2*(sw**4)*OP(2,lambda4,4)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda3,3)+2*SW2*OP(2,lambda4,4)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda3,3)-2*(sw**4)*OP(2,lambda3,3)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda4,4)+2*SW2*OP(2,lambda3,3)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda4,4)+2*(sw**4)*pp(1,4)*pp(2,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-2*SW2*pp(1,4)*pp(2,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+2*(-1+SW2)*SW2*OP(1,lambda4,4)*OP(4,lambda1,1)*(OP(2,lambda3,3)*OP(3,lambda2,2)-pp(2,3)*O(lambda2,2,lambda3,3))+2*(sw**4)*pp(1,3)*pp(2,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)-2*SW2*pp(1,3)*pp(2,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+2*(-1+SW2)*SW2*OP(1,lambda3,3)*OP(3,lambda1,1)*(OP(2,lambda4,4)*OP(4,lambda2,2)-pp(2,4)*O(lambda2,2,lambda4,4))-(sw**4)*OP(1,lambda2,2)*OP(2,lambda1,1)*pp(3,4)*O(lambda3,3,lambda4,4)-((WMASS**4)*OP(1,lambda2,2)*OP(2,lambda1,1)*pp(3,4)*O(lambda3,3,lambda4,4))/(ZMASS**4)+(sw**4)*pp(1,2)*pp(3,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+((WMASS**4)*pp(1,2)*pp(3,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4))/(ZMASS**4)))/(WMASS**2)

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( ft6.ne.0d0 ) then  !ft6 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(ft6*(gw**4)*SW2*(ZMASS**2)*
     &             (-(OP(2,lambda4,4)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda3,3))
     &             -OP(2,lambda3,3)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda4,4)
     &             +pp(1,4)*pp(2,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)
     &             +pp(1,3)*pp(2,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)
     &             +(1 - 2*SW2)**2*(OP(1,lambda4,4)*OP(4,lambda1,1)
     &             *(OP(2,lambda3,3)*OP(3,lambda2,2)-pp(2,3)*O(lambda2,2,lambda3,3))
     &             +OP(1,lambda3,3)*OP(3,lambda1,1)
     &             *(OP(2,lambda4,4)*OP(4,lambda2,2)-pp(2,4)*O(lambda2,2,lambda4,4)))
     &             +4*(sw**4)
     &             *(OP(1,lambda2,2)*OP(2,lambda1,1)*OP(3,lambda4,4)*OP(4,lambda3,3)
     &             -OP(3,lambda4,4)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda2,2)
     &             -OP(2,lambda4,4)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda3,3)
     &             -OP(2,lambda3,3)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda4,4)
     &             +pp(1,4)*pp(2,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)
     &             +pp(1,3)*pp(2,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)
     &             -OP(1,lambda2,2)*OP(2,lambda1,1)*pp(3,4)*O(lambda3,3,lambda4,4)
     &             +pp(1,2)*pp(3,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4))
     &             -4*SW2
     &             *(-(OP(3,lambda4,4)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda2,2))
     &             -OP(2,lambda4,4)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda3,3)
     &             -OP(2,lambda3,3)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda4,4)
     &             +pp(1,4)*pp(2,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)
     &             +pp(1,3)*pp(2,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)
     &             +pp(1,2)*pp(3,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)
     &             +OP(1,lambda2,2)*OP(2,lambda1,1)
     &             *(OP(3,lambda4,4)*OP(4,lambda3,3)-pp(3,4)*O(lambda3,3,lambda4,4)))))/(WMASS**2)
               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( ft7.ne.0d0 ) then  !ft7 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(ft7*(gw**4)*SW2*(ZMASS**2)*(OP(1,lambda4,4)*OP(2,lambda1,1)*OP(3,lambda2,2)*OP(4,lambda3,3)-OP(2,lambda4,4)*OP(4,lambda3,3)*pp(1,3)*O(lambda1,1,lambda2,2)-OP(2,lambda3,3)*OP(3,lambda4,4)*pp(1,4)*O(lambda1,1,lambda2,2)-OP(1,lambda4,4)*OP(4,lambda3,3)*pp(2,3)*O(lambda1,1,lambda2,2)+OP(1,lambda4,4)*OP(2,lambda3,3)*pp(3,4)*O(lambda1,1,lambda2,2)-OP(3,lambda4,4)*OP(4,lambda2,2)*pp(1,2)*O(lambda1,1,lambda3,3)-OP(3,lambda2,2)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda4,4)-OP(3,lambda4,4)*OP(4,lambda1,1)*pp(1,2)*O(lambda2,2,lambda3,3)+OP(2,lambda1,1)*OP(3,lambda4,4)*pp(1,4)*O(lambda2,2,lambda3,3)-OP(1,lambda4,4)*OP(2,lambda1,1)*pp(3,4)*O(lambda2,2,lambda3,3)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-OP(3,lambda1,1)*OP(4,lambda3,3)*pp(1,2)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*OP(4,lambda3,3)*pp(1,3)*O(lambda2,2,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(-(OP(3,lambda4,4)*pp(2,4)*O(lambda1,1,lambda2,2))+OP(2,lambda4,4)*(4*(-1+SW2)*SW2*OP(3,lambda2,2)*OP(4,lambda1,1)+(1-2*SW2)**2*pp(3,4)*O(lambda1,1,lambda2,2))+4*SW2*(-(OP(4,lambda2,2)*pp(2,3)*O(lambda1,1,lambda4,4))+pp(2,4)*(OP(3,lambda4,4)*O(lambda1,1,lambda2,2)+OP(3,lambda2,2)*O(lambda1,1,lambda4,4))+OP(4,lambda1,1)*pp(2,3)*O(lambda2,2,lambda4,4))+(1-2*SW2)**2*OP(2,lambda1,1)*(OP(3,lambda4,4)*OP(4,lambda2,2)-pp(3,4)*O(lambda2,2,lambda4,4))+(sw**4)*(4*OP(4,lambda2,2)*pp(2,3)*O(lambda1,1,lambda4,4)-4*(pp(2,4)*(OP(3,lambda4,4)*O(lambda1,1,lambda2,2)+OP(3,lambda2,2)*O(lambda1,1,lambda4,4))+OP(4,lambda1,1)*pp(2,3)*O(lambda2,2,lambda4,4))))+OP(3,lambda2,2)*OP(4,lambda1,1)*pp(1,2)*O(lambda3,3,lambda4,4)+OP(3,lambda1,1)*OP(4,lambda2,2)*pp(1,2)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(4,lambda2,2)*pp(1,3)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(3,lambda2,2)*pp(1,4)*O(lambda3,3,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+(1-2*SW2)**2*OP(1,lambda2,2)*(OP(3,lambda4,4)*pp(2,4)*O(lambda1,1,lambda3,3)+OP(2,lambda4,4)*(OP(3,lambda1,1)*OP(4,lambda3,3)-pp(3,4)*O(lambda1,1,lambda3,3))+OP(4,lambda3,3)*pp(2,3)*O(lambda1,1,lambda4,4)+OP(2,lambda3,3)*(OP(3,lambda4,4)*OP(4,lambda1,1)-pp(3,4)*O(lambda1,1,lambda4,4))-OP(4,lambda1,1)*pp(2,3)*O(lambda3,3,lambda4,4)-OP(3,lambda1,1)*pp(2,4)*O(lambda3,3,lambda4,4))+4*(sw**4)*(-(OP(2,lambda4,4)*OP(4,lambda3,3)*pp(1,3)*O(lambda1,1,lambda2,2))-OP(2,lambda4,4)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda3,3)-OP(3,lambda2,2)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda4,4)-OP(2,lambda3,3)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda4,4)-OP(2,lambda4,4)*OP(4,lambda1,1)*pp(1,3)*O(lambda2,2,lambda3,3)+OP(2,lambda4,4)*OP(3,lambda1,1)*pp(1,4)*O(lambda2,2,lambda3,3)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-OP(3,lambda4,4)*(OP(4,lambda2,2)*pp(1,2)*O(lambda1,1,lambda3,3)+OP(4,lambda1,1)*pp(1,2)*O(lambda2,2,lambda3,3)+pp(1,4)*(OP(2,lambda3,3)*O(lambda1,1,lambda2,2)-OP(2,lambda1,1)*O(lambda2,2,lambda3,3)))+OP(1,lambda4,4)*(-(OP(4,lambda3,3)*pp(2,3)*O(lambda1,1,lambda2,2))+OP(2,lambda3,3)*(OP(3,lambda1,1)*OP(4,lambda2,2)+pp(3,4)*O(lambda1,1,lambda2,2))-OP(4,lambda2,2)*pp(2,3)*O(lambda1,1,lambda3,3)+OP(3,lambda2,2)*pp(2,4)*O(lambda1,1,lambda3,3)-OP(3,lambda1,1)*pp(2,4)*O(lambda2,2,lambda3,3)+OP(2,lambda1,1)*(OP(3,lambda2,2)*OP(4,lambda3,3)-pp(3,4)*O(lambda2,2,lambda3,3)))-OP(3,lambda1,1)*OP(4,lambda3,3)*pp(1,2)*O(lambda2,2,lambda4,4)+OP(2,lambda3,3)*OP(4,lambda1,1)*pp(1,3)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*OP(4,lambda3,3)*pp(1,3)*O(lambda2,2,lambda4,4)-OP(2,lambda3,3)*OP(3,lambda1,1)*pp(1,4)*O(lambda2,2,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(3,lambda2,2)*OP(4,lambda1,1)*pp(1,2)*O(lambda3,3,lambda4,4)+OP(3,lambda1,1)*OP(4,lambda2,2)*pp(1,2)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(4,lambda2,2)*pp(1,3)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(3,lambda2,2)*pp(1,4)*O(lambda3,3,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4))-4*SW2*(-(OP(2,lambda4,4)*OP(4,lambda3,3)*pp(1,3)*O(lambda1,1,lambda2,2))-OP(2,lambda4,4)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda3,3)-OP(3,lambda2,2)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda4,4)-OP(2,lambda3,3)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda4,4)-OP(2,lambda4,4)*OP(4,lambda1,1)*pp(1,3)*O(lambda2,2,lambda3,3)+OP(2,lambda4,4)*OP(3,lambda1,1)*pp(1,4)*O(lambda2,2,lambda3,3)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-OP(3,lambda4,4)*(OP(4,lambda2,2)*pp(1,2)*O(lambda1,1,lambda3,3)+OP(4,lambda1,1)*pp(1,2)*O(lambda2,2,lambda3,3)+pp(1,4)*(OP(2,lambda3,3)*O(lambda1,1,lambda2,2)-OP(2,lambda1,1)*O(lambda2,2,lambda3,3)))+OP(1,lambda4,4)*(-(OP(4,lambda3,3)*pp(2,3)*O(lambda1,1,lambda2,2))+OP(2,lambda3,3)*(OP(3,lambda1,1)*OP(4,lambda2,2)+pp(3,4)*O(lambda1,1,lambda2,2))-OP(4,lambda2,2)*pp(2,3)*O(lambda1,1,lambda3,3)+OP(3,lambda2,2)*pp(2,4)*O(lambda1,1,lambda3,3)-OP(3,lambda1,1)*pp(2,4)*O(lambda2,2,lambda3,3)+OP(2,lambda1,1)*(OP(3,lambda2,2)*OP(4,lambda3,3)-pp(3,4)*O(lambda2,2,lambda3,3)))-OP(3,lambda1,1)*OP(4,lambda3,3)*pp(1,2)*O(lambda2,2,lambda4,4)+OP(2,lambda3,3)*OP(4,lambda1,1)*pp(1,3)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*OP(4,lambda3,3)*pp(1,3)*O(lambda2,2,lambda4,4)-OP(2,lambda3,3)*OP(3,lambda1,1)*pp(1,4)*O(lambda2,2,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(3,lambda2,2)*OP(4,lambda1,1)*pp(1,2)*O(lambda3,3,lambda4,4)+OP(3,lambda1,1)*OP(4,lambda2,2)*pp(1,2)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(4,lambda2,2)*pp(1,3)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(3,lambda2,2)*pp(1,4)*O(lambda3,3,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4))))/(4*(WMASS**2))

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( ft8.ne.0d0 ) then  !ft8 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(-2*ft8*(gw**4)*(sw**4)*(-1+SW2)*SW2*(ZMASS**4)*(OP(1,lambda2,2)*OP(2,lambda1,1)*OP(3,lambda4,4)*OP(4,lambda3,3)-OP(3,lambda4,4)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda2,2)-OP(2,lambda4,4)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda3,3)-OP(2,lambda3,3)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*OP(4,lambda1,1)*(OP(2,lambda3,3)*OP(3,lambda2,2)-pp(2,3)*O(lambda2,2,lambda3,3))+pp(1,3)*pp(2,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*OP(3,lambda1,1)*(OP(2,lambda4,4)*OP(4,lambda2,2)-pp(2,4)*O(lambda2,2,lambda4,4))-OP(1,lambda2,2)*OP(2,lambda1,1)*pp(3,4)*O(lambda3,3,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)))/(WMASS**4)

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( ft9.ne.0d0 ) then  !ft9 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) -(ft9*(gw**4)*(sw**4)*(-1+SW2)*SW2*(ZMASS**4)*(OP(1,lambda4,4)*OP(2,lambda3,3)*OP(3,lambda1,1)*OP(4,lambda2,2)+OP(1,lambda4,4)*OP(2,lambda1,1)*OP(3,lambda2,2)*OP(4,lambda3,3)-OP(2,lambda4,4)*OP(4,lambda3,3)*pp(1,3)*O(lambda1,1,lambda2,2)-OP(2,lambda3,3)*OP(3,lambda4,4)*pp(1,4)*O(lambda1,1,lambda2,2)-OP(1,lambda4,4)*OP(4,lambda3,3)*pp(2,3)*O(lambda1,1,lambda2,2)+OP(1,lambda4,4)*OP(2,lambda3,3)*pp(3,4)*O(lambda1,1,lambda2,2)-OP(3,lambda4,4)*OP(4,lambda2,2)*pp(1,2)*O(lambda1,1,lambda3,3)-OP(2,lambda4,4)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda3,3)-OP(1,lambda4,4)*OP(4,lambda2,2)*pp(2,3)*O(lambda1,1,lambda3,3)+OP(1,lambda4,4)*OP(3,lambda2,2)*pp(2,4)*O(lambda1,1,lambda3,3)-OP(3,lambda2,2)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda4,4)-OP(2,lambda3,3)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda4,4)-OP(3,lambda4,4)*OP(4,lambda1,1)*pp(1,2)*O(lambda2,2,lambda3,3)-OP(2,lambda4,4)*OP(4,lambda1,1)*pp(1,3)*O(lambda2,2,lambda3,3)+OP(2,lambda4,4)*OP(3,lambda1,1)*pp(1,4)*O(lambda2,2,lambda3,3)+OP(2,lambda1,1)*OP(3,lambda4,4)*pp(1,4)*O(lambda2,2,lambda3,3)-OP(1,lambda4,4)*OP(3,lambda1,1)*pp(2,4)*O(lambda2,2,lambda3,3)-OP(1,lambda4,4)*OP(2,lambda1,1)*pp(3,4)*O(lambda2,2,lambda3,3)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-OP(3,lambda1,1)*OP(4,lambda3,3)*pp(1,2)*O(lambda2,2,lambda4,4)+OP(2,lambda3,3)*OP(4,lambda1,1)*pp(1,3)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*OP(4,lambda3,3)*pp(1,3)*O(lambda2,2,lambda4,4)-OP(2,lambda3,3)*OP(3,lambda1,1)*pp(1,4)*O(lambda2,2,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(-(OP(3,lambda4,4)*pp(2,4)*O(lambda1,1,lambda2,2))+OP(2,lambda4,4)*(OP(3,lambda2,2)*OP(4,lambda1,1)+pp(3,4)*O(lambda1,1,lambda2,2))+OP(4,lambda2,2)*pp(2,3)*O(lambda1,1,lambda4,4)-OP(3,lambda2,2)*pp(2,4)*O(lambda1,1,lambda4,4)-OP(4,lambda1,1)*pp(2,3)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*(OP(3,lambda4,4)*OP(4,lambda2,2)-pp(3,4)*O(lambda2,2,lambda4,4)))+OP(3,lambda2,2)*OP(4,lambda1,1)*pp(1,2)*O(lambda3,3,lambda4,4)+OP(3,lambda1,1)*OP(4,lambda2,2)*pp(1,2)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(4,lambda2,2)*pp(1,3)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(3,lambda2,2)*pp(1,4)*O(lambda3,3,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+OP(1,lambda2,2)*(OP(3,lambda4,4)*pp(2,4)*O(lambda1,1,lambda3,3)+OP(2,lambda4,4)*(OP(3,lambda1,1)*OP(4,lambda3,3)-pp(3,4)*O(lambda1,1,lambda3,3))+OP(4,lambda3,3)*pp(2,3)*O(lambda1,1,lambda4,4)+OP(2,lambda3,3)*(OP(3,lambda4,4)*OP(4,lambda1,1)-pp(3,4)*O(lambda1,1,lambda4,4))-OP(4,lambda1,1)*pp(2,3)*O(lambda3,3,lambda4,4)-OP(3,lambda1,1)*pp(2,4)*O(lambda3,3,lambda4,4))))/(2*(WMASS**4))

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fm0.ne.0d0 ) then  !fm0 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) -2*fm0*(gw**2)*(sw**2)*(ZMASS**2)*(OP(1,lambda2,2)*OP(2,lambda1,1)-pp(1,2)*O(lambda1,1,lambda2,2))*O(lambda3,3,lambda4,4)

               enddo                      
              enddo
             enddo 
            enddo
         endif


         if( fm1.ne.0d0 ) then  !fm1 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(fm1*(gw**2)*(sw**2)*(ZMASS**2)*(OP(1,lambda2,2)*OP(2,lambda4,4)*O(lambda1,1,lambda3,3)+OP(1,lambda2,2)*OP(2,lambda3,3)*O(lambda1,1,lambda4,4)-pp(1,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*(-(OP(2,lambda3,3)*O(lambda1,1,lambda2,2))+OP(2,lambda1,1)*O(lambda2,2,lambda3,3))-pp(1,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(-(OP(2,lambda4,4)*O(lambda1,1,lambda2,2))+OP(2,lambda1,1)*O(lambda2,2,lambda4,4))))/2

               enddo                      
              enddo
             enddo 
            enddo
         endif


         if( fm2.ne.0d0 ) then  !fm2 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) -(fm2*(gw**2)*(sw**2)*(ZMASS**2)*(OP(1,lambda2,2)*OP(2,lambda1,1)-pp(1,2)*O(lambda1,1,lambda2,2))*O(lambda3,3,lambda4,4))

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fm3.ne.0d0 ) then  !fm3 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(fm3*(gw**2)*(sw**2)*(ZMASS**2)*(OP(1,lambda2,2)*OP(2,lambda4,4)*O(lambda1,1,lambda3,3)+OP(1,lambda2,2)*OP(2,lambda3,3)*O(lambda1,1,lambda4,4)-pp(1,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*(-(OP(2,lambda3,3)*O(lambda1,1,lambda2,2))+OP(2,lambda1,1)*O(lambda2,2,lambda3,3))-pp(1,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(-(OP(2,lambda4,4)*O(lambda1,1,lambda2,2))+OP(2,lambda1,1)*O(lambda2,2,lambda4,4))))/4

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fm4.ne.0d0 ) then  !fm4 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +fm4*(gw**2)*(sw**2)*(ZMASS**2)*(OP(1,lambda2,2)*OP(2,lambda1,1)-pp(1,2)*O(lambda1,1,lambda2,2))*O(lambda3,3,lambda4,4)

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fm5.ne.0d0 ) then  !fm5 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(fm5*(gw**2)*(sw**2)*(ZMASS**2)*(OP(1,lambda2,2)*OP(2,lambda4,4)*O(lambda1,1,lambda3,3)+OP(1,lambda2,2)*OP(2,lambda3,3)*O(lambda1,1,lambda4,4)-pp(1,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*(-(OP(2,lambda3,3)*O(lambda1,1,lambda2,2))+OP(2,lambda1,1)*O(lambda2,2,lambda3,3))-pp(1,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(-(OP(2,lambda4,4)*O(lambda1,1,lambda2,2))+OP(2,lambda1,1)*O(lambda2,2,lambda4,4))))/4

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fm5hc.ne.0d0 ) then  !fm5hc coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(fm5hc*(gw**2)*(sw**2)*(ZMASS**2)*(OP(1,lambda2,2)*OP(2,lambda4,4)*O(lambda1,1,lambda3,3)+OP(1,lambda2,2)*OP(2,lambda3,3)*O(lambda1,1,lambda4,4)-pp(1,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*(-(OP(2,lambda3,3)*O(lambda1,1,lambda2,2))+OP(2,lambda1,1)*O(lambda2,2,lambda3,3))-pp(1,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(-(OP(2,lambda4,4)*O(lambda1,1,lambda2,2))+OP(2,lambda1,1)*O(lambda2,2,lambda4,4))))/4

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fm6.ne.0d0 ) then  !fm6 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) -(fm6*(gw**2)*(sw**2)*(ZMASS**2)*(OP(1,lambda2,2)*OP(2,lambda1,1)-pp(1,2)*O(lambda1,1,lambda2,2))*O(lambda3,3,lambda4,4))

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fm7.ne.0d0 ) then  !fm7 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(fm7*(gw**2)*(sw**2)*(ZMASS**2)*(-(OP(1,lambda2,2)*OP(2,lambda4,4)*O(lambda1,1,lambda3,3))-OP(1,lambda2,2)*OP(2,lambda3,3)*O(lambda1,1,lambda4,4)+pp(1,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*(OP(2,lambda3,3)*O(lambda1,1,lambda2,2)-OP(2,lambda1,1)*O(lambda2,2,lambda3,3))+pp(1,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(OP(2,lambda4,4)*O(lambda1,1,lambda2,2)-OP(2,lambda1,1)*O(lambda2,2,lambda4,4))))/4

               enddo                      
              enddo
             enddo 
            enddo
         endif         

      end








ccccccccccccccc
c aawmwp-vertex c
ccccccccccccccc

      subroutine aawmwp_ONanom4(a1,a2,wm1,wp1,M)
c
c This subroutine calculates the vvvv-vertex with anomalous couplings
c and the on-shell polarization vectors for vector bosons
c
c IN:
c     polv1, polv2      : the incoming polvector wz(bosons in the CM frame)
c     polv3, polv4      : the outgoing polvector wz(bosons in the CM frame)
c OUT:
c     M(lambda1,lambda2) 
c
      IMPLICIT NONE

      INTEGER NGRAPHS
      PARAMETER( NGRAPHS=20 )
      DOUBLE PRECISION rZERO, rONE, rHalf
      PARAMETER( rZERO=0.0D0, rONE=1.0D0, rHALF = 0.5D0 )
      DOUBLE COMPLEX cIMAG, cZERO, cONE
      PARAMETER( cIMAG=(0.0D0,1.0D0), cZERO=(0.0D0,0.0D0), cONE=(1.0D0,0.0D0))
      DOUBLE COMPLEX a1(6,-1:1), a2(6,-1:1), wm1(6,-1:1), wp1(6,-1:1)      
      DOUBLE COMPLEX M(-1:1,-1:1,-1:1,-1:1)
      
      DOUBLE COMPLEX OFFeps(6,-1:1,NGRAPHS), OFFepsC(6,-1:1,NGRAPHS)
      DOUBLE PRECISION K(0:3,4)      
      DOUBLE COMPLEX O(-1:1,4,-1:1,4), OP(4,-1:1,4)
      DOUBLE COMPLEXPP(1:4,1:4)
    
      INTEGER I,J,N,MU,NU
      INTEGER LAMBDA1, LAMBDA2, LAMBDA3, LAMBDA4, LAMBDA
      INTEGER L1,L2

      DOUBLE COMPLEX DOTCC, DOTRC,contract_Tjj
      DOUBLE PRECISION DOTRR
      EXTERNAL DOTCC, DOTRC, DOTRR, contract_Tjj
      double precision sw2, sw
#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "gmunu.inc"
#include "transm.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"
      sw2 = sin2w
      sw = sqrt(sin2w)
      OFFeps(6,-1:1,4) = cZero
      K(0:3,1:4) = rZero
      M(-1:1,-1:1,-1:1,-1:1) = cZero
      PP(1:4,1:4) = rZero
      O(-1:1,1:4,-1:1,1:4) = cZero
      OP(1:4,-1:1,1:4) = cZero  
      
      OFFeps(1:6,-1:1,1) = a1(1:6,-1:1)
      OFFeps(1:6,-1:1,2) = a2(1:6,-1:1)
      OFFeps(1:6,-1:1,3) = wm1(1:6,-1:1)
      OFFeps(1:6,-1:1,4) = wp1(1:6,-1:1) 

      do i = 1,4
        K(0,i) = dble(OFFeps(5,-1,i))*(-rONE) !K1(0:3)
        K(1,i) = dble(OFFeps(6,-1,i))*(-rONE) !K2(0:3)
        K(2,i) = dimag(OFFeps(6,-1,i))*(-rONE)!K3(0:3)
        K(3,i) = dimag(OFFeps(5,-1,i))*(-rONE)!K4(0:3)
      enddo
      



      do i = 1,4
        do j = 1,4
          pp(i,j) = dotrr(k(0,i),k(0,j))
        enddo
      enddo

      do L1 = -1,1
        do L2 = -1,1 
          do i = 1,4
            do j = 1,4
              O(L1,i,L2,j) = dotcc(OFFeps(1,L1,i),OFFeps(1,L2,j))
            enddo
          enddo
        enddo
      enddo

      do LAMBDA = -1,1
        do i = 1,4 
          do j = 1,4
            OP(j,lambda,i) = dotrc(k(0,j),OFFeps(1,lambda,i))
          enddo
        enddo
      enddo


!!!! ===============================================================      
!!!! COUPLINGS CALCULATION
!!!! ===============================================================



                 

         if( ft0.ne.0d0 ) then  !ft0 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +8*ft0*(gw**4)*SW2*(OP(1,lambda2,2)*OP(2,lambda1,1)-pp(1,2)*O(lambda1,1,lambda2,2))*(OP(3,lambda4,4)*OP(4,lambda3,3)-pp(3,4)*O(lambda3,3,lambda4,4))

               enddo                      
              enddo
             enddo 
            enddo
         endif


         if( ft1.ne.0d0 ) then  !ft1 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +4*ft1*(gw**4)*SW2*(-(OP(2,lambda4,4)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda3,3))-OP(2,lambda3,3)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*OP(4,lambda1,1)*(OP(2,lambda3,3)*OP(3,lambda2,2)-pp(2,3)*O(lambda2,2,lambda3,3))+pp(1,3)*pp(2,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*OP(3,lambda1,1)*(OP(2,lambda4,4)*OP(4,lambda2,2)-pp(2,4)*O(lambda2,2,lambda4,4)))

               enddo                      
              enddo
             enddo 
            enddo
         endif


         if( ft2.ne.0d0 ) then  !ft2 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +ft2*(gw**4)*SW2*(OP(1,lambda4,4)*OP(2,lambda1,1)*OP(3,lambda2,2)*OP(4,lambda3,3)-OP(2,lambda4,4)*OP(4,lambda3,3)*pp(1,3)*O(lambda1,1,lambda2,2)-OP(2,lambda3,3)*OP(3,lambda4,4)*pp(1,4)*O(lambda1,1,lambda2,2)-OP(1,lambda4,4)*OP(4,lambda3,3)*pp(2,3)*O(lambda1,1,lambda2,2)+OP(1,lambda4,4)*OP(2,lambda3,3)*pp(3,4)*O(lambda1,1,lambda2,2)-OP(3,lambda4,4)*OP(4,lambda2,2)*pp(1,2)*O(lambda1,1,lambda3,3)-OP(3,lambda2,2)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda4,4)-OP(3,lambda4,4)*OP(4,lambda1,1)*pp(1,2)*O(lambda2,2,lambda3,3)+OP(2,lambda1,1)*OP(3,lambda4,4)*pp(1,4)*O(lambda2,2,lambda3,3)-OP(1,lambda4,4)*OP(2,lambda1,1)*pp(3,4)*O(lambda2,2,lambda3,3)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-OP(3,lambda1,1)*OP(4,lambda3,3)*pp(1,2)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*OP(4,lambda3,3)*pp(1,3)*O(lambda2,2,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*((-(OP(3,lambda4,4)*pp(2,4))+OP(2,lambda4,4)*pp(3,4))*O(lambda1,1,lambda2,2)+OP(2,lambda1,1)*(OP(3,lambda4,4)*OP(4,lambda2,2)-pp(3,4)*O(lambda2,2,lambda4,4)))+OP(3,lambda2,2)*OP(4,lambda1,1)*pp(1,2)*O(lambda3,3,lambda4,4)+OP(3,lambda1,1)*OP(4,lambda2,2)*pp(1,2)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(4,lambda2,2)*pp(1,3)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(3,lambda2,2)*pp(1,4)*O(lambda3,3,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+OP(1,lambda2,2)*(OP(3,lambda4,4)*pp(2,4)*O(lambda1,1,lambda3,3)+OP(2,lambda4,4)*(OP(3,lambda1,1)*OP(4,lambda3,3)-pp(3,4)*O(lambda1,1,lambda3,3))+OP(4,lambda3,3)*pp(2,3)*O(lambda1,1,lambda4,4)+OP(2,lambda3,3)*(OP(3,lambda4,4)*OP(4,lambda1,1)-pp(3,4)*O(lambda1,1,lambda4,4))-OP(4,lambda1,1)*pp(2,3)*O(lambda3,3,lambda4,4)-OP(3,lambda1,1)*pp(2,4)*O(lambda3,3,lambda4,4)))

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( ft3.ne.0d0 ) then  !ft3 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +2*ft3*gw**4*sw**2*(pp(1,4)*OP(3,lambda2,2)*(-OP(2,lambda4,4))*O(lambda1,1,lambda3,3)-pp(1,3)*OP(4,lambda2,2)*OP(2,lambda3,3)*O(lambda1,1,lambda4,4)-pp(1,3)*OP(4,lambda1,1)*OP(2,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,4)*OP(3,lambda1,1)*OP(2,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*(-pp(2,3)*OP(4,lambda2,2)*O(lambda1,1,lambda3,3)+pp(2,4)*(OP(3,lambda2,2)*O(lambda1,1,lambda3,3)-OP(3,lambda1,1)*O(lambda2,2,lambda3,3))+OP(3,lambda1,1)*OP(4,lambda2,2)*OP(2,lambda3,3))+pp(1,3)*OP(4,lambda1,1)*OP(2,lambda3,3)*O(lambda2,2,lambda4,4)-pp(1,4)*OP(3,lambda1,1)*OP(2,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(-pp(2,3)*OP(4,lambda1,1)*O(lambda2,2,lambda4,4)+pp(2,3)*OP(4,lambda2,2)*O(lambda1,1,lambda4,4)-pp(2,4)*OP(3,lambda2,2)*O(lambda1,1,lambda4,4)+OP(4,lambda1,1)*OP(3,lambda2,2)*OP(2,lambda4,4))+pp(1,4)*pp(2,3)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3))

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( ft4.ne.0d0 ) then  !ft4 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(gw**4*ft4*sw**2*(-(pp(1, 4)*OP(2, lambda4, 4)*OP(3, lambda2, 2)*O(lambda1, 1, lambda3, 3)) - pp(1, 3)*OP(2, lambda3, 3)*OP(4, lambda2, 2)*O(lambda1, 1, lambda4, 4) + pp(1, 4)*OP(2, lambda4, 4)*OP(3, lambda1, 1)*O(lambda2, 2, lambda3, 3) - pp(1, 3)*OP(2, lambda4, 4)*OP(4, lambda1, 1)*O(lambda2, 2, lambda3, 3) + pp(1, 3)*pp(2, 4)*O(lambda1, 1, lambda4, 4)*O(lambda2, 2, lambda3, 3) + OP(1, lambda4, 4)*(OP(2, lambda3, 3)*OP(3, lambda1, 1)*OP(4, lambda2, 2) - pp(2, 3)*OP(4, lambda2, 2)*O(lambda1, 1, lambda3, 3) +   pp(2, 4)*(OP(3, lambda2, 2)*O(lambda1, 1, lambda3, 3) - OP(3, lambda1, 1)*O(lambda2, 2, lambda3, 3))) - pp(1, 4)*OP(2, lambda3, 3)*OP(3, lambda1, 1)*O(lambda2, 2, lambda4, 4) + pp(1, 3)*OP(2, lambda3, 3)*OP(4, lambda1, 1)*O(lambda2, 2, lambda4, 4) + pp(1, 4)*pp(2, 3)*O(lambda1, 1, lambda3, 3)*O(lambda2, 2, lambda4, 4) + OP(1, lambda3, 3)*(OP(2, lambda4, 4)*OP(3, lambda2, 2)*OP(4, lambda1, 1) - pp(2, 4)*OP(3, lambda2, 2)*O(lambda1, 1, lambda4, 4) + pp(2, 3)*OP(4, lambda2, 2)*O(lambda1, 1, lambda4, 4) -   pp(2, 3)*OP(4, lambda1, 1)*O(lambda2, 2, lambda4, 4))))/2

               enddo                      
              enddo
             enddo 
            enddo
         endif


         if( ft5.ne.0d0 ) then  !ft5 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(-2*ft5*(gw**4)*SW2*(-1+SW2)*(ZMASS**2)*(OP(1,lambda2,2)*OP(2,lambda1,1)-pp(1,2)*O(lambda1,1,lambda2,2))*(OP(3,lambda4,4)*OP(4,lambda3,3)-pp(3,4)*O(lambda3,3,lambda4,4)))/(WMASS**2)

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( ft6.ne.0d0 ) then  !ft6 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(ft6*(gw**4)*SW2*(-1+SW2)*(ZMASS**2)*(OP(2,lambda4,4)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda3,3)+OP(2,lambda3,3)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda4,4)-pp(1,4)*pp(2,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*OP(4,lambda1,1)*(-(OP(2,lambda3,3)*OP(3,lambda2,2))+pp(2,3)*O(lambda2,2,lambda3,3))-pp(1,3)*pp(2,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*OP(3,lambda1,1)*(-(OP(2,lambda4,4)*OP(4,lambda2,2))+pp(2,4)*O(lambda2,2,lambda4,4))))/(WMASS**2)

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( ft7.ne.0d0 ) then  !ft7 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) -(ft7*(gw**4)*SW2*(-1+SW2)*(ZMASS**2)*(OP(1,lambda4,4)*OP(2,lambda1,1)*OP(3,lambda2,2)*OP(4,lambda3,3)-OP(2,lambda4,4)*OP(4,lambda3,3)*pp(1,3)*O(lambda1,1,lambda2,2)-OP(2,lambda3,3)*OP(3,lambda4,4)*pp(1,4)*O(lambda1,1,lambda2,2)-OP(1,lambda4,4)*OP(4,lambda3,3)*pp(2,3)*O(lambda1,1,lambda2,2)+OP(1,lambda4,4)*OP(2,lambda3,3)*pp(3,4)*O(lambda1,1,lambda2,2)-OP(3,lambda4,4)*OP(4,lambda2,2)*pp(1,2)*O(lambda1,1,lambda3,3)-OP(3,lambda2,2)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda4,4)-OP(3,lambda4,4)*OP(4,lambda1,1)*pp(1,2)*O(lambda2,2,lambda3,3)+OP(2,lambda1,1)*OP(3,lambda4,4)*pp(1,4)*O(lambda2,2,lambda3,3)-OP(1,lambda4,4)*OP(2,lambda1,1)*pp(3,4)*O(lambda2,2,lambda3,3)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-OP(3,lambda1,1)*OP(4,lambda3,3)*pp(1,2)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*OP(4,lambda3,3)*pp(1,3)*O(lambda2,2,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*((-(OP(3,lambda4,4)*pp(2,4))+OP(2,lambda4,4)*pp(3,4))*O(lambda1,1,lambda2,2)+OP(2,lambda1,1)*(OP(3,lambda4,4)*OP(4,lambda2,2)-pp(3,4)*O(lambda2,2,lambda4,4)))+OP(3,lambda2,2)*OP(4,lambda1,1)*pp(1,2)*O(lambda3,3,lambda4,4)+OP(3,lambda1,1)*OP(4,lambda2,2)*pp(1,2)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(4,lambda2,2)*pp(1,3)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(3,lambda2,2)*pp(1,4)*O(lambda3,3,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+OP(1,lambda2,2)*(OP(3,lambda4,4)*pp(2,4)*O(lambda1,1,lambda3,3)+OP(2,lambda4,4)*(OP(3,lambda1,1)*OP(4,lambda3,3)-pp(3,4)*O(lambda1,1,lambda3,3))+OP(4,lambda3,3)*pp(2,3)*O(lambda1,1,lambda4,4)+OP(2,lambda3,3)*(OP(3,lambda4,4)*OP(4,lambda1,1)-pp(3,4)*O(lambda1,1,lambda4,4))-OP(4,lambda1,1)*pp(2,3)*O(lambda3,3,lambda4,4)-OP(3,lambda1,1)*pp(2,4)*O(lambda3,3,lambda4,4))))/(4*(WMASS**2))

               enddo                      
              enddo
             enddo 
            enddo
         endif


         if( fm0.ne.0d0 ) then  !fm0 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) -2*fm0*(gw**2)*(sw**2)*(WMASS**2)*(OP(1,lambda2,2)*OP(2,lambda1,1)-pp(1,2)*O(lambda1,1,lambda2,2))*O(lambda3,3,lambda4,4)

               enddo                      
              enddo
             enddo 
            enddo
         endif


         if( fm1.ne.0d0 ) then  !fm1 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(fm1*(gw**2)*(sw**2)*(WMASS**2)*(OP(1,lambda2,2)*OP(2,lambda4,4)*O(lambda1,1,lambda3,3)+OP(1,lambda2,2)*OP(2,lambda3,3)*O(lambda1,1,lambda4,4)-pp(1,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*(-(OP(2,lambda3,3)*O(lambda1,1,lambda2,2))+OP(2,lambda1,1)*O(lambda2,2,lambda3,3))-pp(1,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(-(OP(2,lambda4,4)*O(lambda1,1,lambda2,2))+OP(2,lambda1,1)*O(lambda2,2,lambda4,4))))/2

               enddo                      
              enddo
             enddo 
            enddo
         endif


         if( fm2.ne.0d0 ) then  !fm2 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +fm2*(gw**2)*(sw**2)*(-1+SW2)*(ZMASS**2)*(OP(1,lambda2,2)*OP(2,lambda1,1)-pp(1,2)*O(lambda1,1,lambda2,2))*O(lambda3,3,lambda4,4)

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fm3.ne.0d0 ) then  !fm3 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(fm3*(gw**2)*(sw**2)*(-1+SW2)*(ZMASS**2)*(-(OP(1,lambda2,2)*OP(2,lambda4,4)*O(lambda1,1,lambda3,3))-OP(1,lambda2,2)*OP(2,lambda3,3)*O(lambda1,1,lambda4,4)+pp(1,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*(OP(2,lambda3,3)*O(lambda1,1,lambda2,2)-OP(2,lambda1,1)*O(lambda2,2,lambda3,3))+pp(1,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(OP(2,lambda4,4)*O(lambda1,1,lambda2,2)-OP(2,lambda1,1)*O(lambda2,2,lambda4,4))))/4

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fm4.ne.0d0 ) then  !fm4 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) -(fm4*(gw**2)*(sw**2)*(WMASS**2)*(OP(1,lambda2,2)*OP(2,lambda1,1)-pp(1,2)*O(lambda1,1,lambda2,2))*O(lambda3,3,lambda4,4))

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fm5.ne.0d0 ) then  !fm5 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(fm5*(gw**2)*(sw**2)*(WMASS**2)*(-(OP(1,lambda2,2)*OP(2,lambda4,4)*O(lambda1,1,lambda3,3))-OP(1,lambda2,2)*OP(2,lambda3,3)*O(lambda1,1,lambda4,4)+pp(1,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*(OP(2,lambda3,3)*O(lambda1,1,lambda2,2)-OP(2,lambda1,1)*O(lambda2,2,lambda3,3))+pp(1,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(OP(2,lambda4,4)*O(lambda1,1,lambda2,2)-OP(2,lambda1,1)*O(lambda2,2,lambda4,4))))/4

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fm5hc.ne.0d0 ) then  !fm5hc coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(fm5hc*(gw**2)*(sw**2)*(WMASS**2)*(-(OP(1,lambda2,2)*OP(2,lambda4,4)*O(lambda1,1,lambda3,3))-OP(1,lambda2,2)*OP(2,lambda3,3)*O(lambda1,1,lambda4,4)+pp(1,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*(OP(2,lambda3,3)*O(lambda1,1,lambda2,2)-OP(2,lambda1,1)*O(lambda2,2,lambda3,3))+pp(1,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(OP(2,lambda4,4)*O(lambda1,1,lambda2,2)-OP(2,lambda1,1)*O(lambda2,2,lambda4,4))))/4

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fm6.ne.0d0 ) then  !fm6 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) -(fm6*(gw**2)*(sw**2)*(WMASS**2)*(OP(1,lambda2,2)*OP(2,lambda1,1)-pp(1,2)*O(lambda1,1,lambda2,2))*O(lambda3,3,lambda4,4))

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fm7.ne.0d0 ) then  !fm7 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(fm7*(gw**2)*(sw**2)*(WMASS**2)*(-(OP(1,lambda2,2)*OP(2,lambda4,4)*O(lambda1,1,lambda3,3))-OP(1,lambda2,2)*OP(2,lambda3,3)*O(lambda1,1,lambda4,4)+pp(1,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*(OP(2,lambda3,3)*O(lambda1,1,lambda2,2)-OP(2,lambda1,1)*O(lambda2,2,lambda3,3))+pp(1,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(OP(2,lambda4,4)*O(lambda1,1,lambda2,2)-OP(2,lambda1,1)*O(lambda2,2,lambda4,4))))/4

               enddo                      
              enddo
             enddo 
            enddo
         endif         

      end






ccccccccccccccc
c ahhh-vertex c
ccccccccccccccc

      subroutine ahhh_ONanom4(a1,h1,h2,h3,M)
c
c This subroutine calculates the vvvv-vertex with anomalous couplings
c and the on-shell polarization vectors for vector bosons
c
c IN:
c     polv1, polv2      : the incoming polvector wz(bosons in the CM frame)
c     polv3, polv4      : the outgoing polvector wz(bosons in the CM frame)
c OUT:
c     M(lambda1,lambda2) 
c
      IMPLICIT NONE

      INTEGER NGRAPHS
      PARAMETER( NGRAPHS=20 )
      DOUBLE PRECISION rZERO, rONE, rHalf
      PARAMETER( rZERO=0.0D0, rONE=1.0D0, rHALF = 0.5D0 )
      DOUBLE COMPLEX cIMAG, cZERO, cONE
      PARAMETER( cIMAG=(0.0D0,1.0D0), cZERO=(0.0D0,0.0D0), cONE=(1.0D0,0.0D0))
      DOUBLE COMPLEX a1(6,-1:1), h1(6,-1:1), h2(6,-1:1), h3(6,-1:1)      
      DOUBLE COMPLEX M(-1:1,-1:1,-1:1,-1:1)
      
      DOUBLE COMPLEX OFFeps(6,-1:1,NGRAPHS), OFFepsC(6,-1:1,NGRAPHS)
      DOUBLE PRECISION K(0:3,4)      
      DOUBLE COMPLEX O(-1:1,4,-1:1,4), OP(4,-1:1,4)
      DOUBLE COMPLEXPP(1:4,1:4)
    
      INTEGER I,J,N,MU,NU
      INTEGER LAMBDA1, LAMBDA2, LAMBDA3, LAMBDA4, LAMBDA
      INTEGER L1,L2

      DOUBLE COMPLEX DOTCC, DOTRC,contract_Tjj
      DOUBLE PRECISION DOTRR
      EXTERNAL DOTCC, DOTRC, DOTRR, contract_Tjj
      double precision sw2, sw
#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "gmunu.inc"
#include "transm.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"
      sw2 = sin2w
      sw = sqrt(sin2w)
      OFFeps(6,-1:1,4) = cZero
      K(0:3,1:4) = rZero
      M(-1:1,-1:1,-1:1,-1:1) = cZero
      PP(1:4,1:4) = rZero
      O(-1:1,1:4,-1:1,1:4) = cZero
      OP(1:4,-1:1,1:4) = cZero  
      
      OFFeps(1:6,-1:1,1) = a1(1:6,-1:1)
      OFFeps(1:6,-1:1,2) = h1(1:6,-1:1)
      OFFeps(1:6,-1:1,3) = h2(1:6,-1:1)
      OFFeps(1:6,-1:1,4) = h3(1:6,-1:1) 

      do i = 1,4
        K(0,i) = dble(OFFeps(5,-1,i))*(-rONE) !K1(0:3)
        K(1,i) = dble(OFFeps(6,-1,i))*(-rONE) !K2(0:3)
        K(2,i) = dimag(OFFeps(6,-1,i))*(-rONE)!K3(0:3)
        K(3,i) = dimag(OFFeps(5,-1,i))*(-rONE)!K4(0:3)
      enddo
      



      do i = 1,4
        do j = 1,4
          pp(i,j) = dotrr(k(0,i),k(0,j))
        enddo
      enddo

      do L1 = -1,1
        do L2 = -1,1 
          do i = 1,4
            do j = 1,4
              O(L1,i,L2,j) = dotcc(OFFeps(1,L1,i),OFFeps(1,L2,j))
            enddo
          enddo
        enddo
      enddo

      do LAMBDA = -1,1
        do i = 1,4 
          do j = 1,4
            OP(j,lambda,i) = dotrc(k(0,j),OFFeps(1,lambda,i))
          enddo
        enddo
      enddo


!!!! ===============================================================      
!!!! COUPLINGS CALCULATION
!!!! ===============================================================


        

      end


ccccccccccccccc
c ahhz-vertex c
ccccccccccccccc

      subroutine ahhz_ONanom4(a1,h1,h2,z,M)
c
c This subroutine calculates the vvvv-vertex with anomalous couplings
c and the on-shell polarization vectors for vector bosons
c
c IN:
c     polv1, polv2      : the incoming polvector wz(bosons in the CM frame)
c     polv3, polv4      : the outgoing polvector wz(bosons in the CM frame)
c OUT:
c     M(lambda1,lambda2) 
c
      IMPLICIT NONE

      INTEGER NGRAPHS
      PARAMETER( NGRAPHS=20 )
      DOUBLE PRECISION rZERO, rONE, rHalf
      PARAMETER( rZERO=0.0D0, rONE=1.0D0, rHALF = 0.5D0 )
      DOUBLE COMPLEX cIMAG, cZERO, cONE
      PARAMETER( cIMAG=(0.0D0,1.0D0), cZERO=(0.0D0,0.0D0), cONE=(1.0D0,0.0D0))
      DOUBLE COMPLEX a1(6,-1:1), h1(6,-1:1), h2(6,-1:1), z(6,-1:1)      
      DOUBLE COMPLEX M(-1:1,-1:1,-1:1,-1:1)
      
      DOUBLE COMPLEX OFFeps(6,-1:1,NGRAPHS), OFFepsC(6,-1:1,NGRAPHS)
      DOUBLE PRECISION K(0:3,4)      
      DOUBLE COMPLEX O(-1:1,4,-1:1,4), OP(4,-1:1,4)
      DOUBLE COMPLEXPP(1:4,1:4)
    
      INTEGER I,J,N,MU,NU
      INTEGER LAMBDA1, LAMBDA2, LAMBDA3, LAMBDA4, LAMBDA
      INTEGER L1,L2

      DOUBLE COMPLEX DOTCC, DOTRC,contract_Tjj
      DOUBLE PRECISION DOTRR
      EXTERNAL DOTCC, DOTRC, DOTRR, contract_Tjj
      double precision sw2, sw
#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "gmunu.inc"
#include "transm.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"
      sw2 = sin2w
      sw = sqrt(sin2w)
      OFFeps(6,-1:1,4) = cZero
      K(0:3,1:4) = rZero
      M(-1:1,-1:1,-1:1,-1:1) = cZero
      PP(1:4,1:4) = rZero
      O(-1:1,1:4,-1:1,1:4) = cZero
      OP(1:4,-1:1,1:4) = cZero  
      
      OFFeps(1:6,-1:1,1) = a1(1:6,-1:1)
      OFFeps(1:6,-1:1,2) = h1(1:6,-1:1)
      OFFeps(1:6,-1:1,3) = h2(1:6,-1:1)
      OFFeps(1:6,-1:1,4) = z(1:6,-1:1) 

      do i = 1,4
        K(0,i) = dble(OFFeps(5,-1,i))*(-rONE) !K1(0:3)
        K(1,i) = dble(OFFeps(6,-1,i))*(-rONE) !K2(0:3)
        K(2,i) = dimag(OFFeps(6,-1,i))*(-rONE)!K3(0:3)
        K(3,i) = dimag(OFFeps(5,-1,i))*(-rONE)!K4(0:3)
      enddo
      



      do i = 1,4
        do j = 1,4
          pp(i,j) = dotrr(k(0,i),k(0,j))
        enddo
      enddo

      do L1 = -1,1
        do L2 = -1,1 
          do i = 1,4
            do j = 1,4
              O(L1,i,L2,j) = dotcc(OFFeps(1,L1,i),OFFeps(1,L2,j))
            enddo
          enddo
        enddo
      enddo

      do LAMBDA = -1,1
        do i = 1,4 
          do j = 1,4
            OP(j,lambda,i) = dotrc(k(0,j),OFFeps(1,lambda,i))
          enddo
        enddo
      enddo


!!!! ===============================================================      
!!!! COUPLINGS CALCULATION
!!!! ===============================================================



         if( fm0.ne.0d0 ) then  !fm0 coupling
            do lambda1 = -1,1
!             do lambda2 = -1,1
!              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,0,0,lambda4) = M(lambda1,0,0,lambda4) +(2*fm0*(gw**2)*sw*WMASS*1.0*1.0*pp(2,3)*(OP(1,lambda4,4)*OP(4,lambda1,1)-pp(1,4)*O(lambda1,1,lambda4,4)))/ZMASS

               enddo                      
!              enddo
!             enddo 
            enddo
         endif


         if( fm1.ne.0d0 ) then  !fm1 coupling
            do lambda1 = -1,1
!             do lambda2 = -1,1
!              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,0,0,lambda4) = M(lambda1,0,0,lambda4) -(fm1*(gw**2)*sw*WMASS*1.0*1.0*(OP(3,lambda4,4)*(OP(4,lambda1,1)*pp(1,2)-OP(2,lambda1,1)*pp(1,4))+OP(2,lambda4,4)*(OP(4,lambda1,1)*pp(1,3)-OP(3,lambda1,1)*pp(1,4))+OP(1,lambda4,4)*OP(3,lambda1,1)*pp(2,4)+OP(1,lambda4,4)*OP(2,lambda1,1)*pp(3,4)-pp(1,3)*pp(2,4)*O(lambda1,1,lambda4,4)-pp(1,2)*pp(3,4)*O(lambda1,1,lambda4,4)))/(2*ZMASS)

               enddo                      
!              enddo
!             enddo 
            enddo
         endif


         if( fm2.ne.0d0 ) then  !fm2 coupling
            do lambda1 = -1,1
!             do lambda2 = -1,1
!              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,0,0,lambda4) = M(lambda1,0,0,lambda4) -((fm2*(gw**2)*(sw**3)*ZMASS*1.0*1.0*pp(2,3)*(OP(1,lambda4,4)*OP(4,lambda1,1)-pp(1,4)*O(lambda1,1,lambda4,4)))/WMASS)

               enddo                      
!              enddo
!             enddo 
            enddo
         endif

         if( fm3.ne.0d0 ) then  !fm3 coupling
            do lambda1 = -1,1
!             do lambda2 = -1,1
!              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,0,0,lambda4) = M(lambda1,0,0,lambda4) +(fm3*(gw**2)*(sw**3)*ZMASS*1.0*1.0*(OP(3,lambda4,4)*(OP(4,lambda1,1)*pp(1,2)-OP(2,lambda1,1)*pp(1,4))+OP(2,lambda4,4)*(OP(4,lambda1,1)*pp(1,3)-OP(3,lambda1,1)*pp(1,4))+OP(1,lambda4,4)*OP(3,lambda1,1)*pp(2,4)+OP(1,lambda4,4)*OP(2,lambda1,1)*pp(3,4)-pp(1,3)*pp(2,4)*O(lambda1,1,lambda4,4)-pp(1,2)*pp(3,4)*O(lambda1,1,lambda4,4)))/(4*WMASS)

               enddo                      
!              enddo
!             enddo 
            enddo
         endif

         if( fm4.ne.0d0 ) then  !fm4 coupling
            do lambda1 = -1,1
!             do lambda2 = -1,1
!              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,0,0,lambda4) = M(lambda1,0,0,lambda4) +(fm4*(gw**2)*sw*(-1+2*SW2)*ZMASS*1.0*1.0*pp(2,3)*(OP(1,lambda4,4)*OP(4,lambda1,1)-pp(1,4)*O(lambda1,1,lambda4,4)))/(2*WMASS)

               enddo                      
!              enddo
!             enddo 
            enddo
         endif

         if( fm5.ne.0d0 ) then  !fm5 coupling
            do lambda1 = -1,1
!             do lambda2 = -1,1
!              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,0,0,lambda4) = M(lambda1,0,0,lambda4) +(fm5*(gw**2)*sw*(-1+2*SW2)*ZMASS*1.0*1.0*(OP(3,lambda4,4)*(OP(4,lambda1,1)*pp(1,2)-OP(2,lambda1,1)*pp(1,4))+OP(2,lambda4,4)*(OP(4,lambda1,1)*pp(1,3)-OP(3,lambda1,1)*pp(1,4))+OP(1,lambda4,4)*OP(3,lambda1,1)*pp(2,4)+OP(1,lambda4,4)*OP(2,lambda1,1)*pp(3,4)-pp(1,3)*pp(2,4)*O(lambda1,1,lambda4,4)-pp(1,2)*pp(3,4)*O(lambda1,1,lambda4,4)))/(8*WMASS)

               enddo                      
!              enddo
!             enddo 
            enddo
         endif

         if( fm5hc.ne.0d0 ) then  !fm5hc coupling
            do lambda1 = -1,1
!             do lambda2 = -1,1
!              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,0,0,lambda4) = M(lambda1,0,0,lambda4) +(fm5hc*(gw**2)*sw*(-1+2*SW2)*ZMASS*1.0*1.0*(OP(3,lambda4,4)*(OP(4,lambda1,1)*pp(1,2)-OP(2,lambda1,1)*pp(1,4))+OP(2,lambda4,4)*(OP(4,lambda1,1)*pp(1,3)-OP(3,lambda1,1)*pp(1,4))+OP(1,lambda4,4)*OP(3,lambda1,1)*pp(2,4)+OP(1,lambda4,4)*OP(2,lambda1,1)*pp(3,4)-pp(1,3)*pp(2,4)*O(lambda1,1,lambda4,4)-pp(1,2)*pp(3,4)*O(lambda1,1,lambda4,4)))/(8*WMASS)

               enddo                      
!              enddo
!             enddo 
            enddo
         endif

         if( fm6.ne.0d0 ) then  !fm6 coupling
            do lambda1 = -1,1
!             do lambda2 = -1,1
!              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,0,0,lambda4) = M(lambda1,0,0,lambda4) +(fm6*(gw**2)*sw*WMASS*1.0*1.0*pp(2,3)*(OP(1,lambda4,4)*OP(4,lambda1,1)-pp(1,4)*O(lambda1,1,lambda4,4)))/ZMASS

               enddo                      
!              enddo
!             enddo 
            enddo
         endif

         if( fm7.ne.0d0 ) then  !fm7 coupling
            do lambda1 = -1,1
!             do lambda2 = -1,1
!              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,0,0,lambda4) = M(lambda1,0,0,lambda4) +(fm7*(gw**2)*sw*WMASS*1.0*1.0*(OP(3,lambda4,4)*(OP(4,lambda1,1)*pp(1,2)-OP(2,lambda1,1)*pp(1,4))+OP(2,lambda4,4)*(OP(4,lambda1,1)*pp(1,3)-OP(3,lambda1,1)*pp(1,4))+OP(1,lambda4,4)*OP(3,lambda1,1)*pp(2,4)+OP(1,lambda4,4)*OP(2,lambda1,1)*pp(3,4)-pp(1,3)*pp(2,4)*O(lambda1,1,lambda4,4)-pp(1,2)*pp(3,4)*O(lambda1,1,lambda4,4)))/(4*ZMASS)

               enddo                      
!              enddo
!             enddo 
            enddo
         endif         

      end



ccccccccccccccc
c ahzz-vertex c
ccccccccccccccc

      subroutine ahzz_ONanom4(a1,h,z1,z2,M)
c
c This subroutine calculates the vvvv-vertex with anomalous couplings
c and the on-shell polarization vectors for vector bosons
c
c IN:
c     polv1, polv2      : the incoming polvector wz(bosons in the CM frame)
c     polv3, polv4      : the outgoing polvector wz(bosons in the CM frame)
c OUT:
c     M(lambda1,lambda2) 
c
      IMPLICIT NONE

      INTEGER NGRAPHS
      PARAMETER( NGRAPHS=20 )
      DOUBLE PRECISION rZERO, rONE, rHalf
      PARAMETER( rZERO=0.0D0, rONE=1.0D0, rHALF = 0.5D0 )
      DOUBLE COMPLEX cIMAG, cZERO, cONE
      PARAMETER( cIMAG=(0.0D0,1.0D0), cZERO=(0.0D0,0.0D0), cONE=(1.0D0,0.0D0))
      DOUBLE COMPLEX a1(6,-1:1), h(6,-1:1), z1(6,-1:1), z2(6,-1:1)      
      DOUBLE COMPLEX M(-1:1,-1:1,-1:1,-1:1)
      
      DOUBLE COMPLEX OFFeps(6,-1:1,NGRAPHS), OFFepsC(6,-1:1,NGRAPHS)
      DOUBLE PRECISION K(0:3,4)      
      DOUBLE COMPLEX O(-1:1,4,-1:1,4), OP(4,-1:1,4)
      DOUBLE COMPLEXPP(1:4,1:4)
    
      INTEGER I,J,N,MU,NU
      INTEGER LAMBDA1, LAMBDA2, LAMBDA3, LAMBDA4, LAMBDA
      INTEGER L1,L2

      DOUBLE COMPLEX DOTCC, DOTRC,contract_Tjj
      DOUBLE PRECISION DOTRR
      EXTERNAL DOTCC, DOTRC, DOTRR, contract_Tjj
      double precision sw2, sw
#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "gmunu.inc"
#include "transm.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"
      sw2 = sin2w
      sw = sqrt(sin2w)
      OFFeps(6,-1:1,4) = cZero
      K(0:3,1:4) = rZero
      M(-1:1,-1:1,-1:1,-1:1) = cZero
      PP(1:4,1:4) = rZero
      O(-1:1,1:4,-1:1,1:4) = cZero
      OP(1:4,-1:1,1:4) = cZero  
      
      OFFeps(1:6,-1:1,1) = a1(1:6,-1:1)
      OFFeps(1:6,-1:1,2) = h(1:6,-1:1)
      OFFeps(1:6,-1:1,3) = z1(1:6,-1:1)
      OFFeps(1:6,-1:1,4) = z2(1:6,-1:1) 

      do i = 1,4
        K(0,i) = dble(OFFeps(5,-1,i))*(-rONE) !K1(0:3)
        K(1,i) = dble(OFFeps(6,-1,i))*(-rONE) !K2(0:3)
        K(2,i) = dimag(OFFeps(6,-1,i))*(-rONE)!K3(0:3)
        K(3,i) = dimag(OFFeps(5,-1,i))*(-rONE)!K4(0:3)
      enddo
      



      do i = 1,4
        do j = 1,4
          pp(i,j) = dotrr(k(0,i),k(0,j))
        enddo
      enddo

      do L1 = -1,1
        do L2 = -1,1 
          do i = 1,4
            do j = 1,4
              O(L1,i,L2,j) = dotcc(OFFeps(1,L1,i),OFFeps(1,L2,j))
            enddo
          enddo
        enddo
      enddo

      do LAMBDA = -1,1
        do i = 1,4 
          do j = 1,4
            OP(j,lambda,i) = dotrc(k(0,j),OFFeps(1,lambda,i))
          enddo
        enddo
      enddo


!!!! ===============================================================      
!!!! COUPLINGS CALCULATION
!!!! ===============================================================



         if( fm5.ne.0d0 ) then  !fm5 coupling
            do lambda1 = -1,1
!             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,0,lambda3,lambda4) = M(lambda1,0,lambda3,lambda4) -(fm5*(gw**2)*sw*(ZMASS**2)*1.0*(OP(3,lambda4,4)*pp(1,2)*O(lambda1,1,lambda3,3)-OP(2,lambda4,4)*pp(1,4)*O(lambda1,1,lambda3,3)+OP(1,lambda4,4)*(OP(2,lambda3,3)*OP(3,lambda1,1)-OP(2,lambda1,1)*OP(4,lambda3,3)+(-pp(2,3)+pp(2,4))*O(lambda1,1,lambda3,3))+OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda4,4)-OP(2,lambda3,3)*pp(1,3)*O(lambda1,1,lambda4,4)+OP(1,lambda3,3)*(-(OP(2,lambda1,1)*OP(3,lambda4,4))+OP(2,lambda4,4)*OP(4,lambda1,1)+(pp(2,3)-pp(2,4))*O(lambda1,1,lambda4,4))-OP(3,lambda1,1)*pp(1,2)*O(lambda3,3,lambda4,4)-OP(4,lambda1,1)*pp(1,2)*O(lambda3,3,lambda4,4)+OP(2,lambda1,1)*pp(1,3)*O(lambda3,3,lambda4,4)+OP(2,lambda1,1)*pp(1,4)*O(lambda3,3,lambda4,4)))/(8*WMASS)

               enddo                      
              enddo
!             enddo 
            enddo
         endif

         if( fm5hc.ne.0d0 ) then  !fm5hc coupling
            do lambda1 = -1,1
!             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,0,lambda3,lambda4) = M(lambda1,0,lambda3,lambda4)  +(fm5hc*(gw**2)*sw*(ZMASS**2)*1.0*(OP(3,lambda4,4)*pp(1,2)*O(lambda1,1,lambda3,3)-OP(2,lambda4,4)*pp(1,4)*O(lambda1,1,lambda3,3)+OP(1,lambda4,4)*(OP(2,lambda3,3)*OP(3,lambda1,1)-OP(2,lambda1,1)*OP(4,lambda3,3)+(-pp(2,3)+pp(2,4))*O(lambda1,1,lambda3,3))+OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda4,4)-OP(2,lambda3,3)*pp(1,3)*O(lambda1,1,lambda4,4)+OP(1,lambda3,3)*(-(OP(2,lambda1,1)*OP(3,lambda4,4))+OP(2,lambda4,4)*OP(4,lambda1,1)+(pp(2,3)-pp(2,4))*O(lambda1,1,lambda4,4))-OP(3,lambda1,1)*pp(1,2)*O(lambda3,3,lambda4,4)-OP(4,lambda1,1)*pp(1,2)*O(lambda3,3,lambda4,4)+OP(2,lambda1,1)*pp(1,3)*O(lambda3,3,lambda4,4)+OP(2,lambda1,1)*pp(1,4)*O(lambda3,3,lambda4,4)))/(8*WMASS)

               enddo                      
              enddo
!             enddo 
            enddo
         endif



      end





ccccccccccccccc
c ahwmwp-vertex c
ccccccccccccccc

      subroutine ahwmwp_ONanom4(a1,h,wm1,wp1,M)
c
c This subroutine calculates the vvvv-vertex with anomalous couplings
c and the on-shell polarization vectors for vector bosons
c
c IN:
c     polv1, polv2      : the incoming polvector wz(bosons in the CM frame)
c     polv3, polv4      : the outgoing polvector wz(bosons in the CM frame)
c OUT:
c     M(lambda1,lambda2) 
c
      IMPLICIT NONE

      INTEGER NGRAPHS
      PARAMETER( NGRAPHS=20 )
      DOUBLE PRECISION rZERO, rONE, rHalf
      PARAMETER( rZERO=0.0D0, rONE=1.0D0, rHALF = 0.5D0 )
      DOUBLE COMPLEX cIMAG, cZERO, cONE
      PARAMETER( cIMAG=(0.0D0,1.0D0), cZERO=(0.0D0,0.0D0), cONE=(1.0D0,0.0D0))
      DOUBLE COMPLEX a1(6,-1:1), h(6,-1:1), wm1(6,-1:1), wp1(6,-1:1)      
      DOUBLE COMPLEX M(-1:1,-1:1,-1:1,-1:1)
      
      DOUBLE COMPLEX OFFeps(6,-1:1,NGRAPHS), OFFepsC(6,-1:1,NGRAPHS)
      DOUBLE PRECISION K(0:3,4)      
      DOUBLE COMPLEX O(-1:1,4,-1:1,4), OP(4,-1:1,4)
      DOUBLE COMPLEXPP(1:4,1:4)
    
      INTEGER I,J,N,MU,NU
      INTEGER LAMBDA1, LAMBDA2, LAMBDA3, LAMBDA4, LAMBDA
      INTEGER L1,L2

      DOUBLE COMPLEX DOTCC, DOTRC,contract_Tjj
      DOUBLE PRECISION DOTRR
      EXTERNAL DOTCC, DOTRC, DOTRR, contract_Tjj
      double precision sw2, sw
#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "gmunu.inc"
#include "transm.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"
      sw2 = sin2w
      sw = sqrt(sin2w)
      OFFeps(6,-1:1,4) = cZero
      K(0:3,1:4) = rZero
      M(-1:1,-1:1,-1:1,-1:1) = cZero
      PP(1:4,1:4) = rZero
      O(-1:1,1:4,-1:1,1:4) = cZero
      OP(1:4,-1:1,1:4) = cZero  
      
      OFFeps(1:6,-1:1,1) = a1(1:6,-1:1)
      OFFeps(1:6,-1:1,2) = h(1:6,-1:1)
      OFFeps(1:6,-1:1,3) = wm1(1:6,-1:1)
      OFFeps(1:6,-1:1,4) = wp1(1:6,-1:1) 

      do i = 1,4
        K(0,i) = dble(OFFeps(5,-1,i))*(-rONE) !K1(0:3)
        K(1,i) = dble(OFFeps(6,-1,i))*(-rONE) !K2(0:3)
        K(2,i) = dimag(OFFeps(6,-1,i))*(-rONE)!K3(0:3)
        K(3,i) = dimag(OFFeps(5,-1,i))*(-rONE)!K4(0:3)
      enddo
      



      do i = 1,4
        do j = 1,4
          pp(i,j) = dotrr(k(0,i),k(0,j))
        enddo
      enddo

      do L1 = -1,1
        do L2 = -1,1 
          do i = 1,4
            do j = 1,4
              O(L1,i,L2,j) = dotcc(OFFeps(1,L1,i),OFFeps(1,L2,j))
            enddo
          enddo
        enddo
      enddo

      do LAMBDA = -1,1
        do i = 1,4 
          do j = 1,4
            OP(j,lambda,i) = dotrc(k(0,j),OFFeps(1,lambda,i))
          enddo
        enddo
      enddo


!!!! ===============================================================      
!!!! COUPLINGS CALCULATION
!!!! ===============================================================




         if( fm4.ne.0d0 ) then  !fm4 coupling
            do lambda1 = -1,1
!             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,0,lambda3,lambda4) = M(lambda1,0,lambda3,lambda4)  -(fm4*(gw**2)*sw*WMASS*1.0*(OP(1,lambda3,3)*OP(2,lambda4,4)*OP(3,lambda1,1)-OP(1,lambda4,4)*OP(2,lambda3,3)*OP(4,lambda1,1)-OP(2,lambda4,4)*pp(1,3)*O(lambda1,1,lambda3,3)+OP(2,lambda3,3)*pp(1,4)*O(lambda1,1,lambda4,4)))/2

               enddo                      
              enddo
!             enddo 
            enddo
         endif

         if( fm5.ne.0d0 ) then  !fm5 coupling
            do lambda1 = -1,1
!             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,0,lambda3,lambda4) = M(lambda1,0,lambda3,lambda4)  -(fm5*(gw**2)*sw*WMASS*1.0*(OP(1,lambda4,4)*(OP(2,lambda3,3)*OP(3,lambda1,1)-OP(2,lambda1,1)*OP(4,lambda3,3)-pp(2,3)*O(lambda1,1,lambda3,3))+OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda4,4)-OP(2,lambda3,3)*pp(1,3)*O(lambda1,1,lambda4,4)+OP(1,lambda3,3)*pp(2,3)*O(lambda1,1,lambda4,4)-OP(4,lambda1,1)*pp(1,2)*O(lambda3,3,lambda4,4)+OP(2,lambda1,1)*pp(1,4)*O(lambda3,3,lambda4,4)))/4

               enddo                      
              enddo
!             enddo 
            enddo
         endif

         if( fm5hc.ne.0d0 ) then  !fm5hc coupling
            do lambda1 = -1,1
!             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,0,lambda3,lambda4) = M(lambda1,0,lambda3,lambda4)  +(fm5hc*(gw**2)*sw*WMASS*1.0*(OP(3,lambda4,4)*pp(1,2)*O(lambda1,1,lambda3,3)-OP(2,lambda4,4)*pp(1,4)*O(lambda1,1,lambda3,3)+OP(1,lambda4,4)*pp(2,4)*O(lambda1,1,lambda3,3)-OP(1,lambda3,3)*(OP(2,lambda1,1)*OP(3,lambda4,4)-OP(2,lambda4,4)*OP(4,lambda1,1)+pp(2,4)*O(lambda1,1,lambda4,4))-OP(3,lambda1,1)*pp(1,2)*O(lambda3,3,lambda4,4)+OP(2,lambda1,1)*pp(1,3)*O(lambda3,3,lambda4,4)))/4

               enddo                      
              enddo
!             enddo 
            enddo
         endif



         if( fm7.ne.0d0 ) then  !fm7 coupling
            do lambda1 = -1,1
!             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,0,lambda3,lambda4) = M(lambda1,0,lambda3,lambda4)  +(fm7*(gw**2)*sw*WMASS*1.0*(OP(3,lambda4,4)*pp(1,2)*O(lambda1,1,lambda3,3)+OP(2,lambda4,4)*pp(1,4)*O(lambda1,1,lambda3,3)+OP(1,lambda4,4)*(OP(2,lambda3,3)*OP(3,lambda1,1)+OP(2,lambda1,1)*OP(4,lambda3,3)-(pp(2,3)+pp(2,4))*O(lambda1,1,lambda3,3))-OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda4,4)-OP(2,lambda3,3)*pp(1,3)*O(lambda1,1,lambda4,4)+OP(1,lambda3,3)*(-(OP(2,lambda1,1)*OP(3,lambda4,4))-OP(2,lambda4,4)*OP(4,lambda1,1)+(pp(2,3)+pp(2,4))*O(lambda1,1,lambda4,4))-OP(3,lambda1,1)*pp(1,2)*O(lambda3,3,lambda4,4)+OP(4,lambda1,1)*pp(1,2)*O(lambda3,3,lambda4,4)+OP(2,lambda1,1)*pp(1,3)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*pp(1,4)*O(lambda3,3,lambda4,4)))/4

               enddo                      
              enddo
!             enddo 
            enddo
         endif         

      end



ccccccccccccccc
c azzz-vertex c
ccccccccccccccc

      subroutine azzz_ONanom4(a1,z1,z2,z3,M)
c
c This subroutine calculates the vvvv-vertex with anomalous couplings
c and the on-shell polarization vectors for vector bosons
c
c IN:
c     polv1, polv2      : the incoming polvector wz(bosons in the CM frame)
c     polv3, polv4      : the outgoing polvector wz(bosons in the CM frame)
c OUT:
c     M(lambda1,lambda2) 
c
      IMPLICIT NONE

      INTEGER NGRAPHS
      PARAMETER( NGRAPHS=20 )
      DOUBLE PRECISION rZERO, rONE, rHalf
      PARAMETER( rZERO=0.0D0, rONE=1.0D0, rHALF = 0.5D0 )
      DOUBLE COMPLEX cIMAG, cZERO, cONE
      PARAMETER( cIMAG=(0.0D0,1.0D0), cZERO=(0.0D0,0.0D0), cONE=(1.0D0,0.0D0))
      DOUBLE COMPLEX a1(6,-1:1), z1(6,-1:1), z2(6,-1:1), z3(6,-1:1)      
      DOUBLE COMPLEX M(-1:1,-1:1,-1:1,-1:1)
      
      DOUBLE COMPLEX OFFeps(6,-1:1,NGRAPHS), OFFepsC(6,-1:1,NGRAPHS)
      DOUBLE PRECISION K(0:3,4)      
      DOUBLE COMPLEX O(-1:1,4,-1:1,4), OP(4,-1:1,4)
      DOUBLE COMPLEXPP(1:4,1:4)
    
      INTEGER I,J,N,MU,NU
      INTEGER LAMBDA1, LAMBDA2, LAMBDA3, LAMBDA4, LAMBDA
      INTEGER L1,L2

      DOUBLE COMPLEX DOTCC, DOTRC,contract_Tjj
      DOUBLE PRECISION DOTRR
      EXTERNAL DOTCC, DOTRC, DOTRR, contract_Tjj
      double precision sw2, sw
#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "gmunu.inc"
#include "transm.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"
      sw2 = sin2w
      sw = sqrt(sin2w)
      OFFeps(6,-1:1,4) = cZero
      K(0:3,1:4) = rZero
      M(-1:1,-1:1,-1:1,-1:1) = cZero
      PP(1:4,1:4) = rZero
      O(-1:1,1:4,-1:1,1:4) = cZero
      OP(1:4,-1:1,1:4) = cZero  
      
      OFFeps(1:6,-1:1,1) = a1(1:6,-1:1)
      OFFeps(1:6,-1:1,2) = z1(1:6,-1:1)
      OFFeps(1:6,-1:1,3) = z2(1:6,-1:1)
      OFFeps(1:6,-1:1,4) = z3(1:6,-1:1) 

      do i = 1,4
        K(0,i) = dble(OFFeps(5,-1,i))*(-rONE) !K1(0:3)
        K(1,i) = dble(OFFeps(6,-1,i))*(-rONE) !K2(0:3)
        K(2,i) = dimag(OFFeps(6,-1,i))*(-rONE)!K3(0:3)
        K(3,i) = dimag(OFFeps(5,-1,i))*(-rONE)!K4(0:3)
      enddo
      



      do i = 1,4
        do j = 1,4
          pp(i,j) = dotrr(k(0,i),k(0,j))
        enddo
      enddo

      do L1 = -1,1
        do L2 = -1,1 
          do i = 1,4
            do j = 1,4
              O(L1,i,L2,j) = dotcc(OFFeps(1,L1,i),OFFeps(1,L2,j))
            enddo
          enddo
        enddo
      enddo

      do LAMBDA = -1,1
        do i = 1,4 
          do j = 1,4
            OP(j,lambda,i) = dotrc(k(0,j),OFFeps(1,lambda,i))
          enddo
        enddo
      enddo


!!!! ===============================================================      
!!!! COUPLINGS CALCULATION
!!!! ===============================================================


                

         if( ft0.ne.0d0 ) then  !ft0 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(-8*ft0*(gw**4)*sw*(-1+SW2)*WMASS*(OP(1,lambda2,2)*OP(2,lambda1,1)*OP(3,lambda4,4)*OP(4,lambda3,3)-OP(3,lambda4,4)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda2,2)-OP(2,lambda4,4)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda3,3)-OP(2,lambda3,3)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*OP(4,lambda1,1)*(OP(2,lambda3,3)*OP(3,lambda2,2)-pp(2,3)*O(lambda2,2,lambda3,3))+pp(1,3)*pp(2,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*OP(3,lambda1,1)*(OP(2,lambda4,4)*OP(4,lambda2,2)-pp(2,4)*O(lambda2,2,lambda4,4))-OP(1,lambda2,2)*OP(2,lambda1,1)*pp(3,4)*O(lambda3,3,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)))/ZMASS

               enddo                      
              enddo
             enddo 
            enddo
         endif


         if( ft1.ne.0d0 ) then  !ft1 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(-8*ft1*(gw**4)*sw*(-1+SW2)*WMASS*(OP(1,lambda2,2)*OP(2,lambda1,1)*OP(3,lambda4,4)*OP(4,lambda3,3)-OP(3,lambda4,4)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda2,2)-OP(2,lambda4,4)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda3,3)-OP(2,lambda3,3)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*OP(4,lambda1,1)*(OP(2,lambda3,3)*OP(3,lambda2,2)-pp(2,3)*O(lambda2,2,lambda3,3))+pp(1,3)*pp(2,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*OP(3,lambda1,1)*(OP(2,lambda4,4)*OP(4,lambda2,2)-pp(2,4)*O(lambda2,2,lambda4,4))-OP(1,lambda2,2)*OP(2,lambda1,1)*pp(3,4)*O(lambda3,3,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)))/ZMASS

               enddo                      
              enddo
             enddo 
            enddo
         endif


         if( ft2.ne.0d0 ) then  !ft2 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(-2*ft2*(gw**4)*sw*(-1+SW2)*WMASS*(OP(1,lambda4,4)*OP(2,lambda3,3)*OP(3,lambda1,1)*OP(4,lambda2,2)+OP(1,lambda4,4)*OP(2,lambda1,1)*OP(3,lambda2,2)*OP(4,lambda3,3)-OP(2,lambda4,4)*OP(4,lambda3,3)*pp(1,3)*O(lambda1,1,lambda2,2)-OP(2,lambda3,3)*OP(3,lambda4,4)*pp(1,4)*O(lambda1,1,lambda2,2)-OP(1,lambda4,4)*OP(4,lambda3,3)*pp(2,3)*O(lambda1,1,lambda2,2)+OP(1,lambda4,4)*OP(2,lambda3,3)*pp(3,4)*O(lambda1,1,lambda2,2)-OP(3,lambda4,4)*OP(4,lambda2,2)*pp(1,2)*O(lambda1,1,lambda3,3)-OP(2,lambda4,4)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda3,3)-OP(1,lambda4,4)*OP(4,lambda2,2)*pp(2,3)*O(lambda1,1,lambda3,3)+OP(1,lambda4,4)*OP(3,lambda2,2)*pp(2,4)*O(lambda1,1,lambda3,3)-OP(3,lambda2,2)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda4,4)-OP(2,lambda3,3)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda4,4)-OP(3,lambda4,4)*OP(4,lambda1,1)*pp(1,2)*O(lambda2,2,lambda3,3)-OP(2,lambda4,4)*OP(4,lambda1,1)*pp(1,3)*O(lambda2,2,lambda3,3)+OP(2,lambda4,4)*OP(3,lambda1,1)*pp(1,4)*O(lambda2,2,lambda3,3)+OP(2,lambda1,1)*OP(3,lambda4,4)*pp(1,4)*O(lambda2,2,lambda3,3)-OP(1,lambda4,4)*OP(3,lambda1,1)*pp(2,4)*O(lambda2,2,lambda3,3)-OP(1,lambda4,4)*OP(2,lambda1,1)*pp(3,4)*O(lambda2,2,lambda3,3)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-OP(3,lambda1,1)*OP(4,lambda3,3)*pp(1,2)*O(lambda2,2,lambda4,4)+OP(2,lambda3,3)*OP(4,lambda1,1)*pp(1,3)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*OP(4,lambda3,3)*pp(1,3)*O(lambda2,2,lambda4,4)-OP(2,lambda3,3)*OP(3,lambda1,1)*pp(1,4)*O(lambda2,2,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(-(OP(3,lambda4,4)*pp(2,4)*O(lambda1,1,lambda2,2))+OP(2,lambda4,4)*(OP(3,lambda2,2)*OP(4,lambda1,1)+pp(3,4)*O(lambda1,1,lambda2,2))+OP(4,lambda2,2)*pp(2,3)*O(lambda1,1,lambda4,4)-OP(3,lambda2,2)*pp(2,4)*O(lambda1,1,lambda4,4)-OP(4,lambda1,1)*pp(2,3)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*(OP(3,lambda4,4)*OP(4,lambda2,2)-pp(3,4)*O(lambda2,2,lambda4,4)))+OP(3,lambda2,2)*OP(4,lambda1,1)*pp(1,2)*O(lambda3,3,lambda4,4)+OP(3,lambda1,1)*OP(4,lambda2,2)*pp(1,2)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(4,lambda2,2)*pp(1,3)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(3,lambda2,2)*pp(1,4)*O(lambda3,3,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+OP(1,lambda2,2)*(OP(3,lambda4,4)*pp(2,4)*O(lambda1,1,lambda3,3)+OP(2,lambda4,4)*(OP(3,lambda1,1)*OP(4,lambda3,3)-pp(3,4)*O(lambda1,1,lambda3,3))+OP(4,lambda3,3)*pp(2,3)*O(lambda1,1,lambda4,4)+OP(2,lambda3,3)*(OP(3,lambda4,4)*OP(4,lambda1,1)-pp(3,4)*O(lambda1,1,lambda4,4))-OP(4,lambda1,1)*pp(2,3)*O(lambda3,3,lambda4,4)-OP(3,lambda1,1)*pp(2,4)*O(lambda3,3,lambda4,4))))/ZMASS

               enddo                      
              enddo
             enddo 
            enddo
         endif
         
         if( ft3.ne.0d0 ) then  !ft3 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) -(1/ZMASS)*2*gw**4*ft3*sw*(SW2-1)*WMASS*(OP(1,lambda4,4)*OP(2,lambda3,3)*OP(3,lambda1,1)*OP(4,lambda2,2)-OP(3,lambda4,4)*pp(1,2)*O(lambda1,1,lambda3,3)*OP(4,lambda2,2)-OP(1,lambda4,4)*pp(2,3)*O(lambda1,1,lambda3,3)*OP(4,lambda2,2)-OP(2,lambda3,3)*pp(1,3)*O(lambda1,1,lambda4,4)*OP(4,lambda2,2)+OP(3,lambda1,1)*pp(1,2)*O(lambda3,3,lambda4,4)*OP(4,lambda2,2)-OP(2,lambda1,1)*pp(1,3)*O(lambda3,3,lambda4,4)*OP(4,lambda2,2)+OP(1,lambda4,4)*OP(2,lambda1,1)*OP(3,lambda2,2)*OP(4,lambda3,3)-OP(2,lambda4,4)*OP(4,lambda3,3)*pp(1,3)*O(lambda1,1,lambda2,2)-OP(2,lambda3,3)*OP(3,lambda4,4)*pp(1,4)*O(lambda1,1,lambda2,2)-OP(1,lambda4,4)*OP(4,lambda3,3)*pp(2,3)*O(lambda1,1,lambda2,2)+OP(1,lambda4,4)*OP(2,lambda3,3)*pp(3,4)*O(lambda1,1,lambda2,2)-OP(2,lambda4,4)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda3,3)+OP(1,lambda4,4)*OP(3,lambda2,2)*pp(2,4)*O(lambda1,1,lambda3,3)-OP(3,lambda2,2)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda4,4)-OP(3,lambda4,4)*OP(4,lambda1,1)*pp(1,2)*O(lambda2,2,lambda3,3)-OP(2,lambda4,4)*OP(4,lambda1,1)*pp(1,3)*O(lambda2,2,lambda3,3)+OP(2,lambda4,4)*OP(3,lambda1,1)*pp(1,4)*O(lambda2,2,lambda3,3)+OP(2,lambda1,1)*OP(3,lambda4,4)*pp(1,4)*O(lambda2,2,lambda3,3)-OP(1,lambda4,4)*OP(3,lambda1,1)*pp(2,4)*O(lambda2,2,lambda3,3)-OP(1,lambda4,4)*OP(2,lambda1,1)*pp(3,4)*O(lambda2,2,lambda3,3)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-OP(3,lambda1,1)*OP(4,lambda3,3)*pp(1,2)*O(lambda2,2,lambda4,4)+OP(2,lambda3,3)*OP(4,lambda1,1)*pp(1,3)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*OP(4,lambda3,3)*pp(1,3)*O(lambda2,2,lambda4,4)-OP(2,lambda3,3)*OP(3,lambda1,1)*pp(1,4)*O(lambda2,2,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(-OP(3,lambda4,4)*pp(2,4)*O(lambda1,1,lambda2,2)+OP(2,lambda4,4)*(OP(3,lambda2,2)*OP(4,lambda1,1)+pp(3,4)*O(lambda1,1,lambda2,2))+OP(4,lambda2,2)*pp(2,3)*O(lambda1,1,lambda4,4)-OP(3,lambda2,2)*pp(2,4)*O(lambda1,1,lambda4,4)-OP(4,lambda1,1)*pp(2,3)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*(OP(3,lambda4,4)*OP(4,lambda2,2)-pp(3,4)*O(lambda2,2,lambda4,4)))+OP(3,lambda2,2)*OP(4,lambda1,1)*pp(1,2)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(3,lambda2,2)*pp(1,4)*O(lambda3,3,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+OP(1,lambda2,2)*(OP(3,lambda4,4)*pp(2,4)*O(lambda1,1,lambda3,3)+OP(2,lambda4,4)*(OP(3,lambda1,1)*OP(4,lambda3,3)-pp(3,4)*O(lambda1,1,lambda3,3))+OP(4,lambda3,3)*pp(2,3)*O(lambda1,1,lambda4,4)+OP(2,lambda3,3)*(OP(3,lambda4,4)*OP(4,lambda1,1)-pp(3,4)*O(lambda1,1,lambda4,4))-OP(4,lambda1,1)*pp(2,3)*O(lambda3,3,lambda4,4)-OP(3,lambda1,1)*pp(2,4)*O(lambda3,3,lambda4,4)))

               enddo                      
              enddo
             enddo 
            enddo
         endif


         if( ft4.ne.0d0 ) then  !ft4 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4)+(gw**4*ft4*sw**3*(-1 + 2*SW2)*ZMASS*(OP(1, lambda4, 4)*OP(2, lambda3, 3)*OP(3, lambda1, 1)*OP(4, lambda2, 2) + OP(1, lambda4, 4)*OP(2, lambda1, 1)*OP(3, lambda2, 2)*OP(4, lambda3, 3) + pp(3, 4)*OP(1, lambda4, 4)*OP(2, lambda3, 3)*O(lambda1, 1, lambda2, 2) - pp(1, 4)*OP(2, lambda3, 3)*OP(3, lambda4, 4)*O(lambda1, 1, lambda2, 2) - pp(2, 3)*OP(1, lambda4, 4)*OP(4, lambda3, 3)*O(lambda1, 1, lambda2, 2) - pp(1, 3)*OP(2, lambda4, 4)*OP(4, lambda3, 3)*O(lambda1, 1, lambda2, 2) + pp(2, 4)*OP(1, lambda4, 4)*OP(3, lambda2, 2)*O(lambda1, 1, lambda3, 3) - pp(1, 4)*OP(2, lambda4, 4)*OP(3, lambda2, 2)*O(lambda1, 1, lambda3, 3) - pp(2, 3)*OP(1, lambda4, 4)*OP(4, lambda2, 2)*O(lambda1, 1, lambda3, 3) - pp(1, 2)*OP(3, lambda4, 4)*OP(4, lambda2, 2)*O(lambda1, 1, lambda3, 3) - pp(1, 3)*OP(2, lambda3, 3)*OP(4, lambda2, 2)*O(lambda1, 1, lambda4, 4) - pp(1, 2)*OP(3, lambda2, 2)*OP(4, lambda3, 3)*O(lambda1, 1, lambda4, 4) - pp(3, 4)*OP(1, lambda4, 4)*OP(2, lambda1, 1)*O(lambda2, 2, lambda3, 3) - pp(2, 4)*OP(1, lambda4, 4)*OP(3, lambda1, 1)*O(lambda2, 2, lambda3, 3) + pp(1, 4)*OP(2, lambda4, 4)*OP(3, lambda1, 1)*O(lambda2, 2, lambda3, 3) + pp(1, 4)*OP(2, lambda1, 1)*OP(3, lambda4, 4)*O(lambda2, 2, lambda3, 3) - pp(1, 3)*OP(2, lambda4, 4)*OP(4, lambda1, 1)*O(lambda2, 2, lambda3, 3) - pp(1, 2)*OP(3, lambda4, 4)*OP(4, lambda1, 1)*O(lambda2, 2, lambda3, 3) + pp(1, 3)*pp(2, 4)*O(lambda1, 1, lambda4, 4)*O(lambda2, 2, lambda3, 3) + pp(1, 2)*pp(3, 4)*O(lambda1, 1, lambda4, 4)*O(lambda2, 2, lambda3, 3) - pp(1, 4)*OP(2, lambda3, 3)*OP(3, lambda1, 1)*O(lambda2, 2, lambda4, 4) + pp(1, 3)*OP(2, lambda3, 3)*OP(4, lambda1, 1)*O(lambda2, 2, lambda4, 4) + pp(1, 3)*OP(2, lambda1, 1)*OP(4, lambda3, 3)*O(lambda2, 2, lambda4, 4) - pp(1, 2)*OP(3, lambda1, 1)*OP(4, lambda3, 3)*O(lambda2, 2, lambda4, 4) + pp(1, 4)*pp(2, 3)*O(lambda1, 1, lambda3, 3)*O(lambda2, 2, lambda4, 4) + pp(1, 2)*pp(3, 4)*O(lambda1, 1, lambda3, 3)*O(lambda2, 2, lambda4, 4) + OP(1, lambda3, 3)*(-(pp(2, 4)*OP(3, lambda4, 4)*O(lambda1, 1, lambda2, 2)) + OP(2, lambda4, 4)*(OP(3, lambda2, 2)*OP(4, lambda1, 1) + pp(3, 4)*O(lambda1, 1, lambda2, 2)) -   pp(2, 4)*OP(3, lambda2, 2)*O(lambda1, 1, lambda4, 4) + pp(2, 3)*OP(4, lambda2, 2)*O(lambda1, 1, lambda4, 4) - pp(2, 3)*OP(4, lambda1, 1)*O(lambda2, 2, lambda4, 4) +   OP(2, lambda1, 1)*(OP(3, lambda4, 4)*OP(4, lambda2, 2) - pp(3, 4)*O(lambda2, 2, lambda4, 4))) - pp(1, 4)*OP(2, lambda1, 1)*OP(3, lambda2, 2)*O(lambda3, 3, lambda4, 4) + pp(1, 2)*OP(3, lambda2, 2)*OP(4, lambda1, 1)*O(lambda3, 3, lambda4, 4) - pp(1, 3)*OP(2, lambda1, 1)*OP(4, lambda2, 2)*O(lambda3, 3, lambda4, 4) + pp(1, 2)*OP(3, lambda1, 1)*OP(4, lambda2, 2)*O(lambda3, 3, lambda4, 4) + pp(1, 4)*pp(2, 3)*O(lambda1, 1, lambda2, 2)*O(lambda3, 3, lambda4, 4) + pp(1, 3)*pp(2, 4)*O(lambda1, 1, lambda2, 2)*O(lambda3, 3, lambda4, 4) + OP(1, lambda2, 2)*(pp(2, 4)*OP(3, lambda4, 4)*O(lambda1, 1, lambda3, 3) +   OP(2, lambda4, 4)*(OP(3, lambda1, 1)*OP(4, lambda3, 3) - pp(3, 4)*O(lambda1, 1, lambda3, 3)) + pp(2, 3)*OP(4, lambda3, 3)*O(lambda1, 1, lambda4, 4) +   OP(2, lambda3, 3)*(OP(3, lambda4, 4)*OP(4, lambda1, 1) - pp(3, 4)*O(lambda1, 1, lambda4, 4)) - pp(2, 4)*OP(3, lambda1, 1)*O(lambda3, 3, lambda4, 4) -   pp(2, 3)*OP(4, lambda1, 1)*O(lambda3, 3, lambda4, 4))))/(2*WMASS)
 
               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( ft5.ne.0d0 ) then  !ft5 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(2*ft5*(gw**4)*(sw**3)*(-1+2*SW2)*ZMASS*(OP(1,lambda2,2)*OP(2,lambda1,1)*OP(3,lambda4,4)*OP(4,lambda3,3)-OP(3,lambda4,4)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda2,2)-OP(2,lambda4,4)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda3,3)-OP(2,lambda3,3)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*OP(4,lambda1,1)*(OP(2,lambda3,3)*OP(3,lambda2,2)-pp(2,3)*O(lambda2,2,lambda3,3))+pp(1,3)*pp(2,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*OP(3,lambda1,1)*(OP(2,lambda4,4)*OP(4,lambda2,2)-pp(2,4)*O(lambda2,2,lambda4,4))-OP(1,lambda2,2)*OP(2,lambda1,1)*pp(3,4)*O(lambda3,3,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)))/WMASS

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( ft6.ne.0d0 ) then  !ft6 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(2*ft6*(gw**4)*(sw**3)*(-1+2*SW2)*ZMASS*(OP(1,lambda2,2)*OP(2,lambda1,1)*OP(3,lambda4,4)*OP(4,lambda3,3)-OP(3,lambda4,4)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda2,2)-OP(2,lambda4,4)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda3,3)-OP(2,lambda3,3)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*OP(4,lambda1,1)*(OP(2,lambda3,3)*OP(3,lambda2,2)-pp(2,3)*O(lambda2,2,lambda3,3))+pp(1,3)*pp(2,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*OP(3,lambda1,1)*(OP(2,lambda4,4)*OP(4,lambda2,2)-pp(2,4)*O(lambda2,2,lambda4,4))-OP(1,lambda2,2)*OP(2,lambda1,1)*pp(3,4)*O(lambda3,3,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)))/WMASS

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( ft7.ne.0d0 ) then  !ft7 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(ft7*(gw**4)*(sw**3)*(-1+2*SW2)*ZMASS*(OP(1,lambda4,4)*OP(2,lambda3,3)*OP(3,lambda1,1)*OP(4,lambda2,2)+OP(1,lambda4,4)*OP(2,lambda1,1)*OP(3,lambda2,2)*OP(4,lambda3,3)-OP(2,lambda4,4)*OP(4,lambda3,3)*pp(1,3)*O(lambda1,1,lambda2,2)-OP(2,lambda3,3)*OP(3,lambda4,4)*pp(1,4)*O(lambda1,1,lambda2,2)-OP(1,lambda4,4)*OP(4,lambda3,3)*pp(2,3)*O(lambda1,1,lambda2,2)+OP(1,lambda4,4)*OP(2,lambda3,3)*pp(3,4)*O(lambda1,1,lambda2,2)-OP(3,lambda4,4)*OP(4,lambda2,2)*pp(1,2)*O(lambda1,1,lambda3,3)-OP(2,lambda4,4)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda3,3)-OP(1,lambda4,4)*OP(4,lambda2,2)*pp(2,3)*O(lambda1,1,lambda3,3)+OP(1,lambda4,4)*OP(3,lambda2,2)*pp(2,4)*O(lambda1,1,lambda3,3)-OP(3,lambda2,2)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda4,4)-OP(2,lambda3,3)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda4,4)-OP(3,lambda4,4)*OP(4,lambda1,1)*pp(1,2)*O(lambda2,2,lambda3,3)-OP(2,lambda4,4)*OP(4,lambda1,1)*pp(1,3)*O(lambda2,2,lambda3,3)+OP(2,lambda4,4)*OP(3,lambda1,1)*pp(1,4)*O(lambda2,2,lambda3,3)+OP(2,lambda1,1)*OP(3,lambda4,4)*pp(1,4)*O(lambda2,2,lambda3,3)-OP(1,lambda4,4)*OP(3,lambda1,1)*pp(2,4)*O(lambda2,2,lambda3,3)-OP(1,lambda4,4)*OP(2,lambda1,1)*pp(3,4)*O(lambda2,2,lambda3,3)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-OP(3,lambda1,1)*OP(4,lambda3,3)*pp(1,2)*O(lambda2,2,lambda4,4)+OP(2,lambda3,3)*OP(4,lambda1,1)*pp(1,3)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*OP(4,lambda3,3)*pp(1,3)*O(lambda2,2,lambda4,4)-OP(2,lambda3,3)*OP(3,lambda1,1)*pp(1,4)*O(lambda2,2,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(-(OP(3,lambda4,4)*pp(2,4)*O(lambda1,1,lambda2,2))+OP(2,lambda4,4)*(OP(3,lambda2,2)*OP(4,lambda1,1)+pp(3,4)*O(lambda1,1,lambda2,2))+OP(4,lambda2,2)*pp(2,3)*O(lambda1,1,lambda4,4)-OP(3,lambda2,2)*pp(2,4)*O(lambda1,1,lambda4,4)-OP(4,lambda1,1)*pp(2,3)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*(OP(3,lambda4,4)*OP(4,lambda2,2)-pp(3,4)*O(lambda2,2,lambda4,4)))+OP(3,lambda2,2)*OP(4,lambda1,1)*pp(1,2)*O(lambda3,3,lambda4,4)+OP(3,lambda1,1)*OP(4,lambda2,2)*pp(1,2)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(4,lambda2,2)*pp(1,3)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(3,lambda2,2)*pp(1,4)*O(lambda3,3,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+OP(1,lambda2,2)*(OP(3,lambda4,4)*pp(2,4)*O(lambda1,1,lambda3,3)+OP(2,lambda4,4)*(OP(3,lambda1,1)*OP(4,lambda3,3)-pp(3,4)*O(lambda1,1,lambda3,3))+OP(4,lambda3,3)*pp(2,3)*O(lambda1,1,lambda4,4)+OP(2,lambda3,3)*(OP(3,lambda4,4)*OP(4,lambda1,1)-pp(3,4)*O(lambda1,1,lambda4,4))-OP(4,lambda1,1)*pp(2,3)*O(lambda3,3,lambda4,4)-OP(3,lambda1,1)*pp(2,4)*O(lambda3,3,lambda4,4))))/(2*WMASS)

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( ft8.ne.0d0 ) then  !ft8 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(-2*ft8*(gw**4)*(sw**7)*(ZMASS**3)*(OP(1,lambda2,2)*OP(2,lambda1,1)*OP(3,lambda4,4)*OP(4,lambda3,3)-OP(3,lambda4,4)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda2,2)-OP(2,lambda4,4)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda3,3)-OP(2,lambda3,3)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*OP(4,lambda1,1)*(OP(2,lambda3,3)*OP(3,lambda2,2)-pp(2,3)*O(lambda2,2,lambda3,3))+pp(1,3)*pp(2,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*OP(3,lambda1,1)*(OP(2,lambda4,4)*OP(4,lambda2,2)-pp(2,4)*O(lambda2,2,lambda4,4))-OP(1,lambda2,2)*OP(2,lambda1,1)*pp(3,4)*O(lambda3,3,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)))/(WMASS**3)

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( ft9.ne.0d0 ) then  !ft9 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) -(ft9*(gw**4)*(sw**7)*(ZMASS**3)*(OP(1,lambda4,4)*OP(2,lambda3,3)*OP(3,lambda1,1)*OP(4,lambda2,2)+OP(1,lambda4,4)*OP(2,lambda1,1)*OP(3,lambda2,2)*OP(4,lambda3,3)-OP(2,lambda4,4)*OP(4,lambda3,3)*pp(1,3)*O(lambda1,1,lambda2,2)-OP(2,lambda3,3)*OP(3,lambda4,4)*pp(1,4)*O(lambda1,1,lambda2,2)-OP(1,lambda4,4)*OP(4,lambda3,3)*pp(2,3)*O(lambda1,1,lambda2,2)+OP(1,lambda4,4)*OP(2,lambda3,3)*pp(3,4)*O(lambda1,1,lambda2,2)-OP(3,lambda4,4)*OP(4,lambda2,2)*pp(1,2)*O(lambda1,1,lambda3,3)-OP(2,lambda4,4)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda3,3)-OP(1,lambda4,4)*OP(4,lambda2,2)*pp(2,3)*O(lambda1,1,lambda3,3)+OP(1,lambda4,4)*OP(3,lambda2,2)*pp(2,4)*O(lambda1,1,lambda3,3)-OP(3,lambda2,2)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda4,4)-OP(2,lambda3,3)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda4,4)-OP(3,lambda4,4)*OP(4,lambda1,1)*pp(1,2)*O(lambda2,2,lambda3,3)-OP(2,lambda4,4)*OP(4,lambda1,1)*pp(1,3)*O(lambda2,2,lambda3,3)+OP(2,lambda4,4)*OP(3,lambda1,1)*pp(1,4)*O(lambda2,2,lambda3,3)+OP(2,lambda1,1)*OP(3,lambda4,4)*pp(1,4)*O(lambda2,2,lambda3,3)-OP(1,lambda4,4)*OP(3,lambda1,1)*pp(2,4)*O(lambda2,2,lambda3,3)-OP(1,lambda4,4)*OP(2,lambda1,1)*pp(3,4)*O(lambda2,2,lambda3,3)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-OP(3,lambda1,1)*OP(4,lambda3,3)*pp(1,2)*O(lambda2,2,lambda4,4)+OP(2,lambda3,3)*OP(4,lambda1,1)*pp(1,3)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*OP(4,lambda3,3)*pp(1,3)*O(lambda2,2,lambda4,4)-OP(2,lambda3,3)*OP(3,lambda1,1)*pp(1,4)*O(lambda2,2,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(-(OP(3,lambda4,4)*pp(2,4)*O(lambda1,1,lambda2,2))+OP(2,lambda4,4)*(OP(3,lambda2,2)*OP(4,lambda1,1)+pp(3,4)*O(lambda1,1,lambda2,2))+OP(4,lambda2,2)*pp(2,3)*O(lambda1,1,lambda4,4)-OP(3,lambda2,2)*pp(2,4)*O(lambda1,1,lambda4,4)-OP(4,lambda1,1)*pp(2,3)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*(OP(3,lambda4,4)*OP(4,lambda2,2)-pp(3,4)*O(lambda2,2,lambda4,4)))+OP(3,lambda2,2)*OP(4,lambda1,1)*pp(1,2)*O(lambda3,3,lambda4,4)+OP(3,lambda1,1)*OP(4,lambda2,2)*pp(1,2)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(4,lambda2,2)*pp(1,3)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(3,lambda2,2)*pp(1,4)*O(lambda3,3,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+OP(1,lambda2,2)*(OP(3,lambda4,4)*pp(2,4)*O(lambda1,1,lambda3,3)+OP(2,lambda4,4)*(OP(3,lambda1,1)*OP(4,lambda3,3)-pp(3,4)*O(lambda1,1,lambda3,3))+OP(4,lambda3,3)*pp(2,3)*O(lambda1,1,lambda4,4)+OP(2,lambda3,3)*(OP(3,lambda4,4)*OP(4,lambda1,1)-pp(3,4)*O(lambda1,1,lambda4,4))-OP(4,lambda1,1)*pp(2,3)*O(lambda3,3,lambda4,4)-OP(3,lambda1,1)*pp(2,4)*O(lambda3,3,lambda4,4))))/(2*(WMASS**3))

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fm0.ne.0d0 ) then  !fm0 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) -2*fm0*(gw**2)*sw*WMASS*ZMASS*(OP(1,lambda4,4)*OP(4,lambda1,1)*O(lambda2,2,lambda3,3)-pp(1,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda3,3)*OP(3,lambda1,1)*O(lambda2,2,lambda4,4)-pp(1,3)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda2,2)*OP(2,lambda1,1)*O(lambda3,3,lambda4,4)-pp(1,2)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4))

               enddo                      
              enddo
             enddo 
            enddo
         endif


         if( fm1.ne.0d0 ) then  !fm1 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(fm1*(gw**2)*sw*WMASS*ZMASS*(OP(1,lambda2,2)*OP(2,lambda4,4)*O(lambda1,1,lambda3,3)-OP(1,lambda2,2)*OP(3,lambda4,4)*O(lambda1,1,lambda3,3)+OP(1,lambda2,2)*OP(2,lambda3,3)*O(lambda1,1,lambda4,4)-OP(1,lambda2,2)*OP(4,lambda3,3)*O(lambda1,1,lambda4,4)-pp(1,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-pp(1,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*(-(OP(2,lambda3,3)*O(lambda1,1,lambda2,2))+OP(4,lambda3,3)*O(lambda1,1,lambda2,2)-OP(3,lambda2,2)*O(lambda1,1,lambda3,3)+OP(4,lambda2,2)*O(lambda1,1,lambda3,3)+OP(2,lambda1,1)*O(lambda2,2,lambda3,3)+OP(3,lambda1,1)*O(lambda2,2,lambda3,3))-pp(1,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)-pp(1,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(-(OP(2,lambda4,4)*O(lambda1,1,lambda2,2))+OP(3,lambda4,4)*O(lambda1,1,lambda2,2)+OP(3,lambda2,2)*O(lambda1,1,lambda4,4)-OP(4,lambda2,2)*O(lambda1,1,lambda4,4)+OP(2,lambda1,1)*O(lambda2,2,lambda4,4)+OP(4,lambda1,1)*O(lambda2,2,lambda4,4))+OP(1,lambda2,2)*OP(3,lambda1,1)*O(lambda3,3,lambda4,4)+OP(1,lambda2,2)*OP(4,lambda1,1)*O(lambda3,3,lambda4,4)-pp(1,3)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)-pp(1,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)))/2

               enddo                      
              enddo
             enddo 
            enddo
         endif


         if( fm2.ne.0d0 ) then  !fm2 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(fm2*(gw**2)*(sw**3)*(ZMASS**3)*(OP(1,lambda4,4)*OP(4,lambda1,1)*O(lambda2,2,lambda3,3)-pp(1,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda3,3)*OP(3,lambda1,1)*O(lambda2,2,lambda4,4)-pp(1,3)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda2,2)*OP(2,lambda1,1)*O(lambda3,3,lambda4,4)-pp(1,2)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)))/WMASS

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fm3.ne.0d0 ) then  !fm3 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(fm3*(gw**2)*(sw**3)*(ZMASS**3)*(-(OP(1,lambda2,2)*OP(2,lambda4,4)*O(lambda1,1,lambda3,3))+OP(1,lambda2,2)*OP(3,lambda4,4)*O(lambda1,1,lambda3,3)-OP(1,lambda2,2)*OP(2,lambda3,3)*O(lambda1,1,lambda4,4)+OP(1,lambda2,2)*OP(4,lambda3,3)*O(lambda1,1,lambda4,4)+pp(1,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*(OP(2,lambda3,3)*O(lambda1,1,lambda2,2)-OP(4,lambda3,3)*O(lambda1,1,lambda2,2)+OP(3,lambda2,2)*O(lambda1,1,lambda3,3)-OP(4,lambda2,2)*O(lambda1,1,lambda3,3)-OP(2,lambda1,1)*O(lambda2,2,lambda3,3)-OP(3,lambda1,1)*O(lambda2,2,lambda3,3))+pp(1,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(OP(2,lambda4,4)*O(lambda1,1,lambda2,2)-OP(3,lambda4,4)*O(lambda1,1,lambda2,2)-OP(3,lambda2,2)*O(lambda1,1,lambda4,4)+OP(4,lambda2,2)*O(lambda1,1,lambda4,4)-OP(2,lambda1,1)*O(lambda2,2,lambda4,4)-OP(4,lambda1,1)*O(lambda2,2,lambda4,4))-OP(1,lambda2,2)*OP(3,lambda1,1)*O(lambda3,3,lambda4,4)-OP(1,lambda2,2)*OP(4,lambda1,1)*O(lambda3,3,lambda4,4)+pp(1,3)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(1,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)))/(4*WMASS)

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fm4.ne.0d0 ) then  !fm4 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) -(fm4*(gw**2)*sw*(-1+2*SW2)*(ZMASS**3)*(OP(1,lambda4,4)*OP(4,lambda1,1)*O(lambda2,2,lambda3,3)-pp(1,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda3,3)*OP(3,lambda1,1)*O(lambda2,2,lambda4,4)-pp(1,3)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda2,2)*OP(2,lambda1,1)*O(lambda3,3,lambda4,4)-pp(1,2)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)))/(2*WMASS)

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fm5.ne.0d0 ) then  !fm5 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(fm5*(gw**2)*sw*(-1+2*SW2)*(ZMASS**3)*(-(OP(1,lambda2,2)*OP(2,lambda4,4)*O(lambda1,1,lambda3,3))+OP(1,lambda2,2)*OP(3,lambda4,4)*O(lambda1,1,lambda3,3)-OP(1,lambda2,2)*OP(2,lambda3,3)*O(lambda1,1,lambda4,4)+OP(1,lambda2,2)*OP(4,lambda3,3)*O(lambda1,1,lambda4,4)+pp(1,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*(OP(2,lambda3,3)*O(lambda1,1,lambda2,2)-OP(4,lambda3,3)*O(lambda1,1,lambda2,2)+OP(3,lambda2,2)*O(lambda1,1,lambda3,3)-OP(4,lambda2,2)*O(lambda1,1,lambda3,3)-OP(2,lambda1,1)*O(lambda2,2,lambda3,3)-OP(3,lambda1,1)*O(lambda2,2,lambda3,3))+pp(1,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(OP(2,lambda4,4)*O(lambda1,1,lambda2,2)-OP(3,lambda4,4)*O(lambda1,1,lambda2,2)-OP(3,lambda2,2)*O(lambda1,1,lambda4,4)+OP(4,lambda2,2)*O(lambda1,1,lambda4,4)-OP(2,lambda1,1)*O(lambda2,2,lambda4,4)-OP(4,lambda1,1)*O(lambda2,2,lambda4,4))-OP(1,lambda2,2)*OP(3,lambda1,1)*O(lambda3,3,lambda4,4)-OP(1,lambda2,2)*OP(4,lambda1,1)*O(lambda3,3,lambda4,4)+pp(1,3)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(1,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)))/(8*WMASS)

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fm5hc.ne.0d0 ) then  !fm5hc coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(fm5hc*(gw**2)*sw*(-1+2*SW2)*(ZMASS**3)*(-(OP(1,lambda2,2)*OP(2,lambda4,4)*O(lambda1,1,lambda3,3))+OP(1,lambda2,2)*OP(3,lambda4,4)*O(lambda1,1,lambda3,3)-OP(1,lambda2,2)*OP(2,lambda3,3)*O(lambda1,1,lambda4,4)+OP(1,lambda2,2)*OP(4,lambda3,3)*O(lambda1,1,lambda4,4)+pp(1,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*(OP(2,lambda3,3)*O(lambda1,1,lambda2,2)-OP(4,lambda3,3)*O(lambda1,1,lambda2,2)+OP(3,lambda2,2)*O(lambda1,1,lambda3,3)-OP(4,lambda2,2)*O(lambda1,1,lambda3,3)-OP(2,lambda1,1)*O(lambda2,2,lambda3,3)-OP(3,lambda1,1)*O(lambda2,2,lambda3,3))+pp(1,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(OP(2,lambda4,4)*O(lambda1,1,lambda2,2)-OP(3,lambda4,4)*O(lambda1,1,lambda2,2)-OP(3,lambda2,2)*O(lambda1,1,lambda4,4)+OP(4,lambda2,2)*O(lambda1,1,lambda4,4)-OP(2,lambda1,1)*O(lambda2,2,lambda4,4)-OP(4,lambda1,1)*O(lambda2,2,lambda4,4))-OP(1,lambda2,2)*OP(3,lambda1,1)*O(lambda3,3,lambda4,4)-OP(1,lambda2,2)*OP(4,lambda1,1)*O(lambda3,3,lambda4,4)+pp(1,3)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(1,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)))/(8*WMASS)

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fm6.ne.0d0 ) then  !fm6 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) -(fm6*(gw**2)*sw*WMASS*ZMASS*(OP(1,lambda4,4)*OP(4,lambda1,1)*O(lambda2,2,lambda3,3)-pp(1,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda3,3)*OP(3,lambda1,1)*O(lambda2,2,lambda4,4)-pp(1,3)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda2,2)*OP(2,lambda1,1)*O(lambda3,3,lambda4,4)-pp(1,2)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)))

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fm7.ne.0d0 ) then  !fm7 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(fm7*(gw**2)*sw*WMASS*ZMASS*(-(OP(1,lambda2,2)*OP(2,lambda4,4)*O(lambda1,1,lambda3,3))+OP(1,lambda2,2)*OP(3,lambda4,4)*O(lambda1,1,lambda3,3)-OP(1,lambda2,2)*OP(2,lambda3,3)*O(lambda1,1,lambda4,4)+OP(1,lambda2,2)*OP(4,lambda3,3)*O(lambda1,1,lambda4,4)+pp(1,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*(OP(2,lambda3,3)*O(lambda1,1,lambda2,2)-OP(4,lambda3,3)*O(lambda1,1,lambda2,2)+OP(3,lambda2,2)*O(lambda1,1,lambda3,3)-OP(4,lambda2,2)*O(lambda1,1,lambda3,3)-OP(2,lambda1,1)*O(lambda2,2,lambda3,3)-OP(3,lambda1,1)*O(lambda2,2,lambda3,3))+pp(1,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(OP(2,lambda4,4)*O(lambda1,1,lambda2,2)-OP(3,lambda4,4)*O(lambda1,1,lambda2,2)-OP(3,lambda2,2)*O(lambda1,1,lambda4,4)+OP(4,lambda2,2)*O(lambda1,1,lambda4,4)-OP(2,lambda1,1)*O(lambda2,2,lambda4,4)-OP(4,lambda1,1)*O(lambda2,2,lambda4,4))-OP(1,lambda2,2)*OP(3,lambda1,1)*O(lambda3,3,lambda4,4)-OP(1,lambda2,2)*OP(4,lambda1,1)*O(lambda3,3,lambda4,4)+pp(1,3)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(1,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)))/4

               enddo                      
              enddo
             enddo 
            enddo
         endif         

      end




ccccccccccccccc
c azwmwp-vertex c
ccccccccccccccc

      subroutine azwmwp_ONanom4(a1,z,wm1,wp1,M)
c
c This subroutine calculates the vvvv-vertex with anomalous couplings
c and the on-shell polarization vectors for vector bosons
c
c IN:
c     polv1, polv2      : the incoming polvector wz(bosons in the CM frame)
c     polv3, polv4      : the outgoing polvector wz(bosons in the CM frame)
c OUT:
c     M(lambda1,lambda2) 
c
      IMPLICIT NONE

      INTEGER NGRAPHS
      PARAMETER( NGRAPHS=20 )
      DOUBLE PRECISION rZERO, rONE, rHalf
      PARAMETER( rZERO=0.0D0, rONE=1.0D0, rHALF = 0.5D0 )
      DOUBLE COMPLEX cIMAG, cZERO, cONE
      PARAMETER( cIMAG=(0.0D0,1.0D0), cZERO=(0.0D0,0.0D0), cONE=(1.0D0,0.0D0))
      DOUBLE COMPLEX a1(6,-1:1), z(6,-1:1), wm1(6,-1:1), wp1(6,-1:1)      
      DOUBLE COMPLEX M(-1:1,-1:1,-1:1,-1:1)
      
      DOUBLE COMPLEX OFFeps(6,-1:1,NGRAPHS), OFFepsC(6,-1:1,NGRAPHS)
      DOUBLE PRECISION K(0:3,4)      
      DOUBLE COMPLEX O(-1:1,4,-1:1,4), OP(4,-1:1,4)
      DOUBLE COMPLEXPP(1:4,1:4)
    
      INTEGER I,J,N,MU,NU
      INTEGER LAMBDA1, LAMBDA2, LAMBDA3, LAMBDA4, LAMBDA
      INTEGER L1,L2

      DOUBLE COMPLEX DOTCC, DOTRC,contract_Tjj
      DOUBLE PRECISION DOTRR
      EXTERNAL DOTCC, DOTRC, DOTRR, contract_Tjj
      double precision sw2, sw
#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "gmunu.inc"
#include "transm.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"
      sw2 = sin2w
      sw = sqrt(sin2w)
      OFFeps(6,-1:1,4) = cZero
      K(0:3,1:4) = rZero
      M(-1:1,-1:1,-1:1,-1:1) = cZero
      PP(1:4,1:4) = rZero
      O(-1:1,1:4,-1:1,1:4) = cZero
      OP(1:4,-1:1,1:4) = cZero  
      
      OFFeps(1:6,-1:1,1) = a1(1:6,-1:1)
      OFFeps(1:6,-1:1,2) = z(1:6,-1:1)
      OFFeps(1:6,-1:1,3) = wm1(1:6,-1:1)
      OFFeps(1:6,-1:1,4) = wp1(1:6,-1:1) 

      do i = 1,4
        K(0,i) = dble(OFFeps(5,-1,i))*(-rONE) !K1(0:3)
        K(1,i) = dble(OFFeps(6,-1,i))*(-rONE) !K2(0:3)
        K(2,i) = dimag(OFFeps(6,-1,i))*(-rONE)!K3(0:3)
        K(3,i) = dimag(OFFeps(5,-1,i))*(-rONE)!K4(0:3)
      enddo
      



      do i = 1,4
        do j = 1,4
          pp(i,j) = dotrr(k(0,i),k(0,j))
        enddo
      enddo

      do L1 = -1,1
        do L2 = -1,1 
          do i = 1,4
            do j = 1,4
              O(L1,i,L2,j) = dotcc(OFFeps(1,L1,i),OFFeps(1,L2,j))
            enddo
          enddo
        enddo
      enddo

      do LAMBDA = -1,1
        do i = 1,4 
          do j = 1,4
            OP(j,lambda,i) = dotrc(k(0,j),OFFeps(1,lambda,i))
          enddo
        enddo
      enddo


!!!! ===============================================================      
!!!! COUPLINGS CALCULATION
!!!! ===============================================================

                 

         if( ft0.ne.0d0 ) then  !ft0 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(8*ft0*(gw**4)*sw*WMASS*(OP(1,lambda2,2)*OP(2,lambda1,1)-pp(1,2)*O(lambda1,1,lambda2,2))*(OP(3,lambda4,4)*OP(4,lambda3,3)-pp(3,4)*O(lambda3,3,lambda4,4)))/ZMASS

               enddo                      
              enddo
             enddo 
            enddo
         endif


         if( ft1.ne.0d0 ) then  !ft1 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(4*ft1*(gw**4)*sw*WMASS*(-(OP(2,lambda4,4)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda3,3))-OP(2,lambda3,3)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda4,4)+pp(1,4)*pp(3,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*OP(4,lambda1,1)*(OP(2,lambda3,3)*OP(3,lambda2,2)-pp(3,2)*O(lambda2,2,lambda3,3))+pp(1,3)*pp(4,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*OP(3,lambda1,1)*(OP(2,lambda4,4)*OP(4,lambda2,2)-pp(4,2)*O(lambda2,2,lambda4,4))))/ZMASS

               enddo                      
              enddo
             enddo 
            enddo
         endif


         if( ft2.ne.0d0 ) then  !ft2 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(ft2*(gw**4)*sw*WMASS*(OP(1,lambda2,2)*OP(2,lambda3,3)*OP(3,lambda4,4)*OP(4,lambda1,1)+OP(1,lambda2,2)*OP(2,lambda4,4)*OP(3,lambda1,1)*OP(4,lambda3,3)-OP(2,lambda4,4)*OP(4,lambda3,3)*pp(1,3)*O(lambda1,1,lambda2,2)-OP(2,lambda3,3)*OP(3,lambda4,4)*pp(1,4)*O(lambda1,1,lambda2,2)-OP(3,lambda4,4)*OP(4,lambda2,2)*pp(1,2)*O(lambda1,1,lambda3,3)-OP(1,lambda2,2)*OP(2,lambda4,4)*pp(3,4)*O(lambda1,1,lambda3,3)+OP(1,lambda2,2)*OP(3,lambda4,4)*pp(4,2)*O(lambda1,1,lambda3,3)-OP(3,lambda2,2)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda4,4)+OP(1,lambda2,2)*OP(4,lambda3,3)*pp(3,2)*O(lambda1,1,lambda4,4)-OP(1,lambda2,2)*OP(2,lambda3,3)*pp(3,4)*O(lambda1,1,lambda4,4)-OP(3,lambda4,4)*OP(4,lambda1,1)*pp(1,2)*O(lambda2,2,lambda3,3)+OP(2,lambda1,1)*OP(3,lambda4,4)*pp(1,4)*O(lambda2,2,lambda3,3)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*(OP(2,lambda1,1)*OP(3,lambda2,2)*OP(4,lambda3,3)-OP(4,lambda3,3)*pp(3,2)*O(lambda1,1,lambda2,2)+OP(2,lambda3,3)*pp(3,4)*O(lambda1,1,lambda2,2)-OP(2,lambda1,1)*pp(3,4)*O(lambda2,2,lambda3,3))-OP(3,lambda1,1)*OP(4,lambda3,3)*pp(1,2)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*OP(4,lambda3,3)*pp(1,3)*O(lambda2,2,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(OP(3,lambda4,4)*(OP(2,lambda1,1)*OP(4,lambda2,2)-pp(4,2)*O(lambda1,1,lambda2,2))+pp(3,4)*(OP(2,lambda4,4)*O(lambda1,1,lambda2,2)-OP(2,lambda1,1)*O(lambda2,2,lambda4,4)))+OP(3,lambda2,2)*OP(4,lambda1,1)*pp(1,2)*O(lambda3,3,lambda4,4)+OP(3,lambda1,1)*OP(4,lambda2,2)*pp(1,2)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(4,lambda2,2)*pp(1,3)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(3,lambda2,2)*pp(1,4)*O(lambda3,3,lambda4,4)-OP(1,lambda2,2)*OP(4,lambda1,1)*pp(3,2)*O(lambda3,3,lambda4,4)-OP(1,lambda2,2)*OP(3,lambda1,1)*pp(4,2)*O(lambda3,3,lambda4,4)+pp(1,4)*pp(3,2)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(1,3)*pp(4,2)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)))/ZMASS

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( ft3.ne.0d0 ) then  !ft3 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(1/ZMASS)*2*ft3*WMASS*gw**4*sw*(pp(1,4)*OP(3,lambda2,2)*(-OP(2,lambda4,4))*O(lambda1,1,lambda3,3)-pp(1,3)*OP(4,lambda2,2)*OP(2,lambda3,3)*O(lambda1,1,lambda4,4)-pp(1,3)*OP(4,lambda1,1)*OP(2,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,4)*OP(3,lambda1,1)*OP(2,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*(O(lambda1,1,lambda3,3)*(pp(4,2)*OP(3,lambda2,2)-pp(3,2)*OP(4,lambda2,2))+OP(3,lambda1,1)*(OP(4,lambda2,2)*OP(2,lambda3,3)-pp(4,2)*O(lambda2,2,lambda3,3)))+pp(1,3)*OP(4,lambda1,1)*OP(2,lambda3,3)*O(lambda2,2,lambda4,4)-pp(1,4)*OP(3,lambda1,1)*OP(2,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(OP(3,lambda2,2)*(OP(4,lambda1,1)*OP(2,lambda4,4)-pp(4,2)*O(lambda1,1,lambda4,4))+pp(3,2)*(OP(4,lambda2,2)*O(lambda1,1,lambda4,4)-OP(4,lambda1,1)*O(lambda2,2,lambda4,4)))+pp(1,4)*pp(3,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,3)*pp(4,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3))
                  
               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( ft4.ne.0d0 ) then  !ft4 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(gw**4*ft4*sw**3*ZMASS*(pp(1, 4)*OP(2, lambda4, 4)*OP(3, lambda2, 2)*O(lambda1, 1, lambda3, 3) + pp(1, 3)*OP(2, lambda3, 3)*OP(4, lambda2, 2)*O(lambda1, 1, lambda4, 4) - pp(1, 4)*OP(2, lambda4, 4)*OP(3, lambda1, 1)*O(lambda2, 2, lambda3, 3) + pp(1, 3)*OP(2, lambda4, 4)*OP(4, lambda1, 1)*O(lambda2, 2, lambda3, 3) - pp(1, 3)*pp(4, 2)*O(lambda1, 1, lambda4, 4)*O(lambda2, 2, lambda3, 3) + OP(1, lambda4, 4)*((-(pp(4, 2)*OP(3, lambda2, 2)) + pp(3, 2)*OP(4, lambda2, 2))*O(lambda1, 1, lambda3, 3) +   OP(3, lambda1, 1)*(-(OP(2, lambda3, 3)*OP(4, lambda2, 2)) + pp(4, 2)*O(lambda2, 2, lambda3, 3))) + pp(1, 4)*OP(2, lambda3, 3)*OP(3, lambda1, 1)*O(lambda2, 2, lambda4, 4) - pp(1, 3)*OP(2, lambda3, 3)*OP(4, lambda1, 1)*O(lambda2, 2, lambda4, 4) - pp(1, 4)*pp(3, 2)*O(lambda1, 1, lambda3, 3)*O(lambda2, 2, lambda4, 4) + OP(1, lambda3, 3)*(OP(3, lambda2, 2)*(-(OP(2, lambda4, 4)*OP(4, lambda1, 1)) + pp(4, 2)*O(lambda1, 1, lambda4, 4)) +   pp(3, 2)*(-(OP(4, lambda2, 2)*O(lambda1, 1, lambda4, 4)) + OP(4, lambda1, 1)*O(lambda2, 2, lambda4, 4)))))/(2*WMASS)
               enddo                      
              enddo
             enddo 
            enddo
         endif 

         
         if( ft5.ne.0d0 ) then  !ft5 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(-2*ft5*(gw**4)*(sw**3)*ZMASS*(OP(1,lambda2,2)*OP(2,lambda1,1)-pp(1,2)*O(lambda1,1,lambda2,2))*(OP(3,lambda4,4)*OP(4,lambda3,3)-pp(3,4)*O(lambda3,3,lambda4,4)))/WMASS

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( ft6.ne.0d0 ) then  !ft6 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(ft6*(gw**4)*(sw**3)*ZMASS*(OP(2,lambda4,4)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda3,3)+OP(2,lambda3,3)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda4,4)-pp(1,4)*pp(3,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*OP(4,lambda1,1)*(-(OP(2,lambda3,3)*OP(3,lambda2,2))+pp(3,2)*O(lambda2,2,lambda3,3))-pp(1,3)*pp(4,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*OP(3,lambda1,1)*(-(OP(2,lambda4,4)*OP(4,lambda2,2))+pp(4,2)*O(lambda2,2,lambda4,4))))/WMASS

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( ft7.ne.0d0 ) then  !ft7 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) -(ft7*(gw**4)*(sw**3)*ZMASS*(OP(1,lambda2,2)*OP(2,lambda3,3)*OP(3,lambda4,4)*OP(4,lambda1,1)+OP(1,lambda2,2)*OP(2,lambda4,4)*OP(3,lambda1,1)*OP(4,lambda3,3)-OP(2,lambda4,4)*OP(4,lambda3,3)*pp(1,3)*O(lambda1,1,lambda2,2)-OP(2,lambda3,3)*OP(3,lambda4,4)*pp(1,4)*O(lambda1,1,lambda2,2)-OP(3,lambda4,4)*OP(4,lambda2,2)*pp(1,2)*O(lambda1,1,lambda3,3)-OP(1,lambda2,2)*OP(2,lambda4,4)*pp(3,4)*O(lambda1,1,lambda3,3)+OP(1,lambda2,2)*OP(3,lambda4,4)*pp(4,2)*O(lambda1,1,lambda3,3)-OP(3,lambda2,2)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda4,4)+OP(1,lambda2,2)*OP(4,lambda3,3)*pp(3,2)*O(lambda1,1,lambda4,4)-OP(1,lambda2,2)*OP(2,lambda3,3)*pp(3,4)*O(lambda1,1,lambda4,4)-OP(3,lambda4,4)*OP(4,lambda1,1)*pp(1,2)*O(lambda2,2,lambda3,3)+OP(2,lambda1,1)*OP(3,lambda4,4)*pp(1,4)*O(lambda2,2,lambda3,3)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*(OP(2,lambda1,1)*OP(3,lambda2,2)*OP(4,lambda3,3)-OP(4,lambda3,3)*pp(3,2)*O(lambda1,1,lambda2,2)+OP(2,lambda3,3)*pp(3,4)*O(lambda1,1,lambda2,2)-OP(2,lambda1,1)*pp(3,4)*O(lambda2,2,lambda3,3))-OP(3,lambda1,1)*OP(4,lambda3,3)*pp(1,2)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*OP(4,lambda3,3)*pp(1,3)*O(lambda2,2,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(OP(3,lambda4,4)*(OP(2,lambda1,1)*OP(4,lambda2,2)-pp(4,2)*O(lambda1,1,lambda2,2))+pp(3,4)*(OP(2,lambda4,4)*O(lambda1,1,lambda2,2)-OP(2,lambda1,1)*O(lambda2,2,lambda4,4)))+OP(3,lambda2,2)*OP(4,lambda1,1)*pp(1,2)*O(lambda3,3,lambda4,4)+OP(3,lambda1,1)*OP(4,lambda2,2)*pp(1,2)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(4,lambda2,2)*pp(1,3)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(3,lambda2,2)*pp(1,4)*O(lambda3,3,lambda4,4)-OP(1,lambda2,2)*OP(4,lambda1,1)*pp(3,2)*O(lambda3,3,lambda4,4)-OP(1,lambda2,2)*OP(3,lambda1,1)*pp(4,2)*O(lambda3,3,lambda4,4)+pp(1,4)*pp(3,2)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(1,3)*pp(4,2)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)))/(4*WMASS)

               enddo                      
              enddo
             enddo 
            enddo
         endif



         if( fm0.ne.0d0 ) then  !fm0 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(-2*fm0*(gw**2)*sw*(WMASS**3)*(OP(1,lambda2,2)*OP(2,lambda1,1)-pp(1,2)*O(lambda1,1,lambda2,2))*O(lambda3,3,lambda4,4))/ZMASS

               enddo                      
              enddo
             enddo 
            enddo
         endif


         if( fm1.ne.0d0 ) then  !fm1 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(fm1*(gw**2)*sw*(WMASS**3)*(-(OP(1,lambda3,3)*OP(2,lambda4,4)*O(lambda1,1,lambda2,2))+OP(1,lambda2,2)*(OP(2,lambda4,4)*O(lambda1,1,lambda3,3)+OP(2,lambda3,3)*O(lambda1,1,lambda4,4))-pp(1,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*(-(OP(2,lambda3,3)*O(lambda1,1,lambda2,2))+OP(2,lambda1,1)*O(lambda2,2,lambda3,3))+OP(1,lambda3,3)*OP(2,lambda1,1)*O(lambda2,2,lambda4,4)-pp(1,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)))/(2*ZMASS)

               enddo                      
              enddo
             enddo 
            enddo
         endif


         if( fm2.ne.0d0 ) then  !fm2 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +fm2*(gw**2)*(sw**3)*WMASS*ZMASS*(OP(1,lambda2,2)*OP(2,lambda1,1)-pp(1,2)*O(lambda1,1,lambda2,2))*O(lambda3,3,lambda4,4)

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fm3.ne.0d0 ) then  !fm3 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) -(fm3*(gw**2)*(sw**3)*WMASS*ZMASS*(-(OP(1,lambda3,3)*OP(2,lambda4,4)*O(lambda1,1,lambda2,2))+OP(1,lambda2,2)*(OP(2,lambda4,4)*O(lambda1,1,lambda3,3)+OP(2,lambda3,3)*O(lambda1,1,lambda4,4))-pp(1,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*(-(OP(2,lambda3,3)*O(lambda1,1,lambda2,2))+OP(2,lambda1,1)*O(lambda2,2,lambda3,3))+OP(1,lambda3,3)*OP(2,lambda1,1)*O(lambda2,2,lambda4,4)-pp(1,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)))/4

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fm4.ne.0d0 ) then  !fm4 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(fm4*(gw**2)*sw*WMASS*ZMASS*(OP(1,lambda4,4)*OP(4,lambda1,1)*O(lambda2,2,lambda3,3)-pp(1,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda3,3)*OP(3,lambda1,1)*O(lambda2,2,lambda4,4)-pp(1,3)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+(-1+2*SW2)*OP(1,lambda2,2)*OP(2,lambda1,1)*O(lambda3,3,lambda4,4)+pp(1,2)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)-2*SW2*pp(1,2)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)))/2

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fm5.ne.0d0 ) then  !fm5 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) -(fm5*(gw**2)*sw*WMASS*ZMASS*(-(OP(1,lambda3,3)*OP(3,lambda2,2)*O(lambda1,1,lambda4,4))-pp(1,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*(-(OP(2,lambda3,3)*O(lambda1,1,lambda2,2))-OP(4,lambda3,3)*O(lambda1,1,lambda2,2)+OP(3,lambda2,2)*O(lambda1,1,lambda3,3)+OP(2,lambda1,1)*O(lambda2,2,lambda3,3)-OP(3,lambda1,1)*O(lambda2,2,lambda3,3)+SW2*(OP(2,lambda3,3)*O(lambda1,1,lambda2,2)-OP(2,lambda1,1)*O(lambda2,2,lambda3,3)))+SW2*(OP(1,lambda3,3)*(OP(2,lambda4,4)*O(lambda1,1,lambda2,2)-OP(2,lambda1,1)*O(lambda2,2,lambda4,4))+pp(1,2)*(O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)))+pp(1,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)-OP(1,lambda2,2)*(-(OP(2,lambda3,3)*O(lambda1,1,lambda4,4))-OP(4,lambda3,3)*O(lambda1,1,lambda4,4)+SW2*(OP(2,lambda4,4)*O(lambda1,1,lambda3,3)+OP(2,lambda3,3)*O(lambda1,1,lambda4,4))+OP(4,lambda1,1)*O(lambda3,3,lambda4,4))))/4

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fm5hc.ne.0d0 ) then  !fm5hc coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) -(fm5hc*(gw**2)*sw*WMASS*ZMASS*(-(OP(1,lambda3,3)*OP(2,lambda4,4)*O(lambda1,1,lambda2,2))-OP(1,lambda3,3)*OP(3,lambda4,4)*O(lambda1,1,lambda2,2)+OP(1,lambda3,3)*OP(4,lambda2,2)*O(lambda1,1,lambda4,4)-OP(1,lambda4,4)*(OP(4,lambda2,2)*O(lambda1,1,lambda3,3)+SW2*(-(OP(2,lambda3,3)*O(lambda1,1,lambda2,2))+OP(2,lambda1,1)*O(lambda2,2,lambda3,3)))+OP(1,lambda3,3)*OP(2,lambda1,1)*O(lambda2,2,lambda4,4)-OP(1,lambda3,3)*OP(4,lambda1,1)*O(lambda2,2,lambda4,4)-pp(1,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+SW2*(OP(1,lambda3,3)*(OP(2,lambda4,4)*O(lambda1,1,lambda2,2)-OP(2,lambda1,1)*O(lambda2,2,lambda4,4))+pp(1,2)*(O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)))+pp(1,3)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+OP(1,lambda2,2)*(-((-1+SW2)*OP(2,lambda4,4)*O(lambda1,1,lambda3,3))+OP(3,lambda4,4)*O(lambda1,1,lambda3,3)-SW2*OP(2,lambda3,3)*O(lambda1,1,lambda4,4)-OP(3,lambda1,1)*O(lambda3,3,lambda4,4))))/4

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fm6.ne.0d0 ) then  !fm6 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) -((fm6*(gw**2)*sw*(WMASS**3)*(OP(1,lambda2,2)*OP(2,lambda1,1)-pp(1,2)*O(lambda1,1,lambda2,2))*O(lambda3,3,lambda4,4))/ZMASS)

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fm7.ne.0d0 ) then  !fm7 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) -(fm7*(gw**2)*sw*WMASS*ZMASS*(-(OP(1,lambda3,3)*OP(2,lambda4,4)*O(lambda1,1,lambda2,2))-OP(1,lambda3,3)*OP(3,lambda4,4)*O(lambda1,1,lambda2,2)+OP(1,lambda3,3)*OP(3,lambda2,2)*O(lambda1,1,lambda4,4)-OP(1,lambda3,3)*OP(4,lambda2,2)*O(lambda1,1,lambda4,4)-pp(1,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-pp(1,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*(-(OP(2,lambda3,3)*O(lambda1,1,lambda2,2))-OP(4,lambda3,3)*O(lambda1,1,lambda2,2)-OP(3,lambda2,2)*O(lambda1,1,lambda3,3)+OP(4,lambda2,2)*O(lambda1,1,lambda3,3)+OP(2,lambda1,1)*O(lambda2,2,lambda3,3)+OP(3,lambda1,1)*O(lambda2,2,lambda3,3)+SW2*(OP(2,lambda3,3)*O(lambda1,1,lambda2,2)-OP(2,lambda1,1)*O(lambda2,2,lambda3,3)))+OP(1,lambda3,3)*OP(2,lambda1,1)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*OP(4,lambda1,1)*O(lambda2,2,lambda4,4)-pp(1,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)-pp(1,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+SW2*(OP(1,lambda3,3)*(OP(2,lambda4,4)*O(lambda1,1,lambda2,2)-OP(2,lambda1,1)*O(lambda2,2,lambda4,4))+pp(1,2)*(O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)))+pp(1,3)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(1,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+OP(1,lambda2,2)*(-((-1+SW2)*OP(2,lambda4,4)*O(lambda1,1,lambda3,3))+OP(3,lambda4,4)*O(lambda1,1,lambda3,3)+OP(2,lambda3,3)*O(lambda1,1,lambda4,4)-SW2*OP(2,lambda3,3)*O(lambda1,1,lambda4,4)+OP(4,lambda3,3)*O(lambda1,1,lambda4,4)-OP(3,lambda1,1)*O(lambda3,3,lambda4,4)-OP(4,lambda1,1)*O(lambda3,3,lambda4,4))))/4

               enddo                      
              enddo
             enddo 
            enddo
         endif         

      end




ccccccccccccccc
c hhhh-vertex c
ccccccccccccccc

      subroutine hhhh_ONanom4(h1,h2,h3,h4,M)
c
c This subroutine calculates the vvvv-vertex with anomalous couplings
c and the on-shell polarization vectors for vector bosons
c
c IN:
c     polv1, polv2      : the incoming polvector wz(bosons in the CM frame)
c     polv3, polv4      : the outgoing polvector wz(bosons in the CM frame)
c OUT:
c     M(lambda1,lambda2) 
c
      IMPLICIT NONE

      INTEGER NGRAPHS
      PARAMETER( NGRAPHS=20 )
      DOUBLE PRECISION rZERO, rONE, rHalf
      PARAMETER( rZERO=0.0D0, rONE=1.0D0, rHALF = 0.5D0 )
      DOUBLE COMPLEX cIMAG, cZERO, cONE
      PARAMETER( cIMAG=(0.0D0,1.0D0), cZERO=(0.0D0,0.0D0), cONE=(1.0D0,0.0D0))
      DOUBLE COMPLEX h1(6,-1:1), h2(6,-1:1), h3(6,-1:1), h4(6,-1:1)      
      DOUBLE COMPLEX M(-1:1,-1:1,-1:1,-1:1)
      
      DOUBLE COMPLEX OFFeps(6,-1:1,NGRAPHS), OFFepsC(6,-1:1,NGRAPHS)
      DOUBLE PRECISION K(0:3,4)      
      DOUBLE COMPLEX O(-1:1,4,-1:1,4), OP(4,-1:1,4)
      DOUBLE COMPLEXPP(1:4,1:4)
    
      INTEGER I,J,N,MU,NU
      INTEGER LAMBDA1, LAMBDA2, LAMBDA3, LAMBDA4, LAMBDA
      INTEGER L1,L2

      DOUBLE COMPLEX DOTCC, DOTRC,contract_Tjj
      DOUBLE PRECISION DOTRR
      EXTERNAL DOTCC, DOTRC, DOTRR, contract_Tjj
      double precision sw2, sw
#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "gmunu.inc"
#include "transm.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"
      sw2 = sin2w
      sw = sqrt(sin2w)
      OFFeps(6,-1:1,4) = cZero
      K(0:3,1:4) = rZero
      M(-1:1,-1:1,-1:1,-1:1) = cZero
      PP(1:4,1:4) = rZero
      O(-1:1,1:4,-1:1,1:4) = cZero
      OP(1:4,-1:1,1:4) = cZero  
      
      OFFeps(1:6,-1:1,1) = h1(1:6,-1:1)
      OFFeps(1:6,-1:1,2) = h2(1:6,-1:1)
      OFFeps(1:6,-1:1,3) = h3(1:6,-1:1)
      OFFeps(1:6,-1:1,4) = h4(1:6,-1:1) 

      do i = 1,4
        K(0,i) = dble(OFFeps(5,-1,i))*(-rONE) !K1(0:3)
        K(1,i) = dble(OFFeps(6,-1,i))*(-rONE) !K2(0:3)
        K(2,i) = dimag(OFFeps(6,-1,i))*(-rONE)!K3(0:3)
        K(3,i) = dimag(OFFeps(5,-1,i))*(-rONE)!K4(0:3)
      enddo
      



      do i = 1,4
        do j = 1,4
          pp(i,j) = dotrr(k(0,i),k(0,j))
        enddo
      enddo

      do L1 = -1,1
        do L2 = -1,1 
          do i = 1,4
            do j = 1,4
              O(L1,i,L2,j) = dotcc(OFFeps(1,L1,i),OFFeps(1,L2,j))
            enddo
          enddo
        enddo
      enddo

      do LAMBDA = -1,1
        do i = 1,4 
          do j = 1,4
            OP(j,lambda,i) = dotrc(k(0,j),OFFeps(1,lambda,i))
          enddo
        enddo
      enddo


!!!! ===============================================================      
!!!! COUPLINGS CALCULATION
!!!! ===============================================================



         if( fs0.ne.0d0 ) then  !fs0 coupling
!            do lambda1 = -1,1
!             do lambda2 = -1,1
!              do lambda3 = -1,1
!               do lambda4 = -1,1
                  M(0,0,0,0) = M(0,0,0,0) +2*(fs0*(pp(1,4)*pp(2,3)+pp(1,3)*pp(2,4)+pp(1,2)*pp(3,4)))
!               enddo                      
!              enddo
!             enddo 
!            enddo
         endif

         if( fs1.ne.0d0 ) then  !fs1 coupling
!            do lambda1 = -1,1
!             do lambda2 = -1,1
!              do lambda3 = -1,1
!               do lambda4 = -1,1
                  M(0,0,0,0) = M(0,0,0,0) +2*(fs1*(pp(1,4)*pp(2,3)+pp(1,3)*pp(2,4)+pp(1,2)*pp(3,4)))
!               enddo                      
!              enddo
!             enddo 
!            enddo
         endif

         if( fs2.ne.0d0 ) then  !fs2 coupling
!            do lambda1 = -1,1
!             do lambda2 = -1,1
!              do lambda3 = -1,1
!               do lambda4 = -1,1
                  M(0,0,0,0) = M(0,0,0,0) +2*(fs2*(pp(1,4)*pp(2,3)+pp(1,3)*pp(2,4)+pp(1,2)*pp(3,4)))

!               enddo                      
!              enddo
!             enddo 
!            enddo
         endif                  

     

      end




ccccccccccccccc
c hhhz-vertex c
ccccccccccccccc

      subroutine hhhz_ONanom4(h1,h2,h3,z,M)
c
c This subroutine calculates the vvvv-vertex with anomalous couplings
c and the on-shell polarization vectors for vector bosons
c
c IN:
c     polv1, polv2      : the incoming polvector wz(bosons in the CM frame)
c     polv3, polv4      : the outgoing polvector wz(bosons in the CM frame)
c OUT:
c     M(lambda1,lambda2) 
c
      IMPLICIT NONE

      INTEGER NGRAPHS
      PARAMETER( NGRAPHS=20 )
      DOUBLE PRECISION rZERO, rONE, rHalf
      PARAMETER( rZERO=0.0D0, rONE=1.0D0, rHALF = 0.5D0 )
      DOUBLE COMPLEX cIMAG, cZERO, cONE
      PARAMETER( cIMAG=(0.0D0,1.0D0), cZERO=(0.0D0,0.0D0), cONE=(1.0D0,0.0D0))
      DOUBLE COMPLEX h1(6,-1:1), h2(6,-1:1), h3(6,-1:1), z(6,-1:1)      
      DOUBLE COMPLEX M(-1:1,-1:1,-1:1,-1:1)
      
      DOUBLE COMPLEX OFFeps(6,-1:1,NGRAPHS), OFFepsC(6,-1:1,NGRAPHS)
      DOUBLE PRECISION K(0:3,4)      
      DOUBLE COMPLEX O(-1:1,4,-1:1,4), OP(4,-1:1,4)
      DOUBLE COMPLEXPP(1:4,1:4)
    
      INTEGER I,J,N,MU,NU
      INTEGER LAMBDA1, LAMBDA2, LAMBDA3, LAMBDA4, LAMBDA
      INTEGER L1,L2

      DOUBLE COMPLEX DOTCC, DOTRC,contract_Tjj
      DOUBLE PRECISION DOTRR
      EXTERNAL DOTCC, DOTRC, DOTRR, contract_Tjj
      double precision sw2, sw
#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "gmunu.inc"
#include "transm.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"
      sw2 = sin2w
      sw = sqrt(sin2w)
      OFFeps(6,-1:1,4) = cZero
      K(0:3,1:4) = rZero
      M(-1:1,-1:1,-1:1,-1:1) = cZero
      PP(1:4,1:4) = rZero
      O(-1:1,1:4,-1:1,1:4) = cZero
      OP(1:4,-1:1,1:4) = cZero  
      
      OFFeps(1:6,-1:1,1) = h1(1:6,-1:1)
      OFFeps(1:6,-1:1,2) = h2(1:6,-1:1)
      OFFeps(1:6,-1:1,3) = h3(1:6,-1:1)
      OFFeps(1:6,-1:1,4) = z(1:6,-1:1) 

      do i = 1,4
        K(0,i) = dble(OFFeps(5,-1,i))*(-rONE) !K1(0:3)
        K(1,i) = dble(OFFeps(6,-1,i))*(-rONE) !K2(0:3)
        K(2,i) = dimag(OFFeps(6,-1,i))*(-rONE)!K3(0:3)
        K(3,i) = dimag(OFFeps(5,-1,i))*(-rONE)!K4(0:3)
      enddo
      



      do i = 1,4
        do j = 1,4
          pp(i,j) = dotrr(k(0,i),k(0,j))
        enddo
      enddo

      do L1 = -1,1
        do L2 = -1,1 
          do i = 1,4
            do j = 1,4
              O(L1,i,L2,j) = dotcc(OFFeps(1,L1,i),OFFeps(1,L2,j))
            enddo
          enddo
        enddo
      enddo

      do LAMBDA = -1,1
        do i = 1,4 
          do j = 1,4
            OP(j,lambda,i) = dotrc(k(0,j),OFFeps(1,lambda,i))
          enddo
        enddo
      enddo


!!!! ===============================================================      
!!!! COUPLINGS CALCULATION
!!!! ===============================================================





      end


ccccccccccccccc
c hhzz-vertex c
ccccccccccccccc

      subroutine hhzz_ONanom4(h1,h2,z1,z2,M)
c
c This subroutine calculates the vvvv-vertex with anomalous couplings
c and the on-shell polarization vectors for vector bosons
c
c IN:
c     polv1, polv2      : the incoming polvector wz(bosons in the CM frame)
c     polv3, polv4      : the outgoing polvector wz(bosons in the CM frame)
c OUT:
c     M(lambda1,lambda2) 
c
      IMPLICIT NONE

      INTEGER NGRAPHS
      PARAMETER( NGRAPHS=20 )
      DOUBLE PRECISION rZERO, rONE, rHalf
      PARAMETER( rZERO=0.0D0, rONE=1.0D0, rHALF = 0.5D0 )
      DOUBLE COMPLEX cIMAG, cZERO, cONE
      PARAMETER( cIMAG=(0.0D0,1.0D0), cZERO=(0.0D0,0.0D0), cONE=(1.0D0,0.0D0))
      DOUBLE COMPLEX h1(6,-1:1), h2(6,-1:1), z1(6,-1:1), z2(6,-1:1)      
      DOUBLE COMPLEX M(-1:1,-1:1,-1:1,-1:1)
      
      DOUBLE COMPLEX OFFeps(6,-1:1,NGRAPHS), OFFepsC(6,-1:1,NGRAPHS)
      DOUBLE PRECISION K(0:3,4)      
      DOUBLE COMPLEX O(-1:1,4,-1:1,4), OP(4,-1:1,4)
      DOUBLE COMPLEXPP(1:4,1:4)
    
      INTEGER I,J,N,MU,NU
      INTEGER LAMBDA1, LAMBDA2, LAMBDA3, LAMBDA4, LAMBDA
      INTEGER L1,L2

      DOUBLE COMPLEX DOTCC, DOTRC,contract_Tjj
      DOUBLE PRECISION DOTRR
      EXTERNAL DOTCC, DOTRC, DOTRR, contract_Tjj
      double precision sw2, sw
#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "gmunu.inc"
#include "transm.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"
      sw2 = sin2w
      sw = sqrt(sin2w)
      OFFeps(6,-1:1,4) = cZero
      K(0:3,1:4) = rZero
      M(-1:1,-1:1,-1:1,-1:1) = cZero
      PP(1:4,1:4) = rZero
      O(-1:1,1:4,-1:1,1:4) = cZero
      OP(1:4,-1:1,1:4) = cZero  
      
      OFFeps(1:6,-1:1,1) = h1(1:6,-1:1)
      OFFeps(1:6,-1:1,2) = h2(1:6,-1:1)
      OFFeps(1:6,-1:1,3) = z1(1:6,-1:1)
      OFFeps(1:6,-1:1,4) = z2(1:6,-1:1) 

      do i = 1,4
        K(0,i) = dble(OFFeps(5,-1,i))*(-rONE) !K1(0:3)
        K(1,i) = dble(OFFeps(6,-1,i))*(-rONE) !K2(0:3)
        K(2,i) = dimag(OFFeps(6,-1,i))*(-rONE)!K3(0:3)
        K(3,i) = dimag(OFFeps(5,-1,i))*(-rONE)!K4(0:3)
      enddo
      



      do i = 1,4
        do j = 1,4
          pp(i,j) = dotrr(k(0,i),k(0,j))
        enddo
      enddo

      do L1 = -1,1
        do L2 = -1,1 
          do i = 1,4
            do j = 1,4
              O(L1,i,L2,j) = dotcc(OFFeps(1,L1,i),OFFeps(1,L2,j))
            enddo
          enddo
        enddo
      enddo

      do LAMBDA = -1,1
        do i = 1,4 
          do j = 1,4
            OP(j,lambda,i) = dotrc(k(0,j),OFFeps(1,lambda,i))
          enddo
        enddo
      enddo


!!!! ===============================================================      
!!!! COUPLINGS CALCULATION
!!!! ===============================================================



         if( fs0.ne.0d0 ) then  !fs0 coupling
!            do lambda1 = -1,1
!             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(0,0,lambda3,lambda4) = M(0,0,lambda3,lambda4) -2*fs0*(ZMASS**2)*1.0*1.0*(OP(1,lambda4,4)*OP(2,lambda3,3)+OP(1,lambda3,3)*OP(2,lambda4,4)-pp(1,2)*O(lambda3,3,lambda4,4))

               enddo                      
              enddo
!             enddo 
!            enddo
         endif

         if( fs1.ne.0d0 ) then  !fs1 coupling
!            do lambda1 = -1,1
!             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(0,0,lambda3,lambda4) = M(0,0,lambda3,lambda4) -2*fs1*(ZMASS**2)*1.0*1.0*pp(1,2)*O(lambda3,3,lambda4,4)

               enddo                      
              enddo
!             enddo 
!            enddo
         endif

         if( fs2.ne.0d0 ) then  !fs2 coupling
!            do lambda1 = -1,1
!             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(0,0,lambda3,lambda4) = M(0,0,lambda3,lambda4) -2*fs2*(ZMASS**2)*1.0*1.0*pp(1,2)*O(lambda3,3,lambda4,4)

               enddo                      
              enddo
!             enddo 
!            enddo
         endif                  


         if( fm0.ne.0d0 ) then  !fm0 coupling
!            do lambda1 = -1,1
!             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(0,0,lambda3,lambda4) = M(0,0,lambda3,lambda4) -2*fm0*(gw**2)*(-1+SW2)*1.0*1.0*pp(1,2)*(OP(3,lambda4,4)*OP(4,lambda3,3)-pp(3,4)*O(lambda3,3,lambda4,4))

               enddo                      
              enddo
!             enddo 
!            enddo
         endif


         if( fm1.ne.0d0 ) then  !fm1 coupling
!            do lambda1 = -1,1
!             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(0,0,lambda3,lambda4) = M(0,0,lambda3,lambda4) +(fm1*(gw**2)*(-1+SW2)*1.0*1.0*(OP(1,lambda4,4)*OP(4,lambda3,3)*pp(2,3)+OP(1,lambda3,3)*OP(3,lambda4,4)*pp(2,4)+OP(2,lambda4,4)*(OP(4,lambda3,3)*pp(1,3)-OP(1,lambda3,3)*pp(3,4))+OP(2,lambda3,3)*(OP(3,lambda4,4)*pp(1,4)-OP(1,lambda4,4)*pp(3,4))-pp(1,4)*pp(2,3)*O(lambda3,3,lambda4,4)-pp(1,3)*pp(2,4)*O(lambda3,3,lambda4,4)))/2

               enddo                      
              enddo
!             enddo 
!            enddo
         endif


         if( fm2.ne.0d0 ) then  !fm2 coupling
!            do lambda1 = -1,1
!             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(0,0,lambda3,lambda4) = M(0,0,lambda3,lambda4) +(fm2*(gw**2)*(sw**2)*SW2*(ZMASS**2)*1.0*1.0*pp(1,2)*(OP(3,lambda4,4)*OP(4,lambda3,3)-pp(3,4)*O(lambda3,3,lambda4,4)))/(WMASS**2)

               enddo                      
              enddo
!             enddo 
!            enddo
         endif

         if( fm3.ne.0d0 ) then  !fm3 coupling
!            do lambda1 = -1,1
!             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(0,0,lambda3,lambda4) = M(0,0,lambda3,lambda4) -(fm3*(gw**2)*(sw**2)*SW2*(ZMASS**2)*1.0*1.0*(OP(1,lambda4,4)*OP(4,lambda3,3)*pp(2,3)+OP(1,lambda3,3)*OP(3,lambda4,4)*pp(2,4)+OP(2,lambda4,4)*(OP(4,lambda3,3)*pp(1,3)-OP(1,lambda3,3)*pp(3,4))+OP(2,lambda3,3)*(OP(3,lambda4,4)*pp(1,4)-OP(1,lambda4,4)*pp(3,4))-pp(1,4)*pp(2,3)*O(lambda3,3,lambda4,4)-pp(1,3)*pp(2,4)*O(lambda3,3,lambda4,4)))/(4*(WMASS**2))

               enddo                      
              enddo
!             enddo 
!            enddo
         endif

         if( fm4.ne.0d0 ) then  !fm4 coupling
!            do lambda1 = -1,1
!             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(0,0,lambda3,lambda4) = M(0,0,lambda3,lambda4) +fm4*(gw**2)*(sw**2)*1.0*1.0*pp(1,2)*(OP(3,lambda4,4)*OP(4,lambda3,3)-pp(3,4)*O(lambda3,3,lambda4,4))

               enddo                      
              enddo
!             enddo 
!            enddo
         endif

         if( fm5.ne.0d0 ) then  !fm5 coupling
!            do lambda1 = -1,1
!             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(0,0,lambda3,lambda4) = M(0,0,lambda3,lambda4) +(fm5*(gw**2)*(sw**2)*1.0*1.0*(OP(1,lambda4,4)*OP(4,lambda3,3)*pp(2,3)+OP(1,lambda3,3)*OP(3,lambda4,4)*pp(2,4)+OP(2,lambda4,4)*(OP(4,lambda3,3)*pp(1,3)-OP(1,lambda3,3)*pp(3,4))+OP(2,lambda3,3)*(OP(3,lambda4,4)*pp(1,4)-OP(1,lambda4,4)*pp(3,4))-pp(1,4)*pp(2,3)*O(lambda3,3,lambda4,4)-pp(1,3)*pp(2,4)*O(lambda3,3,lambda4,4)))/4

               enddo                      
              enddo
!             enddo 
!            enddo
         endif

         if( fm5hc.ne.0d0 ) then  !fm5hc coupling
!            do lambda1 = -1,1
!             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(0,0,lambda3,lambda4) = M(0,0,lambda3,lambda4) +(fm5hc*(gw**2)*(sw**2)*1.0*1.0*(OP(1,lambda4,4)*OP(4,lambda3,3)*pp(2,3)+OP(1,lambda3,3)*OP(3,lambda4,4)*pp(2,4)+OP(2,lambda4,4)*(OP(4,lambda3,3)*pp(1,3)-OP(1,lambda3,3)*pp(3,4))+OP(2,lambda3,3)*(OP(3,lambda4,4)*pp(1,4)-OP(1,lambda4,4)*pp(3,4))-pp(1,4)*pp(2,3)*O(lambda3,3,lambda4,4)-pp(1,3)*pp(2,4)*O(lambda3,3,lambda4,4)))/4

               enddo                      
              enddo
!             enddo 
!            enddo
         endif

         if( fm6.ne.0d0 ) then  !fm6 coupling
!            do lambda1 = -1,1
!             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(0,0,lambda3,lambda4) = M(0,0,lambda3,lambda4) -(fm6*(gw**2)*(-1+SW2)*1.0*1.0*pp(1,2)*(OP(3,lambda4,4)*OP(4,lambda3,3)-pp(3,4)*O(lambda3,3,lambda4,4)))

               enddo                      
              enddo
!             enddo 
!            enddo
         endif

         if( fm7.ne.0d0 ) then  !fm7 coupling
!            do lambda1 = -1,1
!             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(0,0,lambda3,lambda4) = M(0,0,lambda3,lambda4) -(fm7*(gw**2)*(-1+SW2)*1.0*1.0*(OP(1,lambda4,4)*OP(4,lambda3,3)*pp(2,3)+OP(1,lambda3,3)*OP(3,lambda4,4)*pp(2,4)+OP(2,lambda4,4)*(OP(4,lambda3,3)*pp(1,3)-OP(1,lambda3,3)*pp(3,4))+OP(2,lambda3,3)*(OP(3,lambda4,4)*pp(1,4)-OP(1,lambda4,4)*pp(3,4))-pp(1,4)*pp(2,3)*O(lambda3,3,lambda4,4)-pp(1,3)*pp(2,4)*O(lambda3,3,lambda4,4)))/4

               enddo                      
              enddo
!             enddo 
!            enddo
         endif         

      end





ccccccccccccccc
c hhwmwp-vertex c
ccccccccccccccc

      subroutine hhwmwp_ONanom4(h1,h2,wm1,wp1,M)
c
c This subroutine calculates the vvvv-vertex with anomalous couplings
c and the on-shell polarization vectors for vector bosons
c
c IN:
c     polv1, polv2      : the incoming polvector wz(bosons in the CM frame)
c     polv3, polv4      : the outgoing polvector wz(bosons in the CM frame)
c OUT:
c     M(lambda1,lambda2) 
c
      IMPLICIT NONE

      INTEGER NGRAPHS
      PARAMETER( NGRAPHS=20 )
      DOUBLE PRECISION rZERO, rONE, rHalf
      PARAMETER( rZERO=0.0D0, rONE=1.0D0, rHALF = 0.5D0 )
      DOUBLE COMPLEX cIMAG, cZERO, cONE
      PARAMETER( cIMAG=(0.0D0,1.0D0), cZERO=(0.0D0,0.0D0), cONE=(1.0D0,0.0D0))
      DOUBLE COMPLEX h1(6,-1:1), h2(6,-1:1), wm1(6,-1:1), wp1(6,-1:1)      
      DOUBLE COMPLEX M(-1:1,-1:1,-1:1,-1:1)
      
      DOUBLE COMPLEX OFFeps(6,-1:1,NGRAPHS), OFFepsC(6,-1:1,NGRAPHS)
      DOUBLE PRECISION K(0:3,4)      
      DOUBLE COMPLEX O(-1:1,4,-1:1,4), OP(4,-1:1,4)
      DOUBLE COMPLEXPP(1:4,1:4)
    
      INTEGER I,J,N,MU,NU
      INTEGER LAMBDA1, LAMBDA2, LAMBDA3, LAMBDA4, LAMBDA
      INTEGER L1,L2

      DOUBLE COMPLEX DOTCC, DOTRC,contract_Tjj
      DOUBLE PRECISION DOTRR
      EXTERNAL DOTCC, DOTRC, DOTRR, contract_Tjj
      double precision sw2, sw
#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "gmunu.inc"
#include "transm.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"
      sw2 = sin2w
      sw = sqrt(sin2w)
      OFFeps(6,-1:1,4) = cZero
      K(0:3,1:4) = rZero
      M(-1:1,-1:1,-1:1,-1:1) = cZero
      PP(1:4,1:4) = rZero
      O(-1:1,1:4,-1:1,1:4) = cZero
      OP(1:4,-1:1,1:4) = cZero  
      
      OFFeps(1:6,-1:1,1) = h1(1:6,-1:1)
      OFFeps(1:6,-1:1,2) = h2(1:6,-1:1)
      OFFeps(1:6,-1:1,3) = wm1(1:6,-1:1)
      OFFeps(1:6,-1:1,4) = wp1(1:6,-1:1) 

      do i = 1,4
        K(0,i) = dble(OFFeps(5,-1,i))*(-rONE) !K1(0:3)
        K(1,i) = dble(OFFeps(6,-1,i))*(-rONE) !K2(0:3)
        K(2,i) = dimag(OFFeps(6,-1,i))*(-rONE)!K3(0:3)
        K(3,i) = dimag(OFFeps(5,-1,i))*(-rONE)!K4(0:3)
      enddo
      



      do i = 1,4
        do j = 1,4
          pp(i,j) = dotrr(k(0,i),k(0,j))
        enddo
      enddo

      do L1 = -1,1
        do L2 = -1,1 
          do i = 1,4
            do j = 1,4
              O(L1,i,L2,j) = dotcc(OFFeps(1,L1,i),OFFeps(1,L2,j))
            enddo
          enddo
        enddo
      enddo

      do LAMBDA = -1,1
        do i = 1,4 
          do j = 1,4
            OP(j,lambda,i) = dotrc(k(0,j),OFFeps(1,lambda,i))
          enddo
        enddo
      enddo


!!!! ===============================================================      
!!!! COUPLINGS CALCULATION
!!!! ===============================================================



         if( fs0.ne.0d0 ) then  !fs0 coupling
!            do lambda1 = -1,1
!             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(0,0,lambda3,lambda4) = M(0,0,lambda3,lambda4) -(fs0*(WMASS**2)*(OP(1,lambda4,4)*OP(2,lambda3,3)+OP(1,lambda3,3)*OP(2,lambda4,4))*1.0*1.0)

               enddo                      
              enddo
!             enddo 
!            enddo
         endif

         if( fs1.ne.0d0 ) then  !fs1 coupling
!            do lambda1 = -1,1
!             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(0,0,lambda3,lambda4) = M(0,0,lambda3,lambda4) -2*fs1*(WMASS**2)*1.0*1.0*pp(1,2)*O(lambda3,3,lambda4,4)

               enddo                      
              enddo
!             enddo 
!            enddo
         endif

         if( fs2.ne.0d0 ) then  !fs2 coupling
!            do lambda1 = -1,1
!             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(0,0,lambda3,lambda4) = M(0,0,lambda3,lambda4) -(fs2*(WMASS**2)*(OP(1,lambda4,4)*OP(2,lambda3,3)+OP(1,lambda3,3)*OP(2,lambda4,4))*1.0*1.0)

               enddo                      
              enddo
!             enddo 
!            enddo
         endif                  



         if( fm0.ne.0d0 ) then  !fm0 coupling
!            do lambda1 = -1,1
!             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(0,0,lambda3,lambda4) = M(0,0,lambda3,lambda4) +2*fm0*(gw**2)*1.0*1.0*pp(1,2)*(OP(3,lambda4,4)*OP(4,lambda3,3)-pp(3,4)*O(lambda3,3,lambda4,4))

               enddo                      
              enddo
!             enddo 
!            enddo
         endif


         if( fm1.ne.0d0 ) then  !fm1 coupling
!            do lambda1 = -1,1
!             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(0,0,lambda3,lambda4) = M(0,0,lambda3,lambda4) -(fm1*(gw**2)*1.0*1.0*(OP(1,lambda4,4)*OP(4,lambda3,3)*pp(2,3)+OP(1,lambda3,3)*OP(3,lambda4,4)*pp(2,4)+OP(2,lambda4,4)*(OP(4,lambda3,3)*pp(1,3)-OP(1,lambda3,3)*pp(3,4))+OP(2,lambda3,3)*(OP(3,lambda4,4)*pp(1,4)-OP(1,lambda4,4)*pp(3,4))-pp(1,4)*pp(2,3)*O(lambda3,3,lambda4,4)-pp(1,3)*pp(2,4)*O(lambda3,3,lambda4,4)))/2

               enddo                      
              enddo
!             enddo 
!            enddo
         endif




         if( fm6.ne.0d0 ) then  !fm6 coupling
!            do lambda1 = -1,1
!             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(0,0,lambda3,lambda4) = M(0,0,lambda3,lambda4) +fm6*(gw**2)*1.0*1.0*pp(1,2)*(OP(3,lambda4,4)*OP(4,lambda3,3)-pp(3,4)*O(lambda3,3,lambda4,4))

               enddo                      
              enddo
!             enddo 
!            enddo
         endif

         if( fm7.ne.0d0 ) then  !fm7 coupling
!            do lambda1 = -1,1
!             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(0,0,lambda3,lambda4) = M(0,0,lambda3,lambda4) +(fm7*(gw**2)*1.0*1.0*(OP(1,lambda4,4)*OP(4,lambda3,3)*pp(2,3)+OP(1,lambda3,3)*OP(3,lambda4,4)*pp(2,4)+OP(2,lambda4,4)*(OP(4,lambda3,3)*pp(1,3)-OP(1,lambda3,3)*pp(3,4))+OP(2,lambda3,3)*(OP(3,lambda4,4)*pp(1,4)-OP(1,lambda4,4)*pp(3,4))-pp(1,4)*pp(2,3)*O(lambda3,3,lambda4,4)-pp(1,3)*pp(2,4)*O(lambda3,3,lambda4,4)))/4

               enddo                      
              enddo
!             enddo 
!            enddo
         endif         

      end






ccccccccccccccc
c hzzz-vertex c
ccccccccccccccc

      subroutine hzzz_ONanom4(h1,z1,z2,z3,M)
c
c This subroutine calculates the vvvv-vertex with anomalous couplings
c and the on-shell polarization vectors for vector bosons
c
c IN:
c     polv1, polv2      : the incoming polvector wz(bosons in the CM frame)
c     polv3, polv4      : the outgoing polvector wz(bosons in the CM frame)
c OUT:
c     M(lambda1,lambda2) 
c
      IMPLICIT NONE

      INTEGER NGRAPHS
      PARAMETER( NGRAPHS=20 )
      DOUBLE PRECISION rZERO, rONE, rHalf
      PARAMETER( rZERO=0.0D0, rONE=1.0D0, rHALF = 0.5D0 )
      DOUBLE COMPLEX cIMAG, cZERO, cONE
      PARAMETER( cIMAG=(0.0D0,1.0D0), cZERO=(0.0D0,0.0D0), cONE=(1.0D0,0.0D0))
      DOUBLE COMPLEX h1(6,-1:1), z1(6,-1:1), z2(6,-1:1), z3(6,-1:1)      
      DOUBLE COMPLEX M(-1:1,-1:1,-1:1,-1:1)
      
      DOUBLE COMPLEX OFFeps(6,-1:1,NGRAPHS), OFFepsC(6,-1:1,NGRAPHS)
      DOUBLE PRECISION K(0:3,4)      
      DOUBLE COMPLEX O(-1:1,4,-1:1,4), OP(4,-1:1,4)
      DOUBLE COMPLEXPP(1:4,1:4)
    
      INTEGER I,J,N,MU,NU
      INTEGER LAMBDA1, LAMBDA2, LAMBDA3, LAMBDA4, LAMBDA
      INTEGER L1,L2

      DOUBLE COMPLEX DOTCC, DOTRC,contract_Tjj
      DOUBLE PRECISION DOTRR
      EXTERNAL DOTCC, DOTRC, DOTRR, contract_Tjj
      double precision sw2, sw
#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "gmunu.inc"
#include "transm.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"
      sw2 = sin2w
      sw = sqrt(sin2w)
      OFFeps(6,-1:1,4) = cZero
      K(0:3,1:4) = rZero
      M(-1:1,-1:1,-1:1,-1:1) = cZero
      PP(1:4,1:4) = rZero
      O(-1:1,1:4,-1:1,1:4) = cZero
      OP(1:4,-1:1,1:4) = cZero  
      
      OFFeps(1:6,-1:1,1) = h1(1:6,-1:1)
      OFFeps(1:6,-1:1,2) = z1(1:6,-1:1)
      OFFeps(1:6,-1:1,3) = z2(1:6,-1:1)
      OFFeps(1:6,-1:1,4) = z3(1:6,-1:1) 

      do i = 1,4
        K(0,i) = dble(OFFeps(5,-1,i))*(-rONE) !K1(0:3)
        K(1,i) = dble(OFFeps(6,-1,i))*(-rONE) !K2(0:3)
        K(2,i) = dimag(OFFeps(6,-1,i))*(-rONE)!K3(0:3)
        K(3,i) = dimag(OFFeps(5,-1,i))*(-rONE)!K4(0:3)
      enddo
      



      do i = 1,4
        do j = 1,4
          pp(i,j) = dotrr(k(0,i),k(0,j))
        enddo
      enddo

      do L1 = -1,1
        do L2 = -1,1 
          do i = 1,4
            do j = 1,4
              O(L1,i,L2,j) = dotcc(OFFeps(1,L1,i),OFFeps(1,L2,j))
            enddo
          enddo
        enddo
      enddo

      do LAMBDA = -1,1
        do i = 1,4 
          do j = 1,4
            OP(j,lambda,i) = dotrc(k(0,j),OFFeps(1,lambda,i))
          enddo
        enddo
      enddo


!!!! ===============================================================      
!!!! COUPLINGS CALCULATION
!!!! ===============================================================





      end

ccccccccccccccc
c hzwmwp-vertex c
ccccccccccccccc

      subroutine hzwmwp_ONanom4(h,z,wm1,wp1,M)
c
c This subroutine calculates the vvvv-vertex with anomalous couplings
c and the on-shell polarization vectors for vector bosons
c
c IN:
c     polv1, polv2      : the incoming polvector wz(bosons in the CM frame)
c     polv3, polv4      : the outgoing polvector wz(bosons in the CM frame)
c OUT:
c     M(lambda1,lambda2) 
c
      IMPLICIT NONE

      INTEGER NGRAPHS
      PARAMETER( NGRAPHS=20 )
      DOUBLE PRECISION rZERO, rONE, rHalf
      PARAMETER( rZERO=0.0D0, rONE=1.0D0, rHALF = 0.5D0 )
      DOUBLE COMPLEX cIMAG, cZERO, cONE
      PARAMETER( cIMAG=(0.0D0,1.0D0), cZERO=(0.0D0,0.0D0), cONE=(1.0D0,0.0D0))
      DOUBLE COMPLEX h(6,-1:1), z(6,-1:1), wm1(6,-1:1), wp1(6,-1:1)      
      DOUBLE COMPLEX M(-1:1,-1:1,-1:1,-1:1)
      
      DOUBLE COMPLEX OFFeps(6,-1:1,NGRAPHS), OFFepsC(6,-1:1,NGRAPHS)
      DOUBLE PRECISION K(0:3,4)      
      DOUBLE COMPLEX O(-1:1,4,-1:1,4), OP(4,-1:1,4)
      DOUBLE COMPLEXPP(1:4,1:4)
    
      INTEGER I,J,N,MU,NU
      INTEGER LAMBDA1, LAMBDA2, LAMBDA3, LAMBDA4, LAMBDA
      INTEGER L1,L2

      DOUBLE COMPLEX DOTCC, DOTRC,contract_Tjj
      DOUBLE PRECISION DOTRR
      EXTERNAL DOTCC, DOTRC, DOTRR, contract_Tjj
      double precision sw2, sw
#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "gmunu.inc"
#include "transm.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"
      sw2 = sin2w
      sw = sqrt(sin2w)
      OFFeps(6,-1:1,4) = cZero
      K(0:3,1:4) = rZero
      M(-1:1,-1:1,-1:1,-1:1) = cZero
      PP(1:4,1:4) = rZero
      O(-1:1,1:4,-1:1,1:4) = cZero
      OP(1:4,-1:1,1:4) = cZero  
      
      OFFeps(1:6,-1:1,1) = h(1:6,-1:1)
      OFFeps(1:6,-1:1,2) = z(1:6,-1:1)
      OFFeps(1:6,-1:1,3) = wm1(1:6,-1:1)
      OFFeps(1:6,-1:1,4) = wp1(1:6,-1:1) 

      do i = 1,4
        K(0,i) = dble(OFFeps(5,-1,i))*(-rONE) !K1(0:3)
        K(1,i) = dble(OFFeps(6,-1,i))*(-rONE) !K2(0:3)
        K(2,i) = dimag(OFFeps(6,-1,i))*(-rONE)!K3(0:3)
        K(3,i) = dimag(OFFeps(5,-1,i))*(-rONE)!K4(0:3)
      enddo
      



      do i = 1,4
        do j = 1,4
          pp(i,j) = dotrr(k(0,i),k(0,j))
        enddo
      enddo

      do L1 = -1,1
        do L2 = -1,1 
          do i = 1,4
            do j = 1,4
              O(L1,i,L2,j) = dotcc(OFFeps(1,L1,i),OFFeps(1,L2,j))
            enddo
          enddo
        enddo
      enddo

      do LAMBDA = -1,1
        do i = 1,4 
          do j = 1,4
            OP(j,lambda,i) = dotrc(k(0,j),OFFeps(1,lambda,i))
          enddo
        enddo
      enddo


!!!! ===============================================================      
!!!! COUPLINGS CALCULATION
!!!! ===============================================================



         if( fs0.ne.0d0 ) then  !fs0 coupling
!            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(0,lambda2,lambda3,lambda4) = M(0,lambda2,lambda3,lambda4) +fs0*(WMASS**2)*ZMASS*1.0*(OP(1,lambda4,4)*O(lambda2,2,lambda3,3)-OP(1,lambda3,3)*O(lambda2,2,lambda4,4))

               enddo                      
              enddo
             enddo 
!            enddo
         endif

         if( fs2.ne.0d0 ) then  !fs2 coupling
!            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(0,lambda2,lambda3,lambda4) = M(0,lambda2,lambda3,lambda4) -fs2*(WMASS**2)*ZMASS*1.0*(OP(1,lambda4,4)*O(lambda2,2,lambda3,3)-OP(1,lambda3,3)*O(lambda2,2,lambda4,4))

               enddo                      
              enddo
             enddo 
!            enddo
         endif                  



         if( fm4.ne.0d0 ) then  !fm4 coupling
!            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(0,lambda2,lambda3,lambda4) = M(0,lambda2,lambda3,lambda4) +(fm4*(gw**2)*(sw**2)*ZMASS*1.0*(OP(1,lambda4,4)*(OP(2,lambda3,3)*OP(3,lambda2,2)-pp(3,2)*O(lambda2,2,lambda3,3))+OP(1,lambda3,3)*(-(OP(2,lambda4,4)*OP(4,lambda2,2))+pp(4,2)*O(lambda2,2,lambda4,4))))/2

               enddo                      
              enddo
             enddo 
!            enddo
         endif

         if( fm5.ne.0d0 ) then  !fm5 coupling
!            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(0,lambda2,lambda3,lambda4) = M(0,lambda2,lambda3,lambda4) +(fm5*(gw**2)*(sw**2)*ZMASS*1.0*(-(OP(2,lambda4,4)*pp(1,3)*O(lambda2,2,lambda3,3))+OP(4,lambda3,3)*pp(1,2)*O(lambda2,2,lambda4,4)+OP(2,lambda3,3)*pp(1,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(OP(2,lambda4,4)*OP(3,lambda2,2)-pp(3,2)*O(lambda2,2,lambda4,4))-OP(4,lambda2,2)*pp(1,2)*O(lambda3,3,lambda4,4)+OP(1,lambda2,2)*(-(OP(2,lambda4,4)*OP(4,lambda3,3))+pp(4,2)*O(lambda3,3,lambda4,4))))/4

               enddo                      
              enddo
             enddo 
!            enddo
         endif

         if( fm5hc.ne.0d0 ) then  !fm5hc coupling
!            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(0,lambda2,lambda3,lambda4) = M(0,lambda2,lambda3,lambda4) +(fm5hc*(gw**2)*(sw**2)*ZMASS*1.0*(-(OP(3,lambda4,4)*pp(1,2)*O(lambda2,2,lambda3,3))-OP(2,lambda4,4)*pp(1,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*(-(OP(2,lambda3,3)*OP(4,lambda2,2))+pp(4,2)*O(lambda2,2,lambda3,3))+OP(2,lambda3,3)*pp(1,4)*O(lambda2,2,lambda4,4)+OP(3,lambda2,2)*pp(1,2)*O(lambda3,3,lambda4,4)+OP(1,lambda2,2)*(OP(2,lambda3,3)*OP(3,lambda4,4)-pp(3,2)*O(lambda3,3,lambda4,4))))/4

               enddo                      
              enddo
             enddo 
!            enddo
         endif


         if( fm7.ne.0d0 ) then  !fm7 coupling
!            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(0,lambda2,lambda3,lambda4) = M(0,lambda2,lambda3,lambda4) +(fm7*(gw**2)*ZMASS*1.0*(-(OP(1,lambda2,2)*OP(2,lambda3,3)*OP(3,lambda4,4))+OP(1,lambda2,2)*OP(2,lambda4,4)*OP(4,lambda3,3)+OP(3,lambda4,4)*pp(1,2)*O(lambda2,2,lambda3,3)-OP(2,lambda4,4)*pp(1,3)*O(lambda2,2,lambda3,3)-OP(2,lambda4,4)*pp(1,4)*O(lambda2,2,lambda3,3)+OP(3,lambda4,4)*pp(1,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*((-1+SW2)*OP(2,lambda3,3)*OP(4,lambda2,2)+OP(3,lambda2,2)*OP(4,lambda3,3)-(pp(3,4)+(-1+SW2)*pp(4,2))*O(lambda2,2,lambda3,3))-OP(4,lambda3,3)*pp(1,2)*O(lambda2,2,lambda4,4)+OP(2,lambda3,3)*pp(1,3)*O(lambda2,2,lambda4,4)-OP(4,lambda3,3)*pp(1,3)*O(lambda2,2,lambda4,4)+OP(2,lambda3,3)*pp(1,4)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(-((-1+SW2)*OP(2,lambda4,4)*OP(3,lambda2,2))-OP(3,lambda4,4)*OP(4,lambda2,2)+((-1+SW2)*pp(3,2)+pp(3,4))*O(lambda2,2,lambda4,4))-OP(3,lambda2,2)*pp(1,2)*O(lambda3,3,lambda4,4)+OP(4,lambda2,2)*pp(1,2)*O(lambda3,3,lambda4,4)+OP(4,lambda2,2)*pp(1,3)*O(lambda3,3,lambda4,4)-OP(3,lambda2,2)*pp(1,4)*O(lambda3,3,lambda4,4)+OP(1,lambda2,2)*pp(3,2)*O(lambda3,3,lambda4,4)-OP(1,lambda2,2)*pp(4,2)*O(lambda3,3,lambda4,4)+SW2*(-(OP(3,lambda4,4)*pp(1,2)*O(lambda2,2,lambda3,3))+OP(2,lambda4,4)*(pp(1,3)+pp(1,4))*O(lambda2,2,lambda3,3)+OP(4,lambda3,3)*pp(1,2)*O(lambda2,2,lambda4,4)-OP(2,lambda3,3)*pp(1,3)*O(lambda2,2,lambda4,4)-OP(2,lambda3,3)*pp(1,4)*O(lambda2,2,lambda4,4)+OP(3,lambda2,2)*pp(1,2)*O(lambda3,3,lambda4,4)-OP(4,lambda2,2)*pp(1,2)*O(lambda3,3,lambda4,4)+OP(1,lambda2,2)*(OP(2,lambda3,3)*OP(3,lambda4,4)-OP(2,lambda4,4)*OP(4,lambda3,3)+(-pp(3,2)+pp(4,2))*O(lambda3,3,lambda4,4)))))/4

               enddo                      
              enddo
             enddo 
!            enddo
         endif         

      end




ccccccccccccccc
c zzzz-vertex c
ccccccccccccccc

      subroutine zzzz_ONanom4(z1,z2,z3,z4,M)
c
c This subroutine calculates the vvvv-vertex with anomalous couplings
c and the on-shell polarization vectors for vector bosons
c
c IN:
c     polv1, polv2      : the incoming polvector wz(bosons in the CM frame)
c     polv3, polv4      : the outgoing polvector wz(bosons in the CM frame)
c OUT:
c     M(lambda1,lambda2) 
c
      IMPLICIT NONE

      INTEGER NGRAPHS
      PARAMETER( NGRAPHS=20 )
      DOUBLE PRECISION rZERO, rONE, rHalf
      PARAMETER( rZERO=0.0D0, rONE=1.0D0, rHALF = 0.5D0 )
      DOUBLE COMPLEX cIMAG, cZERO, cONE
      PARAMETER( cIMAG=(0.0D0,1.0D0), cZERO=(0.0D0,0.0D0), cONE=(1.0D0,0.0D0))
      DOUBLE COMPLEX z1(6,-1:1), z2(6,-1:1), z3(6,-1:1), z4(6,-1:1)      
      DOUBLE COMPLEX M(-1:1,-1:1,-1:1,-1:1)
      
      DOUBLE COMPLEX OFFeps(6,-1:1,NGRAPHS), OFFepsC(6,-1:1,NGRAPHS)
      DOUBLE PRECISION K(0:3,4)      
      DOUBLE COMPLEX O(-1:1,4,-1:1,4), OP(4,-1:1,4)
      DOUBLE COMPLEXPP(1:4,1:4)
    
      INTEGER I,J,N,MU,NU
      INTEGER LAMBDA1, LAMBDA2, LAMBDA3, LAMBDA4, LAMBDA
      INTEGER L1,L2

      DOUBLE COMPLEX DOTCC, DOTRC,contract_Tjj
      DOUBLE PRECISION DOTRR
      EXTERNAL DOTCC, DOTRC, DOTRR, contract_Tjj
      double precision sw2, sw
#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "gmunu.inc"
#include "transm.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"
      sw2 = sin2w
      sw = sqrt(sin2w)
      OFFeps(6,-1:1,4) = cZero
      K(0:3,1:4) = rZero
      M(-1:1,-1:1,-1:1,-1:1) = cZero
      PP(1:4,1:4) = rZero
      O(-1:1,1:4,-1:1,1:4) = cZero
      OP(1:4,-1:1,1:4) = cZero  
      
      OFFeps(1:6,-1:1,1) = z1(1:6,-1:1)
      OFFeps(1:6,-1:1,2) = z2(1:6,-1:1)
      OFFeps(1:6,-1:1,3) = z3(1:6,-1:1)
      OFFeps(1:6,-1:1,4) = z4(1:6,-1:1) 

      do i = 1,4
        K(0,i) = dble(OFFeps(5,-1,i))*(-rONE) !K1(0:3)
        K(1,i) = dble(OFFeps(6,-1,i))*(-rONE) !K2(0:3)
        K(2,i) = dimag(OFFeps(6,-1,i))*(-rONE)!K3(0:3)
        K(3,i) = dimag(OFFeps(5,-1,i))*(-rONE)!K4(0:3)
      enddo
      



      do i = 1,4
        do j = 1,4
          pp(i,j) = dotrr(k(0,i),k(0,j))
        enddo
      enddo

      do L1 = -1,1
        do L2 = -1,1 
          do i = 1,4
            do j = 1,4
              O(L1,i,L2,j) = dotcc(OFFeps(1,L1,i),OFFeps(1,L2,j))
            enddo
          enddo
        enddo
      enddo

      do LAMBDA = -1,1
        do i = 1,4 
          do j = 1,4
            OP(j,lambda,i) = dotrc(k(0,j),OFFeps(1,lambda,i))
          enddo
        enddo
      enddo


!!!! ===============================================================      
!!!! COUPLINGS CALCULATION
!!!! ===============================================================
!

         if( fs0.ne.0d0 ) then  !fs0 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +2*fs0*(ZMASS**4)*(O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4))

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fs1.ne.0d0 ) then  !fs1 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +2*fs1*(ZMASS**4)*(O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4))

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fs2.ne.0d0 ) then  !fs2 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +2*fs2*(ZMASS**4)*(O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4))

               enddo                      
              enddo
             enddo 
            enddo
         endif                  

         if( ft0.ne.0d0 ) then  !ft0 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +8*ft0*(gw**4)*(-1+SW2)**2*(OP(1,lambda2,2)*OP(2,lambda1,1)*OP(3,lambda4,4)*OP(4,lambda3,3)-OP(3,lambda4,4)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda2,2)-OP(2,lambda4,4)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda3,3)-OP(2,lambda3,3)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*OP(4,lambda1,1)*(OP(2,lambda3,3)*OP(3,lambda2,2)-pp(2,3)*O(lambda2,2,lambda3,3))+pp(1,3)*pp(2,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*OP(3,lambda1,1)*(OP(2,lambda4,4)*OP(4,lambda2,2)-pp(2,4)*O(lambda2,2,lambda4,4))-OP(1,lambda2,2)*OP(2,lambda1,1)*pp(3,4)*O(lambda3,3,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4))

               enddo                      
              enddo
             enddo 
            enddo
         endif


         if( ft1.ne.0d0 ) then  !ft1 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +8*ft1*(gw**4)*(-1+SW2)**2*(OP(1,lambda2,2)*OP(2,lambda1,1)*OP(3,lambda4,4)*OP(4,lambda3,3)-OP(3,lambda4,4)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda2,2)-OP(2,lambda4,4)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda3,3)-OP(2,lambda3,3)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*OP(4,lambda1,1)*(OP(2,lambda3,3)*OP(3,lambda2,2)-pp(2,3)*O(lambda2,2,lambda3,3))+pp(1,3)*pp(2,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*OP(3,lambda1,1)*(OP(2,lambda4,4)*OP(4,lambda2,2)-pp(2,4)*O(lambda2,2,lambda4,4))-OP(1,lambda2,2)*OP(2,lambda1,1)*pp(3,4)*O(lambda3,3,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4))

               enddo                      
              enddo
             enddo 
            enddo
         endif


         if( ft2.ne.0d0 ) then  !ft2 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +2*ft2*(gw**4)*(-1+SW2)**2*(OP(1,lambda4,4)*OP(2,lambda3,3)*OP(3,lambda1,1)*OP(4,lambda2,2)+OP(1,lambda4,4)*OP(2,lambda1,1)*OP(3,lambda2,2)*OP(4,lambda3,3)-OP(2,lambda4,4)*OP(4,lambda3,3)*pp(1,3)*O(lambda1,1,lambda2,2)-OP(2,lambda3,3)*OP(3,lambda4,4)*pp(1,4)*O(lambda1,1,lambda2,2)-OP(1,lambda4,4)*OP(4,lambda3,3)*pp(2,3)*O(lambda1,1,lambda2,2)+OP(1,lambda4,4)*OP(2,lambda3,3)*pp(3,4)*O(lambda1,1,lambda2,2)-OP(3,lambda4,4)*OP(4,lambda2,2)*pp(1,2)*O(lambda1,1,lambda3,3)-OP(2,lambda4,4)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda3,3)-OP(1,lambda4,4)*OP(4,lambda2,2)*pp(2,3)*O(lambda1,1,lambda3,3)+OP(1,lambda4,4)*OP(3,lambda2,2)*pp(2,4)*O(lambda1,1,lambda3,3)-OP(3,lambda2,2)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda4,4)-OP(2,lambda3,3)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda4,4)-OP(3,lambda4,4)*OP(4,lambda1,1)*pp(1,2)*O(lambda2,2,lambda3,3)-OP(2,lambda4,4)*OP(4,lambda1,1)*pp(1,3)*O(lambda2,2,lambda3,3)+OP(2,lambda4,4)*OP(3,lambda1,1)*pp(1,4)*O(lambda2,2,lambda3,3)+OP(2,lambda1,1)*OP(3,lambda4,4)*pp(1,4)*O(lambda2,2,lambda3,3)-OP(1,lambda4,4)*OP(3,lambda1,1)*pp(2,4)*O(lambda2,2,lambda3,3)-OP(1,lambda4,4)*OP(2,lambda1,1)*pp(3,4)*O(lambda2,2,lambda3,3)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-OP(3,lambda1,1)*OP(4,lambda3,3)*pp(1,2)*O(lambda2,2,lambda4,4)+OP(2,lambda3,3)*OP(4,lambda1,1)*pp(1,3)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*OP(4,lambda3,3)*pp(1,3)*O(lambda2,2,lambda4,4)-OP(2,lambda3,3)*OP(3,lambda1,1)*pp(1,4)*O(lambda2,2,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(-(OP(3,lambda4,4)*pp(2,4)*O(lambda1,1,lambda2,2))+OP(2,lambda4,4)*(OP(3,lambda2,2)*OP(4,lambda1,1)+pp(3,4)*O(lambda1,1,lambda2,2))+OP(4,lambda2,2)*pp(2,3)*O(lambda1,1,lambda4,4)-OP(3,lambda2,2)*pp(2,4)*O(lambda1,1,lambda4,4)-OP(4,lambda1,1)*pp(2,3)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*(OP(3,lambda4,4)*OP(4,lambda2,2)-pp(3,4)*O(lambda2,2,lambda4,4)))+OP(3,lambda2,2)*OP(4,lambda1,1)*pp(1,2)*O(lambda3,3,lambda4,4)+OP(3,lambda1,1)*OP(4,lambda2,2)*pp(1,2)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(4,lambda2,2)*pp(1,3)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(3,lambda2,2)*pp(1,4)*O(lambda3,3,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+OP(1,lambda2,2)*(OP(3,lambda4,4)*pp(2,4)*O(lambda1,1,lambda3,3)+OP(2,lambda4,4)*(OP(3,lambda1,1)*OP(4,lambda3,3)-pp(3,4)*O(lambda1,1,lambda3,3))+OP(4,lambda3,3)*pp(2,3)*O(lambda1,1,lambda4,4)+OP(2,lambda3,3)*(OP(3,lambda4,4)*OP(4,lambda1,1)-pp(3,4)*O(lambda1,1,lambda4,4))-OP(4,lambda1,1)*pp(2,3)*O(lambda3,3,lambda4,4)-OP(3,lambda1,1)*pp(2,4)*O(lambda3,3,lambda4,4)))

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( ft3.ne.0d0 ) then  !ft3 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +2*ft3*(SW2-1)**2*(gw**4)*(-pp(1,2)*OP(4,lambda2,2)*OP(3,lambda4,4)*O(lambda1,1,lambda3,3)-pp(2,3)*OP(4,lambda2,2)*OP(1,lambda4,4)*O(lambda1,1,lambda3,3)-pp(1,3)*OP(4,lambda2,2)*OP(2,lambda3,3)*O(lambda1,1,lambda4,4)+pp(1,2)*OP(3,lambda1,1)*OP(4,lambda2,2)*O(lambda3,3,lambda4,4)-pp(1,3)*OP(2,lambda1,1)*OP(4,lambda2,2)*O(lambda3,3,lambda4,4)-pp(1,3)*OP(4,lambda3,3)*OP(2,lambda4,4)*O(lambda1,1,lambda2,2)-pp(1,4)*OP(2,lambda3,3)*OP(3,lambda4,4)*O(lambda1,1,lambda2,2)-pp(2,3)*OP(4,lambda3,3)*OP(1,lambda4,4)*O(lambda1,1,lambda2,2)+pp(3,4)*OP(2,lambda3,3)*OP(1,lambda4,4)*O(lambda1,1,lambda2,2)-pp(1,4)*OP(3,lambda2,2)*OP(2,lambda4,4)*O(lambda1,1,lambda3,3)+pp(2,4)*OP(3,lambda2,2)*OP(1,lambda4,4)*O(lambda1,1,lambda3,3)-pp(1,2)*OP(3,lambda2,2)*OP(4,lambda3,3)*O(lambda1,1,lambda4,4)-pp(1,2)*OP(4,lambda1,1)*OP(3,lambda4,4)*O(lambda2,2,lambda3,3)-pp(1,3)*OP(4,lambda1,1)*OP(2,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,4)*OP(3,lambda1,1)*OP(2,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,4)*OP(2,lambda1,1)*OP(3,lambda4,4)*O(lambda2,2,lambda3,3)-pp(2,4)*OP(3,lambda1,1)*OP(1,lambda4,4)*O(lambda2,2,lambda3,3)-pp(3,4)*OP(2,lambda1,1)*OP(1,lambda4,4)*O(lambda2,2,lambda3,3)-pp(1,2)*OP(3,lambda1,1)*OP(4,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,3)*OP(4,lambda1,1)*OP(2,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,3)*OP(2,lambda1,1)*OP(4,lambda3,3)*O(lambda2,2,lambda4,4)-pp(1,4)*OP(3,lambda1,1)*OP(2,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(pp(2,4)*(-OP(3,lambda4,4))*O(lambda1,1,lambda2,2)+OP(2,lambda4,4)*(OP(4,lambda1,1)*OP(3,lambda2,2)+pp(3,4)*O(lambda1,1,lambda2,2))+pp(2,3)*OP(4,lambda2,2)*O(lambda1,1,lambda4,4)-pp(2,4)*OP(3,lambda2,2)*O(lambda1,1,lambda4,4)-pp(2,3)*OP(4,lambda1,1)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*(OP(4,lambda2,2)*OP(3,lambda4,4)-pp(3,4)*O(lambda2,2,lambda4,4)))+pp(1,2)*OP(4,lambda1,1)*OP(3,lambda2,2)*O(lambda3,3,lambda4,4)-pp(1,4)*OP(2,lambda1,1)*OP(3,lambda2,2)*O(lambda3,3,lambda4,4)+OP(1,lambda2,2)*(pp(2,4)*OP(3,lambda4,4)*O(lambda1,1,lambda3,3)+OP(2,lambda4,4)*(OP(3,lambda1,1)*OP(4,lambda3,3)-pp(3,4)*O(lambda1,1,lambda3,3))+pp(2,3)*OP(4,lambda3,3)*O(lambda1,1,lambda4,4)+OP(2,lambda3,3)*(OP(4,lambda1,1)*OP(3,lambda4,4)-pp(3,4)*O(lambda1,1,lambda4,4))-pp(2,3)*OP(4,lambda1,1)*O(lambda3,3,lambda4,4)-pp(2,4)*OP(3,lambda1,1)*O(lambda3,3,lambda4,4))+OP(3,lambda1,1)*OP(4,lambda2,2)*OP(2,lambda3,3)*OP(1,lambda4,4)+OP(2,lambda1,1)*OP(3,lambda2,2)*OP(4,lambda3,3)*OP(1,lambda4,4)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4))
               enddo                      
              enddo
             enddo 
            enddo
         endif


         if( ft4.ne.0d0 ) then  !ft4 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +gw**4*ft4*sw**4*(OP(1, lambda4, 4)*OP(2, lambda3, 3)*OP(3, lambda1, 1)*OP(4, lambda2, 2) + OP(1, lambda4, 4)*OP(2, lambda1, 1)*OP(3, lambda2, 2)*OP(4, lambda3, 3) + pp(3, 4)*OP(1, lambda4, 4)*OP(2, lambda3, 3)*O(lambda1, 1, lambda2, 2) - pp(1, 4)*OP(2, lambda3, 3)*OP(3, lambda4, 4)*O(lambda1, 1, lambda2, 2) - pp(2, 3)*OP(1, lambda4, 4)*OP(4, lambda3, 3)*O(lambda1, 1, lambda2, 2) - pp(1, 3)*OP(2, lambda4, 4)*OP(4, lambda3, 3)*O(lambda1, 1, lambda2, 2) + pp(2, 4)*OP(1, lambda4, 4)*OP(3, lambda2, 2)*O(lambda1, 1, lambda3, 3) - pp(1, 4)*OP(2, lambda4, 4)*OP(3, lambda2, 2)*O(lambda1, 1, lambda3, 3) - pp(2, 3)*OP(1, lambda4, 4)*OP(4, lambda2, 2)*O(lambda1, 1, lambda3, 3) - pp(1, 2)*OP(3, lambda4, 4)*OP(4, lambda2, 2)*O(lambda1, 1, lambda3, 3) - pp(1, 3)*OP(2, lambda3, 3)*OP(4, lambda2, 2)*O(lambda1, 1, lambda4, 4) - pp(1, 2)*OP(3, lambda2, 2)*OP(4, lambda3, 3)*O(lambda1, 1, lambda4, 4) - pp(3, 4)*OP(1, lambda4, 4)*OP(2, lambda1, 1)*O(lambda2, 2, lambda3, 3) - pp(2, 4)*OP(1, lambda4, 4)*OP(3, lambda1, 1)*O(lambda2, 2, lambda3, 3) + pp(1, 4)*OP(2, lambda4, 4)*OP(3, lambda1, 1)*O(lambda2, 2, lambda3, 3) + pp(1, 4)*OP(2, lambda1, 1)*OP(3, lambda4, 4)*O(lambda2, 2, lambda3, 3) - pp(1, 3)*OP(2, lambda4, 4)*OP(4, lambda1, 1)*O(lambda2, 2, lambda3, 3) - pp(1, 2)*OP(3, lambda4, 4)*OP(4, lambda1, 1)*O(lambda2, 2, lambda3, 3) + pp(1, 3)*pp(2, 4)*O(lambda1, 1, lambda4, 4)*O(lambda2, 2, lambda3, 3) + pp(1, 2)*pp(3, 4)*O(lambda1, 1, lambda4, 4)*O(lambda2, 2, lambda3, 3) - pp(1, 4)*OP(2, lambda3, 3)*OP(3, lambda1, 1)*O(lambda2, 2, lambda4, 4) + pp(1, 3)*OP(2, lambda3, 3)*OP(4, lambda1, 1)*O(lambda2, 2, lambda4, 4) + pp(1, 3)*OP(2, lambda1, 1)*OP(4, lambda3, 3)*O(lambda2, 2, lambda4, 4) - pp(1, 2)*OP(3, lambda1, 1)*OP(4, lambda3, 3)*O(lambda2, 2, lambda4, 4) + pp(1, 4)*pp(2, 3)*O(lambda1, 1, lambda3, 3)*O(lambda2, 2, lambda4, 4) + pp(1, 2)*pp(3, 4)*O(lambda1, 1, lambda3, 3)*O(lambda2, 2, lambda4, 4) + OP(1, lambda3, 3)*(-(pp(2, 4)*OP(3, lambda4, 4)*O(lambda1, 1, lambda2, 2)) + OP(2, lambda4, 4)*(OP(3, lambda2, 2)*OP(4, lambda1, 1) + pp(3, 4)*O(lambda1, 1, lambda2, 2)) -   pp(2, 4)*OP(3, lambda2, 2)*O(lambda1, 1, lambda4, 4) + pp(2, 3)*OP(4, lambda2, 2)*O(lambda1, 1, lambda4, 4) - pp(2, 3)*OP(4, lambda1, 1)*O(lambda2, 2, lambda4, 4) +   OP(2, lambda1, 1)*(OP(3, lambda4, 4)*OP(4, lambda2, 2) - pp(3, 4)*O(lambda2, 2, lambda4, 4))) - pp(1, 4)*OP(2, lambda1, 1)*OP(3, lambda2, 2)*O(lambda3, 3, lambda4, 4) + pp(1, 2)*OP(3, lambda2, 2)*OP(4, lambda1, 1)*O(lambda3, 3, lambda4, 4) - pp(1, 3)*OP(2, lambda1, 1)*OP(4, lambda2, 2)*O(lambda3, 3, lambda4, 4) + pp(1, 2)*OP(3, lambda1, 1)*OP(4, lambda2, 2)*O(lambda3, 3, lambda4, 4) + pp(1, 4)*pp(2, 3)*O(lambda1, 1, lambda2, 2)*O(lambda3, 3, lambda4, 4) + pp(1, 3)*pp(2, 4)*O(lambda1, 1, lambda2, 2)*O(lambda3, 3, lambda4, 4) + OP(1, lambda2, 2)*(pp(2, 4)*OP(3, lambda4, 4)*O(lambda1, 1, lambda3, 3) +   OP(2, lambda4, 4)*(OP(3, lambda1, 1)*OP(4, lambda3, 3) - pp(3, 4)*O(lambda1, 1, lambda3, 3)) + pp(2, 3)*OP(4, lambda3, 3)*O(lambda1, 1, lambda4, 4) +   OP(2, lambda3, 3)*(OP(3, lambda4, 4)*OP(4, lambda1, 1) - pp(3, 4)*O(lambda1, 1, lambda4, 4)) - pp(2, 4)*OP(3, lambda1, 1)*O(lambda3, 3, lambda4, 4) -   pp(2, 3)*OP(4, lambda1, 1)*O(lambda3, 3, lambda4, 4)))

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( ft5.ne.0d0 ) then  !ft5 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(-4*ft5*(gw**4)*(sw**4)*(-1+SW2)*(ZMASS**2)*(OP(1,lambda2,2)*OP(2,lambda1,1)*OP(3,lambda4,4)*OP(4,lambda3,3)-OP(3,lambda4,4)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda2,2)-OP(2,lambda4,4)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda3,3)-OP(2,lambda3,3)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*OP(4,lambda1,1)*(OP(2,lambda3,3)*OP(3,lambda2,2)-pp(2,3)*O(lambda2,2,lambda3,3))+pp(1,3)*pp(2,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*OP(3,lambda1,1)*(OP(2,lambda4,4)*OP(4,lambda2,2)-pp(2,4)*O(lambda2,2,lambda4,4))-OP(1,lambda2,2)*OP(2,lambda1,1)*pp(3,4)*O(lambda3,3,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)))/(WMASS**2)

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( ft6.ne.0d0 ) then  !ft6 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(-4*ft6*(gw**4)*(sw**4)*(-1+SW2)*(ZMASS**2)*(OP(1,lambda2,2)*OP(2,lambda1,1)*OP(3,lambda4,4)*OP(4,lambda3,3)-OP(3,lambda4,4)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda2,2)-OP(2,lambda4,4)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda3,3)-OP(2,lambda3,3)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*OP(4,lambda1,1)*(OP(2,lambda3,3)*OP(3,lambda2,2)-pp(2,3)*O(lambda2,2,lambda3,3))+pp(1,3)*pp(2,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*OP(3,lambda1,1)*(OP(2,lambda4,4)*OP(4,lambda2,2)-pp(2,4)*O(lambda2,2,lambda4,4))-OP(1,lambda2,2)*OP(2,lambda1,1)*pp(3,4)*O(lambda3,3,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)))/(WMASS**2)

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( ft7.ne.0d0 ) then  !ft7 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) -((ft7*(gw**4)*(sw**4)*(-1+SW2)*(ZMASS**2)*(OP(1,lambda4,4)*OP(2,lambda3,3)*OP(3,lambda1,1)*OP(4,lambda2,2)+OP(1,lambda4,4)*OP(2,lambda1,1)*OP(3,lambda2,2)*OP(4,lambda3,3)-OP(2,lambda4,4)*OP(4,lambda3,3)*pp(1,3)*O(lambda1,1,lambda2,2)-OP(2,lambda3,3)*OP(3,lambda4,4)*pp(1,4)*O(lambda1,1,lambda2,2)-OP(1,lambda4,4)*OP(4,lambda3,3)*pp(2,3)*O(lambda1,1,lambda2,2)+OP(1,lambda4,4)*OP(2,lambda3,3)*pp(3,4)*O(lambda1,1,lambda2,2)-OP(3,lambda4,4)*OP(4,lambda2,2)*pp(1,2)*O(lambda1,1,lambda3,3)-OP(2,lambda4,4)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda3,3)-OP(1,lambda4,4)*OP(4,lambda2,2)*pp(2,3)*O(lambda1,1,lambda3,3)+OP(1,lambda4,4)*OP(3,lambda2,2)*pp(2,4)*O(lambda1,1,lambda3,3)-OP(3,lambda2,2)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda4,4)-OP(2,lambda3,3)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda4,4)-OP(3,lambda4,4)*OP(4,lambda1,1)*pp(1,2)*O(lambda2,2,lambda3,3)-OP(2,lambda4,4)*OP(4,lambda1,1)*pp(1,3)*O(lambda2,2,lambda3,3)+OP(2,lambda4,4)*OP(3,lambda1,1)*pp(1,4)*O(lambda2,2,lambda3,3)+OP(2,lambda1,1)*OP(3,lambda4,4)*pp(1,4)*O(lambda2,2,lambda3,3)-OP(1,lambda4,4)*OP(3,lambda1,1)*pp(2,4)*O(lambda2,2,lambda3,3)-OP(1,lambda4,4)*OP(2,lambda1,1)*pp(3,4)*O(lambda2,2,lambda3,3)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-OP(3,lambda1,1)*OP(4,lambda3,3)*pp(1,2)*O(lambda2,2,lambda4,4)+OP(2,lambda3,3)*OP(4,lambda1,1)*pp(1,3)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*OP(4,lambda3,3)*pp(1,3)*O(lambda2,2,lambda4,4)-OP(2,lambda3,3)*OP(3,lambda1,1)*pp(1,4)*O(lambda2,2,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(-(OP(3,lambda4,4)*pp(2,4)*O(lambda1,1,lambda2,2))+OP(2,lambda4,4)*(OP(3,lambda2,2)*OP(4,lambda1,1)+pp(3,4)*O(lambda1,1,lambda2,2))+OP(4,lambda2,2)*pp(2,3)*O(lambda1,1,lambda4,4)-OP(3,lambda2,2)*pp(2,4)*O(lambda1,1,lambda4,4)-OP(4,lambda1,1)*pp(2,3)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*(OP(3,lambda4,4)*OP(4,lambda2,2)-pp(3,4)*O(lambda2,2,lambda4,4)))+OP(3,lambda2,2)*OP(4,lambda1,1)*pp(1,2)*O(lambda3,3,lambda4,4)+OP(3,lambda1,1)*OP(4,lambda2,2)*pp(1,2)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(4,lambda2,2)*pp(1,3)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(3,lambda2,2)*pp(1,4)*O(lambda3,3,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+OP(1,lambda2,2)*(OP(3,lambda4,4)*pp(2,4)*O(lambda1,1,lambda3,3)+OP(2,lambda4,4)*(OP(3,lambda1,1)*OP(4,lambda3,3)-pp(3,4)*O(lambda1,1,lambda3,3))+OP(4,lambda3,3)*pp(2,3)*O(lambda1,1,lambda4,4)+OP(2,lambda3,3)*(OP(3,lambda4,4)*OP(4,lambda1,1)-pp(3,4)*O(lambda1,1,lambda4,4))-OP(4,lambda1,1)*pp(2,3)*O(lambda3,3,lambda4,4)-OP(3,lambda1,1)*pp(2,4)*O(lambda3,3,lambda4,4))))/(WMASS**2))

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( ft8.ne.0d0 ) then  !ft8 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(2*ft8*(gw**4)*(sw**8)*(ZMASS**4)*(OP(1,lambda2,2)*OP(2,lambda1,1)*OP(3,lambda4,4)*OP(4,lambda3,3)-OP(3,lambda4,4)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda2,2)-OP(2,lambda4,4)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda3,3)-OP(2,lambda3,3)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*OP(4,lambda1,1)*(OP(2,lambda3,3)*OP(3,lambda2,2)-pp(2,3)*O(lambda2,2,lambda3,3))+pp(1,3)*pp(2,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*OP(3,lambda1,1)*(OP(2,lambda4,4)*OP(4,lambda2,2)-pp(2,4)*O(lambda2,2,lambda4,4))-OP(1,lambda2,2)*OP(2,lambda1,1)*pp(3,4)*O(lambda3,3,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)))/(WMASS**4)

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( ft9.ne.0d0 ) then  !ft9 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(ft9*(gw**4)*(sw**8)*(ZMASS**4)*(OP(1,lambda4,4)*OP(2,lambda3,3)*OP(3,lambda1,1)*OP(4,lambda2,2)+OP(1,lambda4,4)*OP(2,lambda1,1)*OP(3,lambda2,2)*OP(4,lambda3,3)-OP(2,lambda4,4)*OP(4,lambda3,3)*pp(1,3)*O(lambda1,1,lambda2,2)-OP(2,lambda3,3)*OP(3,lambda4,4)*pp(1,4)*O(lambda1,1,lambda2,2)-OP(1,lambda4,4)*OP(4,lambda3,3)*pp(2,3)*O(lambda1,1,lambda2,2)+OP(1,lambda4,4)*OP(2,lambda3,3)*pp(3,4)*O(lambda1,1,lambda2,2)-OP(3,lambda4,4)*OP(4,lambda2,2)*pp(1,2)*O(lambda1,1,lambda3,3)-OP(2,lambda4,4)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda3,3)-OP(1,lambda4,4)*OP(4,lambda2,2)*pp(2,3)*O(lambda1,1,lambda3,3)+OP(1,lambda4,4)*OP(3,lambda2,2)*pp(2,4)*O(lambda1,1,lambda3,3)-OP(3,lambda2,2)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda4,4)-OP(2,lambda3,3)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda4,4)-OP(3,lambda4,4)*OP(4,lambda1,1)*pp(1,2)*O(lambda2,2,lambda3,3)-OP(2,lambda4,4)*OP(4,lambda1,1)*pp(1,3)*O(lambda2,2,lambda3,3)+OP(2,lambda4,4)*OP(3,lambda1,1)*pp(1,4)*O(lambda2,2,lambda3,3)+OP(2,lambda1,1)*OP(3,lambda4,4)*pp(1,4)*O(lambda2,2,lambda3,3)-OP(1,lambda4,4)*OP(3,lambda1,1)*pp(2,4)*O(lambda2,2,lambda3,3)-OP(1,lambda4,4)*OP(2,lambda1,1)*pp(3,4)*O(lambda2,2,lambda3,3)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-OP(3,lambda1,1)*OP(4,lambda3,3)*pp(1,2)*O(lambda2,2,lambda4,4)+OP(2,lambda3,3)*OP(4,lambda1,1)*pp(1,3)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*OP(4,lambda3,3)*pp(1,3)*O(lambda2,2,lambda4,4)-OP(2,lambda3,3)*OP(3,lambda1,1)*pp(1,4)*O(lambda2,2,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(-(OP(3,lambda4,4)*pp(2,4)*O(lambda1,1,lambda2,2))+OP(2,lambda4,4)*(OP(3,lambda2,2)*OP(4,lambda1,1)+pp(3,4)*O(lambda1,1,lambda2,2))+OP(4,lambda2,2)*pp(2,3)*O(lambda1,1,lambda4,4)-OP(3,lambda2,2)*pp(2,4)*O(lambda1,1,lambda4,4)-OP(4,lambda1,1)*pp(2,3)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*(OP(3,lambda4,4)*OP(4,lambda2,2)-pp(3,4)*O(lambda2,2,lambda4,4)))+OP(3,lambda2,2)*OP(4,lambda1,1)*pp(1,2)*O(lambda3,3,lambda4,4)+OP(3,lambda1,1)*OP(4,lambda2,2)*pp(1,2)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(4,lambda2,2)*pp(1,3)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(3,lambda2,2)*pp(1,4)*O(lambda3,3,lambda4,4)+pp(1,4)*pp(2,3)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+OP(1,lambda2,2)*(OP(3,lambda4,4)*pp(2,4)*O(lambda1,1,lambda3,3)+OP(2,lambda4,4)*(OP(3,lambda1,1)*OP(4,lambda3,3)-pp(3,4)*O(lambda1,1,lambda3,3))+OP(4,lambda3,3)*pp(2,3)*O(lambda1,1,lambda4,4)+OP(2,lambda3,3)*(OP(3,lambda4,4)*OP(4,lambda1,1)-pp(3,4)*O(lambda1,1,lambda4,4))-OP(4,lambda1,1)*pp(2,3)*O(lambda3,3,lambda4,4)-OP(3,lambda1,1)*pp(2,4)*O(lambda3,3,lambda4,4))))/(2*(WMASS**4))

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fm0.ne.0d0 ) then  !fm0 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) -2*fm0*(gw**2)*(WMASS**2)*(OP(3,lambda4,4)*OP(4,lambda3,3)*O(lambda1,1,lambda2,2)+OP(2,lambda4,4)*OP(4,lambda2,2)*O(lambda1,1,lambda3,3)+OP(2,lambda3,3)*OP(3,lambda2,2)*O(lambda1,1,lambda4,4)+OP(1,lambda4,4)*OP(4,lambda1,1)*O(lambda2,2,lambda3,3)-pp(1,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-pp(2,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda3,3)*OP(3,lambda1,1)*O(lambda2,2,lambda4,4)-pp(1,3)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)-pp(2,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda2,2)*OP(2,lambda1,1)*O(lambda3,3,lambda4,4)-pp(1,2)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)-pp(3,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4))

               enddo                      
              enddo
             enddo 
            enddo
         endif


         if( fm1.ne.0d0 ) then  !fm1 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(fm1*(gw**2)*(WMASS**2)*(OP(2,lambda3,3)*OP(3,lambda4,4)*O(lambda1,1,lambda2,2)+OP(2,lambda4,4)*OP(4,lambda3,3)*O(lambda1,1,lambda2,2)+OP(1,lambda2,2)*OP(2,lambda4,4)*O(lambda1,1,lambda3,3)+OP(2,lambda4,4)*OP(3,lambda2,2)*O(lambda1,1,lambda3,3)-OP(1,lambda2,2)*OP(3,lambda4,4)*O(lambda1,1,lambda3,3)+OP(3,lambda4,4)*OP(4,lambda2,2)*O(lambda1,1,lambda3,3)+OP(1,lambda2,2)*OP(2,lambda3,3)*O(lambda1,1,lambda4,4)+OP(2,lambda3,3)*OP(4,lambda2,2)*O(lambda1,1,lambda4,4)-OP(1,lambda2,2)*OP(4,lambda3,3)*O(lambda1,1,lambda4,4)+OP(3,lambda2,2)*OP(4,lambda3,3)*O(lambda1,1,lambda4,4)-OP(2,lambda4,4)*OP(3,lambda1,1)*O(lambda2,2,lambda3,3)-OP(2,lambda1,1)*OP(3,lambda4,4)*O(lambda2,2,lambda3,3)+OP(2,lambda4,4)*OP(4,lambda1,1)*O(lambda2,2,lambda3,3)+OP(3,lambda4,4)*OP(4,lambda1,1)*O(lambda2,2,lambda3,3)-pp(1,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-pp(1,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-pp(2,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-pp(3,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*(-(OP(2,lambda3,3)*O(lambda1,1,lambda2,2))+OP(4,lambda3,3)*O(lambda1,1,lambda2,2)-OP(3,lambda2,2)*O(lambda1,1,lambda3,3)+OP(4,lambda2,2)*O(lambda1,1,lambda3,3)+OP(2,lambda1,1)*O(lambda2,2,lambda3,3)+OP(3,lambda1,1)*O(lambda2,2,lambda3,3))+OP(2,lambda3,3)*OP(3,lambda1,1)*O(lambda2,2,lambda4,4)-OP(2,lambda3,3)*OP(4,lambda1,1)*O(lambda2,2,lambda4,4)-OP(2,lambda1,1)*OP(4,lambda3,3)*O(lambda2,2,lambda4,4)+OP(3,lambda1,1)*OP(4,lambda3,3)*O(lambda2,2,lambda4,4)-pp(1,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)-pp(1,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)-pp(2,3)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)-pp(3,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(-(OP(2,lambda4,4)*O(lambda1,1,lambda2,2))+OP(3,lambda4,4)*O(lambda1,1,lambda2,2)+OP(3,lambda2,2)*O(lambda1,1,lambda4,4)-OP(4,lambda2,2)*O(lambda1,1,lambda4,4)+OP(2,lambda1,1)*O(lambda2,2,lambda4,4)+OP(4,lambda1,1)*O(lambda2,2,lambda4,4))+OP(1,lambda2,2)*OP(3,lambda1,1)*O(lambda3,3,lambda4,4)+OP(2,lambda1,1)*OP(3,lambda2,2)*O(lambda3,3,lambda4,4)+OP(1,lambda2,2)*OP(4,lambda1,1)*O(lambda3,3,lambda4,4)-OP(3,lambda2,2)*OP(4,lambda1,1)*O(lambda3,3,lambda4,4)+OP(2,lambda1,1)*OP(4,lambda2,2)*O(lambda3,3,lambda4,4)-OP(3,lambda1,1)*OP(4,lambda2,2)*O(lambda3,3,lambda4,4)-pp(1,3)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)-pp(1,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)-pp(2,3)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)-pp(2,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)))/2

               enddo                      
              enddo
             enddo 
            enddo
         endif


         if( fm2.ne.0d0 ) then  !fm2 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) -((fm2*(gw**2)*(sw**4)*(ZMASS**4)*(OP(3,lambda4,4)*OP(4,lambda3,3)*O(lambda1,1,lambda2,2)+OP(2,lambda4,4)*OP(4,lambda2,2)*O(lambda1,1,lambda3,3)+OP(2,lambda3,3)*OP(3,lambda2,2)*O(lambda1,1,lambda4,4)+OP(1,lambda4,4)*OP(4,lambda1,1)*O(lambda2,2,lambda3,3)-pp(1,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-pp(2,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda3,3)*OP(3,lambda1,1)*O(lambda2,2,lambda4,4)-pp(1,3)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)-pp(2,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda2,2)*OP(2,lambda1,1)*O(lambda3,3,lambda4,4)-pp(1,2)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)-pp(3,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)))/(WMASS**2))

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fm3.ne.0d0 ) then  !fm3 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(fm3*(gw**2)*(sw**4)*(ZMASS**4)*(OP(2,lambda3,3)*OP(3,lambda4,4)*O(lambda1,1,lambda2,2)+OP(2,lambda4,4)*OP(4,lambda3,3)*O(lambda1,1,lambda2,2)+OP(1,lambda2,2)*OP(2,lambda4,4)*O(lambda1,1,lambda3,3)+OP(2,lambda4,4)*OP(3,lambda2,2)*O(lambda1,1,lambda3,3)-OP(1,lambda2,2)*OP(3,lambda4,4)*O(lambda1,1,lambda3,3)+OP(3,lambda4,4)*OP(4,lambda2,2)*O(lambda1,1,lambda3,3)+OP(1,lambda2,2)*OP(2,lambda3,3)*O(lambda1,1,lambda4,4)+OP(2,lambda3,3)*OP(4,lambda2,2)*O(lambda1,1,lambda4,4)-OP(1,lambda2,2)*OP(4,lambda3,3)*O(lambda1,1,lambda4,4)+OP(3,lambda2,2)*OP(4,lambda3,3)*O(lambda1,1,lambda4,4)-OP(2,lambda4,4)*OP(3,lambda1,1)*O(lambda2,2,lambda3,3)-OP(2,lambda1,1)*OP(3,lambda4,4)*O(lambda2,2,lambda3,3)+OP(2,lambda4,4)*OP(4,lambda1,1)*O(lambda2,2,lambda3,3)+OP(3,lambda4,4)*OP(4,lambda1,1)*O(lambda2,2,lambda3,3)-pp(1,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-pp(1,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-pp(2,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-pp(3,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*(-(OP(2,lambda3,3)*O(lambda1,1,lambda2,2))+OP(4,lambda3,3)*O(lambda1,1,lambda2,2)-OP(3,lambda2,2)*O(lambda1,1,lambda3,3)+OP(4,lambda2,2)*O(lambda1,1,lambda3,3)+OP(2,lambda1,1)*O(lambda2,2,lambda3,3)+OP(3,lambda1,1)*O(lambda2,2,lambda3,3))+OP(2,lambda3,3)*OP(3,lambda1,1)*O(lambda2,2,lambda4,4)-OP(2,lambda3,3)*OP(4,lambda1,1)*O(lambda2,2,lambda4,4)-OP(2,lambda1,1)*OP(4,lambda3,3)*O(lambda2,2,lambda4,4)+OP(3,lambda1,1)*OP(4,lambda3,3)*O(lambda2,2,lambda4,4)-pp(1,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)-pp(1,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)-pp(2,3)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)-pp(3,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(-(OP(2,lambda4,4)*O(lambda1,1,lambda2,2))+OP(3,lambda4,4)*O(lambda1,1,lambda2,2)+OP(3,lambda2,2)*O(lambda1,1,lambda4,4)-OP(4,lambda2,2)*O(lambda1,1,lambda4,4)+OP(2,lambda1,1)*O(lambda2,2,lambda4,4)+OP(4,lambda1,1)*O(lambda2,2,lambda4,4))+OP(1,lambda2,2)*OP(3,lambda1,1)*O(lambda3,3,lambda4,4)+OP(2,lambda1,1)*OP(3,lambda2,2)*O(lambda3,3,lambda4,4)+OP(1,lambda2,2)*OP(4,lambda1,1)*O(lambda3,3,lambda4,4)-OP(3,lambda2,2)*OP(4,lambda1,1)*O(lambda3,3,lambda4,4)+OP(2,lambda1,1)*OP(4,lambda2,2)*O(lambda3,3,lambda4,4)-OP(3,lambda1,1)*OP(4,lambda2,2)*O(lambda3,3,lambda4,4)-pp(1,3)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)-pp(1,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)-pp(2,3)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)-pp(2,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)))/(4*(WMASS**2))

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fm4.ne.0d0 ) then  !fm4 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) -(fm4*(gw**2)*(sw**2)*(ZMASS**2)*(OP(3,lambda4,4)*OP(4,lambda3,3)*O(lambda1,1,lambda2,2)+OP(2,lambda4,4)*OP(4,lambda2,2)*O(lambda1,1,lambda3,3)+OP(2,lambda3,3)*OP(3,lambda2,2)*O(lambda1,1,lambda4,4)+OP(1,lambda4,4)*OP(4,lambda1,1)*O(lambda2,2,lambda3,3)-pp(1,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-pp(2,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda3,3)*OP(3,lambda1,1)*O(lambda2,2,lambda4,4)-pp(1,3)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)-pp(2,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda2,2)*OP(2,lambda1,1)*O(lambda3,3,lambda4,4)-pp(1,2)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)-pp(3,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)))

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fm5.ne.0d0 ) then  !fm5 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(fm5*(gw**2)*(sw**2)*(ZMASS**2)*(-(OP(2,lambda3,3)*OP(3,lambda4,4)*O(lambda1,1,lambda2,2))-OP(2,lambda4,4)*OP(4,lambda3,3)*O(lambda1,1,lambda2,2)-OP(1,lambda2,2)*OP(2,lambda4,4)*O(lambda1,1,lambda3,3)-OP(2,lambda4,4)*OP(3,lambda2,2)*O(lambda1,1,lambda3,3)+OP(1,lambda2,2)*OP(3,lambda4,4)*O(lambda1,1,lambda3,3)-OP(3,lambda4,4)*OP(4,lambda2,2)*O(lambda1,1,lambda3,3)-OP(1,lambda2,2)*OP(2,lambda3,3)*O(lambda1,1,lambda4,4)-OP(2,lambda3,3)*OP(4,lambda2,2)*O(lambda1,1,lambda4,4)+OP(1,lambda2,2)*OP(4,lambda3,3)*O(lambda1,1,lambda4,4)-OP(3,lambda2,2)*OP(4,lambda3,3)*O(lambda1,1,lambda4,4)+OP(2,lambda4,4)*OP(3,lambda1,1)*O(lambda2,2,lambda3,3)+OP(2,lambda1,1)*OP(3,lambda4,4)*O(lambda2,2,lambda3,3)-OP(2,lambda4,4)*OP(4,lambda1,1)*O(lambda2,2,lambda3,3)-OP(3,lambda4,4)*OP(4,lambda1,1)*O(lambda2,2,lambda3,3)+pp(1,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(2,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(3,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*(OP(2,lambda3,3)*O(lambda1,1,lambda2,2)-OP(4,lambda3,3)*O(lambda1,1,lambda2,2)+OP(3,lambda2,2)*O(lambda1,1,lambda3,3)-OP(4,lambda2,2)*O(lambda1,1,lambda3,3)-OP(2,lambda1,1)*O(lambda2,2,lambda3,3)-OP(3,lambda1,1)*O(lambda2,2,lambda3,3))-OP(2,lambda3,3)*OP(3,lambda1,1)*O(lambda2,2,lambda4,4)+OP(2,lambda3,3)*OP(4,lambda1,1)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*OP(4,lambda3,3)*O(lambda2,2,lambda4,4)-OP(3,lambda1,1)*OP(4,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(2,3)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(3,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(OP(2,lambda4,4)*O(lambda1,1,lambda2,2)-OP(3,lambda4,4)*O(lambda1,1,lambda2,2)-OP(3,lambda2,2)*O(lambda1,1,lambda4,4)+OP(4,lambda2,2)*O(lambda1,1,lambda4,4)-OP(2,lambda1,1)*O(lambda2,2,lambda4,4)-OP(4,lambda1,1)*O(lambda2,2,lambda4,4))-OP(1,lambda2,2)*OP(3,lambda1,1)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(3,lambda2,2)*O(lambda3,3,lambda4,4)-OP(1,lambda2,2)*OP(4,lambda1,1)*O(lambda3,3,lambda4,4)+OP(3,lambda2,2)*OP(4,lambda1,1)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(4,lambda2,2)*O(lambda3,3,lambda4,4)+OP(3,lambda1,1)*OP(4,lambda2,2)*O(lambda3,3,lambda4,4)+pp(1,3)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(1,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(2,3)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(2,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)))/4

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fm5hc.ne.0d0 ) then  !fm5hc coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(fm5hc*(gw**2)*(sw**2)*(ZMASS**2)*(-(OP(2,lambda3,3)*OP(3,lambda4,4)*O(lambda1,1,lambda2,2))-OP(2,lambda4,4)*OP(4,lambda3,3)*O(lambda1,1,lambda2,2)-OP(1,lambda2,2)*OP(2,lambda4,4)*O(lambda1,1,lambda3,3)-OP(2,lambda4,4)*OP(3,lambda2,2)*O(lambda1,1,lambda3,3)+OP(1,lambda2,2)*OP(3,lambda4,4)*O(lambda1,1,lambda3,3)-OP(3,lambda4,4)*OP(4,lambda2,2)*O(lambda1,1,lambda3,3)-OP(1,lambda2,2)*OP(2,lambda3,3)*O(lambda1,1,lambda4,4)-OP(2,lambda3,3)*OP(4,lambda2,2)*O(lambda1,1,lambda4,4)+OP(1,lambda2,2)*OP(4,lambda3,3)*O(lambda1,1,lambda4,4)-OP(3,lambda2,2)*OP(4,lambda3,3)*O(lambda1,1,lambda4,4)+OP(2,lambda4,4)*OP(3,lambda1,1)*O(lambda2,2,lambda3,3)+OP(2,lambda1,1)*OP(3,lambda4,4)*O(lambda2,2,lambda3,3)-OP(2,lambda4,4)*OP(4,lambda1,1)*O(lambda2,2,lambda3,3)-OP(3,lambda4,4)*OP(4,lambda1,1)*O(lambda2,2,lambda3,3)+pp(1,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(2,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(3,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*(OP(2,lambda3,3)*O(lambda1,1,lambda2,2)-OP(4,lambda3,3)*O(lambda1,1,lambda2,2)+OP(3,lambda2,2)*O(lambda1,1,lambda3,3)-OP(4,lambda2,2)*O(lambda1,1,lambda3,3)-OP(2,lambda1,1)*O(lambda2,2,lambda3,3)-OP(3,lambda1,1)*O(lambda2,2,lambda3,3))-OP(2,lambda3,3)*OP(3,lambda1,1)*O(lambda2,2,lambda4,4)+OP(2,lambda3,3)*OP(4,lambda1,1)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*OP(4,lambda3,3)*O(lambda2,2,lambda4,4)-OP(3,lambda1,1)*OP(4,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(2,3)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(3,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(OP(2,lambda4,4)*O(lambda1,1,lambda2,2)-OP(3,lambda4,4)*O(lambda1,1,lambda2,2)-OP(3,lambda2,2)*O(lambda1,1,lambda4,4)+OP(4,lambda2,2)*O(lambda1,1,lambda4,4)-OP(2,lambda1,1)*O(lambda2,2,lambda4,4)-OP(4,lambda1,1)*O(lambda2,2,lambda4,4))-OP(1,lambda2,2)*OP(3,lambda1,1)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(3,lambda2,2)*O(lambda3,3,lambda4,4)-OP(1,lambda2,2)*OP(4,lambda1,1)*O(lambda3,3,lambda4,4)+OP(3,lambda2,2)*OP(4,lambda1,1)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(4,lambda2,2)*O(lambda3,3,lambda4,4)+OP(3,lambda1,1)*OP(4,lambda2,2)*O(lambda3,3,lambda4,4)+pp(1,3)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(1,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(2,3)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(2,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)))/4

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fm6.ne.0d0 ) then  !fm6 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) -(fm6*(gw**2)*(WMASS**2)*(OP(3,lambda4,4)*OP(4,lambda3,3)*O(lambda1,1,lambda2,2)+OP(2,lambda4,4)*OP(4,lambda2,2)*O(lambda1,1,lambda3,3)+OP(2,lambda3,3)*OP(3,lambda2,2)*O(lambda1,1,lambda4,4)+OP(1,lambda4,4)*OP(4,lambda1,1)*O(lambda2,2,lambda3,3)-pp(1,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-pp(2,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda3,3)*OP(3,lambda1,1)*O(lambda2,2,lambda4,4)-pp(1,3)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)-pp(2,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda2,2)*OP(2,lambda1,1)*O(lambda3,3,lambda4,4)-pp(1,2)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)-pp(3,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)))

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fm7.ne.0d0 ) then  !fm7 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(fm7*(gw**2)*(WMASS**2)*(-(OP(2,lambda3,3)*OP(3,lambda4,4)*O(lambda1,1,lambda2,2))-OP(2,lambda4,4)*OP(4,lambda3,3)*O(lambda1,1,lambda2,2)-OP(1,lambda2,2)*OP(2,lambda4,4)*O(lambda1,1,lambda3,3)-OP(2,lambda4,4)*OP(3,lambda2,2)*O(lambda1,1,lambda3,3)+OP(1,lambda2,2)*OP(3,lambda4,4)*O(lambda1,1,lambda3,3)-OP(3,lambda4,4)*OP(4,lambda2,2)*O(lambda1,1,lambda3,3)-OP(1,lambda2,2)*OP(2,lambda3,3)*O(lambda1,1,lambda4,4)-OP(2,lambda3,3)*OP(4,lambda2,2)*O(lambda1,1,lambda4,4)+OP(1,lambda2,2)*OP(4,lambda3,3)*O(lambda1,1,lambda4,4)-OP(3,lambda2,2)*OP(4,lambda3,3)*O(lambda1,1,lambda4,4)+OP(2,lambda4,4)*OP(3,lambda1,1)*O(lambda2,2,lambda3,3)+OP(2,lambda1,1)*OP(3,lambda4,4)*O(lambda2,2,lambda3,3)-OP(2,lambda4,4)*OP(4,lambda1,1)*O(lambda2,2,lambda3,3)-OP(3,lambda4,4)*OP(4,lambda1,1)*O(lambda2,2,lambda3,3)+pp(1,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,3)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(2,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(3,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*(OP(2,lambda3,3)*O(lambda1,1,lambda2,2)-OP(4,lambda3,3)*O(lambda1,1,lambda2,2)+OP(3,lambda2,2)*O(lambda1,1,lambda3,3)-OP(4,lambda2,2)*O(lambda1,1,lambda3,3)-OP(2,lambda1,1)*O(lambda2,2,lambda3,3)-OP(3,lambda1,1)*O(lambda2,2,lambda3,3))-OP(2,lambda3,3)*OP(3,lambda1,1)*O(lambda2,2,lambda4,4)+OP(2,lambda3,3)*OP(4,lambda1,1)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*OP(4,lambda3,3)*O(lambda2,2,lambda4,4)-OP(3,lambda1,1)*OP(4,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(2,3)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(3,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(OP(2,lambda4,4)*O(lambda1,1,lambda2,2)-OP(3,lambda4,4)*O(lambda1,1,lambda2,2)-OP(3,lambda2,2)*O(lambda1,1,lambda4,4)+OP(4,lambda2,2)*O(lambda1,1,lambda4,4)-OP(2,lambda1,1)*O(lambda2,2,lambda4,4)-OP(4,lambda1,1)*O(lambda2,2,lambda4,4))-OP(1,lambda2,2)*OP(3,lambda1,1)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(3,lambda2,2)*O(lambda3,3,lambda4,4)-OP(1,lambda2,2)*OP(4,lambda1,1)*O(lambda3,3,lambda4,4)+OP(3,lambda2,2)*OP(4,lambda1,1)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(4,lambda2,2)*O(lambda3,3,lambda4,4)+OP(3,lambda1,1)*OP(4,lambda2,2)*O(lambda3,3,lambda4,4)+pp(1,3)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(1,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(2,3)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(2,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)))/4

               enddo                      
              enddo
             enddo 
            enddo
         endif         

      end





ccccccccccccccc
c zzwmwp-vertex c
ccccccccccccccc

      subroutine zzwmwp_ONanom4(z1,z2,wm1,wp1,M)
c
c This subroutine calculates the vvvv-vertex with anomalous couplings
c and the on-shell polarization vectors for vector bosons
c
c IN:
c     polv1, polv2      : the incoming polvector wz(bosons in the CM frame)
c     polv3, polv4      : the outgoing polvector wz(bosons in the CM frame)
c OUT:
c     M(lambda1,lambda2) 
c
      IMPLICIT NONE

      INTEGER NGRAPHS
      PARAMETER( NGRAPHS=20 )
      DOUBLE PRECISION rZERO, rONE, rHalf
      PARAMETER( rZERO=0.0D0, rONE=1.0D0, rHALF = 0.5D0 )
      DOUBLE COMPLEX cIMAG, cZERO, cONE
      PARAMETER( cIMAG=(0.0D0,1.0D0), cZERO=(0.0D0,0.0D0), cONE=(1.0D0,0.0D0))
      DOUBLE COMPLEX z1(6,-1:1), z2(6,-1:1), wm1(6,-1:1), wp1(6,-1:1)      
      DOUBLE COMPLEX M(-1:1,-1:1,-1:1,-1:1)
      
      DOUBLE COMPLEX OFFeps(6,-1:1,NGRAPHS), OFFepsC(6,-1:1,NGRAPHS)
      DOUBLE PRECISION K(0:3,4)      
      DOUBLE COMPLEX O(-1:1,4,-1:1,4), OP(4,-1:1,4)
      DOUBLE COMPLEXPP(1:4,1:4)
    
      INTEGER I,J,N,MU,NU
      INTEGER LAMBDA1, LAMBDA2, LAMBDA3, LAMBDA4, LAMBDA
      INTEGER L1,L2

      DOUBLE COMPLEX DOTCC, DOTRC,contract_Tjj
      DOUBLE PRECISION DOTRR
      EXTERNAL DOTCC, DOTRC, DOTRR, contract_Tjj
      double precision sw2, sw
#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "gmunu.inc"
#include "transm.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"
      sw2 = sin2w
      sw = sqrt(sin2w)
      OFFeps(6,-1:1,4) = cZero
      K(0:3,1:4) = rZero
      M(-1:1,-1:1,-1:1,-1:1) = cZero
      PP(1:4,1:4) = rZero
      O(-1:1,1:4,-1:1,1:4) = cZero
      OP(1:4,-1:1,1:4) = cZero  
      
      OFFeps(1:6,-1:1,1) = z1(1:6,-1:1)
      OFFeps(1:6,-1:1,2) = z2(1:6,-1:1)
      OFFeps(1:6,-1:1,3) = wm1(1:6,-1:1)
      OFFeps(1:6,-1:1,4) = wp1(1:6,-1:1) 

      do i = 1,4
        K(0,i) = dble(OFFeps(5,-1,i))*(-rONE) !K1(0:3)
        K(1,i) = dble(OFFeps(6,-1,i))*(-rONE) !K2(0:3)
        K(2,i) = dimag(OFFeps(6,-1,i))*(-rONE)!K3(0:3)
        K(3,i) = dimag(OFFeps(5,-1,i))*(-rONE)!K4(0:3)
      enddo




      do i = 1,4
        do j = 1,4
          pp(i,j) = dotrr(k(0,i),k(0,j))
        enddo
      enddo

      do L1 = -1,1
        do L2 = -1,1 
          do i = 1,4
            do j = 1,4
              O(L1,i,L2,j) = dotcc(OFFeps(1,L1,i),OFFeps(1,L2,j))
            enddo
          enddo
        enddo
      enddo

      do LAMBDA = -1,1
        do i = 1,4 
          do j = 1,4
            OP(j,lambda,i) = dotrc(k(0,j),OFFeps(1,lambda,i))
          enddo
        enddo
      enddo


!!!! ===============================================================      
!!!! COUPLINGS CALCULATION
!!!! ===============================================================
 


         if( fs0.ne.0d0 ) then  !fs0 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +fs0*(WMASS**2)*(ZMASS**2)*(O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4))

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fs1.ne.0d0 ) then  !fs1 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +2*fs1*(WMASS**2)*(ZMASS**2)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fs2.ne.0d0 ) then  !fs2 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +fs2*(WMASS**2)*(ZMASS**2)*(O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4))

               enddo                      
              enddo
             enddo 
            enddo
         endif                  

         if( ft0.ne.0d0 ) then  !ft0 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) -8*ft0*(gw**4)*(-1+SW2)*(OP(1,lambda2,2)*OP(2,lambda1,1)-pp(1,2)*O(lambda1,1,lambda2,2))*(OP(3,lambda4,4)*OP(4,lambda3,3)-pp(3,4)*O(lambda3,3,lambda4,4))

               enddo                      
              enddo
             enddo 
            enddo
         endif


         if( ft1.ne.0d0 ) then  !ft1 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +4*ft1*(gw**4)*(-1+SW2)*(OP(2,lambda4,4)*OP(4,lambda2,2)*pp(3,1)*O(lambda1,1,lambda3,3)+OP(2,lambda3,3)*OP(3,lambda2,2)*(-(OP(1,lambda4,4)*OP(4,lambda1,1))+pp(4,1)*O(lambda1,1,lambda4,4))+OP(1,lambda4,4)*OP(4,lambda1,1)*pp(3,2)*O(lambda2,2,lambda3,3)-pp(3,2)*pp(4,1)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-pp(3,1)*pp(4,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*OP(3,lambda1,1)*(-(OP(2,lambda4,4)*OP(4,lambda2,2))+pp(4,2)*O(lambda2,2,lambda4,4)))

               enddo                      
              enddo
             enddo 
            enddo
         endif


         if( ft2.ne.0d0 ) then  !ft2 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) -(ft2*(gw**4)*(-1+SW2)*(OP(1,lambda4,4)*OP(2,lambda1,1)*OP(3,lambda2,2)*OP(4,lambda3,3)-OP(2,lambda4,4)*OP(4,lambda3,3)*pp(3,1)*O(lambda1,1,lambda2,2)-OP(1,lambda4,4)*OP(4,lambda3,3)*pp(3,2)*O(lambda1,1,lambda2,2)+OP(1,lambda4,4)*OP(2,lambda3,3)*pp(3,4)*O(lambda1,1,lambda2,2)+OP(1,lambda3,3)*OP(2,lambda4,4)*pp(3,4)*O(lambda1,1,lambda2,2)-OP(1,lambda2,2)*OP(2,lambda4,4)*pp(3,4)*O(lambda1,1,lambda3,3)-OP(3,lambda2,2)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda4,4)+OP(1,lambda2,2)*OP(4,lambda3,3)*pp(3,2)*O(lambda1,1,lambda4,4)-OP(1,lambda2,2)*OP(2,lambda3,3)*pp(3,4)*O(lambda1,1,lambda4,4)-OP(1,lambda4,4)*OP(2,lambda1,1)*pp(3,4)*O(lambda2,2,lambda3,3)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(3,lambda4,4)*(-(OP(2,lambda3,3)*pp(4,1)*O(lambda1,1,lambda2,2))-OP(1,lambda3,3)*pp(4,2)*O(lambda1,1,lambda2,2)+OP(1,lambda2,2)*pp(4,2)*O(lambda1,1,lambda3,3)+OP(4,lambda2,2)*(OP(1,lambda3,3)*OP(2,lambda1,1)-pp(1,2)*O(lambda1,1,lambda3,3))+OP(2,lambda1,1)*pp(4,1)*O(lambda2,2,lambda3,3)+OP(4,lambda1,1)*(OP(1,lambda2,2)*OP(2,lambda3,3)-pp(1,2)*O(lambda2,2,lambda3,3)))+OP(2,lambda1,1)*OP(4,lambda3,3)*pp(3,1)*O(lambda2,2,lambda4,4)-OP(1,lambda3,3)*OP(2,lambda1,1)*pp(3,4)*O(lambda2,2,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(3,lambda2,2)*OP(4,lambda1,1)*pp(1,2)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(4,lambda2,2)*pp(3,1)*O(lambda3,3,lambda4,4)-OP(1,lambda2,2)*OP(4,lambda1,1)*pp(3,2)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(3,lambda2,2)*pp(4,1)*O(lambda3,3,lambda4,4)+pp(3,2)*pp(4,1)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(3,1)*pp(4,2)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+OP(3,lambda1,1)*(OP(4,lambda3,3)*(OP(1,lambda2,2)*OP(2,lambda4,4)-pp(1,2)*O(lambda2,2,lambda4,4))+(OP(4,lambda2,2)*pp(1,2)-OP(1,lambda2,2)*pp(4,2))*O(lambda3,3,lambda4,4))))

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( ft3.ne.0d0 ) then  !ft3 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +2*ft3*(SW2-1)*(gw**4)*(pp(3,2)*OP(4,lambda2,2)*OP(1,lambda4,4)*O(lambda1,1,lambda3,3)+pp(3,1)*OP(4,lambda2,2)*OP(2,lambda3,3)*O(lambda1,1,lambda4,4)-pp(3,2)*OP(4,lambda2,2)*OP(1,lambda3,3)*O(lambda1,1,lambda4,4)+OP(3,lambda2,2)*(pp(4,1)*OP(2,lambda4,4)*O(lambda1,1,lambda3,3)+pp(4,2)*(OP(1,lambda3,3)*O(lambda1,1,lambda4,4)-OP(1,lambda4,4)*O(lambda1,1,lambda3,3))+OP(4,lambda1,1)*(-OP(1,lambda3,3))*OP(2,lambda4,4))+pp(3,1)*OP(4,lambda1,1)*OP(2,lambda4,4)*O(lambda2,2,lambda3,3)-pp(3,1)*OP(4,lambda1,1)*OP(2,lambda3,3)*O(lambda2,2,lambda4,4)+pp(3,2)*OP(4,lambda1,1)*OP(1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(3,lambda1,1)*(-pp(4,1)*OP(2,lambda4,4)*O(lambda2,2,lambda3,3)+pp(4,2)*OP(1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(4,1)*OP(2,lambda3,3)*O(lambda2,2,lambda4,4)+OP(4,lambda2,2)*OP(2,lambda3,3)*(-OP(1,lambda4,4)))-pp(3,2)*pp(4,1)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)-pp(3,1)*pp(4,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3))
               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( ft4.ne.0d0 ) then  !ft4 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) + (gw**4*ft4*sw**4*ZMASS**2*(-(pp(3, 2)*OP(1, lambda4, 4)*OP(4, lambda2, 2)*O(lambda1, 1, lambda3, 3)) + pp(3, 2)*OP(1, lambda3, 3)*OP(4, lambda2, 2)*O(lambda1, 1, lambda4, 4) - pp(3, 1)*OP(2, lambda3, 3)*OP(4, lambda2, 2)*O(lambda1, 1, lambda4, 4) + OP(3, lambda2, 2)*(OP(1, lambda3, 3)*OP(2, lambda4, 4)*OP(4, lambda1, 1) - pp(4, 1)*OP(2, lambda4, 4)*O(lambda1, 1, lambda3, 3) + pp(4, 2)*(OP(1, lambda4, 4)*O(lambda1, 1, lambda3, 3) - OP(1, lambda3, 3)*O(lambda1, 1, lambda4, 4))) - pp(3, 1)*OP(2, lambda4, 4)*OP(4, lambda1, 1)*O(lambda2, 2, lambda3, 3) + pp(3, 1)*pp(4, 2)*O(lambda1, 1, lambda4, 4)*O(lambda2, 2, lambda3, 3) - pp(3, 2)*OP(1, lambda3, 3)*OP(4, lambda1, 1)*O(lambda2, 2, lambda4, 4) + pp(3, 1)*OP(2, lambda3, 3)*OP(4, lambda1, 1)*O(lambda2, 2, lambda4, 4) + pp(3, 2)*pp(4, 1)*O(lambda1, 1, lambda3, 3)*O(lambda2, 2, lambda4, 4) + OP(3, lambda1, 1)*(OP(1, lambda4, 4)*OP(2, lambda3, 3)*OP(4, lambda2, 2) - pp(4, 2)*OP(1, lambda4, 4)*O(lambda2, 2, lambda3, 3) + pp(4, 1)*OP(2, lambda4, 4)*O(lambda2, 2, lambda3, 3) - pp(4, 1)*OP(2, lambda3, 3)*O(lambda2, 2, lambda4, 4))))/(2*WMASS**2)
               enddo                      
              enddo
             enddo 
            enddo
         endif


         if( ft5.ne.0d0 ) then  !ft5 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(2*ft5*(gw**4)*(sw**4)*(ZMASS**2)*(OP(1,lambda2,2)*OP(2,lambda1,1)-pp(1,2)*O(lambda1,1,lambda2,2))*(OP(3,lambda4,4)*OP(4,lambda3,3)-pp(3,4)*O(lambda3,3,lambda4,4)))/(WMASS**2)

               enddo                      
              enddo
             enddo 
            enddo
         endif


         if( ft6.ne.0d0 ) then  !ft6 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(ft6*(gw**4)*(sw**4)*(ZMASS**2)*(-(OP(2,lambda4,4)*OP(4,lambda2,2)*pp(3,1)*O(lambda1,1,lambda3,3))+OP(2,lambda3,3)*OP(3,lambda2,2)*(OP(1,lambda4,4)*OP(4,lambda1,1)-pp(4,1)*O(lambda1,1,lambda4,4))-OP(1,lambda4,4)*OP(4,lambda1,1)*pp(3,2)*O(lambda2,2,lambda3,3)+pp(3,2)*pp(4,1)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(3,1)*pp(4,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*OP(3,lambda1,1)*(OP(2,lambda4,4)*OP(4,lambda2,2)-pp(4,2)*O(lambda2,2,lambda4,4))))/(WMASS**2)

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( ft7.ne.0d0 ) then  !ft7 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(ft7*(gw**4)*(sw**4)*(ZMASS**2)*(OP(1,lambda4,4)*OP(2,lambda1,1)*OP(3,lambda2,2)*OP(4,lambda3,3)-OP(2,lambda4,4)*OP(4,lambda3,3)*pp(3,1)*O(lambda1,1,lambda2,2)-OP(1,lambda4,4)*OP(4,lambda3,3)*pp(3,2)*O(lambda1,1,lambda2,2)+OP(1,lambda4,4)*OP(2,lambda3,3)*pp(3,4)*O(lambda1,1,lambda2,2)+OP(1,lambda3,3)*OP(2,lambda4,4)*pp(3,4)*O(lambda1,1,lambda2,2)-OP(1,lambda2,2)*OP(2,lambda4,4)*pp(3,4)*O(lambda1,1,lambda3,3)-OP(3,lambda2,2)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda4,4)+OP(1,lambda2,2)*OP(4,lambda3,3)*pp(3,2)*O(lambda1,1,lambda4,4)-OP(1,lambda2,2)*OP(2,lambda3,3)*pp(3,4)*O(lambda1,1,lambda4,4)-OP(1,lambda4,4)*OP(2,lambda1,1)*pp(3,4)*O(lambda2,2,lambda3,3)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(3,lambda4,4)*(-(OP(2,lambda3,3)*pp(4,1)*O(lambda1,1,lambda2,2))-OP(1,lambda3,3)*pp(4,2)*O(lambda1,1,lambda2,2)+OP(1,lambda2,2)*pp(4,2)*O(lambda1,1,lambda3,3)+OP(4,lambda2,2)*(OP(1,lambda3,3)*OP(2,lambda1,1)-pp(1,2)*O(lambda1,1,lambda3,3))+OP(2,lambda1,1)*pp(4,1)*O(lambda2,2,lambda3,3)+OP(4,lambda1,1)*(OP(1,lambda2,2)*OP(2,lambda3,3)-pp(1,2)*O(lambda2,2,lambda3,3)))+OP(2,lambda1,1)*OP(4,lambda3,3)*pp(3,1)*O(lambda2,2,lambda4,4)-OP(1,lambda3,3)*OP(2,lambda1,1)*pp(3,4)*O(lambda2,2,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(3,lambda2,2)*OP(4,lambda1,1)*pp(1,2)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(4,lambda2,2)*pp(3,1)*O(lambda3,3,lambda4,4)-OP(1,lambda2,2)*OP(4,lambda1,1)*pp(3,2)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(3,lambda2,2)*pp(4,1)*O(lambda3,3,lambda4,4)+pp(3,2)*pp(4,1)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(3,1)*pp(4,2)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+OP(3,lambda1,1)*(OP(4,lambda3,3)*(OP(1,lambda2,2)*OP(2,lambda4,4)-pp(1,2)*O(lambda2,2,lambda4,4))+(OP(4,lambda2,2)*pp(1,2)-OP(1,lambda2,2)*pp(4,2))*O(lambda3,3,lambda4,4))))/(4*(WMASS**2))

               enddo                      
              enddo
             enddo 
            enddo
         endif



         if( fm0.ne.0d0 ) then  !fm0 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) -2*fm0*(gw**2)*(ZMASS**2)*(((-1+SW2)**2)*OP(1,lambda2,2)*OP(2,lambda1,1)*O(lambda3,3,lambda4,4)+O(lambda1,1,lambda2,2)*(OP(3,lambda4,4)*OP(4,lambda3,3)-(((-1+SW2)**2)*pp(1,2)+pp(3,4))*O(lambda3,3,lambda4,4)))

               enddo                      
              enddo
             enddo 
            enddo
         endif


         if( fm1.ne.0d0 ) then  !fm1 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) -(fm1*(gw**2)*(ZMASS**2)*(OP(1,lambda4,4)*OP(2,lambda3,3)*O(lambda1,1,lambda2,2)+OP(1,lambda3,3)*OP(2,lambda4,4)*O(lambda1,1,lambda2,2)-OP(1,lambda2,2)*OP(2,lambda4,4)*O(lambda1,1,lambda3,3)-OP(3,lambda4,4)*OP(4,lambda2,2)*O(lambda1,1,lambda3,3)-OP(1,lambda2,2)*OP(2,lambda3,3)*O(lambda1,1,lambda4,4)-OP(3,lambda2,2)*OP(4,lambda3,3)*O(lambda1,1,lambda4,4)-OP(1,lambda4,4)*OP(2,lambda1,1)*O(lambda2,2,lambda3,3)-OP(3,lambda4,4)*OP(4,lambda1,1)*O(lambda2,2,lambda3,3)+pp(1,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(3,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-OP(1,lambda3,3)*OP(2,lambda1,1)*O(lambda2,2,lambda4,4)-OP(3,lambda1,1)*OP(4,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(3,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+2*SW2*(-(OP(1,lambda3,3)*OP(2,lambda4,4)*O(lambda1,1,lambda2,2))+OP(1,lambda2,2)*(OP(2,lambda4,4)*O(lambda1,1,lambda3,3)+OP(2,lambda3,3)*O(lambda1,1,lambda4,4))-pp(1,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*(-(OP(2,lambda3,3)*O(lambda1,1,lambda2,2))+OP(2,lambda1,1)*O(lambda2,2,lambda3,3))+OP(1,lambda3,3)*OP(2,lambda1,1)*O(lambda2,2,lambda4,4)-pp(1,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4))+(sw**4)*(OP(1,lambda3,3)*OP(2,lambda4,4)*O(lambda1,1,lambda2,2)-OP(1,lambda2,2)*(OP(2,lambda4,4)*O(lambda1,1,lambda3,3)+OP(2,lambda3,3)*O(lambda1,1,lambda4,4))+pp(1,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*(OP(2,lambda3,3)*O(lambda1,1,lambda2,2)-OP(2,lambda1,1)*O(lambda2,2,lambda3,3))-OP(1,lambda3,3)*OP(2,lambda1,1)*O(lambda2,2,lambda4,4)+pp(1,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4))+OP(3,lambda2,2)*OP(4,lambda1,1)*O(lambda3,3,lambda4,4)+OP(3,lambda1,1)*OP(4,lambda2,2)*O(lambda3,3,lambda4,4)))/2

               enddo                      
              enddo
             enddo 
            enddo
         endif


         if( fm2.ne.0d0 ) then  !fm2 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) -(fm2*(gw**2)*(sw**4)*(ZMASS**2)*(OP(1,lambda2,2)*OP(2,lambda1,1)-pp(1,2)*O(lambda1,1,lambda2,2))*O(lambda3,3,lambda4,4))

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fm3.ne.0d0 ) then  !fm3 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(fm3*(gw**2)*(sw**4)*(ZMASS**2)*(-(OP(1,lambda3,3)*OP(2,lambda4,4)*O(lambda1,1,lambda2,2))+OP(1,lambda2,2)*(OP(2,lambda4,4)*O(lambda1,1,lambda3,3)+OP(2,lambda3,3)*O(lambda1,1,lambda4,4))-pp(1,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*(-(OP(2,lambda3,3)*O(lambda1,1,lambda2,2))+OP(2,lambda1,1)*O(lambda2,2,lambda3,3))+OP(1,lambda3,3)*OP(2,lambda1,1)*O(lambda2,2,lambda4,4)-pp(1,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)))/4

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fm4.ne.0d0 ) then  !fm4 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) -(fm4*(gw**2)*(sw**2)*(ZMASS**2)*(OP(2,lambda4,4)*OP(4,lambda2,2)*O(lambda1,1,lambda3,3)+OP(2,lambda3,3)*OP(3,lambda2,2)*O(lambda1,1,lambda4,4)+OP(1,lambda4,4)*OP(4,lambda1,1)*O(lambda2,2,lambda3,3)-pp(3,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-pp(4,1)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda3,3)*OP(3,lambda1,1)*O(lambda2,2,lambda4,4)-pp(3,1)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)-pp(4,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+2*(-1+SW2)*OP(1,lambda2,2)*OP(2,lambda1,1)*O(lambda3,3,lambda4,4)+2*pp(1,2)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)-2*SW2*pp(1,2)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)))/2

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fm5.ne.0d0 ) then  !fm5 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) -(fm5*(gw**2)*(sw**2)*(ZMASS**2)*(OP(1,lambda4,4)*OP(2,lambda3,3)*O(lambda1,1,lambda2,2)+OP(1,lambda3,3)*OP(2,lambda4,4)*O(lambda1,1,lambda2,2)+OP(1,lambda4,4)*OP(4,lambda3,3)*O(lambda1,1,lambda2,2)+OP(2,lambda4,4)*OP(4,lambda3,3)*O(lambda1,1,lambda2,2)-OP(1,lambda2,2)*OP(2,lambda4,4)*O(lambda1,1,lambda3,3)-OP(1,lambda4,4)*OP(3,lambda2,2)*O(lambda1,1,lambda3,3)+OP(2,lambda4,4)*OP(3,lambda2,2)*O(lambda1,1,lambda3,3)-OP(1,lambda2,2)*OP(2,lambda3,3)*O(lambda1,1,lambda4,4)+OP(1,lambda3,3)*OP(3,lambda2,2)*O(lambda1,1,lambda4,4)-OP(1,lambda2,2)*OP(4,lambda3,3)*O(lambda1,1,lambda4,4)-OP(1,lambda4,4)*OP(2,lambda1,1)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*OP(3,lambda1,1)*O(lambda2,2,lambda3,3)-OP(2,lambda4,4)*OP(3,lambda1,1)*O(lambda2,2,lambda3,3)+pp(1,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-pp(3,1)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-OP(1,lambda3,3)*OP(2,lambda1,1)*O(lambda2,2,lambda4,4)+OP(2,lambda3,3)*OP(3,lambda1,1)*O(lambda2,2,lambda4,4)-OP(2,lambda1,1)*OP(4,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)-pp(3,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+SW2*(-(OP(1,lambda3,3)*OP(2,lambda4,4)*O(lambda1,1,lambda2,2))+OP(1,lambda2,2)*(OP(2,lambda4,4)*O(lambda1,1,lambda3,3)+OP(2,lambda3,3)*O(lambda1,1,lambda4,4))-pp(1,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*(-(OP(2,lambda3,3)*O(lambda1,1,lambda2,2))+OP(2,lambda1,1)*O(lambda2,2,lambda3,3))+OP(1,lambda3,3)*OP(2,lambda1,1)*O(lambda2,2,lambda4,4)-pp(1,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4))+OP(1,lambda2,2)*OP(4,lambda1,1)*O(lambda3,3,lambda4,4)+OP(2,lambda1,1)*OP(4,lambda2,2)*O(lambda3,3,lambda4,4)-pp(4,1)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)-pp(4,2)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)))/4

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fm5hc.ne.0d0 ) then  !fm5hc coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) -(fm5hc*(gw**2)*(sw**2)*(ZMASS**2)*(OP(1,lambda4,4)*OP(2,lambda3,3)*O(lambda1,1,lambda2,2)+OP(1,lambda3,3)*OP(2,lambda4,4)*O(lambda1,1,lambda2,2)+OP(1,lambda3,3)*OP(3,lambda4,4)*O(lambda1,1,lambda2,2)+OP(2,lambda3,3)*OP(3,lambda4,4)*O(lambda1,1,lambda2,2)-OP(1,lambda2,2)*OP(2,lambda4,4)*O(lambda1,1,lambda3,3)-OP(1,lambda2,2)*OP(3,lambda4,4)*O(lambda1,1,lambda3,3)+OP(1,lambda4,4)*OP(4,lambda2,2)*O(lambda1,1,lambda3,3)-OP(1,lambda2,2)*OP(2,lambda3,3)*O(lambda1,1,lambda4,4)-OP(1,lambda3,3)*OP(4,lambda2,2)*O(lambda1,1,lambda4,4)+OP(2,lambda3,3)*OP(4,lambda2,2)*O(lambda1,1,lambda4,4)-OP(1,lambda4,4)*OP(2,lambda1,1)*O(lambda2,2,lambda3,3)-OP(2,lambda1,1)*OP(3,lambda4,4)*O(lambda2,2,lambda3,3)+OP(2,lambda4,4)*OP(4,lambda1,1)*O(lambda2,2,lambda3,3)+pp(1,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-pp(4,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-OP(1,lambda3,3)*OP(2,lambda1,1)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*OP(4,lambda1,1)*O(lambda2,2,lambda4,4)-OP(2,lambda3,3)*OP(4,lambda1,1)*O(lambda2,2,lambda4,4)+pp(1,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)-pp(4,1)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+SW2*(-(OP(1,lambda3,3)*OP(2,lambda4,4)*O(lambda1,1,lambda2,2))+OP(1,lambda2,2)*(OP(2,lambda4,4)*O(lambda1,1,lambda3,3)+OP(2,lambda3,3)*O(lambda1,1,lambda4,4))-pp(1,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*(-(OP(2,lambda3,3)*O(lambda1,1,lambda2,2))+OP(2,lambda1,1)*O(lambda2,2,lambda3,3))+OP(1,lambda3,3)*OP(2,lambda1,1)*O(lambda2,2,lambda4,4)-pp(1,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4))+OP(1,lambda2,2)*OP(3,lambda1,1)*O(lambda3,3,lambda4,4)+OP(2,lambda1,1)*OP(3,lambda2,2)*O(lambda3,3,lambda4,4)-pp(3,1)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)-pp(3,2)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)))/4

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fm6.ne.0d0 ) then  !fm6 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) -(fm6*(gw**2)*(ZMASS**2)*(((-1+SW2)**2)*OP(1,lambda2,2)*OP(2,lambda1,1)*O(lambda3,3,lambda4,4)+O(lambda1,1,lambda2,2)*(OP(3,lambda4,4)*OP(4,lambda3,3)-(((-1+SW2)**2)*pp(1,2)+pp(3,4))*O(lambda3,3,lambda4,4))))

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fm7.ne.0d0 ) then  !fm7 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(fm7*(gw**2)*(ZMASS**2)*(OP(1,lambda4,4)*OP(2,lambda3,3)*O(lambda1,1,lambda2,2)+OP(1,lambda3,3)*OP(2,lambda4,4)*O(lambda1,1,lambda2,2)+OP(1,lambda3,3)*OP(3,lambda4,4)*O(lambda1,1,lambda2,2)+OP(2,lambda3,3)*OP(3,lambda4,4)*O(lambda1,1,lambda2,2)+OP(1,lambda4,4)*OP(4,lambda3,3)*O(lambda1,1,lambda2,2)+OP(2,lambda4,4)*OP(4,lambda3,3)*O(lambda1,1,lambda2,2)-OP(1,lambda2,2)*OP(2,lambda4,4)*O(lambda1,1,lambda3,3)-OP(1,lambda2,2)*OP(3,lambda4,4)*O(lambda1,1,lambda3,3)-OP(1,lambda4,4)*OP(4,lambda2,2)*O(lambda1,1,lambda3,3)-OP(3,lambda4,4)*OP(4,lambda2,2)*O(lambda1,1,lambda3,3)-OP(1,lambda2,2)*OP(2,lambda3,3)*O(lambda1,1,lambda4,4)+OP(1,lambda3,3)*OP(4,lambda2,2)*O(lambda1,1,lambda4,4)-OP(2,lambda3,3)*OP(4,lambda2,2)*O(lambda1,1,lambda4,4)-OP(1,lambda2,2)*OP(4,lambda3,3)*O(lambda1,1,lambda4,4)-OP(1,lambda4,4)*OP(2,lambda1,1)*O(lambda2,2,lambda3,3)-OP(2,lambda1,1)*OP(3,lambda4,4)*O(lambda2,2,lambda3,3)-OP(2,lambda4,4)*OP(4,lambda1,1)*O(lambda2,2,lambda3,3)-OP(3,lambda4,4)*OP(4,lambda1,1)*O(lambda2,2,lambda3,3)+pp(1,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(3,1)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(3,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(4,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-OP(1,lambda3,3)*OP(2,lambda1,1)*O(lambda2,2,lambda4,4)-OP(1,lambda3,3)*OP(4,lambda1,1)*O(lambda2,2,lambda4,4)+OP(2,lambda3,3)*OP(4,lambda1,1)*O(lambda2,2,lambda4,4)-OP(2,lambda1,1)*OP(4,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(3,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(3,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(4,1)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+(sw**4)*(OP(1,lambda3,3)*OP(2,lambda4,4)*O(lambda1,1,lambda2,2)-OP(1,lambda2,2)*(OP(2,lambda4,4)*O(lambda1,1,lambda3,3)+OP(2,lambda3,3)*O(lambda1,1,lambda4,4))+pp(1,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*(OP(2,lambda3,3)*O(lambda1,1,lambda2,2)-OP(2,lambda1,1)*O(lambda2,2,lambda3,3))-OP(1,lambda3,3)*OP(2,lambda1,1)*O(lambda2,2,lambda4,4)+pp(1,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4))+OP(1,lambda2,2)*OP(4,lambda1,1)*O(lambda3,3,lambda4,4)+OP(2,lambda1,1)*OP(4,lambda2,2)*O(lambda3,3,lambda4,4)-pp(3,1)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)-pp(3,2)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)-pp(4,1)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)-pp(4,2)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+OP(3,lambda2,2)*(OP(1,lambda4,4)*O(lambda1,1,lambda3,3)-OP(2,lambda4,4)*O(lambda1,1,lambda3,3)-OP(1,lambda3,3)*O(lambda1,1,lambda4,4)-OP(4,lambda3,3)*O(lambda1,1,lambda4,4)+SW2*(-(OP(1,lambda4,4)*O(lambda1,1,lambda3,3))+OP(2,lambda4,4)*O(lambda1,1,lambda3,3)+OP(1,lambda3,3)*O(lambda1,1,lambda4,4))-(-1+SW2)*OP(2,lambda1,1)*O(lambda3,3,lambda4,4)+OP(4,lambda1,1)*O(lambda3,3,lambda4,4))+OP(3,lambda1,1)*(-(OP(1,lambda4,4)*O(lambda2,2,lambda3,3))+OP(2,lambda4,4)*O(lambda2,2,lambda3,3)-OP(2,lambda3,3)*O(lambda2,2,lambda4,4)-OP(4,lambda3,3)*O(lambda2,2,lambda4,4)+SW2*(OP(1,lambda4,4)*O(lambda2,2,lambda3,3)-OP(2,lambda4,4)*O(lambda2,2,lambda3,3)+OP(2,lambda3,3)*O(lambda2,2,lambda4,4))-(-1+SW2)*OP(1,lambda2,2)*O(lambda3,3,lambda4,4)+OP(4,lambda2,2)*O(lambda3,3,lambda4,4))+SW2*(-2*OP(1,lambda4,4)*OP(2,lambda3,3)*O(lambda1,1,lambda2,2)-2*OP(1,lambda3,3)*OP(2,lambda4,4)*O(lambda1,1,lambda2,2)-OP(1,lambda3,3)*OP(3,lambda4,4)*O(lambda1,1,lambda2,2)-OP(2,lambda3,3)*OP(3,lambda4,4)*O(lambda1,1,lambda2,2)-OP(1,lambda4,4)*OP(4,lambda3,3)*O(lambda1,1,lambda2,2)-OP(2,lambda4,4)*OP(4,lambda3,3)*O(lambda1,1,lambda2,2)+2*OP(1,lambda2,2)*OP(2,lambda4,4)*O(lambda1,1,lambda3,3)+OP(1,lambda2,2)*OP(3,lambda4,4)*O(lambda1,1,lambda3,3)+2*OP(1,lambda2,2)*OP(2,lambda3,3)*O(lambda1,1,lambda4,4)+OP(1,lambda2,2)*OP(4,lambda3,3)*O(lambda1,1,lambda4,4)+2*OP(1,lambda4,4)*OP(2,lambda1,1)*O(lambda2,2,lambda3,3)+OP(2,lambda1,1)*OP(3,lambda4,4)*O(lambda2,2,lambda3,3)-2*pp(1,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-pp(3,1)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-pp(4,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+2*OP(1,lambda3,3)*OP(2,lambda1,1)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*OP(4,lambda3,3)*O(lambda2,2,lambda4,4)-2*pp(1,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)-pp(3,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)-pp(4,1)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(3,1)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(3,2)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(4,1)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+pp(4,2)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+OP(4,lambda1,1)*(OP(2,lambda4,4)*O(lambda2,2,lambda3,3)+(OP(1,lambda3,3)-OP(2,lambda3,3))*O(lambda2,2,lambda4,4)-OP(1,lambda2,2)*O(lambda3,3,lambda4,4))+OP(4,lambda2,2)*(OP(1,lambda4,4)*O(lambda1,1,lambda3,3)+(-OP(1,lambda3,3)+OP(2,lambda3,3))*O(lambda1,1,lambda4,4)-OP(2,lambda1,1)*O(lambda3,3,lambda4,4)))))/4

               enddo                      
              enddo
             enddo 
            enddo
         endif         

      end





ccccccccccccccc
c wpwmwmwp-vertex c
ccccccccccccccc

      subroutine wpwmwmwp_ONanom4(wp1,wm1,wm2,wp2,M)
c
c This subroutine calculates the vvvv-vertex with anomalous couplings
c and the on-shell polarization vectors for vector bosons
c
c IN:
c     polv1, polv2      : the incoming polvector wz(bosons in the CM frame)
c     polv3, polv4      : the outgoing polvector wz(bosons in the CM frame)
c OUT:
c     M(lambda1,lambda2) 
c
      IMPLICIT NONE

      INTEGER NGRAPHS
      PARAMETER( NGRAPHS=20 )
      DOUBLE PRECISION rZERO, rONE, rHalf
      PARAMETER( rZERO=0.0D0, rONE=1.0D0, rHALF = 0.5D0 )
      DOUBLE COMPLEX cIMAG, cZERO, cONE
      PARAMETER( cIMAG=(0.0D0,1.0D0), cZERO=(0.0D0,0.0D0), cONE=(1.0D0,0.0D0))
      DOUBLE COMPLEX wp1(6,-1:1), wm1(6,-1:1), wm2(6,-1:1), wp2(6,-1:1)      
      DOUBLE COMPLEX M(-1:1,-1:1,-1:1,-1:1)
      
      DOUBLE COMPLEX OFFeps(6,-1:1,NGRAPHS), OFFepsC(6,-1:1,NGRAPHS)
      DOUBLE PRECISION K(0:3,4)      
      DOUBLE COMPLEX O(-1:1,4,-1:1,4), OP(4,-1:1,4)
      DOUBLE COMPLEXPP(1:4,1:4)
    
      INTEGER I,J,N,MU,NU
      INTEGER LAMBDA1, LAMBDA2, LAMBDA3, LAMBDA4, LAMBDA
      INTEGER L1,L2

      DOUBLE COMPLEX DOTCC, DOTRC,contract_Tjj
      DOUBLE PRECISION DOTRR
      EXTERNAL DOTCC, DOTRC, DOTRR, contract_Tjj
      double precision sw2, sw
#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "gmunu.inc"
#include "transm.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"
      sw2 = sin2w
      sw = sqrt(sin2w)
      OFFeps(6,-1:1,4) = cZero
      K(0:3,1:4) = rZero
      M(-1:1,-1:1,-1:1,-1:1) = cZero
      PP(1:4,1:4) = rZero
      O(-1:1,1:4,-1:1,1:4) = cZero
      OP(1:4,-1:1,1:4) = cZero  
      
      OFFeps(1:6,-1:1,1) = wp1(1:6,-1:1)
      OFFeps(1:6,-1:1,2) = wm1(1:6,-1:1)
      OFFeps(1:6,-1:1,3) = wm2(1:6,-1:1)
      OFFeps(1:6,-1:1,4) = wp2(1:6,-1:1) 

      do i = 1,4
        K(0,i) = dble(OFFeps(5,-1,i))*(-rONE) !K1(0:3)
        K(1,i) = dble(OFFeps(6,-1,i))*(-rONE) !K2(0:3)
        K(2,i) = dimag(OFFeps(6,-1,i))*(-rONE)!K3(0:3)
        K(3,i) = dimag(OFFeps(5,-1,i))*(-rONE)!K4(0:3)
      enddo
      



      do i = 1,4
        do j = 1,4
          pp(i,j) = dotrr(k(0,i),k(0,j))
        enddo
      enddo

      do L1 = -1,1
        do L2 = -1,1 
          do i = 1,4
            do j = 1,4
              O(L1,i,L2,j) = dotcc(OFFeps(1,L1,i),OFFeps(1,L2,j))
            enddo
          enddo
        enddo
      enddo

      do LAMBDA = -1,1
        do i = 1,4 
          do j = 1,4
            OP(j,lambda,i) = dotrc(k(0,j),OFFeps(1,lambda,i))
          enddo
        enddo
      enddo


!!!! ===============================================================      
!!!! COUPLINGS CALCULATION
!!!! ===============================================================





         if( fs0.ne.0d0 ) then  !fs0 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +4*fs0*(WMASS**4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fs1.ne.0d0 ) then  !fs1 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +2*fs1*(WMASS**4)*(O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4))

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fs2.ne.0d0 ) then  !fs2 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +2*fs2*(WMASS**4)*(O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4))

               enddo                      
              enddo
             enddo 
            enddo
         endif                  

         if( ft0.ne.0d0 ) then  !ft0 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +8*ft0*(gw**4)*(-(OP(3,lambda4,4)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda2,2))-OP(2,lambda3,3)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda4,4)+pp(1,4)*pp(3,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*OP(4,lambda1,1)*(OP(2,lambda3,3)*OP(3,lambda2,2)-pp(3,2)*O(lambda2,2,lambda3,3))+pp(1,2)*pp(3,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+OP(1,lambda2,2)*OP(2,lambda1,1)*(OP(3,lambda4,4)*OP(4,lambda3,3)-pp(3,4)*O(lambda3,3,lambda4,4)))

               enddo                      
              enddo
             enddo 
            enddo
         endif


         if( ft1.ne.0d0 ) then  !ft1 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +4*ft1*(gw**4)*(2*OP(1,lambda3,3)*OP(2,lambda4,4)*OP(3,lambda1,1)*OP(4,lambda2,2)-OP(3,lambda4,4)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda2,2)-2*OP(2,lambda4,4)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda3,3)-OP(2,lambda3,3)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda4,4)+pp(1,4)*pp(3,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda4,4)*OP(4,lambda1,1)*(OP(2,lambda3,3)*OP(3,lambda2,2)-pp(3,2)*O(lambda2,2,lambda3,3))-2*OP(1,lambda3,3)*OP(3,lambda1,1)*pp(2,4)*O(lambda2,2,lambda4,4)+2*pp(1,3)*pp(2,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+OP(1,lambda2,2)*OP(2,lambda1,1)*(OP(3,lambda4,4)*OP(4,lambda3,3)-pp(3,4)*O(lambda3,3,lambda4,4)))

               enddo                      
              enddo
             enddo 
            enddo
         endif


         if( ft2.ne.0d0 ) then  !ft2 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +ft2*(gw**4)*(OP(1,lambda4,4)*OP(2,lambda3,3)*OP(3,lambda1,1)*OP(4,lambda2,2)+2*OP(1,lambda4,4)*OP(2,lambda1,1)*OP(3,lambda2,2)*OP(4,lambda3,3)-OP(2,lambda4,4)*OP(4,lambda3,3)*pp(1,3)*O(lambda1,1,lambda2,2)-2*OP(2,lambda3,3)*OP(3,lambda4,4)*pp(1,4)*O(lambda1,1,lambda2,2)-2*OP(1,lambda4,4)*OP(4,lambda3,3)*pp(3,2)*O(lambda1,1,lambda2,2)+2*OP(1,lambda4,4)*OP(2,lambda3,3)*pp(3,4)*O(lambda1,1,lambda2,2)-OP(3,lambda4,4)*OP(4,lambda2,2)*pp(1,2)*O(lambda1,1,lambda3,3)-OP(2,lambda4,4)*OP(3,lambda2,2)*pp(1,4)*O(lambda1,1,lambda3,3)+OP(1,lambda4,4)*OP(3,lambda2,2)*pp(2,4)*O(lambda1,1,lambda3,3)-OP(1,lambda4,4)*OP(4,lambda2,2)*pp(3,2)*O(lambda1,1,lambda3,3)-2*OP(3,lambda2,2)*OP(4,lambda3,3)*pp(1,2)*O(lambda1,1,lambda4,4)-OP(2,lambda3,3)*OP(4,lambda2,2)*pp(1,3)*O(lambda1,1,lambda4,4)-2*OP(3,lambda4,4)*OP(4,lambda1,1)*pp(1,2)*O(lambda2,2,lambda3,3)-OP(2,lambda4,4)*OP(4,lambda1,1)*pp(1,3)*O(lambda2,2,lambda3,3)+OP(2,lambda4,4)*OP(3,lambda1,1)*pp(1,4)*O(lambda2,2,lambda3,3)+2*OP(2,lambda1,1)*OP(3,lambda4,4)*pp(1,4)*O(lambda2,2,lambda3,3)-OP(1,lambda4,4)*OP(3,lambda1,1)*pp(2,4)*O(lambda2,2,lambda3,3)-2*OP(1,lambda4,4)*OP(2,lambda1,1)*pp(3,4)*O(lambda2,2,lambda3,3)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+2*pp(1,2)*pp(3,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-OP(3,lambda1,1)*OP(4,lambda3,3)*pp(1,2)*O(lambda2,2,lambda4,4)+OP(2,lambda3,3)*OP(4,lambda1,1)*pp(1,3)*O(lambda2,2,lambda4,4)+OP(2,lambda1,1)*OP(4,lambda3,3)*pp(1,3)*O(lambda2,2,lambda4,4)-OP(2,lambda3,3)*OP(3,lambda1,1)*pp(1,4)*O(lambda2,2,lambda4,4)+pp(1,4)*pp(3,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(OP(2,lambda4,4)*pp(3,4)*O(lambda1,1,lambda2,2)+OP(3,lambda4,4)*(OP(2,lambda1,1)*OP(4,lambda2,2)-pp(2,4)*O(lambda1,1,lambda2,2))+OP(4,lambda2,2)*pp(3,2)*O(lambda1,1,lambda4,4)+OP(3,lambda2,2)*(OP(2,lambda4,4)*OP(4,lambda1,1)-pp(2,4)*O(lambda1,1,lambda4,4))-OP(4,lambda1,1)*pp(3,2)*O(lambda2,2,lambda4,4)-OP(2,lambda1,1)*pp(3,4)*O(lambda2,2,lambda4,4))+2*OP(3,lambda2,2)*OP(4,lambda1,1)*pp(1,2)*O(lambda3,3,lambda4,4)+OP(3,lambda1,1)*OP(4,lambda2,2)*pp(1,2)*O(lambda3,3,lambda4,4)-OP(2,lambda1,1)*OP(4,lambda2,2)*pp(1,3)*O(lambda3,3,lambda4,4)-2*OP(2,lambda1,1)*OP(3,lambda2,2)*pp(1,4)*O(lambda3,3,lambda4,4)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+2*pp(1,4)*pp(3,2)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+OP(1,lambda2,2)*(-(OP(2,lambda4,4)*pp(3,4)*O(lambda1,1,lambda3,3))+OP(3,lambda4,4)*(2*OP(2,lambda3,3)*OP(4,lambda1,1)+pp(2,4)*O(lambda1,1,lambda3,3))+2*OP(4,lambda3,3)*pp(3,2)*O(lambda1,1,lambda4,4)-2*OP(2,lambda3,3)*pp(3,4)*O(lambda1,1,lambda4,4)-2*OP(4,lambda1,1)*pp(3,2)*O(lambda3,3,lambda4,4)+OP(3,lambda1,1)*(OP(2,lambda4,4)*OP(4,lambda3,3)-pp(2,4)*O(lambda3,3,lambda4,4))))

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( ft3.ne.0d0 ) then  !ft3 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +2.0d0*ft3*(gw**4)*(-pp(1,2)*OP(4,lambda2,2)*OP(3,lambda4,4)*O(lambda1,1,lambda3,3)-pp(3,2)*OP(4,lambda2,2)*OP(1,lambda4,4)*O(lambda1,1,lambda3,3)-pp(1,3)*OP(4,lambda2,2)*OP(2,lambda3,3)*O(lambda1,1,lambda4,4)+pp(1,2)*OP(3,lambda1,1)*OP(4,lambda2,2)*O(lambda3,3,lambda4,4)-pp(1,3)*OP(2,lambda1,1)*OP(4,lambda2,2)*O(lambda3,3,lambda4,4)-pp(1,3)*OP(4,lambda3,3)*OP(2,lambda4,4)*O(lambda1,1,lambda2,2)-pp(1,4)*OP(3,lambda2,2)*OP(2,lambda4,4)*O(lambda1,1,lambda3,3)+pp(2,4)*OP(3,lambda2,2)*OP(1,lambda4,4)*O(lambda1,1,lambda3,3)-pp(1,3)*OP(4,lambda1,1)*OP(2,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,4)*OP(3,lambda1,1)*OP(2,lambda4,4)*O(lambda2,2,lambda3,3)-pp(2,4)*OP(3,lambda1,1)*OP(1,lambda4,4)*O(lambda2,2,lambda3,3)-pp(1,2)*OP(3,lambda1,1)*OP(4,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,3)*OP(4,lambda1,1)*OP(2,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,3)*OP(2,lambda1,1)*OP(4,lambda3,3)*O(lambda2,2,lambda4,4)-pp(1,4)*OP(3,lambda1,1)*OP(2,lambda3,3)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*(pp(3,4)*OP(2,lambda4,4)*O(lambda1,1,lambda2,2)+OP(3,lambda4,4)*(OP(2,lambda1,1)*OP(4,lambda2,2)-pp(2,4)*O(lambda1,1,lambda2,2))+pp(3,2)*OP(4,lambda2,2)*O(lambda1,1,lambda4,4)+OP(3,lambda2,2)*(OP(4,lambda1,1)*OP(2,lambda4,4)-pp(2,4)*O(lambda1,1,lambda4,4))-pp(3,2)*OP(4,lambda1,1)*O(lambda2,2,lambda4,4)-pp(3,4)*OP(2,lambda1,1)*O(lambda2,2,lambda4,4))+OP(1,lambda2,2)*(O(lambda1,1,lambda3,3)*(pp(2,4)*OP(3,lambda4,4)-pp(3,4)*OP(2,lambda4,4))+OP(3,lambda1,1)*(OP(4,lambda3,3)*OP(2,lambda4,4)-pp(2,4)*O(lambda3,3,lambda4,4)))+OP(3,lambda1,1)*OP(4,lambda2,2)*OP(2,lambda3,3)*OP(1,lambda4,4)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+pp(1,4)*pp(3,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,2)*pp(3,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)+pp(1,3)*pp(2,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4))

               enddo                      
              enddo
             enddo 
            enddo
         endif
         
         if( fm0.ne.0d0 ) then  !fm0 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) -2*fm0*(gw**2)*(WMASS**2)*(OP(3,lambda4,4)*OP(4,lambda3,3)*O(lambda1,1,lambda2,2)+OP(2,lambda3,3)*OP(3,lambda2,2)*O(lambda1,1,lambda4,4)+OP(1,lambda4,4)*OP(4,lambda1,1)*O(lambda2,2,lambda3,3)-pp(1,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-pp(3,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda2,2)*OP(2,lambda1,1)*O(lambda3,3,lambda4,4)-pp(1,2)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)-pp(3,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4))

               enddo                      
              enddo
             enddo 
            enddo
         endif


         if( fm1.ne.0d0 ) then  !fm1 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) +(fm1*(gw**2)*(WMASS**2)*(-(OP(1,lambda4,4)*OP(2,lambda3,3)*O(lambda1,1,lambda2,2))-OP(1,lambda3,3)*OP(2,lambda4,4)*O(lambda1,1,lambda2,2)+OP(2,lambda3,3)*OP(3,lambda4,4)*O(lambda1,1,lambda2,2)+OP(1,lambda4,4)*OP(4,lambda3,3)*O(lambda1,1,lambda2,2)+OP(1,lambda4,4)*OP(4,lambda2,2)*O(lambda1,1,lambda3,3)+OP(3,lambda4,4)*OP(4,lambda2,2)*O(lambda1,1,lambda3,3)-OP(1,lambda3,3)*OP(4,lambda2,2)*O(lambda1,1,lambda4,4)+OP(1,lambda4,4)*OP(2,lambda1,1)*O(lambda2,2,lambda3,3)-OP(2,lambda4,4)*OP(3,lambda1,1)*O(lambda2,2,lambda3,3)-OP(2,lambda1,1)*OP(3,lambda4,4)*O(lambda2,2,lambda3,3)+OP(3,lambda4,4)*OP(4,lambda1,1)*O(lambda2,2,lambda3,3)-pp(1,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-pp(3,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda3,3)*OP(2,lambda1,1)*O(lambda2,2,lambda4,4)+OP(2,lambda3,3)*OP(3,lambda1,1)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*OP(4,lambda1,1)*O(lambda2,2,lambda4,4)+OP(3,lambda1,1)*OP(4,lambda3,3)*O(lambda2,2,lambda4,4)-pp(1,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)-pp(1,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)-pp(3,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)-pp(3,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)-OP(3,lambda1,1)*OP(4,lambda2,2)*O(lambda3,3,lambda4,4)-pp(1,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)-pp(3,2)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)+OP(3,lambda2,2)*(OP(2,lambda4,4)*O(lambda1,1,lambda3,3)+OP(4,lambda3,3)*O(lambda1,1,lambda4,4)+(OP(2,lambda1,1)-OP(4,lambda1,1))*O(lambda3,3,lambda4,4))+OP(1,lambda2,2)*(OP(2,lambda4,4)*O(lambda1,1,lambda3,3)+OP(2,lambda3,3)*O(lambda1,1,lambda4,4)-OP(4,lambda3,3)*O(lambda1,1,lambda4,4)+OP(4,lambda1,1)*O(lambda3,3,lambda4,4))))/2

               enddo                      
              enddo
             enddo 
            enddo
         endif



         if( fm6.ne.0d0 ) then  !fm6 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) -(fm6*(gw**2)*(WMASS**2)*(OP(3,lambda4,4)*OP(4,lambda3,3)*O(lambda1,1,lambda2,2)+OP(2,lambda3,3)*OP(3,lambda2,2)*O(lambda1,1,lambda4,4)+OP(1,lambda4,4)*OP(4,lambda1,1)*O(lambda2,2,lambda3,3)-pp(1,4)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)-pp(3,2)*O(lambda1,1,lambda4,4)*O(lambda2,2,lambda3,3)+OP(1,lambda2,2)*OP(2,lambda1,1)*O(lambda3,3,lambda4,4)-pp(1,2)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)-pp(3,4)*O(lambda1,1,lambda2,2)*O(lambda3,3,lambda4,4)))

               enddo                      
              enddo
             enddo 
            enddo
         endif

         if( fm7.ne.0d0 ) then  !fm7 coupling
            do lambda1 = -1,1
             do lambda2 = -1,1
              do lambda3 = -1,1
               do lambda4 = -1,1
                  M(lambda1,lambda2,lambda3,lambda4) = M(lambda1,lambda2,lambda3,lambda4) -(fm7*(gw**2)*(WMASS**2)*(-(OP(1,lambda3,3)*OP(2,lambda4,4)*O(lambda1,1,lambda2,2))+OP(1,lambda2,2)*OP(2,lambda4,4)*O(lambda1,1,lambda3,3)+OP(2,lambda4,4)*OP(3,lambda2,2)*O(lambda1,1,lambda3,3)+OP(1,lambda4,4)*OP(4,lambda2,2)*O(lambda1,1,lambda3,3)+OP(3,lambda4,4)*OP(4,lambda2,2)*O(lambda1,1,lambda3,3)-OP(1,lambda3,3)*OP(4,lambda2,2)*O(lambda1,1,lambda4,4)-OP(2,lambda4,4)*OP(3,lambda1,1)*O(lambda2,2,lambda3,3)+OP(1,lambda3,3)*OP(2,lambda1,1)*O(lambda2,2,lambda4,4)+OP(2,lambda3,3)*OP(3,lambda1,1)*O(lambda2,2,lambda4,4)+OP(1,lambda3,3)*OP(4,lambda1,1)*O(lambda2,2,lambda4,4)+OP(3,lambda1,1)*OP(4,lambda3,3)*O(lambda2,2,lambda4,4)-pp(1,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)-pp(1,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)-pp(3,2)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)-pp(3,4)*O(lambda1,1,lambda3,3)*O(lambda2,2,lambda4,4)-OP(3,lambda1,1)*OP(4,lambda2,2)*O(lambda3,3,lambda4,4)))/2

               enddo                      
              enddo
             enddo 
            enddo
         endif         

      end
