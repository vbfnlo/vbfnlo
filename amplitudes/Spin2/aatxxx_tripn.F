      subroutine aatxxx_tripn(a1,a2,tc , vertex)
c
c This subroutine computes an amplitude of the three-point coupling of
c two photons and a neutral spin-2 triplet tensor boson.
c
c input:
c       complex a1(6)          : Photon 1                          a1
c       complex a2(6)          : Photon 2                          a2
c       complex tc(18)         : input  tensor                     T
c
c output:
c       complex vertex         : amplitude                gamma(a1,a2,T)
c
c- by Jessica Frank
c     
      implicit none
      double complex a1(6), a2(6), tc(18), vertex

      double complex ft(6,4), mf(4)
      double precision p1(4), p2(4), p1p2
      double precision MET(4,4), T1(4,4,4,4)

      integer i, j, k, l

      double complex cZero
      double precision rZero, rTwo, rOne
      parameter( rZero = 0.0d0, rTwo = 2.0d0, rOne = 1.0d0 )
      parameter( cZero = ( 0.0d0, 0.0d0 ) )

      double precision cos2w, cosw, sinw, gw1, v

#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/spin2coupl.inc"
   
      cos2w = 1.0d0 - sin2w
      cosw = sqrt(cos2w)
      sinw = sqrt(sin2w)
      gw1 = ee / cosw
      v = rTwo*wmass/gw

      MET(1,1) =  1.0d0
      MET(2,2) = -1.0d0
      MET(3,3) = -1.0d0
      MET(4,4) = -1.0d0
      MET(1,2) = 0.0d0
      MET(1,3) = 0.0d0
      MET(1,4) = 0.0d0
      MET(2,1) = 0.0d0
      MET(2,3) = 0.0d0
      MET(2,4) = 0.0d0
      MET(3,1) = 0.0d0
      MET(3,2) = 0.0d0
      MET(3,4) = 0.0d0
      MET(4,1) = 0.0d0
      MET(4,2) = 0.0d0
      MET(4,3) = 0.0d0
    
 
c      multiplication factor in order to get the right signs for the summation
c      with upper and lower indices:

       mf(1) = 1.0d0
       mf(2) = -1.0d0
       mf(3) = -1.0d0
       mf(4) = -1.0d0
 

      ft(1,1) = tc(1)
      ft(1,2) = tc(2)
      ft(1,3) = tc(3)
      ft(1,4) = tc(4)
      ft(2,1) = tc(5)
      ft(2,2) = tc(6)
      ft(2,3) = tc(7)
      ft(2,4) = tc(8)
      ft(3,1) = tc(9)
      ft(3,2) = tc(10)
      ft(3,3) = tc(11)
      ft(3,4) = tc(12)
      ft(4,1) = tc(13)
      ft(4,2) = tc(14)
      ft(4,3) = tc(15)
      ft(4,4) = tc(16)
      ft(5,1) = tc(17)
      ft(6,1) = tc(18)

      p1(1) = dreal(a1(5))
      p1(2) = dreal(a1(6))
      p1(3) = dimag(a1(6))
      p1(4) = dimag(a1(5))
      p2(1) = dreal(a2(5))
      p2(2) = dreal(a2(6))
      p2(3) = dimag(a2(6))
      p2(4) = dimag(a2(5))

      p1p2 = p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)-p1(4)*p2(4)

      do i = 1,4
         do j = 1,4
            do k = 1,4
               do l = 1,4
 
      T1(i,j,k,l) = p1(l)*p2(k)*MET(i,j)-p1(j)*p2(l)*MET(i,k)
     &             -p2(i)*p1(l)*MET(j,k)+p1p2*MET(i,l)*MET(j,k)

                 enddo
            enddo
         enddo
      enddo

      vertex = rZero

      do i = 1,4
         do j = 1,4
            do k = 1,4
               do l = 1,4
                
      vertex = vertex+rOne/lambda*ft(k,l)*a1(i)*a2(j)*
     &         rTwo*f7*cosw*sinw*T1(i,j,k,l)
     &         *mf(i)*mf(j)*mf(k)*mf(l) 

               enddo
            enddo
         enddo
      enddo

      return
      end
