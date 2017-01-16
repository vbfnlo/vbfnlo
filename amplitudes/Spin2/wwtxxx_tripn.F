      subroutine wwtxxx_tripn(wp,wm,tc , vertex)
c
c This subroutine computes an amplitude of the three-point coupling of
c a W plus,a W minus and a neutral spin-2 triplet tensor boson.
c
c input:
c       complex wp(6)          : W plus                               wp
c       complex wm(6)          : W minus                              wm
c       complex tc(18)         : input  tensor                         T
c
c output:
c       complex vertex         : amplitude                gamma(wp,wm,T)
c
c- by Jessica Frank
c     
      implicit none
      double complex wp(6), wm(6), tc(18) , vertex

      double complex ft(6,4), mf(4)

      double precision MET(4,4), T3(4,4,4,4)

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


      do i = 1,4
         do j = 1,4
            do k = 1,4
               do l = 1,4

      T3(i,j,k,l) = MET(i,l)*MET(j,k)

                 enddo
            enddo
         enddo
      enddo

      vertex = rZero

      do i = 1,4
         do j = 1,4
            do k = 1,4
               do l = 1,4
                
      vertex = vertex+(rOne/lambda*ft(k,l)*wp(i)*wm(j)*
     &         f6*gw**2*v**2/4.0d0*T3(i,j,k,l))
     &         *mf(i)*mf(j)*mf(k)*mf(l) 
      
               enddo
            enddo
         enddo
      enddo
 

      return
      end
