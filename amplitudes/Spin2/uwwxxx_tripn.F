      subroutine uwwxxx_tripn(wp,wm , uww)
c
c This subroutine computes an off-shell tensor current 
c from the coupling of a W+, a W- and a neutral spin-2 triplet tensor boson.
c
c input:
c       complex wp(6)          : W+                             wp
c       complex wm(6)          : W-                             wm

c
c output:
c       complex uww(18)        : tensor current (neutral spin-2 triplet)   j^mu^nu(T:wp,wm)
c     
c- by Jessica Frank
c
      implicit none
      double complex wp(6), wm(6) 

      integer i, j, k, l, m, n

      double complex yww(6,4), uww(18)
      double precision pT(4)
      double precision MET(4,4), B(4,4,4,4), mf(4) 
      double complex cZero, d
      double precision rZero, rTwo, rOne
      double precision pT2, T3(4,4,4,4)

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

      
      yww(5,1) = wp(5)+wm(5)
      yww(6,1) = wp(6)+wm(6)

      pT(1) = dreal(yww(5,1))
      pT(2) = dreal(yww(6,1))
      pT(3) = dimag(yww(6,1))
      pT(4) = dimag(yww(5,1))
      
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


      do i = 1,4
         do j = 1,4
            do k = 1,4
               do l = 1,4


      T3(i,j,k,l) = MET(i,l)*MET(j,k)


      pT2  = pT(1)**2-pT(2)**2-pT(3)**2-pT(4)**2


      if ( msp2tripn.gt.rZero ) then
         d = rOne/dcmplx( pT2-msp2tripn**2, msp2tripn*sp2tripnwidth )
      else
         d = rOne/dcmplx( pT2, rZero )
      end if


      B(i,j,k,l) =
     &   (rOne/rTwo*(MET(i,k)*MET(j,l)+MET(i,l)*MET(j,k)
     &  -MET(i,j)*MET(k,l))-rOne/(rTwo*(msp2tripn**2))*(MET(i,k)*pT(j)
     &   *pT(l)+MET(j,l)*pT(i)*pT(k)+MET(i,l)*pT(j)*pT(k)
     &   +MET(j,k)*pT(i)*pT(l))+rOne/6.0d0*(MET(i,j)
     &   +rTwo/(msp2tripn**2)*pT(i)*pT(j))*(MET(k,l)
     &   +rTwo/(msp2tripn**2)*pT(k)*pT(l)))


              enddo
           enddo
         enddo
      enddo


      do i = 1,4
         do j = 1,4  

      yww(i,j) = cZero       

      do k = 1,4
         do l = 1,4
            do m = 1,4
               do n = 1,4


       yww(i,j) = yww(i,j)-(B(k,l,i,j)*d*wp(m)*wm(n)*
     &         rOne/lambda*
     &         f6*gw**2*v**2/4.0d0*T3(m,n,k,l))
     &         *mf(k)*mf(l)*mf(m)*mf(n) 


                     enddo
                  enddo
               enddo
            enddo

         enddo
      enddo


      uww(1) = yww(1,1)
      uww(2) = yww(1,2)
      uww(3) = yww(1,3)
      uww(4) = yww(1,4)
      uww(5) = yww(2,1)
      uww(6) = yww(2,2)
      uww(7) = yww(2,3)
      uww(8) = yww(2,4)
      uww(9) = yww(3,1)
      uww(10) = yww(3,2)
      uww(11) = yww(3,3)
      uww(12) = yww(3,4)
      uww(13) = yww(4,1)
      uww(14) = yww(4,2)
      uww(15) = yww(4,3)
      uww(16) = yww(4,4)
      uww(17) = yww(5,1)
      uww(18) = yww(6,1)


      return
      end
