      subroutine uwaxxx_trippm(w1,a2 , uwa)
c
c This subroutine computes an off-shell tensor current from 
c the coupling of a W boson, a photon and a charged spin-2 triplet tensor boson
c
c input:
c       complex w1(6)          : W                                      w1
c       complex a2(6)          : photon                                 a2

c output:
c       complex uwa(18)        : tensor current (charged spin-2 particle)   j^mu^nu(T:w1,a2)
c     
c- by Jessica Frank
c
      implicit none
      double complex w1(6), a2(6), uwa(18)

      integer i, j, k, l, m, n

      double complex ywa(6,4)
      double precision p1(4), p2(4), pT(4), p1p2, pT2
      double precision MET(4,4), B(4,4,4,4), T1(4,4,4,4), mf(4) 
      double complex cZero, d
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
 
      p1(1) = dreal(w1(5))
      p1(2) = dreal(w1(6))
      p1(3) = dimag(w1(6))
      p1(4) = dimag(w1(5))
      p2(1) = dreal(a2(5))
      p2(2) = dreal(a2(6))
      p2(3) = dimag(a2(6))
      p2(4) = dimag(a2(5))

      
      ywa(5,1) = w1(5)+a2(5)
      ywa(6,1) = w1(6)+a2(6)

      pT(1) = dreal(ywa(5,1))
      pT(2) = dreal(ywa(6,1))
      pT(3) = dimag(ywa(6,1))
      pT(4) = dimag(ywa(5,1))
      
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


      p1p2 = p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)-p1(4)*p2(4)


      do i = 1,4
         do j = 1,4
            do k = 1,4
               do l = 1,4

 
      T1(i,j,k,l) = p1(l)*p2(k)*MET(i,j)-p1(j)*p2(l)*MET(i,k)
     &             -p2(i)*p1(l)*MET(j,k)+p1p2*MET(i,l)*MET(j,k)


      pT2  = pT(1)**2-pT(2)**2-pT(3)**2-pT(4)**2

 
      if ( msp2trippm.gt.rZero ) then
       d = rOne/dcmplx( pT2-msp2trippm**2, msp2trippm*sp2trippmwidth )
      else
         d = rOne/dcmplx( pT2, rZero )
      end if


      B(i,j,k,l) =
     &   (rOne/rTwo*(MET(i,k)*MET(j,l)+MET(i,l)*MET(j,k)
     &   -MET(i,j)*MET(k,l))-rOne/(rTwo*msp2trippm**2)*(MET(i,k)*pT(j)
     &   *pT(l)+MET(j,l)*pT(i)*pT(k)+MET(i,l)*pT(j)*pT(k)
     &   +MET(j,k)*pT(i)*pT(l))+rOne/6.0d0*(MET(i,j)
     &   +rTwo/msp2trippm**2*pT(i)*pT(j))*(MET(k,l)
     &   +rTwo/msp2trippm**2*pT(k)*pT(l)))


              enddo
           enddo
         enddo
      enddo

      do i = 1,4
         do j = 1,4  

      ywa(i,j) = cZero       

      do k = 1,4
         do l = 1,4
            do m = 1,4
               do n = 1,4

      ywa(i,j) = ywa(i,j)-(B(k,l,i,j)*d*w1(m)*a2(n)*
     &         rOne/lambda*
     &         (f7*cosw*T1(m,n,k,l))) 
     &         *mf(k)*mf(l)*mf(m)*mf(n) 

               enddo
            enddo
         enddo
        enddo

       enddo
      enddo

      uwa(1) = ywa(1,1)
      uwa(2) = ywa(1,2)
      uwa(3) = ywa(1,3)
      uwa(4) = ywa(1,4)
      uwa(5) = ywa(2,1)
      uwa(6) = ywa(2,2)
      uwa(7) = ywa(2,3)
      uwa(8) = ywa(2,4)
      uwa(9) = ywa(3,1)
      uwa(10) = ywa(3,2)
      uwa(11) = ywa(3,3)
      uwa(12) = ywa(3,4)
      uwa(13) = ywa(4,1)
      uwa(14) = ywa(4,2)
      uwa(15) = ywa(4,3)
      uwa(16) = ywa(4,4)
      uwa(17) = ywa(5,1)
      uwa(18) = ywa(6,1)

      return
      end
