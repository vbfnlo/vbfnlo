      subroutine uzzxxx_sing(z1,z2 , uzz)
c
c This subroutine computes an off-shell tensor current for a spin-2 singlet tensor boson
c from the z-z-tensor boson coupling.
c
c input:
c       complex z1(6)          : Z boson 1                              z1
c       complex z2(6)          : Z boson 2                              z2

c output:
c       complex uzz(18)        : tensor current         j^mu^nu(T:z1,z2)
c     
c- by Jessica Frank
c
      implicit none
      double complex z1(6), z2(6), uzz(18)

      integer i, j, k, l, m, n, a, s

      double complex yzz(6,4)
      double precision p1(4), p2(4), pT(4), p1p2, Leci(4,4,4,4)
      double precision MET(4,4), B(4,4,4,4), T1(4,4,4,4), mf(4) 
      double complex cZero, d
      double precision rZero, rTwo, rOne
      double precision pT2, T2(4,4,4,4), T3(4,4,4,4)
      double precision T2part1(4,4,4,4), T2part2(4,4,4,4)

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
 
      p1(1) = dreal(z1(5))
      p1(2) = dreal(z1(6))
      p1(3) = dimag(z1(6))
      p1(4) = dimag(z1(5))
      p2(1) = dreal(z2(5))
      p2(2) = dreal(z2(6))
      p2(3) = dimag(z2(6))
      p2(4) = dimag(z2(5))

      
      yzz(5,1) = z1(5)+z2(5)
      yzz(6,1) = z1(6)+z2(6)

      pT(1) = dreal(yzz(5,1))
      pT(2) = dreal(yzz(6,1))
      pT(3) = dimag(yzz(6,1))
      pT(4) = dimag(yzz(5,1))
      
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


       if ((f3.ne.0.0d0) .or. (f4.ne.0.0d0)) then

        do i = 1,4
         do j = 1,4
            do k = 1,4
               do l = 1,4
          Leci(i,j,k,l) = rZero
              enddo
           enddo
         enddo
      enddo


      Leci(1,2,3,4) = Leci(1,2,3,4)+rOne
      Leci(1,2,4,3) = Leci(1,2,4,3)-rOne
      Leci(1,3,2,4) = Leci(1,3,2,4)-rOne
      Leci(1,3,4,2) = Leci(1,3,4,2)+rOne
      Leci(1,4,2,3) = Leci(1,4,2,3)+rOne  
      Leci(1,4,3,2) = Leci(1,4,3,2)-rOne 
      Leci(2,1,3,4) = Leci(2,1,3,4)-rOne 
      Leci(2,1,4,3) = Leci(2,1,4,3)+rOne 
      Leci(2,3,1,4) = Leci(2,3,1,4)+rOne 
      Leci(2,3,4,1) = Leci(2,3,4,1)-rOne 
      Leci(2,4,1,3) = Leci(2,4,1,3)-rOne 
      Leci(2,4,3,1) = Leci(2,4,3,1)+rOne 
      Leci(3,1,2,4) = Leci(3,1,2,4)+rOne 
      Leci(3,1,4,2) = Leci(3,1,4,2)-rOne 
      Leci(3,2,1,4) = Leci(3,2,1,4)-rOne 
      Leci(3,2,4,1) = Leci(3,2,4,1)+rOne 
      Leci(3,4,1,2) = Leci(3,4,1,2)+rOne  
      Leci(3,4,2,1) = Leci(3,4,2,1)-rOne 
      Leci(4,1,2,3) = Leci(4,1,2,3)-rOne 
      Leci(4,1,3,2) = Leci(4,1,3,2)+rOne 
      Leci(4,2,1,3) = Leci(4,2,1,3)+rOne 
      Leci(4,2,3,1) = Leci(4,2,3,1)-rOne   
      Leci(4,3,1,2) = Leci(4,3,1,2)-rOne 
      Leci(4,3,2,1) = Leci(4,3,2,1)+rOne 

       endif

      p1p2 = p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)-p1(4)*p2(4)


      do i = 1,4
         do j = 1,4
            do k = 1,4
               do l = 1,4

 
      T1(i,j,k,l) = p1(l)*p2(k)*MET(i,j)-p1(j)*p2(l)*MET(i,k)
     &             -p2(i)*p1(l)*MET(j,k)+p1p2*MET(i,l)*MET(j,k)


       if ((f3.ne.0.0d0) .or. (f4.ne.0.0d0)) then

      T2part1(i,j,k,l) = rZero

      do a = 1,4
         do s = 1,4

      T2part1(i,j,k,l) = T2part1(i,j,k,l)+(p2(a)*p1(s)*MET(i,k)
     &    *Leci(j,l,s,a)+p1(a)*p2(s)*MET(j,k)*Leci(i,l,s,a))
     &    *mf(a)*mf(s)

         enddo
      enddo

      T2part2(i,j,k,l) = rZero

      do a = 1,4

      T2part2(i,j,k,l) = T2part2(i,j,k,l)+((p1(a)*p2(k)
     &                    -p2(a)*p1(k))*Leci(i,j,l,a))*mf(a)

      enddo

      T2(i,j,k,l) = T2part1(i,j,k,l)+T2part2(i,j,k,l)

       endif


      T3(i,j,k,l) = MET(i,l)*MET(j,k)


      pT2  = pT(1)**2-pT(2)**2-pT(3)**2-pT(4)**2
 
      if ( sp2mass.gt.rZero ) then
         d = rOne/dcmplx( pT2-sp2mass**2, sp2mass*sp2width )
      else
         d = rOne/dcmplx( pT2, rZero )
      end if


      B(i,j,k,l) =
     &   (rOne/rTwo*(MET(i,k)*MET(j,l)+MET(i,l)*MET(j,k)
     &   -MET(i,j)*MET(k,l))-rOne/(rTwo*sp2mass**2)*(MET(i,k)*pT(j)
     &   *pT(l)+MET(j,l)*pT(i)*pT(k)+MET(i,l)*pT(j)*pT(k)
     &   +MET(j,k)*pT(i)*pT(l))+rOne/6.0d0*(MET(i,j)
     &   +rTwo/sp2mass**2*pT(i)*pT(j))*(MET(k,l)
     &   +rTwo/sp2mass**2*pT(k)*pT(l)))

              enddo
           enddo
         enddo
      enddo

      do i = 1,4
         do j = 1,4  

      yzz(i,j) = cZero       

      do k = 1,4
         do l = 1,4
            do m = 1,4
               do n = 1,4

       if ((f3.EQ.0.0d0) .and. (f4.EQ.0.0d0)) then

      yzz(i,j) = yzz(i,j)-(B(k,l,i,j)*d*z1(m)*z2(n)*
     &         rOne/lambda*
     &         (rTwo*(f2*cos2w+f1*sin2w)*T1(m,n,k,l)
     &         +f5*(gw**2+gw1**2)*v**2/rTwo*T3(m,n,k,l))) 
     &         *mf(k)*mf(l)*mf(m)*mf(n) 

       else

      yzz(i,j) = yzz(i,j)-(B(k,l,i,j)*d*z1(m)*z2(n)*
     &         rOne/lambda*
     &         (rTwo*(f2*cos2w+f1*sin2w)*T1(m,n,k,l)
     &         +(f4*cos2w+f3*sin2w)*T2(m,n,k,l)
     &         +f5*(gw**2+gw1**2)*v**2/rTwo*T3(m,n,k,l))) 
     &         *mf(k)*mf(l)*mf(m)*mf(n) 

       endif

               enddo
            enddo
         enddo
        enddo

       enddo
      enddo

      uzz(1) = yzz(1,1)
      uzz(2) = yzz(1,2)
      uzz(3) = yzz(1,3)
      uzz(4) = yzz(1,4)
      uzz(5) = yzz(2,1)
      uzz(6) = yzz(2,2)
      uzz(7) = yzz(2,3)
      uzz(8) = yzz(2,4)
      uzz(9) = yzz(3,1)
      uzz(10) = yzz(3,2)
      uzz(11) = yzz(3,3)
      uzz(12) = yzz(3,4)
      uzz(13) = yzz(4,1)
      uzz(14) = yzz(4,2)
      uzz(15) = yzz(4,3)
      uzz(16) = yzz(4,4)
      uzz(17) = yzz(5,1)
      uzz(18) = yzz(6,1)

      return
      end
