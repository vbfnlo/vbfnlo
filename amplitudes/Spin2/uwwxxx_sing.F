      subroutine uwwxxx_sing(wp,wm , uww)
c
c This subroutine computes an off-shell tensor current for a spin-2 singlet tensor boson
c from the w+, w-, tensor boson coupling.
c
c input:
c       complex wp(6)          : W+                              wp
c       complex wm(6)          : W-                              wm

c
c output:
c       complex uww(18)        : tensor current         j^mu^nu(T:wp,wm)
c     
c- by Jessica Frank
c
      implicit none
      double complex wp(6), wm(6) 

      integer i, j, k, l, m, n, a, s

      double complex yww(6,4), uww(18)
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
      
      p1(1) = dreal(wp(5))
      p1(2) = dreal(wp(6))
      p1(3) = dimag(wp(6))
      p1(4) = dimag(wp(5))
      p2(1) = dreal(wm(5))
      p2(2) = dreal(wm(6))
      p2(3) = dimag(wm(6))
      p2(4) = dimag(wm(5))

      
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
     &   -MET(i,j)*MET(k,l))-rOne/(rTwo*(sp2mass**2))*(MET(i,k)*pT(j)
     &   *pT(l)+MET(j,l)*pT(i)*pT(k)+MET(i,l)*pT(j)*pT(k)
     &   +MET(j,k)*pT(i)*pT(l))+rOne/6.0d0*(MET(i,j)
     &   +rTwo/(sp2mass**2)*pT(i)*pT(j))*(MET(k,l)
     &   +rTwo/(sp2mass**2)*pT(k)*pT(l)))


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

       if ((f3.EQ.0.0d0) .and. (f4.EQ.0.0d0)) then

       yww(i,j) = yww(i,j)-(B(k,l,i,j)*d*wp(m)*wm(n)*
     &         rOne/lambda*(rTwo*f2*T1(m,n,k,l)
     &         +f5*gw**2*v**2/rTwo*T3(m,n,k,l)))
     &         *mf(k)*mf(l)*mf(m)*mf(n) 

       else

       yww(i,j) = yww(i,j)-(B(k,l,i,j)*d*wp(m)*wm(n)*
     &         rOne/lambda*
     &         (rTwo*f2*T1(m,n,k,l)+f4*T2(m,n,k,l)
     &         +f5*gw**2*v**2/rTwo*T3(m,n,k,l)))
     &         *mf(k)*mf(l)*mf(m)*mf(n) 

       endif
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
