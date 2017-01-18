module monaco_rng_sob
   implicit none
!   This is a stripped down version of Algorithm 659
!   as appeared in ACM Trans. Math. Software, vol.14, no.1
!   March 1988, pgs. 88-100

   save
   double precision recipd
   integer*8 v(40,30), s, maxcol, counter
   integer*8 :: atmost
   integer*8 :: poly(40) 
   public :: imonso, monsob
   
contains
   subroutine imonso( ndim )
      implicit none
      integer*8 i, j, j2, k, l, m, newv
      integer*8 temp
      integer*4, intent(in) :: ndim
      logical  includ(8)

      atmost = 1073741823
      poly(1:40) = (/1, 3, 7, 11, 13, 19, 25, 37, 59, 47, 61, 55, 41, &
              67, 97, 91, 109, 103, 115, 131, 193, 137, 145, 143, &
              241, 157, 185, 167, 229, 171, 213, 191, 253, 203, 211, &
              239, 247, 285, 369, 299/)

      v(1:40,1) = 1
      v(3:40,2) = (/1, 3, 1, 3, 1, 3, 3, 1, 3, 1, 3, 1, 3, &
                      1, 1, 3, 1, 3, 1, 3, 1, 3, 3, 1, 3, 1, &
                      3, 1, 3, 1, 1, 3, 1, 3, 1, 3, 1, 3/)
      v(4:40,3) = (/7, 5, 1, 3, 3, 7, 5, 5, 7, 7, 1, 3, 3, &
                      7, 5, 1, 1, 5, 3, 3, 1, 7, 5, 1, 3, 3, &
                      7, 5, 1, 1, 5, 7, 7, 5, 1, 3, 3/)
      v(6:40,4) = (/1, 7, 9, 13, 11, 1, 3, 7, 9, 5, 13, 13, &
                      11, 3, 15, 5, 3, 15, 7, 9, 13, 9, 1, 11, &
                      7, 5, 15, 1, 15, 11, 5, 3, 1, 7, 9/)
      v(8:40,5) = (/9, 3, 27, 15, 29, 21, 23, 19, 11, 25, 7, &
                      13, 17, 1, 25, 29, 3, 31, 11, 5, 23, 27, &
                      19, 21, 5, 1, 17, 13, 7, 15, 9, 31, 9/)
      v(14:40,6) = (/37, 33, 7, 5, 11, 39, 63, 27, 17, 15, &
                       23, 29, 3, 21, 13, 31, 25, 9, 49, 33, &
                       19, 29, 11, 19, 27, 15, 25/)
      v(20:40,7) = (/13, 33, 115, 41, 79, 17, 29, 119, 75, &
                       73, 105, 7, 59, 65, 21, 3, 113, 61, 89, &
                       45, 107/)
      v(38:40,8) = (/7, 23, 39/)


!      data tau /0, 0, 1, 3, 5, 8, 11, 15, 19, 23, 27, 31, 35/

      s = ndim
      i = atmost
      maxcol = 0
      do while (i > 0)
         maxcol = maxcol + 1
         i = i / 2
      enddo

      v(1,1:maxcol) = 1

      do i = 2, s
         j = poly(i)
         m = 0
         do while(j > 1)
            j = j / 2
            m = m + 1
         enddo

         j = poly(i)
         do k = m, 1, -1
            j2 = j / 2
            includ(k) = (j .ne. (2 * j2))
            j = j2
         end do

         do j = m + 1, maxcol
            newv = v(i,j-m)
            l = 1
            do k = 1, m
               l = 2 * l
               if ( includ(k) ) then
                  temp = l * v(i,j-k)
                  newv = ieor(newv, temp)
               end if
            end do
            v(i,j) = newv
         end do
      end do

      l = 1
      do j = maxcol - 1, 1, -1
         l = 2 * l
         v(1:s,j) = v(1:s,j) * l
      end do

      recipd = 1.0d0 / (2 * l)

      counter = 0
      poly(1:s) = 0
   end subroutine

   subroutine monsob( rvec )
      implicit none
      integer*8 i, i2, l
      double precision rvec(40)

      l = 0
      i = counter
   10   l = l + 1 !TODO rewrite for
      i2 = i / 2
      if ( i .ne. (2 * i2) ) then
         i = i2
         goto 10
      end if

      if ( l .gt. maxcol ) then
         stop "MONACO:  Sobol generator - too many calls."  !not too likely
      end if

      do i = 1, s
         poly(i) = ieor(poly(i), v(i,l))
      end do
      rvec(1:s) = poly(1:s) * recipd

      counter = counter + 1
   end subroutine

end module

