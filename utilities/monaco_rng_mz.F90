module monaco_rng_mz
   type :: monaco_state
      integer*4 ndim
      integer*8 :: i=97, j=33
      double precision u(97), c, cd, cm
      integer*8 :: ncalls
   end type monaco_state

   type(monaco_state), save :: s ! rng state
   type(monaco_state), save :: saved_s ! rng state

   private
   public :: monran, monran_set, monaco_state, imonrn
   public :: monran_print_state
   public :: get_next_rn

contains
   double precision function get_next_rn() result(uni)
      implicit none
      uni = s%u(s%i) - s%u(s%j)
      if ( uni .lt. 0.0d0 ) then
         uni = uni + 1.0d0
      end if
      s%u(s%i) = uni
      s%i = s%i - 1
      if ( s%i .le. 0 ) then
         s%i = 97
      end if
      s%j = s%j - 1
      if ( s%j .le. 0 ) then
         s%j = 97
      end if
      s%c = s%c - s%cd
      if ( s%c .lt. 0.0d0 ) then
         s%c = s%c + s%cm
      end if
      uni = uni - s%c
      if ( uni .lt. 0.0d0 ) then
         uni = uni + 1.0d0
      end if
      s%ncalls = s%ncalls + 1
   end function 

   subroutine monran( rvec )
   implicit none
   double precision, dimension(s%ndim), intent(out) :: rvec
   integer ivec

   do ivec = 1, s%ndim
      rvec(ivec) = get_next_rn()
   enddo
   end subroutine

!--------- allow for restart of random numbers at predetermined point
   subroutine monran_set(id)
       use globalvars, only: lglobalprint
   implicit none
   integer*4 id

   if (id.eq.1) then   ! save to monrncsave
      saved_s = s
      if (lglobalprint) print*," random number setting in monaco saved "
   elseif(id.eq.2) then   ! restore monrnc from  monrncsave
      s = saved_s
      print*," random number setting in monaco restored to saved "
   endif
   end subroutine

   type(monaco_state) function monran_getstate() result(res)
      res = s
   end function 

   subroutine monran_setstate(new_state)
      type(monaco_state) :: new_state
      s=new_state
   end subroutine 

   subroutine monran_print_state()
      print*, 'ncall, i, j', s%ncalls, s%i, s%j
      print*, 'u(1:3), c, cd, cm', s%u(1:3), s%c, s%cd, s%cm
   end subroutine

!   This subroutine is a "universal" random number generator
!   as proposed by Marsaglia and Zaman in report FSU-SCRI-87-50
!   Slightly modified by F. James, 1988, to generate a vector
!   of pseudorandom numbers "rvec" of length "len", and also
!   by A. Duff, 1991.

!   To get the values in the Marsaglia-Zaman paper, put
!   ij = 1802, kl = 9373

   subroutine imonrn( ndim, seed1, seed2)

   implicit none
   integer*4, intent(in) :: seed1, seed2
   integer*4, intent(in) :: ndim
   integer*8 i, j, k, l, m, ij, kl, ii, jj
   double precision si, t 

   ij = 1802 + seed1
   kl = 9373 + seed2

   s%ndim = ndim
   i = mod( ij / 177, 177 ) + 2
   j = mod( ij, 177 ) + 2
   k = mod( kl / 169, 178 ) + 1
   l = mod( kl, 169 )

   do ii = 1, 97
      si = 0.0d0
      t = 0.5d0
      do jj = 1, 24
         m = mod( mod( i * j, 179 ) * k, 179 )
         i = j
         j = k
         k = m
         l = mod( 53 * l + 1, 169 )
         if ( mod( l * m, 64 ) .ge. 32 ) then
            si = si + t
         end if
         t = 0.5d0 * t
      end do
      s%u(ii) = si
   end do

   s%c = 362436.0d0 / 16777216.0d0
   s%cd = 7654321.0d0 / 16777216.0d0
   s%cm = 16777213.0d0 / 16777216.0d0

   s%ncalls = 0

   end subroutine

end module
