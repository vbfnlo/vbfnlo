
!*************************************************************************
      SUBROUTINE InitRandomNumbers()
!*************************************************************************
!     Select and initialize random number generator.
!*************************************************************************
          use globalvars, only: lglobalprint, ldoblha, seed
          use readinput
          use monaco, only: rtype
      implicit none

!     So far we have only one random number generator here: The one from
!     monaco. If it is meaningful, other generators might be included later, 
!     being interfaced by this routine and the 'RandomNumber' function.

      integer seed1, seed2
      
      if (.not.ldoblha) then
        call loadfile("random.dat",.false.)

        ! don't read seed if already set via command line argument
        if (seed == 0) then
            call read_int("SEED",seed,0)
        endif
        call read_int("RTYPE", rtype, 0, .true.)

        call closefile

        if (lglobalprint) print *," Random number generator initialized with seed = ",seed
      endif
      
! set the two seeds
      if (mod(seed,2).eq.0) then
         seed1 = seed / 2
         seed2 = seed / 2
      else
         seed1 = (seed + 1) / 2
         seed2 = (seed - 1) / 2
      endif

      call iranmr( seed1, seed2 )
      
      END


!-------------------------------------------------------------
!
!	This subroutine is a "universal" random number generator
!	as proposed by Marsaglia and Zaman in report FSU-SCRI-87-50
!	Slightly modified by F. James, 1988, to generate a vector
!	of pseudorandom numbers "rvec" of length "len", and also
!	by A. Duff, 1991.
!
!	To get the values in the Marsaglia-Zaman paper, put
!	ij = 1802, kl = 9373
!
      subroutine iranmr( seed1, seed2 )
!
      implicit none
      integer*4 seed1, seed2
      integer*8 i, j, k, l, m, ij, kl, ii, jj
      double precision s, t, u(97), c, cd, cm
!
      common /comrmr/ u, c, cd, cm
      save /comrmr/
!
!     data ij /1802/, kl /9373/
      ij = 1802 + seed2
      kl = 9373 + seed1

      i = mod( ij / 177, 177 ) + 2
      j = mod( ij, 177 ) + 2
      k = mod( kl / 169, 178 ) + 1
      l = mod( kl, 169 )
!
      do ii = 1, 97
         s = 0.0d0
         t = 0.5d0
         do jj = 1, 24
            m = mod( mod( i * j, 179 ) * k, 179 )
            i = j
            j = k
            k = m
            l = mod( 53 * l + 1, 169 )
            if ( mod( l * m, 64 ) .ge. 32 ) then
               s = s + t
            end if
            t = 0.5d0 * t
         end do
         u(ii) = s
      end do
!
      c = 362436.0d0 / 16777216.0d0
      cd = 7654321.0d0 / 16777216.0d0
      cm = 16777213.0d0 / 16777216.0d0
!
      return
      end


!*************************************************************************
      FUNCTION RandomNumber()
!*************************************************************************
!     Return a random number
!*************************************************************************
      implicit none
#include "global.inc"
      real*8 RandomNumber

      integer*8 i, j, ivec
      integer*4 ndim
      double precision r, uni
      double precision u(97), c, cd, cm

      common /comrmr/ u, c, cd, cm
      save /comrmr/

      data i /97/, j /33/

      if (lwritedata.or.lreaddata) then
        RandomNumber = 0
        return
      endif 

      uni = u(i) - u(j)
      if ( uni .lt. 0.0d0 ) then
          uni = uni + 1.0d0
      end if
      u(i) = uni
      i = i - 1
      if ( i .le. 0 ) then
          i = 97
      end if
      j = j - 1
      if ( j .le. 0 ) then
          j = 97
      end if
      c = c - cd
      if ( c .lt. 0.0d0 ) then
          c = c + cm
      end if
      uni = uni - c
      if ( uni .lt. 0.0d0 ) then
          uni = uni + 1.0d0
      end if
      r = uni

      RandomNumber = r
      end function

