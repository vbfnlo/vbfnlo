
c*************************************************************************
      SUBROUTINE InitRandomNumbers(seed)
c*************************************************************************
c     Select and initialize random number generator.
c*************************************************************************
      implicit none

c     So far we have only one random number generator here: The one from
c     monaco. If it is meaningful, other generators might be included later, 
c     being interfaced by this routine and the 'RandomNumber' function.

      integer seed, seed1, seed2
      
      call loadfile("random.dat",.false.)
      call read_int("SEED",seed,0)
      call closefile

      print *," Random number generator initialized with seed = ",seed
      
c set the two seeds
      if (mod(seed,2).eq.0) then
         seed1 = seed / 2
         seed2 = seed / 2
      else
         seed1 = (seed + 1) / 2
         seed2 = (seed - 1) / 2
      endif

      call iranmr( 1, seed1, seed2 )
      
      RETURN
      END

c*************************************************************************
      FUNCTION RandomNumber()
c*************************************************************************
c     Return a random number
c*************************************************************************
      implicit none
#include "global.inc"
      real*8 R, RandomNumber
      
      call ranmar( R )

      RandomNumber=R

      RETURN
      END
