c************************************************************************


      subroutine Wardtest_Pent(psi, p1, p2, minus_qz1, minus_qz2, 
     1                            minus_qz3, z1, z2, z3, scale, alfas,
     2                            j5contracted, mborn, gaugetest)

c****************************************************************
c
c       Vera Hankele, <vera@particle.uni-karlsruhe.de>
c	Initial version:  2008
c	Last modified:  March 2009
c
c this subroutine computes the born and virtual matrix elements with 
c Z3, A3 replaced by q_3 and compares it with the analytic result
c which consists of a sum of diagrams with 2 vector bosons
c in the final state.
c
c***************************************************************

      implicit none

c  input / output
      double complex psi(2,-1:1,2)
      double precision p1(0:3), p2(0:3)
      double precision minus_qz1(0:3), minus_qz2(0:3), minus_qz3(0:3)
      double complex z1(6), z2(6), z3(6)
      double precision scale
      double precision alfas
      double complex j5contracted(3), mborn
      integer gaugetest

c  local variables
      integer mu, number
      double complex mvcheck(1:3)
      double precision pi, colorfactor, factor, temp
      double precision couplfactor
      double precision fac(3), deltaP(3)
      parameter (colorfactor = 4.0d0/3.0d0)
      parameter (couplfactor = 0.01d0)
      parameter (pi = 3.141592653589793d0)

      double precision minus_qz1z2(0:3), minus_qz2z3(0:3)
      double complex gaugebox1(2),gaugebox1_born(2)

      double complex result(-1:1,8),result_born(-1:1,8)

      integer icount, bcount, lbox
      data lbox /0/
      data bcount /0/
      data icount /0/
      save icount, bcount, lbox

      double complex zero
      parameter (zero=(0d0,0d0))

      logical ldebug
      parameter (ldebug = .false.)

c  external functions
      double complex s1c, dotcc, dotrc
      external s1c, dotcc, dotrc

c--------------------------- Begin Code ----------------------------
c eps_mu ~ q_mu * fac, fac = |eps(0)|/q(0)
      fac(1) = abs(z1(1))/abs(minus_qz1(0))
      fac(2) = abs(z2(1))/abs(minus_qz2(0))
      fac(3) = abs(z3(1))/abs(minus_qz3(0))

      do mu = 0,3
         minus_qz1z2(mu) = minus_qz1(mu) + minus_qz2(mu)
         minus_qz2z3(mu) = minus_qz2(mu) + minus_qz3(mu)
      enddo

      call boxlinemm(p1,minus_qz1z2,minus_qz3,p2,
     1   psi(1,-1,2),psi(1,-1,1),z2(1), z3(1),
     2   -1,scale,9,1,gaugebox1,result(-1,1),gaugebox1_born,
     3   result_born(-1,1))

      call boxlinemm(p1,minus_qz1z2,minus_qz3,p2,
     1   psi(1,-1,2),psi(1,-1,1),z1(1), z3(1),
     2   -1,scale,9,-1,gaugebox1,result(-1,2),gaugebox1_born,
     3   result_born(-1,2))

      call boxlinemm(p1,minus_qz1,minus_qz2z3,p2,
     1   psi(1,-1,2),psi(1,-1,1),z1(1), z3(1),
     2   -1,scale,9,1,gaugebox1,result(-1,3),gaugebox1_born,
     3   result_born(-1,3))

      call boxlinemm(p1,minus_qz1,minus_qz2z3,p2,
     1   psi(1,-1,2),psi(1,-1,1),z1(1), z2(1),
     2   -1,scale,9,-1,gaugebox1,result(-1,4),gaugebox1_born,
     3   result_born(-1,4))

      mvcheck(1) = -result(-1,1)
      mvcheck(2) = result(-1,2)-result(-1,3)
      mvcheck(3) = result(-1,4)

      deltaP(1) = abs(fac(1)*(j5contracted(1)-mvcheck(1)))
      deltaP(2) = abs(fac(2)*(j5contracted(2)-mvcheck(2)))
      deltaP(3) = abs(fac(3)*(j5contracted(3)-mvcheck(3)))

c select the largest pentagon contribution
      temp = max(deltaP(1), deltaP(2), deltaP(3))
      if (deltaP(1).eq.temp) number = 1
      if (deltaP(2).eq.temp) number = 2
      if (deltaP(3).eq.temp) number = 3

      factor = colorfactor*couplfactor*fac(number)*alfas/2.0d0/pi

      if( abs((j5contracted(number)-mvcheck(number))/mborn*factor)
     1                          .gt.1.0d-3)then
         icount = icount+1
         gaugetest = 1
      else
         bcount = bcount+1
         gaugetest = 0
      endif      

c      do i = 10,10000000,1000
c      if (bcount.eq.1d2*i) print*,"ratio= ",icount*1.0d0/(bcount+icount),
c     1                              lbox*1.0d0/(bcount+icount)
c      enddo

      return
      end

c************************************************************************


      subroutine Wardtest_Box(z1, qz1, z2, qz2, boxcontr, borncontr, 
     &                                          mborn,gaugetest,alfas)

c****************************************************************
c
c       Vera Hankele, <vera@particle.uni-karlsruhe.de>
c	Initial version:  February 2009
c	Last modified:  February 2009
c
c this subroutine computes the born and virtual matrix elements with 
c the external polarization vectors replaced by their momentum
c and compares it with the analytic result.
c
c***************************************************************

      implicit none

c  input / output
      double precision qz1(0:4), qz2(0:4)
      double complex z1(6), z2(6)
      double complex boxcontr(2), borncontr(2), mborn
      integer gaugetest
      double precision alfas

c  local variables
      integer i ,number
      double complex mvcheck(2), cte
      double precision pi, colorfactor, factor, temp
      double precision couplfactor
      double precision fac(2), deltaP(2)
      parameter (colorfactor = 4.0d0/3.0d0)
      parameter (couplfactor = 0.1d0)
      parameter (pi = 3.141592653589793d0)

c      integer icount, bcount
c      data bcount /0/
c      data icount /0/
c      save icount, bcount

      double complex zero
      parameter (zero=(0d0,0d0))

c--------------------------- Begin Code ----------------------------
c eps_mu ~ q_mu * fac, fac = |eps(0)|/q(0)
      fac(1) = abs(z1(1))/qz1(0)
      fac(2) = abs(z2(1))/qz2(0)

      cte = (3*pi*(0d0,1d0)-4.0d0*pi**2/3.0d0+7.0d0)

c Box contribution contracted with external momentum is proportional
c to Born contribution contracted with external momentum
      mvcheck(1) = -cte * borncontr(1) 
      mvcheck(2) = -cte * borncontr(2)

      deltaP(1) = abs(fac(1)*(boxcontr(1)-mvcheck(1)))
      deltaP(2) = abs(fac(2)*(boxcontr(2)-mvcheck(2)))

c select the largest pentagon contribution
            temp = max(deltaP(1), deltaP(2))
            if (deltaP(1).eq.temp) number = 1
            if (deltaP(2).eq.temp) number = 2

            factor = colorfactor*couplfactor*fac(number)*alfas/2.0d0/pi

         if( abs((boxcontr(number)-mvcheck(number))/mborn*factor)
     1                          .gt.1.0d-3)then
c            icount = icount+1
            gaugetest = 1
         else
            gaugetest = 0
         endif      

c         bcount = bcount+1

c      do i = 10,10000000,1000
c      if (bcount.eq.1d2*i) print*,"ratio= ",icount*1.0d0/bcount, icount, bcount
c      enddo

      return
      end

