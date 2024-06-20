****************************************************************************
*   This file contains functions and subroutines to calculate useful
*   kinematical variables 
c***************************************************************************
c   LIST OF ALL FUNCTIONS AND SUBROUTINES IN THIS FILE:
c
c     SUBROUTINE BOOSTN(P,R,Q)
c     SUBROUTINE boostl
c     SUBROUTINE legoy(p,pt,y,phi)
c     SUBROUTINE lego(p,pt,y,phi)
c     FUNCTION RSEPS(Y1,PHI1,Y2,PHI2)
c     FUNCTION rsep(p1,p2)
c     SUBROUTINE order2(x,y,low,high)
c     SUBROUTINE order3(x,y,z,low,mid,high)
c     FUNCTION rsepf(p1,p2)
c     SUBROUTINE switchmom 
c     FUNCTION srt(x)
c
c***************************************************************************
*   Last modified: 18.07.2006
****************************************************************************


C********** BOOSTN ****************************************************
C
      SUBROUTINE BOOSTN(P,R,Q)
C
C
C     The four vector P is assumed to be given in the rest frame of R,
C     which must be a timelike vector.
C     output Q is the vector P boosted to the frame in which R is given.
C                                              Compare Jackson, p.517
C                                              D. Zeppenfeld (28.6.1985)
C     New version with energy stored in zeroth component of four vector
C     arrays. Checked on May 27, 1988.
C
      implicit none

      REAL*8 P(0:3),R(0:3),Q(0:3),R_SQ
      REAL*8 BETA(3), X, Y, GAMMA
      INTEGER I
      X = 0D0
      Y = 0D0
      DO I = 1,3
         BETA(I) = R(I)/R(0)
         X = X + BETA(I)**2
         Y = Y + BETA(I)*P(I)
      ENDDO
      IF (X.GT.1D-16.AND.X.LT.(1D0-1D-12)) THEN
         GAMMA = 1D0/DSQRT(1D0-X)
         DO I = 1,3
            Q(I) = P(I)+BETA(I)*(Y*(GAMMA-1D0)/X + GAMMA*P(0))
         ENDDO
         Q(0) = GAMMA*(P(0) + Y)
      ELSE
         DO I = 0,3
            Q(I) = P(I)
         ENDDO
         R_SQ = R(0)**2-R(1)**2-R(2)**2-R(3)**2
         ! IF(X.GE.(1D0-1D-12)) 
         IF(R_SQ.le.0.d0) 
     *      WRITE(6,1000) R,R_SQ
      ENDIF
 1000 FORMAT (" The reference vector ",4G12.3," is not timelike."/
     1        " R**2 = ",G12.3)
      END
 

C********** BOOSTL ****************************************************
      SUBROUTINE boostl(y, p)
c
c     boost the vector p along the z-axis 
c
      implicit none
      double precision y, p(0:3), tmp(0:3)
      double precision coshy, sinhy

      sinhy = sinh(y)
      coshy = sqrt(1+sinhy**2)

      tmp(0) = coshy*p(0) - sinhy*p(3)
      tmp(3) = coshy*p(3) - sinhy*p(0)

      p(0) = tmp(0)
      p(3) = tmp(3)

      end


c-----------------------------------------------------------------------------
c	This subroutine converts a four-momentum into
c	pT, rapidity and phi variables.
c
      subroutine legoy(
     &              p,		!in:  particle four-momentum
     &              pt,		!out:  corresponding pT (transverse momentum)
     &              y,		!out:  corresponding y ( true rapidity)
     &              phi		!out:  corresponding phi (azimuthal angle)
     &               )
      implicit none
c
c declare input/output variables
c
      real*8 p(0:3), pt, y, phi
c
c declare local variables
c
      real*8 r_zero, r_half, r_pi, r_large
      parameter( r_zero=0.0d0, r_half=0.5d0 )
      parameter( r_pi=3.141592653589793238d0 )
      parameter( r_large=1.0d+3 )
c
      real*8 pt2, mperp2, arg
c
      pt2 = p(1)**2 + p(2)**2
      if ( pt2 .gt. r_zero ) then
         pt = sqrt( pt2 )
c
         mperp2 = p(0)**2 - p(3)**2
         if ( mperp2 .gt. r_zero ) then
            if ( p(3) .ge. r_zero ) then
               arg = (p(0) + p(3))**2 / mperp2
            else
               arg = mperp2 / (p(0) - p(3))**2
            end if
            y = r_half * log( arg )
         else
            if ( p(3) .ge. r_zero ) then
               y = r_large
            else
               y = -r_large
            end if
         end if
         phi = atan2( p(2), p(1) )
c
      else
         if ( p(3) .gt. r_zero ) then
            pt = r_zero
            y = r_large
            phi = r_zero
         else if ( p(3) .lt. r_zero ) then
            pt = r_zero
            y = -r_large
            phi = r_zero
         else
            pt = r_zero
            y = r_zero
            phi = r_zero
         end if
      end if
c
      return
      end
c
c------------------------------------------------------------------------------
c
c	This subroutine converts a four-momentum into
c	pT, rapidity and phi variables.
c
      subroutine lego(
     &              p,		!in:  particle four-momentum
     &              pt,		!out:  corresponding pT (transverse momentum)
     &              y,		!out:  corresponding eta (pseudorapidity)
     &              phi		!out:  corresponding phi (azimuthal angle)
     &               )
      implicit none
c
c declare input/output variables
c
      real*8 p(0:3), pt, y, phi
c
c declare local variables
c
      real*8 r_zero, r_half, r_pi, r_large
      parameter( r_zero=0.0d0, r_half=0.5d0 )
      parameter( r_pi=3.141592653589793238d0 )
      parameter( r_large=1.0d+3 )
c
      real*8 pt2, pabs, arg
c
      pt2 = p(1)**2 + p(2)**2
      if ( pt2 .gt. r_zero ) then
         pt = sqrt( pt2 )
         pabs = sqrt(pt2 + p(3)**2)
c
         if ( p(3) .ge. r_zero ) then
            arg = (pabs + p(3)) / pt
         else
            arg = pt / (pabs - p(3))
         end if
         y = log( arg )
         phi = atan2( p(2), p(1) )
c
      else
         if ( p(3) .gt. r_zero ) then
            pt = r_zero
            y = r_large
            phi = r_zero
         else if ( p(3) .lt. r_zero ) then
            pt = r_zero
            y = -r_large
            phi = r_zero
         else
            pt = r_zero
            y = r_zero
            phi = r_zero
         end if
      end if
c
      return
      end
C
C**********************  RSEPS  ****************************************
C
      FUNCTION RSEPS(Y1,PHI1,Y2,PHI2)
      IMPLICIT NONE
C
C   calculate the separation in the lego plot between two momenta with 
C   pseudorapidities Y1 and Y2  and azimuthal angles PHI1 and PHI2
C
      DOUBLE PRECISION RSEPS, Y1,PHI1,Y2,PHI2, DELPHI,PI
      PARAMETER (PI=3.14159 26535 89793238d0)
      DELPHI = ABS(PHI1-PHI2)
      IF (DELPHI.GT.PI) THEN
         DELPHI = 2*PI-DELPHI
      ENDIF
      IF (DELPHI.LT.0 .OR. DELPHI.GT.PI) THEN
         PRINT*," Problem in RSEPN. DELPHI = ",DELPHI
      ENDIF
      RSEPS = SQRT( (Y1-Y2)**2 + DELPHI**2 + 1D-30 )
      END
c
c----------------------------------------------------------------------------
c
c	This subroutine calculates the seperation in
c	R = ( deltaPhi^2 + deltaY^2 ) between two
c	four momenta.  Phi is the azimuthal angle, while
c	Y is the rapidity.
c
      real*8 function rsep(
     &                   p1,	!in:  four momentum 1
     &                   p2	!in:  four momentum 2
     &                    )
      implicit none
c
c declare input/output variables
c
      real*8 p1(0:3), p2(0:3)
c
c declare local variables
c
      real*8 r_two, r_pi
      parameter( r_two=2.0d0, r_pi=3.141592653589793238d0 )
c
      real*8 pt1, y1, phi1, pt2, y2, phi2
      real*8 dphi, dy
c
c determine "lego" plot parameters
c
      call lego( p1, pt1, y1, phi1 )
      call lego( p2, pt2, y2, phi2 )
c
      dy = y2 - y1
      dphi = abs( phi2 - phi1 )
c
      if ( dphi .gt. r_pi ) then
         rsep = sqrt( (r_two * r_pi - dphi)**2 + dy**2 )
      else
         rsep = sqrt( dphi**2 + dy**2 )
      end if
c
c done
c
      return
      end
c
c-------------------------------------------------------------------------------
c
c	This subroutine orders two given values
c
      subroutine order2(
     &              x,		!in:  input value 1
     &              y,		!in:  input value 2
     &              low,	!out:  lowest value
     &              high	!out:  highest value
     &                 )
      implicit none
c
c declare input/output variables
c
      real*8 x, y, low, high
c
      if ( x .gt. y ) then
         high = x
         low = y
      else
         high = y
         low = x
      end if
c
      return
      end
c
c-------------------------------------------------------------------------------
c
c	This subroutine orders three given values
c
      subroutine order3(
     &              x,		!in:  input value 1
     &              y,		!in:  input value 2
     &              z,		!in:  input value 3
     &              low,	!out:  lowest value
     &              mid,	!out:  middle value
     &              high	!out:  highest value
     &                 )
      implicit none
c
c declare input/output variables
c
      real*8 x, y, z, low, mid, high
c
c do simple tree based sort
c
      if ( (x .gt. y) .and. (x .gt. z) ) then
         if ( y .gt. z ) then
            high = x
            mid = y
            low = z
         else
            high = x
            mid = z
            low = y
         end if
      else if ( (y .gt. x) .and. (y .gt. z) ) then
         if ( x .gt. z ) then
            high = y
            mid = x
            low = z
         else
            high = y
            mid = z
            low = x
         end if
      else
         if ( x .gt. y ) then
            high = z
            mid = x
            low = y
         else
            high = z
            mid = y
            low = x
         end if
      end if
c
c done
c
      return
      end 


c
c----------------------------------------------------------------------------
c
      real*8 function rsepf(
     &                   p1,                  !in:  momentum 1
     &                   p2                   !in:  momentum 2
     &                    )
      implicit none
c
c declare input/output variables
c
      real*8 p1(0:7), p2(0:7)
c
c declare local variables
c
      real*8 r_pi, twopi
      parameter( r_pi=3.141592653589793238d0, twopi=2d0*r_pi )
c
      real*8 dphi, dy
c
      dy = p2(6) - p1(6)
      dphi = abs( p2(7) - p1(7) )
c
      if ( dphi .gt. r_pi ) then
         rsepf = sqrt( (twopi - dphi)**2 + dy**2 )
      else
         rsepf = sqrt( dphi**2 + dy**2 )
      end if
c
c done
c
      return
      end

*****************************************************
      subroutine switchmom (p1, p2)
c
c     switches momenta p1 <-> p2
c
      implicit none
      double precision p1(0:3), p2(0:3), switch(0:3,2)

      integer mu

      do mu = 0, 3
c         print*, p1(mu), p2(mu)
         switch(mu,1) = p1(mu)
         switch(mu,2) = p2(mu)
         p1(mu) = switch(mu,2)
         p2(mu) = switch(mu,1)
c         print*, p1(mu), p2(mu)
      enddo
c      print*, "---------------------"

      return
      end

*****************************************************
      real*8 function srt(x)
*****************************************************
      implicit none
      real*8 x
      if (x.gt.0) then
         srt = sqrt(x)
      elseif(x.lt.0) then
         srt = - sqrt(-x)
      else
         srt = 0
      endif
      end


********************************************************************************
********************************************************************************

*** The process pp -> W-W-jj is obtained by calculating pbar pbar -> W+W+jj 
*** In order to get the parity-odd distributions correct, the momenta must be
*** flipped.  This subroutine does that and re-calculates the rapidity.

      subroutine momswitch(p,v,jets,leptons,invisible)

      implicit none

#include "global.inc"

      real*8 p(0:3,max_p), v(0:3,max_v)
      real*8 jets(0:7,max_jets), leptons(0:8,max_v)
      real*8 invisible(0:8,max_v)

** do-loop variables for flipping momenta 
      integer p1, p2

** dummy variables for re-calculating rapidity.  
      double precision arg, pabs


      do p1 = 1, 3
         do p2 = 1, max_v
            v(p1,p2) = -v(p1,p2)
         end do
         do p2 = 1, max_p
            p(p1,p2) = -p(p1,p2)
         end do
      end do   
      do p1 = 1, 7
         if (p1 .ne. 5) then
            do p2 = 1, max_jets
               if (jets(p1,p2) .ne. 0d0) jets(p1,p2) = -jets(p1,p2)
            end do
         end if
      end do  
      do p1 = 1, 8
         if (p1 .ne. 5) then
            do p2 = 1, max_v
               if (leptons(p1,p2) .ne. 0d0) leptons(p1,p2) = 
     -              -leptons(p1,p2)
               if (invisible(p1,p2) .ne. 0D0) invisible(p1,p2) = 
     -              -invisible(p1,p2)
            end do
         end if
      end do
      do p2 = 1, max_v
         if (leptons(5,p2)**2 .gt. 0d0) then
            pabs = leptons(5,p2)**2 + leptons(3,p2)**2
            pabs = sqrt(pabs)
            if (leptons(3,p2) .ge. 0d0) then
               arg = (pabs + leptons(3,p2)) / leptons(5,p2)
            else
               arg = leptons(5,p2) / (pabs - leptons(3,p2))
            end if
            leptons(6,p2) = log(arg)
         end if
         
         if (invisible(5,p2)**2 .gt. 0d0) then
            pabs = invisible(5,p2)**2 + invisible(3,p2)**2
            pabs = sqrt(pabs)
            if (invisible(3,p2) .ge. 0d0) then
               arg = (pabs + invisible(3,p2)) / invisible(5,p2)
            else
               arg = invisible(5,p2) / (pabs - invisible(3,p2))
            end if
            invisible(6,p2) = log(arg)
         end if
      end do
      do p2 = 1, max_jets
         if (jets(5,p2)**2 .gt. 0d0) then
            pabs = jets(5,p2)**2 + jets(3,p2)**2
            pabs = sqrt(pabs)
            if (jets(3,p2) .ge. 0d0) then
               arg = (pabs + jets(3,p2)) / jets(5,p2)
            else
               arg = jets(5,p2) / (pabs - jets(3,p2))
            end if
            jets(6,p2) = log(arg)
         end if
      end do
      

      end


********************************************************************************
