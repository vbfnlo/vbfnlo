c
c This file contains routines that substitute the helas calls for the triple-vector-boson
c couplings with anomalous couplings.
c
c   LIST OF ALL SUBROUTINES IN THIS FILE:
c
c subroutine wwa_anomal3(wm, wp, a, vertex)
c subroutine jwpawm_anomal3(a, wm, jvv)
c subroutine jwmwpa_anomal3(wp, a, jvv)
c subroutine jaww_anomal3(wm, wp, jvv)
c
c subroutine wwz_anomal3(wm, wp, z, vertex)
c subroutine jwpzwm_anomal3(z, wm, jvv)
c subroutine jwmwpz_anomal3(wp, z, jvv)
c subroutine jzww_anomal3(wm, wp, jvv)
c
c subroutine zzh_anomal3(v1, v2, sc, vertex)
c subroutine jhzz_anomal3(v1,v2 , hvv)
c subroutine azh_anomal3(a, z, sc, vertex)
c subroutine jhaz_anomal3(a, z, hvv)
c subroutine aah_anomal3(v1, v2, sc, vertex)
c subroutine jhaa_anomal3(v1,v2 , hvv)
c subroutine wwh_anomal3(v1, v2, sc, vertex)
c subroutine jhww_anomal3(v1,v2 , hvv)


CCCCCCCCCCCCCCCC
CCCCC WWA CCCCCC
CCCCCCCCCCCCCCCC

      subroutine wwa_anomal3(wm, wp, a, vertex)
c
c This subroutine computes an amplitude of the three-point coupling of
c the gauge bosons w+ w- a with anomalous couplings.
c
c input:
c       complex wm(6)          : vector               flow-out W-
c       complex wp(6)          : vector               flow-out W+
c       complex a(6)           : vector               flow-out A
c output:
c       complex vertex         : amplitude               gamma(wm,wp,a)
c


      implicit none
      double complex wm(6),wp(6),a(6),vertex, dotcc, dotrc, epscccr, epscrrr, epsccrr
      double precision pwm(0:3),pwp(0:3),pa(0:3),dotrr
      external dotcc, dotrc, dotrr, epscccr, epscrrr, epsccrr

#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"

      pwm(0) = dble( wm(5))
      pwm(1) = dble( wm(6))
      pwm(2) = dimag(wm(6))
      pwm(3) = dimag(wm(5))
      pwp(0) = dble( wp(5))
      pwp(1) = dble( wp(6))
      pwp(2) = dimag(wp(6))
      pwp(3) = dimag(wp(5))
      pa(0) = dble( a(5))
      pa(1) = dble( a(6))
      pa(2) = dimag(a(6))
      pa(3) = dimag(a(5))

      vertex = - ( dotcc(wm,wp)*(dotrc(pwm,a )-dotrc(pwp,a ))
     &           + dotcc(wp,a )*(dotrc(pwp,wm)-dotrc(pa ,wm))
     &           + dotcc(a ,wm)*(dotrc(pa ,wp)-dotrc(pwm,wp))  )*gwwa

      if (fwww .ne. 0) then
        vertex = vertex + fwww * 3./2. * gw**3 * sqrt(sin2w) * ( dotrc(pa,wm) * dotrc(pwp,a) * dotrc(pwm,wp)
     &      - dotrc(pa,wp) * dotrc(pwp,wm) * dotrc(pwm,a)
     &      + dotrr(pa,pwm) * ( dotrc(pwp,wm) * dotcc(wp,a) - dotrc(pwp,a) * dotcc(wp,wm) )
     &      + dotrr(pwm,pwp) * ( dotrc(pa,wp) * dotcc(wm,a) - dotrc(pa,wm) * dotcc(wp,a) )
     &      + dotrr(pwp,pa) * ( dotrc(pwm,a) * dotcc(wp,wm) - dotrc(pwm,wp) * dotcc(wm,a) ) )
      endif

      if (fw .ne. 0) then
        vertex = vertex + fw * gw * wmass**2 / 2. * sqrt(sin2w) * ( dotrc(pa,wm)*dotcc(wp,a) - dotrc(pa,wp)*dotcc(wm,a) )
      endif

      if (fb .ne. 0) then
        vertex = vertex + fb * gw * wmass**2 / 2. * sqrt(sin2w) * ( dotrc(pa,wm)*dotcc(wp,a) - dotrc(pa,wp)*dotcc(wm,a) )
      endif

      if (fbw .ne. 0) then
        vertex = vertex + fbw * gw * wmass**2 * sqrt(sin2w) * ( dotrc(pa,wp)*dotcc(wm,a) - dotrc(pa,wm)*dotcc(wp,a) )
      endif

      if (fdw .ne. 0) then
        vertex = vertex - fdw * gw**3 * 2. * sqrt(sin2w) * (
     &    - dotcc(wm,a) * ( dotrc(pa,wp) * ( 2.*dotrr(pa,pwm) + dotrr(pa,pwp) )
     &    - dotrc(pwm,wp) * ( 2.*dotrr(pa,pwm) + dotrr(pwm,pwp) ) )
     &    + dotcc(wp,a) * ( dotrc(pa,wm) * ( dotrr(pa,pwm) + 2.*dotrr(pa,pwp) )
     &    - dotrc(pwp,wm) * ( 2.*dotrr(pa,pwp) + dotrr(pwm,pwp) ) )
     &    - dotcc(wp,wm) * ( dotrc(pwm,a) * ( dotrr(pa,pwm) + 2.*dotrr(pwm,pwp) )
     &    - dotrc(pwp,a) * ( dotrr(pa,pwp) + 2.*dotrr(pwm,pwp) ) )
     &    + dotrc(pa,wp)*dotrc(pa,wm)*dotrc(pwm,a)
     &    - dotrc(pwm,wp)*dotrc(pa,wm)*dotrc(pwm,a)
     &    - dotrc(pa,wp)*dotrc(pa,wm)*dotrc(pwp,a)
     &    + dotrc(pa,wp)*dotrc(pwp,wm)*dotrc(pwp,a)
     &    + dotrc(pwm,wp)*dotrc(pwp,wm)*dotrc(pwm,a)
     &    - dotrc(pwm,wp)*dotrc(pwp,wm)*dotrc(pwp,a) )
      endif

      if (fbtilde .ne. 0) then
        vertex = vertex + fbtilde * gw * wmass**2 / 2. * sqrt(sin2w) * epscccr(wp,wm,a,pa)
      endif

c zero due to momentum conservation
c      if (fwwtilde .ne. 0) then
c        vertex = vertex + fwwtilde * 2. * gw * wmass**2 * sqrt(sin2w) * (
c     &      epscccr(wp,wm,a,pa) + epscccr(wp,wm,a,pwm) + epscccr(wp,wm,a,pwp) )
c      endif

      if (fwwwtilde .ne. 0) then
        vertex = vertex + fwwwtilde * gw**3 / 2. * sqrt(sin2w) * ( dotcc(wp,wm)*epscrrr(a,pa,pwm,pwp)
     &      + dotcc(wp,a)*epscrrr(wm,pa,pwm,pwp) + dotcc(wm,a)*epscrrr(wp,pa,pwm,pwp)
     &      - dotrc(pa,wp)*epsccrr(wm,a,pwm,pwp) - dotrc(pwm,wp)*epsccrr(wm,a,pa,pwp)
     &      - dotrc(pa,wm)*epsccrr(wp,a,pwm,pwp) + dotrc(pwp,wm)*epsccrr(wp,a,pa,pwm)
     &      + dotrc(pwm,a)*epsccrr(wp,wm,pa,pwp) + dotrc(pwp,a)*epsccrr(wp,wm,pa,pwm)
     &      - dotrr(pa,pwm)*epscccr(wp,wm,a,pwp) - dotrr(pa,pwp)*epscccr(wp,wm,a,pwm)
     &      - dotrr(pwm,pwp)*epscccr(wp,wm,a,pa) )
      endif

      if (fdwtilde .ne. 0) then
        vertex = vertex + fdwtilde * 2. * gw**3 * sqrt(sin2w) * (
     &     - ( dotrr(pa,pwm) + dotrr(pa,pwp) ) * epscccr(wp,wm,a,pa)
     &     - ( dotrr(pa,pwm) + dotrr(pwm,pwp) ) * epscccr(wp,wm,a,pwm)
     &     - ( dotrr(pa,pwp) + dotrr(pwm,pwp) ) * epscccr(wp,wm,a,pwp)
     &     + ( dotrc(pa,wp) - dotrc(pwm,wp) ) * epsccrr(wm,a,pa,pwm)
     &     + ( dotrc(pwp,wm) - dotrc(pa,wm) ) * epsccrr(wp,a,pa,pwp)
     &     + ( dotrc(pwm,a) - dotrc(pwp,a) ) * epsccrr(wp,wm,pwm,pwp) )
      endif

      if (fwtilde .ne. 0) then
        vertex = vertex + fwtilde / 2. * gw * wmass**2 * sqrt(sin2w) * epscccr(wp,wm,a,pa)
      endif

c Owtilde as linear combination. Same result but much slower.
c      if (fwtilde .ne. 0) then
c        vertex = vertex - fwtilde / 2. * gw * wmass**2 * sqrt(sin2w) *
c     &      ( epscccr(wp,wm,a,pa) + 2. * ( epscccr(wp,wm,a,pwm) + epscccr(wp,wm,a,pwp) ) )
c      endif

      if (fbwtilde .ne. 0) then
        vertex = vertex - fbwtilde * gw * wmass**2 * sqrt(sin2w) * epscccr(wp,wm,a,pa)
      endif

      return
      end


      subroutine jwpawm_anomal3(a, wm, jvv)
c
c This subroutine computes an off-shell vector current attached to a W+
c from the three-point gauge boson coupling w+ w- a with anomalous
c couplings.  The vector propagator is given in Feynman gauge for a
c massless vector and in unitary gauge for a massive vector.
c
c input:
c       complex a(6)           : first  vector                        v1
c       complex wm(6)          : second vector                        wm
c output:
c       complex jvv(6)         : vector current            j^mu(v:a,wm)
c
c where all the bosons are defined by the flowing-OUT quantum number.

      implicit none
      double complex a(6), wm(6), jvv(6), j12(4), dg, dotcc, dotrc
      double complex epsccrmu, epscrrmu, epsccrr, epscrrr
      double complex cm2        ! mass**2- I Gamma mass (Fabio)
      double precision pa(4),pwm(4),q(4),s
      double precision vm2, dotrr, epsrrrmu
      integer i
      external dotcc, dotrc, dotrr, epsccrmu, epsrrrmu, epscrrmu, epscrrr, epsccrr

#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"

      jvv(5) = a(5)+wm(5)
      jvv(6) = a(6)+wm(6)

      pa(1) =  dble( a(5))
      pa(2) =  dble( a(6))
      pa(3) =  dimag(a(6))
      pa(4) =  dimag(a(5))
      pwm(1) =  dble( wm(5))
      pwm(2) =  dble( wm(6))
      pwm(3) =  dimag(wm(6))
      pwm(4) =  dimag(wm(5))
      q(1)  = -dble( jvv(5))
      q(2)  = -dble( jvv(6))
      q(3)  = -dimag(jvv(6))
      q(4)  = -dimag(jvv(5))
      s = dotrr(q,q)
      vm2 = wmass**2

c     Vertex
      do i=1,4
            j12(i) = -( (pa(i)-pwm(i))*dotcc(a,wm)
     &                      + (dotrc(pwm,a)-dotrc(q,a))*wm(i)
     &                      + (dotrc(q,wm)-dotrc(pa,wm))*a(i) )* gwwa

        if (fwww .ne. 0) then
          j12(i) = j12(i) + fwww * 3./2. * gw**3 * sqrt(sin2w) *
     &        ( dotrc(pa,wm) * dotrc(q,a) * pwm(i) - pa(i) * dotrc(q,wm) * dotrc(pwm,a)
     &        + dotrr(pa,pwm) * ( dotrc(q,wm) * a(i) - dotrc(q,a) * wm(i) )
     &        + dotrr(pwm,q) * ( pa(i) * dotcc(wm,a) - dotrc(pa,wm) * a(i) )
     &        + dotrr(q,pa) * ( dotrc(pwm,a) * wm(i) - pwm(i) * dotcc(wm,a) ) )
        endif

        if (fw .ne. 0) then
          j12(i) = j12(i) + fw * gw * wmass**2 / 2. * sqrt(sin2w) *
     &                           ( dotrc(pa,wm)*a(i) - pa(i)*dotcc(wm,a) )
        endif

        if (fb .ne. 0) then
          j12(i) = j12(i) + fb * gw * wmass**2 / 2. * sqrt(sin2w) *
     &                           ( dotrc(pa,wm)*a(i) - pa(i)*dotcc(wm,a) )
        endif

        if (fbw .ne. 0) then
          j12(i) = j12(i) + fbw * gw * wmass**2 * sqrt(sin2w) * ( pa(i)*dotcc(wm,a) - dotrc(pa,wm)*a(i) )
        endif

        if (fdw .ne. 0) then
          j12(i) = j12(i) - fdw * gw**3 * 2. * sqrt(sin2w) * (
     &      - dotcc(wm,a) * ( pa(i) * ( 2.*dotrr(pa,pwm) + dotrr(pa,q) ) - pwm(i) * ( 2.*dotrr(pa,pwm) + dotrr(pwm,q) ) )
     &      + a(i) * ( dotrc(pa,wm) * ( dotrr(pa,pwm) + 2.*dotrr(pa,q) ) - dotrc(q,wm) * ( 2.*dotrr(pa,q) + dotrr(pwm,q) ) )
     &      - wm(i) * ( dotrc(pwm,a) * ( dotrr(pa,pwm) + 2.*dotrr(pwm,q) ) - dotrc(q,a) * ( dotrr(pa,q) + 2.*dotrr(pwm,q) ) )
     &      + pa(i)*dotrc(pa,wm)*dotrc(pwm,a) - pwm(i)*dotrc(pa,wm)*dotrc(pwm,a) - pa(i)*dotrc(pa,wm)*dotrc(q,a)
     &      + pa(i)*dotrc(q,wm)*dotrc(q,a) + pwm(i)*dotrc(q,wm)*dotrc(pwm,a) - pwm(i)*dotrc(q,wm)*dotrc(q,a) )
        endif

        if (fbtilde .ne. 0) then
          j12(i) = j12(i) - fbtilde * gw * wmass**2 / 2. * sqrt(sin2w) * epsccrmu(wm,a,pa,i-1)
        endif

c zero due to momentum conservation
c        if (fwwtilde .ne. 0) then
c          j12(i) = j12(i) - fwwtilde * 2. * gw * wmass**2 * sqrt(sin2w) * (
c     &      epsccrmu(wm,a,pa,i-1) + epsccrmu(wm,a,pwm,i-1) + epsccrmu(wm,a,q,i-1) )
c        endif

        if (fwwwtilde .ne. 0) then
          j12(i) = j12(i) + fwwwtilde * gw**3 / 2. * sqrt(sin2w) * ( wm(i)*epscrrr(a,pa,pwm,q)
     &        + a(i)*epscrrr(wm,pa,pwm,q) - dotcc(wm,a)*epsrrrmu(pa,pwm,q,i-1)
     &        - pa(i)*epsccrr(wm,a,pwm,q) - pwm(i)*epsccrr(wm,a,pa,q)
     &        + dotrc(pa,wm)*epscrrmu(a,pwm,q,i-1) - dotrc(q,wm)*epscrrmu(a,pa,pwm,i-1)
     &        - dotrc(pwm,a)*epscrrmu(wm,pa,q,i-1) - dotrc(q,a)*epscrrmu(wm,pa,pwm,i-1)
     &        + dotrr(pa,pwm)*epsccrmu(wm,a,q,i-1) + dotrr(pa,q)*epsccrmu(wm,a,pwm,i-1)
     &        + dotrr(pwm,q)*epsccrmu(wm,a,pa,i-1) )
        endif

        if (fdwtilde .ne. 0) then
          j12(i) = j12(i) + fdwtilde * 2. * gw**3 * sqrt(sin2w) * (
     &     + ( dotrr(pa,pwm) + dotrr(pa,q) ) * epsccrmu(wm,a,pa,i-1)
     &     + ( dotrr(pa,pwm) + dotrr(pwm,q) ) * epsccrmu(wm,a,pwm,i-1)
     &     + ( dotrr(pa,q) + dotrr(pwm,q) ) * epsccrmu(wm,a,q,i-1)
     &     + ( pa(i) - pwm(i) ) * epsccrr(wm,a,pa,pwm)
     &     - ( dotrc(q,wm) - dotrc(pa,wm) ) * epscrrmu(a,pa,q,i-1)
     &     - ( dotrc(pwm,a) - dotrc(q,a) ) * epscrrmu(wm,pwm,q,i-1) )
        endif

        if (fwtilde .ne. 0) then
          j12(i) = j12(i) - fwtilde / 2. * gw * wmass**2 * sqrt(sin2w) * epsccrmu(wm,a,pa,i-1)
        endif

c        if (fwtilde .ne. 0) then
c          j12(i) = j12(i) + fwtilde / 2. * gw * wmass**2 * sqrt(sin2w) *
c     &      ( epsccrmu(wm,a,pa,i-1) + 2. * ( epsccrmu(wm,a,pwm,i-1) + epsccrmu(wm,a,q,i-1) ) )
c        endif

        if (fbwtilde .ne. 0) then
          j12(i) = j12(i) + fbwtilde * gw * wmass**2 * sqrt(sin2w) * epsccrmu(wm,a,pa,i-1)
        endif

      enddo


c     Propagator
      cm2=dcmplx( vm2, -wmass*wwidth )       !Fabio's implementation of the fixed width
      dg = 1.d0 / dcmplx( s-vm2, wmass*wwidth )       !factor of i^2 here
      do i=1,4
         jvv(i) = dg * ( j12(i) - q(i)*dotrc(q,j12)/cm2 )
      enddo

      return
      end


      subroutine jwmwpa_anomal3(wp, a, jvv)
c
c This subroutine computes an off-shell vector current attached to a W-
c from the three-point gauge boson coupling w+ w- a with anomalous
c couplings.  The vector propagator is given in Feynman gauge for a
c massless vector and in unitary gauge for a massive vector.
c
c input:
c       complex wp(6)          : first  vector                        wp
c       complex a(6)           : second vector                        a
c output:
c       complex jvv(6)         : vector current            j^mu(v:wp,a)
c
c where all the bosons are defined by the flowing-OUT quantum number.

      implicit none
      double complex wp(6), a(6), jvv(6), j12(4), dg, dotcc, dotrc
      double complex epsccrmu, epscrrmu, epscrrr, epsccrr
      double complex cm2        ! mass**2- I Gamma mass (Fabio)
      double precision pwp(4),pa(4),q(4),s
      double precision vm2, dotrr, epsrrrmu
      integer i
      external dotcc, dotrc, dotrr, epsccrmu, epsrrrmu, epscrrmu, epscrrr, epsccrr

#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"

      jvv(5) = wp(5)+a(5)
      jvv(6) = wp(6)+a(6)

      pwp(1) =  dble( wp(5))
      pwp(2) =  dble( wp(6))
      pwp(3) =  dimag(wp(6))
      pwp(4) =  dimag(wp(5))
      pa(1) =  dble( a(5))
      pa(2) =  dble( a(6))
      pa(3) =  dimag(a(6))
      pa(4) =  dimag(a(5))
      q(1)  = -dble( jvv(5))
      q(2)  = -dble( jvv(6))
      q(3)  = -dimag(jvv(6))
      q(4)  = -dimag(jvv(5))
      s = dotrr(q,q)
      vm2 = wmass**2

c     Vertex
      do i=1,4
            j12(i) = -gwwa * ( (pwp(i)-pa(i))*dotcc(wp,a)
     &                                     + (dotrc(pa,wp)-dotrc(q,wp))*a(i)
     &                                     + (dotrc(q,a)-dotrc(pwp,a))*wp(i) )

        if (fwww .ne. 0) then
          j12(i) = j12(i) + fwww * 3./2. * gw**3 * sqrt(sin2w) *
     &        ( pa(i) * dotrc(pwp,a) * dotrc(q,wp) - dotrc(pa,wp) * pwp(i) * dotrc(q,a)
     &        + dotrr(pa,q) * ( pwp(i) * dotcc(wp,a) - dotrc(pwp,a) * wp(i) )
     &        + dotrr(q,pwp) * ( dotrc(pa,wp) * a(i) - pa(i) * dotcc(wp,a) )
     &        + dotrr(pwp,pa) * ( dotrc(q,a) * wp(i) - dotrc(q,wp) * a(i) ) )
        endif

        if (fw .ne. 0) then
          j12(i) = j12(i) + fw * gw * wmass**2 / 2. * sqrt(sin2w) *
     &        ( pa(i)*dotcc(wp,a) - dotrc(pa,wp)*a(i) )
        endif

        if (fb .ne. 0) then
          j12(i) = j12(i) + fb * gw * wmass**2 / 2. * sqrt(sin2w) *
     &        ( pa(i)*dotcc(wp,a) - dotrc(pa,wp)*a(i) )
        endif

        if (fbw .ne. 0) then
          j12(i) = j12(i) + fbw * gw * wmass**2 * sqrt(sin2w) * ( dotrc(pa,wp)*a(i) - pa(i)*dotcc(wp,a) )
        endif

        if (fdw .ne. 0) then
          j12(i) = j12(i) - fdw * gw**3 * 2. * sqrt(sin2w) * (
     &      - a(i) * ( dotrc(pa,wp) * ( 2.*dotrr(pa,q) + dotrr(pa,pwp) ) - dotrc(q,wp) * ( 2.*dotrr(pa,q) + dotrr(q,pwp) ) )
     &      + dotcc(wp,a) * ( pa(i) * ( dotrr(pa,q) + 2.*dotrr(pa,pwp) ) - pwp(i) * ( 2.*dotrr(pa,pwp) + dotrr(q,pwp) ) )
     &      - wp(i) * ( dotrc(q,a) * ( dotrr(pa,q) + 2.*dotrr(q,pwp) ) - dotrc(pwp,a) * ( dotrr(pa,pwp) + 2.*dotrr(q,pwp) ) )
     &      + dotrc(pa,wp)*pa(i)*dotrc(q,a) - dotrc(q,wp)*pa(i)*dotrc(q,a) - dotrc(pa,wp)*pa(i)*dotrc(pwp,a)
     &      + dotrc(pa,wp)*pwp(i)*dotrc(pwp,a) + dotrc(q,wp)*pwp(i)*dotrc(q,a) - dotrc(q,wp)*pwp(i)*dotrc(pwp,a) )
        endif

        if (fbtilde .ne. 0) then
          j12(i) = j12(i) + fbtilde * gw * wmass**2 / 2. * sqrt(sin2w) * epsccrmu(wp,a,pa,i-1)
        endif

c zero due to momentum conservation
c        if (fwwtilde .ne. 0) then
c          j12(i) = j12(i) + fwwtilde * 2. * gw * wmass**2 * sqrt(sin2w) * (
c     &      epsccrmu(wp,a,pa,i-1) + epsccrmu(wp,a,q,i-1) + epsccrmu(wp,a,pwp,i-1) )
c        endif

        if (fwwwtilde .ne. 0) then
          j12(i) = j12(i) + fwwwtilde * gw**3 / 2. * sqrt(sin2w) * ( wp(i)*epscrrr(a,pa,q,pwp)
     &      - dotcc(wp,a)*epsrrrmu(pa,q,pwp,i-1) + a(i)*epscrrr(wp,pa,q,pwp)
     &      + dotrc(pa,wp)*epscrrmu(a,q,pwp,i-1) + dotrc(q,wp)*epscrrmu(a,pa,pwp,i-1)
     &      - pa(i)*epsccrr(wp,a,q,pwp) + pwp(i)*epsccrr(wp,a,pa,q)
     &      + dotrc(q,a)*epscrrmu(wp,pa,pwp,i-1) + dotrc(pwp,a)*epscrrmu(wp,pa,q,i-1)
     &      - dotrr(pa,q)*epsccrmu(wp,a,pwp,i-1) - dotrr(pa,pwp)*epsccrmu(wp,a,q,i-1)
     &      - dotrr(q,pwp)*epsccrmu(wp,a,pa,i-1) )
        endif

        if (fdwtilde .ne. 0) then
          j12(i) = j12(i) + fdwtilde * 2. * gw**3 * sqrt(sin2w) * (
     &     - ( dotrr(pa,q) + dotrr(pa,pwp) ) * epsccrmu(wp,a,pa,i-1)
     &     - ( dotrr(pa,q) + dotrr(q,pwp) ) * epsccrmu(wp,a,q,i-1)
     &     - ( dotrr(pa,pwp) + dotrr(q,pwp) ) * epsccrmu(wp,a,pwp,i-1)
     &     - ( dotrc(pa,wp) - dotrc(q,wp) ) * epscrrmu(a,pa,q,i-1)
     &     + ( pwp(i) - pa(i) ) * epsccrr(wp,a,pa,pwp)
     &     + ( dotrc(q,a) - dotrc(pwp,a) ) * epscrrmu(wp,q,pwp,i-1) )
        endif

        if (fwtilde .ne. 0) then
          j12(i) = j12(i) + fwtilde / 2. * gw * wmass**2 * sqrt(sin2w) * epsccrmu(wp,a,pa,i-1)
        endif

c        if (fwtilde .ne. 0) then
c          j12(i) = j12(i) - fwtilde / 2. * gw * wmass**2 * sqrt(sin2w) *
c     &      ( epsccrmu(wp,a,pa,i-1) + 2. * ( epsccrmu(wp,a,q,i-1) + epsccrmu(wp,a,pwp,i-1) ) )
c        endif

        if (fbwtilde .ne. 0) then
          j12(i) = j12(i) - fbwtilde * gw * wmass**2 * sqrt(sin2w) * epsccrmu(wp,a,pa,i-1)
        endif

      enddo


c     Propagator
      cm2=dcmplx( vm2, -wmass*wwidth )       !Fabio's implementation of the fixed width
      dg = 1.d0 / dcmplx( s-vm2, wmass*wwidth )
      do i=1,4
         jvv(i) = dg * ( j12(i) - q(i)*dotrc(q,j12)/cm2 )
      enddo

      return
      end


      subroutine jaww_anomal3(wm, wp, jvv)
c
c This subroutine computes an off-shell vector current attached to an A
c from the three-point gauge boson coupling w+ w- a with anomalous
c couplings.  The vector propagator is given in Feynman gauge for a
c massless vector and in unitary gauge for a massive vector.
c
c input:
c       complex wm(6)          : vector               flow-out W-
c       complex wp(6)          : vector               flow-out W+
c where all the bosons are defined by the flowing-OUT quantum number.
c
c output:
c       complex jvv(6)         : vector current A          j^mu(v:wp,wp)
c
c where all the bosons are defined by the flowing-OUT quantum number.

      implicit none
      double complex wm(6), wp(6), jvv(6), j12(4), dg, dotcc, dotrc
      double complex epsccrmu, epscrrmu, epscrrr, epsccrr
      double precision pwm(4),pwp(4),q(4),s
      double precision dotrr, epsrrrmu
      integer i
      external dotcc, dotrc, dotrr, epsccrmu, epscrrmu, epscrrr, epsccrr, epsrrrmu

#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"

      jvv(5) = wm(5)+wp(5)
      jvv(6) = wm(6)+wp(6)

      pwm(1) =  dble( wm(5))
      pwm(2) =  dble( wm(6))
      pwm(3) =  dimag(wm(6))
      pwm(4) =  dimag(wm(5))
      pwp(1) =  dble( wp(5))
      pwp(2) =  dble( wp(6))
      pwp(3) =  dimag(wp(6))
      pwp(4) =  dimag(wp(5))
      q(1)  = -dble( jvv(5))
      q(2)  = -dble( jvv(6))
      q(3)  = -dimag(jvv(6))
      q(4)  = -dimag(jvv(5))
      s = dotrr(q,q)


c     Vertex
      do i=1,4
            j12(i) = -gwwa * (  dotcc(wm,wp)   *pwm(i)
     &                        - dotcc(wm,wp)   *pwp(i)
     &                        + dotrc(pwp,wm)  *wp(i)
     &                        - dotrc(q,wm)    *wp(i)
     &                        + dotrc(q,wp)    *wm(i)
     &                        - dotrc(pwm,wp)  *wm(i) )

        if (fwww .ne. 0) then
          j12(i) = j12(i) + fwww * 3./2. * gw**3 * sqrt(sin2w) *
     &                       ( dotrc(q,wm)    * dotrc(pwm,wp) * pwp(i)
     &                       - dotrr(q,pwm)   * dotcc(wm,wp)  * pwp(i)
     &                       + dotrr(pwp,q)   * dotcc(wm,wp)  * pwm(i)
     &                       - dotrc(q,wp)    * dotrc(pwp,wm) * pwm(i)
     &                       + dotrr(q,pwm)   * dotrc(pwp,wm) * wp(i)
     &                       - dotrr(pwm,pwp) * dotrc(q,wm)   * wp(i)
     &                       + dotrr(pwm,pwp) * dotrc(q,wp)   * wm(i)
     &                       - dotrr(pwp,q)   * dotrc(pwm,wp) * wm(i)  )

        endif

        if (fw .ne. 0) then
          j12(i) = j12(i) + fw * gw * wmass**2 / 2. * sqrt(sin2w) *
     &        ( dotrc(q,wm)*wp(i) - dotrc(q,wp)*wm(i) )
        endif

        if (fb .ne. 0) then
          j12(i) = j12(i) + fb * gw * wmass**2 / 2. * sqrt(sin2w) *
     &        ( dotrc(q,wm)*wp(i) - dotrc(q,wp)*wm(i) )
        endif

        if (fbw .ne. 0) then
          j12(i) = j12(i) + fbw * gw * wmass**2 * sqrt(sin2w) * ( dotrc(q,wp)*wm(i) - dotrc(q,wm)*wp(i) )
        endif

        if (fdw .ne. 0) then
          j12(i) = j12(i) - fdw * gw**3 * 2. * sqrt(sin2w) * (
     &      - wm(i) * ( dotrc(q,wp) * ( 2.*dotrr(q,pwm) + dotrr(q,pwp) ) - dotrc(pwm,wp) * ( 2.*dotrr(q,pwm) + dotrr(pwm,pwp) ) )
     &      + wp(i) * ( dotrc(q,wm) * ( dotrr(q,pwm) + 2.*dotrr(q,pwp) ) - dotrc(pwp,wm) * ( 2.*dotrr(q,pwp) + dotrr(pwm,pwp) ) )
     &      - dotcc(wp,wm) * ( pwm(i) * ( dotrr(q,pwm) + 2.*dotrr(pwm,pwp) ) - pwp(i) * ( dotrr(q,pwp) + 2.*dotrr(pwm,pwp) ) )
     &      + dotrc(q,wp)*dotrc(q,wm)*pwm(i) - dotrc(pwm,wp)*dotrc(q,wm)*pwm(i) - dotrc(q,wp)*dotrc(q,wm)*pwp(i)
     &      + dotrc(q,wp)*dotrc(pwp,wm)*pwp(i) + dotrc(pwm,wp)*dotrc(pwp,wm)*pwm(i) - dotrc(pwm,wp)*dotrc(pwp,wm)*pwp(i) )
        endif

        if (fbtilde .ne. 0) then
          j12(i) = j12(i) - fbtilde * gw * wmass**2 / 2. * sqrt(sin2w) * epsccrmu(wp,wm,q,i-1)
        endif

c zero due to momentum conservation
c        if (fwwtilde .ne. 0) then
c          j12(i) = j12(i) - fwwtilde * 2. * gw * wmass**2 * sqrt(sin2w) * (
c     &      epsccrmu(wp,wm,q,i-1) + epsccrmu(wp,wm,pwm,i-1) + epsccrmu(wp,wm,pwp,i-1) )
c        endif

        if (fwwwtilde .ne. 0) then
          j12(i) = j12(i) + fwwwtilde * gw**3 / 2. * sqrt(sin2w) * ( - dotcc(wp,wm)*epsrrrmu(q,pwm,pwp,i-1)
     &      + wp(i)*epscrrr(wm,q,pwm,pwp) + wm(i)*epscrrr(wp,q,pwm,pwp)
     &      - dotrc(q,wp)*epscrrmu(wm,pwm,pwp,i-1) - dotrc(pwm,wp)*epscrrmu(wm,q,pwp,i-1)
     &      - dotrc(q,wm)*epscrrmu(wp,pwm,pwp,i-1) + dotrc(pwp,wm)*epscrrmu(wp,q,pwm,i-1)
     &      + pwm(i)*epsccrr(wp,wm,q,pwp) + pwp(i)*epsccrr(wp,wm,q,pwm)
     &      + dotrr(q,pwm)*epsccrmu(wp,wm,pwp,i-1) + dotrr(q,pwp)*epsccrmu(wp,wm,pwm,i-1)
     &      + dotrr(pwm,pwp)*epsccrmu(wp,wm,q,i-1) )
        endif

        if (fdwtilde .ne. 0) then
          j12(i) = j12(i) + fdwtilde * 2. * gw**3 * sqrt(sin2w) * (
     &     + ( dotrr(q,pwm) + dotrr(q,pwp) ) * epsccrmu(wp,wm,q,i-1)
     &     + ( dotrr(q,pwm) + dotrr(pwm,pwp) ) * epsccrmu(wp,wm,pwm,i-1)
     &     + ( dotrr(q,pwp) + dotrr(pwm,pwp) ) * epsccrmu(wp,wm,pwp,i-1)
     &     + ( dotrc(q,wp) - dotrc(pwm,wp) ) * epscrrmu(wm,q,pwm,i-1)
     &     + ( dotrc(pwp,wm) - dotrc(q,wm) ) * epscrrmu(wp,q,pwp,i-1)
     &     + ( pwm(i) - pwp(i) ) * epsccrr(wp,wm,pwm,pwp) )
        endif

        if (fwtilde .ne. 0) then
          j12(i) = j12(i) - fwtilde / 2. * gw * wmass**2 * sqrt(sin2w) * epsccrmu(wp,wm,q,i-1)
        endif

c        if (fwtilde .ne. 0) then
c          j12(i) = j12(i) + fwtilde / 2. * gw * wmass**2 * sqrt(sin2w) *
c     &      ( epsccrmu(wp,wm,q,i-1) + 2. * ( epsccrmu(wp,wm,pwm,i-1) + epsccrmu(wp,wm,pwp,i-1) ) )
c        endif

        if (fbwtilde .ne. 0) then
          j12(i) = j12(i) + fbwtilde * gw * wmass**2 * sqrt(sin2w) * epsccrmu(wp,wm,q,i-1)
        endif

      enddo

c     Propagator
      dg = 1.d0 / s
      do i=1,4
         jvv(i) = dg*j12(i)
      enddo

      return
      end


CCCCCCCCCCCCCCCCC
CCCCCC WWZ CCCCCC
CCCCCCCCCCCCCCCCC

      subroutine wwz_anomal3(wm, wp, z, vertex)
c
c This subroutine computes an amplitude of the three-point coupling of
c the gauge bosons w+ w- z with anomalous couplings.
c
c input:
c       complex wm(6)          : vector               flow-out W-
c       complex wp(6)          : vector               flow-out W+
c       complex z(6)           : vector               flow-out Z
c output:
c       complex vertex         : amplitude               gamma(wm,wp,z)
c
      implicit none
      double complex wm(6),wp(6),z(6),vertex
      double complex dotcc, dotrc, epscccr, epscrrr, epsccrr
      double precision pwm(0:3),pwp(0:3),pz(0:3), dotrr, cos2w
      external dotcc, dotrc, dotrr, epscccr, epscrrr, epsccrr

#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"

      cos2w=1.0d0 - sin2w

      pwm(0) = dble( wm(5))
      pwm(1) = dble( wm(6))
      pwm(2) = dimag(wm(6))
      pwm(3) = dimag(wm(5))
      pwp(0) = dble( wp(5))
      pwp(1) = dble( wp(6))
      pwp(2) = dimag(wp(6))
      pwp(3) = dimag(wp(5))
      pz(0) = dble( z(5))
      pz(1) = dble( z(6))
      pz(2) = dimag(z(6))
      pz(3) = dimag(z(5))

      vertex = - ( dotcc(wm,wp)*(dotrc(pwm,z )-dotrc(pwp,z ))
     &           + dotcc(wp,z )*(dotrc(pwp,wm)-dotrc(pz ,wm))
     &           + dotcc(z ,wm)*(dotrc(pz ,wp)-dotrc(pwm,wp))  )*gwwz

      if (fwww .ne. 0) then
        vertex = vertex + fwww * 3./2. * gw**3 * sqrt(cos2w) * ( dotrc(pz,wm) * dotrc(pwp,z) * dotrc(pwm,wp)
     &      - dotrc(pz,wp) * dotrc(pwp,wm) * dotrc(pwm,z)
     &      + dotrr(pz,pwm) * ( dotrc(pwp,wm) * dotcc(wp,z) - dotrc(pwp,z) * dotcc(wp,wm) )
     &      + dotrr(pwm,pwp) * ( dotrc(pz,wp) * dotcc(wm,z) - dotrc(pz,wm) * dotcc(wp,z) )
     &      + dotrr(pwp,pz) * ( dotrc(pwm,z) * dotcc(wp,wm) - dotrc(pwm,wp) * dotcc(wm,z) ) )
      endif

      if (fw .ne. 0) then
        vertex = vertex + fw * gw * wmass**2 / ( 2. * sqrt(cos2w) ) * ( (dotrc(pwm,wp)-dotrc(pz,wp)) * dotcc(wm,z)
     &      + (dotrc(pz,wm)-dotrc(pwp,wm)) * dotcc(wp,z) + (dotrc(pwp,z)-dotrc(pwm,z)) * dotcc(wp,wm)
     &      - sin2w * ( dotrc(pz,wm)*dotcc(wp,z) - dotrc(pz,wp)*dotcc(wm,z) ) )
      endif

      if (fb .ne. 0) then
        vertex = vertex - fb * gw * wmass**2 / ( 2. * sqrt(cos2w) ) *
     &       sin2w * ( dotrc(pz,wm)*dotcc(wp,z) - dotrc(pz,wp)*dotcc(wm,z) )
      endif

      if (fbw .ne. 0) then
        vertex = vertex - fbw * gw * wmass**2 * sin2w / sqrt(cos2w) *
     &      ( dotrc(pz,wp)*dotcc(wm,z) - dotrc(pz,wm)*dotcc(wp,z) )
      endif

      if (fdw .ne. 0) then
        vertex = vertex - fdw * gw**3 * 2. * sqrt(cos2w) * (
     &    - dotcc(wm,z) * ( dotrc(pz,wp) * ( 2.*dotrr(pz,pwm) + dotrr(pz,pwp) )
     &    - dotrc(pwm,wp) * ( 2.*dotrr(pz,pwm) + dotrr(pwm,pwp) ) )
     &    + dotcc(wp,z) * ( dotrc(pz,wm) * ( dotrr(pz,pwm) + 2.*dotrr(pz,pwp) )
     &    - dotrc(pwp,wm) * ( 2.*dotrr(pz,pwp) + dotrr(pwm,pwp) ) )
     &    - dotcc(wp,wm) * ( dotrc(pwm,z) * ( dotrr(pz,pwm) + 2.*dotrr(pwm,pwp) )
     &    - dotrc(pwp,z) * ( dotrr(pz,pwp) + 2.*dotrr(pwm,pwp) ) )
     &    + dotrc(pz,wp)*dotrc(pz,wm)*dotrc(pwm,z)
     &    - dotrc(pwm,wp)*dotrc(pz,wm)*dotrc(pwm,z)
     &    - dotrc(pz,wp)*dotrc(pz,wm)*dotrc(pwp,z)
     &    + dotrc(pz,wp)*dotrc(pwp,wm)*dotrc(pwp,z)
     &    + dotrc(pwm,wp)*dotrc(pwp,wm)*dotrc(pwm,z)
     &    - dotrc(pwm,wp)*dotrc(pwp,wm)*dotrc(pwp,z) )
      endif

      if (fbtilde .ne. 0) then
        vertex = vertex - fbtilde * gw * wmass**2 / 2. * sin2w / sqrt(cos2w) * epscccr(wp,wm,z,pz)
      endif

c zero due to momentum conservation
c      if (fwwtilde .ne. 0) then
c        vertex = vertex + fwwtilde * 2. * gw * wmass**2 * sqrt(cos2w) * (
c     &      epscccr(wp,wm,z,pwm) + epscccr(wp,wm,z,pwp) + epscccr(wp,wm,z,pz) )
c      endif

      if (fwwwtilde .ne. 0) then
        vertex = vertex + fwwwtilde * gw**3 / 2. * sqrt(cos2w) * (dotcc(wp,wm)*epscrrr(z,pwm,pwp,pz)
     &      + dotcc(wp,z)*epscrrr(wm,pwm,pwp,pz) + dotcc(wm,z)*epscrrr(wp,pwm,pwp,pz)
     &      + dotrc(pwm,wp)*epsccrr(wm,z,pwp,pz) - dotrc(pz,wp)*epsccrr(wm,z,pwm,pwp)
     &      - dotrc(pwp,wm)*epsccrr(wp,z,pwm,pz) - dotrc(pz,wm)*epsccrr(wp,z,pwm,pwp)
     &      - dotrc(pwm,z)*epsccrr(wp,wm,pwp,pz) - dotrc(pwp,z)*epsccrr(wp,wm,pwm,pz)
     &      - dotrr(pwm,pwp)*epscccr(wp,wm,z,pz) - dotrr(pwm,pz)*epscccr(wp,wm,z,pwp)
     &      - dotrr(pwp,pz)*epscccr(wp,wm,z,pwm) )
      endif

      if (fdwtilde .ne. 0) then
        vertex = vertex + fdwtilde * 2. * gw**3 * sqrt(cos2w) * (
     &     - ( dotrr(pz,pwm) + dotrr(pz,pwp) ) * epscccr(wp,wm,z,pz)
     &     - ( dotrr(pz,pwm) + dotrr(pwm,pwp) ) * epscccr(wp,wm,z,pwm)
     &     - ( dotrr(pz,pwp) + dotrr(pwm,pwp) ) * epscccr(wp,wm,z,pwp)
     &     + ( dotrc(pz,wp) - dotrc(pwm,wp) ) * epsccrr(wm,z,pz,pwm)
     &     + ( dotrc(pwp,wm) - dotrc(pz,wm) ) * epsccrr(wp,z,pz,pwp)
     &     + ( dotrc(pwm,z) - dotrc(pwp,z) ) * epsccrr(wp,wm,pwm,pwp) )
      endif

      if (fwtilde .ne. 0) then
        vertex = vertex + fwtilde / 2. * gw * wmass**2 / sqrt(cos2w) *
     &      ( epscccr(wp,wm,z,pwm) + epscccr(wp,wm,z,pwp) + cos2w * epscccr(wp,wm,z,pz) )
      endif

c Owtilde as linear combination. Same result but much slower.
c      if (fwtilde .ne. 0) then
c        vertex = vertex - fwtilde / 2. * gw * wmass**2 *
c     &      ( 2. * sqrt(cos2w) * ( epscccr(wp,wm,z,pwm) + epscccr(wp,wm,z,pwp) )
c     &      + ( sqrt(cos2w) + 1. / sqrt(cos2w) ) * epscccr(wp,wm,z,pz) )
c      endif

      if (fbwtilde .ne. 0) then
        vertex = vertex + fbwtilde * gw * wmass**2 * sin2w / sqrt(cos2w) * epscccr(wp,wm,z,pz)
      endif

      return
      end


      subroutine jwpzwm_anomal3(z, wm, jvv)
c
c This subroutine computes an off-shell vector current attached to a W+
c from the three-point gauge boson coupling w+ w- z with anomalous
c couplings.  The vector propagator is given in Feynman gauge for a
c massless vector and in unitary gauge for a massive vector.
c
c input:
c       complex z(6)           : first  vector                        z
c       complex wm(6)          : second vector                        wm
c output:
c       complex jvv(6)         : vector current            j^mu(v:z,wm)
c
c where all the bosons are defined by the flowing-OUT quantum number.

      implicit none
      double complex z(6), wm(6), jvv(6), j12(4), dg, dotcc, dotrc
      double complex epsccrmu, epscrrmu, epsccrr, epscrrr
      double complex cm2        ! mass**2- I Gamma mass (Fabio)
      double precision pz(4),pwm(4),q(4),s
      double precision vm2, dotrr, cos2w, epsrrrmu
      integer i
      external dotcc, dotrc, dotrr, epsccrmu, epsrrrmu, epscrrmu, epsccrr, epscrrr

#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"

      cos2w=1.0d0 - sin2w

      jvv(5) = z(5)+wm(5)
      jvv(6) = z(6)+wm(6)

      pz(1) =  dble( z(5))
      pz(2) =  dble( z(6))
      pz(3) =  dimag(z(6))
      pz(4) =  dimag(z(5))
      pwm(1) =  dble( wm(5))
      pwm(2) =  dble( wm(6))
      pwm(3) =  dimag(wm(6))
      pwm(4) =  dimag(wm(5))
      q(1)  = -dble( jvv(5))
      q(2)  = -dble( jvv(6))
      q(3)  = -dimag(jvv(6))
      q(4)  = -dimag(jvv(5))
      s = dotrr(q,q)
      vm2 = wmass**2

c     Vertex
      do i=1,4
            j12(i) = -gwwz * ( (pz(i)-pwm(i))*dotcc(z,wm)
     &                                     + (dotrc(pwm,z)-dotrc(q,z))*wm(i)
     &                                     + (dotrc(q,wm)-dotrc(pz,wm))*z(i) )

        if (fwww .ne. 0) then
          j12(i) = j12(i) + fwww * 3./2. * gw**3 * sqrt(cos2w) *
     &        ( dotrc(pz,wm) * dotrc(q,z) * pwm(i) - pz(i) * dotrc(q,wm) * dotrc(pwm,z)
     &        + dotrr(pz,pwm) * ( dotrc(q,wm) * z(i) - dotrc(q,z) * wm(i) )
     &        + dotrr(pwm,q) * ( pz(i) * dotcc(wm,z) - dotrc(pz,wm) * z(i) )
     &        + dotrr(q,pz) * ( dotrc(pwm,z) * wm(i) - pwm(i) * dotcc(wm,z) ) )
        endif

        if (fw .ne. 0) then
          j12(i) = j12(i) + fw * gw * wmass**2 / ( 2. * sqrt(cos2w) ) *
     &        ( (pwm(i)-pz(i)) * dotcc(wm,z)
     &        + (dotrc(pz,wm)-dotrc(q,wm)) * z(i) + (dotrc(q,z)-dotrc(pwm,z)) * wm(i)
     &        - sin2w * ( dotrc(pz,wm)*z(i) - pz(i)*dotcc(wm,z) ) )
        endif

        if (fb .ne. 0) then
          j12(i) = j12(i) - fb * gw * wmass**2 / ( 2. * sqrt(cos2w) ) *
     &        sin2w * ( dotrc(pz,wm)*z(i) - pz(i)*dotcc(wm,z) )
        endif

        if (fbw .ne. 0) then
          j12(i) = j12(i) - fbw * gw * wmass**2 * sin2w / sqrt(cos2w) *
     &      ( pz(i)*dotcc(wm,z) - dotrc(pz,wm)*z(i) )
        endif

        if (fdw .ne. 0) then
          j12(i) = j12(i) - fdw * gw**3 * 2. * sqrt(cos2w) * (
     &      - dotcc(wm,z) * ( pz(i) * ( 2.*dotrr(pz,pwm) + dotrr(pz,q) ) - pwm(i) * ( 2.*dotrr(pz,pwm) + dotrr(pwm,q) ) )
     &      + z(i) * ( dotrc(pz,wm) * ( dotrr(pz,pwm) + 2.*dotrr(pz,q) ) - dotrc(q,wm) * ( 2.*dotrr(pz,q) + dotrr(pwm,q) ) )
     &      - wm(i) * ( dotrc(pwm,z) * ( dotrr(pz,pwm) + 2.*dotrr(pwm,q) ) - dotrc(q,z) * ( dotrr(pz,q) + 2.*dotrr(pwm,q) ) )
     &      + pz(i)*dotrc(pz,wm)*dotrc(pwm,z) - pwm(i)*dotrc(pz,wm)*dotrc(pwm,z) - pz(i)*dotrc(pz,wm)*dotrc(q,z)
     &      + pz(i)*dotrc(q,wm)*dotrc(q,z) + pwm(i)*dotrc(q,wm)*dotrc(pwm,z) - pwm(i)*dotrc(q,wm)*dotrc(q,z) )
        endif

        if (fbtilde .ne. 0) then
          j12(i) = j12(i) + fbtilde * gw * wmass**2 / 2. * sin2w / sqrt(cos2w) * epsccrmu(wm,z,pz,i-1)
        endif

c zero due to momentum conservation
c        if (fwwtilde .ne. 0) then
c          j12(i) = j12(i) - fwwtilde * 2. * gw * wmass**2 * sqrt(cos2w) * (
c     &      epsccrmu(wm,z,pwm,i-1) + epsccrmu(wm,z,q,i-1) + epsccrmu(wm,z,pz,i-1) )
c        endif

        if (fwwwtilde .ne. 0) then
          j12(i) = j12(i) + fwwwtilde * gw**3 / 2. * sqrt(cos2w) * (wm(i)*epscrrr(z,pwm,q,pz)
     &      + z(i)*epscrrr(wm,pwm,q,pz) - dotcc(wm,z)*epsrrrmu(pwm,q,pz,i-1)
     &      + pwm(i)*epsccrr(wm,z,q,pz) - pz(i)*epsccrr(wm,z,pwm,q)
     &      + dotrc(q,wm)*epscrrmu(z,pwm,pz,i-1) + dotrc(pz,wm)*epscrrmu(z,pwm,q,i-1)
     &      + dotrc(pwm,z)*epscrrmu(wm,q,pz,i-1) + dotrc(q,z)*epscrrmu(wm,pwm,pz,i-1)
     &      + dotrr(pwm,q)*epsccrmu(wm,z,pz,i-1) + dotrr(pwm,pz)*epsccrmu(wm,z,q,i-1)
     &      + dotrr(q,pz)*epsccrmu(wm,z,pwm,i-1) )
        endif

        if (fdwtilde .ne. 0) then
          j12(i) = j12(i) + fdwtilde * 2. * gw**3 * sqrt(cos2w) * (
     &     + ( dotrr(pz,pwm) + dotrr(pz,q) ) * epsccrmu(wm,z,pz,i-1)
     &     + ( dotrr(pz,pwm) + dotrr(pwm,q) ) * epsccrmu(wm,z,pwm,i-1)
     &     + ( dotrr(pz,q) + dotrr(pwm,q) ) * epsccrmu(wm,z,q,i-1)
     &     + ( pz(i) - pwm(i) ) * epsccrr(wm,z,pz,pwm)
     &     - ( dotrc(q,wm) - dotrc(pz,wm) ) * epscrrmu(z,pz,q,i-1)
     &     - ( dotrc(pwm,z) - dotrc(q,z) ) * epscrrmu(wm,pwm,q,i-1) )
        endif

        if (fwtilde .ne. 0) then
          j12(i) = j12(i) - fwtilde / 2. * gw * wmass**2 / sqrt(cos2w) *
     &      ( epsccrmu(wm,z,pwm,i-1) + epsccrmu(wm,z,q,i-1) + cos2w * epsccrmu(wm,z,pz,i-1) )
        endif

c        if (fwtilde .ne. 0) then
c          j12(i) = j12(i) + fwtilde / 2. * gw * wmass**2 *
c     &      ( 2. * sqrt(cos2w) * ( epsccrmu(wm,z,pwm,i-1) + epsccrmu(wm,z,q,i-1) )
c     &      + ( sqrt(cos2w) + 1. / sqrt(cos2w) ) * epsccrmu(wm,z,pz,i-1) )
c        endif

        if (fbwtilde .ne. 0) then
          j12(i) = j12(i) - fbwtilde * gw * wmass**2 * sin2w / sqrt(cos2w) * epsccrmu(wm,z,pz,i-1)
        endif

      enddo


c     Propagator
      cm2=dcmplx( vm2, -wmass*wwidth )       !Fabio's implementation of the fixed width
      dg = 1.d0 / dcmplx( s-vm2, wmass*wwidth )
      do i=1,4
         jvv(i) = dg * ( j12(i) - q(i)*dotrc(q,j12)/cm2 )
      enddo

      return
      end


      subroutine jwmwpz_anomal3(wp, z, jvv)
c
c This subroutine computes an off-shell vector current attached to a W-
c from the three-point gauge boson coupling w+ w- z with anomalous
c couplings.  The vector propagator is given in Feynman gauge for a
c massless vector and in unitary gauge for a massive vector.
c
c input:
c       complex wp(6)          : first  vector                        wp
c       complex z(6)           : second vector                        z
c output:
c       complex jvv(6)         : vector current wm         j^mu(v:wp,z)
c
c where all the bosons are defined by the flowing-OUT quantum number.

      implicit none
      double complex wp(6), z(6), jvv(6), j12(4), dg, dotcc, dotrc
      double complex epsccrmu, epscrrmu, epsccrr, epscrrr
      double complex cm2        ! mass**2- I Gamma mass (Fabio)
      double precision pwp(4),pz(4),q(4),s
      double precision vm2, dotrr, cos2w, epsrrrmu
      integer i
      external dotcc, dotrc, dotrr, epsccrmu, epsrrrmu, epscrrmu, epsccrr, epscrrr

#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"

      cos2w=1.0d0 - sin2w

      jvv(5) = wp(5)+z(5)
      jvv(6) = wp(6)+z(6)

      pwp(1) =  dble( wp(5))
      pwp(2) =  dble( wp(6))
      pwp(3) =  dimag(wp(6))
      pwp(4) =  dimag(wp(5))
      pz(1) =  dble( z(5))
      pz(2) =  dble( z(6))
      pz(3) =  dimag(z(6))
      pz(4) =  dimag(z(5))
      q(1)  = -dble( jvv(5))
      q(2)  = -dble( jvv(6))
      q(3)  = -dimag(jvv(6))
      q(4)  = -dimag(jvv(5))
      s = dotrr(q,q)
      vm2 = wmass**2

c     Vertex
      do i=1,4
            j12(i) = -gwwz * ( (pwp(i)-pz(i))*dotcc(wp,z)
     &                                     + (dotrc(pz,wp)-dotrc(q,wp))*z(i)
     &                                     + (dotrc(q,z)-dotrc(pwp,z))*wp(i) )

        if (fwww .ne. 0) then
          j12(i) = j12(i) + fwww * 3./2. * gw**3 * sqrt(cos2w) *
     &        ( pz(i) * dotrc(pwp,z) * dotrc(q,wp) - dotrc(pz,wp) * pwp(i) * dotrc(q,z)
     &        + dotrr(pz,q) * ( pwp(i) * dotcc(wp,z) - dotrc(pwp,z) * wp(i) )
     &        + dotrr(q,pwp) * ( dotrc(pz,wp) * z(i) - pz(i) * dotcc(wp,z) )
     &        + dotrr(pwp,pz) * ( dotrc(q,z) * wp(i) - dotrc(q,wp) * z(i) ) )
        endif

        if (fw .ne. 0) then
          j12(i) = j12(i) + fw * gw * wmass**2 / ( 2. * sqrt(cos2w) ) *
     &        ( (dotrc(q,wp)-dotrc(pz,wp)) * z(i)
     &        + (pz(i)-pwp(i)) * dotcc(wp,z) + (dotrc(pwp,z)-dotrc(q,z)) * wp(i)
     &        - sin2w * ( pz(i)*dotcc(wp,z) - dotrc(pz,wp)*z(i) ) )
        endif

        if (fb .ne. 0) then
          j12(i) = j12(i) - fb * gw * wmass**2 / ( 2. * sqrt(cos2w) ) *
     &        sin2w * ( pz(i)*dotcc(wp,z) - dotrc(pz,wp)*z(i) )
        endif

        if (fbw .ne. 0) then
          j12(i) = j12(i) - fbw * gw * wmass**2 * sin2w / sqrt(cos2w) *
     &      ( dotrc(pz,wp)*z(i) - pz(i)*dotcc(wp,z) )
        endif

        if (fdw .ne. 0) then
          j12(i) = j12(i) - fdw * gw**3 * 2. * sqrt(cos2w) * (
     &      - z(i) * ( dotrc(pz,wp) * ( 2.*dotrr(pz,q) + dotrr(pz,pwp) ) - dotrc(q,wp) * ( 2.*dotrr(pz,q) + dotrr(q,pwp) ) )
     &      + dotcc(wp,z) * ( pz(i) * ( dotrr(pz,q) + 2.*dotrr(pz,pwp) ) - pwp(i) * ( 2.*dotrr(pz,pwp) + dotrr(q,pwp) ) )
     &      - wp(i) * ( dotrc(q,z) * ( dotrr(pz,q) + 2.*dotrr(q,pwp) ) - dotrc(pwp,z) * ( dotrr(pz,pwp) + 2.*dotrr(q,pwp) ) )
     &      + dotrc(pz,wp)*pz(i)*dotrc(q,z) - dotrc(q,wp)*pz(i)*dotrc(q,z) - dotrc(pz,wp)*pz(i)*dotrc(pwp,z)
     &      + dotrc(pz,wp)*pwp(i)*dotrc(pwp,z) + dotrc(q,wp)*pwp(i)*dotrc(q,z) - dotrc(q,wp)*pwp(i)*dotrc(pwp,z) )
        endif

        if (fbtilde .ne. 0) then
          j12(i) = j12(i) - fbtilde * gw * wmass**2 / 2. * sin2w / sqrt(cos2w) * epsccrmu(wp,z,pz,i-1)
        endif

c zero due to momentum conservation
c        if (fwwtilde .ne. 0) then
c          j12(i) = j12(i) + fwwtilde * 2. * gw * wmass**2 * sqrt(cos2w) * (
c     &      epsccrmu(wp,z,q,i-1) + epsccrmu(wp,z,pwp,i-1) + epsccrmu(wp,z,pz,i-1) )
c        endif

        if (fwwwtilde .ne. 0) then
          j12(i) = j12(i) + fwwwtilde * gw**3 / 2. * sqrt(cos2w) * ( wp(i)*epscrrr(z,q,pwp,pz)
     &      - dotcc(wp,z)*epsrrrmu(q,pwp,pz,i-1) + z(i)*epscrrr(wp,q,pwp,pz)
     &      - dotrc(q,wp)*epscrrmu(z,pwp,pz,i-1) + dotrc(pz,wp)*epscrrmu(z,q,pwp,i-1)
     &      - pwp(i)*epsccrr(wp,z,q,pz) - pz(i)*epsccrr(wp,z,q,pwp)
     &      - dotrc(q,z)*epscrrmu(wp,pwp,pz,i-1) - dotrc(pwp,z)*epscrrmu(wp,q,pz,i-1)
     &      - dotrr(q,pwp)*epsccrmu(wp,z,pz,i-1) - dotrr(q,pz)*epsccrmu(wp,z,pwp,i-1)
     &      - dotrr(pwp,pz)*epsccrmu(wp,z,q,i-1) )
        endif

        if (fdwtilde .ne. 0) then
          j12(i) = j12(i) + fdwtilde * 2. * gw**3 * sqrt(cos2w) * (
     &     - ( dotrr(pz,q) + dotrr(pz,pwp) ) * epsccrmu(wp,z,pz,i-1)
     &     - ( dotrr(pz,q) + dotrr(q,pwp) ) * epsccrmu(wp,z,q,i-1)
     &     - ( dotrr(pz,pwp) + dotrr(q,pwp) ) * epsccrmu(wp,z,pwp,i-1)
     &     - ( dotrc(pz,wp) - dotrc(q,wp) ) * epscrrmu(z,pz,q,i-1)
     &     + ( pwp(i) - pz(i) ) * epsccrr(wp,z,pz,pwp)
     &     + ( dotrc(q,z) - dotrc(pwp,z) ) * epscrrmu(wp,q,pwp,i-1) )
        endif

        if (fwtilde .ne. 0) then
          j12(i) = j12(i) + fwtilde / 2. * gw * wmass**2 / sqrt(cos2w) *
     &      ( epsccrmu(wp,z,q,i-1) + epsccrmu(wp,z,pwp,i-1) + cos2w * epsccrmu(wp,z,pz,i-1) )
        endif

c        if (fwtilde .ne. 0) then
c          j12(i) = j12(i) - fwtilde / 2. * gw * wmass**2 *
c     &      ( 2. * sqrt(cos2w) * ( epsccrmu(wp,z,q,i-1) + epsccrmu(wp,z,pwp,i-1) )
c     &      + ( sqrt(cos2w) + 1. / sqrt(cos2w) ) * epsccrmu(wp,z,pz,i-1) )
c        endif

        if (fbwtilde .ne. 0) then
          j12(i) = j12(i) + fbwtilde * gw * wmass**2 * sin2w / sqrt(cos2w) * epsccrmu(wp,z,pz,i-1)
        endif

      enddo


c     Propagator
      cm2=dcmplx( vm2, -wmass*wwidth )       !Fabio's implementation of the fixed width
      dg = 1.d0 / dcmplx( s-vm2, wmass*wwidth )
      do i=1,4
         jvv(i) = dg * ( j12(i) - q(i)*dotrc(q,j12)/cm2 )
      enddo

      return
      end


      subroutine jzww_anomal3(wm, wp, jvv)
c
c This subroutine computes an off-shell vector current attached to a Z
c from the three-point gauge boson coupling w+ w- z with anomalous
c couplings.  The vector propagator is given in Feynman gauge for a
c massless vector and in unitary gauge for a massive vector.
c
c input:
c       complex wm(6)          : vector               flow-out W-
c       complex wp(6)          : vector               flow-out W+
c where all the bosons are defined by the flowing-OUT quantum number.
c
c output:
c       complex jvv(6)         : vector current Z          j^mu(v:wm,wp)
c
c where all the bosons are defined by the flowing-OUT quantum number.

      implicit none
      double complex wm(6), wp(6), jvv(6), j12(4), dg, dotcc, dotrc
      double complex epsccrmu, epscrrmu, epsccrr, epscrrr
      double complex cm2        ! mass**2- I Gamma mass (Fabio)
      double precision pwm(4),pwp(4),q(4),s
      double precision vm2, dotrr, cos2w, epsrrrmu
      integer i
      external dotcc, dotrc, dotrr, epsccrmu, epsrrrmu, epscrrmu, epsccrr, epscrrr

#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"

      cos2w=1.0d0 - sin2w

      jvv(5) = wm(5)+wp(5)
      jvv(6) = wm(6)+wp(6)

      pwm(1) =  dble( wm(5))
      pwm(2) =  dble( wm(6))
      pwm(3) =  dimag(wm(6))
      pwm(4) =  dimag(wm(5))
      pwp(1) =  dble( wp(5))
      pwp(2) =  dble( wp(6))
      pwp(3) =  dimag(wp(6))
      pwp(4) =  dimag(wp(5))
      q(1)  = -dble( jvv(5))
      q(2)  = -dble( jvv(6))
      q(3)  = -dimag(jvv(6))
      q(4)  = -dimag(jvv(5))
      s = dotrr(q,q)
      vm2 = zmass**2

c     Vertex
      do i=1,4
            j12(i) = -gwwz * ( dotcc(wm,wp)* pwm(i)
     &                       - dotcc(wm,wp) * pwp(i)
     &                       + dotrc(pwp,wm)* wp(i)
     &                       - dotrc(q,wm)  * wp(i)
     &                       + dotrc(q,wp)  * wm(i)
     &                       - dotrc(pwm,wp)* wm(i) )

        if (fwww .ne. 0) then
          j12(i) = j12(i) + fwww * 3./2. * gw**3 * sqrt(cos2w) *
     &        ( dotrc(q,wm)    * dotrc(pwm,wp) * pwp(i)
     &        - dotrr(q,pwm)   * dotcc(wm,wp)  * pwp(i)
     &        + dotrr(pwp,q)   * dotcc(wm,wp)  * pwm(i)
     &        - dotrc(q,wp)    * dotrc(pwp,wm) * pwm(i)
     &        + dotrr(q,pwm)   * dotrc(pwp,wm) * wp(i)
     &        - dotrr(pwm,pwp) * dotrc(q,wm)   * wp(i)
     &        + dotrr(pwm,pwp) * dotrc(q,wp)   * wm(i)
     &        - dotrr(pwp,q)   * dotrc(pwm,wp) * wm(i)  )
        endif

        if (fw .ne. 0) then
          j12(i) = j12(i) + fw * gw * wmass**2 / ( 2. * sqrt(cos2w) ) *
     &        (  dotrc(pwm,wp)     * wm(i)       !g1
     &         - dotrc(q,wp)       * wm(i)       !deltakappa
     &         + dotrc(q,wm)       * wp(i)       !deltakappa
     &         - dotrc(pwp,wm)     * wp(i)       !g1
     &         + dotcc(wm,wp)      *pwp(i)       !g1
     &         - dotcc(wm,wp)      *pwm(i)       !g1
     &         + sin2w*dotrc(q,wp) * wm(i)       !deltakappa
     &         - sin2w*dotrc(q,wm) * wp(i)       !deltakappa
     &          )
        endif

        if (fb .ne. 0) then
          j12(i) = j12(i) + fb * gw * wmass**2 / ( 2. * sqrt(cos2w) ) *
     &        sin2w *      (   + dotrc(q,wp)*wm(i)   !deltakappa
     &                         - dotrc(q,wm)*wp(i))  !deltakappa
        endif

        if (fbw .ne. 0) then
          j12(i) = j12(i) - fbw * gw * wmass**2 * sin2w / sqrt(cos2w) *
     &      ( dotrc(q,wp)*wm(i) - dotrc(q,wm)*wp(i) )
        endif

        if (fdw .ne. 0) then
          j12(i) = j12(i) - fdw * gw**3 * 2. * sqrt(cos2w) * (
     &      - wm(i) * ( dotrc(q,wp) * ( 2.*dotrr(q,pwm) + dotrr(q,pwp) ) - dotrc(pwm,wp) * ( 2.*dotrr(q,pwm) + dotrr(pwm,pwp) ) )
     &      + wp(i) * ( dotrc(q,wm) * ( dotrr(q,pwm) + 2.*dotrr(q,pwp) ) - dotrc(pwp,wm) * ( 2.*dotrr(q,pwp) + dotrr(pwm,pwp) ) )
     &      - dotcc(wp,wm) * ( pwm(i) * ( dotrr(q,pwm) + 2.*dotrr(pwm,pwp) ) - pwp(i) * ( dotrr(q,pwp) + 2.*dotrr(pwm,pwp) ) )
     &      + dotrc(q,wp)*dotrc(q,wm)*pwm(i) - dotrc(pwm,wp)*dotrc(q,wm)*pwm(i) - dotrc(q,wp)*dotrc(q,wm)*pwp(i)
     &      + dotrc(q,wp)*dotrc(pwp,wm)*pwp(i) + dotrc(pwm,wp)*dotrc(pwp,wm)*pwm(i) - dotrc(pwm,wp)*dotrc(pwp,wm)*pwp(i) )
        endif

        if (fbtilde .ne. 0) then
          j12(i) = j12(i) + fbtilde * gw * wmass**2 / 2. * sin2w / sqrt(cos2w) * epsccrmu(wp,wm,q,i-1)
        endif

c zero due to momentum conservation
c        if (fwwtilde .ne. 0) then
c          j12(i) = j12(i) - fwwtilde * 2. * gw * wmass**2 * sqrt(cos2w) * (
c     &      epsccrmu(wp,wm,pwm,i-1) + epsccrmu(wp,wm,pwp,i-1) + epsccrmu(wp,wm,q,i-1) )
c        endif

        if (fwwwtilde .ne. 0) then
          j12(i) = j12(i) + fwwwtilde * gw**3 / 2. * sqrt(cos2w) * ( - dotcc(wp,wm)*epsrrrmu(pwm,pwp,q,i-1)
     &      + wp(i)*epscrrr(wm,pwm,pwp,q) + wm(i)*epscrrr(wp,pwm,pwp,q)
     &      + dotrc(pwm,wp)*epscrrmu(wm,pwp,q,i-1) - dotrc(q,wp)*epscrrmu(wm,pwm,pwp,i-1)
     &      - dotrc(pwp,wm)*epscrrmu(wp,pwm,q,i-1) - dotrc(q,wm)*epscrrmu(wp,pwm,pwp,i-1)
     &      - pwm(i)*epsccrr(wp,wm,pwp,q) - pwp(i)*epsccrr(wp,wm,pwm,q)
     &      + dotrr(pwm,pwp)*epsccrmu(wp,wm,q,i-1) + dotrr(pwm,q)*epsccrmu(wp,wm,pwp,i-1)
     &      + dotrr(pwp,q)*epsccrmu(wp,wm,pwm,i-1) )
        endif

        if (fdwtilde .ne. 0) then
          j12(i) = j12(i) + fdwtilde * 2. * gw**3 * sqrt(cos2w) * (
     &     + ( dotrr(q,pwm) + dotrr(q,pwp) ) * epsccrmu(wp,wm,q,i-1)
     &     + ( dotrr(q,pwm) + dotrr(pwm,pwp) ) * epsccrmu(wp,wm,pwm,i-1)
     &     + ( dotrr(q,pwp) + dotrr(pwm,pwp) ) * epsccrmu(wp,wm,pwp,i-1)
     &     + ( dotrc(q,wp) - dotrc(pwm,wp) ) * epscrrmu(wm,q,pwm,i-1)
     &     + ( dotrc(pwp,wm) - dotrc(q,wm) ) * epscrrmu(wp,q,pwp,i-1)
     &     + ( pwm(i) - pwp(i) ) * epsccrr(wp,wm,pwm,pwp) )
        endif

        if (fwtilde .ne. 0) then
          j12(i) = j12(i) - fwtilde / 2. * gw * wmass**2 / sqrt(cos2w) *
     &      ( epsccrmu(wp,wm,pwm,i-1) + epsccrmu(wp,wm,pwp,i-1) + cos2w * epsccrmu(wp,wm,q,i-1) )
        endif

c        if (fwtilde .ne. 0) then
c          j12(i) = j12(i) + fwtilde / 2. * gw * wmass**2 *
c     &      ( 2. * sqrt(cos2w) * ( epsccrmu(wp,wm,pwm,i-1) + epsccrmu(wp,wm,pwp,i-1) )
c     &      + ( sqrt(cos2w) + 1. / sqrt(cos2w) ) * epsccrmu(wp,wm,q,i-1) )
c        endif

        if (fbwtilde .ne. 0) then
          j12(i) = j12(i) - fbwtilde * gw * wmass**2 * sin2w / sqrt(cos2w) * epsccrmu(wp,wm,q,i-1)
        endif

      enddo


c     Propagator
      cm2=dcmplx( vm2, -zmass*zwidth )       !Fabio's implementation of the fixed width
      dg = 1.d0 / dcmplx( s-vm2, zmass*zwidth )
      do i=1,4
         jvv(i) = dg * ( j12(i) - q(i)*dotrc(q,j12)/cm2 )
      enddo


      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCC
C Vertices with Higgs C
CCCCCCCCCCCCCCCCCCCCCCC


      subroutine zzh_anomal3(v1, v2, sc, vertex)
c
c This subroutine computes an amplitude of the vector-vector-scalar
c coupling z z h with anomalous couplings.
c
c input:
c       complex z1(6)          : first  vector                        z1
c       complex z2(6)          : second vector                        z2
c       complex sc(3)          : input  scalar                        s
c
c output:
c       complex vertex         : amplitude                gamma(z1,z2,s)
c
c where all the bosons are defined by the flowing-OUT quantum number.

      implicit none
      double complex v1(6),v2(6),sc(3),vertex, dotcc, dotrc, epsccrr
      double precision p1(4),p2(4),q(4), dotrr, cos2w
      external dotcc, dotrc, dotrr, epsccrr

#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"

      cos2w=1.0d0 - sin2w

      p1(1) = dble( v1(5))
      p1(2) = dble( v1(6))
      p1(3) = dimag(v1(6))
      p1(4) = dimag(v1(5))
      p2(1) = dble( v2(5))
      p2(2) = dble( v2(6))
      p2(3) = dimag(v2(6))
      p2(4) = dimag(v2(5))
      q(1)  = dble( sc(2))
      q(2)  = dble( sc(3))
      q(3)  = dimag(sc(3))
      q(4)  = dimag(sc(2))

      vertex = gw*zmass/sqrt(cos2w)*sc(1)*dotcc(v1,v2)
      if (with_secondhiggs .and. use_secondhiggs) then
        vertex = vertex*gcba
      elseif (with_secondhiggs) then
        vertex = vertex*gsba
      endif

      if (fphi1 .ne. 0) then
         vertex = vertex + fphi1 * 4. * wmass**3 / gw / cos2w * sc(1) * dotcc(v1,v2)
      endif

      if (fw .ne. 0) then
         vertex = vertex + fw * gw * wmass / 2. * sc(1) * ( (dotrr(q,p1) + dotrr(q,p2)) * dotcc(v1,v2)
     &      - dotrc(q,v1)*dotrc(p1,v2) - dotrc(q,v2)*dotrc(p2,v1) )
      endif

      if (fb .ne. 0) then
         vertex = vertex + fb * gw * sin2w / cos2w * wmass / 2. * sc(1) * ( (dotrr(q,p1) + dotrr(q,p2)) * dotcc(v1,v2)
     &      - dotrc(q,v1)*dotrc(p1,v2) - dotrc(q,v2)*dotrc(p2,v1) )
      endif

      if (fww .ne. 0) then
         vertex = vertex - fww * 2. * gw * wmass * cos2w  * sc(1) * ( - dotrr(p1,p2) * dotcc(v1,v2)
     &      + dotrc(p2,v1)*dotrc(p1,v2) )
      endif

      if (fbb .ne. 0) then
         vertex = vertex - fbb * 2. * gw * wmass * sin2w**2. / cos2w  * sc(1) * ( - dotrr(p1,p2) * dotcc(v1,v2)
     &      + dotrc(p2,v1)*dotrc(p1,v2) )
      endif

      if (fbw .ne. 0) then
         vertex = vertex - fbw * 2. * gw * wmass * sin2w * sc(1) * ( - dotrr(p1,p2) * dotcc(v1,v2)
     &      + dotrc(p2,v1)*dotrc(p1,v2) )
      endif

      if (fbtilde .ne. 0) then
         vertex = vertex + fbtilde * gw * wmass / 2. * sin2w / cos2w * sc(1) * ( epsccrr(v1,v2,q,p1)
     &      - epsccrr(v1,v2,q,p2) )
      endif

      if (fwwtilde .ne. 0) then
         vertex = vertex - fwwtilde * gw * wmass * 2. * cos2w * sc(1) * epsccrr(v1,v2,p1,p2)
      endif

      if (fbbtilde .ne. 0) then
         vertex = vertex - fbbtilde * gw * wmass * 2. * sin2w**2 / cos2w * sc(1) * epsccrr(v1,v2,p1,p2)
      endif

      if (fwtilde .ne. 0) then
         vertex = vertex + fwtilde * gw * wmass / 2. * sc(1) * ( epsccrr(v1,v2,q,p1) - epsccrr(v1,v2,q,p2) )
      endif

c Owtilde as linear combination. Same result but much slower.
c      if (fwtilde .ne. 0) then
c         vertex = vertex + fwtilde * gw * wmass / 2. / cos2w * sc(1) *
c     &      ( sin2w * ( epsccrr(v1,v2,q,p1) - epsccrr(v1,v2,q,p2) - 2. * epsccrr(v1,v2,p1,p2) )
c     &      + 2. * cos2w * epsccrr(v1,v2,p1,p2) )
c      endif

      if (fbwtilde .ne. 0) then
         vertex = vertex - fbwtilde * gw * wmass * 2. * sin2w * sc(1) * epsccrr(v1,v2,p1,p2)
      endif

      return
      end


      subroutine jhzz_anomal3(v1,v2 , hvv)
c
c This subroutine computes an off-shell scalar current from the vector-
c vector-scalar coupling z z h with anomalous couplings.
c
c input:
c       complex v1(6)          : first  vector                        z1
c       complex v2(6)          : second vector                        z2
c
c output:
c       complex hvv(3)         : off-shell scalar current     j(s:z1,z2)
c
c where all the bosons are defined by the flowing-OUT quantum number.

      implicit none
      double complex v1(6),v2(6),hvv(3), dotrc, dotcc, epsccrr
      double precision p1(4),p2(4),q(4),q2, dotrr, cos2w
      external dotcc, dotrc, dotrr, epsccrr

#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"

      cos2w=1.0d0 - sin2w

      hvv(2) = v1(5)+v2(5)
      hvv(3) = v1(6)+v2(6)

      p1(1) = dble( v1(5))
      p1(2) = dble( v1(6))
      p1(3) = dimag(v1(6))
      p1(4) = dimag(v1(5))
      p2(1) = dble( v2(5))
      p2(2) = dble( v2(6))
      p2(3) = dimag(v2(6))
      p2(4) = dimag(v2(5))
      q(1) = -dble( hvv(2))
      q(2) = -dble( hvv(3))
      q(3) = -dimag(hvv(3))
      q(4) = -dimag(hvv(2))
      q2 = dotrr(q,q)

c     Vertex:
      hvv(1) = gw*zmass/sqrt(cos2w)*dotcc(v1,v2)
      if (with_secondhiggs .and. use_secondhiggs) then
        hvv(1) = hvv(1)*gcba
      elseif (with_secondhiggs) then
        hvv(1) = hvv(1)*gsba
      endif

      if (fphi1 .ne. 0) then
         hvv(1) = hvv(1) + fphi1 * 4. * wmass**3 / gw / cos2w * dotcc(v1,v2)
      endif

      if (fw .ne. 0) then
         hvv(1) = hvv(1) + fw * gw * wmass / 2. * ( (dotrr(q,p1) + dotrr(q,p2)) * dotcc(v1,v2)
     &      - dotrc(q,v1)*dotrc(p1,v2) - dotrc(q,v2)*dotrc(p2,v1) )
      endif

      if (fb .ne. 0) then
         hvv(1) = hvv(1) + fb * gw * sin2w / cos2w * wmass / 2. * ( (dotrr(q,p1) + dotrr(q,p2)) * dotcc(v1,v2)
     &      - dotrc(q,v1)*dotrc(p1,v2) - dotrc(q,v2)*dotrc(p2,v1) )
      endif

      if (fww .ne. 0) then
         hvv(1) = hvv(1) - fww * 2. * gw * wmass * cos2w  * ( - dotrr(p1,p2) * dotcc(v1,v2)
     &      + dotrc(p2,v1)*dotrc(p1,v2) )
      endif

      if (fbb .ne. 0) then
         hvv(1) = hvv(1) - fbb * 2. * gw * wmass * sin2w**2. / cos2w  * ( - dotrr(p1,p2) * dotcc(v1,v2)
     &      + dotrc(p2,v1)*dotrc(p1,v2) )
      endif

      if (fbw .ne. 0) then
         hvv(1) = hvv(1) - fbw * 2. * gw * wmass * sin2w * ( - dotrr(p1,p2) * dotcc(v1,v2)
     &      + dotrc(p2,v1)*dotrc(p1,v2) )
      endif

      if (fbtilde .ne. 0) then
         hvv(1) = hvv(1) + fbtilde * gw * wmass / 2. * sin2w / cos2w * ( epsccrr(v1,v2,q,p1)
     &      - epsccrr(v1,v2,q,p2) )
      endif

      if (fwwtilde .ne. 0) then
         hvv(1) = hvv(1) - fwwtilde * gw * wmass * 2. * cos2w * epsccrr(v1,v2,p1,p2)
      endif

      if (fbbtilde .ne. 0) then
         hvv(1) = hvv(1) - fbbtilde * gw * wmass * 2. * sin2w**2 / cos2w * epsccrr(v1,v2,p1,p2)
      endif

      if (fwtilde .ne. 0) then
         hvv(1) = hvv(1) + fwtilde * gw * wmass / 2. * ( epsccrr(v1,v2,q,p1) - epsccrr(v1,v2,q,p2) )
      endif

c Owtilde as linear combination. Same result but much slower.
c      if (fwtilde .ne. 0) then
c         hvv(1) = hvv(1) + fwtilde * gw * wmass / 2. / cos2w *
c     &      ( sin2w * ( epsccrr(v1,v2,q,p1) - epsccrr(v1,v2,q,p2) - 2. * epsccrr(v1,v2,p1,p2) )
c     &      + 2. * cos2w * epsccrr(v1,v2,p1,p2) )
c      endif

      if (fbwtilde .ne. 0) then
         hvv(1) = hvv(1) - fbwtilde * gw * wmass * 2. * sin2w * epsccrr(v1,v2,p1,p2)
      endif


c     Propagator:
      if (with_secondhiggs .and. use_secondhiggs) then
        hvv(1) = - hvv(1) / dcmplx( q2-h2mass**2, h2mass*h2width )
      else
        hvv(1) = - hvv(1) / dcmplx( q2-hmass**2, hmass*hwidth )
      endif

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCcc

      subroutine azh_anomal3(a, z, sc, vertex)
c
c This subroutine computes an amplitude of the vector-vector-scalar
c coupling a z h with anomalous couplings.
c
c input:
c       complex a(6)          : photon
c       complex z(6)          : z-boson
c       complex sc(3)          : higgs
c
c output:
c       complex vertex         : amplitude                gamma(a,z,s)
c
c where all the bosons are defined by the flowing-OUT quantum number.

      implicit none
      double complex a(6),z(6),sc(3),vertex, dotcc, dotrc, epsccrr
      double precision pa(4),pz(4),q(4), dotrr, cos2w
      external dotcc, dotrc, dotrr, epsccrr

#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"

      cos2w=1.0d0 - sin2w

      pa(1) = dble( a(5))
      pa(2) = dble( a(6))
      pa(3) = dimag(a(6))
      pa(4) = dimag(a(5))
      pz(1) = dble( z(5))
      pz(2) = dble( z(6))
      pz(3) = dimag(z(6))
      pz(4) = dimag(z(5))
      q(1)  = dble( sc(2))
      q(2)  = dble( sc(3))
      q(3)  = dimag(sc(3))
      q(4)  = dimag(sc(2))

      vertex = dcmplx(0.,0.)

      if (fw .ne. 0) then
         vertex = vertex + fw * gw * wmass * sqrt(sin2w/cos2w) / 2. * sc(1) * ( dotrr(pa,q) * dotcc(a,z)
     &      - dotrc(q,a)*dotrc(pa,z) )
      endif

      if (fb .ne. 0) then
         vertex = vertex + fb * gw * wmass * sqrt(sin2w/cos2w) / 2. * sc(1) * ( - dotrr(pa,q) * dotcc(a,z)
     &      + dotrc(q,a)*dotrc(pa,z) )
      endif

      if (fww .ne. 0) then
         vertex = vertex - fww * 2. * gw * wmass * sqrt(sin2w*cos2w)  * sc(1) * ( - dotrr(pa,pz) * dotcc(a,z)
     &      + dotrc(pz,a)*dotrc(pa,z) )
      endif

      if (fbb .ne. 0) then
         vertex = vertex + fbb * 2. * gw * wmass * sqrt(sin2w/cos2w) * sin2w  * sc(1) * ( - dotrr(pa,pz) * dotcc(a,z)
     &      + dotrc(pz,a)*dotrc(pa,z) )
      endif

      if (fbw .ne. 0) then
         vertex = vertex + fbw * gw * wmass * sqrt(sin2w/cos2w) * sc(1) *
     &      ( - dotcc(a,z)*dotrr(pa,pz) * ( cos2w - sin2w )
     &      + ( cos2w - sin2w ) * dotrc(pz,a)*dotrc(pa,z) )
      endif

      if (fbtilde .ne. 0) then
         vertex = vertex + fbtilde * gw * wmass / 2. * sqrt(sin2w/cos2w) * sc(1) * epsccrr(a,z,pa,q)
      endif

      if (fwwtilde .ne. 0) then
         vertex = vertex - fwwtilde * gw * wmass * 2. * sqrt(sin2w*cos2w) * sc(1) * epsccrr(a,z,pa,pz)
      endif

      if (fbbtilde .ne. 0) then
         vertex = vertex + fbbtilde * gw * wmass * 2. * sin2w * sqrt(sin2w/cos2w) * sc(1) * epsccrr(a,z,pa,pz)
      endif

      if (fwtilde .ne. 0) then
         vertex = vertex - fwtilde * gw * wmass / 2. * sqrt(sin2w/cos2w) * sc(1) * epsccrr(a,z,pa,q)
      endif

c Owtilde as linear combination. Same result but much slower.
c      if (fwtilde .ne. 0) then
c         vertex = vertex + fwtilde * gw * wmass / 2. * sqrt(sin2w/cos2w) * sc(1) *
c     &       ( epsccrr(a,z,pa,q) + 2. * epsccrr(a,z,pa,pz) )
c      endif

      if (fbwtilde .ne. 0) then
         vertex = vertex + fbwtilde * gw * wmass * sqrt(sin2w) * sc(1) *
     &      ( sqrt(cos2w) - sin2w / sqrt(cos2w) ) * epsccrr(a,z,pa,pz)
      endif

      return
      end


      subroutine jhaz_anomal3(a, z, hvv)
c
c This subroutine computes an off-shell scalar current from the vector-
c vector-scalar coupling a z h with anomalous couplings.
c
c input:
c       complex a(6)          : first  vector                        a
c       complex z(6)          : second vector                        z
c
c output:
c       complex hvv(3)         : off-shell scalar current     j(s:a,z)
c
c where all the bosons are defined by the flowing-OUT quantum number.

      implicit none
      double complex a(6),z(6),hvv(3), dotrc, dotcc, epsccrr
      double precision pa(4),pz(4),q(4),q2, dotrr, cos2w
      external dotcc, dotrc, dotrr, epsccrr

#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"

      cos2w=1.0d0 - sin2w

      hvv(2) = a(5)+z(5)
      hvv(3) = a(6)+z(6)

      pa(1) = dble( a(5))
      pa(2) = dble( a(6))
      pa(3) = dimag(a(6))
      pa(4) = dimag(a(5))
      pz(1) = dble( z(5))
      pz(2) = dble( z(6))
      pz(3) = dimag(z(6))
      pz(4) = dimag(z(5))
      q(1) = -dble( hvv(2))
      q(2) = -dble( hvv(3))
      q(3) = -dimag(hvv(3))
      q(4) = -dimag(hvv(2))
      q2 = dotrr(q,q)

c     Vertex:
      hvv(1) = dcmplx(0.,0.)


      if (fw .ne. 0) then
         hvv(1) = hvv(1) + fw * gw * wmass * sqrt(sin2w/cos2w) / 2. * ( dotrr(pa,q) * dotcc(a,z)
     &      - dotrc(q,a)*dotrc(pa,z) )
      endif

      if (fb .ne. 0) then
         hvv(1) = hvv(1) + fb * gw * wmass * sqrt(sin2w/cos2w) / 2. * ( - dotrr(pa,q) * dotcc(a,z)
     &      + dotrc(q,a)*dotrc(pa,z) )
      endif

      if (fww .ne. 0) then
         hvv(1) = hvv(1) - fww * 2. * gw * wmass * sqrt(sin2w*cos2w)  * ( - dotrr(pa,pz) * dotcc(a,z)
     &      + dotrc(pz,a)*dotrc(pa,z) )
      endif

      if (fbb .ne. 0) then
         hvv(1) = hvv(1) + fbb * 2. * gw * wmass * sqrt(sin2w/cos2w) * sin2w  * ( - dotrr(pa,pz) * dotcc(a,z)
     &      + dotrc(pz,a)*dotrc(pa,z) )
      endif

      if (fbw .ne. 0) then
         hvv(1) = hvv(1) + fbw * gw * wmass * sqrt(sin2w/cos2w) *
     &      ( - dotcc(a,z)*dotrr(pa,pz) * ( cos2w - sin2w )
     &      + ( cos2w - sin2w ) * dotrc(pz,a)*dotrc(pa,z) )
      endif

      if (fbtilde .ne. 0) then
         hvv(1) = hvv(1) + fbtilde * gw * wmass / 2. * sqrt(sin2w/cos2w) * epsccrr(a,z,pa,q)
      endif

      if (fwwtilde .ne. 0) then
         hvv(1) = hvv(1) - fwwtilde * gw * wmass * 2. * sqrt(sin2w*cos2w) * epsccrr(a,z,pa,pz)
      endif

      if (fbbtilde .ne. 0) then
         hvv(1) = hvv(1) + fbbtilde * gw * wmass * 2. * sin2w * sqrt(sin2w/cos2w) * epsccrr(a,z,pa,pz)
      endif

      if (fwtilde .ne. 0) then
         hvv(1) = hvv(1) - fwtilde * gw * wmass / 2. * sqrt(sin2w/cos2w) * epsccrr(a,z,pa,q)
      endif

c Owtilde as linear combination. Same result but much slower.
c      if (fwtilde .ne. 0) then
c         hvv(1) = hvv(1) + fwtilde * gw * wmass / 2. * sqrt(sin2w/cos2w) *
c     &       ( epsccrr(a,z,pa,q) + 2. * epsccrr(a,z,pa,pz) )
c      endif

      if (fbwtilde .ne. 0) then
         hvv(1) = hvv(1) + fbwtilde * gw * wmass * sqrt(sin2w) *
     &      ( sqrt(cos2w) - sin2w / sqrt(cos2w) ) * epsccrr(a,z,pa,pz)
      endif


c     Propagator:
      if (with_secondhiggs .and. use_secondhiggs) then
        hvv(1) = - hvv(1) / dcmplx( q2-h2mass**2, h2mass*h2width )
      else
        hvv(1) = - hvv(1) / dcmplx( q2-hmass**2, hmass*hwidth )
      endif

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine aah_anomal3(v1, v2, sc, vertex)
c
c This subroutine computes an amplitude of the vector-vector-scalar
c coupling a a h with anomalous couplings.
c
c input:
c       complex a1(6)          : first  vector                        a1
c       complex a2(6)          : second vector                        a2
c       complex sc(3)          : input  scalar                        s
c
c output:
c       complex vertex         : amplitude                gamma(a1,a2,s)
c
c where all the bosons are defined by the flowing-OUT quantum number.

      implicit none
      double complex v1(6),v2(6),sc(3),vertex, dotcc, dotrc, epsccrr
      double precision p1(4),p2(4),q(4), dotrr, cos2w
      external dotcc, dotrc, dotrr, epsccrr

#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"

      cos2w=1.0d0 - sin2w

      p1(1) = dble( v1(5))
      p1(2) = dble( v1(6))
      p1(3) = dimag(v1(6))
      p1(4) = dimag(v1(5))
      p2(1) = dble( v2(5))
      p2(2) = dble( v2(6))
      p2(3) = dimag(v2(6))
      p2(4) = dimag(v2(5))
      q(1)  = dble( sc(2))
      q(2)  = dble( sc(3))
      q(3)  = dimag(sc(3))
      q(4)  = dimag(sc(2))

      vertex = dcmplx(0.,0.)

      if (fww .ne. 0) then
         vertex = vertex - fww * 2. * gw * wmass * sin2w  * sc(1) * ( - dotrr(p1,p2) * dotcc(v1,v2)
     &      + dotrc(p2,v1)*dotrc(p1,v2) )
      endif

      if (fbb .ne. 0) then
         vertex = vertex - fbb * 2. * gw * wmass * sin2w  * sc(1) * ( - dotrr(p1,p2) * dotcc(v1,v2)
     &      + dotrc(p2,v1)*dotrc(p1,v2) )
      endif

      if (fbw .ne. 0) then
         vertex = vertex + fbw * 2. * gw * wmass * sin2w  * sc(1) * ( - dotrr(p1,p2) * dotcc(v1,v2)
     &      + dotrc(p2,v1)*dotrc(p1,v2) )
      endif

      if (fwwtilde .ne. 0) then
         vertex = vertex - fwwtilde * gw * wmass * 2. * sin2w * sc(1) * epsccrr(v1,v2,p1,p2)
      endif

      if (fbbtilde .ne. 0) then
         vertex = vertex - fbbtilde * gw * wmass * 2. * sin2w * sc(1) * epsccrr(v1,v2,p1,p2)
      endif

      if (fbwtilde .ne. 0) then
         vertex = vertex + fbwtilde * gw * wmass * 2. * sin2w * sc(1) * epsccrr(v1,v2,p1,p2)
      endif

      return
      end


      subroutine jhaa_anomal3(v1,v2 , hvv)
c
c This subroutine computes an off-shell scalar current from the vector-
c vector-scalar coupling a a h with anomalous couplings.
c
c input:
c       complex v1(6)          : first  vector                        a1
c       complex v2(6)          : second vector                        a2
c
c output:
c       complex hvv(3)         : off-shell scalar current     j(s:a1,a2)
c
c where all the bosons are defined by the flowing-OUT quantum number.

      implicit none
      double complex v1(6),v2(6),hvv(3), dotrc, dotcc, epsccrr
      double precision p1(4),p2(4),q(4),q2, dotrr, cos2w
      external dotcc, dotrc, dotrr, epsccrr

#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"

      cos2w=1.0d0 - sin2w

      hvv(2) = v1(5)+v2(5)
      hvv(3) = v1(6)+v2(6)

      p1(1) = dble( v1(5))
      p1(2) = dble( v1(6))
      p1(3) = dimag(v1(6))
      p1(4) = dimag(v1(5))
      p2(1) = dble( v2(5))
      p2(2) = dble( v2(6))
      p2(3) = dimag(v2(6))
      p2(4) = dimag(v2(5))
      q(1) = -dble( hvv(2))
      q(2) = -dble( hvv(3))
      q(3) = -dimag(hvv(3))
      q(4) = -dimag(hvv(2))
      q2 = dotrr(q,q)

c     Vertex:
      hvv(1) = dcmplx(0.,0.)

      if (fww .ne. 0) then
         hvv(1) = hvv(1) - fww * 2. * gw * wmass * sin2w  * ( - dotrr(p1,p2) * dotcc(v1,v2)
     &      + dotrc(p2,v1)*dotrc(p1,v2) )
      endif

      if (fbb .ne. 0) then
         hvv(1) = hvv(1) - fbb * 2. * gw * wmass * sin2w  * ( - dotrr(p1,p2) * dotcc(v1,v2)
     &      + dotrc(p2,v1)*dotrc(p1,v2) )
      endif

      if (fbw .ne. 0) then
         hvv(1) = hvv(1) + fbw * 2. * gw * wmass * sin2w  * ( - dotrr(p1,p2) * dotcc(v1,v2)
     &      + dotrc(p2,v1)*dotrc(p1,v2) )
      endif

      if (fwwtilde .ne. 0) then
         hvv(1) = hvv(1) - fwwtilde * gw * wmass * 2. * sin2w * epsccrr(v1,v2,p1,p2)
      endif

      if (fbbtilde .ne. 0) then
         hvv(1) = hvv(1) - fbbtilde * gw * wmass * 2. * sin2w * epsccrr(v1,v2,p1,p2)
      endif

      if (fbwtilde .ne. 0) then
         hvv(1) = hvv(1) + fbwtilde * gw * wmass * 2. * sin2w * epsccrr(v1,v2,p1,p2)
      endif


c     Propagator:
      if (with_secondhiggs .and. use_secondhiggs) then
        hvv(1) = - hvv(1) / dcmplx( q2-h2mass**2, h2mass*h2width )
      else
        hvv(1) = - hvv(1) / dcmplx( q2-hmass**2, hmass*hwidth )
      endif

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine wwh_anomal3(v1, v2, sc, vertex)
c
c This subroutine computes an amplitude of the vector-vector-scalar
c coupling w+ w- h with anomalous couplings.
c
c input:
c       complex wm(6)          : first  vector                        wm
c       complex wp(6)          : second vector                        wp
c       complex sc(3)          : input  scalar                        s
c
c output:
c       complex vertex         : amplitude                gamma(wm,wp,s)
c
c where all the bosons are defined by the flowing-OUT quantum number.

      implicit none
      double complex v1(6),v2(6),sc(3),vertex, dotcc, dotrc, epsccrr
      double precision p1(4),p2(4),q(4), dotrr, cos2w
      external dotcc, dotrc, dotrr, epsccrr

#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"

      cos2w=1.0d0 - sin2w

      p1(1) = dble( v1(5))
      p1(2) = dble( v1(6))
      p1(3) = dimag(v1(6))
      p1(4) = dimag(v1(5))
      p2(1) = dble( v2(5))
      p2(2) = dble( v2(6))
      p2(3) = dimag(v2(6))
      p2(4) = dimag(v2(5))
      q(1)  = dble( sc(2))
      q(2)  = dble( sc(3))
      q(3)  = dimag(sc(3))
      q(4)  = dimag(sc(2))

      vertex = gw*wmass*sc(1)*dotcc(v1,v2)
      if (with_secondhiggs .and. use_secondhiggs) then
        vertex = vertex*gcba
      elseif (with_secondhiggs) then
        vertex = vertex*gsba
      endif

      if (fw .ne. 0) then
         vertex = vertex + fw * gw * wmass / 2. * sc(1) * ( (dotrr(q,p1) + dotrr(q,p2)) * dotcc(v1,v2)
     &      - dotrc(q,v1)*dotrc(p1,v2) - dotrc(q,v2)*dotrc(p2,v1) )
      endif

      if (fww .ne. 0) then
         vertex = vertex - fww * 2.* gw * wmass * sc(1) * ( - dotrr(p1,p2) * dotcc(v1,v2)
     &      + dotrc(p1,v2)*dotrc(p2,v1) )
      endif

      if (fwwtilde .ne. 0) then
         vertex = vertex + fwwtilde * gw * wmass * 2. * sc(1) * epsccrr(v2,v1,p1,p2)
      endif

      if (fwtilde .ne. 0) then
         vertex = vertex - fwtilde * gw * wmass / 2. * sc(1) * ( epsccrr(v2,v1,q,p1) - epsccrr(v2,v1,q,p2) )
      endif

c Owtilde as linear combination. Same result but much slower.
c      if (fwtilde .ne. 0) then
c         vertex = vertex - fwtilde * gw * wmass * sc(1) * epsccrr(v2,v1,p1,p2)
c      endif

      return
      end


      subroutine jhww_anomal3(v1,v2 , hvv)
c
c This subroutine computes an off-shell scalar current from the vector-
c vector-scalar coupling w w h with anomalous couplings.
c
c input:
c       complex wm(6)          : first  vector                        wm
c       complex wp(6)          : second vector                        wp
c
c output:
c       complex hvv(3)         : off-shell scalar current     j(s:wm,wp)
c
c where all the bosons are defined by the flowing-OUT quantum number.

      implicit none
      double complex v1(6),v2(6),hvv(3), dotrc, dotcc, epsccrr
      double precision p1(4),p2(4),q(4),q2, dotrr
      external dotcc, dotrc, dotrr, epsccrr

#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"

      hvv(2) = v1(5)+v2(5)
      hvv(3) = v1(6)+v2(6)

      p1(1) = dble( v1(5))
      p1(2) = dble( v1(6))
      p1(3) = dimag(v1(6))
      p1(4) = dimag(v1(5))
      p2(1) = dble( v2(5))
      p2(2) = dble( v2(6))
      p2(3) = dimag(v2(6))
      p2(4) = dimag(v2(5))
      q(1) = -dble( hvv(2))
      q(2) = -dble( hvv(3))
      q(3) = -dimag(hvv(3))
      q(4) = -dimag(hvv(2))
      q2 = dotrr(q,q)

c     Vertex:
      hvv(1) = gw*wmass*dotcc(v1,v2)
      if (with_secondhiggs .and. use_secondhiggs) then
        hvv(1) = hvv(1)*gcba
      elseif (with_secondhiggs) then
        hvv(1) = hvv(1)*gsba
      endif

      if (fw .ne. 0) then
         hvv(1) = hvv(1) + fw * gw * wmass / 2. * ( (dotrr(q,p1) + dotrr(q,p2)) * dotcc(v1,v2)
     &      - dotrc(q,v1)*dotrc(p1,v2) - dotrc(q,v2)*dotrc(p2,v1) )
      endif

      if (fww .ne. 0) then
         hvv(1) = hvv(1) - fww * 2. * gw * wmass * ( - dotrr(p1,p2) * dotcc(v1,v2)
     &      + dotrc(p1,v2)*dotrc(p2,v1) )
      endif

      if (fwwtilde .ne. 0) then
         hvv(1) = hvv(1) + fwwtilde * gw * wmass * 2. * epsccrr(v2,v1,p1,p2)
      endif

      if (fwtilde .ne. 0) then
         hvv(1) = hvv(1) - fwtilde * gw * wmass / 2. * ( epsccrr(v2,v1,q,p1) - epsccrr(v2,v1,q,p2) )
      endif

c      if (fwtilde .ne. 0) then
c         hvv(1) = hvv(1) - fwtilde * gw * wmass * epsccrr(v2,v1,p1,p2)
c      endif

c     Propagator:
      if (with_secondhiggs .and. use_secondhiggs) then
        hvv(1) = - hvv(1) / dcmplx( q2-h2mass**2, h2mass*h2width )
      else
        hvv(1) = - hvv(1) / dcmplx( q2-hmass**2, hmass*hwidth )
      endif

      return
      end
