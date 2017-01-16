c
c This file contains routines that substitute the helas calls for the quartic-vector-boson
c couplings with anomalous couplings.
c
c   LIST OF ALL SUBROUTINES IN THIS FILE:
c
c subroutine wwaa_anomal4(wm,a1,wp,a2, wwaa)
c subroutine wwza_anomal4(wm,z,wp,a, wwza)
c subroutine wwzz_anomal4(wm,z1,wp,z2, wwzz)
c subroutine wwww_anomal4(wm1,wp1,wm2,wp2, wwww)
c
c subroutine zzzz_anomal4(z1,z2,z3,z4, zzzz)
c subroutine azzz_anomal4(a,z1,z2,z3, azzz)
c subroutine aazz_anomal4(a1,a2,z1,z2, aazz)
c subroutine aaaz_anomal4(a1,a2,a3,z, aaaz)
c subroutine aaaa_anomal4(a1,a2,a3,a4, aaaa)
c
c subroutine anomal_formfactor(p1,p2,p3,p4,p5)


CCCCCCCCCCCCCCC
C WWAA-vertex C
CCCCCCCCCCCCCCC

      subroutine wwaa_anomal4(wm,a1,wp,a2, wwaa)
c
c This subroutine computes the wwaa-vertex with anomalous couplings
c where all the bosons are defined by the flowing-OUT quantum number.
c     
      implicit none
      double complex wp(6),wm(6),a1(6),a2(6),wwaa 
      double complex dotcc, dotrc, epsccrr, epscccr, epscccc
      double precision pwp(0:3),pwm(0:3),pa1(0:3),pa2(0:3), dotrr
      double precision cos2w, g1, cosw, sinw

      external dotcc, dotrc, dotrr, epsccrr, epscccr, epscccc

#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"

      cos2w = 1.0d0 - sin2w
      cosw = sqrt(cos2w)
      sinw = sqrt(sin2w)
      g1 = ee / cosw

      pwm(0) = dble( wm(5))
      pwm(1) = dble( wm(6))
      pwm(2) = dimag(wm(6))
      pwm(3) = dimag(wm(5))
      pwp(0) = dble( wp(5))
      pwp(1) = dble( wp(6))
      pwp(2) = dimag(wp(6))
      pwp(3) = dimag(wp(5))
      pa1(0) = dble( a1(5))
      pa1(1) = dble( a1(6))
      pa1(2) = dimag(a1(6))
      pa1(3) = dimag(a1(5))
      pa2(0) = dble( a2(5))
      pa2(1) = dble( a2(6))
      pa2(2) = dimag(a2(6))
      pa2(3) = dimag(a2(5))

      wwaa = gwwa**2 * (dotcc(wm,a1)*dotcc(wp,a2) + dotcc(wm,a2)*dotcc(wp,a1) - 2.0d0*dotcc(wm,wp)*dotcc(a1,a2))

      if (ft0 .ne. 0) then
        wwaa = wwaa + ft0 * 8. * gw**4 * sin2w * ( dotrc(pwp,wm)*dotrc(pwm,wp) - dotcc(wp,wm)*dotrr(pwp,pwm) ) * 
     &   ( dotrc(pa1,a2)*dotrc(pa2,a1) - dotcc(a1,a2)*dotrr(pa1,pa2) )
      endif

      if (ft1 .ne. 0) then
        wwaa = wwaa + ft1 * 4. * gw**4 * sin2w *
     &   ( ( dotrc(pwp,a1)*dotrc(pa1,wp) - dotcc(a1,wp)*dotrr(pa1,pwp) ) 
     &      * ( dotrc(pwm,a2)*dotrc(pa2,wm) - dotcc(a2,wm)*dotrr(pa2,pwm) )
     &   + ( dotrc(pwm,a1)*dotrc(pa1,wm) - dotcc(a1,wm)*dotrr(pa1,pwm) ) 
     &      * ( dotrc(pwp,a2)*dotrc(pa2,wp) - dotcc(a2,wp)*dotrr(pa2,pwp) ) )
      endif

      if (ft2 .ne. 0) then
        wwaa = wwaa + ft2 * gw**4 * sin2w *
     &    ( dotcc(a1,a2)*dotcc(wm,wp) * ( dotrr(pa1,pwp)*dotrr(pa2,pwm) + dotrr(pa1,pwm)*dotrr(pa2,pwp) )
     &    + dotcc(wm,wp) * ( - dotrc(pa2,a1) * ( dotrc(pwp,a2)*dotrr(pa1,pwm) + dotrc(pwm,a2)*dotrr(pa1,pwp) )
     &                       + dotrc(pwp,a1) * ( dotrc(pwm,a2)*dotrr(pa1,pa2) - dotrc(pa1,a2)*dotrr(pa2,pwm) )
     &                       + dotrc(pwm,a1) * ( dotrc(pwp,a2)*dotrr(pa1,pa2) - dotrc(pa1,a2)*dotrr(pa2,pwp) ) )
     &    + dotcc(a1,wp)*dotcc(a2,wm)*dotrr(pa1,pa2)*dotrr(pwm,pwp)
     &    + dotcc(a1,wm)*dotcc(a2,wp)*dotrr(pa1,pa2)*dotrr(pwm,pwp)
     &    - dotcc(a1,a2) * ( + dotrc(pwp,wm) * ( dotrc(pa2,wp)*dotrr(pa1,pwm) + dotrc(pa1,wp)*dotrr(pa2,pwm) )
     &                       + dotrc(pwm,wp) * ( dotrc(pa2,wm)*dotrr(pa1,pwp) + dotrc(pa1,wm)*dotrr(pa2,pwp) )
     &                       - dotrr(pwm,pwp) * ( dotrc(pa1,wp)*dotrc(pa2,wm) + dotrc(pa1,wm)*dotrc(pa2,wp) ) )
     &    - dotcc(a2,wp) * ( dotrc(pwm,a1)*dotrc(pwp,wm)*dotrr(pa1,pa2) 
     &                     + dotrc(pa2,a1) * ( dotrc(pa1,wm)*dotrr(pwm,pwp) - dotrc(pwp,wm)*dotrr(pa1,pwm) ) )
     &    - dotcc(a1,wp) * ( dotrc(pwm,a2)*dotrc(pwp,wm)*dotrr(pa1,pa2) 
     &                     + dotrc(pa1,a2) * ( dotrc(pa2,wm)*dotrr(pwm,pwp) - dotrc(pwp,wm)*dotrr(pa2,pwm) ) )
     &    - dotcc(a2,wm) * ( dotrc(pwp,a1)*dotrc(pwm,wp)*dotrr(pa1,pa2) 
     &                     + dotrc(pa2,a1) * ( dotrc(pa1,wp)*dotrr(pwm,pwp) - dotrc(pwm,wp)*dotrr(pa1,pwp) ) )
     &    - dotcc(a1,wm) * ( dotrc(pwp,a2)*dotrc(pwm,wp)*dotrr(pa1,pa2) 
     &                     + dotrc(pa1,a2) * ( dotrc(pa2,wp)*dotrr(pwm,pwp) - dotrc(pwm,wp)*dotrr(pa2,pwp) ) )
     &   + dotrc(pa2,a1)*dotrc(pwm,a2)*dotrc(pa1,wp)*dotrc(pwp,wm)
     &   + dotrc(pa2,wp)*dotrc(pwm,a1)*dotrc(pa1,a2)*dotrc(pwp,wm)
     &   + dotrc(pa2,a1)*dotrc(pwm,wp)*dotrc(pa1,wm)*dotrc(pwp,a2)
     &   + dotrc(pa2,wm)*dotrc(pwm,wp)*dotrc(pa1,a2)*dotrc(pwp,a1) )
      endif

      if (ft5 .ne. 0) then
        wwaa = wwaa + ft5 * 2. * gw**4 * sin2w * ( dotrc(pwp,wm)*dotrc(pwm,wp) - dotcc(wp,wm)*dotrr(pwp,pwm) ) * 
     &   ( dotrc(pa1,a2)*dotrc(pa2,a1) - dotcc(a1,a2)*dotrr(pa1,pa2) )
      endif

      if (ft6 .ne. 0) then
        wwaa = wwaa + ft6 * gw**4 * sin2w *
     &    ( ( dotrc(pwp,a1)*dotrc(pa1,wp) - dotcc(a1,wp)*dotrr(pa1,pwp) ) 
     &      * ( dotrc(pwm,a2)*dotrc(pa2,wm) - dotcc(a2,wm)*dotrr(pa2,pwm) )
     &    + ( dotrc(pwm,a1)*dotrc(pa1,wm) - dotcc(a1,wm)*dotrr(pa1,pwm) ) 
     &      * ( dotrc(pwp,a2)*dotrc(pa2,wp) - dotcc(a2,wp)*dotrr(pa2,pwp) ) )
      endif

      if (ft7 .ne. 0) then
        wwaa = wwaa + ft7 / 4. * gw**4 * sin2w *
     &    ( dotcc(a1,a2)*dotcc(wm,wp) * ( dotrr(pa1,pwp)*dotrr(pa2,pwm) + dotrr(pa1,pwm)*dotrr(pa2,pwp) )
     &    + dotcc(wm,wp) * ( - dotrc(pa2,a1) * ( dotrc(pwp,a2)*dotrr(pa1,pwm) + dotrc(pwm,a2)*dotrr(pa1,pwp) )
     &                       + dotrc(pwp,a1) * ( dotrc(pwm,a2)*dotrr(pa1,pa2) - dotrc(pa1,a2)*dotrr(pa2,pwm) )
     &                       + dotrc(pwm,a1) * ( dotrc(pwp,a2)*dotrr(pa1,pa2) - dotrc(pa1,a2)*dotrr(pa2,pwp) ) )
     &    + dotcc(a1,wp)*dotcc(a2,wm)*dotrr(pa1,pa2)*dotrr(pwm,pwp)
     &    + dotcc(a1,wm)*dotcc(a2,wp)*dotrr(pa1,pa2)*dotrr(pwm,pwp)
     &    - dotcc(a1,a2) * ( + dotrc(pwp,wm) * ( dotrc(pa2,wp)*dotrr(pa1,pwm) + dotrc(pa1,wp)*dotrr(pa2,pwm) )
     &                       + dotrc(pwm,wp) * ( dotrc(pa2,wm)*dotrr(pa1,pwp) + dotrc(pa1,wm)*dotrr(pa2,pwp) )
     &                       - dotrr(pwm,pwp) * ( dotrc(pa1,wp)*dotrc(pa2,wm) + dotrc(pa1,wm)*dotrc(pa2,wp) ) )
     &    - dotcc(a2,wp) * ( dotrc(pwm,a1)*dotrc(pwp,wm)*dotrr(pa1,pa2)
     &                     + dotrc(pa2,a1) * (dotrc(pa1,wm)*dotrr(pwm,pwp) - dotrc(pwp,wm)*dotrr(pa1,pwm) ) )
     &    - dotcc(a1,wp) * ( dotrc(pwm,a2)*dotrc(pwp,wm)*dotrr(pa1,pa2)
     &                     + dotrc(pa1,a2) * (dotrc(pa2,wm)*dotrr(pwm,pwp) - dotrc(pwp,wm)*dotrr(pa2,pwm) ) )
     &    - dotcc(a2,wm) * ( dotrc(pwp,a1)*dotrc(pwm,wp)*dotrr(pa1,pa2)
     &                     + dotrc(pa2,a1) * (dotrc(pa1,wp)*dotrr(pwm,pwp) - dotrc(pwm,wp)*dotrr(pa1,pwp) ) )
     &    - dotcc(a1,wm) * ( dotrc(pwp,a2)*dotrc(pwm,wp)*dotrr(pa1,pa2)
     &                     + dotrc(pa1,a2) * (dotrc(pa2,wp)*dotrr(pwm,pwp) - dotrc(pwm,wp)*dotrr(pa2,pwp) ) )
     &    + dotrc(pa2,a1)*dotrc(pwm,a2)*dotrc(pa1,wp)*dotrc(pwp,wm)
     &    + dotrc(pwm,a1)*dotrc(pa1,a2)*dotrc(pa2,wp)*dotrc(pwp,wm)
     &    + dotrc(pa2,a1)*dotrc(pwp,a2)*dotrc(pa1,wm)*dotrc(pwm,wp)
     &    + dotrc(pwp,a1)*dotrc(pa1,a2)*dotrc(pa2,wm)*dotrc(pwm,wp) )
      endif

      if (fm0 .ne. 0) then
        wwaa = wwaa - fm0 * 2. * gw**2 * wmass**2 * sin2w *
     &      dotcc(wm,wp) * ( dotrc(pa2,a1)*dotrc(pa1,a2) - dotcc(a1,a2)*dotrr(pa1,pa2) )
      endif

      if (fm1 .ne. 0) then
        wwaa = wwaa + fm1 / 2. * gw**2 * wmass**2 * sin2w *
     &      ( dotrc(pa2,a1) * ( dotrc(pa1,wm)*dotcc(a2,wp) + dotrc(pa1,wp)*dotcc(a2,wm) )
     &      + dotrc(pa1,a2) * ( dotrc(pa2,wm)*dotcc(a1,wp) + dotrc(pa2,wp)*dotcc(a1,wm) )
     &      - dotcc(a1,a2) * ( dotrc(pa1,wp)*dotrc(pa2,wm) + dotrc(pa1,wm)*dotrc(pa2,wp) )
     &      - dotrr(pa1,pa2) * (dotcc(a1,wp)*dotcc(a2,wm) + dotcc(a1,wm)*dotcc(a2,wp) ) )
      endif

      if (fm2 .ne. 0) then
        wwaa = wwaa - fm2 * gw**2 * wmass**2 * sin2w * dotcc(wp,wm) * 
     &   ( dotrc(pa1,a2) * dotrc(pa2,a1) - dotcc(a1,a2) * dotrr(pa1,pa2) )
      endif

      if (fm3 .ne. 0) then
        wwaa = wwaa + fm3 / 4. * gw**2 * wmass**2 * sin2w *
     &   ( dotrc(pa2,a1) * ( dotrc(pa1,wm)*dotcc(a2,wp) + dotrc(pa1,wp)*dotcc(a2,wm) )
     &   + dotrc(pa1,a2) * ( dotrc(pa2,wm)*dotcc(a1,wp) + dotrc(pa2,wp)*dotcc(a1,wm) )
     &   - dotcc(a1,a2) * ( dotrc(pa1,wp)*dotrc(pa2,wm) + dotrc(pa1,wm)*dotrc(pa2,wp) )
     &   - dotrr(pa1,pa2) * ( dotcc(a1,wp)*dotcc(a2,wm) + dotcc(a1,wm)*dotcc(a2,wp) ) )
      endif

      if (fm4 .ne. 0) then
        wwaa = wwaa - fm4 * gw**2 * wmass**2 * sin2w *
     &   dotcc(wm,wp) * ( dotrc(pa2,a1)*dotrc(pa1,a2) - dotcc(a1,a2)*dotrr(pa1,pa2) )
      endif

      if (fm5 .ne. 0) then
        wwaa = wwaa + fm5 / 4. * gw**2 * wmass**2 * sin2w *
     &   ( - dotrc(pa2,a1) * ( dotrc(pa1,wm)*dotcc(a2,wp) + dotrc(pa1,wp)*dotcc(a2,wm) )
     &     - dotrc(pa1,a2) * ( dotrc(pa2,wm)*dotcc(a1,wp) + dotrc(pa2,wp)*dotcc(a1,wm) )
     &     + dotcc(a1,a2) * ( dotrc(pa1,wp)*dotrc(pa2,wm) + dotrc(pa1,wm)*dotrc(pa2,wp) )
     &     + dotrr(pa1,pa2) * ( dotcc(a1,wp)*dotcc(a2,wm) + dotcc(a1,wm)*dotcc(a2,wp) ) )
      endif

      if (fm6 .ne. 0) then
        wwaa = wwaa - fm6 * gw**2 * wmass**2 * sin2w *
     &      dotcc(wp,wm) * ( dotrc(pa1,a2) * dotrc(pa2,a1) - dotcc(a1,a2) * dotrr(pa1,pa2) )
      endif

      if (fm7 .ne. 0) then
        wwaa = wwaa + fm7 / 4. * gw**2 * wmass**2 * sin2w *
     &   ( - dotrc(pa2,a1) * ( dotrc(pa1,wm)*dotcc(a2,wp) + dotrc(pa1,wp)*dotcc(a2,wm) )
     &     - dotrc(pa1,a2) * ( dotrc(pa2,wm)*dotcc(a1,wp) + dotrc(pa2,wp)*dotcc(a1,wm) )
     &     + dotcc(a1,a2) * ( dotrc(pa1,wp)*dotrc(pa2,wm) + dotrc(pa1,wm)*dotrc(pa2,wp) )
     &     + dotrr(pa1,pa2) * ( dotcc(a1,wp)*dotcc(a2,wm) + dotcc(a1,wm)*dotcc(a2,wp) ) )
      endif

      if (fwww .ne. 0) then
        wwaa = wwaa - fwww * 3./2. * gw**4 * sin2w *
     &      ( (dotrc(pa2,a1) * (dotrc(pwm,a2) + dotrc(pwp,a2)) + dotrc(pa1,a2) * (dotrc(pwm,a1) + dotrc(pwp,a1)) ) * dotcc(wp,wm)
     &      - dotcc(a1,a2) * dotcc(wp,wm) * (dotrr(pwm,pa1) + dotrr(pwm,pa2) + dotrr(pwp,pa1) + dotrr(pwp,pa2))
     &      - dotcc(a2,wm) * ( dotrc(pa2,a1)*dotrc(pwm,wp) + dotrc(pwm,a1)*(dotrc(pa1,wp)-dotrc(pa2,wp)) 
     &                       + dotrc(pwp,a1)*dotrc(pa2,wp) )
     &      - dotcc(a1,wm) * ( dotrc(pa1,a2)*dotrc(pwm,wp) + dotrc(pwp,a2)*dotrc(pa1,wp)
     &                       + dotrc(pwm,a2) * (dotrc(pa2,wp)-dotrc(pa1,wp)) )
     &      - dotcc(a2,wp) * ( dotrc(pa2,a1)*dotrc(pwp,wm) + dotrc(pwp,a1) * (dotrc(pa1,wm)-dotrc(pa2,wm))
     &                       + dotrc(pwm,a1)*dotrc(pa2,wm) )
     &      - dotcc(a1,wp) * ( dotrc(pa1,a2)*dotrc(pwp,wm) + dotrc(pwm,a2)*dotrc(pa1,wm)
     &                       + dotrc(pwp,a2) * (dotrc(pa2,wm) - dotrc(pa1,wm)) )
     &      + dotcc(a1,a2) * ( (dotrc(pa1,wp)+dotrc(pa2,wp))*dotrc(pwp,wm) + dotrc(pwm,wp)*(dotrc(pa1,wm)+dotrc(pa2,wm)) )
     &      + dotcc(a1,wm) * dotcc(a2,wp) * ( dotrr(pwm,pa2) + dotrr(pwp,pa1) )
     &      + dotcc(a1,wp) * dotcc(a2,wm) * ( dotrr(pwm,pa1) + dotrr(pwp,pa2) ) )
      endif

      if (fdw .ne. 0) then
        wwaa = wwaa - fdw * gw**4 * 2. * sin2w * (
     &      + dotcc(a1,a2) * ( - dotrc(pwm,wp) * ( dotrc(pa1,wm) + dotrc(pa2,wm) - 2.*dotrc(pwp,wm) )
     &                         + dotrc(pa1,wp) * ( 2.*dotrc(pa1,wm) - 2.*dotrc(pa2,wm) - dotrc(pwp,wm) )
     &                         - dotrc(pa2,wp) * ( 2.*dotrc(pa1,wm) - 2.*dotrc(pa2,wm) + dotrc(pwp,wm) ) )
     &      + dotcc(wp,a2) * ( + dotrc(pa1,wm) * ( dotrc(pa2,a1) - dotrc(pwm,a1) )
     &                         + dotrc(pwp,wm) * ( 2.*dotrc(pa2,a1) + dotrc(pwm,a1) - dotrc(pwp,a1) )
     &                         - dotrc(pa2,wm) * dotrc(pa2,a1) )
     &      + dotcc(wm,a2) * ( + dotrc(pa1,wp) * ( dotrc(pa2,a1) - dotrc(pwp,a1) )
     &                         + dotrc(pwm,wp) * ( 2.*dotrc(pa2,a1) - dotrc(pwm,a1) + dotrc(pwp,a1) )
     &                         - dotrc(pa2,wp) * dotrc(pa2,a1) )
     &      - dotcc(wm,a1) * ( + dotrc(pa2,wp) * ( dotrc(pwp,a2) - dotrc(pa1,a2) )
     &                         + dotrc(pwm,wp) * ( - 2.*dotrc(pa1,a2) + dotrc(pwm,a2) - dotrc(pwp,a2) )
     &                         + dotrc(pa1,wp) * dotrc(pa1,a2) )
     &      - dotcc(wp,a1) * ( + dotrc(pa2,wm) * ( dotrc(pwm,a2) - dotrc(pa1,a2) )
     &                         + dotrc(pwp,wm) * ( - 2.*dotrc(pa1,a2) - dotrc(pwm,a2) + dotrc(pwp,a2) )
     &                         + dotrc(pa1,wm) * dotrc(pa1,a2) )
     &      + dotcc(wp,wm) * ( + dotrc(pa2,a1) * ( 2.*dotrc(pa1,a2) - dotrc(pwm,a2) - dotrc(pwp,a2) )
     &                         - dotrc(pwp,a1) * ( dotrc(pa1,a2) + 2.*dotrc(pwm,a2) - 2.*dotrc(pwp,a2) )
     &                         - dotrc(pwm,a1) * ( dotrc(pa1,a2) - 2.*dotrc(pwm,a2) + 2.*dotrc(pwp,a2) ) )
     &      - dotcc(wp,wm)*dotcc(a1,a2) * ( 4.*dotrr(pa1,pa2) + dotrr(pa1,pwm) 
     &                                      + dotrr(pa1,pwp) + dotrr(pa2,pwm) + dotrr(pa2,pwp) + 4.*dotrr(pwm,pwp) )
     &      + dotcc(wp,a1)*dotcc(wm,a2) * ( dotrr(pa1,pa2) + 2.*dotrr(pa1,pwp) + 2.*dotrr(pa2,pwm) + dotrr(pwm,pwp) )
     &      + dotcc(wp,a2)*dotcc(wm,a1) * ( dotrr(pa1,pa2) + 2.*dotrr(pa1,pwm) + 2.*dotrr(pa2,pwp) + dotrr(pwm,pwp) ) )
      endif

      if (fwwwtilde .ne. 0) then
        wwaa = wwaa - fwwwtilde / 2. * gw**4 * sin2w *
     &      ( dotcc(wp,wm) * ( epsccrr(a1,a2,pa1,pwm) + epsccrr(a1,a2,pa1,pwp) - epsccrr(a1,a2,pa2,pwm) - epsccrr(a1,a2,pa2,pwp) )
     &      + dotcc(wp,a1) * ( epsccrr(wm,a2,pa1,pwm) + epsccrr(wm,a2,pa1,pwp) + epsccrr(wm,a2,pa2,pwp) )
     &      + dotcc(wp,a2) * ( epsccrr(wm,a1,pa1,pwp) + epsccrr(wm,a1,pa2,pwm) + epsccrr(wm,a1,pa2,pwp) )
     &      + dotcc(wm,a1) * ( epsccrr(wp,a2,pa1,pwm) + epsccrr(wp,a2,pa1,pwp) + epsccrr(wp,a2,pa2,pwm) )
     &      + dotcc(wm,a2) * ( epsccrr(wp,a1,pa1,pwm) + epsccrr(wp,a1,pa2,pwm) + epsccrr(wp,a1,pa2,pwp) )
     &      + dotcc(a1,a2) * ( epsccrr(wp,wm,pa1,pwm) - epsccrr(wp,wm,pa1,pwp) + epsccrr(wp,wm,pa2,pwm) - epsccrr(wp,wm,pa2,pwp) )
     &      + ( dotrc(pa1,wp) - dotrc(pa2,wp) ) * ( epscccr(wm,a1,a2,pwm) + epscccr(wm,a1,a2,pwp) )
     &      + ( dotrc(pa1,wm) - dotrc(pa2,wm) ) * ( epscccr(wp,a1,a2,pwm) + epscccr(wp,a1,a2,pwp) )
     &      - ( dotrc(pwm,a1) - dotrc(pwp,a1) ) * ( epscccr(wp,wm,a2,pa1) + epscccr(wp,wm,a2,pa2) )
     &      - ( dotrc(pwm,a2) - dotrc(pwp,a2) ) * ( epscccr(wp,wm,a1,pa1) + epscccr(wp,wm,a1,pa2) )
     &      + epscccc(wp,wm,a1,a2) * ( - dotrr(pa1,pwm) + dotrr(pa1,pwp) + dotrr(pa2,pwm) - dotrr(pa2,pwp) )
     &      + dotrc(pwm,wp) * ( epscccr(wm,a1,a2,pa1) - epscccr(wm,a1,a2,pa2) )
     &      + dotrc(pwp,wm) * ( epscccr(wp,a1,a2,pa1) - epscccr(wp,a1,a2,pa2) )
     &      + dotrc(pa1,a2) * ( epscccr(wp,wm,a1,pwp) - epscccr(wp,wm,a1,pwm) )
     &      + dotrc(pa2,a1) * ( epscccr(wp,wm,a2,pwp) - epscccr(wp,wm,a2,pwm) ) )
      endif

      if (fdwtilde .ne. 0) then
        wwaa = wwaa - fdwtilde * 2. * gw**4 * sin2w * (
     &      + 2.*dotcc(wp,wm)*epsccrr(a1,a2,pa1,pa2) + dotcc(wp,a2)*epsccrr(wm,a1,pa1,pwm)
     &      + dotcc(wm,a2)*epsccrr(wp,a1,pa1,pwp) + dotcc(wp,a1)*epsccrr(wm,a2,pa2,pwm)
     &      + dotcc(wm,a1)*epsccrr(wp,a2,pa2,pwp) - 2.*dotcc(a1,a2)*epsccrr(wp,wm,pwm,pwp)
     &      + epscccc(wp,wm,a1,a2) * ( dotrr(pa1,pwm) - dotrr(pa1,pwp) - dotrr(pa2,pwm) + dotrr(pa2,pwp) )
     &      + ( dotrc(pa1,wp) - dotrc(pa2,wp) - dotrc(pwm,wp) ) * epscccr(wm,a1,a2,pa1)
     &      + ( dotrc(pa1,wp) - dotrc(pa2,wp) + dotrc(pwm,wp) ) * epscccr(wm,a1,a2,pa2)
     &      + ( dotrc(pa1,wm) - dotrc(pa2,wm) - dotrc(pwp,wm) ) * epscccr(wp,a1,a2,pa1)
     &      + ( dotrc(pa1,wm) - dotrc(pa2,wm) + dotrc(pwp,wm) ) * epscccr(wp,a1,a2,pa2)
     &      - ( dotrc(pa1,a2) + dotrc(pwm,a2) - dotrc(pwp,a2) ) * epscccr(wp,wm,a1,pwp)
     &      + ( dotrc(pa1,a2) - dotrc(pwm,a2) + dotrc(pwp,a2) ) * epscccr(wp,wm,a1,pwm)
     &      - ( dotrc(pa2,a1) + dotrc(pwm,a1) - dotrc(pwp,a1) ) * epscccr(wp,wm,a2,pwp)
     &      + ( dotrc(pa2,a1) - dotrc(pwm,a1) + dotrc(pwp,a1) ) * epscccr(wp,wm,a2,pwm) )
      endif

      return
      end


*********************************************************************************************************
CCCCCCCCCCCCCCC
C WWZA-vertex C
CCCCCCCCCCCCCCC

      subroutine wwza_anomal4(wm,z,wp,a, wwza)
c
c This subroutine computes the wwza-vertex with anomalous couplings
c where all the bosons are defined by the flowing-OUT quantum number.
c     
      implicit none
      double complex wp(6),wm(6),z(6),a(6),wwza
      double complex dotcc, dotrc, epsccrr, epscccr, epscccc
      double precision pwp(0:3),pwm(0:3),pz(0:3),pa(0:3), dotrr
      double precision cos2w, g1, cosw, sinw

      external dotcc, dotrc, dotrr, epsccrr, epscccr, epscccc

#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"  

      cos2w = 1.0d0 - sin2w
      cosw = sqrt(cos2w)
      sinw = sqrt(sin2w)
      g1 = ee / cosw

      pwm(0) = dble( wm(5))
      pwm(1) = dble( wm(6))
      pwm(2) = dimag(wm(6))
      pwm(3) = dimag(wm(5))
      pwp(0) = dble( wp(5))
      pwp(1) = dble( wp(6))
      pwp(2) = dimag(wp(6))
      pwp(3) = dimag(wp(5))
      pz(0)  = dble( z(5))
      pz(1)  = dble( z(6))
      pz(2)  = dimag(z(6))
      pz(3)  = dimag(z(5))
      pa(0)  = dble( a(5))
      pa(1)  = dble( a(6))
      pa(2)  = dimag(a(6))
      pa(3)  = dimag(a(5))

      wwza = gwwz * gwwa * (dotcc(wm,z)*dotcc(wp,a) + dotcc(wm,a)*dotcc(wp,z) - 2.0d0*dotcc(wm,wp)*dotcc(z,a))

      if (ft0 .ne. 0) then
        wwza = wwza + ft0 * 8. * gw**4 * sqrt(cos2w * sin2w) * 
     &   ( dotrc(pwp,wm)*dotrc(pwm,wp) - dotcc(wp,wm)*dotrr(pwp,pwm) ) * ( dotrc(pz,a)*dotrc(pa,z) - dotcc(z,a)*dotrr(pz,pa) )
      endif

      if (ft1 .ne. 0) then
        wwza = wwza + ft1 * 4. * gw**4 * sqrt(cos2w * sin2w) * 
     &   ( ( dotrc(pwp,a)*dotrc(pa,wp) - dotcc(a,wp)*dotrr(pa,pwp) ) * ( dotrc(pwm,z)*dotrc(pz,wm) - dotcc(wm,z)*dotrr(pwm,pz) )
     &   + ( dotrc(pwm,a)*dotrc(pa,wm) - dotcc(a,wm)*dotrr(pa,pwm) ) * ( dotrc(pwp,z)*dotrc(pz,wp) - dotcc(wp,z)*dotrr(pwp,pz) ) )
      endif

      if (ft2 .ne. 0) then
        wwza = wwza + ft2 * gw**4 * sinw * cosw *
     &    ( dotcc(a,wp)*dotcc(wm,z)*dotrr(pa,pz)*dotrr(pwm,pwp)
     &    + dotcc(a,wm)*dotcc(wp,z)*dotrr(pa,pz)*dotrr(pwm,pwp)
     &    + dotcc(wp,z) * ( dotrc(pz,a) * ( dotrc(pwp,wm)*dotrr(pa,pwm) - dotrc(pa,wm)*dotrr(pwm,pwp) )
     &                    - dotrc(pwm,a)*dotrc(pwp,wm)*dotrr(pa,pz) )
     &    + dotcc(wm,z) * ( dotrc(pz,a) * ( dotrc(pwm,wp)*dotrr(pa,pwp) - dotrc(pa,wp)*dotrr(pwm,pwp) )
     &                    - dotrc(pwp,a)*dotrc(pwm,wp)*dotrr(pa,pz) )
     &    - dotcc(a,wp) * ( dotrc(pwp,wm) * ( dotrc(pwm,z)*dotrr(pa,pz) - dotrc(pa,z)*dotrr(pwm,pz) )
     &                    + dotrc(pa,z)*dotrc(pz,wm)*dotrr(pwm,pwp) )
     &    + dotcc(a,z)*dotcc(wm,wp) * ( dotrr(pa,pwp)*dotrr(pwm,pz) + dotrr(pa,pwm)*dotrr(pwp,pz) )
     &    - dotcc(a,z) * ( dotrc(pz,wm) * ( dotrc(pwm,wp)*dotrr(pa,pwp) - dotrc(pa,wp)*dotrr(pwm,pwp) )
     &                   + dotrc(pwp,wm) * ( dotrc(pz,wp)*dotrr(pa,pwm) + dotrc(pa,wp)*dotrr(pwm,pz) )
     &                   + dotrc(pa,wm) * ( dotrc(pwm,wp)*dotrr(pwp,pz) - dotrc(pz,wp)*dotrr(pwm,pwp) ) )
     &    - dotcc(a,wm) * ( dotrc(pa,z)*dotrc(pz,wp)*dotrr(pwm,pwp) 
     &                    + dotrc(pwm,wp) * ( dotrc(pwp,z)*dotrr(pa,pz) - dotrc(pa,z)*dotrr(pwp,pz) ) )
     &    - dotcc(wm,wp) * ( dotrc(pz,a) * ( dotrc(pwp,z)*dotrr(pa,pwm) + dotrc(pwm,z)*dotrr(pa,pwp) )
     &                   + dotrc(pwp,a) * ( dotrc(pa,z)*dotrr(pwm,pz) - dotrc(pwm,z)*dotrr(pa,pz) )
     &                   + dotrc(pwm,a) * ( dotrc(pa,z)*dotrr(pwp,pz) - dotrc(pwp,z)*dotrr(pa,pz) ) )
     &    + dotrc(pwp,a)*dotrc(pa,z)*dotrc(pwm,wp)*dotrc(pz,wm)
     &    + dotrc(pwp,wm)*dotrc(pa,z)*dotrc(pwm,a)*dotrc(pz,wp)
     &    + dotrc(pwp,wm)*dotrc(pa,wp)*dotrc(pwm,z)*dotrc(pz,a)
     &    + dotrc(pwp,z)*dotrc(pa,wm)*dotrc(pwm,wp)*dotrc(pz,a) )
      endif

      if (ft5 .ne. 0) then
        wwza = wwza - ft5 * 2. * gw**4 * sin2w * sqrt(sin2w / cos2w) *
     &   ( dotrc(pwp,wm)*dotrc(pwm,wp) - dotcc(wp,wm)*dotrr(pwp,pwm) ) * ( dotrc(pz,a)*dotrc(pa,z) - dotcc(z,a)*dotrr(pz,pa) )
      endif

      if (ft6 .ne. 0) then
        wwza = wwza - ft6 * gw**4 * sin2w * sqrt(sin2w / cos2w) *
     &   ( ( dotrc(pwp,a)*dotrc(pa,wp) - dotcc(a,wp)*dotrr(pa,pwp) ) * ( dotrc(pwm,z)*dotrc(pz,wm) - dotcc(z,wm)*dotrr(pz,pwm) )
     &   + ( dotrc(pwm,a)*dotrc(pa,wm) - dotcc(a,wm)*dotrr(pa,pwm) ) * ( dotrc(pwp,z)*dotrc(pz,wp) - dotcc(z,wp)*dotrr(pz,pwp) ) )
      endif

      if (ft7 .ne. 0) then
        wwza = wwza - ft7 / 4. * gw**4 * sin2w * sqrt(sin2w / cos2w) *
     &    ( dotcc(a,wp)*dotcc(wm,z)*dotrr(pa,pz)*dotrr(pwm,pwp)
     &    + dotcc(a,wm)*dotcc(wp,z)*dotrr(pa,pz)*dotrr(pwm,pwp)
     &    + dotcc(wp,z) * ( dotrc(pz,a) * ( dotrc(pwp,wm)*dotrr(pa,pwm) - dotrc(pa,wm)*dotrr(pwm,pwp) )
     &                    - dotrc(pwm,a)*dotrc(pwp,wm)*dotrr(pa,pz) )
     &    + dotcc(wm,z) * ( dotrc(pz,a) * ( dotrc(pwm,wp)*dotrr(pa,pwp) - dotrc(pa,wp)*dotrr(pwm,pwp) )
     &                    - dotrc(pwp,a)*dotrc(pwm,wp)*dotrr(pa,pz) )
     &    - dotcc(wp,a) * ( dotrc(pwp,wm) * ( dotrc(pwm,z)*dotrr(pa,pz) - dotrc(pa,z)*dotrr(pwm,pz) )
     &                    + dotrc(pa,z)*dotrc(pz,wm)*dotrr(pwm,pwp) )
     &    + dotcc(a,z)*dotcc(wm,wp) * ( dotrr(pa,pwp)*dotrr(pwm,pz) + dotrr(pa,pwm)*dotrr(pwp,pz) )
     &    - dotcc(a,z) * ( dotrc(pz,wm) * ( dotrc(pwm,wp)*dotrr(pa,pwp) - dotrc(pa,wp)*dotrr(pwm,pwp) )
     &                   + dotrc(pwp,wm) * ( dotrc(pz,wp)*dotrr(pa,pwm) + dotrc(pa,wp)*dotrr(pwm,pz) )
     &                   + dotrc(pa,wm) * ( dotrc(pwm,wp)*dotrr(pwp,pz) - dotrc(pz,wp)*dotrr(pwm,pwp) ) )
     &    - dotcc(wm,a) * ( dotrc(pwm,wp) * ( dotrc(pwp,z)*dotrr(pa,pz) - dotrc(pa,z)*dotrr(pwp,pz) )
     &                    + dotrc(pa,z)*dotrc(pz,wp)*dotrr(pwm,pwp) )
     &    - dotcc(wm,wp) * ( dotrc(pz,a) * ( dotrc(pwm,z)*dotrr(pa,pwp) + dotrc(pwp,z)*dotrr(pa,pwm) )
     &                   + dotrc(pwp,a) * ( - dotrc(pwm,z)*dotrr(pa,pz) + dotrc(pa,z)*dotrr(pwm,pz) )
     &                   + dotrc(pwm,a) * ( dotrc(pa,z)*dotrr(pwp,pz) - dotrc(pwp,z)*dotrr(pa,pz) ) )
     &    + dotrc(pwp,a)*dotrc(pa,z)*dotrc(pwm,wp)*dotrc(pz,wm)
     &    + dotrc(pwp,wm)*dotrc(pa,z)*dotrc(pwm,a)*dotrc(pz,wp)
     &    + dotrc(pwp,wm)*dotrc(pa,wp)*dotrc(pwm,z)*dotrc(pz,a)
     &    + dotrc(pwp,z)*dotrc(pa,wm)*dotrc(pwm,wp)*dotrc(pz,a) )
      endif

      if (fm0 .ne. 0) then
        wwza = wwza - fm0 * 2. * gw**2 * wmass**2 * sinw * cosw *
     &      dotcc(wm,wp) * ( dotrc(pz,a)*dotrc(pa,z) - dotcc(a,z)*dotrr(pa,pz) )
      endif

      if (fm1 .ne. 0) then
        wwza = wwza + fm1 / 2. * gw**2 * wmass**2 * sinw * cosw *
     &      ( dotrc(pz,a) * ( dotrc(pa,wm)*dotcc(wp,z) + dotrc(pa,wp)*dotcc(wm,z) )
     &      - dotcc(a,z) * ( dotrc(pa,wp)*dotrc(pz,wm) + dotrc(pa,wm)*dotrc(pz,wp) )
     &      + dotrc(pa,z) * ( dotrc(pz,wm)*dotcc(a,wp) + dotrc(pz,wp)*dotcc(a,wm) )
     &      - dotrr(pa,pz) * ( dotcc(a,wp)*dotcc(wm,z) + dotcc(a,wm)*dotcc(wp,z) ) )
      endif

      if (fm2 .ne. 0) then
        wwza = wwza + fm2 * gw**2 * wmass**2 * sin2w * sinw / cosw * dotcc(wp,wm) * 
     &   ( dotrc(pa,z) * dotrc(pz,a) - dotcc(a,z) * dotrr(pa,pz) )
      endif

      if (fm3 .ne. 0) then
        wwza = wwza - fm3 / 4. * gw**2 * wmass**2 * sin2w * sinw / cosw *
     &   ( dotrc(pz,a) * ( dotrc(pa,wm)*dotcc(wp,z) + dotrc(pa,wp)*dotcc(wm,z) )
     &   - dotcc(a,z) * ( dotrc(pa,wp)*dotrc(pz,wm) + dotrc(pa,wm)*dotrc(pz,wp) ) 
     &   + dotrc(pa,z) * ( dotrc(pz,wm)*dotcc(a,wp) + dotrc(pz,wp)*dotcc(a,wm) )
     &   - dotrr(pa,pz) * ( dotcc(a,wp)*dotcc(wm,z) + dotcc(a,wm)*dotcc(wp,z) ) )
      endif

      if (fm4 .ne. 0) then
        wwza = wwza + fm4 / 2. * gw**2 * wmass**2 * sinw / cosw *
     &   ( dotcc(wp,z) * ( dotrc(pwm,a)*dotrc(pa,wm) - dotcc(a,wm)*dotrr(pa,pwm) )
     &   + dotcc(wm,z) * ( dotrc(pwp,a)*dotrc(pa,wp) - dotcc(a,wp)*dotrr(pa,pwp) )
     &   + ( cos2w - sin2w ) * dotcc(wm,wp) * ( dotcc(a,z)*dotrr(pa,pz) - dotrc(pz,a)*dotrc(pa,z) ) )
      endif

      if (fm5 .ne. 0) then
        wwza = wwza + fm5 / 8. * gw**2 * wmass**2 * sinw / cosw *
     &   ( 2.*dotrc(pa,wm)*dotrc(pwm,wp)*dotcc(a,z) + 2.*dotrc(pa,z)*dotrc(pwm,a)*dotcc(wm,wp)
     &   - 2.*dotrc(pa,z)*dotrc(pwm,wp)*dotcc(a,wm) - 2.*dotcc(a,z)*dotcc(wm,wp)*dotrr(pa,pwm)
     &   + 2.*dotrc(pwp,a)*dotrc(pa,wm)*dotcc(wp,z) - 2.*dotrc(pa,wm)*dotrc(pwp,z)*dotcc(a,wp)
     &   + 2.*dotrc(pa,wp)*dotrc(pwp,z)*dotcc(a,wm) - 2.*dotcc(a,wm)*dotcc(wp,z)*dotrr(pa,pwp)
     &   + dotrc(pz,wm) * ( dotrc(pa,z)*dotcc(a,wp) - 2.*sin2w*dotrc(pa,wp)*dotcc(a,z) )
     &   + ( cos2w - sin2w ) * ( dotrc(pa,wm)*dotrc(pz,wp)*dotcc(a,z)
     &                         - dotrc(pa,z) * ( dotrc(pz,wm)*dotcc(a,wp) + dotrc(pz,wp)*dotcc(a,wm) )
     &                         + dotrr(pa,pz) * ( dotcc(a,wp)*dotcc(wm,z) + dotcc(a,wm)*dotcc(wp,z) ) )
     &   + dotrc(pz,a) * ( 2.*sin2w*dotrc(pa,wp)*dotcc(wm,z) - 2.*cos2w*dotrc(pa,wm)*dotcc(wp,z) )
     &   + dotrc(pa,wm)*dotrc(pz,wp)*dotcc(a,z) - dotrc(pa,z)*dotrc(pz,wp)*dotcc(a,wm)
     &   + dotcc(a,wm)*dotcc(wp,z)*dotrr(pa,pz) - dotcc(a,wp)*dotcc(wm,z)*dotrr(pa,pz) )
      endif

      if (fm6 .ne. 0) then
        wwza = wwza - fm6 * gw**2 * wmass**2 * sinw * cosw *
     &    dotcc(wm,wp) * ( dotrc(pz,a)*dotrc(pa,z) - dotcc(a,z)*dotrr(pa,pz) ) 
      endif

      if (fm7 .ne. 0) then
        wwza = wwza + fm7 / 8. * gw**2 * wmass**2 * sinw / cosw *
     &   ( 2.*dotrc(pa,wm)*dotrc(pwm,wp)*dotcc(a,z) - 2.*dotrc(pa,z)*dotrc(pwm,wp)*dotcc(a,wm)
     &   + 2.*dotrc(pwm,a) * ( dotrc(pa,z)*dotcc(wm,wp) - dotrc(pa,wp)*dotcc(wm,z) )
     &   - 2.*dotrc(pa,wm)*dotrc(pwm,z)*dotcc(a,wp) + 2.*dotrc(pa,wp)*dotrc(pwm,z)*dotcc(a,wm)
     &   - 2.*dotcc(a,z)*dotcc(wm,wp)*dotrr(pa,pwm) + 2.*dotcc(a,wp)*dotcc(wm,z)*dotrr(pa,pwm)
     &   + 2.*dotrc(pa,wp)*dotrc(pwp,wm)*dotcc(a,z) - 2.*dotrc(pa,z)*dotrc(pwp,wm)*dotcc(a,wp)
     &   + 2.*dotrc(pwp,a) * ( dotrc(pa,z)*dotcc(wm,wp) - dotrc(pa,wm)*dotcc(wp,z) )
     &   + 2.*dotrc(pa,wm)*dotrc(pwp,z)*dotcc(a,wp) - 2.*dotrc(pa,wp)*dotrc(pwp,z)*dotcc(a,wm)
     &   - 2.*dotcc(a,z)*dotcc(wm,wp)*dotrr(pa,pwp) + 2.*dotcc(a,wm)*dotcc(wp,z)*dotrr(pa,pwp)
     &   - 2.*cos2w*dotrc(pz,a) * ( dotrc(pa,wm)*dotcc(wp,z) + dotrc(pa,wp)*dotcc(wm,z) )
     &   + ( cos2w - sin2w ) * ( dotrc(pa,wp)*dotrc(pz,wm)*dotcc(a,z) + dotrc(pa,wm)*dotrc(pz,wp)*dotcc(a,z)
     &                         - dotrc(pa,z)*dotrc(pz,wm)*dotcc(a,wp) - dotrc(pa,z)*dotrc(pz,wp)*dotcc(a,wm) )
     &   + 2.*cos2w*dotrr(pa,pz) * ( dotcc(a,wp)*dotcc(wm,z) + dotcc(a,wm)*dotcc(wp,z) )
     &   + dotrc(pa,wp)*dotrc(pz,wm)*dotcc(a,z) + dotrc(pa,wm)*dotrc(pz,wp)*dotcc(a,z) 
     &   - dotrc(pa,z)*dotrc(pz,wm)*dotcc(a,wp) - dotrc(pa,z)*dotrc(pz,wp)*dotcc(a,wm) )
      endif

      if (fwww .ne. 0) then
        wwza = wwza - fwww * 3./2. * gw**4 * sqrt(sin2w*cos2w) * 
     &      ( (dotrc(pz,a)*(dotrc(pwm,z)+dotrc(pwp,z)) + dotrc(pa,z)*(dotrc(pwm,a)+dotrc(pwp,a))) * dotcc(wp,wm)
     &      - dotcc(a,z) * dotcc(wp,wm) * ( dotrr(pwm,pa) + dotrr(pwm,pz) + dotrr(pwp,pa) + dotrr(pwp,pz) )
     &      - dotcc(z,wm) * ( dotrc(pz,a)*dotrc(pwm,wp) + dotrc(pwm,a)*(dotrc(pa,wp)-dotrc(pz,wp))
     &                      + dotrc(pwp,a)*dotrc(pz,wp) )
     &      - dotcc(a,wm) * ( dotrc(pa,z)*dotrc(pwm,wp) + dotrc(pwp,z)*dotrc(pa,wp)
     &                      + dotrc(pwm,z)*(dotrc(pz,wp)-dotrc(pa,wp)) )
     &      - dotcc(z,wp) * ( dotrc(pz,a)*dotrc(pwp,wm) + dotrc(pwp,a)*(dotrc(pa,wm)-dotrc(pz,wm))
     &                      + dotrc(pwm,a)*dotrc(pz,wm) )
     &      - dotcc(a,wp) * ( dotrc(pa,z)*dotrc(pwp,wm) + dotrc(pwm,z)*dotrc(pa,wm)
     &                      + dotrc(pwp,z)*(dotrc(pz,wm)-dotrc(pa,wm)) )
     &      + dotcc(a,z) * ( (dotrc(pa,wp)+dotrc(pz,wp))*dotrc(pwp,wm) + dotrc(pwm,wp)*(dotrc(pa,wm)+dotrc(pz,wm)) )
     &      + dotcc(a,wm) * dotcc(z,wp) * ( dotrr(pwm,pz) + dotrr(pwp,pa) )
     &      + dotcc(a,wp) * dotcc(z,wm) * ( dotrr(pwm,pa) + dotrr(pwp,pz) ) )
      endif

      if (fw .ne. 0) then
        wwza = wwza + fw * gw**2 * wmass**2 * sqrt(sin2w/cos2w) / 2. * ( - 2.*dotcc(a,z)*dotcc(wm,wp)
     &      + dotcc(a,wp)*dotcc(z,wm) + dotcc(a,wm)*dotcc(z,wp) )
      endif

      if (fdw .ne. 0) then
        wwza = wwza - fdw * gw**4 * 2. * sinw * cosw * (
     &      + dotcc(a,z) * ( - dotrc(pwm,wp) * ( dotrc(pa,wm) + dotrc(pz,wm) - 2.*dotrc(pwp,wm) )
     &                         + dotrc(pa,wp) * ( 2.*dotrc(pa,wm) - 2.*dotrc(pz,wm) - dotrc(pwp,wm) )
     &                         - dotrc(pz,wp) * ( 2.*dotrc(pa,wm) - 2.*dotrc(pz,wm) + dotrc(pwp,wm) ) )
     &      + dotcc(wp,z) * ( + dotrc(pa,wm) * ( dotrc(pz,a) - dotrc(pwm,a) )
     &                         + dotrc(pwp,wm) * ( 2.*dotrc(pz,a) + dotrc(pwm,a) - dotrc(pwp,a) )
     &                         - dotrc(pz,wm) * dotrc(pz,a) )
     &      + dotcc(wm,z) * ( + dotrc(pa,wp) * ( dotrc(pz,a) - dotrc(pwp,a) )
     &                         + dotrc(pwm,wp) * ( 2.*dotrc(pz,a) - dotrc(pwm,a) + dotrc(pwp,a) )
     &                         - dotrc(pz,wp) * dotrc(pz,a) )
     &      - dotcc(wm,a) * ( + dotrc(pz,wp) * ( dotrc(pwp,z) - dotrc(pa,z) )
     &                         + dotrc(pwm,wp) * ( - 2.*dotrc(pa,z) + dotrc(pwm,z) - dotrc(pwp,z) )
     &                         + dotrc(pa,wp) * dotrc(pa,z) )
     &      - dotcc(wp,a) * ( + dotrc(pz,wm) * ( dotrc(pwm,z) - dotrc(pa,z) )
     &                         + dotrc(pwp,wm) * ( - 2.*dotrc(pa,z) - dotrc(pwm,z) + dotrc(pwp,z) )
     &                         + dotrc(pa,wm) * dotrc(pa,z) )
     &      + dotcc(wp,wm) * ( + dotrc(pz,a) * ( 2.*dotrc(pa,z) - dotrc(pwm,z) - dotrc(pwp,z) )
     &                         - dotrc(pwp,a) * ( dotrc(pa,z) + 2.*dotrc(pwm,z) - 2.*dotrc(pwp,z) )
     &                         - dotrc(pwm,a) * ( dotrc(pa,z) - 2.*dotrc(pwm,z) + 2.*dotrc(pwp,z) ) )
     &      - dotcc(wp,wm)*dotcc(a,z) * ( 4.*dotrr(pa,pz) + dotrr(pa,pwm) + dotrr(pa,pwp) 
     &                                    + dotrr(pz,pwm) + dotrr(pz,pwp) + 4.*dotrr(pwm,pwp) )
     &      + dotcc(wp,a)*dotcc(wm,z) * ( dotrr(pa,pz) + 2.*dotrr(pa,pwp) + 2.*dotrr(pz,pwm) + dotrr(pwm,pwp) )
     &      + dotcc(wp,z)*dotcc(wm,a) * ( dotrr(pa,pz) + 2.*dotrr(pa,pwm) + 2.*dotrr(pz,pwp) + dotrr(pwm,pwp) ) )
      endif

      if (fwwwtilde .ne. 0) then
        wwza = wwza - fwwwtilde / 2. * gw**4 * sinw * cosw *
     &     ( dotcc(wp,wm) * ( epsccrr(a,z,pa,pwm) + epsccrr(a,z,pa,pwp) + epsccrr(a,z,pwm,pz) + epsccrr(a,z,pwp,pz) )
     &     + dotcc(wp,a) * ( epsccrr(wm,z,pa,pwm) + epsccrr(wm,z,pa,pwp) - epsccrr(wm,z,pwp,pz) )
     &     + dotcc(wp,z) * ( epsccrr(wm,a,pa,pwp) - epsccrr(wm,a,pwm,pz) - epsccrr(wm,a,pwp,pz) )
     &     + dotcc(wm,a) * ( epsccrr(wp,z,pa,pwm) + epsccrr(wp,z,pa,pwp) - epsccrr(wp,z,pwm,pz) )
     &     + dotcc(wm,z) * ( epsccrr(wp,a,pa,pwm) - epsccrr(wp,a,pwm,pz) - epsccrr(wp,a,pwp,pz) )
     &     + dotcc(a,z) * ( epsccrr(wp,wm,pa,pwm) - epsccrr(wp,wm,pa,pwp) - epsccrr(wp,wm,pwm,pz) + epsccrr(wp,wm,pwp,pz) )
     &     + ( dotrc(pa,wp) - dotrc(pz,wp) ) * ( epscccr(wm,a,z,pwm) + epscccr(wm,a,z,pwp) )
     &     + ( dotrc(pa,wm) - dotrc(pz,wm) ) * ( epscccr(wp,a,z,pwm) + epscccr(wp,a,z,pwp) )
     &     - ( dotrc(pwm,a) - dotrc(pwp,a) ) * ( epscccr(wp,wm,z,pa) + epscccr(wp,wm,z,pz) )
     &     - ( dotrc(pwm,z) - dotrc(pwp,z) ) * ( epscccr(wp,wm,a,pa) + epscccr(wp,wm,a,pz) )
     &     + epscccc(wp,wm,a,z) * ( - dotrr(pa,pwm) + dotrr(pa,pwp) + dotrr(pwm,pz) - dotrr(pwp,pz) )
     &     + dotrc(pa,z) * ( epscccr(wp,wm,a,pwp) - epscccr(wp,wm,a,pwm) )
     &     + dotrc(pwm,wp) * ( epscccr(wm,a,z,pa) - epscccr(wm,a,z,pz) )
     &     + dotrc(pwp,wm) * ( epscccr(wp,a,z,pa) - epscccr(wp,a,z,pz) )
     &     + dotrc(pz,a) * ( epscccr(wp,wm,z,pwp) - epscccr(wp,wm,z,pwm) ) )
      endif

      if (fdwtilde .ne. 0) then
        wwza = wwza - fdwtilde * gw**4 * 2. * sinw * cosw * (
     &      + 2.*dotcc(wp,wm)*epsccrr(a,z,pa,pz) + dotcc(wp,z)*epsccrr(wm,a,pa,pwm)
     &      + dotcc(wm,z)*epsccrr(wp,a,pa,pwp) + dotcc(wp,a)*epsccrr(wm,z,pz,pwm)
     &      + dotcc(wm,a)*epsccrr(wp,z,pz,pwp) - 2.*dotcc(a,z)*epsccrr(wp,wm,pwm,pwp)
     &      + epscccc(wp,wm,a,z) * ( dotrr(pa,pwm) - dotrr(pa,pwp) - dotrr(pz,pwm) +dotrr(pz,pwp) )
     &      + ( dotrc(pa,wp) - dotrc(pz,wp) - dotrc(pwm,wp) ) * epscccr(wm,a,z,pa)
     &      + ( dotrc(pa,wp) - dotrc(pz,wp) + dotrc(pwm,wp) ) * epscccr(wm,a,z,pz)
     &      + ( dotrc(pa,wm) - dotrc(pz,wm) - dotrc(pwp,wm) ) * epscccr(wp,a,z,pa)
     &      + ( dotrc(pa,wm) - dotrc(pz,wm) + dotrc(pwp,wm) ) * epscccr(wp,a,z,pz)
     &      - ( dotrc(pa,z) + dotrc(pwm,z) - dotrc(pwp,z) ) * epscccr(wp,wm,a,pwp)
     &      + ( dotrc(pa,z) - dotrc(pwm,z) + dotrc(pwp,z) ) * epscccr(wp,wm,a,pwm)
     &      - ( dotrc(pz,a) + dotrc(pwm,a) - dotrc(pwp,a) ) * epscccr(wp,wm,z,pwp)
     &      + ( dotrc(pz,a) - dotrc(pwm,a) + dotrc(pwp,a) ) * epscccr(wp,wm,z,pwm) )
      endif

      return
      end


*********************************************************************************************************
CCCCCCCCCCCCCCC
C WWZZ-vertex C
CCCCCCCCCCCCCCC

      subroutine wwzz_anomal4(wm,z1,wp,z2, wwzz)
c
c This subroutine computes the wwzz-vertex with anomalous couplings
c where all the bosons are defined by the flowing-OUT quantum number.
c     
      implicit none
      double complex wp(6),wm(6),z1(6),z2(6),wwzz 
      double complex dotcc, dotrc, epsccrr, epscccr, epscccc
      double precision pwp(0:3),pwm(0:3),pz1(0:3),pz2(0:3), dotrr
      double precision cos2w, g1, cosw, sinw

      external dotcc, dotrc, dotrr, epsccrr, epscccr, epscccc

#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"   

      cos2w=1.0d0 - sin2w
      cosw = sqrt(cos2w)
      sinw = sqrt(sin2w)
      g1 = ee / cosw

      pwm(0) = dble( wm(5))
      pwm(1) = dble( wm(6))
      pwm(2) = dimag(wm(6))
      pwm(3) = dimag(wm(5))
      pwp(0) = dble( wp(5))
      pwp(1) = dble( wp(6))
      pwp(2) = dimag(wp(6))
      pwp(3) = dimag(wp(5))
      pz1(0) = dble( z1(5))
      pz1(1) = dble( z1(6))
      pz1(2) = dimag(z1(6))
      pz1(3) = dimag(z1(5))
      pz2(0) = dble( z2(5))
      pz2(1) = dble( z2(6))
      pz2(2) = dimag(z2(6))
      pz2(3) = dimag(z2(5))

      wwzz = gwwz**2 * (dotcc(wm,z1)*dotcc(wp,z2) + dotcc(wm,z2)*dotcc(wp,z1) - 2.0d0*dotcc(wm,wp)*dotcc(z1,z2))


      if (kmatrix) then
           ! This is the fs0/fs1-amplitude separated in s/t/u-channels and unitarized
           ! following the kmatrix prescription (see subroutine anomal_formfactor)
        if (fs2 .ne. fs0) then
          write(*,*) "Sorry, K-matrix is only available for fs0 == fs2"
          stop
        endif
        if (fs0 .ne. 0 .or. fs1 .ne. 0 .or. fs2 .ne. 0) then
         wwzz = wwzz 
     &         + 2. * wmass**4 / cos2w * fswmwp * dotcc(wm,wp)*dotcc(z1,z2)
     &         +      wmass**4 / cos2w * fswmz1 * dotcc(wm,z1)*dotcc(wp,z2)
     &         +      wmass**4 / cos2w * fswmz2 * dotcc(wm,z2)*dotcc(wp,z1)
        
        endif
      else ! This is the amplitude without k-matrix unitarisation separated in fs1
           ! and fs0 part
      if (fs0 .ne. 0) then
        wwzz = wwzz + fs0 * wmass**4 / cos2w * (dotcc(wm,z1)*dotcc(wp,z2) + dotcc(wm,z2)*dotcc(wp,z1))
      endif

      if (fs1 .ne. 0) then
        wwzz = wwzz + fs1 * wmass**4 / cos2w * 2.0d0 * dotcc(wm,wp)*dotcc(z1,z2)
      endif

      if (fs2 .ne. 0) then
        wwzz = wwzz + fs2 * wmass**4 / cos2w * (dotcc(wm,z1)*dotcc(wp,z2) + dotcc(wm,z2)*dotcc(wp,z1))
      endif
      endif

      if (ft0 .ne. 0) then
        wwzz = wwzz + ft0 * 8. * gw**4 * cos2w * ( dotrc(pwp,wm)*dotrc(pwm,wp) - dotcc(wp,wm)*dotrr(pwp,pwm) )
     &   * ( dotrc(pz1,z2)*dotrc(pz2,z1) - dotcc(z1,z2)*dotrr(pz1,pz2) )
      endif

      if (ft1 .ne. 0) then
        wwzz = wwzz + ft1 * 4. * gw**4 * cos2w * 
     &    ( ( dotrc(pwm,z2)*dotrc(pz2,wm) - dotcc(wm,z2)*dotrr(pwm,pz2) ) 
     &      * ( dotrc(pwp,z1)*dotrc(pz1,wp) - dotcc(wp,z1)*dotrr(pwp,pz1) )
     &    + ( dotrc(pwm,z1)*dotrc(pz1,wm) - dotcc(wm,z1)*dotrr(pwm,pz1) ) 
     &      * ( dotrc(pwp,z2)*dotrc(pz2,wp) - dotcc(wp,z2)*dotrr(pwp,pz2) ) )
      endif

      if (ft2 .ne. 0) then
        wwzz = wwzz + ft2 * gw**4 * cos2w * 
     &    ( dotcc(wm,wp)*dotcc(z1,z2) * ( dotrr(pwm,pz2)*dotrr(pwp,pz1) + dotrr(pwm,pz1)*dotrr(pwp,pz2) )
     &    + dotcc(z1,z2) * ( - dotrc(pwp,wm) * ( dotrc(pz2,wp)*dotrr(pwm,pz1) + dotrc(pz1,wp)*dotrr(pwm,pz2) )
     &                       + dotrc(pz2,wm) * ( dotrc(pz1,wp)*dotrr(pwm,pwp) - dotrc(pwm,wp)*dotrr(pwp,pz1) )
     &                       + dotrc(pz1,wm) * ( dotrc(pz2,wp)*dotrr(pwm,pwp) - dotrc(pwm,wp)*dotrr(pwp,pz2) ) )
     &    + dotcc(wm,z2)*dotcc(wp,z1)*dotrr(pwm,pwp)*dotrr(pz1,pz2)
     &    + dotcc(wm,z1)*dotcc(wp,z2)*dotrr(pwm,pwp)*dotrr(pz1,pz2)
     &    - dotcc(wm,wp) * ( + dotrc(pz2,z1) * ( dotrc(pwp,z2)*dotrr(pwm,pz1) + dotrc(pwm,z2)*dotrr(pwp,pz1) )
     &                       + dotrc(pz1,z2) * ( dotrc(pwp,z1)*dotrr(pwm,pz2) + dotrc(pwm,z1)*dotrr(pwp,pz2) )
     &                       - dotrr(pz1,pz2) * ( dotrc(pwm,z2)*dotrc(pwp,z1) + dotrc(pwm,z1)*dotrc(pwp,z2) ) )
     &    - dotcc(wp,z2) * ( dotrc(pz1,wm)*dotrc(pz2,z1)*dotrr(pwm,pwp)
     &                     + dotrc(pwp,wm) * ( dotrc(pwm,z1)*dotrr(pz1,pz2) - dotrc(pz2,z1)*dotrr(pwm,pz1) ) )
     &    - dotcc(wm,z2) * ( dotrc(pz1,wp)*dotrc(pz2,z1)*dotrr(pwm,pwp)
     &                     + dotrc(pwm,wp) * ( dotrc(pwp,z1)*dotrr(pz1,pz2) - dotrc(pz2,z1)*dotrr(pwp,pz1) ) )
     &    - dotcc(wp,z1) * ( dotrc(pz2,wm)*dotrc(pz1,z2)*dotrr(pwm,pwp)
     &                     + dotrc(pwp,wm) * ( dotrc(pwm,z2)*dotrr(pz1,pz2) - dotrc(pz1,z2)*dotrr(pwm,pz2) ) )
     &    - dotcc(wm,z1) * ( dotrc(pz2,wp)*dotrc(pz1,z2)*dotrr(pwm,pwp)
     &                     + dotrc(pwm,wp) * ( dotrc(pwp,z2)*dotrr(pz1,pz2) - dotrc(pz1,z2)*dotrr(pwp,pz2) ) )
     &    + dotrc(pwm,z2)*dotrc(pwp,wm)*dotrc(pz1,wp)*dotrc(pz2,z1)
     &    + dotrc(pwm,wp)*dotrc(pwp,z2)*dotrc(pz1,wm)*dotrc(pz2,z1)
     &    + dotrc(pwm,z1)*dotrc(pwp,wm)*dotrc(pz1,z2)*dotrc(pz2,wp)
     &    + dotrc(pwm,wp)*dotrc(pwp,z1)*dotrc(pz1,z2)*dotrc(pz2,wm) )
      endif

      if (ft5 .ne. 0) then
        wwzz = wwzz + ft5 * 2. * gw**4 * sin2w**2 / cos2w * ( dotrc(pwp,wm)*dotrc(pwm,wp) - dotcc(wp,wm)*dotrr(pwp,pwm) )
     &   * ( dotrc(pz1,z2)*dotrc(pz2,z1) - dotcc(z1,z2)*dotrr(pz1,pz2) )
      endif

      if (ft6 .ne. 0) then
        wwzz = wwzz + ft6 * gw**4 * sin2w**2 / cos2w *
     &    ( ( dotrc(pwp,z1)*dotrc(pz1,wp) - dotcc(z1,wp)*dotrr(pz1,pwp) ) 
     &      * ( dotrc(pwm,z2)*dotrc(pz2,wm) - dotcc(z2,wm)*dotrr(pz2,pwm) )
     &    + ( dotrc(pwm,z1)*dotrc(pz1,wm) - dotcc(z1,wm)*dotrr(pz1,pwm) ) 
     &      * ( dotrc(pwp,z2)*dotrc(pz2,wp) - dotcc(z2,wp)*dotrr(pz2,pwp) ) )
      endif

      if (ft7 .ne. 0) then
        wwzz = wwzz + ft7 / 4. * gw**4 * sin2w**2 / cos2w *
     &    ( dotcc(wm,wp)*dotcc(z1,z2) * ( dotrr(pwm,pz2)*dotrr(pwp,pz1) + dotrr(pwm,pz1)*dotrr(pwp,pz2) )
     &    + dotcc(z1,z2) * ( - dotrc(pwp,wm) * ( dotrc(pz2,wp)*dotrr(pwm,pz1) + dotrc(pz1,wp)*dotrr(pwm,pz2) )
     &                       + dotrc(pz2,wm) * ( dotrc(pz1,wp)*dotrr(pwm,pwp) - dotrc(pwm,wp)*dotrr(pwp,pz1) )
     &                       + dotrc(pz1,wm) * ( dotrc(pz2,wp)*dotrr(pwm,pwp) - dotrc(pwm,wp)*dotrr(pwp,pz2) ) )
     &    + dotcc(wm,z2)*dotcc(wp,z1)*dotrr(pwm,pwp)*dotrr(pz1,pz2)
     &    + dotcc(wm,z1)*dotcc(wp,z2)*dotrr(pwm,pwp)*dotrr(pz1,pz2)
     &    - dotcc(wm,wp) * ( + dotrc(pz2,z1) * ( dotrc(pwp,z2)*dotrr(pwm,pz1) + dotrc(pwm,z2)*dotrr(pwp,pz1) )
     &                       + dotrc(pz1,z2) * ( dotrc(pwp,z1)*dotrr(pwm,pz2) + dotrc(pwm,z1)*dotrr(pwp,pz2) )
     &                       - dotrr(pz1,pz2) * ( dotrc(pwm,z2)*dotrc(pwp,z1) + dotrc(pwm,z1)*dotrc(pwp,z2) ) )
     &    - dotcc(wp,z2) * ( dotrc(pz1,wm)*dotrc(pz2,z1)*dotrr(pwm,pwp)
     &                     + dotrc(pwp,wm) * ( dotrc(pwm,z1)*dotrr(pz1,pz2) - dotrc(pz2,z1)*dotrr(pwm,pz1) ) )
     &    - dotcc(wm,z2) * ( dotrc(pz1,wp)*dotrc(pz2,z1)*dotrr(pwm,pwp)
     &                     + dotrc(pwm,wp) * ( dotrc(pwp,z1)*dotrr(pz1,pz2) - dotrc(pz2,z1)*dotrr(pwp,pz1) ) )
     &    - dotcc(wp,z1) * ( dotrc(pz1,z2)*dotrc(pz2,wm)*dotrr(pwm,pwp)
     &                     + dotrc(pwp,wm) * ( dotrc(pwm,z2)*dotrr(pz1,pz2) - dotrc(pz1,z2)*dotrr(pwm,pz2) ) )
     &    - dotcc(wm,z1) * ( dotrc(pz1,z2)*dotrc(pz2,wp)*dotrr(pwm,pwp)
     &                     + dotrc(pwm,wp) * ( dotrc(pwp,z2)*dotrr(pz1,pz2) - dotrc(pz1,z2)*dotrr(pwp,pz2) ) )
     &    + dotrc(pwm,z2)*dotrc(pwp,wm)*dotrc(pz1,wp)*dotrc(pz2,z1)
     &    + dotrc(pwm,wp)*dotrc(pwp,z2)*dotrc(pz1,wm)*dotrc(pz2,z1)
     &    + dotrc(pwm,z1)*dotrc(pwp,wm)*dotrc(pz1,z2)*dotrc(pz2,wp)
     &    + dotrc(pwm,wp)*dotrc(pwp,z1)*dotrc(pz1,z2)*dotrc(pz2,wm) )
      endif

      if (fm0 .ne. 0) then
        wwzz = wwzz - fm0 * 2. * gw**2 * wmass**2 *
     &      ( 1. / cos2w * dotcc(z1,z2) * ( dotrc(pwm,wp)*dotrc(pwp,wm) - dotcc(wm,wp)*dotrr(pwm,pwp) )
     &      + cos2w * dotcc(wm,wp) * ( dotrc(pz1,z2)*dotrc(pz2,z1) - dotcc(z1,z2)*dotrr(pz1,pz2) ) )
      endif

      if (fm1 .ne. 0) then
        wwzz = wwzz + fm1 / 16. * gw**2 * wmass**2 / cos2w * 
     &        ( 8. * ( dotrc(pwp,wm) * ( dotrc(pwm,z1)*dotcc(wp,z2) + dotrc(pwm,z2)*dotcc(wp,z1) )
     &               + dotrc(pwm,wp) * ( dotrc(pwp,z1)*dotcc(wm,z2) + dotrc(pwp,z2)*dotcc(wm,z1) )
     &               - dotcc(wm,wp) * ( dotrc(pwm,z2)*dotrc(pwp,z1) + dotrc(pwm,z1)*dotrc(pwp,z2) ) )
     &        - 8. * dotcc(wm,z2)*dotcc(wp,z1)*dotrr(pwm,pwp) - 8. * dotcc(wm,z1)*dotcc(wp,z2)*dotrr(pwm,pwp)
     &        + 8. * cos2w**2 * dotrc(pz1,wm) * ( dotrc(pz2,z1)*dotcc(wp,z2) - dotrc(pz2,wp)*dotcc(z1,z2) )
     &        + 8. * cos2w**2 * dotrc(pz1,z2)*dotrc(pz2,wp)*dotcc(wm,z1)
     &        + 8. * cos2w**2 * dotrc(pz2,wm) * ( dotrc(pz1,z2)*dotcc(wp,z1) - dotrc(pz1,wp)*dotcc(z1,z2) )
     &        + ( sin2w**2 - 4.*sin2w + cos2w**2 + 4.*cos2w - 6.*sin2w*cos2w ) * dotrc(pz1,wp)*dotrc(pz2,z1)*dotcc(wm,z2)
     &        - 8. * cos2w**2 * dotrr(pz1,pz2) * ( dotcc(wm,z2)*dotcc(wp,z1) + dotcc(wm,z1)*dotcc(wp,z2) )
     &        + 3. * dotrc(pz1,wp)*dotrc(pz2,z1)*dotcc(wm,z2) )
      endif

      if (fm2 .ne. 0) then
        wwzz = wwzz - fm2 * gw**2 * wmass**2 * sin2w**2 / cos2w * dotcc(wp,wm) * 
     &   ( dotrc(pz1,z2) * dotrc(pz2,z1) - dotcc(z1,z2) * dotrr(pz1,pz2) )
      endif

      if (fm3 .ne. 0) then
        wwzz = wwzz + fm3 / 4. * gw**2 * wmass**2 * sin2w**2 / cos2w *
     &    ( dotrc(pz1,wp) * dotrc(pz2,z1) * dotcc(wm,z2)
     &    + dotrc(pz1,wm) * ( dotrc(pz2,z1)*dotcc(wp,z2) - dotrc(pz2,wp)*dotcc(z1,z2) )
     &    + dotrc(pz1,z2) * dotrc(pz2,wp) * dotcc(wm,z1)
     &    + dotrc(pz2,wm) * ( dotrc(pz1,z2)*dotcc(wp,z1) - dotrc(pz1,wp)*dotcc(z1,z2) )
     &    - dotrr(pz1,pz2) * ( dotcc(wm,z2)*dotcc(wp,z1) + dotcc(wm,z1)*dotcc(wp,z2) ) )
      endif

      if (fm4 .ne. 0) then
        wwzz = wwzz - fm4 / 2. * gw**2 * wmass**2 * sin2w / cos2w *
     &   ( dotcc(wm,z2) * ( dotrc(pwp,z1)*dotrc(pz1,wp) - dotcc(wp,z1) * ( dotrr(pwm,pz2) + dotrr(pwp,pz1) ) )
     &   + dotcc(wm,z1) * ( dotrc(pwp,z2)*dotrc(pz2,wp) - dotcc(wp,z2) * ( dotrr(pwm,pz1) + dotrr(pwp,pz2) ) )
     &   + dotrc(pwm,z1)*dotrc(pz1,wm)*dotcc(wp,z2) + dotrc(pwm,z2)*dotrc(pz2,wm)*dotcc(wp,z1)
     &   + 2. * cos2w * dotcc(wm,wp) * ( dotcc(z1,z2)*dotrr(pz1,pz2) - dotrc(pz1,z2)*dotrc(pz2,z1) ) )
      endif

      if (fm5 .ne. 0) then
        wwzz = wwzz - fm5 / 8. * gw**2 * wmass**2 * sin2w / cos2w *
     &    ( 2.*dotrc(pz1,wm) * ( dotcc(z1,z2) * ( dotrc(pwm,wp) + cos2w*dotrc(pz2,wp) )
     &                         + dotrc(pwp,z1)*dotcc(wp,z2) - dotrc(pwp,z2)*dotcc(wp,z1)
     &                         - cos2w*dotrc(pz2,z1)*dotcc(wp,z2) )
     &    + 2.*dotrc(pz2,wm) * ( dotcc(z1,z2) * ( dotrc(pwm,wp) + cos2w*dotrc(pz1,wp) )
     &                         - dotrc(pwp,z1)*dotcc(wp,z2) + dotrc(pwp,z2)*dotcc(wp,z1)
     &                         - cos2w*dotrc(pz1,z2)*dotcc(wp,z1) )
     &    - 2.*dotrc(pwm,wp)*dotrc(pz1,z2)*dotcc(wm,z1) + 2.*dotrc(pwm,z1)*dotrc(pz1,z2)*dotcc(wm,wp)
     &    - 2.*dotcc(wm,wp)*dotcc(z1,z2)*dotrr(pwm,pz1) - 2.*dotrc(pwm,wp)*dotrc(pz2,z1)*dotcc(wm,z2)
     &    + 2.*dotrc(pwm,z2)*dotrc(pz2,z1)*dotcc(wm,wp) - 2.*dotcc(wm,wp)*dotcc(z1,z2)*dotrr(pwm,pz2)
     &    + 2.*dotrc(pwp,z2)*dotrc(pz1,wp)*dotcc(wm,z1) - 2.*dotcc(wm,z1)*dotcc(wp,z2)*dotrr(pwp,pz1)
     &    + 2.*dotrc(pwp,z1)*dotrc(pz2,wp)*dotcc(wm,z2) - 2.*dotcc(wm,z2)*dotcc(wp,z1)*dotrr(pwp,pz2)
     &    - ( cos2w - sin2w ) * dotrc(pz1,wp)*dotrc(pz2,z1)*dotcc(wm,z2)
     &    - ( cos2w - sin2w ) * dotrc(pz1,z2)*dotrc(pz2,wp)*dotcc(wm,z1) 
     &    + 2.*cos2w*dotrr(pz1,pz2) * ( dotcc(wm,z2)*dotcc(wp,z1) + dotcc(wm,z1)*dotcc(wp,z2) )
     &    - dotrc(pz1,wp)*dotrc(pz2,z1)*dotcc(wm,z2) - dotrc(pz1,z2)*dotrc(pz2,wp)*dotcc(wm,z1) )
      endif

      if (fm6 .ne. 0) then
        wwzz = wwzz - fm6 * gw**2 * wmass**2 *
     &   ( 1. / cos2w * dotcc(z1,z2) * ( dotrc(pwm,wp)*dotrc(pwp,wm) - dotcc(wm,wp)*dotrr(pwm,pwp) )
     &   + cos2w * dotcc(wp,wm) * ( dotrc(pz1,z2)*dotrc(pz2,z1) - dotcc(z1,z2)*dotrr(pz1,pz2) ) )
      endif

      if (fm7 .ne. 0) then
        wwzz = wwzz + fm7 / 8. * gw**2 * wmass**2 * (
     &   - 2. / cos2w * ( dotrc(pwp,wm) * ( dotcc(wp,z2)*dotrc(pwm,z1) + dotcc(wp,z1)*dotrc(pwm,z2) )
     &                  + dotrc(pwm,wp) * ( dotcc(wm,z2)*dotrc(pwp,z1) + dotcc(wm,z1)*dotrc(pwp,z2) )
     &                  - dotcc(wm,wp) * ( dotrc(pwp,z1)*dotrc(pwm,z2) + dotrc(pwm,z1)*dotrc(pwp,z2) )
     &                  - dotrr(pwm,pwp) * ( dotcc(wm,z2)*dotcc(wp,z1) + dotcc(wm,z1)*dotcc(wp,z2) ) )
     &  + 2.*dotrc(pwp,wm)*dotrc(pz1,wp)*dotcc(z1,z2) + 2.*dotrc(pwp,wm)*dotrc(pz2,wp)*dotcc(z1,z2)
     &  - 2.*dotcc(wm,z2)*dotrc(pz1,wp)*dotrc(pwm,z1) + 2.*dotcc(wm,z2)*dotrc(pz2,wp)*dotrc(pwm,z1)
     &  - 2.*dotcc(wm,z2)*dotrc(pz2,wp)*dotrc(pwp,z1) - 2.*dotcc(wp,z2)*dotrc(pz2,z1)*dotrc(pwp,wm)
     &  - 2.*dotcc(wm,z2)*dotrc(pwm,wp)*dotrc(pz2,z1) - dotcc(wm,z2)*dotrc(pz1,wp)*dotrc(pz2,z1)
     &  + 2.*dotcc(wm,z1)*dotrc(pz1,wp)*dotrc(pwm,z2) - 2.*dotcc(wm,z1)*dotrc(pz2,wp)*dotrc(pwm,z2)
     &  + 2.*dotcc(wm,wp)*dotrc(pz2,z1)*dotrc(pwm,z2) - 2.*dotcc(wm,z1)*dotrc(pz1,wp)*dotrc(pwp,z2)
     &  + 2.*dotcc(wm,wp)*dotrc(pz2,z1)*dotrc(pwp,z2)
     &  + 2.*dotrc(pz1,wm) * ( cos2w*dotrc(pz2,wp)*dotcc(z1,z2) + dotrc(pwm,wp)*dotcc(z1,z2)
     &                       - dotcc(wp,z2) * ( cos2w*dotrc(pz2,z1) + dotrc(pwp,z1) )
     &                       + dotcc(wp,z1) * ( dotrc(pwp,z2) - dotrc(pwm,z2) ) )
     &  - 2.*dotrc(pwp,wm)*dotrc(pz1,z2)*dotcc(wp,z1) - 2.*dotcc(wm,z1)*dotrc(pwm,wp)*dotrc(pz1,z2)
     &  - dotcc(wm,z1)*dotrc(pz2,wp)*dotrc(pz1,z2) + 2.*dotcc(wm,wp)*dotrc(pwm,z1)*dotrc(pz1,z2)
     &  + 2.*dotcc(wm,wp)*dotrc(pwp,z1)*dotrc(pz1,z2)
     &  + 2.*dotrc(pz2,wm) * ( dotcc(z1,z2) * ( cos2w*dotrc(pz1,wp) + dotrc(pwm,wp) )
     &                       - dotcc(wp,z2) * ( dotrc(pwm,z1) - dotrc(pwp,z1) )
     &                       - dotcc(wp,z1) * ( cos2w*dotrc(pz1,z2) + dotrc(pwp,z2) ) )
     &  + 2.*dotcc(wm,z2)*dotcc(wp,z1)*dotrr(pwm,pz1) - 2.*dotcc(wm,wp)*dotcc(z1,z2)*dotrr(pwm,pz1)
     &  + 2.*dotcc(wm,z1)*dotcc(wp,z2)*dotrr(pwm,pz2) - 2.*dotcc(wm,wp)*dotcc(z1,z2)*dotrr(pwm,pz2)
     &  + 2.*dotcc(wm,z1)*dotcc(wp,z2)*dotrr(pwp,pz1) - 2.*dotcc(wm,wp)*dotcc(z1,z2)*dotrr(pwp,pz1)
     &  + 2.*dotcc(wm,z2)*dotcc(wp,z1)*dotrr(pwp,pz2) - 2.*dotcc(wm,wp)*dotcc(z1,z2)*dotrr(pwp,pz2)
     &  + dotcc(wm,z2)*dotcc(wp,z1)*dotrr(pz1,pz2) + dotcc(wm,z1)*dotcc(wp,z2)*dotrr(pz1,pz2)
     &  + ( cos2w - sin2w ) * ( dotcc(wm,z2) * ( dotcc(wp,z1)*dotrr(pz1,pz2) - dotrc(pz1,wp)*dotrc(pz2,z1) )
     &                        + dotcc(wm,z1) * ( dotcc(wp,z2)*dotrr(pz1,pz2) - dotrc(pz2,wp)*dotrc(pz1,z2) ) ) )
      endif

      if (fwww .ne. 0) then
        wwzz = wwzz - fwww * 3./2. * gw**4 * cos2w * 
     &      ( ( dotrc(pz2,z1)*(dotrc(pwm,z2)+dotrc(pwp,z2)) + dotrc(pz1,z2)*(dotrc(pwm,z1)+dotrc(pwp,z1)) ) * dotcc(wp,wm)
     &      - dotcc(z1,z2) * dotcc(wm,wp) * ( dotrr(pwm,pz1) + dotrr(pwm,pz2) + dotrr(pwp,pz1) + dotrr(pwp,pz2) )
     &      - dotcc(z2,wm) * ( dotrc(pz2,z1)*dotrc(pwm,wp) + dotrc(pwm,z1)*(dotrc(pz1,wp)-dotrc(pz2,wp))
     &                       + dotrc(pwp,z1)*dotrc(pz2,wp) )
     &      - dotcc(z1,wm) * ( dotrc(pz1,z2)*dotrc(pwm,wp) + dotrc(pwp,z2)*dotrc(pz1,wp)
     &                       + dotrc(pwm,z2)*(dotrc(pz2,wp)-dotrc(pz1,wp)) )
     &      - dotcc(z2,wp) * ( dotrc(pz2,z1)*dotrc(pwp,wm) + dotrc(pwp,z1)*(dotrc(pz1,wm)-dotrc(pz2,wm))
     &                       + dotrc(pwm,z1)*dotrc(pz2,wm) )
     &      - dotcc(z1,wp) * ( dotrc(pz1,z2)*dotrc(pwp,wm) + dotrc(pwm,z2)*dotrc(pz1,wm)
     &                       + dotrc(pwp,z2)*(dotrc(pz2,wm)-dotrc(pz1,wm)) )
     &      + dotcc(z1,z2) * ( (dotrc(pz1,wp)+dotrc(pz2,wp))*dotrc(pwp,wm) + dotrc(pwm,wp)*(dotrc(pz1,wm)+dotrc(pz2,wm)) )
     &      + dotcc(z1,wm) * dotcc(z2,wp) * ( dotrr(pwm,pz2) + dotrr(pwp,pz1) )
     &      + dotcc(z1,wp) * dotcc(z2,wm) * ( dotrr(pwm,pz1) + dotrr(pwp,pz2) ) )
      endif

      if (fw .ne. 0) then
        wwzz = wwzz + fw * gw**2 * wmass**2 * ( - 2.*dotcc(z1,z2)*dotcc(wm,wp) + dotcc(z1,wp)*dotcc(z2,wm)
     &      + dotcc(z1,wm)*dotcc(z2,wp) )
      endif

      if (fdw .ne. 0) then
        wwzz = wwzz - fdw * gw**4 * 2. * cos2w * (
     &      + dotcc(z1,z2) * ( - dotrc(pwm,wp) * ( dotrc(pz1,wm) + dotrc(pz2,wm) - 2.*dotrc(pwp,wm) )
     &                         + dotrc(pz1,wp) * ( 2.*dotrc(pz1,wm) - 2.*dotrc(pz2,wm) - dotrc(pwp,wm) )
     &                         - dotrc(pz2,wp) * ( 2.*dotrc(pz1,wm) - 2.*dotrc(pz2,wm) + dotrc(pwp,wm) ) )
     &      + dotcc(wp,z2) * ( + dotrc(pz1,wm) * ( dotrc(pz2,z1) - dotrc(pwm,z1) )
     &                         + dotrc(pwp,wm) * ( 2.*dotrc(pz2,z1) + dotrc(pwm,z1) - dotrc(pwp,z1) )
     &                         - dotrc(pz2,wm) * dotrc(pz2,z1) )
     &      + dotcc(wm,z2) * ( + dotrc(pz1,wp) * ( dotrc(pz2,z1) - dotrc(pwp,z1) )
     &                         + dotrc(pwm,wp) * ( 2.*dotrc(pz2,z1) - dotrc(pwm,z1) + dotrc(pwp,z1) )
     &                         - dotrc(pz2,wp) * dotrc(pz2,z1) )
     &      - dotcc(wm,z1) * ( + dotrc(pz2,wp) * ( dotrc(pwp,z2) - dotrc(pz1,z2) )
     &                         + dotrc(pwm,wp) * ( - 2.*dotrc(pz1,z2) + dotrc(pwm,z2) - dotrc(pwp,z2) )
     &                         + dotrc(pz1,wp) * dotrc(pz1,z2) )
     &      - dotcc(wp,z1) * ( + dotrc(pz2,wm) * ( dotrc(pwm,z2) - dotrc(pz1,z2) )
     &                         + dotrc(pwp,wm) * ( - 2.*dotrc(pz1,z2) - dotrc(pwm,z2) + dotrc(pwp,z2) )
     &                         + dotrc(pz1,wm) * dotrc(pz1,z2) )
     &      + dotcc(wp,wm) * ( + dotrc(pz2,z1) * ( 2.*dotrc(pz1,z2) - dotrc(pwm,z2) - dotrc(pwp,z2) )
     &                         - dotrc(pwp,z1) * ( dotrc(pz1,z2) + 2.*dotrc(pwm,z2) - 2.*dotrc(pwp,z2) )
     &                         - dotrc(pwm,z1) * ( dotrc(pz1,z2) - 2.*dotrc(pwm,z2) + 2.*dotrc(pwp,z2) ) )
     &      - dotcc(wp,wm)*dotcc(z1,z2) * ( 4.*dotrr(pz1,pz2) + dotrr(pz1,pwm) + dotrr(pz1,pwp) 
     &                                      + dotrr(pz2,pwm) + dotrr(pz2,pwp) + 4.*dotrr(pwm,pwp) )
     &      + dotcc(wp,z1)*dotcc(wm,z2) * ( dotrr(pz1,pz2) + 2.*dotrr(pz1,pwp) + 2.*dotrr(pz2,pwm) + dotrr(pwm,pwp) )
     &      + dotcc(wp,z2)*dotcc(wm,z1) * ( dotrr(pz1,pz2) + 2.*dotrr(pz1,pwm) + 2.*dotrr(pz2,pwp) + dotrr(pwm,pwp) ) )
      endif

      if (fwwwtilde .ne. 0) then
        wwzz = wwzz + fwwwtilde / 2. * gw**4 * cos2w *
     &    ( dotcc(wp,wm) * ( epsccrr(z1,z2,pwm,pz1) - epsccrr(z1,z2,pwm,pz2) + epsccrr(z1,z2,pwp,pz1) - epsccrr(z1,z2,pwp,pz2) )
     &    + dotcc(wp,z1) * ( epsccrr(wm,z2,pwm,pz1) + epsccrr(wm,z2,pwp,pz1) + epsccrr(wm,z2,pwp,pz2) )
     &    + dotcc(wp,z2) * ( epsccrr(wm,z1,pwm,pz2) + epsccrr(wm,z1,pwp,pz1) + epsccrr(wm,z1,pwp,pz2) )
     &    + dotcc(wm,z1) * ( epsccrr(wp,z2,pwm,pz1) + epsccrr(wp,z2,pwm,pz2) + epsccrr(wp,z2,pwp,pz1) )
     &    + dotcc(wm,z2) * ( epsccrr(wp,z1,pwm,pz1) + epsccrr(wp,z1,pwm,pz2) + epsccrr(wp,z1,pwp,pz2) )
     &    + dotcc(z1,z2) * ( epsccrr(wp,wm,pwm,pz1) + epsccrr(wp,wm,pwm,pz2) - epsccrr(wp,wm,pwp,pz1) - epsccrr(wp,wm,pwp,pz2) )
     &    - ( dotrc(pz1,wp) - dotrc(pz2,wp) ) * ( epscccr(wm,z1,z2,pwm) + epscccr(wm,z1,z2,pwp) )
     &    - ( dotrc(pz1,wm) - dotrc(pz2,wm) ) * ( epscccr(wp,z1,z2,pwm) + epscccr(wp,z1,z2,pwp) )
     &    + ( dotrc(pwm,z1) - dotrc(pwp,z1) ) * ( epscccr(wp,wm,z2,pz1) + epscccr(wp,wm,z2,pz2) )
     &    + ( dotrc(pwm,z2) - dotrc(pwp,z2) ) * ( epscccr(wp,wm,z1,pz1) + epscccr(wp,wm,z1,pz2) )
     &    + epscccc(wp,wm,z1,z2) * ( dotrr(pwm,pz1) - dotrr(pwm,pz2) - dotrr(pwp,pz1) + dotrr(pwp,pz2) )
     &    + dotrc(pz1,z2) * ( epscccr(wp,wm,z1,pwm) - epscccr(wp,wm,z1,pwp) )
     &    + dotrc(pz2,z1) * ( epscccr(wp,wm,z2,pwm) - epscccr(wp,wm,z2,pwp) )
     &    + dotrc(pwm,wp) * ( epscccr(wm,z1,z2,pz2) - epscccr(wm,z1,z2,pz1) )
     &    + dotrc(pwp,wm) * ( epscccr(wp,z1,z2,pz2) - epscccr(wp,z1,z2,pz1) ) )
      endif

      if (fdwtilde .ne. 0) then
        wwzz = wwzz - fdwtilde * gw**4 * 2. * cos2w * (
     &      + 2.*dotcc(wp,wm)*epsccrr(z1,z2,pz1,pz2) + dotcc(wp,z2)*epsccrr(wm,z1,pz1,pwm)
     &      + dotcc(wm,z2)*epsccrr(wp,z1,pz1,pwp) + dotcc(wp,z1)*epsccrr(wm,z2,pz2,pwm)
     &      + dotcc(wm,z1)*epsccrr(wp,z2,pz2,pwp) - 2.*dotcc(z1,z2)*epsccrr(wp,wm,pwm,pwp)
     &      + epscccc(wp,wm,z1,z2) * ( dotrr(pwm,pz1) - dotrr(pwm,pz2) - dotrr(pwp,pz1) + dotrr(pwp,pz2) )
     &      + ( dotrc(pz1,wp) - dotrc(pz2,wp) - dotrc(pwm,wp) ) * epscccr(wm,z1,z2,pz1)
     &      + ( dotrc(pz1,wp) - dotrc(pz2,wp) + dotrc(pwm,wp) ) * epscccr(wm,z1,z2,pz2)
     &      + ( dotrc(pz1,wm) - dotrc(pz2,wm) - dotrc(pwp,wm) ) * epscccr(wp,z1,z2,pz1)
     &      + ( dotrc(pz1,wm) - dotrc(pz2,wm) + dotrc(pwp,wm) ) * epscccr(wp,z1,z2,pz2)
     &      - ( dotrc(pz1,z2) + dotrc(pwm,z2) - dotrc(pwp,z2) ) * epscccr(wp,wm,z1,pwp)
     &      + ( dotrc(pz1,z2) - dotrc(pwm,z2) + dotrc(pwp,z2) ) * epscccr(wp,wm,z1,pwm)
     &      - ( dotrc(pz2,z1) + dotrc(pwm,z1) - dotrc(pwp,z1) ) * epscccr(wp,wm,z2,pwp)
     &      + ( dotrc(pz2,z1) - dotrc(pwm,z1) + dotrc(pwp,z1) ) * epscccr(wp,wm,z2,pwm) )
      endif

      return
      end


*********************************************************************************************************
CCCCCCCCCCCCCCC
C WWWW-vertex C
CCCCCCCCCCCCCCC

      subroutine wwww_anomal4(wm1,wp1,wm2,wp2, wwww)
c
c This subroutine computes the wwww-vertex with anomalous couplings
c where all the bosons are defined by the flowing-OUT quantum number.
c     
      implicit none
      double complex wp1(6),wm1(6),wp2(6),wm2(6),wwww 
      double complex dotcc, dotrc, epsccrr, epscccr, epscccc
      double precision pwp1(0:3),pwm1(0:3),pwp2(0:3),pwm2(0:3), dotrr
      double precision cos2w, g1, cosw, sinw

      external dotcc, dotrc, dotrr, epsccrr, epscccr, epscccc

#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"   

      cos2w = 1.0d0 - sin2w
      cosw = sqrt(cos2w)
      sinw = sqrt(sin2w)
      g1 = ee / cosw

      pwm1(0) = dble( wm1(5))
      pwm1(1) = dble( wm1(6))
      pwm1(2) = dimag(wm1(6))
      pwm1(3) = dimag(wm1(5))
      pwp1(0) = dble( wp1(5))
      pwp1(1) = dble( wp1(6))
      pwp1(2) = dimag(wp1(6))
      pwp1(3) = dimag(wp1(5))
      pwm2(0) = dble( wm2(5))
      pwm2(1) = dble( wm2(6))
      pwm2(2) = dimag(wm2(6))
      pwm2(3) = dimag(wm2(5))
      pwp2(0) = dble( wp2(5))
      pwp2(1) = dble( wp2(6))
      pwp2(2) = dimag(wp2(6))
      pwp2(3) = dimag(wp2(5))

      wwww = - ee**2 / sin2w * (dotcc(wm1,wp1)*dotcc(wm2,wp2) + dotcc(wm1,wp2)*dotcc(wm2,wp1) 
     &                        - 2.0d0*dotcc(wm1,wm2)*dotcc(wp1,wp2))


      if (kmatrix) then
           ! This is the fs0/fs1-amplitude separated in s/t/u-channels and unitarized
           ! following the kmatrix prescription (see subroutine anomal_formfactor)
        if (fs2 .ne. fs0) then
          write(*,*) "Sorry, K-matrix is only available for fs0 == fs2"
          stop
        endif
        if (fs0 .ne. 0 .or. fs1 .ne. 0 .or. fs2 .ne. 0) then
         wwww = wwww 
     &         + 4. * wmass**4 * fswp1wp2 * dotcc(wp1,wp2) * dotcc(wm1,wm2) 
     &         + 2. * wmass**4 * fswp1wm2 * dotcc(wp1,wm2) * dotcc(wp2,wm1)
     &         + 2. * wmass**4 * fswp1wm1 * dotcc(wp1,wm1) * dotcc(wp2,wm2)
        endif
      else ! This is the amplitude without k-matrix unitarisation separated in fs1
           ! and fs0 part
       if (fs0 .ne. 0) then
         wwww = wwww + 4. * fs0 * wmass**4 
     &         * dotcc(wp1,wp2) * dotcc(wm1,wm2)
       endif
 
       if (fs1 .ne. 0) then
         wwww = wwww + 2. * fs1 * wmass**4 
     &      * (dotcc(wp1,wm2) * dotcc(wp2,wm1)
     &      +  dotcc(wp1,wm1) * dotcc(wp2,wm2))
       endif

       if (fs2 .ne. 0) then
         wwww = wwww + 2. * fs2 * wmass**4 
     &      * (dotcc(wp1,wm2) * dotcc(wp2,wm1)
     &      +  dotcc(wp1,wm1) * dotcc(wp2,wm2))
       endif
      endif

      if (ft0 .ne. 0) then
        wwww = wwww + 8. * ft0 * gw**4 * (
     &      ( dotrc(pwm1,wp2)*dotrc(pwp2,wm1) - dotcc(wm1,wp2)*dotrr(pwm1,pwp2) )
     &      * ( dotrc(pwm2,wp1)*dotrc(pwp1,wm2) - dotcc(wm2,wp1)*dotrr(pwm2,pwp1) )
     &    + ( dotrc(pwm1,wp1)*dotrc(pwp1,wm1) - dotcc(wm1,wp1)*dotrr(pwm1,pwp1) )
     &      * ( dotrc(pwm2,wp2)*dotrc(pwp2,wm2) - dotcc(wm2,wp2)*dotrr(pwm2,pwp2) ) )
      endif

      if (ft1 .ne. 0) then
        wwww = wwww + 4. * ft1 * gw**4 *
     &    ( ( dotrc(pwm1,wp2)*dotrc(pwp2,wm1) - dotcc(wm1,wp2)*dotrr(pwm1,pwp2) ) 
     &            * ( dotrc(pwm2,wp1)*dotrc(pwp1,wm2) - dotcc(wm2,wp1)*dotrr(pwm2,pwp1) )
     &    + ( dotrc(pwm1,wp1)*dotrc(pwp1,wm1) - dotcc(wm1,wp1)*dotrr(pwm1,pwp1) ) 
     &            * ( dotrc(pwm2,wp2)*dotrc(pwp2,wm2) - dotcc(wm2,wp2)*dotrr(pwm2,pwp2) )
     &    + 2. * ( dotrc(pwm1,wm2)*dotrc(pwm2,wm1) - dotcc(wm1,wm2)*dotrr(pwm1,pwm2) ) 
     &            * ( dotrc(pwp1,wp2)*dotrc(pwp2,wp1) - dotcc(wp1,wp2)*dotrr(pwp1,pwp2) ) )
      endif

      if (ft2 .ne. 0) then
        wwww = wwww + ft2 * gw**4 *
     &    ( 2. * dotrc(pwp1,wm1) * dotrc(pwp2,wm2) * dotrc(pwm2,wp1) * dotrc(pwm1,wp2)
     &    + dotrc(pwm2,wm1) * dotrc(pwp1,wm2) * dotrc(pwp2,wp1) * dotrc(pwm1,wp2)
     &    + 2. * dotrc(pwp2,wm1) * dotrc(pwp1,wm2) * dotrc(pwm1,wp1) * dotrc(pwm2,wp2)
     &    + dotrc(pwp1,wm1) * dotrc(pwm1,wm2) * dotrc(pwp2,wp1) * dotrc(pwm2,wp2)
     &    + dotrc(pwm2,wm1) * dotrc(pwp2,wm2) * dotrc(pwm1,wp1) * dotrc(pwp1,wp2)
     &    + dotrc(pwp2,wm1) * dotrc(pwm1,wm2) * dotrc(pwm2,wp1) * dotrc(pwp1,wp2)
     &    + dotcc(wm1,wm2)*dotcc(wp1,wp2) * ( dotrr(pwm1,pwp2)*dotrr(pwm2,pwp1) + dotrr(pwm1,pwp1)*dotrr(pwm2,pwp2) )
     &    + dotcc(wp1,wp2) * ( - dotrc(pwm2,wm1) * ( dotrc(pwp2,wm2)*dotrr(pwm1,pwp1) + dotrc(pwp1,wm2)*dotrr(pwm1,pwp2) )
     &                         + dotrc(pwp2,wm1) * ( dotrc(pwp1,wm2)*dotrr(pwm1,pwm2) - dotrc(pwm1,wm2)*dotrr(pwm2,pwp1) )
     &                         + dotrc(pwp1,wm1) * ( dotrc(pwp2,wm2)*dotrr(pwm1,pwm2) - dotrc(pwm1,wm2)*dotrr(pwm2,pwp2) ) )
     &    - dotcc(wm1,wm2) * ( + dotrc(pwp2,wp1) * ( dotrc(pwm2,wp2)*dotrr(pwm1,pwp1) + dotrc(pwm1,wp2)*dotrr(pwm2,pwp1) )
     &                         + dotrc(pwp1,wp2) * ( dotrc(pwm2,wp1)*dotrr(pwm1,pwp2) + dotrc(pwm1,wp1)*dotrr(pwm2,pwp2) )
     &                         - dotrr(pwp1,pwp2) * ( dotrc(pwm2,wp1)*dotrc(pwm1,wp2) + dotrc(pwm1,wp1)*dotrc(pwm2,wp2) ) )
     &    + dotcc(wm1,wp1)*dotcc(wm2,wp2) * ( 2.*dotrr(pwm1,pwp2)*dotrr(pwm2,pwp1) + dotrr(pwm1,pwm2)*dotrr(pwp1,pwp2) )
     &    + dotcc(wm1,wp2)*dotcc(wm2,wp1) * ( 2.*dotrr(pwm1,pwp1)*dotrr(pwm2,pwp2) + dotrr(pwm1,pwm2)*dotrr(pwp1,pwp2) )
     &    - dotcc(wm2,wp2) * ( + dotrc(pwp1,wm1) * ( dotrc(pwp2,wp1)*dotrr(pwm1,pwm2) + 2.*dotrc(pwm2,wp1)*dotrr(pwm1,pwp2) )
     &                         + dotrc(pwp2,wm1) * ( 2.*dotrc(pwm1,wp1)*dotrr(pwm2,pwp1) - 2.*dotrc(pwm2,wp1)*dotrr(pwm1,pwp1) )
     &                         + dotrc(pwm2,wm1) * ( dotrc(pwm1,wp1)*dotrr(pwp1,pwp2) - dotrc(pwp2,wp1)*dotrr(pwm1,pwp1) ) )
     &    - dotcc(wm1,wp2) * ( + 2.*dotrc(pwp2,wm2) * ( dotrc(pwm2,wp1)*dotrr(pwm1,pwp1) - dotrc(pwm1,wp1)*dotrr(pwm2,pwp1) )
     &                         + dotrc(pwp1,wm2) * ( dotrc(pwp2,wp1)*dotrr(pwm1,pwm2) + 2.*dotrc(pwm1,wp1)*dotrr(pwm2,pwp2) )
     &                         + dotrc(pwm1,wm2) * ( dotrc(pwm2,wp1)*dotrr(pwp1,pwp2) - dotrc(pwp2,wp1)*dotrr(pwm2,pwp1) ) )
     &    - dotcc(wm2,wp1) * ( + dotrc(pwp2,wm1) * ( dotrc(pwp1,wp2)*dotrr(pwm1,pwm2) + 2.*dotrc(pwm2,wp2)*dotrr(pwm1,pwp1) )
     &                         + dotrc(pwp1,wm1) * ( 2.*dotrc(pwm1,wp2)*dotrr(pwm2,pwp2) - 2.*dotrc(pwm2,wp2)*dotrr(pwm1,pwp2) )
     &                         + dotrc(pwm2,wm1) * ( dotrc(pwm1,wp2)*dotrr(pwp1,pwp2) - dotrc(pwp1,wp2)*dotrr(pwm1,pwp2) ) )
     &    - dotcc(wm1,wp1) * ( + dotrc(pwp2,wm2) * ( dotrc(pwp1,wp2)*dotrr(pwm1,pwm2) + 2.*dotrc(pwm1,wp2)*dotrr(pwm2,pwp1) )
     &                         + dotrc(pwp1,wm2) * ( 2.*dotrc(pwm2,wp2)*dotrr(pwm1,pwp2) - 2.*dotrc(pwm1,wp2)*dotrr(pwm2,pwp2) )
     &                         + dotrc(pwm1,wm2) * ( dotrc(pwm2,wp2)*dotrr(pwp1,pwp2) - dotrc(pwp1,wp2)*dotrr(pwm2,pwp2) ) ) )
      endif

      if (fm0 .ne. 0) then
        wwww = wwww - fm0 * 2. * gw**2 * wmass**2 *
     &      ( dotcc(wm1,wp2) * ( dotrc(pwm2,wp1)*dotrc(pwp1,wm2) - dotcc(wm2,wp1) * ( dotrr(pwm1,pwp2) + dotrr(pwm2,pwp1) ) )
     &      + dotcc(wm1,wp1) * ( dotrc(pwm2,wp2)*dotrc(pwp2,wm2) - dotcc(wm2,wp2) * ( dotrr(pwm1,pwp1) + dotrr(pwm2,pwp2) ) )
     &      + dotrc(pwm1,wp1) * dotrc(pwp1,wm1) * dotcc(wm2,wp2) + dotrc(pwm1,wp2) * dotrc(pwp2,wm1) * dotcc(wm2,wp1) )
      endif

      if (fm1 .ne. 0) then
        wwww = wwww + fm1 / 2. * gw**2 *wmass**2 *
     &        ( dotcc(wp1,wp2) * ( dotrc(pwm1,wm2) * ( dotrc(pwp1,wm1) + dotrc(pwp2,wm1) )
     &                           + dotrc(pwm2,wm1) * ( dotrc(pwp1,wm2) + dotrc(pwp2,wm2) ) )
     &        + dotcc(wm2,wp2) * ( dotrc(pwp2,wm1) * ( dotrc(pwm1,wp1) - dotrc(pwm2,wp1) )
     &                           + dotrc(pwm2,wp1)*dotrc(pwp1,wm1) - dotrc(pwm2,wm1)*dotrc(pwp2,wp1) )
     &        + dotcc(wm1,wp2) * ( dotrc(pwp2,wm2) * ( dotrc(pwm2,wp1) - dotrc(pwm1,wp1) )
     &                           + dotrc(pwm1,wp1)*dotrc(pwp1,wm2) - dotrc(pwm1,wm2)*dotrc(pwp2,wp1) )
     &        + dotcc(wm2,wp1) * ( dotrc(pwp1,wm1) * ( dotrc(pwm1,wp2) - dotrc(pwm2,wp2) )
     &                           - dotrc(pwm2,wm1)*dotrc(pwp1,wp2) + dotrc(pwm2,wp2)*dotrc(pwp2,wm1) )
     &        - dotcc(wm1,wp1) * ( dotrc(pwp1,wm2) * ( dotrc(pwm1,wp2) - dotrc(pwm2,wp2) )
     &                           + dotrc(pwm1,wm2)*dotrc(pwp1,wp2) - dotrc(pwm1,wp2)*dotrc(pwp2,wm2) )
     &        + dotcc(wm1,wm2) * ( dotrc(pwp1,wp2) * ( dotrc(pwm1,wp1) + dotrc(pwm2,wp1) )
     &                           + dotrc(pwp2,wp1) * ( dotrc(pwm1,wp2) + dotrc(pwm2,wp2) ) )
     &        - dotcc(wm1,wm2)*dotcc(wp1,wp2) * ( dotrr(pwm1,pwp1) + dotrr(pwm1,pwp2) + dotrr(pwm2,pwp1) + dotrr(pwm2,pwp2) )
     &        - dotcc(wm1,wp1)*dotcc(wm2,wp2) * ( dotrr(pwm1,pwp2) + dotrr(pwm2,pwp1) )
     &        - dotcc(wm1,wp2)*dotcc(wm2,wp1) * ( dotrr(pwm1,pwp1) + dotrr(pwm2,pwp2) ) )
      endif

      if (fm6 .ne. 0) then
        wwww = wwww - fm6 * gw**2 * wmass**2 *
     &   ( dotcc(wm1,wp2) * ( dotrc(pwm2,wp1)*dotrc(pwp1,wm2) - dotcc(wm2,wp1) * ( dotrr(pwm1,pwp2) + dotrr(pwm2,pwp1) ) )
     &   + dotcc(wm1,wp1) * ( dotrc(pwm2,wp2)*dotrc(pwp2,wm2) - dotcc(wm2,wp2) * ( dotrr(pwm1,pwp1) + dotrr(pwm2,pwp2) ) )
     &   + dotrc(pwm1,wp1)*dotrc(pwp1,wm1)*dotcc(wm2,wp2) + dotrc(pwm1,wp2)*dotrc(pwp2,wm1)*dotcc(wm2,wp1) )

      endif

      if (fm7 .ne. 0) then
        wwww = wwww - fm7 / 2. * gw**2 * wmass**2 * 
     &   ( dotcc(wp1,wp2) * ( dotrc(pwm1,wm2) * ( dotrc(pwp1,wm1) + dotrc(pwp2,wm1) )
     &                      + dotrc(pwm2,wm1) * ( dotrc(pwp1,wm2) + dotrc(pwp2,wm2) ) )
     &   + dotcc(wm1,wm2) * ( dotrc(pwp1,wp2) * ( dotrc(pwm1,wp1) + dotrc(pwm2,wp1) )
     &                      + dotrc(pwp2,wp1) * ( dotrc(pwm1,wp2) + dotrc(pwm2,wp2) ) )
     &   - dotcc(wm1,wm2)*dotcc(wp1,wp2) * ( dotrr(pwm1,pwp1) + dotrr(pwm1,pwp2) + dotrr(pwm2,pwp1) + dotrr(pwm2,pwp2) )
     &   - dotrc(pwm1,wm2)*dotrc(pwp1,wp2)*dotcc(wm1,wp1) - dotrc(pwm1,wm2)*dotrc(pwp2,wp1)*dotcc(wm1,wp2)
     &   - dotrc(pwm2,wm1)*dotrc(pwp1,wp2)*dotcc(wm2,wp1) - dotrc(pwm2,wm1)*dotrc(pwp2,wp1)*dotcc(wm2,wp2) )
      endif

      if (fwww .ne. 0) then
        wwww = wwww - fwww * 3./2. * gw**4 * 
     &      ( ( dotrc(pwm2,wp1) * (dotrc(pwp1,wm1) - dotrc(pwp2,wm1)) + dotrc(pwm1,wp1) * dotrc(pwp2,wm1)
     &            + dotrc(pwp2,wp1) * dotrc(pwm2,wm1) ) * dotcc(wp2,wm2)
     &      - dotcc(wp1,wm1) * dotcc(wp2,wm2) * (dotrr(pwp1,pwm2) + dotrr(pwm1,pwp2))
     &      + dotcc(wp1,wm2) * ( dotrc(pwm2,wm1) * dotrc(pwp1,wp2) + dotrc(pwp1,wm1) *
     &         (dotrc(pwm1,wp2) - dotrc(pwm2,wp2)) + dotrc(pwp2,wm1) * dotrc(pwm2,wp2) )
     &      - dotcc(wm1,wm2) * ( (dotrc(pwm1,wp1) + dotrc(pwm2,wp1)) * dotrc(pwp1,wp2)
     &         + dotrc(pwp2,wp1) * (dotrc(pwm1,wp2) + dotrc(pwm2,wp2)) )
     &      + dotcc(wm1,wp2) * ( dotrc(pwp2,wp1) * dotrc(pwm1,wm2) + dotrc(pwm1,wp1) *
     &         (dotrc(pwp1,wm2) - dotrc(pwp2,wm2)) + dotrc(pwm2,wp1) * dotrc(pwp2,wm2) )
     &      + dotcc(wp1,wm1) * ( dotrc(pwm2,wp2) * dotrc(pwp1,wm2) + dotrc(pwp1,wp2) * dotrc(pwm1,wm2)
     &         + dotrc(pwm1,wp2) * (dotrc(pwp2,wm2) - dotrc(pwp1,wm2)) )
     &      - dotcc(wp1,wp2) * ( (dotrc(pwp1,wm1) + dotrc(pwp2,wm1)) * dotrc(pwm1,wm2)
     &         + dotrc(pwm2,wm1) * (dotrc(pwp1,wm2) + dotrc(pwp2,wm2)) )
     &      - dotcc(wp1,wm2) * dotcc(wm1,wp2) * (dotrr(pwp1,pwm1) + dotrr(pwp2,pwm2))
     &      + dotcc(wp1,wp2) * dotcc(wm1,wm2) * ( dotrr(pwp1,pwm1) + dotrr(pwp1,pwm2)
     &         + dotrr(pwm1,pwp2) + dotrr(pwp2,pwm2) ) )
      endif

      if (fw .ne. 0) then
        wwww = wwww + fw * gw**2 * wmass**2 * ( 2.*dotcc(wp1,wp2)*dotcc(wm1,wm2) - dotcc(wp1,wm1)*dotcc(wp2,wm2)
     &      - dotcc(wp1,wm2)*dotcc(wp2,wm1) )
      endif

      if (fdw .ne. 0) then
        wwww = wwww + fdw * gw**4 * 2. * (
     &      + dotcc(wp2,wm2) * ( + dotrc(pwm1,wp1) * ( dotrc(pwm2,wm1) - dotrc(pwp1,wm1) )
     &                           + dotrc(pwp2,wp1) * ( 2.*dotrc(pwm2,wm1) + dotrc(pwp1,wm1) - dotrc(pwp2,wm1) )
     &                           - dotrc(pwm2,wp1)*dotrc(pwm2,wm1) )
     &      + dotcc(wm1,wm2) * ( - dotrc(pwp2,wp1) * ( dotrc(pwm1,wp2) + dotrc(pwm2,wp2) - 2.*dotrc(pwp1,wp2) )
     &                           + dotrc(pwm1,wp1) * ( 2.*dotrc(pwm1,wp2) - 2.*dotrc(pwm2,wp2) - dotrc(pwp1,wp2) )
     &                           - dotrc(pwm2,wp1) * ( 2.*dotrc(pwm1,wp2) - 2.*dotrc(pwm2,wp2) + dotrc(pwp1,wp2) ) )
     &      + dotcc(wp1,wm2) * ( + dotrc(pwm2,wm1) * ( dotrc(pwm1,wp2) - dotrc(pwm2,wp2) + 2.*dotrc(pwp1,wp2) )
     &                           + dotrc(pwp2,wm1) * ( dotrc(pwp1,wp2) - dotrc(pwm1,wp2) )
     &                           - dotrc(pwp1,wm1)*dotrc(pwp1,wp2) )
     &      - dotcc(wp1,wm1) * ( + dotrc(pwm2,wp2) * ( dotrc(pwp2,wm2) - dotrc(pwm1,wm2) )
     &                           + dotrc(pwp1,wp2) * ( - 2.*dotrc(pwm1,wm2) + dotrc(pwp1,wm2) - dotrc(pwp2,wm2) )
     &                           + dotrc(pwm1,wp2)*dotrc(pwm1,wm2) )
     &      - dotcc(wm1,wp2) * ( + dotrc(pwm2,wp1) * ( dotrc(pwp1,wm2) - dotrc(pwm1,wm2) )
     &                           + dotrc(pwp2,wp1) * ( - 2.*dotrc(pwm1,wm2) - dotrc(pwp1,wm2) + dotrc(pwp2,wm2) )
     &                           + dotrc(pwm1,wp1)*dotrc(pwm1,wm2) )
     &      + dotcc(wp1,wp2) * ( + dotrc(pwm2,wm1) * ( 2.*dotrc(pwm1,wm2) - dotrc(pwp1,wm2) - dotrc(pwp2,wm2) )
     &                           - dotrc(pwp2,wm1) * ( dotrc(pwm1,wm2) + 2.*dotrc(pwp1,wm2) - 2.*dotrc(pwp2,wm2) )
     &                           - dotrc(pwp1,wm1) * ( dotrc(pwm1,wm2) - 2.*dotrc(pwp1,wm2) + 2.*dotrc(pwp2,wm2) ) )
     &      + dotcc(wp1,wm1)*dotcc(wp2,wm2) * ( dotrr(pwm1,pwm2) + 2.*dotrr(pwm1,pwp1) + 2.*dotrr(pwm2,pwp2) + dotrr(pwp1,pwp2) )
     &      + dotcc(wp1,wm2)*dotcc(wm1,wp2) * ( dotrr(pwm1,pwm2) + 2.*dotrr(pwm1,pwp2) + 2.*dotrr(pwm2,pwp1) + dotrr(pwp1,pwp2) )
     &      - dotcc(wp1,wp2)*dotcc(wm1,wm2) * ( 4.*dotrr(pwm1,pwm2) + dotrr(pwm1,pwp1) + dotrr(pwm1,pwp2) 
     &                                          + dotrr(pwm2,pwp1) + dotrr(pwm2,pwp2) + 4.*dotrr(pwp1,pwp2) ) )
      endif

      if (fwwwtilde .ne. 0) then
        wwww = wwww + fwwwtilde / 2. * gw**4 * 
     &      ( dotcc(wp1,wm1) * ( epsccrr(wp2,wm2,pwm1,pwp1) + epsccrr(wp2,wm2,pwm1,pwp2) + epsccrr(wp2,wm2,pwm2,pwp1) )
     &      + dotcc(wp1,wp2) * ( epsccrr(wm1,wm2,pwm1,pwp1) + epsccrr(wm1,wm2,pwm1,pwp2) 
     &                           - epsccrr(wm1,wm2,pwm2,pwp1) - epsccrr(wm1,wm2,pwm2,pwp2) )
     &      - dotcc(wp1,wm2) * ( epsccrr(wm1,wp2,pwm1,pwp1) + epsccrr(wm1,wp2,pwm2,pwp1) + epsccrr(wm1,wp2,pwm2,pwp2) )
     &      + dotcc(wm1,wp2) * ( epsccrr(wp1,wm2,pwm1,pwp1) + epsccrr(wp1,wm2,pwm1,pwp2) + epsccrr(wp1,wm2,pwm2,pwp2) )
     &      + dotcc(wm1,wm2) * ( - epsccrr(wp1,wp2,pwm1,pwp1) + epsccrr(wp1,wp2,pwm1,pwp2) 
     &                           - epsccrr(wp1,wp2,pwm2,pwp1) + epsccrr(wp1,wp2,pwm2,pwp2) )
     &      + dotcc(wp2,wm2) * ( epsccrr(wp1,wm1,pwm1,pwp2) + epsccrr(wp1,wm1,pwm2,pwp1) + epsccrr(wp1,wm1,pwm2,pwp2) )
     &      - ( dotrc(pwm1,wp1) - dotrc(pwm2,wp1) ) * ( epscccr(wm1,wp2,wm2,pwp1) + epscccr(wm1,wp2,wm2,pwp2) )
     &      + ( dotrc(pwp1,wm1) - dotrc(pwp2,wm1) ) * ( epscccr(wp1,wp2,wm2,pwm1) + epscccr(wp1,wp2,wm2,pwm2) )
     &      + ( dotrc(pwm1,wp2) - dotrc(pwm2,wp2) ) * ( epscccr(wp1,wm1,wm2,pwp1) + epscccr(wp1,wm1,wm2,pwp2) )
     &      - ( dotrc(pwp1,wm2) - dotrc(pwp2,wm2) ) * ( epscccr(wp1,wm1,wp2,pwm1) + epscccr(wp1,wm1,wp2,pwm2) )
     &      + epscccc(wp1,wm1,wp2,wm2) * ( - dotrr(pwm1,pwp1) + dotrr(pwm1,pwp2) + dotrr(pwm2,pwp1) - dotrr(pwm2,pwp2) )
     &      + dotrc(pwp1,wp2) * ( epscccr(wp1,wm1,wm2,pwm1) - epscccr(wp1,wm1,wm2,pwm2) )
     &      + dotrc(pwp2,wp1) * ( epscccr(wm1,wp2,wm2,pwm2) - epscccr(wm1,wp2,wm2,pwm1) )
     &      + dotrc(pwm1,wm2) * ( epscccr(wp1,wm1,wp2,pwp2) - epscccr(wp1,wm1,wp2,pwp1) )
     &      + dotrc(pwm2,wm1) * ( epscccr(wp1,wp2,wm2,pwp1) - epscccr(wp1,wp2,wm2,pwp2) ) )
      endif

      if (fdwtilde .ne. 0) then
        wwww = wwww + fdwtilde * gw**4 * 2. * (
     &      + 2.*dotcc(wp1,wp2)*epsccrr(wm1,wm2,pwm1,pwm2) + dotcc(wp2,wm2)*epsccrr(wp1,wm1,pwm1,pwp1)
     &      - dotcc(wp1,wm2)*epsccrr(wm1,wp2,pwm1,pwp2) + dotcc(wm1,wp2)*epsccrr(wp1,wm2,pwm2,pwp1)
     &      + dotcc(wp1,wm1)*epsccrr(wp2,wm2,pwm2,pwp2) + 2.*dotcc(wm1,wm2)*epsccrr(wp1,wp2,pwp1,pwp2)
     &      + epscccc(wp1,wm1,wp2,wm2) * ( dotrr(pwm1,pwp1) - dotrr(pwm1,pwp2) - dotrr(pwm2,pwp1) + dotrr(pwm2,pwp2) )
     &      + ( dotrc(pwm1,wp2) - dotrc(pwm2,wp2) - dotrc(pwp1,wp2) ) * epscccr(wp1,wm1,wm2,pwm1)
     &      + ( dotrc(pwm1,wp2) - dotrc(pwm2,wp2) + dotrc(pwp1,wp2) ) * epscccr(wp1,wm1,wm2,pwm2)
     &      - ( dotrc(pwm1,wp1) - dotrc(pwm2,wp1) + dotrc(pwp2,wp1) ) * epscccr(wm1,wp2,wm2,pwm2)
     &      + ( - dotrc(pwm1,wp1) + dotrc(pwm2,wp1) + dotrc(pwp2,wp1) ) * epscccr(wm1,wp2,wm2,pwm1)
     &      - ( dotrc(pwm1,wm2) + dotrc(pwp1,wm2) - dotrc(pwp2,wm2) ) * epscccr(wp1,wm1,wp2,pwp2)
     &      + ( dotrc(pwm1,wm2) - dotrc(pwp1,wm2) + dotrc(pwp2,wm2) ) * epscccr(wp1,wm1,wp2,pwp1)
     &      + ( dotrc(pwm2,wm1) + dotrc(pwp1,wm1) - dotrc(pwp2,wm1) ) * epscccr(wp1,wp2,wm2,pwp2)
     &      - ( dotrc(pwm2,wm1) - dotrc(pwp1,wm1) + dotrc(pwp2,wm1) ) * epscccr(wp1,wp2,wm2,pwp1) )
      endif

      return
      end


*********************************************************************************************************
*********************************************************************************************************
CCCCCCCCCCCCCCC
C ZZZZ-vertex C
CCCCCCCCCCCCCCC

      subroutine zzzz_anomal4(z1,z2,z3,z4, zzzz)
c
c This subroutine computes the zzzz-vertex with anomalous couplings
c where all the bosons are defined by the flowing-OUT quantum number.
c     
      implicit none
      double complex z1(6),z2(6),z3(6),z4(6),zzzz 
      double complex dotcc, dotrc, epsccrr, epscccr, epscccc
      double precision pz1(0:3),pz2(0:3),pz3(0:3),pz4(0:3), dotrr
      double precision cos2w, cosw, sinw

      external dotcc, dotrc, dotrr, epsccrr, epscccr, epscccc

#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"   

      cos2w = 1.0d0 - sin2w
      cosw = sqrt(cos2w)
      sinw = sqrt(sin2w)


      pz1(0) = dble( z1(5))
      pz1(1) = dble( z1(6))
      pz1(2) = dimag(z1(6))
      pz1(3) = dimag(z1(5))
      pz2(0) = dble( z2(5))
      pz2(1) = dble( z2(6))
      pz2(2) = dimag(z2(6))
      pz2(3) = dimag(z2(5))
      pz3(0) = dble( z3(5))
      pz3(1) = dble( z3(6))
      pz3(2) = dimag(z3(6))
      pz3(3) = dimag(z3(5))
      pz4(0) = dble( z4(5))
      pz4(1) = dble( z4(6))
      pz4(2) = dimag(z4(6))
      pz4(3) = dimag(z4(5))


      zzzz = dcmplx(0.,0.)


      if (kmatrix) then        
        if (fs2 .ne. fs0) then
          write(*,*) "Sorry, K-matrix is only available for fs0 == fs2"
          stop
        endif
       if (fs0 .ne. 0 .or. fs1 .ne. 0 .or. fs2 .ne. 0) then
        zzzz = zzzz 
     &         + 2. * wmass**4 / cos2w**2 * fsz1z4 * dotcc(z1,z4)*dotcc(z2,z3)
     &         + 2. * wmass**4 / cos2w**2 * fsz1z3 * dotcc(z1,z3)*dotcc(z2,z4)
     &         + 2. * wmass**4 / cos2w**2 * fsz1z2 * dotcc(z1,z2)*dotcc(z3,z4)   
       endif

      else
       if (fs0 .ne. 0) then
        zzzz = zzzz + 2. * fs0 * wmass**4 / cos2w**2 * 
     &      ( dotcc(z1,z4)*dotcc(z2,z3) + dotcc(z1,z3)*dotcc(z2,z4) + dotcc(z1,z2)*dotcc(z3,z4) )
       endif

       if (fs1 .ne. 0) then     
         zzzz = zzzz + 2. * fs1 * wmass**4 / cos2w**2 
     &      * ( dotcc(z1,z4)*dotcc(z2,z3) + dotcc(z1,z3)*dotcc(z2,z4) + dotcc(z1,z2)*dotcc(z3,z4) )
       endif

       if (fs2 .ne. 0) then
        zzzz = zzzz + 2. * fs2 * wmass**4 / cos2w**2 * 
     &      ( dotcc(z1,z4)*dotcc(z2,z3) + dotcc(z1,z3)*dotcc(z2,z4) + dotcc(z1,z2)*dotcc(z3,z4) )
       endif
      endif

      if (fm0 .ne. 0) then
        zzzz = zzzz + 2. * fm0 * gw**2 * wmass**2 *
     &    ( dotcc(z1,z4) * ( dotcc(z2,z3) * ( dotrr(pz1,pz4) + dotrr(pz2,pz3) ) - dotrc(pz2,z3)*dotrc(pz3,z2) )
     &    + dotcc(z1,z3) * ( dotcc(z2,z4) * ( dotrr(pz1,pz3) + dotrr(pz2,pz4) ) - dotrc(pz2,z4)*dotrc(pz4,z2) )
     &    + dotcc(z1,z2) * ( dotcc(z3,z4) * ( dotrr(pz1,pz2) + dotrr(pz3,pz4) ) - dotrc(pz3,z4)*dotrc(pz4,z3) )
     &    - dotrc(pz1,z2)*dotrc(pz2,z1)*dotcc(z3,z4)
     &    - dotrc(pz1,z3)*dotrc(pz3,z1)*dotcc(z2,z4)
     &    - dotrc(pz1,z4)*dotrc(pz4,z1)*dotcc(z2,z3) )
      endif

      if (fm1 .ne. 0) then
        zzzz = zzzz + fm1 / 2. * gw**2 * wmass**2 *
     &    ( dotrc(pz3,z1) * ( dotcc(z2,z3) * ( dotrc(pz1,z4) - dotrc(pz2,z4) )
     &                      + dotcc(z3,z4) * ( dotrc(pz1,z2) - dotrc(pz4,z2) )
     &                      + dotcc(z2,z4) * ( dotrc(pz2,z3) + dotrc(pz4,z3) ) )
     &    + dotrc(pz4,z1) * ( dotcc(z2,z4) * ( dotrc(pz1,z3) - dotrc(pz2,z3) )
     &                      + dotcc(z3,z4) * ( dotrc(pz1,z2) - dotrc(pz3,z2) )
     &                      + dotcc(z2,z3) * ( dotrc(pz2,z4) + dotrc(pz3,z4) ) )
     &   + dotrc(pz1,z3) * dotrc(pz2,z1) * dotcc(z2,z4)
     &   + dotrc(pz1,z2) * dotrc(pz2,z3) * dotcc(z1,z4)
     &   + dotrc(pz1,z4) * dotrc(pz2,z1) * dotcc(z2,z3)
     &   - dotrc(pz1,z4) * dotrc(pz2,z3) * dotcc(z1,z2)
     &   + dotrc(pz1,z2) * dotrc(pz2,z4) * dotcc(z1,z3)
     &   - dotrc(pz1,z3) * dotrc(pz2,z4) * dotcc(z1,z2)
     &   - dotcc(z1,z4) * dotcc(z2,z3) * dotrr(pz1,pz2)
     &   - dotcc(z1,z3) * dotcc(z2,z4) * dotrr(pz1,pz2)
     &   + dotrc(pz1,z3) * dotrc(pz3,z2) * dotcc(z1,z4)
     &   - dotrc(pz1,z4) * dotrc(pz3,z2) * dotcc(z1,z3)
     &   - dotrc(pz1,z2) * dotrc(pz3,z4) * dotcc(z1,z3)
     &   + dotrc(pz1,z3) * dotrc(pz3,z4) * dotcc(z1,z2)
     &   - dotcc(z1,z2) * dotcc(z3,z4) * dotrr(pz1,pz3)
     &   - dotcc(z1,z4) * dotcc(z2,z3) * dotrr(pz1,pz3)
     &   - dotrc(pz1,z3) * dotrc(pz4,z2) * dotcc(z1,z4)
     &   - dotrc(pz1,z2) * dotrc(pz4,z3) * dotcc(z1,z4)
     &   + dotrc(pz1,z4) * dotrc(pz4,z2) * dotcc(z1,z3)
     &   + dotrc(pz1,z4) * dotrc(pz4,z3) * dotcc(z1,z2)
     &   - dotcc(z1,z2) * dotcc(z3,z4) * dotrr(pz1,pz4)
     &   - dotcc(z1,z3) * dotcc(z2,z4) * dotrr(pz1,pz4)
     &   + dotrc(pz2,z1) * dotrc(pz3,z2) * dotcc(z3,z4)
     &   + dotrc(pz2,z4) * dotrc(pz3,z2) * dotcc(z1,z3)
     &   - dotrc(pz2,z1) * dotrc(pz3,z4) * dotcc(z2,z3)
     &   + dotrc(pz2,z3) * dotrc(pz3,z4) * dotcc(z1,z2)
     &   - dotcc(z1,z2) * dotcc(z3,z4) * dotrr(pz2,pz3)
     &   - dotcc(z1,z3) * dotcc(z2,z4) * dotrr(pz2,pz3)
     &   + dotrc(pz2,z1) * dotrc(pz4,z2) * dotcc(z3,z4)
     &   + dotrc(pz2,z3) * dotrc(pz4,z2) * dotcc(z1,z4)
     &   - dotrc(pz2,z1) * dotrc(pz4,z3) * dotcc(z2,z4)
     &   + dotrc(pz2,z4) * dotrc(pz4,z3) * dotcc(z1,z2)
     &   - dotcc(z1,z2) * dotcc(z3,z4) * dotrr(pz2,pz4)
     &   - dotcc(z1,z4) * dotcc(z2,z3) * dotrr(pz2,pz4)
     &   + dotrc(pz3,z2) * dotrc(pz4,z3) * dotcc(z1,z4)
     &   + dotrc(pz3,z4) * dotrc(pz4,z2) * dotcc(z1,z3)
     &   - dotrr(pz3,pz4) * ( dotcc(z1,z4)*dotcc(z2,z3) + dotcc(z1,z3)*dotcc(z2,z4) ) )
      endif

      if (fm2 .ne. 0) then
        zzzz = zzzz - fm2 * gw**2 * wmass**2 * sin2w**2 / cos2w**2 *
     &    ( dotcc(z1,z4) * ( dotrc(pz2,z3)*dotrc(pz3,z2) - dotcc(z2,z3) * ( dotrr(pz1,pz4) + dotrr(pz2,pz3) ) )
     &    + dotcc(z1,z3) * ( dotrc(pz2,z4)*dotrc(pz4,z2) - dotcc(z2,z4) * ( dotrr(pz1,pz3) + dotrr(pz2,pz4) ) )  
     &    + dotcc(z1,z2) * ( dotrc(pz3,z4)*dotrc(pz4,z3) - dotcc(z3,z4) * ( dotrr(pz1,pz2) + dotrr(pz3,pz4) ) )  
     &    + dotrc(pz1,z2) * dotrc(pz2,z1) * dotcc(z3,z4)
     &    + dotrc(pz1,z3) * dotrc(pz3,z1) * dotcc(z2,z4)
     &    + dotrc(pz1,z4) * dotrc(pz4,z1) * dotcc(z2,z3) )
      endif

      if (fm3 .ne. 0) then
        zzzz = zzzz + fm3 / 4. * gw**2 * wmass**2 * sin2w**2 / cos2w**2 *
     &    ( dotrc(pz3,z1) * ( dotcc(z2,z3) * ( dotrc(pz1,z4) - dotrc(pz2,z4) )
     &                      + dotcc(z3,z4) * ( dotrc(pz1,z2) - dotrc(pz4,z2) )
     &                      + dotcc(z2,z4) * ( dotrc(pz2,z3) + dotrc(pz4,z3) ) )
     &    + dotrc(pz4,z1) * ( dotcc(z2,z4) * ( dotrc(pz1,z3) - dotrc(pz2,z3) )
     &                      + dotcc(z3,z4) * ( dotrc(pz1,z2) - dotrc(pz3,z2) )
     &                      + dotcc(z2,z3) * ( dotrc(pz2,z4) + dotrc(pz3,z4) ) )
     &   + dotrc(pz1,z3) * dotrc(pz2,z1) * dotcc(z2,z4)
     &   + dotrc(pz1,z2) * dotrc(pz2,z3) * dotcc(z1,z4)
     &   + dotrc(pz1,z4) * dotrc(pz2,z1) * dotcc(z2,z3)
     &   - dotrc(pz1,z4) * dotrc(pz2,z3) * dotcc(z1,z2)
     &   + dotrc(pz1,z2) * dotrc(pz2,z4) * dotcc(z1,z3)
     &   - dotrc(pz1,z3) * dotrc(pz2,z4) * dotcc(z1,z2)
     &   - dotcc(z1,z4) * dotcc(z2,z3) * dotrr(pz1,pz2)
     &   - dotcc(z1,z3) * dotcc(z2,z4) * dotrr(pz1,pz2)
     &   + dotrc(pz1,z3) * dotrc(pz3,z2) * dotcc(z1,z4)
     &   - dotrc(pz1,z4) * dotrc(pz3,z2) * dotcc(z1,z3)
     &   - dotrc(pz1,z2) * dotrc(pz3,z4) * dotcc(z1,z3)
     &   + dotrc(pz1,z3) * dotrc(pz3,z4) * dotcc(z1,z2)
     &   - dotcc(z1,z2) * dotcc(z3,z4) * dotrr(pz1,pz3)
     &   - dotcc(z1,z4) * dotcc(z2,z3) * dotrr(pz1,pz3)
     &   - dotrc(pz1,z3) * dotrc(pz4,z2) * dotcc(z1,z4)
     &   - dotrc(pz1,z2) * dotrc(pz4,z3) * dotcc(z1,z4)
     &   + dotrc(pz1,z4) * dotrc(pz4,z2) * dotcc(z1,z3)
     &   + dotrc(pz1,z4) * dotrc(pz4,z3) * dotcc(z1,z2)
     &   - dotcc(z1,z2) * dotcc(z3,z4) * dotrr(pz1,pz4)
     &   - dotcc(z1,z3) * dotcc(z2,z4) * dotrr(pz1,pz4)
     &   + dotrc(pz2,z1) * dotrc(pz3,z2) * dotcc(z3,z4)
     &   + dotrc(pz2,z4) * dotrc(pz3,z2) * dotcc(z1,z3)
     &   - dotrc(pz2,z1) * dotrc(pz3,z4) * dotcc(z2,z3)
     &   + dotrc(pz2,z3) * dotrc(pz3,z4) * dotcc(z1,z2)
     &   - dotcc(z1,z2) * dotcc(z3,z4) * dotrr(pz2,pz3)
     &   - dotcc(z1,z3) * dotcc(z2,z4) * dotrr(pz2,pz3)
     &   + dotrc(pz2,z1) * dotrc(pz4,z2) * dotcc(z3,z4)
     &   + dotrc(pz2,z3) * dotrc(pz4,z2) * dotcc(z1,z4)
     &   - dotrc(pz2,z1) * dotrc(pz4,z3) * dotcc(z2,z4)
     &   + dotrc(pz2,z4) * dotrc(pz4,z3) * dotcc(z1,z2)
     &   - dotcc(z1,z2) * dotcc(z3,z4) * dotrr(pz2,pz4)
     &   - dotcc(z1,z4) * dotcc(z2,z3) * dotrr(pz2,pz4)
     &   + dotrc(pz3,z2) * dotrc(pz4,z3) * dotcc(z1,z4)
     &   + dotrc(pz3,z4) * dotrc(pz4,z2) * dotcc(z1,z3)
     &   - dotrr(pz3,pz4) * ( dotcc(z1,z4)*dotcc(z2,z3) + dotcc(z1,z3)*dotcc(z2,z4) ) )
      endif

      if (fm4 .ne. 0) then
        zzzz = zzzz - fm4 * gw**2 * wmass**2 * sin2w / cos2w *
     &    ( dotcc(z1,z4) * ( dotrc(pz2,z3)*dotrc(pz3,z2) - dotcc(z2,z3) * ( dotrr(pz1,pz4) + dotrr(pz2,pz3) ) )
     &    + dotcc(z1,z3) * ( dotrc(pz2,z4)*dotrc(pz4,z2) - dotcc(z2,z4) * ( dotrr(pz1,pz3) + dotrr(pz2,pz4) ) )  
     &    + dotcc(z1,z2) * ( dotrc(pz3,z4)*dotrc(pz4,z3) - dotcc(z3,z4) * ( dotrr(pz1,pz2) + dotrr(pz3,pz4) ) )  
     &    + dotrc(pz1,z2) * dotrc(pz2,z1) * dotcc(z3,z4)
     &    + dotrc(pz1,z3) * dotrc(pz3,z1) * dotcc(z2,z4)
     &    + dotrc(pz1,z4) * dotrc(pz4,z1) * dotcc(z2,z3) )
      endif

      if (fm5 .ne. 0) then
        zzzz = zzzz - fm5 / 4. * gw**2 * wmass**2 * sin2w / cos2w *
     &    ( dotrc(pz3,z1) * ( dotcc(z2,z3) * ( dotrc(pz1,z4) - dotrc(pz2,z4) )
     &                      + dotcc(z3,z4) * ( dotrc(pz1,z2) - dotrc(pz4,z2) )
     &                      + dotcc(z2,z4) * ( dotrc(pz2,z3) + dotrc(pz4,z3) ) )
     &    + dotrc(pz4,z1) * ( dotcc(z2,z4) * ( dotrc(pz1,z3) - dotrc(pz2,z3) )
     &                      + dotcc(z3,z4) * ( dotrc(pz1,z2) - dotrc(pz3,z2) )
     &                      + dotcc(z2,z3) * ( dotrc(pz2,z4) + dotrc(pz3,z4) ) )
     &   + dotrc(pz1,z3) * dotrc(pz2,z1) * dotcc(z2,z4)
     &   + dotrc(pz1,z2) * dotrc(pz2,z3) * dotcc(z1,z4)
     &   + dotrc(pz1,z4) * dotrc(pz2,z1) * dotcc(z2,z3)
     &   - dotrc(pz1,z4) * dotrc(pz2,z3) * dotcc(z1,z2)
     &   + dotrc(pz1,z2) * dotrc(pz2,z4) * dotcc(z1,z3)
     &   - dotrc(pz1,z3) * dotrc(pz2,z4) * dotcc(z1,z2)
     &   - dotcc(z1,z4) * dotcc(z2,z3) * dotrr(pz1,pz2)
     &   - dotcc(z1,z3) * dotcc(z2,z4) * dotrr(pz1,pz2)
     &   + dotrc(pz1,z3) * dotrc(pz3,z2) * dotcc(z1,z4)
     &   - dotrc(pz1,z4) * dotrc(pz3,z2) * dotcc(z1,z3)
     &   - dotrc(pz1,z2) * dotrc(pz3,z4) * dotcc(z1,z3)
     &   + dotrc(pz1,z3) * dotrc(pz3,z4) * dotcc(z1,z2)
     &   - dotcc(z1,z2) * dotcc(z3,z4) * dotrr(pz1,pz3)
     &   - dotcc(z1,z4) * dotcc(z2,z3) * dotrr(pz1,pz3)
     &   - dotrc(pz1,z3) * dotrc(pz4,z2) * dotcc(z1,z4)
     &   - dotrc(pz1,z2) * dotrc(pz4,z3) * dotcc(z1,z4)
     &   + dotrc(pz1,z4) * dotrc(pz4,z2) * dotcc(z1,z3)
     &   + dotrc(pz1,z4) * dotrc(pz4,z3) * dotcc(z1,z2)
     &   - dotcc(z1,z2) * dotcc(z3,z4) * dotrr(pz1,pz4)
     &   - dotcc(z1,z3) * dotcc(z2,z4) * dotrr(pz1,pz4)
     &   + dotrc(pz2,z1) * dotrc(pz3,z2) * dotcc(z3,z4)
     &   + dotrc(pz2,z4) * dotrc(pz3,z2) * dotcc(z1,z3)
     &   - dotrc(pz2,z1) * dotrc(pz3,z4) * dotcc(z2,z3)
     &   + dotrc(pz2,z3) * dotrc(pz3,z4) * dotcc(z1,z2)
     &   - dotcc(z1,z2) * dotcc(z3,z4) * dotrr(pz2,pz3)
     &   - dotcc(z1,z3) * dotcc(z2,z4) * dotrr(pz2,pz3)
     &   + dotrc(pz2,z1) * dotrc(pz4,z2) * dotcc(z3,z4)
     &   + dotrc(pz2,z3) * dotrc(pz4,z2) * dotcc(z1,z4)
     &   - dotrc(pz2,z1) * dotrc(pz4,z3) * dotcc(z2,z4)
     &   + dotrc(pz2,z4) * dotrc(pz4,z3) * dotcc(z1,z2)
     &   - dotcc(z1,z2) * dotcc(z3,z4) * dotrr(pz2,pz4)
     &   - dotcc(z1,z4) * dotcc(z2,z3) * dotrr(pz2,pz4)
     &   + dotrc(pz3,z2) * dotrc(pz4,z3) * dotcc(z1,z4)
     &   + dotrc(pz3,z4) * dotrc(pz4,z2) * dotcc(z1,z3)
     &   - dotrr(pz3,pz4) * ( dotcc(z1,z4)*dotcc(z2,z3) + dotcc(z1,z3)*dotcc(z2,z4) ) )
      endif

      if (fm6 .ne. 0) then
        zzzz = zzzz + fm6 * gw**2 * wmass**2 * 
     &    ( dotcc(z1,z4) * ( dotcc(z2,z3) * ( dotrr(pz1,pz4) + dotrr(pz2,pz3) ) - dotrc(pz2,z3)*dotrc(pz3,z2) )
     &    + dotcc(z1,z3) * ( dotcc(z2,z4) * ( dotrr(pz1,pz3) + dotrr(pz2,pz4) ) - dotrc(pz2,z4)*dotrc(pz4,z2) )  
     &    + dotcc(z1,z2) * ( dotcc(z3,z4) * ( dotrr(pz1,pz2) + dotrr(pz3,pz4) ) - dotrc(pz3,z4)*dotrc(pz4,z3) )  
     &    - dotrc(pz1,z2) * dotrc(pz2,z1) * dotcc(z3,z4)
     &    - dotrc(pz1,z3) * dotrc(pz3,z1) * dotcc(z2,z4)
     &    - dotrc(pz1,z4) * dotrc(pz4,z1) * dotcc(z2,z3) )
      endif

      if (fm7 .ne. 0) then
        zzzz = zzzz + fm7 / 4. * gw**2 * wmass**2 * 
     &    ( - dotrc(pz3,z1) * ( dotcc(z2,z3) * ( dotrc(pz1,z4) - dotrc(pz2,z4) )
     &                        + dotcc(z3,z4) * ( dotrc(pz1,z2) - dotrc(pz4,z2) )
     &                        + dotcc(z2,z4) * ( dotrc(pz2,z3) + dotrc(pz4,z3) ) )
     &    - dotrc(pz4,z1) * ( dotcc(z2,z4) * ( dotrc(pz1,z3) - dotrc(pz2,z3) )
     &                      + dotcc(z3,z4) * ( dotrc(pz1,z2) - dotrc(pz3,z2) )
     &                      + dotcc(z2,z3) * ( dotrc(pz2,z4) + dotrc(pz3,z4) ) )
     &   - dotrc(pz1,z3) * dotrc(pz2,z1) * dotcc(z2,z4)
     &   - dotrc(pz1,z2) * dotrc(pz2,z3) * dotcc(z1,z4)
     &   - dotrc(pz1,z4) * dotrc(pz2,z1) * dotcc(z2,z3)
     &   + dotrc(pz1,z4) * dotrc(pz2,z3) * dotcc(z1,z2)
     &   - dotrc(pz1,z2) * dotrc(pz2,z4) * dotcc(z1,z3)
     &   + dotrc(pz1,z3) * dotrc(pz2,z4) * dotcc(z1,z2)
     &   + dotcc(z1,z4) * dotcc(z2,z3) * dotrr(pz1,pz2)
     &   + dotcc(z1,z3) * dotcc(z2,z4) * dotrr(pz1,pz2)
     &   - dotrc(pz1,z3) * dotrc(pz3,z2) * dotcc(z1,z4)
     &   + dotrc(pz1,z4) * dotrc(pz3,z2) * dotcc(z1,z3)
     &   + dotrc(pz1,z2) * dotrc(pz3,z4) * dotcc(z1,z3)
     &   - dotrc(pz1,z3) * dotrc(pz3,z4) * dotcc(z1,z2)
     &   + dotcc(z1,z2) * dotcc(z3,z4) * dotrr(pz1,pz3)
     &   + dotcc(z1,z4) * dotcc(z2,z3) * dotrr(pz1,pz3)
     &   + dotrc(pz1,z3) * dotrc(pz4,z2) * dotcc(z1,z4)
     &   + dotrc(pz1,z2) * dotrc(pz4,z3) * dotcc(z1,z4)
     &   - dotrc(pz1,z4) * dotrc(pz4,z2) * dotcc(z1,z3)
     &   - dotrc(pz1,z4) * dotrc(pz4,z3) * dotcc(z1,z2)
     &   + dotcc(z1,z2) * dotcc(z3,z4) * dotrr(pz1,pz4)
     &   + dotcc(z1,z3) * dotcc(z2,z4) * dotrr(pz1,pz4)
     &   - dotrc(pz2,z1) * dotrc(pz3,z2) * dotcc(z3,z4)
     &   - dotrc(pz2,z4) * dotrc(pz3,z2) * dotcc(z1,z3)
     &   + dotrc(pz2,z1) * dotrc(pz3,z4) * dotcc(z2,z3)
     &   - dotrc(pz2,z3) * dotrc(pz3,z4) * dotcc(z1,z2)
     &   + dotcc(z1,z2) * dotcc(z3,z4) * dotrr(pz2,pz3)
     &   + dotcc(z1,z3) * dotcc(z2,z4) * dotrr(pz2,pz3)
     &   - dotrc(pz2,z1) * dotrc(pz4,z2) * dotcc(z3,z4)
     &   - dotrc(pz2,z3) * dotrc(pz4,z2) * dotcc(z1,z4)
     &   + dotrc(pz2,z1) * dotrc(pz4,z3) * dotcc(z2,z4)
     &   - dotrc(pz2,z4) * dotrc(pz4,z3) * dotcc(z1,z2)
     &   + dotcc(z1,z2) * dotcc(z3,z4) * dotrr(pz2,pz4)
     &   + dotcc(z1,z4) * dotcc(z2,z3) * dotrr(pz2,pz4)
     &   - dotrc(pz3,z2) * dotrc(pz4,z3) * dotcc(z1,z4)
     &   - dotrc(pz3,z4) * dotrc(pz4,z2) * dotcc(z1,z3)
     &   + dotrr(pz3,pz4) * ( dotcc(z1,z4)*dotcc(z2,z3) + dotcc(z1,z3)*dotcc(z2,z4) ) )
      endif

      if (ft0 .ne. 0) then
        zzzz = zzzz + ft0 * 8. * gw**4 * cos2w**2 *
     &    ( ( dotrc(pz1,z4)*dotrc(pz4,z1) - dotcc(z1,z4)*dotrr(pz1,pz4) ) 
     &            * ( dotrc(pz2,z3)*dotrc(pz3,z2) - dotcc(z2,z3)*dotrr(pz2,pz3) )
     &    + ( dotrc(pz1,z3)*dotrc(pz3,z1) - dotcc(z1,z3)*dotrr(pz1,pz3) ) 
     &            * ( dotrc(pz2,z4)*dotrc(pz4,z2) - dotcc(z2,z4)*dotrr(pz2,pz4) )  
     &    + ( dotrc(pz1,z2)*dotrc(pz2,z1) - dotcc(z1,z2)*dotrr(pz1,pz2) ) 
     &            * ( dotrc(pz3,z4)*dotrc(pz4,z3) - dotcc(z3,z4)*dotrr(pz3,pz4) ) )
      endif

      if (ft1 .ne. 0) then
        zzzz = zzzz + ft1 * 8. * gw**4 * cos2w**2 *
     &    ( ( dotrc(pz1,z4)*dotrc(pz4,z1) - dotcc(z1,z4)*dotrr(pz1,pz4) ) 
     &            * ( dotrc(pz2,z3)*dotrc(pz3,z2) - dotcc(z2,z3)*dotrr(pz2,pz3) )
     &    + ( dotrc(pz1,z3)*dotrc(pz3,z1) - dotcc(z1,z3)*dotrr(pz1,pz3) ) 
     &            * ( dotrc(pz2,z4)*dotrc(pz4,z2) - dotcc(z2,z4)*dotrr(pz2,pz4) )  
     &    + ( dotrc(pz1,z2)*dotrc(pz2,z1) - dotcc(z1,z2)*dotrr(pz1,pz2) ) 
     &            * ( dotrc(pz3,z4)*dotrc(pz4,z3) - dotcc(z3,z4)*dotrr(pz3,pz4) ) )
      endif

      if (ft2 .ne. 0) then
        zzzz = zzzz + ft2 * 2. * gw**4 * cos2w**2 *
     &    ( dotrc(pz3,z1)*dotrc(pz4,z2)*dotrc(pz2,z3)*dotrc(pz1,z4)
     &    + dotrc(pz2,z1)*dotrc(pz3,z2)*dotrc(pz4,z3)*dotrc(pz1,z4)
     &    + dotrc(pz4,z1)*dotrc(pz3,z2)*dotrc(pz1,z3)*dotrc(pz2,z4)
     &    + dotrc(pz3,z1)*dotrc(pz1,z2)*dotrc(pz4,z3)*dotrc(pz2,z4)
     &    + dotrc(pz2,z1)*dotrc(pz4,z2)*dotrc(pz1,z3)*dotrc(pz3,z4)
     &    + dotrc(pz4,z1)*dotrc(pz1,z2)*dotrc(pz2,z3)*dotrc(pz3,z4)
     &    + dotcc(z1,z2)*dotcc(z3,z4) * ( dotrr(pz1,pz4)*dotrr(pz2,pz3) + dotrr(pz1,pz3)*dotrr(pz2,pz4) )
     &    + dotcc(z3,z4) * ( - dotrc(pz2,z1) * ( dotrc(pz4,z2)*dotrr(pz1,pz3) + dotrc(pz3,z2)*dotrr(pz1,pz4) )
     &                       + dotrc(pz4,z1) * ( dotrc(pz3,z2)*dotrr(pz1,pz2) - dotrc(pz1,z2)*dotrr(pz2,pz3) )
     &                       + dotrc(pz3,z1) * ( dotrc(pz4,z2)*dotrr(pz1,pz2) - dotrc(pz1,z2)*dotrr(pz2,pz4) ) )
     &    - dotcc(z1,z2) * ( dotrc(pz4,z3) * ( dotrc(pz2,z4)*dotrr(pz1,pz3) + dotrc(pz1,z4)*dotrr(pz2,pz3) )
     &                     + dotrc(pz3,z4) * ( dotrc(pz2,z3)*dotrr(pz1,pz4) + dotrc(pz1,z3)*dotrr(pz2,pz4) )
     &                     - dotrr(pz3,pz4) * ( dotrc(pz2,z3)*dotrc(pz1,z4) + dotrc(pz1,z3)*dotrc(pz2,z4) ) )
     &    + dotcc(z1,z3)*dotcc(z2,z4) * ( dotrr(pz1,pz4)*dotrr(pz2,pz3) + dotrr(pz1,pz2)*dotrr(pz3,pz4) )
     &    + dotcc(z1,z4)*dotcc(z2,z3) * ( dotrr(pz1,pz3)*dotrr(pz2,pz4) + dotrr(pz1,pz2)*dotrr(pz3,pz4) )
     &    - dotcc(z2,z4) * ( + dotrc(pz3,z1) * ( dotrc(pz4,z3)*dotrr(pz1,pz2) + dotrc(pz2,z3)*dotrr(pz1,pz4) )
     &                       + dotrc(pz4,z1) * ( dotrc(pz1,z3)*dotrr(pz2,pz3) - dotrc(pz2,z3)*dotrr(pz1,pz3) )
     &                       + dotrc(pz2,z1) * ( dotrc(pz1,z3)*dotrr(pz3,pz4) - dotrc(pz4,z3)*dotrr(pz1,pz3) ) )
     &    - dotcc(z1,z4) * ( + dotrc(pz4,z2) * ( dotrc(pz2,z3)*dotrr(pz1,pz3) - dotrc(pz1,z3)*dotrr(pz2,pz3) )
     &                       + dotrc(pz3,z2) * ( dotrc(pz4,z3)*dotrr(pz1,pz2) + dotrc(pz1,z3)*dotrr(pz2,pz4) )
     &                       + dotrc(pz1,z2) * ( dotrc(pz2,z3)*dotrr(pz3,pz4) - dotrc(pz4,z3)*dotrr(pz2,pz3) ) )
     &    - dotcc(z2,z3) * ( + dotrc(pz4,z1) * ( dotrc(pz3,z4)*dotrr(pz1,pz2) + dotrc(pz2,z4)*dotrr(pz1,pz3) )
     &                       + dotrc(pz3,z1) * ( dotrc(pz1,z4)*dotrr(pz2,pz4) - dotrc(pz2,z4)*dotrr(pz1,pz4) )
     &                       + dotrc(pz2,z1) * ( dotrc(pz1,z4)*dotrr(pz3,pz4) - dotrc(pz3,z4)*dotrr(pz1,pz4) ) )
     &    - dotcc(z1,z3) * ( + dotrc(pz4,z2) * ( dotrc(pz3,z4)*dotrr(pz1,pz2) + dotrc(pz1,z4)*dotrr(pz2,pz3) )
     &                       + dotrc(pz3,z2) * ( dotrc(pz2,z4)*dotrr(pz1,pz4) - dotrc(pz1,z4)*dotrr(pz2,pz4) )
     &                       + dotrc(pz1,z2) * ( dotrc(pz2,z4)*dotrr(pz3,pz4) - dotrc(pz3,z4)*dotrr(pz2,pz4) ) ) )
      endif

      if (ft5 .ne. 0) then
        zzzz = zzzz + ft5 * 4. * gw**4 * sin2w**2 *
     &   ( ( dotrc(pz1,z4)*dotrc(pz4,z1) - dotcc(z1,z4)*dotrr(pz1,pz4) ) 
     &            * ( dotrc(pz2,z3)*dotrc(pz3,z2) - dotcc(z2,z3)*dotrr(pz2,pz3) )
     &   + ( dotrc(pz1,z3)*dotrc(pz3,z1) - dotcc(z1,z3)*dotrr(pz1,pz3) ) 
     &            * ( dotrc(pz2,z4)*dotrc(pz4,z2) - dotcc(z2,z4)*dotrr(pz2,pz4) )
     &   + ( dotrc(pz1,z2)*dotrc(pz2,z1) - dotcc(z1,z2)*dotrr(pz1,pz2) ) 
     &            * ( dotrc(pz3,z4)*dotrc(pz4,z3) - dotcc(z3,z4)*dotrr(pz3,pz4) ) )
      endif

      if (ft6 .ne. 0) then
        zzzz = zzzz + ft6 * 4. * gw**4 * sin2w**2 *
     &   ( ( dotrc(pz1,z4)*dotrc(pz4,z1) - dotcc(z1,z4)*dotrr(pz1,pz4) ) 
     &            * ( dotrc(pz2,z3)*dotrc(pz3,z2) - dotcc(z2,z3)*dotrr(pz2,pz3) )
     &   + ( dotrc(pz1,z3)*dotrc(pz3,z1) - dotcc(z1,z3)*dotrr(pz1,pz3) ) 
     &            * ( dotrc(pz2,z4)*dotrc(pz4,z2) - dotcc(z2,z4)*dotrr(pz2,pz4) )
     &   + ( dotrc(pz1,z2)*dotrc(pz2,z1) - dotcc(z1,z2)*dotrr(pz1,pz2) ) 
     &            * ( dotrc(pz3,z4)*dotrc(pz4,z3) - dotcc(z3,z4)*dotrr(pz3,pz4) ) )
      endif

      if (ft7 .ne. 0) then
        zzzz = zzzz + ft7 * gw**4 * sin2w**2 *
     &    ( dotrc(pz3,z1) * dotrc(pz4,z2) * dotrc(pz2,z3) * dotrc(pz1,z4)
     &    + dotrc(pz2,z1) * dotrc(pz3,z2) * dotrc(pz4,z3) * dotrc(pz1,z4)
     &    + dotrc(pz4,z1) * dotrc(pz3,z2) * dotrc(pz1,z3) * dotrc(pz2,z4)
     &    + dotrc(pz3,z1) * dotrc(pz1,z2) * dotrc(pz4,z3) * dotrc(pz2,z4)
     &    + dotrc(pz2,z1) * dotrc(pz4,z2) * dotrc(pz1,z3) * dotrc(pz3,z4)
     &    + dotrc(pz4,z1) * dotrc(pz1,z2) * dotrc(pz2,z3) * dotrc(pz3,z4)
     &    + dotcc(z1,z2)*dotcc(z3,z4) * ( dotrr(pz1,pz4)*dotrr(pz2,pz3) + dotrr(pz1,pz3)*dotrr(pz2,pz4) )
     &    + dotcc(z3,z4) * ( - dotrc(pz2,z1) * ( dotrc(pz4,z2)*dotrr(pz1,pz3) + dotrc(pz3,z2)*dotrr(pz1,pz4) )
     &                       + dotrc(pz4,z1) * ( dotrc(pz3,z2)*dotrr(pz1,pz2) - dotrc(pz1,z2)*dotrr(pz2,pz3) )
     &                       + dotrc(pz3,z1) * ( dotrc(pz4,z2)*dotrr(pz1,pz2) - dotrc(pz1,z2)*dotrr(pz2,pz4) ) )
     &    - dotcc(z1,z2) * ( dotrc(pz4,z3) * ( dotrc(pz2,z4)*dotrr(pz1,pz3) + dotrc(pz1,z4)*dotrr(pz2,pz3) )
     &                     + dotrc(pz3,z4) * ( dotrc(pz2,z3)*dotrr(pz1,pz4) + dotrc(pz1,z3)*dotrr(pz2,pz4) )
     &                     - dotrr(pz3,pz4) * ( dotrc(pz2,z3)*dotrc(pz1,z4) + dotrc(pz1,z3)*dotrc(pz2,z4) ) )
     &    + dotcc(z1,z3)*dotcc(z2,z4) * ( dotrr(pz1,pz4)*dotrr(pz2,pz3) + dotrr(pz1,pz2)*dotrr(pz3,pz4) ) 
     &    + dotcc(z1,z4)*dotcc(z2,z3) * ( dotrr(pz1,pz3)*dotrr(pz2,pz4) + dotrr(pz1,pz2)*dotrr(pz3,pz4) ) 
     &    - dotcc(z2,z4) * ( dotrc(pz3,z1) * ( dotrc(pz4,z3)*dotrr(pz1,pz2) + dotrc(pz2,z3)*dotrr(pz1,pz4) )
     &                     + dotrc(pz4,z1) * ( dotrc(pz1,z3)*dotrr(pz2,pz3) - dotrc(pz2,z3)*dotrr(pz1,pz3) )
     &                     + dotrc(pz2,z1) * ( dotrc(pz1,z3)*dotrr(pz3,pz4) - dotrc(pz4,z3)*dotrr(pz1,pz3) ) )
     &    - dotcc(z1,z4) * ( dotrc(pz4,z2) * ( dotrc(pz2,z3)*dotrr(pz1,pz3) - dotrc(pz1,z3)*dotrr(pz2,pz3) )
     &                     + dotrc(pz3,z2) * ( dotrc(pz4,z3)*dotrr(pz1,pz2) + dotrc(pz1,z3)*dotrr(pz2,pz4) )
     &                     + dotrc(pz1,z2) * ( dotrc(pz2,z3)*dotrr(pz3,pz4) - dotrc(pz4,z3)*dotrr(pz2,pz3) ) )
     &    - dotcc(z2,z3) * ( dotrc(pz4,z1) * ( dotrc(pz3,z4)*dotrr(pz1,pz2) + dotrc(pz2,z4)*dotrr(pz1,pz3) )
     &                     + dotrc(pz3,z1) * ( dotrc(pz1,z4)*dotrr(pz2,pz4) - dotrc(pz2,z4)*dotrr(pz1,pz4) )
     &                     + dotrc(pz2,z1) * ( dotrc(pz1,z4)*dotrr(pz3,pz4) - dotrc(pz3,z4)*dotrr(pz1,pz4) ) )
     &    - dotcc(z1,z3) * ( dotrc(pz4,z2) * ( dotrc(pz3,z4)*dotrr(pz1,pz2) + dotrc(pz1,z4)*dotrr(pz2,pz3) )
     &                     + dotrc(pz3,z2) * ( dotrc(pz2,z4)*dotrr(pz1,pz4) - dotrc(pz1,z4)*dotrr(pz2,pz4) )
     &                     + dotrc(pz1,z2) * ( dotrc(pz2,z4)*dotrr(pz3,pz4) - dotrc(pz3,z4)*dotrr(pz2,pz4) ) ) )
      endif

      if (ft8 .ne. 0) then
        zzzz = zzzz + ft8 * 2. * gw**4 * sin2w**4 / cos2w**2 *
     &   ( ( dotrc(pz1,z4)*dotrc(pz4,z1) - dotcc(z1,z4)*dotrr(pz1,pz4) ) 
     &            * ( dotrc(pz2,z3)*dotrc(pz3,z2) - dotcc(z2,z3)*dotrr(pz2,pz3) )
     &   + ( dotrc(pz1,z3)*dotrc(pz3,z1) - dotcc(z1,z3)*dotrr(pz1,pz3) ) 
     &            * ( dotrc(pz2,z4)*dotrc(pz4,z2) - dotcc(z2,z4)*dotrr(pz2,pz4) )
     &   + ( dotrc(pz1,z2)*dotrc(pz2,z1) - dotcc(z1,z2)*dotrr(pz1,pz2) ) 
     &            * ( dotrc(pz3,z4)*dotrc(pz4,z3) - dotcc(z3,z4)*dotrr(pz3,pz4) ) )
      endif

      if (ft9 .ne. 0) then
        zzzz = zzzz + ft9 / 2. * gw**4 * sin2w**4 / cos2w**2 * 
     &    ( dotrc(pz3,z1) * dotrc(pz4,z2) * dotrc(pz2,z3) * dotrc(pz1,z4)
     &    + dotrc(pz2,z1) * dotrc(pz3,z2) * dotrc(pz4,z3) * dotrc(pz1,z4)
     &    + dotrc(pz4,z1) * dotrc(pz3,z2) * dotrc(pz1,z3) * dotrc(pz2,z4)
     &    + dotrc(pz3,z1) * dotrc(pz1,z2) * dotrc(pz4,z3) * dotrc(pz2,z4)
     &    + dotrc(pz2,z1) * dotrc(pz4,z2) * dotrc(pz1,z3) * dotrc(pz3,z4)
     &    + dotrc(pz4,z1) * dotrc(pz1,z2) * dotrc(pz2,z3) * dotrc(pz3,z4)
     &    + dotcc(z1,z2)*dotcc(z3,z4) * ( dotrr(pz1,pz4)*dotrr(pz2,pz3) + dotrr(pz1,pz3)*dotrr(pz2,pz4) )
     &    + dotcc(z3,z4) * ( - dotrc(pz2,z1) * ( dotrc(pz4,z2)*dotrr(pz1,pz3) + dotrc(pz3,z2)*dotrr(pz1,pz4) )
     &                       + dotrc(pz4,z1) * ( dotrc(pz3,z2)*dotrr(pz1,pz2) - dotrc(pz1,z2)*dotrr(pz2,pz3) )
     &                       + dotrc(pz3,z1) * ( dotrc(pz4,z2)*dotrr(pz1,pz2) - dotrc(pz1,z2)*dotrr(pz2,pz4) ) )
     &    - dotcc(z1,z2) * ( dotrc(pz4,z3) * ( dotrc(pz2,z4)*dotrr(pz1,pz3) + dotrc(pz1,z4)*dotrr(pz2,pz3) )
     &                     + dotrc(pz3,z4) * ( dotrc(pz2,z3)*dotrr(pz1,pz4) + dotrc(pz1,z3)*dotrr(pz2,pz4) )
     &                     - dotrr(pz3,pz4) * ( dotrc(pz2,z3)*dotrc(pz1,z4) + dotrc(pz1,z3)*dotrc(pz2,z4) ) )
     &    + dotcc(z1,z3)*dotcc(z2,z4) * ( dotrr(pz1,pz4)*dotrr(pz2,pz3) + dotrr(pz1,pz2)*dotrr(pz3,pz4) ) 
     &    + dotcc(z1,z4)*dotcc(z2,z3) * ( dotrr(pz1,pz3)*dotrr(pz2,pz4) + dotrr(pz1,pz2)*dotrr(pz3,pz4) ) 
     &    - dotcc(z2,z4) * ( dotrc(pz3,z1) * ( dotrc(pz4,z3)*dotrr(pz1,pz2) + dotrc(pz2,z3)*dotrr(pz1,pz4) )
     &                     + dotrc(pz4,z1) * ( dotrc(pz1,z3)*dotrr(pz2,pz3) - dotrc(pz2,z3)*dotrr(pz1,pz3) )
     &                     + dotrc(pz2,z1) * ( dotrc(pz1,z3)*dotrr(pz3,pz4) - dotrc(pz4,z3)*dotrr(pz1,pz3) ) )
     &    - dotcc(z1,z4) * ( dotrc(pz4,z2) * ( dotrc(pz2,z3)*dotrr(pz1,pz3) - dotrc(pz1,z3)*dotrr(pz2,pz3) )
     &                     + dotrc(pz3,z2) * ( dotrc(pz4,z3)*dotrr(pz1,pz2) + dotrc(pz1,z3)*dotrr(pz2,pz4) )
     &                     + dotrc(pz1,z2) * ( dotrc(pz2,z3)*dotrr(pz3,pz4) - dotrc(pz4,z3)*dotrr(pz2,pz3) ) )
     &    - dotcc(z2,z3) * ( dotrc(pz4,z1) * ( dotrc(pz3,z4)*dotrr(pz1,pz2) + dotrc(pz2,z4)*dotrr(pz1,pz3) )
     &                     + dotrc(pz3,z1) * ( dotrc(pz1,z4)*dotrr(pz2,pz4) - dotrc(pz2,z4)*dotrr(pz1,pz4) )
     &                     + dotrc(pz2,z1) * ( dotrc(pz1,z4)*dotrr(pz3,pz4) - dotrc(pz3,z4)*dotrr(pz1,pz4) ) )
     &    - dotcc(z1,z3) * ( dotrc(pz4,z2) * ( dotrc(pz3,z4)*dotrr(pz1,pz2) + dotrc(pz1,z4)*dotrr(pz2,pz3) )
     &                     + dotrc(pz3,z2) * ( dotrc(pz2,z4)*dotrr(pz1,pz4) - dotrc(pz1,z4)*dotrr(pz2,pz4) )
     &                     + dotrc(pz1,z2) * ( dotrc(pz2,z4)*dotrr(pz3,pz4) - dotrc(pz3,z4)*dotrr(pz2,pz4) ) ) )
      endif

      return
      end

*********************************************************************************************************
CCCCCCCCCCCCCCC
C AZZZ-vertex C
CCCCCCCCCCCCCCC

      subroutine azzz_anomal4(a,z1,z2,z3, azzz)
c
c This subroutine computes the azzz-vertex with anomalous couplings
c where all the bosons are defined by the flowing-OUT quantum number.
c     
      implicit none
      double complex z1(6),z2(6),z3(6),a(6),azzz 
      double complex dotcc, dotrc, epsccrr, epscccr, epscccc
      double precision pz1(0:3),pz2(0:3),pz3(0:3),pa(0:3), dotrr
      double precision cos2w, cosw, sinw

      external dotcc, dotrc, dotrr, epsccrr, epscccr, epscccc

#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"   

      cos2w = 1.0d0 - sin2w
      cosw = sqrt(cos2w)
      sinw = sqrt(sin2w)


      pz1(0) = dble( z1(5))
      pz1(1) = dble( z1(6))
      pz1(2) = dimag(z1(6))
      pz1(3) = dimag(z1(5))
      pz2(0) = dble( z2(5))
      pz2(1) = dble( z2(6))
      pz2(2) = dimag(z2(6))
      pz2(3) = dimag(z2(5))
      pz3(0) = dble( z3(5))
      pz3(1) = dble( z3(6))
      pz3(2) = dimag(z3(6))
      pz3(3) = dimag(z3(5))
      pa(0) = dble( a(5))
      pa(1) = dble( a(6))
      pa(2) = dimag(a(6))
      pa(3) = dimag(a(5))


      azzz = dcmplx(0.,0.)

      if (fm0 .ne. 0) then
        azzz = azzz - 2. * fm0 * gw**2 * wmass**2 * sinw / cosw *
     &    ( dotcc(z2,z3) * ( dotrc(pz1,a)*dotrc(pa,z1) - dotcc(a,z1)*dotrr(pa,pz1) )
     &    + dotcc(z1,z3) * ( dotrc(pz2,a)*dotrc(pa,z2) - dotcc(a,z2)*dotrr(pa,pz2) )
     &    + dotcc(z1,z2) * ( dotrc(pz3,a)*dotrc(pa,z3) - dotcc(a,z3)*dotrr(pa,pz3) ) )
      endif

      if (fm1 .ne. 0) then
        azzz = azzz + fm1 / 2. * gw**2 * wmass**2 * sinw / cosw *
     &    ( - dotrc(pa,z2) * dotrc(pz1,z3) * dotcc(a,z1)
     &      + dotrc(pa,z1) * dotrc(pz1,z2) * dotcc(a,z3)
     &      - dotrc(pa,z3) * dotrc(pz1,z2) * dotcc(a,z1)
     &      + dotrc(pz1,a) * ( dotrc(pa,z2) * dotcc(z1,z3) + dotrc(pa,z3) * dotcc(z1,z2) )
     &      + dotrc(pa,z1) * dotrc(pz1,z3) * dotcc(a,z2)
     &      - dotcc(a,z3) * dotcc(z1,z2) * dotrr(pa,pz1)
     &      - dotcc(a,z2) * dotcc(z1,z3) * dotrr(pa,pz1)
     &      + dotrc(pa,z2) * dotrc(pz2,z1) * dotcc(a,z3)
     &      + dotrc(pa,z2) * dotrc(pz2,z3) * dotcc(a,z1)
     &      - dotrc(pa,z3) * dotrc(pz2,z1) * dotcc(a,z2)
     &      + dotrc(pz2,a) * ( dotrc(pa,z1) * dotcc(z2,z3) + dotrc(pa,z3) * dotcc(z1,z2) )
     &      - dotrc(pa,z1) * dotrc(pz2,z3) * dotcc(a,z2)
     &      - dotcc(a,z3) * dotcc(z1,z2) * dotrr(pa,pz2)
     &      - dotcc(a,z1) * dotcc(z2,z3) * dotrr(pa,pz2)
     &      - dotrc(pa,z2) * dotrc(pz3,z1) * dotcc(a,z3)
     &      + dotrc(pz3,a) * ( dotrc(pa,z1) * dotcc(z2,z3) + dotrc(pa,z2) * dotcc(z1,z3) )
     &      - dotrc(pa,z1) * dotrc(pz3,z2) * dotcc(a,z3)
     &      + dotrc(pa,z3) * dotrc(pz3,z1) * dotcc(a,z2)
     &      + dotrc(pa,z3) * dotrc(pz3,z2) * dotcc(a,z1)
     &      - dotrr(pa,pz3) * ( dotcc(a,z2) * dotcc(z1,z3) + dotcc(a,z1) * dotcc(z2,z3) ) )
      endif

      if (fm2 .ne. 0) then
        azzz = azzz + fm2 * gw**2 * wmass**2 * sinw**3 / cosw**3 *
     &    ( dotcc(z2,z3) * ( dotrc(pz1,a)*dotrc(pa,z1) - dotcc(a,z1)*dotrr(pa,pz1) )
     &    + dotcc(z1,z3) * ( dotrc(pz2,a)*dotrc(pa,z2) - dotcc(a,z2)*dotrr(pa,pz2) )
     &    + dotcc(z1,z2) * ( dotrc(pz3,a)*dotrc(pa,z3) - dotcc(a,z3)*dotrr(pa,pz3) ) )
      endif

      if (fm3 .ne. 0) then
        azzz = azzz - fm3 / 4. * gw**2 * wmass**2 * sinw**3 / cosw**3 *
     &    ( - dotrc(pa,z2) * dotrc(pz1,z3) * dotcc(a,z1)
     &      + dotrc(pa,z1) * dotrc(pz1,z2) * dotcc(a,z3)
     &      - dotrc(pa,z3) * dotrc(pz1,z2) * dotcc(a,z1)
     &      + dotrc(pz1,a) * ( dotrc(pa,z2) * dotcc(z1,z3) + dotrc(pa,z3) * dotcc(z1,z2) )
     &      + dotrc(pa,z1) * dotrc(pz1,z3) * dotcc(a,z2)
     &      - dotcc(a,z3) * dotcc(z1,z2) * dotrr(pa,pz1)
     &      - dotcc(a,z2) * dotcc(z1,z3) * dotrr(pa,pz1)
     &      + dotrc(pa,z2) * dotrc(pz2,z1) * dotcc(a,z3)
     &      + dotrc(pa,z2) * dotrc(pz2,z3) * dotcc(a,z1)
     &      - dotrc(pa,z3) * dotrc(pz2,z1) * dotcc(a,z2)
     &      + dotrc(pz2,a) * ( dotrc(pa,z1) * dotcc(z2,z3) + dotrc(pa,z3) * dotcc(z1,z2) )
     &      - dotrc(pa,z1) * dotrc(pz2,z3) * dotcc(a,z2)
     &      - dotcc(a,z3) * dotcc(z1,z2) * dotrr(pa,pz2)
     &      - dotcc(a,z1) * dotcc(z2,z3) * dotrr(pa,pz2)
     &      - dotrc(pa,z2) * dotrc(pz3,z1) * dotcc(a,z3)
     &      + dotrc(pz3,a) * ( dotrc(pa,z1) * dotcc(z2,z3) + dotrc(pa,z2) * dotcc(z1,z3) )
     &      - dotrc(pa,z1) * dotrc(pz3,z2) * dotcc(a,z3)
     &      + dotrc(pa,z3) * dotrc(pz3,z1) * dotcc(a,z2)
     &      + dotrc(pa,z3) * dotrc(pz3,z2) * dotcc(a,z1)
     &      - dotrr(pa,pz3) * ( dotcc(a,z2) * dotcc(z1,z3) + dotcc(a,z1) * dotcc(z2,z3) ) )
      endif

      if (fm4 .ne. 0) then
        azzz = azzz - fm4 / 4. * gw**2 * wmass**2 * ( sinw**3/cosw**3 + sinw/cosw**3 - 3.*sinw/cosw ) *
     &    ( dotcc(z2,z3) * ( dotrc(pz1,a)*dotrc(pa,z1) - dotcc(a,z1)*dotrr(pa,pz1) )
     &    + dotcc(z1,z3) * ( dotrc(pz2,a)*dotrc(pa,z2) - dotcc(a,z2)*dotrr(pa,pz2) )
     &    + dotcc(z1,z2) * ( dotrc(pz3,a)*dotrc(pa,z3) - dotcc(a,z3)*dotrr(pa,pz3) ) )
      endif

      if (fm5 .ne. 0) then
        azzz = azzz - fm5 / 16. * gw**2 * wmass**2 * ( sinw**3/cosw**3 + sinw/cosw**3 - 3.*sinw/cosw ) *
     &    ( - dotrc(pa,z2) * dotrc(pz1,z3) * dotcc(a,z1)
     &      + dotrc(pa,z1) * dotrc(pz1,z2) * dotcc(a,z3)
     &      - dotrc(pa,z3) * dotrc(pz1,z2) * dotcc(a,z1)
     &      + dotrc(pz1,a) * ( dotrc(pa,z2) * dotcc(z1,z3) + dotrc(pa,z3) * dotcc(z1,z2) )
     &      + dotrc(pa,z1) * dotrc(pz1,z3) * dotcc(a,z2)
     &      - dotcc(a,z3) * dotcc(z1,z2) * dotrr(pa,pz1)
     &      - dotcc(a,z2) * dotcc(z1,z3) * dotrr(pa,pz1)
     &      + dotrc(pa,z2) * dotrc(pz2,z1) * dotcc(a,z3)
     &      + dotrc(pa,z2) * dotrc(pz2,z3) * dotcc(a,z1)
     &      - dotrc(pa,z3) * dotrc(pz2,z1) * dotcc(a,z2)
     &      + dotrc(pz2,a) * ( dotrc(pa,z1) * dotcc(z2,z3) + dotrc(pa,z3) * dotcc(z1,z2) )
     &      - dotrc(pa,z1) * dotrc(pz2,z3) * dotcc(a,z2)
     &      - dotcc(a,z3) * dotcc(z1,z2) * dotrr(pa,pz2)
     &      - dotcc(a,z1) * dotcc(z2,z3) * dotrr(pa,pz2)
     &      - dotrc(pa,z2) * dotrc(pz3,z1) * dotcc(a,z3)
     &      + dotrc(pz3,a) * ( dotrc(pa,z1) * dotcc(z2,z3) + dotrc(pa,z2) * dotcc(z1,z3) )
     &      - dotrc(pa,z1) * dotrc(pz3,z2) * dotcc(a,z3)
     &      + dotrc(pa,z3) * dotrc(pz3,z1) * dotcc(a,z2)
     &      + dotrc(pa,z3) * dotrc(pz3,z2) * dotcc(a,z1)
     &      - dotrr(pa,pz3) * ( dotcc(a,z2) * dotcc(z1,z3) + dotcc(a,z1) * dotcc(z2,z3) ) )
      endif

      if (fm6 .ne. 0) then
        azzz = azzz - fm6 * gw**2 * wmass**2 * sinw / cosw *
     &    ( dotcc(z2,z3) * ( dotrc(pz1,a)*dotrc(pa,z1) - dotcc(a,z1)*dotrr(pa,pz1) )
     &    + dotcc(z1,z3) * ( dotrc(pz2,a)*dotrc(pa,z2) - dotcc(a,z2)*dotrr(pa,pz2) )
     &    + dotcc(z1,z2) * ( dotrc(pz3,a)*dotrc(pa,z3) - dotcc(a,z3)*dotrr(pa,pz3) ) )
      endif

      if (fm7 .ne. 0) then
        azzz = azzz - fm7 / 4. * gw**2 * wmass**2 * sinw / cosw *
     &    ( - dotrc(pa,z2) * dotrc(pz1,z3) * dotcc(a,z1)
     &      + dotrc(pa,z1) * dotrc(pz1,z2) * dotcc(a,z3)
     &      - dotrc(pa,z3) * dotrc(pz1,z2) * dotcc(a,z1)
     &      + dotrc(pz1,a) * ( dotrc(pa,z2) * dotcc(z1,z3) + dotrc(pa,z3) * dotcc(z1,z2) )
     &      + dotrc(pa,z1) * dotrc(pz1,z3) * dotcc(a,z2)
     &      - dotcc(a,z3) * dotcc(z1,z2) * dotrr(pa,pz1)
     &      - dotcc(a,z2) * dotcc(z1,z3) * dotrr(pa,pz1)
     &      + dotrc(pa,z2) * dotrc(pz2,z1) * dotcc(a,z3)
     &      + dotrc(pa,z2) * dotrc(pz2,z3) * dotcc(a,z1)
     &      - dotrc(pa,z3) * dotrc(pz2,z1) * dotcc(a,z2)
     &      + dotrc(pz2,a) * ( dotrc(pa,z1) * dotcc(z2,z3) + dotrc(pa,z3) * dotcc(z1,z2) )
     &      - dotrc(pa,z1) * dotrc(pz2,z3) * dotcc(a,z2)
     &      - dotcc(a,z3) * dotcc(z1,z2) * dotrr(pa,pz2)
     &      - dotcc(a,z1) * dotcc(z2,z3) * dotrr(pa,pz2)
     &      - dotrc(pa,z2) * dotrc(pz3,z1) * dotcc(a,z3)
     &      + dotrc(pz3,a) * ( dotrc(pa,z1) * dotcc(z2,z3) + dotrc(pa,z2) * dotcc(z1,z3) )
     &      - dotrc(pa,z1) * dotrc(pz3,z2) * dotcc(a,z3)
     &      + dotrc(pa,z3) * dotrc(pz3,z1) * dotcc(a,z2)
     &      + dotrc(pa,z3) * dotrc(pz3,z2) * dotcc(a,z1)
     &      - dotrr(pa,pz3) * ( dotcc(a,z2) * dotcc(z1,z3) + dotcc(a,z1) * dotcc(z2,z3) ) )
      endif

      if (ft0 .ne. 0) then
        azzz = azzz + ft0 * 8. * gw**4 * sinw * cosw**3 *
     &    ( ( dotrc(pz3,a)*dotrc(pa,z3) - dotcc(a,z3)*dotrr(pa,pz3) ) 
     &            * ( dotrc(pz1,z2)*dotrc(pz2,z1) - dotcc(z1,z2)*dotrr(pz1,pz2) )
     &    + ( dotrc(pz2,a)*dotrc(pa,z2) - dotcc(a,z2)*dotrr(pa,pz2) ) 
     &            * ( dotrc(pz1,z3)*dotrc(pz3,z1) - dotcc(z1,z3)*dotrr(pz1,pz3) )  
     &    + ( dotrc(pz1,a)*dotrc(pa,z1) - dotcc(a,z1)*dotrr(pa,pz1) ) 
     &            * ( dotrc(pz2,z3)*dotrc(pz3,z2) - dotcc(z2,z3)*dotrr(pz2,pz3) ) )
      endif

      if (ft1 .ne. 0) then
        azzz = azzz + ft1 * 8. * gw**4 * sinw * cosw**3 *
     &    ( ( dotrc(pz3,a)*dotrc(pa,z3) - dotcc(a,z3)*dotrr(pa,pz3) ) 
     &            * ( dotrc(pz1,z2)*dotrc(pz2,z1) - dotcc(z1,z2)*dotrr(pz1,pz2) )
     &    + ( dotrc(pz2,a)*dotrc(pa,z2) - dotcc(a,z2)*dotrr(pa,pz2) ) 
     &            * ( dotrc(pz1,z3)*dotrc(pz3,z1) - dotcc(z1,z3)*dotrr(pz1,pz3) )  
     &    + ( dotrc(pz1,a)*dotrc(pa,z1) - dotcc(a,z1)*dotrr(pa,pz1) ) 
     &            * ( dotrc(pz2,z3)*dotrc(pz3,z2) - dotcc(z2,z3)*dotrr(pz2,pz3) ) )
      endif

      if (ft2 .ne. 0) then
        azzz = azzz + ft2 * 2. * gw**4 * sinw * cosw**3 *
     &    ( dotrc(pz2,a)*dotrc(pz3,z1)*dotrc(pz1,z2)*dotrc(pa,z3)
     &    + dotrc(pz1,a)*dotrc(pz2,z1)*dotrc(pz3,z2)*dotrc(pa,z3)
     &    + dotrc(pz3,a)*dotrc(pz2,z1)*dotrc(pa,z2)*dotrc(pz1,z3)
     &    + dotrc(pz2,a)*dotrc(pa,z1)*dotrc(pz3,z2)*dotrc(pz1,z3)
     &    + dotrc(pz1,a)*dotrc(pz3,z1)*dotrc(pa,z2)*dotrc(pz2,z3)
     &    + dotrc(pz3,a)*dotrc(pa,z1)*dotrc(pz1,z2)*dotrc(pz2,z3)
     &    + dotcc(a,z1)*dotcc(z2,z3) * ( dotrr(pa,pz3)*dotrr(pz1,pz2) + dotrr(pa,pz2)*dotrr(pz1,pz3) )
     &    + dotcc(z2,z3) * ( - dotrc(pz1,a) * ( dotrc(pz3,z1)*dotrr(pa,pz2) + dotrc(pz2,z1)*dotrr(pa,pz3) )
     &                       + dotrc(pz3,a) * ( dotrc(pz2,z1)*dotrr(pa,pz1) - dotrc(pa,z1)*dotrr(pz1,pz2) )
     &                       + dotrc(pz2,a) * ( dotrc(pz3,z1)*dotrr(pa,pz1) - dotrc(pa,z1)*dotrr(pz1,pz3) ) )
     &    - dotcc(a,z1) * ( dotrc(pz3,z2) * ( dotrc(pz1,z3)*dotrr(pa,pz2) + dotrc(pa,z3)*dotrr(pz1,pz2) )
     &                     + dotrc(pz2,z3) * ( dotrc(pz1,z2)*dotrr(pa,pz3) + dotrc(pa,z2)*dotrr(pz1,pz3) )
     &                     - dotrr(pz2,pz3) * ( dotrc(pz1,z2)*dotrc(pa,z3) + dotrc(pa,z2)*dotrc(pz1,z3) ) )
     &    + dotcc(a,z2)*dotcc(z1,z3) * ( dotrr(pa,pz3)*dotrr(pz1,pz2) + dotrr(pa,pz1)*dotrr(pz2,pz3) )
     &    + dotcc(a,z3)*dotcc(z1,z2) * ( dotrr(pa,pz2)*dotrr(pz1,pz3) + dotrr(pa,pz1)*dotrr(pz2,pz3) )
     &    - dotcc(z1,z3) * ( + dotrc(pz2,a) * ( dotrc(pz3,z2)*dotrr(pa,pz1) + dotrc(pz1,z2)*dotrr(pa,pz3) )
     &                       + dotrc(pz3,a) * ( dotrc(pa,z2)*dotrr(pz1,pz2) - dotrc(pz1,z2)*dotrr(pa,pz2) )
     &                       + dotrc(pz1,a) * ( dotrc(pa,z2)*dotrr(pz2,pz3) - dotrc(pz3,z2)*dotrr(pa,pz2) ) )
     &    - dotcc(a,z3) * ( + dotrc(pz3,z1) * ( dotrc(pz1,z2)*dotrr(pa,pz2) - dotrc(pa,z2)*dotrr(pz1,pz2) )
     &                      + dotrc(pz2,z1) * ( dotrc(pz3,z2)*dotrr(pa,pz1) + dotrc(pa,z2)*dotrr(pz1,pz3) )
     &                      + dotrc(pa,z1) * ( dotrc(pz1,z2)*dotrr(pz2,pz3) - dotrc(pz3,z2)*dotrr(pz1,pz2) ) )
     &    - dotcc(z1,z2) * ( + dotrc(pz3,a) * ( dotrc(pz2,z3)*dotrr(pa,pz1) + dotrc(pz1,z3)*dotrr(pa,pz2) )
     &                       + dotrc(pz2,a) * ( dotrc(pa,z3)*dotrr(pz1,pz3) - dotrc(pz1,z3)*dotrr(pa,pz3) )
     &                       + dotrc(pz1,a) * ( dotrc(pa,z3)*dotrr(pz2,pz3) - dotrc(pz2,z3)*dotrr(pa,pz3) ) )
     &    - dotcc(a,z2) * ( + dotrc(pz3,z1) * ( dotrc(pz2,z3)*dotrr(pa,pz1) + dotrc(pa,z3)*dotrr(pz1,pz2) )
     &                      + dotrc(pz2,z1) * ( dotrc(pz1,z3)*dotrr(pa,pz3) - dotrc(pa,z3)*dotrr(pz1,pz3) )
     &                      + dotrc(pa,z1) * ( dotrc(pz1,z3)*dotrr(pz2,pz3) - dotrc(pz2,z3)*dotrr(pz1,pz3) ) ) )
      endif

      if (ft5 .ne. 0) then
        azzz = azzz - ft5 * 2. * gw**4 * ( sinw**3 * cosw - sinw**5 / cosw ) *
     &   ( ( dotrc(pz3,a)*dotrc(pa,z3) - dotcc(a,z3)*dotrr(pa,pz3) ) 
     &            * ( dotrc(pz1,z2)*dotrc(pz2,z1) - dotcc(z1,z2)*dotrr(pz1,pz2) )
     &   + ( dotrc(pz2,a)*dotrc(pa,z2) - dotcc(a,z2)*dotrr(pa,pz2) ) 
     &            * ( dotrc(pz1,z3)*dotrc(pz3,z1) - dotcc(z1,z3)*dotrr(pz1,pz3) )
     &   + ( dotrc(pz1,a)*dotrc(pa,z1) - dotcc(a,z1)*dotrr(pa,pz1) ) 
     &            * ( dotrc(pz2,z3)*dotrc(pz3,z2) - dotcc(z2,z3)*dotrr(pz2,pz3) ) )
      endif

      if (ft6 .ne. 0) then
        azzz = azzz - ft6 * 2. * gw**4 * ( sinw**3 * cosw - sinw**5 / cosw ) *
     &   ( ( dotrc(pz3,a)*dotrc(pa,z3) - dotcc(a,z3)*dotrr(pa,pz3) ) 
     &            * ( dotrc(pz1,z2)*dotrc(pz2,z1) - dotcc(z1,z2)*dotrr(pz1,pz2) )
     &   + ( dotrc(pz2,a)*dotrc(pa,z2) - dotcc(a,z2)*dotrr(pa,pz2) ) 
     &            * ( dotrc(pz1,z3)*dotrc(pz3,z1) - dotcc(z1,z3)*dotrr(pz1,pz3) )
     &   + ( dotrc(pz1,a)*dotrc(pa,z1) - dotcc(a,z1)*dotrr(pa,pz1) ) 
     &            * ( dotrc(pz2,z3)*dotrc(pz3,z2) - dotcc(z2,z3)*dotrr(pz2,pz3) ) )
      endif

      if (ft7 .ne. 0) then
        azzz = azzz - ft7 / 2. * gw**4 * ( sinw**3 * cosw - sinw**5 / cosw ) *
     &    ( dotrc(pz3,a) * dotrc(pz2,z1) * dotrc(pa,z2) * dotrc(pz1,z3)
     &    - dotcc(z1,z2) * dotrc(pz3,a) * dotrc(pz1,z3) * dotrr(pa,pz2)
     &    - dotcc(a,z1) * dotrc(pz3,z2) * dotrc(pz1,z3) * dotrr(pa,pz2)
     &    - dotcc(a,z2) * dotrc(pz2,z1) * dotrc(pz1,z3) * dotrr(pa,pz3)
     &    + dotrc(pz3,a) * dotrc(pa,z1) * dotrc(pz1,z2) * dotrc(pz2,z3)
     &    + dotcc(z2,z3) * dotrc(pz3,a) * dotrc(pz2,z1) * dotrr(pa,pz1)
     &    - dotcc(a,z3) * dotrc(pz2,z1) * dotrc(pz3,z2) * dotrr(pa,pz1)
     &    - dotcc(z1,z2) * dotrc(pz3,a) * dotrc(pz2,z3) * dotrr(pa,pz1)
     &    - dotcc(a,z2) * dotrc(pz3,z1) * dotrc(pz2,z3) * dotrr(pa,pz1)
     &    + dotcc(z1,z3) * dotrc(pz3,a) * dotrc(pz1,z2) * dotrr(pa,pz2)
     &    - dotcc(a,z3) * dotrc(pz3,z1) * dotrc(pz1,z2) * dotrr(pa,pz2)
     &    - dotcc(a,z1) * dotrc(pz1,z2) * dotrc(pz2,z3) * dotrr(pa,pz3)
     &    - dotcc(z2,z3) * dotrc(pz3,a) * dotrc(pa,z1) * dotrr(pz1,pz2)
     &    - dotcc(z1,z3) * dotrc(pz3,a) * dotrc(pa,z2) * dotrr(pz1,pz2)
     &    + dotcc(a,z3) * dotrc(pz3,z1) * dotrc(pa,z2) * dotrr(pz1,pz2)
     &    + dotcc(a,z3) * dotrc(pa,z1) * dotrc(pz3,z2) * dotrr(pz1,pz2)
     &    - dotcc(a,z2) * dotrc(pz3,z1) * dotrc(pa,z3) * dotrr(pz1,pz2)
     &    - dotcc(a,z1) * dotrc(pz3,z2) * dotrc(pa,z3) * dotrr(pz1,pz2)
     &    + dotcc(a,z2) * dotcc(z1,z3) * dotrr(pa,pz3) * dotrr(pz1,pz2)
     &    + dotcc(a,z1) * dotcc(z2,z3) * dotrr(pa,pz3) * dotrr(pz1,pz2)
     &    - dotcc(a,z3) * dotrc(pz2,z1) * dotrc(pa,z2) * dotrr(pz1,pz3)
     &    + dotcc(a,z2) * dotrc(pz2,z1) * dotrc(pa,z3) * dotrr(pz1,pz3)
     &    + dotcc(a,z2) * dotrc(pa,z1) * dotrc(pz2,z3) * dotrr(pz1,pz3)
     &    - dotcc(a,z1) * dotrc(pa,z2) * dotrc(pz2,z3) * dotrr(pz1,pz3)
     &    + dotcc(a,z3) * dotcc(z1,z2) * dotrr(pa,pz2) * dotrr(pz1,pz3)
     &    + dotcc(a,z1) * dotcc(z2,z3) * dotrr(pa,pz2) * dotrr(pz1,pz3)
     &    + dotrc(pz2,a) * ( dotrc(pa,z1)*dotrc(pz3,z2)*dotrc(pz1,z3)
     &                     + dotcc(z1,z2)*dotrc(pz1,z3)*dotrr(pa,pz3)
     &                     + dotrc(pz3,z1) * ( dotrc(pz1,z2)*dotrc(pa,z3) + dotcc(z2,z3)*dotrr(pa,pz1) )
     &                     - dotcc(z1,z3) * ( dotrc(pz3,z2)*dotrr(pa,pz1) + dotrc(pz1,z2)*dotrr(pa,pz3) )
     &                     - dotrr(pz1,pz3) * ( dotrc(pa,z1)*dotcc(z2,z3) + dotrc(pa,z3)*dotcc(z1,z2) ) )
     &    + dotrr(pz2,pz3) * ( dotcc(a,z1) * ( dotrc(pz1,z2)*dotrc(pa,z3) + dotrc(pa,z2)*dotrc(pz1,z3) )
     &                       + dotcc(a,z3) * ( dotcc(z1,z2)*dotrr(pa,pz1) - dotrc(pa,z1)*dotrc(pz1,z2) )
     &                       + dotcc(a,z2) * ( dotcc(z1,z3)*dotrr(pa,pz1) - dotrc(pa,z1)*dotrc(pz1,z3) ) )
     &    + dotrc(pz1,a) * ( dotrc(pz3,z1) * ( dotrc(pa,z2)*dotrc(pz2,z3) - dotcc(z2,z3)*dotrr(pa,pz2) )
     &                     + dotrc(pz2,z1) * ( dotrc(pz3,z2)*dotrc(pa,z3) - dotcc(z2,z3)*dotrr(pa,pz3) )
     &                     + dotcc(z1,z3) * ( dotrc(pz3,z2)*dotrr(pa,pz2) - dotrc(pa,z2)*dotrr(pz2,pz3) )
     &                     + dotcc(z1,z2) * ( dotrc(pz2,z3)*dotrr(pa,pz3) - dotrc(pa,z3)*dotrr(pz2,pz3) ) ) )
      endif

      if (ft8 .ne. 0) then
        azzz = azzz - ft8 * 2. * gw**4 * sinw**7 / cosw**3 *
     &   ( ( dotrc(pz3,a)*dotrc(pa,z3) - dotcc(a,z3)*dotrr(pa,pz3) ) 
     &            * ( dotrc(pz1,z2)*dotrc(pz2,z1) - dotcc(z1,z2)*dotrr(pz1,pz2) )
     &   + ( dotrc(pz2,a)*dotrc(pa,z2) - dotcc(a,z2)*dotrr(pa,pz2) ) 
     &            * ( dotrc(pz1,z3)*dotrc(pz3,z1) - dotcc(z1,z3)*dotrr(pz1,pz3) )
     &   + ( dotrc(pz1,a)*dotrc(pa,z1) - dotcc(a,z1)*dotrr(pa,pz1) ) 
     &            * ( dotrc(pz2,z3)*dotrc(pz3,z2) - dotcc(z2,z3)*dotrr(pz2,pz3) ) )
      endif

      if (ft9 .ne. 0) then
        azzz = azzz - ft9 / 2. * gw**4 * sinw**7 / cosw**3 *
     &    ( dotrc(pz3,a) * dotrc(pz2,z1) * dotrc(pa,z2) * dotrc(pz1,z3)
     &    - dotcc(z1,z2) * dotrc(pz3,a) * dotrc(pz1,z3) * dotrr(pa,pz2)
     &    - dotcc(a,z1) * dotrc(pz3,z2) * dotrc(pz1,z3) * dotrr(pa,pz2)
     &    - dotcc(a,z2) * dotrc(pz2,z1) * dotrc(pz1,z3) * dotrr(pa,pz3)
     &    + dotrc(pz3,a) * dotrc(pa,z1) * dotrc(pz1,z2) * dotrc(pz2,z3)
     &    + dotcc(z2,z3) * dotrc(pz3,a) * dotrc(pz2,z1) * dotrr(pa,pz1)
     &    - dotcc(a,z3) * dotrc(pz2,z1) * dotrc(pz3,z2) * dotrr(pa,pz1)
     &    - dotcc(z1,z2) * dotrc(pz3,a) * dotrc(pz2,z3) * dotrr(pa,pz1)
     &    - dotcc(a,z2) * dotrc(pz3,z1) * dotrc(pz2,z3) * dotrr(pa,pz1)
     &    + dotcc(z1,z3) * dotrc(pz3,a) * dotrc(pz1,z2) * dotrr(pa,pz2)
     &    - dotcc(a,z3) * dotrc(pz3,z1) * dotrc(pz1,z2) * dotrr(pa,pz2)
     &    - dotcc(a,z1) * dotrc(pz1,z2) * dotrc(pz2,z3) * dotrr(pa,pz3)
     &    - dotcc(z2,z3) * dotrc(pz3,a) * dotrc(pa,z1) * dotrr(pz1,pz2)
     &    - dotcc(z1,z3) * dotrc(pz3,a) * dotrc(pa,z2) * dotrr(pz1,pz2)
     &    + dotcc(a,z3) * dotrc(pz3,z1) * dotrc(pa,z2) * dotrr(pz1,pz2)
     &    + dotcc(a,z3) * dotrc(pa,z1) * dotrc(pz3,z2) * dotrr(pz1,pz2)
     &    - dotcc(a,z2) * dotrc(pz3,z1) * dotrc(pa,z3) * dotrr(pz1,pz2)
     &    - dotcc(a,z1) * dotrc(pz3,z2) * dotrc(pa,z3) * dotrr(pz1,pz2)
     &    + dotcc(a,z2) * dotcc(z1,z3) * dotrr(pa,pz3) * dotrr(pz1,pz2)
     &    + dotcc(a,z1) * dotcc(z2,z3) * dotrr(pa,pz3) * dotrr(pz1,pz2)
     &    - dotcc(a,z3) * dotrc(pz2,z1) * dotrc(pa,z2) * dotrr(pz1,pz3)
     &    + dotcc(a,z2) * dotrc(pz2,z1) * dotrc(pa,z3) * dotrr(pz1,pz3)
     &    + dotcc(a,z2) * dotrc(pa,z1) * dotrc(pz2,z3) * dotrr(pz1,pz3)
     &    - dotcc(a,z1) * dotrc(pa,z2) * dotrc(pz2,z3) * dotrr(pz1,pz3)
     &    + dotcc(a,z3) * dotcc(z1,z2) * dotrr(pa,pz2) * dotrr(pz1,pz3)
     &    + dotcc(a,z1) * dotcc(z2,z3) * dotrr(pa,pz2) * dotrr(pz1,pz3)
     &    + dotrc(pz2,a) * ( dotrc(pa,z1)*dotrc(pz3,z2)*dotrc(pz1,z3)
     &                     + dotcc(z1,z2)*dotrc(pz1,z3)*dotrr(pa,pz3)
     &                     + dotrc(pz3,z1) * ( dotrc(pz1,z2)*dotrc(pa,z3) + dotcc(z2,z3)*dotrr(pa,pz1) )
     &                     - dotcc(z1,z3) * ( dotrc(pz3,z2)*dotrr(pa,pz1) + dotrc(pz1,z2)*dotrr(pa,pz3) )
     &                     - dotrr(pz1,pz3) * ( dotrc(pa,z1)*dotcc(z2,z3) + dotrc(pa,z3)*dotcc(z1,z2) ) )
     &    + dotrr(pz2,pz3) * ( dotcc(a,z1) * ( dotrc(pz1,z2)*dotrc(pa,z3) + dotrc(pa,z2)*dotrc(pz1,z3) )
     &                       + dotcc(a,z3) * ( dotcc(z1,z2)*dotrr(pa,pz1) - dotrc(pa,z1)*dotrc(pz1,z2) )
     &                       + dotcc(a,z2) * ( dotcc(z1,z3)*dotrr(pa,pz1) - dotrc(pa,z1)*dotrc(pz1,z3) ) )
     &    + dotrc(pz1,a) * ( dotrc(pz3,z1) * ( dotrc(pa,z2)*dotrc(pz2,z3) - dotcc(z2,z3)*dotrr(pa,pz2) )
     &                     + dotrc(pz2,z1) * ( dotrc(pz3,z2)*dotrc(pa,z3) - dotcc(z2,z3)*dotrr(pa,pz3) )
     &                     + dotcc(z1,z3) * ( dotrc(pz3,z2)*dotrr(pa,pz2) - dotrc(pa,z2)*dotrr(pz2,pz3) )
     &                     + dotcc(z1,z2) * ( dotrc(pz2,z3)*dotrr(pa,pz3) - dotrc(pa,z3)*dotrr(pz2,pz3) ) ) )
      endif

      return
      end

*********************************************************************************************************
CCCCCCCCCCCCCCC
C AAZZ-vertex C
CCCCCCCCCCCCCCC

      subroutine aazz_anomal4(a1,a2,z1,z2, aazz)
c
c This subroutine computes the aazz-vertex with anomalous couplings
c where all the bosons are defined by the flowing-OUT quantum number.
c     
      implicit none
      double complex z1(6),z2(6),a1(6),a2(6),aazz 
      double complex dotcc, dotrc, epsccrr, epscccr, epscccc
      double precision pz1(0:3),pz2(0:3),pa1(0:3),pa2(0:3), dotrr
      double precision cos2w, cosw, sinw

      external dotcc, dotrc, dotrr, epsccrr, epscccr, epscccc

#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"   

      cos2w = 1.0d0 - sin2w
      cosw = sqrt(cos2w)
      sinw = sqrt(sin2w)


      pz1(0) = dble( z1(5))
      pz1(1) = dble( z1(6))
      pz1(2) = dimag(z1(6))
      pz1(3) = dimag(z1(5))
      pz2(0) = dble( z2(5))
      pz2(1) = dble( z2(6))
      pz2(2) = dimag(z2(6))
      pz2(3) = dimag(z2(5))
      pa1(0) = dble( a1(5))
      pa1(1) = dble( a1(6))
      pa1(2) = dimag(a1(6))
      pa1(3) = dimag(a1(5))
      pa2(0) = dble( a2(5))
      pa2(1) = dble( a2(6))
      pa2(2) = dimag(a2(6))
      pa2(3) = dimag(a2(5))


      aazz = dcmplx(0.,0.)

      if (fm0 .ne. 0) then
        aazz = aazz - 2. * fm0 * gw**2 * wmass**2 * sin2w / cos2w *
     &    dotcc(z1,z2) * ( dotrc(pa2,a1)*dotrc(pa1,a2) - dotcc(a1,a2)*dotrr(pa1,pa2) )
      endif

      if (fm1 .ne. 0) then
        aazz = aazz + fm1 / 2. * gw**2 * wmass**2 * sin2w / cos2w *
     &    ( dotrc(pa2,a1) * ( dotrc(pa1,z1)*dotcc(a2,z2) + dotrc(pa1,z2)*dotcc(a2,z1) )
     &    + dotrc(pa1,a2) * ( dotrc(pa2,z1)*dotcc(a1,z2) + dotrc(pa2,z2)*dotcc(a1,z1) )
     &    - dotcc(a1,a2) * ( dotrc(pa1,z2)*dotrc(pa2,z1) + dotrc(pa1,z1)*dotrc(pa2,z2) )
     &    - dotrr(pa1,pa2) * ( dotcc(a1,z2)*dotcc(a2,z1) + dotcc(a1,z1)*dotcc(a2,z2) ) )
      endif

      if (fm2 .ne. 0) then
        aazz = aazz + fm2 * gw**2 * wmass**2 * sin2w / cos2w *
     &    dotcc(z1,z2) * ( dotcc(a1,a2)*dotrr(pa1,pa2) - dotrc(pa2,a1)*dotrc(pa1,a2) )
      endif

      if (fm3 .ne. 0) then
        aazz = aazz + fm3 / 4. * gw**2 * wmass**2 * sin2w / cos2w *
     &    ( dotrc(pa2,a1) * ( dotrc(pa1,z1)*dotcc(a2,z2) + dotrc(pa1,z2)*dotcc(a2,z1) )
     &    + dotrc(pa1,a2) * ( dotrc(pa2,z1)*dotcc(a1,z2) + dotrc(pa2,z2)*dotcc(a1,z1) )
     &    - dotcc(a1,a2) * ( dotrc(pa1,z2)*dotrc(pa2,z1) + dotrc(pa1,z1)*dotrc(pa2,z2) )
     &    - dotrr(pa1,pa2) * ( dotcc(a1,z2)*dotcc(a2,z1) + dotcc(a1,z1)*dotcc(a2,z2) ) )
      endif

      if (fm4 .ne. 0) then
        aazz = aazz + fm4 * gw**2 * wmass**2 * sin2w / cos2w *
     &    dotcc(z1,z2) * (  dotrc(pa2,a1)*dotrc(pa1,a2) - dotcc(a1,a2)*dotrr(pa1,pa2) )
      endif

      if (fm5 .ne. 0) then
        aazz = aazz + fm5 / 4. * gw**2 * wmass**2 * sin2w / cos2w *
     &    ( dotrc(pa2,a1) * ( dotrc(pa1,z1)*dotcc(a2,z2) + dotrc(pa1,z2)*dotcc(a2,z1) )
     &    + dotrc(pa1,a2) * ( dotrc(pa2,z1)*dotcc(a1,z2) + dotrc(pa2,z2)*dotcc(a1,z1) )
     &    - dotcc(a1,a2) * ( dotrc(pa1,z2)*dotrc(pa2,z1) + dotrc(pa1,z1)*dotrc(pa2,z2) )
     &    - dotrr(pa1,pa2) * ( dotcc(a1,z2)*dotcc(a2,z1) + dotcc(a1,z1)*dotcc(a2,z2) ) )
      endif

      if (fm6 .ne. 0) then
        aazz = aazz - fm6 * gw**2 * wmass**2 * sin2w / cos2w *
     &    dotcc(z1,z2) * ( dotrc(pa2,a1)*dotrc(pa1,a2) - dotcc(a1,a2)*dotrr(pa1,pa2) )
      endif

      if (fm7 .ne. 0) then
        aazz = aazz - fm7 / 4. * gw**2 * wmass**2 * sin2w / cos2w *
     &    ( dotrc(pa2,a1) * ( dotrc(pa1,z1)*dotcc(a2,z2) + dotrc(pa1,z2)*dotcc(a2,z1) )
     &    + dotrc(pa1,a2) * ( dotrc(pa2,z1)*dotcc(a1,z2) + dotrc(pa2,z2)*dotcc(a1,z1) )
     &    - dotcc(a1,a2) * ( dotrc(pa1,z2)*dotrc(pa2,z1) + dotrc(pa1,z1)*dotrc(pa2,z2) )
     &    - dotrr(pa1,pa2) * ( dotcc(a1,z2)*dotcc(a2,z1) + dotcc(a1,z1)*dotcc(a2,z2) ) )
      endif

      if (ft0 .ne. 0) then
        aazz = aazz + ft0 * 8. * gw**4 * sin2w * cos2w *
     &    ( ( dotrc(pz2,a1)*dotrc(pa1,z2) - dotcc(a1,z2)*dotrr(pa1,pz2) ) 
     &            * ( dotrc(pz1,a2)*dotrc(pa2,z1) - dotcc(z1,a2)*dotrr(pz1,pa2) )
     &    + ( dotrc(pz1,a1)*dotrc(pa1,z1) - dotcc(a1,z1)*dotrr(pa1,pz1) ) 
     &            * ( dotrc(pz2,a2)*dotrc(pa2,z2) - dotcc(a2,z2)*dotrr(pa2,pz2) )  
     &    + ( dotrc(pa2,a1)*dotrc(pa1,a2) - dotcc(a1,a2)*dotrr(pa1,pa2) ) 
     &            * ( dotrc(pz1,z2)*dotrc(pz2,z1) - dotcc(z1,z2)*dotrr(pz1,pz2) ) )
      endif

      if (ft1 .ne. 0) then
        aazz = aazz + ft1 * 8. * gw**4 * sin2w * cos2w *
     &    ( ( dotrc(pz2,a1)*dotrc(pa1,z2) - dotcc(a1,z2)*dotrr(pa1,pz2) ) 
     &            * ( dotrc(pz1,a2)*dotrc(pa2,z1) - dotcc(z1,a2)*dotrr(pz1,pa2) )
     &    + ( dotrc(pz1,a1)*dotrc(pa1,z1) - dotcc(a1,z1)*dotrr(pa1,pz1) ) 
     &            * ( dotrc(pz2,a2)*dotrc(pa2,z2) - dotcc(a2,z2)*dotrr(pa2,pz2) )


     &    + ( dotrc(pa2,a1)*dotrc(pa1,a2) - dotcc(a1,a2)*dotrr(pa1,pa2) ) 
     &            * ( dotrc(pz1,z2)*dotrc(pz2,z1) - dotcc(z1,z2)*dotrr(pz1,pz2) ) )
      endif

      if (ft2 .ne. 0) then
        aazz = aazz + ft2 * 2. * gw**4 * sin2w * cos2w *
     &    ( dotrc(pz1,a1)*dotrc(pz2,a2)*dotrc(pa2,z1)*dotrc(pa1,z2)
     &    + dotrc(pa2,a1)*dotrc(pz1,a2)*dotrc(pz2,z1)*dotrc(pa1,z2)
     &    + dotrc(pz2,a1)*dotrc(pz1,a2)*dotrc(pa1,z1)*dotrc(pa2,z2)
     &    + dotrc(pz1,a1)*dotrc(pa1,a2)*dotrc(pz2,z1)*dotrc(pa2,z2)
     &    + dotrc(pa2,a1)*dotrc(pz2,a2)*dotrc(pa1,z1)*dotrc(pz1,z2)
     &    + dotrc(pz2,a1)*dotrc(pa1,a2)*dotrc(pa2,z1)*dotrc(pz1,z2)
     &    + dotcc(a1,a2)*dotcc(z1,z2) * ( dotrr(pa1,pz2)*dotrr(pa2,pz1) + dotrr(pa1,pz1)*dotrr(pa2,pz2) )
     &    + dotcc(z1,z2) * ( - dotrc(pa2,a1) * ( dotrc(pz2,a2)*dotrr(pa1,pz1) + dotrc(pz1,a2)*dotrr(pa1,pz2) )
     &                       + dotrc(pz2,a1) * ( dotrc(pz1,a2)*dotrr(pa1,pa2) - dotrc(pa1,a2)*dotrr(pa2,pz1) )
     &                       + dotrc(pz1,a1) * ( dotrc(pz2,a2)*dotrr(pa1,pa2) - dotrc(pa1,a2)*dotrr(pa2,pz2) ) )
     &    - dotcc(a1,a2) * ( dotrc(pz2,z1) * ( dotrc(pa2,z2)*dotrr(pa1,pz1) + dotrc(pa1,z2)*dotrr(pa2,pz1) )
     &                     + dotrc(pz1,z2) * ( dotrc(pa2,z1)*dotrr(pa1,pz2) + dotrc(pa1,z1)*dotrr(pa2,pz2) )
     &                     - dotrr(pz1,pz2) * ( dotrc(pa2,z1)*dotrc(pa1,z2) + dotrc(pa1,z1)*dotrc(pa2,z2) ) )
     &    + dotcc(a1,z1)*dotcc(a2,z2) * ( dotrr(pa1,pz2)*dotrr(pa2,pz1) + dotrr(pa1,pa2)*dotrr(pz1,pz2) )
     &    + dotcc(a1,z2)*dotcc(a2,z1) * ( dotrr(pa1,pz1)*dotrr(pa2,pz2) + dotrr(pa1,pa2)*dotrr(pz1,pz2) )
     &    - dotcc(a2,z2) * ( + dotrc(pz1,a1) * ( dotrc(pz2,z1)*dotrr(pa1,pa2) + dotrc(pa2,z1)*dotrr(pa1,pz2) )
     &                       + dotrc(pz2,a1) * ( dotrc(pa1,z1)*dotrr(pa2,pz1) - dotrc(pa2,z1)*dotrr(pa1,pz1) )
     &                       + dotrc(pa2,a1) * ( dotrc(pa1,z1)*dotrr(pz1,pz2) - dotrc(pz2,z1)*dotrr(pa1,pz1) ) )
     &    - dotcc(a1,z2) * ( + dotrc(pz2,a2) * ( dotrc(pa2,z1)*dotrr(pa1,pz1) - dotrc(pa1,z1)*dotrr(pa2,pz1) )
     &                       + dotrc(pz1,a2) * ( dotrc(pz2,z1)*dotrr(pa1,pa2) + dotrc(pa1,z1)*dotrr(pa2,pz2) )
     &                       + dotrc(pa1,a2) * ( dotrc(pa2,z1)*dotrr(pz1,pz2) - dotrc(pz2,z1)*dotrr(pa2,pz1) ) )
     &    - dotcc(a2,z1) * ( + dotrc(pz2,a1) * ( dotrc(pz1,z2)*dotrr(pa1,pa2) + dotrc(pa2,z2)*dotrr(pa1,pz1) )
     &                       + dotrc(pz1,a1) * ( dotrc(pa1,z2)*dotrr(pa2,pz2) - dotrc(pa2,z2)*dotrr(pa1,pz2) )
     &                       + dotrc(pa2,a1) * ( dotrc(pa1,z2)*dotrr(pz1,pz2) - dotrc(pz1,z2)*dotrr(pa1,pz2) ) )
     &    - dotcc(a1,z1) * ( + dotrc(pz2,a2) * ( dotrc(pz1,z2)*dotrr(pa1,pa2) + dotrc(pa1,z2)*dotrr(pa2,pz1) )
     &                       + dotrc(pz1,a2) * ( dotrc(pa2,z2)*dotrr(pa1,pz2) - dotrc(pa1,z2)*dotrr(pa2,pz2) )
     &                       + dotrc(pa1,a2) * ( dotrc(pa2,z2)*dotrr(pz1,pz2) - dotrc(pz1,z2)*dotrr(pa2,pz2) ) ) )
      endif

      if (ft5 .ne. 0) then
        aazz = aazz - ft5 / 2. * gw**4 *
     &   ( 4.*sin2w**2 * ( 2. * ( dotrc(pz2,a1)*dotrc(pa1,z2) - dotcc(a1,z2)*dotrr(pa1,pz2) ) 
     &                        * ( dotrc(pz1,a2)*dotrc(pa2,z1) - dotcc(z1,a2)*dotrr(pz1,pa2) )
     &                   + 2. * ( dotrc(pz1,a1)*dotrc(pa1,z1) - dotcc(a1,z1)*dotrr(pa1,pz1) ) 
     &                        * ( dotrc(pz2,a2)*dotrc(pa2,z2) - dotcc(z2,a2)*dotrr(pz2,pa2) )
     &                   - sin2w / cos2w * dotrc(pa2,a1)*dotrc(pa1,a2)*dotrc(pz1,z2)*dotrc(pz2,z1) )
     &   + ( 4.*sin2w/cos2w - 8.*sin2w**2 ) * dotcc(a1,a2)*dotrr(pa1,pa2) 
     &            * ( dotrc(pz1,z2)*dotrc(pz2,z1) - dotcc(z1,z2)*dotrr(pz1,pz2) )
     &   + ( 4.*sin2w/cos2w - 8.*sin2w**2 ) * dotrc(pa2,a1)*dotrc(pa1,a2)*dotcc(z1,z2)*dotrr(pz1,pz2)
     &   - 4.*sin2w*cos2w * dotrc(pa2,a1)*dotrc(pa1,a2)*dotrc(pz1,z2)*dotrc(pz2,z1) )
      endif

      if (ft6 .ne. 0) then
        aazz = aazz - ft6 / 4. * gw**4 *
     &    ( 16. * sin2w**2 * dotrc(pa2,a1)*dotrc(pa1,a2) * ( dotrc(pz1,z2)*dotrc(pz2,z1) - dotcc(z1,z2)*dotrr(pz1,pz2) )
     &    + 2. * sin2w/cos2w * ( 2. * (cos2w-sin2w)**2 * dotrc(pz2,a1)*dotrc(pa1,z2) 
     &                         * ( dotcc(a2,z1)*dotrr(pa2,pz1) - dotrc(pz1,a2)*dotrc(pa2,z1) )
     &                         + ( 2. - 8.*sin2w*cos2w ) 
     &                           * ( dotcc(a1,z2)*dotrr(pa1,pz2) * ( dotrc(pz1,a2)*dotrc(pa2,z1) - dotcc(a2,z1)*dotrr(pa2,pz1) )
     &                             - ( dotrc(pz1,a1)*dotrc(pa1,z1) - dotcc(a1,z1)*dotrr(pa1,pz1) ) 
     &                             * ( dotrc(pz2,a2)*dotrc(pa2,z2) - dotcc(a2,z2)*dotrr(pa2,pz2) ) )
     &                         - ( 1. - 8.*sin2w*cos2w ) * dotcc(a1,a2)*dotcc(z1,z2)*dotrr(pa1,pa2)*dotrr(pz1,pz2)
     &                         + dotcc(a1,a2)*dotrr(pa1,pa2) 
     &                           * ( dotcc(z1,z2)*dotrr(pz1,pz2) -  8.*sin2w*cos2w*dotrc(pz1,z2)*dotrc(pz2,z1) ) ) )
      endif

      if (ft7 .ne. 0) then
        aazz = aazz + ft7 / 8. * gw**4 * sin2w / cos2w *
     &    ( - 2.*(1.-4.*cos2w*sin2w) * dotrc(pz2,a1) * dotrc(pa1,a2) * dotcc(z1,z2) * dotrr(pa2,pz1)
     &    + 2.*(1.-4.*cos2w*sin2w) * dotrr(pz1,pz2) * ( dotcc(a1,a2) * ( dotrc(pa2,z1)*dotrc(pa1,z2) 
     &                                                                  + dotrc(pa1,z1)*dotrc(pa2,z2) )
     &                                                + dotcc(a1,z2) * ( dotcc(a2,z1)*dotrr(pa1,pa2) 
     &                                                                  - dotrc(pa1,a2)*dotrc(pa2,z1) )
     &                                                + dotcc(a1,z1) * ( dotcc(a2,z2)*dotrr(pa1,pa2) 
     &                                                                  - dotrc(pa1,a2)*dotrc(pa2,z2) ) )
     &    - 8.*sin2w*cos2w * dotrc(pz2,a1) * dotrc(pz1,a2) * dotrc(pa1,z1) * dotrc(pa2,z2)
     &    + (2.-8.*sin2w*cos2w) * dotrc(pz2,a1) * dotrc(pa1,a2) * dotrc(pa2,z1) * dotrc(pz1,z2)
     &    + (2.-8.*sin2w*cos2w) * dotrc(pz2,a1) * dotrc(pz1,a2) * dotcc(z1,z2) * dotrr(pa1,pa2)
     &    - (2.-8.*sin2w*cos2w) * dotcc(a1,z2) * dotrc(pz1,a2) * dotrc(pz2,z1) * dotrr(pa1,pa2)
     &    - (2.-8.*sin2w*cos2w) * dotcc(a2,z1) * dotrc(pz1,z2) * dotrc(pz2,a1) * dotrr(pa1,pa2)
     &    - (2.-8.*sin2w*cos2w) * dotcc(a1,z1) * dotrc(pz1,z2) * dotrc(pz2,a2) * dotrr(pa1,pa2)
     &    - 8.*sin2w*cos2w * dotcc(a2,z2) * dotrc(pz2,a1) * dotrc(pa2,z1) * dotrr(pa1,pz1)
     &    + 8.*sin2w*cos2w * dotcc(a1,z2) * dotrc(pz2,a2) * dotrc(pa2,z1) * dotrr(pa1,pz1)
     &    + 8.*sin2w*cos2w * dotcc(a2,z1) * dotrc(pz2,a1) * dotrc(pa2,z2) * dotrr(pa1,pz1)
     &    - (2.-8.*sin2w*cos2w) * dotcc(a1,a2) * dotrc(pz2,z1) * dotrc(pa2,z2) * dotrr(pa1,pz1)
     &    + 8.*sin2w*cos2w * dotcc(a1,z1) * dotrc(pz1,a2) * dotrc(pa2,z2) * dotrr(pa1,pz2)
     &    - (2.-8.*sin2w*cos2w) * dotcc(a1,a2) * dotrc(pa2,z1) * dotrc(pz1,z2) * dotrr(pa1,pz2)
     &    + 8.*sin2w*cos2w * dotcc(a2,z2) * dotrc(pz2,a1) * dotrc(pa1,z1) * dotrr(pa2,pz1)
     &    - 8.*sin2w*cos2w * dotcc(a1,z2) * dotrc(pz2,a2) * dotrc(pa1,z1) * dotrr(pa2,pz1)
     &    + (2.-8.*sin2w*cos2w) * dotrc(pa1,a2) * dotrc(pz2,z1) * dotcc(a1,z2) * dotrr(pa2,pz1)
     &    + 8.*sin2w*cos2w * dotcc(a1,z1) * dotrc(pz2,a2) * dotrc(pa1,z2) * dotrr(pa2,pz1)
     &    - (2.-8.*sin2w*cos2w) * dotcc(a1,a2) * dotrc(pz2,z1) * dotrc(pa1,z2) * dotrr(pa2,pz1)
     &    - 8.*sin2w*cos2w * dotcc(a1,z1) * dotcc(a2,z2) * dotrr(pa1,pz2) * dotrr(pa2,pz1)
     &    + (2.-8.*sin2w*cos2w) * dotcc(a1,a2) * dotcc(z1,z2) * dotrr(pa1,pz2) * dotrr(pa2,pz1)
     &    + 8.*sin2w*cos2w * dotcc(a1,z2) * dotrc(pz1,a2) * dotrc(pa1,z1) * dotrr(pa2,pz2)
     &    - 8.*sin2w*cos2w * dotcc(a1,z1) * dotrc(pz1,a2) * dotrc(pa1,z2) * dotrr(pa2,pz2)
     &    + (2.-8.*sin2w*cos2w) * dotrc(pa1,a2) * dotrc(pz1,z2) * dotcc(a1,z1) * dotrr(pa2,pz2)
     &    - (2.-8.*sin2w*cos2w) * dotcc(a1,a2) * dotrc(pz1,z2) * dotrc(pa1,z1) * dotrr(pa2,pz2)
     &    - 8.*sin2w*cos2w * dotcc(a1,z2) * dotcc(a2,z1) * dotrr(pa1,pz1) * dotrr(pa2,pz2)
     &    + (2.-8.*sin2w*cos2w) * dotcc(a1,a2) * dotcc(z1,z2) * dotrr(pa1,pz1) * dotrr(pa2,pz2)
     &    + (2.-8.*sin2w*cos2w) * dotrc(pa2,a1) * ( dotrc(pz2,a2) * ( dotrc(pa1,z1)*dotrc(pz1,z2) - dotcc(z1,z2)*dotrr(pa1,pz1) )
     &                                            + dotrc(pz1,a2) * ( dotrc(pz2,z1)*dotrc(pa1,z2) - dotcc(z1,z2)*dotrr(pa1,pz2) )
     &                                            + dotcc(a2,z2) * ( dotrc(pz2,z1)*dotrr(pa1,pz1) - dotrc(pa1,z1)*dotrr(pz1,pz2) )
     &                                            + dotcc(a2,z1) * ( dotrc(pz1,z2)*dotrr(pa1,pz2) - dotrc(pa1,z2)*dotrr(pz1,pz2)))
     &    + dotrc(pz1,a1) * ( 8.*sin2w*cos2w * dotcc(a2,z1)*dotrc(pa1,z2)*dotrr(pa2,pz2)
     &                      - (2.-8.*sin2w*cos2w) * dotcc(a2,z2) * dotrc(pz2,z1) * dotrr(pa1,pa2)   
     &                      + dotrc(pz2,a2) * ( - 8.*sin2w*cos2w * dotrc(pa2,z1)*dotrc(pa1,z2) 
     &                                          + (2.-8.*sin2w*cos2w) * dotcc(z1,z2)*dotrr(pa1,pa2) )
     &                      + 8.*sin2w*cos2w * dotcc(a2,z2) * dotrc(pa2,z1) * dotrr(pa1,pz2)
     &                      - 8.*sin2w*cos2w * dotcc(a2,z1) * dotrc(pa2,z2) * dotrr(pa1,pz2)
     &                      + (2.-8.*sin2w*cos2w) * dotrc(pa1,a2) 
     &                              * ( dotrc(pz2,z1)*dotrc(pa2,z2) - dotcc(z1,z2)*dotrr(pa2,pz2) ) ) )
      endif

      if (ft8 .ne. 0) then
        aazz = aazz + ft8 * 2. * gw**4 * sin2w**3 / cos2w *
     &    ( ( dotrc(pz2,a1)*dotrc(pa1,z2) - dotcc(a1,z2)*dotrr(pa1,pz2) ) 
     &            * ( dotrc(pz1,a2)*dotrc(pa2,z1) - dotcc(z1,a2)*dotrr(pz1,pa2) )
     &    + ( dotrc(pz1,a1)*dotrc(pa1,z1) - dotcc(a1,z1)*dotrr(pa1,pz1) ) 
     &            * ( dotrc(pz2,a2)*dotrc(pa2,z2) - dotcc(a2,z2)*dotrr(pa2,pz2) )  
     &    + ( dotrc(pa2,a1)*dotrc(pa1,a2) - dotcc(a1,a2)*dotrr(pa1,pa2) ) 
     &            * ( dotrc(pz1,z2)*dotrc(pz2,z1) - dotcc(z1,z2)*dotrr(pz1,pz2) ) )
      endif

      if (ft9 .ne. 0) then
        aazz = aazz + ft9 / 2. * gw**4 * sin2w**3 / cos2w *
     &    ( dotrc(pz1,a1) * dotrc(pz2,a2) * dotrc(pa2,z1) * dotrc(pa1,z2)
     &    + dotrc(pa2,a1) * dotrc(pz1,a2) * dotrc(pz2,z1) * dotrc(pa1,z2)
     &    + dotrc(pz2,a1) * dotrc(pz1,a2) * dotrc(pa1,z1) * dotrc(pa2,z2)
     &    + dotrc(pz1,a1) * dotrc(pa1,a2) * dotrc(pz2,z1) * dotrc(pa2,z2)
     &    + dotrc(pa2,a1) * dotrc(pz2,a2) * dotrc(pa1,z1) * dotrc(pz1,z2)
     &    + dotrc(pz2,a1) * dotrc(pa1,a2) * dotrc(pa2,z1) * dotrc(pz1,z2)
     &    + dotcc(a1,a2)*dotcc(z1,z2) * ( dotrr(pa1,pz2)*dotrr(pa2,pz1) + dotrr(pa1,pz1)*dotrr(pa2,pz2) )
     &    + dotcc(z1,z2) * ( - dotrc(pa2,a1) * ( dotrc(pz2,a2)*dotrr(pa1,pz1) + dotrc(pz1,a2)*dotrr(pa1,pz2) )
     &                       + dotrc(pz2,a1) * ( dotrc(pz1,a2)*dotrr(pa1,pa2) - dotrc(pa1,a2)*dotrr(pa2,pz1) )
     &                       + dotrc(pz1,a1) * ( dotrc(pz2,a2)*dotrr(pa1,pa2) - dotrc(pa1,a2)*dotrr(pa2,pz2) ) )
     &    - dotcc(a1,a2) * ( dotrc(pz2,z1) * ( dotrc(pa2,z2)*dotrr(pa1,pz1) + dotrc(pa1,z2)*dotrr(pa2,pz1) )
     &                     + dotrc(pz1,z2) * ( dotrc(pa2,z1)*dotrr(pa1,pz2) + dotrc(pa1,z1)*dotrr(pa2,pz2) )
     &                     - dotrr(pz1,pz2) * ( dotrc(pa2,z1)*dotrc(pa1,z2) + dotrc(pa1,z1)*dotrc(pa2,z2) ) )
     &    + dotcc(a1,z1)*dotcc(a2,z2) * ( dotrr(pa1,pz2)*dotrr(pa2,pz1) + dotrr(pa1,pa2)*dotrr(pz1,pz2) ) 
     &    + dotcc(a1,z2)*dotcc(a2,z1) * ( dotrr(pa1,pz1)*dotrr(pa2,pz2) + dotrr(pa1,pa2)*dotrr(pz1,pz2) ) 
     &    - dotcc(a2,z2) * ( dotrc(pz1,a1) * ( dotrc(pz2,z1)*dotrr(pa1,pa2) + dotrc(pa2,z1)*dotrr(pa1,pz2) )
     &                     + dotrc(pz2,a1) * ( dotrc(pa1,z1)*dotrr(pa2,pz1) - dotrc(pa2,z1)*dotrr(pa1,pz1) )
     &                     + dotrc(pa2,a1) * ( dotrc(pa1,z1)*dotrr(pz1,pz2) - dotrc(pz2,z1)*dotrr(pa1,pz1) ) )
     &    - dotcc(a1,z2) * ( dotrc(pz2,a2) * ( dotrc(pa2,z1)*dotrr(pa1,pz1) - dotrc(pa1,z1)*dotrr(pa2,pz1) )
     &                     + dotrc(pz1,a2) * ( dotrc(pz2,z1)*dotrr(pa1,pa2) + dotrc(pa1,z1)*dotrr(pa2,pz2) )
     &                     + dotrc(pa1,a2) * ( dotrc(pa2,z1)*dotrr(pz1,pz2) - dotrc(pz2,z1)*dotrr(pa2,pz1) ) )
     &    - dotcc(a2,z1) * ( dotrc(pz2,a1) * ( dotrc(pz1,z2)*dotrr(pa1,pa2) + dotrc(pa2,z2)*dotrr(pa1,pz1) )
     &                     + dotrc(pz1,a1) * ( dotrc(pa1,z2)*dotrr(pa2,pz2) - dotrc(pa2,z2)*dotrr(pa1,pz2) )
     &                     + dotrc(pa2,a1) * ( dotrc(pa1,z2)*dotrr(pz1,pz2) - dotrc(pz1,z2)*dotrr(pa1,pz2) ) )
     &    - dotcc(a1,z1) * ( dotrc(pz2,a2) * ( dotrc(pz1,z2)*dotrr(pa1,pa2) + dotrc(pa1,z2)*dotrr(pa2,pz1) )
     &                     + dotrc(pz1,a2) * ( dotrc(pa2,z2)*dotrr(pa1,pz2) - dotrc(pa1,z2)*dotrr(pa2,pz2) )
     &                     + dotrc(pa1,a2) * ( dotrc(pa2,z2)*dotrr(pz1,pz2) - dotrc(pz1,z2)*dotrr(pa2,pz2) ) ) )
      endif

      return
      end


*********************************************************************************************************
CCCCCCCCCCCCCCC
C AAAZ-vertex C
CCCCCCCCCCCCCCC

      subroutine aaaz_anomal4(a1,a2,a3,z, aaaz)
c
c This subroutine computes the aaaz-vertex with anomalous couplings
c where all the bosons are defined by the flowing-OUT quantum number.
c     
      implicit none
      double complex z(6),a1(6),a2(6),a3(6),aaaz 
      double complex dotcc, dotrc, epsccrr, epscccr, epscccc
      double precision pz(0:3),pa1(0:3),pa2(0:3),pa3(0:3), dotrr
      double precision cos2w, cosw, sinw

      external dotcc, dotrc, dotrr, epsccrr, epscccr, epscccc

#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"   

      cos2w = 1.0d0 - sin2w
      cosw = sqrt(cos2w)
      sinw = sqrt(sin2w)


      pz(0) = dble( z(5))
      pz(1) = dble( z(6))
      pz(2) = dimag(z(6))
      pz(3) = dimag(z(5))
      pa1(0) = dble( a1(5))
      pa1(1) = dble( a1(6))
      pa1(2) = dimag(a1(6))
      pa1(3) = dimag(a1(5))
      pa2(0) = dble( a2(5))
      pa2(1) = dble( a2(6))
      pa2(2) = dimag(a2(6))
      pa2(3) = dimag(a2(5))
      pa3(0) = dble( a3(5))
      pa3(1) = dble( a3(6))
      pa3(2) = dimag(a3(6))
      pa3(3) = dimag(a3(5))

      aaaz = dcmplx(0.,0.)

      if (ft0 .ne. 0) then
        aaaz = aaaz + ft0 * 8. * gw**4 * sinw**3 * cosw *
     &    ( ( dotrc(pz,a1)*dotrc(pa1,z) - dotcc(a1,z)*dotrr(pa1,pz) ) 
     &            * ( dotrc(pa3,a2)*dotrc(pa2,a3) - dotcc(a2,a3)*dotrr(pa2,pa3) )
     &    + ( dotrc(pa3,a1)*dotrc(pa1,a3) - dotcc(a1,a3)*dotrr(pa1,pa3) ) 
     &            * ( dotrc(pz,a2)*dotrc(pa2,z) - dotcc(a2,z)*dotrr(pa2,pz) )  
     &    + ( dotrc(pa2,a1)*dotrc(pa1,a2) - dotcc(a1,a2)*dotrr(pa1,pa2) ) 
     &            * ( dotrc(pz,a3)*dotrc(pa3,z) - dotcc(a3,z)*dotrr(pa3,pz) ) )
      endif

      if (ft1 .ne. 0) then
        aaaz = aaaz + ft1 * 8. * gw**4 * sinw**3 * cosw *
     &    ( ( dotrc(pz,a1)*dotrc(pa1,z) - dotcc(a1,z)*dotrr(pa1,pz) ) 
     &            * ( dotrc(pa3,a2)*dotrc(pa2,a3) - dotcc(a2,a3)*dotrr(pa2,pa3) )
     &    + ( dotrc(pa3,a1)*dotrc(pa1,a3) - dotcc(a1,a3)*dotrr(pa1,pa3) ) 
     &            * ( dotrc(pz,a2)*dotrc(pa2,z) - dotcc(a2,z)*dotrr(pa2,pz) )  
     &    + ( dotrc(pa2,a1)*dotrc(pa1,a2) - dotcc(a1,a2)*dotrr(pa1,pa2) ) 
     &            * ( dotrc(pz,a3)*dotrc(pa3,z) - dotcc(a3,z)*dotrr(pa3,pz) ) )
      endif

      if (ft2 .ne. 0) then
        aaaz = aaaz + ft2 * 2. * gw**4 * sinw**3 * cosw *
     &    ( dotrc(pa3,a1)*dotrc(pz,a2)*dotrc(pa2,a3)*dotrc(pa1,z)
     &    + dotrc(pa2,a1)*dotrc(pa3,a2)*dotrc(pz,a3)*dotrc(pa1,z)
     &    + dotrc(pz,a1)*dotrc(pa3,a2)*dotrc(pa1,a3)*dotrc(pa2,z)
     &    + dotrc(pa3,a1)*dotrc(pa1,a2)*dotrc(pz,a3)*dotrc(pa2,z)
     &    + dotrc(pa2,a1)*dotrc(pz,a2)*dotrc(pa1,a3)*dotrc(pa3,z)
     &    + dotrc(pz,a1)*dotrc(pa1,a2)*dotrc(pa2,a3)*dotrc(pa3,z)
     &    + dotcc(a1,a2)*dotcc(a3,z) * ( dotrr(pa1,pz)*dotrr(pa2,pa3) + dotrr(pa1,pa3)*dotrr(pa2,pz) )
     &    + dotcc(a3,z) * ( - dotrc(pa2,a1) * ( dotrc(pz,a2)*dotrr(pa1,pa3) + dotrc(pa3,a2)*dotrr(pa1,pz) )
     &                      + dotrc(pz,a1) * ( dotrc(pa3,a2)*dotrr(pa1,pa2) - dotrc(pa1,a2)*dotrr(pa2,pa3) )
     &                      + dotrc(pa3,a1) * ( dotrc(pz,a2)*dotrr(pa1,pa2) - dotrc(pa1,a2)*dotrr(pa2,pz) ) )
     &    - dotcc(a1,a2) * ( dotrc(pz,a3) * ( dotrc(pa2,z)*dotrr(pa1,pa3) + dotrc(pa1,z)*dotrr(pa2,pa3) )
     &                     + dotrc(pa3,z) * ( dotrc(pa2,a3)*dotrr(pa1,pz) + dotrc(pa1,a3)*dotrr(pa2,pz) )
     &                     - dotrr(pa3,pz) * ( dotrc(pa2,a3)*dotrc(pa1,z) + dotrc(pa1,a3)*dotrc(pa2,z) ) )
     &    + dotcc(a1,a3)*dotcc(a2,z) * ( dotrr(pa1,pz)*dotrr(pa2,pa3) + dotrr(pa1,pa2)*dotrr(pa3,pz) )
     &    + dotcc(a1,z)*dotcc(a2,a3) * ( dotrr(pa1,pa3)*dotrr(pa2,pz) + dotrr(pa1,pa2)*dotrr(pa3,pz) )
     &    - dotcc(a2,z) * ( + dotrc(pa3,a1) * ( dotrc(pz,a3)*dotrr(pa1,pa2) + dotrc(pa2,a3)*dotrr(pa1,pz) )
     &                      + dotrc(pz,a1) * ( dotrc(pa1,a3)*dotrr(pa2,pa3) - dotrc(pa2,a3)*dotrr(pa1,pa3) )
     &                      + dotrc(pa2,a1) * ( dotrc(pa1,a3)*dotrr(pa3,pz) - dotrc(pz,a3)*dotrr(pa1,pa3) ) )
     &    - dotcc(a1,z) * ( + dotrc(pz,a2) * ( dotrc(pa2,a3)*dotrr(pa1,pa3) - dotrc(pa1,a3)*dotrr(pa2,pa3) )
     &                      + dotrc(pa3,a2) * ( dotrc(pz,a3)*dotrr(pa1,pa2) + dotrc(pa1,a3)*dotrr(pa2,pz) )
     &                      + dotrc(pa1,a2) * ( dotrc(pa2,a3)*dotrr(pa3,pz) - dotrc(pz,a3)*dotrr(pa2,pa3) ) )
     &    - dotcc(a2,a3) * ( + dotrc(pz,a1) * ( dotrc(pa3,z)*dotrr(pa1,pa2) + dotrc(pa2,z)*dotrr(pa1,pa3) )
     &                       + dotrc(pa3,a1) * ( dotrc(pa1,z)*dotrr(pa2,pz) - dotrc(pa2,z)*dotrr(pa1,pz) )
     &                       + dotrc(pa2,a1) * ( dotrc(pa1,z)*dotrr(pa3,pz) - dotrc(pa3,z)*dotrr(pa1,pz) ) )
     &    - dotcc(a1,a3) * ( + dotrc(pz,a2) * ( dotrc(pa3,z)*dotrr(pa1,pa2) + dotrc(pa1,z)*dotrr(pa2,pa3) )
     &                       + dotrc(pa3,a2) * ( dotrc(pa2,z)*dotrr(pa1,pz) - dotrc(pa1,z)*dotrr(pa2,pz) )
     &                       + dotrc(pa1,a2) * ( dotrc(pa2,z)*dotrr(pa3,pz) - dotrc(pa3,z)*dotrr(pa2,pz) ) ) )
      endif

      if (ft5 .ne. 0) then
        aaaz = aaaz + ft5 * 2. * gw**4 * ( sinw**3 * cosw - sinw**5 / cosw ) *
     &   ( ( dotrc(pz,a1)*dotrc(pa1,z) - dotcc(a1,z)*dotrr(pa1,pz) ) 
     &            * ( dotrc(pa3,a2)*dotrc(pa2,a3) - dotcc(a2,a3)*dotrr(pa2,pa3) )
     &   + ( dotrc(pa3,a1)*dotrc(pa1,a3) - dotcc(a1,a3)*dotrr(pa1,pa3) ) 
     &            * ( dotrc(pz,a2)*dotrc(pa2,z) - dotcc(a2,z)*dotrr(pa2,pz) )
     &   + ( dotrc(pa2,a1)*dotrc(pa1,a2) - dotcc(a1,a2)*dotrr(pa1,pa2) ) 
     &            * ( dotrc(pz,a3)*dotrc(pa3,z) - dotcc(a3,z)*dotrr(pa3,pz) ) )
      endif

      if (ft6 .ne. 0) then
        aaaz = aaaz + ft6 * 2. * gw**4 * ( sinw**3 * cosw - sinw**5 / cosw ) *
     &   ( ( dotrc(pz,a1)*dotrc(pa1,z) - dotcc(a1,z)*dotrr(pa1,pz) ) 
     &            * ( dotrc(pa3,a2)*dotrc(pa2,a3) - dotcc(a2,a3)*dotrr(pa2,pa3) )
     &   + ( dotrc(pa3,a1)*dotrc(pa1,a3) - dotcc(a1,a3)*dotrr(pa1,pa3) ) 
     &            * ( dotrc(pz,a2)*dotrc(pa2,z) - dotcc(a2,z)*dotrr(pa2,pz) )
     &   + ( dotrc(pa2,a1)*dotrc(pa1,a2) - dotcc(a1,a2)*dotrr(pa1,pa2) ) 
     &            * ( dotrc(pz,a3)*dotrc(pa3,z) - dotcc(a3,z)*dotrr(pa3,pz) ) )
      endif

      if (ft7 .ne. 0) then
        aaaz = aaaz + ft7 / 2. * gw**4 * ( sinw**3 * cosw - sinw**5 / cosw ) *
     &    ( dotrc(pz,a1) * dotrc(pa3,a2) * dotrc(pa1,a3) * dotrc(pa2,z)
     &    - dotcc(a2,a3) * dotrc(pz,a1) * dotrc(pa2,z) * dotrr(pa1,pa3)
     &    - dotcc(a1,a2) * dotrc(pz,a3) * dotrc(pa2,z) * dotrr(pa1,pa3)
     &    - dotcc(a1,a3) * dotrc(pa3,a2) * dotrc(pa2,z) * dotrr(pa1,pz)
     &    + dotrc(pz,a1) * dotrc(pa1,a2) * dotrc(pa2,a3) * dotrc(pa3,z)
     &    + dotcc(a3,z) * dotrc(pz,a1) * dotrc(pa3,a2) * dotrr(pa1,pa2)
     &    - dotcc(a1,z) * dotrc(pa3,a2) * dotrc(pz,a3) * dotrr(pa1,pa2)
     &    - dotcc(a2,a3) * dotrc(pz,a1) * dotrc(pa3,z) * dotrr(pa1,pa2)
     &    - dotcc(a1,a3) * dotrc(pz,a2) * dotrc(pa3,z) * dotrr(pa1,pa2)
     &    + dotcc(a2,z) * dotrc(pz,a1) * dotrc(pa2,a3) * dotrr(pa1,pa3)
     &    - dotcc(a1,z) * dotrc(pz,a2) * dotrc(pa2,a3) * dotrr(pa1,pa3)
     &    - dotcc(a1,a2) * dotrc(pa2,a3) * dotrc(pa3,z) * dotrr(pa1,pz)
     &    - dotcc(a3,z) * dotrc(pz,a1) * dotrc(pa1,a2) * dotrr(pa2,pa3)
     &    - dotcc(a2,z) * dotrc(pz,a1) * dotrc(pa1,a3) * dotrr(pa2,pa3)
     &    + dotcc(a1,z) * dotrc(pz,a2) * dotrc(pa1,a3) * dotrr(pa2,pa3)
     &    + dotcc(a1,z) * dotrc(pa1,a2) * dotrc(pz,a3) * dotrr(pa2,pa3)
     &    - dotcc(a1,a3) * dotrc(pz,a2) * dotrc(pa1,z) * dotrr(pa2,pa3)
     &    - dotcc(a1,a2) * dotrc(pz,a3) * dotrc(pa1,z) * dotrr(pa2,pa3)
     &    + dotcc(a1,a3) * dotcc(a2,z) * dotrr(pa1,pz) * dotrr(pa2,pa3)
     &    + dotcc(a1,a2) * dotcc(a3,z) * dotrr(pa1,pz) * dotrr(pa2,pa3)
     &    - dotcc(a1,z) * dotrc(pa3,a2) * dotrc(pa1,a3) * dotrr(pa2,pz)
     &    + dotcc(a1,a3) * dotrc(pa3,a2) * dotrc(pa1,z) * dotrr(pa2,pz)
     &    + dotcc(a1,a3) * dotrc(pa1,a2) * dotrc(pa3,z) * dotrr(pa2,pz)
     &    - dotcc(a1,a2) * dotrc(pa1,a3) * dotrc(pa3,z) * dotrr(pa2,pz)
     &    + dotcc(a1,z) * dotcc(a2,a3) * dotrr(pa1,pa3) * dotrr(pa2,pz)
     &    + dotcc(a1,a2) * dotcc(a3,z) * dotrr(pa1,pa3) * dotrr(pa2,pz)
     &    + dotrc(pa3,a1) * ( dotrc(pa1,a2)*dotrc(pz,a3)*dotrc(pa2,z)
     &                      + dotcc(a2,a3)*dotrc(pa2,z)*dotrr(pa1,pz)
     &                      + dotrc(pz,a2) * ( dotrc(pa2,a3)*dotrc(pa1,z) + dotcc(a3,z)*dotrr(pa1,pa2) )
     &                      - dotcc(a2,z) * ( dotrc(pz,a3)*dotrr(pa1,pa2) + dotrc(pa2,a3)*dotrr(pa1,pz) )
     &                      - dotrr(pa2,pz) * ( dotrc(pa1,a2)*dotcc(a3,z) + dotrc(pa1,z)*dotcc(a2,a3) ) )
     &    + dotrr(pa3,pz) * ( dotcc(a1,a2) * ( dotrc(pa2,a3)*dotrc(pa1,z) + dotrc(pa1,a3)*dotrc(pa2,z) )
     &                      + dotcc(a1,z) * ( dotcc(a2,a3)*dotrr(pa1,pa2) - dotrc(pa1,a2)*dotrc(pa2,a3) )
     &                      + dotcc(a1,a3) * ( dotcc(a2,z)*dotrr(pa1,pa2) - dotrc(pa1,a2)*dotrc(pa2,z) ) )
     &    + dotrc(pa2,a1) * ( dotrc(pz,a2) * ( dotrc(pa1,a3)*dotrc(pa3,z) - dotcc(a3,z)*dotrr(pa1,pa3) )
     &                      + dotrc(pa3,a2) * ( dotrc(pz,a3)*dotrc(pa1,z) - dotcc(a3,z)*dotrr(pa1,pz) )
     &                      + dotcc(a2,z) * ( dotrc(pz,a3)*dotrr(pa1,pa3) - dotrc(pa1,a3)*dotrr(pa3,pz) )
     &                      + dotcc(a2,a3) * ( dotrc(pa3,z)*dotrr(pa1,pz) - dotrc(pa1,z)*dotrr(pa3,pz) ) ) )
      endif

      if (ft8 .ne. 0) then
        aaaz = aaaz - ft8 * 2. * gw**4 * sinw**5 / cosw *
     &   ( ( dotrc(pz,a1)*dotrc(pa1,z) - dotcc(a1,z)*dotrr(pa1,pz) ) 
     &            * ( dotrc(pa3,a2)*dotrc(pa2,a3) - dotcc(a2,a3)*dotrr(pa2,pa3) )
     &   + ( dotrc(pa3,a1)*dotrc(pa1,a3) - dotcc(a1,a3)*dotrr(pa1,pa3) ) 
     &            * ( dotrc(pz,a2)*dotrc(pa2,z) - dotcc(a2,z)*dotrr(pa2,pz) )
     &   + ( dotrc(pa2,a1)*dotrc(pa1,a2) - dotcc(a1,a2)*dotrr(pa1,pa2) ) 
     &            * ( dotrc(pz,a3)*dotrc(pa3,z) - dotcc(a3,z)*dotrr(pa3,pz) ) )
      endif

      if (ft9 .ne. 0) then
        aaaz = aaaz - ft9 / 2. * gw**4 * sinw**5 / cosw *
     &    ( dotrc(pz,a1) * dotrc(pa3,a2) * dotrc(pa1,a3) * dotrc(pa2,z)
     &    - dotcc(a2,a3) * dotrc(pz,a1) * dotrc(pa2,z) * dotrr(pa1,pa3)
     &    - dotcc(a1,a2) * dotrc(pz,a3) * dotrc(pa2,z) * dotrr(pa1,pa3)
     &    - dotcc(a1,a3) * dotrc(pa3,a2) * dotrc(pa2,z) * dotrr(pa1,pz)
     &    + dotrc(pz,a1) * dotrc(pa1,a2) * dotrc(pa2,a3) * dotrc(pa3,z)
     &    + dotcc(a3,z) * dotrc(pz,a1) * dotrc(pa3,a2) * dotrr(pa1,pa2)
     &    - dotcc(a1,z) * dotrc(pa3,a2) * dotrc(pz,a3) * dotrr(pa1,pa2)
     &    - dotcc(a2,a3) * dotrc(pz,a1) * dotrc(pa3,z) * dotrr(pa1,pa2)
     &    - dotcc(a1,a3) * dotrc(pz,a2) * dotrc(pa3,z) * dotrr(pa1,pa2)
     &    + dotcc(a2,z) * dotrc(pz,a1) * dotrc(pa2,a3) * dotrr(pa1,pa3)
     &    - dotcc(a1,z) * dotrc(pz,a2) * dotrc(pa2,a3) * dotrr(pa1,pa3)
     &    - dotcc(a1,a2) * dotrc(pa2,a3) * dotrc(pa3,z) * dotrr(pa1,pz)
     &    - dotcc(a3,z) * dotrc(pz,a1) * dotrc(pa1,a2) * dotrr(pa2,pa3)
     &    - dotcc(a2,z) * dotrc(pz,a1) * dotrc(pa1,a3) * dotrr(pa2,pa3)
     &    + dotcc(a1,z) * dotrc(pz,a2) * dotrc(pa1,a3) * dotrr(pa2,pa3)
     &    + dotcc(a1,z) * dotrc(pa1,a2) * dotrc(pz,a3) * dotrr(pa2,pa3)
     &    - dotcc(a1,a3) * dotrc(pz,a2) * dotrc(pa1,z) * dotrr(pa2,pa3)
     &    - dotcc(a1,a2) * dotrc(pz,a3) * dotrc(pa1,z) * dotrr(pa2,pa3)
     &    + dotcc(a1,a3) * dotcc(a2,z) * dotrr(pa1,pz) * dotrr(pa2,pa3)
     &    + dotcc(a1,a2) * dotcc(a3,z) * dotrr(pa1,pz) * dotrr(pa2,pa3)
     &    - dotcc(a1,z) * dotrc(pa3,a2) * dotrc(pa1,a3) * dotrr(pa2,pz)
     &    + dotcc(a1,a3) * dotrc(pa3,a2) * dotrc(pa1,z) * dotrr(pa2,pz)
     &    + dotcc(a1,a3) * dotrc(pa1,a2) * dotrc(pa3,z) * dotrr(pa2,pz)
     &    - dotcc(a1,a2) * dotrc(pa1,a3) * dotrc(pa3,z) * dotrr(pa2,pz)
     &    + dotcc(a1,z) * dotcc(a2,a3) * dotrr(pa1,pa3) * dotrr(pa2,pz)
     &    + dotcc(a1,a2) * dotcc(a3,z) * dotrr(pa1,pa3) * dotrr(pa2,pz)
     &    + dotrc(pa3,a1) * ( dotrc(pa1,a2)*dotrc(pz,a3)*dotrc(pa2,z)
     &                      + dotcc(a2,a3)*dotrc(pa2,z)*dotrr(pa1,pz)
     &                      + dotrc(pz,a2) * ( dotrc(pa2,a3)*dotrc(pa1,z) + dotcc(a3,z)*dotrr(pa1,pa2) )
     &                      - dotcc(a2,z) * ( dotrc(pz,a3)*dotrr(pa1,pa2) + dotrc(pa2,a3)*dotrr(pa1,pz) )
     &                      - dotrr(pa2,pz) * ( dotrc(pa1,a2)*dotcc(a3,z) + dotrc(pa1,z)*dotcc(a2,a3) ) )
     &    + dotrr(pa3,pz) * ( dotcc(a1,a2) * ( dotrc(pa2,a3)*dotrc(pa1,z) + dotrc(pa1,a3)*dotrc(pa2,z) )
     &                      + dotcc(a1,z) * ( dotcc(a2,a3)*dotrr(pa1,pa2) - dotrc(pa1,a2)*dotrc(pa2,a3) )
     &                      + dotcc(a1,a3) * ( dotcc(a2,z)*dotrr(pa1,pa2) - dotrc(pa1,a2)*dotrc(pa2,z) ) )
     &    + dotrc(pa2,a1) * ( dotrc(pz,a2) * ( dotrc(pa1,a3)*dotrc(pa3,z) - dotcc(a3,z)*dotrr(pa1,pa3) )
     &                      + dotrc(pa3,a2) * ( dotrc(pz,a3)*dotrc(pa1,z) - dotcc(a3,z)*dotrr(pa1,pz) )
     &                      + dotcc(a2,z) * ( dotrc(pz,a3)*dotrr(pa1,pa3) - dotrc(pa1,a3)*dotrr(pa3,pz) )
     &                      + dotcc(a2,a3) * ( dotrc(pa3,z)*dotrr(pa1,pz) - dotrc(pa1,z)*dotrr(pa3,pz) ) ) )
      endif

      return
      end


*********************************************************************************************************
CCCCCCCCCCCCCCC
C AAAA-vertex C
CCCCCCCCCCCCCCC

      subroutine aaaa_anomal4(a1,a2,a3,a4, aaaa)
c
c This subroutine computes the aaaa-vertex with anomalous couplings
c where all the bosons are defined by the flowing-OUT quantum number.
c     
      implicit none
      double complex a1(6),a2(6),a3(6),a4(6),aaaa 
      double complex dotcc, dotrc, epsccrr, epscccr, epscccc
      double precision pa1(0:3),pa2(0:3),pa3(0:3),pa4(0:3), dotrr
      double precision cos2w, cosw, sinw

      external dotcc, dotrc, dotrr, epsccrr, epscccr, epscccc

#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"   

      cos2w = 1.0d0 - sin2w
      cosw = sqrt(cos2w)
      sinw = sqrt(sin2w)


      pa1(0) = dble( a1(5))
      pa1(1) = dble( a1(6))
      pa1(2) = dimag(a1(6))
      pa1(3) = dimag(a1(5))
      pa2(0) = dble( a2(5))
      pa2(1) = dble( a2(6))
      pa2(2) = dimag(a2(6))
      pa2(3) = dimag(a2(5))
      pa3(0) = dble( a3(5))
      pa3(1) = dble( a3(6))
      pa3(2) = dimag(a3(6))
      pa3(3) = dimag(a3(5))
      pa4(0) = dble( a4(5))
      pa4(1) = dble( a4(6))
      pa4(2) = dimag(a4(6))
      pa4(3) = dimag(a4(5))


      aaaa = dcmplx(0.,0.)

      if (ft0 .ne. 0) then
        aaaa = aaaa + ft0 * 8. * gw**4 * sin2w**2 *
     &    ( ( dotrc(pa4,a1)*dotrc(pa1,a4) - dotcc(a1,a4)*dotrr(pa1,pa4) ) 
     &            * ( dotrc(pa3,a2)*dotrc(pa2,a3) - dotcc(a2,a3)*dotrr(pa2,pa3) )
     &    + ( dotrc(pa3,a1)*dotrc(pa1,a3) - dotcc(a1,a3)*dotrr(pa1,pa3) ) 
     &            * ( dotrc(pa4,a2)*dotrc(pa2,a4) - dotcc(a2,a4)*dotrr(pa2,pa4) )  
     &    + ( dotrc(pa2,a1)*dotrc(pa1,a2) - dotcc(a1,a2)*dotrr(pa1,pa2) ) 
     &            * ( dotrc(pa4,a3)*dotrc(pa3,a4) - dotcc(a3,a4)*dotrr(pa3,pa4) ) )
      endif

      if (ft1 .ne. 0) then
        aaaa = aaaa + ft1 * 8. * gw**4 * sin2w**2 *
     &    ( ( dotrc(pa4,a1)*dotrc(pa1,a4) - dotcc(a1,a4)*dotrr(pa1,pa4) ) 
     &            * ( dotrc(pa3,a2)*dotrc(pa2,a3) - dotcc(a2,a3)*dotrr(pa2,pa3) )
     &    + ( dotrc(pa3,a1)*dotrc(pa1,a3) - dotcc(a1,a3)*dotrr(pa1,pa3) ) 
     &            * ( dotrc(pa4,a2)*dotrc(pa2,a4) - dotcc(a2,a4)*dotrr(pa2,pa4) )  
     &    + ( dotrc(pa2,a1)*dotrc(pa1,a2) - dotcc(a1,a2)*dotrr(pa1,pa2) ) 
     &            * ( dotrc(pa4,a3)*dotrc(pa3,a4) - dotcc(a3,a4)*dotrr(pa3,pa4) ) )
      endif

      if (ft2 .ne. 0) then
        aaaa = aaaa + ft2 * 2. * gw**4 * sin2w**2 *
     &    ( dotrc(pa3,a1)*dotrc(pa4,a2)*dotrc(pa2,a3)*dotrc(pa1,a4)
     &    + dotrc(pa2,a1)*dotrc(pa3,a2)*dotrc(pa4,a3)*dotrc(pa1,a4)
     &    + dotrc(pa4,a1)*dotrc(pa3,a2)*dotrc(pa1,a3)*dotrc(pa2,a4)
     &    + dotrc(pa3,a1)*dotrc(pa1,a2)*dotrc(pa4,a3)*dotrc(pa2,a4)
     &    + dotrc(pa2,a1)*dotrc(pa4,a2)*dotrc(pa1,a3)*dotrc(pa3,a4)
     &    + dotrc(pa4,a1)*dotrc(pa1,a2)*dotrc(pa2,a3)*dotrc(pa3,a4)
     &    + dotcc(a1,a2)*dotcc(a3,a4) * ( dotrr(pa1,pa4)*dotrr(pa2,pa3) + dotrr(pa1,pa3)*dotrr(pa2,pa4) )
     &    + dotcc(a3,a4) * ( - dotrc(pa2,a1) * ( dotrc(pa4,a2)*dotrr(pa1,pa3) + dotrc(pa3,a2)*dotrr(pa1,pa4) )
     &                       + dotrc(pa4,a1) * ( dotrc(pa3,a2)*dotrr(pa1,pa2) - dotrc(pa1,a2)*dotrr(pa2,pa3) )
     &                       + dotrc(pa3,a1) * ( dotrc(pa4,a2)*dotrr(pa1,pa2) - dotrc(pa1,a2)*dotrr(pa2,pa4) ) )
     &    - dotcc(a1,a2) * ( dotrc(pa4,a3) * ( dotrc(pa2,a4)*dotrr(pa1,pa3) + dotrc(pa1,a4)*dotrr(pa2,pa3) )
     &                     + dotrc(pa3,a4) * ( dotrc(pa2,a3)*dotrr(pa1,pa4) + dotrc(pa1,a3)*dotrr(pa2,pa4) )
     &                     - dotrr(pa3,pa4) * ( dotrc(pa2,a3)*dotrc(pa1,a4) + dotrc(pa1,a3)*dotrc(pa2,a4) ) )
     &    + dotcc(a1,a3)*dotcc(a2,a4) * ( dotrr(pa1,pa4)*dotrr(pa2,pa3) + dotrr(pa1,pa2)*dotrr(pa3,pa4) )
     &    + dotcc(a1,a4)*dotcc(a2,a3) * ( dotrr(pa1,pa3)*dotrr(pa2,pa4) + dotrr(pa1,pa2)*dotrr(pa3,pa4) )
     &    - dotcc(a2,a4) * ( + dotrc(pa3,a1) * ( dotrc(pa4,a3)*dotrr(pa1,pa2) + dotrc(pa2,a3)*dotrr(pa1,pa4) )
     &                       + dotrc(pa4,a1) * ( dotrc(pa1,a3)*dotrr(pa2,pa3) - dotrc(pa2,a3)*dotrr(pa1,pa3) )
     &                       + dotrc(pa2,a1) * ( dotrc(pa1,a3)*dotrr(pa3,pa4) - dotrc(pa4,a3)*dotrr(pa1,pa3) ) )
     &    - dotcc(a1,a4) * ( + dotrc(pa4,a2) * ( dotrc(pa2,a3)*dotrr(pa1,pa3) - dotrc(pa1,a3)*dotrr(pa2,pa3) )
     &                       + dotrc(pa3,a2) * ( dotrc(pa4,a3)*dotrr(pa1,pa2) + dotrc(pa1,a3)*dotrr(pa2,pa4) )
     &                       + dotrc(pa1,a2) * ( dotrc(pa2,a3)*dotrr(pa3,pa4) - dotrc(pa4,a3)*dotrr(pa2,pa3) ) )
     &    - dotcc(a2,a3) * ( + dotrc(pa4,a1) * ( dotrc(pa3,a4)*dotrr(pa1,pa2) + dotrc(pa2,a4)*dotrr(pa1,pa3) )
     &                       + dotrc(pa3,a1) * ( dotrc(pa1,a4)*dotrr(pa2,pa4) - dotrc(pa2,a4)*dotrr(pa1,pa4) )
     &                       + dotrc(pa2,a1) * ( dotrc(pa1,a4)*dotrr(pa3,pa4) - dotrc(pa3,a4)*dotrr(pa1,pa4) ) )
     &    - dotcc(a1,a3) * ( + dotrc(pa4,a2) * ( dotrc(pa3,a4)*dotrr(pa1,pa2) + dotrc(pa1,a4)*dotrr(pa2,pa3) )
     &                       + dotrc(pa3,a2) * ( dotrc(pa2,a4)*dotrr(pa1,pa4) - dotrc(pa1,a4)*dotrr(pa2,pa4) )
     &                       + dotrc(pa1,a2) * ( dotrc(pa2,a4)*dotrr(pa3,pa4) - dotrc(pa3,a4)*dotrr(pa2,pa4) ) ) )
      endif

      if (ft5 .ne. 0) then
        aaaa = aaaa + ft5 * 4. * gw**4 * sin2w**2 *
     &   ( ( dotrc(pa1,a4)*dotrc(pa4,a1) - dotcc(a1,a4)*dotrr(pa1,pa4) ) 
     &            * ( dotrc(pa2,a3)*dotrc(pa3,a2) - dotcc(a2,a3)*dotrr(pa2,pa3) )
     &   + ( dotrc(pa1,a3)*dotrc(pa3,a1) - dotcc(a1,a3)*dotrr(pa1,pa3) ) 
     &            * ( dotrc(pa2,a4)*dotrc(pa4,a2) - dotcc(a2,a4)*dotrr(pa2,pa4) )
     &   + ( dotrc(pa1,a2)*dotrc(pa2,a1) - dotcc(a1,a2)*dotrr(pa1,pa2) ) 
     &            * ( dotrc(pa3,a4)*dotrc(pa4,a3) - dotcc(a3,a4)*dotrr(pa3,pa4) ) )
      endif

      if (ft6 .ne. 0) then
        aaaa = aaaa + ft6 * 4. * gw**4 * sin2w**2 *
     &   ( ( dotrc(pa1,a4)*dotrc(pa4,a1) - dotcc(a1,a4)*dotrr(pa1,pa4) ) 
     &            * ( dotrc(pa2,a3)*dotrc(pa3,a2) - dotcc(a2,a3)*dotrr(pa2,pa3) )
     &   + ( dotrc(pa1,a3)*dotrc(pa3,a1) - dotcc(a1,a3)*dotrr(pa1,pa3) ) 
     &            * ( dotrc(pa2,a4)*dotrc(pa4,a2) - dotcc(a2,a4)*dotrr(pa2,pa4) )
     &   + ( dotrc(pa1,a2)*dotrc(pa2,a1) - dotcc(a1,a2)*dotrr(pa1,pa2) ) 
     &            * ( dotrc(pa3,a4)*dotrc(pa4,a3) - dotcc(a3,a4)*dotrr(pa3,pa4) ) )
      endif

      if (ft7 .ne. 0) then
        aaaa = aaaa + ft7 * gw**4 * sin2w**2 *
     &    ( dotrc(pa3,a1) * dotrc(pa4,a2) * dotrc(pa2,a3) * dotrc(pa1,a4)
     &    + dotrc(pa2,a1) * dotrc(pa3,a2) * dotrc(pa4,a3) * dotrc(pa1,a4)
     &    + dotrc(pa4,a1) * dotrc(pa3,a2) * dotrc(pa1,a3) * dotrc(pa2,a4)
     &    + dotrc(pa3,a1) * dotrc(pa1,a2) * dotrc(pa4,a3) * dotrc(pa2,a4)
     &    + dotrc(pa2,a1) * dotrc(pa4,a2) * dotrc(pa1,a3) * dotrc(pa3,a4)
     &    + dotrc(pa4,a1) * dotrc(pa1,a2) * dotrc(pa2,a3) * dotrc(pa3,a4)
     &    + dotcc(a1,a2)*dotcc(a3,a4) * ( dotrr(pa1,pa4)*dotrr(pa2,pa3) + dotrr(pa1,pa3)*dotrr(pa2,pa4) )
     &    + dotcc(a3,a4) * ( - dotrc(pa2,a1) * ( dotrc(pa4,a2)*dotrr(pa1,pa3) + dotrc(pa3,a2)*dotrr(pa1,pa4) )
     &                       + dotrc(pa4,a1) * ( dotrc(pa3,a2)*dotrr(pa1,pa2) - dotrc(pa1,a2)*dotrr(pa2,pa3) )
     &                       + dotrc(pa3,a1) * ( dotrc(pa4,a2)*dotrr(pa1,pa2) - dotrc(pa1,a2)*dotrr(pa2,pa4) ) )
     &    - dotcc(a1,a2) * ( dotrc(pa4,a3) * ( dotrc(pa2,a4)*dotrr(pa1,pa3) + dotrc(pa1,a4)*dotrr(pa2,pa3) )
     &                     + dotrc(pa3,a4) * ( dotrc(pa2,a3)*dotrr(pa1,pa4) + dotrc(pa1,a3)*dotrr(pa2,pa4) )
     &                     - dotrr(pa3,pa4) * ( dotrc(pa2,a3)*dotrc(pa1,a4) + dotrc(pa1,a3)*dotrc(pa2,a4) ) )
     &    + dotcc(a1,a3)*dotcc(a2,a4) * ( dotrr(pa1,pa4)*dotrr(pa2,pa3) + dotrr(pa1,pa2)*dotrr(pa3,pa4) ) 
     &    + dotcc(a1,a4)*dotcc(a2,a3) * ( dotrr(pa1,pa3)*dotrr(pa2,pa4) + dotrr(pa1,pa2)*dotrr(pa3,pa4) ) 
     &    - dotcc(a2,a4) * ( dotrc(pa3,a1) * ( dotrc(pa4,a3)*dotrr(pa1,pa2) + dotrc(pa2,a3)*dotrr(pa1,pa4) )
     &                     + dotrc(pa4,a1) * ( dotrc(pa1,a3)*dotrr(pa2,pa3) - dotrc(pa2,a3)*dotrr(pa1,pa3) )
     &                     + dotrc(pa2,a1) * ( dotrc(pa1,a3)*dotrr(pa3,pa4) - dotrc(pa4,a3)*dotrr(pa1,pa3) ) )
     &    - dotcc(a1,a4) * ( dotrc(pa4,a2) * ( dotrc(pa2,a3)*dotrr(pa1,pa3) - dotrc(pa1,a3)*dotrr(pa2,pa3) )
     &                     + dotrc(pa3,a2) * ( dotrc(pa4,a3)*dotrr(pa1,pa2) + dotrc(pa1,a3)*dotrr(pa2,pa4) )
     &                     + dotrc(pa1,a2) * ( dotrc(pa2,a3)*dotrr(pa3,pa4) - dotrc(pa4,a3)*dotrr(pa2,pa3) ) )
     &    - dotcc(a2,a3) * ( dotrc(pa4,a1) * ( dotrc(pa3,a4)*dotrr(pa1,pa2) + dotrc(pa2,a4)*dotrr(pa1,pa3) )
     &                     + dotrc(pa3,a1) * ( dotrc(pa1,a4)*dotrr(pa2,pa4) - dotrc(pa2,a4)*dotrr(pa1,pa4) )
     &                     + dotrc(pa2,a1) * ( dotrc(pa1,a4)*dotrr(pa3,pa4) - dotrc(pa3,a4)*dotrr(pa1,pa4) ) )
     &    - dotcc(a1,a3) * ( dotrc(pa4,a2) * ( dotrc(pa3,a4)*dotrr(pa1,pa2) + dotrc(pa1,a4)*dotrr(pa2,pa3) )
     &                     + dotrc(pa3,a2) * ( dotrc(pa2,a4)*dotrr(pa1,pa4) - dotrc(pa1,a4)*dotrr(pa2,pa4) )
     &                     + dotrc(pa1,a2) * ( dotrc(pa2,a4)*dotrr(pa3,pa4) - dotrc(pa3,a4)*dotrr(pa2,pa4) ) ) )
      endif

      if (ft8 .ne. 0) then
        aaaa = aaaa + ft8 * 2. * gw**4 * sin2w**2 *
     &   ( ( dotrc(pa1,a4)*dotrc(pa4,a1) - dotcc(a1,a4)*dotrr(pa1,pa4) ) 
     &            * ( dotrc(pa2,a3)*dotrc(pa3,a2) - dotcc(a2,a3)*dotrr(pa2,pa3) )
     &   + ( dotrc(pa1,a3)*dotrc(pa3,a1) - dotcc(a1,a3)*dotrr(pa1,pa3) ) 
     &            * ( dotrc(pa2,a4)*dotrc(pa4,a2) - dotcc(a2,a4)*dotrr(pa2,pa4) )
     &   + ( dotrc(pa1,a2)*dotrc(pa2,a1) - dotcc(a1,a2)*dotrr(pa1,pa2) ) 
     &            * ( dotrc(pa3,a4)*dotrc(pa4,a3) - dotcc(a3,a4)*dotrr(pa3,pa4) ) )
      endif

      if (ft9 .ne. 0) then
        aaaa = aaaa + ft9 / 2. * gw**4 * sin2w**2 *
     &    ( dotrc(pa3,a1) * dotrc(pa4,a2) * dotrc(pa2,a3) * dotrc(pa1,a4)
     &    + dotrc(pa2,a1) * dotrc(pa3,a2) * dotrc(pa4,a3) * dotrc(pa1,a4)
     &    + dotrc(pa4,a1) * dotrc(pa3,a2) * dotrc(pa1,a3) * dotrc(pa2,a4)
     &    + dotrc(pa3,a1) * dotrc(pa1,a2) * dotrc(pa4,a3) * dotrc(pa2,a4)
     &    + dotrc(pa2,a1) * dotrc(pa4,a2) * dotrc(pa1,a3) * dotrc(pa3,a4)
     &    + dotrc(pa4,a1) * dotrc(pa1,a2) * dotrc(pa2,a3) * dotrc(pa3,a4)
     &    + dotcc(a1,a2)*dotcc(a3,a4) * ( dotrr(pa1,pa4)*dotrr(pa2,pa3) + dotrr(pa1,pa3)*dotrr(pa2,pa4) )
     &    + dotcc(a3,a4) * ( - dotrc(pa2,a1) * ( dotrc(pa4,a2)*dotrr(pa1,pa3) + dotrc(pa3,a2)*dotrr(pa1,pa4) )
     &                       + dotrc(pa4,a1) * ( dotrc(pa3,a2)*dotrr(pa1,pa2) - dotrc(pa1,a2)*dotrr(pa2,pa3) )
     &                       + dotrc(pa3,a1) * ( dotrc(pa4,a2)*dotrr(pa1,pa2) - dotrc(pa1,a2)*dotrr(pa2,pa4) ) )
     &    - dotcc(a1,a2) * ( dotrc(pa4,a3) * ( dotrc(pa2,a4)*dotrr(pa1,pa3) + dotrc(pa1,a4)*dotrr(pa2,pa3) )
     &                     + dotrc(pa3,a4) * ( dotrc(pa2,a3)*dotrr(pa1,pa4) + dotrc(pa1,a3)*dotrr(pa2,pa4) )
     &                     - dotrr(pa3,pa4) * ( dotrc(pa2,a3)*dotrc(pa1,a4) + dotrc(pa1,a3)*dotrc(pa2,a4) ) )
     &    + dotcc(a1,a3)*dotcc(a2,a4) * ( dotrr(pa1,pa4)*dotrr(pa2,pa3) + dotrr(pa1,pa2)*dotrr(pa3,pa4) ) 
     &    + dotcc(a1,a4)*dotcc(a2,a3) * ( dotrr(pa1,pa3)*dotrr(pa2,pa4) + dotrr(pa1,pa2)*dotrr(pa3,pa4) ) 
     &    - dotcc(a2,a4) * ( dotrc(pa3,a1) * ( dotrc(pa4,a3)*dotrr(pa1,pa2) + dotrc(pa2,a3)*dotrr(pa1,pa4) )
     &                     + dotrc(pa4,a1) * ( dotrc(pa1,a3)*dotrr(pa2,pa3) - dotrc(pa2,a3)*dotrr(pa1,pa3) )
     &                     + dotrc(pa2,a1) * ( dotrc(pa1,a3)*dotrr(pa3,pa4) - dotrc(pa4,a3)*dotrr(pa1,pa3) ) )
     &    - dotcc(a1,a4) * ( dotrc(pa4,a2) * ( dotrc(pa2,a3)*dotrr(pa1,pa3) - dotrc(pa1,a3)*dotrr(pa2,pa3) )

     &                     + dotrc(pa3,a2) * ( dotrc(pa4,a3)*dotrr(pa1,pa2) + dotrc(pa1,a3)*dotrr(pa2,pa4) )
     &                     + dotrc(pa1,a2) * ( dotrc(pa2,a3)*dotrr(pa3,pa4) - dotrc(pa4,a3)*dotrr(pa2,pa3) ) )
     &    - dotcc(a2,a3) * ( dotrc(pa4,a1) * ( dotrc(pa3,a4)*dotrr(pa1,pa2) + dotrc(pa2,a4)*dotrr(pa1,pa3) )
     &                     + dotrc(pa3,a1) * ( dotrc(pa1,a4)*dotrr(pa2,pa4) - dotrc(pa2,a4)*dotrr(pa1,pa4) )
     &                     + dotrc(pa2,a1) * ( dotrc(pa1,a4)*dotrr(pa3,pa4) - dotrc(pa3,a4)*dotrr(pa1,pa4) ) )
     &    - dotcc(a1,a3) * ( dotrc(pa4,a2) * ( dotrc(pa3,a4)*dotrr(pa1,pa2) + dotrc(pa1,a4)*dotrr(pa2,pa3) )
     &                     + dotrc(pa3,a2) * ( dotrc(pa2,a4)*dotrr(pa1,pa4) - dotrc(pa1,a4)*dotrr(pa2,pa4) )
     &                     + dotrc(pa1,a2) * ( dotrc(pa2,a4)*dotrr(pa3,pa4) - dotrc(pa3,a4)*dotrr(pa2,pa4) ) ) )
      endif

      return
      end


********************************************************************************
********************************************************************************

** This subroutine applies the formfactor for the dim-8 and dim-6 operators
** for anomalous gauge boson couplings.
** Input via anomV.dat.
** The coefficients FB, FW and FWWW 
** (which parametrise the WWZ and WWgamma couplings) can have individual 
** formfactors

      subroutine anomal_formfactor(p1,p2,p3,p4)

      implicit none

#include "VBFNLO/utilities/an_couplings.inc"   
*#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/mssm.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/process.inc"

** Input variables: p1,p2,p3,p4 are the boson momenta
      real*8 p1(0:3),p2(0:3),p3(0:3),p4(0:3)

      double precision ffvalue_fwww, ffvalue_fw, ffvalue_fb, ffvalue, ffvalue_KZ, ffvalue_KA

** K-Matrix variables
      double complex DA00,DA11,DA02,DA20,DA22

      if (formfact .OR. kmatrix) then  
                            ! In this part the center of mass energy
                            !ss of the electroweak systems gets calculated. This is needed both for
                            !formfactors and k-matrix unitarisation.
** energy scale: 
        SELECT CASE(procID) ! different cases for VVjj via VBF, diboson and triboson processes
                            ! here also different form factors depending on the vector boson momenta could be implemented
                            ! VVj and Vjj in VBF processes use own implementation for anomalous couplings and form factors!
        ! VVjj in VBF:
        ! p1, p2 are momenta of incoming bosons, p3, p4 are momenta of outgoing bosons
        CASE(WPWMjj,ZZjj_ll,ZZjj_lnu,WPZjj,WMZjj,WPWPjj,WMWMjj,WPAjj,WMAjj,
     &       WPhadWMjj, WPWMhadjj, WPhadZjj, WPZhadjj, WMhadZjj, WMZhadjj, ZZhadjj, WPhadWPjj, WMhadWMjj,
     &       ZAjj,ZAjj_n)

           ss=(p3(0)+p4(0))**2-(p3(1)+p4(1))**2-(p3(2)+p4(2))**2-(p3(3)+p4(3))**2 ! two outgoing boson momenta

        ! Vjj in VBF:
        ! p1, p2 are momenta of incoming bosons, p3, p4 dummy variables
        CASE(Zjj_l,Zjj_nu,WPjj,WMjj)

           ss=(p1(0)+p2(0))**2-(p1(1)+p2(1))**2-(p1(2)+p2(2))**2-(p1(3)+p2(3))**2 ! two incoming boson momenta

        ! Triboson
        ! p1 is momentum of incoming boson, p2, p3, p4 are momenta of outgoing bosons
        CASE(WWZ,ZZZ,ZZWP,ZZWM,WWWP,WWWM,WWA,ZZA,ZZnA,WPZA,WMZA,WPAA,WMAA,ZAA,ZnAA,AAA,
     &       WPAAj,WMAAj,
     &       WPhadWMZ, WPWMhadZ, WWZhad, ZZhadWP, ZZWPhad, ZZhadWM, ZZWMhad,
     &       WPhadWMWP, WPWMhadWP, WMhadWPWM, WMWPhadWM, ZZZhad,
     &       WPhadWMA, WPWMhadA, ZZhadA, WPhadZA, WPZhadA, WMhadZA, WMZhadA)

           ss=(p1(0))**2-(p1(1))**2-(p1(2))**2-(p1(3))**2 ! one incoming boson momentum

        ! Diboson
        ! p1 is momentum of incoming boson, p2, p3 are momenta of outgoing bosons, p4 is dummy variable
        CASE(WW,WPZ,WMZ,GFWW,
     &       WPhadWMlep,WPlepWMhad,WPhadZ,WPZhad,WMhadZ,WMZhad,GFWPhadWMlep,GFWPlepWMhad,
     &       WPH,WMH,WPHJ,WMHJ)

           ss=(p1(0))**2-(p1(1))**2-(p1(2))**2-(p1(3))**2 ! one incoming boson momentum

        CASE DEFAULT
           print *, "form factor not available for this process!"
        END SELECT

        if (kmatrix) then !K- Matrix unitarisation for L_S0 and L_S1

         DA00 = -((7.*(fs0+fs2) + 11.*fS1)**2*ss**4)/(6.*((0,192)*Pi + (7.*(fs0+fs2) + 11.*fS1)*ss**2)) 
         DA02 = -((2.*(fs0+fs2) + fS1)**2*ss**4)/(30.*((0,960)*Pi + (2.*(fs0+fs2) + fS1)*ss**2))
         DA11 = -(((fs0+fs2) - 2.*fS1)**2*ss**4)/(12.*((0,384)*Pi + ((fs0+fs2) - 2.*fS1)*ss**2))
         DA20 = -((2.*(fs0+fs2) + fS1)**2*ss**4)/(3.*((0,96)*Pi + (2.*(fs0+fs2) + fS1)*ss**2))
         DA22 = -(((fs0+fs2) + 2.*fS1)**2*ss**4)/(60.*((0,1920)*Pi + ((fs0+fs2) + 2.*fS1)*ss**2))
 
        SELECT CASE(procID)

!       ** W+.W- **
        CASE(WPWMjj, WPhadWMjj, WPWMhadjj) !these processes can have z-bosons in
            !the initial state. Therefore both wwzz and wwww need unitarisation.

         fswp1wm1 = fs1 + fs2 +  1./(3. * ss**2) * (2. * DA00 + DA20)
     &                        - 10./(3. * ss**2) * (2. * DA02 + DA22)
         fswp1wp2 = fs0 +        1./(2. * ss**2) * (10. * DA02 + 3. * DA11 + 5. * DA22)
         fswp1wm2 = fs1 + fs2 +  1./(1. * ss**2) * (10. * DA02 - 3. * DA11 + 5. * DA22)

         fswmwp = fs1 + 2./(3.*ss**2) * (DA00 - DA20) - 20./(3.*ss**2) * (DA02 - DA22)
         fswmz1 = fS0 + fS2 + 20./(ss**2) * (DA02 - DA22)
         fswmz2 = fswmz1        

!       ** W+.W+ **
        CASE(WPWPjj, WMWMjj, WPhadWPjj, WMhadWMjj)
         
         fswp1wp2 = fs0 +        1./ss**2 * (DA20 - 10. * DA22) 
         fswp1wm2 = fs1 + fs2 + 30./ss**2 * DA22
         fswp1wm1 = fswp1wm2

!        ** W.Z **         
        CASE(WPZjj, WMZjj, WPhadZjj, WPZhadjj, WMhadZjj, WMZhadjj)

         fswmwp = fs1 -       3./ss**2 * DA11 + 15./ss**2 * DA22
         fswmz1 = fs0 + fs2 + 2./ss**2 * DA20 - 20./ss**2 * DA22
         fswmz2 = fs0 + fs2 + 6./ss**2 * DA11 + 30./ss**2 * DA22

!       ** Z.Z **  
        CASE(ZZjj_ll, ZZjj_lnu, ZZhadjj) !these processes can have w-bosons in
            !the initial state. Therefore both wwzz and zzzz need
            !unitarisation, just like in the W+.W- case.
         
         fswmwp = fs1 + 2./(3.*ss**2) * (DA00 - DA20) - 20./(3.*ss**2) * (DA02 - DA22)
         fswmz1 = fS0 + fS2 + 20./(ss**2) * (DA02 - DA22)
         fswmz2 = fswmz1        

         fsz1z2 = fS0 + fS1 + fs2 + 2./(3.*ss**2) * (DA00 + 2.*DA20) - 20./(3.*ss**2) * (DA02 + 2.*DA22)
         fsz1z3 = fS0 + fS1 + fs2 + 10./ss**2 * (DA02 + 2. * DA22) 
         fsz1z4 = fsz1z3
     
        CASE DEFAULT
           print *, "K- matrix unitarisation not available for this process!"
        END SELECT
        endif  !if (kmatrix)

        if (formfact) then 
** 'universal' formfactors:
         ffvalue=(1./(1.+ss/ffmassscale2))**DBLE(ffexponent)

** individual formfactors for fw,fb and fwww
         if (formfac_ind) then
           ffvalue_fwww = (1./(1. + ss/massscale2FWWW))**DBLE(ffexpFWWW)
           ffvalue_fw = (1./(1. + ss/massscale2FW))**DBLE(ffexpFW)
           if (trianom .eq. 1) then
              ffvalue_fb = (1./(1. + ss/massscale2FB))**DBLE(ffexpFB)
           else if (fb_0 .ne. 0d0) then
              ffvalue_fb = 2d0*aDkappa0_0/(1d0 + ss/
     &             massscale2KA)**dble(ffexpKA) - MW2*fw_0*ffvalue_fw
              ffvalue_fb = ffvalue_fb/(2d0*aDkappa0_0 - MW2*fw_0)
           else
              ffvalue_fb = 1d0
           end if
         else
           ffvalue_fb = ffvalue
           ffvalue_fw = ffvalue
           ffvalue_fwww = ffvalue
         end if


**  dim-8 operators:
         fs0=fs0_0*ffvalue
         fs1=fs1_0*ffvalue
         fs2=fs2_0*ffvalue
         ft0=ft0_0*ffvalue
         ft1=ft1_0*ffvalue
         ft2=ft2_0*ffvalue
         ft5=ft5_0*ffvalue
         ft6=ft6_0*ffvalue
         ft7=ft7_0*ffvalue
         ft8=ft8_0*ffvalue
         ft9=ft9_0*ffvalue
         fm0=fm0_0*ffvalue
         fm1=fm1_0*ffvalue
         fm2=fm2_0*ffvalue
         fm3=fm3_0*ffvalue
         fm4=fm4_0*ffvalue
         fm5=fm5_0*ffvalue
         fm6=fm6_0*ffvalue
         fm7=fm7_0*ffvalue
 
**  dim-6 operators:
 
         fwww = fwww_0*ffvalue_fwww
 
         fw = fw_0*ffvalue_fw
 
         fb = fb_0*ffvalue_fb
 
         fbw=fbw_0*ffvalue
         fdw=fdw_0*ffvalue
         fww=fww_0*ffvalue
         fbb=fbb_0*ffvalue
         fwwtilde=fwwtilde_0*ffvalue
         fbwtilde=fbwtilde_0*ffvalue
         fbbtilde=fbbtilde_0*ffvalue
         fwtilde=fwtilde_0*ffvalue
         fbtilde=fbtilde_0*ffvalue
         fwwwtilde=fwwwtilde_0*ffvalue
         fdwtilde=fdwtilde_0*ffvalue
         fphi1=fphi1_0*ffvalue


**  Now moving to the trianom=2 parametrisation (needed for Vjj in VBF)
          if (trianom .eq. 2) then
             ffvalue_KA = (1./(1.+ss/massscale2KA))**DBLE(ffexpKA)
             ffvalue_KZ = (1./(1.+ss/massscale2KZ))**DBLE(ffexpKZ)
          else
             ! numerator is multiplied by appropriate formfactors
             if ( (fw_0 + fb_0) .ne. 0d0 ) then
               ffvalue_KA = (fw + fb)/(fw_0 + fb_0)
             else
               ffvalue_KA = 0
             endif
             if ( (fw_0*CW2 - SW2*fb_0) .ne. 0d0 ) then
               ffvalue_KZ = (CW2*fw - SW2*fb)/(fw_0*CW2 - SW2*fb_0)
             else
               ffvalue_KZ = 0
             endif
          end if
          if (lambda0 .ne. 0D0) then
             lambda0 = lambda0_0
             lambda0 = lambda0*ffvalue_fwww
          end if
          if (zDg0 .ne. 0D0) then
             zDg0 = zDg0_0
             zDg0 = zDg0*ffvalue_fw
          end if
          if ((aDkappa0 .ne. 0d0) .or. (zDkappa0 .ne. 0d0)) then
             aDkappa0 = aDkappa0_0
             aDkappa0 = aDkappa0*ffvalue_KA
             zDkappa0 = zDkappa0_0
             zDkappa0 = zDkappa0*ffvalue_KZ
          end if

       endif 
      endif


      end
