c  for V1V2 --> spin-2 particle --> V3V4

c by Dieter Zeppenfeld and Jessica Frank

      subroutine vvsp2tovv(idfsvv,p1,p2,q1,q2,eps1,eps2,
     1           wwsp2vv,wwsp2vvx,zzsp2vv,azsp2vv,zasp2vv,aasp2vv)

c  This subroutine calculates leptonic tensors for different V3V4 
c  or V3V4 -> 4 leptons (or neutrinos) final states via a spin-2 resonance
c  (only resonant spin-2 contributions (singlet or triplet), no SM or non-resonant spin-2 diagrams)
c
c  Input:   idfsvv:    ID of final-state V3,V4:
c                      1: gamma gamma (-> 4l) 
c                      2: WW (-> lnulnu)
c                      3: ZZ (-> 4l/nu)
c                      4: gamma Z (-> 4l/nu) 
c
c           p1,p2:     incoming momenta of V1,V2,
c           q1,q2:     outgoing momenta of V3,V4,
c           eps1,eps2: V3,V4 polarisation vectors (in case of final state photons) 
c                      or V3,V4 currents from final-state leptons (in case of V3V4 -> 4 leptons (or neutrinos))
c
c  Output: wwsp2vv(mu,nu) etc.: leptonic tensor for V1V2 = WW etc.
C 
      IMPLICIT NONE

C  input and output variables
      INTEGER idfsvv
      REAL*8  p1(0:3), p2(0:3), q1(0:3), q2(0:3)
      complex*16 eps1(0:3), eps2(0:3)
      double complex wwsp2vv(0:3,0:3), zzsp2vv(0:3,0:3), azsp2vv(0:3,0:3), 
     1           wwsp2vvx(0:3,0:3), zasp2vv(0:3,0:3), aasp2vv(0:3,0:3)

c for histograms
      REAL*8 p1inimom0, p1inimom1, p1inimom2, p1inimom3
      REAL*8 p2inimom0, p2inimom1, p2inimom2, p2inimom3

      COMMON /histpinimom/ p1inimom0, p1inimom1, p1inimom2, p1inimom3,
     &                     p2inimom0, p2inimom1, p2inimom2, p2inimom3

C  local variables
      INTEGER  mu,nu,i
      real*8 k(0:3), q1mq2(0:3), p1mp2(0:3) 
      real*8 p1mp2sq, q1mq2sq, q12p12, kp1mp2, kq1mq2, 
     1          ax1, ax5, ax6, ax7, mtsq, sbar, p1sq,p2sq, 
     2       spin2massgen, spin2widthgen
      complex*16 t11(0:3,0:3), t21(0:3,0:3), t12(0:3,0:3), t22(0:3,0:3)
      complex*16 e1p1mp2, e2p1mp2, keps1, keps2, e1e2, ax4, 
     1          resmunu, reskk, reskmu1, reskmu2, reskmu3, reskmu4, 
     2          res12, resknu1, resknu2, resknu3, resknu4, 
     3          reseps1mua, reseps1mub, reseps1nua, reseps1nub,  
     4          reseps2mua, reseps2mub, reseps2nua, reseps2nub,
     5          respp, respq, resqp, resqq, formfac
      REAL*8     f1ww, f2ww, f1zz,f2zz,f1az,f1aa
      double precision cos2w, cosw, sinw, gw1, v
      complex*16 zmw2, zmz2
      logical lreset, lk11, lk21, lk12, lk22
      data lreset /.true./
      save  zmw2,zmz2,f1ww, f2ww, f1zz,f2zz,f1az,f1aa
      save  lreset, lk11, lk21, lk12, lk22

#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/spin2coupl.inc"

      p1inimom0=p1(0)
      p1inimom1=p1(1)
      p1inimom2=p1(2)
      p1inimom3=p1(3)
      p2inimom0=p2(0)
      p2inimom1=p2(1)
      p2inimom2=p2(2)
      p2inimom3=p2(3)

      if (lreset) then
         cos2w = 1.0d0 - sin2w
         cosw = sqrt(cos2w)
         sinw = sqrt(sin2w)
         gw1 = ee / cosw
         v = 2.0d0*wmass/gw

         f1ww = 2.0d0*f2/lambda
         f2ww = f5*gw**2*v**2/(2.0d0*lambda) + f6*gw**2*v**2/(4.0d0*lambda)
         f1zz = 2.0d0*(f2*cos2w+f1*sin2w)/lambda - 2.0d0*f7*cosw*sinw/lambda
         f2zz = f5*(gw**2+gw1**2)*v**2/(2.0d0*lambda) - f6*(gw**2+gw1**2)*v**2/(4.0d0*lambda)
         f1aa = 2.0d0*(f1*cos2w+f2*sin2w)/lambda + 2.0d0*f7*cosw*sinw/lambda
         f1az = 2.0d0*cosw*sinw*(f2-f1)/lambda + f7*(cos2w-sin2w)/lambda

         zmw2 = dcmplx(wmass**2, -wmass*wwidth)
         zmz2 = dcmplx(zmass**2, -zmass*zwidth)
         lreset = .false.

       if ((idfsvv.eq.1).or.(idfsvv.eq.4)) then  ! V3V4 = gamma gamma or gamma Z
         lk11 = f1ww.ne.0d0 .or. f1zz.ne.0d0 .or. 
     1          f1az.ne.0d0 .or. f1aa.ne.0d0
         lk21 = f2ww.ne.0d0 .or. f2zz.ne.0d0
         lk12 = .false.
         lk22 = .false.
       elseif ((idfsvv.eq.2).or.(idfsvv.eq.3)) then  ! V3V4 = WW or ZZ
         lk11 = f1ww.ne.0d0 .or. f1zz.ne.0d0 .or. 
     1          f1az.ne.0d0 .or. f1aa.ne.0d0
         lk21 = f2ww.ne.0d0 .or. f2zz.ne.0d0
         lk12 = f2ww.ne.0d0 .or. f1ww.ne.0d0 .or. f1zz.ne.0d0 .or.
     1          f1az.ne.0d0 .or. f1aa.ne.0d0
         lk22 = f2ww.ne.0d0 .or. f2zz.ne.0d0
       else
         print*,"final state not implemented in vvsp2tovv or called with wrong argument idfsvv"
         stop
       endif

      endif

      do mu = 0,3
         k(mu) = p1(mu) + p2(mu)
         p1mp2(mu) = p1(mu) - p2(mu)
         q1mq2(mu) = q1(mu) - q2(mu)
      enddo

      if ((f6.EQ.0.0d0) .and. (f7.EQ.0.0d0)) then      !singlet case
      spin2massgen=sp2mass
      spin2widthgen=sp2width
      mtsq = sp2mass**2
      elseif ((f1.EQ.0.0d0) .and. (f2.EQ.0.0d0) .and. (f3.EQ.0.0d0) .and.
     &    (f4.EQ.0.0d0) .and. (f5.EQ.0.0d0)) then      !triplet case
      spin2massgen=msp2tripn
      spin2widthgen=sp2tripnwidth
      mtsq = msp2tripn**2
      else
         print*,"spin-2 couplings not valid -> choose spin-2 singlet or triplet case and set the other couplings to zero"
         stop
      endif

      sbar = k(0)**2 - k(1)**2 - k(2)**2 - k(3)**2 - mtsq
      p1sq = p1(0)**2 - p1(1)**2 - p1(2)**2 - p1(3)**2
      p2sq = p2(0)**2 - p2(1)**2 - p2(2)**2 - p2(3)**2
      p1mp2sq = p1mp2(0)**2 - p1mp2(1)**2 - p1mp2(2)**2 - p1mp2(3)**2  
      q1mq2sq = q1mq2(0)**2 - q1mq2(1)**2 - q1mq2(2)**2 - q1mq2(3)**2  
      q12p12 = q1mq2(0)*p1mp2(0) - q1mq2(1)*p1mp2(1) - 
     1               q1mq2(2)*p1mp2(2) - q1mq2(3)*p1mp2(3)
      kp1mp2 = k(0)*p1mp2(0) - k(1)*p1mp2(1) - 
     1                k(2)*p1mp2(2) - k(3)*p1mp2(3) 
      kq1mq2 = k(0)*q1mq2(0) - k(1)*q1mq2(1) - 
     1                k(2)*q1mq2(2) - k(3)*q1mq2(3) 

      e1p1mp2 = eps1(0)*p1mp2(0) - eps1(1)*p1mp2(1) - 
     1                   eps1(2)*p1mp2(2) - eps1(3)*p1mp2(3) 
      e2p1mp2 = eps2(0)*p1mp2(0) - eps2(1)*p1mp2(1) - 
     1                   eps2(2)*p1mp2(2) - eps2(3)*p1mp2(3) 
      keps1 = k(0)*eps1(0) - k(1)*eps1(1) - k(2)*eps1(2) - k(3)*eps1(3)
      keps2 = k(0)*eps2(0) - k(1)*eps2(1) - k(2)*eps2(2) - k(3)*eps2(3)
      e1e2 = eps1(0)*eps2(0) - eps1(1)*eps2(1) - eps1(2)*eps2(2) - 
     1                                         eps1(3)*eps2(3) 

C set parameters for K1 * K1 interaction
      if (lk11) then
      ax1= p1mp2sq - mtsq + sbar
      ax4 = 3d0*mtsq*(sbar - mtsq + q1mq2sq)*
     1                          (keps1*e2p1mp2+keps2*e1p1mp2)
      ax5 = q1mq2sq*2d0*mtsq-mtsq**2-2d0*sbar**2
      ax6 =  -2d0*(mtsq+q1mq2sq+3d0*sbar)
      ax7 = mtsq*(9d0*sbar - 5d0*mtsq + q1mq2sq )

      res12 = (sbar + mtsq - p1mp2sq)*(sbar + mtsq - q1mq2sq) 
     1            * 3d0 * mtsq**2

      reskmu1 = 2d0*e1e2*(- 2d0*kp1mp2 + mtsq + 2d0*sbar)* kq1mq2**2 +
     1 6d0*mtsq*(e1e2*q12p12 + keps1*e2p1mp2 - keps2*e1p1mp2) * kq1mq2 +
     2 ax4 + 2d0*e1e2*(ax5*( -kp1mp2 -mtsq +sbar) - 6d0*mtsq*sbar**2) +
     3 2d0*keps1*keps2*( ax6*(kp1mp2-sbar) + ax7 )

      resknu1 = 2d0*e1e2*(- 2d0*kp1mp2 - mtsq - 2d0*sbar)* kq1mq2**2 +
     1 6d0*mtsq*(e1e2*q12p12 + keps1*e2p1mp2 - keps2*e1p1mp2) * kq1mq2 +
     2 ax4 + 2d0*e1e2*(ax5*( -kp1mp2 +mtsq -sbar) + 6d0*mtsq*sbar**2) +
     3 2d0*keps1*keps2*( ax6*(kp1mp2+sbar) - ax7 )

      reskk = 2d0*(2d0*ax1*e1e2*kq1mq2**2-ax1*e1e2*mtsq**2+2d0*ax1*
     . e1e2*mtsq*q1mq2sq-2d0*ax1*e1e2*sbar**2+2d0*ax1*keps1*keps2*
     . mtsq+2d0*ax1*keps1*keps2*q1mq2sq+6d0*ax1*keps1*keps2*sbar-6d0*
     . e1e2*mtsq*sbar**2+12d0*keps1*keps2*mtsq*sbar)

      reskmu2 = 
     1 6d0*mtsq*(-e1e2*kp1mp2*kq1mq2+e1e2*kq1mq2*mtsq-e1e2*kq1mq2*
     2 p1mp2sq+e1e2*mtsq*q12p12-e1p1mp2*keps2*mtsq+e2p1mp2*keps1*mtsq)
 
      reskmu3 = 3d0*mtsq*(-e2p1mp2*mtsq**2+e2p1mp2*mtsq*q1mq2sq-e2p1mp2*
     . mtsq*sbar+2d0*keps2*kp1mp2*kq1mq2+keps2*kp1mp2*mtsq-keps2*
     . kp1mp2*q1mq2sq-keps2*kp1mp2*sbar-2d0*keps2*kq1mq2*mtsq+2d0*
     . keps2*kq1mq2*p1mp2sq-keps2*mtsq**2+keps2*mtsq*p1mp2sq-2d0*
     . keps2*mtsq*q12p12+keps2*mtsq*q1mq2sq-keps2*mtsq*sbar-keps2*
     . p1mp2sq*q1mq2sq-keps2*p1mp2sq*sbar)
 
      reskmu4=3d0*mtsq*(-e1p1mp2*mtsq**2+e1p1mp2*mtsq*q1mq2sq-e1p1mp2*
     . mtsq*sbar-2d0*keps1*kp1mp2*kq1mq2+keps1*kp1mp2*mtsq-keps1*
     . kp1mp2*q1mq2sq-keps1*kp1mp2*sbar+2d0*keps1*kq1mq2*mtsq-2d0*
     . keps1*kq1mq2*p1mp2sq-keps1*mtsq**2+keps1*mtsq*p1mp2sq+2d0*
     . keps1*mtsq*q12p12+keps1*mtsq*q1mq2sq-keps1*mtsq*sbar-keps1*
     . p1mp2sq*q1mq2sq-keps1*p1mp2sq*sbar)
 
      resknu2=6d0*mtsq*(e1e2*kp1mp2*kq1mq2+e1e2*kq1mq2*mtsq-e1e2*kq1mq2*
     . p1mp2sq-e1e2*mtsq*q12p12+e1p1mp2*keps2*mtsq-e2p1mp2*keps1*mtsq)

      resknu3=3d0*mtsq*(e2p1mp2*mtsq**2-e2p1mp2*mtsq*q1mq2sq+e2p1mp2*
     . mtsq*sbar-2d0*keps2*kp1mp2*kq1mq2-keps2*kp1mp2*mtsq+keps2*kp1mp2*
     . q1mq2sq+keps2*kp1mp2*sbar-2d0*keps2*kq1mq2*mtsq+2d0*keps2*
     . kq1mq2*p1mp2sq-keps2*mtsq**2+keps2*mtsq*p1mp2sq+2d0*keps2*mtsq
     . *q12p12+keps2*mtsq*q1mq2sq-keps2*mtsq*sbar-keps2*p1mp2sq*
     . q1mq2sq-keps2*p1mp2sq*sbar)
 
      resknu4=3d0*mtsq*(e1p1mp2*mtsq**2-e1p1mp2*mtsq*q1mq2sq+e1p1mp2*
     . mtsq*sbar+2d0*keps1*kp1mp2*kq1mq2-keps1*kp1mp2*mtsq+keps1*kp1mp2*
     . q1mq2sq+keps1*kp1mp2*sbar+2d0*keps1*kq1mq2*mtsq-2d0*keps1*
     . kq1mq2*p1mp2sq-keps1*mtsq**2+keps1*mtsq*p1mp2sq-2d0*keps1*mtsq
     . *q12p12+keps1*mtsq*q1mq2sq-keps1*mtsq*sbar-keps1*p1mp2sq*
     . q1mq2sq-keps1*p1mp2sq*sbar)
 
      reseps1mua=3d0*mtsq*(e2p1mp2*mtsq**2-e2p1mp2*mtsq*q1mq2sq+
     .e2p1mp2*mtsq
     . *sbar-2d0*keps2*kp1mp2*kq1mq2-keps2*kp1mp2*mtsq+keps2*kp1mp2*
     . q1mq2sq+keps2*kp1mp2*sbar+2d0*keps2*kq1mq2*sbar+2d0*keps2*mtsq
     . *q12p12-keps2*mtsq*sbar-keps2*q1mq2sq*sbar-keps2*sbar**2)
 
      reseps1mub=6d0*keps2*mtsq**2*(mtsq-p1mp2sq+sbar)
 
      reseps1nua=3d0*mtsq*(e2p1mp2*mtsq**2-e2p1mp2*mtsq*q1mq2sq+e2p1mp2*
     . mtsq*sbar-2d0*keps2*kp1mp2*kq1mq2-keps2*kp1mp2*mtsq+keps2*kp1mp2*
     . q1mq2sq+keps2*kp1mp2*sbar-2d0*keps2*kq1mq2*sbar+2d0*keps2*mtsq
     . *q12p12+keps2*mtsq*sbar+keps2*q1mq2sq*sbar+keps2*sbar**2)
 
      reseps1nub=6d0*keps2*mtsq**2*(mtsq-p1mp2sq+sbar)
  
      reseps2mua=3d0*mtsq*(e1p1mp2*mtsq**2-e1p1mp2*mtsq*q1mq2sq+e1p1mp2*
     . mtsq*sbar+2d0*keps1*kp1mp2*kq1mq2-keps1*kp1mp2*mtsq+keps1*kp1mp2*
     . q1mq2sq+keps1*kp1mp2*sbar-2d0*keps1*kq1mq2*sbar-2d0*keps1*mtsq
     . *q12p12-keps1*mtsq*sbar-keps1*q1mq2sq*sbar-keps1*sbar**2)
 
      reseps2mub=6d0*keps1*mtsq**2*(-mtsq+p1mp2sq-sbar)
 
      reseps2nua=3d0*mtsq*(e1p1mp2*mtsq**2-e1p1mp2*mtsq*q1mq2sq+e1p1mp2*
     . mtsq*sbar+2d0*keps1*kp1mp2*kq1mq2-keps1*kp1mp2*mtsq+keps1*kp1mp2*
     . q1mq2sq+keps1*kp1mp2*sbar+2d0*keps1*kq1mq2*sbar-2d0*keps1*mtsq
     . *q12p12+keps1*mtsq*sbar+keps1*q1mq2sq*sbar+keps1*sbar**2)
 
      reseps2nub=6d0*keps1*mtsq**2*(-mtsq+p1mp2sq-sbar)

      respp=4d0*mtsq*(-e1e2*kq1mq2**2-e1e2*mtsq**2+2d0*e1e2*mtsq*
     . q1mq2sq+e1e2*sbar**2+5d0*keps1*keps2*mtsq-keps1*keps2*q1mq2sq-
     . 3d0*keps1*keps2*sbar)
 
      respq=6d0*mtsq*(e1e2*kp1mp2*kq1mq2+e1e2*kq1mq2*sbar-e1e2*mtsq*
     . q12p12+e1p1mp2*keps2*mtsq-e2p1mp2*keps1*mtsq)
 
      resqp=6d0*mtsq*(e1e2*kp1mp2*kq1mq2-e1e2*kq1mq2*sbar-e1e2*mtsq*
     . q12p12+e1p1mp2*keps2*mtsq-e2p1mp2*keps1*mtsq)
 
      resqq=6d0*e1e2*mtsq**2*(-mtsq+p1mp2sq-sbar)

      resmunu=-(reskk*(-kp1mp2+mtsq+sbar)+reskmu1*(kp1mp2-p1mp2sq)+
     . reskmu2*(kq1mq2-q12p12)+reskmu3*(-e1p1mp2+keps1)+
     . reskmu4*(-e2p1mp2+keps2))

      do mu = 0,3
         do nu = 0,3
            t11(mu,nu) = reskk*k(mu)*k(nu) + k(mu)*(reskmu1*p1mp2(nu)
     1           +reskmu2*q1mq2(nu)+reskmu3*eps1(nu)+reskmu4*eps2(nu)) +
     2  k(nu)*(resknu1*p1mp2(mu)+resknu2*q1mq2(mu)+resknu3*eps1(mu)+
     3             resknu4*eps2(mu)) +
     4             res12 * (eps1(mu)*eps2(nu) + eps1(nu)*eps2(mu)) +
     5  eps1(mu)*(reseps1mua*p1mp2(nu) + reseps1mub*q1mq2(nu)) +
     6  eps1(nu)*(reseps1nua*p1mp2(mu) + reseps1nub*q1mq2(mu)) +
     7  eps2(mu)*(reseps2mua*p1mp2(nu) + reseps2mub*q1mq2(nu)) +
     8  eps2(nu)*(reseps2nua*p1mp2(mu) + reseps2nub*q1mq2(mu)) +
     9  respp*p1mp2(mu)*p1mp2(nu) + respq*p1mp2(mu)*q1mq2(nu) +
     1          resqp*q1mq2(mu)*p1mp2(nu) + resqq*q1mq2(mu)*q1mq2(nu)
         enddo
      enddo
      t11(0,0) = t11(0,0) + resmunu
      do i = 1,3
         t11(i,i) = t11(i,i) - resmunu
      enddo
      else
         do mu = 0,3
            do nu = 0,3 
               t11(mu,nu) = 0d0
            enddo
         enddo
      endif

      if (lk21) then
C prepare K2*K1 interaction
         resmunu=2d0*mtsq*(-e1e2*kq1mq2**2-e1e2*mtsq**2+2d0*e1e2*mtsq*
     . q1mq2sq+e1e2*sbar**2+5d0*keps1*keps2*mtsq-keps1*keps2*q1mq2sq-
     . 3d0*keps1*keps2*sbar)

         reskk=
     .  2d0*(-2d0*e1e2*kq1mq2**2+e1e2*mtsq**2-2d0*e1e2*mtsq*q1mq2sq
     . +2d0*e1e2*sbar**2-2d0*keps1*keps2*mtsq-2d0*keps1*keps2*q1mq2sq
     . -6d0*keps1*keps2*sbar)

         res12=3d0*mtsq**2*(mtsq-q1mq2sq+sbar)

         reskmu1=0d0
 
         reskmu2=6d0*e1e2*kq1mq2*mtsq
 
         reskmu3=3d0*keps2*mtsq*(-2d0*kq1mq2-mtsq+q1mq2sq+sbar)
 
         reskmu4=3d0*keps1*mtsq*(2d0*kq1mq2-mtsq+q1mq2sq+sbar)
 
         resknu1=0d0
 
         resknu2=6d0*e1e2*kq1mq2*mtsq

         resknu3=3d0*keps2*mtsq*(-2d0*kq1mq2-mtsq+q1mq2sq+sbar)
 
         resknu4=3d0*keps1*mtsq*(2d0*kq1mq2-mtsq+q1mq2sq+sbar)
 
         reseps1mua=0d0
 
         reseps1mub=6d0*keps2*mtsq**2
 
         reseps1nua=0d0
 
         reseps1nub=6d0*keps2*mtsq**2
  
         reseps2mua=0d0
 
         reseps2mub=-6d0*keps1*mtsq**2
 
         reseps2nua=0d0
 
         reseps2nub=-6d0*keps1*mtsq**2

         respp=0d0
         respq=0d0
         resqp=0d0
         resqq=-6d0*e1e2*mtsq**2

         do mu = 0,3
            do nu = 0,3
            t21(mu,nu) = reskk*k(mu)*k(nu) + k(mu)*(reskmu1*p1mp2(nu)
     1           +reskmu2*q1mq2(nu)+reskmu3*eps1(nu)+reskmu4*eps2(nu)) +
     2  k(nu)*(resknu1*p1mp2(mu)+resknu2*q1mq2(mu)+resknu3*eps1(mu)+
     3             resknu4*eps2(mu)) +
     4             res12 * (eps1(mu)*eps2(nu) + eps1(nu)*eps2(mu)) +
     5  eps1(mu)*(reseps1mua*p1mp2(nu) + reseps1mub*q1mq2(nu)) +
     6  eps1(nu)*(reseps1nua*p1mp2(mu) + reseps1nub*q1mq2(mu)) +
     7  eps2(mu)*(reseps2mua*p1mp2(nu) + reseps2mub*q1mq2(nu)) +
     8  eps2(nu)*(reseps2nua*p1mp2(mu) + reseps2nub*q1mq2(mu)) +
     9  respp*p1mp2(mu)*p1mp2(nu) + respq*p1mp2(mu)*q1mq2(nu) +
     1          resqp*q1mq2(mu)*p1mp2(nu) + resqq*q1mq2(mu)*q1mq2(nu)
            enddo
         enddo
         t21(0,0) = t21(0,0) + resmunu
         do i = 1,3
            t21(i,i) = t21(i,i) - resmunu
         enddo
      else
         do mu = 0,3
            do nu = 0,3 
               t21(mu,nu) = 0d0
            enddo
         enddo
      endif

      if (lk12) then
C prepare K1*K2 interaction
         resmunu=2d0*e1e2*mtsq*(-kp1mp2**2-mtsq**2+2d0*mtsq*p1mp2sq+sbar**2)
     . +2d0*keps1*keps2*(-2d0*kp1mp2**2+mtsq**2-2d0*mtsq*p1mp2sq+2d0*
     . sbar**2)+6d0*keps1*e2p1mp2*kp1mp2*mtsq+6d0*keps2*e1p1mp2*
     . kp1mp2*mtsq-6d0*e1p1mp2*e2p1mp2*mtsq**2

         reskk=2d0*e1e2*mtsq*(mtsq-p1mp2sq-sbar)+4d0*keps1*keps2*(mtsq-
     . p1mp2sq-sbar)

         res12=3d0*mtsq**2*(mtsq-p1mp2sq+sbar)

         reskmu1=2d0*e1e2*mtsq*(kp1mp2+mtsq-sbar)+2d0*keps1*keps2*(2d0*
     . kp1mp2-mtsq-2d0*sbar)-3d0*keps1*e2p1mp2*mtsq-
     . 3d0*keps2*e1p1mp2*mtsq
 
         reskmu2=0d0
 
        reskmu3=3d0*keps2*mtsq*(kp1mp2-mtsq+p1mp2sq)-3d0*e2p1mp2*mtsq**2
 
        reskmu4=3d0*keps1*mtsq*(kp1mp2-mtsq+p1mp2sq)-3d0*e1p1mp2*mtsq**2
 
         resknu1=2d0*e1e2*mtsq*(kp1mp2-mtsq+sbar)+2d0*keps1*keps2*(2d0*
     . kp1mp2+mtsq+2d0*sbar)-3d0*keps1*e2p1mp2*mtsq-3d0*keps2*e1p1mp2
     . *mtsq
 
         resknu2=0d0

       resknu3=3d0*keps2*mtsq*(-kp1mp2-mtsq+p1mp2sq)+3d0*e2p1mp2*mtsq**2
 
       resknu4=3d0*keps1*mtsq*(-kp1mp2-mtsq+p1mp2sq)+3d0*e1p1mp2*mtsq**2

         reseps1mua=3d0*keps2*mtsq*(-kp1mp2+sbar)+3d0*e2p1mp2*mtsq**2
 
         reseps1mub = 0d0

         reseps1nua=-3d0*keps2*mtsq*(kp1mp2+sbar)+3d0*e2p1mp2*mtsq**2
 
         reseps1nub = 0d0
  
         reseps2mua = 3d0*keps1*mtsq*(-kp1mp2+sbar)+3d0*e1p1mp2*mtsq**2
 
         reseps2mub =0d0
 
         reseps2nua =-3d0*keps1*mtsq*(kp1mp2+sbar)+3d0*e1p1mp2*mtsq**2
 
         reseps2nub=0d0

         respp=-4d0*e1e2*mtsq**2+4d0*keps1*keps2*mtsq
 
         respq=0d0
         resqp=0d0
         resqq=0d0

         do mu = 0,3
            do nu = 0,3
            t12(mu,nu) = reskk*k(mu)*k(nu) + k(mu)*(reskmu1*p1mp2(nu)
     1           +reskmu3*eps1(nu)+reskmu4*eps2(nu)) +
     2  k(nu)*(resknu1*p1mp2(mu)+resknu3*eps1(mu)+resknu4*eps2(mu)) +
     3             res12 * (eps1(mu)*eps2(nu) + eps1(nu)*eps2(mu)) +
     4  eps1(mu)*reseps1mua*p1mp2(nu) + eps1(nu)*reseps1nua*p1mp2(mu) +
     5  eps2(mu)*reseps2mua*p1mp2(nu) + eps2(nu)*reseps2nua*p1mp2(mu) +
     6  respp*p1mp2(mu)*p1mp2(nu) 
            enddo
         enddo
         t12(0,0) = t12(0,0) + resmunu
         do i = 1,3
            t12(i,i) = t12(i,i) - resmunu
         enddo
      else
         do mu = 0,3
            do nu = 0,3 
               t12(mu,nu) = 0d0
            enddo
         enddo
      endif

      if (lk22) then
C prepare K2*K2 interaction
         resmunu=-2d0*e1e2*mtsq**2+2d0*keps1*keps2*mtsq

         reskk=2d0*e1e2*mtsq+4d0*keps1*keps2

         res12=3d0*mtsq**2

         reskmu3=-3d0*keps2*mtsq
 
         reskmu4=-3d0*keps1*mtsq
 
         resknu3=-3d0*keps2*mtsq
 
         resknu4=-3d0*keps1*mtsq
 

         do mu = 0,3
            do nu = 0,3
               t22(mu,nu) = reskk*k(mu)*k(nu) + 
     1              k(mu)*(reskmu3*eps1(nu)+reskmu4*eps2(nu)) +
     2              k(nu)*(resknu3*eps1(mu)+resknu4*eps2(mu)) +
     4              res12 * (eps1(mu)*eps2(nu) + eps1(nu)*eps2(mu)) 
            enddo
         enddo
         t22(0,0) = t22(0,0) + resmunu
         do i = 1,3
            t22(i,i) = t22(i,i) - resmunu
         enddo
      else
         do mu = 0,3
            do nu = 0,3 
               t22(mu,nu) = 0d0
            enddo
         enddo
      endif

      formfac=(lambdaff**2/(abs(sbar+mtsq)+lambdaff**2)*
     &         lambdaff**2/(abs(p1sq)+lambdaff**2)*
     &         lambdaff**2/(abs(p2sq)+lambdaff**2))**nff


      do mu = 0,3
         do nu = 0,3
            t11(mu,nu) = t11(mu,nu)*formfac/( 96d0*mtsq**2 
     1           * dcmplx(sbar,spin2massgen*spin2widthgen) )
            t21(mu,nu) = t21(mu,nu)*formfac/( 24d0*mtsq**2 
     1           * dcmplx(sbar,spin2massgen*spin2widthgen) )
            t12(mu,nu) = t12(mu,nu)*formfac/( 24d0*mtsq**2 
     1           * dcmplx(sbar,spin2massgen*spin2widthgen) )
            t22(mu,nu) = t22(mu,nu)*formfac/( 6d0*mtsq**2 
     1           * dcmplx(sbar,spin2massgen*spin2widthgen) )
         enddo
      enddo


      do mu = 0,3
       do nu = 0,3

        if (idfsvv.eq.1) then  ! V3V4 = gamma gamma

            wwsp2vv(mu,nu) = 
     1          (t11(mu,nu)*f1ww*f1aa + t21(mu,nu)*f2ww*f1aa) /
     2          ( (p1sq - zmw2) * (p2sq - zmw2) )
            wwsp2vvx(mu,nu) =
     1          (t11(nu,mu)*f1ww*f1aa + t21(nu,mu)*f2ww*f1aa) /
     2          ( (p1sq - zmw2) * (p2sq - zmw2) )
            zzsp2vv(mu,nu) = 
     1          (t11(mu,nu)*f1zz*f1aa + t21(mu,nu)*f2zz*f1aa) /
     2          ( (p1sq - zmz2) * (p2sq - zmz2) )
            azsp2vv(mu,nu) = t11(mu,nu)*f1az*f1aa/
     1          ( p1sq  * (p2sq - zmz2) )
            zasp2vv(mu,nu) = t11(nu,mu)*f1az*f1aa/
     1          ( p2sq  * (p1sq - zmz2) )
            aasp2vv(mu,nu) = t11(mu,nu)*f1aa*f1aa/
     1          ( p1sq  * p2sq )

        elseif (idfsvv.eq.2) then  ! V3V4 = WW

            wwsp2vv(mu,nu) = 
     1          (t11(mu,nu)*f1ww*f1ww + (t21(mu,nu)+t12(mu,nu))*f2ww*f1ww
     2            + t22(mu,nu)*f2ww*f2ww) /
     3          ( (p1sq - zmw2) * (p2sq - zmw2) )
            wwsp2vvx(mu,nu) =
     1          (t11(nu,mu)*f1ww*f1ww + (t21(nu,mu)+t12(nu,mu))*f2ww*f1ww
     2            + t22(nu,mu)*f2ww*f2ww) /
     3          ( (p1sq - zmw2) * (p2sq - zmw2) )
            zzsp2vv(mu,nu) = 
     1          (t11(mu,nu)*f1zz*f1ww + t21(mu,nu)*f2zz*f1ww + 
     2           t12(mu,nu)*f1zz*f2ww + t22(mu,nu)*f2zz*f2ww) /
     3          ( (p1sq - zmz2) * (p2sq - zmz2) )
            azsp2vv(mu,nu) = (t11(mu,nu)*f1az*f1ww + t12(mu,nu)*f1az*f2ww) /
     1          ( p1sq  * (p2sq - zmz2) )
            zasp2vv(mu,nu) = (t11(nu,mu)*f1az*f1ww + t12(nu,mu)*f1az*f2ww) /
     1          ( p2sq  * (p1sq - zmz2) )
            aasp2vv(mu,nu) = (t11(mu,nu)*f1aa*f1ww + t12(mu,nu)*f1aa*f2ww) /
     1          ( p1sq  * p2sq )

        elseif (idfsvv.eq.3) then  ! V3V4 = ZZ

            wwsp2vv(mu,nu) = 
     1          (t11(mu,nu)*f1ww*f1zz + t12(mu,nu)*f2zz*f1ww
     2           + t21(mu,nu)*f2ww*f1zz + t22(mu,nu)*f2ww*f2zz) /
     3          ( (p1sq - zmw2) * (p2sq - zmw2) )
            wwsp2vvx(mu,nu) =
     1          (t11(nu,mu)*f1ww*f1zz + t12(nu,mu)*f2zz*f1ww
     2           + t21(nu,mu)*f2ww*f1zz + t22(nu,mu)*f2ww*f2zz) /
     3          ( (p1sq - zmw2) * (p2sq - zmw2) )
            zzsp2vv(mu,nu) = 
     1          (t11(mu,nu)*f1zz*f1zz + t21(mu,nu)*f2zz*f1zz + 
     2           t12(mu,nu)*f1zz*f2zz + t22(mu,nu)*f2zz*f2zz) /
     3          ( (p1sq - zmz2) * (p2sq - zmz2) )
            azsp2vv(mu,nu) = (t11(mu,nu)*f1az*f1zz + t12(mu,nu)*f1az*f2zz) /
     1          ( p1sq  * (p2sq - zmz2) )
            zasp2vv(mu,nu) = (t11(nu,mu)*f1az*f1zz + t12(nu,mu)*f1az*f2zz) /
     1          ( p2sq  * (p1sq - zmz2) )
            aasp2vv(mu,nu) = (t11(mu,nu)*f1aa*f1zz + t12(mu,nu)*f1aa*f2zz) /
     1          ( p1sq  * p2sq )

        elseif (idfsvv.eq.4) then  ! V3V4 = gamma Z

            wwsp2vv(mu,nu) = 
     1          (t11(mu,nu)*f1ww*f1az + t21(mu,nu)*f2ww*f1az) /
     2          ( (p1sq - zmw2) * (p2sq - zmw2) )
            wwsp2vvx(mu,nu) =
     1          (t11(nu,mu)*f1ww*f1az + t21(nu,mu)*f2ww*f1az) /
     2          ( (p1sq - zmw2) * (p2sq - zmw2) )
            zzsp2vv(mu,nu) = 
     1          (t11(mu,nu)*f1zz*f1az + t21(mu,nu)*f2zz*f1az) /
     2          ( (p1sq - zmz2) * (p2sq - zmz2) )
            azsp2vv(mu,nu) = t11(mu,nu)*f1az*f1az/
     1          ( p1sq  * (p2sq - zmz2) )
            zasp2vv(mu,nu) = t11(nu,mu)*f1az*f1az/
     1          ( p2sq  * (p1sq - zmz2) )
            aasp2vv(mu,nu) = t11(mu,nu)*f1aa*f1az/
     1          ( p1sq  * p2sq )

        endif

       enddo
      enddo

!          do mu=0,3
!           do nu=0,3
!         write(*,*) "mu", mu
!         write(*,*) "nu", nu
!         write(*,*) "wwsp2aa(mu,nu)", wwsp2vv(mu,nu)
!         write(*,*) "wwsp2aax(mu,nu)", wwsp2vvx(mu,nu)
!         write(*,*) "zzsp2aa(mu,nu)", zzsp2vv(mu,nu)
!         write(*,*) "azsp2aa(mu,nu)", azsp2vv(mu,nu)
!         write(*,*) "zasp2aa(mu,nu)", zasp2vv(mu,nu)
!         write(*,*) "aasp2aa(mu,nu)", aasp2vv(mu,nu)
!         pause
!           enddo
!          enddo

      end 
