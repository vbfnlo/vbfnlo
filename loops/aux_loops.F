
c
c------------- C(m,q1^2,q2^2,P^2)  3-point fuction ---------------------------
c
c      complex*16 function I3point(m,q1sq,q2sq,Psq)
      complex*16 function I3point(mt,q1sqt,q2sqt,Psqt)
      implicit none
      double precision m, q1sq, q2sq, psq 
      double precision mt,q1sqt,q2sqt,psqt,eps1

c evaluate scalar 3-point function for equal masses m on propagators 
c  
c  I3 = 1/(i*pi^2) * Int d^4k [k^2-m^2]^-1 [(k+q1)^2-m^2]^-1 [(k-q2)^2-m^2]^-1
c
c     = i/pi^2 C("t Hooft,Veltman) = -C0(Passarino,Veltman)
c
c in terms of 12 dilogarithms. P = q1+q2 is assumed to be space-like. q1 
c and q2 may be space-, time- or lightlike
c
c	Dieter Zeppenfeld, <dieter@pheno.physics.wisc.edu>
c	Initial version:  2000 November 12
c	Last modified: Michael Kubocz
c       date         : 20.01.2011

      double precision lambda, rtlam, q1q2, q1p, q2p, fmsq
      double precision rp(3), rm(3), eps, precise, pi

      complex*16 sp(3), sm(3), fmt2, ieps, fac, z, di, vli2
c      external vli2
      parameter (eps=1d-18, ieps = eps*(0,1), precise=1d-3 )
      parameter (pi=3.14159 26535 89793238d0, fac=-1d0/2d0)

      eps1=1d-7
      if(abs(mt).le.eps1) then
         m=0d0
      else
         m=mt
      endif
      if(abs(q1sqt).le.eps1) then
         q1sq=0d0
      else
         q1sq=q1sqt
      endif
      if(abs(q2sqt).le.eps1) then
         q2sq=0d0
      else
         q2sq=q2sqt
      endif
      if(abs(Psqt).le.eps1) then
         Psq=0d0
      else
         Psq=Psqt
      endif

c determine dot products of input momenta
      q1q2 = 0.5d0*(psq - q1sq - q2sq)
      q1p  = 0.5d0*(psq + q1sq - q2sq)
      q2p  = 0.5d0*(psq - q1sq + q2sq)
      lambda = q1q2**2-q1sq*q2sq
      
      if (lambda.le.0d0) then
         write(*,*) " singular lambda in 3-point function. Reset to 0 "
         I3point = 0
         return
      else
         rtlam = dsqrt(lambda)
      endif
      
c determine factors for call of spence functions
      fmsq = 4*m**2
      fmt2 = fmsq - ieps
      di = q1sq*q2sq*psq + fmt2*lambda 
      if ( abs(dreal(di)).gt.eps )then
         di = 1d0/di
      else
         print*," singular point called in I3point "
         stop
      endif

      rp(1) = q2p + rtlam
      rm(1) = q2p - rtlam
      if (abs(rm(1)).lt.abs(rp(1))*precise ) then
         rm(1) = q2sq*psq/rp(1)
      elseif (abs(rp(1)).lt.abs(rm(1))*precise ) then
         rp(1) = q2sq*psq/rm(1)
      endif

      rp(2) = q1q2 + rtlam
      rm(2) = q1q2 - rtlam
      if (abs(rm(2)).lt.abs(rp(2))*precise ) then
         rm(2) = q1sq*q2sq/rp(2)
      elseif (abs(rp(2)).lt.abs(rm(2))*precise ) then
         rp(2) = q1sq*q2sq/rm(2)
      endif

      rp(3) = q1p + rtlam
      rm(3) = q1p - rtlam
      if (abs(rm(3)).lt.abs(rp(3))*precise ) then
         rm(3) = q1sq*psq/rp(3)
      elseif (abs(rp(3)).lt.abs(rm(3))*precise ) then
         rp(3) = q1sq*psq/rm(3)
      endif

      if (q1sq.ne.0d0) then
         z = rtlam*sqrt((q1sq-fmt2)*q1sq) 
         sp(1) = di*( q2p*q1sq + z )
         sm(1) = di*( q2p*q1sq - z )
      else
         sp(1) = 0
         sm(1) = 0
      endif

      if (psq.ne.0d0) then
         z = rtlam*sqrt((psq-fmt2)*psq) 
         sp(2) = di*( q1q2*psq + z )
         sm(2) = di*( q1q2*psq - z )
      else
         sp(2) = 0
         sm(2) = 0
      endif

      if (q2sq.ne.0d0) then
         z = rtlam*sqrt((q2sq-fmt2)*q2sq) 
         sp(3) = di*( q1p*q2sq + z )
         sm(3) = di*( q1p*q2sq - z )
      else
         sp(3) = 0
         sm(3) = 0
      endif

      z = 0
      if (q1sq.lt.0d0 .or. q1sq.gt.fmsq) then
         z = z + vli2(rm(1)*sp(1)) + vli2(rm(1)*sm(1)) - 
     &           vli2(rp(1)*sm(1)) - vli2(rp(1)*sp(1))
      elseif (q1sq.gt.0d0) then
         z = z + 2*dreal( vli2(rm(1)*sp(1)) ) - 
     &           2*dreal( vli2(rp(1)*sp(1)) )
      endif
      if (psq.lt.0d0 .or. psq.gt.fmsq) then
         z = z - vli2(rm(2)*sp(2)) - vli2(rm(2)*sm(2)) + 
     &           vli2(rp(2)*sm(2)) + vli2(rp(2)*sp(2))
      elseif (psq.gt.0d0) then
         z = z - 2*dreal( vli2(rm(2)*sp(2)) ) + 
     &           2*dreal( vli2(rp(2)*sp(2)) )
      endif
      if (q2sq.lt.0d0 .or. q2sq.gt.fmsq) then
         z = z + vli2(rm(3)*sp(3)) + vli2(rm(3)*sm(3)) - 
     &           vli2(rp(3)*sm(3)) - vli2(rp(3)*sp(3))
      elseif (q2sq.gt.0d0) then
         z = z + 2*dreal( vli2(rm(3)*sp(3)) ) - 
     &           2*dreal( vli2(rp(3)*sp(3)) )
      endif

      I3point = z*fac/rtlam
      return
      end


********************************************************************************
********************************************************************************


c
c-------------  spence function  vli2(z) ------------------------
c
      complex*16 function vli2(zin)
      implicit none
      complex*16 zin, z, u, u2, unpo, ans, zext
      double precision r, r2, r2n, fac
c
c determine the value of the dilogarithm 
c
c    vli2(z) = - int_0^1  log(1-zt)/t dt  with cut along the positive 
c                                        real axis, z>1
c
c	Dieter Zeppenfeld, <dieter@pheno.physics.wisc.edu>
c	Initial version:  2000 November 6
c	Last modified:    2000 November 12
c  
      integer i
      double precision c0,c1,c2,c4,c6,c8,c10,c12,c14,c16,c18,c20,c22
      double precision b0,b1,b2,b4,b6,b8,b10,b12,b14,b16,b18,b20,b22
      double precision d0,d1,d2,d4,d6,d8,d10,d12,d14,d16,d18,d20,d22
      parameter (b0=1d0,            d0 =1d0,      c0= b0/d0)
      parameter (b1=-1d0/2d0,       d1 =d0*2d0,   c1= b1/d1)
      parameter (b2= 1d0/6d0,       d2 =d1*3d0,   c2= b2/d2)
      parameter (b4=-1d0/30d0,      d4 =d2*20d0,  c4= b4/d4)
      parameter (b6=1d0/42d0,       d6 =d4*42d0,  c6= b6/d6)
      parameter (b8=-1d0/30d0,      d8 =d6*72d0,  c8= b8/d8)
      parameter (b10=5d0/66d0,      d10=d8*110d0, c10=b10/d10)
      parameter (b12=-691d0/2730d0, d12=d10*156d0,c12=b12/d12)
      parameter (b14=7d0/6d0,       d14=d12*210d0,c14=b14/d14)
      parameter (b16=-3617d0/510d0, d16=d14*272d0,c16=b16/d16)
      parameter (b18=43867d0/798d0, d18=d16*342d0,c18=b18/d18)
      parameter (b20=-174611d0/330d0,d20=d18*420d0,c20=b20/d20)
      parameter (b22=854513d0/138d0,d22=d20*506d0,c22=b22/d22)
      double precision eps, epst, pi, pi2o6
      parameter (eps=1d-16, epst=1d-3)
      parameter (pi=3.14159 26535 89793238d0, pi2o6=pi**2/6d0)
c
c debug information
      logical ldebug
      parameter (ldebug=.false.)

      z = zin
c      print*," vli2 call with z = ",z
      u = z**2
      r2 = dreal(z)**2+dimag(z)**2 
      if (r2.lt.eps) then
         vli2 = z + u/4d0
         return
      elseif (r2.lt.epst) then
         ans = z + u/4d0
         do i = 3,11
            u = u*z
            ans = ans + u/i**2
         enddo
         vli2 = ans
         return
      endif
      if (dreal(z).ge.1d0 .and. dimag(z).eq.0 ) then
         z = z + (0d0,1d0)*eps
      endif
c
c use z-->1/z and z--> 1-z mappings of the spence function to restrict 
c agument to unit circle in the complex plane with Re(z) <= 0.5
c
      zext = (0d0,0d0)
      fac = 1
      if (r2.gt.1d0) then     ! map z ---> 1/z
         fac = -fac
         zext = -pi2o6 - 0.5d0*(log(-z))**2
         z = 1d0/z
      endif
      if (dreal(z).gt.0.5d0) then     ! map new z ---> 1-z
         zext = zext + fac*(pi2o6-log(z)*log(1-z))
         fac = -fac
         z = 1-z
      endif
c
c now use t = 1 - exp(-u) mapping to write Li(z) in terms of Bernoulli 
c numbers
c
      u = - log(1-z)
      r2 = abs(u)**2
      u2 = u*u
      ans = u*(c0 + u*(c1+c2*u))
      r2n = r2*r2       !r^4

      unpo = u2*u2*u
      ans = ans + c4*unpo

      unpo = unpo*u2
      ans = ans + c6*unpo

      r = r2n*r2n       !r^8
      unpo = unpo*u2
      ans = ans + c8*unpo

      r2n = r*r2        !r^10 
      if ((r2n*c10).gt.eps) then
         unpo = unpo*u2
         ans = ans + c10*unpo
      else
         vli2 = fac * ans + zext
         if (ldebug) print*," exit li2s at n=8 "
         return
      endif

      unpo = unpo*u2
      ans = ans + c12*unpo

      unpo = unpo*u2
      ans = ans + c14*unpo

      unpo = unpo*u2
      ans = ans + c16*unpo

      r2n = r2n*r
      if ((r2n*c18).gt.eps) then
         unpo = unpo*u2
         ans = ans + c18*unpo
      else
         vli2 = fac * ans + zext
         if (ldebug) print*," exit li2s at n=16 "
         return
      endif

      unpo = unpo*u2
      ans = ans + c20*unpo

      unpo = unpo*u2
      ans = ans + c22*unpo

      vli2 = fac * ans + zext
      end
