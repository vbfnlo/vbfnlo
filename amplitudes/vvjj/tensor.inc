c Dieters version from Apr. 12, 2005, extended for real emission contributions
c
c  precalculated lepton tensors forr WW --> l nubar lbar nu
c
c  lepton spinors and W+- polarization vectors
      double complex wep, wve, wmu, wvm, wp, wm, wpp, wmp, xp, xm
      double precision qp, qm
      common /wpmspinor/ wep(6), wve(6), wmu(6), wvm(6),
     1                   wp(6), wm(6),
     2                   wpp(6), wmp(6), xp, xm, qp(0:4), qm(0:4)
      double complex wel, wvu, we, wu
      double complex xxew,xxuw,wet,wut
      double precision qww
      common /wpmspinor/ wel(6), wvu(6),
     1                   we(6), wu(6), qww(0:4),
     2			 wet(6),wut(6), xxew,xxuw

c  leptonic tensors
      double complex aww, zww, aaww, azww, zaww, zzww, 
     1           NCwpa, NCwpz, CCwpa, CCwpz, wwww5, 
     2           NCwma, NCwmz, CCwma, CCwmz, wwww6
      logical lfs
      common /tensww/ aww(0:5), zww(0:5), aaww(0:3,0:3,3), 
     1                azww(0:3,0:3,3), zaww(0:3,0:3,3), 
     2                zzww(0:3,0:3,3), wwww5(0:3,0:3,3), 
     3                wwww6(0:3,0:3,3),
     4                NCwpa(0:3,0:3,2,3), NCwpz(0:3,0:3,2,3), 
     5                CCwpa(0:3,0:3,2,3), CCwpz(0:3,0:3,2,3),
     6                NCwma(0:3,0:3,2,3), NCwmz(0:3,0:3,2,3), 
     7                CCwma(0:3,0:3,2,3), CCwmz(0:3,0:3,2,3), lfs(8)
c  lfs remembers whether or not an fsign combination needs to be calculated 
c     lfs(1) = .false.      fsign(1)=fsign(2) = -1 has been calcd and stored
c     lfs(2) = .false.      fsign(1)=fsign(2) = +1 has been calcd and stored
c     lfs(3) = .false.      fsign(3)=fsign(4) = -1 has been calcd and stored
c     lfs(4) = .false.      fsign(3)=fsign(4) = +1 has been calcd and stored
c lfs(1,+1) to lfs(4,+1) are for gsign=+1 (outgoing gluon) in case of real emission
c for incoming gluons use lfs(1,-1) to lfs(4,-1):
c     lfs(5) = .false.      fsign(1) = -1 fsign(2) = +1 has been calcd ...
c     lfs(6) = .false.      fsign(1) = +1 fsign(2) = -1 has been calcd ...
c     lfs(7) = .false.      fsign(3) = -1 fsign(4) = +1 has been calcd ...
c     lfs(8) = .false.      fsign(3) = +1 fsign(4) = -1 has been calcd ...
c for testing purposes
      double complex tamp,tampr
      integer ip,ipr,is1,is3
      common /wwtest/ tamp(200,-1:1,-1:1,6), ip
      common /wwrtest/ tampr(900,-1:1,-1:1,-1:1,6), ipr


c  leptonic tensors
      double complex wppwpp,wmmwmm, 
     #		       awen,zwen,waen,wzen,
     #		       awmu,zwmu,wamu,wzmu     
      logical lfp,lfm
      common /tenswpm/ wppwpp(0:3,0:3,3),wmmwmm(0:3,0:3,3),
     #		       awen(0:3,0:3,2,3),zwen(0:3,0:3,2,3),
     #		       waen(0:3,0:3,2,3),wzen(0:3,0:3,2,3),
     #		       awmu(0:3,0:3,2,3),zwmu(0:3,0:3,2,3),
     #		       wamu(0:3,0:3,2,3),wzmu(0:3,0:3,2,3),	
     #		       lfp(4),lfm(4)

c for testing purposes:
      double complex tamwp,tamwpr
      integer ipp,ip1,ip3,ipr1,ipr3
      common /wptest/ tamwp(100,-1:1,-1:1,4), ipp
      common /wprtest/ tamwpr(500,-1:1,-1:1,-1:1)
      
      double complex tamwm,tamwmr
      integer im1,im3,imr1,imr3
c     integer is1,is3
      common /wmtest/ tamwm(100,-1:1,-1:1)
      common /wmrtest/ tamwmr(500,-1:1,-1:1,-1:1)
