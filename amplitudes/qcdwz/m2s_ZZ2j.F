c ************************************************************************************
c Author: Matthias Kerner
c Date: 2012/13
c ************************************************************************************

!#define NINH_COMPARE
! #define CMPLX_MASS_SCHEME
      real*8 function m2s_QCDzz2j(
     &                   xi,  !in:  Feynman x parameters of incoming fermions
     &                   p,   !in:  momenta associated with external fermions
     &                   v,     !in:  momenta associated with Z/H dec. fermions
     &                   rn4, ! 1: lha subprocess, 2,3: lephel, 4: gaugeTest V or x in fin coll
     &                   Nbos, ! 22 = ZZ_ll, 21 = ZZ_lnu, 2=Z_l, -2=Z_nu
     &                   nlo
     &                        )
      implicit none
c
c declare input/output variables
c
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/cuts.inc"
#include "VBFNLO/amplitudes/vvjj/tensorz.inc"
#include "VBFNLO/loops/Utilities/GaugeTest.inc"
#include "cacheZZ.inc"
#include "common.inc"
      integer nbos,nlo
      real*8 xi(nx), p(0:3,max_p,max_kin), v(0:3,max_v,max_kin), m2s,rn4(4)
c
c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"
      real*8  alphas, xmt, alfa, xmz, xmw,s2w, xmh
      common /bkopin/   alphas, xmt, alfa, xmz, xmw, s2w, xmh
c
c declare local variables
c
      real*8 q_sf

      integer init/0/, I, mu,q,h
      save init
c
c declare variables for summation process
c
      integer FSIGN(4)
      integer id(4) ! id(feyn index) = ps index
c
c store contributions from subprocess ip in res(ip) where
      real*8 res(maxNumSubProc)
      real*8 weight
c
c declare parton distribution variables
c
      real*8 pdf(-6:6,2,2)


      double precision XGW, XGZ, XGH, GAMT,
     & BWNE,BWUD,BWTB,BZNN,BZEE,BZUU,BZDD,BZTT,
     & BHWW,BHZZ,BHGG,BHTT,BHBB,BHCC,BHTAU,BHMU,BHGAM,BHGAMZ
      COMMON /BRANCH/ BWNE,BWUD,BWTB,BZNN,BZEE,BZUU,BZDD,BZTT,
     &                BHWW,BHZZ,BHGG,BHTT,BHBB,BHCC,BHTAU,BHMU,
     &                BHGAM, BHGAMZ, XGW, XGZ, XGH, GAMT


      DOUBLE PRECISION CLR,XM2,XMG,B,VV,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),VV(4,5),A(4,5)


c local variables
      integer FL4qcdVV, ires_min,ires_max
      save ires_min, ires_max
      real*8 mass2,mjj2
      external FL4qcdVV,mass2,mjj2,dotcc
      logical ldebug,allsub
      parameter (ldebug = .false., allsub=.false.)
      complex*16 phase1,phase2,dotcc
      real*8 RNcom
      common /QCDewPhase / RNcom
      

! #include "temp.inc" 
      integer hel/1/, nhel
      save hel

      real*8 xiFinColl(2,2)


      complex*16 ec V currents


c V1/V2/V12,
c    for qhel=0:
c       AZ=1,2 : A and Z current wo coupling
c    qhel = -1:1:
c       AZ=3,4 : V coupling to u/d quark with qhel inclouding couplings
      complex*16 epsV,w4l(0:4) ! (mu, qhel, 12=AZ 34=including coupl to ud, V1 V2 V12, hel->0=rnd hel, kin)
      real*8 momV
      common /QCDZZcurrents/ epsV(0:3,-1:1,4,3,0:4,7),w4l,momV(0:4,3,7)
      logical lAZ(0:2,3) ! V.or.A A V,    V1 V2 V12
      common /QCDZZcacheL/ lAZ
      
!       if(init.ge.1.and.nlo.le.3) then
!         m2s_QCDzz2j = 0d0
!         return
!         endif


      if ( init .lt. 1 ) then

         write(6,*) " "
        if(nbos .eq. 22) then
            write(6,*) "QCD ZZ jj -> 4l amplitude square information:"
         elseif (nbos.eq.21) then
            write(6,*) "QCD ZZ jj -> 2l 2v amplitude square information:"
         elseif (nbos.eq.20) then
            write(6,*) "QCD ZA jj -> 2l A amplitude square information:"
         elseif (nbos.eq.23) then
            write(6,*) "QCD ZA jj -> 2v A amplitude square information:"
        elseif(nbos.eq.2) then
            write(6,*) "QCD Z jj -> 2l amplitude square information:"
        elseif(nbos.eq.-2) then
            write(6,*) "QCD Z jj -> 2v amplitude square information:"
         elseif (nbos.eq.24) then
            write(6,*) "QCD AA jj amplitude square information:"
        endif
         if(nfl.eq.5) then
          write(6,*) "external b quarks are included"
         else
          write(6,*) "external b quarks are excluded"
         endif
         write(6,*) "-----------------------------------------------"
         write(6,*) " "
      endif


!       if(nlo.ne.0) then
          GaugeAccuracyBoxL=GaugeAccuracy
          GaugeAccuracyBoxLnA=GaugeAccuracy
          GaugeAccuracyPenL=GaugeAccuracy
          GaugeAccuracyPenLnA=GaugeAccuracy
          GaugeAccuracyHexL=GaugeAccuracy
          GaugeAccuracyHexLnA=GaugeAccuracy
!       endif

      if ( ldebug ) then
 10      format( " p(", i1, ") = ", 4(f10.3, 2x) )
 20      format( " v(", i1, ") = ", 4(f10.3, 2x) )
c
         write(6,*) " "
         do i = 1, n_v
            write(6,20) i, v(0,i,1), v(1,i,1), v(2,i,1), v(3,i,1)
         end do
         print*," parton momenta "
         do i = 1, n_p
            write(6,10) i, p(0,i,1), p(1,i,1), p(2,i,1), p(3,i,1)
         end do
         print*,"Z momenta"
         do i=1,2
          write(6,10) i, (v(mu,2*i-1,1)+v(mu,2*i,1),mu=0,3)
         enddo
c
         write(6,*) " "
         write(6,*) "xi(1), xi(2) =", xi
         write(6,*) "mu_f =",sqrt(mufsq(1,1)),sqrt(mufsq(2,1))
      end if


c call PDF subroutine in order to determine parton
c distributions in the incoming (anti)protons.
c

      q_sf = sqrt(mufsq(1,1))
      call pdfproton( xi(1), q_sf, pdf(-6,1,1) )        ! f_a(x1)
      pdf(:,1,1) = pdf(:,1,1)/xi(1)

      q_sf = sqrt(mufsq(2,1))
      call pdfproton( xi(2), q_sf, pdf(-6,2,1) )        ! f_a(x2)
      pdf(:,2,1) = pdf(:,2,1)/xi(2)
      
      if(nlo.eq.3) then
        xiFinColl(:,1) = xi
        xiFinColl(:,2) = xi/((1d0-xi)*rn4(2)+xi) ! transform rn el [0;1] to x el [z;1], get xa = xxa/x, multiply jacobian to pdf

        q_sf = sqrt(mufsq(1,1))
        call pdfproton( xiFinColl(1,2), q_sf, pdf(-6,1,2) )        ! f_a(x1)
        pdf(:,1,2) = pdf(:,1,2)/xiFinColl(1,2)

        q_sf = sqrt(mufsq(2,1))
        call pdfproton( xiFinColl(2,2), q_sf, pdf(-6,2,2) )        ! f_a(x2)
        pdf(:,2,2) = pdf(:,2,2)/xiFinColl(2,2)

        call CScolAmpQZZReset
        call CScolAmpGNewReset
      endif


#ifdef NINH_COMPARE
        als(1,1) = 0.129808134118326d0
        als(2,1) = 0.129808134118326d0
        alphas=als(1,1)
        G = sqrt(4d0*pi*alphas)
        gg(1) = -G
        gg(2) = -G
        WMASS = 80.398d0
        ZMASS = 91.1876d0
        ZWIDTH= 2.508419926825088d0
        WWIDTH= 2.097672973437829d0
        HWIDTH=0d0
        mursq = ZMASS**2 
        mufsq = ZMASS**2
        
        gwf(1) = -sqrt( 1.16637d-5/sqrt(2d0))*80.398d0*2d0
        pdf = 1d0
! 
        do i=-6,6
          pdf(i,1,1) =1d0+0.01d0*i
          pdf(i,2,1) =1d0+0.02d0*i
        enddo

! ! ! ! ! !         p(0:3,1,1) = (/ 500.000000000000d0, 0.000000000000000d0,0.000000000000000d0,500.000000000000d0/)
! ! ! ! ! !         p(0:3,2,1) = (/ 500.000000000000d0, 0.000000000000000d0,0.000000000000000d0,-500.000000000000d0/)
! ! ! ! ! !         momW(0:3,1,1) = (/ 216.232776975442d0, -110.114342023960d0,14.5729785186567d0,-167.198192012463d0/)
! ! ! ! ! !         momW(0:3,2,1) = (/ 337.997244776963d0,-59.5891662343044d0,-150.880208615759d0,285.416525429737d0/)
! ! ! ! ! !         p(0:3,3,1) = (/158.506206722148d0,-60.0863508547814d0,136.307230097103d0,54.1681367009778d0/)
! ! ! ! ! !         p(0:3,4,1) = (/287.263771525447d0,229.789859113046d0,0.000000000000000d0,-172.386470118252d0/)
! ! ! ! ! !         p(0:3,5,1) = momW(0:3,1,1)
! ! ! ! ! !         p(0:3,6,1) = momW(0:3,2,1)
! ! ! ! ! !         momW(4,1,1) = mass2(momW(0,1,1))
! ! ! ! ! !         momW(4,2,1) = mass2(momW(0,2,1))

! ! ! ! ! !         p(0:3,1,1) = (/ 1000.000000000000d0, 0.000000000000000d0,0.000000000000000d0,1000.000000000000d0/)
! ! ! ! ! !         p(0:3,2,1) = (/ 1000.000000000000d0, 0.000000000000000d0,0.000000000000000d0,-1000.000000000000d0/)
! ! ! ! ! !         momW(0:3,1,1) = (/217.091616120748d0, -150.418551440425d0, 109.940119952077d0,77.1515443358419d0  /)
! ! ! ! ! !         momW(0:3,2,1) = (/ 973.657135114260d0,858.210753856590d0, 0.000000000000000d0,-452.789885399190d0 /)
! ! ! ! ! !         p(0:3,3,1) = (/744.278366210229d0,-648.780827758112d0, -107.384826245447d0,348.571689949576d0 /)
! ! ! ! ! !         p(0:3,4,1) = (/64.9728825547623d0, -59.0113746580532d0, -2.55529370662999d0,27.0666511137723d0 /)
! ! ! ! ! !         p(0:3,5,1) = momW(0:3,1,1)
! ! ! ! ! !         p(0:3,6,1) = momW(0:3,2,1)
! ! ! ! !         
! ! ! ! ! ! c Paper
! ! ! ! ! !     p(0:3,1,1) = (/500.000000000000d0,500.000000000000d0,0.00000000000000d0,0.00000000000000d0/)
! ! ! ! ! !     p(0:3,2,1) =  (/500.000000000000d0,-500.000000000000d0,0.00000000000000d0,0.00000000000000d0/)
! ! ! ! ! !     p(0:3,3,1) =  (/54.2314070117999d0,-31.1330162081798d0,-7.92796656791140d0,43.6912823611163d0/)
! ! ! ! ! !     p(0:3,4,1) = (/214.488870161418d0,-27.0607980217775d0,-98.5198083786150d0,188.592247959949d0/)
! ! ! ! ! !     v(0:3,1,1) = (/85.5312248384887d0,-8.22193223977868d0,36.1637837682033d0,-77.0725048002414d0/)
! ! ! ! ! !     v(0:3,2,1) = (/181.428811610043d0,-57.8599829481937d0,-171.863734086635d0,-5.61185898481311d0/)
! ! ! ! ! !     v(0:3,3,1) = (/82.8493010774356d0,-65.9095476235891d0,-49.8952157196287d0,5.51413360058664d0/)
! ! ! ! ! !     v(0:3,4,1) = (/381.470385300815d0,190.185277041519d0,292.042940984587d0,-155.113300136598d0/)
! ! ! ! ! ! 
! ! ! ! ! c 2->6
!         p(0:3,1,1) = (/ 500.000000000000d0, 0.000000000000000d0,0.000000000000000d0,500.000000000000d0/)
!         p(0:3,2,1) = (/ 500.000000000000d0, 0.000000000000000d0,0.000000000000000d0,-500.000000000000d0/)
!         p(0:3,3,1) = (/ 88.1907917270800d0, -71.8338387342512d0,-26.8500535430029d0,-43.5498562901036d0/) 
!         p(0:3,4,1) = (/ 457.182458695087d0,  404.255345805865d0, 0.000000000000000d0,-213.526148108059d0/)
!         v(0:3,1,1) = (/ 100.788423998349d0, -45.8888076188476d0, 49.1583762549171d0, 75.0731496033797d0/)
!         v(0:3,2,1) = (/    199.936105394529d0,-187.698928314700d0,-25.1678128082751d0,64.1103716125017d0/)
!         v(0:3,3,1) = (/149.590646659259d0,-96.4973730388405d0,2.67112519191297d0,114.273722503239d0/)
!         v(0:3,4,1) = (/4.31157352569544d0,-2.33639809922604d0,0.188364904448005d0,3.61876067904234d0/)

c point 2
!       v(:,1,1) = (/   128.28517384798186640d0,     12.81205307979342223d0,    -69.19348313631607539d0,    107.26229074820207643d0 /)
!       v(:,2,1) = (/   556.02626579167394993d0,   -111.18103223107226540d0,   -276.99462728273675793d0,    469.12467722212636545d0 /)
!       v(:,3,1) = (/   181.93738317191755982d0,      2.47586265168042274d0,    181.73731791511434608d0,      8.16264519600288097d0 /)
!       v(:,4,1) = (/    50.32502148264003949d0,     10.99910781031214846d0,     26.51191270377822917d0,     41.33697980492176782d0 /)
! c     parton momenta 
!       p(:,1,1) = (/   817.04583297571480216d0,      0.00000000000000000d0,      0.00000000000000000d0,    817.04583297571480216d0 /)
!       p(:,2,1) = (/   507.68698916035094726d0,      0.00000000000000000d0,      0.00000000000000000d0,   -507.68698916035094726d0 /)
!       p(:,3,1) = (/   154.60391740584950071d0,    -51.30711878197600839d0,      9.48511674338542221d0,   -145.53344426588444094d0 /)
!       p(:,4,1) = (/   253.55506043600280464d0,    136.20112747126228214d0,    128.45376305677484652d0,   -170.99430489000485522d0 /)

! ! ! ! ! ! 
c paper
      als=1d0
      alphas=1d0
      G=sqrt(4d0*pi)
      gg(1)=-G
      gg(2)=-G
      WMASS=80.385d0
      ZMASS=91.1876d0
      print*,"GWF",GWF
      GWF(1)=-sqrt(4d0*pi)/sqrt(1d0-(WMASS/ZMASS)**2)/sqrt(2d0)
      print*,"GWF",GWF
      GWF(2)=0d0
!       WWIDTH=2.097672973437829d0
!       ZWIDTH=2.508419926825088d0
      WWIDTH=2.09761d0
      ZWIDTH=2.508905d0 
      HWIDTH=4.195d-3
!       GWWA=sqrt(4d0*pi)
      print*,"GWWA",GWWA
      print*,sqrt(4d0*pi)
      GAL=GWWA
      CLR(3,1,-1:1:2) = +GWWA*2d0/3d0
      CLR(4,1,-1:1:2) = -GWWA*1d0/3d0
!       print*,"GWF",GWF
      print*,"GWWA",GWWA
      print*,GAL ,'gal'  
      print*,GAL**2/(4d0*pi) ,'alpha'  
      print*,"couplings",clr(3:4,1,-1:1:2)
      GAU(1:2) = -clr(3,1,1)
      GAD(1:2) = -clr(4,1,1)
!       call KOPPLN(1,G,G*s,s,c,z,w,q,g)
      p(0:3,1,1) = (/18.3459102072588d0, 0.0d0, 0.0d0, 18.3459102072588d0/)
      p(0:3,2,1) = (/4853.43796816526d0, 0.0d0, 0.0d0, -4853.43796816526d0/)
      p(0:3,3,1) = (/235.795970274883d0, -57.9468743482139d0, -7.096445419113396d-15, -228.564869022223d0/)
      p(0:3,4,1) = (/141.477229270568d0, -45.5048903376581d0, -65.9221967646567d0, -116.616359620580d0/)
      v(0:3,1,1) = (/276.004829895761d0, 31.4878768361538d0, -8.65306166938040d0, -274.066240646098d0/)
      v(0:3,2,1) = (/1909.28515244344d0, 29.6334571080402d0, 40.1409467910328d0, -1908.63311192893d0/)
      v(0:3,3,1) = (/2241.46026948104d0, 28.1723094714198d0, 30.2470561132914d0, -2241.07910976778d0/)
      v(0:3,4,1) = (/67.7604270068059d0, 14.1581212702582d0, 4.18725552971283d0, -66.1323669723852d0/)
! 

      p(0:3,1,1)=  (/  32.0772251055223d0, 0d0, 0d0, 32.0772251055223d0/)
      p(0:3,2,1) =   (/  2801.69305619768d0, 0d0, 0d0, -2801.69305619768d0/)
      p(0:3,3,1) =   (/  226.525314156010d0, -10.2177083492279d0, -1.251308382450315D-15, -226.294755550298d0/)
      p(0:3,4,1) =   (/  327.281588297290d0, -6.48554750244653d0, -10.1061447270513d0, -327.061219882068d0/)
      v(0:3,1,1) =   (/  646.824307052136d0, 36.0746355875450d0, -26.0379256562231d0, -645.292438579767d0/)
      v(0:3,2,1) =   (/  1598.85193997112d0, -2.88431497177613d0, 24.4490976584709d0, -1598.66239347157d0/)
      v(0:3,3,1) =   (/  34.2871318266438d0, -16.4870647640944d0, 11.6949727248035d0, 27.6949763915464d0/)
      v(0:3,4,1) = 0d0

#endif


      if(nbos.eq.24) then
        momV(0:3,1:2,1) = v(:,1:2,1)
        momV(0:3,3,1) = momV(0:3,1,1)+momV(0:3,2,1)
        momV(4,1:2,1) = 0d0
        momV(4,3,1) = mass2(momV(0,3,1))
      else
        momV(0:3,1,1) = v(0:3,1,1)+v(0:3,2,1)
        momV(0:3,2,1) = v(0:3,3,1)+v(0:3,4,1)
        momV(0:3,3,1) = momV(0:3,1,1)+momV(0:3,2,1)

        momV(4,1,1) = mass2(momV(0,1,1))
        momV(4,2,1) = mass2(momV(0,2,1))
        momV(4,3,1) = mass2(momV(0,3,1))
      endif

#ifdef CMPLX_MASS_SCHEME
c cmplx mass for e coupling
        clrc=clr
        csin2w=1d0-(WMASS**2)/(ZMASS**2)
        !csin2w=1d0-(WMASS**2-(0d0,1d0)*WMASS*WWIDTH)/(ZMASS**2-(0d0,1d0)*ZMASS*ZWIDTH)
        e=gAl(1)

        qf=-1; t3=-0.5d0
        gV=t3/2d0-qf*csin2w
        gA=-t3/2d0
        gzlc(1)=-e/sqrt(csin2w)/sqrt(1d0-csin2w)*(gV-gA)
        qf=-1; t3=0d0
        gV=t3/2d0-qf*csin2w
        gA=-t3/2d0
        gzlc(2)=-e/sqrt(csin2w)/sqrt(1d0-csin2w)*(gV+gA)
!         print*,"gzl",gzl
!         print*,"gzlc","gzlc
!         print*
        gzl=gzlc
        
c cmplx mass for u/d coupling
        qf=2d0/3d0; t3=0.5d0
        gV=t3/2d0-qf*csin2w
        gA=-t3/2d0
        clrc(3,2,-1)=e/sqrt(csin2w)/sqrt(1d0-csin2w)*(gV-gA)
        qf=2d0/3d0; t3=0d0
        gV=t3/2d0-qf*csin2w
        gA=-t3/2d0
        clrc(3,2, 1)=e/sqrt(csin2w)/sqrt(1d0-csin2w)*(gV+gA)
!         print*,clrc(3,2,-1:1:2)
!         print*,clr(3,2,-1:1:2)
        qf=-1d0/3d0; t3=-0.5d0
        gV=t3/2d0-qf*csin2w
        gA=-t3/2d0
        clrc(4,2,-1)=e/sqrt(csin2w)/sqrt(1d0-csin2w)*(gV-gA)
        qf=-1d0/3d0; t3=0d0
        gV=t3/2d0-qf*csin2w
        gA=-t3/2d0
        clrc(4,2, 1)=e/sqrt(csin2w)/sqrt(1d0-csin2w)*(gV+gA)
!         print*,clrc(4,2,-1:1:2)
!         print*,clr(4,2,-1:1:2)
!         pause
#endif


c calc Z currents
      select case(nbos)
      case(-2) ! Zn
        nhel = 1 
      case(2,21,23)
        nhel = 2
      case(22,20,24)
        nhel = 4
      case default 
        print*,__FILE__,__LINE__
        stop
      end select

      do hel=1,nhel
        ie = (-1)**(hel+1) ! hel leppair 1
        iu = sign(1,2-hel) ! hel leppair 2
c epsV(:,qhel, q,i,lhel,kin)
c   i= 1: leppair 1; 2: leppair 2; 3: V->4lep
c   lhel: 1..4,  V1 same for i=13/24; V2 same for i=12/34
c   qhel!=0: epsV including coupling to quark q
c   qhel=0: q=3: photon without coupling
c   qhel=0: q=4: Z without coupling

cmkza Leptonic tensor defined

      CALL IXXXXX(v(0,1,1),0d0 ,-ie,-1,lep)         !e+       
      CALL OXXXXX(v(0,2,1),0d0 ,+ie, 1,lem)         !e- 

      CALL JIOXXX(lep,lem,GZL ,ZMASS,ZWIDTH,ze)     !Zl
      CALL JIOXXX(lep,lem,GAL ,0d0,0d0,ae)      !Al
      epsV(:,0,3,1,hel,1) = ae(1:4)
      epsV(:,0,4,1,hel,1) = ze(1:4)
      lAZ(:,1) = .true. ! A and Z cotribution for V1
      
      if(nbos.eq.22) then
        CALL IXXXXX(v(0,3,1),0d0,-iu,-1,lup)         !mu+      
        CALL OXXXXX(v(0,4,1),0d0 ,+iu, 1,lum)         !mu- 

        CALL JIOXXX(lup,lum,GZL ,ZMASS,ZWIDTH,zu)     !Zl
        CALL JIOXXX(lup,lum,GAL ,0d0,0d0,au)      !Al
#ifdef WITH_VBF
        call vto4l(v,hel,azz,zzztens)
!         CALL HVVXXX(ze,zu,GZZH ,HMASS   ,HWIDTH  ,lep)
!         print*,lep(1)
!         vh = 2.d0*sqrt(xm2(3)*(1-xm2(3)/xm2(2)))/(-clr(2,1,-1)) ! 2*MW*SW/e
!         print*,gzzh*dotcc(zu(1),ze(1))/dcmplx(momv(4,3,1)-HMASS**2,HMASS*HWIDTH) 
!         CALL SXXXXX(momv(0,3,1),-1,lup)   
!         call VVHAXX(zu,ze,lup,2,2,lep(1)) 
!         print*,lep(1)/dcmplx(momv(4,3,1)-HMASS**2,HMASS*HWIDTH) 
!         print*
! 
!         print*,-B(2,2,6)*dotcc(zu(1),ze(1))
!         print*,2d0*xm2(2)/vh
!         print*,GZZH
!         print*,B(2,2,6)
!         pause
!         print*,lep(1)
!         w4l(hel) = -B(2,2,6)*dotcc(zu(1),ze(1))
        w4l(hel) = gzzh*dotcc(zu(1),ze(1))*(-1)/dcmplx(momv(4,3,1)-HMASS**2,HMASS*HWIDTH)  ! eff H 
!         print*,w4l(hel)/dcmplx(momv(4,3,1)-HMASS**2,HMASS*HWIDTH) 
!         
!         print*,w4l(hel)/dcmplx(momv(4,3,1)-HMASS**2,HMASS*HWIDTH) /lep(1)
!         pause
!         return
#else
        print*,'please compile with vbf processes enabled'
        stop
#endif
        epsV(:,0,3,2,hel,1) = au(1:4)
        epsV(:,0,4,2,hel,1) = zu(1:4)
        epsV(:,0,3,3,hel,1) = azz(0:3)
        epsV(:,0,4,3,hel,1) = zzztens(0:3)
        lAZ(:,2:3) = .true. ! A and Z cotribution for V2
      elseif(nbos.eq.21) then   
        CALL OXXXXX(v(0,3,1),0d0 ,-1, 1,lum)         !vm
        CALL IXXXXX(v(0,4,1),0d0 ,+1,-1,lup)         ! vm~   

        CALL JIOXXX(lup,lum,GZN ,ZMASS,ZWIDTH,zu)     !Zl
#ifdef WITH_VBF
        call vto4ln(v,hel,azz,zzztens)
        w4l(hel) = -B(2,2,6)*dotcc(zu,ze) ! - = Vertex*Prop
#else
        print*,'please compile with vbf processes enabled'
        stop
#endif
        epsV(:,0,3,2,hel,1) = 0d0
        epsV(:,0,4,2,hel,1) = zu(1:4)
        epsV(:,0,3,3,hel,1) = azz(0:3)
        epsV(:,0,4,3,hel,1) = zzztens(0:3)
        lAZ(:,2:3) = .true. 
c        lAZ(2,2) = .false. ! no Z for V2
      elseif(nbos.eq.20) then   ! ZA
        call VXXXXX(v(0,3,1),0d0,iu,1,au) ! A

        epsV(:,0,3,2,hel,1) = au(1:4)
        epsV(:,0,4,2,hel,1) = 0d0
        
        CALL FVIXXX(lep,au,GAL ,0d0    ,0d0    ,lup)
        CALL FVOXXX(lem,au,GAL ,0d0    ,0d0,lum)

        CALL JIOXXX(lup,lem,GZL ,ZMASS,ZWIDTH,zu)     !Zl
        epsV(:,0,4,3,hel,1) = zu(1:4)
        CALL JIOXXX(lep,lum,GZL ,ZMASS,ZWIDTH,zu)     !Zl
        epsV(:,0,4,3,hel,1) = epsV(:,0,4,3,hel,1)+zu(1:4)

        CALL JIOXXX(lup,lem,GAL ,0d0,0d0,zu)     !Zl
        epsV(:,0,3,3,hel,1) = zu(1:4)
        CALL JIOXXX(lep,lum,GAL ,0d0,0d0,zu)     !Zl
        epsV(:,0,3,3,hel,1) = epsV(:,0,3,3,hel,1)+zu(1:4)

!         epsV(:,0,3:4,1:2,hel,1) = 0d0

        lAZ(:,2:3) = .true. 
        lAZ(1,2) = .true. ! no photon for V2

      elseif(nbos.eq.23) then   
        CALL OXXXXX(v(0,1,1),0d0 ,-1, 1,lem)         !vm
        CALL IXXXXX(v(0,2,1),0d0 ,+1,-1,lep)         ! vm~   
        CALL JIOXXX(lep,lem,GZN ,ZMASS,ZWIDTH,ze)     !Zn
        epsV(:,0,3,1,hel,1) = 0d0
        epsV(:,0,4,1,hel,1) = ze(1:4)
        
        call VXXXXX(v(0,3,1),0d0,ie,1,au) ! A
        epsV(:,0,3,2,hel,1) = au(1:4)
        epsV(:,0,4,2,hel,1) = 0d0
        
        epsV(:,0,3:4,3,hel,1) = 0d0 ! V->4l
        
        lAZ = .true.
c======
c try if this is working
!         lAZ=.false.
!         lAZ(2,1) = .true. ! Z contribution for V1
!         lAZ(1,2) = .true. ! A contribution for V2
c try if this is working
c======
        
      elseif(nbos.eq.24) then
        call VXXXXX(v(0,1,1),0d0,ie,1,ae) 
        call VXXXXX(v(0,2,1),0d0,iu,1,au) 
        
        epsV(:,0,3,1,hel,1) = ae(1:4)
        epsV(:,0,3,2,hel,1) = au(1:4)
        epsV(:,0,3,3,hel,1) = 0d0
        
        
        lAZ = .true.
!         lAZ(0:1,1:2) = .true.
        

      elseif(abs(nbos).eq.2) then
        if(nbos.lt.0) then
        CALL OXXXXX(v(0,1,1),0d0 ,-1, 1,lem)         !vm
        CALL IXXXXX(v(0,2,1),0d0 ,+1,-1,lep)         ! vm~   
          CALL JIOXXX(lep,lem,GZN ,ZMASS,ZWIDTH,ze) ! Z-> nu nubar
          ae = 0d0
        endif
        epsV(:,0,3,3,hel,1) = ae(1:4)
        epsV(:,0,4,3,hel,1) = ze(1:4)
        lAZ(:,3) = .true. ! A and Z cotribution for V3, don't use V1 and V2
        lAZ(:,1:2) = .false.
      endif
c add A and Z current including couplings to u and d quarks
      do q=3,4
      do h=-1,1,2
      do i=1,3
        if(lAZ(0,i)) then
          epsV(:,h,q,i,hel,1) = 0d0
#ifdef CMPLX_MASS_SCHEME
          if(lAZ(1,i)) epsV(:,h,q,i,hel,1) = epsV(:,h,q,i,hel,1) + epsV(:,0,3,i,hel,1)*CLRc(q,1,h)
          if(lAZ(2,i)) epsV(:,h,q,i,hel,1) = epsV(:,h,q,i,hel,1) + epsV(:,0,4,i,hel,1)*CLRc(q,2,h) 
#else
          if(lAZ(1,i)) epsV(:,h,q,i,hel,1) = epsV(:,h,q,i,hel,1) + epsV(:,0,3,i,hel,1)*CLR(q,1,h)
          if(lAZ(2,i)) epsV(:,h,q,i,hel,1) = epsV(:,h,q,i,hel,1) + epsV(:,0,4,i,hel,1)*CLR(q,2,h) 
#endif
        endif
      enddo;enddo;enddo
! c copy photon contribution to index 1 and 2 for convenience
      epsV(:,0,1,:,:,1) = epsV(:,0,3,:,:,1) 
      epsV(:,0,2,:,:,1) = epsV(:,0,4,:,:,1) 
      
      enddo
      
      phase1=exp((0d0,1d0)*pi*(rn4(3)))
      phase2=exp((0d0,1d0)*pi*(rn4(4)))
      if(nhel.eq.4) then
        epsV(:,:,:,1,0,1) = phase1*epsV(:,:,:,1,1,1)+conjg(phase1)*epsV(:,:,:,1,2,1)
        epsV(:,:,:,2,0,1) = phase2*epsV(:,:,:,2,1,1)+conjg(phase2)*epsV(:,:,:,2,3,1)
        epsV(:,:,:,3,0,1) = phase1*(phase2*(epsV(:,:,:,3,1,1))+conjg(phase2)*epsV(:,:,:,3,3,1))
     &                     +conjg(phase1)*(phase2*(epsV(:,:,:,3,2,1))+conjg(phase2)*epsV(:,:,:,3,4,1))
        w4l(0) = phase1*(phase2*(w4l(1))+conjg(phase2)*w4l(3))
     &                     +conjg(phase1)*(phase2*(w4l(2))+conjg(phase2)*w4l(4))
      elseif(nhel.eq.2) then
         if(nbos.eq.21.or.nbos.eq.2) then
        epsV(:,:,:,1,0,1) = phase1*epsV(:,:,:,1,1,1)+conjg(phase1)*epsV(:,:,:,1,2,1)
        epsV(:,:,:,2,0,1) = epsV(:,:,:,2,2,1)
        epsV(:,:,:,3,0,1) = phase1*epsV(:,:,:,3,1,1)+conjg(phase1)*epsV(:,:,:,3,2,1)
        elseif(nbos.eq.23) then
        epsV(:,:,:,2,0,1) = phase1*epsV(:,:,:,2,1,1)+conjg(phase1)*epsV(:,:,:,2,2,1)
        epsV(:,:,:,1,0,1) = epsV(:,:,:,1,2,1)
        epsV(:,:,:,3,0,1) = 0d0
        else
          print*,__FILE__,__LINE__
          stop
        endif
      else ! Zn
        epsV(:,:,:,3,0,1) = epsV(:,:,:,3,1,1)
      endif


c reset subprocess counter
      ires_min = FL4QCDVV(id,id,-sub_number) +1


! c different family
!       id(1) = 1
!       id(2) = 2
!       id(3) = 3
!       id(4) = 4
!       
!       fsign(1:4) = 1
!       print*,'my diff family'
!       pMG(0,1)=0d0
!       call amp4qzzBorn(P,id,fsign,1,amp)
!       pmG(0,1)=pMG(0,1)+amp(0,3,3,1)
!       call amp4qzzBorn(P,id,fsign,2,amp)
!       pmG(0,1)=pMG(0,1)+amp(0,3,3,1)
!       call amp4qzzBorn(P,id,fsign,3,amp)
!       pmG(0,1)=pMG(0,1)+amp(0,3,3,1)
!       call amp4qzzBorn(P,id,fsign,4,amp)
!       pmG(0,1)=pMG(0,1)+amp(0,3,3,1)
!       print*,pMG(0,1)
! 
!       id(1)=2
!       id(2)=1
!       pMG(0,1)=0d0
!       call amp4qzzBorn(P,id,fsign,1,amp)
!       pmG(0,1)=pMG(0,1)+amp(0,3,3,2)
!       call amp4qzzBorn(P,id,fsign,2,amp)
!       pmG(0,1)=pMG(0,1)+amp(0,3,3,2)
!       call amp4qzzBorn(P,id,fsign,3,amp)
!       pmG(0,1)=pMG(0,1)+amp(0,3,3,2)
!       call amp4qzzBorn(P,id,fsign,4,amp)
!       pmG(0,1)=pMG(0,1)+amp(0,3,3,2)
!       print*,pMG(0,1)
! 
! 
! c same family
!       id(1) = 1
!       id(2) = 2
!       id(3) = 3
!       id(4) = 4
!       lokt = .true.
!       fsign(1:4) = 1
!       print*,'my same family'
!       pMG(0,1)=0d0
!       call precalcQCDZZjj(p,lokt,1,nbos)
!       call amp4qzzBorn(P,id,fsign,1,amp)
!       pmG(0,1)=pMG(0,1)+amp(0,3,3,0)
!       call precalcQCDZZjj(p,lokt,2,nbos)
!       call amp4qzzBorn(P,id,fsign,2,amp)
!       pmG(0,1)=pMG(0,1)+amp(0,3,3,0)
!       call precalcQCDZZjj(p,lokt,3,nbos)
!       call amp4qzzBorn(P,id,fsign,3,amp)
!       pmG(0,1)=pMG(0,1)+amp(0,3,3,0)
!       call precalcQCDZZjj(p,lokt,4,nbos)
!       call amp4qzzBorn(P,id,fsign,4,amp)
!       pmG(0,1)=pMG(0,1)+amp(0,3,3,0)
!       print*,pMG(0,1)
! !       
!       pMG(:,:4) = p(:,:4,1)
!       pMG(:,5:8) = v(:,:4,1)
!       print*,'mg'
!       call smatrixuuzz(pMG,amp(0,3,3,0))
!       print*,amp(0,3,3,0)
!       pause
! ! 
!       id(1) = 1
!       id(2) = 2
!       id(3) = 3
!       id(4) = 4
!       fsign=1
! ! !       
!       call amp2q2gZZborn(p,id,fsign,1,amp(0,3,3,0))
!       call amp2q2gZZborn(p,id,fsign,2,amp(0,3,3,1))
!       call amp2q2gZZborn(p,id,fsign,3,amp(0,3,3,2))
!       call amp2q2gZZborn(p,id,fsign,4,amp(0,3,4,0))
!       amp(0,3,3,0) = amp(0,3,3,0)+amp(0,3,3,1)+amp(0,3,3,2)+amp(0,3,4,0)
!       print*,amp(0,3,3,0)
! ! 
!       pMG(:,:4) = p(:,:4,1)
!       pMG(:,5:8) = v(:,:4,1)
!       print*,'mg'
!       call smatrixugugzz(pMG,amp(0,3,3,0))
!       print*,amp(0,3,3,0)
!       pause


! ! ! ZA
! !       id(1) = 1
! !       id(2) = 2
! !       id(3) = 3
! !       id(4) = 4
! !       fsign=1
! ! ! !       
! !       call amp2q2gZZborn(p,id,fsign,1,amp(0,3,3,0))
! !       call amp2q2gZZborn(p,id,fsign,2,amp(0,3,3,1))
! !       call amp2q2gZZborn(p,id,fsign,3,amp(0,3,3,2))
! !       call amp2q2gZZborn(p,id,fsign,4,amp(0,3,4,0))
! !       amp(0,3,3,0) = amp(0,3,3,0)+amp(0,3,3,1)+amp(0,3,3,2)+amp(0,3,4,0)
! !       print*,amp(0,3,3,0)
! ! ! 
! !       pMG(:,:4) = p(:,:4,1)
! !       pMG(:,5:8) = v(:,:4,1)
! !       print*,'mg'
! !       call smatrixugugza(pMG,amp(0,3,3,0))
! !       print*,amp(0,3,3,0)
! !       pause

! ! ! paper
! ! !       id(1) = 3
! ! !       id(2) = 1
! ! !       id(3) = 4
! ! !       id(4) = 2
! ! !       fsign = -1
! ! !       fsign(3) = 1
! ! !       
! ! !       ampg=0d0
! ! !       do i=-1,1,2; do j=1,2; do k=1,2; do l=1,4
! ! !         call amp2q2gZZvirt(p,id,fsign,2,nbos,i,j,k,l,0.3d0, ampgtemp)
! ! !         ampg=ampg+ampgtemp
! ! !       enddo; enddo; enddo;enddo
! ! !       print*,ampg(0,:)/2/2/8/8
! ! !       pause
! ! 
! ! !       id(1) = 1
! ! !       id(2) = 2
! ! !       id(3) = 3
! ! !       id(4) = 4
! ! !       fsign = 1
! ! !       
! ! !       amp=0d0
! ! !       do l=1,4
! ! !         call amp4qzzVirt(P,id,fsign,2,l,amptemp(0,3,3,0))
! ! !         amp=amp+amptemp
! ! !       enddo
! ! !       print*,amp(0,3,3,0)/4/9/2
! ! !       print*,amp(0,3,3,1)/4/9
! ! !       print*,amp(0,3,4,1)/4/9
! ! !       print*,amp(0,4,4,0)/4/9/2
! ! !       print*,amp(0,4,4,1)/4/9
! ! !       pause
 
c   q1 -------------- q3
c            9
c            6
c            9
c   g2 96969696969696 g4

c*****  2.  ********
c**Diagram**********  q1 g2 --->  q3  g4 W Z   **********************
c**Diagram**********  g2 q1 --->  q3  g4 W Z   **********************
c**Phasespace*******  1  2         3   4

      id(1) = 1
      id(2) = 2
      id(3) = 3
      id(4) = 4

      fsign(1) =  1
      fsign(2) = -1
      fsign(3) =  1
      fsign(4) =  1
      
      if(sub_number.eq.1.or.allsub) call getRes2Q2GZZ(p,xiFinColl,id,fsign,pdf,nbos,nlo,rn4(1),res)

! switch initial state
      id(1)=2
      id(2)=1
      
      if(sub_number.eq.1.or.allsub) call getRes2Q2GZZ(p,xiFinColl,id,fsign,pdf,nbos,nlo,rn4(1),res)

c*****  4.  ********
c**Diagram**********  qb3 g2 ---> qb1 g4 W Z   **********************
c**Diagram**********  g2 qb3 ---> qb1 g4 W Z   **********************
c**Phasespace*******  1   2        3   4

      id(1) = 3
      id(2) = 2
      id(3) = 1
      id(4) = 4

      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = -1
      fsign(4) =  1

      if(sub_number.eq.2.or.allsub) call getRes2Q2GZZ(p,xiFinColl,id,fsign,pdf,nbos,nlo,rn4(1),res)

! switch initial state
      id(2)=1
      id(3)=2

      if(sub_number.eq.2.or.allsub) call getRes2Q2GZZ(p,xiFinColl,id,fsign,pdf,nbos,nlo,rn4(1),res)

c*****  3.  ********
c**Diagram**********  q1 qb3 ---> g2  g4 W Z   **********************
c**Diagram**********  qb3 q1 ---> g2  g4 W Z   **********************
c**Phasespace*******  1   2        3   4

      id(1) = 1
      id(2) = 3
      id(3) = 2
      id(4) = 4

      fsign(1) =  1
      fsign(2) =  1
      fsign(3) =  -1
      fsign(4) =  1

      if(sub_number.eq.3.or.allsub) call getRes2Q2GZZ(p,xiFinColl,id,fsign,pdf,nbos,nlo,rn4(1),res)

! switch initial state
      id(1)=2
      id(3)=1

      if(sub_number.eq.3.or.allsub) call getRes2Q2GZZ(p,xiFinColl,id,fsign,pdf,nbos,nlo,rn4(1),res)

c*****  1.  ********
c**Diagram**********  g2 g4 ---> qbar1 Q3 W Z   **********************
c**Phasespace*******  1  2         3    4

      id(1) = 3
      id(2) = 1
      id(3) = 4
      id(4) = 2

      fsign(1) = -1
      fsign(2) = -1
      fsign(3) =  1
      fsign(4) = -1

      if(sub_number.eq.4.or.allsub) call getRes2Q2GZZ(p,xiFinColl,id,fsign,pdf,nbos,nlo,rn4(1),res)
      
c end of 2q2g processes
c begin 4q processes
c
c   q1 -------------- q2
c            9
c            6
c            9
c   q3 -------------- q4
c
c

c*****  1.  ********
c**Diagram**********  q1 q2 ---> q3 q4 W Z   **********************
c**Phasespace*******  1   2      3   4

      id(1) = 1
      id(2) = 2
      id(3) = 3
      id(4) = 4
      
      fsign = 1


      if(allsub.or.sub_number.eq.5) call getRes4QZZ(res,p,xifincoll,id,fsign,pdf,nbos,nlo )

c*****  2.  ********
c**Diagram**********  q1 qb3 ---> qb2 q4 W Z   **********************
c**Diagram**********  qb3 q1 ---> qb2 q4 W Z   **********************
c**Phasespace*******  1   2        3   4

      id(1) = 1
      id(2) = 3
      id(3) = 4
      id(4) = 2
      
      fsign(1) = 1
      fsign(2) = -1
      fsign(3) = 1
      fsign(4) = -1

      if(allsub.or.sub_number.eq.6) call getRes4QZZ(res,p,xifincoll,id,fsign,pdf,nbos,nlo )
      
      id(4) = 1
      id(1) = 2
      if(allsub.or.sub_number.eq.6) call getRes4QZZ(res,p,xifincoll,id,fsign,pdf,nbos,nlo )

c*****  3.  ********
c**Diagram**********  qb3 qb4 ---> qb1 qb2 W Z   **********************
c**Phasespace*******  1   2        3   4

      id(1) = 3
      id(2) = 4
      id(3) = 1
      id(4) = 2
      
      fsign(1:4) = -1

      if(allsub.or.sub_number.eq.7) call getRes4QZZ(res,p,xifincoll,id,fsign,pdf,nbos,nlo )
c
c  end of subprocess evaluation
c
c

      ires_max = numcolor
      if(init.lt.sub_number) then
         if(ldebug) print*, "Number of Subprocesses: ", numcolor
!          init =init+1
          init = max(1,sub_number)
      endif

      m2s = 0d0
      do i=ires_min,ires_max
         m2s=m2s + res(i)
      enddo

c select subprocess
      iprocess =0
      weight = 0d0
      i=ires_min
      do while((m2s.gt.0d0).and.(i.le.ires_max).and.(iprocess.eq.0))
         weight = weight+res(i)/m2s
         if(weight.ge.rn4(1)) then
            iprocess = i
         endif
         i = i+1
      enddo

      m2s = m2s*(als(1,1)/alphas)**2

      m2s_QCDzz2j =m2s
!       print*,m2s
!       pause

      RETURN
      END


      subroutine getRes2Q2GZZ(p,xiFinColl,id,fsign,pdf,Nbos,nlo,rn2,res)
      implicit none
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/global.inc"
      real*8 p(0:3,max_p,max_kin), pdf(-6:6,2,2), xiFinColl(2,2),rn2(2)! for gauge test
      integer id(4), fsign(4),nbos,nlo,hel
      real*8 res(maxNumSubProc)  !out: subprocess cross-sections

      integer FL4qcdVV
      external FL4qcdVV
      
      real*8 amp(0:2,3:4), polcol
      integer colstr,q,i,k,j,qt
      LOGICAL ldebug
      parameter (ldebug=.false.)
      integer iflav(1:4)
      real*8 QCDVVjjFinCollBornKin

      if(nlo.le.1) then

#ifdef NINH_COMPARE
        amp=0d0
        do i=-1,1,2; do k=1,2; do j=1,2; do q=1,4
        call amp2q2gZZvirt(p,id,fsign,2,nbos,i,j,k,q,rn2(2), amptemp(0,3,1))
        amp=amp+amptemp(:,:,1)
        enddo;enddo;enddo;enddo
        print*,"amplitude",amp/8/8/2/2
        pause
#else
!       res=0d0
!       return
        call amp2q2gZZborn(p,id,fsign,0,amp)

!         hel=8d0*rn2(1)
!         if(hel.gt.7) then
!           print*,'error in ',__FILE__,__LINE__
!           stop
!         endif
!         k=mod(hel,2)+1; hel=hel/2
!         j=mod(hel,2)+1; hel=hel/2
!         i=1-2*hel
!         call amp2q2gZZvirt(p,id,fsign,4,nbos,i,j,k,0,rn2(2), amp)
!         amp=amp*8      
! ! !         hel=8*rn2(1)
! ! !         if(hel.gt.7) then
! ! !           print*,'error in ',__FILE__,__LINE__
! ! !           stop
! ! !         endif
! ! !         k=mod(hel,2)+1; hel=hel/2
! ! !         j=mod(hel,2)+1; hel=hel/2
! ! !         i=1-2*hel
! ! !         call amp2q2gZZvirt(p,id,fsign,4,nbos,i,j,k,0,rn2(2), amp)
! ! !         amp=amp*8      
#endif
      
      elseif(nlo.ne.3) then ! use fincoll to calc nlo.eq.3
        hel=8d0*rn2(1)
        if(hel.gt.7) then
          print*,'error in ',__FILE__,__LINE__
          stop
        endif
        k=mod(hel,2)+1; hel=hel/2
        j=mod(hel,2)+1; hel=hel/2
        i=1-2*hel
        call amp2q2gZZvirt(p,id,fsign,nlo,nbos,i,j,k,0,rn2(2), amp)
        amp=amp*8      
      endif


c get pol. sum/avg.
      polcol = 0.25d0
      if (id(2).le.2) then
        polcol = polcol/8d0
      else
        polcol = polcol/3d0
      endif
      if (id(4).le.2) then
        polcol = polcol/8d0
      else
        polcol = polcol/3d0
      endif
      if(id(4)+id(2).eq.7) polcol=polcol/2d0 ! final statge g g
      

c sum over q flavors and colors
      do q = 1,nfl ! quark
        iflav(id(1)) = q*fsign(1)
        iflav(id(3)) = q*fsign(3)
        qt=3+mod(q,2)
        do colstr = 1,2
          iflav(id(2)) = 21
          iflav(id(4)) = 21
          k = FL4qcdVV(iflav,id,colstr)
          iflav(id(2)) = 0
          iflav(id(4)) = 0
          if(nlo.eq.3) then
            res(k) = 0d0
            if(colstr.eq.1) res(k) = QCDVVjjFinCollBornKin(p,xiFinColl,id,fsign,iflav,0,nbos,pdf) ! born is included here
          else
            res(k) = amp(colstr,qt)*POLCOL*
     &               pdf(sign1*iflav(1),1,1)*pdf(sign2*iflav(2),2,1)
          endif
          if(ldebug .and.colstr.eq.1)  print*,k, iflav,amp(0,qt)*Polcol,amp(0,qt)*POLCOL*
     &               pdf(sign1*iflav(1),1,1)*pdf(sign2*iflav(2),2,1)
        enddo
      enddo
      end

      subroutine getRes4QZZ(res,p,xifincoll,id,fsign,pdf,Nbos,nlo)
      implicit none
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/global.inc"
#include "common.inc"
      real*8 res(1:maxNumSubProc) !out: subprocess cross-sections
      real*8 p(0:3,max_p,max_kin)
      integer fsign(4), id(4),nlo
      real*8 pdf(-6:6,2,2)
      integer nbos

      integer FL4qcdVV
      external FL4qcdVV

      real*8 POLCOL, polcol0    ! ave over intial polarizations and colors
      real*8 xiFinColl(2,2)
      integer j,i,k
      integer iflav(1:4),q1,q2
      integer temp
      logical ldebug,idFl
      parameter (ldebug = .false.)
      real*8 amp(0:2,3:4,3:4,0:2),QCDVVjjFinCollBornKin

      external QCDVVjjFinCollBornKin
      
      polcol0 = 1.0d0/(4.0d0 * 9.0d0)

!       j=RandomNumber()*4+1
      if(nlo.eq.0) then
!       res=0d0
!       return
        call amp4qzzBorn(P,id,fsign,0,amp)
!           call amp4qzzVirt(P,id,fsign,4,0,amp)
!         call amp4qzzBorn(P,id,fsign,1,amptemp(0,3,3,0,1))
!         call amp4qzzBorn(P,id,fsign,2,amptemp(0,3,3,0,2))
!         amp=amptemp(:,:,:,:,1)+amptemp(:,:,:,:,2)
!       print*,amp
#ifdef NINH_COMPARE
!         call amp4qzzBorn(P,id,fsign,1,amptemp(0,3,3,0,1))
!         call amp4qzzBorn(P,id,fsign,2,amptemp(0,3,3,0,2))
!         call amp4qzzBorn(P,id,fsign,3,amptemp(0,3,3,0,3))
!         call amp4qzzBorn(P,id,fsign,4,amptemp(0,3,3,0,4))
!         amp=amptemp(:,:,:,:,1)+amptemp(:,:,:,:,2)+amptemp(:,:,:,:,3)+amptemp(:,:,:,:,4)


          call amp4qzzVirt(P,id,fsign,nlo,1,amptemp(0,3,3,0,1))
          call amp4qzzVirt(P,id,fsign,nlo,2,amptemp(0,3,3,0,2))
          call amp4qzzVirt(P,id,fsign,nlo,3,amptemp(0,3,3,0,3))
          call amp4qzzVirt(P,id,fsign,nlo,4,amptemp(0,3,3,0,4))
          amp=amptemp(:,:,:,:,1)+amptemp(:,:,:,:,2)+amptemp(:,:,:,:,3)+amptemp(:,:,:,:,4)
          print*,amp(0,3,3,0)/4/9/2
          print*,amp(0,4,4,0)/4/9/2
          print*,amp(1,3,3,1)/4/9
          print*,amp(1,3,4,1)/4/9
          print*,amp(1,4,3,1)/4/9
          print*,amp(1,4,4,1)/4/9
          pause
#endif
      else
        if(nlo.ne.3) then
#ifdef NINH_COMPARE
          call amp4qzzVirt(P,id,fsign,nlo,1,amptemp(0,3,3,0,1))
          call amp4qzzVirt(P,id,fsign,nlo,2,amptemp(0,3,3,0,2))
          call amp4qzzVirt(P,id,fsign,nlo,3,amptemp(0,3,3,0,3))
          call amp4qzzVirt(P,id,fsign,nlo,4,amptemp(0,3,3,0,4))
          amp=amptemp(:,:,:,:,1)+amptemp(:,:,:,:,2)+amptemp(:,:,:,:,3)+amptemp(:,:,:,:,4)
          print*,amp(0,3,3,0)/4/9/2
          pause
#else
          call amp4qzzVirt(P,id,fsign,nlo,0,amp)
#endif

        endif
      endif
      
c sum over all quark combinations
#ifdef NINH_COMPARE
      do j = 1,4!nfl  ! quark family 1 and 2
      do k = 1,4!nfl
#else
      do j = 1,nfl  ! quark family 1 and 2
      do k = 1,nfl
#endif
        IdFl = j.eq.k
        q1=3+mod(j,2)
        q2=3+mod(k,2)
c fill flavor information for pdf and lha output
        iflav(id(1)) = j*fsign(1)
        iflav(id(3)) = j*fsign(3)
        iflav(id(2)) = k*fsign(2)
        iflav(id(4)) = k*fsign(4)

c  pick polcol
        if(iflav(3).eq.iflav(4)) then
            polcol = polcol0/2.0d0
        else
            polcol = polcol0
        endif

c get m2s*pdf
        if(idFl) then  ! same q families
          i = FL4qcdVV(iflav,id,1)
          if(nlo.eq.3) then
            res(i) = QCDVVjjFinCollBornKin(p,xiFinColl,id,fsign,iflav,0,nbos,pdf) ! born is included here
          else
            res(i) = POLCOL * amp(0,q1,q2,0)
     $                *pdf(sign1*iflav(1),1,1)*pdf(sign2*iflav(2),2,1)
          endif
          if(ldebug) print*,i, iflav,POLCOL * amp(0,q1,q2,0),POLCOL * amp(0,q1,q2,0)
     $                *pdf(sign1*iflav(1),1,1)*pdf(sign2*iflav(2),2,1)
!             print*,iflav
        else
          i = FL4qcdVV(iflav,id,1)
          if(nlo.eq.3) then
            res(i) = QCDVVjjFinCollBornKin(p,xiFinColl,id,fsign,iflav,1,nbos,pdf) ! born is included here
          else
            res(i) = POLCOL * amp(0,q1,q2,1)
     $                *pdf(sign1*iflav(1),1,1)*pdf(sign2*iflav(2),2,1)
          endif
          if(ldebug) print*,i, iflav,POLCOL * amp(0,q1,q2,1),POLCOL * amp(0,q1,q2,1)
     $                *pdf(sign1*iflav(1),1,1)*pdf(sign2*iflav(2),2,1)
!             print*,iflav
          if(iflav(1)*iflav(2).lt.0) then ! include colstr 2 for q qb initial states
            ! swap particle 1 and 2; change sign if needed
            temp         = sign(iflav(id(1)),fsign(2))
            iflav(id(1)) = sign(iflav(id(2)),fsign(1))
            iflav(id(2)) = temp
            i = FL4qcdVV(iflav,id,2)
            if(nlo.eq.3) then
              res(i) = QCDVVjjFinCollBornKin(p,xiFinColl,id,fsign,iflav,2,nbos,pdf) ! born is included here
            else
              res(i) = POLCOL * amp(0,q1,q2,2)
     $                *pdf(sign1*iflav(1),1,1)*pdf(sign2*iflav(2),2,1)
            endif
            if(ldebug) print*,i, iflav,POLCOL * amp(0,q1,q2,2),POLCOL * amp(0,q1,q2,2)
     $                *pdf(sign1*iflav(1),1,1)*pdf(sign2*iflav(2),2,1)
!             print*,iflav
          endif
        endif
      enddo
      enddo
      return
      end

       
c 4q WW Amplitude
c    q1 >>> #### >>> q3 
c           ####          , WW
c           ####
c    q2 >>> #### >>> q4  
c   res(colstr, q1,q2,flav config=0,2)
c      colstr=0..2: for catani seymour: <|>,  <|T1T2|>,  <|T1T3|>
c                   for lo: col summed, and the 2 different color flows
c      q1/2=3,4 for u and d:
c      flavour: 0: same family (t and u channel)
c      flavour: 1: diff family (only t channel), q1 corresponds to qline q1>q3
c      flavour: 2: diff family (only u channel), q1 corresponds to qline q2>q3
      subroutine amp4qzzCS(P,kin,id,fsign,lephel,res)
      implicit none
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "cacheZZ.inc"
      real*8 P(0:3,max_p,max_kin),res(0:2,3:4,3:4,0:2)
      integer kin, id(4),fsign(4),lephel
      logical lnlo,luseprecalc
      
      real*8 pbar(0:3,4),pdia(0:3,4),momTemp(0:4)
      complex*16 psitemp(2,-1:1,4)
      complex*16 currij(0:3,-1:1,1:2,3:4), currVij(0:3,-1:1,3:4,1:2,3:4,3)
      complex*16 amp(2)

      complex*16 tempCurr(0:3,-1:1,4)
      complex*16 dotcc

      real*8 qjj2,qsq(1:2,3:4,0:1)
      real*8 squared(2,2)

      
      integer mu,tempi4(4), lepkin,colstr,h1,h2,q,q1,q2
      integer ikin, in,out,iV,hel,iV2,in2,out2

#include "VBFNLO/utilities/scales.inc"
      real*8  alphas, xmt, alfa, xmz, xmw,s2w, xmh
      common /bkopin/   alphas, xmt, alfa, xmz, xmw, s2w, xmh

c V currents
c V1/V2/V12,
c    for qhel=0:
c       AZ=1,2 : A and Z current wo coupling
c    qhel = -1:1:
c       AZ=3,4 : V coupling to u/d quark with qhel inclouding couplings
      complex*16 epsV,w4l(0:4) ! (mu, qhel, 12=AZ 34=including coupl to ud, V1 V2 V12, hel->0=rnd hel, kin)
      real*8 momV
      common /QCDZZcurrents/ epsV(0:3,-1:1,4,3,0:4,7),w4l,momV(0:4,3,7)
      logical lAZ(0:2,3) ! V.or.A A V,    V1 V2 V12
      common /QCDZZcacheL/ lAZ
      

      lnlo = .true.
      ikin = kin
      lusePreCalc=kin.ne.1 ! don't use precalc for finColl
      goto 11
      
      entry amp4qzzBorn(P,id,fsign,lephel,res)
      ikin = 1
      lnlo= .false.
      lusePreCalc=.false.
      
 11   continue
      lepkin = ikin
      if(ikin .gt.7) lepkin = 1

c in   q1 >>> #### >>> q3  out
c             ####          , WZ
c             ####
c in2  q2 >>> #### >>> q4  out2

      do mu=0,3
        pbar(mu,1)=p(mu,id(1),ikin) ! change convetnion for pbar for call of psi0m
        pbar(mu,2)=p(mu,id(3),ikin)
        pbar(mu,3)=p(mu,id(2),ikin)
        pbar(mu,4)=p(mu,id(4),ikin)
        pdia(mu,1)=pbar(mu,1)*fsign(1)
        pdia(mu,2)=pbar(mu,3)*fsign(2)
        pdia(mu,3)=pbar(mu,2)*fsign(3)
        pdia(mu,4)=pbar(mu,4)*fsign(4)
      enddo

      if(.not.lusePreCalc) then
      tempi4(1) = fsign(1)
      tempi4(2) = fsign(3)
      tempi4(3) = fsign(2)
      tempi4(4) = fsign(4)
c get bra and ket 
      call psi0m(4,pbar,tempi4,psitemp)
      psi(:,:,id(1),ikin) = psitemp(:,:,1)
      psi(:,:,id(2),ikin) = psitemp(:,:,3)
      psi(:,:,-id(3),ikin) = psitemp(:,:,2)
      psi(:,:,-id(4),ikin) = psitemp(:,:,4)
      
!       psiV = 0d0

c contract q with V
      do in=1,2 ! incoming quarks
      do iV=3,1,-1 ! V1, V2, V12
      do q=3,4 ! u/d
      do hel=-1,1,2
      if(lAZ(0,iV)) then
        call ket2c(psi(1,hel,id(in),ikin),.true.,pdia(0,in),hel,momV(0,iV,lepkin),epsV(0,hel,q,iV,lephel,lepkin),
     $             psiV(1,hel,q,id(in),iV,ikin),mompsiV(0,id(in),iV,ikin))
        out = in+2
        call bra2c(psi(1,hel,-id(out),ikin),.true.,pdia(0,out),hel,momV(0,iV,lepkin),epsV(0,hel,q,iV,lephel,lepkin),
     $             psiV(1,hel,q,-id(out),iV,ikin),mompsiV(0,-id(out),iV,ikin))
c contract with 2nd V
        if(iV.le.2) then
          iV2= 3-iV
          if(lAZ(0,iV2)) then
            call ket2c(psiV(1,hel,q,id(in),iV,ikin),.false.,mompsiV(0,id(in),iV,ikin),hel,momV(0,iV2,lepkin),
     $                 epsV(0,hel,q,iV2,lephel,lepkin),psitemp(1,hel,1),momtemp)
            psiV(1:2,hel,q,id(in),3,ikin)= psiV(1:2,hel,q,id(in),3,ikin)+psitemp(:,hel,1)
            call bra2c(psiV(1,hel,q,-id(out),iV,ikin),.false.,mompsiV(0,-id(out),iV,ikin),hel,momV(0,iV2,lepkin),
     $                 epsV(0,hel,q,iV2,lephel,lepkin),psitemp(1,hel,1),momtemp)
            psiV(1:2,hel,q,-id(out),3,ikin)= psiV(1:2,hel,q,-id(out),3,ikin)+psitemp(:,hel,1)
          endif
        endif
      endif
      enddo;enddo;enddo;enddo
!       print*,psiV(1:2,-1:1:2,:,id(1),1:2,1)
!       print*,psiV(:,-1,3,id(1),2,ikin)

      endif ! usePreCalc
      
c build currents
      do in=1,2
      do out=3,4
c <q|mu|q>
        call curr(1,psi(1,-1,-id(out),ikin),psi(1,-1,id(in),ikin),currij(0,-1,in,out))
        do q=3,4
c <q|mu|Vq> + <qV|mu|q>
          do iV=1,2
          if(lAZ(0,iV)) then
            call curr(1,psi(1,-1,-id(out),ikin),psiV(1,-1,q,id(in),iV,ikin),tempCurr(0,-1,1))
            call curr(1,psiV(1,-1,q,-id(out),iV,ikin),psi(1,-1,id(in),ikin),tempCurr(0,-1,2))
            currVij(:,:,q,in,out,iV) = tempCurr(:,:,1) + tempCurr(:,:,2)
          endif
          enddo
c <q|mu|VVq> + <qVV|mu|q> + <qV|mu|Vq>
!           if(lAZ(0,3)) then ! should always be true
            call curr(1,psi(1,-1,-id(out),ikin),psiV(1,-1,q,id(in),3,ikin),tempCurr(0,-1,1))
            call curr(1,psiV(1,-1,q,-id(out),3,ikin),psi(1,-1,id(in),ikin),tempCurr(0,-1,2))
            currVij(:,:,q,in,out,3) = tempCurr(:,:,1) + tempCurr(:,:,2)
!           else
!             currVij(:,:,q,in,out,3) = 0d0
!           endif
          if(lAZ(0,1).and.lAZ(0,2)) then
            call curr(1,psiV(1,-1,q,-id(out),1,ikin),psiV(1,-1,q,id(in),2,ikin),tempCurr(0,-1,3))
            call curr(1,psiV(1,-1,q,-id(out),2,ikin),psiV(1,-1,q,id(in),1,ikin),tempCurr(0,-1,4))
            currVij(:,:,q,in,out,3) = currVij(:,:,q,in,out,3) + tempCurr(:,:,3) + tempCurr(:,:,4)
          endif
        enddo
        qsq(in,out,0) = qjj2(pdia(0,in),pdia(0,out))
        qsq(in,out,1) = qjj2(pdia(0,in),mompsiV(0,-id(out),1,ikin))
      enddo
      enddo
!       print*,currij(:,-1,3,2)/qsq(3,2,0)*G
!       print*,dotcc(currij(0,-1,3,2),currVij(0,-1,3,1,4,3))/qsq(3,2,0)*G**2
        
c calc amplitudes

      res=0d0

      do q1=3,4
      do q2=3,4
c sum over q helicities
      do h1=-1,1,2
      do h2=-1,1,2
      
      amp = 0d0
      out=3 ! fix final state
      out2=4
      do in=1,2
        in2=3-in
        colstr = in
!         print*,currVij(0,h2,q2,in2,out2,3)
        
        if(lAZ(0,1)) then
c in>V1>out, in2>V2>out2
          amp(colstr) = amp(colstr) + dotcc(currVij(0,h1,q1,in,out,1),currVij(0,h2,q2,in2,out2,2))/qsq(in,out,1)
c in>V2>out, in2>V1>out2
          amp(colstr) = amp(colstr) + dotcc(currVij(0,h1,q1,in,out,2),currVij(0,h2,q2,in2,out2,1))/qsq(in2,out2,1)
        endif
!         if(lAZ(0,3)) then ! should always be true
c in>V12>out, in2>out2
          amp(colstr) = amp(colstr) + dotcc(currVij(0,h1,q1,in,out,3),currij(0,h2,in2,out2))/qsq(in2,out2,0)
c in>out, in2>V12>out2
          amp(colstr) = amp(colstr) + dotcc(currij(0,h1,in,out),currVij(0,h2,q2,in2,out2,3))/qsq(in,out,0)
!         endif
      enddo

c square amplitudes
      amp = amp*G**2
      
      amp(2) = -amp(2)
      squared(1,1) = dble(amp(1)*dconjg(amp(1)))
      squared(2,2) = dble(amp(2)*dconjg(amp(2)))
      if(h1.eq.h2) then
        squared(1,2) = dble(amp(1)*dconjg(amp(2)))
      else
        squared(1,2) = 0d0
      endif
!       print*,amp

      res(0,q1,q2,1) = res(0,q1,q2,1)+2d0*squared(1,1)
      res(0,q1,q2,2) = res(0,q1,q2,2)+2d0*squared(2,2)
      if(q1.eq.q2) then
        res(0,q1,q2,0) = res(0,q1,q2,0)+2d0*(squared(1,1)+squared(2,2)-2d0/3d0*squared(1,2))
      endif

      if(lnlo) then ! CS color factors
        if(q1.eq.q2) then
          res(1,q1,q2,0) = res(1,q1,q2,0)-2d0/3d0*(squared(1,1)+squared(2,2)-10d0/3d0*squared(1,2))
          res(2,q1,q2,0) = res(2,q1,q2,0)+1d0/3d0*(squared(1,1)-7d0*squared(2,2)-2d0/3d0*squared(1,2))
        endif
        res(1,q1,q2,1) = res(1,q1,q2,1)-2d0/3d0*squared(1,1)
        res(2,q1,q2,1) = res(2,q1,q2,1)+1d0/3d0*squared(1,1)
        res(1,q1,q2,2) = res(1,q1,q2,2)-2d0/3d0*squared(2,2)
        res(2,q1,q2,2) = res(2,q1,q2,2)-7d0/3d0*squared(2,2)
      else ! get colorstr weighted amp for lha output
        res(1,q1,q2,1) = res(0,q1,q2,1)*squared(1,1)/(squared(1,1)+squared(2,2))
      endif
      
      enddo; enddo ! hel
      enddo; enddo ! q1, q2

      end       


c q1 q3 g2 g4 WZ Amplitude
c    ug2dgwzCS
c       Amplitude as needed for Catani Seymor with specified helicities (see arXiv:0905.0883)
c       resDiag: real Matrix elements <l|l>, <l|T1T2|l>, <l|T1T3|l>
c       resOffd: complex Matrix elements <-|+>,<-|T1T2|+>,<-|T1T3|+>, where the gluon, that is involved in dipoles has hel. index 0
c       all hel. configurations will summed if ghel_i =0 
c    ug2dgwzBorn
c       Amplitude as needed for Born process
c       res(0:2): amplitude for specified hel. configuration or summed if ghel_i = 0, index 1 and 2 for colstr for LHA output
      SUBROUTINE amp2q2gZZcs(P,kin,id,fsign,lephel, usePreCalc,resDiag) ! sum over glu hel if ghel ==0
      implicit none
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "cacheZZ.inc"
      real*8 P(0:3,max_p,max_kin), resDiag(0:2,0:2,0:2,3:4), res(0:2,3:4)
      integer kin, id(4),fsign(4),lephel, ikin
      logical lnlo,usePreCalc,lusePreCalc
      
      real*8 pbar(0:3,2),pdia(0:4,4),mpdia(0:4,4),momtemp(0:4),epsGG(0:3,2,2)

      complex*16  amp(2,2,2),psitemp(2,-1:1,2)
      complex*16 s1c,s1r
      real*8 dotrr,mjj2
      
      complex*16 psiVGlocal(2,-1:1,3:4,2,5,2)
      integer i,j,mu, pol1,pol2,pol(2), tempi2(2), lepkin,idg1,idg2
      integer glu1,glu2,iV,iV2, q, hel

c V currents
c V1/V2/V12,
c    for qhel=0:
c       AZ=1,2 : A and Z current wo coupling
c    qhel = -1:1:
c       AZ=3,4 : V coupling to u/d quark with qhel inclouding couplings
      complex*16 epsV,w4l(0:4) ! (mu, qhel, 12=AZ 34=including coupl to ud, V1 V2 V12, hel->0=rnd hel, kin)
      real*8 momV
      common /QCDZZcurrents/ epsV(0:3,-1:1,4,3,0:4,7),w4l,momV(0:4,3,7)
      logical lAZ(0:2,3) ! V.or.A A V,    V1 V2 V12
      common /QCDZZcacheL/ lAZ


      lnlo = .true.
      ikin = kin
      lusePreCalc = usePreCalc
      goto 11

      entry amp2q2gZZborn(P,id,fsign,lephel, res)
      lusePreCalc = .false.
      ikin = 1
      lnlo= .false.
      
      res=0d0
 11   continue 

      if(id(2).ge.id(4)) then
        print*,__FILE__,__LINE__
        stop
      endif


      lepkin = ikin
      if(ikin .gt.7) lepkin = 1


c use pdia for psi, -pdia for psibar
c -pdia for bosons !!!!!
      pdia(:3,1:2) = p(:,1:2,ikin) ! psi in initial state has sign 1
      pdia(:3,3:4) =-p(:,3:4,ikin) ! psi in final state has sign -1
      pdia(4,:) = 0d0
      mpdia = -pdia

      if(.not. lusePreCalc) then


c get bra and ket 
      tempi2(1) = merge(1,-1,id(1).le.2)
      tempi2(2) = merge(1,-1,id(3).ge.3)
      pbar(:,1) = p(:,id(1),ikin)
      pbar(:,2) = p(:,id(3),ikin)
      call psi0m(2,pbar,tempi2,Psi(1,-1,id(1),ikin))
      psi(:,:,-id(3),ikin) = psi(:,:,id(1)+1,ikin)

      do i=2,4,2
        call polvec(p(0,id(i),ikin),1,epsG(0,id(i),1,ikin))
        call polvec(p(0,id(i),ikin),2,epsG(0,id(i),2,ikin))
      enddo

!         call helvec(p(0,id(2),ikin),-1,-1,epsG(0,id(2),1,ikin))
!         call helvec(p(0,id(4),ikin),1, -1,epsG(0,id(4),1,ikin))
!         call helvec(p(0,id(2),ikin),-1, 1,epsG(0,id(2),2,ikin))
!         call helvec(p(0,id(4),ikin),1,  1,epsG(0,id(4),2,ikin))

c psiG
      do j=2,4,2
      do hel=-1,1,2
        call bra2r(psi(1,hel,-id(3),ikin),.true.,mpdia(0,id(3)),hel,mpdia(0,id(j)),epsG(0,id(j),1,ikin)
     &            ,psiG(1,hel,-id(3),id(j),1,ikin),momPsiG(0,-id(3),id(j),ikin))
        call bra2r(psi(1,hel,-id(3),ikin),.true.,mpdia(0,id(3)),hel,mpdia(0,id(j)),epsG(0,id(j),2,ikin)
     &            ,psiG(1,hel,-id(3),id(j),2,ikin),momPsiG(0,-id(3),id(j),ikin))

        call ket2r(psi(1,hel, id(1),ikin),.true.,pdia(0,id(1)),hel,mpdia(0,id(j)),epsG(0,id(j),1,ikin)
     &            ,psiG(1,hel, id(1),id(j),1,ikin),momPsiG(0,id(1),id(j),ikin))
        call ket2r(psi(1,hel, id(1),ikin),.true.,pdia(0,id(1)),hel,mpdia(0,id(j)),epsG(0,id(j),2,ikin)
     &            ,psiG(1,hel, id(1),id(j),2,ikin),momPsiG(0,id(1),id(j),ikin))
      enddo; enddo
      
      psiV=0d0

      do iV=3,1,-1 ! V1, V2, V12
      do q=3,4 ! u/d
      do hel=-1,1,2
      if(lAZ(0,iV)) then
        call ket2c(psi(1,hel,id(1),ikin),.true.,pdia(0,id(1)),hel,momV(0,iV,lepkin),epsV(0,hel,q,iV,lephel,lepkin),
     $             psiV(1,hel,q,id(1),iV,ikin),mompsiV(0,id(1),iV,ikin))
        call bra2c(psi(1,hel,-id(3),ikin),.true.,mpdia(0,id(3)),hel,momV(0,iV,lepkin),epsV(0,hel,q,iV,lephel,lepkin),
     $             psiV(1,hel,q,-id(3),iV,ikin),mompsiV(0,-id(3),iV,ikin))
        if(iV.eq.3) then ! this is one contribution to qVV
          psiV(:,hel,q, id(1),3,ikin) =  psiV(:,hel,q, id(1),3,ikin)
          psiV(:,hel,q,-id(3),3,ikin) =  psiV(:,hel,q,-id(3),3,ikin)
        endif
c contract with 2nd V
        if(iV.le.2) then
          iV2= 3-iV
          if(lAZ(0,iV2)) then
            call ket2c(psiV(1,hel,q,id(1),iV,ikin),.false.,mompsiV(0,id(1),iV,ikin),hel,momV(0,iV2,lepkin),
     $                 epsV(0,hel,q,iV2,lephel,lepkin),psitemp(1,hel,1),momtemp)
            psiV(1:2,hel,q,id(1),3,ikin)= psiV(1:2,hel,q,id(1),3,ikin)+psitemp(:,hel,1)
            call bra2c(psiV(1,hel,q,-id(3),iV,ikin),.false.,mompsiV(0,-id(3),iV,ikin),hel,momV(0,iV2,lepkin),
     $                 epsV(0,hel,q,iV2,lephel,lepkin),psitemp(1,hel,1),momtemp)
            psiV(1:2,hel,q,-id(3),3,ikin)= psiV(1:2,hel,q,-id(3),3,ikin)+psitemp(1:2,hel,1)
          endif
        endif
      endif
      enddo;enddo;enddo

      endif ! .not. usePrecalc

c 3glu vert
      do pol1=1,2; do pol2=1,2
      do mu = 0,3
c gluon copplings
        epsGG(mu,pol1,pol2)= (2d0*dotrr(pdia(0,id(4)),epsG(0,id(2),pol1,ikin))*epsG(mu,id(4),pol2,ikin)
     &              -2d0*dotrr(pdia(0,id(2)),epsG(0,id(4),pol2,ikin))*epsG(mu,id(2),pol1,ikin)
     &        +dotrr(epsG(0,id(2),pol1,ikin),epsG(0,id(4),pol2,ikin))*(pdia(mu,id(2))-pdia(mu,id(4))))
     &        /mjj2(pdia(0,id(2)),pdia(0,id(4)))
      enddo
      enddo; enddo


      do iV=2,1,-1 ! V1, V2
      do q=3,4 ! u/d
      do hel=-1,1,2
      if(lAZ(0,iV)) then
c contract with g
          do j=2,4,2; do pol1=1,2
            call bra2r(psiV(1,hel,q,-id(3),iV,ikin),.false.,mompsiV(0,-id(3),iV,ikin),hel,mpdia(0,id(j)),
     $                 epsG(0,id(j),pol1,ikin),psitemp(1,hel,1),momtemp)
            call bra2c(psiG(1,hel,-id(3),id(j),pol1,ikin),.false.,momPsiG(0,-id(3),id(j),ikin),hel,momV(0,iV,lepkin),
     $                 epsV(0,hel,q,iV,lephel,lepkin),psitemp(1,hel,2),momtemp)
            psiVGlocal(1:2,hel,q,iV,id(j),pol1)=psitemp(1:2,hel,1)+psitemp(1:2,hel,2)
          enddo; enddo
      endif
      enddo;enddo;enddo

      if(lnlo) then
        do q=3,4
        do pol1=0,2
        do pol2=0,2
          resDiag(0,pol1,pol2,q) = 0d0
          resDiag(1,pol1,pol2,q) = 0d0
          resDiag(2,pol1,pol2,q) = 0d0
        enddo; enddo; enddo
      endif

      do q=3,4
      do hel=-1,1,2

      do pol1 = 1,2
      do pol2 = 1,2
        pol(1) = pol1
        pol(2) = pol2
        amp(:,pol1,pol2)=0d0

c g(g1,g2) amplitudes
        amp(1,pol1,pol2) = +s1r(psi(1,hel,-id(3),ikin),epsGG(0,pol1,pol2),.true.,hel,PsiV(1,hel,q,id(1),3,ikin))
     &           +s1r(psiV(1,hel,q,-id(3),3,ikin),epsGG(0,pol1,pol2),.true.,hel,Psi(1,hel,id(1),ikin))
!         print*,amp(1,pol1,pol2)
        if(lAZ(0,1)) then
          amp(1,pol1,pol2) = amp(1,pol1,pol2)
     &           +s1r(psiV(1,hel,q,-id(3),1,ikin),epsGG(0,pol1,pol2),.true.,hel,PsiV(1,hel,q,id(1),2,ikin))
     &           +s1r(psiV(1,hel,q,-id(3),2,ikin),epsGG(0,pol1,pol2),.true.,hel,PsiV(1,hel,q,id(1),1,ikin))
        endif
!         print*,PsiV(:,hel,q,:,1,ikin)
!         print*,PsiV(:,hel,q,:,2,ikin)
!         print*,PsiV(:,hel,q,:,3,ikin)
!         print*,amp(1,pol1,pol2)
!         pause
        amp(2,pol1,pol2) = -amp(1,pol1,pol2)

        do glu1=1,2
          glu2 = 3-glu1
          idg1=id(2*glu1)
          idg2=id(2*glu2)
c VVgg, ggVV, g(WV)g amplitudes
          amp(glu1,pol1,pol2) = amp(glu1,pol1,pol2)
     &     + s1r(PsiG(1,hel,-id(3),idg1,pol(glu1),ikin),epsG(0,idg2,pol(glu2),ikin),.false.,hel,PsiV(1,hel,q,id(1),3,ikin))
     &     + s1c(PsiG(1,hel,-id(3),idg1,pol(glu1),ikin),epsV(0,hel,q,3,lephel,lepkin),.true.,hel,
     &           PsiG(1,hel,id(1),idg2,pol(glu2),ikin))
     &     + s1r(PsiV(1,hel,q,-id(3),3,ikin),epsG(0,idg1,pol(glu1),ikin),.false.,hel,PsiG(1,hel,id(1),idg2,pol(glu2),ikin))
! ! c all other amplitudes
        if(lAZ(0,1)) then
        do iV=1,2
        iV2=3-iV
          amp(glu1,pol1,pol2) = amp(glu1,pol1,pol2)
     &     + s1r(PsiVGlocal(1,hel,q,iV,idg1,pol(glu1)),epsG(0,idg2,pol(glu2),ikin),.false.,hel,PsiV(1,hel,q,id(1),iV2,ikin))
     &     + s1c(PsiVGlocal(1,hel,q,iV,idg1,pol(glu1)),epsV(0,hel,q,iV2,lephel,lepkin),.true.,hel,
     &                           PsiG(1,hel,id(1),idg2,pol(glu2),ikin))
        enddo !iV
        endif
        enddo !glu1
        
        
c multiply coupl. const. and get |amp|^2
        amp(1,pol1,pol2) = amp(1,pol1,pol2)*GG(1)**2
        amp(2,pol1,pol2) = amp(2,pol1,pol2)*GG(1)**2
!         print*,amp(:,pol1,pol2)
!       print*,amp(1,pol1,pol2),amp(2,pol1,pol2)
        if (lnlo) then
          resDiag(0,pol1,pol2,q) = resDiag(0,pol1,pol2,q) 
     &                   + dble(16d0*(amp(1,pol1,pol2)*dconjg(amp(1,pol1,pol2))+amp(2,pol1,pol2)*dconjg(amp(2,pol1,pol2)))
     &                   -4d0*amp(1,pol1,pol2)*dconjg(amp(2,pol1,pol2)))/3d0
          resDiag(1,pol1,pol2,q) = resDiag(1,pol1,pol2,q) 
     &                   + dble(amp(1,pol1,pol2)*dconjg(amp(1,pol1,pol2))-8d0*amp(2,pol1,pol2)*dconjg(amp(2,pol1,pol2))
     &                   +2d0*amp(1,pol1,pol2)*dconjg(amp(2,pol1,pol2)))
          resDiag(2,pol1,pol2,q) = resDiag(2,pol1,pol2,q) 
     &                   + dble(-amp(1,pol1,pol2)*dconjg(amp(1,pol1,pol2))-amp(2,pol1,pol2)*dconjg(amp(2,pol1,pol2))
     &                   -20d0*amp(1,pol1,pol2)*dconjg(amp(2,pol1,pol2)))/9d0
        else
          res(0,q) = res(0,q)+ dble(16d0*(amp(1,pol1,pol2)*dconjg(amp(1,pol1,pol2))+amp(2,pol1,pol2)*dconjg(amp(2,pol1,pol2)))
     &                              -4d0*amp(1,pol1,pol2)*dconjg(amp(2,pol1,pol2)))/3d0
          res(1,q) = res(1,q)+dble(amp(1,pol1,pol2)*dconjg(amp(1,pol1,pol2)))
          res(2,q) = res(2,q)+dble(amp(2,pol1,pol2)*dconjg(amp(2,pol1,pol2)))
        endif
      enddo; enddo ! pol summation
      ! return <1|Tij|2>+<2|Tij|1>, pol1 fixed
      if (lnlo) then
        do pol1 = 1,2
          resDiag(0,pol1,0,q) = resDiag(0,pol1,0,q)
     &                              +2d0*dble(16d0*(amp(1,pol1,1)*dconjg(amp(1,pol1,2))+amp(2,pol1,1)*dconjg(amp(2,pol1,2)))
     &                              -2d0*amp(1,pol1,1)*dconjg(amp(2,pol1,2))-2d0*amp(2,pol1,1)*dconjg(amp(1,pol1,2)))/3d0
          resDiag(1,pol1,0,q) = resDiag(1,pol1,0,q)
     &                              + 2d0*dble(amp(1,pol1,1)*dconjg(amp(1,pol1,2))-8d0*amp(2,pol1,1)*dconjg(amp(2,pol1,2))
     &                              +amp(1,pol1,1)*dconjg(amp(2,pol1,2))+amp(2,pol1,1)*dconjg(amp(1,pol1,2)))
          resDiag(2,pol1,0,q) = resDiag(2,pol1,0,q)
     &                              + 2d0*dble(-amp(1,pol1,1)*dconjg(amp(1,pol1,2))-amp(2,pol1,1)*dconjg(amp(2,pol1,2))
     &                              -10d0*amp(1,pol1,1)*dconjg(amp(2,pol1,2))-10d0*amp(2,pol1,1)*dconjg(amp(1,pol1,2)))/9d0
          resDiag(0,0,pol1,q) = resDiag(0,0,pol1,q)
     &                              + 2d0*dble(16d0*(amp(1,1,pol1)*dconjg(amp(1,2,pol1))+amp(2,1,pol1)*dconjg(amp(2,2,pol1)))
     &                              -2d0*amp(1,1,pol1)*dconjg(amp(2,2,pol1))-2d0*amp(2,1,pol1)*dconjg(amp(1,2,pol1)))/3d0
          resDiag(1,0,pol1,q) = resDiag(1,0,pol1,q)
     &                              + 2d0*dble(amp(1,1,pol1)*dconjg(amp(1,2,pol1))-8d0*amp(2,1,pol1)*dconjg(amp(2,2,pol1))
     &                              +amp(1,1,pol1)*dconjg(amp(2,2,pol1))+amp(2,1,pol1)*dconjg(amp(1,2,pol1)))
          resDiag(2,0,pol1,q) = resDiag(2,0,pol1,q)
     &                              + 2d0*dble(-amp(1,1,pol1)*dconjg(amp(1,2,pol1))-amp(2,1,pol1)*dconjg(amp(2,2,pol1))
     &                              -10d0*amp(1,1,pol1)*dconjg(amp(2,2,pol1))-10d0*amp(2,1,pol1)*dconjg(amp(1,2,pol1)))/9d0
        enddo
      endif

      enddo; enddo ! q, hel

      do q=3,4
      if(lnlo) then
        do pol1=1,2; do pol2=1,2
            resDiag(0,0,0,q)=resDiag(0,0,0,q)+resDiag(0,pol1,pol2,q)
            resDiag(1,0,0,q)=resDiag(1,0,0,q)+resDiag(1,pol1,pol2,q)
            resDiag(2,0,0,q)=resDiag(2,0,0,q)+resDiag(2,pol1,pol2,q)
        enddo; enddo
      else
        res(1,q) = res(1,q)/(res(1,q)+res(2,q))*res(0,q)
        res(2,q) = res(0,q)-res(1,q)
      endif
      enddo

      end


