c ************************************************************************************
c Author: Matthias Kerner
c Date: 2012/13
c ************************************************************************************

! #define NINH_COMPARE

      real*8 function m2s_QCDwz2j(
     &                   xi,  !in:  Feynman x parameters of incoming fermions
     &                   p,   !in:  momenta associated with external fermions
     &                   v,     !in:  momenta associated with Z/H dec. fermions
     &                   rn3,
     &                   Nbos, ! 32 = WpZ, 42 = WmZ
     &                   nlo
     &                        )

      use globalvars, only: ldoblha

      implicit none
c
c declare input/output variables
c
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/loops/Utilities/GaugeTest.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"
#include "common.inc"
      integer nbos
      real*8 xi(nx), p(0:3,max_p,max_kin), v(0:3,max_v,max_kin), m2s,rn3(3)
      integer nlo
c
c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"
      real*8  alphas, xmt, alfa, xmz, xmw,s2w, xmh
      common /bkopin/   alphas, xmt, alfa, xmz, xmw, s2w, xmh


c
c declare local variables
c
      real*8 q_sf

      integer init/0/, I
      save init
c
c declare variables for summation process
c
      integer FSIGN(4),hel(2)
      integer id(4) ! id(feyn index) = ps index
c
c store contributions from subprocess ip in res(ip) where
      real*8 res(maxNumSubProc)
      real*8 weight, rn 
c
c declare parton distribution variables
c
      real*8 pdf(-6:6,2,2)
      real*8 xiFinColl(2,2)
      
!       real*8 ampborntemp(0:2,3:4,0:2)

!       integer temphel(9)
!       logical lokt(16)


      complex*16 epsW(0:3,7),epsWV(0:3,-1:1,7) 
      complex*16 epsV(0:3,-1:1,3:4,-1:1,7) ! ( , mu-hel, u/d , quark-hel), coupl. const to quarks included
      complex*16 epsA(0:3,-1:1),epsZ(0:3,-1:1)
      real*8 momV(0:4,7),momW(0:4,7),momWV(0:4,7)
      common /qcdwznew/ epsW, epsWV, epsV,epsA,epsZ,momV,momW,momWV
 

c local variables
      integer FL4qcdVV, ires_min,ires_max
      save ires_min, ires_max
      external FL4qcdVV
      logical ldebug,allsub
      parameter (ldebug = .false., allsub=.false.)


      double precision  clr, xm2, xmg, b
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6)


      if ( (init .eq. 0) .and. (.not. ldoblha) ) then

         write(6,*) " "
        if(nbos .eq. 32) then
         write(6,*) "W+Z jj amplitude square information:"
        elseif(nbos.eq.42) then
         write(6,*) "W-Z jj amplitude square information:"
        elseif(nbos.eq.31) then
         write(6,*) "W+A jj amplitude square information:"
        elseif(nbos.eq.41) then
         write(6,*) "W-A jj amplitude square information:"
	elseif(nbos.eq.3) then
         write(6,*) "W+ jj amplitude square information:"
	elseif(nbos.eq.4) then
         write(6,*) "W- jj amplitude square information:"
        endif
         if(nfl.eq.5) then
          write(6,*) "external b quarks are included"
         else
          write(6,*) "external b quarks are excluded"
         endif
         write(6,*) "alpha parameter of dipoles: ",alphaDip
         write(6,*) "precision for gauge test:   ",GaugeAccuracy
         write(6,*) "-----------------------------------------------"
         
      else if ( (init .eq. 0) .and. ldoblha ) then
        init = 1
      endif
      
      if(nlo.ne.0) then
          GaugeAccuracyBoxL=GaugeAccuracy
          GaugeAccuracyBoxLnA=GaugeAccuracy
          GaugeAccuracyPenL=GaugeAccuracy
          GaugeAccuracyPenLnA=GaugeAccuracy
          GaugeAccuracyHexL=GaugeAccuracy
          GaugeAccuracyHexLnA=GaugeAccuracy
      endif

      if ( ldebug ) then
 10      format( " p(", i1, ") = (/ ", 3(f23.17,"d0,",2x), f23.17, "d0 /)" )
 20      format( " v(", i1, ") = (/ ", 3(f23.17,"d0,",2x), f23.17, "d0 /)" )
c
         write(6,*) " "
         do i = 1, n_v
            write(6,20) i, v(0,i,1), v(1,i,1), v(2,i,1), v(3,i,1)
         end do
         print*," parton momenta "
         do i = 1, n_p
            write(6,10) i, p(0,i,1), p(1,i,1), p(2,i,1), p(3,i,1)
         end do
         print*, "momentum conservation"
!          print*,(sum(p(mu,1:2,1))-sum(p(mu,3:4,1))-sum(v(mu,1:4,1)),mu=0,3)
c
         write(6,*) " "
         write(6,*) "xi(1), xi(2) =", xi
         write(6,*) "mu_f =",sqrt(mufsq(1,1)),sqrt(mufsq(2,1))
      end if
      

!        p(:,1,1) = (/514.089409873829    d0,   0.000000000000000d+000 , 0.000000000000000d+000 ,514.089409873829d0    /)
!        p(:,2,1) = (/   984.032214635719 d0,      0.000000000000000d+000 , 0.000000000000000d+000 ,-984.032214635719     d0 /)
!        p(:,3,1) = (/  500.670758477960  d0,     -158.175146947933   d0,     67.7843228513764  d0,   470.167115883784  d0 /)   
!        p(:,4,1) = (/  237.450379690683  d0,      2.60428596690189   d0,     25.2218642991637  d0,   -236.092689575782 d0 /)    
!        v(:,1,1) = (/  257.848330711510  d0,     -31.1971844337297   d0,    -36.4614559575744  d0,   -253.343757696087 d0 /)    
!        v(:,2,1) = (/  79.0088043361969  d0,      28.7776467124131   d0,     11.3233402353580  d0,   -72.7050216851994 d0 /)    
!        v(:,3,1) = (/  247.603696525318  d0,      60.7609367826238   d0,    -61.7947919972497  d0,   -231.942024601658 d0 /)    
!        v(:,4,1) = (/  175.539654767881  d0,      97.2294619197244   d0,    -6.07327943107390  d0,   -146.026427086948 d0 /)    
!        rn3(1)= 0.905948162078857d0
!        rn3(2)=0.876931428909302d0
!        rn3(3)=5.415290594100952E-002


      if (.not.ldoblha) then
c     
c call PDF subroutine in order to determine parton
c distributions in the incoming (anti)protons.
c
      q_sf = sqrt(mufsq(1,1))
      call pdfproton( xi(1), q_sf, pdf(-6,1,1) )        ! f_a(x1)

      q_sf = sqrt(mufsq(2,1))
      call pdfproton( xi(2), q_sf, pdf(-6,2,1) )        ! f_a(x2)
      do i=-5,5
        pdf(i,1,1) = pdf(i,1,1)/xi(1)
        pdf(i,2,1) = pdf(i,2,1)/xi(2)
      enddo

c pdfs for fin. coll.
      if(nlo.eq.3) then
	xiFinColl(1,1) = xi(1)
        xiFinColl(2,1) = xi(2)
	xiFinColl(1,2) = xi(1)/((1d0-xi(1))*rn3(3)+xi(1)) ! transform rn el [0;1] to x el [z;1], get xa = xxa/x, multiply jacobian to pdf
        xiFinColl(2,2) = xi(2)/((1d0-xi(2))*rn3(3)+xi(2))

	q_sf = sqrt(mufsq(1,1))
	call pdfproton( xiFinColl(1,2), q_sf, pdf(-6,1,2) )        ! f_a(x1)

	q_sf = sqrt(mufsq(2,1))
	call pdfproton( xiFinColl(2,2), q_sf, pdf(-6,2,2) )        ! f_a(x2)
        do i=-5,5
          pdf(i,1,2) = pdf(i,1,2)/xiFinColl(1,2)
          pdf(i,2,2) = pdf(i,2,2)/xiFinColl(2,2)
        enddo
        call CScolAmpQNewReset
        call CScolAmpGNewReset
      endif

      endif !ldoblha 

!       pdf=1d0
!       als=1d0
!       print*,rn3
!       rn3=0.5d0
!       rn3(3)=0.995772600173950d0

      call calcEWparts(v,nbos,rn3(2))

c reset subprocess counter
      ires_min = FL4QCDVV(id,id,-sub_number) +1
      
c random helicity
      rn=rn3(1) ! don't change rn3, so that one can reevaluate the complete ps point 
      rn3(1) = 2d0*rn3(1)
      hel(1) = (-1)**INT(rn3(1))
      rn3(1) = rn3(1) - INT(rn3(1))
      rn3(1) = 2d0*rn3(1)
      hel(2) = (-1)**INT(rn3(1))
      rn3(1) = rn


#ifdef NINH_COMPARE


!       als(1,1) = 0.129808134118326d0
!       als(2,1) = 0.129808134118326d0
      GG(1)= sqrt(als(1,1)*4d0*pi)
      GG(2) = GG(1)
      G = -GG(1)
      alphas = als(1,1)

c point 1
!        p(:,1,1) = (/ 500d0,0d0,0d0, 500d0 /)
!        p(:,2,1) = (/ 500d0,0d0,0d0,-500d0 /)
!        p(:,3,1) = (/ 88.1907917270800d0,-71.8338387342512d0,-26.8500535430029d0,-43.5498562901036d0 /)
!        p(:,4,1) = (/ 457.182458695087d0, 404.255345805865d0, 0.00000000000000d0,-213.526148108059d0 /)
!        v(:,1,1) = (/ 100.788423998349d0,-45.8888076188476d0, 49.1583762549171d0, 75.0731496033797d0 /)
!        v(:,2,1) = (/ 199.936105394529d0,-187.698928314700d0,-25.1678128082751d0, 64.1103716125017d0 /)
!        v(:,3,1) = (/ 149.590646659259d0,-96.4973730388405d0, 2.67112519191297d0, 114.273722503239d0 /)
!        v(:,4,1) = (/ 4.31157352569544d0,-2.33639809922604d0, 0.18836490444800d0, 3.61876067904234d0 /)

c point 2
!       v(:,1,1) = (/   128.28517384798186640d0,     12.81205307979342223d0,    -69.19348313631607539d0,    107.26229074820207643d0 /)
!       v(:,2,1) = (/   556.02626579167394993d0,   -111.18103223107226540d0,   -276.99462728273675793d0,    469.12467722212636545d0 /)
!       v(:,3,1) = (/   181.93738317191755982d0,      2.47586265168042274d0,    181.73731791511434608d0,      8.16264519600288097d0 /)
!       v(:,4,1) = (/    50.32502148264003949d0,     10.99910781031214846d0,     26.51191270377822917d0,     41.33697980492176782d0 /)
! c 	parton momenta 
!       p(:,1,1) = (/   817.04583297571480216d0,      0.00000000000000000d0,      0.00000000000000000d0,    817.04583297571480216d0 /)
!       p(:,2,1) = (/   507.68698916035094726d0,      0.00000000000000000d0,      0.00000000000000000d0,   -507.68698916035094726d0 /)
!       p(:,3,1) = (/   154.60391740584950071d0,    -51.30711878197600839d0,      9.48511674338542221d0,   -145.53344426588444094d0 /)
!       p(:,4,1) = (/   253.55506043600280464d0,    136.20112747126228214d0,    128.45376305677484652d0,   -170.99430489000485522d0 /)

c WA paper
!       print*,GWF
!       print*,GWWA
!       print*,GAL
!       print*,clr(3:4,1,-1:1:2)
!       print*
      als=1d0
      alphas=1d0
      G=sqrt(4d0*pi)
      gg(1)=-G
      gg(2)=-G
      WMASS=80.385d0
      ZMASS=91.1876d0
      GWF(1)=-sqrt(4d0*pi)/sqrt(1d0-(WMASS/ZMASS)**2)/sqrt(2d0)
      GWF(2)=0d0
      WWIDTH=2.097672973437829d0
      ZWIDTH=2.508419926825088d0
!       WWIDTH=2.09761d0
!       ZWIDTH=2.5089d0 
      GWWA=sqrt(4d0*pi)
      GAL=GWWA
      CLR(3,1,-1:1:2) = +GWWA*2d0/3d0
      CLR(4,1,-1:1:2) = -GWWA*1d0/3d0
!       print*,GWF
!       print*,GWWA
!       print*,GAL
!       print*,clr(3:4,1,-1:1:2)
      GAU(1:2) = -clr(3,1,1)
      GAD(1:2) = -clr(4,1,1)
      p(0:3,1,1)=  (/  32.0772251055223d0, 0d0, 0d0, 32.0772251055223d0/)
      p(0:3,2,1) =   (/  2801.69305619768d0, 0d0, 0d0, -2801.69305619768d0/)
      p(0:3,3,1) =   (/  226.525314156010d0, -10.2177083492279d0, -1.251308382450315D-15, -226.294755550298d0/)
      p(0:3,4,1) =   (/  327.281588297290d0, -6.48554750244653d0, -10.1061447270513d0, -327.061219882068d0/)
      v(0:3,1,1) =   (/  646.824307052136d0, 36.0746355875450d0, -26.0379256562231d0, -645.292438579767d0/)
      v(0:3,2,1) =   (/  1598.85193997112d0, -2.88431497177613d0, 24.4490976584709d0, -1598.66239347157d0/)
      v(0:3,3,1) =   (/  34.2871318266438d0, -16.4870647640944d0, 11.6949727248035d0, 27.6949763915464d0/)

!       print*,(sum(p(mu,3:4,1))+sum(v(mu,1:3,1))-sum(p(mu,1:2,1)),mu=0,3)

      xi(1)= 0.468022827150478d0
      xi(2)= 0.786655185286678d0
      q_sf = sqrt(mufsq(1,1))
      call pdfproton( xi(1), q_sf, pdf(-6,1,1) )        ! f_a(x1)
      pdf(:,1,1) = pdf(:,1,1)/xi(1)

      q_sf = sqrt(mufsq(2,1))
      call pdfproton( xi(2), q_sf, pdf(-6,2,1) )        ! f_a(x2)
      pdf(:,2,1) = pdf(:,2,1)/xi(2)

        do i=-6,6
          pdf(i,1,1) =1d0!+0.01d0*i
          pdf(i,2,1) =1d0!+0.02d0*i
        enddo      
        
      
!       if(init.eq.0) print*,(p(mu,3,1)+p(mu,4,1)+sum(v(mu,1:4,1))-sum(p(mu,1:2,1)),mu=0,3)
        
      call calcEWparts(v,nbos,rn3(2))
!       call calcEWpartsOld(v(0,1,1))

!       nlo=1
c gg ->ub d

      id(1) = 3
      id(3) = 4
      id(2) =1
      id(4) = 2
      fsign=-1
      fsign(3)=1
      call getRes2Q2GWZ(p,xiFinColl,id,fsign,pdf,nbos,2,hel,rn3(3),res) 
      
c u db -> cb c
      id(1) = 1
      id(2) = 3
      id(3) = 2
      id(4) = 4
      fsign=1
      fsign(2:3) = -1
      call getRes4QWZ(res,p,xiFinColl,fsign,id,pdf,nbos,2,rn3(3) )

      print*,'mg'

      pMG(:,1:4) = p(:,1:4,1)
      pMG(:,5:7) = v(:,1:3,1)
      call Smatrixggudenua(pMG,rn3(3))
      print*,rn3(3)
      return

c u db -> uu
      
      
!       call getRes4QWZ(res,p,fsign,id,pdf,nbos,nlo,rn3(3))
!       goto 123

!       id(1) = 3
!       id(2) = 1
!       id(3) = 4
!       id(4) = 2
! 
!       fsign(1) = -1
!       fsign(2) = -1
!       fsign(3) =  1
!       fsign(4) = -1
!       call getRes2Q2GWZ(p,id,fsign,pdf,nbos,nlo,hel,rn3(3),res)
!       res = res/2d0
!       goto 123

#endif
! !       id(1) = 1
! !       id(2) = 2
! !       id(3) = 3
! !       id(4) = 4
! ! 
! !       fsign(1) =  1
! !       fsign(2) =  -1
! !       fsign(3) =  1
! !       fsign(4) =  1
! ! 
! !       lokt=.false.
! !       lokt(1) = .true.
! !       call precalcQCDVVjj(p,lokt,0,nbos)
! !       call ug2dgwzCSNew(P,1,id,fsign,nbos,0)
      
! !       pMG(:,1:4) = p(:,1:4,1)
! !       pMG(:,5:6) = v(:,1:2,1)
! ! !       call SMATRIX_UUUDEN(pMG,res(1))
! !       call SMATRIX_UDDDen(pMG,res(1))
! !       print*,res(1)
! !       call amp4qwBorn(p,id,fsign,nbos,0,ampborntemp)
! !       print*,ampborntemp(0,3,0),ampborntemp(1,3,1),ampborntemp(2,3,2)
! ! !       id(3) = 4
! ! !       id(4) = 3
! !       id(1) = 2
! !       id(2) = 1
! !       call amp4qwBorn(p,id,fsign,nbos,0,ampborntemp)
! !       print*,ampborntemp(0,3,0),ampborntemp(2,3,2),ampborntemp(1,3,1)
! !       call amp4qwVirt(P,id,fsign,nbos,0,1,rn3(3),ampborntemp)
! !       print*,ampborntemp(0,3,0),ampborntemp(2,3,2),ampborntemp(1,3,1)
! !       
! !       pause

      
! begin 2q 2g subprocesses
! Feynam diagram labels:

c   q1 -------------- q3
c            9
c            6
c            9
c   g2 96969696969696 g4

       if (ldoblha) then
c id is diagtoPhys in standard VBFNLO terminology
         id(blha_physdiagmap(1,blha_cursubproc,blha_curproc))=1
         id(blha_physdiagmap(2,blha_cursubproc,blha_curproc))=2
         id(blha_physdiagmap(3,blha_cursubproc,blha_curproc))=3
         id(blha_physdiagmap(4,blha_cursubproc,blha_curproc))=4
         fsign(1) = blha_fsign(1,blha_cursubproc,blha_curproc)
         fsign(2) = blha_fsign(2,blha_cursubproc,blha_curproc)
         fsign(3) = blha_fsign(3,blha_cursubproc,blha_curproc)
         fsign(4) = blha_fsign(4,blha_cursubproc,blha_curproc)

         if ((blha_idsubproc(blha_cursubproc,blha_curproc)/100) .eq. 1) then
           call getRes2Q2GWZ(p,xiFinColl,id,fsign,pdf,nbos,nlo,hel,rn3(3),res)
         else if ((blha_idsubproc(blha_cursubproc,blha_curproc)/100) .eq. 2) then
           call getRes4QWZ(res,p,xiFinColl,fsign,id,pdf,nbos,nlo,rn3(3) )
         endif
         return
       endif

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

      if(sub_number.eq.1.or.allsub) call getRes2Q2GWZ(p,xiFinColl,id,fsign,pdf,nbos,nlo,hel,rn3(3),res)

! switch initial state
      id(1)=2
      id(2)=1
      
      if(sub_number.eq.1.or.allsub) call getRes2Q2GWZ(p,xiFinColl,id,fsign,pdf,nbos,nlo,hel,rn3(3),res)

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

      if(sub_number.eq.2.or.allsub) call getRes2Q2GWZ(p,xiFinColl,id,fsign,pdf,nbos,nlo,hel,rn3(3),res)

! switch initial state
      id(2)=1
      id(3)=2

      if(sub_number.eq.2.or.allsub) call getRes2Q2GWZ(p,xiFinColl,id,fsign,pdf,nbos,nlo,hel,rn3(3),res)

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

      if(sub_number.eq.3.or.allsub) call getRes2Q2GWZ(p,xiFinColl,id,fsign,pdf,nbos,nlo,hel,rn3(3),res)

! switch initial state
      id(1)=2
      id(3)=1

      if(sub_number.eq.3.or.allsub) call getRes2Q2GWZ(p,xiFinColl,id,fsign,pdf,nbos,nlo,hel,rn3(3),res)

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
      if(sub_number.eq.4.or.allsub) call getRes2Q2GWZ(p,xiFinColl,id,fsign,pdf,nbos,nlo,hel,rn3(3),res)
      
c end of 2q2g processes
c begin 4q processes
c
c   q1 -------------- q3
c            9
c            6
c            9
c   q2 -------------- q4
c
c  W is radiated of 13 quark line
c
c*****  1.  ********
c**Diagram**********  q1 q2 ---> q3 q4 W Z   ********************** 
c**Diagram**********  q2 q1 ---> q3 q4 W Z   **********************
c**Phasespace*******  1   2        3   4       **********************

      id(1) = 1
      id(2) = 2
      id(3) = 4
      id(4) = 3

      fsign(1) =  1
      fsign(2) =  1
      fsign(3) =  1
      fsign(4) =  1
   

      if(sub_number.eq.5.or.allsub) call getRes4QWZ(res,p,xiFinColl,fsign,id,pdf,nbos,nlo,rn3(3) )
! switch initial state
      id(1) = 2
      id(2) = 1
      
      if(sub_number.eq.5.or.allsub) call getRes4QWZ(res,p,xiFinColl,fsign,id,pdf,nbos,nlo,rn3(3) )

! 
! c*****  2.  ********
! c**Diagram**********  q1 qb4 ---> q3 qb2 W Z   **********************
! c**Diagram**********  qb4 q1 ---> q3 qb2 W Z   **********************
! c**Phasespace*******  1   2        3   4
! 
      id(1) = 1
      id(2) = 3
      id(3) = 4
      id(4) = 2

      fsign(1) =  1
      fsign(2) = -1
      fsign(3) =  1
      fsign(4) = -1
      

      if(sub_number.eq.6.or.allsub) call getRes4QWZ(res,p,xiFinColl,fsign,id,pdf,nbos,nlo,rn3(3) )
! switch initial state
      id(4) = 1
      id(1) = 2

      if(sub_number.eq.6.or.allsub) call getRes4QWZ(res,p,xiFinColl,fsign,id,pdf,nbos,nlo,rn3(3) )
! 
! 
! c*****  3.  ********
! c**Diagram**********  qb3 q2 ---> qb1 q4 W Z   **********************
! c**Diagram**********  q2 qb3 ---> qb1 q4 W Z   **********************
! c**Phasespace*******  1   2      3   4
! 
      id(1) = 3
      id(2) = 2
      id(3) = 1
      id(4) = 4

      fsign(1) = -1
      fsign(2) =  1
      fsign(3) = -1
      fsign(4) =  1 

      if(sub_number.eq.7.or.allsub) call getRes4QWZ(res,p,xiFinColl,fsign,id,pdf,nbos,nlo,rn3(3) )
! switch initial state
      id(2) = 1
      id(3) = 2
      
      if(sub_number.eq.7.or.allsub) call getRes4QWZ(res,p,xiFinColl,fsign,id,pdf,nbos,nlo,rn3(3) )

c*****  4.  ********
c**Diagram**********  qb3 qb4 ---> qb1 qb2 W Z   **********************
c**Diagram**********  qb4 qb3 ---> qb1 qb2 W Z   **********************
c**Phasespace*******  1   2        3   4

      id(1) = 3
      id(2) = 4
      id(3) = 1
      id(4) = 2

      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = -1
      fsign(4) = -1

      if(sub_number.eq.8.or.allsub) call getRes4QWZ(res,p,xiFinColl,fsign,id,pdf,nbos,nlo,rn3(3) )
! switch initial state
      id(4) = 1
      id(3) = 2

      if(sub_number.eq.8.or.allsub) call getRes4QWZ(res,p,xiFinColl,fsign,id,pdf,nbos,nlo,rn3(3) )
c
c
c  end of subprocess evaluation
c
c
      continue

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

c select subprocess for lha output
      iprocess =0
      weight = 0d0
      i=ires_min
      do while((i.le.ires_max).and.(iprocess.eq.0))
         weight = weight+res(i)/m2s
         if(weight.ge.rn3(3)) then
            iprocess = i
         endif
         i = i+1
      enddo
      
*     scale for running alphas
      m2s = m2s*(als(1,1)/alphas)**2
      

      m2s_QCDwz2j =m2s
#ifdef NINH_COMPARE      
      print*,m2s,m2s/2d0
#endif

!       print*,m2s, hel,rn3(3)
!       pause
      END


      subroutine getRes2Q2GWZ(p,xiFinColl,id,fsign,pdf,Nbos,nlo,hel,rn,res)

      use globalvars, only: ldoblha

      implicit none
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"
      real*8 p(0:3,max_p,max_kin), pdf(-6:6,2,2), xiFinColl(2,2),rn! for gauge test
      integer id(4), fsign(4),nbos,nlo,hel(2)
      real*8 res(maxNumSubProc)  !out: subprocess cross-sections

      integer FL4qcdVV
      external FL4qcdVV
      
      real*8 amp(0:2), polcol, ampDiag(0:2,0:2,0:2)
      integer colstr,d,u,i,k
      LOGICAL ldebug
      parameter (ldebug=.false.)
      integer iflav(1:4)

      real*8 QCDVVjjFinCollBornKin
      external randomnumber


      complex*16 epsW(0:3,7),epsWV(0:3,-1:1,7) 
      complex*16 epsV(0:3,-1:1,3:4,-1:1,7) ! ( , mu-hel, u/d , quark-hel), coupl. const to quarks included
      complex*16 epsA(0:3,-1:1),epsZ(0:3,-1:1)
      real*8 momV(0:4,7),momW(0:4,7),momWV(0:4,7)
      common /qcdwznew/ epsW, epsWV, epsV,epsA,epsZ,momV,momW,momWV


      if(nlo.eq.0) then
!         hel sum
	if(nbos.gt.10) then
#ifdef WITH_QCDVV
          if (ldoblha) then
            call ug2dgwzCSnew(p,1,id,fsign,nbos,0,.false.,ampDiag)
          else
            call ug2dgwzBornnew(p,id,fsign,nbos,0,amp)
          endif
#endif
	else
          if (ldoblha) then
            call ug2dgwCS(p,1,id,fsign,nbos,ampDiag)
          else
	    call ug2dgwBorn(p,id,fsign,nbos,amp)
	  endif
	endif


#ifdef NINH_COMPARE
	temp2= 0d0
	call ug2dgwzBornnew(p,id,fsign,nbos,1,amp)
	temp2 = amp(0)
	call ug2dgwzBornnew(p,id,fsign,nbos,-1,amp)
	temp2 = amp(0)+temp2
	amp(0) = temp2
	amp(1) = temp2
	amp(2) = 0d0
        print*,id
        print*,fsign
        print*,nbos
        print*,amp
#endif

      elseif(nlo.le.2) then ! use fincoll to calc nlo.eq.1
	
! ! ! ! ! 	  print*,hel,rn,id,fsign,nbos
! ! ! ! ! 	  call ug2dgwVirt(P,id,fsign,nbos,2,hel(1),hel(2),rn,amp)
! ! ! ! ! 	  print*,amp
! ! ! ! ! 	  call ug2dgwzVirt(P,id,fsign,nbos,2,0,hel(1),hel(2),rn,amp)
! ! ! ! ! 	  print*,amp
! ! ! ! ! 	  pause

#ifdef NINH_COMPARE
	temp2= 0d0
	call ug2dgwzVirt(P,id,fsign,nbos,nlo,-1,-1,-1,rn,amp)
	temp2 = temp2+amp(0)
	call ug2dgwzVirt(P,id,fsign,nbos,nlo,-1,-1,+1,rn,amp)
	temp2 = temp2+amp(0)
	call ug2dgwzVirt(P,id,fsign,nbos,nlo,-1,+1,-1,rn,amp)
	temp2 = temp2+amp(0)
	call ug2dgwzVirt(P,id,fsign,nbos,nlo,-1,+1,+1,rn,amp)
	temp2 = temp2+amp(0)
	call ug2dgwzVirt(P,id,fsign,nbos,nlo, 1,-1,-1,rn,amp)
	temp2 = temp2+amp(0)
	call ug2dgwzVirt(P,id,fsign,nbos,nlo, 1,-1,+1,rn,amp)
	temp2 = temp2+amp(0)
	call ug2dgwzVirt(P,id,fsign,nbos,nlo, 1,+1,-1,rn,amp)
	temp2 = temp2+amp(0)
	call ug2dgwzVirt(P,id,fsign,nbos,nlo, 1,+1,+1,rn,amp)
	temp2 = temp2+amp(0)
	amp(0) = temp2
	amp(1) = temp2
	amp(2) = 0d0
#else
      ! rnd hel
#ifdef WITH_NLO
	if(nbos.gt.10) then
#ifdef WITH_QCDVV
          if (ldoblha) then
            do i=-1,1,2
            do k=-1,1,2
              call ug2dgwzVirt(P,id,fsign,nbos,1,0,i,k,rn,amp)
            enddo
            enddo
          else
	    call ug2dgwzVirt(P,id,fsign,nbos,nlo,0,hel(1),hel(2),rn,amp)
          endif
          do i=0,2
            amp(i)=amp(i)*4
          enddo

#endif
	else
          if (ldoblha) then
            do i=-1,1,2
            do k=-1,1,2
              call ug2dgwVirt(P,id,fsign,nbos,1,i,k,amp)
            enddo
            enddo
          else
	    call ug2dgwVirt(P,id,fsign,nbos,nlo,hel(1),hel(2),amp)
          endif
          do i=0,2
            amp(i)=amp(i)*4
          enddo
	endif
#endif
#endif

      endif

      if (ldoblha) then
        return
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
      do i = 1,2 ! family
	u = 2*i             !d,s quarks
	d = u - 1           !u,c quarks
	if(Nbos.eq.42.or.nbos.eq.4.or.nbos.eq.41) then ! W- case
	    iflav(id(1)) = d*fsign(1)
	    iflav(id(3)) = u*fsign(3)
	elseif(Nbos.eq.32.or.nbos.eq.3.or.nbos.eq.31) then ! w+ case
	    iflav(id(1)) = u*fsign(1)
	    iflav(id(3)) = d*fsign(3)
	else
	    stop
	endif
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
	    res(k) = amp(colstr)*POLCOL*
     &		     pdf(sign1*iflav(1),1,1)*pdf(sign2*iflav(2),2,1)
	  endif
	  if(ldebug .and.colstr.eq.1)  print*,k, iflav
	enddo
	
      enddo
      end


      subroutine getRes4QWZ(res,p,xiFinColl,fsign,id,pdf,Nbos,nlo,rn)

      use globalvars, only: ldoblha

      implicit none
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"
#include "common.inc"
      real*8 res(1:maxNumSubProc) !out: subprocess cross-sections
      real*8 p(0:3,max_p,max_kin), xiFinColl(2,2),rn
      integer fsign(4), id(4)
      real*8 pdf(-6:6,2,2)
      integer nbos,nlo
      
      integer FL4qcdVV,di
      external FL4qcdVV,di

      real*8 POLCOL, polcol0    ! ave over intial polarizations and colors
      integer colstr,d,u,q,j,k,temp
      integer iflav(1:4)
      integer family,flav
      logical ldebug
      parameter (ldebug = .false.)
      real*8 amp(0:2,3:4,0:2),QCDVVjjFinCollBornKin


      complex*16 epsW(0:3,7),epsWV(0:3,-1:1,7) 
      complex*16 epsV(0:3,-1:1,3:4,-1:1,7) ! ( , mu-hel, u/d , quark-hel), coupl. const to quarks included
      complex*16 epsA(0:3,-1:1),epsZ(0:3,-1:1)
      real*8 momV(0:4,7),momW(0:4,7),momWV(0:4,7)
      common /qcdwznew/ epsW, epsWV, epsV,epsA,epsZ,momV,momW,momWV


      polcol0 = 1.0d0/(4.0d0 * 9.0d0)
      
      if(nlo.eq.0) then
	if(nbos.gt.10) then
#ifdef WITH_QCDVV
         if (ldoblha) then
          call amp4qwzCS(p,1,id,fsign,nbos,0,amp)
         else
	  call amp4qwzBorn(p,id,fsign,nbos,0,0,amp)
         endif
#ifdef NINH_COMPARE
!           call amp4qwzBorn(p,id,fsign,nbos,0,1,amp)
!           call amp4qwzBorn(p,id,fsign,nbos,0,-1,amptemp)
          call amp4qwzVirt(P,id,fsign,nbos,0,1,nlo,rn,amp)
          call amp4qwzVirt(P,id,fsign,nbos,0,-1,nlo,rn,amptemp)
          amp=amp+amptemp
#endif
#endif
	else ! W+ or W-
         if (ldoblha) then
          call amp4qwCS(p,1,id,fsign,nbos,0,amp)
         else
	  call amp4qwBorn(p,id,fsign,nbos,0,amp)
	 endif
	endif
      elseif(nlo.le.2) then
#ifdef WITH_NLO
	if(nbos.gt.10) then
#ifdef WITH_QCDVV
	  call amp4qwzVirt(P,id,fsign,nbos,0,0,nlo,rn,amp)
#endif
	else
	  call amp4qwVirt(P,id,fsign,nbos,nlo,rn,amp)
	endif
#endif

! ! ! ! ! 	  call amp4qwzBorn(p,id,fsign,nbos,0,0,amp)

#ifdef NINH_COMPARE
	call amp4qwzVirt(P,id,fsign,nbos,0,1,nlo,rn,amp)
	call amp4qwzVirt(P,id,fsign,nbos,0,-1,nlo,rn,amptemp)
	amp = amp+amptemp
#endif
      endif
	
      if (ldoblha) then
        return
      endif

c sum over all quark combinations
      do j = 1,2  ! W emitting Quark family
	u = 2*j              !u,c quarks
	d = u - 1            !d,s quarks
	do flav = 1,2   ! flavor of quark line not emitting W
          do family =0,1+nfl/5*(-flav+2) ! 2 families + b Quark if nf == 5
            q = flav + 2*family  ! u,d,s,c quarks


c fill flavor information for pdf and lha output
            if(Nbos.eq.42.or.nbos.eq.4.or.nbos.eq.41) then  ! w- case
               iflav(id(1)) = d*fsign(1)
               iflav(id(3)) = u*fsign(3)
            elseif(Nbos.eq.32.or.nbos.eq.3.or.nbos.eq.31) then ! w+ case
               iflav(id(1)) = u*fsign(1)
               iflav(id(3)) = d*fsign(3)
            else
               stop
            endif

            iflav(id(2)) = q * fsign(2)
            iflav(id(4)) = q * fsign(4)

c  pick polcol
            if(iflav(3).eq.iflav(4)) then
               polcol = polcol0/2.0d0
            else
               polcol = polcol0
            endif
            

c get m2s*pdf
            if(j-1.eq.family) then  ! same q families
	      if(iflav(1).ne.iflav(2).or.di(id,1).lt.di(id,2)) then ! see my table of subprocesses
                do colstr = 1,2  ! two color flows
                  k = FL4qcdVV(iflav,id,colstr)
		  if(nlo.eq.3) then
		    res(k) = 0d0
		    if(colstr.eq.1) res(k) = QCDVVjjFinCollBornKin(p,xiFinColl,id,fsign,iflav,0,nbos,pdf) ! born is included here
		  else
		    res(k) = POLCOL * amp(colstr,5-flav,0)*
     $                     pdf(sign1*iflav(1),1,1)*pdf(sign2*iflav(2),2,1)
		  endif
		  if(ldebug.and.colstr.eq.1) print*,k, iflav, res(k),amp(1,5-flav,0),amp(2,5-flav,0),amp(0,5-flav,0)
               enddo
	      endif
	    else
	      k = FL4qcdVV(iflav,id,1)
	      if(nlo.eq.3) then
		res(k) = QCDVVjjFinCollBornKin(p,xiFinColl,id,fsign,iflav,1,nbos,pdf) ! born is included here
	      else
		res(k) = POLCOL * amp(1,5-flav,1)
     $		      *pdf(sign1*iflav(1),1,1)*pdf(sign2*iflav(2),2,1)
	      endif
	      if(ldebug) print*,k, iflav,res(k)

	      if(iflav(1)*iflav(2).lt.0) then ! include colstr 2 for q qb initial states
		if(mod(iflav(id(1))+iflav(id(2)),2).eq.0) then ! swap particle 1 and 2; change sign if needed
		  temp         = sign(iflav(id(1)),fsign(2))
		  iflav(id(1)) = sign(iflav(id(2)),fsign(1))
		  iflav(id(2)) = temp
		else
		  temp         = sign(iflav(id(3)),fsign(4))
		  iflav(id(3)) = sign(iflav(id(4)),fsign(3))
		  iflav(id(4)) = temp
		endif
		k = FL4qcdVV(iflav,id,2)
		if(nlo.eq.3) then
		  res(k) = QCDVVjjFinCollBornKin(p,xiFinColl,id,fsign,iflav,2,nbos,pdf) ! born is included here
		else
		  res(k) = POLCOL * amp(2,5-flav,2)
     $		      *pdf(sign1*iflav(1),1,1)*pdf(sign2*iflav(2),2,1)
		endif
		if(ldebug) print*,k, iflav,res(k)
	      endif

            endif
	  enddo
	enddo
      enddo
      return
      end

      integer function FL4qcdVV(iflav,id,colstr)
      Implicit none
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/lha.inc"
c     COMMOM BLOCK FOR DIAGRAMMATIC LABELS
      integer id(4),iflav(4),colstr
      integer listp,nmin
      save listp,nmin
      logical lfill
      integer lastsub
      save lastsub
c      integer,save :: listp ! this is the list position
c      integer,save :: numresets
      data lastsub /0/
      data listp /0/

      if (colstr.le.-1) then  ! we are resetting, use -sub_number to identify process range
        lfill=colstr.ne.lastsub ! fill parton information for new subprocesses
        if(lfill) then
          lastsub=colstr
          nmin=listp
        endif
        listp=nmin ! set position to last subprocess
      else
         listp = listp + 1
         if(lfill.and.(lha.or.hepmc)) then
             call fillColoredPartonsQCDVV(iflav,id,colstr,listp)
             helicity(1)=9
             helicity(2)=9
             helicity(3)=9
             helicity(4)=9
         endif                  ! numresets.eq.0
      endif
!       print*,'fl:',listp, iflav,colstr,"min,fill",nmin,lfill
      numcolor = listp ! common block info
      FL4qcdVV   = listp ! outputed function
c
      end


	subroutine fillColoredPartonsQCDVV(iflav,id,colstr,listp1)
	implicit none
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/process.inc"

c          COMMOM BLOCK FOR DIAGRAMMATIC LABELS
	integer id(4)
	integer i,colstr,listp1,listp
	integer ngluon,npart
	parameter (npart = 4)
	integer iflav(npart),pi(1:2)
	integer col(1:2)
c
	listp = listp1 
	lnup(listp) = 4+3*N_lepton/2
c
c         lidprup(listposition)=listposition so we don't need to use it.
c         scalup is same for all subprocesses
c         lwgtup(listposition) is already stored in res(fl(....)) so
c         there's
c         no need to worry about it here.
c           write(*,*) 'writing to lhepup'
c

	listup(1,listp)=-1   !incoming quark or gluon
	listup(2,listp)=-1
	listup(3,listp)=1    !outgoin quark or gluon
	listup(4,listp)=1

c
c         The particle codes:
      lidup(1,listp)=iflav(1)
      lidup(2,listp)=iflav(2)
      lidup(3,listp)=iflav(3)
      lidup(4,listp)=iflav(4)

c
c         mothers:
      do i = 3,npart
	  lmothup(1,i,listp) = 1
	  lmothup(2,i,listp) = 2
      enddo
c
c     4Q and 2Q2G are the possibilities.
c     First, count the number of gluons to decide whether
c     4Q or 2Q2G.
      ngluon = 0
      do i = 1, npart
	  if(iflav(i).eq.21) then
	    ngluon = ngluon + 1
	  endif
      enddo
      
c
      do i = 1,2
	  col(i) = 0
      enddo
      if(ngluon.eq.0) then ! 4Q case
	  col(colstr) = 1 ! selects color structure
c
c
	  if(id(1).lt.3) then ! 1 and 2 are incoming pythia labels
	    licolup(1,id(1),listp) = 501
	    licolup(2,id(1),listp) = 0 ! anti-color
	  else
	    licolup(1,id(1),listp) = 0 ! color
	    licolup(2,id(1),listp) =  501 !anti color
	  endif
c
	  if(id(2).lt.3) then
	    licolup(1,id(2),listp) = 502
	    licolup(2,id(2),listp) = 0 !anti color
	  else
	    licolup(1,id(2),listp) = 0
	    licolup(2,id(2),listp) = 502 ! anti-color
	  endif
c
	  if(id(3).lt.3) then
	    licolup(1,id(3),listp) = 0
	    licolup(2,id(3),listp) = 502*col(1)+501*col(2) ! anti-color
	  else
	    licolup(1,id(3),listp) = 502*col(1)+501*col(2)
	    licolup(2,id(3),listp) = 0
	  endif
c
	  if(id(4).lt.3) then
	    licolup(1,id(4),listp) = 0 ! color
	    licolup(2,id(4),listp) = 501*col(1) +
     $               502*col(2)   ! anti-color
	  else
	    licolup(1,id(4),listp) = 501*col(1) +
     $               502*col(2)  ! color
	    licolup(2,id(4),listp) = 0 ! anti-color
	  endif
c
c
      else  ! 2Q2G case
c     The next case goes here.
c
c
c    q1 ----->--------------->--------q2
c             0             0
c             0             0
c             0             0
c           g_pi(1)       g_pi(2)
c
c        pi = (1),(12)   See HZ
c    colstr =  1,   2
c
	  select case(colstr)
	  case(1)
	    pi(1) = 4
	    pi(2) = 2
	  case(2)
	    pi(1) = 2
	    pi(2) = 4
	  endselect
c
	  if(id(1).lt.3) then  ! 1,2 are incoming
	    licolup(1,id(1),listp) = 501 ! color
	    licolup(2,id(1),listp) = 0 !anti-color
	  else
	    licolup(1,id(1),listp) = 0 ! color
	    licolup(2,id(1),listp) = 501 !anti-color
	  endif
c
	  if(id(3).lt.3) then  ! 1,2 are incoming
	    licolup(1,id(3),listp) = 0 ! color
	    licolup(2,id(3),listp) = 503 !anti-color
	  else
	    licolup(1,id(3),listp) = 503 ! color
	    licolup(2,id(3),listp) = 0 !anti-color
	  endif
c
	  if(id(pi(1)).lt.3) then
	    licolup(1,id(pi(1)),listp) = 502 ! color
	    licolup(2,id(pi(1)),listp) = 501 !anti-color
	  else
	    licolup(1,id(pi(1)),listp) = 501 ! color
	    licolup(2,id(pi(1)),listp) = 502 !anti-color
	  endif
	  if(id(pi(2)).lt.3) then
	    licolup(1,id(pi(2)),listp) = 503 ! color
	    licolup(2,id(pi(2)),listp) = 502 !anti-color
	  else
	    licolup(1,id(pi(2)),listp) = 502 ! color
	    licolup(2,id(pi(2)),listp) = 503 !anti-color
	  endif
c
      endif
      end ! fillColoredPartons
! 
      SUBROUTINE calcEWparts(v,bos,rn)
c     calculates currents
c     WW:  W+ > e+ ve
c     WZ:  Z  > mu+ mu-
c     WA:  A  > mu+ mu-
c     WWV: W+ > e+ ve mu+ mu-
c     2. Index = spin of mu+, mu-

      use globalvars, only: ldoblha

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"
#include "common.inc"

      real*8 v(0:3,max_v,max_kin)
      logical lokt(max_kin)
      complex*16 WA(6), WZ(6), WW(6)
      complex*16 WWeA(6), WWeZ(6),WWnZ(6),WWWA(6), WWWZ(6), WWnW(6)
      complex*16 EL(6), VE(6), MUP(6),MUM(6)
      complex*16 TEMP(6)

      integer hel,i,j,k, kin, mu,maxkin,bos!,sig
      real*8 mass2
      external mass2
      real*8 zero
      logical ldebug,lWA
      parameter (ldebug=.true.)
      parameter (zero = 0d0)

      complex*16 epsW(0:3,7),epsWV(0:3,-1:1,7) 
      complex*16 epsV(0:3,-1:1,3:4,-1:1,7) ! ( , mu-hel, u/d , quark-hel), coupl. const to quarks included
      complex*16 epsA(0:3,-1:1),epsZ(0:3,-1:1)
      real*8 momV(0:4,7),momW(0:4,7),momWV(0:4,7)
      common /qcdwznew/ epsW, epsWV, epsV,epsA,epsZ,momV,momW,momWV


      real*8 rn
      complex*16 phase

      double precision  clr, xm2, xmg, b
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6)
      maxkin = 1
      goto 12
      
      entry calcEWpartsNEW(v, lokt,bos,rn)
      maxkin = 7
      
 12   continue

      lWA=mod(bos,10).eq.1
!       sig = (37-bos)/5 ! +- 1
      do kin=1,maxkin 
      if(kin.eq.1.or.lokt(kin).or.ldebug)then
      
      if(bos.gt.10) then ! WpZ, WmZ
      phase = exp((0d0,2d0)*pi*rn)

c     Wave Functions for e+, ne, mu+, mu-
      do hel = -1,1,2
	if(bos/10.eq.3) then
	  CALL IXXXXX(v(0,1,kin),ZERO ,1,-1,EL(1))              
	  CALL OXXXXX(v(0,2,kin),ZERO ,-1,1,VE(1))

          if(.not.lWA) then
            CALL IXXXXX(v(0,3,kin),ZERO ,-hel,-1,MUP(1))        
            CALL OXXXXX(v(0,4,kin),ZERO , hel, 1,MUM(1))
          endif
	else
	  CALL OXXXXX(v(0,1,kin),ZERO ,-1,1,EL(1))              
	  CALL IXXXXX(v(0,2,kin),ZERO ,1,-1,VE(1))

          if(.not.lWA) then
            CALL OXXXXX(v(0,3,kin),ZERO , hel, 1,MUM(1))        
            CALL IXXXXX(v(0,4,kin),ZERO ,-hel,-1,MUP(1))
          endif
	endif
	if (bos/10.eq.3) then
c     calculate A, Z, W, W3 current
          CALL JIOXXX(EL ,ve ,GWF ,WMASS   ,WWIDTH  ,WW)

          if(lWA) then
            call VXXXXX(v(0,3,kin),ZERO,hel,1,WA)
          else
            CALL JIOXXX(mup,mum,GAL ,ZERO    ,AWIDTH  ,WA)
            CALL JIOXXX(mup,mum,GZL ,ZMASS   ,ZWIDTH  ,WZ)
c     e -> Z e
            CALL FVIXXX(EL,wz,GZL ,ZERO    ,ZERO    ,temp)
            CALL JIOXXX(temp,ve,GWF ,WMASS   ,WWIDTH  ,WWeZ)
          endif

c     e -> A e
	  CALL FVIXXX(EL,wa,GAL ,ZERO    ,ZERO    ,temp)
	  CALL JIOXXX(temp,ve,GWF ,WMASS   ,WWIDTH  ,WWeA)

          if(.not.lWA) then
c     ve -> Z ve
            CALL FVOXXX(ve,wz,GZN ,ZERO    ,ZERO    ,temp)
            CALL JIOXXX(EL,temp,GWF ,WMASS   ,WWIDTH  ,WWnZ)

c     vm -> W+ mu-
            CALL FVOXXX(mum,ww,GWF ,ZERO    ,ZERO    ,temp)
            CALL JIOXXX(mup,temp,GWF ,WMASS   ,WWIDTH  ,WWnW)

c     W -> W Z, W A
            CALL JVVXXX(ww,wz,GWWZ ,WMASS   ,WWIDTH  ,WWWZ)
          endif
	  CALL JVVXXX(ww,wa,GWWA ,WMASS   ,WWIDTH  ,WWWA)

	else
c     calculate A, Z, W, W3 current
          CALL JIOXXX(ve,el ,GWF ,WMASS   ,WWIDTH  ,WW)
          if(lWA) then
            call VXXXXX(v(0,3,kin),ZERO,hel,1,WA)
          else
            CALL JIOXXX(mup,mum,GAL ,ZERO    ,AWIDTH  ,WA)
            CALL JIOXXX(mup,mum,GZL ,ZMASS   ,ZWIDTH  ,WZ)
c     e -> Z e
            CALL FVOXXX(EL,wz,GZL ,ZERO    ,ZERO    ,temp)
            CALL JIOXXX(ve,temp,GWF ,WMASS   ,WWIDTH  ,WWeZ)
          endif
  ! 	print*,WA
  ! 	print*,WZ

c     e -> A e
	  CALL FVOXXX(EL,wa,GAL ,ZERO    ,ZERO    ,temp)
	  CALL JIOXXX(ve,temp,GWF ,WMASS   ,WWIDTH  ,WWeA)

          if(.not.lWA) then
c     ve -> Z ve
            CALL FVIXXX(ve,wz,GZN ,ZERO    ,ZERO    ,temp)
            CALL JIOXXX(temp,el,GWF ,WMASS   ,WWIDTH  ,WWnZ)

c     vm -> W+ mu-
            CALL FVIXXX(mup,ww,GWF ,ZERO    ,ZERO    ,temp)
            CALL JIOXXX(temp,mum,GWF ,WMASS   ,WWIDTH  ,WWnW)   
          endif

c     W -> W Z, W A
	  CALL JVVXXX(wz,ww,GWWZ ,WMASS   ,WWIDTH  ,WWWZ) ! ordering of wa, ww is important !
	  CALL JVVXXX(wa,ww,GWWA ,WMASS   ,WWIDTH  ,WWWA)

	endif

	do i=1,4
          if(lWA) then
            epsWV(i-1,hel,kin)=GWF(1)*(+WWWA(i)+WWeA(i))
          else
            epsWV(i-1,hel,kin)=GWF(1)*(+WWWZ(i)+WWWA(i)+WWeA(i)+WWeZ(i)+WWnZ(i)+WWNW(i))
          endif
	  
	  do k=-1,1,2 ! quark hel
	  do j=3,4 ! quark type
            if(lWA) then
              epsV(i-1,hel,j,k,kin) = - WA(i)*CLR(j,1,k) ! going from D_mu = d_mu + i g A  ->  D_mu = d_mu - igA, see 1209.6213 (*sign to account for conventions of CLR() and GWA ...*)
	    else
              epsV(i-1,hel,j,k,kin) = -WZ(i)*CLR(j,2,k) - WA(i)*CLR(j,1,k) ! going from D_mu = d_mu + i g A  ->  D_mu = d_mu - igA, see 1209.6213 (*sign to account for conventions of CLR() and GWA ...*)
            endif
	  enddo
	  enddo
	  epsW(i-1,kin) = WW(i)*GWF(1)  
	  if(kin.eq.1) then
	    epsA(i-1,hel) = -WA(i)
	    epsZ(i-1,hel) = -WZ(i)
            if(lWA) epsZ(i-1,hel) = 0d0
	  endif
	enddo
      enddo ! hel

      do mu=0,3 ! multiply with random phase and calculate momenta
	epsWV(mu,0,kin) = epsWV(mu,-1,kin)*phase +  epsWV(mu,1,kin)*CONJG(phase)
        do k=-1,1,2
        do j=3,4
	  epsV(mu,0,j,k,kin) = phase*epsV(mu,-1,j,k,kin) +epsV(mu,1,j,k,kin) *conjg(phase)
        enddo
        enddo
	if(kin.eq.1) then ! only needed for fermion loop contributions
	  epsA(mu,0) = phase*epsA(mu,-1)+conjg(phase)*epsA(mu,1)
	  epsZ(mu,0) = phase*epsZ(mu,-1)+conjg(phase)*epsZ(mu,1)
	endif
        if (ldoblha) then
          epsWV(mu,0,kin) = (2-blha_ranhelcomb)*epsWV(mu,+1,kin) + (blha_ranhelcomb-1)*epsWV(mu,-1,kin)
          do k=-1,1,2
          do j=3,4
            epsV(mu,0,j,k,kin) = (2-blha_ranhelcomb)*epsV(mu,+1,j,k,kin) + (blha_ranhelcomb-1)*epsV(mu,-1,j,k,kin)
          enddo
          enddo
          if(kin.eq.1) then ! only needed for fermion loop contributions
            epsA(mu,0) = (2-blha_ranhelcomb)*epsA(mu,+1)+(blha_ranhelcomb-1)*epsA(mu,-1)
            epsZ(mu,0) = (2-blha_ranhelcomb)*epsZ(mu,+1)+(blha_ranhelcomb-1)*epsZ(mu,-1)
          endif
        endif
	momV (mu,kin) = v(mu,3,kin)+v(mu,4,kin)
	momW (mu,kin) = v(mu,1,kin)+v(mu,2,kin)
	momWV(mu,kin) = momV(mu,kin)+momW(mu,kin)
      enddo
      
      
      else ! Wp, Wm
	if(bos.eq.3) then
	  CALL IXXXXX(v(0,1,kin),ZERO ,1,-1,EL(1))              
	  CALL OXXXXX(v(0,2,kin),ZERO ,-1,1,VE(1))
	  CALL JIOXXX(EL ,ve ,GWF ,WMASS   ,WWIDTH  ,WW)
	else
	  CALL OXXXXX(v(0,1,kin),ZERO ,-1,1,EL(1))              
	  CALL IXXXXX(v(0,2,kin),ZERO ,1,-1,VE(1))
	  CALL JIOXXX(ve,el ,GWF ,WMASS   ,WWIDTH  ,WW)
	endif
	do mu=0,3
	  epsW(mu,kin) = WW(mu+1)*GWF(1)  
	  momW (mu,kin) = v(mu,1,kin)+v(mu,2,kin)
	enddo
      endif
            
      momV(4,kin) = mass2(momV(0,kin))
      momW(4,kin) = mass2(momW(0,kin))
      momWV(4,kin) = mass2(momWV(0,kin))
      endif !lokt
      enddo ! kin
      end
