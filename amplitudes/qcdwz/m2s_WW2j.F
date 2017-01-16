c ************************************************************************************
c Author: Matthias Kerner
c Date: 2012/13
c ************************************************************************************

! #define NINH_COMPARE
! #define POWHEG_COMPARE
! #define CALC_QCD_EW_INTERF
      real*8 function m2s_QCDww2j(
     &                   xi,  !in:  Feynman x parameters of incoming fermions
     &                   p,   !in:  momenta associated with external fermions
     &                   v,     !in:  momenta associated with Z/H dec. fermions
     &                   rn2,
     &                   Nbos, ! 33 = WpWp, 44 = WmWm
     &                   nlo
     &                        )
      implicit none
c
c declare input/output variables
c
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/loops/Utilities/GaugeTest.inc"
#include "common.inc"
      integer nbos,nlo
      real*8 xi(nx), p(0:3,max_p,max_kin), v(0:3,max_v,max_kin), m2s,rn2(2)
c
c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"
      real*8  alphas, xmt, alfa, xmz, xmw,s2w, xmh
      common /bkopin/   alphas, xmt, alfa, xmz, xmw, s2w, xmh
c
c declare local variables
c
      real*8 q_sf
      integer init/0/, I, J, mu
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


c local variables
      integer FL4qcdVV, ires_min,ires_max
      save ires_min, ires_max
      real*8 mass2,mjj2
      external FL4qcdVV,mass2,mjj2
      logical ldebug,allsub
      parameter (ldebug = .false., allsub=.false.)

#ifdef NINH_COMPARE
      integer hel1,hel2
      real*8 debugsum
#endif
      real*8 xiFinColl(2,2)
      

      complex*16 lepton(6,3)

#ifdef CALC_QCD_EW_INTERF
c for EW contributions
      real*8 vtemp
      common /wwInterferences/ vtemp(0:3,4)
#endif

      complex*16 epsW
      real*8 momW
      common /QCDWWcurrents/ epsW(0:3,2,7),momW(0:4,2,7)
c EW/QCD interferences
      integer WWJJ_interf
      common /cWWJJ_interf/ WWJJ_interf

#ifdef POWHEG_COMPARE
         wmass = 80.419d0
         wwidth= 2.141d0
         GWF(1) = - sqrt(0.43907954221934414d0/2d0)
#endif
!       wwidth =  2.097672973437829d0
!       als= 0.129808134118326d0
!       als = 0.117981587488294d0
!       wwidth=2.085d0

#ifdef NINH_COMPARE
!       wwidth=2.097672973437829d0
#endif

      if ( init .lt. 1 ) then

         write(6,*) " "
        if(nbos .eq. 33) then
         write(6,*) "W+W+ jj amplitude square information:"
        else
         write(6,*) "W-W- jj amplitude square information:"
        endif
         write(6,*) "-----------------------------------------------"
         write(6,*) " "
         write(6,*) "alpha parameter of dipoles: ",alphaDip
         write(6,*) "precision of gauge test:   ",GaugeAccuracy
         write(6,*) "-----------------------------------------------"
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
         print*,"W momenta"
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

      q_sf = sqrt(mufsq(2,1))
      call pdfproton( xi(2), q_sf, pdf(-6,2,1) )        ! f_a(x2)
      do i=-5,5
        pdf(i,1,1) = pdf(i,1,1)/xi(1)
        pdf(i,2,1) = pdf(i,2,1)/xi(2)
      enddo
      
      if(nlo.eq.3) then
        xiFinColl(1,1) = xi(1)
        xiFinColl(2,1) = xi(2)
        xiFinColl(1,2) = xi(1)/((1d0-xi(1))*rn2(2)+xi(1)) ! transform rn el [0;1] to x el [z;1], get xa = xxa/x, multiply jacobian to pdf
        xiFinColl(2,2) = xi(2)/((1d0-xi(2))*rn2(2)+xi(2)) 

        q_sf = sqrt(mufsq(1,1))
        call pdfproton( xiFinColl(1,2), q_sf, pdf(-6,1,2) )        ! f_a(x1)

        q_sf = sqrt(mufsq(2,1))
        call pdfproton( xiFinColl(2,2), q_sf, pdf(-6,2,2) )        ! f_a(x2)
        
        do i=-5,5
          pdf(i,1,2) = pdf(i,1,2)/xiFinColl(1,2)
          pdf(i,2,2) = pdf(i,2,2)/xiFinColl(2,2)
        enddo
        
        call CScolAmpQNewReset
      endif


! ! ! ! ! ! ! ! ! ! ! ! c Paper
! ! ! ! ! !       p(0:3,1,1) = (/500.000000000000d0,500.000000000000d0,0.00000000000000d0,0.00000000000000d0/)
! ! ! ! ! !       p(0:3,2,1) =  (/500.000000000000d0,-500.000000000000d0,0.00000000000000d0,0.00000000000000d0/)
! ! ! ! ! !       p(0:3,3,1) =  (/54.2314070117999d0,-31.1330162081798d0,-7.92796656791140d0,43.6912823611163d0/)
! ! ! ! ! !       p(0:3,4,1) = (/214.488870161418d0,-27.0607980217775d0,-98.5198083786150d0,188.592247959949d0/)
! ! ! ! ! !       v(0:3,1,1) = (/85.5312248384887d0,-8.22193223977868d0,36.1637837682033d0,-77.0725048002414d0/)
! ! ! ! ! !       v(0:3,2,1) = (/181.428811610043d0,-57.8599829481937d0,-171.863734086635d0,-5.61185898481311d0/)
! ! ! ! ! !       v(0:3,3,1) = (/82.8493010774356d0,-65.9095476235891d0,-49.8952157196287d0,5.51413360058664d0/)
! ! ! ! ! !       v(0:3,4,1) = (/381.470385300815d0,190.185277041519d0,292.042940984587d0,-155.113300136598d0/)
!

!       als=1d0
!       alphas=1d0
!       G=sqrt(4d0*pi)
!       gg(1)=-G
!       gg(2)=-G
! !       GWF(1)=4d0*pi/137d0/WMASS**2*sqrt(2d0)/8d0
! !       print*,GWF(1)
!       GWF(1)=sqrt(4d0*pi)/sqrt(1d0-(WMASS/ZMASS)**2)/sqrt(2d0)
!       print*,GWF(1)
!       WWIDTH=2.09761d0
!       GWF(2)=0d0
!       p(0:3,1,1) = (/18.3459102072588d0, 0.0d0, 0.0d0, 18.3459102072588d0/)
!       p(0:3,2,1) = (/4853.43796816526d0, 0.0d0, 0.0d0, -4853.43796816526d0/)
!       p(0:3,3,1) = (/235.795970274883d0, -57.9468743482139d0, -7.096445419113396d-15, -228.564869022223d0/)
!       p(0:3,4,1) = (/141.477229270568d0, -45.5048903376581d0, -65.9221967646567d0, -116.616359620580d0/)
!       v(0:3,1,1) = (/276.004829895761d0, 31.4878768361538d0, -8.65306166938040d0, -274.066240646098d0/)
!       v(0:3,2,1) = (/1909.28515244344d0, 29.6334571080402d0, 40.1409467910328d0, -1908.63311192893d0/)
!       v(0:3,3,1) = (/2241.46026948104d0, 28.1723094714198d0, 30.2470561132914d0, -2241.07910976778d0/)
!       v(0:3,4,1) = (/67.7604270068059d0, 14.1581212702582d0, 4.18725552971283d0, -66.1323669723852d0/)
! 
!       print*,(sum(p(mu,3:4,1))+sum(v(mu,1:4,1))-sum(p(mu,1:2,1)),mu=0,3)


      do mu=0,3
        momW(mu,1,1) = v(mu,1,1)+v(mu,2,1)
        momW(mu,2,1) = v(mu,3,1)+v(mu,4,1)
      enddo

      momW(4,1,1) = mass2(momW(0,1,1))
      momW(4,2,1) = mass2(momW(0,2,1))

!       print*,sqrt(momW(4,1,1))
!       print*,sqrt(momW(4,2,1))

c calc W currents
      do j=1,2
	if(nbos.eq.33) then
	  CALL IXXXXX(v(0,2*j-1,1),0d0,1,-1,lepton(1,1))        
	  CALL OXXXXX(v(0,2*j  ,1),0d0,-1,1,lepton(1,2))
	  CALL JIOXXX(lepton(1,1) ,lepton(1,2) ,GWF ,WMASS   ,WWIDTH  ,lepton(1,3))
	else
	  CALL OXXXXX(v(0,2*j-1,1),0d0,-1,1,lepton(1,1))        
	  CALL IXXXXX(v(0,2*j  ,1),0d0,1,-1,lepton(1,2))
	  CALL JIOXXX(lepton(1,2) ,lepton(1,1) ,GWF ,WMASS   ,WWIDTH  ,lepton(1,3))
	endif
	do mu=0,3
	  epsW(mu,j,1)=lepton(mu+1,3)*GWF(1)
	  momW(mu,j,1)=v(mu,2*j-1,1)+v(mu,2*j,1)
	enddo
	momW(4,j,1) = mass2(momW(0,j,1))
      enddo
      
#ifdef CALC_QCD_EW_INTERF
      vtemp = v(:,:4,1)
#endif


#ifdef NINH_COMPARE
!         als(1,1) = 0.122058009253390d0
!         als(2,1) = 0.122058009253390d0
        als(1,1) = 1d0
        als(2,1) = 1d0

        alphas=als(1,1)
        G = sqrt(4d0*pi*alphas)
        gg(1) = -G
        gg(2) = -G
        WMASS = 80.398d0
        mursq(1,1) = WMASS**2
	gwf(1) = -sqrt( 1.16637d-5/sqrt(2d0))*80.398d0*2d0
!         GWF(1) = -1.166370000000000d-5
!         do i=-6,6
!           pdf(i,1,1) =1d0+0.01d0*i
!           pdf(i,2,1) =1d0+0.02d0*i
!         enddo
      pdf=1d0

        p(0:3,1,1) = (/ 500.000000000000d0, 0.000000000000000d0,0.000000000000000d0,500.000000000000d0/)
        p(0:3,2,1) = (/ 500.000000000000d0, 0.000000000000000d0,0.000000000000000d0,-500.000000000000d0/)
        momW(0:3,1,1) = (/ 216.232776975442d0, -110.114342023960d0,14.5729785186567d0,-167.198192012463d0/)
        momW(0:3,2,1) = (/ 337.997244776963d0,-59.5891662343044d0,-150.880208615759d0,285.416525429737d0/)
        p(0:3,3,1) = (/158.506206722148d0,-60.0863508547814d0,136.307230097103d0,54.1681367009778d0/)
        p(0:3,4,1) = (/287.263771525447d0,229.789859113046d0,0.000000000000000d0,-172.386470118252d0/)
        p(0:3,5,1) = momW(0:3,1,1)
        p(0:3,6,1) = momW(0:3,2,1)
        momW(4,1,1) = mass2(momW(0,1,1))
        momW(4,2,1) = mass2(momW(0,2,1))


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
! ! ! ! ! ! 	p(0:3,1,1) = (/500.000000000000d0,500.000000000000d0,0.00000000000000d0,0.00000000000000d0/)
! ! ! ! ! ! 	p(0:3,2,1) =  (/500.000000000000d0,-500.000000000000d0,0.00000000000000d0,0.00000000000000d0/)
! ! ! ! ! ! 	p(0:3,3,1) =  (/54.2314070117999d0,-31.1330162081798d0,-7.92796656791140d0,43.6912823611163d0/)
! ! ! ! ! ! 	p(0:3,4,1) = (/214.488870161418d0,-27.0607980217775d0,-98.5198083786150d0,188.592247959949d0/)
! ! ! ! ! ! 	v(0:3,1,1) = (/85.5312248384887d0,-8.22193223977868d0,36.1637837682033d0,-77.0725048002414d0/)
! ! ! ! ! ! 	v(0:3,2,1) = (/181.428811610043d0,-57.8599829481937d0,-171.863734086635d0,-5.61185898481311d0/)
! ! ! ! ! ! 	v(0:3,3,1) = (/82.8493010774356d0,-65.9095476235891d0,-49.8952157196287d0,5.51413360058664d0/)
! ! ! ! ! ! 	v(0:3,4,1) = (/381.470385300815d0,190.185277041519d0,292.042940984587d0,-155.113300136598d0/)
! ! ! ! ! ! 
! ! ! ! ! ! c 2->6
! ! ! ! ! ! !         p(0:3,1,1) = (/ 500.000000000000d0, 0.000000000000000d0,0.000000000000000d0,500.000000000000d0/)
! ! ! ! ! ! !         p(0:3,2,1) = (/ 500.000000000000d0, 0.000000000000000d0,0.000000000000000d0,-500.000000000000d0/)
! ! ! ! ! ! !         p(0:3,3,1) = (/ 88.1907917270800d0, -71.8338387342512d0,-26.8500535430029d0,-43.5498562901036d0/) 
! ! ! ! ! ! !         p(0:3,4,1) = (/ 457.182458695087d0,  404.255345805865d0, 0.000000000000000d0,-213.526148108059d0/)
! ! ! ! ! ! !         v(0:3,1,1) = (/ 100.788423998349d0, -45.8888076188476d0, 49.1583762549171d0, 75.0731496033797d0/)
! ! ! ! ! ! !         v(0:3,2,1) = (/    199.936105394529d0,-187.698928314700d0,-25.1678128082751d0,64.1103716125017d0/)
! ! ! ! ! ! !         v(0:3,3,1) = (/149.590646659259d0,-96.4973730388405d0,2.67112519191297d0,114.273722503239d0/)
! ! ! ! ! ! !         v(0:3,4,1) = (/4.31157352569544d0,-2.33639809922604d0,0.188364904448005d0,3.61876067904234d0/)
! ! ! ! ! ! 

      debugsum = 0d0
      do hel1=-1,1
      do hel2=-1,1
!         print*,'helicities: ',hel1,hel2
        call VXXXXX(momW(0,1,1),WMASS,hel1,1,lepton(1,1))
        epsW(0:3,1,1)=lepton(1:4,1)
        call VXXXXX(momW(0,2,1),WMASS,hel2,1,lepton(1,1))
        epsW(0:3,2,1)=lepton(1:4,1)

        epsW(0:3,1,1) = epsW(0:3,1,1)*GWF(1)/sqrt(2d0) ! Wff coupling and identical factor
        epsW(0:3,2,1) = epsW(0:3,2,1)*GWF(1)
        res=0d0
#endif
 ! ! !    comparison paper   
!       id(1) = 1
!       id(2) = 3
!       id(3) = 2
!       id(4) = 4
!       
!       fsign(1) = 1
!       fsign(2) = -1
!       fsign(3) = -1
!       fsign(4) = 1
!       
!       print*,'a'
!       call amp4qwwVirt(P,id,fsign,2,res)
!       print*,'b'
!       pause

c reset subprocess counter
      ires_min = FL4QCDVV(id,id,-sub_number) +1 

 
c begin 4q processes
c
c   q1 -------------- q3
c            9
c            6
c            9
c   q2 -------------- q4
c
c
#ifndef NINH_COMPARE
c*****  1.  ********
c**Diagram**********  q1 q2 ---> q3 q4 W Z   **********************
c**Diagram**********  q2 q1 ---> q3 q4 W Z   **********************
c**Phasespace*******  1   2      3   4

      id(1) = 1
      id(2) = 2
      id(3) = 3
      id(4) = 4

      fsign(1) =  1
      fsign(2) =  1
      fsign(3) =  1
      fsign(4) =  1 

      if(allsub.or.sub_number.eq.1) call getRes4QWW(res,p,xifincoll,fsign,id,pdf,nbos,nlo )

c*****  2.  ********
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

      if(allsub.or.sub_number.eq.3) call getRes4QWW(res,p,xifincoll,fsign,id,pdf,nbos,nlo )
#endif
c*****  3.  ********
c**Diagram**********  q1 qb3 ---> qb2 q4 W Z   **********************
c**Phasespace*******  1   2        3   4       **********************

      id(1) = 1
      id(2) = 3
      id(3) = 4
      id(4) = 2

      fsign(1) =  1
      fsign(2) = -1
      fsign(3) =  1
      fsign(4) = -1
      

      if(allsub.or.sub_number.eq.2) call getRes4QWW(res,p,xifincoll,fsign,id,pdf,nbos,nlo ) 

! flip initial state

      id(1) = 3
      id(2) = 2
      id(3) = 1
      id(4) = 4

      fsign(1) = -1
      fsign(2) =  1
      fsign(3) = -1
      fsign(4) =  1
      
      
      if(allsub.or.sub_number.eq.2) call getRes4QWW(res,p,xifincoll,fsign,id,pdf,nbos,nlo ) 

c
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
      do while((i.le.ires_max).and.(iprocess.eq.0))
         weight = weight+res(i)/m2s
         if(weight.ge.rn2(1)) then
            iprocess = i
         endif
         i = i+1
      enddo

      
*     scale for running alphas
#ifndef CALC_QCD_EW_INTERF
      m2s = m2s*(als(1,1)/alphas)**2
#else
      if(nlo.ne.0) then
        print*,__FILE__,__LINE__
        stop
      endif
#endif

      m2s_QCDww2j =m2s

#ifdef NINH_COMPARE
!       print*,m2s
!       pause
      debugsum=debugsum+m2s
      enddo
      enddo
!       m2s_QCDww2j = debugsum
! ! ! ! !       m2s = debugsum
      print*,'sum', debugsum
#endif
!       pause

      RETURN
      END


      subroutine getRes4QWW(res,p,xifincoll,fsign,id,pdf,Nbos,nlo)
      implicit none
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/global.inc"
      real*8 res(1:maxNumSubProc) !out: subprocess cross-sections
      real*8 p(0:3,max_p,max_kin)
      integer fsign(4), id(4),nlo
      real*8 pdf(-6:6,2,2)
      integer nbos

      integer FL4qcdVV
      external FL4qcdVV

      real*8 POLCOL, polcol0    ! ave over intial polarizations and colors
      real*8 xiFinColl(2,2)
      integer colstr,j,i,k
      integer iflav(1:4)
      integer temp
      logical ldebug,idFl
      parameter (ldebug = .false.)
      real*8 amp(0:2,0:2),QCDVVjjFinCollBornKin

      external QCDVVjjFinCollBornKin
c EW/QCD interferences
      integer WWJJ_interf
      common /cWWJJ_interf/ WWJJ_interf
      
      polcol0 = 1.0d0/(4.0d0 * 9.0d0)

      if(nlo.eq.0) then
	call amp4qwwBorn(P,id,fsign,amp)
      else
        if(nlo.le.2) then
#ifdef WITH_NLO
          call amp4qwwVirt(P,id,fsign,nlo,amp)
#endif
        endif
      endif
      
c sum over all quark combinations
      do j = 0,2,2  ! quark family 1 and 2
!       print*,"same family"
      do k = 0,2,2
        IdFl = j.eq.k
c fill flavor information for pdf and lha output
	if(Nbos.eq.33) then  ! w+ case
	    iflav(id(1)) = (2+j)*fsign(1)
	    iflav(id(3)) = (1+j)*fsign(3)
	    iflav(id(2)) = (2+k)*fsign(2)
	    iflav(id(4)) = (1+k)*fsign(4)
	elseif(Nbos.eq.44) then ! w- case
	    iflav(id(1)) = (1+j)*fsign(1)
	    iflav(id(3)) = (2+j)*fsign(3)
	    iflav(id(2)) = (1+k)*fsign(2)
	    iflav(id(4)) = (2+k)*fsign(4)
	else
	    stop
	endif

c  pick polcol
	if(iflav(3).eq.iflav(4)) then
	    polcol = polcol0/2.0d0
	else
	    polcol = polcol0
	endif
#ifdef CALC_QCD_EW_INTERF 
        ! only t-channel EW * u/s channel QCD -> compensate missing factor 2 here, 
        if(WWJJ_interf.eq.1) polcol=polcol0
#endif

c get m2s*pdf
	if(idFl) then  ! same q families
          do colstr=1,2
            i = FL4qcdVV(iflav,id,colstr)
            if(nlo.eq.3) then
              if (colstr.eq.1) then
                res(i) = QCDVVjjFinCollBornKin(p,xiFinColl,id,fsign,iflav,0,nbos,pdf) ! born is included here
              else
                res(i) = 0d0
              endif
            else
              res(i) = POLCOL * amp(colstr,0)
     $		      *pdf(sign1*iflav(1),1,1)*pdf(sign2*iflav(2),2,1)
            endif
            if(ldebug) print*,i, iflav!,res(i)
          enddo
	else
	  i = FL4qcdVV(iflav,id,1)
          if(nlo.eq.3) then
            res(i) = QCDVVjjFinCollBornKin(p,xiFinColl,id,fsign,iflav,1,nbos,pdf) ! born is included here
          else
            res(i) = POLCOL * amp(0,1)
     $		      *pdf(sign1*iflav(1),1,1)*pdf(sign2*iflav(2),2,1)
          endif
	  if(ldebug) print*,i, iflav!,res(i)
          if(WWJJ_interf.ne.1) then ! not for VBF contribution
	  if(iflav(1)*iflav(2).lt.0) then ! include colstr 2 for q qb initial states
	    ! swap particle 1 and 2; change sign if needed
	    temp         = sign(iflav(id(1)),fsign(2))
	    iflav(id(1)) = sign(iflav(id(2)),fsign(1))
	    iflav(id(2)) = temp
	    i = FL4qcdVV(iflav,id,2)
            if(nlo.eq.3) then
              res(i) = QCDVVjjFinCollBornKin(p,xiFinColl,id,fsign,iflav,2,nbos,pdf) ! born is included here
            else
              res(i) = POLCOL * amp(0,2)
     $		      *pdf(sign1*iflav(1),1,1)*pdf(sign2*iflav(2),2,1)
            endif
	    if(ldebug) print*,i, iflav 
	  endif
          endif
	endif
      enddo
      enddo
      return
      end

       
c 4q WW Amplitude
c    q1 >>> #### >>> q3 
c	    ####          , WW
c	    ####
c    q2 >>> #### >>> q4  
c   res(colstr, flav config=0,2)
c      colstr=0..2: for catani seymour: <|>,  <|T1T2|>,  <|T1T3|>
c		    for lo: col summed, and the 2 different color flows
c      flavour: 0: same family (t and u channel)
c      flavour: 1: diff family (only t channel)
c      flavour: 2: diff family (only u channel)
      subroutine amp4qwwCS(P,kin,id,fsign,res)
      implicit none
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
      real*8 P(0:3,max_p,max_kin),res(0:2,0:2)
      integer kin, id(4),fsign(4)
      logical lnlo
      
      real*8 pbar(0:3,4),pdia(0:3,4),vtemp
      common /wwInterferences/ vtemp(0:3,4)
      real*8 mompsiW(0:4,4,2)
      complex*16 psi(2,-1:1,4), psiW(2,4,2)
      complex*16 amp(2)
      complex*16 gluWij(0:3,1:3,2:4,2)
      complex*16 tempCurr(0:3,-1:1,4)
      complex*16 dotcc
      real*8 qjj2
      real*8 squared(2,2)
      
      integer i,j,mu,tempi4(4), lepkin
      integer ikin, in,out

#include "VBFNLO/utilities/scales.inc"
      real*8  alphas, xmt, alfa, xmz, xmw,s2w, xmh
      common /bkopin/   alphas, xmt, alfa, xmz, xmw, s2w, xmh

      complex*16 epsW
      real*8 momW
      common /QCDWWcurrents/ epsW(0:3,2,7),momW(0:4,2,7)


c EW/QCD interferences
      integer WWJJ_interf
      common /cWWJJ_interf/ WWJJ_interf


      lnlo = .true.
      ikin = kin
      goto 11
      
      entry amp4qwwBorn(P,id,fsign,res)
      ikin = 1
      lnlo= .false.
      
 11   continue
      lepkin = ikin
      if(ikin .gt.7) lepkin = 1

c new numbering convention
c in   q1 >>> #### >>> q2  out
c             ####          , WZ
c	      ####
c in2  q3 >>> #### >>> q4  out2

      do mu=0,3
	pbar(mu,1)=p(mu,id(1),ikin)
	pbar(mu,2)=p(mu,id(3),ikin)
	pbar(mu,3)=p(mu,id(2),ikin)
	pbar(mu,4)=p(mu,id(4),ikin)
	pdia(mu,1)=pbar(mu,1)*fsign(1)
	pdia(mu,2)=pbar(mu,2)*fsign(3)
	pdia(mu,3)=pbar(mu,3)*fsign(2)
	pdia(mu,4)=pbar(mu,4)*fsign(4)
      enddo      
      

c get bra and ket 
      tempi4(1) = fsign(1)
      tempi4(2) = fsign(3)
      tempi4(3) = fsign(2)
      tempi4(4) = fsign(4)
      call psi0m(4,pbar,tempi4,psi)


!       print*,psi(1:2,-1,1)
c new convention
c in   >W>>>>W>  out
c          S         , WZ
c	   S
c in2  >>>>>>>>  out2

c W emission
      do in=1,3,2
      do j=1,2
	call ket2c(psi(1,-1,in),.true.,pdia(0,in),-1,momW(0,j,lepkin),epsW(0,j,lepkin),psiW(1,in,j),mompsiW(0,in,j))
      enddo
      enddo

      do out=2,4,2
      do j=1,2
	call bra2c(psi(1,-1,out),.true.,pdia(0,out),-1,momW(0,j,lepkin),epsW(0,j,lepkin),psiW(1,out,j),mompsiW(0,out,j))
      enddo
      enddo

c calculate glu currents
      do in=1,3,2
      do out=2,4,2
      do i=1,2
	call curr(-1,psi(1,-1,out),psiW(1,in,i),tempCurr(0,-1,1))
	call curr(-1,psiW(1,out,i),psi(1,-1,in),tempCurr(0,-1,2))
	do mu=0,3
	  gluWij(mu,in,out,i) = tempCurr(mu,-1,1)+tempCurr(mu,-1,2)
	enddo
      enddo
      enddo
      enddo
c calc amplitudes      
      amp(1) = ( dotcc(gluWij(0,1,2,1),gluWij(0,3,4,2))/qjj2(pdia(0,1),momPsiW(0,2,1))
     1	        +dotcc(gluWij(0,1,2,2),gluWij(0,3,4,1))/qjj2(pdia(0,1),momPsiW(0,2,2)))*G**2
      amp(2) =-( dotcc(gluWij(0,1,4,1),gluWij(0,3,2,2))/qjj2(pdia(0,1),momPsiW(0,4,1))
     1    	+dotcc(gluWij(0,1,4,2),gluWij(0,3,2,1))/qjj2(pdia(0,1),momPsiW(0,4,2)))*G**2

c square amplitudes
      do j=0,2
	res(j,0) = 0d0
	res(j,1) = 0d0
	res(j,2) = 0d0
      enddo


      squared(1,1) = dble(amp(1)*dconjg(amp(1)))
      squared(2,2) = dble(amp(2)*dconjg(amp(2)))
      squared(1,2) = dble(amp(1)*dconjg(amp(2)))

      res(0,1) = 2d0*squared(1,1)
      res(0,2) = 2d0*squared(2,2)
      res(0,0) = 2d0*(squared(1,1)+squared(2,2)-2d0/3d0*squared(1,2))

      if(lnlo) then ! CS color factors
	res(1,0) = -2d0/3d0*(squared(1,1)+squared(2,2)-10d0/3d0*squared(1,2))
	res(2,0) = +1d0/3d0*(squared(1,1)-7d0*squared(2,2)-2d0/3d0*squared(1,2))
	res(1,1) = -2d0/3d0*squared(1,1)
	res(2,1) = +1d0/3d0*squared(1,1)
	res(1,2) = -2d0/3d0*squared(2,2)
	res(2,2) = -7d0/3d0*squared(2,2)
      else ! get colorstr weighted amp for lha output
	res(1,0) = res(0,0)*squared(1,1)/(squared(1,1)+squared(2,2))
	res(2,0) = res(0,0)-res(1,0)
      endif


      end       


           subroutine fillColoredPartonsWZ(if1,if2,if3,if4,id,colstr,listp1)
           implicit none
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/lha.inc"

c          COMMOM BLOCK FOR DIAGRAMMATIC LABELS
           integer id(4)
           integer i,if1,if2,if3,if4,colstr,listp1,listp
           integer ngluon,npart
           parameter (npart = 4)
           integer iflav(npart),pi(1:2)
           integer col(1:2)
c
           listp = listp1 + numdecay
           lnup(listp) = numParticles+2
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
          lidup(1,listp)=if1
          lidup(2,listp)=if2
          lidup(3,listp)=if3
          lidup(4,listp)=if4

c         lSPINUP goes here:
c          do i = 1,4                   ! quarks and gluons are unpolarized
c             lspinup(i,listp) = 9      ! handle them in fillhepup.f
c          enddo
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
          iflav(1) = if1
          iflav(2) = if2
          iflav(3) = if3
          iflav(4) = if4
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
                licolup(1,id(1),listp) = 501*col(1) +
     $               501*col(2) ! color
                licolup(2,id(1),listp) = 0 ! anti-color
             else
                licolup(1,id(1),listp) = 0 ! color
                licolup(2,id(1),listp) =  501*col(1) +
     $               501*col(2) ! anti-color
             endif
c
             if(id(2).lt.3) then
                licolup(1,id(2),listp) = 502*col(1) +! color
     $               502*col(2)
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
                pi(1) = 2
                pi(2) = 4
             case(2)
                pi(1) = 4
                pi(2) = 2
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



