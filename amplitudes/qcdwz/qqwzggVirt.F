c ************************************************************************************
c Author: Matthias Kerner
c Date: 2012/13
c ************************************************************************************

! ! #define CALC_WJJ      
      SUBROUTINE ug2dgwzVirt(P,id,fsign,nbos,nlo,lephel,ghel2,ghel4,rn,res) ! sum over glu hel if ghel ==0

      use globalvars, only: ldoblha

      implicit none
#include "common.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"
      real*8 P(0:3,max_p,max_kin), res(0:2),rn ! for gauge test
      integer id(4),fsign(4),nbos,nlo,lephel,ghel2 ,ghel4 

      
      complex*16 epsW(0:3,7),epsWV(0:3,-1:1,7) 
      complex*16 epsV(0:3,-1:1,3:4,-1:1,7) ! ( , mu-hel, u/d , quark-hel), coupl. const to quarks included
      complex*16 epsA(0:3,-1:1),epsZ(0:3,-1:1)
      real*8 momV(0:4,7),momW(0:4,7),momWV(0:4,7)
      common /qcdwznew/ epsW, epsWV, epsV,epsA,epsZ,momV,momW,momWV
 

      real*8 pbar(0:3,2),qbar(0:4,2),pdia(0:3,2),qdia(0:4,3),kdia(0:3,3),momTemp(0:4)
      complex*16 qdiac(0:4,3),kdiac(0:3,3)
      real*8 momPsiW(0:4,2),momPsiV(0:4,2),momPsiWV(0:4,2),momPsiG(0:4,2,2)
      complex*16 Psi(2,-1:1,2), PsiW(2,2),PsiV(2,2),PsiWV(2,2), PsiG(2,2,2), PsiWG(2,2,2)
      complex*16 epsG(0:3,3), PsiTemp(2,-1:1),PsiVV(2,2), amp(3),glutemp(0:3)
      complex*16 gluW(0:3), gluWV(0:3), gluWG(0:3,2,2), tempcurr(0:3,-1:1,4)
      complex*16 ampb(3),ampv(3),ampi(3), tempborn, tempvirt,tempvirt2,tempvirt3,tempvirt4
      complex*16 ortb(3),ortv(3),orti(3),ort(3)
      complex*16 s1c, dotcc, dotrc,tempc
      real*8 mass2,musq
      integer i,j,mu,tempi2(2),qtype(2)
      integer glu1,glu2,div
      integer cpl(3),nfi(3),nGaugeAbe,nGaugeNonAbe
      real*8 mf(3)
      logical lsimple,lborn,lbox,lpent,lhex,lCSI,lnloWV,lnloW,lV,lVV,lG,lGG,ldebug,tempbool
      complex*16 hexamp(6,3),hexres(13,30),hexgauge(4,30),hexborn(30),hexgaugeb(4,30)
      complex*16 boxgauge(2,3),boxres(5,3),boxgaugeb(2,3),boxborn(3),nothing(9),checkg(5,3,2),checkgb(2,2)
      complex*16 pengauge(3,8),penres(8,8),pengaugeb(3,8),penborn(8),floopres(6)
      real*8 gaugeReal(2,4),floopgauge(6),logmuosij(3,4)
      logical wardidtest(0:22)
      real*8 td,gaugeprecision,mjj2,dotrr
      external mjj2,dotrr
      real*8 CF,CA,CFCA,TR

      parameter( CF = 4d0/3d0, CA = 3d0, CFCA=-1d0/6d0, TR=0.5d0)

      logical ltop,printall
      logical printGaugeFail,IsNaN
      external isnan
!       parameter ()

      double precision  clr, xm2, xmg, b
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6)


      logical printGaugePrec, GaugeOK(2),return0ifGaugeFails
      logical recalculatePSpoint
      common /comRecalcPSpoint/ recalculatePSpoint
!       parameter ()
      real*8 comGaugePrec
      common /VVJJGaugePrec/ comGaugePrec
cfc
      common/ hexgaugeValue/gaugehexABEabsvalue
      real*8 gaugehexABEabsvalue
      logical wardidtest_Q(0:22),wardpAbe(2),wardpNoAbe(2),wardhAbe(2)
      logical wardhNoAbe(2)
      logical wardfb(2)
      common/wardpfermionbox/wardfb
      common/wardpenAbe/wardpAbe
      common/wardpenNoAbe/wardpNoAbe
      common/wardhexAbe/wardhAbe
      common/wardHexNoAbe/wardhNoAbe
      ldebug = .false.
! !       printall=.true.
      printGaugePrec = .false.
      printGaugeFail = .false.
      return0ifGaugeFails = .false.
! !  333  continue 
      ltop = NFl.ge.5

!       if(nloContrib.eq.0) then
	lborn  = nlo.eq.1
	
	lsimple= nlo.eq.2 ! fermion loops, glu triangle, self energies, counterterms
	lbox   = nlo.eq.2
	lpent  = nlo.eq.2
	lhex   = nlo.eq.2
	lCSI   = nlo.eq.2
	
	lnloW  = nlo.eq.2
	lnloWV = nlo.eq.2
	if(nlo.eq.3) then
          res(0) = 0d0; res(1) = 0d0; res(2) = 0d0
	  return
	endif
      if (ldoblha) then
        mursq(1,1) = blha_scale**2
        lborn=  .true. 
        lsimple=.true. 
        lbox =  .true. 
        lpent=  .true. 
        lhex =  .true. 
        lCSI =  .false. 
	lnloW = .true.
	lnloWV= .true.
      endif
        
!       else
! ! 	lborn  = nloContrib.eq.1
! ! 	lsimple= nloContrib.eq.2
! ! 	lbox   = nloContrib.eq.3
! ! 	lpent  = nloContrib.eq.4.or.nloContrib.eq.5
! ! 	lhex   = nloContrib.eq.6
! ! 	lCSI   = nloContrib.eq.7
! 
!         lborn  = nloContrib.eq.1
!         lsimple= nloContrib.eq.2
!         lbox   = nloContrib.eq.2
!         lpent  = nloContrib.eq.2
!         lhex   = nloContrib.eq.2
!         lCSI   = nloContrib.eq.2
!         lnloW=.true.
!         lnloWV=.true.
!       endif
! 	
! 	lnloW   = nloContrib.ne.5
! 	lnloWV   = nloContrib.ne.4
!       endif
	
      
!       lnloWV = .true. ! only inlcude virtuals with both w and Z coupling to quarks
!       lnloW  = .true. ! only inlcude virtuals with both only W coupling to quarks, W->4l
! 
! 
!       lborn  = .false.
!       lsimple= .true. ! self energies triangles, floops
!       lbox   = .true.
!       lpent  = .true.
!       lhex   = .true.
!       lCSI   = .false.
!       lnloW = .true.
!       lnloWV= .true.

!       lborn  = .false.
!       lsimple= .true. ! self energies triangles, floops
!       lbox   = .false.
!       lpent  = .false.
!       lhex   = .false.
!       lCSI   = .false.
!       lnloW = .true.
!       lnloWV= .true.
      
      div =0
      
      ldebug=ldebug.or.div.ne.0

      musq = mursq(1,1)
!       musq = 100000d0
      if(musq.le.0d0) then
	print*,"unreasonable musq"
      endif
      
      lV  = .true.
      lVV = .true.
      lG  = .true.
      lGG = .true.

      do i=1,3
	res(i-1) = 0d0
	amp(i) = 0d0
	ampb(i) = 0d0
	ampv(i) = 0d0
	ampi(i) = 0d0
      enddo
      
      if(ldebug) print*, id

      do mu=0,3
	pbar(mu,1)=p(mu,id(1),1)
	pbar(mu,2)=p(mu,id(3),1)
	qbar(mu,1)=p(mu,id(2),1)
	qbar(mu,2)=p(mu,id(4),1)
	pdia(mu,1)=pbar(mu,1)*fsign(1)
	pdia(mu,2)=pbar(mu,2)*fsign(3)
	qdia(mu,1)=qbar(mu,1)*fsign(2)
	qdia(mu,2)=qbar(mu,2)*fsign(4)
	qdia(mu,3)=qdia(mu,1)+qdia(mu,2)
      enddo      
      qbar(4,1) = 0d0
      qbar(4,2) = 0d0
      qdia(4,1) = 0d0
      qdia(4,2) = 0d0
      qdia(4,3) = mass2(qdia(0,3))
      
      qtype(1) = nbos/10 ! (W+Z)? u type : d type
      qtype(2) = 7-qtype(1)


#ifdef CALC_WJJ
      lV=.true.
      lVV=.false.
      lnloWV=.false.
      epsWV(:,lephel,:)=epsW
      momWV=momW
      momW=momW/2d0
      momV=momW
      qtype(1) = 3
      qtype(2) = 4
#endif

c get bra and ket 
      tempi2(1) = fsign(1)
      tempi2(2) = fsign(3)
      call psi0m(2,pbar,tempi2,Psi)
      
c get <q,glu| and |glu,q>

      call helvec(qbar(0,1),fsign(2),ghel2*fsign(2),epsG(0,1))
      call helvec(qbar(0,2),fsign(4),ghel4*fsign(4),epsG(0,2))

!       pause
!       epsG(0:3,1) = qbar(0:3,1)
!       epsG(0:3,2) = qbar(0:3,2)
!       epsWV(0:3,lephel,1) = momWV(0:3,1)
      
      do mu = 0,3
c       momenta are defined as outgoing
	  epsG(mu,3) = (2d0*dotrc(qdia(0,2),epsG(0,1))*epsG(mu,2)
     &		      -2d0*dotrc(qdia(0,1),epsG(0,2))*epsG(mu,1)
     &		      +dotcc(epsG(0,1),epsG(0,2))*(qdia(mu,1)-qdia(mu,2)))/qdia(4,3) ! * [ta,tb]

      enddo
      
      do mu=0,3
        if(.not. lG) epsg(mu,3) = 0d0
        if(.not. lGG) epsg(mu,1) = 0d0
        if(.not. lGG) epsg(mu,2) = 0d0
        if(.not. lV) epsWV(mu,lephel,1) = 0d0
        if(.not. lVV) epsW(mu,1) = 0d0
        if(.not. lVV) epsV(mu,lephel,3,-1,1) = 0d0
        if(.not. lVV) epsV(mu,lephel,4,-1,1) = 0d0
        if(.not. lVV) epsA(mu,-1) = 0d0
        if(.not. lVV) epsZ(mu,-1) = 0d0
      enddo


c ---------------------------------------
c $$$$$$$$$$ B O R N $$$$$$$$$$$$$$$$$$$$
c ---------------------------------------

c get <q,V| and |V,q>
      call ket2c(psi(1,-1,1),.true.,pdia(0,1),-1,momW(0,1), epsW(0,1),PsiW(1,1),momPsiW(0,1))
      call bra2c(psi(1,-1,2),.true.,pdia(0,2),-1,momW(0,1), epsW(0,1),PsiW(1,2),momPsiW(0,2))
      call ket2c(psi(1,-1,1),.true.,pdia(0,1),-1,momV(0,1), epsV(0,lephel,qtype(1),-1,1),PsiV(1,1),momPsiV(0,1)) ! u for W+, d for W-
      call bra2c(psi(1,-1,2),.true.,pdia(0,2),-1,momV(0,1), epsV(0,lephel,qtype(2),-1,1),PsiV(1,2),momPsiV(0,2)) ! d for W+, u for W-     
      call ket2c(psi(1,-1,1),.true.,pdia(0,1),-1,momWV(0,1), epsWV(0,lephel,1),PsiWV(1,1),momPsiWV(0,1))
      call bra2c(psi(1,-1,2),.true.,pdia(0,2),-1,momWV(0,1), epsWV(0,lephel,1),PsiWV(1,2),momPsiWV(0,2))
      do i=1,2
	call ket2c(psi(1,-1,1),.true.,pdia(0,1),-1,qdia(0,i), epsG(0,i),PsiG(1,1,i),momPsiG(0,1,i))
	call bra2c(psi(1,-1,2),.true.,pdia(0,2),-1,qdia(0,i), epsG(0,i),PsiG(1,2,i),momPsiG(0,2,i))
      enddo

c get <q,VV| and |VV,q>
      call ket2c(psiV(1,1),.false.,momPsiV(0,1),-1,momW(0,1), epsW(0,1),PsiTemp(1,-1),momTemp)
      call ket2c(psiW(1,1),.false.,momPsiW(0,1),-1,momV(0,1), epsV(0,lephel,qtype(2),-1,1),PsiTemp(1,1),momTemp)
      PsiVV(1,1) = Psitemp(1,-1) + Psitemp(1,1) - PsiWV(1,1)
      PsiVV(2,1) = Psitemp(2,-1) + Psitemp(2,1) - PsiWV(2,1)
      call bra2c(psiV(1,2),.false.,momPsiV(0,2),-1,momW(0,1), epsW(0,1),PsiTemp(1,-1),momTemp)
      call bra2c(psiW(1,2),.false.,momPsiW(0,2),-1,momV(0,1), epsV(0,lephel,qtype(1),-1,1),PsiTemp(1,1),momTemp)
      PsiVV(1,2) = Psitemp(1,-1) + Psitemp(1,1) - PsiWV(1,2)
      PsiVV(2,2) = Psitemp(2,-1) + Psitemp(2,1) - PsiWV(2,2)
      
      do i=1,2
	call ket2c(psiW(1,1),.false.,momPsiW(0,1),-1,qdia(0,i),epsG(0,i),Psitemp(1,-1),momtemp)
	call ket2c(psiG(1,1,i),.false.,momPsiG(0,1,i),-1,momW(0,1),epsW(0,1),PsiTemp(1,1),momtemp)
	PsiWG(1,1,i) = Psitemp(1,-1) + PsiTemp(1,1)
	PsiWG(2,1,i) = Psitemp(2,-1) + PsiTemp(2,1)
	call bra2c(psiW(1,2),.false.,momPsiW(0,2),-1,qdia(0,i),epsG(0,i),Psitemp(1,-1),momtemp)
	call bra2c(psiG(1,2,i),.false.,momPsiG(0,2,i),-1,momW(0,1),epsW(0,1),PsiTemp(1,1),momtemp)
	PsiWG(1,2,i) = Psitemp(1,-1) + PsiTemp(1,1)
	PsiWG(2,2,i) = Psitemp(2,-1) + PsiTemp(2,1)
      enddo
      
c build currents
      call curr(-1,Psi(1,-1,2),PsiVV(1,1),tempcurr(0,-1,1))
      call curr(-1,PsiVV(1,2),Psi(1,-1,1),tempcurr(0,-1,2))
      call curr(-1,PsiW(1,2),PsiV(1,1),tempcurr(0,-1,3))
      call curr(-1,PsiV(1,2),PsiW(1,1),tempcurr(0,-1,4))
      do mu=0,3
        gluWV(mu) = +tempcurr(mu,-1,1)+tempcurr(mu,-1,2)+tempcurr(mu,-1,3)+tempcurr(mu,-1,4)
      enddo
      
      call curr(-1,Psi(1,-1,2),PsiW(1,1),tempcurr(0,-1,1))
      call curr(-1,PsiW(1,2),Psi(1,-1,1),tempcurr(0,-1,2))
      do mu=0,3
        gluW(mu) = -tempcurr(mu,-1,1)-tempcurr(mu,-1,2)
      enddo

c g(g1,g2) amplitudesl4*fsign(4),epsG(0,2))
! !       print*,qdia(4,:)
! !       qdia(0:3,1:2) = -qdia(0:3,1:2)
! !       call Vert3G(qdia(0,1),qdia(0,2),qdia(0,3),epsG(0,1),epsG(0,2),gluWV,musq,amp(1),-1)
! !       qdia(0:3,1:2) = -qdia(0:3,1:2)
! !       amp(1) = -amp(1)/qdia(4,3)
! !       amp(2) = -amp(1)
! !       print*,amp

      amp(1) = -dotcc(epsG(0,3),gluWV) ! col1 = in>g1>g2>out = tbta , but epsG3 has col [ta,tb]
      amp(2) = -amp(1)

c reevaluate current with some contributions disabled to split up nlo code
      if(lsimple) then
	call curr(-1,Psi(1,-1,2),PsiWV(1,1),tempcurr(0,-1,1)) 
	call curr(-1,PsiWV(1,2),Psi(1,-1,1),tempcurr(0,-1,2))
        do mu=0,3
          if(.not.lnloW) then
            gluWV(mu) = gluWV(mu)+tempcurr(mu,-1,1)+tempcurr(mu,-1,2) ! remove the W->4l contribution ( - (-curr))
          else if(.not.lnloWV) then
            gluWV(mu) = -tempcurr(mu,-1,1)-tempcurr(mu,-1,2) ! only use this contribution
          endif
        enddo
      endif

 
      do glu1=1,2
	glu2 = 3-glu1

c gluWG(:,glu,pos) currents
	call curr(-1,Psi(1,-1,2),PsiWG(1,1,glu1),tempcurr(0,-1,1))
	call curr(-1,PsiW(1,2),PsiG(1,1,glu1),tempcurr(0,-1,2))
        do mu=0,3
          gluWG(mu,glu1,1) = tempcurr(mu,-1,1)+tempcurr(mu,-1,2)
        enddo
	call curr(-1,PsiWG(1,2,glu1),psi(1,-1,1),tempcurr(0,-1,1))
	call curr(-1,PsiG(1,2,glu1),PsiW(1,1),tempcurr(0,-1,2))
        do mu=0,3
          gluWG(mu,glu1,2) = tempcurr(mu,-1,1)+tempcurr(mu,-1,2)
        enddo
      
	tempbool=.true. ! this is set to true, so that you can aplly gauge test for the born, should be changed to false in the end
c V-> WV amplitudes
	amp(glu1) = amp(glu1) + s1c(PsiG(1,2,glu2),epsWV(0,lephel,1),.true.,-1,PsiG(1,1,glu1))
c VVgg, ggVV amplitudes
	amp(glu1) = amp(glu1) 
     &	       - s1c(psiG(1,2,glu2),epsG(0,glu1),tempbool,-1,PsiVV(1,1))
     &	       - s1c(psiVV(1,2),epsG(0,glu2),tempbool,-1,PsiG(1,1,glu1))
c all other amplitudes
	amp(glu1) = amp(glu1)
     &		    - s1c(psiV(1,2),epsG(0,glu2),tempbool,-1,PsiWG(1,1,glu1))
     &		    - s1c(psiG(1,2,glu2),epsV(0,lephel,qtype(2),-1,1),.true.,-1,PsiWG(1,1,glu1))

	amp(glu1) = amp(glu1)
     &		    - s1c(psiWG(1,2,glu2),epsG(0,glu1),tempbool,-1,PsiV(1,1))
     &		    - s1c(psiWG(1,2,glu2),epsV(0,lephel,qtype(1),-1,1),.true.,-1,PsiG(1,1,glu1))
      enddo
	
c multiply coupl. const. and get |amp|^2
      amp(1) = amp(1)*GG(1)**2
      amp(2) = amp(2)*GG(1)**2
! 	print*,amp(1,hel1,hel2),amp(2,hel1,hel2)
      if (lborn) then
      res(0) = res(0)+ dble(16d0*(amp(1)*dconjg(amp(1))+amp(2)*dconjg(amp(2)))
     &                              -4d0*amp(1)*dconjg(amp(2)))/3d0
      endif ! born
      
      ort(1) = (0d0,0.5d0)*(amp(2)-amp(1))
      ort(2) = 0.5d0*(amp(1)+amp(2))
      ort(3) = (amp(1)+amp(2))/6d0
!       pause
!       return


c---------------------------------
c	VIRTUALS start here
c---------------------------------


c  change momentum convention-> all momenta flowing in
      do mu=0,3
	pdia(mu,2) = -pdia(mu,2)
	qdia(mu,1) = -qdia(mu,1) ! qdia: gluon
	qdia(mu,2) = -qdia(mu,2)
	qdia(mu,3) = -qdia(mu,3)
	kdia(mu,1) = -momW(mu,1) ! kdia: weak boson
	kdia(mu,2) = -momV(mu,1)
	kdia(mu,3) = -momWV(mu,1)
	qdiac(mu,1) = qdia(mu,1)
	qdiac(mu,2) = qdia(mu,2)
	qdiac(mu,3) = qdia(mu,3)
	kdiac(mu,1) = kdia(mu,1)
	kdiac(mu,2) = kdia(mu,2)
	kdiac(mu,3) = kdia(mu,3)
      enddo
      if(ldebug) then
	print*,'mom cons.', (pdia(mu,1)+pdia(mu,2)+qdia(mu,1)+qdia(mu,2)+kdia(mu,1)+kdia(mu,2),mu=0,3)
	print*,'mom cons.', (pdia(mu,1)+pdia(mu,2)+qdia(mu,3)+kdia(mu,3),mu=0,3)
	gaugeprecision = 0d0
      endif
      

c ---------------------------------------
c $$$$$$$$$$ SIMPLE LOOPS $$$$$$$$$$$$$$$
c ---------------------------------------

      if(lsimple) then
        gaugeOK(1) = .true.
      cpl(1) = 3 ! light up quarks
      cpl(2) = 4 ! light down quarks
      cpl(3) = 3 ! top
      
      nfi(1) = 2 ! n light up quarks
      if(NFl .eq.4) then
	nfi(2) = 2
      else
	nfi(2) = 3
      endif
      nfi(3) = 1

      mf(1) = 0d0
      mf(2) = 0d0
      mf(3) = tmass

      if(.not.(lnloW.or.lnloWV)) then
	print*,'check lnloW nlnoWV switch'
	stop
      endif

c glu self enregie
      tempborn=dotcc(epsG(0,3),gluWV)*GG(1)**2 ! 1/q^2 included in epsG
      call gluonSelf(qdia(4,3),musq,NFl,tempc,div)
      tempvirt = tempc*tempborn

      if(ltop) then
	call gluonSelfMass(qdia(4,3),musq,tmass,tempc,div)
	tempvirt = tempvirt + tempc*tempborn
      endif

      ampv(1) = ampv(1)-tempvirt
      ampv(2) = ampv(2)+tempvirt
      
      
c top loop for external gluons
c use on shell renormalisation -> =0
cc      if(ltop) then
cc	tempvirt =0d0
cc	call gluonSelfMass(qdia(4,1),musq,tmass,tempc,div)
cc	tempvirt = tempvirt + tempc
cc	call gluonSelfMass(qdia(4,2),musq,tmass,tempc,div)
cc	tempvirt = tempvirt + tempc
cc	tempvirt=tempvirt/2d0
cc	ampv(1) = ampv(1) + amp(1)*tempvirt
cc	ampv(2) = ampv(2) + amp(2)*tempvirt
cc      endif
	
      
c V- and U-loop contributions to 3 gluon vertex
      do mu=0,3
        momtemp(mu) = -qdia(mu,3)
      enddo
c!       call Vert3G(qdia(0,1),qdia(0,2),momtemp(0),epsG(0,1),epsG(0,2),gluWV,musq,tempborn,-1)
c!       tempborn = tempborn*GG(1)**2/qdia(4,3)
c!       print*,'test',tempborn
      
      call Vert3G(qdia(0,1),qdia(0,2),momtemp(0),epsG(0,1),epsG(0,2),gluWV,musq,tempvirt,div)
      tempvirt=tempvirt*GG(1)**2/qdia(4,3)

      ampv(1) = ampv(1)-tempvirt
      ampv(2) = ampv(2)+tempvirt


c glu glu glu* ! 2 directions have opposite sign but different color structure
      call FermionVertCouplingT(1d0,1d0,1d0,1d0,1d0,1d0,
     $		0d0,
     $		qdia(0,1),qdia(0,2),momtemp(0),
     $		epsG(0,1),epsG(0,2),gluWV,
     $		musq,1,1,floopres(1),floopgauge(1),div)
      tempvirt = -floopres(1)/2d0/qdia(4,3)*G**2
c see comments in FermionVertEvenTemp for minus sign
c color structure is Tr[ta,tb,tc]tc=tatb/2-delta_ab/12
c routine with glu 1 and 2 interchanged will have different sign and color structure tbta/2-delta_ab/8
c the delta_ab color structure will cancel
      ampv(2) = ampv(2)+tempvirt*NFl ! 5 light quarks
      ampv(1) = ampv(1)-tempvirt*NFl

      if(floopgauge(1).gt.GaugeAccuracy) then
        if(printGaugeFail) print*, 'gauge test failed in qqgg, line',__LINE__,floopgauge(1)
        gaugeOK(1) = .false.
      endif

      if(ltop) then
        call FermionVertCouplingT(1d0,1d0,1d0,1d0,1d0,1d0,
     $		tmass,
     $		qdia(0,1),qdia(0,2),momtemp(0),
     $		epsG(0,1),epsG(0,2),gluWV,
     $		musq,1,1,floopres(1),floopgauge(1),div)

        tempvirt = -floopres(1)/2d0/qdia(4,3)*G**2
        ampv(2) = ampv(2)+tempvirt
        ampv(1) = ampv(1)-tempvirt

        if(floopgauge(1).gt.GaugeAccuracy) then
          if(printGaugeFail) print*, 'gauge test failed in qqgg, line',__LINE__,floopgauge(1)
          gaugeOK(1) = .false.
        endif
      endif

      if(lnloWV) then
c glu3 glu*,Z
      do mu=0,3
        momtemp(mu) = -qdia(mu,3)-kdia(mu,2)
      enddo
      momtemp(4) = mass2(momtemp(0))
      if(NFl.eq.5) then
      do i=2,3 ! only b and t contribute
        call FermionVertCouplingT(1d0,1d0,clr(cpl(i),2,+1),1d0,1d0,clr(cpl(i),2,-1),
     $		mf(i),
     $		momtemp,qdia(0,3),kdia(0,2),
     $		gluW,epsG(0,3),epsZ(0,lephel),
     $		musq,1,1,floopres(1),floopgauge(1),div)
        call FermionVertCouplingT(1d0,1d0,clr(cpl(i),2,+1),1d0,1d0,clr(cpl(i),2,-1),
     $		mf(i),
     $		qdia(0,3),momtemp,kdia(0,2),
     $		epsG(0,3),gluW,epsZ(0,lephel),
     $		musq,1,1,floopres(2),floopgauge(2),div)
        tempvirt =  -(floopres(1)+floopres(2))*TR/momtemp(4)*G**2  ! should have a minus sign (see comment in FermionVertEventtempCoupling.F)
        ampv(2) = ampv(2) + tempvirt ! epsG has colstr [ta,tb]
        ampv(1) = ampv(1) - tempvirt

        if(floopgauge(1).gt.GaugeAccuracy.or.floopgauge(2).gt.GaugeAccuracy) then
          if(printGaugeFail) print*, 'gauge test failed in qqgg, line',__LINE__,floopgauge(1),floopgauge(4)
          gaugeOK(1) = .false.
        endif
      enddo
      endif
	  

c glu Z glu* triangles
c (gZg_mu).(qqWg)^mu,  (gZg_mu).glu.(qqW)^nu
      do glu1=1,2
        glu2=3-glu1
	
c	combine gluW and glu2
      do mu =0,3
        momtemp(mu) = pdia(mu,1)+pdia(mu,2)+kdia(mu,1)
      enddo
      momtemp(4) = mass2(momtemp(0))
      do mu = 0,3
c       momenta are defined as incoming -> different sign compared to epsG3
        glutemp(mu) = -(2d0*dotrc(qdia(0,glu2),gluW)*epsG(mu,glu2)
     &		      -2d0*dotrc(momtemp(0),epsG(0,glu2))*gluW(mu)
     &		      +dotcc(gluW(0),epsG(0,glu2))*(momtemp(mu)-qdia(mu,glu2)))/momtemp(4) ! * [t_glu2,t_glu1]
      enddo
      do mu=0,3
        momtemp(mu) =momtemp(mu) + qdia(mu,glu2)
      enddo
      momtemp(4) = mass2(momtemp(0))

      if(NFl.eq.5) then
      do i=2,3 ! only b and t contributions
c 	triangle(Z,glu1) 888 gluWG(glu2)
        call FermionVertCouplingT(1d0,1d0,clr(cpl(i),2,+1),1d0,1d0,clr(cpl(i),2,-1),
     $		mf(i),
     $		qdia(0,glu1),momtemp(0),kdia(0,2),
     $		epsG(0,glu1),gluWG(0,glu2,1),epsZ(0,lephel),
     $		musq,1,1,floopres(1),floopgauge(1),div)

        call FermionVertCouplingT(1d0,1d0,clr(cpl(i),2,+1),1d0,1d0,clr(cpl(i),2,-1),
     $		mf(i),
     $		qdia(0,glu1),momtemp(0),kdia(0,2),
     $		epsG(0,glu1),gluWG(0,glu2,2),epsZ(0,lephel),
     $		musq,0,0,floopres(2),floopgauge(2),div)

c 	triangle(Z,glu1) 888 3gluVert(glu2) 888 gluW
        call FermionVertCouplingT(1d0,1d0,clr(cpl(i),2,+1),1d0,1d0,clr(cpl(i),2,-1),
     $		mf(i),
     $		qdia(0,glu1),momtemp(0),kdia(0,2),
     $		epsG(0,glu1),glutemp(0),epsZ(0,lephel),
     $		musq,0,0,floopres(3),floopgauge(3),div)

c 	triangle(Z,glu1) 888 gluWG(glu2)
        call FermionVertCouplingT(1d0,1d0,clr(cpl(i),2,+1),1d0,1d0,clr(cpl(i),2,-1),
     $		mf(i),
     $		momtemp(0),qdia(0,glu1),kdia(0,2),
     $		gluWG(0,glu2,1),epsG(0,glu1),epsZ(0,lephel),
     $		musq,1,1,floopres(4),floopgauge(4),div)

        call FermionVertCouplingT(1d0,1d0,clr(cpl(i),2,+1),1d0,1d0,clr(cpl(i),2,-1),
     $		mf(i),
     $		momtemp(0),qdia(0,glu1),kdia(0,2),
     $		gluWG(0,glu2,2),epsG(0,glu1),epsZ(0,lephel),
     $		musq,0,0,floopres(5),floopgauge(5),div)

c 	triangle(Z,glu1) 888 3gluVert(glu2) 888 gluW
        call FermionVertCouplingT(1d0,1d0,clr(cpl(i),2,+1),1d0,1d0,clr(cpl(i),2,-1),
     $		mf(i),
     $		momtemp(0),qdia(0,glu1),kdia(0,2),
     $		glutemp(0),epsG(0,glu1),epsZ(0,lephel),
     $		musq,0,0,floopres(6),floopgauge(6),div)

        ampv(glu2) = ampv(glu2) - (floopres(1)+floopres(4)-floopres(3)-floopres(6))*TR/momtemp(4)*G**2 ! minus sign -> see comments in FermVertCoupTemp.F
        ampv(glu1) = ampv(glu1) - (floopres(2)+floopres(5)+floopres(3)+floopres(6))*TR/momtemp(4)*G**2

        if(floopgauge(1).gt.GaugeAccuracy.or.floopgauge(4).gt.GaugeAccuracy) then
          if(printGaugeFail) print*, 'gauge test failed in qqgg, line',__LINE__,floopgauge(1),floopgauge(4)
          gaugeOK(1) = .false.
        endif
      enddo
      endif
      enddo !glu1

c glu glu glu* Z 

      do i=1,3
      if(i.ne.3.or.ltop) then

        do mu=0,3
          momtemp(mu) = -qdia(mu,3)-kdia(mu,2)
        enddo
        momtemp(4) = mass2(momtemp(0))

	wardfb(1)=.false.
        wardfb(2)=.false.
c        gaugeOK(1)=.true.
        gaugeOK(2)=.true.

        call FermionBoxCouplingT(1d0,1d0,1d0,clr(cpl(i),1,+1),1d0,1d0,1d0,clr(cpl(i),1,-1),
     $		mf(i),
     $		qdia(0,1),qdia(0,2),momtemp,kdia(0,2),
     $		epsG(0,1),epsG(0,2),gluW,epsA(0,lephel),
     $		musq,1,1,floopres(1),floopgauge(1),div)
        gaugeOK(1)= wardfb(1).and.gaugeOK(1)
        gaugeOK(2)= wardfb(2).and.gaugeOK(2)

        call FermionBoxCouplingT(1d0,1d0,1d0,clr(cpl(i),2,+1),1d0,1d0,1d0,clr(cpl(i),2,-1),
     $		mf(i),
     $		qdia(0,1),qdia(0,2),momtemp,kdia(0,2),
     $		epsG(0,1),epsG(0,2),gluW,epsZ(0,lephel),
     $		musq,0,0,floopres(2),floopgauge(2),div)

        call FermionBoxCouplingT(1d0,1d0,clr(cpl(i),1,+1),1d0,1d0,1d0,clr(cpl(i),1,-1),1d0,
     $		mf(i),
     $		qdia(0,1),qdia(0,2),kdia(0,2),momtemp,
     $		epsG(0,1),epsG(0,2),epsA(0,lephel),gluW,
     $		musq,1,1,floopres(3),floopgauge(3),div)
        gaugeOK(1)= wardfb(1).and.gaugeOK(1)
        gaugeOK(2)= wardfb(2).and.gaugeOK(2)

        call FermionBoxCouplingT(1d0,1d0,clr(cpl(i),2,+1),1d0,1d0,1d0,clr(cpl(i),2,-1),1d0,
     $		mf(i),
     $		qdia(0,1),qdia(0,2),kdia(0,2),momtemp,
     $		epsG(0,1),epsG(0,2),epsZ(0,lephel),gluW,
     $		musq,0,0,floopres(4),floopgauge(4),div)

        call FermionBoxCouplingT(1d0,clr(cpl(i),1,+1),1d0,1d0,1d0,clr(cpl(i),1,-1),1d0,1d0,
     $		mf(i),
     $		qdia(0,1),kdia(0,2),qdia(0,2),momtemp,
     $		epsG(0,1),epsA(0,lephel),epsG(0,2),gluW,
     $		musq,1,1,floopres(5),floopgauge(5),div)
        gaugeOK(1)= wardfb(1).and.gaugeOK(1)
        gaugeOK(2)= wardfb(2).and.gaugeOK(2)

        call FermionBoxCouplingT(1d0,clr(cpl(i),2,+1),1d0,1d0,1d0,clr(cpl(i),2,-1),1d0,1d0,
     $		mf(i),
     $		qdia(0,1),kdia(0,2),qdia(0,2),momtemp,
     $		epsG(0,1),epsZ(0,lephel),epsG(0,2),gluW,
     $		musq,0,0,floopres(6),floopgauge(6),div)

        if(ldebug) gaugeprecision = max(gaugeprecision,floopgauge(1))
        if(ldebug) gaugeprecision = max(gaugeprecision,floopgauge(3))
        if(ldebug) gaugeprecision = max(gaugeprecision,floopgauge(5))
        if(printgaugeprec) print*,"Fermion Box Gauge: ", floopgauge(1),floopgauge(3),floopgauge(5)
        tempvirt = (floopres(1)+floopres(2)+floopres(3)+floopres(4)+floopres(5)+floopres(6))/momtemp(4)*G**2 
cccccccccccccccccccccccc
cccccccccccccccccccccccc
c   GGGA test
!       call FermionBoxGGGA(qdia(0,1),qdia(0,2),momtemp,kdia(0,2),epsG(0,1),epsG(0,2),gluW,epsA(0,lephel),musq,cpl(i),mf(i),tempvirt) ! GGGA
!       tempvirt = tempvirt*G / momtemp(4) ! 2nd G is included in FBox      call FermionBoxGGGA(qdia(0,1),qdia(0,2),momtemp,kdia(0,2),epsG(0,1),epsG(0,2),gluW,epsA(0,lephel),musq,cpl(i),tempvirt) ! GGGA
! 
!       call FermionBoxGGGZ(qdia(0,1),qdia(0,2),momtemp,kdia(0,2),epsG(0,1),epsG(0,2),gluW,epsZ(0,lephel),musq,cpl(i),mf(i),tempvirt) ! GGGA
!       tempvirt = tempvirt*G / momtemp(4) ! 2nd G is included in FBox
c   GGGA test
cccccccccccccccccccccccc
cccccccccccccccccccccccc
c color structure is Tr[ta,tb,tc]tc=tatb/2-delta_ab/12
        ampv(2) = ampv(2)+tempvirt*nfi(i)/2d0
        ampv(3) = ampv(3)-tempvirt*nfi(i)/12d0

        if(max(floopgauge(1),floopgauge(3),floopgauge(5)).gt.GaugeAccuracy) then
          if(printGaugeFail) print*, 'gauge test failed in qqgg, line',__LINE__,floopgauge(1),floopgauge(3),floopgauge(5)
          gaugeOK(1) = .false.
        endif

c change glu1 and glu2
        call FermionBoxCouplingT(1d0,1d0,1d0,clr(cpl(i),1,+1),1d0,1d0,1d0,clr(cpl(i),1,-1),
     $		mf(i),
     $		qdia(0,2),qdia(0,1),momtemp,kdia(0,2),
     $		epsG(0,2),epsG(0,1),gluW,epsA(0,lephel),
     $		musq,1,1,floopres(1),floopgauge(1),div)

        gaugeOK(1)= wardfb(1).and.gaugeOK(1)
        gaugeOK(2)= wardfb(2).and.gaugeOK(2)

        call FermionBoxCouplingT(1d0,1d0,1d0,clr(cpl(i),2,+1),1d0,1d0,1d0,clr(cpl(i),2,-1),
     $		mf(i),
     $		qdia(0,2),qdia(0,1),momtemp,kdia(0,2),
     $		epsG(0,2),epsG(0,1),gluW,epsZ(0,lephel),
     $		musq,0,0,floopres(2),floopgauge(2),div)

        call FermionBoxCouplingT(1d0,1d0,clr(cpl(i),1,+1),1d0,1d0,1d0,clr(cpl(i),1,-1),1d0,
     $		mf(i),
     $		qdia(0,2),qdia(0,1),kdia(0,2),momtemp,
     $		epsG(0,2),epsG(0,1),epsA(0,lephel),gluW,
     $		musq,1,1,floopres(3),floopgauge(3),div)
        gaugeOK(1)= wardfb(1).and.gaugeOK(1)
        gaugeOK(2)= wardfb(2).and.gaugeOK(2)

        call FermionBoxCouplingT(1d0,1d0,clr(cpl(i),2,+1),1d0,1d0,1d0,clr(cpl(i),2,-1),1d0,
     $		mf(i),
     $		qdia(0,2),qdia(0,1),kdia(0,2),momtemp,
     $		epsG(0,2),epsG(0,1),epsZ(0,lephel),gluW,
     $		musq,0,0,floopres(4),floopgauge(4),div)

        call FermionBoxCouplingT(1d0,clr(cpl(i),1,+1),1d0,1d0,1d0,clr(cpl(i),1,-1),1d0,1d0,
     $		mf(i),
     $		qdia(0,2),kdia(0,2),qdia(0,1),momtemp,
     $		epsG(0,2),epsA(0,lephel),epsG(0,1),gluW,
     $		musq,1,1,floopres(5),floopgauge(5),div)
        gaugeOK(1)= wardfb(1).and.gaugeOK(1)
        gaugeOK(2)= wardfb(2).and.gaugeOK(2)

        call FermionBoxCouplingT(1d0,clr(cpl(i),2,+1),1d0,1d0,1d0,clr(cpl(i),2,-1),1d0,1d0,
     $		mf(i),
     $		qdia(0,2),kdia(0,2),qdia(0,1),momtemp,
     $		epsG(0,2),epsZ(0,lephel),epsG(0,1),gluW,
     $		musq,0,0,floopres(6),floopgauge(6),div)

        if(ldebug) gaugeprecision = max(gaugeprecision,floopgauge(1))
        if(ldebug) gaugeprecision = max(gaugeprecision,floopgauge(3))
        if(ldebug) gaugeprecision = max(gaugeprecision,floopgauge(5))
        if(printgaugeprec) print*,"Fermion Box Gauge: ", floopgauge(1),floopgauge(3),floopgauge(5)
        tempvirt = (floopres(1)+floopres(2)+floopres(3)+floopres(4)+floopres(5)+floopres(6))/momtemp(4)*G**2

cccccccccccccccccccccccc
cccccccccccccccccccccccc
c   GGGA test
!       call FermionBoxGGGA(qdia(0,2),qdia(0,1),momtemp,kdia(0,2),epsG(0,2),epsG(0,1),gluW,epsA(0,lephel),musq,cpl(i),mf(i),tempvirt) ! GGGA
!       tempvirt = tempvirt*G / momtemp(4) ! 2nd G is included in FBox
! 
!       call FermionBoxGGGZ(qdia(0,2),qdia(0,1),momtemp,kdia(0,2),epsG(0,2),epsG(0,1),gluW,epsZ(0,lephel),musq,cpl(i),mf(i),tempvirt) ! GGGA
!       tempvirt = tempvirt*G / momtemp(4) ! 2nd G is included in FBox
c   GGGA test
cccccccccccccccccccccccc
cccccccccccccccccccccccc
c color structure is Tr[tb,ta,tc]tc=tbta/2-delta_ab/12
        ampv(1) = ampv(1)+tempvirt*nfi(i)/2d0
        ampv(3) = ampv(3)-tempvirt*nfi(i)/12d0

        if(max(floopgauge(1),floopgauge(3),floopgauge(5)).gt.GaugeAccuracy) then
          if(printGaugeFail) print*, 'gauge test failed in qqgg, line',__LINE__,floopgauge(1),floopgauge(3),floopgauge(5)
          gaugeOK(1) = .false.
        endif
       endif
       enddo
      endif ! lnloWV
      
      box2count2m = box2count2m + 1
      box2count2sub(sub_number) = box2count2sub(sub_number) + 1
      if(.not.gaugeOK(1)) then
        box2countm   = box2countm   + 1
        box2countsub  (sub_number) = box2countsub  (sub_number) + 1
      endif
      if(.not.gaugeOK(2)) then
        box2countm_Q = box2countm_Q + 1
        box2countsub_q(sub_number) = box2countsub_q(sub_number) + 1
        if(div.eq.0.and.return0ifGaugeFails) then
          res(0) = 0d0; res(1) = 0d0; res(2) = 0d0
          return
        else
          recalculatePSpoint = .true.
        endif
      endif
      

      endif ! lsimple
      
      if(printall) then
	print*,'ampv,simple',ampv
      endif
     
c ---------------------------------------
c $$$$$$$$$$ B O X E S $$$$$$$$$$$$$$$$$$
c ---------------------------------------

      if(lbox.and.lnloW) then

      if(ldebug) then
        print*
        print*,'box'
      endif
      boxborn(1) = 0d0
      boxborn(2) = 0d0
      boxborn(3) = 0d0

c q> GG > WV> q
      call boxlineABETotal(pdia(0,1),qdia(0,3),kdia(0,3),pdia(0,2),
     &            psi(1,-1,2),psi(1,-1,1),epsG(0,3),epsWV(0,lephel,1),
     &            -1,musq,1,2,3,1,boxgauge(1,1),boxres(1,1),boxgaugeb(1,1),boxborn(1),div)

c nA gauge test
      call boxlineABETotal(pdia(0,1),qdia(0,3),kdia(0,3),pdia(0,2),
     &            psi(1,-1,2),psi(1,-1,1),qdiac(0,3),epsWV(0,lephel,1),
     &            -1,musq,1,2,9,-1,nothing(1),checkg(1,1,1),nothing(1),checkgb(1,1),div)

c q> WV> GG> q
      call boxlineABETotal(pdia(0,1),kdia(0,3),qdia(0,3),pdia(0,2),
     &            psi(1,-1,2),psi(1,-1,1),epsWV(0,lephel,1),epsG(0,3),
     &            -1,musq,1,3,3,1,boxgauge(1,2),boxres(1,2),boxgaugeb(1,2),boxborn(2),div)

c nA gauge test
      call boxlineABETotal(pdia(0,1),kdia(0,3),qdia(0,3),pdia(0,2),
     &            psi(1,-1,2),psi(1,-1,1),epsWV(0,lephel,1),qdiac(0,3),
     &            -1,musq,1,3,9,-1,nothing(1),checkg(1,2,1),nothing(1),checkgb(2,1),div)

c non Abelian
      call BoxlineNoAbeTotal(pdia(0,1),kdia(0,3),pdia(0,2),qdia(0,3),
     &            psi(1,-1,2),psi(1,-1,1),epsWV(0,lephel,1),epsG(0,3),
     &            -1,musq,2,1,boxgauge(1,3),boxres(1,3),div)
      
c gauge tests 
      call box_gauge(pdia(0,1),qdia(0,3),kdia(0,3),pdia(0,2),musq,
     2    boxgauge(1,1),boxgaugeb(1,1),wardidtest(1))
      call box_gauge(pdia(0,1),kdia(0,3),qdia(0,3),pdia(0,2),musq,
     2    boxgauge(1,2),boxgaugeb(1,2),wardidtest(2))

c      check epsG -> pG
      if (abs(boxgauge(1,3)).lt.1d-6) then
        td=abs(-(checkg(2,2,1)+checkg(2,1,1))/2d0+boxgauge(1,3))
      else
        td=abs(-(checkg(2,2,1)+checkg(2,1,1))/2d0/boxgauge(1,3)+1d0)
      endif
      wardidtest(3) = td.lt.GaugeAccuracy

      if(ldebug) gaugeprecision = max(gaugeprecision,td)
      if(printgaugeprec) print*,"BoxLine Gauge: ", td
      if(ldebug) print*,wardidtest(1),wardidtest(2),wardidtest(3),td,checkg(1,3,1)
c sum
      wardidtest(0) = wardidtest(1).and.wardidtest(2).and.wardidtest(3)
      wardidtest(0) = wardidtest(0).or.ldebug.or.div.ne.0

      if(printAll) then
	print*,'box'
	print*,wardidtest
	print*,boxres
	print*,boxborn
      endif

      boxcount2m = boxcount2m + 1
      boxcount2sub(sub_number) = boxcount2sub(sub_number) + 1

      tempborn = boxborn(1)+boxborn(2)
      tempvirt = CF*(boxres(1,1)+boxres(1,2))+CFCA*(boxres(2,1)+boxres(2,2))+CA*boxres(1,3)
      tempborn = -tempborn*G**2 ! sign boxline routine differs from Peskin-S convetnion
      tempvirt = -tempvirt*G**2
      
      ampb(1) = ampb(1) - tempborn ! [ta,tb] = -col1 + col2
      ampb(2) = ampb(2) + tempborn
      if(div.eq.0) then
        ampv(1) = ampv(1)+CF*tempborn ! rational term
        ampv(2) = ampv(2)-CF*tempborn ! rational term
      endif
      ampv(1) = ampv(1) - tempvirt
      ampv(2) = ampv(2) + tempvirt

      if(.not.wardidtest(0)) then
        boxcountm   = boxcountm   + 1
        boxcountm_q = boxcountm_q + 1
        boxcountsub  (sub_number) = boxcountsub  (sub_number) + 1
        boxcountsub_q(sub_number) = boxcountsub_q(sub_number) + 1
        if(printGaugeFail) print*, 'gauge test failed in qqgg, line',__LINE__,wardidtest,td
        if(return0ifGaugeFails) then
          res(0) = 0d0; res(1) = 0d0; res(2) = 0d0
          return
        else
          recalculatePSpoint=.true.
        endif
      endif

      
      endif ! lbox
      
c ---------------------------------------
c $$$$$$$$$$ PENTAGONS $$$$$$$$$$$$$$$$$$
c ---------------------------------------

      if(lpent) then
        gaugeOK(1) = .true.
        gaugeOK(2) = .true.
      if(ldebug) then
        print*
        print*,'pent'
      endif

      if(lnloWV) then

      wardpAbe(1)=.false.
      wardpAbe(2)=.false.
      wardpNoAbe(1)=.false.
      wardpNoAbe(2)=.false.
      wardhAbe(1)=.false.
      wardhAbe(2)=.false.
      wardhNoAbe(1)=.false.
      wardhNoAbe(2)=.false.

c gg, W, V type      

c q > gg> W> V> q
      call penlineAbeTotalQ(pdia(0,1),qdia(0,3),kdia(0,1),kdia(0,2),pdia(0,2),
     &		psi(1,-1,2),psi(1,-1,1),epsG(0,3),epsW(0,1),epsV(0,lephel,qtype(2),-1,1),
     &		-1,musq,1,2,4,1,pengauge(1,1),penres(1,1),pengaugeb(1,1),penborn(1),div)!,wardidtest(1))
      wardidtest(1)=wardpAbe(1)
      wardidtest_Q(1)=wardpAbe(2)    
c q > W> gg> V> q
      call penlineAbeTotalQ(pdia(0,1),kdia(0,1),qdia(0,3),kdia(0,2),pdia(0,2),
     &		psi(1,-1,2),psi(1,-1,1),epsW(0,1),epsG(0,3),epsV(0,lephel,qtype(2),-1,1),
     &		-1,musq,1,3,4,1,pengauge(1,2),penres(1,2),pengaugeb(1,2),penborn(2),div)!,wardidtest(2))
      wardidtest(2)=wardpAbe(1)
      wardidtest_Q(2)=wardpAbe(2)    
c q > W> V> gg> q
      call penlineAbeTotalQ(pdia(0,1),kdia(0,1),kdia(0,2),qdia(0,3),pdia(0,2),
     &		psi(1,-1,2),psi(1,-1,1),epsW(0,1),epsV(0,lephel,qtype(2),-1,1),epsG(0,3),
     &		-1,musq,1,4,4,1,pengauge(1,3),penres(1,3),pengaugeb(1,3),penborn(3),div)!,wardidtest(3))
      wardidtest(3)=wardpAbe(1)
      wardidtest_Q(3)=wardpAbe(2)    
c q > gg> V> W> q
      call penlineAbeTotalQ(pdia(0,1),qdia(0,3),kdia(0,2),kdia(0,1),pdia(0,2),
     &		psi(1,-1,2),psi(1,-1,1),epsG(0,3),epsV(0,lephel,qtype(1),-1,1),epsW(0,1),
     &		-1,musq,1,2,4,1,pengauge(1,4),penres(1,4),pengaugeb(1,4),penborn(4),div)!,wardidtest(4))
      wardidtest(4)=wardpAbe(1)
      wardidtest_Q(4)=wardpAbe(2)    
c q > V> gg> W> q
      call penlineAbeTotalQ(pdia(0,1),kdia(0,2),qdia(0,3),kdia(0,1),pdia(0,2),
     &		psi(1,-1,2),psi(1,-1,1),epsV(0,lephel,qtype(1),-1,1),epsG(0,3),epsW(0,1),
     &		-1,musq,1,3,4,1,pengauge(1,5),penres(1,5),pengaugeb(1,5),penborn(5),div)!,wardidtest(5))
      wardidtest(5)=wardpAbe(1)
      wardidtest_Q(5)=wardpAbe(2)    
c q > V> W> gg> q
      call penlineAbeTotalQ(pdia(0,1),kdia(0,2),kdia(0,1),qdia(0,3),pdia(0,2),
     &		psi(1,-1,2),psi(1,-1,1),epsV(0,lephel,qtype(1),-1,1),epsW(0,1),epsG(0,3),
     &		-1,musq,1,4,4,1,pengauge(1,6),penres(1,6),pengaugeb(1,6),penborn(6),div)!,wardidtest(6))
      wardidtest(6)=wardpAbe(1)
      wardidtest_Q(6)=wardpAbe(2)    

      wardpAbe(1)=.false.
      wardpAbe(2)=.false.
      wardpNoAbe(1)=.false.
      wardpNoAbe(2)=.false.
      wardhAbe(1)=.false.
      wardhAbe(2)=.false.
      wardhNoAbe(1)=.false.
      wardhNoAbe(2)=.false.


c nonAbelian
c q > W> V> q
      call penlineNoAbeTotalQ(pdia(0,1),kdia(0,1),kdia(0,2),pdia(0,2),qdia(0,3),
     &          psi(1,-1,2),psi(1,-1,1),epsW(0,1),epsV(0,lephel,qtype(2),-1,1),epsG(0,3),
     &          -1,musq,4,1,pengauge(1,7),penres(1,7),div)!,wardidtest(7))
      wardidtest(7)=wardpNoAbe(1)  .or.wardpNoAbe(2)    
      wardidtest_Q(7)=wardpNoAbe(2)  
c q > V> W> q
      call penlineNoAbeTotalQ(pdia(0,1),kdia(0,2),kdia(0,1),pdia(0,2),qdia(0,3),
     &          psi(1,-1,2),psi(1,-1,1),epsV(0,lephel,qtype(1),-1,1),epsW(0,1),epsG(0,3),
     &          -1,musq,4,1,pengauge(1,8),penres(1,8),div)!,wardidtest(8))
      wardidtest(8)=wardpNoAbe(1)  .or.wardpNoAbe(2)    
      wardidtest_Q(8)=wardpNoAbe(2)  

c      if(ldebug) gaugeprecision = max(gaugeprecision,td)


      wardidtest(0)=wardidtest(1).and.wardidtest(2).and.wardidtest(3).and.wardidtest(4).and.wardidtest(5)
     &              .and.wardidtest(6).and.wardidtest(7).and.wardidtest(8)
      wardidtest_q(0)=wardidtest_q(1).and.wardidtest_q(2).and.wardidtest_q(3).and.wardidtest_q(4).and.wardidtest_q(5)
     &              .and.wardidtest_q(6).and.wardidtest_q(7).and.wardidtest_q(8)


      if(ldebug) print*,wardidtest
      if(ldebug) print*,wardidtest_q
      wardidtest_q(0) = wardidtest_q(0).or.ldebug.or.div.ne.0

      if(printAll) then
	print*,'pen',0
	print*,wardidtest
	print*,wardidtest_q
	print*,penres
	print*,penborn
      endif

      gaugeOK(1) = gaugeOK(1).and.wardidtest(0)
      gaugeOK(2) = gaugeOK(2).and.wardidtest_q(0)
c sum
      tempborn = penborn(1)+penborn(2)+penborn(3)+penborn(4)+penborn(5)+penborn(6)
      tempvirt =   CF*(penres(1,1)+penres(1,2)+penres(1,3)+penres(1,4)+penres(1,5)+penres(1,6))
     &  	          +CFCA*(penres(2,1)+penres(2,2)+penres(2,3)+penres(2,4)+penres(2,5)+penres(2,6))
     &  	          +CA*(penres(1,7)+penres(1,8))

      tempborn = -tempborn*G**2 ! sign penline routine differs from Peskin-S convetnion
      tempvirt = -tempvirt*G**2
      
      ampb(1) = ampb(1) - tempborn ! [ta,tb] is col2-col1
      ampb(2) = ampb(2) + tempborn
      if(div.eq.0) then
        ampv(1) = ampv(1)+CF*tempborn ! rational term
        ampv(2) = ampv(2)-CF*tempborn ! rational term
      endif
      ampv(1) = ampv(1) - tempvirt
      ampv(2) = ampv(2) + tempvirt

      if(.not.wardidtest_Q(0)) then
        if(printGaugeFail) print*, 'gauge test failed in qqgg, line',__LINE__,wardidtest,td
      endif
      
      endif ! lnloWV

      if(lnloW) then

      tempvirt=0d0
      tempvirt2=0d0
      tempvirt3=0d0
      tempvirt4=0d0
      tempborn=0d0
c g,g,WV type    
      do glu1=1,2
      glu2=3-glu1
      
      do i=1,8; do j=1,8
        penres(i,j)=0d0
      enddo; enddo

      wardpAbe(1)=.false.
      wardpAbe(2)=.false.
      wardpNoAbe(1)=.false.
      wardpNoAbe(2)=.false.
      wardhAbe(1)=.false.
      wardhAbe(2)=.false.
      wardhNoAbe(1)=.false.
      wardhNoAbe(2)=.false.

      
c glu1,glu2,WV
      call penlineAbeTotalQ(pdia(0,1),qdia(0,glu1),qdia(0,glu2),kdia(0,3),pdia(0,2),
     &          psi(1,-1,2),psi(1,-1,1),epsG(0,glu1),epsG(0,glu2),epsWV(0,lephel,1),
     &          -1,musq,2,4,4,1,pengauge(1,1),penres(1,1),pengaugeb(1,1),penborn(1),div)!,wardidtest(1))
      wardidtest(1)=wardpAbe(1)
      wardidtest_Q(1)=wardpAbe(2)  
c glu1,WV,glu2

      call penlineAbeTotalQ(pdia(0,1),qdia(0,glu1),kdia(0,3),qdia(0,glu2),pdia(0,2),
     &          psi(1,-1,2),psi(1,-1,1),epsG(0,glu1),epsWV(0,lephel,1),epsG(0,glu2),
     &          -1,musq,2,3,4,1,pengauge(1,2),penres(1,2),pengaugeb(1,2),penborn(2),div)!,wardidtest(2))
      wardidtest(2)=wardpAbe(1)
      wardidtest_Q(2)=wardpAbe(2)  
c WV,glu1,glu2
      call penlineAbeTotalQ(pdia(0,1),kdia(0,3),qdia(0,glu1),qdia(0,glu2),pdia(0,2),
     &          psi(1,-1,2),psi(1,-1,1),epsWV(0,lephel,1),epsG(0,glu1),epsG(0,glu2),
     &          -1,musq,2,2,4,1,pengauge(1,3),penres(1,3),pengaugeb(1,3),penborn(3),div)!,wardidtest(3))
      wardidtest(3)=wardpAbe(1)
      wardidtest_Q(3)=wardpAbe(2)  

      wardpAbe(1)=.false.
      wardpAbe(2)=.false.
      wardpNoAbe(1)=.false.
      wardpNoAbe(2)=.false.
      wardhAbe(1)=.false.
      wardhAbe(2)=.false.
      wardhNoAbe(1)=.false.
      wardhNoAbe(2)=.false.


c nonAbe glu2
c glu1,WV
      call penlineNoAbe_pTQ(pdia(0,1),qdia(0,glu1),kdia(0,3),pdia(0,2),qdia(0,glu2),
     &          psi(1,-1,2),psi(1,-1,1),epsG(0,glu1),epsWV(0,lephel,1),epsG(0,glu2),
     &          -1,musq,2,2,4,1,pengauge(1,4),penres(1,4),div)!,wardidtest(4))
      wardidtest(4)=wardpNoAbe(1)
      wardidtest_Q(4)=wardpNoAbe(2)  
c WV,glu1
      call penlineNoAbe_pTQ(pdia(0,1),kdia(0,3),qdia(0,glu1),pdia(0,2),qdia(0,glu2),
     &          psi(1,-1,2),psi(1,-1,1),epsWV(0,lephel,1),epsG(0,glu1),epsG(0,glu2),
     &          -1,musq,2,3,4,1,pengauge(1,5),penres(1,5),div)!,wardidtest(5))
      wardidtest(5)=wardpNoAbe(1)
      wardidtest_Q(5)=wardpNoAbe(2)  

      wardpAbe(1)=.false.
      wardpAbe(2)=.false.
      wardpNoAbe(1)=.false.
      wardpNoAbe(2)=.false.
      wardhAbe(1)=.false.
      wardhAbe(2)=.false.
      wardhNoAbe(1)=.false.
      wardhNoAbe(2)=.false.


c 2 triple vertices
      call NoAbe2gaPenLineTotalQ(pdia(0,1),kdia(0,3),pdia(0,2),qdia(0,glu2),qdia(0,glu1),
     &		psi(1,-1,2),psi(1,-1,1),epsWV(0,lephel,1),epsG(0,glu2),epsG(0,glu1),
     &		-1,musq,2,0,1,nothing(1),penres(1,6),gaugereal(1,1),div,1)
      wardidtest(6)=wardpNoAbe(1)
      wardidtest_Q(6)=wardpNoAbe(2)  

      wardpAbe(1)=.false.
      wardpAbe(2)=.false.
      wardpNoAbe(1)=.false.
      wardpNoAbe(2)=.false.
      wardhAbe(1)=.false.
      wardhAbe(2)=.false.
      wardhNoAbe(1)=.false.
      wardhNoAbe(2)=.false.

! c quad vertex
      call NoAbe4vertBoxLineTotalQ(pdia(0,1),kdia(0,3),pdia(0,2),qdia(0,glu2),qdia(0,glu1),
     &		psi(1,-1,2),psi(1,-1,1),epsWV(0,lephel,1),epsG(0,glu2),epsG(0,glu1),
     &		-1,musq,2,0,1,nothing(1),penres(1,7),gaugereal(1,2),div,1)

      wardidtest(7)=wardpNoAbe(1)
      wardidtest_Q(7)=wardpNoAbe(2)  
! 

      do i=1,3
        penres(i,7) = penres(i,7)/2d0 ! we calculate this contribution twice, but only need it once
      enddo

      wardidtest(0)=wardidtest(1).and.wardidtest(2).and.wardidtest(3).and.wardidtest(4).and.wardidtest(5)
     &              .and.wardidtest(6).and.wardidtest(7)
      wardidtest_q(0)=wardidtest_q(1).and.wardidtest_q(2).and.wardidtest_q(3).and.wardidtest_q(4).and.wardidtest_q(5)
     &              .and.wardidtest_q(6).and.wardidtest_q(7)


      if(ldebug) print*,wardidtest
      wardidtest_q(0) = wardidtest_q(0).or.ldebug.or.div.ne.0


      if(printAll) then
	print*,'pen',glu1
	print*,wardidtest
	print*,wardidtest_q
	print*,penres
	print*,penborn
      endif

      gaugeOK(1) = gaugeOK(1).and.wardidtest(0)
      gaugeOK(2) = gaugeOK(2).and.wardidtest_q(0)
c sum
      tempborn = -penborn(1)-penborn(2)-penborn(3) ! penlineAbe and penlineNoAbe have different sign convention then Hex*
      tempvirt = -penres(1,1)-penres(1,2)-penres(1,3)-(0,-1.5d0)*
     &          (penres(2,4)+penres(2,5))+( -1.5d0*penres(1,6)
     &          + 1.5d0*(penres(2,7)+penres(3,7)-2d0*penres(1,7)))
      tempvirt3= -penres(2,1)-penres(2,2)-penres(2,3)-(0,0.25d0)
     &          *(penres(1,4)+penres(1,5))+(-0.25d0*penres(1,6)
     &          +0.25d0*(2d0*penres(2,7)-penres(1,7)-penres(3,7)))
      tempvirt2=                                     -(0,-1.5d0)*(penres(3,4)+penres(3,5))+
     &                                1.5d0*(penres(1,7)+penres(2,7)-2d0*penres(3,7))

      tempborn  = tempborn*G**2 
      tempvirt  = tempvirt*G**2
      tempvirt2 = tempvirt2*G**2
      tempvirt3 = tempvirt3*G**2
      
      ampb(glu1) = ampb(glu1) + tempborn ! >1> 2> is colstr 1
      if(div.eq.0) ampv(glu1) = ampv(glu1)-CF*tempborn ! rational term
      ampv(glu1) = ampv(glu1) + tempvirt
      ampv(glu2) = ampv(glu2) + tempvirt2
      ampv(3)    = ampv(3)    + tempvirt3

      if(.not.wardidtest_q(0)) then
        if(printGaugeFail) print*, 'gauge test failed in qqgg, line',__LINE__,(wardidtest(i),i=1,7),td
      endif
      enddo ! glu

      endif ! nlnoW

      pentcount2m = pentcount2m + 1
      pentcount2sub(sub_number) = pentcount2sub(sub_number) + 1
      if(.not.gaugeOK(1)) then
        pentcountm = pentcountm + 1
        pentcountsub(sub_number) = pentcountsub(sub_number) + 1
      endif
      if(.not.gaugeOK(2)) then
        pentcountm_Q = pentcountm_q + 1
        pentcountsub_Q(sub_number) = pentcountsub_Q(sub_number) + 1
        if(return0ifGaugeFails) then
          res(0) = 0d0; res(1) = 0d0; res(2) = 0d0
          return
        else
          recalculatePSpoint=.true.
        endif
      endif
      endif ! lpent


c ---------------------------------------
c $$$$$$$$$$ HEXAGONS $$$$$$$$$$$$$$$$$$$
c ---------------------------------------

      if(lhex.and.lnloWV) then
        gaugeOK(1) = .true.
        gaugeOK(2) = .true.
      if(ldebug) then
        print*
        print*,'hexagons'
      endif
      tempvirt = 0d0
      tempvirt2=0d0
      tempvirt3=0d0
      tempvirt4=0d0

      ngaugeAbe=4*rn+1
      if(ngaugeAbe.eq.5) then
        ngaugeAbe=4
      endif
      ngaugeNonAbe=3*rn+1
      if(ngaugeNonAbe.eq.4) then
        ngaugeNonAbe=3
      endif
      
      do glu1 = 1,2
      glu2 = 3-glu1
      i=1
      
!       print*
      if(ldebug) print*, "glu1 = ", glu1

      wardpAbe(1)=.false.
      wardpAbe(2)=.false.
      wardpNoAbe(1)=.false.
      wardpNoAbe(2)=.false.
      wardhAbe(1)=.false.
      wardhAbe(2)=.false.
      wardhNoAbe(1)=.false.
      wardhNoAbe(2)=.false.
!       print*,'a'

c glu1 glu2 W Z      
      call HexlineAbeTotalQ(pdia(0,1),qdia(0,glu1),qdia(0,glu2),kdia(0,1),kdia(0,2),pdia(0,2),
     &            psi(1,-1,2),psi(1,-1,1),epsG(0,glu1),epsG(0,glu2),epsW(0,1),epsV(0,lephel,qtype(2),-1,1),
     &            -1,musq,2,23,1,hexres(1,1),hexborn(1),hexgauge(1,1),hexgaugeb(1,1),div,nGaugeAbe)!,wardidtest(1))
      wardidtest(1)=wardhAbe(1)
      wardidtest_Q(1)=wardhAbe(2)  
!       print*,'b'

c glu1 W glu2 Z      
      call HexlineAbeTotalQ(pdia(0,1),qdia(0,glu1),kdia(0,1),qdia(0,glu2),kdia(0,2),pdia(0,2),
     &            psi(1,-1,2),psi(1,-1,1),epsG(0,glu1),epsW(0,1),epsG(0,glu2),epsV(0,lephel,qtype(2),-1,1),
     &            -1,musq,2,24,1,hexres(1,2),hexborn(2),hexgauge(1,2),hexgaugeb(1,2),div,nGaugeAbe)!,wardidtest(2))
      wardidtest(2)=wardhAbe(1)
      wardidtest_Q(2)=wardhAbe(2)  
!       print*,'c'
    
c glu1 W Z glu2
      call HexlineAbeTotalQ(pdia(0,1),qdia(0,glu1),kdia(0,1),kdia(0,2),qdia(0,glu2),pdia(0,2),
     &            psi(1,-1,2),psi(1,-1,1),epsG(0,glu1),epsW(0,1),epsV(0,lephel,qtype(2),-1,1),epsG(0,glu2),
     &            -1,musq,2,25,1,hexres(1,3),hexborn(3),hexgauge(1,3),hexgaugeb(1,3),div,nGaugeAbe)!,wardidtest(3))
      wardidtest(3)=wardhAbe(1)
      wardidtest_Q(3)=wardhAbe(2)  
!       print*,'d'

     
c glu1 glu2 Z W
      call HexlineAbeTotalQ(pdia(0,1),qdia(0,glu1),qdia(0,glu2),kdia(0,2),kdia(0,1),pdia(0,2),
     &            psi(1,-1,2),psi(1,-1,1),epsG(0,glu1),epsG(0,glu2),epsV(0,lephel,qtype(1),-1,1),epsW(0,1),
     &            -1,musq,2,23,1,hexres(1,4),hexborn(4),hexgauge(1,4),hexgaugeb(1,4),div,nGaugeAbe)!,wardidtest(4))
      wardidtest(4)=wardhAbe(1)
      wardidtest_Q(4)=wardhAbe(2)  
!       print*,'e'

     
c glu1 Z glu2 W      
      call HexlineAbeTotalQ(pdia(0,1),qdia(0,glu1),kdia(0,2),qdia(0,glu2),kdia(0,1),pdia(0,2),
     &            psi(1,-1,2),psi(1,-1,1),epsG(0,glu1),epsV(0,lephel,qtype(1),-1,1),epsG(0,glu2),epsW(0,1),
     &            -1,musq,2,24,1,hexres(1,5),hexborn(5),hexgauge(1,5),hexgaugeb(1,5),div,nGaugeAbe)!,wardidtest(5))
      wardidtest(5)=wardhAbe(1)
      wardidtest_Q(5)=wardhAbe(2)  
!       print*,'f'

     
c glu1 Z W glu2
      call HexlineAbeTotalQ(pdia(0,1),qdia(0,glu1),kdia(0,2),kdia(0,1),qdia(0,glu2),pdia(0,2),
     &            psi(1,-1,2),psi(1,-1,1),epsG(0,glu1),epsV(0,lephel,qtype(1),-1,1),epsW(0,1),epsG(0,glu2),
     &            -1,musq,2,25,1,hexres(1,6),hexborn(6),hexgauge(1,6),hexgaugeb(1,6),div,nGaugeAbe)!,wardidtest(6))
      wardidtest(6)=wardhAbe(1)
      wardidtest_Q(6)=wardhAbe(2)  
!       print*,'g'

      
c W glu1 glu2 Z      
      call HexlineAbeTotalQ(pdia(0,1),kdia(0,1),qdia(0,glu1),qdia(0,glu2),kdia(0,2),pdia(0,2),
     &            psi(1,-1,2),psi(1,-1,1),epsW(0,1),epsG(0,glu1),epsG(0,glu2),epsV(0,lephel,qtype(2),-1,1),
     &            -1,musq,2,34,1,hexres(1,7),hexborn(7),hexgauge(1,7),hexgaugeb(1,7),div,nGaugeAbe)!,wardidtest(7))     
      wardidtest(7)=wardhAbe(1)
      wardidtest_Q(7)=wardhAbe(2)  
!       print*,'h'


c W glu1 Z glu2
      call HexlineAbeTotalQ(pdia(0,1),kdia(0,1),qdia(0,glu1),kdia(0,2),qdia(0,glu2),pdia(0,2),
     &            psi(1,-1,2),psi(1,-1,1),epsW(0,1),epsG(0,glu1),epsV(0,lephel,qtype(2),-1,1),epsG(0,glu2),
     &            -1,musq,2,35,1,hexres(1,8),hexborn(8),hexgauge(1,8),hexgaugeb(1,8),div,nGaugeAbe)!,wardidtest(8))     
      wardidtest(8)=wardhAbe(1)
      wardidtest_Q(8)=wardhAbe(2)  
!       print*,'i'

     
c W Z glu1 glu2
      call HexlineAbeTotalQ(pdia(0,1),kdia(0,1),kdia(0,2),qdia(0,glu1),qdia(0,glu2),pdia(0,2),
     &            psi(1,-1,2),psi(1,-1,1),epsW(0,1),epsV(0,lephel,qtype(2),-1,1),epsG(0,glu1),epsG(0,glu2),
     &            -1,musq,2,45,1,hexres(1,9),hexborn(9),hexgauge(1,9),hexgaugeb(1,9),div,nGaugeAbe)!,wardidtest(9))
      wardidtest(9)=wardhAbe(1)
      wardidtest_Q(9)=wardhAbe(2)  
!       print*,'j'

     
c Z glu1 glu2 W
      call HexlineAbeTotalQ(pdia(0,1),kdia(0,2),qdia(0,glu1),qdia(0,glu2),kdia(0,1),pdia(0,2),
     &            psi(1,-1,2),psi(1,-1,1),epsV(0,lephel,qtype(1),-1,1),epsG(0,glu1),epsG(0,glu2),epsW(0,1),
     &            -1,musq,2,34,1,hexres(1,10),hexborn(10),hexgauge(1,10),hexgaugeb(1,10),div,nGaugeAbe)!,wardidtest(10))     
      wardidtest(10)=wardhAbe(1)
      wardidtest_Q(10)=wardhAbe(2)  
!       print*,'k'


c Z glu1 W glu2
      call HexlineAbeTotalQ(pdia(0,1),kdia(0,2),qdia(0,glu1),kdia(0,1),qdia(0,glu2),pdia(0,2),
     &            psi(1,-1,2),psi(1,-1,1),epsV(0,lephel,qtype(1),-1,1),epsG(0,glu1),epsW(0,1),epsG(0,glu2),
     &            -1,musq,2,35,1,hexres(1,11),hexborn(11),hexgauge(1,11),hexgaugeb(1,11),div,nGaugeAbe)!,wardidtest(11))     
      wardidtest(11)=wardhAbe(1)
      wardidtest_Q(11)=wardhAbe(2)  
!       print*,'l'

     
c Z W glu1 glu2
      call HexlineAbeTotalQ(pdia(0,1),kdia(0,2),kdia(0,1),qdia(0,glu1),qdia(0,glu2),pdia(0,2),
     &            psi(1,-1,2),psi(1,-1,1),epsV(0,lephel,qtype(1),-1,1),epsW(0,1),epsG(0,glu1),epsG(0,glu2),
     &            -1,musq,2,45,1,hexres(1,12),hexborn(12),hexgauge(1,12),hexgaugeb(1,12),div,nGaugeAbe)!,wardidtest(12))
      wardidtest(12)=wardhAbe(1)
      wardidtest_Q(12)=wardhAbe(2)  
!       print*,'m'


      wardpAbe(1)=.false.
      wardpAbe(2)=.false.
      wardpNoAbe(1)=.false.
      wardpNoAbe(2)=.false.
      wardhAbe(1)=.false.
      wardhAbe(2)=.false.
      wardhNoAbe(1)=.false.
      wardhNoAbe(2)=.false.


c no Abe glu2
c glu1 W Z      
      call HexlineNoAbe1gTotalQ(pdia(0,1),qdia(0,glu1),kdia(0,1),kdia(0,2),pdia(0,2),qdia(0,glu2),
     &            psi(1,-1,2),psi(1,-1,1),epsG(0,glu1),epsW(0,1),epsV(0,lephel,qtype(2),-1,1),epsG(0,glu2),
     &            -1,musq,2,2,1,hexres(1,13),hexgauge(1,13),div,nGaugeNonAbe)!,wardidtest(13))
      wardidtest(13)=wardhNoAbe(1) 
      wardidtest_Q(13)=wardhNoAbe(2)  


c glu1 Z W 
      call HexlineNoAbe1gTotalQ(pdia(0,1),qdia(0,glu1),kdia(0,2),kdia(0,1),pdia(0,2),qdia(0,glu2),
     &            psi(1,-1,2),psi(1,-1,1),epsG(0,glu1),epsV(0,lephel,qtype(1),-1,1),epsW(0,1),epsG(0,glu2),
     &            -1,musq,2,2,1,hexres(1,14),hexgauge(1,14),div,nGaugeNonAbe)!,wardidtest(14))
      wardidtest(14)=wardhNoAbe(1)
      wardidtest_Q(14)=wardhNoAbe(2)  


c W glu1 Z      
      call HexlineNoAbe1gTotalQ(pdia(0,1),kdia(0,1),qdia(0,glu1),kdia(0,2),pdia(0,2),qdia(0,glu2),
     &            psi(1,-1,2),psi(1,-1,1),epsW(0,1),epsG(0,glu1),epsV(0,lephel,qtype(2),-1,1),epsG(0,glu2),
     &            -1,musq,2,3,1,hexres(1,15),hexgauge(1,15),div,nGaugeNonAbe)!,wardidtest(15))
      wardidtest(15)=wardhNoAbe(1)
      wardidtest_Q(15)=wardhNoAbe(2)  


c W Z glu1
      call HexlineNoAbe1gTotalQ(pdia(0,1),kdia(0,1),kdia(0,2),qdia(0,glu1),pdia(0,2),qdia(0,glu2),
     &            psi(1,-1,2),psi(1,-1,1),epsW(0,1),epsV(0,lephel,qtype(2),-1,1),epsG(0,glu1),epsG(0,glu2),
     &            -1,musq,2,4,1,hexres(1,16),hexgauge(1,16),div,nGaugeNonAbe)!,wardidtest(16))
      wardidtest(16)=wardhNoAbe(1)
      wardidtest_Q(16)=wardhNoAbe(2)  


c Z glu1 W 
      call HexlineNoAbe1gTotalQ(pdia(0,1),kdia(0,2),qdia(0,glu1),kdia(0,1),pdia(0,2),qdia(0,glu2),
     &            psi(1,-1,2),psi(1,-1,1),epsV(0,lephel,qtype(1),-1,1),epsG(0,glu1),epsW(0,1),epsG(0,glu2),
     &            -1,musq,2,3,1,hexres(1,17),hexgauge(1,17),div,nGaugeNonAbe)!,wardidtest(17))
      wardidtest(17)=wardhNoAbe(1)
      wardidtest_Q(17)=wardhNoAbe(2)  


c Z W glu1
      call HexlineNoAbe1gTotalQ(pdia(0,1),kdia(0,2),kdia(0,1),qdia(0,glu1),pdia(0,2),qdia(0,glu2),
     &            psi(1,-1,2),psi(1,-1,1),epsV(0,lephel,qtype(1),-1,1),epsW(0,1),epsG(0,glu1),epsG(0,glu2),
     &            -1,musq,2,4,1,hexres(1,18),hexgauge(1,18),div,nGaugeNonAbe)!,wardidtest(18))
      wardidtest(18)=wardhNoAbe(1)
      wardidtest_Q(18)=wardhNoAbe(2)  


      wardpAbe(1)=.false.
      wardpAbe(2)=.false.
      wardpNoAbe(1)=.false.
      wardpNoAbe(2)=.false.
      wardhAbe(1)=.false.
      wardhAbe(2)=.false.
      wardhNoAbe(1)=.false.
      wardhNoAbe(2)=.false.


c 2 3glu vertices
c W V
      call NoAbe2gaHexLineTotalQ(pdia(0,1),kdia(0,1),kdia(0,2),pdia(0,2),qdia(0,glu2),qdia(0,glu1),
     &             psi(1,-1,2),psi(1,-1,1),epsW(0,1),epsV(0,lephel,qtype(2),-1,1),epsG(0,glu2),epsG(0,glu1),
     &		   -1,musq,2,0,1,hexamp(1,1),hexres(1,19),gaugereal(1,1),div,3)
      wardidtest(19)=wardhNoAbe(1)
      wardidtest_Q(19)=wardhNoAbe(2)  

c V W
      call NoAbe2gaHexLineTotalQ(pdia(0,1),kdia(0,2),kdia(0,1),pdia(0,2),qdia(0,glu2),qdia(0,glu1),
     &             psi(1,-1,2),psi(1,-1,1),epsV(0,lephel,qtype(1),-1,1),epsW(0,1),epsG(0,glu2),epsG(0,glu1),
     &		   -1,musq,2,0,1,hexamp(1,1),hexres(1,20),gaugereal(1,2),div,3)
      wardidtest(20)=wardhNoAbe(1)
      wardidtest_Q(20)=wardhNoAbe(2)  


      wardpAbe(1)=.false.
      wardpAbe(2)=.false.
      wardpNoAbe(1)=.false.
      wardpNoAbe(2)=.false.
      wardhAbe(1)=.false.
      wardhAbe(2)=.false.
      wardhNoAbe(1)=.false.
      wardhNoAbe(2)=.false.

c 4glu vertex
c W V
      call NoAbe4vertHexLineTotalQ(pdia(0,1),kdia(0,1),kdia(0,2),pdia(0,2),qdia(0,glu2),qdia(0,glu1),
     &             psi(1,-1,2),psi(1,-1,1),epsW(0,1),epsV(0,lephel,qtype(2),-1,1),epsG(0,glu2),epsG(0,glu1),
     &		   -1,musq,2,0,1,hexamp(1,1),hexres(1,21),gaugereal(1,3),div,3)
      wardidtest(21)   = wardhNoAbe(1)
      wardidtest_Q(21) = wardhNoAbe(2)

c V W
      call NoAbe4vertHexLineTotalQ(pdia(0,1),kdia(0,2),kdia(0,1),pdia(0,2),qdia(0,glu2),qdia(0,glu1),
     &             psi(1,-1,2),psi(1,-1,1),epsV(0,lephel,qtype(1),-1,1),epsW(0,1),epsG(0,glu2),epsG(0,glu1),
     &		   -1,musq,2,0,1,hexamp(1,1),hexres(1,22),gaugereal(1,4),div,3)

      wardidtest(22)   = wardhNoAbe(1)
      wardidtest_Q(22) = wardhNoAbe(2)

     
      wardidtest(0) = wardidtest(1).and.wardidtest(2).and.wardidtest(3)
     &          .and.wardidtest(4).and.wardidtest(5).and.wardidtest(6).and.
     &                wardidtest(7).and.wardidtest(8).and.wardidtest(9)
     &          .and.wardidtest(10).and.wardidtest(11).and.wardidtest(12).and.
     &                wardidtest(13).and.wardidtest(14).and.wardidtest(15)
     &          .and.wardidtest(16).and.wardidtest(17).and.wardidtest(18).and.
     &                wardidtest(19).and.wardidtest(20).and.wardidtest(21).and.wardidtest(22)
      wardidtest_q(0) = wardidtest_q(1).and.wardidtest_q(2).and.wardidtest_q(3)
     &          .and.wardidtest_q(4).and.wardidtest_q(5).and.wardidtest_q(6).and.
     &                wardidtest_q(7).and.wardidtest_q(8).and.wardidtest_q(9)
     &          .and.wardidtest_q(10).and.wardidtest_q(11).and.wardidtest_q(12).and.
     &                wardidtest_q(13).and.wardidtest_q(14).and.wardidtest_q(15)
     &          .and.wardidtest_q(16).and.wardidtest_q(17).and.wardidtest_q(18).and.
     &                wardidtest_q(19).and.wardidtest_q(20).and.wardidtest_q(21).and.wardidtest_q(22)
     
      wardidtest_Q(0) = wardidtest_Q(0).or.ldebug
!       if(.not.wardidtest_Q(0)) then
!         print*,wardidtest(:)
!         print*,wardidtest_Q(:)
!         pause
!       endif

      gaugeOK(1) = gaugeOK(1) .and. wardidtest(0)
      gaugeOK(2) = gaugeOK(2) .and. wardidtest_Q(0)

      gaugereal(1,1) = max(gaugereal(1,1),gaugereal(1,2),gaugereal(1,3)
     &          ,gaugereal(1,4),gaugereal(2,1),gaugereal(2,2)
     &          ,gaugereal(2,3),gaugereal(2,4))
      if(printGaugePrec.and.gaugereal(1,1).gt.comGaugePrec) then
        print*,'GaugePrec',gaugereal(1,1)
        comGaugePrec=gaugereal(1,1)
      endif

      do i=1,3
        hexres(i,21) = hexres(i,21)/2d0 ! we calculate this contribution twice, but only need it once
        hexres(i,22) = hexres(i,22)/2d0
      enddo
      
      if(printAll) then
	print*,'Hex',glu1
	print*,wardidtest
	print*,wardidtest_q
	print*,hexres
	print*,hexborn
      endif
      
      
      if(ldebug) print*,wardidtest
      if(ldebug) print*,wardidtest_q
      
c abelian
c sign convention of hexagons is as Peskin-S
	

        tempborn = hexborn(1)+hexborn(2)+hexborn(3)+hexborn(4)+hexborn(5)+hexborn(6)+hexborn(7)+hexborn(8)+
     &             hexborn(9)+hexborn(10)+hexborn(11)+hexborn(12)
     
        tempvirt = hexres(1,1)+hexres(1,2)+hexres(1,3)+hexres(1,4)+hexres(1,5)+hexres(1,6)+hexres(1,7)+hexres(1,8)+
     &             hexres(1,9)+hexres(1,10)+hexres(1,11)+hexres(1,12)
     
        tempvirt2 = 0d0

        tempvirt3= hexres(2,1)+hexres(2,2)+hexres(2,3)+hexres(2,4)+hexres(2,5)+hexres(2,6)+hexres(2,7)+hexres(2,8)+
     &             hexres(2,9)+hexres(2,10)+hexres(2,11)+hexres(2,12)

c nAbelian
        tempvirt3= tempvirt3+(0d0,0.25d0)*(hexres(1,13)+hexres(1,14)+hexres(1,15)+hexres(1,16)+hexres(1,17)+hexres(1,18))
        tempvirt = tempvirt +(0d0,-1.5d0)*(hexres(2,13)+hexres(2,14)+hexres(2,15)+hexres(2,16)+hexres(2,17)+hexres(2,18))
        tempvirt2= tempvirt2+(0d0,-1.5d0)*(hexres(3,13)+hexres(3,14)+hexres(3,15)+hexres(3,16)+hexres(3,17)+hexres(3,18))

c 2 3glu vertices
        tempvirt = tempvirt  -1.5d0 *(hexres(1,19)+hexres(1,20))
        tempvirt3= tempvirt3 -0.25d0*(hexres(1,19)+hexres(1,20))

c 4 glu vertex
        tempvirt = tempvirt +1.5d0*(hexres(3,21)+hexres(3,22)+hexres(2,21)+hexres(2,22)-2d0*(hexres(1,21)+hexres(1,22)))
        tempvirt2= tempvirt2+1.5d0*(hexres(1,21)+hexres(1,22)+hexres(2,21)+hexres(2,22)-2d0*(hexres(3,21)+hexres(3,22)))
        tempvirt3= tempvirt3+0.25d0*(2d0*(hexres(2,21)+hexres(2,22))-hexres(1,21)-hexres(1,22)-hexres(3,21)-hexres(3,22))
            
        tempborn = tempborn*G**2
        tempvirt = tempvirt*G**2
        tempvirt2= tempvirt2*G**2
        tempvirt3= tempvirt3*G**2
        
	
        ampb(glu1) = ampb(glu1) + tempborn
        if(div.eq.0) ampv(glu1) = ampv(glu1)-CF*tempborn ! rational term
        ampv(glu1) = ampv(glu1) + tempvirt
        ampv(glu2) = ampv(glu2) + tempvirt2
        ampv(3) = ampv(3) + tempvirt3
      if(.not.wardidtest_Q(0)) then
c	if(printGaugeFail) then
        If(printGaugeFail) print*, 'gauge test failed in qqgg, line',__LINE__,wardidtest,gaugereal
      endif


      enddo
     
      hexcount2m = hexcount2m + 1
      hexcount2sub(sub_number) = hexcount2sub(sub_number) + 1
      if(.not.gaugeOK(1)) then
        hexcountm = hexcountm + 1
        hexcountsub(sub_number) = hexcountsub(sub_number) + 1
      endif
      if(.not.gaugeOK(2)) then
        hexcountm_Q = hexcountm_q + 1
        hexcountsub_Q(sub_number) = hexcountsub_Q(sub_number) + 1
        if(return0ifGaugeFails) then
          res(0) = 0d0; res(1) = 0d0; res(2) = 0d0
          return
        else
          recalculatePSpoint=.true.
        endif
      endif
! 
      endif ! hex

c---------------------------------
c	Counter terms, see 0908.4124
c---------------------------------

      if(div.eq.1) then
	i=NFl
	ampv(1) = ampv(1) -ampb(1)*2d0*(11d0/6d0*CA-2d0/3d0*TR*NFl)
	ampv(2) = ampv(2) -ampb(2)*2d0*(11d0/6d0*CA-2d0/3d0*TR*NFl)
      endif
      if(div.eq.0.and.lsimple) then ! 
        i=NFl
	ampv(1) = ampv(1) -ampb(1)*2d0*(11d0/6d0*CA-2d0/3d0*TR*NFl)*log(musq/mursq(1,1))
	ampv(2) = ampv(2) -ampb(2)*2d0*(11d0/6d0*CA-2d0/3d0*TR*NFl)*log(musq/mursq(1,1))
      endif 


c---------------------------------
c	I Operator
c---------------------------------
c use amplitude calculated with loop routines ampb
c to split up into W, and WV contributions
c amp still contains all contributions
!       if(nlocontrib.eq.7) then
!         ampb(1) = amp(1)
!         ampb(2) = amp(2)
!         ampb(3) = amp(3)
!       endif

      if(lCSI.or.ldoblha) then
      logmuosij(1,2) = log(musq/dotrr(pbar(0,1),pbar(0,2))/2d0)
      logmuosij(1,3) = log(musq/dotrr(pbar(0,1),qbar(0,1))/2d0)
      logmuosij(1,4) = log(musq/dotrr(pbar(0,1),qbar(0,2))/2d0)
      logmuosij(2,3) = log(musq/dotrr(pbar(0,2),qbar(0,1))/2d0)
      logmuosij(2,4) = log(musq/dotrr(pbar(0,2),qbar(0,2))/2d0)
      logmuosij(3,4) = log(musq/dotrr(qbar(0,1),qbar(0,2))/2d0)
      ampb(1) = amp(1)
      ampb(2) = amp(2)
      ampb(3) = amp(3)
      endif
      
      if(lCSI) then
      if (ldebug) then
	print*,'pbar_ij'
	print*,dotrr(pbar(0,1),pbar(0,2))
	print*,dotrr(pbar(0,1),qbar(0,1))
	print*,dotrr(pbar(0,1),qbar(0,2))
	print*,dotrr(pbar(0,2),qbar(0,1))
	print*,dotrr(pbar(0,2),qbar(0,2))
	print*,dotrr(qbar(0,1),qbar(0,2))
      endif
      
      if(div.eq.2) then
	ampi(1) = ampb(1)*26d0/3d0
	ampi(2) = ampb(2)*26d0/3d0
	ampi(3) = 0d0
      elseif(div.eq.1) then
	ampi(1) = -(ampb(1)*(-45d0 + 2d0*NFL + logmuosij(1,2) - 9d0*logmuosij(1,3) - 9d0*logmuosij(2,4) - 9d0*logmuosij(3,4)))/3.
	ampi(2) = -(ampb(2)*(-45d0 + 2d0*NFL + logmuosij(1,2) - 9d0*logmuosij(1,4) - 9d0*logmuosij(2,3) - 9d0*logmuosij(3,4)))/3.
	ampi(3) = ((ampb(1) + ampb(2))*logmuosij(1,2) - ampb(2)*logmuosij(1,3) - ampb(1)*logmuosij(1,4) - ampb(1)*logmuosij(2,3)
     &              - ampb(2)*logmuosij(2,4) + ampb(1)*logmuosij(3,4) + ampb(2)*logmuosij(3,4))/2.
      else
        if(alphaDip.eq.1d0) then
	  ampi(1) = -(ampb(1)*(-840 + 32*NFL + 104*Pi**2 + 9*logmuosij(1,2) + 3*logmuosij(1,2)**2 + 3*(-30 + NFL)*logmuosij(1,3) 
     &               - 27*logmuosij(1,3)**2 - 90*logmuosij(2,4) + 3*NFL*logmuosij(2,4) - 27*logmuosij(2,4)**2 - 99*logmuosij(3,4)
     &               + 6*NFL*logmuosij(3,4) - 27*logmuosij(3,4)**2))/18.
	  ampi(2) = -(ampb(2)*(-840 + 32*NFL + 104*Pi**2 + 9*logmuosij(1,2) + 3*logmuosij(1,2)**2 + 3*(-30 + NFL)*logmuosij(1,4) -
     &                 27*logmuosij(1,4)**2 - 90*logmuosij(2,3) + 3*NFL*logmuosij(2,3) - 27*logmuosij(2,3)**2 - 99*logmuosij(3,4)
     &                  + 6*NFL*logmuosij(3,4) - 27*logmuosij(3,4)**2))/18.
	  ampi(3) = (27*(ampb(1) + ampb(2))*logmuosij(1,2) + 9*(ampb(1) + ampb(2))*logmuosij(1,2)**2 + 
     &               (-30 + NFL)*ampb(2)*logmuosij(1,3) - 9*ampb(2)*logmuosij(1,3)**2 - 30*ampb(1)*logmuosij(1,4) + 
     &               NFL*ampb(1)*logmuosij(1,4) - 9*ampb(1)*logmuosij(1,4)**2 - 30*ampb(1)*logmuosij(2,3) + 
     &               NFL*ampb(1)*logmuosij(2,3) - 9*ampb(1)*logmuosij(2,3)**2 - 30*ampb(2)*logmuosij(2,4) + 
     &               NFL*ampb(2)*logmuosij(2,4) - 9*ampb(2)*logmuosij(2,4)**2 + 33*ampb(1)*logmuosij(3,4) - 
     &               2*NFL*ampb(1)*logmuosij(3,4) + 33*ampb(2)*logmuosij(3,4) - 2*NFL*ampb(2)*logmuosij(3,4) + 
     &               9*ampb(1)*logmuosij(3,4)**2 + 9*ampb(2)*logmuosij(3,4)**2)/36.
	else
	  ampi(1) = -(ampb(1)*(-570 - 270*alphaDip + 20*NFL + 12*alphaDip*NFL + 104*Pi**2 - 6*(-45 + 2*NFL)*Log(alphaDip) + 
     &                156*Log(alphaDip)**2 + 9*logmuosij(1,2) + 3*logmuosij(1,2)**2 - 90*logmuosij(1,3) + 3*NFL*logmuosij(1,3) - 
     &                27*logmuosij(1,3)**2 - 90*logmuosij(2,4) + 3*NFL*logmuosij(2,4) - 27*logmuosij(2,4)**2 - 99*logmuosij(3,4) +
     &                6*NFL*logmuosij(3,4) - 27*logmuosij(3,4)**2))/18.
	  ampi(2) = -(ampb(2)*(-570 - 270*alphaDip + 20*NFL + 12*alphaDip*NFL + 104*Pi**2 - 6*(-45 + 2*NFL)*Log(alphaDip) + 
     &                156*Log(alphaDip)**2 + 9*logmuosij(1,2) + 3*logmuosij(1,2)**2 - 90*logmuosij(1,4) + 3*NFL*logmuosij(1,4) - 
     &                27*logmuosij(1,4)**2 - 90*logmuosij(2,3) + 3*NFL*logmuosij(2,3) - 27*logmuosij(2,3)**2 - 99*logmuosij(3,4) + 
     &                6*NFL*logmuosij(3,4) - 27*logmuosij(3,4)**2))/18.
	  ampi(3) = (27*(ampb(1) + ampb(2))*logmuosij(1,2) + 9*(ampb(1) + ampb(2))*logmuosij(1,2)**2 + 
     &                (-30 + NFL)*ampb(2)*logmuosij(1,3) - 9*ampb(2)*logmuosij(1,3)**2 - 30*ampb(1)*logmuosij(1,4) + 
     &                NFL*ampb(1)*logmuosij(1,4) - 9*ampb(1)*logmuosij(1,4)**2 - 30*ampb(1)*logmuosij(2,3) + 
     &                NFL*ampb(1)*logmuosij(2,3) - 9*ampb(1)*logmuosij(2,3)**2 - 30*ampb(2)*logmuosij(2,4) + 
     &                NFL*ampb(2)*logmuosij(2,4) - 9*ampb(2)*logmuosij(2,4)**2 + 33*ampb(1)*logmuosij(3,4) - 
     &                2*NFL*ampb(1)*logmuosij(3,4) + 33*ampb(2)*logmuosij(3,4) - 2*NFL*ampb(2)*logmuosij(3,4) + 
     &                9*ampb(1)*logmuosij(3,4)**2 + 9*ampb(2)*logmuosij(3,4)**2)/36.
	endif
      endif
      endif ! lCSI
      if (ldoblha) then
c remember for the conversions that we take alphas/(4*pi)*2*Re(conjg(ampb)*ampv) later
c convert from Gamma(1+epsilon) to 1/Gamma(1-epsilon) convention : -alphas*pi/6*(CF+CA)*Born
        ampv(1) = ampv(1) - Pi**2/3d0*(blha_CF+blha_CA)*ampb(1)
        ampv(2) = ampv(2) - Pi**2/3d0*(blha_CF+blha_CA)*ampb(2)
c add dim red contribution if blha_tgamma{Quark.Gluon} != 0: -alphas/pi*(gammatilde_q+gammatilde_g)(with incl.-)*Born
        ampv(1) = ampv(1) - 2d0*(blha_tgammaQuark+blha_tgammaGluon)*ampb(1)
        ampv(2) = ampv(2) - 2d0*(blha_tgammaQuark+blha_tgammaGluon)*ampb(2)

c set divergent epsilon^-2 and epsilon^-1 terms
        if (blha_amptype(blha_curproc).eq.1) then ! loop
            blha_amp(1)=blha_amp(1)
     &        -als(1,1)*(blha_CA+blha_CF)/Pi*
     &          dble(blha_CF**2*blha_CA*(ampb(1)*dconjg(ampb(1))
     &                              +ampb(2)*dconjg(ampb(2)))
     &              -blha_CF*ampb(1)*dconjg(ampb(2)))
            blha_amp(2)=blha_amp(2)
     &        -als(1,1)*blha_CF/blha_CA/4d0/Pi*
     &          ( ( logmuosij(1,2) - logmuosij(1,4) - logmuosij(2,3) 
     &              + blha_CA**2*(4*blha_CF*(blha_gammaGluon + blha_gammaQuark) 
     &              - logmuosij(1,3) + logmuosij(1,4) + logmuosij(2,3) 
     &              - logmuosij(2,4) + 4*CF**2*(logmuosij(1,4) 
     &              + logmuosij(2,3) - logmuosij(3,4))) + logmuosij(3,4) 
     &              + 4*blha_CA**3*blha_CF*logmuosij(3,4)
     &            )*dble(ampb(2)*dconjg(ampb(2)))
     &           +( logmuosij(1,2) - logmuosij(1,4) - logmuosij(2,3) 
     &              + 4*blha_CA**2*blha_CF*(blha_gammaGluon + blha_gammaQuark 
     &              + blha_CF*(logmuosij(1,4) + logmuosij(2,3) 
     &              - logmuosij(3,4))) + logmuosij(3,4) 
     &              + 2*blha_CA**3*blha_CF*(logmuosij(1,3) - logmuosij(1,4) 
     &              - logmuosij(2,3) + logmuosij(2,4) + 2*logmuosij(3,4))
     &            )*dble(ampb(1)*dconjg(ampb(1)))
     &           +( -2*blha_gammaGluon - 2*blha_gammaQuark 
     &              + 2*blha_CA*logmuosij(1,2) - 2*blha_CF*logmuosij(1,2) 
     &              - blha_CA*logmuosij(1,3) - blha_CA*logmuosij(1,4) 
     &              - blha_CA*logmuosij(2,3) - blha_CA*logmuosij(2,4)
     &            )*dble(ampb(2)*dconjg(ampb(1)))*2d0*blha_CA
     &          )
        endif
      endif


c square amplitudes

      if(ldebug) then
	ortb(1) = (0d0,0.5d0)*(ampb(2)-ampb(1))
	ortb(2) = 0.5d0*(ampb(1)+ampb(2))
	ortb(3) = (ampb(1)+ampb(2))/6d0
      endif

      ortv(1) = (0d0,0.5d0)*(ampv(2)-ampv(1))
      ortv(2) = 0.5d0*(ampv(1)+ampv(2))
      ortv(3) = ampv(3)+(ampv(1)+ampv(2))/6d0

      orti(1) = (0d0,0.5d0)*(ampi(2)-ampi(1))
      orti(2) = 0.5d0*(ampi(1)+ampi(2))
      orti(3) = ampi(3)+(ampi(1)+ampi(2))/6d0


      if(div.ne.0) then
	print*,id
	print*
	print*,'final'
	print*,'div',div
	print*,'amp born  ',amp
	print*,'amp bornV ',ampb
	print*,'ampi      ',ampi
	print*,'ampv      ',ampv
	print*,'ampv/ampb ',(ampv(i)/ampb(i),i=1,3)
	print*,'-ampv/ampi',(-ampv(i)/ampi(i),i=1,3)

	print*
	print*,'final'
	print*,'div',div
	print*,'ort born  ',ort
	print*,'ort bornV ',ortb
	print*,'orti      ',orti
	print*,'ortv      ',ortv
	print*,'ortv/ortb ',(ortv(i)/ortb(i),i=1,3)
	print*,'-ortv/orti',(-ortv(i)/orti(i),i=1,3)
      endif
      
      res(2) = dble(12d0*ort(1)*dconjg(ort (1))+20d0/3d0*ort(2)*dconjg(ort (2))+24d0*ort(3)*dconjg(ort (3)))
      res(1) = dble(12d0*ort(1)*dconjg(orti(1))+20d0/3d0*ort(2)*dconjg(orti(2))+24d0*ort(3)*dconjg(orti(3)))
      res(0) = dble(12d0*ort(1)*dconjg(ortv(1))+20d0/3d0*ort(2)*dconjg(ortv(2))+24d0*ort(3)*dconjg(ortv(3)))

      if(.not.lborn) res(2) = 0d0
      
      if(ldebug) then
	print*,'Re(virt*born)'
	print*,'scale   ',sqrt(musq)
	print*,'res Born',res(2)
	print*,'res I*B2',res(1)
	print*,'res V*B ',res(0)
	print*,'(I+V)*B)',res(0)+res(1)
	print*,'-V/I    ',-res(0)/res(1)
	print*,'canell  ',abs((res(1)+res(0))/res(1))
	print*,(res(0)-res(1))/res(2)
	print*,'gaugeprec. wo Hex', gaugeprecision
!  	pause
      endif
      res(1) = res(1)*als(1,1)/2d0/pi ! 2Re *als/4pi
      res(0) = res(0)*als(1,1)/2d0/pi 

      if(IsNan(res,0d0,3,0)) then
        if(ldebug.or..not.recalculatePSpoint) print*,"NaN in "__FILE__,__LINE__
        do i=0,2
          res(i) = 0d0
        enddo
      endif
      if(recalculatePSpoint) then
        recalculatePSpoint=.false.
        res=0d0
      endif

      if (ldoblha) then
        if (blha_amptype(blha_curproc).eq.1) then ! loop
          blha_amp(3)= blha_amp(3)+
     &      res(0)
          blha_amp(4)= blha_amp(4)+
     &      res(2)
        else
          call BLHA_amptypeerror(blha_amptype(blha_curproc),__FILE__,__LINE__)
        endif
      endif

      res(1) = res(0)+res(1) ! V+I
      res(0) = res(1)+res(2) ! V+I+B

      end
