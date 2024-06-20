c ************************************************************************************
c Author: Matthias Kerner
c Date: 2012/13
c ************************************************************************************

c 4q WW Amplitude
      subroutine amp4qwwVirt(P,id,fsign,nlo,res)
      implicit none
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/scales.inc"
#include "common.inc"
      real*8 P(0:3,max_p,max_kin),res(0:2,0:2)
      integer id(4),fsign(4),nlo
      
      real*8 pbar(0:3,4),pdia(0:3,4),kdia(0:3,2),momTemp(0:4)
      real*8 resB(0:2),resI(0:2),resV(0:2)
      complex*16 kdiac(0:3,2),momtempc(0:3)
      real*8 mompsiW(0:4,4,2)
      complex*16 psi(2,-1:1,4), psiW(2,4,2),psi4(4,4)
      complex*16 amp(2)
      complex*16 ampv(2),ampb(2),ampvcross1(2)
      complex*16 gluWij(0:3,1:3,2:4,2)
      complex*16 tempCurr(0:3,-1:1,4)
      complex*16 dotcc
      real*8 mass2,qjj2,dotrr
      real*8 squared(2,2), tempR, musq,logmuosij(4,4)

      
      integer i,j,mu,tempi4(4),cross,col(2),w1,w2,crosssign
      integer in,out,in2,out2

      complex*16 epsW
      real*8 momW
      common /QCDWWcurrents/ epsW(0:3,2,7),momW(0:4,2,7)
      
      real*8 deltaCVirt
      parameter (deltaCVIRT = 0d0)


      real*8 pi2,theta
      logical IsNaN
      external theta,IsNaN
      parameter (pi2 =9.86960440108935861883449099987615114d0)
      logical recalculatePSpoint
      common /comRecalcPSpoint/ recalculatePSpoint
      
      logical lhex,lbox, lborn,lcsi, wardidtest(3)
      complex*16 hexamp(9),hexres,hexborn
      complex*16 boxgauge(2,3),boxres(5,3),boxgaugeb(2,3),boxborn(3),nothing(3),checkg(5,3,2),checkgb(2,2),boxresRen(3)

      real*8 td,gaugeprecision
      real*8 CF,CA,CFCA,TR
      parameter( CF = 4d0/3d0, CA = 3d0, CFCA=-1d0/6d0, TR=0.5d0)
      integer div
      logical ldebug,shiftaround,ltop
      real*8 K_cont
      
      ltop = NFl.eq.5
      
      ldebug = .false.
      shiftaround = .false.
      K_cont = 50d0
      if(.not.shiftaround) K_cont=0d0
      
!       if(nloContrib.eq.0) then
	lborn = nlo.eq.1
	lCSI = nlo.eq.2
	lbox = nlo.eq.2
	lhex = nlo.eq.2
!       else
! 	lborn = nloContrib.eq.1
! 	lCSI  = nloContrib.eq.2
! 	lbox  = nloContrib.eq.3
! 	lhex  = nloContrib.eq.4
!       endif
      if(.not.(lborn.or.lcsi.or.lbox.or.lhex)) print*, "nloflag", nlo
!        lhex = .true.
!       lbox = .true.
!       lborn= .false.
!       lCSI = .false.
      

      div = 0


c new numbering convention
c in   q1 >>> #### >>> q2  out
c             ####          , WZ
c	      ####
c in2  q3 >>> #### >>> q4  out2

      do mu=0,3
	pbar(mu,1)=p(mu,id(1),1)
	pbar(mu,2)=p(mu,id(3),1)
	pbar(mu,3)=p(mu,id(2),1)
	pbar(mu,4)=p(mu,id(4),1)
	pdia(mu,1)=pbar(mu,1)*fsign(1)
	pdia(mu,2)=pbar(mu,2)*fsign(3)
	pdia(mu,3)=pbar(mu,3)*fsign(2)
	pdia(mu,4)=pbar(mu,4)*fsign(4)
      enddo      
      
      do mu=0,2
	resB(mu) = 0d0
	resI(mu) = 0d0
	resV(mu) = 0d0
      enddo

c get bra and ket 
      tempi4(1) = fsign(1)
      tempi4(2) = fsign(3)
      tempi4(3) = fsign(2)
      tempi4(4) = fsign(4)
      call psi0m(4,pbar,tempi4,psi)
!       pause
      
      
!       print*,psi(1:2,-1,1)
c new convention
c in   >W>>>>W>  out
c          S         , WZ
c	   S
c in2  >>>>>>>>  out2

c W emission
      do in=1,3,2
      do j=1,2
	call ket2c(psi(1,-1,in),.true.,pdia(0,in),-1,momW(0,j,1),epsW(0,j,1),psiW(1,in,j),mompsiW(0,in,j))
      enddo
      enddo

      do out=2,4,2
      do j=1,2
	call bra2c(psi(1,-1,out),.true.,pdia(0,out),-1,momW(0,j,1),epsW(0,j,1),psiW(1,out,j),mompsiW(0,out,j))
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
	  gluWij(mu,in,out,i) = -gluWij(mu,in,out,i)/qjj2(pdia(0,in),momPsiW(0,out,i))
	enddo
      enddo
      enddo
      enddo
c calc amplitudes      
      amp(1) = ( dotcc(gluWij(0,1,2,1),gluWij(0,3,4,2))*qjj2(pdia(0,1),momPsiW(0,2,1)) ! both Wij have one propagator factor included
     1	        +dotcc(gluWij(0,1,2,2),gluWij(0,3,4,1))*qjj2(pdia(0,1),momPsiW(0,2,2)))*G**2
      amp(2) =-( dotcc(gluWij(0,1,4,1),gluWij(0,3,2,2))*qjj2(pdia(0,1),momPsiW(0,4,1))
     1    	+dotcc(gluWij(0,1,4,2),gluWij(0,3,2,1))*qjj2(pdia(0,1),momPsiW(0,4,2)))*G**2

      
      squared(1,1) = dble(amp(1)*dconjg(amp(1)))
      squared(2,2) = dble(amp(2)*dconjg(amp(2)))
      squared(1,2) = dble(amp(1)*dconjg(amp(2)))
      

      if(lborn) then
	resB(1) = resB(1) + 2d0*squared(1,1)*(1d0+K_cont*als(1,1)/(4d0*Pi))
	resB(2) = resB(2) + 2d0*squared(2,2)*(1d0+K_cont*als(1,1)/(4d0*Pi))
	resB(0) = resB(0) + 2d0*(squared(1,1)+squared(2,2)-2d0/3d0*squared(1,2))*(1d0+K_cont*als(1,1)/(4d0*Pi))
      endif

!       print*, amp(1)/2d0/G**2/GWF(1)**4*((momW(4,1,1)**2-WMASS**2 + (0d0,1d0)*WWIDTH*WMASS)*(momW(4,2,1)**2-WMASS**2 +(0d0,1d0)*WWIDTH*WMASS))/momW(4,1,1)**2/momw(4,2,1)**2
!       pause
c---------------------------------
c	VIRTUALS start here
c---------------------------------


c  change momentum convention-> all momenta flowing in
      do mu=0,3
	pdia(mu,2)=-pdia(mu,2)
	pdia(mu,4)=-pdia(mu,4)
	kdia(mu,1)=-momW(mu,1,1)
	kdia(mu,2)=-momW(mu,2,1)
	kdiac(mu,1) = kdia(mu,1)
	kdiac(mu,2) = kdia(mu,2)
      enddo    
      
      ampb(1) = 0d0
      ampb(2) = 0d0
      ampv(1) = 0d0
      ampv(2) = 0d0
      
      if(ldebug) then
	print*
	print*,"diagtophys", id(1), id(3), id(2), id(4)
	print*,"fsign     ",fsign(1),fsign(3),fsign(2),fsign(4)
	print*,'pdia:'
	do i=1,4
	  print*,(pdia(mu,i),mu=0,3)
	enddo
	do i=1,2
	  print*,(kdia(mu,i),mu=0,3)
	enddo
! 	print*,"mom conservation"
! 	print*,( sum( (/ kdia(mu,1:2),pdia(mu,1:4) /) ), mu=0,3)
	
	gaugeprecision = 0d0
      endif
      
      musq = mursq(1,1) ! 80d0**2
      if(musq.le.0d0) then
	print*,"unreasonable musq"
      endif
!       musq = 80.398d0**2

c  Catani Seymour I Operator
c  with factorized  Gamma(1+ep)
c  difference red-reg =2 Born CF (Alpha/2pi)

      if(lCSI) then
	do i=1,3
	do j=i+1,4
	  logmuosij(i,j) = log(musq/dotrr(pbar(0,i),pbar(0,j))/2d0)
	enddo
	enddo
	if(div.eq.0) then
c in dim. reg.
          if(alphaDip.eq.1d0) then
	  resI(0) = + (als(1,1)*((-32*(-15 + 2*Pi**2)*(3*squared(1,1) - 2*squared(1,2) + 3*squared(2,2)))/27. + 
     -      (-((logmuosij(1,2)**2 - logmuosij(1,4)**2 - logmuosij(2,3)**2 + logmuosij(3,4)**2)*
     -            (3*squared(1,1) - 2*squared(1,2) - 21*squared(2,2))) + 
     -         2*(logmuosij(1,3)**2 - logmuosij(1,4)**2 - logmuosij(2,3)**2 + logmuosij(2,4)**2)*
     -          (3*squared(1,1) - 10*squared(1,2) + 3*squared(2,2)) + 
     -         8*(logmuosij(1,4)**2 + logmuosij(2,3)**2)*(3*squared(1,1) - 2*squared(1,2) + 3*squared(2,2)))/9. + 
     -      3*(-((logmuosij(1,2) - logmuosij(1,4) - logmuosij(2,3) + logmuosij(3,4))
     &          *(3*squared(1,1) - 2*squared(1,2) - 21*squared(2,2)))/9. + 
     -         (8*(logmuosij(1,4) + logmuosij(2,3))*(3*squared(1,1) - 2*squared(1,2) + 3*squared(2,2)))/9. + 
     -         (-logmuosij(1,3) + logmuosij(1,4) + logmuosij(2,3) - logmuosij(2,4))*
     -          ((20*squared(1,2))/9. - (2*(squared(1,1) + squared(2,2)))/3.))))/(2.*Pi)
	  resI(1) = -(als(1,1)*(-480 + 64*Pi**2 + 9*logmuosij(1,2) + 3*logmuosij(1,2)**2 
     & - 18*logmuosij(1,3) - 6*logmuosij(1,3)**2 - 63*logmuosij(1,4) - 
     -       21*logmuosij(1,4)**2 - 63*logmuosij(2,3) - 21*logmuosij(2,3)**2
     &  - 18*logmuosij(2,4) - 6*logmuosij(2,4)**2 + 9*logmuosij(3,4) + 
     -       3*logmuosij(3,4)**2)*squared(1,1))/(18.*Pi)
	  resI(2) = -(als(1,1)*(-480 + 64*Pi**2 - 63*logmuosij(1,2) - 21*logmuosij(1,2)**2
     &  - 18*logmuosij(1,3) - 6*logmuosij(1,3)**2 + 9*logmuosij(1,4) + 
     -       3*logmuosij(1,4)**2 + 9*logmuosij(2,3) + 3*logmuosij(2,3)**2 - 18*logmuosij(2,4)
     &  - 6*logmuosij(2,4)**2 - 63*logmuosij(3,4) - 
     -       21*logmuosij(3,4)**2)*squared(2,2))/(18.*Pi)

          else
            resI(0) = resI(0) +         (als(1,1)*((-16*(-21 + 4*Pi**2 - 9*alphaDip + 9*Log(alphaDip)
     &           + 6*Log(alphaDip)**2)*(3*squared(1,1) - 2*squared(1,2) + 3*squared(2,2)))/27. + 
     -      (-((logmuosij(1,2)**2 - logmuosij(1,4)**2 - logmuosij(2,3)**2 + logmuosij(3,4)**2)*
     -            (3*squared(1,1) - 2*squared(1,2) - 21*squared(2,2))) + 
     -         2*(logmuosij(1,3)**2 - logmuosij(1,4)**2 - logmuosij(2,3)**2 + logmuosij(2,4)**2)*
     -          (3*squared(1,1) - 10*squared(1,2) + 3*squared(2,2)) + 
     -         8*(logmuosij(1,4)**2 + logmuosij(2,3)**2)*(3*squared(1,1)
     &           - 2*squared(1,2) + 3*squared(2,2)))/9. + 
     -      3*(-((logmuosij(1,2) - logmuosij(1,4) - logmuosij(2,3)
     &           + logmuosij(3,4))*(3*squared(1,1) - 2*squared(1,2) - 21*squared(2,2)))/9. + 
     -         (8*(logmuosij(1,4) + logmuosij(2,3))*(3*squared(1,1) - 2*squared(1,2) + 3*squared(2,2)))/9. + 
     -         (-logmuosij(1,3) + logmuosij(1,4) + logmuosij(2,3) - logmuosij(2,4))*
     -          ((20*squared(1,2))/9. - (2*(squared(1,1) + squared(2,2)))/3.))))/(2.*Pi)
            resI(1) = resI(1) + (als(1,1)*(336 + 144*alphaDip - 64*Pi**2 - 144*Log(alphaDip)
     &           - 96*Log(alphaDip)**2 - 9*logmuosij(1,2) - 3*logmuosij(1,2)**2 + 
     -         18*logmuosij(1,3) + 6*logmuosij(1,3)**2 + 63*logmuosij(1,4) + 21*logmuosij(1,4)**2
     &           + 63*logmuosij(2,3) + 21*logmuosij(2,3)**2 + 
     -         18*logmuosij(2,4) + 6*logmuosij(2,4)**2 - 9*logmuosij(3,4) 
     &          - 3*logmuosij(3,4)**2)*squared(1,1))/(18.*Pi)
            resI(2) = resI(2) + (als(1,1)*(336 + 144*alphaDip - 64*Pi**2
     &           - 144*Log(alphaDip) - 96*Log(alphaDip)**2 + 63*logmuosij(1,2) + 21*logmuosij(1,2)**2 + 
     -         18*logmuosij(1,3) + 6*logmuosij(1,3)**2 - 9*logmuosij(1,4) -
     &           3*logmuosij(1,4)**2 - 9*logmuosij(2,3) - 3*logmuosij(2,3)**2 + 
     -         18*logmuosij(2,4) + 6*logmuosij(2,4)**2 + 63*logmuosij(3,4) 
     &          + 21*logmuosij(3,4)**2)*squared(2,2))/(18.*Pi)
          endif


	  resI(0) = resI(0)+deltaCVIRT*2d0*(squared(1,1)+squared(2,2)-2d0/3d0*squared(1,2))*als(1,1)/2d0/pi
	  resI(1) = resI(1)+deltaCVIRT*2d0*(squared(1,1))*als(1,1)/2d0/pi
	  resI(2) = resI(2)+deltaCVIRT*2d0*(squared(2,2))*als(1,1)/2d0/pi

! ! ! c factorize 1/Gamma(1-eps)  instead of Gamma(1+eps)
! ! ! 	  res(0,0) = res(0,0) + (16*als(1,1)*squared(1,1))/(3.*Pi)*pi**2/6d0
! ! ! 	  res(0,1) = res(0,1) + (16*als(1,1)*(3*squared(1,1) - 2*squared(1,2) + 3*squared(2,2)))/(9.*Pi)*pi**2/6d0

	else if(div.eq.2) then
	  resI(1) = (16*als(1,1)*squared(1,1))/(3.*Pi)
	  resI(2) = (16*als(1,1)*squared(2,2))/(3.*Pi)
	  resI(0) = (16*als(1,1)*(3*squared(1,1) - 2*squared(1,2) + 3*squared(2,2)))/(9.*Pi)
	else if(div.eq.1) then
	  resI(0) = (als(1,1)*((16*(3*squared(1,1) - 2*squared(1,2) + 3*squared(2,2)))/3. + 
     -      2*(-((logmuosij(1,2) - logmuosij(1,4) - logmuosij(2,3)
     &  + logmuosij(3,4))*(3*squared(1,1) - 2*squared(1,2) - 21*squared(2,2)))/9. + 
     -         (8*(logmuosij(1,4) + logmuosij(2,3))*(3*squared(1,1) 
     & - 2*squared(1,2) + 3*squared(2,2)))/9. + 
     -         (-logmuosij(1,3) + logmuosij(1,4) + logmuosij(2,3) - logmuosij(2,4))*
     -          ((20*squared(1,2))/9. - (2*(squared(1,1) + squared(2,2)))/3.))))/(2.*Pi)
	  resI(1) = -(als(1,1)*(-24 + logmuosij(1,2) - 2*logmuosij(1,3) 
     & - 7*logmuosij(1,4) - 7*logmuosij(2,3) - 2*logmuosij(2,4) + logmuosij(3,4))*
     -     squared(1,1))/(3.*Pi)
	  resI(2) = +(als(1,1)*(24 + 7*logmuosij(1,2) + 2*logmuosij(1,3) 
     & - logmuosij(1,4) - logmuosij(2,3) + 2*logmuosij(2,4) + 7*logmuosij(3,4))*
     -    squared(2,2))/(3.*Pi)
	endif
      endif

c      add Counter term
      if(div.eq.1.and.lbox) then
	resv(1) = resv(1) - 4d0*als(1,1)/(4d0*pi)*(11d0/6d0*CA-2d0/3d0*TR*NFl)*2d0*squared(1,1)  ! 2 *CT_Paco,because there are CT coming from upper and lower line boxes; ohter factors 2 coming from 2Re[M M] and color
	resv(2) = resv(2) - 4d0*als(1,1)/(4d0*pi)*(11d0/6d0*CA-2d0/3d0*TR*NFl)*2d0*squared(2,2)  ! 2 *CT_Paco,because there are CT coming from upper and lower line boxes; ohter factors 2 coming from 2Re[M M] and color
	resv(0) = resv(0) - 4d0*als(1,1)/(4d0*pi)*(11d0/6d0*CA
     &          -2d0/3d0*TR*NFl)*2d0*(squared(1,1)+squared(2,2)-2d0/3d0*squared(1,2))
	if(ltop) then
	  resv(1) = resv(1) - 4d0*als(1,1)/(4d0*pi)*(-2d0/3d0*TR)*2d0*squared(1,1)  ! 2 *CT_Paco,because there are CT coming from upper and lower line boxes; ohter factors 2 coming from 2Re[M M] and color
	  resv(2) = resv(2) - 4d0*als(1,1)/(4d0*pi)*(-2d0/3d0*TR)*2d0*squared(2,2)  ! 2 *CT_Paco,because there are CT coming from upper and lower line boxes; ohter factors 2 coming from 2Re[M M] and color
	  resv(0) = resv(0) - 4d0*als(1,1)/(4d0*pi)*(
     &          -2d0/3d0*TR)*2d0*(squared(1,1)+squared(2,2)-2d0/3d0*squared(1,2))
	endif
      endif
      if(div.eq.0.and.lbox.and.ltop) then !finite part coming from onshell renormalization of top quark
	resv(1) = resv(1) + als(1,1)/pi*(-(11d0/6d0*CA-2d0/3d0*TR*NFl)*log(musq/mursq(1,1))+log(musq/tmass/tmass)/3d0)*2d0*squared(1,1)
	resv(2) = resv(2) + als(1,1)/pi*(-(11d0/6d0*CA-2d0/3d0*TR*NFl)*log(musq/mursq(1,1))+log(musq/tmass/tmass)/3d0)*2d0*squared(2,2)
	resv(0) = resv(0) + als(1,1)/pi*(-(11d0/6d0*CA-2d0/3d0*TR*NFl)*log(musq/mursq(1,1))+log(musq/tmass/tmass)/3d0)
     &          *2d0*(squared(1,1)+squared(2,2)-2d0/3d0*squared(1,2))
      endif
      

      do cross =0,1  ! cross fermions
      crosssign = 1-2*cross

      col(1) = 1+cross ! change col str for fermions flipped
      col(2) = 3-col(1)
      
      in = 1
      in2 =3
      if(cross.eq.0) then
	out= 2
	out2=4
      else
	out = 4
	out2= 2
      endif      
      
!       Lm = Log(musq)
!       Ls = mjj2(pdia(0,in),pdia(0,out))
!       ths=theta(ls)
!       ls = log(abs(ls))
!       cte(1) = -3d0*(Lm - Ls + (0d0,1d0)*pi*ths) 
!      -  - (Lm - Ls + (0d0,1d0)*pi*ths)**2
! c   substracted the c_virt
!      -  -7d0 + 1d0/3d0*pi2
! 
!       Ls = mjj2(pdia(0,in2),pdia(0,out2))
!       ths=theta(ls)
!       ls = log(abs(ls))
!       cte(2) = -3d0*(Lm - Ls + (0d0,1d0)*pi*ths) 
!      -  - (Lm - Ls + (0d0,1d0)*pi*ths)**2
! c   substracted the c_virt
!      -  -7d0 + 1d0/3d0*pi2
! 
!       cte(:) = -2d0
      
      if(lbox) then

      if(ldebug) then
	print*
	print*,'box'
      endif

      do w1 = 1,2 ! W coupling to upper line
      w2 = 3-w1
	
c upper line boxes
      do mu=0,3
	momTemp(mu) = pdia(mu,in2)+pdia(mu,out2)+kdia(mu,w2)
	momtempc(mu) = momtemp(mu)
      enddo
      
      if(ldebug) print*,"upper line, w1=",w1

c abelian
c g W
      call BoxlineAbeTotal(pdia(0,in),momtemp(0),kdia(0,w1),pdia(0,out),
     &             psi(1,-1,out),psi(1,-1,in),gluWij(0,in2,out2,w2),epsW(0,w1,1),
     &             -1,musq,1,2,3,1,boxgauge(1,1),boxres(1,1),boxgaugeb(1,1),boxborn(1),div)
c for nonAbe gauge test
      call BoxlineAbeTotal(pdia(0,in),momtemp(0),kdia(0,w1),pdia(0,out),
     &             psi(1,-1,out),psi(1,-1,in),momtempc(0),epsW(0,w1,1),
     &             -1,musq,1,2,9,-1,nothing(1),checkg(1,1,1),nothing(1),checkgb(1,1),div)

	
c Wg
      call BoxlineAbeTotal(pdia(0,in),kdia(0,w1),momtemp(0),pdia(0,out),
     &             psi(1,-1,out),psi(1,-1,in),epsW(0,w1,1),gluWij(0,in2,out2,w2),
     &             -1,musq,1,3,3,1,boxgauge(1,2),boxres(1,2),boxgaugeb(1,2),boxborn(2),div)
      call BoxlineAbeTotal(pdia(0,in),kdia(0,w1),momtemp(0),pdia(0,out),
     &             psi(1,-1,out),psi(1,-1,in),epsW(0,w1,1),momtempc(0),
     &             -1,musq,1,3,9,-1,nothing(1),checkg(1,2,1),nothing(1),checkgb(2,1),div)

! c non abelian
      call BoxlineNoAbeTotal(pdia(0,in),kdia(0,w1),pdia(0,out),momtemp(0),
     &             psi(1,-1,out),psi(1,-1,in),epsW(0,w1,1),gluWij(0,in2,out2,w2),
     &             -1,musq,3,1,boxgauge(1,3),boxres(1,3),div)
      call BoxlineNoAbeTotal(pdia(0,in),kdia(0,w1),pdia(0,out),momtemp(0),
     &             psi(1,-1,out),psi(1,-1,in),epsW(0,w1,1),momtempc(0),
     &             -1,musq,9,-1,nothing,checkg(1,3,1),div)

! ! ! ! !       call BoxlineNoAbeTotal_p(pdia(0,in),kdia(0,w1),pdia(0,out),momtemp(0),
! ! ! ! !      &             psi(1,-1,out),psi(1,-1,in),epsW(0,w1,1),gluWij(0,in2,out2,w2),
! ! ! ! !      &             -1,musq,1,3,1,nothing,boxres(1,3),div)
! ! ! ! ! 
! ! ! ! !       boxres(1,3) = boxres(2,3)


c gluon self energies     
      boxresren(1) = 0d0
      boxresren(2) = 0d0
      boxresren(3) = 0d0
      call gluonSelf(mass2(momtemp(0)),musq,NFl,boxresRen(1),div) ! has to be multiplied with Mborn
      if(ltop) call gluonSelfMass(mass2(momtemp(0)),musq,tmass,boxresRen(2),div)

c gauge checks
      call box_gauge(pdia(0,in),momtemp(0),kdia(0,w1),pdia(0,out),musq,
     2    boxgauge(1,1),boxgaugeb(1,1),wardidtest(1))
      call box_gauge(pdia(0,in),kdia(0,w1),momtemp(0),pdia(0,out),musq,
     2    boxgauge(1,2),boxgaugeb(1,2),wardidtest(2))
	
c      check epsG -> pG
      if (abs(checkg(1,3,1)).lt.1d-3) then
	td=abs(-(checkg(2,2,1)+checkg(2,1,1))/2d0+checkg(1,3,1))
      else
	td=abs(-(checkg(2,2,1)+checkg(2,1,1))/2d0/checkg(1,3,1)+1d0)
      endif
      wardidtest(3) = td.lt.GaugeAccuracy
      if(ldebug) gaugeprecision = max(gaugeprecision,td)

c counter for instatbilites, increase number of boxes oinly once
      boxcount2m = boxcount2m+1
      boxcount2sub(sub_number) = boxcount2sub(sub_number) + 1

c add to ampvv
      
	ampv(col(1)) = ampv(col(1)) - crosssign*G**2* ( CF*( boxres(1,1)+boxres(1,2)) ! multiply coupling and color factors, 2nd alphas will be multiplied at the end
     1				+ CFCA*(boxres(2,1)+boxres(2,2)) + CA*boxres(1,3)
     2				+ (boxresren(1)+boxresren(2))*(boxborn(1)+boxborn(2))  )! renbox
	if(div.eq.0) ampv(col(1)) = ampv(col(1)) + CF*crosssign*G**2*(boxborn(1)+boxborn(2)) ! dim. reg.
      if(.not.(wardidtest(1).and.wardidtest(2).and.wardidtest(3).or.ldebug.or.div.ne.0)) then
        recalculatePSpoint=.true.
        boxcountm   = boxcountm   + 1
        boxcountm_q = boxcountm_q + 1
        boxcountsub  (sub_number) = boxcountsub  (sub_number) + 1
        boxcountsub_q(sub_number) = boxcountsub_q(sub_number) + 1
      endif
      ampb(col(1)) = ampb(col(1)) - (boxborn(1)+boxborn(2))*crosssign*G**2 ! this is only for testing; when calculating boxes for both quark lines -> add born twice


!       ampb(1) = (boxborn(1)+boxborn(2))*crosssign*G**2
!       ampv(1) = ampv(col(1)) + crosssign*G**2* ( CF*( boxres(1,1)+boxres(1,2)) ! multiply coupling and color factors, 2nd alphas will be multiplied at the end
!      1				+ CFCA*(boxres(2,1)+boxres(2,2)) + CA*boxres(1,3))
!       print*,'box',boxres(1:2,1)/boxborn(1)
!       print*,'box',boxres(1,3)/(boxborn(1)+boxborn(2))
!       
!       print*,ampb(1)*cte(col(1)),cte(col(1))
!       print*,ampv(1)
!       print*,ampv(1)-ampb(1)*cte(col(1))
!       print*,boxborn
!       print*,boxres(:,1)
!       print*,boxres(:,2)
!       print*,boxres(:,3)
!       pause

c boxline corrections to 2nd quark line

      do mu=0,3
	momTemp(mu) = pdia(mu,in)+pdia(mu,out)+kdia(mu,w1)
	momtempc(mu) = momtemp(mu)
      enddo
      
      if(ldebug) print*,"lower line, w1=",w1

c abelian
c g W
      call BoxlineAbeTotal(pdia(0,in2),momtemp(0),kdia(0,w2),pdia(0,out2),
     &             psi(1,-1,out2),psi(1,-1,in2),gluWij(0,in,out,w1),epsW(0,w2,1),
     &             -1,musq,1,2,3,1,boxgauge(1,1),boxres(1,1),boxgaugeb(1,1),boxborn(1),div)
c for nonAbe gauge test
      call BoxlineAbeTotal(pdia(0,in2),momtemp(0),kdia(0,w2),pdia(0,out2),
     &             psi(1,-1,out2),psi(1,-1,in2),momtempc(0),epsW(0,w2,1),
     &             -1,musq,1,2,9,-1,nothing(1),checkg(1,1,1),nothing(1),checkgb(1,1),div)
      call BoxlineAbeTotal(pdia(0,in2),momtemp(0),kdia(0,w2),pdia(0,out2),
     &             psi(1,-1,out2),psi(1,-1,in2),gluWij(0,in,out,w1),kdiac(0,w2),
     &             -1,musq,1,2,9,-1,nothing(1),checkg(1,1,2),nothing(1),checkgb(1,2),div)
	
c Wg
      call BoxlineAbeTotal(pdia(0,in2),kdia(0,w2),momtemp(0),pdia(0,out2),
     &             psi(1,-1,out2),psi(1,-1,in2),epsW(0,w2,1),gluWij(0,in,out,w1),
     &             -1,musq,1,3,3,1,boxgauge(1,2),boxres(1,2),boxgaugeb(1,2),boxborn(2),div)
      call BoxlineAbeTotal(pdia(0,in2),kdia(0,w2),momtemp(0),pdia(0,out2),
     &             psi(1,-1,out2),psi(1,-1,in2),kdiac(0,w2),gluWij(0,in,out,w1),
     &             -1,musq,1,3,9,-1,nothing(1),checkg(1,2,2),nothing(1),checkgb(2,2),div)
      call BoxlineAbeTotal(pdia(0,in2),kdia(0,w2),momtemp(0),pdia(0,out2),
     &             psi(1,-1,out2),psi(1,-1,in2),epsW(0,w2,1),momtempc(0),
     &             -1,musq,1,3,9,-1,nothing(1),checkg(1,2,1),nothing(1),checkgb(2,1),div)

c non abelian
      call BoxlineNoAbeTotal(pdia(0,in2),kdia(0,w2),pdia(0,out2),momtemp(0),
     &             psi(1,-1,out2),psi(1,-1,in2),epsW(0,w2,1),gluWij(0,in,out,w1),
     &             -1,musq,3,1,boxgauge(1,3),boxres(1,3),div)
      call BoxlineNoAbeTotal(pdia(0,in2),kdia(0,w2),pdia(0,out2),momtemp(0),
     &             psi(1,-1,out2),psi(1,-1,in2),kdiac(0,w2),gluWij(0,in,out,w1),
     &             -1,musq,9,-1,nothing,checkg(1,3,2),div)
      call BoxlineNoAbeTotal(pdia(0,in2),kdia(0,w2),pdia(0,out2),momtemp(0),
     &             psi(1,-1,out2),psi(1,-1,in2),epsW(0,w2,1),momtempc(0),
     &             -1,musq,9,-1,nothing,checkg(1,3,1),div)

c gauge checks
      call box_gauge(pdia(0,in2),momtemp(0),kdia(0,w2),pdia(0,out2),musq,
     2    boxgauge(1,1),boxgaugeb(1,1),wardidtest(1))
      call box_gauge(pdia(0,in2),kdia(0,w2),momtemp(0),pdia(0,out2),musq,
     2    boxgauge(1,2),boxgaugeb(1,2),wardidtest(2))
	
c     check epsG -> pG
      if (abs(checkg(1,3,1)).lt.1d-3) then
	td=abs(-(checkg(2,2,1)+checkg(2,1,1))/2d0+checkg(1,3,1))
      else
	td=abs(-(checkg(2,2,1)+checkg(2,1,1))/2d0/checkg(1,3,1)+1d0)
      endif
      wardidtest(3) = td.lt.GaugeAccuracy

      if(ldebug) gaugeprecision = max(gaugeprecision,td)

c add to ampvv
      
	ampv(col(1)) = ampv(col(1)) - crosssign*G**2* ( CF*( boxres(1,1)+boxres(1,2)) ! multiply coupling and color factors
     1				+ CFCA*(boxres(2,1)+boxres(2,2)) + CA*boxres(1,3))
	if(div.eq.0) ampv(col(1)) = ampv(col(1)) + CF*crosssign*G**2*(boxborn(1)+boxborn(2)) ! dim. reg.  
      if(.not.(wardidtest(1).and.wardidtest(2).and.wardidtest(3).or.ldebug.or.div.ne.0)) then
        boxcountm   = boxcountm   + 1
        boxcountm_q = boxcountm_q + 1
        boxcountsub  (sub_number) = boxcountsub  (sub_number) + 1
        boxcountsub_q(sub_number) = boxcountsub_q(sub_number) + 1
	recalculatePSpoint=.true.
      endif
      ampb(col(1)) = ampb(col(1)) - (boxborn(1)+boxborn(2))*crosssign*G**2 ! this is only for testing; when calculating boxes for both quark lines -> add born twice
      

      enddo !w1


      if(ldebug) then
	call check("identical born amp", ampb(col(1))/amp(col(1))/2) ! ampb counted twice due to boxes of both quark lines
      endif

      if(shiftaround) then
	ampv(col(1)) = ampv(col(1)) - amp(col(1))*K_cont/2d0 ! 2 compensating 2 Re
      endif
!       pause
      endif ! lbox


      if(lhex) then
      
      psi4(1,1) = psi(1,-1,1)
      psi4(2,1) = psi(2,-1,1)
      psi4(3,1) = 0d0
      psi4(4,1) = 0d0

      psi4(1,2) = 0d0
      psi4(2,2) = 0d0
      psi4(3,2) = psi(1,-1,2)
      psi4(4,2) = psi(2,-1,2)

      psi4(1,3) = psi(1,-1,3)
      psi4(2,3) = psi(2,-1,3)
      psi4(3,3) = 0d0
      psi4(4,3) = 0d0

      psi4(1,4) = 0d0
      psi4(2,4) = 0d0
      psi4(3,4) = psi(1,-1,4)
      psi4(4,4) = psi(2,-1,4)
      
      ampb(1) = 0d0
      ampb(2) = 0d0

c counter of hexagons for instatilities
      hexcount2m = hexcount2m+1
      hexcount2sub(sub_number) = hexcount2sub(sub_number)+1
      
      do w1=1,2
      w2 = 3-w1


c in  -->---$$$$$---<-- in2
c          |     |
c          |     |
c out --<---$$$$$--->-- out2

c in mu w1 nu out ; in2 mu w2 nu out2
      if(ldebug) then
	print*
	print*, 'hex Dire, w1=',w1
      endif

c do gauge test first
c w2 -> mom_w2
      call HexBox277T(pdia(0,out),pdia(0,out2),kdia(0,w2),pdia(0,in2),pdia(0,in),kdia(0,w1),
     &               psi4(1,out),psi4(1,in),psi4(1,out2),psi4(1,in2),kdiac(0,w2),epsW(0,w1,1),
     &               musq,1,0,0,hexamp,hexres,hexborn,div)
      if(ldebug) then
	print*, "Gauge test w2"
! 	print*, "gauge born", hexborn
! 	print*, "true born, but summed w=1,2", ampb(col(1))
! 	print*, "hexres    ", hexres
! 	print*, "hexamp",hexamp
	call check('r127',(hexamp(1)+hexamp(2))/(-hexamp(7)))
	call check('r348',(hexamp(3)+hexamp(4))/(-hexamp(8)))
	call check('r569',(hexamp(5)+hexamp(6))/(-hexamp(9)))
	if(ldebug) gaugeprecision = max(gaugeprecision,abs(1d0-(hexamp(1)+hexamp(2))/(-hexamp(7))))
	if(ldebug) gaugeprecision = max(gaugeprecision,abs(1d0-(hexamp(3)+hexamp(4))/(-hexamp(8))))
	if(ldebug) gaugeprecision = max(gaugeprecision,abs(1d0-(hexamp(5)+hexamp(6))/(-hexamp(9))))
	
      endif
      
      wardidtest(1) = abs((hexamp(1)+hexamp(2))/(-hexamp(7))-1d0) .lt. GaugeAccuracy
     &            .and.abs((hexamp(3)+hexamp(4))/(-hexamp(8))-1d0) .lt. GaugeAccuracy
     &            .and.abs((hexamp(5)+hexamp(6))/(-hexamp(9))-1d0) .lt. GaugeAccuracy

c w1 -> mom_w1
      call HexBox277T(pdia(0,out),pdia(0,out2),kdia(0,w2),pdia(0,in2),pdia(0,in),kdia(0,w1),
     &               psi4(1,out),psi4(1,in),psi4(1,out2),psi4(1,in2),epsW(0,w2,1),kdiac(0,w1),
     &               musq,-1,0,0,hexamp,hexres,hexborn,div)
      if(ldebug) then
	print*, "Gauge test w1"
! 	print*, "gauge born", hexborn
! 	print*, "true born, but summed w=1,2", ampb(col(1))
! 	print*, "hexres    ", hexres
! 	print*, "hexamp",hexamp
	call check('r135',(hexamp(1)+hexamp(3))/(-hexamp(5)))
	call check('r246',(hexamp(2)+hexamp(4))/(-hexamp(6)))
	call check('r789',(hexamp(7)+hexamp(8))/(-hexamp(9)))
	if(ldebug) gaugeprecision = max(gaugeprecision,abs(1d0-(hexamp(1)+hexamp(3))/(-hexamp(5))))
	if(ldebug) gaugeprecision = max(gaugeprecision,abs(1d0-(hexamp(2)+hexamp(4))/(-hexamp(6))))
	if(ldebug) gaugeprecision = max(gaugeprecision,abs(1d0-(hexamp(7)+hexamp(8))/(-hexamp(9))))
      endif

      wardidtest(2) = abs((hexamp(1)+hexamp(3))/(-hexamp(5))-1d0) .lt. GaugeAccuracy
     &            .and.abs((hexamp(2)+hexamp(4))/(-hexamp(6))-1d0) .lt. GaugeAccuracy
     &            .and.abs((hexamp(7)+hexamp(8))/(-hexamp(9))-1d0) .lt. GaugeAccuracy

c calc amplitude
      call HexBox277T(pdia(0,out),pdia(0,out2),kdia(0,w2),pdia(0,in2),pdia(0,in),kdia(0,w1),
     &               psi4(1,out),psi4(1,in),psi4(1,out2),psi4(1,in2),epsW(0,w2,1),epsW(0,w1,1),
     &               musq,-1,0,0,hexamp,hexres,hexborn,div)

c add to ampv with color color factors
      
	ampv(col(1)) = ampv(col(1)) - hexres*crosssign*G**2/6d0 ! convert to born col str
	ampv(col(2)) = ampv(col(2)) + hexres*crosssign*G**2/2d0
      if(.not.(wardidtest(1).and.wardidtest(2).or.div.ne.0.or.ldebug)) then
        hexcountm = hexcountm + 1
        hexcountsub(sub_number) = hexcountsub(sub_number) + 1
        hexcountm_Q = hexcountm_q + 1
        hexcountsub_Q(sub_number) = hexcountsub_Q(sub_number) + 1
        recalculatePSpoint=.true.
      endif
      ampb(col(1)) = ampb(col(1)) + hexborn*crosssign*G**2
!             if(.not.(wardidtest(1).and.wardidtest(2))) then
! 	      print*,wardidtest
! 	      pause
! 	    endif

c out --<---$$$$$---<-- in2
c          |     |
c          |     |
c in  -->---$$$$$--->-- out2

c in mu w1 nu out ; in2 nu w2 mu out2
      if(ldebug) then
	print*
	print*,'hexCross, w1 = ',w1
      endif

c do gauge test first
c w2 -> mom_w2
      call HexBoxCross277T(pdia(0,in),pdia(0,out2),kdia(0,w2),pdia(0,in2),pdia(0,out),kdia(0,w1),
     &               psi4(1,in),psi4(1,out),psi4(1,out2),psi4(1,in2),kdiac(0,w2),epsW(0,w1,1),
     &               musq,1,0,0,hexamp,hexres,hexborn,div)
      if(ldebug) then
	print*, "Gauge test w2"
! 	print*, "gauge born", hexborn
! 	print*, "true born, but summed w=1,2", ampb(col(1))
! 	print*, "hexres    ", hexres
! 	print*, "hexamp",hexamp
	call check('r127',(hexamp(1)+hexamp(2))/(-hexamp(7)))
	call check('r348',(hexamp(3)+hexamp(4))/(-hexamp(8)))
	call check('r569',(hexamp(5)+hexamp(6))/(-hexamp(9)))
	if(ldebug) gaugeprecision = max(gaugeprecision,abs(1d0-(hexamp(1)+hexamp(2))/(-hexamp(7))))
	if(ldebug) gaugeprecision = max(gaugeprecision,abs(1d0-(hexamp(3)+hexamp(4))/(-hexamp(8))))
	if(ldebug) gaugeprecision = max(gaugeprecision,abs(1d0-(hexamp(5)+hexamp(6))/(-hexamp(9))))
      endif
      
      wardidtest(1) = abs((hexamp(1)+hexamp(2))/(-hexamp(7))-1d0) .lt. GaugeAccuracy
     &            .and.abs((hexamp(3)+hexamp(4))/(-hexamp(8))-1d0) .lt. GaugeAccuracy
     &            .and.abs((hexamp(5)+hexamp(6))/(-hexamp(9))-1d0) .lt. GaugeAccuracy

c w1 -> mom_w1
      call HexBoxCross277T(pdia(0,in),pdia(0,out2),kdia(0,w2),pdia(0,in2),pdia(0,out),kdia(0,w1),
     &               psi4(1,in),psi4(1,out),psi4(1,out2),psi4(1,in2),epsW(0,w2,1),kdiac(0,w1),
     &               musq,-1,0,0,hexamp,hexres,hexborn,div)
      if(ldebug) then
	print*, "Gauge test w1"
! 	print*, "gauge born", hexborn
! 	print*, "true born, but summed w=1,2", ampb(col(1))
! 	print*, "hexres    ", hexres
! 	print*, "hexamp",hexamp
	call check('r135',(hexamp(1)+hexamp(3))/(-hexamp(5)))
	call check('r246',(hexamp(2)+hexamp(4))/(-hexamp(6)))
	call check('r789',(hexamp(7)+hexamp(8))/(-hexamp(9)))
	if(ldebug) gaugeprecision = max(gaugeprecision,abs(1d0-(hexamp(1)+hexamp(3))/(-hexamp(5))))
	if(ldebug) gaugeprecision = max(gaugeprecision,abs(1d0-(hexamp(2)+hexamp(4))/(-hexamp(6))))
	if(ldebug) gaugeprecision = max(gaugeprecision,abs(1d0-(hexamp(7)+hexamp(8))/(-hexamp(9))))
      endif

      wardidtest(2) = abs((hexamp(1)+hexamp(3))/(-hexamp(5))-1d0) .lt. GaugeAccuracy
     &            .and.abs((hexamp(2)+hexamp(4))/(-hexamp(6))-1d0) .lt. GaugeAccuracy
     &            .and.abs((hexamp(7)+hexamp(8))/(-hexamp(9))-1d0) .lt. GaugeAccuracy

c calc amplitude
      call HexBoxCross277T(pdia(0,in),pdia(0,out2),kdia(0,w2),pdia(0,in2),pdia(0,out),kdia(0,w1),
     &               psi4(1,in),psi4(1,out),psi4(1,out2),psi4(1,in2),epsW(0,w2,1),epsW(0,w1,1),
     &               musq,-1,0,0,hexamp,hexres,hexborn,div)


c add to ampv with color color factors
      
	ampv(col(1)) = ampv(col(1)) + hexres*crosssign*G**2*4d0/3d0 ! convert to born col str
	ampv(col(2)) = ampv(col(2)) + hexres*crosssign*G**2/2d0
      if(.not.(wardidtest(1).and.wardidtest(2).or.div.ne.0.or.ldebug)) then
        hexcountm = hexcountm + 1
        hexcountsub(sub_number) = hexcountsub(sub_number) + 1
        hexcountm_Q = hexcountm_q + 1
        hexcountsub_Q(sub_number) = hexcountsub_Q(sub_number) + 1
        recalculatePSpoint=.true.
      endif
      ampb(col(1)) = ampb(col(1)) + hexborn*crosssign*G**2


      enddo ! w1
      if(ldebug) then
	print*
	print*, "Hex", hexres
	print*, "born from HEX", ampb(col(1))/2 ! we calculated this in Dire and in Cross
	print*, "matthias born",amp(col(1))
	call check("hexborn/born",ampb(col(1))/amp(col(1))/2)
      endif

      endif ! lhex
      

      squared(1,1) = 2d0*dble(ampv(1)*dconjg(amp(1)))
      squared(2,1) = 2d0*dble(ampv(2)*dconjg(amp(1)))
      squared(2,2) = 2d0*dble(ampv(2)*dconjg(amp(2)))
      squared(1,2) = 2d0*dble(ampv(1)*dconjg(amp(2)))
      
c this is the result for different flavours |----|**2 + |--++|**2
      if(cross.eq.0) then
	tempR = (2d0*squared(1,1)-2d0/3d0*squared(2,1))*als(1,1)/(4d0* Pi) ! alphas/4pi
	resv(1) = resv(1)+tempR
      else
	tempR = (2d0*squared(2,2)-2d0/3d0*squared(1,2))*als(1,1)/(4d0* Pi) ! alphas/4pi
	resv(2) = resv(2)+tempR
      endif
      
c different flavours
c add up both crossing, before we can square the amplitude
      if(cross.eq.0) then
	ampvcross1(1) = ampv(1) ! store this contribution
	ampvcross1(2) = ampv(2)
      else
	ampv(1) = ampv(1) + ampvcross1(1)
	ampv(2) = ampv(2) + ampvcross1(2)
      
	squared(1,1) = 2d0*dble(ampv(1)*dconjg(amp(1)))
	squared(2,1) = 2d0*dble(ampv(2)*dconjg(amp(1)))
	squared(2,2) = 2d0*dble(ampv(2)*dconjg(amp(2)))
	squared(1,2) = 2d0*dble(ampv(1)*dconjg(amp(2)))
      
	tempR = (2d0*(squared(1,1)+squared(2,2))-2d0/3d0*(squared(1,2)+squared(2,1)))*als(1,1)/(4d0* Pi) ! alphas/4pi
	resv(0) = resv(0)+tempR
      endif

      ampv(1) = 0d0 ! reset amp for cross==1 
      ampv(2) = 0d0
      
      enddo ! cross

      if(div.ne.0) then
	print*,'div ',div
	print*,'born ',resB
	print*,'I    ',resI
	print*,'virt ',resV
	print*,'ratio',(resI(i)/resV(i),i=0,2)
	print*,'canel',((resI(i)+resV(i))/resV(i),i=0,2)
        print*,'V/B  ',ampv/ampb
      endif
!       print*,recalculatePSpoint
!         print*,'V/B  ',resV(1)/resB(1)/2
!         print*,'V/B  ',(resV(1)/2-16d0/3d0*resB(1)*pi**2/6)/resB(1)

      res(0,0) = resB(0)+resI(0)+resV(0)
      res(0,1) = resB(1)+resI(1)+resV(1)
      res(0,2) = resB(2)+resI(2)+resV(2) 
      
      res(1,0) = res(0,0)
      res(1,1) = res(0,1)
      res(1,2) = res(0,2)
      res(2,0) = 0d0
      res(2,1) = 0d0
      res(2,2) = 0d0
      
      
      if(IsNan(res,0d0,9,0)) then
        if(ldebug.or..not.recalculatePSpoint) print*,"NaN in "__FILE__,__LINE__
        do i=0,2; do j=0,2
          res(i,j) = 0d0
        enddo; enddo
      endif
      if(recalculatePSpoint) then
        recalculatePSpoint=.false.
        res=0d0
      endif

      if(ldebug) then
	print*,"res: ",res(0,0),res(0,1),res(0,2)
      endif
      end       

      subroutine check(text,value)
      implicit none
      character*(*) text
      complex*16 value
      
      write(*,"('check ',A,': ',T35,2(F21.16))") text,value
      if(abs(value-1d0).gt.1d-6) then
	print*,'bad precision'
! 	pause
      endif
      
      end

