c ************************************************************************************
c Author: Matthias Kerner
c Date: 2012/13
c ************************************************************************************

! #define CALC_WJJ
c 4q WZ Amplitude
      subroutine amp4qwzVirt(P,id,fsign,nbos,qhela,lephel,nlo,rn,res)

      use globalvars, only: ldoblha

      implicit none
#include "common.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"
      real*8 P(0:3,max_p,max_kin),rn,res(0:2,3:4,0:2)
      integer id(4),fsign(4),nbos,qhel,lephel,nlo,qhela
      
      real*8 pbar(0:3,4),pdia(0:3,4),kdia(0:3,3),momTemp(0:4)
      complex*16 kdiac(0:3,3)
      real*8 mompsiW(0:4,4),mompsiV(0:4,4),mompsiWV(0:4,4)
      real*8 gluijMom(0:4,1:3,2:4),gluWijMom(0:4,1:3,2:4),gluVijMom(0:4,1:3,2:4)
      complex*16 gluijMomc(0:4,1:3,2:4),gluWijMomc(0:4,1:3,2:4),gluVijMomc(0:4,1:3,2:4)
      complex*16 psi(2,-1:1,4), psiW(2,4),psi4(4,4),psiV(2,-1:1,3:4,4),psiWV(2,4)
      complex*16 psiTemp(2,-1:1),psiVV(2,4), amp(2,3:4,-1:1) ! colstr, 2nd, hel
      complex*16 ampv(2,3:4,-1:1),ampb(2,3:4,-1:1),ampvcross1(2,3:4,-1:1)
      complex*16 gluij(0:3,-1:1,1:3,2:4),gluWij(0:3,1:3,2:4),gluWVij(0:3,1:3,2:4),tempvirt
      complex*16 gluVij(0:3,-1:1,3:4,1:3,2:4),tempCurr(0:3,-1:1,4)
      complex*16 s1c, dotcc
      real*8 mass2,qjj2,dotrr
      real*8 squared(2,2), tempR, musq,logmuosij(4,4)
      real*8 resv(0:2,3:4),resI(0:2,3:4),resb(0:2,3:4)


      double precision  clr, xm2, xmg, b
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6)

      complex*16 epsW(0:3,7),epsWV(0:3,-1:1,7) 
      complex*16 epsV(0:3,-1:1,3:4,-1:1,7) ! ( , mu-hel, u/d , quark-hel), coupl. const to quarks included
      complex*16 epsA(0:3,-1:1),epsZ(0:3,-1:1)
      real*8 momV(0:4,7),momW(0:4,7),momWV(0:4,7)
      common /qcdwznew/ epsW, epsWV, epsV,epsA,epsZ,momV,momW,momWV

      
      integer i,j,k,h,q,mu,tempi4(4),cross,col(2),crosssign,colstr,helcfg,comp,ngaugeh(-4:4)
      integer minhel, maxhel, in,out,in2,out2,hel(4),qtype(3),rndhel,gtesthel,gtestq,temphel ! fermion line 1:  qtype(1) > W > qtype(2), if cross==1: only qtype(3) needs to be calculated for 2nd qline

      real*8 deltaCVirt
      parameter (deltaCVIRT = -0d0)


      logical eval(4,-1:1),return0ifGaugeFails,isnan
      external isNaN
      logical recalculatePSpoint
      common /comRecalcPSpoint/ recalculatePSpoint
      
      complex*16 vertres(2)
      real*8 vertgauge(6)
      logical lhex,lpenl,lpenb,lbox, lborn,lborn2,lcsi, wardidtest(0:8),lhel(-1:1),lpasshtest(3)
      complex*16 hexamp(9),hexres,hexborn,tempborn,gluren(3),tempHexRes(2,-1:1),tempHexResG(2,2),tempHexResB(2,-1:1)
      complex*16 penbamp(3),penbres,penbborn
      complex*16 boxgauge(2,3),boxres(5,3,3:4,-1:1),boxgaugeb(2,3),boxborn(3,3:4,-1:1),nothing(8),checkg(5,3,2),checkgb(2,2)
      complex*16 pengauge(3,8),penres(8,8,-1:1),pengaugeb(3,8),penborn(8,-1:1) 

      real*8 td,gaugeprecision
      real*8 CF,CA,CFCA,TR
      parameter( CF = 4d0/3d0, CA = 3d0, CFCA=-1d0/6d0, TR=0.5d0)
      logical ltop, gaugeOK(2)
      integer div
      logical ldebug, lW, lWV1, lWV2
      

      logical wardidtest_Q(0:8),wardpAbe(2),wardpNoAbe(2)
      common/wardpenAbe/wardpAbe
      common/wardpenNoAbe/wardpNoAbe
      

      ltop = NFL.ge.5
      ldebug = .false.
      return0ifGaugeFails = .false.
      
      do i=0,2; do j=3,4
        resB(i,j)=0d0
        resV(i,j)=0d0
        resI(i,j)=0d0
      enddo; enddo
      
      lborn= nlo.eq.1
      lCSI = nlo.eq.2
      lbox = nlo.eq.2
      lpenl= nlo.eq.2
      lpenb= nlo.eq.2
      lhex = nlo.eq.2
      lborn2=nlo.eq.2
      if (ldoblha) then
        mursq(1,1) = blha_scale**2
        lborn=  .true. 
        lCSI =  .false. 
        lbox =  .true. 
        lpenl=  .true. 
        lpenb=  .true. 
        lhex =  .true. 
        lborn2= .true. 
      endif

!       if(.not.(lborn.or.lcsi.or.lbox.or.lhex)) print*, "nloflag", nlo
!       lborn = .false.
!       lhex  = .false.
!       lpenb = .false.
!       lpenl = .true.
!       lbox  = .true.
!       lborn2= .false. ! vertline and self energies
!       lCSI  = .false.

!       lhex  = .true.
!       lpenb = .true.
!       lpenl = .true.
!       lbox  = .true.
!       lborn = .false.
!       lborn2= .true.
!       lCSI  = .true.

!       lhex  = .true.
!       lpenb = .true.
!       lpenl = .true.
!       lbox  = .true.
!       lborn = .false.
!       lborn2= .true.
!       lCSI  = .true.
      
      qhel=0
      
      
      lW = .true. ! in>W>out  in2>V>out2
      lWV1 = .true. ! in>WV>out
      lWV2 = .true. ! in>W>V>out
      
      if(.not.lWV1) then
        do mu=0,3
          epsWV(mu,lephel,1) = 0d0
        enddo
      endif

      div = 0
      
      ldebug = ldebug.or.div.ne.0
      
c select a random helicity for gauge checks
      tempR=rn ! don't change rn, so that ps point can be reevaluated
      rn = 2d0*rn
      rndhel = (-1)**INT(rn)
      rn = rn-INT(rn)
      rn = 2d0*rn ! 
      rndhel = rndhel*(3+INT(rn)) ! also rnd number specifying q type
      rn = rn-INT(rn)
      rn=tempR

      musq = mursq(1,1)
c      musq =10000000d0
      if(musq.le.0d0) then
	print*,"unreasonable musq"
      endif


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
      
      qtype(1) = nbos/10 ! (W+Z)? u type : d type
      qtype(2) = 7-qtype(1)

c get bra and ket 
      tempi4(1) = fsign(1)
      tempi4(2) = fsign(3)
      tempi4(3) = fsign(2)
      tempi4(4) = fsign(4)
      call psi0m(4,pbar,tempi4,psi)

#ifdef CALC_WJJ
      lW=.false.
      lWV1=.true.
      lWV2=.false.
      epsWV(:,lephel,:)=epsW
      momWV=momW
      momW=momW/2d0
      momV=momW
      qtype(1) = 3
      qtype(2) = 4
#endif

!       do i=-1,1
! 	CALL IXXXXX(P(0,id(1),1 ),0d0 ,i*fsign(1),fsign(1),W_MG(1,i,1   ))
! 	CALL IXXXXX(P(0,id(2),1   ),0d0 ,i*fsign(2),fsign(2),W_MG(1,i,3   ))
! 	CALL OXXXXX(P(0,id(3) ,1  ),0d0 ,i*fsign(3),fsign(3),W_MG(1,i,2   ))
! 	CALL OXXXXX(P(0,id(4)  ,1 ),0d0 ,i*fsign(4),fsign(4),W_MG(1,i,4   ))   
! 	print*,i,fsign
! 	print*,W_MG(1:4,i,1)
! 	print*,W_MG(1:4,i,2)
! 	print*,W_MG(1:4,i,3)
! 	print*,W_MG(1:4,i,4)
!       enddo
!       print*,psi(:,1,:)
!       psi(1:2,:,1)=W_MG(1:2,:,1)
!       psi(1:2,:,3)=W_MG(1:2,:,3)
!       psi(1:2,:,2)=W_MG(1:2,:,2)
!       psi(1:2,:,4)=W_MG(1:2,:,4)
!       print*,psi(:,1,:)
!       pause


      do i=1,4
      do j=-1,1,2
	eval(i,j) = .true. ! to avoid reevaluating parts multiple times
      enddo;enddo
      gluWij(0,1,2) = 0d0
      gluWij(0,1,4) = 0d0
      gluWij(0,3,2) = 0d0
      gluWij(0,3,4) = 0d0
      
      do i=-1,1,2
	amp(1,3,i) = 0d0
	amp(2,3,i) = 0d0
	amp(1,4,i) = 0d0
	amp(2,4,i) = 0d0
      enddo

      minhel = max(qhel,1)
      maxhel = qhel
      if (maxhel.eq.0) maxhel = 3

      do in = 1,3,2  ! determine q line coupling to W
      do out= 2,4,2
      do helcfg = minhel,maxhel
c new convention
c in   >W>>>>W>  out
c          S         , WZ
c	   S
c in2  >>>>>>>>  out2

      in2=4-in
      out2=6-out
      colstr= (in+out)/2
c test for reasonable combination of helcfg, in, out
      lpasshtest(helcfg) = .not.(colstr.eq.3  .or.  helcfg+colstr.eq.4)

      if(lpasshtest(helcfg)) then
	hel(in)  = -1
	hel(out) = -1
	if(helcfg.eq.1) then
	  hel(in2)  = -1
	  hel(out2) = -1
	else
	  hel(in2)  = 1
	  hel(out2) = 1
	endif
	
	if(eval(in,-1)) then ! calc |W,in>, |Z,in>, |WZ,in>
	  call ket2c(psi(1,-1,in),.true.,pdia(0,in),-1,momW(0,1),epsW(0,1),psiW(1,in),mompsiW(0,in))
	  call ket2c(psi(1,-1,in),.true.,pdia(0,in),-1,momWV(0,1),epsWV(0,lephel,1),psiWV(1,in),mompsiWV(0,in))
	  call ket2c(psi(1,-1,in),.true.,pdia(0,in),-1,momV(0,1),epsV(0,lephel,qtype(1),-1,1),psiV(1,-1,qtype(1),in),mompsiV(0,in))
	  call ket2c(psi(1,-1,in),.true.,pdia(0,in),-1,momV(0,1),epsV(0,lephel,qtype(2),-1,1),psiV(1,-1,qtype(2),in),mompsiV(0,in)) ! speed up: calc only for spec. helcfg
	  call ket2c(psiW(1,in),.false.,mompsiW(0,in),-1,momV(0,1),epsV(0,lephel,qtype(2),-1,1),psitemp(1,-1),momtemp)
	  call ket2c(psiV(1,-1,qtype(1),in),.false.,mompsiV(0,in),-1,momW(0,1),epsW(0,1),psitemp(1,1),momtemp)
	  psiVV(1,in) = psiTemp(1,-1)+psiTemp(1,1)-psiWV(1,in)
	  psiVV(2,in) = psiTemp(2,-1)+psiTemp(2,1)-psiWV(2,in)
	  eval(in,-1) = .false.
	endif
	if(eval(out,-1)) then ! calc <out,W|, <out,V|, <out,WV|
	  call bra2c(psi(1,-1,out),.true.,pdia(0,out),-1,momW(0,1),epsW(0,1),psiW(1,out),mompsiW(0,out))
	  call bra2c(psi(1,-1,out),.true.,pdia(0,out),-1,momWV(0,1),epsWV(0,lephel,1),psiWV(1,out),mompsiWV(0,out))
	  call bra2c(psi(1,-1,out),.true.,pdia(0,out),-1,momV(0,1),epsV(0,lephel,qtype(2),-1,1),psiV(1,-1,qtype(2),out),mompsiV(0,out))
	  call bra2c(psi(1,-1,out),.true.,pdia(0,out),-1,momV(0,1),epsV(0,lephel,qtype(1),-1,1),psiV(1,-1,qtype(1),out),mompsiV(0,out)) ! speed up: calc only for spec. helcfg
	  call bra2c(psiW(1,out),.false.,mompsiW(0,out),-1,momV(0,1),epsV(0,lephel,qtype(1),-1,1),psitemp(1,-1),momtemp)
	  call bra2c(psiV(1,-1,qtype(2),out),.false.,mompsiV(0,out),-1,momW(0,1),epsW(0,1),psitemp(1,1),momtemp)
	  psiVV(1,out) = psiTemp(1,-1)+psiTemp(1,1)-psiWV(1,out)
	  psiVV(2,out) = psiTemp(2,-1)+psiTemp(2,1)-psiWV(2,out)
	  eval(out,-1) = .false.
	endif
	if(eval(in2,hel(in2))) then ! calc |V,in2> (and |W,in2>, |WV,in2> if needed for another comb. of in/out/helcfg)
	  call ket2c(psi(1,hel(in2),in2),.true.,pdia(0,in2),hel(in2), momV(0,1),epsV(0,lephel,3,hel(in2),1),psiV(1,hel(in2),3,in2),mompsiV(0,in2))
	  call ket2c(psi(1,hel(in2),in2),.true.,pdia(0,in2),hel(in2), momV(0,1),epsV(0,lephel,4,hel(in2),1),psiV(1,hel(in2),4,in2),mompsiV(0,in2)) ! speed up possible if spec. helcfg
	  if(hel(in2).eq.-1) then
	    call ket2c(psi(1,-1,in2),.true.,pdia(0,in2),-1,momW(0,1),epsW(0,1),psiW(1,in2),mompsiW(0,in2))
	    call ket2c(psi(1,-1,in2),.true.,pdia(0,in2),-1,momWV(0,1),epsWV(0,lephel,1),psiWV(1,in2),mompsiWV(0,in2))
	    call ket2c(psiW(1,in2),.false.,mompsiW(0,in2),-1,momV(0,1),epsV(0,lephel,qtype(2),-1,1),psitemp(1,-1),momtemp)
	    call ket2c(psiV(1,-1,qtype(1),in2),.false.,mompsiV(0,in2),-1,momW(0,1),epsW(0,1),psitemp(1,1),momtemp)
	    psiVV(1,in2) = psiTemp(1,-1)+psiTemp(1,1)-psiWV(1,in2)
	    psiVV(2,in2) = psiTemp(2,-1)+psiTemp(2,1)-psiWV(2,in2)
	  endif
	  eval(in2,hel(in2)) = .false.
	endif
	if(eval(out2,hel(out2))) then ! calc |V,out2> (and |W,out2>, |WV,out2> if needed for another comb. of in/out/helcfg)
	  call bra2c(psi(1,hel(out2),out2),.true.,pdia(0,out2),hel(out2), momV(0,1),epsV(0,lephel,3,hel(out2),1),psiV(1,hel(out2),3,out2),mompsiV(0,out2))
	  call bra2c(psi(1,hel(out2),out2),.true.,pdia(0,out2),hel(out2), momV(0,1),epsV(0,lephel,4,hel(out2),1),psiV(1,hel(out2),4,out2),mompsiV(0,out2)) ! speed up possible if spec. helcfg
	  if(hel(out2).eq.-1) then
	    call bra2c(psi(1,-1,out2),.true.,pdia(0,out2),-1,momW(0,1),epsW(0,1),psiW(1,out2),mompsiW(0,out2))
	    call bra2c(psi(1,-1,out2),.true.,pdia(0,out2),-1,momWV(0,1),epsWV(0,lephel,1),psiWV(1,out2),mompsiWV(0,out2))
	    call bra2c(psiW(1,out2),.false.,mompsiW(0,out2),-1,momV(0,1),epsV(0,lephel,qtype(1),-1,1),psitemp(1,-1),momtemp)
	    call bra2c(psiV(1,-1,qtype(2),out2),.false.,mompsiV(0,out2),-1,momW(0,1),epsW(0,1),psitemp(1,1),momtemp)
	    psiVV(1,out2) = psiTemp(1,-1)+psiTemp(1,1)-psiWV(1,out2)
	    psiVV(2,out2) = psiTemp(2,-1)+psiTemp(2,1)-psiWV(2,out2)
	  endif
	  eval(out2,hel(out2)) = .false.
	endif
	
      endif ! passtest
      enddo ! helcfg
      
      
      if(.not.lWV2) then
        do i=1,4
          psiVV(1,i)=-psiWV(1,i)
          psiVV(2,i)=-psiWV(2,i)
        enddo
      endif
	
c calculate glu currents
      if(gluWij(0,in,out).eq.0d0) then
	call curr(-1,psi(1,-1,out),psiVV(1,in),tempCurr(0,-1,1))
	call curr(-1,psiVV(1,out),psi(1,-1,in),tempCurr(0,-1,2))
	call curr(-1,psiW(1,out),psiV(1,-1,qtype(1),in),tempCurr(0,-1,3))
	call curr(-1,psiV(1,-1,qtype(2),out),psiW(1,in),tempCurr(0,-1,4))
	if(.not.lwv2) then
          do mu=0,3; do h=-1,1,2
            tempcurr(mu,h,3) = 0d0
            tempcurr(mu,h,4) = 0d0
          enddo; enddo
        endif
	do mu=0,3
	  gluWVij(mu,in,out) = tempCurr(mu,-1,1)+tempCurr(mu,-1,2)+tempCurr(mu,-1,3)+tempCurr(mu,-1,4)
	enddo
	call curr(-1,psi(1,-1,out),psiW(1,in),tempCurr(0,-1,1))
	call curr(-1,psiW(1,out),psi(1,-1,in),tempCurr(0,-1,2))
	do mu=0,3
	  gluWij(mu,in,out) = -tempCurr(mu,-1,1)-tempCurr(mu,-1,2)
	enddo
	
c calculate lower line currents
	call curr(1,psi(1,-1,out2),psi(1,-1,in2),gluij(0,-1,in2,out2))
	
	call curr(1,psi(1,-1,out2),psiV(1,-1,3,in2),tempCurr(0,-1,1))
	call curr(1,psiV(1,-1,3,out2),psi(1,-1,in2),tempCurr(0,-1,2))
	call curr(1,psi(1,-1,out2),psiV(1,-1,4,in2),tempCurr(0,-1,3))
	call curr(1,psiV(1,-1,4,out2),psi(1,-1,in2),tempCurr(0,-1,4))
	do mu=0,3; do i=-1,1,2
	  gluVij(mu,i,3,in2,out2) = -tempCurr(mu,i,1)-tempCurr(mu,i,2)
	  gluVij(mu,i,4,in2,out2) = -tempCurr(mu,i,3)-tempCurr(mu,i,4)
	enddo;enddo
	
      endif
      
c eventually calc amplitudes
      do helcfg = minhel,maxhel
      if(lpasshtest(helcfg)) then
        hel(in) = -1
        hel(out)= -1
        if(helcfg.eq.1) then
          hel(in2) = -1
          hel(out2)= -1
	else
	  hel(in2) = 1
          hel(out2)= 1
	endif	
	tempvirt = GG(1)**2/qjj2(pdia(0,in2),pdia(0,out2))
     &		* s1c(psi(1,hel(out2),out2),gluWVij(0,in,out),.true.,hel(out2),psi(1,hel(in2),in2))
		
	if(lW) then
	  i = 1
	else
	  i=0
	endif
	
        if(colstr.eq.1.or.in.eq.3) then ! 2nd = qtype(1)
	    amp(colstr,qtype(1),hel(in2)) = tempvirt- i* GG(1)**2/qjj2(pdia(0,out),mompsiW(0,in)) * (
     &		  s1c(psi(1,hel(out2),out2),gluWij(0,in,out),.true.,hel(out2),psiV(1,hel(in2),qtype(1),in2))
     &		 +s1c(psiV(1,hel(out2),qtype(1),out2),gluWij(0,in,out),.true.,hel(out2),psi(1,hel(in2),in2)))
	endif
	

        if(colstr.eq.1.or.out.eq.4) then ! 2nd = qtype(2)
	    amp(colstr,qtype(2),hel(in2)) = tempvirt- i* GG(1)**2/qjj2(pdia(0,out),mompsiW(0,in)) * (
     &		  s1c(psi(1,hel(out2),out2),gluWij(0,in,out),.true.,hel(out2),psiV(1,hel(in2),qtype(2),in2))
     &		 +s1c(psiV(1,hel(out2),qtype(2),out2),gluWij(0,in,out),.true.,hel(out2),psi(1,hel(in2),in2)))
	endif
	
	if(in.eq.3) amp(colstr,qtype(1),hel(in2)) = -amp(colstr,qtype(1),hel(in2)) ! flipped fermions
	if(out.eq.4)amp(colstr,qtype(2),hel(in2)) = -amp(colstr,qtype(2),hel(in2))
	

      endif ! passhtest
      enddo ! helcfg
      enddo; enddo ! in and out|
      
      if(lCSI.or.ldoblha) then
	do i=1,3
	do j=i+1,4
	  logmuosij(i,j) = log(musq/dotrr(pbar(0,i),pbar(0,j))/2d0)
	enddo
	enddo
      endif
      
c square amplitudes
      do i = 3,4 ! 2nd quark
	do j=0,2
	  res(j,i,0) = 0d0
	  res(j,i,1) = 0d0
	  resv(j,i) = 0d0
	  resI(j,i) = 0d0
	enddo
	do h=-1,1,2
	  
	  squared(1,1) = dble(amp(1,i,h)*dconjg(amp(1,i,h)))
	  squared(2,2) = dble(amp(2,i,h)*dconjg(amp(2,i,h)))
	  if(h.eq.-1) then
	    squared(1,2) = dble(amp(1,i,h)*dconjg(amp(2,i,h)))
	  else
	    squared(1,2) = 0d0 ! only one colstr contributes for --++ and -+-+
	  endif
c	  |in out in2 out2|
c	  |----|**2 = T_ij M_i M_j
c	  |--++|**2 = T_11 M_1 M_1
c	  |-+-+|**2 = T_22 M_2 M_2

	  if(lborn) then
	    resB(1,i) = resB(1,i) + 2d0*squared(1,1)
	    resB(2,i) = resB(2,i) + 2d0*squared(2,2)
	    resB(0,i) = resB(0,i) + 2d0*(squared(1,1)+squared(2,2)-2d0/3d0*squared(1,2))
	  endif
	  

c  Catani Seymour I Operator
c  with factorized  Gamma(1+ep)
c  difference red-reg =2 Born CF (Alpha/2pi)

      if(lCSI) then
	if(div.eq.0) then
c in dim. reg.
	  if(alphaDip.eq.1d0) then
	    resI(0,i) = resI(0,i) + (als(1,1)*((-32*(-15 + 2*Pi**2)*(3*squared(1,1) - 2*squared(1,2) + 3*squared(2,2)))/27. + 
     -      (-((logmuosij(1,2)**2 - logmuosij(1,4)**2 - logmuosij(2,3)**2 + logmuosij(3,4)**2)*
     -            (3*squared(1,1) - 2*squared(1,2) - 21*squared(2,2))) + 
     -         2*(logmuosij(1,3)**2 - logmuosij(1,4)**2 - logmuosij(2,3)**2 + logmuosij(2,4)**2)*
     -          (3*squared(1,1) - 10*squared(1,2) + 3*squared(2,2)) + 
     -         8*(logmuosij(1,4)**2 + logmuosij(2,3)**2)*(3*squared(1,1) - 2*squared(1,2) + 3*squared(2,2)))/9. + 
     -      3*(-((logmuosij(1,2) - logmuosij(1,4) - logmuosij(2,3) + logmuosij(3,4))*(3*squared(1,1) 
     &           - 2*squared(1,2) - 21*squared(2,2)))/9. + 
     -         (8*(logmuosij(1,4) + logmuosij(2,3))*(3*squared(1,1) - 2*squared(1,2) + 3*squared(2,2)))/9. + 
     -         (-logmuosij(1,3) + logmuosij(1,4) + logmuosij(2,3) - logmuosij(2,4))*
     -          ((20*squared(1,2))/9. - (2*(squared(1,1) + squared(2,2)))/3.))))/(2.*Pi)

	  resI(1,i) = resI(1,i)-(als(1,1)*(-480 + 64*Pi**2 + 9*logmuosij(1,2) + 3*logmuosij(1,2)**2
     &           - 18*logmuosij(1,3) - 6*logmuosij(1,3)**2 - 63*logmuosij(1,4) - 
     -       21*logmuosij(1,4)**2 - 63*logmuosij(2,3) - 21*logmuosij(2,3)**2 - 18*logmuosij(2,4) 
     &          - 6*logmuosij(2,4)**2 + 9*logmuosij(3,4) + 
     -       3*logmuosij(3,4)**2)*squared(1,1))/(18.*Pi)
	  resI(2,i) = resI(2,i)-(als(1,1)*(-480 + 64*Pi**2 - 63*logmuosij(1,2) - 21*logmuosij(1,2)**2
     &           - 18*logmuosij(1,3) - 6*logmuosij(1,3)**2 + 9*logmuosij(1,4) + 
     -       3*logmuosij(1,4)**2 + 9*logmuosij(2,3) + 3*logmuosij(2,3)**2 - 18*logmuosij(2,4)
     &           - 6*logmuosij(2,4)**2 - 63*logmuosij(3,4) - 
     -       21*logmuosij(3,4)**2)*squared(2,2))/(18.*Pi)
	  else
	    resI(0,i) = resI(0,i) +         (als(1,1)*((-16*(-21 + 4*Pi**2 - 9*alphaDip + 9*Log(alphaDip)
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
	    resI(1,i) = resI(1,i) + (als(1,1)*(336 + 144*alphaDip - 64*Pi**2 - 144*Log(alphaDip)
     &           - 96*Log(alphaDip)**2 - 9*logmuosij(1,2) - 3*logmuosij(1,2)**2 + 
     -         18*logmuosij(1,3) + 6*logmuosij(1,3)**2 + 63*logmuosij(1,4) + 21*logmuosij(1,4)**2
     &           + 63*logmuosij(2,3) + 21*logmuosij(2,3)**2 + 
     -         18*logmuosij(2,4) + 6*logmuosij(2,4)**2 - 9*logmuosij(3,4) 
     &          - 3*logmuosij(3,4)**2)*squared(1,1))/(18.*Pi)
	    resI(2,i) = resI(2,i) + (als(1,1)*(336 + 144*alphaDip - 64*Pi**2
     &           - 144*Log(alphaDip) - 96*Log(alphaDip)**2 + 63*logmuosij(1,2) + 21*logmuosij(1,2)**2 + 
     -         18*logmuosij(1,3) + 6*logmuosij(1,3)**2 - 9*logmuosij(1,4) -
     &           3*logmuosij(1,4)**2 - 9*logmuosij(2,3) - 3*logmuosij(2,3)**2 + 
     -         18*logmuosij(2,4) + 6*logmuosij(2,4)**2 + 63*logmuosij(3,4) 
     &          + 21*logmuosij(3,4)**2)*squared(2,2))/(18.*Pi)
	  endif

	  resI(0,i) = resI(0,i)+deltaCVIRT*2d0*(squared(1,1)+squared(2,2)-2d0/3d0*squared(1,2))*als(1,1)/2d0/pi
	  resI(1,i) = resI(1,i)+deltaCVIRT*2d0*(squared(1,1))*als(1,1)/2d0/pi
	  resI(2,i) = resI(2,i)+deltaCVIRT*2d0*(squared(2,2))*als(1,1)/2d0/pi

! ! ! c factorize 1/Gamma(1-eps)  instead of Gamma(1+eps)
! ! ! 	  res(1,0) = res(1,0) + (16*als(1,1)*squared(1,1))/(3.*Pi)*pi**2/6d0
! ! ! 	  res(1,1) = res(1,1) + (16*als(1,1)*(3*squared(1,1) - 2*squared(1,2) + 3*squared(2,2)))/(9.*Pi)*pi**2/6d0

	else if(div.eq.2) then
	  resI(1,i) = resI(1,i) + (16*als(1,1)*squared(1,1))/(3.*Pi)
	  resI(2,i) = resI(2,i) + (16*als(1,1)*squared(2,2))/(3.*Pi)
	  resI(0,i) = resI(0,i) + (16*als(1,1)*(3*squared(1,1) - 2*squared(1,2) + 3*squared(2,2)))/(9.*Pi)
	else if(div.eq.1) then
	  resI(0,i) = resI(0,i) + (als(1,1)*((16*(3*squared(1,1) - 2*squared(1,2) + 3*squared(2,2)))/3. + 
     -      2*(-((logmuosij(1,2) - logmuosij(1,4) - logmuosij(2,3)
     &           + logmuosij(3,4))*(3*squared(1,1) - 2*squared(1,2) - 21*squared(2,2)))/9. + 
     -         (8*(logmuosij(1,4) + logmuosij(2,3))*(3*squared(1,1) - 2*squared(1,2) + 3*squared(2,2)))/9. + 
     -         (-logmuosij(1,3) + logmuosij(1,4) + logmuosij(2,3) - logmuosij(2,4))*
     -          ((20*squared(1,2))/9. - (2*(squared(1,1) + squared(2,2)))/3.))))/(2.*Pi)
	  resI(1,i) = resI(1,i) -(als(1,1)*(-24 + logmuosij(1,2) - 2*logmuosij(1,3) 
     &          - 7*logmuosij(1,4) - 7*logmuosij(2,3) - 2*logmuosij(2,4) + logmuosij(3,4))*
     -     squared(1,1))/(3.*Pi)
	  resI(2,i) = resI(2,i) +(als(1,1)*(24 + 7*logmuosij(1,2) + 2*logmuosij(1,3)
     &           - logmuosij(1,4) - logmuosij(2,3) + 2*logmuosij(2,4) + 7*logmuosij(3,4))*
     -    squared(2,2))/(3.*Pi)
	endif
      endif
      if (ldoblha) then
c convert from Gamma(1+epsilon) to 1/Gamma(1-epsilon) convention : -alphas*pi/3*CF*Born
        resV(0,i) = resV(0,i) - als(1,1)*Pi/3d0*blha_CF*
     &    (blha_CF*blha_CA/2d0*(squared(1,1)+squared(2,2))-blha_CF*squared(1,2))
        resV(1,i) = resV(1,i) - als(1,1)*Pi/3d0*blha_CF*
     &    (blha_CF*blha_CA/2d0*squared(1,1))
        resV(2,i) = resV(2,i) - als(1,1)*Pi/3d0*blha_CF*
     &    (blha_CF*blha_CA/2d0*squared(2,2))
c add dim red contribution if blha_tgammaQuark != 0: -2*alphas/pi*gammatilde(with incl.-)*Born
        resV(0,i) = resV(0,i) -2d0*als(1,1)/Pi*blha_tgammaQuark*
     &    (blha_CF*blha_CA/2d0*(squared(1,1)+squared(2,2))-blha_CF*squared(1,2))
        resV(1,i) = resV(1,i) -2d0*als(1,1)/Pi*blha_tgammaQuark*
     &    (blha_CF*blha_CA/2d0*squared(1,1))
        resV(2,i) = resV(2,i) -2d0*als(1,1)/Pi*blha_tgammaQuark*
     &    (blha_CF*blha_CA/2d0*squared(2,2))

c set divergent epsilon^-2 and epsilon^-1 terms
        if ( (blha_amptype(blha_curproc).eq.1) .and.  ! loop
     &       (i.eq.mod(blha_idsubproc(blha_cursubproc,blha_curproc),100)/10) ) then
          if (mod(blha_idsubproc(blha_cursubproc,blha_curproc),10).eq.0) then
            blha_amp(1)=blha_amp(1)
     &        -2d0*als(1,1)*blha_CF/Pi*
     &          (blha_CF*blha_CA/2d0*(squared(1,1)+squared(2,2))-blha_CF*squared(1,2))
            blha_amp(2)=blha_amp(2)
     &        -als(1,1)*blha_CF/4d0/Pi*
     &          ( ( -logmuosij(1,2) + 2*logmuosij(1,3) - logmuosij(1,4) 
     &              - logmuosij(2,3) + 2*blha_CA*(2*blha_gammaQuark 
     &              + blha_CF*(logmuosij(1,4) + logmuosij(2,3))) 
     &              + 2*logmuosij(2,4) - logmuosij(3,4)
     &            )*squared(1,1)
     &           +( 2*(logmuosij(1,3) - logmuosij(1,4) - logmuosij(2,3) 
     &              + logmuosij(2,4)) + blha_CA**2*(-logmuosij(1,2) 
     &              + logmuosij(1,4) + logmuosij(2,3) - logmuosij(3,4)) 
     &              + blha_CA*(4*blha_gammaQuark + 2*blha_CF*(2*logmuosij(1,2) 
     &              - logmuosij(1,4) - logmuosij(2,3) + 2*logmuosij(3,4)))
     &            )*squared(2,2)
     &           +( -2*(-logmuosij(1,2) + logmuosij(1,4) + logmuosij(2,3) 
     &              + 2*blha_CA**2*(logmuosij(1,3) - logmuosij(1,4) 
     &              - logmuosij(2,3) + logmuosij(2,4)) + blha_CA*(4*blha_gammaQuark 
     &              - 2*blha_CF*(logmuosij(1,3) - 2*logmuosij(1,4) 
     &              - 2*logmuosij(2,3) + logmuosij(2,4))) - logmuosij(3,4))
     &            )/blha_CA*squared(1,2)
     &          )
          else if (mod(blha_idsubproc(blha_cursubproc,blha_curproc),10).eq.1) then
            blha_amp(1)=blha_amp(1)
     &        -2d0*als(1,1)*blha_CF/Pi*
     &          blha_CF*blha_CA/2d0*(squared(1,1))
            blha_amp(2)=blha_amp(2)
     &        -als(1,1)*blha_CF/4d0/Pi*
     &          ( ( -logmuosij(1,2) + 2*logmuosij(1,3) - logmuosij(1,4) 
     &              - logmuosij(2,3) + 2*blha_CA*(2*blha_gammaQuark 
     &              + blha_CF*(logmuosij(1,4) + logmuosij(2,3))) 
     &              + 2*logmuosij(2,4) - logmuosij(3,4)
     &            )*squared(1,1)
     &          )
          else if (mod(blha_idsubproc(blha_cursubproc,blha_curproc),10).eq.2) then
            blha_amp(1)=blha_amp(1)
     &        -2d0*als(1,1)*blha_CF/Pi*
     &          blha_CF*blha_CA/2d0*(squared(2,2))
            blha_amp(2)=blha_amp(2)
     &        -als(1,1)*blha_CF/4d0/Pi*
     &          ( ( 2*(logmuosij(1,3) - logmuosij(1,4) - logmuosij(2,3) 
     &              + logmuosij(2,4)) + blha_CA**2*(-logmuosij(1,2) 
     &              + logmuosij(1,4) + logmuosij(2,3) - logmuosij(3,4)) 
     &              + blha_CA*(4*blha_gammaQuark + 2*blha_CF*(2*logmuosij(1,2) 
     &              - logmuosij(1,4) - logmuosij(2,3) + 2*logmuosij(3,4)))
     &            )*squared(2,2)
     &          )
          endif
        endif
      endif

c      add Counter term
      if(div.eq.1.and.lbox) then
	resv(1,i) = resv(1,i) - 4d0*als(1,1)/(4d0*pi)*(11d0/6d0*CA-2d0/3d0*TR*NFL)*2d0*squared(1,1)  ! 2 *CT_Paco,because there are CT coming from upper and lower line boxes; ohter factors 2 coming from 2Re[M M] and color
	resv(2,i) = resv(2,i) - 4d0*als(1,1)/(4d0*pi)*(11d0/6d0*CA-2d0/3d0*TR*NFL)*2d0*squared(2,2)  ! 2 *CT_Paco,because there are CT coming from upper and lower line boxes; ohter factors 2 coming from 2Re[M M] and color
	resv(0,i) = resv(0,i) - 4d0*als(1,1)/(4d0*pi)*(11d0/6d0*CA
     &          -2d0/3d0*TR*NFL)*2d0*(squared(1,1)+squared(2,2)-2d0/3d0*squared(1,2))
	if(ltop) then
	  resv(1,i) = resv(1,i) - 4d0*als(1,1)/(4d0*pi)*(-2d0/3d0*TR)*2d0*squared(1,1)  ! 2 *CT_Paco,because there are CT coming from upper and lower line boxes; ohter factors 2 coming from 2Re[M M] and color
	  resv(2,i) = resv(2,i) - 4d0*als(1,1)/(4d0*pi)*(-2d0/3d0*TR)*2d0*squared(2,2)  ! 2 *CT_Paco,because there are CT coming from upper and lower line boxes; ohter factors 2 coming from 2Re[M M] and color
	  resv(0,i) = resv(0,i) - 4d0*als(1,1)/(4d0*pi)*(
     &          -2d0/3d0*TR)*2d0*(squared(1,1)+squared(2,2)-2d0/3d0*squared(1,2))
	endif
      endif
      if(div.eq.0.and.lbox.and.ltop) then !finite part coming from onshell renormalization of top quark
	resv(1,i) = resv(1,i) + als(1,1)/pi*(-(11d0/6d0*CA-2d0/3d0*TR*NFL)*log(musq/mursq(1,1))+log(musq/tmass/tmass)/3d0)*2d0*squared(1,1)
	resv(2,i) = resv(2,i) + als(1,1)/pi*(-(11d0/6d0*CA-2d0/3d0*TR*NFL)*log(musq/mursq(1,1))+log(musq/tmass/tmass)/3d0)*2d0*squared(2,2)
	resv(0,i) = resv(0,i) + als(1,1)/pi*(-(11d0/6d0*CA-2d0/3d0*TR*NFL)*log(musq/mursq(1,1))+log(musq/tmass/tmass)/3d0)
     &          *2d0*(squared(1,1)+squared(2,2)-2d0/3d0*squared(1,2))
      endif

	enddo
      enddo

c      print*,'a',nlo,lborn,lcsi,lborn2,lbox,lpenl,lpenb,lhex
      
c---------------------------------
c	VIRTUALS start here
c---------------------------------

c  change momentum convention-> all momenta flowing in
      do mu=0,3
        pdia(mu,2)=-pdia(mu,2)
        pdia(mu,4)=-pdia(mu,4)
        kdia(mu,1)=-momW(mu,1)
        kdia(mu,2)=-momV(mu,1)
        kdia(mu,3)=kdia(mu,1)+kdia(mu,2)
        kdiac(mu,1) = kdia(mu,1)
        kdiac(mu,2) = kdia(mu,2)
        kdiac(mu,3) = kdia(mu,3)
      enddo    
      
      do i=3,4
      do j=-1,1,2
        ampb(1,i,j) = 0d0
        ampb(2,i,j) = 0d0
        ampv(1,i,j) = 0d0
        ampv(2,i,j) = 0d0
      enddo
      enddo
      
      do in=1,3,2
      do out=2,4,2
        do mu=0,3
          gluijMom(mu,in,out) = pdia(mu,in)+pdia(mu,out)
          gluWijMom(mu,in,out) = gluijMom(mu,in,out) + kdia(mu,1)
          gluVijMom(mu,in,out) = gluijMom(mu,in,out) + kdia(mu,2)
          gluijMomc(mu,in,out) = gluijMom(mu,in,out)
          gluWijMomc(mu,in,out) = gluWijMom(mu,in,out)
          gluVijMomc(mu,in,out) = gluVijMom(mu,in,out)
        enddo
        gluijMom(4,in,out) = mass2(gluijMom(0,in,out))
        gluVijMom(4,in,out) = mass2(gluVijMom(0,in,out))
        gluWijMom(4,in,out) = mass2(gluWijMom(0,in,out))
        gluijMomc(4,in,out) = gluijMom(4,in,out)
        gluWijMomc(4,in,out) = gluWijMom(4,in,out)
        gluVijMomc(4,in,out) = gluVijMom(4,in,out)
      enddo; enddo
      
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
        print*,"mom conservation"
!         print*,( sum( (/ kdia(mu,1:2),pdia(mu,1:4) /) ), mu=0,3)
        
        gaugeprecision = 0d0
      endif
    
      do in = 1,3,2  ! determine q line coupling to W
      do out= 2,4,2
      
      
      cross = (in+out)/2-1
      if (cross.gt.1) goto 40
c determine which helicity of the 2nd qline has to be calculated
      if(cross.eq.0) then
        lhel(-1) = qhel.eq.0.or.qhel.eq.1
        lhel( 1) = qhel.eq.0.or.qhel.eq.2
      else 
        lhel(-1) = qhel.eq.0.or.qhel.eq.1
        lhel( 1) = qhel.eq.0.or.qhel.eq.3
        if(in.eq.3) then ! incoming quarks have to be the same > fix 2nd quark
          qtype(3) = qtype(1)
        else
          qtype(3) = qtype(2)
        endif
      endif
      
      in2 = 4-in
      out2= 6-out
      crosssign= 1-2*cross
      helcfg = cross + 2 ! helcfg for positive helicity
      colstr = cross+1
      col(1) = colstr
      col(2) = 3-colstr

      gtesthel = sign(1,rndhel) ! rnd. helicity for gauge tests * qtype for gtest
      gtestq = abs(rndhel)
      if(.not.lhel(gtesthel)) gtesthel = -gtesthel
      if(cross.eq.1.and.qtype(3).ne.abs(gtestq)) gtestq = 7-gtestq

      
      if(ldebug) print*,'gtest hel,q',gtesthel,gtestq

c ------------------------------------------
c calculate vertline and gluon self energies     
c ------------------------------------------
      if(lborn2) then

      gaugeOK(1) = .true.
      do h =-1,1,2
      if(lhel(h)) then

c self energies for in>W>out in2>V>out2
      call gluonSelf(gluVijMom(4,in2,out2),musq,NFL,gluren(1),div)
      if(ltop) then
        call gluonSelfMass(gluVijMom(4,in2,out2),musq,tmass,gluren(2),div)
      else
        gluren(2) = 0d0
      endif
      do q=3,4
        if (cross.eq.0.or.qtype(3).eq.q) then
          tempborn = dotcc(gluWij(0,in,out),gluVij(0,h,q,in2,out2))/gluVijMom(4,in2,out2)*G**2*crosssign
          
          ampb(colstr,q,h) = ampb(colstr,q,h) + tempborn
          ampv(colstr,q,h) = ampv(colstr,q,h) + tempborn*(gluren(1)+gluren(2))
        endif
      enddo
	  
c self energies for in>WV>out in2>out2
      tempborn = dotcc(gluWVij(0,in,out),gluij(0,h,in2,out2))/gluijMom(4,in2,out2)*G**2*crosssign
      call gluonSelf(gluijMom(4,in2,out2),musq,NFL,gluren(1),div)
      if(ltop) then 
        call gluonSelfMass(gluijMom(4,in2,out2),musq,tmass,gluren(2),div)
      else
        gluren(2) = 0d0
      endif
c vertline
      call vertline(gluijMom(4,in2,out2),musq,1,gluren(3),div)

      do q=3,4
        if (cross.eq.0.or.qtype(3).eq.q) then
          ampb(colstr,q,h) = ampb(colstr,q,h) + tempborn*2d0 ! vertbox and self energy
          ampv(colstr,q,h) = ampv(colstr,q,h) + tempborn*(gluren(1)+gluren(2)+gluren(3))
        endif
      enddo

c Fermion loop triangles
      tempvirt = 0d0
c massless up type
c u and d contributions cancel for same mass quarks
c so eather calc 3*d + 2*u +1*t, or just calc 1*d + 1*t
c 	  call FermionVertCouplingT(1d0,1d0,clr(3,2,+1),1d0,1d0,clr(3,2,-1),
c      $		0d0,
c      $		gluWijMom(0,in,out),gluijMom(0,in2,out2),kdia(0,2),
c      $		gluWij(0,in,out),gluij(0,h,in2,out2),epsZ(0,lephel),
c      $		musq,1,1,vertres(1),vertgauge(1),div)
c 	  call FermionVertCouplingT(1d0,1d0,clr(3,2,+1),1d0,1d0,clr(3,2,-1),
c      $		0d0,
c      $		gluijMom(0,in2,out2),gluWijMom(0,in,out),kdia(0,2),
c      $		gluij(0,h,in2,out2),gluWij(0,in,out),epsZ(0,lephel),
c      $		musq,1,1,vertres(2),vertgauge(2),div)
c 	  tempvirt = tempvirt + 1d0*(vertres(1)+vertres(2))
c 	  print*,vertres(1:2),vertres(1)+vertres(2)


      if(NFL.eq.5) then
c massless bottom quark
        call FermionVertCouplingT(1d0,1d0,clr(4,2,+1),1d0,1d0,clr(4,2,-1),
     $		0d0,
     $		gluWijMom(0,in,out),gluijMom(0,in2,out2),kdia(0,2),
     $		gluWij(0,in,out),gluij(0,h,in2,out2),epsZ(0,lephel),
     $		musq,1,1,vertres(1),vertgauge(1),div)
        call FermionVertCouplingT(1d0,1d0,clr(4,2,+1),1d0,1d0,clr(4,2,-1),
     $		0d0,
     $		gluijMom(0,in2,out2),gluWijMom(0,in,out),kdia(0,2),
     $		gluij(0,h,in2,out2),gluWij(0,in,out),epsZ(0,lephel),
     $		musq,1,1,vertres(2),vertgauge(2),div)
        tempvirt = tempvirt + 1d0*(vertres(1)+vertres(2))

c top quark
        call FermionVertCouplingT(1d0,1d0,clr(3,2,+1),1d0,1d0,clr(3,2,-1),
     $		tmass,
     $		gluWijMom(0,in,out),gluijMom(0,in2,out2),kdia(0,2),
     $		gluWij(0,in,out),gluij(0,h,in2,out2),epsZ(0,lephel),
     $		musq,1,1,vertres(1),vertgauge(1),div)
        call FermionVertCouplingT(1d0,1d0,clr(3,2,+1),1d0,1d0,clr(3,2,-1),
     $		tmass,
     $		gluijMom(0,in2,out2),gluWijMom(0,in,out),kdia(0,2),
     $		gluij(0,h,in2,out2),gluWij(0,in,out),epsZ(0,lephel),
     $		musq,1,1,vertres(2),vertgauge(2),div)
        tempvirt = tempvirt + vertres(1)+vertres(2)
      else
        vertgauge(5) = 1d-19
        vertgauge(6) = 1d-19
      endif
	  
c 	  print*,tempvirt
c 	  call FermionTriangleGGZ(gluWijMom(0,in,out),gluijMom(0,in2,out2),kdia(0,2),
c      $		gluWij(0,in,out),gluij(0,h,in2,out2),epsZ(0,lephel),musq,tmass,tempborn)
c 	  print*,tempborn*2d0,tempborn*2d0/tempvirt
c 	  pause

      vertgauge(1) = max(vertgauge(1),vertgauge(2),vertgauge(3),vertgauge(4),vertgauge(5),vertgauge(6))
      do q=3,4
        if (cross.eq.0.or.qtype(3).eq.q) then
          ampv(colstr,q,h) = ampv(colstr,q,h) - tempvirt/gluWijMom(4,in,out)/gluijMom(4,in2,out2)*G**2*crosssign*TR
          ! see note in fermionVertCouplTemp.F for - sign
          endif
      enddo
      if(vertgauge(1).Gt.GaugeAccuracy) then
        gaugeOK(1) = .false.
      endif

      endif ! hel
      enddo ! hel
      
      box2count2m = box2count2m + 1
      box2count2sub(sub_number) = box2count2sub(sub_number) + 1
      if(.not.gaugeOK(1)) then
        box2countm   = box2countm   + 1
        box2countm_Q = box2countm_Q + 1
        box2countsub  (sub_number) = box2countsub  (sub_number) + 1
        box2countsub_q(sub_number) = box2countsub_q(sub_number) + 1
        if (return0ifGaugeFails) then
          do i=0,2; do j=3,4; do k=0,2
            res(i,j,k) = 0d0
          enddo; enddo; enddo
          return
        else
          recalculatePSpoint=.true.
        endif
      endif
      
      endif ! lborn
      
c ---------------------------------------
c $$$$$$$$$$ B O X E S $$$$$$$$$$$$$$$$$$
c ---------------------------------------
      
      if(lbox) then
        gaugeOK(1) = .true.

      if(ldebug) then
        print*
        print*,'box'
      endif

      do i=-4,4
        ngaugeh(i)=9
      enddo
      ngaugeh(gtesthel) = 3
      ngaugeh(gtesthel*gtestq) = 3

      do i=1,3
        boxborn(i,3,-1) = 0d0
        boxborn(i,4,-1) = 0d0
        boxborn(i,3, 1) = 0d0
        boxborn(i,4, 1) = 0d0
      enddo

c upper line boxes

c in> G> WV> out, in2>out2 
!       print*,'box 1a'
      comp=1
      if(lhel(-1)) then
        call boxlineABETotal(pdia(0,in),gluijMom(0,in2,out2),kdia(0,3),pdia(0,out),
     &            psi(1,-1,out),psi(1,-1,in),gluij(0,-1,in2,out2),epsWV(0,lephel,1),
     &            -1,musq,1,2,ngaugeh(-1),comp,boxgauge(1,1),boxres(1,1,3,-1),boxgaugeb(1,1),boxborn(1,3,-1),div)
        comp =-1
      endif
      
      if(lhel(1)) then
        call boxlineABETotal(pdia(0,in),gluijMom(0,in2,out2),kdia(0,3),pdia(0,out),
     &            psi(1,-1,out),psi(1,-1,in),gluij(0, 1,in2,out2),epsWV(0,lephel,1),
     &            -1,musq,1,2,ngaugeh(1),comp,boxgauge(1,1),boxres(1,1,3,1),boxgaugeb(1,1),boxborn(1,3,1),div)
        comp =-1
      endif

c	for gauge test of nonAbe boxes
      call boxlineABETotal(pdia(0,in),gluijMom(0,in2,out2),kdia(0,3),pdia(0,out),
     &            psi(1,-1,out),psi(1,-1,in),gluijMomc(0,in2,out2),epsWV(0,lephel,1),
     &            -1,musq,1,2,9,-1,nothing(1),checkg(1,1,1),nothing(1),checkgb(1,1),div)
      

c in> WV> G> out, in2>out2 
!       print*,'box 1b'
      comp=1
      if(lhel(-1)) then
        call boxlineABETotal(pdia(0,in),kdia(0,3),gluijMom(0,in2,out2),pdia(0,out),
     &            psi(1,-1,out),psi(1,-1,in),epsWV(0,lephel,1),gluij(0,-1,in2,out2),
     &            -1,musq,1,3,ngaugeh(-1),comp,boxgauge(1,2),boxres(1,2,3,-1),boxgaugeb(1,2),boxborn(2,3,-1),div)
        comp =-1
      endif
      
      if(lhel(1)) then
        call boxlineABETotal(pdia(0,in),kdia(0,3),gluijMom(0,in2,out2),pdia(0,out),
     &            psi(1,-1,out),psi(1,-1,in),epsWV(0,lephel,1),gluij(0, 1,in2,out2),
     &            -1,musq,1,3,ngaugeh(1),comp,boxgauge(1,2),boxres(1,2,3,1),boxgaugeb(1,2),boxborn(2,3,1),div)
        comp =-1
      endif

c	for gauge test of nonAbe boxes
      call boxlineABETotal(pdia(0,in),kdia(0,3),gluijMom(0,in2,out2),pdia(0,out),
     &            psi(1,-1,out),psi(1,-1,in),epsWV(0,lephel,1),gluijMomc(0,in2,out2),
     &            -1,musq,1,3,9,-1,nothing(1),checkg(1,2,1),nothing(1),checkgb(2,1),div)

c non abelian
      comp=1
!       print*,'box 1c'
      if(lhel(-1)) then
        call BoxlineNoAbeTotal(pdia(0,in),kdia(0,3),pdia(0,out),gluijMom(0,in2,out2),
     &             psi(1,-1,out),psi(1,-1,in),epsWV(0,lephel,1),gluij(0,-1,in2,out2),
     &             -1,musq,ngaugeh(-1),comp,boxgauge(1,3),boxres(1,3,3,-1),div)
        comp = -1
      endif

      if(lhel(1)) then
        call BoxlineNoAbeTotal(pdia(0,in),kdia(0,3),pdia(0,out),gluijMom(0,in2,out2),
     &             psi(1,-1,out),psi(1,-1,in),epsWV(0,lephel,1),gluij(0,1,in2,out2),
     &             -1,musq,ngaugeh(1),comp,boxgauge(1,3),boxres(1,3,3,1),div)
        comp = -1
      endif

      call BoxlineNoAbeTotal(pdia(0,in),kdia(0,3),pdia(0,out),gluijMom(0,in2,out2),
     &             psi(1,-1,out),psi(1,-1,in),epsWV(0,lephel,1),gluijMomc(0,in2,out2),
     &             -1,musq,9,-1,nothing,checkg(1,3,1),div)

c gauge tests 
      call box_gauge(pdia(0,in),gluijMom(0,in2,out2),kdia(0,3),pdia(0,out),musq,
     2    boxgauge(1,1),boxgaugeb(1,1),wardidtest(1))
      call box_gauge(pdia(0,in),kdia(0,3),gluijMom(0,in2,out2),pdia(0,out),musq,
     2    boxgauge(1,2),boxgaugeb(1,2),wardidtest(2))
!       print*,abs(boxgaugeb(:,1:2))

c      check epsG -> pG
      if (abs(checkg(1,3,1)).lt.1d-6) then
        td=abs(-(checkg(2,2,1)+checkg(2,1,1))/2d0+checkg(1,3,1))
      else
        td=abs(-(checkg(2,2,1)+checkg(2,1,1))/2d0/checkg(1,3,1)+1d0)
      endif
      wardidtest(3) = td.lt.GaugeAccuracy

      if(ldebug) print*,wardidtest(1),wardidtest(2),wardidtest(3),td,checkg(1,3,1)
!       pause

      if(ldebug) gaugeprecision = max(gaugeprecision,td)
      
c sum
      wardidtest(0) = wardidtest(1).and.wardidtest(2).and.wardidtest(3).or.ldebug
!       print*,td,abs(checkg(1,3,1))
!       print*,wardidtest(1:3)
!       if(.not.wardidtest(0)) pause
      gaugeOK(1) = gaugeOK(1).and.wardidtest(0)

      do h=-1,1,2
        if(lhel(h)) then
          tempborn = boxborn(1,3,h)+boxborn(2,3,h)
          tempvirt = CF*(boxres(1,1,3,h)+boxres(1,2,3,h))+CFCA*(boxres(2,1,3,h)+boxres(2,2,3,h))+CA*boxres(1,3,3,h)
          tempborn = tempborn*crosssign*G**2/gluijMom(4,in2,out2)
          tempvirt = tempvirt*crosssign*G**2/gluijMom(4,in2,out2)	  
          
          if(cross.eq.0.or.qtype(3).eq.3) then
            ampb(colstr,3,h) = ampb(colstr,3,h) - tempborn
            ampv(colstr,3,h) = ampv(colstr,3,h) - tempvirt
            if(div.eq.0) ampv(colstr,3,h) = ampv(colstr,3,h) + CF*tempborn ! dim. reg.  
          endif
          if(cross.eq.0.or.qtype(3).eq.4) then
            ampb(colstr,4,h) = ampb(colstr,4,h) - tempborn
            ampv(colstr,4,h) = ampv(colstr,4,h) - tempvirt
            if(div.eq.0) ampv(colstr,4,h) = ampv(colstr,4,h) + CF*tempborn ! dim. reg.  
          endif
        endif
      enddo


c in> G> W> out, in2>V>out2 
!       print*,'box 2a'
      comp=1
      if(lhel(-1).and.(cross.eq.0.or.qtype(3).eq.3)) then
        call boxlineABETotal(pdia(0,in),gluVijMom(0,in2,out2),kdia(0,1),pdia(0,out),
     &            psi(1,-1,out),psi(1,-1,in),gluVij(0,-1,3,in2,out2),epsW(0,1),
     &            -1,musq,1,2,ngaugeh(-3),comp,boxgauge(1,1),boxres(1,1,3,-1),boxgaugeb(1,1),boxborn(1,3,-1),div)
        comp =-1
      endif
      
      if(lhel(-1).and.(cross.eq.0.or.qtype(3).eq.4)) then
        call boxlineABETotal(pdia(0,in),gluVijMom(0,in2,out2),kdia(0,1),pdia(0,out),
     &            psi(1,-1,out),psi(1,-1,in),gluVij(0,-1,4,in2,out2),epsW(0,1),
     &            -1,musq,1,2,ngaugeh(-4),comp,boxgauge(1,1),boxres(1,1,4,-1),boxgaugeb(1,1),boxborn(1,4,-1),div)
        comp =-1
      endif

      if(lhel(1).and.(cross.eq.0.or.qtype(3).eq.3)) then
        call boxlineABETotal(pdia(0,in),gluVijMom(0,in2,out2),kdia(0,1),pdia(0,out),
     &            psi(1,-1,out),psi(1,-1,in),gluVij(0,1,3,in2,out2),epsW(0,1),
     &            -1,musq,1,2,ngaugeh(3),comp,boxgauge(1,1),boxres(1,1,3,1),boxgaugeb(1,1),boxborn(1,3,1),div)
        comp =-1
      endif
      
      if(lhel(1).and.(cross.eq.0.or.qtype(3).eq.4)) then
        call boxlineABETotal(pdia(0,in),gluVijMom(0,in2,out2),kdia(0,1),pdia(0,out),
     &            psi(1,-1,out),psi(1,-1,in),gluVij(0,1,4,in2,out2),epsW(0,1),
     &            -1,musq,1,2,ngaugeh(4),comp,boxgauge(1,1),boxres(1,1,4,1),boxgaugeb(1,1),boxborn(1,4,1),div)
        comp =-1
      endif

c        gauge test of nAbe boxes
      call boxlineABETotal(pdia(0,in),gluVijMom(0,in2,out2),kdia(0,1),pdia(0,out),
     &            psi(1,-1,out),psi(1,-1,in),gluVijMomc(0,in2,out2),epsW(0,1),
     &            -1,musq,1,2,9,-1,nothing,checkg(1,1,1),nothing,checkgb(1,1),div)
	
      
c in> W> G> out, in2>V>out2 
!       print*,'box 2b'
      comp=1
      if(lhel(-1).and.(cross.eq.0.or.qtype(3).eq.3)) then
        call boxlineABETotal(pdia(0,in),kdia(0,1),gluVijMom(0,in2,out2),pdia(0,out),
     &            psi(1,-1,out),psi(1,-1,in),epsW(0,1),gluVij(0,-1,3,in2,out2),
     &            -1,musq,1,3,ngaugeh(-3),comp,boxgauge(1,2),boxres(1,2,3,-1),boxgaugeb(1,2),boxborn(2,3,-1),div)
        comp =-1
      endif

      if(lhel(-1).and.(cross.eq.0.or.qtype(3).eq.4)) then
        call boxlineABETotal(pdia(0,in),kdia(0,1),gluVijMom(0,in2,out2),pdia(0,out),
     &            psi(1,-1,out),psi(1,-1,in),epsW(0,1),gluVij(0,-1,4,in2,out2),
     &            -1,musq,1,3,ngaugeh(-4),comp,boxgauge(1,2),boxres(1,2,4,-1),boxgaugeb(1,2),boxborn(2,4,-1),div)
        comp =-1
      endif

      if(lhel(1).and.(cross.eq.0.or.qtype(3).eq.3)) then
        call boxlineABETotal(pdia(0,in),kdia(0,1),gluVijMom(0,in2,out2),pdia(0,out),
     &            psi(1,-1,out),psi(1,-1,in),epsW(0,1),gluVij(0,1,3,in2,out2),
     &            -1,musq,1,3,ngaugeh(3),comp,boxgauge(1,2),boxres(1,2,3,1),boxgaugeb(1,2),boxborn(2,3,1),div)
        comp =-1
      endif

      if(lhel(1).and.(cross.eq.0.or.qtype(3).eq.4)) then
        call boxlineABETotal(pdia(0,in),kdia(0,1),gluVijMom(0,in2,out2),pdia(0,out),
     &            psi(1,-1,out),psi(1,-1,in),epsW(0,1),gluVij(0,1,4,in2,out2),
     &            -1,musq,1,3,ngaugeh(4),comp,boxgauge(1,2),boxres(1,2,4,1),boxgaugeb(1,2),boxborn(2,4,1),div)
        comp =-1	
      endif

c        gauge test of nAbe boxes
      call boxlineABETotal(pdia(0,in),kdia(0,1),gluVijMom(0,in2,out2),pdia(0,out),
     &            psi(1,-1,out),psi(1,-1,in),epsW(0,1),gluVijMomc(0,in2,out2),
     &            -1,musq,1,3,9,-1,nothing,checkg(1,2,1),nothing,checkgb(2,1),div)
	

c non abelian
!       print*,'box 2c'
      comp=1
      if(lhel(-1).and.(cross.eq.0.or.qtype(3).eq.3)) then
        call BoxlineNoAbeTotal(pdia(0,in),kdia(0,1),pdia(0,out),gluVijMom(0,in2,out2),
     &             psi(1,-1,out),psi(1,-1,in),epsW(0,1),gluVij(0,-1,3,in2,out2),
     &             -1,musq,ngaugeh(-3),comp,boxgauge(1,3),boxres(1,3,3,-1),div)
        comp = -1
      endif

      if(lhel(-1).and.(cross.eq.0.or.qtype(3).eq.4)) then
        call BoxlineNoAbeTotal(pdia(0,in),kdia(0,1),pdia(0,out),gluVijMom(0,in2,out2),
     &             psi(1,-1,out),psi(1,-1,in),epsW(0,1),gluVij(0,-1,4,in2,out2),
     &             -1,musq,ngaugeh(-4),comp,boxgauge(1,3),boxres(1,3,4,-1),div)
        comp = -1
      endif

      if(lhel(1).and.(cross.eq.0.or.qtype(3).eq.3)) then
        call BoxlineNoAbeTotal(pdia(0,in),kdia(0,1),pdia(0,out),gluVijMom(0,in2,out2),
     &             psi(1,-1,out),psi(1,-1,in),epsW(0,1),gluVij(0,1,3,in2,out2),
     &             -1,musq,ngaugeh(3),comp,boxgauge(1,3),boxres(1,3,3,1),div)
        comp = -1
      endif

      if(lhel(1).and.(cross.eq.0.or.qtype(3).eq.4)) then
        call BoxlineNoAbeTotal(pdia(0,in),kdia(0,1),pdia(0,out),gluVijMom(0,in2,out2),
     &             psi(1,-1,out),psi(1,-1,in),epsW(0,1),gluVij(0,1,4,in2,out2),
     &             -1,musq,ngaugeh(4),comp,boxgauge(1,3),boxres(1,3,4,1),div)
        comp = -1
      endif

      call BoxlineNoAbeTotal(pdia(0,in),kdia(0,1),pdia(0,out),gluVijMom(0,in2,out2),
     &             psi(1,-1,out),psi(1,-1,in),epsW(0,1),gluVijmomc(0,in2,out2),
     &             -1,musq,9,-1,nothing,checkg(1,3,1),div)

c gauge tests 
      call box_gauge(pdia(0,in),gluijMom(0,in2,out2),kdia(0,3),pdia(0,out),musq,
     2    boxgauge(1,1),boxgaugeb(1,1),wardidtest(1))
      call box_gauge(pdia(0,in),kdia(0,3),gluijMom(0,in2,out2),pdia(0,out),musq,
     2    boxgauge(1,2),boxgaugeb(1,2),wardidtest(2))

c      check epsG -> pG
      if (abs(checkg(1,3,1)).lt.1d-4) then
        td=abs(-(checkg(2,2,1)+checkg(2,1,1))/2d0+checkg(1,3,1))
      else
        td=abs(-(checkg(2,2,1)+checkg(2,1,1))/2d0/checkg(1,3,1)+1d0)
      endif
      wardidtest(3) = td.lt.GaugeAccuracy
      if(ldebug) print*,wardidtest(1),wardidtest(2),wardidtest(3),td,checkg(1,3,1)
      

      if(ldebug) gaugeprecision = max(gaugeprecision,td)
      
c sum
      wardidtest(0) = wardidtest(1).and.wardidtest(2).and.wardidtest(3).or.ldebug
      gaugeOK(1) = gaugeOK(1).and.wardidtest(0)
!             print*,td,abs(checkg(1,3,1))
!       print*,wardidtest(1:3)
!       if(.not.wardidtest(0)) pause
      do h=-1,1,2
      do q=3,4
        if(lhel(h).and.(cross.eq.0.or.qtype(3).eq.q)) then
          tempborn = boxborn(1,q,h)+boxborn(2,q,h)
          tempvirt = CF*(boxres(1,1,q,h)+boxres(1,2,q,h))+CFCA*(boxres(2,1,q,h)+boxres(2,2,q,h))+CA*boxres(1,3,q,h)
          tempborn = tempborn*crosssign*G**2/gluVijMom(4,in2,out2)
          tempvirt = tempvirt*crosssign*G**2/gluVijMom(4,in2,out2)	  
          
          ampb(colstr,q,h) = ampb(colstr,q,h) - tempborn
          ampv(colstr,q,h) = ampv(colstr,q,h) - tempvirt
          if(div.eq.0) ampv(colstr,q,h) = ampv(colstr,q,h) + CF*tempborn ! dim. reg.  
        endif
      enddo   
      enddo   

c lower line boxes
c in> W> out, in2> G> V>out2 
!       print*,'box 3a'
      comp=1
      if(lhel(-1).and.(cross.eq.0.or.qtype(3).eq.3)) then
        call boxlineABETotal(pdia(0,in2),gluWijMom(0,in,out),kdia(0,2),pdia(0,out2),
     &            psi(1,-1,out2),psi(1,-1,in2),gluWij(0,in,out),epsV(0,lephel,3,-1,1),
     &            -1,musq,1,2,ngaugeh(-3),comp,boxgauge(1,1),boxres(1,1,3,-1),boxgaugeb(1,1),boxborn(1,3,-1),div)
        comp =-1
      endif
      
      if(lhel(-1).and.(cross.eq.0.or.qtype(3).eq.4)) then
        call boxlineABETotal(pdia(0,in2),gluWijMom(0,in,out),kdia(0,2),pdia(0,out2),
     &            psi(1,-1,out2),psi(1,-1,in2),gluWij(0,in,out),epsV(0,lephel,4,-1,1),
     &            -1,musq,1,2,ngaugeh(-4),comp,boxgauge(1,1),boxres(1,1,4,-1),boxgaugeb(1,1),boxborn(1,4,-1),div)
        comp =-1
      endif
      
      if(lhel(1).and.(cross.eq.0.or.qtype(3).eq.3)) then
        call boxlineABETotal(pdia(0,in2),gluWijMom(0,in,out),kdia(0,2),pdia(0,out2),
     &            psi(1,1,out2),psi(1,1,in2),gluWij(0,in,out),epsV(0,lephel,3,1,1),
     &            1,musq,1,2,ngaugeh(3),comp,boxgauge(1,1),boxres(1,1,3,1),boxgaugeb(1,1),boxborn(1,3,1),div)
        comp =-1
      endif
      
      if(lhel(1).and.(cross.eq.0.or.qtype(3).eq.4)) then
        call boxlineABETotal(pdia(0,in2),gluWijMom(0,in,out),kdia(0,2),pdia(0,out2),
     &            psi(1,1,out2),psi(1,1,in2),gluWij(0,in,out),epsV(0,lephel,4,1,1),
     &            1,musq,1,2,ngaugeh(4),comp,boxgauge(1,1),boxres(1,1,4,1),boxgaugeb(1,1),boxborn(1,4,1),div)
        comp =-1
      endif

c        gauge test of nAbe boxes
      call boxlineABETotal(pdia(0,in2),gluWijMom(0,in,out),kdia(0,2),pdia(0,out2),
     &            psi(1,gtesthel,out2),psi(1,gtesthel,in2),gluWijMomc(0,in,out),epsV(0,lephel,gtestq,gtesthel,1),
     &            gtesthel,musq,1,2,9,-1,nothing,checkg(1,1,1),nothing,checkgb(1,1),div)

      
c in> W> out, in2> V> G>out2 
!       print*,'box 3b'
      comp=1
      if(lhel(-1).and.(cross.eq.0.or.qtype(3).eq.3)) then
        call boxlineABETotal(pdia(0,in2),kdia(0,2),gluWijMom(0,in,out),pdia(0,out2),
     &            psi(1,-1,out2),psi(1,-1,in2),epsV(0,lephel,3,-1,1),gluWij(0,in,out),
     &            -1,musq,1,3,ngaugeh(-3),comp,boxgauge(1,2),boxres(1,2,3,-1),boxgaugeb(1,2),boxborn(2,3,-1),div)
        comp =-1	
      endif

      if(lhel(-1).and.(cross.eq.0.or.qtype(3).eq.4)) then
        call boxlineABETotal(pdia(0,in2),kdia(0,2),gluWijMom(0,in,out),pdia(0,out2),
     &            psi(1,-1,out2),psi(1,-1,in2),epsV(0,lephel,4,-1,1),gluWij(0,in,out),
     &            -1,musq,1,3,ngaugeh(-4),comp,boxgauge(1,2),boxres(1,2,4,-1),boxgaugeb(1,2),boxborn(2,4,-1),div)
        comp =-1
      endif
      
      if(lhel(1).and.(cross.eq.0.or.qtype(3).eq.3)) then
        call boxlineABETotal(pdia(0,in2),kdia(0,2),gluWijMom(0,in,out),pdia(0,out2),
     &            psi(1,1,out2),psi(1,1,in2),epsV(0,lephel,3,1,1),gluWij(0,in,out),
     &            1,musq,1,3,ngaugeh(3),comp,boxgauge(1,2),boxres(1,2,3,1),boxgaugeb(1,2),boxborn(2,3,1),div)
        comp =-1
      endif
      
      if(lhel(1).and.(cross.eq.0.or.qtype(3).eq.4)) then
        call boxlineABETotal(pdia(0,in2),kdia(0,2),gluWijMom(0,in,out),pdia(0,out2),
     &            psi(1,1,out2),psi(1,1,in2),epsV(0,lephel,4,1,1),gluWij(0,in,out),
     &            1,musq,1,3,ngaugeh(4),comp,boxgauge(1,2),boxres(1,2,4,1),boxgaugeb(1,2),boxborn(2,4,1),div)
        comp =-1
      endif

c        gauge test of nAbe boxes
      call boxlineABETotal(pdia(0,in2),kdia(0,2),gluWijMom(0,in,out),pdia(0,out2),
     &            psi(1,gtesthel,out2),psi(1,gtesthel,in2),epsV(0,lephel,gtestq,gtesthel,1),gluWijMomc(0,in,out),
     &            gtesthel,musq,1,3,9,-1,nothing,checkg(1,2,1),nothing,checkgb(2,1),div)
      
c non abelian
!       print*,'box 3c'
      comp=1
      if(lhel(-1).and.(cross.eq.0.or.qtype(3).eq.3)) then
        call BoxlineNoAbeTotal(pdia(0,in2),kdia(0,2),pdia(0,out2),gluWijMom(0,in,out),
     &             psi(1,-1,out2),psi(1,-1,in2),epsV(0,lephel,3,-1,1),gluWij(0,in,out),
     &             -1,musq,ngaugeh(-3),comp,boxgauge(1,3),boxres(1,3,3,-1),div)
        comp = -1
      endif

      if(lhel(-1).and.(cross.eq.0.or.qtype(3).eq.4)) then
        call BoxlineNoAbeTotal(pdia(0,in2),kdia(0,2),pdia(0,out2),gluWijMom(0,in,out),
     &             psi(1,-1,out2),psi(1,-1,in2),epsV(0,lephel,4,-1,1),gluWij(0,in,out),
     &             -1,musq,ngaugeh(-4),comp,boxgauge(1,3),boxres(1,3,4,-1),div)
        comp = -1
      endif

      if(lhel(1).and.(cross.eq.0.or.qtype(3).eq.3)) then
        call BoxlineNoAbeTotal(pdia(0,in2),kdia(0,2),pdia(0,out2),gluWijMom(0,in,out),
     &             psi(1,1,out2),psi(1,1,in2),epsV(0,lephel,3,1,1),gluWij(0,in,out),
     &             1,musq,ngaugeh(3),comp,boxgauge(1,3),boxres(1,3,3,1),div)
        comp = -1
      endif

      if(lhel(1).and.(cross.eq.0.or.qtype(3).eq.4)) then
        call BoxlineNoAbeTotal(pdia(0,in2),kdia(0,2),pdia(0,out2),gluWijMom(0,in,out),
     &             psi(1,1,out2),psi(1,1,in2),epsV(0,lephel,4,1,1),gluWij(0,in,out),
     &             1,musq,ngaugeh(4),comp,boxgauge(1,3),boxres(1,3,4,1),div)
        comp = -1
      endif

      call BoxlineNoAbeTotal(pdia(0,in2),kdia(0,2),pdia(0,out2),gluWijMom(0,in,out),
     &             psi(1,gtesthel,out2),psi(1,gtesthel,in2),epsV(0,lephel,gtestq,gtesthel,1),gluWijmomc(0,in,out),
     &             gtesthel,musq,9,-1,nothing,checkg(1,3,1),div)


c gauge tests 
      call box_gauge(pdia(0,in2),gluWijMom(0,in,out),kdia(0,2),pdia(0,out2),musq,
     2    boxgauge(1,1),boxgaugeb(1,1),wardidtest(1))
      call box_gauge(pdia(0,in2),kdia(0,2),gluWijMom(0,in,out),pdia(0,out2),musq,
     2    boxgauge(1,2),boxgaugeb(1,2),wardidtest(2))

c      check epsG -> pG
      if (abs(checkg(1,3,1)).lt.1d-4) then
        td=abs(-(checkg(2,2,1)+checkg(2,1,1))/2d0+checkg(1,3,1))
      else
        td=abs(-(checkg(2,2,1)+checkg(2,1,1))/2d0/checkg(1,3,1)+1d0)
      endif
      wardidtest(3) = td.lt.GaugeAccuracy
      if(ldebug) print*,wardidtest(1),wardidtest(2),wardidtest(3),td,checkg(1,3,1)
!       pause

      if(ldebug) gaugeprecision = max(gaugeprecision,td)
      
c sum
      wardidtest(0) = wardidtest(1).and.wardidtest(2).and.wardidtest(3).or.ldebug
      gaugeOK(1) = gaugeOK(1).and.wardidtest(0)
!             print*,td,abs(checkg(1,3,1))
!       print*,wardidtest(1:3)
!       if(.not.wardidtest(0)) pause
      do h=-1,1,2
      do q=3,4
        if(lhel(h).and.(cross.eq.0.or.qtype(3).eq.q)) then
          tempborn = boxborn(1,q,h)+boxborn(2,q,h)
          tempvirt = CF*(boxres(1,1,q,h)+boxres(1,2,q,h))+CFCA*(boxres(2,1,q,h)+boxres(2,2,q,h))+CA*boxres(1,3,q,h)
          tempborn = tempborn*crosssign*G**2/gluVijMom(4,in2,out2)
          tempvirt = tempvirt*crosssign*G**2/gluVijMom(4,in2,out2)	  
          
          ampb(colstr,q,h) = ampb(colstr,q,h) - tempborn
          ampv(colstr,q,h) = ampv(colstr,q,h) - tempvirt
          if(div.eq.0) ampv(colstr,q,h) = ampv(colstr,q,h) + CF*tempborn ! dim. reg.  
        endif
      enddo   
      enddo   

      boxcount2m = boxcount2m + 1
      boxcount2sub(sub_number) = boxcount2sub(sub_number) + 1

      if(.not.gaugeOK(1)) then
        boxcountm   = boxcountm   + 1
        boxcountm_q = boxcountm_q + 1
        boxcountsub  (sub_number) = boxcountsub  (sub_number) + 1
        boxcountsub_q(sub_number) = boxcountsub_q(sub_number) + 1
        if(return0ifGaugeFails) then
          do i=0,2; do j=3,4; do k=0,2
            res(i,j,k) = 0d0
          enddo; enddo; enddo
          return
        else
          recalculatePSpoint=.true.
        endif
      endif

      endif ! lbox


c ---------------------------------------
c $$$$$$$$$$ PENTAGONS $$$$$$$$$$$$$$$$$$
c ---------------------------------------


      if(lpenl.and.lWV2) then
        gaugeOK(1) = .true.
        gaugeOK(2) = .true.

      if(ldebug) then
        print*
        print*,'pent'
      endif
      
      ngaugeh(gtesthel) = 4
      ngaugeh(-gtesthel) = 9


      do i=1,8
        penborn(i,-1)=0d0
        penborn(i, 1)=0d0
      enddo


      if(lhel(-1)) then
         ngaugeh(-1)=4
         ngaugeh(1)=9
      endif


c in> G> W> V> out, in2>out2 
      comp=1
      if(lhel(-1)) then
        call penlineABETotalQ(pdia(0,in),gluijMom(0,in2,out2),kdia(0,1),kdia(0,2),pdia(0,out),
     &            psi(1,-1,out),psi(1,-1,in),gluij(0,-1,in2,out2),epsW(0,1),epsV(0,lephel,qtype(2),-1,1),
     &            -1,musq,1,2,ngaugeh(-1),comp,pengauge(1,1),penres(1,1,-1),pengaugeb(1,1),penborn(1,-1),div)!,wardidtest(1))
        comp =-1
      endif
      
      if(lhel(1)) then
        call penlineABETotalQ(pdia(0,in),gluijMom(0,in2,out2),kdia(0,1),kdia(0,2),pdia(0,out),
     &            psi(1,-1,out),psi(1,-1,in),gluij(0, 1,in2,out2),epsW(0,1),epsV(0,lephel,qtype(2),-1,1),
     &            -1,musq,1,2,ngaugeh(1),comp,pengauge(1,1),penres(1,1,1),pengaugeb(1,1),penborn(1,1),div)!,wardidtest(1))
        comp =-1
      endif
      wardidtest(1)=wardpAbe(1)
      wardidtest_Q(1)=wardpAbe(2)     

c in> W> G> V> out, in2>out2 
      comp=1
      if(lhel(-1)) then
        call penlineABETotalQ(pdia(0,in),kdia(0,1),gluijMom(0,in2,out2),kdia(0,2),pdia(0,out),
     &            psi(1,-1,out),psi(1,-1,in),epsW(0,1),gluij(0,-1,in2,out2),epsV(0,lephel,qtype(2),-1,1),
     &            -1,musq,1,3,ngaugeh(-1),comp,pengauge(1,2),penres(1,2,-1),pengaugeb(1,2),penborn(2,-1),div)!,wardidtest(2))
        comp =-1
      endif
      
      if(lhel(1)) then
        call penlineABETotalQ(pdia(0,in),kdia(0,1),gluijMom(0,in2,out2),kdia(0,2),pdia(0,out),
     &            psi(1,-1,out),psi(1,-1,in),epsW(0,1),gluij(0, 1,in2,out2),epsV(0,lephel,qtype(2),-1,1),
     &            -1,musq,1,3,ngaugeh(1),comp,pengauge(1,2),penres(1,2,1),pengaugeb(1,2),penborn(2,1),div)!,wardidtest(2))
        comp =-1
      endif
      wardidtest(2)=wardpAbe(1)
      wardidtest_Q(2)=wardpAbe(2)     
      
c in> W> V> G> out, in2>out2 
      comp=1
      if(lhel(-1)) then
        call penlineABETotalQ(pdia(0,in),kdia(0,1),kdia(0,2),gluijMom(0,in2,out2),pdia(0,out),
     &            psi(1,-1,out),psi(1,-1,in),epsW(0,1),epsV(0,lephel,qtype(2),-1,1),gluij(0,-1,in2,out2),
     &            -1,musq,1,4,ngaugeh(-1),comp,pengauge(1,3),penres(1,3,-1),pengaugeb(1,3),penborn(3,-1),div)!,wardidtest(3))
        comp =-1
      endif
      
      if(lhel(1)) then
        call penlineABETotalQ(pdia(0,in),kdia(0,1),kdia(0,2),gluijMom(0,in2,out2),pdia(0,out),
     &            psi(1,-1,out),psi(1,-1,in),epsW(0,1),epsV(0,lephel,qtype(2),-1,1),gluij(0, 1,in2,out2),
     &            -1,musq,1,4,ngaugeh(1),comp,pengauge(1,3),penres(1,3,1),pengaugeb(1,3),penborn(3,1),div)!,wardidtest(3))
        comp =-1
      endif
      wardidtest(3)=wardpAbe(1)
      wardidtest_Q(3)=wardpAbe(2)     

c in> G> V> W> out, in2>out2 
      comp=1
      if(lhel(-1)) then
        call penlineABETotalQ(pdia(0,in),gluijMom(0,in2,out2),kdia(0,2),kdia(0,1),pdia(0,out),
     &            psi(1,-1,out),psi(1,-1,in),gluij(0,-1,in2,out2),epsV(0,lephel,qtype(1),-1,1),epsW(0,1),
     &            -1,musq,1,2,ngaugeh(-1),comp,pengauge(1,4),penres(1,4,-1),pengaugeb(1,4),penborn(4,-1),div)!,wardidtest(4))
        comp =-1
      endif
      
      if(lhel(1)) then
        call penlineABETotalQ(pdia(0,in),gluijMom(0,in2,out2),kdia(0,2),kdia(0,1),pdia(0,out),
     &            psi(1,-1,out),psi(1,-1,in),gluij(0, 1,in2,out2),epsV(0,lephel,qtype(1),-1,1),epsW(0,1),
     &            -1,musq,1,2,ngaugeh(1),comp,pengauge(1,4),penres(1,4,1),pengaugeb(1,4),penborn(4,1),div)!,wardidtest(4))
        comp =-1
      endif
      wardidtest(4)=wardpAbe(1)
      wardidtest_Q(4)=wardpAbe(2)     

c in> V> G> W> out, in2>out2 
      comp=1
      if(lhel(-1)) then
        call penlineABETotalQ(pdia(0,in),kdia(0,2),gluijMom(0,in2,out2),kdia(0,1),pdia(0,out),
     &            psi(1,-1,out),psi(1,-1,in),epsV(0,lephel,qtype(1),-1,1),gluij(0,-1,in2,out2),epsW(0,1),
     &            -1,musq,1,3,ngaugeh(-1),comp,pengauge(1,5),penres(1,5,-1),pengaugeb(1,5),penborn(5,-1),div)!,wardidtest(5))
        comp =-1
      endif
      
      if(lhel(1)) then
        call penlineABETotalQ(pdia(0,in),kdia(0,2),gluijMom(0,in2,out2),kdia(0,1),pdia(0,out),
     &            psi(1,-1,out),psi(1,-1,in),epsV(0,lephel,qtype(1),-1,1),gluij(0, 1,in2,out2),epsW(0,1),
     &            -1,musq,1,3,ngaugeh(1),comp,pengauge(1,5),penres(1,5,1),pengaugeb(1,5),penborn(5,1),div)!,wardidtest(5))
        comp =-1
      endif
      wardidtest(5)=wardpAbe(1)
      wardidtest_Q(5)=wardpAbe(2)     

c in> V> W> G> out, in2>out2 
      comp=1
      if(lhel(-1)) then
        call penlineABETotalQ(pdia(0,in),kdia(0,2),kdia(0,1),gluijMom(0,in2,out2),pdia(0,out),
     &            psi(1,-1,out),psi(1,-1,in),epsV(0,lephel,qtype(1),-1,1),epsW(0,1),gluij(0,-1,in2,out2),
     &            -1,musq,1,4,ngaugeh(-1),comp,pengauge(1,6),penres(1,6,-1),pengaugeb(1,6),penborn(6,-1),div)!,wardidtest(6))
        comp =-1
      endif
      
      if(lhel(1)) then
        call penlineABETotalQ(pdia(0,in),kdia(0,2),kdia(0,1),gluijMom(0,in2,out2),pdia(0,out),
     &            psi(1,-1,out),psi(1,-1,in),epsV(0,lephel,qtype(1),-1,1),epsW(0,1),gluij(0, 1,in2,out2),
     &            -1,musq,1,4,ngaugeh(1),comp,pengauge(1,6),penres(1,6,1),pengaugeb(1,6),penborn(6,1),div)!,wardidtest(6))
        comp =-1
      endif
      wardidtest(6)=wardpAbe(1)
      wardidtest_Q(6)=wardpAbe(2)     
c
c  no Abe
c
c      ngaugeh(gtesthel) = 4
c      ngaugeh(-gtesthel) = 9

c in> W> V> out
      comp=1
      if(lhel(-1)) then
        call penlineNoABETotalQ(pdia(0,in),kdia(0,1),kdia(0,2),pdia(0,out),gluijMom(0,in2,out2),
     &            psi(1,-1,out),psi(1,-1,in),epsW(0,1),epsV(0,lephel,qtype(2),-1,1),gluij(0,-1,in2,out2),
     &            -1,musq,ngaugeh(-1),comp,pengauge(1,7),penres(1,7,-1),div)
        comp =-1
      endif

      if(lhel(1)) then
        call penlineNoABETotalQ(pdia(0,in),kdia(0,1),kdia(0,2),pdia(0,out),gluijMom(0,in2,out2),
     &            psi(1,-1,out),psi(1,-1,in),epsW(0,1),epsV(0,lephel,qtype(2),-1,1),gluij(0,1,in2,out2),
     &            -1,musq,ngaugeh(1),comp,pengauge(1,7),penres(1,7,1),div)
        comp =-1
      endif
      wardidtest(7)=wardpNoAbe(1)
      wardidtest_Q(7)=wardpNoAbe(2)     

c in> V> W> out
      comp=1
      if(lhel(-1)) then
        call penlineNoABETotalQ(pdia(0,in),kdia(0,2),kdia(0,1),pdia(0,out),gluijMom(0,in2,out2),
     &            psi(1,-1,out),psi(1,-1,in),epsV(0,lephel,qtype(1),-1,1),epsW(0,1),gluij(0,-1,in2,out2),
     &            -1,musq,ngaugeh(-1),comp,pengauge(1,8),penres(1,8,-1),div)
        comp =-1
      endif

      if(lhel(1)) then
        call penlineNoABETotalQ(pdia(0,in),kdia(0,2),kdia(0,1),pdia(0,out),gluijMom(0,in2,out2),
     &            psi(1,-1,out),psi(1,-1,in),epsV(0,lephel,qtype(1),-1,1),epsW(0,1),gluij(0,1,in2,out2),
     &            -1,musq,ngaugeh(1),comp,pengauge(1,8),penres(1,8,1),div)
        comp =-1
      endif
      wardidtest(8)=wardpNoabe(1)
      wardidtest_Q(8)=wardpNoAbe(2)     

      wardidtest(0)=wardidtest(1).and.wardidtest(2).and.wardidtest(3).and.wardidtest(4).and.wardidtest(5)
     &              .and.wardidtest(6).and.wardidtest(7).and.wardidtest(8)
      wardidtest_q(0)=wardidtest_q(1).and.wardidtest_q(2).and.wardidtest_q(3).and.wardidtest_q(4).and.wardidtest_q(5)
     &              .and.wardidtest_q(6).and.wardidtest_q(7).and.wardidtest_q(8)

      wardidtest_q(0) = wardidtest_q(0).or.ldebug
      gaugeOK(1) = gaugeOK(1) .and.wardidtest(0)
      gaugeOK(2) = gaugeOK(2) .and.wardidtest_Q(0)
      
      if(ldebug) print*,wardidtest
      
c sum
      do h=-1,1,2
        if(lhel(h)) then
          tempborn = penborn(1,h)+penborn(2,h)+penborn(3,h)+penborn(4,h)+penborn(5,h)+penborn(6,h)
          tempvirt = CF*(penres(1,1,h)+penres(1,2,h)+penres(1,3,h)+penres(1,4,h)+penres(1,5,h)+penres(1,6,h))
     &		    +CFCA*(penres(2,1,h)+penres(2,2,h)+penres(2,3,h)+penres(2,4,h)+penres(2,5,h)+penres(2,6,h))
     &		    +CA*(penres(1,7,h)+penres(1,8,h))
          tempborn = tempborn*crosssign*G**2/gluijMom(4,in2,out2)
          tempvirt = tempvirt*crosssign*G**2/gluijMom(4,in2,out2)	  
          
          if(cross.eq.0.or.qtype(3).eq.3) then
            ampb(colstr,3,h) = ampb(colstr,3,h) - tempborn
            ampv(colstr,3,h) = ampv(colstr,3,h) - tempvirt
            if(div.eq.0) ampv(colstr,3,h) = ampv(colstr,3,h) + CF*tempborn ! dim. reg.  
          endif
          if(cross.eq.0.or.qtype(3).eq.4) then
            ampb(colstr,4,h) = ampb(colstr,4,h) - tempborn
            ampv(colstr,4,h) = ampv(colstr,4,h) - tempvirt
            if(div.eq.0) ampv(colstr,4,h) = ampv(colstr,4,h) + CF*tempborn ! dim. reg.  
          endif
        endif
      enddo

      
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
          do i=0,2; do j=3,4; do k=0,2
            res(i,j,k) = 0d0
          enddo; enddo; enddo
          return
        else
          recalculatePSpoint=.true.
        endif
      endif

      endif ! lpenline

      
c ---------------------------------------
c $$$$$$$$$$ PENBOX  $$$$$$$$$$$$$$$$$$$$
c ---------------------------------------
      
      if(lpenb) then
        gaugeOK(1) = .true.
      
      psi4(1,1) = psi(1,-1,1)
      psi4(2,1) = psi(2,-1,1)
      psi4(3,1) = psi(1, 1,1)
      psi4(4,1) = psi(2, 1,1)

      psi4(3,2) = psi(1,-1,2)
      psi4(4,2) = psi(2,-1,2)
      psi4(1,2) = psi(1, 1,2)
      psi4(2,2) = psi(2, 1,2)

      psi4(1,3) = psi(1,-1,3)
      psi4(2,3) = psi(2,-1,3)
      psi4(3,3) = psi(1, 1,3)
      psi4(4,3) = psi(2, 1,3)

      psi4(3,4) = psi(1,-1,4)
      psi4(4,4) = psi(2,-1,4)
      psi4(1,4) = psi(1, 1,4)
      psi4(2,4) = psi(2, 1,4)


c in2  -->---$$$$$---<-- in
c           |     |__W
c           |     |
c out2 --<---$$$$$--->-- out


c penbox  Dire
      if(ldebug) then
        print*
        print*, 'penbox Dire, in,out',in,out
      endif
      
c in > W > Z > out
c start with gauge test
      temphel = 67
      if(gtesthel.eq.-1) temphel = 77
      call PenBoxT(pdia(0,out2),pdia(0,out),kdia(0,3),pdia(0,in),pdia(0,in2),
     &               psi4(1,out2),psi4(1,in2),psi4(1,out),psi4(1,in),kdiac(0,3),
     &               musq,1,0,0,penbamp,penbres,penbborn,div,temphel)
     
      if (abs(penbamp(1)).gt.abs(penbamp(2))) then
        td = abs((penbamp(2)+penbamp(3))/penbamp(1)+1d0)
      else
        td = abs((penbamp(1)+penbamp(3))/penbamp(2)+1d0)
      endif
!       if(abs(penbamp(1)).lt.1d-6.or.abs(penbamp(2)).lt.1d-6) print*,__LINE__,penbamp
      
      wardidtest(0) = td.lt. GaugeAccuracy.or.ldebug
      GaugeOK(1) = GaugeOK(1).and.wardidtest(0)
      if(ldebug) print*, wardidtest(0),td
      if(ldebug) gaugeprecision = max(gaugeprecision,td)
      
c passed test > calc amp
      do h=-1,1,2
      if(lhel(h)) then
        temphel = 67
        if(h.eq.-1) temphel = 77
        call PenBoxT(pdia(0,out2),pdia(0,out),kdia(0,3),pdia(0,in),pdia(0,in2),
     &               psi4(1,out2),psi4(1,in2),psi4(1,out),psi4(1,in),epsWV(0,lephel,1),
     &               musq,-1,0,0,penbamp,penbres,penbborn,div,temphel)
        penbborn = penbborn*GG(1)**2*crosssign
        penbres = penbres*GG(1)**2*crosssign
        if(cross.eq.0.or.qtype(3).eq.3) then
          ampb(col(1),3,h) = ampb(col(1),3,h) + penbborn
          ampv(col(1),3,h) = ampv(col(1),3,h) - penbres/6d0
          ampv(col(2),3,h) = ampv(col(2),3,h) + penbres/2d0
        endif
        if(cross.eq.0.or.qtype(3).eq.4) then
          ampb(col(1),4,h) = ampb(col(1),4,h) + penbborn
          ampv(col(1),4,h) = ampv(col(1),4,h) - penbres/6d0
          ampv(col(2),4,h) = ampv(col(2),4,h) + penbres/2d0
        endif
      endif
      enddo


c penbox  Cross
      if(ldebug) then
        print*
        print*, 'penbox Cross, in,out',in,out
      endif
      
c in > W > Z > out
c start with gauge test
      temphel = 67
      if(gtesthel.eq.-1) temphel = 77
      call PenBoxCrossT(pdia(0,in2),pdia(0,out),kdia(0,3),pdia(0,in),pdia(0,out2),
     &               psi4(1,in2),psi4(1,out2),psi4(1,out),psi4(1,in),kdiac(0,3),
     &               musq,1,0,0,penbamp,penbres,penbborn,div,temphel)

      if (abs(penbamp(1)).gt.abs(penbamp(2))) then
        td = abs((penbamp(2)+penbamp(3))/penbamp(1)+1d0)
      else
        td = abs((penbamp(1)+penbamp(3))/penbamp(2)+1d0)
      endif
!       if(abs(penbamp(1)).lt.1d-6.or.abs(penbamp(2)).lt.1d-6) print*,__LINE__,penbamp
      
      wardidtest(0) = td.lt.GaugeAccuracy.or.ldebug
      GaugeOK(1) = GaugeOK(1).and.wardidtest(0)
      if(ldebug) print*, wardidtest(0),td
      if(ldebug) gaugeprecision = max(gaugeprecision,td)
      
c passed test > calc amp
      do h=-1,1,2
      if(lhel(h)) then
        temphel = 67
        if(h.eq.-1) temphel = 77
        call PenBoxCrossT(pdia(0,in2),pdia(0,out),kdia(0,3),pdia(0,in),pdia(0,out2),
     &               psi4(1,in2),psi4(1,out2),psi4(1,out),psi4(1,in),epsWV(0,lephel,1),
     &               musq,-1,0,0,penbamp,penbres,penbborn,div,temphel)
        penbborn = penbborn*GG(1)**2*crosssign
        penbres = penbres*GG(1)**2*crosssign
        if(cross.eq.0.or.qtype(3).eq.3) then
          ampb(col(1),3,h) = ampb(col(1),3,h) + penbborn
          ampv(col(1),3,h) = ampv(col(1),3,h) + penbres*4d0/3d0
          ampv(col(2),3,h) = ampv(col(2),3,h) + penbres/2d0
        endif
        if(cross.eq.0.or.qtype(3).eq.4) then
          ampb(col(1),4,h) = ampb(col(1),4,h) + penbborn
          ampv(col(1),4,h) = ampv(col(1),4,h) + penbres*4d0/3d0
          ampv(col(2),4,h) = ampv(col(2),4,h) + penbres/2d0
        endif
      endif
      enddo

      box2count2m = box2count2m + 1
      box2count2sub(sub_number) = box2count2sub(sub_number) + 1
      if(.not.gaugeOK(1)) then
        box2countm   = box2countm   + 1
        box2countm_Q = box2countm_Q + 1
        box2countsub  (sub_number) = box2countsub  (sub_number) + 1
        box2countsub_q(sub_number) = box2countsub_q(sub_number) + 1
        if(return0ifGaugeFails) then
          do i=0,2; do j=3,4; do k=0,2
            res(i,j,k) = 0d0
          enddo; enddo; enddo
          return
        else
          recalculatePSpoint=.true.
        endif
      endif
      
      endif ! lpenbox

c ---------------------------------------
c $$$$$$$$$$ HEXAGONS  $$$$$$$$$$$$$$$$$$
c ---------------------------------------
      
      if(lhex.and.(lW.or.lWV2)) then
        gaugeOK(1)=.true.
      
      psi4(1,1) = psi(1,-1,1)
      psi4(2,1) = psi(2,-1,1)
      psi4(3,1) = psi(1, 1,1)
      psi4(4,1) = psi(2, 1,1)

      psi4(3,2) = psi(1,-1,2)
      psi4(4,2) = psi(2,-1,2)
      psi4(1,2) = psi(1, 1,2)
      psi4(2,2) = psi(2, 1,2)

      psi4(1,3) = psi(1,-1,3)
      psi4(2,3) = psi(2,-1,3)
      psi4(3,3) = psi(1, 1,3)
      psi4(4,3) = psi(2, 1,3)

      psi4(3,4) = psi(1,-1,4)
      psi4(4,4) = psi(2,-1,4)
      psi4(1,4) = psi(1, 1,4)
      psi4(2,4) = psi(2, 1,4)
      
   
c in2  -->---$$$$$---<-- in
c           |     |
c           |     |
c out2 --<---$$$$$--->-- out


c hexbox 1 Dire
      if(ldebug) then
        print*
        print*, 'hex1 Dire, in,out',in,out
      endif
      
      do i=1,2
      do h=-1,1,2
      tempHexRes(i,h) = 0d0
      tempHexResB(i,h) = 0d0
      enddo
      tempHexResG(i,1) = 0d0
      tempHexResG(i,2) = 0d0
      enddo

c in > W > Z > out
c start with gauge test, by now only V-> pV
      temphel = 67
      if(gtesthel.eq.-1) temphel = 77

      call HexBoxT(pdia(0,out2),pdia(0,out),kdia(0,2),kdia(0,1),pdia(0,in),pdia(0,in2),
     &               psi4(1,out2),psi4(1,in2),psi4(1,out),psi4(1,in),kdiac(0,2),epsW(0,1),
     &               musq,1,0,0,hexamp,hexres,hexborn,div,temphel)
      tempHexResG(1,1) = hexres

      call HexBoxT(pdia(0,out2),pdia(0,out),kdia(0,2),kdia(0,1),pdia(0,in),pdia(0,in2),
     &               psi4(1,out2),psi4(1,in2),psi4(1,out),psi4(1,in),epsV(0,lephel,qtype(2),-1,1),kdiac(0,1),
     &               musq,-1,0,0,hexamp,hexres,hexborn,div,temphel)
      tempHexResG(1,2) = hexres

c hel2 = -1
      if(lhel(-1)) then
        call HexBox77T(pdia(0,out2),pdia(0,out),kdia(0,2),kdia(0,1),pdia(0,in),pdia(0,in2),
     &               psi4(1,out2),psi4(1,in2),psi4(1,out),psi4(1,in),epsV(0,lephel,qtype(2),-1,1),epsW(0,1),
     &               musq,-1,0,0,hexamp,hexres,hexborn,div)
        tempHexRes(1,-1) = hexres
        tempHexResB(1,-1) = hexborn
      endif

c hel2= +1
      if(lhel(1)) then
        call HexBox67T(pdia(0,out2),pdia(0,out),kdia(0,2),kdia(0,1),pdia(0,in),pdia(0,in2),
     &               psi4(1,out2),psi4(1,in2),psi4(1,out),psi4(1,in),epsV(0,lephel,qtype(2),-1,1),epsW(0,1),
     &               musq,-1,0,0,hexamp,hexres,hexborn,div)
        tempHexRes(1,1) = hexres
        tempHexResB(1,1) = hexborn
      endif

c in > Z > W > out,  hel2=-1
c start with gauge test, by now only V-> pV
        call HexBoxT(pdia(0,out2),pdia(0,out),kdia(0,1),kdia(0,2),pdia(0,in),pdia(0,in2),
     &               psi4(1,out2),psi4(1,in2),psi4(1,out),psi4(1,in),epsW(0,1),kdiac(0,2),
     &               musq,1,0,0,hexamp,hexres,hexborn,div,temphel)
        tempHexResG(2,1) = hexres

        call HexBoxT(pdia(0,out2),pdia(0,out),kdia(0,1),kdia(0,2),pdia(0,in),pdia(0,in2),
     &               psi4(1,out2),psi4(1,in2),psi4(1,out),psi4(1,in),kdiac(0,1),epsV(0,lephel,qtype(2),-1,1),
     &               musq,-1,0,0,hexamp,hexres,hexborn,div,temphel)
        tempHexResG(2,2) = hexres

c hel2 = -1
      if(lhel(-1)) then
        call HexBox77T(pdia(0,out2),pdia(0,out),kdia(0,1),kdia(0,2),pdia(0,in),pdia(0,in2),
     &               psi4(1,out2),psi4(1,in2),psi4(1,out),psi4(1,in),epsW(0,1),epsV(0,lephel,qtype(1),-1,1),
     &               musq,-1,0,0,hexamp,hexres,hexborn,div)
        tempHexRes(2,-1) = hexres
        tempHexResB(2,-1) = hexborn
      endif

c hel2= 1
      if(lhel(1)) then
        call HexBox67T(pdia(0,out2),pdia(0,out),kdia(0,1),kdia(0,2),pdia(0,in),pdia(0,in2),
     &               psi4(1,out2),psi4(1,in2),psi4(1,out),psi4(1,in),epsW(0,1),epsV(0,lephel,qtype(1),-1,1),
     &               musq,-1,0,0,hexamp,hexres,hexborn,div)
        tempHexRes(2,1) = hexres
        tempHexResB(2,1) = hexborn
      endif

      if(ldebug) then
       gaugeprecision=max(gaugeprecision,abs((tempHexResG(1,1)+tempHexResG(2,1))/max(abs(tempHexResG(1,1)),abs(tempHexResG(2,1)))))
       gaugeprecision=max(gaugeprecision,abs((tempHexResG(1,2)+tempHexResG(2,2))/max(abs(tempHexResG(1,2)),abs(tempHexResG(2,2)))))
       call check('gtest HexB1',-tempHexResG(1,1)/tempHexResG(2,1))
      endif
      
      wardidtest(1) = abs((tempHexResG(1,1)+tempHexResG(2,1))/tempHexResG(1,1)).lt.GaugeAccuracy
      wardidtest(2) = abs((tempHexResG(1,2)+tempHexResG(2,2))/tempHexResG(1,2)).lt.GaugeAccuracy
!       if(abs(tempHexResG(1,1)).lt.1d-6.or.abs(tempHexResG(1,1)).lt.1d-6) print*,__LINE__,tempHexResG
      wardidtest(0) = wardidtest(1).and.wardidtest(2).or.ldebug
      do q=3,4
      do h=-1,1,2
      if(lhel(h).and.(cross.eq.0.or.qtype(3).eq.q)) then
        ampb(col(1),q,h) = ampb(col(1),q,h) + (tempHexResB(1,h)+tempHexResB(2,h))*GG(1)**2*crosssign
        ampv(col(1),q,h) = ampv(col(1),q,h) - (tempHexRes(1,h)+tempHexRes(2,h))*GG(1)**2*crosssign/6d0
        ampv(col(2),q,h) = ampv(col(2),q,h) + (tempHexRes(1,h)+tempHexRes(2,h))*GG(1)**2*crosssign/2d0
      endif
      enddo;enddo
      if(.NOT.wardidtest(0)) then
        gaugeOK(1) = .false.
      endif ! wardidtest


c hexbox 1 Cross
      if(ldebug) then
        print*
        print*, 'hex1 Cross, in,out',in,out
      endif
      
      do i=1,2
      do h=-1,1,2
      tempHexRes(i,h) = 0d0
      tempHexResB(i,h) = 0d0
      enddo
      tempHexResG(i,1) = 0d0
      tempHexResG(i,2) = 0d0
      enddo
      
c in > W > Z > out
c star with gauge test, by now only V-> pV
      call HexBoxCross1T(pdia(0,in2),pdia(0,out),kdia(0,2),kdia(0,1),pdia(0,in),pdia(0,out2),
     &               psi4(1,in2),psi4(1,out2),psi4(1,out),psi4(1,in),kdiac(0,2),epsW(0,1),
     &               musq,1,0,0,hexamp,hexres,hexborn,div,temphel)
      tempHexResG(1,1) = hexres

      call HexBoxCross1T(pdia(0,in2),pdia(0,out),kdia(0,2),kdia(0,1),pdia(0,in),pdia(0,out2),
     &               psi4(1,in2),psi4(1,out2),psi4(1,out),psi4(1,in),epsV(0,lephel,qtype(2),-1,1),kdiac(0,1),
     &               musq,-1,0,0,hexamp,hexres,hexborn,div,temphel)
      tempHexResG(1,2) = hexres

c hel2=-1
      if(lhel(-1)) then
        call HexBoxCross177T(pdia(0,in2),pdia(0,out),kdia(0,2),kdia(0,1),pdia(0,in),pdia(0,out2),
     &               psi4(1,in2),psi4(1,out2),psi4(1,out),psi4(1,in),epsV(0,lephel,qtype(2),-1,1),epsW(0,1),
     &               musq,-1,0,0,hexamp,hexres,hexborn,div)
        tempHexRes(1,-1) = hexres
        tempHexResB(1,-1) = hexborn
      endif

c hel2=+1
      if(lhel(1)) then
        call HexBoxCross167T(pdia(0,in2),pdia(0,out),kdia(0,2),kdia(0,1),pdia(0,in),pdia(0,out2),
     &               psi4(1,in2),psi4(1,out2),psi4(1,out),psi4(1,in),epsV(0,lephel,qtype(2),-1,1),epsW(0,1),
     &               musq,-1,0,0,hexamp,hexres,hexborn,div)
        tempHexRes(1,1) = hexres
        tempHexResB(1,1) = hexborn
      endif


c in > Z > W > out
c star with gauge test, by now only V-> pV
        call HexBoxCross1T(pdia(0,in2),pdia(0,out),kdia(0,1),kdia(0,2),pdia(0,in),pdia(0,out2),
     &               psi4(1,in2),psi4(1,out2),psi4(1,out),psi4(1,in),epsW(0,1),kdiac(0,2),
     &               musq,1,0,0,hexamp,hexres,hexborn,div,temphel)
        tempHexResG(2,1) = hexres

        call HexBoxCross1T(pdia(0,in2),pdia(0,out),kdia(0,1),kdia(0,2),pdia(0,in),pdia(0,out2),
     &               psi4(1,in2),psi4(1,out2),psi4(1,out),psi4(1,in),kdiac(0,1),epsV(0,lephel,qtype(2),-1,1),
     &               musq,1,0,0,hexamp,hexres,hexborn,div,temphel)
        tempHexResG(2,2) = hexres

c hel2=-1
      if(lhel(-1)) then
        call HexBoxCross177T(pdia(0,in2),pdia(0,out),kdia(0,1),kdia(0,2),pdia(0,in),pdia(0,out2),
     &               psi4(1,in2),psi4(1,out2),psi4(1,out),psi4(1,in),epsW(0,1),epsV(0,lephel,qtype(1),-1,1),
     &               musq,-1,0,0,hexamp,hexres,hexborn,div)
        tempHexRes(2,-1) = hexres
        tempHexResB(2,-1) = hexborn
      endif

c hel2=+1
      if(lhel(1)) then
        call HexBoxCross167T(pdia(0,in2),pdia(0,out),kdia(0,1),kdia(0,2),pdia(0,in),pdia(0,out2),
     &               psi4(1,in2),psi4(1,out2),psi4(1,out),psi4(1,in),epsW(0,1),epsV(0,lephel,qtype(1),-1,1),
     &               musq,-1,0,0,hexamp,hexres,hexborn,div)
        tempHexRes(2,1) = hexres
        tempHexResB(2,1) = hexborn
      endif

      if(ldebug) then
       gaugeprecision=max(gaugeprecision,abs((tempHexResG(1,1)+tempHexResG(2,1))/max(abs(tempHexResG(1,1)),abs(tempHexResG(2,1)))))
       gaugeprecision=max(gaugeprecision,abs((tempHexResG(1,2)+tempHexResG(2,2))/max(abs(tempHexResG(1,2)),abs(tempHexResG(2,2)))))
       call check('gtest HexB1',-tempHexResG(1,1)/tempHexResG(2,1))
      endif
      
      wardidtest(1) = abs((tempHexResG(1,1)+tempHexResG(2,1))/tempHexResG(1,1)).lt.GaugeAccuracy
      wardidtest(2) = abs((tempHexResG(1,2)+tempHexResG(2,2))/tempHexResG(1,2)).lt.GaugeAccuracy
!       if(abs(tempHexResG(1,1)).lt.1d-6.or.abs(tempHexResG(1,1)).lt.1d-6) print*,__LINE__,tempHexResG
      wardidtest(0) = wardidtest(1).and.wardidtest(2).or.ldebug
      do q=3,4
      do h=-1,1,2
      if(lhel(h).and.(cross.eq.0.or.qtype(3).eq.q)) then
        ampb(col(1),q,h) = ampb(col(1),q,h) + (tempHexResB(1,h)+tempHexResB(2,h))*GG(1)**2*crosssign
        ampv(col(1),q,h) = ampv(col(1),q,h) + (tempHexRes(1,h)+tempHexRes(2,h))*GG(1)**2*crosssign*4d0/3d0
        ampv(col(2),q,h) = ampv(col(2),q,h) + (tempHexRes(1,h)+tempHexRes(2,h))*GG(1)**2*crosssign/2d0
      endif
      enddo;enddo
      if(.not.wardidtest(0)) then
        gaugeOK(1) = .false.
      endif ! wardidtest


c in w  out ; in2 v out2
      if(ldebug) then
        print*
        print*, 'hex2 Dire, in,out',in,out
      endif
      
c start with gtest
      if (gtesthel.eq.-1) then
        call HexBox277T(pdia(0,out),pdia(0,out2),kdia(0,2),pdia(0,in2),pdia(0,in),kdia(0,1),
     &               psi4(1,out),psi4(1,in),psi4(1,out2),psi4(1,in2),kdiac(0,2),epsW(0,1),
     &               musq,1,0,0,hexamp,hexres,hexborn,div)
      else
        call HexBox276T(pdia(0,out),pdia(0,out2),kdia(0,2),pdia(0,in2),pdia(0,in),kdia(0,1),
     &               psi4(1,out),psi4(1,in),psi4(1,out2),psi4(1,in2),kdiac(0,2),epsW(0,1),
     &               musq,1,0,0,hexamp,hexres,hexborn,div)
      endif

      if(ldebug) then
        print*, "Gauge test HexBox"
        call check('r127',(hexamp(1)+hexamp(2))/(-hexamp(7)))
        call check('r348',(hexamp(3)+hexamp(4))/(-hexamp(8)))
        call check('r569',(hexamp(5)+hexamp(6))/(-hexamp(9)))
        if(ldebug) gaugeprecision = max(gaugeprecision,abs(1d0-(hexamp(1)+hexamp(2))/(-hexamp(7))))
        if(ldebug) gaugeprecision = max(gaugeprecision,abs(1d0-(hexamp(3)+hexamp(4))/(-hexamp(8))))
        if(ldebug) gaugeprecision = max(gaugeprecision,abs(1d0-(hexamp(5)+hexamp(6))/(-hexamp(9))))  
      endif
	
      wardidtest(1) = abs((hexamp(1)+hexamp(2))/(-hexamp(7))-1d0) .lt. GaugeAccuracy
     &           .and.abs((hexamp(3)+hexamp(4))/(-hexamp(8))-1d0) .lt. GaugeAccuracy
     &           .and.abs((hexamp(5)+hexamp(6))/(-hexamp(9))-1d0) .lt. GaugeAccuracy
      
!       if(min(abs(hexamp(7)),abs(hexamp(8)),abs(hexamp(9))).lt.1d-7) print*,__LINE__,hexamp,'\n'
      if (gtesthel.eq.-1) then
        call HexBox277T(pdia(0,out),pdia(0,out2),kdia(0,2),pdia(0,in2),pdia(0,in),kdia(0,1),
     &               psi4(1,out),psi4(1,in),psi4(1,out2),psi4(1,in2),epsV(0,lephel,gtestq,-1,1),kdiac(0,1),
     &               musq,-1,0,0,hexamp,hexres,hexborn,div)
      else
        call HexBox276T(pdia(0,out),pdia(0,out2),kdia(0,2),pdia(0,in2),pdia(0,in),kdia(0,1),
     &               psi4(1,out),psi4(1,in),psi4(1,out2),psi4(1,in2),epsV(0,lephel,gtestq,1,1),kdiac(0,1),
     &               musq,-1,0,0,hexamp,hexres,hexborn,div)
      endif

      if(ldebug) then
        print*, "Gauge test HexBox"
        call check('r135',(hexamp(1)+hexamp(3))/(-hexamp(5)))
        call check('r246',(hexamp(2)+hexamp(4))/(-hexamp(6)))
        call check('r789',(hexamp(7)+hexamp(8))/(-hexamp(9)))
        if(ldebug) gaugeprecision = max(gaugeprecision,abs(1d0-(hexamp(1)+hexamp(3))/(-hexamp(5))))
        if(ldebug) gaugeprecision = max(gaugeprecision,abs(1d0-(hexamp(2)+hexamp(4))/(-hexamp(6))))
        if(ldebug) gaugeprecision = max(gaugeprecision,abs(1d0-(hexamp(7)+hexamp(8))/(-hexamp(9))))
      endif

      wardidtest(2) = abs((hexamp(1)+hexamp(3))/(-hexamp(5))-1d0) .lt. GaugeAccuracy
     &           .and.abs((hexamp(2)+hexamp(4))/(-hexamp(6))-1d0) .lt. GaugeAccuracy
     &           .and.abs((hexamp(7)+hexamp(8))/(-hexamp(9))-1d0) .lt. GaugeAccuracy
      
!       if(min(abs(hexamp(5)),abs(hexamp(6)),abs(hexamp(9))).lt.1d-7) print*,__LINE__,hexamp,'\n'
      wardidtest(0) = wardidtest(1).and.wardidtest(2).or.ldebug.or.div.ne.0
      gaugeOK(1) = gaugeOK(1) .and. wardidtest(0)


! c wardtest ok? calc amplitude
!       if(wardidtest(0)) then

c qflav2 = d, hel2=-1
      if(lhel(-1).and.(cross.eq.0.or.qtype(3).eq.3)) then
        call HexBox277T(pdia(0,out),pdia(0,out2),kdia(0,2),pdia(0,in2),pdia(0,in),kdia(0,1),
     &               psi4(1,out),psi4(1,in),psi4(1,out2),psi4(1,in2),epsV(0,lephel,3,-1,1),epsW(0,1),
     &               musq,-1,0,0,hexamp,hexres,hexborn,div)
        ampb(col(1),3,-1) = ampb(col(1),3,-1) + hexborn*GG(1)**2*crosssign
        ampv(col(1),3,-1) = ampv(col(1),3,-1) - hexres*crosssign*G**2/6d0 ! convert to born col str
        ampv(col(2),3,-1) = ampv(col(2),3,-1) + hexres*crosssign*G**2/2d0
      endif

c qflav2 = u, hel2=-1
      if(lhel(-1).and.(cross.eq.0.or.qtype(3).eq.4)) then
        call HexBox277T(pdia(0,out),pdia(0,out2),kdia(0,2),pdia(0,in2),pdia(0,in),kdia(0,1),
     &               psi4(1,out),psi4(1,in),psi4(1,out2),psi4(1,in2),epsV(0,lephel,4,-1,1),epsW(0,1),
     &               musq,-1,0,0,hexamp,hexres,hexborn,div)
        ampb(col(1),4,-1) = ampb(col(1),4,-1) + hexborn*GG(1)**2*crosssign
        ampv(col(1),4,-1) = ampv(col(1),4,-1) - hexres*crosssign*G**2/6d0 ! convert to born col str
        ampv(col(2),4,-1) = ampv(col(2),4,-1) + hexres*crosssign*G**2/2d0
      endif

c qflav2 = d, hel2=+1
      if(lhel(-1).and.(cross.eq.0.or.qtype(3).eq.3)) then
        call HexBox276T(pdia(0,out),pdia(0,out2),kdia(0,2),pdia(0,in2),pdia(0,in),kdia(0,1),
     &               psi4(1,out),psi4(1,in),psi4(1,out2),psi4(1,in2),epsV(0,lephel,3,1,1),epsW(0,1),
     &               musq,-1,0,0,hexamp,hexres,hexborn,div)
        ampb(col(1),3,1) = ampb(col(1),3,1) + hexborn*GG(1)**2*crosssign
        ampv(col(1),3,1) = ampv(col(1),3,1) - hexres*crosssign*G**2/6d0 ! convert to born col str
        ampv(col(2),3,1) = ampv(col(2),3,1) + hexres*crosssign*G**2/2d0
      endif

c qflav2 = u, hel2=1
      if(lhel(-1).and.(cross.eq.0.or.qtype(3).eq.4)) then
        call HexBox276T(pdia(0,out),pdia(0,out2),kdia(0,2),pdia(0,in2),pdia(0,in),kdia(0,1),
     &               psi4(1,out),psi4(1,in),psi4(1,out2),psi4(1,in2),epsV(0,lephel,4,1,1),epsW(0,1),
     &               musq,-1,0,0,hexamp,hexres,hexborn,div)
        ampb(col(1),4,1) = ampb(col(1),4,1) + hexborn*GG(1)**2*crosssign
        ampv(col(1),4,1) = ampv(col(1),4,1) - hexres*crosssign*G**2/6d0 ! convert to born col str
        ampv(col(2),4,1) = ampv(col(2),4,1) + hexres*crosssign*G**2/2d0
      endif
      
!       endif ! wardidtest


c out --<---$$$$$---<-- in2
c          |     |
c          |     |
c in  -->---$$$$$--->-- out2

c in mu w1 nu out ; in2 nu w2 mu out2
      if(ldebug) then
        print*
        print*,'hex2 Cross'
      endif

c start with gtest
      if (gtesthel.eq.-1) then
        call HexBoxCross277T(pdia(0,in),pdia(0,out2),kdia(0,2),pdia(0,in2),pdia(0,out),kdia(0,1),
     &               psi4(1,in),psi4(1,out),psi4(1,out2),psi4(1,in2),kdiac(0,2),epsW(0,1),
     &               musq,1,0,0,hexamp,hexres,hexborn,div)
      else
        call HexBoxCross276T(pdia(0,in),pdia(0,out2),kdia(0,2),pdia(0,in2),pdia(0,out),kdia(0,1),
     &               psi4(1,in),psi4(1,out),psi4(1,out2),psi4(1,in2),kdiac(0,2),epsW(0,1),
     &               musq,1,0,0,hexamp,hexres,hexborn,div)
      endif

      if(ldebug) then
        print*, "Gauge test HexBox"
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

!       if(min(abs(hexamp(7)),abs(hexamp(8)),abs(hexamp(9))).lt.1d-7) print*,__LINE__,hexamp,'\n'
      if(gtesthel.eq.-1) then
        call HexBoxCross277T(pdia(0,in),pdia(0,out2),kdia(0,2),pdia(0,in2),pdia(0,out),kdia(0,1),
     &               psi4(1,in),psi4(1,out),psi4(1,out2),psi4(1,in2),epsV(0,lephel,gtestq,gtesthel,1),kdiac(0,1),
     &               musq,-1,0,0,hexamp,hexres,hexborn,div)
      else
        call HexBoxCross276T(pdia(0,in),pdia(0,out2),kdia(0,2),pdia(0,in2),pdia(0,out),kdia(0,1),
     &               psi4(1,in),psi4(1,out),psi4(1,out2),psi4(1,in2),epsV(0,lephel,gtestq,gtesthel,1),kdiac(0,1),
     &               musq,-1,0,0,hexamp,hexres,hexborn,div)
      endif

      if(ldebug) then
        print*, "Gauge test HexBox"
        call check('r135',(hexamp(1)+hexamp(3))/(-hexamp(5)))
        call check('r246',(hexamp(2)+hexamp(4))/(-hexamp(6)))
        call check('r789',(hexamp(7)+hexamp(8))/(-hexamp(9)))
        if(ldebug) gaugeprecision = max(gaugeprecision,abs(1d0-(hexamp(1)+hexamp(3))/(-hexamp(5))))
        if(ldebug) gaugeprecision = max(gaugeprecision,abs(1d0-(hexamp(2)+hexamp(4))/(-hexamp(6))))
        if(ldebug) gaugeprecision = max(gaugeprecision,abs(1d0-(hexamp(7)+hexamp(8))/(-hexamp(9))))
      endif

      wardidtest(2) = abs((hexamp(1)+hexamp(3))/(-hexamp(5))-1d0) .lt. GaugeAccuracy
     &           .and.abs((hexamp(2)+hexamp(4))/(-hexamp(6))-1d0) .lt. GaugeAccuracy
     &           .and.abs((hexamp(7)+hexamp(8))/(-hexamp(9))-1d0) .lt. GaugeAccuracy

!       if(min(abs(hexamp(5)),abs(hexamp(6)),abs(hexamp(9))).lt.1d-7) print*,__LINE__,hexamp,'\n'
      wardidtest(0) = wardidtest(1).and.wardidtest(2).or.ldebug.or.div.ne.0
      gaugeOK(1) = gaugeOK(1) .and. wardidtest(0)


! c wardtest ok? calc amplitude
!       if(wardidtest(0)) then

c qflav2 = d, hel2=-1
      if(lhel(-1).and.(cross.eq.0.or.qtype(3).eq.3)) then
        call HexBoxCross277T(pdia(0,in),pdia(0,out2),kdia(0,2),pdia(0,in2),pdia(0,out),kdia(0,1),
     &               psi4(1,in),psi4(1,out),psi4(1,out2),psi4(1,in2),epsV(0,lephel,3,-1,1),epsW(0,1),
     &               musq,comp,0,0,hexamp,hexres,hexborn,div)
        ampb(col(1),3,-1) = ampb(col(1),3,-1) + hexborn*GG(1)**2*crosssign
        ampv(col(1),3,-1) = ampv(col(1),3,-1) + hexres*crosssign*G**2*4d0/3d0 ! convert to born col str
        ampv(col(2),3,-1) = ampv(col(2),3,-1) + hexres*crosssign*G**2/2d0
      endif

c qflav2 = u, hel2=-1
      if(lhel(-1).and.(cross.eq.0.or.qtype(3).eq.4)) then
        call HexBoxCross277T(pdia(0,in),pdia(0,out2),kdia(0,2),pdia(0,in2),pdia(0,out),kdia(0,1),
     &               psi4(1,in),psi4(1,out),psi4(1,out2),psi4(1,in2),epsV(0,lephel,4,-1,1),epsW(0,1),
     &               musq,comp,0,0,hexamp,hexres,hexborn,div)
        ampb(col(1),4,-1) = ampb(col(1),4,-1) + hexborn*GG(1)**2*crosssign
        ampv(col(1),4,-1) = ampv(col(1),4,-1) + hexres*crosssign*G**2*4d0/3d0 ! convert to born col str
        ampv(col(2),4,-1) = ampv(col(2),4,-1) + hexres*crosssign*G**2/2d0
      endif

c qflav2 = d, hel2=+1
      if(lhel(-1).and.(cross.eq.0.or.qtype(3).eq.3)) then
        call HexBoxCross276T(pdia(0,in),pdia(0,out2),kdia(0,2),pdia(0,in2),pdia(0,out),kdia(0,1),
     &               psi4(1,in),psi4(1,out),psi4(1,out2),psi4(1,in2),epsV(0,lephel,3,1,1),epsW(0,1),
     &               musq,comp,0,0,hexamp,hexres,hexborn,div)
        ampb(col(1),3,1) = ampb(col(1),3,1) + hexborn*GG(1)**2*crosssign
        ampv(col(1),3,1) = ampv(col(1),3,1) + hexres*crosssign*G**2*4d0/3d0 ! convert to born col str
        ampv(col(2),3,1) = ampv(col(2),3,1) + hexres*crosssign*G**2/2d0
      endif

c qflav2 = u, hel2=1
      if(lhel(-1).and.(cross.eq.0.or.qtype(3).eq.4)) then
        call HexBoxCross276T(pdia(0,in),pdia(0,out2),kdia(0,2),pdia(0,in2),pdia(0,out),kdia(0,1),
     &               psi4(1,in),psi4(1,out),psi4(1,out2),psi4(1,in2),epsV(0,lephel,4,1,1),epsW(0,1),
     &               musq,comp,0,0,hexamp,hexres,hexborn,div)
        ampb(col(1),4,1) = ampb(col(1),4,1) + hexborn*GG(1)**2*crosssign
        ampv(col(1),4,1) = ampv(col(1),4,1) + hexres*crosssign*G**2*4d0/3d0 ! convert to born col str
        ampv(col(2),4,1) = ampv(col(2),4,1) + hexres*crosssign*G**2/2d0
      endif

!       endif ! wardidtest

      hexcount2m = hexcount2m + 1
      hexcount2sub(sub_number) = hexcount2sub(sub_number) + 1
      if(.not.gaugeOK(1)) then
        hexcountm = hexcountm + 1
        hexcountsub(sub_number) = hexcountsub(sub_number) + 1
        hexcountm_Q = hexcountm_q + 1
        hexcountsub_Q(sub_number) = hexcountsub_Q(sub_number) + 1
        if(return0ifGaugeFails) then
          do i=0,2; do j=3,4; do k=0,2
            res(i,j,k) = 0d0
          enddo; enddo; enddo
          return
        else
          recalculatePSpoint=.true.
        endif
      endif

      endif ! lhex


      if(ldebug) then
	print*,'same born amp?'
	print*,((ampb(i,j,-1)/amp(i,j,-1),i=1,2),j=3,4)
	print*,((ampb(i,j,1)/amp(i,j,1),i=1,2),j=3,4)
      endif

c 
      do q=3,4
      do h=-1,1,2
	squared(1,1) = 2d0*dble(ampv(1,q,h)*dconjg(amp(1,q,h)))
	squared(2,1) = 2d0*dble(ampv(2,q,h)*dconjg(amp(1,q,h)))
	squared(2,2) = 2d0*dble(ampv(2,q,h)*dconjg(amp(2,q,h)))
	squared(1,2) = 2d0*dble(ampv(1,q,h)*dconjg(amp(2,q,h)))
	
c this is the result for different flavours |----|**2 + |--++|**2
	if(cross.eq.0) then
	  tempR = (2d0*squared(1,1)-2d0/3d0*squared(2,1))*als(1,1)/(4d0* Pi) ! alphas/4pi
	  resv(1,q) = resv(1,q)+tempR
	else
	  if(qtype(3).eq.q) then
	    tempR = (2d0*squared(2,2)-2d0/3d0*squared(1,2))*als(1,1)/(4d0* Pi) ! alphas/4pi
	    resv(2,q) = resv(2,q)+tempR
	  endif
	endif
	
c different flavours
	if(h.eq.1) then 
c                 cross.eq.0: |--++|**2 = MV_i*T_i1*MB_1
c                 cross.eq.1: |-+-+|**2 = MV_i*T_i2*MB_2
	    tempR = (2d0*(squared(col(1),col(1)))-2d0/3d0*(squared(col(2),col(1))))*als(1,1)/(4d0* Pi) ! alphas/4pi
	    resv(0,q) = resv(0,q)+tempR

	else
c here we have to add up both crossing, before we can square the amplitude
	  if(cross.eq.0) then
	    ampvcross1(1,q,h) = ampv(1,q,h) ! store this contribution
	    ampvcross1(2,q,h) = ampv(2,q,h)
	  else
	    if(q.eq.qtype(3)) then
	      ampv(1,q,h) = ampv(1,q,h) + ampvcross1(1,q,h)
	      ampv(2,q,h) = ampv(2,q,h) + ampvcross1(2,q,h)
	      
	    
	      squared(1,1) = 2d0*dble(ampv(1,q,h)*dconjg(amp(1,q,h)))
	      squared(2,1) = 2d0*dble(ampv(2,q,h)*dconjg(amp(1,q,h)))
	      squared(2,2) = 2d0*dble(ampv(2,q,h)*dconjg(amp(2,q,h)))
	      squared(1,2) = 2d0*dble(ampv(1,q,h)*dconjg(amp(2,q,h)))
	    
	      tempR = (2d0*(squared(1,1)+squared(2,2))-2d0/3d0*(squared(1,2)+squared(2,1)))*als(1,1)/(4d0* Pi) ! alphas/4pi
	      resv(0,q) = resv(0,q)+tempR
	    endif
	  endif
	endif
	ampv(1,q,h) = 0d0 ! reset amp for cross==1 
	ampv(2,q,h) = 0d0
      enddo
      enddo 
	      
      
 40   continue      
      enddo ! out
      enddo ! in
      
      if(ldebug) print*,"gaugeprecision ",gaugeprecision
      
      if(div.ne.0) then ! check poles
	print*
	print*,"poles div = ", div
	print*,'id',id
	print*,'qhel',qhel
	print*
	print*,"resb           ",resb
	print*,"res I Operator ",resI
	print*,"res virt+CT    ",resV
	print*,'sum I+virt+CT  ',((resV(i,j)+resI(i,j),i=0,2),j=3,4)
	print*,'cancellation   ',((abs((resI(i,j)+resV(i,j))/resI(i,j)),i=0,2),j=3,4)
      endif

      if(ldebug) then
	print*,"res Born: ",resB
	print*,"res Virt: ",((resV(i,j)+resI(i,j),i=0,2),j=3,4)
! 	pause
      endif
      
c add all cotributions
c same family
      res(0,3,0) = resB(0,3)+resV(0,3)+resI(0,3)
      res(0,4,0) = resB(0,4)+resV(0,4)+resI(0,4)
      res(1,3,0) = res(0,3,0)
      res(1,4,0) = res(0,4,0)
      res(2,3,0) = 0d0
      res(2,4,0) = 0d0
c different family
      res(1,4,1) = resB(1,4)+resV(1,4)+resI(1,4)
      res(2,4,2) = resB(2,4)+resV(2,4)+resI(2,4)
      res(1,3,1) = resB(1,3)+resV(1,3)+resI(1,3)
      res(2,3,2) = resB(2,3)+resV(2,3)+resI(2,3)
      res(0,4,1) = res(1,4,1)
      res(0,4,2) = res(2,4,2)
      res(0,3,1) = res(1,3,1)
      res(0,3,2) = res(2,3,2)
!       pause
      
      if(recalculatePSpoint) then
        recalculatePSpoint=.false.
        res=0d0
        resV=0d0
      endif
      
      if (ldoblha) then
        if (blha_amptype(blha_curproc).eq.1) then ! loop
          blha_amp(3)= blha_amp(3)+
     &      resV(mod(blha_idsubproc(blha_cursubproc,blha_curproc),10),
     &           mod(blha_idsubproc(blha_cursubproc,blha_curproc),100)/10)
          blha_amp(4)= blha_amp(4)+
     &      resB(mod(blha_idsubproc(blha_cursubproc,blha_curproc),10),
     &           mod(blha_idsubproc(blha_cursubproc,blha_curproc),100)/10)
        else
          call BLHA_amptypeerror(blha_amptype(blha_curproc),__FILE__,__LINE__)
        endif
        return
      endif 

      if(IsNan(res,0d0,18,0)) then
        if(ldebug.or..not.recalculatePSpoint) print*,"NaN in "__FILE__,__LINE__
        do i=0,2; do q=3,4; do j=0,2
          res(i,q,j) = 0d0
        enddo; enddo; enddo
      endif


      end       
