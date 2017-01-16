c ************************************************************************************
c Author: Matthias Kerner
c Date: 2012/13
c ************************************************************************************


c 4q W Amplitude
c    q1 >>> #### >>> q3 
c	    ####          , WZ
c	    ####
c    q2 >>> #### >>> q4  
c   qhel: el. (1,2,3) correspond to helicity configurations
c		1: all -
c		2: h1=h3=-1 , h2=h4= + 
c		3: h1!= h3 (and h2, h4 specified by charge flow of the specific subprocess)
c	  qhel= 0: sum hel. confg. 1-3
c   res(colstr,2nd qtype =3,4  ;  same qflavor=0,1)
c      colstr=0..2: for catani seymour: <|>,  <|T1T2|>,  <|T1T3|>
c		    for lo: col summed, and the 2 different color flows
      subroutine amp4qwCS(P,kin,id,fsign,nbos,qhel,res)

      use globalvars, only: ldoblha

      implicit none
#include "common.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"
      real*8 P(0:3,max_p,max_kin),res(0:2,3:4,0:2)
      integer kin, id(4),fsign(4),nbos,qhel
      logical lnlo

      complex*16 epsW(0:3,7),epsWV(0:3,-1:1,7) 
      complex*16 epsV(0:3,-1:1,3:4,-1:1,7) ! ( , mu-hel, u/d , quark-hel), coupl. const to quarks included
      complex*16 epsA(0:3,-1:1),epsZ(0:3,-1:1)
      real*8 momV(0:4,7),momW(0:4,7),momWV(0:4,7)
      common /qcdwznew/ epsW, epsWV, epsV,epsA,epsZ,momV,momW,momWV

      
      real*8 pbar(0:3,4),pdia(0:3,4)
      real*8 mompsiW(0:4,4)
      complex*16 psi(2,-1:1,4), psiW(2,4)
      complex*16 amp(2,3:4,-1:1) ! colstr, 2nd, hel
      complex*16 gluWij(0:3,1:3,2:4),gluij(0:3,-1:1,1:3,2:4)
      complex*16 tempCurr(0:3,-1:1,4)
      complex*16 dotcc
      real*8 qjj2
      real*8 squared(2,2), tempR
      
      integer i,j,mu, tempi4(4), lepkin,colstr,qtype(2)
      integer  ikin, in,out,in2,out2,helcfg

      lnlo = .true.
      ikin = kin
      goto 11
      
      entry amp4qwBorn(P,id,fsign,nbos,qhel,res)
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
      
      qtype(1) = nbos
      qtype(2) = 7-qtype(1)
      
c get bra and ket 
      tempi4(1) = fsign(1)
      tempi4(2) = fsign(3)
      tempi4(3) = fsign(2)
      tempi4(4) = fsign(4)
      call psi0m(4,pbar,tempi4,psi)
      
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

c new convention
c in   >W>>>>W>  out
c          S         , WZ
c	   S
c in2  >>>>>>>>  out2
c get <qW| and |Wq>
      do in = 1,3,2
	  call ket2c(psi(1,-1,in),.true.,pdia(0,in),-1,momW(0,lepkin),epsW(0,lepkin),psiW(1,in),mompsiW(0,in))
	  out = in+1
	  call bra2c(psi(1,-1,out),.true.,pdia(0,out),-1,momW(0,lepkin),epsW(0,lepkin),psiW(1,out),mompsiW(0,out))
      enddo

      do in = 1,3,2  ! determine q line coupling to W
      do out= 2,4,2

      colstr=(in+out)/2
!       if(colstr.ne.3) then
      if(   colstr.eq.1 .and.qhel.le.2
     $	.or.colstr.eq.2 .and.qhel.ne.2) then
	
	in2 = 4-in  ! 2nd q index
	out2 =6-out

c calculate glu currents
	call curr(-1,psi(1,-1,out),psiW(1,in),tempCurr(0,-1,1))
	call curr(-1,psiW(1,out),psi(1,-1,in),tempCurr(0,-1,2))
	do mu=0,3
	  gluWij(mu,in,out) =  tempCurr(mu,-1,1)+tempCurr(mu,-1,2)
	enddo
	call curr(1,psi(1,-1,out2),psi(1,-1,in2),gluij(0,-1,in2,out2))

      
c eventually calc amplitudes
	tempR = dble(GG(1)**2/qjj2(pdia(0,in2),pdia(0,out2)))
	i = out/2 ! swapped initial or final state? only needed for colstr 2
	amp(colstr,qtype(i),-1) = tempR * dotcc(gluWij(0,in,out),gluij(0,-1,in2,out2))
	amp(colstr,qtype(i), 1) = tempR * dotcc(gluWij(0,in,out),gluij(0, 1,in2,out2))
	if(colstr.eq.1) then
	  amp(1,7-qtype(i),-1) = amp(1,qtype(i),-1)
	  amp(1,7-qtype(i), 1) = amp(1,qtype(i), 1)
	endif
	
      endif
      enddo; enddo ! in and out

c square amplitudes
      do i = 3,4 ! 2nd quark
	do j=0,2
	  res(j,i,0) = 0d0
	  res(j,i,1) = 0d0
	  res(j,i,2) = 0d0
	enddo
	
	do helcfg=-1,1,2
! artificially get back the qhel flag
	 if(  helcfg.eq.-1.and.qhel.le.1
     $	  .or.helcfg.eq.1.and.qhel.ne.1) then
!           if(qhel.le.2) amp(4-qhel,i,1) = 0d0
! end qhel

	  amp(2,i,helcfg) = - amp(2,i,helcfg) ! change of fermion indices
	  
	  squared(1,1) = dble(amp(1,i,helcfg)*dconjg(amp(1,i,helcfg)))
	  squared(2,2) = dble(amp(2,i,helcfg)*dconjg(amp(2,i,helcfg)))
	  squared(1,2) = dble(amp(1,i,helcfg)*dconjg(amp(2,i,helcfg)))
	  if(helcfg.eq.1) squared(1,2) = 0d0

	  if(lnlo) then ! CS color factors
            if (ldoblha) then
              res(0,i,0) = res(0,i,0) + blha_CF*blha_CA/2d0*(squared(1,1)+squared(2,2))
     &                                - blha_CF*blha_CA*(blha_CA-2d0*blha_CF)*squared(1,2)
              res(0,i,1) = res(0,i,1) + blha_CF*blha_CA/2d0*(squared(1,1))
              res(0,i,2) = res(0,i,2) + blha_CF*blha_CA/2d0*(squared(2,2))
              res(1,i,0) = res(1,i,0) - blha_CF*blha_CA/2d0*(blha_CA-2d0*blha_CF)
     &                                  *(squared(1,1)+squared(2,2)-2d0*(blha_CA-blha_CF)*squared(1,2))
              res(1,i,1) = res(1,i,1) - blha_CF*blha_CA/2d0*(blha_CA-2d0*blha_CF)*squared(1,1)
              res(1,i,2) = res(1,i,2) - blha_CF*blha_CA/2d0*(blha_CA-2d0*blha_CF)*squared(2,2)
              res(2,i,0) = res(2,i,0) + blha_CF*blha_CA/4d0*(blha_CA-2d0*blha_CF)*squared(1,1)
     &                                + blha_CF*blha_CA/4d0*(blha_CA-4d0*blha_CF)*squared(2,2)
     &                                - blha_CF/2d0*(blha_CA-2d0*blha_CF)*squared(1,2)
              res(2,i,1) = res(2,i,1) + blha_CF*blha_CA/4d0*(blha_CA-2d0*blha_CF)*squared(1,1)
              res(2,i,2) = res(2,i,2) + blha_CF*blha_CA/4d0*(blha_CA-4d0*blha_CF)*squared(2,2)
            else 
	      res(0,i,0) = res(0,i,0) + 2d0*(squared(1,1)+squared(2,2)-2d0/3d0*squared(1,2))
	      res(0,i,1) = res(0,i,1) + 2d0*(squared(1,1))
	      res(0,i,2) = res(0,i,2) + 2d0*(squared(2,2))
	      res(1,i,0) = res(1,i,0) - 2d0/3d0*(squared(1,1)+squared(2,2)-10d0/3d0*squared(1,2))
	      res(1,i,1) = res(1,i,1) - 2d0/3d0*squared(1,1)
	      res(1,i,2) = res(1,i,2) - 2d0/3d0*squared(2,2)
	      res(2,i,0) = res(2,i,0) + 1d0/3d0*(squared(1,1)-7d0*squared(2,2)-2d0/3d0*squared(1,2))
	      res(2,i,1) = res(2,i,1) + 1d0/3d0*squared(1,1)
	      res(2,i,2) = res(2,i,2) - 7d0/3d0*squared(2,2)
            endif
	  else ! get colorstr weighted amp for lha output

	    res(1,i,1) = res(1,i,1) + 2d0*squared(1,1)
	    res(2,i,2) = res(2,i,2) + 2d0*squared(2,2)
	    res(0,i,0) = res(0,i,0) + 2d0*(squared(1,1)+squared(2,2)-2d0/3d0*squared(1,2))

	    res(1,i,0) = res(1,i,0)+squared(1,1)*2d0*(1d0-2d0/3d0*squared(1,2)/(squared(1,1)+squared(2,2)))
	    res(2,i,0) = res(0,i,0) - res(1,i,0)
	  endif
	 endif
	enddo
! ! ! ! ! ! ! ! ! ! ! ! ! 	if(lnlo) then 
! ! ! ! ! ! ! ! ! ! ! ! ! c 	    we swapped particle 1 and 2 --> res(::,2) is born, T21 and T23, --> we have to modify (2,:,2) again
! ! ! ! ! ! ! ! ! ! ! ! ! c           or                  3 and 4 --> res(::,2) is born, T12 and T14
! ! ! ! ! ! ! ! ! ! ! ! ! 	    res(2,i,2) = -4d0/3d0*res(0,i,2)-res(1,i,2)-res(2,i,2)
! ! ! ! ! ! ! ! ! ! ! ! ! 	endif
      enddo

      if (ldoblha) then
        if (blha_amptype(blha_curproc).eq.0) then ! tree
          blha_amp(1)=blha_amp(1)+
     &      res(0,
     &          mod(blha_idsubproc(blha_cursubproc,blha_curproc),100)/10,
     &          mod(blha_idsubproc(blha_cursubproc,blha_curproc),10)    )

        else if (blha_amptype(blha_curproc).eq.2) then ! cctree
          ! 1-2
          call BLHA_cctree(1,2,
     &      + res(1,
     &            mod(blha_idsubproc(blha_cursubproc,blha_curproc),100)/10,
     &            mod(blha_idsubproc(blha_cursubproc,blha_curproc),10)    )
     &    )
          ! 1-3
          call BLHA_cctree(1,3,
     &      + res(2,
     &            mod(blha_idsubproc(blha_cursubproc,blha_curproc),100)/10,
     &            mod(blha_idsubproc(blha_cursubproc,blha_curproc),10)    )
     &    )
          ! 1-4 = -(CF<> + <12> + <13>)
          call BLHA_cctree(1,4,
     &      - ( blha_CF*res(0,
     &                      mod(blha_idsubproc(blha_cursubproc,blha_curproc),100)/10,
     &                      mod(blha_idsubproc(blha_cursubproc,blha_curproc),10))
     &          + res(1,
     &                mod(blha_idsubproc(blha_cursubproc,blha_curproc),100)/10,
     &                mod(blha_idsubproc(blha_cursubproc,blha_curproc),10))
     &          + res(2,
     &                mod(blha_idsubproc(blha_cursubproc,blha_curproc),100)/10,
     &                mod(blha_idsubproc(blha_cursubproc,blha_curproc),10))   )
     &    )
          ! 2-3 = -(CF<> + <12> + <13>)
          call BLHA_cctree(2,3,
     &      - ( blha_CF*res(0,
     &                      mod(blha_idsubproc(blha_cursubproc,blha_curproc),100)/10,
     &                      mod(blha_idsubproc(blha_cursubproc,blha_curproc),10))
     &          + res(1,
     &                mod(blha_idsubproc(blha_cursubproc,blha_curproc),100)/10,
     &                mod(blha_idsubproc(blha_cursubproc,blha_curproc),10))
     &          + res(2,
     &                mod(blha_idsubproc(blha_cursubproc,blha_curproc),100)/10,
     &                mod(blha_idsubproc(blha_cursubproc,blha_curproc),10))   )
     &    )
          ! 2-4 = <13>
          call BLHA_cctree(2,4,
     &      + res(2,
     &            mod(blha_idsubproc(blha_cursubproc,blha_curproc),100)/10,
     &            mod(blha_idsubproc(blha_cursubproc,blha_curproc),10)    )
     &    )
          ! 3-4 = <12>
          call BLHA_cctree(3,4,
     &      + res(1,
     &            mod(blha_idsubproc(blha_cursubproc,blha_curproc),100)/10,
     &            mod(blha_idsubproc(blha_cursubproc,blha_curproc),10)    )
     &    )

        else 
          call BLHA_amptypeerror(blha_amptype(blha_curproc),__FILE__,__LINE__)
        endif
      endif

      end


c q1 q3 g2 g4 W Amplitude
c    ug2dgwzCS
c	Amplitude as needed for Catani Seymor with specified helicities (see arXiv:0905.0883)
c	resDiag: real Matrix elements <l|l>, <l|T1T2|l>, <l|T1T3|l>
c	resOffd: complex Matrix elements <-|+>,<-|T1T2|+>,<-|T1T3|+>, where the gluon, that is involved in dipoles has hel. index 0
c       all hel. configurations will summed if ghel_i =0 
c    ug2dgwzBorn
c	Amplitude as needed for Born process
c	res(0:2): amplitude for specified hel. configuration or summed if ghel_i = 0, index 1 and 2 for colstr for LHA output
      SUBROUTINE ug2dgwCS(P,kin,id,fsign,nbos, resDiag) ! sum over glu hel if ghel ==0

      use globalvars, only: ldoblha

      implicit none
#include "common.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"
      real*8 P(0:3,max_p,max_kin), resDiag(0:2,0:2,0:2), res(0:2)
      integer kin, id(4),fsign(4),nbos , ikin
      logical lnlo

      complex*16 epsW(0:3,7),epsWV(0:3,-1:1,7) 
      complex*16 epsV(0:3,-1:1,3:4,-1:1,7) ! ( , mu-hel, u/d , quark-hel), coupl. const to quarks included
      complex*16 epsA(0:3,-1:1),epsZ(0:3,-1:1)
      real*8 momV(0:4,7),momW(0:4,7),momWV(0:4,7)
      common /qcdwznew/ epsW, epsWV, epsV,epsA,epsZ,momV,momW,momWV

      
      real*8 pbar(0:3,2),kbar(0:4,2),pdia(0:3,2),kdia(0:4,3)
      real*8 momPsiW(0:4,2),momPsiG(0:4,2,2)
      complex*16 Psi(2,-1:1,2), PsiW(2,2), PsiG(2,2,2,2) ! (,ferm,hel,glu)
      complex*16 amp(2,0:2,0:2)
      complex*16  scres(0:2,2), amph(2,-1:1,-1:1), im
      real*8 epsG(0:3,2,2), epsGG(0:4)
      complex*16 s1c,s1r
      real*8 mass2,dotrr
      parameter (im = (0d0,1d0))
      
      integer i,mu, pol1,pol2,hel(2), tempi2(2), lepkin,qtype(2)
      integer glu1,glu2
      lnlo = .true.
      ikin = kin
      goto 11
      entry ug2dgwBorn(P,id,fsign,nbos, res)

      ikin = 1
      lnlo= .false.
      res(0) = 0d0
      res(1) = 0d0
      res(2) = 0d0
      
 11   continue
      lepkin = ikin
      if(ikin .gt.7) lepkin = 1

      if(lnlo) then
	resDiag(1,0,0) = 0d0
	resDiag(2,0,0) = 0d0
	resDiag(0,0,0) = 0d0
      endif

      fsign(1)=1
      fsign(2)=1
      fsign(3)=1
      fsign(4)=1
      if(id(2).le.2) fsign(2) = -1
      if(id(4).le.2) fsign(4) = -1
      if(id(1).ge.3) fsign(1) = -1
      if(id(3).le.2) fsign(3) = -1
      do mu=0,3
	pbar(mu,1)=p(mu,id(1),ikin)
	pbar(mu,2)=p(mu,id(3),ikin)
	kbar(mu,1)=p(mu,id(2),ikin)
	kbar(mu,2)=p(mu,id(4),ikin)
	pdia(mu,1)=pbar(mu,1)*fsign(1)
	pdia(mu,2)=pbar(mu,2)*fsign(3)
	kdia(mu,1)=kbar(mu,1)*fsign(2)
	kdia(mu,2)=kbar(mu,2)*fsign(4)
	kdia(mu,3)=kdia(mu,1)+kdia(mu,2)
      enddo      
      kbar(4,1) = 0d0
      kbar(4,2) = 0d0
      kdia(4,1) = 0d0
      kdia(4,2) = 0d0
      kdia(4,3) = mass2(kdia(0,3))
      
      qtype(1) = nbos ! (W+)? u type : d type
      qtype(2) = 7-qtype(1)
      
c get bra and ket 
      tempi2(1) = fsign(1)
      tempi2(2) = fsign(3)
      call psi0m(2,pbar,tempi2,Psi)
c get <q,glu| and |glu,q>
      do pol1 = 1,2
      	call polvec(kbar(0,1),pol1,epsG(0,1,pol1))
        call polvec(kbar(0,2),pol1,epsG(0,2,pol1))
      do i=1,2
        call ket2r(psi(1,-1,1),.true.,pdia(0,1),-1,kdia(0,i), epsG(0,i,pol1),PsiG(1,1,pol1,i),momPsiG(0,1,i))
        call bra2r(psi(1,-1,2),.true.,pdia(0,2),-1,kdia(0,i), epsG(0,i,pol1),PsiG(1,2,pol1,i),momPsiG(0,2,i))
      enddo
      enddo
c get <q,V| and |V,q>
      call ket2c(psi(1,-1,1),.true.,pdia(0,1),-1,momW(0,lepkin), epsW(0,lepkin),PsiW(1,1),momPsiW(0,1))
      call bra2c(psi(1,-1,2),.true.,pdia(0,2),-1,momW(0,lepkin), epsW(0,lepkin),PsiW(1,2),momPsiW(0,2))
!       print*,(PsiV(i,1),i=1,2)
!       print*,(epsG(i,2,gpol2),i=0,3)
      do pol1 = 1,2
      do pol2 = 1,2
        hel(1) = pol1
        hel(2) = pol2
        do mu = 0,3
	  epsGG(mu) = (2d0*dotrr(kdia(0,2),epsG(0,1,pol1))*epsG(mu,2,pol2)
     &		      -2d0*dotrr(kdia(0,1),epsG(0,2,pol2))*epsG(mu,1,pol1)
     &		      +dotrr(epsG(0,1,pol1),epsG(0,2,pol2))*(kdia(mu,1)-kdia(mu,2)))/kdia(4,3)
        enddo
c g(g1,g2) amplitudes
        amp(1,pol1,pol2) = -s1r(psi(1,-1,2),epsGG,.true.,-1,PsiW(1,1))
     &           -s1r(psiW(1,2),epsGG,.true.,-1,Psi(1,-1,1))
	amp(2,pol1,pol2) = -amp(1,pol1,pol2)
! 	print*,s1c(psi(1,-1,2),epsGG,.true.,-1,PsiWV(1,1))
	do glu1=1,2
	  glu2 = 3-glu1
c V-> WV amplitudes
	  amp(glu1,pol1,pol2) = amp(glu1,pol1,pol2)
     &		     + s1r(PsiG(1,2,hel(glu1),glu1),epsG(0,glu2,hel(glu2)),.false.,-1,PsiW(1,1))
     &		     + s1c(PsiG(1,2,hel(glu1),glu1),epsW(0,lepkin),.true.,-1,PsiG(1,1,hel(glu2),glu2))
     &		     + s1r(PsiW(1,2),epsG(0,glu1,hel(glu1)),.false.,-1,PsiG(1,1,hel(glu2),glu2))
	enddo
	
c multiply coupl. const. and get |amp|^2
	amp(1,pol1,pol2) = amp(1,pol1,pol2)*GG(1)**2
	amp(2,pol1,pol2) = amp(2,pol1,pol2)*GG(1)**2
! 	print*,amp(1,pol1,pol2),amp(2,pol1,pol2)
	if (lnlo) then
          if (ldoblha) then
            resDiag(0,pol1,pol2) = dble(
     &        blha_CF**2*blha_CA*                  (amp(1,pol1,pol2)*dconjg(amp(1,pol1,pol2))
     &                                             +amp(2,pol1,pol2)*dconjg(amp(2,pol1,pol2)))
     &       -blha_CF*blha_CA*(blha_CA-2d0*blha_CF)*amp(1,pol1,pol2)*dconjg(amp(2,pol1,pol2)))
            resDiag(1,pol1,pol2) = dble(
     &        blha_CF*blha_CA**2*(blha_CA-2d0*blha_CF)/4d0*amp(1,pol1,pol2)*dconjg(amp(1,pol1,pol2))
     &       -blha_CF**2*blha_CA**2/2d0*                   amp(2,pol1,pol2)*dconjg(amp(2,pol1,pol2))
     &       +blha_CF*blha_CA**2*(blha_CA-2d0*blha_CF)/2d0*amp(1,pol1,pol2)*dconjg(amp(2,pol1,pol2)))
            resDiag(2,pol1,pol2) = dble(
     &       -blha_CF*(blha_CA-2d0*blha_CF)/4d0*                     (amp(1,pol1,pol2)*dconjg(amp(1,pol1,pol2))
     &                                                               +amp(2,pol1,pol2)*dconjg(amp(2,pol1,pol2)))
     &       -blha_CF*blha_CA*(blha_CA-blha_CF)*(blha_CA-2d0*blha_CF)*amp(1,pol1,pol2)*dconjg(amp(2,pol1,pol2)))
          else
  	    resDiag(0,pol1,pol2) = dble(16d0*(amp(1,pol1,pol2)*dconjg(amp(1,pol1,pol2))+amp(2,pol1,pol2)*dconjg(amp(2,pol1,pol2)))
     &                                -4d0*amp(1,pol1,pol2)*dconjg(amp(2,pol1,pol2)))/3d0
            resDiag(1,pol1,pol2) = dble(amp(1,pol1,pol2)*dconjg(amp(1,pol1,pol2))-8d0*amp(2,pol1,pol2)*dconjg(amp(2,pol1,pol2))
     &                                +2d0*amp(1,pol1,pol2)*dconjg(amp(2,pol1,pol2)))
            resDiag(2,pol1,pol2) = dble(-amp(1,pol1,pol2)*dconjg(amp(1,pol1,pol2))-amp(2,pol1,pol2)*dconjg(amp(2,pol1,pol2))
     &                                -20d0*amp(1,pol1,pol2)*dconjg(amp(2,pol1,pol2)))/9d0
          endif
	  resDiag(0,0,0)=resDiag(0,0,0)+resDiag(0,pol1,pol2)
	  resDiag(1,0,0)=resDiag(1,0,0)+resDiag(1,pol1,pol2)
	  resDiag(2,0,0)=resDiag(2,0,0)+resDiag(2,pol1,pol2)
	else
  	  res(0) = res(0)+ dble(16d0*(amp(1,pol1,pol2)*dconjg(amp(1,pol1,pol2))+amp(2,pol1,pol2)*dconjg(amp(2,pol1,pol2)))
     &                              -4d0*amp(1,pol1,pol2)*dconjg(amp(2,pol1,pol2)))/3d0
          res(1) = res(1)+dble(amp(1,pol1,pol2)*dconjg(amp(1,pol1,pol2)))
          res(2) = res(2)+dble(amp(2,pol1,pol2)*dconjg(amp(2,pol1,pol2)))
        endif
      enddo; enddo ! hel summation
      if (lnlo) then
        if (ldoblha) then ! <M,i-|Tij|M,i+>
c convert cartesian to helicity basis
          do glu1=1,2
            amph(glu1,+1,+1) = 
     &        ( +                     amp(glu1,1,1)
     &          -im*fsign(2)         *amp(glu1,2,1)
     &          -im*         fsign(4)*amp(glu1,1,2)
     &          -   fsign(2)*fsign(4)*amp(glu1,2,2) )/2d0
            amph(glu1,-1,+1) = 
     &        ( +                     amp(glu1,1,1)
     &          +im*fsign(2)*         amp(glu1,2,1)
     &          -im*         fsign(4)*amp(glu1,1,2)
     &          +   fsign(2)*fsign(4)*amp(glu1,2,2) )/2d0
            amph(glu1,+1,-1) = 
     &        ( +                     amp(glu1,1,1)
     &          -im*fsign(2)*         amp(glu1,2,1)
     &          +im*         fsign(4)*amp(glu1,1,2)
     &          +   fsign(2)*fsign(4)*amp(glu1,2,2) )/2d0
            amph(glu1,-1,-1) = 
     &        ( +                     amp(glu1,1,1)
     &          +im*fsign(2)*         amp(glu1,2,1)
     &          +im*         fsign(4)*amp(glu1,1,2)
     &          -   fsign(2)*fsign(4)*amp(glu1,2,2) )/2d0
          enddo

          scres(:,:) = 0d0
          do pol1 = -1,1,2
            scres(0,1) = scres(0,1) 
     &        + blha_CF**2*blha_CA*(amph(1,+1,pol1)*dconjg(amph(1,-1,pol1))
     &                             +amph(2,+1,pol1)*dconjg(amph(2,-1,pol1)))
     &        - blha_CF/2d0       *(amph(1,+1,pol1)*dconjg(amph(2,-1,pol1))
     &                             +amph(2,+1,pol1)*dconjg(amph(1,-1,pol1)))
            scres(1,1) = scres(1,1) 
     &        + blha_CF*blha_CA/4d0*       amph(1,+1,pol1)*dconjg(amph(1,-1,pol1))
     &        - blha_CF**2*blha_CA**2/2d0* amph(2,+1,pol1)*dconjg(amph(2,-1,pol1))
     &        + blha_CF*blha_CA/4d0*      (amph(1,+1,pol1)*dconjg(amph(2,-1,pol1))
     &                                    +amph(2,+1,pol1)*dconjg(amph(1,-1,pol1)))
            scres(2,1) = scres(2,1) 
     &        - blha_CF/blha_CA/4d0*          (amph(1,+1,pol1)*dconjg(amph(1,-1,pol1))
     &                                        +amph(2,+1,pol1)*dconjg(amph(2,-1,pol1)))
     &        - blha_CF*(blha_CA-blha_CF)/2d0*(amph(1,+1,pol1)*dconjg(amph(2,-1,pol1))
     &                                        +amph(2,+1,pol1)*dconjg(amph(1,-1,pol1)))
            scres(0,2) = scres(0,2) 
     &        + blha_CF**2*blha_CA*(amph(1,pol1,+1)*dconjg(amph(1,pol1,-1))
     &                             +amph(2,pol1,+1)*dconjg(amph(2,pol1,-1)))
     &        - blha_CF/2d0       *(amph(1,pol1,+1)*dconjg(amph(2,pol1,-1))
     &                             +amph(2,pol1,+1)*dconjg(amph(1,pol1,-1)))
            scres(1,2) = scres(1,2) 
     &        + blha_CF*blha_CA/4d0*       amph(1,pol1,+1)*dconjg(amph(1,pol1,-1))
     &        - blha_CF**2*blha_CA**2/2d0* amph(2,pol1,+1)*dconjg(amph(2,pol1,-1))
     &        + blha_CF*blha_CA/4d0*      (amph(1,pol1,+1)*dconjg(amph(2,pol1,-1))
     &                                    +amph(2,pol1,+1)*dconjg(amph(1,pol1,-1)))
            scres(2,2) = scres(2,2) 
     &        - blha_CF/blha_CA/4d0*          (amph(1,pol1,+1)*dconjg(amph(1,pol1,-1))
     &                                        +amph(2,pol1,+1)*dconjg(amph(2,pol1,-1)))
     &        - blha_CF*(blha_CA-blha_CF)/2d0*(amph(1,pol1,+1)*dconjg(amph(2,pol1,-1))
     &                                        +amph(2,pol1,+1)*dconjg(amph(1,pol1,-1)))
          enddo 
        else
          do pol1 = 1,2
  	    resDiag(0,pol1,0) = 2d0*dble(16d0*(amp(1,pol1,1)*dconjg(amp(1,pol1,2))+amp(2,pol1,1)*dconjg(amp(2,pol1,2)))
     &                                -2d0*amp(1,pol1,1)*dconjg(amp(2,pol1,2))-2d0*amp(2,pol1,1)*dconjg(amp(1,pol1,2)))/3d0
            resDiag(1,pol1,0) = 2d0*dble(amp(1,pol1,1)*dconjg(amp(1,pol1,2))-8d0*amp(2,pol1,1)*dconjg(amp(2,pol1,2))
     &                                +amp(1,pol1,1)*dconjg(amp(2,pol1,2))+amp(2,pol1,1)*dconjg(amp(1,pol1,2)))
            resDiag(2,pol1,0) = 2d0*dble(-amp(1,pol1,1)*dconjg(amp(1,pol1,2))-amp(2,pol1,1)*dconjg(amp(2,pol1,2))
     &                                -10d0*amp(1,pol1,1)*dconjg(amp(2,pol1,2))-10d0*amp(2,pol1,1)*dconjg(amp(1,pol1,2)))/9d0
  	    resDiag(0,0,pol1) = 2d0*dble(16d0*(amp(1,1,pol1)*dconjg(amp(1,2,pol1))+amp(2,1,pol1)*dconjg(amp(2,2,pol1)))
     &                                -2d0*amp(1,1,pol1)*dconjg(amp(2,2,pol1))-2d0*amp(2,1,pol1)*dconjg(amp(1,2,pol1)))/3d0
            resDiag(1,0,pol1) = 2d0*dble(amp(1,1,pol1)*dconjg(amp(1,2,pol1))-8d0*amp(2,1,pol1)*dconjg(amp(2,2,pol1))
     &                                +amp(1,1,pol1)*dconjg(amp(2,2,pol1))+amp(2,1,pol1)*dconjg(amp(1,2,pol1)))
            resDiag(2,0,pol1) = 2d0*dble(-amp(1,1,pol1)*dconjg(amp(1,2,pol1))-amp(2,1,pol1)*dconjg(amp(2,2,pol1))
     &                                -10d0*amp(1,1,pol1)*dconjg(amp(2,2,pol1))-10d0*amp(2,1,pol1)*dconjg(amp(1,2,pol1)))/9d0
          enddo
        endif
      else
        res(1) = res(1)/(res(1)+res(2))*res(0)
	res(2) = res(0)-res(1)
      endif
      if (ldoblha) then
        if (blha_amptype(blha_curproc).eq.0) then ! tree
          blha_amp(1)=blha_amp(1)+
     &      resDiag(0,0,0)
        else if (blha_amptype(blha_curproc).eq.2) then ! cctree
          ! 1-2
          call BLHA_cctree(1,2,
     &      + resDiag(1,0,0)
     &    )
          ! 1-3
          call BLHA_cctree(1,3,
     &      + resDiag(2,0,0)
     &    )
          ! 1-4 = -(CF<> + <12> + <13>)
          call BLHA_cctree(1,4,
     &      - ( blha_CF*resDiag(0,0,0)
     &          + resDiag(1,0,0)
     &          + resDiag(2,0,0) )
     &    )
          ! 2-3 = -(CF<> + <12> + <13>)
          call BLHA_cctree(2,3,
     &      - ( blha_CF*resDiag(0,0,0)
     &          + resDiag(1,0,0)
     &          + resDiag(2,0,0) )
     &    )
          ! 2-4 = (CF-CA)<>+<13>
          call BLHA_cctree(2,4,
     &      ( + (blha_CF-blha_CA)*resDiag(0,0,0)
     &        + resDiag(2,0,0) )
     &    )
          ! 3-4 = <12>
          call BLHA_cctree(3,4,
     &      + resDiag(1,0,0)
     &    )
        else if (blha_amptype(blha_curproc).eq.3) then ! sctree
          ! 2-1 = <2|12|2>
          call BLHA_sctree(2,1,
     &      scres(1,1)
     &    )
          ! 2-2 = CA<2|2>
          call BLHA_sctree(2,2,
     &      blha_CA*scres(0,1)
     &    )
          ! 2-3 = -(CF<2|2> + <2|12|2> + <2|13|2>)
          call BLHA_sctree(2,3,
     &      - ( blha_CF*scres(0,1)
     &          + scres(1,1)
     &          + scres(2,1) )
     &    )
          ! 2-4 = (CF-CA)<2|2> + <2|13|2>
          call BLHA_sctree(2,4,
     &      ( + (blha_CF-blha_CA)*scres(0,1)
     &        + scres(2,1) )
     &    )
          ! 4-1 = -(CF<4|4> + <4|12|4> + <4|13|4>)
          call BLHA_sctree(4,1,
     &      - ( blha_CF*scres(0,2)
     &          + scres(1,2)
     &          + scres(2,2) )
     &    )
          ! 4-2 = (CF-CA)<4|4> + <4|13|4>
          call BLHA_sctree(4,2,
     &      ( + (blha_CF-blha_CA)*scres(0,2)
     &        + scres(2,2) )
     &    )
          ! 4-3 = <4|12|4>
          call BLHA_sctree(4,3,
     &      scres(1,2)
     &    )
          ! 4-4 = CA<4|4>
          call BLHA_sctree(4,4,
     &      blha_CA*scres(0,2)
     &    )

        else 
          call BLHA_amptypeerror(blha_amptype(blha_curproc),__FILE__,__LINE__)
        endif
      endif
      end
      
      
c 4q WZ Amplitude
      subroutine amp4qwVirt(P,id,fsign,nbos,nlo,rn,res)

      use globalvars, only: ldoblha

      implicit none
#include "common.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"
      real*8 P(0:3,max_p,max_kin),rn,res(0:2,3:4,0:2)
      integer id(4),fsign(4),nbos,qhel,nlo
!       logical skipidfl(2)
      
      real*8 pbar(0:3,4),pdia(0:3,4),kdia(0:3,1)
      complex*16 kdiac(0:3,1)
      real*8 mompsiW(0:4,4)
      real*8 gluijMom(0:4,1:3,2:4),gluWijMom(0:4,1:3,2:4)
      complex*16 gluijMomc(0:4,1:3,2:4),gluWijMomc(0:4,1:3,2:4)
      complex*16 psi(2,-1:1,4), psiW(2,4),psi4(4,4)
      complex*16 amp(2,3:4,-1:1) ! colstr, 2nd, hel
      complex*16 ampv(2,3:4,-1:1),ampb(2,3:4,-1:1),ampvcross1(2,3:4,-1:1)
      complex*16 gluij(0:3,-1:1,1:3,2:4),gluWij(0:3,1:3,2:4),tempvirt
      complex*16 tempCurr(0:3,-1:1,4)
      complex*16 dotcc
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

      
      integer i,j,h,q,mu,tempi4(4),cross,col(2),crosssign,colstr,helcfg,comp,ngaugeh(-4:4)
      integer minhel, maxhel, in,out,in2,out2,hel(4),qtype(3),rndhel,gtesthel,gtestq,temphel ! fermion line 1:  qtype(1) > W > qtype(2), if cross==1: only qtype(3) needs to be calculated for 2nd qline

      real*8 deltaCVirt
      parameter (deltaCVIRT = -0d0)


      logical isNaN
      external isnan
      logical recalculatePSpoint
      common /comRecalcPSpoint/ recalculatePSpoint


      logical lpenb,lbox, lborn,lborn2,lcsi, wardidtest(0:8),lhel(-1:1),lpasshtest(3)
      complex*16 tempborn,gluren(3)
      complex*16 penbamp(3),penbres,penbborn
      complex*16 boxgauge(2,3),boxres(5,3,3:4,-1:1),boxgaugeb(2,3),boxborn(3,3:4,-1:1),nothing(8),checkg(5,3,2),checkgb(2,2)


      real*8 td,gaugeprecision
      real*8 CF,CA,CFCA,TR
      parameter( CF = 4d0/3d0, CA = 3d0, CFCA=-1d0/6d0, TR=0.5d0)
      logical ltop
      integer div
      logical ldebug
      
      
      ltop = NFL.eq.5

      ldebug = .false.
      
      do i=0,2; do j=3,4
        resB(i,j)=0d0
        resV(i,j)=0d0
        resI(i,j)=0d0
      enddo; enddo
      
      lborn= nlo.eq.1
      lCSI = nlo.eq.2
      lbox = nlo.eq.2
      lpenb= nlo.eq.2
      lborn2=nlo.eq.2
      if (ldoblha) then
        mursq(1,1) = blha_scale**2
        lborn=  .true. 
        lCSI =  .false. 
        lbox =  .true. 
        lpenb=  .true. 
        lborn2= .true. 
      endif
      
!       lborn = .true.
!       lCSI  = .true.
!       lbox  = .true.
!       lpenb = .true.
!       lborn2= .true.
     
      qhel=0
      div = 0
      
c select a random helicity for gauge checks
      rn = 2d0*rn
      rndhel = (-1)**INT(rn)

      musq = mursq(1,1)
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
      
      qtype(1) = nbos ! (W+)? u type : d type
      qtype(2) = 7-qtype(1)

c get bra and ket 
      tempi4(1) = fsign(1)
      tempi4(2) = fsign(3)
      tempi4(3) = fsign(2)
      tempi4(4) = fsign(4)
      call psi0m(4,pbar,tempi4,psi)
      
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

c new convention
c in   >W>>>>W>  out
c          S         , WZ
c	   S
c in2  >>>>>>>>  out2

c get <qW| and |Wq>
      do in = 1,3,2
	  call ket2c(psi(1,-1,in),.true.,pdia(0,in),-1,momW(0,1),epsW(0,1),psiW(1,in),mompsiW(0,in))
	  out = in+1
	  call bra2c(psi(1,-1,out),.true.,pdia(0,out),-1,momW(0,1),epsW(0,1),psiW(1,out),mompsiW(0,out))
      enddo
      
      minhel = max(qhel,1)
      maxhel = qhel
      if (maxhel.eq.0) maxhel = 3

      do in = 1,3,2  ! determine q line coupling to W
      do out= 2,4,2
      
      in2=4-in
      out2=6-out
      colstr= (in+out)/2
      
c calculate glu currents
      if(gluWij(0,in,out).eq.0d0) then
	call curr(-1,psi(1,-1,out),psiW(1,in),tempCurr(0,-1,1))
	call curr(-1,psiW(1,out),psi(1,-1,in),tempCurr(0,-1,2))
	do mu=0,3
	  gluWij(mu,in,out) = -tempCurr(mu,-1,1)-tempCurr(mu,-1,2)
	enddo
	
c calculate lower line currents
	call curr(1,psi(1,-1,out2),psi(1,-1,in2),gluij(0,-1,in2,out2))
      endif

c eventually calc amplitudes
      do helcfg = minhel,maxhel
c test for reasonable combination of helcfg, in, out
      lpasshtest(helcfg) = .not.(colstr.eq.3  .or.  helcfg+colstr.eq.4)
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
     &		* dotcc(gluWij(0,in,out),gluij(0,hel(in2),in2,out2))
	if(colstr.eq.1.or.in .eq.3) amp(colstr,qtype(1),hel(in2)) = tempvirt
	if(colstr.eq.1.or.out.eq.4) amp(colstr,qtype(2),hel(in2)) = tempvirt
	
	if(in.eq.3) amp(colstr,3,hel(in2)) = -amp(colstr,3,hel(in2)) ! flipped fermions
	if(out.eq.4)amp(colstr,4,hel(in2)) = -amp(colstr,4,hel(in2))

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
     -      3*(-((logmuosij(1,2) - logmuosij(1,4) - logmuosij(2,3) 
     &          + logmuosij(3,4))*(3*squared(1,1) - 2*squared(1,2) - 21*squared(2,2)))/9. + 
     -         (8*(logmuosij(1,4) + logmuosij(2,3))*(3*squared(1,1) 
     &          - 2*squared(1,2) + 3*squared(2,2)))/9. + 
     -         (-logmuosij(1,3) + logmuosij(1,4) + logmuosij(2,3) - logmuosij(2,4))*
     -          ((20*squared(1,2))/9. - (2*(squared(1,1) + squared(2,2)))/3.))))/(2.*Pi)
	  resI(1,i) = resI(1,i)-(als(1,1)*(-480 + 64*Pi**2 + 9*logmuosij(1,2)
     &           + 3*logmuosij(1,2)**2 - 18*logmuosij(1,3) - 6*logmuosij(1,3)**2 - 63*logmuosij(1,4) - 
     -       21*logmuosij(1,4)**2 - 63*logmuosij(2,3) - 21*logmuosij(2,3)**2 
     &          - 18*logmuosij(2,4) - 6*logmuosij(2,4)**2 + 9*logmuosij(3,4) + 
     -       3*logmuosij(3,4)**2)*squared(1,1))/(18.*Pi)
	  resI(2,i) = resI(2,i)-(als(1,1)*(-480 + 64*Pi**2 - 63*logmuosij(1,2) 
     &          - 21*logmuosij(1,2)**2 - 18*logmuosij(1,3) - 6*logmuosij(1,3)**2 + 9*logmuosij(1,4) + 
     -       3*logmuosij(1,4)**2 + 9*logmuosij(2,3) + 3*logmuosij(2,3)**2 
     &          - 18*logmuosij(2,4) - 6*logmuosij(2,4)**2 - 63*logmuosij(3,4) - 
     -       21*logmuosij(3,4)**2)*squared(2,2))/(18.*Pi) 
	  else 
	    resI(0,i) = resI(0,i) +(als(1,1)*((16*(21 + 9*alphaDip - 4*Pi**2
     &           - 9*Log(alphaDip) - 6*Log(alphaDip)**2)*(3*squared(1,1)
     &           - 2*squared(1,2) + 3*squared(2,2)))/
     -       27. + (-((logmuosij(1,2)**2 - logmuosij(1,4)**2 - logmuosij(2,3)**2 + logmuosij(3,4)**2)*
     -            (3*squared(1,1) - 2*squared(1,2) - 21*squared(2,2))) + 
     -         2*(logmuosij(1,3)**2 - logmuosij(1,4)**2 - logmuosij(2,3)**2 + logmuosij(2,4)**2)*
     -          (3*squared(1,1) - 10*squared(1,2) + 3*squared(2,2)) + 
     -         8*(logmuosij(1,4)**2 + logmuosij(2,3)**2)*(3*squared(1,1)
     &           - 2*squared(1,2) + 3*squared(2,2)))/9. + 
     -      3*(-((logmuosij(1,2) - logmuosij(1,4) - logmuosij(2,3) 
     &          + logmuosij(3,4))*(3*squared(1,1) - 2*squared(1,2) - 21*squared(2,2)))/9. + 
     -         (8*(logmuosij(1,4) + logmuosij(2,3))*(3*squared(1,1) 
     &          - 2*squared(1,2) + 3*squared(2,2)))/9. + 
     -         (-logmuosij(1,3) + logmuosij(1,4) + logmuosij(2,3) 
     &          - logmuosij(2,4))*
     -          ((20*squared(1,2))/9. - (2*(squared(1,1) + squared(2,2)))/3.))))/(2.*Pi)
	    resI(1,i) = resI(1,i) + (als(1,1)*(336 + 144*alphaDip - 64*Pi**2
     &           - 144*Log(alphaDip) - 96*Log(alphaDip)**2 - 9*logmuosij(1,2) - 3*logmuosij(1,2)**2 + 
     -         18*logmuosij(1,3) + 6*logmuosij(1,3)**2 + 63*logmuosij(1,4) 
     &          + 21*logmuosij(1,4)**2 + 63*logmuosij(2,3) + 21*logmuosij(2,3)**2 + 
     -         18*logmuosij(2,4) + 6*logmuosij(2,4)**2 - 9*logmuosij(3,4) 
     &          - 3*logmuosij(3,4)**2)*squared(1,1))/(18.*Pi)
	    resI(2,i) = resI(2,i) + (als(1,1)*(336 + 144*alphaDip - 64*Pi**2 
     &          - 144*Log(alphaDip) - 96*Log(alphaDip)**2 + 63*logmuosij(1,2) + 21*logmuosij(1,2)**2 + 
     -         18*logmuosij(1,3) + 6*logmuosij(1,3)**2 - 9*logmuosij(1,4) 
     &          - 3*logmuosij(1,4)**2 - 9*logmuosij(2,3) - 3*logmuosij(2,3)**2 + 
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
     -         (8*(logmuosij(1,4) + logmuosij(2,3))*(3*squared(1,1) 
     &          - 2*squared(1,2) + 3*squared(2,2)))/9. + 
     -         (-logmuosij(1,3) + logmuosij(1,4) + logmuosij(2,3) - logmuosij(2,4))*
     -          ((20*squared(1,2))/9. - (2*(squared(1,1) + squared(2,2)))/3.))))/(2.*Pi)
	  resI(1,i) = resI(1,i) -(als(1,1)*(-24 + logmuosij(1,2) 
     &          - 2*logmuosij(1,3) - 7*logmuosij(1,4) - 7*logmuosij(2,3) - 2*logmuosij(2,4) + logmuosij(3,4))*
     -     squared(1,1))/(3.*Pi)
	  resI(2,i) = resI(2,i) +(als(1,1)*(24 + 7*logmuosij(1,2)
     &           + 2*logmuosij(1,3) - logmuosij(1,4) - logmuosij(2,3) + 2*logmuosij(2,4) + 7*logmuosij(3,4))*
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
	resv(1,i) = resv(1,i) + als(1,1)/pi*(-(11d0/6d0*CA-2d0/3d0*TR*NFL)*log(musq/mursq(1,1))+log(musq/tmass/tmass)/3d0)
     $                          *2d0*squared(1,1)
	resv(2,i) = resv(2,i) + als(1,1)/pi*(-(11d0/6d0*CA-2d0/3d0*TR*NFL)*log(musq/mursq(1,1))+log(musq/tmass/tmass)/3d0)
     $                          *2d0*squared(2,2)
	resv(0,i) = resv(0,i) + als(1,1)/pi*(-(11d0/6d0*CA-2d0/3d0*TR*NFL)*log(musq/mursq(1,1))+log(musq/tmass/tmass)/3d0)
     &          *2d0*(squared(1,1)+squared(2,2)-2d0/3d0*squared(1,2))
      endif

	enddo
      enddo


c---------------------------------
c	VIRTUALS start here
c---------------------------------

c  change momentum convention-> all momenta flowing in
      do mu=0,3
	pdia(mu,2)=-pdia(mu,2)
	pdia(mu,4)=-pdia(mu,4)
	kdia(mu,1)=-momW(mu,1)
	kdiac(mu,1) = kdia(mu,1)
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
	  gluijMomc(mu,in,out) = gluijMom(mu,in,out)
	  gluWijMomc(mu,in,out) = gluWijMom(mu,in,out)
	enddo
	gluijMom(4,in,out) = mass2(gluijMom(0,in,out))
	gluWijMom(4,in,out) = mass2(gluWijMom(0,in,out))
	gluijMomc(4,in,out) = gluijMom(4,in,out)
	gluWijMomc(4,in,out) = gluWijMom(4,in,out)
      enddo; enddo
      
      if(ldebug) then
	print*
	print*,"diagtophys", id(1), id(3), id(2), id(4)
	print*,"fsign     ",fsign(1),fsign(3),fsign(2),fsign(4)
	print*,'pdia:'
	do i=1,4
	  print*,(pdia(mu,i),mu=0,3)
	enddo
	do i=1,1
	  print*,(kdia(mu,i),mu=0,3)
	enddo
	print*,"mom conservation"
! 	print*,( sum( (/ kdia(mu,1:1),pdia(mu,1:4) /) ), mu=0,3)
	
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

      do h =-1,1,2
	if(lhel(h)) then
	  
c self energies for in>W>out in2>out2
	  tempborn = dotcc(gluWij(0,in,out),gluij(0,h,in2,out2))/gluijMom(4,in2,out2)*G**2*crosssign
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
	endif
      enddo
      
      endif ! lborn2
      
c ---------------------------------------
c $$$$$$$$$$ B O X E S $$$$$$$$$$$$$$$$$$
c ---------------------------------------
      
      if(lbox) then

      if(ldebug) then
	print*
	print*,'box'
      endif

      do i=-1,1,2
        ngaugeh(i)=9
      enddo
      ngaugeh(gtesthel) = 3


      do i=1,3; do j=3,4; do h=-1,1,2
        boxborn(i,j,h) = 0d0
      enddo; enddo; enddo

c upper line boxes

c in> G> W> out, in2>out2 
      comp=1
      if(lhel(-1)) then
	call boxlineABETotal(pdia(0,in),gluijMom(0,in2,out2),kdia(0,1),pdia(0,out),
     &            psi(1,-1,out),psi(1,-1,in),gluij(0,-1,in2,out2),epsW(0,1),
     &            -1,musq,1,2,ngaugeh(-1),comp,boxgauge(1,1),boxres(1,1,3,-1),boxgaugeb(1,1),boxborn(1,3,-1),div)
	comp =-1
      endif
      
      if(lhel(1)) then
	call boxlineABETotal(pdia(0,in),gluijMom(0,in2,out2),kdia(0,1),pdia(0,out),
     &            psi(1,-1,out),psi(1,-1,in),gluij(0, 1,in2,out2),epsW(0,1),
     &            -1,musq,1,2,ngaugeh(1),comp,boxgauge(1,1),boxres(1,1,3,1),boxgaugeb(1,1),boxborn(1,3,1),div)
	comp =-1
      endif

c	for gauge test of nonAbe boxes
      call boxlineABETotal(pdia(0,in),gluijMom(0,in2,out2),kdia(0,1),pdia(0,out),
     &            psi(1,-1,out),psi(1,-1,in),gluijMomc(0,in2,out2),epsW(0,1),
     &            -1,musq,1,2,9,-1,nothing(1),checkg(1,1,1),nothing(1),checkgb(1,1),div)
      

c in> W> G> out, in2>out2 
      comp=1
      if(lhel(-1)) then
	call boxlineABETotal(pdia(0,in),kdia(0,1),gluijMom(0,in2,out2),pdia(0,out),
     &            psi(1,-1,out),psi(1,-1,in),epsW(0,1),gluij(0,-1,in2,out2),
     &            -1,musq,1,3,ngaugeh(-1),comp,boxgauge(1,2),boxres(1,2,3,-1),boxgaugeb(1,2),boxborn(2,3,-1),div)
	comp =-1
      endif
      
      if(lhel(1)) then
	call boxlineABETotal(pdia(0,in),kdia(0,1),gluijMom(0,in2,out2),pdia(0,out),
     &            psi(1,-1,out),psi(1,-1,in),epsW(0,1),gluij(0, 1,in2,out2),
     &            -1,musq,1,3,ngaugeh(1),comp,boxgauge(1,2),boxres(1,2,3,1),boxgaugeb(1,2),boxborn(2,3,1),div)
	comp =-1
      endif

c	for gauge test of nonAbe boxes
      call boxlineABETotal(pdia(0,in),kdia(0,1),gluijMom(0,in2,out2),pdia(0,out),
     &            psi(1,-1,out),psi(1,-1,in),epsW(0,1),gluijMomc(0,in2,out2),
     &            -1,musq,1,3,9,-1,nothing(1),checkg(1,2,1),nothing(1),checkgb(2,1),div)

c non abelian
      comp=1
      if(lhel(-1)) then
	call BoxlineNoAbeTotal(pdia(0,in),kdia(0,1),pdia(0,out),gluijMom(0,in2,out2),
     &             psi(1,-1,out),psi(1,-1,in),epsW(0,1),gluij(0,-1,in2,out2),
     &             -1,musq,ngaugeh(-1),comp,boxgauge(1,3),boxres(1,3,3,-1),div)
	comp = -1
      endif

      if(lhel(1)) then
	call BoxlineNoAbeTotal(pdia(0,in),kdia(0,1),pdia(0,out),gluijMom(0,in2,out2),
     &             psi(1,-1,out),psi(1,-1,in),epsW(0,1),gluij(0,1,in2,out2),
     &             -1,musq,ngaugeh(1),comp,boxgauge(1,3),boxres(1,3,3,1),div)
	comp = -1
      endif

      call BoxlineNoAbeTotal(pdia(0,in),kdia(0,1),pdia(0,out),gluijMom(0,in2,out2),
     &             psi(1,-1,out),psi(1,-1,in),epsW(0,1),gluijMomc(0,in2,out2),
     &             -1,musq,9,-1,nothing,checkg(1,3,1),div)

c gauge tests 
      call box_gauge(pdia(0,in),gluijMom(0,in2,out2),kdia(0,1),pdia(0,out),musq,
     2    boxgauge(1,1),boxgaugeb(1,1),wardidtest(1))
      call box_gauge(pdia(0,in),kdia(0,1),gluijMom(0,in2,out2),pdia(0,out),musq,
     2    boxgauge(1,2),boxgaugeb(1,2),wardidtest(2))

c      check epsG -> pG
      if (abs(checkg(1,3,1)).lt.1d-3) then
	td=abs(-(checkg(2,2,1)+checkg(2,1,1))/2d0-checkg(1,3,1))
      else
	td=abs(-(checkg(2,2,1)+checkg(2,1,1))/2d0/checkg(1,3,1)+1d0)
      endif
      wardidtest(3) = td.lt.1d-2

      if(ldebug) print*,wardidtest(1),wardidtest(2),wardidtest(3),td,checkg(1,3,1)
!       pause

      if(ldebug) gaugeprecision = max(gaugeprecision,td)
      
c sum
      wardidtest(0) = wardidtest(1).and.wardidtest(2).and.wardidtest(3).or.ldebug
	boxcount2m = boxcount2m+1
	boxcount2sub(sub_number) = boxcount2sub(sub_number)+1
      if(.not.wardidtest(0)) then
	boxcountm = boxcountm+1
	boxcountsub(sub_number) = boxcountsub(sub_number)+1
        recalculatePSpoint= .true.
        
!         do i=0,2; do j=3,4; do k=0,2
!           res (i,j,k) = 0d0
!         enddo; enddo; enddo
! 	return
      endif

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

!       pause

      endif ! lbox

c ---------------------------------------
c $$$$$$$$$$ PENBOX  $$$$$$$$$$$$$$$$$$$$
c ---------------------------------------
      
      if(lpenb) then
      
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
      
c in > W > out
c start with gauge test
      temphel = 67
      if(gtesthel.eq.-1) temphel = 77
      call PenBoxT(pdia(0,out2),pdia(0,out),kdia(0,1),pdia(0,in),pdia(0,in2),
     &               psi4(1,out2),psi4(1,in2),psi4(1,out),psi4(1,in),kdiac(0,1),
     &               musq,1,0,0,penbamp,penbres,penbborn,div,temphel)
     
      if (abs(penbamp(1)).gt.1d-2) then
	td = abs((penbamp(2)+penbamp(3))/penbamp(1)+1d0)
      else
	td = abs(penbamp(1)+penbamp(2)+penbamp(3))
      endif
      
      wardidtest(0) = td.lt.1d-2 .or. ldebug
      if(ldebug) print*, wardidtest(0),td
      if(ldebug) gaugeprecision = max(gaugeprecision,td)
      
c passed test > calc amp
      do h=-1,1,2
      if(lhel(h)) then
	temphel = 67
	if(h.eq.-1) temphel = 77
	call PenBoxT(pdia(0,out2),pdia(0,out),kdia(0,1),pdia(0,in),pdia(0,in2),
     &               psi4(1,out2),psi4(1,in2),psi4(1,out),psi4(1,in),epsW(0,1),
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

      
	box2count2m = box2count2m + 1
	box2count2sub(sub_number) = box2count2sub(sub_number) + 1
      if(.not.wardidtest(0)) then
	box2countm = box2countm + 1
	box2countsub(sub_number) = box2countsub(sub_number) + 1
        recalculatePSpoint=.true.
      endif


c penbox  Cross
      if(ldebug) then
	print*
	print*, 'penbox Cross, in,out',in,out
      endif
      
c in > W > Z > out
c start with gauge test
      temphel = 67
      if(gtesthel.eq.-1) temphel = 77
      call PenBoxCrossT(pdia(0,in2),pdia(0,out),kdia(0,1),pdia(0,in),pdia(0,out2),
     &               psi4(1,in2),psi4(1,out2),psi4(1,out),psi4(1,in),kdiac(0,1),
     &               musq,1,0,0,penbamp,penbres,penbborn,div,temphel)

      if (abs(penbamp(1)).gt.1d-2) then
	td = abs((penbamp(2)+penbamp(3))/penbamp(1)+1d0)
      else
	td = abs(penbamp(1)+penbamp(2)+penbamp(3))
      endif
      
      wardidtest(0) = td.lt.1d-2 .or. ldebug
      if(ldebug) print*, wardidtest(0),td
      if(ldebug) gaugeprecision = max(gaugeprecision,td)
      
c passed test > calc amp
      do h=-1,1,2
      if(lhel(h)) then
	temphel = 67
	if(h.eq.-1) temphel = 77
	call PenBoxCrossT(pdia(0,in2),pdia(0,out),kdia(0,1),pdia(0,in),pdia(0,out2),
     &               psi4(1,in2),psi4(1,out2),psi4(1,out),psi4(1,in),epsW(0,1),
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
      if(.not.wardidtest(0)) then
	box2countm = box2countm + 1
	box2countsub(sub_number) = box2countsub(sub_number) + 1
        recalculatePSpoint=.true.
      endif

      endif ! lpenbox


      if(ldebug) then
	print*,'same born amp?'
        do i=1,2; do j=3,4
          print*,ampb(i,j,-1)/amp(i,j,-1)
          print*,ampb(i,j,1)/amp(i,j,1)
        enddo ; enddo
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
	
c same flavours
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
      
!       if(skipidFL(1)) res(:,qtype(1),1) = 0d0
!       if(skipidFL(2)) res(:,qtype(2),1) = 0d0
      
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
	print*,'sum I+virt+CT  ',((resV(i,j)+resI(i,j),i=1,2),j=3,4)
	print*,'cancellation   ',((abs((resI(i,j)+resV(i,j))/resI(i,j)),i=1,2),j=3,4)
! 	pause
      endif

      if(ldebug) then
	print*,"res Born: ",resB
	print*,"res Virt: ",((resV(i,j)+resI(i,j),i=1,2),j=3,4)
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

      if(recalculatePSpoint) then
        recalculatePSpoint=.false.
        do i=0,2; do q=3,4; do j=0,2
          res(i,q,j) = 0d0
        enddo; enddo; enddo
        resB = 0
        resV = 0
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
      endif 

      if(IsNan(res,0d0,18,0)) then
        if(ldebug.or..not.recalculatePSpoint) print*,"NaN in "__FILE__,__LINE__
        do i=0,2; do q=3,4; do j=0,2
          res(i,q,j) = 0d0
        enddo; enddo; enddo
      endif

      end       


      SUBROUTINE ug2dgwVirt(P,id,fsign,nbos,nlo,ghel2,ghel4,res) ! sum over glu hel if ghel ==0

      use globalvars, only: ldoblha

      implicit none
#include "common.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"
      real*8 P(0:3,max_p,max_kin), res(0:2) ! for gauge test
      integer id(4),fsign(4),nbos,nlo,lephel,ghel2 ,ghel4 
      logical lnlo
      
      complex*16 epsW(0:3,7),epsWV(0:3,-1:1,7) 
      complex*16 epsV(0:3,-1:1,3:4,-1:1,7) ! ( , mu-hel, u/d , quark-hel), coupl. const to quarks included
      complex*16 epsA(0:3,-1:1),epsZ(0:3,-1:1)
      real*8 momV(0:4,7),momW(0:4,7),momWV(0:4,7)
      common /qcdwznew/ epsW, epsWV, epsV,epsA,epsZ,momV,momW,momWV


      real*8 pbar(0:3,2),qbar(0:4,2),pdia(0:3,2),qdia(0:4,3),kdia(0:3,1),momTemp(0:4)
      complex*16 qdiac(0:4,3),kdiac(0:3,1)
      real*8 momPsiW(0:4,2),momPsiG(0:4,2,2)
      complex*16 Psi(2,-1:1,2), PsiW(2,2),PsiG(2,2,2)
      complex*16 epsG(0:3,3), amp(3)
      complex*16 gluW(0:3), tempcurr(0:3,-1:1,4)
      complex*16 ampb(3),ampv(3),ampi(3), tempborn, tempvirt,tempvirt2,tempvirt3,tempvirt4
      complex*16 ortb(3),ortv(3),orti(3),ort(3)
      complex*16 s1c, dotcc, dotrc,tempc
      real*8 mass2,musq,pi2
      parameter (pi2 =9.86960440108935861883449099987615114d0)
      integer i,j,mu,tempi2(2),qtype(2)
      integer glu1,glu2,div

      logical recalculatePSpoint
      common /comRecalcPSpoint/ recalculatePSpoint


      logical lsimple,lborn,lbox,lpent,lCSI,ldebug

      complex*16 boxgauge(2,3),boxres(5,3),boxgaugeb(2,3),boxborn(3),nothing(9),checkg(5,3,2),checkgb(2,2)
      complex*16 pengauge(3,8),penres(8,8),pengaugeb(3,8),penborn(8),floopres(6)
      real*8 gaugeReal(2,4),floopgauge(6),logmuosij(3,4)
      logical wardidtest(0:22),isNan
      external isNaN
      real*8 td,gaugeprecision,mjj2,dotrr
      external mjj2,dotrr
      real*8 CF,CA,CFCA,TR
      parameter( CF = 4d0/3d0, CA = 3d0, CFCA=-1d0/6d0, TR=0.5d0)
      logical ltop
      double precision  clr, xm2, xmg, b
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6)
      
      ltop = NFL.eq.5
      ldebug = .false.

	lborn  = nlo.eq.1
	
	lsimple= nlo.eq.2 ! fermion loops, glu triangle, self energies, counterterms
	lbox   = nlo.eq.2
	lpent  = nlo.eq.2
	lCSI   = nlo.eq.2 
      if (ldoblha) then
        mursq(1,1) = blha_scale**2
        lborn=  .true. 
        lsimple=.true. 
        lbox =  .true. 
        lpent=  .true. 
        lCSI =  .false. 
      endif
	
	if(nlo.eq.3) then
	  res(0) = 0d0
          res(1) = 0d0
          res(2) = 0d0
	  return
	endif


!       lborn  = .true.
!       lsimple= .true. ! self energies triangles, floops
!       lbox   = .true.
!       lpent  = .true.
!       lCSI   = .true.

      div = 0

      musq = mursq(1,1)
      if(musq.le.0d0) then
	print*,"unreasonable musq"
      endif
      
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
      
      qtype(1) = nbos ! (W+Z)? u type : d type
      qtype(2) = 7-qtype(1)
      

c get bra and ket 
      tempi2(1) = fsign(1)
      tempi2(2) = fsign(3)
      call psi0m(2,pbar,tempi2,Psi)
      
c get <q,glu| and |glu,q>

      call helvec(qbar(0,1),fsign(2),ghel2,epsG(0,1))
      call helvec(qbar(0,2),fsign(4),ghel4,epsG(0,2))
      
      do mu = 0,3
	  epsG(mu,3) = (2d0*dotrc(qdia(0,2),epsG(0,1))*epsG(mu,2)
     &		      -2d0*dotrc(qdia(0,1),epsG(0,2))*epsG(mu,1)
     &		      +dotcc(epsG(0,1),epsG(0,2))*(qdia(mu,1)-qdia(mu,2)))/qdia(4,3) ! * [ta,tb]

      enddo


c ---------------------------------------
c $$$$$$$$$$ B O R N $$$$$$$$$$$$$$$$$$$$
c ---------------------------------------

c get <q,V| and |V,q>
      call ket2c(psi(1,-1,1),.true.,pdia(0,1),-1,momW(0,1), epsW(0,1),PsiW(1,1),momPsiW(0,1))
      call bra2c(psi(1,-1,2),.true.,pdia(0,2),-1,momW(0,1), epsW(0,1),PsiW(1,2),momPsiW(0,2))
      do i=1,2
	call ket2c(psi(1,-1,1),.true.,pdia(0,1),-1,qdia(0,i), epsG(0,i),PsiG(1,1,i),momPsiG(0,1,i))
	call bra2c(psi(1,-1,2),.true.,pdia(0,2),-1,qdia(0,i), epsG(0,i),PsiG(1,2,i),momPsiG(0,2,i))
      enddo
      
c build currents
      call curr(-1,Psi(1,-1,2),PsiW(1,1),tempcurr(0,-1,1))
      call curr(-1,PsiW(1,2),Psi(1,-1,1),tempcurr(0,-1,2))
      do mu=0,3
        gluW(mu) = -tempcurr(mu,-1,1)-tempcurr(mu,-1,2)
      enddo

      amp(1) = -dotcc(epsG(0,3),gluW) ! col1 = in>g1>g2>out = tbta 
      amp(2) = -amp(1)


      do glu1=1,2
	glu2 = 3-glu1

	amp(glu1) = amp(glu1) + s1c(PsiG(1,2,glu2),epsW(0,1),.true.,-1,PsiG(1,1,glu1))
	amp(glu1) = amp(glu1) + s1c(PsiG(1,2,glu2),epsG(0,glu1),.false.,-1,PsiW(1,1))
	amp(glu1) = amp(glu1) + s1c(PsiW(1,2),epsG(0,glu2),.false.,-1,PsiG(1,1,glu1))
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
	qdiac(mu,1) = qdia(mu,1)
	qdiac(mu,2) = qdia(mu,2)
	qdiac(mu,3) = qdia(mu,3)
	kdiac(mu,1) = kdia(mu,1)
      enddo
      
      if(ldebug) then
	print*,'mom cons.', (pdia(mu,1)+pdia(mu,2)+qdia(mu,1)+qdia(mu,2)+kdia(mu,1),mu=0,3)
	gaugeprecision = 0d0
      endif
      
c ---------------------------------------
c $$$$$$$$$$ SIMPLE LOOPS $$$$$$$$$$$$$$$
c ---------------------------------------

      if(lsimple) then

c glu self enregie
      tempborn=dotcc(epsG(0,3),gluW)*GG(1)**2 ! 1/q^2 included in epsG
      call gluonSelf(qdia(4,3),musq,NFL,tempc,div)
      tempvirt = tempc*tempborn

      if(ltop) then
	call gluonSelfMass(qdia(4,3),musq,tmass,tempc,div)
	tempvirt = tempvirt + tempc*tempborn
      endif

      ampv(1) = ampv(1)-tempvirt
      ampv(2) = ampv(2)+tempvirt
      
      
c top loop for external gluons
c use onshell renormalisation
c      if(ltop) thenc
c	tempvirt =0d0
c	call gluonSelfMass(qdia(4,1),musq,tmass,tempc,div)
c	tempvirt = tempvirt + tempc
c	call gluonSelfMass(qdia(4,2),musq,tmass,tempc,div)
c	tempvirt = tempvirt + tempc
c	tempvirt=tempvirt/2d0
c	ampv(1) = ampv(1) + amp(1)*tempvirt
c	ampv(2) = ampv(2) + amp(2)*tempvirt
c      endif
	
      
c V- and U-loop contributions to 3 gluon vertex
      do mu=0,3
        momtemp(mu) = -qdia(mu,3)
      enddo
!       call Vert3G(qdia(0,1),qdia(0,2),momtemp(0),epsG(0,1),epsG(0,2),gluWV,musq,tempborn,-1)
!       tempborn = tempborn*GG(1)**2/qdia(4,3)
!       print*,'test',tempborn
      
      call Vert3G(qdia(0,1),qdia(0,2),momtemp(0),epsG(0,1),epsG(0,2),gluW,musq,tempvirt,div)
      tempvirt=tempvirt*GG(1)**2/qdia(4,3)

      ampv(1) = ampv(1)-tempvirt
      ampv(2) = ampv(2)+tempvirt
      
c glu glu glu* ! 2 directions have opposite sign but different color structure
      call FermionVertCouplingT(1d0,1d0,1d0,1d0,1d0,1d0,
     $		0d0,
     $		qdia(0,1),qdia(0,2),momtemp(0),
     $		epsG(0,1),epsG(0,2),gluW,
     $		musq,1,1,floopres(1),floopgauge(1),div)
      tempvirt = -floopres(1)/2d0/qdia(4,3)*G**2
c see comments in FermionVertEvenTemp for minus sign
c color structure is Tr[ta,tb,tc]tc=tatb/2-delta_ab/8
c routine with glu 1 and 2 interchanged will have different sign and color structure tbta/2-delta_ab/8
      
	ampv(2) = ampv(2)+tempvirt*NFL ! 5 light quarks
	ampv(1) = ampv(1)-tempvirt*NFL

	box2count2m = box2count2m + 1
	box2count2sub(sub_number) = box2count2sub(sub_number) + 1
      if(floopgauge(1).gt.1d-2) then
	box2countm = box2countm + 1
	box2countsub(sub_number) = box2countsub(sub_number) + 1
        recalculatePSpoint=.true.
      endif

      if(ltop) then
	call FermionVertCouplingT(1d0,1d0,1d0,1d0,1d0,1d0,
     $		tmass,
     $		qdia(0,1),qdia(0,2),momtemp(0),
     $		epsG(0,1),epsG(0,2),gluW,
     $		musq,1,1,floopres(1),floopgauge(1),div)
	tempvirt = -floopres(1)/2d0/qdia(4,3)*G**2
	
	  ampv(2) = ampv(2)+tempvirt
	  ampv(1) = ampv(1)-tempvirt
	if(floopgauge(1).gt.1d-2) then
          recalculatePSpoint=.true.
	endif
      endif

      
      endif ! lsimple

c ---------------------------------------
c $$$$$$$$$$ B O X E S $$$$$$$$$$$$$$$$$$
c ---------------------------------------

      if(lbox) then

      if(ldebug) then
	print*
	print*,'box'
      endif
      boxborn(1) = 0d0
      boxborn(2) = 0d0
      boxborn(3) = 0d0

c q> GG > WV> q
      call boxlineABETotal(pdia(0,1),qdia(0,3),kdia(0,1),pdia(0,2),
     &            psi(1,-1,2),psi(1,-1,1),epsG(0,3),epsW(0,1),
     &            -1,musq,1,2,3,1,boxgauge(1,1),boxres(1,1),boxgaugeb(1,1),boxborn(1),div)

c nA gauge test
      call boxlineABETotal(pdia(0,1),qdia(0,3),kdia(0,1),pdia(0,2),
     &            psi(1,-1,2),psi(1,-1,1),qdiac(0,3),epsW(0,1),
     &            -1,musq,1,2,9,-1,nothing(1),checkg(1,1,1),nothing(1),checkgb(1,1),div)

c q> WV> GG> q
      call boxlineABETotal(pdia(0,1),kdia(0,1),qdia(0,3),pdia(0,2),
     &            psi(1,-1,2),psi(1,-1,1),epsW(0,1),epsG(0,3),
     &            -1,musq,1,3,3,1,boxgauge(1,2),boxres(1,2),boxgaugeb(1,2),boxborn(2),div)

c nA gauge test
      call boxlineABETotal(pdia(0,1),kdia(0,1),qdia(0,3),pdia(0,2),
     &            psi(1,-1,2),psi(1,-1,1),epsW(0,1),qdiac(0,3),
     &            -1,musq,1,3,9,-1,nothing(1),checkg(1,2,1),nothing(1),checkgb(2,1),div)

c non Abelian
      call BoxlineNoAbeTotal(pdia(0,1),kdia(0,1),pdia(0,2),qdia(0,3),
     &            psi(1,-1,2),psi(1,-1,1),epsW(0,1),epsG(0,3),
     &            -1,musq,2,1,boxgauge(1,3),boxres(1,3),div)
      
c gauge tests 
      call box_gauge(pdia(0,1),qdia(0,3),kdia(0,1),pdia(0,2),musq,
     2    boxgauge(1,1),boxgaugeb(1,1),wardidtest(1))
      call box_gauge(pdia(0,1),kdia(0,1),qdia(0,3),pdia(0,2),musq,
     2    boxgauge(1,2),boxgaugeb(1,2),wardidtest(2))

c      check epsG -> pG
      if (abs(boxgauge(1,3)).lt.1d-6) then
	td=abs(-(checkg(2,2,1)+checkg(2,1,1))/2d0+boxgauge(1,3))
      else
	td=abs(-(checkg(2,2,1)+checkg(2,1,1))/2d0/boxgauge(1,3)+1d0)
      endif
      wardidtest(3) = td.lt.1d-3


      if(ldebug) gaugeprecision = max(gaugeprecision,td)
      if(ldebug) print*,wardidtest(1),wardidtest(2),wardidtest(3),td,checkg(1,3,1)
c sum
      wardidtest(0) = wardidtest(1).and.wardidtest(2).and.wardidtest(3)
      wardidtest(0) = wardidtest(0).or.ldebug.or.div.ne.0

      
	boxcount2m = boxcount2m+1
	boxcount2sub(sub_number) = boxcount2sub(sub_number) + 1
      if(.not.wardidtest(0)) then
	boxcountm = boxcountm+1
	boxcountsub(sub_number) = boxcountsub(sub_number) + 1
        recalculatePSpoint=.true.
      endif

      
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

      
      endif ! lbox
      
c ---------------------------------------
c $$$$$$$$$$ PENTAGONS $$$$$$$$$$$$$$$$$$
c ---------------------------------------

      if(lpent) then

      if(ldebug) then
	print*
	print*,'pent'
      endif

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
      
c glu1,glu2,WV
      call penlineAbeTotal(pdia(0,1),qdia(0,glu1),qdia(0,glu2),kdia(0,1),pdia(0,2),
     &          psi(1,-1,2),psi(1,-1,1),epsG(0,glu1),epsG(0,glu2),epsW(0,1),
     &          -1,musq,2,4,4,1,pengauge(1,1),penres(1,1),pengaugeb(1,1),penborn(1),div)
c glu1,WV,glu2

      call penlineAbeTotal(pdia(0,1),qdia(0,glu1),kdia(0,1),qdia(0,glu2),pdia(0,2),
     &          psi(1,-1,2),psi(1,-1,1),epsG(0,glu1),epsW(0,1),epsG(0,glu2),
     &          -1,musq,2,3,4,1,pengauge(1,2),penres(1,2),pengaugeb(1,2),penborn(2),div)

c WV,glu1,glu2
      call penlineAbeTotal(pdia(0,1),kdia(0,1),qdia(0,glu1),qdia(0,glu2),pdia(0,2),
     &          psi(1,-1,2),psi(1,-1,1),epsW(0,1),epsG(0,glu1),epsG(0,glu2),
     &          -1,musq,2,2,4,1,pengauge(1,3),penres(1,3),pengaugeb(1,3),penborn(3),div)

c nonAbe glu2
c glu1,WV
      call penlineNoAbe_pT(pdia(0,1),qdia(0,glu1),kdia(0,1),pdia(0,2),qdia(0,glu2),
     &          psi(1,-1,2),psi(1,-1,1),epsG(0,glu1),epsW(0,1),epsG(0,glu2),
     &          -1,musq,2,2,4,1,pengauge(1,4),penres(1,4),div)

c WV,glu1
      call penlineNoAbe_pT(pdia(0,1),kdia(0,1),qdia(0,glu1),pdia(0,2),qdia(0,glu2),
     &          psi(1,-1,2),psi(1,-1,1),epsW(0,1),epsG(0,glu1),epsG(0,glu2),
     &          -1,musq,2,3,4,1,pengauge(1,5),penres(1,5),div)

c 2 triple vertices
      call NoAbe2gaPenLineTotal(pdia(0,1),kdia(0,1),pdia(0,2),qdia(0,glu2),qdia(0,glu1),
     &		psi(1,-1,2),psi(1,-1,1),epsW(0,1),epsG(0,glu2),epsG(0,glu1),
     &		-1,musq,2,0,1,nothing(1),penres(1,6),gaugereal(1,1),div,1)


! c quad vertex
      call NoAbe4vertBoxLineTotal(pdia(0,1),kdia(0,1),pdia(0,2),qdia(0,glu2),qdia(0,glu1),
     &		psi(1,-1,2),psi(1,-1,1),epsW(0,1),epsG(0,glu2),epsG(0,glu1),
     &		-1,musq,2,0,1,nothing(1),penres(1,7),gaugereal(1,2),div,1)

! 
c gauge tests
c glu1,glu2,WV
      call pent_gauge(pdia(0,1),qdia(0,glu1),qdia(0,glu2),kdia(0,1),pdia(0,2),
     &          psi(1,-1,2),psi(1,-1,1),epsG(0,glu1),epsG(0,glu2),epsW(0,1),
     &          -1,musq,pengauge(1,1),pengaugeb(1,1),wardidtest(1))
! 
c glu1,WV,glu2
      call pent_gauge(pdia(0,1),qdia(0,glu1),kdia(0,1),qdia(0,glu2),pdia(0,2),
     &          psi(1,-1,2),psi(1,-1,1),epsG(0,glu1),epsW(0,1),epsG(0,glu2),
     &          -1,musq,pengauge(1,2),pengaugeb(1,2),wardidtest(2))

c WV,glu1,glu2
      call pent_gauge(pdia(0,1),kdia(0,1),qdia(0,glu1),qdia(0,glu2),pdia(0,2),
     &          psi(1,-1,2),psi(1,-1,1),epsW(0,1),epsG(0,glu1),epsG(0,glu2),
     &          -1,musq,pengauge(1,3),pengaugeb(1,3),wardidtest(3))
! 
! c glu1,WV
      call pentNoAbe_gauge(1d-1,pdia(0,1),qdia(0,glu1),kdia(0,1),pdia(0,2),qdia(0,glu2),
     &          psi(1,-1,2),psi(1,-1,1),epsG(0,glu1),epsW(0,1),epsG(0,glu2),
     &          -1,musq,pengauge(1,4),td,wardidtest(4))
     
c WV,glu1
      call pentNoAbe_gauge(1d-1,pdia(0,1),kdia(0,1),qdia(0,glu1),pdia(0,2),qdia(0,glu2),
     &          psi(1,-1,2),psi(1,-1,1),epsW(0,1),epsG(0,glu1),epsG(0,glu2),
     &          -1,musq,pengauge(1,5),td,wardidtest(5))

      penres(1,7) = penres(1,7)/2d0 ! we calculate this contribution twice, but only need it once
      penres(2,7) = penres(2,7)/2d0 
      penres(3,7) = penres(3,7)/2d0 

      wardidtest(6) = abs(gaugereal(1,1)).lt.1d-1
      wardidtest(7) = abs(gaugereal(1,2)).lt.1d-1
      wardidtest(0)=wardidtest(1).and.wardidtest(2).and.wardidtest(3).and.wardidtest(4).and.wardidtest(5)
     &              .and.wardidtest(6).and.wardidtest(7)
      if(ldebug) print*,(wardidtest(i),i=0,7)
      wardidtest(0) = wardidtest(0).or.ldebug.or.div.ne.0

c sum
	pentcount2m = pentcount2m+1
	pentcount2sub(sub_number) = pentcount2sub(sub_number) + 1
      if(.not.wardidtest(0)) then
	pentcountm = pentcountm+1
	pentcountsub(sub_number) = pentcountsub(sub_number) + 1
        recalculatePSpoint=.true.
      endif

	tempborn = -penborn(1)-penborn(2)-penborn(3) ! penlineAbe and penlineNoAbe have different sign convention
	tempvirt = -penres(1,1)-penres(1,2)-penres(1,3)-(0,-1.5d0)*(penres(2,4)
     &          +penres(2,5))+( -1.5d0*penres(1,6)+ 1.5d0*(penres(2,7)
     &          +penres(3,7)-2d0*penres(1,7)))
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
      enddo ! glu

      endif ! lpent


c---------------------------------
c	I Operator
c---------------------------------
c use amplitude calculated with loop routines ampb
c to split up into W, and WV contributions
c amp still contains all contributions

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
	ampi(3) = ((ampb(1) + ampb(2))*logmuosij(1,2) - ampb(2)*logmuosij(1,3) - ampb(1)*logmuosij(1,4) - ampb(1)*logmuosij(2,3)- 
     $              ampb(2)*logmuosij(2,4) + ampb(1)*logmuosij(3,4) + ampb(2)*logmuosij(3,4))/2.
      else
        if(alphaDip.eq.1d0) then
	  ampi(1) = -(ampb(1)*(-840 + 32*NFL + 104*Pi**2 + 9*logmuosij(1,2) + 3*logmuosij(1,2)**2 + 3*(-30 + NFL)*logmuosij(1,3) - 
     $                27*logmuosij(1,3)**2 - 90*logmuosij(2,4) + 3*NFL*logmuosij(2,4) - 27*logmuosij(2,4)**2 - 99*logmuosij(3,4) +
     $                 6*NFL*logmuosij(3,4) - 27*logmuosij(3,4)**2))/18.
	  ampi(2) = -(ampb(2)*(-840 + 32*NFL + 104*Pi**2 + 9*logmuosij(1,2) + 3*logmuosij(1,2)**2 + 3*(-30 + NFL)*logmuosij(1,4) - 
     $                27*logmuosij(1,4)**2 - 90*logmuosij(2,3) + 3*NFL*logmuosij(2,3) - 27*logmuosij(2,3)**2 - 99*logmuosij(3,4) + 
     $                6*NFL*logmuosij(3,4) - 27*logmuosij(3,4)**2))/18.
	  ampi(3) = (27*(ampb(1) + ampb(2))*logmuosij(1,2) + 9*(ampb(1) + ampb(2))*logmuosij(1,2)**2 + 
     $              (-30 + NFL)*ampb(2)*logmuosij(1,3) - 9*ampb(2)*logmuosij(1,3)**2 - 30*ampb(1)*logmuosij(1,4) + 
     $              NFL*ampb(1)*logmuosij(1,4) - 9*ampb(1)*logmuosij(1,4)**2 - 30*ampb(1)*logmuosij(2,3) + 
     $              NFL*ampb(1)*logmuosij(2,3) - 9*ampb(1)*logmuosij(2,3)**2 - 30*ampb(2)*logmuosij(2,4) + 
     $              NFL*ampb(2)*logmuosij(2,4) - 9*ampb(2)*logmuosij(2,4)**2 + 33*ampb(1)*logmuosij(3,4) - 
     $              2*NFL*ampb(1)*logmuosij(3,4) + 33*ampb(2)*logmuosij(3,4) - 2*NFL*ampb(2)*logmuosij(3,4) + 
     $              9*ampb(1)*logmuosij(3,4)**2 + 9*ampb(2)*logmuosij(3,4)**2)/36.
	else
	  ampi(1) = -(ampb(1)*(-570 - 270*alphaDip + 20*NFL + 12*alphaDip*NFL + 104*Pi**2 - 6*(-45 + 2*NFL)*Log(alphaDip) + 1
     $              56*Log(alphaDip)**2 + 9*logmuosij(1,2) + 3*logmuosij(1,2)**2 - 90*logmuosij(1,3) + 3*NFL*logmuosij(1,3) - 
     $              27*logmuosij(1,3)**2 - 90*logmuosij(2,4) + 3*NFL*logmuosij(2,4) - 27*logmuosij(2,4)**2 - 99*logmuosij(3,4) + 
     $              6*NFL*logmuosij(3,4) - 27*logmuosij(3,4)**2))/18.
	  ampi(2) = -(ampb(2)*(-570 - 270*alphaDip + 20*NFL + 12*alphaDip*NFL + 104*Pi**2 - 6*(-45 + 2*NFL)*Log(alphaDip) + 
     $              156*Log(alphaDip)**2 + 9*logmuosij(1,2) + 3*logmuosij(1,2)**2 - 90*logmuosij(1,4) + 3*NFL*logmuosij(1,4) - 
     $              27*logmuosij(1,4)**2 - 90*logmuosij(2,3) + 3*NFL*logmuosij(2,3) - 27*logmuosij(2,3)**2 - 99*logmuosij(3,4) + 
     $              6*NFL*logmuosij(3,4) - 27*logmuosij(3,4)**2))/18.
	  ampi(3) = (27*(ampb(1) + ampb(2))*logmuosij(1,2) + 9*(ampb(1) + ampb(2))*logmuosij(1,2)**2 + (-30 + NFL)*
     $              ampb(2)*logmuosij(1,3) - 9*ampb(2)*logmuosij(1,3)**2 - 30*ampb(1)*logmuosij(1,4) + NFL*ampb(1)*logmuosij(1,4) - 
     $              9*ampb(1)*logmuosij(1,4)**2 - 30*ampb(1)*logmuosij(2,3) + NFL*ampb(1)*logmuosij(2,3) - 
     $              9*ampb(1)*logmuosij(2,3)**2 - 30*ampb(2)*logmuosij(2,4) + NFL*ampb(2)*logmuosij(2,4) - 
     $              9*ampb(2)*logmuosij(2,4)**2 + 33*ampb(1)*logmuosij(3,4) - 2*NFL*ampb(1)*logmuosij(3,4) + 
     $              33*ampb(2)*logmuosij(3,4) - 2*NFL*ampb(2)*logmuosij(3,4) + 9*ampb(1)*logmuosij(3,4)**2 + 
     $              9*ampb(2)*logmuosij(3,4)**2)/36.
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

c---------------------------------
c	Counter terms
c---------------------------------
      if(div.eq.1) then
	ampv(1) = ampv(1) -ampb(1)*2d0*(11d0/6d0*CA-2d0/3d0*TR*NFL)
	ampv(2) = ampv(2) -ampb(2)*2d0*(11d0/6d0*CA-2d0/3d0*TR*NFL)
      endif
      if(div.eq.0.and.ltop.and.lsimple) then 
	ampv(1) = ampv(1) -ampb(1)*2d0*(11d0/6d0*CA-2d0/3d0*TR*NFL)*log(musq/mursq(1,1))
	ampv(2) = ampv(2) -ampb(2)*2d0*(11d0/6d0*CA-2d0/3d0*TR*NFL)*log(musq/mursq(1,1))
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
! 	pause
      endif
      res(1) = res(1)*als(1,1)/2d0/pi ! 2Re *als/4pi
      res(0) = res(0)*als(1,1)/2d0/pi 

      if(recalculatePSpoint) then
        recalculatePSpoint=.false.
        res(0)=0d0; res(1)=0d0; res(2)=0d0
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
      
      
      if(IsNan(res,0d0,3,0)) then
        if(ldebug.or..not.recalculatePSpoint) print*,"NaN in "__FILE__,__LINE__
        do i=0,2
          res(i) = 0d0
        enddo
      endif
      
      end

      SUBROUTINE ug2dggwOwnAllHel(P,id,fsign,nbos, res) ! sum over glu hel if ghel ==0

      use globalvars, only: ldoblha

      implicit none 
#include "common.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"
      real*8 P(0:3,max_p,max_kin), res(0:6)
      integer id(5),fsign(5),nbos 

      complex*16 epsW(0:3,7),epsWV(0:3,-1:1,7) 
      complex*16 epsV(0:3,-1:1,3:4,-1:1,7) ! ( , mu-hel, u/d , quark-hel), coupl. const to quarks included
      complex*16 epsA(0:3,-1:1),epsZ(0:3,-1:1)
      real*8 momV(0:4,7),momW(0:4,7),momWV(0:4,7)
      common /qcdwznew/ epsW, epsWV, epsV,epsA,epsZ,momV,momW,momWV

      
      real*8 pbar(0:3,2),kbar(0:4,3),pdia(0:3,2),kdia(0:4,7),momTemp(0:4) ! kdia for glu: 1 2 3 12 13 23 123
      real*8 momPsiW(0:4,2),momPsiG(0:4,2,7)
      real*8 momPsiWG(0:4,2,3)
      real*8 triD(0:4,3,2,2,2),triV(0:4,3,2,2,2) ! 3-gluon polarications (see Nucl Phys B313 560, (5.13-14) indices: triD(mu,i) with j<k, triV(j)
      complex*16 Psi(2,-1:1,2), PsiW(2,2), PsiG(2,2,7,2) ! (,ferm,hel,glu,pol)
      complex*16 PsiWG(2,2,3,2)
      complex*16 PsiTemp(2,-1:1), amp(6)
      complex*16 s1r,temp
      real*8 epsG(0:3,3,2),eps3G(0:3,4:6,0:2,0:2,0:2),epsTemp(0:3)
      real*8 mass2, colsumsq,dotrr
      integer denom(6)
      double precision CF(6,6)
      
      integer i,j,mu, tempi2(2), lepkin,qtype(2),ggsign(3:5),pol1,pol2,pol3,pol(0:3)
      integer ig(3),perm, permutations(3,6)

c permutations convenient with Madgraph u g -> e+ ve d g mu+ mu- g  
      data (permutations(i,1),i=1,3) / 1, 2, 3/
      data (permutations(i,2),i=1,3) / 1, 3, 2/
      data (permutations(i,3),i=1,3) / 3, 1, 2/
      data (permutations(i,4),i=1,3) / 2, 1, 3/
      data (permutations(i,5),i=1,3) / 2, 3, 1/
      data (permutations(i,6),i=1,3) / 3, 2, 1/

c color factors from MG
      DATA Denom(1  )/            9/                                       
      DATA (CF(i,1  ),i=1  ,6  ) /    64,   -8,    1,   -8,    1,   10/    
C               T[ 5, 1, 2, 6, 9]                                          
      DATA Denom(2  )/            9/                                       
      DATA (CF(i,2  ),i=1  ,6  ) /    -8,   64,   -8,    1,   10,    1/    
C               T[ 5, 1, 2, 9, 6]                                          
      DATA Denom(3  )/            9/                                       
      DATA (CF(i,3  ),i=1  ,6  ) /     1,   -8,   64,   10,    1,   -8/    
C               T[ 5, 1, 9, 2, 6]                                          
      DATA Denom(4  )/            9/                                       
      DATA (CF(i,4  ),i=1  ,6  ) /    -8,    1,   10,   64,   -8,    1/    
C               T[ 5, 1, 6, 2, 9]                                          
      DATA Denom(5  )/            9/                                       
      DATA (CF(i,5  ),i=1  ,6  ) /     1,   10,    1,   -8,   64,   -8/    
C               T[ 5, 1, 6, 9, 2]                                          
      DATA Denom(6  )/            9/                                       
      DATA (CF(i,6  ),i=1  ,6  ) /    10,    1,   -8,    1,   -8,   64/    

      lepkin = 1
      do mu=0,3
	pbar(mu,1)=p(mu,id(1),1)
	pbar(mu,2)=p(mu,id(3),1)
	kbar(mu,1)=p(mu,id(2),1)
	kbar(mu,2)=p(mu,id(4),1)
	kbar(mu,3)=p(mu,id(5),1)
	
	pdia(mu,1)=pbar(mu,1)*fsign(1)
	pdia(mu,2)=pbar(mu,2)*fsign(3)
	kdia(mu,1)=kbar(mu,1)*fsign(2)
	kdia(mu,2)=kbar(mu,2)*fsign(4)
	kdia(mu,3)=kbar(mu,3)*fsign(5)
	kdia(mu,4)=kdia(mu,1)+kdia(mu,2)
	kdia(mu,5)=kdia(mu,1)+kdia(mu,3)
	kdia(mu,6)=kdia(mu,2)+kdia(mu,3)
	kdia(mu,7)=kdia(mu,1)+kdia(mu,6)
      enddo      
      kbar(4,1) = 0d0
      kbar(4,2) = 0d0
      kbar(4,3) = 0d0
      kdia(4,1) = 0d0
      kdia(4,2) = 0d0
      kdia(4,3) = 0d0
      kdia(4,4) = mass2(kdia(0,4))
      kdia(4,5) = mass2(kdia(0,5))
      kdia(4,6) = mass2(kdia(0,6))
      kdia(4,7) = mass2(kdia(0,7))
      
      qtype(1) = nbos ! (W+Z)? u type : d type
      qtype(2) = 7-qtype(1)
      
c get bra and ket 
      tempi2(1) = fsign(1)
      tempi2(2) = fsign(3)
      call psi0m(2,pbar,tempi2,Psi)
c get eps glu
      do pol1 = 1,2
	call polvec(kbar(0,1),pol1,epsG(0,1,pol1))
	call polvec(kbar(0,2),pol1,epsG(0,2,pol1))
	call polvec(kbar(0,3),pol1,epsG(0,3,pol1))
      enddo
      
      do pol1=1,2; do pol2=1,2
      do mu = 0,3
	eps3G(mu,4,pol1,pol2,0)= (2d0*dotrr(kdia(0,2),epsG(0,1,pol1))*epsG(mu,2,pol2)
     &		    -2d0*dotrr(kdia(0,1),epsG(0,2,pol2))*epsG(mu,1,pol1)
     &		    +dotrr(epsG(0,1,pol1),epsG(0,2,pol2))*(kdia(mu,1)-kdia(mu,2)))/kdia(4,4)
	eps3G(mu,5,pol1,0,pol2)= (2d0*dotrr(kdia(0,3),epsG(0,1,pol1))*epsG(mu,3,pol2)
     &		    -2d0*dotrr(kdia(0,1),epsG(0,3,pol2))*epsG(mu,1,pol1)
     &		    +dotrr(epsG(0,1,pol1),epsG(0,3,pol2))*(kdia(mu,1)-kdia(mu,3)))/kdia(4,5)
	eps3G(mu,6,0,pol1,pol2)= (2d0*dotrr(kdia(0,3),epsG(0,2,pol1))*epsG(mu,3,pol2)
     &		    -2d0*dotrr(kdia(0,2),epsG(0,3,pol2))*epsG(mu,2,pol1)
     &		    +dotrr(epsG(0,2,pol1),epsG(0,3,pol2))*(kdia(mu,2)-kdia(mu,3)))/kdia(4,6)
      enddo
      enddo;enddo

      do pol1=1,2; do pol2=1,2; do pol3=1,2
      do mu=0,3
	triD(mu,1,pol1,pol2,pol3) = (2d0*(dotrr(epsG(0,1,pol1),kdia(0,6))*eps3G(mu,6,0,pol2,pol3)
     &		          -dotrr(kdia(0,1),eps3G(0,6,0,pol2,pol3))*epsG(mu,1,pol1))
     &		       +   dotrr(epsG(0,1,pol1),eps3G(0,6,0,pol2,pol3))*(kdia(mu,1)-kdia(mu,6)))/kdia(4,7)
	triD(mu,2,pol1,pol2,pol3) = (2d0*(dotrr(epsG(0,2,pol2),kdia(0,5))*eps3G(mu,5,pol1,0,pol3)
     &		          -dotrr(kdia(0,2),eps3G(0,5,pol1,0,pol3))*epsG(mu,2,pol2))
     &		       +   dotrr(epsG(0,2,pol2),eps3G(0,5,pol1,0,pol3))*(kdia(mu,2)-kdia(mu,5)))/kdia(4,7)
	triD(mu,3,pol1,pol2,pol3) = (2d0*(dotrr(epsG(0,3,pol3),kdia(0,4))*eps3G(mu,4,pol1,pol2,0)
     &		          -dotrr(kdia(0,3),eps3G(0,4,pol1,pol2,0))*epsG(mu,3,pol3))
     &		       +   dotrr(epsG(0,3,pol3),eps3G(0,4,pol1,pol2,0))*(kdia(mu,3)-kdia(mu,4)))/kdia(4,7)
	triV(mu,1,pol1,pol2,pol3) = (2d0*epsG(mu,1,pol1)*dotrr(epsG(0,2,pol2),epsG(0,3,pol3))
     &		        -  epsG(mu,2,pol2)*dotrr(epsG(0,1,pol1),epsG(0,3,pol3))
     &		        -  epsG(mu,3,pol3)*dotrr(epsG(0,1,pol1),epsG(0,2,pol2)))/kdia(4,7)
	triV(mu,2,pol1,pol2,pol3) = (2d0*epsG(mu,2,pol2)*dotrr(epsG(0,1,pol1),epsG(0,3,pol3))
     &		        -  epsG(mu,1,pol1)*dotrr(epsG(0,2,pol2),epsG(0,3,pol3))
     &		        -  epsG(mu,3,pol3)*dotrr(epsG(0,1,pol1),epsG(0,2,pol2)))/kdia(4,7)
	triV(mu,3,pol1,pol2,pol3) = -triV(mu,1,pol1,pol2,pol3)-triV(mu,2,pol1,pol2,pol3)
      enddo
      enddo;enddo;enddo

c get <q,glu| and |glu,q>
      do i=1,3
      do pol1=1,2
	call ket2r(psi(1,-1,1),.true.,pdia(0,1),-1,kdia(0,i), epsG(0,i,pol1),PsiG(1,1,i,pol1),momPsiG(0,1,i))
	call bra2r(psi(1,-1,2),.true.,pdia(0,2),-1,kdia(0,i), epsG(0,i,pol1),PsiG(1,2,i,pol1),momPsiG(0,2,i))
      enddo; enddo
c get <q,V| and |V,q>
      call ket2c(psi(1,-1,1),.true.,pdia(0,1),-1,momW(0,lepkin), epsW(0,lepkin),PsiW(1,1),momPsiW(0,1))
      call bra2c(psi(1,-1,2),.true.,pdia(0,2),-1,momW(0,lepkin), epsW(0,lepkin),PsiW(1,2),momPsiW(0,2))
c get <q,Vg| and |Vg,q>, 
      do i=1,3
      do pol1=1,2
	call ket2r(psiW(1,1),.false.,momPsiW(0,1),-1,kdia(0,i), epsG(0,i,pol1),PsiTemp(1,-1),momPsiWG(0,1,i))
	call ket2c(psiG(1,1,i,pol1),.false.,momPsiG(0,1,i),-1,momW(0,lepkin), epsW(0,lepkin),PsiTemp(1,1),momTemp)
	PsiWG(1,1,i,pol1) = Psitemp(1,-1) + Psitemp(1,1)
	PsiWG(2,1,i,pol1) = Psitemp(2,-1) + Psitemp(2,1)
	call bra2c(psiG(1,2,i,pol1),.false.,momPsiG(0,2,i),-1,momW(0,lepkin), epsW(0,lepkin),PsiTemp(1,-1),momPsiWG(0,2,i))
	call bra2r(psiW(1,2),.false.,momPsiW(0,2),-1,kdia(0,i), epsG(0,i,pol1),PsiTemp(1,1),momTemp)
	PsiWG(1,2,i,pol1) = Psitemp(1,-1) + Psitemp(1,1)
	PsiWG(2,2,i,pol1) = Psitemp(2,-1) + Psitemp(2,1)
      enddo
      enddo

      do i=0,6
	res(i) = 0d0
      enddo
      colsumsq = 0d0
c get amps ans ampsq for all polarisations
      do pol1=1,2; do pol2=1,2; do pol3=1,2

      do perm = 1,6
	ig(1) = permutations(1,perm)
	ig(2) = permutations(2,perm)
	ig(3) = permutations(3,perm)
	! 12 13 and 23 ordered
	if(perm.le.3) then; ggsign(3) = 1; else; ggsign(3) = -1; endif
	if(perm.le.2.or.perm.eq.4) then; ggsign(4) = 1; else; ggsign(4) = -1; endif
	if(perm/2.eq.2.or.perm.eq.1) then; ggsign(5) = 1; else; ggsign(5) = -1; endif

c 1 gluon splitting amplitudes
	i = ig(2)+ig(3)+1
	pol(1) = pol1; pol(2) = pol2; pol(3) = pol3
	pol(0) = pol(ig(1)); pol(ig(1)) = 0 ! swap to get convention of epsG and psi*G
	temp = s1r(psiWG(1,2,ig(1),pol(0)),eps3G(0,i,pol(1),pol(2),pol(3)),.true.,-1,Psi(1,-1,1))
     &        +s1r(psiG(1,2,ig(1),pol(0)),eps3G(0,i,pol(1),pol(2),pol(3)),.true.,-1,PsiW(1,1))
	amp(perm) = -temp*ggsign(i-1)

	i = ig(1)+ig(2)+1
	pol(ig(1)) = pol(0); pol(0) = pol(ig(3)); pol(ig(3)) = 0 ! swap back and swap for next calls
	temp = s1r(psiW(1,2),eps3G(0,i,pol(1),pol(2),pol(3)),.true.,-1,PsiG(1,1,ig(3),pol(0)))
     &        +s1r(psi(1,-1,2),eps3G(0,i,pol(1),pol(2),pol(3)),.true.,-1,PsiWG(1,1,ig(3),pol(0)))
	amp(perm) = amp(perm)-temp*ggsign(i-1)
	pol(ig(3)) = pol(0) ! swap back
	
c g(g1,g2,g3) amplitudes	
	do mu=0,3
	  epsTemp(mu) = ggsign(ig(2)+ig(3))*triD(mu,ig(1),pol1,pol2,pol3)-ggsign(ig(1)+ig(2))*triD(mu,ig(3),pol1,pol2,pol3)+
     $                        triV(mu,ig(2),pol1,pol2,pol3)
	enddo
        amp(perm)= amp(perm)
     &           +s1r(psi(1,-1,2),epsTemp,.true.,-1,PsiW(1,1))
     &           +s1r(psiW(1,2),epsTemp,.true.,-1,Psi(1,-1,1))

c amplitudes without glu splitting
	amp(perm) = amp(perm)
     &             +s1r(psiWG(1,2,ig(1),pol(ig(1))),epsG(0,ig(2),pol(ig(2))),.false.,-1,PsiG(1,1,ig(3),pol(ig(3))))
     &             +s1r(psiG(1,2,ig(1),pol(ig(1))),epsG(0,ig(2),pol(ig(2))),.false.,-1,PsiWG(1,1,ig(3),pol(ig(3))))
	amp(perm) = amp(perm) * G**3
      enddo ! perm

c BLHA colour factors
      if (ldoblha) then
        denom(:) = 1d0
        do i=1,6
          CF(i,i            ) =  blha_CA*blha_CF**3
          CF(i,mod(i+ 1,6)+1) =  blha_CF*(blha_CA-2d0*blha_CF)/4d0
          CF(i,mod(i+ 3,6)+1) =  blha_CF*(blha_CA-2d0*blha_CF)/4d0
          CF(i,mod( 8-i,6)+1) = -blha_CA*blha_CF**2*(blha_CA-2d0*blha_CF)/2d0
          CF(i,mod(10-i,6)+1) = -blha_CA*blha_CF**2*(blha_CA-2d0*blha_CF)/2d0
          CF(i,7-i          ) =  blha_CA*blha_CF*(blha_CA-2d0*blha_CF)*(blha_CA-blha_CF)/2d0
        enddo
      endif

c square amp
      DO I = 1, 6
          TEMP = (0.D0,0.D0)
          DO J = 1, 6
              TEMP = TEMP + CF(J,I)*AMP(J)
          ENDDO
          res(0) =res(0)+dble(TEMP*DCONJG(AMP(I))/DENOM(I))
          res(i) =res(i)+ dble(AMP(I)*DCONJG(AMP(I))/DENOM(I))
      ENDDO
      enddo; enddo;enddo ! pol
      
      colsumsq = 0d0
      do i=1,6
	colsumsq = colsumsq+res(i)
      enddo
      DO i=1,6
         res(i)=res(i)/colsumsq*res(0)
      ENDDO
      
      if (ldoblha) then
        if (blha_amptype(blha_curproc).eq.0) then ! tree
          blha_amp(1)=blha_amp(1)+
     &      res(0)
        else 
          call BLHA_amptypeerror(blha_amptype(blha_curproc),__FILE__,__LINE__)
        endif
      endif

      end

      subroutine amp4q1gw(P,id,fsign,nbos,qhel,ghel,res)

      use globalvars, only: ldoblha

      implicit none
#include "common.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"
      real*8 P(0:3,max_p,max_kin),res(0:4,3:4,0:1)
      integer id(5),fsign(5),nbos,qhel,ghel

      complex*16 epsW(0:3,7),epsWV(0:3,-1:1,7) 
      complex*16 epsV(0:3,-1:1,3:4,-1:1,7) ! ( , mu-hel, u/d , quark-hel), coupl. const to quarks included
      complex*16 epsA(0:3,-1:1),epsZ(0:3,-1:1)
      real*8 momV(0:4,7),momW(0:4,7),momWV(0:4,7)
      common /qcdwznew/ epsW, epsWV, epsV,epsA,epsZ,momV,momW,momWV

      
      real*8 pbar(0:3,5),pdia(0:3,5),momTemp(0:4)
      real*8 mompsiW(0:4,4)
      real*8 mompsiWG(0:4,4),mompsiG(0:4,4)
      complex*16 psi(2,-1:1,4), psiW(2,4)
      complex*16 psiG(2,-1:1,4,-1:1), psiWG(2,4,-1:1)
      complex*16 psiTemp(2,-1:1), amp(4,2,3:4,3,-1:1) ! numbering as in MG u u -> e+ ve u d mu+ mu- g (1<>3, 2<>4 interchanged!), colstr, 2nd, hel
      complex*16 epsG(0:3,-1:1)
      complex*16 gluWij(0:3,1:3,2:4),temp
      complex*16 gluWGij(0:3,1:3,2:4,-1:1,2)
      complex*16 gluij(0:3,-1:1,1:3,2:4),gluGij(0:3,-1:1,1:3,2:4,-1:1,2)

      complex*16 tempCurr(0:3,-1:1,4),amp_ort
      complex*16 s1c, dotrc
      complex*16 dotcc
      external s1c, dotcc, dotrc,dotrr,qjj2,mass2
      real*8 mass2,qjj2
      real*8 tempR,colsq(0:4,2)
      real*8 dotrr
      logical skip(3:4)
      
      integer i,j,k,l,m,mu,hel(4), tempi4(4),colstr,qtype(2),gh
      integer minhel, maxhel, in,out,in2,out2,helcfg
      integer minhelq,maxhelq
      logical eval(4,-1:1), lpasstest(3)
      integer color_Z(2,2,4),color_C(2,2)
      data (color_Z(2,2,i),i=1,4) / 2, 2, 2, 2 /
      data (color_Z(1,2,i),i=1,4) / 2,-4,-2, 4 /
      data (color_Z(2,1,i),i=1,4) /-4, 2, 4,-2 /
      data (color_Z(1,1,i),i=1,4) /-4,-4,-4,-4 /
      data color_C(2,2), color_C(1,2), color_C(2,1), color_C(1,1) / 5, 3, 3, 1 /
      double precision color_m(8,8)

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
	pbar(mu,5)=p(mu,id(5),1)
	pdia(mu,1)=pbar(mu,1)*fsign(1)
	pdia(mu,2)=pbar(mu,2)*fsign(3)
	pdia(mu,3)=pbar(mu,3)*fsign(2)
	pdia(mu,4)=pbar(mu,4)*fsign(4)
	pdia(mu,5)=pbar(mu,5)*fsign(5)
      enddo      
      

      qtype(1) = nbos ! (W+Z)? u type : d type
      qtype(2) = 7-qtype(1)
      
c get bra and ket 
      tempi4(1) = fsign(1)
      tempi4(2) = fsign(3)
      tempi4(3) = fsign(2)
      tempi4(4) = fsign(4)
      call psi0m(4,pbar,tempi4,psi)
      do i = -1,1,2
      	call helvec(pbar(0,5),fsign(5),i,epsG(0,i))
      enddo
      do i=1,4
      do j=-1,1,2
	eval(i,j) = .true. ! to avoid reevaluating parts multiple times
      enddo;enddo
      gluWij(0,1,2) = 0d0
      gluWij(0,1,4) = 0d0
      gluWij(0,3,2) = 0d0
      gluWij(0,3,4) = 0d0
      
      do i=0,4
        res(i,3,0) = 0d0
        res(i,3,1) = 0d0
        res(i,4,0) = 0d0
        res(i,4,1) = 0d0
      enddo

      minhelq = max(qhel,1)
      maxhelq = qhel
      if (maxhelq.eq.0) maxhelq = 3

      if(ghel.eq.0) then
        minhel=-1
        maxhel=1
      else
        minhel=ghel
        maxhel=ghel
      endif

      do in = 1,3,2  ! determine q line coupling to W
      do out= 2,4,2
      do helcfg = minhelq,maxhelq
c new convention
c in   >W>>>>W>  out
c          S         , WZ
c	   S
c in2  >>>>>>>>  out2

c test for reasonable combination of helcfg, in, out
      in2 = 4-in  ! 2nd q index
      out2 =6-out
      lpasstest(helcfg) = .true.
      colstr = (in+out)/2
      if(colstr.eq.3  .or.  helcfg+colstr.eq.4)   lpasstest(helcfg) = .false.

      if(lpasstest(helcfg)) then
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
	  do gh = minhel,maxhel,2
	    call ket2c(psi(1,-1,in),.true.,pdia(0,in),-1,pDia(0,5),epsG(0,gh),psiG(1,-1,in,gh),mompsiG(0,in))

	    call ket2c(psiG(1,-1,in,gh),.false.,mompsiG(0,in),-1,momW(0,1),epsW(0,1),psiWG(1,in,gh),mompsiWG(0,in))
	    call ket2c(psiW(1,in),.false.,mompsiW(0,in),-1,pDia(0,5),epsG(0,gh),psiTemp(1,1),momTemp)
	    do i=1,2; psiWG(i,in,gh) = psiWG(i,in,gh) + psiTemp(i,1); enddo
	  enddo
	  eval(in,-1) = .false.
	endif

	if(eval(out,-1)) then ! calc <out,W|, <out,V|, <out,WV|
	  call bra2c(psi(1,-1,out),.true.,pdia(0,out),-1,momW(0,1),epsW(0,1),psiW(1,out),mompsiW(0,out))
	  do gh = minhel,maxhel,2
! 	   print*,'a'
	    call bra2c(psi(1,-1,out),.true.,pdia(0,out),-1,pDia(0,5),epsG(0,gh),psiG(1,-1,out,gh),mompsiG(0,out))

	    call bra2c(psiG(1,-1,out,gh),.false.,momPsiG(0,out),-1,momW(0,1),epsW(0,1),psiWG(1,out,gh),mompsiWG(0,out))
	    call bra2c(psiW(1,out),.false.,mompsiW(0,out),-1,pDia(0,5),epsG(0,gh),psiTemp(1,1),momTemp)
	    do i=1,2; psiWG(i,out,gh) = psiWG(i,out,gh) + psiTemp(i,1); enddo
	  enddo
	endif

	if(eval(in2,hel(in2))) then ! calc |V,in2> (and |W,in2>, |WV,in2> if needed for another comb. of in/out/helcfg)
         if(hel(in2).eq.-1) then
	    call ket2c(psi(1,-1,in2),.true.,pdia(0,in2),-1,momW(0,1),epsW(0,1),psiW(1,in2),mompsiW(0,in2))
	 endif
	 do gh = minhel,maxhel,2
          call ket2c(psi(1,hel(in2),in2),.true.,pdia(0,in2),hel(in2), pDia(0,5),epsG(0,gh),PsiG(1,hel(in2),in2,gh),momPsiG(0,in2))

	  if(hel(in2).eq.-1) then
	      call ket2c(psiG(1,-1,in2,gh),.false.,momPsiG(0,in2),-1,momW(0,1),epsW(0,1),psiWG(1,in2,gh),mompsiWG(0,in2))

	      call ket2c(psiW(1,in2),.false.,mompsiW(0,in2),-1,pDia(0,5),epsG(0,gh),psiTemp(1,1),momTemp)
	      do i=1,2; psiWG(i,in2,gh) = psiWG(i,in2,gh) + psiTemp(i,1); enddo
          endif
	 enddo
	  eval(in2,hel(in2)) = .false.
	endif

	if(eval(out2,hel(out2))) then ! calc |V,out2> (and |W,out2>, |WV,out2> if needed for another comb. of in/out/helcfg)
	  if(hel(out2).eq.-1) then
	    call bra2c(psi(1,-1,out2),.true.,pdia(0,out2),-1,momW(0,1),epsW(0,1),psiW(1,out2),mompsiW(0,out2))
	  endif
	  
	  do gh=minhel,maxhel,2
	    call bra2c(psi(1,hel(out2),out2),.true.,pdia(0,out2),hel(out2),pDia(0,5),epsG(0,gh),psiG(1,hel(out2),out2,gh),
     $                 mompsiG(0,out2))

	    if(hel(out2).eq.-1) then
	      call bra2c(psiG(1,-1,out2,gh),.false.,momPsiG(0,out2),-1,momW(0,1),epsW(0,1),psiWG(1,out2,gh),mompsiWG(0,out2))

	      call bra2c(psiW(1,out2),.false.,mompsiW(0,out2),-1,pDia(0,5),epsG(0,gh),psiTemp(1,1),momTemp)
	      do i=1,2; psiWG(i,out2,gh) = psiWG(i,out2,gh) + psiTemp(i,1); enddo
	          
	    endif
	  enddo
	endif
	endif !passtest
	enddo ! helcfg

c calculate glu currents
	if(gluWij(0,in,out).eq.0d0) then
c    W
	  call curr(-1,psi(1,-1,out),psiW(1,in),tempCurr(0,-1,1))
	  call curr(-1,psiW(1,out),psi(1,-1,in),tempCurr(0,-1,2))
	  do mu=0,3
	    gluWij(mu,in,out) =  tempCurr(mu,-1,1)+tempCurr(mu,-1,2)
	  enddo

	  do gh = minhel, maxhel,2
c    Wg
	    call curr(-1,psiG(1,-1,out,gh),psiW(1,in),tempCurr(0,-1,1))
	    call curr(-1,psiWG(1,out,gh),psi(1,-1,in),tempCurr(0,-1,2))
	    call curr(-1,psi(1,-1,out),psiWG(1,in,gh),tempCurr(0,-1,3))
	    call curr(-1,psiW(1,out),psiG(1,-1,in,gh),tempCurr(0,-1,4))
	    do mu=0,3
	      gluWGij(mu,in,out,gh,1) =  tempCurr(mu,-1,1)+tempCurr(mu,-1,2)
	      gluWGij(mu,in,out,gh,2) =  tempCurr(mu,-1,3)+tempCurr(mu,-1,4)
	    enddo
	  enddo

c lower line
c    q qbar
	  call curr(1,psi(1,-1,out2),psi(1,-1,in2),gluij(0,-1,in2,out2))

c    g
	  do gh =minhel,maxhel,2
	    call curr(1,psiG(1,-1,out2,gh),psi(1,-1,in2),gluGij(0,-1,in2,out2,gh,1))
	    call curr(1,psi(1,-1,out2),psiG(1,-1,in2,gh),gluGij(0,-1,in2,out2,gh,2))
	  enddo
	endif

c eventually calc amplitudes
      do helcfg = minhelq,maxhelq
      if(lpasstest(helcfg)) then
        hel(in) = -1
        hel(out)= -1
        if(helcfg.eq.1) then
          hel(in2) = -1
          hel(out2)= -1
	else
	  hel(in2) = 1
          hel(out2)= 1
	endif
      do gh=minhel,maxhel,2

	! crossed qlines only for specific in/out
	skip(qtype(1)) = out.eq.4
	skip(qtype(2)) =  in.eq.3

	! <WV.G> <.>   (. stands for uncontracted glu index (relative to ext. glu)
	temp = dotcc(gluij(0,hel(out2),in2,out2),gluWGij(0,in,out,gh,1))*G**3/qjj2(pdia(0,in2),pdia(0,out2))
! 	print*,temp
	if(.not.skip(3)) amp(out,colstr,3,helcfg,gh) = temp
	if(.not.skip(4)) amp(out,colstr,4,helcfg,gh) = temp
	! <WVG.> <.>   (. stands for uncontracted glu index (relative to ext. glu)
	temp = dotcc(gluij(0,hel(out2),in2,out2),gluWGij(0,in,out,gh,2))*G**3/qjj2(pdia(0,in2),pdia(0,out2))
! 	print*,temp
	if(.not.skip(3)) amp(in,colstr,3,helcfg,gh) = temp
	if(.not.skip(4)) amp(in,colstr,4,helcfg,gh) = temp
	! <WV.> <.G>   ( . stands for uncontracted glu index (relative to ext. glu)
	temp = G**3/qjj2(momPsiG(0,in2),pdia(0,out2))
	if(.not.skip(3)) amp(out2,colstr,3,helcfg,gh) = 0d0 + temp*dotcc(gluGij(0,hel(out2),in2,out2,gh,1),gluWij(0,in,out))
	if(.not.skip(4)) amp(out2,colstr,4,helcfg,gh) = 0d0 + temp*dotcc(gluGij(0,hel(out2),in2,out2,gh,1),gluWij(0,in,out))
! 	print*,temp*dotcc(gluGij(0,hel(out2),in2,out2,gh,1),gluWVij(0,in,out))
	! <WV.> <G.>   ( . stands for uncontracted glu index (relative to ext. glu)
	if(.not.skip(3)) amp(in2, colstr,3,helcfg,gh) = 0d0 + temp*dotcc(gluGij(0,hel(out2),in2,out2,gh,2),gluWij(0,in,out))
	if(.not.skip(4)) amp(in2, colstr,4,helcfg,gh) = 0d0 + temp*dotcc(gluGij(0,hel(out2),in2,out2,gh,2),gluWij(0,in,out))

c ggg Vertex contributions
	! <out|WV|in> glu <out2|in2>
	do mu = 0,3
	  momtemp(mu)= pdia(mu,out2)-pdia(mu,in2)
	enddo
	temp = G**3/mass2(momtemp)/qjj2(pdia(0,out),momPsiW(0,in))  
     &          *( dotrc(pdia(0,5),gluWij(0,in,out))*dotcc(epsG(0,gh),gluij(0,hel(out2),in2,out2))
     &            +dotrc(momtemp,  epsG(0,gh))*dotcc(gluWij(0,in,out),gluij(0,hel(out2),in2,out2))
     &            -dotrc(pdia(0,5),gluij(0,hel(out2),in2,out2))*dotcc(epsG(0,gh),gluWij(0,in,out)))
	if(.not.skip(3)) amp(out,colstr,3,helcfg,gh) = amp(out,colstr,3,helcfg,gh) - 2d0*temp
        if(.not.skip(4)) amp(out,colstr,4,helcfg,gh) = amp(out,colstr,4,helcfg,gh) - 2d0*temp
        if(.not.skip(3)) amp(in ,colstr,3,helcfg,gh) = amp(in, colstr,3,helcfg,gh) + 2d0*temp
        if(.not.skip(4)) amp(in ,colstr,4,helcfg,gh) = amp(in, colstr,4,helcfg,gh) + 2d0*temp
!         print*,temp

      enddo ! gluhel
      else ! not passtes
	if(colstr.le.2) then
	do i=1,4
! 	  print*,i,colstr,helcfg
	  amp(i,colstr,3,helcfg,-1) = 0d0
	  amp(i,colstr,4,helcfg,-1) = 0d0
	  amp(i,colstr,3,helcfg, 1) = 0d0
	  amp(i,colstr,4,helcfg, 1) = 0d0
	enddo
	endif
      endif ! not passtest
      enddo ! helcfg
      enddo; enddo ! in and out
c square amplitudes
      do i = 3,4 ! 2nd quark
	do j=0,4
	  res(j,i,0) = 0d0
	  res(j,i,1) = 0d0
	enddo
	do helcfg = minhelq,maxhelq
	do gh = minhel,maxhel,2
          if (ldoblha) then
            do m=1,4
c (1:4,1:4)
              color_m(m  ,m            ) =  blha_CA*blha_CF**2/2d0
              color_m(m  ,mod( 6-m,4)+1) = -blha_CA*blha_CF*(blha_CA-2d0*blha_CF)/4d0
              color_m(m  ,mod(m+ 1,4)+1) = -blha_CA*blha_CF*(blha_CA-2d0*blha_CF)/2d0
              color_m(m  ,5-m          ) = -blha_CA*blha_CF*(blha_CA-4d0*blha_CF)/4d0
c (5:8,5:8)
              color_m(m+4,m+4          ) =  blha_CA*blha_CF**2/2d0
              color_m(m+4,mod( 6-m,4)+5) = -blha_CA*blha_CF*(blha_CA-4d0*blha_CF)/4d0
              color_m(m+4,mod(m+ 1,4)+5) = -blha_CA*blha_CF*(blha_CA-2d0*blha_CF)/2d0
              color_m(m+4,9-m          ) = -blha_CA*blha_CF*(blha_CA-2d0*blha_CF)/4d0
c (1:4,5:8)
              color_m(m  ,m+4          ) =  blha_CA*blha_CF**2*(blha_CA-2d0*blha_CF)/2d0
              color_m(m  ,mod( 6-m,4)+5) = -blha_CF*(blha_CA-2d0*blha_CF)/4d0
              color_m(m  ,mod(m+ 1,4)+5) = -blha_CA*blha_CF*(blha_CA-2d0*blha_CF)*(blha_CA-blha_CF)/2d0
              color_m(m  ,9-m          ) = -blha_CF*(blha_CA-2d0*blha_CF)/4d0
            enddo
c upper right and lower left 4x4 blocks are identical
            color_m(5:8,1:4) = color_m(1:4,5:8)
            do k=1,4 
            do l=1,4
              res(0,i,1) = res(0,i,1) 
     &          + dreal( color_m(k,l)*
     &              amp(k,1,i,helcfg,gh)*dconjg(amp(l,1,i,helcfg,gh))
     &            )
            enddo
            enddo
            do k=1,8 
            do l=1,8
              res(0,i,0) = res(0,i,0) 
     &          + dreal( color_m(k,l)*
     &              amp(mod(k-1,4)+1,(k-1)/4+1,i,helcfg,gh)*dconjg(amp(mod(l-1,4)+1,(l-1)/4+1,i,helcfg,gh))
     &            )
            enddo
            enddo
          else
c for lhe output
          colsq(0,1) = 0d0
          colsq(0,2) = 0d0
          do k=1,4; do l=1,2
            colsq(k,l) = dble(amp(k,l,i,helcfg,gh)*dconjg(amp(k,l,i,helcfg,gh)))
          enddo; enddo
c    diff q families
	  do k=1,2; do l=1,2
c start with col 2
! ! ! ! ! 	    amp_ort = 0d0
! ! ! ! ! 	    do m=1,4
! ! ! ! ! 	      amp_ort = amp_ort + amp(m,2,i,helcfg,gh)*color_Z(k,l,m)
! ! ! ! ! 	    enddo
! ! ! ! ! 	    res(2,i,1) = res(2,i,1) + dble(amp_ort*dconjg(amp_ort))*color_c(k,l)/36d0
c reset and do col 1
	    amp_ort = 0d0
	    do m=1,4
	      amp_ort = amp_ort + amp(m,1,i,helcfg,gh)*color_Z(k,l,m)
	    enddo
	    tempR=dble(amp_ort*dconjg(amp_ort))*color_c(k,l)/36d0
            res(0,i,1) = res(0,i,1) + tempR
            colsq(0,1) = colsq(0,1) + tempR ! for lhe output
c    same q family
	    if(nbos.eq.i) then ! same initial q type
	      temp = (k*2)-3
	      amp_ort = amp_ort - temp*amp(1,2,i,helcfg,gh)*color_Z(k,l,3) ! change of initial state changes nnumbering of amplitudes!
	      amp_ort = amp_ort - temp*amp(2,2,i,helcfg,gh)*color_Z(k,l,2)
	      amp_ort = amp_ort - temp*amp(3,2,i,helcfg,gh)*color_Z(k,l,1)
	      amp_ort = amp_ort - temp*amp(4,2,i,helcfg,gh)*color_Z(k,l,4)
	    else ! same final q type
	      temp = (l*2)-3
	      amp_ort = amp_ort - temp*amp(1,2,i,helcfg,gh)*color_Z(k,l,1)
	      amp_ort = amp_ort - temp*amp(2,2,i,helcfg,gh)*color_Z(k,l,4)
	      amp_ort = amp_ort - temp*amp(3,2,i,helcfg,gh)*color_Z(k,l,3)
	      amp_ort = amp_ort - temp*amp(4,2,i,helcfg,gh)*color_Z(k,l,2)
	    endif
            tempR=dble(amp_ort*dconjg(amp_ort))*color_c(k,l)/36d0
            res(0,i,0) = res(0,i,0) + tempR
            colsq(0,2) = colsq(0,2) + tempR ! for lhe output
 	  enddo;enddo
c lhe output
          tempR=colsq(1,1)+colsq(2,1)+colsq(3,1)+colsq(4,1)
          if(tempR.ne.0d0) then
            res(1,i,1) = res(1,i,1)+ colsq(0,1)*(colsq(1,1)+colsq(4,1))/tempR
            res(2,i,1) = res(2,i,1)+ colsq(0,1)*(colsq(2,1)+colsq(3,1))/tempR
          endif
          tempR=tempR+colsq(1,2)+colsq(2,2)+colsq(3,2)+colsq(4,2)
          res(1,i,0) = res(1,i,0)+ colsq(0,2)*(colsq(1,1)+colsq(4,1))/tempR
          res(2,i,0) = res(2,i,0)+ colsq(0,2)*(colsq(2,1)+colsq(3,1))/tempR
          if(i.eq.3) then ! change particle 1,3 for colstr 2
            res(3,i,0)= res(3,i,0)+ colsq(0,2)*(colsq(3,2)+colsq(4,2))/tempR
            res(4,i,0)= res(4,i,0)+ colsq(0,2)*(colsq(1,2)+colsq(2,2))/tempR
          else ! change particle 1,3 for colstr 2
            res(3,i,0)= res(3,i,0)+ colsq(0,2)*(colsq(1,2)+colsq(2,2))/tempR
            res(4,i,0)= res(4,i,0)+ colsq(0,2)*(colsq(3,2)+colsq(4,2))/tempR
          endif
          endif !ldoblha
	enddo;enddo ! hel
      enddo
      if (ldoblha) then
        if (blha_amptype(blha_curproc).eq.0) then ! tree
          blha_amp(1)=blha_amp(1)+
     &      res(0,
     &          mod(blha_idsubproc(blha_cursubproc,blha_curproc),100)/10,
     &          mod(blha_idsubproc(blha_cursubproc,blha_curproc),10)    )
        else
          call BLHA_amptypeerror(blha_amptype(blha_curproc),__FILE__,__LINE__)
        endif
      endif

      end

