c ************************************************************************************
c Author: Matthias Kerner
c Date: 2012/13
c ************************************************************************************

c #define NINH_COMPARE 

c 4q WZ Amplitude
c    q1 >>> #### >>> q3 
c	    ####          , WZ
c	    ####
c    q2 >>> #### >>> q4  
c   qhel: el. (1,2,3) correspond to helicity configurations
c		1: all -
c		2: h1=h3=-1 , h2=h4= + 
c		3: h1!= h3 (and h2, h4 specified by charge flow of the specific subprocess)
c	  qhel= 0: sum hel. confg. 1-3
c   res(colstr,2nd qtype =3,4  ;  diff qflavor=0,1)
c      colstr=0..2: for catani seymour: <|>,  <|T1T2|>,  <|T1T3|>
c		    for lo: col summed, and the 2 different color flows
      subroutine amp4qwzCS(P,kin,id,fsign,nbos,qhel,res)

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

      
      real*8 pbar(0:3,4),pdia(0:3,4),momTemp(0:4)
      real*8 mompsiW(0:4,4),mompsiV(0:4,4),mompsiWV(0:4,4)
      complex*16 psi(2,-1:1,4), psiW(2,4),psiV(2,-1:1,3:4,4),psiWV(2,4)
      complex*16 psiTemp(2,-1:1),psiVV(2,4), amp(2,3:4,3) ! colstr, 2nd, hel
      complex*16 gluWij(0:3,1:3,2:4), gluWVij(0:3,1:3,2:4),temp
      complex*16 tempCurr(0:3,-1:1,4)
      complex*16 s1c
      real*8 qjj2
      real*8 squared(2,2)
      integer helswitch
      common/helswitchTemp/ helswitch

      integer i,j,mu,hel(4), tempi4(4), lepkin,colstr,qtype(2)
      integer minhel, maxhel, ikin, in,out,in2,out2,helcfg, lephel,lephel0
      logical eval(4,-1:1), lpasstest
      lnlo = .true.
      ikin = kin
      lephel = 0
      goto 11
      
      entry amp4qwzBorn(P,id,fsign,nbos,qhel,lephel0,res)
      lephel = lephel0
      ikin = 1
      lnlo= .false.
      
 11   continue
      lepkin = ikin
      if(ikin .gt.7) lepkin = 1

#ifdef NINH_COMPARE
      lephel = helswitch
cc      restemp = 0d0
cc      do lephel=-1,1,2
#endif
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
      
      if(nbos.gt.10) then
	qtype(1) = nbos/10 ! (W+Z)? u type : d type
	qtype(2) = 7-qtype(1)
      else
	qtype(1) = nbos
	qtype(2) = 7-qtype(1)
      endif
      
c get bra and ket 
      tempi4(1) = fsign(1)
      tempi4(2) = fsign(3)
      tempi4(3) = fsign(2)
      tempi4(4) = fsign(4)
      call psi0m(4,pbar,tempi4,psi)
      
      do i=1,4
      do j=-1,1,2
	eval(i,j) = .true. ! to avoid reevaluating parts multiple times
      enddo;enddo
      gluWij(0,1,2) = 0d0
      gluWij(0,1,4) = 0d0
      gluWij(0,3,2) = 0d0
      gluWij(0,3,4) = 0d0
      
      do i=1,3
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

c test for reasonable combination of helcfg, in, out
      lpasstest = .true.
      colstr = (in+out)/2
      if(helcfg .eq.1.and.(in+out.ne.7)) then
	hel(1) = -1; hel(2) = -1; hel(3)=-1; hel(4) = -1
      elseif(helcfg.eq.2.and.in+out.eq.3) then
	hel(1) = -1; hel(2) = -1; hel(3)= 1; hel(4) =  1
      elseif(helcfg.eq.3.and.(in+out.eq.5)) then
	hel(in) =-1; hel(out)=-1; hel(4-in)=1; hel(6-out) =1
      else
        lpasstest = .false.
      endif


      if(lpasstest) then
	in2 = 4-in  ! 2nd q index
	out2 =6-out
	
	if(eval(in,-1)) then ! calc |W,in>, |Z,in>, |WZ,in>
	  call ket2c(psi(1,-1,in),.true.,pdia(0,in),-1,momW(0,lepkin),epsW(0,lepkin),psiW(1,in),mompsiW(0,in))
	  call ket2c(psi(1,-1,in),.true.,pdia(0,in),-1,momWV(0,lepkin),epsWV(0,lephel,lepkin),psiWV(1,in),mompsiWV(0,in))
	  call ket2c(psi(1,-1,in),.true.,pdia(0,in),-1,momV(0,lepkin),epsV(0,lephel,qtype(1),-1,lepkin),psiV(1,-1,qtype(1),in),mompsiV(0,in))
	  call ket2c(psi(1,-1,in),.true.,pdia(0,in),-1,momV(0,lepkin),epsV(0,lephel,qtype(2),-1,lepkin),psiV(1,-1,qtype(2),in),mompsiV(0,in)) ! speed up: calc only for spec. helcfg
	  call ket2c(psiW(1,in),.false.,mompsiW(0,in),-1,momV(0,lepkin),epsV(0,lephel,qtype(2),-1,lepkin),psitemp(1,-1),momtemp)
	  call ket2c(psiV(1,-1,qtype(1),in),.false.,mompsiV(0,in),-1,momW(0,lepkin),epsW(0,lepkin),psitemp(1,1),momtemp)
	  psiVV(1,in) = psiTemp(1,-1)+psiTemp(1,1)-psiWV(1,in)
	  psiVV(2,in) = psiTemp(2,-1)+psiTemp(2,1)-psiWV(2,in)
	  eval(in,-1) = .false.
	endif
	if(eval(out,-1)) then ! calc <out,W|, <out,V|, <out,WV|
	  call bra2c(psi(1,-1,out),.true.,pdia(0,out),-1,momW(0,lepkin),epsW(0,lepkin),psiW(1,out),mompsiW(0,out))
	  call bra2c(psi(1,-1,out),.true.,pdia(0,out),-1,momWV(0,lepkin),epsWV(0,lephel,lepkin),psiWV(1,out),mompsiWV(0,out))
	  call bra2c(psi(1,-1,out),.true.,pdia(0,out),-1,momV(0,lepkin),epsV(0,lephel,qtype(2),-1,lepkin),psiV(1,-1,qtype(2),out),mompsiV(0,out))
	  call bra2c(psi(1,-1,out),.true.,pdia(0,out),-1,momV(0,lepkin),epsV(0,lephel,qtype(1),-1,lepkin),psiV(1,-1,qtype(1),out),mompsiV(0,out)) ! speed up: calc only for spec. helcfg
	  call bra2c(psiW(1,out),.false.,mompsiW(0,out),-1,momV(0,lepkin),epsV(0,lephel,qtype(1),-1,lepkin),psitemp(1,-1),momtemp)
	  call bra2c(psiV(1,-1,qtype(2),out),.false.,mompsiV(0,out),-1,momW(0,lepkin),epsW(0,lepkin),psitemp(1,1),momtemp)
	  psiVV(1,out) = psiTemp(1,-1)+psiTemp(1,1)-psiWV(1,out)
	  psiVV(2,out) = psiTemp(2,-1)+psiTemp(2,1)-psiWV(2,out)
	  eval(out,-1) = .false.
	endif
	if(eval(in2,hel(in2))) then ! calc |V,in2> (and |W,in2>, |WV,in2> if needed for another comb. of in/out/helcfg)
	  call ket2c(psi(1,hel(in2),in2),.true.,pdia(0,in2),hel(in2), momV(0,lepkin),epsV(0,lephel,3,hel(in2),lepkin),psiV(1,hel(in2),3,in2),mompsiV(0,in2))
	  call ket2c(psi(1,hel(in2),in2),.true.,pdia(0,in2),hel(in2), momV(0,lepkin),epsV(0,lephel,4,hel(in2),lepkin),psiV(1,hel(in2),4,in2),mompsiV(0,in2)) ! speed up possible if spec. helcfg
	  if(hel(in2).eq.-1) then
	    call ket2c(psi(1,-1,in2),.true.,pdia(0,in2),-1,momW(0,lepkin),epsW(0,lepkin),psiW(1,in2),mompsiW(0,in2))
	    call ket2c(psi(1,-1,in2),.true.,pdia(0,in2),-1,momWV(0,lepkin),epsWV(0,lephel,lepkin),psiWV(1,in2),mompsiWV(0,in2))
	    call ket2c(psiW(1,in2),.false.,mompsiW(0,in2),-1,momV(0,lepkin),epsV(0,lephel,qtype(2),-1,lepkin),psitemp(1,-1),momtemp)
	    call ket2c(psiV(1,-1,qtype(1),in2),.false.,mompsiV(0,in2),-1,momW(0,lepkin),epsW(0,lepkin),psitemp(1,1),momtemp)
	    psiVV(1,in2) = psiTemp(1,-1)+psiTemp(1,1)-psiWV(1,in2)
	    psiVV(2,in2) = psiTemp(2,-1)+psiTemp(2,1)-psiWV(2,in2)
	  endif
	  eval(in2,hel(in2)) = .false.
	endif
	if(eval(out2,hel(out2))) then ! calc |V,out2> (and |W,out2>, |WV,out2> if needed for another comb. of in/out/helcfg)
	  call bra2c(psi(1,hel(out2),out2),.true.,pdia(0,out2),hel(out2), momV(0,lepkin),epsV(0,lephel,3,hel(out2),lepkin),psiV(1,hel(out2),3,out2),mompsiV(0,out2))
	  call bra2c(psi(1,hel(out2),out2),.true.,pdia(0,out2),hel(out2), momV(0,lepkin),epsV(0,lephel,4,hel(out2),lepkin),psiV(1,hel(out2),4,out2),mompsiV(0,out2)) ! speed up possible if spec. helcfg
	  if(hel(out2).eq.-1) then
	    call bra2c(psi(1,-1,out2),.true.,pdia(0,out2),-1,momW(0,lepkin),epsW(0,lepkin),psiW(1,out2),mompsiW(0,out2))
	    call bra2c(psi(1,-1,out2),.true.,pdia(0,out2),-1,momWV(0,lepkin),epsWV(0,lephel,lepkin),psiWV(1,out2),mompsiWV(0,out2))
	    call bra2c(psiW(1,out2),.false.,mompsiW(0,out2),-1,momV(0,lepkin),epsV(0,lephel,qtype(1),-1,lepkin),psitemp(1,-1),momtemp)
	    call bra2c(psiV(1,-1,qtype(2),out2),.false.,mompsiV(0,out2),-1,momW(0,lepkin),epsW(0,lepkin),psitemp(1,1),momtemp)
	    psiVV(1,out2) = psiTemp(1,-1)+psiTemp(1,1)-psiWV(1,out2)
	    psiVV(2,out2) = psiTemp(2,-1)+psiTemp(2,1)-psiWV(2,out2)
	  endif
	  eval(out2,hel(out2)) = .false.
	endif

c calculate glu currents
	if(gluWij(0,in,out).eq.0d0) then
	  call curr(-1,psi(1,-1,out),psiVV(1,in),tempCurr(0,-1,1))
	  call curr(-1,psiVV(1,out),psi(1,-1,in),tempCurr(0,-1,2))
	  call curr(-1,psiW(1,out),psiV(1,-1,qtype(1),in),tempCurr(0,-1,3))
	  call curr(-1,psiV(1,-1,qtype(2),out),psiW(1,in),tempCurr(0,-1,4))
	  do mu=0,3
	    gluWVij(mu,in,out) = tempCurr(mu,-1,1)+tempCurr(mu,-1,2)+tempCurr(mu,-1,3)+tempCurr(mu,-1,4)
	  enddo
	  call curr(-1,psi(1,-1,out),psiW(1,in),tempCurr(0,-1,1))
	  call curr(-1,psiW(1,out),psi(1,-1,in),tempCurr(0,-1,2))
	  do mu=0,3
	    gluWij(mu,in,out) =  tempCurr(mu,-1,1)+tempCurr(mu,-1,2)
	  enddo
	endif
	
c eventually calc amplitudes
	temp = GG(1)**2/qjj2(pdia(0,in2),pdia(0,out2))
     &		* s1c(psi(1,hel(out2),out2),gluWVij(0,in,out),.true.,hel(out2),psi(1,hel(in2),in2))
        if(colstr.eq.1.or.in.eq.3) then ! 2nd = qtype(1)
	    amp(colstr,qtype(1),helcfg) = temp+ GG(1)**2/qjj2(pdia(0,out),mompsiW(0,in)) * (
     &		  s1c(psi(1,hel(out2),out2),gluWij(0,in,out),.true.,hel(out2),psiV(1,hel(in2),qtype(1),in2))
     &		 +s1c(psiV(1,hel(out2),qtype(1),out2),gluWij(0,in,out),.true.,hel(out2),psi(1,hel(in2),in2)))
	endif

        if(colstr.eq.1.or.out.eq.4) then ! 2nd = qtype(2)
	    amp(colstr,qtype(2),helcfg) = temp+ GG(1)**2/qjj2(pdia(0,out),mompsiW(0,in)) * (
     &		  s1c(psi(1,hel(out2),out2),gluWij(0,in,out),.true.,hel(out2),psiV(1,hel(in2),qtype(2),in2))
     &		 +s1c(psiV(1,hel(out2),qtype(2),out2),gluWij(0,in,out),.true.,hel(out2),psi(1,hel(in2),in2)))
	endif

! 	print*,amp(1,3,2)
      endif ! passtest
      enddo ! helcfg
      enddo; enddo ! in and out

c square amplitudes
      do i = 3,4 ! 2nd quark
	do j=0,2
	  res(j,i,0) = 0d0
	  res(j,i,1) = 0d0
	  res(j,i,2) = 0d0
	enddo
	do helcfg = minhel,maxhel
	  amp(2,i,helcfg) = - amp(2,i,helcfg) ! change of fermion indices
	  
	  squared(1,1) = dble(amp(1,i,helcfg)*dconjg(amp(1,i,helcfg)))
	  squared(2,2) = dble(amp(2,i,helcfg)*dconjg(amp(2,i,helcfg)))
	  squared(1,2) = dble(amp(1,i,helcfg)*dconjg(amp(2,i,helcfg)))

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

	    res(1,i,0) = res(1,i,0)  + squared(1,1)*2d0*(1d0-2d0/3d0*squared(1,2)/(squared(1,1)+squared(2,2)))
	    res(2,i,0) = res(0,i,0) - res(1,i,0)
	  endif
	enddo
	if(lnlo) then 
c 	    we swapped particle 1 and 2 --> res(::,2) is born, T21 and T23, --> we have to modify (2,:,2) again
c           or                  3 and 4 --> res(::,2) is born, T12 and T14
! 	    res(2,i,2) = -4d0/3d0*res(0,i,2)-res(1,i,2)-res(2,i,2)
	endif
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

#ifdef NINH_COMPARE
c      restemp=restemp+res
c      enddo
c      res=restemp
#endif
!       pause
      end


      subroutine amp4q1gwzNew(P,id,fsign,nbos,qhel,ghel,res)

      use globalvars, only: ldoblha

      implicit none
#include "common.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"
#include "brakets.inc"
      real*8 P(0:3,max_p,max_kin),res(0:4,3:4,0:1)
      integer id(5),fsign(5),nbos,qhel,ghel

      complex*16 epsW(0:3,7),epsWV(0:3,-1:1,7) 
      complex*16 epsV(0:3,-1:1,3:4,-1:1,7) ! ( , mu-hel, u/d , quark-hel), coupl. const to quarks included
      complex*16 epsA(0:3,-1:1),epsZ(0:3,-1:1)
      real*8 momV(0:4,7),momW(0:4,7),momWV(0:4,7)
      common /qcdwznew/ epsW, epsWV, epsV,epsA,epsZ,momV,momW,momWV

      
      real*8 pbar(0:3,5),pdia(0:3,5),momTemp(0:4)
!       real*8 mompsiW(0:4,4),mompsiV(0:4,4),mompsiWV(0:4,4)
!       real*8 mompsiWG(0:4,4),
      real*8 mompsiVG(0:4,4),mompsiWVG(0:4,4) !,mompsiGG(0:4,2,2),mompsiG(0:4,4)
!       complex*16 psi(2,-1:1,4), psiW(2,4),psiV(2,-1:1,3:4,4),psiWV(2,4)
!       complex*16 psiG(2,-1:1,4,-1:1), psiWG(2,4,-1:1),
      complex*16 psiVG(2,-1:1,3:4,4,2),psiWVG(2,4,2)
      complex*16 psiTemp(2,-1:1),psiVVG(2,4,2), amp(4,2,3:4,3,2) ! numbering as in MG u u -> e+ ve u d mu+ mu- g (1<>3, 2<>4 interchanged!), colstr, 2nd, hel
!       complex*16 epsG(0:3,-1:1)
      complex*16 gluWij(0:3,1:3,2:4), gluWVij(0:3,1:3,2:4),temp,temp2
      complex*16 gluWGij(0:3,1:3,2:4,2,2), gluWVGij(0:3,1:3,2:4,2,2)
      complex*16 gluij(0:3,-1:1,1:3,2:4),gluGij(0:3,-1:1,1:3,2:4,2,2)
      complex*16 gluVij(0:3,-1:1,3:4,1:3,2:4),gluVGij(0:3,-1:1,3:4,1:3,2:4,2,2)
      complex*16 tempCurr(0:3,-1:1,4),amp_ort
      complex*16 s1c, dotrc
      complex*16 dotcc
      external s1c, dotcc, dotrc,dotrr,qjj2,mass2, qjj2Q
      real*8 mass2,qjj2,qjj2Q
      real*8 tempR,colsq(0:4,2)
      real*8 dotrr
      logical skip(3:4)
      integer helswitch
      common/helswitchTemp/ helswitch

#ifdef NINH_COMPARE
      real*8 restemp(0:4,3:4,0:1)
#endif
      
      integer i,j,k,l,m,mu,hel(4),colstr,qtype(2),gh
      integer in,out,in2,out2,helcfg, lephel
      integer minhelq,maxhelq,idi,idi2,ido,ido2,idtemp(4)
      logical eval(4,-1:1), lpasstest(3)
      integer color_Z(2,2,4),color_C(2,2)
      data (color_Z(2,2,i),i=1,4) / 2, 2, 2, 2 /
      data (color_Z(1,2,i),i=1,4) / 2,-4,-2, 4 /
      data (color_Z(2,1,i),i=1,4) /-4, 2, 4,-2 /
      data (color_Z(1,1,i),i=1,4) /-4,-4,-4,-4 /
      data color_C(2,2), color_C(1,2), color_C(2,1), color_C(1,1) / 5, 3, 3, 1 /
      double precision color_m(8,8)

      lephel = 0

#ifdef NINH_COMPARE
c      restemp = 0d0
c      do lephel = -1,1,2
      lephel=helswitch
#endif

c new numbering convention
c in   q1 >>> #### >>> q2  out
c             ####          , WZ
c             ####
c in2  q3 >>> #### >>> q4  out2
!       do i=1,5
!       print*,(p(mu,i,1),mu=0,3)
!       enddo
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
      

      qtype(1) = nbos/10 ! (W+Z)? u type : d type
      qtype(2) = 7-qtype(1)
      
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
      
      idtemp(1)=id(1)
      idtemp(3)=id(2)
      idtemp(2)=id(3)
      idtemp(4)=id(4)

      minhelq = max(qhel,1)
      maxhelq = qhel
      if (maxhelq.eq.0) maxhelq = 3


      do in = 1,3,2  ! determine q line coupling to W
      do out= 2,4,2
      do helcfg = minhelq,maxhelq
c new convention
c in   >W>>>>W>  out
c          S         , WZ
c          S
c in2  >>>>>>>>  out2

c test for reasonable combination of helcfg, in, out
      in2 = 4-in  ! 2nd q index
      out2 =6-out
      idi=idtemp(in)
      ido=-idtemp(out)
      idi2=idtemp(in2)
      ido2=-idtemp(out2)
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
          do gh = 1,2
              !  |V,G,in>
            call ket2c(psiG(1,-1,idi,id(5),gh,1),.false.,mompsiG(0,idi,id(5),1),-1,momV(0,1),epsV(0,lephel,qtype(1),-1,1),
     $  psiVG(1,-1,qtype(1),in,gh),mompsiVG(0,in))
            call ket2c(psiG(1,-1,idi,id(5),gh,1),.false.,mompsiG(0,idi,id(5),1),-1,momV(0,1),epsV(0,lephel,qtype(2),-1,1),
     $   psiVG(1,-1,qtype(2),in,gh),mompsiVG(0,in)) ! speed up: calc only for spec. helcfg
            call ket2c(psiG(1,-1,idi,id(5),gh,1),.false.,mompsiG(0,idi,id(5),1),-1,momWV(0,1),epsWV(0,lephel,1),
     $   psiWVG(1,in,gh),mompsiWVG(0,in))

             !   |G,V,in> and add them
            call ket2r(psiV(1,-1,idi,qtype(1),1),.false.,mompsiV(0,idi,2,1),-1,pDia(0,5),epsG(0,id(5),gh,1),psiTemp(1,1),momTemp)
              do i=1,2; psiVG(i,-1,qtype(1),in,gh) = psiVG(i,-1,qtype(1),in,gh) + psiTemp(i,1); enddo
            call ket2r(psiV(1,-1,idi,qtype(2),1),.false.,mompsiV(0,idi,2,1),-1,pDia(0,5),epsG(0,id(5),gh,1),psiTemp(1,1),momTemp)
              do i=1,2; psiVG(i,-1,qtype(2),in,gh) = psiVG(i,-1,qtype(2),in,gh) + psiTemp(i,1); enddo

            call ket2c(psiWG(1,idi,id(5),gh,1),.false.,mompsiWG(0,idi,id(5),1),-1,momV(0,1),epsV(0,lephel,qtype(2),-1,1),
     $   psitemp(1,-1),momtemp)
            call ket2c(psiVG(1,-1,qtype(1),in,gh),.false.,mompsiVG(0,in),-1,momW(0,1),epsW(0,1),psitemp(1,1),momtemp)
            call ket2r(psiWV(1,idi,1),.false.,mompsiV(0,idi,3,1),-1,pDia(0,5),epsG(0,id(5),gh,1),psiTemp(1,0),momTemp)
              do i=1,2; psiVVG(i,in,gh) = psiTemp(i,-1)+psiTemp(i,0) + psiTemp(i,1)-psiWVG(i,in,gh); enddo          
          enddo
          eval(in,-1) = .false.
        endif
! 
        if(eval(out,-1)) then ! calc <out,W|, <out,V|, <out,WV|
          do gh = 1,2
            call bra2c(psiG(1,-1,ido,id(5),gh,1),.false.,momPsiG(0,ido,id(5),1),-1,momWV(0,1),epsWV(0,lephel,1),psiWVG(1,out,gh),
     $   mompsiWVG(0,out))
            call bra2c(psiG(1,-1,ido,id(5),gh,1),.false.,momPsiG(0,ido,id(5),1),-1,momV(0,1),epsV(0,lephel,qtype(2),-1,1),
     $   psiVG(1,-1,qtype(2),out,gh),mompsiVG(0,out))
            call bra2c(psiG(1,-1,ido,id(5),gh,1),.false.,momPsiG(0,ido,id(5),1),-1,momV(0,1),epsV(0,lephel,qtype(1),-1,1),
     $   psiVG(1,-1,qtype(1),out,gh),mompsiVG(0,out)) ! speed up: calc only for spec. helcfg
! ! 
            call bra2r(psiV(1,-1,ido,qtype(1),1),.false.,mompsiV(0,ido,2,1),-1,pDia(0,5),epsG(0,id(5),gh,1),psiTemp(1,1),momTemp)
              do i=1,2; psiVG(i,-1,qtype(1),out,gh) = psiVG(i,-1,qtype(1),out,gh) + psiTemp(i,1); enddo
            call bra2r(psiV(1,-1,ido,qtype(2),1),.false.,mompsiV(0,ido,2,1),-1,pDia(0,5),epsG(0,id(5),gh,1),psiTemp(1,1),momTemp)
              do i=1,2; psiVG(i,-1,qtype(2),out,gh) = psiVG(i,-1,qtype(2),out,gh) + psiTemp(i,1); enddo
! 
            call bra2c(psiWG(1,ido,id(5),gh,1),.false.,mompsiWG(0,ido,id(5),1),-1,momV(0,1),epsV(0,lephel,qtype(1),-1,1),
     $   psitemp(1,-1),momtemp)
            call bra2c(psiVG(1,-1,qtype(2),out,gh),.false.,mompsiVG(0,out),-1,momW(0,1),epsW(0,1),psitemp(1,1),momtemp)
            call bra2r(psiWV(1,ido,1),.false.,mompsiV(0,ido,3,1),-1,pDia(0,5),epsG(0,id(5),gh,1),psiTemp(1,0),momTemp)
              do i=1,2; psiVVG(i,out,gh) = psiTemp(i,-1)+psiTemp(i,0) + psiTemp(i,1)-psiWVG(i,out,gh); enddo        
          enddo
          eval(out,-1)=.false.
        endif
! 
        if(eval(in2,hel(in2))) then ! calc |V,in2> (and |W,in2>, |WV,in2> if needed for another comb. of in/out/helcfg)
          do gh = 1,2
            call ket2c(psiG(1,hel(in2),idi2,id(5),gh,1),.false.,momPsiG(0,idi2,id(5),1),hel(in2), momV(0,1),
     $   epsV(0,lephel,3,hel(in2),1),psiVG(1,hel(in2),3,in2,gh),mompsiVG(0,in2))
            call ket2c(psiG(1,hel(in2),idi2,id(5),gh,1),.false.,momPsiG(0,idi2,id(5),1),hel(in2), momV(0,1),
     $   epsV(0,lephel,4,hel(in2),1),psiVG(1,hel(in2),4,in2,gh),mompsiVG(0,in2)) ! speed up possible if spec. helcfg
! 
            call ket2r(psiV(1,hel(in2),idi2,3,1),.false.,mompsiV(0,idi2,2,1),hel(in2),pDia(0,5),epsG(0,id(5),gh,1),
     $   psiTemp(1,1),momTemp)
              do i=1,2; psiVG(i,hel(in2),3,in2,gh) = psiVG(i,hel(in2),3,in2,gh) + psiTemp(i,1); enddo
            call ket2r(psiV(1,hel(in2),idi2,4,1),.false.,mompsiV(0,idi2,2,1),hel(in2),pDia(0,5),epsG(0,id(5),gh,1),
     $   psiTemp(1,1),momTemp)
              do i=1,2; psiVG(i,hel(in2),4,in2,gh) = psiVG(i,hel(in2),4,in2,gh) + psiTemp(i,1); enddo

            if(hel(in2).eq.-1) then
              call ket2c(psiG(1,-1,idi2,id(5),gh,1),.false.,momPsiG(0,idi2,id(5),1),-1,momWV(0,1),epsWV(0,lephel,1),
     $   psiWVG(1,in2,gh),mompsiWVG(0,in2))
! 
              call ket2c(psiWG(1,idi2,id(5),gh,1),.false.,mompsiWG(0,idi2,id(5),1),-1,momV(0,1),epsV(0,lephel,qtype(2),-1,1),
     $   psitemp(1,-1),momtemp)
              call ket2c(psiVG(1,-1,qtype(1),in2,gh),.false.,mompsiVG(0,in2),-1,momW(0,1),epsW(0,1),psitemp(1,1),momtemp)
              call ket2r(psiWV(1,idi2,1),.false.,mompsiV(0,idi2,3,1),-1,pDia(0,5),epsG(0,id(5),gh,1),psiTemp(1,0),momtemp)
                do i=1,2; psiVVG(i,in2,gh) = psiTemp(i,-1)+psiTemp(i,0) + psiTemp(i,1)-psiWVG(i,in2,gh); enddo  
            endif
          enddo
          eval(in2,hel(in2)) = .false.
        endif
        if(eval(out2,hel(out2))) then ! calc |V,out2> (and |W,out2>, |WV,out2> if needed for another comb. of in/out/helcfg)
          do gh=1,2
            call bra2c(psiG(1,hel(out2),ido2,id(5),gh,1),.false.,momPsiG(0,ido2,id(5),1),hel(out2), momV(0,1),
     $   epsV(0,lephel,3,hel(out2),1),psiVG(1,hel(out2),3,out2,gh),mompsiVG(0,out2))
            call bra2c(psiG(1,hel(out2),ido2,id(5),gh,1),.false.,momPsiG(0,ido2,id(5),1),hel(out2), momV(0,1),
     $   epsV(0,lephel,4,hel(out2),1),psiVG(1,hel(out2),4,out2,gh),mompsiVG(0,out2)) ! speed up possible if spec. helcfg
            
            call bra2r(psiV(1,hel(out2),ido2,3,1),.false.,mompsiV(0,ido2,2,1),hel(out2),pDia(0,5),epsG(0,id(5),gh,1),
     $   psiTemp(1,1),momTemp)
              do i=1,2; psiVG(i,hel(out2),3,out2,gh) = psiVG(i,hel(out2),3,out2,gh) + psiTemp(i,1); enddo
            call bra2r(psiV(1,hel(out2),ido2,4,1),.false.,mompsiV(0,ido2,2,1),hel(out2),pDia(0,5),epsG(0,id(5),gh,1),
     $   psiTemp(1,1),momTemp)
              do i=1,2; psiVG(i,hel(out2),4,out2,gh) = psiVG(i,hel(out2),4,out2,gh) + psiTemp(i,1); enddo           
            
            if(hel(out2).eq.-1) then
              call bra2c(psiG(1,-1,ido2,id(5),gh,1),.false.,momPsiG(0,ido2,id(5),1),-1,momWV(0,1),epsWV(0,lephel,1),
     $   psiWVG(1,out2,gh),mompsiWVG(0,out2))
              
              call bra2c(psiWG(1,ido2,id(5),gh,1),.false.,mompsiWG(0,ido2,id(5),1),-1,momV(0,1),epsV(0,lephel,qtype(1),-1,1),
     $   psitemp(1,-1),momtemp)
              call bra2c(psiVG(1,-1,qtype(2),out2,gh),.false.,mompsiVG(0,out2),-1,momW(0,1),epsW(0,1),psitemp(1,1),momtemp)
              call bra2r(psiWV(1,ido2,1),.false.,mompsiV(0,ido2,3,1),-1,pDia(0,5),epsG(0,id(5),gh,1),psiTemp(1,0),momTemp)
                do i=1,2; psiVVG(i,out2,gh) = psiTemp(i,-1)+psiTemp(i,0) + psiTemp(i,1)-psiWVG(i,out2,gh); enddo            
            endif
          enddo
          eval(out2,hel(out2))=.false.
        endif
        endif !passtest
        enddo ! helcfg

c calculate glu currents
        if(gluWij(0,in,out).eq.0d0) then
c    WV upper line
          call curr(-1,psi(1,-1,ido,1),psiWV(1,idi,1),tempCurr(0,-1,1))
          call curr(-1,psiWV(1,ido,1),psi(1,-1,idi,1),tempCurr(0,-1,2))
          call curr(-1,psiW(1,ido,1),psiV(1,-1,idi,qtype(1),1),tempCurr(0,-1,3))
          call curr(-1,psiV(1,-1,ido,qtype(2),1),psiW(1,idi,1),tempCurr(0,-1,4))
          do mu=0,3
            gluWVij(mu,in,out) = tempCurr(mu,-1,1)+tempCurr(mu,-1,2)+tempCurr(mu,-1,3)+tempCurr(mu,-1,4)
          enddo

c    W
          call curr(-1,psi(1,-1,ido,1),psiW(1,idi,1),tempCurr(0,-1,1))
          call curr(-1,psiW(1,ido,1),psi(1,-1,idi,1),tempCurr(0,-1,2))
          do mu=0,3
            gluWij(mu,in,out) =  tempCurr(mu,-1,1)+tempCurr(mu,-1,2)
          enddo
c    WVg
          do gh = 1,2
c        g on the left
            call curr(-1,psiG(1,-1,ido,id(5),gh,1),psiWV(1,idi,1),tempCurr(0,-1,1))
            call curr(-1,psiVVG(1,out,gh),psi(1,-1,idi,1),tempCurr(0,-1,2))
            call curr(-1,psiWG(1,ido,id(5),gh,1),psiV(1,-1,idi,qtype(1),1),tempCurr(0,-1,3))
            call curr(-1,psiVG(1,-1,qtype(2),out,gh),psiW(1,idi,1),tempCurr(0,-1,4))
            do mu=0,3
              gluWVGij(mu,in,out,gh,1) = tempCurr(mu,-1,1)+tempCurr(mu,-1,2)+tempCurr(mu,-1,3)+tempCurr(mu,-1,4)
            enddo
! c        g on the right
            call curr(-1,psi(1,-1,ido,1),psiVVG(1,in,gh),tempCurr(0,-1,1))
            call curr(-1,psiWV(1,ido,1),psiG(1,-1,idi,id(5),gh,1),tempCurr(0,-1,2))
            call curr(-1,psiW(1,ido,1),psiVG(1,-1,qtype(1),in,gh),tempCurr(0,-1,3))
            call curr(-1,psiV(1,-1,ido,qtype(2),1),psiWG(1,idi,id(5),gh,1),tempCurr(0,-1,4))
            do mu=0,3
              gluWVGij(mu,in,out,gh,2) = tempCurr(mu,-1,1)+tempCurr(mu,-1,2)+tempCurr(mu,-1,3)+tempCurr(mu,-1,4)
            enddo 
c    Wg
            call curr(-1,psiG(1,-1,ido,id(5),gh,1),psiW(1,idi,1),tempCurr(0,-1,1))
            call curr(-1,psiWG(1,ido,id(5),gh,1),psi(1,-1,idi,1),tempCurr(0,-1,2))
            call curr(-1,psi(1,-1,ido,1),psiWG(1,idi,id(5),gh,1),tempCurr(0,-1,3))
            call curr(-1,psiW(1,ido,1),psiG(1,-1,idi,id(5),gh,1),tempCurr(0,-1,4))
            do mu=0,3
              gluWGij(mu,in,out,gh,1) =  tempCurr(mu,-1,1)+tempCurr(mu,-1,2)
              gluWGij(mu,in,out,gh,2) =  tempCurr(mu,-1,3)+tempCurr(mu,-1,4)
            enddo
          enddo
          
          
c lower line
c    q qbar
          call curr(1,psi(1,-1,ido2,1),psi(1,-1,idi2,1),gluij(0,-1,in2,out2))
c    V
          call curr(1,psi(1,-1,ido2,1),psiV(1,-1,idi2,3,1),tempCurr(0,-1,1))
          call curr(1,psiV(1,-1,ido2,3,1),psi(1,-1,idi2,1),tempCurr(0,-1,2))
          call curr(1,psi(1,-1,ido2,1),psiV(1,-1,idi2,4,1),tempCurr(0,-1,3))
          call curr(1,psiV(1,-1,ido2,4,1),psi(1,-1,idi2,1),tempCurr(0,-1,4))
          do mu=0,3; do i=-1,1,2
            gluVij(mu,i,3,in2,out2) = tempCurr(mu,i,1)+tempCurr(mu,i,2)
            gluVij(mu,i,4,in2,out2) = tempCurr(mu,i,3)+tempCurr(mu,i,4)
          enddo;enddo
c    g
          do gh =1,2
            call curr(1,psiG(1,-1,ido2,id(5),gh,1),psi(1,-1,idi2,1),gluGij(0,-1,in2,out2,gh,1))
            call curr(1,psi(1,-1,ido2,1),psiG(1,-1,idi2,id(5),gh,1),gluGij(0,-1,in2,out2,gh,2))
c    Vg     
            call curr(1,psiVG(1,-1,3,out2,gh),psi(1,-1,idi2,1),tempCurr(0,-1,1))
            call curr(1,psiG(1,-1,ido2,id(5),gh,1),psiV(1,-1,idi2,3,1),tempCurr(0,-1,2))
            call curr(1,psiVG(1,-1,4,out2,gh),psi(1,-1,idi2,1),tempCurr(0,-1,3))
            call curr(1,psiG(1,-1,ido2,id(5),gh,1),psiV(1,-1,idi2,4,1),tempCurr(0,-1,4))
            do mu=0,3; do i=-1,1,2
              gluVGij(mu,i,3,in2,out2,gh,1) = tempCurr(mu,i,1)+tempCurr(mu,i,2)
              gluVGij(mu,i,4,in2,out2,gh,1) = tempCurr(mu,i,3)+tempCurr(mu,i,4)
            enddo;enddo
            call curr(1,psiV(1,-1,ido2,3,1),psiG(1,-1,idi2,id(5),gh,1),tempCurr(0,-1,1))
            call curr(1,psi(1,-1,ido2,1),psiVG(1,-1,3,in2,gh),tempCurr(0,-1,2))
            call curr(1,psiV(1,-1,ido2,4,1),psiG(1,-1,idi2,id(5),gh,1),tempCurr(0,-1,3))
            call curr(1,psi(1,-1,ido2,1),psiVG(1,-1,4,in2,gh),tempCurr(0,-1,4))
            do mu=0,3; do i=-1,1,2
              gluVGij(mu,i,3,in2,out2,gh,2) = tempCurr(mu,i,1)+tempCurr(mu,i,2)
              gluVGij(mu,i,4,in2,out2,gh,2) = tempCurr(mu,i,3)+tempCurr(mu,i,4)
            enddo;enddo
          enddo
        endif
!       print*,'own'
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
      do gh=1,2

        ! crossed qlines only for specific in/out
        skip(qtype(1)) = out.eq.4
        skip(qtype(2)) =  in.eq.3

        ! <WV.G> <.>   (. stands for uncontracted glu index (relative to ext. glu)
        temp = dotcc(gluij(0,hel(out2),in2,out2),gluWVGij(0,in,out,gh,1))*G**3/qjj2(pdia(0,in2),pdia(0,out2))
!       print*,temp
        if(.not.skip(3)) amp(out,colstr,3,helcfg,gh) = temp
        if(.not.skip(4)) amp(out,colstr,4,helcfg,gh) = temp
        ! <WVG.> <.>   (. stands for uncontracted glu index (relative to ext. glu)
        temp = dotcc(gluij(0,hel(out2),in2,out2),gluWVGij(0,in,out,gh,2))*G**3/qjj2(pdia(0,in2),pdia(0,out2))
!       print*,temp
        if(.not.skip(3)) amp(in,colstr,3,helcfg,gh) = temp
        if(.not.skip(4)) amp(in,colstr,4,helcfg,gh) = temp
        ! <W.G> <.V>   ( . stands for uncontracted glu index (relative to ext. glu)
        temp = G**3/qjj2(momPsiV(0,idi2,2,1),pdia(0,out2))
        if(.not.skip(3)) amp(out,colstr,3,helcfg,gh) = amp(out,colstr,3,helcfg,gh) + temp*dotcc(gluVij(0,hel(out2),3,in2,out2),
     $   gluWGij(0,in,out,gh,1))
        if(.not.skip(4)) amp(out,colstr,4,helcfg,gh) = amp(out,colstr,4,helcfg,gh) + temp*dotcc(gluVij(0,hel(out2),4,in2,out2),
     $   gluWGij(0,in,out,gh,1))
!       print*,temp*dotcc(gluVij(0,hel(out2),3,in2,out2),gluWGij(0,in,out,gh,1))
        ! <WG.> <.V>   ( . stands for uncontracted glu index (relative to ext. glu)
        if(.not.skip(3)) amp(in, colstr,3,helcfg,gh) = amp(in ,colstr,3,helcfg,gh) + temp*dotcc(gluVij(0,hel(out2),3,in2,out2),
     $   gluWGij(0,in,out,gh,2))
        if(.not.skip(4)) amp(in, colstr,4,helcfg,gh) = amp(in ,colstr,4,helcfg,gh) + temp*dotcc(gluVij(0,hel(out2),4,in2,out2),
     $   gluWGij(0,in,out,gh,2))
!       print*,temp*dotcc(gluVij(0,hel(out2),4,in2,out2),gluWGij(0,in,out,gh,2))
        ! <WV.> <.G>   ( . stands for uncontracted glu index (relative to ext. glu)
        temp = G**3/qjj2(momPsiG(0,idi2,id(5),1),pdia(0,out2))
        if(.not.skip(3)) amp(out2,colstr,3,helcfg,gh) = 0d0 + temp*dotcc(gluGij(0,hel(out2),in2,out2,gh,1),gluWVij(0,in,out))
        if(.not.skip(4)) amp(out2,colstr,4,helcfg,gh) = 0d0 + temp*dotcc(gluGij(0,hel(out2),in2,out2,gh,1),gluWVij(0,in,out))
! !       print*,temp*dotcc(gluGij(0,hel(out2),in2,out2,gh,1),gluWVij(0,in,out))
        ! <WV.> <G.>   ( . stands for uncontracted glu index (relative to ext. glu)
        if(.not.skip(3)) amp(in2, colstr,3,helcfg,gh) = 0d0 + temp*dotcc(gluGij(0,hel(out2),in2,out2,gh,2),gluWVij(0,in,out))
        if(.not.skip(4)) amp(in2, colstr,4,helcfg,gh) = 0d0 + temp*dotcc(gluGij(0,hel(out2),in2,out2,gh,2),gluWVij(0,in,out))
! !       print*,temp*dotcc(gluGij(0,hel(out2),in2,out2,gh,2),gluWVij(0,in,out))
        ! <W.> <V.G>   ( . stands for uncontracted glu index (relative to ext. glu)
        temp = G**3/qjj2(momPsiVG(0,in2),pdia(0,out2))
        if(.not.skip(3)) amp(out2,colstr,3,helcfg,gh) = amp(out2,colstr,3,helcfg,gh) + 
     $   temp*dotcc(gluVGij(0,hel(out2),3,in2,out2,gh,1),gluWij(0,in,out))
        if(.not.skip(4)) amp(out2,colstr,4,helcfg,gh) = amp(out2,colstr,4,helcfg,gh) + 
     $   temp*dotcc(gluVGij(0,hel(out2),4,in2,out2,gh,1),gluWij(0,in,out))
!       print*,temp*dotcc(gluVGij(0,hel(out2),4,in2,out2,gh,1),gluWij(0,in,out))
        ! <W.> <VG.>   ( . stands for uncontracted glu index (relative to ext. glu)
        if(.not.skip(3)) amp(in2 ,colstr,3,helcfg,gh) = amp(in2 ,colstr,3,helcfg,gh) + 
     $   temp*dotcc(gluVGij(0,hel(out2),3,in2,out2,gh,2),gluWij(0,in,out))
        if(.not.skip(4)) amp(in2 ,colstr,4,helcfg,gh) = amp(in2 ,colstr,4,helcfg,gh) + 
     $   temp*dotcc(gluVGij(0,hel(out2),4,in2,out2,gh,2),gluWij(0,in,out))
!       print*,temp 
c ggg Vertex contributions
        ! <out|WV|in> glu <out2|in2>
        do mu = 0,3
          momtemp(mu)= pdia(mu,out2)-pdia(mu,in2)
        enddo
        temp = G**3/mass2(momtemp)/qjj2(momPsiV(0,ido,2,1),momPsiV(0,idi,1,1))  
     &          *( dotrc(pdia(0,5),gluWVij(0,in,out))*dotrc(epsG(0,id(5),gh,1),gluij(0,hel(out2),in2,out2))
     &            +dotrr(momtemp,  epsG(0,id(5),gh,1))*dotcc(gluWVij(0,in,out),gluij(0,hel(out2),in2,out2))
     &            -dotrc(pdia(0,5),gluij(0,hel(out2),in2,out2))*dotrc(epsG(0,id(5),gh,1),gluWVij(0,in,out)))
        if(.not.skip(3)) amp(out,colstr,3,helcfg,gh) = amp(out,colstr,3,helcfg,gh) - 2d0*temp
        if(.not.skip(4)) amp(out,colstr,4,helcfg,gh) = amp(out,colstr,4,helcfg,gh) - 2d0*temp
        if(.not.skip(3)) amp(in ,colstr,3,helcfg,gh) = amp(in, colstr,3,helcfg,gh) + 2d0*temp
        if(.not.skip(4)) amp(in ,colstr,4,helcfg,gh) = amp(in, colstr,4,helcfg,gh) + 2d0*temp
!         print*,temp
!         print*,temp
        ! <W.> glu <V.>
        do mu = 0,3
          momtemp(mu)= momPsiV(mu,ido2,2,1)-pdia(mu,in2)
        enddo
        temp = G**3/mass2(momtemp)/qjj2(pdia(0,out),momPsiV(0,idi,1,1))
        temp2 = temp*(
     &             dotrc(pdia(0,5),gluWij(0,in,out))*dotrc(epsG(0,id(5),gh,1),gluVij(0,hel(out2),3,in2,out2))
     &            +dotrr(momtemp,  epsG(0,id(5),gh,1))*dotcc(gluWij(0,in,out),gluVij(0,hel(out2),3,in2,out2))
     &            -dotrc(pdia(0,5),gluVij(0,hel(out2),3,in2,out2))*dotrc(epsG(0,id(5),gh,1),gluWij(0,in,out)))
        if(.not.skip(3)) amp(out,colstr,3,helcfg,gh) = amp(out,colstr,3,helcfg,gh) -2d0*temp2
        if(.not.skip(3)) amp(in, colstr,3,helcfg,gh) = amp(in, colstr,3,helcfg,gh) +2d0*temp2
!       print*,temp2
        temp2 = temp*(
     &             dotrc(pdia(0,5),gluWij(0,in,out))*dotrc(epsG(0,id(5),gh,1),gluVij(0,hel(out2),4,in2,out2))
     &            +dotrr(momtemp,  epsG(0,id(5),gh,1))*dotcc(gluWij(0,in,out),gluVij(0,hel(out2),4,in2,out2))
     &            -dotrc(pdia(0,5),gluVij(0,hel(out2),4,in2,out2))*dotrc(epsG(0,id(5),gh,1),gluWij(0,in,out)))
        if(.not.skip(4)) amp(out,colstr,4,helcfg,gh) = amp(out,colstr,4,helcfg,gh) -2d0*temp2
        if(.not.skip(4)) amp(in, colstr,4,helcfg,gh) = amp(in, colstr,4,helcfg,gh) +2d0*temp2

      enddo ! gluhel
      else ! not passtes
        if(colstr.le.2) then
        do i=1,4
!         print*,i,colstr,helcfg
          amp(i,colstr,3,helcfg,2) = 0d0
          amp(i,colstr,4,helcfg,2) = 0d0
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
        do gh = 1,2
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
            colsq(k,l) = amp(k,l,i,helcfg,gh)*dconjg(amp(k,l,i,helcfg,gh))
          enddo; enddo
c    diff q families
          do k=1,2; do l=1,2
! ! ! ! c start with col 2
! ! ! !             amp_ort = 0d0
! ! ! !             do m=1,4
! ! ! !               amp_ort = amp_ort + amp(m,2,i,helcfg,gh)*color_Z(k,l,m)
! ! ! !             enddo
! ! ! !             res(2,i,1) = res(2,i,1) + dble(amp_ort*dconjg(amp_ort))*color_c(k,l)/36d0
c reset and do col 1
            amp_ort = 0d0
            do m=1,4
              amp_ort = amp_ort + amp(m,1,i,helcfg,gh)*color_Z(k,l,m)
            enddo
            tempR=dble(amp_ort*dconjg(amp_ort))*color_c(k,l)/36d0
            res(0,i,1) = res(0,i,1) + tempR
            colsq(0,1) = colsq(0,1) + tempR ! for lhe output
c    same q family
            if(nbos/10.eq.i) then ! same initial q type
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

#ifdef NINH_COMPARE
cc      restemp = restemp + res
cc      enddo
cc      res = restemp
#endif
      end

       
