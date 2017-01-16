c ************************************************************************************
c Author: Matthias Kerner
c Date: 2013
c ************************************************************************************

! #define NINH_COMPARE
! #define CMPLX_MASS_SCHEME

      subroutine m2s_qcdzz3j(
     &                   xi,  !in:  Feynman x parameters of incoming fermions
     &                   p,   !in:  momenta associated with external fermions
     &                   v,     !in:  momenta associated with Z/H dec. fermions
     &                   nbos,
     &                   lnlo,
     &                   lokt,
     &                   xyuz,
     &                   rn,
     &                   m2s
     &                        )
      implicit none
c 
c declare input/output variables
c
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/amplitudes/vvjj/tensorz.inc"
#include "common.inc"
#include "cacheZZ.inc"
      integer nbos
      real*8 xi(nx), p(0:3,max_p,max_kin), v(0:3,max_v,max_kin), m2s(0:max_kin)
      real*8 xyuz(2,2:max_kin),rn(3)
      logical lokt(max_kin), lnlo
c
c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"
      double precision  alphas, xmt, alfa, xmz, xmw,s2w, xmh
      common /bkopin/   alphas, xmt, alfa, xmz, xmw, s2w, xmh
c
c declare local variables
c
      real*8 q_sf

      integer init/0/, I, J, k
      save init
c
c declare variables for summation process
c
      integer FSIGN(5+max_v) ! diagram index
      integer id(5) ! id(diagram index) = ps index
c
c store contributions from subprocess ip in res(ip) where
      real*8 res(maxNumSubProc)
      real*8 weight
c
c declare parton distribution variables
c
      real*8 pdf(-6:6,2,16)
c local varriables
      integer FL5qcdWZ, ires_min, ires_max, kin
      save ires_min, ires_max
      data ires_min,ires_max / 1, maxnumsubproc/
      external FL5qcdWZ

      logical ldebug,lallsub
!       parameter (ldebug=.false.)
      integer nhel
      real*8 mjj2, mass2
      integer hel, h, q
      integer qcdVVjjkin
      complex*16 phase1, phase2

      double precision  clr, xm2, xmg, b
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6)
      

      integer helswitch
      common/helswitchTemp/ helswitch


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
      

      ldebug = .false.

      if ( init .eq. 0 ) then
        if(nbos .eq. 22) then
         write(6,*) "QCD ZZ jjj -> 4l amplitude square information:"
        elseif(nbos.eq.21) then
         write(6,*) "QCD ZZ jjj -> 2l 2v amplitude square information:"
        elseif(nbos.eq.20) then
         write(6,*) "QCD ZA jjj -> 2l A amplitude square information:"
        elseif(nbos.eq.23) then
         write(6,*) "QCD ZA jjj -> 2v A amplitude square information:"
        elseif(nbos.eq.2) then
         write(6,*) "QCD Z jjj -> 2l amplitude square information:"
        elseif(nbos.eq.-2) then
         write(6,*) "QCD Z jjj -> 2v amplitude square information:"
        elseif (nbos.eq.24) then
         write(6,*) "QCD AA jjj amplitude square information:"
        endif
         if(nfl.eq.5) then
          write(6,*) "external b quarks are included"
         else
          write(6,*) "external b quarks are excluded"
         endif
         write(6,*) "alpha parameter of dipoles: ",alphaDip
         write(6,*) "-----------------------------------------------"
!          call CSGeneralSetAlpha(alphaDip)
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
         id(1)=1
         id(2)=2
         id(3)=3
         id(4)=4
         id(5)=5
         print*,'mjj(i,j) kins'
	do i = 1, 4
	do j = i+1,5
            print '(2I3, F20.9, 5I3)',i,j, sqrt(mjj2(p(0,i,1),p(0,j,1))), (qcdVVjjKin(id,i,j,k),k=1,5)
         end do; enddo
c
         write(6,*) " "
         write(6,*) "xi(1), xi(2) =", xi
         write(6,*) "mu_f =",sqrt(mufsq(1,1)),sqrt(mufsq(2,1))
         write(6,'(A,16I2)') ' lokt',(i,i=1,16)
         write(6,*) 'lokt',lokt
      endif

      lallsub = .true.

      do i=ires_min,ires_max
        res(i) = 0d0
      enddo


#ifdef NINH_COMPARE

      GG(1)= sqrt(als(1,1)*4d0*pi)
      GG(2) = GG(1)
      G = -GG(1)
      alphas = als(1,1)

      lokt=.true.
      p(0:3,1,1) = (/ 6730.08280683459d0       ,0.000000000000000d0  ,0.000000000000000d0, 6730.08280683459d0 /)
      p(0:3,2,1) = (/ 3947.20327739161d0     ,  0.000000000000000d0 , 0.000000000000000d0,  -3947.20327739161d0   /)
 
      p(0:3,3,1) = (/ 248.686154232159d0       ,-242.461245361708d0      ,-2.969293880494970d-14,   55.2932889637688d0     /)
      p(0:3,4,1) = (/ 800.241550046922d0       ,-8.75872413497419d0       ,-2.14219053547030d0     ,-800.190748629811d0     /)
      p(0:3,5,1) = (/ 6524.68199680939d0        ,57.5418907510555d0        ,71.7218259837599d0     ,   6524.03403347767d0     /)
      v(0:3,1,1) = (/401.085745435270d0        ,156.975925900645d0        ,99.8515562313491d0    ,-355.328018311502d0     /)
      v(0:3,2,1) = (/205.854205368964d0        ,75.7095603734415d0        ,103.784076996598d0     ,-160.850494864368d0     /)
      v(0:3,3,1) = (/303.513926899382d0       ,-7.87454107262725d0       ,-59.5491280591892d0     ,-297.510666652759d0     /)
      v(0:3,4,1) = (/2193.22250543412d0       ,-31.1328664558320d0       ,-213.666140617048d0     ,-2182.56786454002d0     /)
      xi(1) = p(0,1,1)/7000d0
      xi(2) = p(0,2,1)/7000d0
      do i=1,4
      do j=i,5
      print*,i,j,mjj2(p(0,i,1),p(0,j,1))/2d0
      enddo
      enddo
        call qcdVVJJptildemake(p,xyuz,v,lokt)
      pdf = 1d0
      
      do i=-4,4
        pdf(i,1,:) = 1d0+1d-2*i
        pdf(i,2,:) = 1d0+2d-2*i
      enddo

      xi(1) = 0.938549583219634d0
      xi(2) = 0.740985563038073d0

      q_sf = sqrt(mufsq(1,1))
      call pdfproton( xi(1), q_sf, pdf(-6,1,1) )        ! f_a(x1)

      q_sf = sqrt(mufsq(2,1))
      call pdfproton( xi(2), q_sf, pdf(-6,2,1) )        ! f_a(x2)

      pdf(-6:6,1,1) = pdf(-6:6,1,1)/xi(1)
      pdf(-6:6,2,1) = pdf(-6:6,2,1)/xi(2)

        do i =2,16
          if(lokt(i)) then
            if(.not.lokt(1).or.mufsq(1,1).ne.mufsq(1,i).or.mufsq(2,1).ne.mufsq(2,i)) then
              q_sf = sqrt(mufsq(1,i))
              call pdfproton( xi(1), q_sf, pdf(-6,1,i) )        ! f_a(x1)
              q_sf = sqrt(mufsq(2,i))
              call pdfproton( xi(2), q_sf, pdf(-6,2,i) )        ! f_a(x2)
              pdf(-6:6,1,i) = pdf(-6:6,1,i)/xi(1)
              pdf(-6:6,2,i) = pdf(-6:6,2,i)/xi(2)
            else
              do j = -6,6
                pdf(j,1,i) = pdf(j,1,1)
                pdf(j,2,i) = pdf(j,2,1)
              enddo
            endif
          endif
        enddo
#endif


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

      if(lnlo) then
	do i =2,16
	  if(lokt(i)) then
	    if(.not.lokt(1).or.mufsq(1,1).ne.mufsq(1,i).or.mufsq(2,1).ne.mufsq(2,i)) then
	      q_sf = sqrt(mufsq(1,i))
	      call pdfproton( xi(1), q_sf, pdf(-6,1,i) )        ! f_a(x1)
	      q_sf = sqrt(mufsq(2,i))
	      call pdfproton( xi(2), q_sf, pdf(-6,2,i) )        ! f_a(x2)
              do j=-5,5
                pdf(j,1,i) = pdf(j,1,i)/xi(1)
                pdf(j,2,i) = pdf(j,2,i)/xi(2)
              enddo
	    else
	      do j = -5,5
		pdf(j,1,i) = pdf(j,1,1)
		pdf(j,2,i) = pdf(j,2,1)
	      enddo
	    endif
	  endif
	enddo
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
!         print*,gzl
!         print*,gzlc
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
      

c calculate leptonic decay
      do kin=1,7
      if(kin.eq.1.or.lokt(kin)) then
      
      if(nbos.eq.24) then
        momV(0:3,1:2,kin) = v(:,1:2,kin)
        momV(0:3,3,kin) = momV(0:3,1,kin)+momV(0:3,2,kin)
        momV(4,1:2,kin) = 0d0
        momV(4,3,kin) = mass2(momV(0,3,kin))
      else
        momV(0:3,1,kin) = v(0:3,1,kin)+v(0:3,2,kin)
        momV(0:3,2,kin) = v(0:3,3,kin)+v(0:3,4,kin)
        momV(0:3,3,kin) = momV(0:3,1,kin)+momV(0:3,2,kin)

        momV(4,1,kin) = mass2(momV(0,1,kin))
        momV(4,2,kin) = mass2(momV(0,2,kin))
        momV(4,3,kin) = mass2(momV(0,3,kin))
      endif


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
c   lhel: 1..4,  V1 same for i=12/34; V2 same for i=13/24
c   qhel!=0: epsV including coupling to quark q
c   qhel=0: q=3: photon without coupling
c   qhel=0: q=4: Z without coupling

      CALL IXXXXX(v(0,1,kin),0d0 ,-ie,-1,lep)         !e+       
      CALL OXXXXX(v(0,2,kin),0d0 ,+ie, 1,lem)         !e- 

      CALL JIOXXX(lep,lem,GZL ,ZMASS,ZWIDTH,ze)     !Zl
      CALL JIOXXX(lep,lem,GAL ,0d0,0d0,ae)      !Al
      epsV(:,0,3,1,hel,kin) = ae(1:4)
      epsV(:,0,4,1,hel,kin) = ze(1:4)
      lAZ(:,1) = .true. ! A and Z cotribution for V1
      
      if(nbos.eq.22) then
        CALL IXXXXX(v(0,3,kin),0d0,-iu,-1,lup)         !mu+      
        CALL OXXXXX(v(0,4,kin),0d0 ,+iu, 1,lum)    

        CALL JIOXXX(lup,lum,GZL ,ZMASS,ZWIDTH,zu)     !Zl
        CALL JIOXXX(lup,lum,GAL ,0d0,0d0,au)      !Al
#ifdef WITH_VBF
        call vto4l(v(0,1,kin),hel,azz,zzztens)
#else
        print*,'please compile with vbf processes enabled'
        stop
#endif
        epsV(:,0,3,2,hel,kin) = au(1:4)
        epsV(:,0,4,2,hel,kin) = zu(1:4)
        epsV(:,0,3,3,hel,kin) = azz(0:3)
        epsV(:,0,4,3,hel,kin) = zzztens(0:3)
        lAZ(:,2:3) = .true. ! A and Z cotribution for V2
      elseif(nbos.eq.21) then
        CALL OXXXXX(v(0,3,kin),0d0 ,-1, 1,lum)         !vm
        CALL IXXXXX(v(0,4,kin),0d0 ,+1,-1,lup)         ! vm~   

        CALL JIOXXX(lup,lum,GZN ,ZMASS,ZWIDTH,zu)     !Zl
#ifdef WITH_VBF
        call vto4ln(v(0,1,kin),hel,azz,zzztens)
#else
        print*,'please compile with vbf processes enabled'
        stop
#endif
        epsV(:,0,3,2,hel,kin) = 0d0
        epsV(:,0,4,2,hel,kin) = zu(1:4)
        epsV(:,0,3,3,hel,kin) = azz(0:3)
        epsV(:,0,4,3,hel,kin) = zzztens(0:3)
        lAZ(:,2:3) = .true. 
        lAZ(1,2) = .false. ! no photon for V2
      elseif(nbos.eq.20) then   
        call VXXXXX(v(0,3,kin),0d0,iu,1,au) ! A

        epsV(:,0,3,2,hel,kin) = au(1:4)
        epsV(:,0,4,2,hel,kin) = 0d0
        
        CALL FVIXXX(lep,au,GAL ,0d0    ,0d0    ,lup)
        CALL FVOXXX(lem,au,GAL ,0d0    ,0d0,lum)

        CALL JIOXXX(lup,lem,GZL ,ZMASS,ZWIDTH,zu)     !Zl
        epsV(:,0,4,3,hel,kin) = zu(1:4)
        CALL JIOXXX(lep,lum,GZL ,ZMASS,ZWIDTH,zu)     !Zl
        epsV(:,0,4,3,hel,kin) = epsV(:,0,4,3,hel,kin)+zu(1:4)

        CALL JIOXXX(lup,lem,GAL ,0d0,0d0,zu)     !Zl
        epsV(:,0,3,3,hel,kin) = zu(1:4)
        CALL JIOXXX(lep,lum,GAL ,0d0,0d0,zu)     !Zl
        epsV(:,0,3,3,hel,kin) = epsV(:,0,3,3,hel,kin)+zu(1:4)

!         epsV(:,0,3:4,1:2,hel,1) = 0d0

        lAZ(:,2:3) = .true. 
        lAZ(2,2) = .false. ! no Z for V2
      elseif(nbos.eq.23) then   
        CALL OXXXXX(v(0,1,kin),0d0 ,-1, 1,lem)         !vm
        CALL IXXXXX(v(0,2,kin),0d0 ,+1,-1,lep)         ! vm~   
        CALL JIOXXX(lep,lem,GZN ,ZMASS,ZWIDTH,ze)     !Zn
        epsV(:,0,3,1,hel,kin) = 0d0
        epsV(:,0,4,1,hel,kin) = ze(1:4)
        
        call VXXXXX(v(0,3,kin),0d0,ie,1,au) ! A
        epsV(:,0,3,2,hel,kin) = au(1:4)
        epsV(:,0,4,2,hel,kin) = 0d0
        
        epsV(:,0,3:4,3,hel,kin) = 0d0 ! V->4l
        
        lAZ = .true.
        
c======
c try if this is working
!         lAZ=.false.
!         lAZ(2,1) = .true. ! Z contribution for V1
!         lAZ(1,2) = .true. ! A contribution for V2
c try if this is working
c======

      elseif(nbos.eq.24) then
        call VXXXXX(v(0,1,kin),0d0,ie,1,ae) 
        call VXXXXX(v(0,2,kin),0d0,iu,1,au) 
        
        epsV(:,0,3,1,hel,kin) = ae(1:4)
        epsV(:,0,3,2,hel,kin) = au(1:4)
        epsV(:,0,3,3,hel,kin) = 0d0
        
        lAZ = .true.

      elseif(abs(nbos).eq.2) then
        if(nbos.lt.0) then
        CALL OXXXXX(v(0,1,kin),0d0 ,-1, 1,lem)         !vm
        CALL IXXXXX(v(0,2,kin),0d0 ,+1,-1,lep)         ! vm~   
          CALL JIOXXX(lep,lem,GZN ,ZMASS,ZWIDTH,ze) ! Z-> nu nubar
          ae = 0d0
        endif
        epsV(:,0,3,3,hel,kin) = ae(1:4)
        epsV(:,0,4,3,hel,kin) = ze(1:4)
        lAZ(:,3) = .true. ! A and Z cotribution for V3, don't use V1 and V2
        lAZ(:,1:2) = .false.
      endif
c add A and Z current including couplings to u and d quarks
      do q=3,4
      do h=-1,1,2
      do i=1,3
        if(lAZ(0,i)) then
          epsV(:,h,q,i,hel,kin) = 0d0
#ifdef CMPLX_MASS_SCHEME
          if(lAZ(1,i)) epsV(:,h,q,i,hel,kin) = epsV(:,h,q,i,hel,kin) + epsV(:,0,3,i,hel,kin)*CLRc(q,1,h)
          if(lAZ(2,i)) epsV(:,h,q,i,hel,kin) = epsV(:,h,q,i,hel,kin) + epsV(:,0,4,i,hel,kin)*CLRc(q,2,h)
#else
          if(lAZ(1,i)) epsV(:,h,q,i,hel,kin) = epsV(:,h,q,i,hel,kin) + epsV(:,0,3,i,hel,kin)*CLR(q,1,h)
          if(lAZ(2,i)) epsV(:,h,q,i,hel,kin) = epsV(:,h,q,i,hel,kin) + epsV(:,0,4,i,hel,kin)*CLR(q,2,h)  
#endif
        endif
      enddo;enddo;enddo

      enddo
      
      phase1=exp((0d0,1d0)*pi*(rn(2)))
      phase2=exp((0d0,1d0)*pi*(rn(3)))

      if(nhel.eq.4) then
        epsV(:,:,:,1,0,kin) = phase1*epsV(:,:,:,1,1,kin)+conjg(phase1)*epsV(:,:,:,1,2,kin)
        epsV(:,:,:,2,0,kin) = phase2*epsV(:,:,:,2,1,kin)+conjg(phase2)*epsV(:,:,:,2,3,kin)
        epsV(:,:,:,3,0,kin) = phase1*(phase2*(epsV(:,:,:,3,1,kin))+conjg(phase2)*epsV(:,:,:,3,3,kin))
     &                     +conjg(phase1)*(phase2*(epsV(:,:,:,3,2,kin))+conjg(phase2)*epsV(:,:,:,3,4,kin))
      elseif(nhel.eq.2) then
         if(nbos.eq.21.or.nbos.eq.2) then
        epsV(:,:,:,1,0,kin) = phase1*epsV(:,:,:,1,1,kin)+conjg(phase1)*epsV(:,:,:,1,2,kin)
        epsV(:,:,:,2,0,kin) = epsV(:,:,:,2,2,kin)
        epsV(:,:,:,3,0,kin) = phase1*epsV(:,:,:,3,1,kin)+conjg(phase1)*epsV(:,:,:,3,2,kin)
        elseif(nbos.eq.23) then
        epsV(:,:,:,2,0,kin) = phase1*epsV(:,:,:,2,1,kin)+conjg(phase1)*epsV(:,:,:,2,2,kin)
        epsV(:,:,:,1,0,kin) = epsV(:,:,:,1,2,kin)
        epsV(:,:,:,3,0,kin) = 0d0
        else
          print*,__FILE__,__LINE__
          stop
        endif
      else ! Zn
        epsV(:,:,:,3,0,kin) = epsV(:,:,:,3,1,kin)
      endif


      endif
      enddo

cccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccc
ccc  Mad Graph
c     q1 g2 -> q3 g4 g5 wz  diagram
c     p1 p2    p3 p4 p5     phasespace
! 
!       lokt(1)=.true.
!       id(1) = 1
!       id(2) = 2
!       id(3) = 3
!       id(4) = 4
!       id(5) = 5
! 
!       fsign(1) =  1
!       fsign(2) = -1
!       fsign(3) =  1
!       fsign(4) =  1
!       fsign(5) =  1
!       
!       pMG(:,1:5) = p(:,1:5,1)
!       pMG(:,6:9) = v(:,1:4,1)
!       call SMATRIXUGUGG(pMG,res(1))
!       print*,res(1)
!       call precalcQCDZZjj(p,lokt,1,nbos)
!       call qg2qggZZ(p, id, fsign,1, amptemp)
! !       print*,amptemp(0,3),amptemp(0,3)/2
! !       pause
! !       return
!       res(1) = amptemp(0,3)
!       call precalcQCDZZjj(p,lokt,2,nbos)
!       call qg2qggZZ(p, id, fsign,2, amptemp)
!       res(1) = res(1)+amptemp(0,3)
!       call precalcQCDZZjj(p,lokt,3,nbos)
!       call qg2qggZZ(p, id, fsign,3, amptemp)
!       res(1) = res(1)+amptemp(0,3)
!       call precalcQCDZZjj(p,lokt,4,nbos)
!       call qg2qggZZ(p, id, fsign,4, amptemp)
!       res(1) = res(1)+amptemp(0,3)
!       print*,res(1)
!       pause
!       return
!


!       lokt(1)=.true.
!       id(1) = 1
!       id(2) = 2
!       id(3) = 3
!       id(4) = 4
!       id(5) = 5
! 
!       fsign(1) =  1
!       fsign(2) =  1
!       fsign(3) =  1
!       fsign(4) =  1
!       fsign(5) =  1
!       
!       pMG(:,1:5) = p(:,1:5,1)
!       pMG(:,6:9) = v(:,1:4,1)
! !       call SMATRIXUDUDGZZ(pMG,res(1))
!       call SMATRIXUUUUGZZ(pMG,res(1))
!       print*,res(1)
!       call precalcQCDZZjj(p,lokt,1,nbos)
!       call amp4q1gzz(p, id, fsign,1, amptemp4q)
! !       print*,amptemp4q(0,3,4,1)
! !       pause
! !       return
!       res(1) = amptemp4q(0,3,3,0)
!       call precalcQCDZZjj(p,lokt,2,nbos)
!       call amp4q1gzz(p, id, fsign,2, amptemp4q)
!       res(1) = res(1) + amptemp4q(0,3,3,0)
!       call precalcQCDZZjj(p,lokt,3,nbos)
!       call amp4q1gzz(p, id, fsign,3, amptemp4q)
!       res(1) = res(1) + amptemp4q(0,3,3,0)
!       call precalcQCDZZjj(p,lokt,4,nbos)
!       call amp4q1gzz(p, id, fsign,4, amptemp4q)
!       res(1) = res(1) + amptemp4q(0,3,3,0)
!       print*,res(1)
!       pause
!       return
      

!       lokt(1)=.true.
!       id(1) = 1
!       id(2) = 2
!       id(3) = 3
!       id(4) = 4
!       id(5) = 5
! 
!       fsign(1) =  1
!       fsign(2) =  1
!       fsign(3) =  1
!       fsign(4) =  1
!       fsign(5) =  1
!       call precalcQCDZZjj(p,lokt,0,nbos)
!       call amp4q1gzz(p, id, fsign,0, amptemp4q)
!       print*,amptemp4q(0,3,3,1),amptemp4q(0,4,4,1),amptemp4q(0,4,3,1)
!       p(:,:,2) = p(:,:,1); p(:,1,1)=p(:,2,2);p(:,2,1)=p(:,1,2) ! switch p1 p2
!       call precalcQCDZZjj(p,lokt,0,nbos)
!       call amp4q1gzz(p, id, fsign,0, amptemp4q)
!       print*,amptemp4q(0,3,3,2),amptemp4q(0,4,4,2),amptemp4q(0,4,3,2)
!       p(:,:,2) = p(:,:,1); p(:,3,1)=p(:,4,2);p(:,4,1)=p(:,3,2) ! switch p3 p4
!       call precalcQCDZZjj(p,lokt,0,nbos)
!       call amp4q1gzz(p, id, fsign,0, amptemp4q)
!       print*,amptemp4q(0,3,3,1),amptemp4q(0,4,4,1),amptemp4q(0,3,4,1)
!       
!       ! same quark families
!       call precalcQCDZZjj(p,lokt,0,nbos)
!       call amp4q1gzz(p, id, fsign,0, amptemp4q)
!       print*,amptemp4q(0,3,3,0),amptemp4q(0,4,4,0),amptemp4q(0,4,3,0)
!       p(:,:,2) = p(:,:,1); p(:,1,1)=p(:,2,2);p(:,2,1)=p(:,1,2) ! switch p1 p2
!       call precalcQCDZZjj(p,lokt,0,nbos)
!       call amp4q1gzz(p, id, fsign,0, amptemp4q)
!       print*,amptemp4q(0,3,3,0),amptemp4q(0,4,4,0),amptemp4q(0,4,3,0)
!       p(:,:,2) = p(:,:,1); p(:,3,1)=p(:,4,2);p(:,4,1)=p(:,3,2) ! switch p3 p4
!       call precalcQCDZZjj(p,lokt,0,nbos)
!       call amp4q1gzz(p, id, fsign,0, amptemp4q)
!       print*,amptemp4q(0,3,3,0),amptemp4q(0,4,4,0),amptemp4q(0,3,4,0)
!       pause
!       return


!  
!       if(sub_number.eq.2.or.lallsub) call getRes2Q3GZZ(p,id,fsign,pdf,nbos,lnlo,lokt,xyuz,res,v)
ccc  Mad Graph
cccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccc

      call precalcQCDZZjj(p,lokt,0)
      
      ires_min = FL5qcdWZ(id,id, -100000) +1

      
      call CSColAmpGNewReset()
      call CSColAmpQZZReset()

C--------begin subprocess evaluation-----------------------------------
C begin 2Q 3G Processes
c--crossing 1----------------------------------------------------------
c     g4 g2 -> qbar1 q3 g5 wz  diagram
c     p1 p2    p3   p4  p5     phasespace

      id(1) = 3
      id(2) = 2
      id(3) = 4
      id(4) = 1
      id(5) = 5

      fsign(1) = -1
      fsign(2) = -1
      fsign(3) =  1
      fsign(4) = -1
      fsign(5) =  1
      
      if(sub_number.eq.1.or.lallsub) call getRes2Q3GZZ(p,id,fsign,pdf,nbos,lnlo,lokt,xyuz,res)

c--crossing 2----------------------------------------------------------
c     q1 g2 -> q3 g4 g5 wz  diagram
c     p1 p2    p3 p4 p5     phasespace

      id(1) = 1
      id(2) = 2
      id(3) = 3
      id(4) = 4
      id(5) = 5

      fsign(1) =  1
      fsign(2) = -1
      fsign(3) =  1
      fsign(4) =  1
      fsign(5) =  1
 
      if(sub_number.eq.2.or.lallsub) call getRes2Q3GZZ(p,id,fsign,pdf,nbos,lnlo,lokt,xyuz,res)

c switch initial state
      id(1) = 2
      id(2) = 1

      if(sub_number.eq.3.or.lallsub) call getRes2Q3GZZ(p,id,fsign,pdf,nbos,lnlo,lokt,xyuz,res)

c--crossing 3----------------------------------------------------------
c     qb3 g2 -> qb1 g4 g5 wz  diagram
c     p1  p2     p3 p4 p5     phasespace

      id(1) = 3
      id(2) = 2
      id(3) = 1
      id(4) = 4
      id(5) = 5

      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = -1
      fsign(4) =  1
      fsign(5) =  1

      if(sub_number.eq.4.or.lallsub) call getRes2Q3GZZ(p,id,fsign,pdf,nbos,lnlo,lokt,xyuz,res)

c switch initial state
      id(3) = 2
      id(2) = 1

      if(sub_number.eq.5.or.lallsub) call getRes2Q3GZZ(p,id,fsign,pdf,nbos,lnlo,lokt,xyuz,res)

c--crossing 4----------------------------------------------------------
c     q1 qb3 -> g2 g4 g5 wz  diagram
c     p1  p2    p3 p4 p5     phasespace

      id(1) = 1
      id(2) = 3
      id(3) = 2
      id(4) = 4
      id(5) = 5

      fsign(1) =  1
      fsign(2) =  1
      fsign(3) = -1
      fsign(4) =  1
      fsign(5) =  1

      if(sub_number.eq.6.or.lallsub) call getRes2Q3GZZ(p,id,fsign,pdf,nbos,lnlo,lokt,xyuz,res)

c switch initial state
      id(1) = 2
      id(3) = 1

      if(sub_number.eq.7.or.lallsub) call getRes2Q3GZZ(p,id,fsign,pdf,nbos,lnlo,lokt,xyuz,res)

c
c
c  end of 2q 3g processes
c  begin 4q 1g processes
c
c
c--crossing 1----------------------------------------------------------
c     q1 q2 -> q3 q4 g zz  diagram
c     p1 p2    p3 p4 p5     phasespace

      id(1) = 1
      id(2) = 2
      id(3) = 3
      id(4) = 4
      id(5) = 5

      fsign(1) =  1
      fsign(2) =  1
      fsign(3) =  1
      fsign(4) =  1
      fsign(5) =  1
      
      if(sub_number.eq.8.or.lallsub) call getRes4Q1GZZ(p,id,fsign,pdf,nbos,lnlo,lokt,xyuz,res)
c
c--crossing 2----------------------------------------------------------
c     q1 qb3 ---> qb2 q4 g zz  diagram
c     p1 p2         p3 p4 p5     phasespace

      id(1) = 1
      id(2) = 3
      id(3) = 4
      id(4) = 2
      
      fsign(1) = 1
      fsign(2) = -1
      fsign(3) = 1
      fsign(4) = -1
      
      if(sub_number.eq.9.or.lallsub) call getRes4Q1GZZ(p,id,fsign,pdf,nbos,lnlo,lokt,xyuz,res)

      id(4) = 1
      id(1) = 2
      if(sub_number.eq.10.or.lallsub) call getRes4Q1GZZ(p,id,fsign,pdf,nbos,lnlo,lokt,xyuz,res)

c--crossing 3----------------------------------------------------------
c     qb3 qb4 ---> qb1 qb2 g wz  diagram
c     p1  p2    p3 p4 p5     phasespace

      id(1) = 3
      id(2) = 4
      id(3) = 1
      id(4) = 2
      
      fsign(1:4) = -1

      if(sub_number.eq.11.or.lallsub) call getRes4Q1GZZ(p,id,fsign,pdf,nbos,lnlo,lokt,xyuz,res)

c--crossing 4----------------------------------------------------------
c     g5 q1 -> qb2 q4 q3 wz  diagram
c     p1 p2    p3 p4  p5     phasespace

      id(1) = 2
      id(2) = 3
      id(3) = 5
      id(4) = 4
      id(5) = 1

c  Change momnentum convention:
      fsign(1) =  1
      fsign(2) = -1
      fsign(3) =  1
      fsign(4) =  1
      fsign(5) = -1

      if(sub_number.eq.12.or.lallsub) call getRes4Q1GZZ(p,id,fsign,pdf,nbos,lnlo,lokt,xyuz,res)

      id(5) = 2
      id(1) = 1

      if(sub_number.eq.13.or.lallsub) call getRes4Q1GZZ(p,id,fsign,pdf,nbos,lnlo,lokt,xyuz,res)

c--crossing 5----------------------------------------------------------
c     g5 qb3 -> qb2 qb1 q4 wz  diagram
c     p1  p2    p3  p4  p5     phasespace

      id(1) = 4
      id(2) = 3
      id(3) = 2
      id(4) = 5
      id(5) = 1

      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = -1
      fsign(4) =  1
      fsign(5) = -1

      if(sub_number.eq.14.or.lallsub) call getRes4Q1GZZ(p,id,fsign,pdf,nbos,lnlo,lokt,xyuz,res)

      id(5) = 2
      id(3) = 1

      if(sub_number.eq.15.or.lallsub) call getRes4Q1GZZ(p,id,fsign,pdf,nbos,lnlo,lokt,xyuz,res)
c
c
c---end of subprocess evaluation-------------------------------
c
c

      if(lnlo) then
	ires_max = 16
	ires_min = 1
      else
	ires_max = numcolor
      endif
      if(init.eq.0) then
	 if(.not.lnlo) print*, "Number of (color) subprocesses: ", ires_max
	init =1
      endif
      
      m2s(0) = 0d0
      do i=ires_min,ires_max
	if(lnlo) then
	  res(i) = res(i)*(als(1,i)/alphas)**3
	  m2s(i) = res(i)
	else
	  res(i) = res(i)*(als(1,1)/alphas)**3
	endif
	m2s(0)=m2s(0) + res(i)
      enddo
!       print*,res
c select subprocess
      iprocess =0
      weight = 0d0
       
      i=ires_min
      do while((i.le.ires_max).and.(iprocess.eq.0))
         weight = weight+res(i)/m2s(0)
         if(weight.ge.rn(1)) then
            iprocess = i
         endif
         i = i+1
      enddo

      if(iprocess.eq.0..and..not.lnlo) then
!          print*, 'error in subprocess selection'
      endif
!       ldebug = ldebug.or. abs(max(m2s(0),m2s(1))/min(m2s(0),m2s(1))-1d0).lt.1d-1
      
      if(ldebug ) then
        write(6,*) " "
         do i = 1, n_v
            write(6,20) i, v(0,i,1), v(1,i,1), v(2,i,1), v(3,i,1)
         end do
         print*," parton momenta "
         do i = 1, n_p
            write(6,10) i, p(0,i,1), p(1,i,1), p(2,i,1), p(3,i,1)
         end do
         id(1)=1
         id(2)=2
         id(3)=3
         id(4)=4
         id(5)=5
         print*,'mjj2(i,j) kins'
	do i = 1, 4
	do j = i+1,5
            print*,i,j, mjj2(p(0,i,1),p(0,j,1)), (qcdVVjjKin(id,i,j,k),k=1,5)
         end do; enddo
c
         write(6,*) " "
         write(6,*) "xi(1), xi(2) =", xi
         write(6,*) "mu_f =",sqrt(mufsq(1,1)),sqrt(mufsq(2,1))
	do i=0,16
	 print*,'m2s( ',i,' ) ', m2s(i), lokt(max(i,1))
	enddo
      endif
!       pause

#ifdef NINH_COMPARE
      m2s = m2s
      print*,m2s(1),m2s(0)-m2s(1),m2s(0)
#endif

      end

      subroutine getRes2Q3GZZ(p,id,fsign,pdf,Nbos,lnlo,lokt,xyuz, res)
      implicit none
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
      real*8 res(1:maxnumsubproc) ! lnlo ? res(kin) : res(subproc)
      real*8 p(0:3,max_p,max_kin),pdf(-6:6,2,16),xyuz(2,2:max_kin)
      integer fsign(5), id(5), Nbos, idb(4),fsignB(4)
      logical lnlo, lokt(16)
      
      real*8 amp(0:6,3:4) ! m^2: lnlo ? (Born, T1T2, T1T3) : les houches color flow
      integer flav(5), k, colstr,q,qt
      integer FL5QCDWZ
      real*8 polcol(maxnumsubproc), dip(3:4,4,5,5),dipSum(16,3:4)
      common / polcolz / polcol
      integer qcdVVjjKin
      integer emit1, emit2, spect, kin, i,j,spectB,emitB

      
      logical ldebug
      parameter (ldebug = .false.)
      external CScolAmpGNewUsePrecalc
      integer di

      double precision  alphas, xmt, alfa, xmz, xmw,s2w, xmh
      common /bkopin/   alphas, xmt, alfa, xmz, xmw, s2w, xmh
      
#ifdef NINH_COMPARE
      integer helswitch
      common/helswitchTemp/ helswitch
#endif

   
      flav(id(2)) = 0 ! Feyn. Daigram labels  2, 4, 5 are gluons
      flav(id(4)) = 0
      flav(id(5)) = 0
      
c get 3 jet amplitude
      if(lokt(1)) then
        call qg2qggZZ(p, id, fsign,0, amp)
      else
	  amp(0,:) = 0d0
      endif
#ifdef WITH_NLO
      if(lnlo) then
c calculate dipoles
      do i=1,16
        dipsum(i,3)=0d0
        dipsum(i,4)=0d0
      enddo
c qg splitting
      do emit1 = 1,3,2 ! qg splitting
      do emit2 = 2,5
      if(emit2.ne.3) then
      do spect  = 1,5
      if(spect.ne.emit1.and.spect.ne.emit2) then
        kin = qcdVVjjKin(id,emit1,emit2,spect)
        if(kin.gt.1)then
        if(lokt(kin)) then
          if(ldebug) print*, 'calc Dipole (',id(emit1),id(emit2),')',id(spect),'with kin',kin
          ! set up born
          do i=1,4
            idB(i) = id(i)
          enddo
          spectB=spect
          if(id(emit2).le.2) call swapInt(idB(emit1),idB(emit2)) ! i g->qq splitting, swap q with one g
          if(idB(2).ge.idB(4)) then
             call swapInt(idB(2),idB(4)) ! sort glu for borncache
             if(mod(spectB,2).eq.0) spectB=6-spectB
          endif
          if (spectB.eq.5) spectB=4
          if(id(spect).le.2.and.id(emit1).le.3.and.id(emit2).le.3) then
            ! for ii kin: order of final state fixed by RE order; by
            ! now, we removed particle 5, instead of 3; swapping 3 and 4
            ! will fix this
            i=di(idB,4)
            j=di(idB,3)
            idB(i)=3
            idB(j)=4
          endif

          if(id(5).ne.5) print*,__FILE__,__LINE__
          ! get dipole
          call CSgeneral(p,kin,xyuz,"qg",id(emit1),id(emit2),id(spect),0,0,
     &                     CScolAmpGNewUsePrecalc,idB,fsignB,emit1,spectB,nbos,2,dip(3,emit1,emit2,spect))
          dipsum(kin,3)=dipsum(kin,3)+dip(3,emit1,emit2,spect)
          dipsum(kin,4)=dipsum(kin,4)+dip(4,emit1,emit2,spect)
          if(ldebug) print*, 'Dipole: ',dip(3:4,emit1,emit2,spect)*alphas
        endif
        endif
      endif
      enddo
      endif
      enddo;enddo
c gg splitting
      do emit1=2,4,2
      do emit2=4,5
      do spect=1,5
      if(emit1.ne.emit2.and.emit1.ne.spect.and.emit2.ne.spect) then
        kin = qcdVVjjKin(id,emit1,emit2,spect)
        if(kin.gt.1) then
        if(lokt(kin)) then
          if(ldebug) print*, 'calc Dipole (',id(emit1),id(emit2),')',id(spect),'with kin',kin
          ! set up born
          do i=1,4
            idB(i) = id(i)
          enddo
          spectB=spect
          emitB=emit1
          if(idB(2).ge.idB(4)) then
             call swapInt(idB(2),idB(4)) ! sort glu for borncache
             if(mod(spectB,2).eq.0) spectB=6-spectB
             emitB=6-emitB
          endif
          if (spectB.eq.5) spectB=4
          call CSgeneral(p,kin,xyuz,"gg",id(emit1),id(emit2),id(spect),0,0,
     &                     CScolAmpGNewUsePrecalc,idB,fsignB,emitB,spectB,nbos,2,dip(3,emit1,emit2,spect))
          dipsum(kin,3)=dipsum(kin,3)+dip(3,emit1,emit2,spect)
          dipsum(kin,4)=dipsum(kin,4)+dip(4,emit1,emit2,spect)
          if(ldebug) print*, 'Dipole: ',dip(:,emit1,emit2,spect)*alphas
        endif
        endif
      endif
      enddo;enddo;enddo

      endif
#endif

c calculate subprocess
	do q = 1,nfl ! sum over quarks
          flav(id(1)) = q*fsign(1)  ! up type ( or anti up)
          flav(id(3)) = q*fsign(3) 
          
          qt = 3+mod(q,2) ! 3,4 for u,d

	  if (ldebug) then
	    print*
	    print*, 'new subprocess'
	    print*, 'partons:    ',flav
	    print*, 'di          ',(di(id,i),i=1,5)
	    print*, 'id          ',id
	    print*, 'lokt        ',lokt
	  endif
c NLO CS Subtraction
	  if(lnlo) then
#ifdef WITH_NLO
        k = FL5QCDWZ(flav,id,colstr) 
        if(lokt(1)) then
          res(1) = res(1) + amp(0,qt)*polcol(k)*
     1		pdf(sign1*flav(1),1,1)*pdf(sign2*flav(2),2,1)
          if(ldebug) print*,"|M_m+1|^2", amp(0,qt)
        endif
        do kin=2,16
          res(kin) = res(kin) - dipsum(kin,qt)*polcol(k)*alphas*
     1		pdf(sign1*flav(1),1,kin)*pdf(sign2*flav(2),2,kin)
          if(ldebug) print*,"Dipole   ", dipsum(kin,qt)*polcol(k)*alphas*
     1          pdf(sign1*flav(1),1,kin)*pdf(sign2*flav(2),2,kin)
        enddo
#endif
	  else
c LO subprocesses
	    do colstr = 1,6
	      k = FL5QCDWZ(flav,id,colstr) ! get process number, calc polcol and fill les houches data
! 	      print*,pdf(sign1*flav(1),1,1),pdf(sign2*flav(2),2,1)
	      res(k) = amp(colstr,qt)*POLCOL(k)*
     1		pdf(sign1*flav(1),1,1)*pdf(sign2*flav(2),2,1)
	    enddo
!             print*,flav
	  endif
	enddo ! Quark family
      end


      subroutine getRes4Q1GZZ(p,id,fsign,pdf,Nbos,lnlo,lokt,xyuz,res)
c hel summation, if rn < 0
      implicit none
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/global.inc"
#include "common.inc"
      real*8 res(1:maxNumSubProc) !out: subprocess cross-sections
      real*8 p(0:3,max_p,max_kin)
      integer fsign(5), id(5),fsignB(5)
      real*8 pdf(-6:6,2,16), xyuz(2,2:max_kin)
      logical lnlo,lokt(16)

      integer idB(4),emitB,spectB
      integer nbos

      integer FL5qcdWZ, di


      double precision  alphas, xmt, alfa, xmz, xmw,s2w, xmh
      common /bkopin/   alphas, xmt, alfa, xmz, xmw, s2w, xmh

      integer j,i,k,n,q1,q2,q,tempi
      integer iflav(1:5)

      integer kin,emit1,emit2,spect, qcdVVjjkin
      logical ldebug, idFl
      parameter (ldebug = .false.)
      real*8 amp(0:4,3:4,3:4,0:2),dip(3:4,3:4,0:2),dipsum(3:4,3:4,0:2,2:16),dipq(3:4)
      external CScolAmpGNewUsePrecalc,CScolAmpQZZ

      double precision polcol(maxnumsubproc)

      common /polcolz/ polcol

      if(ldebug) then
	print*,'new subprocess'
	print*,(di(id,i),i=1,5)
      endif
c precalc amplitude and dipoles
      if(lokt(1)) then
        call amp4q1gzz(p,id,fsign,0,amp)
      else
        do i=0,4; do j=3,4; do k=3,4
          amp(i,j,k,0) = 0d0
          amp(i,j,k,1) = 0d0
          amp(i,j,k,2) = 0d0
        enddo; enddo; enddo
      endif
      dipSum=0d0

#ifdef WITH_NLO
c qq dipoles
      if(lnlo) then
        do emit1=1,2; do emit2=3,4; do spect=1,5
          kin = qcdVVjjKin(id,emit1,emit2,spect)
          if(kin.gt.0) then
            if(lokt(kin)) then
              idB(1)=id(3-emit1)
              idB(3)=id(7-emit2)
              idB(2)=min(id(emit1),id(emit2))

              i=max(id(emit1),id(emit2))
              if(idB(1).gt.i) idB(1)=idB(1)-1
              if(idB(3).gt.i) idB(3)=idB(3)-1

              idB(4)=10-idB(1)-idB(2)-idB(3)
              emitB=2
              if(spect.le.2) then
                spectB=1
              elseif(spect.le.4) then
                spectB=3
              else
                spectB=4
              endif
              if(idB(2).ge.idB(4)) then
                call swapInt(idB(2),idB(4))
                if(spectB.eq.4) spectB=2
                emitB = 4
              endif
              call CSGeneral(p,kin,xyuz,"qq",id(emit1),id(emit2),id(spect),0,0,
     $                       CScolAmpGNewUsePrecalc,idB,fsignB,emitB,spectB,nbos,2,dipq(3))
              if(emit1.eq.1.and.emit2.eq.3) then
                do q=3,4; do i=0,1
                  dipsum(q,3,i,kin) = dipsum(q,3,i,kin) + alphas*dipq(3)
                  dipsum(q,4,i,kin) = dipsum(q,4,i,kin) + alphas*dipq(4)
                enddo; enddo
              elseif(emit1.eq.2.and.emit2.eq.4) then
                do q=3,4; do i=0,1
                  dipsum(3,q,i,kin) = dipsum(3,q,i,kin) + alphas*dipq(3)
                  dipsum(4,q,i,kin) = dipsum(4,q,i,kin) + alphas*dipq(4)
                enddo; enddo
              elseif(emit1.eq.1.and.emit2.eq.4) then
                do q=3,4; do i=0,2,2
                  dipsum(3,q,i,kin) = dipsum(3,q,i,kin) + alphas*dipq(3)
                  dipsum(4,q,i,kin) = dipsum(4,q,i,kin) + alphas*dipq(4)
                enddo; enddo
              elseif(emit1.eq.2.and.emit2.eq.3) then
                do q=3,4; do i=0,2,2
                  dipsum(q,3,i,kin) = dipsum(q,3,i,kin) + alphas*dipq(3)
                  dipsum(q,4,i,kin) = dipsum(q,4,i,kin) + alphas*dipq(4)
                enddo; enddo
              endif
            endif
          endif
        enddo;enddo; enddo

c qg dipoles         
        do emit1=1,4
        do spect=1,4
        if(emit1.ne.spect) then
          kin = qcdVVjjKin(id,emit1,5,spect)
          if(kin.gt.0) then
          if(lokt(kin)) then
            do i=1,4
              idB(i)=id(i)
            enddo
            if(id(5).le.2) then
              idB(emit1) = id(5)
              do i=1,4
                if(idb(i).gt.id(emit1)) idB(i) = idb(i) -1
              enddo
            endif
            call CSgeneral(p,kin,xyuz,"qg",id(emit1),id(5),id(spect),0,0,
     &            CScolAmpQZZ,idB,fsign(1),emit1,spect,nbos,12,dip(3,3,0))
            dipSum(:,:,:,kin)=dipSum(:,:,:,kin)+dip*alphas
          endif;endif
        endif
        enddo;enddo       
      endif
#endif
      
      iflav(id(5)) = 0 ! gluon

      do j=1,nfl ! quark 1
      do k=1,nfl ! quark 2
        IdFl = j.eq.k
        q1=3+mod(j,2)
        q2=3+mod(k,2)
c fill flavor information for pdf and lha output
        iflav(id(1)) = j*fsign(1)
        iflav(id(3)) = j*fsign(3)
        iflav(id(2)) = k*fsign(2)
        iflav(id(4)) = k*fsign(4)
        
        if(lnlo) then
          i = FL5qcdWZ(iflav,id,1)
          n = merge(0,1,idFl)
          if(lokt(1)) res(1) = res(1)  + POLCOL(i) * amp(0,q1,q2,n)*pdf(sign1*iflav(1),1,1)*pdf(sign2*iflav(2),2,1)
          do kin=2,16; if(lokt(kin)) then
            res(kin) = res(kin) - POLCOL(i) * dipSum(q1,q2,n,kin)*pdf(sign1*iflav(1),1,kin)*pdf(sign2*iflav(2),2,kin)
          endif; enddo
          if(iflav(1)*iflav(2).lt.0.and..not.idFl) then  ! include colstr 2 for q qb initial state
            tempi        = sign(iflav(id(1)),fsign(2))
            iflav(id(1)) = sign(iflav(id(2)),fsign(1))
            iflav(id(2)) = tempi
            i = FL5qcdWZ(iflav,id,2)
            if(lokt(1)) res(1) = res(1)  + POLCOL(i) * amp(0,q1,q2,2)*pdf(sign1*iflav(1),1,1)*pdf(sign2*iflav(2),2,1)
            do kin=2,16; if(lokt(kin)) then
              res(kin) = res(kin) - POLCOL(i) * dipSum(q1,q2,2,kin)*pdf(sign1*iflav(1),1,kin)*pdf(sign2*iflav(2),2,kin)
            endif; enddo
          endif
        else
          i = FL5qcdWZ(iflav,id,1)
          n = merge(0,1,idFl)
          res(i) = POLCOL(i) * amp(0,q1,q2,n)*pdf(sign1*iflav(1),1,1)*pdf(sign2*iflav(2),2,1)
!           print*,iflav,'a',i!,res(i)
          if(iflav(1)*iflav(2).lt.0.and..not.idFl) then  ! include colstr 2 for q qb initial state
            
            tempi        = sign(iflav(id(1)),fsign(2))
            iflav(id(1)) = sign(iflav(id(2)),fsign(1))
            iflav(id(2)) = tempi
            i = FL5qcdWZ(iflav,id,2)
            res(i) = POLCOL(i) * amp(0,q1,q2,2)*pdf(sign1*iflav(1),1,1)*pdf(sign2*iflav(2),2,1)
!             print*,iflav,'b',q1,q2,i!,res(i)
          endif
        endif
          
      enddo
      enddo
      
      end


      SUBROUTINE qg2qggZZ(P,id,fsign,lephel, res) ! sum over glu hel if ghel ==0
      implicit none 
#include "cacheZZ.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
      real*8 P(0:3,max_p,max_kin), res(0:6,3:4)
      integer kin, id(5),fsign(5),lephel 

      real*8 pbar(0:3,2),kbar(0:4,3),pdia(0:3,2),kdia(0:4,7) ! kdia for glu: 1 2 3 12 13 23 123
! !       real*8 momPsiW(0:4,2),momPsiV(0:4,2),momPsiWV(0:4,2),momPsiG(0:4,2,7)
      real*8 triD(0:4,3,2,2,2),triV(0:4,3,2,2,2) ! 3-gluon polarications (see Nucl Phys B313 560, (5.13-14) indices: triD(mu,i) with j<k, triV(j)
!       complex*16 PsiVVG(2,-1:1,3:4,2,3,2)
      complex*16 amp(6)
      complex*16 s1r,temp
      real*8 eps3G(0:3,4:6,0:2,0:2,0:2),epsTemp(0:3)
      real*8 mass2, colsumsq,dotrr
      integer denom(6),CF(6,6),idg(3),qt
      
      integer i,j,mu,qh, lepkin,ggsign(3:5),pol1,pol2,pol3,pol(0:3)
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
      

      lepkin = 1
      kin=1
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

      idg(1) = id(2)
      idg(2) = id(4)
      idg(3) = id(5)
!       print*,psiG(:,1,-id(3),id(2),1,1)*GG
      
      do pol1=1,2; do pol2=1,2
      do mu = 0,3
        eps3G(mu,4,pol1,pol2,0)= (2d0*dotrr(kdia(0,2),epsG(0,id(2),pol1,1))*epsG(mu,id(4),pol2,1)
     &              -2d0*dotrr(kdia(0,1),epsG(0,id(4),pol2,1))*epsG(mu,id(2),pol1,1)
     &              +dotrr(epsG(0,id(2),pol1,1),epsG(0,id(4),pol2,1))*(kdia(mu,1)-kdia(mu,2)))/kdia(4,4)
        eps3G(mu,5,pol1,0,pol2)= (2d0*dotrr(kdia(0,3),epsG(0,id(2),pol1,1))*epsG(mu,id(5),pol2,1)
     &              -2d0*dotrr(kdia(0,1),epsG(0,id(5),pol2,1))*epsG(mu,id(2),pol1,1)
     &              +dotrr(epsG(0,id(2),pol1,1),epsG(0,id(5),pol2,1))*(kdia(mu,1)-kdia(mu,3)))/kdia(4,5)
        eps3G(mu,6,0,pol1,pol2)= (2d0*dotrr(kdia(0,3),epsG(0,id(4),pol1,1))*epsG(mu,id(5),pol2,1)
     &              -2d0*dotrr(kdia(0,2),epsG(0,id(5),pol2,1))*epsG(mu,id(4),pol1,1)
     &              +dotrr(epsG(0,id(4),pol1,1),epsG(0,id(5),pol2,1))*(kdia(mu,2)-kdia(mu,3)))/kdia(4,6)
      enddo
      enddo;enddo

      do pol1=1,2; do pol2=1,2; do pol3=1,2
      do mu=0,3
        triD(mu,1,pol1,pol2,pol3) = (2d0*(dotrr(epsG(0,id(2),pol1,1),kdia(0,6))*eps3G(mu,6,0,pol2,pol3)
     &                    -dotrr(kdia(0,1),eps3G(0,6,0,pol2,pol3))*epsG(mu,id(2),pol1,1))
     &                 +   dotrr(epsG(0,id(2),pol1,1),eps3G(0,6,0,pol2,pol3))*(kdia(mu,1)-kdia(mu,6)))/kdia(4,7)
        triD(mu,2,pol1,pol2,pol3) = (2d0*(dotrr(epsG(0,id(4),pol2,1),kdia(0,5))*eps3G(mu,5,pol1,0,pol3)
     &                    -dotrr(kdia(0,2),eps3G(0,5,pol1,0,pol3))*epsG(mu,id(4),pol2,1))
     &                 +   dotrr(epsG(0,id(4),pol2,1),eps3G(0,5,pol1,0,pol3))*(kdia(mu,2)-kdia(mu,5)))/kdia(4,7)
        triD(mu,3,pol1,pol2,pol3) = (2d0*(dotrr(epsG(0,id(5),pol3,1),kdia(0,4))*eps3G(mu,4,pol1,pol2,0)
     &                    -dotrr(kdia(0,3),eps3G(0,4,pol1,pol2,0))*epsG(mu,id(5),pol3,1))
     &                 +   dotrr(epsG(0,id(5),pol3,1),eps3G(0,4,pol1,pol2,0))*(kdia(mu,3)-kdia(mu,4)))/kdia(4,7)
        triV(mu,1,pol1,pol2,pol3) = (2d0*epsG(mu,id(2),pol1,1)*dotrr(epsG(0,id(4),pol2,1),epsG(0,id(5),pol3,1))
     &                  -  epsG(mu,id(4),pol2,1)*dotrr(epsG(0,id(2),pol1,1),epsG(0,id(5),pol3,1))
     &                  -  epsG(mu,id(5),pol3,1)*dotrr(epsG(0,id(2),pol1,1),epsG(0,id(4),pol2,1)))/kdia(4,7)
        triV(mu,2,pol1,pol2,pol3) = (2d0*epsG(mu,id(4),pol2,1)*dotrr(epsG(0,id(2),pol1,1),epsG(0,id(5),pol3,1))
     &                  -  epsG(mu,id(2),pol1,1)*dotrr(epsG(0,id(4),pol2,1),epsG(0,id(5),pol3,1))
     &                  -  epsG(mu,id(5),pol3,1)*dotrr(epsG(0,id(2),pol1,1),epsG(0,id(4),pol2,1)))/kdia(4,7)
        triV(mu,3,pol1,pol2,pol3) = -triV(mu,1,pol1,pol2,pol3)-triV(mu,2,pol1,pol2,pol3)
      enddo
      enddo;enddo;enddo


c get amplitudes
      do qt=3,4
      do i=0,6
        res(i,qt) = 0d0
      enddo
      colsumsq = 0d0

c get amps ans ampsq for all polarisations
      
      do qh=-1,1,2
      do pol1=1,2; do pol2=1,2; do pol3=1,2

      do perm = 1,6
        ig(1) = permutations(1,perm)
        ig(2) = permutations(2,perm)
        ig(3) = permutations(3,perm)
        ! 12 13 and 23 ordered
        if(perm.le.3) then; ggsign(3) = 1; else; ggsign(3) = -1; endif
        if(perm.le.2.or.perm.eq.4) then; ggsign(4) = 1; else; ggsign(4) = -1; endif
        if(perm/2.eq.2.or.perm.eq.1) then; ggsign(5) = 1; else; ggsign(5) = -1; endif
        amp(perm)=0d0
c 1 gluon splitting amplitudes
        i = ig(2)+ig(3)+1
        pol(1) = pol1; pol(2) = pol2; pol(3) = pol3
        pol(0) = pol(ig(1)); pol(ig(1)) = 0 ! swap to get convention of epsG and psi*G
        temp = s1r(psiVG(1,qh,qt,-id(3),3,idg(ig(1)),pol(0)),eps3G(0,i,pol(1),pol(2),pol(3)),.true.,qh,Psi(1,qh,id(1),1))
     &        +s1r(psiG(1,qh,-id(3),idg(ig(1)),pol(0),1),eps3G(0,i,pol(1),pol(2),pol(3)),.true.,qh,PsiV(1,qh,qt,id(1),3,1))
        if(lAZ(0,1)) then
          temp=temp
     &        +s1r(psiVG(1,qh,qt,-id(3),1,idg(ig(1)),pol(0)),eps3G(0,i,pol(1),pol(2),pol(3)),.true.,qh,PsiV(1,qh,qt,id(1),2,1))
     &        +s1r(psiVG(1,qh,qt,-id(3),2,idg(ig(1)),pol(0)),eps3G(0,i,pol(1),pol(2),pol(3)),.true.,qh,PsiV(1,qh,qt,id(1),1,1))
        endif
        amp(perm) = -temp*ggsign(i-1)

        i = ig(1)+ig(2)+1
        pol(ig(1)) = pol(0); pol(0) = pol(ig(3)); pol(ig(3)) = 0 ! swap back and swap for next calls
        temp = s1r(psiV(1,qh,qt,-id(3),3,1),eps3G(0,i,pol(1),pol(2),pol(3)),.true.,qh,PsiG(1,qh,id(1),idg(ig(3)),pol(0),1))
     &        +s1r(psi(1,qh,-id(3),1),eps3G(0,i,pol(1),pol(2),pol(3)),.true.,qh,PsiVG(1,qh,qt,id(1),3,idg(ig(3)),pol(0)))
        if(lAZ(0,1)) then
          temp=temp
     &        +s1r(psiV(1,qh,qt,-id(3),1,1),eps3G(0,i,pol(1),pol(2),pol(3)),.true.,qh,psiVG(1,qh,qt,id(1),2,idg(ig(3)),pol(0)))
     &        +s1r(psiV(1,qh,qt,-id(3),2,1),eps3G(0,i,pol(1),pol(2),pol(3)),.true.,qh,psiVG(1,qh,qt,id(1),1,idg(ig(3)),pol(0)))
        endif
        amp(perm) = amp(perm)-temp*ggsign(i-1)
        pol(ig(3)) = pol(0) ! swap back
        
c g(g1,g2,g3) amplitudes        
        do mu=0,3
          epsTemp(mu) = ggsign(ig(2)+ig(3))*triD(mu,ig(1),pol1,pol2,pol3)-ggsign(ig(1)+ig(2))*triD(mu,ig(3),pol1,pol2,pol3)+
     $                         triV(mu,ig(2),pol1,pol2,pol3)
        enddo
        amp(perm)= amp(perm)
     &           +s1r(psi(1,qh,-id(3),1),epsTemp,.true.,qh,PsiV(1,qh,qt,id(1),3,1))
     &           +s1r(psiV(1,qh,qt,-id(3),3,1),epsTemp,.true.,qh,Psi(1,qh,id(1),1))
        if(lAZ(0,1)) then
          amp(perm)= amp(perm)
     &           +s1r(psiV(1,qh,qt,-id(3),1,1),epsTemp,.true.,qh,PsiV(1,qh,qt,id(1),2,1))
     &           +s1r(psiV(1,qh,qt,-id(3),2,1),epsTemp,.true.,qh,PsiV(1,qh,qt,id(1),1,1))
        endif

c amplitudes without glu splitting
        amp(perm) = amp(perm)
     &             +s1r(psiVG(1,qh,qt,-id(3),3,idg(ig(1)),pol(ig(1))),epsG(0,idg(ig(2)),pol(ig(2)),1),.false.,qh,
     &                  PsiG(1,qh,id(1),idg(ig(3)),pol(ig(3)),1))
     &             +s1r(psiG(1,qh,-id(3),idg(ig(1)),pol(ig(1)),1),epsG(0,idg(ig(2)),pol(ig(2)),1),.false.,qh,
     &                  PsiVG(1,qh,qt,id(1),3,idg(ig(3)),pol(ig(3))))
        if(lAZ(0,1)) then
          amp(perm) = amp(perm)
     &             +s1r(psiVG(1,qh,qt,-id(3),1,idg(ig(1)),pol(ig(1))),epsG(0,idg(ig(2)),pol(ig(2)),1),.false.,qh,
     &                  psiVG(1,qh,qt,id(1),2,idg(ig(3)),pol(ig(3))))
     &             +s1r(psiVG(1,qh,qt,-id(3),2,idg(ig(1)),pol(ig(1))),epsG(0,idg(ig(2)),pol(ig(2)),1),.false.,qh,
     &                  psiVG(1,qh,qt,id(1),1,idg(ig(3)),pol(ig(3))))
        endif
        amp(perm) = amp(perm) * G**3
!         print*,amp(perm)
      enddo ! perm
!       print*,amp

c square amp
      DO I = 1, 6
          TEMP = (0.D0,0.D0)
          DO J = 1, 6
              TEMP = TEMP + CF(J,I)*AMP(J)
          ENDDO
          res(0,qt) =res(0,qt)+TEMP*DCONJG(AMP(I))/DENOM(I)
          res(i,qt) =res(i,qt)+ AMP(I)*DCONJG(AMP(I))/DENOM(I)
      ENDDO
      enddo; enddo;enddo ! pol
      enddo ! qhel
      
      do i=1,6
        colsumsq = colsumsq+res(i,qt)
      enddo
      DO i=1,6
         res(i,qt) = res(i,qt)/colsumsq*res(0,qt)
      ENDDO
      enddo !qt

      end


      subroutine amp4q1gzz(P,id,fsign,lephel,res)
      implicit none
#include "common.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "cacheZZ.inc"
      real*8 P(0:3,max_p,max_kin),res(0:4,3:4,3:4,0:2)
      integer id(5),fsign(5),lephel

      
      real*8 pbar(0:3,5),pdia(0:3,5)
      complex*16 currij(0:3,-1:1,1:2,3:4),currGij(0:3,-1:1,1:2,3:4,2,2)
      complex*16 currVij(0:3,-1:1,3:4,1:2,3:4,3),currVGij1(0:3,-1:1,3:4,1:2,3:4,3,2),currVGij2(0:3,-1:1,3:4,1:2,3:4,3,2)
      complex*16 tempCurr(0:3,-1:1,4),amp_ort(0:2),amp(4,2),temp
      complex*16 s1c, dotrc
      complex*16 dotcc
      external s1c, dotcc, dotrc,dotrr,qjj2,mass2, qjj2Q
      real*8 mass2,qjj2,qjj2Q,qVij(0:4,1:2,3:4,0:3)
      real*8 tempR
      real*8 dotrr

      integer helswitch
      common/helswitchTemp/ helswitch

      integer i,j,k,l,m,mu,colstr,qt,iV,pol,q1,q2,h1,h2
      integer in,out,in2,out2
      integer idi,ido

      integer color_Z(2,2,4),color_C(2,2)
      data (color_Z(2,2,i),i=1,4) / 2, 2, 2, 2 /
      data (color_Z(1,2,i),i=1,4) / 2,-2,-4, 4 /
      data (color_Z(2,1,i),i=1,4) /-4, 4, 2,-2 /
      data (color_Z(1,1,i),i=1,4) /-4,-4,-4,-4 /
      data color_C(2,2), color_C(1,2), color_C(2,1), color_C(1,1) / 5, 3, 3, 1 /

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
      

c numbering convention
c in   q1 >>> #### >>> q3  out
c             ####          , WZ
c             ####
c in2  q2 >>> #### >>> q4  out2
!       do i=1,5
!       print*,(p(mu,i,1),mu=0,3)
!       enddo
      do mu=0,3
        pbar(mu,1)=p(mu,id(1),1)
        pbar(mu,2)=p(mu,id(2),1)
        pbar(mu,3)=p(mu,id(3),1)
        pbar(mu,4)=p(mu,id(4),1)
        pbar(mu,5)=p(mu,id(5),1)
        pdia(mu,1)=pbar(mu,1)*fsign(1)
        pdia(mu,2)=pbar(mu,2)*fsign(2)
        pdia(mu,3)=pbar(mu,3)*fsign(3)
        pdia(mu,4)=pbar(mu,4)*fsign(4)
        pdia(mu,5)=pbar(mu,5)*fsign(5)
      enddo      
      
      
      do in = 1,2  ! determine q line coupling to W
      do out= 3,4
        idi=id(in)
        ido=-id(out)
c qbar q
        call curr(1,psi(1,-1,ido,1),psi(1,-1,idi,1),currij(0,-1,in,out))
c qbar V q
        do iV=1,3
        if(lAZ(0,iV)) then
        do qt=3,4
          call curr(1,psiV(1,-1,qt,ido,iV,1),psi(1,-1,idi,1),tempCurr(0,-1,1))
          call curr(1,psi(1,-1,ido,1),psiV(1,-1,qt,idi,iV,1),tempCurr(0,-1,2))
          currVij(:,:,qt,in,out,iV) = tempCurr(:,:,1)+tempCurr(:,:,2)
          if(iV.eq.3.and.lAZ(0,1)) then
            call curr(1,psiV(1,-1,qt,ido,1,1),psiV(1,-1,qt,idi,2,1),tempCurr(0,-1,1))
            call curr(1,psiV(1,-1,qt,ido,2,1),psiV(1,-1,qt,idi,1,1),tempCurr(0,-1,2))
            currVij(:,:,qt,in,out,iV) = currVij(:,:,qt,in,out,iV) + tempCurr(:,:,1)+tempCurr(:,:,2)
          endif
        enddo;endif;enddo
c qbar g q
        do pol=1,2
c g on left side
          call curr(1,psiG(1,-1,ido,id(5),pol,1),psi(1,-1,idi,1),currGij(0,-1,in,out,pol,1))
c g on right side
          call curr(1,psi(1,-1,ido,1),psiG(1,-1,idi,id(5),pol,1),currGij(0,-1,in,out,pol,2))
c qqbar Vg q
          do iV=1,3
          if(lAZ(0,iV)) then
          do qt=3,4
c g on left side
            call curr(1,psiG(1,-1,ido,id(5),pol,1),psiV(1,-1,qt,idi,iV,1),tempCurr(0,-1,1))
            call curr(1,psiVG(1,-1,qt,ido,iV,id(5),pol),psi(1,-1,idi,1),tempCurr(0,-1,2))
            currVGij1(:,:,qt,in,out,iV,pol)=tempCurr(:,:,1)+tempCurr(:,:,2)
            if(iV.eq.3.and.lAZ(0,1)) then
              call curr(1,psiVG(1,-1,qt,ido,1,id(5),pol),psiV(1,-1,qt,idi,2,1),tempCurr(0,-1,1))
              call curr(1,psiVG(1,-1,qt,ido,2,id(5),pol),psiV(1,-1,qt,idi,1,1),tempCurr(0,-1,2))
              currVGij1(:,:,qt,in,out,iV,pol)=currVGij1(:,:,qt,in,out,iV,pol)+tempCurr(:,:,1)+tempCurr(:,:,2)
            endif
c g on right side
            call curr(1,psi(1,-1,ido,1),psiVG(1,-1,qt,idi,iV,id(5),pol),tempCurr(0,-1,1))
            call curr(1,psiV(1,-1,qt,ido,iV,1),psiG(1,-1,idi,id(5),pol,1),tempCurr(0,-1,2))
            currVGij2(:,:,qt,in,out,iV,pol)=tempCurr(:,:,1)+tempCurr(:,:,2)
            if(iV.eq.3.and.lAZ(0,1)) then
              call curr(1,psiV(1,-1,qt,ido,1,1),psiVG(1,-1,qt,idi,2,id(5),pol),tempCurr(0,-1,1))
              call curr(1,psiV(1,-1,qt,ido,2,1),psiVG(1,-1,qt,idi,1,id(5),pol),tempCurr(0,-1,2))
              currVGij2(:,:,qt,in,out,iV,pol)=currVGij2(:,:,qt,in,out,iV,pol)+tempCurr(:,:,1)+tempCurr(:,:,2)
            endif
          enddo;endif;enddo
        enddo ! pol
        do mu=0,3
          qVij(mu,in,out,0) = pdia(mu,in)-pdia(mu,out)
          qVij(mu,in,out,1) = momPsiV(mu,id(in),1,1)-pdia(mu,out)
          qVij(mu,in,out,2) = momPsiV(mu,id(in),2,1)-pdia(mu,out)
          qVij(mu,in,out,3) = momPsiV(mu,id(in),3,1)-pdia(mu,out)
        enddo
        qVij(4,in,out,0) = mass2(qVij(0,in,out,0))
        qVij(4,in,out,1) = mass2(qVij(0,in,out,1))
        qVij(4,in,out,2) = mass2(qVij(0,in,out,2))
        qVij(4,in,out,3) = mass2(qVij(0,in,out,3))

      enddo; enddo ! in, out
!       print*,psiG(:,-1,id(1),id(5),2,1)*G
!       print*,currGij(:,-1,1,3,2,2)*G**2/qsq(4,2,4,3)
!       print*
      
c calculate amplitude
      do q1=3,4 ! qt of 2 quark lines
      do q2=3,4
        do i=0,4; do j=0,2
          res(i,q1,q2,j) = 0d0
        enddo; enddo
      do h1=-1,1,2
      do h2=-1,1,2
      do pol=1,2

      amp = 0d0
      out=3 ! fix final state
      out2=4
      do in=1,2
        in2=3-in
        colstr = in

        if(lAZ(0,1)) then
c in>V1>out, in2>V2 >out2
         amp(out2,colstr)=amp(out2,colstr) + dotcc(currVij(0,h1,q1,in,out,1),currVGij1(0,h2,q2,in2,out2,2,pol))/qVij(4,in,out,1)
         amp( in2,colstr)=amp( in2,colstr) + dotcc(currVij(0,h1,q1,in,out,1),currVGij2(0,h2,q2,in2,out2,2,pol))/qVij(4,in,out,1)
         amp(out ,colstr)=amp(out ,colstr) + dotcc(currVGij1(0,h1,q1,in,out,1,pol),currVij(0,h2,q2,in2,out2,2))/qVij(4,in2,out2,2)
         amp( in ,colstr)=amp( in ,colstr) + dotcc(currVGij2(0,h1,q1,in,out,1,pol),currVij(0,h2,q2,in2,out2,2))/qVij(4,in2,out2,2)
!       3 glu vert
        temp = 1d0/qVij(4,in,out,1)/qVij(4,in2,out2,2)
     &          *( dotrc(pdia(0,5),currVij(0,h1,q1,in,out,1))*dotrc(epsG(0,id(5),pol,1),currVij(0,h2,q2,in2,out2,2))
     &            +dotrr(qVij(0,in,out,1),epsG(0,id(5),pol,1))*dotcc(currVij(0,h1,q1,in,out,1),currVij(0,h2,q2,in2,out2,2))
     &            -dotrc(pdia(0,5),currVij(0,h2,q2,in2,out2,2))*dotrc(epsG(0,id(5),pol,1),currVij(0,h1,q1,in,out,1)))
        amp(out,colstr) = amp(out,colstr) - 2d0*temp
        amp(in ,colstr) = amp(in, colstr) + 2d0*temp

! c in>V2>out, in2>V1>out2
         amp(out2,colstr)=amp(out2,colstr) + dotcc(currVij(0,h1,q1,in,out,2),currVGij1(0,h2,q2,in2,out2,1,pol))/qVij(4,in,out,2)
         amp( in2,colstr)=amp( in2,colstr) + dotcc(currVij(0,h1,q1,in,out,2),currVGij2(0,h2,q2,in2,out2,1,pol))/qVij(4,in,out,2)
         amp(out ,colstr)=amp(out ,colstr) + dotcc(currVGij1(0,h1,q1,in,out,2,pol),currVij(0,h2,q2,in2,out2,1))/qVij(4,in2,out2,1)
         amp( in ,colstr)=amp( in ,colstr) + dotcc(currVGij2(0,h1,q1,in,out,2,pol),currVij(0,h2,q2,in2,out2,1))/qVij(4,in2,out2,1)
!       3 glu vert
        temp = 1d0/qVij(4,in,out,2)/qVij(4,in2,out2,1)
     &          *( dotrc(pdia(0,5),currVij(0,h1,q1,in,out,2))*dotrc(epsG(0,id(5),pol,1),currVij(0,h2,q2,in2,out2,1))
     &            +dotrr(qVij(0,in,out,2),epsG(0,id(5),pol,1))*dotcc(currVij(0,h1,q1,in,out,2),currVij(0,h2,q2,in2,out2,1))
     &            -dotrc(pdia(0,5),currVij(0,h2,q2,in2,out2,1))*dotrc(epsG(0,id(5),pol,1),currVij(0,h1,q1,in,out,2)))
        amp(out,colstr) = amp(out,colstr) - 2d0*temp
        amp(in ,colstr) = amp(in, colstr) + 2d0*temp

        endif
!         if(lAZ(0,3)) then ! should always be true
c in>V12>out, in2>out2
        amp(out2,colstr) = amp(out2,colstr) + dotcc(currVij(0,h1,q1,in,out,3),currGij(0,h2,in2,out2,pol,1))/qVij(4,in,out,3)
        amp( in2,colstr) = amp( in2,colstr) + dotcc(currVij(0,h1,q1,in,out,3),currGij(0,h2,in2,out2,pol,2))/qVij(4,in,out,3)
        amp(out ,colstr) = amp(out ,colstr) + dotcc(currVGij1(0,h1,q1,in,out,3,pol),currij(0,h2,in2,out2))/qVij(4,in2,out2,0)
        amp( in ,colstr) = amp( in ,colstr) + dotcc(currVGij2(0,h1,q1,in,out,3,pol),currij(0,h2,in2,out2))/qVij(4,in2,out2,0)
!       3 glu vert
        temp = 1d0/qVij(4,in,out,3)/qVij(4,in2,out2,0)
     &          *( dotrc(pdia(0,5),currVij(0,h1,q1,in,out,3))*dotrc(epsG(0,id(5),pol,1),currij(0,h2,in2,out2))
     &            +dotrr(qVij(0,in,out,3),epsG(0,id(5),pol,1))*dotcc(currVij(0,h1,q1,in,out,3),currij(0,h2,in2,out2))
     &            -dotrc(pdia(0,5),currij(0,h2,in2,out2))*dotrc(epsG(0,id(5),pol,1),currVij(0,h1,q1,in,out,3)))
        amp(out,colstr) = amp(out,colstr) - 2d0*temp
        amp(in ,colstr) = amp(in, colstr) + 2d0*temp

c in>out, in2>V12>out2
        amp(out2,colstr) = amp(out2,colstr) + dotcc(currij(0,h1,in,out),currVGij1(0,h2,q2,in2,out2,3,pol))/qVij(4,in,out,0)
        amp( in2,colstr) = amp( in2,colstr) + dotcc(currij(0,h1,in,out),currVGij2(0,h2,q2,in2,out2,3,pol))/qVij(4,in,out,0)
        amp(out ,colstr) = amp(out ,colstr) + dotcc(currGij(0,h1,in,out,pol,1),currVij(0,h2,q2,in2,out2,3))/qVij(4,in2,out2,3)
        amp( in ,colstr) = amp( in ,colstr) + dotcc(currGij(0,h1,in,out,pol,2),currVij(0,h2,q2,in2,out2,3))/qVij(4,in2,out2,3)
!       3 glu vert
        temp = 1d0/qVij(4,in,out,0)/qVij(4,in2,out2,3)
     &          *( dotrc(pdia(0,5),currij(0,h1,in,out))*dotrc(epsG(0,id(5),pol,1),currVij(0,h2,q2,in2,out2,3))
     &            +dotrr(qVij(0,in,out,0),epsG(0,id(5),pol,1))*dotcc(currij(0,h1,in,out),currVij(0,h2,q2,in2,out2,3))
     &            -dotrc(pdia(0,5),currVij(0,h2,q2,in2,out2,3))*dotrc(epsG(0,id(5),pol,1),currij(0,h1,in,out)))
        amp(out,colstr) = amp(out,colstr) - 2d0*temp
        amp(in ,colstr) = amp(in, colstr) + 2d0*temp
      enddo

c square amplitudes
      amp = amp*G**3
      do k=1,2; do l=1,2
        amp_ort = 0d0
        do m=1,4
          amp_ort(1) = amp_ort(1) + amp(m,1)*color_Z(k,l,m)
        enddo
        temp = (k*2)-3
        amp_ort(2) = amp_ort(2) - temp*amp(1,2)*color_Z(k,l,2) ! change of initial state changes nnumbering of amplitudes!
        amp_ort(2) = amp_ort(2) - temp*amp(2,2)*color_Z(k,l,1)
        amp_ort(2) = amp_ort(2) - temp*amp(3,2)*color_Z(k,l,3)
        amp_ort(2) = amp_ort(2) - temp*amp(4,2)*color_Z(k,l,4)
        
        amp_ort(0) = amp_ort(1) + amp_ort(2)
        
c different quarks ( if q1==q2 -> different families)
        tempR=dble(amp_ort(1)*dconjg(amp_ort(1)))*color_c(k,l)/36d0
        res(0,q1,q2,1) = res(0,q1,q2,1) + tempR
        tempR=dble(amp_ort(2)*dconjg(amp_ort(2)))*color_c(k,l)/36d0
        res(0,q1,q2,2) = res(0,q1,q2,2) + tempR
c same quarks, interferences for equal helicities
        if(q1.eq.q2) then
          if(h1.eq.h2) then 
            tempR=dble(amp_ort(0)*dconjg(amp_ort(0)))*color_c(k,l)/36d0
            res(0,q1,q2,0) = res(0,q1,q2,0) + tempR
          else
            tempR=dble(amp_ort(1)*dconjg(amp_ort(1)))*color_c(k,l)/36d0
            res(0,q1,q2,0) = res(0,q1,q2,0) + tempR
            tempR=dble(amp_ort(2)*dconjg(amp_ort(2)))*color_c(k,l)/36d0
            res(0,q1,q2,0) = res(0,q1,q2,0) + tempR
          endif
        endif
      enddo; enddo ! k,l
            
      enddo ! pol
      enddo; enddo ! hel
      enddo; enddo ! qt
#ifdef NINH_COMPARE
cc      restemp = restemp + res
cc      enddo
cc      res = restemp
#endif
      end

       
