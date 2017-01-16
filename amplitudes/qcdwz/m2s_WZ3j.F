c ************************************************************************************
c Author: Matthias Kerner
c Date: 2012/13
c ************************************************************************************

! #define NINH_COMPARE

      subroutine m2s_qcdwz3j(
     &                   xi,  !in:  Feynman x parameters of incoming fermions
     &                   p,   !in:  momenta associated with external fermions
     &                   v,     !in:  momenta associated with Z/H dec. fermions
     &                   nbos, ! 32 = WPZ, 42 = WMZ
     &                   lnlo,
     &                   lokt,
     &                   xyuz,
     &                   rn,
     &                   m2s
     &                        )

      use globalvars, only: ldoblha

      implicit none
c 
c declare input/output variables
c
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"
#include "common.inc"
      integer nbos
      real*8 xi(nx), p(0:3,max_p,max_kin), v(0:3,max_v,max_kin), m2s(0:max_kin)
      real*8 xyuz(2,2:max_kin),rn(2)
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
      real*8 pdf(-6:6,2,22)
c local varriables
      integer FL5qcdWZ, ires_min, ires_max
      save ires_min, ires_max
      data ires_min,ires_max / 1, maxnumsubproc/
      external FL5qcdWZ
      logical FLskip(2)
      logical ldebug,lallsub
!       parameter (ldebug=.false.)
      real*8 mjj2
      integer qcdVVjjkin

      double precision  clr, xm2, xmg, b
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6)
      
      integer helswitch
      common/helswitchTemp/ helswitch


      ldebug = .false.

      if ( (init .eq. 0) .and. (.not. ldoblha) ) then
        if(nbos .eq. 32) then
         write(6,*) "QCD W+Z jjj amplitude square information:"
        elseif(nbos.eq.42) then
         write(6,*) "QCD W-Z jjj amplitude square information:"
        elseif(nbos.eq.31) then
         write(6,*) "QCD W+A jjj amplitude square information:"
        elseif(nbos.eq.41) then
         write(6,*) "QCD W-A jjj amplitude square information:"
	elseif(nbos.eq.3) then
	 write(6,*) "QCD W+ jjj amplitude square information:"
	else
	 write(6,*) "QCD W- jjj amplitude square information:"
        endif
         if(nfl.eq.5) then
          write(6,*) "external b quarks are included"
         else
          write(6,*) "external b quarks are excluded"
         endif
         write(6,*) "alpha parameter of dipoles: ",alphaDip
         write(6,*) "-----------------------------------------------"
!         call CSGeneralSetAlpha(alphaDip)
      else if ( (init .eq. 0) .and. ldoblha ) then
        init = 1
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

      lallsub=.true.

      do i=ires_min,ires_max
	res(i) = 0d0
      enddo
      
      
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
      
      endif !ldoblha 


      call calcEWpartsNew(v,lokt,nbos,rn(1))

      ires_min = FL5qcdWZ(id,id, -100000) +1


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


      call calcEWpartsNew(v,lokt,nbos,rn(1))
	rn(1) = -1d0
	

      pMG(:,1) = p(:,1,1)
      pMG(:,2) = p(:,2,1)
      pMG(:,3) = p(:,4,1)
      pMG(:,4) = p(:,3,1)
      pMG(:,5) = p(:,5,1)
      pMG(:,6) = v(:,1,1)
      pMG(:,7) = v(:,2,1)
      pMG(:,8) = v(:,3,1)
      pMG(:,9) = v(:,4,1)
!       call SMATRIXMD(pMG,m2s(0))
!       call DIPOLSUM(pMG,xi,1d0,m2s(1))
!       print*,m2s(0:1)
!       pause
#endif
      
#ifdef NINH_COMPARE
      if(nbos.gt.10) then
        call precalcQCDWZjj(p,lokt,1,nbos)
      endif
      helswitch=1
      goto 24
   23 continue
      call precalcQCDWZjj(p,lokt,-1,nbos)
      helswitch=-1
   24 continue
#else
      call precalcQCDWZjj(p,lokt,0,nbos)
#endif

c reset cache of born amplitudes
      call CSColAmpGNewReset()
      call CSColAmpQNewReset()

       if (ldoblha) then
c id is diagtoPhys in standard VBFNLO terminology
         id(blha_physdiagmap(1,blha_cursubproc,blha_curproc))=1
         id(blha_physdiagmap(2,blha_cursubproc,blha_curproc))=2
         id(blha_physdiagmap(3,blha_cursubproc,blha_curproc))=3
         id(blha_physdiagmap(4,blha_cursubproc,blha_curproc))=4
         id(blha_physdiagmap(5,blha_cursubproc,blha_curproc))=5
         fsign(1) = blha_fsign(1,blha_cursubproc,blha_curproc)
         fsign(2) = blha_fsign(2,blha_cursubproc,blha_curproc)
         fsign(3) = blha_fsign(3,blha_cursubproc,blha_curproc)
         fsign(4) = blha_fsign(4,blha_cursubproc,blha_curproc)
         fsign(5) = blha_fsign(5,blha_cursubproc,blha_curproc)

         if ((blha_idsubproc(blha_cursubproc,blha_curproc)/100) .eq. 1) then
           call getRes2Q3GWZ(p,id,fsign,pdf,nbos,lnlo,lokt,xyuz,res)
         else if ((blha_idsubproc(blha_cursubproc,blha_curproc)/100) .eq. 2) then
           call getRes4Q1GWZ(p,id,fsign,FLskip,pdf,nbos,lnlo,lokt,xyuz,res)
         endif
         return
       endif

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
      
      if(sub_number.eq.1.or.lallsub) call getRes2Q3GWZ(p,id,fsign,pdf,nbos,lnlo,lokt,xyuz,res)

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
 
      if(sub_number.eq.2.or.lallsub) call getRes2Q3GWZ(p,id,fsign,pdf,nbos,lnlo,lokt,xyuz,res)

c switch initial state
      id(1) = 2
      id(2) = 1

      if(sub_number.eq.3.or.lallsub) call getRes2Q3GWZ(p,id,fsign,pdf,nbos,lnlo,lokt,xyuz,res)

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

      if(sub_number.eq.4.or.lallsub)        call getRes2Q3GWZ(p,id,fsign,pdf,nbos,lnlo,lokt,xyuz,res)

c switch initial state
      id(3) = 2
      id(2) = 1

      if(sub_number.eq.5.or.lallsub)        call getRes2Q3GWZ(p,id,fsign,pdf,nbos,lnlo,lokt,xyuz,res)

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

      if(sub_number.eq.6.or.lallsub)        call getRes2Q3GWZ(p,id,fsign,pdf,nbos,lnlo,lokt,xyuz,res)

c switch initial state
      id(1) = 2
      id(3) = 1

      if(sub_number.eq.7.or.lallsub)        call getRes2Q3GWZ(p,id,fsign,pdf,nbos,lnlo,lokt,xyuz,res)

c
c
c  end of 2q 3g processes
c  begin 4q 1g processes
c
c
c--crossing 1----------------------------------------------------------
c     q1 qb3 -> qb2 q4 g wz  diagram
c     p1  p2    p3 p4 p5     phasespace

      id(1) = 1
      id(2) = 3
      id(3) = 2
      id(4) = 4
      id(5) = 5

      fsign(1) =  1
      fsign(2) = -1
      fsign(3) = -1
      fsign(4) =  1
      fsign(5) =  1

      FLskip(1) = .false.
      FLskip(2) = .false.

      if(sub_number.eq.8.or.lallsub)        call getRes4Q1GWZ(p,id,fsign,FLskip,pdf,nbos,lnlo,lokt,xyuz,res)

      id(1) = 2
      id(3) = 1

      if(sub_number.eq.9.or.lallsub)        call getRes4Q1GWZ(p,id,fsign,FLskip,pdf,nbos,lnlo,lokt,xyuz,res)

c--crossing 2----------------------------------------------------------
c     q2 qb4 -> qb1 q3 g wz  diagram
c     p1  p2    p3 p4 p5     phasespace

      id(1) = 3
      id(2) = 1
      id(3) = 4
      id(4) = 2
      id(5) = 5

c  Change momnentum convention:
      fsign(1) = -1
      fsign(2) =  1
      fsign(3) =  1
      fsign(4) = -1
      fsign(5) =  1

      if(sub_number.eq.10.or.lallsub)        call getRes4Q1GWZ(p,id,fsign,FLskip,pdf,nbos,lnlo,lokt,xyuz,res)

      id(2) = 2
      id(4) = 1

      if(sub_number.eq.11.or.lallsub)        call getRes4Q1GWZ(p,id,fsign,FLskip,pdf,nbos,lnlo,lokt,xyuz,res)

c--crossing 3----------------------------------------------------------
c     q1 q2 -> q3 q4 g wz  diagram
c     p1 p2    p3 p4 p5     phasespace

      id(1) = 1
      id(2) = 2
      id(3) = 4
      id(4) = 3
      id(5) = 5

c  Change momnentum convention:
      fsign(1) =  1
      fsign(2) =  1
      fsign(3) =  1
      fsign(4) =  1
      fsign(5) =  1
      
      flskip(1) = .false.
      flskip(2) = .false.
      
      if(sub_number.eq.12.or.lallsub)        call getRes4Q1GWZ(p,id,fsign,FLskip,pdf,nbos,lnlo,lokt,xyuz,res)

      id(1) = 2
      id(2) = 1

      FLskip(1) = .true.
      if(sub_number.eq.13.or.lallsub)        call getRes4Q1GWZ(p,id,fsign,FLskip,pdf,nbos,lnlo,lokt,xyuz,res)

c--crossing 4----------------------------------------------------------
c     q1 qb4 -> qb2 q3 g wz  diagram
c     p1  p2    p3 p4 p5     phasespace

      id(1) = 1
      id(2) = 3
      id(3) = 4
      id(4) = 2
      id(5) = 5

c  Change momnentum convention:
      fsign(1) =  1
      fsign(2) = -1
      fsign(3) =  1
      fsign(4) = -1
      fsign(5) =  1

      FLskip(1) = .true.
      FLskip(2) = .true.

      if(sub_number.eq.14.or.lallsub)        call getRes4Q1GWZ(p,id,fsign,FLskip,pdf,nbos,lnlo,lokt,xyuz,res)

      id(1) = 2
      id(4) = 1

      if(sub_number.eq.15.or.lallsub)        call getRes4Q1GWZ(p,id,fsign,FLskip,pdf,nbos,lnlo,lokt,xyuz,res)

c--crossing 5----------------------------------------------------------
c     qb3 q2 -> qb1 q4 g wz  diagram
c     p1  p2    p3 p4 p5     phasespace

      id(1) = 3
      id(2) = 2
      id(3) = 1
      id(4) = 4
      id(5) = 5

c  Change momnentum convention:
      fsign(1) = -1
      fsign(2) =  1
      fsign(3) = -1
      fsign(4) =  1
      fsign(5) =  1

      FLskip(1) = .true.
      FLskip(2) = .true.

      if(sub_number.eq.16.or.lallsub)        call getRes4Q1GWZ(p,id,fsign,FLskip,pdf,nbos,lnlo,lokt,xyuz,res)

      id(3) = 2
      id(2) = 1

      if(sub_number.eq.17.or.lallsub)        call getRes4Q1GWZ(p,id,fsign,FLskip,pdf,nbos,lnlo,lokt,xyuz,res)

c--crossing 6----------------------------------------------------------
c     qb3 qb4 -> qb1 qb2 g wz  diagram
c     p1  p2     p3  p4 p5     phasespace

      id(1) = 3
      id(2) = 4
      id(3) = 1
      id(4) = 2
      id(5) = 5

c  Change momnentum convention:
      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = -1
      fsign(4) = -1
      fsign(5) =  1

      FLskip(1) = .false.
      FLskip(2) = .false.

      if(sub_number.eq.18.or.lallsub)        call getRes4Q1GWZ(p,id,fsign,FLskip,pdf,nbos,lnlo,lokt,xyuz,res)

      id(3) = 2
      id(4) = 1

      FLskip(2) = .true.
      if(sub_number.eq.19.or.lallsub)        call getRes4Q1GWZ(p,id,fsign,FLskip,pdf,nbos,lnlo,lokt,xyuz,res)

c--crossing 7----------------------------------------------------------
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

      FLskip(1) = .false.
      FLskip(2) = .false.

      if(sub_number.eq.20.or.lallsub)        call getRes4Q1GWZ(p,id,fsign,FLskip,pdf,nbos,lnlo,lokt,xyuz,res)

      id(5) = 2
      id(1) = 1

      if(sub_number.eq.21.or.lallsub)        call getRes4Q1GWZ(p,id,fsign,FLskip,pdf,nbos,lnlo,lokt,xyuz,res)

c--crossing 8----------------------------------------------------------
c     g5 q2 -> qb1 q3 q4 wz  diagram
c     p1  p2    p3 p4 p5     phasespace

      id(1) = 3
      id(2) = 2
      id(3) = 5
      id(4) = 4
      id(5) = 1

      fsign(1) = -1
      fsign(2) =  1
      fsign(3) =  1
      fsign(4) =  1
      fsign(5) = -1

      FLskip(1) = .true.

      if(sub_number.eq.22.or.lallsub)        call getRes4Q1GWZ(p,id,fsign,FLskip,pdf,nbos,lnlo,lokt,xyuz,res)

      id(5) = 2
      id(2) = 1

      if(sub_number.eq.23.or.lallsub)        call getRes4Q1GWZ(p,id,fsign,FLskip,pdf,nbos,lnlo,lokt,xyuz,res)

c--crossing 9----------------------------------------------------------
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

      FLskip(1) = .false.
      FLskip(2) = .false.

      if(sub_number.eq.24.or.lallsub)        call getRes4Q1GWZ(p,id,fsign,FLskip,pdf,nbos,lnlo,lokt,xyuz,res)

      id(5) = 2
      id(3) = 1

      if(sub_number.eq.25.or.lallsub)        call getRes4Q1GWZ(p,id,fsign,FLskip,pdf,nbos,lnlo,lokt,xyuz,res)

c--crossing 10----------------------------------------------------------
c     g5 qb4 -> qb2 qb1 q3 wz  diagram
c     p1  p2    p3  p4  p5     phasespace
! 
      id(1) = 4
      id(2) = 3
      id(3) = 5
      id(4) = 2
      id(5) = 1

      fsign(1) = -1
      fsign(2) = -1
      fsign(3) =  1
      fsign(4) = -1
      fsign(5) = -1

      FLskip(1) = .false.
      FLskip(2) = .true.

      if(sub_number.eq.26.or.lallsub)        call getRes4Q1GWZ(p,id,fsign,FLskip,pdf,nbos,lnlo,lokt,xyuz,res)

      id(5) = 2
      id(4) = 1

      if(sub_number.eq.27.or.lallsub)        call getRes4Q1GWZ(p,id,fsign,FLskip,pdf,nbos,lnlo,lokt,xyuz,res)

      FLskip(2) = .false.

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
	 if(ldebug) print*, "Number of (color) subprocesses: ", ires_max
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
         if(weight.ge.rn(2)) then
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

#ifdef NINH_COMPARE
      if(helswitch.eq.1) then
        goto 23
      endif
      m2s = m2s
      print*,m2s(1),m2s(0)-m2s(1),m2s(0)
#endif

      end

      subroutine getRes2Q3GWZ(p,id,fsign,pdf,Nbos,lnlo,lokt,xyuz, res)

      use globalvars, only: ldoblha

      implicit none
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
      real*8 res(1:maxnumsubproc) ! lnlo ? res(kin) : res(subproc)
      real*8 p(0:3,max_p,max_kin),pdf(-6:6,2,22),xyuz(2,2:max_kin)
      integer fsign(5), id(5), Nbos,idb(4),fsignB(4)
      logical lnlo, lokt(16)
      
      real*8 amp(0:6)

      integer family, flav(5), k, colstr
      integer FL5QCDWZ
      real*8 polcol(2000), dip(4,5,5),dipSum(16)
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
	if(nbos.gt.10) then
#ifdef WITH_QCDVV
#ifdef NINH_COMPARE
	  call ug2dggwzOwnAllHelNew(p, id, fsign,nbos,helswitch, amp)
#else
	  call ug2dggwzOwnAllHelNew(p, id, fsign,nbos,0, amp)
#endif
#endif
	else 
	  call ug2dggwOwnAllHel(p, id, fsign,nbos, amp)
	endif
      else
	  amp(0) = 0d0
      endif
      if (ldoblha) then
        return
      endif
#ifdef WITH_NLO
      if(lnlo) then
c calculate dipoles
      do i=1,16
        dipsum(i)=0d0
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
     &                     CScolAmpGNewUsePrecalc,idB,fsignB,emit1,spectB,nbos,1,dip(emit1,emit2,spect))
          dipsum(kin)=dipsum(kin)+dip(emit1,emit2,spect)
          if(ldebug) print*, 'Dipole: ',dip(emit1,emit2,spect)*alphas
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
     &                     CScolAmpGNewUsePrecalc,idB,fsignB,emitB,spectB,nbos,1,dip(emit1,emit2,spect))
          dipsum(kin)=dipsum(kin)+dip(emit1,emit2,spect)
          if(ldebug) print*, 'Dipole: ',dip(emit1,emit2,spect)*alphas
        endif
        endif
      endif
      enddo;enddo;enddo

      endif
#endif

c calculate subprocess
	do family = 1,2 ! sum over quarks
	  if (Nbos.eq. 32.or.nbos.eq.3.or.nbos.eq.31) then
	    flav(id(1)) = 2*family*fsign(1)  ! up type ( or anti up)
	    flav(id(3)) = (2*family-1)*fsign(3) 
	  else 
            flav(id(3)) = 2*family*fsign(3)  ! up type ( or anti up)
            flav(id(1)) = (2*family-1)*fsign(1) 	    
	  endif

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
        k = FL5QCDWZ(flav,id,1) 
        if(lokt(1)) then
          res(1) = res(1) + amp(0)*polcol(k)*
     1		pdf(sign1*flav(1),1,1)*pdf(sign2*flav(2),2,1)
          if(ldebug) print*,"|M_m+1|^2", amp(0)
        endif
        do kin=2,16
          if (lokt(kin)) res(kin) = res(kin) - dipsum(kin)*polcol(k)*alphas*
     1		pdf(sign1*flav(1),1,kin)*pdf(sign2*flav(2),2,kin)
        enddo
#endif
	  else
c LO subprocesses
	    do colstr = 1,6
	      k = FL5QCDWZ(flav,id,colstr) ! get process number, calc polcol and fill les houches data
! 	      print*,pdf(sign1*flav(1),1,1),pdf(sign2*flav(2),2,1)
	      res(k) = amp(colstr)*POLCOL(k)*pdf(sign1*flav(1),1,1)*pdf(sign2*flav(2),2,1)
	    enddo
	  endif
	enddo ! Quark family
c reduce process counter for next helicity configuration


      end


      subroutine getRes4Q1GWZ(p,id,fsign,SkipIdFl,pdf,Nbos,lnlo,lokt,xyuz,res)

      use globalvars, only: ldoblha

      implicit none
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/global.inc"
#include "common.inc"
      real*8 res(1:maxNumSubProc) !out: subprocess cross-sections
      real*8 p(0:3,max_p,max_kin)
      integer fsign(5), id(5),fsignB(5)
      real*8 pdf(-6:6,2,22), xyuz(2,2:max_kin)
      logical SkipIdFl(1:2),lnlo,lokt(16)
      logical calc12Id,calc34Id
      integer diffFamily,idB(4),emitB,spectB
      integer nbos

      integer FL5qcdWZ, di

      double precision  alphas, xmt, alfa, xmz, xmw,s2w, xmh
      common /bkopin/   alphas, xmt, alfa, xmz, xmw, s2w, xmh

      integer colstr,d,u,q,j,i,k
      integer iflav(1:5)
      integer family,flav
      integer kin,emit1,emit2,spect, qcdVVjjkin
      logical ldebug
      parameter (ldebug = .false.)
      real*8 amp(0:4,3:4,0:1),dip(3:4,0:2,4,5,5)
      external CScolAmpGNewUsePrecalc,CScolAmpQNew

      double precision polcol(2000)

      common /polcolz/ polcol

      if(ldebug) then
	print*,'new subprocess'
	print*,(di(id,i),i=1,5)
      endif
c precalc amplitude and dipoles
      if(lokt(1)) then
	if(nbos.gt.10) then
#ifdef WITH_QCDVV
          call amp4q1gwzNew(p,id,fsign,nbos,0,0,amp)
!           print*,amp(0,:,:)
!           call amp4q1gwzQuad(P,id,fsign,nbos,0,0,amp)
!           print*,amp(0,:,:)
!           pause
#endif
	else
	  call amp4q1gw (p,id,fsign,nbos,0,0,amp)
	endif
      else
        do i=0,4; do j=3,4
          amp(i,j,0) = 0d0
          amp(i,j,1) = 0d0
        enddo; enddo
      endif
      if (ldoblha) then
        return
      endif

#ifdef WITH_NLO
c qq dipoles
      do emit1=1,2; do emit2=3,4; do spect=1,5
        if(emit1+emit2.eq.6.or.
     $     emit1+emit2.eq.5.and.(
     $             emit1.eq.1.and..not.SkipIdFl(1)
     $         .or.emit2.eq.3.and..not.SkipIdFl(2))) then
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
     $                       CScolAmpGNewUsePrecalc,idB,fsignB,emitB,spectB,nbos,1,dip(3,0,emit1,emit2,spect))
              dip(3,0,emit1,emit2,spect)=dip(3,0,emit1,emit2,spect)*alphas
              dip(4,0,emit1,emit2,spect)=dip(3,0,emit1,emit2,spect)
              dip(4,1,emit1,emit2,spect)=dip(3,0,emit1,emit2,spect)
              dip(3,1,emit1,emit2,spect)=dip(3,0,emit1,emit2,spect)
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
     &            CScolAmpQNew,idB,fsign(1),emit1,spect,nbos,6,dip(3,0,emit1,5,spect))
          do i=3,4; do j=0,1
            dip(i,j,emit1,5,spect)=dip(i,j,emit1,5,spect)*alphas
          enddo; enddo
        endif;endif
      endif
      enddo;enddo       
#endif
        
      do j = 1,2  ! W emitting Quark family
	u = 2*j              !u,c quarks
	d = u - 1            !d,s quarks
	do flav = 1,2!1,2   ! flavor of quark line not emitting W
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
            iflav(id(5)) = 0
            
            if(ldebug) print*,'partons: ',iflav
            diffFamily = 1
            if( j-1.eq.family) diffFamily = 0
            calc12Id = (abs(iflav(id(1))).eq.q).and..not.SkipIdFl(1)
            calc34Id = (abs(iflav(id(3))).eq.q).and..not.SkipIdFl(2)
c get m2s*pdf
c   different q families
	    if(lnlo) then ! CS subtraction
	      k = FL5qcdWZ(iflav,id,1)
	      if(diffFamily.eq.1.or.calc12Id.or.calc34Id) then
		if(lokt(1)) then ! real emision ME
		  res(1) = res(1) + polcol(k)*amp(0,5-flav,diffFamily)*pdf(sign1*iflav(1),1,1)*pdf(sign2*iflav(2),2,1)
		  if(ldebug) print*,'ME ',1,amp(0,5-flav,diffFamily),polcol(k)*amp(0,5-flav,diffFamily)
     $			          *pdf(sign1*iflav(1),1,1)*pdf(sign2*iflav(2),2,1)

		endif

		do emit1=1,4 ! glu emision dipoles
		do spect=1,4
		  if(emit1.ne.spect) then
		    kin = qcdVVjjkin(id,emit1,5,spect)
		    if(kin.gt.0) then
                      if(lokt(kin)) then
			
                        res(kin) = res(kin) - polcol(k)*dip(5-flav,diffFamily,emit1,5,spect)
     $			          *pdf(sign1*iflav(1),1,kin)*pdf(sign2*iflav(2),2,kin)
                        if(ldebug) print*,'Dip',id(emit1)*100+id(5)*10+id(spect),dip(5-flav,diffFamily,emit1,5,spect),
     $                                     - polcol(k)*dip(5-flav,diffFamily,emit1,5,spect)
     $			          *pdf(sign1*iflav(1),1,kin)*pdf(sign2*iflav(2),2,kin)
                      endif
		    endif
		  endif
		enddo; enddo
		! calc q2q4 Dipole
                do spect = 1,5,2
                  kin = qcdVVjjkin(id,2,4,spect)
                  if(kin.gt.0) then
                    if(lokt(kin)) then
                      res(kin) = res(kin) - polcol(k)*dip(5-flav,diffFamily,2,4,spect)
     $			          *pdf(sign1*iflav(1),1,kin)*pdf(sign2*iflav(2),2,kin)
                      if(ldebug) print*,'Dip',id(2)*100+id(4)*10+id(spect),dip(5-flav,diffFamily,2,4,spect),
     $					- polcol(k)*dip(5-flav,diffFamily,2,4,spect)
     $			          *pdf(sign1*iflav(1),1,kin)*pdf(sign2*iflav(2),2,kin)
                    endif
                  endif
                enddo
	      endif
	      if(diffFamily.eq.0) then ! calculate dipoles q1q4 and q2q3
		emit1=0
		if(calc12Id) then
		  emit1=1
		  emit2=4
		elseif(calc34Id) then
		  emit1=2
		  emit2=3
		endif
		if(emit1.ne.0) then
		  do spect = 1,5
		    if(spect.ne.emit1.and.spect.ne.emit2) then
		      kin = qcdVVjjkin(id,emit1,emit2,spect)
		      if(kin.gt.0) then 
                        if(lokt(kin)) then
			res(kin) = res(kin) - polcol(k)*dip(5-flav,diffFamily,emit1,emit2,spect)
     $			          *pdf(sign1*iflav(1),1,kin)*pdf(sign2*iflav(2),2,kin)
			if(ldebug) print*,'Dip',id(emit1)*100+id(emit2)*10+id(spect),dip(5-flav,diffFamily,emit1,emit2,spect),
     $			                  - polcol(k)*dip(5-flav,diffFamily,emit1,emit2,spect)
     $			          *pdf(sign1*iflav(1),1,kin)*pdf(sign2*iflav(2),2,kin)
                        endif
                      endif
		    endif
		  enddo
		endif
	      endif
	    else ! born process
c       born
	      if(diffFamily.eq.1.or.calc12Id.or.calc34Id) then
		do colstr = 1,4
		  k = FL5QCDWZ(iflav,id,colstr)
		  res(k) = res(k) + polcol(k)*amp(colstr,5-flav,diffFamily)
     $			          *pdf(sign1*iflav(1),1,1)*pdf(sign2*iflav(2),2,1)
		enddo
	      endif
	    endif
	  enddo
	enddo
      enddo
      
      end


      integer function FL5qcdWZ(flav,id,colstr)
c     adopted from Terrance Figy, June 12,2001
      Implicit none
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/lha.inc"
c     COMMOM BLOCK FOR DIAGRAMMATIC LABELS
      integer id(5)
c     This commom block must be in m2s.
c
c     common block for the number of colored processes
      double precision polcol(2000)
      common /polcolz/ polcol
      integer flav(5),colstr
      integer listp
      save listp
      integer numresets
      save numresets
c      integer,save :: listp ! this is the list position
c      integer,save :: numresets
      data numresets /-1/
      data listp /0/
      external calpolcol

      if (colstr.le.-1) then  ! we are resetting
         if (colstr.gt.-10000) then
            listp = -colstr
         else
            listp = 0
         endif
         numresets = numresets+1
         numcolor = 0
      else
c increment the counter regardless of whether or not it's the 1st
c time throught
         listp = listp + 1 
         if(numcolor.lt.listp.and.numresets.eq.0) then
             call fillColoredPartonsVV3j(flav,id,colstr,listp)
             helicity(1)=9
             helicity(2)=9
             helicity(3)=9
             helicity(4)=9
             helicity(5)=9
             call calpolcol(flav,polcol(listp))
             if(listp.gt.2000) then 
               print*,"Error in ",__FILE__,__LINE__
               stop
             endif
         endif                  ! numresets.eq.0
      endif
      numcolor = max(listp,numcolor) ! common block info
      FL5qcdWZ   = listp ! outputed function
c
      end

      subroutine fillColoredPartonsVV3j(flav,id, colstr,listp1)
      implicit none
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/process.inc"

c          COMMOM BLOCK FOR DIAGRAMMATIC LABELS
      integer id(5)
      integer i,if(5),flav(5),colstr,listp
      integer ngluon,listp1
      integer pi(1:3)
      integer col(1:4)
      
      if(colstr.le.0) return

      listp = listp1 + numdecay
      lnup(listp) = 5+3*N_lepton/2
c         lidprup(listposition)=listposition so we don't need to use it.
c         scalup is same for all subprocesses
c         lwgtup(listposition) is already stored in res(fl(....)) so
c         there's
c         no need to worry about it here.
c           write(*,*) 'writing to lhepup'
      listup(1,listp)=-1        !incoming quark or gluon
      listup(2,listp)=-1
      listup(3,listp)=1         !outgoin quark or gluon
      listup(4,listp)=1
      listup(5,listp)=1

c
      do i=1,5
	if(i) = flav(i)
	if(if(i).eq.0) if(i)= 21 
      enddo
c         The particle codes:
      lidup(1,listp)=if(1)
      lidup(2,listp)=if(2)
      lidup(3,listp)=if(3)
      lidup(4,listp)=if(4)
      lidup(5,listp)=if(5)
c
c         mothers:
      do i = 3,numParticles      ! i=3,np   should be sufficient
         lmothup(1,i,listp) = 1
         lmothup(2,i,listp) = 2
      enddo
c
c     Color flow information goes here:
c
c     4Q1G and 2Q3G are the possibilities.
c     First, count the number of gluons to decide whether
c     4Q1G or 2Q3G.
      ngluon = 0
      do i = 1,5
         if(if(i).eq.21) then
            ngluon = ngluon + 1
         endif
      enddo
C reset color structure
      do i = 1,4
         col(i) = 0
      enddo

      if(ngluon.eq.1) then      ! 4Qg case
c
c
c     For 4Q1G there are 4 color structures.
c
c
c
c           0  g1                              q1---->------->------q3
c           0                                            0
c           0                                            0
c   q1  ---->-------->--------q3                         0
c                0                              q2--->------->------q4
c                0                                    0
c                0                                    0
c   q2  ---->--------->-------q4                      0  g1
c        color structure 1                        color structure 2
c
c
c           0  g1                              q1---->------->------q4
c           0                                            0
c           0                                            0
c   q1  ---->-------->--------q4                         0
c                0                              q2--->------->------q3
c                0                                    0
c                0                                    0
c   q2  ---->--------->-------q3                      0  g1
c        color structure 3                        color structure 4

c
c    The momentum flow of the gluon is outgoing.
c
         col(colstr) = 1        ! selects color structure
c
         if(id(1).lt.3) then   ! 1 and 2 are incoming pythia labels
            licolup(1,id(1),listp) = 501 ! color
            licolup(2,id(1),listp) = 0 ! anti-color
         else
            licolup(1,id(1),listp) = 0 ! color
            licolup(2,id(1),listp) =  501 ! anti-color
         endif

         if(id(2).lt.3) then
            licolup(1,id(2),listp) = 502 ! color
            licolup(2,id(2),listp) = 0 ! anti-color
         else
            licolup(1,id(2),listp) = 0 ! color
            licolup(2,id(2),listp) =  502 ! anti-color
         endif
c
         if(id(3).lt.3) then
            licolup(1,id(3),listp) = 0 ! color
            licolup(2,id(3),listp) = 502*col(1)+503*(col(2)+col(3))+501*col(4)  ! anti-color
         else
            licolup(1,id(3),listp) =502*col(1)+503*(col(2)+col(3))+501*col(4)  ! color
            licolup(2,id(3),listp) = 0 ! anti-color
         endif
c
c
         if(id(4).lt.3) then
            licolup(1,id(4),listp) = 0 ! color
            licolup(2,id(4),listp) = 503*(col(1)+col(4))+501*col(2)+502*col(3)! anti-color
         else
            licolup(1,id(4),listp) = 503*(col(1)+col(4))+501*col(2)+502*col(3)! color
            licolup(2,id(4),listp) = 0 ! anti-color
         endif
c
         if(id(5).lt.3) then
            licolup(1,id(5),listp) = 503  ! color
            licolup(2,id(5),listp) = 501*(col(1)+col(3))+502*(col(2)+col(4))! anti-color
         else
            licolup(1,id(5),listp) = 501*(col(1)+col(3))+502*(col(2)+col(4))! color
            licolup(2,id(5),listp) = 503 ! anti-color
         endif

c
c
      else                      ! 2Q3G case
c     The next case goes here.
c
c
c    q1 ----->---------------->-------q3
c             0      0       0
c             0      0       0
c             0      0       0
c           g_pi(1)  g_pi(2) g_pi(3)
c
c
         select case(colstr)
      case(1)
         pi(1) = 5
         pi(2) = 4
         pi(3) = 2
      case(2)
         pi(1) = 4
         pi(2) = 5
         pi(3) = 2
      case(3)
         pi(1) = 4
         pi(2) = 2
         pi(3) = 5
      case(4)
         pi(1) = 5
         pi(2) = 2
         pi(3) = 4
      case(5)
         pi(1) = 2
         pi(2) = 5
         pi(3) = 4
      case(6)
         pi(1) = 2
         pi(2) = 4
         pi(3) = 5
         endselect
c
         if(id(1).lt.3) then   ! 1,2 are incoming
            licolup(1,id(1),listp) = 501 ! color
            licolup(2,id(1),listp) = 0 !anti-color
         else
            licolup(1,id(1),listp) = 0 ! color
            licolup(2,id(1),listp) = 501 !anti-color
         endif
c
         if(id(3).lt.3) then   ! 1,2 are incoming
            licolup(1,id(3),listp) = 0 ! color
            licolup(2,id(3),listp) = 504 !anti-color
         else
            licolup(1,id(3),listp) = 504 ! color
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
         if(id(pi(3)).lt.3) then
            licolup(1,id(pi(3)),listp) = 504 ! color
            licolup(2,id(pi(3)),listp) = 503 !anti-color
         else
            licolup(1,id(pi(3)),listp) = 503 ! color
            licolup(2,id(pi(3)),listp) = 504 !anti-color
         endif
      endif                     ! 2Q3g or 4Q1g
      end                       ! fillColoredPartons
      integer function qcdVVjjKin(id,emit1,emit2,spect)
      implicit none
      integer emit1,emit2,spect,id(5), e1,e2
      
      if(id(emit1).gt.id(emit2)) then
	e1 = id(emit2)
	e2 = id(emit1)
      else
	e1 = id(emit1)
	e2 = id(emit2)
      endif

      if(e2.lt.3.or.spect.eq.emit1.or.spect.eq.emit2) then
	qcdVVjjKin = -1;
	return
      endif
      if(e1.lt.3) then
	if(id(spect).lt.3) then
	  qcdVVjjKin = 3*e1+e2-4
	else
	  qcdVVjjKin = 3*e1+e2+id(spect)+1
	endif
      else
	if(id(spect).lt.3) then
	  qcdVVjjKin = 3*id(spect)+e1+e2+1
	else
	  qcdVVjjKin = id(spect)+5
	endif
      endif
      end

      integer function di(id,i)
c inverse of id(4) or id(5)
      implicit none
      integer id(4),i,j
      do j=1,4
	if(id(j).eq.i) then
	  di = j
	  return
	endif
      enddo
      di=5
      end
      
      subroutine swapInt(x,y)
      implicit none
      integer x,y,a
      a=x
      x=y
      y=a
      end


ccccccc begin subroutine calpolcol ccccccccccc
      subroutine calpolcol(flavor,factor)
      Implicit none
#include "VBFNLO/utilities/global.inc"
      integer flavor(n_p)        ! in: flavors
      real*8 factor             ! out: polcol
c local variables
      integer idfl(n_p),i,j
      real*8 lfactor
      logical ldebug
      parameter (ldebug=.false.)
c
      if (ldebug) then
         print*,' '
         print*,' Enter calpolcol '
         print*,' flavor = ',flavor
      endif
      lfactor = 0.250d0         ! ave over initial state helicities
      do i = nx+1,n_p            ! initialize identical flavor counter
         idfl(i)=1
          enddo
      do i = nx+2,n_p
         do j = nx+1,i-1
            if ( flavor(i).eq.flavor(j) ) then
               idfl(j)=idfl(j)+1
             endif
          enddo
          enddo
      if (ldebug) print*,' idfl = ',idfl
      do i=nx+1,n_p-1
         lfactor = lfactor/idfl(i)
                enddo
      if (ldebug) print*,' final S factor ',0.25d0/lfactor
c The color average factor 
      do i = 1,2
         if(flavor(i).eq.21.or.flavor(i).eq.0) then
            lfactor = lfactor/8.0d0 ! gluon in i.s.
         elseif (abs(flavor(i)).le.6 ) then
            lfactor = lfactor/3.0d0 ! quark in i.s.
             endif
      enddo
      if (ldebug) then
         print*,' polcol = 1/',1/lfactor
!          pause
          endif
      factor = lfactor
      return
      end 
ccccccccccccc end subroutine ccccccccccccccccccccccccccccc

