c ************************************************************************************
c Author: Matthias Kerner
c Date: 2012/13
c ************************************************************************************

! #define NINH_COMPARE
! #define POWHEG_COMPARE
      subroutine m2s_qcdww3j(
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
      implicit none
c 
c declare input/output variables
c
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/lha.inc"
#include "common.inc"
      integer nbos
      real*8 xi(nx), p(0:3,max_p,max_kin), v(0:3,max_v,max_kin), m2s(0:max_kin)
      real*8 xyuz(2,2:max_kin)
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

      integer init/0/, I, J, k, mu
      save init
c
c declare variables for summation process
c
      integer FSIGN(5) ! diagram index
      integer id(5) ! id(diagram index) = ps index
c
c store contributions from subprocess ip in res(ip) where
      real*8 res(maxNumSubProc)
      real*8 weight, rn

c
c declare parton distribution variables
c
      real*8 pdf(-6:6,2,22)
c local varriables
      integer FL5QCDWZ, ires_min, ires_max, kin
      save ires_min, ires_max
      data ires_min,ires_max / 1, maxnumsubproc/
      real*8 mass2
      external FL5QCDWZ,mass2
      logical SkipIdFl
      logical ldebug
!       parameter (ldebug=.false.)
      integer nhel(9) ! hel with diagram index
      logical skip (0:2,2) 
      real*8 mjj2
      integer qcdVVjjkin
      
      real*8 temp3 /0d0/
      integer NN /0/


      complex*16 lepton(6,3)
      logical allsub

      complex*16 epsW
      real*8 momW
      common /QCDWWcurrents/ epsW(0:3,2,7),momW(0:4,2,7)


!       integer countProc
!       common /countP/ countProc(5)
!       do i=1,5
! 	countProc(i) =0
!       enddo
!       lokt=.true.
!       dt(1) = dtime(timearray)
!       do timeiter= 1,1000

!       ldebug = .false.


#ifdef POWHEG_COMPARE
         wmass = 80.419d0
         wwidth= 2.141d0
         GWF(1) = - sqrt(0.43907954221934414d0/2d0)
#endif


      if ( init .eq. 0 ) then
        if(nbos .eq. 33) then
         write(6,*) "QCD W+W+ jjj amplitude square information:"
        else
         write(6,*) "QCD W-W- jjj amplitude square information:"
        endif
         write(6,*) "alpha parameter of dipoles: ",alphaDip
         write(6,*) "-----------------------------------------------"
         write(6,*) " "
c         call CSGeneralSetAlpha(alphaDip)
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
      end if

      do i=ires_min,ires_max
	res(i) = 0d0
      enddo
      
      allsub=.true.

	
c 2->7
#ifdef NINH_COMPARE
      lokt=.true.
      lokt(8:16) = .true.
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
      print*,'a'
#endif
c     
c call PDF subroutine in order to determine parton
c distributions in the incoming (anti)protons.
c
      q_sf = sqrt(mufsq(1,1))
      call pdfproton( xi(1), q_sf, pdf(-6,1,1) )        ! x1 f_a(x1)

      q_sf = sqrt(mufsq(2,1))
      call pdfproton( xi(2), q_sf, pdf(-6,2,1) )        ! x2 f_a(x2)

      do i=-6,6
	pdf(i,1,1)=pdf(i,1,1)/xi(1)
	pdf(i,2,1)=pdf(i,2,1)/xi(2)
      enddo
      if(lnlo) then
	do i =2,16
	  if(lokt(i)) then
	    if(mufsq(1,1).ne.mufsq(1,i).or.mufsq(2,1).ne.mufsq(2,i)) then
	      q_sf = sqrt(mufsq(1,i))
	      call pdfproton( xi(1), q_sf, pdf(-6,1,i) )        ! x1 f_a(x1)
	      q_sf = sqrt(mufsq(2,i))
	      call pdfproton( xi(2), q_sf, pdf(-6,2,i) )        ! x2 f_a(x2)
              do j=-5,5
                pdf(j,1,i) = pdf(j,1,i)/xi(1)
                pdf(j,2,i) = pdf(j,2,i)/xi(2)
              enddo
	    else
	      do j = -6,6
		pdf(j,1,i) = pdf(j,1,1)
		pdf(j,2,i) = pdf(j,2,1)
	      enddo
	    endif
	  endif
	enddo
      endif


c calc W currents
      do kin=1,7
      if(lokt(kin).or.kin.eq.1) then
	do j=1,2
	  if(nbos.eq.33) then
	    CALL IXXXXX(v(0,2*j-1,kin),0d0,1,-1,lepton(1,1))        
	    CALL OXXXXX(v(0,2*j  ,kin),0d0,-1,1,lepton(1,2))
	    CALL JIOXXX(lepton(1,1) ,lepton(1,2) ,GWF ,WMASS   ,WWIDTH  ,lepton(1,3))
	  else
	    CALL OXXXXX(v(0,2*j-1,kin),0d0,-1,1,lepton(1,1))        
	    CALL IXXXXX(v(0,2*j  ,kin),0d0,1,-1,lepton(1,2))
	    CALL JIOXXX(lepton(1,2) ,lepton(1,1) ,GWF ,WMASS   ,WWIDTH  ,lepton(1,3))
	  endif
	  do mu=0,3
	    epsW(mu,j,kin)=lepton(mu+1,3)*GWF(1)
	    momW(mu,j,kin)=v(mu,2*j-1,kin)+v(mu,2*j,kin)
	  enddo
	  momW(4,j,kin) = mass2(momW(0,j,kin))
	enddo
      endif
      enddo
      
      call CSColAmpQNewReset()

#ifdef NINH_COMPARE
        alphas=als(1,1)

        G = sqrt(4d0*pi*alphas)
        gg(1) = -G
        gg(2) = -G
        do i=-6,6
          pdf(i,1,:) =1d0+0.01d0*i
          pdf(i,2,:) =1d0+0.02d0*i
        enddo
        xi(1) =1d0
        xi(2) =1d0        
#endif
c 2->5  

!         p(0:3,1,1) = (/ 500.000000000000d0, 0.000000000000000d0,0.000000000000000d0,500.000000000000d0/)
!         p(0:3,2,1) = (/ 500.000000000000d0, 0.000000000000000d0,0.000000000000000d0,-500.000000000000d0/)
!         momW(0:3,1,1) = (/ 397.866392966475d0,-297.523578137646d0,32.1947053257073d0, 249.553700932918d0/)
!         momW(0:3,2,1) = (/ 122.262841651602d0,-83.8081688340796d0,-30.1075145326538d0,23.5391684234863d0/)
!         p(0:3,3,1) = (/ 10.0696722001341d0,5.72361968659522d0,-2.06452467121156d0, 8.02347890776386d0/)
!         p(0:3,4,1) = (/0.485024026870584d0,0.146631598351091d0,-2.266612184189497d-2,  0.461772376747906d0/)
!         p(0:3,5,1) = (/ 469.316069154918d0,375.461495686779d0,0d0,  -281.578120640916d0/)
! 	v(:,1,1) = momw(0:3,1,1)
! 	v(:,2,1) = momw(0:3,2,1)
! 	v(:,3:4,1) = 0d0
! 
! 	call qcdVVJJptildemake(p,xyuz,v)
! 	momw(0:3,1:2,1:7)=v(:,1:2,1:7)
! 
!         do i=1,7
! 	  momw(4,1,i) = sqrt(mass2(momw(0,1,i)))
! 	  momw(4,2,i) = sqrt(mass2(momw(0,2,i)))
! 	enddo
! 
!          print*," parton momenta "
!          do i = 1, 5
!             write(6,10) i, p(0,i,1), p(1,i,1), p(2,i,1), p(3,i,1)
!          end do
!          print*," boson momenta "
!          do i = 1, 2
!             write(6,10) i, momw(0,i,1), momw(1,i,1), momw(2,i,1), momw(3,i,1)
!          end do
! 	do kin=2,16
! 	print*
! 	print*,kin
! 	         print*," parton momenta "
!          do i = 1, 4
!             write(6,10) i, p(0,i,kin), p(1,i,kin), p(2,i,kin), p(3,i,kin)
!          end do
!          print*," boson momenta "
!          do i = 1, 2
!             write(6,10) i, momw(0,i,kin), momw(1,i,kin), momw(2,i,kin), momw(3,i,kin)
!          end do
!          enddo
!          print*,'sum momenta',(sum( (/ p(mu,3:5,7),v(mu,1:4,7),-p(mu,1:2,7) /) ),mu=0,3)
!       endif
! 
!       debugsum = 0d0
!       do hel1=-1,1
!       do hel2=-1,1
! !         print*,'helicities: ',hel1,hel2
!         do kin=1,7
! 	  call VXXXXX(momW(0,1,kin),WMASS,hel1,1,lepton(1,1))
! 	  epsW(0:3,1,kin)=lepton(1:4,1)
! 	  call VXXXXX(momW(0,2,kin),WMASS,hel2,1,lepton(1,1))
! 	  epsW(0:3,2,kin)=lepton(1:4,1)
!   !         call helvec(momW(0,1,1),1,hel1,epsW(0,1,1)) ! mom, outgoing=1, hel, epsilon
!   !         call helvec(momW(0,2,1),1,hel2,epsW(0,2,1))
!   !         print*,epsW(0:3,1,1)
! 	  epsW(0:3,1,kin) = epsW(0:3,1,kin)*GWF(1)/sqrt(2d0) ! Wff coupling and identical factor
! 	  epsW(0:3,2,kin) = epsW(0:3,2,kin)*GWF(1)
! 	enddo
! 	m2s=0d0
! 	res=0d0
! 	borncache4q=-1d0
      
      ires_min = FL5QCDWZ(id,id, -100000) +1
 
C--------begin subprocess evaluation-----------------------------------

c*****  1.  ********
c**Diagram**********  q1 q2 ---> q3 q4 g W Z   **********************
c**Phasespace*******  1   2      3   4 g

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


      SkipIdFl = .false.

      if(sub_number.eq.1.or.allsub) call getRes4Q1GWW(p,xi,id,fsign,SkipIdFl,pdf,nbos,lnlo,lokt,xyuz,res )
c*****  2.  ********
c**Diagram**********  qb3 qb4 ---> qb1 qb2 g W Z   **********************
c**Phasespace*******  1   2        3   4

      id(1) = 4
      id(2) = 3
      id(3) = 1
      id(4) = 2
      id(5) = 5

      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = -1
      fsign(4) = -1
      fsign(5) = 1

      if(sub_number.eq.2.or.allsub) call getRes4Q1GWW(p,xi,id,fsign,SkipIdFl,pdf,nbos,lnlo,lokt,xyuz,res )
      

c*****  3.  ********
c**Diagram**********  q1 qb3 ---> qb2 q4 g W Z   **********************
c**Diagram**********  qb3 q1 ---> qb2 q4 g W Z   **********************
c**Phasespace*******  1   2        3   4       **********************

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

c reset borncache

      if(sub_number.eq.3.or.allsub) call getRes4Q1GWW(p,xi,id,fsign,SkipIdFl,pdf,nbos,lnlo,lokt,xyuz,res )
! switch initial state
      id(1) = 2
      id(3) = 1
      
      if(sub_number.eq.3.or.allsub) call getRes4Q1GWW(p,xi,id,fsign,SkipIdFl,pdf,nbos,lnlo,lokt,xyuz,res )


c*****  4.  ********
c**Diagram**********  q1 qb4 ---> q3 qb2 g W Z   **********************
c**Diagram**********  qb4 q1 ---> q3 qb2 g W Z   **********************
c**Phasespace*******  1   2        3   4

      id(1) = 1
      id(2) = 3
      id(3) = 4
      id(4) = 2
      id(5) = 5

c  same as crossing 2 for identical flavours
      SkipIdFl=.true.

      fsign(1) =  1
      fsign(2) = -1
      fsign(3) =  1
      fsign(4) = -1
      fsign(5) =  1

      if(sub_number.eq.3.or.allsub) call getRes4Q1GWW(p,xi,id,fsign,SkipIdFl,pdf,nbos,lnlo,lokt,xyuz,res )
! switch initial state
      id(1) = 2
      id(4) = 1

      if(sub_number.eq.3.or.allsub) call getRes4Q1GWW(p,xi,id,fsign,SkipIdFl,pdf,nbos,lnlo,lokt,xyuz,res )

c*****  5.  ********
c**Diagram**********  q1 g ---> q3 qb2 q4 W Z   **********************
c**Diagram**********  g q1 ---> q3 qb2 q4 W Z   **********************
c**Phasespace*******  1   2     3   4  5

      id(1) = 1
      id(2) = 3
      id(3) = 4
      id(4) = 5
      id(5) = 2

      SkipIdFl=.false.

      fsign(1) =  1
      fsign(2) = -1
      fsign(3) =  1
      fsign(4) =  1
      fsign(5) = -1

      if(sub_number.eq.4.or.allsub) call getRes4Q1GWW(p,xi,id,fsign,SkipIdFl,pdf,nbos,lnlo,lokt,xyuz,res )
! switch initial state
      id(1) = 2
      id(5) = 1

      if(sub_number.eq.4.or.allsub) call getRes4Q1GWW(p,xi,id,fsign,SkipIdFl,pdf,nbos,lnlo,lokt,xyuz,res )

c*****  6.  ********
c**Diagram**********  qb3 g  ---> qb1 qb2 q4 W Z   **********************
c**Diagram**********  g  qb3 ---> qb1 qb2 q4 W Z   **********************
c**Phasespace*******  1   2        3   4  5

      id(1) = 3
      id(2) = 4
      id(3) = 1
      id(4) = 5
      id(5) = 2

c  same as crossing 2 for identical flavours
      SkipIdFl=.false.

      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = -1
      fsign(4) =  1
      fsign(5) = -1

      if(sub_number.eq.5.or.allsub) call getRes4Q1GWW(p,xi,id,fsign,SkipIdFl,pdf,nbos,lnlo,lokt,xyuz,res )
! switch initial state
      id(3) = 2
      id(5) = 1

      if(sub_number.eq.5.or.allsub) call getRes4Q1GWW(p,xi,id,fsign,SkipIdFl,pdf,nbos,lnlo,lokt,xyuz,res )

c
c
c---end of subprocess evaluation-------------------------------
c
c

	if(lnlo) then
	  ires_max = 17
	  ires_min = 1
	else
          ires_max = numcolor
	endif
      if(init.eq.0) then
! 	if(ldebug) print*, "Number of (color) subprocesses: ", ires_max
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
         if(weight.ge.rn) then
            iprocess = i
         endif
         i = i+1
      enddo

!       if(iprocess.eq.0..and..not.lnlo) then
!          print*, 'error in subprocess selection'
!       endif
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
	do i=1,17
	 print*,'m2s( ',i,' ) ', m2s(i), lokt(i)
	enddo
      endif

!       print*,m2s(0)-m2s(1), m2s(1)
!       debugsum=debugsum+m2s(0)
!       enddo
!       enddo
!       m2s(0) = debugsum
!       print*,'sum', debugsum
!       print*,'disabled born cache?'
!       pause

!       enddo !timeiter
!       dt(1) = dtime(timearray)
!       print*,dt(1)
!       pause
      end


      subroutine getRes4Q1GWW(p,xi,id,fsign,SkipIdFl,pdf,Nbos,lnlo,lokt,xyuz,res)
c hel summation, if rn < 0
      implicit none
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/global.inc"
      real*8 res(1:maxNumSubProc) !out: subprocess cross-sections
      real*8 p(0:3,max_p,max_kin),xi(2)
      integer fsign(5), id(5)
      real*8 pdf(-6:6,2,22), xyuz(2,2:max_kin)
      logical lnlo,lokt(max_kin)
      logical idFl, skipIdFl
      integer nbos,difffamily

      integer FL5qcdWZ, di,idb(4)


      double precision  alphas, xmt, alfa, xmz, xmw,s2w, xmh
      common /bkopin/   alphas, xmt, alfa, xmz, xmw, s2w, xmh

      integer j,i,k,colstr
      integer iflav(5)
      integer kin,emit,spect, qcdVVjjkin
      logical ldebug
      parameter (ldebug = .false.)
      real*8 amp(0:4,0:2),dip(3:4,0:2,4,5,5)
      double precision polcol(2000)

      common /polcolz/ polcol
      external CScolAmpQNew

c  nloContrib 11: only calc RE+Dip
c  nloContrib 12: only calc finColl
      
      
      if(ldebug) then
	print*,'new subprocess'
	print*,(di(id,i),i=1,5)
      endif

c precalc amplitude and dipoles
      if(lokt(1)) then
          call amp4q1gwwHelSum(p,id,fsign,amp)
      else
        do i=0,4; do j=0,2
          amp(i,j) = 0d0
        enddo; enddo
      endif
      
!       if(ldebug.and.lokt(1)) print*,'amp',amp
            
      
      if(lnlo) then
      do emit=1,4
      do spect=1,4
        if(emit.ne.spect) then
          kin = qcdVVjjKin(id,emit,5,spect)
          if(kin.gt.0) then
            if(lokt(kin)) then
!               if(ldebug) print*,'emit,spect, kin',emit,5,spect,kin
              do i=1,4
                idB(i) = id(i)
              enddo
              if(id(5).le.2) then
                idB(emit) = id(5)
                do i=1,4
                  if(idb(i).gt.id(emit)) idB(i) = idb(i) -1
                enddo
              endif
              call CSgeneral(p,kin,xyuz,"qg",id(emit),id(5),id(spect),0,0,
     &            CScolAmpQNew,idB,fsign(1),emit,spect,nbos,6,dip(3,0,emit,5,spect))
              dip(3,0,emit,5,spect) = dip(3,0,emit,5,spect)*alphas
              dip(3,1,emit,5,spect) = dip(3,1,emit,5,spect)*alphas
!               if(ldebug) print *,dip(3,0,emit,5,spect),dip(3,1,emit,5,spect)
            endif
          endif
        endif
      enddo	  
      enddo
      endif

      do j=0,2,2
      do k=0,2,2
	IdFl = j.eq.k
	difffamily = 1
	if(idFl) difffamily = 0

	if(.not.(skipIdFl.and.IdFl)) then
c fill flavor information for pdf and lha output
	  if(Nbos.eq.33) then  ! w+ case
	      iflav(id(1)) = (2+j)*fsign(1)
	      iflav(id(3)) = (1+j)*fsign(3)
	      iflav(id(2)) = (2+k)*fsign(2)
	      iflav(id(4)) = (1+k)*fsign(4)
	      iflav(id(5)) = 0
	  elseif(Nbos.eq.44) then ! w- case
	      iflav(id(1)) = (1+j)*fsign(1)
	      iflav(id(3)) = (2+j)*fsign(3)
	      iflav(id(2)) = (1+k)*fsign(2)
	      iflav(id(4)) = (2+k)*fsign(4)
	      iflav(id(5)) = 0
	  else
	      stop
	  endif
	  
	  if(lnlo) then 
	    i = FL5qcdWZ(iflav,id,1)
	    if(lokt(1)) then 
	      res(1) = res(1) +polcol(i)*amp(0,difffamily)
     $		      *pdf(sign1*iflav(1),1,1)*pdf(sign2*iflav(2),2,1)
	    endif
            
	      
	    do emit=1,4
	    do spect=1,4
	      if(emit.ne.spect) then
		kin = qcdVVjjkin(id,emit,5,spect)
		if(kin.gt.0) then
		  if(lokt(kin)) then
		    res(kin) = res(kin) - polcol(i)*dip(3,difffamily,emit,5,spect)
     $			          *pdf(sign1*iflav(1),1,kin)*pdf(sign2*iflav(2),2,kin)
		  endif
		endif
	      endif
	    enddo;enddo
	  else
            do colstr=1,4
              i = FL5qcdWZ(iflav,id,colstr)
              if(idFl) then  ! same q families
                res(i) = POLCOL(i) * amp(colstr,0)
     $		      *pdf(sign1*iflav(1),1,1)*pdf(sign2*iflav(2),2,1)
                if(ldebug) print*,i, iflav,res(i)
              else
                res(i) = POLCOL(i) * amp(colstr,1)
     $		      *pdf(sign1*iflav(1),1,1)*pdf(sign2*iflav(2),2,1)
                if(ldebug) print*,i, iflav,res(i)
              endif	  
            enddo
	  endif ! LO
! 	  print*,iflav,res(1),sum(res(2:17))
	endif ! id FL
      enddo
      enddo
      return
      end


      subroutine amp4q1gwwHelSum(P,id,fsign,res)
      implicit none
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
      real*8 P(0:3,max_p,max_kin),res(0:4,0:2)
      integer id(5),fsign(5),pol
      
      real*8 pbar(0:3,5),pdia(0:3,5),momTemp(0:4)
      real*8 mompsiW(0:4,4,2),mompsiG(0:4,4),momCurrij(0:4,1:3,2:4,2)
      complex*16 psi(2,-1:1,4), psiW(2,4,2),psiG(2,4,2),psiWG(2,4,2,2)
      complex*16 psiTemp(2,-1:1),amp(4,2)
      real*8 epsG(0:3,2)
      complex*16 gluWij(0:3,1:3,2:4,2),gluWGij(0:3,1:3,2:4,2,2,2),temp
      complex*16 tempCurr(0:3,-1:1,4),amp_ort
      complex*16 s1c, dotcc, dotrc
      real*8 mass2,qjj2,dotrr
      
      real*8 tempR,colsq(0:4,2)
      
      integer i,j,k,l,m,mu,tempi4(4),in,out
      integer minhel, maxhel
      external s1c,dotcc,dotrc,mass2,qjj2,minhel,maxhel,dotrr

      complex*16 epsW
      real*8 momW
      common /QCDWWcurrents/ epsW(0:3,2,7),momW(0:4,2,7)

      integer color_Z(2,2,4),color_C(2,2)
      data (color_Z(2,2,i),i=1,4) / 2, 2, 2, 2 /
      data (color_Z(1,2,i),i=1,4) / 2,-4,-2, 4 /
      data (color_Z(2,1,i),i=1,4) /-4, 2, 4,-2 /
      data (color_Z(1,1,i),i=1,4) /-4,-4,-4,-4 /
      data color_C(2,2), color_C(1,2), color_C(2,1), color_C(1,1) / 5, 3, 3, 1 /
  
!       print*,sum(epsw),sum(momW)!sum(P(0:3,1:5,1)),id,fsign,pol,res
c new numbering convention
c in   q1 >>> #### >>> q2  out
c             ####          , WZ
c             ####
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
      
c get bra and ket 
      tempi4(1) = fsign(1)
      tempi4(2) = fsign(3)
      tempi4(3) = fsign(2)
      tempi4(4) = fsign(4)
      call psi0m(4,pbar,tempi4,psi)

      call polvec(pbar(0,5),1,epsG(0,1))
      call polvec(pbar(0,5),2,epsG(0,2))
      
c W and glu emissions emission
      do in=1,3,2
      do pol=1,2
        call ket2r(psi(1,-1,in),.true.,pdia(0,in),-1,pdia(0,5),epsG(0,pol),psiG(1,in,pol),mompsiG(0,in))
      enddo
      do j=1,2

        call ket2c(psi(1,-1,in),.true.,pdia(0,in),-1,momW(0,j,1),epsW(0,j,1),psiW(1,in,j),mompsiW(0,in,j))

        do pol=1,2
          call ket2c(psiG(1,in,pol),.false.,momPsiG(0,in),-1,momW(0,j,1),epsW(0,j,1),psiWG(1,in,j,pol),momtemp)
          call ket2r(psiW(1,in,j),.false.,momPsiW(0,in,j),-1,pdia(0,5),epsG(0,pol),psitemp(1,1),momtemp)
          psiWG(1,in,j,pol) = psiWG(1,in,j,pol)+psitemp(1,1)
          psiWG(2,in,j,pol) = psiWG(2,in,j,pol)+psitemp(2,1)
        enddo
      enddo
      enddo

      do out=2,4,2
      do pol=1,2
        call bra2r(psi(1,-1,out),.true.,pdia(0,out),-1,pdia(0,5),epsG(0,pol),psiG(1,out,pol),mompsiG(0,out))
      enddo
      do j=1,2
        call bra2c(psi(1,-1,out),.true.,pdia(0,out),-1,momW(0,j,1),epsW(0,j,1),psiW(1,out,j),mompsiW(0,out,j))

        do pol=1,2
          call bra2c(psiG(1,out,pol),.false.,momPsiG(0,out),-1,momW(0,j,1),epsW(0,j,1),psiWG(1,out,j,pol),momtemp)
          call bra2r(psiW(1,out,j),.false.,momPsiW(0,out,j),-1,pdia(0,5),epsG(0,pol),psitemp(1,1),momtemp)
          psiWG(1,out,j,pol) = psiWG(1,out,j,pol)+psitemp(1,1)
          psiWG(2,out,j,pol) = psiWG(2,out,j,pol)+psitemp(2,1)
        enddo
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
        do pol=1,2
c G attached to out
          call curr(-1,psiG(1,out,pol),psiW(1,in,i),tempCurr(0,-1,1))
          call curr(-1,psiWG(1,out,i,pol),psi(1,-1,in),tempCurr(0,-1,2))
          do mu=0,3
            gluWGij(mu,in,out,i,pol,2) = tempCurr(mu,-1,1)+tempCurr(mu,-1,2)
          enddo   
c G attachaed to in
          call curr(-1,psi(1,-1,out),psiWG(1,in,i,pol),tempCurr(0,-1,1))
          call curr(-1,psiW(1,out,i),psiG(1,in,pol),tempCurr(0,-1,2))
          do mu=0,3
            gluWGij(mu,in,out,i,pol,1) = tempCurr(mu,-1,1)+tempCurr(mu,-1,2)
          enddo
        enddo
      enddo
      enddo
      enddo
   
      do i=1,3,2
      do j=2,4,2
      do k=1,2
        do mu=0,3
          momCurrij(mu,i,j,k) = mompsiW(mu,j,k)-pdia(mu,i)
        enddo; 
        momCurrij(4,i,j,k) = mass2(momCurrij(0,i,j,k))
      enddo; enddo; enddo
      
      do i=0,4; do j=0,2
        res(i,j) = 0d0
      enddo; enddo

      do pol=1,2

c calc amplitudes      
      do i=1,2 ! pos of G
        ! <1|G|2><3|4>
        amp(i,1)   = dotcc(gluWGij(0,1,2,1,pol,i),gluWij(0,3,4,2))/momCurrij(4,3,4,2)
     1              +dotcc(gluWGij(0,1,2,2,pol,i),gluWij(0,3,4,1))/momCurrij(4,3,4,1)
        ! <1|2><3|G|4>
        amp(i+2,1) = dotcc(gluWGij(0,3,4,1,pol,i),gluWij(0,1,2,2))/momCurrij(4,1,2,2)
     1              +dotcc(gluWGij(0,3,4,2,pol,i),gluWij(0,1,2,1))/momCurrij(4,1,2,1)
        ! <1|G|4><3|2>
        j=i
        if(j.eq.2) j=4
        amp(j,2)   = dotcc(gluWGij(0,1,4,1,pol,i),gluWij(0,3,2,2))/momCurrij(4,3,2,2)
     1              +dotcc(gluWGij(0,1,4,2,pol,i),gluWij(0,3,2,1))/momCurrij(4,3,2,1)
        ! <1|4><3|G|2>
        j=4-i
        amp(j,2)   = dotcc(gluWGij(0,3,2,1,pol,i),gluWij(0,1,4,2))/momCurrij(4,1,4,2)
     1              +dotcc(gluWGij(0,3,2,2,pol,i),gluWij(0,1,4,1))/momCurrij(4,1,4,1)
      enddo

c ggg Vertex contributions
      temp=0d0
      do i=1,2
        j=3-i
        temp = temp+ 1d0/momCurrij(4,1,2,i)/momCurrij(4,3,4,j)
     &          *( dotrc(pdia(0,5),gluWij(0,1,2,i))*dotrc(epsG(0,pol),gluWij(0,3,4,j))
     &            +dotrr(momCurrij(0,3,4,j),epsG(0,pol))*dotcc(gluWij(0,1,2,i),gluWij(0,3,4,j))
     &            -dotrc(pdia(0,5),gluWij(0,3,4,j))*dotrc(epsG(0,pol),gluWij(0,1,2,i)))
      enddo
      amp(2,1) = amp(2,1) - 2d0*temp
      amp(1,1) = amp(1,1) + 2d0*temp

      temp=0d0
      do i=1,2
        j=3-i
        temp = temp+ 1d0/momCurrij(4,1,4,i)/momCurrij(4,3,2,j)
     &          *(-dotrc(pdia(0,5),gluWij(0,1,4,i))*dotrc(epsG(0,pol),gluWij(0,3,2,j))
     &            +dotrr(momCurrij(0,1,4,i),  epsG(0,pol))*dotcc(gluWij(0,1,4,i),gluWij(0,3,2,j))
     &            +dotrc(pdia(0,5),gluWij(0,3,2,j))*dotrc(epsG(0,pol),gluWij(0,1,4,i)))
      enddo
      amp(2,2) = amp(2,2) - 2d0*temp
      amp(3,2) = amp(3,2) + 2d0*temp
c square amplitudes
!       print*,amp*gg(1)**3
c for lhe output
      colsq(0,1) = 0d0
      colsq(0,2) = 0d0
      do k=1,4; do l=1,2
        colsq(k,l) = dble(amp(k,l)*dconjg(amp(k,l)))
      enddo; enddo
c    diff q families
      do k=1,2; do l=1,2
        amp_ort = 0d0
        do m=1,4
          amp_ort = amp_ort + amp(m,1)*color_Z(k,l,m)
        enddo
        tempR=dble(amp_ort*dconjg(amp_ort))*color_c(k,l)/36d0
        res(0,1) = res(0,1) + tempR
        colsq(0,1) = colsq(0,1) + tempR ! for lhe output
c    same q family
          temp = (k*2)-3
          amp_ort = amp_ort - temp*amp(1,2)*color_Z(k,l,3) ! change of initial state changes nnumbering of amplitudes!
          amp_ort = amp_ort - temp*amp(2,2)*color_Z(k,l,2)
          amp_ort = amp_ort - temp*amp(3,2)*color_Z(k,l,1)
          amp_ort = amp_ort - temp*amp(4,2)*color_Z(k,l,4)
          tempR=dble(amp_ort*dconjg(amp_ort))*color_c(k,l)/36d0
          res(0,0) = res(0,0) + tempR
          colsq(0,2) = colsq(0,2) + tempR ! for lhe output
      enddo;enddo
c lhe output
          tempR=colsq(1,1)+colsq(2,1)+colsq(3,1)+colsq(4,1)
          if(tempR.ne.0d0) then
            res(1,1) = res(1,1)+ colsq(0,1)*(colsq(1,1)+colsq(4,1))/tempR
            res(2,1) = res(2,1)+ colsq(0,1)*(colsq(2,1)+colsq(3,1))/tempR
          endif
          tempR=tempR+colsq(1,2)+colsq(2,2)+colsq(3,2)+colsq(4,2)
          res(1,0) = res(1,0)+ colsq(0,2)*(colsq(1,1)+colsq(4,1))/tempR
          res(2,0) = res(2,0)+ colsq(0,2)*(colsq(2,1)+colsq(3,1))/tempR
          res(3,0) = res(3,0)+ colsq(0,2)*(colsq(1,2)+colsq(2,2))/tempR
          res(4,0) = res(4,0)+ colsq(0,2)*(colsq(3,2)+colsq(4,2))/tempR
      enddo ! pol
      do i=0,4; do j=0,2
        res(i,j) = res(i,j) *G**6
      enddo; enddo

      end
