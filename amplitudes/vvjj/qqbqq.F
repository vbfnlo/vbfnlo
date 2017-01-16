C********************  qqbqqi ***********************************

      subroutine qqbqqi(bos,lflavr)
          use globalvars, only: ldoblha
      implicit none
      integer bos, lflavr(5:6)
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

C     QQBQQi initializes the calculation of matrix elements for 
c     electroweak single weak boson production by quark quark scattering
c     i.e. weak boson fusion
c
C        q1 q3    ---->   q2 q4 V (gluon) ,   V ---> f5-bar f6
c
C     and crossing related processes. Pauli interference terms for
c     identical fermions are neglected. In particular, only the
c     t-channel exchange of elctroweak bosons is considered.
c     s-channel production of weak boson pairs is NOT implemented.
C     
C     QQBQQI must be called first to initialize the desired
c     production process. This subroutine works with both
c
C   subroutine qqZqq              for the lowest order 
c              qqwpqq             processes
c              qqwmqq
c
c   and 
c
C   subroutine qqZqqj             for the real emission corrections
c              qqwpqqj            to these processes, i.e. for the above
c              qqwmqqj            processes with one additional gluon 
c
C     INPUT:
C     BOS = 2              Z,(gamma)  production
C     BOS = 3              W+  production
C     BOS = 4              W-  production
C     LFLAVR(5:6)          the flavor of the BOS decay products in
c                            the standard notation, i.e.the fermion
c                            indices are:  
C                            1: neutrino of electron (muon,tau)
C                            2: electron (muon,tau)
C                            3: up quark (charm,top)
c                            4: down quark (strange,bottom)
C
C     The main task of QQBQQI is to set up the products of coupling
c     constants needed in Feynman graphs a, ..., g which are stored
c     in fcpl(4,4,8,6,7) and which are calculated from output of
c     KOPPLN

C
C  Here fcpl(vl,vh,sig,i,j) contains the coupling constant factors. 
C     vl  = 1,..,4    the lower one of the boson indices i of the V_i
C     vh  = 1,..,4    the higher one of the boson indices i of the V_i
C                     vl = 2 e.g. means V_lower is a Z-boson
c		      Note:  the subscript on vi is determined by its 
c		      origin, not its identity.  V1 is the produced
c                     boson (connected to the 5-6 fermion line),
c		      V2 is the scattered boson connected to the 1-2
c                     fermion line and V3 is the scattered boson
c                     connected to the 3-4 fermion line.	
C     sig = 1,..,8    identifies the (maximally) 8 different helicity
C                     combinations
C     i   = 1,..,6    numbers the possible subprocesses (4 for W
c                     production and 6 for Z production)
C     j   = 1,..,7    the 7 Feynman graphs
C  
      double precision fcpl(4,4,8,6,7)
      integer vlmin(6,7), vlmax(6,7), vhmin(6,7), vhmax(6,7), 
     $     nprmax, sigmax, sig1(8,4), sig3(8,4), sig5(8,4)
      common /cqqbqq/ fcpl, nprmax, sigmax,
     1     vlmin, vlmax, vhmin, vhmax, sig1, sig3, sig5
      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)
      
      integer b2min(6,4), b2max(6,4), b3min(6,4), b3max(6,4),
     $     b1min(4), b1max(4), fl1(6,4), fl2(6,4),
     $     fl3(6,4), fl4(6,4), zsig1(8,2), wsig1(3,4), 
     $     zsig3(8,2), wsig3(3,4), zsig5(8,2), wsig5(3,4), 
     $     i, imin, j, vl, vh, v1, v2, v3, iv,  
     $     isig, isig1, isig3, isig5
      save b2min, b2max,b3min,b3max,fl1,fl2,fl3,fl4,
     1     zsig1,zsig3,zsig5,wsig1,wsig3,wsig5
      logical rejgr, ldebug
      parameter (ldebug=.false.)
      
      data b1min /1,1,4,3/, b1max /1,2,4,3/
      
      data b2min /1,1,1,1,3,4,  1,1,1,1,3,4,  1,1,3,3,9,9,  1,1,4,4,9,9/
      data b2max /2,2,2,2,3,4,  2,2,2,2,3,4,  2,2,3,3,0,0,  2,2,4,4,0,0/
      data b3min /1,1,1,1,4,3,  1,1,1,1,4,3,  3,3,1,1,9,9,  4,4,1,1,9,9/
      data b3max /2,2,2,2,4,3,  2,2,2,2,4,3,  3,3,2,2,0,0,  4,4,2,2,0,0/
c                             |<--gam,Z--->|  |<---W+--->|  |<---W- --->|
      
      data fl1 /3,3,4,4,3,4,  3,3,4,4,3,4,  3,4,3,3,0,0,  3,4,4,4,0,0/
      data fl2 /3,3,4,4,4,3,  3,3,4,4,4,3,  3,4,4,4,0,0,  3,4,3,3,0,0/
      data fl3 /3,4,3,4,4,3,  3,4,3,4,4,3,  3,3,3,4,0,0,  4,4,3,4,0,0/
      data fl4 /3,4,3,4,3,4,  3,4,3,4,3,4,  4,4,3,4,0,0,  3,3,3,4,0,0/
c	                    |<--gam,Z--->|  |<---W+--->|  |<---W- --->|
      
      data zsig1 /-1,-1,-1,-1, 1, 1, 1, 1,    -1,-1,6*0/
      data zsig3 /-1,-1, 1, 1,-1,-1, 1, 1,    -1,-1,6*0/
      data zsig5 /-1, 1,-1, 1,-1, 1,-1, 1,    -1, 1,6*0/
C Feynman graphs        1,...,6                    7
c                   No 3 boson vertex		3 boson vertex
      
      data wsig1 /-1, 1,0,   -1,0,-1,   -1, 1,-1,    -1, 1,-1/
      data wsig3 /-1,-1,0,   -1,0, 1,   -1,-1, 1,    -1,-1, 1/
      data wsig5 /-1,-1,0,   -1,0,-1,   -1,-1,-1,    -1,-1,-1/
C Feynman graphs    1,2         3,4        5,6           7
      
c     Note: when determining the possible values of sigi in the above 
c     data statements, sigi = +,- unless it only couples to a W, then
c     sigi = - only.  Zeros occur when certain graphs do not exist,
c     such as there is no 3 boson vertex except 2 W"s and 1 gamma or
c     Z.

      if (bos.eq.2.or.bos.eq.1) then        ! Z or A production 
         nprmax = 6             ! # of flavor combinations, uudd, etc.
         sigmax = 8             ! max # of combinations
                                !                  (sig1,sig2,sig3)
         
C     6 subprocesses for Z production. Fix sig1,3,5 ranges:
         do i = 1,8             ! 8 helicity combinations
            do j = 1,3          ! makes 3 copies of zsig1,2,3(i,1):
               sig1(i,j) = zsig1(i,1) ! 	j = 1 for graphs 1, 2
               sig3(i,j) = zsig3(i,1) ! 	j = 2 for graphs 3, 4
               sig5(i,j) = zsig5(i,1) ! 	j = 3 for graphs 5, 6
            enddo			
            sig1(i,4) = zsig1(i,2) ! Make 1 copy of zsig1,2,3(i,2):
            sig3(i,4) = zsig3(i,2) ! 	j = 4 for graph 7
            sig5(i,4) = zsig5(i,2)
         enddo
         
      elseif (bos.eq.3 .or.  bos.eq.4) then
         nprmax = 4
         sigmax = 3
C  4 subprocesses for W+ or W- production. Fix sig1,3,5 ranges:
         do j = 1,4
            do i = 1,8          ! re-zero all entries since some
               sig1(i,j) = 0    ! will not be written over
               sig3(i,j) = 0
               sig5(i,j) = 0
            enddo
            do i = 1,3
               sig1(i,j) = wsig1(i,j)
                               ! j = 1 for graphs 1, 2.  j = 2 for
               sig3(i,j) = wsig3(i,j)
                               ! graphs 3, 4.  j = 3 for graphs 5, 6.
               sig5(i,j) = wsig5(i,j)
                               ! j = 4 for graph 7.  (same as Z case)
            enddo			
         enddo
         
      else
         print*," Unacceptable BOS value in QQBQQI. BOS = ",bos
         print*," Recompile with BOS = 1,2,3, or 4 "
         stop
      endif

C     next set up the range of boson identifiers for the lower index
c     (VLMIN/MAX) and the higher index (VHMIN/MAX) exchanged boson
c     for the 7 Feynman graphs
      do j = 1, nprmax
         do i = 1, 2            ! graphs 1, 2
            vlmin(j,i) = b1min(bos) ! V1/bos is produced boson 
            vlmax(j,i) = b1max(bos) !	denoted b1
            vhmin(j,i) = b2min(j,bos) ! no triple vertex so 
            vhmax(j,i) = b2max(j,bos) !     	b2/V2 = b3/V3
         enddo
         
         do i = 3, 4            ! graphs 3, 4
            vlmin(j,i) = b1min(bos) ! V1/bos is produced boson
            vlmax(j,i) = b1max(bos) ! 	denoted b1
            vhmin(j,i) = b3min(j,bos) ! no triple vertex so 
            vhmax(j,i) = b3max(j,bos) ! 	b2/V2 = b3/V3
         enddo
         
         do i = 5, 7            ! graphs 5, 6.  No V1 produced
            vlmin(j,i) = b2min(j,bos) ! Vl = V2/b2, Vh = V3/b3.  In 
            vlmax(j,i) = b2max(j,bos) ! general b2 .ne. b3 .
            vhmin(j,i) = b3min(j,bos) ! graph 7 ==> 3 boson vertex,
            vhmax(j,i) = b3max(j,bos) ! use all 3 bosons, V1, VL, VH
         enddo
         
C     No 3 boson vertex for processes 1,...,4 in Z production.  Z
c     production ==> a W from both the 1-2 and 3-4 fermion lines.  
         if (bos.eq.2 .and. j.le.4) then
            vlmin(j,7) = 0      ! normal range is 1,...,6
            vlmax(j,7) = -1     ! 0, -1 correspond to 
            vhmin(j,7) = 0      !	nothing
            vhmax(j,7) = -1
         endif
      enddo
      
C  Reset the coupling factors
      do j = 1,7             ! all graphs
         do i = 1,6          ! all flavors
            do isig = 1,8    ! all helicity combinations
               do vh = 1,4
                  do vl = 1,4
                     fcpl(vl,vh,isig,i,j) = 0d0
                  enddo
               enddo
            enddo
         enddo
      enddo

c     This huge loop sets up the products of couplings for graphs
c     1,...,6.  The strong coupling for the gluon vertex will be
c     added later in the code.
      do j = 1,6 ! The first 6 graphs
         do isig = 1,sigmax ! all helicity combinations
            isig1 = sig1(isig,(j+1)/2)         ! (j+1)/2 = integer ==> same
            isig3 = sig3(isig,(j+1)/2)         ! values for graphs 1 & 2,
            isig5 = sig5(isig,(j+1)/2)         ! 3 & 4, 5 & 6.
            do i = 1,nprmax                  ! all flavor combinations
               do vh = vhmin(i,j),vhmax(i,j)
                  do vl = vlmin(i,j),vlmax(i,j)
                     if (j.eq.1) then
                        rejgr = (bos.eq.3 .and. i.eq.3) .or.
     $                          (bos.eq.4 .and. i.eq.4)
                        if (.not.rejgr)
     $                       fcpl(vl,vh,isig,i,1) =
     $                       clr(fl1(i,bos),vh,isig1)
     $                       *clr(fl3(i,bos),vh,isig3)
     $                       *clr(fl4(i,bos),vl,isig3)*
     $                       clr(lflavr(5),vl,isig5)
                        if(bos.eq.1)
     $                       fcpl(vl,vh,isig,i,1) =
     $                       clr(fl1(i,bos),vh,isig1)
     $                       *clr(fl3(i,bos),vh,isig3)
     $                       *clr(fl4(i,bos),vl,isig3)
c insert non-executed print of fcpl to avoid g77 assigning 0 to fcpl
                        if (vl.eq.-2) then
                           print*," coupling for graph 1 ",
     $                            fcpl(vl,vh,isig,i,1)
                           print*," lwarn = ",lwarn
                        endif
                     elseif (j.eq.2) then
                        rejgr = (bos.eq.3 .and. i.eq.4) .or.
     $                          (bos.eq.4 .and. i.eq.3)
                        if (.not.rejgr)
     $                       fcpl(vl,vh,isig,i,2) =
     $                       clr(fl1(i,bos),vh,isig1)
     $                       *clr(fl4(i,bos),vh,isig3)
     $                       *clr(fl3(i,bos),vl,isig3)
     $                       *clr(lflavr(5),vl,isig5)
                        if (bos.eq.1)
     $                       fcpl(vl,vh,isig,i,2) =
     $                       clr(fl1(i,bos),vh,isig1)
     $                       *clr(fl4(i,bos),vh,isig3)
     $                       *clr(fl3(i,bos),vl,isig3)
                     elseif (j.eq.3) then
                        rejgr = (bos.eq.3 .and. i.eq.1) .or.
     $                          (bos.eq.4 .and. i.eq.2)
                        if (.not.rejgr)
     $                       fcpl(vl,vh,isig,i,3) =
     $                       clr(fl1(i,bos),vh,isig1)
     $                       *clr(fl4(i,bos),vh,isig3)*
     $                       clr(fl2(i,bos),vl,isig1)
     $                       *clr(lflavr(5),vl,isig5)
                        if (bos.eq.1)
     $                       fcpl(vl,vh,isig,i,3) =
     $                       clr(fl1(i,bos),vh,isig1)
     $                       *clr(fl4(i,bos),vh,isig3)*
     $                       clr(fl2(i,bos),vl,isig1)
                     elseif (j.eq.4) then
                        rejgr = (bos.eq.3 .and. i.eq.2) .or.
     $                          (bos.eq.4 .and. i.eq.1)
                        if (.not.rejgr)
     $                       fcpl(vl,vh,isig,i,4) =
     $                       clr(fl2(i,bos),vh,isig1)
     $                       *clr(fl4(i,bos),vh,isig3)
     $                       *clr(fl1(i,bos),vl,isig1)
     $                       *clr(lflavr(5),vl,isig5)
                        if (bos.eq.1)
     $                       fcpl(vl,vh,isig,i,4) =
     $                       clr(fl2(i,bos),vh,isig1)
     $                       *clr(fl4(i,bos),vh,isig3)
     $                       *clr(fl1(i,bos),vl,isig1)
                     elseif (j.eq.5) then
                        rejgr = bos.eq.2 .and. i.ge.5 .and.
     $                          mod(lflavr(5),2).eq.mod(i,2)
                        if (.not.rejgr)
     $                       fcpl(vl,vh,isig,i,5) =
     $                       clr(fl1(i,bos),vl,isig1)
     $                       *clr(lflavr(5),vl,isig5)
     $                       *clr(lflavr(6),vh,isig5)
     $                       *clr(fl3(i,bos),vh,isig3)
                        if(bos.eq.1)
     $                       fcpl(vl,vh,isig,i,5) = 0d0
                     elseif (j.eq.6) then
                        rejgr = bos.eq.2 .and. i.ge.5 .and.
     $                          mod(lflavr(5),2).ne.mod(i,2)
                        if (.not.rejgr)
     $                       fcpl(vl,vh,isig,i,6) =
     $                       clr(fl1(i,bos),vl,isig1)
     $                       *clr(lflavr(6),vl,isig5)
     $                       *clr(lflavr(5),vh,isig5)
     $                       *clr(fl3(i,bos),vh,isig3)
                        if (bos.eq.1)
     $                       fcpl(vl,vh,isig,i,6) = 0d0

                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo
c      stop
c     assign the products of the couplings for the triple boson
c     vertex/graph 7.
      if (bos.eq.2.or.bos.eq.1) then
         imin = 5
      else
         imin = 1
      endif
      do isig = 1,3
         isig1 = sig1(isig,4)
         isig3 = sig3(isig,4)
         isig5 = sig5(isig,4)
         do i = imin,nprmax
            do iv = 1,2
               if (vlmin(i,7).eq.vlmax(i,7)) then
                  v2 = vlmin(i,7)
                  if (vhmin(i,7).eq. vhmax(i,7)) then
                     v3 = vhmin(i,7)
                     v1 = iv
                  else
                     v3 = iv
                     v1 = 7-v2
                  endif
               else
                  v2 = iv
                  v3 = vhmin(i,7)
                  v1 = 7-v3
               endif
               rejgr = (bos.eq.2 .and. isig.eq.3).or.(bos.eq.1 .and. isig.eq.3)
               if (.not.rejgr) fcpl(iv,1,isig,i,7) =
     $              clr(fl1(i,bos),v2,isig1)
     $              *clr(fl3(i,bos),v3,isig3)
     $              *clr(lflavr(5),v1,isig5)
     $              *b(iv,3,4)
               if (.not.rejgr.and.bos.eq.1) fcpl(iv,1,isig,i,7) =
     $              clr(fl1(i,bos),v2,isig1)
     $              *clr(fl3(i,bos),v3,isig3)
!     $              *clr(lflavr(5),v1,isig5)
     $              *b(iv,3,4)
            enddo
         enddo
      enddo
      call tbvini(.true.)
      if (ldebug) then
         open (12,file="zqqv",status="unknown")
         write (12,*) fcpl
         close(12)
      endif
      if (.not.ldoblha) then
        print*," "
        print*," Setup of couplings in qqbqqi successfully",
     &         " completed "
        print*," W/Z production with finite width implemented via"
        print*," "
        print*,"            igauge = ",igauge
        print*," "
        if (igauge.eq.1) then
           print*," i.e. the overall factor scheme is used "
        elseif (igauge.eq.2) then
           print*," i.e. the complex mass scheme is used "
        elseif (igauge.eq.3 .and. bos.eq.2) then
           print*," i.e. the naive Z Breit-Wigner scheme is used "
        elseif (igauge.eq.-1) then
           print*," i.e. the W/Z width is set to 0 for testing purposes "
        else
           print*," which is illegal. Reset it in global.inc "
           stop
        endif
! moved to cuts.F
!      if (qsqAmin.gt.0d0) then
!         print*," "
!         print 171," Photon exchange with virtuality Q^2 < ",qsqAmin,
!     1   " GeV^2 is eliminated. "
!         print*," Recover by calculating  q gamma ---> q g V   etc."
!      endif
! 171  format(a,f10.1,a)
        print*," "
        print*,"---------------------------------------------------------"
      endif
      end


C
C********************  qqZqqMod ***********************************
C
      subroutine qqZqq(pbar,sign,NLO, 
     1                 uucc,uuss,ddcc,ddss,udsc,ducs)
          use globalvars, only: ldoblha
      implicit none
      double precision uucc,uuss,ddcc,ddss,
     1       udsc,ducs,          pbar(0:3,6),
     2       uucs,ddcs,udcc,udss,
     3       uusc,ddsc,ducc,duss
      integer sign(6), NLO
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"
C  This code is modified to allow for virtual corrections, more precisely
C  the interference of Born with the finite part of virtual diagrams
C  for 
c
c  INPUT:  NLO = 1       return uucc = |M_born|^2 + 2Re(M_Born^* M_virt)
c          NLO = 0       return uucc = |M_born|^2   etc.
c
c
c  There are two options for selecting helicties: 1) do loop summation
c  to provide unpolarized cross sections, and 2) calculation of 1 chirality
c  combination only. See common block /chelsum/ below. uucc etc. are the 
c  1) polarization summed cross section or 2) the result for the selected
c  helicity combination
c
c	Dieter Zeppenfeld, <dieter@pheno.physics.wisc.edu>
c	Initial version:  2002 October 24
c	Last modified:    2003 January 3 
C
C  QQBQQ calculates the matrix elements**2 for electroweak single
c  weak boson production by quark quark scattering
C
C        q1 q3    ---->   q2 q4 V,   V ---> f5-bar f6
C
C  and crossing related processes. Pauli interference terms for
c  identical fermions are neglected. In particular, only the
c  t-channel exchange of elctroweak bosons is considered. s-channel
c  production of weak boson pairs is NOT implemented.
C
C  QQBQQI must be called first to initialize the desired production
c  process. The main task of QQBQQI is to set up the products of coupling
c  constants needed in Feynman graphs a, ..., g which are stored in
c  fcpl
      double precision fcpl(4,4,8,6,7)
      integer vlmin(6,7), vlmax(6,7), vhmin(6,7), vhmax(6,7), 
     1     nprmax, sigmax, sig1(8,4), sig3(8,4), sig5(8,4)
      common /cqqbqq/ fcpl, nprmax, sigmax,
     1     vlmin, vlmax, vhmin, vhmax, sig1, sig3, sig5
      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)
C  Here fcpl(vl,vh,sig,i,j) contains the coupling constant factors. 
C     vl  = 1,..,4    the lower one of the boson indices i of the V_i
C     vh  = 1,..,4    the higher one of the boson indices i of the V_i
C                     vl = 2 e.g. means V_lower is a Z-boson
C     sig = 1,..,8    identifies the (maximally) 8 different helicity
C                     combinations
C     i   = 1,..,6    numbers the possible subprocesses
c                     (4 for W production and 6 for Z production)
C     j   = 1,..,7    the 7 Feynman graphs
C
c  helicity selection
      integer jsig, jsig1, jsig3, jsig5
      common /chelsum/ jsig,jsig1,jsig3,jsig5
C Select helicity combination by setting jsig externally. This program will
c then fill jsig1,jsig3,jsig5 to correspond to the chiralities of the 1-2,
c 3-4 and 5-6 fermion lines for external use, e.g. in calculating polarized
c cross sections. This is option 2 and it requires that jsig lies within
c the allowed sig range [1,sigmax] for the process. If jsig is outside this
c range the sum, NOT average!, over all polarization states will be calculated
c
c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"

C  local variables
      integer i,ibos,j,jj,k,mu,vl,vh,ic,
     1        isig,isig1,isig3,isig5, jsigmin, jsigmax
      double precision res(6,-1:1),resv(6,-1:1), p(0:3,6),
     1       p21(0:4),p43(0:4),p65(0:4), dummy(0:4)
      double precision tree(6)
      integer id(4)
      double precision pi2o3,c2,c2o4pi
      parameter (pi2o3=pi**2/3d0)
      parameter (c2=4d0/3d0, c2o4pi=c2/4d0/pi)
      logical lnlo
      double precision lnmusqp1p2, lnmusqp3p4

#include "VBFNLO/utilities/koppln.inc"

      double complex psi(2,-1:1,6), j21(0:3,-1:1), j43(0:3,-1:1),
     1           j65(0:3,-1:1), prop65(4), propl(4,3), proph(4,3), 
     2           corr56, braket(2), matr, matr7(2), 
     3           fac(6,8,7),mm(6,8), s1c, mvirt, fac_virt(6,8,4),
     4           mv12(6,8), mv34(6,8)
      save fac,fac_virt
      double complex prop21(4), prop43(4)
      external s1c
      double complex m1, je1(0:5),je2(0:5)
      DOUBLE COMPLEX resultgauge(2),resultb,resultgaugeb(2)
      integer jd, isPhoton, isNoPhoton
      logical ldebug, lflagv
      parameter (ldebug=.false.)
      common /bnumerr/lflagv           ! if true: numerical problems in

c      double precision td


! calculation of M_virt
      double precision scale(2),pcall(0:3,6),cte!, resultgauge(2), resultgaugeb(2)
      integer igaugetest,icomp
      parameter(igaugetest=3, cte=-pi**2/3d0+7d0)
      double precision dotrr
      external dotrr
      lnlo = NLO.eq.1
      lflagv = .false.

C  Entry for Z production
      ibos = 2
      goto 40

C  Entry for W+ production
      entry qqWPqq(pbar,sign,NLO, uucs,ddcs,udcc,udss)
      ibos = 3
      goto 40

C  Entry for W- production
      entry qqWMqq(pbar,sign,NLO, uusc,ddsc,ducc,duss)
      ibos = 4
      goto 40

C  Entry for A production
      entry qqAqq(pbar,sign,NLO, uucc,uuss,ddcc,ddss,udsc,ducs)
      ibos = 1

40    continue  

      lnlo = NLO.eq.1
      lflagv = .false.

c assign chirality summation range
      if (jsig.ge.1 .and. jsig.le.sigmax) then !fixed chirality selection
         jsigmin = jsig
         jsigmax = jsig
         jsig1 = sig1(jsig,3)
         jsig3 = sig3(jsig,3)
         jsig5 = sig5(jsig,3)
      else                                     !do loop over chiralities
         jsigmin = 1
         jsigmax = sigmax
         jsig1 = 0
         jsig3 = 0
         jsig5 = 0
      endif
c reset summation arrays fac and fac_virt to zero
      do j = 1,7                  ! all graphs
         do k = 1,nprmax          ! all flavors
            do isig = 1,sigmax    ! all helicity combinations
               fac(k,isig,j) = (0d0,0d0)
               if (j.le.4) fac_virt(k,isig,j) = (0d0,0d0)
            enddo
         enddo
      enddo
C  Define the internal momenta
      do mu = 0,3
         do i = 1,6
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo
         p21(mu) = p(mu,2)-p(mu,1)
         p43(mu) = p(mu,4)-p(mu,3)
         p65(mu) = p(mu,6)-p(mu,5)
      enddo
      if (lwarn) then
         do mu=0,3
            dummy(mu) = p(mu,1)+p(mu,3)+p(mu,5)-
     1                 (p(mu,2)+p(mu,4)+p(mu,6))
         enddo
         if (abs(dummy(0)).gt.1d-8 .or.
     1       abs(dummy(1)).gt.1d-8 .or.
     2       abs(dummy(2)).gt.1d-8 .or.
     3       abs(dummy(3)).gt.1d-8 ) then
            print*," momentum mismatch in qqbqq "
            print 102," dummy = ",(dummy(mu),mu=0,3)
            print 103," sign  = ",sign
 101        format( " p(", i1, ") = ", 4(f10.3, 2x) )
 102        format(a,4f10.3)
 103        format(a,6i5)
            do i = 1,6
               write(6,101) i, p(0,i), p(1,i), p(2,i), p(3,i)
            end do
         endif
      endif

      p21(4) = p21(0)**2-p21(1)**2-p21(2)**2-p21(3)**2
      p43(4) = p43(0)**2-p43(1)**2-p43(2)**2-p43(3)**2
      p65(4) = p65(0)**2-p65(1)**2-p65(2)**2-p65(3)**2

      
      scale(1)=abs(2d0*dotrr(p(0,1),p(0,2)))
      scale(2)=abs(2d0*dotrr(p(0,3),p(0,4)))

c igauge allows to select different ways to insure EW gauge invariance
c in the face of finite W and Z width effects
c  
c   igauge = 1       overall factor scheme
c            2       complex mass scheme, m**2-->m**2-i*m*Gamma
c                    in propagators
c            3       for ibos=2 only: naive Z Breit Wigner
c           -1       width = 0, for testing purposes only
c
c igauge is set in global.inc

C  Get the vector boson propagator factors
      prop21(1) = 1/p21(4)
      if (abs(igauge).eq.1 .or. igauge.eq.3) then
         prop21(2) = 1/(p21(4)-xm2(2))
         prop21(3) = 1/(p21(4)-xm2(3))
      elseif (igauge.eq.2) then
         prop21(2) = 1/dcmplx(p21(4)-xm2(2),xmg(2))
         prop21(3) = 1/dcmplx(p21(4)-xm2(3),xmg(3))
      endif
      prop21(4) = prop21(3)

      prop43(1) = 1/p43(4)
      if (abs(igauge).eq.1 .or. igauge.eq.3) then
         prop43(2) = 1/(p43(4)-xm2(2))
         prop43(3) = 1/(p43(4)-xm2(3))
      elseif (igauge.eq.2) then
         prop43(2) = 1/dcmplx(p43(4)-xm2(2),xmg(2))
         prop43(3) = 1/dcmplx(p43(4)-xm2(3),xmg(3))
      endif
      prop43(4) = prop43(3)

      if (ibos.eq.2) then
         prop65(1) = 1/p65(4)
         if (igauge.ne.-1) then
            prop65(2) = 1/dcmplx(p65(4)-xm2(2),xmg(2))
         else
            prop65(2) = 1/(p65(4)-xm2(2))
         endif
         prop65(3) = 0
         prop65(4) = 0
C  correction factor for diagrams 5,6 and s-channel photon 
         if (igauge.eq.1) then
            corr56 = (p65(4)-xm2(2))*prop65(2)
            prop65(1) = prop65(1) * corr56
         else
            corr56 = 1
         endif
      elseif (ibos.eq.1) then
         prop65(1) = 1d0
         prop65(2) = 0d0
         prop65(3) = 0d0
         prop65(4) = 0d0
C  correction factor for diagrams 5,6 and s-channel photon 
         !diagrams 5 and 6 are not existent, so set
         corr56=0d0
      else
         if (igauge.ne.-1) then
            prop65(3) = 1/dcmplx(p65(4)-xm2(3),xmg(3))
         else
            prop65(3) = 1/(p65(4)-xm2(3))
         endif
         prop65(4) = prop65(3)
         prop65(1) = 0
         prop65(2) = 0
C  correction factor for diagrams 5,6
         if (igauge.eq.1) then
            corr56 = (p65(4)-xm2(3))*prop65(3)
         else
            corr56 = 1
         endif
      endif
C
C  Assign the possible propagators for the three classes of Feynman
c  graphs (excluding the 3 boson vertex)
      do vh = 1,4
         propl(vh,1) = prop65(vh)
         proph(vh,1) = prop21(vh)
         propl(vh,2) = prop65(vh)
         proph(vh,2) = prop43(vh)
         propl(vh,3) = prop21(vh)
         proph(vh,3) = prop43(vh)
      enddo
C
C  Get the external spinors (including factor sqrt(2E) )
      if (ibos.eq.1) then
        call psi0m(4,pbar,sign,psi)
      else
        call psi0m(6,pbar,sign,psi)
      endif

C  get the f-fbar currents J21^mu, J43^mu, J65^mu
      call curr(1,psi(1,-1,2),psi(1,-1,1),j21)
      call curr(1,psi(1,-1,4),psi(1,-1,3),j43)
      if (ibos.eq.2) then
         call curr(1,psi(1,-1,6),psi(1,-1,5),j65)
      elseif (ibos.eq.1) then
         call helvec(p65,1,-1,j65(0,-1)) ! helicity=-1
         j65(0,0)=0d0
         j65(1,0)=0d0
         j65(2,0)=0d0
         j65(3,0)=0d0
         call helvec(p65,1, 1,j65(0, 1)) ! helicity=+1
      else
         call curr(-1,psi(1,-1,6),psi(1,-1,5),j65)
      endif

c  fix ID for boxline call for first 4 Feynman graphs
      id(1) = 6-sign(3)         ! ID = 5,7 for graph 1
      id(2) = 7-sign(3)         ! ID = 6,8 for graph 2
      id(3) = 2-sign(1)         ! ID = 1,3 for graph 3
      id(4) = 3-sign(1)         ! ID = 2,4 for graph 4

      if(lnlo.and.ibos.eq.1)then
         do mu=0,3
            pcall(mu,1)=-p(mu,1)
            pcall(mu,2)=p(mu,2)
            pcall(mu,4)=p(mu,4)
            pcall(mu,3)=-p(mu,3)
         enddo
      endif

      if (ibos.eq.1) then
         isPhoton = 1
         isNoPhoton = 0
      else
         isPhoton = 0
         isNoPhoton = 1
      endif

C  get the amplitudes of the first 6 Feynman graphs for all allowed
c  helicity combinations
      do 50 j = 1,6             ! Feynman graphs
         jj = (j+1)/2
         icomp = 1
         do 50 isig =(1-nlo)*jsigmin*isNoPhoton+isPhoton,jsigmax,
     $        (1-nlo)+nlo*jsigmin
                                          ! helicity combinations
            if (isig.eq.0) then
               if (j.ge.5) goto 50
               isig1 = sig1(1,3)
               isig3 = sig3(1,3)
               isig5 = 0
               do ic = 0,3
                  j65(ic,isig5) = dcmplx(p65(ic),0d0)
               enddo
            else
               isig1 = sig1(isig,jj)
               if (isig1.eq.0) goto 50
               isig3 = sig3(isig,jj)
               isig5 = sig5(isig,jj)
            endif
C
C  graphs a and b (1,2): 2 bosons attached to 34 line
            if (j.eq.1) then
               call bra2c(psi(1,isig3,4),.true.,p(0,4),isig3,
     &              p65,j65(0,isig5),braket,dummy)


               matr = -s1c(braket,j21(0,isig1),.true.,isig3,
     &              psi(1,isig3,3))


#ifdef WITH_NLO 
               if (lnlo) then
                  if (ibos.ne.1) then
                     call boxline(id(j), 
     1                    psi(1,isig3,3),psi(1,isig3,4),p(0,3),p(0,4),isig3,
     2                    j21(0,isig1),j65(0,isig5),p21,p65,
     3                    matr,   mvirt)
                  else
                     CALL BOXLINEMM(pcall(0,3),p21,p65,pcall(0,4),
     $                    psi(1,isig3,4),psi(1,isig3,3),
     $                    j21(0,isig1),j65(0,isig5),isig3,scale(2),
     $                    igaugetest,icomp,resultgauge,mvirt,
     $                    resultgaugeb,resultb)
                     icomp=-1
                  endif
               endif
                  
cdebug
               if (ldebug.and.ibos.ne.1) then
                  do mu = 0,3
                     je1(mu) = j21(mu,isig1)
                     je2(mu) = j65(mu,isig5)
                  enddo
                  je1(4) = dcmplx(p21(0),p21(3))
                  je1(5) = dcmplx(p21(1),p21(2))
                  je2(4) = dcmplx(p65(0),p65(3))
                  je2(5) = dcmplx(p65(1),p65(2))
                  call boxline6(jd, 
     1              psi(1,isig3,3),psi(1,isig3,4),p(0,3),p(0,4),isig3,
     1              je1, je2,
     3              matr,   m1)
                  if (isig.eq.0) then
                     print*," mborn = ",matr,mvirt,m1
                  else
                     print*," j=1, isig3 = ",isig3,":",m1/mvirt,mvirt
                  endif
!                  read(*,"(a)") text
               endif
#endif
cdebug
            elseif (j.eq.2) then
               call ket2c(psi(1,isig3,3),.true.,p(0,3),isig3,
     &              p65,j65(0,isig5),braket,dummy)
               matr = -s1c(psi(1,isig3,4),j21(0,isig1),.true.,
     &              isig3,braket)
#ifdef WITH_NLO 
               if (lnlo) then
                  if (ibos.ne.1) then 
                     call boxline(id(j), 
     1                    psi(1,isig3,3),psi(1,isig3,4),p(0,3),p(0,4),isig3,
     2                    j65(0,isig5),j21(0,isig1),p65,p21,
     3                    matr,   mvirt)
                  else
                     CALL BOXLINEMM(pcall(0,3),p65,p21,pcall(0,4),
     $                    psi(1,isig3,4),psi(1,isig3,3),
     $                    j65(0,isig5),j21(0,isig1),isig3,scale(2),
     $                    igaugetest,icomp,resultgauge,mvirt,
     $                    resultgaugeb,resultb)
                     icomp=-1
                  endif
               endif
#endif
C
C  graphs c and d (3,4)
            elseif (j.eq.3) then
               call bra2c(psi(1,isig1,2),.true.,p(0,2),isig1,
     &              p65,j65(0,isig5),braket,dummy)
               matr = -s1c(braket,j43(0,isig3),.true.,isig1,
     &              psi(1,isig1,1))
#ifdef WITH_NLO 
               if (lnlo) then
                  if (ibos.ne.1) then 
                     call boxline(id(j), 
     1                    psi(1,isig1,1),psi(1,isig1,2),p(0,1),p(0,2),isig1,
     2                    j43(0,isig3),j65(0,isig5),p43,p65,
     3                    matr,   mvirt)
                  else
                     CALL BOXLINEMM(pcall(0,1),p43,p65,pcall(0,2),
     $                    psi(1,isig1,2),psi(1,isig1,1),
     $                    j43(0,isig3),j65(0,isig5),isig1,scale(1),
     $                    igaugetest,icomp,resultgauge,mvirt,
     $                    resultgaugeb,resultb)
                     icomp=-1
                  endif
               endif
#endif

            elseif (j.eq.4) then
               call ket2c(psi(1,isig1,1),.true.,p(0,1),isig1,
     &              p65,j65(0,isig5),braket,dummy)
               matr = -s1c(psi(1,isig1,2),j43(0,isig3),.true.,
     &              isig1,braket)
#ifdef WITH_NLO 
               if (lnlo) then
                  if (ibos.ne.1) then 
                     call boxline(id(j), 
     1                    psi(1,isig1,1),psi(1,isig1,2),p(0,1),p(0,2),isig1,
     2                    j65(0,isig5),j43(0,isig3),p65,p43,
     3                    matr,   mvirt)
                  else
                     CALL BOXLINEMM(pcall(0,1),p65,p43,pcall(0,2),
     $                    psi(1,isig1,2),psi(1,isig1,1),
     $                    j65(0,isig5),j43(0,isig3),isig1,scale(1),
     $                    igaugetest,icomp,resultgauge,mvirt,
     $                    resultgaugeb,resultb)
                     icomp=-1
                  endif
               endif
#endif

C
C     graphs e and f (5,6): final state leptons: no NLO QCD corrections
            elseif (j.eq.5.and.ibos.ne.1) then
               call bra2c(psi(1,isig5,6),.true.,p(0,6),isig5,
     &              p43,j43(0,isig3),braket,dummy)
               matr = -s1c(braket,j21(0,isig1),.true.,isig5,
     &              psi(1,isig5,5))
               matr = matr*corr56

            elseif (j.eq.6.and.ibos.ne.1) then
               call ket2c(psi(1,isig5,5),.true.,p(0,5),isig5,
     &              p43,j43(0,isig3),
     1              braket,dummy)
               matr = -s1c(psi(1,isig5,6),j21(0,isig1),.true.,
     &              isig5,braket)
               matr = matr*corr56

            endif
c
c for isig=0 and nlo=1 check that the virtual boxline contribution is small
c compared to the Born amplitude
            if (isig.eq.0 .and.ibos.ne.1) then
               if (abs(mvirt).gt.0.1*abs(matr)) then
c set flag at 1d-1. This should be compared to pi^2/3-7 = 3.71 times M_born
c to be added below. Both factors will be multiplied by alphas/3pi at the end!
c out of a run of 10^6 events this gave 80 correct flags and 210 wrong ones 
c based on recalculating the |M|^2 at a boosted phase space point
                  lflagv = .true.
                  if (ldebug) then
                     print*," WARNING: numerical problems for event "
                     print*," Feynman graph ",j
                     print*," mvirt = ",mvirt
                     print*," mborn = ",matr
                     print*," |ratio| = ",abs(mvirt/matr)
                  endif
               endif
               goto 50
            endif
C
C  now get the coupling*propagator factors for subprocess k, helicity
C  combination isig, Feynman graph j
            do k = 1,nprmax
               fac(k,isig,j) = 0
               do vh = vhmin(k,j),vhmax(k,j)
                  do vl = vlmin(k,j),vlmax(k,j)
                     fac(k,isig,j) = fac(k,isig,j) + 
     1               fcpl(vl,vh,isig,k,j)*propl(vl,jj)*proph(vh,jj)
                  enddo
               enddo

               if (ibos.ne.1) then
                  if (j.le.4) then
                     if (.not.lflagv .and. (mvirt.eq.mvirt)) then
                        fac_virt(k,isig,j) = fac(k,isig,j)*mvirt
                     else
                        fac_virt(k,isig,j) = 0
                     endif
                  endif
               else
                  if (j.le.4) then
                     if (.not.lflagv .and. (mvirt.eq.mvirt)) then
                        fac_virt(k,isig,j) = fac(k,isig,j)*(mvirt+cte*matr)
                     else
                        fac_virt(k,isig,j) = 0
                     endif
                  endif
               endif
               fac(k,isig,j) = fac(k,isig,j)*matr
               
               if(j.gt.4.and.ibos.eq.1) fac(k,isig,j) = 0d0 
            enddo
50    continue
C
C  and finally the three boson vertex. An extra factor of -1 is
c  included because OUTGOING momenta p21,p43,p65 are used in the call
c  of tbv The correspondence of the boson identifiers V_i and the
c  currents J_ij is
C          V_1  =  J_65
C          V_2  =  J_21
C          V_3  =  J_43
C 
      do 60 isig = 1,8
         do 60 k = 1,nprmax
            fac(k,isig,7) = 0
60    continue
        if (ibos.eq.2) then
         do isig = max(jsigmin,1),min(jsigmax,2),1
            isig1 = sig1(isig,4)
            isig3 = sig3(isig,4)
            isig5 = sig5(isig,4)
C process 5 has V2 = W+, V3 = W-
            call tbv(j21(0,isig1),j43(0,isig3),j65(0,isig5),
     &           p21,p43,p65,matr7)
            fac(5,isig,7) = matr7(1)*prop65(1)*fcpl(1,1,isig,5,7) +
     1                      matr7(2)*prop65(2)*fcpl(2,1,isig,5,7)
C process 6 has V2 = W-, V3 = W+
            call tbv(j43(0,isig3),j21(0,isig1),j65(0,isig5),p43,
     &           p21,p65,matr7)
            fac(6,isig,7) = matr7(1)*prop65(1)*fcpl(1,1,isig,6,7) +
     1                      matr7(2)*prop65(2)*fcpl(2,1,isig,6,7)
            do k = 5,6
               fac(k,isig,7) = fac(k,isig,7)*prop21(3)*prop43(4)
            enddo
         enddo
      elseif (ibos.eq.3) then
         do 70 k = 1,4
            do 70 isig = jsigmin,jsigmax
               isig1 = sig1(isig,4)
               isig3 = sig3(isig,4)
               isig5 = sig5(isig,4)
               if (k.le.2 .and. isig.ne.3) then
C processes 1 and 2 have V1 = W-, V2 = Z/gamma, V3 = W+
                  call tbv(j43(0,isig3),j65(0,isig5),j21(0,isig1),
     &                     p43,p65,p21,matr7)
                  matr7(1) = matr7(1)*prop21(1)*prop43(3)
                  matr7(2) = matr7(2)*prop21(2)*prop43(3)
               elseif (k.gt.2 .and. isig.ne.2) then
C processes 3 and 4 have V1 = W-, V2 = W+, V3 = Z/gamma
                  call tbv(j21(0,isig1),j65(0,isig5),j43(0,isig3),
     &                     p21,p65,p43,matr7)
                  matr7(1) = matr7(1)*prop43(1)*prop21(3)
                  matr7(2) = matr7(2)*prop43(2)*prop21(3)
               else
                  goto 70
               endif
               fac(k,isig,7) = ( matr7(1)*fcpl(1,1,isig,k,7)+
     1                           matr7(2)*fcpl(2,1,isig,k,7)  )
     &                       *prop65(4)
70       continue 
      elseif (ibos.eq.4) then
         do 80 k = 1,4
            do 80 isig = jsigmin,jsigmax
               isig1 = sig1(isig,4)
               isig3 = sig3(isig,4)
               isig5 = sig5(isig,4)
               if (k.le.2 .and. isig.ne.3) then
C processes 1 and 2 have V1 = W+, V2 = Z/gamma, V3 = W-
                  call tbv(j65(0,isig5),j43(0,isig3),j21(0,isig1),
     &                     p65,p43,p21,matr7)
                  matr7(1) = matr7(1)*prop21(1)*prop43(4)
                  matr7(2) = matr7(2)*prop21(2)*prop43(4)
               elseif (k.gt.2 .and. isig.ne.2) then
C processes 3 and 4 have V1 = W+, V2 = W-, V3 = Z/gamma
                  call tbv(j65(0,isig5),j21(0,isig1),j43(0,isig3),
     &                     p65,p21,p43,matr7)
                  matr7(1) = matr7(1)*prop43(1)*prop21(4)
                  matr7(2) = matr7(2)*prop43(2)*prop21(4)
               else
                  goto 80
               endif
               fac(k,isig,7) = ( matr7(1)*fcpl(1,1,isig,k,7)+
     1                           matr7(2)*fcpl(2,1,isig,k,7)  )
     &                       *prop65(3)
80       continue 
      elseif (ibos.eq.1) then    
         do isig = max(jsigmin,1),min(jsigmax,2),1
            isig1 = sig1(isig,4)
            isig3 = sig3(isig,4)
            isig5 = sig5(isig,4)
C     process 5 has V2 = W+, V3 = W-
            call tbv(j21(0,isig1),j43(0,isig3),j65(0,isig5),
     &           p21,p43,p65,matr7)
            fac(5,isig,7) = matr7(1)*prop65(1)*fcpl(1,1,isig,5,7) +
     1           matr7(2)*prop65(2)*fcpl(2,1,isig,5,7)
C     process 6 has V2 = W-, V3 = W+
            call tbv(j43(0,isig3),j21(0,isig1),j65(0,isig5),p43,
     &           p21,p65,matr7)
            fac(6,isig,7) = matr7(1)*prop65(1)*fcpl(1,1,isig,6,7) +
     1           matr7(2)*prop65(2)*fcpl(2,1,isig,6,7)
            do k = 5,6
               fac(k,isig,7) = fac(k,isig,7)*prop21(3)*prop43(4)
            enddo
         enddo
      endif
C
C  Now sum the contributions from the seven Feynman graphs for all
c  subprocesses
c      if (ldebug) print*," sigmax = ",sigmax
      do k = 1,nprmax
         do isig5 = -1,1
            res(k,isig5) = 0
            resv(k,isig5) = 0
            tree(k) = 0
         enddo
         do isig = jsigmin,jsigmax
            isig5 = sig5(isig,3)
            if (ldebug) then
               if( (sig5(isig,1).ne.isig5 .and. sig5(isig,1).ne.0) .or.
     1             (sig5(isig,2).ne.isig5 .and. sig5(isig,2).ne.0) .or.
     2             (sig5(isig,4).ne.isig5 .and. sig5(isig,4).ne.0) 
     3           ) then
                  print*," check isig5 assignment for isig = ",isig
                  print*," sig5 = ",(sig5(isig,ic),ic=1,4)
               endif
            endif
            mm(k,isig) = 0
            do j = 1,7
               mm(k,isig) = mm(k,isig) + fac(k,isig,j)
            enddo
            res(k,isig5) = res(k,isig5) + dreal(mm(k,isig))**2
     &                                  + dimag(mm(k,isig))**2
            if (lnlo) then
               mv12(k,isig) = (0d0,0d0)         ! boxline correction to 12 line
               mv34(k,isig) = (0d0,0d0)         ! boxline correction to 34 line
               do j = 1,2
                  mv34(k,isig) = mv34(k,isig) + fac_virt(k,isig,j)
                  mv12(k,isig) = mv12(k,isig) + fac_virt(k,isig,j+2)
               enddo
c  add Born type term and multiply by F_q = alphas*C_2/4pi
c  the factor pi^2/3+9/2 for the born term is after adding the subtraction term
c  and the counter term for the renormalization of the pdfs
               mv12(k,isig) = als(1,1)*c2o4pi*( mv12(k,isig) + 
     &                                          mm(k,isig)*cvirt )
               mv34(k,isig) = als(2,1)*c2o4pi*( mv34(k,isig) + 
     &                                          mm(k,isig)*cvirt )
               resv(k,isig5) = resv(k,isig5) +
     &                   2*dreal(mm(k,isig)*conjg(mv12(k,isig))) +
     &                   2*dreal(mm(k,isig)*conjg(mv34(k,isig)))
            endif
         enddo 
         if (lnlo) then 
            res(k,-1) = res(k,-1)+ resv(k,-1)
            res(k,1)  = res(k,1) + resv(k,1)
         endif
         res(k,0) = res(k,-1)+res(k,1)
         resv(k,0) = resv(k,-1)+resv(k,1)
         if (lnlo) then
           tree(k) = res(k,0) - resv(k,0)
         else
           tree(k) = res(k,0) 
         endif
      enddo

C  Assign the final results for the subprocesses, include color factor 9
c  for sum over external quark colors

      if (ibos.eq.2) then
c eliminate processes with photon virtuality below cutoff
         if ( abs(p21(4)).lt.qsqAmin .or.  
     1        abs(p43(4)).lt.qsqAmin ) then
            do i = 1,4
               res(i,0) = 0
               tree(i) = 0
            enddo
         endif
         uucc = res(1,0)*9
         uuss = res(2,0)*9
         ddcc = res(3,0)*9
         ddss = res(4,0)*9
         udsc = res(5,0)*9
         ducs = res(6,0)*9
      elseif (ibos.eq.3) then
c eliminate processes with photon virtuality below cutoff
         if ( abs(p21(4)).lt.qsqAmin ) then
            res(1,0) = 0
            res(2,0) = 0
            tree(1) = 0
            tree(2) = 0
         elseif ( abs(p43(4)).lt.qsqAmin ) then
            res(3,0) = 0
            res(4,0) = 0
            tree(3) = 0
            tree(4) = 0
         endif
         uucs = res(1,0)*9
         ddcs = res(2,0)*9
         udcc = res(3,0)*9
         udss = res(4,0)*9
      elseif (ibos.eq.4) then
c eliminate processes with photon virtuality below cutoff
         if ( abs(p21(4)).lt.qsqAmin ) then
            res(1,0) = 0
            res(2,0) = 0
            tree(1) = 0
            tree(2) = 0
         elseif ( abs(p43(4)).lt.qsqAmin ) then
            res(3,0) = 0
            res(4,0) = 0
            tree(3) = 0
            tree(4) = 0
         endif
         uusc = res(1,0)*9
         ddsc = res(2,0)*9
         ducc = res(3,0)*9
         duss = res(4,0)*9
      elseif (ibos.eq.1) then
c eliminate processes with photon virtuality below cutoff
         if ( abs(p21(4)).lt.qsqAmin .or.  
     $        abs(p43(4)).lt.qsqAmin ) then
            do i = 1,4
               res(i,0) = 0
               tree(i) = 0
            enddo
         endif
         uucc = res(1,0)*9
         uuss = res(2,0)*9
         ddcc = res(3,0)*9
         ddss = res(4,0)*9
         udsc = res(5,0)*9
         ducs = res(6,0)*9
      endif

      if (ldoblha) then
c add colour factor 
        do i = 1,6
          res(i,0) = res(i,0)*9
          tree(i)  = tree(i) *9
        enddo

        if (blha_amptype(blha_curproc).eq.0) then ! tree
          blha_amp(1)=blha_amp(1)+
     &      res(blha_idsubproc(blha_cursubproc,blha_curproc),0)

        else if (blha_amptype(blha_curproc).eq.1) then ! loop
          lnmusqp1p2 = log(blha_scale**2/(-p21(4)))
          lnmusqp3p4 = log(blha_scale**2/(-p43(4)))

          blha_amp(1)=blha_amp(1)+
     &      blha_CF/(2d0*pi)*2d0*
     &       (als(1,1)+als(2,1))*(-1d0)*
     &       tree(blha_idsubproc(blha_cursubproc,blha_curproc))
          blha_amp(2)=blha_amp(2)+
     &      1d0/(2d0*pi)*2d0*
     &       ( als(1,1)*(-blha_gammaQuark-blha_CF*lnmusqp1p2)
     &        +als(2,1)*(-blha_gammaQuark-blha_CF*lnmusqp3p4) )*
     &       tree(blha_idsubproc(blha_cursubproc,blha_curproc))
c this is more complicated than necessary, 
c but should work for other VBF processes as well
          blha_amp(3)= blha_amp(3)+
     &      (res(blha_idsubproc(blha_cursubproc,blha_curproc),0)   ! tree + loop + I
     &       -tree(blha_idsubproc(blha_cursubproc,blha_curproc)))  ! -tree
     &       * blha_CF/c2                                          ! correct for Nc ?= inf
     &      + (-(als(1,1)+als(2,1))/2d0/pi*                        ! -I
     &          (blha_CF*cvirt+
     &            2d0*(blha_KQuark+blha_gammaQuark+
     &                 blha_tgammaQuark+(pi**2/6d0-1)*blha_CF) )
     &         +als(1,1)/pi*                                       ! +(mu^2/2p1.p2)^epsilon
     &         (-blha_gammaQuark*lnmusqp1p2-blha_CF/2d0*lnmusqp1p2**2)
     &         +als(2,1)/pi*                                       ! +(mu^2/2p3.p4)^epsilon
     &         (-blha_gammaQuark*lnmusqp3p4-blha_CF/2d0*lnmusqp3p4**2)
     &        )*tree(blha_idsubproc(blha_cursubproc,blha_curproc))
          blha_amp(4)= blha_amp(4)+
     &       tree(blha_idsubproc(blha_cursubproc,blha_curproc))

        else if (blha_amptype(blha_curproc).eq.2) then ! cctree
          ! 1-2
          call BLHA_cctree(1,2,
     &      -blha_CF*
     &       tree(blha_idsubproc(blha_cursubproc,blha_curproc))
     &    )
          ! 3-4
          call BLHA_cctree(3,4,
     &      -blha_CF*
     &       tree(blha_idsubproc(blha_cursubproc,blha_curproc))
     &    )

        else 
          call BLHA_amptypeerror(blha_amptype(blha_curproc),__FILE__,__LINE__)
        endif
      endif
      
      return
      end


#ifdef WITH_NLO 
c
c----------  boxline(id, ... ) ----------------------------------------
c
      subroutine boxline(id, psi1,psi2,k1,k2, isig,
     1                   eps1,eps2,q1,q2, 
     2                   mborn,    mvirt)
      implicit none
      integer id, isig
      double precision k1(0:3),k2(0:3), q1(0:4),q2(0:4)
      double complex psi1(2),psi2(2),eps1(0:3),eps2(0:3),mborn,mvirt
c 
c  determine the finite virtual corrections to the fundamental box-type
c  graph, i.e. the sum of all virtual corrections to the Born diagram
c
c    psi1(k1)    ---->---------->-------------->--  psi2(k2)
c                         S              S
c                         S              S
c                         S              S
c                     q1  V           q2 V
c                         S              S
c                         S              S
c
c                      eps1(q1)       eps2(q2)
c
c   note: k1 = k2+q1+q2, k1.k1 = 0 = k2.k2, q1.q1,q2.q2 arbitrary
c
c  INPUT:      id = 1...8      identifier for graph to determine 
c                              correct D_ij function etc., see BCD_fill
c              psi1(2)         the incoming Weyl spinor
c              psi2(2)         the outgoing Weyl spinor, psi^+(k2)
c       k1(0:3),k2(0:3)        the corresponding momenta in direction of
c                              the fermion line
c              isig            chirality of Weyl fermion (=tau in HZ notation)
c       eps1,eps2(0:3)         the "polarization vectors/currents" 
c                              corresponding to the 2 external EW bosons
c       q1(0:4),q2(0:4)        outflowing momenta of EW bosons 
c              mborn           the amplitude for the corresponding Born diagram
c              alfas           value of the QCD coupling at the appropiate
c                              renormalization scale
c  OUTPUT:     mvirt           the amplitude for the finite part of all 
c                              virtual 1-loop QCD corrections
c
Common block for D_ij etc. information, calculated in BCD_fill
      double complex D0v(8), Dijv(3,13,8), Teps1(8), Teps2(8), Tborn(8),
     1               Tg1(8), Tg2(8)
      double precision ss(8), tt(8), q1sq(8), q2sq(8)
      double precision aq1(0:3,8), aq2(0:3,8), ak1(0:3,8), ak2(0:3,8)
      common /bcd_qqv/ D0v,Dijv, Teps1,Teps2,Tborn,Tg1,Tg2,
     1                 ss,tt,q1sq,q2sq
      common /bcd_qqv_aux/ aq1,aq2,ak1,ak2
Common block bcd_qqv_aux is used for debugging purposes only
      logical ldebug, gaugecheck, loutput
      parameter (ldebug=.false., loutput = .false.)
common block for debugging via bcd_fill_c
      double complex D0c, Dijc(3,13), Teps1c, Teps2c, Tbornc,Tg1c,Tg2c
      double precision sc, tc, q1sqc, q2sqc
      common /bcd_qqv_c/ D0c,Dijc,Teps1c,Teps2c,Tbornc,Tg1c,Tg2c, 
     1                   sc,tc,q1sqc,q2sqc
C   local variables
      integer mu
      double precision s,t,mom(0:3),dot0p
      external dot0p
      double complex me1, me2, mq, d0t, s1c, s1r, sc3, z
      external d0t, s1c, s1r, sc3
      double complex ce1,ce2,cq,cborn
      double complex d0, d11,d12,d13, d21,d22,d23,d24,d25,d26,d27
      double complex d31,d32,d33,d34,d35,d36,d37,d38,d39,d310,
     1               d311,d312,d313
      double complex eps1k1, eps1k2, eps1q2, eps1eps2, 
     1               eps2k1, eps2q1, eps2k2, eps1q1, eps2q2,
     2               gce1, gce2, gcq

Check consistency of kinematics
      if (ldebug) then
         print*," "
         print*," comparison of input to BOXLINE for id = ",id
         print("(a,3f16.4)")," q1sq = ",q1sq(id),q1(4),q1sq(id)/q1(4)
         print("(a,3f16.4)")," q2sq = ",q2sq(id),q2(4),q2sq(id)/q2(4)
         s = -2*dot0p(k1,k2)
         t =  2*dot0p(k2,q2)+q2(4)
         print("(a,3f16.4)")," s = ",ss(id),s,ss(id)/s
         print("(a,3f16.4)")," t = ",tt(id),t,tt(id)/t
c         do mu = 0,3
c            if (abs(k1(mu)-ak1(mu,id)).gt.(1d-10*k1(0)) .or.
c     1          abs(k2(mu)-ak2(mu,id)).gt.(1d-10*k2(0)) .or.
c     1          abs(q1(mu)-aq1(mu,id)).gt.(1d-10*k2(0)) .or.
c     1          abs(q2(mu)-aq2(mu,id)).gt.(1d-10*k2(0)) ) then
c               print*," WARNING: momentum mismatch in boxline, mu =",mu
c               print*," k1: ",k1(mu)/(ak1(mu,id)+1d-80)
c               print*," k2: ",k2(mu)/(ak2(mu,id)+1d-80)
c               print*," q1: ",q1(mu)/(aq1(mu,id)+1d-80)
c               print*," q2: ",q2(mu)/(aq2(mu,id)+1d-80)
c            endif
c         enddo
         print*," input psi1 ",psi1
         print*," input psi2 ",psi2
         print*," input eps1 ",eps1
         print*," input eps2 ",eps2
      endif

c  get matrix elements 
c
c  me1 = psibar(2) eps1slash psi(1)
c  me2 = psibar(2) eps2slash psi(1)
c  mq  = psibar(2) (q1-q2)slash psi(1)
      me1 = -s1c(psi2,eps1,.true.,isig,psi1)
      me2 = -s1c(psi2,eps2,.true.,isig,psi1)
      mq  = -2*s1r(psi2,q1,.true.,isig,psi1)
      if (ldebug) then
         print*," check M_q: ",mq/s1r(psi2,q2,.true.,isig,psi1),
     1          " = 2?"
         do mu = 0,3
            mom(mu) = k2(mu)+q2(mu)
         enddo
         t = mom(0)**2-mom(1)**2-mom(2)**2-mom(3)**2
         z = -sc3(psi2,eps2,mom,eps1,psi1,isig)/t
         print*," M_born:    ",mborn/z," = 1? Abs(ratio) = ",
     &                     abs(mborn/z)
         print*," D0 = ",d0v(id)
      endif
c
c prepare the coefficients
      s = ss(id)
      t = tt(id)

      eps1k1 = eps1(0)*k1(0)-eps1(1)*k1(1)-eps1(2)*k1(2)-eps1(3)*k1(3)
      eps1k2 = eps1(0)*k2(0)-eps1(1)*k2(1)-eps1(2)*k2(2)-eps1(3)*k2(3)
      eps1q2 = eps1(0)*q2(0)-eps1(1)*q2(1)-eps1(2)*q2(2)-eps1(3)*q2(3)
      eps1eps2 = eps1(0)*eps2(0)-eps1(1)*eps2(1)-
     1           eps1(2)*eps2(2)-eps1(3)*eps2(3)
      eps2k1 = eps2(0)*k1(0)-eps2(1)*k1(1)-eps2(2)*k1(2)-eps2(3)*k1(3)
      eps2k2 = eps2(0)*k2(0)-eps2(1)*k2(1)-eps2(2)*k2(2)-eps2(3)*k2(3)
      eps2q1 = eps2(0)*q1(0)-eps2(1)*q1(1)-eps2(2)*q1(2)-eps2(3)*q1(3)

      eps1q1 = eps1(0)*q1(0)-eps1(1)*q1(1)-eps1(2)*q1(2)-eps1(3)*q1(3)
      eps2q2 = eps2(0)*q2(0)-eps2(1)*q2(1)-eps2(2)*q2(2)-eps2(3)*q2(3)

      gaugecheck = abs(eps1q1).gt.1d-10*abs(eps1q2) .or.
     1             abs(eps2q2).gt.1d-10*abs(eps2q1)
      if (ldebug) print*," gaugecheck = ",gaugecheck
      gaugecheck = gaugecheck .or. ldebug
c
      d0 = d0v(id)

      d11 = dijv(1,1,id)
      d12 = dijv(1,2,id)
      d13 = dijv(1,3,id)

      d21 = dijv(2,1,id)
      d22 = dijv(2,2,id)
      d23 = dijv(2,3,id)
      d24 = dijv(2,4,id)
      d25 = dijv(2,5,id)
      d26 = dijv(2,6,id)
      d27 = dijv(2,7,id)

      d31 = dijv(3,1,id)
      d32 = dijv(3,2,id)
      d33 = dijv(3,3,id)
      d34 = dijv(3,4,id)
      d35 = dijv(3,5,id)
      d36 = dijv(3,6,id)
      d37 = dijv(3,7,id)
      d38 = dijv(3,8,id)
      d39 = dijv(3,9,id)
      d310 = dijv(3,10,id)
      d311 = dijv(3,11,id)
      d312 = dijv(3,12,id)
      d313 = dijv(3,13,id)
      if (ldebug) then
         call BCD_fill_c(k1,k2,q1,q2)
         print*," D0: ",abs(d0/d0c-1)
         print("(a,3g10.2)")," D1j ",abs(d11/dijc(1,1)-1),
     1                               abs(d12/dijc(1,2)-1),
     2                               abs(d13/dijc(1,3)-1)
         print("(a,7g10.2)")," D2j ",abs(d21/dijc(2,1)-1),
     1                               abs(d22/dijc(2,2)-1),
     1                               abs(d23/dijc(2,3)-1),
     1                               abs(d24/dijc(2,4)-1),
     1                               abs(d25/dijc(2,5)-1),
     1                               abs(d26/dijc(2,6)-1),
     1                               abs(d27/dijc(2,7)-1)
         print("(a,7g10.2)")," D3j ",abs(d31/dijc(3,1)-1),
     1                               abs(d32/dijc(3,2)-1),
     1                               abs(d33/dijc(3,3)-1),
     1                               abs(d34/dijc(3,4)-1),
     1                               abs(d35/dijc(3,5)-1),
     1                               abs(d36/dijc(3,6)-1),
     1                               abs(d37/dijc(3,7)-1)
         print("(a,7g10.2)"),"     ",abs(d38/dijc(3,8)-1),
     1                               abs(d39/dijc(3,9)-1),
     1                               abs(d310/dijc(3,10)-1),
     1                               abs(d311/dijc(3,11)-1),
     1                               abs(d312/dijc(3,12)-1),
     1                               abs(d313/dijc(3,13)-1)
      endif
c
      ce1 = -2*(d12+d24)*eps2k2*(q1sq(id)+q2sq(id)-3*s-4*t) - 
     1      2*(d12-d24)*eps2q1*(q2sq(id)-t) 
     2  +4*(- d11*eps2k2*s - d12*eps2k1*t + d13*eps2k2*
     3  (q2sq(id)-s-t) + d13*eps2q1*q2sq(id) - d21*eps2k2*s 
     4 - d22*eps2k2*t - d22*eps2q1*q2sq(id) + d23*eps2k2*q1sq(id)
     5 + d25*eps2k2*(q2sq(id)-s-2*t) 
     6 - d26*eps2k2*(q2sq(id)-s-t)+d26*eps2q1*t+2*d27*eps2q1
     7 - d32*eps2k2*q2sq(id) - d34*eps2k2*(q2sq(id)-t) 
     8 + d35*eps2k2*(q2sq(id)-s-t)+d36*eps2k2*(2*q2sq(id)-t) 
     9 + d37*eps2k2*q1sq(id) + d38*eps2k2*(q1sq(id)+q2sq(id)-s)
     & - d39*eps2k2*q1sq(id) 
     1 - d310*eps2k2*(q1sq(id)+2*q2sq(id)-2*s-t)
     2 - 4*d311*eps2k2 + 6*d312*eps2k2 + 2*d313*eps2q1)

      ce2= -2*(d12+d24)*(eps1k2*(q1sq(id)+q2sq(id)-s-2*t)
     1  +eps1q2*(q2sq(id)-s-3*t)) + 
     2  4*(d13*eps1k2*(q1sq(id)-2*s-t) - 
     3  d13*eps1q2*(2*s+t) + d22*eps1k2*q2sq(id)
     4  + d22*eps1q2*q2sq(id) - d23*eps1k2*t + 
     5    d23*eps1q2*(q1sq(id)-t) - d24*eps1k1*q2sq(id) 
     6  + d25*eps1k2*q1sq(id) + d25*eps1k1*(q2sq(id)-2*s-t) 
     7  + d26*eps1k2*t - d26*eps1k1*(q1sq(id)-s) - 
     8  2*d27*eps1q2 + d33*eps1k2*q1sq(id) + d33*eps1q2*q1sq(id)
     9  + d37*eps1k1*(q2sq(id)-s-t) + d38*eps1k1*q2sq(id)
     &  - d39*eps1k1*(q1sq(id)+q2sq(id)-s) - d310*eps1k1*(q2sq(id)-t)
     1  +2*d311*eps1k2+2*d312*eps1q2-6*d313*eps1k1)

C NOTE: cq = 1/2 * cq1 in my reduce code
      cq = (d12+d24)*eps1eps2*s + 2*(4*d12*eps2k2*eps1k2 + 
     1 3*d12*eps2k2*eps1q2+d12*eps2q1*eps1k2-4*d13*eps2k2*eps1k2-
     2 2*d13*eps2k2*eps1q2-2*d13*eps2q1*eps1k2-d13*eps1eps2*s+
     3 2*d22*eps2k2*eps1q2-d22*eps1eps2*t-2*d23*eps2q1*eps1k2-
     4 2*d23*eps2q1*eps1q2-d23*eps1eps2*t+6*d24*eps2k2*eps1k2+
     5 3*d24*eps2k2*eps1q2+d24*eps2q1*eps1k2-6*d25*eps2k2*eps1k2-
     6 2*d25*eps2k2*eps1q2-2*d25*eps2q1*eps1k2-d25*eps1eps2*s-
     7 4*d26*eps2k2*eps1q2+4*d26*eps2q1*eps1k2+2*d26*eps2q1*eps1q2+ 
     8 d26*eps1eps2*(s+2*t)-d32*eps1eps2*q2sq(id)+
     9 d33*eps1eps2*q1sq(id)+2*d34*eps2k2*eps1k2-2*d35*eps2k2*eps1k2+
     & 2*d36*eps2k2*eps1q2+d36*eps1eps2*(q2sq(id)-t)-
     1 2*d37*eps2q1*eps1k2+d37*eps1eps2*(q2sq(id)-s-t)+
     2 2*d38*eps2q1*eps1q2+d38*eps1eps2*(q1sq(id)+2*q2sq(id)-s) -
     3 2*d39*eps2q1*eps1q2-d39*eps1eps2*(2*q1sq(id)+q2sq(id)-s) -
     4 2*d310*eps2k2*eps1q2+2*d310*eps2q1*eps1k2-d310*eps1eps2*
     5 (2*q2sq(id)-s-2*t)+4*d312*eps1eps2-4*d313*eps1eps2)

      cborn = -2*((d36+d37-2*d310)*(q2sq(id)-t) + 
     1    d38*(q1sq(id)+2*q2sq(id)) - d39*(2*q1sq(id)+q2sq(id)) )
     2    -2*(d0+d11+d12-2*d13+d24-d25+d26-d37-d38+d39+d310)*s
     3    -2*(-(d22+d23-2*d26)*t + 2*d27 - d32*q2sq(id)
     4    +d33*q1sq(id)+ 6*(d312-d313))

      if (gaugecheck) then
         gce1 = -2*eps2q2*( d12*(2*q1sq(id)+q2sq(id)-s-3*t) - 
     1       2*(d13-d22)*(q1sq(id)+q2sq(id)-s-t)
     2     - d24*(q2sq(id)-s-t)-2*d26*q1sq(id)-4*d27-4*d312)

         gce2 = 2*eps1q1*(d12*(q2sq(id)+t)-2*d13*(q2sq(id)+s) + 
     1      2*d22*q2sq(id) + 2*d23*(q1sq(id)+q2sq(id)-s-2*t)
     2     +d24*(q2sq(id)+t)+2*d25*(-q2sq(id)+s+t)-2*d26*(q2sq(id)-t)
     3     -4*d27+2*d33*q1sq(id) + 4*d313)

         gcq = 4*eps1q1*eps2q2*(d12-d13+d22-d23+d38-d39) +
     1        2*eps1q1*eps2k2*(d12-4*d23+d24+2*d26-2*d37+2*d310) -
     2        4*eps1q1*eps2q1*(d23-d26+d33-d39) +
     3        2*eps2q2*eps1k2*(3*d12-4*d13+4*d22+d24-2*d25-
     4                                     2*d26+2*d36-2*d310) +
     5        4*eps2q2*eps1q2*(d12-d13+2*d22-2*d26+d32-d38)
         ce1 = ce1 + gce1
         ce2 = ce2 + gce2
         cq = cq + gcq
      endif
c
c  now add the terms from the vertex and propagator corrections
      ce1 = ce1 + 2*eps2k2*Teps2(id)
      ce2 = ce2 + 2*eps1k1*Teps1(id)
      if (gaugecheck) then
         ce1 = ce1 + 2*eps2q2*Tg2(id)
         ce2 = ce2 - 2*eps1q1*Tg1(id)
      endif
      cborn = cborn*t - Tborn(id)
c  and the final result
      mvirt = me1*ce1 + me2*ce2 + mq*cq + mborn*cborn
      if (ldebug) then
         print*," mvirt/mborn ",mvirt/mborn
         if (abs(eps1(0)/q1(0)-1).lt.1d-10) then
            print*," mborn type ",(mborn*cborn+ me2*ce2)/mborn
            print*," mq type    ",(me1*ce1+ mq*cq)/mborn
         else
            print*," mborn type ",(mborn*cborn+ me1*ce1)/mborn
            print*," mq type    ",(me2*ce2+ mq*cq)/mborn
         endif
      endif

      return
      end


#endif
