!*************************************************************************  
!
!     This file contains all functions needed for the
!     momentum mapping of the 
!     Binoth Les Houches Accords 
!     (arXiv:1001.1307 and arXiv:1308.3462)
!
!*************************************************************************  
!
!     Author: Michael Rauch <michael.rauch@kit.edu>
!     Initial version: Sep 2013
!     Last changed: Dec 2013 by Michael Rauch
!
!*************************************************************************  
!   LIST OF ALL FUNCTIONS AND SUBROUTINES IN THIS FILE:
!*************************************************************************
!     SUBROUTINE MomMapping_VBF5(pdgparton, mapparton)
!     SUBROUTINE MomMapping_VBF4(pdgparton, mapparton)
!     SUBROUTINE MomMapping_EW(nparton,nelweak,pdgelweak,mapelweak,oQCD,oQED)
!     SUBROUTINE makefermionpairs(nf,listf,lflip,npairs,listpairs,nother)
!     LOGICAL FUNCTION permutations(nl,list)
!     LOGICAL FUNCTION isparton(pdg)
!     LOGICAL FUNCTION samegen(pdg1,pdg2)
!*************************************************************************  

      SUBROUTINE MomMapping_VBF5(pdgparton, mapparton)
!*************************************************************************
!     Momentum mapping for VBF+j
!*************************************************************************
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

      integer pdgparton(5)
      integer mapparton(5)

! local variables
      integer i,j,k,l,nsp,gsign
      integer npairs, listpairs(5,2), nother
! inverse function to listpairs(:,i), BLHA->diag map
      integer invlistpairs(5)
      integer invmap(5)
! external functions
      integer findpairpartner
      external findpairpartner

      call makefermionpairs(5,pdgparton,.true.,npairs,listpairs,nother)

! exactly one gluon allowed
      if (nother.ne.1) then
        blha_numsubproc(blha_numproc) = -1
        return
      endif

      nsp=0
      do i=1,npairs

! check whether s-channel diagram: both initial-state connected
        if ( (((mapparton(listpairs(1,i))-1)/2.eq.0).and.   &
!           as (12)
              ((mapparton(listpairs(2,i))-1)/2.eq.0)) .or.  &
             (((mapparton(listpairs(3,i))-1)/2.eq.0).and.   &
!           or as (34)
              ((mapparton(listpairs(4,i))-1)/2.eq.0)) ) cycle

        nsp=nsp+1

! inverse listpairs
        do j=1,5
          invlistpairs(listpairs(j,i)) = j
        enddo

! particlemap and fsign
! set gluon to position 5 -- gets overwritten later if this is not true
        blha_particlemap(5,nsp,blha_numproc) = mapparton(listpairs(5,i))
        blha_physdiagmap(5,nsp,blha_numproc) = 5
        blha_fsign(5,nsp,blha_numproc) = +1

        do j=1,2
          blha_particlemap(j,nsp,blha_numproc) = mapparton(j)
! find partner 
          if (pdgparton(j) .eq. 21) then ! gluon
            blha_physdiagmap(j,nsp,blha_numproc) = 5
            blha_fsign(5,nsp,blha_numproc) = -1
! f-s fermion on pair position -- f-s anti-fermion on 5
! -- take the pair not containing the i-s particle =/= j
            l = 2*(2-((invlistpairs(3-j)-1)/2))
            blha_particlemap(j+2,nsp,blha_numproc) = mapparton(listpairs(l,i))
            blha_physdiagmap(j+2,nsp,blha_numproc) = 2*j
            blha_fsign(2*j  ,nsp,blha_numproc) = +1
            blha_particlemap(5  ,nsp,blha_numproc) = mapparton(listpairs(l-1,i))
            blha_physdiagmap(5  ,nsp,blha_numproc) = 2*j-1
            blha_fsign(2*j-1,nsp,blha_numproc) = -1
          else 
! find partner -- "4*((i-1)/2)+3-i" finds the corresponding integer for pair groups (12),(34),...
            l = findpairpartner(invlistpairs(j))
            blha_particlemap(j+2,nsp,blha_numproc) = mapparton(listpairs(l,i))
            if (pdgparton(j) .gt. 0) then
              blha_physdiagmap(j,nsp,blha_numproc) = 2*j-1
              blha_fsign(2*j-1,nsp,blha_numproc) = +1
              k = findpairpartner(2*j-1)
              blha_physdiagmap(j+2,nsp,blha_numproc) = k
              blha_fsign(k,nsp,blha_numproc) = +1
            else 
              blha_physdiagmap(j,nsp,blha_numproc) = 2*j
              blha_fsign(2*j  ,nsp,blha_numproc) = -1
              k = findpairpartner(2*j)
              blha_physdiagmap(j+2,nsp,blha_numproc) = k
              blha_fsign(k,nsp,blha_numproc) = -1
            endif
          endif
        enddo

! create inverse map
        do j=1,5
          do k=1,5
            if (mapparton(k) .eq. blha_particlemap(j,nsp,blha_numproc)) then
              invmap(blha_physdiagmap(j,nsp,blha_numproc)) = k
            endif
          enddo
        enddo

! choose right subproc
        if (abs(pdgparton(invmap(1))).eq. &
            abs(pdgparton(invmap(2)))) then ! upper-Z
          if (abs(pdgparton(invmap(3))).eq. &
              abs(pdgparton(invmap(4)))) then ! NC-ZZ
            if (mod(abs(pdgparton(invmap(1))),2) .eq. 0) then
              if (mod(abs(pdgparton(invmap(3))),2) .eq. 0) then
                blha_idsubproc(nsp,blha_numproc) = 1 ! uucc
              else
                blha_idsubproc(nsp,blha_numproc) = 2 ! uuss
              endif
            else
              if (mod(abs(pdgparton(invmap(3))),2) .eq. 0) then
                blha_idsubproc(nsp,blha_numproc) = 3 ! ddcc
              else
                blha_idsubproc(nsp,blha_numproc) = 4 ! ddss
              endif
            endif
          else ! CC-ZW
            if (mod(abs(pdgparton(invmap(1))),2) .eq. 0) then
              blha_idsubproc(nsp,blha_numproc) = 1 ! uucs/uusc
            else
              blha_idsubproc(nsp,blha_numproc) = 2 ! ddcs/ddsc
            endif
          endif
        else ! upper-W
          if (abs(pdgparton(invmap(3))).eq. &
              abs(pdgparton(invmap(4)))) then ! CC-WZ
            if (mod(abs(pdgparton(invmap(3))),2) .eq. 0) then
              blha_idsubproc(nsp,blha_numproc) = 3 ! udcc/ducc
            else
              blha_idsubproc(nsp,blha_numproc) = 4 ! udss/duss
            endif
          else ! NC-WW
            if (mod(abs(pdgparton(invmap(1))),2) .eq. 0) then
              blha_idsubproc(nsp,blha_numproc) = 5 ! udsc
            else
              blha_idsubproc(nsp,blha_numproc) = 6 ! ducs
            endif
          endif
        endif

! choose where gluon can be attached
        if (blha_fsign(5,nsp,blha_numproc) .eq. +1) then   !final-state gluon -> both ok
          blha_idsubproc(nsp,blha_numproc) =  &
            blha_idsubproc(nsp,blha_numproc) + 10
        else if ( (blha_fsign(1,nsp,blha_numproc) .eq.  &
                  -blha_fsign(2,nsp,blha_numproc) ) .and. &
                  (blha_fsign(3,nsp,blha_numproc) .eq.  &
                   blha_fsign(4,nsp,blha_numproc) ) ) then ! only to (12)
          blha_idsubproc(nsp,blha_numproc) = &
            blha_idsubproc(nsp,blha_numproc) + 20
        else if ( (blha_fsign(1,nsp,blha_numproc) .eq.  &
                   blha_fsign(2,nsp,blha_numproc) ) .and. &
                  (blha_fsign(3,nsp,blha_numproc) .eq.  &
                  -blha_fsign(4,nsp,blha_numproc) ) ) then ! only to (34)
          blha_idsubproc(nsp,blha_numproc) =  &
            blha_idsubproc(nsp,blha_numproc) + 30
        else
! something went wrong here -- should not happen
          blha_idsubproc(nsp,blha_numproc) = 0
        endif
      enddo

      blha_numsubproc(blha_numproc) = nsp



      return
      end
      
!*************************************************************************  

      SUBROUTINE MomMapping_VBF4(pdgparton, mapparton)
!*************************************************************************
!     Momentum mapping for VBF
!*************************************************************************
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

      integer pdgparton(4)
      integer mapparton(4)

! local variables
      integer i,j,k,l,nsp
      integer npairs, listpairs(4,2), nother
! inverse function to listpairs(:,i), BLHA->diag map
      integer invlistpairs(4)
      integer invmap(4)
! external functions
      integer findpairpartner
      external findpairpartner

      call makefermionpairs(4,pdgparton,.true.,npairs,listpairs,nother)

! only all-quarks allowed
      if (nother.ne.0) then
        blha_numsubproc(blha_numproc) = -1
        return
      endif

      nsp=0
      do i=1,npairs

! check whether s-channel diagram: both initial-state connected
        if ( (((mapparton(listpairs(1,i))-1)/2.eq.0).and.   &
!           as (12)
              ((mapparton(listpairs(2,i))-1)/2.eq.0)) .or.  &
             (((mapparton(listpairs(3,i))-1)/2.eq.0).and.   &
!           or as (34)
              ((mapparton(listpairs(4,i))-1)/2.eq.0)) ) cycle

        nsp=nsp+1

! inverse listpairs
        do j=1,4
          invlistpairs(listpairs(j,i)) = j
        enddo

! particlemap and fsign
        do j=1,2
          blha_particlemap(j,nsp,blha_numproc) = mapparton(j)
! find partner -- "4*((i-1)/2)+3-i" finds the corresponding integer for pair groups (12),(34),...
          l = findpairpartner(invlistpairs(j))
          blha_particlemap(j+2,nsp,blha_numproc) = mapparton(listpairs(l,i))
          if (pdgparton(j) .gt. 0) then
            blha_physdiagmap(j,nsp,blha_numproc) = 2*j-1
            blha_fsign(2*j-1,nsp,blha_numproc) = +1
            k = findpairpartner(2*j-1)
            blha_physdiagmap(j+2,nsp,blha_numproc) = k
            blha_fsign(k,nsp,blha_numproc) = +1
          else 
            blha_physdiagmap(j,nsp,blha_numproc) = 2*j
            blha_fsign(2*j  ,nsp,blha_numproc) = -1
            k = findpairpartner(2*j)
            blha_physdiagmap(j+2,nsp,blha_numproc) = k
            blha_fsign(k,nsp,blha_numproc) = -1
          endif
        enddo

! create inverse map
        do j=1,4
          do k=1,4
            if (mapparton(k) .eq. blha_particlemap(j,nsp,blha_numproc)) then
              invmap(blha_physdiagmap(j,nsp,blha_numproc)) = k
            endif
          enddo
        enddo

! choose right subproc
        if (abs(pdgparton(invmap(1))).eq. &
            abs(pdgparton(invmap(2)))) then ! upper-Z
          if (abs(pdgparton(invmap(3))).eq. &
              abs(pdgparton(invmap(4)))) then ! NC-ZZ
            if (mod(abs(pdgparton(invmap(1))),2) .eq. 0) then
              if (mod(abs(pdgparton(invmap(3))),2) .eq. 0) then
                blha_idsubproc(nsp,blha_numproc) = 1 ! uucc
              else
                blha_idsubproc(nsp,blha_numproc) = 2 ! uuss
              endif
            else
              if (mod(abs(pdgparton(invmap(3))),2) .eq. 0) then
                blha_idsubproc(nsp,blha_numproc) = 3 ! ddcc
              else
                blha_idsubproc(nsp,blha_numproc) = 4 ! ddss
              endif
            endif
          else ! CC-ZW
            if (mod(abs(pdgparton(invmap(1))),2) .eq. 0) then
              blha_idsubproc(nsp,blha_numproc) = 1 ! uucs/uusc
            else
              blha_idsubproc(nsp,blha_numproc) = 2 ! ddcs/ddsc
            endif
          endif
        else ! upper-W
          if (abs(pdgparton(invmap(3))).eq. &
              abs(pdgparton(invmap(4)))) then ! CC-WZ
            if (mod(abs(pdgparton(invmap(3))),2) .eq. 0) then
              blha_idsubproc(nsp,blha_numproc) = 3 ! udcc/ducc
            else
              blha_idsubproc(nsp,blha_numproc) = 4 ! udss/duss
            endif
          else ! NC-WW
            if (mod(abs(pdgparton(invmap(1))),2) .eq. 0) then
              blha_idsubproc(nsp,blha_numproc) = 5 ! udsc
            else
              blha_idsubproc(nsp,blha_numproc) = 6 ! ducs
            endif
          endif
        endif
      enddo

      blha_numsubproc(blha_numproc) = nsp

      return
      end

!*************************************************************************  

      SUBROUTINE MomMapping_EW(nparton,nelweak,pdgelweak,mapelweak, &
                               oQCD,oQED)
!*************************************************************************
!     Map electro-weak part of VBF processes to process IDs
!*************************************************************************
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/process.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

      integer nparton,nelweak
      integer pdgelweak(nelweak)
      integer mapelweak(nelweak)
      integer oQCD,oQED
! local variables
      integer i,j,k,l,m
      integer nsp
      integer npairs, listpairs(nelweak,6), nother
      integer nwp, nwm, nzl, nzn, nh, ngamma, nhel, nid
      logical lVBFlepdec, lVBFhaddec, lQCDlepdec

      nsp = blha_numsubproc(blha_numproc)
      if (blha_numsubproc(blha_numproc) .le. 0) return

      call makefermionpairs(nelweak,pdgelweak,.false., &
                            npairs,listpairs,nother)

      lVBFlepdec = (oQCD .eq. nparton-4) .and. (oQED .eq. nelweak+2 )
      lVBFhaddec = (oQCD .eq. nparton-6) .and. (oQED .eq. nelweak+4 )

      blha_alphasorder(blha_numproc) = oQCD
      blha_alphaorder(blha_numproc)  = oQED

      do i=1,npairs
! get number of W,Z,H,gamma
        nwp = 0
        nwm = 0
        nzl = 0
        nzn = 0
        nh  = 0
        ngamma = 0
        nhel = 1
        do j=1,nelweak-nother,2
          if (pdgelweak(listpairs(j,i))+ &
              pdgelweak(listpairs(j+1,i)).eq.0) then
            if (mod(abs(pdgelweak(listpairs(j,i))),2).eq.1) then
              nzl = nzl + 1
              nhel = nhel * 2
            else
              nzn = nzn + 1
            endif
          else if (pdgelweak(listpairs(j,i))+ &
                   pdgelweak(listpairs(j+1,i)).eq.+1) then
            nwp = nwp + 1
          else if (pdgelweak(listpairs(j,i))+ &
                   pdgelweak(listpairs(j+1,i)).eq.-1) then
            nwm = nwm + 1
          endif
          finallep(j)   = pdgelweak(listpairs(j,i))
          finallep(j+1) = pdgelweak(listpairs(j+1,i))
        enddo
        do j=nelweak-nother+1,nelweak
          if (pdgelweak(listpairs(j,i)).eq.25) then
            nh = nh + 1
          else if (pdgelweak(listpairs(j,i)).eq.22) then
            ngamma = ngamma + 1
            nhel = nhel * 2
          endif
        enddo
        nid = nwp   *10**5 &
             +nwm   *10**4 &
             +nzl   *10**3 &
             +nzn   *10**2 &
             +nh    *10**1 &
             +ngamma*10**0

        blha_multsubproc((i-1)*nsp+1,blha_numproc) = nhel

        SELECT CASE (nid)
        CASE(10) ! H
          if ((nparton.eq.4).and.lVBFlepdec) then
            blha_procsubproc((i-1)*nsp+1,blha_numproc) = Hjj
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
          else if ((nparton.eq.5).and.lVBFlepdec) then
            blha_procsubproc((i-1)*nsp+1,blha_numproc) = Hjjj
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
          else if ((nparton.eq.6).and.lVBFlepdec) then
            blha_procsubproc((i-1)*nsp+1,blha_numproc) = Hjjj
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .true.
          else 
            blha_numsubproc(blha_numproc) = -1
            return
          endif
        CASE(100000) ! W+
          if (lVBFlepdec) then
            blha_procsubproc((i-1)*nsp+1,blha_numproc) = Wpjj
            if (nparton.eq.4) then
              blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
            else if (nparton.eq.5) then
              blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .true.
            endif
          else 
            blha_numsubproc(blha_numproc) = -1
            return
          endif
        CASE(10000) ! W-
          if (lVBFlepdec) then
            blha_procsubproc((i-1)*nsp+1,blha_numproc) = Wmjj
            if (nparton.eq.4) then
              blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
            else if (nparton.eq.5) then
              blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .true.
            endif
          else 
            blha_numsubproc(blha_numproc) = -1
            return
          endif
        CASE(1000) ! Z_l
          blha_procsubproc((i-1)*nsp+1,blha_numproc) = Zjj_l
          if ((nparton.eq.4).and.lVBFlepdec) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
            blha_multsubproc((i-1)*nsp+1,blha_numproc) = 1      ! helicity summation handled internally in this process
          else if ((nparton.eq.5).and.lVBFlepdec) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .true.
            blha_multsubproc((i-1)*nsp+1,blha_numproc) = 1      ! helicity summation handled internally in this process
          else 
            blha_numsubproc(blha_numproc) = -1
            return
          endif
        CASE(100) ! Z_nu
          blha_procsubproc((i-1)*nsp+1,blha_numproc) = Zjj_nu
          if ((nparton.eq.4).and.lVBFlepdec) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
          else if ((nparton.eq.5).and.lVBFlepdec) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .true.
          else 
            blha_numsubproc(blha_numproc) = -1
            return
          endif
        CASE(1) ! gamma
          blha_procsubproc((i-1)*nsp+1,blha_numproc) = Ajj
          if ((nparton.eq.4).and.lVBFlepdec) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
            blha_multsubproc((i-1)*nsp+1,blha_numproc) = 1      ! helicity summation handled internally in this process
          else if ((nparton.eq.5).and.lVBFlepdec) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .true.
            blha_multsubproc((i-1)*nsp+1,blha_numproc) = 1      ! helicity summation handled internally in this process
          else 
            blha_numsubproc(blha_numproc) = -1
            return
          endif
        CASE(20) ! HH
          blha_procsubproc((i-1)*nsp+1,blha_numproc) = HHjj
          if ((nparton.eq.4).and.lVBFlepdec) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
          else if ((nparton.eq.5).and.lVBFlepdec) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .true.
          else 
            blha_numsubproc(blha_numproc) = -1
            return
          endif
        CASE(11) ! H gamma
          blha_procsubproc((i-1)*nsp+1,blha_numproc) = HAjj
          if ((nparton.eq.4).and.lVBFlepdec) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
            blha_multsubproc((i-1)*nsp+1,blha_numproc) = 1      ! helicity summation handled internally in this process
          else if ((nparton.eq.5).and.lVBFlepdec) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .true.
            blha_multsubproc((i-1)*nsp+1,blha_numproc) = 1      ! helicity summation handled internally in this process
          else 
            blha_numsubproc(blha_numproc) = -1
            return
          endif
        CASE(110000) ! W+ W-
          blha_procsubproc((i-1)*nsp+1,blha_numproc) = WpWmjj
          if ((nparton.eq.4).and.lVBFlepdec) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
          else if ((nparton.eq.5).and.lVBFlepdec) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .true.
          else 
            blha_numsubproc(blha_numproc) = -1
            return
          endif
        CASE(101000) ! W+ Z_l
          if (lVBFlepdec) then
            blha_procsubproc((i-1)*nsp+1,blha_numproc) = WpZjj
            if (nparton.eq.4) then
              blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
            else if (nparton.eq.5) then
              blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .true.
            endif
          else 
            blha_numsubproc(blha_numproc) = -1
            return
          endif
        CASE(11000) ! W- Z_l
          if (lVBFlepdec) then
            blha_procsubproc((i-1)*nsp+1,blha_numproc) = WmZjj
            if (nparton.eq.4) then
              blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
            else if (nparton.eq.5) then
              blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .true.
            endif
          else 
            blha_numsubproc(blha_numproc) = -1
            return
          endif
        CASE(2000) ! Z_l Z_l
          blha_procsubproc((i-1)*nsp+1,blha_numproc) = ZZjj_ll
          if ((nparton.eq.4).and.lVBFlepdec) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
          else if ((nparton.eq.5).and.lVBFlepdec) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .true.
          else 
            blha_numsubproc(blha_numproc) = -1
            return
          endif
        CASE(1100) ! Z_l Z_nu
          blha_procsubproc((i-1)*nsp+1,blha_numproc) = ZZjj_lnu
          if ((nparton.eq.4).and.lVBFlepdec) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
          else if ((nparton.eq.5).and.lVBFlepdec) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .true.
          else 
            blha_numsubproc(blha_numproc) = -1
            return
          endif
        CASE(100001) ! W+ A
          if (lVBFlepdec) then
            blha_procsubproc((i-1)*nsp+1,blha_numproc) = WpAjj
            if (nparton.eq.4) then
              blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
            else if (nparton.eq.5) then
              blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .true.
            endif
          else 
            blha_numsubproc(blha_numproc) = -1
            return
          endif
        CASE(10001) ! W- A
          if (lVBFlepdec) then
            blha_procsubproc((i-1)*nsp+1,blha_numproc) = WmAjj
            if (nparton.eq.4) then
              blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
            else if (nparton.eq.5) then
              blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .true.
            endif
          else 
            blha_numsubproc(blha_numproc) = -1
            return
          endif
        CASE(1001) ! Z_l A
          if (lVBFlepdec) then
          blha_procsubproc((i-1)*nsp+1,blha_numproc) = ZAjj
            if (nparton.eq.4) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
            else if (nparton.eq.5) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .true.
            endif
          else 
            blha_numsubproc(blha_numproc) = -1
            return
          endif
        CASE(101) ! Z_nu A
          if (lVBFlepdec) then
          blha_procsubproc((i-1)*nsp+1,blha_numproc) = ZAjj_n
            if (nparton.eq.4) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
            else if (nparton.eq.5) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .true.
            endif
          else 
            blha_numsubproc(blha_numproc) = -1
            return
          endif
        CASE(200000) ! W+ W+
          blha_procsubproc((i-1)*nsp+1,blha_numproc) = WpWpjj
          if ((nparton.eq.4).and.lVBFlepdec) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
          else if ((nparton.eq.5).and.lVBFlepdec) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .true.
          else 
            blha_numsubproc(blha_numproc) = -1
            return
          endif
        CASE(20000) ! W- W-
          blha_procsubproc((i-1)*nsp+1,blha_numproc) = WmWmjj
          if ((nparton.eq.4).and.lVBFlepdec) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
          else if ((nparton.eq.5).and.lVBFlepdec) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .true.
          else 
            blha_numsubproc(blha_numproc) = -1
            return
          endif
        CASE DEFAULT
! not yet supported
          blha_numsubproc(blha_numproc) = -1
          return
        END SELECT

! get correct electro-weak order
        procid = blha_procsubproc((i-1)*nsp+1,blha_numproc)
        LOplusJet = blha_lojsubproc((i-1)*nsp+1,blha_numproc)
        call proc_assignment
        call lepton_assignment
        do j=1,nelweak-nother,2                      ! lepton pairs
          k=1
          l=0
          do while (l.eq.0)                          ! found first match
            if (listpairs(k,i).gt.0) then
              l = pdgelweak(listpairs(k,i))
            else
              l = 0
            endif
            if (finallep(j) .ne. l) then 
              k = k+1
              l = 0
            endif
          enddo
          if (mod(k,2).eq.0) then                    ! get partner: next or previous
            l = k-1
          else 
            l = k+1
          endif
          blha_particlemap(nparton+j,(i-1)*nsp+1,blha_numproc) =  &
            mapelweak(listpairs(k,i))
          blha_fsign(nparton+j,(i-1)*nsp+1,blha_numproc)       =  &
            sign(1,pdgelweak(listpairs(k,i)))
          blha_particlemap(nparton+j+1,(i-1)*nsp+1,blha_numproc) = &
            mapelweak(listpairs(l,i))
          blha_fsign(nparton+j+1,(i-1)*nsp+1,blha_numproc)       = &
            sign(1,pdgelweak(listpairs(l,i)))
          listpairs(k,i) = 0                         ! used: clear entries 
          listpairs(l,i) = 0
        enddo
        j=nelweak-nother+1
        do k=nelweak-nother+1,nelweak              ! now add all Higgs 
          if (listpairs(k,i) .gt. 0) then
            if (pdgelweak(listpairs(k,i)).eq.25) then
              blha_particlemap(nparton+j,(i-1)*nsp+1,blha_numproc) =  &
                mapelweak(listpairs(k,i))
              blha_fsign(nparton+j,(i-1)*nsp+1,blha_numproc)       = +1
              listpairs(k,i) = 0                     ! used: clear entries 
              j=j+1
            endif
          endif
        enddo
        do k=nelweak-nother+1,nelweak              ! now add all photons 
          if (listpairs(k,i) .gt. 0) then
            if (pdgelweak(listpairs(k,i)).eq.22) then
              blha_particlemap(nparton+j,(i-1)*nsp+1,blha_numproc) =  &
                mapelweak(listpairs(k,i))
              blha_fsign(nparton+j,(i-1)*nsp+1,blha_numproc)       = +1
              listpairs(k,i) = 0                     ! used: clear entries 
              j=j+1
            endif
          endif
        enddo
! copy to other subprocesses
        do m=1,nsp ! over all partonic subprocesses
          blha_procsubproc((i-1)*nsp+m,blha_numproc) =  &
            blha_procsubproc((i-1)*nsp+1,blha_numproc) 
          blha_lojsubproc((i-1)*nsp+m,blha_numproc) =  &
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) 
          blha_multsubproc((i-1)*nsp+m,blha_numproc) =  &
            blha_multsubproc((i-1)*nsp+1,blha_numproc) 
        enddo
! copy to other parts from corresponding first entry
        do m=1,nsp ! over all partonic subprocesses
          do j=1,nparton
            blha_particlemap(j,(i-1)*nsp+m,blha_numproc) =  &
              blha_particlemap(j,m,blha_numproc)
            blha_fsign(j,(i-1)*nsp+m,blha_numproc)       =  &
              blha_fsign(j,m,blha_numproc) 
            blha_idsubproc((i-1)*nsp+m,blha_numproc)     =  &
              blha_idsubproc(m,blha_numproc) 
          enddo
          do j=nparton+1,nelweak+nparton
            blha_particlemap(j,(i-1)*nsp+m,blha_numproc) =  &
              blha_particlemap(j,(i-1)*nsp+1,blha_numproc)
            blha_fsign(j,(i-1)*nsp+m,blha_numproc)       =  &
              blha_fsign(j,(i-1)*nsp+1,blha_numproc) 
          enddo
        enddo

      enddo

      blha_numsubproc(blha_numproc) = nsp*npairs

      return
      end
!*************************************************************************  

      LOGICAL FUNCTION isparton(pdg)
!*************************************************************************
!     Whether pdg denotes a parton
!*************************************************************************
      implicit none

      integer pdg

      isparton = (abs(pdg)>0) .and. ((abs(pdg)<6) .or. (pdg.eq.21))
      return
      end

!*************************************************************************  

      LOGICAL FUNCTION samegen(pdg1,pdg2)
!*************************************************************************
!     Whether pdg denotes a parton
!*************************************************************************
      implicit none

      integer pdg1, pdg2
      integer apdg1, apdg2

      apdg1 = abs(pdg1)
      apdg2 = abs(pdg2)
      
      if ( ((apdg1.eq.0) .and. (apdg2.ne.0)) .or. &
           ((apdg1.ne.0) .and. (apdg2.eq.0)) ) then
        samegen = .false.
      else
        samegen = (apdg1+1)/2 .eq. (apdg2+1)/2
      endif

      return
      end

!*************************************************************************  

      SUBROUTINE makefermionpairs(nf,listf,lflip,npairs,listpairs,nother)
!*************************************************************************
!     Initialize variables
!*************************************************************************
      implicit none

      integer nf !number of fermions
      integer npairs,nother
      integer listf(nf)
      integer listpairs(nf,*)
      logical lflip

! local variables
      integer i,j
      integer mylistf(nf)
      integer nbra, listbra(nf)
      integer nket, listket(nf)
      integer listother(nf)
      integer nperm(nf)
      logical lpairok, lperm
! functions
      logical permutations, samegen
      external permutations, samegen

! reminder: bra: final-state fermion, initial-state anti-fermion
!           ket: initial-state fermion, final-state anti-fermion

! convert listf to all out-going
      do i=1,nf
        if ((i.le.2) .and. lflip) then
          mylistf(i) = -listf(i)
        else
          mylistf(i) = listf(i)
        endif
      enddo

! sort into bras, kets and other
      nbra=0
      nket=0
      nother=0
      do i=1,nf
        if ((mylistf(i).gt.0) .and. (mylistf(i).lt.20)) then
          nbra=nbra+1
          listbra(nbra)=i
        else if ((mylistf(i).lt.0) .and. (mylistf(i).gt.-20)) then
          nket=nket+1
          listket(nket)=i
        else
          nother=nother+1
          listother(nother)=i
        endif
      enddo

      npairs=0
! number of bras and kets does not match -> invalid
      if (nbra.ne.nket) return 

! kets are taken as ordered; loop over all bra permutations
      do i=1,nbra
        nperm(i) = i
      enddo
      lperm=.true.
      do while(lperm)
! check if pairs are valid
        lpairok=.true.
        do j=1,nket
          if (.not. samegen( &
              listf(listket(j)), &
              listf(listbra(abs(nperm(j)))) &
            )) lpairok=.false.
        enddo
        if (lpairok) then
          npairs=npairs+1
          do j=1,nket
            listpairs(2*j-1,npairs)=listket(j)
            listpairs(2*j,npairs)=listbra(abs(nperm(j)))
          enddo
          do j=1,nother
            listpairs(nbra+nket+j,npairs)=listother(j)
          enddo
        endif

! get next permutation
        lperm = permutations(nbra,nperm)
      enddo

      return
      end

!*************************************************************************  

      LOGICAL FUNCTION permutations(nl,list)
!*************************************************************************
!     Generates all permutations of 1,2,...,n
!     Steinhaus-Johnson-Trotter-Algorithm
!     Even's speedup implemented using signs
!     -> use abs on returned list
!*************************************************************************
      implicit none
      integer nl
      integer list(nl)
      logical linit
! local variables
      integer i,n,tmp
      logical lmobile

! find largest mobile integer
      lmobile = .false.
      do n=nl,1,-1
        do i=1,nl
          if (abs(list(i)).eq.n) then
! check if it's mobile
            if (list(i).gt.0) then
              if (i.gt.1) then
                if (abs(list(i-1)).lt.n) then
                  lmobile=.true.
                  tmp = list(i-1)
                  list(i-1) = list(i)
                  list(i) = tmp
                endif
              endif
            else
              if (i.lt.nl) then
                if (abs(list(i+1)).lt.n) then
                  lmobile=.true.
                  tmp = list(i+1)
                  list(i+1) = list(i)
                  list(i) = tmp
                endif
              endif
            endif
          endif
          if (lmobile) exit
        enddo
        if (lmobile) exit
      enddo
! if mobile is non-largest, swap all directions of larger ones
      if (lmobile.and.(n.ne.nl)) then
        do i=1,nl
          if (abs(list(i)).gt.n) &
            list(i) = -list(i)
        enddo
      endif
      
! if all non-mobile found last one
      if (lmobile) then
        permutations = .true.
      else
        permutations = .false.
      endif

      return
      end

!*************************************************************************

      INTEGER FUNCTION findpairpartner(i)
!*************************************************************************
!     Find corresponding partner to i in pairs
!     (1<->2),(3<->4),...
!*************************************************************************
      implicit none
      integer i

      findpairpartner = 4*((i-1)/2)+3 - i

      return
      end

!*************************************************************************

