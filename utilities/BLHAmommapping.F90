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
!     Last changed: Mar 2017 by Michael Rauch
!
!*************************************************************************  
!   LIST OF ALL FUNCTIONS AND SUBROUTINES IN THIS FILE:
!*************************************************************************
!     SUBROUTINE MomMapping_VBF6(pdgparton, mapparton)
!     SUBROUTINE MomMapping_VBF5(pdgparton, mapparton)
!     SUBROUTINE MomMapping_VBF4(pdgparton, mapparton)
!     SUBROUTINE MomMapping_EW(nparton,nelweak,pdgelweak,mapelweak,oQCD,oQED)
!     SUBROUTINE makefermionpairs(nf,listf,lflip,npairs,listpairs,nother)
!     LOGICAL FUNCTION permutations(nl,list)
!     LOGICAL FUNCTION isparton(pdg)
!     LOGICAL FUNCTION samegen(pdg1,pdg2)
!*************************************************************************  

      SUBROUTINE MomMapping_VBF6(pdgparton, mapparton)
!*************************************************************************
!     Momentum mapping for VBF+2j
!*************************************************************************
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

      integer pdgparton(6)
      integer mapparton(6)

! local variables
      integer i,j,k,l,nsp
      integer npairs, listpairs(6,6), nother
! inverse function to listpairs(:,i), BLHA->diag map
      integer invlistpairs(6)
      integer invmap(6)
! external functions
      integer findpairpartner
      logical samegen
      external findpairpartner, samegen

      call makefermionpairs(6,pdgparton,.true.,npairs,listpairs,nother)

      nsp=0

      if (nother.eq.2) then

        do i=1,npairs
          ! check whether s-channel diagram: both initial-state connected
          !! as (12)
          if ( (((mapparton(listpairs(1,i))-1)/2.eq.0).and.   &
                ((mapparton(listpairs(2,i))-1)/2.eq.0)) .or.  &
          !! or as (34)
               (((mapparton(listpairs(3,i))-1)/2.eq.0).and.   &
                ((mapparton(listpairs(4,i))-1)/2.eq.0)) ) cycle

! inverse listpairs
          do j=1,6
            invlistpairs(listpairs(j,i)) = j
          enddo

          nsp=nsp+1

          ! 1 and 2 stay at positions
          blha_particlemap(1,nsp,blha_numproc) = mapparton(1)
          blha_particlemap(2,nsp,blha_numproc) = mapparton(2)

          ! set gluons to 5 and 6, overwrite later if initial-state
          blha_particlemap(5,nsp,blha_numproc) = mapparton(listpairs(5,i))
          blha_physdiagmap(5,nsp,blha_numproc) = 5
          blha_fsign(5,nsp,blha_numproc) = +1
          blha_particlemap(6,nsp,blha_numproc) = mapparton(listpairs(6,i))
          blha_physdiagmap(6,nsp,blha_numproc) = 6
          blha_fsign(6,nsp,blha_numproc) = +1

          if ( pdgparton(1) .eq. 21 ) then
            ! initial-state gluon - map to 5
            blha_physdiagmap(1,nsp,blha_numproc) = 5
            blha_fsign(5,nsp,blha_numproc) = -1
            if ( pdgparton(2) .eq. 21 ) then
              ! 2 initial-state gluons
              blha_physdiagmap(2,nsp,blha_numproc) = 6
              blha_fsign(6,nsp,blha_numproc) = -1
              ! fermions in pairs
              do j=1,4
                blha_particlemap(j+2,nsp,blha_numproc) = &
                  mapparton(listpairs(j,i))
                blha_physdiagmap(j+2,nsp,blha_numproc) = j
                blha_fsign(j,nsp,blha_numproc) = (-1)**j
              enddo
            else
              ! line containing parton 2
              k = invlistpairs(2)
              l = listpairs(findpairpartner(k),i)
              blha_particlemap(4,nsp,blha_numproc) = mapparton(l)
              if (mod(k,2).eq.1) then
                blha_physdiagmap(2,nsp,blha_numproc) = 3
                blha_physdiagmap(4,nsp,blha_numproc) = 4
                blha_fsign(3,nsp,blha_numproc) = +1
                blha_fsign(4,nsp,blha_numproc) = +1
                l = listpairs(4-k,i) ! ket of other line
              else
                blha_physdiagmap(2,nsp,blha_numproc) = 4
                blha_physdiagmap(4,nsp,blha_numproc) = 3
                blha_fsign(3,nsp,blha_numproc) = -1
                blha_fsign(4,nsp,blha_numproc) = -1
                l = listpairs(5-k,i) ! ket of other line
              endif
              ! final-final line
              k = listpairs(findpairpartner(invlistpairs(l)),i)
              blha_particlemap(5,nsp,blha_numproc) = mapparton(l)
              blha_particlemap(3,nsp,blha_numproc) = mapparton(k)
              blha_physdiagmap(5,nsp,blha_numproc) = 1
              blha_fsign(1,nsp,blha_numproc) = -1
              blha_physdiagmap(3,nsp,blha_numproc) = 2
              blha_fsign(2,nsp,blha_numproc) = +1
            endif
          else if ( pdgparton(2) .eq. 21 ) then
            ! initial-state gluon - map to 5
            blha_physdiagmap(2,nsp,blha_numproc) = 5
            blha_fsign(5,nsp,blha_numproc) = -1
            ! line containing parton 1
            k = invlistpairs(1)
            l = listpairs(findpairpartner(k),i)
            blha_particlemap(3,nsp,blha_numproc) = mapparton(l)
            if (mod(k,2).eq.1) then
              blha_physdiagmap(1,nsp,blha_numproc) = 1
              blha_physdiagmap(3,nsp,blha_numproc) = 2
              blha_fsign(1,nsp,blha_numproc) = +1
              blha_fsign(2,nsp,blha_numproc) = +1
              l = listpairs(4-k,i) ! ket of other line
            else
              blha_physdiagmap(1,nsp,blha_numproc) = 2
              blha_physdiagmap(3,nsp,blha_numproc) = 1
              blha_fsign(1,nsp,blha_numproc) = -1
              blha_fsign(2,nsp,blha_numproc) = -1
              l = listpairs(5-k,i) ! ket of other line
            endif
            ! final-final line
            k = listpairs(findpairpartner(invlistpairs(l)),i)
            blha_particlemap(5,nsp,blha_numproc) = mapparton(l)
            blha_particlemap(4,nsp,blha_numproc) = mapparton(k)
            blha_physdiagmap(5,nsp,blha_numproc) = 3
            blha_fsign(3,nsp,blha_numproc) = -1
            blha_physdiagmap(4,nsp,blha_numproc) = 4
            blha_fsign(4,nsp,blha_numproc) = +1
          else
            do j=1,2
              blha_particlemap(j,nsp,blha_numproc) = mapparton(j)
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
          endif

          ! create inverse map
          do j=1,6
            do k=1,6
              if (mapparton(k) .eq. blha_particlemap(j,nsp,blha_numproc)) then
                invmap(blha_physdiagmap(j,nsp,blha_numproc)) = k
              endif
            enddo
          enddo

          ! choose right subproc
          blha_idsubproc(nsp,blha_numproc) = 100001
          if (abs(pdgparton(invmap(1))).eq. &
              abs(pdgparton(invmap(2)))) then ! upper-Z
            if (abs(pdgparton(invmap(3))).eq. &
                abs(pdgparton(invmap(4)))) then ! NC-ZZ
              if (mod(abs(pdgparton(invmap(1))),2) .eq. 1) then
                blha_idsubproc(nsp,blha_numproc) = &
                  blha_idsubproc(nsp,blha_numproc) + 1000
              endif
              if (mod(abs(pdgparton(invmap(3))),2) .eq. 1) then
                blha_idsubproc(nsp,blha_numproc) = &
                  blha_idsubproc(nsp,blha_numproc) + 100
              endif
            else ! CC-ZW
              if (mod(abs(pdgparton(invmap(1))),2) .eq. 0) then
                blha_idsubproc(nsp,blha_numproc) = 0 ! uucs/uusc
              else
                blha_idsubproc(nsp,blha_numproc) = 0 ! ddcs/ddsc
              endif
            endif
          else ! upper-W
            if (abs(pdgparton(invmap(3))).eq. &
                abs(pdgparton(invmap(4)))) then ! CC-WZ
              if (mod(abs(pdgparton(invmap(3))),2) .eq. 0) then
                blha_idsubproc(nsp,blha_numproc) = 0 ! udcc/ducc
              else
                blha_idsubproc(nsp,blha_numproc) = 0 ! udss/duss
              endif
            else ! NC-WW
              blha_idsubproc(nsp,blha_numproc) = &
                blha_idsubproc(nsp,blha_numproc) + 10000
              if (mod(abs(pdgparton(invmap(1))),2) .eq. 1) then
                blha_idsubproc(nsp,blha_numproc) = &
                  blha_idsubproc(nsp,blha_numproc) + 1000
              endif
            endif
          endif

!          ! check if we've seen this process already -> symmetry factor
!          if ( blha_idsubproc(nsp,blha_numproc) .ne. 0 ) then
!            do j=1,nsp-1
!              if (blha_idsubproc(j,blha_numproc)/10 .eq. &
!                  blha_idsubproc(nsp,blha_numproc)/10 ) then ! ignore last digit with multiplicity information
!                blha_idsubproc(j,blha_numproc) = blha_idsubproc(j,blha_numproc) + 1
!                nsp = nsp-1
!                exit
!              endif
!            enddo
!          endif

        enddo

      else if (nother.eq.0) then

        do i=1,npairs

          ! check whether s-channel diagram: both initial-state connected
          !! as (12)
          if ( (((mapparton(listpairs(1,i))-1)/2.eq.0).and.   &
                ((mapparton(listpairs(2,i))-1)/2.eq.0)) .or.  &
          !! or as (34)
               (((mapparton(listpairs(3,i))-1)/2.eq.0).and.   &
                ((mapparton(listpairs(4,i))-1)/2.eq.0)) .or.  &
          !! or as (56)
               (((mapparton(listpairs(5,i))-1)/2.eq.0).and.   &
                ((mapparton(listpairs(6,i))-1)/2.eq.0)) ) cycle

          nsp=nsp+1

          ! inverse listpairs
          do j=1,6
            invlistpairs(listpairs(j,i)) = j
          enddo

          do j=1,2
            blha_particlemap(j,nsp,blha_numproc) = mapparton(j)
          ! find partner 
            l = findpairpartner(invlistpairs(j))
            blha_particlemap(j+2,nsp,blha_numproc) = mapparton(listpairs(l,i))
            listpairs(invlistpairs(j),i)=0
            listpairs(l,i)=0
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
          ! remaining pair
          do j=1,5,2
            if (listpairs(j,i).eq.0) cycle
            blha_particlemap(5,nsp,blha_numproc) = mapparton(listpairs(j,i))
            blha_particlemap(6,nsp,blha_numproc) = mapparton(listpairs(j+1,i))
            exit
          enddo
          blha_physdiagmap(5,nsp,blha_numproc) = 5
          blha_physdiagmap(6,nsp,blha_numproc) = 6
          blha_fsign(5,nsp,blha_numproc) = -1
          blha_fsign(6,nsp,blha_numproc) = +1

          ! create inverse map
          do j=1,6
            do k=1,6
              if (mapparton(k) .eq. blha_particlemap(j,nsp,blha_numproc)) then
                invmap(blha_physdiagmap(j,nsp,blha_numproc)) = k
              endif
            enddo
          enddo

          ! choose right subproc
          blha_idsubproc(nsp,blha_numproc) = 200000
          if (abs(pdgparton(invmap(1))).eq. &
              abs(pdgparton(invmap(2)))) then ! upper-Z
            if (abs(pdgparton(invmap(3))).eq. &
                abs(pdgparton(invmap(4)))) then ! lower-ZZ
              if (abs(pdgparton(invmap(5))).eq. &
                  abs(pdgparton(invmap(6)))) then ! ZZZ
                if (mod(abs(pdgparton(invmap(1))),2) .eq. 1) then
                  blha_idsubproc(nsp,blha_numproc) = &
                    blha_idsubproc(nsp,blha_numproc) + 1000
                endif
                if (mod(abs(pdgparton(invmap(3))),2) .eq. 1) then
                  blha_idsubproc(nsp,blha_numproc) = &
                    blha_idsubproc(nsp,blha_numproc) + 100
                endif
                if (mod(abs(pdgparton(invmap(5))),2) .eq. 1) then
                  blha_idsubproc(nsp,blha_numproc) = &
                    blha_idsubproc(nsp,blha_numproc) + 10
                endif
                if (abs(pdgparton(invmap(1))).eq. &
                    abs(pdgparton(invmap(5)))) then ! upper-interf
                  blha_idsubproc(nsp,blha_numproc) = &
                    blha_idsubproc(nsp,blha_numproc) + 1
                endif
                if (abs(pdgparton(invmap(3))).eq. &
                    abs(pdgparton(invmap(5)))) then ! lower-interf
                  blha_idsubproc(nsp,blha_numproc) = &
                    blha_idsubproc(nsp,blha_numproc) + 2
                endif
                if (pdgparton(invmap(1)).eq. &
                    pdgparton(invmap(3))) then ! non-existent-interf
                  blha_idsubproc(nsp,blha_numproc) = &
                    blha_idsubproc(nsp,blha_numproc) + 4
                endif
              else !ZZW
                blha_idsubproc(nsp,blha_numproc) = 0
              endif
            else ! lower-ZW
              if (abs(pdgparton(invmap(5))).eq. &
                  abs(pdgparton(invmap(6)))) then ! ZWZ
                blha_idsubproc(nsp,blha_numproc) = 0
              else ! ZWW
                call ID_VBF6_6qW(nsp,pdgparton,invmap)
              endif
            endif
          else ! upper-W
            if (abs(pdgparton(invmap(3))).eq. &
                abs(pdgparton(invmap(4)))) then ! lower-WZ
              if (abs(pdgparton(invmap(5))).eq. &
                  abs(pdgparton(invmap(6)))) then ! WZZ
                blha_idsubproc(nsp,blha_numproc) = 0
              else ! WZW
                call ID_VBF6_6qW(nsp,pdgparton,invmap)
              endif
            else ! lower-WW
              if (abs(pdgparton(invmap(5))).eq. &
                  abs(pdgparton(invmap(6)))) then ! WWZ
                call ID_VBF6_6qW(nsp,pdgparton,invmap)
              else ! WWW
                blha_idsubproc(nsp,blha_numproc) = 0
              endif
            endif
          endif

          if ( blha_idsubproc(nsp,blha_numproc) .lt. 200000 ) then
            ! not a valid process
            nsp = nsp-1
!          else
!            ! check if we've seen this process already -> ignore as
!            ! this is an already included interference term
!            do j=1,nsp-1
!              if (blha_idsubproc(j,blha_numproc) .eq. &
!                  blha_idsubproc(nsp,blha_numproc) ) then
!                nsp = nsp-1
!                exit
!              endif
!            enddo
          endif

        enddo

      endif

      blha_numsubproc(blha_numproc) = nsp

      return
      end
      
!*************************************************************************  

      SUBROUTINE ID_VBF6_6qW(nsp,pdgparton,invmap)
!*************************************************************************
!     Get the ID for the 6 quark with W exchange process
!*************************************************************************
      implicit none

      integer, intent(in) :: nsp,pdgparton(6),invmap(6)

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

! local variables
      integer i

! external 
      logical samegen
      external samegen

      ! W exchange process
      blha_idsubproc(nsp,blha_numproc) = &
        blha_idsubproc(nsp,blha_numproc) + 10000

      if ( ( pdgparton(1) .gt. 0 .and. &
             mod(abs(pdgparton(1)),2) .eq. 1 ) .or. &
           ( pdgparton(1) .lt. 0 .and. &
             mod(abs(pdgparton(1)),2) .eq. 0 ) ) then
        blha_idsubproc(nsp,blha_numproc) = &
          blha_idsubproc(nsp,blha_numproc) + 1000
      endif
      if ( ( pdgparton(2) .gt. 0 .and. &
             mod(abs(pdgparton(2)),2) .eq. 1 ) .or. &
           ( pdgparton(2) .lt. 0 .and. &
             mod(abs(pdgparton(2)),2) .eq. 0 ) ) then
        blha_idsubproc(nsp,blha_numproc) = &
          blha_idsubproc(nsp,blha_numproc) + 100
      endif
      if (mod(abs(pdgparton(invmap(5))),2) .eq. 1) then
        blha_idsubproc(nsp,blha_numproc) = &
          blha_idsubproc(nsp,blha_numproc) + 10
      endif

      if (pdgparton(invmap(1)) .eq. pdgparton(invmap(2)) ) then ! only col1/col3
        if ( ( pdgparton(1) .gt. 0 .and. &
                 mod(abs(pdgparton(1)),2) .eq. &
                   mod(abs(pdgparton(invmap(5))),2) ) .or. &
             ( pdgparton(1) .lt. 0 .and. &
                 mod(abs(pdgparton(1)),2) .ne. &
                   mod(abs(pdgparton(invmap(5))),2) ) ) then
          if ( samegen(pdgparton(invmap(1)), &
                       pdgparton(invmap(5))) ) then ! upper-interf
            ! interference part of another contribution - clear
            blha_idsubproc(nsp,blha_numproc) = 0
          else
            ! col3
            blha_idsubproc(nsp,blha_numproc) = &
              blha_idsubproc(nsp,blha_numproc) + 5
            if ( pdgparton(1) .lt. 0 ) then
              ! correct flavour of 56 line if 1 is anti-quark
              if (mod(abs(pdgparton(invmap(5))),2) .eq. 1) then
                blha_idsubproc(nsp,blha_numproc) = &
                  blha_idsubproc(nsp,blha_numproc) - 10
              else
                blha_idsubproc(nsp,blha_numproc) = &
                  blha_idsubproc(nsp,blha_numproc) + 10
              endif
              ! correct momentum mapping as this contribution comes from a
              ! flipped diagram: 1 <-> 5
              do i=1,6
                if (blha_physdiagmap(i,nsp,blha_numproc) .eq. 1) then
                  blha_physdiagmap(i,nsp,blha_numproc) = 5
                  blha_physdiagmap(5,nsp,blha_numproc) = 1
                  exit
                endif
              enddo
            else
              ! correct momentum mapping as this contribution comes from a
              ! flipped diagram: 2 <-> 6
              do i=1,6
                if (blha_physdiagmap(i,nsp,blha_numproc) .eq. 2) then
                  blha_physdiagmap(i,nsp,blha_numproc) = 6
                  blha_physdiagmap(6,nsp,blha_numproc) = 2
                  exit
                endif
              enddo
            endif
          endif
        else
          ! col1
          blha_idsubproc(nsp,blha_numproc) = &
            blha_idsubproc(nsp,blha_numproc) + 1
          if ( samegen(pdgparton(invmap(1)), &
                       pdgparton(invmap(5))) ) then ! upper-interf
            blha_idsubproc(nsp,blha_numproc) = &
              blha_idsubproc(nsp,blha_numproc) + 2000
          endif
        endif
      else if (pdgparton(invmap(3)) .eq. pdgparton(invmap(4)) ) then ! only col2/col4
        if ( ( pdgparton(2) .gt. 0 .and. &
                 mod(abs(pdgparton(2)),2) .eq. &
                   mod(abs(pdgparton(invmap(5))),2) ) .or. &
             ( pdgparton(2) .lt. 0 .and. &
                 mod(abs(pdgparton(2)),2) .ne. &
                   mod(abs(pdgparton(invmap(5))),2) ) ) then
          if ( samegen(pdgparton(invmap(3)), &
                       pdgparton(invmap(5))) ) then ! lower-interf
            ! interference part of another contribution - clear
            blha_idsubproc(nsp,blha_numproc) = 0
          else
            ! col4
            blha_idsubproc(nsp,blha_numproc) = &
              blha_idsubproc(nsp,blha_numproc) + 6
            ! correct flavour of 56 line if 2 is anti-quark
            if ( pdgparton(2) .lt. 0 ) then
              if (mod(abs(pdgparton(invmap(5))),2) .eq. 1) then
                blha_idsubproc(nsp,blha_numproc) = &
                  blha_idsubproc(nsp,blha_numproc) - 10
              else
                blha_idsubproc(nsp,blha_numproc) = &
                  blha_idsubproc(nsp,blha_numproc) + 10
              endif
              ! correct momentum mapping as this contribution comes from a
              ! flipped diagram: 3 <-> 5
              do i=1,6
                if (blha_physdiagmap(i,nsp,blha_numproc) .eq. 3) then
                  blha_physdiagmap(i,nsp,blha_numproc) = 5
                  blha_physdiagmap(5,nsp,blha_numproc) = 3
                  exit
                endif
              enddo
            else
              ! correct momentum mapping as this contribution comes from a
              ! flipped diagram: 4 <-> 6
              do i=1,6
                if (blha_physdiagmap(i,nsp,blha_numproc) .eq. 4) then
                  blha_physdiagmap(i,nsp,blha_numproc) = 6
                  blha_physdiagmap(6,nsp,blha_numproc) = 4
                  exit
                endif
              enddo
            endif
          endif
        else
          ! col2
          blha_idsubproc(nsp,blha_numproc) = &
            blha_idsubproc(nsp,blha_numproc) + 2
          if ( samegen(pdgparton(invmap(3)), &
                       pdgparton(invmap(5))) ) then ! lower-interf
            blha_idsubproc(nsp,blha_numproc) = &
              blha_idsubproc(nsp,blha_numproc) + 200
          endif
        endif
      else
        ! col 1 and 2
          blha_idsubproc(nsp,blha_numproc) = &
            blha_idsubproc(nsp,blha_numproc) + 3
        if ( samegen(pdgparton(invmap(1)), &
                     pdgparton(invmap(5))) ) then ! upper-interf
          blha_idsubproc(nsp,blha_numproc) = &
            blha_idsubproc(nsp,blha_numproc) + 2000
        endif
        if ( samegen(pdgparton(invmap(3)), &
                     pdgparton(invmap(5))) ) then ! lower-interf
          blha_idsubproc(nsp,blha_numproc) = &
            blha_idsubproc(nsp,blha_numproc) + 200
        endif
      endif

      ! check for multiplicity factor
      if ( pdgparton(blha_particlemap(3,nsp,blha_numproc)) .eq. &
             pdgparton(blha_particlemap(4,nsp,blha_numproc)) .or. &
           pdgparton(blha_particlemap(3,nsp,blha_numproc)) .eq. &
             pdgparton(blha_particlemap(5,nsp,blha_numproc)) .or. &
           pdgparton(blha_particlemap(3,nsp,blha_numproc)) .eq. &
             pdgparton(blha_particlemap(6,nsp,blha_numproc)) ) then
        blha_idsubproc(nsp,blha_numproc) = &
          blha_idsubproc(nsp,blha_numproc) + 4000
      endif
      if ( pdgparton(blha_particlemap(4,nsp,blha_numproc)) .eq. &
             pdgparton(blha_particlemap(3,nsp,blha_numproc)) .or. &
           pdgparton(blha_particlemap(4,nsp,blha_numproc)) .eq. &
             pdgparton(blha_particlemap(5,nsp,blha_numproc)) .or. &
           pdgparton(blha_particlemap(4,nsp,blha_numproc)) .eq. &
             pdgparton(blha_particlemap(6,nsp,blha_numproc)) ) then
        blha_idsubproc(nsp,blha_numproc) = &
          blha_idsubproc(nsp,blha_numproc) + 400
      endif

      return
      end

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
      integer i,j,k,l,nsp
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

      SUBROUTINE MomMapping_QCD5(pdgparton, mapparton)
!*************************************************************************
!     Momentum mapping for QCD+3j
!*************************************************************************
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

      integer pdgparton(5)
      integer mapparton(5)

! local variables
      integer i,j,k,l,nsp
      integer npairs, listpairs(5,2), nother
      integer tmp
! inverse function to listpairs(:,i), BLHA->diag map
      integer invlistpairs(5)
      integer invmap(5)
! external functions
      integer findpairpartner
      external findpairpartner

      call makefermionpairs(5,pdgparton,.true.,npairs,listpairs,nother)

      nsp=1

! inverse listpairs
      do j=1,5
        invlistpairs(listpairs(j,1)) = j
      enddo

      if (nother.eq.3) then 
! 2q 3g case

! particlemap and fsign
        do j=1,5
          blha_particlemap(j,nsp,blha_numproc) = mapparton(j)
        enddo
! q -> 1
        blha_physdiagmap(listpairs(1,1),nsp,blha_numproc) = 1
        blha_fsign(1,nsp,blha_numproc) = int(sign(1.,2.5-listpairs(1,1)))
! qbar -> 3
        blha_physdiagmap(listpairs(2,1),nsp,blha_numproc) = 3
        blha_fsign(3,nsp,blha_numproc) = int(sign(1.,listpairs(2,1)-2.5))
! g -> 2,4,5
        blha_physdiagmap(listpairs(3,1),nsp,blha_numproc) = 2
        blha_fsign(2,nsp,blha_numproc) = int(sign(1.,listpairs(3,1)-2.5))
        blha_physdiagmap(listpairs(4,1),nsp,blha_numproc) = 4
        blha_fsign(4,nsp,blha_numproc) = int(sign(1.,listpairs(4,1)-2.5))
        blha_physdiagmap(listpairs(5,1),nsp,blha_numproc) = 5
        blha_fsign(5,nsp,blha_numproc) = int(sign(1.,listpairs(5,1)-2.5))

        blha_idsubproc(nsp,blha_numproc) = 100

      else if (nother.eq.1) then 
! 4q1g case

! particlemap and fsign
        do j=1,5
          blha_particlemap(j,nsp,blha_numproc) = mapparton(j)
        enddo
        blha_physdiagmap(listpairs(1,1),nsp,blha_numproc) = 1
        blha_fsign(1,nsp,blha_numproc) = sign(1,pdgparton(listpairs(1,1)))
        blha_physdiagmap(listpairs(3,1),nsp,blha_numproc) = 2
        blha_fsign(2,nsp,blha_numproc) = sign(1,pdgparton(listpairs(3,1)))
        blha_physdiagmap(listpairs(2,1),nsp,blha_numproc) = 3
        blha_fsign(3,nsp,blha_numproc) = sign(1,pdgparton(listpairs(2,1)))
        blha_physdiagmap(listpairs(4,1),nsp,blha_numproc) = 4
        blha_fsign(4,nsp,blha_numproc) = sign(1,pdgparton(listpairs(4,1)))
        blha_physdiagmap(listpairs(5,1),nsp,blha_numproc) = 5
        blha_fsign(5,nsp,blha_numproc) = int(sign(1.,listpairs(5,1)-2.5))

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
            abs(pdgparton(invmap(3)))) then ! upper-Z
          if (abs(pdgparton(invmap(2))).eq. &
              abs(pdgparton(invmap(4)))) then ! NC-ZZ
            if (mod(abs(pdgparton(invmap(1))),2) .eq. 0) then
              if (mod(abs(pdgparton(invmap(2))),2) .eq. 0) then
                blha_idsubproc(nsp,blha_numproc) = 0 ! XXX - uucc
              else
                blha_idsubproc(nsp,blha_numproc) = 0 ! XXX - uuss
              endif
            else
              if (mod(abs(pdgparton(invmap(2))),2) .eq. 0) then
                blha_idsubproc(nsp,blha_numproc) = 0 ! XXX - ddcc
              else
                blha_idsubproc(nsp,blha_numproc) = 0 ! XXX - ddss
              endif
            endif
          else ! CC-ZW
! W line in W/WZ/WA must be 13
! W on 24 line -> exchange
            blha_physdiagmap(listpairs(3,1),nsp,blha_numproc) = 1
            blha_physdiagmap(listpairs(1,1),nsp,blha_numproc) = 2
            blha_physdiagmap(listpairs(4,1),nsp,blha_numproc) = 3
            blha_physdiagmap(listpairs(2,1),nsp,blha_numproc) = 4
            tmp = blha_fsign(1,nsp,blha_numproc) 
            blha_fsign(1,nsp,blha_numproc) = blha_fsign(2,nsp,blha_numproc)
            blha_fsign(2,nsp,blha_numproc) = tmp
            tmp = blha_fsign(3,nsp,blha_numproc) 
            blha_fsign(3,nsp,blha_numproc) = blha_fsign(4,nsp,blha_numproc)
            blha_fsign(4,nsp,blha_numproc) = tmp
            if (mod(abs(pdgparton(invmap(1))),2) .eq. 0) then
              blha_idsubproc(nsp,blha_numproc) = 232 - npairs ! uucs/uusc
            else
              blha_idsubproc(nsp,blha_numproc) = 242 - npairs ! ddcs/ddsc
            endif
          endif
        else ! upper-W
          if (abs(pdgparton(invmap(2))).eq. &
              abs(pdgparton(invmap(4)))) then ! CC-WZ
            if (mod(abs(pdgparton(invmap(2))),2) .eq. 0) then
              blha_idsubproc(nsp,blha_numproc) = 232 - npairs ! udcc/ducc
            else
              blha_idsubproc(nsp,blha_numproc) = 242 - npairs ! udss/duss
            endif
          else ! NC-WW
            if (mod(abs(pdgparton(invmap(1))),2) .eq. 0) then
              blha_idsubproc(nsp,blha_numproc) = 0 ! XXX - udsc
            else
              blha_idsubproc(nsp,blha_numproc) = 0 ! XXX - ducs
            endif
          endif
        endif

      endif

      blha_numsubproc(blha_numproc) = nsp

      return
      end

!*************************************************************************  
      SUBROUTINE MomMapping_QCD(pdgparton, mapparton, np)
!*************************************************************************
!     Momentum mapping for QCD+0j, QCD+1j
!     corresponds to MomMapping_QCD2 and MomMapping_QCD3         
!*************************************************************************
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

      integer, intent(in) :: np ! number of partons
      integer, dimension(np), intent(in) :: pdgparton
      integer, dimension(np), intent(in) :: mapparton

! local variables
      integer i,j,k,l,nsp
      integer npairs, listpairs(np,2), nother
      integer tmp
! inverse function to listpairs(:,i), BLHA->diag map
      integer invlistpairs(np) 
      integer invmap(np)
! external functions
      integer findpairpartner
      external findpairpartner

      call makefermionpairs(np,pdgparton,.true.,npairs,listpairs,nother)

      nsp=1

! inverse listpairs
      do j=1,np
        invlistpairs(listpairs(j,1)) = j
      enddo

      if (nother > np-2) then
         print *, "wrong nother=", nother, ' for np=', np
      endif

! particlemap and fsign
      do j=1,np
        blha_particlemap(j,nsp,blha_numproc) = mapparton(j)
      enddo

! q -> 1
      blha_physdiagmap(listpairs(1,1),nsp,blha_numproc) = 1
      blha_fsign(1,nsp,blha_numproc) = sign(1,pdgparton(listpairs(1,1)))
! qbar -> 2
      blha_physdiagmap(listpairs(2,1),nsp,blha_numproc) = 2
      blha_fsign(2,nsp,blha_numproc) = sign(1,pdgparton(listpairs(2,1)))

      if (np == 3 .and. nother == 1) then
         blha_physdiagmap(listpairs(3,1),nsp,blha_numproc) = 3
         blha_fsign(3,nsp,blha_numproc) = int(sign(1.,listpairs(3,1)-2.5))
      endif
      blha_idsubproc(nsp,blha_numproc) = 1

! create inverse map
      do j=1,np
        do k=1,np
          if (mapparton(k) .eq. blha_particlemap(j,nsp,blha_numproc)) then
            invmap(blha_physdiagmap(j,nsp,blha_numproc)) = k
          endif
        enddo
      enddo

! choose right subproc
       ! TODO: add here for ZZ, Z, ...
        ! if (abs(pdgparton(invmap(1))).eq.
       ! &      abs(pdgparton(invmap(3)))) then ! upper-Z
        !   if (abs(pdgparton(invmap(2))).eq.
       ! &        abs(pdgparton(invmap(4)))) then ! NC-ZZ
        !     if (mod(abs(pdgparton(invmap(1))),2) .eq. 0) then
        !       if (mod(abs(pdgparton(invmap(2))),2) .eq. 0) then
        !         blha_idsubproc(nsp,blha_numproc) = 0 ! XXX - uucc
        !       else
        !         blha_idsubproc(nsp,blha_numproc) = 0 ! XXX - uuss
        !       endif
        !     else
        !       if (mod(abs(pdgparton(invmap(2))),2) .eq. 0) then
        !         blha_idsubproc(nsp,blha_numproc) = 0 ! XXX - ddcc
        !       else
        !         blha_idsubproc(nsp,blha_numproc) = 0 ! XXX - ddss
        !       endif
        !     endif
        !   else ! CC-ZW
! c W line in W/WZ/WA must be 13
! c W on 24 line -> exchange
        !     blha_physdiagmap(listpairs(3,1),nsp,blha_numproc) = 1
        !     blha_physdiagmap(listpairs(1,1),nsp,blha_numproc) = 2
        !     blha_physdiagmap(listpairs(4,1),nsp,blha_numproc) = 3
        !     blha_physdiagmap(listpairs(2,1),nsp,blha_numproc) = 4
        !     tmp = blha_fsign(1,nsp,blha_numproc) 
        !     blha_fsign(1,nsp,blha_numproc) = blha_fsign(2,nsp,blha_numproc)
        !     blha_fsign(2,nsp,blha_numproc) = tmp
        !     tmp = blha_fsign(3,nsp,blha_numproc) 
        !     blha_fsign(3,nsp,blha_numproc) = blha_fsign(4,nsp,blha_numproc)
        !     blha_fsign(4,nsp,blha_numproc) = tmp
        !     if (mod(abs(pdgparton(invmap(1))),2) .eq. 0) then
        !       blha_idsubproc(nsp,blha_numproc) = 232 - npairs ! uucs/uusc
        !     else
        !       blha_idsubproc(nsp,blha_numproc) = 242 - npairs ! ddcs/ddsc
        !     endif
        !   endif
        ! else ! upper-W
        !   if (abs(pdgparton(invmap(2))).eq.
        ! &        abs(pdgparton(invmap(4)))) then ! CC-WZ
        !     if (mod(abs(pdgparton(invmap(2))),2) .eq. 0) then
        !       blha_idsubproc(nsp,blha_numproc) = 232 - npairs ! udcc/ducc
        !     else
        !       blha_idsubproc(nsp,blha_numproc) = 242 - npairs ! udss/duss
        !     endif
        !   else ! NC-WW
        !     if (mod(abs(pdgparton(invmap(1))),2) .eq. 0) then
        !       blha_idsubproc(nsp,blha_numproc) = 0 ! XXX - udsc
        !     else
        !       blha_idsubproc(nsp,blha_numproc) = 0 ! XXX - ducs
        !     endif
        !   endif
        ! endif

      blha_numsubproc(blha_numproc) = nsp

      return
      end


      SUBROUTINE MomMapping_QCD4(pdgparton, mapparton)
!*************************************************************************
!     Momentum mapping for QCD+2j
!*************************************************************************
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

      integer pdgparton(4)
      integer mapparton(4)

! local variables
      integer i,j,k,l,nsp
      integer npairs, listpairs(4,2), nother
      integer tmp
! inverse function to listpairs(:,i), BLHA->diag map
      integer invlistpairs(4)
      integer invmap(4)
! external functions
      integer findpairpartner
      external findpairpartner

      call makefermionpairs(4,pdgparton,.true.,npairs,listpairs,nother)

      nsp=1

! inverse listpairs
      do j=1,4
        invlistpairs(listpairs(j,1)) = j
      enddo

      if (nother.eq.2) then 
! 2q 2g case

! particlemap and fsign
        do j=1,4
          blha_particlemap(j,nsp,blha_numproc) = mapparton(j)
        enddo
! q -> 1
        blha_physdiagmap(listpairs(1,1),nsp,blha_numproc) = 1
        blha_fsign(1,nsp,blha_numproc) = int(sign(1.,2.5-listpairs(1,1))) 
        ! equivalent: sign(1,pdgparton(listpairs(1,1)))
! qbar -> 3
        blha_physdiagmap(listpairs(2,1),nsp,blha_numproc) = 3
        blha_fsign(3,nsp,blha_numproc) = int(sign(1.,listpairs(2,1)-2.5))
! g -> 2,4
        blha_physdiagmap(listpairs(3,1),nsp,blha_numproc) = 2
        blha_fsign(2,nsp,blha_numproc) = int(sign(1.,listpairs(3,1)-2.5))
        blha_physdiagmap(listpairs(4,1),nsp,blha_numproc) = 4
        blha_fsign(4,nsp,blha_numproc) = int(sign(1.,listpairs(4,1)-2.5))

        blha_idsubproc(nsp,blha_numproc) = 100

      else if (nother.eq.0) then 
! 4q case

! particlemap and fsign
        do j=1,4
          blha_particlemap(j,nsp,blha_numproc) = mapparton(j)
        enddo
        blha_physdiagmap(listpairs(1,1),nsp,blha_numproc) = 1
        blha_fsign(1,nsp,blha_numproc) = sign(1,pdgparton(listpairs(1,1)))
        blha_physdiagmap(listpairs(3,1),nsp,blha_numproc) = 2
        blha_fsign(2,nsp,blha_numproc) = sign(1,pdgparton(listpairs(3,1)))
        blha_physdiagmap(listpairs(2,1),nsp,blha_numproc) = 3
        blha_fsign(3,nsp,blha_numproc) = sign(1,pdgparton(listpairs(2,1)))
        blha_physdiagmap(listpairs(4,1),nsp,blha_numproc) = 4
        blha_fsign(4,nsp,blha_numproc) = sign(1,pdgparton(listpairs(4,1)))

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
            abs(pdgparton(invmap(3)))) then ! upper-Z
          if (abs(pdgparton(invmap(2))).eq. &
              abs(pdgparton(invmap(4)))) then ! NC-ZZ
            if (mod(abs(pdgparton(invmap(1))),2) .eq. 0) then
              if (mod(abs(pdgparton(invmap(2))),2) .eq. 0) then
                blha_idsubproc(nsp,blha_numproc) = 0 ! XXX - uucc
              else
                blha_idsubproc(nsp,blha_numproc) = 0 ! XXX - uuss
              endif
            else
              if (mod(abs(pdgparton(invmap(2))),2) .eq. 0) then
                blha_idsubproc(nsp,blha_numproc) = 0 ! XXX - ddcc
              else
                blha_idsubproc(nsp,blha_numproc) = 0 ! XXX - ddss
              endif
            endif
          else ! CC-ZW
! W line in W/WZ/WA must be 13
! W on 24 line -> exchange
            blha_physdiagmap(listpairs(3,1),nsp,blha_numproc) = 1
            blha_physdiagmap(listpairs(1,1),nsp,blha_numproc) = 2
            blha_physdiagmap(listpairs(4,1),nsp,blha_numproc) = 3
            blha_physdiagmap(listpairs(2,1),nsp,blha_numproc) = 4
            tmp = blha_fsign(1,nsp,blha_numproc) 
            blha_fsign(1,nsp,blha_numproc) = blha_fsign(2,nsp,blha_numproc)
            blha_fsign(2,nsp,blha_numproc) = tmp
            tmp = blha_fsign(3,nsp,blha_numproc) 
            blha_fsign(3,nsp,blha_numproc) = blha_fsign(4,nsp,blha_numproc)
            blha_fsign(4,nsp,blha_numproc) = tmp
            if (mod(abs(pdgparton(invmap(1))),2) .eq. 0) then
              blha_idsubproc(nsp,blha_numproc) = 232 - npairs ! uucs/uusc
            else
              blha_idsubproc(nsp,blha_numproc) = 242 - npairs ! ddcs/ddsc
            endif
          endif
        else ! upper-W
          if (abs(pdgparton(invmap(2))).eq. &
              abs(pdgparton(invmap(4)))) then ! CC-WZ
            if (mod(abs(pdgparton(invmap(2))),2) .eq. 0) then
              blha_idsubproc(nsp,blha_numproc) = 232 - npairs ! udcc/ducc
            else
              blha_idsubproc(nsp,blha_numproc) = 242 - npairs ! udss/duss
            endif
          else ! NC-WW
            if (mod(abs(pdgparton(invmap(1))),2) .eq. 0) then
              blha_idsubproc(nsp,blha_numproc) = 0 ! XXX - udsc
            else
              blha_idsubproc(nsp,blha_numproc) = 0 ! XXX - ducs
            endif
          endif
        endif

      endif !qqgg

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

      ! these ids are used to generate a process id
      ! their value is only used in this function
      ! the value is arbitray as long as the resulting ID is unique
      integer IDwp, IDwm, IDzl, IDzn, IDh, IDgamma
      parameter (IDwp   =100000)
      parameter (IDwm   =10000)
      parameter (IDzl   =1000)
      parameter (IDzn   =100)
      parameter (IDh    =10)
      parameter (IDgamma=1)

      nsp = blha_numsubproc(blha_numproc)
      if (blha_numsubproc(blha_numproc) .le. 0) return

      call makefermionpairs(nelweak,pdgelweak,.false., &
                            npairs,listpairs,nother)

      lVBFlepdec = (oQCD .eq. nparton-4) .and. (oQED .eq. nelweak+2 )
      lVBFhaddec = (oQCD .eq. nparton-6) .and. (oQED .eq. nelweak+4 )
      lQCDlepdec = (oQCD .eq. nparton-2) .and. (oQED .eq. nelweak )

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
        ! create id to map out processes
        nid = nwp   *IDwp    &
             +nwm   *IDwm    &
             +nzl   *IDzl    &
             +nzn   *IDzn    &
             +nh    *IDh     &
             +ngamma*IDgamma

        blha_multsubproc((i-1)*nsp+1,blha_numproc) = nhel

        SELECT CASE (nid)
        CASE(IDh) ! H
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
        CASE(IDwp) ! W+
          if (lVBFlepdec) then
            blha_procsubproc((i-1)*nsp+1,blha_numproc) = Wpjj
            if (nparton.eq.4) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
            else if (nparton.eq.5) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .true.
            endif
          else if (lQCDlepdec) then
            if (nparton.eq.4) then
              blha_procsubproc((i-1)*nsp+1,blha_numproc) = QCDWPjj
              blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
            else if (nparton.eq.5) then
              blha_procsubproc((i-1)*nsp+1,blha_numproc) = QCDWPjj
              blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .true.
            endif
          else 
            blha_numsubproc(blha_numproc) = -1
            return
          endif
        CASE(IDWm) ! W-
          if (lVBFlepdec) then
            blha_procsubproc((i-1)*nsp+1,blha_numproc) = Wmjj
            if (nparton.eq.4) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
            else if (nparton.eq.5) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .true.
            endif
          else if (lQCDlepdec) then
            if (nparton.eq.4) then
              blha_procsubproc((i-1)*nsp+1,blha_numproc) = QCDWMjj
              blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
            else if (nparton.eq.5) then
              blha_procsubproc((i-1)*nsp+1,blha_numproc) = QCDWMjj
              blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .true.
            endif
          else 
            blha_numsubproc(blha_numproc) = -1
            return
          endif
        CASE(IDzl) ! Z_l
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
        CASE(IDzn) ! Z_nu
          blha_procsubproc((i-1)*nsp+1,blha_numproc) = Zjj_nu
          if ((nparton.eq.4).and.lVBFlepdec) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
          else if ((nparton.eq.5).and.lVBFlepdec) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .true.
          else 
            blha_numsubproc(blha_numproc) = -1
            return
          endif
        CASE(IDgamma) ! gamma
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
        CASE(2*IDh) ! HH
          blha_procsubproc((i-1)*nsp+1,blha_numproc) = HHjj
          if ((nparton.eq.4).and.lVBFlepdec) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
          else if ((nparton.eq.5).and.lVBFlepdec) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .true.
          else 
            blha_numsubproc(blha_numproc) = -1
            return
          endif
        CASE(IDh + IDgamma) ! H gamma
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
        CASE(IDwp + IDwm) ! W+ W-
          blha_procsubproc((i-1)*nsp+1,blha_numproc) = WpWmjj
          if ((nparton.eq.4).and.lVBFlepdec) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
          else if ((nparton.eq.5).and.lVBFlepdec) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .true.
          else 
            blha_numsubproc(blha_numproc) = -1
            return
          endif
        CASE(IDWp + IDzl) ! W+ Z_l
          if (lVBFlepdec) then
            blha_procsubproc((i-1)*nsp+1,blha_numproc) = WpZjj
            if (nparton.eq.4) then
              blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
            else if (nparton.eq.5) then
              blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .true.
            endif
          else if (lQCDlepdec) then
!            if (nparton.eq.2) then
!              blha_procsubproc((i-1)*nsp+1,blha_numproc) = WPZ
!              blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
!           else if (nparton.eq.3) then 
!              blha_procsubproc((i-1)*nsp+1,blha_numproc) = WPZj
!              blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
!           else if (nparton.eq.4) then
            if (nparton.eq.4) then
              blha_procsubproc((i-1)*nsp+1,blha_numproc) = QCDWPZjj
              blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
            else if (nparton.eq.5) then
              blha_procsubproc((i-1)*nsp+1,blha_numproc) = QCDWPZjj
              blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .true.
            endif
          else 
            blha_numsubproc(blha_numproc) = -1
            return
          endif
        CASE(IDwm + IDzl) ! W- Z_l
          if (lVBFlepdec) then
          blha_procsubproc((i-1)*nsp+1,blha_numproc) = WmZjj
            if (nparton.eq.4) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
            else if (nparton.eq.5) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .true.
            endif
          else if (lQCDlepdec) then
            if (nparton.eq.4) then
              blha_procsubproc((i-1)*nsp+1,blha_numproc) = QCDWmZjj
              blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
            else if (nparton.eq.5) then
              blha_procsubproc((i-1)*nsp+1,blha_numproc) = QCDWmZjj
              blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .true.
            endif
          else 
            blha_numsubproc(blha_numproc) = -1
            return
          endif
        CASE(2*IDzl) ! Z_l Z_l
          blha_procsubproc((i-1)*nsp+1,blha_numproc) = ZZjj_ll
          if ((nparton.eq.4).and.lVBFlepdec) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
          else if ((nparton.eq.5).and.lVBFlepdec) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .true.
          else 
            blha_numsubproc(blha_numproc) = -1
            return
          endif
        CASE(IDzl + IDzn) ! Z_l Z_nu
          blha_procsubproc((i-1)*nsp+1,blha_numproc) = ZZjj_lnu
          if ((nparton.eq.4).and.lVBFlepdec) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
          else if ((nparton.eq.5).and.lVBFlepdec) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .true.
          else 
            blha_numsubproc(blha_numproc) = -1
            return
          endif
        CASE(IDwp + IDgamma) ! W+ A
          if (lVBFlepdec) then
          blha_procsubproc((i-1)*nsp+1,blha_numproc) = WpAjj
            if (nparton.eq.4) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
            else if (nparton.eq.5) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .true.
            endif
          else if (lQCDlepdec) then
            if (nparton.eq.4) then
              blha_procsubproc((i-1)*nsp+1,blha_numproc) = QCDWPAjj
              blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
            else if (nparton.eq.5) then
              blha_procsubproc((i-1)*nsp+1,blha_numproc) = QCDWPAjj
              blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .true.
            endif
          else 
            blha_numsubproc(blha_numproc) = -1
            return
          endif
        CASE(IDwm + IDgamma) ! W- A
          if (lVBFlepdec) then
            blha_procsubproc((i-1)*nsp+1,blha_numproc) = WmAjj
            if (nparton.eq.4) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
            else if (nparton.eq.5) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .true.
            endif
          else if (lQCDlepdec) then
            if (nparton.eq.4) then
              blha_procsubproc((i-1)*nsp+1,blha_numproc) = QCDWMAjj
              blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
            else if (nparton.eq.5) then
              blha_procsubproc((i-1)*nsp+1,blha_numproc) = QCDWMAjj
              blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .true.
            endif
          else 
            blha_numsubproc(blha_numproc) = -1
            return
          endif
        CASE(IDzl + IDgamma) ! Z_l A
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
        CASE(IDzn + IDgamma) ! Z_nu A
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
        CASE(2*IDwp) ! W+ W+
          blha_procsubproc((i-1)*nsp+1,blha_numproc) = WpWpjj
          if ((nparton.eq.4).and.lVBFlepdec) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
          else if ((nparton.eq.5).and.lVBFlepdec) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .true.
          else 
            blha_numsubproc(blha_numproc) = -1
            return
          endif
        CASE(2*IDwm) ! W- W-
          blha_procsubproc((i-1)*nsp+1,blha_numproc) = WmWmjj
          if ((nparton.eq.4).and.lVBFlepdec) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
          else if ((nparton.eq.5).and.lVBFlepdec) then
            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .true.
          else 
            blha_numsubproc(blha_numproc) = -1
            return
          endif
!       CASE(IDh + IDwp) !W+ H
!          if (nparton.eq.2) then
!             blha_procsubproc((i-1)*nsp+1,blha_numproc) = WPH
!             blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
!          else if (nparton.eq.3) then 
!             blha_procsubproc((i-1)*nsp+1,blha_numproc) = WPHj
!             blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
!          else if (nparton.eq.4) then
!             blha_procsubproc((i-1)*nsp+1,blha_numproc) = WPHj
!            blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .true.
!          else 
!            blha_procsubproc((i-1)*nsp+1,blha_numproc) = 0
!            return
!          endif
!       CASE(IDh + IDwm) !W- H
!          if (nparton.eq.2) then
!             blha_procsubproc((i-1)*nsp+1,blha_numproc) = WMH
!             blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
!          else if (nparton.eq.3) then 
!             blha_procsubproc((i-1)*nsp+1,blha_numproc) = WMHj
!             blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .false.
!          else if (nparton.eq.4) then
!             blha_procsubproc((i-1)*nsp+1,blha_numproc) = WMHj
!             blha_lojsubproc((i-1)*nsp+1,blha_numproc) = .true.
!          else
!            blha_numsubproc(blha_numproc) = -1
!            return
!          endif
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
          k=0
          l=1
          do while (k.eq.0)                          ! found first match
            if (listpairs(l,i).gt.0) then
              k = pdgelweak(listpairs(l,i))
            else
              k = 0
            endif
            if (finallep(j) .ne. k) then
              l = l+1
              k = 0
            else
              if (mod(l,2).eq.0) then                ! get partner: next or previous
                m = l-1
              else
                m = l+1
              endif
              if (finallep(j+1) .ne. pdgelweak(listpairs(m,i))) then  ! check pdg id of partner
                l = l+1
                k = 0
              endif
            endif
          enddo
          blha_particlemap(nparton+j,(i-1)*nsp+1,blha_numproc) =  &
            mapelweak(listpairs(l,i))
          blha_fsign(nparton+j,(i-1)*nsp+1,blha_numproc)       =  &
            sign(1,pdgelweak(listpairs(l,i)))
          blha_particlemap(nparton+j+1,(i-1)*nsp+1,blha_numproc) = &
            mapelweak(listpairs(m,i))
          blha_fsign(nparton+j+1,(i-1)*nsp+1,blha_numproc)       = &
            sign(1,pdgelweak(listpairs(m,i)))
          if ( abs(pdgelweak(listpairs(l,i))+ &
                   pdgelweak(listpairs(m,i))) .eq. 0) then
            blha_bosons((j+1)/2,(i-1)*nsp+1,blha_numproc) = 23
          else if ( abs(pdgelweak(listpairs(l,i))+ &
                         pdgelweak(listpairs(m,i))) .eq. 1) then
            blha_bosons((j+1)/2,(i-1)*nsp+1,blha_numproc) = 24
          endif
          listpairs(l,i) = 0                         ! used: clear entries
          listpairs(m,i) = 0
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
            blha_physdiagmap(j,(i-1)*nsp+m,blha_numproc) =  &
              blha_physdiagmap(j,m,blha_numproc)
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
      blha_pssubstep(blha_numproc) = nsp
      blha_numbosons(blha_numproc) = (nelweak-nother)/2

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

      integer, intent(in) :: nf !number of fermions
      integer, intent(in) :: listf(nf)
      logical, intent(in) :: lflip
      integer, intent(out) :: npairs,nother
      integer, intent(out) :: listpairs(nf,*)

! local variables
      integer i,j,k
      integer mylistf(nf)
      integer nbra, listbra(nf)
      integer nket, listket(nf)
      integer listother(nf)
      integer nperm(nf)
      integer temppair(2)
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

! reorder bra-ket pairs so that fermion pairs of the same type stay at the same position
      do i=2,npairs
        do j=1,nbra ! loop over fermion pairs in listpairs(:,1)
          k=1
          do while (k.le.nbra) ! loop over fermion pairs in listpairs(:,i)
            if (k.eq.j) then
              k=k+1
              cycle
            endif
            if ( ( (mylistf(listpairs(2*k-1,i)).eq.mylistf(listpairs(2*j-1,1))) .and.    &
                   (mylistf(listpairs(2*k  ,i)).eq.mylistf(listpairs(2*j  ,1))) )        &
                 .and..not.                                            &
                 ( (mylistf(listpairs(2*k-1,i)).eq.mylistf(listpairs(2*k-1,1))) .and.    &
                   (mylistf(listpairs(2*k  ,i)).eq.mylistf(listpairs(2*k  ,1))) )        &
                 .and..not.                                            &
                 ( (mylistf(listpairs(2*j-1,i)).eq.mylistf(listpairs(2*j-1,1))) .and.    &
                   (mylistf(listpairs(2*j  ,i)).eq.mylistf(listpairs(2*j  ,1))) )        &
               ) then
               ! different-index pair matches but same-index one doesn't -> swap
              temppair(1) = listpairs(2*k-1,i)
              temppair(2) = listpairs(2*k  ,i)
              listpairs(2*k-1,i) = listpairs(2*j-1,i)
              listpairs(2*k  ,i) = listpairs(2*j  ,i)
              listpairs(2*j-1,i) = temppair(1)
              listpairs(2*j  ,i) = temppair(2)
              k=0 ! restart
            endif
            k=k+1
          enddo
        enddo
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

