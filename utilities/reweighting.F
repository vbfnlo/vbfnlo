****************************************************************************
*     This file contains routines for the unweighting of events.
****************************************************************************

****************************************************************************
      Subroutine ResetWeights(ps_number,sub_number)
****************************************************************************
*     call this subroutine at the beginning of each iteration
****************************************************************************
      implicit none

      integer ps_number, sub_number

#include "lha.inc"

      ! for multi-PS, except for multiChannelUnweighting: maxweight of subsequent PS should not be smaller
      ! sub_number > 1 for processes where subprocesses are integrated individually (currently WAAj)
      ! for multiChannelUnweighting treat all PS / subprocesses independently
      if ((ps_number.eq.1 .and. sub_number.eq.1) .or. multiChannelUnweighting) maxweight=0d0

      maxweight_PS = 0d0
      sumwgt = 0d0

      nevraw = 0

      estimatedEvents = 0
      nevovr=0
      nevnewmax=0
      nevacc=0
      nevall=0
      prenevacc=0

      return
      end

****************************************************************************
      Subroutine ReweightingInfo(unweighted)
****************************************************************************
      implicit none

      logical unweighted
      integer getEstimatedEvents
      external getEstimatedEvents

#include "lha.inc"

      if (unweighted) then

            if(lha.or.hepmc)then
               if(unweighting)then
                  if (partialUnweight) then
                     if (.not.multiChannelUnweighting) print*,"maximal weight =",maxweight_PS
                     print*,"number of events with weight = 1 after partial unweighting =", nevacc
                     print*,"number of events with weight > 1 after partial unweighting =", nevovr
                  else
                     if (.not.multiChannelUnweighting) print*,"maximal weight =",wtmax
                     print*,"number of unweighted events =",nevacc
                  endif
               else
                  print*,"maximal weight =",wtmax
               endif
            endif
            print*,""

      else

            print*,"number of events after cuts in the last iteration =", nevraw
            print*,"maximal weight in the last iteration =", maxweight_PS
            print*,"average weight in the last iteration =", sumwgt / nevraw
            print*,"ratio largest / average weight =",maxweight_PS / (sumwgt / nevraw)
            if (partialUnweight) then
              if (nevovr.ge.1) then
                print*, "estimated number of fully unweighted events:", int(nevraw / (maxweight_PS / (sumwgt / nevraw)))
              else
                print*, "estimated number of fully unweighted events:", getEstimatedEvents()
              endif
            else
               if (desiredEventCount.eq.0)
     &            print*, "estimated number of unweighted events:", int(nevraw / (maxweight_PS / (sumwgt / nevraw)))
            endif
            print*,""

      endif
      
      return
      end


****************************************************************************
      Subroutine Reweight(evnweight, maxWeightFixed)
****************************************************************************

      implicit none
      logical maxWeightFixed
      real*8 evnweight,RandomNumber,rn
      external RandomNumber

#include "lha.inc"

      nevraw = nevraw + 1                 ! total events that passed cuts so far in this PS
      sumwgt = sumwgt + evnweight         ! sum of eventweights so far in this PS

      if (evnweight.le.0d0) then
         evnweight = 0d0
         return
      endif

      if(evnweight.ge.maxweight_PS) maxweight_PS = evnweight   ! maximal weight in this PS

      if(unweighting)then
         if(evnweight.ge.maxweight) then
            if (maxWeightFixed .and. partialUnweight) then
               ! print*, "event weight > 1 :", evnweight / maxweight
               nevovr    = nevovr + 1             ! estimating number of events with larger weights than used for unweighting
               evnweight = -evnweight             ! in case of option partialUnweight = .true.
               ! estimatedEvents = estimatedEvents + 1            ! don't count those towards the desired event count: They will end up additionally in event file
            else
               estimatedEvents = int(estimatedEvents * ( dble(maxweight) / dble(evnweight))) + 1
               nevnewmax = nevnewmax + 1           ! number of pre-accepted events with new maxweight
               maxweight = evnweight               ! increase maxweight to current maximum
               maxweight_PS = maxweight
            endif
         else
            rn=RandomNumber()
            if(rn*maxweight.lt.evnweight) then ! event is to be kept
               prenevacc = prenevacc+1         ! number of pre-accepted events
               estimatedEvents = estimatedEvents + 1
               evnweight = maxweight
            else
               evnweight = 0d0
            endif
         endif
      else
         if(evnweight.ge.maxweight) maxweight = evnweight
      endif

      return
      end


******************************************************
      integer function getEstimatedEvents()
******************************************************

      implicit none

#include "lha.inc"

      ! return estimated number of unweighted events in the current phase space
      getEstimatedEvents = estimatedEvents

      ! print*, prenevacc, nevnewmax, nevraw, maxweight, nevovr, estimatedEvents

      return
      end


