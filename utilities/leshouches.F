c------------------------------------------------------------------------------
      subroutine lha_file(imode)

      implicit none

#include "global.inc"
#include "lha.inc"
#include "scales.inc"
#include "process.inc"

      double precision rn, RandomNumber
      double precision unwSumXS, probabilityKeepOverweight, rescaleWeight5Sigma
      integer k

      integer imode
      integer i1,istop
      integer pdflib, valueLHA

      character*250 procfile
      character*250 sampleFile
      character*250 sampleFilePartial
      character*250 sampleinfoFile
      character*7 istr 
      character*250 line
      logical lglobalprint
      
      real*8 xwgtup_max_partial
      real*8 contribOverweightSum(maxContribs), preAcceptedNormalEvents(maxContribs)

      character*250 GetInputPath
      external GetInputPath

      COMMON/PDFparameters/pdflib,valueLHA

      if (hepmc) then
          sampleFile = trim(hepmcFile)//".tmp"
          sampleFilePartial = trim(hepmcFile)//"_partial.tmp"
          sampleinfoFile = trim(hepmcFile)//".info"
      endif
      if (lha) then
          sampleFile = trim(lhaFile)//".tmp"
          sampleFilePartial = trim(lhaFile)//"_partial.tmp"
          sampleinfoFile = trim(lhaFile)//".info"
      endif

      if(imode.eq.-1) then    
         open(31,file=sampleFile,form="unformatted")
         open(32,file=sampleinfoFile)
      elseif(imode.eq.0) then
c...finding the maximal weight 
         close(31)
         open(31,file=sampleFile,form="unformatted")
         istop=1
         wtmax=0d0
         do while(istop.eq.1)
            if (multiChannelUnweighting) then
              ! contributionNumber is the "channel index" for multichannel processes
              read(31,end=200) contributionNumber,nup,idprup,xwgtup,scalup,aqedup,aqcdup
            else
              read(31,end=200) nup,idprup,xwgtup,scalup,aqedup,aqcdup
            endif
            do i1=1,nup
               read(31),idup(i1),istup(i1),mothup(1,i1),
     &              mothup(2,i1),icolup(1,i1),icolup(2,i1), 
     &              pup(1,i1),pup(2,i1),pup(3,i1),pup(4,i1),
     &              pup(5,i1),vtimup(i1),spinup(i1)
            enddo   
            if(xwgtup.gt.wtmax) wtmax=xwgtup
            if (multiChannelUnweighting) then
               if(xwgtup.gt.contribMaxWeight(contributionNumber)) contribMaxWeight(contributionNumber) = xwgtup
            endif
         enddo
 200     continue
         close(31)
         ! for multichannel processes with fixed event number: calculate XS fraction of each channel
         if (multiChannelUnweighting) then 
            unwSumXS = 0.
            do k=1,numContribs
              unwSumXS = unwSumXS + contribXS(k)
            enddo
            do k=1,numContribs
              contribXSFraction(k) = contribXS(k) / unwSumXS
              if (k.gt.1) contribXSFraction(k) = contribXSFraction(k) + contribXSFraction(k-1)
              ! print*, k, contribXS(k), contribXSFraction(k), contribMaxWeight(k)
            enddo
         endif
         open(32,file=sampleinfoFile)
c... PDG codes of the two incoming beam particles
         SELECT CASE(procID)
         CASE(WMWMjj, WMhadWMjj)    !W-W-jj is obtained from W+W+jj, flip IDs
            idbmup(1)=-2212*sign1      
            idbmup(2)=-2212*sign2
         CASE DEFAULT
            idbmup(1)=2212*sign1      
            idbmup(2)=2212*sign2
         END SELECT
c...energies in Gev of two incoming beam particles
         ebmup(1)=ecm/2          
         ebmup(2)=ecm/2
c...number of different external processes        
         nprup=1 
c...master switch dictaiting for events weights and cross sections
c...maximum event weight
         if(unweighting)then
            idwtup=3 
            xmaxup(1)=1d0
         else
            idwtup=2
            xmaxup(1)=wtmax
         endif
c...integer identifier of each external process           
         lprup(1)=procID  
c...not t o clash with pdflib
         pdfgup(1)=0
         pdfgup(2)=0
c...CTEQ6L1 LHA (LHAGLUE)                
         select case(pdflib)
         case(0)   ! cteq6ll - hardwired
         pdfsup(1)=10042
         pdfsup(2)=10042
         case(1)   ! lhapdf
         pdfsup(1)=valueLHA !10042
         pdfsup(2)=valueLHA !10042
         case(2)   ! mrst2004qed
         pdfsup(1)=20460
         pdfsup(2)=20460
         case(3)   ! MSTW2008
         pdfsup(1)=21000
         pdfsup(2)=21000

         end select
c...cross section for each external process in pb
         xsecup(1)=xsection/1000
c...statistical error on the cross section in pb
         xerrup(1)=sdev/1000
         write(32,5100)idbmup(1),idbmup(2),ebmup(1),ebmup(2),
     &        pdfgup(1),pdfgup(2),pdfsup(1),pdfsup(2),idwtup,nprup
         write(32,5200)xsecup(1),xerrup(1),xmaxup(1),lprup(1)         
c         close(31)
         close(32)
      elseif(imode.eq.1) then
         open(31,file=sampleFile,form="unformatted")
         open(32,file=sampleinfoFile)
        if (lha) then 
         open(33,file=lhaFile)
         write(33,"(A)") '<LesHouchesEvents version="1.0">'
         write(33,"(A)") "<!--"
         write(33,"(A)") 
     &   "File generated with VBFNLO - parton level MC program (LO)" 
        endif
#ifdef WITH_HEPMC
        if (hepmc) then 
         call HepMCHeader(hepmcFILE);
        endif
#endif
         
         lglobalprint = .false.
c...open procinfo file
         procfile = "procinfo.dat"
         call GetInputFile(55,procfile,*77)
         
         if (LOplusjet) then
            write (istr,"(a3,I4.4)") "###",procID
         else
            write (istr,"(a2,I4.4)") "##",procID
         end if

 50      continue
         read(55,"(a100)",end=99,err=99) line

         if (index(line,"##").eq.1) lglobalprint = .false.

         if (index(line,istr).eq.1) then 
            lglobalprint = .true.
            goto 50
         endif
         if (lglobalprint.and.lha) write(33,"(A)") trim(line)
#ifdef WITH_HEPMC
         if (lglobalprint.and.hepmc) call HepMCHeader2(line);
#endif
         goto 50
 99      continue      
         close(55)
 77      continue

         if (lha) then
           write(33,"(A)") "-->"
           write(33,"(A)") "<init>"
         endif
         read(32,*)idbmup(1),idbmup(2),ebmup(1),ebmup(2),
     &        pdfgup(1),pdfgup(2),pdfsup(1),pdfsup(2),idwtup,nprup
         read(32,*)xsecup(1),xerrup(1),xmaxup(1),lprup(1)         

         ! determine events with weight > 1 to be stored to file
         if (multiChannelUnweighting) then
            do k=1,numContribs
               contribOverweightSum(k)=0
               preAcceptedNormalEvents(k) = 0
            enddo
         else
            preAcceptedNormalEvents(1) = 0
         endif
         if (partialunweight) then
            open(34,file=sampleFilePartial,form="unformatted")
            nevovr = 0
            xwgtup_max_partial = 0d0
            istop=1
            do while(istop.eq.1) 
               if (multiChannelUnweighting) then
                 read(31,end=101) contributionNumber,nup,idprup,xwgtup,scalup,aqedup,aqcdup
               else
                 read(31,end=101) nup,idprup,xwgtup,scalup,aqedup,aqcdup
               endif
               do i1=1,nup
                  read(31),idup(i1),istup(i1),mothup(1,i1),
     &                 mothup(2,i1),icolup(1,i1),icolup(2,i1), 
     &                 pup(1,i1),pup(2,i1),pup(3,i1),pup(4,i1),
     &                 pup(5,i1),vtimup(i1),spinup(i1)
               enddo
               if (xwgtup.lt.0d0) then
                  ! unweighting in processes with more than one subchannel / phasespace and requesting specific number of events
                  if (multiChannelUnweighting) then
                     wtmax = contribMaxWeight(contributionNumber)
                     probabilityKeepOverweight = ( (-xwgtup) / wtmax * contribXS(contributionNumber) * dble(desiredEventCount) 
     &                                 / unwSumXS / dble(contribEstimatedEvents(contributionNumber)) )
                     rescaleWeight5Sigma = ( contribXS(contributionNumber) * dble(desiredEventCount) / unwSumXS 
     &                                 / dble(contribEstimatedEvents(contributionNumber)) )
                  else
                     probabilityKeepOverweight = ( (-xwgtup) / wtmax * dble(desiredEventCount)
     &                                 / dble(contribEstimatedEvents(1)) )
                     rescaleWeight5Sigma = ( dble(desiredEventCount) / dble(contribEstimatedEvents(1)) )
                  endif
                  ! print*, xwgtup, wtmax, desiredEventCount, contribEstimatedEvents(1)
                  ! print*, probabilityKeepOverweight
                  rn=RandomNumber()
                  if( rn .lt. probabilityKeepOverweight ) then ! event is to be kept
                     xwgtup = -xwgtup / wtmax
                     xwgtup = max( xwgtup * rescaleWeight5Sigma, 1d0 )        ! compensate for larger probability of events with
                                                                              ! large weights to end up in event file
                     if (xwgtup.gt.1d0) then                                  ! only treat events specially with final weight > 1
                        nevovr = nevovr + 1
                        if (xwgtup .gt. xwgtup_max_partial) xwgtup_max_partial = xwgtup       ! largest weight for event file header
                        if (multiChannelUnweighting) then     ! store contributions from overweight events to get XS fractions right
                           contribOverweightSum(contributionNumber) = contribOverweightSum(contributionNumber) + xwgtup
                        endif
                     else
                        if (multichannelunweighting) then
                           preAcceptedNormalEvents(contributionNumber) = preAcceptedNormalEvents(contributionNumber) + 1
                        else
                           preAcceptedNormalEvents(1) = preAcceptedNormalEvents(1) + 1
                        endif
                     endif
                     ! write the information of large weights to a temporary event file
                     if (multiChannelUnweighting) then
                        ! contributionNumber is the "channel index" for multichannel processes
                        write(34)contributionNumber,nup,idprup,xwgtup,scalup,aqedup,aqcdup
                     else
                        write(34)nup,idprup,xwgtup,scalup,aqedup,aqcdup
                     endif
                     do i1=1,nup
                        write(34),idup(i1),istup(i1),mothup(1,i1),
     &                  mothup(2,i1),icolup(1,i1),icolup(2,i1), 
     &                  pup(1,i1),pup(2,i1),pup(3,i1),pup(4,i1),
     &                  pup(5,i1),vtimup(i1),spinup(i1)
                     enddo
                  endif
               endif
            enddo
 101        continue
            if (nevovr.gt.0) then
               idwtup=2
               xmaxup(1) = xwgtup_max_partial
            endif
            if (multiChannelUnweighting) then
               do k=2,numContribs
                  contribOverweightSum(k) = contribOverweightSum(k) + contribOverweightSum(k-1)
               enddo
            endif
            close(34)   ! temporary file for events with large weights
            close(31)
            open(31,file=sampleFile,form="unformatted")
         endif   ! partialUnweight

         ! write out header of event files
         if (lha) then
           write(33,5100)idbmup(1),idbmup(2),ebmup(1),ebmup(2),
     &          pdfgup(1),pdfgup(2),pdfsup(1),pdfsup(2),idwtup,nprup
           write(33,5200)xsecup(1),xerrup(1),xmaxup(1),lprup(1)
           write(33,"(A)") "</init>"
         endif
#ifdef WITH_HEPMC
         if (hepmc) then
           call HepMCHeader3(idbmup(1),idbmup(2),ebmup(1),ebmup(2),
     &          pdfgup(1),pdfgup(2),pdfsup(1),pdfsup(2),idwtup,nprup,
     &          xsecup(1),xerrup(1),xmaxup(1),lprup(1));
         endif
#endif

         ! store events with weight > 1 to file
         if (partialunweight) then
            open(34,file=sampleFilePartial,form="unformatted")
            istop=1
            do while(istop.eq.1) 
               if (multiChannelUnweighting) then
                 read(34,end=102) contributionNumber,nup,idprup,xwgtup,scalup,aqedup,aqcdup
               else
                 read(34,end=102) nup,idprup,xwgtup,scalup,aqedup,aqcdup
               endif
               do i1=1,nup
                  read(34),idup(i1),istup(i1),mothup(1,i1),
     &                 mothup(2,i1),icolup(1,i1),icolup(2,i1), 
     &                 pup(1,i1),pup(2,i1),pup(3,i1),pup(4,i1),
     &                 pup(5,i1),vtimup(i1),spinup(i1)
               enddo
               ! write events to final event file
               if (lha) then
                  write(33,"(A)") "<event>"
                  write(33,5300)nup,idprup,xwgtup,scalup,aqedup,aqcdup
                  do i1=1,nup
                     write(33,5400),idup(i1),istup(i1),mothup(1,i1),
     &                    mothup(2,i1),icolup(1,i1),icolup(2,i1), 
     &                    pup(1,i1),pup(2,i1),pup(3,i1),pup(4,i1),
     &                    pup(5,i1),vtimup(i1),spinup(i1)
                  enddo   
                  write(33,"(A)") "</event>"
               endif
#ifdef WITH_HEPMC
               if (hepmc) then
                 call HepMCEvent(scalup,aqcdup,aqedup,xwgtup)
                 do i1=1,nup
                   call HepMCEventP(idup(i1),istup(i1),
     &                 pup(1,i1),pup(2,i1),pup(3,i1),pup(4,i1),
     &                 pup(5,i1))
                 enddo
                 call HepMCEventE()
               endif
#endif
            enddo
 102        continue
            close(34,status="delete")   ! temporary file for events with large weights
         endif   ! partialUnweight

         ! write "normal" events to file
         istop=1
         do while(istop.eq.1) 
            nevall=nevall+1
            if (multiChannelUnweighting) then
              read(31,end=100) contributionNumber,nup,idprup,xwgtup,scalup,aqedup,aqcdup
            else
              read(31,end=100) nup,idprup,xwgtup,scalup,aqedup,aqcdup
            endif
            do i1=1,nup
               read(31),idup(i1),istup(i1),mothup(1,i1),
     &              mothup(2,i1),icolup(1,i1),icolup(2,i1), 
     &              pup(1,i1),pup(2,i1),pup(3,i1),pup(4,i1),
     &              pup(5,i1),vtimup(i1),spinup(i1)

            enddo
            if(unweighting)then
               ! unweighting in processes with more than one subchannel / phasespace and requesting specific number of events
               if (multiChannelUnweighting) then
                  if (preAcceptedNormalEvents(contributionNumber).gt.0) then
                     nevacc = nevacc + preAcceptedNormalEvents(contributionNumber)
                     preAcceptedNormalEvents(contributionNumber) = 0
                  endif
                  wtmax = contribMaxWeight(contributionNumber)
!                  if ( nevacc.ge.desiredEventCount*contribXSFraction(contributionNumber) ) then
                  if ( (nevacc + contribOverweightSum(contributionNumber)) .ge. 
     &                  (desiredEventCount + contribOverweightSum(numContribs)) * contribXSFraction(contributionNumber) ) then
                     ! enough events in the current contribution
                     xwgtup = 0
                  endif
                  if (contributionNumber.gt.1) then
!                     if (nevacc.lt.desiredEventCount*contribXSFraction(contributionNumber-1)) then
                     if ( (nevacc + contribOverweightSum(contributionNumber-1)) .lt. 
     &                  (desiredEventCount + contribOverweightSum(numContribs)) * contribXSFraction(contributionNumber-1) ) then
                        ! too few events in the last contribution
                        istop = 0
                     endif
                  endif
               else
                  if (preAcceptedNormalEvents(1).gt.0) then
                     nevacc = nevacc + preAcceptedNormalEvents(1)
                     preAcceptedNormalEvents(1) = 0
                  endif
                  if (nevacc.ge.desiredEventCount .and. desiredEventCount.gt.0) then      ! requesting specific number of events
                     ! enough events found
                     istop = 0
                     xwgtup = 0
                  endif
               endif
               rn=RandomNumber()
               if(rn*wtmax.lt.xwgtup) then ! event is to be kept
                  nevacc=nevacc+1
                  xwgtup=1d0
                  if (lha) then
                  write(33,"(A)") "<event>"
                  write(33,5300)nup,idprup,xwgtup,scalup,aqedup,aqcdup
                  do i1=1,nup
                     write(33,5400),idup(i1),istup(i1),mothup(1,i1),
     &                    mothup(2,i1),icolup(1,i1),icolup(2,i1), 
     &                    pup(1,i1),pup(2,i1),pup(3,i1),pup(4,i1),
     &                    pup(5,i1),vtimup(i1),spinup(i1)
                  enddo   
                  write(33,"(A)") "</event>"
                  endif
#ifdef WITH_HEPMC
                  if (hepmc) then
                    call HepMCEvent(scalup,aqcdup,aqedup,xwgtup)
                    do i1=1,nup
                      call HepMCEventP(idup(i1),istup(i1),
     &                    pup(1,i1),pup(2,i1),pup(3,i1),pup(4,i1),
     &                    pup(5,i1))
                    enddo
                    call HepMCEventE()
                  endif
#endif
               endif
            else
               if (lha) then
               write(33,"(A)") "<event>"
               write(33,5300)nup,idprup,xwgtup,scalup,aqedup,aqcdup
               do i1=1,nup
                  write(33,5400),idup(i1),istup(i1),mothup(1,i1),
     &                 mothup(2,i1),icolup(1,i1),icolup(2,i1), 
     &                 pup(1,i1),pup(2,i1),pup(3,i1),pup(4,i1),
     &                 pup(5,i1),vtimup(i1),spinup(i1)
               enddo   
               write(33,"(A)") "</event>"               
               endif
#ifdef WITH_HEPMC
               if (hepmc) then
                 call HepMCEvent(scalup,aqcdup,aqedup,xwgtup)
                 do i1=1,nup
                   call HepMCEventP(idup(i1),istup(i1),
     &                 pup(1,i1),pup(2,i1),pup(3,i1),pup(4,i1),
     &                 pup(5,i1))
                 enddo
                 call HepMCEventE()
               endif
#endif
            endif
         enddo
 100     continue
         nevall=nevall-1
         
         ! possible extension for the future: randomize order of events for multichannel PS
         ! instead of writing everything directly to HepMC / LHE, write to temporary file
         ! with random access.
         ! Sort array (i,rand(i)) by second column
         ! Read events sorted(i) from file and write to final event file
         ! (seen in MadEvent)
         ! For:
         ! Unweighted events
         ! Partially unweighted events? (maybe only events with weight 1? )
         ! Weighted events? (massive additional IO)
         ! Include switch in lha.inc / global.inc?
         
         if (lha) then
           write(33,"(A)") "</LesHouchesEvents>"
           close(33)
         endif
#ifdef WITH_HEPMC
         if (hepmc) call HepMCEnd();
#endif
         close(31,status="delete")
         close(32,status="delete")
         
         if (multiChannelUnweighting) then
            if (nevacc.lt.desiredEventCount) then
                print*, ""
                print*, "Could not get the desired number of events for at least one subchannel."
                print*, "Therefore we can't guarantee a correct sampling of the events and the"
                print*, "event file will be removed."
                print*, ""
                nevacc = 0
                if (lha) then
                  OPEN (33, FILE=lhaFile, STATUS="OLD")
                  CLOSE (33, STATUS="DELETE")
                endif
                if (hepmc) then
                  OPEN (33, FILE=hepmcFile, STATUS="OLD")
                  CLOSE (33, STATUS="DELETE")
                endif
            endif
         else
            if (unweighting .and. nevacc.lt.desiredEventCount) then
                print*, ""
                print*, "Could not get the desired number of ", desiredEventCount, "events."
                print*, "Writing out all events that have been obtained so far."
                print*, ""
            endif
         endif
      endif
 5100 format(1p,2i8,2e14.6,6i8)
 5200 format(1p,3e14.6,i6)
 5300 format(1p,2i6,4e14.6)
 5400 format(1p,i8,5i5,5e18.10,e9.1,e9.1)  
      end

c------------------------------------------------------------------------------

      subroutine lha_put(p,v)

      implicit none

#include "global.inc"
#include "lha.inc"
#include "scales.inc"
#include "process.inc"
#include "mssm.inc"

      double precision pinv, RandomNumber

      integer iput,gnlo  
      integer i1,i11,i2
      integer idflav(6),listposition

      double precision p(0:3,max_p,max_kin)
      double precision v(0:3,max_v,max_kin)
      double precision vv(4),pp1(4),pp2(4),ncms1(3),ncms2(3)
      double precision lpp1, lpp2, mtau, pp, rn
      double precision  alfas, xmt, alfa, xmz, xmw, sin2w, xmh, gf, alfas_lo
      common /bkopin/ alfas, xmt, alfa, xmz, xmw, sin2w, xmh, gf, alfas_lo


      mtau = ML

 
      if(doNLO)then
         gnlo=1
      else
         gnlo=0
      endif
c...number of particle entries in the current parton level event 
      nup=lnup(iprocess)
c...identity of the current process
      idprup=procID
c...event weight
      xwgtup=eventweight        !/1000
c...scale Q of event
      scalup=max((sqrt(mufsq(1,1))),(sqrt(mufsq(2,1))))  
c...QED coupling
      aqedup=alfa
c...QCD coupling
      aqcdup=max(als(1,1),als(2,1))
c...number of boson(s)
      if(n_v.eq.1)then
         n_b=1
      elseif(n_v.eq.3)then
         n_b=2
      elseif(n_v.eq.5)then
         n_b=3
      else
         n_b=n_v/2
      endif
c------------------------------------
c...Quarks & gluon
c------------------------------------

      do i1=1,n_p 

c...particle identity code as in PDG               
         idup(i1)=lidup(i1,iprocess)
c...status code of particle (=-1 incoming particle =1 outgoing)
         istup(i1)=listup(i1,iprocess)
c...momenta of partons     
         pup(1,i1)=p(1,i1,1)
         pup(2,i1)=p(2,i1,1)
         pup(3,i1)=p(3,i1,1)
         pup(4,i1)=p(0,i1,1)
         pup(5,i1)=0d0
c...orgin information
         mothup(1,i1)=lmothup(1,i1,iprocess)
         mothup(2,i1)=lmothup(2,i1,iprocess)
c...colour information
         icolup(1,i1)=licolup(1,i1,iprocess)
         icolup(2,i1)=licolup(2,i1,iprocess)
c...secondary vertex information (invariant lifetime)
         vtimup(i1)=0d0
c...spin information
         spinup(i1)=9d0
      enddo
c...particle identity code for vector boson decay products as in PDG 

c------------------------------------
c...Decay products of vector bosons
c------------------------------------

      call fillColorless(idflav,listposition)

      do i1=1,n_v
         i11=n_p+n_b+i1
         
         if(i11.le.nup)then
c...particle identity code
            idup(i11)=lidup(i11,1)
c...status code of particle (=-1 incoming particle =1 outgoing)
            istup(i11)=listup(i11,1)
c...colour information
            icolup(1,i11)=licolup(1,i11,1)
            icolup(2,i11)=licolup(2,i11,1)
c...orgin information
            mothup(1,i11)=lmothup(1,i11,1)
            mothup(2,i11)=lmothup(2,i11,1)
c...life time and spin information
            vtimup(i11)=0d0
            spinup(i11)=helicity(i1)

c...momenta of vector boson(s) decay products
            pup(1,i11)=v(1,i1,1)
            pup(2,i11)=v(2,i1,1)
            pup(3,i11)=v(3,i1,1)
            pup(4,i11)=v(0,i1,1)
            select case(procID)
            case(HHjj)
               spinup(i11)=9d0
               pup(5,i11)=dsqrt(pup(4,i11)**2-pup(3,i11)**2
     &              -pup(2,i11)**2-pup(1,i11)**2)
            case DEFAULT
               pup(5,i11)=0d0
            end select
            if(taumass_switch) then
               if(abs(idflav(i1)).eq.15) pup(5,i11)=mtau
            endif

         endif
      enddo  

c...helicity configurations for some processes (decay of Higgs boson only)
      select case(procID)
      ! AA and bb final states: unpolarized
      case(Hjj_AA, Hjj_bbar, Hjjj_AA,Hjjj_bbar, HAjj_AA, HAjj_bbar,
     -     GFHjj_AA, GFHjj_bbar, HHjj_bbAA)
         do i1=1,n_v
            i11=n_p+n_b+i1
            spinup(i11)=9d0
         enddo
      ! tautau and mumu final states: coinflip
      case(Hjj_mu, Hjj_tau, HAjj_mu, HAjj_tau,
     -        Hjjj_mu, Hjjj_tau,
     -        GFHjj_mu, GFHjj_tau)
         do i1=1,n_v,2
            i11=n_p+n_b+i1
            rn=RandomNumber()
            if (rn.gt.0.5d0) then
               spinup(i11)   =  1
               spinup(i11+1) =  1
            else
               spinup(i11)   = -1
               spinup(i11+1) = -1
            endif
         enddo
      case(HHjj_bbtau) ! special case HHjj --> b bbar tau tau
         do i1=1,n_v
            i11=n_p+n_b+i1
            if(i1.le.2) then
               spinup(i11)=9d0
            else
               if(i1.eq.3) then
                  rn=RandomNumber()
                  if (rn.gt.0.5d0) then
                     spinup(i11)   =  1
                     spinup(i11+1) =  1
                  else
                     spinup(i11)   = -1
                     spinup(i11+1) = -1
                  endif
               endif
            endif
         enddo
      ! WW and ZZ final states: helicity variable is set in subroutine m2s_VVsum
      end select

c...mass effects for tau leptons
      if(taumass_switch)then

         do i1=1,n_v-1,2
            i11=n_p+n_b+i1

            ! assumption: photons at the end of v => decay pairs start at i1=1,3,...
            if ( abs(idflav(i1)).eq.15 .and. abs(idflav(i1+1)).eq.16 .or.           ! W decay
     &           abs(idflav(i1)).eq.16 .and. abs(idflav(i1+1)).eq.15 ) then

               if (mothup(1,i11).ne.mothup(1,i11+1) .or. mothup(2,i11).ne.mothup(2,i11+1)) then
                  print*, "Fatal error in Les Houches event output for the inclusion of tau masses:"
                  print*, "The particle ordering for this process is not as expected."
                  print*, "Please report this problem to the authors."
                  stop
               endif

               if ( abs(idflav(i1)).eq.15 ) then      ! W -> tau + neutrino
                  do i2=1,3
                     pp1(i2)=v(i2,i1,1)   ! tau+ or tau-
                     pp2(i2)=v(i2,i1+1,1) ! nu_tau or nu_tau~
                  enddo
                  pp1(4)=v(0,i1,1)
                  pp2(4)=v(0,i1+1,1)
               else                                   ! W -> neutrino + tau
                  do i2=1,3
                     pp1(i2)=v(i2,i1+1,1) ! tau+ or tau-
                     pp2(i2)=v(i2,i1,1)   ! nu_tau or nu_tau~
                  enddo
                  pp1(4)=v(0,i1+1,1)
                  pp2(4)=v(0,i1,1)
               endif

               do i2=1,4
                  vv(i2)=pp1(i2)+pp2(i2)
               enddo
               pinv=dsqrt(vv(4)**2-vv(1)**2-vv(2)**2-vv(3)**2)

               ! dismiss this event if not enough energy available for inclusion of tau mass
               if (pinv .le. mtau) return
               
               call bostdq(1,vv,pp1,pp1)
               call bostdq(1,vv,pp2,pp2)
               
               lpp1=dsqrt(pp1(1)**2+pp1(2)**2+pp1(3)**2)
               lpp2=dsqrt(pp2(1)**2+pp2(2)**2+pp2(3)**2)
               do i2=1,3
                  ncms1(i2)=pp1(i2)/lpp1
                  ncms2(i2)=pp2(i2)/lpp2
               enddo
               pp1(4)=(pinv**2+mtau**2)/(2d0*pinv)
               pp2(4)=(pinv**2-mtau**2)/(2d0*pinv)
               pp=dsqrt(pp1(4)**2-mtau**2)            
               do i2=1,3
                  pp1(i2)=pp*ncms1(i2)
                  pp2(i2)=pp2(4)*ncms2(i2)
               enddo 

               call bostdq(-1,vv,pp1,pp1)
               call bostdq(-1,vv,pp2,pp2)

               do i2=1,4
                  pup(i2,n_p+n_b+i1)=pp1(i2)
                  pup(i2,n_p+n_b+i1+1)=pp2(i2)
               enddo

               if ( abs(idflav(i1)).eq.15 ) then      ! W -> tau + neutrino
                  do i2=1,4
                     pup(i2,n_p+n_b+i1)   = pp1(i2)   ! tau
                     pup(i2,n_p+n_b+i1+1) = pp2(i2)   ! neutrino
                  enddo
               else                                   ! W -> neutrino + tau
                  do i2=1,4
                     pup(i2,n_p+n_b+i1+1) = pp1(i2)   ! tau
                     pup(i2,n_p+n_b+i1)   = pp2(i2)   ! neutrino
                  enddo
               endif

            elseif ( abs(idflav(i1)).eq.15 .and. abs(idflav(i1+1)).eq.15 ) then     ! Z decay

               if (mothup(1,i11).ne.mothup(1,i11+1) .or. mothup(2,i11).ne.mothup(2,i11+1)) then
                  print*, "Fatal error in Les Houches event output for the inclusion of tau masses:"
                  print*, "The particle ordering for this process is not as expected."
                  print*, "Please report this problem to the authors."
                  stop
               endif

               do i2=1,3
                  pp1(i2)=v(i2,i1,1)
                  pp2(i2)=v(i2,i1+1,1)
               enddo
               pp1(4)=v(0,i1,1)
               pp2(4)=v(0,i1+1,1)

               do i2=1,4
                  vv(i2)=pp1(i2)+pp2(i2)
               enddo
               pinv=dsqrt(vv(4)**2-vv(1)**2-vv(2)**2-vv(3)**2)

               ! dismiss this event if not enough energy available for inclusion of tau mass
               if (pinv .le. 2*mtau) return

               call bostdq(1,vv,pp1,pp1)
               call bostdq(1,vv,pp2,pp2)

               lpp1=dsqrt(pp1(1)**2+pp1(2)**2+pp1(3)**2)
               lpp2=dsqrt(pp2(1)**2+pp2(2)**2+pp2(3)**2)
               do i2=1,3
                  ncms1(i2)=pp1(i2)/lpp1
                  ncms2(i2)=pp2(i2)/lpp2
               enddo
               pp=dsqrt(pp1(4)**2-mtau**2)
               do i2=1,3
                  pp1(i2)=pp*ncms1(i2)
                  pp2(i2)=pp*ncms2(i2)
               enddo

               call bostdq(-1,vv,pp1,pp1)
               call bostdq(-1,vv,pp2,pp2)

               do i2=1,4
                  pup(i2,n_p+n_b+i1)   = pp1(i2)
                  pup(i2,n_p+n_b+i1+1) = pp2(i2)
               enddo

            endif

         enddo

      endif

c------------------------------------
c...Vector bosons
c------------------------------------

c...life time and spin information
      do i1=1,n_b
         i11=n_p+i1

c...particle identity code
         idup(i11)=lidup(i11,1)
c...status code of particle (=-1 incoming particle =1 outgoing)
         istup(i11)=listup(i11,1)
c...colour information
         icolup(1,i11)=licolup(1,i11,1)
         icolup(2,i11)=licolup(2,i11,1)
c...orgin information
         mothup(1,i11)=lmothup(1,i11,1)
         mothup(2,i11)=lmothup(2,i11,1)
         vtimup(i11)=0d0
         spinup(i11)=9d0
      enddo
      iput=0
      do i1=1,n_b*2,2
         iput=iput+1
	 i11=n_p+iput		
c...momenta of vector boson(s)
         pup(1,i11)=v(1,i1,1)+v(1,i1+1,1)
         pup(2,i11)=v(2,i1,1)+v(2,i1+1,1)
         pup(3,i11)=v(3,i1,1)+v(3,i1+1,1)
         pup(4,i11)=v(0,i1,1)+v(0,i1+1,1)
c...invariant mass of vector boson(s)
         if (idup(i11).ne.22) then
            pup(5,i11)=dsqrt(pup(4,i11)**2-pup(3,i11)**2
     &           -pup(2,i11)**2-pup(1,i11)**2)
         else 
            pup(5,i11)=0d0
         endif
      enddo

c...write the information to the temporary event file
      if (multiChannelUnweighting) then
         ! contributionNumber is the "channel index" for multichannel processes
         write(31)contributionNumber,nup,idprup,xwgtup,scalup,aqedup,aqcdup
      else
         write(31)nup,idprup,xwgtup,scalup,aqedup,aqcdup
      endif
      do i1=1,nup
         write(31),idup(i1),istup(i1),mothup(1,i1),
     &        mothup(2,i1),icolup(1,i1),icolup(2,i1), 
     &        pup(1,i1),pup(2,i1),pup(3,i1),pup(4,i1),
     &        pup(5,i1),vtimup(i1),spinup(i1)
      enddo
  
      end
      
c------------------------------------------------------------------------------

      subroutine kf_leptons(idflav)

      implicit none

#include "global.inc"
#include "process.inc"

      integer idflav(6), i

      select case(procID)
      case(Hjj_tau, HAjj_tau, Hjjj_tau, GFHjj_tau) 
         idflav(1)=-15
         idflav(2)= 15  
      case(Hjj_AA, HAjj_AA, Hjjj_AA, AAjj, GFHjj_AA,QCDAAjj) 
         idflav(1)= 22
         idflav(2)= 22   
      case(Hjj_bbar, HAjj_bbar, Hjjj_bbar, GFHjj_bbar)
         idflav(1)= 5
         idflav(2)=-5
      case(Hjj_mu, HAjj_mu, Hjjj_mu, GFHjj_mu)
         idflav(1)=-13
         idflav(2)= 13
      case(Hjj_WW, Hjjj_WW,WPWMjj, Hjj_ZZ_ll, Hjjj_ZZ_ll,ZZjj_ll, 
     1     Hjj_ZZ_lnu,Hjjj_ZZ_lnu, 
     2     ZZjj_lnu,Zjj_l, Zjj_nu, WPjj, WMjj, QCDZjj_l, QCDZjj_nu,
     3     WWZ, ZZWP, ZZWM, ZZZ, WW, GFWW, WpZ, WmZ, ZZ, GFZZ, 
     &     WWj, GFWWj, ZZj, GFZZj,
     &     WWWP, WWWM, WPZjj, WMZjj,
     4     HAjj_WW, HAjj_ZZ_ll, HAjj_ZZ_lnu,
     5     GFHjj_WW, GFHjj_ZZ_ll, GFHjj_ZZ_lnu,
     6     WPWPjj, WMWMjj,
     7     Sp2jj_WW, Sp2jj_ZZ_ll, Sp2jj_ZZ_lnu,
     8     QCDWPjj,QCDWMjj,QCDWPWPjj,QCDWMWMjj,QCDWPZjj,QCDWMZjj,
     9     QCDZZjj_ll,QCDZZjj_lnu)
         do i = 1, N_lepton
            idflav(i) = finallep(i)
         enddo
      ! first boson decays hadronically
      CASE(WPhadWMlep, GFWPhadWMlep, ZZhad, GFZZhad,
     &     WPhadZ, WMhadZ,
     &     WPhadWMZ, ZZhadWP, ZZhadWM, WPhadWMWP, WMhadWPWM, ZZZhad,
     &     WPhadWMjj, ZZhadjj, WPhadZjj, WMhadZjj,
     &     Hjj_WPhadWM, Hjj_ZZhad, WPhadWPjj, WMhadWMjj)
         do i = 1, N_lepton + N_decayquarks
           idflav(i) = Particle_IDs(i)
         enddo
         if (particle_ids(1).eq.93 .or. particle_ids(1).eq.94) then
           if (procid.eq.WPhadZjj .or. procid.eq.ZZhadjj .or. procid.eq.WPhadWPjj .or. procid.eq.WPhadZj) then
             ! switch to antiquark - quark ordering
             idflav(1) = finalquarks_psp(2)
             idflav(2) = finalquarks_psp(1)
           else
             idflav(1) = finalquarks_psp(1)
             idflav(2) = finalquarks_psp(2)
           endif
         endif
      ! second boson decays hadronically
      CASE(WPlepWMhad, GFWPlepWMhad,
     &     WPZhad, WMZhad,
     &     WPWMhadZ, WPWMhadWP, WMWPhadWM,
     &     WPWMhadjj, WPZhadjj, WMZhadjj,
     &     Hjj_WPWMhad)
         do i = 1, N_lepton + N_decayquarks
           idflav(i) = Particle_IDs(i)
         enddo
         if (particle_ids(3).eq.93 .or. particle_ids(3).eq.94) then
           if (procid.eq.WPWMhadjj .or. procid.eq.WPZhadjj .or. procid.eq.WMZhadjj .or. procid.eq.WPZhadj) then
             ! switch to antiquark - quark ordering
             idflav(3) = finalquarks_psp(2)
             idflav(4) = finalquarks_psp(1)
           else
             idflav(3) = finalquarks_psp(1)
             idflav(4) = finalquarks_psp(2)
           endif
         endif
      ! third boson decays hadronically
      CASE(WWZhad, ZZWPhad, ZZWMhad)
         do i = 1, N_lepton + N_decayquarks
           idflav(i) = Particle_IDs(i)
         enddo
         if (particle_ids(5).eq.93 .or. particle_ids(5).eq.94) then
           idflav(5) = finalquarks_psp(1)
           idflav(6) = finalquarks_psp(2)
         endif
      ! first boson in triboson with photon processes decays hadronically
      CASE(WPhadWMA, WPhadZA, WMhadZA, ZZhadA)
         do i = 1, N_lepton + N_decayquarks
           idflav(i) = Particle_IDs(i)
         enddo
         if (particle_ids(1).eq.93 .or. particle_ids(1).eq.94) then
           idflav(1) = finalquarks_psp(1)
           idflav(2) = finalquarks_psp(2)
         endif
         idflav(N_lepton+N_decayquarks+1) = 22
      ! second boson in triboson with photon processes decays hadronically
      CASE(WPWMhadA, WPZhadA, WMZhadA)
         do i = 1, N_lepton + N_decayquarks
           idflav(i) = Particle_IDs(i)
         enddo
         if (particle_ids(3).eq.93 .or. particle_ids(3).eq.94) then
           idflav(3) = finalquarks_psp(1)
           idflav(4) = finalquarks_psp(2)
         endif
         idflav(N_lepton+N_decayquarks+1) = 22
      case(WpA, WmA, ZA, GFZA, WWA, ZZA, ZZnA, WPZA, WMZA,WPAjj,WMAjj,
     &      QCDWPAjj,QCDWMAjj,QCDZAjj_l,QCDZAjj_n,ZAjj,ZAjj_n)
         do i = 1, N_lepton
            idflav(i) = finallep(i)
         enddo
         idflav(N_lepton+1) = 22
      case(WPH, WMH, WPHj, WMHj) 
         do i = 1, N_lepton
            idflav(i) = finallep(i)
         enddo
         idflav(N_lepton+1) = 25
      case(WP_only, WM_only, WPj, WMj) 
         do i = 1, N_lepton
            idflav(i) = finallep(i)
         enddo
      case(AA, GFAA, WPAA, WMAA, ZAA, ZnAA,WPAAj,WMAAj)
         do i = 1, N_lepton
            idflav(i) = finallep(i)
         enddo
         idflav(N_lepton+1) = 22
         idflav(N_lepton+2) = 22
      case(AAA)
         idflav(1) = 22
         idflav(2) = 22
         idflav(3) = 22
      case(HHjj_bbAA)
         idflav(1) = 5
         idflav(2) =-5
         idflav(3) = 22
         idflav(4) = 22
      case(HHjj_bbtau)
         idflav(1) = 5
         idflav(2) =-5
         idflav(3) =-15
         idflav(4) = 15
      case(HHjj)
         idflav(1) = 25
         idflav(2) = 25
      case default
c         print*,"no numbers for the final state particles assigned 
c     1            in leshouches.f"
      end select

      end
      
c------------------------------------------------------------------------------

      subroutine bostdq(idir,vv,pp,q)

      implicit none

c------------------------------------------------------------------------------
c...boost along arbitrary vector vv
c...four-vector pp is boosted from an actual frame to the rest frame
c...of the four-vector vv (for idir=1) or back (for idir=-1)
c...q is a resulting four-vector
c...note: vv must be time-like, pp may be arbitrary
c------------------------------------------------------------------------------
      
      double precision v(4),p(4),q(4),pp(4),vv(4)
      integer i,idir
      double precision amv,wsp

      do i=1,4
         v(i)=vv(i)
         p(i)=pp(i)
      enddo

      amv=(v(4)**2-v(1)**2-v(2)**2-v(3)**2)
      if (amv.le.0d0) then
        print*,"bosstv: warning amv**2=",amv
      endif
      amv=dsqrt(abs(amv))
      if (idir.eq.-1) then
        q(4)=( p(1)*v(1)+p(2)*v(2)+p(3)*v(3)+p(4)*v(4))/amv
        wsp =(q(4)+p(4))/(v(4)+amv)
      elseif (idir.eq.1) then
        q(4)=(-p(1)*v(1)-p(2)*v(2)-p(3)*v(3)+p(4)*v(4))/amv
        wsp =-(q(4)+p(4))/(v(4)+amv)
      else
        print*,"boostv: wrong value of idir = ",idir
      endif

      q(1)=p(1)+wsp*v(1)
      q(2)=p(2)+wsp*v(2)
      q(3)=p(3)+wsp*v(3)

      end

c------------------------------------------------------------------------------

      subroutine fillColorless(idflav,listposition)

      implicit none

#include "global.inc"
#include "lha.inc"
#include "process.inc"

      integer i1,i11,iput
      integer idflav(6),listposition,id1,id2,id3,id4,id5,id6
 
      listposition=1

c...particle identity code
      select case(procID)
      case(Hjj_WW, Hjjj_WW,WPWMjj,WPhadWMjj,WPWMhadjj,GFHjj_WW,
     &     Hjj_WPhadWM, Hjj_WPWMhad, Sp2jj_WW) !W+W-
         lidup(n_p+1,listposition)= 24
         lidup(n_p+2,listposition)=-24
      case(WPWPjj, WPhadWPjj,QCDWPWPjj) !W+W+
         lidup(n_p+1,listposition)= 24
         lidup(n_p+2,listposition)= 24
      case(WMWMjj, WMhadWMjj,QCDWMWMjj) !W-W-
         lidup(n_p+1,listposition)= -24
         lidup(n_p+2,listposition)= -24
      case(Hjj_ZZ_ll, Hjj_ZZ_lnu, Hjjj_ZZ_ll, Hjjj_ZZ_lnu, 
     1        ZZjj_ll, ZZjj_lnu, GFHjj_ZZ_ll, GFHjj_ZZ_lnu, 
     2        ZZhadjj, Hjj_ZZhad,
     3        Sp2jj_ZZ_ll, Sp2jj_ZZ_lnu) !ZZ
         lidup(n_p+1,listposition)= 23
         lidup(n_p+2,listposition)= 23
      case(Hjj, Hjj_AA, Hjj_mu, Hjj_tau, Hjj_bbar,
     -        GFHjj_AA, GFHjj_mu, GFHjj_tau, GFHjj_bbar)!H
         lidup(n_p+1,listposition)= 25
      case(AAjj)  ! spin-2 -> gamma gamma
         lidup(n_p+1,listposition)= 39
      case(Hjjj, Hjjj_AA, Hjjj_mu, Hjjj_tau, Hjjj_bbar) !H
         lidup(n_p+1,listposition)= 25
      case(HHjj)
         n_b = 0
         lidup(n_p+1,listposition)= 25
         lidup(n_p+2,listposition)= 25
      case(HHjj_bbtau, HHjj_bbAA)
         lidup(n_p+1,listposition)= 25
         lidup(n_p+2,listposition)= 25
      case(HAjj,  HAjj_AA, HAjj_mu, HAjj_tau, HAjj_bbar)
         lidup(n_p+1,listposition)= 25
         lidup(n_p+2,listposition)= 22
      case(HAjj_WW)
         lidup(n_p+1,listposition)= 24
         lidup(n_p+2,listposition)=-24
         lidup(n_p+3,listposition)= 22
      case(HAjj_ZZ_ll, HAjj_ZZ_lnu)
         lidup(n_p+1,listposition)= 23
         lidup(n_p+2,listposition)= 23
         lidup(n_p+3,listposition)= 22
      case(Ajj)
         lidup(n_p+1,listposition)= 22
      case(Zjj_l, Zjj_nu, QCDZjj_l, QCDZjj_nu)   !Z
         lidup(n_p+1,listposition)= 23
      case(WPjj,QCDWPjj)   !W+
         lidup(n_p+1,listposition)= 24 
      case(WMjj,QCDWMjj)   !W-
         lidup(n_p+1,listposition)=-24 
      case(WWZ, WPhadWMZ, WPWMhadZ, WWZhad) !WWZ
         lidup(n_p+1,listposition)= 24
         lidup(n_p+2,listposition)=-24
         lidup(n_p+3,listposition)= 23
      case(ZZWP, ZZhadWP, ZZWPhad) !ZZW+
         lidup(n_p+1,listposition)= 23
         lidup(n_p+2,listposition)= 23
         lidup(n_p+3,listposition)= 24
      case(ZZWM, ZZhadWM, ZZWMhad) !ZZW-
         lidup(n_p+1,listposition)= 23
         lidup(n_p+2,listposition)= 23
         lidup(n_p+3,listposition)=-24
      case(WWWP, WPhadWMWP, WPWMhadWP) !W+W-W+
         lidup(n_p+1,listposition)= 24
         lidup(n_p+2,listposition)=-24
         lidup(n_p+3,listposition)= 24
      case(WWWM, WMhadWPWM, WMWPhadWM) !W-W+W-
         lidup(n_p+1,listposition)=-24
         lidup(n_p+2,listposition)= 24
         lidup(n_p+3,listposition)=-24
      case(ZZZ, ZZZhad) !ZZZ
         lidup(n_p+1,listposition)= 23
         lidup(n_p+2,listposition)= 23
         lidup(n_p+3,listposition)= 23
      case(WWA, WPhadWMA, WPWMhadA)
         n_b = 2
         lidup(n_p+1,listposition)= 24
         lidup(n_p+2,listposition)=-24
         lidup(n_p+3,listposition)= 22
      case(ZZA, ZZhadA, ZZnA)
         n_b = 2
         lidup(n_p+1,listposition)= 23
         lidup(n_p+2,listposition)= 23
         lidup(n_p+3,listposition)= 22
      case(WPZA, WPhadZA, WPZhadA)
         n_b = 2
         lidup(n_p+1,listposition)= 24
         lidup(n_p+2,listposition)= 23
         lidup(n_p+3,listposition)= 22
      case(WMZA, WMhadZA, WMZhadA)
         n_b = 2
         lidup(n_p+1,listposition)=-24
         lidup(n_p+2,listposition)= 23
         lidup(n_p+3,listposition)= 22
      case(WPAA,WPAAj)
         n_b = 1
         lidup(n_p+1,listposition)= 24
         lidup(n_p+2,listposition)= 22
         lidup(n_p+3,listposition)= 22
      case(WMAA,WMAAj)
         n_b = 1
         lidup(n_p+1,listposition)=-24
         lidup(n_p+2,listposition)= 22
         lidup(n_p+3,listposition)= 22
      case(ZAA,ZnAA)
         n_b = 1
         lidup(n_p+1,listposition)= 23
         lidup(n_p+2,listposition)= 22
         lidup(n_p+3,listposition)= 22
      case(AAA)
         n_b = 0
         lidup(n_p+1,listposition)= 22
         lidup(n_p+2,listposition)= 22
         lidup(n_p+3,listposition)= 22
      case(WW,GFWW,
     &     WPhadWMlep,WPlepWMhad,GFWPhadWMlep,GFWPlepWMhad) !WW
         lidup(n_p+1,listposition)= 24
         lidup(n_p+2,listposition)=-24
      case(WpZ,WPhadZ,WPZhad) !W+Z
         lidup(n_p+1,listposition)= 24
         lidup(n_p+2,listposition)= 23
      case(WmZ,WMhadZ,WMZhad) !W-Z
         lidup(n_p+1,listposition)=-24
         lidup(n_p+2,listposition)= 23
      case(ZZ,GFZZ, !ZZ
     &     ZZhad,GFZZhad,
     &     QCDZZjj_ll,QCDZZjj_lnu)
         lidup(n_p+1,listposition)= 23
         lidup(n_p+2,listposition)= 23
      case(WpA) !W+A
         n_b = 1
         lidup(n_p+1,listposition)= 24
         lidup(n_p+2,listposition)= 22
      case(WmA) !W-A
         n_b = 1
         lidup(n_p+1,listposition)=-24
         lidup(n_p+2,listposition)= 22
      case(ZA,GFZA) !ZA
         n_b = 1
         lidup(n_p+1,listposition)= 23
         lidup(n_p+2,listposition)= 22
      case(AA,GFAA) !AA
         n_b = 0
         lidup(n_p+1,listposition)= 22
         lidup(n_p+2,listposition)= 22
      case(WPH, WPHj) !W+H
         n_b = 2
         lidup(n_p+1,listposition)= 24
         lidup(n_p+2,listposition)= 25 
         ! higgsdecay lhe output not yet implmented
      case(WMH, WMHj) !W-H
         n_b = 2
         lidup(n_p+1,listposition)=-24
         lidup(n_p+2,listposition)= 25 
         ! higgsdecay lhe output not yet implmented
      case(WP_only, WPj) !W+
         n_b = 1
         lidup(n_p+1,listposition)= 24
      case(WM_only, WMj) !W-
         n_b = 1
         lidup(n_p+1,listposition)=-24
      case(WPZjj, WPhadZjj, WPZhadjj,QCDWPZjj) !W+Zjj
         lidup(n_p+1,listposition)= 24
         lidup(n_p+2,listposition)= 23
      case(WMZjj, WMhadZjj, WMZhadjj,QCDWMZjj) !W-Zjj
         lidup(n_p+1,listposition)=-24
         lidup(n_p+2,listposition)= 23
      case(WPAjj,QCDWPAjj) !W+Ajj
         n_b=1
         lidup(n_p+1,listposition)= 24
         lidup(n_p+2,listposition)= 22
      case(WMAjj,QCDWMAjj) !W-Ajj
         n_b=1
         lidup(n_p+1,listposition)=-24
         lidup(n_p+2,listposition)= 22
      case(ZAjj,ZAjj_n) !ZAjj
         n_b=1
         lidup(n_p+1,listposition)= 23
         lidup(n_p+2,listposition)= 22
      case(QCDZAjj_l,QCDZAjj_n) !ZAjj
         n_b=1
         lidup(n_p+1,listposition)= 23
         lidup(n_p+2,listposition)= 22
      end select
 
      do i1=1,n_b
         i11=n_p+i1
c...colour information
         licolup(1,i11,listposition)=0
         licolup(2,i11,listposition)=0
         
c...orgin information of vector boson(s)
         lmothup(1,i11,listposition)=1
         lmothup(2,i11,listposition)=2
c...status code of particle (=-1 incoming particle =1 outgoing)
         listup(i11,listposition)=2

c...set higgs/photon as final particle in some processes with stable higgs
         select case(procID)
         case(HAjj,Hjj,GFHjj,Hjjj,Ajj)  
            listup(i11,listposition)=1
         case(WPH, WMH, WPHj, WMHj)
            if (i1.eq.2) listup(i11,listposition)=1
         end select
      enddo

c...setting photon as final particle in VBF Higgs + photon production
      select case(process)
      case(HAjj_AA,HAjj_mu,HAjj_tau,HAjj_bbar) 
         listup(6,listposition)=1
      case(HAjj_WW,HAjj_ZZ_ll,HAjj_ZZ_lnu)
         listup(7,listposition)=1
      case(HAjjjLO_AA,HAjjjLO_mu,HAjjjLO_tau,HAjjjLO_bbar) 
         listup(7,listposition)=1
      case(HAjjjLO_WW,HAjjjLO_ZZ_ll,HAjjjLO_ZZ_lnu)
         listup(8,listposition)=1
      end select

c...setting spin and mothers of photon(s) in diboson and triboson processes
      select case(procid)
      case(WPA,WMA,ZA,WWA,ZZA,ZZnA,WPZA,WMZA,WPAjj,WMAjj,QCDWPAjj,QCDWMAjj,
     &     QCDZAjj_l,QCDZAjj_n,ZAjj,ZAjj_n,
     &     WPhadWMA, WPWMhadA, ZZhadA, WPhadZA, WPZhadA, WMhadZA, WMZhadA)  !  one final photon
         lmothup(1,nup,listposition) = 1
         lmothup(2,nup,listposition) = 2
      case(AA,WPAA,WMAA,ZAA,ZnAA)  !  two final photons
         lmothup(1,nup,listposition) = 1
         lmothup(2,nup,listposition) = 2
         lmothup(1,nup-1,listposition) = 1
         lmothup(2,nup-1,listposition) = 2
      case(AAA)  ! three final photons
         lmothup(1,nup,listposition) = 1
         lmothup(2,nup,listposition) = 2
         lmothup(1,nup-1,listposition) = 1
         lmothup(2,nup-1,listposition) = 2
         lmothup(1,nup-2,listposition) = 1
         lmothup(2,nup-2,listposition) = 2
      case(HHjj)
         lmothup(1,nup,listposition) = 1
         lmothup(2,nup,listposition) = 2
         lmothup(1,nup-1,listposition) = 1
         lmothup(2,nup-1,listposition) = 2
      end select

      call kf_leptons(idflav)

      id1=idflav(1)
      id2=idflav(2)
      id3=idflav(3)
      id4=idflav(4)
      id5=idflav(5)
      id6=idflav(6)

      do i1=1,n_v
         i11=n_p+n_b+i1
         if(i11.le.nup)then
            lidup(i11,listposition)=idflav(i1) 
            listup(i11,listposition)= 1
            licolup(1,i11,listposition)= 0
            licolup(2,i11,listposition)= 0
            if(idflav(i1).eq.5)then
              if(process.eq.Hjj_bbar.or.process.eq.HAjj_bbar
     $                .or.process.eq.HHjj_bbtau
     $                .or.process.eq.HHjj_bbAA)then
                 licolup(1,i11,listposition)= 503
                 licolup(2,i11,listposition)= 0
              elseif(process.eq.HjjjLO_bbar.or.process.eq.HAjjjLO_bbar
     $                .or.process.eq.Hjjj_bbar
     $                .or.process.eq.HHjjjLO_bbtau
     $                .or.process.eq.HHjjjLO_bbAA)then
                 licolup(1,i11,listposition)= 504
                 licolup(2,i11,listposition)= 0
              elseif(process.eq.HjjjjLO_bbar)then
                 licolup(1,i11,listposition)= 505
                 licolup(2,i11,listposition)= 0
              else if(process.eq.GFHjj_bbar)then
                 licolup(1,i11,listposition)= 506
                 licolup(2,i11,listposition)= 0
              endif
            elseif(idflav(i1).eq.-5) then
              if(process.eq.Hjj_bbar.or.process.eq.HAjj_bbar
     $                .or.process.eq.HHjj_bbtau
     $                .or.process.eq.HHjj_bbAA)then            
                 licolup(1,i11,listposition)= 0
                 licolup(2,i11,listposition)= 503               
              elseif(process.eq.HjjjLO_bbar.or.process.eq.HAjjjLO_bbar
     $                .or.process.eq.Hjjj_bbar
     $                .or.process.eq.HHjjjLO_bbtau
     $                .or.process.eq.HHjjjLO_bbAA)then
                 licolup(1,i11,listposition)= 0
                 licolup(2,i11,listposition)= 504
              elseif(process.eq.HjjjjLO_bbar)then
                 licolup(1,i11,listposition)= 0
                 licolup(2,i11,listposition)= 505
              else if(process.eq.GFHjj_bbar)then
                 licolup(1,i11,listposition)= 0
                 licolup(2,i11,listposition)= 506
              endif
            endif

            if(idflav(i1).le.5 .and. idflav(i1).ge.1) then
              select case (process)
              case (WPhadWMlep,WPlepWMhad,GFWPhadWMlep,GFWPlepWMhad,
     &              ZZhad,GFZZhad,WPhadZ,WPZhad,WMhadZ,WMZhad,
     &              WPhadWMZ, WPWMhadZ, WWZhad, ZZhadWP, ZZWPhad, ZZhadWM, ZZWMhad,
     &              WPhadWMWP, WPWMhadWP, WMhadWPWM, WMWPhadWM, ZZZhad,
     &              WPhadWMA, WPWMhadA, ZZhadA, WPhadZA, WPZhadA, WMhadZA, WMZhadA)
                 licolup(1,i11,listposition)= 502
                 licolup(2,i11,listposition)= 0
              case (WPhadWMlepjLO,WPlepWMhadjLO,ZZhadjLO,
     &              WPhadZjLO,WPZhadjLO,WMhadZjLO,WMZhadjLO,
     &              WPhadWMZjLO, WPWMhadZjLO, WWZhadjLO, ZZhadWPjLO, ZZWPhadjLO, ZZhadWMjLO, ZZWMhadjLO,
     &              WPhadWMWPjLO, WPWMhadWPjLO, WMhadWPWMjLO, WMWPhadWMjLO, ZZZhadjLO,
     &              WPhadWMAjLO, WPWMhadAjLO, ZZhadAjLO, WPhadZAjLO, WPZhadAjLO, WMhadZAjLO, WMZhadAjLO,
     &              WPhadWMjj, WPWMhadjj, WPhadZjj, WPZhadjj, WMhadZjj, WMZhadjj, ZZhadjj,
     &              Hjj_WPhadWM, Hjj_WPWMhad, Hjj_ZZhad, WPhadWPjj, WMhadWMjj)
                 licolup(1,i11,listposition)= 503
                 licolup(2,i11,listposition)= 0
              case (WPhadWMjjjLO, WPWMhadjjjLO, WPhadZjjjLO, WPZhadjjjLO, WMhadZjjjLO, WMZhadjjjLO, ZZhadjjjLO,
     &              HjjjLO_WPhadWM, HjjjLO_WPWMhad, HjjjLO_ZZhad, WPhadWPjjjLO, WMhadWMjjjLO)
                 licolup(1,i11,listposition)= 504
                 licolup(2,i11,listposition)= 0
              end select
            elseif(idflav(i1).le.-1 .and. idflav(i1).ge.-5) then
              select case (process)
              case (WPhadWMlep,WPlepWMhad,GFWPhadWMlep,GFWPlepWMhad,
     &              ZZhad,GFZZhad,WPhadZ,WPZhad,WMhadZ,WMZhad,
     &              WPhadWMZ, WPWMhadZ, WWZhad, ZZhadWP, ZZWPhad, ZZhadWM, ZZWMhad,
     &              WPhadWMWP, WPWMhadWP, WMhadWPWM, WMWPhadWM, ZZZhad,
     &              WPhadWMA, WPWMhadA, ZZhadA, WPhadZA, WPZhadA, WMhadZA, WMZhadA)
                 licolup(1,i11,listposition)= 0
                 licolup(2,i11,listposition)= 502
              case (WPhadWMlepjLO,WPlepWMhadjLO,ZZhadjLO,
     &              WPhadZjLO,WPZhadjLO,WMhadZjLO,WMZhadjLO,
     &              WPhadWMZjLO, WPWMhadZjLO, WWZhadjLO, ZZhadWPjLO, ZZWPhadjLO, ZZhadWMjLO, ZZWMhadjLO,
     &              WPhadWMWPjLO, WPWMhadWPjLO, WMhadWPWMjLO, WMWPhadWMjLO, ZZZhadjLO,
     &              WPhadWMAjLO, WPWMhadAjLO, ZZhadAjLO, WPhadZAjLO, WPZhadAjLO, WMhadZAjLO, WMZhadAjLO,
     &              WPhadWMjj, WPWMhadjj, WPhadZjj, WPZhadjj, WMhadZjj, WMZhadjj, ZZhadjj,
     &              Hjj_WPhadWM, Hjj_WPWMhad, Hjj_ZZhad, WPhadWPjj, WMhadWMjj)
                 licolup(1,i11,listposition)= 0
                 licolup(2,i11,listposition)= 503
              case (WPhadWMjjjLO, WPWMhadjjjLO, WPhadZjjjLO, WPZhadjjjLO, WMhadZjjjLO, WMZhadjjjLO, ZZhadjjjLO,
     &              HjjjLO_WPhadWM, HjjjLO_WPWMhad, HjjjLO_ZZhad, WPhadWPjjjLO, WMhadWMjjjLO)
                 licolup(1,i11,listposition)= 0
                 licolup(2,i11,listposition)= 504
              end select
            endif
         endif
      enddo   

      iput=0
      do i1=1,n_b*2,2
         iput=iput+1
         i11=n_p+n_b+i1
         if(i11.le.nup)then
            lmothup(1,i11,listposition)=n_p+iput
            lmothup(2,i11,listposition)=0
            lmothup(1,i11+1,listposition)=n_p+iput
            lmothup(2,i11+1,listposition)=0
         endif
      enddo      

      end
