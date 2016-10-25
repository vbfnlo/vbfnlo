***************************************************************************
*
*  This is the GGFLO main program. 
*
*  This is basically a duplicate of VBFNLO, but without the NLO part 
*  and some different initialization routines.
*
****************************************************************************

****************************************************************************
*  Last modified: 29.07.2008
****************************************************************************

      subroutine doggflo
      implicit none

c include global variables via include-files, e.g. "procID"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/tau_hel.inc"
#include "VBFNLO/utilities/process.inc"

c define local variables

c     variables for the integration loops
      integer iteration, seed
      integer*8 Ncall, i, two
      parameter (two = 2)

c... control bad GramDets
      Logical Singular
      integer counter(10)
      COMMON/DetCount/ counter,iteration,Singular

c     digit characters for numbering the grid output file
      character*1 digit(61)
      data digit /'1','2','3','4','5','6','7','8','9','A',
     1            'B','C','D','E','F','G','H','I','J','K',
     1            'L','M','N','O','P','Q','R','S','T','U',
     1            'V','W','X','Y','Z','a','b','c','d','e',
     1            'f','g','h','i','j','k','l','m','n','o',
     1            'p','q','r','s','t','u','v','w','x','y',
     1            'z'/
c     the filenames for the grid output
      character*250 file_name_out, file_name_in

c     random number array controlled by monaco
      real*8 rand(max_PS_dim),rn

c     4-momenta of the partons involved in the basic QCD 2->2(+1)+X process.
c     The partons are the incoming quarks/gluons + the 2 tagging jets + 
c     the real emission (for NLO only)
      real*8 p(0:3,max_p,max_kin) 

c     Additional particles like a higgs boson, or any decay products 
c     (also hadronic ones). in short: everything which is not a parton 
c     involved in the basic QCD 2->2(+1)+X process.
      real*8 v(0:3,max_v,max_kin) 

c     Jets after jet-defining function
      real*8 jets(0:7,max_jets,max_kin) 
      integer n_jets(max_kin)
c     Charged leptons, i.e. electrons, muons, taus but not neutrinos
      real*8 leptons(0:7,max_v,max_kin)
      integer n_leptons(max_kin)
c     photons
      real*8 photons(0:7,max_v,max_kin)
      integer n_photons(max_kin)
c     invisible particles, i.e. neutrinos, LSP
      real*8 invisible(0:8,max_v,max_kin)
      integer n_invisible(max_kin)          
      
c     event weights, jacobi-factors, amplitude-squares
      real*8 weight, dps, m2s(0:3), dsig(0:3)

c     the result of the integration: crossection, error, chisquare
      real*8 chi2
      real*8 NLOxsection(-2:5), NLOsdev(-2:5), NLOchi2(-2:5)

c     other variables
      real*8 x(nx)   ! feynman x parameters
      logical cut_ok ! whether event passes cuts
      real*8 xuz(2,2:3)
      integer gnlo, L, mu
      logical lokt(3)
      common /nlovariable/ gnlo
      logical noEventOutput

c     variables for unweighting in LHE output
      integer estimatedEventNumber, estimatedEventNumber_old, getEstimatedEvents
      integer desiredPlus5Sigma
      real*8 increaseNcallRatio
      external getEstimatedEvents
      logical extendForEvents
      integer*8 Ncall_LO

c     external functions
      logical Cuts     ! general cut-function
      real*8 amplitude ! general squared-amplitude function
      double precision reweight_grid      ! modify grid generation for PS-dependent statistics variation in histograms
      external reweight_grid
      external Cuts, amplitude

      real*8 RandomNumber
      external RandomNumber     

** Functions and variables to fully define the input gridname
      character*250 GetInputPath
      character*250 path
      external GetInputPath    

c     debugging and testing variables and switches
      logical showtiming,GramCheck
      parameter (showtiming=.false.)

      call PrintVBFNLOVersion(2)
      
** LoopTools initialisation
#ifdef WITH_LT
      call initLT      
#endif

      if (showtiming) call mtime

c     first, do all the initialization
      
      call InitGlobalParameters_GGFLO

      call PrintProcInfo(.false.)

      call InitPDFs(1)
      
      call InitCouplings

*      call InitScales

      call InitProcess

      call InitCuts
      
      call InitPhaseSpace
    
      call InitRandomNumbers(seed)

      call InitRandomNumbers(seed)

      call InitHistograms


** Setting input path
      path = GetInputPath()

     

c...Les Houches interface
      noEventOutput = .false.
      if ((lha.or.hepmc).and..not.doNLO) then 
         call lha_file(-1)
      endif


c... set default value for "Singular"
      Singular = .false.


c==========================
c    The LO Calculation
c==========================

      ! no multi-channel processes in GGFLO
      multiChannelUnweighting = .false.
    
      if (doBorn) then

** FF initialisation
         call initFF
         call initGaugeTest

*         call InitPDFs(0)

         gnlo = 0

         file_name_in  = path(1:len_trim(path))//"/"//
     &        gridname1(1)
     &        (1:INDEX(gridname1(1),'  ')-1)
         file_name_in = file_name_in(1:len_trim(file_name_in))


         DO iteration = 1, N_iterations1
            Ncall = two**( N_points(0) - N_iterations1 + iteration )
           
            file_name_out = gridname1(1)(1:INDEX(gridname1(1),'  ')-1)//
     &                   '.out.'//digit(iteration)
            if ( iteration .eq. 1 ) then
               call monaco_init( PS_dimension, Ncall, seed )
               call monaco_read(file_name_in)
            else
               call monaco_init2( PS_dimension, Ncall )
            end if       
            if (iteration.eq.N_iterations1) then
               call monran_set(1) ! save random number info
            endif

            call ResetWeights(1,1)

            ! run in last iteration can be extended when running with LHE output and the user
            ! requested a fixed number of unweighted events
            extendForEvents = .false.
            estimatedEventNumber_old = 0
            i = 1
            Ncall_LO = Ncall
            do while (i .le. Ncall_LO)

c... set default value for "Singular"
               Singular = .false.

c     retrieve a set of random numbers and the corresponding weight
               call monaco_get( rand, weight)

c     transform the random numbers to a physical phase space point
               call phasespace(rand, p, x, v, 1, dps  )

               if ( dps .gt. 0d0 ) then

                  if (alllep_flag) call lepton_gen

c     define jets, leptons, photons
                  call defineJets(p, v, 1, jets, n_jets)
                  call DefineLeptons(v, 1, leptons, n_leptons)
                  call DefineInvisible(v,1, invisible, n_invisible)
                  call DefinePhotons(v, 1, photons, n_photons)

c     apply cuts
                  cut_ok = Cuts(jets, n_jets, leptons, n_leptons,
     1                            invisible,n_invisible,photons, n_photons, 1)

                  if (cut_ok) then
c     calculate fac. and ren. scale for this phase space point
                     call Scales(p, v, jets, n_jets, leptons, 
     1                       n_leptons, photons, n_photons, 1, 0)
c     call the matrix element code and get an amplitude square
                     m2s(0) = Amplitude(rand, p, x, v, 0) 
                  else
                     dps = 0d0
                     m2s(0) = 0d0
                  end if 

                  dsig(0) = m2s(0)*dps

                  if (cut_ok) then
                     eventweight=dsig(0)*weight
c     in the last iteration : 
                     if (iteration.eq.N_iterations1) then

c     pass event information to the histogram routine
                        if (.not.extendForEvents) call HistogramEvent(weight,dsig(0),dps,
     &                                rand, p(0,1,1), x, v(0,1,1), 
     &                                jets(0,1,1), n_jets(1), 
     &                                leptons(0,1,1), n_leptons(1), 
     &                                invisible(0,1,1),n_invisible(1),
     &                                photons(0,1,1), n_photons(1), 0)

c     the unweighting of events : 
                        
c                        if (unweighting) then 
                           call Reweight(eventweight, extendForEvents) ! reweight
c                        endif

                        if (eventweight.ne.0d0) then ! event is to be kept

                           if((lha.or.hepmc).and..not.doNLO) call lha_put(p,v) !  Les Houches interface

                        endif
                     endif
                     
                  endif ! cutok

               else     ! dps > 0d0
                  m2s(0) = 0d0
                  dsig(0) = 0d0
               end if

c     pass weight to monaco
               if (.not.extendForEvents) call monaco_put(rand, weight, dsig(0), reweight_grid(p(0,1,1),v(0,1,1)))

               if (unweighting .and. iteration.eq.N_iterations1 .and. desiredEventCount.gt.0 .and. i.eq.NCall_LO) then
                  estimatedEventNumber = getEstimatedEvents()
                  desiredPlus5Sigma = int(dble(desiredEventCount) + 5.0*sqrt(dble(desiredEventCount)))
                  increaseNcallRatio = (dble(desiredPlus5Sigma) / dble(estimatedEventNumber))
                  if (increaseNcallRatio.gt.2d0) increaseNcallRatio = 2d0
                  if (increaseNcallRatio.lt.1.01d0) increaseNcallRatio = 1.01d0
                  print*, ""
                  if (partialUnweight) then
                     print*, "Estimated number of partially unweighted events so far:", estimatedEventNumber
                  else
                     print*, "Estimated number of unweighted events so far:", estimatedEventNumber
                  endif
                  if (estimatedEventNumber.lt.desiredPlus5Sigma) then
                     if (int((dble(Ncall_LO) * increaseNcallRatio)) / Ncall .lt. 65) then
                        Ncall_LO = int(dble(Ncall_LO) * increaseNcallRatio)
                        if (.not.extendForEvents) then
                           print*, ""
                           print*, ""
                           print*, "*****************************************************************************"
                           print*, ""
                           print*, "Cross section calculation finished, starting to generate additional events..."
                           print*, ""
                           print*, "*****************************************************************************"
                           print*, ""
                           print*, ""
                        endif
                        extendForEvents = .true.
                        print*, "Extending this iteration by a factor of", increaseNcallRatio,
     &                          "to get the desired number of events."
                        if (estimatedEventNumber.lt.estimatedEventNumber_old) then
                           print*, ""
                           print*, "Warning:"
                           print*, "The estimated number of unweighted events decreased"
                           print*, "due to a newly found larger maximal event weight."
                           print*, "Please consider increasing LO_POINTS if this occurs"
                           print*, "several times in a row!"
                           print*, ""
                        endif
                        estimatedEventNumber_old = estimatedEventNumber
                     else
                        print*, ""
                        print*, "Warning:"
                        print*, "Could not get enough unweighted events with the current"
                        print*, "precision on the cross section."
                        print*, "Please consider increasing LO_POINTS and/or LO_ITERATIONS."
                        print*, "Writing out all events obtained so far..."
                        print*, ""
                     endif
                  endif
               endif

               i = i + 1      ! increment loop variable

            end do
        
            call monaco_write( file_name_out )

         END DO

** closing FF
         call closeFF

         call monaco_result(xsection, sdev, chi2)
     
c     print out the result. 
         print *,""
         print *," result (LO): ",xsection," +- ",sdev," fb  "
     &        ,sdev/xsection*100d0,'%'
         print *,""
      else
         print *," No LO calculation performed."
      endif

c...Les Houches interface
      if((lha.or.hepmc).and..not.doNLO) then
         call lha_file(0)
         call lha_file(1)
      endif

      if(.not.doNLO) then
         print*,'number of events in the last iteration =',ncall_lo
         call reweightinginfo(.false.)
         call reweightinginfo(.true.)
      endif

      call WriteHistograms


** closing LoopTools
#ifdef WITH_LT
      call closeLT
#endif


* event output is commented out, because it's buggy and not standardized
*      call CloseEventOutput
     
c print out results in a file also
      open (20,
     &      file="xsection.out",
     &      status='REPLACE')
      write(20,*) xsection, sdev
      write(20,*) NLOxsection(-2), NLOsdev(-2)
      close(20)

      if (showtiming) call mtime

c...Les Houches interface
      if (lha.and..not.doNLO) then 
         print*,''
         print*,'created LHA event file for the LO calculation : ',
     &        lhafile
         print*,''
      endif
      if (hepmc.and..not.doNLO) then 
         print*,''
         print*,'created HepMC event file for the LO calculation : ',
     &        hepmcfile
         print*,''
      endif
      if (noEventOutput) then 
         print*,''
         print*,'Sorry, event output is not available for this process! '
      endif
      if ((lha.or.hepmc).and..not.doNLO) then
         if (unweighting) then
           if (partialUnweight) then
             write(*,'(A,I8)')         ' number of events with weight = 1 written to file                      : ', nevacc
             if (nevovr.gt.0) then
                write(*,'(A,G11.4,A,I8)') ' number of events with weight between 1 and',xmaxup(1),' written to file : ', nevovr
             else
                print*, "No events with weight > 1 occured ==> events are fully unweighted!"
             endif
           else
             print*,'number of unweighted events written to file : ', nevacc
           endif
         else
           print*,'number of weighted events written to file : ', nevall
         endif
         print*,''
      endif

c.. Printing summary of results

      write(*,*)' '
      write(*,*)' '
      print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' //
     &     '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

      write(*,*)'    VBFNLO ', vbfnloversionstring(1:len_trim(vbfnloversionstring))

      call printFinalProc

      if (doBorn) then
         print *," TOTAL result (LO):       ",
     &        xsection," +- ",sdev," fb  ",
     &        sdev/xsection*100d0,'%'
      end if

      print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' //
     &     '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(*,*)' '
      write(*,*)' '
   
      STOP
      END
c     END OF GGFLO
    
