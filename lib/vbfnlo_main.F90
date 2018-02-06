!**************************************************************************
!
! This is the VBFNLO main program.
!
!***************************************************************************

module vbfnlo_main
   use monaco
   use globalvars, only: lglobalprint
   use VBFNLOVersion, only: setVersion, printVersion
   implicit none

   ! include global variables via include-files, e.g. "procID"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/tau_hel.inc"
#include "VBFNLO/utilities/process.inc"
#include "VBFNLO/utilities/scales.inc"

   ! define local variables
   !--------------------------------------------------------
   ! variables for the integration loops
   integer ps_number
   integer iteration
   integer*8 Ncall, i, two
   parameter (two = 2)
   common/Ncall1/Ncall,iteration

   ! digit characters for numbering the grid output file
   character*1 digit(61)
   data digit /'1','2','3','4','5','6','7','8','9','A', &
               'B','C','D','E','F','G','H','I','J','K', &
               'L','M','N','O','P','Q','R','S','T','U', &
               'V','W','X','Y','Z','a','b','c','d','e', &
               'f','g','h','i','j','k','l','m','n','o', &
               'p','q','r','s','t','u','v','w','x','y', &
               'z'/
   ! the filenames for the grid output
   character*250 file_name_out, file_name_in

   ! random number array controlled by monaco
   real*8 rand(max_PS_dim)

   ! 4-momenta of the partons involved in the basic QCD 2->2(+1)+X process.
   ! The partons are the incoming quarks/gluons + the 2 tagging jets +
   ! the real emission (for NLO only)
   real*8 p(0:3,max_p,max_kin)

   ! Additional particles like a higgs boson, or any decay products
   ! (also hadronic ones). in short: everything which is not a parton
   ! involved in the basic QCD 2->2(+1)+X process.
   real*8 v(0:3,max_v,max_kin)

   ! Jets after jet-defining function
   real*8 jets(0:7,max_jets,max_kin)
   integer n_jets(max_kin)
   ! Charged leptons, i.e. electrons, muons, taus but not neutrinos
   real*8 leptons(0:8,max_v,max_kin)
   integer n_leptons(max_kin)
   ! photons
   real*8 photons(0:7,max_v,max_kin)
   integer n_photons(max_kin)
   ! invisible particles, i.e. neutrinos, LSP
   real*8 invisible(0:8,max_v,max_kin)
   integer n_invisible(max_kin)

   ! event weights, jacobi-factors, amplitude-squares
   real*8 weight, dps, m2s(0:max_kin), dsig(0:max_kin)

   ! the result of the integration: chisquare
   real*8 chi2

   ! other variables
   real*8 x(nx)   ! feynman x parameters
   logical cut_ok ! whether event passes cuts
   logical ps_ok  ! phase space generator is suitable for this ps point
   real*8 xuz(2,2:max_kin)
   integer gnlo, L, nd
   logical lokt(max_kin),realcont
   common /nlovariable/ gnlo
   logical noEventOutput

   ! variables for unweighting in LHE output
   integer estimatedEventNumber, estimatedEventNumber_old, getEstimatedEvents
   integer desiredPlus5Sigma, maxNcallFactor, k, N_iterations1_orig
   real*8 increaseNcallRatio
   double precision unwSumXS
   external getEstimatedEvents
   logical extendForEvents, multiUnweightingFirstRun, multiUnweightingSecondRun
   integer*8 Ncall_LO

   ! external functions
   logical Cuts       ! general cut-function
   logical Choose_PS  ! general PS-function
   real*8 amplitude   ! general squared-amplitude function
   external Cuts, amplitude, Choose_PS

   real*8 RandomNumber
   external RandomNumber

   !* Functions and variables to fully define the input gridname
   character*250 GetInputPath
   character*250 path
   external GetInputPath

   ! debugging and testing variables and switches
   real time,time0

   ! debug variable to output dipole test data
   logical ldiptest
   parameter (ldiptest=.false.)

   logical Finalprinting
   parameter(Finalprinting=.true.)

!... control bad GramDets
      Logical Singular
      integer counter(10)
      COMMON/DetCount/ Singular

#ifdef WITH_MPI
   ! MPI variables
   integer :: ierr
#endif


contains
   subroutine dovbfnlo
#ifdef WITH_MPI
      use vbfnlo_mpi
#endif
      use globalvars, only: isggflo, seed
      use stopwatch, only: starttimer, stoptimer, printtime, printduration

      lglobalprint = .true.
#ifdef WITH_MPI
      call InitMPI()
#endif

      call SetVersion
      if (lglobalprint) call PrintVersion
      !* LoopTools initialisation
#ifdef WITH_LT
      call initLT
#endif



      call InitGlobalParameters
      if (isggflo) infloops = .true.

      if (lglobalprint) call PrintProcInfo(.false.)

      call InitPDFs(1)

      call InitCouplings

      ! call InitScales

      call InitProcess
      call InitCuts
      call InitPhaseSpace
      call InitRandomNumbers(seed)
      call InitHistograms


      ! make grid directory
      IF((Loops_sub_LO.gt.1).or.(Loops_sub_NLO.gt.1)) then
         CALL system("mkdir -p GRID")
      endif

      if (lwritedata) call OpenDataFile(seed)


      !* Setting input path
      path = GetInputPath()


!... set default value for "Singular"
      Singular = .false.


      call lha_output

      boxcountm  = 0
      boxcount2m = 0
      box2countm  = 0
      box2count2m = 0
      pentcountm  = 0
      pentcount2m = 0
      hexcountm  = 0
      hexcount2m = 0

      !==========================
      ! The LO Calculation
      !==========================

      final_xsec = 0d0
      final_sdev2= 0d0
      final_time = 0d0

      ! initialize unweighting of events for multichannel processes
      multiUnweightingFirstRun = .false.
      multiUnweightingSecondRun = .false.
      multiChannelUnweighting = .false.

      if (showtiming .and. lglobalprint) then
          call printtime('Starttime: ')
          call starttimer
      endif

      if (doBorn) then
         call born_loop
      else
         print *," No LO calculation performed."
      endif

      call check_xsec(final_xsec(0,0,-1), 'LO')

      !...Les Houches interface
      if((lha.or.hepmc).and..not.doNLO) then
         xsection = final_xsec(0,0,-1)
         sdev = sqrt(final_sdev2(0,0,-1))
         call lha_file(0)
         call lha_file(1)
      endif

      if(.not.doNLO) call reweightinginfo(.true.)


#ifdef WITH_NLO
      !* FF initialisation
      if (.not. isggflo) then
          call initFF

          if ((doFloops).and.hasFLoops.and.(floops.ne.0)) then
             call fermion_loop
          else
             if(hasFLoops) print *," No NLO fermion loop calculation performed."
          endif
          call check_xsec(final_xsec(0,0,0), 'NLO fermion loop')

          if ((doVirtuals.and.hasNLO) .or. ewcor_switch) then
             call virtuals_loop
          endif
          call check_xsec(final_xsec(0,0,0), 'NLO virtual')

          !* closing FF
          call closeFF

          if (doEmissions.and.hasNLO) then
             call reals_loop
          endif

          if (ewcor_switch .and. (sector .ge. 3)) then
             call real_photons_loop
          end if
      endif

      !* closing LoopTools
#ifdef WITH_LT
      call closeLT
#endif

      if (((dovirtuals.and.doEmissions.and.hasNLO).or.ewcor_switch) .and. lglobalprint) then
         print *," "
         print *," final result at NLO "
         print *," sigma = ",final_xsec(0,0,0)," +- ",sqrt(final_sdev2(0,0,0))," fb" &
            ,sqrt(final_sdev2(0,0,0))/final_xsec(0,0,0)*100d0,'%'
      endif

#endif
      if ((doVirtuals.and.hasNLO) .or. ewcor_switch .or. isggflo) call Final_Instabilities
      ! if(FinalPrinting .and. (NLO_Loops .gt. 1) .and. dovirtuals)
      if(FinalPrinting .and. lglobalprint) call Final_Printing

      if (lwritedata .or. lreaddata) call CloseDataFile()

      call final_lha_out

#ifdef WITH_MPI
      call combine_histograms
      if (vbfnlo_mpi_myid == vbfnlo_mpi_masterid) &
#endif
         call WriteHistograms

      ! print out results in a file also
      open (20, &
         file=xsecfile, &
         status='REPLACE')
      write(20,*) final_xsec(0,0,-1),sqrt(final_sdev2(0,0,-1))
      write(20,*) final_xsec(0,0,0),sqrt(final_sdev2(0,0,0))
      close(20)

#ifdef WITH_MPI
      call FinalMPI
#endif
      !.. Printing summary of results

      if (lglobalprint) then
         write(*,*)' '
         write(*,*)' '
         print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' // &
         '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

         call printFinalProc
         call print_final_xsection
         call print_virtual_instabilities

         if (doNLO .or. ewcor_switch) then
            write(*,*) " K-Factor:               ",  final_xsec(0,0,0)/final_xsec(0,0,-1)
         end if
         print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' // &
         '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write(*,*)' '
         write(*,*)' '
      endif

      if (showtiming .and. lglobalprint) then
          call printtime('  Endtime: ')
          call stoptimer
          call printduration
      endif

   end subroutine

   subroutine check_xsec(xsec, typestr)
      character(len=*) typestr
      real*8 xsec
      if (IsNaN(xsec)) then
         print*, 'Sorry, something has gone wrong with the'
         print*, trim(typestr) // 'cross section!'
         stop
      endif
   endsubroutine

   subroutine born_loop
#ifdef WITH_MPI
      use vbfnlo_mpi
#endif
      use globalvars, only: isggflo

      if (isggflo) then
!* FF initialisation
         call initFF
         call initGaugeTest

         call Init_Instabilities(1)
     endif
      
      ! call InitPDFs(0)
      gnlo = 0
      do     ! loop for multichannel processes and definite number of unweighted events

         call born_unweighting_runs_init

         do ps_number = 1, PS_Loops !multi-channel
            if (PS_Loops.gt.1 .and. lglobalprint) then
               write(*,*)" "
               write(*,*)"PHASE SPACE", ps_number
               write(*,*)" "
            endif

            do sub_number = 1, Loops_sub_LO ! multi-processes
               call cpu_time(time0)

               if (Loops_sub_LO.gt.1 .and. lglobalprint) then
                  write(*,*)" "
                  write(*,*)"Subprocess", sub_number
                  write(*,*)" "
               endif

               if (multiUnweightingSecondRun) then
                  ! only one iteration for second run, reuse last grid from first run
                  if (Loops_sub_LO.gt.1) then
                     if (N_iterations1_orig.eq.1) then
                        file_name_in = get_filename_in(sub_number, ps_number, prefix_grid=.false., append_out=.false.)
                     else
                        file_name_in = get_filename_in(sub_number, ps_number, prefix_grid=.true., append_out=.true.)
                     endif
                  else
                     if (N_iterations1_orig.eq.1) then
                        file_name_in = get_filename_in(-1, ps_number, prefix_grid=.false., append_out=.false.)
                     else
                        file_name_in = get_filename_in(-1, ps_number, prefix_grid=.false., append_out=.true.)
                     endif
                  end if
               else
                  if (Loops_sub_LO.gt.1) then
                     file_name_in = get_filename_in(sub_number, ps_number, prefix_grid=.false., append_out=.false.)
                  else
                     file_name_in = get_filename_in(-1, ps_number, prefix_grid=.false., append_out=.false.)
                  end if
               endif

               xsection=0d0
               sdev=0d0

               if (lreaddata) then
                  call numdata(Ncall,PS_dimension,1,PS_number)
                  N_iterations1 = 1
               endif

               do iteration = 1, N_iterations1
                  ! Adjust statistics for different phase space
                  call statistics_proc(ps_number)

                  if (Loops_sub_LO.gt.1) then
                     file_name_out = 'GRID/'//gridname1(ps_number) &
                        (1:INDEX(gridname1(ps_number),'  ')-1)// &
                        '_'//digit(sub_number)//'.out.'//digit(iteration)
                  else
                     file_name_out = gridname1(ps_number) &
                        (1:INDEX(gridname1(ps_number),'  ')-1)// &
                        '.out.'//digit(iteration)
                  end if

                  if ( iteration .eq. 1 ) then
                     call monaco_init( PS_dimension, Ncall)
                     call monaco_read(file_name_in)
                  else
                     call monaco_init2( PS_dimension, Ncall )
                  end if
                  ! if (iteration.eq.N_iterations1) then
                  !    call monran_set(1) ! save random number info, TODO: remove
                  ! endif
                  call ResetWeights(ps_number,sub_number)

                  ! run in last iteration can be extended when running with LHE output and the user
                  ! requested a fixed number of unweighted events
                  extendForEvents = .false.
                  estimatedEventNumber_old = 0
                  ! i = 0
                  Ncall_LO = Ncall
#ifdef WITH_MPI
                  call mpi_prepare_loop(iteration, Ncall_LO)
                  do i = 1, ncall_thisproc
#else
                  do i = 1, NCall_LO
#endif
                     call born_phasespacepoint
                  enddo
#ifdef WITH_MPI
                  call mpi_after_loop(NCall_LO)
#else
         call monaco_end_of_iteration()
#endif

#ifdef WITH_MPI
                  if (vbfnlo_mpi_myid == vbfnlo_mpi_masterid) then
#endif
                  if (.not.multiUnweightingSecondRun) call monaco_write( file_name_out )
#ifdef WITH_MPI
                  endif
#endif
               enddo

               if (.not.multiUnweightingSecondRun) call monaco_result(xsection, sdev, chi2)

               if (multiChannelUnweighting) then
                  ! save max weight / XS of this part
                  if (multiUnweightingFirstRun) then
                     contribXS(contributionNumber) = xsection
                     numContribs = contributionNumber
                     if (numContribs.gt.maxContribs) then
                        print*, "Too many different phase spaces / subprocesses.  &
                           &Please increase parameter maxContribs in utilities/lha.inc!"
                        stop
                     endif
                  else
                     contribMaxWeight(contributionNumber) = 0.       ! init for lha_file(0)
                     contribEstimatedEvents(contributionNumber) = estimatedEventNumber
                  endif
                  contributionNumber = contributionNumber + 1
               else
                  contribEstimatedEvents(1) = estimatedEventNumber
               endif

               if (.not.multiUnweightingSecondRun) then
                  ! print out the result.
                  if (lglobalprint) then
                     print *,""
                     if(Loops_sub_LO.eq.1) then
                        print *," result (LO): ",xsection," +- ",sdev," fb  " &
                           ,sdev/xsection*100d0,'%'
                     else
                        print *," result sub (LO): ",xsection," +- ",sdev, &
                           " fb  ",sdev/xsection*100d0,'%'
                     end if
                     print *,""
                  endif

                  call cpu_time(time)
                  call AddToXsec(xsection,sdev**2,time-time0,sub_number,ps_number,-1)

               endif

               if(.not.doNLO .and. .not.multiUnweightingFirstRun) then
                  if (Loops_sub_LO.eq.1) then
                     if (PS_Loops.gt.1) print*,'in this phasespace:'
                  else
                     print*,'in this phasespace / subprocess:'
                  endif
                  print*,'number of events in the last iteration =',ncall_lo
                  call reweightinginfo(.false.)
               endif


            enddo ! end sub-processes

            if(Loops_sub_LO.gt.1 .and. .not.multiUnweightingSecondRun) then
               print *,""
               print *," result (LO): ", &
                  final_xsec(0,ps_number,-1)," +- ",sqrt(final_sdev2(0,ps_number,-1))," fb  ", &
                  sqrt(final_sdev2(0,ps_number,-1))/final_xsec(0,ps_number,-1)*100d0,'%'
               print *," "
            endif

         enddo                  !multi-channel

         if (isggflo) then
             call closeFF
         endif



         if (ps_loops.gt.1 .and. .not.multiUnweightingSecondRun) then
            print *,""
            print *," Total (multi-channel) result (LO): ", &
               final_xsec(0,0,-1)," +- ",sqrt(final_sdev2(0,0,-1))," fb  ", &
               sqrt(final_sdev2(0,0,-1))/final_xsec(0,0,-1)*100d0,'%'
            print *," "
         endif

         if (.not.multiUnweightingFirstRun) exit
      enddo

   endsubroutine

   subroutine born_phasespacepoint
      use globalvars, only: isggflo

      if (isggflo) then
          !... set default value for "Singular"
          Singular = .false.
      endif
    

      ! retrieve a set of random numbers and the corresponding weight
      call monaco_get( rand, weight)

      if (lreaddata) call readdata(PS_dimension,rand,weight,1,m2s(0))
      ! transform the random numbers to a physical phase space point
      call phasespace(rand, p, x, v, ps_number, dps  )

      if ( dps .gt. 0d0 ) then
         if (alllep_flag) call lepton_gen

         call DefineJets(p, v, 1, jets, n_jets)
         call DefineLeptons(v, 1, leptons, n_leptons)
         call DefineInvisible(v,1, invisible, n_invisible)
         call DefinePhotons(v, 1, photons, n_photons)

         ! in case of multi-channel: choose phase space:
         PS_OK = Choose_PS(ps_number,v, 1)
         cut_ok = Cuts(jets, n_jets, leptons, n_leptons,  &
            invisible,n_invisible,photons, n_photons, 1)

         if (cut_ok .and. PS_ok) then
            ! calculate fac. and ren. scale for this phase space point
            call Scales(p, v, jets, n_jets, leptons,  &
               n_leptons, photons, n_photons, 1, 0)
            ! call the matrix element code and get an amplitude square
            if (.not.lreaddata) &
               m2s(0) =  Amplitude(rand, p, x, v, 0, ps_number)

            if ((iteration.eq.N_iterations1).and. &
               lwritedata) call writedata(PS_dimension, &
               rand,weight,1,m2s(0),PS_number)
         else
            dps = 0d0
            m2s(0) = 0d0
         endif
         dsig(0) = m2s(0)*dps

         if (cut_ok .and. PS_ok) then
            eventweight=dsig(0)*weight

            ! in the last iteration :
            if (iteration.eq.N_iterations1) then

               ! pass event information to the histogram routine
               if (.not.extendForEvents .and. .not.multiUnweightingSecondRun)  &
                  call HistogramEvent(weight,dsig(0),dps, &
                  rand, p(0,1,1), x, v(0,1,1),  &
                  jets(0,1,1), n_jets(1),  &
                  leptons(0,1,1), n_leptons(1),  &
                  invisible(0,1,1),n_invisible(1), &
                  photons(0,1,1), n_photons(1), 0)

               call Reweight(eventweight, extendForEvents)        ! reweight

               if (eventweight.ne.0d0 .and. .not.multiUnweightingFirstRun) then ! event is to be kept
                  if((lha.or.hepmc).and..not.doNLO) call lha_put(p,v)
               endif
            endif
         endif   ! cutok
      else       ! dps > 0d0
         m2s(0) = 0d0
         dsig(0) = 0d0
      endif

      ! pass weight to monaco
      if (.not.extendForEvents .and. .not.multiUnweightingSecondRun)  &
         call monaco_put(rand, weight, dsig(0), &
              reweight_grid(p(0,1,1),v(0,1,1)))
      call born_extend_unweighting

   endsubroutine


   subroutine born_extend_unweighting
      ! run until the desired amount of unweighted events is achieved
      if (unweighting .and. iteration.eq.N_iterations1 .and. desiredEventCount.gt.0  &
         .and. i.eq.NCall_LO .and. .not.multiUnweightingFirstRun) then
      estimatedEventNumber = getEstimatedEvents()
      print*, ""
      if (multiChannelUnweighting) then
         desiredPlus5Sigma = int(dble(desiredEventCount)*contribXSFraction(contributionNumber)  &
            + 10.0*sqrt(dble(desiredEventCount)*contribXSFraction(contributionNumber)))
         if (partialUnweight) then
            print*, "Estimated number of partially unweighted events in this subchannel so far:" &
               , estimatedEventNumber
         else
            print*, "Estimated number of unweighted events in this subchannel so far:", estimatedEventNumber
         endif
      else
         desiredPlus5Sigma = int(dble(desiredEventCount) + 5.0*sqrt(dble(desiredEventCount)))
         if (partialUnweight) then
            print*, "Estimated number of partially unweighted events so far:", estimatedEventNumber
         else
            print*, "Estimated number of unweighted events so far:", estimatedEventNumber
         endif
      endif
      increaseNcallRatio = (dble(desiredPlus5Sigma) / dble(estimatedEventNumber))
      if (increaseNcallRatio.gt.2d0) increaseNcallRatio = 2d0
      if (increaseNcallRatio.lt.1.01d0) increaseNcallRatio = 1.01d0
      if (estimatedEventNumber.lt.desiredPlus5Sigma) then
         if ((dble(Ncall_LO) * increaseNcallRatio) / dble(Ncall) .lt. maxNcallFactor) then
            Ncall_LO = dble(Ncall_LO) * increaseNcallRatio
            if (.not.extendForEvents .and. .not.multiChannelUnweighting) then
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
            print*, "Extending this iteration by a factor of", increaseNcallRatio, &
               "to get the desired number of events."
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
            if (.not.multiChannelUnweighting) print*, "Writing out all events obtained so far..."
            print*, ""
         endif
      endif
   endif
end subroutine

function get_filename_in(sub_number_in, ps_number_in, prefix_grid, append_out) result(fname)
   implicit none
   integer, intent(in) :: sub_number_in
   integer, intent(in) :: ps_number_in
   logical, intent(in) :: prefix_grid
   logical, intent(in) :: append_out

   include 'VBFNLO/utilities/global.inc' !
   character*250 :: fname
   if (prefix_grid) then
      fname = 'GRID/'
   else
      if (append_out) then
         fname = ''
      else
         fname  = trim(path)//"/" 
      endif
   endif

   fname = trim(fname) // trim(gridname1(ps_number_in))
   if (sub_number_in > -1) then
      fname = trim(fname) //'_'//digit(sub_number_in)
   endif
   if (append_out) then
      fname = trim(fname) // '.out.'//digit(N_iterations1_orig)
   endif
end function


subroutine born_unweighting_runs_init
   ! settings concerning definite number of unweighted events
   if (multiUnweightingFirstRun) then
      multiUnweightingFirstRun = .false.      ! reset for 2nd run
      multiUnweightingSecondRun = .true.
      unwSumXS = 0.
      do k=1,numContribs
         unwSumXS = unwSumXS + contribXS(k)
      enddo
      do k=1,numContribs
         contribXSFraction(k) = contribXS(k) / unwSumXS
      enddo
      contributionNumber = 1
      N_iterations1_orig = N_iterations1
      N_iterations1 = 1           ! only one iteration in the 2nd run
      print*, ""
      print*, "****************************************************"
      print*, ""
      print*, "  Starting second run to get the unweighted events  "
      print*, ""
      print*, "****************************************************"
      print*, ""
   elseif (unweighting .and. desiredEventCount.gt.0 .and. (PS_Loops.gt.1 .or. Loops_sub_LO.gt.1)) then
      multiChannelUnweighting = .true.
      multiUnweightingFirstRun = .true.
      contributionNumber = 1
      maxNcallFactor = 129
      print*, ""
      print*, "********************************************"
      print*, ""
      print*, "  Starting first run to get cross sections  "
      print*, ""
      print*, "********************************************"
      print*, ""
   else
      maxNcallFactor = 65
   endif
end subroutine




subroutine fermion_loop
#ifdef WITH_MPI
      use vbfnlo_mpi
#endif

   !c This bit is here because we need the LO definitions for that
   !c and we save reinitializing things again

   ! should be calculate with LO PDFs
   if (.not.doBorn) call InitPDFs(0)
   ! and different PS
   infloops = .true.
   call InitPhaseSpace
   call InitGaugeTest

   call InitFLoop
   ! leave enough room for multiple loop contributions
   gnlo = 11

   do ps_number = 1, PS_Loops !multi-channel
      if (PS_Loops.gt.1 .and. lglobalprint) then
         write(*,*)" "
         write(*,*)"PHASE SPACE", ps_number
         write(*,*)" "
      endif


      do sub_number = 1, Loops_sub_LO ! multi-processes
         call cpu_time(time0)
         if (Loops_sub_LO.gt.1 .and. lglobalprint) then
            write(*,*)" "
            write(*,*)"Subprocess", sub_number
            write(*,*)" "
         endif

         if (Loops_sub_LO.gt.1) then
            file_name_in  = trim(path)//"/"// &
               gridname4(ps_number) &
               (1:INDEX(gridname4(ps_number),'  ')-1)
            file_name_in = trim(file_name_in) &
               //'_'//digit(sub_number)
         else
            file_name_in  = trim(path)//"/"// &
               gridname4(ps_number) &
               (1:INDEX(gridname4(ps_number),'  ')-1)
            file_name_in = trim(file_name_in)
         end if

         xsection=0d0
         sdev=0d0

         if (lreaddata) then
            call numdata(Ncall,PS_dimension,1,200+PS_number)
            N_iterations1 = 1
         endif

         do iteration = 1, N_iterations1

            ! Adjust statistics for different phase space
            call statistics_proc_FLoops(ps_number)

            if (Loops_sub_LO.gt.1) then
               file_name_out = 'GRID/'//gridname4(ps_number) &
                  (1:INDEX(gridname4(ps_number),'  ')-1)// &
                  '_'//digit(sub_number)//'.out.'//digit(iteration)
            else
               file_name_out = gridname4(ps_number) &
                  (1:INDEX(gridname4(ps_number),'  ')-1)// &
                  '.out.'//digit(iteration)
            end if

            if ( iteration .eq. 1 ) then
               call monaco_init( PS_dimension, Ncall)
               call monaco_read(file_name_in)
            else
               call monaco_init2( PS_dimension, Ncall )
            end if
            ! if (iteration.eq.N_iterations1) then
            !    call monran_set(1) ! save random number info, TODO: remove
            ! endif
#ifdef WITH_MPI
            call mpi_prepare_loop(iteration, Ncall)
            do i = 1, ncall_thisproc
#else
            do i = 1, NCall
#endif
               call fermion_phasespacepoint
            enddo
#ifdef WITH_MPI
            call mpi_after_loop(NCall)
#else
         call monaco_end_of_iteration()
#endif

#ifdef WITH_MPI
            if (vbfnlo_mpi_myid == vbfnlo_mpi_masterid) &
#endif
            call monaco_write( file_name_out )

         enddo

         call monaco_result(xsection, sdev, chi2)

         ! print out the result.
         if (lglobalprint) then
         print *,""
         if(Loops_sub_LO.eq.1) then
            print *," result (NLO fermion loop): ",xsection, &
               " +- ",sdev," fb  ",sdev/xsection*100d0,'%'
         else
            print *," result sub (NLO fermion loop): ",xsection, &
               " +- ",sdev," fb  ",sdev/xsection*100d0,'%'
         end if
         print *,""
         endif

         call cpu_time(time)
         call AddToXsec(xsection,sdev**2,time-time0,sub_number,ps_number,gnlo)

      end do ! end sub-processes

      if (lglobalprint) then
      if(Loops_sub_LO.gt.1) then
         print *,""
         print *," result (NLO fermion loop): ",final_xsec(0,ps_number,gnlo)," +- " &
            ,sqrt(final_sdev2(0,ps_number,gnlo))," fb  ",sqrt(final_sdev2(0,ps_number,gnlo))/ &
            final_xsec(0,ps_number,gnlo)*100d0,'%'
         print *," "
      endif
      endif
   enddo                  !multi-channel


   if (lglobalprint .and. ps_loops.gt.1) then
      print *,""
      print *," Total (multi-channel) result (NLO fermion loop): ", &
         final_xsec(0,0,gnlo)," +- ",sqrt(final_sdev2(0,0,gnlo))," fb  ", &
         sqrt(final_sdev2(0,0,gnlo))/final_xsec(0,0,gnlo)*100d0,'%'
      print *," "
   endif

   infloops = .false.
   call InitPhaseSpace
end subroutine


subroutine fermion_phasespacepoint
   ! retrieve a set of random numbers and the corresponding weight
   call monaco_get( rand, weight)

   if (lreaddata) call readdata(PS_dimension,rand,weight,1,m2s(0))
   ! transform the random numbers to a physical phase space point
   call phasespace(rand, p, x, v, ps_number, dps  )

   if ( dps .gt. 0d0 ) then

      if (alllep_flag) call lepton_gen

      ! define jets, leptons, photons
      call defineJets(p, v, 1, jets, n_jets)
      call DefineLeptons(v, 1, leptons, n_leptons)
      call DefineInvisible(v,1, invisible, n_invisible)
      call DefinePhotons(v, 1, photons, n_photons)

      ! in case of multi-channel: choose phase space:
      PS_OK = Choose_PS(ps_number,v, 1)
      ! apply cuts
      cut_ok = Cuts(jets, n_jets, leptons, n_leptons,  &
         invisible,n_invisible,photons, n_photons, 1)

      if (cut_ok .and. PS_ok) then
         ! calculate fac. and ren. scale for this phase space point (LO scale)
         call Scales(p, v, jets, n_jets, leptons,  &
            n_leptons, photons, n_photons, 1, 0)
         ! call the matrix element code and get an amplitude square
         if (.not.lreaddata) m2s(0) =  &
            Amplitude(rand, p, x, v, gnlo, ps_number)

         if ((iteration.eq.N_iterations1).and. &
            lwritedata) call writedata(PS_dimension, &
            rand,weight,1,m2s(0),200+PS_number)
      else
         dps = 0d0
         m2s(0) = 0d0
      end if
      dsig(0) = m2s(0)*dps

      if (cut_ok .and. PS_ok) then
         ! in the last iteration :
         if (iteration.eq.N_iterations1) then

            ! pass event information to the histogram routine
            call HistogramEvent(weight,dsig(0),dps, &
               rand, p(0,1,1), x, v(0,1,1),  &
               jets(0,1,1), n_jets(1),  &
               leptons(0,1,1), n_leptons(1),  &
               invisible(0,1,1),n_invisible(1), &
               photons(0,1,1), n_photons(1), -1)
         endif

      endif   ! cutok

   else       ! dps > 0d0
      m2s(0) = 0d0
      dsig(0) = 0d0
   end if

   ! pass weight to monaco
   call monaco_put(rand, weight, dsig(0), reweight_grid(p(0,1,1),v(0,1,1)))

end subroutine



subroutine virtuals_loop
#ifdef WITH_MPI
      use vbfnlo_mpi
#endif
    use stdoutredirect, only: normalstdout, redirectstdout
   ! Only one iteration, using the grid from the LO run.
   iteration = 1

   call InitPDFs(1)
   call InitGaugeTest

   do ps_number = 1,PS_Loops !multi-channel

      if (PS_Loops.gt.1 .and. lglobalprint) then
         print*," "
         print*,"PHASE SPACE ", ps_number
         print*," "
      endif
      ! Initialize to ZERO the instabilities and other variables

      call Init_Instabilities(ps_number)

      do sub_number = 1, Loops_sub_LO ! multi-processes

         if (Loops_sub_LO.gt.1 .and. lglobalprint) then
            write(*,*)" "
            write(*,*)"Subprocess", sub_number
            write(*,*)" "
         endif

         if (Loops_sub_LO.gt.1) then
            if (N_iterations1.eq.1) then
               file_name_in  = trim(path)//"/"//trim(gridname1(ps_number))
               file_name_in = trim(file_name_in)//'_'//digit(sub_number)
            else
               file_name_in  = 'GRID/'//trim(gridname1(ps_number))// &
                  '_'//digit(sub_number)//'.out.'//digit(N_iterations1)
            endif
         else
            if (N_iterations1.eq.1) then
               file_name_in  = trim(path)//"/"//trim(gridname1(ps_number))
               file_name_in = trim(file_name_in)
            else
               file_name_in  = trim(gridname1(ps_number)) // &
                  '.out.'//digit(N_iterations1)
            endif
         end if

         do gnlo = 1, NLO_Loops ! loop over different virtual contns

            if (lreaddata) then
               call numdata(Ncall,PS_dimension,1,10*gnlo+PS_number)
            endif

            ! Adjust statistics for different phase space regions
            call statistics_proc_NLO(ps_number)

            call cpu_time(time0)

            if (Loops_sub_LO.gt.1) then
               file_name_out = 'GRID/'//trim(gridname1(ps_number)) // &
                  '_'//digit(sub_number)//'.vout.'//digit(1)
            else
               file_name_out = trim(gridname1(ps_number)) // &
                  '.vout.'//digit(gnlo)
            end if

            call monaco_init( PS_dimension, Ncall)
            call monaco_read(file_name_in)
            ! call monran_set(2)

            select case(process)
            case(Hjj, Hjj_AA, Hjj_mu, Hjj_tau, Hjj_bbar)
                if (ewcor_switch .and. gnlo == 2) then
                    call redirectstdout()
                endif
            end select

#ifdef WITH_MPI
            call mpi_prepare_loop(iteration, ncall)
            do i = 1, ncall_thisproc
#else
            do i = 1, ncall
#endif
               call virtual_phasespacepoint
            enddo  ! Ncalls

            select case(process)
            case(Hjj, Hjj_AA, Hjj_mu, Hjj_tau, Hjj_bbar)
                if (ewcor_switch .and. gnlo == 2) then
                    call normalstdout()
                endif
            end select

#ifdef WITH_MPI
            call mpi_after_loop(ncall)
#else
         call monaco_end_of_iteration()
#endif

#ifdef WITH_MPI
            if (vbfnlo_mpi_myid == vbfnlo_mpi_masterid) &
#endif
            call monaco_write(file_name_out)

            call monaco_result(xsection, sdev, chi2)

            if (lglobalprint) then
               print *," "
               print *," result for virtual contribution, no. ",gnlo
               print *," sigma = ",xsection," +- ", sdev," fb  "
            endif

            ! printing stabilities
            if (lglobalprint) call Instabilities(ps_number)

            call cpu_time(time)
            call AddToXsec(xsection,sdev**2,time-time0,sub_number,ps_number,gnlo)

            if (lglobalprint) then
               print *," "
               print *,'total result for virtual contns so far '
               print *,' sigma = ',sum(final_xsec(sub_number,ps_number,1:NLO_Loops)), &
                       ' +- ',sqrt(sum(final_sdev2(sub_number,ps_number,1:NLO_Loops))), &
                       ' fb'
               print *,' '
            endif

         enddo   ! different virtual contributions

      enddo               ! end subprocesses

      if(Loops_sub_LO.gt.1)then
         print *,' total result for virtual PS ', ps_number
         print *,' sigma = ',final_xsec(0,ps_number,0),' +- ',sqrt(final_sdev2(0,ps_number,0)), &
            ' fb'
         print *,' '
      endif

   enddo                  !multi-channel

   if (lglobalprint .and. ps_loops.gt.1) then
      print *,""
      print *," Total (multi-channel) result (NLO virtual): " &
         ,final_xsec(0,0,0)," +- ",sqrt(final_sdev2(0,0,0))," fb  " &
         ,sqrt(final_sdev2(0,0,0))/final_xsec(0,0,0)*100d0,'%'
      print *,""
   endif
end subroutine

subroutine virtual_phasespacepoint
   call monaco_get( rand, weight)

   if (lreaddata) call readdata(PS_dimension,rand, &
      weight,1,m2s(0))

   call phasespace(rand, p, x, v, ps_number, dps  )

   if ( dps .gt. 0d0 ) then

      if (alllep_flag) call lepton_gen

      call defineJets(p, v, 1, jets, n_jets)
      call DefineLeptons(v, 1, leptons, n_leptons)
      call DefinePhotons(v, 1, photons, n_photons)
      call DefineInvisible(v,1, invisible, n_invisible)

      PS_OK = Choose_PS(PS_number, v, 1)

      cut_ok = Cuts(jets, n_jets, leptons, n_leptons,  &
         invisible,n_invisible,photons, n_photons, 1)

      if (cut_ok.and.ps_ok) then
         call Scales(p, v, jets, n_jets, leptons,  &
            n_leptons, photons, n_photons, 1, 1)
         if (.not.lreaddata) m2s(0) = Amplitude(rand,  &
            p, x, v, gnlo, ps_number)
         if (lwritedata) call writedata(PS_dimension, &
            rand,weight,1,m2s(0),10*gnlo+PS_number)
#ifdef WITH_LT
         call clearLT
#endif
      else
         dps = 0d0
         m2s(0) = 0d0
      end if

      dsig(0) = m2s(0)*dps

      call HistogramEvent(weight,dsig(0),dps, &
         rand, p(0,1,1), x, v(0,1,1),  &
         jets(0,1,1), n_jets(1),  &
         leptons(0,1,1), n_leptons(1),  &
         invisible(0,1,1),n_invisible(1),      &
         photons(0,1,1), n_photons(1), -1)

   else
      m2s(0) = 0d0
      dsig(0) = 0d0
   end if

   call monaco_put(rand, weight, dsig(0), reweight_grid(p(0,1,1),v(0,1,1)))
end subroutine


subroutine reals_loop
#ifdef WITH_MPI
      use vbfnlo_mpi
#endif
   if (.not. doVirtuals) call InitPDFs(1)
   call InitRealEmissions
   if (ldiptest) call InitDipTest

   do ps_number = 1,PS_Loops

      if (ps_loops.gt.1 .and. lglobalprint) then
         print*," "
         print*,"PHASE SPACE ", ps_number
         print*," "
      endif

      do sub_number = 1, Loops_sub_NLO ! multi-processes

         if (Loops_sub_NLO.gt.1 .and. lglobalprint) then
            write(*,*)" "
            write(*,*)"Subprocess", sub_number
            write(*,*)" "
         endif
         call cpu_time(time0)

         file_name_in  = trim(path)//"/"// &
            gridname2(ps_number) &
            (1:INDEX(gridname2(ps_number),'  ')-1)
         if (Loops_sub_NLO.gt.1) then
            file_name_in = trim(file_name_in) &
               //'_'//digit(sub_number)
         else
            file_name_in = trim(file_name_in)
         end if

         if (lreaddata) then
            call numdata(Ncall,PS_dimension,n_kin,100+PS_number)
            N_iterations2 = 1
         endif

         ! for histogram error estimation
         call InitRealHist(n_kin)

         DO iteration = 1, N_iterations2
            ! Control the statistics for different phase space
            call statistics_proc_Real(ps_number)

            if (Loops_sub_NLO.gt.1) then
               file_name_out = 'GRID/'//gridname2(ps_number) &
                  (1:INDEX(gridname2(ps_number),'  ')-1)// &
                  '_'//digit(sub_number)// &
                  '.out.'//digit(iteration)
            else
               file_name_out = gridname2(ps_number) &
                  (1:INDEX(gridname2(ps_number),'  ')-1)// &
                  '.out.'//digit(iteration)
            end if

            if ( iteration .eq. 1 ) then
               call monaco_init( PS_dimension, Ncall)
               call monaco_read(file_name_in)
            else
               call monaco_init2( PS_dimension, Ncall )
            end if

#ifdef WITH_MPI
            call mpi_prepare_loop(iteration, ncall)
            do i = 1, ncall_thisproc
#else
            do i = 1, Ncall
#endif
               call reals_phasespacepoint
            end do

#ifdef WITH_MPI
            call mpi_after_loop(ncall)
#else
         call monaco_end_of_iteration()
#endif

#ifdef WITH_MPI
            if (vbfnlo_mpi_myid == vbfnlo_mpi_masterid) &
#endif
               call monaco_write( file_name_out )

         enddo

         call monaco_result(xsection, sdev, chi2)

         if (lglobalprint) then
             print *,""
             if(Loops_sub_NLO.eq.1) then
                print *," result for real emissions "
             else
                print *," result subprocess for real emissions "
             endif
             print *," sigma = ",xsection," +- ",sdev," fb"
             print *," "
         endif

         call cpu_time(time)
         call AddToXsec(xsection,sdev**2,time-time0,sub_number,ps_number,NLO_Loops+1)

      enddo               ! endsuprocess

      if(Loops_sub_NLO.gt.1) then
         print *," result for real emissions "
         print *," sigma = ",final_xsec(0,ps_number,NLO_Loops+1)," +- ",sqrt(final_sdev2(0,ps_number,NLO_Loops+1)) &
            ," fb"
         print *," "
      endif

   enddo                  !multi-channel

   if (ldiptest) call FinishDipTest

   if (ps_loops.gt.1) then
      print *,""
      print *," Total (multi-channel) result (NLO real): " &
         ,final_xsec(0,0,NLO_Loops+1)," +- ",sqrt(final_sdev2(0,0,NLO_Loops+1))," fb  " &
         ,sqrt(final_sdev2(0,0,NLO_Loops+1))/final_xsec(0,0,NLO_Loops+1)*100d0,'%'
      print *,""
   endif
end subroutine

subroutine reals_phasespacepoint
   call monaco_get( rand, weight)

   if (lreaddata) call readdata(PS_dimension,rand, &
      weight,n_kin,m2s(0))

   call phasespace(rand, p, x, v, ps_number, dps  )

   if ( dps .gt. 0d0 ) then

      if (alllep_flag) call lepton_gen

      call ptilde(p(0,1,1),xuz(1,2),v(0,1,1))
      call defineJets(p, v, n_kin, jets, n_jets)
      call DefineLeptons(v, n_kin, leptons, n_leptons)
      call DefinePhotons(v, n_kin, photons, n_photons)
      call DefineInvisible(v,n_kin, invisible, n_invisible)

      ps_ok = Choose_PS(PS_number, v, 1)


      lokt(1) = Cuts(jets, n_jets, leptons, n_leptons, &
         invisible,n_invisible,photons,n_photons, 1) &
         .and.ps_ok

      cut_ok = lokt(1)

      do nd = 2, n_kin
         ps_ok = Choose_PS(PS_number, v, nd)
         lokt(nd) = Cuts(jets, n_jets, leptons, n_leptons, &
            invisible,n_invisible, &
            photons, n_photons, nd).and.ps_ok
         cut_ok = cut_ok .or. lokt(nd)
      enddo

      if (cut_ok) then
         call Scales(p, v, jets, n_jets, leptons,  &
            n_leptons, photons, n_photons, n_kin, 1)

         if (.not.lreaddata) call RE_Amplitude(rand,  &
            p, x, v, lokt, xuz, m2s,ps_number)
         if ((iteration.eq.N_iterations2).and. &
            lwritedata) call writedata(PS_dimension, &
            rand,weight,n_kin,m2s(1),100+PS_number)

         if (ldiptest .and. m2s(0).ne.0d0) call WriteDipTestEvent(p(0,1,1),m2s,dps,weight,lokt,iteration)

      else
         dps = 0d0
         do L = 0,n_kin
            m2s(L) = 0d0
         enddo
      end if

      do L = 0,n_kin
         dsig(L) = m2s(L)*dps
      enddo


      if (iteration.eq.N_iterations2) then

         realcont=.false.
         do L = 1,n_kin
            if (dsig(L).ne.0d0) then
               call HistogramEvent(weight,dsig(L),dps, &
                  rand, p(0,1,L), x, v(0,1,L),  &
                  jets(0,1,L), n_jets(L),  &
                  leptons(0,1,L), n_leptons(L), &
                  invisible(0,1,L),n_invisible(L),         &
                  photons(0,1,L), n_photons(L), L)
               realcont=.true.
            endif
         enddo
         if (realcont) call SaveRealHist

      endif

   else
      do L = 0,n_kin
         m2s(L) = 0d0
         dsig(L) = 0d0
      enddo
   end if

   call monaco_put(rand, weight, dsig(0), reweight_grid(p(0,1,1),v(0,1,1)))

endsubroutine


subroutine lha_output
#ifdef WITH_MPI
      use vbfnlo_mpi

   if ((lha .or. hepmc) .and. vbfnlo_mpi_nprocs > 1) then
       print*, "Event output is not available in combination with MPI."
       print*, "Please use MPI to generate an optimized grid and then generated events"
       print*, "using single-core runs."
       call FinalMPI()
       stop
   endif
#endif
   !...Les Houches interface
   noEventOutput = .false.
   if ((lha.or.hepmc).and..not.doNLO) then
      Select Case(procID)
      CASE(WMAJ, WPAJ, WPZJ, WMZJ, WPAAJ, WMAAJ, &
            WMAJLO, WPAJLO, WPZJLO, WMZJLO, &
            WPHJ_AA,WPHJ_mu,WPHJ_tau,WPHJ_bbar,WPHJ_WW,WPHJ_ZZ_ll,WPHJ_ZZ_lnu, & ! decay channels do not yet have
            WMHJ_AA,WMHJ_mu,WMHJ_tau,WMHJ_bbar,WMHJ_WW,WMHJ_ZZ_ll,WMHJ_ZZ_lnu, & ! a working event output
            WPH_AA,WPH_mu,WPH_tau,WPH_bbar,WPH_WW,WPH_ZZ_ll,WPH_ZZ_lnu, &
            WMH_AA,WMH_mu,WMH_tau,WMH_bbar,WMH_WW,WMH_ZZ_ll,WMH_ZZ_lnu &
            )
         ! in WH: Higgs decay LHA output not yet implemented
         lha=.false.
         hepmc=.false.
         noEventOutput = .true.
         print*,''
         print*,'Sorry, event file is not available for the process', &
            procID
         print*,''
      CASE(WPHJ, WMHJ, WPJ, WMJ)
         if (procid.ne.process) then
            ! in W(H)JjLO: color flow not yet implemented
            lha=.false.
            hepmc=.false.
            noEventOutput = .true.
            print*,''
            print*,'Sorry, event file is not available for the process', procID
            print*,''
         else
            call lha_file(-1)
         endif
      CASE DEFAULT
         call lha_file(-1)
      END SELECT
   endif
end subroutine


subroutine real_photons_loop
#ifdef WITH_MPI
      use vbfnlo_mpi
#endif
   !sophy: calculating the real photon emission.  Adapted from code written by
   ! Terrance Figy
   print *,""
   print *,"Starting real photon emissions"

   ! if (.not. ewcor_switch) call InitPDFs(1)

   call InitQEDEmissions

   do ps_number = 1,PS_Loops

      if (ps_loops.gt.1 .and. lglobalprint) then
         print*," "
         print*,"PHASE SPACE ", ps_number
         print*," "
      endif

      call cpu_time(time0)

      file_name_in  = trim(path)//"/"// &
         gridname3(ps_number) &
         (1:INDEX(gridname3(ps_number),'  ')-1)
      file_name_in = trim(file_name_in)

      ! for histogram error estimation
      call InitRealHist(n_qed)

      DO iteration = 1, N_iterations2
         Ncall = 2**( N_points(-1) - N_iterations2 + iteration )

         file_name_out = gridname3(ps_number) &
            (1:INDEX(gridname3(ps_number),'  ')-1)// &
            '.out.'//digit(iteration)

         if ( iteration .eq. 1 ) then
            call monaco_init( PS_dimension, Ncall )
            call monaco_read(file_name_in)
         else
            call monaco_init2( PS_dimension, Ncall )
         end if

#ifdef WITH_MPI
         call mpi_prepare_loop(iteration, ncall)
         do i = 1, ncall_thisproc
#else
         do i = 1, Ncall
#endif
         call real_photons_phasespacepoint
         end do

#ifdef WITH_MPI
         call mpi_after_loop(ncall)
#else
         call monaco_end_of_iteration()
#endif

#ifdef WITH_MPI
            if (vbfnlo_mpi_myid == vbfnlo_mpi_masterid) &
#endif
         call monaco_write( file_name_out )

      END DO

      call monaco_result(xsection, sdev, chi2)

      call cpu_time(time)
      call AddToXsec(xsection,sdev**2,time-time0,1,ps_number,NLO_Loops+2)

      ! print out the result.
      print *,""
      print *," result for QED real emissions "
      print *," sigma = ",xsection," +- ",sdev," fb"
      print *," "

   enddo                  !multi-channel

   if (ps_loops.gt.1) then
      print *,""
      print *," Total (multi-channel) result (QED real): " &
         ,final_xsec(0,0,NLO_loops+2)," +- ",sqrt(final_sdev2(0,0,NLO_loops+2))," fb  " &
         ,sqrt(final_sdev2(0,0,NLO_loops+2))/final_xsec(0,0,NLO_loops+2)*100d0,'%'
      print *,""
   endif

end subroutine

subroutine real_photons_phasespacepoint
   call monaco_get( rand, weight)

   call phasespace(rand, p, x, v, ps_number, dps  )

   if ( dps .gt. 0d0 ) then

      if (alllep_flag) call lepton_gen

      call ptilde_qed(p(0,1,1),xuz(1,2),v(0,1,1))
      call defineJets(p, v, n_qed, jets, n_jets)
      call DefineLeptons(v, n_qed, leptons, n_leptons)
      call DefinePhotons(v, n_qed, photons, n_photons)
      call DefineInvisible(v,n_qed, invisible, n_invisible)

      ps_ok = Choose_PS(PS_number, v, 1)

      lokt(1) = Cuts(jets, n_jets, leptons, n_leptons, &
         invisible,n_invisible, photons,n_photons, 1) &
         .and.ps_ok
      cut_ok = lokt(1)
      do nd = 2, n_qed
         ps_ok = Choose_PS(PS_number, v, nd)
         lokt(nd) = Cuts(jets, n_jets, leptons, n_leptons, &
            invisible,n_invisible, &
            photons, n_photons, nd) .and.ps_ok
         cut_ok = cut_ok .or. lokt(nd)
      enddo
      if (cut_ok) then
         call Scales(p, v, jets, n_jets, leptons,  &
            n_leptons, photons, n_photons, n_kin, 1)
         call QED_Amplitude(rand, p, x, v, lokt, xuz, m2s, &
            ps_number)
      else
         dps = 0d0
         m2s(0:n_qed) = 0d0
      end if

      do L = 0,n_qed
         dsig(L) = m2s(L)*dps
      enddo

      if (iteration.eq.N_iterations2) then

         realcont=.false.
         do L = 1,n_qed
            if (dsig(L).ne.0d0) then
               call HistogramEvent(weight,dsig(L),dps, &
                  rand, p(0,1,L), x, v(0,1,L),  &
                  jets(0,1,L), n_jets(L),  &
                  leptons(0,1,L), n_leptons(L), &
                  invisible(0,1,L),n_invisible(L),         &
                  photons(0,1,L), n_photons(L), L)
               realcont=.true.
            endif
         enddo
         if (realcont) call SaveRealHist

      endif

   else
      do L = 0,n_qed
         m2s(L) = 0d0
         dsig(L) = 0d0
      enddo
   end if


   call monaco_put(rand, weight, dsig(0), reweight_grid(p(0,1,1),v(0,1,1)))


end subroutine


subroutine final_lha_out
   !...Les Houches / HepMC interface
   if (lha.and..not.doNLO) then
      print*,''
      print*,'created LHA event file for the LO calculation : ', &
         lhafile
   endif
   if (hepmc.and..not.doNLO) then
      print*,''
      print*,'created HepMC event file for the LO calculation : ', &
         hepmcfile
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
   if ((ps_loops.gt.1 .or. Loops_sub_LO.gt.1) .and. .not.doNLO .and. (hepmc .or. lha)) then
      print*, ""
      print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' // &
      '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      print*, ""
      print*, "Important note:"
      if (Loops_sub_LO.eq.1) then
         print*, "The events of different phase spaces are written block-wise into the event file."
      else
         print*, "The events of different phase spaces and subprocesses are written block-wise into the event file."
      endif
      print*, "Therefore the event file should always be used completely, otherwise some parts of phase"
      print*, "space are underrepresentated. Using only parts of the event file gives only correct results"
      print*, "if the events are taken randomly from the whole file."
      print*, ""
      print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' // &
      '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      print*, ""
   endif

end subroutine

subroutine print_final_xsection
    character(len=*), parameter :: xsecformat = &
        '(" ", A23, f20.14, "  +-", f20.14, " fb", f20.14, " %")'

   if (doBorn) then
      print xsecformat,"TOTAL result (LO):    ", &
         final_xsec(0,0,-1),sqrt(final_sdev2(0,0,-1)), &
         sqrt(final_sdev2(0,0,-1))/final_xsec(0,0,-1)*100d0
   end if
   if ((doVirtuals.and.hasNLO) .or. ewcor_switch) then
      print xsecformat,"NLO virtual result:   ", &
         sum(final_xsec(0,0,1:NLO_Loops)),sqrt(sum(final_sdev2(0,0,1:NLO_Loops))), &
         sqrt(sum(final_sdev2(0,0,1:NLO_Loops)))/sum(final_xsec(0,0,1:NLO_Loops))*100d0
   end if
   if (doFLoops.and.hasFloops.and.(floops.ne.0)) then
      print xsecformat,"NLO fermion loop result: ", &
         final_xsec(0,0,11),sqrt(final_sdev2(0,0,11)), &
         sqrt(final_sdev2(0,0,11))/final_xsec(0,0,11)*100d0
   end if
   if (doEmissions.and.hasNLO) then
      print xsecformat,"QCD real emission:    ", &
         final_xsec(0,0,NLO_Loops+1),sqrt(final_sdev2(0,0,NLO_Loops+1)), &
         sqrt(final_sdev2(0,0,NLO_Loops+1))/final_xsec(0,0,NLO_Loops+1)*100d0
   end if
   if (ewcor_switch .and. (sector .ge. 3)) then
      print xsecformat,"QED real emission:    ", &
         final_xsec(0,0,NLO_Loops+2),sqrt(final_sdev2(0,0,NLO_Loops+2)), &
         sqrt(final_sdev2(0,0,NLO_Loops+2))/final_xsec(0,0,NLO_Loops+2)*100d0
   end if
   if ((doEmissions.and.hasNLO) .or. (ewcor_switch .and. (sector .ge. 3))) then
   print xsecformat,"TOTAL result (NLO):   ", &
      final_xsec(0,0,0),sqrt(final_sdev2(0,0,0)), &
      sqrt(final_sdev2(0,0,0))/final_xsec(0,0,0)*100d0
end if
   end subroutine

   subroutine print_virtual_instabilities
      if(doVirtuals) then
         select case (procID)
         case(QCDWPZjj, QCDWMZjj,QCDWPWPjj,QCDWMWMjj,QCDWPAjj,QCDWMAjj, QCDWPjj,QCDWMjj &
               ,QCDZAjj_l,QCDZAjj_n,QCDZZjj_ll,QCDZZjj_lnu,QCDZjj_l, QCDZjj_nu,QCDAAjj &
               )
            print*," Approx. error due to instable points: ",instab_error/final_xsec(0,0,0)*100d0," %"
#ifndef WITH_QUAD
            print*," This error can be reduced by using a compiler that supports quadruple precision"
#endif
         end select
      endif
   end subroutine


end module

