! This file can be used to generate a set of phase space points in VBFNLO, write them to a file and read them in again.
! The generated files can be used to compare amplitude or cross section level results to other MC programs for specific
! phase space points, which can be useful for debugging or accuracy comparisons.
!
! Usage (check amplitudes/hjjj/m2s_qqh3j.F for an example): 
! 1) Set PSwrite = .true. to generate a set of Phase space points in PSfilename. The easiest way to get the correct 
!    set of momenta for some process is to insert the following in the beginning of the respective amplitudes/YYY/m2s_XXX.F file:
! right after "subroutine m2s_XXX":
!
!        use fixedPS
!
! after the declaration of variables:
!
!        call writePSpoint(p,v)
!
! 2) Read in phase space points from file setting PSread = .true, PSwrite = .false. and by adding:
!
!        call readPSpoint(p,v)
!
! 3) At the very end of a routine, write the results wanted to Resfilename by setting PSread=.true., Reswrite =.true. and  with:
!
!      if(Reswrite) write(ResUnit,*) m2s
!
! If you wish to add other files, e.g. to read in data from another program, check the fixedPS bits in lib/vbfnlo_main.F90 and apply
! changes accordingly.
        
        module fixedPS
          use globalvars
          implicit none
          
          integer, parameter :: PSpoints = 1000 ! Number of PS points to be generated or read in
          integer, save :: printCount = 1  
          integer, save :: readCount = 1  
          
          character*250, save :: PSfilename = 'PS.input'
          integer, parameter :: PSunit = 9999
          character*250, save :: Resfilename = 'res.output'
          integer, parameter :: Resunit = 8888  
          
          ! Switches for reading and writing 
          logical, parameter :: PSwrite   =.false.
          logical, parameter :: PSread    =.false.
          logical, parameter :: Reswrite  =.false.
contains

          subroutine startPSwrite
            implicit none

            open(PSunit, file=PSfilename, status='replace')
            print*, "-----------------------------------------------"
            print*, "Phase space points will be generated and saved to ", PSfilename
            print*, "Number of points to be generated: ", PSpoints
            print*, "-----------------------------------------------"
          end subroutine

          subroutine startPSread
            implicit none

            open(PSunit, file=PSfilename, status='old', action='read')
            print*, "-----------------------------------------------"
            print*, "Fixed set of phase space points will be read in from: ", PSfilename
            if(Reswrite) then
              open(Resunit, file=Resfilename, status='replace')
              print*, "VBFNLO results for fixed phase space points will be written to:  ", Resfilename
            endif
            print*, "-----------------------------------------------"
          end subroutine

          subroutine stopPSwrite
            implicit none 

            if((printcount.gt.(PSpoints))) then
              print*, ""
              print*, "-----------------------------------------------"
              print*, "Generation of fixed phase space points finished."
              print*, "Number of points generated: ", PSpoints
              print*, "Output file: ", PSfilename
              print*, "-----------------------------------------------"
              stop
            endif
          end subroutine
             
          subroutine writePSpoint(p,v)
            implicit none
#include "VBFNLO/utilities/global.inc"
            real*8  p(0:3,max_p,max_kin), v(0:3,max_v,max_kin)
            integer mu, i

            if (PSwrite.and.printCount.le.PSpoints) then
             do i=1,n_p
               write(PSunit,*) (p(mu,i,1),mu=0,3)
             enddo
             do i=1,n_v
               write(PSunit,*) (v(mu,i,1),mu=0,3)
             enddo
             printCount = printCount + 1
            endif
          end subroutine

          subroutine readPSpoint(p,v)
            implicit none
#include "VBFNLO/utilities/global.inc"
            real*8  p(0:3,max_p,max_kin), v(0:3,max_v,max_kin)
            integer mu, i

            if (PSread) then
             readCount = readCount + 1
             if(readCount.gt.PSpoints+1) stop
             do i=1,n_p
               read(PSunit,*) (p(mu,i,1),mu=0,3)
             enddo
             do i=1,n_v
               read(PSunit,*) (v(mu,i,1),mu=0,3)
             enddo
            endif
          end subroutine

        end module fixedPS

