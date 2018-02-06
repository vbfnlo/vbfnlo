!*************************************************************************  
!
!     Redirect stdout to file
!
!*************************************************************************  
!
!     Author: Michael Rauch <michael.rauch@kit.edu>
!     Initial version: Jul 2015
!     Last changed: Jul 2015 by Michael Rauch
!
!*************************************************************************  
!   LIST OF ALL FUNCTIONS AND SUBROUTINES IN THIS FILE:
!*************************************************************************
!     SUBROUTINE findstdout()
!     SUBROUTINE redirectstdout(outfile)
!     SUBROUTINE normalstdout()
!*************************************************************************  

! make this more flexible if we encounter an odd Fortran compiler at some point
MODULE stdoutredirect
!*************************************************************************  
!  Module containing definitions for redirecting stdout to a file
!*************************************************************************  
! make this more flexible if we encounter an odd Fortran compiler at some point
    integer, parameter :: stdout=6
    character(len=255) :: stdoutfile=""
    character(len=32) :: redirectfile="LoopTools.out"

contains
    SUBROUTINE redirectstdout()
!*************************************************************************
!     redirect stdout
!*************************************************************************
      implicit none

      logical :: lexist

      if ( .not. findstdout() ) return   ! no known way to get stdout back - we are screwed

      close(unit=stdout)
      inquire(file=redirectfile, exist=lexist)
      if (lexist) then
          open(unit=stdout,file=redirectfile,status="old",position='append')
      else
          open(unit=stdout,file=redirectfile,status="new")
      endif

  end subroutine

!*************************************************************************  

    SUBROUTINE normalstdout()
!*************************************************************************
!     stdout back to normal destination
!*************************************************************************
        implicit none

        if ( .not. findstdout() ) return   ! no known way to get stdout back 

        close(unit=stdout)
        open(unit=stdout,file=stdoutfile,status="old",position='append',err=99)
        return 

        ! retry without append
99      open(unit=stdout,file=stdoutfile,status="old")

    end subroutine

!*************************************************************************  

    logical FUNCTION findstdout()
!*************************************************************************
!     find how to get back stdout
!*************************************************************************
        implicit none

        logical :: lexist

        findstdout = .true.

        if (stdoutfile == "not found") then
            findstdout = .false.
            return
        endif
        if (stdoutfile /= "") return

        stdoutfile = "/dev/stdout"
        inquire(file=stdoutfile, exist=lexist)
        if ( lexist ) return 

        stdoutfile = "/proc/self/fd/1"
        inquire(file=stdoutfile, exist=lexist)
        if ( lexist ) return 

        stdoutfile = "not found"
        findstdout = .false.

    end function

!*************************************************************************  

END MODULE
