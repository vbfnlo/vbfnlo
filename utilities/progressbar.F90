!*************************************************************************  
!     This file contains a subroutine for printing a progress bar 
!     while integrating
!*************************************************************************  

      module progressbar
        implicit none

        logical :: progressprint

        integer*8, private :: npoints
        integer, private :: progress

        character(len=53), private :: header="|+----+----+----+----+----+----+----+----+----+----+|"
        character(len=58), private :: bar

        contains
          subroutine initprogress(n)
            implicit none
            integer*8, intent(in) :: n

            npoints = n
            progress = -1
            bar="|                                                   | ???%"
            write(unit=6,fmt="(a53)") header

          end subroutine

          subroutine printprogress(n)
            implicit none
            integer*8, intent(in) :: n
            integer currprogress
            integer k

            currprogress = (100*n)/npoints

            if (currprogress > progress) then
              progress = currprogress
              write(unit=bar(55:57),fmt="(i3)") progress
              do k=0,progress/2
                bar(2+k:2+k)="*"
              enddo
              ! print the progress bar.
              write(unit=6,fmt="(a1,a58)",advance="no") char(13), bar
              if (progress.eq.100) then
                write(unit=6,fmt=*)
              else
                flush(unit=6)
              endif
            endif

          end subroutine

      end module




