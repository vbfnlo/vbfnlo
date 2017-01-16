************************************************************************************
* This file contains routines for checking the CS dipoles of the NLO calculations. *
************************************************************************************


      Subroutine InitDipTest()

      implicit none

#include "global.inc"
#include "cuts.inc"

      ! open files, write gnuplot command file

      CALL system("mkdir -p Dipole")

      open(90,FILE='./Dipole/dipoletest.gp',STATUS='REPLACE')
      write (90,"(A)") 'set output "hist_dipoletest.ps"'
      write (90,"(A)") "set terminal postscript color"
      write (90,"(A)") "set style data points"
      write (90,"(A)") "set key off"
      write (90,"(A)") "set logscale x"
      write (90,"(A)") "set logscale y"
      write (90,"(A)") "set yrange[1e-7:]"
      write (90,"(A)") ""

      open(90,FILE='./Dipole/p1p3.dat',STATUS='REPLACE')
      write (90,"(A)") 'set title "(real-dipoles)/real / dm_p1_p3"'
      write (90,"(A)") 'plot "p1p3.dat" using 1:(abs($2/$3))'
      open(92,FILE='./Dipole/p2p3.dat',STATUS='REPLACE')
      write (90,"(A)") 'set title "(real-dipoles)/real / dm_p2_p3"'
      write (90,"(A)") 'plot "p2p3.dat" using 1:(abs($2/$3))'
      open(93,FILE='./Dipole/soft3.dat',STATUS='REPLACE')
      write (90,"(A)") 'set title "(real-dipoles)/real / dE_parton_3"'
      write (90,"(A)") 'plot "soft3.dat" using 1:(abs($2/$3))'
      if (n_p-nx .gt. 1) then
         open(94,FILE='./Dipole/p1p4.dat',STATUS='REPLACE')
         write (90,"(A)") 'set title "(real-dipoles)/real / dm_p1_p4"'
         write (90,"(A)") 'plot "p1p4.dat" using 1:(abs($2/$3))'
         open(95,FILE='./Dipole/p2p4.dat',STATUS='REPLACE')
         write (90,"(A)") 'set title "(real-dipoles)/real / dm_p2_p4"'
         write (90,"(A)") 'plot "p2p4.dat" using 1:(abs($2/$3))'
         open(96,FILE='./Dipole/p3p4.dat',STATUS='REPLACE')
         write (90,"(A)") 'set title "(real-dipoles)/real / dm_p3_p4"'
         write (90,"(A)") 'plot "p3p4.dat" using 1:(abs($2/$3))'
         open(97,FILE='./Dipole/soft4.dat',STATUS='REPLACE')
         write (90,"(A)") 'set title "(real-dipoles)/real / dE_parton_4"'
         write (90,"(A)") 'plot "soft4.dat" using 1:(abs($2/$3))'
      endif
      if (n_p-nx .gt. 2) then
         open(98,FILE='./Dipole/p1p5.dat',STATUS='REPLACE')
         write (90,"(A)") 'set title "(real-dipoles)/real / dm_p1_p5"'
         write (90,"(A)") 'plot "p1p5.dat" using 1:(abs($2/$3))'
         open(99,FILE='./Dipole/p2p5.dat',STATUS='REPLACE')
         write (90,"(A)") 'set title "(real-dipoles)/real / dm_p2_p5"'
         write (90,"(A)") 'plot "p2p5.dat" using 1:(abs($2/$3))'
         open(100,FILE='./Dipole/p3p5.dat',STATUS='REPLACE')
         write (90,"(A)") 'set title "(real-dipoles)/real / dm_p3_p5"'
         write (90,"(A)") 'plot "p3p5.dat" using 1:(abs($2/$3))'
         open(101,FILE='./Dipole/p4p5.dat',STATUS='REPLACE')
         write (90,"(A)") 'set title "(real-dipoles)/real / dm_p4_p5"'
         write (90,"(A)") 'plot "p4p5.dat" using 1:(abs($2/$3))'
         open(102,FILE='./Dipole/soft5.dat',STATUS='REPLACE')
         write (90,"(A)") 'set title "(real-dipoles)/real / dE_parton_5"'
         write (90,"(A)") 'plot "soft5.dat" using 1:(abs($2/$3))'
      endif
      if (n_p-nx .gt. 3) then
         open(103,FILE='./Dipole/p1p6.dat',STATUS='REPLACE')
         write (90,"(A)") 'set title "(real-dipoles)/real / dm_p1_p6"'
         write (90,"(A)") 'plot "p1p6.dat" using 1:(abs($2/$3))'
         open(104,FILE='./Dipole/p2p6.dat',STATUS='REPLACE')
         write (90,"(A)") 'set title "(real-dipoles)/real / dm_p2_p6"'
         write (90,"(A)") 'plot "p2p6.dat" using 1:(abs($2/$3))'
         open(105,FILE='./Dipole/p3p6.dat',STATUS='REPLACE')
         write (90,"(A)") 'set title "(real-dipoles)/real / dm_p3_p6"'
         write (90,"(A)") 'plot "p3p6.dat" using 1:(abs($2/$3))'
         open(106,FILE='./Dipole/p4p6.dat',STATUS='REPLACE')
         write (90,"(A)") 'set title "(real-dipoles)/real / dm_p4_p6"'
         write (90,"(A)") 'plot "p4p6.dat" using 1:(abs($2/$3))'
         open(107,FILE='./Dipole/p5p6.dat',STATUS='REPLACE')
         write (90,"(A)") 'set title "(real-dipoles)/real / dm_p5_p6"'
         write (90,"(A)") 'plot "p5p6.dat" using 1:(abs($2/$3))'
         open(108,FILE='./Dipole/soft6.dat',STATUS='REPLACE')
         write (90,"(A)") 'set title "(real-dipoles)/real / dE_parton_6"'
         write (90,"(A)") 'plot "soft6.dat" using 1:(abs($2/$3))'
      endif

      close(90)

      return
      end


      Subroutine WriteDipTestEvent(p,m2s,dps,weight,lokt, iteration)

      implicit none

#include "global.inc"
#include "cuts.inc"

      ! write histogram data

      real*8 p(0:3,max_p) 
      real*8 weight, dps, m2s(0:max_kin)
      logical lokt(max_kin), loktdip
      integer i
      integer iteration
            
      real*8 mjj2
      external mjj2

      logical alliterations               ! should all iterations of the real emission calculation be used or only the last one?
      parameter ( alliterations = .false. )
      
      double precision cutoff_invmass, cutoff_soft
      cutoff_invmass = 1d0                ! max value for invariant mass of the collinear partons that is logged
      cutoff_soft = 1d0                   ! max value for energy of the soft partons that is logged


      ! values written out to file:
      ! * inv mass / energy
      ! * sum of M^2 for real emission kinematics and all dipoles
      ! * M^2 for real emission kinematics
      ! * sum of M^2 for real emission kinematics and all dipoles
      ! * list of booleans: for all kinematics: are cuts fulfilled?


      if (.not. alliterations .and. iteration.lt.N_iterations2) return

      loktdip=.false.
      do i = 2, n_kin
            loktdip = loktdip .or. lokt(i)
      enddo
      if ((.not.loktdip) .or. .not.lokt(1)) return      ! the RE or all of the dipoles missed the cuts
                                                        ! this line can be commented out if events in which the dipoles or the RE
                                                        ! miss the cuts are of interest as well...


            if(mjj2(p(0,1),p(0,3)).lt.cutoff_invmass) write(91,*) sqrt(mjj2(p(0,1),p(0,3))), m2s(0),m2s(1),
     &            m2s(0)*dps*weight,(lokt(i), i=1,n_kin)
            if(mjj2(p(0,2),p(0,3)).lt.cutoff_invmass) write(92,*) sqrt(mjj2(p(0,2),p(0,3))), m2s(0),m2s(1),
     &            m2s(0)*dps*weight,(lokt(i), i=1,n_kin)
            if(p(0,3).lt.cutoff_soft) write(93,*) p(0,3), m2s(0),m2s(1),m2s(0)*dps*weight,(lokt(i), i=1,n_kin)

            if (n_p-nx .gt. 1) then
                  if(mjj2(p(0,1),p(0,4)).lt.cutoff_invmass) write(94,*) sqrt(mjj2(p(0,1),p(0,4))), m2s(0),m2s(1),
     &                  m2s(0)*dps*weight,(lokt(i), i=1,n_kin)
                  if(mjj2(p(0,2),p(0,4)).lt.cutoff_invmass) write(95,*) sqrt(mjj2(p(0,2),p(0,4))), m2s(0),m2s(1),
     &                  m2s(0)*dps*weight,(lokt(i), i=1,n_kin)
                  if(mjj2(p(0,3),p(0,4)).lt.cutoff_invmass) write(96,*) sqrt(mjj2(p(0,3),p(0,4))), m2s(0),m2s(1),
     &                  m2s(0)*dps*weight,(lokt(i), i=1,n_kin)
                  if(p(0,4).lt.cutoff_soft) write(97,*) p(0,4), m2s(0),m2s(1),m2s(0)*dps*weight,(lokt(i), i=1,n_kin)
            endif

            if (n_p-nx .gt. 2) then
                  if(mjj2(p(0,1),p(0,5)).lt.cutoff_invmass) write(98,*) sqrt(mjj2(p(0,1),p(0,5))), m2s(0),m2s(1),
     &                  m2s(0)*dps*weight,(lokt(i), i=1,n_kin)
                  if(mjj2(p(0,2),p(0,5)).lt.cutoff_invmass) write(99,*) sqrt(mjj2(p(0,2),p(0,5))), m2s(0),m2s(1),
     &                  m2s(0)*dps*weight,(lokt(i), i=1,n_kin)
                  if(mjj2(p(0,3),p(0,5)).lt.cutoff_invmass) write(100,*) sqrt(mjj2(p(0,3),p(0,5))), m2s(0),m2s(1),
     &                  m2s(0)*dps*weight,(lokt(i), i=1,n_kin)
                  if(mjj2(p(0,4),p(0,5)).lt.cutoff_invmass) write(101,*) sqrt(mjj2(p(0,4),p(0,5))), m2s(0),m2s(1),
     &                  m2s(0)*dps*weight,(lokt(i), i=1,n_kin)
                  if(p(0,5).lt.cutoff_soft) write(102,*)  p(0,5), m2s(0),m2s(1),m2s(0)*dps*weight,(lokt(i), i=1,n_kin)
            endif

            if (n_p-nx .gt. 3) then       ! Hjjj
                  if(mjj2(p(0,1),p(0,6)).lt.cutoff_invmass) write(103,*) sqrt(mjj2(p(0,1),p(0,6))), m2s(0),m2s(1),
     &                  m2s(0)*dps*weight,(lokt(i), i=1,n_kin)
                  if(mjj2(p(0,2),p(0,6)).lt.cutoff_invmass) write(104,*) sqrt(mjj2(p(0,2),p(0,6))), m2s(0),m2s(1),
     &                  m2s(0)*dps*weight,(lokt(i), i=1,n_kin)
                  if(mjj2(p(0,3),p(0,6)).lt.cutoff_invmass) write(105,*) sqrt(mjj2(p(0,3),p(0,6))), m2s(0),m2s(1),
     &                  m2s(0)*dps*weight,(lokt(i), i=1,n_kin)
                  if(mjj2(p(0,4),p(0,6)).lt.cutoff_invmass) write(106,*) sqrt(mjj2(p(0,4),p(0,6))), m2s(0),m2s(1),
     &                  m2s(0)*dps*weight,(lokt(i), i=1,n_kin)
                  if(mjj2(p(0,5),p(0,6)).lt.cutoff_invmass) write(107,*) sqrt(mjj2(p(0,5),p(0,6))), m2s(0),m2s(1),
     &                  m2s(0)*dps*weight,(lokt(i), i=1,n_kin)
                  if(p(0,6).lt.cutoff_soft) write(108,*) p(0,6), m2s(0),m2s(1),m2s(0)*dps*weight,(lokt(i), i=1,n_kin)
            endif

      
      return
      end


      Subroutine FinishDipTest()

      implicit none

#include "global.inc"
#include "cuts.inc"

      integer maxfilenumber, i

      ! close files

      maxfilenumber=93
      if (n_p-nx .gt. 1) maxfilenumber=97
      if (n_p-nx .gt. 2) maxfilenumber=102
      if (n_p-nx .gt. 3) maxfilenumber=108

      do i = 91,maxfilenumber
         close(i)
      enddo
      print*,'Createted Dipole test files'

      return
      end

