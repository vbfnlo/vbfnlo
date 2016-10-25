********************************************************************************
********************************************************************************
*** FFini.F                                                                  ***
*** 23rd May 2012                                                            ***
***                                                                          ***
*** This folder contains the basic FF: initialisation, closing etc.          ***
*** used by the main vbfnlo and ggflo codes                                  ***
***                                                                          ***
********************************************************************************
********************************************************************************

*** We start with the FF initialisation

      subroutine initFF

      implicit none

      write(*,*)'  '
      write(*,*)'Initialising FF ..'
      write(*,*)'  '

      call vbfffini
#ifdef WITH_QUAD
      call vbfquadffini
#endif

      end



********************************************************************************
********************************************************************************

*** This subroutine closes FF and checks errors

      subroutine closeFF

      implicit none

      write(*,*)'  '
      write(*,*)'Summary for FF ..'
      write(*,*)'  '

      call vbfffexi
#ifdef WITH_QUAD
      call vbfquadffexi
#endif

      end



********************************************************************************
