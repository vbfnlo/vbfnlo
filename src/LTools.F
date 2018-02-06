********************************************************************************
********************************************************************************
*** LTools.F                                                                 ***
*** 9th February 2011                                                        ***
***                                                                          ***
*** This folder contains the basic LoopTools: initialisation, closing etc.   ***
*** used by the main vbfnlo and ggflo codes                                  ***
***                                                                          ***
********************************************************************************
********************************************************************************

*** We start with the LoopTools initialisation

      subroutine initLT

      implicit none

#include "looptools.h"


      write(*,*)'  '
      write(*,*)'Initialising LoopTools ..'
      write(*,*)'  '

#ifdef LT_NEW
      call ltini
#endif
#ifdef LT_OLD
      call ffini
#endif


      end


********************************************************************************
********************************************************************************

*** This subroutine closes LoopTools and checks errors

      subroutine closeLT

      implicit none

#include "looptools.h"


#ifdef LT_NEW
      call ltexi
#endif
#ifdef LT_OLD
      call ffexi
#endif


      end


********************************************************************************
********************************************************************************

*** This subroutine clears the LoopTools cache

      subroutine clearLT

      implicit none

#include "looptools.h"


      call clearcache


      end


********************************************************************************
