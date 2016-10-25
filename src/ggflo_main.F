***************************************************************************
*
*  This is the GGFLO main program. 
*
*  This is basically a duplicate of VBFNLO, but without the NLO part 
*  and some different initialization routines.
*
****************************************************************************

      PROGRAM ggflo
      implicit none
#include "VBFNLO/utilities/global.inc"

      call ParseOptions
      ldoscales = .true.
      lglobalprint = .true.
      ldoblha = .false.
      call doggflo

      end
