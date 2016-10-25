***************************************************************************
*
*  This is the VBFNLO main program. 
*
****************************************************************************

****************************************************************************
*  Last modified: 03.05.2011
****************************************************************************

      PROGRAM vbfnlo
      implicit none
#include "VBFNLO/utilities/global.inc"

      call ParseOptions
      ldoscales = .true.
      lglobalprint = .true.
      ldoblha = .false.
      call dovbfnlo

      end
