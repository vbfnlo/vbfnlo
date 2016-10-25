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

      call ParseOptions
      call doggflo

      end
