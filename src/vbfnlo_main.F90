!**************************************************************************
!  This is the VBFNLO main program. 
!***************************************************************************

      PROGRAM vbfnlo
         use vbfnlo_main, only: dovbfnlo
         use globalvars, only: ldoscales, ldoblha
         use cmd_args, only: ParseOptions
      implicit none

      call ParseOptions
      ldoscales = .true.
      ldoblha = .false.
      call dovbfnlo

      end
