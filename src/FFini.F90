!** We start with the FF initialisation
      subroutine initFF
       use globalvars, only: lglobalprint
      implicit none

      if (lglobalprint) then
         write(*,*)'  '
         write(*,*)'Initialising FF ..'
         write(*,*)'  '
      endif

      call vbfffini
#ifdef WITH_QUAD
      call vbfquadffini
#endif
      end


!** This subroutine closes FF and checks errors
      subroutine closeFF
       use globalvars, only: lglobalprint
          
      implicit none

      if (lglobalprint) then
         write(*,*)'  '
         write(*,*)'Summary for FF ..'
         write(*,*)'  '
      endif

      call vbfffexi
#ifdef WITH_QUAD
      call vbfquadffexi
#endif
      end
