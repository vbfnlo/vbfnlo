!!! THIS IS A FORTRAN 90 FILE: CONTINUATION LINES AND COMMENTS HAVE A DIFFERENT SYNTAX !!!

! this file provides the three functions GetInputPath, GetInputFile and GetPdfsetsPath
! in a .F90 file in order to use free-form source format which allows longer @PATHS@ .


!*************************************************************************  
      character*250 FUNCTION GetInputPath()
          use readinput, only: inputpath
!*************************************************************************  
!     returns the path to input files
!     it is taken from the command line option input=my_input_path
!     or, if no command line option is given, from the 
!     environment variable VBF_INPUT_PATH
!*************************************************************************  
      implicit none
      character*250  s
      
! first, check for command line option
      if (len_trim(inputpath).gt.0) then
         GetInputPath = inputpath
         RETURN
      endif

! next, check for environment variable
      call GETENV( "VBF_INPUT_PATH" , s )
      if (len_trim(s).gt.0) then
         GetInputPath = s
         RETURN
      endif

! else, current directory
      GetInputPath = "."      
      END


!*************************************************************************  
      SUBROUTINE GetInputFile(inunit, filename, *)
          use readinput, only: inputpath
!*************************************************************************  
!     returns the path to input files
!     it is taken from the command line option input=my_input_path
!     or, if no command line option is given, from the 
!     environment variable VBF_INPUT_PATH
!*************************************************************************  
      implicit none
      integer inunit
      character*(*) filename

      character*250 pathname, s
      integer ier
      
! check with same order as GetInputPath(), 
! but continue to further options if not found 

! first, check for command line option
      if (len_trim(inputpath).gt.0) then
         pathname = trim(inputpath)//"/"//filename
         open (unit=inunit, file=pathname, status="old", iostat=ier)
         if (ier.eq.0) return
      endif

! next, check for environment variable
      call GETENV( "VBF_INPUT_PATH" , s )
      if (len_trim(s).gt.0) then
         pathname = trim(s)//"/"//filename
         open (unit=inunit, file=pathname, status="old", iostat=ier)
         if (ier.eq.0) return
      endif

! next, current working directory
      pathname = "./"//filename
      open (unit=inunit, file=pathname, status="old", iostat=ier)
      if (ier.eq.0) return

! compile time installation directory 
      pathname = INPUTPATH //"/"//filename
      open (unit=inunit, file=pathname, status="old", iostat=ier)
      if (ier.eq.0) return

      RETURN 1
      END


!*************************************************************************  
      character*250 FUNCTION GetPdfsetsPath()
          use readinput, only: pdfpath
!*************************************************************************  
!     Returns the path to the VBFNLO-internal pdfsets files.
!     The default is hardwired to the share directory
!     of the installation.
!*************************************************************************  
      implicit none

      character*250  s

! first, check for command line option
      if (len_trim(pdfpath) > 0) then
         GetPdfsetsPath = pdfpath
         RETURN
      endif

! check for environment variable
      call GETENV("VBF_PDFSETS_PATH", s)
      if (len_trim(s).gt.0) then
         GetPdfsetsPath = s
         RETURN
      endif

! default to installation, set when building
      GetPdfsetsPath = PDFSETSPATH
      END


