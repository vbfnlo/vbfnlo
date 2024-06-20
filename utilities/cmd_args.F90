
module cmd_args

      contains
!*************************************************************************  
      subroutine ParseOptions
!*************************************************************************  
! go through list of command-line options
! handle everything and warn about all unknown options
!*************************************************************************  
        use VBFNLOVersion, only: setVersion, vbfnloversionstring
        use readinput, only: inputpath, pdfpath
        use globalvars, only: seed
      implicit none
      character*250 s, sorig
      integer narg,i,j,k,le
      
#include "global.inc"

      seed=0
      inputpath=""
      pdfpath=""

! go through list
      narg = IARGC()
      if (narg.gt.0)then
         do i=1,narg
            call GETARG(i , s)
            sorig=s
! upper case everything in front of a "="
            le = index(s,"=")
! if there is no "=" upper case everything
            if (le.eq.0) le = len_trim(s)+1
            do j=1,le-1
               k=ichar(s(j:j))
               if(k.ge.97.and.k.le.122) then
                  k=ichar(s(j:j))-32   
                  s(j:j)=char(k)
               endif
            enddo
! remove leading dashes
            do while (s(1:1).eq.'-')
              s = s(2:len(s))
            enddo
            le = index(s,"=")
            if (le.eq.0) le=len_trim(s)+1
            if (s(1:le-1).eq."INPUT") then
                inputpath=s(le+1:len(s))
            elseif (s(1:le-1).eq."PDFSETS") then
                pdfpath=s(le+1:len(s))
            elseif (s(1:le-1).eq."SEED") then
                read(s(le+1:len(s)), *) seed
            elseif (s(1:le-1).eq."VERSION") then
                call setVersion
                print *, "VBFNLO ", trim(vbfnloversionstring)
                stop
            else
                print *, "Warning: Unknown option ", trim(sorig)
            endif
         enddo
      endif

      END subroutine

!*************************************************************************  
      function VBFNLOCapabilities(option)
        use VBFNLOVersion, only: vbfnloversionnumber
!*************************************************************************  
! test all options VBFNLO can handle
!*************************************************************************  
      implicit none

#include "global.inc"

      integer VBFNLOCapabilities
      character(*) option


! general
      IF (trim(option).eq."VERSION") then
      VBFNLOCapabilities = vbfnloversionnumber

! processes
#ifdef WITH_VBF
      ELSEIF (trim(option).eq."VBF") then
      VBFNLOCapabilities = 1
#endif
#ifdef WITH_HJJJ
      ELSEIF (trim(option).eq."HJJJ") then
      VBFNLOCapabilities = 1
#endif
#ifdef WITH_DIBOSON
      ELSEIF (trim(option).eq."DIBOSON") then
      VBFNLOCapabilities = 1
#endif
#ifdef WITH_DIBOSONJET
      ELSEIF (trim(option).eq."DIBOSONJET") then
      VBFNLOCapabilities = 1
#endif
#ifdef WITH_TRIBOSON
      ELSEIF (trim(option).eq."TRIBOSON") then
      VBFNLOCapabilities = 1
#endif
#ifdef WITH_TRIBOSONJET
      ELSEIF (trim(option).eq."TRIBOSONJET") then
      VBFNLOCapabilities = 1
#endif
#ifdef WITH_GGF
      ELSEIF (trim(option).eq."GGF") then
      VBFNLOCapabilities = 1
#endif
#ifdef WITH_QCDV
      ELSEIF (trim(option).eq."QCDV") then
      VBFNLOCapabilities = 1
#endif
#ifdef WITH_QCDVV
      ELSEIF (trim(option).eq."QCDVV") then
      VBFNLOCapabilities = 1
#endif

! options
#ifdef WITH_KK
      ELSEIF (trim(option).eq."KK") then
      VBFNLOCapabilities = 1
#endif
#ifdef WITH_NLO
      ELSEIF (trim(option).eq."NLO") then
      VBFNLOCapabilities = 1
#endif

! programs
#ifdef WITH_FH
      ELSEIF (trim(option).eq."FH") then
      VBFNLOCapabilities = FHVERSION
#endif
#ifdef WITH_HEPMC
      ELSEIF (trim(option).eq."HEPMC") then
      VBFNLOCapabilities = 1
#endif
#ifdef WITH_LHA
      ELSEIF (trim(option).eq."LHA") then
      VBFNLOCapabilities = 1
#endif
#ifdef WITH_LT
      ELSEIF (trim(option).eq."LT") then
      VBFNLOCapabilities = 1
#endif
#ifdef WITH_ROOT
      ELSEIF (trim(option).eq."ROOT") then
      VBFNLOCapabilities = 1
#endif

      ELSE
        VBFNLOCapabilities = 0

      ENDIF

      end function

      end module
