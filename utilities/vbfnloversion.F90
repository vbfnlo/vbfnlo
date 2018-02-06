module VBFNLOVersion

  integer vbfnloversionnumber
  character*50 vbfnloversionstring
  character*255 vbfnloreference(3)

contains

  subroutine setVersion
      implicit none

      integer progtype
      integer major, minor, patch, subnum
      character*20 substr, strsubnum
      character*12 strgithash
      character*2 strmajor, strminor, strpatch

      major = 3               ! first digit
      minor = 0               ! second digit
      patch = 0               ! third digit
      substr = "beta"         ! beta/rc
      subnum = 5              ! number for beta/rc

      vbfnloreference(1) =  'Comput. Phys. Commun. 180 (2009) 1661[arXiv:0811.4559 [hep-ph]]'
      vbfnloreference(2) = 'arXiV:1107.4038, arXiV:1207.4975, arXiv:1404.3940'
      vbfnloreference(3) = 'http://www.itp.kit.edu/vbfnlo'

      ! this gets filled by code in Makefile.am
#include "vbfnloversion.inc"
      if(len_trim(strgithash) .le. 4) then
        strgithash = ""
      endif

      vbfnloversionnumber=major*10000+minor*100+patch

      if (major.lt.10) then
        write(strmajor,'(I1)') major
      else
        write(strmajor,'(I2)') major
      endif
      if (minor.lt.10) then
        write(strminor,'(I1)') minor
      else
        write(strminor,'(I2)') minor
      endif
      if (patch.lt.10) then
        write(strpatch,'(I1)') patch
      else
        write(strpatch,'(I2)') patch
      endif
      if (trim(substr) .ne. "") then
        if (subnum.lt.10) then
          write(strsubnum,'(X,A,I1)') trim(substr), subnum
        else
          write(strsubnum,'(X,A,I2)') trim(substr), subnum
        endif
      else 
        strsubnum = ""
      endif

      vbfnloversionstring = &
        trim(strmajor)//"."// trim(strminor)//"."// trim(strpatch)// trim(strsubnum)// &
        " "// trim(strgithash) 
      end subroutine

      subroutine printVersion
          use readinput

#include "global.inc"
#include "process.inc"


        write(*,*)'  '
        write(*,*)'-----------------------------------------------------------------'
        write(*,*)'VBFNLO, version ', trim(vbfnloversionstring)
        write(*,*) trim(vbfnloreference(1))
        write(*,*) trim(vbfnloreference(2))
        write(*,*) trim(vbfnloreference(3))

        ! get process ID from vbfnlo.dat and 
        call loadfile("vbfnlo.dat",.false.)
        call read_int("PROCESS",procID,100,.true.)
        call read_logical("LOPROCESS_PLUS_JET",LOplusjet,.false.,.true.)
        call PrintProcInfo(.true.)
        call closefile

        write(*,*)'-----------------------------------------------------------------'
        write(*,*)'  '
        write(*,*)'  '

    end subroutine
end module

