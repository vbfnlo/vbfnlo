      subroutine mtime
      implicit none
      character(len=13) cd0, ct0, cz0
      integer iv0(8)
      CALL DATE_AND_TIME(cd0,ct0,cz0,iv0)
      print*,""
      print*,"TIME=",iv0
      return
      end

