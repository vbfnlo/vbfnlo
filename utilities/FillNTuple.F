      SUBROUTINE FillNTuple(weight, dsig, p, v, NLO)
      implicit none

#include "global.inc"
      integer*8 nevent
      common/Nevent/nevent

      real*8 weight, dsig
      integer NLO
      real*8 p(0:3,max_p), v(0:3,max_v)
      real*8 pfill(4*max_p+4*max_v+2)
      real*8 dw

      integer imax_vmom, imax_pmom, k, i, l

      character*100 GetInputPath
      character*100 path,filename
      external GetInputPath

      if (NLO .lt. -1) return


      if (dsig.eq.0d0) return

      dw = dsig * weight
c Here we write out n tuples
c      write(*,*) "nevent=",nevent
c Determine the nonzero entries in p and v.

      imax_vmom = 0
      do k = 1, max_v
        if((v(0,k)+v(1,k)+v(2,k)+v(3,k)).ne.0d0) imax_vmom = imax_vmom + 1
      enddo

      imax_pmom = 2
      do k = 3, max_p
        if((p(0,k)+p(1,k)+p(2,k)+p(3,k)).ne.0d0) imax_pmom = imax_pmom + 1
      enddo

c      write(*,*) 'imax_pmom=',imax_pmom
c      write(*,*) 'imax_vmom=',imax_vmom
c      write(*,*) "nevent=",nevent

c Do we need the event number?

c Fill the partonic momenta.
      do i = 1, imax_pmom
        pfill(4*(i-1)+1) = p(1,i) ! px_i
        pfill(4*(i-1)+2) = p(2,i) ! py_i
        pfill(4*(i-1)+3) = p(3,i) ! pz_i
        pfill(4*(i-1)+4) = p(0,i) ! E_i
      enddo
c Fill the decay momenta.
      do i = 1, imax_vmom
        pfill(4*imax_pmom + 4*(i-1)+1) = v(1,i) ! px_i
        pfill(4*imax_pmom + 4*(i-1)+2) = v(2,i) ! py_i
        pfill(4*imax_pmom + 4*(i-1)+3) = v(3,i) ! pz_i
        pfill(4*imax_pmom + 4*(i-1)+4) = v(0,i) ! E_i
      enddo
c Fill the event weight
      pfill(4*imax_pmom + 4*imax_vmom + 2) = dw
c Fill the event weight
      pfill(4*imax_pmom + 4*imax_vmom + 1) = nevent
c Write to file

      path = GetInputPath()
c      write(*,*) "pfill=",pfill
      if(NLO.eq.0) then
          filename="ntupleLO.dat"
          path = path(1:len_trim(path))//"/"//filename
          open(unit=40,file=path(1:len_trim(path)),position="append",status="old",action="write")
          l = 4*imax_pmom + 4*imax_vmom + 2
c          write(*,*) "l=",l 
          write(40,41) (pfill(i),i=1,l)
 41       format(*(G0,3X))        
          close(40)
      elseif (NLO.ne.0) then
          filename="ntupleNLO.dat"
          path = path(1:len_trim(path))//"/"//filename
          open(unit=40,file=path(1:len_trim(path)),position="append",status="old",action="write")
          
          write(40,42) (pfill(i),i=1,4*imax_pmom + 4*imax_vmom + 2)
 42       format(*(G0,3X))
          close(40) 
      else
      write(*,*) "There is problem."
      stop
      endif


      END

