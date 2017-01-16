c --------------------------------------
c   A0finG2=A0finG(M,musq)                                 
c --------------------------------------
      complex*16 function A0fin(M,musq)
      implicit none
      real*8 pi
      complex*16 Ipi, Ieps 
      parameter(pi=3.14159265358979324d0,
     & Ipi=(0d0,3.14159265358979324d0),Ieps=(0d0,1d-38))

      double precision M,musq
      A0fin = M*M*(-dLog(M*M/musq)+1.d0)
      end

c --------------------------------------
c   A0finG2=A0finG(M,musq)                                 
c --------------------------------------
      complex*16 function A0finDiv(M,i)
      implicit none
      double precision M
      integer i
      if(i.eq.1) then
      A0finDiv = M*M
      else
      A0finDiv=0d0
      endif
      end
