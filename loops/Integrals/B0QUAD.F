c   Bubble correction: fin_QUADite part.           
c
c    musq = mu^2 e' una scala dimensionale esterna!!
c 
c  int d^dk/(2 pi)^d 1/(k^2)/(k+p1)^2 = 
c          N_ep * B0fin_QUAD(p1sq,musq)
c          N_ep = i/(4 pi)^2 (4 pi)^ep Gamma(1+ep) (musq)^(-ep)

c      complex*32 function B0fin_QUAD(p1sq,musq)
      subroutine B0fin_QUAD(p1sq,musq,B0_fin_QUAD) 
      implicit none
      complex*32 B0_fin_QUAD
      real*16  p1sq,musq
c     ris(1) = finite part, ris(2) = coeff of 1/ep, 
c     ris(3) = coeff of 1/ep^2      
      complex*32 ris(3)
      complex*32 l1
      complex*32 I
      parameter (I=(0,1q0))
      real*16 pi
      parameter (pi=3.14159265358979323846264338327950288q0)
      complex*32 ipi
      parameter (ipi=(0,3.14159265358979323846264338327950288q0))


c      complex*32 B0
c      logical debugB0C0D0
c      common/debug_B0C0D0/debugB0C0D0
      real*16  tiny
      parameter (tiny=1q-13)
      
c      if (debugB0C0D0) then
c         B0fin_QUAD = B0(p1sq,musq)
c         return
c      endif

      if (musq.lt.0) then
      write(*,*) 
     # 'ERROR in B0fin_QUAD: mu^2 MUST be a positive number'
      stop
      endif

c     if (p1sq.eq.0) then
      if (abs(p1sq).lt.tiny) then
c BUBBLE WITH ZERO EXTERNAL INVARIANTS
c     write(*,*) 'Warning: B0fin_QUAD called with external'//
c     #       ' invariant equal to zero'
         ris(3) = 0q0
         ris(2) = 0q0
         ris(1) = 0q0
         B0_fin_QUAD = ris(1)
         return
      endif
         
      l1 =  Log(abs(p1sq/musq))
      if (p1sq.gt.0q0) then
         l1 = l1 -ipi
      endif
      ris(3) = 0q0
      ris(2) = 1q0
      ris(1) = 2q0-l1
      B0_fin_QUAD = ris(1)
      return
      end


c --------------------------------------
c   B0finG2=B0finG(M1,0,s,musq)                                 
c --------------------------------------
      subroutine B0fin_QUAD_Div(s,musq,id,B0finG2Div)
      implicit none 
      real*16 s,musq
      integer id
      complex*32 B0finG2Div
      if(id.eq.1) then
         if(s.gt.1d-7) then
            B0finG2Div = 1.q0
         else
            B0finG2Div = 0q0
         endif
      else
      B0finG2Div = 0q0
      endif 
      end
