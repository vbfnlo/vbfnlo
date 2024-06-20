ccccccccccccccccccccccccccccccccccccccccccccccccc
c Parameters used for the frame transformation
ccccccccccccccccccccccccccccccccccccccccccccccccc

c  transformation matrices
      real*8 TT, TTinv
      common /transm/ TT(0:3,0:3,5), TTinv(0:3,0:3,5)

      real*8 ecm_real
      common /transm/ ecm_real
 
      real*8 sc
      common /cms/ sc
