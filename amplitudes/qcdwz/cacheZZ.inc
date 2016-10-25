      complex*16 psi(2,-1:1,-5:5,16) !neg id for outgoing (bra), pos id for incoming (ket)
      complex*16 psiG(2,-1:1,-5:5,5,2,16) ! ,helq,qid,gid,kin
      complex*16 psiV(2,-1:1,3:4,-5:5,3,16)
      complex*16 psiVG(2,-1:1,3:4,-5:5,3,5,2)
      real*8 epsG(0:3,5,2,16)!,epsGG(0:3,4,5,2,2,16)
      real*8 momPsiG(0:4,-5:5,5,16),momPsiV(0:4,-5:5,3,16)
      real*8 momPsiVG(0:4,-5:5,3,5,16)
      common /QCDZZcacheC/ psi,psiG,psiV,psiVG
      common /QCDZZcacheR/ epsG,momPsiG,momPsiV,momPsiVG



