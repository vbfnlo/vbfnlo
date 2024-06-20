      complex*16 psi(2,-1:1,-5:5,16) !neg id for outgoing (bra), pos id for incoming (ket)
      complex*16 psiG(2,-1:1,-5:5,5,2,16) ! ,helq,qid,gid,kin
      complex*16 psiW(2,-5:5,16),psiWV(2,-5:5,16),psiV(2,-1:1,-5:5,3:4,16)
      complex*16 psiWG(2,-5:5,5,2,16)  !,psiVG(2,-1:1,-5:5,3:4,5,2,16)
      real*8 epsG(0:3,5,2,16),epsGG(0:3,4,5,2,2,16)
      real*8 momPsiG(0:4,-5:5,5,16),momPsiV(0:4,-5:5,3,16),momPsiWG(0:4,-5:5,5,16)
      common /QCDbraketsC/ psi,psiG,psiW,psiWV,psiV,psiWG
      common /QCDbraketsR/ epsG,epsGG,momPsiG,momPsiV,momPsiWG

