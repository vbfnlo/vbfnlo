c ************************************************************************************
c Author: Matthias Kerner
c Date: 2012/13
c ************************************************************************************

c q1 q3 g2 g4 WZ Amplitude
c    ug2dgwzCS
c       Amplitude as needed for Catani Seymor
c       res(colcorr, pol1,pol2): collcorr=0..2 =  <l|l>, <l|T1T2|l>, <l|T1T3|l>
c          pol1,2=0,0: summed over all pol.
c          one pol_i=0: 2*Re(<1_i|TT|2_i>), second pol. fixed
c    ug2dgwzBorn
c       Amplitude as needed for Born process
c       res(0:2): amplitude for specified hel, index 1 and 2 for colstr for LHA output
      SUBROUTINE ug2dgwzCSNew(P,kin,id,fsign,nbos,lephel, usePreCalc,resDiag) ! sum over glu hel if ghel ==0

      use globalvars, only: ldoblha

      implicit none
#include "common.inc"
#include "brakets.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"
      real*8 P(0:3,max_p,max_kin), resDiag(0:2,0:2,0:2), res(0:2)
      integer kin, id(4),fsign(4),nbos,lephel, ikin
      logical lnlo,usePreCalc,lusePreCalc
      
      real*8 pbar(0:3,2),pdia(0:4,4),mpdia(0:4,4),momtemp(0:4)

      complex*16 epsW(0:3,7),epsWV(0:3,-1:1,7) 
      complex*16 epsV(0:3,-1:1,3:4,-1:1,7) ! ( , mu-hel, u/d , quark-hel), coupl. const to quarks included
      complex*16 epsA(0:3,-1:1),epsZ(0:3,-1:1)
      real*8 momV(0:4,7),momW(0:4,7),momWV(0:4,7)
      common /qcdwznew/ epsW, epsWV, epsV,epsA,epsZ,momV,momW,momWV

      complex*16 scres(0:2,2), amph(2,-1:1,-1:1), im
      parameter (im = (0d0,1d0))
      complex*16  amp(2,2,2),psitemp(2,-1:1)
      complex*16 s1c,s1r
      real*8 dotrr,mjj2
      
      integer i,j,h,mu, pol1,pol2,pol(2), tempi2(2), lepkin,qtype(2),idg1,idg2
      integer glu1,glu2
      lnlo = .true.
      ikin = kin
      lusePreCalc = usePreCalc
      goto 11      
      entry ug2dgwzBornNew(P,id,fsign,nbos,lephel, res)
      lusePreCalc = .false.
      ikin = 1
      lnlo= .false.
      
      res(0) = 0d0
      res(1) = 0d0
      res(2) = 0d0
 11   continue ! this part can also be calculated with precalcQCDVVjj

      if(id(2).ge.id(4)) then
        print*,__FILE__,__LINE__
        stop
      endif
      if(.not. lusePreCalc) then

      lepkin = ikin
      if(ikin .gt.7) lepkin = 1
      
      qtype(1) = nbos/10 ! (W+Z)? u type : d type
      qtype(2) = 7-qtype(1)

c use pdia for psi, -pdia for psibar
c -pdia for bosons !!!!!
      do mu=0,3
        pdia(mu,1) = p(mu,1,ikin) ! psi in initial state has sign 1
        pdia(mu,2) = p(mu,2,ikin) 
        pdia(mu,3) =-p(mu,3,ikin) ! psi in final state has sign -1
        pdia(mu,4) =-p(mu,4,ikin)
        pdia(4,mu+1) = 0d0
        pbar(mu,1) = p(mu,id(1),ikin)
        pbar(mu,2) = p(mu,id(3),ikin)
        do i=1,4
          mpdia(mu,i) = -pdia(mu,i)
        enddo
      enddo

c get bra and ket 
      if(id(1).le.2) then; tempi2(1) = 1; else; tempi2(1) = -1; endif
      if(id(3).ge.3) then; tempi2(2) = 1; else; tempi2(2) = -1; endif
      
      call psi0m(2,pbar,tempi2,Psi(1,-1,id(1),ikin))
      do i=1,2; do h=-1,1,2
        psi(i,h,-id(3),ikin) = psi(i,h,id(1)+1,ikin)
      enddo; enddo

      do i=2,4,2
        call polvec(p(0,id(i),ikin),1,epsG(0,id(i),1,ikin))
        call polvec(p(0,id(i),ikin),2,epsG(0,id(i),2,ikin))
      enddo
      do pol1=1,2; do pol2=1,2
      do mu = 0,3
c gluon copplings
c add a minus sign because of pdia convention 
        epsGG(mu,id(2),id(4),pol1,pol2,ikin)= -(2d0*dotrr(pdia(0,id(4)),epsG(0,id(2),pol1,ikin))*epsG(mu,id(4),pol2,ikin)
     &              -2d0*dotrr(pdia(0,id(2)),epsG(0,id(4),pol2,ikin))*epsG(mu,id(2),pol1,ikin)
     &        +dotrr(epsG(0,id(2),pol1,ikin),epsG(0,id(4),pol2,ikin))*(pdia(mu,id(2))-pdia(mu,id(4))))
     &        /mjj2(pdia(0,id(2)),pdia(0,id(4)))
      enddo
      enddo; enddo
c psiG
      do j=2,4,2
        call bra2r(psi(1,-1,-id(3),ikin),.true.,mpdia(0,id(3)),-1,mpdia(0,id(j)),epsG(0,id(j),1,ikin),psiG(1,-1,-id(3),id(j),1,ikin)
     &            ,momPsiG(0,-id(3),id(j),ikin))
        call bra2r(psi(1,-1,-id(3),ikin),.true.,mpdia(0,id(3)),-1,mpdia(0,id(j)),epsG(0,id(j),2,ikin),psiG(1,-1,-id(3),id(j),2,ikin)
     &            ,momPsiG(0,-id(3),id(j),ikin))

        call ket2r(psi(1,-1, id(1),ikin),.true.,pdia(0,id(1)),-1,mpdia(0,id(j)),epsG(0,id(j),1,ikin),psiG(1,-1, id(1),id(j),1,ikin)
     &            ,momPsiG(0,id(1),id(j),ikin))
        call ket2r(psi(1,-1, id(1),ikin),.true.,pdia(0,id(1)),-1,mpdia(0,id(j)),epsG(0,id(j),2,ikin),psiG(1,-1, id(1),id(j),2,ikin)
     &            ,momPsiG(0,id(1),id(j),ikin))
      enddo

c psiV, psiW
        call bra2c(psi(1,-1,-id(3),ikin),.true.,mpdia(0,id(3)),-1,momW(0,lepkin),epsW(0,lepkin),
     $             psiW(1,-id(3),ikin),momPsiV(0,-id(3),1,ikin))
        call ket2c(psi(1,-1, id(1),ikin),.true., pdia(0,id(1)),-1,momW(0,lepkin),epsW(0,lepkin),
     $             psiW(1, id(1),ikin),momPsiV(0,id(1),1,ikin))
      
        do j=3,4
          call bra2c(psi(1,-1,-id(3),ikin),.true.,mpdia(0,id(3)),-1,momV(0,lepkin),epsV(0,lephel,qtype(2),-1,lepkin),
     $               psiV(1,-1,-id(3),qtype(2),ikin),momPsiV(0,-id(3),2,ikin))
          call ket2c(psi(1,-1, id(1),ikin),.true.,pdia(0,id(1)),-1,momV(0,lepkin),epsV(0,lephel,qtype(1),-1,lepkin),
     $               psiV(1,-1, id(1),qtype(1),ikin),momPsiV(0,id(1),2,ikin))
        enddo

c psiVV
        call bra2c(psi(1,-1,-id(3),ikin),.true.,mpdia(0,id(3)),-1,momWV(0,lepkin),epsWV(0,lephel,lepkin),
     $             psiWV(1,-id(3),ikin),momPsiV(0,-id(3),3,ikin))
        call ket2c(psi(1,-1, id(1),ikin),.true.,pdia(0,id(1)),-1,momWV(0,lepkin),epsWV(0,lephel,lepkin),psiWV(1, id(1),ikin),
     $             momPsiV(0,id(1),3,ikin))

        call ket2c(psiV(1,-1,id(1),qtype(1),ikin),.false.,momPsiV(0,id(1),2,ikin),-1,momW(0,lepkin), epsW(0,lepkin),
     $             PsiTemp(1,-1),momTemp)
        call ket2c(psiW(1,id(1),ikin),.false.,momPsiV(0,id(1),1,ikin),-1,momV(0,lepkin), epsV(0,lephel,qtype(2),-1,lepkin),
     $             PsiTemp(1,1),momTemp)
        psiWV(1,id(1),ikin) = psitemp(1,-1) + psitemp(1,1) - psiWV(1,id(1),ikin)
        psiWV(2,id(1),ikin) = psitemp(2,-1) + psitemp(2,1) - psiWV(2,id(1),ikin)

        call bra2c(psiV(1,-1,-id(3),qtype(2),ikin),.false.,momPsiV(0,-id(3),2,ikin),-1,momW(0,lepkin), epsW(0,lepkin),
     $             PsiTemp(1,-1),momTemp)
        call bra2c(psiW(1,-id(3),ikin),.false.,momPsiV(0,-id(3),1,ikin),-1,momV(0,lepkin), epsV(0,lephel,qtype(1),-1,lepkin),
     $             PsiTemp(1,1),momTemp)
        psiWV(1,-id(3),ikin) = psitemp(1,-1) + psitemp(1,1) - psiWV(1,-id(3),ikin)
        psiWV(2,-id(3),ikin) = psitemp(2,-1) + psitemp(2,1) - psiWV(2,-id(3),ikin)

c psiWG
        do j=2,4,2
          do pol1 = 1,2
            call ket2c(psiG(1,-1,id(1),id(j),pol1,ikin),.false.,momPsiG(0,id(1),id(j),ikin),-1,momW(0,lepkin),epsW(0,lepkin),
     $                 PsiTemp(1,-1),momtemp)
            call ket2r(psiW(1,id(1),ikin),.false.,momPsiV(0,id(1),1,ikin),-1,mpdia(0,id(j)),epsG(0,id(j),pol1,ikin),
     $                 PsiTemp(1,1),momPsiWG(0,id(1),id(j),ikin))
            PsiWG(1,id(1),id(j),pol1,ikin) = PsiTemp(1,-1)+PsiTemp(1,1)
            PsiWG(2,id(1),id(j),pol1,ikin) = PsiTemp(2,-1)+PsiTemp(2,1)

            call bra2c(psiG(1,-1,-id(3),id(j),pol1,ikin),.false.,momPsiG(0,-id(3),id(j),ikin),-1,momW(0,lepkin),epsW(0,lepkin),
     $                 PsiTemp(1,-1),momtemp)
            call bra2r(psiW(1,-id(3),ikin),.false.,momPsiV(0,-id(3),1,ikin),-1,mpdia(0,id(j)),epsG(0,id(j),pol1,ikin),
     $                 PsiTemp(1,1),momPsiWG(0,-id(3),id(j),ikin))
            PsiWG(1,-id(3),id(j),pol1,ikin) = PsiTemp(1,-1)+PsiTemp(1,1)
            PsiWG(2,-id(3),id(j),pol1,ikin) = PsiTemp(2,-1)+PsiTemp(2,1)
          enddo
        enddo

      else ! use precalc
        lepkin = ikin
        if(ikin .gt.7) lepkin = 1
        
        qtype(1) = nbos/10 ! (W+Z)? u type : d type
        qtype(2) = 7-qtype(1)
      endif

      if(lnlo) then
        resDiag(1,0,0) = 0d0
        resDiag(2,0,0) = 0d0
        resDiag(0,0,0) = 0d0
      endif
      
      do pol1 = 1,2
      do pol2 = 1,2
        pol(1) = pol1
        pol(2) = pol2

c g(g1,g2) amplitudes
        amp(1,pol1,pol2) = -s1r(psi(1,-1,-id(3),ikin),epsGG(0,id(2),id(4),pol1,pol2,ikin),.true.,-1,PsiWV(1,id(1),ikin))
     &           -s1r(psiWV(1,-id(3),ikin),epsGG(0,id(2),id(4),pol1,pol2,ikin),.true.,-1,Psi(1,-1,id(1),ikin))
     &           -s1r(psiW(1,-id(3),ikin),epsGG(0,id(2),id(4),pol1,pol2,ikin),.true.,-1,PsiV(1,-1,id(1),qtype(1),ikin))
     &           -s1r(psiV(1,-1,-id(3),qtype(2),ikin),epsGG(0,id(2),id(4),pol1,pol2,ikin),.true.,-1,PsiW(1,id(1),ikin))
        amp(2,pol1,pol2) = -amp(1,pol1,pol2)

        do glu1=1,2
          glu2 = 3-glu1
          idg1=id(2*glu1)
          idg2=id(2*glu2)
c VVgg, ggVV, g(WV)g amplitudes
          amp(glu1,pol1,pol2) = amp(glu1,pol1,pol2)
     &     + s1r(PsiG(1,-1,-id(3),idg1,pol(glu1),ikin),epsG(0,idg2,pol(glu2),ikin),.false.,-1,PsiWV(1,id(1),ikin))
     &     - s1c(PsiG(1,-1,-id(3),idg1,pol(glu1),ikin),epsWV(0,lephel,lepkin),.true.,-1,PsiG(1,-1,id(1),idg2,pol(glu2),ikin))
     &     + s1r(PsiWV(1,-id(3),ikin),epsG(0,idg1,pol(glu1),ikin),.false.,-1,PsiG(1,-1,id(1),idg2,pol(glu2),ikin))
c all other amplitudes
          amp(glu1,pol1,pol2) = amp(glu1,pol1,pol2)
     &     + s1r(PsiWG(1,-id(3),idg1,pol(glu1),ikin),epsG(0,idg2,pol(glu2),ikin),.false.,-1,PsiV(1,-1,id(1),qtype(1),ikin))
     &     + s1c(PsiWG(1,-id(3),idg1,pol(glu1),ikin),epsV(0,lephel,qtype(1),-1,lepkin),.true.,-1,
     &                           PsiG(1,-1,id(1),idg2,pol(glu2),ikin))
     &     + s1r(PsiV(1,-1,-id(3),qtype(2),ikin),epsG(0,idg1,pol(glu1),ikin),.false.,-1,PsiWG(1,id(1),idg2,pol(glu2),ikin))
     &     + s1c(PsiG(1,-1,-id(3),idg1,pol(glu1),ikin),epsV(0,lephel,qtype(2),-1,lepkin),.true.,-1,
     &                           PsiWG(1,id(1),idg2,pol(glu2),ikin))
        enddo
        
c multiply coupl. const. and get |amp|^2
        amp(1,pol1,pol2) = amp(1,pol1,pol2)*GG(1)**2
        amp(2,pol1,pol2) = amp(2,pol1,pol2)*GG(1)**2
!       print*,amp(1,pol1,pol2),amp(2,pol1,pol2)
        if (lnlo) then
          if (ldoblha) then
            resDiag(0,pol1,pol2) = dble(
     &        blha_CF**2*blha_CA*                  (amp(1,pol1,pol2)*dconjg(amp(1,pol1,pol2))
     &                                             +amp(2,pol1,pol2)*dconjg(amp(2,pol1,pol2)))
     &       -blha_CF*blha_CA*(blha_CA-2d0*blha_CF)*amp(1,pol1,pol2)*dconjg(amp(2,pol1,pol2)))
            resDiag(1,pol1,pol2) = dble(
     &        blha_CF*blha_CA**2*(blha_CA-2d0*blha_CF)/4d0*amp(1,pol1,pol2)*dconjg(amp(1,pol1,pol2))
     &       -blha_CF**2*blha_CA**2/2d0*                   amp(2,pol1,pol2)*dconjg(amp(2,pol1,pol2))
     &       +blha_CF*blha_CA**2*(blha_CA-2d0*blha_CF)/2d0*amp(1,pol1,pol2)*dconjg(amp(2,pol1,pol2)))
            resDiag(2,pol1,pol2) = dble(
     &       -blha_CF*(blha_CA-2d0*blha_CF)/4d0*                     (amp(1,pol1,pol2)*dconjg(amp(1,pol1,pol2))
     &                                                               +amp(2,pol1,pol2)*dconjg(amp(2,pol1,pol2)))
     &       -blha_CF*blha_CA*(blha_CA-blha_CF)*(blha_CA-2d0*blha_CF)*amp(1,pol1,pol2)*dconjg(amp(2,pol1,pol2)))
          else
            resDiag(0,pol1,pol2) = dble(16d0*(amp(1,pol1,pol2)*dconjg(amp(1,pol1,pol2))+amp(2,pol1,pol2)*dconjg(amp(2,pol1,pol2)))
     &                                -4d0*amp(1,pol1,pol2)*dconjg(amp(2,pol1,pol2)))/3d0
            resDiag(1,pol1,pol2) = dble(amp(1,pol1,pol2)*dconjg(amp(1,pol1,pol2))-8d0*amp(2,pol1,pol2)*dconjg(amp(2,pol1,pol2))
     &                                +2d0*amp(1,pol1,pol2)*dconjg(amp(2,pol1,pol2)))
            resDiag(2,pol1,pol2) = dble(-amp(1,pol1,pol2)*dconjg(amp(1,pol1,pol2))-amp(2,pol1,pol2)*dconjg(amp(2,pol1,pol2))
     &                                -20d0*amp(1,pol1,pol2)*dconjg(amp(2,pol1,pol2)))/9d0
          endif
          resDiag(0,0,0)=resDiag(0,0,0)+resDiag(0,pol1,pol2)
          resDiag(1,0,0)=resDiag(1,0,0)+resDiag(1,pol1,pol2)
          resDiag(2,0,0)=resDiag(2,0,0)+resDiag(2,pol1,pol2)
        else
          res(0) = res(0)+ dble(16d0*(amp(1,pol1,pol2)*dconjg(amp(1,pol1,pol2))+amp(2,pol1,pol2)*dconjg(amp(2,pol1,pol2)))
     &                              -4d0*amp(1,pol1,pol2)*dconjg(amp(2,pol1,pol2)))/3d0
          res(1) = res(1)+dble(amp(1,pol1,pol2)*dconjg(amp(1,pol1,pol2)))
          res(2) = res(2)+dble(amp(2,pol1,pol2)*dconjg(amp(2,pol1,pol2)))
        endif
      enddo; enddo ! pol summation
      ! return <1|Tij|2>+<2|Tij|1>, pol1 fixed
      if (lnlo) then
        if (ldoblha) then ! <M,i-|Tij|M,i+>
c convert cartesian to helicity basis
          do glu1=1,2
            amph(glu1,+1,+1) = 
     &        ( +                     amp(glu1,1,1)
     &          -im*fsign(2)         *amp(glu1,2,1)
     &          -im*         fsign(4)*amp(glu1,1,2)
     &          -   fsign(2)*fsign(4)*amp(glu1,2,2) )/2d0
            amph(glu1,-1,+1) = 
     &        ( +                     amp(glu1,1,1)
     &          +im*fsign(2)*         amp(glu1,2,1)
     &          -im*         fsign(4)*amp(glu1,1,2)
     &          +   fsign(2)*fsign(4)*amp(glu1,2,2) )/2d0
            amph(glu1,+1,-1) = 
     &        ( +                     amp(glu1,1,1)
     &          -im*fsign(2)*         amp(glu1,2,1)
     &          +im*         fsign(4)*amp(glu1,1,2)
     &          +   fsign(2)*fsign(4)*amp(glu1,2,2) )/2d0
            amph(glu1,-1,-1) = 
     &        ( +                     amp(glu1,1,1)
     &          +im*fsign(2)*         amp(glu1,2,1)
     &          +im*         fsign(4)*amp(glu1,1,2)
     &          -   fsign(2)*fsign(4)*amp(glu1,2,2) )/2d0
          enddo

          scres(:,:) = 0d0
          do pol1 = -1,1,2
            scres(0,1) = scres(0,1) 
     &        + blha_CF**2*blha_CA*(amph(1,+1,pol1)*dconjg(amph(1,-1,pol1))
     &                             +amph(2,+1,pol1)*dconjg(amph(2,-1,pol1)))
     &        - blha_CF/2d0       *(amph(1,+1,pol1)*dconjg(amph(2,-1,pol1))
     &                             +amph(2,+1,pol1)*dconjg(amph(1,-1,pol1)))
            scres(1,1) = scres(1,1) 
     &        + blha_CF*blha_CA/4d0*       amph(1,+1,pol1)*dconjg(amph(1,-1,pol1))
     &        - blha_CF**2*blha_CA**2/2d0* amph(2,+1,pol1)*dconjg(amph(2,-1,pol1))
     &        + blha_CF*blha_CA/4d0*      (amph(1,+1,pol1)*dconjg(amph(2,-1,pol1))
     &                                    +amph(2,+1,pol1)*dconjg(amph(1,-1,pol1)))
            scres(2,1) = scres(2,1) 
     &        - blha_CF/blha_CA/4d0*          (amph(1,+1,pol1)*dconjg(amph(1,-1,pol1))
     &                                        +amph(2,+1,pol1)*dconjg(amph(2,-1,pol1)))
     &        - blha_CF*(blha_CA-blha_CF)/2d0*(amph(1,+1,pol1)*dconjg(amph(2,-1,pol1))
     &                                        +amph(2,+1,pol1)*dconjg(amph(1,-1,pol1)))
            scres(0,2) = scres(0,2) 
     &        + blha_CF**2*blha_CA*(amph(1,pol1,+1)*dconjg(amph(1,pol1,-1))
     &                             +amph(2,pol1,+1)*dconjg(amph(2,pol1,-1)))
     &        - blha_CF/2d0       *(amph(1,pol1,+1)*dconjg(amph(2,pol1,-1))
     &                             +amph(2,pol1,+1)*dconjg(amph(1,pol1,-1)))
            scres(1,2) = scres(1,2) 
     &        + blha_CF*blha_CA/4d0*       amph(1,pol1,+1)*dconjg(amph(1,pol1,-1))
     &        - blha_CF**2*blha_CA**2/2d0* amph(2,pol1,+1)*dconjg(amph(2,pol1,-1))
     &        + blha_CF*blha_CA/4d0*      (amph(1,pol1,+1)*dconjg(amph(2,pol1,-1))
     &                                    +amph(2,pol1,+1)*dconjg(amph(1,pol1,-1)))
            scres(2,2) = scres(2,2) 
     &        - blha_CF/blha_CA/4d0*          (amph(1,pol1,+1)*dconjg(amph(1,pol1,-1))
     &                                        +amph(2,pol1,+1)*dconjg(amph(2,pol1,-1)))
     &        - blha_CF*(blha_CA-blha_CF)/2d0*(amph(1,pol1,+1)*dconjg(amph(2,pol1,-1))
     &                                        +amph(2,pol1,+1)*dconjg(amph(1,pol1,-1)))
          enddo 
        else
          do pol1 = 1,2
            resDiag(0,pol1,0) = 2d0*dble(16d0*(amp(1,pol1,1)*dconjg(amp(1,pol1,2))+amp(2,pol1,1)*dconjg(amp(2,pol1,2)))
     &                                -2d0*amp(1,pol1,1)*dconjg(amp(2,pol1,2))-2d0*amp(2,pol1,1)*dconjg(amp(1,pol1,2)))/3d0
            resDiag(1,pol1,0) = 2d0*dble(amp(1,pol1,1)*dconjg(amp(1,pol1,2))-8d0*amp(2,pol1,1)*dconjg(amp(2,pol1,2))
     &                                +amp(1,pol1,1)*dconjg(amp(2,pol1,2))+amp(2,pol1,1)*dconjg(amp(1,pol1,2)))
            resDiag(2,pol1,0) = 2d0*dble(-amp(1,pol1,1)*dconjg(amp(1,pol1,2))-amp(2,pol1,1)*dconjg(amp(2,pol1,2))
     &                                -10d0*amp(1,pol1,1)*dconjg(amp(2,pol1,2))-10d0*amp(2,pol1,1)*dconjg(amp(1,pol1,2)))/9d0
            resDiag(0,0,pol1) = 2d0*dble(16d0*(amp(1,1,pol1)*dconjg(amp(1,2,pol1))+amp(2,1,pol1)*dconjg(amp(2,2,pol1)))
     &                                -2d0*amp(1,1,pol1)*dconjg(amp(2,2,pol1))-2d0*amp(2,1,pol1)*dconjg(amp(1,2,pol1)))/3d0
            resDiag(1,0,pol1) = 2d0*dble(amp(1,1,pol1)*dconjg(amp(1,2,pol1))-8d0*amp(2,1,pol1)*dconjg(amp(2,2,pol1))
     &                                +amp(1,1,pol1)*dconjg(amp(2,2,pol1))+amp(2,1,pol1)*dconjg(amp(1,2,pol1)))
            resDiag(2,0,pol1) = 2d0*dble(-amp(1,1,pol1)*dconjg(amp(1,2,pol1))-amp(2,1,pol1)*dconjg(amp(2,2,pol1))
     &                                -10d0*amp(1,1,pol1)*dconjg(amp(2,2,pol1))-10d0*amp(2,1,pol1)*dconjg(amp(1,2,pol1)))/9d0
          enddo
        endif
      else
        res(1) = res(1)/(res(1)+res(2))*res(0)
        res(2) = res(0)-res(1)
      endif
      if (ldoblha) then
        if (blha_amptype(blha_curproc).eq.0) then ! tree
          blha_amp(1)=blha_amp(1)+
     &      resDiag(0,0,0)
        else if (blha_amptype(blha_curproc).eq.2) then ! cctree
          ! 1-2
          call BLHA_cctree(1,2,
     &      + resDiag(1,0,0)
     &    )
          ! 1-3
          call BLHA_cctree(1,3,
     &      + resDiag(2,0,0)
     &    )
          ! 1-4 = -(CF<> + <12> + <13>)
          call BLHA_cctree(1,4,
     &      - ( blha_CF*resDiag(0,0,0)
     &          + resDiag(1,0,0)
     &          + resDiag(2,0,0) )
     &    )
          ! 2-3 = -(CF<> + <12> + <13>)
          call BLHA_cctree(2,3,
     &      - ( blha_CF*resDiag(0,0,0)
     &          + resDiag(1,0,0)
     &          + resDiag(2,0,0) )
     &    )
          ! 2-4 = (CF-CA)<>+<13>
          call BLHA_cctree(2,4,
     &      ( + (blha_CF-blha_CA)*resDiag(0,0,0)
     &        + resDiag(2,0,0) )
     &    )
          ! 3-4 = <12>
          call BLHA_cctree(3,4,
     &      + resDiag(1,0,0)
     &    )
        else if (blha_amptype(blha_curproc).eq.3) then ! sctree
          ! 2-1 = <2|12|2>
          call BLHA_sctree(2,1,
     &      scres(1,1)
     &    )
          ! 2-2 = CA<2|2>
          call BLHA_sctree(2,2,
     &      blha_CA*scres(0,1)
     &    )
          ! 2-3 = -(CF<2|2> + <2|12|2> + <2|13|2>)
          call BLHA_sctree(2,3,
     &      - ( blha_CF*scres(0,1)
     &          + scres(1,1)
     &          + scres(2,1) )
     &    )
          ! 2-4 = (CF-CA)<2|2> + <2|13|2>
          call BLHA_sctree(2,4,
     &      ( + (blha_CF-blha_CA)*scres(0,1)
     &        + scres(2,1) )
     &    )
          ! 4-1 = -(CF<4|4> + <4|12|4> + <4|13|4>)
          call BLHA_sctree(4,1,
     &      - ( blha_CF*scres(0,2)
     &          + scres(1,2)
     &          + scres(2,2) )
     &    )
          ! 4-2 = (CF-CA)<4|4> + <4|13|4>
          call BLHA_sctree(4,2,
     &      ( + (blha_CF-blha_CA)*scres(0,2)
     &        + scres(2,2) )
     &    )
          ! 4-3 = <4|12|4>
          call BLHA_sctree(4,3,
     &      scres(1,2)
     &    )
          ! 4-4 = CA<4|4>
          call BLHA_sctree(4,4,
     &      blha_CA*scres(0,2)
     &    )

        else 
          call BLHA_amptypeerror(blha_amptype(blha_curproc),__FILE__,__LINE__)
        endif
      endif
      end


      SUBROUTINE ug2dggwzOwnAllHelNew(P,id,fsign,nbos,lephel, res) ! sum over glu hel if ghel ==0

      use globalvars, only: ldoblha

      implicit none 
#include "common.inc"
#include "brakets.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"
      real*8 P(0:3,max_p,max_kin), res(0:6)
      integer id(5),fsign(5),nbos,lephel 

      complex*16 epsW(0:3,7),epsWV(0:3,-1:1,7) 
      complex*16 epsV(0:3,-1:1,3:4,-1:1,7) ! ( , mu-hel, u/d , quark-hel), coupl. const to quarks included
      complex*16 epsA(0:3,-1:1),epsZ(0:3,-1:1)
      real*8 momV(0:4,7),momW(0:4,7),momWV(0:4,7)
      common /qcdwznew/ epsW, epsWV, epsV,epsA,epsZ,momV,momW,momWV

      
      real*8 pbar(0:3,2),kbar(0:4,3),pdia(0:3,2),kdia(0:4,7),momTemp(0:4) ! kdia for glu: 1 2 3 12 13 23 123
! !       real*8 momPsiW(0:4,2),momPsiV(0:4,2),momPsiWV(0:4,2),momPsiG(0:4,2,7)
      real*8 momPsiVG(0:4,2,3)!,momPsiWG(0:4,2,3)
      real*8 triD(0:4,3,2,2,2),triV(0:4,3,2,2,2) ! 3-gluon polarications (see Nucl Phys B313 560, (5.13-14) indices: triD(mu,i) with j<k, triV(j)
      complex*16 PsiVVG(2,2,3,2),PsiVG(2,2,3,2)
      complex*16 PsiTemp(2,-1:1), amp(6)
      complex*16 s1r,temp
      real*8 eps3G(0:3,4:6,0:2,0:2,0:2),epsTemp(0:3)
      real*8 mass2, colsumsq,dotrr
      integer denom(6),idg(3)
      double precision CF(6,6)
      
      integer i,j,mu, lepkin,qtype(2),ggsign(3:5),pol1,pol2,pol3,pol(0:3)
      integer ig(3),perm, permutations(3,6)

c permutations convenient with Madgraph u g -> e+ ve d g mu+ mu- g  
      data (permutations(i,1),i=1,3) / 1, 2, 3/
      data (permutations(i,2),i=1,3) / 1, 3, 2/
      data (permutations(i,3),i=1,3) / 3, 1, 2/
      data (permutations(i,4),i=1,3) / 2, 1, 3/
      data (permutations(i,5),i=1,3) / 2, 3, 1/
      data (permutations(i,6),i=1,3) / 3, 2, 1/

c color factors from MG
      DATA Denom(1  )/            9/                                       
      DATA (CF(i,1  ),i=1  ,6  ) /    64,   -8,    1,   -8,    1,   10/    
C               T[ 5, 1, 2, 6, 9]                                          
      DATA Denom(2  )/            9/                                       
      DATA (CF(i,2  ),i=1  ,6  ) /    -8,   64,   -8,    1,   10,    1/    
C               T[ 5, 1, 2, 9, 6]                                          
      DATA Denom(3  )/            9/                                       
      DATA (CF(i,3  ),i=1  ,6  ) /     1,   -8,   64,   10,    1,   -8/    
C               T[ 5, 1, 9, 2, 6]                                          
      DATA Denom(4  )/            9/                                       
      DATA (CF(i,4  ),i=1  ,6  ) /    -8,    1,   10,   64,   -8,    1/    
C               T[ 5, 1, 6, 2, 9]                                          
      DATA Denom(5  )/            9/                                       
      DATA (CF(i,5  ),i=1  ,6  ) /     1,   10,    1,   -8,   64,   -8/    
C               T[ 5, 1, 6, 9, 2]                                          
      DATA Denom(6  )/            9/                                       
      DATA (CF(i,6  ),i=1  ,6  ) /    10,    1,   -8,    1,   -8,   64/    

      lepkin = 1
      do mu=0,3
        pbar(mu,1)=p(mu,id(1),1)
        pbar(mu,2)=p(mu,id(3),1)
        kbar(mu,1)=p(mu,id(2),1)
        kbar(mu,2)=p(mu,id(4),1)
        kbar(mu,3)=p(mu,id(5),1)
        
        pdia(mu,1)=pbar(mu,1)*fsign(1)
        pdia(mu,2)=pbar(mu,2)*fsign(3)
        kdia(mu,1)=kbar(mu,1)*fsign(2)
        kdia(mu,2)=kbar(mu,2)*fsign(4)
        kdia(mu,3)=kbar(mu,3)*fsign(5)

        kdia(mu,4)=kdia(mu,1)+kdia(mu,2)
        kdia(mu,5)=kdia(mu,1)+kdia(mu,3)
        kdia(mu,6)=kdia(mu,2)+kdia(mu,3)
        kdia(mu,7)=kdia(mu,1)+kdia(mu,6)
      enddo      
      kbar(4,1) = 0d0
      kbar(4,2) = 0d0
      kbar(4,3) = 0d0
      kdia(4,1) = 0d0
      kdia(4,2) = 0d0
      kdia(4,3) = 0d0

      kdia(4,4) = mass2(kdia(0,4))
      kdia(4,5) = mass2(kdia(0,5))
      kdia(4,6) = mass2(kdia(0,6))

      kdia(4,7) = mass2(kdia(0,7))
      
      qtype(1) = nbos/10 ! (W+Z)? u type : d type
      qtype(2) = 7-qtype(1)

      idg(1) = id(2)
      idg(2) = id(4)
      idg(3) = id(5)
      
      do pol1=1,2; do pol2=1,2
      do mu = 0,3
        eps3G(mu,4,pol1,pol2,0)= (2d0*dotrr(kdia(0,2),epsG(0,id(2),pol1,1))*epsG(mu,id(4),pol2,1)
     &              -2d0*dotrr(kdia(0,1),epsG(0,id(4),pol2,1))*epsG(mu,id(2),pol1,1)
     &              +dotrr(epsG(0,id(2),pol1,1),epsG(0,id(4),pol2,1))*(kdia(mu,1)-kdia(mu,2)))/kdia(4,4)
        eps3G(mu,5,pol1,0,pol2)= (2d0*dotrr(kdia(0,3),epsG(0,id(2),pol1,1))*epsG(mu,id(5),pol2,1)
     &              -2d0*dotrr(kdia(0,1),epsG(0,id(5),pol2,1))*epsG(mu,id(2),pol1,1)
     &              +dotrr(epsG(0,id(2),pol1,1),epsG(0,id(5),pol2,1))*(kdia(mu,1)-kdia(mu,3)))/kdia(4,5)
        eps3G(mu,6,0,pol1,pol2)= (2d0*dotrr(kdia(0,3),epsG(0,id(4),pol1,1))*epsG(mu,id(5),pol2,1)
     &              -2d0*dotrr(kdia(0,2),epsG(0,id(5),pol2,1))*epsG(mu,id(4),pol1,1)
     &              +dotrr(epsG(0,id(4),pol1,1),epsG(0,id(5),pol2,1))*(kdia(mu,2)-kdia(mu,3)))/kdia(4,6)
      enddo
      enddo;enddo

      do pol1=1,2; do pol2=1,2; do pol3=1,2
      do mu=0,3
        triD(mu,1,pol1,pol2,pol3) = (2d0*(dotrr(epsG(0,id(2),pol1,1),kdia(0,6))*eps3G(mu,6,0,pol2,pol3)
     &                    -dotrr(kdia(0,1),eps3G(0,6,0,pol2,pol3))*epsG(mu,id(2),pol1,1))
     &                 +   dotrr(epsG(0,id(2),pol1,1),eps3G(0,6,0,pol2,pol3))*(kdia(mu,1)-kdia(mu,6)))/kdia(4,7)
        triD(mu,2,pol1,pol2,pol3) = (2d0*(dotrr(epsG(0,id(4),pol2,1),kdia(0,5))*eps3G(mu,5,pol1,0,pol3)
     &                    -dotrr(kdia(0,2),eps3G(0,5,pol1,0,pol3))*epsG(mu,id(4),pol2,1))
     &                 +   dotrr(epsG(0,id(4),pol2,1),eps3G(0,5,pol1,0,pol3))*(kdia(mu,2)-kdia(mu,5)))/kdia(4,7)
        triD(mu,3,pol1,pol2,pol3) = (2d0*(dotrr(epsG(0,id(5),pol3,1),kdia(0,4))*eps3G(mu,4,pol1,pol2,0)
     &                    -dotrr(kdia(0,3),eps3G(0,4,pol1,pol2,0))*epsG(mu,id(5),pol3,1))
     &                 +   dotrr(epsG(0,id(5),pol3,1),eps3G(0,4,pol1,pol2,0))*(kdia(mu,3)-kdia(mu,4)))/kdia(4,7)
        triV(mu,1,pol1,pol2,pol3) = (2d0*epsG(mu,id(2),pol1,1)*dotrr(epsG(0,id(4),pol2,1),epsG(0,id(5),pol3,1))
     &                  -  epsG(mu,id(4),pol2,1)*dotrr(epsG(0,id(2),pol1,1),epsG(0,id(5),pol3,1))
     &                  -  epsG(mu,id(5),pol3,1)*dotrr(epsG(0,id(2),pol1,1),epsG(0,id(4),pol2,1)))/kdia(4,7)
        triV(mu,2,pol1,pol2,pol3) = (2d0*epsG(mu,id(4),pol2,1)*dotrr(epsG(0,id(2),pol1,1),epsG(0,id(5),pol3,1))
     &                  -  epsG(mu,id(2),pol1,1)*dotrr(epsG(0,id(4),pol2,1),epsG(0,id(5),pol3,1))
     &                  -  epsG(mu,id(5),pol3,1)*dotrr(epsG(0,id(2),pol1,1),epsG(0,id(4),pol2,1)))/kdia(4,7)
        triV(mu,3,pol1,pol2,pol3) = -triV(mu,1,pol1,pol2,pol3)-triV(mu,2,pol1,pol2,pol3)
      enddo
      enddo;enddo;enddo

c get <q,Vg| and |Vg,q>, 
      do i=1,3
      do pol1=1,2
        call ket2r(psiV(1,-1,id(1),qtype(1),1),.false.,momPsiV(0,id(1),2,1),-1,kdia(0,i), epsG(0,idg(i),pol1,1),
     $             PsiTemp(1,-1),momPsiVG(0,1,i))
        call ket2c(psiG(1,-1,id(1),idg(i),pol1,1),.false.,momPsiG(0,id(1),idg(i),1),-1,momV(0,lepkin), 
     $             epsV(0,lephel,qtype(1),-1,lepkin),PsiTemp(1,1),momTemp)
        PsiVG(1,1,i,pol1) = Psitemp(1,-1) + Psitemp(1,1)
        PsiVG(2,1,i,pol1) = Psitemp(2,-1) + Psitemp(2,1)
        call bra2c(psiG(1,-1,-id(3),idg(i),pol1,1),.false.,momPsiG(0,-id(3),idg(i),1),-1,momV(0,lepkin), 
     $             epsV(0,lephel,qtype(2),-1,lepkin),PsiTemp(1,-1),momPsiVG(0,2,i))
        call bra2r(psiV(1,-1,-id(3),qtype(2),1),.false.,momPsiV(0,-id(3),2,1),-1,kdia(0,i), epsG(0,idg(i),pol1,1),
     $             PsiTemp(1,1),momTemp)
        PsiVG(1,2,i,pol1) = Psitemp(1,-1) + Psitemp(1,1)
        PsiVG(2,2,i,pol1) = Psitemp(2,-1) + Psitemp(2,1)
c get <q,VVg| and |VVg,q>
        call ket2r(psiWV(1,id(1),1),.false.,momPsiV(0,id(1),3,1),-1,kdia(0,i), epsG(0,idg(i),pol1,1),PsiTemp(1,-1),momTemp)
        call ket2c(psiVG(1,1,i,pol1),.false.,momPsiVG(0,1,i),-1,momW(0,lepkin), epsW(0,lepkin),PsiTemp(1,0),momTemp)
        call ket2c(psiWG(1,id(1),idg(i),pol1,1),.false.,momPsiWG(0,id(1),idg(i),1),-1,momV(0,lepkin), 
     $             epsV(0,lephel,qtype(2),-1,lepkin),PsiTemp(1,1),momTemp)
        PsiVVG(1,1,i,pol1) = Psitemp(1,-1) + Psitemp(1,0) + Psitemp(1,1)
        PsiVVG(2,1,i,pol1) = Psitemp(2,-1) + Psitemp(2,0) + Psitemp(2,1)
        call ket2c(psiG(1,-1,id(1),idg(i),pol1,1),.false.,momPsiG(0,id(1),idg(i),1),-1,momWV(0,lepkin), epsWV(0,lephel,lepkin),
     $             PsiTemp(1,1),momTemp)
        PsiVVG(1,1,i,pol1) = PsiVVG(1,1,i,pol1) - Psitemp(1,1)
        PsiVVG(2,1,i,pol1) = PsiVVG(2,1,i,pol1) - Psitemp(2,1)
        call bra2r(psiWV(1,-id(3),1),.false.,momPsiV(0,-id(3),3,1),-1,kdia(0,i), epsG(0,idg(i),pol1,1),PsiTemp(1,-1),momTemp)
        call bra2c(psiVG(1,2,i,pol1),.false.,momPsiVG(0,2,i),-1,momW(0,lepkin), epsW(0,lepkin),PsiTemp(1,0),momTemp)
        call bra2c(psiWG(1,-id(3),idg(i),pol1,1),.false.,momPsiWG(0,-id(3),idg(i),1),-1,momV(0,lepkin), 
     $             epsV(0,lephel,qtype(1),-1,lepkin),PsiTemp(1,1),momTemp)
        PsiVVG(1,2,i,pol1) = Psitemp(1,-1) + Psitemp(1,1) + Psitemp(1,0)
        PsiVVG(2,2,i,pol1) = Psitemp(2,-1) + Psitemp(2,1) + Psitemp(2,0)
        call bra2c(psiG(1,-1,-id(3),idg(i),pol1,1),.false.,momPsiG(0,-id(3),idg(i),1),-1,momWV(0,lepkin), epsWV(0,lephel,lepkin),
     $             PsiTemp(1,1),momTemp)
        PsiVVG(1,2,i,pol1) = PsiVVG(1,2,i,pol1) - Psitemp(1,1)
        PsiVVG(2,2,i,pol1) = PsiVVG(2,2,i,pol1) - Psitemp(2,1)
      enddo
      enddo

      do i=0,6
        res(i) = 0d0
      enddo
      colsumsq = 0d0
c get amps ans ampsq for all polarisations
      do pol1=1,2; do pol2=1,2; do pol3=1,2

      do perm = 1,6
        ig(1) = permutations(1,perm)
        ig(2) = permutations(2,perm)
        ig(3) = permutations(3,perm)
        ! 12 13 and 23 ordered
        if(perm.le.3) then; ggsign(3) = 1; else; ggsign(3) = -1; endif
        if(perm.le.2.or.perm.eq.4) then; ggsign(4) = 1; else; ggsign(4) = -1; endif
        if(perm/2.eq.2.or.perm.eq.1) then; ggsign(5) = 1; else; ggsign(5) = -1; endif

c 1 gluon splitting amplitudes
        i = ig(2)+ig(3)+1
        pol(1) = pol1; pol(2) = pol2; pol(3) = pol3
        pol(0) = pol(ig(1)); pol(ig(1)) = 0 ! swap to get convention of epsG and psi*G
        temp = s1r(psiVVG(1,2,ig(1),pol(0)),eps3G(0,i,pol(1),pol(2),pol(3)),.true.,-1,Psi(1,-1,id(1),1))
     &        +s1r(psiVG(1,2,ig(1),pol(0)),eps3G(0,i,pol(1),pol(2),pol(3)),.true.,-1,PsiW(1,id(1),1))
     &        +s1r(psiWG(1,-id(3),idg(ig(1)),pol(0),1),eps3G(0,i,pol(1),pol(2),pol(3)),.true.,-1,PsiV(1,-1,id(1),qtype(1),1))
     &        +s1r(psiG(1,-1,-id(3),idg(ig(1)),pol(0),1),eps3G(0,i,pol(1),pol(2),pol(3)),.true.,-1,PsiWV(1,id(1),1))
        amp(perm) = -temp*ggsign(i-1)

        i = ig(1)+ig(2)+1
        pol(ig(1)) = pol(0); pol(0) = pol(ig(3)); pol(ig(3)) = 0 ! swap back and swap for next calls
        temp = s1r(psiWV(1,-id(3),1),eps3G(0,i,pol(1),pol(2),pol(3)),.true.,-1,PsiG(1,-1,id(1),idg(ig(3)),pol(0),1))
     &        +s1r(psiV(1,-1,-id(3),qtype(2),1),eps3G(0,i,pol(1),pol(2),pol(3)),.true.,-1,PsiWG(1,id(1),idg(ig(3)),pol(0),1))
     &        +s1r(psiW(1,-id(3),1),eps3G(0,i,pol(1),pol(2),pol(3)),.true.,-1,PsiVG(1,1,ig(3),pol(0)))
     &        +s1r(psi(1,-1,-id(3),1),eps3G(0,i,pol(1),pol(2),pol(3)),.true.,-1,PsiVVG(1,1,ig(3),pol(0)))
        amp(perm) = amp(perm)-temp*ggsign(i-1)
        pol(ig(3)) = pol(0) ! swap back
        
c g(g1,g2,g3) amplitudes        
        do mu=0,3
          epsTemp(mu) = ggsign(ig(2)+ig(3))*triD(mu,ig(1),pol1,pol2,pol3)-ggsign(ig(1)+ig(2))*triD(mu,ig(3),pol1,pol2,pol3)+
     $                         triV(mu,ig(2),pol1,pol2,pol3)
        enddo
        amp(perm)= amp(perm)
     &           +s1r(psi(1,-1,-id(3),1),epsTemp,.true.,-1,PsiWV(1,id(1),1))
     &           +s1r(psiWV(1,-id(3),1),epsTemp,.true.,-1,Psi(1,-1,id(1),1))
     &           +s1r(psiW(1,-id(3),1),epsTemp,.true.,-1,PsiV(1,-1,id(1),qtype(1),1))
     &           +s1r(psiV(1,-1,-id(3),qtype(2),1),epsTemp,.true.,-1,PsiW(1,id(1),1))

! c amplitudes without glu splitting
        amp(perm) = amp(perm)
     &             +s1r(psiVVG(1,2,ig(1),pol(ig(1))),epsG(0,idg(ig(2)),pol(ig(2)),1),.false.,-1,
     &                  PsiG(1,-1,id(1),idg(ig(3)),pol(ig(3)),1))
     &             +s1r(psiVG(1,2,ig(1),pol(ig(1))),epsG(0,idg(ig(2)),pol(ig(2)),1),.false.,-1,
     &                  PsiWG(1,id(1),idg(ig(3)),pol(ig(3)),1))
     &             +s1r(psiWG(1,-id(3),idg(ig(1)),pol(ig(1)),1),epsG(0,idg(ig(2)),pol(ig(2)),1),.false.,-1,
     &                  PsiVG(1,1,ig(3),pol(ig(3))))
     &             +s1r(psiG(1,-1,-id(3),idg(ig(1)),pol(ig(1)),1),epsG(0,idg(ig(2)),pol(ig(2)),1),.false.,-1,
     &                  PsiVVG(1,1,ig(3),pol(ig(3))))
        amp(perm) = amp(perm) * G**3
      enddo ! perm

c BLHA colour factors
      if (ldoblha) then
        denom(:) = 1d0
        do i=1,6
          CF(i,i            ) =  blha_CA*blha_CF**3
          CF(i,mod(i+ 1,6)+1) =  blha_CF*(blha_CA-2d0*blha_CF)/4d0
          CF(i,mod(i+ 3,6)+1) =  blha_CF*(blha_CA-2d0*blha_CF)/4d0
          CF(i,mod( 8-i,6)+1) = -blha_CA*blha_CF**2*(blha_CA-2d0*blha_CF)/2d0
          CF(i,mod(10-i,6)+1) = -blha_CA*blha_CF**2*(blha_CA-2d0*blha_CF)/2d0
          CF(i,7-i          ) =  blha_CA*blha_CF*(blha_CA-2d0*blha_CF)*(blha_CA-blha_CF)/2d0
        enddo
      endif

c square amp
      DO I = 1, 6
          TEMP = (0.D0,0.D0)
          DO J = 1, 6
              TEMP = TEMP + CF(J,I)*AMP(J)
          ENDDO
          res(0) =res(0)+TEMP*DCONJG(AMP(I))/DENOM(I)
          res(i) =res(i)+ AMP(I)*DCONJG(AMP(I))/DENOM(I)
      ENDDO
      enddo; enddo;enddo ! pol
      
      colsumsq = 0d0
      do i=1,6
        colsumsq = colsumsq+res(i)
      enddo
      DO i=1,6
         res(i)=res(i)/colsumsq*res(0)
      ENDDO

      if (ldoblha) then
        if (blha_amptype(blha_curproc).eq.0) then ! tree
          blha_amp(1)=blha_amp(1)+
     &      res(0)
        else 
          call BLHA_amptypeerror(blha_amptype(blha_curproc),__FILE__,__LINE__)
        endif
      endif

      end


