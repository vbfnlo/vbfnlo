      SUBROUTINE amp2q2gZZVirt(P,id,fsign,nlo,bos,qhel0,ghel10,ghel20,lephel0,rn, res) ! sum over glu hel if ghel ==0
      implicit none
#include "common.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/scales.inc"
      real*8 P(0:3,max_p,max_kin), res(0:2,3:4),rn
      integer id(4),fsign(4),lephel,qhel,ghel1,ghel2,nlo, qhel0,ghel10,ghel20,lephel0,bos
      
      real*8 pbar(0:3,2),qdia(0:3,2),gdia(0:4,3),epsGr(0:3,3)
      complex*16 epsG(0:3,3)

      complex*16  amp(2),psitemp(2,-1:1,2),ampv(3,3:4),ampb(2,3:4),ampborn(2,3:4),ampi(3)
      complex*16 ortb(3),ortv(3),orti(3),ort(3)
      complex*16 s1c, dotcc, dotrc,s1r
      real*8 mass2,dotrr,mjj2
      
      complex*16 psiVGlocal(2,-1:1,3:4,2,2,2)
      integer i,j,mu, tempi2(2)
      integer glu1,glu2,iV1,iV2, q
      

      real*8 momPsiG(0:4,2,2),mompsiV(0:4,2,3),momtemp(0:4),td
      complex*16 psi(2,-1:1,2) !neg id for outgoing (bra), pos id for incoming (ket)
      complex*16 psiG(2,-1:1,2,2) ! ,helq,qid,gid,kin
      complex*16 psiV(2,-1:1,3:4,2,3),momtempc(0:4)
      complex*16 gluV(0:3,-1:1,3:4,0:3),checkg(5,3),checkgb(2)
      
      logical ltop, lborn,lbox,lpent, lhex, lfloop,lcsi,lsimple, ldebug,isNaN,lgaugeBorn,lZZ,lZ4l,lhiggs

      external isNaN
      integer div, idia, ngauge, nGaugeAbe, nGaugeNoAbe, comp,flsplit
      real*8 musq, kdia(0:4,3)
      real*8 CF,CA,CFCA,TR
      parameter( CF = 4d0/3d0, CA = 3d0, CFCA=-1d0/6d0, TR=0.5d0)
      
      real*8 gaugereal(2,2),logmuosij(3,4),floopgauge(10)
      complex*16 tempres(13),tempresb,hexgauge(4,30),hexgaugeb(30),tempc,tempborn
      complex*16 pengauge(3,4),pengaugeb(3,4)
      complex*16 tempvirt(0:3,3:4),ctemp(18),tempcurr(0:3,-1:1,2),nothing(9)


      logical wardidtest(0:22),wardidtest_q(0:22),GaugeOK(2),gaugetestfailed
      logical wardpAbe(2),wardpNoAbe(2),wardhAbe(2)
      logical wardhNoAbe(2)
      logical wardfb(2)
      common/wardpfermionbox/wardfb
      common/wardpenAbe/wardpAbe
      common/wardpenNoAbe/wardpNoAbe
      common/wardhexAbe/wardhAbe
      common/wardHexNoAbe/wardhNoAbe
      logical GaugeE,GaugeO
      common/FloopPenggzz/GaugeE,GaugeO     
c V currents
c V1/V2/V12,
c    for qhel=0:
c       AZ=1,2 : A and Z current wo coupling
c    qhel = -1:1:
c       AZ=3,4 : V coupling to u/d quark with qhel inclouding couplings
      complex*16 epsV,w4l(0:4) ! (mu, qhel, 12=AZ 34=including coupl to ud, V1 V2 V12, hel->0=rnd hel, kin)
      real*8 momV
      common /QCDZZcurrents/ epsV(0:3,-1:1,4,3,0:4,7),w4l,momV(0:4,3,7)
      logical lAZ(0:2,3) ! V.or.A A V,    V1 V2 V12
      common /QCDZZcacheL/ lAZ

      GaugeE=.false.
      GaugeO=.false.
      
      qhel = qhel0
      ghel1= ghel10
      ghel2= ghel20
      lephel = lephel0

      ltop = NFl.ge.5
! ! !       ltop=.true.
      ldebug=.false.
      lgaugeBorn=.true. 
      
      div=0
      
      if(ldebug) then
        print*,'new amp'
        print*,id
        print*,qhel,lephel,ghel1,ghel2
      endif
!       pause
      musq=mursq(1,1)

!       if(nloContrib.eq.0) then
      lborn  = nlo.eq.1
      lfloop = nlo.ge.4
      lsimple= nlo.eq.2 ! fermion loops, glu triangle, self energies, counterterms
      lbox   = nlo.eq.2
      lpent  = nlo.eq.2
      lhex   = nlo.eq.2
      lCSI   = nlo.eq.2

      flsplit=0  ! inlcude all floops

! !       flsplit = nlo ! pentagon = 5, other floop=4
      
!        lborn  = .false.
!        lfloop = .true.
!        lsimple= .true.
!        lbox   = .true.
!        lpent  = .true.
!        lhex   = .true.
!        lCSI   = .true.

      lZZ =  lAZ(0,1)
      lZ4l = lAZ(0,3)
      lhiggs=nfl.ge.5!.true.
!       
!       lborn  = .false.
!       lfloop = .true.
!       lsimple= .false.
!       lbox   = .false.
!       lpent  = .false.
!       lhex   = .false.
!       lCSI   = .false.


      lfloop=lfloop.and.div.eq.0
      res=0d0
      
      ampb=0d0
      ampv=0d0
      ampi=0d0
      ampborn=0d0
      
!       epsV(:,:,:,1:2,:,1) = 0d0 
!       w4l=0d0
      
!       epsV(:,:,:,3,:,1) = 0d0 
!       epsV(:,:,4,3,:,1) = 0d0 
      epsGr=0d0

      if(id(2).ge.id(4)) then
        print*,__FILE__,__LINE__
        stop
      endif

c use pdia for psi, -pdia for psibar
c -pdia for bosons !!!!!
      qdia(:3,1) = fsign(1)*p(:,id(1),1) ! psi in initial state has sign 1
      qdia(:3,2) = fsign(3)*p(:,id(3),1) ! psi in final state has sign -1
      gdia(:3,1) = fsign(2)*p(:,id(2),1)
      gdia(:3,2) = fsign(4)*p(:,id(4),1)
      gdia(:3,3) = gdia(:3,1) + gdia(:3,2)
      gdia(4,:2) = 0d0
      gdia(4,3) = mass2(gdia(0,3))

c get bra and ket 
      tempi2(1) = merge(1,-1,id(1).le.2)
      tempi2(2) = merge(1,-1,id(3).ge.3)
      pbar(:,1) = p(:,id(1),1)
      pbar(:,2) = p(:,id(3),1)
      call psi0m(2,pbar,tempi2,psi(1,-1,1))

      do i=2,4,2
        call polvec(p(0,id(2),1),ghel1,epsGr(0,1))
        call polvec(p(0,id(4),1),ghel2,epsGr(0,2))
      enddo
!       epsGr(:3,2) = gdia(:3,2)
      epsG= epsGr

      
!       if(lborn) then

c psiG
      do j=1,2
        call bra2r(psi(1,qhel,2),.true.,qdia(0,2),qhel,gdia(0,j),epsGr(0,j),psiG(1,qhel,2,j),momPsiG(0,2,j))
        call ket2r(psi(1,qhel,1),.true.,qdia(0,1),qhel,gdia(0,j),epsGr(0,j),psiG(1,qhel,1,j),momPsiG(0,1,j))
      enddo
      

      psiV=0d0

      do iV1=3,1,-1 ! V1, V2, V12
      do q=3,4 ! u/d
      if(lAZ(0,iV1)) then
        call ket2c(psi(1,qhel,1),.true.,qdia(0,1),qhel,momV(0,iV1,1),epsV(0,qhel,q,iV1,lephel,1),
     $             psiV(1,qhel,q,1,iV1),mompsiV(0,1,iV1))
        call bra2c(psi(1,qhel,2),.true.,qdia(0,2),qhel,momV(0,iV1,1),epsV(0,qhel,q,iV1,lephel,1),
     $             psiV(1,qhel,q,2,iV1),mompsiV(0,2,iV1))

c contract with 2nd V
        if(iV1.le.2) then
          iV2= 3-iV1
          if(lAZ(0,iV2)) then
            call ket2c(psiV(1,qhel,q,1,iV1),.false.,mompsiV(0,1,iV1),qhel,momV(0,iV2,1),
     $                 epsV(0,qhel,q,iV2,lephel,1),psitemp(1,qhel,1),momtemp)
            psiV(1:2,qhel,q,1,3)= psiV(1:2,qhel,q,1,3)+psitemp(:,qhel,1)
            call bra2c(psiV(1,qhel,q,2,iV1),.false.,mompsiV(0,2,iV1),qhel,momV(0,iV2,1),
     $                 epsV(0,qhel,q,iV2,lephel,1),psitemp(1,qhel,1),momtemp)
            psiV(1:2,qhel,q,2,3)= psiV(1:2,qhel,q,2,3)+psitemp(1:2,qhel,1)
          endif
        endif
      endif
      enddo;enddo

c 3glu vert
      do mu = 0,3
c gluon copplings
c momenta are flowing out -> extra sign
        epsGr(mu,3)= -(2d0*dotrr(gdia(0,2),epsGr(0,1))*epsGr(mu,2)
     &              -2d0*dotrr(gdia(0,1),epsGr(0,2))*epsGr(mu,1)
     &        +dotrr(epsGr(0,1),epsGr(0,2))*(gdia(mu,1)-gdia(mu,2)))
     &        /mjj2(gdia(0,1),gdia(0,2)) ! * [ta.tb]
      enddo
      epsG(:,3) = epsGr(:,3)


      do iV1=2,1,-1 ! V1, V2
      do q=3,4 ! u/d
      if(lAZ(0,iV1)) then
c contract with g
          do j=1,2
            call bra2r(psiV(1,qhel,q,2,iV1),.false.,mompsiV(0,2,iV1),qhel,gdia(0,j),
     $                 epsGr(0,j),psitemp(1,qhel,1),momtemp)
            call bra2c(psiG(1,qhel,2,j),.false.,momPsiG(0,2,j),qhel,momV(0,iV1,1),
     $                 epsV(0,qhel,q,iV1,lephel,1),psitemp(1,qhel,2),momtemp)
            psiVGlocal(1:2,qhel,q,2,iV1,j)=psitemp(1:2,qhel,1)+psitemp(1:2,qhel,2)

            call ket2r(psiV(1,qhel,q,1,iV1),.false.,mompsiV(0,1,iV1),qhel,gdia(0,j),
     $                 epsGr(0,j),psitemp(1,qhel,1),momtemp)
            call ket2c(psiG(1,qhel,1,j),.false.,momPsiG(0,1,j),qhel,momV(0,iV1,1),
     $                 epsV(0,qhel,q,iV1,lephel,1),psitemp(1,qhel,2),momtemp)
            psiVGlocal(1:2,qhel,q,1,iV1,j)=psitemp(1:2,qhel,1)+psitemp(1:2,qhel,2)
          enddo
      endif
      enddo;enddo

c build currents and calc amplitudes
      call curr(qhel,Psi(1,-1,2),Psi(1,-1,1),gluV(0,-1,3,0))
      gluV(:,-1:1:2,4,0) = gluV(:,-1:1:2,3,0)
      do q=3,4
        do iV1=1,3
          if(lAZ(0,iV1)) then
            call curr(qhel,Psi(1,-1,2),PsiV(1,-1,q,1,iV1),tempcurr(0,-1,1))
            call curr(qhel,PsiV(1,-1,q,2,iV1),Psi(1,-1,1),tempcurr(0,-1,2))
            gluV(:,-1:1:2,q,iV1) = tempcurr(:,-1:1:2,1) + tempcurr(:,-1:1:2,2)
          endif
        enddo
        if(lAZ(0,1)) then
          call curr(qhel,PsiV(1,-1,q,2,2),PsiV(1,-1,q,1,1),tempcurr(0,-1,1))
          call curr(qhel,PsiV(1,-1,q,2,1),PsiV(1,-1,q,1,2),tempcurr(0,-1,2))
          gluV(:,-1:1:2,q,3) = gluV(:,-1:1:2,q,3)+tempcurr(:,-1:1:2,1) + tempcurr(:,-1:1:2,2)
        endif


        amp(:)=0d0
c amp(1) = tb ta with Peskin-S. conventions of factors i

c g(g1,g2) amplitudes
!         amp(2) = +s1r(psi(1,qhel,2),epsGr(0,3),.true.,qhel,PsiV(1,qhel,q,1,3))
!      &           +s1r(psiV(1,qhel,q,2,3),epsGr(0,3),.true.,qhel,Psi(1,qhel,1))
!         if(lAZ(0,1)) then
!           amp(2) = amp(2)
!      &           +s1r(psiV(1,qhel,q,2,1),epsGr(0,3),.true.,qhel,PsiV(1,qhel,q,1,2))
!      &           +s1r(psiV(1,qhel,q,2,2),epsGr(0,3),.true.,qhel,PsiV(1,qhel,q,1,1))
!         endif
        ! Vertex * (iM) = (-I)(-I)
        amp(2) = -dotrc(epsGr(0,3),gluV(0,qhel,q,3))
        amp(1) = -amp(2)

        do glu1=1,2
          glu2 = 3-glu1
c VVgg, ggVV, g(WV)g amplitudes
          amp(glu1) = amp(glu1)
     &     - s1r(PsiG(1,qhel,2,glu2),epsGr(0,glu1),lgaugeBorn,qhel,PsiV(1,qhel,q,1,3))
     &     - s1c(PsiG(1,qhel,2,glu2),epsV(0,qhel,q,3,lephel,1),.true.,qhel,PsiG(1,qhel,1,glu1))
     &     - s1r(PsiV(1,qhel,q,2,3),epsGr(0,glu2),lgaugeBorn,qhel,PsiG(1,qhel,1,glu1))

! c all other amplitudes
        if(lAZ(0,1)) then
        do iV1=1,2
        iV2=3-iV1
          amp(glu1) = amp(glu1)
     &     - s1r(PsiVGlocal(1,qhel,q,2,iV1,glu2),epsGr(0,glu1),lgaugeBorn,qhel,PsiV(1,qhel,q,1,iV2))
     &     - s1c(PsiVGlocal(1,qhel,q,2,iV1,glu2),epsV(0,qhel,q,iV2,lephel,1),.true.,qhel,PsiG(1,qhel,1,glu1))
        enddo !iV1
        endif
        enddo !glu1

        ampborn(1,q) = amp(1)
        ampborn(2,q) = amp(2)
        
c multiply coupl. const. and get |amp|^2
        amp(1) = amp(1)*GG(1)**2
        amp(2) = amp(2)*GG(1)**2

        res(0,q) = res(0,q)+ dble(16d0*(amp(1)*dconjg(amp(1))+amp(2)*dconjg(amp(2)))
     &                              -4d0*amp(1)*dconjg(amp(2)))/3d0
        res(1,q) = res(1,q)+dble(amp(1)*dconjg(amp(1)))
        res(2,q) = res(2,q)+dble(amp(2)*dconjg(amp(2)))

        res(1,q) = res(1,q)/(res(1,q)+res(2,q))*res(0,q)
        res(2,q) = res(0,q)-res(1,q)
      enddo
      
      
!       endif ! born
      
c---------------------------------
c---------------------------------
c---------------------------------
c       V I R T U A L S 
c---------------------------------
c---------------------------------
c---------------------------------
      gaugetestfailed = .false.

c change momentum convention -> all momenta flowing in
      gdia(:3,:) = -gdia(:3,:)
      qdia(:,2) = -qdia(:,2)
      kdia(:3,:) = -momV(:3,:,1)
      kdia(4,:)=momv(4,:,1)

      
c simple diagrams
      if(lsimple) then
      
      GaugeOK(1) = .true.

c glu self enregie
      do q=3,4
        tempborn=-dotcc(epsG(0,3),gluV(0,qhel,q,3)) ! 1/q^2 included in epsG ! Vertex*(iM) = (-I) (-I)
        call gluonSelf(gdia(4,3),musq,nfl,tempc,div)
        tempres(1) = tempc*tempborn

        if(ltop) then
          call gluonSelfMass(gdia(4,3),musq,tmass,tempc,div)
          tempres(1) = tempres(1) + tempc*tempborn
        endif

        ampv(1,q) = ampv(1,q)-tempres(1)
        ampv(2,q) = ampv(2,q)+tempres(1)

c V- and U-loop contributions to 3 gluon vertex
      do mu=0,3
        momtemp(mu) = -gdia(mu,3)
      enddo
!       print*,'test',tempborn
!        call Vert3G(gdia(0,1),gdia(0,2),momtemp(0),epsG(0,1),epsG(0,2),gluV(0,qhel,q,3),musq,tempborn,-1)
!        tempborn = tempborn/gdia(4,3) ! gluV*prop: (-I)(-I), Vert3G uses Peskin convention -> additional -1
!        print*,'test',tempborn
!        pause
      
      call Vert3G(gdia(0,1),gdia(0,2),momtemp(0),epsG(0,1),epsG(0,2),gluV(0,qhel,q,3),musq,tempc,div)
      tempc=+tempc/gdia(4,3)! gluV*prop: (-I)(-I), Vert3G uses Peskin convention -> additional -1

      ampv(1,q) = ampv(1,q)-tempc
      ampv(2,q) = ampv(2,q)+tempc
      enddo ! q

c glu glu glu* ! 2 directions have opposite sign but different color structure
      comp=1
      do q=3,4
        call FermionVertCouplingT(1d0,1d0,1d0,1d0,1d0,1d0,
     $          0d0,
     $          gdia(0,1),gdia(0,2),momtemp(0),
     $          epsG(0,1),epsG(0,2),gluV(0,qhel,q,3),
     $          musq,comp,comp,tempres(1),floopgauge(1),div)
      tempres(1) = -tempres(1)/2d0/gdia(4,3) ! gluV*prop: (-I)(-I)
      comp=-1
c color structure is Tr[ta,tb,tc]tc=tatb/2-delta_ab/12
c routine with glu 1 and 2 interchanged will have different sign and color structure tbta/2-delta_ab/8
c the delta_ab color structure will cancel
      ampv(2,q) = ampv(2,q)+tempres(1)*NFl ! 5 light quarks
      ampv(1,q) = ampv(1,q)-tempres(1)*NFl

      if(comp.eq.1.and.floopgauge(1).gt.GaugeAccuracy) then
        gaugeOK(1) = .false.
      endif
      enddo ! q

      if(ltop) then
      comp=1
      do q=3,4
        call FermionVertCouplingT(1d0,1d0,1d0,1d0,1d0,1d0,
     $          tmass,
     $          gdia(0,1),gdia(0,2),momtemp(0),
     $          epsG(0,1),epsG(0,2),gluV(0,qhel,q,3),
     $          musq,comp,comp,tempres(1),floopgauge(1),div)
      comp=-1
      tempres(1) = -tempres(1)/2d0/gdia(4,3)
      ampv(2,q) = ampv(2,q)+tempres(1)
      ampv(1,q) = ampv(1,q)-tempres(1)

        if(comp.eq.1.and.floopgauge(1).gt.GaugeAccuracy) then
          gaugeOK(1) = .false.
        endif
      enddo ! q

      if(.not.lbox) then
        boxcount2m = boxcount2m + 1
        boxcount2sub(sub_number) = boxcount2sub(sub_number) + 1
      endif
      if(.not.gaugeOK(1)) then
        boxcountm   = boxcountm   + 1
        boxcountm_Q = boxcountm_Q + 1
        boxcountsub  (sub_number) = boxcountsub  (sub_number) + 1
        boxcountsub_q(sub_number) = boxcountsub_q(sub_number) + 1
        gaugetestfailed=.true.
      endif

      endif


      endif ! lsimple


cccccccccccccccccccccc
c Boxes
cccccccccccccccccccccc
      if(lbox.and.lZ4l) then
      
      GaugeOK(1) = .true.

        tempvirt(0:3,3:4) = 0d0
c q > gg> VV> q
        ngauge=3 !ngaugeAbe
        comp=1
        do q=3,4
          call boxlineABETotal(qdia(0,1),gdia(0,3),kdia(0,3),qdia(0,2),
     &          psi(1,qhel,2),psi(1,qhel,1),epsG(0,3),epsV(0,qhel,q,3,lephel,1),
     &          qhel,musq,1,2,ngauge,comp,pengauge(1,1),tempres,pengaugeb(1,1),tempresb,div)
          tempvirt(0,q) = tempvirt(0,q) + tempresb
          tempvirt(1,q) = tempvirt(1,q) + tempres(1)
          tempvirt(2,q) = tempvirt(2,q) + tempres(2)
          ngauge=9
          comp=-1
        enddo

c abe gauge test
        call box_gauge(qdia(0,1),gdia(0,3),kdia(0,3),qdia(0,2),
     &         musq,pengauge(1,1),pengaugeb(1,1),wardidtest(1))

c for nAbe gaugetest
        momtempc(:3) = gdia(0:3,3)
        call boxlineABETotal(qdia(0,1),gdia(0,3),kdia(0,3),qdia(0,2),
     &          psi(1,qhel,2),psi(1,qhel,1),momtempc(0),epsV(0,qhel,3,3,lephel,1),
     &          qhel,musq,1,2,9,-1,nothing(1),checkg(1,1),nothing(1),checkgb(1),div)
c q > VV> gg> q
        ngauge=3 !ngaugeAbe
        comp=1
        do q=3,4
          call boxlineABETotal(qdia(0,1),kdia(0,3),gdia(0,3),qdia(0,2),
     &          psi(1,qhel,2),psi(1,qhel,1),epsV(0,qhel,q,3,lephel,1),epsG(0,3),
     &          qhel,musq,1,3,ngauge,comp,pengauge(1,2),tempres,pengaugeb(1,2),tempresb,div)
          tempvirt(0,q) = tempvirt(0,q) + tempresb
          tempvirt(1,q) = tempvirt(1,q) + tempres(1)
          tempvirt(2,q) = tempvirt(2,q) + tempres(2)
          ngauge=9
          comp=-1
        enddo

c abe gauge test
        call box_gauge(qdia(0,1),kdia(0,3),gdia(0,3),qdia(0,2),
     &         musq,pengauge(1,2),pengaugeb(1,2),wardidtest(2))

c for nAbe gaugetest
        call boxlineABETotal(qdia(0,1),kdia(0,3),gdia(0,3),qdia(0,2),
     &          psi(1,qhel,2),psi(1,qhel,1),epsV(0,qhel,3,3,lephel,1),momtempc(0),
     &          qhel,musq,1,3,9,-1,nothing(1),checkg(1,2),nothing(1),checkgb(2),div)


c sum up boxline abe
c colstr, epsG(:,3) = [ta,tb] = colstr2
c signconvention boxline:-Pestkin convention
        do q=3,4
          ampb(1,q) = ampb(1,q) + tempvirt(0,q)
          ampv(1,q) = ampv(1,q) + tempvirt(1,q)*CF
          ampv(1,q) = ampv(1,q) + tempvirt(2,q)*CFCA
          if(div.eq.0) ampv(1,q) = ampv(1,q) - tempvirt(0,q)*CF
          ampb(2,q) = ampb(2,q) - tempvirt(0,q)
          ampv(2,q) = ampv(2,q) - tempvirt(1,q)*CF
          ampv(2,q) = ampv(2,q) - tempvirt(2,q)*CFCA
          if(div.eq.0) ampv(2,q) = ampv(2,q) + tempvirt(0,q)*CF
        enddo

c non Abe Boxline
c q > VV> q
        ngauge=ngaugeNoAbe
        comp=1
        ngauge=2
        do q=3,4
          call BoxlineNoAbeTotal(qdia(0,1),kdia(0,3),qdia(0,2),gdia(0,3),
     &          psi(1,qhel,2),psi(1,qhel,1),epsV(0,qhel,q,3,lephel,1),epsG(0,3),
     &          qhel,musq,ngauge,comp,pengauge(1,3),tempres,div)
!             print*,q,checkg(1,3)
          ampv(1,q) = ampv(1,q) + tempres(1)*CA
          ampv(2,q) = ampv(2,q) - tempres(1)*CA
          ngauge=9
          comp=-1
        enddo
c gauge Test
          call BoxlineNoAbeTotal(qdia(0,1),kdia(0,3),qdia(0,2),gdia(0,3),
     &       psi(1,qhel,2),psi(1,qhel,1),epsV(0,qhel,3,3,lephel,1),momtempc(0),
     &          qhel,musq,9,-1,nothing(1),checkg(1,3),div)

        if (abs(pengauge(1,3)).lt.1d-8) then
          td=abs(-(checkg(2,2)+checkg(2,1))/2d0+pengauge(1,3))
        else
          td=abs(-(checkg(2,2)+checkg(2,1))/2d0/pengauge(1,3)+1d0)
        endif
        wardidtest(3) = td.lt.GaugeAccuracy

        gaugeOK(1) = gaugeOK(1).and.wardidtest(1).and.wardidtest(2).and.wardidtest(3).or.ldebug

      boxcount2m = boxcount2m + 1
      boxcount2sub(sub_number) = boxcount2sub(sub_number) + 1
      if(.not.gaugeOK(1)) then
        gaugetestfailed=.true.
        boxcountm   = boxcountm   + 1
        boxcountm_q = boxcountm_q + 1
        boxcountsub  (sub_number) = boxcountsub  (sub_number) + 1
        boxcountsub_q(sub_number) = boxcountsub_q(sub_number) + 1
      endif

      endif ! box

cccccccccccccccccccccc
c Pentagons
cccccccccccccccccccccc
      if(lpent) then


      nGaugeAbe=4
      nGaugeNoAbe=4
      
      GaugeOK(1) = .true.
      GaugeOK(2) = .true.

      if(lAZ(0,1).and.lAZ(0,2).and.lZZ) then

      do iV1=1,2
        iV2 =3-iV1

        tempvirt(0:3,3:4) = 0d0
        idia=1
c q > gg> V> V> q
        ngauge=ngaugeAbe
        comp=1
        do q=3,4
          call penlineAbeTotalQ(qdia(0,1),gdia(0,3),kdia(0,iV1),kdia(0,iV2),qdia(0,2),
     &          psi(1,qhel,2),psi(1,qhel,1),epsG(0,3),epsV(0,qhel,q,iV1,lephel,1),epsV(0,qhel,q,iV2,lephel,1),
     &          qhel,musq,1,2,ngauge,comp,pengauge(1,idia),tempres,pengaugeb(1,idia),tempresb,div)
          tempvirt(0,q) = tempvirt(0,q) + tempresb
          tempvirt(1,q) = tempvirt(1,q) + tempres(1)
          tempvirt(2,q) = tempvirt(2,q) + tempres(2)
          wardidtest(idia)=wardpAbe(1)
          wardidtest_Q(idia)=wardpAbe(2)    
          ngauge=9
          comp=-1
        enddo
        idia = idia+1

c q > V> gg> V> q
        ngauge=ngaugeAbe
        comp=1
        do q=3,4
          call penlineAbeTotalQ(qdia(0,1),kdia(0,iV1),gdia(0,3),kdia(0,iV2),qdia(0,2),
     &          psi(1,qhel,2),psi(1,qhel,1),epsV(0,qhel,q,iV1,lephel,1),epsG(0,3),epsV(0,qhel,q,iV2,lephel,1),
     &          qhel,musq,1,3,ngauge,comp,pengauge(1,idia),tempres,pengaugeb(1,idia),tempresb,div)
          tempvirt(0,q) = tempvirt(0,q) + tempresb
          tempvirt(1,q) = tempvirt(1,q) + tempres(1)
          tempvirt(2,q) = tempvirt(2,q) + tempres(2)
          wardidtest(idia)=wardpAbe(1)
          wardidtest_Q(idia)=wardpAbe(2)    
          ngauge=9
          comp=-1
        enddo
        idia = idia+1

c q > V> V> gg> q
        ngauge=ngaugeAbe
        comp=1
        do q=3,4
          call penlineAbeTotalQ(qdia(0,1),kdia(0,iV1),kdia(0,iV2),gdia(0,3),qdia(0,2),
     &          psi(1,qhel,2),psi(1,qhel,1),epsV(0,qhel,q,iV1,lephel,1),epsV(0,qhel,q,iV2,lephel,1),epsG(0,3),
     &          qhel,musq,1,4,ngauge,comp,pengauge(1,idia),tempres,pengaugeb(1,idia),tempresb,div)
          tempvirt(0,q) = tempvirt(0,q) + tempresb
          tempvirt(1,q) = tempvirt(1,q) + tempres(1)
          tempvirt(2,q) = tempvirt(2,q) + tempres(2)
          wardidtest(idia)=wardpAbe(1)
          wardidtest_Q(idia)=wardpAbe(2)    
          ngauge=9
          comp=-1
        enddo
        idia = idia+1


c sum up penline abe
c colstr, epsG(:,3) = [ta,tb] = colstr2
c signconvention penline:-Pestkin convention
        do q=3,4
          ampb(1,q) = ampb(1,q) - tempvirt(0,q)
          ampv(1,q) = ampv(1,q) - tempvirt(1,q)*CF
          ampv(1,q) = ampv(1,q) - tempvirt(2,q)*CFCA
          if(div.eq.0) ampv(1,q) = ampv(1,q) + tempvirt(0,q)*CF
          ampb(2,q) = ampb(2,q) + tempvirt(0,q)
          ampv(2,q) = ampv(2,q) + tempvirt(1,q)*CF
          ampv(2,q) = ampv(2,q) + tempvirt(2,q)*CFCA
          if(div.eq.0) ampv(2,q) = ampv(2,q) - tempvirt(0,q)*CF
        enddo

        tempvirt = 0d0

c non Abe Penline
c q > V> V> q
        ngauge=ngaugeNoAbe
        comp=1
        do q=3,4
          call penlineNoAbeTotalQ(qdia(0,1),kdia(0,iV1),kdia(0,iV2),qdia(0,2),gdia(0,3),
     &          psi(1,qhel,2),psi(1,qhel,1),epsV(0,qhel,q,iV1,lephel,1),epsV(0,qhel,q,iV2,lephel,1),epsG(0,3),
     &          qhel,musq,4,comp,pengauge(1,idia),tempres,div)
          ampv(1,q) = ampv(1,q) - tempres(1)*CA
          ampv(2,q) = ampv(2,q) + tempres(1)*CA
          wardidtest(idia)=wardpNoAbe(1)
          wardidtest_Q(idia)=wardpNoAbe(2)  
          ngauge=9
          comp=-1
        enddo
        idia = idia+1

        GaugeOK(1) =GaugeOK(1).and.  wardidtest(1).and.wardidtest(2).and.wardidtest(3).and.wardidtest(4).or.ldebug
        GaugeOK(2) =GaugeOK(2).and.  wardidtest_Q(1).and.wardidtest_Q(2).and.wardidtest_Q(3).and.wardidtest_Q(4).or.ldebug

        if(ldebug) then
          print*,'Gauge Test PenLine, glu1,iV1=',glu1,iV1
          print*,wardidtest(1:4)
          print*,wardidtest_Q(1:4)
        endif


      enddo ! iV1

      endif ! lAZ(0,1:2)
      
      if(lZ4l) then
      
      do glu1=1,2
      glu2=3-glu1
      
      tempvirt=0d0


      idia=1
c q > g> g > V> q
        ngauge=ngaugeAbe
        comp=1
        do q=3,4
          call penlineAbeTotalQ(qdia(0,1),gdia(0,glu1),gdia(0,glu2),kdia(0,3),qdia(0,2),
     &          psi(1,qhel,2),psi(1,qhel,1),epsG(0,glu1),epsG(0,glu2),epsV(0,qhel,q,3,lephel,1),
     &          qhel,musq,2,4,ngauge,comp,pengauge(1,idia),tempres,pengaugeb(1,idia),tempresb,div)
          tempvirt(0,q) = tempvirt(0,q) + tempresb
          tempvirt(1,q) = tempvirt(1,q) + tempres(1)
          tempvirt(2,q) = tempvirt(2,q) + tempres(2)
          wardidtest(idia)=wardpAbe(1)
          wardidtest_Q(idia)=wardpAbe(2)    
          ngauge=9
          comp=-1
        enddo
        idia = idia+1

c q > g> V > g> q
        ngauge=ngaugeAbe
        comp=1
        do q=3,4
          call penlineAbeTotalQ(qdia(0,1),gdia(0,glu1),kdia(0,3),gdia(0,glu2),qdia(0,2),
     &          psi(1,qhel,2),psi(1,qhel,1),epsG(0,glu1),epsV(0,qhel,q,3,lephel,1),epsG(0,glu2),
     &          qhel,musq,2,3,ngauge,comp,pengauge(1,idia),tempres,pengaugeb(1,idia),tempresb,div)
          tempvirt(0,q) = tempvirt(0,q) + tempresb
          tempvirt(1,q) = tempvirt(1,q) + tempres(1)
          tempvirt(2,q) = tempvirt(2,q) + tempres(2)
          wardidtest(idia)=wardpAbe(1)
          wardidtest_Q(idia)=wardpAbe(2)    
          ngauge=9
          comp=-1
        enddo
        idia = idia+1

c q > V> g> g > q
        ngauge=ngaugeAbe
        comp=1
        do q=3,4
          call penlineAbeTotalQ(qdia(0,1),kdia(0,3),gdia(0,glu1),gdia(0,glu2),qdia(0,2),
     &          psi(1,qhel,2),psi(1,qhel,1),epsV(0,qhel,q,3,lephel,1),epsG(0,glu1),epsG(0,glu2),
     &          qhel,musq,2,2,ngauge,comp,pengauge(1,idia),tempres,pengaugeb(1,idia),tempresb,div)
          tempvirt(0,q) = tempvirt(0,q) + tempresb
          tempvirt(1,q) = tempvirt(1,q) + tempres(1)
          tempvirt(2,q) = tempvirt(2,q) + tempres(2)
          wardidtest(idia)=wardpAbe(1)
          wardidtest_Q(idia)=wardpAbe(2)    
          ngauge=9
          comp=-1
        enddo
        idia = idia+1

c sum up penline abe
c signconvention penline:-Pestkin convention
        do q=3,4
          ampb(glu1,q) = ampb(glu1,q) + tempvirt(0,q)
          if(div.eq.0) ampv(glu1,q) = ampv(glu1,q) - tempvirt(0,q)*CF
          ampv(glu1,q) = ampv(glu1,q) + tempvirt(1,q)
          ampv(3,q) = ampv(3,q) + tempvirt(2,q)
        enddo
        
        
        GaugeOK(1) =GaugeOK(1).and.  wardidtest(1).and.wardidtest(2).and.wardidtest(3).or.ldebug
        GaugeOK(2) =GaugeOK(2).and.  wardidtest_Q(1).and.wardidtest_Q(2).and.wardidtest_Q(3).or.ldebug

        tempvirt = 0d0

      idia=1
c non Abe Penline
c q > V> g> q
        ngauge=ngaugeNoAbe
        comp=1
        do q=3,4
          call penlineNoAbe_pTQ(qdia(0,1),kdia(0,3),gdia(0,glu1),qdia(0,2),gdia(0,glu2),
     &          psi(1,qhel,2),psi(1,qhel,1),epsV(0,qhel,q,3,lephel,1),epsG(0,glu1),epsG(0,glu2),
     &          qhel,musq,2,3,ngauge,comp,pengauge(1,idia),tempres,div)
          tempvirt(1,q) = tempvirt(1,q) + tempres(1)
          tempvirt(2,q) = tempvirt(2,q) + tempres(2)
          tempvirt(3,q) = tempvirt(3,q) + tempres(3)
          wardidtest(idia)=wardpNoAbe(1)
          wardidtest_Q(idia)=wardpNoAbe(2)  
          ngauge=9
          comp=-1
        enddo
        idia = idia+1

c non Abe Penline
c q > g> V> q
        ngauge=ngaugeNoAbe
        comp=1
        do q=3,4
          call penlineNoAbe_pTQ(qdia(0,1),gdia(0,glu1),kdia(0,3),qdia(0,2),gdia(0,glu2),
     &          psi(1,qhel,2),psi(1,qhel,1),epsG(0,glu1),epsV(0,qhel,q,3,lephel,1),epsG(0,glu2),
     &          qhel,musq,2,2,ngauge,comp,pengauge(1,idia),tempres,div)
          tempvirt(1,q) = tempvirt(1,q) + tempres(1)
          tempvirt(2,q) = tempvirt(2,q) + tempres(2)
          tempvirt(3,q) = tempvirt(3,q) + tempres(3)
          wardidtest(idia)=wardpNoAbe(1)
          wardidtest_Q(idia)=wardpNoAbe(2)  
          ngauge=9
          comp=-1
        enddo
        idia = idia+1

c sum up penline NoAbe
c signconvention penline:-Pestkin convention
        do q=3,4
          ampv(glu1,q) = ampv(glu1,q) + tempvirt(2,q)*(0d0,-1.5d0)
          ampv(glu2,q) = ampv(glu2,q) + tempvirt(3,q)*(0d0,-1.5d0)
          ampv(3,q) = ampv(3,q) + tempvirt(1,q)*(0d0,0.25d0)
        enddo

c quad vertex Penline
        ngauge=1
        comp=1
        do q=3,4
          call NoAbe4vertBoxLineTotalQ(qdia(0,1),kdia(0,3),qdia(0,2),gdia(0,glu2),gdia(0,glu1),
     &          psi(1,qhel,2),psi(1,qhel,1),epsV(0,qhel,q,3,lephel,1),epsG(0,glu2),epsG(0,glu1),
     &          qhel,musq,2,0,comp,nothing(1),tempres(1),gaugereal(1,2),div,ngauge)
          ampv(glu1,q) = ampv(glu1,q) - 1.5d0*(tempres(2)+tempres(3)-2d0*tempres(1))/2d0
          ampv(glu2,q) = ampv(glu2,q) - 1.5d0*(tempres(1)+tempres(2)-2d0*tempres(3))/2d0
          ampv(3,q) = ampv(3,q)  -0.25d0*(2d0*tempres(2)-tempres(1)-tempres(3))/2d0
          wardidtest(idia)=wardpNoAbe(1)
          wardidtest_Q(idia)=wardpNoAbe(2)  
          ngauge=9
          comp=-1
        enddo
        idia = idia+1

c 2 3glu vert Penline
        ngauge=1
        comp=1
        do q=3,4
          call NoAbe2gaPenLineTotalQ(qdia(0,1),kdia(0,3),qdia(0,2),gdia(0,glu2),gdia(0,glu1),
     &          psi(1,qhel,2),psi(1,qhel,1),epsV(0,qhel,q,3,lephel,1),epsG(0,glu2),epsG(0,glu1),
     &          qhel,musq,2,0,comp,nothing(1),tempres(1),gaugereal(1,2),div,ngauge)
          ampv(glu1,q) = ampv(glu1,q) + 1.5d0*tempres(1)
          ampv(3,q) = ampv(3,q) + 0.25d0*tempres(1)
          wardidtest(idia)=wardpNoAbe(1)
          wardidtest_Q(idia)=wardpNoAbe(2)  
          ngauge=9
          comp=-1
        enddo

        GaugeOK(1) =GaugeOK(1).and.  wardidtest(1).and.wardidtest(2).and.wardidtest(3).and.wardidtest(4).or.ldebug
        GaugeOK(2) =GaugeOK(2).and.  wardidtest_Q(1).and.wardidtest_Q(2).and.wardidtest_Q(3).and.wardidtest_Q(4).or.ldebug
! 
!         idia = idia+1
! 
! 
      enddo !glu1

      endif ! z4l


      pentcount2m = pentcount2m + 1
      pentcount2sub(sub_number) = pentcount2sub(sub_number) + 1
      if(.not.gaugeOK(1)) then
        pentcountm = pentcountm + 1
        pentcountsub(sub_number) = pentcountsub(sub_number) + 1
      endif
      if(.not.gaugeOK(2)) then
        pentcountm_Q = pentcountm_q + 1
        pentcountsub_Q(sub_number) = pentcountsub_Q(sub_number) + 1
        gaugetestfailed = .true.
      endif

      endif ! pent
!       lhex=.false.


c ---------------------------------------
c $$$$$$$$$$ HEXAGONS $$$$$$$$$$$$$$$$$$$
c ---------------------------------------
!       epsV(0:3,qhel,3,2,lephel,1) = kdia(0:3,2)
!       epsV(0:3,qhel,4,2,lephel,1) = kdia(0:3,2)
!       ampv=0d0

      if(lhex.and.lAZ(0,1).and.lAZ(0,2).and.lzz) then ! only include if we have V1 and V2
!       if(.true.) then

      ngaugeAbe=4*rn+1
      if(ngaugeAbe.eq.5) then
        ngaugeAbe=4
      endif
      ngaugeNoAbe=3*rn+1
      if(ngaugeNoAbe.eq.4) then
        ngaugeNoAbe=3
      endif
      
      GaugeOK(1) = .true.
      GaugeOK(2) = .true.
      
!       epsG(0:3,1) = gdia(0:3,1)
!       ampv=0d0
      
      
      do glu1=1,2
      do iV1=1,2
!       ampv=0d0
!       ampb=0d0

        iV2 = 3-iV1
        glu2= 3-glu1

        tempvirt(0:3,3:4) = 0d0
        idia=1
        
c glu1 glu2 V V      
        ngauge=ngaugeAbe
        comp=1
        do q=3,4
          call HexlineAbeTotalQ(qdia(0,1),gdia(0,glu1),gdia(0,glu2),kdia(0,iV1),kdia(0,iV2),qdia(0,2),
     &         psi(1,qhel,2),psi(1,qhel,1),epsG(0,glu1),epsG(0,glu2),epsV(0,qhel,q,iV1,lephel,1),epsV(0,qhel,q,iV2,lephel,1),
     &         qhel,musq,2,23,comp,tempres,tempresb,hexgauge(1,idia),hexgaugeb(idia),div,ngauge)
          tempvirt(0,q) = tempvirt(0,q) + tempresb
          tempvirt(1,q) = tempvirt(1,q) + tempres(1)
          tempvirt(2,q) = tempvirt(2,q) + tempres(2)
          wardidtest(idia)=wardhAbe(1)
          wardidtest_Q(idia)=wardhAbe(2)  
          ngauge=9
          comp=-1
        enddo
        idia=idia +1
        
c glu1 V glu2 V      
        ngauge=ngaugeAbe
        comp=1
        do q=3,4
          call HexlineAbeTotalQ(qdia(0,1),gdia(0,glu1),kdia(0,iV1),gdia(0,glu2),kdia(0,iV2),qdia(0,2),
     &         psi(1,qhel,2),psi(1,qhel,1),epsG(0,glu1),epsV(0,qhel,q,iV1,lephel,1),epsG(0,glu2),epsV(0,qhel,q,iV2,lephel,1),
     &         qhel,musq,2,24,comp,tempres,tempresb,hexgauge(1,idia),hexgaugeb(idia),div,ngauge)
          tempvirt(0,q) = tempvirt(0,q) + tempresb
          tempvirt(1,q) = tempvirt(1,q) + tempres(1)
          tempvirt(2,q) = tempvirt(2,q) + tempres(2)
          wardidtest(idia)=wardhAbe(1)
          wardidtest_Q(idia)=wardhAbe(2)  
          ngauge=9
          comp=-1
        enddo
        idia=idia +1
        
c glu1  V V glu2      
        ngauge=ngaugeAbe
        comp=1
        do q=3,4
          call HexlineAbeTotalQ(qdia(0,1),gdia(0,glu1),kdia(0,iV1),kdia(0,iV2),gdia(0,glu2),qdia(0,2),
     &         psi(1,qhel,2),psi(1,qhel,1),epsG(0,glu1),epsV(0,qhel,q,iV1,lephel,1),epsV(0,qhel,q,iV2,lephel,1),epsG(0,glu2),
     &         qhel,musq,2,25,comp,tempres,tempresb,hexgauge(1,idia),hexgaugeb(idia),div,ngauge)
          tempvirt(0,q) = tempvirt(0,q) + tempresb
          tempvirt(1,q) = tempvirt(1,q) + tempres(1)
          tempvirt(2,q) = tempvirt(2,q) + tempres(2)
          wardidtest(idia)=wardhAbe(1)
          wardidtest_Q(idia)=wardhAbe(2)  
          ngauge=9
          comp=-1
        enddo
        idia=idia +1
        
c V glu1 glu2 V
        ngauge=ngaugeAbe
        comp=1
        do q=3,4
          call HexlineAbeTotalQ(qdia(0,1),kdia(0,iV1),gdia(0,glu1),gdia(0,glu2),kdia(0,iV2),qdia(0,2),
     &         psi(1,qhel,2),psi(1,qhel,1),epsV(0,qhel,q,iV1,lephel,1),epsG(0,glu1),epsG(0,glu2),epsV(0,qhel,q,iV2,lephel,1),
     &         qhel,musq,2,34,comp,tempres,tempresb,hexgauge(1,idia),hexgaugeb(idia),div,ngauge)
          tempvirt(0,q) = tempvirt(0,q) + tempresb
          tempvirt(1,q) = tempvirt(1,q) + tempres(1)
          tempvirt(2,q) = tempvirt(2,q) + tempres(2)
          wardidtest(idia)=wardhAbe(1)
          wardidtest_Q(idia)=wardhAbe(2)  
          ngauge=9
          comp=-1
        enddo
        idia=idia +1
        
c V glu1 V glu2     
        ngauge=ngaugeAbe
        comp=1
        do q=3,4
          call HexlineAbeTotalQ(qdia(0,1),kdia(0,iV1),gdia(0,glu1),kdia(0,iV2),gdia(0,glu2),qdia(0,2),
     &         psi(1,qhel,2),psi(1,qhel,1),epsV(0,qhel,q,iV1,lephel,1),epsG(0,glu1),epsV(0,qhel,q,iV2,lephel,1),epsG(0,glu2),
     &         qhel,musq,2,35,comp,tempres,tempresb,hexgauge(1,idia),hexgaugeb(idia),div,ngauge)
          tempvirt(0,q) = tempvirt(0,q) + tempresb
          tempvirt(1,q) = tempvirt(1,q) + tempres(1)
          tempvirt(2,q) = tempvirt(2,q) + tempres(2)
          wardidtest(idia)=wardhAbe(1)
          wardidtest_Q(idia)=wardhAbe(2)  
          ngauge=9
          comp=-1
        enddo
        idia=idia +1
        
c V V glu1 glu2
        ngauge=ngaugeAbe
        comp=1
        do q=3,4
          call HexlineAbeTotalQ(qdia(0,1),kdia(0,iV1),kdia(0,iV2),gdia(0,glu1),gdia(0,glu2),qdia(0,2),
     &         psi(1,qhel,2),psi(1,qhel,1),epsV(0,qhel,q,iV1,lephel,1),epsV(0,qhel,q,iV2,lephel,1),epsG(0,glu1),epsG(0,glu2),
     &         qhel,musq,2,45,comp,tempres,tempresb,hexgauge(1,idia),hexgaugeb(idia),div,ngauge)
          tempvirt(0,q) = tempvirt(0,q) + tempresb
          tempvirt(1,q) = tempvirt(1,q) + tempres(1)
          tempvirt(2,q) = tempvirt(2,q) + tempres(2)
          wardidtest(idia)=wardhAbe(1)
          wardidtest_Q(idia)=wardhAbe(2)  
          ngauge=9
          comp=-1
        enddo
        idia=idia +1

c sum up hexline abe
        do q=3,4
          ampb(glu1,q) = ampb(glu1,q) + tempvirt(0,q)
          ampv(glu1,q) = ampv(glu1,q) + tempvirt(1,q)
          if(div.eq.0) ampv(glu1,q) = ampv(glu1,q) - CF*tempvirt(0,q) ! rational term in MSbar
          ampv(3,q) = ampv(3,q) + tempvirt(2,q) ! delta_ab colstr
        enddo
        

c nAbe Hexline
        tempvirt = 0d0
      
c glu1 V V
        ngauge=ngaugeNoAbe
        comp=1
        do q=3,4
          call HexlineNoAbe1gTotalQ(qdia(0,1),gdia(0,glu1),kdia(0,iV1),kdia(0,iV2),qdia(0,2),gdia(0,glu2),
     &            psi(1,qhel,2),psi(1,qhel,1),epsG(0,glu1),epsV(0,qhel,q,iV1,lephel,1),epsV(0,qhel,q,iV2,lephel,1),epsG(0,glu2),
     &            qhel,musq,2,2,comp,tempres,hexgauge(1,idia),div,ngauge)
          tempvirt(1,q) = tempvirt(1,q) + tempres(1)
          tempvirt(2,q) = tempvirt(2,q) + tempres(2)
          tempvirt(3,q) = tempvirt(3,q) + tempres(3)
          wardidtest(idia)=wardhNoAbe(1) 
          wardidtest_Q(idia)=wardhNoAbe(2)  
          ngauge=9
          comp=-1
        enddo
        idia=idia+1
      
c glu1 V V
        ngauge=ngaugeNoAbe
        comp=1
        do q=3,4
          call HexlineNoAbe1gTotalQ(qdia(0,1),kdia(0,iV1),gdia(0,glu1),kdia(0,iV2),qdia(0,2),gdia(0,glu2),
     &            psi(1,qhel,2),psi(1,qhel,1),epsV(0,qhel,q,iV1,lephel,1),epsG(0,glu1),epsV(0,qhel,q,iV2,lephel,1),epsG(0,glu2),
     &            qhel,musq,2,3,comp,tempres,hexgauge(1,idia),div,ngauge)
          tempvirt(1,q) = tempvirt(1,q) + tempres(1)
          tempvirt(2,q) = tempvirt(2,q) + tempres(2)
          tempvirt(3,q) = tempvirt(3,q) + tempres(3)
          wardidtest(idia)=wardhNoAbe(1) 
          wardidtest_Q(idia)=wardhNoAbe(2)  
          ngauge=9
          comp=-1
        enddo
        idia=idia+1
      
c glu1 V V
        ngauge=ngaugeNoAbe
        comp=1
        do q=3,4
          call HexlineNoAbe1gTotalQ(qdia(0,1),kdia(0,iV1),kdia(0,iV2),gdia(0,glu1),qdia(0,2),gdia(0,glu2),
     &            psi(1,qhel,2),psi(1,qhel,1),epsV(0,qhel,q,iV1,lephel,1),epsV(0,qhel,q,iV2,lephel,1),epsG(0,glu1),epsG(0,glu2),
     &            qhel,musq,2,4,comp,tempres,hexgauge(1,idia),div,ngauge)
          tempvirt(1,q) = tempvirt(1,q) + tempres(1)
          tempvirt(2,q) = tempvirt(2,q) + tempres(2)
          tempvirt(3,q) = tempvirt(3,q) + tempres(3)
          wardidtest(idia)=wardhNoAbe(1) 
          wardidtest_Q(idia)=wardhNoAbe(2)  
          ngauge=9
          comp=-1
        enddo
        idia=idia+1

      
c sum up HexLine noAbe
        do q=3,4
          ampv(3,q)    = ampv(3,q)    + (0d0,0.25d0)*tempvirt(1,q)
          ampv(glu1,q) = ampv(glu1,q) + (0d0,-1.5d0)*tempvirt(2,q)
          ampv(glu2,q) = ampv(glu2,q) + (0d0,-1.5d0)*tempvirt(3,q)
        enddo


c 2 3glu vert
        ngauge=3
        comp=1
        do q=3,4
          call NoAbe2gaHexLineTotalQ(qdia(0,1),kdia(0,iV1),kdia(0,iV2),qdia(0,2),gdia(0,glu2),gdia(0,glu1),
     &             psi(1,qhel,2),psi(1,qhel,1),epsV(0,qhel,q,iV1,lephel,1),epsV(0,qhel,q,iV2,lephel,1),epsG(0,glu2),epsG(0,glu1),
     &             qhel,musq,2,0,comp,ctemp,tempres,gaugereal(1,1),div,ngauge)
          wardidtest(idia)=wardhNoAbe(1)
          wardidtest_Q(idia)=wardhNoAbe(2)  
          ngauge=9
          comp=-1
          ampv(glu1,q) = ampv(glu1,q)  - 1.5d0*tempres(1)
          ampv(3,q)    = ampv(3,q)     -0.25d0*tempres(1)
        enddo
        idia=idia+1

c 4glu vert
        ngauge=3
        comp=1
        do q=3,4
          call NoAbe4vertHexLineTotalQ(qdia(0,1),kdia(0,iV1),kdia(0,iV2),qdia(0,2),gdia(0,glu2),gdia(0,glu1),
     &             psi(1,qhel,2),psi(1,qhel,1),epsV(0,qhel,q,iV1,lephel,1),epsV(0,qhel,q,iV2,lephel,1),epsG(0,glu2),epsG(0,glu1),
     &             qhel,musq,2,0,comp,ctemp,tempres,gaugereal(1,2),div,ngauge)
          wardidtest(idia)=wardhNoAbe(1)
          wardidtest_Q(idia)=wardhNoAbe(2)  
          
!           ngauge=0
!           comp=-1
          
          ampv(glu1,q) = ampv(glu1,q) + 1.5d0*(tempres(3)+tempres(2)-2d0*tempres(1))/2d0 ! factor 1/2 -> diagram calculated twice
          ampv(glu2,q) = ampv(glu2,q) + 1.5d0*(tempres(1)+tempres(2)-2d0*tempres(3))/2d0 ! factor 1/2 -> diagram calculated twice
          ampv(3,q)    = ampv(3,q)    +0.25d0*(2d0*tempres(2)-tempres(1)-tempres(3))/2d0
        enddo
        idia=idia+1
          
        if(ldebug) then
          print*,'Gauge Test HexLine, glu1,iV1=',glu1,iV1
          print*,wardidtest(1:11)
          print*,wardidtest_Q(1:11)
        endif

!         print*,glu1,iV1
!         print*,ampv(:,3)
!         pause
        
        GaugeOK(1) = GaugeOK(1).and.wardidtest(1).and.wardidtest(2).and.wardidtest(3).and.wardidtest(4).and.wardidtest(5).and.
     $             wardidtest(6).and.wardidtest(7).and.wardidtest(8).and.wardidtest(9).and.wardidtest(10).and.wardidtest(11)   
        GaugeOK(2) = GaugeOK(2).and.wardidtest_Q(1).and.wardidtest_Q(2).and.wardidtest_Q(3).and.wardidtest_Q(4)
     $               .and.wardidtest_Q(5).and.wardidtest_Q(6).and.wardidtest_Q(7).and.wardidtest_Q(8)
     $               .and.wardidtest_Q(9).and.wardidtest_Q(10).and.wardidtest_Q(11)


      enddo;enddo ! iV1 glu

      hexcount2m = hexcount2m + 1
      hexcount2sub(sub_number) = hexcount2sub(sub_number) + 1
      if(.not.gaugeOK(1)) then
        hexcountm = hexcountm + 1
        hexcountsub(sub_number) = hexcountsub(sub_number) + 1
      endif
      if(.not.gaugeOK(2)) then
        hexcountm_Q = hexcountm_q + 1
        hexcountsub_Q(sub_number) = hexcountsub_Q(sub_number) + 1
        gaugetestfailed = .true.
      endif

      endif !lHex

c---------------------------------
c       Counter terms, see 0908.4124
c---------------------------------
      res=0d0
      do q=3,4
      
      if(ldebug) print*,'results, q=',q

      if(div.eq.1) then
        i=NFl
        ampv(1,q) = ampv(1,q) -ampb(1,q)*2d0*(11d0/6d0*CA-2d0/3d0*TR*NFl)
        ampv(2,q) = ampv(2,q) -ampb(2,q)*2d0*(11d0/6d0*CA-2d0/3d0*TR*NFl)
      endif
      if(div.eq.0.and.lsimple) then ! 
        i=NFl
        ampv(1,q) = ampv(1,q) -ampb(1,q)*2d0*(11d0/6d0*CA-2d0/3d0*TR*NFl)*log(musq/mursq(1,1))
        ampv(2,q) = ampv(2,q) -ampb(2,q)*2d0*(11d0/6d0*CA-2d0/3d0*TR*NFl)*log(musq/mursq(1,1))
      endif 


c---------------------------------
c       I Operator
c---------------------------------
c qqgg=1234
      if(lCSI) then
      logmuosij(1,2) = log(musq/dotrr(p(0,id(1),1),p(0,id(3),1))/2d0)
      logmuosij(1,3) = log(musq/dotrr(p(0,id(1),1),p(0,id(2),1))/2d0)
      logmuosij(1,4) = log(musq/dotrr(p(0,id(1),1),p(0,id(4),1))/2d0)
      logmuosij(2,3) = log(musq/dotrr(p(0,id(3),1),p(0,id(2),1))/2d0)
      logmuosij(2,4) = log(musq/dotrr(p(0,id(3),1),p(0,id(4),1))/2d0)
      logmuosij(3,4) = log(musq/dotrr(p(0,id(2),1),p(0,id(4),1))/2d0)
      
      
       ampb(1,q) = ampborn(1,q)
       ampb(2,q) = ampborn(2,q)

      if(div.eq.2) then
        ampi(1) = ampb(1,q)*26d0/3d0
        ampi(2) = ampb(2,q)*26d0/3d0
        ampi(3) = 0d0
      elseif(div.eq.1) then
        ampi(1) = -(ampb(1,q)*(-45d0 + 2d0*NFL + logmuosij(1,2) - 9d0*logmuosij(1,3) - 9d0*logmuosij(2,4) -9d0*logmuosij(3,4)))/3d0
        ampi(2) = -(ampb(2,q)*(-45d0 + 2d0*NFL + logmuosij(1,2) - 9d0*logmuosij(1,4) - 9d0*logmuosij(2,3) -9d0*logmuosij(3,4)))/3d0
        ampi(3) = ((ampb(1,q) + ampb(2,q))*logmuosij(1,2) - ampb(2,q)*logmuosij(1,3) - ampb(1,q)*logmuosij(1,4) - 
     &              ampb(1,q)*logmuosij(2,3) - ampb(2,q)*logmuosij(2,4) + ampb(1,q)*logmuosij(3,4) + ampb(2,q)*logmuosij(3,4))/2d0
      else
        if(alphaDip.eq.1d0) then
          ampi(1) = -(ampb(1,q)*(-840 + 32*NFL + 104*Pi**2 + 9*logmuosij(1,2) + 3*logmuosij(1,2)**2 + 3*(-30 + NFL)*logmuosij(1,3) 
     &               - 27*logmuosij(1,3)**2 - 90*logmuosij(2,4) + 3*NFL*logmuosij(2,4) - 27*logmuosij(2,4)**2 - 99*logmuosij(3,4)
     &               + 6*NFL*logmuosij(3,4) - 27*logmuosij(3,4)**2))/18.
          ampi(2) = -(ampb(2,q)*(-840 + 32*NFL + 104*Pi**2 + 9*logmuosij(1,2) + 3*logmuosij(1,2)**2 + 3*(-30 + NFL)*logmuosij(1,4)-
     &                 27*logmuosij(1,4)**2 - 90*logmuosij(2,3) + 3*NFL*logmuosij(2,3) - 27*logmuosij(2,3)**2 - 99*logmuosij(3,4)
     &                  + 6*NFL*logmuosij(3,4) - 27*logmuosij(3,4)**2))/18.
          ampi(3) = (27*(ampb(1,q) + ampb(2,q))*logmuosij(1,2) + 9*(ampb(1,q) + ampb(2,q))*logmuosij(1,2)**2 + 
     &               (-30 + NFL)*ampb(2,q)*logmuosij(1,3) - 9*ampb(2,q)*logmuosij(1,3)**2 - 30*ampb(1,q)*logmuosij(1,4) + 
     &               NFL*ampb(1,q)*logmuosij(1,4) - 9*ampb(1,q)*logmuosij(1,4)**2 - 30*ampb(1,q)*logmuosij(2,3) + 
     &               NFL*ampb(1,q)*logmuosij(2,3) - 9*ampb(1,q)*logmuosij(2,3)**2 - 30*ampb(2,q)*logmuosij(2,4) + 
     &               NFL*ampb(2,q)*logmuosij(2,4) - 9*ampb(2,q)*logmuosij(2,4)**2 + 33*ampb(1,q)*logmuosij(3,4) - 
     &               2*NFL*ampb(1,q)*logmuosij(3,4) + 33*ampb(2,q)*logmuosij(3,4) - 2*NFL*ampb(2,q)*logmuosij(3,4) + 
     &               9*ampb(1,q)*logmuosij(3,4)**2 + 9*ampb(2,q)*logmuosij(3,4)**2)/36.
        else
          ampi(1) = -(ampb(1,q)*(-570 - 270*alphaDip + 20*NFL + 12*alphaDip*NFL + 104*Pi**2 - 6*(-45 + 2*NFL)*Log(alphaDip) + 
     &                156*Log(alphaDip)**2 + 9*logmuosij(1,2) + 3*logmuosij(1,2)**2 - 90*logmuosij(1,3) + 3*NFL*logmuosij(1,3) - 
     &                27*logmuosij(1,3)**2 - 90*logmuosij(2,4) + 3*NFL*logmuosij(2,4) - 27*logmuosij(2,4)**2 - 99*logmuosij(3,4) +
     &                6*NFL*logmuosij(3,4) - 27*logmuosij(3,4)**2))/18.
          ampi(2) = -(ampb(2,q)*(-570 - 270*alphaDip + 20*NFL + 12*alphaDip*NFL + 104*Pi**2 - 6*(-45 + 2*NFL)*Log(alphaDip) + 
     &                156*Log(alphaDip)**2 + 9*logmuosij(1,2) + 3*logmuosij(1,2)**2 - 90*logmuosij(1,4) + 3*NFL*logmuosij(1,4) - 
     &                27*logmuosij(1,4)**2 - 90*logmuosij(2,3) + 3*NFL*logmuosij(2,3) - 27*logmuosij(2,3)**2 - 99*logmuosij(3,4) + 
     &                6*NFL*logmuosij(3,4) - 27*logmuosij(3,4)**2))/18.
          ampi(3) = (27*(ampb(1,q) + ampb(2,q))*logmuosij(1,2) + 9*(ampb(1,q) + ampb(2,q))*logmuosij(1,2)**2 + 
     &                (-30 + NFL)*ampb(2,q)*logmuosij(1,3) - 9*ampb(2,q)*logmuosij(1,3)**2 - 30*ampb(1,q)*logmuosij(1,4) + 
     &                NFL*ampb(1,q)*logmuosij(1,4) - 9*ampb(1,q)*logmuosij(1,4)**2 - 30*ampb(1,q)*logmuosij(2,3) + 
     &                NFL*ampb(1,q)*logmuosij(2,3) - 9*ampb(1,q)*logmuosij(2,3)**2 - 30*ampb(2,q)*logmuosij(2,4) + 
     &                NFL*ampb(2,q)*logmuosij(2,4) - 9*ampb(2,q)*logmuosij(2,4)**2 + 33*ampb(1,q)*logmuosij(3,4) - 
     &                2*NFL*ampb(1,q)*logmuosij(3,4) + 33*ampb(2,q)*logmuosij(3,4) - 2*NFL*ampb(2,q)*logmuosij(3,4) + 
     &                9*ampb(1,q)*logmuosij(3,4)**2 + 9*ampb(2,q)*logmuosij(3,4)**2)/36.
        endif
      endif
!      else
!        ampi=0d0
      endif ! lCSI


c square amplitudes

      if(ldebug) then
        ortb(1) = (0d0,0.5d0)*(ampb(2,q)-ampb(1,q))
        ortb(2) = 0.5d0*(ampb(1,q)+ampb(2,q))
        ortb(3) = (ampb(1,q)+ampb(2,q))/6d0
      endif


      ort(1) = (0d0,0.5d0)*(ampborn(2,q)-ampborn(1,q))
      ort(2) = 0.5d0*(ampborn(1,q)+ampborn(2,q))
      ort(3) = (ampborn(1,q)+ampborn(2,q))/6d0

      ortv(1) = (0d0,0.5d0)*(ampv(2,q)-ampv(1,q))
      ortv(2) = 0.5d0*(ampv(1,q)+ampv(2,q))
      ortv(3) = ampv(3,q)+(ampv(1,q)+ampv(2,q))/6d0

      orti(1) = (0d0,0.5d0)*(ampi(2)-ampi(1))
      orti(2) = 0.5d0*(ampi(1)+ampi(2))
      orti(3) = ampi(3)+(ampi(1)+ampi(2))/6d0


      if(div.ne.0) then
        print*,id
        print*
        print*,'final'
        print*,'div',div
        print*,'amp born  ',ampborn(:,q)
        print*,'amp bornV ',ampb(:,q)
        print*,'ampi      ',ampi(:)
        print*,'ampv      ',ampv(:,q)
        print*,'ampv/ampb ',(ampv(i,q)/ampb(i,q),i=1,2),ampv(3,q)
        print*,'-ampv/ampi',(-ampv(i,q)/ampi(i),i=1,3)

        print*
        print*,'final'
        print*,'div',div
        print*,'ort born  ',ort
        print*,'ort bornV ',ortb
        print*,'orti      ',orti
        print*,'ortv      ',ortv
        print*,'ortv/ortb ',(ortv(i)/ortb(i),i=1,3)
        print*,'-ortv/orti',(-ortv(i)/orti(i),i=1,3)
      endif

      res(2,q) = dble(12d0*ort(1)*dconjg(ort (1))+20d0/3d0*ort(2)*dconjg(ort (2))+24d0*ort(3)*dconjg(ort (3)))
      res(1,q) = dble(12d0*ort(1)*dconjg(orti(1))+20d0/3d0*ort(2)*dconjg(orti(2))+24d0*ort(3)*dconjg(orti(3)))
      res(0,q) = dble(12d0*ort(1)*dconjg(ortv(1))+20d0/3d0*ort(2)*dconjg(ortv(2))+24d0*ort(3)*dconjg(ortv(3)))

      if(.not.lborn) res(2,q) = 0d0
      
      if(ldebug) then
        print*,'Re(virt*born)'
        print*,'scale   ',sqrt(musq)
        print*,'res Born',res(2,q)
        print*,'res I*B2',res(1,q)
        print*,'res V*B ',res(0,q)
        print*,'(I+V)*B)',res(0,q)+res(1,q)
        print*,'-V/I    ',-res(0,q)/res(1,q)
        print*,'canell  ',abs((res(1,q)+res(0,q))/res(1,q))
        print*,(res(0,q)-res(1,q))/res(2,q)
!       pause
      endif
      res(2,q) = res(2,q)*G**4
      res(1,q) = res(1,q)*G**4*als(1,1)/2d0/pi ! 2Re *als/4pi
      res(0,q) = res(0,q)*G**4*als(1,1)/2d0/pi 

      res(1,q) = res(0,q)+res(1,q) ! V+I
      res(0,q) = res(1,q)+res(2,q) ! V+I+B
      
      enddo


      if(IsNan(res,0d0,6,0)) then
        if(ldebug.or..not.gaugetestfailed) print*,"NaN in "__FILE__,__LINE__
        do i=0,2
          res(i,3) = 0d0
          res(i,4) = 0d0
        enddo
      endif
      if(gaugetestfailed) then
        res=0d0
      endif

      if(ldebug) then
        print*,'gaugetestfailed',gaugetestfailed
        print*,"hel",qhel,lephel,ghel1,ghel2
      endif
c      pause

      end


