
       
c 4q WW Amplitude
c    q1 >>> #### >>> q3 
c           ####          , WW
c           ####
c    q2 >>> #### >>> q4  
c   res(colstr, q1,q2,flav config=0,2)
c      colstr=0..2: for catani seymour: <|>,  <|T1T2|>,  <|T1T3|>
c                   for lo: col summed, and the 2 different color flows
c      q1/2=3,4 for u and d:
c      flavour: 0: same family (t and u channel)
c      flavour: 1: diff family (only t channel), q1 corresponds to qline q1>q3
c      flavour: 2: diff family (only u channel), q1 corresponds to qline q2>q3
      subroutine amp4qzzVirt(P,id,fsign,nlo,lephel,res)
      implicit none
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "cacheZZ.inc"
#include "common.inc"
      real*8 P(0:3,max_p,max_kin),res(0:2,3:4,3:4,0:2)
      integer id(4),fsign(4),lephel,nlo
      logical luseprecalc
      
      real*8 pbar(0:3,4),pdia(0:3,4),momTemp(0:4),kdia(0:3,3)
      complex*16 psitemp(2,-1:1,4),kdiac(0:3,3)
      complex*16 currij(0:3,-1:1,1:2,3:4), currVij(0:3,-1:1,3:4,1:2,3:4,3)
      complex*16 amp(2),ampborn(2,-1:1,-1:1,3:4,3:4),ampb(2,-1:1,-1:1,3:4,3:4)
      complex*16 ampv(2,-1:1,-1:1,3:4,3:4),ampvCross1(2,-1:1,-1:1,3:4,3:4)
      real*8 resb(3:4,3:4,0:2),resi(3:4,3:4,0:2),resv(3:4,3:4,0:2)
      complex*16 temp
      complex*16 tempCurr(0:3,-1:1,4),psi4(4,4)
      complex*16 dotcc,momtempc(0:4)

      complex*16 boxgauge(2),boxres(5),boxgaugeb(2),boxborn,nothing(3),checkg(5,3),checkgb(2)
      complex*16 penbamp(3),penbres,penbborn,tempborn,tempvirt,hexamp(9),hexres,hexborn,hexresG(2,2)
      complex*16 pengauge(3),penres(8),pengaugeb(3),penborn 
      real*8 qjj2,qsq(1:2,3:4,0:3),logmuosij(4,4),musq,dotrr,td
      real*8 squared(2,2),gaugeprecision
      logical lborn,lcsi,lhex,lpenl,lpenb,lboxl,lsimple,lfloop,ldebug
      logical GaugeOK(0:2),wardidtest(0:8),wardidtest_q(0:8),ltop,lhiggs,lZ4l,lZZ
      integer div, col(2),comp,ngauge
      real*8 CF,CA,CFCA,TR
      parameter( CF = 4d0/3d0, CA = 3d0, CFCA=-1d0/6d0, TR=0.5d0)

      
      integer i,j,mu,tempi4(4), lepkin,colstr,h1,h2,q,q1,q2
      integer ikin, in,out,iV,hel,iV2,in2,out2
      integer qin,qin2,qout,qout2,qline, bos

      logical wardpAbe(2),wardpNoAbe(2)
      common/wardpenAbe/wardpAbe
      common/wardpenNoAbe/wardpNoAbe

#include "VBFNLO/utilities/scales.inc"
      real*8  alphas, xmt, alfa, xmz, xmw,s2w, xmh
      common /bkopin/   alphas, xmt, alfa, xmz, xmw, s2w, xmh

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
      logical GaugeE,GaugeO
      common/FloopPenggzz/GaugeE,GaugeO     

      GaugeE=.false.
      GaugeO=.false.
      
      ikin=1
      lepkin = ikin
      bos=22
      if(ikin .gt.7) lepkin = 1
      
      ldebug=.false.
      ltop=nfl.ge.5
      
      div= 0
!       epsV(:,:,:,1:2,:,1) = 0d0 
!       epsV(:,:,:,3,:,1) = 0d0  
     
      lborn =nlo.eq.1
      lCSI = nlo.eq.2
      lsimple=nlo.eq.2
      lfloop=nlo.eq.4
      lboxl= nlo.eq.2
      lpenl= nlo.eq.2
      lpenb= nlo.eq.2
      lhex = nlo.eq.2
      
cMK       lborn=.false.
cMK       lCSI = .false.
cMK       lsimple=.false.
cMK       lfloop=.true.
cMK       lboxl= .false.
cMK       lpenl= .false.
cMK       lpenb= .false.
cMK       lhex = .false.


!       lborn=.false.
!       lCSI = .false.
!       lsimple=.false.
!       lfloop=.false.
!       lboxl= .true.
!       lpenl= .false.
!       lpenb= .false.
!       lhex = .false.

      lZZ =  lAZ(0,1)
      lZ4l = lAZ(0,3)
      lhiggs=nfl.ge.5
      
      musq = mursq(1,1)
      if(musq.le.0d0) then
        print*,"unreasonable musq"
      endif
      
      ampborn=0d0
      ampb=0d0
      ampv=0d0
      ampvCross1=0d0
      resb=0d0
      resi=0d0
      resv=0d0

c in   q1 >>> #### >>> q3  out
c             ####          , WZ
c             ####
c in2  q2 >>> #### >>> q4  out2

      do mu=0,3
        pbar(mu,1)=p(mu,id(1),ikin) ! change convetnion for pbar for call of psi0m
        pbar(mu,2)=p(mu,id(3),ikin)
        pbar(mu,3)=p(mu,id(2),ikin)
        pbar(mu,4)=p(mu,id(4),ikin)
        pdia(mu,1)=pbar(mu,1)*fsign(1)
        pdia(mu,2)=pbar(mu,3)*fsign(2)
        pdia(mu,3)=pbar(mu,2)*fsign(3)
        pdia(mu,4)=pbar(mu,4)*fsign(4)
      enddo

      if(.not.lusePreCalc) then
      tempi4(1) = fsign(1)
      tempi4(2) = fsign(3)
      tempi4(3) = fsign(2)
      tempi4(4) = fsign(4)
c get bra and ket 
      call psi0m(4,pbar,tempi4,psitemp)
      psi(:,:,id(1),ikin) = psitemp(:,:,1)
      psi(:,:,id(2),ikin) = psitemp(:,:,3)
      psi(:,:,-id(3),ikin) = psitemp(:,:,2)
      psi(:,:,-id(4),ikin) = psitemp(:,:,4)
      
!       psiV = 0d0

c contract q with V
      do in=1,2 ! incoming quarks
      do iV=3,1,-1 ! V1, V2, V12
      do q=3,4 ! u/d
      do hel=-1,1,2
      if(lAZ(0,iV)) then
        call ket2c(psi(1,hel,id(in),ikin),.true.,pdia(0,in),hel,momV(0,iV,lepkin),epsV(0,hel,q,iV,lephel,lepkin),
     $             psiV(1,hel,q,id(in),iV,ikin),mompsiV(0,id(in),iV,ikin))
        out = in+2
        call bra2c(psi(1,hel,-id(out),ikin),.true.,pdia(0,out),hel,momV(0,iV,lepkin),epsV(0,hel,q,iV,lephel,lepkin),
     $             psiV(1,hel,q,-id(out),iV,ikin),mompsiV(0,-id(out),iV,ikin))
c contract with 2nd V
        if(iV.le.2) then
          iV2= 3-iV
          if(lAZ(0,iV2)) then
            call ket2c(psiV(1,hel,q,id(in),iV,ikin),.false.,mompsiV(0,id(in),iV,ikin),hel,momV(0,iV2,lepkin),
     $                 epsV(0,hel,q,iV2,lephel,lepkin),psitemp(1,hel,1),momtemp)
            psiV(1:2,hel,q,id(in),3,ikin)= psiV(1:2,hel,q,id(in),3,ikin)+psitemp(:,hel,1)
            call bra2c(psiV(1,hel,q,-id(out),iV,ikin),.false.,mompsiV(0,-id(out),iV,ikin),hel,momV(0,iV2,lepkin),
     $                 epsV(0,hel,q,iV2,lephel,lepkin),psitemp(1,hel,1),momtemp)
            psiV(1:2,hel,q,-id(out),3,ikin)= psiV(1:2,hel,q,-id(out),3,ikin)+psitemp(:,hel,1)
          endif
        endif
      endif
      enddo;enddo;enddo;enddo
!       print*,psiV(1:2,-1:1:2,:,id(1),1:2,1)
!       print*,psiV(:,-1,3,id(1),2,ikin)

      endif ! usePreCalc
      
c build currents
      do in=1,2
      do out=3,4
c <q|mu|q>
        call curr(1,psi(1,-1,-id(out),ikin),psi(1,-1,id(in),ikin),currij(0,-1,in,out))
        do q=3,4
c <q|mu|Vq> + <qV|mu|q>
          do iV=1,2
          if(lAZ(0,iV)) then
            call curr(1,psi(1,-1,-id(out),ikin),psiV(1,-1,q,id(in),iV,ikin),tempCurr(0,-1,1))
            call curr(1,psiV(1,-1,q,-id(out),iV,ikin),psi(1,-1,id(in),ikin),tempCurr(0,-1,2))
            currVij(:,:,q,in,out,iV) = tempCurr(:,:,1) + tempCurr(:,:,2)
          endif
          enddo
c <q|mu|VVq> + <qVV|mu|q> + <qV|mu|Vq>
!           if(lAZ(0,3)) then ! should always be true
            call curr(1,psi(1,-1,-id(out),ikin),psiV(1,-1,q,id(in),3,ikin),tempCurr(0,-1,1))
            call curr(1,psiV(1,-1,q,-id(out),3,ikin),psi(1,-1,id(in),ikin),tempCurr(0,-1,2))
            currVij(:,:,q,in,out,3) = tempCurr(:,:,1) + tempCurr(:,:,2)
!           else
!             currVij(:,:,q,in,out,3) = 0d0
!           endif
          if(lAZ(0,1).and.lAZ(0,2)) then
            call curr(1,psiV(1,-1,q,-id(out),1,ikin),psiV(1,-1,q,id(in),2,ikin),tempCurr(0,-1,3))
            call curr(1,psiV(1,-1,q,-id(out),2,ikin),psiV(1,-1,q,id(in),1,ikin),tempCurr(0,-1,4))
            currVij(:,:,q,in,out,3) = currVij(:,:,q,in,out,3) + tempCurr(:,:,3) + tempCurr(:,:,4)
          endif
        enddo
        qsq(in,out,0) = qjj2(pdia(0,in),pdia(0,out))
        qsq(in,out,1) = qjj2(pdia(0,in),mompsiV(0,-id(out),1,ikin))
        qsq(in,out,2) = qjj2(pdia(0,in),mompsiV(0,-id(out),2,ikin))
        qsq(in,out,3) = qjj2(pdia(0,in),mompsiV(0,-id(out),3,ikin))
      enddo
      enddo
!       print*,currij(:,-1,3,2)/qsq(3,2,0)*G
!       print*,dotcc(currij(0,-1,3,2),currVij(0,-1,3,1,4,3))/qsq(3,2,0)*G**2
        
c calc amplitudes

      res=0d0

      do q1=3,4
      do q2=3,4
c sum over q helicities
      do h1=-1,1,2
      do h2=-1,1,2
      
      amp = 0d0
      out=3 ! fix final state
      out2=4
      do in=1,2
        in2=3-in
        colstr = in
!         print*,currVij(0,h2,q2,in2,out2,3)
        
        if(lAZ(0,1)) then
c in>V1>out, in2>V2>out2
          amp(colstr) = amp(colstr) + dotcc(currVij(0,h1,q1,in,out,1),currVij(0,h2,q2,in2,out2,2))/qsq(in,out,1)
c in>V2>out, in2>V1>out2
          amp(colstr) = amp(colstr) + dotcc(currVij(0,h1,q1,in,out,2),currVij(0,h2,q2,in2,out2,1))/qsq(in2,out2,1)
        endif
!         if(lAZ(0,3)) then ! should always be true
c in>V12>out, in2>out2
          amp(colstr) = amp(colstr) + dotcc(currVij(0,h1,q1,in,out,3),currij(0,h2,in2,out2))/qsq(in2,out2,0)
c in>out, in2>V12>out2
          amp(colstr) = amp(colstr) + dotcc(currij(0,h1,in,out),currVij(0,h2,q2,in2,out2,3))/qsq(in,out,0)
!         endif
        ampborn(colstr,h1,h2,q1,q2) = amp(colstr)
      enddo

c square amplitudes
      amp=amp*G**2

      amp(2) = -amp(2) ! antisymmetrize
      squared(1,1) = dble(amp(1)*dconjg(amp(1)))
      squared(2,2) = dble(amp(2)*dconjg(amp(2)))
      if(h1.eq.h2) then
        squared(1,2) = dble(amp(1)*dconjg(amp(2)))
      else
        squared(1,2) = 0d0
      endif
!       print*,amp

      res(0,q1,q2,1) = res(0,q1,q2,1)+2d0*squared(1,1)
      res(0,q1,q2,2) = res(0,q1,q2,2)+2d0*squared(2,2)
      if(q1.eq.q2) then
        res(0,q1,q2,0) = res(0,q1,q2,0)+2d0*(squared(1,1)+squared(2,2)-2d0/3d0*squared(1,2))
      endif


      resb(q1,q2,1) = resb(q1,q2,1) + 2d0*squared(1,1)
      resb(q1,q2,2) = resb(q1,q2,2) + 2d0*squared(2,2)
      if(q1.eq.q2) then
        resb(q1,q2,0) = resb(q1,q2,0) + 2d0*(squared(1,1)+squared(2,2)-2d0/3d0*squared(1,2))
      endif
      
      if (lCSI) then ! I operator

        do i=1,3
        do j=i+1,4
          logmuosij(i,j) = log(musq/dotrr(pbar(0,i),pbar(0,j))/2d0)
        enddo
        enddo

        if(div.eq.0) then
          resi(q1,q2,0) =resi(q1,q2,0)+ (als(1,1)*((-32*(-15 + 2*Pi**2)*(3*squared(1,1) - 2*squared(1,2) + 3*squared(2,2)))/27. + 
     -      (-((logmuosij(1,2)**2 - logmuosij(1,4)**2 - logmuosij(2,3)**2 + logmuosij(3,4)**2)*
     -            (3*squared(1,1) - 2*squared(1,2) - 21*squared(2,2))) + 
     -         2*(logmuosij(1,3)**2 - logmuosij(1,4)**2 - logmuosij(2,3)**2 + logmuosij(2,4)**2)*
     -          (3*squared(1,1) - 10*squared(1,2) + 3*squared(2,2)) + 
     -         8*(logmuosij(1,4)**2 + logmuosij(2,3)**2)*(3*squared(1,1) - 2*squared(1,2) + 3*squared(2,2)))/9. + 
     -      3*(-((logmuosij(1,2) - logmuosij(1,4) - logmuosij(2,3) + logmuosij(3,4))*(3*squared(1,1) 
     &           - 2*squared(1,2) - 21*squared(2,2)))/9. + 
     -         (8*(logmuosij(1,4) + logmuosij(2,3))*(3*squared(1,1) - 2*squared(1,2) + 3*squared(2,2)))/9. + 
     -         (-logmuosij(1,3) + logmuosij(1,4) + logmuosij(2,3) - logmuosij(2,4))*
     -          ((20*squared(1,2))/9. - (2*(squared(1,1) + squared(2,2)))/3.))))/(2.*Pi)
          resi(q1,q2,1) = resi(q1,q2,1) -(als(1,1)*(-480 + 64*Pi**2 + 9*logmuosij(1,2) + 3*logmuosij(1,2)**2
     &           - 18*logmuosij(1,3) - 6*logmuosij(1,3)**2 - 63*logmuosij(1,4) - 
     -       21*logmuosij(1,4)**2 - 63*logmuosij(2,3) - 21*logmuosij(2,3)**2 - 18*logmuosij(2,4) 
     &          - 6*logmuosij(2,4)**2 + 9*logmuosij(3,4) + 
     -       3*logmuosij(3,4)**2)*squared(1,1))/(18.*Pi)
          resi(q1,q2,2) = resi(q1,q2,2) -(als(1,1)*(-480 + 64*Pi**2 - 63*logmuosij(1,2) - 21*logmuosij(1,2)**2
     &           - 18*logmuosij(1,3) - 6*logmuosij(1,3)**2 + 9*logmuosij(1,4) + 
     -       3*logmuosij(1,4)**2 + 9*logmuosij(2,3) + 3*logmuosij(2,3)**2 - 18*logmuosij(2,4)
     &           - 6*logmuosij(2,4)**2 - 63*logmuosij(3,4) - 
     -       21*logmuosij(3,4)**2)*squared(2,2))/(18.*Pi)
        else if (div.eq.1) then
          resi(q1,q2,0) = resi(q1,q2,0)+ (als(1,1)*((16*(3*squared(1,1) - 2*squared(1,2) + 3*squared(2,2)))/3. + 
     -      2*(-((logmuosij(1,2) - logmuosij(1,4) - logmuosij(2,3)
     &           + logmuosij(3,4))*(3*squared(1,1) - 2*squared(1,2) - 21*squared(2,2)))/9. + 
     -         (8*(logmuosij(1,4) + logmuosij(2,3))*(3*squared(1,1) - 2*squared(1,2) + 3*squared(2,2)))/9. + 
     -         (-logmuosij(1,3) + logmuosij(1,4) + logmuosij(2,3) - logmuosij(2,4))*
     -          ((20*squared(1,2))/9. - (2*(squared(1,1) + squared(2,2)))/3.))))/(2.*Pi)
          resi(q1,q2,1) = resi(q1,q2,1) -(als(1,1)*(-24 + logmuosij(1,2) - 2*logmuosij(1,3) 
     &          - 7*logmuosij(1,4) - 7*logmuosij(2,3) - 2*logmuosij(2,4) + logmuosij(3,4))*
     -     squared(1,1))/(3.*Pi)
          resi(q1,q2,2) = resi(q1,q2,2)+ (als(1,1)*(24 + 7*logmuosij(1,2) + 2*logmuosij(1,3)
     &           - logmuosij(1,4) - logmuosij(2,3) + 2*logmuosij(2,4) + 7*logmuosij(3,4))*
     -    squared(2,2))/(3.*Pi)
        else if (div.eq.2) then
          resi(q1,q2,0) =resi(q1,q2,0)+  (16*als(1,1)*(3*squared(1,1) - 2*squared(1,2) + 3*squared(2,2)))/(9.*Pi)
          resi(q1,q2,1) =resi(q1,q2,1)+  (16*als(1,1)*squared(1,1))/(3.*Pi)
          resi(q1,q2,2) =resi(q1,q2,2)+ (16*als(1,1)*squared(2,2))/(3.*Pi)
        else
          print*,"wrong div"
          stop
        endif
      endif ! I operator
          
c      add Counter term
      if(div.eq.1.and.lboxl) then 
        resv(q1,q2,1) = resv(q1,q2,1) - 4d0*als(1,1)/(4d0*pi)*(11d0/6d0*CA-2d0/3d0*TR*NFL)*2d0*squared(1,1)  ! 2 *CT_Paco,because there are CT coming from upper and lower line boxes; ohter factors 2 coming from 2Re[M M] and color
        resv(q1,q2,2) = resv(q1,q2,2) - 4d0*als(1,1)/(4d0*pi)*(11d0/6d0*CA-2d0/3d0*TR*NFL)*2d0*squared(2,2)  ! 2 *CT_Paco,because there are CT coming from upper and lower line boxes; ohter factors 2 coming from 2Re[M M] and color
        resv(q1,q2,0) = resv(q1,q2,0) - 4d0*als(1,1)/(4d0*pi)*(11d0/6d0*CA
     &          -2d0/3d0*TR*NFL)*2d0*(squared(1,1)+squared(2,2)-2d0/3d0*squared(1,2))
        if(ltop) then
          resv(q1,q2,1) = resv(q1,q2,1) - 4d0*als(1,1)/(4d0*pi)*(-2d0/3d0*TR)*2d0*squared(1,1)  ! 2 *CT_Paco,because there are CT coming from upper and lower line boxes; ohter factors 2 coming from 2Re[M M] and color
          resv(q1,q2,2) = resv(q1,q2,2) - 4d0*als(1,1)/(4d0*pi)*(-2d0/3d0*TR)*2d0*squared(2,2)  ! 2 *CT_Paco,because there are CT coming from upper and lower line boxes; ohter factors 2 coming from 2Re[M M] and color
          resv(q1,q2,0) = resv(q1,q2,0) - 4d0*als(1,1)/(4d0*pi)*(
     &          -2d0/3d0*TR)*2d0*(squared(1,1)+squared(2,2)-2d0/3d0*squared(1,2))
        endif
      endif
      if(div.eq.0.and.lboxl.and.ltop) then !finite part coming from onshell renormalization of top quark
        resv(q1,q2,1) = resv(q1,q2,1) + als(1,1)/pi*
     &                    (-(11d0/6d0*CA-2d0/3d0*TR*NFL)*log(musq/mursq(1,1))+log(musq/tmass/tmass)/3d0)*2d0*squared(1,1)
        resv(q1,q2,2) = resv(q1,q2,2) + als(1,1)/pi*
     &                    (-(11d0/6d0*CA-2d0/3d0*TR*NFL)*log(musq/mursq(1,1))+log(musq/tmass/tmass)/3d0)*2d0*squared(2,2)
        resv(q1,q2,0) = resv(q1,q2,0) + als(1,1)/pi*(-(11d0/6d0*CA-2d0/3d0*TR*NFL)*log(musq/mursq(1,1))+log(musq/tmass/tmass)/3d0)
     &          *2d0*(squared(1,1)+squared(2,2)-2d0/3d0*squared(1,2))
      endif

      enddo; enddo ! hel
      enddo; enddo ! q1, q2


ccccccccccccccccccccccccccccccccccc
cccccc V I R T U A L S
ccccccccccccccccccccccccccccccccccc
      GaugeOK(0) = .true. 

c  change momentum convention -> all particles flowing in
      pdia(:,1:2)=pdia(:,1:2)
      pdia(:,3:4)=-pdia(:,3:4)
      
      psitemp(:,:,1)=psi(:,:,id(1),ikin)
      psitemp(:,:,2)=psi(:,:,id(2),ikin)
      psitemp(:,:,3)=psi(:,:,-id(3),ikin)
      psitemp(:,:,4)=psi(:,:,-id(4),ikin)

      kdia=-momV(:3,:,lepkin)
      kdiac = kdia

c include propagator in currij, currVij
      do in=1,2
      do out=3,4
      do h1=-1,1,2
        currij(:,h1,in,out)=-currij(:,h1,in,out)/qsq(in,out,0)
      do q=3,4
        if (lAZ(0,1)) currVij(:,h1,q,in,out,1)=-currVij(:,h1,q,in,out,1)/qsq(in,out,1)
        if (lAZ(0,2)) currVij(:,h1,q,in,out,2)=-currVij(:,h1,q,in,out,2)/qsq(in,out,2)
        if (lAZ(0,3)) currVij(:,h1,q,in,out,3)=-currVij(:,h1,q,in,out,3)/qsq(in,out,3)
      enddo;enddo;enddo;enddo

      out=3 ! fix final state
      out2=4
      do in=1,2
      in2= 3-in

      ampb=0d0
      ampv=0d0
      
      col(1) = in ! born color structure
      col(2) =  in2


c ---------------------------------------
c $$$$$$$$$$ S I M P L E $$$$$$$$$$$$$$$$
c ---------------------------------------

      if (lsimple) then
        
      do qline=1,2 ! corrections to both qlines

        if(qline.eq.1) then
          qin=in
          qout=out
          qin2=in2
          qout2=out2
        else
          qin=in2
          qout=out2
          qin2=in
          qout2=out
        endif
        
        do q1=3,4
        do h1=-1,1,2
        do h2=-1,1,2
c both Z on one qline
c vertline
          tempborn=dotcc(currVij(0,h1,q1,qin,qout,3),currij(0,h2,qin2,qout2))*qsq(qin,qout,3) ! don't include prop. twice
          call vertline(qsq(qin,qout,3),musq,1,tempvirt,div)
          tempvirt=tempborn*tempvirt
c glu self energy
          call gluonSelf(qsq(qin,qout,3),musq,NFL,temp,div)
          tempvirt=tempvirt+tempborn*temp
          if(ltop) then 
            call gluonSelfMass(qsq(qin,qout,3),musq,tmass,temp,div)
            tempvirt=tempvirt+tempborn*temp
          endif
          
          if(qline.eq.1) then
            ampb(col(1),h1,h2,q1,3) = ampb(col(1),h1,h2,q1,3) +tempborn*2d0 ! vertline and selfenergy
            ampv(col(1),h1,h2,q1,3) = ampv(col(1),h1,h2,q1,3) +tempvirt

            ampb(col(1),h1,h2,q1,4) = ampb(col(1),h1,h2,q1,4) +tempborn*2d0
            ampv(col(1),h1,h2,q1,4) = ampv(col(1),h1,h2,q1,4) +tempvirt
          else
            ampb(col(1),h2,h1,3,q1) = ampb(col(1),h2,h1,3,q1)  +tempborn*2d0
            ampv(col(1),h2,h1,3,q1) = ampv(col(1),h2,h1,3,q1)  +tempvirt

            ampb(col(1),h2,h1,4,q1) = ampb(col(1),h2,h1,4,q1) + tempborn*2d0
            ampv(col(1),h2,h1,4,q1) = ampv(col(1),h2,h1,4,q1) + tempvirt
          endif

c Z on different qline
         if(lZZ) then
          do q2=3,4
c glu self energy
          tempborn=dotcc(currVij(0,h1,q1,qin,qout,1),currVij(0,h2,q2,qin2,qout2,2))*qsq(qin,qout,1) ! don't include prop. twice
          call gluonSelf(qsq(qin,qout,1),musq,NFL,temp,div)
          tempvirt=tempborn*temp
          if(ltop) then 
            call gluonSelfMass(qsq(qin,qout,1),musq,tmass,temp,div)
            tempvirt=tempvirt+tempborn*temp
          endif
          
          if(qline.eq.1) then
            ampb(col(1),h1,h2,q1,q2) = ampb(col(1),h1,h2,q1,q2) +tempborn
            ampv(col(1),h1,h2,q1,q2) = ampv(col(1),h1,h2,q1,q2) +tempvirt
          else
            ampb(col(1),h2,h1,q2,q1) = ampb(col(1),h2,h1,q2,q1)  +tempborn
            ampv(col(1),h2,h1,q2,q1) = ampv(col(1),h2,h1,q2,q1)  +tempvirt
          endif
          enddo ! q2
         endif
        enddo; enddo; enddo
      enddo
      endif

c ---------------------------------------
c $$$$$$$$$$ B O X E S $$$$$$$$$$$$$$$$$$
c ---------------------------------------

      if (lboxl) then
        
      gaugeOK(1) = .true.
      do qline=1,2 ! corrections to both qlines
        if(ldebug) then
          print*,'boxline qline=',qline
        endif

        if(qline.eq.1) then
          qin=in
          qout=out
          qin2=in2
          qout2=out2
        else
          qin=in2
          qout=out2
          qin2=in
          qout2=out
        endif
        
        if(lZ4l) then
        
c q > Gij > Z4l > q
        comp=1
        ngauge=3
        momtemp(:3) = -pdia(:,qin)-pdia(:,qout)-kdia(:3,3)
        momtempc=momTemp
        if(ldebug) print*,qin,qout,qin2,qout2

        do q1=3,4
        do h1=-1,1,2
        do h2=-1,1,2
          call boxlineABETotal(pdia(0,qin),momtemp(0),kdia(0,3),pdia(0,qout),
     &            psitemp(1,h1,qout),psitemp(1,h1,qin),currij(0,h2,qin2,qout2),epsV(0,h1,q1,3,lephel,1),
     &            h1,musq,1,2,ngauge,comp,boxgauge(1),boxres(1),boxgaugeb(1),boxborn,div)
          comp =-1
          ngauge=9
          
          if(qline.eq.1) then
            ampb(col(1),h1,h2,q1,3) = ampb(col(1),h1,h2,q1,3) - boxborn
            ampv(col(1),h1,h2,q1,3) = ampv(col(1),h1,h2,q1,3) - boxres(1)*CF - boxres(2)*CFCA
            if( div.eq.0) ampv(col(1),h1,h2,q1,3) = ampv(col(1),h1,h2,q1,3) + boxborn*CF

            ampb(col(1),h1,h2,q1,4) = ampb(col(1),h1,h2,q1,4) - boxborn
            ampv(col(1),h1,h2,q1,4) = ampv(col(1),h1,h2,q1,4) - boxres(1)*CF - boxres(2)*CFCA
            if( div.eq.0) ampv(col(1),h1,h2,q1,4) = ampv(col(1),h1,h2,q1,4) + boxborn*CF
          else
            ampb(col(1),h2,h1,3,q1) = ampb(col(1),h2,h1,3,q1)  - boxborn
            ampv(col(1),h2,h1,3,q1) = ampv(col(1),h2,h1,3,q1)  - boxres(1)*CF - boxres(2)*CFCA
            if( div.eq.0) ampv(col(1),h2,h1,3,q1) = ampv(col(1),h2,h1,3,q1) + boxborn*CF

            ampb(col(1),h2,h1,4,q1) = ampb(col(1),h2,h1,4,q1) - boxborn
            ampv(col(1),h2,h1,4,q1) = ampv(col(1),h2,h1,4,q1) - boxres(1)*CF - boxres(2)*CFCA
            if( div.eq.0) ampv(col(1),h2,h1,4,q1) = ampv(col(1),h2,h1,4,q1) + boxborn*CF
          endif
          
        enddo; enddo; enddo

c for nAbe gauge test
      call boxlineABETotal(pdia(0,qin),momtemp(0),kdia(0,3),pdia(0,qout),
     &            psitemp(1,-1,qout),psitemp(1,-1,qin),momtempc(0),epsV(0,-1,3,3,lephel,1),
     &            -1,musq,1,2,9,-1,nothing(1),checkg(1,1),nothing(1),checkgb(1),div)
c gauge tests Abe
      call box_gauge(pdia(0,qin),momtemp(0),kdia(0,3),pdia(0,qout),musq,
     2    boxgauge(1),boxgaugeb(1),wardidtest(1))
      

c q > Z4l > currij > q
        comp=1
        ngauge=3

        do q1=3,4
        do h1=-1,1,2
        do h2=-1,1,2
          call boxlineABETotal(pdia(0,qin),kdia(0,3),momtemp(0),pdia(0,qout),
     &            psitemp(1,h1,qout),psitemp(1,h1,qin),epsV(0,h1,q1,3,lephel,1),currij(0,h2,qin2,qout2),
     &            h1,musq,1,3,ngauge,comp,boxgauge(1),boxres(1),boxgaugeb(1),boxborn,div)
          comp =-1
          ngauge=9

          if(qline.eq.1) then
            ampb(col(1),h1,h2,q1,3) = ampb(col(1),h1,h2,q1,3) - boxborn
            ampv(col(1),h1,h2,q1,3) = ampv(col(1),h1,h2,q1,3) - boxres(1)*CF - boxres(2)*CFCA
            if( div.eq.0) ampv(col(1),h1,h2,q1,3) = ampv(col(1),h1,h2,q1,3) + boxborn*CF

            ampb(col(1),h1,h2,q1,4) = ampb(col(1),h1,h2,q1,4) - boxborn 
            ampv(col(1),h1,h2,q1,4) = ampv(col(1),h1,h2,q1,4) - boxres(1)*CF  - boxres(2)*CFCA
            if( div.eq.0) ampv(col(1),h1,h2,q1,4) = ampv(col(1),h1,h2,q1,4)  + boxborn*CF
          else
            ampb(col(1),h2,h1,3,q1) = ampb(col(1),h2,h1,3,q1)  - boxborn 
            ampv(col(1),h2,h1,3,q1) = ampv(col(1),h2,h1,3,q1)  - boxres(1)*CF - boxres(2)*CFCA
            if( div.eq.0) ampv(col(1),h2,h1,3,q1) = ampv(col(1),h2,h1,3,q1) + boxborn*CF

            ampb(col(1),h2,h1,4,q1) = ampb(col(1),h2,h1,4,q1) - boxborn
            ampv(col(1),h2,h1,4,q1) = ampv(col(1),h2,h1,4,q1) - boxres(1)*CF - boxres(2)*CFCA
            if( div.eq.0) ampv(col(1),h2,h1,4,q1) = ampv(col(1),h2,h1,4,q1) + boxborn*CF
          endif
          
        enddo; enddo; enddo

c for nAbe gauge test
        call boxlineABETotal(pdia(0,qin),kdia(0,3),momtemp(0),pdia(0,qout),
     &            psitemp(1,-1,qout),psitemp(1,-1,qin),epsV(0,-1,3,3,lephel,1),momtempc(0),
     &            -1,musq,1,3,9,-1,nothing(1),checkg(1,2),nothing(1),checkgb(2),div)
c gauge tests Abe
        call box_gauge(pdia(0,qin),kdia(0,3),momtemp(0),pdia(0,qout),musq,
     2            boxgauge(1),boxgaugeb(1),wardidtest(2))
          
c boxline nAbe
c start with gauge test
        call BoxlineNoAbeTotal(pdia(0,qin),kdia(0,3),pdia(0,qout),momtemp(0),
     &             psitemp(1,-1,qout),psitemp(1,-1,qin),epsV(0,-1,3,3,lephel,1),momtempc(0),
     &             -1,musq,9,1,nothing,checkg(1,3),div)
        if (abs(checkg(1,3)).lt.1d-8) then
          td=abs(-(checkg(2,2)+checkg(2,1))/2d0+checkg(1,3))
        else
          td=abs(-(checkg(2,2)+checkg(2,1))/2d0/checkg(1,3)+1d0)
        endif
        wardidtest(3) = td.lt.GaugeAccuracy
          
        do q1=3,4
        do h1=-1,1,2
        do h2=-1,1,2
          call BoxlineNoAbeTotal(pdia(0,qin),kdia(0,3),pdia(0,qout),momtemp(0),
     &             psitemp(1,h1,qout),psitemp(1,h1,qin),epsV(0,h1,q1,3,lephel,1),currij(0,h2,qin2,qout2),
     &             h1,musq,9,-1,boxgauge(1),boxres(1),div)

          if(qline.eq.1) then
            ampv(col(1),h1,h2,q1,3) = ampv(col(1),h1,h2,q1,3) - boxres(1)*CA 
            ampv(col(1),h1,h2,q1,4) = ampv(col(1),h1,h2,q1,4) - boxres(1)*CA 
          else
            ampv(col(1),h2,h1,3,q1) = ampv(col(1),h2,h1,3,q1)  - boxres(1)*CA 
            ampv(col(1),h2,h1,4,q1) = ampv(col(1),h2,h1,4,q1) - boxres(1)*CA
          endif
          
        enddo; enddo; enddo
        gaugeOK(1) = gaugeOK(1).and.wardidtest(1).and.wardidtest(2).and.wardidtest(3).or.ldebug

        endif ! lZ4l
c q>V>q  q>V>q boxline corretions
        if(ldebug) print*,'boxline in>iV>out',qin,iV,qout
        if(lZZ) then
        do iV=1,2
        iV2=3-iV

        comp=1
        ngauge=3
        momtemp(:3) = -pdia(:,qin)-pdia(:,qout)-kdia(:3,iV)
        momtempc=momTemp

c qin > currVij(iV2) > iV1 > qout
        do q1=3,4
        do q2=3,4
        do h1=-1,1,2
        do h2=-1,1,2
          call boxlineABETotal(pdia(0,qin),momtemp(0),kdia(0,iV),pdia(0,qout),
     &            psitemp(1,h1,qout),psitemp(1,h1,qin),currVij(0,h2,q2,qin2,qout2,iV2),epsV(0,h1,q1,iV,lephel,1),
     &            h1,musq,1,2,ngauge,comp,boxgauge(1),boxres(1),boxgaugeb(1),boxborn,div)
          comp =-1
          ngauge=9
          
          if(qline.eq.1) then
            ampb(col(1),h1,h2,q1,q2) = ampb(col(1),h1,h2,q1,q2) - boxborn
            ampv(col(1),h1,h2,q1,q2) = ampv(col(1),h1,h2,q1,q2) - boxres(1)*CF - boxres(2)*CFCA
            if( div.eq.0) ampv(col(1),h1,h2,q1,q2) = ampv(col(1),h1,h2,q1,q2) + boxborn*CF
          else
            ampb(col(1),h2,h1,q2,q1) = ampb(col(1),h2,h1,q2,q1)  - boxborn
            ampv(col(1),h2,h1,q2,q1) = ampv(col(1),h2,h1,q2,q1)  - boxres(1)*CF - boxres(2)*CFCA
            if( div.eq.0) ampv(col(1),h2,h1,q2,q1) = ampv(col(1),h2,h1,q2,q1) + boxborn*CF
          endif
          
        enddo; enddo; enddo; enddo

c for nAbe gauge test
        call boxlineABETotal(pdia(0,qin),momtemp(0),kdia(0,iV),pdia(0,qout),
     &            psitemp(1,-1,qout),psitemp(1,-1,qin),momtempc(0),epsV(0,-1,3,iV,lephel,1),
     &            -1,musq,1,2,9,-1,nothing(1),checkg(1,1),nothing(1),checkgb(1),div)
c gauge tests Abe
        call box_gauge(pdia(0,qin),momtemp(0),kdia(0,iV),pdia(0,qout),musq,
     2    boxgauge(1),boxgaugeb(1),wardidtest(1))
      
c qin > iV1 >currVij(iV2) >  qout
        comp =1
        ngauge=3

        do q1=3,4
        do q2=3,4
        do h1=-1,1,2
        do h2=-1,1,2
          call boxlineABETotal(pdia(0,qin),kdia(0,iV),momtemp(0),pdia(0,qout),
     &            psitemp(1,h1,qout),psitemp(1,h1,qin),epsV(0,h1,q1,iV,lephel,1),currVij(0,h2,q2,qin2,qout2,iV2),
     &            h1,musq,1,3,ngauge,comp,boxgauge(1),boxres(1),boxgaugeb(1),boxborn,div)
          comp =-1
          ngauge=9
          
          if(qline.eq.1) then
            ampb(col(1),h1,h2,q1,q2) = ampb(col(1),h1,h2,q1,q2) - boxborn
            ampv(col(1),h1,h2,q1,q2) = ampv(col(1),h1,h2,q1,q2) - boxres(1)*CF - boxres(2)*CFCA
            if( div.eq.0) ampv(col(1),h1,h2,q1,q2) = ampv(col(1),h1,h2,q1,q2) + boxborn*CF
          else
            ampb(col(1),h2,h1,q2,q1) = ampb(col(1),h2,h1,q2,q1)  - boxborn
            ampv(col(1),h2,h1,q2,q1) = ampv(col(1),h2,h1,q2,q1)  - boxres(1)*CF - boxres(2)*CFCA
            if( div.eq.0) ampv(col(1),h2,h1,q2,q1) = ampv(col(1),h2,h1,q2,q1) + boxborn*CF
          endif
          
        enddo; enddo; enddo; enddo

c for nAbe gauge test
      call boxlineABETotal(pdia(0,qin),kdia(0,iV),momtemp(0),pdia(0,qout),
     &            psitemp(1,-1,qout),psitemp(1,-1,qin),epsV(0,-1,3,iV,lephel,1),momtempc(0),
     &            -1,musq,1,3,9,-1,nothing(1),checkg(1,2),nothing(1),checkgb(2),div)
c gauge tests Abe
      call box_gauge(pdia(0,qin),kdia(0,iV),momtemp(0),pdia(0,qout),musq,
     2    boxgauge(1),boxgaugeb(1),wardidtest(2))

c boxline nAbe
c start with gauge test
        call BoxlineNoAbeTotal(pdia(0,qin),kdia(0,iV),pdia(0,qout),momtemp(0),
     &             psitemp(1,-1,qout),psitemp(1,-1,qin),epsV(0,-1,3,iV,lephel,1),momtempc(0),
     &             -1,musq,9,1,nothing,checkg(1,3),div)

        if (abs(checkg(1,3)).lt.1d-8) then
          td=abs(-(checkg(2,2)+checkg(2,1))/2d0+checkg(1,3))
        else
          td=abs(-(checkg(2,2)+checkg(2,1))/2d0/checkg(1,3)+1d0)
        endif
        wardidtest(3) = td.lt.GaugeAccuracy
          
        do q1=3,4
        do q2=3,4
        do h1=-1,1,2
        do h2=-1,1,2
          call BoxlineNoAbeTotal(pdia(0,qin),kdia(0,iV),pdia(0,qout),momtemp(0),
     &             psitemp(1,h1,qout),psitemp(1,h1,qin),epsV(0,h1,q1,iV,lephel,1),currVij(0,h2,q2,qin2,qout2,iV2),
     &             h1,musq,9,-1,boxgauge(1),boxres(1),div)

          if(qline.eq.1) then
            ampv(col(1),h1,h2,q1,q2) = ampv(col(1),h1,h2,q1,q2) - boxres(1)*CA 
          else
            ampv(col(1),h2,h1,q2,q1) = ampv(col(1),h2,h1,q2,q1)  - boxres(1)*CA 
          endif
          
        enddo; enddo; enddo; enddo
        gaugeOK(1) = gaugeOK(1).and.wardidtest(1).and.wardidtest(2).and.wardidtest(3).or.ldebug
        if(ldebug) print*,'wardidtest',wardidtest(1:3)
      enddo ! iV
      endif ! lZZ
      enddo ! qline

      boxcount2m = boxcount2m + 1
      boxcount2sub(sub_number) = boxcount2sub(sub_number) + 1

      if(.not.gaugeOK(1)) then
        GaugeOK(0)=.false.
        boxcountm   = boxcountm   + 1
        boxcountm_q = boxcountm_q + 1
        boxcountsub  (sub_number) = boxcountsub  (sub_number) + 1
        boxcountsub_q(sub_number) = boxcountsub_q(sub_number) + 1
      endif

      endif ! boxl


c ---------------------------------------
c $$$$$$$$$$ Penline $$$$$$$$$$$$$$$$$$$$
c ---------------------------------------

      if(lpenl.and.lZZ) then
        
      gaugeOK(1) = .true.
      gaugeOK(2) = .true.
      do qline=1,2 ! corrections to both qlines
        if(ldebug) then
          print*,'penline qline=',qline
        endif

        if(qline.eq.1) then
          qin=in
          qout=out
          qin2=in2
          qout2=out2
        else
          qin=in2
          qout=out2
          qin2=in
          qout2=out
        endif
        momtemp(:3) = -pdia(:,qin)-pdia(:,qout)-kdia(:3,3)
        momtempc=momTemp
        
        do iV=1,2
        iV2=3-iV

c q > Gij > iV > iV2 > q
        comp=1
        ngauge=4
        if(ldebug) print*,qin,qout,qin2,qout2

        do q1=3,4
        do h1=-1,1,2
        do h2=-1,1,2
          call penlineABETotalQ(pdia(0,qin),momtemp(0),kdia(0,iV),kdia(0,iV2),pdia(0,qout),
     &            psitemp(1,h1,qout),psitemp(1,h1,qin),currij(0,h2,qin2,qout2),
     &            epsV(0,h1,q1,iV,lephel,1),epsV(0,h1,q1,iV2,lephel,1),
     &            h1,musq,1,2,ngauge,comp,pengauge(1),penres(1),pengaugeb(1),penborn,div)
          comp =-1
          ngauge=9
          
          if(qline.eq.1) then
            ampb(col(1),h1,h2,q1,3) = ampb(col(1),h1,h2,q1,3) + penborn
            ampv(col(1),h1,h2,q1,3) = ampv(col(1),h1,h2,q1,3) + penres(1)*CF + penres(2)*CFCA
            if( div.eq.0) ampv(col(1),h1,h2,q1,3) = ampv(col(1),h1,h2,q1,3) - penborn*CF

            ampb(col(1),h1,h2,q1,4) = ampb(col(1),h1,h2,q1,4) + penborn
            ampv(col(1),h1,h2,q1,4) = ampv(col(1),h1,h2,q1,4) + penres(1)*CF + penres(2)*CFCA
            if( div.eq.0) ampv(col(1),h1,h2,q1,4) = ampv(col(1),h1,h2,q1,4) - penborn*CF
          else
            ampb(col(1),h2,h1,3,q1) = ampb(col(1),h2,h1,3,q1) + penborn
            ampv(col(1),h2,h1,3,q1) = ampv(col(1),h2,h1,3,q1) + penres(1)*CF + penres(2)*CFCA
            if( div.eq.0) ampv(col(1),h2,h1,3,q1) = ampv(col(1),h2,h1,3,q1) - penborn*CF

            ampb(col(1),h2,h1,4,q1) = ampb(col(1),h2,h1,4,q1) + penborn
            ampv(col(1),h2,h1,4,q1) = ampv(col(1),h2,h1,4,q1) + penres(1)*CF + penres(2)*CFCA
            if( div.eq.0) ampv(col(1),h2,h1,4,q1) = ampv(col(1),h2,h1,4,q1) - penborn*CF
          endif
        enddo; enddo; enddo
        wardidtest(1)=wardpAbe(1)
        wardidtest_Q(1)=wardpAbe(2)   

c q > iV > Gij > iV2 > q
        comp=1
        ngauge=4
        do q1=3,4
        do h1=-1,1,2
        do h2=-1,1,2
          call penlineABETotalQ(pdia(0,qin),kdia(0,iV),momtemp(0),kdia(0,iV2),pdia(0,qout),
     &            psitemp(1,h1,qout),psitemp(1,h1,qin),
     &            epsV(0,h1,q1,iV,lephel,1),currij(0,h2,qin2,qout2),epsV(0,h1,q1,iV2,lephel,1),
     &            h1,musq,1,3,ngauge,comp,pengauge(1),penres(1),pengaugeb(1),penborn,div)
          comp =-1
          ngauge=9
          
          if(qline.eq.1) then
            ampb(col(1),h1,h2,q1,3) = ampb(col(1),h1,h2,q1,3) + penborn
            ampv(col(1),h1,h2,q1,3) = ampv(col(1),h1,h2,q1,3) + penres(1)*CF + penres(2)*CFCA
            if( div.eq.0) ampv(col(1),h1,h2,q1,3) = ampv(col(1),h1,h2,q1,3) - penborn*CF

            ampb(col(1),h1,h2,q1,4) = ampb(col(1),h1,h2,q1,4) + penborn
            ampv(col(1),h1,h2,q1,4) = ampv(col(1),h1,h2,q1,4) + penres(1)*CF + penres(2)*CFCA
            if( div.eq.0) ampv(col(1),h1,h2,q1,4) = ampv(col(1),h1,h2,q1,4) - penborn*CF
          else
            ampb(col(1),h2,h1,3,q1) = ampb(col(1),h2,h1,3,q1) + penborn
            ampv(col(1),h2,h1,3,q1) = ampv(col(1),h2,h1,3,q1) + penres(1)*CF + penres(2)*CFCA
            if( div.eq.0) ampv(col(1),h2,h1,3,q1) = ampv(col(1),h2,h1,3,q1) - penborn*CF

            ampb(col(1),h2,h1,4,q1) = ampb(col(1),h2,h1,4,q1) + penborn
            ampv(col(1),h2,h1,4,q1) = ampv(col(1),h2,h1,4,q1) + penres(1)*CF + penres(2)*CFCA
            if( div.eq.0) ampv(col(1),h2,h1,4,q1) = ampv(col(1),h2,h1,4,q1) - penborn*CF
          endif
        enddo; enddo; enddo
        wardidtest(2)=wardpAbe(1)
        wardidtest_Q(2)=wardpAbe(2)   

c q > iV > iV2 > Gij > q
        comp=1
        ngauge=4
        do q1=3,4
        do h1=-1,1,2
        do h2=-1,1,2
          call penlineABETotalQ(pdia(0,qin),kdia(0,iV),kdia(0,iV2),momtemp(0),pdia(0,qout),
     &            psitemp(1,h1,qout),psitemp(1,h1,qin),
     &            epsV(0,h1,q1,iV,lephel,1),epsV(0,h1,q1,iV2,lephel,1),currij(0,h2,qin2,qout2),
     &            h1,musq,1,4,ngauge,comp,pengauge(1),penres(1),pengaugeb(1),penborn,div)
          comp =-1
          ngauge=9
          
          if(qline.eq.1) then
            ampb(col(1),h1,h2,q1,3) = ampb(col(1),h1,h2,q1,3) + penborn
            ampv(col(1),h1,h2,q1,3) = ampv(col(1),h1,h2,q1,3) + penres(1)*CF + penres(2)*CFCA
            if( div.eq.0) ampv(col(1),h1,h2,q1,3) = ampv(col(1),h1,h2,q1,3) - penborn*CF

            ampb(col(1),h1,h2,q1,4) = ampb(col(1),h1,h2,q1,4) + penborn
            ampv(col(1),h1,h2,q1,4) = ampv(col(1),h1,h2,q1,4) + penres(1)*CF + penres(2)*CFCA
            if( div.eq.0) ampv(col(1),h1,h2,q1,4) = ampv(col(1),h1,h2,q1,4) - penborn*CF
          else
            ampb(col(1),h2,h1,3,q1) = ampb(col(1),h2,h1,3,q1) + penborn
            ampv(col(1),h2,h1,3,q1) = ampv(col(1),h2,h1,3,q1) + penres(1)*CF + penres(2)*CFCA
            if( div.eq.0) ampv(col(1),h2,h1,3,q1) = ampv(col(1),h2,h1,3,q1) - penborn*CF

            ampb(col(1),h2,h1,4,q1) = ampb(col(1),h2,h1,4,q1) + penborn
            ampv(col(1),h2,h1,4,q1) = ampv(col(1),h2,h1,4,q1) + penres(1)*CF + penres(2)*CFCA
            if( div.eq.0) ampv(col(1),h2,h1,4,q1) = ampv(col(1),h2,h1,4,q1) - penborn*CF
          endif
        enddo; enddo; enddo
        wardidtest(3)=wardpAbe(1)
        wardidtest_Q(3)=wardpAbe(2)   

c penline nonAbe
c q > iV > iV2 > q
        comp=1
        ngauge=4
        do q1=3,4
        do h1=-1,1,2
        do h2=-1,1,2
          call penlineNoABETotalQ(pdia(0,qin),kdia(0,iV),kdia(0,iV2),pdia(0,qout),momtemp(0),
     &            psitemp(1,h1,qout),psitemp(1,h1,qin),
     &            epsV(0,h1,q1,iV,lephel,1),epsV(0,h1,q1,iV2,lephel,1),currij(0,h2,qin2,qout2),
     &            h1,musq,ngauge,comp,pengauge(1),penres(1),div)
          comp =-1
          ngauge=9
          
          if(qline.eq.1) then
            ampv(col(1),h1,h2,q1,3) = ampv(col(1),h1,h2,q1,3) + penres(1)*CA
            ampv(col(1),h1,h2,q1,4) = ampv(col(1),h1,h2,q1,4) + penres(1)*CA
          else
            ampv(col(1),h2,h1,3,q1) = ampv(col(1),h2,h1,3,q1) + penres(1)*CA
            ampv(col(1),h2,h1,4,q1) = ampv(col(1),h2,h1,4,q1) + penres(1)*CA
          endif
        enddo; enddo; enddo
        wardidtest(4)=wardpNoAbe(1)
        wardidtest_Q(4)=wardpNoAbe(2)     
        GaugeOK(1) = GaugeOK(1).and.wardidtest(1).and.wardidtest(2).and.wardidtest(3).and.wardidtest(4)
     &                  .or.ldebug
        GaugeOK(2) = GaugeOK(2).and.wardidtest_q(1).and.wardidtest_q(2).and.wardidtest_q(3).and.wardidtest_q(4)
     &                  .or.ldebug
        enddo ! iV

      enddo ! qline


      pentcount2m = pentcount2m + 1
      pentcount2sub(sub_number) = pentcount2sub(sub_number) + 1
      if(.not.gaugeOK(1)) then
        pentcountm = pentcountm + 1
        pentcountsub(sub_number) = pentcountsub(sub_number) + 1
      endif
      if(.not.gaugeOK(2)) then
        GaugeOK(0) = .false.
        pentcountm_Q = pentcountm_q + 1
        pentcountsub_Q(sub_number) = pentcountsub_Q(sub_number) + 1
      endif

      endif ! penl

c ---------------------------------------
c $$$$$$$$$$ PenBox $$$$$$$$$$$$$$$$$$$$$
c ---------------------------------------

      if(lpenb.and.lZ4l) then

      if(ldebug) then
        print*
        print*, 'penbox Dire, in,out',in,out
      endif

      psi4(1,in) = psi(1,-1,id(in),ikin)
      psi4(2,in) = psi(2,-1,id(in),ikin)
      psi4(3,in) = psi(1, 1,id(in),ikin)
      psi4(4,in) = psi(2, 1,id(in),ikin)

      psi4(3,out) = psi(1,-1,-id(out),ikin)
      psi4(4,out) = psi(2,-1,-id(out),ikin)
      psi4(1,out) = psi(1, 1,-id(out),ikin)
      psi4(2,out) = psi(2, 1,-id(out),ikin)

      psi4(1,in2) = psi(1,-1,id(in2),ikin)
      psi4(2,in2) = psi(2,-1,id(in2),ikin)
      psi4(3,in2) = psi(1, 1,id(in2),ikin)
      psi4(4,in2) = psi(2, 1,id(in2),ikin)

      psi4(3,out2) = psi(1,-1,-id(out2),ikin)
      psi4(4,out2) = psi(2,-1,-id(out2),ikin)
      psi4(1,out2) = psi(1, 1,-id(out2),ikin)
      psi4(2,out2) = psi(2, 1,-id(out2),ikin)

      gaugeOK(1) = .true.
c in2  -->---$$$$$---<-- in
c           |     |__Z
c           |     |
c out2 --<---$$$$$--->-- out

c start with gauge test
c use hel=77
      call PenBoxT(pdia(0,out2),pdia(0,out),kdia(0,3),pdia(0,in),pdia(0,in2),
     &         psi4(1,out2),psi4(1,in2),psi4(1,out),psi4(1,in),kdiac(0,3),
     &         musq,1,0,0,penbamp,penbres,penbborn,div,77)
     
      if (abs(penbamp(1)).gt.abs(penbamp(2))) then
        td = abs((penbamp(2)+penbamp(3))/penbamp(1)+1d0)
      else
        td = abs((penbamp(1)+penbamp(3))/penbamp(2)+1d0)
      endif
      if(ldebug) print*,'penbox gauge 1',td
      GaugeOK(1) = GaugeOK(1).and.(td.le.Gaugeaccuracy.or.ldebug)
      if(ldebug) gaugeprecision = max(gaugeprecision,td)
      
      if(GaugeOK(1)) then
        do h1=-1,1,2
        do h2=-1,1,2
        do q=3,4
          call PenBoxT(pdia(0,out2),pdia(0,out),kdia(0,3),pdia(0,in),pdia(0,in2),
     &               psi4(1,out2),psi4(1,in2),psi4(1,out),psi4(1,in),epsV(0,h1,q,3,lephel,lepkin),
     &               musq,-1,0,0,penbamp,penbres,penbborn,div,2*h2+h1)
          ampb(col(1),h1,h2,q,3) = ampb(col(1),h1,h2,q,3) - penbborn
          ampv(col(1),h1,h2,q,3) = ampv(col(1),h1,h2,q,3) + penbres/6d0
          ampv(col(2),h1,h2,q,3) = ampv(col(2),h1,h2,q,3) - penbres/2d0

          ampb(col(1),h1,h2,q,4) = ampb(col(1),h1,h2,q,4) - penbborn
          ampv(col(1),h1,h2,q,4) = ampv(col(1),h1,h2,q,4) + penbres/6d0
          ampv(col(2),h1,h2,q,4) = ampv(col(2),h1,h2,q,4) - penbres/2d0
        enddo;enddo;enddo
        
      endif

c in  -->---$$$$$---<-- in2
c           |     |__Z
c           |     |
c out --<---$$$$$--->-- out2

c start with gauge test
c use hel=77
      call PenBoxT(pdia(0,out),pdia(0,out2),kdia(0,3),pdia(0,in2),pdia(0,in),
     &         psi4(1,out),psi4(1,in),psi4(1,out2),psi4(1,in2),kdiac(0,3),
     &         musq,1,0,0,penbamp,penbres,penbborn,div,77)
     
      if (abs(penbamp(1)).gt.abs(penbamp(2))) then
        td = abs((penbamp(2)+penbamp(3))/penbamp(1)+1d0)
      else
        td = abs((penbamp(1)+penbamp(3))/penbamp(2)+1d0)
      endif
      if(ldebug) print*,'penbox gauge 2',td
      GaugeOK(1) = GaugeOK(1).and.(td.le.Gaugeaccuracy.or.ldebug)
      if(ldebug) gaugeprecision = max(gaugeprecision,td)
      
      if(GaugeOK(1)) then
        do h1=-1,1,2
        do h2=-1,1,2
        do q=3,4
          call PenBoxT(pdia(0,out),pdia(0,out2),kdia(0,3),pdia(0,in2),pdia(0,in),
     &               psi4(1,out),psi4(1,in),psi4(1,out2),psi4(1,in2),epsV(0,h2,q,3,lephel,lepkin),
     &               musq,-1,0,0,penbamp,penbres,penbborn,div,2*h1+h2)
          ampb(col(1),h1,h2,3,q) = ampb(col(1),h1,h2,3,q) - penbborn
          ampv(col(1),h1,h2,3,q) = ampv(col(1),h1,h2,3,q) + penbres/6d0
          ampv(col(2),h1,h2,3,q) = ampv(col(2),h1,h2,3,q) - penbres/2d0

          ampb(col(1),h1,h2,4,q) = ampb(col(1),h1,h2,4,q) - penbborn
          ampv(col(1),h1,h2,4,q) = ampv(col(1),h1,h2,4,q) + penbres/6d0
          ampv(col(2),h1,h2,4,q) = ampv(col(2),h1,h2,4,q) - penbres/2d0
        enddo;enddo;enddo
        
      endif

      if(ldebug) then
        print*
        print*, 'penbox Cross, in,out',in,out
      endif

c out2 --<---$$$$$---<-- in
c           |     |__Z
c           |     |
c in2  -->---$$$$$--->-- out

c start with gauge test
c use hel=77
      call PenBoxCrossT(pdia(0,in2),pdia(0,out),kdia(0,3),pdia(0,in),pdia(0,out2),
     &         psi4(1,in2),psi4(1,out2),psi4(1,out),psi4(1,in),kdiac(0,3),
     &         musq,1,0,0,penbamp,penbres,penbborn,div,77)
     
      if (abs(penbamp(1)).gt.abs(penbamp(2))) then
        td = abs((penbamp(2)+penbamp(3))/penbamp(1)+1d0)
      else
        td = abs((penbamp(1)+penbamp(3))/penbamp(2)+1d0)
      endif
      if(ldebug) print*,'penbox gauge 1',td
      GaugeOK(1) = GaugeOK(1).and.(td.le.Gaugeaccuracy.or.ldebug)
      if(ldebug) gaugeprecision = max(gaugeprecision,td)
      
      if(GaugeOK(1)) then
        do h1=-1,1,2
        do h2=-1,1,2
        do q=3,4
          call PenBoxCrossT(pdia(0,in2),pdia(0,out),kdia(0,3),pdia(0,in),pdia(0,out2),
     &            psi4(1,in2),psi4(1,out2),psi4(1,out),psi4(1,in),epsV(0,h1,q,3,lephel,lepkin),
     &               musq,-1,0,0,penbamp,penbres,penbborn,div,2*h2+h1)
          ampb(col(1),h1,h2,q,3) = ampb(col(1),h1,h2,q,3) - penbborn
          ampv(col(1),h1,h2,q,3) = ampv(col(1),h1,h2,q,3) - penbres*4d0/3d0
          ampv(col(2),h1,h2,q,3) = ampv(col(2),h1,h2,q,3) - penbres/2d0

          ampb(col(1),h1,h2,q,4) = ampb(col(1),h1,h2,q,4) - penbborn
          ampv(col(1),h1,h2,q,4) = ampv(col(1),h1,h2,q,4) - penbres*4d0/3d0
          ampv(col(2),h1,h2,q,4) = ampv(col(2),h1,h2,q,4) - penbres/2d0
        enddo;enddo;enddo
        
      endif

c out --<---$$$$$---<-- in2
c           |     |__Z
c           |     |
c in  -->---$$$$$--->-- out2

c start with gauge test
c use hel=77
      call PenBoxCrossT(pdia(0,in),pdia(0,out2),kdia(0,3),pdia(0,in2),pdia(0,out),
     &         psi4(1,in),psi4(1,out),psi4(1,out2),psi4(1,in2),kdiac(0,3),
     &         musq,1,0,0,penbamp,penbres,penbborn,div,77)
     
      if (abs(penbamp(1)).gt.abs(penbamp(2))) then
        td = abs((penbamp(2)+penbamp(3))/penbamp(1)+1d0)
      else
        td = abs((penbamp(1)+penbamp(3))/penbamp(2)+1d0)
      endif
      if(ldebug) print*,'penbox gauge 2',td
      GaugeOK(1) = GaugeOK(1).and.(td.le.Gaugeaccuracy.or.ldebug)
      if(ldebug) gaugeprecision = max(gaugeprecision,td)
      
      if(GaugeOK(1)) then
        do h1=-1,1,2
        do h2=-1,1,2
        do q=3,4
          call PenBoxCrossT(pdia(0,in),pdia(0,out2),kdia(0,3),pdia(0,in2),pdia(0,out),
     &         psi4(1,in),psi4(1,out),psi4(1,out2),psi4(1,in2),epsV(0,h2,q,3,lephel,lepkin),
     &               musq,-1,0,0,penbamp,penbres,penbborn,div,2*h1+h2)
          ampb(col(1),h1,h2,3,q) = ampb(col(1),h1,h2,3,q) - penbborn
          ampv(col(1),h1,h2,3,q) = ampv(col(1),h1,h2,3,q) - penbres*4d0/3d0
          ampv(col(2),h1,h2,3,q) = ampv(col(2),h1,h2,3,q) - penbres/2d0

          ampb(col(1),h1,h2,4,q) = ampb(col(1),h1,h2,4,q) - penbborn
          ampv(col(1),h1,h2,4,q) = ampv(col(1),h1,h2,4,q) - penbres*4d0/3d0
          ampv(col(2),h1,h2,4,q) = ampv(col(2),h1,h2,4,q) - penbres/2d0
        enddo;enddo;enddo
        
      endif

      box2count2m = box2count2m + 1
      box2count2sub(sub_number) = box2count2sub(sub_number) + 1
      if(.not.gaugeOK(1)) then
        GaugeOK(0) = .false.
        box2countm   = box2countm   + 1
        box2countm_Q = box2countm_Q + 1
        box2countsub  (sub_number) = box2countsub  (sub_number) + 1
        box2countsub_q(sub_number) = box2countsub_q(sub_number) + 1
      endif
      endif ! penb
      
c ---------------------------------------
c $$$$$$$$$$ HexBox $$$$$$$$$$$$$$$$$$$$$
c ---------------------------------------

      if(lhex.and.lZZ) then

      if(ldebug) then
        print*
        print*, 'HexBox Dire, in,out',in,out
      endif

      psi4(1,in) = psi(1,-1,id(in),ikin)
      psi4(2,in) = psi(2,-1,id(in),ikin)
      psi4(3,in) = psi(1, 1,id(in),ikin)
      psi4(4,in) = psi(2, 1,id(in),ikin)

      psi4(3,out) = psi(1,-1,-id(out),ikin)
      psi4(4,out) = psi(2,-1,-id(out),ikin)
      psi4(1,out) = psi(1, 1,-id(out),ikin)
      psi4(2,out) = psi(2, 1,-id(out),ikin)

      psi4(1,in2) = psi(1,-1,id(in2),ikin)
      psi4(2,in2) = psi(2,-1,id(in2),ikin)
      psi4(3,in2) = psi(1, 1,id(in2),ikin)
      psi4(4,in2) = psi(2, 1,id(in2),ikin)

      psi4(3,out2) = psi(1,-1,-id(out2),ikin)
      psi4(4,out2) = psi(2,-1,-id(out2),ikin)
      psi4(1,out2) = psi(1, 1,-id(out2),ikin)
      psi4(2,out2) = psi(2, 1,-id(out2),ikin)

      gaugeOK(1) = .true.

      do iV=1,2
      iV2=3-iV
c in2  -->---$$$$$---<-- in
c           |     |__Z
c           |     |__Z
c out2 --<---$$$$$--->-- out

c start with gauge test
c use hel=77
      call HexBoxT(pdia(0,out2),pdia(0,out),kdia(0,iV),kdia(0,iV2),pdia(0,in),pdia(0,in2),
     &         psi4(1,out2),psi4(1,in2),psi4(1,out),psi4(1,in),kdiac(0,iV),epsV(0,-1,3,iV2,lephel,lepkin),
     &         musq,1,0,0,hexamp,hexresG(iV,1),hexborn,div,77)

      call HexBoxT(pdia(0,out2),pdia(0,out),kdia(0,iV),kdia(0,iV2),pdia(0,in),pdia(0,in2),
     &         psi4(1,out2),psi4(1,in2),psi4(1,out),psi4(1,in),epsV(0,-1,3,iV,lephel,lepkin),kdiac(0,iV2),
     &         musq,-1,0,0,hexamp,hexresG(iV2,2),hexborn,div,77)
     
      do h1=-1,1,2
      do h2=-1,1,2
      do q=3,4
        call HexBoxT(pdia(0,out2),pdia(0,out),kdia(0,iV),kdia(0,iV2),pdia(0,in),pdia(0,in2),
     &          psi4(1,out2),psi4(1,in2),psi4(1,out),psi4(1,in),epsV(0,h1,q,iV,lephel,lepkin),epsV(0,h1,q,iV2,lephel,lepkin),
     &          musq,-1,0,0,hexamp,hexres,hexborn,div,2*h2+h1)
          ampb(col(1),h1,h2,q,3) = ampb(col(1),h1,h2,q,3) + hexborn
          ampv(col(1),h1,h2,q,3) = ampv(col(1),h1,h2,q,3) - hexres/6d0
          ampv(col(2),h1,h2,q,3) = ampv(col(2),h1,h2,q,3) + hexres/2d0

          ampb(col(1),h1,h2,q,4) = ampb(col(1),h1,h2,q,4) + hexborn
          ampv(col(1),h1,h2,q,4) = ampv(col(1),h1,h2,q,4) - hexres/6d0
          ampv(col(2),h1,h2,q,4) = ampv(col(2),h1,h2,q,4) + hexres/2d0
      enddo;enddo;enddo
      enddo ! iV

      if(ldebug) then
        print*,'HexBox 1 gauge', (abs(sum(hexresg(i,:))/hexresg(i,1)),i=1,2)
      endif

      GaugeOK(1) = GaugeOK(1).and.abs(sum(hexresg(1,:))/hexresg(1,1)) .le. GaugeAccuracy
     $                       .and.abs(sum(hexresg(2,:))/hexresg(2,1)) .le. GaugeAccuracy

      do iV=1,2
      iV2=3-iV
c in  -->---$$$$$---<-- in2
c           |     |__Z
c           |     |__Z
c out --<---$$$$$--->-- out2

c start with gauge test
c use hel=77
      call HexBoxT(pdia(0,out),pdia(0,out2),kdia(0,iV),kdia(0,iV2),pdia(0,in2),pdia(0,in),
     &         psi4(1,out),psi4(1,in),psi4(1,out2),psi4(1,in2),kdiac(0,iV),epsV(0,-1,3,iV2,lephel,lepkin),
     &         musq,1,0,0,hexamp,hexresG(iV,1),hexborn,div,77)

      call HexBoxT(pdia(0,out),pdia(0,out2),kdia(0,iV),kdia(0,iV2),pdia(0,in2),pdia(0,in),
     &         psi4(1,out),psi4(1,in),psi4(1,out2),psi4(1,in2),epsV(0,-1,3,iV,lephel,lepkin),kdiac(0,iV2),
     &         musq,-1,0,0,hexamp,hexresG(iV2,2),hexborn,div,77)
     
      do h1=-1,1,2
      do h2=-1,1,2
      do q=3,4
        call HexBoxT(pdia(0,out),pdia(0,out2),kdia(0,iV),kdia(0,iV2),pdia(0,in2),pdia(0,in),
     &          psi4(1,out),psi4(1,in),psi4(1,out2),psi4(1,in2),epsV(0,h2,q,iV,lephel,lepkin),epsV(0,h2,q,iV2,lephel,lepkin),
     &          musq,-1,0,0,hexamp,hexres,hexborn,div,2*h1+h2)
          ampb(col(1),h1,h2,3,q) = ampb(col(1),h1,h2,3,q) + hexborn
          ampv(col(1),h1,h2,3,q) = ampv(col(1),h1,h2,3,q) - hexres/6d0
          ampv(col(2),h1,h2,3,q) = ampv(col(2),h1,h2,3,q) + hexres/2d0
                                                      
          ampb(col(1),h1,h2,4,q) = ampb(col(1),h1,h2,4,q) + hexborn
          ampv(col(1),h1,h2,4,q) = ampv(col(1),h1,h2,4,q) - hexres/6d0
          ampv(col(2),h1,h2,4,q) = ampv(col(2),h1,h2,4,q) + hexres/2d0
      enddo;enddo;enddo
      enddo ! iV

      if(ldebug) then
        print*,'HexBox 1 gauge', (abs(sum(hexresg(i,:))/hexresg(i,1)),i=1,2)
      endif

      GaugeOK(1) = GaugeOK(1).and.abs(sum(hexresg(1,:))/hexresg(1,1)) .le. GaugeAccuracy
     $                       .and.abs(sum(hexresg(2,:))/hexresg(2,1)) .le. GaugeAccuracy

      if(ldebug) then
        print*, 'HexBox Cross 1'
      endif

      do iV=1,2
      iV2=3-iV
c out2 --<---$$$$$---<-- in
c           |     |__Z
c           |     |__Z
c in2  -->---$$$$$--->-- out

c start with gauge test
c use hel=77
      call HexBoxCross1T(pdia(0,in2),pdia(0,out),kdia(0,iV),kdia(0,iV2),pdia(0,in),pdia(0,out2),
     &         psi4(1,in2),psi4(1,out2),psi4(1,out),psi4(1,in),kdiac(0,iV),epsV(0,-1,3,iV2,lephel,lepkin),
     &         musq,1,0,0,hexamp,hexresG(iV,1),hexborn,div,77)

      call HexBoxCross1T(pdia(0,in2),pdia(0,out),kdia(0,iV),kdia(0,iV2),pdia(0,in),pdia(0,out2),
     &         psi4(1,in2),psi4(1,out2),psi4(1,out),psi4(1,in),epsV(0,-1,3,iV,lephel,lepkin),kdiac(0,iV2),
     &         musq,-1,0,0,hexamp,hexresG(iV2,2),hexborn,div,77)
     
      do h1=-1,1,2
      do h2=-1,1,2
      do q=3,4
        call HexBoxCross1T(pdia(0,in2),pdia(0,out),kdia(0,iV),kdia(0,iV2),pdia(0,in),pdia(0,out2),
     &          psi4(1,in2),psi4(1,out2),psi4(1,out),psi4(1,in),epsV(0,h1,q,iV,lephel,lepkin),epsV(0,h1,q,iV2,lephel,lepkin),
     &          musq,-1,0,0,hexamp,hexres,hexborn,div,2*h2+h1)
          ampb(col(1),h1,h2,q,3) = ampb(col(1),h1,h2,q,3) + hexborn
          ampv(col(1),h1,h2,q,3) = ampv(col(1),h1,h2,q,3) + hexres*CF
          ampv(col(2),h1,h2,q,3) = ampv(col(2),h1,h2,q,3) + hexres/2d0

          ampb(col(1),h1,h2,q,4) = ampb(col(1),h1,h2,q,4) + hexborn
          ampv(col(1),h1,h2,q,4) = ampv(col(1),h1,h2,q,4) + hexres*CF
          ampv(col(2),h1,h2,q,4) = ampv(col(2),h1,h2,q,4) + hexres/2d0
      enddo;enddo;enddo
      enddo ! iV

      if(ldebug) then
        print*,'HexBox Cross 1 gauge', (abs(sum(hexresg(i,:))/hexresg(i,1)),i=1,2)
      endif

      GaugeOK(1) = GaugeOK(1).and.abs(sum(hexresg(1,:))/hexresg(1,1)) .le. GaugeAccuracy
     $                       .and.abs(sum(hexresg(2,:))/hexresg(2,1)) .le. GaugeAccuracy

      do iV=1,2
      iV2=3-iV
c out -->---$$$$$---<-- in2
c           |     |__Z
c           |     |__Z
c in  --<---$$$$$--->-- out2

c start with gauge test
c use hel=77
      call HexBoxCross1T(pdia(0,in),pdia(0,out2),kdia(0,iV),kdia(0,iV2),pdia(0,in2),pdia(0,out),
     &         psi4(1,in),psi4(1,out),psi4(1,out2),psi4(1,in2),kdiac(0,iV),epsV(0,-1,3,iV2,lephel,lepkin),
     &         musq,1,0,0,hexamp,hexresG(iV,1),hexborn,div,77)

      call HexBoxCross1T(pdia(0,in),pdia(0,out2),kdia(0,iV),kdia(0,iV2),pdia(0,in2),pdia(0,out),
     &         psi4(1,in),psi4(1,out),psi4(1,out2),psi4(1,in2),epsV(0,-1,3,iV,lephel,lepkin),kdiac(0,iV2),
     &         musq,-1,0,0,hexamp,hexresG(iV2,2),hexborn,div,77)
     
      do h1=-1,1,2
      do h2=-1,1,2
      do q=3,4
        call HexBoxCross1T(pdia(0,in),pdia(0,out2),kdia(0,iV),kdia(0,iV2),pdia(0,in2),pdia(0,out),
     &          psi4(1,in),psi4(1,out),psi4(1,out2),psi4(1,in2),epsV(0,h2,q,iV,lephel,lepkin),epsV(0,h2,q,iV2,lephel,lepkin),
     &          musq,-1,0,0,hexamp,hexres,hexborn,div,2*h1+h2)
          ampb(col(1),h1,h2,3,q) = ampb(col(1),h1,h2,3,q) + hexborn
          ampv(col(1),h1,h2,3,q) = ampv(col(1),h1,h2,3,q) + hexres*CF
          ampv(col(2),h1,h2,3,q) = ampv(col(2),h1,h2,3,q) + hexres/2d0
                                                      
          ampb(col(1),h1,h2,4,q) = ampb(col(1),h1,h2,4,q) + hexborn
          ampv(col(1),h1,h2,4,q) = ampv(col(1),h1,h2,4,q) + hexres*CF
          ampv(col(2),h1,h2,4,q) = ampv(col(2),h1,h2,4,q) + hexres/2d0
      enddo;enddo;enddo
      enddo ! iV

      if(ldebug) then
        print*,'HexBox Cross 1 gauge', (abs(sum(hexresg(i,:))/hexresg(i,1)),i=1,2)
      endif

      GaugeOK(1) = GaugeOK(1).and.abs(sum(hexresg(1,:))/hexresg(1,1)) .le. GaugeAccuracy
     $                       .and.abs(sum(hexresg(2,:))/hexresg(2,1)) .le. GaugeAccuracy


      if(ldebug) print*,'HexBox2'
c in  -->---$$$$$---<-- in2
c        Z__|     |__Z
c           |     |
c out --<---$$$$$--->-- out2
      do iV=1,2
      iV2=3-iV

c start with gauge test
c use hel=77
      call HexBox2T(pdia(0,out),pdia(0,out2),kdia(0,iV),pdia(0,in2),pdia(0,in),kdia(0,iV2),
     &         psi4(1,out),psi4(1,in),psi4(1,out2),psi4(1,in2),kdiac(0,iV),epsV(0,-1,3,iV2,lephel,lepkin),
     &         musq,1,0,0,hexamp,hexres,hexborn,div,77)
      if(ldebug) then
        print*, "Gauge test HexBox 2"
        call check('r127',(hexamp(1)+hexamp(2))/(-hexamp(7)))
        call check('r348',(hexamp(3)+hexamp(4))/(-hexamp(8)))
        call check('r569',(hexamp(5)+hexamp(6))/(-hexamp(9)))
      endif
      GaugeOK(1) =GaugeOK(1).and. abs((hexamp(1)+hexamp(2))/(-hexamp(7))-1d0) .lt. GaugeAccuracy
     &           .and.abs((hexamp(3)+hexamp(4))/(-hexamp(8))-1d0) .lt. GaugeAccuracy
     &           .and.abs((hexamp(5)+hexamp(6))/(-hexamp(9))-1d0) .lt. GaugeAccuracy
     
      do h1=-1,1,2
      do h2=-1,1,2
      do q1=3,4
      do q2=3,4
        call HexBox2T(pdia(0,out),pdia(0,out2),kdia(0,iV),pdia(0,in2),pdia(0,in),kdia(0,iV2),
     &         psi4(1,out),psi4(1,in),psi4(1,out2),psi4(1,in2),epsV(0,h2,q2,iV,lephel,lepkin),epsV(0,h1,q1,iV2,lephel,lepkin),
     &         musq,-1,0,0,hexamp,hexres,hexborn,div,2*h1+h2)
          ampb(col(1),h1,h2,q1,q2) = ampb(col(1),h1,h2,q1,q2) + hexborn
          ampv(col(1),h1,h2,q1,q2) = ampv(col(1),h1,h2,q1,q2) - hexres/6d0
          ampv(col(2),h1,h2,q1,q2) = ampv(col(2),h1,h2,q1,q2) + hexres/2d0
      enddo;enddo;enddo; enddo
      enddo ! iV

c out --<---$$$$$---<-- in2
c        Z__|     |__Z
c           |     |
c in  -->---$$$$$--->-- out2
      do iV=1,2
      iV2=3-iV

c start with gauge test
c use hel=77
      call HexBoxCross2T(pdia(0,in),pdia(0,out2),kdia(0,iV),pdia(0,in2),pdia(0,out),kdia(0,iV2),
     &         psi4(1,in),psi4(1,out),psi4(1,out2),psi4(1,in2),kdiac(0,iV),epsV(0,-1,3,iV2,lephel,lepkin),
     &         musq,1,0,0,hexamp,hexres,hexborn,div,77)
      if(ldebug) then
        print*, "Gauge test HexBox Cross 2"
        call check('r127',(hexamp(1)+hexamp(2))/(-hexamp(7)))
        call check('r348',(hexamp(3)+hexamp(4))/(-hexamp(8)))
        call check('r569',(hexamp(5)+hexamp(6))/(-hexamp(9)))
      endif
      GaugeOK(1) =GaugeOK(1).and. abs((hexamp(1)+hexamp(2))/(-hexamp(7))-1d0) .lt. GaugeAccuracy
     &           .and.abs((hexamp(3)+hexamp(4))/(-hexamp(8))-1d0) .lt. GaugeAccuracy
     &           .and.abs((hexamp(5)+hexamp(6))/(-hexamp(9))-1d0) .lt. GaugeAccuracy
     
      do h1=-1,1,2
      do h2=-1,1,2
      do q1=3,4
      do q2=3,4
        call HexBoxCross2T(pdia(0,in),pdia(0,out2),kdia(0,iV),pdia(0,in2),pdia(0,out),kdia(0,iV2),
     &         psi4(1,in),psi4(1,out),psi4(1,out2),psi4(1,in2),epsV(0,h2,q2,iV,lephel,lepkin),epsV(0,h1,q1,iV2,lephel,lepkin),
     &         musq,-1,0,0,hexamp,hexres,hexborn,div,2*h1+h2)
          ampb(col(1),h1,h2,q1,q2) = ampb(col(1),h1,h2,q1,q2) + hexborn
          ampv(col(1),h1,h2,q1,q2) = ampv(col(1),h1,h2,q1,q2) + hexres*CF
          ampv(col(2),h1,h2,q1,q2) = ampv(col(2),h1,h2,q1,q2) + hexres/2d0
      enddo;enddo;enddo; enddo
      enddo ! iV

      hexcount2m = hexcount2m + 1
      hexcount2sub(sub_number) = hexcount2sub(sub_number) + 1
      if(.not.gaugeOK(1)) then
        GaugeOK(0) = .false.
        hexcountm = hexcountm + 1
        hexcountsub(sub_number) = hexcountsub(sub_number) + 1
        hexcountm_Q = hexcountm_q + 1
        hexcountsub_Q(sub_number) = hexcountsub_Q(sub_number) + 1
      endif

      endif ! hexbox


c ---------------------------------------
c $$$$$$$$$$ Sum up everything $$$$$$$$$$
c ---------------------------------------

      
      if(ldebug) then
        h1=1
        h2=-1
        print*,'showing results for q,hel1,2:',h1,h2
      do q1=3,4
      do q2=3,4
        print*,'same born amp?'
        print*,ampb(col(1),h1,h2,q1,q2)/ampborn(col(1),h1,h2,q1,q2)
      enddo
      enddo
      print*,ampb(col(1),1:1:2,-1:1:2,:,:)/ampborn(col(1),1:1:2,-1:1:2,:,:)
      endif

      do q1=3,4
      do q2=3,4
      do h1=-1,1,2
      do h2=-1,1,2
        squared(1,1) = 2d0*dble(ampv(1,h1,h2,q1,q2)*dconjg(ampborn(1,h1,h2,q1,q2)))
        squared(2,1) = 2d0*dble(ampv(2,h1,h2,q1,q2)*dconjg(ampborn(1,h1,h2,q1,q2)))
        squared(2,2) = 2d0*dble(ampv(2,h1,h2,q1,q2)*dconjg(ampborn(2,h1,h2,q1,q2)))
        squared(1,2) = 2d0*dble(ampv(1,h1,h2,q1,q2)*dconjg(ampborn(2,h1,h2,q1,q2)))

c different flavor result (for q1.eq.q2 the generations can be different)
        resv(q1,q2,col(1)) = resv(q1,q2,col(1)) + (2d0*squared(col(1),col(1))-2d0/3d0*squared(col(2),col(1)))
     $                          *G**4*als(1,1)/(4d0* Pi) ! alphas/4pi
c same flavor
        if(q1.eq.q2) then
          if(h1.ne.h2) then
            resv(q1,q2,0) = resv(q1,q2,0) + (2d0*squared(col(1),col(1))-2d0/3d0*squared(col(2),col(1)))
     $                          *G**4*als(1,1)/(4d0* Pi) ! alphas/4pi
          else
            if(in.eq.1) then
              ampvCross1=ampv ! store
            else ! antisymmetrize
              ampv(:,h1,h2,q1,q2) =  ampvCross1(:,h1,h2,q1,q2) - ampv(:,h1,h2,q1,q2) 
              ampborn(2,h1,h2,q1,q2) = -ampborn(2,h1,h2,q1,q2)
              squared(1,1) = 2d0*dble(ampv(1,h1,h2,q1,q2)*dconjg(ampborn(1,h1,h2,q1,q2)))
              squared(2,1) = 2d0*dble(ampv(2,h1,h2,q1,q2)*dconjg(ampborn(1,h1,h2,q1,q2)))
              squared(2,2) = 2d0*dble(ampv(2,h1,h2,q1,q2)*dconjg(ampborn(2,h1,h2,q1,q2)))
              squared(1,2) = 2d0*dble(ampv(1,h1,h2,q1,q2)*dconjg(ampborn(2,h1,h2,q1,q2)))
            resv(q1,q2,0) = resv(q1,q2,0) + (2d0*(squared(1,1)+squared(2,2))-2d0/3d0*(squared(2,1)+squared(1,2)))
     $                            *G**4*als(1,1)/(4d0* Pi) ! alphas/4pi
            endif
          endif
        endif
          
      enddo; enddo ! helicities
      enddo; enddo

      enddo ! in


      if (ldebug) then
        print*,'GaugeOK',GaugeOK(0)
        do q1=3,4
        do q2=3,4
          print*,'q1,q2=' , q1,q2

          if(div.ne.0) then
            print*,'results for same flavor and both crossings'
            print*,'resB',resB(q1,q2,:)
            print*,'resV',resV(q1,q2,:)
            print*,'resI',resI(q1,q2,:)
            print*,'canc',(resV(q1,q2,:)+resI(q1,q2,:))/resI(q1,q2,:)
          endif
        enddo
        enddo

      endif

      if(GaugeOK(0)) then ! all gauge test ok
        res(0,:,:,:)=resV+resI
        if(lborn) res(0,:,:,:)=res(0,:,:,:) + resB
        res(1,:,:,1)=res(0,:,:,1)
        res(2,:,:,2)=res(0,:,:,2)
      else
        res=0d0
      endif

      end       
