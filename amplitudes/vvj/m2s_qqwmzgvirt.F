! This subroutine returns the virtual correction  p p --> e- ve mu+ mu- j

      subroutine m2s_ppwzjvirt(xi,rn,pdfchange,p,v,sumresult,nlo)

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/scales.inc"

      integer k,mu,photpolin,photin,glupolin,gluin,signspinor(2),map(2),nlo,
     &   incross,mapspin1(2),mapspin2(2)

      integer badpointspent(6),badpointspentnoab(6),badpointsbox(6),badpoints

     
      common /badpointsWMZG/ badpointspent,badpointsbox,badpoints

      common /badpointsnoabWMZG/ badpointspentnoab


      double precision photmom(0:4),glumom(0:4),wmom(0:3),pentmom(0:3,5),
     &  wdecmom(0:3,2),spinormom(0:3,2),mur,m2s(6),
     &  boxmom(0:3,4),CF,CA,NF,TR,m2sborn(6),m2sbornwopdf(6),m2stest(6),sumresult,boxmomtilde(0:3,4),
     &  sumresultalt,finiteubd,s,t,u,mursq1,rn,
     &  finiteubg,finiteubgcr,finiteqg,s15,eulergamma,inv2,inv3,musq,tf,finiteubdf,bornelement,
     &  mt2,effglumom(0:3),quarkmom(0:3),triangleresult(2,3),quarkmom2(0:3),trianglemom(0:3,3),fermres
      double complex photpol(0:3,2), glupol(0:3), wpol(0:3),psi(2,-1:1,2),
     &  resultgauge(3),result(6),resultgaugeb(3),temp(6),
     &  resultb(6),virt,ampborn,cte,penresult(8),cfc,cac,cfcac,boxresult(5),test(11),testb(11),
     &  ancont,storephotpol(0:3,-1:1,2),storewpol(0:3),storewphotpol(0:3,-1:1),testgauge1,
     &  testgauge2,pic2,im,cacgraph,fact,factfigy,wavea(6,-1:1),
     &  wavez(6,-1:1),effglucurr(0:3),ampbornubg(2,-1:1,-1:1),ampborndg(2,-1:1,-1:1),
     &  ampbornubd(2,-1:1,-1:1),cptri,resultgauge2(3),resultgauge3(3)

      double complex AMPUBDGEMVEBMUMMU,AMPUBGDBEMVEBMUMM,AMPDGUBEMVEBMUMMU
      external AMPUBDGEMVEBMUMMU,AMPUBGDBEMVEBMUMM,AMPDGUBEMVEBMUMMU
      integer NHEL(7)

      double precision accuracy,accuracyres,zero

      parameter (accuracy=0.1d0,zero=0d0)

      logical wardresult,ldebug,borndebug,fincolcalc,pdfchange,calctriangle,ubd(2),ubg(2),dg(2),badgauge
      parameter (ldebug=.true.,borndebug=.false.,fincolcalc=.false.,calctriangle=.true.)
     
      parameter (NF=5d0,inv2=1d0/2d0,inv3=1d0/3d0)
      parameter (CA=3d0)
      parameter (TR=1d0/2d0,tf=1d0/2d0)
      parameter (CF=4d0/3d0)
      parameter (eulergamma=5.772156649015328D-1)

      double precision p(0:3,max_p,max_kin), v(0:3,max_v,max_kin),xi(nx),pdf(-6:6,2)
     &   

      double precision helastestmom(0:3,7),testh,m2s_qqwzglu_hel

      double precision dotrr,theta
      external dotrr,theta,m2s_qqwzglu_hel

      double precision wemnugamom(0:3)
      double complex currlepwga(0:3,-1:1),currlepwga2(0:3,-1:1),currlepz(0:3,-1:1),currlepa(0:3,-1:1),
     &  wem(6),wve(6),wmcurr(6),wmum(6),wmup(6),bornubd(2,-1:1,-1:1),bornubg(2,-1:1,-1:1),borndg(2,-1:1,-1:1)

      double precision tmp(0:3)
      double complex xp1,xp2,xw,xg,dotcr,wpoltilde(0:3),photpoltilde(0:3,2),glupoltilde(0:3),
     &   testpb1,testpb2
      external dotcr
!      logical nlopent, nlobox

** Used for gauge test of non-abelian boxes
      double complex nothing(5),gaugeCA1(5),gaugeCA2(5),polgluontemp(0:3)

      INTEGER init/0/
      SAVE init


      if (.not.pdfchange) then

      if ( init .eq. 0 ) then
         write(6,*) ' '
         write(6,*) 'W-Zj virtual amplitude square information:'
         write(6,*) '-----------------------------------------------'
         write(6,*) ' '
         write(6,*) ' '
         write(6,*) '-----------------------------------------------'
      endif
      endif
      init=1

      if (nlo.eq.1) then
         sumresult=m2s_qqwzglu_hel(xi,p(0,1,1),v(0,1,1),nlo,pdfchange)
      else

      call getyourscalesready(xi,p(0,1,1),v(0,1,1),pdfchange,pdf,nlo,1)!.true.)

      do k=1,2
      ubd(k)=.true.
      ubg(k)=.true.
      dg(k)= .true.
      enddo

!       call calcWgajScales(p(0,1,1),v(0,1,1),1)
      mur=mursq(1,1)
      musq=mur

! running coupling
!       alfas  =  als(1,1) 
!       G      =  dsqrt(4d0*pi*alfas)
!       GG(1)  =  -G
!       GG(2)  =  -G


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do incross=1,2
      do glupolin=-1,1,2
      do photpolin=-1,1,2
      bornubd(incross,glupolin,photpolin)=0d0
      NHEL(1)=1
      NHEL(2)=-1
      NHEL(3)=glupolin
      NHEL(4)=-1
      NHEL(5)=1
      if (photpolin.eq.-1) then
      NHEL(6)=-1
      NHEL(7)=1
      else
      NHEL(6)=1
      NHEL(7)=-1
      endif
      do k=0,3
         if (incross.eq.1) then
         helastestmom(k,1)=p(k,1,1)
         helastestmom(k,2)=p(k,2,1)
         else
         helastestmom(k,2)=p(k,1,1)
         helastestmom(k,1)=p(k,2,1)
         endif
         helastestmom(k,3)=p(k,3,1)
         helastestmom(k,4)=v(k,1,1)
         helastestmom(k,5)=v(k,2,1)
         helastestmom(k,6)=v(k,3,1)
         helastestmom(k,7)=v(k,4,1)
      enddo
      bornubd(incross,glupolin,photpolin)=AMPUBDGEMVEBMUMMU(helastestmom,NHEL)
      enddo
      enddo
      enddo

      do incross=1,2
      do glupolin=-1,1,2
      do photpolin=-1,1,2
      bornubg(incross,glupolin,photpolin)=0d0
      NHEL(1)=1
      NHEL(2)=glupolin
      NHEL(3)=1
      NHEL(4)=-1
      NHEL(5)=1
      if (photpolin.eq.-1) then
      NHEL(6)=-1
      NHEL(7)=1
      else
      NHEL(6)=1
      NHEL(7)=-1
      endif
      do k=0,3
         if (incross.eq.1) then
         helastestmom(k,1)=p(k,1,1)
         helastestmom(k,2)=p(k,2,1)
         else
         helastestmom(k,2)=p(k,1,1)
         helastestmom(k,1)=p(k,2,1)
         endif
         helastestmom(k,3)=p(k,3,1)
         helastestmom(k,4)=v(k,1,1)
         helastestmom(k,5)=v(k,2,1)
         helastestmom(k,6)=v(k,3,1)
         helastestmom(k,7)=v(k,4,1)
      enddo

      if (incross.eq.1) then
          bornubg(incross,glupolin,photpolin)=AMPUBGDBEMVEBMUMM(helastestmom,NHEL)
      elseif (incross.eq.2) then
          bornubg(incross,glupolin,photpolin)=-AMPUBGDBEMVEBMUMM(helastestmom,NHEL)
      endif
      enddo
      enddo
      enddo

      do incross=1,2
      do glupolin=-1,1,2
      do photpolin=-1,1,2
      borndg(incross,glupolin,photpolin)=0d0 
      NHEL(1)=-1
      NHEL(2)=glupolin
      NHEL(3)=-1
      NHEL(4)=-1
      NHEL(5)=1
      if (photpolin.eq.-1) then
      NHEL(6)=-1
      NHEL(7)=1
      else
      NHEL(6)=1
      NHEL(7)=-1
      endif
      do k=0,3
         if (incross.eq.1) then
         helastestmom(k,1)=p(k,1,1)
         helastestmom(k,2)=p(k,2,1)
         else
         helastestmom(k,2)=p(k,1,1)
         helastestmom(k,1)=p(k,2,1)
         endif
         helastestmom(k,3)=p(k,3,1)
         helastestmom(k,4)=v(k,1,1)
         helastestmom(k,5)=v(k,2,1)
         helastestmom(k,6)=v(k,3,1)
         helastestmom(k,7)=v(k,4,1)
      enddo
      if (incross.eq.1) then
        borndg(incross,glupolin,photpolin)=-AMPDGUBEMVEBMUMMU(helastestmom,NHEL)
      elseif (incross.eq.2) then
        borndg(incross,glupolin,photpolin)=AMPDGUBEMVEBMUMMU(helastestmom,NHEL)
      endif
      enddo
      enddo
      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      print*, bornubd
!      print*, ''
!      print*, bornubg
!      print*, ''
!      print*, borndg
!      STOP


      im=dcmplx(0d0,1d0)
      cfc=dcmplx(cf) 	
      cac=dcmplx(ca)
      cacgraph=cac
      cfcac=(cfc-cac/2d0)
      pic2=dcmplx(pi**2)

      sumresult=0d0
      sumresultalt=0d0
      do mu=1,6
      m2s(mu)=0d0
      m2sborn(mu)=0d0
      m2sbornwopdf(mu)=0d0
      m2stest(mu)=0d0
      enddo
      mapspin1(1)=1
      mapspin1(2)=2
      mapspin2(1)=2
      mapspin2(2)=1

!       print*, mursq(1,1)
!       print*, mufsq(1,1)
!       print*, mufsq(2,1)

! ! pdfs
!       q_sf = sqrt(mufsq(1,1))
!       call pdfproton( xi(1), q_sf, pdf(-6,1) )        ! f_a(x1)
! 
!       q_sf = sqrt(mufsq(1,1))
!       call pdfproton( xi(2), q_sf, pdf(-6,2) )        ! f_b(x2)
! 
!       do k=-6,6,1
!         pdftmp(k,1)=pdf(k,1)/xi(1)
!         pdftmp(k,2)=pdf(k,2)/xi(2)
!       enddo
!       do k=-6,6,1
!         pdf(k,1)=pdftmp(sign1*k,1)
!         pdf(k,2)=pdftmp(sign2*k,2)
!       enddo
! 
! ! change pdfs for W+ case
!       if (pdfchange) then
!       do k=-6,6
!         pdftmp(k,1)=pdf(k,1)
!         pdftmp(k,2)=pdf(k,2)
!       enddo
!       do k=-6,6
!         pdf(k,1)=pdftmp(-k,1)
!         pdf(k,2)=pdftmp(-k,2)
!       enddo
!       endif

      map(1)=-1
      map(2)=1

      if (borndebug) then
      print*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      print*,"Born Debug Mode"
      print*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   ub d -> w gamma g
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      photin=1
      gluin=1
      do mu=0,3
      photmom(mu)=v(mu,3,1)+v(mu,4,1)
      glumom(mu)=p(mu,3,1)
      wdecmom(mu,1)=v(mu,2,1) !neutrino to become ket
      wdecmom(mu,2)=v(mu,1,1) !electron to become bra
      wmom(mu)=wdecmom(mu,1)+wdecmom(mu,2)
      wemnugamom(mu)=wmom(mu)+photmom(mu)
      enddo
      glumom(4)=0d0
      photmom(4)=dotrr(photmom,photmom)


      if (dotrr(v(0,3,1),v(0,3,1)).ge.1e-6) then
      print*, "photon", dotrr(v(0,3,1),v(0,3,1))
      endif
      if (dotrr(p(0,3,1),p(0,3,1)).ge.1e-6) then
      print*, "gluon", dotrr(p(0,3,1),p(0,3,1))
      endif


      do photpolin=-1,1,2
      do glupolin=-1,1,2
!
!  Polarizations using HELAS
!
!  effective photon/Z polarization
!
!
!  This is the Z->leptons current
!
      if (photpolin.eq.-1) then
      CALL OXXXXX(v(0,3,1   ),zero ,-1, 1, wmum ) !  mu-
      CALL IXXXXX(v(0,4,1   ),zero , 1,-1, wmup ) !  mu+
      elseif (photpolin.eq.1) then
      CALL OXXXXX(v(0,3,1   ),zero , 1, 1, wmum ) !  mu-
      CALL IXXXXX(v(0,4,1   ),zero ,-1,-1, wmup ) !  mu+
      endif

      do mu=0,3
      CALL VCARTX(photmom,ZMASS , Zwidth ,mu,-1, wavez(1,-1))  
      CALL IOVXXX(wmup,wmum,wavez(1,-1),GZL ,currlepz(mu,photpolin))
      enddo
!
!  This is the A->leptons current
!
      do mu=0,3
      CALL VCARTX(photmom,zero ,zero ,mu,-1, wavea(1,-1))
      CALL IOVXXX(wmup,wmum,wavea(1,-1),GAL ,currlepa(mu,photpolin))
      enddo
      
      do mu=0,3  ! only left-handed contributions
      storephotpol(mu,photpolin,1)=currlepa(mu,photpolin)+currlepz(mu,photpolin)/gau(1)*gzu(1) 
      storephotpol(mu,photpolin,2)=currlepa(mu,photpolin)+currlepz(mu,photpolin)/gad(1)*gzd(1)
      enddo

!
!   finally get the right effective photon polarization
!

      do mu=0,3
      photpol(mu,1)= storephotpol(mu,photpolin,1)
      photpol(mu,2)= storephotpol(mu,photpolin,2)
      enddo

!
!   effective WA current
!
      do mu=1,6
      enddo
      CALL OXXXXX(v(0,1,1   ),0d0 ,-1, 1, wem ) !  e-
      CALL IXXXXX(v(0,2,1   ),0d0 , 1,-1, wve ) !  ve~
      CALL WTOEVEMUPMUMvirt(photpolin,wem,wve,wmum,wmup,wemnugamom,currlepwga2(0,photpolin))

      do mu=0,3 
         currlepwga(mu,photpolin)=(currlepwga2(mu,photpolin)-
     &        (currlepwga2(0,photpolin)*wemnugamom(0)-currlepwga2(1,photpolin)*wemnugamom(1)
     &        -currlepwga2(2,photpolin)*wemnugamom(2)-currlepwga2(3,photpolin)*wemnugamom(3))
     &        *dcmplx(wemnugamom(mu))/dcmplx(wmass**2,-wmass*wwidth))
      enddo


      do mu=0,3
      storewphotpol(mu,photpolin)=currlepwga(mu,photpolin)
      enddo
       
!
!   gluon polarization
!
      CALL VXXXXX(glumom,0d0 ,glupolin,gluin, temp)
      do mu=0,3
      glupol(mu)=temp(mu+1)
      enddo
!       do mu=0,3                          !!!!! GAUGE CHECK
!       glupol(mu)=dcmplx(glumom(mu))
!       enddo
!
!   W current
!
      Call JIOXXX(wve,wem,gwf,wmass,wwidth,wmcurr)
      do mu=0,3
      wpol(mu)=wmcurr(mu+1)
      storewpol(mu)=wpol(mu)
      enddo
!
!     initial state exchanged
!
      do incross=1,2

      ampborn=dcmplx(0d0,0d0)
      ampbornubd(incross,glupolin,photpolin)=dcmplx(0d0,0d0)
      virt=dcmplx(0d0,0d0)
      fact=dcmplx(0d0,0d0)
      do mu=1,6
      result(mu)=dcmplx(0d0,0d0)
      resultb(mu)=dcmplx(0d0,0d0)
      enddo 
      do mu=1,8
      penresult(mu)=dcmplx(0d0,0d0)
      enddo
      do mu=1,5
      boxresult(mu)=dcmplx(0d0,0d0)
      enddo

!
!     SPINORS
!
      do mu=0,3
      spinormom(mu,1)= p(mu,mapspin2(incross),1) ! d to become ket
      spinormom(mu,2)= p(mu,mapspin1(incross),1) ! ub to become bra
      enddo
      signspinor(1) = 1
      signspinor(2) =-1
      CALL PSI0M(2,spinormom(0,1),signspinor,psi)
!
!   Pentagon Contribution
!
      do mu=0,3
      pentmom(mu,1)=1d0*p(mu,mapspin1(incross),1)   !ub incoming
      pentmom(mu,2)=1d0*p(mu,mapspin2(incross),1)   !d incoming
      pentmom(mu,3)=-1d0*glumom(mu)  !g outgoing
      pentmom(mu,4)=-1d0*photmom(mu) !gamma outgoing
      pentmom(mu,5)=-1d0*wmom(mu)    !w- outgoing
      enddo

      do mu=0,3
      tmp(mu)=wmom(mu)+photmom(mu)
      enddo
      xp2=dotcr(photpol(0,2),tmp)/dcmplx(dotrr(photmom,tmp))
      xp1=dotcr(photpol(0,1),tmp)/dcmplx(dotrr(photmom,tmp))
      xw=dotcr(wpol,tmp)/dcmplx(dotrr(wmom,tmp))
      xg=0d0!dotcr(glupol,tmp)/dcmplx(dotrr(glumom,tmp))
      do mu=0,3
         photpoltilde(mu,2)=photpol(mu,2)-xp2*dcmplx(photmom(mu),0d0)
         photpoltilde(mu,1)=photpol(mu,1)-xp1*dcmplx(photmom(mu),0d0)
         wpoltilde(mu)   =wpol(mu)-xw*dcmplx(wmom(mu),0d0)
         glupoltilde(mu) =glupol(mu)!-xg*dcmplx(glumom(mu),0d0)
      enddo

      s15=2d0*dotrr(pentmom(0,1),pentmom(0,2))
      cte= -(-3d0*(Log(mur) - Log(abs(s15)) + (0d0,1d0)*pi*theta(s15)) 
     &  - (Log(mur) - Log(abs(s15)) + (0d0,1d0)*pi*theta(s15))**2
c   substracted the c_virt
     &  -7d0 + 1d0/3d0*pi**2)
!
!   abelian
!
!       CALL  penlineABEmmm(pentmom(0,2),pentmom(0,3),pentmom(0,4),pentmom(0,5),pentmom(0,1),
!      &      psi(1,-1,2),psi(1,-1,1),glupol,photpol(0,2),wpol,
!      &      -1,    mur,    1,        2  ,    4      , 1 ,resultgauge(1),
!      &      penresult,resultgaugeb(1),resultb(1))
! 
!       result(1)=(cfc*penresult(1)+cfcac*penresult(2))*gad(1)
!       resultb(1)=resultb(1)*gad(1)  ! gives +1*Madgraph
!       test(1)=result(1)!+cte*resultb(1)
!       testb(1)=resultb(1) 
! 
!       call pent_gauge(pentmom(0,2),pentmom(0,3),pentmom(0,4),pentmom(0,5),pentmom(0,1),
!      &      psi(1,-1,2),psi(1,-1,1),glupol,photpol(0,2),wpol,
!      &      -1,    mur ,resultgauge,resultgaugeb,wardresult)
!       print*, result(1),resultb(1)
     

      CALL  penlineABEmmm(pentmom(0,2),pentmom(0,3),pentmom(0,4),pentmom(0,5),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),glupoltilde,photpoltilde(0,2),wpoltilde,
     &      -1,    mur,    1,        2  ,    4      , 1 ,resultgauge(1),
     &      penresult,resultgaugeb(1),resultb(1))

      call pent_to_boxnoab(pentmom(0,2),pentmom(0,3),pentmom(0,4),pentmom(0,5),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),glupol,photpol(0,2),wpol,xg,xp2,xw,-1,mur,2,testpb1,testpb2)

      call pent_gauge(pentmom(0,2),pentmom(0,3),pentmom(0,4),pentmom(0,5),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),glupoltilde,photpoltilde(0,2),wpoltilde,
     &      -1,    mur ,resultgauge,resultgaugeb,wardresult)

!       if (wardresult.and.nlopent) then
!       result(1)=(cfc*penresult(1)+cfcac*penresult(2))*gad(1)
!       endif
!       if (nlobox) then
!       result(1)=testpb1*gad(1)
!       endif

      if (wardresult) then
      result(1)=(cfc*penresult(1)+cfcac*penresult(2)+testpb1)*gad(1)
      resultb(1)=(resultb(1)+testpb2)*gad(1)  ! gives +1*Madgraph
      else
      result(1)=(testpb1)*gad(1)
      resultb(1)=(resultb(1)+testpb2)*gad(1)  ! gives +1*Madgraph
      endif

      if (ldebug) then !.and.nlopent
!       if (abs(testpb1).eq.0d0) then
!       ubd(incross)=.false.
!       endif
      if (.not.wardresult) then
      badpointspent(incross)=badpointspent(incross)+1
      ubd(incross)=.false.
      endif
      endif

!       CALL  penlineABEmmm(pentmom(0,2),pentmom(0,4),pentmom(0,5),pentmom(0,3),pentmom(0,1),
!      &      psi(1,-1,2),psi(1,-1,1),photpol(0,2),wpol,glupol,
!      &      -1,mur,1,4,4,1,resultgauge(1),penresult,resultgaugeb(1),resultb(2))
! 
!       result(2)=(cfc*penresult(1)+cfcac*penresult(2))*gad(1)
!       resultb(2)=resultb(2)*gad(1)
!       test(2)=result(2)!+cte*resultb(2)
!       testb(2)=resultb(2)!+cte*resultb(2)
! 
!       call pent_gauge(pentmom(0,2),pentmom(0,4),pentmom(0,5),pentmom(0,3),pentmom(0,1),
!      &      psi(1,-1,2),psi(1,-1,1),photpol(0,2),wpol,glupol,
!      &      -1,mur,resultgauge,resultgaugeb,wardresult)
!       print*, result(2),resultb(2)

      CALL penlineABEmmm(pentmom(0,2),pentmom(0,4),pentmom(0,5),pentmom(0,3),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),photpoltilde(0,2),wpoltilde,glupoltilde,
     &      -1,mur,1,4,4,1,resultgauge(1),penresult,resultgaugeb(1),resultb(2))

      CALL pent_to_boxnoab(pentmom(0,2),pentmom(0,4),pentmom(0,5),pentmom(0,3),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),photpol(0,2),wpol,glupol,xp2,xw,xg,-1,mur,3,testpb1,testpb2)

      call pent_gauge(pentmom(0,2),pentmom(0,4),pentmom(0,5),pentmom(0,3),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),photpoltilde(0,2),wpoltilde,glupoltilde,
     &      -1,mur,resultgauge,resultgaugeb,wardresult)

!       if (wardresult.and.nlopent) then
!       result(2)=(cfc*penresult(1)+cfcac*penresult(2))*gad(1)
!       endif
!       if (nlobox) then
!       result(2)=testpb1*gad(1)
!       endif

      if (wardresult) then 
      result(2)=(cfc*penresult(1)+cfcac*penresult(2)+testpb1)*gad(1)
      resultb(2)=(resultb(2)+testpb2)*gad(1)  ! gives +1*Madgraph
      else
      result(2)=(testpb1)*gad(1)
      resultb(2)=(resultb(2)+testpb2)*gad(1)  ! gives +1*Madgraph
      endif

      if (ldebug) then !.and.nlopent
!       if (abs(testpb1).eq.0d0) then
!       ubd(incross)=.false.
!       endif
      if (.not.wardresult) then
      badpointspent(incross)=badpointspent(incross)+1 
      ubd(incross)=.false.
      endif
      endif

!       CALL  penlineABEmmm(pentmom(0,2),pentmom(0,5),pentmom(0,3),pentmom(0,4),pentmom(0,1),
!      &      psi(1,-1,2),psi(1,-1,1),wpol,glupol,photpol(0,1),
!      &      -1,mur,1,3,4,1,resultgauge(1),penresult,resultgaugeb(1),resultb(3))
! 
!       result(3)=(cfc*penresult(1)+cfcac*penresult(2))*gau(1)
!       resultb(3)=resultb(3)*gau(1)    ! gives +1*HELAS
!       test(3)=result(3)!+cte*resultb(3)
!       testb(3)=resultb(3)
!       CALL pent_gauge(pentmom(0,2),pentmom(0,5),pentmom(0,3),pentmom(0,4),pentmom(0,1),
!      &      psi(1,-1,2),psi(1,-1,1),wpol,glupol,photpol(0,1),
!      &      -1,    mur,resultgauge,resultgaugeb,wardresult)
!       print*, result(3),resultb(3)
     

      CALL  penlineABEmmm(pentmom(0,2),pentmom(0,5),pentmom(0,3),pentmom(0,4),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),wpoltilde,glupoltilde,photpoltilde(0,1),
     &      -1,mur,1,3,4,1,resultgauge(1),penresult,resultgaugeb(1),resultb(3))

      CALL pent_to_boxnoab(pentmom(0,2),pentmom(0,5),pentmom(0,3),pentmom(0,4),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),wpol,glupol,photpol(0,1),xw,xg,xp1,-1,mur,6,testpb1,testpb2)

      CALL pent_gauge(pentmom(0,2),pentmom(0,5),pentmom(0,3),pentmom(0,4),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),wpoltilde,glupoltilde,photpoltilde(0,1),
     &      -1,    mur,resultgauge,resultgaugeb,wardresult)

!       if (wardresult.and.nlopent) then
!       result(3)=(cfc*penresult(1)+cfcac*penresult(2))*gau(1)
!       endif
!       if (nlobox) then
!       result(3)=(testpb1)*gau(1)
!       endif

      if (wardresult) then
      result(3)=(cfc*penresult(1)+cfcac*penresult(2)+testpb1)*gau(1)
      resultb(3)=(resultb(3)+testpb2)*gau(1)  ! gives +1*Madgraph
      else
      result(3)=(testpb1)*gau(1)
      resultb(3)=(resultb(3)+testpb2)*gau(1)  ! gives +1*Madgraph
      endif

      if (ldebug) then !.and.nlopent
!       if (abs(testpb1).eq.0d0) then
!       ubd(incross)=.false.
!       endif
      if (.not.wardresult) then
      badpointspent(incross)=badpointspent(incross)+1 
      ubd(incross)=.false.
      endif
      endif

!       CALL penlineABEmmm(pentmom(0,2),pentmom(0,4),pentmom(0,3),pentmom(0,5),pentmom(0,1),
!      &      psi(1,-1,2),psi(1,-1,1),photpol(0,2),glupol,wpol,
!      &      -1,mur,1,3,4,1,resultgauge(1),penresult,resultgaugeb(1),resultb(4))
! 
!       result(4)=(cfc*penresult(1)+cfcac*penresult(2))*gad(1)
!       resultb(4)=resultb(4)*gad(1)
!       test(4)=result(4)!+cte*resultb(4)
!       testb(4)=resultb(4)
!       CALL pent_gauge(pentmom(0,2),pentmom(0,4),pentmom(0,3),pentmom(0,5),pentmom(0,1),
!      &      psi(1,-1,2),psi(1,-1,1),photpol(0,2),glupol,wpol,
!      &      -1,    mur ,resultgauge,resultgaugeb,wardresult)
!       print*, result(4),resultb(4)

      CALL penlineABEmmm(pentmom(0,2),pentmom(0,4),pentmom(0,3),pentmom(0,5),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),photpoltilde(0,2),glupoltilde,wpoltilde,
     &      -1,mur,1,3,4,1,resultgauge(1),penresult,resultgaugeb(1),resultb(4))

      CALL pent_to_boxnoab(pentmom(0,2),pentmom(0,4),pentmom(0,3),pentmom(0,5),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),photpol(0,2),glupol,wpol,xp2,xg,xw,-1,mur,6,testpb1,testpb2)

      CALL pent_gauge(pentmom(0,2),pentmom(0,4),pentmom(0,3),pentmom(0,5),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),photpoltilde(0,2),glupoltilde,wpoltilde,
     &      -1,    mur ,resultgauge,resultgaugeb,wardresult)

      if (wardresult) then
      result(4)=(cfc*penresult(1)+cfcac*penresult(2)+testpb1)*gad(1)
      resultb(4)=(resultb(4)+testpb2)*gad(1)  ! gives +1*Madgraph
      else
      result(4)=(testpb1)*gad(1)
      resultb(4)=(resultb(4)+testpb2)*gad(1)  ! gives +1*Madgraph
      endif

!       if (wardresult.and.nlopent) then
!       result(4)=(cfc*penresult(1)+cfcac*penresult(2))*gad(1)
!       endif
!       if (nlobox) then
!       result(4)=(testpb1)*gad(1)
!       endif

      if (ldebug) then !.and.nlopent
!       if (abs(testpb1).eq.0d0) then
!       ubd(incross)=.false.
!       endif
      if (.not.wardresult) then
      badpointspent(incross)=badpointspent(incross)+1
      ubd(incross)=.false.
      endif
      endif

!       CALL penlineABEmmm(pentmom(0,2),pentmom(0,3),pentmom(0,5),pentmom(0,4),pentmom(0,1),
!      &      psi(1,-1,2),psi(1,-1,1),glupol,wpol,photpol(0,1),
!      &      -1,mur,1,2,4,1,resultgauge(1),penresult,resultgaugeb(1),resultb(5))
! 
!       result(5)=(cfc*penresult(1)+cfcac*penresult(2))*gau(1)
!       resultb(5)=resultb(5)*gau(1)
!       test(5)=result(5)!+cte*resultb(5)
!       testb(5)=resultb(5)
!       call pent_gauge(pentmom(0,2),pentmom(0,3),pentmom(0,5),pentmom(0,4),pentmom(0,1),
!      &      psi(1,-1,2),psi(1,-1,1),glupol,wpol,photpol(0,1),
!      &      -1,    mur,resultgauge,resultgaugeb,wardresult)
!       print*, result(5),resultb(5)

      CALL penlineABEmmm(pentmom(0,2),pentmom(0,3),pentmom(0,5),pentmom(0,4),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),glupoltilde,wpoltilde,photpoltilde(0,1),
     &      -1,mur,1,2,4,1,resultgauge(1),penresult,resultgaugeb(1),resultb(5))

      CALL pent_to_boxnoab(pentmom(0,2),pentmom(0,3),pentmom(0,5),pentmom(0,4),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),glupol,wpol,photpol(0,1),xg,xw,xp1,-1,mur,2,testpb1,testpb2)

      call pent_gauge(pentmom(0,2),pentmom(0,3),pentmom(0,5),pentmom(0,4),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),glupoltilde,wpoltilde,photpoltilde(0,1),
     &      -1,    mur,resultgauge,resultgaugeb,wardresult)

      if (wardresult) then
      result(5)=(cfc*penresult(1)+cfcac*penresult(2)+testpb1)*gau(1)
      resultb(5)=(resultb(5)+testpb2)*gau(1)  ! gives +1*Madgraph
      else
      result(5)=(testpb1)*gau(1)
      resultb(5)=(resultb(5)+testpb2)*gau(1)  ! gives +1*Madgraph
      endif

!       if (wardresult.and.nlopent) then
!       result(5)=(cfc*penresult(1)+cfcac*penresult(2))*gau(1)
!       endif
!       if (nlobox) then
!       result(5)=(testpb1)*gau(1)
!       endif

      if (ldebug) then  !.and.nlopent
!       if (abs(testpb1).eq.0d0) then
!       ubd(incross)=.false.
!       endif
      if (.not.wardresult) then
      badpointspent(incross)=badpointspent(incross)+1 
      ubd(incross)=.false.
      endif
      endif
!                             p1           p2           p3           p4           p5    
!       CALL  penlineABEmmm(pentmom(0,2),pentmom(0,5),pentmom(0,4),pentmom(0,3),pentmom(0,1),
!      &      psi(1,-1,2),psi(1,-1,1),wpol,photpol(0,1),glupol,
!      &      -1,mur,1,4,4,1,resultgauge(1),penresult,resultgaugeb(1),resultb(6))
! 
!       result(6)=(cfc*penresult(1)+cfcac*penresult(2))*gau(1)
!       resultb(6)=resultb(6)*gau(1)
!       test(6)=result(6)!+cte*resultb(6)
!       testb(6)=resultb(6)
!       CALL pent_gauge(pentmom(0,2),pentmom(0,5),pentmom(0,4),pentmom(0,3),pentmom(0,1),
!      &      psi(1,-1,2),psi(1,-1,1),wpol,photpol(0,1),glupol,
!      &      -1,    mur,resultgauge,resultgaugeb,wardresult)
!       print*, result(6),resultb(6)

      CALL  penlineABEmmm(pentmom(0,2),pentmom(0,5),pentmom(0,4),pentmom(0,3),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),wpoltilde,photpoltilde(0,1),glupoltilde,
     &      -1,mur,1,4,4,1,resultgauge(1),penresult,resultgaugeb(1),resultb(6))

      CALL pent_to_boxnoab(pentmom(0,2),pentmom(0,5),pentmom(0,4),pentmom(0,3),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),wpol,photpol(0,1),glupol,xw,xp1,xg,-1,mur,3,testpb1,testpb2)

      CALL pent_gauge(pentmom(0,2),pentmom(0,5),pentmom(0,4),pentmom(0,3),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),wpoltilde,photpoltilde(0,1),glupoltilde,
     &      -1,    mur,resultgauge,resultgaugeb,wardresult)

      if (wardresult) then
      result(6)=(cfc*penresult(1)+cfcac*penresult(2)+testpb1)*gau(1)
      resultb(6)=(resultb(6)+testpb2)*gau(1)  ! gives +1*Madgraph
      else
      result(6)=(testpb1)*gau(1)
      resultb(6)=(resultb(6)+testpb2)*gau(1)  ! gives +1*Madgraph
      endif

!       if (wardresult.and.nlopent) then
!       result(6)=(cfc*penresult(1)+cfcac*penresult(2))*gau(1)
!       endif
!       if (nlobox) then
!       result(6)=(testpb1)*gau(1)
!       endif

      if (ldebug) then !.and.nlopent
!       if (abs(testpb1).eq.0d0) then
!       ubd(incross)=.false.
!       endif
      if (.not.wardresult) then
      badpointspent(incross)=badpointspent(incross)+1 
      ubd(incross)=.false.
      endif
      endif

!  sum the abelian pentagon contributions
      virt=virt+(result(1)+result(2)+result(3)+result(4)+result(5)+result(6))*gg(1)**3*gwf(1)
!  sum the born contributions
      ampborn=ampborn+(resultb(1)+resultb(2)+resultb(3)+resultb(4)+resultb(5)+resultb(6))*gg(1)*gwf(1)

!
!  non-abelian
!
      CALL  penlineNoAbemmm(pentmom(0,2),pentmom(0,4),pentmom(0,5),pentmom(0,1),pentmom(0,3),
     &      psi(1,-1,2),psi(1,-1,1),photpol(0,2),wpol,glupol,
     &      -1,    mur,     4      , 1 ,resultgauge(1),result(1))

      CALL pentNoAbe_gauge(accuracy,pentmom(0,2),pentmom(0,4),pentmom(0,5),pentmom(0,1),pentmom(0,3),
     &      psi(1,-1,2),psi(1,-1,1),photpol(0,2),wpol,glupol,
     &      -1,    mur, resultgauge(1), accuracyres,Wardresult)

      if (ldebug) then
      if (.not.wardresult) then
      badpointspentnoab(incross)=badpointspentnoab(incross)+1 
!      sumresult=0d0
!      return
      ubd(incross)=.false.
      endif
      endif

      result(1)=cacgraph*result(1)*gad(1)
      test(7)=result(1)

      CALL  penlineNoAbemmm(pentmom(0,2),pentmom(0,5),pentmom(0,4),pentmom(0,1),pentmom(0,3),
     &      psi(1,-1,2),psi(1,-1,1),wpol,photpol(0,1),glupol,
     &      -1,    mur,     4      , 1 ,resultgauge(1),result(2))

      CALL pentNoAbe_gauge(accuracy,pentmom(0,2),pentmom(0,5),pentmom(0,4),pentmom(0,1),pentmom(0,3),
     &      psi(1,-1,2),psi(1,-1,1),wpol,photpol(0,1),glupol,
     &      -1,    mur, resultgauge(1), accuracyres,Wardresult)

      if (ldebug) then
      if (.not.wardresult) then
      badpointspentnoab(incross)=badpointspentnoab(incross)+1 
!      sumresult=0d0
!      return
      ubd(incross)=.false.
      endif
      endif

      result(2)=cacgraph*result(2)*gau(1)
      test(8)=result(2)

!  sum the non-abelian pentagon contributions
      virt=virt+(result(1)+result(2))*gg(1)**3*gwf(1)
!
!   Box Contribution
!
      do mu=0,3
      boxmom(mu,1)=pentmom(mu,1)   !ub incoming
      boxmom(mu,2)=pentmom(mu,2)   !d incoming
      boxmom(mu,3)=-1d0*glumom(mu)  !g outgoing
      boxmom(mu,4)=-1d0*wemnugamom(mu) !Wgamma outgoing
      polgluontemp(mu)=boxmom(mu,3)
      enddo

!
!   abelian
!
!       if (nlobox) then
!    note: the second call is needed for the gauge test of the non-abelian boxes
      CALL boxlineABEmm(boxmom(0,2),boxmom(0,3),boxmom(0,4),boxmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),glupol,currlepwga(0,photpolin),
     &      -1, mur, 1, 2 ,3,1,resultgauge(1),boxresult,
     &      resultgaugeb(1),resultb(1))
      CALL boxlineABEmm(boxmom(0,2),boxmom(0,3),boxmom(0,4),boxmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),polgluontemp,currlepwga(0,photpolin),
     &      -1, mur, 1, 2 ,3,-1,nothing,gaugeCA1,
     &      nothing,nothing(1))

      result(1)=-1d0*(cfc*boxresult(1)+cfcac*boxresult(2))
      resultb(1)=-1d0*resultb(1) ! boxlineABEmm gives -1*HELAS for Born
      test(9)=result(1)!+cte*resultb(1)
      testb(9)=resultb(1)
!
!   gauge test for boxes
!
      if (abs(resultgaugeb(1)).le.1e-10) then
      testgauge1=resultgauge(1)+cte*resultgaugeb(1)
      else
      testgauge1=resultgauge(1)/(-cte*resultgaugeb(1))-1.
      endif
      if (abs(resultgaugeb(2)).le.1e-10) then
      testgauge2=resultgauge(2)+cte*resultgaugeb(2)
      else
      testgauge2=resultgauge(2)/(-cte*resultgaugeb(2))-1.
      endif
      if (ldebug) then
      if (abs(testgauge1).ge.1e-2) then
!      print*, "here ub d, box 1, g1", "here???",abs(testgauge1)
      STOP
      badpointsbox(incross)=badpointsbox(incross)+1
!      sumresult=0d0
!      return
      ubd(incross)=.false.
      endif
      if (abs(testgauge2).ge.1e-2) then
!      print*, "here ub d, box 1, g2"
      badpointsbox(incross)=badpointsbox(incross)+1
!      sumresult=0d0
!      return
      ubd(incross)=.false.
      endif
      endif

!    note: the second call is needed for the gauge test of the non-abelian boxes
      CALL boxlineABEmm(boxmom(0,2),boxmom(0,4),boxmom(0,3),boxmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),currlepwga(0,photpolin),glupol,
     &      -1, mur, 1, 3 ,3,1,resultgauge2(1),boxresult,
     &      resultgaugeb(1),resultb(2))
      CALL boxlineABEmm(boxmom(0,2),boxmom(0,4),boxmom(0,3),boxmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),currlepwga(0,photpolin),polgluontemp,
     &      -1, mur, 1, 3 ,3,-1,nothing,gaugeCA2,
     &      nothing,nothing(1))
      
      result(2)=-1d0*(cfc*boxresult(1)+cfcac*boxresult(2))
      resultb(2)=-1d0*resultb(2) ! boxlineABEmm gives -1*HELAS for Born
      test(10)=result(2)!+cte*resultb(2)
      testb(10)=resultb(2)

!
!   gauge test for boxes
!
      if (ldebug) then
      if (abs(resultgaugeb(1)).le.1e-10) then
      testgauge1=resultgauge2(1)+cte*resultgaugeb(1)
      else
      testgauge1=resultgauge2(1)/(-cte*resultgaugeb(1))-1.
      endif
      if (abs(resultgaugeb(2)).le.1e-10) then
      testgauge2=resultgauge2(2)+cte*resultgaugeb(2)
      else
      testgauge2=resultgauge2(2)/(-cte*resultgaugeb(2))-1.
      endif

      if (abs(testgauge1).ge.1e-2) then
!      print*, "here ub d, box 2, g1"
      badpointsbox(incross)=badpointsbox(incross)+1
!      sumresult=0d0
!      return
      endif
      if (abs(testgauge2).ge.1e-2) then
!      print*, "here ub d, box 2, g2"
      badpointsbox(incross)=badpointsbox(incross)+1
!      sumresult=0d0
!      return
      ubd(incross)=.false.
      endif
      endif

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Non-abelian Part
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C W
      CALL boxlineNoAbemm(boxmom(0,2),boxmom(0,4),boxmom(0,1),boxmom(0,3),
     &      psi(1,-1,2),psi(1,-1,1),currlepwga(0,photpolin),glupol,
     &      -1, mur, 3, 1,resultgauge3(1),result(3))

** gauge test:      
      if (abs(resultgauge3(2)).lt.1d-7) then
           testgauge2=abs(-(gaugeCA1(2)+ gaugeCA2(2))/2d0- resultgauge3(2) )
      else
           testgauge2=abs(-(gaugeCA1(2)+ gaugeCA2(2))/2d0/resultgauge3(2)+1d0)
      endif
      if (abs(testgauge2).ge.1d-1) then
         badpointsbox(incross)=badpointsbox(incross)+1
         ubd(incross)=.false.
      endif

      result(3)=-1d0*result(3)*cacgraph
      test(11)=result(3)

      virt=virt+(result(1)+result(2)+result(3))*gwf(1)*gg(1)**3
      ampborn=ampborn+(resultb(1)+resultb(2))*gwf(1)*gg(1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            Finite Terms from Pacos expansion     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      s=2.*dotrr(pentmom(0,1),pentmom(0,2))
      t=2.*dotrr(pentmom(0,1),pentmom(0,3))
      u=2.*dotrr(pentmom(0,2),pentmom(0,3))

       factfigy=-2*EulerGamma*Inv3*nf*TF + 2*Inv3*nf*TF*Log(4*Pi) + 
     -     CF*(-1 + 3*EulerGamma - EulerGamma**2 - Inv2*Inv3*Pi**2 + 
     -     Log(musq)*(2*EulerGamma - 2*Log(4*Pi)) - 
     -     Log(4*Pi)**2 - 2*EulerGamma*Log(Abs(s)) + (0,2)*
     -     EulerGamma*Pi*theta(s) + Log(4*Pi)*(-3 + 2*EulerGamma + 2*
     -     Log(Abs(s)) - (0,2)*Pi*theta(s))) + 
     -  CA*((11 - 3*EulerGamma)*EulerGamma*Inv2*Inv3 - Inv2**2*Inv3*
     -     Pi**2 + Log(musq)*(EulerGamma - Log(4*Pi)) - Inv2*
     -     Log(4*Pi)**2 + EulerGamma*Log(Abs(s)) - 
     -     EulerGamma*Log(Abs(t)) - EulerGamma*Log(Abs(u)) + (0,1)*
     -     EulerGamma*Pi*(-theta(s) + theta(t) + theta(u)) + 
     -     Log(4*Pi)*(EulerGamma - 11*Inv2*Inv3 - Log(Abs(s)) + 
     -     Log(Abs(t)) + Log(Abs(u)) - (0,1)*Pi*(-theta(s) + theta(t) + 
     -     theta(u))))

      virt=virt
     &      +factfigy*ampborn*gg(1)**2


!       endif !!! BOX contributions

      ampbornubd(incross,glupolin,photpolin)=bornubd(incross,glupolin,
     -     photpolin)!ampborn
      ampborn=bornubd(incross,glupolin,photpolin)
      m2s(incross)= m2s(incross)+2d0*dreal(virt*conjg(ampborn))/
     -     (4d0*pi)**2/9d0
      m2sborn(incross)=m2sborn(incross)+abs(ampborn)**2/9d0 
      m2stest(incross)=m2stest(incross)+2d0*dreal(((test(1)+test(2)+
     -     test(3)+test(4)+test(5)
     &    +test(6)+test(7)+test(8)+test(9)+test(10)+test(11))*gg(1)**3*
     -     gwf(1))*conjg(ampborn))/9d0


      enddo  !loop over initial state exchanged
      enddo  !loop over gluon polarization
      enddo  !loop over photon polarization

      m2s(1)=m2s(1)*(pdf(-2,1)*pdf(1,2)+pdf(-4,1)*pdf(3,2))
      m2s(2)=m2s(2)*(pdf(-2,2)*pdf(1,1)+pdf(-4,2)*pdf(3,1))
      m2stest(1)=m2stest(1)*(pdf(-2,1)*pdf(1,2)+pdf(-4,1)*pdf(3,2))
      m2stest(2)=m2stest(2)*(pdf(-2,2)*pdf(1,1)+pdf(-4,2)*pdf(3,1))
      m2sbornwopdf(1)=m2sborn(1)
      m2sbornwopdf(2)=m2sborn(2)
      m2sborn(1)=m2sborn(1)*(pdf(-2,1)*pdf(1,2)+pdf(-4,1)*pdf(3,2))
      m2sborn(2)=m2sborn(2)*(pdf(-2,2)*pdf(1,1)+pdf(-4,2)*pdf(3,1))

      if (borndebug) then

      do k=0,3
         helastestmom(k,1)=p(k,1,1)
         helastestmom(k,2)=p(k,2,1)
         helastestmom(k,3)=p(k,3,1)
         helastestmom(k,4)=v(k,1,1)
         helastestmom(k,5)=v(k,2,1)
         helastestmom(k,6)=v(k,3,1)
         helastestmom(k,7)=v(k,4,1)
      enddo

! u~ d -> g a e- ve~ 

      Call SUBDGEMVEBMUMMUP(helastestmom,testh)
      testh=testh*(pdf(-2,1)*pdf(1,2)+pdf(-4,1)*pdf(3,2))
      print*,  "born from virtuals/HELAS = ",m2sborn(1)/testh, "  =1 ?"


      do k=0,3
         helastestmom(k,1)=p(k,2,1)
         helastestmom(k,2)=p(k,1,1)
         helastestmom(k,3)=p(k,3,1)
         helastestmom(k,4)=v(k,1,1)
         helastestmom(k,5)=v(k,2,1)
         helastestmom(k,6)=v(k,3,1)
         helastestmom(k,7)=v(k,4,1)
      enddo

      Call SUBDGEMVEBMUMMUP(helastestmom,testh)
      testh=testh*(pdf(-2,2)*pdf(1,1)+pdf(-4,2)*pdf(3,1))

      do mu=0,3
      boxmomtilde(mu,1)=boxmom(mu,1) !INCOMING
      boxmomtilde(mu,2)=-boxmom(mu,2) !outgoing
      boxmomtilde(mu,3)=-boxmom(mu,3) !outgoing
      boxmomtilde(mu,4)=-boxmom(mu,4) !outgoing
      enddo


      print*,  "born from virtuals/HELAS = ",m2sborn(2)/testh, "  =1 ?"
      print*,  mur
      print*, "-s23",-2*dotrr(boxmomtilde(0,2),boxmomtilde(0,3))  !
      print*, "-t", 2*dotrr(boxmomtilde(0,1),boxmomtilde(0,3))   !  additional minus due to incoming conv
      print*, "-t12",2*dotrr(boxmomtilde(0,1),boxmomtilde(0,2))  !
      print*, "finished ub d process"

      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   ub g -> w gamma db global minus compared to ub d/Madgraph
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!       cte=(-1.0d0*pic2/3.0d0+7.0d0)
      photin=1
      gluin=-1
!
!     initial state exchanged
!
      do incross=1,2

      do mu=0,3
      glumom(mu)=p(mu,mapspin2(incross),1)
      enddo

      do photpolin=-1,1,2
      do glupolin=-1,1,2

!
!  Polarizations using HELAS
!
!   photon polarization
!
      do mu=0,3
      photpol(mu,1)=storephotpol(mu,photpolin,1)
      photpol(mu,2)=storephotpol(mu,photpolin,2)
      currlepwga(mu,photpolin)=storewphotpol(mu,photpolin)
      wpol(mu)=storewpol(mu)
      enddo

!
!   gluon polarization
!
      CALL VXXXXX(glumom,0d0 ,glupolin,gluin, temp)
      do mu=0,3
      glupol(mu)=temp(mu+1)
      enddo

      ampborn=dcmplx(0d0,0d0)
      ampbornubg(incross,glupolin,photpolin)=dcmplx(0d0,0d0)
      ampborn=dcmplx(0d0,0d0)
      virt=dcmplx(0d0,0d0)
      fact=dcmplx(0d0,0d0)
      do mu=1,6
      result(mu)=dcmplx(0d0,0d0)
      resultb(mu)=dcmplx(0d0,0d0)
      enddo 
      do mu=1,8
      penresult(mu)=dcmplx(0d0,0d0)
      enddo
      do mu=1,5
      boxresult(mu)=dcmplx(0d0,0d0)
      enddo
!
!     SPINORS
!
      do mu=0,3
      spinormom(mu,1)= p(mu,3,1)  ! db to become ket
      spinormom(mu,2)= p(mu,mapspin1(incross),1) ! ub to become bra
      enddo
      signspinor(1) =-1
      signspinor(2) =-1
      CALL PSI0M(2,spinormom(0,1),signspinor,psi)

!
!   Pentagon Contribution
!
      do mu=0,3
      pentmom(mu,1)=1d0*p(mu,mapspin1(incross),1)   !ub incoming
      pentmom(mu,2)=-1d0*p(mu,3,1)   !db outgouing
      pentmom(mu,3)=1d0*glumom(mu)   !g incoming
      pentmom(mu,4)=-1d0*photmom(mu) !gamma outgoing
      pentmom(mu,5)=-1d0*wmom(mu)    !w- outgoing
      enddo

      do mu=0,3
      tmp(mu)=wmom(mu)+photmom(mu)
      enddo
      xp1=dotcr(photpol(0,1),tmp)/dcmplx(dotrr(photmom,tmp))
      xp2=dotcr(photpol(0,2),tmp)/dcmplx(dotrr(photmom,tmp))
      xw=dotcr(wpol,tmp)/dcmplx(dotrr(wmom,tmp))
      xg=0d0!dotcr(glupol,tmp)/dcmplx(dotrr(glumom,tmp))
      do mu=0,3
         photpoltilde(mu,2)=photpol(mu,2)-xp2*dcmplx(photmom(mu),0d0)
         photpoltilde(mu,1)=photpol(mu,1)-xp1*dcmplx(photmom(mu),0d0)
         wpoltilde(mu)   =wpol(mu)-xw*dcmplx(wmom(mu),0d0)
         glupoltilde(mu) =glupol(mu)!-xg*dcmplx(glumom(mu),0d0)
      enddo

      s15=2d0*dotrr(pentmom(0,1),pentmom(0,2))
      cte= -(-3d0*(Log(mur) - Log(abs(s15)) + (0d0,1d0)*pi*theta(s15)) 
     &  - (Log(mur) - Log(abs(s15)) + (0d0,1d0)*pi*theta(s15))**2
c   substracted the c_virt
     &  -7d0 + 1d0/3d0*pi**2)
!
!   abelian
!
!       CALL  penlineABEmmm(pentmom(0,2),pentmom(0,3),pentmom(0,4),pentmom(0,5),pentmom(0,1),
!      &      psi(1,-1,2),psi(1,-1,1),glupol,photpol(0,2),wpol,
!      &      -1,    mur,    1,        2  ,    4      , 1 ,resultgauge(1),
!      &      penresult,resultgaugeb(1),resultb(1))
! 
!       result(1)=(cfc*penresult(1)+cfcac*penresult(2))*gad(1)
!       resultb(1)=resultb(1)*gad(1)  ! gives -1*Madgraph
!       test(1)=result(1)+cte*resultb(1)
! 
!       call pent_gauge(pentmom(0,2),pentmom(0,3),pentmom(0,4),pentmom(0,5),pentmom(0,1),
!      &      psi(1,-1,2),psi(1,-1,1),glupol,photpol(0,2),wpol,
!      &      -1,    mur ,resultgauge,resultgaugeb,wardresult)
!       print*, "here 1",result(1),resultb(1), wardresult

      CALL  penlineABEmmm(pentmom(0,2),pentmom(0,3),pentmom(0,4),p
     -     entmom(0,5),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),glupoltilde,photpoltilde(0,2),wpoltilde,
     &      -1,    mur,    1,        2  ,    4      , 1 ,resultgauge(1),
     &      penresult,resultgaugeb(1),resultb(1))

      call pent_to_boxnoab(pentmom(0,2),pentmom(0,3),pentmom(0,4),
     -     pentmom(0,5),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),glupol,photpol(0,2),wpol,xg,xp2,
     -     xw,-1,mur,2,testpb1,testpb2)

      call pent_gauge(pentmom(0,2),pentmom(0,3),pentmom(0,4),
     -     pentmom(0,5),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),glupoltilde,photpoltilde(0,2),wpoltilde,
     &      -1,    mur ,resultgauge,resultgaugeb,wardresult)

      if (wardresult) then
      result(1)=(cfc*penresult(1)+cfcac*penresult(2)+testpb1)*gad(1)
      resultb(1)=(resultb(1)+testpb2)*gad(1)  ! gives +1*Madgraph
      else
      result(1)=(testpb1)*gad(1)
      resultb(1)=(resultb(1)+testpb2)*gad(1)  ! gives +1*Madgraph
      endif

      if (ldebug) then
!       if (abs(testpb1).eq.0d0) then
!       ubg(incross)=.false.
!       endif
!       print*, "penline gauge check gluon-photon-W", wardresult
      if (.not.wardresult) then
!       print*, "bad point pent 1"
      badpointspent(incross+2)=badpointspent(incross+2)+1
!      sumresult=0d0
!      return
      ubg(incross)=.false.
      endif
      endif

!       CALL  penlineABEmmm(pentmom(0,2),pentmom(0,4),pentmom(0,5),pentmom(0,3),pentmom(0,1),
!      &      psi(1,-1,2),psi(1,-1,1),photpol(0,2),wpol,glupol,
!      &      -1,mur,1,4,4,1,resultgauge(1),penresult,resultgaugeb(1),resultb(2))
! 
!       result(2)=(cfc*penresult(1)+cfcac*penresult(2))*gad(1)
!       resultb(2)=resultb(2)*gad(1)
!       test(2)=result(2)+cte*resultb(2)
! 
!       call pent_gauge(pentmom(0,2),pentmom(0,4),pentmom(0,5),pentmom(0,3),pentmom(0,1),
!      &      psi(1,-1,2),psi(1,-1,1),photpol(0,2),wpol,glupol,
!      &      -1,mur,resultgauge,resultgaugeb,wardresult)
!       print*, result(2),resultb(2)

      CALL penlineABEmmm(pentmom(0,2),pentmom(0,4),pentmom(0,5),pentmom(0,3),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),photpoltilde(0,2),wpoltilde,glupoltilde,
     &      -1,mur,1,4,4,1,resultgauge(1),penresult,resultgaugeb(1),resultb(2))

      CALL pent_to_boxnoab(pentmom(0,2),pentmom(0,4),pentmom(0,5),pentmom(0,3),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),photpol(0,2),wpol,glupol,xp2,xw,xg,-1,mur,3,testpb1,testpb2)

      call pent_gauge(pentmom(0,2),pentmom(0,4),pentmom(0,5),pentmom(0,3),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),photpoltilde(0,2),wpoltilde,glupoltilde,
     &      -1,mur,resultgauge,resultgaugeb,wardresult)

      if (wardresult) then
      result(2)=(cfc*penresult(1)+cfcac*penresult(2)+testpb1)*gad(1)
      resultb(2)=(resultb(2)+testpb2)*gad(1)  ! gives +1*Madgraph
      else
      result(2)=(testpb1)*gad(1)
      resultb(2)=(resultb(2)+testpb2)*gad(1)  ! gives +1*Madgraph
      endif

      if (ldebug) then
!       if (abs(testpb1).eq.0d0) then
!       ubg(incross)=.false.
!       endif
!       print*, "penline gauge check photon-W-gluon", wardresult
      if (.not.wardresult) then
!       print*, "bad point pent 2"
      badpointspent(incross+2)=badpointspent(incross+2)+1 
!      sumresult=0d0
!      return
      ubg(incross)=.false.
      endif
      endif

!       CALL  penlineABEmmm(pentmom(0,2),pentmom(0,5),pentmom(0,3),pentmom(0,4),pentmom(0,1),
!      &      psi(1,-1,2),psi(1,-1,1),wpol,glupol,photpol(0,1),
!      &      -1,mur,1,3,4,1,resultgauge(1),penresult,resultgaugeb(1),resultb(3))
! 
!       result(3)=(cfc*penresult(1)+cfcac*penresult(2))*gau(1)
!       resultb(3)=resultb(3)*gau(1)    ! gives +1*HELAS
!       test(3)=result(3)+cte*resultb(3)
! 
!       CALL pent_gauge(pentmom(0,2),pentmom(0,5),pentmom(0,3),pentmom(0,4),pentmom(0,1),
!      &      psi(1,-1,2),psi(1,-1,1),wpol,glupol,photpol(0,1),
!      &      -1,    mur,resultgauge,resultgaugeb,wardresult)
!       print*, result(3),resultb(3)

      CALL  penlineABEmmm(pentmom(0,2),pentmom(0,5),pentmom(0,3),pentmom(0,4),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),wpoltilde,glupoltilde,photpoltilde(0,1),
     &      -1,mur,1,3,4,1,resultgauge(1),penresult,resultgaugeb(1),resultb(3))

      CALL pent_to_boxnoab(pentmom(0,2),pentmom(0,5),pentmom(0,3),pentmom(0,4),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),wpol,glupol,photpol(0,1),xw,xg,xp1,-1,mur,6,testpb1,testpb2)

      CALL pent_gauge(pentmom(0,2),pentmom(0,5),pentmom(0,3),pentmom(0,4),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),wpoltilde,glupoltilde,photpoltilde(0,1),
     &      -1,    mur,resultgauge,resultgaugeb,wardresult)

      if (wardresult) then
      result(3)=(cfc*penresult(1)+cfcac*penresult(2)+testpb1)*gau(1)
      resultb(3)=(resultb(3)+testpb2)*gau(1)  ! gives +1*Madgraph
      else
      result(3)=(testpb1)*gau(1)
      resultb(3)=(resultb(3)+testpb2)*gau(1)  ! gives +1*Madgraph
      endif

      if (ldebug) then
!       if (abs(testpb1).eq.0d0) then
!       ubg(incross)=.false.
!       endif
!       print*, "penline gauge check W-gluon-photon", wardresult
      if (.not.wardresult) then
!       print*, "bad point pent 3"
      badpointspent(incross+2)=badpointspent(incross+2)+1 
!      sumresult=0d0
!      return
      ubg(incross)=.false.
      endif
      endif

!       CALL penlineABEmmm(pentmom(0,2),pentmom(0,4),pentmom(0,3),pentmom(0,5),pentmom(0,1),
!      &      psi(1,-1,2),psi(1,-1,1),photpol(0,2),glupol,wpol,
!      &      -1,mur,1,3,4,1,resultgauge(1),penresult,resultgaugeb(1),resultb(4))
! 
!       result(4)=(cfc*penresult(1)+cfcac*penresult(2))*gad(1)
!       resultb(4)=resultb(4)*gad(1)
!       test(4)=result(4)+cte*resultb(4)
! 
!       CALL pent_gauge(pentmom(0,2),pentmom(0,4),pentmom(0,3),pentmom(0,5),pentmom(0,1),
!      &      psi(1,-1,2),psi(1,-1,1),photpol(0,2),glupol,wpol,
!      &      -1,    mur ,resultgauge,resultgaugeb,wardresult)
!       print*, result(4),resultb(4)

      CALL penlineABEmmm(pentmom(0,2),pentmom(0,4),pentmom(0,3),pentmom(0,5),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),photpoltilde(0,2),glupoltilde,wpoltilde,
     &      -1,mur,1,3,4,1,resultgauge(1),penresult,resultgaugeb(1),resultb(4))

      CALL pent_to_boxnoab(pentmom(0,2),pentmom(0,4),pentmom(0,3),pentmom(0,5),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),photpol(0,2),glupol,wpol,xp2,xg,xw,-1,mur,6,testpb1,testpb2)

      CALL pent_gauge(pentmom(0,2),pentmom(0,4),pentmom(0,3),pentmom(0,5),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),photpoltilde(0,2),glupoltilde,wpoltilde,
     &      -1,    mur ,resultgauge,resultgaugeb,wardresult)

      if (wardresult) then
      result(4)=(cfc*penresult(1)+cfcac*penresult(2)+testpb1)*gad(1)
      resultb(4)=(resultb(4)+testpb2)*gad(1)  ! gives +1*Madgraph
      else
      result(4)=(testpb1)*gad(1)
      resultb(4)=(resultb(4)+testpb2)*gad(1)  ! gives +1*Madgraph
      endif

      if (ldebug) then
!       if (abs(testpb1).eq.0d0) then
!       ubg(incross)=.false.
!       endif
!       print*, "penline gauge check photon-gluon-W", wardresult
      if (.not.wardresult) then
!       print*, "bad point pent 4"
      badpointspent(incross+2)=badpointspent(incross+2)+1
!      sumresult=0d0
!      return 
      ubg(incross)=.false.
      endif
      endif

!       CALL penlineABEmmm(pentmom(0,2),pentmom(0,3),pentmom(0,5),pentmom(0,4),pentmom(0,1),
!      &      psi(1,-1,2),psi(1,-1,1),glupol,wpol,photpol(0,1),
!      &      -1,mur,1,2,4,1,resultgauge(1),penresult,resultgaugeb(1),resultb(5))
! 
!       result(5)=(cfc*penresult(1)+cfcac*penresult(2))*gau(1)
!       resultb(5)=resultb(5)*gau(1)
!       test(5)=result(5)+cte*resultb(5)
! 
!       call pent_gauge(pentmom(0,2),pentmom(0,3),pentmom(0,5),pentmom(0,4),pentmom(0,1),
!      &      psi(1,-1,2),psi(1,-1,1),glupol,wpol,photpol(0,1),
!      &      -1,    mur,resultgauge,resultgaugeb,wardresult)
!       print*, result(5),resultb(5)

      CALL penlineABEmmm(pentmom(0,2),pentmom(0,3),pentmom(0,5),pentmom(0,4),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),glupoltilde,wpoltilde,photpoltilde(0,1),
     &      -1,mur,1,2,4,1,resultgauge(1),penresult,resultgaugeb(1),resultb(5))

      CALL pent_to_boxnoab(pentmom(0,2),pentmom(0,3),pentmom(0,5),pentmom(0,4),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),glupol,wpol,photpol(0,1),xg,xw,xp1,-1,mur,2,testpb1,testpb2)

      call pent_gauge(pentmom(0,2),pentmom(0,3),pentmom(0,5),pentmom(0,4),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),glupoltilde,wpoltilde,photpoltilde(0,1),
     &      -1,    mur,resultgauge,resultgaugeb,wardresult)

      if (wardresult) then
      result(5)=(cfc*penresult(1)+cfcac*penresult(2)+testpb1)*gau(1)
      resultb(5)=(resultb(5)+testpb2)*gau(1)  ! gives +1*Madgraph
      else
      result(5)=(testpb1)*gau(1)
      resultb(5)=(resultb(5)+testpb2)*gau(1)  ! gives +1*Madgraph
      endif

      if (ldebug) then
!       if (abs(testpb1).eq.0d0) then
!       ubg(incross)=.false.
!       endif
!       print*, "penline gauge check photon-W-gluon", wardresult
      if (.not.wardresult) then
!       print*, "bad point pent 5"
      badpointspent(incross+2)=badpointspent(incross+2)+1 
!      sumresult=0d0
!      return
      ubg(incross)=.false.
      endif
      endif
!                             p1           p2           p3           p4           p5    
!       CALL  penlineABEmmm(pentmom(0,2),pentmom(0,5),pentmom(0,4),pentmom(0,3),pentmom(0,1),
!      &      psi(1,-1,2),psi(1,-1,1),wpol,photpol(0,1),glupol,
!      &      -1,mur,1,4,4,1,resultgauge(1),penresult,resultgaugeb(1),resultb(6))
! 
!       result(6)=(cfc*penresult(1)+cfcac*penresult(2))*gau(1)
!       resultb(6)=resultb(6)*gau(1)
!       test(6)=result(6)+cte*resultb(6)
! 
!       CALL pent_gauge(pentmom(0,2),pentmom(0,5),pentmom(0,4),pentmom(0,3),pentmom(0,1),
!      &      psi(1,-1,2),psi(1,-1,1),wpol,photpol(0,1),glupol,
!      &      -1,    mur,resultgauge,resultgaugeb,wardresult)
!       print*, result(6),resultb(6)

      CALL  penlineABEmmm(pentmom(0,2),pentmom(0,5),pentmom(0,4),pentmom(0,3),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),wpoltilde,photpoltilde(0,1),glupoltilde,
     &      -1,mur,1,4,4,1,resultgauge(1),penresult,resultgaugeb(1),resultb(6))

      CALL pent_to_boxnoab(pentmom(0,2),pentmom(0,5),pentmom(0,4),pentmom(0,3),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),wpol,photpol(0,1),glupol,xw,xp1,xg,-1,mur,3,testpb1,testpb2)

      CALL pent_gauge(pentmom(0,2),pentmom(0,5),pentmom(0,4),pentmom(0,3),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),wpoltilde,photpoltilde(0,1),glupoltilde,
     &      -1,    mur,resultgauge,resultgaugeb,wardresult)

      if (wardresult) then
      result(6)=(cfc*penresult(1)+cfcac*penresult(2)+testpb1)*gau(1)
      resultb(6)=(resultb(6)+testpb2)*gau(1)  ! gives +1*Madgraph
      else
      result(6)=(testpb1)*gau(1)
      resultb(6)=(resultb(6)+testpb2)*gau(1)  ! gives +1*Madgraph
      endif

      if (ldebug) then
!       if (abs(testpb1).eq.0d0) then
!       ubg(incross)=.false.
!       endif
!       print*, "penline gauge check W-photon-gluon", wardresult
      if (.not.wardresult) then
!       print*, "bad point pent 6"
      badpointspent(incross+2)=badpointspent(incross+2)+1 
!      sumresult=0d0
!      return
      ubg(incross)=.false.
      endif
      endif

!  sum the abelian pentagon contributions
      virt=virt+(result(1)+result(2)+result(3)+result(4)+result(5)+result(6))*gg(1)**3*gwf(1)
!  sum the born contributions
      ampborn=ampborn+(resultb(1)+resultb(2)+resultb(3)+resultb(4)+resultb(5)+resultb(6))*gg(1)*gwf(1)
!
!  non-abelian
!
      CALL  penlineNoAbemmm(pentmom(0,2),pentmom(0,4),pentmom(0,5),pentmom(0,1),pentmom(0,3),
     &      psi(1,-1,2),psi(1,-1,1),photpol(0,2),wpol,glupol,
     &      -1,    mur,     4      , 1 ,resultgauge(1),result(1))

      CALL pentNoAbe_gauge(accuracy,pentmom(0,2),pentmom(0,4),pentmom(0,5),pentmom(0,1),pentmom(0,3),
     &      psi(1,-1,2),psi(1,-1,1),photpol(0,2),wpol,glupol,
     &      -1,    mur, resultgauge(1), accuracyres,Wardresult)

      if (ldebug) then
      if (.not.wardresult) then
      badpointspentnoab(incross+2)=badpointspentnoab(incross+2)+1 
!      sumresult=0d0
!      return
      ubg(incross)=.false.
      endif
      endif

      result(1)=cacgraph*result(1)*gad(1)
      test(7)=result(1)

      CALL  penlineNoAbemmm(pentmom(0,2),pentmom(0,5),pentmom(0,4),pentmom(0,1),pentmom(0,3),
     &      psi(1,-1,2),psi(1,-1,1),wpol,photpol(0,1),glupol,
     &      -1,    mur,     4      , 1 ,resultgauge(1),result(2))

      CALL pentNoAbe_gauge(accuracy,pentmom(0,2),pentmom(0,5),pentmom(0,4),pentmom(0,1),pentmom(0,3),
     &      psi(1,-1,2),psi(1,-1,1),wpol,photpol(0,1),glupol,
     &      -1,    mur, resultgauge(1), accuracyres,Wardresult)

      if (ldebug) then
      if (.not.wardresult) then
      badpointspentnoab(incross+2)=badpointspentnoab(incross+2)+1 
!      sumresult=0d0
!      return
      ubg(incross)=.false.
      endif
      endif


      result(2)=cacgraph*result(2)*gau(1)
      test(8)=result(2)

!  sum the non-abelian pentagon contributions
      virt=virt+(result(1)+result(2))*gg(1)**3*gwf(1)
!
!   Box Contribution
!
      do mu=0,3
      boxmom(mu,1)=pentmom(mu,1)   !ub incoming
      boxmom(mu,2)=pentmom(mu,2)   !db outgoing
      boxmom(mu,3)=glumom(mu)  !g incoming
      boxmom(mu,4)=-1d0*wemnugamom(mu) !Wgamma outgoing
      polgluontemp(mu)=boxmom(mu,3)
      enddo

!
!   abelian
!    note: the second call is needed for the gauge test of the non-abelian boxes
      CALL boxlineABEmm(boxmom(0,2),boxmom(0,3),boxmom(0,4),boxmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),glupol,currlepwga(0,photpolin),
     &      -1, mur, 1, 2 ,3,1,resultgauge(1),boxresult,
     &      resultgaugeb(1),resultb(1))
      CALL boxlineABEmm(boxmom(0,2),boxmom(0,3),boxmom(0,4),boxmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),polgluontemp,currlepwga(0,photpolin),
     &      -1, mur, 1, 2 ,3,-1,nothing,gaugeCA1,
     &      nothing,nothing(1))


      result(1)=-1d0*(cfc*boxresult(1)+cfcac*boxresult(2))
      resultb(1)=-1d0*resultb(1) ! boxlineABEmm gives -1*HELAS for Born
      test(9)=result(1)+cte*resultb(1)

!
!   gauge test for boxes
!
      if (ldebug) then

      if (abs(resultgaugeb(1)).le.1e-10) then
      testgauge1=resultgauge(1)+cte*resultgaugeb(1)
      else
      testgauge1=resultgauge(1)/(-cte*resultgaugeb(1))-1.
      endif
      if (abs(resultgaugeb(2)).le.1e-10) then
      testgauge2=resultgauge(2)+cte*resultgaugeb(2)
      else
      testgauge2=resultgauge(2)/(-cte*resultgaugeb(2))-1.
      endif

      if (abs(testgauge1).ge.1e-2) then
!      print*, "here ub g box 1, g1", resultgauge(1)/(-cte*resultgaugeb(1)), resultgauge(1),(resultgaugeb(1))
      badpointsbox(incross+2)=badpointsbox(incross+2)+1
!      sumresult=0d0
!       STOP
!      return
      ubg(incross)=.false.
      endif
      if (abs(testgauge2).ge.1e-2) then
!      print*, "here ub g box 1, g2", resultgauge(1),resultgaugeb(1)
      badpointsbox(incross+2)=badpointsbox(incross+2)+1
!      sumresult=0d0
!       STOP
!      return
      ubg(incross)=.false.
      endif
      endif

!    note: the second call is needed for the gauge test of the non-abelian boxes
      CALL boxlineABEmm(boxmom(0,2),boxmom(0,4),boxmom(0,3),boxmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),currlepwga(0,photpolin),glupol,
     &      -1, mur, 1, 3 ,3,1,resultgauge2(1),boxresult,
     &      resultgaugeb(1),resultb(2))
      CALL boxlineABEmm(boxmom(0,2),boxmom(0,4),boxmom(0,3),boxmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),currlepwga(0,photpolin),polgluontemp,
     &      -1, mur, 1, 3 ,3,-1,nothing,gaugeCA2,
     &      nothing,nothing(1))
      
      result(2)=-1d0*(cfc*boxresult(1)+cfcac*boxresult(2))
      resultb(2)=-1d0*resultb(2) ! boxlineABEmm gives -1*HELAS for Born
      test(10)=result(2)+cte*resultb(2)
!
!   gauge test for boxes
!
      if (ldebug) then

      if (abs(resultgaugeb(1)).le.1e-10) then
      testgauge1=resultgauge2(1)+cte*resultgaugeb(1)
      else
      testgauge1=resultgauge2(1)/(-cte*resultgaugeb(1))-1.
      endif
      if (abs(resultgaugeb(2)).le.1e-10) then
      testgauge2=resultgauge2(2)+cte*resultgaugeb(2)
      else
      testgauge2=resultgauge2(2)/(-cte*resultgaugeb(2))-1.
      endif

      if (abs(testgauge1).ge.1e-2) then
!      print*, "here ub g box 2, g1"
      badpointsbox(incross+2)=badpointsbox(incross+2)+1
!      sumresult=0d0
!       STOP
!      return
      ubg(incross)=.false.
      endif
      if (abs(testgauge2).ge.1e-2) then
!      print*, "here ub g box 2, g2"
      badpointsbox(incross+2)=badpointsbox(incross+2)+1
!      sumresult=0d0
!       STOP
!      return
      ubg(incross)=.false.
      endif
      endif

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Non-abelian Part
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C W
      CALL boxlineNoAbemm(boxmom(0,2),boxmom(0,4),boxmom(0,1),boxmom(0,3),
     &      psi(1,-1,2),psi(1,-1,1),currlepwga(0,photpolin),glupol,
     &      -1, mur, 3, 1,resultgauge3(1),result(3))

** gauge test:
      if (abs(resultgauge3(2)).lt.1d-7) then
           testgauge2=abs(-(gaugeCA1(2)+ gaugeCA2(2))/2d0- resultgauge3(2) )
      else
           testgauge2=abs(-(gaugeCA1(2)+ gaugeCA2(2))/2d0/resultgauge3(2)    +1d0)
      endif   
      if (abs(testgauge2).ge.1d-1) then
         badpointsbox(incross+2)=badpointsbox(incross+2)+1
         ubg(incross)=.false.
      endif

      result(3)=-1d0*result(3)*cacgraph
      test(11)=result(3)

      virt=virt+(result(1)+result(2)+result(3))*gwf(1)*gg(1)**3
      ampborn=ampborn+(resultb(1)+resultb(2))*gwf(1)*gg(1)

      ancont=cac/2.*pic2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            Finite Terms from Pacos expansion     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      s=2.*dotrr(pentmom(0,1),pentmom(0,2))
!      t=2.*dotrr(pentmom(0,2),pentmom(0,3))
!      u=2.*dotrr(pentmom(0,1),pentmom(0,3))

      s=2.*dotrr(pentmom(0,1),pentmom(0,2))
      t=2.*dotrr(pentmom(0,1),pentmom(0,3))
      u=2.*dotrr(pentmom(0,2),pentmom(0,3))


       factfigy=-2*EulerGamma*Inv3*nf*TF + 2*Inv3*nf*TF*Log(4*Pi) + 
     -     CF*(-1 + 3*EulerGamma - EulerGamma**2 - Inv2*Inv3*Pi**2 + 
     -     Log(musq)*(2*EulerGamma - 2*Log(4*Pi)) - 
     -     Log(4*Pi)**2 - 2*EulerGamma*Log(Abs(s)) + (0,2)*EulerGamma*
     -     Pi*theta(s) + Log(4*Pi)*(-3 + 2*EulerGamma + 2*Log(Abs(s)) - 
     -     (0,2)*Pi*theta(s))) + 
     -  CA*((11 - 3*EulerGamma)*EulerGamma*Inv2*Inv3 - Inv2**2*Inv3*
     -     Pi**2 + Log(musq)*(EulerGamma - Log(4*Pi)) - Inv2*
     -     Log(4*Pi)**2 + EulerGamma*Log(Abs(s)) - 
     -     EulerGamma*Log(Abs(t)) - EulerGamma*Log(Abs(u)) + (0,1)*
     -     EulerGamma*Pi*(-theta(s) + theta(t) + theta(u)) + 
     -     Log(4*Pi)*(EulerGamma - 11*Inv2*Inv3 - Log(Abs(s)) + 
     -     Log(Abs(t)) + Log(Abs(u)) - (0,1)*Pi*(-theta(s) + theta(t) + 
     -     theta(u))))


      virt=virt
     &     +factfigy*ampborn*gg(1)**2

      ampbornubg(incross,glupolin,photpolin)=bornubg(incross,glupolin,
     -     photpolin)
      ampborn=bornubg(incross,glupolin,photpolin)

      m2s(incross+2)= m2s(incross+2)+2d0*dreal(virt*conjg(ampborn))/
     -     (4d0*pi)**2/24d0
      m2sborn(incross+2)=m2sborn(incross+2)+abs(ampborn)**2/24d0 
      m2stest(incross+2)=m2stest(incross+2)+2.*dreal(((test(1)+test(2)+
     -     test(3)+test(4)+test(5)
     &    +test(6)+test(7)+test(8)+test(9)+test(10)+test(11))*gg(1)**3*
     -     gwf(1))*conjg(ampborn))/9.

      enddo  !loop over initial state exchanged


      enddo  !loop over gluon polarization
      enddo  !loop over photon polarization

      m2s(3)=m2s(3)*(pdf(-2,1)*pdf(0,2)+pdf(-4,1)*pdf(0,2))
      m2s(4)=m2s(4)*(pdf(-2,2)*pdf(0,1)+pdf(-4,2)*pdf(0,1))
      m2stest(3)=m2stest(3)*(pdf(-2,1)*pdf(0,2)+pdf(-4,1)*pdf(0,2))
      m2stest(4)=m2stest(4)*(pdf(-2,2)*pdf(0,1)+pdf(-4,1)*pdf(0,1))
      m2sbornwopdf(3)=m2sborn(3)
      m2sbornwopdf(4)=m2sborn(4)
      m2sborn(3)=m2sborn(3)*(pdf(-2,1)*pdf(0,2)+pdf(-4,1)*pdf(0,2))
      m2sborn(4)=m2sborn(4)*(pdf(-2,2)*pdf(0,1)+pdf(-4,2)*pdf(0,1))

      if (borndebug) then
      do k=0,3
         helastestmom(k,1)=p(k,1,1)
         helastestmom(k,2)=p(k,2,1)
         helastestmom(k,3)=p(k,3,1)
         helastestmom(k,4)=v(k,1,1)
         helastestmom(k,5)=v(k,2,1)
         helastestmom(k,6)=v(k,3,1)
         helastestmom(k,7)=v(k,4,1)
      enddo

! u~ g -> d a e- ve~ 

      Call SUBGDBEMVEBMUMMUP(helastestmom,testh)
      testh=testh*(pdf(-2,1)*pdf(0,2)+pdf(-4,1)*pdf(0,2))

      print*,  "born from virtuals/HELAS = ",m2sborn(3)/testh, "  =1 ?"

      do k=0,3
         helastestmom(k,1)=p(k,2,1)
         helastestmom(k,2)=p(k,1,1)
         helastestmom(k,3)=p(k,3,1)
         helastestmom(k,4)=v(k,1,1)
         helastestmom(k,5)=v(k,2,1)
         helastestmom(k,6)=v(k,3,1)
         helastestmom(k,7)=v(k,4,1)
      enddo
      Call SUBGDBEMVEBMUMMUP(helastestmom,testh)
      testh=testh*(pdf(-2,2)*pdf(0,1)+pdf(-4,2)*pdf(0,1))

      print*,  "born from virtuals/HELAS = ",m2sborn(4)/testh, "  =1 ? "

      do mu=0,3
      boxmomtilde(mu,1)=boxmom(mu,1) !INCOMING
      boxmomtilde(mu,2)=-boxmom(mu,2) !outgoing
      boxmomtilde(mu,3)=-boxmom(mu,3) !outgoing
      boxmomtilde(mu,4)=boxmom(mu,4) !outgoing
      enddo

      print*,  mur
      print*, "-s23",-2*dotrr(boxmomtilde(0,2),boxmomtilde(0,3)) !
      print*, "-t", 2*dotrr(boxmomtilde(0,1),boxmomtilde(0,3))   !  additional minus due to incoming conv
      print*, "-t12",2*dotrr(boxmomtilde(0,1),boxmomtilde(0,2))  !
      print*, "finished ub g process"


      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   d g -> w gamma u
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!       cte=(-1.0d0*pic2/3.0d0+7.0d0)
!
!     initial state exchanged
!
      do incross=1,2

      do mu=0,3
      glumom(mu)=p(mu,mapspin2(incross),1)
      enddo

      do photpolin=-1,1,2
      do glupolin=-1,1,2
!
!  Polarizations using HELAS
!
!   photon polarization
!
      do mu=0,3
      photpol(mu,1)=storephotpol(mu,photpolin,1)
      photpol(mu,2)=storephotpol(mu,photpolin,2)
      currlepwga(mu,photpolin)=storewphotpol(mu,photpolin)
      wpol(mu)=storewpol(mu)
      enddo
!
!   gluon polarization
!
      CALL VXXXXX(glumom,0d0 ,glupolin,gluin, temp)
      do mu=0,3
      glupol(mu)=temp(mu+1)
      enddo

      ampborn=dcmplx(0d0,0d0)
      ampborndg(incross,glupolin,photpolin)=dcmplx(0d0,0d0)
      virt=dcmplx(0d0,0d0)
      fact=dcmplx(0d0,0d0)
      do mu=1,6
      result(mu)=dcmplx(0d0,0d0)
      resultb(mu)=dcmplx(0d0,0d0)
      enddo 
      do mu=1,8
      penresult(mu)=dcmplx(0d0,0d0)
      enddo
      do mu=1,5
      boxresult(mu)=dcmplx(0d0,0d0)
      enddo
!
!     SPINORS
!
      do mu=0,3
      spinormom(mu,1)= p(mu,mapspin1(incross),1)  ! d to become ket
      spinormom(mu,2)= p(mu,3,1) ! u to become bra
      enddo
      signspinor(1) =1
      signspinor(2) =1
      CALL PSI0M(2,spinormom(0,1),signspinor,psi)

!
!   Pentagon Contribution
!
      do mu=0,3
      pentmom(mu,1)=-1d0*p(mu,3,1)   !u outgoing
      pentmom(mu,2)=1d0*p(mu,mapspin1(incross),1)  !d incoming
      pentmom(mu,3)=1d0*glumom(mu)   !g incoming
      pentmom(mu,4)=-1d0*photmom(mu) !gamma outgoing
      pentmom(mu,5)=-1d0*wmom(mu)    !w- outgoing
      enddo

      do mu=0,3
      tmp(mu)=wmom(mu)+photmom(mu)
      enddo
      xp2=dotcr(photpol(0,2),tmp)/dcmplx(dotrr(photmom,tmp))
      xp1=dotcr(photpol(0,2),tmp)/dcmplx(dotrr(photmom,tmp))
      xw=dotcr(wpol,tmp)/dcmplx(dotrr(wmom,tmp))
      xg=0d0!dotcr(glupol,tmp)/dcmplx(dotrr(glumom,tmp))
      do mu=0,3
         photpoltilde(mu,2)=photpol(mu,2)-xp2*dcmplx(photmom(mu),0d0)
         photpoltilde(mu,1)=photpol(mu,1)-xp1*dcmplx(photmom(mu),0d0)
         wpoltilde(mu)   =wpol(mu)-xw*dcmplx(wmom(mu),0d0)
         glupoltilde(mu) =glupol(mu)!-xg*dcmplx(glumom(mu),0d0)
      enddo

      s15=2d0*dotrr(pentmom(0,1),pentmom(0,2))
      cte= -(-3d0*(Log(mur) - Log(abs(s15)) + (0d0,1d0)*pi*theta(s15)) 
     &  - (Log(mur) - Log(abs(s15)) + (0d0,1d0)*pi*theta(s15))**2
c   substracted the c_virt
     &  -7d0 + 1d0/3d0*pi**2)
!
!   abelian
!
!       CALL  penlineABEmmm(pentmom(0,2),pentmom(0,3),pentmom(0,4),pentmom(0,5),pentmom(0,1),
!      &      psi(1,-1,2),psi(1,-1,1),glupol,photpol(0,2),wpol,
!      &      -1,    mur,    1,        2  ,    4      , 1 ,resultgauge(1),
!      &      penresult,resultgaugeb(1),resultb(1))
! 
!       result(1)=(cfc*penresult(1)+cfcac*penresult(2))*gad(1)
!       resultb(1)=resultb(1)*gad(1)  ! gives -1*Madgraph?????
!       test(1)=result(1)+cte*resultb(1)
! 
!       call pent_gauge(pentmom(0,2),pentmom(0,3),pentmom(0,4),pentmom(0,5),pentmom(0,1),
!      &      psi(1,-1,2),psi(1,-1,1),glupol,photpol(0,2),wpol,
!      &      -1,    mur ,resultgauge,resultgaugeb,wardresult)
!       print*, result(1),resultb(1)

      CALL  penlineABEmmm(pentmom(0,2),pentmom(0,3),pentmom(0,4),
     -     pentmom(0,5),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),glupoltilde,photpoltilde(0,2),wpoltilde,
     &      -1,    mur,    1,        2  ,    4      , 1 ,resultgauge(1),
     &      penresult,resultgaugeb(1),resultb(1))

      call pent_to_boxnoab(pentmom(0,2),pentmom(0,3),pentmom(0,4),
     -     pentmom(0,5),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),glupol,photpol(0,2),wpol,xg,xp2,
     -     xw,-1,mur,2,testpb1,testpb2)

      call pent_gauge(pentmom(0,2),pentmom(0,3),pentmom(0,4),
     -     pentmom(0,5),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),glupoltilde,photpoltilde(0,2),wpoltilde,
     &      -1,    mur ,resultgauge,resultgaugeb,wardresult)

      if (wardresult) then
      result(1)=(cfc*penresult(1)+cfcac*penresult(2)+testpb1)*gad(1)
      resultb(1)=(resultb(1)+testpb2)*gad(1)  ! gives +1*Madgraph
      else
      result(1)=(testpb1)*gad(1)
      resultb(1)=(resultb(1)+testpb2)*gad(1)  ! gives +1*Madgraph
      endif

      if (ldebug) then
!       if (abs(testpb1).eq.0d0) then
!       dg(incross)=.false.
!       endif
!       print*, "penline gauge check gluon-photon-W", wardresult
      if (.not.wardresult) then
!       print*, "bad point pent 1"
      badpointspent(incross+4)=badpointspent(incross+4)+1
!      sumresult=0d0
!       STOP
!      return
      dg(incross)=.false.
      endif
      endif

!       CALL  penlineABEmmm(pentmom(0,2),pentmom(0,4),pentmom(0,5),pentmom(0,3),pentmom(0,1),
!      &      psi(1,-1,2),psi(1,-1,1),photpol(0,2),wpol,glupol,
!      &      -1,mur,1,4,4,1,resultgauge(1),penresult,resultgaugeb(1),resultb(2))
! 
!       result(2)=(cfc*penresult(1)+cfcac*penresult(2))*gad(1)
!       resultb(2)=resultb(2)*gad(1)
!       test(2)=result(2)+cte*resultb(2)
! 
!       call pent_gauge(pentmom(0,2),pentmom(0,4),pentmom(0,5),pentmom(0,3),pentmom(0,1),
!      &      psi(1,-1,2),psi(1,-1,1),photpol(0,2),wpol,glupol,
!      &      -1,mur,resultgauge,resultgaugeb,wardresult)
!       print*, result(2),resultb(2)

      CALL penlineABEmmm(pentmom(0,2),pentmom(0,4),pentmom(0,5),
     -     pentmom(0,3),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),photpoltilde(0,2),wpoltilde,glupoltilde,
     &      -1,mur,1,4,4,1,resultgauge(1),penresult,resultgaugeb(1),resultb(2))

      CALL pent_to_boxnoab(pentmom(0,2),pentmom(0,4),pentmom(0,5),
     -     pentmom(0,3),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),photpol(0,2),wpol,glupol,xp2,xw,
     -     xg,-1,mur,3,testpb1,testpb2)

      call pent_gauge(pentmom(0,2),pentmom(0,4),pentmom(0,5),
     -     pentmom(0,3),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),photpoltilde(0,2),wpoltilde,glupoltilde,
     &      -1,mur,resultgauge,resultgaugeb,wardresult)

      if (wardresult) then
      result(2)=(cfc*penresult(1)+cfcac*penresult(2)+testpb1)*gad(1)
      resultb(2)=(resultb(2)+testpb2)*gad(1)  ! gives +1*Madgraph
      else
      result(2)=(testpb1)*gad(1)
      resultb(2)=(resultb(2)+testpb2)*gad(1)  ! gives +1*Madgraph
      endif

      if (ldebug) then
!       if (abs(testpb1).eq.0d0) then
!       dg(incross)=.false.
!       endif
!       print*, "penline gauge check photon-W-gluon", wardresult
      if (.not.wardresult) then
!       print*, "bad point pent 2"
      badpointspent(incross+4)=badpointspent(incross+4)+1 
!      sumresult=0d0
!       STOP
!      return
      dg(incross)=.false.
      endif
      endif

!       CALL  penlineABEmmm(pentmom(0,2),pentmom(0,5),pentmom(0,3),pentmom(0,4),pentmom(0,1),
!      &      psi(1,-1,2),psi(1,-1,1),wpol,glupol,photpol(0,1),
!      &      -1,mur,1,3,4,1,resultgauge(1),penresult,resultgaugeb(1),resultb(3))
! 
!       result(3)=(cfc*penresult(1)+cfcac*penresult(2))*gau(1)
!       resultb(3)=resultb(3)*gau(1)    ! gives +1*HELAS
!       test(3)=result(3)+cte*resultb(3)
! 
!       CALL pent_gauge(pentmom(0,2),pentmom(0,5),pentmom(0,3),pentmom(0,4),pentmom(0,1),
!      &      psi(1,-1,2),psi(1,-1,1),wpol,glupol,photpol(0,1),
!      &      -1,    mur,resultgauge,resultgaugeb,wardresult)
!       print*, result(3),resultb(3)

      CALL  penlineABEmmm(pentmom(0,2),pentmom(0,5),pentmom(0,3),
     -     pentmom(0,4),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),wpoltilde,glupoltilde,photpoltilde(0,1),
     &      -1,mur,1,3,4,1,resultgauge(1),penresult,resultgaugeb(1),
     -     resultb(3))

      CALL pent_to_boxnoab(pentmom(0,2),pentmom(0,5),pentmom(0,3),
     -     pentmom(0,4),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),wpol,glupol,photpol(0,1),xw,xg,xp1,-1,mur,6,testpb1,testpb2)

      CALL pent_gauge(pentmom(0,2),pentmom(0,5),pentmom(0,3),
     -     pentmom(0,4),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),wpoltilde,glupoltilde,photpoltilde(0,1),
     &      -1,    mur,resultgauge,resultgaugeb,wardresult)

      if (wardresult) then
      result(3)=(cfc*penresult(1)+cfcac*penresult(2)+testpb1)*gau(1)
      resultb(3)=(resultb(3)+testpb2)*gau(1)  ! gives +1*Madgraph
      else
      result(3)=(testpb1)*gau(1)
      resultb(3)=(resultb(3)+testpb2)*gau(1)  ! gives +1*Madgraph
      endif

      if (ldebug) then
!       if (abs(testpb1).eq.0d0) then
!       dg(incross)=.false.
!       endif
!       print*, "penline gauge check W-gluon-photon", wardresult
      if (.not.wardresult) then
!       print*, "bad point pent 3"
      badpointspent(incross+4)=badpointspent(incross+4)+1 
!      sumresult=0d0
!       STOP
!      return
      dg(incross)=.false.
      endif
      endif

!       CALL penlineABEmmm(pentmom(0,2),pentmom(0,4),pentmom(0,3),pentmom(0,5),pentmom(0,1),
!      &      psi(1,-1,2),psi(1,-1,1),photpol(0,2),glupol,wpol,
!      &      -1,mur,1,3,4,1,resultgauge(1),penresult,resultgaugeb(1),resultb(4))
! 
!       result(4)=(cfc*penresult(1)+cfcac*penresult(2))*gad(1)
!       resultb(4)=resultb(4)*gad(1)
!       test(4)=result(4)+cte*resultb(4)
! 
!       CALL pent_gauge(pentmom(0,2),pentmom(0,4),pentmom(0,3),pentmom(0,5),pentmom(0,1),
!      &      psi(1,-1,2),psi(1,-1,1),photpol(0,2),glupol,wpol,
!      &      -1,    mur ,resultgauge,resultgaugeb,wardresult)
!       print*, result(4),resultb(4)

      CALL penlineABEmmm(pentmom(0,2),pentmom(0,4),pentmom(0,3),
     -     pentmom(0,5),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),photpoltilde(0,2),glupoltilde,wpoltilde,
     &      -1,mur,1,3,4,1,resultgauge(1),penresult,resultgaugeb(1),
     -     resultb(4))

      CALL pent_to_boxnoab(pentmom(0,2),pentmom(0,4),pentmom(0,3),pentmom(0,5),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),photpol(0,2),glupol,wpol,xp2,xg,xw,-1,mur,6,testpb1,testpb2)

      CALL pent_gauge(pentmom(0,2),pentmom(0,4),pentmom(0,3),pentmom(0,5),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),photpoltilde(0,2),glupoltilde,wpoltilde,
     &      -1,    mur ,resultgauge,resultgaugeb,wardresult)

      if (wardresult) then
      result(4)=(cfc*penresult(1)+cfcac*penresult(2)+testpb1)*gad(1)
      resultb(4)=(resultb(4)+testpb2)*gad(1)  ! gives +1*Madgraph
      else
      result(4)=(testpb1)*gad(1)
      resultb(4)=(resultb(4)+testpb2)*gad(1)  ! gives +1*Madgraph
      endif

      if (ldebug) then
!       if (abs(testpb1).eq.0d0) then
!       dg(incross)=.false.
!       endif
!       print*, "penline gauge check photon-gluon-W", wardresult
      if (.not.wardresult) then
!       print*, "bad point pent 4"
      badpointspent(incross+4)=badpointspent(incross+4)+1
!      sumresult=0d0
!       STOP
!      return 
      dg(incross)=.false.
      endif
      endif

!       CALL penlineABEmmm(pentmom(0,2),pentmom(0,3),pentmom(0,5),pentmom(0,4),pentmom(0,1),
!      &      psi(1,-1,2),psi(1,-1,1),glupol,wpol,photpol(0,1),
!      &      -1,mur,1,2,4,1,resultgauge(1),penresult,resultgaugeb(1),resultb(5))
! 
!       result(5)=(cfc*penresult(1)+cfcac*penresult(2))*gau(1)
!       resultb(5)=resultb(5)*gau(1)
!       test(5)=result(5)+cte*resultb(5)
! 
!       call pent_gauge(pentmom(0,2),pentmom(0,3),pentmom(0,5),pentmom(0,4),pentmom(0,1),
!      &      psi(1,-1,2),psi(1,-1,1),glupol,wpol,photpol(0,1),
!      &      -1,    mur,resultgauge,resultgaugeb,wardresult)
!       print*, result(5), resultb(5)

      CALL penlineABEmmm(pentmom(0,2),pentmom(0,3),pentmom(0,5),pentmom(0,4),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),glupoltilde,wpoltilde,photpoltilde(0,1),
     &      -1,mur,1,2,4,1,resultgauge(1),penresult,resultgaugeb(1),resultb(5))

      CALL pent_to_boxnoab(pentmom(0,2),pentmom(0,3),pentmom(0,5),pentmom(0,4),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),glupol,wpol,photpol(0,1),xg,xw,xp1,-1,mur,2,testpb1,testpb2)

      call pent_gauge(pentmom(0,2),pentmom(0,3),pentmom(0,5),pentmom(0,4),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),glupoltilde,wpoltilde,photpoltilde(0,1),
     &      -1,    mur,resultgauge,resultgaugeb,wardresult)
 
      if (wardresult) then
      result(5)=(cfc*penresult(1)+cfcac*penresult(2)+testpb1)*gau(1)
      resultb(5)=(resultb(5)+testpb2)*gau(1)  ! gives +1*Madgraph
      else
      result(5)=(testpb1)*gau(1)
      resultb(5)=(resultb(5)+testpb2)*gau(1)  ! gives +1*Madgraph
      endif

      if (ldebug) then
!       if (abs(testpb1).eq.0d0) then
!       dg(incross)=.false.
!       endif
!       print*, "penline gauge check photon-W-gluon", wardresult
      if (.not.wardresult) then
!       print*, "bad point pent 5"
      badpointspent(incross+4)=badpointspent(incross+4)+1 
!      sumresult=0d0
!       STOP
!      return
      dg(incross)=.false.
      endif
      endif

!                             p1           p2           p3           p4           p5    
!      CALL  penlineABEmmm(pentmom(0,2),pentmom(0,5),pentmom(0,4),pentmom(0,3),pentmom(0,1),
!     &      psi(1,-1,2),psi(1,-1,1),wpol,photpol(0,1),glupol,
!     &      -1,mur,1,4,4,1,resultgauge(1),penresult,resultgaugeb(1),resultb(6))
!
!      result(6)=(cfc*penresult(1)+cfcac*penresult(2))*gau(1)
!      resultb(6)=resultb(6)*gau(1)
!      test(6)=result(6)+cte*resultb(6)
!
!      CALL pent_gauge(pentmom(0,2),pentmom(0,5),pentmom(0,4),pentmom(0,3),pentmom(0,1),
!     &      psi(1,-1,2),psi(1,-1,1),wpol,photpol(0,1),glupol,
!     &      -1,    mur,resultgauge,resultgaugeb,wardresult)

      CALL  penlineABEmmm(pentmom(0,2),pentmom(0,5),pentmom(0,4),pentmom(0,3),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),wpoltilde,photpoltilde(0,1),glupoltilde,
     &      -1,mur,1,4,4,1,resultgauge(1),penresult,resultgaugeb(1),resultb(6))

      CALL pent_to_boxnoab(pentmom(0,2),pentmom(0,5),pentmom(0,4),pentmom(0,3),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),wpol,photpol(0,1),glupol,xw,xp1,xg,-1,mur,3,testpb1,testpb2)

      CALL pent_gauge(pentmom(0,2),pentmom(0,5),pentmom(0,4),pentmom(0,3),pentmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),wpoltilde,photpoltilde(0,1),glupoltilde,
     &      -1,    mur,resultgauge,resultgaugeb,wardresult)

      if (wardresult) then
      result(6)=(cfc*penresult(1)+cfcac*penresult(2)+testpb1)*gau(1)
      resultb(6)=(resultb(6)+testpb2)*gau(1)  ! gives +1*Madgraph
      else
      result(6)=(testpb1)*gau(1)
      resultb(6)=(resultb(6)+testpb2)*gau(1)  ! gives +1*Madgraph
      endif

      if (ldebug) then
!       if (abs(testpb1).eq.0d0) then
!       dg(incross)=.false.
!       endif
!       print*, "penline gauge check W-photon-gluon", wardresult
      if (.not.wardresult) then
!       print*, "bad point pent 6"
      badpointspent(incross+4)=badpointspent(incross+4)+1 
!      sumresult=0d0
!       STOP
!      return
      dg(incross)=.false.
      endif
      endif

!  sum the abelian pentagon contributions
      virt=virt+(result(1)+result(2)+result(3)+result(4)+result(5)+result(6))*gg(1)**3*gwf(1)
!  sum the born contributions
      ampborn=ampborn+(resultb(1)+resultb(2)+resultb(3)+resultb(4)+resultb(5)+resultb(6))*gg(1)*gwf(1)
!
!  non-abelian
!
      CALL  penlineNoAbemmm(pentmom(0,2),pentmom(0,4),pentmom(0,5),pentmom(0,1),pentmom(0,3),
     &      psi(1,-1,2),psi(1,-1,1),photpol(0,2),wpol,glupol,
     &      -1,    mur,     4      , 1 ,resultgauge(1),result(1))

      CALL pentNoAbe_gauge(accuracy,pentmom(0,2),pentmom(0,4),pentmom(0,5),pentmom(0,1),pentmom(0,3),
     &      psi(1,-1,2),psi(1,-1,1),photpol(0,2),wpol,glupol,
     &      -1,    mur, resultgauge(1), accuracyres,Wardresult)

      if (ldebug) then
      if (.not.wardresult) then
      badpointspentnoab(incross+4)=badpointspentnoab(incross+4)+1 
!      sumresult=0d0
!      return
      dg(incross)=.false.
      endif
      endif

      result(1)=cacgraph*result(1)*gad(1)
      test(7)=result(1)

      CALL  penlineNoAbemmm(pentmom(0,2),pentmom(0,5),pentmom(0,4),pentmom(0,1),pentmom(0,3),
     &      psi(1,-1,2),psi(1,-1,1),wpol,photpol(0,1),glupol,
     &      -1,    mur,     4      , 1 ,resultgauge(1),result(2))

      CALL pentNoAbe_gauge(accuracy,pentmom(0,2),pentmom(0,5),pentmom(0,4),pentmom(0,1),pentmom(0,3),
     &      psi(1,-1,2),psi(1,-1,1),wpol,photpol(0,1),glupol,
     &      -1,    mur, resultgauge(1), accuracyres,Wardresult)

      if (ldebug) then
      if (.not.wardresult) then
      badpointspentnoab(incross+4)=badpointspentnoab(incross+4)+1 
!      sumresult=0d0
!      return
      dg(incross)=.false.
      endif
      endif

      result(2)=cacgraph*result(2)*gau(1)
      test(8)=result(2)

!  sum the non-abelian pentagon contributions
      virt=virt+(result(1)+result(2))*gg(1)**3*gwf(1)
!
!   Box Contribution
!
      do mu=0,3
      boxmom(mu,1)=pentmom(mu,1)   !u outgoing
      boxmom(mu,2)=pentmom(mu,2)   !d incoming
      boxmom(mu,3)=glumom(mu)  !g incoming
      boxmom(mu,4)=-1d0*wemnugamom(mu) !Wgamma outgoing
      polgluontemp(mu)=boxmom(mu,3)
      enddo
!
!   abelian
!    note: the second call is needed for the gauge test of the non-abelian boxes
      CALL boxlineABEmm(boxmom(0,2),boxmom(0,3),boxmom(0,4),boxmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),glupol,currlepwga(0,photpolin),
     &      -1, mur, 1, 2 ,3,1,resultgauge(1),boxresult,
     &      resultgaugeb(1),resultb(1))
      CALL boxlineABEmm(boxmom(0,2),boxmom(0,3),boxmom(0,4),boxmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),polgluontemp,currlepwga(0,photpolin),
     &      -1, mur, 1, 2 ,3,-1,nothing,gaugeCA1,
     &      nothing,nothing(1))

      result(1)=-1d0*(cfc*boxresult(1)+cfcac*boxresult(2))
      resultb(1)=-1d0*resultb(1) ! boxlineABEmm gives -1*HELAS for Born
      test(9)=result(1)+cte*resultb(1)

!
!   gauge test for boxes
!
      if (ldebug) then

      if (abs(resultgaugeb(1)).le.1e-10) then
      testgauge1=resultgauge(1)+cte*resultgaugeb(1)
      else
      testgauge1=resultgauge(1)/(-cte*resultgaugeb(1))-1.
      endif
      if (abs(resultgaugeb(2)).le.1e-10) then
      testgauge2=resultgauge(2)+cte*resultgaugeb(2)
      else
      testgauge2=resultgauge(2)/(-cte*resultgaugeb(2))-1.
      endif

      if (abs(testgauge1).ge.1e-2) then
      badpointsbox(incross+4)=badpointsbox(incross+4)+1
!      print*, "here d g box 1, g1" , resultgauge(1)/resultgaugeb(1)
!      sumresult=0d0
!       STOP
!      return
      dg(incross)=.false.
      endif
      if (abs(testgauge2).ge.1e-2) then
!      print*, "here d g box 1, g2"
      badpointsbox(incross+4)=badpointsbox(incross+4)+1
!      sumresult=0d0
!       STOP
!      return
      dg(incross)=.false.
      endif
      endif

!    note: the second call is needed for the gauge test of the non-abelian boxes
      CALL boxlineABEmm(boxmom(0,2),boxmom(0,4),boxmom(0,3),boxmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),currlepwga(0,photpolin),glupol,
     &      -1, mur, 1, 3 ,3,1,resultgauge2(1),boxresult,
     &      resultgaugeb(1),resultb(2))      
      CALL boxlineABEmm(boxmom(0,2),boxmom(0,4),boxmom(0,3),boxmom(0,1),
     &      psi(1,-1,2),psi(1,-1,1),currlepwga(0,photpolin),polgluontemp,
     &      -1, mur, 1, 3 ,3,-1,nothing,gaugeCA2,
     &      nothing,nothing(1))
      
      result(2)=-1d0*(cfc*boxresult(1)+cfcac*boxresult(2))
      resultb(2)=-1d0*resultb(2) ! boxlineABEmm gives -1*HELAS for Born
      test(10)=result(2)+cte*resultb(2)
!
!   gauge test for boxes
!
      if (ldebug) then
      if (abs(resultgaugeb(1)).le.1e-10) then
      testgauge1=resultgauge2(1)+cte*resultgaugeb(1)
      else
      testgauge1=resultgauge2(1)/(-cte*resultgaugeb(1))-1.
      endif
      if (abs(resultgaugeb(2)).le.1e-10) then
      testgauge2=resultgauge2(2)+cte*resultgaugeb(2)
      else
      testgauge2=resultgauge2(2)/(-cte*resultgaugeb(2))-1.
      endif
      if (abs(testgauge1).ge.1e-2) then
!      print*, "here d g box 2, g1"
      badpointsbox(incross+4)=badpointsbox(incross+4)+1
!      sumresult=0d0
!       STOP
!      return
      dg(incross)=.false.
      endif
      if (abs(testgauge2).ge.1e-2) then
!      print*, "here  d g box 2, g2"
      badpointsbox(incross+4)=badpointsbox(incross+4)+1
!      sumresult=0d0
!       STOP
!      return
      dg(incross)=.false.
      endif
      endif

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Non-abelian Part
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C W
      CALL boxlineNoAbemm(boxmom(0,2),boxmom(0,4),boxmom(0,1),boxmom(0,3),
     &      psi(1,-1,2),psi(1,-1,1),currlepwga(0,photpolin),glupol,
     &      -1, mur, 3, 1,resultgauge3(1),result(3))

** gauge test:
      if (abs(resultgauge3(2)).lt.1d-7) then
           testgauge2=abs(-(gaugeCA1(2)+ gaugeCA2(2))/2d0- resultgauge3(2) )
      else
           testgauge2=abs(-(gaugeCA1(2)+ gaugeCA2(2))/2d0/resultgauge3(2)    +1d0)
      endif    
      if (abs(testgauge2).ge.1d-1) then
         badpointsbox(incross+4)=badpointsbox(incross+4)+1
         dg(incross)=.false.
      endif

      result(3)=-1d0*result(3)*cacgraph
      test(11)=result(3)

      virt=virt+(result(1)+result(2)+result(3))*gwf(1)*gg(1)**3
      ampborn=ampborn+(resultb(1)+resultb(2))*gwf(1)*gg(1)

      ancont=cac/2.*pic2
!       virt=virt+ancont*ampborn*gg(1)**2/(4.*pi)**2
!       virt=virt+(cte*cfc+ancont)*ampborn*gg(1)**2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            Finite Terms from Pacos expansion     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      s=2.*dotrr(pentmom(0,1),pentmom(0,2))
!      t=2.*dotrr(pentmom(0,2),pentmom(0,3))
!      u=2.*dotrr(pentmom(0,1),pentmom(0,3))


      s=2.*dotrr(pentmom(0,1),pentmom(0,2))
      t=2.*dotrr(pentmom(0,1),pentmom(0,3))
      u=2.*dotrr(pentmom(0,2),pentmom(0,3))


       factfigy=-2*EulerGamma*Inv3*nf*TF + 2*Inv3*nf*TF*Log(4*Pi) + 
     -     CF*(-1 + 3*EulerGamma - EulerGamma**2 - Inv2*Inv3*Pi**2 + 
     -     Log(musq)*(2*EulerGamma - 2*Log(4*Pi)) - 
     -     Log(4*Pi)**2 - 2*EulerGamma*Log(Abs(s)) + (0,2)*EulerGamma*
     -     Pi*theta(s) + Log(4*Pi)*(-3 + 2*EulerGamma + 2*Log(Abs(s)) - 
     -     (0,2)*Pi*theta(s))) + 
     -  CA*((11 - 3*EulerGamma)*EulerGamma*Inv2*Inv3 - Inv2**2*Inv3*
     -     Pi**2 + Log(musq)*(EulerGamma - Log(4*Pi)) - Inv2*
     -     Log(4*Pi)**2 + EulerGamma*Log(Abs(s)) - 
     -     EulerGamma*Log(Abs(t)) - EulerGamma*Log(Abs(u)) + (0,1)*
     -     EulerGamma*Pi*(-theta(s) + theta(t) + theta(u)) + 
     -     Log(4*Pi)*(EulerGamma - 11*Inv2*Inv3 - Log(Abs(s)) + 
     -     Log(Abs(t)) + Log(Abs(u)) - (0,1)*Pi*(-theta(s) + theta(t) + 
     -     theta(u))))

      virt=virt
     &  +factfigy*ampborn*gg(1)**2

!       print*, "proc3",virt,mur
      ampborndg(incross,glupolin,photpolin)=borndg(incross,glupolin,photpolin)
      ampborn=borndg(incross,glupolin,photpolin)

      m2s(incross+4)= m2s(incross+4)+2d0*dreal(virt*conjg(ampborn))/(4d0*pi)**2/24d0
      m2sborn(incross+4)=m2sborn(incross+4)+abs(ampborn)**2/24d0 
      m2stest(incross+4)=m2stest(incross+4)+2.*dreal(((test(1)+test(2)+test(3)+test(4)+test(5)
     &    +test(6)+test(7)+test(8)+test(9)+test(10)+test(11))*gg(1)**3*gwf(1))*conjg(ampborn))/9.

      enddo  !loop over initial state exchanged
      enddo  !loop over gluon polarization
      enddo  !loop over photon polarization

      m2s(5)=m2s(5)*(pdf(1,1)*pdf(0,2)+pdf(3,1)*pdf(0,2))
      m2s(6)=m2s(6)*(pdf(1,2)*pdf(0,1)+pdf(3,2)*pdf(0,1))
      m2stest(5)=m2stest(5)*(pdf(1,1)*pdf(0,2)+pdf(3,1)*pdf(0,2))
      m2stest(6)=m2stest(6)*(pdf(0,1)*pdf(1,2)+pdf(0,1)*pdf(3,2))
      m2sbornwopdf(5)=m2sborn(5)
      m2sbornwopdf(6)=m2sborn(6)
      m2sborn(5)=m2sborn(5)*(pdf(1,1)*pdf(0,2)+pdf(3,1)*pdf(0,2))
      m2sborn(6)=m2sborn(6)*(pdf(1,2)*pdf(0,1)+pdf(3,2)*pdf(0,1))

      if (borndebug) then

      do k=0,3
         helastestmom(k,1)=p(k,1,1)
         helastestmom(k,2)=p(k,2,1)
         helastestmom(k,3)=p(k,3,1)
         helastestmom(k,4)=v(k,1,1)
         helastestmom(k,5)=v(k,2,1)
         helastestmom(k,6)=v(k,3,1)
         helastestmom(k,7)=v(k,4,1)
      enddo

! d g -> u a e- ve~ 

      Call SDGUBEMVEBMUMMUP(helastestmom,testh)
      testh=testh*(pdf(1,1)*pdf(0,2)+pdf(3,1)*pdf(0,2))

      print*,  "born from virtuals/HELAS = ",m2sborn(5)/testh, "  =1 ?"

      do k=0,3
         helastestmom(k,1)=p(k,2,1)
         helastestmom(k,2)=p(k,1,1)
         helastestmom(k,3)=p(k,3,1)
         helastestmom(k,4)=v(k,1,1)
         helastestmom(k,5)=v(k,2,1)
         helastestmom(k,6)=v(k,3,1)
         helastestmom(k,7)=v(k,4,1)
      enddo
      Call SDGUBEMVEBMUMMUP(helastestmom,testh)
      testh=testh*(pdf(0,1)*pdf(1,2)+pdf(0,1)*pdf(3,2))

      print*,  "born from virtuals/HELAS = ",m2sborn(6)/testh, "  =1 ?"


      do mu=0,3
      boxmomtilde(mu,1)=-boxmom(mu,1) !outgoing
      boxmomtilde(mu,2)=boxmom(mu,2) !incoming
      boxmomtilde(mu,3)=boxmom(mu,3) !incoming
      boxmomtilde(mu,4)=-boxmom(mu,4) !outgoing
      enddo

      print*,  mur
      print*, "-s23",-2*dotrr(boxmomtilde(0,2),boxmomtilde(0,3))  !
      print*, "-t", 2*dotrr(boxmomtilde(0,1),boxmomtilde(0,3))   !  additional minus due to incoming conv
      print*, "-t12",2*dotrr(boxmomtilde(0,1),boxmomtilde(0,2))  !
      print*, "finished d g process"
      STOp
      endif

!   sum the virtual contributions and born

      bornelement=m2sborn(1)+m2sborn(2)+m2sborn(3)+m2sborn(4)+
     -     m2sborn(5)+m2sborn(6)
!       sumresult=m2s(1)+m2s(2)+m2s(3)+m2s(4)+m2s(5)+m2s(6)!!+bornelement
!      sumresultalt=m2stest(1)+m2stest(2)+m2stest(3)+m2stest(4)+m2stest(5)+m2stest(6)

!      do k=1,6
!         m2s(k)=m2sborn(k)
!      enddo

      sumresult=0d0
      do k=1,2
      if (ubd(k)) then
      sumresult=sumresult+m2s(k)
      endif
      if (ubg(k)) then
      sumresult=sumresult+m2s(k+2)
      endif
      if (dg(k)) then
      sumresult=sumresult+m2s(k+4)
      endif
      enddo

      badgauge=ubd(1).and.ubd(2).and.ubg(1).and.ubg(2).and.
     &   dg(1).and.dg(2)
      if (.not.badgauge) badpoints=badpoints+1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      finite after CATANI SEYMOUR I OPERATOR    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      s=2d0*dotrr(p(0,1,1),p(0,2,1))
      t=-2d0*dotrr(p(0,1,1),p(0,3,1))
      u=-2d0*dotrr(p(0,2,1),p(0,3,1))
      mursq1=4.*pi*mursq(1,1)

      finiteubd=(25*alfas*CA)/(9D0*Pi) + (5*alfas*CF)/Pi - 
     -     (11*alfas*CA*EulerGamma)/(12d0*Pi) - (3*alfas*CF*EulerGamma)/
     -     (2D0*Pi) + (alfas*CA*EulerGamma**2)/(4D0*Pi) + 
     -  (alfas*CF*EulerGamma**2)/(2D0*Pi) - (7*alfas*CA*Pi)/24d0 - 
     -     (7*alfas*CF*Pi)/12d0 - (8*alfas*Nf*TR)/(9D0*Pi) + (alfas*
     -     EulerGamma*Nf*TR)/(3D0*Pi) - 
     -  (3*alfas*CA*Log(mursq1/S))/(4D0*Pi) + (3*alfas*CF*
     -     Log(mursq1/S))/(2D0*Pi) + (alfas*CA*EulerGamma*
     -     Log(mursq1/S))/(2D0*Pi) - 
     -  (alfas*CF*EulerGamma*Log(mursq1/S))/Pi - (alfas*CA*
     -     Log(mursq1/S)**2)/(4D0*Pi) + (alfas*CF*Log(mursq1/S)**2)/(2D0*Pi) + 
     -  (5*alfas*CA*Log(-(mursq1/T)))/(6D0*Pi) - (alfas*CA*
     -     EulerGamma*Log(-(mursq1/T)))/(2D0*Pi) - (alfas*Nf*TR*
     -     Log(-(mursq1/T)))/(6D0*Pi) + 
     -  (alfas*CA*Log(-(mursq1/T))**2)/(4D0*Pi) + (5*alfas*CA*
     -     Log(-(mursq1/U)))/(6D0*Pi) - (alfas*CA*EulerGamma*
     -     Log(-(mursq1/U)))/(2D0*Pi) - 
     -  (alfas*Nf*TR*Log(-(mursq1/U)))/(6D0*Pi) + (alfas*CA*
     -     Log(-(mursq1/U))**2)/(4D0*Pi)

      finiteubg=(25*alfas*CA)/(9D0*Pi) + (5*alfas*CF)/Pi - 
     -     (11*alfas*CA*EulerGamma)/(12d0*Pi) - (3*alfas*CF*EulerGamma)/
     -     (2D0*Pi) + (alfas*CA*EulerGamma**2)/(4D0*Pi) + 
     -  (alfas*CF*EulerGamma**2)/(2D0*Pi) - (7*alfas*CA*Pi)/24d0 - 
     -     (7*alfas*CF*Pi)/12d0 - (8*alfas*Nf*TR)/(9D0*Pi) + (alfas*
     -     EulerGamma*Nf*TR)/(3D0*Pi) + 
     -  (5*alfas*CA*Log(mursq1/S))/(6D0*Pi) - (alfas*CA*EulerGamma*
     -     Log(mursq1/S))/(2D0*Pi) - (alfas*Nf*TR*Log(mursq1/S))/(6D0*Pi) + 
     -  (alfas*CA*Log(mursq1/S)**2)/(4D0*Pi) - (3*alfas*CA*
     -     Log(-(mursq1/T)))/(4D0*Pi) + (3*alfas*CF*Log(-(mursq1/T)))/
     -     (2D0*Pi) + 
     -  (alfas*CA*EulerGamma*Log(-(mursq1/T)))/(2D0*Pi) - (alfas*CF*
     -     EulerGamma*Log(-(mursq1/T)))/Pi - (alfas*CA*
     -     Log(-(mursq1/T))**2)/(4D0*Pi) + 
     -  (alfas*CF*Log(-(mursq1/T))**2)/(2D0*Pi) + (5*alfas*CA*
     -     Log(-(mursq1/U)))/(6D0*Pi) - (alfas*CA*EulerGamma*
     -     Log(-(mursq1/U)))/(2D0*Pi) - 
     -  (alfas*Nf*TR*Log(-(mursq1/U)))/(6D0*Pi) + (alfas*CA*
     -     Log(-(mursq1/U))**2)/(4D0*Pi)


      u=-2d0*dotrr(p(0,1,1),p(0,3,1))
      t=-2d0*dotrr(p(0,2,1),p(0,3,1))


      finiteubgcr=(25*alfas*CA)/(9D0*Pi) + (5*alfas*CF)/Pi - 
     -     (11*alfas*CA*EulerGamma)/(12d0*Pi) - (3*alfas*CF*EulerGamma)/
     -     (2D0*Pi) + (alfas*CA*EulerGamma**2)/(4D0*Pi) + 
     -  (alfas*CF*EulerGamma**2)/(2D0*Pi) - (7*alfas*CA*Pi)/24d0 - 
     -     (7*alfas*CF*Pi)/12d0 - (8*alfas*Nf*TR)/(9D0*Pi) + (alfas*
     -     EulerGamma*Nf*TR)/(3D0*Pi) + 
     -  (5*alfas*CA*Log(mursq1/S))/(6D0*Pi) - (alfas*CA*EulerGamma*
     -     Log(mursq1/S))/(2D0*Pi) - (alfas*Nf*TR*Log(mursq1/S))/(6D0*Pi) + 
     -  (alfas*CA*Log(mursq1/S)**2)/(4D0*Pi) - (3*alfas*CA*
     -     Log(-(mursq1/T)))/(4D0*Pi) + (3*alfas*CF*Log(-(mursq1/T)))/
     -     (2D0*Pi) + 
     -  (alfas*CA*EulerGamma*Log(-(mursq1/T)))/(2D0*Pi) - (alfas*CF*
     -     EulerGamma*Log(-(mursq1/T)))/Pi - (alfas*CA*
     -     Log(-(mursq1/T))**2)/(4D0*Pi) + 
     -  (alfas*CF*Log(-(mursq1/T))**2)/(2D0*Pi) + (5*alfas*CA*
     -     Log(-(mursq1/U)))/(6D0*Pi) - (alfas*CA*EulerGamma*
     -     Log(-(mursq1/U)))/(2D0*Pi) - 
     -  (alfas*Nf*TR*Log(-(mursq1/U)))/(6D0*Pi) + (alfas*CA*
     -     Log(-(mursq1/U))**2)/(4D0*Pi)


!       print*, m2s(1),finiteubd*(pdf(-2,1)*pdf(1,2)+pdf(-4,1)*pdf(3,2))*m2sbornwopdf(1)
!       print*, m2s(2),finiteubd*(pdf(-2,2)*pdf(1,1)+pdf(-4,2)*pdf(3,1))*m2sbornwopdf(2)
!       print*, m2s(3),finiteubg*  (pdf(-2,1)*pdf(0,2)+pdf(-4,1)*pdf(0,2))*m2sbornwopdf(3)
!       print*, m2s(4),finiteubgcr*(pdf(-2,2)*pdf(0,1)+pdf(-4,2)*pdf(0,1))*m2sbornwopdf(4)
!       print*, m2s(5),finiteubg*  (pdf( 1,1)*pdf(0,2)+pdf( 3,1)*pdf(0,2))*m2sbornwopdf(5)
!       print*, m2s(6),finiteubgcr*(pdf( 1,2)*pdf(0,1)+pdf( 3,2)*pdf(0,1))*m2sbornwopdf(6)
!       STOP

!       finiteubdf=finiteubd*(pdf(-2,1)*pdf(1,2)+pdf(-4,1)*pdf(3,2))*m2sbornwopdf(1)
!      &          +finiteubd*(pdf(-2,2)*pdf(1,1)+pdf(-4,2)*pdf(3,1))*m2sbornwopdf(2)
!       finiteqg=finiteubg*  (pdf(-2,1)*pdf(0,2)+pdf(-4,1)*pdf(0,2))*m2sbornwopdf(3)
!      &        +finiteubgcr*(pdf(-2,2)*pdf(0,1)+pdf(-4,2)*pdf(0,1))*m2sbornwopdf(4)
!      &        +finiteubg*  (pdf( 1,1)*pdf(0,2)+pdf( 3,1)*pdf(0,2))*m2sbornwopdf(5)
!      &        +finiteubgcr*(pdf( 1,2)*pdf(0,1)+pdf( 3,2)*pdf(0,1))*m2sbornwopdf(6)


      finiteubdf=0d0
      finiteqg=0d0
!       if (nlobox) then
      if (ubd(1)) then
      finiteubdf=finiteubdf+finiteubd*(pdf(-2,1)*pdf(1,2)+pdf(-4,1)*pdf(3,2))*m2sbornwopdf(1) 
      endif
      if (ubd(2)) then
      finiteubdf=finiteubdf+finiteubd*(pdf(-2,2)*pdf(1,1)+pdf(-4,2)*pdf(3,1))*m2sbornwopdf(2)     
      endif
      if (ubg(1)) then
      finiteqg=finiteqg+finiteubg*(pdf(-2,1)*pdf(0,2)+pdf(-4,1)*pdf(0,2))*m2sbornwopdf(3)  
      endif
      if (ubg(2)) then
      finiteqg=finiteqg+finiteubgcr*(pdf(-2,2)*pdf(0,1)+pdf(-4,2)*pdf(0,1))*m2sbornwopdf(4)      
      endif
      if (dg(1)) then
      finiteqg=finiteqg+finiteubg*(pdf( 1,1)*pdf(0,2)+pdf( 3,1)*pdf(0,2))*m2sbornwopdf(5)  
      endif
      if (dg(2)) then
      finiteqg=finiteqg+finiteubgcr*(pdf( 1,2)*pdf(0,1)+pdf( 3,2)*pdf(0,1))*m2sbornwopdf(6)  
      endif
!       endif


      triangleresult(1,1)=0d0
      triangleresult(1,2)=0d0
      triangleresult(1,3)=0d0
      triangleresult(2,1)=0d0
      triangleresult(2,2)=0d0
      triangleresult(2,3)=0d0
      fermres=0d0

      if (calctriangle) then

! Triangle contribution, prelimary
!
!  ub g processes
!
      mt2=tmass**2
      cptri=dcmplx(TR*gwwa/(dsqrt(sin2w*(1d0-sin2w)))*alfas/(4d0*pi),0d0)

      do incross=1,2

      do mu=0,3
      glumom(mu)=p(mu,mapspin1(incross),1)
      quarkmom(mu)=p(mu,mapspin2(incross),1)
      effglumom(mu)=quarkmom(mu)-p(mu,3,1)-wmom(mu)
      enddo
      do mu=0,3
      trianglemom(mu,1)=-1d0*glumom(mu)   ! gluon with m2=0d0 incoming
      trianglemom(mu,2)=-1d0*effglumom(mu) ! effective gluon incoming
      trianglemom(mu,3)=photmom(mu)        ! z boson  outgoing
      enddo

!      do wpolin=-1,1
!      CALL VXXXXX(wmom,wmass ,wpolin,win, wmcurr)
      do mu=0,3
         wmcurr(mu+1)=storewpol(mu)
      enddo

      CALL QBTOQBWG(quarkmom,p(0,3,1),wmcurr,effglumom,effglucurr)

      do glupolin=-1,1,2

      CALL VXXXXX(glumom,0d0 ,glupolin,-1, temp)
      do mu=0,3
      glupol(mu)=temp(mu+1)
      enddo

      do photpolin=-1,1,2
      result(1)=0d0

      CALL WZTRIANGLE(trianglemom(0,1),trianglemom(0,2),trianglemom(0,3)
     &   ,glupol,effglucurr,currlepz(0,photpolin),mt2,musq,result(1))


      triangleresult(incross,1)=triangleresult(incross,1)+2d0*dimag(result(1)
     &    *conjg(ampbornubg(mapspin2(incross),glupolin,photpolin)))*cptri/24d0

!   incross=1 interferes with gub 
!   incross=2 interferes with ubg
!   structure ampbornqg (incross,glupolin,photpolin)
!
      enddo !photpolin
      enddo !glupolin
      enddo !incross

!
!  d g processes
!
      do incross=1,2

      do mu=0,3
      glumom(mu)=p(mu,mapspin1(incross),1)
      quarkmom(mu)=p(mu,mapspin2(incross),1)
      effglumom(mu)=quarkmom(mu)-p(mu,3,1)-wmom(mu)
      enddo
      do mu=0,3
      trianglemom(mu,1)=-1d0*glumom(mu)   ! gluon with m2=0d0 incoming
      trianglemom(mu,2)=-1d0*effglumom(mu) ! effective gluon incoming
      trianglemom(mu,3)=photmom(mu)        ! z boson  outgoing
      enddo

!      do wpolin=-1,1
!      CALL VXXXXX(wmom,wmass ,wpolin,win, wmcurr)
      do mu=0,3
         wmcurr(mu+1)=storewpol(mu)
      enddo

      CALL QTOQWG(quarkmom,p(0,3,1),wmcurr,effglumom,effglucurr)

      do glupolin=-1,1,2

      CALL VXXXXX(glumom,0d0 ,glupolin,-1, temp)
      do mu=0,3
      glupol(mu)=temp(mu+1)
      enddo

      do photpolin=-1,1,2
      result(1)=0d0

      CALL WZTRIANGLE(trianglemom(0,1),trianglemom(0,2),trianglemom(0,3)
     &   ,glupol,effglucurr,currlepz(0,photpolin),mt2,musq,result(1))

      triangleresult(incross,2)=triangleresult(incross,2)-2d0*dimag(result(1)
     &    *conjg(ampborndg(mapspin2(incross),glupolin,photpolin)))*cptri/24d0

!   incross=1 interferes with gd
!   incross=2 interferes with dg
!   structure ampbornqg (incross,glupolin,photpolin)
!
      enddo !photpolin
      enddo !glupolin
      enddo !incross

!
!  ub d processes
!
      do incross=1,2

      do mu=0,3
      glumom(mu)=p(mu,3,1) 
      quarkmom(mu)=p(mu,mapspin1(incross),1)
      quarkmom2(mu)=p(mu,mapspin2(incross),1)
      effglumom(mu)=quarkmom(mu)+quarkmom2(mu)-wmom(mu)
      enddo
      do mu=0,3
      trianglemom(mu,1)=glumom(mu)  ! m2=0 outgoing
      trianglemom(mu,2)=-1d0*effglumom(mu)  !  incoming
      trianglemom(mu,3)=photmom(mu)  ! zboson outgoing
      enddo

!       do mu=0,3
!       print*, trianglemom(mu,1)+ trianglemom(mu,2)+ trianglemom(mu,3)
!       enddo
!       STOP

!      do wpolin=-1,1
      do mu=0,3
         wmcurr(mu+1)=storewpol(mu)
      enddo

      CALL QBQTOWG(quarkmom,quarkmom2,wmcurr,effglumom,effglucurr)

      do glupolin=-1,1,2

      CALL VXXXXX(glumom,0d0 ,glupolin,1, temp)
      do mu=0,3
      glupol(mu)=temp(mu+1)
      enddo

      do photpolin=-1,1,2
      result(1)=0d0

      CALL WZTRIANGLE(trianglemom(0,1),trianglemom(0,2),trianglemom(0,3)
     &   ,glupol,effglucurr,currlepz(0,photpolin),mt2,musq,result(1))

      triangleresult(incross,3)=triangleresult(incross,3)+2d0*dimag(result(1)
     &    *conjg(ampbornubd(incross,glupolin,photpolin)))*cptri/9d0

!   incross=1 interferes with ubd
!   incross=2 interferes with dub
!   structure ampbornqbq (incross,glupolin,photpolin)

      enddo !photpolin
      enddo !glupolin
      enddo !incross

      triangleresult(1,1)=triangleresult(1,1)*(pdf(0,1)*pdf(-2,2)+pdf(0,1)*pdf(-4,2)) !g ub
      triangleresult(2,1)=triangleresult(2,1)*(pdf(0,2)*pdf(-2,1)+pdf(0,2)*pdf(-4,1)) !ub g
      triangleresult(1,2)=triangleresult(1,2)*(pdf(0,1)*pdf( 1,2)+pdf(0,1)*pdf( 3,2)) !g d
      triangleresult(2,2)=triangleresult(2,2)*(pdf(0,2)*pdf( 1,1)+pdf(0,2)*pdf( 3,1)) !d g
      triangleresult(1,3)=triangleresult(1,3)*(pdf(-2,1)*pdf(1,2)+pdf(-4,1)*pdf(3,2)) !ub d
      triangleresult(2,3)=triangleresult(2,3)*(pdf(-2,2)*pdf(1,1)+pdf(-4,2)*pdf(3,1)) !d ub

      fermres=triangleresult(1,1)+triangleresult(2,1)+triangleresult(1,2)
     &  +triangleresult(2,2)+triangleresult(1,3)+triangleresult(2,3)

!      print*, triangleresult
!      print*, " "
!      print*, fermres
!      STOP
      
      sumresult=fermres+sumresult

      endif


      sumresult=sumresult+finiteqg+finiteubdf

      endif


      end


      subroutine WZTRIANGLE(k1,k2,k3,e1,e2,e3,mt2,musq,result)

      implicit none

      integer mu
      double precision k1(0:3),k2(0:3),k3(0:3),mt2,musq,mt
      double complex e1(0:3),e2(0:3),e3(0:3),result

      logical ldebug
      parameter (ldebug=.false.)

!
!  local variables
!
      double complex pair1, pair2, pair3, pair4, pair5,
     &  eps1,eps2,eps3,eps4,eps5,abb1,abb2,abb3,abb4,abb5,
     &  abbsum1,abbsum2,abbsum3,abbsum4,mg2,mz2

      double complex dotcr,epscccr,epsccrr,b0ce
      double precision dotrr,psqt

      external dotcr,dotrr,epscccr,epsccrr,b0ce

      complex*16 i3point
      external i3point

      mt=dsqrt(mt2)
      mg2=dcmplx(dotrr(k2,k2))
      mz2=dcmplx(dotrr(k3,k3))
      psqt=dotrr(k2,k2)+2d0*dotrr(k2,k3)+dotrr(k3,k3)

      Pair1 = 0d0!dotcr(e1, k1)
      Pair3 = dotcr(e1, k2)
      Pair2 = dotcr(e2, k1)
      Pair4 = dotcr(e2, k2)
      Pair5 = dotcr(e3, k2)
      Eps1 = Epscccr(e1, e2, e3, k1)
      Eps2 = Epscccr(e1, e2, e3, k2)
      Eps5 = Epsccrr(e1, e2, k1, k2)
      Eps4 = Epsccrr(e1, e3, k1, k2)
      Eps3 = Epsccrr(e2, e3, k1, k2)
      Abb1 = 0d0!Eps3*Pair1
      Abb2 = Eps4*Pair2
      Abb3 = Eps3*Pair3
      Abb4 = Eps4*Pair4
      Abb5 = Eps5*Pair5 
      AbbSum2 = Abb2 - Abb3
      AbbSum3 = Abb3 + Abb4 + Abb5
      AbbSum4 = Eps1 - Eps2
      AbbSum1 = Eps1 + Eps2

      result=(-Eps2 - (2*AbbSum2*mg2)/(mg2 - mz2)**2 + (-2*Abb4 + 
     -     2*AbbSum3 + AbbSum4*mg2)/(-mg2 + mz2))*
     -   B0ce(mg2,0d0,0d0,musq) + (Eps2 + (2*AbbSum2*mg2)/(mg2 - mz2)**2 + 
     -     (-2*Abb4 + 2*AbbSum3 + AbbSum4*mg2)/(mg2 - mz2))*
     -     B0ce(mg2,mt2,mt2,musq) + 
     -  (-AbbSum1 - AbbSum4 + 2*Eps1 + Eps2 + (2*AbbSum2*mg2)/(mg2 - mz2)**2 + 
     -     (-2*(Abb4 - AbbSum3) + AbbSum4*mg2)/(mg2 - mz2))*
     -     B0ce(mz2,0d0,0d0,musq) + 
     -  (AbbSum1 + AbbSum4 - 2*Eps1 - Eps2 - (2*AbbSum2*mg2)/(mg2 - mz2)**2 + 
     -     (-2*(Abb4 - AbbSum3) + AbbSum4*mg2)/(-mg2 + mz2))*
     -     B0ce(mz2,mt2,mt2,musq) + 
     -  (2*mt2*(2*AbbSum2 + AbbSum4*(mg2 - mz2))*
     -     i3point(mt,mg2,mz2,psqt))/(mg2 - mz2)


      if (ldebug) then
      print*, "input"
      print*," "
      print*,"k1",k1
      print*,"k2",k2
      print*,"k3",k3
      print*,"e1",e1
      print*,"e2",e2
      print*,"e3",e3
      print*,"mt2",mt2
      print*,"musq", musq
      print*," "
      print*,"momentum conservation"
      do mu=0,3
      print*,k1(mu)+k2(mu)+k3(mu)
      enddo
      print*, " "
      print*,"output "
      print*," "
      print*,Pair1, "= dotcr(e1, k1)"
      print*,Pair3, " = dotcr(e1, k2)"
      print*,Pair2, " = dotcr(e2, k1) "
      print*,Pair4, " = dotcr(e2, k2)"
      print*,Pair5, " = dotcr(e3, k2)"
      print*,Eps1, " = Epscccr(e1, e2, e3, k1)"
      print*,Eps2, " = Epscccr(e1, e2, e3, k2)"
      print*,Eps5, " = Epsccrr(e1, e2, k1, k2)"
      print*,Eps4, " = Epsccrr(e1, e3, k1, k2)"
      print*,Eps3, " = Epsccrr(e2, e3, k1, k2)"
      print*,Abb1, " = Eps3*Pair1 "
      print*,Abb2, " = Eps4*Pair2"
      print*,Abb3, " = Eps3*Pair3"
      print*,Abb4, " = Eps4*Pair4"
      print*,Abb5, " = Eps5*Pair5"
      print*,AbbSum2, " = Abb2 - Abb3"
      print*,AbbSum3, " = Abb3 + Abb4 + Abb5"
      print*,AbbSum4, " = Eps1 - Eps2"
      print*,AbbSum1, " = Eps1 + Eps2"
      print*," "
      print*, "scalar function B0ce(mg2,0d0,0d0,musq) =",B0ce(mg2,0d0,0d0,musq)
      print*, "scalar function B0ce(mg2,mt2,mt2,musq) =",B0ce(mg2,mt2,mt2,musq)
      print*, "scalar function B0ce(mz2,mt2,mt2,musq) =",B0ce(mz2,mt2,mt2,musq)
      print*, "scalar function C0ce(mg2,mz2,0d0,mt2,mt2,mt2,musq) =", i3point(dsqrt(mt2),mg2,mz2,psqt)
      print*," "
      print*,"result", result
!       STOP
      endif


      end


      double complex function b0ce(p1sq,m1s,m2s,musq)

      implicit none

#include "VBFNLO/utilities/global.inc"

      double precision m1s,m2s,musq
      double complex p1sq
      double precision eulergamma
      parameter (eulergamma=0.5772156649015329d0)
      double complex qlI2,b0eval(0:2)
      external qlI2

! this is the Ellis-Zanderighi B Function 
      b0eval(0) = qlI2(p1sq,m1s,m2s,musq,0)
! this is the conversion factor to LoopTools
      b0eval(1) =(-eulergamma+log(4d0*pi))*qlI2(p1sq,m1s,m2s,musq,-1)
!      b0eval(2) =(-eulergamma+log(4d0*pi))*qlI2(p1sq,m1s,m2s,musq,-2) !=0

      b0ce=b0eval(0)+b0eval(1)!+b0eval(2)
      return
      end


      subroutine finitecollubd(z,x,pa,pb,pin,mufac,quark,result)

      implicit none

#include "VBFNLO/utilities/global.inc"

      integer quark

      double precision z,x,pa(0:3),pb(0:3),pin(0:3),result,dotrr,zox,pdfgzox,
     &   pdfqzox, pdfqz,pdftmp(-6:6),mufac,mufacsq,NF,
     &   paprimepb2,paprimepi2

      logical ifxmzge0

      external dotrr

      NF=5d0
      mufacsq=mufac**2
      zox=z/x

      if ((x-z).ge.0) then
      ifxmzge0=.true.
      call pdfproton( zox, mufac, pdftmp(-6) )
      pdfgzox=pdftmp(0)/zox
      pdfqzox=pdftmp(quark)/zox
      else
      ifxmzge0=.false.
      pdfgzox=0d0
      pdfqzox=0d0
      endif
      paprimepb2=2d0*dotrr(pa,pb)
      paprimepi2=2d0*dotrr(pa,pin)

      call pdfproton( z, mufac, pdftmp(-6) )
      pdfqz=pdftmp(quark)/z


      result=(pdfqz*(-113 + 2*Nf + 6*Pi**2 - 18*Log(mufacsq/paprimepi2)*(3 + 4*Log(1 - z)) + 
     -      Log(mufacsq/paprimepb2)*(6 + 8*Log(1 - z)) + 2*Log(1 -z)*(-33 + 2*Nf + 14*Log(1 - z))))/24.

      if (ifxmzge0) then
      result=result+
     - ((3*pdfgzox*(-1 + x)*(1 + 2*(-1 + x)*x) - 8*(pdfqzox - 2*pdfqz*x + pdfqzox*x**2))*Log(mufacsq/paprimepb2) - 
     -    9*(3*pdfgzox*(-1 + x)*(1 + 2*(-1 + x)*x) - 8*(pdfqzox - 2*pdfqz*x + pdfqzox*x**2))*Log(mufacsq/paprimepi2) + 
     -    7*(3*pdfgzox*(-1 + x)*(1 + 2*(-1 + x)*x) - 8*(pdfqzox - 2*pdfqz*x + pdfqzox*x**2))*Log(1 - x) + 
     -    4*(((-33 + 2*Nf)*pdfqz - 12*pdfgzox*(-1 + x)**2)*x + pdfqzox*(17 - 2*Nf - 16*(-2 + x)*x) + 
     -       2*(-1 + x)*(-8*pdfqzox*(1 + x) + pdfgzox*(3 + 6*(-1 + x)*x))*Log(1./x) + 32*pdfqzox*Log(x)))/
     -  (48.*(-1 + x)*x)
      endif

!      print*, result, pdfqzox, pdfgzox
!      STOP

      end


      subroutine finitecollgq(z,x,pa,pb,pin,mufac,result)

      implicit none

#include "VBFNLO/utilities/global.inc"

      integer i

      double precision z,x,pa(0:3),pb(0:3),pin(0:3),result,dotrr,zox,pdfgzox,
     &   pdfqzox, papi,papb,pdftmp(-6:6),mufacsq,pdfgz,NF,mufac,paprimepi2,
     &   paprimepb2,pdfqbzox

      logical ifxmzge0

      external dotrr

      pdfqzox=0d0
      pdfqbzox=0d0
      result=0d0

      NF=5d0

      mufacsq=mufac**2  ! q_sf is SQRT[mufsq]

      zox=z/x

      if ((x-z).ge.0) then

      ifxmzge0=.true.
      call pdfproton( zox, mufac, pdftmp(-6) )
      pdfgzox=pdftmp(0)/zox
!
!     Sum over the Quark->Gluon-Splittings
! 
        do i=1,4
         pdfqzox=pdfqzox+pdftmp(i)/zox
        enddo
!
!     Sum over the Antiquark->Gluon-Splittings
! 
        do i=-4,-1
         pdfqzox=pdfqzox+pdftmp(i)/zox
        enddo
!
      else

      ifxmzge0=.false.
      pdfgzox=0d0
      pdfqzox=0d0

      endif

      papb=dotrr(pa,pb)
      papi=dotrr(pa,pin)
      paprimepi2=2d0*papi
      paprimepb2=2d0*papb
!
!    Get gluon-pdf
!
      call pdfproton( z, mufac, pdftmp(-6) )
      pdfgz=pdftmp(0)/z

      result=(pdfgz*(-681 + 32*NF + 18*Pi**2 + 3*Log(mufacsq/paprimepb2)*(-33 + 2*NF - 72*Log(1 - z)) + 
     -      3*Log(mufacsq/paprimepi2)*(-33 + 2*NF - 72*Log(1 - z)) + 162*Log(1 - z)*(-1 + 2*Log(1 - z))))/72.


      if (ifxmzge0) then
      result=result+
     -  (-4*(9*pdfgz*x**2 + 2*pdfqzox*(-1 + x)*(2 + (-2 + x)*x) - 9*pdfgzox*(1 + (-1 + x)*x)**2)*
     -     Log(mufacsq/paprimepb2) - 4*(9*pdfgz*x**2 + 2*pdfqzox*(-1 + x)*(2 + (-2 + x)*x) - 
     -       9*pdfgzox*(1 + (-1 + x)*x)**2)*Log(mufacsq/paprimepi2) + 
     -    12*(9*pdfgz*x**2 + 2*pdfqzox*(-1 + x)*(2 + (-2 + x)*x) - 9*pdfgzox*(1 + (-1 + x)*x)**2)*Log(1 - x) - 
     -    8*(-1 + x)*(-2*pdfqzox*(2 + (-2 + x)*x) + 9*pdfgzox*(-1 + x*(2 + (-1 + x)*x)))*Log(1/x) + 
     -    x*(27*pdfgzox + (-27*pdfgz + 16*pdfqzox*(-1 + x))*x + 72*pdfgzox*Log(x)))/(12.*(-1 + x)*x**2)
      endif

      end


      subroutine finitecollqqg(z,x,pa,pb,pin,mufac,quark,result)

      implicit none

#include "VBFNLO/utilities/global.inc"

      integer quark

      double precision z,x,pa(0:3),pb(0:3),pin(0:3),result,dotrr,zox,pdfgzox,
     &   pdfqzox, pdfqz,paprimepi2,paprimepb2,pdftmp(-6:6),mufacsq,NF,pdfgz,mufac

      logical ifxmzge0

      external dotrr

      result=0d0
      NF=5d0
      mufacsq=mufac**2

      zox=z/x

      if ((x-z).ge.0) then
      ifxmzge0=.true.
      call pdfproton( zox, mufac, pdftmp(-6) )
      pdfgzox=pdftmp(0)/zox
      pdfqzox=pdftmp(quark)/zox
      else
      ifxmzge0=.false.
      pdfgzox=0d0
      pdfqzox=0d0
      endif
      paprimepb2=2d0*dotrr(pa,pb)
      paprimepi2=2d0*dotrr(pa,pin)

      call pdfproton( z, mufac, pdftmp(-6) )
      pdfqz=pdftmp(quark)/z
      pdfgz=pdftmp(0)/z


      result=(pdfqz*(-231 - 2*Pi**2 - 54*Log(mufacsq/paprimepb2)*(3 + 
     -     4*Log(1 - z)) + 
     -      6*Log(mufacsq/paprimepi2)*(3 + 4*Log(1 - z)) + 6*Log(1 - z)*
     -     (3 + 34*Log(1 - z))))/72.

      if(ifxmzge0) then
      result=result+
     -   (-9*(3*pdfgzox*(-1 + x)*(1 + 2*(-1 + x)*x) - 8*(pdfqzox - 
     -     2*pdfqz*x + pdfqzox*x**2))*Log(mufacsq/paprimepb2) + 
     -    (3*pdfgzox*(-1 + x)*(1 + 2*(-1 + x)*x) - 8*(pdfqzox - 
     -     2*pdfqz*x + pdfqzox*x**2))*Log(mufacsq/paprimepi2) + 
     -    17*(3*pdfgzox*(-1 + x)*(1 + 2*(-1 + x)*x) - 8*(pdfqzox - 
     -     2*pdfqz*x + pdfqzox*x**2))*Log(1 - x) + 
     -    4*(3*(pdfqz - 4*pdfgzox*(-1 + x)**2)*x + pdfqzox*(-19 - 
     -     16*(-2 + x)*x) + 
     -       2*(-1 + x)*(-8*pdfqzox*(1 + x) + 
     -     pdfgzox*(3 + 6*(-1 + x)*x))*Log(1./x) + 32*pdfqzox*Log(x)))/
     -  (48.*(-1 + x)*x)

      endif


      end


      SUBROUTINE QBTOQBWG(pup,pdown,wmwave,gmom,currg)
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : u~ -> d~ w- g  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   2,NEIGEN=  1) 
#include "genps.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   6, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL),pup(0:3),pdown(0:3),gmom(0:3)
C  
C LOCAL VARIABLES 
C  
      INTEGER mu
      COMPLEX*16 AMP(NGRAPHS),currg(0:3),wmwave(6)
      COMPLEX*16 W(6,NWAVEFUNCS)

#include "VBFNLO/utilities/coupl.inc"
C ----------
C BEGIN CODE
C ----------
      do mu=0,3
      p(mu,1)=pup(mu)
      p(mu,2)=pdown(mu)
      enddo

      do mu=1,6
      w(mu,3)=wmwave(mu)
      enddo

      CALL OXXXXX(P(0,1   ),ZERO ,1,-1,W(1,1   ))        
      CALL IXXXXX(P(0,2   ),ZERO ,1,-1,W(1,2   ))        

      do mu=0,3

      CALL VCARTX(gmom,zero ,zero,mu,1,W(1,4)) 
!       CALL VXXXXX(P(0,3   ),WMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
!       CALL VXXXXX(P(0,4   ),ZERO ,NHEL(4   ),+1*IC(4   ),W(1,4   ))        
      CALL FVIXXX(W(1,2   ),W(1,3   ),GWF ,ZERO    ,ZERO    ,W(1,5   ))    
      CALL IOVXXX(W(1,5   ),W(1,1   ),W(1,4   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,1   ),W(1,3   ),GWF ,ZERO    ,ZERO    ,W(1,6   ))    
      CALL IOVXXX(W(1,2   ),W(1,6   ),W(1,4   ),GG ,AMP(2   ))             
      currg(mu) = -AMP(   1)-AMP(   2)

      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END


      SUBROUTINE QTOQWG(pdown,pup,wmwave,gmom,currg)
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : d -> u w- g  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   2,NEIGEN=  1) 
#include "genps.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   6, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL),pup(0:3),pdown(0:3),gmom(0:3)
C  
C LOCAL VARIABLES 
C  
      INTEGER mu
      COMPLEX*16 AMP(NGRAPHS),currg(0:3),wmwave(6)
      COMPLEX*16 W(6,NWAVEFUNCS)

#include "VBFNLO/utilities/coupl.inc"
C ----------
C BEGIN CODE
C ----------

      do mu=0,3
      p(mu,1)=pdown(mu)
      p(mu,2)=pup(mu)
      enddo

      do mu=1,6
      w(mu,3)=wmwave(mu)
      enddo

      CALL IXXXXX(P(0,1   ),ZERO ,-1 ,+1,W(1,1   ))        
      CALL OXXXXX(P(0,2   ),ZERO ,-1 ,+1,W(1,2   ))        


      do mu=0,3

      CALL VCARTX(gmom,zero ,zero,mu,1,W(1,4)) 

!       CALL VXXXXX(P(0,3   ),WMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
!       CALL VXXXXX(P(0,4   ),ZERO ,NHEL(4   ),+1*IC(4   ),W(1,4   ))   
     
      CALL FVOXXX(W(1,2   ),W(1,3   ),GWF ,ZERO    ,ZERO    ,W(1,5   ))    
      CALL IOVXXX(W(1,1   ),W(1,5   ),W(1,4   ),GG ,AMP(1   ))             
      CALL FVIXXX(W(1,1   ),W(1,3   ),GWF ,ZERO    ,ZERO    ,W(1,6   ))    
      CALL IOVXXX(W(1,6   ),W(1,2   ),W(1,4   ),GG ,AMP(2   ))             
      currg(mu) = -AMP(   1)-AMP(   2)
      Enddo
      END


      SUBROUTINE QBQTOWG(pup,pdown,wmwave,gmom,currg)
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : u~ d -> w- g  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   2,NEIGEN=  1) 
#include "genps.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   6, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL),pup(0:3),pdown(0:3),gmom(0:3)
C  
C LOCAL VARIABLES 
C  
      INTEGER mu
      COMPLEX*16 AMP(NGRAPHS),currg(0:3),wmwave(6)
      COMPLEX*16 W(6,NWAVEFUNCS)

#include "VBFNLO/utilities/coupl.inc"
C ----------
C BEGIN CODE
C ----------
      do mu=0,3
      p(mu,1)=pup(mu)
      p(mu,2)=pdown(mu)
      enddo

      do mu=1,6
      w(mu,3)=wmwave(mu)
      enddo


      CALL OXXXXX(P(0,1   ),ZERO ,1,-1,W(1,1   ))        
      CALL IXXXXX(P(0,2   ),ZERO ,-1,+1,W(1,2   ))

      do mu=0,3

      CALL VCARTX(gmom,zero ,zero,mu,1,W(1,4)) 

!      CALL VXXXXX(P(0,3   ),WMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
!      CALL VXXXXX(P(0,4   ),ZERO ,NHEL(4   ),+1*IC(4   ),W(1,4   ))        
      CALL FVIXXX(W(1,2   ),W(1,3   ),GWF ,ZERO    ,ZERO    ,W(1,5   ))    
      CALL IOVXXX(W(1,5   ),W(1,1   ),W(1,4   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,1   ),W(1,3   ),GWF ,ZERO    ,ZERO    ,W(1,6   ))    
      CALL IOVXXX(W(1,2   ),W(1,6   ),W(1,4   ),GG ,AMP(2   ))             
      currg(mu) = -AMP(   1)-AMP(   2)
      enddo
      END


      double complex function qlI2(p1,m1,m2,mu2,ep) 
      implicit none

      double precision pi,pisq,pisqo6
      parameter(pi=3.14159265358979d0,pisq=pi*pi,pisqo6=pisq/6d0)

      double precision zip,half,one,two,three,four,eight
      parameter(zip=0d0,half=0.5d0,one=1d0,two=2d0)
      parameter(three=3d0,four=4d0,eight=8d0)

      double complex im,impi,czip,chalf,cone,ctwo,c2ipi
      parameter(im=(0d0,1d0),impi=(0d0,3.14159265358979d0),
     . czip=(0d0,0d0),chalf=(0.5d0,0d0),cone=(1d0,0d0),ctwo=(2d0,0d0),
     . c2ipi=2d0*impi)


      double precision m1,m2,mu2
      double complex p1
      integer ep
C     p1=p1(1) is the squared four-momentum of the external particle i 
C     mi=m(i)^2, i=1,2 are the squares of the mass of the propagator i 
C     mu2 is the square of the scale mu
C     ep=-2,-1,0 chooses the coefficient in the Laurent series.
      double precision p1o,m1o,m2o,mu2o,pp1,mm1,mm2,newmu2,scalefac
      logical qlzero 
      double complex Ival(-2:0),qlI2fin
      data p1o,m1o,m2o,mu2o/3*0d0,-1d0/
      save Ival,p1o,m1o,m2o,mu2o

C--If we have already calculated, use the saved value
C--else setup the arrays
      if    ((dble(p1) .eq. p1o)
     . .and. (m1 .eq. m1o)
     . .and. (m2 .eq. m2o)
     . .and. (mu2 .eq. mu2o)) then
      qlI2=Ival(ep)
      return
      else
      p1o=p1
      m1o=m1
      m2o=m2
      mu2o=mu2
      endif
      if (mu2 .le. 0d0) then
      write(6,*) 'stopping because mu2 .le. 0d0 in qlI2, mu2=',mu2
      write(6,*) 'Rerun with positive mu2'
      stop
      endif

      scalefac=max(abs(p1),abs(m1),abs(m2),abs(mu2))
      pp1=p1/scalefac
      mm1=m1/scalefac
      mm2=m2/scalefac
      newmu2=mu2/scalefac

      Ival(-2)=czip
      Ival(-1)=czip
      Ival(0)=czip
      if ((qlzero(pp1)).and.(qlzero(mm1)).and.(qlzero(mm2))) then
      qlI2=Ival(ep)
      else
      Ival(-1)=cone
      Ival(0)=qlI2fin(pp1,mm1,mm2,newmu2)
      qlI2=Ival(ep)
      endif

      return
      end


      double complex function qlI2fin(p1sq,m0s,m1s,musq) 
C---- Implementation of the formulae of Denner and Dittmaier
C----%\cite{Denner:2005nn}
C----\bibitem{Denner:2005nn}
C----  A.~Denner and S.~Dittmaier,
C----  %``Reduction schemes for one-loop tensor integrals,''
C----  Nucl.\ Phys.\  B {\bf 734}, 62 (2006)
C----  [arXiv:hep-ph/0509141].
C----  %%CITATION = NUPHA,B734,62;%%
      implicit none

      double precision pi,pisq,pisqo6
      parameter(pi=3.14159265358979d0,pisq=pi*pi,pisqo6=pisq/6d0)

      double precision zip,half,one,two,three,four,eight
      parameter(zip=0d0,half=0.5d0,one=1d0,two=2d0)
      parameter(three=3d0,four=4d0,eight=8d0)

      double complex im,impi,czip,chalf,cone,ctwo,c2ipi
      parameter(im=(0d0,1d0),impi=(0d0,3.14159265358979d0),
     . czip=(0d0,0d0),chalf=(0.5d0,0d0),cone=(1d0,0d0),ctwo=(2d0,0d0),
     . c2ipi=2d0*impi)

      double precision p1sq,m0s,m1s,m0sq,m1sq,musq
      double complex xp,xm,b,rt,arg,arg1,qlfndd,Vcln
      logical qlzero

      m0sq=min(m0s,m1s)
      m1sq=max(m0s,m1s)

      if ((qlzero(abs(p1sq/musq)))
     . .and. (qlzero(abs(m0sq/musq)))
     . .and. (qlzero(abs(m1sq/musq)))) then
      write(6,*) 'setting psq=m0sq=m1sq=0  self-energy to zero'
      write(6,*) 'p1sq,m0sq,m1sq=',p1sq,m0sq,m1sq
      qlI2fin=czip
      return

      elseif (qlzero(m0sq/musq)) then
      arg=dcmplx(1d0-m1sq/p1sq)
      arg1=dcmplx((m1sq-p1sq)/musq)

C---deal with special cases for m0sq=0
C----- (a,0,a)  p1sq=m1sq, DD(4.13)
      if (qlzero(abs(arg1))) then
      qlI2fin=dcmplx(log(musq/m1sq))+ctwo
C----- (0,0,a) 
      elseif (qlzero(abs(p1sq/musq))) then
      qlI2fin=dcmplx(log(musq/m1sq))+cone
C----- (a,0,0)
      elseif (qlzero(abs(m1sq/musq))) then
      qlI2fin=-Vcln(arg1,-1d0)+ctwo
      else
C----- (a,0,c)
      qlI2fin=-Vcln(arg1,-1d0)+cone-qlfndd(0,arg,1d0)
      endif
      return
      elseif (qlzero(abs(p1sq/musq))) then
C---deal with special case, p1sq=0
      if (qlzero(abs((m1sq-m0sq)/musq))) then    ! (m1sq = m0sq)
      qlI2fin=dcmplx(log(musq/m0sq))
      else
      xp=dcmplx(m0sq/(m0sq-m1sq))  ! other root is formally infinite
      qlI2fin=dcmplx(log(musq/m0sq))-qlfndd(0,xp,1d0)
      endif     
      else
C----general case, DD (4.8)
      b=dcmplx(m1sq-m0sq-p1sq)
      rt=sqrt(dcmplx((m1sq-m0sq-p1sq)**2-4d0*p1sq*m0sq))
      xp=0.5d0*(-b+rt)/p1sq
      xm=0.5d0*(-b-rt)/p1sq

      qlI2fin=dcmplx(log(musq/m0sq))-qlfndd(0,xp,1d0)-qlfndd(0,xm,-1d0)

      endif
      return 

      end


      logical function qlzero(psq)
      implicit none
      
      double precision qlonshellcutoff
      parameter(qlonshellcutoff=1d-10)

      double precision psq
      if (abs(psq) .lt. qlonshellcutoff) then
      qlzero=.true.
      else
      qlzero=.false.
      endif
      return
      end


      double complex function Vcln(z, isig)
      implicit none
      double precision pi
      parameter(pi=3.14159265358979d0)

      double complex z
      double precision isig

      if( DIMAG(z) .eq. 0d0 .and. DBLE(z) .le. 0d0 ) then
         Vcln = log(-z) + DCMPLX(0D0, sign(pi, isig))
      else
         Vcln = log(z)
      endif
      end

      double complex function qlfndd(n,x,iep)
C----Implementation of DD Eq. 4.11
C----%\cite{Denner:2005nn}
C----\bibitem{Denner:2005nn}
C----  A.~Denner and S.~Dittmaier,
C----  %``Reduction schemes for one-loop tensor integrals,''
C----  Nucl.\ Phys.\  B {\bf 734}, 62 (2006)
C----  [arXiv:hep-ph/0509141].
C----  %%CITATION = NUPHA,B734,62;%%

      implicit none

      double precision pi,pisq,pisqo6
      parameter(pi=3.14159265358979d0,pisq=pi*pi,pisqo6=pisq/6d0)

      double precision zip,half,one,two,three,four,eight
      parameter(zip=0d0,half=0.5d0,one=1d0,two=2d0)
      parameter(three=3d0,four=4d0,eight=8d0)

      double complex im,impi,czip,chalf,cone,ctwo,c2ipi
      parameter(im=(0d0,1d0),impi=(0d0,3.14159265358979d0),
     . czip=(0d0,0d0),chalf=(0.5d0,0d0),cone=(1d0,0d0),ctwo=(2d0,0d0),
     . c2ipi=2d0*impi)

      integer j,n,infty
      double complex xm1,x,Vcln
      double precision iep
      logical qlzero
      parameter(infty=16) ! number of terms in sum
      
      xm1=x-cone
      if (abs(x) .lt. 10d0) then
        if (qlzero(abs(x-cone))) then
             qlfndd=czip
             else
          qlfndd=(cone-dcmplx(x**(n+1)))*(Vcln(xm1,iep)-Vcln(x,iep))
          endif
        do j=0,n
          qlfndd=qlfndd-dcmplx(x**(n-j))/dfloat(j+1)
        enddo
      elseif (abs(x) .ge. 10d0) then
        qlfndd=Vcln(cone-cone/x,iep)
        do j=n+1,n+infty
          qlfndd=qlfndd+dcmplx(x**(n-j))/dfloat(j+1)
        enddo
      endif
           
      return
      end
