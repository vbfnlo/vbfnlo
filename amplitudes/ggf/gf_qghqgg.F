***********************************************************************
      SUBROUTINE Higgs3jQG(PBAR,FSIGN,FPOL,KBAR,GSIGN,PHBAR,HSIGN,RHsum)
*----------------------------------------------------------------------
c... author: Michael Kubocz
c... date  : 15.02.2009
c
c  modified: Francisco Campanario 
c  add:  Ward identities for pentagons and hexagons
c       + control of instabilities in Dble and QUAD
***********************************************************************
      IMPLICIT NONE

#include "ggf_process.inc"
#include "ggf_dotproducts.inc"
#include "ggfloMG.inc" 
      
      logical eva
      real*8 phc(0:3)
      common/phc1/phc,eva
      
      logical SingularC(2)
      common/singluonver/SingularC
      
      integer InterpolCall(4),ICallF(4),SingularF(4)
      COMMON/ICall/InterpolCall,ICallF,SingularF
      integer GCcounter(4)
      logical GaugeFail
      COMMON/GaugeCHeck/GCcounter,GaugeFail

c... Information on coupling constants from KOPPELN
      real*8 ALFAS,XMT,ALFA,XMZ,XMW,SIN2W,XMH
      real*8 clr, xm2, xmg, b, v, a
      COMMON /BKOPIN/ ALFAS,XMT,ALFA,XMZ,XMW,SIN2W,XMH
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)

      character*17 wr2
      CHARACTER*50 wr3,wr4,wr1
      Parameter(wr4='(A10,I1,A1,I2,A2,E28.20,E28.20)',
     &          wr1='(A14,E28.20,E28.20)',
     &          wr2='(A10,E28.20)',wr3='(A5,I1,A2,E28.20)')

c... Momenta of all particles (2 quarks, 3 gluons, Higgs)
      real*8 PBAR(0:3,2),KBAR(0:4,3),PHBAR(0:4)
      real*8 K(0:4,3),P(0:4,2),PH(0:4)
c... 4-dim momenta variables
      real*8 Kp(0:3,3),Pp(0:3,2),PHp(0:3)
      real*8 DUMMYp(0:3),Kpenp(0:3,4)
c...
      integer FSIGN(2),FPOL(2),SIGMAF(2),GSIGN(3),HSIGN

c... momenta combinations
      real*8 p21(0:4),p21K1(0:4),p21K2(0:4),p21K3(0:4),
     &       p21K12(0:4),p21K13(0:4),p21K23(0:4),
     &       p21K1H(0:4),p21K2H(0:4),p21K3H(0:4),
     &       K12(0:4),K13(0:4),K23(0:4),K123(0:4),
     &       PHK1(0:4),PHK2(0:4),PHK3(0:4),PHp21(0:4),
     &       PHK13(0:4),PHK23(0:4),PHK12(0:4)
c... 4-dim momenta variables
      real*8 p21p(0:3),p21K1p(0:3),p21K2p(0:3),p21K3p(0:3),
     &       p21K12p(0:3),p21K13p(0:3),p21K23p(0:3),
     &       p21K1Hp(0:3),p21K2Hp(0:3),p21K3Hp(0:3),
     &       K12p(0:3),K13p(0:3),K23p(0:3),K123p(0:3),
     &       PHK1p(0:3),PHK2p(0:3),PHK3p(0:3),PHp21p(0:3),
     &       PHK13p(0:3),PHK23p(0:3),PHK12p(0:3)

c... Polarization vectors
      real*8 epsg(0:3,2,3)
      complex*16 epsgC(0:3,2,3)
      integer i1pol,i2pol,i3pol,i1,i2,i3,i4,
     &     i1polmax,i2polmax,i3polmax

c... Effective polarization vectors (with Higgs-vertex)
      complex*16 epsHeff(0:3,2,3),j21Heff(0:3),
     &           epsHtop(0:3,2,3),j21Htop(0:3)

c... Spinors and quark currents with and without gluon emision
      complex*16 psi(2,-1:1,2),braketg(2,-1:1,12,2),
     &         j21aux(0:3,-1:1),j21auxi(0:3,-1:1,2),j21auxf(0:3,-1:1,2),
     &         j21(0:3),j21i(0:3,2,3),j21f(0:3,2,3),j21if(0:3,-1:1,2,2,6)
      integer isig,sigma

c... Quark currents with gluon and Higgs emission
      complex*16 j21Hefff(0:3,2,3),j21Heffi(0:3,2,3),
     &           j21Heffif(0:3,-1:1,2,2,6),j21iHefff(0:3,-1:1,2,2,6),
     &           j21Htopf(0:3,2,3),j21Htopi(0:3,2,3),
     &           j21Htopif(0:3,-1:1,2,2,6),j21iHtopf(0:3,-1:1,2,2,6)


c... Auxiliary "current" variables
      complex*16 GGG12(0:4,2,2),GGG13(0:4,2,2),GGG23(0:4,2,2),
     &           GGG21(0:4,2,2),GGG31(0:4,2,2),GGG32(0:4,2,2),
     &           GGGj211(0:4,2),GGG2j21(0:4,2),GGG3j21(0:4,2),
     &           GGG123(0:4,2,2,2),
     &           GGG132a(0:4,2,2,2),GGG132b(0:4,2,2,2)

      complex*16 GGGPH12(0:4,2,2),GGGPH13(0:4,2,2),GGGPH23(0:4,2,2),
     &           GGGPH32(0:4,2,2),
     &           GGG1PH2(0:4,2,2),GGG1PH3(0:4,2,2),
     &           GGGPA12(0:4,2,2),GGGPA13(0:4,2,2),GGGPA23(0:4,2,2),
     &           GGGPA32(0:4,2,2),
     &           GGG1PA2(0:4,2,2),GGG1PA3(0:4,2,2)

      complex*16 GGGPHtop12(0:4,2,2),GGGPHtop13(0:4,2,2),
     &           GGGPHtop23(0:4,2,2),GGGPHtop32(0:4,2,2),
     &           GGG1PHtop2(0:4,2,2),GGG1PHtop3(0:4,2,2)

c... Effective vertices and loop functions     
      complex*16 H2G_Loop,H3G_Loop
      complex*16 H2G_effective,H3G_effective,GGG_vertex,H4G_effective2
      complex*16 qpentH_top1234(4),qpentH_top1342(4),qpentH_top1423(4)

c... Color matrix
      real*8 RGB(6,6)

c... Further local variables
      complex*16 color(6),colorT,colorPen(6),ampQG(87)
      complex*16 Heff3j(87),Htop3j(87),Hbottom3j(87),
     &           Aeff3j(87),Atop3j(87),Abottom3j(87)
      integer ienter,mu,i,j,l
      complex*16 auxPen,auxRest
      real*8 pi,fac,vev,DUMMY(0:4),Kpen(0:4,4),RH,RHpenSq,RHpenI,RHpenI2,
     &       RHsum,RHsumFi
      complex*16 EPSpen(0:3,15)
      parameter (pi=3.14159 26535 89793d0)
      save fac
      data ienter /0/
      save ienter
c... Debbuging variables
      real*8 polar(0:3,4), phiggs(0:3),ratioVal
      complex*16 polarC(0:3,4)
      complex*16 H3G_LoopGC,H4G_Loop,pentP,epsChelas(1:6,-1:1,3)
      real*8 aux1(0:3,4),sum(0:3),auxGauge,auxGauge2,ratioPB
      complex*16 GCepsg(0:3,4),Pen2134,Box1,Box2,GCres,aux1C(0:3,4)
      complex*16 pent1,pent2,Box3,BoxDiff,pentMK
c.... Check with Paco's code
      real*8 pC(0:3,6),pp23(0:3),pp45(0:3),pp234(0:3),pCsum(0:3),m,
     &       p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
      real*8 p1sq,p1p2,p1p3,p1p4,p1p5,p1p6,p2sq,p2p3,p2p4,p2p5,p2p6,p3sq,
     &       p3p4,p3p5,p3p6,p4sq,p4p5,p4p6,p5sq,p5p6,p6sq
      complex*16 mup(0:3,5),mup4test(0:3),
     &           mup1(0:3),mup2(0:3),mup3(0:3),mup4(0:3),mup5(0:3)
      complex*16 p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p2mup2,p2mup3,
     &           p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup5,p4mup1,
     &           p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5mup4,
     &           p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,
     &           mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5
      complex*16 resulttemp1,resulttemp,result,test4g(0:3)

c... External functions
      complex*16 dotrc,dotcc
      real*8 dotrr
      EXTERNAL dotrr,dotrc,dotcc,
     &         H2G_Loop,H3G_Loop,!H4G_LoopDD2,
     &         H2G_effective,H3G_effective,
     &         GGG_vertex,
     &         H3G_LoopGC,H4G_Loop,H4G_LoopDD
c... FC
       complex*16 FDeffPenLoop
       external FDeffPenLoop

      Common/Forpent/Kpen,epsgPen
      complex*16  epsgPen(0:3,4)

      logical SelectABottom,SelectATop,SelectHBottom,SelectHTop
      common/select/SelectABottom,SelectATop,SelectHBottom,SelectHTop
      
      SelectAbottom=   .false.
      SelectATop=      .false.
      SelectHbottom=   .false.
      SelectHTop=      .false.     

ccccccccccccccccccccccccccccccccc

      verdbleSin=.false.
      verquadSin=.false.
      boxdbleSin=.false.
      boxquadSin=.false.
      pendbleSin=.false.
      penquadSin=.false.

      boxSingular(1)=.false.
      boxSingular(2)=.false.
      verSingular(1)=.false.
      verSingular(2)=.false.
      
      SingularC(1)=.false.
      SingularC(2)=.false.
         

        PenSingular(1)=.false.
        PenSingular(2)=.false.
        


      
      if (ienter.eq.0) then
         ienter = ienter+1
         vev = 2*xmw/(v(1,3)*dsqrt(8d0))
c... Common overall factor
         fac=64d0*pi*alfas**5/(vev*vev)
      endif

      RGB(1,2)=1d0
      RGB(1,3)=-8d0
      RGB(1,4)=-8d0
      RGB(1,5)=1d0
      RGB(1,6)=10d0
      RGB(2,3)=-8d0
      RGB(2,4)=10d0
      RGB(2,5)=1d0
      RGB(2,6)=-8d0
      RGB(3,4)=1d0
      RGB(3,5)=10d0
      RGB(3,6)=1d0
      RGB(4,5)=-8d0
      RGB(4,6)=1d0
      RGB(5,6)=-8d0
      do i=1,6
         RGB(i,i)=64d0
         do j=1,6
            if(i.gt.j) then
               RGB(i,j)=RGB(j,i)
            endif
         enddo
      enddo

c...  Assign internal momentum variables
      SIGMAF(1) = FSIGN(1)*FPOL(1)
      SIGMAF(2) = FSIGN(2)*FPOL(2)
      IF (SIGMAF(1).NE.SIGMAF(2)) THEN
        WRITE(*,*) ' Helicity nonconserving amplitude called in HQQGG '
        RETURN
      ENDIF
      SIGMA = SIGMAF(1)

      if(MGdebug) then
         SIGMA= -1  
      endif

c... The K1,K2 (gluon momenta) are outgoing, not ingoing !!!!
c... This is necessary for using the bra-ket subroutines
      do mu = 0,3
         P(mu,1) = FSIGN(1)*PBAR(mu,1) 
         P(mu,2) = FSIGN(2)*PBAR(mu,2) 
         K(mu,1) = GSIGN(1)*KBAR(mu,1)
         K(mu,2) = GSIGN(2)*KBAR(mu,2)
         K(mu,3) = GSIGN(3)*KBAR(mu,3)
         PH(mu)  = hsign*phbar(mu)
         p21(mu) = P(mu,2)-P(mu,1)
         p21(mu)=-(K(mu,1)+K(mu,2)+K(mu,3)+PH(mu))
         K12(mu) = K(mu,1)+K(mu,2)
         K13(mu) = K(mu,1)+K(mu,3)
         K23(mu) = K(mu,2)+K(mu,3)
         p21K1(mu) = p21(mu)+K(mu,1)
         p21K2(mu) = p21(mu)+K(mu,2)
         p21K3(mu) = p21(mu)+K(mu,3)
         p21K1H(mu) = p21K1(mu)+PH(mu)
         p21K2H(mu) = p21K2(mu)+PH(mu)
         p21K3H(mu) = p21K3(mu)+PH(mu)
         p21K12(mu) = p21(mu)+K12(mu)
         p21K13(mu) = p21(mu)+K13(mu)
         p21K23(mu) = p21(mu)+K23(mu)
         PHK12(mu) = PH(mu)+K12(mu)
         PHK13(mu) = PH(mu)+K13(mu)
         PHK23(mu) = PH(mu)+K23(mu)
         K123(mu) = K12(mu)+K(mu,3)
c... 4-dim momenta variables
         Pp(mu,1) = P(mu,1)
         Pp(mu,2) = P(mu,2)
         Kp(mu,1) = K(mu,1)
         Kp(mu,2) = K(mu,2)
         Kp(mu,3) = K(mu,3)
         PHp(mu)  = PH(mu)
         p21p(mu) = p21(mu)
         K12p(mu) = K12(mu)
         K13p(mu) = K13(mu)
         K23p(mu) = K23(mu)
         p21K1p(mu) = p21K1(mu)
         p21K2p(mu) = p21K2(mu)
         p21K3p(mu) = p21K3(mu)
         p21K1Hp(mu) = p21K1H(mu)
         p21K2Hp(mu) = p21K2H(mu)
         p21K3Hp(mu) = p21K3H(mu)
         p21K12p(mu) = p21K12(mu)
         p21K13p(mu) = p21K13(mu)
         p21K23p(mu) = p21K23(mu)
         PHK12p(mu) = PHK12(mu)
         PHK13p(mu) = PHK13(mu)
         PHK23p(mu) = PHK23(mu)
         K123p(mu) = K123(mu)
      enddo
      P(4,1)=0D0
      P(4,2)=0D0
      K(4,1)=0D0
      K(4,2)=0D0
      K(4,3)=0D0
      p21K1(4)=dotrr(p21K1p,p21K1p) 
      p21K2(4)=dotrr(p21K2p,p21K2p) 
      p21K3(4)=dotrr(p21K3p,p21K3p) 
      p21K1H(4)=dotrr(p21K1Hp,p21K1Hp)
      p21K2H(4)=dotrr(p21K2Hp,p21K2Hp)
      p21K3H(4)=dotrr(p21K3Hp,p21K3Hp)
      p21(4)=dotrr(p21p,p21p) 
      K12(4)=dotrr(K12p,K12p)
      K13(4)=dotrr(K13p,K13p)
      K23(4)=dotrr(K23p,K23p)
      p21K12(4)=dotrr(p21K12p,p21K12p)
      p21K13(4)=dotrr(p21K13p,p21K13p)
      p21K23(4)=dotrr(p21K23p,p21K23p)
      PHK12(4)=dotrr(PHK12p,PHK12p)
      PHK13(4)=dotrr(PHK13p,PHK13p)
      PHK23(4)=dotrr(PHK23p,PHK23p)
      K123(4)=dotrr(K123p,K123p)


c... get the external spinors (including factor sqrt(2E) )
      CALL psi0m(2,pbar(0,1),fsign,psi)

c... get the f-fbar current J21^mu
      call curr(1,psi(1,-1,2),psi(1,-1,1),j21aux)
      do mu=0,3
         j21(mu)=j21aux(mu,SIGMA)/p21(4) !with attached propagator !!!
      enddo

      if(MGdebug) then
         CALL vxxxxx(KBAR(0,1),0d0,-1,-1 ,epsChelas(1,-1,1))
         CALL vxxxxx(KBAR(0,1),0d0, 1,-1 ,epsChelas(1, 1,1))
         CALL vxxxxx(KBAR(0,2),0d0,-1,-1 ,epsChelas(1,-1,2))
         CALL vxxxxx(KBAR(0,2),0d0, 1,-1 ,epsChelas(1, 1,2))
         CALL vxxxxx(KBAR(0,3),0d0,-1,-1 ,epsChelas(1,-1,3))
         CALL vxxxxx(KBAR(0,3),0d0, 1,-1 ,epsChelas(1, 1,3))
         do mu=0,3
            epsgC(mu,1,1)=epsChelas(mu+1, 1,1)
            epsgC(mu,2,1)=epsChelas(mu+1,-1,1)
            epsgC(mu,1,2)=epsChelas(mu+1, 1,2)
            epsgC(mu,2,2)=epsChelas(mu+1,-1,2)
            epsgC(mu,1,3)=epsChelas(mu+1, 1,3)
            epsgC(mu,2,3)=epsChelas(mu+1,-1,3)
         enddo
         
         if(lgauge2) then
            do mu=0,3
               do i1pol=1,2
                  epsgC(mu,i1pol,1)=K(mu,1)/K(0,1)
               enddo 
            enddo
         endif
c... alternative conversion of polarization vectors to the MadGraph scheme
c      do mu=0,3
c         epsgC(mu,1,1)=1d0/dsqrt(2d0)*DCMPLX( epsg(mu,1,1),epsg(mu,2,1))
c         epsgC(mu,2,1)=1d0/dsqrt(2d0)*DCMPLX(-epsg(mu,1,1),epsg(mu,2,1))
c         epsgC(mu,1,2)=1d0/dsqrt(2d0)*DCMPLX( epsg(mu,1,2),epsg(mu,2,2))
c         epsgC(mu,2,2)=1d0/dsqrt(2d0)*DCMPLX(-epsg(mu,1,2),epsg(mu,2,2))
c         epsgC(mu,1,3)=1d0/dsqrt(2d0)*DCMPLX( epsg(mu,1,3),epsg(mu,2,3))
c         epsgC(mu,2,3)=1d0/dsqrt(2d0)*DCMPLX(-epsg(mu,1,3),epsg(mu,2,3))
c      enddo
      endif

c... get the gluon polarization vector and the gluon emission spinors
      do i1pol = 1,2	! 2 gluon polarizations
         if(.not.MGdebug) then
         CALL polvec(kbar(0,1),i1pol,epsg(0,i1pol,1)) ! get gluon pol.vectors
         CALL polvec(kbar(0,2),i1pol,epsg(0,i1pol,2)) ! get gluon pol.vectors
         CALL polvec(kbar(0,3),i1pol,epsg(0,i1pol,3)) ! get gluon pol.vectors
         do mu=0,3
            epsgC(mu,i1pol,1)=epsg(mu,i1pol,1) ! complex pol. vector
            epsgC(mu,i1pol,2)=epsg(mu,i1pol,2) ! complex pol. vector
            epsgC(mu,i1pol,3)=epsg(mu,i1pol,3) ! complex pol. vector
         enddo
         endif
         if(lgauge2) then
            do mu=0,3
               epsgC(mu,i1pol,1)=Kbar(mu,1)/K(0,1)
            enddo
         endif

         do isig = -1,1,2
            CALL ket2c(psi(1,isig,1),.false.,p(0,1),isig,K(0,1),
     &       epsgC(0,i1pol,1),braketg(1,isig,1,i1pol),DUMMY)  ! |K1,1>_i1pol,isig
            CALL bra2c(psi(1,isig,2),.false.,p(0,2),isig,K(0,1),
     &       epsgC(0,i1pol,1),braketg(1,isig,2,i1pol),DUMMY)  ! <2,K1|_i1pol,isig
            CALL ket2c(psi(1,isig,1),.false.,p(0,1),isig,K(0,2),
     &       epsgC(0,i1pol,2),braketg(1,isig,3,i1pol),DUMMY)  ! |K2,1>_i1pol,isig
            CALL bra2c(psi(1,isig,2),.false.,p(0,2),isig,K(0,2),
     &       epsgC(0,i1pol,2),braketg(1,isig,4,i1pol),DUMMY)  ! <2,K2|_i1pol,isig
            CALL ket2c(psi(1,isig,1),.false.,p(0,1),isig,K(0,3),
     &       epsgC(0,i1pol,3),braketg(1,isig,5,i1pol),DUMMY)  ! |K3,1>_i1pol,isig
            CALL bra2c(psi(1,isig,2),.false.,p(0,2),isig,K(0,3),
     &       epsgC(0,i1pol,3),braketg(1,isig,6,i1pol),DUMMY)  ! <2,K3|_i1pol,isig
         enddo
      enddo

      do i1pol = 1,2 ! gluon polarizations
         CALL curr(1,psi(1,-1,2),braketg(1,-1,1,i1pol),
     &             j21auxi(0,-1,i1pol))           !=  <2|vertex|K1,1>_i1pol
         CALL curr(1,braketg(1,-1,2,i1pol),psi(1,-1,1),
     &             j21auxf(0,-1,i1pol))           !=  <2,K1|vertex|1>_i1pol
         do mu=0,3 !with attached propagator !!!
            j21i(mu,i1pol,1)=j21auxi(mu,SIGMA,i1pol)/p21K1(4)
            j21f(mu,i1pol,1)=j21auxf(mu,SIGMA,i1pol)/p21K1(4)
         enddo
         CALL curr(1,psi(1,-1,2),braketg(1,-1,3,i1pol),
     &             j21auxi(0,-1,i1pol))           !=  <2|vertex|K2,1>_i1pol
         CALL curr(1,braketg(1,-1,4,i1pol),psi(1,-1,1),
     &             j21auxf(0,-1,i1pol))           !=  <2,K2|vertex|1>_i1pol
         do mu=0,3 !with attached propagator !!!
            j21i(mu,i1pol,2)=j21auxi(mu,SIGMA,i1pol)/p21K2(4)
            j21f(mu,i1pol,2)=j21auxf(mu,SIGMA,i1pol)/p21K2(4)
         enddo
         CALL curr(1,psi(1,-1,2),braketg(1,-1,5,i1pol),
     &             j21auxi(0,-1,i1pol))	          !=  <2|vertex|K3,1>_i1pol
         CALL curr(1,braketg(1,-1,6,i1pol),psi(1,-1,1),
     &             j21auxf(0,-1,i1pol))  	  !=  <2,K3|vertex|1>_i1pol
         do mu=0,3 !with attached propagator !!!
            j21i(mu,i1pol,3)=j21auxi(mu,SIGMA,i1pol)/p21K3(4)
            j21f(mu,i1pol,3)=j21auxf(mu,SIGMA,i1pol)/p21K3(4)
         enddo
      enddo

      do i1pol=1,2
         do i2pol=1,2
            CALL curr(1,braketg(1,-1,4,i1pol),braketg(1,-1,1,i2pol),
     &                j21if(0,-1,i1pol,i2pol,1))            != <2,K2|vertex|K1,1>
            CALL curr(1,braketg(1,-1,2,i1pol),braketg(1,-1,3,i2pol),
     &                j21if(0,-1,i1pol,i2pol,2))            != <2,K1|vertex|K2,1>
            CALL curr(1,braketg(1,-1,6,i1pol),braketg(1,-1,1,i2pol),
     &                j21if(0,-1,i1pol,i2pol,3))            != <2,K3|vertex|K1,1>
            CALL curr(1,braketg(1,-1,2,i1pol),braketg(1,-1,5,i2pol),
     &                j21if(0,-1,i1pol,i2pol,4))            != <2,K1|vertex|K3,1>
            CALL curr(1,braketg(1,-1,6,i1pol),braketg(1,-1,3,i2pol),
     &                j21if(0,-1,i1pol,i2pol,5))            != <2,K3|vertex|K2,1>
            CALL curr(1,braketg(1,-1,4,i1pol),braketg(1,-1,5,i2pol),
     &                j21if(0,-1,i1pol,i2pol,6))            != <2,K2|vertex|K3,1>
         enddo
      enddo

      do i=1,2
c... J21-1; p21K1: 49
         CALL  GGG_current(p21p,j21,Kp(0,1),epsgC(0,i,1),DUMMYp
     &                    ,GGGj211(0,i))
c... 2-J21; p21K2: 32,35
         CALL  GGG_current(Kp(0,2),epsgC(0,i,2),p21p,j21,DUMMYp
     &                    ,GGG2j21(0,i))
c... 3-J21; p21K3: 60,64
         CALL  GGG_current(Kp(0,3),epsgC(0,i,3),p21p,j21,DUMMYp
     &                    ,GGG3j21(0,i))
         do j=1,2
c... K12 : 24,41,58,61,64
            CALL  GGG_current(Kp(0,1),epsgC(0,i,1),Kp(0,2),
     &                        epsgC(0,j,2),DUMMYp,GGG12(0,i,j))
c... K13 : 8,27,31,33,35,73
            CALL  GGG_current(Kp(0,1),epsgC(0,i,1),Kp(0,3),
     &                        epsgC(0,j,3),DUMMYp,GGG13(0,i,j))
c... K23 : 16
            CALL  GGG_current(Kp(0,2),epsgC(0,i,2),Kp(0,3),
     &                        epsgC(0,j,3),DUMMYp,GGG23(0,i,j))
c... K32 : 42,48,81
            CALL  GGG_current(Kp(0,3),epsgC(0,i,3),Kp(0,2),
     &                        epsgC(0,j,2),DUMMYp,GGG32(0,i,j))
         enddo
      enddo
      do i=1,2
         do j=1,2
            do l=1,2
c... K123; K123 : 62
               CALL  GGG_current(K12p,GGG12(0,i,j),Kp(0,3),
     &                          epsgC(0,l,3),DUMMYp,GGG123(0,i,j,l))
c... K132a; K123 : 33
               CALL  GGG_current(K13p,GGG13(0,i,j),Kp(0,2),
     &                         epsgC(0,l,2),DUMMYp,GGG132a(0,i,j,l))
c... K132b; K123 : 50
               CALL  GGG_current(Kp(0,1),epsgC(0,i,1),K23p,
     &                         GGG32(0,j,l),DUMMYp,GGG132b(0,i,j,l))
            enddo
         enddo
      enddo

      do i=1,6
         color(i) =DCMPLX(0d0,0d0)
      enddo
         colorT=DCMPLX(0d0,0d0)
      do i=1,87
         Heff3j(i)=DCMPLX(0d0,0d0)
         Htop3j(i)=DCMPLX(0d0,0d0)
         Hbottom3j(i)=DCMPLX(0d0,0d0)
         Aeff3j(i)=DCMPLX(0d0,0d0)
         Atop3j(i)=DCMPLX(0d0,0d0)
         Abottom3j(i)=DCMPLX(0d0,0d0)
         ampQG(i) =DCMPLX(0d0,0d0) 
      enddo
      RH      = 0.0d0
      RHPenSq = 0.0d0
      RHpenI  = 0.0d0
      RHsum   = 0.0d0

*********************************************************************
******************** EFFECTIVE CP-EVEN HIGGS ************************
*********************************************************************
      if(Heff) then
         CALL EpsH_effective(p21p,.true.,j21,PHp,PHp21p,j21Heff)
         PHp21(4)=dotrr(PHp21p,PHp21p) ! both gluons are virtual !!!
         do i1pol=1,2           ! gluon polarisations
            CALL EpsH_effective(Kp(0,1),.false.,epsgC(0,i1pol,1),PHp,
     &           PHK1p,epsHeff(0,i1pol,1))
            PHK1(4)=dotrr(PHK1p,PHK1p)
            CALL EpsH_effective(Kp(0,2),.false.,epsgC(0,i1pol,2),PHp,
     &           PHK2p,epsHeff(0,i1pol,2))
            PHK2(4)=dotrr(PHK2p,PHK2p)
            CALL EpsH_effective(Kp(0,3),.false.,epsgC(0,i1pol,3),PHp,
     &           PHK3p,epsHeff(0,i1pol,3))
            PHK3(4)=dotrr(PHK3p,PHK3p)
c... 4/5-dim momentum variables
            do mu=0,3
               PHp21(mu)=PHp21p(mu)
               PHK1(mu)=PHK1p(mu)
               PHK2(mu)=PHK2p(mu)
               PHK3(mu)=PHK3p(mu)
            enddo
c  with effective polarization vector 
         do isig = -1,1,2
            CALL ket2c(psi(1,isig,1),.false.,p(0,1),isig,PHK1,
     &           epsHeff(0,i1pol,1),braketg(1,isig,7,i1pol),DUMMY)  ! |H,K1,1>_i1pol,SIGMA1
            CALL bra2c(psi(1,isig,2),.false.,p(0,2),isig,PHK1,
     &           epsHeff(0,i1pol,1),braketg(1,isig,8,i1pol),DUMMY)  ! <2,H,K1|_i1pol,SIGMA2
            CALL ket2c(psi(1,isig,1),.false.,p(0,1),isig,PHK2,
     &           epsHeff(0,i1pol,2),braketg(1,isig,9,i1pol),DUMMY)  ! |H,K2,1>_2i1pol,SIGMA1
            CALL bra2c(psi(1,isig,2),.false.,p(0,2),isig,PHK2,
     &           epsHeff(0,i1pol,2),braketg(1,isig,10,i1pol),DUMMY) ! <2,H,K2|_i1pol,SIGMA2
            CALL ket2c(psi(1,isig,1),.false.,p(0,1),isig,PHK3,
     &           epsHeff(0,i1pol,3),braketg(1,isig,11,i1pol),DUMMY) ! |H,K3,1>_2i1pol,SIGMA1
            CALL bra2c(psi(1,isig,2),.false.,p(0,2),isig,PHK3,
     &           epsHeff(0,i1pol,3),braketg(1,isig,12,i1pol),DUMMY) ! <2,H,K3|_i1pol,SIGMA2
         enddo

            CALL curr(1,psi(1,-1,2),braketg(1,-1,7,i1pol),
     &                j21auxi(0,-1,i1pol))        !=  <2|vertex|H,K1,1>_i1pol,isig1
            CALL curr(1,braketg(1,-1,8,i1pol),psi(1,-1,1),
     &                j21auxf(0,-1,i1pol))        !=  <2,H,K1|vertex|1>_i1pol,isig1
            do mu=0,3 !with attached propagator !!!
               j21Heffi(mu,i1pol,1)=j21auxi(mu,SIGMA,i1pol)/p21K1H(4)
               j21Hefff(mu,i1pol,1)=j21auxf(mu,SIGMA,i1pol)/p21K1H(4)
            enddo
            CALL curr(1,psi(1,-1,2),braketg(1,-1,9,i1pol),
     &                j21auxi(0,-1,i1pol))        !=  <2|vertex|H,K2,1>_i1pol,isig1
            CALL curr(1,braketg(1,-1,10,i1pol),psi(1,-1,1),
     &                j21auxf(0,-1,i1pol))        !=  <2,H,K2|vertex|1>_i1pol,isig1
            do mu=0,3 !with attached propagator !!!
               j21Heffi(mu,i1pol,2)=j21auxi(mu,SIGMA,i1pol)/p21K2H(4)
               j21Hefff(mu,i1pol,2)=j21auxf(mu,SIGMA,i1pol)/p21K2H(4)
            enddo
            CALL curr(1,psi(1,-1,2),braketg(1,-1,11,i1pol),
     &                j21auxi(0,-1,i1pol))        !=  <2|vertex|H,K3,1>_i1pol,isig1
            CALL curr(1,braketg(1,-1,12,i1pol),psi(1,-1,1),
     &                j21auxf(0,-1,i1pol))        !=  <2,H,K3|vertex|1>_i1pol,isig1
            do mu=0,3 !with attached propagator !!!
               j21Heffi(mu,i1pol,3)=j21auxi(mu,SIGMA,i1pol)/p21K3H(4)
               j21Hefff(mu,i1pol,3)=j21auxf(mu,SIGMA,i1pol)/p21K3H(4)
            enddo
      enddo

      do i1pol=1,2              ! gluon polarizations
         do i2pol=1,2           ! gluon polarizations
c  Higgs production in initial state
           CALL curr(1,braketg(1,-1,4,i1pol),braketg(1,-1,7,i2pol),
     &           j21Heffif(0,-1,i1pol,i2pol,1))        !=  <2,K2|vertex|H,K1,1>

            CALL curr(1,braketg(1,-1,2,i1pol),braketg(1,-1,9,i2pol),
     &           j21Heffif(0,-1,i1pol,i2pol,2))        !=  <2,K1|vertex|H,K2,1>

            CALL curr(1,braketg(1,-1,2,i1pol),braketg(1,-1,11,i2pol),
     &           j21Heffif(0,-1,i1pol,i2pol,3))        !=  <2,K1|vertex|H,K3,1>

            CALL curr(1,braketg(1,-1,6,i1pol),braketg(1,-1,7,i2pol),
     &           j21Heffif(0,-1,i1pol,i2pol,4))        !=  <2,K3|vertex|H,K1,1>

            CALL curr(1,braketg(1,-1,6,i1pol),braketg(1,-1,9,i2pol),
     &           j21Heffif(0,-1,i1pol,i2pol,5))        !=  <2,K3|vertex|H,K2,1>

            CALL curr(1,braketg(1,-1,4,i1pol),braketg(1,-1,11,i2pol),
     &           j21Heffif(0,-1,i1pol,i2pol,6))        !=  <2,K2|vertex|H,K3,1>

c  Higgs production in final state
            CALL curr(1,braketg(1,-1,10,i1pol),braketg(1,-1,1,i2pol),
     &           j21iHefff(0,-1,i1pol,i2pol,1))        !=  <2,H,K2|vertex|K1,1>

            CALL curr(1,braketg(1,-1,8,i1pol),braketg(1,-1,3,i2pol),
     &           j21iHefff(0,-1,i1pol,i2pol,2))        !=  <2,H,K1|vertex|K2,1>

            CALL curr(1,braketg(1,-1,8,i1pol),braketg(1,-1,5,i2pol),
     &           j21iHefff(0,-1,i1pol,i2pol,3))        !=  <2,H,K1|vertex|K3,1>

            CALL curr(1,braketg(1,-1,12,i1pol),braketg(1,-1,1,i2pol),
     &           j21iHefff(0,-1,i1pol,i2pol,4))        !=  <2,H,K3|vertex|K1,1>

            CALL curr(1,braketg(1,-1,10,i1pol),braketg(1,-1,5,i2pol),
     &           j21iHefff(0,-1,i1pol,i2pol,5))        !=  <2,H,K2|vertex|K3,1>

            CALL curr(1,braketg(1,-1,12,i1pol),braketg(1,-1,3,i2pol),
     &           j21iHefff(0,-1,i1pol,i2pol,6))        !=  <2,H,K3|vertex|K2,1>
         enddo
      enddo
      do i=1,2
         do j=1,2
            CALL  GGG_current(PHK1p,epsHeff(0,i,1),Kp(0,2),epsgC(0,j,2),
     &                        DUMMYp,GGGPH12(0,i,j)) !54
            CALL  GGG_current(PHK1p,epsHeff(0,i,1),Kp(0,3),epsgC(0,j,3),
     &                        DUMMYp,GGGPH13(0,i,j)) !25
            CALL  GGG_current(PHK2p,epsHeff(0,i,2),Kp(0,3),epsgC(0,j,3),
     &                        DUMMYp,GGGPH23(0,i,j)) !44
            CALL  GGG_current(PHK3p,epsHeff(0,i,3),Kp(0,2),epsgC(0,j,2),
     &                        DUMMYp,GGGPH32(0,i,j)) !46
            CALL  GGG_current(Kp(0,1),epsgC(0,i,1),PHK2p,epsHeff(0,j,2),
     &                        DUMMYp,GGG1PH2(0,i,j)) !56
            CALL  GGG_current(Kp(0,1),epsgC(0,i,1),PHK3p,epsHeff(0,j,3),
     &                        DUMMYp,GGG1PH3(0,i,j)) !29
         enddo
      enddo
      endif

*********************************************************************
****************** CP-EVEN HIGGS WITH TOP-LOOP **********************
*********************************************************************
      if (HloopTop) then
         do mu=0,3
            Kpen(mu,1)=K(mu,1)
            Kpen(mu,2)=K(mu,2)
            Kpen(mu,3)=K(mu,3)
            Kpen(mu,4)=p21(mu)
         enddo
         Kpen(4,1)=K(4,1)
         Kpen(4,2)=K(4,2)
         Kpen(4,3)=K(4,3)
         Kpen(4,4)=p21(4)
         do i=1,4
            do j=i,4
               if (i.eq.j) then
                  q_dot_q(i,j)=Kpen(4,i)
               else
                  q_dot_q(i,j)=dotrr(Kpen(0,i),Kpen(0,j))
                  q_dot_q(j,i) = q_dot_q(i,j)
               endif
            enddo
         enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
cFC Is this necessary???
c
cfc         CALL calculate_E_functions() !depend on q_dot_q(i,j) only
c
cFC
cFC        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         CALL EpsH_Loop(p21,.true.,j21,PHp,PHp21p,j21Htop,m_top)
         PHp21(4)=dotrr(PHp21p,PHp21p) ! both gluons are virtual !!!
         do i1pol=1,2           ! gluon polarisations
            CALL EpsH_Loop(K(0,1),.false.,epsgC(0,i1pol,1),PHp,PHK1p,
     &                          epsHtop(0,i1pol,1),m_top)
            PHK1(4)=dotrr(PHK1p,PHK1p)
            CALL EpsH_Loop(K(0,2),.false.,epsgC(0,i1pol,2),PHp,PHK2p,
     &                          epsHtop(0,i1pol,2),m_top)
            PHK2(4)=dotrr(PHK2p,PHK2p)
            CALL EpsH_Loop(K(0,3),.false.,epsgC(0,i1pol,3),PHp,PHK3p,
     &                       epsHtop(0,i1pol,3),m_top)
            PHK3(4)=dotrr(PHK3p,PHK3p)
c... 4/5-dim momentum variables
            do mu=0,3
               PHp21(mu)=PHp21p(mu)
               PHK1(mu)=PHK1p(mu)
               PHK2(mu)=PHK2p(mu)
               PHK3(mu)=PHK3p(mu)
            enddo
c  with effective polarization vector 
         do isig = -1,1,2
            CALL ket2c(psi(1,isig,1),.false.,p(0,1),isig,PHK1,
     &           epsHtop(0,i1pol,1),braketg(1,isig,7,i1pol),DUMMY)  ! |H,K1,1>_i1pol,SIGMA1
            CALL bra2c(psi(1,isig,2),.false.,p(0,2),isig,PHK1,
     &           epsHtop(0,i1pol,1),braketg(1,isig,8,i1pol),DUMMY)  ! <2,H,K1|_i1pol,SIGMA2
            CALL ket2c(psi(1,isig,1),.false.,p(0,1),isig,PHK2,
     &           epsHtop(0,i1pol,2),braketg(1,isig,9,i1pol),DUMMY)  ! |H,K2,1>_2i1pol,SIGMA1
            CALL bra2c(psi(1,isig,2),.false.,p(0,2),isig,PHK2,
     &           epsHtop(0,i1pol,2),braketg(1,isig,10,i1pol),DUMMY) ! <2,H,K2|_i1pol,SIGMA2
            CALL ket2c(psi(1,isig,1),.false.,p(0,1),isig,PHK3,
     &           epsHtop(0,i1pol,3),braketg(1,isig,11,i1pol),DUMMY) ! |H,K3,1>_2i1pol,SIGMA1
            CALL bra2c(psi(1,isig,2),.false.,p(0,2),isig,PHK3,
     &           epsHtop(0,i1pol,3),braketg(1,isig,12,i1pol),DUMMY) ! <2,H,K3|_i1pol,SIGMA2
         enddo

            CALL curr(1,psi(1,-1,2),braketg(1,-1,7,i1pol),
     &                j21auxi(0,-1,i1pol))        !=  <2|vertex|H,K1,1>_i1pol,isig1
            CALL curr(1,braketg(1,-1,8,i1pol),psi(1,-1,1),
     &                j21auxf(0,-1,i1pol))        !=  <2,H,K1|vertex|1>_i1pol,isig1
            do mu=0,3 !with attached propagator !!!
               j21Htopi(mu,i1pol,1)=j21auxi(mu,SIGMA,i1pol)/p21K1H(4)
               j21Htopf(mu,i1pol,1)=j21auxf(mu,SIGMA,i1pol)/p21K1H(4)
            enddo
            CALL curr(1,psi(1,-1,2),braketg(1,-1,9,i1pol),
     &                j21auxi(0,-1,i1pol))        !=  <2|vertex|H,K2,1>_i1pol,isig1
            CALL curr(1,braketg(1,-1,10,i1pol),psi(1,-1,1),
     &                j21auxf(0,-1,i1pol))        !=  <2,H,K2|vertex|1>_i1pol,isig1
            do mu=0,3 !with attached propagator !!!
               j21Htopi(mu,i1pol,2)=j21auxi(mu,SIGMA,i1pol)/p21K2H(4)
               j21Htopf(mu,i1pol,2)=j21auxf(mu,SIGMA,i1pol)/p21K2H(4)
            enddo
            CALL curr(1,psi(1,-1,2),braketg(1,-1,11,i1pol),
     &                j21auxi(0,-1,i1pol))        !=  <2|vertex|H,K3,1>_i1pol,isig1
            CALL curr(1,braketg(1,-1,12,i1pol),psi(1,-1,1),
     &                j21auxf(0,-1,i1pol))        !=  <2,H,K3|vertex|1>_i1pol,isig1
            do mu=0,3 !with attached propagator !!!
               j21Htopi(mu,i1pol,3)=j21auxi(mu,SIGMA,i1pol)/p21K3H(4)
               j21Htopf(mu,i1pol,3)=j21auxf(mu,SIGMA,i1pol)/p21K3H(4)
            enddo
      enddo

      do i1pol=1,2              ! gluon polarizations
         do i2pol=1,2           ! gluon polarizations
c  Higgs production in initial state
            CALL curr(1,braketg(1,-1,4,i1pol),braketg(1,-1,7,i2pol),
     &           j21Htopif(0,-1,i1pol,i2pol,1))        !=  <2,K2|vertex|H,K1,1>

            CALL curr(1,braketg(1,-1,2,i1pol),braketg(1,-1,9,i2pol),
     &           j21Htopif(0,-1,i1pol,i2pol,2))        !=  <2,K1|vertex|H,K2,1>

            CALL curr(1,braketg(1,-1,2,i1pol),braketg(1,-1,11,i2pol),
     &           j21Htopif(0,-1,i1pol,i2pol,3))        !=  <2,K1|vertex|H,K3,1>

            CALL curr(1,braketg(1,-1,6,i1pol),braketg(1,-1,7,i2pol),
     &           j21Htopif(0,-1,i1pol,i2pol,4))        !=  <2,K3|vertex|H,K1,1>

            CALL curr(1,braketg(1,-1,6,i1pol),braketg(1,-1,9,i2pol),
     &           j21Htopif(0,-1,i1pol,i2pol,5))        !=  <2,K3|vertex|H,K2,1>

            CALL curr(1,braketg(1,-1,4,i1pol),braketg(1,-1,11,i2pol),
     &           j21Htopif(0,-1,i1pol,i2pol,6))        !=  <2,K2|vertex|H,K3,1>

c  Higgs production in final state
            CALL curr(1,braketg(1,-1,10,i1pol),braketg(1,-1,1,i2pol),
     &           j21iHtopf(0,-1,i1pol,i2pol,1))        !=  <2,H,K2|vertex|K1,1>

            CALL curr(1,braketg(1,-1,8,i1pol),braketg(1,-1,3,i2pol),
     &           j21iHtopf(0,-1,i1pol,i2pol,2))        !=  <2,H,K1|vertex|K2,1>

            CALL curr(1,braketg(1,-1,8,i1pol),braketg(1,-1,5,i2pol),
     &           j21iHtopf(0,-1,i1pol,i2pol,3))        !=  <2,H,K1|vertex|K3,1>

            CALL curr(1,braketg(1,-1,12,i1pol),braketg(1,-1,1,i2pol),
     &           j21iHtopf(0,-1,i1pol,i2pol,4))        !=  <2,H,K3|vertex|K1,1>

            CALL curr(1,braketg(1,-1,10,i1pol),braketg(1,-1,5,i2pol),
     &           j21iHtopf(0,-1,i1pol,i2pol,5))        !=  <2,H,K2|vertex|K3,1>

            CALL curr(1,braketg(1,-1,12,i1pol),braketg(1,-1,3,i2pol),
     &           j21iHtopf(0,-1,i1pol,i2pol,6))        !=  <2,H,K3|vertex|K2,1>
         enddo
      enddo
      do i=1,2
         do j=1,2
            CALL  GGG_current(PHK1p,epsHtop(0,i,1),Kp(0,2),epsgC(0,j,2),
     &                        DUMMYp,GGGPHtop12(0,i,j)) !54
            CALL  GGG_current(PHK1p,epsHtop(0,i,1),Kp(0,3),epsgC(0,j,3),
     &                        DUMMYp,GGGPHtop13(0,i,j)) !25
            CALL  GGG_current(PHK2p,epsHtop(0,i,2),Kp(0,3),epsgC(0,j,3),
     &                        DUMMYp,GGGPHtop23(0,i,j)) !44
            CALL  GGG_current(PHK3p,epsHtop(0,i,3),Kp(0,2),epsgC(0,j,2),
     &                        DUMMYp,GGGPHtop32(0,i,j)) !46
            CALL  GGG_current(Kp(0,1),epsgC(0,i,1),PHK2p,epsHtop(0,j,2),
     &                        DUMMYp,GGG1PHtop2(0,i,j)) !56
            CALL  GGG_current(Kp(0,1),epsgC(0,i,1),PHK3p,epsHtop(0,j,3),
     &                        DUMMYp,GGG1PHtop3(0,i,j)) !29
         enddo
      enddo
      endif

c... MGdebug
      if(MGdebug) then
         i1polmax=1
         i2polmax=1
         i3polmax=1
      else
         i1polmax=2
         i2polmax=2
         i3polmax=2
      endif

      do i1pol=1,i1polmax
      do i2pol=1,i2polmax
      do i3pol=1,i3polmax
*********************************************************************
******************** EFFECTIVE CP-EVEN HIGGS ************************
*********************************************************************
         if(Heff) then
         Heff3j(1)=cplModHeff                                       !1
     &        *dotcc(j21iHefff(0,SIGMA,i1pol,i3pol,3),epsgC(0,i2pol,2))
*----------------------------------------------------------------------
         Heff3j(2)=cplModHeff                                       !2
     &        *dotcc(j21if(0,SIGMA,i1pol,i3pol,4),epsHeff(0,i2pol,2)) 
*----------------------------------------------------------------------         
         Heff3j(3)=cplModHeff                                       !3
     &        *dotcc(j21Heffif(0,SIGMA,i1pol,i3pol,3),epsgC(0,i2pol,2))
*----------------------------------------------------------------------
         Heff3j(4)=cplModHeff*GGG_vertex(Kp(0,3),p21K2p,PHK1p,         !4 
     &        epsgC(0,i3pol,3),j21i(0,i2pol,2),epsHeff(0,i1pol,1))
*----------------------------------------------------------------------
         Heff3j(5)=-cplModHeff*GGG_vertex(Kp(0,3),p21K2Hp,Kp(0,1),     !5 minus
     &        epsgC(0,i3pol,3),j21Heffi(0,i2pol,2),epsgC(0,i1pol,1))
*----------------------------------------------------------------------
         Heff3j(6)=cplModHeff*GGG_vertex(PHK3p,p21K2p,Kp(0,1),         !6 
     &        epsHeff(0,i3pol,3),j21i(0,i2pol,2),epsgC(0,i1pol,1))
*----------------------------------------------------------------------
         Heff3j(7)=-cplModHeff                                      
     &      *H3G_effective(Kp(0,3),p21K2p,Kp(0,1),epsgC(0,i3pol,3),    !7
     &                             j21i(0,i2pol,2),epsgC(0,i1pol,1))
*----------------------------------------------------------------------
         Heff3j(8)=cplModHeff                                       !8 
     &    *H2G_effective(K13p,p21K2p,GGG13(0,i1pol,i3pol),j21i(0,i2pol,2))
*----------------------------------------------------------------------
         Heff3j(9)=cplModHeff                                       !9
     &     *dotcc(j21if(0,SIGMA,i3pol,i2pol,5),epsHeff(0,i1pol,1)) 
*----------------------------------------------------------------------
         Heff3j(10)=cplModHeff                                      !10
     &     *dotcc(j21Heffif(0,SIGMA,i3pol,i2pol,5),epsgC(0,i1pol,1))
*----------------------------------------------------------------------
         Heff3j(11)=cplModHeff                                      !11
     &     *dotcc(j21iHefff(0,SIGMA,i3pol,i2pol,6),epsgC(0,i1pol,1))
*----------------------------------------------------------------------
         Heff3j(12)=-cplModHeff*GGG_vertex(p21K1Hp,Kp(0,2),Kp(0,3),    !12 minus
     &        j21Hefff(0,i1pol,1),epsgC(0,i2pol,2),epsgC(0,i3pol,3))
*----------------------------------------------------------------------
         Heff3j(13)=-cplModHeff*GGG_vertex(p21K1p,Kp(0,3),PHK2p,       !13 minus
     &        j21f(0,i1pol,1),epsgC(0,i3pol,3),epsHeff(0,i2pol,2))
*----------------------------------------------------------------------
         Heff3j(14)=cplModHeff*GGG_vertex(p21K1p,Kp(0,2),PHK3p,        !14  
     &        j21f(0,i1pol,1),epsgC(0,i2pol,2),epsHeff(0,i3pol,3))
*----------------------------------------------------------------------
         Heff3j(15)=cplModHeff                                      !15
     &              *H3G_effective(p21K1p,Kp(0,3),Kp(0,2),j21f(0,i1pol,1),
     &                             epsgC(0,i3pol,3),epsgC(0,i2pol,2))
*----------------------------------------------------------------------
         Heff3j(16)=cplModHeff                                      !16 
     &    *H2G_effective(p21K1p,K23p,j21f(0,i1pol,1),GGG23(0,i2pol,i3pol))
*----------------------------------------------------------------------
         Heff3j(17)=cplModHeff                                     !17
     &               *dotcc(j21iHefff(0,SIGMA,i1pol,i2pol,2),epsgC(0,i3pol,3))
*----------------------------------------------------------------------
         Heff3j(18)=cplModHeff                                     !18
     &               *dotcc(j21Heffif(0,SIGMA,i1pol,i2pol,2),epsgC(0,i3pol,3))
*----------------------------------------------------------------------
         Heff3j(19)=cplModHeff                                     !19
     &               *dotcc(j21if(0,SIGMA,i1pol,i2pol,2),epsHeff(0,i3pol,3)) 
*----------------------------------------------------------------------
         Heff3j(20)=cplModHeff*GGG_vertex(p21K3p,PHK1p,Kp(0,2),       !20 +
     &                j21i(0,i3pol,3),epsHeff(0,i1pol,1),epsgC(0,i2pol,2))
*----------------------------------------------------------------------
         Heff3j(21)=cplModHeff*GGG_vertex(p21K3p,Kp(0,1),PHK2p,       !21 +
     &                j21i(0,i3pol,3),epsgC(0,i1pol,1),epsHeff(0,i2pol,2))
*----------------------------------------------------------------------
         Heff3j(22)=-cplModHeff*GGG_vertex(p21K3Hp,Kp(0,1),Kp(0,2),   !22 minus
     &                j21Heffi(0,i3pol,3),epsgC(0,i1pol,1),epsgC(0,i2pol,2))
*----------------------------------------------------------------------
         Heff3j(23)=-cplModHeff                                    !23 minus
     &              *H3G_effective(p21K3p,Kp(0,1),Kp(0,2),j21i(0,i3pol,3),
     &                             epsgC(0,i1pol,1),epsgC(0,i2pol,2)) 
*----------------------------------------------------------------------
         Heff3j(24)=cplModHeff                                     !24
     &           *H2G_effective(p21K3p,K12p,j21i(0,i3pol,3),GGG12(0,i1pol,i2pol))
*----------------------------------------------------------------------
         Heff3j(25)=-cplModHeff*GGG_vertex(Kp(0,2),p21p,PHK13p,       !25 minus
     &                epsgC(0,i2pol,2),j21,GGGPH13(0,i1pol,i3pol))
*----------------------------------------------------------------------
         Heff3j(26)=-cplModHeff*H4G_effective2(epsHeff(0,i1pol,1), !26 minus
     &                 epsgC(0,i3pol,3),epsgC(0,i2pol,2),j21)
*----------------------------------------------------------------------
         Heff3j(27)=cplModHeff*GGG_vertex(PHK2p,p21p,K13p,            !27  
     &                epsHeff(0,i2pol,2),j21,GGG13(0,i1pol,i3pol))
*----------------------------------------------------------------------
         Heff3j(28)=cplModHeff*H4G_effective2(epsgC(0,i1pol,1),    !28
     &                 epsgC(0,i3pol,3),epsHeff(0,i2pol,2),j21)
*----------------------------------------------------------------------
         Heff3j(29)=-cplModHeff*GGG_vertex(Kp(0,2),p21p,PHK13p,       !29 minus
     &                epsgC(0,i2pol,2),j21,GGG1PH3(0,i1pol,i3pol))
*----------------------------------------------------------------------
         Heff3j(30)=cplModHeff*H4G_effective2(epsgC(0,i1pol,1),    !30
     &                 epsHeff(0,i3pol,3),epsgC(0,i2pol,2),j21)  
*----------------------------------------------------------------------
         Heff3j(31)=-cplModHeff                                    
     &              *H3G_effective(Kp(0,2),p21p,K13p,epsgC(0,i2pol,2),!31 minus
     &                             j21,GGG13(0,i1pol,i3pol)) 
*----------------------------------------------------------------------
         Heff3j(32)=cplModHeff                                    
     &              *H3G_effective(Kp(0,3),p21K2p,Kp(0,1),            !32
     &                epsgC(0,i3pol,3),GGG2j21(0,i2pol),epsgC(0,i1pol,1)) 
*----------------------------------------------------------------------
         Heff3j(33)=cplModHeff                                     !33 VZ
     &           *H2G_effective(p21p,K123p,j21,GGG132a(0,i1pol,i3pol,i2pol))
*----------------------------------------------------------------------
         Heff3j(34)=-cplModHeff*H4G_effective2(epsgC(0,i1pol,1),
     &                  epsgC(0,i3pol,3),epsgC(0,i2pol,2),j21Heff) !34 minus
*----------------------------------------------------------------------
         Heff3j(35)=-cplModHeff                                    !35
     &     *H2G_effective(p21K2p,K13p,GGG2j21(0,i2pol),GGG13(0,i1pol,i3pol))
*----------------------------------------------------------------------
         Heff3j(37)=cplModHeff*GGG_vertex(p21K3p,PHK1p,Kp(0,2),       !37 
     &           j21f(0,i3pol,3),epsHeff(0,i1pol,1),epsgC(0,i2pol,2))
*----------------------------------------------------------------------
         Heff3j(38)=cplModHeff*GGG_vertex(p21K3p,Kp(0,1),PHK2p,       !38 
     &           j21f(0,i3pol,3),epsgC(0,i1pol,1),epsHeff(0,i2pol,2))
*----------------------------------------------------------------------
         Heff3j(39)=-cplModHeff*GGG_vertex(p21K3Hp,Kp(0,1),Kp(0,2),   !39 minus
     &           j21Hefff(0,i3pol,3),epsgC(0,i1pol,1),epsgC(0,i2pol,2))
*----------------------------------------------------------------------
         Heff3j(40)=-cplModHeff                                    !40 minus
     &              *H3G_effective(Kp(0,2),p21K3p,Kp(0,1),epsgC(0,i2pol,2),
     &                j21f(0,i3pol,3),epsgC(0,i1pol,1)) 
*----------------------------------------------------------------------
         Heff3j(41)=cplModHeff                                     !41
     &    *H2G_effective(p21K3p,K12p,j21f(0,i3pol,3),GGG12(0,i1pol,i2pol))
*----------------------------------------------------------------------
         Heff3j(42)=-cplModHeff*GGG_vertex(p21p,PHK1p,K23p,           !42 +
     &                j21,epsHeff(0,i1pol,1),GGG32(0,i3pol,i2pol))
*----------------------------------------------------------------------
         Heff3j(43)=cplModHeff*H4G_effective2(epsHeff(0,i1pol,1),  !43
     &                       j21,epsgC(0,i2pol,2),epsgC(0,i3pol,3)) 
*----------------------------------------------------------------------
         Heff3j(44)=-cplModHeff*GGG_vertex(p21p,Kp(0,1),PHK23p,       !44 +
     &                j21,epsgC(0,i1pol,1),GGGPH23(0,i2pol,i3pol))
*----------------------------------------------------------------------
         Heff3j(45)=cplModHeff*H4G_effective2(epsgC(0,i1pol,1),    !45 
     &                       j21,epsHeff(0,i2pol,2),epsgC(0,i3pol,3)) 
*----------------------------------------------------------------------
         Heff3j(46)=cplModHeff*GGG_vertex(p21p,Kp(0,1),PHK23p,        !46 
     &                j21,epsgC(0,i1pol,1),GGGPH32(0,i3pol,i2pol))
*----------------------------------------------------------------------
         Heff3j(47)=-cplModHeff*H4G_effective2(epsgC(0,i1pol,1),   !47 minus
     &                       j21,epsgC(0,i2pol,2),epsHeff(0,i3pol,3)) 
*----------------------------------------------------------------------
         Heff3j(48)=cplModHeff                                  
     &              *H3G_effective(p21p,Kp(0,1),K23p,                 !48
     &                    j21,epsgC(0,i1pol,1),GGG32(0,i3pol,i2pol)) 
*----------------------------------------------------------------------
         Heff3j(49)=-cplModHeff                                    !49 minus
     &              *H3G_effective(Kp(0,3),Kp(0,2),p21K1p,
     &                epsgC(0,i3pol,3),epsgC(0,i2pol,2),GGGj211(0,i1pol)) 
*----------------------------------------------------------------------
         Heff3j(50)=-cplModHeff                                    !50 minus
     &     *H2G_effective(p21p,K123p,j21,GGG132b(0,i1pol,i3pol,i2pol))
*----------------------------------------------------------------------
         Heff3j(51)=-cplModHeff*H4G_effective2(epsgC(0,i1pol,1),
     &                  j21Heff,epsgC(0,i2pol,2),epsgC(0,i3pol,3)) !51 minus
*----------------------------------------------------------------------
         Heff3j(52)=cplModHeff                                     !52 
     &    *H2G_effective(p21K1p,K23p,GGGj211(0,i1pol),GGG32(0,i3pol,i2pol))
*----------------------------------------------------------------------
         Heff3j(54)=-cplModHeff*GGG_vertex(p21p,PHK12p,Kp(0,3),       !54 minus
     &                 j21,GGGPH12(0,i1pol,i2pol),epsgC(0,i3pol,3))
*----------------------------------------------------------------------
         Heff3j(55)=-cplModHeff*H4G_effective2(epsHeff(0,i1pol,1), !55 minus
     &                        epsgC(0,i2pol,2),epsgC(0,i3pol,3),j21) 
*----------------------------------------------------------------------
         Heff3j(56)=-cplModHeff*GGG_vertex(p21p,PHK12p,Kp(0,3),       !56 minus
     &                 j21,GGG1PH2(0,i1pol,i2pol),epsgC(0,i3pol,3))
*----------------------------------------------------------------------
         Heff3j(57)=cplModHeff*H4G_effective2(epsgC(0,i1pol,1),    !57
     &                        epsHeff(0,i2pol,2),epsgC(0,i3pol,3),j21) 
*---------------------------------------------------------------------- 
         Heff3j(58)=-cplModHeff*GGG_vertex(p21p,K12p,PHK3p,           !58 minus
     &                 j21,GGG12(0,i1pol,i2pol),epsHeff(0,i3pol,3))
*----------------------------------------------------------------------
         Heff3j(59)=cplModHeff*H4G_effective2(epsgC(0,i1pol,1),    !43
     &                        epsgC(0,i2pol,2),epsHeff(0,i3pol,3),j21) 
*----------------------------------------------------------------------
         Heff3j(60)=cplModHeff                                     !60
     &              *H3G_effective(Kp(0,2),p21K3p,Kp(0,1),
     &                      epsgC(0,i2pol,2),GGG3j21(0,i3pol),epsgC(0,i1pol,1)) 
*----------------------------------------------------------------------
         Heff3j(61)=cplModHeff                                     
     &              *H3G_effective(p21p,K12p,Kp(0,3),                 !61
     &                      j21,GGG12(0,i1pol,i2pol),epsgC(0,i3pol,3)) 
*----------------------------------------------------------------------
         Heff3j(62)=cplModHeff                                     !62
     &     *H2G_effective(p21p,K123p,j21,GGG123(0,i1pol,i2pol,i3pol))
*----------------------------------------------------------------------
         Heff3j(63)=-cplModHeff*H4G_effective2(epsgC(0,i1pol,1),
     &                  epsgC(0,i2pol,2),epsgC(0,i3pol,3),j21Heff) !63 minus
*----------------------------------------------------------------------
         Heff3j(64)=-cplModHeff                                    !64 minus
     &     *H2G_effective(p21K3p,K12p,GGG3j21(0,i3pol),GGG12(0,i1pol,i2pol))
*----------------------------------------------------------------------
         Heff3j(66)=cplModHeff                                     !66
     &               *dotcc(j21if(0,SIGMA,i2pol,i3pol,6),epsHeff(0,i1pol,1)) 
*----------------------------------------------------------------------
         Heff3j(67)=cplModHeff                                     !67
     &        *dotcc(j21iHefff(0,SIGMA,i2pol,i3pol,5),epsgC(0,i1pol,1)) 
*----------------------------------------------------------------------
         Heff3j(68)=cplModHeff                                     !68
     &           *dotcc(j21Heffif(0,SIGMA,i2pol,i3pol,6),epsgC(0,i1pol,1)) 
*----------------------------------------------------------------------
         Heff3j(69)=cplModHeff*GGG_vertex(p21K2p,PHK1p,Kp(0,3),       !69 
     &                j21f(0,i2pol,2),epsHeff(0,i1pol,1),epsgC(0,i3pol,3))
*----------------------------------------------------------------------
         Heff3j(70)=-cplModHeff*GGG_vertex(p21K2Hp,Kp(0,1),Kp(0,3),   !70 minus
     &           j21Hefff(0,i2pol,2),epsgC(0,i1pol,1),epsgC(0,i3pol,3))
*----------------------------------------------------------------------
         Heff3j(71)=cplModHeff*GGG_vertex(p21K2p,Kp(0,1),PHK3p,       !71 
     &                j21f(0,i2pol,2),epsgC(0,i1pol,1),epsHeff(0,i3pol,3))
*----------------------------------------------------------------------
         Heff3j(72)=-cplModHeff                                    !72 minus
     &              *H3G_effective(p21K2p,Kp(0,1),Kp(0,3),
     &                j21f(0,i2pol,2),epsgC(0,i1pol,1),epsgC(0,i3pol,3)) 
*----------------------------------------------------------------------
         Heff3j(73)=cplModHeff                                     !73
     &     *H2G_effective(p21K2p,K13p,j21f(0,i2pol,2),GGG13(0,i1pol,i3pol))
*----------------------------------------------------------------------
         Heff3j(74)=cplModHeff                                     !74
     &        *dotcc(j21Heffif(0,SIGMA,i3pol,i1pol,4),epsgC(0,i2pol,2))
*----------------------------------------------------------------------
         Heff3j(75)=cplModHeff                                     !75
     &        *dotcc(j21if(0,SIGMA,i3pol,i1pol,3),epsHeff(0,i2pol,2))
*----------------------------------------------------------------------
         Heff3j(76)=cplModHeff                                     !76
     &        *dotcc(j21iHefff(0,SIGMA,i3pol,i1pol,4),epsgC(0,i2pol,2))
*----------------------------------------------------------------------
         Heff3j(77)=-cplModHeff*GGG_vertex(Kp(0,2),Kp(0,3),p21K1Hp,   !77 minus
     &                epsgC(0,i2pol,2),epsgC(0,i3pol,3),j21Heffi(0,i1pol,1))
*----------------------------------------------------------------------
         Heff3j(78)=cplModHeff*GGG_vertex(PHK2p,Kp(0,3),p21K1p,       !78
     &                epsHeff(0,i2pol,2),epsgC(0,i3pol,3),j21i(0,i1pol,1))
*----------------------------------------------------------------------
         Heff3j(79)=-cplModHeff*GGG_vertex(Kp(0,2),p21K1p,PHK3p,      !79 minus
     &                epsgC(0,i2pol,2),j21i(0,i1pol,1),epsHeff(0,i3pol,3))
*----------------------------------------------------------------------
         Heff3j(80)=cplModHeff                                     !80 
     &              *H3G_effective(Kp(0,3),Kp(0,2),p21K1p,epsgC(0,i3pol,3),
     &                             epsgC(0,i2pol,2),j21i(0,i1pol,1))
*----------------------------------------------------------------------
         Heff3j(81)=-cplModHeff                                    !81 minus
     &     *H2G_effective(p21K1p,K23p,j21i(0,i1pol,1),GGG32(0,i3pol,i2pol))
*----------------------------------------------------------------------
         Heff3j(82)=cplModHeff                                     !82
     &           *dotcc(j21Heffif(0,SIGMA,i2pol,i1pol,1),epsgC(0,i3pol,3)) 
*----------------------------------------------------------------------
         Heff3j(83)=cplModHeff                                     !83
     &           *dotcc(j21iHefff(0,SIGMA,i2pol,i1pol,1),epsgC(0,i3pol,3)) 
*----------------------------------------------------------------------
         Heff3j(84)=cplModHeff                                     !84
     &           *dotcc(j21if(0,SIGMA,i2pol,i1pol,1),epsHeff(0,i3pol,3)) 
*----------------------------------------------------------------------
         Heff3j(36)=cplModHeff
     &    *H4G_effective2(epsgC(0,i1pol,1),epsgC(0,i3pol,3),       !36 Pent
     &                    epsgC(0,i2pol,2),j21)
*----------------------------------------------------------------------
         Heff3j(53)=cplModHeff
     &    *H4G_effective2(j21,epsgC(0,i1pol,1),epsgC(0,i2pol,2),   !53 Pent
     &                        epsgC(0,i3pol,3))
*----------------------------------------------------------------------
         Heff3j(65)=cplModHeff
     &    *H4G_effective2(epsgC(0,i1pol,1),epsgC(0,i2pol,2),       !65 Pent
     &                    epsgC(0,i3pol,3),j21)
*----------------------------------------------------------------------

         if(MGdebug) then
            do i=1,84
               GGFLOamp2(i1pol,i2pol,i3pol,i)=Heff3j(i)*dsqrt(fac)
c               print*,'GGFLOamp2(',i,')',GGFLOamp2(i1pol,i2pol,i3pol,i)
            enddo
         endif
      endif                     ! end effevtive CP-even

*********************************************************************
****************** CP-EVEN HIGGS WITH TOP-LOOP **********************
*********************************************************************
      if (HloopTop) then
         do mu=0,3
            EPSpen(mu,1)=epsgC(mu,i1pol,1)
            EPSpen(mu,2)=epsgC(mu,i2pol,2)
            EPSpen(mu,3)=epsgC(mu,i3pol,3)
            EPSpen(mu,4)=j21(mu)
ccc NEW for pentagons 
            EPSgpen(mu,1)=epsgC(mu,i1pol,1)
            EPSgpen(mu,2)=epsgC(mu,i2pol,2)
            EPSgpen(mu,3)=epsgC(mu,i3pol,3)
            EPSgpen(mu,4)=j21(mu)
             PHc(mu)=ph(mu)
         enddo
         
         if(peneva) then
c        pendbleSin=.false.
c        penquadSin=.false.
         eva=.true.
      PenSingular(1)=.false.
      PenSingular(2)=.false.
cfc         CALL auxiliary_calculation(EPSpen(0,1),Kpen(0,1))
cfc         CALL H4G_LoopDD(4,3,2,1,m_top,qpentH_top1234(1))
cfc         CALL H4G_LoopDD(2,3,4,1,m_top,qpentH_top1234(2)) 
cfc         CALL H4G_LoopDD(2,1,4,3,m_top,qpentH_top1234(3)) 
cfc         CALL H4G_LoopDD(3,2,1,4,m_top,qpentH_top1234(4))
cfc        Htop3j(85)=cplModHtop*(qpentH_top1234(1)+qpentH_top1234(2)
cfc     &                         +qpentH_top1234(3)+qpentH_top1234(4))
cFC c
c         print*, " Htop (85) before",  Htop3j(85)
c        Htop3j(85)=Htop3j(85)/(FDeffPenLoop(4,3,2,1)+FDeffPenLoop(2,3,4,1)+
c     &  FDeffPenLoop(2,1,4,3) +   FDeffPenLoop(3,2,1,4) )

      selectHtop=.true.
         Htop3j(85)=(FDeffPenLoop(4,3,2,1)+FDeffPenLoop(2,3,4,1)+
     &  FDeffPenLoop(2,1,4,3) +   FDeffPenLoop(3,2,1,4) )
      IF(PenSingular(1)) then
         PendbleSin=.true.
      endif
      IF(PenSingular(2)) then
         PenquadSin=.true.
      endif
c         print*, " Htop (85) after",  abs(Htop3j(85)-1d0)

cfc         CALL H4G_LoopDD(2,1,3,4,m_top,qpentH_top1342(1))
cfc         CALL H4G_LoopDD(3,4,2,1,m_top,qpentH_top1342(2))
cfc         CALL H4G_LoopDD(4,2,1,3,m_top,qpentH_top1342(3))
cfc         CALL H4G_LoopDD(2,4,3,1,m_top,qpentH_top1342(4))
cfc         
cfc         Htop3j(86)=cplModHtop*(qpentH_top1342(1)+qpentH_top1342(2)
cfc     &                         +qpentH_top1342(3)+qpentH_top1342(4))
cfc
c        print*, " Htop (86) before",  Htop3j(86)
c        Htop3j(86)=Htop3j(86)/(FDeffPenLoop(2,1,3,4)+FDeffPenLoop(3,4,2,1)+
c     &               FDeffPenLoop(4,2,1,3)+FDeffPenLoop(2,4,3,1))
      PenSingular(1)=.false.
      PenSingular(2)=.false.
         Htop3j(86)=(FDeffPenLoop(2,1,3,4)+FDeffPenLoop(3,4,2,1)+
     &               FDeffPenLoop(4,2,1,3)+FDeffPenLoop(2,4,3,1))
      IF(PenSingular(1)) then
         PendbleSin=.true.
      endif
      IF(PenSingular(2)) then
         PenquadSin=.true.
      endif
c         print*, " Htop (86) after", abs( Htop3j(86)-1d0)
cfc
cfc         CALL H4G_LoopDD(2,3,1,4,m_top,qpentH_top1423(1))
cfc         CALL H4G_LoopDD(3,1,4,2,m_top,qpentH_top1423(2))
cfc         CALL H4G_LoopDD(4,2,3,1,m_top,qpentH_top1423(3))
cfc         CALL H4G_LoopDD(3,2,4,1,m_top,qpentH_top1423(4))
cfc         
cfc         Htop3j(87)=cplModHtop*(qpentH_top1423(1)+qpentH_top1423(2)
cfc     &                         +qpentH_top1423(3)+qpentH_top1423(4))
c         print*, " Htop (87) before",  Htop3j(87)

c        Htop3j(87)=Htop3j(87)/(FDeffPenLoop(2,3,1,4)+FDeffPenLoop(3,1,4,2)+
c     &               FDeffPenLoop(4,2,3,1)+FDeffPenLoop(3,2,4,1))
      PenSingular(1)=.false.
      PenSingular(2)=.false.
         Htop3j(87)=(FDeffPenLoop(2,3,1,4)+FDeffPenLoop(3,1,4,2)+
     &     FDeffPenLoop(4,2,3,1)+FDeffPenLoop(3,2,4,1))

      selectHtop=.false.

      IF(PenSingular(1)) then
          PendbleSin=.true.
      endif
      IF(PenSingular(2)) then
         PenquadSin=.true.
      endif
c         print*, " Htop (87) after",abs(  Htop3j(87)-1d0)
c         pause

         eva=.false.
         endif

         Htop3j(1)=cplModHtop                                     !1
     &               *dotcc(j21iHtopf(0,SIGMA,i1pol,i3pol,3),epsgC(0,i2pol,2))
*----------------------------------------------------------------------
         Htop3j(2)=cplModHtop                                     !2
     &               *dotcc(j21if(0,SIGMA,i1pol,i3pol,4),epsHtop(0,i2pol,2)) 
*----------------------------------------------------------------------  
         Htop3j(3)=cplModHtop                                     !3
     &               *dotcc(j21Htopif(0,SIGMA,i1pol,i3pol,3),epsgC(0,i2pol,2))
*----------------------------------------------------------------------
         Htop3j(4)=cplModHtop*GGG_vertex(Kp(0,3),p21K2p,PHK1p,       !4 
     &                epsgC(0,i3pol,3),j21i(0,i2pol,2),epsHtop(0,i1pol,1))
*----------------------------------------------------------------------
         Htop3j(5)=-cplModHtop*GGG_vertex(Kp(0,3),p21K2Hp,Kp(0,1),   !5 minus
     &                epsgC(0,i3pol,3),j21Htopi(0,i2pol,2),epsgC(0,i1pol,1))
*----------------------------------------------------------------------
         Htop3j(6)=cplModHtop*GGG_vertex(PHK3p,p21K2p,Kp(0,1),       !6 
     &                epsHtop(0,i3pol,3),j21i(0,i2pol,2),epsgC(0,i1pol,1))
*----------------------------------------------------------------------
         Htop3j(7)=-cplModHtop                                    !7
     &      *H3G_Loop(K(0,3),p21K2,K(0,1),epsgC(0,i3pol,3),
     &                             j21i(0,i2pol,2),epsgC(0,i1pol,1),m_top)
*----------------------------------------------------------------------
         Htop3j(8)=cplModHtop                                     !8 
     &    *H2G_Loop(K13,p21K2,GGG13(0,i1pol,i3pol),j21i(0,i2pol,2),m_top)
*----------------------------------------------------------------------
         Htop3j(9)=cplModHtop                                     !9
     &     *dotcc(j21if(0,SIGMA,i3pol,i2pol,5),epsHtop(0,i1pol,1)) 
*----------------------------------------------------------------------
         Htop3j(10)=cplModHtop                                    !10
     &     *dotcc(j21Htopif(0,SIGMA,i3pol,i2pol,5),epsgC(0,i1pol,1))
*----------------------------------------------------------------------
         Htop3j(11)=cplModHtop                                    !11
     &               *dotcc(j21iHtopf(0,SIGMA,i3pol,i2pol,6),epsgC(0,i1pol,1))
*----------------------------------------------------------------------
         Htop3j(12)=-cplModHtop*GGG_vertex(p21K1Hp,Kp(0,2),Kp(0,3),  !12 minus
     &                j21Htopf(0,i1pol,1),epsgC(0,i2pol,2),epsgC(0,i3pol,3))
*----------------------------------------------------------------------
         Htop3j(13)=-cplModHtop*GGG_vertex(p21K1p,Kp(0,3),PHK2p,     !13 minus
     &                j21f(0,i1pol,1),epsgC(0,i3pol,3),epsHtop(0,i2pol,2))
*----------------------------------------------------------------------
         Htop3j(14)=cplModHtop*GGG_vertex(p21K1p,Kp(0,2),PHK3p,      !14  
     &                j21f(0,i1pol,1),epsgC(0,i2pol,2),epsHtop(0,i3pol,3))
*----------------------------------------------------------------------
         Htop3j(15)=cplModHtop                                    !15
     &              *H3G_Loop(p21K1,K(0,3),K(0,2),j21f(0,i1pol,1),
     &                             epsgC(0,i3pol,3),epsgC(0,i2pol,2),m_top)
*----------------------------------------------------------------------
         Htop3j(16)=cplModHtop                                    !16 
     &    *H2G_Loop(p21K1,K23,j21f(0,i1pol,1),GGG23(0,i2pol,i3pol),m_top)
*----------------------------------------------------------------------
         Htop3j(17)=cplModHtop                                   !17
     &               *dotcc(j21iHtopf(0,SIGMA,i1pol,i2pol,2),epsgC(0,i3pol,3))
*----------------------------------------------------------------------
         Htop3j(18)=cplModHtop                                   !18
     &               *dotcc(j21Htopif(0,SIGMA,i1pol,i2pol,2),epsgC(0,i3pol,3))
*----------------------------------------------------------------------
         Htop3j(19)=cplModHtop                                   !19
     &               *dotcc(j21if(0,SIGMA,i1pol,i2pol,2),epsHtop(0,i3pol,3)) 
*----------------------------------------------------------------------
         Htop3j(20)=cplModHtop*GGG_vertex(p21K3p,PHK1p,Kp(0,2),     !20 +
     &                j21i(0,i3pol,3),epsHtop(0,i1pol,1),epsgC(0,i2pol,2))
*----------------------------------------------------------------------
         Htop3j(21)=cplModHtop*GGG_vertex(p21K3p,Kp(0,1),PHK2p,     !21 +
     &                j21i(0,i3pol,3),epsgC(0,i1pol,1),epsHtop(0,i2pol,2))
*----------------------------------------------------------------------
         Htop3j(22)=-cplModHtop*GGG_vertex(p21K3Hp,Kp(0,1),Kp(0,2), !22 minus
     &                j21Htopi(0,i3pol,3),epsgC(0,i1pol,1),epsgC(0,i2pol,2))
*----------------------------------------------------------------------
         Htop3j(23)=-cplModHtop                                  !23 minus
     &              *H3G_Loop(p21K3,K(0,1),K(0,2),j21i(0,i3pol,3),
     &                             epsgC(0,i1pol,1),epsgC(0,i2pol,2),m_top) 
*----------------------------------------------------------------------
         Htop3j(24)=cplModHtop                                   !24
     &           *H2G_Loop(p21K3,K12,j21i(0,i3pol,3),GGG12(0,i1pol,i2pol),m_top)
*----------------------------------------------------------------------
         Htop3j(25)=-cplModHtop*GGG_vertex(Kp(0,2),p21p,PHK13p,    !25 minus
     &                epsgC(0,i2pol,2),j21,GGGPHtop13(0,i1pol,i3pol))
*----------------------------------------------------------------------
         Htop3j(26)=-cplModHtop*H4G_effective2(epsHtop(0,i1pol,1),!26 minus
     &                 epsgC(0,i3pol,3),epsgC(0,i2pol,2),j21)
*----------------------------------------------------------------------
         Htop3j(27)=cplModHtop*GGG_vertex(PHK2p,p21p,K13p,         !27  
     &                epsHtop(0,i2pol,2),j21,GGG13(0,i1pol,i3pol))
*----------------------------------------------------------------------
         Htop3j(28)=cplModHtop*H4G_effective2(epsgC(0,i1pol,1),  !28
     &                 epsgC(0,i3pol,3),epsHtop(0,i2pol,2),j21)
*----------------------------------------------------------------------
         Htop3j(29)=-cplModHtop*GGG_vertex(Kp(0,2),p21p,PHK13p,    !29 minus
     &                epsgC(0,i2pol,2),j21,GGG1PHtop3(0,i1pol,i3pol))
*----------------------------------------------------------------------
         Htop3j(30)=cplModHtop*H4G_effective2(epsgC(0,i1pol,1),  !30
     &                 epsHtop(0,i3pol,3),epsgC(0,i2pol,2),j21)  
*----------------------------------------------------------------------
         Htop3j(31)=-cplModHtop                                  !31 minus
     &              *H3G_Loop(K(0,2),p21,K13,epsgC(0,i2pol,2),
     &                             j21,GGG13(0,i1pol,i3pol),m_top) 
*----------------------------------------------------------------------
         Htop3j(32)=cplModHtop                                   !32
     &        *H3G_Loop(K(0,3),p21K2,K(0,1),
     &           epsgC(0,i3pol,3),GGG2j21(0,i2pol),epsgC(0,i1pol,1),m_top) 
*----------------------------------------------------------------------
         Htop3j(33)=cplModHtop                                   !33 VZ
     &           *H2G_Loop(p21,K123,j21,GGG132a(0,i1pol,i3pol,i2pol),m_top)
*----------------------------------------------------------------------
c         Htop3j(34)=-cplModHtop*H4G_effective2(epsgC(0,i1pol,1),
c     &                  epsgC(0,i3pol,3),epsgC(0,i2pol,2),j21Htop)  !34 minus
         CALL Vertex4g(epsgC(0,i1pol,1), epsgC(0,i3pol,3),
     &        epsgC(0,i2pol,2),test4g)
         Htop3j(34)=-cplModHtop*H2G_Loop(p21,K123,j21,test4g,m_top)/(K123(4))
*----------------------------------------------------------------------
         Htop3j(35)=-cplModHtop     
     &     *H2G_Loop(p21K2,K13,GGG2j21(0,i2pol),GGG13(0,i1pol,i3pol),m_top)
*----------------------------------------------------------------------
         Htop3j(36)=DCMPLX(0d0,0d0) !---> see Pentagons (85-87)
*----------------------------------------------------------------------
         Htop3j(37)=cplModHtop*GGG_vertex(p21K3p,PHK1p,Kp(0,2),     !37 
     &           j21f(0,i3pol,3),epsHtop(0,i1pol,1),epsgC(0,i2pol,2))
*----------------------------------------------------------------------
         Htop3j(38)=cplModHtop*GGG_vertex(p21K3p,Kp(0,1),PHK2p,     !38 
     &           j21f(0,i3pol,3),epsgC(0,i1pol,1),epsHtop(0,i2pol,2))
*----------------------------------------------------------------------
         Htop3j(39)=-cplModHtop*GGG_vertex(p21K3Hp,Kp(0,1),Kp(0,2), !39 minus
     &           j21Htopf(0,i3pol,3),epsgC(0,i1pol,1),epsgC(0,i2pol,2))
*----------------------------------------------------------------------
         Htop3j(40)=-cplModHtop                                  !40 minus
     &              *H3G_Loop(K(0,2),p21K3,K(0,1),epsgC(0,i2pol,2),
     &                j21f(0,i3pol,3),epsgC(0,i1pol,1),m_top) 
*----------------------------------------------------------------------
         Htop3j(41)=cplModHtop                                   !41
     &    *H2G_Loop(p21K3,K12,j21f(0,i3pol,3),GGG12(0,i1pol,i2pol),m_top)
*----------------------------------------------------------------------
         Htop3j(42)=-cplModHtop*GGG_vertex(p21p,PHK1p,K23p,     !42 +
     &                j21,epsHtop(0,i1pol,1),GGG32(0,i3pol,i2pol))
*----------------------------------------------------------------------
         Htop3j(43)=cplModHtop*H4G_effective2(epsHtop(0,i1pol,1),!43
     &                       j21,epsgC(0,i2pol,2),epsgC(0,i3pol,3)) 
*----------------------------------------------------------------------
         Htop3j(44)=-cplModHtop*GGG_vertex(p21p,Kp(0,1),PHK23p,   !44 +
     &                j21,epsgC(0,i1pol,1),GGGPHtop23(0,i2pol,i3pol))
*----------------------------------------------------------------------
         Htop3j(45)=cplModHtop*H4G_effective2(epsgC(0,i1pol,1),!45 
     &                       j21,epsHtop(0,i2pol,2),epsgC(0,i3pol,3)) 
*----------------------------------------------------------------------
         Htop3j(46)=cplModHtop*GGG_vertex(p21p,Kp(0,1),PHK23p,    !46 
     &                j21,epsgC(0,i1pol,1),GGGPHtop32(0,i3pol,i2pol))
*----------------------------------------------------------------------
         Htop3j(47)=-cplModHtop*H4G_effective2(epsgC(0,i1pol,1),!47 minus
     &                       j21,epsgC(0,i2pol,2),epsHtop(0,i3pol,3)) 
*----------------------------------------------------------------------
         Htop3j(48)=cplModHtop                                   !48
     &              *H3G_Loop(p21,K(0,1),K23,
     &                    j21,epsgC(0,i1pol,1),GGG32(0,i3pol,i2pol),m_top) 
*----------------------------------------------------------------------
         Htop3j(49)=-cplModHtop                                  !49 minus
     &              *H3G_Loop(K(0,3),K(0,2),p21K1,
     &                epsgC(0,i3pol,3),epsgC(0,i2pol,2),GGGj211(0,i1pol),m_top) 
*----------------------------------------------------------------------
         Htop3j(50)=-cplModHtop                                  !50 minus
     &     *H2G_Loop(p21,K123,j21,GGG132b(0,i1pol,i3pol,i2pol),m_top)
*----------------------------------------------------------------------
c         Htop3j(51)=-cplModHtop*H4G_effective2(epsgC(0,i1pol,1),
c     &                  j21Htop,epsgC(0,i2pol,2),epsgC(0,i3pol,3))  !51 minus
         CALL Vertex4g(epsgC(0,i3pol,3), epsgC(0,i2pol,2),
     &        epsgC(0,i1pol,1),test4g)
         Htop3j(51)=cplModHtop*H2G_Loop(p21,K123,j21,test4g,m_top)/(K123(4))
*----------------------------------------------------------------------
         Htop3j(52)=cplModHtop                                   !52 
     &    *H2G_Loop(p21K1,K23,GGGj211(0,i1pol),GGG32(0,i3pol,i2pol),m_top)
*----------------------------------------------------------------------
         Htop3j(53)=DCMPLX(0d0,0d0) !---> see Pentagons (85-87)
*----------------------------------------------------------------------
         Htop3j(54)=-cplModHtop*GGG_vertex(p21p,PHK12p,Kp(0,3),    !54 minus
     &                 j21,GGGPHtop12(0,i1pol,i2pol),epsgC(0,i3pol,3))
*----------------------------------------------------------------------
         Htop3j(55)=-cplModHtop*H4G_effective2(epsHtop(0,i1pol,1),!55 minus
     &                        epsgC(0,i2pol,2),epsgC(0,i3pol,3),j21) 
*----------------------------------------------------------------------
         Htop3j(56)=-cplModHtop*GGG_vertex(p21p,PHK12p,Kp(0,3),    !56 minus
     &                 j21,GGG1PHtop2(0,i1pol,i2pol),epsgC(0,i3pol,3))
*----------------------------------------------------------------------
         Htop3j(57)=cplModHtop*H4G_effective2(epsgC(0,i1pol,1),!43
     &                        epsHtop(0,i2pol,2),epsgC(0,i3pol,3),j21) 
*---------------------------------------------------------------------- 
         Htop3j(58)=-cplModHtop*GGG_vertex(p21p,K12p,PHK3p,      !58 minus
     &                 j21,GGG12(0,i1pol,i2pol),epsHtop(0,i3pol,3))
*----------------------------------------------------------------------
         Htop3j(59)=cplModHtop*H4G_effective2(epsgC(0,i1pol,1),!43
     &                        epsgC(0,i2pol,2),epsHtop(0,i3pol,3),j21) 
*----------------------------------------------------------------------
         Htop3j(60)=cplModHtop                                   !60
     &              *H3G_Loop(K(0,2),p21K3,K(0,1),
     &        epsgC(0,i2pol,2),GGG3j21(0,i3pol),epsgC(0,i1pol,1),m_top) 
*----------------------------------------------------------------------
         Htop3j(61)=cplModHtop                                   !61 
     &              *H3G_Loop(p21,K12,K(0,3),
     &                      j21,GGG12(0,i1pol,i2pol),epsgC(0,i3pol,3),m_top) 
*----------------------------------------------------------------------
         Htop3j(62)=cplModHtop                                   !62
     &     *H2G_Loop(p21,K123,j21,GGG123(0,i1pol,i2pol,i3pol),m_top)
*----------------------------------------------------------------------
c         Htop3j(63)=-cplModHtop*H4G_effective2(epsgC(0,i1pol,1),
c     &                  epsgC(0,i2pol,2),epsgC(0,i3pol,3),j21Htop)  !63 minus
         CALL Vertex4g(epsgC(0,i1pol,1), epsgC(0,i2pol,2),
     &        epsgC(0,i3pol,3),test4g)
         Htop3j(63)=-cplModHtop*H2G_Loop(p21,K123,j21,test4g,m_top)/(K123(4))
*----------------------------------------------------------------------
         Htop3j(64)=-cplModHtop                                  !64 minus
     &     *H2G_Loop(p21K3,K12,GGG3j21(0,i3pol),GGG12(0,i1pol,i2pol),m_top)
*----------------------------------------------------------------------
         Htop3j(65)=DCMPLX(0d0,0d0) !---> see Pentagons (85-87)
*----------------------------------------------------------------------
         Htop3j(66)=cplModHtop                                   !66
     &               *dotcc(j21if(0,SIGMA,i2pol,i3pol,6),epsHtop(0,i1pol,1)) 
*----------------------------------------------------------------------
         Htop3j(67)=cplModHtop                                   !67
     &           *dotcc(j21iHtopf(0,SIGMA,i2pol,i3pol,5),epsgC(0,i1pol,1)) 
*----------------------------------------------------------------------
         Htop3j(68)=cplModHtop                                   !68
     &           *dotcc(j21Htopif(0,SIGMA,i2pol,i3pol,6),epsgC(0,i1pol,1)) 
*----------------------------------------------------------------------
         Htop3j(69)=cplModHtop*GGG_vertex(p21K2p,PHK1p,Kp(0,3),     !69 
     &                j21f(0,i2pol,2),epsHtop(0,i1pol,1),epsgC(0,i3pol,3))
*----------------------------------------------------------------------
         Htop3j(70)=-cplModHtop*GGG_vertex(p21K2Hp,Kp(0,1),Kp(0,3), !70 minus
     &           j21Htopf(0,i2pol,2),epsgC(0,i1pol,1),epsgC(0,i3pol,3))
*----------------------------------------------------------------------
         Htop3j(71)=cplModHtop*GGG_vertex(p21K2p,Kp(0,1),PHK3p,     !71 
     &                j21f(0,i2pol,2),epsgC(0,i1pol,1),epsHtop(0,i3pol,3))
*----------------------------------------------------------------------
         Htop3j(72)=-cplModHtop                                  !72 minus
     &              *H3G_Loop(p21K2,K(0,1),K(0,3),
     &                j21f(0,i2pol,2),epsgC(0,i1pol,1),epsgC(0,i3pol,3),m_top) 
*----------------------------------------------------------------------
         Htop3j(73)=cplModHtop                                   !73
     &     *H2G_Loop(p21K2,K13,j21f(0,i2pol,2),GGG13(0,i1pol,i3pol),m_top)
*----------------------------------------------------------------------
         Htop3j(74)=cplModHtop                                   !74
     &        *dotcc(j21Htopif(0,SIGMA,i3pol,i1pol,4),epsgC(0,i2pol,2))
*----------------------------------------------------------------------
         Htop3j(75)=cplModHtop                                   !75
     &        *dotcc(j21if(0,SIGMA,i3pol,i1pol,3),epsHtop(0,i2pol,2))
*----------------------------------------------------------------------
         Htop3j(76)=cplModHtop                                   !76
     &        *dotcc(j21iHtopf(0,SIGMA,i3pol,i1pol,4),epsgC(0,i2pol,2))
*----------------------------------------------------------------------
         Htop3j(77)=-cplModHtop*GGG_vertex(Kp(0,2),Kp(0,3),p21K1Hp, !77 minus
     &                epsgC(0,i2pol,2),epsgC(0,i3pol,3),j21Htopi(0,i1pol,1))
*----------------------------------------------------------------------
         Htop3j(78)=cplModHtop*GGG_vertex(PHK2p,Kp(0,3),p21K1p,     !78
     &                epsHtop(0,i2pol,2),epsgC(0,i3pol,3),j21i(0,i1pol,1))
*----------------------------------------------------------------------
         Htop3j(79)=-cplModHtop*GGG_vertex(Kp(0,2),p21K1p,PHK3p,    !79 minus
     &                epsgC(0,i2pol,2),j21i(0,i1pol,1),epsHtop(0,i3pol,3))
*----------------------------------------------------------------------
         Htop3j(80)=cplModHtop                                   !80 
     &              *H3G_Loop(K(0,3),K(0,2),p21K1,epsgC(0,i3pol,3),
     &                             epsgC(0,i2pol,2),j21i(0,i1pol,1),m_top)
*----------------------------------------------------------------------
         Htop3j(81)=-cplModHtop                                  !81 minus
     &     *H2G_Loop(p21K1,K23,j21i(0,i1pol,1),GGG32(0,i3pol,i2pol),m_top)
*----------------------------------------------------------------------
         Htop3j(82)=cplModHtop                                   !82
     &           *dotcc(j21Htopif(0,SIGMA,i2pol,i1pol,1),epsgC(0,i3pol,3)) 
*----------------------------------------------------------------------
         Htop3j(83)=cplModHtop                                   !83
     &           *dotcc(j21iHtopf(0,SIGMA,i2pol,i1pol,1),epsgC(0,i3pol,3)) 
*----------------------------------------------------------------------
         Htop3j(84)=cplModHtop                                   !84
     &           *dotcc(j21if(0,SIGMA,i2pol,i1pol,1),epsHtop(0,i3pol,3)) 
         
         if (MGdebug) then
            do i=1,87
               GGFLOamp2(i1pol,i2pol,i3pol,i)=Htop3j(i)*dsqrt(fac)
c               print*,'GGFLOamp2(',i,')',GGFLOamp2(i1pol,i2pol,i3pol,i)
            enddo
         endif
      endif                    ! end CP-even with top-loop

cc start check
      if (HloopBottom) then
         CALL  HQG_bottom(PBAR,FSIGN,FPOL,KBAR,GSIGN,PHBAR,HSIGN,
     &                    i1pol,i2pol,i3pol,Hbottom3j)
      endif
      if (Aeff) then
         CALL  AQG_eff(PBAR,FSIGN,FPOL,KBAR,GSIGN,PHBAR,HSIGN,
     &                 i1pol,i2pol,i3pol,Aeff3j)
      endif
      if (AloopTop) then
         CALL  AQG_top(PBAR,FSIGN,FPOL,KBAR,GSIGN,PHBAR,HSIGN,
     &                 i1pol,i2pol,i3pol,Atop3j)
      endif
      if (AloopBottom) then
         CALL  AQG_bottom(PBAR,FSIGN,FPOL,KBAR,GSIGN,PHBAR,HSIGN,
     &                    i1pol,i2pol,i3pol,Abottom3j)
      endif

cc end check

      ampQG(1) =Heff3j(1) +Htop3j(1) +Hbottom3j(1)
     &         +Aeff3j(1) +Atop3j(1) +Abottom3j(1)
      ampQG(2) =Heff3j(2) +Htop3j(2) +Hbottom3j(2)
     &         +Aeff3j(2) +Atop3j(2) +Abottom3j(2)
      ampQG(3) =Heff3j(3) +Htop3j(3) +Hbottom3j(3)
     &         +Aeff3j(3) +Atop3j(3) +Abottom3j(3)
      ampQG(4) =Heff3j(4) +Htop3j(4) +Hbottom3j(4)
     &         +Aeff3j(4) +Atop3j(4) +Abottom3j(4)
      ampQG(5) =Heff3j(5) +Htop3j(5) +Hbottom3j(5)
     &         +Aeff3j(5) +Atop3j(5) +Abottom3j(5)
      ampQG(6) =Heff3j(6) +Htop3j(6) +Hbottom3j(6)
     &         +Aeff3j(6) +Atop3j(6) +Abottom3j(6)
      ampQG(7) =Heff3j(7) +Htop3j(7) +Hbottom3j(7)
     &         +Aeff3j(7) +Atop3j(7) +Abottom3j(7)
      ampQG(8) =Heff3j(8) +Htop3j(8) +Hbottom3j(8)
     &         +Aeff3j(8) +Atop3j(8) +Abottom3j(8)
      ampQG(9) =Heff3j(9) +Htop3j(9) +Hbottom3j(9)
     &         +Aeff3j(9) +Atop3j(9) +Abottom3j(9)
      ampQG(10)=Heff3j(10)+Htop3j(10)+Hbottom3j(10)
     &         +Aeff3j(10)+Atop3j(10)+Abottom3j(10)
                                                                
      ampQG(11)=Heff3j(11)+Htop3j(11)+Hbottom3j(11)
     &         +Aeff3j(11)+Atop3j(11)+Abottom3j(11)
      ampQG(12)=Heff3j(12)+Htop3j(12)+Hbottom3j(12)
     &         +Aeff3j(12)+Atop3j(12)+Abottom3j(12)
      ampQG(13)=Heff3j(13)+Htop3j(13)+Hbottom3j(13)
     &         +Aeff3j(13)+Atop3j(13)+Abottom3j(13)
      ampQG(14)=Heff3j(14)+Htop3j(14)+Hbottom3j(14)
     &         +Aeff3j(14)+Atop3j(14)+Abottom3j(14)
      ampQG(15)=Heff3j(15)+Htop3j(15)+Hbottom3j(15)
     &         +Aeff3j(15)+Atop3j(15)+Abottom3j(15)
      ampQG(16)=Heff3j(16)+Htop3j(16)+Hbottom3j(16)
     &         +Aeff3j(16)+Atop3j(16)+Abottom3j(16)
      ampQG(17)=Heff3j(17)+Htop3j(17)+Hbottom3j(17)
     &         +Aeff3j(17)+Atop3j(17)+Abottom3j(17)
      ampQG(18)=Heff3j(18)+Htop3j(18)+Hbottom3j(18)
     &         +Aeff3j(18)+Atop3j(18)+Abottom3j(18)
      ampQG(19)=Heff3j(19)+Htop3j(19)+Hbottom3j(19)
     &         +Aeff3j(19)+Atop3j(19)+Abottom3j(19)
                                                                
      ampQG(20)=Heff3j(20)+Htop3j(20)+Hbottom3j(20)
     &         +Aeff3j(20)+Atop3j(20)+Abottom3j(20)
      ampQG(21)=Heff3j(21)+Htop3j(21)+Hbottom3j(21)
     &         +Aeff3j(21)+Atop3j(21)+Abottom3j(21)
      ampQG(22)=Heff3j(22)+Htop3j(22)+Hbottom3j(22)
     &         +Aeff3j(22)+Atop3j(22)+Abottom3j(22)
      ampQG(23)=Heff3j(23)+Htop3j(23)+Hbottom3j(23)
     &         +Aeff3j(23)+Atop3j(23)+Abottom3j(23)
      ampQG(24)=Heff3j(24)+Htop3j(24)+Hbottom3j(24)
     &         +Aeff3j(24)+Atop3j(24)+Abottom3j(24)
      ampQG(25)=Heff3j(25)+Htop3j(25)+Hbottom3j(25)
     &         +Aeff3j(25)+Atop3j(25)+Abottom3j(25)
      ampQG(26)=Heff3j(26)+Htop3j(26)+Hbottom3j(26)
     &         +Aeff3j(26)+Atop3j(26)+Abottom3j(26)
      ampQG(27)=Heff3j(27)+Htop3j(27)+Hbottom3j(27)
     &         +Aeff3j(27)+Atop3j(27)+Abottom3j(27)
      ampQG(28)=Heff3j(28)+Htop3j(28)+Hbottom3j(28)
     &         +Aeff3j(28)+Atop3j(28)+Abottom3j(28)
      ampQG(29)=Heff3j(29)+Htop3j(29)+Hbottom3j(29)
     &         +Aeff3j(29)+Atop3j(29)+Abottom3j(29)
                                                                
      ampQG(30)=Heff3j(30)+Htop3j(30)+Hbottom3j(30)
     &         +Aeff3j(30)+Atop3j(30)+Abottom3j(30)
      ampQG(31)=Heff3j(31)+Htop3j(31)+Hbottom3j(31)
     &         +Aeff3j(31)+Atop3j(31)+Abottom3j(31)
      ampQG(32)=Heff3j(32)+Htop3j(32)+Hbottom3j(32)
     &         +Aeff3j(32)+Atop3j(32)+Abottom3j(32)
      ampQG(33)=Heff3j(33)+Htop3j(33)+Hbottom3j(33)
     &         +Aeff3j(33)+Atop3j(33)+Abottom3j(33)
      ampQG(34)=Heff3j(34)+Htop3j(34)+Hbottom3j(34)
     &         +Aeff3j(34)+Atop3j(34)+Abottom3j(34)
      ampQG(35)=Heff3j(35)+Htop3j(35)+Hbottom3j(35)
     &         +Aeff3j(35)+Atop3j(35)+Abottom3j(35)
      ampQG(36)=Heff3j(36)+Htop3j(36)+Hbottom3j(36)
     &         +Aeff3j(36)+Atop3j(36)+Abottom3j(36)
      ampQG(37)=Heff3j(37)+Htop3j(37)+Hbottom3j(37)
     &         +Aeff3j(37)+Atop3j(37)+Abottom3j(37)
      ampQG(38)=Heff3j(38)+Htop3j(38)+Hbottom3j(38)
     &         +Aeff3j(38)+Atop3j(38)+Abottom3j(38)
      ampQG(39)=Heff3j(39)+Htop3j(39)+Hbottom3j(39)
     &         +Aeff3j(39)+Atop3j(39)+Abottom3j(39)
                                                                
      ampQG(40)=Heff3j(40)+Htop3j(40)+Hbottom3j(40)
     &         +Aeff3j(40)+Atop3j(40)+Abottom3j(40)
      ampQG(41)=Heff3j(41)+Htop3j(41)+Hbottom3j(41)
     &         +Aeff3j(41)+Atop3j(41)+Abottom3j(41)
      ampQG(42)=Heff3j(42)+Htop3j(42)+Hbottom3j(42)
     &         +Aeff3j(42)+Atop3j(42)+Abottom3j(42)
      ampQG(43)=Heff3j(43)+Htop3j(43)+Hbottom3j(43)
     &         +Aeff3j(43)+Atop3j(43)+Abottom3j(43)
      ampQG(44)=Heff3j(44)+Htop3j(44)+Hbottom3j(44)
     &         +Aeff3j(44)+Atop3j(44)+Abottom3j(44)
      ampQG(45)=Heff3j(45)+Htop3j(45)+Hbottom3j(45)
     &         +Aeff3j(45)+Atop3j(45)+Abottom3j(45)
      ampQG(46)=Heff3j(46)+Htop3j(46)+Hbottom3j(46)
     &         +Aeff3j(46)+Atop3j(46)+Abottom3j(46)
      ampQG(47)=Heff3j(47)+Htop3j(47)+Hbottom3j(47)
     &         +Aeff3j(47)+Atop3j(47)+Abottom3j(47)
      ampQG(48)=Heff3j(48)+Htop3j(48)+Hbottom3j(48)
     &         +Aeff3j(48)+Atop3j(48)+Abottom3j(48)
      ampQG(49)=Heff3j(49)+Htop3j(49)+Hbottom3j(49)
     &         +Aeff3j(49)+Atop3j(49)+Abottom3j(49)
                                                                
      ampQG(50)=Heff3j(50)+Htop3j(50)+Hbottom3j(50)
     &         +Aeff3j(50)+Atop3j(50)+Abottom3j(50)
      ampQG(51)=Heff3j(51)+Htop3j(51)+Hbottom3j(51)
     &         +Aeff3j(51)+Atop3j(51)+Abottom3j(51)
      ampQG(52)=Heff3j(52)+Htop3j(52)+Hbottom3j(52)
     &         +Aeff3j(52)+Atop3j(52)+Abottom3j(52)
      ampQG(53)=Heff3j(53)+Htop3j(53)+Hbottom3j(53)
     &         +Aeff3j(53)+Atop3j(53)+Abottom3j(53)
      ampQG(54)=Heff3j(54)+Htop3j(54)+Hbottom3j(54)
     &         +Aeff3j(54)+Atop3j(54)+Abottom3j(54)
      ampQG(55)=Heff3j(55)+Htop3j(55)+Hbottom3j(55)
     &         +Aeff3j(55)+Atop3j(55)+Abottom3j(55)
      ampQG(56)=Heff3j(56)+Htop3j(56)+Hbottom3j(56)
     &         +Aeff3j(56)+Atop3j(56)+Abottom3j(56)
      ampQG(57)=Heff3j(57)+Htop3j(57)+Hbottom3j(57)
     &         +Aeff3j(57)+Atop3j(57)+Abottom3j(57)
      ampQG(58)=Heff3j(58)+Htop3j(58)+Hbottom3j(58)
     &         +Aeff3j(58)+Atop3j(58)+Abottom3j(58)
      ampQG(59)=Heff3j(59)+Htop3j(59)+Hbottom3j(59)
     &         +Aeff3j(59)+Atop3j(59)+Abottom3j(59)
                                                                
      ampQG(60)=Heff3j(60)+Htop3j(60)+Hbottom3j(60)
     &         +Aeff3j(60)+Atop3j(60)+Abottom3j(60)
      ampQG(61)=Heff3j(61)+Htop3j(61)+Hbottom3j(61)
     &         +Aeff3j(61)+Atop3j(61)+Abottom3j(61)
      ampQG(62)=Heff3j(62)+Htop3j(62)+Hbottom3j(62)
     &         +Aeff3j(62)+Atop3j(62)+Abottom3j(62)
      ampQG(63)=Heff3j(63)+Htop3j(63)+Hbottom3j(63)
     &         +Aeff3j(63)+Atop3j(63)+Abottom3j(63)
      ampQG(64)=Heff3j(64)+Htop3j(64)+Hbottom3j(64)
     &         +Aeff3j(64)+Atop3j(64)+Abottom3j(64)
      ampQG(65)=Heff3j(65)+Htop3j(65)+Hbottom3j(65)
     &         +Aeff3j(65)+Atop3j(65)+Abottom3j(65)
      ampQG(66)=Heff3j(66)+Htop3j(66)+Hbottom3j(66)
     &         +Aeff3j(66)+Atop3j(66)+Abottom3j(66)
      ampQG(67)=Heff3j(67)+Htop3j(67)+Hbottom3j(67)
     &         +Aeff3j(67)+Atop3j(67)+Abottom3j(67)
      ampQG(68)=Heff3j(68)+Htop3j(68)+Hbottom3j(68)
     &         +Aeff3j(68)+Atop3j(68)+Abottom3j(68)
      ampQG(69)=Heff3j(69)+Htop3j(69)+Hbottom3j(69)
     &         +Aeff3j(69)+Atop3j(69)+Abottom3j(69)
                                                                
      ampQG(70)=Heff3j(70)+Htop3j(70)+Hbottom3j(70)
     &         +Aeff3j(70)+Atop3j(70)+Abottom3j(70)
      ampQG(71)=Heff3j(71)+Htop3j(71)+Hbottom3j(71)
     &         +Aeff3j(71)+Atop3j(71)+Abottom3j(71)
      ampQG(72)=Heff3j(72)+Htop3j(72)+Hbottom3j(72)
     &         +Aeff3j(72)+Atop3j(72)+Abottom3j(72)
      ampQG(73)=Heff3j(73)+Htop3j(73)+Hbottom3j(73)
     &         +Aeff3j(73)+Atop3j(73)+Abottom3j(73)
      ampQG(74)=Heff3j(74)+Htop3j(74)+Hbottom3j(74)
     &         +Aeff3j(74)+Atop3j(74)+Abottom3j(74)
      ampQG(75)=Heff3j(75)+Htop3j(75)+Hbottom3j(75)
     &         +Aeff3j(75)+Atop3j(75)+Abottom3j(75)
      ampQG(76)=Heff3j(76)+Htop3j(76)+Hbottom3j(76)
     &         +Aeff3j(76)+Atop3j(76)+Abottom3j(76)
      ampQG(77)=Heff3j(77)+Htop3j(77)+Hbottom3j(77)
     &         +Aeff3j(77)+Atop3j(77)+Abottom3j(77)
      ampQG(78)=Heff3j(78)+Htop3j(78)+Hbottom3j(78)
     &         +Aeff3j(78)+Atop3j(78)+Abottom3j(78)
      ampQG(79)=Heff3j(79)+Htop3j(79)+Hbottom3j(79)
     &         +Aeff3j(79)+Atop3j(79)+Abottom3j(79)
                                                                
      ampQG(80)=Heff3j(80)+Htop3j(80)+Hbottom3j(80)
     &         +Aeff3j(80)+Atop3j(80)+Abottom3j(80)
      ampQG(81)=Heff3j(81)+Htop3j(81)+Hbottom3j(81)
     &         +Aeff3j(81)+Atop3j(81)+Abottom3j(81)
      ampQG(82)=Heff3j(82)+Htop3j(82)+Hbottom3j(82)
     &         +Aeff3j(82)+Atop3j(82)+Abottom3j(82)
      ampQG(83)=Heff3j(83)+Htop3j(83)+Hbottom3j(83)
     &         +Aeff3j(83)+Atop3j(83)+Abottom3j(83)
      ampQG(84)=Heff3j(84)+Htop3j(84)+Hbottom3j(84)
     &         +Aeff3j(84)+Atop3j(84)+Abottom3j(84)

c... pentagon loops
      ampQG(85)=(Htop3j(85)+Hbottom3j(85)
     &          +Atop3j(85)+Abottom3j(85))
      ampQG(86)=(Htop3j(86)+Hbottom3j(86)
     &          +Atop3j(86)+Abottom3j(86))
      ampQG(87)=(Htop3j(87)+Hbottom3j(87)
     &          +Atop3j(87)+Abottom3j(87))

cc start check
c      print*,'ampQG(85)*0.5d0',ampQG(85)*0.5d0
c      print*,'ampQG(86)*0.5d0',ampQG(86)*0.5d0
c      print*,'ampQG(87)*0.5d0',ampQG(87)*0.5d0
c      print*,'Heff3j(36)',Heff3j(36)
c      print*,'Heff3j(53)',Heff3j(53)
c      print*,'Heff3j(65)',Heff3j(65)
c      print*,'c1',-ampQG(53)-ampQG(65)
c      print*,'c2',+ampQG(36)+ampQG(65)
c      print*,'c3',-ampQG(36)+ampQG(53)
c      print*,'c4', ampQG(36)+ampQG(65)
c      print*,'c5',-ampQG(36)+ampQG(53)
c      print*,'c6',-ampQG(53)-ampQG(65)
cc end check

      color(1)=+ampQG(1)+ampQG(2)+ampQG(3)-ampQG(12)+ampQG(13)
     &         +ampQG(14)-ampQG(15)-ampQG(16)+ampQG(20)+ampQG(21)
     &         -ampQG(22)-ampQG(23)-ampQG(24)-ampQG(42)+ampQG(43)
     &         +ampQG(44)+ampQG(45)+ampQG(46)-ampQG(47)+ampQG(48)
     &         -ampQG(49)+ampQG(50)-ampQG(51)-ampQG(52)-ampQG(53)
     &         +ampQG(54)+ampQG(55)+ampQG(56)-ampQG(57)+ampQG(58)
     &         -ampQG(59)-ampQG(60)-ampQG(61)+ampQG(62)+ampQG(63)
     &         -ampQG(64)-ampQG(65)+ampQG(85)*0.5d0

      color(2)=-ampQG(4)+ampQG(5)-ampQG(6)+ampQG(7)+ampQG(8)
     &         +ampQG(9)+ampQG(10)+ampQG(11)-ampQG(25)-ampQG(26)
     &         +ampQG(27)+ampQG(28)-ampQG(29)+ampQG(30)-ampQG(31)
     &         +ampQG(32)-ampQG(33)-ampQG(34)+ampQG(35)+ampQG(36)
     &         +ampQG(37)+ampQG(38)-ampQG(39)-ampQG(40)-ampQG(41)
     &         -ampQG(54)-ampQG(55)-ampQG(56)+ampQG(57)-ampQG(58)
     &         +ampQG(59)+ampQG(60)+ampQG(61)-ampQG(62)-ampQG(63)
     &         +ampQG(64)+ampQG(65)+ampQG(86)*0.5d0

      color(3)=+ampQG(4)-ampQG(5)+ampQG(6)-ampQG(7)-ampQG(8)
     &         +ampQG(12)-ampQG(13)-ampQG(14)+ampQG(15)+ampQG(16)
     &         +ampQG(17)+ampQG(18)+ampQG(19)+ampQG(25)+ampQG(26)
     &         -ampQG(27)-ampQG(28)+ampQG(29)-ampQG(30)+ampQG(31)
     &         -ampQG(32)+ampQG(33)+ampQG(34)-ampQG(35)-ampQG(36)
     &         +ampQG(42)-ampQG(43)-ampQG(44)-ampQG(45)-ampQG(46)
     &         +ampQG(47)-ampQG(48)+ampQG(49)-ampQG(50)+ampQG(51)
     &         +ampQG(52)+ampQG(53)+ampQG(87)*0.5d0

      color(4)=-ampQG(20)-ampQG(21)+ampQG(22)+ampQG(23)+ampQG(24)
     &         -ampQG(25)-ampQG(26)+ampQG(27)+ampQG(28)-ampQG(29)
     &         +ampQG(30)-ampQG(31)+ampQG(32)-ampQG(33)-ampQG(34)
     &         +ampQG(35)+ampQG(36)-ampQG(54)-ampQG(55)-ampQG(56)
     &         +ampQG(57)-ampQG(58)+ampQG(59)+ampQG(60)+ampQG(61)
     &         -ampQG(62)-ampQG(63)+ampQG(64)+ampQG(65)+ampQG(66)
     &         +ampQG(67)+ampQG(68)+ampQG(69)-ampQG(70)+ampQG(71)
     &         -ampQG(72)-ampQG(73)+ampQG(86)*0.5d0

      color(5)=+ampQG(25)+ampQG(26)-ampQG(27)-ampQG(28)+ampQG(29)
     &         -ampQG(30)+ampQG(31)-ampQG(32)+ampQG(33)+ampQG(34)
     &         -ampQG(35)-ampQG(36)+ampQG(42)-ampQG(43)-ampQG(44)
     &         -ampQG(45)-ampQG(46)+ampQG(47)-ampQG(48)+ampQG(49)
     &         -ampQG(50)+ampQG(51)+ampQG(52)+ampQG(53)-ampQG(69)
     &         +ampQG(70)-ampQG(71)+ampQG(72)+ampQG(73)-ampQG(77)
     &         +ampQG(78)+ampQG(79)-ampQG(80)-ampQG(81)+ampQG(82)
     &         +ampQG(83)+ampQG(84)+ampQG(87)*0.5d0

      color(6)=-ampQG(37)-ampQG(38)+ampQG(39)+ampQG(40)+ampQG(41)
     &         -ampQG(42)+ampQG(43)+ampQG(44)+ampQG(45)+ampQG(46)
     &         -ampQG(47)+ampQG(48)-ampQG(49)+ampQG(50)-ampQG(51)
     &         -ampQG(52)-ampQG(53)+ampQG(54)+ampQG(55)+ampQG(56)
     &         -ampQG(57)+ampQG(58)-ampQG(59)-ampQG(60)-ampQG(61)
     &         +ampQG(62)+ampQG(63)-ampQG(64)-ampQG(65)+ampQG(74)
     &         +ampQG(75)+ampQG(76)+ampQG(77)-ampQG(78)-ampQG(79)
     &         +ampQG(80)+ampQG(81)+ampQG(85)*0.5d0

      do i=1,6
         colorT=DCMPLX(0d0,0d0)
         do j=1,6
            colorT=colorT+RGB(i,j)*color(j)
         enddo
         RH=RH+(dreal(colorT)*dreal(color(i))
     &        + dimag(colorT)*dimag(color(i)))/9d0   
      enddo
      RHpenSq=RHPenSq+5d0/18d0*(dreal(ampQG(85)+ampQG(86)+ampQG(87))
     &                         *dreal(ampQG(85)+ampQG(86)+ampQG(87))
     &                         +dimag(ampQG(85)+ampQG(86)+ampQG(87))
     &                         *dimag(ampQG(85)+ampQG(86)+ampQG(87)))
      do i=1,6
         RHpenI=RHpenI - 5d0/9d0
     &         *(dreal(ampQG(85)+ampQG(86)+ampQG(87))*dreal(color(i))
     &          +dimag(ampQG(85)+ampQG(86)+ampQG(87))*dimag(color(i)))
      enddo

c$$$      do i=1,87
c$$$      print*,'ampQG(',i,')',ampQG(i)/(dsqrt(RH+RHpenSq+RHpenI)) !*0.5d0*dsqrt(fac)
c$$$         print*,'ampQG(',i,')',ampQG(i)!/(dsqrt(RH+RHpenSq+RHpenI))!*0.5d0*dsqrt(fac)
c$$$   enddo
c$$$   print*,''
c$$$   do i=1,6
c$$$         print*,'color(',i,')',color(i)!*dsqrt(fac)
c$$$   enddo
c$$$   print*,''

c$$$      print*,'RH               ',RH
c$$$      print*,'RHpenSq          ',RHpenSq
c$$$      print*,'RHpenI           ',RHpenI
c$$$      print*,'RH+RHpenSq+RHpenI',RH+RHpenSq+RHpenI
c$$$      print*,'------------------------------------'

      enddo
      enddo
      enddo

      RHsum=(RH+RHpenSq+RHpenI)*fac
cfc      print*, "RHSUM",RHSUM

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
c 
c         Check whether there is instabilities
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c!                                                       !! 
c!         IF SO ::::-->   RSUM -> 0                     !!
c!                                                       !!
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
c
c         Count instabilites for pen,box,vertices
c         both in quad and dble precision
c           
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c   Check Vertices
c     
c        if(boxeva) then
c   Dble precision
           If(verdbleSin) then
c          print*,"verdbleSi",verdbleSin
c          stop
c bad dble
        verinst(gf_iteration,1,1)=verinst(gf_iteration,1,1)+1
c        print*, " verinst(gf_iteration,1,1)", verinst(gf_iteration,1,1)
       else
c good dble 
        verinst(gf_iteration,2,1)=verinst(gf_iteration,2,1)+1
c        print*, " verinst(gf_iteration,2,1)", verinst(gf_iteration,2,1)
        endif

       If(verQUADSin) then
c bad QUAD
        verinst(gf_iteration,1,2)=verinst(gf_iteration,1,2)+1
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        RHsum=0d0
c
cccccccccccccccccccccccccccccccccccccccccccccccccccc
       else
c good QUAD 
        verinst(gf_iteration,2,2)=verinst(gf_iteration,2,2)+1
        endif
c        endif
c 
c   Check Boxes
c     
        if(boxeva) then
c   Dble precision
       If(boxdbleSin) then
c bad dble
        boxinst(gf_iteration,1,1)=boxinst(gf_iteration,1,1)+1
       else
c good dble 
        boxinst(gf_iteration,2,1)=boxinst(gf_iteration,2,1)+1
        endif

       If(boxQUADSin) then
c bad QUAD
        boxinst(gf_iteration,1,2)=boxinst(gf_iteration,1,2)+1
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        RHsum=0d0
c
cccccccccccccccccccccccccccccccccccccccccccccccccccc
       else
c good QUAD 
        boxinst(gf_iteration,2,2)=boxinst(gf_iteration,2,2)+1
        endif
        endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c   Check Penagons
c     
        if(peneva) then
c   Dble precision
       If(pendbleSin) then
c bad dble
        peninst(gf_iteration,1,1)=peninst(gf_iteration,1,1)+1
       else
c good dble 
        peninst(gf_iteration,2,1)=peninst(gf_iteration,2,1)+1
        endif

       If(penQUADSin) then
c bad QUAD
        peninst(gf_iteration,1,2)=peninst(gf_iteration,1,2)+1
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        RHsum=0d0
c
cccccccccccccccccccccccccccccccccccccccccccccccccccc
       else
c good QUAD 
        peninst(gf_iteration,2,2)=peninst(gf_iteration,2,2)+1
        endif
        endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        if(MGdebug) then
         GGFLOampSq2=RHsum
         do i=1,6
            GGFLOcolor(i)=color(i)*dsqrt(fac)
c            print*,'color(',i,')',color(i)*dsqrt(fac)
c            STOP
         enddo
      endif

      if(lgauge2) then
         auxGauge=0.0d0
         auxGauge2=0.0d0

         print*, "boxes"
         print*
         print*, 7, ampQG(7)         
         print*, 15, ampQG(15)   
         i=23
         print*, i, ampQG(i)   
         i=31
         print*, i, ampQG(i)   
         i=32
         print*, i, ampQG(i)   
         i=40
         print*, i, ampQG(i)   
         i=48
         print*, i, ampQG(i)   
         i=49
         print*, i, ampQG(i)   
         i=60
         print*, i, ampQG(i)   
         i=61
         print*, i, ampQG(i)   
         i=72
         print*, i, ampQG(i)   
         i=80
         print*, i, ampQG(i)   

        print*, "pentagons"
         print*
         i=36
         print*, i, ampQG(i)   

         i=53
         print*, i, ampQG(i)   

         i=65
         print*, i, ampQG(i)   

          print*
         i=85
         print*, i, -2d0*(ampQG(65)+ampQG(53)   )
         i=86
         print*, i, 2d0*(ampQG(36)+ampQG(65)   )
         i=87
         print*, i, 2d0*(ampQG(53)-ampQG(36)   )

         print*
       print*, "vertices"
         print*
         do i=1,6
                print*, i, ampQG(i)
           auxGauge=auxGauge+abs(ampQG(i)*dsqrt(fac))
         enddo
         do i=8,14
                print*, i, ampQG(i)
           auxGauge=auxGauge+abs(ampQG(i)*dsqrt(fac))
         enddo
         do i=16,22
                print*, i, ampQG(i)
           auxGauge=auxGauge+abs(ampQG(i)*dsqrt(fac))
        enddo
        do i=24,30
                print*, i, ampQG(i)
           auxGauge=auxGauge+abs(ampQG(i)*dsqrt(fac))
         enddo
         do i=33,35
                print*, i, ampQG(i)
           auxGauge=auxGauge+abs(ampQG(i)*dsqrt(fac))
         enddo

        do i=37,39
                print*, i, ampQG(i)
           auxGauge=auxGauge+abs(ampQG(i)*dsqrt(fac))
         enddo

         do i=41,47
                print*, i, ampQG(i)
           auxGauge=auxGauge+abs(ampQG(i)*dsqrt(fac))
         enddo
         do i=50,52
                print*, i, ampQG(i)
           auxGauge=auxGauge+abs(ampQG(i)*dsqrt(fac))
         enddo
         do i=54,59
                print*, i, ampQG(i)
           auxGauge=auxGauge+abs(ampQG(i)*dsqrt(fac))
         enddo

         do i=62,64
                print*, i, ampQG(i)
           auxGauge=auxGauge+abs(ampQG(i)*dsqrt(fac))
         enddo
         do i=66,71
                print*, i, ampQG(i)
           auxGauge=auxGauge+abs(ampQG(i)*dsqrt(fac))
         enddo
         do i=73,79
                print*, i, ampQG(i)
           auxGauge=auxGauge+abs(ampQG(i)*dsqrt(fac))
         enddo
         do i=81,100
                print*, i, ampQG(i)
           auxGauge=auxGauge+abs(ampQG(i)*dsqrt(fac))
         enddo

         do i=1,6
            auxGauge2=auxGauge2+abs(color(i)*dsqrt(fac))
         enddo
         print*,'auxGauge              ',auxGauge
         print*,'auxGauge2             ',auxGauge2
         print*,''
         print*,'i1pol,i2pol,i3pol',i1pol,i2pol,i3pol
         do i=1,6
            print*,'color(',i,')         ',color(i)
         enddo
         print*,''
         do i=1,6
            print*,'color(',i,')/auxGauge ',color(i)*dsqrt(fac)/auxGauge
         enddo
         do i=1,6
            print*,'color(',i,')/auxGauge2',color(i)*dsqrt(fac)/auxGauge2
         enddo

         print*,'RHsum / (auxGauge )^2 ',RHsum/(auxGauge**2)
         print*,'RHsum / (auxGauge2)^2 ',RHsum/(auxGauge2**2)
         do i=1,6
            print*,'RHsum / (color(i))**2 ',RHsum/(abs(color(i)
     &           *dsqrt(fac)))**2,i
         enddo
cc start check
         print*,''
         print*,'fac                     ',fac
         print*,'RH                      ',RH
         print*,'RHpenSq                 ',RHpenSq
         print*,'RHpenI                  ',RHpenI
         print*,'RH+RHpenSq+RHpenI       ',RH+RHpenSq+RHpenI
         print*,'RHsum                   ',RHsum
         print*,''
         print*,''
         print*,'RH/  (auxGauge )^2      ',RH*fac/(auxGauge2**2)
         print*,'RHpenSq/  (auxGauge )^2 ',RHpenSq*fac/(auxGauge2**2)
         print*,'RHpenI/  (auxGauge )^2  ',RHpenI*fac/(auxGauge2**2)
         print*,'RH+RHpenSq+RHpenI/(auxGauge )^2',(RH+RHpenSq+RHpenI)*fac/
     &   (auxGauge2**2)
         print*, "(auxGauge2)^2",(auxGauge2**2)
         print*,'RHsum /  (auxGauge2 )^2   ',RHsum/(auxGauge2**2)
         print*,''
         STOP

c         print*,' ' 
c         print*, "polarization GGFLO"
c         print*, ("W(1,2)",epsChelas(mu, 1,1), mu=1,18)
c         print*, ("ggf(1,2)",epsgC(mu,1,1), mu=0,3)
c         print*,'PBAR(mu,2)',(KBAR(mu,1),mu=0,3)
c         print*,' '


cc end check
      endif

cc start check
c$$$         print*,''
c$$$         print*,'RH                      ',RH
c$$$         print*,'RHpenSq                 ',RHpenSq
c$$$         print*,'RHpenI                  ',RHpenI
c$$$         print*,'RH+RHpenSq+RHpenI       ',RH+RHpenSq+RHpenI
c$$$         print*,'RHsum                   ',RHsum
c$$$         print*,''
c$$$         STOP
cc end check



      RETURN

c 100  RHsum=(RH+RHpenSq+RHpenI)*fac
c      print*,'RHsum',RHsum

      END

