***********************************************************************
      SUBROUTINE AQG_bottom(PBAR,FSIGN,FPOL,KBAR,GSIGN,PHBAR,HSIGN,
     &                   i1polext,i2polext,i3polext,Abottom3j)
*----------------------------------------------------------------------
c... author: Michael Kubocz
c... date  : 06.02.2013
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
      integer i1pol,i2pol,i3pol,i1polext,i2polext,i3polext,i1,i2,i3,i4,
     &     i1polmax,i2polmax,i3polmax

c... Effective polarization vectors (with Higgs-vertex)
      complex*16 epsAbottom(0:3,2,3),j21Abottom(0:3)

c... Spinors and quark currents with and without gluon emision
      complex*16 psi(2,-1:1,2),braketg(2,-1:1,12,2),
     &         j21aux(0:3,-1:1),j21auxi(0:3,-1:1,2),j21auxf(0:3,-1:1,2),
     &         j21(0:3),j21i(0:3,2,3),j21f(0:3,2,3),j21if(0:3,-1:1,2,2,6)
      integer isig,sigma

c... Quark currents with gluon and Higgs emission
      complex*16 j21Abottomf(0:3,2,3),j21Abottomi(0:3,2,3),
     &           j21Abottomif(0:3,-1:1,2,2,6),j21iAbottomf(0:3,-1:1,2,2,6)

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

      complex*16 GGGPAbottom12(0:4,2,2),GGGPAbottom13(0:4,2,2),
     &           GGGPAbottom23(0:4,2,2),GGGPAbottom32(0:4,2,2),
     &           GGG1PAbottom2(0:4,2,2),GGG1PAbottom3(0:4,2,2)

c... Effective vertices and loop functions     
      complex*16 A2G_Loop,A3G_Loop
      complex*16 GGG_vertex,H4G_effective2
      complex*16 qpentA_bottom1234(4),qpentA_bottom1342(4),qpentA_bottom1423(4)

c... Color matrix
      real*8 RGB(6,6)

c... Further local variables
      complex*16 color(6),colorT,colorPen(6),ampQG(87)
      complex*16 Abottom3j(87)
      integer ienter,mu,i,j,l
      complex*16 auxPen,auxRest
      real*8 pi,fac,vev,DUMMY(0:4),Kpen(0:4,4),RH,RHpenSq,RHpenI,RHpenI2,
     &     RHsum,RHsumFi
      complex*16 EPSpen(0:3,15)
      parameter (pi=3.14159 26535 89793d0)
      save fac
      data ienter /0/
      save ienter
c... Debbuging variables
      real*8 polar(0:3,4), phiggs(0:3),ratioVal
      complex*16 polarC(0:3,4)
      complex*16 A3G_LoopGC,A4G_Loop,pentP,epsChelas(1:6,-1:1,3)
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
     &         A2G_Loop,A3G_Loop,!A4G_LoopDD2,
     &         GGG_vertex,
     &         A3G_LoopGC,A4G_Loop,A4G_LoopDD
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

cfc      verdbleSin=.false.
cfc      verquadSin=.false.
cfc      boxdbleSin=.false.
cfc      boxquadSin=.false.
cfc      boxSingular(1)=.false.
cfc      boxSingular(2)=.false.
cfc      verSingular(1)=.false.
cfc      verSingular(2)=.false.
cfc      
cfc      SingularC(1)=.false.
cfc      SingularC(2)=.false.

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
         Abottom3j(i)=DCMPLX(0d0,0d0)
      enddo
      RH      = 0.0d0
      RHPenSq = 0.0d0
      RHpenI  = 0.0d0
      RHsum   = 0.0d0

*********************************************************************
**************** CP-ODD HIGGS WITH BOTTOM-LOOP *********************
*********************************************************************
      if (AloopBottom) then
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
         CALL EpsA_Loop(p21,j21,PHp,PHp21p,j21Abottom,m_bottom)
         PHp21(4)=dotrr(PHp21p,PHp21p) ! both gluons are virtual !!!
         do i1pol=1,2           ! gluon polarisations
            CALL EpsA_Loop(Kp(0,1),epsgC(0,i1pol,1),PHp,PHK1p,
     &                          epsAbottom(0,i1pol,1),m_bottom)
            PHK1(4)=dotrr(PHK1p,PHK1p)
            CALL EpsA_Loop(Kp(0,2),epsgC(0,i1pol,2),PHp,PHK2p,
     &                          epsAbottom(0,i1pol,2),m_bottom)
            PHK2(4)=dotrr(PHK2p,PHK2p)
            CALL EpsA_Loop(Kp(0,3),epsgC(0,i1pol,3),PHp,PHK3p,
     &                       epsAbottom(0,i1pol,3),m_bottom)
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
     &           epsAbottom(0,i1pol,1),braketg(1,isig,7,i1pol),DUMMY)  ! |H,K1,1>_i1pol,SIGMA1
            CALL bra2c(psi(1,isig,2),.false.,p(0,2),isig,PHK1,
     &           epsAbottom(0,i1pol,1),braketg(1,isig,8,i1pol),DUMMY)  ! <2,H,K1|_i1pol,SIGMA2
            CALL ket2c(psi(1,isig,1),.false.,p(0,1),isig,PHK2,
     &           epsAbottom(0,i1pol,2),braketg(1,isig,9,i1pol),DUMMY)  ! |H,K2,1>_2i1pol,SIGMA1
            CALL bra2c(psi(1,isig,2),.false.,p(0,2),isig,PHK2,
     &           epsAbottom(0,i1pol,2),braketg(1,isig,10,i1pol),DUMMY) ! <2,H,K2|_i1pol,SIGMA2
            CALL ket2c(psi(1,isig,1),.false.,p(0,1),isig,PHK3,
     &           epsAbottom(0,i1pol,3),braketg(1,isig,11,i1pol),DUMMY) ! |H,K3,1>_2i1pol,SIGMA1
            CALL bra2c(psi(1,isig,2),.false.,p(0,2),isig,PHK3,
     &           epsAbottom(0,i1pol,3),braketg(1,isig,12,i1pol),DUMMY) ! <2,H,K3|_i1pol,SIGMA2
         enddo

            CALL curr(1,psi(1,-1,2),braketg(1,-1,7,i1pol),
     &                j21auxi(0,-1,i1pol))        !=  <2|vertex|H,K1,1>_i1pol,isig1
            CALL curr(1,braketg(1,-1,8,i1pol),psi(1,-1,1),
     &                j21auxf(0,-1,i1pol))        !=  <2,H,K1|vertex|1>_i1pol,isig1
            do mu=0,3 !with attached propagator !!!
               j21Abottomi(mu,i1pol,1)=j21auxi(mu,SIGMA,i1pol)/p21K1H(4)
               j21Abottomf(mu,i1pol,1)=j21auxf(mu,SIGMA,i1pol)/p21K1H(4)
            enddo
            CALL curr(1,psi(1,-1,2),braketg(1,-1,9,i1pol),
     &                j21auxi(0,-1,i1pol))        !=  <2|vertex|H,K2,1>_i1pol,isig1
            CALL curr(1,braketg(1,-1,10,i1pol),psi(1,-1,1),
     &                j21auxf(0,-1,i1pol))        !=  <2,H,K2|vertex|1>_i1pol,isig1
            do mu=0,3 !with attached propagator !!!
               j21Abottomi(mu,i1pol,2)=j21auxi(mu,SIGMA,i1pol)/p21K2H(4)
               j21Abottomf(mu,i1pol,2)=j21auxf(mu,SIGMA,i1pol)/p21K2H(4)
            enddo
            CALL curr(1,psi(1,-1,2),braketg(1,-1,11,i1pol),
     &                j21auxi(0,-1,i1pol))        !=  <2|vertex|H,K3,1>_i1pol,isig1
            CALL curr(1,braketg(1,-1,12,i1pol),psi(1,-1,1),
     &                j21auxf(0,-1,i1pol))        !=  <2,H,K3|vertex|1>_i1pol,isig1
            do mu=0,3 !with attached propagator !!!
               j21Abottomi(mu,i1pol,3)=j21auxi(mu,SIGMA,i1pol)/p21K3H(4)
               j21Abottomf(mu,i1pol,3)=j21auxf(mu,SIGMA,i1pol)/p21K3H(4)
            enddo
      enddo

      do i1pol=1,2              ! gluon polarizations
         do i2pol=1,2           ! gluon polarizations
c  Higgs production in initial state
            CALL curr(1,braketg(1,-1,4,i1pol),braketg(1,-1,7,i2pol),
     &           j21Abottomif(0,-1,i1pol,i2pol,1))        !=  <2,K2|vertex|H,K1,1>

            CALL curr(1,braketg(1,-1,2,i1pol),braketg(1,-1,9,i2pol),
     &           j21Abottomif(0,-1,i1pol,i2pol,2))        !=  <2,K1|vertex|H,K2,1>

            CALL curr(1,braketg(1,-1,2,i1pol),braketg(1,-1,11,i2pol),
     &           j21Abottomif(0,-1,i1pol,i2pol,3))        !=  <2,K1|vertex|H,K3,1>

            CALL curr(1,braketg(1,-1,6,i1pol),braketg(1,-1,7,i2pol),
     &           j21Abottomif(0,-1,i1pol,i2pol,4))        !=  <2,K3|vertex|H,K1,1>

            CALL curr(1,braketg(1,-1,6,i1pol),braketg(1,-1,9,i2pol),
     &           j21Abottomif(0,-1,i1pol,i2pol,5))        !=  <2,K3|vertex|H,K2,1>

            CALL curr(1,braketg(1,-1,4,i1pol),braketg(1,-1,11,i2pol),
     &           j21Abottomif(0,-1,i1pol,i2pol,6))        !=  <2,K2|vertex|H,K3,1>

c  Higgs production in final state
            CALL curr(1,braketg(1,-1,10,i1pol),braketg(1,-1,1,i2pol),
     &           j21iAbottomf(0,-1,i1pol,i2pol,1))        !=  <2,H,K2|vertex|K1,1>

            CALL curr(1,braketg(1,-1,8,i1pol),braketg(1,-1,3,i2pol),
     &           j21iAbottomf(0,-1,i1pol,i2pol,2))        !=  <2,H,K1|vertex|K2,1>

            CALL curr(1,braketg(1,-1,8,i1pol),braketg(1,-1,5,i2pol),
     &           j21iAbottomf(0,-1,i1pol,i2pol,3))        !=  <2,H,K1|vertex|K3,1>

            CALL curr(1,braketg(1,-1,12,i1pol),braketg(1,-1,1,i2pol),
     &           j21iAbottomf(0,-1,i1pol,i2pol,4))        !=  <2,H,K3|vertex|K1,1>

            CALL curr(1,braketg(1,-1,10,i1pol),braketg(1,-1,5,i2pol),
     &           j21iAbottomf(0,-1,i1pol,i2pol,5))        !=  <2,H,K2|vertex|K3,1>

            CALL curr(1,braketg(1,-1,12,i1pol),braketg(1,-1,3,i2pol),
     &           j21iAbottomf(0,-1,i1pol,i2pol,6))        !=  <2,H,K3|vertex|K2,1>
         enddo
      enddo
      do i=1,2
         do j=1,2
            CALL  GGG_current(PHK1p,epsAbottom(0,i,1),Kp(0,2),epsgC(0,j,2),
     &                        DUMMYp,GGGPAbottom12(0,i,j)) !54
            CALL  GGG_current(PHK1p,epsAbottom(0,i,1),Kp(0,3),epsgC(0,j,3),
     &                        DUMMYp,GGGPAbottom13(0,i,j)) !25
            CALL  GGG_current(PHK2p,epsAbottom(0,i,2),Kp(0,3),epsgC(0,j,3),
     &                        DUMMYp,GGGPAbottom23(0,i,j)) !44
            CALL  GGG_current(PHK3p,epsAbottom(0,i,3),Kp(0,2),epsgC(0,j,2),
     &                        DUMMYp,GGGPAbottom32(0,i,j)) !46
            CALL  GGG_current(Kp(0,1),epsgC(0,i,1),PHK2p,epsAbottom(0,j,2),
     &                        DUMMYp,GGG1PAbottom2(0,i,j)) !56
            CALL  GGG_current(Kp(0,1),epsgC(0,i,1),PHK3p,epsAbottom(0,j,3),
     &                        DUMMYp,GGG1PAbottom3(0,i,j)) !29
         enddo
      enddo
      endif

c... MGdebug
c$$$      if(MGdebug) then
c$$$         i1polmax=1
c$$$         i2polmax=1
c$$$         i3polmax=1
c$$$      else
c$$$         i1polmax=2
c$$$         i2polmax=2
c$$$         i3polmax=2
c$$$      endif

*********************************************************************
**************** CP-ODD HIGGS WITH BOTTOM-LOOP *********************
*********************************************************************
      if (AloopBottom) then
         do mu=0,3
            EPSpen(mu,1)=epsgC(mu,i1polext,1)
            EPSpen(mu,2)=epsgC(mu,i2polext,2)
            EPSpen(mu,3)=epsgC(mu,i3polext,3)
            EPSpen(mu,4)=j21(mu)
ccc NEW for pentagons 
            EPSgpen(mu,1)=epsgC(mu,i1polext,1)
            EPSgpen(mu,2)=epsgC(mu,i2polext,2)
            EPSgpen(mu,3)=epsgC(mu,i3polext,3)
            EPSgpen(mu,4)=j21(mu)
             PHc(mu)=ph(mu)
         enddo
         
         if(peneva) then
cfc        pendbleSin=.false.
cfc        penquadSin=.false.
         eva=.true.
      PenSingular(1)=.false.
      PenSingular(2)=.false.
cfc         CALL auxiliary_calculation(EPSpen(0,1),Kpen(0,1))
cfc         CALL A4G_LoopDD(4,3,2,1,m_bottom,qpentA_bottom1234(1))
cfc         CALL A4G_LoopDD(2,3,4,1,m_bottom,qpentA_bottom1234(2)) 
cfc         CALL A4G_LoopDD(2,1,4,3,m_bottom,qpentA_bottom1234(3)) 
cfc         CALL A4G_LoopDD(3,2,1,4,m_bottom,qpentA_bottom1234(4))
cfc        Abottom3j(85)=cplModAbottom*(qpentA_bottom1234(1)+qpentA_bottom1234(2)
cfc     &                         +qpentA_bottom1234(3)+qpentA_bottom1234(4))
cFC c
c         print*, " Abottom (85) before",  Abottom3j(85)
c        Abottom3j(85)=Abottom3j(85)/(FDeffPenLoop(4,3,2,1)+FDeffPenLoop(2,3,4,1)+
c     &  FDeffPenLoop(2,1,4,3) +   FDeffPenLoop(3,2,1,4) )

      selectAbottom=.true.
      
         Abottom3j(85)=(FDeffPenLoop(4,3,2,1)+FDeffPenLoop(2,3,4,1)+
     &  FDeffPenLoop(2,1,4,3) +   FDeffPenLoop(3,2,1,4) )
      IF(PenSingular(1)) then
         PendbleSin=.true.
      endif
      IF(PenSingular(2)) then
         PenquadSin=.true.
      endif
c         print*, " Abottom (85) after",  abs(Abottom3j(85)-1d0)

cfc         CALL A4G_LoopDD(2,1,3,4,m_bottom,qpentA_bottom1342(1))
cfc         CALL A4G_LoopDD(3,4,2,1,m_bottom,qpentA_bottom1342(2))
cfc         CALL A4G_LoopDD(4,2,1,3,m_bottom,qpentA_bottom1342(3))
cfc         CALL A4G_LoopDD(2,4,3,1,m_bottom,qpentA_bottom1342(4))
cfc         
cfc         Abottom3j(86)=cplModAbottom*(qpentA_bottom1342(1)+qpentA_bottom1342(2)
cfc     &                         +qpentA_bottom1342(3)+qpentA_bottom1342(4))
cfc
c        print*, " Abottom (86) before",  Abottom3j(86)
c        Abottom3j(86)=Abottom3j(86)/(FDeffPenLoop(2,1,3,4)+FDeffPenLoop(3,4,2,1)+
c     &               FDeffPenLoop(4,2,1,3)+FDeffPenLoop(2,4,3,1))
      PenSingular(1)=.false.
      PenSingular(2)=.false.
         Abottom3j(86)=(FDeffPenLoop(2,1,3,4)+FDeffPenLoop(3,4,2,1)+
     &               FDeffPenLoop(4,2,1,3)+FDeffPenLoop(2,4,3,1))
      IF(PenSingular(1)) then
         PendbleSin=.true.
      endif
      IF(PenSingular(2)) then
         PenquadSin=.true.
      endif
c         print*, " Abottom (86) after", abs( Abottom3j(86)-1d0)
cfc
cfc         CALL A4G_LoopDD(2,3,1,4,m_bottom,qpentA_bottom1423(1))
cfc         CALL A4G_LoopDD(3,1,4,2,m_bottom,qpentA_bottom1423(2))
cfc         CALL A4G_LoopDD(4,2,3,1,m_bottom,qpentA_bottom1423(3))
cfc         CALL A4G_LoopDD(3,2,4,1,m_bottom,qpentA_bottom1423(4))
cfc         
cfc         Abottom3j(87)=cplModAbottom*(qpentA_bottom1423(1)+qpentA_bottom1423(2)
cfc     &                         +qpentA_bottom1423(3)+qpentA_bottom1423(4))
c         print*, " Abottom (87) before",  Abottom3j(87)

c        Abottom3j(87)=Abottom3j(87)/(FDeffPenLoop(2,3,1,4)+FDeffPenLoop(3,1,4,2)+
c     &               FDeffPenLoop(4,2,3,1)+FDeffPenLoop(3,2,4,1))
      PenSingular(1)=.false.
      PenSingular(2)=.false.
         Abottom3j(87)=(FDeffPenLoop(2,3,1,4)+FDeffPenLoop(3,1,4,2)+
     &     FDeffPenLoop(4,2,3,1)+FDeffPenLoop(3,2,4,1))

      selectAbottom=.false.
         
      IF(PenSingular(1)) then
         PendbleSin=.true.
      endif
      IF(PenSingular(2)) then
         PenquadSin=.true.
      endif
c         print*, " Abottom (87) after",abs(  Abottom3j(87)-1d0)
c         pause

         eva=.false.
         endif

         Abottom3j(1)=cplModAbottom                                     !1
     &               *dotcc(j21iAbottomf(0,SIGMA,i1polext,i3polext,3),epsgC(0,i2polext,2))
*----------------------------------------------------------------------
         Abottom3j(2)=cplModAbottom                                     !2
     &               *dotcc(j21if(0,SIGMA,i1polext,i3polext,4),epsAbottom(0,i2polext,2)) 
*----------------------------------------------------------------------  
         Abottom3j(3)=cplModAbottom                                     !3
     &               *dotcc(j21Abottomif(0,SIGMA,i1polext,i3polext,3),epsgC(0,i2polext,2))
*----------------------------------------------------------------------
         Abottom3j(4)=cplModAbottom*GGG_vertex(Kp(0,3),p21K2p,PHK1p,       !4 
     &                epsgC(0,i3polext,3),j21i(0,i2polext,2),epsAbottom(0,i1polext,1))
*----------------------------------------------------------------------
         Abottom3j(5)=-cplModAbottom*GGG_vertex(Kp(0,3),p21K2Hp,Kp(0,1),   !5 minus
     &                epsgC(0,i3polext,3),j21Abottomi(0,i2polext,2),epsgC(0,i1polext,1))
*----------------------------------------------------------------------
         Abottom3j(6)=cplModAbottom*GGG_vertex(PHK3p,p21K2p,Kp(0,1),       !6 
     &                epsAbottom(0,i3polext,3),j21i(0,i2polext,2),epsgC(0,i1polext,1))
*----------------------------------------------------------------------
         Abottom3j(7)=-cplModAbottom                                    !7
     &      *A3G_Loop(K(0,3),p21K2,K(0,1),epsgC(0,i3polext,3),
     &                             j21i(0,i2polext,2),epsgC(0,i1polext,1),m_bottom)
*----------------------------------------------------------------------
         Abottom3j(8)=cplModAbottom                                     !8 
     &    *A2G_Loop(K13,p21K2,GGG13(0,i1polext,i3polext),j21i(0,i2polext,2),m_bottom)
*----------------------------------------------------------------------
         Abottom3j(9)=cplModAbottom                                     !9
     &     *dotcc(j21if(0,SIGMA,i3polext,i2polext,5),epsAbottom(0,i1polext,1)) 
*----------------------------------------------------------------------
         Abottom3j(10)=cplModAbottom                                    !10
     &     *dotcc(j21Abottomif(0,SIGMA,i3polext,i2polext,5),epsgC(0,i1polext,1))
*----------------------------------------------------------------------
         Abottom3j(11)=cplModAbottom                                    !11
     &               *dotcc(j21iAbottomf(0,SIGMA,i3polext,i2polext,6),epsgC(0,i1polext,1))
*----------------------------------------------------------------------
         Abottom3j(12)=-cplModAbottom*GGG_vertex(p21K1Hp,Kp(0,2),Kp(0,3),  !12 minus
     &                j21Abottomf(0,i1polext,1),epsgC(0,i2polext,2),epsgC(0,i3polext,3))
*----------------------------------------------------------------------
         Abottom3j(13)=-cplModAbottom*GGG_vertex(p21K1p,Kp(0,3),PHK2p,     !13 minus
     &                j21f(0,i1polext,1),epsgC(0,i3polext,3),epsAbottom(0,i2polext,2))
*----------------------------------------------------------------------
         Abottom3j(14)=cplModAbottom*GGG_vertex(p21K1p,Kp(0,2),PHK3p,      !14  
     &                j21f(0,i1polext,1),epsgC(0,i2polext,2),epsAbottom(0,i3polext,3))
*----------------------------------------------------------------------
         Abottom3j(15)=cplModAbottom                                    !15
     &              *A3G_Loop(p21K1,K(0,3),K(0,2),j21f(0,i1polext,1),
     &                             epsgC(0,i3polext,3),epsgC(0,i2polext,2),m_bottom)
*----------------------------------------------------------------------
         Abottom3j(16)=cplModAbottom                                    !16 
     &    *A2G_Loop(p21K1,K23,j21f(0,i1polext,1),GGG23(0,i2polext,i3polext),m_bottom)
*----------------------------------------------------------------------
         Abottom3j(17)=cplModAbottom                                   !17
     &               *dotcc(j21iAbottomf(0,SIGMA,i1polext,i2polext,2),epsgC(0,i3polext,3))
*----------------------------------------------------------------------
         Abottom3j(18)=cplModAbottom                                   !18
     &               *dotcc(j21Abottomif(0,SIGMA,i1polext,i2polext,2),epsgC(0,i3polext,3))
*----------------------------------------------------------------------
         Abottom3j(19)=cplModAbottom                                   !19
     &               *dotcc(j21if(0,SIGMA,i1polext,i2polext,2),epsAbottom(0,i3polext,3)) 
*----------------------------------------------------------------------
         Abottom3j(20)=cplModAbottom*GGG_vertex(p21K3p,PHK1p,Kp(0,2),     !20 +
     &                j21i(0,i3polext,3),epsAbottom(0,i1polext,1),epsgC(0,i2polext,2))
*----------------------------------------------------------------------
         Abottom3j(21)=cplModAbottom*GGG_vertex(p21K3p,Kp(0,1),PHK2p,     !21 +
     &                j21i(0,i3polext,3),epsgC(0,i1polext,1),epsAbottom(0,i2polext,2))
*----------------------------------------------------------------------
         Abottom3j(22)=-cplModAbottom*GGG_vertex(p21K3Hp,Kp(0,1),Kp(0,2), !22 minus
     &                j21Abottomi(0,i3polext,3),epsgC(0,i1polext,1),epsgC(0,i2polext,2))
*----------------------------------------------------------------------
         Abottom3j(23)=-cplModAbottom                                  !23 minus
     &              *A3G_Loop(p21K3,K(0,1),K(0,2),j21i(0,i3polext,3),
     &                             epsgC(0,i1polext,1),epsgC(0,i2polext,2),m_bottom) 
*----------------------------------------------------------------------
         Abottom3j(24)=cplModAbottom                                   !24
     &           *A2G_Loop(p21K3,K12,j21i(0,i3polext,3),GGG12(0,i1polext,i2polext),m_bottom)
*----------------------------------------------------------------------
         Abottom3j(25)=-cplModAbottom*GGG_vertex(Kp(0,2),p21p,PHK13p,    !25 minus
     &                epsgC(0,i2polext,2),j21,GGGPAbottom13(0,i1polext,i3polext))
*----------------------------------------------------------------------
         Abottom3j(26)=-cplModAbottom*H4G_effective2(epsAbottom(0,i1polext,1),!26 minus
     &                 epsgC(0,i3polext,3),epsgC(0,i2polext,2),j21)
*----------------------------------------------------------------------
         Abottom3j(27)=cplModAbottom*GGG_vertex(PHK2p,p21p,K13p,         !27  
     &                epsAbottom(0,i2polext,2),j21,GGG13(0,i1polext,i3polext))
*----------------------------------------------------------------------
         Abottom3j(28)=cplModAbottom*H4G_effective2(epsgC(0,i1polext,1),  !28
     &                 epsgC(0,i3polext,3),epsAbottom(0,i2polext,2),j21)
*----------------------------------------------------------------------
         Abottom3j(29)=-cplModAbottom*GGG_vertex(Kp(0,2),p21p,PHK13p,    !29 minus
     &                epsgC(0,i2polext,2),j21,GGG1PAbottom3(0,i1polext,i3polext))
*----------------------------------------------------------------------
         Abottom3j(30)=cplModAbottom*H4G_effective2(epsgC(0,i1polext,1),  !30
     &                 epsAbottom(0,i3polext,3),epsgC(0,i2polext,2),j21)  
*----------------------------------------------------------------------
         Abottom3j(31)=-cplModAbottom                                  !31 minus
     &              *A3G_Loop(K(0,2),p21,K13,epsgC(0,i2polext,2),
     &                             j21,GGG13(0,i1polext,i3polext),m_bottom) 
*----------------------------------------------------------------------
         Abottom3j(32)=cplModAbottom                                   !32
     &        *A3G_Loop(K(0,3),p21K2,K(0,1),
     &           epsgC(0,i3polext,3),GGG2j21(0,i2polext),epsgC(0,i1polext,1),m_bottom) 
*----------------------------------------------------------------------
         Abottom3j(33)=cplModAbottom                                   !33 VZ
     &           *A2G_Loop(p21,K123,j21,GGG132a(0,i1polext,i3polext,i2polext),m_bottom)
*----------------------------------------------------------------------
c         Abottom3j(34)=-cplModAbottom*H4G_effective2(epsgC(0,i1polext,1),
c     &                  epsgC(0,i3polext,3),epsgC(0,i2polext,2),j21Abottom)  !34 minus
         CALL Vertex4g(epsgC(0,i1polext,1), epsgC(0,i3polext,3),
     &        epsgC(0,i2polext,2),test4g)
         Abottom3j(34)=-cplModAbottom*A2G_Loop(p21,K123,j21,test4g,m_bottom)/(K123(4))
*----------------------------------------------------------------------
         Abottom3j(35)=-cplModAbottom     
     &     *A2G_Loop(p21K2,K13,GGG2j21(0,i2polext),GGG13(0,i1polext,i3polext),m_bottom)
*----------------------------------------------------------------------
         Abottom3j(36)=DCMPLX(0d0,0d0) !---> see Pentagons (85-87)
*----------------------------------------------------------------------
         Abottom3j(37)=cplModAbottom*GGG_vertex(p21K3p,PHK1p,Kp(0,2),     !37 
     &           j21f(0,i3polext,3),epsAbottom(0,i1polext,1),epsgC(0,i2polext,2))
*----------------------------------------------------------------------
         Abottom3j(38)=cplModAbottom*GGG_vertex(p21K3p,Kp(0,1),PHK2p,     !38 
     &           j21f(0,i3polext,3),epsgC(0,i1polext,1),epsAbottom(0,i2polext,2))
*----------------------------------------------------------------------
         Abottom3j(39)=-cplModAbottom*GGG_vertex(p21K3Hp,Kp(0,1),Kp(0,2), !39 minus
     &           j21Abottomf(0,i3polext,3),epsgC(0,i1polext,1),epsgC(0,i2polext,2))
*----------------------------------------------------------------------
         Abottom3j(40)=-cplModAbottom                                  !40 minus
     &              *A3G_Loop(K(0,2),p21K3,K(0,1),epsgC(0,i2polext,2),
     &                j21f(0,i3polext,3),epsgC(0,i1polext,1),m_bottom) 
*----------------------------------------------------------------------
         Abottom3j(41)=cplModAbottom                                   !41
     &    *A2G_Loop(p21K3,K12,j21f(0,i3polext,3),GGG12(0,i1polext,i2polext),m_bottom)
*----------------------------------------------------------------------
         Abottom3j(42)=-cplModAbottom*GGG_vertex(p21p,PHK1p,K23p,     !42 +
     &                j21,epsAbottom(0,i1polext,1),GGG32(0,i3polext,i2polext))
*----------------------------------------------------------------------
         Abottom3j(43)=cplModAbottom*H4G_effective2(epsAbottom(0,i1polext,1),!43
     &                       j21,epsgC(0,i2polext,2),epsgC(0,i3polext,3)) 
*----------------------------------------------------------------------
         Abottom3j(44)=-cplModAbottom*GGG_vertex(p21p,Kp(0,1),PHK23p,   !44 +
     &                j21,epsgC(0,i1polext,1),GGGPAbottom23(0,i2polext,i3polext))
*----------------------------------------------------------------------
         Abottom3j(45)=cplModAbottom*H4G_effective2(epsgC(0,i1polext,1),!45 
     &                       j21,epsAbottom(0,i2polext,2),epsgC(0,i3polext,3)) 
*----------------------------------------------------------------------
         Abottom3j(46)=cplModAbottom*GGG_vertex(p21p,Kp(0,1),PHK23p,    !46 
     &                j21,epsgC(0,i1polext,1),GGGPAbottom32(0,i3polext,i2polext))
*----------------------------------------------------------------------
         Abottom3j(47)=-cplModAbottom*H4G_effective2(epsgC(0,i1polext,1),!47 minus
     &                       j21,epsgC(0,i2polext,2),epsAbottom(0,i3polext,3)) 
*----------------------------------------------------------------------
         Abottom3j(48)=cplModAbottom                                   !48
     &              *A3G_Loop(p21,K(0,1),K23,
     &                    j21,epsgC(0,i1polext,1),GGG32(0,i3polext,i2polext),m_bottom) 
*----------------------------------------------------------------------
         Abottom3j(49)=-cplModAbottom                                  !49 minus
     &              *A3G_Loop(K(0,3),K(0,2),p21K1,
     &                epsgC(0,i3polext,3),epsgC(0,i2polext,2),GGGj211(0,i1polext),m_bottom) 
*----------------------------------------------------------------------
         Abottom3j(50)=-cplModAbottom                                  !50 minus
     &     *A2G_Loop(p21,K123,j21,GGG132b(0,i1polext,i3polext,i2polext),m_bottom)
*----------------------------------------------------------------------
c         Abottom3j(51)=-cplModAbottom*H4G_effective2(epsgC(0,i1polext,1),
c     &                  j21Abottom,epsgC(0,i2polext,2),epsgC(0,i3polext,3))  !51 minus
         CALL Vertex4g(epsgC(0,i3polext,3), epsgC(0,i2polext,2),
     &        epsgC(0,i1polext,1),test4g)
         Abottom3j(51)=cplModAbottom*A2G_Loop(p21,K123,j21,test4g,m_bottom)/(K123(4))
*----------------------------------------------------------------------
         Abottom3j(52)=cplModAbottom                                   !52 
     &    *A2G_Loop(p21K1,K23,GGGj211(0,i1polext),GGG32(0,i3polext,i2polext),m_bottom)
*----------------------------------------------------------------------
         Abottom3j(53)=DCMPLX(0d0,0d0) !---> see Pentagons (85-87)
*----------------------------------------------------------------------
         Abottom3j(54)=-cplModAbottom*GGG_vertex(p21p,PHK12p,Kp(0,3),    !54 minus
     &                 j21,GGGPAbottom12(0,i1polext,i2polext),epsgC(0,i3polext,3))
*----------------------------------------------------------------------
         Abottom3j(55)=-cplModAbottom*H4G_effective2(epsAbottom(0,i1polext,1),!55 minus
     &                        epsgC(0,i2polext,2),epsgC(0,i3polext,3),j21) 
*----------------------------------------------------------------------
         Abottom3j(56)=-cplModAbottom*GGG_vertex(p21p,PHK12p,Kp(0,3),    !56 minus
     &                 j21,GGG1PAbottom2(0,i1polext,i2polext),epsgC(0,i3polext,3))
*----------------------------------------------------------------------
         Abottom3j(57)=cplModAbottom*H4G_effective2(epsgC(0,i1polext,1),!43
     &                        epsAbottom(0,i2polext,2),epsgC(0,i3polext,3),j21) 
*---------------------------------------------------------------------- 
         Abottom3j(58)=-cplModAbottom*GGG_vertex(p21p,K12p,PHK3p,      !58 minus
     &                 j21,GGG12(0,i1polext,i2polext),epsAbottom(0,i3polext,3))
*----------------------------------------------------------------------
         Abottom3j(59)=cplModAbottom*H4G_effective2(epsgC(0,i1polext,1),!43
     &                        epsgC(0,i2polext,2),epsAbottom(0,i3polext,3),j21) 
*----------------------------------------------------------------------
         Abottom3j(60)=cplModAbottom                                   !60
     &              *A3G_Loop(K(0,2),p21K3,K(0,1),
     &        epsgC(0,i2polext,2),GGG3j21(0,i3polext),epsgC(0,i1polext,1),m_bottom) 
*----------------------------------------------------------------------
         Abottom3j(61)=cplModAbottom                                   !61 
     &              *A3G_Loop(p21,K12,K(0,3),
     &                      j21,GGG12(0,i1polext,i2polext),epsgC(0,i3polext,3),m_bottom) 
*----------------------------------------------------------------------
         Abottom3j(62)=cplModAbottom                                   !62
     &     *A2G_Loop(p21,K123,j21,GGG123(0,i1polext,i2polext,i3polext),m_bottom)
*----------------------------------------------------------------------
c         Abottom3j(63)=-cplModAbottom*H4G_effective2(epsgC(0,i1polext,1),
c     &                  epsgC(0,i2polext,2),epsgC(0,i3polext,3),j21Abottom)  !63 minus
         CALL Vertex4g(epsgC(0,i1polext,1), epsgC(0,i2polext,2),
     &        epsgC(0,i3polext,3),test4g)
         Abottom3j(63)=-cplModAbottom*A2G_Loop(p21,K123,j21,test4g,m_bottom)/(K123(4))
*----------------------------------------------------------------------
         Abottom3j(64)=-cplModAbottom                                  !64 minus
     &     *A2G_Loop(p21K3,K12,GGG3j21(0,i3polext),GGG12(0,i1polext,i2polext),m_bottom)
*----------------------------------------------------------------------
         Abottom3j(65)=DCMPLX(0d0,0d0) !---> see Pentagons (85-87)
*----------------------------------------------------------------------
         Abottom3j(66)=cplModAbottom                                   !66
     &               *dotcc(j21if(0,SIGMA,i2polext,i3polext,6),epsAbottom(0,i1polext,1)) 
*----------------------------------------------------------------------
         Abottom3j(67)=cplModAbottom                                   !67
     &           *dotcc(j21iAbottomf(0,SIGMA,i2polext,i3polext,5),epsgC(0,i1polext,1)) 
*----------------------------------------------------------------------
         Abottom3j(68)=cplModAbottom                                   !68
     &           *dotcc(j21Abottomif(0,SIGMA,i2polext,i3polext,6),epsgC(0,i1polext,1)) 
*----------------------------------------------------------------------
         Abottom3j(69)=cplModAbottom*GGG_vertex(p21K2p,PHK1p,Kp(0,3),     !69 
     &                j21f(0,i2polext,2),epsAbottom(0,i1polext,1),epsgC(0,i3polext,3))
*----------------------------------------------------------------------
         Abottom3j(70)=-cplModAbottom*GGG_vertex(p21K2Hp,Kp(0,1),Kp(0,3), !70 minus
     &           j21Abottomf(0,i2polext,2),epsgC(0,i1polext,1),epsgC(0,i3polext,3))
*----------------------------------------------------------------------
         Abottom3j(71)=cplModAbottom*GGG_vertex(p21K2p,Kp(0,1),PHK3p,     !71 
     &                j21f(0,i2polext,2),epsgC(0,i1polext,1),epsAbottom(0,i3polext,3))
*----------------------------------------------------------------------
         Abottom3j(72)=-cplModAbottom                                  !72 minus
     &              *A3G_Loop(p21K2,K(0,1),K(0,3),
     &                j21f(0,i2polext,2),epsgC(0,i1polext,1),epsgC(0,i3polext,3),m_bottom) 
*----------------------------------------------------------------------
         Abottom3j(73)=cplModAbottom                                   !73
     &     *A2G_Loop(p21K2,K13,j21f(0,i2polext,2),GGG13(0,i1polext,i3polext),m_bottom)
*----------------------------------------------------------------------
         Abottom3j(74)=cplModAbottom                                   !74
     &        *dotcc(j21Abottomif(0,SIGMA,i3polext,i1polext,4),epsgC(0,i2polext,2))
*----------------------------------------------------------------------
         Abottom3j(75)=cplModAbottom                                   !75
     &        *dotcc(j21if(0,SIGMA,i3polext,i1polext,3),epsAbottom(0,i2polext,2))
*----------------------------------------------------------------------
         Abottom3j(76)=cplModAbottom                                   !76
     &        *dotcc(j21iAbottomf(0,SIGMA,i3polext,i1polext,4),epsgC(0,i2polext,2))
*----------------------------------------------------------------------
         Abottom3j(77)=-cplModAbottom*GGG_vertex(Kp(0,2),Kp(0,3),p21K1Hp, !77 minus
     &                epsgC(0,i2polext,2),epsgC(0,i3polext,3),j21Abottomi(0,i1polext,1))
*----------------------------------------------------------------------
         Abottom3j(78)=cplModAbottom*GGG_vertex(PHK2p,Kp(0,3),p21K1p,     !78
     &                epsAbottom(0,i2polext,2),epsgC(0,i3polext,3),j21i(0,i1polext,1))
*----------------------------------------------------------------------
         Abottom3j(79)=-cplModAbottom*GGG_vertex(Kp(0,2),p21K1p,PHK3p,    !79 minus
     &                epsgC(0,i2polext,2),j21i(0,i1polext,1),epsAbottom(0,i3polext,3))
*----------------------------------------------------------------------
         Abottom3j(80)=cplModAbottom                                   !80 
     &              *A3G_Loop(K(0,3),K(0,2),p21K1,epsgC(0,i3polext,3),
     &                             epsgC(0,i2polext,2),j21i(0,i1polext,1),m_bottom)
*----------------------------------------------------------------------
         Abottom3j(81)=-cplModAbottom                                  !81 minus
     &     *A2G_Loop(p21K1,K23,j21i(0,i1polext,1),GGG32(0,i3polext,i2polext),m_bottom)
*----------------------------------------------------------------------
         Abottom3j(82)=cplModAbottom                                   !82
     &           *dotcc(j21Abottomif(0,SIGMA,i2polext,i1polext,1),epsgC(0,i3polext,3)) 
*----------------------------------------------------------------------
         Abottom3j(83)=cplModAbottom                                   !83
     &           *dotcc(j21iAbottomf(0,SIGMA,i2polext,i1polext,1),epsgC(0,i3polext,3)) 
*----------------------------------------------------------------------
         Abottom3j(84)=cplModAbottom                                   !84
     &           *dotcc(j21if(0,SIGMA,i2polext,i1polext,1),epsAbottom(0,i3polext,3)) 
         
         if (MGdebug) then
            do i=1,87
               GGFLOamp2(i1polext,i2polext,i3polext,i)=Abottom3j(i)*dsqrt(fac)
c               print*,'GGFLOamp2(',i,')',GGFLOamp2(i1polext,i2polext,i3polext,i)
            enddo
         endif
      endif                    ! end CP-even with bottom-loop




cfcccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
cfcc 
cfcc         Check whether there is instabilities
cfcc
cfcc!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cfcc!                                                       !! 
cfcc!         IF SO ::::-->   RSUM -> 0                     !!
cfcc!                                                       !!
cfcc!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
cfcc
cfcc         Count instabilites for pen,box,vertices
cfcc         both in quad and dble precision
cfcc           
cfccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cfccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cfcc 
cfcc   Check Vertices
cfcc     
cfcc        if(boxeva) then
cfcc   Dble precision
cfc           If(verdbleSin) then
cfcc          print*,"verdbleSi",verdbleSin
cfcc          stop
cfcc bad dble
cfc        verinst(gf_iteration,1,1)=verinst(gf_iteration,1,1)+1
cfcc        print*, " verinst(gf_iteration,1,1)", verinst(gf_iteration,1,1)
cfc       else
cfcc good dble 
cfc        verinst(gf_iteration,2,1)=verinst(gf_iteration,2,1)+1
cfcc        print*, " verinst(gf_iteration,2,1)", verinst(gf_iteration,2,1)
cfc        endif
cfc
cfc       If(verQUADSin) then
cfcc bad QUAD
cfc        verinst(gf_iteration,1,2)=verinst(gf_iteration,1,2)+1
cfccccccccccccccccccccccccccccccccccccccccccccccccccccc
cfcc
cfc        RHsum=0d0
cfcc
cfccccccccccccccccccccccccccccccccccccccccccccccccccccc
cfc       else
cfcc good QUAD 
cfc        verinst(gf_iteration,2,2)=verinst(gf_iteration,2,2)+1
cfc        endif
cfcc        endif
cfcc 
cfcc   Check Boxes
cfcc     
cfc        if(boxeva) then
cfcc   Dble precision
cfc       If(boxdbleSin) then
cfcc bad dble
cfc        boxinst(gf_iteration,1,1)=boxinst(gf_iteration,1,1)+1
cfc       else
cfcc good dble 
cfc        boxinst(gf_iteration,2,1)=boxinst(gf_iteration,2,1)+1
cfc        endif
cfc
cfc       If(boxQUADSin) then
cfcc bad QUAD
cfc        boxinst(gf_iteration,1,2)=boxinst(gf_iteration,1,2)+1
cfccccccccccccccccccccccccccccccccccccccccccccccccccccc
cfcc
cfc        RHsum=0d0
cfcc
cfccccccccccccccccccccccccccccccccccccccccccccccccccccc
cfc       else
cfcc good QUAD 
cfc        boxinst(gf_iteration,2,2)=boxinst(gf_iteration,2,2)+1
cfc        endif
cfc        endif
cfccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cfcc 
cfcc   Check Penagons
cfcc     
cfc        if(peneva) then
cfcc   Dble precision
cfc       If(pendbleSin) then
cfcc bad dble
cfc        peninst(gf_iteration,1,1)=peninst(gf_iteration,1,1)+1
cfc       else
cfcc good dble 
cfc        peninst(gf_iteration,2,1)=peninst(gf_iteration,2,1)+1
cfc        endif
cfc
cfc       If(penQUADSin) then
cfcc bad QUAD
cfc        peninst(gf_iteration,1,2)=peninst(gf_iteration,1,2)+1
cfccccccccccccccccccccccccccccccccccccccccccccccccccccc
cfcc
cfc        RHsum=0d0
cfcc
cfccccccccccccccccccccccccccccccccccccccccccccccccccccc
cfc       else
cfcc good QUAD 
cfc        peninst(gf_iteration,2,2)=peninst(gf_iteration,2,2)+1
cfc        endif
cfc        endif
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
         print*,'i1polext,i2polext,i3polext',i1polext,i2polext,i3polext
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

