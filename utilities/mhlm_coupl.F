c********************************************************
c
c	these subroutines calculate the couplings of Vector bosons
c	and fermions for the 3 site Higgsless model in the limit that
c	the heavy fermions completely decouple.
c
c	Implementation is taken from the package at:
c	https://hep.pa.msu.edu/twiki/bin/view/ThreeSiteModel/MadGraph
c
c*****************************************************
	subroutine MHLM_coupl(wmass,ee,mwp,zmass)
c*****************************************************
c Input:	wmass: mass of SM w boson
c		ee:	 g_wwa coupling
c		mwp:	 mass of the W'
c		zmass: mass of SM z boson

	implicit none

c	SM fermions are treated massless
c	the heavy partner aquire an infinite mass
c	W' doesnt couple to SM fermions through ideal delocalization
	double precision wmass,ee,mwp,zmass,mzp
      double precision GWWA,GWWZ,gOne,GWWpA,GWWpZ,GWWZp,GWpWpA,
     +GWpWpZ,GWWpZp,GWpWpZp,GWWAA,GWWZA,GWWZZ,GWWWW,
     +GWpAWA,GWpZWA,GWpZWZ,GWpWWW,GWZpWA,GWZpWZ,GWpAWpA,GWpZWpA,
     +GWpZWpZ,GWpWWpW,GWpZpWA,GWpZpWZ,GWZpWZp,GWpZpWpA,GW
     +pZpWpZ,GWpWpWpW,GWpZpWZp,GWpZpWpZp,GWpWpWpWp
	double precision gf,gtf,vw0f,vw1f,vwp0f,vwp1f
	double precision vz0f,vz1f,vzp0f,vzp1f, mzpf

#include "kkcoupl.inc"

	external gf,gtf,vw0f,vw1f,vwp0f,vwp1f,vz0f,vz1f,vzp0f
	external vzp1f, mzpf

c---------------------------------
c the Z' mass
c-----------------------------------
	mzp=MZPf(WMASS,ZMASS,MWP)
c*******************************************
c           VVV
c*******************************************

c-------------------------------------------
c       V-light   V-light   V-light
c-------------------------------------------

      gwwa = ee
      gwwz = 
     - gf(ee,WMASS,ZMASS,MWP)
     -	*vW0f(WMASS,MWP)*vW0f(WMASS,MWP)*vZ0f(WMASS,ZMASS,MWP)
     - +gtf(ee,WMASS,ZMASS,MWP)
     -	*vW1f(WMASS,MWP)*vW1f(WMASS,MWP)*vZ1f(WMASS,ZMASS,MWP)


c-------------------------------------------
c       V-heavy   V-light   V-light
c-------------------------------------------

      GWWpA= 0d0

      GWWpZ=
     - gf(ee,WMASS,ZMASS,MWP)
     -	*vW0f(WMASS,MWP)*vWP0f(WMASS,MWP)*vZ0f(WMASS,ZMASS,MWP)
     - +gtf(ee,WMASS,ZMASS,MWP)
     -	*vW1f(WMASS,MWP)*vWP1f(WMASS,MWP)*vZ1f(WMASS,ZMASS,MWP)


      GWWZp=
     - gf(ee,WMASS,ZMASS,MWP)
     -	*vW0f(WMASS,MWP)*vW0f(WMASS,MWP)*vZP0f(WMASS,ZMASS,MWP)
     - +gtf(ee,WMASS,ZMASS,MWP)
     -	*vW1f(WMASS,MWP)*vW1f(WMASS,MWP)*vZP1f(WMASS,ZMASS,MWP)


c-------------------------------------------
c       V-heavy   V-heavy   V-light
c-------------------------------------------
      GWpWpA= ee


      GWpWpZ=
     - gf(ee,WMASS,ZMASS,MWP)
     -	*vWP0f(WMASS,MWP)*vWP0f(WMASS,MWP)*vZ0f(WMASS,ZMASS,MWP)
     - +gtf(ee,WMASS,ZMASS,MWP)
     -	*vWP1f(WMASS,MWP)*vWP1f(WMASS,MWP)*vZ1f(WMASS,ZMASS,MWP)


      GWWpZp=
     - gf(ee,WMASS,ZMASS,MWP)
     -	*vW0f(WMASS,MWP)*vWP0f(WMASS,MWP)*vZP0f(WMASS,ZMASS,MWP)
     - +gtf(ee,WMASS,ZMASS,MWP)
     -	*vW1f(WMASS,MWP)*vWP1f(WMASS,MWP)*vZP1f(WMASS,ZMASS,MWP)


c-------------------------------------------
c       V-heavy   V-heavy   V-heavy
c-------------------------------------------
      GWpWpZp=
     - gf(ee,WMASS,ZMASS,MWP)
     -	*vWP0f(WMASS,MWP)*vWP0f(WMASS,MWP)*vZP0f(WMASS,ZMASS,MWP)
     - +gtf(ee,WMASS,ZMASS,MWP)
     -	*vWP1f(WMASS,MWP)*vWP1f(WMASS,MWP)*vZP1f(WMASS,ZMASS,MWP)


c*******************************************
c           VVVV
c*******************************************
      gOne = 1d0

c-------------------------------------------
c   V-light    V-light   V-light   V-light
c-------------------------------------------
      GWWAA = ee**2

      
      GWWZA = ee*(
     - gf(ee,WMASS,ZMASS,MWP)
     -	*vW0f(WMASS,MWP)*vW0f(WMASS,MWP)
     -	*vZ0f(WMASS,ZMASS,MWP)
     - +gtf(ee,WMASS,ZMASS,MWP)
     -	*vW1f(WMASS,MWP)*vW1f(WMASS,MWP)
     -	*vZ1f(WMASS,ZMASS,MWP)
     - )


      GWWZZ = 
     - gf(ee,WMASS,ZMASS,MWP)**2
     -	*vW0f(WMASS,MWP)*vW0f(WMASS,MWP)
     -	*vZ0f(WMASS,ZMASS,MWP)*vZ0f(WMASS,ZMASS,MWP)
     - +gtf(ee,WMASS,ZMASS,MWP)**2
     -	*vW1f(WMASS,MWP)*vW1f(WMASS,MWP)
     -	*vZ1f(WMASS,ZMASS,MWP)*vZ1f(WMASS,ZMASS,MWP)

c according to the modified helas routine
c https://hep.pa.msu.edu/twiki/bin/view/ThreeSiteModel/MadGraph
c this is the coupling squared of the WWWW vertex
c the output is set to the sqrt of it, so i can be used with normal helas routines.
c however gwwwwp^2 is negative, and it doesnt work for that, but it isnt used in the matrix-
c element calculation anyway
      GWWWW = 
     - gf(ee,WMASS,ZMASS,MWP)**2
     -	*vW0f(WMASS,MWP)*vW0f(WMASS,MWP)
     -	*vW0f(WMASS,MWP)*vW0f(WMASS,MWP)
     - +gtf(ee,WMASS,ZMASS,MWP)**2
     -	*vW1f(WMASS,MWP)*vW1f(WMASS,MWP)
     -	*vW1f(WMASS,MWP)*vW1f(WMASS,MWP)


c-------------------------------------------
c   V-heavy    V-light   V-light   V-light
c-------------------------------------------
      GWpAWA= 0d0


      GWpZWA= ee*(
     - gf(ee,WMASS,ZMASS,MWP)
     -	*vWP0f(WMASS,MWP)*vW0f(WMASS,MWP)
     -	*vZ0f(WMASS,ZMASS,MWP)
     - +gtf(ee,WMASS,ZMASS,MWP)
     -	*vWP1f(WMASS,MWP)*vW1f(WMASS,MWP)
     -	*vZ1f(WMASS,ZMASS,MWP)
     - )


      GWpZWZ= 
     - gf(ee,WMASS,ZMASS,MWP)**2
     -	*vWP0f(WMASS,MWP)*vW0f(WMASS,MWP)
     -	*vZ0f(WMASS,ZMASS,MWP)*vZ0f(WMASS,ZMASS,MWP)
     - +gtf(ee,WMASS,ZMASS,MWP)**2
     -	*vWP1f(WMASS,MWP)*vW1f(WMASS,MWP)
     -	*vZ1f(WMASS,ZMASS,MWP)*vZ1f(WMASS,ZMASS,MWP)

     
      GWpWWW= 
     - gf(ee,WMASS,ZMASS,MWP)**2
     -	*vWP0f(WMASS,MWP)*vW0f(WMASS,MWP)
     -	*vW0f(WMASS,MWP)*vW0f(WMASS,MWP)
     - +gtf(ee,WMASS,ZMASS,MWP)**2
     -	*vWP1f(WMASS,MWP)*vW1f(WMASS,MWP)
     -	*vW1f(WMASS,MWP)*vW1f(WMASS,MWP)


      GWZpWA= ee*(
     - gf(ee,WMASS,ZMASS,MWP)
     -	*vW0f(WMASS,MWP)*vW0f(WMASS,MWP)
     -	*vZP0f(WMASS,ZMASS,MWP)
     - +gtf(ee,WMASS,ZMASS,MWP)
     -	*vW1f(WMASS,MWP)*vW1f(WMASS,MWP)
     -	*vZP1f(WMASS,ZMASS,MWP)
     - )


      GWZpWZ= 
     - gf(ee,WMASS,ZMASS,MWP)**2
     -	*vW0f(WMASS,MWP)*vW0f(WMASS,MWP)
     -	*vZP0f(WMASS,ZMASS,MWP)*vZ0f(WMASS,ZMASS,MWP)
     - +gtf(ee,WMASS,ZMASS,MWP)**2
     -	*vW1f(WMASS,MWP)*vW1f(WMASS,MWP)
     -	*vZP1f(WMASS,ZMASS,MWP)*vZ1f(WMASS,ZMASS,MWP)


c-------------------------------------------
c   V-heavy    V-heavy   V-light   V-light
c-------------------------------------------
      GWpAWpA= ee**2


      GWpZWpA= ee*(
     - gf(ee,WMASS,ZMASS,MWP)
     -	*vWP0f(WMASS,MWP)*vWP0f(WMASS,MWP)
     -	*vZ0f(WMASS,ZMASS,MWP)
     - +gtf(ee,WMASS,ZMASS,MWP)
     -	*vWP1f(WMASS,MWP)*vWP1f(WMASS,MWP)
     -	*vZ1f(WMASS,ZMASS,MWP)
     - )


      GWpZWpZ= 
     - gf(ee,WMASS,ZMASS,MWP)**2
     -	*vWP0f(WMASS,MWP)*vWP0f(WMASS,MWP)
     -	*vZ0f(WMASS,ZMASS,MWP)*vZ0f(WMASS,ZMASS,MWP)
     - +gtf(ee,WMASS,ZMASS,MWP)**2
     -	*vWP1f(WMASS,MWP)*vWP1f(WMASS,MWP)
     -	*vZ1f(WMASS,ZMASS,MWP)*vZ1f(WMASS,ZMASS,MWP)


      GWpWWpW= 
     - gf(ee,WMASS,ZMASS,MWP)**2
     -	*vWP0f(WMASS,MWP)*vWP0f(WMASS,MWP)
     -	*vW0f(WMASS,MWP)*vW0f(WMASS,MWP)
     - +gtf(ee,WMASS,ZMASS,MWP)**2
     -	*vWP1f(WMASS,MWP)*vWP1f(WMASS,MWP)
     -	*vW1f(WMASS,MWP)*vW1f(WMASS,MWP)


      GWpZpWA= ee*(
     - gf(ee,WMASS,ZMASS,MWP)
     -	*vWP0f(WMASS,MWP)*vW0f(WMASS,MWP)
     -	*vZP0f(WMASS,ZMASS,MWP)
     - +gtf(ee,WMASS,ZMASS,MWP)
     -	*vWP1f(WMASS,MWP)*vW1f(WMASS,MWP)
     -	*vZP1f(WMASS,ZMASS,MWP)
     - )


      GWpZpWZ= 
     - gf(ee,WMASS,ZMASS,MWP)**2
     -	*vWP0f(WMASS,MWP)*vW0f(WMASS,MWP)
     -	*vZP0f(WMASS,ZMASS,MWP)*vZ0f(WMASS,ZMASS,MWP)
     - +gtf(ee,WMASS,ZMASS,MWP)**2
     -	*vWP1f(WMASS,MWP)*vW1f(WMASS,MWP)
     -	*vZP1f(WMASS,ZMASS,MWP)*vZ1f(WMASS,ZMASS,MWP)


      GWZpWZp= 
     - gf(ee,WMASS,ZMASS,MWP)**2
     -	*vW0f(WMASS,MWP)*vW0f(WMASS,MWP)
     -	*vZP0f(WMASS,ZMASS,MWP)*vZP0f(WMASS,ZMASS,MWP)
     - +gtf(ee,WMASS,ZMASS,MWP)**2
     -	*vW1f(WMASS,MWP)*vW1f(WMASS,MWP)
     -	*vZP1f(WMASS,ZMASS,MWP)*vZP1f(WMASS,ZMASS,MWP)


c-------------------------------------------
c   V-heavy    V-heavy   V-heavy   V-light
c-------------------------------------------
      GWpZpWpA= ee*(
     - gf(ee,WMASS,ZMASS,MWP)
     -	*vWP0f(WMASS,MWP)*vWP0f(WMASS,MWP)
     -	*vZP0f(WMASS,ZMASS,MWP)
     - +gtf(ee,WMASS,ZMASS,MWP)
     -	*vWP1f(WMASS,MWP)*vWP1f(WMASS,MWP)
     -	*vZP1f(WMASS,ZMASS,MWP)
     - )


      GWpZpWpZ= 
     - gf(ee,WMASS,ZMASS,MWP)**2
     -	*vWP0f(WMASS,MWP)*vWP0f(WMASS,MWP)
     -	*vZP0f(WMASS,ZMASS,MWP)*vZ0f(WMASS,ZMASS,MWP)
     - +gtf(ee,WMASS,ZMASS,MWP)**2
     -	*vWP1f(WMASS,MWP)*vWP1f(WMASS,MWP)
     -	*vZP1f(WMASS,ZMASS,MWP)*vZ1f(WMASS,ZMASS,MWP)


      GWpWpWpW=
     - gf(ee,WMASS,ZMASS,MWP)**2
     -	*vWP0f(WMASS,MWP)*vWP0f(WMASS,MWP)
     -	*vWP0f(WMASS,MWP)*vW0f(WMASS,MWP)
     - +gtf(ee,WMASS,ZMASS,MWP)**2
     -	*vWP1f(WMASS,MWP)*vWP1f(WMASS,MWP)
     -	*vWP1f(WMASS,MWP)*vW1f(WMASS,MWP)


      GWpZpWZp= 
     - gf(ee,WMASS,ZMASS,MWP)**2
     -	*vWP0f(WMASS,MWP)*vW0f(WMASS,MWP)
     -	*vZP0f(WMASS,ZMASS,MWP)*vZP0f(WMASS,ZMASS,MWP)
     - +gtf(ee,WMASS,ZMASS,MWP)**2
     -	*vWP1f(WMASS,MWP)*vW1f(WMASS,MWP)
     -	*vZP1f(WMASS,ZMASS,MWP)*vZP1f(WMASS,ZMASS,MWP)


c-------------------------------------------
c   V-heavy    V-heavy   V-heavy   V-heavy
c-------------------------------------------
      GWpZpWpZp= 
     - gf(ee,WMASS,ZMASS,MWP)**2
     -	*vWP0f(WMASS,MWP)*vWP0f(WMASS,MWP)
     -	*vZP0f(WMASS,ZMASS,MWP)*vZP0f(WMASS,ZMASS,MWP)
     - +gtf(ee,WMASS,ZMASS,MWP)**2
     -	*vWP1f(WMASS,MWP)*vWP1f(WMASS,MWP)
     -	*vZP1f(WMASS,ZMASS,MWP)*vZP1f(WMASS,ZMASS,MWP)


      GWpWpWpWp= 
     - gf(ee,WMASS,ZMASS,MWP)**2
     -	*vWP0f(WMASS,MWP)*vWP0f(WMASS,MWP)
     -	*vWP0f(WMASS,MWP)*vWP0f(WMASS,MWP)
     - +gtf(ee,WMASS,ZMASS,MWP)**2
     -	*vWP1f(WMASS,MWP)*vWP1f(WMASS,MWP)
     -	*vWP1f(WMASS,MWP)*vWP1f(WMASS,MWP)


c--------------------------------------------------------------
c now the output for the kk_coupl_inp routine
c--------------------------------------------------------------
      open (unit=50,file="kk_coupl_inp.dat")
	write(50,*) "KKMAXW =              1"
	write(50,*) "KKMAXZ =              1"
	write(50,*) "KKMAXG =              0"
	write(50,*) "--------------------------------------------"
	write(50,*) "--------------------------------------------"

	write(50,*) "---Kaluza-Klein masses---"
	write(50,*) "KK_MASSW =", mwp
	write(50,*) "KK_MASSZ =", mzp
	write(50,*) "---------------------------------------------"
	write(50,*) "---------------------------------------------"

	write(50,*)"-- RECALL: W0 = SM W, Z0 = SM Z, Z'0 = SM photon"
	write(50,*) "---WWZ coupl---"
	write(50,*) "CPL_W0-W0-ZX =", gwwz, gwwzp, "!W0-W0- Z0 to Z1"
*	write(50,*) gwwzp
	write(50,*) "CPL_W1-W0-ZX =", gwwpz, gwwpzp, "!W1-W0- Z0 to Z1"
*	write(50,*) gwwpzp
	write(50,*) "CPL_W1-W1-ZX =", gwpwpz, gwpwpzp, "!W1-W1- Z0 to Z1"
*	write(50,*) gwpwpzp	
	write(50,*) "---------------------------------------------"

	write(50,*) "--WWZ' coupl----"
	write(50,*) "CPL_W1-W0-GX = ", gwwpa, "!W1-W0- Z'0 to Z'0"
	write(50,*) "CPL_W1-W1-GX = ", gwpwpa, "!W1-W1- Z'0 to Z'0"
	write(50,*) "---------------------------------------------"

	write(50,*) "--WWWW coupl----"
c like written above in the calculation, the output is set to the root of it.
	write(50,*) "CPL_W0-W0-W0-WX =", sqrt(gwwww), "!W0-W0-W0- W0 to W0"
	write(50,*) "CPL_W1-W0-W0-WX =", gwpwww, "!W1-W0-W0- W0 to W0"
	write(50,*) "CPL_W1-W1-W0-WX =", gwpwwpw, "!W1-W1-W0- W0 to W0"
	write(50,*) "CPL_W1-W1-W1-WX =", gwpwpwpw, gwpwpwpwp,"!W1-W1-W1- W0 to W1"
*	write(50,*) gwpwpwpwp
	write(50,*) "---------------------------------------------"

	write(50,*) "--WWZZ coupl----"
	write(50,*) "CPL_W0-W0-Z0-ZX =", gwwzz, "!W0-W0-Z0 -Z0 to Z0"
	write(50,*) "CPL_W0-W0-Z1-ZX =", gwzpwz, gwzpwzp, "!W0-W0-Z1 -Z0 to Z1"
*	write(50,*) gwzpwzp
	write(50,*) "CPL_W1-W0-Z0-ZX =", gwpzwz, "!W1-W0-Z0 -Z0 to Z0"
	write(50,*) "CPL_W1-W0-Z1-ZX =", gwpzpwz, gwpzpwzp, "!W1-W0-Z1 -Z0 to Z1"
*	write(50,*) gwpzpwzp
	write(50,*) "CPL_W1-W1-Z0-ZX =", gwpzwpz, "!W1-W1-Z0 -Z0 to Z0"
	write(50,*) "CPL_W1-W1-Z1-ZX =", gwpzpwpz, gwpzpwpzp, "!W1-W1-Z1 -Z0 to Z1"
*	write(50,*) gwpzpwpzp
	write(50,*) "---------------------------------------------"

	write(50,*)"--WWZ'Z' coupl----"
	write(50,*) "CPL_W0-W0-G0-GX =", gwwaa, "!W0-W0-Z'0 -Z'0 to Z'0"
	write(50,*) "CPL_W1-W0-G0-GX =", gwpawa, "!W1-W0-Z'0 -Z'0 to Z'0"
	write(50,*) "CPL_W1-W1-G0-GX =", gwpawpa, "!W1-W1-Z'0 -Z'0 to Z'0"
	write(50,*) "---------------------------------------------"

	write(50,*) "--WWZ'Z coupl----"
	write(50,*) "CPL_W0-W0-G0-ZX =", gwwza, gwzpwa, "!W0-W0-Z'0- Z0 to Z1"
*	write(50,*) gwzpwa
	write(50,*) "CPL_W1-W0-G0-ZX =", gwpzwa, gwpzpwa, "!W1-W0-Z'0- Z0 to Z1"
*	write(50,*) gwpzpwa
	write(50,*) "CPL_W1-W1-G0-ZX =", gwpzwpa, gwpzpwpa, "!W1-W1-Z'0- Z0 to Z1"
*	write(50,*) gwpzpwpa
	close(50)
	return
	END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	Fermion couplings in 3-site higgsless model            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine mhlm_fermion(mwp)

c	Only interactions of Fermions with the light SM gauge bosons is
c	taken into account. There are no couplings between SM gauge bosons
c	and heavy fermions as well as heavy gauge bosons and light Sm fermions.
c	also the 3rd Generation fermions have the same coupling as all other
c	fermions

	implicit none

#include "coupl.inc"
	double precision mwp
	double precision two,three,zero,pi
	parameter (pi = 3.141592653589793d0)
	parameter (zero=0d0)
	parameter (two=2d0)
	parameter (three=3d0)
	double precision gf,gtf,gpf,vw0f,vw1f
	double precision vz0f,vz1f,vz2f
	double precision vL0f,vL1f,gwfL,gznl,gzll,gzlr
	double precision gzdr,gzdl,gzul,gzur

	DOUBLE PRECISION CLR,XM2,XMG,B,V,A
	COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)
      double precision BWNE,BWUD,BWTB,BZNN,BZEE,BZUU,BZDD,BZTT,
     &                BHWW,BHZZ,BHGG,BHTT,BHBB,BHCC,BHTAU,BHMU,
     &                BHGAM, BHGAMZ, XGW, XGZ, XGH, GAMT
      COMMON /BRANCH/ BWNE,BWUD,BWTB,BZNN,BZEE,BZUU,BZDD,BZTT,
     &                BHWW,BHZZ,BHGG,BHTT,BHBB,BHCC,BHTAU,BHMU,
     &                BHGAM, BHGAMZ, XGW, XGZ, XGH, GAMT

	double precision decw,gwne,gwud,decz,gznn,
     & 			gzee,gzuu,gzdd
	external gf,gtf,gpf,vw0f,vw1f,vz0f,vz1f,vz2f
	external vL0f,vL1f

c-------------------------------------------
c     V-light   F-light    F-light   (SM)
c-------------------------------------------

	
      gal(1) = dcmplx(  ee          , Zero )
      gal(2) = dcmplx(  ee          , Zero )
      gau(1) = dcmplx( -ee*Two/Three, Zero )
      gau(2) = dcmplx( -ee*Two/Three, Zero )
      gad(1) = dcmplx(  ee/Three    , Zero )
      gad(2) = dcmplx(  ee/Three    , Zero )

      gwfL = (gf(-ee,WMASS,ZMASS,MWP)*
     -     vL0f(WMASS,MWP)**2*vW0f(WMASS,MWP) + 
     -    gtf(-ee,WMASS,ZMASS,MWP)*
     -     vL1f(WMASS,MWP)**2*vW1f(WMASS,MWP))/
     -  sqrt(2d0)
      gwf(1) = dcmplx(  gwfL              , Zero )
      gwf(2) = dcmplx(  Zero              , Zero )


      gznL = 
     - 1d0/2d0*gf(-ee,WMASS,ZMASS,MWP)
     -     *vZ0f(WMASS,ZMASS,MWP)
     -     *vL0f(WMASS,MWP)**2
     - +1d0/2d0*gtf(-ee,WMASS,ZMASS,MWP)
     -     *VZ1f(WMASS,ZMASS,MWP)
     -     *vL1f(WMASS,MWP)**2
     - -1d0/2d0*gpf(-ee,WMASS,ZMASS,MWP)
     -     *vZ2f(WMASS,ZMASS,MWP)
      gzn(1) = dcmplx(  gznL                        , Zero )
      gzn(2) = dcmplx(  Zero                        , Zero )

      gzlL = 
     - -1d0/2d0*gf(-ee,WMASS,ZMASS,MWP)
     -     *vZ0f(WMASS,ZMASS,MWP)
     -     *vL0f(WMASS,MWP)**2
     - -1d0/2d0*gtf(-ee,WMASS,ZMASS,MWP)
     -     *vZ1f(WMASS,ZMASS,MWP)
     -     *vL1f(WMASS,MWP)**2
     - -1d0/2d0*gpf(-ee,WMASS,ZMASS,MWP)
     -     *vZ2f(WMASS,ZMASS,MWP)
      gzlR =
     - -gpf(-ee,WMASS,ZMASS,MWP)
     -     *vZ2f(WMASS,ZMASS,MWP)
      gzl(1) = dcmplx(  gzlL                        , Zero )
      gzl(2) = dcmplx(  gzlR                        , Zero )

      gzuL = 
     - 1d0/2d0*gf(-ee,WMASS,ZMASS,MWP)
     - 	*vZ0f(WMASS,ZMASS,MWP)*vL0f(WMASS,MWP)**2
     - +1d0/2d0*gtf(-ee,WMASS,ZMASS,MWP)
     - 	*VZ1f(WMASS,ZMASS,MWP)*vL1f(WMASS,MWP)**2
     - +1d0/6d0*gpf(-ee,WMASS,ZMASS,MWP)
     - 	*vZ2f(WMASS,ZMASS,MWP)
      gzuR = 
     - 2d0/3d0*gpf(-ee,WMASS,ZMASS,MWP)
     - 	*vZ2f(WMASS,ZMASS,MWP)		
      gzu(1) = dcmplx( gzuL                         , Zero )
      gzu(2) = dcmplx( gzuR                         , Zero )

      gzdL = 
     - -1d0/2d0*gf(-ee,WMASS,ZMASS,MWP)
     - 	*vZ0f(WMASS,ZMASS,MWP)*vL0f(WMASS,MWP)**2
     - -1d0/2d0*gtf(-ee,WMASS,ZMASS,MWP)
     - 	*VZ1f(WMASS,ZMASS,MWP)*vL1f(WMASS,MWP)**2
     - +1d0/6d0*gpf(-ee,WMASS,ZMASS,MWP)
     - 	*vZ2f(WMASS,ZMASS,MWP)
      gzdR = 
     - -1d0/3d0*gpf(-ee,WMASS,ZMASS,MWP)
     -	*vZ2f(WMASS,ZMASS,MWP)
      gzd(1) = dcmplx( gzdL                         , Zero )
      gzd(2) = dcmplx( gzdR                         , Zero )

	clr(1,3,-1)=-dreal(gwf(1))	! W+- boson
	clr(1,3, 1)=-dreal(gwf(2))	! zero
	clr(2,3,-1)=-dreal(gwf(1))	! W+- boson
	clr(2,3, 1)=-dreal(gwf(2))	! zero
	clr(3,3,-1)=-dreal(gwf(1))	! W+- boson
	clr(3,3, 1)=-dreal(gwf(2))	! zero
	clr(4,3,-1)=-dreal(gwf(1))	! W+- boson
	clr(4,3, 1)=-dreal(gwf(2))	! zero

	clr(1,4,-1)=-dreal(gwf(1))	! W+- boson
	clr(1,4, 1)=-dreal(gwf(2))	! zero
	clr(2,4,-1)=-dreal(gwf(1))	! W+- boson
	clr(2,4, 1)=-dreal(gwf(2))	! zero
	clr(3,4,-1)=-dreal(gwf(1))	! W+- boson
	clr(3,4, 1)=-dreal(gwf(2))	! zero
	clr(4,4,-1)=-dreal(gwf(1))	! W+- boson
	clr(4,4, 1)=-dreal(gwf(2))	! zero

	clr(1,2,-1)=-dreal(gzn(1))	! Z nu
	clr(1,2, 1)=-dreal(gzn(2))	! zero	
	clr(2,2,-1)=-dreal(gzl(1))	! Z leptons left
	clr(2,2, 1)=-dreal(gzl(2))	! Z leptons right
	clr(3,2,-1)=-dreal(gzu(1))	! Z u left
	clr(3,2, 1)=-dreal(gzu(2))	! Z u right
	clr(4,2,-1)=-dreal(gzd(1))	! Z d left
	clr(4,2, 1)=-dreal(gzd(2))	! Z d right

c recalculate widths of Z and W in mhlm:
      DECW = wmass/12.d0/PI
      GWNE = DECW*( clr(1,3,-1)**2 + clr(1,3,1)**2 )/2d0
      GWUD = DECW*( clr(1,3,-1)**2 + clr(1,3,1)**2 )/2d0*3.d0*(1.d0 + ALFAS/PI)
      XGW  = 3.d0*GWNE + 2.d0*GWUD
      BWNE = GWNE/XGW
      BWUD = GWUD/XGW
      BWTB = 0d0

      DECZ = zmass/12.d0/PI
      GZNN = DECZ*( clr(1,2,-1)**2 + clr(1,2,1)**2 )/2d0
      GZEE = DECZ*( clr(2,2,-1)**2 + clr(2,2,1)**2 )/2d0
      GZUU = DECZ*( clr(3,2,-1)**2 + clr(3,2,1)**2 )/2d0*3.d0*(1.d0 + ALFAS/PI)
      GZDD = DECZ*( clr(4,2,-1)**2 + clr(4,2,1)**2 )/2d0*3.d0*(1.d0 + ALFAS/PI)

      XGZ  = 3.d0*GZNN + 3.d0*GZEE + 2.d0*GZUU + 3.d0*GZDD
      BZNN = GZNN/XGZ
      BZEE = GZEE/XGZ
      BZUU = GZUU/XGZ
      BZDD = GZDD/XGZ
      BZTT = 0d0

      XMG(2) = zmass*XGZ
      XMG(3) = wmass*XGW
      XMG(4) = XMG(3)

	wwidth = xmg(3)/wmass
	zwidth = xmg(2)/zmass

      print*,"*******************************************************"
	print*, "   NEW BRANCHING RATIOS AND WIDTHS FOR W AND Z BOSON"
	print*, ""
	print*, "values above are obsolete due to new fermion couplings"
      print*,"*******************************************************"
      write(*,"(T4,A,T36,A,F10.6,A)") "total widths of Z0",": ",
     &     XGZ," GeV"
      write(*,"(T4,A,T36,A,F10.6,A)") "total widths of W+ or W-",": ",
     &     XGW," GeV"
      write(*,"(T4,A,T36,A,F10.7,A)") "W -> lepton anti-lepton",": ",
     &     BWNE*100.d0," %"
      write(*,"(T4,A,T36,A,F10.7,A)") 
     &     "W -> massless quark anti-quark",": ", BWUD*100.d0," %"
      write(*,"(T4,A,T36,A,F10.7,A)") "W -> top bottom",": ",
     &     BWTB*100.d0," %"
      print *," "
      write(*,"(T4,A,T36,A,F10.7,A)") 
     &     "Z -> neutrino anti-neutrino",": ",BZNN*100.d0," %"
      write(*,"(T4,A,T36,A,F10.7,A)") 
     &     "Z -> charged lepton anti-lepton",": ",BZEE*100.d0," %"
      write(*,"(T4,A,T36,A,F10.7,A)") 
     &     "Z -> up-type quark anti-quark",": ",BZUU*100.d0," %"
      write(*,"(T4,A,T36,A,F10.7,A)") 
     &     "Z -> down-type quark anti-quark",": ",BZDD*100.d0," %"
      write(*,"(T4,A,T36,A,F10.7,A)") "Z -> top anti-top",": ",
     &     BZTT*100.d0," %"


	RETURN
	END


	subroutine setmgtest
	implicit none
#include "VBFNLO/utilities/coupl.inc"
#include "kkcoupl.inc"
#include "mgtest.inc"

      double precision GZpDL, GZpDR, zero1, zero, zpclr(1:5,-1:1)
      double precision GZpUL, GZpUR, GZpLL, GZpLR, xgzp
      double precision GZpNL, GZpNR, gzptl, gzptr, MF
	double precision gf,vZP0f,gtf,VZP1f,VZP2f,gpf, heavyzp
	double precision vL0f,vL1f,vL2f, vRt1f, vRt2f, vLt0f, vLt1f
      double precision deczp, gzpuu, gzpdd, gzpee, gzpnn, gzptt, pi
      double precision bzpuu, bzpdd, bzpee, bzpnn, bzptt, x1,y1,y2, bzpww
      parameter (pi = 3.141592653589793d0)
      logical lmgteston
      parameter (lmgteston = .true.)
	external gf,vZP0f,gtf,VZP1f,VZP2f,gpf
	external vL0f,vL1f,vL2f, vRt1f, vRt2f, vLt0f, vLt1f!, vl

	print*, "RUN USING MADGRAPH MATRIX ELEMENT"
c MF set to fixed value of 4 TeV
      MF=4000.0
c transfer values for couplings, masses, widths to MG matrix element variables      
	zero1=0d0
	zero=0d0
	gwwzp=gwwz_k(0,0,1)
	gwwpz=gwwz_k(1,0,0)
	gwwpzp=gwwz_k(1,0,1)
	gwpwpzp=gwwz_k(1,1,1)
	gwpwpz=gwwz_k(1,1,0)
	mzp=z1mass_k(1)
	mwp=w1mass_k(1)
	wzp=z1width_k(1)
	wwp=w1width_k(1)
	gwwaa=gwwaa_k(0,0,0,0)
	gwzpwz=gwwzz_k(0,0,1,0)
	gwzpwzp=gwwzz_k(0,0,1,1)
	gwwzz=gwwzz_k(0,0,0,0)
	gwzpwa=gwwaz_k(0,0,0,1)
	gwwza=gwwaz_k(0,0,0,0)
	gwwzpa=gwwaz_k(0,0,0,1)
	gwwww=gwwww_k(0,0,0,0)

	if(lmgteston) then
      GZpDL = 
     - -1d0/2d0*gf(-EE,WMASS,ZMASS,MWP)
     - 	*vZP0f(WMASS,ZMASS,MWP)*vL0f(WMASS,MWP)**2
     - -1d0/2d0*gtf(-EE,WMASS,ZMASS,MWP)
     - 	*VZP1f(WMASS,ZMASS,MWP)*vL1f(WMASS,MWP)**2
     - +1d0/6d0*gpf(-EE,WMASS,ZMASS,MWP)
     - 	*vZP2f(WMASS,ZMASS,MWP)
     	GZpDR = 
     - -1d0/3d0*gpf(-EE,WMASS,ZMASS,MWP)
     -	*vZP2f(WMASS,ZMASS,MWP)
      GZpD(1)=dcmplx(GZpDL,Zero)
      GZpD(2)=dcmplx(GZpDR,Zero)
c      GZpD(1)=dcmplx(Zero1,Zero)
c      GZpD(2)=dcmplx(Zero1,Zero)
      GZpUL = 
     - 1d0/2d0*gf(-EE,WMASS,ZMASS,MWP)
     - 	*vZP0f(WMASS,ZMASS,MWP)*vL0f(WMASS,MWP)**2
     - +1d0/2d0*gtf(-EE,WMASS,ZMASS,MWP)
     - 	*VZP1f(WMASS,ZMASS,MWP)*vL1f(WMASS,MWP)**2
     - +1d0/6d0*gpf(-EE,WMASS,ZMASS,MWP)
     - 	*vZP2f(WMASS,ZMASS,MWP)
     	GZpUR = 
     - +2d0/3d0*gpf(-EE,WMASS,ZMASS,MWP)
     - 	*vZP2f(WMASS,ZMASS,MWP)		
      GZpU(1)=dcmplx(GZpUL,Zero)
      GZpU(2)=dcmplx(GZpUR,Zero)
c      GZpU(1)=dcmplx(Zero1,Zero)
c      GZpU(2)=dcmplx(Zero1,Zero)

      GZpLL = 
     - -1d0/2d0*gf(-EE,WMASS,ZMASS,MWP)
     - 	*vZP0f(WMASS,ZMASS,MWP)*vL0f(WMASS,MWP)**2
     - -1d0/2d0*gtf(-EE,WMASS,ZMASS,MWP)
     - 	*VZP1f(WMASS,ZMASS,MWP)*vL1f(WMASS,MWP)**2
     - -1d0/2d0*gpf(-EE,WMASS,ZMASS,MWP)
     - 	*vZP2f(WMASS,ZMASS,MWP)
     	GZpLR = 
     - -1d0*gpf(-EE,WMASS,ZMASS,MWP)
     - 	*vZP2f(WMASS,ZMASS,MWP)
      GZpL(1)=dcmplx(GZpLL,Zero)
      GZpL(2)=dcmplx(GZpLR,Zero)
c      GZpL(1)=dcmplx(Zero1,Zero)
c      GZpL(2)=dcmplx(Zero1,Zero)

      GZpNL = 
     - 1d0/2d0*gf(-EE,WMASS,ZMASS,MWP)
     - 	*vZP0f(WMASS,ZMASS,MWP)*vL0f(WMASS,MWP)**2
     - +1d0/2d0*gtf(-EE,WMASS,ZMASS,MWP)
     - 	*VZP1f(WMASS,ZMASS,MWP)*vL1f(WMASS,MWP)**2
     - -1d0/2d0*gpf(-EE,WMASS,ZMASS,MWP)
     - 	*vZP2f(WMASS,ZMASS,MWP)
     	GZpNR = 0d0
      GZpN(1)=dcmplx(GZpNL,Zero)
      GZpN(2)=dcmplx(GZpNR,Zero)
c      GZpN(1)=dcmplx(Zero1,Zero)
c      GZpN(2)=dcmplx(Zero,Zero)

      GZpTL = 
     - 1d0/2d0*gf(-ee,WMASS,ZMASS,MWP)
     - 	*vZP0f(WMASS,ZMASS,MWP)*vLt0f(WMASS,MWP,TMASS,MF)**2
     - +1d0/2d0*gtf(-ee,WMASS,ZMASS,MWP)
     - 	*VZP1f(WMASS,ZMASS,MWP)*vLt1f(WMASS,MWP,TMASS,MF)**2
     - +1d0/6d0*gpf(-ee,WMASS,ZMASS,MWP)
     - 	*vZP2f(WMASS,ZMASS,MWP)
     	GZpTR = 
     - 1d0/2d0*gtf(-ee,WMASS,ZMASS,MWP)
     - 	*VZP1f(WMASS,ZMASS,MWP)*vRt1f(WMASS,MWP,TMASS,MF)**2
     - +1d0/6d0*gpf(-ee,WMASS,ZMASS,MWP)
     - 	*VZP2f(WMASS,ZMASS,MWP)*vRt1f(WMASS,MWP,TMASS,MF)**2
     - +2d0/3d0*gpf(-ee,WMASS,ZMASS,MWP)
     - 	*vZP2f(WMASS,ZMASS,MWP)*vRt2f(WMASS,MWP,TMASS,MF)**2		
      GZpT(1)=dcmplx(GZpTL,Zero)
      GZpT(2)=dcmplx(GZpTR,Zero)

c calculate new width of Z' including fermion couplings

	zpclr(1,-1)=-dreal(gzpn(1))	! Zp nu left
	zpclr(1, 1)=-dreal(gzpn(2))	! Zp nu right
	zpclr(2,-1)=-dreal(gzpl(1))	! Zp leptons left
	zpclr(2, 1)=-dreal(gzpl(2))	! Zp leptons right
	zpclr(3,-1)=-dreal(gzpu(1))	! Zp u left
	zpclr(3, 1)=-dreal(gzpu(2))	! Zp u right
	zpclr(4,-1)=-dreal(gzpd(1))	! Zp d left
	zpclr(4, 1)=-dreal(gzpd(2))	! Zp d right
	zpclr(5,-1)=-dreal(gzpt(1)) ! Zp t left
	zpclr(5, 1)=-dreal(gzpt(2)) ! Zp t right

      DECZP = z1mass_k(1)/12.d0/PI
      GZPNN = DECZP*( zpclr(1,-1)**2 + zpclr(1,1)**2 )/2d0
      GZPEE = DECZP*( zpclr(2,-1)**2 + zpclr(2,1)**2 )/2d0
      GZPUU = DECZP*( zpclr(3,-1)**2 + zpclr(3,1)**2 )/2d0*3.d0*(1.d0 + ALFAS/PI)
      GZPDD = DECZP*( zpclr(4,-1)**2 + zpclr(4,1)**2 )/2d0*3.d0*(1.d0 + ALFAS/PI)

c now the decay into heavy top pairs, same as in koppln.F
      x1=( zpclr(5,-1)+zpclr(5,1) )/( zpclr(5,1)-zpclr(5,-1) )
      y1=tmass/z1mass_k(1)
      y2=tmass/z1mass_k(1)

      HEAVYZP=   ( 1.d0 - 0.5d0*(Y1**2 + Y2**2) - 
     &           0.5d0*(Y1**2 - Y2**2)**2 +
     &           3.d0*Y1*Y2*((X1**2 - 1.d0)/(X1**2 + 1.d0)) )*
     &           SQRT((1.d0 - Y1**2 - Y2**2)**2 - 4.d0*Y1**2*Y2**2)

      GZPTT = GZPUU*HEAVYZP

c now sum it all up
      XGZP  = 3.d0*GZPNN + 3.d0*GZPEE + 2.d0*GZPUU + 3.d0*GZPDD + GZPTT + z1width_k(1)

      BZPWW = z1width_k(1)/XGZP
      BZPNN = GZPNN/XGZP
      BZPEE = GZPEE/XGZP
      BZPUU = GZPUU/XGZP
      BZPDD = GZPDD/XGZP
      BZPTT = GZPTT/XGZP

      print*,"*******************************************************"
	print*,"        NEW BRANCHING RATIOS FOR Z' BOSON"
      print*,"*******************************************************"
      write(*,"(T4,A,T36,A,F10.6,A)") "total widths of Z'",": ",
     &     XGZP," GeV"
      write(*,"(T4,A,T36,A,F10.7,A)") 
     &     "Z' -> W+ W-",": ",BZPWW*100.d0," %"
      write(*,"(T4,A,T36,A,F10.7,A)") 
     &     "Z' -> neutrino anti-neutrino",": ",BZPNN*100.d0," %"
      write(*,"(T4,A,T36,A,F10.7,A)") 
     &     "Z' -> charged lepton anti-lepton",": ",BZPEE*100.d0," %"
      write(*,"(T4,A,T36,A,F10.7,A)") 
     &     "Z' -> up-type quark anti-quark",": ",BZPUU*100.d0," %"
      write(*,"(T4,A,T36,A,F10.7,A)") 
     &     "Z' -> down-type quark anti-quark",": ",BZPDD*100.d0," %"
      write(*,"(T4,A,T36,A,F10.7,A)") "Z' -> top anti-top",": ",
     &     BZPTT*100.d0," %"
      
      z1width_k(1)=xgzp
      wzp=z1width_k(1)

	endif

	Return
	end


