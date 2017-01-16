C
C******************************  TBV  ********************************
C
C  TBV calculates the WWZ and WWgamma triple boson vertex including 
C  the general 7 anomalous couplings g_1, kappa, lambda, g_4, g_5, 
C  kappa-tilde and lambda-tilde defined by the effective Lagrangian 
C  given in
C
C  K.Hagiwara and D.Zeppenfeld, Nucl. Phys. B282 (1987) 253   ( = HZ)
C
C  INPUT:
C    EPSP(0:3)       complex polarization vector of W^+
C    EPSM(0:3)       complex polarization vector of W^-
C    EPSV(0:3)       complex polarization vector of V = gamma/Z
C                      
C      Q1(0:4)       4-momentum of (virtual) W^+
C      Q2(0:4)       4-momentum of (virtual) W^-
C      Q3(0:4)       4-momentum of (virtual) V
C
C  All 4-momenta are assumed incoming, i.e. Q1 + Q2 + Q3 = 0
C  The fourth components store the square of the Qi:
C  Qi(4) = Qi(0)**2 - Qi(1)**2 - Qi(2)**2 - Qi(3)**2 
C
C  OUTPUT: 
C    RES(2)    complex amplitude for gamma = RES(1)
C              and Z = RES(2) representing 
C              Gamma_V^{\alpha\beta\mu} eps+_\alpha eps-_\beta epsV_\mu
C              as defined in HZ
C
      subroutine tbvini(gaugei)

      implicit none

      double complex  epsp(0:3), epsm(0:3), epsv(0:3), res(2)
      double precision  q1(0:4), q2(0:4), q3(0:4)
      logical gaugei
c
c anomalous couplings:
c here Lambda and nform supply the necessary information on the 
c formfactor which is taken fo Ajj production as 
c
c    1 / ( ( 1 + q1abs/Lambda**2 ) *
c          ( 1 + q2abs/Lambda**2 ) * 
c          ( 1 + q3abs/Lambda**2 ) ) ** nform
c
c where q1abs, q2abs, q3abs are the absolute values of the q**2 of the
c three gauge bosons coupled in the TBV
c The three anomalous inputs lambda, kappa, g can have individual values
c of lambda and nform.  The formfactors and modified parameters for all other processes are 
c calculated using the subroutine anomal_formfactor
c
      double precision  g1(2),kap(2),lam(2)
      double precision  g4(2),g5(2),kaptil(2),lamtil(2)
      double precision  lambda!, lambdaK, lambdaG, lambdaL
      integer  nform!, nformK, nformG, nformL
      common /b_tbv/ g1,kap,lam,g4,g5,kaptil,lamtil,nform
c
c initialize photon couplings
c
      data g1(1)     /1.d0/
      data kap(1)    /1.d0/
      data lam(1)    /0.d0/
      data g4(1)     /0.d0/
      data g5(1)     /0.d0/
      data kaptil(1) /0.d0/
      data lamtil(1) /0.d0/
c
c and Z boson couplings
c
      data g1(2)     /1.d0/
      data kap(2)    /1.d0/
      data lam(2)    /0.d0/
      data g4(2)     /0.d0/
      data g5(2)     /0.d0/
      data kaptil(2) /0.d0/
      data lamtil(2) /0.d0/
c
c and fix the form factor for anomalous couplings
c
!      data lambda    /1.d0/
!      data nform     /1/
c
c other SM couplings
c
#include "VBFNLO/utilities/koppln.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/an_couplings.inc" 
c
c LOCAL variables are
c
      double complex  q1ep,q1em,q1ev, q2ep,q2em,q2ev, q3ep,q3em,q3ev, 
     &                epem,epev,emev, c1,c2,c3,c4, 
     &                spvpm(0:3),spv12(0:3),e1vpm,e2vpm,emv12,epv12
      double precision  q1q2,q1q3, xlam(2),xlamt(2)
      double precision fac
      double precision lamb2i!, lamb2iG, lamb2iK, lamb2iL
      double precision  zero
      parameter( zero = 0.0d0 )
      integer  i
c
c call of TBVINI tests the anomalous couplings supplied by the common 
c block B_TBV for parity violating couplings (to see whether TBV 
c really needs to calculate terms involving the epsilon tensor) and 
c initializes some logical variables to be used in TBV to make the 
c production calls execute faster. In particular for 
c
c GAUGEI = .TRUE.
c
c it is assumed that   EPSP.Q1 = EPSM.Q2 = EPSV.Q3 = 0
c
      double precision  parchc, smc
      logical parity, gauge, sm, ldebug
      parameter (ldebug = .false.)
      data parity/.true./, gauge/.false./, sm/.true./ 
      save parity,gauge,sm,xlam,xlamt


** setting anomalous WWV coupling parameters, converting to tbv notation
      if (with_anom) then
* photon couplings:
         kap(1) = 1D0 + aDkappa0
         lam(1) = lambda0
* Z couplings
         g1(2) = 1d0 + zDg0
         kap(2) = 1d0 + zDkappa0
         lam(2) = lambda0
* formfactor:
         nform = ffexponent
         Lambda = ffmassscale2
      end if
         

      gauge = gaugei
      parchc = zero
      smc = zero

      do i = 1,2
         parchc = parchc + abs(g5(i)) + abs(kaptil(i)) + abs(lamtil(i))
         smc = smc + abs(g1(i)-1.d0) + abs(kap(i)-1.d0) + 
     &               abs(lam(i)) + abs(g4(i))
         xlam(i) = lam(i)/xmw**2
         xlamt(i) = lamtil(i)/xmw**2
      end do

      smc = smc + parchc
      parity = parchc .gt. 1.d-10 
      sm = smc.lt.1.d-10
      if (.not.sm) then
         if (Lambda .gt. 0d0) then
            lamb2i = 1d0/Lambda**2
         else
            lamb2i = 0d0
         end if
         if (ldebug) then
            write(6,*)
            write(6,*) " Values of anomalous couplings:"
            write(6,*) " g1  = ",g1
            write(6,*) " kap = ",kap
            write(6,*) " lam = ",lam
            write(6,*) " g4  = ",g4
            write(6,*) " g5  = ",g5
            write(6,*) " kaptil = ",kaptil
            write(6,*) " lamtil = ",lamtil
            write(6,*)
         end if
      endif
      if (parity) then
         write(6,*) " WARNING in TBV: anomalous coupling section has ",
     &              " not been tested yet for parity odd couplings"
      endif

      if (procID.eq. 150 .and. formfac_ind) then
         print*, "Error: For this process using individual form factors ",
     &           "in anomV.dat is not implemented."
         print*, "Please set FORMFAC_IND to false!"
         stop
      endif

      return

c
c now the main entry of the program
c
      entry tbv(epsp,epsm,epsv,q1,q2,q3,res)

** note: with anom, we have sm=F, gauge=T, parity=F
c
c determine the dot-products
c
      epem = epsp(0)*epsm(0) - epsp(1)*epsm(1) - 
     &       epsp(2)*epsm(2) - epsp(3)*epsm(3)
      epev = epsp(0)*epsv(0) - epsp(1)*epsv(1) - 
     &       epsp(2)*epsv(2) - epsp(3)*epsv(3)
      emev = epsm(0)*epsv(0) - epsm(1)*epsv(1) - 
     &       epsm(2)*epsv(2) - epsm(3)*epsv(3)
c
c first the standard model case
c
      if (sm) then
         q2ev = (q2(0) - q1(0))*epsv(0) - (q2(1)-q1(1))*epsv(1) -
     &          (q2(2) - q1(2))*epsv(2) - (q2(3)-q1(3))*epsv(3)
         q3ep = (q3(0) - q2(0))*epsp(0) - (q3(1)-q2(1))*epsp(1) -
     &          (q3(2) - q2(2))*epsp(2) - (q3(3)-q2(3))*epsp(3)
         q3em = (q3(0) - q1(0))*epsm(0) - (q3(1)-q1(1))*epsm(1) -
     &          (q3(2) - q1(2))*epsm(2) - (q3(3)-q1(3))*epsm(3)
         res(1) = emev*q3ep - epev*q3em + epem*q2ev
         res(2) = res(1)
         return
      endif
c
c next the parity even contributions in the general case


* get formfactor for anomalous couplings 
      if (formfact) then

** formfactor 
         if (procID .eq. 150) then  ! for photon+2jets
            fac = 1.d0/( (1.d0+ abs(q1(4))*lamb2i)*
     &           (1.d0+ abs(q2(4))*lamb2i)*
     &           (1.d0+ abs(q3(4))*lamb2i)  )**nform
         else  ! for massive Vjj
            call anomal_formfactor(q1,q2,q1,q2) ! the last two variables are dummy variables!
            fac = 1D0
         end if

* applying formfactor   
         kap(1) = 1D0 + aDkappa0*fac
         lam(1) = lambda0*fac
         xlam(1) = lam(1)/xmw**2
* Z couplings
         g1(2) = 1d0 + zDg0*fac
         kap(2) = 1d0 + zDkappa0*fac
         lam(2) = lambda0*fac
         xlam(2) = lam(2)/xmw**2
            
      end if


      q1ev = q1(0)*epsv(0)-q1(1)*epsv(1)-q1(2)*epsv(2)-q1(3)*epsv(3)
      q3ep = q3(0)*epsp(0)-q3(1)*epsp(1)-q3(2)*epsp(2)-q3(3)*epsp(3)
      q3em = q3(0)*epsm(0)-q3(1)*epsm(1)-q3(2)*epsm(2)-q3(3)*epsm(3)
      q1q2 = q1(0)*q2(0)  -q1(1)*q2(1)  -q1(2)*q2(2)  -q1(3)*q2(3)

      if (gauge) then
         c1 = -q3em*epev
         c2 = q3ep*emev 
         c3 = q1ev*epem 
         c4 = q1ev*q3ep*q3em

         do i = 1,2
            res(i) = dcmplx(-2.d0+g1(i) + kap(i) + 
     &           xlam(i)*q2(4),-g4(i)) * c1 +
     &           dcmplx(-2.d0 + g1(i) + kap(i) +  
     &           xlam(i)*q1(4), g4(i)) * c2 -
     &           (-2.d0 + 2.d0*g1(i) + xlam(i)*q3(4)) * c3 + 
     &           (2.d0*xlam(i)) * c4
            res(i) = 2.d0*(c1+c2-c3) + res(i)
         end do 

      else
         q2em = q2(0)*epsm(0)-q2(1)*epsm(1)-q2(2)*epsm(2)-q2(3)*epsm(3)
         q3ev = q3(0)*epsv(0)-q3(1)*epsv(1)-q3(2)*epsv(2)-q3(3)*epsv(3)
         q1ep = q1(0)*epsp(0)-q1(1)*epsp(1)-q1(2)*epsp(2)-q1(3)*epsp(3)
         q2ev = -(q1ev+q3ev)
         q2ep = -(q1ep+q3ep)
         q1em = -(q2em+q3em)
         q1q3 = -q1(4) - q1q2 
         c1 = q1em*epev - epem*(q1ev-q2ev) - q2ep*emev
         c2 = q3ep*emev - q3em*epev
         c3 = q1ev*q2ep*q3em-q1em*q2ev*q3ep + 
     &        epem*(q3(4)*q1ev-q1q3*q3ev) -
     &        epev*(q2(4)*q1em-q1q2*q2em) + emev*(q1(4)*q2ep-q1q2*q1ep)
         c4 = q3ep*emev+q3em*epev
     
         do i = 1,2
            res(i) = (g1(i)-1.d0)*c1 + (kap(i)-1.d0)*c2-
     &           xlam(i)*c3 + dcmplx( -g4(i)*dimag(c4),
     &           g4(i)*dreal(c4) )
            res(i) = c1 + c2 + res(i)
         end do
      endif

      if (parity) then
         if (gauge) then
            q1em = -q3em
            q2ep = -q3ep
         endif
         call epsccc(spvpm,epsv,epsp,epsm)
         call epscrr(spv12,epsv,q1,q2)
         e1vpm = q1(0)*spvpm(0)-q1(1)*spvpm(1)-
     &           q1(2)*spvpm(2)-q1(3)*spvpm(3)
         e2vpm = q2(0)*spvpm(0)-q2(1)*spvpm(1)-
     &           q2(2)*spvpm(2)-q2(3)*spvpm(3)
         emv12 = epsm(0)*spv12(0) - epsm(1)*spv12(1) -
     &           epsm(2)*spv12(2) - epsm(3)*spv12(3)
         epv12 = epsp(0)*spv12(0) - epsp(1)*spv12(1) -
     &           epsp(2)*spv12(2) - epsp(3)*spv12(3)
         c1 = q2ep*emv12 + q1em*epv12
         do i = 1,2
            res(i) = res(i) + (xlamt(i)*c1 +
     &           dcmplx(xlamt(i)**q1q2-kaptil(i), 
     &           g5(i))*e1vpm +
     &           dcmplx(xlamt(i)**q1q2-kaptil(i),
     &           -g5(i))*e2vpm)
         end do
      endif
ccc
      return
      end      

