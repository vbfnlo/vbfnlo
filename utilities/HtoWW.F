********************************************************************************
********************************************************************************
*** HtoWW.F                                                                  ***
*** 27 July 2011                                                             ***
*** sophy@particle.uni-karlsruhe.de, rauch@particle.uni-karlsruhe.de         ***
***                                                                          ***
*** This file contains the code needed to compute the H -> VV - 4 lepton     ***
*** matrix squared elements.  These were originally duplicated in the        ***
*** folders amplitudes/vvjj and amplitudes/hjjj, and have now been moved     ***
*** here.  Used for Hjj via VBF and gluon fusion, Hjjj and HAjj.             ***
***                                                                          ***
********************************************************************************
********************************************************************************

c     F1,F2,F3
      subroutine computeFs(qsq1,qsq2,tau1,tau3,flav1,flav3,
     $     bos1,bos2,F1,F2,F3)
      implicit none
      real*8 qsq1,qsq2          !q1^2 and q2^2
      integer tau1,tau3         !helicities of f1 and f3
      integer flav1,flav3       !flavor of fermions
      integer bos1,bos2         !boson type : W+ 3
c                                             W- 4
c                                             Z  2 
c                                         photon 1
c
      double complex prop1(4),prop2(4) !1 = photon propagator
c     2 = Z propagator
      double complex F1,F2,F3
c     koppln commom blocks
      double complex ahvv(3,4,4), ahvvL(3,4,4)
      common/tensorhvv/ ahvv, ahvvL
      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)


** propagator factors:
      prop1(1) = 1.d0/dcmplx(qsq1,0.0d0)
      prop2(1) = 1.d0/dcmplx(qsq2,0.0d0)
      prop1(2) = 1.d0/dcmplx(qsq1-xm2(2),xmg(2))
      prop2(2) = 1.d0/dcmplx(qsq2-xm2(2),xmg(2))
      prop1(3) = 1.d0/dcmplx(qsq1-xm2(3),xmg(3))
      prop2(3) = 1.d0/dcmplx(qsq2-xm2(3),xmg(3))
      prop1(4) = prop1(3)
      prop2(4) = prop2(3)

      if(bos1.eq.2) then  ! H -> ZZ -> llll
         F1 = CLR(flav1,2,tau1)*CLR(flav3,2,tau3)*ahvv(1,2,2)*
     1                              prop1(2)*prop2(2)
         
** NOTE: this includes effective H-gamma-Z and H-gamma-gamma vertices in the SM
         F2 = CLR(flav1,2,tau1)*CLR(flav3,2,tau3)*ahvv(2,2,2)*
     1                              prop1(2)*prop2(2)+
     $       V(flav1,1)*CLR(flav3,2,tau3)*ahvv(2,2,1)*prop1(1)*prop2(2)+
     $       V(flav3,1)*CLR(flav1,2,tau1)*ahvv(2,2,1)*prop1(2)*prop2(1)+
     $       V(flav1,1)*V(flav3,1)*ahvv(2,1,1)*prop1(1)*prop2(1)
c
         F3 = CLR(flav1,2,tau1)*CLR(flav3,2,tau3)*ahvv(3,2,2)*
     1                              prop1(2)*prop2(2)+
     $       V(flav1,1)*CLR(flav3,2,tau3)*ahvv(3,2,1)*prop1(1)*prop2(2)+
     $       V(flav3,1)*CLR(flav1,2,tau1)*ahvv(3,2,1)*prop1(2)*prop2(1)+
     $       V(flav1,1)*V(flav3,1)*ahvv(3,1,1)*prop1(1)*prop2(1)

      else  ! H -> WW -> lvlv

         F1 = CLR(flav1,3,tau1)*CLR(flav3,4,tau3)*ahvv(1,3,4)*
     1                              prop1(3)*prop2(4)
c         
         F2 = CLR(flav1,3,tau1)*CLR(flav3,4,tau3)*ahvv(2,3,4)*
     1                              prop1(3)*prop2(4)
c     
         F3 = CLR(flav1,3,tau1)*CLR(flav3,4,tau3)*ahvv(3,3,4)*
     1                              prop1(3)*prop2(4)

      endif

      return

      end
c


********************************************************************************
********************************************************************************


c     sum over chiralities tau1,tau3 =-1,+1
c     
      subroutine m2s_VVsum_had(L,bos1,bos2,flav1,flav3,m2s)

c     flav1,flav3 are flavors of vector boson decay products
c     possible values:
c     1 = neutrino
c     2 = lepton
c     3 = up-type quark
c     4 = down-type quark
c     5 = down + up - type quark (sum over possibilities)
c
c     possible combinations for hadronic decay:
c     

      implicit none
      integer bos1, bos2, flav1, flav3, i, helicity_up(4)
      real*8 L(0:3,4), m2s, m2s_up, m2s_down, rnumb, rnumb2, RandomNumber, QSQ_had
      external RandomNumber
      real*8 mjj2
      external mjj2
      double precision alphas5, alphas5_hardwired
      external alphas5, alphas5_hardwired

      integer pdflib, valueLHA
      COMMON/PDFparameters/pdflib, valueLHA

#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/process.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"

      m2s = 0.0d0
      m2s_up = 0.0d0
      m2s_down = 0.0d0

      if (flav1.ne.5 .and. flav3.ne.5) then     ! W-decay
         call m2s_VVsum(L,bos1,bos2,flav1,flav3,m2s)
         if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94)
     &        m2s = m2s * 2.d0                  ! factor 2 when both quark combinations selected
      elseif (flav1.eq.5) then                  ! Z-decay
         call m2s_VVsum(L,bos1,bos2,3,flav3,m2s_up)
         if ((lha.or.hepmc) .and. (.not. doNLO)) then
            do i=1,4
              helicity_up(i)= helicity(i)
            enddo
         endif
         if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94)
     &        m2s_up = m2s_up * 2.d0                             ! factor 2 when all quark combinations selected
         call m2s_VVsum(L,bos1,bos2,4,flav3,m2s_down)
         if (finalquarks(1).eq.93) m2s_down = m2s_down * 2.d0    ! factor 2 when 2 quark generations selected
         if (finalquarks(1).eq.94) m2s_down = m2s_down * 3.d0    ! factor 3 when 3 quark generations selected
         m2s = m2s_up + m2s_down
      else
         print*, "In Hjj_ZZhad the first Z has to decay hadronically!"
         stop
      endif


c...Les Houches interface - set decay particles for hadronic decay
      if ((lha.or.hepmc).and..not.doNLO) then
         SELECT CASE(process)
         CASE(Hjj_WPhadWM)
            if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
               rnumb=RandomNumber()
               finalquarks_psp(1) =  2 + 2* INT(rnumb*2)
               finalquarks_psp(2) = -1 - 2* INT(rnumb*2)
            endif
         CASE(Hjj_WPWMhad)
            if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
               rnumb=RandomNumber()
               finalquarks_psp(1) =  1 + 2* INT(rnumb*2)
               finalquarks_psp(2) = -2 - 2* INT(rnumb*2)
            endif
         CASE(Hjj_ZZhad)
            rnumb2 = RandomNumber()
            if (finalquarks(1).eq.93) then
               rnumb=RandomNumber()
               if (m2s_up.le.(rnumb2*m2s)) then    ! up-type
                 finalquarks_psp(1) =  2 + 2* INT(rnumb*2)
                 finalquarks_psp(2) = -2 - 2* INT(rnumb*2)
                 do i=1,4
                   helicity(i)= helicity_up(i)
                 enddo
               else                       ! down-type
                 finalquarks_psp(1) =  1 + 2* INT(rnumb*2)
                 finalquarks_psp(2) = -1 - 2* INT(rnumb*2)
               endif
            elseif (finalquarks(1).eq.94) then
               rnumb=RandomNumber()
               if (m2s_up.le.(rnumb2*m2s)) then    ! up-type
                 finalquarks_psp(1) =  2 + 2* INT(rnumb*2)
                 finalquarks_psp(2) = -2 - 2* INT(rnumb*2)
                 do i=1,4
                   helicity(i)= helicity_up(i)
                 enddo
               else                       ! down-type
                 finalquarks_psp(1) =  1 + 2* INT(rnumb*3)
                 finalquarks_psp(2) = -1 - 2* INT(rnumb*3)
               endif
            endif
         END SELECT
      endif


      ! include NLO decay corrections here
      if (nlo_semilep.gt.0) then
         SELECT CASE(process)
         CASE(Hjj_WPhadWM, Hjj_ZZhad)
           QSQ_had = abs(mjj2(L(0,1),L(0,2)))
         CASE(Hjj_WPWMhad)
           QSQ_had = abs(mjj2(L(0,3),L(0,4)))
         end select
         if (pdflib.eq.1 .and. .not.doNLO) then   ! lhapdf used and LO production
           m2s = m2s * (1.d0 + alphas5_hardwired(max(QSQ_had,1d0),1)/pi)
         else
           m2s = m2s * (1.d0 + alphas5(max(QSQ_had,1d0),1)/pi)
         endif
      endif


      return
      end


********************************************************************************
********************************************************************************

c     sum over chiralities tau1,tau3 =-1,+1
c     
      subroutine m2s_VVsum(L,bos1,bos2,flav1,flav3,m2s)
      implicit none
      integer bos1, bos2, flav1, flav3, i, j, h_index, ie, iu
      real*8 L(0:3,4), m2s, temp_m2s, m2s_contrib(4), rnumb, RandomNumber
      external RandomNumber
#include "VBFNLO/utilities/lha.inc"

      m2s = 0.0d0

      if ((bos1 .eq. 3) .or. (bos1 .eq. 4)) then  ! H -> WW -> 4l

         call m2s_VV(L,bos1,bos2,-1,-1,flav1,flav3,temp_m2s)
         m2s = m2s + temp_m2s
         if ((lha.or.hepmc) .and. (.not. doNLO)) then
            helicity(1)=-1
            helicity(2)= 1
            helicity(3)=-1
            helicity(4)= 1
         end if

      else  ! H -> ZZ -> 4l

         do i = -1,1,2
            do j = -1,1,2
               call m2s_VV(L,bos1,bos2,i,j,flav1,flav3,temp_m2s)
               m2s = m2s + temp_m2s
               h_index = (i+1) + (j+1)/2 + 1 !0,2 + 0,1 + 1 = 1..4
               m2s_contrib(h_index) = temp_m2s
            enddo
         enddo
         if ((lha.or.hepmc) .and. (.not. doNLO)) then
            ! choose random helicity combination, based on xsection as weight
            h_index=0
            temp_m2s=0.d0
            rnumb=RandomNumber()
            do while((h_index.lt.4).and.(temp_m2s.le.rnumb*m2s))
               h_index = h_index + 1
               temp_m2s = temp_m2s + m2s_contrib(h_index)
            enddo
            ie = (-1)**((h_index+1)/2)
            iu = (-1)**(h_index)
            helicity(1)= ie
            helicity(2)=-ie
            helicity(3)= iu
            helicity(4)=-iu
         end if

      end if

      return
      end


********************************************************************************
********************************************************************************

c These subroutines compute the matrix element squared for the decays
c h -> W*W* and h-> (Z/gamma)*(Z/gamma)* for
c anomalous hvv couplings

      subroutine m2s_VV(L,bos1,bos2,tau1,tau3,flav1,flav3,m2s)

      implicit none

      real*8 L(0:3,4)           !lepton 4-mom
      real*8 m2s                !matrix element squared
c     
      real*8 colorfac(4)        !colorfac(flav)
      real*8 Q(0:4,2),p13,p24,p32,p14,p12,p34
      real*8 levi       
      real*8 res(0:3)
c      levi = eps-mu-nu-rho-sigma l(mu,1)l(mu,2)l(mu,3)l(mu,4)
c     Q(mu,i) mu is spacetime index
c     i is either 1 or 2
      complex*16 F1,F2,F3,F1c,F2c,F3c,i_
      parameter(i_ = (0.0d0,1.0d0))

c     tau1,tau3 are +1 or -1 chirality
c     bos1,bos2 is boson type 3 for w and 2 for z/gamma
c     flav1,flav3 are flavors of vector boson decay products
      integer tau1,tau3,flav1,flav3,bos1,bos2,mu,i

      complex*16 temp_m2s

c     set up color factors
      colorfac(1) = 1.0d0
      colorfac(2) = 1.0d0
      colorfac(3) = 3.0d0
      colorfac(4) = 3.0d0
c
c     set up dot products 
      p13 = l(0,1)*l(0,3) - l(1,1)*l(1,3) - l(2,1)*l(2,3)- l(3,1)*l(3,3)
      p24 = l(0,2)*l(0,4) - l(1,2)*l(1,4) - l(2,2)*l(2,4)- l(3,2)*l(3,4)
c
      p32 = l(0,2)*l(0,3) - l(1,2)*l(1,3) - l(2,2)*l(2,3)- l(3,2)*l(3,3)
      p14 = l(0,1)*l(0,4) - l(1,1)*l(1,4) - l(2,1)*l(2,4)- l(3,1)*l(3,4)
      p12 = l(0,2)*l(0,1) - l(1,2)*l(1,1) - l(2,2)*l(2,1)- l(3,2)*l(3,1) 
      p34 = l(0,3)*l(0,4) - l(1,3)*l(1,4) - l(2,3)*l(2,4)- l(3,3)*l(3,4)

c     compute levi-civita symbol
      call EPSRRR(res,l(0,2),l(0,3),l(0,4))
      levi = res(0)*l(0,1)-res(1)*l(1,1)-res(2)*l(2,1)-res(3)*l(3,1)

      do mu = 0,3
         Q(mu,1) = l(mu,1) + l(mu,2)
         Q(mu,2) = l(mu,3) + l(mu,4)
      enddo
c     compute qsq for vector bosons
      do i = 1,2
         Q(4,i) = Q(0,i)**2 -  Q(1,i)**2 - Q(2,i)**2 - Q(3,i)**2
      enddo

      call computeFs(Q(4,1),Q(4,2),tau1,tau3,flav1,flav3,
     $     bos1,bos2,F1,F2,F3)

      F1c = conjg(F1)
      F2c = conjg(F2)
      F3c = conjg(F3)

      temp_m2s =
     &  + f3*f3c * (  - 16*p14**2*p32**2 + 32*p13*p14*p32*p24 - 16*
     &    p13**2*p24**2 + 8*p12*p24**2*p34 + 8*p12*p32**2*p34 + 16*p12*
     &    p14*p32*p34 + 8*p12*p14**2*p34 + 16*p12*p13*p24*p34 + 8*p12*
     &    p13**2*p34 - 16*p12**2*p34**2 + 8*tau1*tau3*p12*p24**2*p34 - 
     &    8*tau1*tau3*p12*p32**2*p34 + 16*tau1*tau3*p12*p14*p32*p34 - 8
     &    *tau1*tau3*p12*p14**2*p34 - 16*tau1*tau3*p12*p13*p24*p34 + 8*
     &    tau1*tau3*p12*p13**2*p34 )

      temp_m2s = temp_m2s + f3*f2c * ( 16*levi*p14*p32 - 16*levi*p13*
     &    p24 - 16*tau1*tau3*levi*p12*p34 - 8*i_*tau3*p12*p24**2*p34 + 
     &    8*i_*tau3*p12*p32**2*p34 - 8*i_*tau3*p12*p14**2*p34 + 8*i_*
     &    tau3*p12*p13**2*p34 - 8*i_*tau1*p12*p24**2*p34 - 8*i_*tau1*
     &    p12*p32**2*p34 + 8*i_*tau1*p12*p14**2*p34 + 8*i_*tau1*p12*
     &    p13**2*p34 )

      temp_m2s = temp_m2s + f3*f1c * (  - 4*levi*p24 + 4*levi*p32 + 4*
     &    levi*p14 - 4*levi*p13 - 4*tau1*tau3*levi*p24 - 4*tau1*tau3*
     &    levi*p32 - 4*tau1*tau3*levi*p14 - 4*tau1*tau3*levi*p13 + 4*i_
     &    *tau3*p14*p32*p24 + 4*i_*tau3*p14*p32**2 - 4*i_*tau3*p14**2*
     &    p32 - 4*i_*tau3*p13*p24**2 - 4*i_*tau3*p13*p32*p24 + 4*i_*
     &    tau3*p13*p14*p24 - 4*i_*tau3*p13*p14*p32 + 4*i_*tau3*p13**2*
     &    p24 - 4*i_*tau3*p12*p24*p34 + 4*i_*tau3*p12*p32*p34 - 4*i_*
     &    tau3*p12*p14*p34 + 4*i_*tau3*p12*p13*p34 + 4*i_*tau1*p14*p32*
     &    p24 - 4*i_*tau1*p14*p32**2 + 4*i_*tau1*p14**2*p32 - 4*i_*tau1
     &    *p13*p24**2 + 4*i_*tau1*p13*p32*p24 - 4*i_*tau1*p13*p14*p24
     &     - 4*i_*tau1*p13*p14*p32 + 4*i_*tau1*p13**2*p24 - 4*i_*tau1*
     &    p12*p24*p34 - 4*i_*tau1*p12*p32*p34 + 4*i_*tau1*p12*p14*p34
     &     + 4*i_*tau1*p12*p13*p34 )

      temp_m2s = temp_m2s + f2*f3c * ( 16*levi*p14*p32 - 16*levi*p13*
     &    p24 - 16*tau1*tau3*levi*p12*p34 + 8*i_*tau3*p12*p24**2*p34 - 
     &    8*i_*tau3*p12*p32**2*p34 + 8*i_*tau3*p12*p14**2*p34 - 8*i_*
     &    tau3*p12*p13**2*p34 + 8*i_*tau1*p12*p24**2*p34 + 8*i_*tau1*
     &    p12*p32**2*p34 - 8*i_*tau1*p12*p14**2*p34 - 8*i_*tau1*p12*
     &    p13**2*p34 )

      temp_m2s = temp_m2s + f2*f2c * ( 16*p14**2*p32**2 - 32*p13*p14*
     &    p32*p24 + 16*p13**2*p24**2 + 8*p12*p24**2*p34 + 8*p12*p32**2*
     &    p34 - 16*p12*p14*p32*p34 + 8*p12*p14**2*p34 - 16*p12*p13*p24*
     &    p34 + 8*p12*p13**2*p34 + 16*p12**2*p34**2 + 8*tau1*tau3*p12*
     &    p24**2*p34 - 8*tau1*tau3*p12*p32**2*p34 - 16*tau1*tau3*p12*
     &    p14*p32*p34 - 8*tau1*tau3*p12*p14**2*p34 + 16*tau1*tau3*p12*
     &    p13*p24*p34 + 8*tau1*tau3*p12*p13**2*p34 )

      temp_m2s = temp_m2s + f2*f1c * (  - 4*p14*p32*p24 + 4*p14*p32**2
     &     + 4*p14**2*p32 + 4*p13*p24**2 - 4*p13*p32*p24 - 4*p13*p14*
     &    p24 - 4*p13*p14*p32 + 4*p13**2*p24 + 4*p12*p24*p34 + 4*p12*
     &    p32*p34 + 4*p12*p14*p34 + 4*p12*p13*p34 - 4*tau1*tau3*p14*p32
     &    *p24 - 4*tau1*tau3*p14*p32**2 - 4*tau1*tau3*p14**2*p32 + 4*
     &    tau1*tau3*p13*p24**2 + 4*tau1*tau3*p13*p32*p24 + 4*tau1*tau3*
     &    p13*p14*p24 - 4*tau1*tau3*p13*p14*p32 + 4*tau1*tau3*p13**2*
     &    p24 + 4*tau1*tau3*p12*p24*p34 - 4*tau1*tau3*p12*p32*p34 - 4*
     &    tau1*tau3*p12*p14*p34 + 4*tau1*tau3*p12*p13*p34 - 4*i_*tau3*
     &    levi*p24 - 4*i_*tau3*levi*p32 + 4*i_*tau3*levi*p14 + 4*i_*
     &    tau3*levi*p13 - 4*i_*tau1*levi*p24 + 4*i_*tau1*levi*p32 - 4*
     &    i_*tau1*levi*p14 + 4*i_*tau1*levi*p13 )

      temp_m2s = temp_m2s + f1*f3c * (  - 4*levi*p24 + 4*levi*p32 + 4*
     &    levi*p14 - 4*levi*p13 - 4*tau1*tau3*levi*p24 - 4*tau1*tau3*
     &    levi*p32 - 4*tau1*tau3*levi*p14 - 4*tau1*tau3*levi*p13 - 4*i_
     &    *tau3*p14*p32*p24 - 4*i_*tau3*p14*p32**2 + 4*i_*tau3*p14**2*
     &    p32 + 4*i_*tau3*p13*p24**2 + 4*i_*tau3*p13*p32*p24 - 4*i_*
     &    tau3*p13*p14*p24 + 4*i_*tau3*p13*p14*p32 - 4*i_*tau3*p13**2*
     &    p24 + 4*i_*tau3*p12*p24*p34 - 4*i_*tau3*p12*p32*p34 + 4*i_*
     &    tau3*p12*p14*p34 - 4*i_*tau3*p12*p13*p34 - 4*i_*tau1*p14*p32*
     &    p24 + 4*i_*tau1*p14*p32**2 - 4*i_*tau1*p14**2*p32 + 4*i_*tau1
     &    *p13*p24**2 - 4*i_*tau1*p13*p32*p24 + 4*i_*tau1*p13*p14*p24
     &     + 4*i_*tau1*p13*p14*p32 - 4*i_*tau1*p13**2*p24 + 4*i_*tau1*
     &    p12*p24*p34 + 4*i_*tau1*p12*p32*p34 - 4*i_*tau1*p12*p14*p34
     &     - 4*i_*tau1*p12*p13*p34 )

      temp_m2s = temp_m2s + f1*f2c * (  - 4*p14*p32*p24 + 4*p14*p32**2
     &     + 4*p14**2*p32 + 4*p13*p24**2 - 4*p13*p32*p24 - 4*p13*p14*
     &    p24 - 4*p13*p14*p32 + 4*p13**2*p24 + 4*p12*p24*p34 + 4*p12*
     &    p32*p34 + 4*p12*p14*p34 + 4*p12*p13*p34 - 4*tau1*tau3*p14*p32
     &    *p24 - 4*tau1*tau3*p14*p32**2 - 4*tau1*tau3*p14**2*p32 + 4*
     &    tau1*tau3*p13*p24**2 + 4*tau1*tau3*p13*p32*p24 + 4*tau1*tau3*
     &    p13*p14*p24 - 4*tau1*tau3*p13*p14*p32 + 4*tau1*tau3*p13**2*
     &    p24 + 4*tau1*tau3*p12*p24*p34 - 4*tau1*tau3*p12*p32*p34 - 4*
     &    tau1*tau3*p12*p14*p34 + 4*tau1*tau3*p12*p13*p34 + 4*i_*tau3*
     &    levi*p24 + 4*i_*tau3*levi*p32 - 4*i_*tau3*levi*p14 - 4*i_*
     &    tau3*levi*p13 + 4*i_*tau1*levi*p24 - 4*i_*tau1*levi*p32 + 4*
     &    i_*tau1*levi*p14 - 4*i_*tau1*levi*p13 )

      ! SM part:
      temp_m2s = temp_m2s + f1*f1c * ( 8*p14*p32 + 8*p13*p24 - 8*tau1*
     &    tau3*p14*p32 + 8*tau1*tau3*p13*p24 )

      m2s = dble(temp_m2s) 
c      print*,"temp_m2s=",temp_m2s
      m2s = m2s * colorfac(flav1)*colorfac(flav3)
c      print*,"m2s=",m2s

      return

      end


********************************************************************************
