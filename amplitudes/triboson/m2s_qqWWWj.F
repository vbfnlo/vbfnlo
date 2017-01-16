*****************************************************************************
c
c    begin subroutine m2s_qqWWWj
c
c*****************************************************************************
c
c       Stefan Prestel, <sprestel@particle.uni-karlsruhe.de>
c	Initial version:  2008
c	Last modified by Vera Hankele:  2009 January 28
c
c       this subroutine computes |ME|^2 for qq -> WWW g
c       and all the crossed processes including the gluon induced
c       subprocesses at LO and the virtual corrections via call of "qqWWWj.f"
c
c INPUT     bos       Boson identifier, 2 = Z, 3,4=W
c           nlo       NLO=1: create subtraction term; LO = 0
c           lok       need to calculate m2s(1:3)? T or F 
c           xi        Feynman x parameters of incoming fermions
c           p         momenta associated with external fermions
c           v         momenta associated with WWZ decay fermions
c           rn        random number array
c           xuz       x_ija for subtraction terms
c           NW        NW=1: Narrow width approximation, NW=0: full calculation
c
c
c OUTPUT    m2s       |M|^2*pdf1*pdf2 
c
c***************************************************************************

      SUBROUTINE m2s_qqWWWj(bos,nlo,lok,xi,p,v,rn,xuz,m2s,NW)

      IMPLICIT NONE

c declare input/output variables
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/process.inc"

      integer N_gen_W
      real*8 fac_W, mjj2
      external mjj2

      integer bos
      double precision xi(nx), p(0:3,max_p,max_kin), v(0:3,max_v,max_kin), 
     1       rn(1), xuz, m2s(0:3)
      logical lok(3)
      integer NW

      double complex zero
      parameter (zero=(0d0,0d0))

c declare external functions
      integer FL_VVVg
      external FL_VVVg

c alfas, scales etc
      double precision x1,x2
      double precision lnQomu(2:3), omxi(2:3), 
     &       Ax(2:3), Bx(2:3), Cx(2:3), Dxi(2:3), tgs2oqsq(2:3), 
     &       ln1mxi, lnrat, z, lnz

c declare local variables
      double precision q_sf
      double precision dotrr
      external dotrr

      integer init/0/, I, J, mu
      save init

c declare variables for summation process
      INTEGER ires, FSIGN(3+n_v), gsign

      integer physToDiag(3), nlo, nmin, nmax
      save nmin, nmax

c store contributions from subprocess ip in res(ip,ID) where
c  ID = 1  : the real emission |M|^2 * pdf
c     = 2,3: subtraction terms for emision off upper or lower line
c     = 0  : subtracted result which drives integration, i.e
c res(*,0) = res(*,1)+res(*,2)+res(*,3)

      double precision res(maxnumsubproc,0:3)

c variables for leshouches interface
      double precision weight,rnumb,RandomNumber

c declare parton distribution variables
      double precision pdf(-6:6,2,3)

c define program switches
      logical ldebug
      parameter (ldebug=.false.)

      double precision swmomenta(0:3,6)

c---------------------------- BEGIN CODE --------------------------------

      if ( init .ne. 0) then
         do i = 1,nmax
            do j=0,3
               res(i,j) = 0
            enddo
         enddo
      else if (init.eq.0) then
            write(6,*) " "
         if (bos.eq.343) then
            write(6,*) "W+ W- W+ j amplitude square information:"
         elseif(bos.eq.434) then
            write(6,*) "W- W+ W- j amplitude square information:"
         else
            write(6,*) "Process not implemented, check boson identifier"
            stop
         endif
         write(6,*) "-----------------------------------------------"
         write(6,*) " "
         call vtowww_reset

         call printnfl(.true.)

         print*," "
         print*," creal = ",crealVV," cvirtual = ",cvirtVV
         print*," "
         init = 1
         do i = 1,maxnumsubproc
            do j = 0,3
               res(i,j) = 0
            enddo
         enddo
      endif

c -------------------------------------------------------------------

c Les Houches interface
      if((lha.or.hepmc).and.(.not.doNLO)) then
         helicity(1) = -1
         helicity(2) =  1
         helicity(3) = -1
         helicity(4) =  1
         helicity(5) = -1
         helicity(6) =  1
      endif

c switch momenta so that electron spinor is always created out of v(0,1,j)
      if(bos.eq.434) then
         do j=1,3
            do mu=0,3
            swmomenta(mu,1) = v(mu,3,j) ! ve
            swmomenta(mu,2) = v(mu,4,j) ! e+
            swmomenta(mu,3) = v(mu,1,j) ! mu-
            swmomenta(mu,4) = v(mu,2,j) ! vm~
            enddo
            do i=1,4
              do mu=0,3
                 v(mu,i,j) = swmomenta(mu,i)
              enddo
            enddo
         enddo
      endif

c lepton spinors and W+, W-, W+/- polarisation vectors
      do J = 1,3

         if(bos.eq.343) then
            CALL OXXXXX(v(0,1,J),ZERO ,-1, 1,wve1(1,J))                  !ve
            CALL IXXXXX(v(0,2,J),ZERO , 1,-1,wep1(1,J))                  !e+
            CALL OXXXXX(v(0,3,J),ZERO ,-1, 1,wmu1(1,J))                  !mu-
            CALL IXXXXX(v(0,4,J),ZERO , 1,-1,wvm1(1,J))                  !vm~
            CALL OXXXXX(v(0,5,J),ZERO ,-1, 1,wvt(1,J))                   !vt
            CALL IXXXXX(v(0,6,J),LMASS, 1,-1,wtp(1,J))                   !ta+
            CALL JIOXXX(wep1(1,J),wve1(1,J),GWF,WMASS,WWIDTH, wp1(1,J))  !W+(1)
            CALL JIOXXX(wvm1(1,J),wmu1(1,J),GWF,WMASS,WWIDTH, w1m(1,J))  !W-(1)
            CALL JIOXXX(wtp(1,J),wvt(1,J),GWF,WMASS,WWIDTH, wp2(1,J))    !W+(2)
            do mu = 0,3
               qp1(mu,J)     = v(mu,1,J)+v(mu,2,J)       !W+ #1
               qp2(mu,J)     = v(mu,5,J)+v(mu,6,J)       !W+ #2
               q1m(mu,J)     = v(mu,3,J)+v(mu,4,J)       !W- #1
               qwp1wm(mu,J)  = qp1(mu,J)+ q1m(mu,J)      !momentum W+#1 + W-#1
               qwp2wm(mu,J)  = qp2(mu,J)+ q1m(mu,J)      !momentum W+#2 + W-#1
               qwpwmwp(mu,J) = qwp1wm(mu,J) + qp2(mu,J)  !momentum W+W-W+
               qwwwp(mu,J)   = qwp1wm(mu,J) + qp2(mu,J) + p(mu,3,J)
            enddo
            qp1(4,J)     = dotrr(qp1(0,J),qp1(0,J))
            qp2(4,J)     = dotrr(qp2(0,J),qp2(0,J))
            q1m(4,J)     = dotrr(q1m(0,J),q1m(0,J))
            qwp1wm(4,J)  = dotrr(qwp1wm(0,J),qwp1wm(0,J))
            qwp2wm(4,J)  = dotrr(qwp2wm(0,J),qwp2wm(0,J))
            qwpwmwp(4,J) = dotrr(qwpwmwp(0,J),qwpwmwp(0,J))
            qwwwp(4,J)   = dotrr(qwwwp(0,J),qwwwp(0,J))

         else if(bos.eq.434) then
            CALL OXXXXX(v(0,1,J),ZERO ,-1, 1,wve2(1,J))                  !ve
            CALL IXXXXX(v(0,2,J),ZERO , 1,-1,wep2(1,J))                  !e+
            CALL OXXXXX(v(0,3,J),ZERO ,-1, 1,wmu2(1,J))                  !mu-
            CALL IXXXXX(v(0,4,J),ZERO , 1,-1,wvm2(1,J))                  !vm~
            CALL OXXXXX(v(0,5,J),LMASS,-1, 1,wtm(1,J))                   !ta-
            CALL IXXXXX(v(0,6,J),ZERO , 1,-1,wvtbar(1,J))                !vt~
            CALL JIOXXX(wep2(1,J),wve2(1,J),GWF,WMASS,WWIDTH, w1p(1,J))  !W+(1)
            CALL JIOXXX(wvm2(1,J),wmu2(1,J),GWF,WMASS,WWIDTH, wm1(1,J))  !W-(1)
            CALL JIOXXX(wvtbar(1,J),wtm(1,J),GWF,WMASS,WWIDTH,wm2(1,J))  !W-(2)
            do mu = 0,3
               q1p(mu,J)     = v(mu,1,J)+v(mu,2,J)       !W+ #1
               qm1(mu,J)     = v(mu,3,J)+v(mu,4,J)       !W- #1
               qm2(mu,J)     = v(mu,5,J)+v(mu,6,J)       !W- #2
               qwpwm1(mu,J)  = q1p(mu,J)+ qm1(mu,J)      !momentum W+#1 + W-#1
               qwpwm2(mu,J)  = q1p(mu,J)+ qm2(mu,J)      !momentum W+#1 + W-#2
               qwpwmwm(mu,J) = qwpwm1(mu,J) + qm2(mu,J) !momentum of W+W-W- 
               qwwwm(mu,J)   = qwpwm1(mu,J) + qm2(mu,J) + p(mu,3,J)
            enddo
            q1p(4,J)     = dotrr(q1p(0,J),q1p(0,J))
            qm1(4,J)     = dotrr(qm1(0,J),qm1(0,J))
            qm2(4,J)     = dotrr(qm2(0,J),qm2(0,J))
            qwpwm1(4,J)  = dotrr(qwpwm1(0,J),qwpwm1(0,J))
            qwpwm2(4,J)  = dotrr(qwpwm2(0,J),qwpwm2(0,J))
            qwpwmwm(4,J) = dotrr(qwpwmwm(0,J),qwpwmwm(0,J))
            qwwwm(4,J)   = dotrr(qwwwm(0,J),qwwwm(0,J))

         else
            write(6,*) "Process not implemented, check boson identifier"
         endif

c calculate leptonic tensors


           SELECT CASE(procid)
           CASE(WPhadWMWP, WPWMhadWP, WMhadWPWM, WMWPhadWM)
               if(bos.eq.343)then
                  call anomal_formfactor(qwpwmwp(0,J),qp1(0,J),q1m(0,J),qp2(0,J))
                  call WWW_atoww1_had(v(0,1,J),J,343,WWW_aww(0,J,1),NW)               !A -> ve  e+  mu-  vm~
                  call WWW_atoww2_had(v(0,1,J),J,343,WWW_aww(0,J,2),NW)               !A -> mu- vm~ vt   tau+
                  call WWW_ztoww1_had(v(0,1,J),J,343,WWW_zww(0,J,1),NW)               !Z -> ve  e+  mu-  vm~
                  call WWW_ztoww2_had(v(0,1,J),J,343,WWW_zww(0,J,2),NW)               !Z -> mu- vm~ vt   tau+
                  call WWW_wptowww_had(v(0,1,J),J,1,wptowww(0,J),NW)                  !W+ -> ve e+ mu- vm~ vt ta+
               else if(bos.eq.434) then
                  call anomal_formfactor(qwpwmwm(0,J),q1p(0,J),qm1(0,J),qm2(0,J))
                  call WWW_atoww1_had(v(0,1,J),J,434,WWW_aww(0,J,1),NW)               !A -> ve  e+  mu-  vm~
                  call WWW_atoww3_had(v(0,1,J),J,434,WWW_aww(0,J,3),NW)               !A -> ve  e+  tau- vt~
                  call WWW_ztoww1_had(v(0,1,J),J,434,WWW_zww(0,J,1),NW)               !Z -> ve  e+  mu-  vm~
                  call WWW_ztoww3_had(v(0,1,J),J,434,WWW_zww(0,J,3),NW)               !Z -> ve  e+  tau- vt~
                  call WWW_wmtowww_had(v(0,1,J),J,1,wmtowww(0,J),NW)                  !W- -> ve e+ mu- vm~ ta- vt~
               else
                  write(6,*) 'Process not implemented,'
                  write(6,*)'check boson identifier'
               endif

           CASE DEFAULT
            if (with_anom) then ! anomalous gauge boson couplings
c           using global form factor for all tensors of one phase space point
c           this ensures proper cancellations for anomalous contributions
c           energy scale is invariant WWW mass
               if(bos.eq.343)then
                  call anomal_formfactor(qwpwmwp(0,J),qp1(0,J),q1m(0,J),qp2(0,J))
                  call WWW_atoww1_anomal(v(0,1,J),J,343,WWW_aww(0,J,1),NW)               !A -> ve  e+  mu-  vm~
                  call WWW_atoww2_anomal(v(0,1,J),J,343,WWW_aww(0,J,2),NW)               !A -> mu- vm~ vt   tau+
                  call WWW_ztoww1_anomal(v(0,1,J),J,343,WWW_zww(0,J,1),NW)               !Z -> ve  e+  mu-  vm~
                  call WWW_ztoww2_anomal(v(0,1,J),J,343,WWW_zww(0,J,2),NW)               !Z -> mu- vm~ vt   tau+
                  call WWW_wptowww_anomal(v(0,1,J),J,1,wptowww(0,J),NW)                  !W+ -> ve e+ mu- vm~ vt ta+
               else if(bos.eq.434) then
                  call anomal_formfactor(qwpwmwm(0,J),q1p(0,J),qm1(0,J),qm2(0,J))
                  call WWW_atoww1_anomal(v(0,1,J),J,434,WWW_aww(0,J,1),NW)               !A -> ve  e+  mu-  vm~
                  call WWW_atoww3_anomal(v(0,1,J),J,434,WWW_aww(0,J,3),NW)               !A -> ve  e+  tau- vt~
                  call WWW_ztoww1_anomal(v(0,1,J),J,434,WWW_zww(0,J,1),NW)               !Z -> ve  e+  mu-  vm~
                  call WWW_ztoww3_anomal(v(0,1,J),J,434,WWW_zww(0,J,3),NW)               !Z -> ve  e+  tau- vt~
                  call WWW_wmtowww_anomal(v(0,1,J),J,1,wmtowww(0,J),NW)                  !W- -> ve e+ mu- vm~ ta- vt~
               else
                  write(6,*) 'Process not implemented,'
                  write(6,*)'check boson identifier'
               endif
            else if (with_kk) then ! Kaluza Klein Scenario
#ifdef WITH_KK
               if(bos.eq.343)then
                  call WWW_atoww1(v(0,1,J),J,343,WWW_aww(0,J,1),NW)               !A -> W+ W- -> ve  e+  mu-  vm~
                  call WWW_atoww2(v(0,1,J),J,343,WWW_aww(0,J,2),NW)               !A -> W+ W- -> mu- vm~ vt   tau+
                  call WWW_ztoww1(v(0,1,J),J,343,WWW_zww(0,J,1),NW)               !Z -> W+ W- -> ve  e+  mu-  vm~
                  call WWW_ztoww2(v(0,1,J),J,343,WWW_zww(0,J,2),NW)               !Z -> W+ W- -> mu- vm~ vt   tau+
                  call WWW_wptowww_kk(v(0,1,J),J,1,wptowww(0,J),NW)               !W+ -> WWW -> ve e+ mu- vm~ vt ta+
               else if(bos.eq.434) then
                  call WWW_atoww1(v(0,1,J),J,434,WWW_aww(0,J,1),NW)               !A -> W+ W- -> ve  e+  mu-  vm~
                  call WWW_atoww3(v(0,1,J),J,434,WWW_aww(0,J,3),NW)               !A -> W+ W- -> ve  e+  tau- vt~
                  call WWW_ztoww1(v(0,1,J),J,434,WWW_zww(0,J,1),NW)               !Z -> W+ W- -> ve  e+  mu-  vm~
                  call WWW_ztoww3(v(0,1,J),J,434,WWW_zww(0,J,3),NW)               !Z -> W+ W- -> ve  e+  tau- vt~
                  call WWW_wmtowww_kk(v(0,1,J),J,1,wmtowww(0,J),NW)               !W- -> WWW -> ve e+ mu- vm~ ta- vt~
               else
                  write(6,*) 'Process not implemented,'
                  write(6,*)'check boson identifier'
               endif
#endif
            else ! SM
               if(bos.eq.343)then
                  call WWW_atoww1(v(0,1,J),J,343,WWW_aww(0,J,1),NW)               !A -> ve  e+  mu-  vm~
                  call WWW_atoww2(v(0,1,J),J,343,WWW_aww(0,J,2),NW)               !A -> mu- vm~ vt   tau+
                  call WWW_ztoww1(v(0,1,J),J,343,WWW_zww(0,J,1),NW)               !Z -> ve  e+  mu-  vm~
                  call WWW_ztoww2(v(0,1,J),J,343,WWW_zww(0,J,2),NW)               !Z -> mu- vm~ vt   tau+
                  call WWW_wptowww(v(0,1,J),J,1,wptowww(0,J),NW)                  !W+ -> ve e+ mu- vm~ vt ta+
               else if(bos.eq.434) then
                  call WWW_atoww1(v(0,1,J),J,434,WWW_aww(0,J,1),NW)               !A -> ve  e+  mu-  vm~
                  call WWW_atoww3(v(0,1,J),J,434,WWW_aww(0,J,3),NW)               !A -> ve  e+  tau- vt~
                  call WWW_ztoww1(v(0,1,J),J,434,WWW_zww(0,J,1),NW)               !Z -> ve  e+  mu-  vm~
                  call WWW_ztoww3(v(0,1,J),J,434,WWW_zww(0,J,3),NW)               !Z -> ve  e+  tau- vt~
                  call WWW_wmtowww(v(0,1,J),J,1,wmtowww(0,J),NW)                  !W- -> ve e+ mu- vm~ ta- vt~
               else
                  write(6,*) 'Process not implemented,'
                  write(6,*)'check boson identifier'
               endif
            endif
           END SELECT

      enddo        !j summation

c  ---------------------------------------------------------------------

c call PDF subroutine in order to determine parton
c distributions in the incoming (anti)protons.
c for the NLO contributions x1 = x*y in my notes with x=xuz, y=xi(i)

      x1 = xi(1)*xuz
      q_sf = sqrt(mufsq(1,1))
      call pdfproton( xi(1), q_sf, pdf(-6,1,1) )     ! f_a(y)=f_a(x1/z)

      if (nlo.gt.0) then
      q_sf = sqrt(mufsq(1,2))                        ! f_a(x1) for upper line 
      call pdfproton( x1, q_sf, pdf(-6,1,2) )        !   NLO correction

      if (mufsq(1,3).ne.mufsq(1,1)) then            ! f_a(x1) for lower line
         q_sf = sqrt(mufsq(1,3))                     !   NLO correction
         call pdfproton( xi(1), q_sf, pdf(-6,1,3) )
      else
         do i = -6,6
            pdf(i,1,3) = pdf(i,1,1)
         enddo
      endif
      endif

      x2 = xi(2)*xuz
      q_sf = sqrt(mufsq(2,1))
      call pdfproton( xi(2), q_sf, pdf(-6,2,1) )     ! f_b(y)=f_a(x2/z)

      if (nlo.gt.0) then
      q_sf = sqrt(mufsq(2,3))                        ! f_b(x2) for lower line 
      call pdfproton( x2, q_sf, pdf(-6,2,3) )        !   NLO correction

      if (mufsq(2,2).ne.mufsq(2,1)) then             ! f_b(x2) for upper line
         q_sf = sqrt(mufsq(2,2))                     !   NLO correction
         call pdfproton( xi(2), q_sf, pdf(-6,2,2) )
      else
         do i = -6,6
            pdf(i,2,2) = pdf(i,2,1)
         enddo
      endif
      endif


      if (nlo.gt.0) then
c and fill the coefficient functions for the finite subtraction terms
C Note that the color factors TR and C2 are NOT!!! included here

      if(bos.eq.343) then
         lnQomu(2) = log(qwwwp(4,2)/mufsq(1,2))
         lnQomu(3) = log(qwwwp(4,3)/mufsq(2,3))
c        factor 4*pi*alfas/(2x pa pb ) = 4*pi*alfas/Q^2
         tgs2oqsq(2) = 4d0*pi*als(1,2)/qwwwp(4,2)
         tgs2oqsq(3) = 4d0*pi*als(2,3)/qwwwp(4,3)
      else
         lnQomu(2) = log(qwwwm(4,2)/mufsq(1,2))
         lnQomu(3) = log(qwwwm(4,3)/mufsq(2,3))
c        factor 4*pi*alfas/(2x pa pb ) = 4*pi*alfas/Q^2
         tgs2oqsq(2) = 4d0*pi*als(1,2)/qwwwm(4,2)
         tgs2oqsq(3) = 4d0*pi*als(2,3)/qwwwm(4,3)
      endif

      omxi(2) = 1d0-x1
      omxi(3) = 1d0-x2

      z = xuz

      do j = 2,3

         ln1mxi = log(omxi(j))

C cut off Int_x1^1 dz log(1-z) at z<1-5E-6. This assures that the relative
C error of Int_x1^1 dz log(1-z) is less than 1E-4/(1-x1). see notes p 28.4

         if (z.lt.0.999995) then
            lnrat = lnQomu(j) + log((1d0-z)**2/z)
            lnz = log(z)
            Ax(j) = 2d0*((z**2+(1-z)**2)*lnrat + 2*z*(1-z))
            Bx(j) = 2d0*(2d0*(lnrat+lnz)/(1d0-z))
            Cx(j) = 2d0*(1-z - 2d0*lnz/(1d0-z) - (1+z)*lnrat)
            Dxi(j) = 2d0*(1.5d0*lnQomu(j) + 2d0*ln1mxi*lnQomu(j)
     1               + 2d0*ln1mxi**2) + crealVV
            do i = 1,5
               pdf(-i,j-1,j)=(pdf(-i,j-1,j)*(Dxi(j)/omxi(j)-Bx(j)) +
     1                     pdf(-i,j-1,1)*(Bx(j)+Cx(j)))*tgs2oqsq(j)
     2                     *1d0/(1d0-z)
               pdf(i,j-1,j)= (pdf(i,j-1,j)*(Dxi(j)/omxi(j)-Bx(j)) +
     1                     pdf(i,j-1,1)*(Bx(j)+Cx(j)))*tgs2oqsq(j)
     2                     *1d0/(1d0-z)

            enddo
            pdf(0,j-1,j)=pdf(0,j-1,1)*Ax(j)*tgs2oqsq(j)
     1                     *1d0/(1d0-z)
         endif
      enddo
      endif

C******************************************************************************

C     Now sum over the subprocesses contributing to W+ W- W production

      nmax = FL_VVVg(fsign,-1)          !reset counter for subprocesses to 0

C*******************  q1 qbar2 ---> W+ W- W+ g   ******************************
c   physToDiag(ext.momentum label) = Feynman diagram label

      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=2
      physToDiag(3)=3    ! gluon

      fsign(1) = 1          !incoming quark
      fsign(2) = -1         !incoming antiquark
      fsign(3) = 1          !outgoing ve
      fsign(4) = -1         !outgoing ep
      fsign(5) = 1          !outgoing mu-
      fsign(6) = -1         !outgoing vm~
      fsign(7) = 1          !outgoing ta-
      fsign(8) = -1         !outgoing ta+
      gsign    = 1          !outgoing gluon

      call qq_WWWj(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     &              pdf,res,nmin,nmax,NW)

C*******************  q2bar q1 ---> W+ W- W+ g   ******************************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=1
      physToDiag(3)=3    ! gluon

      fsign(1) = 1
      fsign(2) = -1
      gsign    = 1

      call qq_WWWj(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     &              pdf,res,nmin,nmax,NW)

C*******************  g q1 ---> W+ W- W+ q2  **********************************

      physToDiag(1)=3
      physToDiag(2)=1
      physToDiag(3)=2

      fsign(1) = 1
      fsign(2) = 1
      gsign    = -1

      call qq_WWWj(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     1              pdf,res,nmin,nmax,NW)

C*******************  qbar2 g ---> W+ W- W+ qbar1   ***************************

      physToDiag(1)=2
      physToDiag(2)=3
      physToDiag(3)=1

      fsign(1) = -1
      fsign(2) = -1
      gsign = -1

      call qq_WWWj(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     1              pdf,res,nmin,nmax,NW)

C*******************  g q2bar ---> W+ W- W+ q1bar  ****************************

      physToDiag(1)=3
      physToDiag(2)=2
      physToDiag(3)=1

      fsign(1) = -1
      fsign(2) = -1
      gsign = -1

      call qq_WWWj(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     1              pdf,res,nmin,nmax,NW)

c*******************  q1 g ---> W+ W- W+ q2   *********************************

      physToDiag(1)=1
      physToDiag(2)=3
      physToDiag(3)=2

      fsign(1) = 1
      fsign(2) = 1
      gsign    = -1

      call qq_WWWj(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     1              pdf,res,nmin,nmax,NW)

c*******************end of process evaluation**********************************

!     reverse momenta switch for W-W+W- (in order for a correct event output etc...)
      if(bos.eq.434) then
         do j=1,3
            do mu=0,3
            swmomenta(mu,1) = v(mu,3,j) ! ve
            swmomenta(mu,2) = v(mu,4,j) ! e+
            swmomenta(mu,3) = v(mu,1,j) ! mu-
            swmomenta(mu,4) = v(mu,2,j) ! vm~
            enddo
            do i=1,4
              do mu=0,3
                 v(mu,i,j) = swmomenta(mu,i)
              enddo
            enddo
         enddo
      endif

      SELECT CASE(procid)
      CASE(WPhadWMWP,WMhadWPWM)
       call hadDecayFactor_W(finalquarks(1), abs(mjj2(v(0,1,1),v(0,2,1))), N_gen_W, fac_W)
       do j=0,3
        DO IRES = 1,NMAX
          res(IRES,j) = res(IRES,j) * fac_W            ! factor fac_W for hadronic decay
        ENDDO
       enddo
      CASE(WPWMhadWP,WMWPhadWM)
       call hadDecayFactor_W(finalquarks(1), abs(mjj2(v(0,3,1),v(0,4,1))), N_gen_W, fac_W)
       do j=0,3
        DO IRES = 1,NMAX
          res(IRES,j) = res(IRES,j) * fac_W            ! factor fac_W for hadronic decay
        ENDDO
       enddo
      END SELECT

      if (init.eq.1) then
         init = init+1
      endif
      do j=0,3
         m2s(j) = 0
         DO IRES = 1,NMAX
            m2s(j) = m2s(j) + RES(IRES,j)
         ENDDO

         if(j.eq.0)then
c...Les Houches interface - the most propable subprocess
            if ((lha.or.hepmc).and..not.doNLO) then
               i=0
               weight=0.d0
               rnumb=RandomNumber()
               do while((i.le.nmax).and.(weight.le.rnumb*m2s(0)))
                  i=i+1
                  weight=weight+res(i,0)
                  iprocess=i
               enddo
               SELECT CASE(procid)
               CASE(WPhadWMWP, WMWPhadWM)
                  if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
                     rnumb=RandomNumber()
                     finalquarks_psp(1) =  2 + 2* INT(rnumb*2)
                     finalquarks_psp(2) = -1 - 2* INT(rnumb*2)
                  endif
               CASE(WMhadWPWM, WPWMhadWP)
                  if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
                     rnumb=RandomNumber()
                     finalquarks_psp(1) =  1 + 2* INT(rnumb*2)
                     finalquarks_psp(2) = -2 - 2* INT(rnumb*2)
                  endif
               END SELECT
            endif
         endif
      enddo

      RETURN
      END

c******************************************************************************
c
c    end subroutine m2s_qqWWWj
c
c******************************************************************************

c******************************************************************************
c
c   begin subroutine qq_WWWj
c
c******************************************************************************
c
c       Stefan Prestel, <sprestel@particle.uni-karlsruhe.de>
c	Initial version:  2008
c	Last modified by Vera Hankele:  2009 January 28
c
c       qq_WWWj calls the amplitude square routines 
c             qqWWWj    for qq --> WWW g --> g 6 leptons
c                           gq --> WWW q --> q 6 leptons 
c       and   qqWWW for the subtraction terms and finite collinear
c       remainders for the NLO cross section calculation
c
c  INPUT:  p(0:3,3,3)      external physical parton momenta
c          v(0:3,nv,3)       Z decay momenta, for Higgs production only
c                          the sum q(mu) = v(mu,1)+...+v(mu,nv) is needed
c          physToDiag(3)   physToDiag(ext.mom. label) = Feynman diagram label
c          fsign,gsign     sign factors for amplitude calls; see qqZqq(j)
c          nlo             nlo = 0: calculate LO only (i.e. no subtraction)
c                          nlo = 1: calculate full NLO subtraction
c          lok(3)          lok(ID)=T means momenta set ID passes acceptance
c                          cuts, i.e res(k,ID) needs to be calculated
c  OUTPUT:
c          udbarWWW(ID)   output in format as in qqZqq(j), but corresponding
c          etc.            to res(*,ID) in m2s_qqZqq
c                          ID = 1  : the real emission |M|^2 * pdf
c                             = 2,3: sutraction terms for emision off 12 or 
c     or   res(k,ID)                 34 line
c                             = 0  : subtracted result which drives 
c                                    integration, i.e 
c                                    res(*,0) = res(*,1)+res(*,2)+res(*,3)
c
c  In and Output
c    nmin, nmax            range of process numbers for this call
c
c
c*****************************************************************************

      SUBROUTINE qq_WWWj(xi,p,v,physToDiag,fsign,gsign,bos, 
     1                    nlo,lok,xuz,
     1                    pdf,res,nmin,nmax,NW)
c
      IMPLICIT NONE

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"

c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/coupl.inc"

      double precision p(0:3,max_p,max_kin), v(0:3,max_v,max_kin), xi(nx), xuz
      double precision pdf(-6:6,2,3), res(maxnumsubproc,0:3)
      integer physToDiag(3), fsign(3+n_v), 
     1        gsign, bos, nlo, nmin, nmax
      logical lok(3)
      integer NW

      double precision C2, TR, N                              ! color factors
      parameter(N=3d0,TR=0.5d0,c2=4d0/3d0)

c Note: Factor 3 from color sum included in qqWWW. I am using color summed
c amplitudes here, while CS use color averages. Thus I later divide by
c 8*3 for an initial qg state and by 3*3 for an initial qq state.

      double precision udbarWWWb(2,2:3), dubarWWWb(2,2:3)
      double precision udbarWWW, dubarWWW, udbarWWWg, dubarWWWg
      double precision udbarWWW2, dubarWWW2, udbarWWWg2, dubarWWWg2

      double precision pbar(0:3,3+n_v), qbar(0:4),q2,v2,sub(2:3),
     1       dotrr, gs2(2:3), polcolq, polcolg
      double precision NCmatrixelt(0:1,3,2), temp, temp2


      integer iflav(3), FL_VVVg
      external dotrr, FL_VVVg

      double complex zero
      parameter (zero=(0d0,0d0))

      save gs2,polcolq,polcolg,udbarWWWb,dubarWWWb
      integer i,if1,j,k,mu

      logical ldebug,ldebugs,ldebug_soco
      parameter (ldebug=.false.,ldebugs = .false.,ldebug_soco = .false.)

      integer stopp   ! variable to stop programm after debugpspoints times
      data stopp /1/
      save stopp
      integer counter ! to count number of M^2 differing to madgraph by epsilon
      data counter/0/
      save counter

      integer debugpspoints    ! number of phasespace points tested
      parameter (debugpspoints = 6)
      double precision epsilon  ! agreement between madgraph, handwritten M^2
      parameter ( epsilon = 1.0d-8)

c------------------------------------------------------------------------------

c initialize output

      udbarWWW  = 0.0d0
      dubarWWW  = 0.0d0
      udbarWWWg = 0.0d0
      dubarWWWg = 0.0d0

      do k =0,1
         do i =1,3
            do j =1,2
                NCmatrixelt(k,i,j) = 0.0d0
            enddo
         enddo
      enddo

c if this is a new phase space point -->
c reset the LO amplitude to 0 to avoid wrong subtraction
      if (physToDiag(1).eq.1.and.physToDiag(2).eq.2) then
         do i = 1,2
            do j = 2,3
               udbarWWWb(i,j) = 0d0
               dubarWWWb(i,j) = 0d0
            enddo
         enddo

c determine strong coupling gs for the two quark lines and factor
c  for polarization and spin average
         gs2(2) = 4d0*pi*als(1,1)
         gs2(3) = 4d0*pi*als(2,1)
         polcolq = 1d0/(4d0*N**2*xi(1)*xi(2))
         polcolg = 1d0/(4d0*N*(N**2-1)*xi(1)*xi(2))
      endif

      nmin = nmax+1
c get the real emission amplitude squared, store it in udbarWWWg, dubarWWWg

      do mu = 0,3
         do i = 1,3
            pbar(mu,physToDiag(i))=p(mu,i,1)
         enddo
         qbar(mu) = pbar(mu,3)
      enddo
      qbar(4)=0d0

      do mu = 0,3               ! kinematics for qq -> WWW -> 6 leptons
         pbar(mu,3) = v(mu,1,1) ! ve
         pbar(mu,4) = v(mu,2,1) ! e+
         pbar(mu,5) = v(mu,3,1) ! mu-
         pbar(mu,6) = v(mu,4,1) ! vmu~
         pbar(mu,7) = v(mu,5,1) ! vt  or ta-
         pbar(mu,8) = v(mu,6,1) ! ta+ or vt~
         pbar(mu,9) = qbar(mu)
      enddo
      fsign(9) = gsign

      if(bos.eq.343) then

         if (lok(1) .or. nlo.eq.0) then

#ifdef WITH_MADGRAPH
            if (ldebug) then
               call SUDBAR_WWWg(pbar,fsign,udbarWWWg2,NW)          !madgraph
            endif
#endif

            call qqWWWpj(pbar,fsign,qbar,gsign,udbarWWWg,NW)

#ifdef WITH_MADGRAPH
            if (ldebug) then
               if (abs(udbarWWWg/udbarWWWg2 - 1) .gt. epsilon) then
c                  print*," "
c                  print*,"udbarWWWg =  ",udbarWWWg
c                  print*,"udbarWWWg2 = ",udbarWWWg2
c                  print*,"udbarWWWg/udbarWWWg2 - 1 = ",udbarWWWg/udbarWWWg2-1
c                  print*," "
                  counter = counter +1
c                  print*,"----------------------------------------"
               endif
               if (stopp .eq. debugpspoints) then
                  print*,counter, " points in ",debugpspoints,
     1               " differ by number greater than ",epsilon
                  print*,"Leaving VBFNLO"
                  stop
               endif
               stopp = stopp +1
            endif
#endif

         endif

      else if(bos.eq.434) then

         if (lok(1) .or. nlo.eq.0) then

#ifdef WITH_MADGRAPH
            if (ldebug) then
               call SDUBAR_WWWg(pbar,fsign,dubarWWWg2,NW)
            endif
#endif

            call qqWWWmj(pbar,fsign,qbar,gsign,dubarWWWg,NW)

#ifdef WITH_MADGRAPH
            if (ldebug) then
               if (abs(dubarWWWg/dubarWWWg2 - 1) .gt. epsilon) then
c                  print*," "
c                  print*,"dubarWWWg =  ",dubarWWWg
c                  print*,"dubarWWWg2 = ",dubarWWWg2
c                  print*,"dubarWWWg/dubarWWWg2 - 1 = ",dubarWWWg/dubarWWWg2-1
c                  print*," "
                  counter = counter +1
c                  print*,"-------------------------------------------"
               endif
               if (stopp .eq. debugpspoints) then
                  print*,counter, " points in ",debugpspoints," 
     1              differ by number greater than ", epsilon
                  print*,"Leaving VBFNLO"
                  stop
               endif
               stopp = stopp +1
            endif
#endif

         endif

      else
         write(6,*) "Process not implemented, check boson identifier"
      endif

      NCmatrixelt(0,2,1)= udbarWWWg
      NCmatrixelt(1,2,1)= dubarWWWg

c for the NLO case get the subtraction terms;
c first the case with a final state gluon (gsign = +1)

c*********************************************
c BEGIN W+W-W+ PROCESS EVALUATION ************
c*********************************************
      if(bos.eq.343) then

c -------begin gsign=1 section--------------
        if (gsign.eq.1) then

#ifdef WITH_NLO
         if(nlo.eq.1) then
            do j=2,3 ! j=2: emission off 1,2 line ; j=3: emission off 3,4 line
               do mu =0,3
                  do i=1,2
                     pbar(mu,PhysToDiag(i)) = p(mu,i,j)
                  enddo
                  do i=1,6
                     pbar(mu,i+2) = v(mu,i,j)
                  enddo
               enddo

               if(lok(j)) then

#ifdef WITH_MADGRAPH
                  if(ldebug) then
                     call sudbar_WWW(pbar,fsign,udbarWWW2,NW) !Madgraph
                  endif
#endif

                  call QQWWWp(pbar,fsign,0,J,udbarWWW,NW)

#ifdef WITH_MADGRAPH
                  if(ldebug) then
                     if(abs(udbarWWW/udbarWWW2-1).gt.1.0d-2) then
                        print*,"(NLO =1) udbarWWW(j=",j,")  : ",udbarWWW
                        print*,"(NLO =1) udbarWWW2(j=",j,") : ",udbarWWW2
                        print*,"(NLO =1) udbarWWW/udbarWWW2(j=",j,")-1:",
     1                                    udbarWWW/udbarWWW2 -1
                     endif
                  endif
#endif

               else
                  udbarWWW = 0.0d0
               endif

               NCmatrixelt(0,j,2) = udbarWWW ! Born M^2 for mom. config. 2, 3
               udbarWWWb(PhysToDiag(1),j) = udbarWWW ! save M^2 for later use

               q2    = 2.0d0*xuz*dotrr(qbar,p(0,j-1,1))
               v2    = 2.0d0*gs2(j)*(2.0d0/(1-xuz) - (1+xuz))
               sub(j)= v2/q2 !M_B^2*sub are dipoles
            enddo
         endif !nlo.eq.1
#endif /* WITH_NLO */

         iflav(3) = 0           ! final state gluon id
         do if1=1,(nfl/2)*2,2
            iflav(1)=( if1 - PhysToDiag(1) + 2 )*fsign(PhysToDiag(1))
            iflav(2)=( if1 - PhysToDiag(2) + 2 )*fsign(PhysToDiag(2))

            if(ldebug) then
               print*,'iflav(1) :  ',iflav(1)
               print*,'iflav(2) :  ',iflav(2)
               print*,'iflav(3) :  ',iflav(3)
            endif

            k=FL_VVVg(iflav,1)

c            if(k.lt.15) print*,'(final gluon ) k = ',k

            if (lok(1)) then
            res(k,1)=pdf(sign1*iflav(1),1,1)
     &              *pdf(sign2*iflav(2),2,1)
     &              *NCmatrixelt(mod(if1+1,2),2,1)*polcolq
            else
               res(k,1) = 0.0d0
            endif

c---------------begin nlo=1 block----------------------------------

            if(nlo.eq.1.and.(lok(2).or.lok(3))) then
c emission from upper line
               res(k,2) =( pdf(sign1*iflav(1),1,2)*pdf(sign2*iflav(2),2,2)!coll
     &                   - pdf(sign1*iflav(1),1,1)*pdf(sign2*iflav(2),2,2)
     &                        *sub(2))  !real emission
     &                   *C2*NCmatrixelt(mod(if1+1,2),2,2)
c emission from lower line
               res(k,3) =( pdf(sign1*iflav(1),1,3)*pdf(sign2*iflav(2),2,3)!coll
     &                   - pdf(sign1*iflav(1),1,3)*pdf(sign2*iflav(2),2,1)
     &                        *sub(3))  !real emission
     &                   *C2*NCmatrixelt(mod(if1+1,2),3,2)

               res(k,2) = res(k,2)*polcolq
               res(k,3) = res(k,3)*polcolq
               res(k,0) = res(k,1) + res(k,2) + res(k,3)

c debugging for collinear and soft divergences -----------------------
                  if(ldebug_soco) then                !begin debug
                     temp = qbar(0)                        !soft
                     if(temp.le.1.0d-4) then
                        temp = abs(res(k,1)/(res(k,2)+res(k,3)))
                        if(temp.le.0.9d0.or.temp.ge.1.1d0) then
                           print*,' '
                           print*,'gluon in the final state '
                           print*,'if1 = ',if1
                           print*,'qbar(0) = ',qbar(0)
                           print*,'res(k,1) = ',res(k,1)
                           print*,'res(k,2)+res(k,3)=',res(k,2)+res(k,3)
                           print*,'res(k,0) = ',res(k,0)
                           print*,'|res(k,1)/res(k,2)| = ',temp
                           print*,' '
                        endif
                     endif

                     temp  = abs(dotrr(qbar,pbar(0,1)))    !collinear
                     temp2 = abs(dotrr(qbar,pbar(0,2)))
                     if(temp.le.1.0d-4.or.temp2.le.1.0d-4) then
                        temp = abs(res(k,1)/(res(k,2)+res(k,3)))
                           print*,' '
                           print*,'gluon in the final state '
                           print*,'if1 = ',if1
                           print*,'qbar*pbar(0,1)=',dotrr(qbar,pbar(0,1))
                           print*,'qbar*pbar(0,2)=',dotrr(qbar,pbar(0,2))
                           print*,'res(k,1) = ',res(k,1)
                           print*,'res(k,2)+res(k,3)=',res(k,2)+res(k,3)
                           print*,'res(k,0) = ',res(k,0)
                           print*,'|res(k,1)/res(k,2)| = ',temp
                           print*,' '
                     endif
                  endif ! end debug
c --------------------------------------------------------------------------

            else
               res(k,0) = res(k,1)
               res(k,2) = 0.0d0
               res(k,3) = 0.0d0

            endif
         enddo

c------------------end of gsign = 1 section-----------------------------

        elseif (gsign.eq.-1) then           !initial gluon section

         j = 0
         if (physToDiag(1).eq.3) then
            j=2                                  ! j=2: emission off 1,2 line
         elseif (physToDiag(2).eq.3) then
            j=3                                  ! j=3: emission off 3,4 line
         endif

         if(nlo.eq.1) then
            q2 = 2.0d0*dotrr(p(0,3,1),p(0,j-1,1))*xuz
            v2 = 2.0d0*gs2(j)*( (1-xuz)**2 + xuz**2 )
            sub(j) = v2/q2                              ! M_B^2*sub are dipoles

            if (physToDiag(1).eq.1.or.PhysToDiag(2).eq.2) then
               NCmatrixelt(0,2,2) = udbarWWWb(1,2)
               NCmatrixelt(0,3,2) = udbarWWWb(1,3)
            else
               NCmatrixelt(0,2,2) = udbarWWWb(2,2)
               NCmatrixelt(0,3,2) = udbarWWWb(2,3)
            endif
         endif

         do if1=1,(nfl/2)*2,2
            iflav(j-1) = 0
            iflav(4-j)=( if1 - PhysToDiag(4-j) + 2 )*fsign(PhysToDiag(4-j))
            iflav(3)  =( if1 - PhysToDiag(3) + 2 )*fsign(physToDiag(3))

            if(ldebug) then
               print*,'iflav(1) :  ',iflav(1)
               print*,'iflav(2) :  ',iflav(2)
               print*,'iflav(3) :  ',iflav(3)
            endif

            k=FL_VVVg(iflav,j)

c            if(k.lt.15) print*,'(initial gluon ) k = ',k

            if (lok(1)) then
               res(k,1)=pdf(sign1*iflav(1),1,1)  
     &              *pdf(sign2*iflav(2),2,1)      
     &              *NCmatrixelt(mod(if1+1,2),2,1)*polcolg
            else
               res(k,1) = 0.0d0
            endif

c ---------------begin nlo=1 block ---------------

            if(nlo.eq.1.and.(lok(2).or.lok(3))) then
               res(k,j) =( pdf(sign1*iflav(1),1,j)*pdf(sign2*iflav(2),2,j)!coll
     &                   - pdf(sign1*iflav(1),1,2*j-3)*
     &                    pdf(sign2*iflav(2),2,4-j)*sub(j)) !real emission
     &                   *C2*NCmatrixelt(mod(if1+1,2),j,2)

               res(k,j)   = res(k,j)*polcolg
               res(k,5-j) = 0.0d0
               res(k,0)   = res(k,1) + res(k,j)

c debugging for collinear and soft divergences -----------------------
                  if (ldebug_soco) then      !begin debug
                     temp = p(0,3,1)                       !soft
                     if(temp.le.1d-2) then
                         temp = abs(res(k,1)/res(k,j))
                        if(temp.le.0.99.or.temp.ge.1.01) then
                           print*,' '
                           print*,' gluon in the initial state'
                           print*,' k = ',k,' if1 = ',if1      
                           print*,'p(0,3,1) = ',p(0,3,1)
                           print*,'res(k,1) = ',res(k,1)
                           print*,'res(k,j) = ',res(k,j)
                           print*,'res(k,-j+5) = ',res(k,-j+5)
                           print*,'res(0)/res(j) = ',temp
                        endif
                     endif

                     temp = abs(dotrr(p(0,3,1),p(0,1,1)))  !collinear
                     temp2 = abs(dotrr(p(0,3,1),p(0,2,1)))
                     if(temp.le.0.0001.or.temp2.le.0.0001) then
                        if(j.eq.3) then
                           temp = abs((res(k,1))+res(k-4,1)/
     &                                 (res(k,2)+res(k-4,2)))
                           if(temp.le.0.99.or.temp.ge.1.01) then
                              print*,' '
                              print*,' gluon in the initial state'
                              print*,' k = ',k,' if1 = ',if1      
                              print*,'pa*p1= ',dotrr(p(0,3,1),p(0,1,1))
                              print*,'pb*p1= ',dotrr(p(0,3,1),p(0,2,1))
                              print*,'res(k,1) = ',res(k,1)
                              print*,'res(k,j) = ',res(k,j)
                              print*,'res(k,-j+5) = ',res(k,-j+5)
                              print*,'res(k,1)+res(k-4,1) = ',
     &                                    res(k,1)+res(k-4,1)
                              print*,'res(k,2)+res(k-4,2) = ',
     &                                    res(k,2)+res(k-4,2)
                              print*,'res(0)/res(j) = ',temp
                           endif
                        endif
                     endif
                  endif !end debug
c----------------------------------------------------------------------
            else
               res(k,0) = res(k,1)
               res(k,2) = 0.0d0
               res(k,3) = 0.0d0

            endif
         enddo
        endif !---------end of gsign=-1 section--------------------------

c*****************************************
c END OF W+W-W+ PROCESS EVALUATION *******
c*****************************************

      elseif(bos.eq.434) then

c----------begin gsign=1 section----------
        if (gsign.eq.1) then

#ifdef WITH_NLO
         if(nlo.eq.1) then
            do j=2,3
               do mu =0,3
                  do i=1,2
                     pbar(mu,PhysToDiag(i)) = p(mu,i,j)
                  enddo
                  do i=1,6
                     pbar(mu,i+2) = v(mu,i,j)
                  enddo
               enddo

               if(lok(j)) then

#ifdef WITH_MADGRAPH
                  if(ldebug) then
                     call sdubar_WWW(pbar,fsign,dubarWWW2,NW) ! Madgraph code
                  endif
#endif

                  call QQWWWm(pbar,fsign,0,J,dubarWWW,NW)

#ifdef WITH_MADGRAPH
                  if(ldebug) then
                     if(abs(dubarWWW/dubarWWW2 -1).gt.1.0d-2) then
                        print*,"(NLO =1) dubarWWW(j=",j,")  : ",dubarWWW
                        print*,"(NLO =1) dubarWWW2(j=",j,") : ",dubarWWW2
                        print*,"(NLO =1) dubarWWW/dubarWWW2(j=",j,") -1",
     &                              dubarWWW/dubarWWW2 -1
                     endif
                  endif
#endif

               else
                  dubarWWW = 0.0d0
               endif

               NCmatrixelt(1,j,2)=dubarWWW  ! Born M^2 for mom. conf. 2 and 3
               dubarWWWb(PhysToDiag(1),j)=dubarWWW  ! save M^2 for later use

               q2    = 2.0d0*xuz*dotrr(qbar,p(0,j-1,1))
               v2    = 2.0d0*gs2(j)*(2.0d0/(1-xuz) - (1+xuz))
               sub(j)= v2/q2 !M_B^2*sub are dipoles
            enddo
         endif !nlo.eq.1
#endif /* WITH_NLO */

         iflav(3) = 0           ! final state gluon id
         do if1=1,(nfl/2)*2,2
            iflav(1)=( if1 + PhysToDiag(1) - 1 )*fsign(PhysToDiag(1))
            iflav(2)=( if1 + PhysToDiag(2) - 1 )*fsign(PhysToDiag(2))

            k=FL_VVVg(iflav,1)
            if(lok(1)) then
               res(k,1) = pdf(sign1*iflav(1),1,1)
     &                   *pdf(sign2*iflav(2),2,1)
     &                   *NCmatrixelt(mod(if1,2),2,1)*polcolq
            else
               res(k,1) = 0.0d0
            endif

c-------------------------begin nlo=1 block-------------------------------

            if(nlo.eq.1.and.(lok(2).or.lok(3))) then
#ifdef WITH_NLO
c  emission from upper line
               res(k,2) =( pdf(sign1*iflav(1),1,2)*pdf(sign2*iflav(2),2,2)!coll
     &                   - pdf(sign1*iflav(1),1,1)*pdf(sign2*iflav(2),2,2)
     &                      *sub(2))  !real emission
     &                   *C2*NCmatrixelt(mod(if1,2),2,2)
c  emission from lower line
               res(k,3) =( pdf(sign1*iflav(1),1,3)*pdf(sign2*iflav(2),2,3)!coll
     &                   - pdf(sign1*iflav(1),1,3)*pdf(sign2*iflav(2),2,1)
     &                      *sub(3))  ! real emission
     &                   *C2*NCmatrixelt(mod(if1,2),3,2)

               res(k,2) = res(k,2)*polcolq
               res(k,3) = res(k,3)*polcolq
               res(k,0) = res(k,1) + res(k,2) + res(k,3)


c debugging for collinear and soft divergences -----------------------
                  if(ldebug_soco) then                    !begin debug
                     temp = qbar(0)                        !soft
                     if(temp.le.1.0d-4) then
                        temp = abs(res(k,1)/(res(k,2)+res(k,3)))
                        if(temp.le.0.9d0.or.temp.ge.1.1d0) then
                           print*,' '
                           print*,'gluon in the final state '
                           print*,'if1 = ',if1
                           print*,'qbar(0) = ',qbar(0)
                           print*,'res(k,1) = ',res(k,1)
                           print*,'res(k,2)+res(k,3)=',res(k,2)+res(k,3)
                           print*,'res(k,0) = ',res(k,0)
                           print*,'|res(k,1)/res(k,2)| = ',temp
                           print*,' '
                        endif
                     endif

                     temp  = abs(dotrr(qbar,pbar(0,1)))    !collinear
                     temp2 = abs(dotrr(qbar,pbar(0,2)))
                     if(temp.le.1.0d-4.or.temp2.le.1.0d-4) then
                        temp = abs(res(k,1)/(res(k,2)+res(k,3)))
                           print*,' '
                           print*,'gluon in the final state '
                           print*,'if1 = ',if1
                           print*,'qbar*pbar(0,1)=',dotrr(qbar,pbar(0,1))
                           print*,'qbar*pbar(0,2)=',dotrr(qbar,pbar(0,2))
                           print*,'res(k,1) = ',res(k,1)
                           print*,'res(k,2)+res(k,3)=',res(k,2)+res(k,3)
                           print*,'res(k,0) = ',res(k,0)
                           print*,'|res(k,1)/res(k,2)| = ',temp
                           print*,' '
                     endif
                  endif ! end debug
c----------------------------------------------------------------------
#endif /* WITH_NLO */
            else
               res(k,0) = res(k,1)
               res(k,2) = 0.0d0
               res(k,3) = 0.0d0

            endif
         enddo
c-------------------------end of gsign=1 section-------------------

        elseif (gsign.eq.-1) then           !initial gluon section
c         if (ldebug) then
c            print*,' initial state gluon section in qq_WWWj '
c         endif

         j = 0
         if (physToDiag(1).eq.3) then
            j=2                                  ! j=2: emission off 1,2 line
         elseif (physToDiag(2).eq.3) then
            j=3                                  ! j=3: emission off 3,4 line
         endif

#ifdef WITH_NLO
         if(nlo.eq.1 ) then
            q2 = 2.0d0*dotrr(p(0,3,1),p(0,j-1,1))*xuz
            v2 = 2.0d0*gs2(j)*( (1-xuz)**2 + xuz**2 )
            sub(j) = v2/q2                              ! M_B^2*sub are dipoles

            if (physToDiag(1).eq.1.or.PhysToDiag(2).eq.2) then
               NCmatrixelt(1,2,2) = dubarWWWb(1,2)
               NCmatrixelt(1,3,2) = dubarWWWb(1,3)
            else
               NCmatrixelt(1,2,2) = dubarWWWb(2,2)
               NCmatrixelt(1,3,2) = dubarWWWb(2,3)
            endif
         endif
#endif

         do if1=1,(nfl/2)*2,2
            iflav(j-1) = 0
            iflav(4-j) = ( if1 + PhysToDiag(4-j) - 1 )*fsign(PhysToDiag(4-j))
            iflav(3) = ( if1 + PhysToDiag(3) - 1 )*fsign(PhysToDiag(3))

            k=FL_VVVg(iflav,j)
            if (lok(1)) then
               res(k,1)=pdf(sign1*iflav(1),1,1)  !1 for j=2;3 for j=3
     &              *pdf(sign2*iflav(2),2,1)       !2 for j=2;1 for j=3
*               res(k,1)=pdf(sign1*iflav(1),1,2*j-3)  !1 for j=2;3 for j=3
*     &              *pdf(sign2*iflav(2),2,4-j)       !2 for j=2;1 for j=3
     &              *NCmatrixelt(mod(if1,2),2,1)*polcolg
            else
               res(k,1) = 0.0d0
            endif

c-------------------begin nlo=1 block---------------------------

            if(nlo.eq.1.and.(lok(2).or.lok(3))) then
#ifdef WITH_NLO
               res(k,j) =( pdf(sign1*iflav(1),1,j)*pdf(sign2*iflav(2),2,j)!coll
     &                   -pdf(sign1*iflav(1),1,2*j-3)*pdf(sign2*iflav(2),2,4-j)
     &                     *sub(j)) !real emission
     &                   *C2*NCmatrixelt(mod(if1,2),j,2)

               res(k,j)   = res(k,j)*polcolg
               res(k,5-j) = 0.0d0
               res(k,0)   = res(k,1) + res(k,j)

c debugging for collinear and soft divergences -----------------------
                  if (ldebug_soco) then             !begin debug
                     temp = p(0,3,1)                       !soft
                     if(temp.le.1d-2) then
                         temp = abs(res(k,1)/res(k,j))
                        if(temp.le.0.99.or.temp.ge.1.01) then
                           print*,' '
                           print*,' gluon in the initial state'
                           print*,' k = ',k,' if1 = ',if1      
                           print*,'p(0,3,1) = ',p(0,3,1)
                           print*,'res(k,1) = ',res(k,1)
                           print*,'res(k,j) = ',res(k,j)
                           print*,'res(k,-j+5) = ',res(k,-j+5)
                           print*,'res(k,1) = ',res(k,1)
                           print*,'res(k,j) = ',res(k,j)
                           print*,'res(0)/res(j) = ',temp
                           endif
                        endif

                     temp = abs(dotrr(p(0,3,1),p(0,1,1)))  !collinear
                     temp2 = abs(dotrr(p(0,3,1),p(0,2,1)))
                     if(temp.le.0.0001.or.temp2.le.0.0001) then
                        if(j.eq.3) then
                           temp = abs((res(k,1))+res(k-4,1)/
     &                                 (res(k,2)+res(k-4,2)))
                           if(temp.le.0.99.or.temp.ge.1.01) then
                              print*,' '
                              print*,' gluon in the initial state'
                              print*,' k = ',k,' if1 = ',if1      
                              print*,'pa*p1= ',dotrr(p(0,3,1),p(0,1,1))
                              print*,'pb*p1= ',dotrr(p(0,3,1),p(0,2,1))
                              print*,'res(k,1) = ',res(k,1)
                              print*,'res(k,j) = ',res(k,j)
                              print*,'res(k,-j+5) = ',res(k,-j+5)
                              print*,'res(k,1)+res(k-4,1) = ',
     &                                   res(k,1)+res(k-4,1)
                              print*,'res(k,2)+res(k-4,2) = ',
     &                                   res(k,2)+res(k-4,2)
                              print*,'res(k,1) = ',res(k,1)
                              print*,'res(k,j) = ',res(k,j)
                              print*,'res(0)/res(j) = ',temp
                           endif
                        endif
                     endif
                  endif !end debug
c----------------------------------------------------------------------
#endif /* WITH_NLO */
            else
               res(k,0) = res(k,1)
               res(k,2) = 0.0d0
               res(k,3) = 0.0d0
            endif !nlo.and.lok
         enddo
        endif  !---------end of gsign=-1 section ----------

c******************* END OF COMPUTATION **********************

      else
         write(6,*) 'Process not implemented, check boson identifier'
      endif

      nmax = FL_VVVg(iflav,-2)

      end

c******************************************************************
c
c   end subroutine qq_WWWj
c
c*****************************************************************
