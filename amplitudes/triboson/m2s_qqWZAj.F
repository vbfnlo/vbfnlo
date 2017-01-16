c*********************  subroutine m2s_qqWZAj*********************************
c
c INPUT     bos       Boson identifier, 1 = A, 3,4=W
c           nlo       NLO=1: create subtraction term; LO = 0
c           lok       need to calculate m2s(1:3)? T or F 
c           xi        Feynman x parameters of incoming fermions
c           p         momenta associated with external fermions
c           v         momenta associated with WZA decay fermions
c           rn        random number array
c           xuz       x_ija for subtraction terms
c          
c
c OUTPUT    m2s       |M|^2*pdf1*pdf2 
c
c       Adapted from m2s_qqWWAj from Guiseppe Bozzi
c       Adapted from m2s_qqVVVj from Vera Hankele
c
c       michael rauch, <rauch@particle.physik.uni-karlsruhe.de>
c       Initial version: December 2008
c       Last modified: January 2009
c
c***************************************************************************

      SUBROUTINE m2s_qqWZAj(bos,nlo,lok,xi,p,v,rn,xuz,m2s,ps_number)

      IMPLICIT NONE

c declare input/output variables

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/process.inc"

c     hadronic decay
      integer N_gen_up, N_gen_down
      real*8 fac_Z_up, fac_Z_down, mjj2
      integer N_gen_W
      real*8 fac_W
      external mjj2

      INTEGER bos, ps_number
      DOUBLE PRECISION xi(nx), p(0:3,max_p,max_kin), v(0:3,max_v,max_kin), 
     1       rn(1), xuz, m2s(0:3)
      logical lok(3)

c  helicity selection

      INTEGER h
      COMMON /hcount / h

      COMPLEX*16 zero
      PARAMETER (zero=(0d0,0d0))

c declare external functions

      integer FL_VVVg
      external FL_VVVg

c alfas, scales etc

      DOUBLE PRECISION x1,x2
      DOUBLE PRECISION lnQomu(2:3), omxi(2:3), 
     &       Ax(2:3), Bx(2:3), Cx(2:3), Dxi(2:3), tgs2oqsq(2:3), 
     &       ln1mxi, lnrat, z, lnz

c declare local variables

      DOUBLE PRECISION q_sf
      DOUBLE PRECISION dotrr
      EXTERNAL dotrr

      INTEGER init/0/, I, J, mu
      SAVE init

c declare variables for summation process

      INTEGER ires, nmaxold, FSIGN(3+n_v), gsign
      SAVE nmaxold

      INTEGER physToDiag(3), nlo, nmin, nmax
      SAVE nmin, nmax

c store contributions from subprocess ip in res(ip,ID) where
c  ID = 1  : the real emission |M|^2 * pdf
c     = 2,3: subtraction terms for emision off upper or lower line
c     = 0  : subtracted result which drives integration, i.e
c res(*,0) = res(*,1)+res(*,2)+res(*,3)

      DOUBLE PRECISION res(maxnumsubproc,0:3)

c declare parton distribution variables

      DOUBLE PRECISION pdf(-6:6,2,3)

c define program switches

      LOGICAL ldebug
      PARAMETER (ldebug=.false.)

c variables for leshouches interface
      DOUBLE PRECISION weight,rnumb,RandomNumber

c---------------------------- BEGIN CODE --------------------------------

      if ( init .ne. 0) then
         do i = 1,nmax
            do j=0,3
               res(i,j) = 0
            enddo
         enddo
      else if (init.eq.0) then
         if (bos.eq.321) then
            write(6,*) " "
            write(6,*) "W+ Z A j amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
         else if (bos.eq.421) then
            write(6,*) " "
            write(6,*) "W- Z A j amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
         endif
         call vtowza_reset
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

      do J = 1,3
         do mu = 0,3
            qp(mu,J) = v(mu,1,J)+v(mu,2,J)   !W+
            qz(mu,J) = v(mu,3,J)+v(mu,4,J)   !Z 
            qa(mu,J) = v(mu,5,J)             !A
            qwz(mu,J) = qp(mu,J) + qz(mu,J)
            qwza(mu,J) = qwz(mu,J) + qa(mu,J) 
         enddo

         qp(4,J) = dotrr(qp(0,J),qp(0,J))
         qz(4,J) = dotrr(qz(0,J),qz(0,J))
         qa(4,J) = dotrr(qa(0,J),qa(0,J))
         qwz(4,J) = dotrr(qwz(0,J),qwz(0,J))
         qwza(4,J) = dotrr(qwza(0,J),qwza(0,J))
      enddo

c select helicity: h = random number for photon helicity

      h = mod(h,4) + 1
      iu = (-1)**(h+1)
      it = (-1)**((h-1)/2)

      if ((lha.or.hepmc) .and. (.not. doNLO)) then
         helicity(1)=-1
         helicity(2)= 1
         helicity(3)=-iu
         helicity(4)= iu
         helicity(5)= it
      end if

c lepton spinors and W+,W-,A polarisation vectors
      do J = 1,3

         CALL OXXXXX(v(0,1,J),ZERO ,-1,1,wve(1,J))   !ve !e-
         CALL IXXXXX(v(0,2,J),ZERO ,1,-1,wep(1,J))   !e+ !ve~
         CALL OXXXXX(v(0,3,J),ZERO ,-iu,1,wmu(1,J))   !mu-      
         CALL IXXXXX(v(0,4,J),ZERO ,iu,-1,wmp(1,J))   !mu+
         CALL VXXXXX(qa(0,J),ZERO ,it,1,wa(1,J))   !A
         CALL JIOXXX(wep(1,J),wve(1,J),GWF,WMASS,WWIDTH, wp(1,J))   !W+


           SELECT CASE(procid)
           CASE(WPZhadA, WMZhadA)
             ! for general output to all combinations: up-type first
             if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
               ! set couplings for Z hadronic decay into uubar
               call setZtouu
             endif
             CALL JIOXXX(wmp(1,J),wmu(1,J),GZ_ZF,ZMASS,ZWIDTH, ztau(1,J)) !Z
             CALL JIOXXX(wmp(1,J),wmu(1,J),GZ_AF,ZERO ,ZERO  , atau(1,J)) !Z

           CASE DEFAULT
             CALL JIOXXX(wmp(1,J),wmu(1,J),GZL,ZMASS,ZWIDTH, ztau(1,J)) !Z
             CALL JIOXXX(wmp(1,J),wmu(1,J),GAL,ZERO ,ZERO  , atau(1,J)) !Z

           END SELECT


         SELECT CASE(procid)
         CASE(WPhadZA, WPZhadA, WMhadZA, WMZhadA)

            call anomal_formfactor(qwza(0,J),qp(0,J),qz(0,J),qa(0,J))
            if (bos.eq.321) then
               call TRIBOSwptowz_had(v(0,1,J),J,h,wpwz(0,J),0)!W -> W+ Z -> ve e+ mu- mu+
               call TRIBOSwptowa_had(v(0,1,J),J,h,wpwa(0,J),0)!W -> W+ A -> ve e+ a
               call wptowza_had(v(0,1,J),J,h,zwwa(0,J),0)    !W -> WZA -> ve e+ mu- mu+ a
            else if (bos.eq.421) then
               call TRIBOSwmtowz_had(v(0,1,J),J,h,wpwz(0,J),0)!W -> W- Z -> ve~ e- mu- mu+
               call TRIBOSwmtowa_had(v(0,1,J),J,h,wpwa(0,J),0)!W -> W- A -> ve~ e- a
               call wmtowza_had(v(0,1,J),J,h,zwwa(0,J),0)    !W -> WZA -> ve~ e- mu- mu+ a
            endif
            call TRIBOSatoza_had(v(0,1,J),J,h,aza(0,J),0)  !A -> Z A -> mu- mu+ a
            call TRIBOSztoza_had(v(0,1,J),J,h,zaz(0,J),0)  !Z -> Z A -> mu- mu+ a


         CASE DEFAULT

         if (with_anom) then ! anomalous gauge boson couplings
c           using global form factor for all tensors of one phase space point
c           this ensures proper cancellations for anomalous contributions
c           energy scale is invariant ZZW mass
            call anomal_formfactor(qwza(0,J),qp(0,J),qz(0,J),qa(0,J))

c for WWA->4l+A precalculate leptonic tensors A->WWA,Z->WWA etc..
c leptonic tensors, last entry defines whether the non resonant diagrams
c contribute ( = 0: all diagrams, = 1: no non resonant diagrams)

            if (bos.eq.321) then
               call TRIBOSwptowz_anomal(v(0,1,J),J,h,wpwz(0,J),0)!W -> W+ Z -> ve e+ mu- mu+
               call TRIBOSwptowa_anomal(v(0,1,J),J,h,wpwa(0,J),0)!W -> W+ A -> ve e+ a
               call wptowza_anomal(v(0,1,J),J,h,zwwa(0,J),0)    !W -> WZA -> ve e+ mu- mu+ a
            else if (bos.eq.421) then
               call TRIBOSwmtowz_anomal(v(0,1,J),J,h,wpwz(0,J),0)!W -> W- Z -> ve~ e- mu- mu+
               call TRIBOSwmtowa_anomal(v(0,1,J),J,h,wpwa(0,J),0)!W -> W- A -> ve~ e- a
               call wmtowza_anomal(v(0,1,J),J,h,zwwa(0,J),0)    !W -> WZA -> ve~ e- mu- mu+ a
            endif
            call TRIBOSatoza_anomal(v(0,1,J),J,h,aza(0,J),0)  !A -> Z A -> mu- mu+ a
            call TRIBOSztoza_anomal(v(0,1,J),J,h,zaz(0,J),0)  !Z -> Z A -> mu- mu+ a


         else !SM

c for WWA->4l+A precalculate leptonic tensors A->WWA,Z->WWA etc..
c leptonic tensors, last entry defines whether the non resonant diagrams
c contribute ( = 0: all diagrams, = 1: no non resonant diagrams)

            if (bos.eq.321) then
               call TRIBOSwptowz(v(0,1,J),J,h,wpwz(0,J),0)!W -> W+ Z -> ve e+ mu- mu+
               call TRIBOSwptowa(v(0,1,J),J,h,wpwa(0,J),0)!W -> W+ A -> ve e+ a
               call wptowza(v(0,1,J),J,h,zwwa(0,J),0)    !W -> WZA -> ve e+ mu- mu+ a
            else if (bos.eq.421) then
               call TRIBOSwmtowz(v(0,1,J),J,h,wpwz(0,J),0)!W -> W- Z -> ve~ e- mu- mu+
               call TRIBOSwmtowa(v(0,1,J),J,h,wpwa(0,J),0)!W -> W- A -> ve~ e- a
               call wmtowza(v(0,1,J),J,h,zwwa(0,J),0)    !W -> WZA -> ve~ e- mu- mu+ a
            endif
            call TRIBOSatoza(v(0,1,J),J,h,aza(0,J),0)  !A -> Z A -> mu- mu+ a
            call TRIBOSztoza(v(0,1,J),J,h,zaz(0,J),0)  !Z -> Z A -> mu- mu+ a

         endif !anom or SM

         END SELECT

      enddo

c  ---------------------------------------------------------------------

c scales and als:

      if (.false.) then
         write(6,*) " "
         write(6,*) "xi(1), xi(2) =", xi
         write(6,*) "mu_f L=1:",sqrt(mufsq(1,1)),sqrt(mufsq(2,1))
         write(6,*) "mu_f L=2:",sqrt(mufsq(1,2)),sqrt(mufsq(2,2))
         write(6,*) "mu_f L=3:",sqrt(mufsq(1,3)),sqrt(mufsq(2,3))
         write(6,*) " "
         write(6,*) "xi(1), xi(2) =", xi
         write(6,*) "alphas L=1:",als(1,1), als(2,1)
         write(6,*) "alphas L=2:",als(1,2), als(2,2)
         write(6,*) "alphas L=3:",als(1,3), als(2,3)
      end if

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

      if (mufsq(1,3).ne.mufsq(1,1)) then             ! f_a(x1) for lower line
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

      lnQomu(2) = log(qwza(4,2)/mufsq(1,2))
      lnQomu(3) = log(qwza(4,3)/mufsq(2,3))
 
      omxi(2) = 1d0-x1
      omxi(3) = 1d0-x2

c factor 4*pi*alfas/(2x pa pb ) = 4*pi*alfas/Q^2

      tgs2oqsq(2) = 4d0*pi*als(1,2)/qwza(4,2)
      tgs2oqsq(3) = 4d0*pi*als(2,3)/qwza(4,3)

      z = xuz

      do j = 2,3

         ln1mxi = log(omxi(j))

C cut off Int_x1^1 dz log(1-z) at z<1-5E-6. This assures that the relative
C error of Int_x1^1 dz log(1-z) is less than 1E-4/(1-x1). see notes p 28.4

         if (z.lt.0.999995) then
            lnrat = lnQomu(j) + log((1d0-z)**2/z)
            lnz = log(z)
            Ax(j) = 2d0*((z**2+(1-z)**2)*lnrat+2*z*(1-z))
            Bx(j) = 2d0*(2d0*(lnrat+lnz)/(1d0-z))
            Cx(j) = 2d0*(1-z-2d0*lnz/(1d0-z)-(1+z)*lnrat)
            Dxi(j) = 2d0*(1.5d0*lnQomu(j)+2d0*ln1mxi*lnQomu(j)+
     1                    2d0*ln1mxi**2)+crealVV
            do i = 1,5
               pdf(-i,j-1,j)=(pdf(-i,j-1,j)*(Dxi(j)/omxi(j)-Bx(j)) +
     1                        pdf(-i,j-1,1)*(Bx(j)+Cx(j)))
     2                        *tgs2oqsq(j)*1d0/(1d0-z)
               pdf(i,j-1,j)= (pdf(i,j-1,j)*(Dxi(j)/omxi(j)-Bx(j)) +
     1                        pdf(i,j-1,1)*(Bx(j)+Cx(j)))
     2                        *tgs2oqsq(j)*1d0/(1d0-z)
            enddo
            pdf(0,j-1,j)=pdf(0,j-1,1)*Ax(j)*tgs2oqsq(j)*1d0/(1d0-z)
         endif
      enddo
      endif

C********************************************************************

C     Now sum over the subprocesses contributing to W+ Z A production

      nmax = FL_VVVg(fsign,-1)!reset counter for subprocesses to 0

C*******************  q1 qbar2 ---> W+ Z A g   **********************
         
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
      fsign(7) = 1          !outgoing A
      gsign    = 1          !outgoing gluon

      call qq_WZAj(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     &             pdf,res,nmin,nmax,ps_number)

C*******************  qbar1 q2 ---> W+ Z A g   **********************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=1
      physToDiag(3)=3    ! gluon

      fsign(1) = 1
      fsign(2) = -1
      gsign    = 1

      call qq_WZAj(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     &             pdf,res,nmin,nmax,ps_number)

C*******************  g q1 ---> W+ Z A q2  ***********************************

      physToDiag(1)=3
      physToDiag(2)=1
      physToDiag(3)=2

      fsign(1) = 1
      fsign(2) = 1
      gsign    = -1

      call qq_WZAj(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     1             pdf,res,nmin,nmax,ps_number)

C*******************  qbar2 g ---> W+ Z A qbar1   **********************

      physToDiag(1)=2
      physToDiag(2)=3
      physToDiag(3)=1

      fsign(1) = -1
      fsign(2) = -1
      gsign    = -1

      call qq_WZAj(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     1             pdf,res,nmin,nmax,ps_number)

C*******************  g q2bar ---> W+ Z A q1bar  **********************

      physToDiag(1)=3
      physToDiag(2)=2
      physToDiag(3)=1

      fsign(1) = -1
      fsign(2) = -1
      gsign    = -1

      call qq_WZAj(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     1             pdf,res,nmin,nmax,ps_number)

C*******************  q1 g ---> W+ Z A q2   **********************

      physToDiag(1)=1
      physToDiag(2)=3
      physToDiag(3)=2

      fsign(1) = 1
      fsign(2) = 1
      gsign    = -1

      call qq_WZAj(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     1             pdf,res,nmin,nmax,ps_number)


c**************  end of process evaluation part one ******************


! now run again for down-type decay products if 93 93 / 94 94 is 
! requested and the z decays hadronically.

      SELECT CASE(procid)
      CASE(WPZhadA, WMZhadA)
       if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
        ! set couplings for Z hadronic decay into ddbar
        call setZtodd

        do J = 1,3
            CALL JIOXXX(wmp(1,J),wmu(1,J),GZ_ZF,ZMASS,ZWIDTH, ztau(1,J)) !Z
            CALL JIOXXX(wmp(1,J),wmu(1,J),GZ_AF,ZERO ,ZERO  , atau(1,J)) !Z
            call anomal_formfactor(qwza(0,J),qp(0,J),qz(0,J),qa(0,J))
            if (bos.eq.321) then
               call TRIBOSwptowz_had(v(0,1,J),J,h,wpwz(0,J),0)!W -> W+ Z -> ve e+ mu- mu+
               call TRIBOSwptowa_had(v(0,1,J),J,h,wpwa(0,J),0)!W -> W+ A -> ve e+ a
               call wptowza_had(v(0,1,J),J,h,zwwa(0,J),0)    !W -> WZA -> ve e+ mu- mu+ a
            else if (bos.eq.421) then
               call TRIBOSwmtowz_had(v(0,1,J),J,h,wpwz(0,J),0)!W -> W- Z -> ve~ e- mu- mu+
               call TRIBOSwmtowa_had(v(0,1,J),J,h,wpwa(0,J),0)!W -> W- A -> ve~ e- a
               call wmtowza_had(v(0,1,J),J,h,zwwa(0,J),0)    !W -> WZA -> ve~ e- mu- mu+ a
            endif
            call TRIBOSatoza_had(v(0,1,J),J,h,aza(0,J),0)  !A -> Z A -> mu- mu+ a
            call TRIBOSztoza_had(v(0,1,J),J,h,zaz(0,J),0)  !Z -> Z A -> mu- mu+ a
        enddo


C*******************  q1 qbar2 ---> W+ Z A g   **********************
         
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
      fsign(7) = 1          !outgoing A
      gsign    = 1          !outgoing gluon

      call qq_WZAj(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     &             pdf,res,nmin,nmax,ps_number)

C*******************  qbar1 q2 ---> W+ Z A g   **********************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=1
      physToDiag(3)=3    ! gluon

      fsign(1) = 1
      fsign(2) = -1
      gsign    = 1

      call qq_WZAj(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     &             pdf,res,nmin,nmax,ps_number)

C*******************  g q1 ---> W+ Z A q2  ***********************************

      physToDiag(1)=3
      physToDiag(2)=1
      physToDiag(3)=2

      fsign(1) = 1
      fsign(2) = 1
      gsign    = -1

      call qq_WZAj(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     1             pdf,res,nmin,nmax,ps_number)

C*******************  qbar2 g ---> W+ Z A qbar1   **********************

      physToDiag(1)=2
      physToDiag(2)=3
      physToDiag(3)=1

      fsign(1) = -1
      fsign(2) = -1
      gsign    = -1

      call qq_WZAj(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     1             pdf,res,nmin,nmax,ps_number)

C*******************  g q2bar ---> W+ Z A q1bar  **********************

      physToDiag(1)=3
      physToDiag(2)=2
      physToDiag(3)=1

      fsign(1) = -1
      fsign(2) = -1
      gsign    = -1

      call qq_WZAj(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     1             pdf,res,nmin,nmax,ps_number)

C*******************  q1 g ---> W+ Z A q2   **********************

      physToDiag(1)=1
      physToDiag(2)=3
      physToDiag(3)=2

      fsign(1) = 1
      fsign(2) = 1
      gsign    = -1

      call qq_WZAj(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     1             pdf,res,nmin,nmax,ps_number)


c***  end of extra process evaluation for hadronic decay  ***

       endif
      end select


c*****************  end of process evaluation  *********************

      if (init.eq.1) then
         init = init+1
      endif

      SELECT CASE(procid)
      CASE(WPhadZA, WMhadZA)
       call hadDecayFactor_W(finalquarks(1), abs(mjj2(v(0,1,1),v(0,2,1))), N_gen_W, fac_W)
       do j=0,3
        DO IRES = 1,NMAX
          res(IRES,j) = res(IRES,j) * fac_W            ! factor fac_W for hadronic decay
        ENDDO
       enddo
      CASE(WPZhadA, WMZhadA)
       call hadDecayFactor_Z(finalquarks(1), abs(mjj2(v(0,3,1),v(0,4,1))), N_gen_up, N_gen_down, fac_Z_up, fac_Z_down)
       do j=0,3
        DO IRES = 1,NMAX
          if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
             if (ires.le.(NMAX/2)) then
                res(IRES,j) = res(IRES,j) * fac_Z_up      ! factor fac_Z_up for up-type with all combinations
             elseif (ires.gt.(NMAX/2)) then
                res(IRES,j) = res(IRES,j) * fac_Z_down    ! factor fac_Z_down for down-type with all combinations
             endif
             if (mod(nmax,2).ne.0) then
               print*, "Something went wrong with the hadronic decays!"
               stop
             endif
          elseif (mod(abs(finalquarks(1)),2).eq.0) then
             res(IRES,j) = res(IRES,j) * fac_Z_up         ! up-type final states
          else
             res(IRES,j) = res(IRES,j) * fac_Z_down       ! down-type final states
          endif
        ENDDO
       enddo
      END SELECT

      do j=0,3
         m2s(j) = 0
         DO IRES = 1,NMAX
            m2s(j) = m2s(j) + RES(IRES,j)
         ENDDO

         if(j.eq.0)then
c...Les Houches interface - the most propable subprocess 3jets at LO  
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
               CASE(WPhadZA)
                  if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
                     rnumb=RandomNumber()
                     finalquarks_psp(1) =  2 + 2* INT(rnumb*2)
                     finalquarks_psp(2) = -1 - 2* INT(rnumb*2)
                  endif
               CASE(WMhadZA)
                  if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
                     rnumb=RandomNumber()
                     finalquarks_psp(1) =  1 + 2* INT(rnumb*2)
                     finalquarks_psp(2) = -2 - 2* INT(rnumb*2)
                  endif
               CASE(WPZhadA, WMZhadA)
                  if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
                     rnumb=RandomNumber()
                     if (i.le.(nmax/2)) then    ! up-type
                       finalquarks_psp(1) =  2 + 2* INT(rnumb*N_gen_up)
                       finalquarks_psp(2) = -2 - 2* INT(rnumb*N_gen_up)
                     else                       ! down-type
                       finalquarks_psp(1) =  1 + 2* INT(rnumb*N_gen_down)
                       finalquarks_psp(2) = -1 - 2* INT(rnumb*N_gen_down)
                     endif
                  endif
               END SELECT
            endif
         endif

         m2s(j) = m2s(j) * 4d0    ! factor 4 for random helicity summation
      enddo


c debugging ---------------------------------------------------------------

      if (ldebug) then
         do j=1,3
            if (lok(j) .and. m2s(j).eq.0 ) then
               print*,j," lok(j) = ",lok(j)," m2s = ",m2s(j)
            endif
         enddo
         if (m2s(0).ne.0) then
            if (abs((m2s(1)+m2s(2)+m2s(3))/m2s(0)-1).gt.1d-5) then
               print*," m2s(0) = ",m2s(0),(m2s(1)+m2s(2)+m2s(3))/m2s(0)
               print*," ratio subtr/real = ",(m2s(2)+m2s(3))/m2s(1)
               print*," m2s(0,1,2,3) = ",m2s
               print*
            endif
         else
            print*," m2s = ",m2s
         endif
         if (nmax.ne.nmaxold) print*," nmax: ",nmaxold,nmax
         nmaxold = nmax
      endif
      RETURN
      END

c*****************************************************************************
c
c    end subroutine m2s_qqWZA1j
c
c*****************************************************************************
c******************************************************************
c
c   begin subroutine qq_WZA1j
c
c*****************************************************************

      SUBROUTINE qq_WZAj(xi,p,v,physToDiag,fsign,gsign,bos, 
     1                   nlo,lok,xuz,
     1                   pdf,res,nmin,nmax,ps_number)
      IMPLICIT NONE

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"

      DOUBLE PRECISION p(0:3,max_p,max_kin), v(0:3,max_v,max_kin), xi(nx), xuz
      DOUBLE PRECISION pdf(-6:6,2,3), res(maxnumsubproc,0:3)
      INTEGER physToDiag(3), fsign(3+n_v), 
     1        gsign, bos, nlo, nmin, nmax, ps_number
      LOGICAL lok(3)
c
c qq_WZAj calls the amplitude square routines 
c             qqWZAj    for qq --> WZA g --> g 4 leptons
c                           gq --> WZA q --> q 4 leptons 
c
c for the subtraction terms for the NLO cross section calculation
c
c  INPUT:  p(0:3,3,3)      external physical parton momenta
c          v(0:3,nv,3)     lepton and photon momenta
c          physToDiag(3)   physToDiag(ext.mom. label) = Feynman diagram label
c          fsign,gsign     sign factors for amplitude calls; see qqZqq(j)
c          nlo             nlo = 0: calculate LO only (i.e. no subtraction)
c                          nlo = 1: calculate full NLO subtraction
c          lok(3)          lok(ID)=T means momenta set ID passes acceptance
c                          cuts, i.e res(k,ID) needs to be calculated
c  OUTPUT:
c                          ID = 1  : the real emission |M|^2 * pdf
c                             = 2,3: sutraction terms for emision off 
c          res(k,ID)                 upper or lower line
c                             = 0  : subtracted result which drives 
c                                    integration, i.e. 
c                                    res(*,0) = res(*,1)+res(*,2)+res(*,3)
c
c  In and Output
c          nmin, nmax      range of process numbers for this call
c
      double precision C2, TR, N                              ! color factors
      parameter(N=3d0,TR=0.5d0,c2=4d0/3d0)

c Note: Factor 3 from color sum included in qqWZA. I am using color summed
c amplitudes here, while CS use color averages. Thus I later divide by
c 8*3 for an initial qg state and by 3*3 for an initial qq state.

c alfas, scales etc

#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/coupl.inc"

c
c  helicity selection

      INTEGER h
      COMMON /hcount / h

      DOUBLE PRECISION udwzab(2,2:3), udwza
      DOUBLE PRECISION udwza2, udwzag
      DOUBLE PRECISION pbar(0:3,3+n_v), qbar(0:4),q2,v2,sub(2:3), 
     1       dotrr, gs2(2:3), polcolq, polcolg
      DOUBLE PRECISION NCmatrixelt(0:1,3,2), temp, temp1, temp2, temp3
      INTEGER iflav(3), diagToPhys(3), FL_VVVg
      EXTERNAL dotrr, FL_VVVg

      COMPLEX*16 zero
      PARAMETER (zero=(0d0,0d0))

      LOGICAL ldebug, lres, ldebugsoft
      SAVE gs2,polcolq,polcolg,udwzab
      INTEGER i,if1,j,k,mu
      
      PARAMETER (ldebug=.false.)
      PARAMETER (ldebugsoft=.false.)
      PARAMETER (lres=.false.)

c--------------------------------------------------------------------

      udwza  = 0d0
      udwzag = 0d0

      if (lres) then
        open(unit=31,file="irfinite.chk",ACCESS='APPEND')
      endif

c  if this is a new phase space point -->
c  reset the LO amplitude to 0 to avoid wrong subtraction
c      if (physToDiag(1).eq.1.and.physToDiag(2).eq.2) then
      if (Fl_VVVg(iflav,-2).eq.0) then
         do i = 1,2
            do j = 2,3
               udwzab(i,j) = 0d0
            enddo
         enddo

c  determine strong coupling gs for the two quark lines and factor for 
c  polarization and spin average
         gs2(2) = 4d0*pi*als(1,1)
         gs2(3) = 4d0*pi*als(2,1)
         polcolq = 1d0/(4d0*N**2*xi(1)*xi(2))
         polcolg = 1d0/(4d0*N*(N**2-1)*xi(1)*xi(2))
      endif

      nmin = nmax+1
      do i = 1,3
         diagToPhys(physToDiag(i)) = i
      enddo
         
c get the real emission amplitude squared, store it in udwza
      do mu = 0,3
         do i = 1,3
            pbar(mu,physToDiag(i))=p(mu,i,1)
         enddo
         qbar(mu) = pbar(mu,3)
      enddo
      qbar(4)=0d0

      do mu = 0,3               ! kinematics for qq -> WZA -> 4 leptons + A
         pbar(mu,3) = v(mu,1,1) ! ve
         pbar(mu,4) = v(mu,2,1) ! e+
         pbar(mu,5) = v(mu,3,1) ! mu-
         pbar(mu,6) = v(mu,4,1) ! mu+
         pbar(mu,7) = v(mu,5,1) ! A
         pbar(mu,8) = qbar(mu)
      enddo
      fsign(8) = gsign

c calculate the matrix elements pa + pb -> ve e+ mu- mu+ a p1
      if (lok(1) .or. nlo.eq.0) then

#ifdef WITH_MADGRAPH
         if (ldebug) then
            call sud_wzag(pbar,fsign,udwzag) !madgraph
         endif
#endif
         
         call qqWZAj(pbar,fsign,qbar,gsign,udwza,bos,ps_number)
         
#ifdef WITH_MADGRAPH
         if (ldebug) then
            if ( (abs(udwza/udwzag - 1) .gt. 1d-4) ) then
               print*,"p(0,1,1) = ",p(0,1,1)
               print*,"udwza  =   ",udwza
               print*,"udwzag =   ",udwzag
               print*,"udwza/udwzag - 1 = ",udwza/udwzag - 1
               print*," "
            endif
            read(*,*)
         endif
#endif

      endif

      NCmatrixelt(0,2,1)=udwza

c for the NLO case get the subtraction terms; 
c first the case with a final state gluon (gsign = +1)

      if (gsign.eq.1) then

#ifdef WITH_NLO
         if (nlo.eq.1) then
            do j = 2,3               ! j=2: emission off 1,2 line
               do mu = 0,3           ! j=3: emission off 3,4 line
                  do i = 1,2
                     pbar(mu,physToDiag(i))=p(mu,i,j)
                  enddo
                  do i = 1,5
                     pbar(mu,i+2)=v(mu,i,j)
                  enddo
               enddo

c calculation of matrix elements for q qb -> W Z A g:

               if(lok(j)) then
#ifdef WITH_MADGRAPH
                  if(ldebug) then
                     call sud_wza(pbar,fsign,udwza2)  !madgraph
                  endif
#endif

                  call QQWZA(pbar,fsign,0,J,udwza,ps_number)
  
#ifdef WITH_MADGRAPH
                  if (ldebug) then
                        print*,"udwza/udwza2 =   ",udwza/udwza2
                  endif
#endif
               else
                  udwza = 0d0
               endif

               NCmatrixelt(0,j,2)=udwza   !Born for momentum
                                          !configuration 2 and 3

c save matrix elements for later use with initial gluons
               udwzab(physToDiag(1),j) = udwza

               q2 = 2d0*xuz*dotrr(qbar,p(0,j-1,1)) !p(mu,j-1,1) is inc.quark
               v2 = 2d0*gs2(j)*                   ! 8 pi alfas
     &             ( 2d0/(1-xuz)-(1+xuz) )
               sub(j) = v2/q2

            enddo
         endif !nlo.eq.1
#endif /* WITH_NLO */

         iflav(3) = 0           ! final state gluon id
         do if1=2,(nfl/2)*2,2
            if (bos.eq.321) then
               iflav(1)=( if1 - physToDiag(1) + 1 )*fsign(physToDiag(1))
               iflav(2)=( if1 - physToDiag(2) + 1 )*fsign(physToDiag(2))
            else if (bos.eq.421) then
               iflav(1)=( if1 + physToDiag(1) - 2 )*fsign(physToDiag(1))
               iflav(2)=( if1 + physToDiag(2) - 2 )*fsign(physToDiag(2))
            endif
            k=FL_VVVg(iflav,2)
            if(lok(1)) then
               res(k,1)=pdf(sign1*iflav(1),1,1)  
     &           *pdf(sign2*iflav(2),2,1)  
     &           *NCmatrixelt(mod(if1,2),2,1)*polcolq
            else
               res(k,1) = 0
            endif

            if ( nlo.eq.1 .and. (lok(2).or.lok(3)) ) then
#ifdef WITH_NLO
               res(k,2) = ( pdf(sign1*iflav(1),1,2) ! collinear
     &              *pdf(sign2*iflav(2),2,2)
     &              -pdf(sign1*iflav(1),1,1) ! real emission
     &              *pdf(sign2*iflav(2),2,2 )*sub(2))*C2
     &              *NCmatrixelt(mod(if1,2),2,2)
               res(k,3) = ( pdf(sign1*iflav(1),1,3) ! collinear
     &              *pdf(sign2*iflav(2),2,3)
     &              -pdf(sign1*iflav(1),1,3) ! real emission
     &              *pdf(sign2*iflav(2),2,1)*sub(3))*C2
     &              *NCmatrixelt(mod(if1,2),3,2)
               
               res(k,2) = res(k,2) * polcolq
               res(k,3) = res(k,3) * polcolq
               res(k,0) = res(k,1) + res(k,2) + res(k,3)

c debugging for collinear and soft divergences -----------------------
               if (lres) then
                 temp1 = sqrt(abs(dotrr(qbar,pbar(0,1)))) ! collinear 1
                 temp2 = sqrt(abs(dotrr(qbar,pbar(0,2)))) ! collinear 2
                 temp3 = abs(res(k,1)/(res(k,2)+res(k,3)))
                 if ( (temp3.gt.0.) .and. (temp3.le.5.) )
     &             write(31,*) "0", 
     &               qbar(0),
     &               temp1,
     &               temp2,
     &               min(temp1,temp2),
     &               abs(res(k,1)/(res(k,2)+res(k,3)))
               endif

               if(ldebugsoft) then !begin debug
                  temp = qbar(0)                     !soft
                  temp1 = abs(dotrr(qbar,pbar(0,1))) !collinear
                  temp2 = abs(dotrr(qbar,pbar(0,2))) !collinear
                  temp3 = abs(res(k,1)/(res(k,2)+res(k,3)))
                  if (abs(1-temp3).ge.1d-3) then
                     if(temp.le.1d-3) then
                        print*," "
                        print*,"soft gluon in the final state"
                        print*,"if1 = ",if1
                        print*,"qbar(0) = ",qbar(0)
                        print*,"test=", sqrt(2d0*dotrr(pbar(0,1),pbar(0,2)))
                        print*,"qbar*pbar(0,1)=",dotrr(qbar,pbar(0,1))
                        print*,"qbar*pbar(0,2)=",dotrr(qbar,pbar(0,2))
                        print*,"res(k,1) = ",res(k,1)
                        print*,"res(k,2) = ",res(k,2)
                        print*,"res(k,3) = ",res(k,3)
                        print*,"res(k,2)+res(k,3)=",res(k,2)+res(k,3)
                        print*,"res(k,0) = ",res(k,0)
                        print*,"|res(k,1)/(res(k,2)+res(k,3))| = ",temp3
                        print*," "
                     endif
                     if(temp1.le.1d-1.or.temp2.le.1d-1) then
                        print*," "
                        print*,"collinear gluon in the final state"
                        print*,"if1 = ",if1
                        print*,"qbar(0) = ",qbar(0)
                        print*,"test=", sqrt(2d0*dotrr(pbar(0,1),pbar(0,2)))
                        print*,"qbar*pbar(0,1)=",dotrr(qbar,pbar(0,1))
                        print*,"qbar*pbar(0,2)=",dotrr(qbar,pbar(0,2))
                        print*,"res(k,1) = ",res(k,1)
                        print*,"res(k,2) = ",res(k,2)
                        print*,"res(k,3) = ",res(k,3)
                        print*,"res(k,2)+res(k,3)=",res(k,2)+res(k,3)
                        print*,"res(k,0) = ",res(k,0)
                        print*,"|res(k,1)/(res(k,2)+res(k,3))| = ",temp3
                        print*," "
                     endif
                  endif
               endif         !end debug
c----------------------------------------------------------------------
#endif /* WITH_NLO */
            else
               res(k,0) = res(k,1)
               res(k,2) = 0
               res(k,3) = 0
            endif !nlo

         enddo !if1
         
      elseif (gsign.eq.-1) then !initial gluon section
         j = 0
         if (nlo.eq.0) then
            if (physToDiag(1).eq.3) then
               j=2           ! j=2: emission off 1,2 line
            elseif (physToDiag(2).eq.3) then
               j=3           ! j=3: emission off 3,4 line
            endif
         elseif (nlo.eq.1) then
#ifdef WITH_NLO
            if (physToDiag(1).eq.3) then
               j=2           ! j=2: emission off 1,2 line
            elseif (physToDiag(2).eq.3) then
               j=3           ! j=3: emission off 3,4 line
            endif
            
            q2 = 2d0 * dotrr(p(0,3,1),p(0,j-1,1))*xuz
            v2 = 2d0*gs2(j) * ( (1-xuz)**2 + xuz**2 )
            sub(j) = v2/q2
            
            if (physToDiag(1).eq.1.or.PhysToDiag(2).eq.2) then
               NCmatrixelt(0,2,2)=udwzab(1,2)
               NCmatrixelt(0,3,2)=udwzab(1,3)
            else
               NCmatrixelt(0,2,2)=udwzab(2,2)
               NCmatrixelt(0,3,2)=udwzab(2,3)
            endif
#endif /* WITH_NLO */
         endif
         
         do if1=2,(nfl/2)*2,2
            iflav(j-1) = 0
            if (bos.eq.321) then
               iflav(4-j)=( if1 - physToDiag(4-j) + 1 )*fsign(physToDiag(4-j))
               iflav(3)  =( if1 - physToDiag(3) + 1 )*fsign(physToDiag(3))
            else if (bos.eq.421) then
               iflav(4-j)=( if1 + physToDiag(4-j) - 2 )*fsign(physToDiag(4-j))
               iflav(3)  =( if1 + physToDiag(3) - 2 )*fsign(physToDiag(3))
            endif
            k=FL_VVVg(iflav,j)
            if (lok(1)) then
               res(k,1)=pdf(sign1*iflav(1),1,1) 
     &              *pdf(sign2*iflav(2),2,1) 
     &              *NCmatrixelt(mod(if1,2),2,1)*polcolg
            else
               res(k,1) = 0
            endif
            
c--------------------------------------------------------------
            if ( nlo.eq.1 .and. (lok(2).or.lok(3)) ) then
#ifdef WITH_NLO                  
               res(k,j) = ( pdf(sign1*iflav(1),1,j) ! collinear
     &              *pdf(sign2*iflav(2),2,j)
     &              -pdf(sign1*iflav(1),1,2*j-3) ! real emission
     &              *pdf(sign2*iflav(2),2,4-j )*sub(j))
     &              *NCmatrixelt(mod(if1,2),j,2)
               
c--------------------------------------------------------------------
               
               res(k,j) = res(k,j)*C2*polcolg
               res(k,-j+5) = 0d0
               
               res(k,0) = res(k,1) + res(k,j)
               
               if (lres) then
                 temp1 = sqrt(abs(dotrr(p(0,j-1,1),p(0,3,1)))) ! collinear 1
                 temp2 = sqrt(abs(dotrr(p(0,j-1,1),p(0,4-j,1)))) ! collinear 2
                 temp3 = abs(res(k,1)/res(k,j))
                 if ( (temp3.gt.0.) .and. (temp3.le.5.) )
     &             write(31,*) "1", 
     &               p(0,j-1,1),
     &               temp1,
     &               temp2,
     &               min(temp1,temp2),
     &               temp3
               endif
               if (ldebugsoft) then !begin debug
                  temp1 = abs(dotrr(p(0,3,1),p(0,1,1)))
                  temp2 = abs(dotrr(p(0,3,1),p(0,2,1)))
                  temp3 = abs(res(k,1)/res(k,j))
                  if (abs(1-temp3).ge.1d-3) then
                    if(temp1.le.1d-1.or.temp2.le.1d-1) then
                       print*," "
                       print*,"collinear quark (initial state gluon)"
                       print*,"k = ",k," if1 = ",if1      
                       print*,"p(0,3,1) =",p(0,3,1)
                       print*,"p(1,3,1)/p(0,3,1)",p(1,3,1)/p(0,3,1)
                       print*,"p(2,3,1)/p(0,3,1)",p(2,3,1)/p(0,3,1)
                       print*,"p(3,3,1)/p(0,3,1)",p(3,3,1)/p(0,3,1)
                       print*,"pa*p1= ",dotrr(p(0,3,1),p(0,1,1))
                       print*,"pb*p1= ",dotrr(p(0,3,1),p(0,2,1))
                       print*,"res(k,1) = ",res(k,1)
                       print*,"res(k,j) = ",res(k,j)
                       print*,"res(1)/res(j) = ",temp3
                    endif
                  endif
               endif         !end debug
#endif /* WITH_NLO */
            else
               res(k,0) = res(k,1)
               res(k,2) = 0
               res(k,3) = 0
            endif
            
         enddo
      endif
      
      nmax = FL_VVVg(iflav,-2)

      if (lres) then
        close(31)
      endif
         
      end
      
c******************************************************************
c     
c   end subroutine qq_WZA1j
c
c*****************************************************************

