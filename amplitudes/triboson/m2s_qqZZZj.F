c*********************  subroutine m2s_qqZZZj*********************************
c
c INPUT     bos       Boson identifier, 2 = Z, 3,4=W
c           nlo       NLO=1: create subtraction term; LO = 0
c           lok       need to calculate m2s(1:3)? T or F 
c           xi        Feynman x parameters of incoming fermions
c           p         momenta associated with external fermions
c           v         momenta associated with WWZ decay fermions
c           rn        random number array
c           xuz       x_ija for subtraction terms
c           NW        = 1: NW Approximation, = 0: Full calculation
c          
c
c OUTPUT    m2s       |M|^2*pdf1*pdf2 
c
c       Adapted from m2s_wbfzhg by D. Zeppenfeld
c
c       Vera Hankele, <vera@particle.uni-karlsruhe.de>
c       Initial version: January 2007
c       Last modified: June 2009
c
c***************************************************************************

      SUBROUTINE m2s_qqZZZj(bos,nlo,lok,xi,p,v,rn,xuz,m2s,NW,ps_number)

      IMPLICIT NONE

c declare input/output variables

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/process.inc"

c     for hadronic decay
      integer N_gen_up, N_gen_down
      real*8 fac_Z_up, fac_Z_down, mjj2
      external mjj2

      INTEGER bos, NW,ps_number
      DOUBLE PRECISION xi(nx), p(0:3,max_p,max_kin), v(0:3,max_v,max_kin), 
     1       rn(1), xuz, m2s(0:3)
      LOGICAL lok(3)

c  helicity selection

      INTEGER h
      COMMON /hcount / h

      DOUBLE COMPLEX zero
      PARAMETER (zero=(0d0,0d0))

c declare external functions

      INTEGER FL_VVVg
      EXTERNAL FL_VVVg

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

      INTEGER ires, FSIGN(3+n_v), gsign

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

      DOUBLE PRECISION weight,rnumb,RandomNumber

c---------------------------- BEGIN CODE --------------------------------

      if ( init .ne. 0) then
         do i = 1,nmax
            do j=0,3
               res(i,j) = 0
            enddo
         enddo
      else if (init.eq.0) then
         if (bos.eq.22 .or. bos.eq.22) then
            write(6,*) " "
            write(6,*) "Z Z Z j amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
            call vtozzz_reset
         endif
         call printnfl(.false.)
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
 
c for ZZZ->6l precalculate A->ZZZ,Z->ZZZ etc..

c lepton spinors and Z/A polarisation vectors
c lepton helicities for Z -> l- l+ not fixed
c -> sum over all possible helicities in |M|**2
 
c select helicity: h ... random number for lepton helicity
c                        combination (h=1:2)

      h = mod(h,8) + 1
      ie = sign(1,4-h)
      iu = (-1)**(h+1)
      it = (-1)**((h-1)/2)

c...Les Houches interface
        if ((lha.or.hepmc).and..not.doNLO) then
           helicity(1)=-ie
           helicity(2)= ie
           helicity(3)=-iu
           helicity(4)= iu
           helicity(5)=-it
           helicity(6)= it
        endif

      do J = 1,3

         CALL OXXXXX(v(0,1,J),ZERO ,-ie,+1,zem(1,J)) !W(1,3))       !e-
         CALL IXXXXX(v(0,2,J),ZERO ,+ie,-1,zep(1,J)) !W(1,4))       !e+
         CALL OXXXXX(v(0,3,J),ZERO ,-iu,+1,zmm(1,J)) !W(1,5))       !mu-
         CALL IXXXXX(v(0,4,J),ZERO ,+iu,-1,zmp(1,J)) !W(1,6))       !mu+
         CALL OXXXXX(v(0,5,J),ZERO ,-it,+1,ztm(1,J)) !W(1,7))       !ta-
         CALL IXXXXX(v(0,6,J),ZERO ,+it,-1,ztp(1,J)) !W(1,8))       !ta+

         SELECT CASE(procid)
         CASE(ZZZhad)
            ! for general output to all combinations: up-type first
            if ((finalquarks(1).eq.93 .or. finalquarks(1).eq.94)) then
              ! set couplings for Z hadronic decay into uubar
              call setZtouu
            endif
           CALL JIOXXX(zep(1,J),zem(1,J),GZ_AF,ZERO , ZERO ,ae(1,J))
           CALL JIOXXX(zep(1,J),zem(1,J),GZ_ZF,ZMASS,ZWIDTH,ze(1,J))

         CASE DEFAULT
           CALL JIOXXX(zep(1,J),zem(1,J),GAL,ZERO , ZERO ,ae(1,J))
           CALL JIOXXX(zep(1,J),zem(1,J),GZL,ZMASS,ZWIDTH,ze(1,J))

         END SELECT

         CALL JIOXXX(zmp(1,J),zmm(1,J),GAL,ZERO ,ZERO  ,amu(1,J))
         CALL JIOXXX(zmp(1,J),zmm(1,J),GZL,ZMASS,ZWIDTH,zmu(1,J))
         CALL JIOXXX(ztp(1,J),ztm(1,J),GAL,ZERO ,ZERO  ,atau(1,J))
         CALL JIOXXX(ztp(1,J),ztm(1,J),GZL,ZMASS,ZWIDTH,ztau(1,J))

         do mu = 0,3
            qze(mu,J) = v(mu,1,J)+v(mu,2,J)   !Z
            qzmu(mu,J) = v(mu,3,J)+v(mu,4,J)   !Z
            qztau(mu,J) = v(mu,5,J)+v(mu,6,J)   !Z
            qzzz(mu,J) = v(mu,1,J)+v(mu,2,J)+v(mu,3,J)+v(mu,4,J)
     1                  +v(mu,5,J)+v(mu,6,J)!+p(mu,3,J)
         enddo

         qze(4,J) = dotrr(qze(0,J),qze(0,J))
         qzmu(4,J) = dotrr(qzmu(0,J),qzmu(0,J))
         qztau(4,J) = dotrr(qztau(0,J),qztau(0,J))
         qzzz(4,J) = dotrr(qzzz(0,J),qzzz(0,J))


         SELECT CASE(procid)
         CASE(ZZZhad)
            call anomal_formfactor(qzzz(0,J),qze(0,J),qzmu(0,J),qztau(0,J))
c leptonic tensors
c all diagrams
            if (ps_number.eq.1) then
               call atozzz_had(v(0,1,J),J,h,azzz(0,J),NW) !A -> Z Z Z -> 6l
               call ztozzz_had(v(0,1,J),J,h,zzzz(0,J),NW) !Z -> Z Z Z -> 6l
               call zzz_atozz_had(v(0,1,J),v(0,2,J),v(0,3,J),v(0,4,J),-ie,-iu,azezm(0,J))
               call zzz_atozz_had(v(0,1,J),v(0,2,J),v(0,5,J),v(0,6,J),-ie,-it,azezt(0,J))
               call atozz_anomal(v(0,3,J),v(0,4,J),v(0,5,J),v(0,6,J),-iu,-it,azmzt(0,J))
               call zzz_ztozz_had(v(0,1,J),v(0,2,J),v(0,3,J),v(0,4,J),-ie,-iu,zzezm(0,J))
               call zzz_ztozz_had(v(0,1,J),v(0,2,J),v(0,5,J),v(0,6,J),-ie,-it,zzezt(0,J))
               call ztozz_anomal(v(0,3,J),v(0,4,J),v(0,5,J),v(0,6,J),-iu,-it,zzmzt(0,J))
            endif
c only higgs contributions
            call ztozzzHiggs_had(v(0,1,J),J,h,zzzz_higgs1(0,J),zzzz_higgs2(0,J), 
     1           zzzz_higgs3(0,J),NW) !Z -> H Z, H -> Z Z for all permutations

         CASE DEFAULT

      if (with_anom) then ! anomalous gauge boson couplings
c       using global form factor for all tensors of one phase space point
c       this ensures proper cancellations for anomalous contributions
c       energy scale is invariant ZZZ mass
         call anomal_formfactor(qzzz(0,J),qze(0,J),qzmu(0,J),qztau(0,J))


c leptonic tensors
         if (NW.eq.1) then
c NW approximation (only Higgs contribution exists in the leptonic tensors)
            call ztozzz_anomal(v(0,1,J),J,h,zzzz(0,J),NW) !Z -> Z Z Z -> 6l
            call atozzz_anomal(v(0,1,J),J,h,azzz(0,J),NW) !A -> Z Z Z -> 6l
c only higgs contributions
         call ztozzzHiggs_anomal(v(0,1,J),J,h,zzzz_higgs1(0,J),zzzz_higgs2(0,J), 
     1           zzzz_higgs3(0,J),NW) !Z -> H Z, H -> Z Z for all permutations
         else
c all diagrams
            if (ps_number.eq.1) then
               call atozzz_anomal(v(0,1,J),J,h,azzz(0,J),NW) !A -> Z Z Z -> 6l
               call ztozzz_anomal(v(0,1,J),J,h,zzzz(0,J),NW) !Z -> Z Z Z -> 6l
               call atozz_anomal(v(0,1,J),v(0,2,J),v(0,3,J),v(0,4,J),-ie,-iu,azezm(0,J))
               call atozz_anomal(v(0,1,J),v(0,2,J),v(0,5,J),v(0,6,J),-ie,-it,azezt(0,J))
               call atozz_anomal(v(0,3,J),v(0,4,J),v(0,5,J),v(0,6,J),-iu,-it,azmzt(0,J))
               call ztozz_anomal(v(0,1,J),v(0,2,J),v(0,3,J),v(0,4,J),-ie,-iu,zzezm(0,J))
               call ztozz_anomal(v(0,1,J),v(0,2,J),v(0,5,J),v(0,6,J),-ie,-it,zzezt(0,J))
               call ztozz_anomal(v(0,3,J),v(0,4,J),v(0,5,J),v(0,6,J),-iu,-it,zzmzt(0,J))
            endif
c only higgs contributions
         call ztozzzHiggs_anomal(v(0,1,J),J,h,zzzz_higgs1(0,J),zzzz_higgs2(0,J), 
     1           zzzz_higgs3(0,J),NW) !Z -> H Z, H -> Z Z for all permutations
         endif


      else ! SM

c leptonic tensors
         if (NW.eq.1) then
c NW approximation (only Higgs contribution exists in the leptonic tensors)
            call ztozzz(v(0,1,J),J,h,zzzz(0,J),NW) !Z -> Z Z Z -> 6l
c only higgs contributions
         call ztozzzHiggs(v(0,1,J),J,h,zzzz_higgs1(0,J),zzzz_higgs2(0,J), 
     1           zzzz_higgs3(0,J),NW) !Z -> H Z, H -> Z Z for all permutations
         else
c all diagrams
            if (ps_number.eq.1) then
               call atozzz(v(0,1,J),J,h,azzz(0,J)) !A -> Z Z Z -> 6l
               call ztozzz(v(0,1,J),J,h,zzzz(0,J),NW) !Z -> Z Z Z -> 6l
               call atozz(v(0,1,J),v(0,2,J),v(0,3,J),v(0,4,J),-ie,-iu,azezm(0,J))
               call atozz(v(0,1,J),v(0,2,J),v(0,5,J),v(0,6,J),-ie,-it,azezt(0,J))
               call atozz(v(0,3,J),v(0,4,J),v(0,5,J),v(0,6,J),-iu,-it,azmzt(0,J))
               call ztozz(v(0,1,J),v(0,2,J),v(0,3,J),v(0,4,J),-ie,-iu,zzezm(0,J))
               call ztozz(v(0,1,J),v(0,2,J),v(0,5,J),v(0,6,J),-ie,-it,zzezt(0,J))
               call ztozz(v(0,3,J),v(0,4,J),v(0,5,J),v(0,6,J),-iu,-it,zzmzt(0,J))
            endif
c only higgs contributions
         call ztozzzHiggs(v(0,1,J),J,h,zzzz_higgs1(0,J),zzzz_higgs2(0,J), 
     1           zzzz_higgs3(0,J),NW) !Z -> H Z, H -> Z Z for all permutations
         endif

      endif !anom couplings or SM

         END SELECT     ! procid

      enddo !J = 1,3

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

      lnQomu(2) = log(qzzz(4,2)/mufsq(1,2))
      lnQomu(3) = log(qzzz(4,3)/mufsq(2,3))
 
      omxi(2) = 1d0-x1
      omxi(3) = 1d0-x2

c      factor 4*pi*alfas/(2x pa pb ) = 4*pi*alfas/Q^2

      tgs2oqsq(2) = 4d0*pi*als(1,2)/qzzz(4,2)
      tgs2oqsq(3) = 4d0*pi*als(2,3)/qzzz(4,3)

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

C********************************************************************

C     Now sum over the subprocesses contributing to Z Z Z production

      nmax = FL_VVVg(fsign,-1)          !reset counter for subprocesses to 0

C*******************  q1 qbar2 ---> Z Z Z g   **********************
         
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

      call qq_ZZZj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     &              pdf,res,nmin,nmax,NW,ps_number)

C*******************  qbar1 q2 ---> Z Z Z g   **********************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=1
      physToDiag(3)=3    ! gluon

      fsign(1) = 1
      fsign(2) = -1

      call qq_ZZZj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     &              pdf,res,nmin,nmax,NW,ps_number)

C*******************  g q1 ---> Z Z Z q2  ***********************************

      physToDiag(1)=3
      physToDiag(2)=1
      physToDiag(3)=2

      fsign(1) = 1
      fsign(2) = 1
      gsign    = -1

      call qq_ZZZj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1              pdf,res,nmin,nmax,NW,ps_number)

C*******************  qbar2 g ---> Z Z Z qbar1   **********************

      physToDiag(1)=2
      physToDiag(2)=3
      physToDiag(3)=1

      fsign(1) = -1
      fsign(2) = -1
      gsign = -1

      call qq_ZZZj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1              pdf,res,nmin,nmax,NW,ps_number)

C*******************  g q2bar ---> Z Z Z q1bar  **********************

      physToDiag(1)=3
      physToDiag(2)=2
      physToDiag(3)=1

      fsign(1) = -1
      fsign(2) = -1
      gsign = -1

      call qq_ZZZj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1              pdf,res,nmin,nmax,NW,ps_number)

C*******************  q1 g ---> Z Z Z q2   **********************

      physToDiag(1)=1
      physToDiag(2)=3
      physToDiag(3)=2

      fsign(1) = 1
      fsign(2) = 1

      call qq_ZZZj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1              pdf,res,nmin,nmax,NW,ps_number)


c**************  end of process evaluation part one ******************


! now run again for down-type decay products if 93 93 / 94 94 is 
! requested and the z decays hadronically.

      SELECT CASE(procid)
      CASE(ZZZhad)
       if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
        ! set couplings for Z hadronic decay into ddbar
        call setZtodd

        do j=1,3
           CALL JIOXXX(zep(1,J),zem(1,J),GZ_AF,ZERO , ZERO ,ae(1,J))
           CALL JIOXXX(zep(1,J),zem(1,J),GZ_ZF,ZMASS,ZWIDTH,ze(1,J))
            if (ps_number.eq.1) then
               call atozzz_had(v(0,1,J),J,h,azzz(0,J),NW) !A -> Z Z Z -> 6l
               call ztozzz_had(v(0,1,J),J,h,zzzz(0,J),NW) !Z -> Z Z Z -> 6l
               call zzz_atozz_had(v(0,1,J),v(0,2,J),v(0,3,J),v(0,4,J),-ie,-iu,azezm(0,J))
               call zzz_atozz_had(v(0,1,J),v(0,2,J),v(0,5,J),v(0,6,J),-ie,-it,azezt(0,J))
               call atozz_anomal(v(0,3,J),v(0,4,J),v(0,5,J),v(0,6,J),-iu,-it,azmzt(0,J))
               call zzz_ztozz_had(v(0,1,J),v(0,2,J),v(0,3,J),v(0,4,J),-ie,-iu,zzezm(0,J))
               call zzz_ztozz_had(v(0,1,J),v(0,2,J),v(0,5,J),v(0,6,J),-ie,-it,zzezt(0,J))
               call ztozz_anomal(v(0,3,J),v(0,4,J),v(0,5,J),v(0,6,J),-iu,-it,zzmzt(0,J))
            endif
c only higgs contributions
            call ztozzzHiggs_had(v(0,1,J),J,h,zzzz_higgs1(0,J),zzzz_higgs2(0,J), 
     1           zzzz_higgs3(0,J),NW) !Z -> H Z, H -> Z Z for all permutations
        enddo


C*******************  q1 qbar2 ---> Z Z Z g   **********************
         
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

      call qq_ZZZj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     &              pdf,res,nmin,nmax,NW,ps_number)

C*******************  qbar1 q2 ---> Z Z Z g   **********************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=1
      physToDiag(3)=3    ! gluon

      fsign(1) = 1
      fsign(2) = -1

      call qq_ZZZj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     &              pdf,res,nmin,nmax,NW,ps_number)

C*******************  g q1 ---> Z Z Z q2  ***********************************

      physToDiag(1)=3
      physToDiag(2)=1
      physToDiag(3)=2

      fsign(1) = 1
      fsign(2) = 1
      gsign    = -1

      call qq_ZZZj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1              pdf,res,nmin,nmax,NW,ps_number)

C*******************  qbar2 g ---> Z Z Z qbar1   **********************

      physToDiag(1)=2
      physToDiag(2)=3
      physToDiag(3)=1

      fsign(1) = -1
      fsign(2) = -1
      gsign = -1

      call qq_ZZZj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1              pdf,res,nmin,nmax,NW,ps_number)

C*******************  g q2bar ---> Z Z Z q1bar  **********************

      physToDiag(1)=3
      physToDiag(2)=2
      physToDiag(3)=1

      fsign(1) = -1
      fsign(2) = -1
      gsign = -1

      call qq_ZZZj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1              pdf,res,nmin,nmax,NW,ps_number)

C*******************  q1 g ---> Z Z Z q2   **********************

      physToDiag(1)=1
      physToDiag(2)=3
      physToDiag(3)=2

      fsign(1) = 1
      fsign(2) = 1

      call qq_ZZZj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1              pdf,res,nmin,nmax,NW,ps_number)


c***  end of extra process evaluation for hadronic decay  ***

       endif
      end select

c*****************  end of process evaluation  *********************


      SELECT CASE(procid)
      CASE(ZZZhad)
       call hadDecayFactor_Z(finalquarks(1), abs(mjj2(v(0,1,1),v(0,2,1))), N_gen_up, N_gen_down, fac_Z_up, fac_Z_down)
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


      if (init.eq.1) then
         init = init+1
      endif
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
               CASE(ZZZhad)
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

         m2s(j) = m2s(j) * 8.0d0     ! factor 2 for random helicity summation
      enddo


      RETURN
      END

c*****************************************************************************
c
c    end subroutine m2s_qqZZZj
c
c*****************************************************************************
c******************************************************************
c
c   begin subroutine qq_ZZZj
c
c*****************************************************************

      subroutine qq_ZZZj(xi,p,v,physToDiag,fsign,gsign,bos, 
     1                    nlo,lok,xuz,
     1                    pdf,res,nmin,nmax, NW,ps_number)
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"

      double precision p(0:3,max_p,max_kin), v(0:3,max_v,max_kin), xi(nx), xuz
      double precision pdf(-6:6,2,3), res(maxnumsubproc,0:3)
      integer physToDiag(3), fsign(3+n_v), 
     1        gsign, bos, nlo, nmin, nmax, NW,ps_number
      logical lok(3)
c
c qq_VVVj calls the amplitude square routines 
c             qqVVVj    for qq --> VVV g --> g 6 leptons
c                           gq --> VVV q --> q 6 leptons 
c
c for the subtraction terms for the NLO cross section calculation
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
c          NW              NW = 1: NW approximation, NW = 0: full calc.
c  OUTPUT:
c          uuZZZ(ID)   output in format as in qqZqq(j), but corresponding
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
      DOUBLE PRECISION C2, TR, N                              ! color factors
      PARAMETER(N=3d0,TR=0.5d0,c2=4d0/3d0)

c Note: Factor 3 from color sum included in qqVVV. I am using color summed
c amplitudes here, while CS use color averages. Thus I later divide by
c 8*3 for an initial qg state and by 3*3 for an initial qq state.

c alfas, scales etc

#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/coupl.inc"

c
c  helicity selection

      INTEGER h
      COMMON /hcount / h

      DOUBLE PRECISION uuzzzb(2,2:3), ddzzzb(2,2:3), uuzzz, ddzzz
      DOUBLE PRECISION uuzzz2, ddzzz2, uuzzzg, ddzzzg
      DOUBLE PRECISION pbar(0:3,3+n_v), qbar(0:4),q2,v2,sub(2:3), 
     1       dotrr, gs2(2:3), polcolq, polcolg
      DOUBLE PRECISION NCmatrixelt(0:1,3,2), temp
      INTEGER iflav(3), diagToPhys(3), FL_VVVg
      EXTERNAL dotrr, FL_VVVg

      DOUBLE COMPLEX zero
      PARAMETER (zero=(0d0,0d0))

      LOGICAL ldebug
      SAVE gs2,polcolq,polcolg,uuzzzb,ddzzzb
      INTEGER i,if1,j,k,mu
      
      PARAMETER (ldebug=.false.)

c--------------------------------------------------------------------

      uuzzz  = 0d0
      ddzzz  = 0d0
      uuzzzg = 0d0
      ddzzzg = 0d0

c  if this is a new phase space point -->
c  reset the LO amplitude to 0 to avoid wrong subtraction
      if (physToDiag(1).eq.1.and.physToDiag(2).eq.2) then
         do i = 1,2
            do j = 2,3
               uuzzzb(i,j) = 0d0
               ddzzzb(i,j) = 0d0
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
         
c get the real emission amplitude squared, store it in uuzzz, ddzzz
      do mu = 0,3
         do i = 1,3
            pbar(mu,physToDiag(i))=p(mu,i,1)
         enddo
         qbar(mu) = pbar(mu,3)
      enddo
      qbar(4)=0d0

      do mu = 0,3               ! kinematics for qq -> ZZZ -> 6 leptons
         pbar(mu,3) = v(mu,1,1)
         pbar(mu,4) = v(mu,2,1)
         pbar(mu,5) = v(mu,3,1)
         pbar(mu,6) = v(mu,4,1)
         pbar(mu,7) = v(mu,5,1)
         pbar(mu,8) = v(mu,6,1)
         pbar(mu,9) = qbar(mu)
      enddo
      fsign(9) = gsign

c calculate the matrix elements pa + pb -> 6l p1
      if (lok(1) .or. nlo.eq.0) then

#ifdef WITH_MADGRAPH
         if (ldebug) then
            call suu_zzzg(pbar,fsign,uuzzzg,NW)          !madgraph
            call sdd_zzzg(pbar,fsign,ddzzzg,NW)
         endif
#endif

         if (ps_number.eq.1) then ! subtract Higgs contr.
            call qqZZZj(pbar,fsign,qbar,gsign,uuzzz,ddzzz,NW)
            call qqZZZj1(pbar,fsign,qbar,gsign,uuzzz2,ddzzz2,NW)
            uuzzz = uuzzz - uuzzz2
            ddzzz = ddzzz - ddzzz2
         else if (ps_number.eq.2) then ! only Higgs contr.
            call qqZZZj1(pbar,fsign,qbar,gsign,uuzzz,ddzzz,NW)
         else ! complete |M|^2
            call qqZZZj(pbar,fsign,qbar,gsign,uuzzz,ddzzz,NW)
         endif

         if (ldebug) then
            if ( (abs(uuzzz/uuzzzg - 1) .gt. 1d-4) .or.
     1           (abs(ddzzz/ddzzzg - 1) .gt. 1d-4) ) then
               print*,"p(0,1,1) = ",p(0,1,1)
               print*,"uuzzz  =   ",uuzzz
               print*,"uuzzzg =  ",uuzzzg
               print*,"uuzzzg/uuzzz - 1 = ",uuzzzg/uuzzz - 1
               print*,"ddzzz =   ",ddzzz
               print*,"ddzzzg =  ",ddzzzg
               print*,"ddzzz/ddzzzg - 1 = ",ddzzz/ddzzzg - 1
               print*," "
            endif
         endif

      endif

      NCmatrixelt(0,2,1)=uuzzz
      NCmatrixelt(1,2,1)=ddzzz

c for the NLO case get the subtraction terms; 
c first the case with a final state gluon (gsign = +1)

      if (gsign.eq.1) then
         if (ldebug) then
            print*," final state gluon section in qq_ZZZj "
            print 101," fsign = ",fsign
 101                format(a,6i5,a,i5)
         endif

#ifdef WITH_NLO
         if (nlo.eq.1) then
            do j = 2,3               ! j=2: emission off 1,2 line
               do mu = 0,3           ! j=3: emission off 3,4 line
                  do i = 1,2
                     pbar(mu,physToDiag(i))=p(mu,i,j)
                  enddo
                  do i = 1,6
                     pbar(mu,i+2)=v(mu,i,j)
                  enddo
               enddo

c calculation of matrix elements for q qb -> Z Z Z g:

               if(lok(j)) then
#ifdef WITH_MADGRAPH
                  if(ldebug) then
                     call suu_zzz(pbar,fsign,uuzzz2,NW)  !madgraph
                     call sdd_zzz(pbar,fsign,ddzzz2,NW)
                  endif
#endif

                  if (ps_number.eq.1) then ! subtract Higgs contr.
                     call qqZZZ(pbar,fsign,0,J,uuzzz,ddzzz,NW)
                     call qqZZZ1(pbar,fsign,0,J,uuzzz2,ddzzz2,NW)
                     uuzzz = uuzzz - uuzzz2
                     ddzzz = ddzzz - ddzzz2
                  else if (ps_number.eq.2) then ! only Higgs contr.
                     call qqZZZ1(pbar,fsign,0,J,uuzzz,ddzzz,NW)
                  else ! complete |M|^2
                     call qqZZZ(pbar,fsign,0,J,uuzzz,ddzzz,NW)
                  endif
  
                  if (ldebug) then
                        print*,"uuzzz/uuzzz2 =   ",uuzzz/uuzzz2
                        print*,"ddzzz/ddzzz2 =   ",ddzzz/ddzzz2
                  endif
               else
                  uuzzz = 0d0
                  ddzzz = 0d0
               endif

               NCmatrixelt(0,j,2)=uuzzz   !Born for momentum
               NCmatrixelt(1,j,2)=ddzzz   !configuration 2 and 3

c save matrix elements for later use with initial gluons
               uuzzzb(physToDiag(1),j) = uuzzz
               ddzzzb(physToDiag(1),j) = ddzzz

               q2 = 2d0*xuz*dotrr(qbar,p(0,j-1,1)) !p(mu,j-1,1) is inc.quark
               v2 = 2d0*gs2(j)*                   ! 8 pi alfas
     &             ( 2d0/(1-xuz)-(1+xuz) )
               sub(j) = v2/q2

            enddo
         endif !nlo.eq.1
#endif /* WITH_NLO */

         iflav(3) = 0           ! final state gluon id
         do if1=1,nfl
            iflav(1)=if1*fsign(physToDiag(1))
            iflav(2)=if1*fsign(physToDiag(2))
               k=FL_VVVg(iflav,2)
               if(lok(1)) then
                  res(k,1)=pdf(sign1*iflav(1),1,1)  
     &                 *pdf(sign2*iflav(2),2,1)  
     &                 *NCmatrixelt(mod(if1,2),2,1)*polcolq
               else
                     res(k,1) = 0
               endif

               if ( nlo.eq.1 .and. (lok(2).or.lok(3)) ) then
#ifdef WITH_NLO
                  res(k,2) = ( pdf(sign1*iflav(1),1,2)     ! collinear
     &                    *pdf(sign2*iflav(2),2,2)
     &                    -pdf(sign1*iflav(1),1,1)  ! real emission
     &                    *pdf(sign2*iflav(2),2,2 )*sub(2))*C2
     &                    *NCmatrixelt(mod(if1,2),2,2)
                  res(k,3) = ( pdf(sign1*iflav(1),1,3)     ! collinear
     &                    *pdf(sign2*iflav(2),2,3)
     &                    -pdf(sign1*iflav(1),1,3)  ! real emission
     &                    *pdf(sign2*iflav(2),2,1)*sub(3))*C2
     &                    *NCmatrixelt(mod(if1,2),3,2)

                  res(k,2) = res(k,2) * polcolq
                  res(k,3) = res(k,3) * polcolq
                  res(k,0) = res(k,1) + res(k,2) + res(k,3)

c debugging for collinear and soft divergences -----------------------
                  if(ldebug) then                              !begin debug
                     temp = qbar(0)              !soft
                     if(temp.le.1.0d-1) then
c                     temp = abs(dotrr(qbar,pbar(0,1)))  !collinear
c                     temp2 = abs(dotrr(qbar,pbar(0,2)))
c                     if(temp.le.0.0001.or.temp2.le.0.0001) then
                        temp = abs(res(k,1)/(res(k,2)+res(k,3)))
                        if(temp.le.0.99d0.or.temp.ge.1.01d0) then
                           print*," "
                           print*,"gluon in the final state "
                           print*,"if1 = ",if1
                           print*,"qbar(0) = ",qbar(0)
                           print*,"qbar*pbar(0,1)=",dotrr(qbar,pbar(0,1))
                           print*,"qbar*pbar(0,2)=",dotrr(qbar,pbar(0,2))
                           print*,"res(k,1) = ",res(k,1)
                           print*,"res(k,2)+res(k,3)=",res(k,2)+res(k,3)
                           print*,"res(k,0) = ",res(k,0)
                           print*,"|res(k,1)/res(k,2)| = ",temp
                           print*," "
                        endif
                     endif
                  endif                                       !end debug
c----------------------------------------------------------------------
#endif /* WITH_NLO */
               else
                     res(k,0) = res(k,1)
                     res(k,2) = 0
                     res(k,3) = 0
               endif
         enddo  

      elseif (gsign.eq.-1) then           !initial gluon section
         j = 0
         if (nlo.eq.0) then
            if (physToDiag(1).eq.3) then
               j=2              ! j=2: emission off 1,2 line
            elseif (physToDiag(2).eq.3) then
               j=3              ! j=3: emission off 3,4 line
            endif
         elseif (nlo.eq.1) then
#ifdef WITH_NLO
            if (physToDiag(1).eq.3) then
               j=2              ! j=2: emission off 1,2 line
            elseif (physToDiag(2).eq.3) then
               j=3              ! j=3: emission off 3,4 line
            endif
          
               q2 = 2d0 * dotrr(p(0,3,1),p(0,j-1,1))*xuz
               v2 = 2d0*gs2(j) * ( (1-xuz)**2 + xuz**2 )
               sub(j) = v2/q2

            if (physToDiag(1).eq.1.or.PhysToDiag(2).eq.2) then
               NCmatrixelt(0,2,2)=uuzzzb(1,2)
               NCmatrixelt(1,2,2)=ddzzzb(1,2)
               NCmatrixelt(0,3,2)=uuzzzb(1,3)
               NCmatrixelt(1,3,2)=ddzzzb(1,3)
            else
               NCmatrixelt(0,2,2)=uuzzzb(2,2)
               NCmatrixelt(1,2,2)=ddzzzb(2,2)
               NCmatrixelt(0,3,2)=uuzzzb(2,3)
               NCmatrixelt(1,3,2)=ddzzzb(2,3)
            endif
#endif
         endif

         do if1=1,nfl
            iflav(j-1) = 0
            iflav(4-j)=if1*fsign(physToDiag(4-j))
            iflav(3)  =if1*fsign(physToDiag(3))
            k=FL_VVVg(iflav,j)
            if (lok(1)) then
               res(k,1)=pdf(sign1*iflav(1),1,1)  
     &              *pdf(sign2*iflav(2),2,1)      
     &              *NCmatrixelt(mod(if1,2),2,1)*polcolg
            else
               res(k,1) = 0
            endif

            if ( nlo.eq.1 .and. (lok(2).or.lok(3)) ) then
#ifdef WITH_NLO
               res(k,j) = ( pdf(sign1*iflav(1),1,j)     ! collinear
     &                 *pdf(sign2*iflav(2),2,j)
     &                 -pdf(sign1*iflav(1),1,2*j-3)  ! real emission
     &                 *pdf(sign2*iflav(2),2,4-j )*sub(j))
     &                 *NCmatrixelt(mod(if1,2),j,2)

               res(k,j) = res(k,j)*C2*polcolg
               res(k,-j+5) = 0d0

               res(k,0) = res(k,1) + res(k,j)

c***********************************************************************
               if (ldebug) then                         !begin debug
                     temp = p(0,3,1)                 !soft
                     if(temp.le.0.1d0) then
c                     temp = abs(dotrr(p(0,3,1),p(0,1,1)))  !collinear
c                     temp2 = abs(dotrr(p(0,3,1),p(0,2,1)))
c                     if(temp.le.0.0001.or.temp2.le.0.0001) then
                     if(j.eq.3) then
                        temp = abs((res(k,1))+res(k-4,1))/
     &                             (res(k,2)+res(k-4,2))
                        if(temp.le.0.9.or.temp.ge.1.1) then
                           print*," "
                           print*," gluon in the initial state"
                           print*," k = ",k," if1 = ",if1      
                           print*,"p(0,3,1) = ",p(0,3,1)
c                           print*,"pa*p1= ",dotrr(p(0,3,1),p(0,1,1))
c                           print*,"pb*p1= ",dotrr(p(0,3,1),p(0,2,1))
                           print*,"res(k,1) = ",res(k,1)
                           print*,"res(k,2) = ",res(k,2)
                     print*,"res(k,1)+res(k-4,1) = ",res(k,1)+res(k-4,1)
                     print*,"res(k,2)+res(k-4,2) = ",res(k,2)+res(k-4,2)
                           endif
                        endif
                     endif
               endif                                   !end debug
c***********************************************************************
#endif /* WITH_NLO */

            else
               res(k,0) = res(k,1)
               res(k,2) = 0
               res(k,3) = 0
            endif

         enddo
      endif

      nmax = FL_VVVg(iflav,-2)

      end

c******************************************************************
c
c   end subroutine qq_ZZZj
c
c*****************************************************************
