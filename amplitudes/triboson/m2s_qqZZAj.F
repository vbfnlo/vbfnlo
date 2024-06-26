c*********************  subroutine m2s_qqZZAj*********************************
c
c INPUT     bos       Boson identifier, 2 = Z, 3,4=W
c           nlo       NLO=1: create subtraction term; LO = 0
c           lok       need to calculate m2s(1:3)? T or F 
c           xi        Feynman x parameters of incoming fermions
c           p         momenta associated with external fermions
c           v         momenta associated with WWZ decay fermions
c           rn        random number array
c           xuz       x_ija for subtraction terms
c          
c
c OUTPUT    m2s       |M|^2*pdf1*pdf2 
c
c       Adapted from m2s_qqZZAj by V. Hankele
c
c       Giuseppe Bozzi, <giuseppe@particle.physik.uni-karlsruhe.de>
c       Last modified: December 2008
c
c***************************************************************************

      subroutine m2s_qqZZAj(bos,nlo,lok,xi,p,v,rn,xuz,m2s)

      use globalvars, only: ldoblha, lglobalprint

      implicit none

c declare input/output variables
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/process.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

c     hadronic decay
      integer N_gen_up, N_gen_down
      real*8 fac_Z_up, fac_Z_down, mjj2
      external mjj2

      integer bos
      double precision xi(nx), p(0:3,max_p,max_kin), v(0:3,max_v,max_kin), 
     1       rn(1), xuz, m2s(0:3)
      logical lok(3)

c  helicity selection
      INTEGER h
      COMMON /hcount / h

      double precision zero
      parameter (zero=0d0)

c declare external functions
      integer Fl_VVVg
      external Fl_VVVg

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

c declare parton distribution variables
      double precision pdf(-6:6,2,3)

c define program switches
      logical ldebug
      parameter (ldebug=.false.)

c variables for leshouches interface
      DOUBLE PRECISION weight,rnumb,RandomNumber

c---------------------------- BEGIN CODE --------------------------------

      if ( init .ne. 0) then
         do i = 1,nmax
            do j=0,3
               res(i,j) = 0
            enddo
         enddo
      else if ((init.eq.0).and.(.not. ldoblha)) then
         if (bos.eq.22 .or. bos.eq.22) then
            if (lglobalprint) then
               write(6,*) " "
               write(6,*) "Z Z A j amplitude square information:"
               write(6,*) "-----------------------------------------------"
               write(6,*) " "
            endif
            call vtozza_reset
         endif
         call printnfl(.false.)

         if (lglobalprint) then
            print*," "
            print*," creal = ",crealVV," cvirtual = ",cvirtVV
            print*," "
         endif
         init = 1
         do i = 1,maxnumsubproc
            do j = 0,3
               res(i,j) = 0
            enddo
         enddo
      else if ((init.eq.0).and.ldoblha) then
        init = 1
        h = 1
      endif

c -------------------------------------------------------------------

        
c for ZZA->4l precalculate A->ZZA,Z->ZZA etc..

c lepton spinors and Z/A polarisation vectors
 
c select helicity: h = random number for lepton+photon helicity (h=1:8)

      if (ldoblha) then
        h = blha_ranhelcomb
      end if

      h = mod(h,8) + 1

      ie = sign(1,4-h)
      if(bos .eq. 221) then
         iu = (-1)**(h+1) !muons
      else if(bos .eq. 251) then
         iu = 1 ! only left-handed neutrinos
      end if
      it = (-1)**((h-1)/2)

c...Les Houches interface
        if ((lha.or.hepmc).and..not.doNLO) then
           helicity(1)=-ie
           helicity(2)= ie
           helicity(3)=-iu
           helicity(4)= iu
           helicity(5) =  it
        endif

      do J = 1,3

         do mu = 0,3
            qze(mu,J) = v(mu,1,J)+v(mu,2,J)  !A,Z(e)
            qzmu(mu,J) = v(mu,3,J)+v(mu,4,J) !A,Z(m)
            qz(mu,J) = v(mu,5,J)             !A
            qzz(mu,J) = qze(mu,J) + qzmu(mu,J)             !ZZ 4-momentum
            qzza(mu,J) = qze(mu,J) + qzmu(mu,J) + qz(mu,J) !ZZA 4-momentum
         enddo

         qze(4,J) = dotrr(qze(0,J),qze(0,J))
         qzmu(4,J) = dotrr(qzmu(0,J),qzmu(0,J))
         qz(4,J) = dotrr(qz(0,J),qz(0,J))
         qzz(4,J) = dotrr(qzz(0,J),qzz(0,J))
         qzza(4,J) = dotrr(qzza(0,J),qzza(0,J))          !s
 
         CALL OXXXXX(v(0,1,J),ZERO ,-ie,+1,zem(1,J)) !e-
         CALL IXXXXX(v(0,2,J),ZERO ,+ie,-1,zep(1,J)) !e+
         CALL OXXXXX(v(0,3,J),ZERO ,-iu,+1,zmm(1,J)) !mu-
         CALL IXXXXX(v(0,4,J),ZERO ,+iu,-1,zmp(1,J)) !mu+
         CALL VXXXXX(qz(0,J),ZERO ,it,1,atau(1,J)) !A

         SELECT CASE(procid)
         CASE(ZZhadA)
           ! for general output to all combinations: up-type first
           if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
             ! set couplings for Z hadronic decay into uubar
             call setZtouu
           endif
           CALL JIOXXX(zep(1,J),zem(1,J),GZ_AF,ZERO , ZERO ,ae(1,J))
           CALL JIOXXX(zep(1,J),zem(1,J),GZ_ZF,ZMASS,ZWIDTH,ze(1,J))

         CASE DEFAULT
           CALL JIOXXX(zep(1,J),zem(1,J),GAL,ZERO , ZERO ,ae(1,J))
           CALL JIOXXX(zep(1,J),zem(1,J),GZL,ZMASS,ZWIDTH,ze(1,J))
         END SELECT
         
         if(bos .eq. 221) then
            CALL JIOXXX(zmp(1,J),zmm(1,J),GAL,ZERO,ZERO, amu(1,J))   !A*m
            CALL JIOXXX(zmp(1,J),zmm(1,J),GZL,ZMASS,ZWIDTH, zmu(1,J))!Z*m
         else if(bos .eq. 251) then
            CALL JIOXXX(zmp(1,J),zmm(1,J),GZN,ZMASS,ZWIDTH, zmu(1,J))!Z*v
            amu(1:4,J) = 0
            amu(5:6,J) = zmu(5:6,J)
         endif        

        SELECT CASE(procid)
        CASE(ZZhadA)
          ! for general output to all combinations: up-type first
          if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
            ! set couplings for Z hadronic decay into uubar
            call setZtouu
          endif

           call anomal_formfactor(qzza(0,J),qze(0,J),qzmu(0,J),qz(0,J))
           call zza_atozz_had(v(0,1,J),J,azz(0,J),0) !A -> Z Z -> e- e+ mu- mu+
           call zza_ztozz_had(v(0,1,J),J,zzztens(0,J),0) !Z -> Z Z -> e- e+ mu- mu+
           call zza_atozae_had(v(0,1,J),J,h,aaze(0,J),0) !A -> Z A -> e- e+ a
           call zza_ztozae_had(v(0,1,J),J,h,zzae(0,J),0) !Z -> Z A -> e- e+ a
           call zza_atozam_anomal(v(0,1,J),J,h,aazm(0,J),0) !A -> Z A -> mu- mu+ a
           call zza_ztozam_anomal(v(0,1,J),J,h,zzam(0,J),0) !Z -> Z A -> mu- mu+ a
           call zza_atozza_had(v(0,1,J),J,h,azza(0,J),0) !A -> ZZA -> e- e+ mu- mu+ a
           call zza_ztozza_had(v(0,1,J),J,h,zzza(0,J),0) !Z -> ZZA -> e- e+ mu- mu+ a

        CASE DEFAULT

        if (with_anom) then ! anomalous gauge boson couplings
c          using global form factor for all tensors of one phase space point
c          this ensures proper cancellations for anomalous contributions
c          energy scale is invariant ZZA mass
           call anomal_formfactor(qzza(0,J),qze(0,J),qzmu(0,J),qz(0,J))

c leptonic tensors
           call zza_atozz_anomal(v(0,1,J),J,azz(0,J),0) !A -> Z Z -> e- e+ mu- mu+
           call zza_ztozz_anomal(v(0,1,J),J,zzztens(0,J),0) !Z -> Z Z -> e- e+ mu- mu+
           call zza_atozae_anomal(v(0,1,J),J,h,aaze(0,J),0) !A -> Z A -> e- e+ a
           call zza_ztozae_anomal(v(0,1,J),J,h,zzae(0,J),0) !Z -> Z A -> e- e+ a
           call zza_atozam_anomal(v(0,1,J),J,h,aazm(0,J),0) !A -> Z A -> mu- mu+ a
           call zza_ztozam_anomal(v(0,1,J),J,h,zzam(0,J),0) !Z -> Z A -> mu- mu+ a
           call zza_atozza_anomal(v(0,1,J),J,h,azza(0,J),0) !A -> ZZA -> e- e+ mu- mu+ a
           call zza_ztozza_anomal(v(0,1,J),J,h,zzza(0,J),0) !Z -> ZZA -> e- e+ mu- mu+ a


        else  ! SM

c leptonic tensors
           if(bos .eq. 221) then !mu- mu+
            call zza_atozz(v(0,1,J),J,azz(0,J),0)      !A -> Z Z -> e- e+ mu- mu+
            call zza_ztozz(v(0,1,J),J,zzztens(0,J),0)  !Z -> Z Z -> e- e+ mu- mu+
            call zza_atozae(v(0,1,J),J,h,aaze(0,J),0)  !A -> Z A -> e- e+ a
            call zza_ztozae(v(0,1,J),J,h,zzae(0,J),0)  !Z -> Z A -> e- e+ a
            call zza_atozam(v(0,1,J),J,h,aazm(0,J),0)  !A -> Z A -> mu- mu+ a
            call zza_ztozam(v(0,1,J),J,h,zzam(0,J),0)  !Z -> Z A -> mu- mu+ a
            call zza_atozza(v(0,1,J),J,h,azza(0,J),0)  !A -> ZZA -> e- e+ mu- mu+ a
            call zza_ztozza(v(0,1,J),J,h,zzza(0,J),0)  !Z -> ZZA -> e- e+ mu- mu+ a
         else if(bos .eq. 251) then !v- v+
            call zzna_atozz(v(0,1,J),J,azz(0,J),0)      !A -> Z Z -> e- e+ v- v+
            call zzna_ztozz(v(0,1,J),J,zzztens(0,J),0)  !Z -> Z Z -> e- e+ v- v+
            call zza_atozae(v(0,1,J),J,h,aaze(0,J),0)  !A -> Z A -> e- e+ a
            call zza_ztozae(v(0,1,J),J,h,zzae(0,J),0)  !Z -> Z A -> e- e+ a
            aazm(0:3,J) = 0
            zzam(0:3,J) = 0
            call zzna_atozza(v(0,1,J),J,h,azza(0,J),0)  !A -> ZZA -> v- v+ mu- mu+ a
            call zzna_ztozza(v(0,1,J),J,h,zzza(0,J),0)  !Z -> ZZA -> v- v+ mu- mu+ a
         endif

        endif ! anom or SM

        end select

      enddo

      if (ldoblha) then
        physToDiag(1)=blha_physdiagmap(1,blha_cursubproc,blha_curproc)
        physToDiag(2)=blha_physdiagmap(2,blha_cursubproc,blha_curproc)
        physToDiag(3)=blha_physdiagmap(3,blha_cursubproc,blha_curproc)

        fsign(1) = blha_fsign(1, blha_cursubproc, blha_curproc)
        fsign(2) = blha_fsign(2, blha_cursubproc, blha_curproc)
        fsign(3) = 1          !outgoing ve
        fsign(4) = -1         !outgoing ep
        fsign(5) = 1          !outgoing mu-
        fsign(6) = -1         !outgoing vm~
        fsign(7) = 1          !outgoing A
        gsign    = blha_fsign(3, blha_cursubproc, blha_curproc)

        call qq_VVV2j(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     &                pdf,res,nmin,nmax)
        return
      end if

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

      lnQomu(2) = log(qzza(4,2)/mufsq(1,2))
      lnQomu(3) = log(qzza(4,3)/mufsq(2,3))
 
      omxi(2) = 1d0-x1
      omxi(3) = 1d0-x2

c      factor 4*pi*alfas/(2x pa pb ) = 4*pi*alfas/Q^2

      tgs2oqsq(2) = 4d0*pi*als(1,2)/qzza(4,2)
      tgs2oqsq(3) = 4d0*pi*als(2,3)/qzza(4,3)

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
     1           + 2d0*ln1mxi**2) + crealVV
            do i = 1,5
               pdf(-i,j-1,j)=(pdf(-i,j-1,j)*(Dxi(j)/omxi(j)-Bx(j)) +
     1                     pdf(-i,j-1,1)*(Bx(j)+Cx(j)))*tgs2oqsq(j)
     2                     *1d0/(1d0-z)
               pdf(i,j-1,j)= (pdf(i,j-1,j)*(Dxi(j)/omxi(j)-Bx(j)) +
     1                     pdf(i,j-1,1)*(Bx(j)+Cx(j)))*tgs2oqsq(j)
     2                     *1d0/(1d0-z)
            enddo
            pdf(0,j-1,j)=pdf(0,j-1,1)*Ax(j)*tgs2oqsq(j)*1d0/(1d0-z)
         endif
      enddo
      endif

C********************************************************************

C     Now sum over the subprocesses contributing to Z Z A production

      nmax = Fl_VVVg(fsign,-1) !reset counter for subprocesses to 0

C*******************  q1 qbar2 ---> Z Z A g   **********************
         
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

      call qq_VVV2j(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     &              pdf,res,nmin,nmax)

C*******************  qbar1 q2 ---> Z Z A g   **********************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=1
      physToDiag(3)=3    ! gluon

      fsign(1) = 1
      fsign(2) = -1
      gsign    = 1

      call qq_VVV2j(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     &              pdf,res,nmin,nmax)

C*******************  g q1 ---> Z Z A q2  ***********************************

      physToDiag(1)=3
      physToDiag(2)=1
      physToDiag(3)=2

      fsign(1) = 1
      fsign(2) = 1
      gsign    = -1

      call qq_VVV2j(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     1              pdf,res,nmin,nmax)

C*******************  qbar2 g ---> Z Z A qbar1   **********************

      physToDiag(1)=2
      physToDiag(2)=3
      physToDiag(3)=1

      fsign(1) = -1
      fsign(2) = -1
      gsign    = -1

      call qq_VVV2j(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     1              pdf,res,nmin,nmax)

C*******************  g q2bar ---> Z Z A q1bar  **********************

      physToDiag(1)=3
      physToDiag(2)=2
      physToDiag(3)=1

      fsign(1) = -1
      fsign(2) = -1
      gsign    = -1

      call qq_VVV2j(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     1              pdf,res,nmin,nmax)

C*******************  q1 g ---> Z Z A q2   **********************

      physToDiag(1)=1
      physToDiag(2)=3
      physToDiag(3)=2

      fsign(1) = 1
      fsign(2) = 1
      gsign    = -1

      call qq_VVV2j(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     1              pdf,res,nmin,nmax)


c**************  end of process evaluation part one ******************


! now run again for down-type decay products if 93 93 / 94 94 is 
! requested and the z decays hadronically.

      SELECT CASE(procid)
      CASE(ZZhadA)
       if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
        ! set couplings for Z hadronic decay into ddbar
        call setZtodd

        do J = 1,3
           CALL JIOXXX(zep(1,J),zem(1,J),GZ_AF,ZERO , ZERO ,ae(1,J))
           CALL JIOXXX(zep(1,J),zem(1,J),GZ_ZF,ZMASS,ZWIDTH,ze(1,J))
           call zza_atozz_had(v(0,1,J),J,azz(0,J),0) !A -> Z Z -> e- e+ mu- mu+
           call zza_ztozz_had(v(0,1,J),J,zzztens(0,J),0) !Z -> Z Z -> e- e+ mu- mu+
           call zza_atozae_had(v(0,1,J),J,h,aaze(0,J),0) !A -> Z A -> e- e+ a
           call zza_ztozae_had(v(0,1,J),J,h,zzae(0,J),0) !Z -> Z A -> e- e+ a
           call zza_atozam_anomal(v(0,1,J),J,h,aazm(0,J),0) !A -> Z A -> mu- mu+ a
           call zza_ztozam_anomal(v(0,1,J),J,h,zzam(0,J),0) !Z -> Z A -> mu- mu+ a
           call zza_atozza_had(v(0,1,J),J,h,azza(0,J),0) !A -> ZZA -> e- e+ mu- mu+ a
           call zza_ztozza_had(v(0,1,J),J,h,zzza(0,J),0) !Z -> ZZA -> e- e+ mu- mu+ a
        enddo

C*******************  q1 qbar2 ---> Z Z A g   **********************
         
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

      call qq_VVV2j(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     &              pdf,res,nmin,nmax)

C*******************  qbar1 q2 ---> Z Z A g   **********************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=1
      physToDiag(3)=3    ! gluon

      fsign(1) = 1
      fsign(2) = -1
      gsign    = 1

      call qq_VVV2j(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     &              pdf,res,nmin,nmax)

C*******************  g q1 ---> Z Z A q2  ***********************************

      physToDiag(1)=3
      physToDiag(2)=1
      physToDiag(3)=2

      fsign(1) = 1
      fsign(2) = 1
      gsign    = -1

      call qq_VVV2j(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     1              pdf,res,nmin,nmax)

C*******************  qbar2 g ---> Z Z A qbar1   **********************

      physToDiag(1)=2
      physToDiag(2)=3
      physToDiag(3)=1

      fsign(1) = -1
      fsign(2) = -1
      gsign    = -1

      call qq_VVV2j(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     1              pdf,res,nmin,nmax)

C*******************  g q2bar ---> Z Z A q1bar  **********************

      physToDiag(1)=3
      physToDiag(2)=2
      physToDiag(3)=1

      fsign(1) = -1
      fsign(2) = -1
      gsign    = -1

      call qq_VVV2j(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     1              pdf,res,nmin,nmax)

C*******************  q1 g ---> Z Z A q2   **********************

      physToDiag(1)=1
      physToDiag(2)=3
      physToDiag(3)=2

      fsign(1) = 1
      fsign(2) = 1
      gsign    = -1

      call qq_VVV2j(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     1              pdf,res,nmin,nmax)


c***  end of extra process evaluation for hadronic decay  ***

       endif
      end select


c*****************  end of process evaluation  *********************

      if (init.eq.1) then
         init = init+1
      endif

      SELECT CASE(procid)
      CASE(ZZhadA)
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
               CASE(ZZhadA)
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

         if (bos .eq. 221) then
           m2s(j) = m2s(j) * 8d0           ! factor 2 for random helicity summation
         else if (bos .eq. 251) then
           m2s(j) = m2s(j) * 4d0           ! factor 2 for random helicity summation - except Z->nu nu~
         endif

      enddo


      RETURN
      END

c*****************************************************************************
c
c    end subroutine m2s_qqZZAj
c
c*****************************************************************************
c******************************************************************
c
c   begin subroutine qq_ZZAj
c
c*****************************************************************

      subroutine qq_VVV2j(xi,p,v,physToDiag,fsign,gsign,bos, 
     1                    nlo,lok,xuz,
     1                    pdf,res,nmin,nmax)

      use globalvars, only: ldoblha

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"

      double precision p(0:3,max_p,max_kin), v(0:3,max_v,max_kin), xi(nx), xuz
      double precision pdf(-6:6,2,3), res(maxnumsubproc,0:3)
      integer physToDiag(3), fsign(3+n_v), 
     1        gsign, bos, nlo, nmin, nmax
      logical lok(3)
c
c qq_VVV2j calls the amplitude square routines 
c             qqZZAj    for qq --> ZZA g --> g 4 leptons
c                           gq --> ZZA q --> q 4 leptons 
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
c  OUTPUT:
c          uuzza(ID)   output in format as in qqZqq(j), but corresponding
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
      double precision C2, TR, N                              ! color factors
      parameter(N=3d0,TR=0.5d0,c2=4d0/3d0)

c Note: Factor 3 from color sum included in qqVVV. I am using color summed
c amplitudes here, while CS use color averages. Thus I later divide by
c 8*3 for an initial qg state and by 3*3 for an initial qq state.

c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/coupl.inc"

c  helicity selection
      INTEGER h
      COMMON /hcount / h

      double precision uuzzab(2,2:3), ddzzab(2,2:3), uuzza, ddzza
      double precision uuzza2, ddzza2, uuzzag, ddzzag
      double precision pbar(0:3,3+n_v), qbar(0:4),q2,v2,sub(2:3), 
     1       dotrr, gs2(2:3), polcolq, polcolg
      double precision NCmatrixelt(0:1,3,2), temp, temp1, temp2, temp3
      integer iflav(3), diagToPhys(3), Fl_Vvvg
      external dotrr, Fl_Vvvg

      double precision zero
      parameter (zero=0d0)

      logical ldebug, ldebugsoft
      save gs2,polcolq,polcolg,uuzzab,ddzzab
      integer i,if1,j,k,mu
      
      parameter (ldebug=.false.)
      parameter (ldebugsoft=.false.)

c--------------------------------------------------------------------

      uuzza  = 0d0
      ddzza  = 0d0
      uuzzag = 0d0
      ddzzag = 0d0

c  if this is a new phase space point -->
c  reset the LO amplitude to 0 to avoid wrong subtraction
      if (physToDiag(1).eq.1.and.physToDiag(2).eq.2) then
         do i = 1,2
            do j = 2,3
               uuzzab(i,j) = 0d0
               ddzzab(i,j) = 0d0
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
         
c get the real emission amplitude squared, store it in uuzza, ddzza
      do mu = 0,3
         do i = 1,3
            pbar(mu,physToDiag(i))=p(mu,i,1)
         enddo
         qbar(mu) = pbar(mu,3)
      enddo
      qbar(4)=0d0

      do mu = 0,3               ! kinematics for qq -> ZZA -> 4 leptons + A
         pbar(mu,3) = v(mu,1,1)
         pbar(mu,4) = v(mu,2,1)
         pbar(mu,5) = v(mu,3,1)
         pbar(mu,6) = v(mu,4,1)
         pbar(mu,7) = v(mu,5,1)
         pbar(mu,8) = qbar(mu)
      enddo
      fsign(8) = gsign

c calculate the matrix elements pa + pb -> 4l a p1
      if (lok(1) .or. nlo.eq.0) then

#ifdef WITH_MADGRAPH
         if (ldebug) then
            call suu_zzag(pbar,fsign,uuzzag) !madgraph
            call sdd_zzag(pbar,fsign,ddzzag)
         endif
#endif
         
         call qqZZAj(pbar,fsign,qbar,gsign,uuzza,ddzza)

#ifdef WITH_MADGRAPH
         if (ldebug) then
            if ( (abs(uuzza/uuzzag - 1) .gt. 1d-4) .or.
     1           (abs(ddzza/ddzzag - 1) .gt. 1d-4) ) then
               print*,"p(0,1,1) = ",p(0,1,1)
               print*,"uuzza  =   ",uuzza
               print*,"uuzzag =  ",uuzzag
               print*,"uuzzag/uuzza - 1 = ",uuzzag/uuzza - 1
               print*,"ddzza =   ",ddzza
               print*,"ddzzag =  ",ddzzag
               print*,"ddzza/ddzzag - 1 = ",ddzzag/ddzza - 1
               print*," "
            endif
         endif
#endif

      endif

      if (ldoblha) then
        return
      end if

      NCmatrixelt(0,2,1)=uuzza
      NCmatrixelt(1,2,1)=ddzza

c for the NLO case get the subtraction terms; 
c first the case with a final state gluon (gsign = +1)
c************************************************************
c
c final gluon section
c
c************************************************************

      if (gsign.eq.1) then

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

c calculation of matrix elements for q qb -> Z Z A g:

               if(lok(j)) then
#ifdef WITH_MADGRAPH
                  if(ldebug) then
                     call suu_zza(pbar,fsign,uuzza2)  !madgraph
                     call sdd_zza(pbar,fsign,ddzza2)
                  endif
#endif

                  call QQZZA(pbar,fsign,0,J,uuzza,ddzza)
  
#ifdef WITH_MADGRAPH
                  if (ldebug) then
                     print*,"ptilde-born"
                     print*,"uuzza/uuzza2 =   ",uuzza/uuzza2
                     print*,"ddzza/ddzza2 =   ",ddzza/ddzza2
                  endif
#endif
               else
                  uuzza = 0d0
                  ddzza = 0d0
               endif

               NCmatrixelt(0,j,2)=uuzza   !Born for momentum
               NCmatrixelt(1,j,2)=ddzza   !configuration 2 and 3

c save matrix elements for later use with initial gluons
               uuzzab(physToDiag(1),j) = uuzza
               ddzzab(physToDiag(1),j) = ddzza

               q2 = 2d0*xuz*dotrr(qbar,p(0,j-1,1)) !p(mu,j-1,1) is inc.quark
               v2 = 2d0*gs2(j)*                   ! 8 pi alfas
     &             ( 2d0/(1-xuz)-(1+xuz) )
               sub(j) = v2/q2

            enddo !j = 2,3
         endif !nlo.eq.1
#endif /* WITH_NLO */

         iflav(3) = 0           ! final state gluon id
         do if1=1,nfl
            iflav(1)=if1*fsign(physToDiag(1))
            iflav(2)=if1*fsign(physToDiag(2))
               k=Fl_Vvvg(iflav,2)
               if(lok(1)) then
                  res(k,1)=pdf(sign1*iflav(1),1,1)  
     &                 *pdf(sign2*iflav(2),2,1)  
     &                 *NCmatrixelt(mod(if1,2),2,1)*polcolq
               else
                     res(k,1) = 0
               endif

               if ( nlo.eq.1 .and. (lok(2).or.lok(3)) ) then
#ifdef WITH_NLO
                  res(k,2) = ( pdf(sign1*iflav(1),1,2) ! collinear
     &                 *pdf(sign2*iflav(2),2,2)
     &                 -pdf(sign1*iflav(1),1,1) ! real emission
     &                 *pdf(sign2*iflav(2),2,2 )*sub(2))*C2
     &                 *NCmatrixelt(mod(if1,2),2,2)
                  res(k,3) = ( pdf(sign1*iflav(1),1,3) ! collinear
     &                 *pdf(sign2*iflav(2),2,3)
     &                 -pdf(sign1*iflav(1),1,3) ! real emission
     &                 *pdf(sign2*iflav(2),2,1)*sub(3))*C2
     &                 *NCmatrixelt(mod(if1,2),3,2)
                  
                  res(k,2) = res(k,2) * polcolq
                  res(k,3) = res(k,3) * polcolq
                  res(k,0) = res(k,1) + res(k,2) + res(k,3)

c debugging for collinear and soft divergences -----------------------
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
c                           read(*,*)
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
c                           read(*,*)
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
         enddo   !if1  

c************************************************************
c
c initial gluon section
c
c************************************************************

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
               NCmatrixelt(0,2,2)=uuzzab(1,2)
               NCmatrixelt(1,2,2)=ddzzab(1,2)
               NCmatrixelt(0,3,2)=uuzzab(1,3)
               NCmatrixelt(1,3,2)=ddzzab(1,3)
            else
               NCmatrixelt(0,2,2)=uuzzab(2,2)
               NCmatrixelt(1,2,2)=ddzzab(2,2)
               NCmatrixelt(0,3,2)=uuzzab(2,3)
               NCmatrixelt(1,3,2)=ddzzab(2,3)
            endif
#endif /* WITH_NLO */
         endif

         do if1=1,nfl
            iflav(j-1) = 0
            iflav(4-j)=if1*fsign(physToDiag(4-j))
            iflav(3)  =if1*fsign(physToDiag(3))
            k=Fl_VVVg(iflav,j)
            if (lok(1)) then
               res(k,1)=pdf(sign1*iflav(1),1,1)  
     &              *pdf(sign2*iflav(2),2,1)      
     &              *NCmatrixelt(mod(if1,2),2,1)*polcolg
            else
               res(k,1) = 0
            endif

            if ( nlo.eq.1 .and. (lok(2).or.lok(3)) ) then
#ifdef WITH_NLO
               res(k,j) = ( pdf(sign1*iflav(1),1,j) ! collinear
     &              *pdf(sign2*iflav(2),2,j)
     &              -pdf(sign1*iflav(1),1,2*j-3) ! real emission
     &              *pdf(sign2*iflav(2),2,4-j )*sub(j))
     &              *NCmatrixelt(mod(if1,2),j,2)

               res(k,j) = res(k,j)*C2*polcolg
               res(k,-j+5) = 0d0
               
               res(k,0) = res(k,1) + res(k,j)
                  if (ldebugsoft) then !begin debug
                     temp1 = abs(dotrr(p(0,3,1),p(0,1,1)))
                     temp2 = abs(dotrr(p(0,3,1),p(0,2,1)))
                     temp3 = abs(res(k,1)/res(k,j))
                     if (abs(1-temp3).ge.1d-3) then
                        if(temp1.le.1d-1.or.temp2.le.1d-1) then
                           print*," "
                           print*,"collinear quark (init state gluon)"
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
c                           read(*,*)
                        endif
                     endif
                  endif         !end debug
#endif /* WITH_NLO */
               
            else
               res(k,0) = res(k,1)
               res(k,2) = 0
               res(k,3) = 0
            endif !nlo .and. lok
            
         enddo !if1
      endif !gsign
      
      nmax = Fl_VVVg(iflav,-2)

      end

c******************************************************************
c
c   end subroutine qq_VVV2j
c
c*****************************************************************
