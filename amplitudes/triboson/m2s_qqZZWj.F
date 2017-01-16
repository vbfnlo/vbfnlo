c*********************  subroutine m2s_qqZZWj*********************************
c
c       Vera Hankele, <vera@particle.uni-karlsruhe.de>
c	Initial version:  2008 February
c	Last modified:  2009 January 28
c
c       this subroutine computes |ME|^2 for qq -> ZZW g
c       and all the crossed processes including the gluon induced
c       subprocesses at LO and the virtual corrections via call of "qqZZWj.f"
c
c	For the decay Z -> l- l+ random helicity summation is employed
c
c INPUT     bos       Boson identifier, 2 = Z, 3,4=W
c           nlo       NLO=1: create subtraction term; LO = 0
c           lok       need to calculate m2s(1:3)? T or F 
c           xi        Feynman x parameters of incoming fermions
c           p         momenta associated with external fermions
c           v         momenta associated with WWZ decay fermions
c           rn        random number array
c           xuz       x_ija for subtraction terms
c           NW        = 1:  NW Approximation, = 0: Full calculation
c          
c
c OUTPUT    m2s       |M|^2*pdf1*pdf2 
c
c***************************************************************************

      SUBROUTINE m2s_qqZZWj(bos,nlo,lok,xi,p,v,rn,xuz,m2s, NW)

      IMPLICIT NONE

c declare input/output variables
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/process.inc"

c     variables for hadronic decay
      integer N_gen_W
      real*8 fac_W, mjj2
      external mjj2
      integer N_gen_up, N_gen_down
      real*8 fac_Z_up, fac_Z_down

      integer bos
      double precision xi(nx), p(0:3,max_p,max_kin), v(0:3,max_v,max_kin), 
     1       rn(1), xuz, m2s(0:3)
      logical lok(3)
      integer NW    !in:  = 1:  NW Approximation, = 0: Full calculation

c  helicity selection
      INTEGER h
      COMMON /hcount / h

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

c declare parton distribution variables
      double precision pdf(-6:6,2,3)

c define program switches
      logical ldebug
      parameter (ldebug=.false.)

      double precision weight,rnumb,RandomNumber

c---------------------------- BEGIN CODE --------------------------------

      if ( init .ne. 0) then
         do i = 1,nmax
            do j=0,3
               res(i,j) = 0
            enddo
         enddo
      else if (init.eq.0) then
         if (bos.eq.223) then
            write(6,*) " "
            write(6,*) "Z Z W+ j amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
         else if (bos.eq.224) then
            write(6,*) " "
            write(6,*) "Z Z W- j amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
         endif
         call vtozzw_reset

         call printnfl(.true.)
         print*," "
         print*," creal = ",crealVV," cvirtual = ",cvirtVV
         print*," "
         h = 1
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
            qze(mu,J) = v(mu,1,J)+v(mu,2,J)   !Z
            qzmu(mu,J) = v(mu,3,J)+v(mu,4,J)  !Z
            qwbos(mu,J) = v(mu,5,J)+v(mu,6,J) !W
            qzzw(mu,J) = qze(mu,J) + qzmu(mu,J) + qwbos(mu,J)
         enddo

         qze(4,J) = dotrr(qze(0,J),qze(0,J))
         qzmu(4,J) = dotrr(qzmu(0,J),qzmu(0,J))
         qwbos(4,J) = dotrr(qwbos(0,J),qwbos(0,J))
         qzzw(4,J) = dotrr(qzzw(0,J),qzzw(0,J))        
      enddo !J = 1,3
 
c for ZZW->6l precalculate A->ZZW,Z->ZZW etc..

c lepton spinors and Z/A polarisation vectors
c lepton helicities for Z -> l- l+ not fixed
c -> sum over all possible helicities in |M|**2
 
c select helicity: h ... random number for lepton helicity
c                        combination (h=1:2)
      h = mod(h,4) + 1
      ie = sign(1,2-h)
      iu = (-1)**(h+1)

c...Les Houches interface
        if ((lha.or.hepmc).and..not.doNLO) then
           helicity(1)=-ie
           helicity(2)= ie
           helicity(3)=-iu
           helicity(4)= iu
           helicity(5)=-1
           helicity(6)= 1
        endif

      do J = 1,3

         CALL OXXXXX(v(0,1,J),ZERO ,-ie,+1,zem(1,J)) !e-
         CALL IXXXXX(v(0,2,J),ZERO ,+ie,-1,zep(1,J)) !e+
         CALL OXXXXX(v(0,3,J),ZERO ,-iu,+1,zmm(1,J)) !mu-
         CALL IXXXXX(v(0,4,J),ZERO ,+iu,-1,zmp(1,J)) !mu+
         CALL OXXXXX(v(0,5,J),ZERO , -1,+1,ztm(1,J)) !vt  !tau-
         CALL IXXXXX(v(0,6,J),ZERO ,  1,-1,ztp(1,J)) !ta+ !vta~

         SELECT CASE(procid)
         CASE(ZZhadWP,ZZhadWM)
            ! for general output to all combinations: up-type first
            if ((finalquarks(1).eq.93 .or. finalquarks(1).eq.94)) then
              ! set couplings for Z hadronic decay into uubar
              call setZtouu
            endif
            CALL JIOXXX(zep(1,J),zem(1,J),GZ_AF,ZERO , ZERO ,ae(1,J))   !A
            CALL JIOXXX(zep(1,J),zem(1,J),GZ_ZF,ZMASS,ZWIDTH,ze(1,J))   !Z

         CASE DEFAULT
            CALL JIOXXX(zep(1,J),zem(1,J),GAL,ZERO , ZERO ,ae(1,J))   !A
            CALL JIOXXX(zep(1,J),zem(1,J),GZL,ZMASS,ZWIDTH,ze(1,J))   !Z

         END SELECT
         CALL JIOXXX(zmp(1,J),zmm(1,J),GAL,ZERO ,ZERO  ,amu(1,J))  !A
         CALL JIOXXX(zmp(1,J),zmm(1,J),GZL,ZMASS,ZWIDTH,zmu(1,J))  !Z
         CALL JIOXXX(ztp(1,J),ztm(1,J),GWF,WMASS,WWIDTH,wbos(1,J)) !W+ !W-

c leptonic tensors

         SELECT CASE(procid)
         CASE(ZZhadWP, ZZhadWM)
               call anomal_formfactor(qzzw(0,J),qze(0,J),qzmu(0,J),qwbos(0,J))
               if (bos.eq.223) then
                 call wptozzw_had(v(0,1,J),J,h,wzzw(0,J), NW) ! W+ -> ZZW
                 call zzw_wptowz_had(v(0,1,J),v(0,2,J),v(0,5,J),v(0,6,J),J,-ie,wzw1(0,J), NW)
                 call zzw_wptowz_anomal(v(0,3,J),v(0,4,J),v(0,5,J),v(0,6,J),J,-iu,wzw2(0,J), NW)       ! no hadronic decay of the second Z
               else if (bos.eq.224) then
                 call wmtozzw_had(v(0,1,J),J,h,wzzw(0,J), NW) ! W- -> ZZW
                 call zzw_wmtowz_had(v(0,1,J),v(0,2,J),v(0,5,J),v(0,6,J),J,-ie,wzw1(0,J), NW)
                 call zzw_wmtowz_anomal(v(0,3,J),v(0,4,J),v(0,5,J),v(0,6,J),J,-iu,wzw2(0,J), NW)       ! no hadronic decay of the second Z
               endif
               call zzw_atozz_had(v(0,1,J),v(0,2,J),v(0,3,J),v(0,4,J),-ie,-iu,azezm(0,J))
               call zzw_ztozz_had(v(0,1,J),v(0,2,J),v(0,3,J),v(0,4,J),-ie,-iu,zzezm(0,J))

         CASE(ZZWPhad, ZZWMhad)
               call anomal_formfactor(qzzw(0,J),qze(0,J),qzmu(0,J),qwbos(0,J))
               if (bos.eq.223) then
                 call wptozzw_had(v(0,1,J),J,h,wzzw(0,J), NW) ! W+ -> ZZW
                 call zzw_wptowz_had(v(0,1,J),v(0,2,J),v(0,5,J),v(0,6,J),J,-ie,wzw1(0,J), NW)
                 call zzw_wptowz_had(v(0,3,J),v(0,4,J),v(0,5,J),v(0,6,J),J,-iu,wzw2(0,J), NW)
               else if (bos.eq.224) then
                 call wmtozzw_had(v(0,1,J),J,h,wzzw(0,J), NW) ! W- -> ZZW
                 call zzw_wmtowz_had(v(0,1,J),v(0,2,J),v(0,5,J),v(0,6,J),J,-ie,wzw1(0,J), NW)
                 call zzw_wmtowz_had(v(0,3,J),v(0,4,J),v(0,5,J),v(0,6,J),J,-iu,wzw2(0,J), NW)
               endif
               call zzw_atozz_had(v(0,1,J),v(0,2,J),v(0,3,J),v(0,4,J),-ie,-iu,azezm(0,J))
               call zzw_ztozz_had(v(0,1,J),v(0,2,J),v(0,3,J),v(0,4,J),-ie,-iu,zzezm(0,J))

         CASE DEFAULT
            if (with_anom) then ! anomalous gauge boson couplings
c              using global form factor for all tensors of one phase space point
c              this ensures proper cancellations for anomalous contributions
c              energy scale is invariant ZZW mass
               call anomal_formfactor(qzzw(0,J),qze(0,J),qzmu(0,J),qwbos(0,J))
               if (bos.eq.223) then
                 call wptozzw_anomal(v(0,1,J),J,h,wzzw(0,J), NW) ! W+ -> ZZW
                 call zzw_wptowz_anomal(v(0,1,J),v(0,2,J),v(0,5,J),v(0,6,J),J,-ie,wzw1(0,J), NW)
                 call zzw_wptowz_anomal(v(0,3,J),v(0,4,J),v(0,5,J),v(0,6,J),J,-iu,wzw2(0,J), NW)
               else if (bos.eq.224) then
                 call wmtozzw_anomal(v(0,1,J),J,h,wzzw(0,J), NW) ! W- -> ZZW
                 call zzw_wmtowz_anomal(v(0,1,J),v(0,2,J),v(0,5,J),v(0,6,J),J,-ie,wzw1(0,J), NW)
                 call zzw_wmtowz_anomal(v(0,3,J),v(0,4,J),v(0,5,J),v(0,6,J),J,-iu,wzw2(0,J), NW)
               endif
               if (NW.eq.0) then
                  call atozz_anomal(v(0,1,J),v(0,2,J),v(0,3,J),v(0,4,J),-ie,-iu,azezm(0,J))
                  call ztozz_anomal(v(0,1,J),v(0,2,J),v(0,3,J),v(0,4,J),-ie,-iu,zzezm(0,J))
               endif
            else if (with_kk) then ! kk scenario
#ifdef WITH_KK
               if (bos.eq.223) then
                  call wptozzw_kk(v(0,1,J),J,h,wzzw(0,J), NW) ! W+ -> ZZW
                  call zzw_wptowz(v(0,1,J),v(0,2,J),v(0,5,J),v(0,6,J),J,-ie,wzw1(0,J), NW)
                  call zzw_wptowz(v(0,3,J),v(0,4,J),v(0,5,J),v(0,6,J),J,-iu,wzw2(0,J), NW)
               else if (bos.eq.224) then
                  call wmtozzw_kk(v(0,1,J),J,h,wzzw(0,J), NW) ! W- -> ZZW
                  call zzw_wmtowz(v(0,1,J),v(0,2,J),v(0,5,J),v(0,6,J),J,-ie,wzw1(0,J), NW)
                  call zzw_wmtowz(v(0,3,J),v(0,4,J),v(0,5,J),v(0,6,J),J,-iu,wzw2(0,J), NW)
               endif
               if (NW.eq.0) then
                  call atozz(v(0,1,J),v(0,2,J),v(0,3,J),v(0,4,J),-ie,-iu,azezm(0,J))
                  call ztozz(v(0,1,J),v(0,2,J),v(0,3,J),v(0,4,J),-ie,-iu,zzezm(0,J))
               endif
#endif
            else
               if (bos.eq.223) then
                 call wptozzw(v(0,1,J),J,h,wzzw(0,J), NW) ! W+ -> ZZW
                 call zzw_wptowz(v(0,1,J),v(0,2,J),v(0,5,J),v(0,6,J),J,-ie,wzw1(0,J), NW)
                 call zzw_wptowz(v(0,3,J),v(0,4,J),v(0,5,J),v(0,6,J),J,-iu,wzw2(0,J), NW)
               else if (bos.eq.224) then
                 call wmtozzw(v(0,1,J),J,h,wzzw(0,J), NW) ! W- -> ZZW
                 call zzw_wmtowz(v(0,1,J),v(0,2,J),v(0,5,J),v(0,6,J),J,-ie,wzw1(0,J), NW)
                 call zzw_wmtowz(v(0,3,J),v(0,4,J),v(0,5,J),v(0,6,J),J,-iu,wzw2(0,J), NW)
               endif
               if (NW.eq.0) then
                  call atozz(v(0,1,J),v(0,2,J),v(0,3,J),v(0,4,J),-ie,-iu,azezm(0,J))
                  call ztozz(v(0,1,J),v(0,2,J),v(0,3,J),v(0,4,J),-ie,-iu,zzezm(0,J))
               endif
            endif
         END SELECT


      enddo  !J = 1,3

c**********************************************************************
c call PDF subroutine in order to determine parton
c distributions in the incoming (anti)protons.
c for the NLO contributions x1 = x*y in my notes with x=xuz, y=xi(i)
c**********************************************************************
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
c*********************************************************************
c and fill the coefficient functions for the finite subtraction terms
C Note that the color factors TR and C2 are NOT!!! included here
c*********************************************************************
      lnQomu(2) = log(qzzw(4,2)/mufsq(1,2))
      lnQomu(3) = log(qzzw(4,3)/mufsq(2,3))
 
      omxi(2) = 1d0-x1
      omxi(3) = 1d0-x2

c factor 4*pi*alfas/(2x pa pb ) = 4*pi*alfas/Q^2
      tgs2oqsq(2) = 4d0*pi*als(1,2)/qzzw(4,2)
      tgs2oqsq(3) = 4d0*pi*als(2,3)/qzzw(4,3)

      z = xuz

      do J = 2,3

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

      enddo !J = 2,3
      endif

c********************************************************************
c
c Now sum over the subprocesses contributing to Z Z W production
c
c********************************************************************
      nmax = FL_VVVg(fsign,-1)   !reset counter for subprocesses to 0

C*******************  q1 qbar2 ---> Z Z W g   **********************
         
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

      call qq_ZZWj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     &              pdf,res,nmin,nmax, NW)

C*******************  qbar1 q2 ---> Z Z W g   **********************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=1
      physToDiag(3)=3    ! gluon

      fsign(1) = 1
      fsign(2) = -1
      gsign    = 1

      call qq_ZZWj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     &              pdf,res,nmin,nmax, NW)

C*******************  g q1 ---> Z Z W q2  ***********************************

      physToDiag(1)=3
      physToDiag(2)=1
      physToDiag(3)=2

      fsign(1) = 1
      fsign(2) = 1
      gsign    = -1

      call qq_ZZWj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1              pdf,res,nmin,nmax, NW)

C*******************  qbar2 g ---> Z Z W qbar1   **********************

      physToDiag(1)=2
      physToDiag(2)=3
      physToDiag(3)=1

      fsign(1) = -1
      fsign(2) = -1
      gsign    = -1

      call qq_ZZWj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1              pdf,res,nmin,nmax, NW)

C*******************  g q2bar ---> Z Z W q1bar  **********************

      physToDiag(1)=3
      physToDiag(2)=2
      physToDiag(3)=1

      fsign(1) = -1
      fsign(2) = -1
      gsign    = -1

      call qq_ZZWj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1              pdf,res,nmin,nmax, NW)

C*******************  q1 g ---> Z Z W q2   **********************

      physToDiag(1)=1
      physToDiag(2)=3
      physToDiag(3)=2

      fsign(1) = 1
      fsign(2) = 1
      gsign    = -1

      call qq_ZZWj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1              pdf,res,nmin,nmax, NW)

c**************  end of process evaluation part one ******************


! now run again for down-type decay products if 93 93 / 94 94 is 
! requested and the z decays hadronically.

      SELECT CASE(procid)
      CASE(ZZhadWP, ZZhadWM)
       if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
        ! set couplings for Z hadronic decay into ddbar
        call setZtodd
        do J=1,3
           CALL JIOXXX(zep(1,J),zem(1,J),GZ_AF,ZERO , ZERO ,ae(1,J))   !A
           CALL JIOXXX(zep(1,J),zem(1,J),GZ_ZF,ZMASS,ZWIDTH,ze(1,J))   !Z
           call anomal_formfactor(qzzw(0,J),qze(0,J),qzmu(0,J),qwbos(0,J))
            if (bos.eq.223) then
              call wptozzw_had(v(0,1,J),J,h,wzzw(0,J), NW) ! W+ -> ZZW
              call zzw_wptowz_had(v(0,1,J),v(0,2,J),v(0,5,J),v(0,6,J),J,-ie,wzw1(0,J), NW)
              call zzw_wptowz_anomal(v(0,3,J),v(0,4,J),v(0,5,J),v(0,6,J),J,-iu,wzw2(0,J), NW)       ! no hadronic decay of the second Z
            else if (bos.eq.224) then
              call wmtozzw_had(v(0,1,J),J,h,wzzw(0,J), NW) ! W- -> ZZW
              call zzw_wmtowz_had(v(0,1,J),v(0,2,J),v(0,5,J),v(0,6,J),J,-ie,wzw1(0,J), NW)
              call zzw_wmtowz_anomal(v(0,3,J),v(0,4,J),v(0,5,J),v(0,6,J),J,-iu,wzw2(0,J), NW)       ! no hadronic decay of the second Z
            endif
            call zzw_atozz_had(v(0,1,J),v(0,2,J),v(0,3,J),v(0,4,J),-ie,-iu,azezm(0,J))
            call zzw_ztozz_had(v(0,1,J),v(0,2,J),v(0,3,J),v(0,4,J),-ie,-iu,zzezm(0,J))
        enddo

C*******************  q1 qbar2 ---> Z Z W g   **********************
         
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

      call qq_ZZWj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     &              pdf,res,nmin,nmax, NW)

C*******************  qbar1 q2 ---> Z Z W g   **********************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=1
      physToDiag(3)=3    ! gluon

      fsign(1) = 1
      fsign(2) = -1
      gsign    = 1

      call qq_ZZWj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     &              pdf,res,nmin,nmax, NW)

C*******************  g q1 ---> Z Z W q2  ***********************************

      physToDiag(1)=3
      physToDiag(2)=1
      physToDiag(3)=2

      fsign(1) = 1
      fsign(2) = 1
      gsign    = -1

      call qq_ZZWj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1              pdf,res,nmin,nmax, NW)

C*******************  qbar2 g ---> Z Z W qbar1   **********************

      physToDiag(1)=2
      physToDiag(2)=3
      physToDiag(3)=1

      fsign(1) = -1
      fsign(2) = -1
      gsign    = -1

      call qq_ZZWj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1              pdf,res,nmin,nmax, NW)

C*******************  g q2bar ---> Z Z W q1bar  **********************

      physToDiag(1)=3
      physToDiag(2)=2
      physToDiag(3)=1

      fsign(1) = -1
      fsign(2) = -1
      gsign    = -1

      call qq_ZZWj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1              pdf,res,nmin,nmax, NW)

C*******************  q1 g ---> Z Z W q2   **********************

      physToDiag(1)=1
      physToDiag(2)=3
      physToDiag(3)=2

      fsign(1) = 1
      fsign(2) = 1
      gsign    = -1

      call qq_ZZWj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1              pdf,res,nmin,nmax, NW)


c***  end of extra process evaluation for hadronic decay  ***

       endif
      end select

c*****************  end of process evaluation  *********************

      SELECT CASE(procid)
      CASE(ZZWPhad, ZZWMhad)
       call hadDecayFactor_W(finalquarks(1), abs(mjj2(v(0,5,1),v(0,6,1))), N_gen_W, fac_W)
       do j=0,3
        DO IRES = 1,NMAX
          res(IRES,j) = res(IRES,j) * fac_W            ! factor fac_W for hadronic decay
        ENDDO
       enddo
      CASE(ZZhadWM, ZZhadWP)
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
               CASE(ZZWPhad)
                  if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
                     rnumb=RandomNumber()
                     finalquarks_psp(1) =  2 + 2* INT(rnumb*2)
                     finalquarks_psp(2) = -1 - 2* INT(rnumb*2)
                  endif
               CASE(ZZWMhad)
                  if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
                     rnumb=RandomNumber()
                     finalquarks_psp(1) =  1 + 2* INT(rnumb*2)
                     finalquarks_psp(2) = -2 - 2* INT(rnumb*2)
                  endif
               CASE(ZZhadWP, ZZhadWM)
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

         m2s(j) = m2s(j) * 4.0d0     ! factor for random helicity summation
      enddo

      RETURN
      END

c*****************************************************************************
c
c    end subroutine m2s_qqZZWj
c
c*****************************************************************************

c******************************************************************
c
c   begin subroutine qq_ZZWj
c
c******************************************************************
c
c       Vera Hankele, <vera@particle.uni-karlsruhe.de>
c	Initial version:  2008 February
c	Last modified:  2009 January 28
c
c       qq_ZZWj calls the amplitude square routines 
c             qqZZWj    for qq --> ZZW g --> g 6 leptons
c                           gq --> ZZW q --> q 6 leptons 
c       and   qqZZW for the subtraction terms and finite collinear
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
c          NW              NW = 1: NW approximation, NW = 0: full calc.
c  OUTPUT:
c          udzzw(ID)   output in format as in qqZqq(j), but corresponding
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
c***********************************************************************

      subroutine qq_ZZWj(xi,p,v,physToDiag,fsign,gsign,bos, 
     1                    nlo,lok,xuz,
     1                    pdf,res,nmin,nmax, NW)
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"

      double precision p(0:3,max_p,max_kin), v(0:3,max_v,max_kin), xi(nx), xuz
      double precision pdf(-6:6,2,3), res(maxnumsubproc,0:3)
      integer physToDiag(3), fsign(3+n_v), 
     1        gsign, bos, nlo, nmin, nmax
      logical lok(3)
      integer NW

      double precision C2, TR, N                              ! color factors
      parameter(N=3d0,TR=0.5d0,c2=4d0/3d0)

c Note: Factor 3 from color sum included in qqZZW. I am using color summed
c amplitudes here, while CS use color averages. Thus I later divide by
c 8*3 for an initial qg state and by 3*3 for an initial qq state.

c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/coupl.inc"

c  helicity selection
      INTEGER h
      COMMON /hcount / h

      double precision udzzwb(2,2:3), duzzwb(2,2:3), udzzw, duzzw
      double precision udzzw2, duzzw2, duzzwg, udzzwg
      double precision pbar(0:3,3+n_v), qbar(0:4),q2,v2,sub(2:3), 
     1       dotrr, gs2(2:3), polcolq, polcolg
      double precision NCmatrixelt(0:1,3,2)
      integer iflav(3), diagToPhys(3), FL_VVVg
      external dotrr, FL_VVVg

      double complex zero
      parameter (zero=(0d0,0d0))

      logical ldebug
      save gs2,polcolq,polcolg,udzzwb,duzzwb
      integer i,if1,j,k,mu
      
      parameter (ldebug=.false.)

c***********************************************************************
c
c Begin Code ***********************************************************
c
c***********************************************************************
      udzzw  = 0d0
      duzzw  = 0d0
      udzzwg = 0d0
      duzzwg = 0d0

c  if this is a new phase space point -->
c  reset the LO amplitude to 0 to avoid wrong subtraction
      if (physToDiag(1).eq.1.and.physToDiag(2).eq.2) then
         do i = 1,2
            do j = 2,3
               udzzwb(i,j) = 0d0
               duzzwb(i,j) = 0d0
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
         
c get the real emission amplitude squared, store it in duzzw, udzzw
      do mu = 0,3
         do i = 1,3
            pbar(mu,physToDiag(i))=p(mu,i,1)
         enddo
         qbar(mu) = pbar(mu,3)
      enddo
      qbar(4)=0d0

      do mu = 0,3               ! kinematics for qq -> ZZW -> 6 leptons
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
            call sud_zzwg(pbar,fsign,udzzwg)          !madgraph
            call sdu_zzwg(pbar,fsign,duzzwg)
         endif
#endif

         if (bos.eq.223) then
            call qqZZWj(pbar,fsign,qbar,gsign,udzzw,bos, NW)
         elseif (bos.eq.224) then
            call qqZZWj(pbar,fsign,qbar,gsign,duzzw,bos, NW)
         endif

#ifdef WITH_MADGRAPH
         if (ldebug) then
            if ( abs(udzzw/udzzwg - 1) .gt. 1d-4) then
               print*,"p(0,1,1) = ",p(0,1,1)
               print*,"udzzw  =   ",udzzw
               print*,"udzzwg =  ",udzzwg
               print*,"udzzwg/udzzw - 1 = ",udzzwg/udzzw - 1
               print*," "
            endif
         endif
#endif

      endif

      NCmatrixelt(0,2,1)=udzzw
      NCmatrixelt(1,2,1)=duzzw

c*************************************************************
c
c for the NLO case get the subtraction terms; 
c first the case with a final state gluon (gsign = +1)
c
c*************************************************************
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

c calculation of matrix elements for q qb -> Z Z W g:

               if(lok(j)) then
#ifdef WITH_MADGRAPH
                  if(ldebug) then
                     call sud_zzw(pbar,fsign,udzzw2)  !madgraph
                     call sdu_zzw(pbar,fsign,duzzw2)
                  endif
#endif

                  if (bos.eq.223) then
                     call qqZZW(pbar,fsign,0,J,udzzw,bos, NW)
                  else if (bos.eq.224) then
                     call qqZZW(pbar,fsign,0,J,duzzw,bos, NW)
                  endif
  
                  if (ldebug) then
                        print*,"udzzw/udzzw2 =   ",udzzw/udzzw2
                        print*,"duzzw/duzzw2 =   ",duzzw/duzzw2
                  endif
               else
                  udzzw = 0d0
                  duzzw = 0d0
               endif

               NCmatrixelt(0,j,2)=udzzw   !Born for momentum
               NCmatrixelt(1,j,2)=duzzw   !configuration 2 and 3

c save matrix elements for later use with initial gluons
               udzzwb(physToDiag(1),j) = udzzw
               duzzwb(physToDiag(1),j) = duzzw

               q2 = 2d0*xuz*dotrr(qbar,p(0,j-1,1)) !p(mu,j-1,1) is inc.quark
               v2 = 2d0*gs2(j)*                    ! 8 pi alfas
     &             ( 2d0/(1-xuz)-(1+xuz) )
               sub(j) = v2/q2

            enddo
         endif !nlo.eq.1
#endif /* WITH_NLO */

         iflav(3) = 0           ! final state gluon id
         do if1=1,(nfl/2)*2,2
            if (bos.eq.223) then
               iflav(1) = (if1-physToDiag(1)+2)*fsign(physToDiag(1))
               iflav(2) = (if1+physToDiag(1)-1)*fsign(physToDiag(2))
            else if (bos.eq.224) then
               iflav(1) = (if1+physToDiag(1)-1)*fsign(physToDiag(1))
               iflav(2) = (if1-physToDiag(1)+2)*fsign(physToDiag(2))
            endif

            k=FL_VVVg(iflav,2)

            if (lok(1)) then
               res(k,1)=pdf(sign1*iflav(1),1,1)  
     &                 *pdf(sign2*iflav(2),2,1)  
     &                 *NCmatrixelt(bos-223,2,1)*polcolq
            else
               res(k,1) = 0
            endif

            if ( nlo.eq.1 .and. (lok(2).or.lok(3)) ) then
#ifdef WITH_NLO
               res(k,2) = ( pdf(sign1*iflav(1),1,2)     ! collinear
     &                    *pdf(sign2*iflav(2),2,2)
     &                    -pdf(sign1*iflav(1),1,1)  ! real emission
     &                    *pdf(sign2*iflav(2),2,2 )*sub(2))*C2
     &                    *NCmatrixelt(bos-223,2,2)
               res(k,3) = ( pdf(sign1*iflav(1),1,3)     ! collinear
     &                    *pdf(sign2*iflav(2),2,3)
     &                    -pdf(sign1*iflav(1),1,3)  ! real emission
     &                    *pdf(sign2*iflav(2),2,1)*sub(3))*C2
     &                    *NCmatrixelt(bos-223,3,2)

               res(k,2) = res(k,2) * polcolq
               res(k,3) = res(k,3) * polcolq
               res(k,0) = res(k,1) + res(k,2) + res(k,3)

c****************************************************************************
c debugging for collinear and soft divergences
c****************************************************************************
c       if(ldebug) then                                          !end debug
c          temp = qbar(0) !soft                                  !debug
cc          if(temp.le.1.0d-1) then                              !debug
cc          temp = abs(dotrr(qbar,pbar(0,1))) !collinear         !debug
cc          temp2 = abs(dotrr(qbar,pbar(0,2)))                   !debug
c          if(temp.le.0.0001.or.temp2.le.0.0001) then            !debug
c             temp = abs(res(k,1)/(res(k,2)+res(k,3)))           !debug
c             if(temp.le.0.99d0.or.temp.ge.1.01d0) then          !debug
c                print*," "                                      !debug
c                print*,"gluon in the final state "              !debug
c                print*,"if1 = ",if1                             !debug
c                print*,"qbar(0) = ",qbar(0)                     !debug
c                print*,"qbar*pbar(0,1)=",dotrr(qbar,pbar(0,1))  !debug
c                print*,"qbar*pbar(0,2)=",dotrr(qbar,pbar(0,2))  !debug
c                print*,"res(k,1) = ",res(k,1)                   !debug
c                print*,"res(k,2)+res(k,3)=",res(k,2)+res(k,3)   !debug
c                print*,"res(k,0) = ",res(k,0)                   !debug
c                print*,"|res(k,1)/res(k,2)| = ",temp            !debug
c                print*," "                                      !debug
c             endif                                              !debug
c          endif                                                 !debug
c       endif                                                    !end debug
c****************************************************************************
#endif /* WITH_NLO */

            else
               res(k,0) = res(k,1)
               res(k,2) = 0
               res(k,3) = 0
            endif
         enddo  

c************************************************************
c
c initial gluon section
c
c************************************************************
      elseif (gsign.eq.-1) then

         if (physToDiag(1).eq.3) then
            j=2              ! j=2: emission off 1,2 line
         elseif (physToDiag(2).eq.3) then
            j=3              ! j=3: emission off 3,4 line
         endif

#ifdef WITH_NLO
         if (nlo.eq.1) then          
            q2 = 2d0 * dotrr(p(0,3,1),p(0,j-1,1))*xuz
            v2 = 2d0*gs2(j) * ( (1-xuz)**2 + xuz**2 )
            sub(j) = v2/q2

            if (physToDiag(1).eq.1.or.PhysToDiag(2).eq.2) then
               NCmatrixelt(0,2,2) = udzzwb(1,2)
               NCmatrixelt(1,2,2) = duzzwb(1,2)
               NCmatrixelt(0,3,2) = udzzwb(1,3)
               NCmatrixelt(1,3,2) = duzzwb(1,3)
            else
               NCmatrixelt(0,2,2) = udzzwb(2,2)
               NCmatrixelt(1,2,2) = duzzwb(2,2)
               NCmatrixelt(0,3,2) = udzzwb(2,3)
               NCmatrixelt(1,3,2) = duzzwb(2,3)
            endif
         endif
#endif /* WITH_NLO */

         iflav(j-1) = 0
         do if1=1,(nfl/2)*2,2
            if (bos.eq.223) then
               iflav(4-j)=(if1+physToDiag(3)-1)*fsign(physToDiag(4-j))
               iflav(3)  =(if1-physToDiag(3)+2)*fsign(physToDiag(3))
            else if (bos.eq.224) then
               iflav(4-j)=(if1-physToDiag(3)+2)*fsign(physToDiag(4-j))
               iflav(3)  =(if1+physToDiag(3)-1)*fsign(physToDiag(3))
            endif
            k=FL_VVVg(iflav,j)

            if (lok(1)) then
               res(k,1)=pdf(sign1*iflav(1),1,1)  
     &              *pdf(sign2*iflav(2),2,1)       
     &              *NCmatrixelt(bos-223,2,1)*polcolg
            else
               res(k,1) = 0
            endif

            if ( nlo.eq.1 .and. (lok(2).or.lok(3)) ) then
#ifdef WITH_NLO
                  res(k,j) = ( pdf(sign1*iflav(1),1,j)     ! collinear
     &                    *pdf(sign2*iflav(2),2,j)
     &                    -pdf(sign1*iflav(1),1,2*j-3)  ! real emission
     &                    *pdf(sign2*iflav(2),2,4-j )*sub(j))
     &                    *NCmatrixelt(bos-223,j,2)

                 res(k,j) = res(k,j)*C2*polcolg
                 res(k,-j+5) = 0d0

                 res(k,0) = res(k,1) + res(k,j)

c***************************************************************************
c debugging of soft and collinear divergencies
c***************************************************************************
c       if (ldebug) then                                       !begin debug
cc          temp = p(0,3,1) !soft                              !debug
cc          if(temp.le.0.1d0) then                             !debug
c          temp = abs(dotrr(p(0,3,1),p(0,1,1)))  !collinear    !debug
c          temp2 = abs(dotrr(p(0,3,1),p(0,2,1)))               !debug
c          if(temp.le.0.0001.or.temp2.le.0.0001) then          !debug
c             if(j.eq.3) then                                  !debug
c                temp = abs((res(k,1))+res(k-4,1))/            !debug
c     &                             (res(k,2)+res(k-4,2))      !debug
c                if(temp.le.0.9.or.temp.ge.1.1) then           !debug
c                   print*," "                                 !debug
c                   print*," gluon in the initial state"       !debug
c                   print*," k = ",k," if1 = ",if1             !debug
c                   print*,"p(0,3,1) = ",p(0,3,1)              !debug
c                   print*,"pa*p1= ",dotrr(p(0,3,1),p(0,1,1))  !debug
c                   print*,"pb*p1= ",dotrr(p(0,3,1),p(0,2,1))  !debug
c                   print*,"res(k,1) = ",res(k,1)              !debug
c                   print*,"res(k,2) = ",res(k,2)              !debug
c                endif                                         !debug
c             endif                                            !debug
c          endif                                               !debug
c       endif                                                  !end debug
c***************************************************************************
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
c   end subroutine qq_ZZWj
c
c*****************************************************************
