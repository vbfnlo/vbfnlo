*****************************************************************************
c
c    begin subroutine m2s_qqWWZj
c
c*****************************************************************************
c
c       Vera Hankele, <vera@particle.uni-karlsruhe.de>
c	Initial version:  2007 January 23
c	Last modified:  2009 January 28
c
c       this subroutine computes |ME|^2 for qq -> ve e+ mu-vm~ ta- ta+ g
c       and all the crossed processes including the gluon induced
c       subprocesses at LO and the virtual corrections via call of "qqWWZj.f"
c
c	For the decay Z -> ta- ta+ random helicity summation is employed
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

      SUBROUTINE m2s_qqWWZj(bos,nlo,lok,xi,p,v,rn,xuz,m2s, NW)

      IMPLICIT NONE

c declare input/output variables
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/process.inc"

      INTEGER bos
      DOUBLE PRECISION xi(nx), p(0:3,max_p,max_kin), v(0:3,max_v,max_kin), 
     1       rn(1), xuz, m2s(0:3)
      LOGICAL lok(3)
      INTEGER NW    !in:  = 1:  NW Approximation, = 0: Full calculation

c  helicity selection
      INTEGER h
      COMMON /hcount / h

      COMPLEX*16 zero
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

c variables for hadronic decays
      integer N_gen_W
      real*8 fac_W, mjj2
      external mjj2
      integer N_gen_up, N_gen_down
      real*8 fac_Z_up, fac_Z_down

c---------------------------- BEGIN CODE --------------------------------

      if ( init .ne. 0) then
         do i = 1,nmax
            do j=0,3
               res(i,j) = 0
            enddo
         enddo
      else if (init.eq.0) then
         if (bos.eq.342 .or. bos.eq.432) then
            write(6,*) " "
            write(6,*) "W+ W- Z j amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
            call vtowwz_reset
         endif

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
            qp(mu,J) = v(mu,1,J)+v(mu,2,J)   !W+
            qm(mu,J) = v(mu,3,J)+v(mu,4,J)   !W-
            qz(mu,J) = v(mu,5,J)+v(mu,6,J)   !Z
            qww(mu,J) = qp(mu,J) + qm(mu,J)
            qwwz(mu,J) = qp(mu,J) + qm(mu,J) + qz(mu,J)! + p(mu,3,J)
         enddo

         qp(4,J) = dotrr(qp(0,J),qp(0,J))
         qm(4,J) = dotrr(qm(0,J),qm(0,J))
         qz(4,J) = dotrr(qz(0,J),qz(0,J))
         qww(4,J) = dotrr(qww(0,J),qww(0,J))
         qwwz(4,J) = dotrr(qwwz(0,J),qwwz(0,J))          !s        
      enddo !J = 1,3
 
c for WWZ->6l precalculate A->WWZ,Z->WWZ etc..

c lepton spinors and W+, W-, Z polarisation vectors
c lepton helicities for Z -> ta- ta+ not fixed
c -> sum over all possible helicities in |M|**2
 
c select helicity: h ... random number for lepton helicity
c                        combination (h=1:2)

       h = mod(h,2) + 1
       it = (-1)**(h+1)

c...Les Houches interface
        if ((lha.or.hepmc).and..not.doNLO) then
           helicity(1)=-1
           helicity(2)= 1
           helicity(3)=-1
           helicity(4)= 1
           helicity(5)=-it
           helicity(6)= it
        endif

      do J = 1,3

         CALL OXXXXX(v(0,1,J),ZERO ,-1,1,wve(1,J))   !ve 
         CALL IXXXXX(v(0,2,J),ZERO ,1,-1,wep(1,J))   !e+        
         CALL OXXXXX(v(0,3,J),ZERO ,-1,1,wmu(1,J))   !mu-      
         CALL IXXXXX(v(0,4,J),ZERO ,1,-1,wvm(1,J))   !vm~
         CALL OXXXXX(v(0,5,J),ZERO ,-it,1,ztm(1,J)) !ta-      
         CALL IXXXXX(v(0,6,J),ZERO ,it,-1,ztp(1,J)) !ta+

         CALL JIOXXX(wep(1,J),wve(1,J),GWF,WMASS,WWIDTH, wp(1,J))  !W+
         CALL JIOXXX(wvm(1,J),wmu(1,J),GWF,WMASS,WWIDTH, wm(1,J))  !W-
         CALL JIOXXX(ztp(1,J),ztm(1,J),GAL,ZERO ,ZERO  ,atau(1,J)) !A
         CALL JIOXXX(ztp(1,J),ztm(1,J),GZL,ZMASS,ZWIDTH,ztau(1,J)) !Z

c leptonic tensors, last entry defines whether the non resonant diagrams
c contribute ( = 0: all diagrams, = 1: no non resonant diagrams)

           SELECT CASE(procid)
           CASE(WPhadWMZ, WPWMhadZ, WWZhad)
             ! for general output to all combinations: up-type first
             if ((finalquarks(1).eq.93 .or. finalquarks(1).eq.94) .and.
     &           (procid.eq.WWZhad)) then
               ! set couplings for Z hadronic decay into uubar
               call setZtouu
             endif
             call anomal_formfactor(qwwz(0,J),qp(0,J),qm(0,J),qz(0,J))
             CALL JIOXXX(ztp(1,J),ztm(1,J),GZ_AF,ZERO ,ZERO  ,atau(1,J)) !A
             CALL JIOXXX(ztp(1,J),ztm(1,J),GZ_ZF,ZMASS,ZWIDTH,ztau(1,J)) !Z
             call TRIBOSatoww_had(v(0,1,J),J,aww(0,J),NW) !A -> W+ W- -> ve e+ mu- vm~
             call TRIBOSztoww_had(v(0,1,J),J,zww(0,J),NW) !Z -> W+ W- -> ve e+ mu- vm~
             call wptowz_had(v(0,1,J),J,h,wpwz(0,J),NW) !W+ -> W+ Z -> ve e+ ta- ta+
             call wmtowz_had(v(0,1,J),J,h,wmwz(0,J),NW) !W- -> W- Z -> mu- vm~ ta- ta+
             call atowwz_had(v(0,1,J),J,h,awwz(0,J),NW)!A -> WWZ ->ve e+ mu- vm~ ta- ta+
             call ztowwz_had(v(0,1,J),J,h,zwwz(0,J),NW)!Z -> WWZ ->ve e+ mu- vm~ ta- ta+
           CASE DEFAULT
               if (with_anom) then ! anomalous gauge boson couplings
c                 using global form factor for all tensors of one phase space point
c                 this ensures proper cancellations for anomalous contributions
c                 energy scale is invariant WWZ mass
                  call anomal_formfactor(qwwz(0,J),qp(0,J),qm(0,J),qz(0,J))
                  call TRIBOSatoww_anomal(v(0,1,J),J,aww(0,J),NW) !A -> W+ W- -> ve e+ mu- vm~
                  call TRIBOSztoww_anomal(v(0,1,J),J,zww(0,J),NW) !Z -> W+ W- -> ve e+ mu- vm~
                  call wptowz_anomal(v(0,1,J),J,h,wpwz(0,J),NW) !W+ -> W+ Z -> ve e+ ta- ta+
                  call wmtowz_anomal(v(0,1,J),J,h,wmwz(0,J),NW) !W- -> W- Z -> mu- vm~ ta- ta+
                  call atowwz_anomal(v(0,1,J),J,h,awwz(0,J),NW)!A -> WWZ ->ve e+ mu- vm~ ta- ta+
                  call ztowwz_anomal(v(0,1,J),J,h,zwwz(0,J),NW)!Z -> WWZ ->ve e+ mu- vm~ ta- ta+
               else if (with_kk) then ! Kaluza Klein Scenario
#ifdef WITH_KK
                  call TRIBOSatoww(v(0,1,J),J,aww(0,J),NW) !A -> W+ W- -> ve e+ mu- vm~
                  call TRIBOSztoww(v(0,1,J),J,zww(0,J),NW) !Z -> W+ W- -> ve e+ mu- vm~
                  call wptowz(v(0,1,J),J,h,wpwz(0,J),NW) !W+ -> W+ Z -> ve e+ ta- ta+
                  call wmtowz(v(0,1,J),J,h,wmwz(0,J),NW) !W- -> W- Z -> mu- vm~ ta- ta+
                  call atowwz_kk(v(0,1,J),J,h,awwz(0,J),NW) !A -> WWZ -> ve e+ mu- vm~ ta- ta+
                  call ztowwz_kk(v(0,1,J),J,h,zwwz(0,J),NW) !Z -> WWZ -> ve e+ mu- vm~ ta- ta+
#endif
               else ! SM scenario
                  call TRIBOSatoww(v(0,1,J),J,aww(0,J),NW) !A -> W+ W- -> ve e+ mu- vm~
                  call TRIBOSztoww(v(0,1,J),J,zww(0,J),NW) !Z -> W+ W- -> ve e+ mu- vm~
                  call wptowz(v(0,1,J),J,h,wpwz(0,J),NW) !W+ -> W+ Z -> ve e+ ta- ta+
                  call wmtowz(v(0,1,J),J,h,wmwz(0,J),NW) !W- -> W- Z -> mu- vm~ ta- ta+
                  call atowwz(v(0,1,J),J,h,awwz(0,J),NW) !A -> WWZ -> ve e+ mu- vm~ ta- ta+
                  call ztowwz(v(0,1,J),J,h,zwwz(0,J),NW) !Z -> WWZ -> ve e+ mu- vm~ ta- ta+
               endif
           END SELECT

      enddo ! J = 1,3

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

      lnQomu(2) = log(qwwz(4,2)/mufsq(1,2))
      lnQomu(3) = log(qwwz(4,3)/mufsq(2,3))
 
      omxi(2) = 1d0-x1
      omxi(3) = 1d0-x2

c      factor 4*pi*alfas/(2x pa pb ) = 4*pi*alfas/Q^2

      tgs2oqsq(2) = 4d0*pi*als(1,2)/qwwz(4,2)
      tgs2oqsq(3) = 4d0*pi*als(2,3)/qwwz(4,3)

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

C     Now sum over the subprocesses contributing to W+ W- Z production

      nmax = FL_VVVg(fsign,-1)          !reset counter for subprocesses to 0

C*******************  q1 qbar2 ---> W+ W- Z g   **********************
         
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

      call qq_WWZj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     &              pdf,res,nmin,nmax, NW)

C*******************  qbar1 q2 ---> W+ W- Z g   **********************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=1
      physToDiag(3)=3    ! gluon

      fsign(1) = 1
      fsign(2) = -1

      call qq_WWZj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     &              pdf,res,nmin,nmax, NW)

C*******************  g q1 ---> W+ W- Z q2  ***********************************

      physToDiag(1)=3
      physToDiag(2)=1
      physToDiag(3)=2

      fsign(1) = 1
      fsign(2) = 1
      gsign    = -1

      call qq_WWZj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1              pdf,res,nmin,nmax, NW)

C*******************  qbar2 g ---> W+ W- Z qbar1   **********************

      physToDiag(1)=2
      physToDiag(2)=3
      physToDiag(3)=1

      fsign(1) = -1
      fsign(2) = -1
      gsign = -1

      call qq_WWZj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1              pdf,res,nmin,nmax, NW)

C*******************  g q2bar ---> W+ W- Z q1bar  **********************

      physToDiag(1)=3
      physToDiag(2)=2
      physToDiag(3)=1

      fsign(1) = -1
      fsign(2) = -1
      gsign = -1

      call qq_WWZj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1              pdf,res,nmin,nmax, NW)

C*******************  q1 g ---> W+ W- Z q2   **********************

      physToDiag(1)=1
      physToDiag(2)=3
      physToDiag(3)=2

      fsign(1) = 1
      fsign(2) = 1

      call qq_WWZj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1              pdf,res,nmin,nmax, NW)


c**************  end of process evaluation part one ******************


! now run again for down-type decay products if 93 93 / 94 94 is 
! requested and the z decays hadronically.

      SELECT CASE(procid)
      CASE(WWZhad)
       if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
        ! set couplings for Z hadronic decay into ddbar
        call setZtodd
        do J = 1,3
          call anomal_formfactor(qwwz(0,J),qp(0,J),qm(0,J),qz(0,J))
          CALL JIOXXX(ztp(1,J),ztm(1,J),GZ_AF,ZERO ,ZERO  ,atau(1,J)) !A
          CALL JIOXXX(ztp(1,J),ztm(1,J),GZ_ZF,ZMASS,ZWIDTH,ztau(1,J)) !Z
          call TRIBOSatoww_had(v(0,1,J),J,aww(0,J),NW) !A -> W+ W- -> ve e+ mu- vm~
          call TRIBOSztoww_had(v(0,1,J),J,zww(0,J),NW) !Z -> W+ W- -> ve e+ mu- vm~
          call wptowz_had(v(0,1,J),J,h,wpwz(0,J),NW) !W+ -> W+ Z -> ve e+ ta- ta+
          call wmtowz_had(v(0,1,J),J,h,wmwz(0,J),NW) !W- -> W- Z -> mu- vm~ ta- ta+
          call atowwz_had(v(0,1,J),J,h,awwz(0,J),NW)!A -> WWZ ->ve e+ mu- vm~ ta- ta+
          call ztowwz_had(v(0,1,J),J,h,zwwz(0,J),NW)!Z -> WWZ ->ve e+ mu- vm~ ta- ta+
        enddo


C*******************  q1 qbar2 ---> W+ W- Z g   **********************
         
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

      call qq_WWZj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     &              pdf,res,nmin,nmax, NW)

C*******************  qbar1 q2 ---> W+ W- Z g   **********************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=1
      physToDiag(3)=3    ! gluon

      fsign(1) = 1
      fsign(2) = -1

      call qq_WWZj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     &              pdf,res,nmin,nmax, NW)

C*******************  g q1 ---> W+ W- Z q2  ***********************************

      physToDiag(1)=3
      physToDiag(2)=1
      physToDiag(3)=2

      fsign(1) = 1
      fsign(2) = 1
      gsign    = -1

      call qq_WWZj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1              pdf,res,nmin,nmax, NW)

C*******************  qbar2 g ---> W+ W- Z qbar1   **********************

      physToDiag(1)=2
      physToDiag(2)=3
      physToDiag(3)=1

      fsign(1) = -1
      fsign(2) = -1
      gsign = -1

      call qq_WWZj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1              pdf,res,nmin,nmax, NW)

C*******************  g q2bar ---> W+ W- Z q1bar  **********************

      physToDiag(1)=3
      physToDiag(2)=2
      physToDiag(3)=1

      fsign(1) = -1
      fsign(2) = -1
      gsign = -1

      call qq_WWZj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1              pdf,res,nmin,nmax, NW)

C*******************  q1 g ---> W+ W- Z q2   **********************

      physToDiag(1)=1
      physToDiag(2)=3
      physToDiag(3)=2

      fsign(1) = 1
      fsign(2) = 1

      call qq_WWZj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1              pdf,res,nmin,nmax, NW)


c***  end of extra process evaluation for hadronic decay  ***

       endif
      end select


c*****************  end of process evaluation  *********************


      SELECT CASE(procid)
      CASE(WPhadWMZ)
       call hadDecayFactor_W(finalquarks(1), abs(mjj2(v(0,1,1),v(0,2,1))), N_gen_W, fac_W)
       do j=0,3
        DO IRES = 1,NMAX
          res(IRES,j) = res(IRES,j) * fac_W            ! factor fac_W for hadronic decay
        ENDDO
       enddo
      CASE(WPWMhadZ)
       call hadDecayFactor_W(finalquarks(1), abs(mjj2(v(0,3,1),v(0,4,1))), N_gen_W, fac_W)
       do j=0,3
        DO IRES = 1,NMAX
          res(IRES,j) = res(IRES,j) * fac_W            ! factor fac_W for hadronic decay
        ENDDO
       enddo
      CASE(WWZhad)
       call hadDecayFactor_Z(finalquarks(1), abs(mjj2(v(0,5,1),v(0,6,1))), N_gen_up, N_gen_down, fac_Z_up, fac_Z_down)
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
               CASE(WPhadWMZ)
                  if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
                     rnumb=RandomNumber()
                     finalquarks_psp(1) =  2 + 2* INT(rnumb*2)
                     finalquarks_psp(2) = -1 - 2* INT(rnumb*2)
                  endif
               CASE(WPWMhadZ)
                  if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
                     rnumb=RandomNumber()
                     finalquarks_psp(1) =  1 + 2* INT(rnumb*2)
                     finalquarks_psp(2) = -2 - 2* INT(rnumb*2)
                  endif
               CASE(WWZhad)
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

         m2s(j) = m2s(j) * 2.0d0     ! factor 2 for random helicity summation
      enddo

      RETURN
      END

c*****************************************************************************
c
c    end subroutine m2s_qqWWZj
c
c*****************************************************************************
c******************************************************************
c
c   begin subroutine qq_WWZj
c
c*****************************************************************
c
c       Vera Hankele, <vera@particle.uni-karlsruhe.de>
c	Initial version:  2006 March 23
c	Last modified:  2009 January 28
c
c       qq_WWZj calls the amplitude square routines 
c             qqWWZj    for qq --> WWZ g --> g 6 leptons
c                           gq --> WWZ q --> q 6 leptons 
c       and   qqWWZ for the subtraction terms and finite collinear
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
c          uuWWZ(ID)   output in format as in qqZqq(j), but corresponding
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
c*****************************************************************************

      SUBROUTINE qq_WWZj(xi,p,v,physToDiag,fsign,gsign,bos, 
     1                    nlo,lok,xuz,
     1                    pdf,res,nmin,nmax, NW)
      IMPLICIT NONE

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"

      DOUBLE PRECISION p(0:3,max_p,max_kin), v(0:3,max_v,max_kin), xi(nx), xuz
      DOUBLE PRECISION pdf(-6:6,2,3), res(maxnumsubproc,0:3)
      INTEGER physToDiag(3), fsign(3+n_v), 
     1        gsign, bos, nlo, nmin, nmax
      LOGICAL lok(3)
      INTEGER NW

c
      DOUBLE PRECISION C2, TR, N                              ! color factors
      PARAMETER(N=3d0,TR=0.5d0,c2=4d0/3d0)

c Note: Factor 3 from color sum included in qqWWZ. I am using color summed
c amplitudes here, while CS use color averages. Thus I later divide by
c 8*3 for an initial qg state and by 3*3 for an initial qq state.

c alfas, scales etc

#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/coupl.inc"

c
c  helicity selection
      INTEGER h
      COMMON /hcount / h

      DOUBLE PRECISION uuwwzb(2,2:3), ddwwzb(2,2:3), uuwwz, ddwwz
      DOUBLE PRECISION uuwwz2, ddwwz2, uuwwzg, ddwwzg
      DOUBLE PRECISION pbar(0:3,3+n_v), qbar(0:4),q2,v2,sub(2:3), 
     1       dotrr, gs2(2:3), polcolq, polcolg
      DOUBLE PRECISION NCmatrixelt(0:1,3,2)
      INTEGER iflav(3), diagToPhys(3), FL_VVVg
      EXTERNAL dotrr, FL_VVVg

      COMPLEX*16 zero
      PARAMETER (zero=(0d0,0d0))

      LOGICAL ldebug
      SAVE gs2,polcolq,polcolg,uuwwzb,ddwwzb
      INTEGER i,if1,j,k,mu
      
      PARAMETER (ldebug=.false.)

c--------------------------------------------------------------------

      uuwwz  = 0d0
      ddwwz  = 0d0
      uuwwzg = 0d0
      ddwwzg = 0d0

c  if this is a new phase space point -->
c  reset the LO amplitude to 0 to avoid wrong subtraction
      if (physToDiag(1).eq.1.and.physToDiag(2).eq.2) then
         do i = 1,2
            do j = 2,3
               uuwwzb(i,j) = 0d0
               ddwwzb(i,j) = 0d0
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
         
c get the real emission amplitude squared, store it in uuwwz, ddwwz
      do mu = 0,3
         do i = 1,3
            pbar(mu,physToDiag(i))=p(mu,i,1)
         enddo
         qbar(mu) = pbar(mu,3)
      enddo
      qbar(4)=0d0

      do mu = 0,3               ! kinematics for qq -> WWZ -> 6 leptons
         pbar(mu,3) = v(mu,1,1) ! ve
         pbar(mu,4) = v(mu,2,1) ! e+
         pbar(mu,5) = v(mu,3,1) ! mu-
         pbar(mu,6) = v(mu,4,1) ! nu_mu~
         pbar(mu,7) = v(mu,5,1) ! ta-
         pbar(mu,8) = v(mu,6,1) ! ta+
         pbar(mu,9) = qbar(mu)
      enddo
      fsign(9) = gsign

c calculate the matrix elements pa + pb -> ve e+ mu- vm~ ta- ta+ p1
      if (lok(1) .or. nlo.eq.0) then

#ifdef WITH_MADGRAPH
         if (ldebug) then
            call suu_wwzg(pbar,fsign,uuwwzg,nw)          !madgraph
            call sdd_wwzg(pbar,fsign,ddwwzg,nw)
         endif
#endif

         call qqWWZj(pbar,fsign,qbar,gsign,uuwwz,ddwwz, NW)

#ifdef WITH_MADGRAPH
         if (ldebug) then
            if ( (abs(uuwwz/uuwwzg - 1) .gt. 1d-4) .or.
     1           (abs(ddwwz/ddwwzg - 1) .gt. 1d-4) ) then
               print*,"p(0,1,1) = ",p(0,1,1)
               print*,"uuwwz  =   ",uuwwz
               print*,"uuwwzg =  ",uuwwzg
               print*,"uuwwz/uuwwzg - 1 = ",uuwwz/uuwwzg - 1
               print*,"ddwwz =   ",ddwwz
               print*,"ddwwzg =  ",ddwwzg
               print*,"ddwwz/ddwwzg - 1 = ",ddwwz/ddwwzg - 1
               print*," "
            endif
         endif
#endif

      endif

      NCmatrixelt(0,2,1)=uuwwz
      NCmatrixelt(1,2,1)=ddwwz

c for the NLO case get the subtraction terms; 
c first the case with a final state gluon (gsign = +1)
c************************************************************
c
c final gluon section
c
c************************************************************

      if (gsign.eq.1) then
         if (ldebug) then
            print*," final state gluon section in qq_WWZj "
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

c calculation of matrix elements for q qb -> W W Z g:

               if(lok(j)) then
#ifdef WITH_MADGRAPH
                  if(ldebug) then
                     call suu_wwz(pbar,fsign,uuwwz2,nw)  !madgraph
                     call sdd_wwz(pbar,fsign,ddwwz2,nw)
                  endif
#endif

                  call QQWWZ(pbar,fsign,0,J,uuwwz,ddwwz, NW)

#ifdef WITH_MADGRAPH
                  if (ldebug) then
                        print*,"uuwwz/uuwwz2 =   ",uuwwz/uuwwz2
                        print*,"ddwwz/ddwwz2 =   ",ddwwz/ddwwz2
                  endif
#endif
               else
                  uuwwz = 0d0
                  ddwwz = 0d0
               endif

               NCmatrixelt(0,j,2)=uuwwz   !Born for momentum
               NCmatrixelt(1,j,2)=ddwwz   !configuration 2 and 3

c save matrix elements for later use with initial gluons
               uuwwzb(physToDiag(1),j) = uuwwz
               ddwwzb(physToDiag(1),j) = ddwwz

               q2 = 2d0*xuz*dotrr(qbar,p(0,j-1,1)) !p(mu,j-1,1) is inc.quark
               v2 = 2d0*gs2(j)*                   ! 8 pi alfas
     &             ( 2d0/(1-xuz)-(1+xuz) )
               sub(j) = v2/q2

            enddo !j = 2,3
         endif !nlo.eq.1
#endif /* WITH_NLO */

         iflav(3) = 0           ! final state gluon id
         do if1=1,(nfl/2)*2
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
     &                       *pdf(sign2*iflav(2),2,2)
     &                       -pdf(sign1*iflav(1),1,1)  ! real emission
     &                       *pdf(sign2*iflav(2),2,2 )*sub(2))*C2
     &                       *NCmatrixelt(mod(if1,2),2,2)
                  res(k,3) = ( pdf(sign1*iflav(1),1,3)     ! collinear
     &                       *pdf(sign2*iflav(2),2,3)
     &                       -pdf(sign1*iflav(1),1,3)  ! real emission
     &                       *pdf(sign2*iflav(2),2,1)*sub(3))*C2
     &                       *NCmatrixelt(mod(if1,2),3,2)

                  res(k,2) = res(k,2) * polcolq
                  res(k,3) = res(k,3) * polcolq
                  res(k,0) = res(k,1) + res(k,2) + res(k,3)

c debugging for collinear and soft divergences -----------------------
c      if(ldebug) then                                            !begin debug
c         temp = qbar(0) !soft                                    !debug
c         if(temp.le.1.0d-1) then                                 !debug
c            temp = abs(dotrr(qbar,pbar(0,1)))  !collinear        !debug
c            temp2 = abs(dotrr(qbar,pbar(0,2)))                   !debug
c            if(temp.le.1.0d-9.or.temp2.le.1.0d-9) then           !debug
c               temp = abs(res(k,1)/(res(k,2)+res(k,3)))          !debug
c               if(temp.le.0.9d0.or.temp.ge.1.1d0) then           !debug
c                  print*," "                                     !debug
c                  print*,"gluon in the final state "             !debug
c                  print*,"if1 = ",if1                            !debug
c                  print*,"qbar(0) = ",qbar(0)                    !debug
c                  print*,"qbar*pbar(0,1)=",dotrr(qbar,pbar(0,1)) !debug
c                  print*,"qbar*pbar(0,2)=",dotrr(qbar,pbar(0,2)) !debug
c                  print*,"res(k,1) = ",res(k,1)                  !debug 
c                  print*,"res(k,2)+res(k,3)=",res(k,2)+res(k,3)  !debug
c                  print*,"res(k,0) = ",res(k,0)                  !debug
c                  print*,"|res(k,1)/res(k,2)| = ",temp           !debug
c                  print*," "                                     !debug
c               endif                                             !debug
c            endif                                                !debug
c         endif                                                   !end debug
c----------------------------------------------------------------------
#endif /* WITH_NLO */
               else
                     res(k,0) = res(k,1)
                     res(k,2) = 0
                     res(k,3) = 0
               endif !nlo
         enddo !if1  

c************************************************************
c
c initial gluon section
c
c************************************************************
      elseif (gsign.eq.-1) then
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
               NCmatrixelt(0,2,2)=uuwwzb(1,2)
               NCmatrixelt(1,2,2)=ddwwzb(1,2)
               NCmatrixelt(0,3,2)=uuwwzb(1,3)
               NCmatrixelt(1,3,2)=ddwwzb(1,3)
            else
               NCmatrixelt(0,2,2)=uuwwzb(2,2)
               NCmatrixelt(1,2,2)=ddwwzb(2,2)
               NCmatrixelt(0,3,2)=uuwwzb(2,3)
               NCmatrixelt(1,3,2)=ddwwzb(2,3)
            endif
#endif /* WITH_NLO */
         endif

         do if1=1,(nfl/2)*2
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
     &                    *pdf(sign2*iflav(2),2,j)
     &                    -pdf(sign1*iflav(1),1,2*j-3)  ! real emission
     &                    *pdf(sign2*iflav(2),2,4-j )*sub(j))
     &                    *NCmatrixelt(mod(if1,2),j,2)

              res(k,j) = res(k,j)*C2*polcolg
              res(k,-j+5) = 0d0

              res(k,0) = res(k,1) + res(k,j)

c debugging for collinear and soft divergences -----------------------
c        if (ldebug) then                                       !begin debug
c           temp = p(0,3,1) !soft                               !debug
cc           if(temp.le.1d-6) then                              !debug
cc           temp = abs(dotrr(p(0,3,1),p(0,1,1))) !collinear    !debug
cc           temp2 = abs(dotrr(p(0,3,1),p(0,2,1)))              !debug
c           if(temp.le.0.0001.or.temp2.le.0.0001) then          !debug
c              if(j.eq.3) then                                  !debug
c                 temp = abs((res(k,1))+res(k-4,1)/             !debug
c     &                             (res(k,2)+res(k-4,2)))      !debug
c                 temp = abs(res(k,1)/res(k,j))                 !debug
c                 if(temp.le.0.99.or.temp.ge.1.01) then         !debug
c                    print*," "                                 !debug
c                    print*," gluon in the initial state"       !debug
c                    print*," k = ",k," if1 = ",if1             !debug
c                    print*,"p(0,3,1) = ",p(0,3,1)              !debug
c                    print*,"res(k,1) = ",res(k,1)              !debug
c                    print*,"res(k,j) = ",res(k,j)              !debug
c                    print*,"res(k,-j+5) = ",res(k,-j+5)        !debug
c                    print*,"res(k,1) = ",res(k,1)              !debug
c                    print*,"res(k,j) = ",res(k,j)              !debug
c                    print*,"res(0)/res(j) = ",temp             !debug
c                 endif                                         !debug
c              endif                                            !debug
c           endif                                               !debug
c        endif                                                  !end debug
c----------------------------------------------------------------------
#endif /* WITH_NLO */

            else
               res(k,0) = res(k,1)
               res(k,2) = 0
               res(k,3) = 0
            endif !nlo .and. lok

         enddo !if1
      endif !gsign

      nmax = FL_VVVg(iflav,-2)

      end

c******************************************************************
c
c   end subroutine qq_WWZj
c
c*****************************************************************

c******************************************************************
c
c   begin function FL_VVVg
c
c*****************************************************************
c
c       Vera Hankele, <vera@particle.uni-karlsruhe.de>
c	Initial version:  2006 March 23
c	Last modified:  2009 January 28
c       
c  fill hepup color and flavor assignments for WBF processes with
c  one attached external gluon and count subprocesses. 
c  There are 2 color structures depending into which fermion line the 
c  gluon is inserted
c
c  colstruc = 2:upper
c             3:lower fermion line
c
c  cases can be distinguished according to which iflav(i)=0, i.e. corresponds
c  to the gluon

c  id1,2 are flavor identifiers for incoming quarks
c  id3,4 are flavor identifiers for outgoing quarks
c  gsign=+1: outgoing gluon, gsign=-1 for incoming gluon
c  colstruc =2,3 determines whether the gluon is coupled to the 
c  uppper or lower fermion line:

c  Note that colstruc is also a flag:
c  if colstruc=-1, then we are resetting.
c
c********************************************************************

      INTEGER FUNCTION FL_VVVg(iflav,colstruc)

      IMPLICIT NONE

      integer iflav(3),colstruc    ! input for color and flavor assignment

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
      
      integer lkup
      common /localkup/ lkup(numParticles,7,maxNumSubProc)
 
      integer listposition
      integer numresets
      save listposition
      save numresets

      data numresets /-1/
      data listposition /0/

      if (colstruc.eq.-1) then  ! we are resetting
c       write(*,*) "we are resetting"
        listposition=0
        numresets=numresets+1
      else if (colstruc.eq.-2) then
c       don"t do anything, just return the number of subprocesses.
      else
c increment the counter regardless of whether or not it"s the 1st time through.
        listposition=listposition+1

c       we fill info for this subprocess,if required
        if(numresets.eq.0 .and. (lha.or.hepmc)) then
          call fillColoredPartons_VVVg(iflav,
     1                                listposition+numdecay,colstruc)
        endif       ! numresets.eq.0
      endif         
      FL_VVVg=listposition
      end
c******************************************************************************
c
c   end function FL_VVVg
c
c******************************************************************************

c*****************************************************************************
c
c    begin  subroutine fillColoredPartons_VVVg
c
c*****************************************************************************
c
c       Vera Hankele, <vera@particle.uni-karlsruhe.de>
c	Initial version:  2006 March 23
c	Last modified:  2009 January 28
c
c  assigns values to the variables in the common block localHEPUP
c  in particular, this subroutine assigns values to those variables that 
c  will be stored in the lookup tables generated by writeHEPUPtable.  
c  As the name suggests, this routine only stores the information for the 
c  colored partons.  Particles without color will be dealt with in the 
c  subroutine fillColorless.
c
c*****************************************************************************

      SUBROUTINE fillColoredPartons_VVVg(iflav,listposition,colstruc)

      IMPLICIT NONE

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/process.inc"
      
      integer iflav(3),id1,id2,id3,listposition,i,i1
      integer iflavour(3),colstruc

      logical ldebug
      parameter (ldebug=.false.)


      do i1=1,3
         iflavour(i1)=iflav(i1) 
         if(iflavour(i1).eq.0) iflavour(i1)=21   
      enddo

c...flavours of 2 quarks and 1 gluon
      id1=iflavour(1)
      id2=iflavour(2)
      id3=iflavour(3)
      
      select case (process)
      case(WWZjLO, ZZZjLO, WPhadWMZjLO, WPWMhadZjLO, WWZhadjLO, ZZZhadjLO)
         lnup(listposition)=numParticles+4 
      case(ZZWPjLO, ZZhadWPjLO, ZZWPhadjLO)
         lnup(listposition)=numParticles+4 
      case(ZZWMjLO, ZZhadWMjLO, ZZWMhadjLO)
         lnup(listposition)=numParticles+4
      case(WWWPjLO,WWWMjLO,WPhadWMWPjLO,WPWMhadWPjLO,WMhadWPWMjLO,WMWPhadWMjLO)
         lnup(listposition)=numParticles+4
      case (WWAjLO, ZZAjLO, WPZAjLO, WMZAjLO,
     &      WPhadWMAjLO, WPWMhadAjLO, ZZhadAjLO, WPhadZAjLO, WPZhadAjLO, WMhadZAjLO, WMZhadAjLO)
         lnup(listposition) = numParticles+2
      case (WPAAjLO, WMAAjLO, ZAAjLO, ZnAAjLO)
         lnup(listposition) = numParticles
      case (AAAjLO)
         lnup(listposition) = numParticles-2
      end select
      
      listup(1,listposition)=-1 !incoming partons
      listup(2,listposition)=-1
      listup(3,listposition)=1  !outgoing partons
 
c...four quarks and gluon
      lidup(1,listposition)=id1    
      lidup(2,listposition)=id2
      lidup(3,listposition)=id3

      do i=3,n_p
         lmothup(1,i,listposition)=1
         lmothup(2,i,listposition)=2
      enddo

c...final state gluon
      if(id3.eq.21)then
         if(id1.gt.0) then
               licolup(1,1,listposition)=501
               licolup(2,1,listposition)=0
               licolup(1,2,listposition)=0
               licolup(2,2,listposition)=502
               licolup(1,3,listposition)=501
               licolup(2,3,listposition)=502
         else
               licolup(1,1,listposition)=0
               licolup(2,1,listposition)=501
               licolup(1,2,listposition)=502
               licolup(2,3,listposition)=0
               licolup(1,3,listposition)=502
               licolup(2,3,listposition)=501               
         endif            
      endif

c...initial particles gluon and quark/antiquark 
      if(id1.eq.21)then 
         if(id2.gt.0)then
            licolup(1,1,listposition)=501
            licolup(2,1,listposition)=502
            licolup(1,2,listposition)=502
            licolup(2,2,listposition)=0
            licolup(1,3,listposition)=501
            licolup(2,3,listposition)=0
         else
            licolup(1,1,listposition)=502
            licolup(2,1,listposition)=501
            licolup(1,2,listposition)=0
            licolup(2,2,listposition)=502
            licolup(1,3,listposition)=0
            licolup(2,3,listposition)=501
         endif 
      endif
      
c...initial particles quark/antiquark and gluon 
      if(id2.eq.21)then 
         if(id1.gt.0)then
            licolup(1,1,listposition)=501
            licolup(2,1,listposition)=0
            licolup(1,2,listposition)=502
            licolup(2,2,listposition)=501
            licolup(1,3,listposition)=502
            licolup(2,3,listposition)=0
         else
            licolup(1,1,listposition)=0
            licolup(2,1,listposition)=501
            licolup(1,2,listposition)=501
            licolup(2,2,listposition)=502
            licolup(1,3,listposition)=0
            licolup(2,3,listposition)=502
         endif 
      endif

      end  ! fillColoredPartons_VVVg
c*****************************************************************************
c
c    end subroutine fillColoredPartons_VVVg
c
c*****************************************************************************
