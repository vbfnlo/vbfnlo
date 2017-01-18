c*********************  subroutine m2s_qqWZj*********************************
c
c INPUT     bos       Boson identifier, 2 = Z, 6=H
c           nlo       NLO=1: create subtraction term; LO = 0
c           lok       need to calculate m2s(1:3)? T or F 
c           xi        Feynman x parameters of incoming fermions
c           p         momenta associated with external fermions
c           v         momenta associated with WW decay fermions
c           rn        random number array
c           xuz       x_ija for subtraction terms
c          
c
c OUTPUT    m2s       |M|^2*pdf1*pdf2 
c
c       Adapted from m2s_VVj by Vera Hankele
c       
c       Johannes Bellm, <johannes@particle.physik.uni-karlsruhe.de>
c	Initial version: Juli 2011
c	Last modified:
c
c***************************************************************************

      subroutine m2s_qqWZj(bos,nlo,lok,xi,p,v,rn,xuz,m2s)

      implicit none

c declare input/output variables

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/process.inc"

      integer bos
      real*8 xi(nx), p(0:3,max_p,max_kin), v(0:3,max_v,max_kin), 
     1       rn(1), xuz, m2s(0:3)
      logical lok(3)

c  helicity selection

      integer jsig, jsig1, jsig3, jsig5
      common /chelsum/ jsig,jsig1,jsig3,jsig5

c  helicity selection
      INTEGER h, it 
      COMMON /hcount / h

      DOUBLE PRECISION zero
      PARAMETER (zero=0d0)

c declare external functions

      integer FL_VVg
      external FL_VVg

c alfas, scales etc

      real*8 x1,x2
      real*8 lnQomu(2:3), omxi(2:3), 
     &       Ax(2:3), Bx(2:3), Cx(2:3), Dxi(2:3), tgs2oqsq(2:3), 
     &       ln1mxi, lnrat, z, lnz, qwpz(0:4,1:3), qwpa(0:4,1:3)

c declare local variables

      real*8 q_sf
      real*8 dotrr
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

      real*8 res(maxnumsubproc,0:3)

c declare parton distribution variables

      real*8 pdf(-6:6,2,3)

c define program switches

      logical ldebug
      parameter (ldebug=.false.)

      real*8 weight,rnumb,RandomNumber

c variables for hadronic decays
      integer N_gen_W
      real*8 fac_W, mjj2
      external mjj2
      integer N_gen_up, N_gen_down
      real*8 fac_Z_up, fac_Z_down

c --------------------------------------------------------------

c if first iteration, output acceptance cuts and fix constant input 
c parameters for call of qqWZj.f

      if ( init .ne. 0) then
         do i = 1,nmax
            do j=0,3
               res(i,j) = 0
            enddo
         enddo
      else if (init.eq.0) then

         if(bos.eq.32) then
            write(6,*) " "
            write(6,*) "W+Z j amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
         else if(bos.eq.42) then
            write(6,*) " "
            write(6,*) "W-Z j amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
         else
            write(6,*) " wrong call of m2s_qqWZ, bos = ",bos
            write(6,*) "-----------------------------------------------"
         endif
         call DIBOSwtowz_reset

         call printnfl(.true.)
         print*," "
         print*," creal = ",crealVV," cvirtual = ",cvirtVV
         init = 1
         do i = 1,maxnumsubproc
            do j = 0,3
               res(i,j) = 0
            enddo
         enddo
      endif

c ------------------------------------------------------------------- 

c for WZ->4l precalculate W->WZ,W->WA  leptonic tensors
c lepton spinors and W+- polarization vectors

c...Les Houches interface

c lepton helicities for Z -> mu- mu+ not fixed
c -> sum over all possible helicities in |M|**2
 
c select helicity: h ... random number for lepton helicity combination (h=1:2)
        h = mod(h,2) + 1
	it = (-1)**(h+1)


        if ((lha.or.hepmc) .and. (.not. doNLO)) then
           helicity(1)=-1
           helicity(2)= 1
           helicity(3)=-it
           helicity(4)= it
        endif

         do J = 1,3
 
           CALL OXXXXX(v(0,1,J),ZERO ,-1,1,wve(1,J)) !W(1,3))     !ve
           CALL IXXXXX(v(0,2,J),ZERO ,1,-1,wep(1,J)) !W(1,4))          !e+   
           CALL OXXXXX(v(0,3,J),ZERO ,-it,1,wmu(1,J)) !W(1,5))          !mu- 
           CALL IXXXXX(v(0,4,J),ZERO ,it,-1,wmubar(1,J)) !W(1,6))          !mu+

c later you need the an: A-> mumu,wp: w+->e+ve and zn: Z->mumu 


           do mu = 0,3
               qp(mu,J) = v(mu,1,J)+v(mu,2,J)   !W+
               qz(mu,J) = v(mu,3,J)+v(mu,4,J)   !Z
               qa(mu,J) = v(mu,3,J)+v(mu,4,J)   !A
               qwpz(mu,J) = qp(mu,J) + qz(mu,J) + p(mu,3,J)
               qwpa(mu,J) = qp(mu,J) + qa(mu,J) + p(mu,3,J)
           enddo
           qp(4,J) = dotrr(qp(0,J),qp(0,J))
           qz(4,J) = dotrr(qz(0,J),qz(0,J))
           qa(4,J) = dotrr(qa(0,J),qa(0,J))
           qwpz(4,J) = dotrr(qwpz(0,J),qwpz(0,J))
           qwpa(4,J) = dotrr(qwpa(0,J),qwpa(0,J))


           SELECT CASE(procid)
           CASE(WPhadZ,WPZhad,WMhadZ,WMZhad)
             ! for general output to all combinations: up-type first
             if ((finalquarks(1).eq.93 .or. finalquarks(1).eq.94) .and.
     &           (procid.eq.WPZhad .or. procid.eq.WMZhad)) then
               ! set couplings for Z hadronic decay into uubar
               call setZtouu
             endif
             ! for hadronic decays anom and SM use the same tensors
             CALL JIOXXX(wep(1,J),wve(1,J),GWF ,WMASS,WWIDTH,wp(1,J)) !!W+
             CALL JIOXXX(wmubar(1,J),wmu(1,J),GZ_ZF ,ZMASS,ZWIDTH,zn(1,J)) !!Z
             CALL JIOXXX(wmubar(1,J),wmu(1,J),GZ_AF ,ZERO,ZERO,an(1,J)) !!A
             if (bos.eq.32) then
                call DIBOSwptowpz_had(v(0,1,J),J,wtowz(0,J))
             else if (bos.eq.42) then
                call DIBOSwmtowmz_had(v(0,1,J),J,wtowz(0,J))
             endif
           CASE DEFAULT
             CALL JIOXXX(wep(1,J),wve(1,J),GWF ,WMASS,WWIDTH,wp(1,J)) !!W+
             CALL JIOXXX(wmubar(1,J),wmu(1,J),GZL ,ZMASS,ZWIDTH,zn(1,J)) !!Z
             CALL JIOXXX(wmubar(1,J),wmu(1,J),GAL ,ZERO,ZERO,an(1,J)) !!A
             if (with_anom) then ! anomalous gauge boson couplings noanomal
                      if (bos.eq.32) then
                            call DIBOSwptowpz_an(v(0,1,J),J,wtowz(0,J))
                      else if (bos.eq.42) then
                            call DIBOSwmtowmz_an(v(0,1,J),J,wtowz(0,J))
                      endif
             else  ! SM
                     if (bos.eq.32) then
                           call DIBOSwptowpz(v(0,1,J),J,wtowz(0,J))
                     else if (bos.eq.42) then
                           call DIBOSwmtowmz(v(0,1,J),J,wtowz(0,J))
                     endif
             endif 
           END SELECT

         enddo

c call PDF subroutine in order to determine parton
c distributions in the incoming (anti)protons.
c for the NLO contributions x1 = x*y in my notes with x=xuz, y=xi(i)

      x1 = xi(1)*xuz
      q_sf = sqrt(mufsq(1,1))
      call pdfproton( xi(1), q_sf, pdf(-6,1,1) )      ! f_a(y)=f_a(x1/z)
      if(nlo.ne.0)then
      q_sf = sqrt(mufsq(1,2))                        ! f_a(x1) for upper line 
      call pdfproton( x1, q_sf, pdf(-6,1,2) )         !   NLO correction

      if (mufsq(1,3).ne.mufsq(1,1)) then              ! f_a(x1) for lower line
         q_sf = sqrt(mufsq(1,3))                      !   NLO correction
         call pdfproton( xi(1), q_sf, pdf(-6,1,3) )
      else
         do i = -6,6
            pdf(i,1,3) = pdf(i,1,1)
         enddo
      endif
      endif


      x2 = xi(2)*xuz
      q_sf = sqrt(mufsq(2,1))
      call pdfproton( xi(2), q_sf, pdf(-6,2,1) )      ! f_b(y)=f_a(x2/z)
      if(nlo.ne.0)then
      q_sf = sqrt(mufsq(2,3))                        ! f_b(x2) for lower line 
      call pdfproton( x2, q_sf, pdf(-6,2,3) )         !   NLO correction

      if (mufsq(2,2).ne.mufsq(2,1)) then              ! f_b(x2) for upper line
         q_sf = sqrt(mufsq(2,2))                      !   NLO correction
         call pdfproton( xi(2), q_sf, pdf(-6,2,2) )
      else
         do i = -6,6
            pdf(i,2,2) = pdf(i,2,1)
         enddo
      endif
      endif


c and fill the coefficient functions for the finite subtraction terms
C Note that the color factors TR and C2 are NOT!!! included here

      lnQomu(2) = log(qwpz(4,2)/mufsq(1,2))
      lnQomu(3) = log(qwpz(4,3)/mufsq(2,3))

      omxi(2) = 1d0-x1
      omxi(3) = 1d0-x2

c      factor 4*pi*alfas/(2x pa pb ) = 4*pi*alfas/Q^2

      tgs2oqsq(2) = 4d0*pi*als(1,2)/qwpz(4,2)
      tgs2oqsq(3) = 4d0*pi*als(2,3)/qwpz(4,3)

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
     1                  + 2d0*ln1mxi**2) + crealVV
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


C********************************************************************

C     Now sum over the subprocesses contributing to W Z production

      nmax = FL_VVg(fsign,-1)          !reset counter for subprocesses to 0

C*******************  q1 qbar2 ---> W Z g   **********************

c   physToDiag(ext.momentum label) = Feynman diagram label

      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=2
      physToDiag(3)=3    ! gluon

      fsign(1) = 1          !incoming quark
      fsign(2) = -1         !incoming antiquark
      fsign(3) = 1          !outgoing neutrino
      fsign(4) = -1         !outgoing lepton
      fsign(5) = 1          !outgoing antilepton
      fsign(6) = -1         !outgoing antineutrino
      gsign    = 1          !outgoing gluon

      call qq_WZj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     &              pdf,res,nmin,nmax)

C*******************  qbar1 q2 ---> W Z g   **********************
 
      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=1
      physToDiag(3)=3    ! gluon

      fsign(1) = 1
      fsign(2) = -1

      call qq_WZj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     &              pdf,res,nmin,nmax)

C*******************  g q1 ---> W Z q2  ***********************************

      physToDiag(1)=3
      physToDiag(2)=1
      physToDiag(3)=2

      fsign(1) = 1
      fsign(2) = 1
      gsign    = -1

      call qq_WZj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1              pdf,res,nmin,nmax)

C*******************  qbar2 g ---> W Z qbar1   **********************

      physToDiag(1)=2
      physToDiag(2)=3
      physToDiag(3)=1

      fsign(1) = -1
      fsign(2) = -1
      gsign = -1

      call qq_WZj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1              pdf,res,nmin,nmax)

C*******************  g q2bar ---> W Z q1bar  **********************

      physToDiag(1)=3
      physToDiag(2)=2
      physToDiag(3)=1

      fsign(1) = -1
      fsign(2) = -1
      gsign = -1

      call qq_WZj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1              pdf,res,nmin,nmax)

C*******************  q1 g ---> W Z q2   **********************

      physToDiag(1)=1
      physToDiag(2)=3
      physToDiag(3)=2

      fsign(1) = 1
      fsign(2) = 1
      gsign    = -1 

      call qq_WZj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1              pdf,res,nmin,nmax)


c**************  end of process evaluation part one ******************


! now run again for down-type decay products if 93 93 / 94 94 is 
! requested and the z decays hadronically.

      SELECT CASE(procid)
      CASE(WPZhad, WMZhad)
       if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
        ! set couplings for Z hadronic decay into ddbar
        call setZtodd
        do J = 1,3
          CALL JIOXXX(wmubar(1,J),wmu(1,J),GZ_ZF ,ZMASS,ZWIDTH,zn(1,J)) !!Z
          CALL JIOXXX(wmubar(1,J),wmu(1,J),GZ_AF ,ZERO,ZERO,an(1,J)) !!A
          if (bos.eq.32) then
             call DIBOSwptowpz_had(v(0,1,J),J,wtowz(0,J))
          else if (bos.eq.42) then
             call DIBOSwmtowmz_had(v(0,1,J),J,wtowz(0,J))
          endif
        enddo

C*******************  q1 qbar2 ---> W Z g   **********************

c   physToDiag(ext.momentum label) = Feynman diagram label

      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=2
      physToDiag(3)=3    ! gluon

      fsign(1) = 1          !incoming quark
      fsign(2) = -1         !incoming antiquark
      fsign(3) = 1          !outgoing neutrino
      fsign(4) = -1         !outgoing lepton
      fsign(5) = 1          !outgoing antilepton
      fsign(6) = -1         !outgoing antineutrino
      gsign    = 1          !outgoing gluon

      call qq_WZj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     &              pdf,res,nmin,nmax)

C*******************  qbar1 q2 ---> W Z g   **********************
 
      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=1
      physToDiag(3)=3    ! gluon

      fsign(1) = 1
      fsign(2) = -1

      call qq_WZj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     &              pdf,res,nmin,nmax)

C*******************  g q1 ---> W Z q2  ***********************************

      physToDiag(1)=3
      physToDiag(2)=1
      physToDiag(3)=2

      fsign(1) = 1
      fsign(2) = 1
      gsign    = -1

      call qq_WZj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1              pdf,res,nmin,nmax)

C*******************  qbar2 g ---> W Z qbar1   **********************

      physToDiag(1)=2
      physToDiag(2)=3
      physToDiag(3)=1

      fsign(1) = -1
      fsign(2) = -1
      gsign = -1

      call qq_WZj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1              pdf,res,nmin,nmax)

C*******************  g q2bar ---> W Z q1bar  **********************

      physToDiag(1)=3
      physToDiag(2)=2
      physToDiag(3)=1

      fsign(1) = -1
      fsign(2) = -1
      gsign = -1

      call qq_WZj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1              pdf,res,nmin,nmax)

C*******************  q1 g ---> W Z q2   **********************

      physToDiag(1)=1
      physToDiag(2)=3
      physToDiag(3)=2

      fsign(1) = 1
      fsign(2) = 1
      gsign    = -1 

      call qq_WZj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1              pdf,res,nmin,nmax)


c***  end of extra process evaluation for hadronic decay  ***

       endif
      end select


c*****************  end of process evaluation  *********************


      SELECT CASE(procid)
      CASE(WPhadZ,WMhadZ)
       call hadDecayFactor_W(finalquarks(1), abs(mjj2(v(0,1,1),v(0,2,1))), N_gen_W, fac_W)
       do j=0,3
        DO IRES = 1,NMAX
          res(IRES,j) = res(IRES,j) * fac_W            ! factor fac_W for hadronic decay
        ENDDO
       enddo
      CASE(WPZhad,WMZhad)
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


      if (init.eq.1) then
         init = init+1
      endif
      do j=0,3
         m2s(j) = 0
         DO IRES = 1,NMAX

            m2s(j) = m2s(j) + RES(IRES,j) 
           
         ENDDO
  
         m2s(j)= m2s(j)*2 ! factor 2 for random helicity summation

         if(j.eq.0)then
c...Les Houches interface - the most propable subprocess 3jets at LO  
            if ((lha.or.hepmc) .and. (.not. doNLO)) then
               i=0
               weight=0.d0
               rnumb=RandomNumber()
               do while((i.le.nmax).and.(weight.le.rnumb*m2s(0)))
                  i=i+1
                  weight=weight+res(i,0)
                  iprocess=i
               enddo
               SELECT CASE(procid)
               CASE(WPhadZ)
                  if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
                     rnumb=RandomNumber()
                     finalquarks_psp(1) =  2 + 2* INT(rnumb*2)
                     finalquarks_psp(2) = -1 - 2* INT(rnumb*2)
                  endif
               CASE(WMhadZ)
                  if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
                     rnumb=RandomNumber()
                     finalquarks_psp(1) =  1 + 2* INT(rnumb*2)
                     finalquarks_psp(2) = -2 - 2* INT(rnumb*2)
                  endif
               CASE(WPZhad,WMZhad)
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

      enddo

      RETURN
      END

c*****************************************************************************
c
c    end subroutine m2s_qqVVj
c
c*****************************************************************************
c******************************************************************
c
c   begin subroutine qq_WPZj
c
c*****************************************************************

      subroutine qq_WZj(xi,p,v,physToDiag,fsign,gsign,bos, 
     1                    nlo,lok,xuz,
     1                    pdf,res,nmin,nmax)
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"

      real*8 p(0:3,max_p,max_kin), v(0:3,max_v,max_kin)
      real*8 xi(nx), xuz
      real*8 pdf(-6:6,2,3), res(maxnumsubproc,0:3)
      integer physToDiag(3), fsign(3+n_v), 
     1        gsign, bos, nlo, nmin, nmax
      logical lok(3)
c
c qq_WPZj calls the amplitude square routines 
c             qqWZj    for qq --> WZ g --> g 4 leptons
c                          gq --> WZ q --> q 4 leptons 
c
c for the subtraction terms for the NLO cross section calculation
c
c  INPUT:  p(0:3,5,3)      external physical parton momenta
c          v(0:3,nv,3)       Z decay momenta, for Higgs production only
c                          the sum q(mu) = v(mu,1)+...+v(mu,nv) is needed
c          physToDiag(3)   physToDiag(ext.mom. label) = Feynman diagram label
c          fsign,gsign     sign factors for amplitude calls; see qqZqq(j)
c          nlo             nlo = 0: calculate LO only (i.e. no subtraction)
c                          nlo = 1: calculate full NLO subtraction
c          lok(3)          lok(ID)=T means momenta set ID passes acceptance
c                          cuts, i.e res(k,ID) needs to be calculated
c  OUTPUT:
c          uucc(ID)   output in format as in qqZqq(j), but corresponding
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

c Note: Factor 3 from color sum included in qqVV. I am using color summed
c amplitudes here, while CS use color averages. Thus I later divide by
c 8*3 for an initial qg state and by 3*3 for an initial qq state.

c alfas, scales etc

#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/coupl.inc"

c
c  helicity selection

      integer jsig, jsig1, jsig3, jsig5
      common /chelsum/ jsig,jsig1,jsig3,jsig5
c
      real*8 udwzb(2,2:3),
     &       udwz, udwzg
      real*8 pbar(0:3,3+n_v), qbar(0:4),q2,v2,sub(2:3), 
     1       e_in(2), dotrr, gs2(2:3), polcolq, polcolg

      double precision NCmatrixelt(3,2), temp, temp2, temp3, temp1
      integer iflav(3), diagToPhys(3), FL_VVg
      external dotrr, FL_VVg

      complex*16 zero
      parameter (zero=(0d0,0d0))


      logical ldebug, lres
      data e_in/0d0,0d0/
      save e_in,gs2,polcolq,polcolg, udwzb
      integer i,if1,j,k,mu
      
      parameter (ldebug=.false.)
      PARAMETER (lres=.false.)

c--------------------------------------------------------------------

      udwz  = 0d0
      udwzg    = 0d0

      if (lres) then
        open(unit=31,file="irfinite_WZ.chk",ACCESS='APPEND')
      endif

c  if this is a new phase space point -->
c  reset the LO amplitude to 0 to avoid wrong subtraction
      if (physToDiag(1).eq.1.and.physToDiag(2).eq.2) then
         do i = 1,2
            do j = 2,3
               udwzb(i,j) = 0d0
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
         
c get the real emission amplitude squared, store it in udwz, ddww
      do mu = 0,3
         do i = 1,3
            pbar(mu,physToDiag(i))=p(mu,i,1)
         enddo
         qbar(mu) = pbar(mu,3)
      enddo
      qbar(4)=0d0

      do mu = 0,3            ! kinematics for WZ -->4 lepton decay
         pbar(mu,3) = v(mu,1,1) ! ve
         pbar(mu,4) = v(mu,2,1) ! e+
         pbar(mu,5) = v(mu,3,1) ! mu-
         pbar(mu,6) = v(mu,4,1) ! nu_mu~
         pbar(mu,7) = qbar(mu)
      enddo
      fsign(7) = gsign

c calculate the matrix elements pa + pb -> ve e+ mu- vm~ p1
      if (lok(1) .or. nlo.eq.0) then

            call qqWZj(pbar,fsign,qbar,gsign,udwz,bos)

c       No debugging here, every subprocess is debugged

      
      endif


      NCmatrixelt(2,1)=udwz
       udwz = 0d0

c for the NLO case get the subtraction terms; 
c first the case with a final state gluon (gsign = +1)

      if (gsign.eq.1) then
!          if (ldebug) then
!             print*," final state gluon section in qq_WPZj "
!             print*," jsig = ",jsig," jsig1,3 = ",jsig1,jsig3
!             print 101," fsign = ",fsign
!  101        format(a,6i5,a,i5)
!          endif

#ifdef WITH_NLO
         if (nlo.eq.1) then
            do j = 2,3               ! j=2: emission off 1,2 line
               do mu = 0,3           ! j=3: emission off 3,4 line
                  do i = 1,2
                     pbar(mu,physToDiag(i))=p(mu,i,j)
                  enddo
                  do i = 1,4
                     pbar(mu,i+2)=v(mu,i,j)
                  enddo
               enddo

c calculation of matrix elemnets for q qb -> W Z g:

               if(lok(j)) then
! #ifdef WITH_MADGRAPH No debugging here because sudwpz(...)is tested in m2s_qqWZ.F with J=1
!                   if(ldebug) then
!                      call sudwpz(pbar,fsign,udwz2)  !madgraph
! 
!                   endif
! #endif
! 
                  call qqWZ(pbar,fsign,0,J,udwz,bos)


               else
                  udwz = 0d0
               endif

               NCmatrixelt(j,2)=udwz   !Born for momentum
                                          !configuration 2 and 3


c save matrix elements for later use with initial gluons
               udwzb(physToDiag(1),j) = udwz

               q2 = 2d0*xuz*dotrr(qbar,p(0,j-1,1)) !p(mu,j-1,1) is inc.quark
               v2 = 2d0*gs2(j)*                   ! 8 pi alfas
     &             ( 2d0/(1-xuz)-(1+xuz) )
               sub(j) = v2/q2

            enddo
         endif !nlo.eq.1
#endif


         iflav(3) = 0           ! final state gluon id
         do if1=2,(nfl/2)*2,2

             if (bos.eq.32) then
                 iflav(1)=( if1 - physToDiag(1) + 1 )*fsign(physToDiag(1))
                 iflav(2)=( if1 - physToDiag(2) + 1 )*fsign(physToDiag(2))
             else if (bos.eq.42) then
                 iflav(1)=( if1 + physToDiag(1) - 2 )*fsign(physToDiag(1))
                 iflav(2)=( if1 + physToDiag(2) - 2 )*fsign(physToDiag(2))
             endif


               k=FL_VVg(iflav,2)
               if(lok(1)) then
                  res(k,1)=pdf(sign1*iflav(1),1,1)  !1 for j=2;3 for j=3
     &                 *pdf(sign2*iflav(2),2,1)  !2 for j=2;1 for j=3
     &                 *NCmatrixelt(2,1)*polcolq
               else
                     res(k,1) = 0
               endif


               if ( (nlo.eq.1) .and. (lok(2).or.lok(3)) ) then
#ifdef WITH_NLO
                  res(k,2) = ( pdf(sign1*iflav(1),1,2)     ! collinear
     &                       *pdf(sign2*iflav(2),2,2)
     &                       -pdf(sign1*iflav(1),1,1)  ! real emission
     &                       *pdf(sign2*iflav(2),2,2 )*sub(2))*C2
     &                       *NCmatrixelt(2,2)

                  res(k,3) = ( pdf(sign1*iflav(1),1,3)     ! collinear
     &                       *pdf(sign2*iflav(2),2,3)
     &                       -pdf(sign1*iflav(1),1,3)  ! real emission
     &                       *pdf(sign2*iflav(2),2,1)*sub(3))*C2
     &                       *NCmatrixelt(3,2)


                  res(k,2) = res(k,2) * polcolq
                  res(k,3) = res(k,3) * polcolq
                  res(k,0) = res(k,1) + res(k,2) + res(k,3)


c debugging for collinear and soft divergences -----------------------

               if (lres) then
                 temp1 = sqrt(abs(dotrr(qbar,pbar(0,1)))) ! collinear 1
                 temp2 = sqrt(abs(dotrr(qbar,pbar(0,2)))) ! collinear 2
                 temp3 = abs(res(k,1)/(res(k,2)+res(k,3)))
                 if ( (temp3.gt.0.) .and. (qbar(0).le.0.2 .or. min(temp1,temp2).le.0.2) )
     &             write(31,*) "0", 
     &               qbar(0),
     &               temp1,
     &               temp2,
     &               min(temp1,temp2),
     &               abs(res(k,1)/(res(k,2)+res(k,3)))
               endif

                  if(ldebug) then                              !begin debug
                     temp = qbar(0)              !soft
                     if(temp.le.1.0d0) then    
!                      print*, lok(1),lok(2),lok(3)
c                     temp = abs(dotrr(qbar,pbar(0,1)))  !collinear
c                     temp2 = abs(dotrr(qbar,pbar(0,2)))
c                     if(temp.le.0.0001.or.temp2.le.0.0001) then
                        temp = abs(res(k,1)/(res(k,2)+res(k,3)))
                        if(temp.le.0.9d0.or.temp.ge.1.1d0) then
                           print*," "
                           print*,"gluon in the final state "
                           print*,"if1 = ",if1
                           print*,"qbar(0) = ",qbar(0)
                           print*,"qbar*pbar(0,1)=",dotrr(qbar,pbar(0,1))
                           print*,"qbar*pbar(0,2)=",dotrr(qbar,pbar(0,2))
                           print*,"res(k,1) = ",res(k,1)
                           print*,"res(k,2) = ",res(k,2)
                           print*,"res(k,3) = ",res(k,3)
                           print*,"res(k,0) = ",res(k,0)
                           print*,"res(k,2)+res(k,3)=",res(k,2)+res(k,3)
                           print*,"|res(k,1)/(res(k,2)+res(k,3))| = ",temp
                           print*," "
                        endif
                     endif
                  endif                                       !end debug
c----------------------------------------------------------------------
#endif

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
               NCmatrixelt(2,2)=udwzb(1,2)
               NCmatrixelt(3,2)=udwzb(1,3)
            else
               NCmatrixelt(2,2)=udwzb(2,2)
               NCmatrixelt(3,2)=udwzb(2,3)
            endif
#endif
         endif !nlo


            iflav(j-1) = 0

         do if1=2,(nfl/2)*2,2


            if (bos.eq.32) then
               iflav(4-j)=( if1 - physToDiag(4-j) + 1 )*fsign(physToDiag(4-j))
               iflav(3)  =( if1 - physToDiag(3) + 1 )*fsign(physToDiag(3))
            else if (bos.eq.42) then
               iflav(4-j)=( if1 + physToDiag(4-j) - 2 )*fsign(physToDiag(4-j))
               iflav(3)  =( if1 + physToDiag(3) - 2 )*fsign(physToDiag(3))
            endif


            k=FL_VVg(iflav,j)


            if (lok(1)) then

               res(k,1)=pdf(sign1*iflav(1),1,1)  !1 for j=2;3 for j=3
     &              *pdf(sign2*iflav(2),2,1)      !2 for j=2;1 for j=3
     &              *NCmatrixelt(2,1)*polcolg       
            else
               res(k,1) = 0
            endif

c--------------------------------------------------------------

            if ( (nlo.eq.1) .and. (lok(2).or.lok(3)) ) then
#ifdef WITH_NLO

               res(k,j) = ( pdf(sign1*iflav(1),1,j)     ! collinear
     &                       *pdf(sign2*iflav(2),2,j)
     &                       -pdf(sign1*iflav(1),1,2*j-3)  ! real emission
     &                       *pdf(sign2*iflav(2),2,4-j )*sub(j))
     &                       *NCmatrixelt(j,2)

              res(k,j) = res(k,j)*C2*polcolg
              res(k,-j+5) = 0d0

              res(k,0) = res(k,1) + res(k,j)

               if (lres) then
                 temp1 = sqrt(abs(dotrr(p(0,j-1,1),p(0,3,1)))) ! collinear 1
                 temp2 = sqrt(abs(dotrr(p(0,j-1,1),p(0,4-j,1)))) ! collinear 2
                 temp3 = abs(res(k,1)/res(k,j))
                 if ( (temp3.gt.0.) .and. (min(temp1,temp2).le.0.2) )
     &             write(31,*) "1", 
     &               p(0,j-1,1),
     &               temp1,
     &               temp2,
     &               min(temp1,temp2),
     &               temp3
               endif

!               if (ldebug) then                         !begin debug
!                  temp = p(0,3,1)                 !soft
!                  if(temp.le.0.001d0) then
! c                     temp = abs(dotrr(p(0,3,1),p(0,1,1)))  !collinear
! c                     temp2 = abs(dotrr(p(0,3,1),p(0,2,1)))
! c                     if(temp.le.0.0001.or.temp2.le.0.0001) then
! c                     if(j.eq.3) then
! c                        temp = abs((res(k,1))+res(k-4,1))/
! c     &                             (res(k,2)+res(k-4,2))
! c                        if(temp.le.0.9.or.temp.ge.1.1) then
! c                           print*," "
! c                           print*," gluon in the initial state"
! c                           print*," k = ",k," if1 = ",if1      
!                         print*,"p(0,3,1) = ",p(0,3,1)
! c                           print*,"pa*p1= ",dotrr(p(0,3,1),p(0,1,1))
! c                           print*,"pb*p1= ",dotrr(p(0,3,1),p(0,2,1))
!                         print*,"res(k,1) = ",res(k,1)
!                         print*,"res(k,",j,") = ",res(k,j)
! c                     print*,"res(k,1)+res(k-4,1) = ",res(k,1)+res(k-4,1)
! c                     print*,"res(k,2)+res(k-4,2) = ",res(k,2)+res(k-4,2)
! c                           endif
! c                        endif
!                  endif
!               endif                                   !end debug
#endif

            else
               res(k,0) = res(k,1)
               res(k,2) = 0
               res(k,3) = 0
            endif

         enddo
      endif

      nmax = FL_VVg(iflav,-2)

      if (lres) then
        close(31)
      endif

      end

c******************************************************************
c
c   end subroutine qq_WPZj
c
c*****************************************************************



