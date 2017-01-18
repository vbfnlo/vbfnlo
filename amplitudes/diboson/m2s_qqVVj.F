c*********************  subroutine m2s_qqVVj*********************************
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
c       Adapted from m2s_wbfzhg by D. Zeppenfeld
c
c       Vera Hankele, <vera@particle.physik.uni-karlsruhe.de>
c	Initial version: April 2006
c	Last modified:
c
c***************************************************************************

      subroutine m2s_qqVVj(bos,nlo,lok,xi,p,v,rn,xuz,m2s,ps_number)

      implicit none

c declare input/output variables

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/process.inc"

      integer bos, ps_number
      real*8 xi(nx), p(0:3,max_p,max_kin), v(0:3,max_v,max_kin), 
     1       rn(1), xuz, m2s(0:3)
      logical lok(3)

c  helicity selection

      integer jsig, jsig1, jsig3, jsig5
      common /chelsum/ jsig,jsig1,jsig3,jsig5

      DOUBLE PRECISION zero
      PARAMETER (zero=0d0)

c declare external functions

      integer FL_VVg
      external FL_VVg

c alfas, scales etc

      real*8 x1,x2
      real*8 lnQomu(2:3), omxi(2:3), 
     &       Ax(2:3), Bx(2:3), Cx(2:3), Dxi(2:3), tgs2oqsq(2:3), 
     &       ln1mxi, lnrat, z, lnz

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

      integer N_gen_W
      real*8 fac_W, mjj2
      external mjj2

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

c --------------------------------------------------------------

c if first iteration, output acceptance cuts and fix constant input 
c parameters for call of qqVVj.f

      if ( init .ne. 0) then
         res(1:nmax,0:3) = 0
      else if (init.eq.0) then
         if (bos.eq.34 .or. bos.eq.43) then
            write(6,*) " "
            write(6,*) "W+W- j amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
            call dibostoww_reset
         endif
         call printnfl(.true.)
         print*," "
         print*," creal = ",crealVV," cvirtual = ",cvirtVV
         init = 1
         res(1:maxnumsubproc,0:3) = 0
      endif

c ------------------------------------------------------------------- 

c for WW->4l precalculate A->WW,Z->WW,AZ->WW etc leptonic tensors
c lepton spinors and W+- polarization vectors

c...Les Houches interface
        if ((lha.or.hepmc).and..not.doNLO) then
           helicity(1)=-1
           helicity(2)= 1
           helicity(3)=-1
           helicity(4)= 1
        endif

        do J = 1,3

         CALL IXXXXX(v(0,2,J),ZERO ,1,-1,wep(1,J)) !W(1,5))          !e+       
         CALL OXXXXX(v(0,1,J),ZERO ,-1,1,wve(1,J)) !W(1,4))          !ve 
         CALL OXXXXX(v(0,3,J),ZERO ,-1,1,wmu(1,J)) !W(1,6))          !mu-      
         CALL IXXXXX(v(0,4,J),ZERO ,1,-1,wvm(1,J)) !W(1,7))          !vm~
         CALL JIOXXX(wep(1,J),wve(1,J),GWF ,WMASS,WWIDTH,wp(1,J)) !W(1,8))!W+
         CALL JIOXXX(wvm(1,J),wmu(1,J),GWF ,WMASS,WWIDTH,wm(1,J)) !W(1,9))!W-

         do mu = 0,3
            qp(mu,J) = v(mu,1,J)+v(mu,2,J)
            qm(mu,J) = v(mu,3,J)+v(mu,4,J)
            qww(mu,J) = qp(mu,J) + qm(mu,J) + p(mu,3,J)
         enddo
         qp(4,J) = dotrr(qp(0,J),qp(0,J))
         qm(4,J) = dotrr(qm(0,J),qm(0,J))
         qww(4,J) = dotrr(qww(0,J),qww(0,J))          !s

c leptonic tensors

         SELECT CASE(process)
         CASE(WPhadWMlep,WPhadWMlepjLO,WPlepWMhad,WPlepWMhadjLO)
            ! for hadronic decays anom and SM use the same tensors
            call anomal_formfactor(qww(0,J),qp(0,J),qm(0,J),qm(0,J))       ! last parameter is "dummy variable"
            call DIBOSatoww_had(v(0,1,J),J,aww(0,J))
            call DIBOSztoww_had(v(0,1,J),J,zww(0,J))
         CASE DEFAULT
            if (with_anom) then
               call anomal_formfactor(qww(0,J),qp(0,J),qm(0,J),qm(0,J))       ! last parameter is "dummy variable"
               call DIBOSatoww_anom(v(0,1,J),J,aww(0,J))
               call DIBOSztoww_anom(v(0,1,J),J,zww(0,J))
            else     ! standard model
               call DIBOSatoww(v(0,1,J),J,aww(0,J))
               call DIBOSztoww(v(0,1,J),J,zww(0,J))
            end if
         END SELECT

        enddo

c call PDF subroutine in order to determine parton
c distributions in the incoming (anti)protons.
c for the NLO contributions x1 = x*y in my notes with x=xuz, y=xi(i)

      x1 = xi(1)*xuz
      q_sf = sqrt(mufsq(1,1))
      call pdfproton( xi(1), q_sf, pdf(-6,1,1) )      ! f_a(y)=f_a(x1/z)

      if (nlo.gt.0) then
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
      else
         do i = -6,6
            pdf(i,1,3) = pdf(i,1,1)
         enddo
      endif

      x2 = xi(2)*xuz
      q_sf = sqrt(mufsq(2,1))
      call pdfproton( xi(2), q_sf, pdf(-6,2,1) )      ! f_b(y)=f_a(x2/z)

      if (nlo.gt.0) then
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
      else
         do i = -6,6
            pdf(i,2,2) = pdf(i,2,1)
         enddo
      endif
      

      if (nlo.gt.0) then
c and fill the coefficient functions for the finite subtraction terms
C Note that the color factors TR and C2 are NOT!!! included here

      lnQomu(2) = log(qww(4,2)/mufsq(1,2))
      lnQomu(3) = log(qww(4,3)/mufsq(2,3))
 
      omxi(2) = 1d0-x1
      omxi(3) = 1d0-x2

c      factor 4*pi*alfas/(2x pa pb ) = 4*pi*alfas/Q^2

      tgs2oqsq(2) = 4d0*pi*als(1,2)/qww(4,2)
      tgs2oqsq(3) = 4d0*pi*als(2,3)/qww(4,3)

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
      endif

C********************************************************************

C     Now sum over the subprocesses contributing to W+ W- production

      nmax = FL_VVg(fsign,-1)          !reset counter for subprocesses to 0

C*******************  q1 qbar2 ---> V V g   **********************
         
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

      call qq_VVj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     &              pdf,res,nmin,nmax,ps_number)

C*******************  qbar1 q2 ---> V V g   **********************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=1
      physToDiag(3)=3    ! gluon

      fsign(1) = 1
      fsign(2) = -1

      call qq_VVj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     &              pdf,res,nmin,nmax,ps_number)

C*******************  g q1 ---> V V q2  ***********************************

      physToDiag(1)=3
      physToDiag(2)=1
      physToDiag(3)=2

      fsign(1) = 1
      fsign(2) = 1
      gsign    = -1

      call qq_VVj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1              pdf,res,nmin,nmax,ps_number)

C*******************  qbar2 g ---> V V qbar1   **********************

      physToDiag(1)=2
      physToDiag(2)=3
      physToDiag(3)=1

      fsign(1) = -1
      fsign(2) = -1
      gsign = -1

      call qq_VVj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1              pdf,res,nmin,nmax,ps_number)

C*******************  g q2bar ---> V V q1bar  **********************

      physToDiag(1)=3
      physToDiag(2)=2
      physToDiag(3)=1

      fsign(1) = -1
      fsign(2) = -1
      gsign = -1

      call qq_VVj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1              pdf,res,nmin,nmax,ps_number)

C*******************  q1 g ---> V V q2   **********************

      physToDiag(1)=1
      physToDiag(2)=3
      physToDiag(3)=2

      fsign(1) = 1
      fsign(2) = 1

      call qq_VVj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1              pdf,res,nmin,nmax,ps_number)

c*****************  end of process evaluation  *********************

      if (init.eq.1) then
         init = init+1
      endif

      SELECT CASE(procid)
      CASE(WPhadWMlep)
       call hadDecayFactor_W(finalquarks(1), abs(mjj2(v(0,1,1),v(0,2,1))), N_gen_W, fac_W)
       do j=0,3
        DO IRES = 1,NMAX
          res(IRES,j) = res(IRES,j) * fac_W            ! factor fac_W for hadronic decay
        ENDDO
       enddo
      CASE(WPlepWMhad)
       call hadDecayFactor_W(finalquarks(1), abs(mjj2(v(0,3,1),v(0,4,1))), N_gen_W, fac_W)
       do j=0,3
        DO IRES = 1,NMAX
          res(IRES,j) = res(IRES,j) * fac_W            ! factor fac_W for hadronic decay
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
               if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
                  rnumb=RandomNumber()
                  SELECT CASE(process)
                  CASE(WPhadWMlepjLO)
                     finalquarks_psp(1) =  2 + 2* INT(rnumb*2)
                     finalquarks_psp(2) = -1 - 2* INT(rnumb*2)
                  CASE(WPlepWMhadjLO)
                     finalquarks_psp(1) =  1 + 2* INT(rnumb*2)
                     finalquarks_psp(2) = -2 - 2* INT(rnumb*2)
                  END SELECT
               endif
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
c   begin subroutine qq_VVj
c
c*****************************************************************

      subroutine qq_VVj(xi,p,v,physToDiag,fsign,gsign,bos, 
     1                    nlo,lok,xuz,
     1                    pdf,res,nmin,nmax,ps_number)
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"

      real*8 p(0:3,max_p,max_kin), v(0:3,max_v,max_kin)
      real*8 xi(nx), xuz
      real*8 pdf(-6:6,2,3), res(maxnumsubproc,0:3)
      integer physToDiag(3), fsign(3+n_v), 
     1        gsign, bos, nlo, nmin, nmax, ps_number
      logical lok(3)
c
c qq_VVj calls the amplitude square routines 
c             qqVVj    for qq --> VV g --> g 4 leptons
c                          gq --> VV q --> q 4 leptons 
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
      real*8 uuwwb(2,2:3), ddwwb(2,2:3),
     &       uuww, ddww, uuwwg, ddwwg
      real*8 uuww2, ddww2
      real*8 pbar(0:3,3+n_v), qbar(0:4),q2,v2,sub(2:3), 
     1       e_in(2), dotrr, gs2(2:3), polcolq, polcolg
      double precision NCmatrixelt(0:1,3,2), temp
      integer iflav(3), diagToPhys(3), FL_VVg
      external dotrr, FL_VVg

      complex*16 zero
      parameter (zero=(0d0,0d0))

      logical ldebug
      data e_in/0d0,0d0/
      save e_in, gs2,polcolq,polcolg, uuwwb,ddwwb
      integer i,if1,j,k,mu
      
      parameter (ldebug=.false.)

c--------------------------------------------------------------------

      uuww  = 0d0
      ddww  = 0d0
      uuwwg    = 0d0
      ddwwg    = 0d0

c  if this is a new phase space point -->
c  reset the LO amplitude to 0 to avoid wrong subtraction
      if (physToDiag(1).eq.1.and.physToDiag(2).eq.2) then
         do i = 1,2
            do j = 2,3
               uuwwb(i,j) = 0d0
               ddwwb(i,j) = 0d0
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
         
c get the real emission amplitude squared, store it in uuww, ddww
      do mu = 0,3
         do i = 1,3
            pbar(mu,physToDiag(i))=p(mu,i,1)
         enddo
         qbar(mu) = pbar(mu,3)
      enddo
      qbar(4)=0d0

      do mu = 0,3            ! kinematics for H-->WW -->4 lepton decay
         pbar(mu,3) = v(mu,1,1) ! ve
         pbar(mu,4) = v(mu,2,1) ! e+
         pbar(mu,5) = v(mu,3,1) ! mu-
         pbar(mu,6) = v(mu,4,1) ! nu_mu~
         pbar(mu,7) = qbar(mu)
      enddo
      fsign(7) = gsign

c calculate the matrix elements pa + pb -> ve e+ mu- vm~ p1
      if (lok(1) .or. nlo.eq.0) then

#ifdef WITH_MADGRAPH
         if (ldebug) then
            call suu_wwg(pbar,fsign,uuwwg)          !madgraph
            call sdd_wwg(pbar,fsign,ddwwg)
         endif
#endif

            call qqVVj(pbar,fsign,qbar,gsign,0,1,uuww,ddww,ps_number)

#ifdef WITH_MADGRAPH
         if (ldebug) then
            print*,"uuwwg =   ",uuwwg
            print*,"uuww = ",uuww
            print*,"uuwwg/uuww = ",uuwwg/uuww
            print*,"ddwwg =   ",ddwwg
            print*,"ddww = ",ddww
            print*,"ddwwg/ddww = ",ddwwg/ddww
         endif
#endif

      endif

      NCmatrixelt(0,2,1)=uuww
      NCmatrixelt(1,2,1)=ddww

c for the NLO case get the subtraction terms; 
c first the case with a final state gluon (gsign = +1)

      if (gsign.eq.1) then
         if (ldebug) then
            print*," final state gluon section in qq_VVj "
            print*," jsig = ",jsig," jsig1,3 = ",jsig1,jsig3
            print 101," fsign = ",fsign
 101        format(a,6i5,a,i5)
         endif

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

c calculation of matrix elemnets for q qb -> W W g:

               if(lok(j)) then
#ifdef WITH_MADGRAPH
                  if(ldebug) then
                     call suu_ww(pbar,fsign,uuww2)  !madgraph
                     call sdd_ww(pbar,fsign,ddww2)
                  endif
#endif

                  call QQVV(pbar,fsign,0,J,uuww,ddww)

#ifdef WITH_MADGRAPH
                  if (ldebug) then
                        print*,"uuww/uuww2 =   ",uuww/uuww2
                        print*,"ddww/ddww2 =   ",ddww/ddww2
                  endif
#endif

               else
                  uuww = 0d0
                  ddww = 0d0
               endif

               NCmatrixelt(0,j,2)=uuww   !Born for momentum
               NCmatrixelt(1,j,2)=ddww   !configuration 2 and 3

c save matrix elements for later use with initial gluons
               uuwwb(physToDiag(1),j) = uuww
               ddwwb(physToDiag(1),j) = ddww

               q2 = 2d0*xuz*dotrr(qbar,p(0,j-1,1)) !p(mu,j-1,1) is inc.quark
               v2 = 2d0*gs2(j)*                   ! 8 pi alfas
     &             ( 2d0/(1-xuz)-(1+xuz) )
               sub(j) = v2/q2

            enddo
         endif !nlo.eq.1
#endif

         iflav(3) = 0           ! final state gluon id
         do if1=1,(nfl/2)*2
            iflav(1)=if1*fsign(physToDiag(1))
            iflav(2)=if1*fsign(physToDiag(2))
               k=FL_VVg(iflav,2)
               if(lok(1)) then
                  res(k,1)=pdf(sign1*iflav(1),1,1) 
     &                 *pdf(sign2*iflav(2),2,1) 
     &                 *NCmatrixelt(mod(if1,2),2,1)*polcolq
               else
                     res(k,1) = 0
               endif

               if ( (nlo.eq.1) .and. (lok(2).or.lok(3)) ) then
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
                  if(ldebug) then                              !begin debug
                     temp = qbar(0)              !soft
                     if(temp.le.1.0d0) then
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
                           print*,"res(k,2)+res(k,3)=",res(k,2)+res(k,3)
                           print*,"res(k,0) = ",res(k,0)
                           print*,"|res(k,1)/res(k,2)| = ",temp
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
               NCmatrixelt(0,2,2)=uuwwb(1,2)
               NCmatrixelt(1,2,2)=ddwwb(1,2)
               NCmatrixelt(0,3,2)=uuwwb(1,3)
               NCmatrixelt(1,3,2)=ddwwb(1,3)
            else
               NCmatrixelt(0,2,2)=uuwwb(2,2)
               NCmatrixelt(1,2,2)=ddwwb(2,2)
               NCmatrixelt(0,3,2)=uuwwb(2,3)
               NCmatrixelt(1,3,2)=ddwwb(2,3)
            endif
#endif
         endif !nlo

         do if1=1,(nfl/2)*2
            iflav(j-1) = 0
            iflav(4-j)=if1*fsign(physToDiag(4-j))
            iflav(3)  =if1*fsign(physToDiag(3))
            k=FL_VVg(iflav,j)
            if (lok(1)) then
               res(k,1)=pdf(sign1*iflav(1),1,1)  
     &              *pdf(sign2*iflav(2),2,1)    
     &              *NCmatrixelt(mod(if1,2),2,1)*polcolg
            else
               res(k,1) = 0
            endif

c--------------------------------------------------------------

            if ( (nlo.eq.1) .and. (lok(2).or.lok(3)) ) then
#ifdef WITH_NLO
               res(k,j) = ( pdf(sign1*iflav(1),1,j)     ! collinear
     &                       *pdf(sign2*iflav(2),2,j))
     &                       *NCmatrixelt(mod(if1,2),j,2)
               res(k,j) = ( pdf(sign1*iflav(1),1,j)     ! collinear
     &                       *pdf(sign2*iflav(2),2,j)
     &                       -pdf(sign1*iflav(1),1,2*j-3)  ! real emission
     &                       *pdf(sign2*iflav(2),2,4-j )*sub(j))
     &                       *NCmatrixelt(mod(if1,2),j,2)

              res(k,j) = res(k,j)*C2*polcolg
              res(k,-j+5) = 0d0

              res(k,0) = res(k,1) + res(k,j)

              if (ldebug) then                         !begin debug
                 temp = p(0,3,1)                 !soft
                 if(temp.le.0.001d0) then
c                     temp = abs(dotrr(p(0,3,1),p(0,1,1)))  !collinear
c                     temp2 = abs(dotrr(p(0,3,1),p(0,2,1)))
c                     if(temp.le.0.0001.or.temp2.le.0.0001) then
c                     if(j.eq.3) then
c                        temp = abs((res(k,1))+res(k-4,1))/
c     &                             (res(k,2)+res(k-4,2))
c                        if(temp.le.0.9.or.temp.ge.1.1) then
c                           print*," "
c                           print*," gluon in the initial state"
c                           print*," k = ",k," if1 = ",if1      
                        print*,"p(0,3,1) = ",p(0,3,1)
c                           print*,"pa*p1= ",dotrr(p(0,3,1),p(0,1,1))
c                           print*,"pb*p1= ",dotrr(p(0,3,1),p(0,2,1))
                        print*,"res(k,1) = ",res(k,1)
                        print*,"res(k,",j,") = ",res(k,j)
c                     print*,"res(k,1)+res(k-4,1) = ",res(k,1)+res(k-4,1)
c                     print*,"res(k,2)+res(k-4,2) = ",res(k,2)+res(k-4,2)
c                           endif
c                        endif
                 endif
              endif                                   !end debug
#endif

            else
               res(k,0) = res(k,1)
               res(k,2) = 0
               res(k,3) = 0
            endif

         enddo
      endif

      nmax = FL_VVg(iflav,-2)

      end

c******************************************************************
c
c   end subroutine qq_VVj
c
c*****************************************************************

c******************************************************************
c
c   begin function FL_VVg
c
c*****************************************************************
      integer function FL_VVg(iflav,colstruc)
      implicit none
      integer iflav(3),colstruc    ! input for color and flavor assignment
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
c  if colstruc=-1, reset listposition=0
c  if colstruc=-2, return listposition without incrementing

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
          call fillColoredPartons_VVg(iflav,
     1                                listposition+numdecay,colstruc)
        endif       ! numresets.eq.0
      endif         
      FL_VVg=listposition
      end
c******************************************************************************
c
c   end function FL_VVg
c
c******************************************************************************

c*****************************************************************************
c
c    begin  subroutine fillColoredPartons_VVg
c
c*****************************************************************************
      subroutine fillColoredPartons_VVg(iflav,listposition,colstruc)
c  assigns values to the variables in the common block localHEPUP
c  in particular, this subroutine assigns values to those variables that 
c  will be stored in the lookup tables generated by writeHEPUPtable.  
c  As the name suggests, this routine only stores the information for the 
c  colored partons.  Particles without color will be dealt with in the 
c  subroutine fillColorless.

      implicit none

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
      case(WWjLO,WPZjLOb,WMZjLOb,ZZjLO,WPAjLOb,WMAjLOb,ZAjLO,AAjLO,
     &     WPhadWMlepjLO,WPlepWMhadjLO,ZZhadjLO,
     &     WPhadZjLO,WPZhadjLO,WMhadZjLO,WMZhadjLO)
         lnup(listposition)=numParticles+1   
      case(WPj,WMj,WPjLOb,WMjLOb)
         lnup(listposition)=6
      case(WPHj,WMHj,WPHjLOb,WMHjLOb)
         lnup(listposition)=7
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
         if(id1.eq.21) then
               licolup(1,1,listposition)=501
               licolup(2,1,listposition)=502
               licolup(1,2,listposition)=502
               licolup(2,2,listposition)=503
               licolup(1,3,listposition)=501
               licolup(2,3,listposition)=503
         else if(id1.lt.0) then
               licolup(1,1,listposition)=0
               licolup(2,1,listposition)=501
               licolup(1,2,listposition)=502
               licolup(2,3,listposition)=0
               licolup(1,3,listposition)=502
               licolup(2,3,listposition)=501               
         else
               licolup(1,1,listposition)=501
               licolup(2,1,listposition)=0
               licolup(1,2,listposition)=0
               licolup(2,2,listposition)=502
               licolup(1,3,listposition)=501
               licolup(2,3,listposition)=502
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

      end  ! fillColoredPartons_VVg
c*****************************************************************************
c
c    end subroutine fillColoredPartons_VVg
c
c*****************************************************************************
