c*********************  subroutine m2s_qqWHj*********************************
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
c  Adapted from m2s_WAj (Johannes Bellm)
c
c  Robin Roth <robin@particle.physik.uni-karlsruhe.de>
c	Initial version: November 2012
c	Last modified: November 2012
c
c***************************************************************************

      subroutine m2s_qqWHj(bos,nlo,lok,xi,p,v,rn,xuz,m2s)

         use globalvars, only: lglobalprint
      implicit none

c declare input/output variables

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/lha.inc"

      integer bos
      real*8 xi(nx), p(0:3,max_p,max_kin), v(0:3,max_v,max_kin), 
     1       rn(1), xuz, m2s(0:3)
      logical lok(3)

c  helicity selection

      integer jsig, jsig1, jsig3, jsig5
      common /chelsum/ jsig,jsig1,jsig3,jsig5

c  helicity selection

      DOUBLE PRECISION zero
      PARAMETER (zero=0d0)

c declare external functions

      integer FL_VVg
      external FL_VVg

c alfas, scales etc

      real*8 x1,x2
      real*8 lnQomu(2:3), omxi(2:3), 
     &       Ax(2:3), Bx(2:3), Cx(2:3), Dxi(2:3), tgs2oqsq(2:3), 
     &       ln1mxi, lnrat, z, lnz, qwph(0:4,1:3)

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
      real*8 pdfckm(-6:6,2,3)

c define program switches

      logical ldebug
      parameter (ldebug=.false.)

      real*8 weight,rnumb,RandomNumber

      ! switch between whj and wj
      logical hashiggs

      ! Higgs propagator factor
      double precision fac
      ! constants needed for fac from common block
      DOUBLE PRECISION CLR,XM2,XMG,B
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6)

c --------------------------------------------------------------

c if first iteration, output acceptance cuts and fix constant input 
c parameters for call of qqWAj.f

      if ( init .ne. 0) then
         do i = 1,nmax
            do j=0,3
               res(i,j) = 0
            enddo
         enddo
      else if (init.eq.0) then
         if (lglobalprint) write(6,*) " "
         if(bos.eq.36) then
            if (lglobalprint) then
               write(6,*) "W+H j amplitude square information:"
            endif
            hashiggs = .true.
         elseif(bos.eq.46) then
            if (lglobalprint) then
               write(6,*) "W-H j amplitude square information:"
            endif
            hashiggs = .true.
         elseif(bos.eq.3) then
            if (lglobalprint) then
               write(6,*) "W+ j amplitude square information:"
            endif
            hashiggs = .false.
         elseif(bos.eq.4) then
            if (lglobalprint) then
               write(6,*) "W- j amplitude square information:"
            endif
            hashiggs = .false.
         else
            write(6,*) " wrong call of m2s_qqWHj, bos = ",bos
            stop
         endif

         if (lglobalprint) then
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
         endif
         
         call wtowh_reset

         call printnfl(.true.)

         if (lglobalprint) then
            print*," "
            print*," creal = ",crealVV," cvirtual = ",cvirtVV
         endif
         init = 1
         do i = 1,maxnumsubproc
            do j = 0,3
               res(i,j) = 0
            enddo
         enddo
      endif

c ------------------------------------------------------------------- 

c for WA->2l+A precalculate W->WA leptonic tensors
c lepton spinors and W+- polarization vectors

c...Les Houches interface
        if ((lha.or.hepmc).and..not.doNLO) then
           helicity(1)=-1
           helicity(2)= 1
        endif

      do J = 1,3
         do mu = 0,3
            if(hashiggs) then
               qp(mu,J) = v(mu,1,J)+v(mu,2,J)   !W+
!               qh(mu,J) = v(mu,3,J) 
               qh(mu,J) = 0d0
               do i=3,n_v ! sum all higgs decay products
                  qh(mu,J) = qh(mu,J) + v(mu,i,J)
               enddo
               qwph(mu,J) = qp(mu,J) + qh(mu,J) + p(mu,3,J)
            else
               qwph(mu,J) = v(mu,1,J)+v(mu,2,J)   !W+
            endif
         enddo
         if (hashiggs) then
            qp(4,J) = dotrr(qp(0,J),qp(0,J))
            qh(4,J) = dotrr(qh(0,J),qh(0,J))
         endif
         qwph(4,J) = dotrr(qwph(0,J),qwph(0,J))


         if (bos.eq.3 .or. bos.eq.36) then
            CALL OXXXXX(v(0,1,J),ZERO ,-1,1,wve(1,J))           !ve 
            CALL IXXXXX(v(0,2,J),ZERO ,1,-1,wep(1,J))           !e+     
         elseif (bos.eq.4 .or. bos.eq.46) then
            CALL IXXXXX(v(0,2,J),ZERO ,1,-1,wve(1,J))           !ve~  
            CALL OXXXXX(v(0,1,J),ZERO ,-1,1,wep(1,J))           !e-    
         endif

         if (hashiggs) then
            CALL SXXXXX(qh(0,J),1,hig(1,J))      
         endif

         !can be moved out of J loop
      if(hashiggs) then
         fac = 1.d0 / ( (qh(4,1)-xm2(6))**2 + xmg(6)**2 ) !propagator^2 only
         fac = fac *16*pi*xmg(6)
      else
         fac = 1d0
      endif

      if(hashiggs) then
         if (with_anom) then ! anomalous gauge boson couplings 
            if (bos.eq.36) then
               call wpmtowh_an(v(0,1,J),J,wtowh(0,J),3)
            else if (bos.eq.46) then
               call wpmtowh_an(v(0,1,J),J,wtowh(0,J),4)
            endif
         else  ! SM
            if (ldebug) print*, "wtowh for bos ", bos
            if (bos.eq.36) then
               call wpmtowh(v(0,1,J),J,wtowh(0,J),3)
            else if (bos.eq.46) then
               call wpmtowh(v(0,1,J),J,wtowh(0,J),4)
            endif
         endif 
      else
         if (bos.eq.3) then
            CALL JIOXXX(wep(1,J),wve(1,J),GWF ,WMASS,WWIDTH,wtowh(0,J)) !!W+
         else if (bos.eq.4) then
            CALL JIOXXX(wve(1,J),wep(1,J),GWF ,WMASS,WWIDTH,wtowh(0,J)) !!W-
         endif
      endif


      enddo

c call PDF subroutine in order to determine parton
c distributions in the incoming (anti)protons.
c for the NLO contributions x1 = x*y in my notes with x=xuz, y=xi(i)

      x1 = xi(1)*xuz
      q_sf = sqrt(mufsq(1,1))
      call pdfproton( xi(1), q_sf, pdf(-6,1,1) )      ! f_a(y)=f_a(x1/z)
      if(nlo.gt.0) then
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


c TODO: combine with above for beam1 using j=1,2 and 4-j
      x2 = xi(2)*xuz
      q_sf = sqrt(mufsq(2,1))
      call pdfproton( xi(2), q_sf, pdf(-6,2,1) )  ! f_b(y)=f_a(x2/z)
      if(nlo.gt.0) then
         q_sf = sqrt(mufsq(2,3))                 ! f_b(x2) for lower line 
         call pdfproton( x2, q_sf, pdf(-6,2,3) ) !   NLO correction

         if (mufsq(2,2).ne.mufsq(2,1)) then      ! f_b(x2) for upper line
            q_sf = sqrt(mufsq(2,2))              !   NLO correction
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

      lnQomu(2) = log(qwph(4,2)/mufsq(1,2))
      lnQomu(3) = log(qwph(4,3)/mufsq(2,3))

      omxi(2) = 1d0-x1
      omxi(3) = 1d0-x2

c      factor 4*pi*alfas/(2x pa pb ) = 4*pi*alfas/Q^2

      tgs2oqsq(2) = 4d0*pi*als(1,2)/qwph(4,2)
      tgs2oqsq(3) = 4d0*pi*als(2,3)/qwph(4,3)

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


      do i=1,3
         call pdfckmtransform(pdf(-6,1,i),pdfckm(-6,1,i))
      enddo


C********************************************************************

C     Now sum over the subprocesses contributing to W+ W- production

      nmax = FL_VVg(fsign,-1)          !reset counter for subprocesses to 0

C*******************  q1 qbar2 ---> W H g   **********************
! here everything done with W+ in mind
! signchange for W- will come in qqWHj      

c   physToDiag(ext.momentum label) = Feynman diagram label
      

      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=2
      physToDiag(3)=3    ! gluon

      fsign(1) = 1          !incoming quark
      fsign(2) = -1         !incoming antiquark
      fsign(3) = 1          !outgoing neutrino
      fsign(4) = -1         !outgoing antilepton

      gsign    = 1          !outgoing gluon

      call qq_WHj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     &              pdf,pdfckm,res,nmin,nmax)

C*******************  qbar1 q2 ---> V V g   **********************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=1
      physToDiag(3)=3    ! gluon

      fsign(1) = 1
      fsign(2) = -1

      call qq_WHj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     &              pdf,pdfckm,res,nmin,nmax)

C*******************  g q1 ---> V V q2  ***********************************
      
      physToDiag(1)=3
      physToDiag(2)=1
      physToDiag(3)=2

      fsign(1) = 1
      fsign(2) = 1
      gsign    = -1

      call qq_WHj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1              pdf,pdfckm,res,nmin,nmax)

C*******************  qbar2 g ---> V V qbar1   **********************

      physToDiag(1)=2
      physToDiag(2)=3
      physToDiag(3)=1

      fsign(1) = -1
      fsign(2) = -1
      gsign = -1

      call qq_WHj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1              pdf,pdfckm,res,nmin,nmax)

C*******************  g q2bar ---> V V q1bar  **********************

      physToDiag(1)=3
      physToDiag(2)=2
      physToDiag(3)=1

      fsign(1) = -1
      fsign(2) = -1
      gsign = -1

      call qq_WHj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1              pdf,pdfckm,res,nmin,nmax)

C*******************  q1 g ---> V V q2   **********************

      physToDiag(1)=1
      physToDiag(2)=3
      physToDiag(3)=2

      fsign(1) = 1
      fsign(2) = 1
      gsign    = -1 

      call qq_WHj(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1              pdf,pdfckm,res,nmin,nmax)

c*****************  end of process evaluation  *********************

      if (init.eq.1) then
         init = init+1
      endif
      do j=0,3
         m2s(j) = 0
         DO IRES = 1,NMAX

            m2s(j) = m2s(j) + RES(IRES,j) 
           
         ENDDO
  
         m2s(j)= m2s(j)*fac  !multipy higgs propagator
c         m2s(j)= m2s(j)*2 ! factor 2 for random helicity summation

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
            endif
         endif

      enddo

      RETURN
      END

c*****************************************************************************
c
c    end subroutine m2s_qqWHj
c
c*****************************************************************************
c******************************************************************
c
c   begin subroutine qq_WHj
c
c*****************************************************************

      subroutine qq_WHj(xi,p,v,physToDiag,fsign,gsign,bos, 
     1                    nlo,lok,xuz,
     1                    pdf,pdfckm,res,nmin,nmax)
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"

      real*8 p(0:3,max_p,max_kin), v(0:3,max_v,max_kin)
      real*8 xi(nx), xuz
      real*8 pdf(-6:6,2,3), pdfckm(-6:6,2,3), res(maxnumsubproc,0:3)
      integer physToDiag(3), fsign(3+n_v), 
     1        gsign, bos, nlo, nmin, nmax
      logical lok(3)
      double precision getpdfckm
c
c qq_WPZj calls the amplitude square routines 
c             qqWAj    for qq --> WA g --> g 2 leptons +A
c                          gq --> WA q --> q 2 leptons +A
c
c for the subtraction terms for the NLO cross section calculation
c
c  INPUT:  p(0:3,5,3)      external physical parton momenta
c          v(0:3,nv,3)     Lepton and Higgs momenta
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

c Note: Factor 3 from color sum included in qqWA. I am using color summed
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
      real*8 udwhb(2,2:3),
     &       udwh, udwhg

      ! 6 = 2 (in) + 1 jet + 3 final state (l, v, h)
      real*8 pbar(0:3,3+n_v), qbar(0:4),q2,v2,sub(2:3), 
     1       dotrr, gs2(2:3), polcolq, polcolg

      double precision NCmatrixelt(3,2), temp
      integer iflav(3), diagToPhys(3), FL_VVg
      ! variables used for including ckm matrix in pdfs
      integer f1,f2

      external dotrr, FL_VVg

      complex*16 zero
      parameter (zero=(0d0,0d0))

      logical ldebug
      save gs2,polcolq,polcolg, udwhb
      integer i,if1,j,k,mu
      
      parameter (ldebug=.false.)

c--------------------------------------------------------------------

      if(ldebug) then 
         print*, "--------- enter m2s_qqWHj  qq_whj-------------------"
         print*, "lok ", lok
         print*, "nlo ", nlo
         print*, "xuz ", xuz
      endif

      udwh  = 0d0
      udwhg = 0d0

c  if this is a new phase space point -->
c  reset the LO amplitude to 0 to avoid wrong subtraction
      if (physToDiag(1).eq.1.and.physToDiag(2).eq.2) then
         do i = 1,2
            do j = 2,3
               udwhb(i,j) = 0d0
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
         
c get the real emission amplitude squared, store it in udwh, ddww
      do mu = 0,3
         do i = 1,3
            pbar(mu,physToDiag(i))=p(mu,i,1)
         enddo
         qbar(mu) = pbar(mu,3)
      enddo
      qbar(4)=0d0 ! massless jet

      do mu = 0,3            ! 
         pbar(mu,3) = v(mu,1,1) ! ve or e-
         pbar(mu,4) = v(mu,2,1) ! e+ or ve~
!         pbar(mu,5) = v(mu,3,1) ! H !TODO: fix for madgraph comparison
!         pbar(mu,6) = qbar(mu)
      enddo
!      fsign(6) = gsign

c calculate the matrix elements pa + pb -> ve e+ H p1
      if (lok(1) .or. nlo.eq.0) then


! fsign should not be flipped for W- case         
! Instead labels on fermionlines move around 
! u->d, dbar -> ubar, ve -> e-, e+ -> ve~

         call qqWHj(pbar,fsign,qbar,gsign,udwh,bos)

c debugging in subprocesses

      endif
      
      NCmatrixelt(2,1)=udwh


c for the NLO case get the subtraction terms; 
c first the case with a final state gluon (gsign = +1)

      if (gsign.eq.1) then
         if (ldebug) then      
c             print*," final state gluon section in qq_WPZj "
c             print*," jsig = ",jsig," jsig1,3 = ",jsig1,jsig3
c             print 101," fsign = ",fsign
c  101        format(a,6i5,a,i5)
         endif
#ifdef WITH_NLO

         if (nlo.eq.1) then
            do j = 2,3               ! j=2: emission off 1,2 line
               do mu = 0,3           ! j=3: emission off 3,4 line
                  do i = 1,2
                     pbar(mu,physToDiag(i))=p(mu,i,j)
                  enddo
                  do i = 1,3
                     pbar(mu,i+2)=v(mu,i,j)
                  enddo
               enddo

c calculation of matrix elemnets for q qb -> W W g:

               if(lok(j)) then

                  call qqWH(pbar,fsign,0,J,udwh,bos)
                  if (ldebug) print *, "j, udwh", j, " ", udwh

c qqWH is debugged in m2s_qqWH.F for J=1 

               else
                  udwh = 0d0
               endif

               NCmatrixelt(j,2)=udwh   !Born for momentum
               !configuration 2 and 3


c save matrix elements for later use with initial gluons
               udwhb(physToDiag(1),j) = udwh

               q2 = 2d0*xuz*dotrr(qbar,p(0,j-1,1)) !p(mu,j-1,1) is inc.quark
               v2 = 2d0*gs2(j)*                   ! 8 pi alfas
     &             ( 2d0/(1-xuz)-(1+xuz) )
               sub(j) = v2/q2
               if (ldebug) print*,"subj ",sub(j)

            enddo
         endif !nlo.eq.1
#endif


         iflav(3) = 0           ! final state gluon id
         do if1=2,(nfl/2)*2,2

            if (bos.eq.3 .or. bos.eq.36) then
               iflav(1)=( if1 - physToDiag(1) + 1 )*fsign(physToDiag(1))
               iflav(2)=( if1 - physToDiag(2) + 1 )*fsign(physToDiag(2))
            else if (bos.eq.4 .or. bos.eq.46) then
               iflav(1)=( if1 + physToDiag(1) - 2 )*fsign(physToDiag(1))
               iflav(2)=( if1 + physToDiag(2) - 2 )*fsign(physToDiag(2))
            endif

            f1=sign1*iflav(1)
            f2=sign2*iflav(2)

            if (ldebug) print*, "iflav ", iflav

            k=FL_VVg(iflav,2)
            if(lok(1)) then 
               res(k,1)= getpdfckm(pdf(-6,1,1),pdfckm(-6,1,1),f1,f2)
     &                 *NCmatrixelt(2,1)*polcolq
            else
               res(k,1) = 0
            endif

            if ( (nlo.eq.1) .and. (lok(2).or.lok(3)) ) then
#ifdef WITH_NLO
               f1=sign1*iflav(1)
               f2=sign2*iflav(2)
               ! collinear - real emission

               res(k,2) = 
     &         (getpdfckm(pdf(-6,1,2),pdfckm(-6,1,2),f1,f2) 
     &   -0.5d0*(pdf(f1,1,1)*pdfckm(f2,2,2)+pdfckm(f1,1,1)*pdf(f2,2,2))*sub(2))
     &                     *C2*NCmatrixelt(2,2)

               res(k,3) = 
     &         (getpdfckm(pdf(-6,1,3),pdfckm(-6,1,3),f1,f2) 
     &   -0.5d0*(pdf(f1,1,3)*pdfckm(f2,2,1)+pdfckm(f1,1,3)*pdf(f2,2,1))*sub(3))
     &                     *C2*NCmatrixelt(3,2)

!              res(k,2) = ( pdf(sign1*iflav(1),1,2)     ! collinear
!    &                     *pdf(sign2*iflav(2),2,2)
!    &                     -pdf(sign1*iflav(1),1,1)  ! real emission
!    &                     *pdf(sign2*iflav(2),2,2 )*sub(2))*C2
!    &                     *NCmatrixelt(2,2)

!              res(k,3) = ( pdf(sign1*iflav(1),1,3)     ! collinear
!    &                     *pdf(sign2*iflav(2),2,3)
!    &                     -pdf(sign1*iflav(1),1,3)  ! real emission
!    &                     *pdf(sign2*iflav(2),2,1)*sub(3))*C2
!    &                     *NCmatrixelt(3,2)
               res(k,2) = res(k,2) * polcolq
               res(k,3) = res(k,3) * polcolq
               res(k,0) = res(k,1) + res(k,2) + res(k,3)


c debugging for collinear and soft divergences -----------------------
               if(ldebug) then                              !begin debug
                  temp = qbar(0)              !soft

                  if(temp.le.1.0d0) then

!                      print*,"jiha",lok(1)
!                      print*,lok(1),lok(2),lok(3)
c                     temp = abs(dotrr(qbar,pbar(0,1)))  !collinear
c                     temp2 = abs(dotrr(qbar,pbar(0,2)))
c                     if(temp.le.0.0001.or.temp2.le.0.0001) then
                     temp = abs(res(k,1)/(res(k,2)+res(k,3)))
!                                              print*,temp
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

            else ! nlo.eq.1 .and. ( lok(2) .or. lok(3))
               res(k,0) = res(k,1)
               res(k,2) = 0
               res(k,3) = 0
            endif
         enddo  

      elseif (gsign.eq.-1) then                !initial gluon section
         j = 0
         if (nlo.eq.0) then
            if (physToDiag(1).eq.3) then
               j=2                             ! j=2: emission off 1,2 line
            elseif (physToDiag(2).eq.3) then
               j=3                             ! j=3: emission off 3,4 line
            endif
         elseif (nlo.eq.1) then
#ifdef WITH_NLO
            if (physToDiag(1).eq.3) then
               j=2                             ! j=2: emission off 1,2 line
            elseif (physToDiag(2).eq.3) then
               j=3                             ! j=3: emission off 3,4 line
            endif

            q2 = 2d0 * dotrr(p(0,3,1),p(0,j-1,1))*xuz
            v2 = 2d0*gs2(j) * ( (1-xuz)**2 + xuz**2 )
            sub(j) = v2/q2

            if (physToDiag(1).eq.1.or.PhysToDiag(2).eq.2) then
               NCmatrixelt(2,2)=udwhb(1,2)
               NCmatrixelt(3,2)=udwhb(1,3)
            else
               NCmatrixelt(2,2)=udwhb(2,2)
               NCmatrixelt(3,2)=udwhb(2,3)
            endif
#endif
         endif !nlo


         iflav(j-1) = 0

         do if1=2,(nfl/2)*2,2


            if (bos.eq.3 .or. bos.eq.36) then
               iflav(4-j)=( if1 - physToDiag(4-j) + 1 )*fsign(physToDiag(4-j))
               iflav(3)  =( if1 - physToDiag(3) + 1 )*fsign(physToDiag(3))
            else if (bos.eq.4 .or. bos.eq.46) then
               iflav(4-j)=( if1 + physToDiag(4-j) - 2 )*fsign(physToDiag(4-j))
               iflav(3)  =( if1 + physToDiag(3) - 2 )*fsign(physToDiag(3))
            endif


            f1=sign1*iflav(1)
            f2=sign2*iflav(2)

            k=FL_VVg(iflav,j)


            if (ldebug) print*, "iflav ", iflav
            if (lok(1)) then
               ! factor 2 removes factor 0.5 from getpdfckm
               ! this 0.5 is not needed as one parton is a gluon and
               !    therefore pdfckm=0 
               ! compiler will optimize this away
               res(k,1)= 2d0*getpdfckm(pdf,pdfckm,f1,f2)
     &                 *NCmatrixelt(2,1)*polcolg
            else
               res(k,1) = 0
            endif

c--------------------------------------------------------------

            if ( (nlo.eq.1) .and. (lok(2).or.lok(3)) ) then
#ifdef WITH_NLO
            f1=sign1*iflav(1)
            f2=sign2*iflav(2)

               ! collinear - real emission
!              print*, 'fristline'
!              print*, pdf(f1,1,j), pdfckm(f2,2,j)
!              print*, pdfckm(f1,1,j), pdf(f2,2,j)
!              print*, pdf(f1,1,2*j-3),pdfckm(f2,2,4-j)
!              print*, pdfckm(f1,1,2*j-3),pdf(f2,2,4-j)

               ! factor 2d0 see above
               res(k,j) = 
     &         ( 2d0*getpdfckm(pdf(-6,1,j),pdfckm(-6,1,j),f1,f2) 
     &         - ( pdfckm(f1,1,2*j-3)*pdf(f2,2,4-j) + pdf(f1,1,2*j-3)*pdfckm(f2,2,4-j) )*sub(j) )
     &                       *NCmatrixelt(j,2)

!              print*, 'new'
!              print*, res(k,j)

!              res(k,j) = ( pdf(sign1*iflav(1),1,j)        ! collinear
!    &                       *pdf(sign2*iflav(2),2,j)
!    &                       -pdf(sign1*iflav(1),1,2*j-3)  ! real emission
!    &                       *pdf(sign2*iflav(2),2,4-j )*sub(j))
!    &                       *NCmatrixelt(j,2)
!              print*, res(k,j)
               if (ldebug) then
                  print*, "pdfs ",
     &                     pdf(sign1*iflav(1),1,j), 
     &                     pdf(sign2*iflav(2),2,j),
     &                     pdf(sign1*iflav(1),1,2*j-3),
     &                     pdf(sign2*iflav(2),2,4-j )
                  print*, "reskj ", res(k,j)
                  print*, "subj ", sub(j)
                  print*, "polcolg ", polcolg
                  print*, "j, Mat ", j, NCmatrixelt(j,2)
               endif

               res(k,j) = res(k,j)*C2*polcolg
               res(k,-j+5) = 0d0

               res(k,0) = res(k,1) + res(k,j)

               if (ldebug) then                            !begin debug
                  temp = p(0,3,1)                          !soft
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
c            if (ldebug) print*, "resk ",res(k,0:3)

         enddo
      endif

      nmax = FL_VVg(iflav,-2)

      end

c******************************************************************
c
c   end subroutine qq_WPZj
c
c*****************************************************************

