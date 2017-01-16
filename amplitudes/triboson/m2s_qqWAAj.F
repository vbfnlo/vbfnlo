c*********************  subroutine m2s_qqWAAj*********************************
c
c INPUT     bos       Boson identifier, 1 = A, 3,4=W
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
c       giuseppe bozzi, <giuseppe@particle.uni-karlsruhe.de>
c       Last modified: December 2009
c
c***************************************************************************

      SUBROUTINE m2s_qqWAAj(bos,nlo,lok,xi,p,v,rn,xuz,m2s,ps_number,NW)

      IMPLICIT NONE

c declare input/output variables

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/lha.inc"

      INTEGER bos, nw,ps_number
      DOUBLE PRECISION xi(nx),p(0:3,max_p,max_kin),v(0:3,max_v,max_kin), 
     1       rn(1), xuz, m2s(0:3)
      logical lok(3)

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
         if (bos.eq.311) then
            write(6,*) " "
            write(6,*) "AAW+j amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
         elseif (bos.eq.411) then
            write(6,*) " "
            write(6,*) "AAW-j amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
         endif
         call vtowaa_reset
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
            qw0(mu,J) = v(mu,1,J)+v(mu,2,J)   !W+-
            qa1(mu,J) = v(mu,3,J)             !A1
            qa2(mu,J) = v(mu,4,J)             !A2
            qaa(mu,J) = qa1(mu,J) + qa2(mu,J)
            qaaw(mu,J) = qa1(mu,J) + qa2(mu,J) + qw0(mu,J)
         enddo

         qw0(4,J) = dotrr(qw0(0,J),qw0(0,J))
         qa1(4,J) = dotrr(qa1(0,J),qa1(0,J))
         qa2(4,J) = dotrr(qa2(0,J),qa2(0,J))
         qaa(4,J) = dotrr(qaa(0,J),qaa(0,J))
         qaaw(4,J) = dotrr(qaaw(0,J),qaaw(0,J))
      enddo

c select helicity: h = random number for photon helicity

      h = mod(h,4) + 1
      ie = (-1)**(h+1)
      iu = (-1)**((h-1)/2)

c Les Houches interface
      if((lha.or.hepmc).and.(.not.doNLO)) then
         helicity(1) = -1
         helicity(2) =  1
         helicity(3) =  iu
         helicity(4) =  it
      endif

c lepton spinors and W+,A polarisation vectors
      do J = 1,3
         CALL OXXXXX(v(0,1,J),ZERO ,-1,1,wve(1,J))   !ve 
         CALL IXXXXX(v(0,2,J),ZERO ,1,-1,wep(1,J))   !e+        
         CALL VXXXXX(qa1(0,J),ZERO ,ie,1,atau1(1,J)) !A
         CALL VXXXXX(qa2(0,J),ZERO ,iu,1,atau2(1,J)) !A
         CALL JIOXXX(wep(1,J),wve(1,J),GWF,WMASS,WWIDTH, w0(1,J))  !W+


         if (with_anom) then ! anomalous gauge boson couplings
c           using global form factor for all tensors of one phase space point
c           this ensures proper cancellations for anomalous contributions
c           energy scale is invariant WAA mass
            call anomal_formfactor(qaaw(0,J),qw0(0,J),qa1(0,J),qa2(0,J))

c leptonic tensors
         if (bos.eq.311) then
            call wptowanew_anomal(v(0,1,J),J,h,wpwa1(0,J),1,NW)!W+ -> W+ A -> ve e+ a
            call wptowanew_anomal(v(0,1,J),J,h,wpwa2(0,J),2,NW)!W+ -> W+ A -> ve e+ a
            call wptowaa_anomal(v(0,1,J),J,h,wpwaa(0,J),0,NW)  !W+ -> ve e+ a a
         elseif (bos.eq.411) then
            call wmtowanew_anomal(v(0,1,J),J,h,wmwa1(0,J),1,NW)!W- -> W- A -> mu- vm~ a
            call wmtowanew_anomal(v(0,1,J),J,h,wmwa2(0,J),2,NW)!W- -> W- A -> mu- vm~ a
            call wmtowaa_anomal(v(0,1,J),J,h,wmwaa(0,J),0,NW)  !W- -> ve e+ a a
         endif


         else !SM

c leptonic tensors
         if (bos.eq.311) then
            call wptowanew(v(0,1,J),J,h,wpwa1(0,J),1,NW)!W+ -> W+ A -> ve e+ a
            call wptowanew(v(0,1,J),J,h,wpwa2(0,J),2,NW)!W+ -> W+ A -> ve e+ a
            call wptowaa(v(0,1,J),J,h,wpwaa(0,J),0,NW)  !W+ -> ve e+ a a
         elseif (bos.eq.411) then
            call wmtowanew(v(0,1,J),J,h,wmwa1(0,J),1,NW)!W- -> W- A -> mu- vm~ a
            call wmtowanew(v(0,1,J),J,h,wmwa2(0,J),2,NW)!W- -> W- A -> mu- vm~ a
            call wmtowaa(v(0,1,J),J,h,wmwaa(0,J),0,NW)  !W- -> ve e+ a a
         endif

         endif ! anom or SM

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

      lnQomu(2) = log(qaaw(4,2)/mufsq(1,2))
      lnQomu(3) = log(qaaw(4,3)/mufsq(2,3))
 
      omxi(2) = 1d0-x1
      omxi(3) = 1d0-x2

c factor 4*pi*alfas/(2x pa pb ) = 4*pi*alfas/Q^2

      tgs2oqsq(2) = 4d0*pi*als(1,2)/qaaw(4,2)
      tgs2oqsq(3) = 4d0*pi*als(2,3)/qaaw(4,3)

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

C     Now sum over the subprocesses contributing to W+ A A production

      nmax = FL_VVVg(fsign,-1)!reset counter for subprocesses to 0

C*******************  q1 qbar2 ---> W+ W- A g   **********************
         
c   physToDiag(ext.momentum label) = Feynman diagram label

      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=2
      physToDiag(3)=3    ! gluon

      fsign(1) = 1          !incoming quark
      fsign(2) = -1         !incoming antiquark
      fsign(3) = 1          !outgoing ve
      fsign(4) = -1         !outgoing ep
      fsign(5) = 1          !outgoing A1
      fsign(6) = 1          !outgoing A2
      gsign    = 1          !outgoing gluon

      call qq_WAAj(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     &              pdf,res,nmin,nmax,ps_number)

C*******************  qbar1 q2 ---> W+ W- A g   **********************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=1
      physToDiag(3)=3    ! gluon

      fsign(1) = 1
      fsign(2) = -1
      gsign    = 1

      call qq_WAAj(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     &              pdf,res,nmin,nmax,ps_number)

C*******************  g q1 ---> W+ W- A q2  ***********************************

      physToDiag(1)=3
      physToDiag(2)=1
      physToDiag(3)=2

      fsign(1) = 1
      fsign(2) = 1
      gsign    = -1

      call qq_WAAj(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     1              pdf,res,nmin,nmax,ps_number)

C*******************  qbar2 g ---> W+ W- A qbar1   **********************

      physToDiag(1)=2
      physToDiag(2)=3
      physToDiag(3)=1

      fsign(1) = -1
      fsign(2) = -1
      gsign    = -1

      call qq_WAAj(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     1              pdf,res,nmin,nmax,ps_number)

C*******************  g q2bar ---> W+ W- A q1bar  **********************

      physToDiag(1)=3
      physToDiag(2)=2
      physToDiag(3)=1

      fsign(1) = -1
      fsign(2) = -1
      gsign    = -1

      call qq_WAAj(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     1              pdf,res,nmin,nmax,ps_number)

C*******************  q1 g ---> W+ W- A q2   **********************

      physToDiag(1)=1
      physToDiag(2)=3
      physToDiag(3)=2

      fsign(1) = 1
      fsign(2) = 1
      gsign    = -1

      call qq_WAAj(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     1              pdf,res,nmin,nmax,ps_number)

c*****************  end of process evaluation  *********************

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
            endif
         endif

         m2s(j) = m2s(j) * 4d0    ! factor for random helicity summation
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
c    end subroutine m2s_qqWWAj
c
c*****************************************************************************
c******************************************************************
c
c   begin subroutine qq_WAAj
c
c*****************************************************************

      SUBROUTINE qq_WAAj(xi,p,v,physToDiag,fsign,gsign,bos, 
     1                    nlo,lok,xuz,
     1                    pdf,res,nmin,nmax,ps_number)
      IMPLICIT NONE

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"

      integer ps_number

      DOUBLE PRECISION p(0:3,max_p,max_kin), v(0:3,max_v,max_kin), xi(nx), xuz
      DOUBLE PRECISION pdf(-6:6,2,3), res(maxnumsubproc,0:3)
      INTEGER physToDiag(3), fsign(3+n_v), 
     1        gsign, bos, nlo, nmin, nmax
      LOGICAL lok(3)
c
c qq_VVV1j calls the amplitude square routines 
c             qqVVVj    for qq --> VVA g --> g 4 leptons
c                           gq --> VVA q --> q 4 leptons 
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

      DOUBLE PRECISION udwaab(2,2:3), duwaab(2,2:3), udwaa, duwaa
      DOUBLE PRECISION udwaa2, duwaa2
      DOUBLE PRECISION pbar(0:3,3+n_v), qbar(0:4),q2,v2,sub(2:3), 
     1       dotrr, gs2(2:3), polcolq, polcolg
      DOUBLE PRECISION NCmatrixelt(0:1,3,2), temp, temp1, temp2, temp3
      INTEGER iflav(3), diagToPhys(3), FL_VVVg
      EXTERNAL dotrr, FL_VVVg

      COMPLEX*16 zero
      PARAMETER (zero=(0d0,0d0))

      LOGICAL ldebug, lres, ldebugsoft
      SAVE gs2,polcolq,polcolg,udwaab,duwaab
      INTEGER i,if1,j,k,mu
      
      PARAMETER (ldebug=.false.)
      PARAMETER (ldebugsoft=.false.)
      PARAMETER (lres=.false.)

c--------------------------------------------------------------------

      udwaa  = 0d0
      duwaa  = 0d0
      udwaa2 = 0d0
      duwaa2 = 0d0

      if (lres) then
        open(unit=31,file="irfinite.chk",ACCESS='APPEND')
      endif

c  if this is a new phase space point -->
c  reset the LO amplitude to 0 to avoid wrong subtraction
      if (physToDiag(1).eq.1.and.physToDiag(2).eq.2) then
         do i = 1,2
            do j = 2,3
               udwaab(i,j) = 0d0
               duwaab(i,j) = 0d0
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
         
c get the real emission amplitude squared, store it in uuwwa, ddwwa
      do mu = 0,3
         do i = 1,3
            pbar(mu,physToDiag(i))=p(mu,i,1)
         enddo
         qbar(mu) = pbar(mu,3)
      enddo
      qbar(4)=0d0

      do mu = 0,3               ! kinematics for qq -> WWA -> 4 leptons + A
         pbar(mu,3) = v(mu,1,1) ! ve
         pbar(mu,4) = v(mu,2,1) ! e+
         pbar(mu,5) = v(mu,3,1) ! A1
         pbar(mu,6) = v(mu,4,1) ! A2
         pbar(mu,7) = qbar(mu)
      enddo
      fsign(7) = gsign

c calculate the matrix elements pa + pb -> ve e+ mu- vm~ a p1
      if (lok(1) .or. nlo.eq.0) then

#ifdef WITH_MADGRAPH
         if (ldebug) then
            if (bos.eq.311) then
               call sud_waag(pbar,fsign,udwaa2)
            else if (bos.eq.411) then
               call sdu_waag(pbar,fsign,duwaa2)
            endif
         endif
#endif
         
         if (bos.eq.311) then
            call qqAAWj(pbar,fsign,qbar,gsign,0,1,udwaa,bos,ps_number)
         elseif (bos.eq.411) then
            call qqAAWj(pbar,fsign,qbar,gsign,0,1,duwaa,bos,ps_number)
         endif
         
         if (ldebug) then
            if (bos.eq.311) then
               print*,"udwaa =   ",udwaa
#ifdef WITH_MADGRAPH
               print*,"udwaa2 =   ",udwaa2
#endif
            elseif (bos.eq.411) then
               print*,"duwaa =   ",duwaa
#ifdef WITH_MADGRAPH
               print*,"duwaa2 =   ",duwaa2
#endif
            endif
            read(*,*)
         endif
      endif

      NCmatrixelt(0,2,1)=udwaa
      NCmatrixelt(1,2,1)=duwaa

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

c calculation of matrix elements for q qb -> W W A g:

               if(lok(j)) then
                  if (bos.eq.311) then
                     call qqAAW(pbar,fsign,0,J,udwaa,bos)
                  elseif (bos.eq.411) then
                     call qqAAW(pbar,fsign,0,J,duwaa,bos)
                  endif
               else
                  udwaa = 0d0
                  duwaa = 0d0
               endif

               NCmatrixelt(0,j,2)=udwaa   !Born for momentum
               NCmatrixelt(1,j,2)=duwaa   !configuration 2 and 3

c save matrix elements for later use with initial gluons
               udwaab(physToDiag(1),j) = udwaa
               duwaab(physToDiag(1),j) = duwaa

               q2 = 2d0*xuz*dotrr(qbar,p(0,j-1,1)) !p(mu,j-1,1) is inc.quark
               v2 = 2d0*gs2(j)*                   ! 8 pi alfas
     &             ( 2d0/(1-xuz)-(1+xuz) )
               sub(j) = v2/q2

            enddo
         endif !nlo.eq.1
#endif /* WITH_NLO */

         iflav(3) = 0           ! final state gluon id
         do if1=1,(nfl/2)*2,2
            if (bos.eq.311) then
               iflav(1)=(if1-physToDiag(1)+2)*fsign(physToDiag(1))
               iflav(2)=(if1+physToDiag(1)-1)*fsign(physToDiag(2))
               k=FL_VVVg(iflav,2)
               if (lok(1)) then
                  res(k,1)=pdf(sign1*iflav(1),1,1)
     &                 *pdf(sign2*iflav(2),2,1)
     &                 *NCmatrixelt(0,2,1)*polcolq
               endif
            elseif (bos.eq.411) then
               iflav(1)=(if1+physToDiag(1)-1)*fsign(physToDiag(1))
               iflav(2)=(if1-physToDiag(1)+2)*fsign(physToDiag(2))
               k=FL_VVVg(iflav,2)
               if (lok(1)) then
                  res(k,1)=pdf(sign1*iflav(1),1,1)
     &              *pdf(sign2*iflav(2),2,1)
     &              *NCmatrixelt(1,2,1)*polcolq
               endif
            else
               res(k,1) = 0d0
            endif
            
            if ( nlo.eq.1 .and. (lok(2).or.lok(3)) ) then
#ifdef WITH_NLO

               if (bos.eq.311) then
                  res(k,2) = ( pdf(sign1*iflav(1),1,2) ! collinear
     &                 *pdf(sign2*iflav(2),2,2)
     &                 -pdf(sign1*iflav(1),1,1) ! real emission
     &                 *pdf(sign2*iflav(2),2,2 )*sub(2))*C2
     &                 *NCmatrixelt(0,2,2)
                  res(k,3) = ( pdf(sign1*iflav(1),1,3) ! collinear
     &                 *pdf(sign2*iflav(2),2,3)
     &                 -pdf(sign1*iflav(1),1,3) ! real emission
     &                 *pdf(sign2*iflav(2),2,1)*sub(3))*C2
     &                 *NCmatrixelt(0,3,2)
                  elseif (bos.eq.411) then
                  res(k,2) = ( pdf(sign1*iflav(1),1,2) ! collinear
     &                 *pdf(sign2*iflav(2),2,2)
     &                 -pdf(sign1*iflav(1),1,1) ! real emission
     &                 *pdf(sign2*iflav(2),2,2 )*sub(2))*C2
     &                 *NCmatrixelt(1,2,2)
                  res(k,3) = ( pdf(sign1*iflav(1),1,3) ! collinear
     &                 *pdf(sign2*iflav(2),2,3)
     &                 -pdf(sign1*iflav(1),1,3) ! real emission
     &                 *pdf(sign2*iflav(2),2,1)*sub(3))*C2
     &                 *NCmatrixelt(1,3,2)
               endif

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
               NCmatrixelt(0,2,2)=udwaab(1,2)
               NCmatrixelt(1,2,2)=duwaab(1,2)
               NCmatrixelt(0,3,2)=udwaab(1,3)
               NCmatrixelt(1,3,2)=duwaab(1,3)
            else
               NCmatrixelt(0,2,2)=udwaab(2,2)
               NCmatrixelt(1,2,2)=duwaab(2,2)
               NCmatrixelt(0,3,2)=udwaab(2,3)
               NCmatrixelt(1,3,2)=duwaab(2,3)
            endif
#endif /* WITH_NLO */
         endif
            
         iflav(j-1) = 0
         do if1=1,(nfl/2)*2,2
            if (bos.eq.311) then
               iflav(4-j)=(if1+PhysToDiag(3)-1)*fsign(physToDiag(4-j))
               iflav(3)  =(if1-PhysToDiag(3)+2)*fsign(physToDiag(3))
               k=FL_VVVg(iflav,j)
               if (lok(1)) then
                  res(k,1)=pdf(sign1*iflav(1),1,1)
     &                 *pdf(sign2*iflav(2),2,1)
     &                 *NCmatrixelt(0,2,1)*polcolg
               else
                  res(k,1) = 0
               endif
            elseif (bos.eq.411) then
               iflav(4-j)=(if1-PhysToDiag(3)+2)*fsign(physToDiag(4-j))
               iflav(3)  =(if1+PhysToDiag(3)-1)*fsign(physToDiag(3))
               k=FL_VVVg(iflav,j)
               if (lok(1)) then
                  res(k,1)=pdf(sign1*iflav(1),1,1)
     &                 *pdf(sign2*iflav(2),2,1)
     &                 *NCmatrixelt(1,2,1)*polcolg
               else
                  res(k,1) = 0
               endif
            endif
               
c--------------------------------------------------------------
            if ( nlo.eq.1 .and. (lok(2).or.lok(3)) ) then
#ifdef WITH_NLO                  
               if (bos.eq.311) then
                  res(k,j) = ( pdf(sign1*iflav(1),1,j) ! collinear
     &                 *pdf(sign2*iflav(2),2,j)
     &                 -pdf(sign1*iflav(1),1,2*j-3) ! real emission
     &                 *pdf(sign2*iflav(2),2,4-j )*sub(j))
     &                 *NCmatrixelt(0,j,2)
               elseif (bos.eq.411) then
                  res(k,j) = ( pdf(sign1*iflav(1),1,j) ! collinear
     &                 *pdf(sign2*iflav(2),2,j)
     &                 -pdf(sign1*iflav(1),1,2*j-3) ! real emission
     &                 *pdf(sign2*iflav(2),2,4-j )*sub(j))
     &                 *NCmatrixelt(1,j,2)
               endif
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
c   end subroutine qq_WWAj
c
c*****************************************************************

