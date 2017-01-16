*****************************************************************************
c
c    begin function m2s_qqZAA
c
c*****************************************************************************
c
c     michael rauch, <rauch@particle.physik.uni-karlsruhe.de>
c     Adapted from m2s_qqZZA
c     Last modified Nov 2010 (giuseppe)
c     
c     this function computes |ME|^2 for qq -> e- e+ A A ("ZAA") 
c     at LO and virtual corrections via call of 'qqZAA.f'
c
c ------------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION m2s_qqZAA(
     &      bos,  !in:  XXX Boson identifier,X=1 photon, X=2 Z, X=3,4 W+-
     &      nlo,  !in:  Born only for nlo=0, +interference for 1
     &      xi,   !in:  Feynman x parameters of incoming fermions
     &      p,    !in:  momenta associated with external fermions
     &      v,    !in:  momenta associated with gauge boson dec. fermions
     &      rn,   !in:  random number array
     &      ps_number
     &                        )
      IMPLICIT NONE

c declare global variables

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/lha.inc"
#include "tensor.inc"

c declare input/output variables
      INTEGER bos, nlo, ps_number
      DOUBLE PRECISION xi(nx), p(0:3,max_p), v(0:3,n_v,3), rn(1), m2s
      DOUBLE PRECISION q34(0:4,3) !check if I really need these

c declare external functions
      DOUBLE PRECISION dotrr
      DOUBLE COMPLEX dotrc
      INTEGER FL_VVV
      EXTERNAL dotrr, dotrc
      EXTERNAL FL_VVV

c  helicity selection
      INTEGER h
      COMMON /hcount / h

      COMPLEX*16 zero
      PARAMETER (zero=(0d0,0d0))

      DOUBLE PRECISION q_sf, qq
      INTEGER init/0/, I, mu
      SAVE init

c declare variables for summation process
      INTEGER ires

      INTEGER FSIGN(2+n_v), physToDiag(2), nmin, nmax
      SAVE nmax

c store contributions from subprocess ip in res(ip) where
      DOUBLE PRECISION res(maxnumsubproc)

c declare parton distribution variables
      DOUBLE PRECISION pdf(-6:6,2)

c define program switches
      LOGICAL ldebug
      PARAMETER (ldebug=.false.)

c variables for leshouches interface
      DOUBLE PRECISION weight,rnumb,RandomNumber

c ------------- BEGIN CODE -------------------

      if ( init .eq. 0 ) then
         
         write(6,*) " "
         if (bos .eq. 211) then
           write(6,*) "ZAA (Z->l+l-) amplitude square information:"
         else if (bos .eq. 511) then
           write(6,*) "ZAA (Z->vlvl~) amplitude square information:"
         else 
           write(6,*) "Unknown boson id in m2s_qqZAA.F: ", bos
           stop
         endif
         write(6,*) "-----------------------------------------------"
         write(6,*) " "
         call vtozaa_reset

         call printnfl(.false.)
         h = 1
         init = 1
         do i = 1,maxnumsubproc
            res(i) = 0
         enddo
      else
         do i = 1,nmax
            res(i) = 0
         enddo
      endif
      
      do mu = 0,3
         qze(mu,1) = v(mu,1,1)+v(mu,2,1)  ! A,Z(e)
         qa1(mu,1) = v(mu,3,1)            ! A1
         qa2(mu,1) = v(mu,4,1)            ! A2
         qaa(mu,1) = qa1(mu,1) + qa2(mu,1)
         qzaa(mu,1) = qze(mu,1) + qa1(mu,1) + qa2(mu,1)
      enddo

      qze(4,1) = dotrr(qze(0,1),qze(0,1))
      qa1(4,1) = dotrr(qa1(0,1),qa1(0,1))
      qa2(4,1) = dotrr(qa2(0,1),qa2(0,1))
      qaa(4,1) = dotrr(qaa(0,1),qaa(0,1))
      qzaa(4,1) = dotrr(qzaa(0,1),qzaa(0,1))

c for ZAA->2lAA precalculate leptonic tensors
c lepton spinors and Z, A polarisation vectors
c select helicity: h = random number for Z,photon helicity (h=1:2)

      if (bos.eq.211) then
        h = mod(h,8) + 1
      elseif (bos.eq.511) then 
        h = mod(h,4) + 1 ! only left-handed neutrinos
      endif
      ie = sign(1,4-h)
      iu = (-1)**(h+1)
      it = (-1)**((h-1)/2)

c...Les Houches interface
      if ((lha.or.hepmc).and..not.doNLO) then
         helicity(1)=-ie
         helicity(2)= ie
         helicity(3)= iu
         helicity(4)= it
      endif
      
      CALL OXXXXX(v(0,1,1),ZERO,-ie,1,zem(1,1))  !e- (v)
      CALL IXXXXX(v(0,2,1),ZERO,+ie,-1,zep(1,1)) !e+ (v~)
      CALL VXXXXX(qa1(0,1),ZERO,iu,1,atau1(1,1)) !A1
      CALL VXXXXX(qa2(0,1),ZERO,it,1,atau2(1,1)) !A2
      if (bos.eq.211) then
         CALL JIOXXX(zep(1,1),zem(1,1),GAL,ZERO,ZERO, ae(1,1)) !A*e
         CALL JIOXXX(zep(1,1),zem(1,1),GZL,ZMASS,ZWIDTH, ze(1,1)) !Z*e
      elseif (bos.eq.511) then
         CALL JIOXXX(zep(1,1),zem(1,1),GZN,ZMASS,ZWIDTH, ze(1,1)) !Z*v
      endif
      
c************ reduction of pentagon contribution *************
#ifdef WITH_NLO
      if (abs(nlo).eq.5 .or. nlo.eq.-4) then
c     reduction to "true" pentagons (see paper)
         qq = qze(0,1)*qaa(0,1) - qze(1,1)*qaa(1,1)-
     &        qze(2,1)*qaa(2,1) - qze(3,1)*qaa(3,1)
         xze = dotrc(qaa(0,1),ze(1,1))/qq
         xae = dotrc(qaa(0,1),ae(1,1))/qq

c------------------------------------------------------------
c         xze = (0.0d0,0.0d0)
c         xae = (0.0d0,0.0d0)
c------------------------------------------------------------

         do mu = 0,3
            zepp(mu+1,1) = ze(mu+1,1) - xze*qze(mu,1)
            aep(mu+1,1) = ae(mu+1,1) - xae*qze(mu,1)
            atau1p(mu+1,1) = atau1(mu+1,1)
            atau2p(mu+1,1) = atau2(mu+1,1)
         enddo
         do mu = 5,6
            zepp(mu,1) = ze(mu,1)
            aep(mu,1) = ae(mu,1)
            atau1p(mu,1) = atau1(mu,1)
            atau2p(mu,1) = atau2(mu,1)
         enddo         

      endif
#endif


      if (with_anom) then ! anomalous gauge boson couplings
c        using global form factor for all tensors of one phase space point
c        this ensures proper cancellations for anomalous contributions
c        energy scale is invariant ZAA mass
         call anomal_formfactor(qzaa(0,1),qze(0,1),qa1(0,1),qa2(0,1))


c leptonic tensors, last entry defines which photon is considered
c (1,2: one of the two photons radiated from the quark line)
c (0: both photons emitted from the final state leptons)
c only needed if bos=211 (e-e+aa)

         if (bos.eq.211) then
            call zaa_atoza_anomal(v(0,1,1),1,h,aaze1(0,1),1) !A -> Z A -> e+ e- a
            call zaa_atoza_anomal(v(0,1,1),1,h,aaze2(0,1),2) !A -> Z A -> e+ e- a
            call zaa_ztoza_anomal(v(0,1,1),1,h,zzae1(0,1),1) !Z -> Z A -> e+ e- a
            call zaa_ztoza_anomal(v(0,1,1),1,h,zzae2(0,1),2) !Z -> Z A -> e+ e- a
            call zaa_atozaa_anomal(v(0,1,1),1,h,azaa(0,1),0) !A -> e+ e- a a
            call zaa_ztozaa_anomal(v(0,1,1),1,h,zzaa(0,1),0) !Z -> e+ e- a a
         elseif (bos.eq.511) then
            call zaa_atoznaa_anomal(v(0,1,1),1,h,azaa(0,1),0) !A -> ve ve~ a a
            call zaa_ztoznaa_anomal(v(0,1,1),1,h,zzaa(0,1),0) !Z -> ve ve~ a a
         endif


      else !SM


c leptonic tensors, last entry defines which photon is considered
c (1,2: one of the two photons radiated from the quark line)
c (0: both photons emitted from the final state leptons)
c only needed if bos=211 (e-e+aa)

         if (bos.eq.211) then
            call zaa_atoza(v(0,1,1),1,h,aaze1(0,1),1) !A -> Z A -> e+ e- a
            call zaa_atoza(v(0,1,1),1,h,aaze2(0,1),2) !A -> Z A -> e+ e- a
            call zaa_ztoza(v(0,1,1),1,h,zzae1(0,1),1) !Z -> Z A -> e+ e- a
            call zaa_ztoza(v(0,1,1),1,h,zzae2(0,1),2) !Z -> Z A -> e+ e- a
            call zaa_atozaa(v(0,1,1),1,h,azaa(0,1),0) !A -> e+ e- a a
            call zaa_ztozaa(v(0,1,1),1,h,zzaa(0,1),0) !Z -> e+ e- a a
         endif

      endif !anom or SM

      
c if desired, print output information

      if ( ldebug ) then
 10      format( " p(", i1, ") = ", 4(f10.3, 2x) )
 20      format( " v(", i1, ") = ", 4(f10.3, 2x) )
c
         write(6,*) " "
         do i = 1, n_v
            write(6,20) i, v(0,i,1), v(1,i,1), v(2,i,1), v(3,i,1)
         end do
         print*," parton momenta "
         do i = 1, 2
            write(6,10) i, p(0,i), p(1,i), p(2,i), p(3,i)
         end do
c
         write(6,*) " "
         write(6,*) "xi(1), xi(2) =", xi
         write(6,*) "mu_f =",sqrt(mufsq(1,1))
      end if
c
c call PDF subroutine in order to determine parton
c distributions in the incoming (anti)protons.
c
      q_sf = sqrt(mufsq(1,1))

      call pdfproton( xi(1), q_sf, pdf(-6,1) )        ! f_a(x1)
      call pdfproton( xi(2), q_sf, pdf(-6,2) )        ! f_a(x2)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     Now sum over the subprocesses contributing to Z production

      nmax = FL_VVV(fsign,-1) !reset counter for subprocesses to 0

C*******************  q1 qbar1 ---> Z A A   **************************

c   physToDiag(ext.momentum label) = Feynman diagram label

      physToDiag(1)=1     !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=2

      fsign(1) = 1
      fsign(2) = -1
      fsign(3) = 1
      fsign(4) = -1
      fsign(5) = 1
      fsign(6) = 1

      call qq_ZAA(xi,p,v,physToDiag,fsign,bos,nlo,pdf,
     &     res,nmin,nmax,ps_number)

C*******************  q1bar q1 ---> Z A A   **********************

      physToDiag(1)=2     !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=1

      fsign(1) = 1
      fsign(2) = -1

      call qq_ZAA(xi,p,v,physToDiag,fsign,bos,nlo,pdf,
     &     res,nmin,nmax,ps_number)

c*****************  end of process evaluation  *********************

      if (init.eq.1) then
         init = init+1
      endif

      m2s = 0
      DO IRES = 1,NMAX
         m2s = m2s + RES(IRES)
      ENDDO

c...Les Houches interface - the most propable subprocess   
      if ((lha.or.hepmc).and..not.doNLO) then
         i=0
         weight=0.d0
         rnumb=RandomNumber()
         do while((i.le.nmax).and.(weight.le.rnumb*m2s))
            i=i+1
            weight=weight+res(i)
            iprocess=i
         enddo
      endif

      if (bos.eq.211) then
        m2s = m2s * 8d0    ! factor 8 for random helicity summation
      elseif (bos.eq.511) then
        m2s = m2s * 4d0    ! factor 4 for random helicity summation
      else 
        m2s = 0d0
      endif
      m2s_qqZAA = m2s
      
      RETURN
      END

c*****************************************************************************
c
c    end function m2s_qqZAA
c
c*****************************************************************************

c******************************************************************
c
c   begin subroutine qq_ZAA
c
c*****************************************************************

      subroutine qq_ZAA(xi,p,v,physToDiag,fsign,bos, nlo,
     1                  pdf,res,nmin,nmax,ps_number)
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"

      double precision p(0:3,n_p), v(0:3,n_v), xi(nx)
      double precision pdf(-6:6,2), res(maxnumsubproc)
      integer physToDiag(2), fsign(2+n_v), bos, nlo, nmin, nmax
      integer ps_number

c  INPUT:  p(0:3,4)        external physical parton momenta
c          v(0:3,nv)       leptons and photon momenta
c          physToDiag(4)   physToDiag(ext.mom. label) = Feynman diagram label
c          fsign           sign factors for amplitude calls; see qqZqq
c          nlo             nlo = 0: calculate LO only
c                          nlo = 1: calculate full NLO (Born + interference)
c  OUTPUT: res(k)
c
c  In and Output:
c          nmin, nmax      range of process numbers for this call
c
      double precision N                              ! color factors
      parameter(N=3d0)
c
c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/coupl.inc"

c local variables
      double precision uuzaa, ddzaa, uuzaa2, ddzaa2
      double precision pbar(0:3,2+n_v), dotrr, polcol
      double precision NCmatrixelt(0:1,0:1)
      double complex zero
      integer iflav(4), diagToPhys(2), FL_VVV
      external dotrr, FL_VVV

      logical ldebug
      parameter (ldebug=.false.)
      parameter (zero = (0d0,0d0) )
      integer i,if1,k,mu

      uuzaa = 0d0
      uuzaa2 = 0d0
      ddzaa = 0d0
      ddzaa2 = 0d0

      nmin = nmax+1
      do i = 1,2
         diagToPhys(physToDiag(i)) = i
      enddo

      polcol = 1d0/(4d0*N**2*xi(1)*xi(2))

      do mu = 0,3
         do i = 1,2
            pbar(mu,physToDiag(i))=p(mu,i)
         enddo
      enddo

      do mu = 0,3             ! kinematics for ZAA
         pbar(mu,3) = v(mu,1) ! e-
         pbar(mu,4) = v(mu,2) ! e+
         pbar(mu,5) = v(mu,3) ! A1
         pbar(mu,6) = v(mu,4) ! A2
      enddo

#ifdef WITH_MADGRAPH
      if (ldebug) then
         if (bos.eq.211) then
            call suu_zaa(pbar,fsign,uuzaa2)
            call sdd_zaa(pbar,fsign,ddzaa2)
         elseif (bos.eq.511) then
            call suu_znaa(pbar,fsign,uuzaa2)
            call sdd_znaa(pbar,fsign,ddzaa2)
         endif
      endif
#endif

      call QQZAA(bos,pbar,fsign,nlo,1,uuzaa,ddzaa,ps_number)

#ifdef WITH_MADGRAPH
      if (ldebug) then
         print*,"uuzaa =   ",uuzaa
         print*,"uuzaa2 =  ",uuzaa2
         print*,"ddzaa =   ",ddzaa
         print*,"ddzaa2 =  ",ddzaa2
         read(*,*)
      endif
#endif

      NCmatrixelt(0,0)=uuzaa
      NCmatrixelt(1,1)=ddzaa

      do if1=1,nfl
            iflav(1)=if1*fsign(physToDiag(1))
            iflav(2)=if1*fsign(physToDiag(2))
            k=FL_VVV(iflav,1)

            res(k)=pdf(sign1*iflav(1),1)
     &            *pdf(sign2*iflav(2),2)
     &            *NCmatrixelt(mod(if1,2),mod(if1,2))*polcol

      enddo 

      nmax = FL_VVV(iflav,-2)

      end
c******************************************************************
c
c   end subroutine qq_ZAA
c
c*****************************************************************


      subroutine vtozaa_reset
      implicit none
#include "tensor.inc"
      complex*16 zero
      parameter (zero=(0d0,0d0))
      integer j,mu

      do mu = 0,5
         do j = 1,3
            aaze1(mu,j) = zero
            zzae1(mu,j) = zero
            aaze2(mu,j) = zero
            zzae2(mu,j) = zero
            azaa(mu,j) = zero
            zzaa(mu,j) = zero
         enddo
      enddo
      print*," precalculated tensors initialized to ",zero 
      return
      end
