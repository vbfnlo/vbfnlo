*****************************************************************************
c
c    begin function m2s_qqWH
c
c*****************************************************************************
c
c       Adapted from m2s_qqWAA by giuseppe bozzi and m2s_qqWA by J. Bellm
c
c       this function computes |ME|^2 for qq-> ve(ve~) e+(-) H ("W+- H")
c	     at LO and virtual corrections via call of 'qqWH.f'
c
c       Robin Roth <robin@particle.physik.uni-karlsruhe.de>
c       Initial version: November 2012
c       Last modified: November 2012
c
c ------------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION m2s_qqWH(
     &      bos,  !in:  XXX Boson identifier,X=1 photon, X=2 Z, X=3,4 W+-, X=6 H
     &      nlo,  !in:  Born only for nlo=0, +interference for 1
     &      xi,   !in:  Feynman x parameters of incoming fermions
     &      p,    !in:  momenta associated with external fermions
     &      v,    !in:  momenta associated with gauge boson dec. fermions
     &      rn    !in:  random number array
     &                        )

      use globalvars, only: ldoblha

      IMPLICIT NONE


c v(mu,1,1), v(mu,2,1) : w+-
c v(mu,3,1), v(mu,4,1) : H decay products (in H-> 2l)
! in H-> 4l also v 5 and v 6


c declare global variables

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/lha.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

c declare input/output variables

      INTEGER bos, nlo !, NW, ps_number
      DOUBLE PRECISION xi(nx), p(0:3,max_p), v(0:3,max_v,max_kin),  ! TODO: use n_v n_p
     1                 rn(1), m2s

c declare external functions

      DOUBLE PRECISION dotrr
      DOUBLE COMPLEX dotrc
      INTEGER FL_VV
      EXTERNAL dotrr, dotrc, FL_VV

c declare local variables

      DOUBLE PRECISION zero
      PARAMETER (zero=0d0)

      DOUBLE PRECISION q_sf 
      INTEGER init/0/, I, mu 
      SAVE init

c declare variables for summation process

      INTEGER ires

      INTEGER FSIGN(4), physToDiag(2), nmin, nmax
      SAVE nmax

c store contributions from subprocess ip in res(ip) where

      DOUBLE PRECISION res(8) ! 8 enough

c declare parton distribution variables

      DOUBLE PRECISION pdf(-6:6,2)
c pdfckm holds pdfs modified with ckm      
      DOUBLE PRECISION pdfckm(-6:6,2)

c switch between W and WH
      logical hashiggs      

c define program switches

      LOGICAL ldebug
      PARAMETER (ldebug=.false.)

c variables for leshouches interface
      DOUBLE PRECISION weight,rnumb,RandomNumber

      
      DOUBLE PRECISION CLR,XM2,XMG,B
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6)
      double precision fac

      double precision mass2
      external mass2
c ------------- BEGIN CODE -------------------

      if (ldebug) then      
         write(6,*) "------------from m2s_qqWH.F----------------------"
         write(6,*) " "
      endif

      if (bos == 36 .or. bos == 46) then
         hashiggs = .true.
      else
         hashiggs = .false.
      endif

      if ( init .eq. 0 .and. .not. ldoblha) then
         write(6,*) " "
         if (bos.eq.36) then
            write(6,*) "W+H amplitude square information"
         elseif(bos.eq.46) then
            write(6,*) "W-H amplitude square information:"
         elseif (bos.eq.3) then
            write(6,*) "W+ amplitude square information"
         elseif(bos.eq.4) then
            write(6,*) "W- amplitude square information:"
         else
            write(6,*) "wrong call of m2s_qqWH.F, bos = ", bos
         endif

         call printnfl(.true.)
c        h = 1
         init = 1
      endif

      res(:) = 0
      call wtowh_reset

c...Les Houches interface
        if ((lha.or.hepmc).and..not.doNLO) then
           helicity(1)=-1
           helicity(2)= 1
        endif

      if ( ldebug ) then
         write(6,*) " m2s_qqWH.F, part1 -------------------------- "
c 10      format( " p(", i1, ") = ", 4(f10.3, 2x) )
c 20      format( " v(", i1, ") = ", 4(f10.3, 2x) )
c
         write(6,*) " "
         do i = 1, n_v
            write(6,20) i, v(0,i,1), v(1,i,1), v(2,i,1), v(3,i,1)
         end do
         print*," parton momenta "
         do i = 1, n_p !max_p -1
            write(6,10) i, p(0,i), p(1,i), p(2,i), p(3,i)
         end do
c
         write(6,*) " "
         write(6,*) "xi(1), xi(2) =", xi
         write(6,*) "mu_f =",sqrt(mufsq(1,1))
      end if !ldebug

      if(hashiggs) then
         qw0(0:3,1) = v(0:3,1,1)+v(0:3,2,1)   !W+-
         qh(0:3,1) = 0d0
         do i=3,n_v ! sum all higgs decay products
            qh(0:3,1) = qh(0:3,1) + v(0:3,i,1)
         enddo
         qwh(0:3,1) = qh(0:3,1)+qw0(0:3,1)
      else
         qwh(0:3,1) = v(0:3,1,1)+v(0:3,2,1)   !W+- 
      endif

      if(hashiggs) then
         qw0(4,1) = dotrr(qw0(0,1),qw0(0,1))
         qh(4,1) = dotrr(qh(0,1),qh(0,1))
      endif
      qwh(4,1) = dotrr(qwh(0,1),qwh(0,1))

      ! lepton spinors (later combinded to boson pol. vectors

      if (bos.eq.3 .or. bos.eq.36) then
         CALL OXXXXX(v(0,1,1),ZERO ,-1,1,wve(1,1))           !ve 
         CALL IXXXXX(v(0,2,1),ZERO ,1,-1,wep(1,1))           !e+     
      elseif ((bos.eq.4 .or. bos.eq.46) .and. .not. ldoblha) then
         CALL OXXXXX(v(0,2,1),ZERO ,-1,1,wep(1,1))            !e-    
         CALL IXXXXX(v(0,1,1),ZERO ,1,-1,wve(1,1))          !ve~  
      elseif ((bos.eq.4 .or. bos.eq.46) .and. ldoblha) then
         CALL OXXXXX(v(0,1,1),ZERO ,-1,1,wep(1,1))            !e-    
         CALL IXXXXX(v(0,2,1),ZERO ,1,-1,wve(1,1))          !ve~  
      endif
      
      if(hashiggs) then
         CALL SXXXXX(qh(0,1),1,hig(1,1))                     !H
      endif

      if (ldebug) then
         print*, 'hashiggs', hashiggs
         print*, 'hig', hig(:,1)
      endif


! from amplitudes/vvjj/qqhqq.F      
c
      if(hashiggs) then
         fac = 1.d0 / ( (qh(4,1)-xm2(6))**2 + xmg(6)**2 ) !propagator^2 only
         fac = fac *16*pi*xmg(6)
      else
         fac = 1d0
      endif
c
c multiply by 16 pi*M Gamma to compensate for the phase space integral in
c LIPSN 
! (here LISVHj)
c
c   Int dq^2/2pi * dOmega/(32 pi^2) * 1/(q^2-M*2)^2+(MGamma)^2) 
c
c   = 1/(16 pi MGamma)
!   
!   assuming m**2 >> MGamma
c
c (where the propagator factor is fac above)


      if(ldebug) print*, 'anom', with_anom

      if(hashiggs) then
         if (with_anom) then ! anomalous gauge boson couplings
            if (bos.eq.36) then
               call wpmtowh_an(v(0,1,1),1,wtowh(0,1),3)  !W+ -> ve e+ h
            elseif (bos.eq.46) then
               call wpmtowh_an(v(0,1,1),1,wtowh(0,1),4)  !W- -> ve e- h
            endif
         else !SM
            if (bos.eq.36) then
               call wpmtowh(v(0,1,1),1,wtowh(0,1),3)  !W+ -> ve e+ h
            elseif (bos.eq.46) then
               call wpmtowh(v(0,1,1),1,wtowh(0,1),4)  !W- -> ve e- h
            endif
         endif !anom or SM
      else !nohiggs
         if (bos.eq.3) then
            CALL JIOXXX(wep(1,1), wve(1,1), GWF, WMASS, WWIDTH, wtowh(0,1)) 
            ! this should better be called w0
         elseif (bos.eq.4) then
            CALL JIOXXX(wve(1,1), wep(1,1), GWF, WMASS, WWIDTH, wtowh(0,1)) 
         endif
      endif !hashiggs

      if (ldoblha) then
         physToDiag(1:2)=blha_physdiagmap(1:2,blha_cursubproc,blha_curproc)
         fsign(1:4) = blha_fsign(1:4,blha_cursubproc,blha_curproc)
         ! if (bos.eq.4 .or. bos.eq.46) then
         !    fsign(3:4)=-fsign(3:4)
         ! endif
         if (ldebug) then
            print*, 'blha in m2s_qqWH'
            print*, physToDiag(1:2)
            print*, 'bos', bos
            print*, 'fsign given', fsign(1:4)
            print*, 'fsign expt ', 1, -1, 1, -1
            print*, 'p1', p(0:3, 1)
            print*, 'p2', p(0:3, 2)
            print*, 'p3', p(0:3, 3)
            print*, 'v1', v(0:3, 1, 1), mass2(v(0:3, 1, 1))
            print*, 'v2', v(0:3, 2, 1), mass2(v(0:3, 2, 1))
            print*, 'v3', v(0:3, 3, 1), mass2(v(0:3, 3, 1))
            print*, 'v4', v(0:3, 4, 1), mass2(v(0:3, 4, 1))
            print*, 'np, nv', n_p, n_v
            print*, 'hmass, hwidth', hmass, hwidth
            print*, 'hmass, hwidth', sqrt(xm2(6)), sqrt(xmg(6))
            print*, 'wmass, wwidth', wmass, wwidth
         endif

         call qq_WH(xi,p,v,physToDiag,fsign,bos,nlo,pdf,pdfckm,res,nmin,nmax)
         return

      endif

      RES(:)=0d0

      if ( ldebug ) then
         write(6,*) " m2s_qqWH.F, part2 -------------------------- "
 10      format( " p(", i1, ") = ", 4(f10.3, 2x) )
 20      format( " v(", i1, ") = ", 4(f10.3, 2x) )
c
         write(6,*) " "
         do i = 1, n_v
            write(6,20) i, v(0,i,1), v(1,i,1), v(2,i,1), v(3,i,1)
         end do
         print*," parton momenta "
         do i = 1, max_p-1
            write(6,10) i, p(0,i), p(1,i), p(2,i), p(3,i)
         end do
c
         write(6,*) " "
         write(6,*) "xi(1), xi(2) =", xi
         write(6,*) "mu_f =",sqrt(mufsq(1,1))
      end if !ldebug


c       call PDF subroutine in order to determine parton
c       distributions in the incoming (anti)protons.

      q_sf = sqrt(mufsq(1,1))

      call pdfproton( xi(1), q_sf, pdf(-6,1) )        ! f_a(x1)
      call pdfproton( xi(2), q_sf, pdf(-6,2) )        ! f_a(x2)

! ckm testing with special pdfs      
!     do i=1,2
!        pdf(-1,i)=pdf(1,i)
!        pdf(2,i)=pdf(1,i)
!        pdf(-2,i)=pdf(1,i)
!        pdf(3,i)=pdf(1,i)
!        pdf(-3,i)=pdf(1,i)
!        pdf(4,i)= 0d0 ! pdf(1,i) !0d0 !
!        pdf(-4,i)=0d0 ! pdf(1,i) !0d0 !
!     enddo

      call pdfckmtransform(pdf,pdfckm)
      if (ldebug) then
         print*, 'pdfs'
         print*, pdf
         print*, pdfckm
      endif


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     Now sum over the subprocesses contributing to Z production

      nmax = FL_VV(fsign,-1) !reset counter for subprocesses to 0

C*******************  q1 qbar2 ---> W H   **************************

c   physToDiag(ext.momentum label) = Feynman diagram label

      physToDiag(1)=1     !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=2

      fsign(1) = 1   ! u
      fsign(2) = -1  ! dbar
      fsign(3) = 1   ! ve
      fsign(4) = -1  ! e- = ebar
!      fsign(5) = 1   ! h

      ! in case of W- (instead of W+) WMH <> WPH
      ! inverse sign (particle/antiparticle) of
      ! u, d, ve, e (i.e. particles 1-4)
      if (bos.eq.4 .or. bos.eq.46) then
         do i=1,4
            fsign(i)=-fsign(i)
         enddo
      endif


      call qq_WH(xi,p,v,physToDiag,fsign,bos,nlo,pdf,pdfckm,res,nmin,nmax)

C*******************  q1bar q2 ---> W H   **********************

      physToDiag(1)=2     !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=1

c      fsign(1) = 1
c      fsign(2) = -1

      call qq_WH(xi,p,v,physToDiag,fsign,bos,nlo,pdf,pdfckm,res,nmin,nmax)

c*****************  end of process evaluation  *********************

      if (init.eq.1) then ! TODO: remove?
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

      ! multiply in factor from higgs propagator
      m2s_qqWH = m2s * fac

      RETURN
      END

c*************************************************************************
c
c    end function m2s_qqWH
c
c*************************************************************************

c******************************************************************
c
c   begin subroutine qq_WH
c
c*****************************************************************

      subroutine qq_WH(xi,p,v,physToDiag,fsign,bos, nlo,
     1                  pdf,pdfckm,res,nmin,nmax)

      use globalvars, only: ldoblha

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

      real*8 p(0:3,max_p-1), v(0:3,n_v), xi(nx)
      real*8 pdf(-6:6,2), pdfckm(-6:6,2), res(8) !res(maxnumsubproc) too big
      integer physToDiag(2), fsign(2+n_v), bos, nlo, nmin, nmax

c wbf_zh calls the amplitude square routines 
c             qqWA    for qq --> VA --> 2 leptons + A
c
c  INPUT:  p(0:3,4)        external physical parton momenta
c          
c          v(0:3,nv)       leptons and higgs decay products
c                          the sum q(mu) = v(mu,1)+...+v(mu,nv) is needed
c                          v(mu,1): ve
c                          v(mu,2): e
c                          v(mu,3), v(mu,4): Higgs decay products
c
c
c          physToDiag(4)   physToDiag(ext.mom. label) = Feynman diagram label
c          fsign           sign factors for amplitude calls; see qqZqq
c          nlo             nlo = 0: calculate LO only
c                          nlo = 1: calculate full NLO (Born + interference)
c  OUTPUT:
c          uuvv(ID)   output in format as in qqZqq(j), but corresponding
c          etc.            to res(*) in m2s_qqZqq
c     or   res(k)   
c
c  In and Output
c    nmin, nmax            range of process numbers for this call
c
      double precision N                              ! color factors
      parameter(N=3d0)
c
c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/coupl.inc"

c local variables
      real*8 udwh  , udwh2
      real*8 pbar(0:3,2+n_v), dotrr, polcol ! 5 = 2in + l v h
      double precision getpdfckm
      double precision NCmatrixelt
      double complex zero
      integer iflav(4), diagToPhys(2), FL_VV
      external dotrr, FL_VV

c      logical ChargedCurrent, sametype, oneAntiparticle
      logical ldebug
      parameter (ldebug=.false.)
      parameter (zero = (0d0,0d0) )
      integer i,if1,k,mu

      udwh = 0d0 
      udwh2 = 0d0 

      nmin = nmax+1
      do i = 1,2
         diagToPhys(physToDiag(i)) = i
      enddo

      polcol = 1d0/(4d0 * N**2 * xi(1) * xi(2))

      do mu = 0,3
         do i = 1,2
            pbar(mu,physToDiag(i))=p(mu,i)
         enddo
      enddo


      !TODO: all not needed
      do mu = 0,3             ! kinematics for WH --> 2 lepton + H decay
         pbar(mu,3) = v(mu,1) ! ve
         pbar(mu,4) = v(mu,2) ! e+
!         pbar(mu,5) = v(mu,3)+v(mu,4) ! H 
      enddo

#ifdef WITH_MADGRAPH
      if(ldebug) then
         if (bos.eq.36) then
            call sudwph(pbar,udwh2)
         else if (bos.eq.46) then
            call sudwmh(pbar,udwh2)
         endif
      endif
#endif

      call QQWH(pbar,fsign,nlo,1,udwh, bos)
      if (ldebug) then
         print*, 'in qq_WH'
         print*, 'pbar', pbar
         print*, 'fsign', fsign
         print*, 'udwh', udwh
      endif


      if (ldoblha) then
         return
      endif

#ifdef WITH_MADGRAPH
      if (ldebug) then
            print*,"-----------------------------------------"
            print*,"m2s_qqWH: comparison with madgraph in qq_WH"
            print*,"-----------------------------------------"
            print*, abs(udwh/udwh2 -1)
c         if(abs(udwh/udwh2 -1).ge.1.0d-10) then !was 1.0d-6
            print*,"udwh =   ",udwh
            print*,"udwh2 =  ",udwh2
c         endif
      endif
#endif


      NCmatrixelt=udwh


c loop over all possible quark flavors (ignoring b)
      do if1=2,(nfl/2)*2,2


c bos 36 or 46 included in fsign         
c see just before call to qq_WH function in this file

c if1 = 2,4 (more or less depending on nfl)
c iflav(1,2) actual incoming quark flavors (e.g -1,2 ...)
           iflav(1)=( if1 - physToDiag(1) + 1 )*fsign(physToDiag(1))
           iflav(2)=( if1 - physToDiag(2) + 1 )*fsign(physToDiag(2))


c originally taken:
c but now 36 and 46 combined           
c from m2s_qqWA:         
c        if (bos.eq.3) then
c          iflav(1)=( if1 - physToDiag(1) + 1 )*fsign(physToDiag(1))
c          iflav(2)=( if1 - physToDiag(2) + 1 )*fsign(physToDiag(2))
c        else if (bos.eq.4) then
c           iflav(1)=( if1 + physToDiag(1) - 2 )*fsign(physToDiag(1))
c           iflav(2)=( if1 + physToDiag(2) - 2 )*fsign(physToDiag(2))
c           iflav(1)= -( if1 - physToDiag(1) + 1 )*fsign(physToDiag(1))
c           iflav(2)= -( if1 - physToDiag(2) + 1 )*fsign(physToDiag(2))
c        endif

c from m2s_qqWAA:
c         if (bos.eq.3) then
c            iflav(1)=( if1 - physToDiag(1) + 2 )*fsign(physToDiag(1))
c            iflav(2)=( if1 - physToDiag(2) - 1 )*fsign(physToDiag(2))
c         else if (bos.eq.4) then
c            iflav(1)=( if1 + physToDiag(1) - 1 )*fsign(physToDiag(1))
c            iflav(2)=( if1 + physToDiag(2) + 2 )*fsign(physToDiag(2))
c         endif

         if (ldebug) then
            print*,"iflav1 ",iflav(1)
            print*,"iflav2 ",iflav(2)
            print*,"pdf1 ",pdf(sign1*iflav(1),1)
            print*,"pdf2 ",pdf(sign2*iflav(2),2)
            print*,"pdf1ckm ",pdfckm(sign1*iflav(1),1)
            print*,"pdf2ckm ",pdfckm(sign2*iflav(2),2)
         endif

         k=FL_VV(iflav,1)

         if (ldebug) then
            print*,"k ",k
            print*,"sign1 ",sign1
            print*,"sign2 ",sign2
            print*,"fsign ",fsign
            print*,"Matrixelement ",NCmatrixelt
            print*,"polcol ",polcol

         endif
         !res(k)=getpdf(pdf,sign1*iflav(1),sign2*iflav(2))
         res(k)=getpdfckm(pdf,pdfckm,sign1*iflav(1),sign2*iflav(2))
     &         *NCmatrixelt*polcol

      enddo 

      nmax = FL_VV(iflav,-2)
      if (ldebug) then
         print*,"nmax ",nmax
         print*,"res ",res
      endif

      end
c******************************************************************
c
c   end subroutine qq_WH
c
c*****************************************************************


      subroutine wtowh_reset
      implicit none
#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
      complex*16 zero
      parameter (zero=(0d0,0d0))
      integer j,mu

      do mu = 0,3
         do j = 1,3
            wtowh(mu,j) = zero            
            wtoonlyw(mu,j) = zero            
         enddo
      enddo
c      print*," W to WH  tensors initialized to ",zero 
      return
      end



