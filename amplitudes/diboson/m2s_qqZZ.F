*****************************************************************************
c
c    begin function m2s_qqZZ
c
c*****************************************************************************
c
c       michael rauch, <rauch@particle.physik.uni-karlsruhe.de>
c       Adapted from m2s_qqAA
c       Last modified Jun 2011 (rauch)
c
c     this function computes |ME|^2 for qq -> Z Z 
c     at LO and virtual corrections via call of 'qqZZ.f'
c
c ------------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION m2s_qqZZ(
     &      bos,  !in:  XXX Boson identifier,X=1 photon, X=2 Z, X=3,4 W+-
     &      nlo,  !in:  Born only for nlo=0, +interference for 1
     &      xi,   !in:  Feynman x parameters of incoming fermions
     &      p,    !in:  momenta associated with external fermions
     &      v,    !in:  momenta associated with gauge boson dec. fermions
     &      rn    !in:  random number array
     &                        )
      IMPLICIT NONE

c declare global variables

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/process.inc"
#include "tensor.inc"

c declare input/output variables
      INTEGER bos, nlo
      DOUBLE PRECISION xi(nx), p(0:3,max_p), v(0:3,n_v,3), rn(1), m2s

c declare external functions
      DOUBLE PRECISION dotrr
      DOUBLE COMPLEX dotrc
      INTEGER FL_VV
      EXTERNAL dotrr, dotrc
      EXTERNAL FL_VV

c  helicity selection
      INTEGER h
      COMMON /hcount / h

      DOUBLE PRECISION zero
      PARAMETER (zero=0d0)

      DOUBLE PRECISION q_sf
      INTEGER init/0/, I, mu
      SAVE init

c declare variables for summation process
      INTEGER ires

      INTEGER FSIGN(2+n_v), physToDiag(2), nmin, nmax
      SAVE nmax

      integer N_gen_up, N_gen_down
      real*8 fac_Z_up, fac_Z_down, mjj2
      external mjj2

c store contributions from subprocess ip in res(ip) where
      DOUBLE PRECISION res(maxnumsubproc)

c declare parton distribution variables
      DOUBLE PRECISION pdf(-6:6,2)

c define program switches
      LOGICAL ldebug, ldebuglt
      PARAMETER (ldebug=.false.)
      PARAMETER (ldebuglt=.false.)

c variables for leshouches interface
      DOUBLE PRECISION weight,rnumb,RandomNumber

c temporary variables

      double complex zm2i(2:3)
      save zm2i
      double precision  clr, xm2, xmg, b
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6)

c ------------- BEGIN CODE -------------------

      if ( init .eq. 0 ) then
         
         write(6,*) " "
         write(6,*) "ZZ amplitude square information:"
         write(6,*) "-----------------------------------------------"
         write(6,*) " "

         call printnfl(.false.)
         h = 1
         init = 1
         do i = 1,maxnumsubproc
            res(i) = 0
         enddo
         zm2i(2) = 1d0/dcmplx(xm2(2),-xmg(2))
         zm2i(3) = 1d0/dcmplx(xm2(3),-xmg(3))
      else
         do i = 1,nmax
            res(i) = 0
         enddo
      endif
      
      do mu = 0,3
         qp(mu,1) = v(mu,1,1)+v(mu,2,1)   !Z_e
         qm(mu,1) = v(mu,3,1)+v(mu,4,1)   !Z_mu
         qaa(mu,1) = qp(mu,1)+qm(mu,1)    !Z_all
      enddo

      qp(4,1) = dotrr(qp(0,1),qp(0,1))
      qm(4,1) = dotrr(qm(0,1),qm(0,1))
      qaa(4,1) = dotrr(qaa(0,1),qaa(0,1))

c photon polarisation vectors
c select helicity: h = random number for photon helicity (h=1:2)
      h = mod(h,4) + 1
      ie = (-1)**(h+1)
      iu = (-1)**((h-1)/2)
      
      if ((lha.or.hepmc) .and. (.not. doNLO)) then
         helicity(1)=-ie
         helicity(2)= ie
         helicity(3)=-iu
         helicity(4)= iu
      end if

      CALL OXXXXX(v(0,1,1),ZERO ,-ie,+1,wve(1,1))                 !e- 
      CALL IXXXXX(v(0,2,1),ZERO ,+ie,-1,wep(1,1))                 !e+        
      CALL OXXXXX(v(0,3,1),ZERO ,-iu,+1,wmu(1,1))                 !mu-      
      CALL IXXXXX(v(0,4,1),ZERO ,+iu,-1,wvm(1,1))                 !mu+
      CALL JIOXXX(wep(1,1),wve(1,1),GZL,ZMASS,ZWIDTH, wp(1,1))    !Z_e
      CALL JIOXXX(wep(1,1),wve(1,1),GAL,ZERO ,ZERO  , atau1(1,1)) !A_e
      CALL JIOXXX(wvm(1,1),wmu(1,1),GZL,ZMASS,ZWIDTH, wm(1,1))    !Z_mu
      CALL JIOXXX(wvm(1,1),wmu(1,1),GAL,ZERO ,ZERO  , atau2(1,1)) !A_mu

      SELECT CASE(procid)
      CASE(ZZhad, GFZZhad)
        ! for general output to all combinations: up-type first
        if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
          ! set couplings for Z hadronic decay into uubar
          call setZtouu
        endif
        CALL JIOXXX(wep(1,1),wve(1,1),GZ_ZF,ZMASS,ZWIDTH, wp(1,1))    !Z_e
        CALL JIOXXX(wep(1,1),wve(1,1),GZ_AF,ZERO ,ZERO  , atau1(1,1)) !A_e
        call zz_ztozz_had(v(0,1,1),1,wz(1,1),0) !Z -> e+ e- mu+ mu-
        call zz_atozz_had(v(0,1,1),1,wa(1,1),0) !A -> e+ e- mu+ mu-
        if ((nlo .eq. -8) .or. (nlo.eq.1)) 
     &    call zz_htozz_had(v(0,1,1),1,w4l(h,1),0) !H -> e+ e- mu+ mu-
      CASE DEFAULT
        call zz_ztozz(v(0,1,1),1,wz(1,1),0) !Z -> e+ e- mu+ mu-
        call zz_atozz(v(0,1,1),1,wa(1,1),0) !A -> e+ e- mu+ mu-
        if ((nlo .eq. -8) .or. (nlo.eq.1)) 
     &    call zz_htozz(v(0,1,1),1,w4l(h,1),0) !H -> e+ e- mu+ mu-
      END SELECT

cc add k^mu k^nu term for Z
c      ctmp1 = -zm2i(2) * dotrc(qp(0,1),wp(1,1))
c      ctmp2 = -zm2i(2) * dotrc(qm(0,1),wm(1,1))
c      ctmp3 = -zm2i(2) * dotrc(qaa(0,1),wz(1,1))
c      do mu=0,3
c        wp(mu+1,1) = wp(mu+1,1) + qp(mu,1) * ctmp1
c        wm(mu+1,1) = wm(mu+1,1) + qm(mu,1) * ctmp2
c        wz(mu+1,1) = wz(mu+1,1) + qaa(mu,1) * ctmp3
c      enddo


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
         do i = 1, max_p-1
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
C     Now sum over the subprocesses contributing to Z Z production

      nmax = FL_VV(fsign,-1) !reset counter for subprocesses to 0

C*******************  q1 qbar1 ---> Z Z   **************************

c   physToDiag(ext.momentum label) = Feynman diagram label

      physToDiag(1)=1     !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=2

      fsign(1) = 1
      fsign(2) = -1
      fsign(3) = 1
      fsign(4) = -1
      fsign(5) = 1
      fsign(6) = -1

      call qq_ZZ(xi,p,v,physToDiag,fsign,bos,nlo,pdf,
     &     res,nmin,nmax)

C*******************  q1bar q1 ---> Z Z   **********************

      if (nlo .ne. -8) then
cc gg is symmetric -> would be double-counting
      physToDiag(1)=2     !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=1

      fsign(1) = 1
      fsign(2) = -1

      call qq_ZZ(xi,p,v,physToDiag,fsign,bos,nlo,pdf,
     &     res,nmin,nmax)
      endif

c**************  end of process evaluation part one ******************


! now run again for down-type decay products if 93 93 / 94 94 is requested.

      SELECT CASE(procid)
      CASE(ZZhad, GFZZhad)
       if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
        ! set couplings for Z hadronic decay into ddbar
        call setZtodd
        CALL JIOXXX(wep(1,1),wve(1,1),GZ_ZF,ZMASS,ZWIDTH, wp(1,1))    !Z_e
        CALL JIOXXX(wep(1,1),wve(1,1),GZ_AF,ZERO ,ZERO  , atau1(1,1)) !A_e
        call zz_ztozz_had(v(0,1,1),1,wz(1,1),0) !Z -> e+ e- mu+ mu-
        call zz_atozz_had(v(0,1,1),1,wa(1,1),0) !A -> e+ e- mu+ mu-
        if ((nlo .eq. -8) .or. (nlo.eq.1)) 
     &    call zz_htozz_had(v(0,1,1),1,w4l(h,1),0) !H -> e+ e- mu+ mu-


C*******************  q1 qbar1 ---> Z Z   **************************

c   physToDiag(ext.momentum label) = Feynman diagram label

      physToDiag(1)=1     !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=2

      fsign(1) = 1
      fsign(2) = -1
      fsign(3) = 1
      fsign(4) = -1
      fsign(5) = 1
      fsign(6) = -1

      call qq_ZZ(xi,p,v,physToDiag,fsign,bos,nlo,pdf,
     &     res,nmin,nmax)

C*******************  q1bar q1 ---> Z Z   **********************

      if (nlo .ne. -8) then
cc gg is symmetric -> would be double-counting
      physToDiag(1)=2     !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=1

      fsign(1) = 1
      fsign(2) = -1

      call qq_ZZ(xi,p,v,physToDiag,fsign,bos,nlo,pdf,
     &     res,nmin,nmax)
      endif

c***  end of extra process evaluation for hadronic decay  ***

       endif
      end select


c*****************  end of process evaluation  *********************


      SELECT CASE(procid)
      CASE(ZZhad,GFZZhad)
        call hadDecayFactor_Z(finalquarks(1), abs(mjj2(v(0,1,1),v(0,2,1))), N_gen_up, N_gen_down, fac_Z_up, fac_Z_down)
        DO IRES = 1,NMAX
          if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
             if (ires.le.(NMAX/2)) then
                res(IRES) = res(IRES) * fac_Z_up      ! factor fac_Z_up for up-type with all combinations
             elseif (ires.gt.(NMAX/2)) then
                res(IRES) = res(IRES) * fac_Z_down    ! factor fac_Z_down for down-type with all combinations
             endif
             if (mod(nmax,2).ne.0) then
               print*, "Something went wrong with the hadronic decays!"
               stop
             endif
          elseif (mod(abs(finalquarks(1)),2).eq.0) then
             res(IRES) = res(IRES) * fac_Z_up         ! up-type final states
          else
             res(IRES) = res(IRES) * fac_Z_down       ! down-type final states
          endif
        ENDDO
      END SELECT

      if (init.eq.1) then
         init = init+1
      endif

      m2s = 0
      DO IRES = 1,NMAX
         m2s = m2s + RES(IRES)
      ENDDO

c...Les Houches interface - the most propable subprocess   
      if ((lha.or.hepmc) .and. (.not. doNLO)) then
         i=0
         weight=0.d0
         rnumb=RandomNumber()
         do while((i.le.nmax).and.(weight.le.rnumb*m2s))
            i=i+1
            weight=weight+res(i)
            iprocess=i
         enddo
         SELECT CASE(process)
         CASE(ZZhad,GFZZhad)
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
      m2s = m2s * 4d0    ! factor 4 for random helicity summation
      m2s_qqZZ = m2s
      RETURN
      END
c*****************************************************************************
c
c    end function m2s_qqZZ
c
c*****************************************************************************
c******************************************************************
c
c   begin subroutine qq_ZZ
c
c*****************************************************************

      subroutine qq_ZZ(xi,p,v,physToDiag,fsign,bos, nlo,
     1                 pdf,res,nmin,nmax)
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"

      double precision p(0:3,n_p), v(0:3,n_v), xi(nx)
      double precision pdf(-6:6,2), res(maxnumsubproc)
      integer physToDiag(2), fsign(2+n_v), bos, nlo, nmin, nmax

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
      double precision uuzz, ddzz, uuzz2, ddzz2
      double precision pbar(0:3,2+n_v), dotrr, polcol, polcolg
      double precision NCmatrixelt(0:1,0:1)
      double complex zero
      integer iflav(4), diagToPhys(2), FL_VV
      external dotrr, FL_VV

      logical ldebug, ldebuglt
      parameter (ldebug=.false.)
      parameter (ldebuglt=.false.)
      parameter (zero = (0d0,0d0) )
      integer i,if1,k,mu

      uuzz = 0d0
      uuzz = 0d0
      ddzz = 0d0
      ddzz = 0d0

      nmin = nmax+1
      do i = 1,2
         diagToPhys(physToDiag(i)) = i
      enddo

      polcol = 1d0/(4d0*N**2*xi(1)*xi(2))
      polcolg = 1d0/(4d0*(N**2-1)**2*xi(1)*xi(2))

      do mu = 0,3
         do i = 1,2
            pbar(mu,physToDiag(i))=p(mu,i)
         enddo
      enddo

      do mu = 0,3             ! kinematics for ZZ
         pbar(mu,3) = v(mu,1) 
         pbar(mu,4) = v(mu,2) 
         pbar(mu,5) = v(mu,3) 
         pbar(mu,6) = v(mu,4) 
      enddo

      if ((nlo .eq. -8).or.((nlo.eq.1).and.(physToDiag(1).eq.1))) then
c gluon-fusion loop contributions
        call GGZZ(pbar,fsign,nlo,1,uuzz)
        iflav(1)=0
        iflav(2)=0
        k=FL_VV(iflav,1)
  
        res(k)=pdf(sign1*iflav(1),1)
     &        *pdf(sign2*iflav(2),2)
     &        *uuzz*polcolg

      endif
      if (nlo .ne. -8) then
#ifdef WITH_MADGRAPH
        if (ldebug) then
           call suu_zz(pbar,fsign,uuzz2)
c           call suu_zz_had_u(pbar,fsign,uuzz2)  ! for hadronic decay into uu~
           call sdd_zz(pbar,fsign,ddzz2)
c           call sdd_zz_had_u(pbar,fsign,ddzz2)  ! for hadronic decay into uu~
        endif
#endif

        call QQZZ(pbar,fsign,nlo,1,uuzz,ddzz)

#ifdef WITH_MADGRAPH
        if (ldebug) then
           print*,"uuzz =   ",uuzz
           print*,"uuzz2 =  ",uuzz2
           print*,"ddzz =   ",ddzz
           print*,"ddzz2 =  ",ddzz2
           read(*,*)
        endif
#endif

        NCmatrixelt(0,0)=uuzz
        NCmatrixelt(1,1)=ddzz

        do if1=1,nfl
              iflav(1)=if1*fsign(physToDiag(1))
              iflav(2)=if1*fsign(physToDiag(2))
              k=FL_VV(iflav,1)
  
              res(k)=pdf(sign1*iflav(1),1)
     &              *pdf(sign2*iflav(2),2)
     &              *NCmatrixelt(mod(if1,2),mod(if1,2))*polcol
  
        enddo 

      endif

      nmax = FL_VV(iflav,-2)

      end
c******************************************************************
c
c   end subroutine qq_ZZ
c
c*****************************************************************

