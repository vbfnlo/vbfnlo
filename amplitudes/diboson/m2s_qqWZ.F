c*********************  subroutine m2s_qqWZ *********************************        
c
c       Adapted from m2s_VV by Vera Hankele
c
c       Johannes Bellm, <johannes@particle.physik.uni-karlsruhe.de>
c	Initial version: Juli 2011
c	Last modified:
c
c***************************************************************************
      real*8 function m2s_qqWZ(
     &      bos,  !in:  XX Boson identifier,X=1 photon, X=2 Z, X=3,4 W+-, X=6 H
     &      nlo,  !in:  Born only for nlo=0, +interference for 1
     &      xi,   !in:  Feynman x parameters of incoming fermions
     &      p,    !in:  momenta associated with external fermions
     &      v,    !in:  momenta associated with gauge boson dec. fermions
     &      rn    !in:  random number array
     &                        )

      use globalvars, only: ldoblha

      implicit none
c
c declare input/output variables
c
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/process.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"
#include "tensor.inc"

      integer bos, nlo
      real*8 xi(nx), p(0:3,max_p), v(0:3,n_v,3), rn(1), m2s

c  helicity selection
c      integer jsig, jsig1, jsig3, jsig5
c      common /chelsum/ jsig,jsig1,jsig3,jsig5


c  helicity selection
      INTEGER h, it 
      COMMON /hcount / h

      double precision zero
      parameter (zero=0d0)

c declare external functions

      real*8 dotrr
      integer FL_VV
      external dotrr
      external FL_VV

c declare local variables

      real*8 q_sf
      integer init/0/, I, mu
      save init

c declare variables for summation process

      INTEGER ires

      integer FSIGN(2+n_v), physToDiag(2), nmin, nmax
      save nmax
c
c store contributions from subprocess ip in res(ip) where
      real*8 res(maxnumsubproc)
c
c declare parton distribution variables
c
      real*8 pdf(-6:6,2)
c
c define program switches
c
      logical ldebug
      parameter (ldebug=.false.)		!output debug information

c variables for leshouches interface
      DOUBLE PRECISION weight,rnumb,RandomNumber

c variables for hadronic decays
      integer N_gen_W
      real*8 fac_W, mjj2
      external mjj2
      integer N_gen_up, N_gen_down
      real*8 fac_Z_up, fac_Z_down


      if ( init .eq. 0  .and. (.not. ldoblha) ) then
         if(bos.eq.32) then
            write(6,*) " "
            write(6,*) "W+Z amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
         else if(bos.eq.42) then
            write(6,*) " "
            write(6,*) "W-Z amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
         else
            write(6,*) " wrong call of m2s_qqWZ, bos = ",bos
            write(6,*) "-----------------------------------------------"
         endif
         call DIBOSwtowz_reset

         call printnfl(.true.)
!          h = 1
         init = 1
         do i = 1,maxnumsubproc
            res(i) = 0
         enddo
      else if ( (init .eq. 0) .and. ldoblha ) then
         init = 1
      else
         do i = 1,nmax
            res(i) = 0
         enddo
      endif


c...Les Houches interface

c lepton helicities for Z -> ta- ta+ not fixed
c -> sum over all possible helicities in |M|**2
 
c select helicity: h ... random number for lepton helicity combination (h=1:2)
      if (ldoblha) then
         h = blha_ranhelcomb
      endif

      h = mod(h,2) + 1
      it = (-1)**(h+1)

      if ((lha.or.hepmc) .and. (.not. doNLO)) then
         helicity(1)=-1
         helicity(2)= 1
         helicity(3)=-it
         helicity(4)= it
      endif

c for WW->4l precalculate A->WW,Z->WW,AZ->WW etc leptonic tensors

      CALL OXXXXX(v(0,1,1),ZERO ,-1,1,wve(1,1)) !W(1,3))     !ve
      CALL IXXXXX(v(0,2,1),ZERO ,1,-1,wep(1,1)) !W(1,4))          !e+          
      CALL OXXXXX(v(0,3,1),ZERO ,-it,1,wmu(1,1)) !W(1,5))          !mu-      
      CALL IXXXXX(v(0,4,1),ZERO ,it,-1,wmubar(1,1)) !W(1,6))          !mu+
      

      do mu = 0,3
         qp(mu,1) = v(mu,1,1)+v(mu,2,1)   !W+
         qz(mu,1) = v(mu,3,1)+v(mu,4,1)   !Z
         qa(mu,1) = v(mu,3,1)+v(mu,4,1)   !A
      enddo
      qp(4,1) = dotrr(qp(0,1),qp(0,1))
      qz(4,1) = dotrr(qz(0,1),qz(0,1))
      qa(4,1) = dotrr(qa(0,1),qa(0,1))

c leptonic tensors
      SELECT CASE(process)
      CASE(WPhadZ,WPZhad,WMhadZ,WMZhad)
        ! for general output to all combinations: up-type first
        if ((finalquarks(1).eq.93 .or. finalquarks(1).eq.94) .and.
     &      (process.eq.WPZhad .or. process.eq.WMZhad)) then
          ! set couplings for Z hadronic decay into uubar
          call setZtouu
        endif
        ! for hadronic decays anom and SM use the same tensors
        if (bos.eq.32) then
           call DIBOSwptowpz_had(v,1,wtowz(0,1))
        else if (bos.eq.42) then
           call DIBOSwmtowmz_had(v,1,wtowz(0,1))
        endif
        ! later you need the an: A-> mumu,wp: w+->e+ve and zn: Z->mumu 
        CALL JIOXXX(wep(1,1),wve(1,1),GWF ,WMASS,WWIDTH,wp(1,1)) !!W+
        CALL JIOXXX(wmubar(1,1),wmu(1,1),GZ_ZF ,ZMASS,ZWIDTH,zn(1,1)) !!Z
        CALL JIOXXX(wmubar(1,1),wmu(1,1),GZ_AF ,ZERO,ZERO,an(1,1)) !!A
      CASE DEFAULT
        if (with_anom) then ! anomalous gauge boson couplings noanomal
            if (bos.eq.32) then
                   call DIBOSwptowpz_an(v,1,wtowz(0,1))
            else if (bos.eq.42) then
                   call DIBOSwmtowmz_an(v,1,wtowz(0,1))
            endif
        else  ! SM
            if (bos.eq.32) then
                  call DIBOSwptowpz(v,1,wtowz(0,1))
            else if (bos.eq.42) then
                  call DIBOSwmtowmz(v,1,wtowz(0,1))
            endif
        endif
        ! later you need the an: A-> mumu,wp: w+->e+ve and zn: Z->mumu 
        CALL JIOXXX(wep(1,1),wve(1,1),GWF ,WMASS,WWIDTH,wp(1,1)) !!W+
        CALL JIOXXX(wmubar(1,1),wmu(1,1),GZL ,ZMASS,ZWIDTH,zn(1,1)) !!Z
        CALL JIOXXX(wmubar(1,1),wmu(1,1),GAL ,ZERO,ZERO,an(1,1)) !!A
      END SELECT


      if (ldoblha) then
         physToDiag(1:2)=blha_physdiagmap(1:2,blha_cursubproc,blha_curproc)
         fsign(1:2) = blha_fsign(1:2,blha_cursubproc,blha_curproc)
         ! fsign for leptions, TODO: check
         fsign(3) = 1
         fsign(4) = -1
         fsign(5) = 1
         fsign(6) = -1 

         call qq_WZ(xi,p,v,physToDiag,fsign,bos,nlo,pdf,res,nmin,nmax) 
         return
      endif

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
         do i = 1, n_p-1
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

        
      DO IRES = 1,NMAX
         RES(IRES)=0d0
      ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     Now sum over the subprocesses contributing to Z production

      nmax = FL_VV(fsign,-1)            !reset counter for subprocesses to 0

C*******************  q1 q2 ---> V V   **************************

c   physToDiag(ext.momentum label) = Feynman diagram label

      physToDiag(1)=1     !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=2

      fsign(1) = 1 
      fsign(2) = -1 
      fsign(3) = 1
      fsign(4) = -1
      fsign(5) = 1
      fsign(6) = -1 

      call qq_WZ(xi,p,v,physToDiag,fsign,bos,nlo,
     1            pdf,res,nmin,nmax) 


C*******************  q1<->q2 ---> V V   **********************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=1

      fsign(1) = 1
      fsign(2) = -1

      call qq_WZ(xi,p,v,physToDiag,fsign,bos, nlo,
     1            pdf,res,nmin,nmax)


c**************  end of process evaluation part one ******************


! now run again for down-type decay products if 93 93 / 94 94 is 
! requested and the z decays hadronically.

      SELECT CASE(procid)
      CASE(WPZhad, WMZhad)
       if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
        ! set couplings for Z hadronic decay into ddbar
        call setZtodd
        if (bos.eq.32) then
           call DIBOSwptowpz_had(v,1,wtowz(0,1))
        else if (bos.eq.42) then
           call DIBOSwmtowmz_had(v,1,wtowz(0,1))
        endif
        CALL JIOXXX(wmubar(1,1),wmu(1,1),GZ_ZF ,ZMASS,ZWIDTH,zn(1,1)) !!Z
        CALL JIOXXX(wmubar(1,1),wmu(1,1),GZ_AF ,ZERO,ZERO,an(1,1)) !!A


C*******************  q1 q2 ---> V V   **************************

c   physToDiag(ext.momentum label) = Feynman diagram label

      physToDiag(1)=1     !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=2

      fsign(1) = 1 
      fsign(2) = -1 
      fsign(3) = 1
      fsign(4) = -1
      fsign(5) = 1
      fsign(6) = -1 

      call qq_WZ(xi,p,v,physToDiag,fsign,bos,nlo,
     1            pdf,res,nmin,nmax) 


C*******************  q1<->q2 ---> V V   **********************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=1

      fsign(1) = 1
      fsign(2) = -1

      call qq_WZ(xi,p,v,physToDiag,fsign,bos, nlo,
     1            pdf,res,nmin,nmax)


c***  end of extra process evaluation for hadronic decay  ***

       endif
      end select


c*****************  end of process evaluation  *********************


      SELECT CASE(procid)
      CASE(WPhadZ,WMhadZ)
        call hadDecayFactor_W(finalquarks(1), abs(mjj2(v(0,1,1),v(0,2,1))), N_gen_W, fac_W)
        DO IRES = 1,NMAX
          res(IRES) = res(IRES) * fac_W            ! factor fac_W for hadronic decay
        ENDDO
      CASE(WPZhad,WMZhad)
        call hadDecayFactor_Z(finalquarks(1), abs(mjj2(v(0,3,1),v(0,4,1))), N_gen_up, N_gen_down, fac_Z_up, fac_Z_down)
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

      m2s = m2s*2d0     ! factor 2 for random helicity summation
!         if(m2s.le.0d0)then
!           print*,"m2s",m2s
!         endif
      m2s_qqWZ = m2s
      RETURN
      END

c*****************************************************************************
c
c    end function m2s_qqWZ
c
c*****************************************************************************

c******************************************************************
c
c   begin subroutine qq_WZ
c
c*****************************************************************

      subroutine qq_WZ(xi,p,v,physToDiag,fsign,bos, nlo,
     1                  pdf,res,nmin,nmax)

      use globalvars, only: ldoblha

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

      real*8 p(0:3,max_p-1), v(0:3,n_v), xi(nx)
      real*8 pdf(-6:6,2), res(maxnumsubproc)
      integer physToDiag(2), fsign(2+n_v), bos, nlo, nmin, nmax
c
c wbf_zh calls the amplitude square routines 
c             qqWZ    for qq --> WZ --> 4 leptons
c
c  INPUT:  p(0:3,4)        external physical parton momenta
c          v(0:3,nv)       Z decay momenta, for Higgs production only
c                          the sum q(mu) = v(mu,1)+...+v(mu,nv) is needed
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
      real*8 udwz, udwz2
      real*8 pbar(0:3,2+n_v), dotrr, polcol
      double precision NCmatrixelt(0:1,0:1)
      double complex zero
      integer iflav(4), diagToPhys(2), FL_VV
      external dotrr, FL_VV


      logical ldebug
      parameter (ldebug=.false.)
      parameter (zero = (0d0,0d0) )
      integer i,if1,k,mu

      udwz = 0d0
      udwz2 = 0d0


      nmin = nmax+1
      do i = 1,2
         diagToPhys(physToDiag(i)) = i
      enddo

      polcol = 1d0/(4d0*N**2*xi(1)*xi(2))

      do i = 1,2
         pbar(0:3,physToDiag(i))=p(0:3,i)
      enddo

      ! kinematics for WZ --> 4 lepton decay
      pbar(0:3,3) = v(0:3,1) ! ve
      pbar(0:3,4) = v(0:3,2) ! e+
      pbar(0:3,5) = v(0:3,3) ! mu-
      pbar(0:3,6) = v(0:3,4) ! mu+

#ifdef WITH_MADGRAPH
         if(ldebug) then

            if (bos.eq.32) then
               ! call sudwpz(pbar,fsign,udwz2)
               call sudwphadz(pbar,fsign,udwz2)         ! hadronic decay. for Z only decay into dd~!
            else if (bos.eq.42) then
               call sduwmz(pbar,fsign,udwz2)              ! here the hadronic stuff is also in the normal routine.
            endif
         endif
#endif

         call QQWZ(pbar,fsign,nlo,1,udwz,bos)

#ifdef WITH_MADGRAPH
         if (ldebug) then
            if(abs(udwz/udwz2 -1).ge.1.0d-5) then
               print*,"udwz =   ",udwz
               print*,"udwz2 =  ",udwz2
               print*, "turn off anomal couplings! or NLO!!" 
            endif
         endif
#endif

      if (ldoblha) then
        return
      endif

      NCmatrixelt(0,0)=udwz


      do if1=2,(nfl/2)*2,2
         if (bos.eq.32) then
            iflav(1)=( if1 - physToDiag(1) + 1 )*fsign(physToDiag(1))
            iflav(2)=( if1 - physToDiag(2) + 1 )*fsign(physToDiag(2))
         else if (bos.eq.42) then
            iflav(1)=( if1 + physToDiag(1) - 2 )*fsign(physToDiag(1))
            iflav(2)=( if1 + physToDiag(2) - 2 )*fsign(physToDiag(2))
         endif
         k=FL_VV(iflav,1)

         res(k)=pdf(sign1*iflav(1),1)
     &         *pdf(sign2*iflav(2),2)
     &         *NCmatrixelt(0,0)*polcol

!         if(res(k).le.0d0)then
!           print*,"res(",k,")",res(k)
!           print*,"pdf(",sign1*iflav(1),",1)",pdf(sign1*iflav(1),1)
!           print*,"pdf(",sign2*iflav(2),",2)",pdf(sign2*iflav(2),2)
!           print*,"NCmatrixelt(0,0)",NCmatrixelt(0,0)
!           print*,"polcol",polcol
! 
!         endif


      enddo 

      nmax = FL_VV(iflav,-2)

      end
c******************************************************************
c
c   end subroutine qq_WZ
c
c*****************************************************************


      subroutine DIBOSwtowz_reset
      implicit none
#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
      complex*16 zero
      parameter (zero=(0d0,0d0))
      integer j,mu

      do mu = 0,3
         do j = 1,3
            wtowz(mu,j) = zero            
         enddo
      enddo
      print*," W to WZ  tensors initialized to ",zero 
      return
      end

      
