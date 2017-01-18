c*****************************************************************************
c
c    begin function m2s_qqVV
c
c*****************************************************************************
c
c       Adapted from m2s_wbfzh by D. Zeppenfeld
c
c
      real*8 function m2s_qqVV(
     &      bos,  !in:  XX Boson identifier,X=1 photon, X=2 Z, X=3,4 W+-, X=6 H
     &      nlo,  !in:  Born only for nlo=0, +interference for 1
     &      xi,   !in:  Feynman x parameters of incoming fermions
     &      p,    !in:  momenta associated with external fermions
     &      v,    !in:  momenta associated with gauge boson dec. fermions
     &      rn    !in:  random number array
     &                        )
      implicit none
c
c declare input/output variables
c
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/process.inc"
#include "tensor.inc"

      integer bos, nlo
      real*8 xi(nx), p(0:3,max_p), v(0:3,n_v,3), rn(1), m2s

c  helicity selection
      INTEGER h
      COMMON /hcount / h

      DOUBLE PRECISION zero
      PARAMETER (zero=0d0)

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
      LOGICAL ldebug, ldebuglt
      parameter (ldebug=.false.)		!output debug information
      PARAMETER (ldebuglt=.false.)

c variables for leshouches interface
      DOUBLE PRECISION weight,rnumb,RandomNumber

      integer N_gen_W
      real*8 fac_W, mjj2
      external mjj2

      if ( init .eq. 0 ) then

         if (bos.eq.34 .or. bos.eq.43) then
            write(6,*) " "
            write(6,*) "W+W- amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
            call DIBOStoww_reset
         endif

         call printnfl(.true.)
         init = 1
         do i = 1,maxnumsubproc
            res(i) = 0
         enddo
      else
         do i = 1,nmax
            res(i) = 0
         enddo
      endif

c...Les Houches interface
        if ((lha.or.hepmc).and..not.doNLO) then
           helicity(1)=-1
           helicity(2)= 1
           helicity(3)=-1
           helicity(4)= 1
        endif

c for WW->4l precalculate A->WW,Z->WW,AZ->WW etc leptonic tensors

      CALL IXXXXX(v(0,2,1),ZERO ,1,-1,wep(1,1)) !W(1,4))          !e+       
      CALL OXXXXX(v(0,1,1),ZERO ,-1,1,wve(1,1)) !W(1,3))          !ve 
      CALL OXXXXX(v(0,3,1),ZERO ,-1,1,wmu(1,1)) !W(1,5))          !mu-      
      CALL IXXXXX(v(0,4,1),ZERO ,1,-1,wvm(1,1)) !W(1,6))          !vm~
      CALL JIOXXX(wep(1,1),wve(1,1),GWF ,WMASS,WWIDTH,wp(1,1)) !W(1,7))!W+
      CALL JIOXXX(wvm(1,1),wmu(1,1),GWF ,WMASS,WWIDTH,wm(1,1)) !W(1,8))!W-

      do mu = 0,3
         qp(mu,1) = v(mu,1,1)+v(mu,2,1)   !W+
         qm(mu,1) = v(mu,3,1)+v(mu,4,1)   !W-
         qww(mu,1) = qp(mu,1)+qm(mu,1)    !Z_all
      enddo
      qp(4,1) = dotrr(qp(0,1),qp(0,1))
      qm(4,1) = dotrr(qm(0,1),qm(0,1))
      qww(4,1) = dotrr(qww(0,1),qww(0,1))

c helicity always 1
      h=1

c leptonic tensors

      SELECT CASE(process)
      CASE(WPhadWMlep,GFWPhadWMlep,WPlepWMhad,GFWPlepWMhad)
        ! for hadronic decays anom and SM use the same tensors
        call anomal_formfactor(qww(0,1),qp(0,1),qm(0,1),qm(0,1))       ! last parameter is "dummy variable"
        call DIBOSatoww_had(v,1,aww(0,1))
        call DIBOSztoww_had(v,1,zww(0,1))
        if ((nlo .eq. -8) .or. (nlo.eq.1)) 
     &    call ww_htoww_had(v(0,1,1),1,w4l(h,1),0) !H -> e+ ve vu~ mu-
      CASE DEFAULT
        if (with_anom) then
           call anomal_formfactor(qww(0,1),qp(0,1),qm(0,1),qm(0,1))       ! last parameter is "dummy variable"
           call DIBOSatoww_anom(v,1,aww(0,1))
           call DIBOSztoww_anom(v,1,zww(0,1))
        else  ! standard model:
           call DIBOSatoww(v,1,aww(0,1))
           call DIBOSztoww(v,1,zww(0,1))
        end if
        if ((nlo .eq. -8) .or. (nlo.eq.1)) 
     &    call ww_htoww(v(0,1,1),1,w4l(h,1),0) !H -> e+ ve vu~ mu-
      END SELECT


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


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     Now sum over the subprocesses contributing to W W production

      nmax = FL_VV(fsign,-1)            !reset counter for subprocesses to 0

C*******************  q1 qbar2 ---> W W   **************************

c   physToDiag(ext.momentum label) = Feynman diagram label

      physToDiag(1)=1     !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=2

      fsign(1) = 1
      fsign(2) = -1
      fsign(3) = 1
      fsign(4) = -1
      fsign(5) = 1
      fsign(6) = -1

      call qq_VV(xi,p,v,physToDiag,fsign,bos,nlo,
     1            pdf,res,nmin,nmax)

C*******************  q1bar q2 ---> W W   **********************

      if (nlo .ne. -8) then
cc gg is symmetric -> would be double-counting
      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=1

      fsign(1) = 1
      fsign(2) = -1

      call qq_VV(xi,p,v,physToDiag,fsign,bos, nlo,
     1            pdf,res,nmin,nmax)
      endif

c*****************  end of process evaluation  *********************

      if (init.eq.1) then
         init = init+1
      endif

      SELECT CASE(process)
      CASE(WPhadWMlep,GFWPhadWMlep)
        call hadDecayFactor_W(finalquarks(1), abs(mjj2(v(0,1,1),v(0,2,1))), N_gen_W, fac_W)
        DO IRES = 1,NMAX
          res(IRES) = res(IRES) * fac_W            ! factor fac_W for hadronic decay
        ENDDO
      CASE(WPlepWMhad,GFWPlepWMhad)
        call hadDecayFactor_W(finalquarks(1), abs(mjj2(v(0,3,1),v(0,4,1))), N_gen_W, fac_W)
        DO IRES = 1,NMAX
          res(IRES) = res(IRES) * fac_W            ! factor fac_W for hadronic decay
        ENDDO
      END SELECT

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
         if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
            rnumb=RandomNumber()
            SELECT CASE(process)
            CASE(WPhadWMlep,GFWPhadWMlep)
               finalquarks_psp(1) =  2 + 2* INT(rnumb*2)
               finalquarks_psp(2) = -1 - 2* INT(rnumb*2)
            CASE(WPlepWMhad,GFWPlepWMhad)
               finalquarks_psp(1) =  1 + 2* INT(rnumb*2)
               finalquarks_psp(2) = -2 - 2* INT(rnumb*2)
            END SELECT
         endif
      endif

      m2s = m2s     ! factor for electrons and muons and tau

      m2s_qqVV = m2s
      
      RETURN
      END

c*****************************************************************************
c
c    end function m2s_qqVV
c
c*****************************************************************************

c******************************************************************
c
c   begin subroutine qq_VV
c
c*****************************************************************

      subroutine qq_VV(xi,p,v,physToDiag,fsign,bos, nlo,
     1                  pdf,res,nmin,nmax)
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"

      real*8 p(0:3,n_p), v(0:3,n_v), xi(nx)
      real*8 pdf(-6:6,2), res(maxnumsubproc)
      integer physToDiag(2), fsign(2+n_v), bos, nlo, nmin, nmax
c
c wbf_zh calls the amplitude square routines 
c             qqVV    for qq --> VV --> 4 leptons
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
      real*8 uuww, ddww, uuww2, ddww2
      real*8 pbar(0:3,2+n_v), dotrr, polcol, polcolg
      double precision NCmatrixelt(0:1,0:1)
      double complex zero
      integer iflav(4), diagToPhys(2), FL_VV
      external dotrr, FL_VV


      logical ldebug, ldebuglt
      parameter (ldebug=.false.)
      parameter (ldebuglt=.false.)
      parameter (zero = (0d0,0d0) )
      integer i,if1,k,mu

      uuww = 0d0
      uuww2 = 0d0
      ddww = 0d0
      ddww2 = 0d0

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


         do mu = 0,3             ! kinematics for WW --> 4 lepton decay
            pbar(mu,3) = v(mu,1) ! ve
            pbar(mu,4) = v(mu,2) ! e+
            pbar(mu,5) = v(mu,3) ! mu-
            pbar(mu,6) = v(mu,4) ! nu_mu~
         enddo

      if ((nlo .eq. -8).or.((nlo.eq.1).and.(physToDiag(1).eq.1))) then
c gluon-fusion loop contributions
        call GGWW(pbar,fsign,nlo,1,uuww)
        iflav(1)=0
        iflav(2)=0
        k=FL_VV(iflav,1)
  
        res(k)=pdf(sign1*iflav(1),1)
     &        *pdf(sign2*iflav(2),2)
     &        *uuww*polcolg

      endif
      if (nlo .ne. -8) then
#ifdef WITH_MADGRAPH
         if(ldebug) then
               call suu_ww(pbar,fsign,uuww2)
               ! call suu_ww_had(pbar,fsign,uuww2)  ! for hadronic decay
               call sdd_ww(pbar,fsign,ddww2)
               ! call sdd_ww_had(pbar,fsign,ddww2)  ! for hadronic decay
	    endif
#endif
       
	   call QQVV(pbar,fsign,nlo,1,uuww,ddww)

#ifdef WITH_MADGRAPH
         if (ldebug) then
            if(abs(uuww/uuww2 - 1d0).ge.1.0d-5) then
               print*,"uuww =   ",uuww
               print*,"uuww2 =  ",uuww2
            endif
            if(abs(ddww/ddww2 - 1d0).ge.1.0d-5) then
               print*,"ddww =   ",ddww
               print*,"ddww2 =  ",ddww2
            endif
         endif
#endif

        NCmatrixelt(0,0)=uuww
        NCmatrixelt(1,1)=ddww


        do if1=1,(nfl/2)*2
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
c   end subroutine qq_VV
c
c*****************************************************************


c******************************************************************
c
c   begin function FL_VV
c
c*****************************************************************
      INTEGER FUNCTION FL_VV(iflav,colstruc)

      IMPLICIT NONE

      INTEGER iflav(4),colstruc    ! input for color and flavor assignment
c  fill hepup color and flavor assignments for WBF processes with
c  two quark lines only

c  id1,2 are flavor identifiers for incoming quarks
c  id3,4 are flavor identifiers for outgoing quarks
c  colstruc = is irrelevant here
c
c  color structure of these interactions is trivial
c  and colstruc was included only to make this program 
c  more consistent with the others.

c  Note that colstruc is also a flag:
c  if colstruc=-1, then we are resetting.

c  The color structure is always the same:
c  color always flows from id1->id3 and id2->id4
c  no gluons are involved.

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
          call fillColoredPartons_VV(iflav,
     1                               listposition+numdecay)
        endif       ! numresets.eq.0
      endif         
      FL_VV=listposition
      end
c******************************************************************************
c
c   end function FL_VV
c
c******************************************************************************

c*****************************************************************************
c
c    begin  subroutine fillColoredPartons_VV
c
c*****************************************************************************

      subroutine fillColoredPartons_VV(iflav,listposition)

c  assigns values to the variables in the common block localHEPUP
c  in particular, this subroutine assigns values to those variables that 
c  will be stored in the lookup tables generated by writeHEPUPtable.  
c  As the name suggests, this routine only stores the information for the 
c  colored partons.  Particles without color will be dealt with in the 
c  subroutine fillColorless.

      IMPLICIT NONE

c  include files
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/process.inc"

c  input parameters
      INTEGER iflav(2),id1,id2,listposition !2 quark flavors and subprocess number

c  local variables

      
      id1=iflav(1)
      id2=iflav(2)
   
      select case(process)
      case(WW,WPZ,WMZ,ZZ,GFWW,GFZZ,
     &     WPhadWMlep,WPlepWMhad,GFWPhadWMlep,GFWPlepWMhad,
     &     ZZhad,GFZZhad,
     &     WPhadZ,WPZhad,WMhadZ,WMZhad)
         lnup(listposition)=numParticles
         ! there is no point in using numParticles as this is just another constant = 8
      case(WP_only,WM_only)
         lnup(listposition)=5
      case(WPH,WMH) 
         lnup(listposition)=6
      case(WPH_AA,WPH_mu,WPH_tau,WPH_bbar,WMH_AA,WMH_mu,WMH_tau,WMH_bbar)
         lnup(listposition)=8
      case(WPA,WMA,ZA,GFZA)
         lnup(listposition)=numParticles-2
      case(AA,GFAA)
         lnup(listposition)=numParticles-4
      case default
         print*,'wrong process number in fillColoredPartons_VV: ',process
         stop
      end select

      listup(1,listposition)=-1   !incoming quarks
      listup(2,listposition)=-1
 
      if (id1.ne.0) then
        lidup(1,listposition)=id1     ! the 2 quarks
      else
        lidup(1,listposition)=21      ! id=0 is gluon, lhe expects ID 21
      endif
      if (id2.ne.0) then
        lidup(2,listposition)=id2
      else
        lidup(2,listposition)=21      ! id=0 is gluon, lhe expects ID 21
      endif

      if (id1.gt.0) then
        licolup(1,1,listposition)=501
        licolup(2,1,listposition)=0
        licolup(1,2,listposition)=0
        licolup(2,2,listposition)=501
      else if (id1.lt.0) then
        licolup(1,1,listposition)=0
        licolup(2,1,listposition)=501
        licolup(1,2,listposition)=501
        licolup(2,2,listposition)=0
      else 
        licolup(1,1,listposition)=501
        licolup(2,1,listposition)=502
        licolup(1,2,listposition)=502
        licolup(2,2,listposition)=501
      endif

      end  ! fillColoredPartons_VV
c*****************************************************************************
c
c    end subroutine fillColoredPartons_VV
c
c*****************************************************************************


      subroutine DIBOStoww_reset
      implicit none
#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
      complex*16 zero
      parameter (zero=(0d0,0d0))
      integer j,mu

      do mu = 0,3
         do j = 1,3
            aww(mu,j) = zero
            zww(mu,j) = zero
         enddo
      enddo
      print*," v(v) to WW tensors initialized to ",zero 
      return
      end
