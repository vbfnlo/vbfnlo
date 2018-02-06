c*****************************************************************************
c
c    begin function m2s_qqVVj2
c
c*****************************************************************************
c
c       Adapted from m2s_qqVV by V. Hankele
c
c
      real*8 function m2s_qqVVj2(
     &      bos,  !in:  XX Boson identifier,X=1 photon, X=2 Z, X=3,4 W+-, X=6 H
     &      nlo,  !in:  Born only for nlo=0, +interference for 1
     &      xi,   !in:  Feynman x parameters of incoming fermions
     &      p,    !in:  momenta associated with external fermions
     &      v,    !in:  momenta associated with gauge boson dec. fermions
     &      rn,   !in:  random number array
     &      ps_number !in:  phase-space number
     &                        )
      implicit none
c
c declare input/output variables
c
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/kkcoupl.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/process.inc"
#include "tensor.inc"

      

      integer bos, nlo, ps_number
      real*8 xi(nx), p(0:3,max_p), v(0:3,n_v,11), rn(1), m2s

c  helicity selection
c      integer jsig, jsig1, jsig3, jsig5
c      common /chelsum/ jsig,jsig1,jsig3,jsig5

      DOUBLE PRECISION zero
      PARAMETER (zero=0d0)

c declare external functions

      real*8 dotrr
      integer FL_VVg
      external dotrr
      external FL_VVg

c declare local variables

      real*8 q_sf
      integer init/0/, I, mu
      save init

c declare variables for summation process

      INTEGER ires

      integer FSIGN(n_p+n_v), gsign, physToDiag(3), nmin, nmax
      save nmax

c declare variables for hadronic decay 
      integer N_gen_W
      real*8 fac_W, mjj2
      external mjj2
c
c store contributions from subprocess ip in res(ip,ID) where
c  ID = 1  : the real emission |M|^2 * pdf
c     = 2,3: subtraction terms for emision off upper or lower line
c     = 0  : subtracted result which drives integration, i.e
c res(*,0) = res(*,1)+res(*,2)+res(*,3)

      real*8 res(maxnumsubproc)
c
c declare parton distribution variables
c
      real*8 pdf(-6:6,2)
c
c define program switches
c
      LOGICAL ldebug
      parameter (ldebug=.false.)          !output debug information

c variables for leshouches interface
      DOUBLE PRECISION weight,rnumb,RandomNumber
c  helicity selection
      INTEGER h
      COMMON /hcount / h

      if ( init .eq. 0 ) then

         if (bos.eq.34 .or. bos.eq.43) then
            write(6,*) " "
            write(6,*) "W+W-j amplitude square information:"
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

      CALL IXXXXX(v(0,2,1),ZERO ,1,-1,wep(1,1)) !W(1,4))       !e+  | e+ | d~ 
      CALL OXXXXX(v(0,1,1),ZERO ,-1,1,wve(1,1)) !W(1,3))       !ve  | ve | u
      CALL OXXXXX(v(0,3,1),ZERO ,-1,1,wmu(1,1)) !W(1,5))       !mu- | d  | e-
      CALL IXXXXX(v(0,4,1),ZERO ,1,-1,wvm(1,1)) !W(1,6))       !vm~ | u~ | ve~
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

      SELECT CASE(process)
      CASE(WPHADWMJ, WPWMHADJ, GFWPHADWMJ, GFWPWMHADJ)
         !hadronic tensors
         call DIBOSZtoQQll(v, 1, zww(0,1))
         call DIBOSAtoQQll(v, 1, aww(0,1))
      CASE DEFAULT
         !leptonic tensors
        if (with_kk) then
#ifdef WITH_KK
           call DIBOSatoww_kk(v,1,aww(0,1))
           call DIBOSztoww_kk(v,1,zww(0,1))
           call DIBOSz1toww_kk(v,1,z1ww(0,1))
           call DIBOSz2toww_kk(v,1,z2ww(0,1))
#endif            
        else
           call DIBOSatoww(v,1,aww(0,1))
           call DIBOSztoww(v,1,zww(0,1))
        end if
      END SELECT

      if ((nlo .eq. -6) .or. (nlo .eq. -8) .or. (nlo.eq.1)) 
     &  call ww_htoww(v(0,1,1),1,w4l(h,1),0) !H -> e+ ve vu~ mu-

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

C     Now sum over the subprocesses contributing to WWj production

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

      call qq_VVj2(xi,p,v,physToDiag,fsign,gsign,bos, nlo,
     &              pdf,res,nmin,nmax,ps_number)

      if (nlo .ne. -8) then
cc gg is symmetric -> would be double-counting
C*******************  qbar1 q2 ---> V V g   **********************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=1
      physToDiag(3)=3    ! gluon

      fsign(1) = 1
      fsign(2) = -1

      call qq_VVj2(xi,p,v,physToDiag,fsign,gsign,bos, nlo,
     &              pdf,res,nmin,nmax,ps_number)

C*******************  g q1 ---> V V q2  ***********************************

      physToDiag(1)=3
      physToDiag(2)=1
      physToDiag(3)=2

      fsign(1) = 1
      fsign(2) = 1
      gsign    = -1

      call qq_VVj2(xi,p,v,physToDiag,fsign,gsign,bos, nlo,
     &              pdf,res,nmin,nmax,ps_number)

C*******************  qbar2 g ---> V V qbar1   **********************

      physToDiag(1)=2
      physToDiag(2)=3
      physToDiag(3)=1

      fsign(1) = -1
      fsign(2) = -1
      gsign = -1

      call qq_VVj2(xi,p,v,physToDiag,fsign,gsign,bos, nlo,
     &              pdf,res,nmin,nmax,ps_number)

C*******************  g q2bar ---> V V q1bar  **********************

      physToDiag(1)=3
      physToDiag(2)=2
      physToDiag(3)=1

      fsign(1) = -1
      fsign(2) = -1
      gsign = -1

      call qq_VVj2(xi,p,v,physToDiag,fsign,gsign,bos, nlo,
     &              pdf,res,nmin,nmax,ps_number)

C*******************  q1 g ---> V V q2   **********************

      physToDiag(1)=1
      physToDiag(2)=3
      physToDiag(3)=2

      fsign(1) = 1
      fsign(2) = 1

      call qq_VVj2(xi,p,v,physToDiag,fsign,gsign,bos, nlo,
     &              pdf,res,nmin,nmax,ps_number)

      endif
c*****************  end of process evaluation  *********************

      SELECT CASE(procid)
      CASE(WPhadWMj,GFWPhadWMj)
       call hadDecayFactor_W(finalquarks(1), abs(mjj2(v(0,1,1),v(0,2,1))), N_gen_W, fac_W)
        DO IRES = 1,NMAX
          res(IRES) = res(IRES) * fac_W            ! factor fac_W for hadronic decay
        ENDDO
      CASE(WPWMhadj,GFWPWMhadj)
       call hadDecayFactor_W(finalquarks(1), abs(mjj2(v(0,3,1),v(0,4,1))), N_gen_W, fac_W)
        DO IRES = 1,NMAX
          res(IRES) = res(IRES) * fac_W            ! factor fac_W for hadronic decay
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

      m2s_qqVVj2 = m2s
      
      RETURN
      END

c*****************************************************************************
c
c    end function m2s_qqVVj2
c
c*****************************************************************************

c******************************************************************
c
c   begin subroutine qq_VVj2
c
c*****************************************************************

      subroutine qq_VVj2(xi,p,v,physToDiag,fsign,gsign, bos, nlo,
     1                  pdf,res,nmin,nmax,ps_number)
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"

      real*8 p(0:3,max_p-1), v(0:3,n_v), xi(nx)
      real*8 pdf(-6:6,2), res(maxnumsubproc)
      integer physToDiag(3), fsign(n_p+n_v), gsign, 
     &        bos, nlo, nmin, nmax, ps_number
c
c qq_VVj2 calls the amplitude square routines 
c             qqVVj    for qq --> VV --> 4 leptons
c
c  INPUT:  p(0:3,3)        external physical parton momenta
c          v(0:3,nv)       decay momenta
c                          the sum q(mu) = v(mu,1)+...+v(mu,nv) is needed
c          physToDiag(3)   physToDiag(ext.mom. label) = Feynman diagram label
c          fsign           sign factors for amplitude calls; see qqZqq
c          nlo             nlo = 0: calculate LO only
c                          nlo = 1: calculate full NLO (Born + interference)
c  OUTPUT:
c          uuvv(ID)   output in format as in qqZqq(j), but corresponding
c          etc.            to res(ID) in m2s_qqZqq
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
      real*8 pbar(0:3,n_p+n_v), qbar(0:4), dotrr, polcolq, polcolg, polcolgg

      double precision NCmatrixelt(0:1,0:1)
      double complex zero
      integer iflav(4), diagToPhys(3), FL_VVg
      external dotrr, FL_VVg


      logical ldebug
      parameter (ldebug=.false.)
      parameter (zero = (0d0,0d0) )
      integer i,if1,j,k,mu

      uuww = 0d0
      uuww2 = 0d0
      ddww = 0d0
      ddww2 = 0d0

      nmin = nmax+1
      do i = 1,n_p
         diagToPhys(physToDiag(i)) = i
      enddo

      polcolq = 1d0/(4d0*N**2*xi(1)*xi(2))
      polcolg = 1d0/(4d0*N*(N**2-1)*xi(1)*xi(2))
      polcolgg= 1d0/(4d0*(N**2-1)**2*xi(1)*xi(2))

      do mu = 0,3
         do i = 1,n_p
            pbar(mu,physToDiag(i))=p(mu,i)
         enddo
         qbar(mu) = pbar(mu,3)
      enddo


      do mu = 0,3             ! kinematics for WW --> 4 lepton decay
         pbar(mu,3) = v(mu,1) ! ve
         pbar(mu,4) = v(mu,2) ! e+
         pbar(mu,5) = v(mu,3) ! mu-
         pbar(mu,6) = v(mu,4) ! nu_mu~
         pbar(mu,7) = qbar(mu)
      enddo
      fsign(7) = gsign

      if ((nlo .eq. -8).or.(
     &  (nlo.eq.1).and.(
     &    (physToDiag(1).eq.1) .and. (physToDiag(2).eq.2)
     & ))) then
c gluon-fusion loop contributions

c instead of anti-fermion here initial/final
        fsign(2)=-fsign(2)
        call GGWWg(pbar,fsign,qbar,gsign,nlo,1,uuww)
c switch back for nlo.eq.1
        fsign(2)=-fsign(2)
        iflav(1)=0
        iflav(2)=0
        k=FL_VVg(iflav,1)
  
        res(k)=pdf(sign1*iflav(1),1)
     &        *pdf(sign2*iflav(2),2)
     &        *uuww*polcolgg

      endif
      if (nlo .ne. -8) then

#ifdef WITH_MADGRAPH
        if(ldebug) then
          call suu_wwg(pbar,fsign,uuww2)
          call sdd_wwg(pbar,fsign,ddww2)
        endif
#endif
      

        call QQVVj(pbar,fsign,qbar,gsign,nlo,1,uuww,ddww,ps_number)

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


        if (gsign .eq. 1) then
           iflav(3) = 0           ! final state gluon id
           do if1=1,(nfl/2)*2
              iflav(1)=if1*fsign(physToDiag(1))
              iflav(2)=if1*fsign(physToDiag(2))
                 k=FL_VVg(iflav,1)
                 res(k)=pdf(sign1*iflav(1),1)
     &                 *pdf(sign2*iflav(2),2)
     &                 *NCmatrixelt(mod(if1,2),mod(if1,2))*polcolq
           enddo
        elseif (gsign.eq.-1) then !initial gluon section
           j = 0
           if (physToDiag(1).eq.3) then
              j=2              ! j=2: emission off 1,2 line
           elseif (physToDiag(2).eq.3) then
              j=3              ! j=3: emission off 3,4 line
           endif
           do if1=1,(nfl/2)*2
              iflav(j-1) = 0
              iflav(4-j)=if1*fsign(physToDiag(4-j))
              iflav(3)  =if1*fsign(physToDiag(3))
              k=FL_VVg(iflav,j)
              res(k)=pdf(sign1*iflav(1),1)
     &              *pdf(sign2*iflav(2),2)
     &              *NCmatrixelt(mod(if1,2),mod(if1,2))*polcolg
           enddo
        endif 
      endif 

      nmax = FL_VVg(iflav,-2)

      end
c******************************************************************
c
c   end subroutine qq_VVj2
c
c*****************************************************************
