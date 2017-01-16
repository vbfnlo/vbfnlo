c*****************************************************************************
c
c    begin function m2s_qqZZj2
c
c*****************************************************************************
c
c       Adapted from m2s_qqVVj2 
c
c
      real*8 function m2s_qqZZj2(
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
      INTEGER h
      COMMON /hcount / h

      complex*16 zero
      parameter (zero=(0d0,0d0))

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
      parameter (ldebug=.false.)          !output debug information

c variables for leshouches interface
      DOUBLE PRECISION weight,rnumb,RandomNumber

cmr
      double precision  clr, xm2, xmg, b2, v2, a2
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B2(6,6,6),
     1                  V2(4,5),A2(4,5)
cmr

      if ( init .eq. 0 ) then

         if (bos.eq.22) then
            write(6,*) " "
            write(6,*) "ZZj amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
         endif

         call printnfl(.false.)
         h=1
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
         qp(mu,1) = v(mu,1,1)+v(mu,2,1)   !Z/A_e
         qm(mu,1) = v(mu,3,1)+v(mu,4,1)   !Z/A_mu
         qaa(mu,1) = qp(mu,1)+qm(mu,1)    !Z_all
      enddo
      qp(4,1) = dotrr(qp(0,1),qp(0,1))
      qm(4,1) = dotrr(qm(0,1),qm(0,1))
      qaa(4,1) = dotrr(qaa(0,1),qaa(0,1))

c polarisation vectors
c select helicity: h = random number for helicity (h=1:2)
      h = mod(h,4) + 1
      ie = (-1)**(h+1)
      iu = (-1)**((h-1)/2)

c...Les Houches interface
        if ((lha.or.hepmc).and..not.doNLO) then
           helicity(1)=-ie
           helicity(2)= ie
           helicity(3)=-iu
           helicity(4)= iu
        endif

c for ZZ->4l precalculate leptonic tensors

      CALL OXXXXX(v(0,1,1),ZERO ,-ie,+1,wve(1,1)) !W(1,3))          !e-
      CALL IXXXXX(v(0,2,1),ZERO ,+ie,-1,wep(1,1)) !W(1,4))          !e+
      CALL OXXXXX(v(0,3,1),ZERO ,-iu,+1,wmu(1,1)) !W(1,5))          !mu-
      CALL IXXXXX(v(0,4,1),ZERO ,+iu,-1,wvm(1,1)) !W(1,6))          !mu+
      CALL JIOXXX(wep(1,1),wve(1,1),GZL ,ZMASS,ZWIDTH,wp(1,1))    !Z_e
      CALL JIOXXX(wep(1,1),wve(1,1),GAL ,ZERO ,ZERO  ,atau1(1,1)) !A_e
      CALL JIOXXX(wvm(1,1),wmu(1,1),GZL ,ZMASS,ZWIDTH,wm(1,1))    !Z_mu
      CALL JIOXXX(wvm(1,1),wmu(1,1),GAL ,ZERO ,ZERO  ,atau2(1,1)) !A_mu
      call zz_ztozz(v(0,1,1),1,wz(1,1),0) !Z -> e+ e- mu+ mu-
      call zz_atozz(v(0,1,1),1,wa(1,1),0) !A -> e+ e- mu+ mu-
      if ((nlo .eq. -8) .or. (nlo.eq.1)) 
     &  call zz_htozz(v(0,1,1),1,w4l(h,1),0) !H -> e+ e- mu+ mu-
      if (narrowwidth.ne.0) then
c no gamma* contributions
        do mu=1,6
          atau1(mu,1) = 0d0
          atau2(mu,1) = 0d0
        enddo
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


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     Now sum over the subprocesses contributing to ZZj production

      nmax = FL_VVg(fsign,-1)          !reset counter for subprocesses to 0

C*******************  q1 qbar2 ---> Z Z g   **********************
         
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

      call qq_ZZj2(xi,p,v,physToDiag,fsign,gsign,bos, nlo,
     &              pdf,res,nmin,nmax,ps_number)

      if (nlo .ne. -8) then
cc gg is symmetric -> would be double-counting
C*******************  qbar1 q2 ---> Z Z g   **********************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=1
      physToDiag(3)=3    ! gluon

      fsign(1) = 1
      fsign(2) = -1

      call qq_ZZj2(xi,p,v,physToDiag,fsign,gsign,bos, nlo,
     &              pdf,res,nmin,nmax,ps_number)

C*******************  g q1 ---> Z Z q2  ***********************************

      physToDiag(1)=3
      physToDiag(2)=1
      physToDiag(3)=2

      fsign(1) = 1
      fsign(2) = 1
      gsign    = -1

      call qq_ZZj2(xi,p,v,physToDiag,fsign,gsign,bos, nlo,
     &              pdf,res,nmin,nmax,ps_number)

C*******************  qbar2 g ---> Z Z qbar1   **********************

      physToDiag(1)=2
      physToDiag(2)=3
      physToDiag(3)=1

      fsign(1) = -1
      fsign(2) = -1
      gsign = -1

      call qq_ZZj2(xi,p,v,physToDiag,fsign,gsign,bos, nlo,
     &              pdf,res,nmin,nmax,ps_number)

C*******************  g q2bar ---> Z Z q1bar  **********************

      physToDiag(1)=3
      physToDiag(2)=2
      physToDiag(3)=1

      fsign(1) = -1
      fsign(2) = -1
      gsign = -1

      call qq_ZZj2(xi,p,v,physToDiag,fsign,gsign,bos, nlo,
     &              pdf,res,nmin,nmax,ps_number)

C*******************  q1 g ---> Z Z q2   **********************

      physToDiag(1)=1
      physToDiag(2)=3
      physToDiag(3)=2

      fsign(1) = 1
      fsign(2) = 1

      call qq_ZZj2(xi,p,v,physToDiag,fsign,gsign,bos, nlo,
     &              pdf,res,nmin,nmax,ps_number)

      endif
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

      m2s = m2s * 4d0      ! factor 4 for random helicity summation

      m2s_qqZZj2 = m2s
      
      RETURN
      END

c*****************************************************************************
c
c    end function m2s_qqZZj2
c
c*****************************************************************************

c******************************************************************
c
c   begin subroutine qq_ZZj2
c
c*****************************************************************

      subroutine qq_ZZj2(xi,p,v,physToDiag,fsign,gsign, bos, nlo,
     1                  pdf,res,nmin,nmax,ps_number)
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"

      real*8 p(0:3,max_p-1), v(0:3,n_v), xi(nx)
      real*8 pdf(-6:6,2), res(maxnumsubproc)
      integer physToDiag(3), fsign(n_p+n_v), gsign, 
     &        bos, nlo, nmin, nmax, ps_number
c
c qq_ZZj2 calls the amplitude square routines 
c             qqZZj    for qq --> ZZ --> 4 leptons
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
      real*8 uuzz, ddzz, uuzz2, ddzz2
      real*8 pbar(0:3,n_p+n_v), qbar(0:4), dotrr
      real*8 polcolq, polcolg, polcolgg

      double precision NCmatrixelt(0:1,0:1)
      double complex zero
      integer iflav(4), diagToPhys(3), FL_VVg
      external dotrr, FL_VVg


      logical ldebug
      parameter (ldebug=.false.)
      parameter (zero = (0d0,0d0) )
      integer i,if1,j,k,mu

      uuzz = 0d0
      uuzz2 = 0d0
      ddzz = 0d0
      ddzz2 = 0d0

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
        call GGZZg(pbar,fsign,qbar,gsign,bos,nlo,1,uuzz)
c switch back for nlo.eq.1
        fsign(2)=-fsign(2)
        iflav(1)=0
        iflav(2)=0
        k=FL_VVg(iflav,1)
  
        res(k)=pdf(sign1*iflav(1),1)
     &        *pdf(sign2*iflav(2),2)
     &        *uuzz*polcolgg

      endif
      if (nlo .ne. -8) then

#ifdef WITH_MADGRAPH
        if(ldebug) then
          call suu_zzg(pbar,fsign,uuzz2)
          call sdd_zzg(pbar,fsign,ddzz2)
        endif
#endif

      
        call QQZZj(pbar,fsign,qbar,gsign,nlo,1,uuzz,ddzz,ps_number)

#ifdef WITH_MADGRAPH
        if (ldebug) then
           if(abs(uuzz/uuzz2 - 1d0).ge.1.0d-5) then
              print*,"uuzz =   ",uuzz
              print*,"uuzz2 =  ",uuzz2
           endif
           if(abs(ddzz/ddzz2 - 1d0).ge.1.0d-5) then
              print*,"ddzz =   ",ddzz
              print*,"ddzz2 =  ",ddzz2
           endif
        endif
#endif

        NCmatrixelt(0,0)=uuzz
        NCmatrixelt(1,1)=ddzz


        if (gsign .eq. 1) then
           iflav(3) = 0           ! final state gluon id
           do if1=1,nfl
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
           do if1=1,nfl
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
c   end subroutine qq_ZZj2
c
c*****************************************************************
