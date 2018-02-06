c---------------------------------------------------------------------------
c	lepton helicities have to be summed over externally 
c		( this routine is for one fixed set of lep. helicities only ) 
c---------------------------------------------------------------------------
      
      subroutine qqhhqq(pbar,fsign, nlo, L,bos,
     1                  uucc,uuss,ddcc,ddss,udsc,ducs)
          use globalvars, only: ldoblha
      implicit none
c
c       Julien Baglio, <baglio@particle.uni-karlsruhe.de>
c       Created on: 21 Mai 2012 (to be cleaned!)  
c       based on the modification of ppzqq.f by
c	Barbara Jaeger, <baj@particle.physik.uni-karlsruhe.de>
c	Initial version:  2005 Sept. 02
c	Last modified:    2005 
C
C  QQHHQQ calculates the matrix elements**2 for
c  Higgs boson pair production by quark quark scattering
C
C        q1 q3    ---->   q2 q4 H H
C
c
C  This code is modified to allow for virtual corrections, more precisely
C  the interference of Born with the finite part of virtual diagrams
C  for 
c
c  INPUT:  NLO = 1       return uucc = |M_born|^2 + 2Re(M_Born^* M_virt)
c          NLO = 0       return uucc = |M_born|^2   etc.
c  integer L (=1,2) is used for the subtraction code to select leptonic 
c  tensors for gluon emission off the upper Line (L=1) or lower quark 
c  line (L=2)
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"
#include "tensorz.inc"
#include "tensorzn.inc"
c for check of individual graphs :
#include "tampz-lists.inc" 
c
c electroweak couplings are taken from KOPPLN
c
      double precision  clr, xm2, xmg, b, v, a
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     1                  V(4,5),A(4,5)
c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"
c
c variables for the main part of the program
c
      double precision  pbar(0:3,4+max_v)
      double precision uucc,uuss,ddcc,ddss,udsc,ducs,res(6),resv(6)
      double precision tree(6)
      double precision  p(0:3,4+max_v), p21(0:4), p43(0:4)

      integer  fsign(4+max_v), nlo, mu, i, k,
     1         isig1, isig3,bos
      integer  ifl(4,6), js1, js3, L, Lold, is1, is3
      double complex prop21(4), prop43(4)
      double complex mat(6,-1:1,-1:1,9), matv(6,-1:1,-1:1,9)
      double complex mm(6,-1:1,-1:1)

      double complex mzz, mzz5, mzz6

      double complex psi(2,-1:1,4), jqq(0:5,-1:1,2) !, eps(0:5)
      double complex jvez(0:5,-1:1,-1:1,4), jvea(0:5,-1:1,-1:1,4),
     1               jvuz(0:5,-1:1,-1:1,4), jvua(0:5,-1:1,-1:1,4)

      double precision dummy(0:4)
      double complex zm2i(2:4)

      double complex contract_Tjj, dotcc, dotrc, dotqj, s1c
      external contract_Tjj, dotcc, dotrc, dotqj, s1c
      logical ldebug, ldebugv, ldebugvf,linit,lerror,lgc(4)
      data linit /.true./, lerror /.false./, ldebug /.false./
      data lgc /4*.false./
      
      logical vvdebug,vbdebug, vlldebug,bbdebug,pdebug
      parameter (vvdebug = .false.,vbdebug = .false.,vlldebug =.false.,
     #		 bbdebug = .false.,pdebug = .false.)
      
      integer*8 icb1, icount1, icb2, icount2
      data icb1/0/,icount1/0/,icb2/0/,icount2/0/,Lold/1/

      parameter (ldebugv=.false., ldebugvf=.false.)
      save ifl, Lold, zm2i, linit, lgc, icb1, icb2, icount1, icount2
      double complex  zero
      parameter (zero = (0d0,0d0) )

!      integer iterations1, iterations2, n2max1, n2max2
!      integer lflavr(5:4+nv)
!      character*50 gridname1, gridname2
!      common /bmainpr/ iterations1, iterations2, n2max1, n2max2,
!     2                 lflavr, gridname1, gridname2
c
cc
c variables for virtual corrections
c
	logical bad_gauge,bad_gauge_sin,bad_gauge_isig
	common /gauge / bad_gauge,bad_gauge_sin,bad_gauge_isig
c
      double precision c2,c2o2pi
      parameter (c2=4d0/3d0, c2o2pi=c2/2d0/pi)
      double precision lnmusqp1p2, lnmusqp3p4

      double precision QCDfac
      logical lnlo
c      parameter (lnlo=.true.)
      lnlo = NLO.ne.0    ! include virtual stuff if T
c
	bad_gauge = .false. ! set F at beginning of each run
c
c
c define flavors of external quarks for the 4 NC and 2 CC subprocesses
c
      if (linit) then
         linit = .false.
         k = 1                  ! uucc
         ifl(1,k) = 3
         ifl(2,k) = 3
         ifl(3,k) = 3
         ifl(4,k) = 3
         k = 2                  ! uuss
         ifl(1,k) = 3
         ifl(2,k) = 3
         ifl(3,k) = 4
         ifl(4,k) = 4
         k = 3                  ! ddcc
         ifl(1,k) = 4
         ifl(2,k) = 4
         ifl(3,k) = 3
         ifl(4,k) = 3
         k = 4                  ! ddss
         ifl(1,k) = 4
         ifl(2,k) = 4
         ifl(3,k) = 4
         ifl(4,k) = 4
         k = 5                  ! udsc
         ifl(1,k) = 3
         ifl(2,k) = 4
         ifl(3,k) = 4
         ifl(4,k) = 3
         k =6                   ! ducs
         ifl(1,k) = 4
         ifl(2,k) = 3
         ifl(3,k) = 3
         ifl(4,k) = 4
         zm2i(2) = 1/dcmplx(xm2(2),-xmg(2))
         zm2i(3) = 1/dcmplx(xm2(3),-xmg(3))
         zm2i(4) = 1/dcmplx(xm2(4),-xmg(4))
      endif

 701  if (ldebug.or.ldebugv) then
         print*," " 
         print*," enter qqhhqq "
      endif
      do k = 1,6
         do isig1 = -1,1,2
            do isig3 = -1,1,2
               do i = 1,9
                  mat(k,isig1,isig3,i) = 0
                  matv(k,isig1,isig3,i) = 0
               enddo
            enddo
         enddo
      enddo
      do k = 0,5
         do isig1 = -1,1,2
            do isig3 = -1,1,2
               do i = 1,4
                  jvez(k,isig1,isig3,i) = 0
                  jvea(k,isig1,isig3,i) = 0
                  jvuz(k,isig1,isig3,i) = 0
                  jvua(k,isig1,isig3,i) = 0
               enddo
            enddo
         enddo
      enddo

c
c identify fermion line fsign factors
c
      is1 = fsign(1)
      is3 = fsign(3)
      js1 = (3+fsign(1))/2       ! 1 for fsign1=-1,2 for fsign1=+1
      js3 = (7+fsign(3))/2       ! 3 for fsign3=-1,4 for fsign3=+1
      if (ldebugv) then
         print*," ffsign = ",(fsign(i),i=1,4)
         print*," is1, is3 = ",is1,is3," js1, js3 = ",js1,js3
         print*," lzs = ",lzs," lzs1 = ",lzs(js1), " lzs3 = ",lzs(js3) 
         if (fsign(1).ne.fsign(2) .or. fsign(3).ne.fsign(4) ) then
            print*," quark fsign factors ",(fsign(i),i=1,4)
            print*," are not implemented yet in qqhhqq"
            stop
         endif
      endif

c
c define the internal momenta
c
      do mu = 0,3
         do i = 1,4+n_v
            p(mu,i) = pbar(mu,i)*fsign(i)
         enddo
	 
           p21(mu) =   p(mu,2) - p(mu,1)
           p43(mu) =   p(mu,4) - p(mu,3)	 
      enddo
      p21(4) = p21(0)**2 - p21(1)**2 - p21(2)**2 - p21(3)**2
      p43(4) = p43(0)**2 - p43(1)**2 - p43(2)**2 - p43(3)**2	 

      if (lwarn) then
         do mu=0,3
            dummy(mu) = p(mu,1)+p(mu,3)+p(mu,5)+p(mu,7)-
     1                 (p(mu,2)+p(mu,4)+p(mu,6)+p(mu,8))
         enddo
         if (abs(dummy(0)).gt.1d-8 .or.
     1       abs(dummy(1)).gt.1d-8 .or.
     2       abs(dummy(2)).gt.1d-8 .or.
     3       abs(dummy(3)).gt.1d-8 ) then
            print*," momentum mismatch in qqhhqq "
            print 102," dummy = ",(dummy(mu),mu=0,3)
            print 103," fsign  = ",fsign
 101        format( " p(", i1, ") = ", 4(f10.3, 2x) )
 102        format(a,4f10.3)
 103        format(a,6i5)
            do i = 1,8
               write(6,101) i, p(0,i), p(1,i), p(2,i), p(3,i)
            end do
         endif 
      endif
c
c get the vector boson propagator factors
c
      prop21(1) = 1/p21(4)
      prop21(2) = 1/dcmplx(p21(4)-xm2(2),xmg(2))
      prop21(3) = 1/dcmplx(p21(4)-xm2(3),xmg(3))
      prop21(4) = prop21(3)

      prop43(1) = 1/p43(4)
      prop43(2) = 1/dcmplx(p43(4)-xm2(2),xmg(2))
      prop43(3) = 1/dcmplx(p43(4)-xm2(3),xmg(3))
      prop43(4) = prop43(3)
c
c get the external quark spinors (including factor sqrt(2E) )
c
      call psi0m(4,pbar(0,1),fsign(1),psi)
c
c get the f-fbar currents J21^mu=jqq(mu,*,1), J43^mu=jqq(mu,*,2) 
c
      call curr6(1,psi(1,-1,2),p(0,2),psi(1,-1,1),p(0,1),jqq(0,-1,1))
      call curr6(1,psi(1,-1,4),p(0,4),psi(1,-1,3),p(0,3),jqq(0,-1,2))
      
c
c -------------------------------------------------------------------
c
c contract with VVtoHH tensors to get Vertex-Vertex scattering diagrams
c
      do isig1 = -1,1,2
         do isig3 = -1,1,2
            mzz =
     4           contract_Tjj(zzhh(0,0,L),jqq(0,isig1,1),jqq(0,isig3,2))
            do k = 1,4
               mat(k,isig1,isig3,1) = 
     4              mzz*clr(ifl(1,k),2,isig1)*clr(ifl(3,k),2,isig3)
               matv(k,isig1,isig3,1) = (0d0,0d0)
            enddo
         enddo !isig3
      enddo !isig1
      
c----------------------      
      
      mzz6 = contract_Tjj(wwhh6(0,0,L),jqq(0,-1,2),jqq(0,-1,1))
      mzz5 = contract_Tjj(wwhh5(0,0,L),jqq(0,-1,1),jqq(0,-1,2))
      mat(5,-1,-1,1) = mzz5*clr(3,3,-1)**2
      mat(6,-1,-1,1) = mzz6*clr(3,3,-1)**2
      matv(5,-1,-1,1) = (0d0,0d0)
      matv(6,-1,-1,1) = (0d0,0d0)
c      
c
c ------------------------------------------------------------------------------
c ------------------------------------------------------------------------------
     
c sum the graphs, square them and map them onto uucc, uuss etc.

c i = 1		VV

c for checking individual pieces only:
      if(.false.) then
         do isig1 = -1,1,2
            do isig3 = -1,1,2
               do  i = 1,8
                  do k = 1,6
                     mat(k,isig1,isig3,i) = 0d0
                  enddo	
               enddo
	    enddo
         enddo 
      endif
      
   
      do k = 1,6
         res(k) = 0
         resv(k) = 0
         tree(k) = 0D0
         do isig1 = -1,1,2
            do isig3 = -1,1,2
               mm(k,isig1,isig3) = 0
               do i = 1,1
                  mm(k,isig1,isig3) = mm(k,isig1,isig3) + 
     1                 mat(k,isig1,isig3,i)
               enddo
	       
               res(k) = res(k) + dreal(mm(k,isig1,isig3))**2
     &                         + dimag(mm(k,isig1,isig3))**2

            enddo
         enddo  

         if (lnlo) then
            QCDfac = ( 1d0+(als(1,1)+als(2,1))*c2o2pi*cvirt )
         else
            QCDfac = 1.d0
         endif
	       
         tree(k) = res(k) * 9D0 
         res(k) = res(k)*QCDfac*9d0   ! 9 is the color sum factor

      enddo

      if (ldebug) then
	print*,"uucc",uucc,res(1)/uucc
	print*,"uuss",uuss,res(2)/uuss
	print*,"ddcc",ddcc,res(3)/ddcc
	print*,"ddss",ddss,res(4)/ddss
	print*,"udsc",udsc,res(5)/udsc
	print*,"ducs",ducs,res(6)/ducs
       print*
      endif
  
      if (lwarn) then
         resv(1) = uucc
         resv(2) = uuss
         resv(3) = ddcc
         resv(4) = ddss
         resv(5) = udsc
         resv(6) = ducs
         do k = 1,6
            if ((res(k)/resv(k)-1).gt.1d-7) then
               if (.not.ldebug) then
                  print*," k =",k,":",res(k)/resv(k)-1
                  ldebug = .true.
                  lerror = .true.
                  goto 701
               endif
            endif
         enddo
         if (lerror) then
            print*," "
            lerror = .false.
            ldebug = .false.
         endif
      endif
c eliminate processes with photon virtuality below cutoff
      if ( abs(p21(4)).lt.qsqAmin .or.  
     1     abs(p43(4)).lt.qsqAmin ) then
         do k = 1,4
            res(k) = 0
            tree(k) = 0
         enddo
      endif


      uucc = res(1)
      uuss = res(2)
      ddcc = res(3)
      ddss = res(4)
      udsc = res(5)
      ducs = res(6)

      if (ldoblha) then
        if (blha_amptype(blha_curproc).eq.0) then ! tree
          blha_amp(1)=blha_amp(1)+
     &      res(blha_idsubproc(blha_cursubproc,blha_curproc))

        else if (blha_amptype(blha_curproc).eq.1) then ! loop
          lnmusqp1p2 = log(blha_scale**2/(-p21(4)))
          lnmusqp3p4 = log(blha_scale**2/(-p43(4)))

          blha_amp(1)=blha_amp(1)+
     &      blha_CF/(2d0*pi)*2d0*
     &       (als(1,1)+als(2,1))*(-1d0)*
     &       tree(blha_idsubproc(blha_cursubproc,blha_curproc))
          blha_amp(2)=blha_amp(2)+
     &      1d0/(2d0*pi)*2d0*
     &       ( als(1,1)*(-blha_gammaQuark-blha_CF*lnmusqp1p2)
     &        +als(2,1)*(-blha_gammaQuark-blha_CF*lnmusqp3p4) )*
     &       tree(blha_idsubproc(blha_cursubproc,blha_curproc))
c this is more complicated than necessary, 
c but should work for other VBF processes as well
          blha_amp(3)= blha_amp(3)+
     &      (res(blha_idsubproc(blha_cursubproc,blha_curproc))     ! tree + loop + I
     &       -tree(blha_idsubproc(blha_cursubproc,blha_curproc)))  ! -tree
     &       * blha_CF/c2                                          ! correct for Nc ?= inf
     &      + (-(als(1,1)+als(2,1))/2d0/pi*                        ! -I
     &          (blha_CF*cvirt+
     &            2d0*(blha_KQuark+blha_gammaQuark+
     &                 blha_tgammaQuark+(pi**2/6d0-1)*blha_CF) )
     &         +als(1,1)/pi*                                       ! +(mu^2/2p1.p2)^epsilon
     &         (-blha_gammaQuark*lnmusqp1p2-blha_CF/2d0*lnmusqp1p2**2)
     &         +als(2,1)/pi*                                       ! +(mu^2/2p3.p4)^epsilon
     &         (-blha_gammaQuark*lnmusqp3p4-blha_CF/2d0*lnmusqp3p4**2)
     &        )*tree(blha_idsubproc(blha_cursubproc,blha_curproc))
          blha_amp(4)= blha_amp(4)+
     &       tree(blha_idsubproc(blha_cursubproc,blha_curproc))

        else if (blha_amptype(blha_curproc).eq.2) then ! cctree
          ! 1-2
          call BLHA_cctree(1,2,
     &      -blha_CF*
     &       tree(blha_idsubproc(blha_cursubproc,blha_curproc))
     &    )
          ! 3-4
          call BLHA_cctree(3,4,
     &      -blha_CF*
     &       tree(blha_idsubproc(blha_cursubproc,blha_curproc))
     &    )

        else 
          call BLHA_amptypeerror(blha_amptype(blha_curproc),__FILE__,__LINE__)
        endif
      endif

      Lold = L

      return
      end


