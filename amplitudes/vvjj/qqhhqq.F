c---------------------------------------------------------------------------
c	lepton helicities have to be summed over externally 
c		( this routine is for one fixed set of lep. helicities only ) 
c---------------------------------------------------------------------------
      
      subroutine qqhhqq(pbar,fsign, nlo, L,bos,
     1                  uucc,uuss,ddcc,ddss,udsc,ducs)
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
      double precision  pbar(0:3,4+max_v), musq
      double precision uucc,uuss,ddcc,ddss,udsc,ducs,res(6),resv(6)
      double precision tree(6)
      double precision  p(0:3,4+max_v), p21(0:4), p43(0:4), pez(0:4),
     1                  puz(0:4), pzz(0:4)
      integer  fsign(4+max_v), nlo, mu, i, j, jj, k, kk, id,
     1         isig, isig1, isig3,bos,h
      integer  ifl(4,6), js1, js3, L, Lold, is1, is3
      double complex prop21(4), prop43(4)
      double complex mat(6,-1:1,-1:1,9), matv(6,-1:1,-1:1,9)
      double complex mm(6,-1:1,-1:1), 
     1               mv12(6,-1:1,-1:1), mv34(6,-1:1,-1:1)
      double complex mau,mzu
      double complex maa, maz, mza, mzz, mzz5, mzz6
      integer id1,id2,id3
      double complex m5e(3,2,2),m5u(3,2,2),
     1 	             m1p(6,2,2),m2p(6,2,2),m3p(6,2,2)
      double complex mv5e(3,2,2),mv5u(3,2,2),
     1 	             mv1p(6,2,2),mv2p(6,2,2),mv3p(6,2,2)
      double complex  m1, m2, ma(2), mz(2),
     1		      m1ae,m1ze,m2ae,m2ze,z1ae,z1ze,z2ae,z2ze,
     1	              m1au,m1zu,m2au,m2zu,z1au,z1zu,z2au,z2zu,
     1		      m3au,m3zu,m4au,m4zu,m3ae,m3ze,m4ae,m4ze,
     1		      z3ae,z3ze,z3zu,z3au,z4ae,z4ze,z4zu,z4au,
     1		      matu(6,-1:1,-1:1,9),mate(6,-1:1,-1:1,9),
     1		      mue,mate12,muu,matu12,mal,mzl,mate34,matu34,
     1		      muae,muze,mlae,mlze
      double complex  mv1ae,mv1ze,mv2ae,mv2ze,mv1au,mv1zu,mv2au,mv2zu,
     1		      mv3au,mv3zu,mv4au,mv4zu,mv3ae,mv3ze,mv4ae,mv4ze
      double complex  matve12,matvu12,matve34,matvu34
      double complex mv1aae,mv1aze,mv1zae,mv1zze,
     1 		     mv2aae,mv2aze,mv2zae,mv2zze,mv1e,mv2e,
     1		     mv1aau,mv1azu,mv1zau,mv1zzu,
     1 		     mv2aau,mv2azu,mv2zau,mv2zzu,mv1u,mv2u,
     1		     mv3aae,mv3aze,mv3zae,mv3zze,
     3 		     mv4aae,mv4aze,mv4zae,mv4zze,mv3e,mv4e,
     3		     mv3aau,mv3azu,mv3zau,mv3zzu,
     3 		     mv4aau,mv4azu,mv4zau,mv4zzu,mv3u,mv4u
      double complex m1aae,m1aze,m1zae,m1zze,z1aze,z1zze,
     1 		     m2aae,m2aze,m2zae,m2zze,m1e,m2e,
     1		     m1aau,m1azu,m1zau,m1zzu,z1azu,z1zzu,
     1 		     m2aau,m2azu,m2zau,m2zzu,m1u,m2u,
     1		     mz1aae,mz1aze,mz1zae,mz1zze,
     1		     mz2aae,mz2aze,mz2zae,mz2zze,
     1		     mz1aau,mz1azu,mz1zau,mz1zzu,
     1		     mz2aau,mz2azu,mz2zau,mz2zzu,
     1		     m3aae,m3aze,m3zae,m3zze,z3aze,z3zze,
     3 		     m4aae,m4aze,m4zae,m4zze,m3e,m4e,
     3		     m3aau,m3azu,m3zau,m3zzu,z3azu,z3zzu,
     3 		     m4aau,m4azu,m4zau,m4zzu,m3u,m4u,
     3		     mz3aae,mz3aze,mz3zae,mz3zze,
     3		     mz4aae,mz4aze,mz4zae,mz4zze,
     3		     mz3aau,mz3azu,mz3zau,mz3zzu,
     3		     mz4aau,mz4azu,mz4zau,mz4zzu
      double complex mezz(4,-1:1,-1:1),meza(4,-1:1,-1:1),
     #		     meaz(4,-1:1,-1:1),meaa(4,-1:1,-1:1),
     # 		     muzz(4,-1:1,-1:1),muza(4,-1:1,-1:1),
     #		     muaz(4,-1:1,-1:1),muaa(4,-1:1,-1:1),
     #		     mvezz(4,2,-1:1,-1:1),mveza(4,2,-1:1,-1:1),
     #		     mveaz(4,2,-1:1,-1:1),mveaa(4,2,-1:1,-1:1),
     # 		     mvuzz(4,2,-1:1,-1:1),mvuza(4,2,-1:1,-1:1),
     #		     mvuaz(4,2,-1:1,-1:1),mvuaa(4,2,-1:1,-1:1),
     #		     zezz(-1:1,-1:1,3),zeza(-1:1,-1:1,3),
     #		     zeaz(-1:1,-1:1,3),zeaa(-1:1,-1:1,3),
     # 		     zuzz(-1:1,-1:1,3),zuza(-1:1,-1:1,3),
     #		     zuaz(-1:1,-1:1,3),zuaa(-1:1,-1:1,3)
      double complex mv1,mv2,mva(2),mvz(2), mv5(3,3:4)
      double complex psi(2,-1:1,4), jqq(0:5,-1:1,2) !, eps(0:5)
      double complex psize(2,-1:1,4),psiae(2,-1:1,4),
     1 		     psizu(2,-1:1,4),psiau(2,-1:1,4),
     1               jez(0:5,-1:1,4),jea(0:5,-1:1,4),
     1		     juz(0:5,-1:1,4),jua(0:5,-1:1,4),
     1               jvez(0:5,-1:1,-1:1,4), jvea(0:5,-1:1,-1:1,4),
     1               jvuz(0:5,-1:1,-1:1,4), jvua(0:5,-1:1,-1:1,4)
      double complex bkjqq(2,-1:1,-1:1,4)
      double complex ja(0:5,-1:1,-1:1,4), jz(0:5,-1:1,-1:1,4), 
     1               jva(0:5,-1:1,-1:1,4), jvz(0:5,-1:1,-1:1,4),
     2               psia(2,-1:1,-1:1,4), psiz(2,-1:1,-1:1,4)
      double complex j5zzue(0:3,-1:1,3,-1:1,2), j5zzeu(0:3,-1:1,3,-1:1,2),
     1 		     j5zaue(0:3,-1:1,3,-1:1,2), j5zaeu(0:3,-1:1,3,-1:1,2),
     1 		     j5azue(0:3,-1:1,3,-1:1,2), j5azeu(0:3,-1:1,3,-1:1,2),
     1 		     j5aaue(0:3,-1:1,3,-1:1,2), j5aaeu(0:3,-1:1,3,-1:1,2)
      double complex epsee(0:3), epsuu(0:3),
     1		     epsCCuu(0:5,-1:1,2),epsCCee(0:5,-1:1,2),
     1		     epsee6(0:3), epsuu6(0:3),
     1		     epsCCuu6(0:5,-1:1,2),epsCCee6(0:5,-1:1,2)	 
      double complex epsaa(0:3),epsaz(0:3),epsza(0:3),epszz(0:3),
     1		     epsaee(0:5,-1:1,3:4,2), epszee(0:5,-1:1,3:4,2),
     1		     epsauu(0:5,-1:1,3:4,2), epszuu(0:5,-1:1,3:4,2),
     1		     qepszee(3:4),qepszuu(3:4)	 
      double precision fqze(0:4,4),fqae(0:4,4),fqzu(0:4,4),fqau(0:4,4)
      double precision fq(0:4,4),qee(0:4), quu(0:4), dummy(0:4)
      double complex zm2i(2:4)
      double complex jj21ae(-1:1),jj21au(-1:1),jj43ae(-1:1),jj43au(-1:1),
     1 		     jj21ze(-1:1),jj21zu(-1:1),jj43ze(-1:1),jj43zu(-1:1)
      double complex z1,zuu(4:5),zee(4:5)
      double complex qepsCCee(2), qepsCCuu(2),qepsCCee6(2), qepsCCuu6(2)
      double complex propt1(-1:1,-1:1,4,2), propt2(-1:1,-1:1,4,2),
     1		     propt(-1:1,-1:1,5:6,2), prop_ee(3), prop_uu(3),
     1		     propbbe, propbbu,
     1		     propbbea(6,-1:1,-1:1), propbbua(6,-1:1,-1:1),
     1		     propbbez(6,-1:1,-1:1), propbbuz(6,-1:1,-1:1)
      double complex fac
      integer idu,idl
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
      double precision xgc1, xgc2
      double complex mv1zz(-1:1),mv2zz(-1:1),mv1za(-1:1),mv2za(-1:1),
     1 		     mv1az(-1:1),mv2az(-1:1),mv1aa(-1:1),mv2aa(-1:1),
     1		     xv1zz(-1:1),xv2zz(-1:1),xv1za(-1:1),xv2za(-1:1),
     1		     xv1az(-1:1),xv2az(-1:1),xv1aa(-1:1),xv2aa(-1:1)
      parameter (ldebugv=.false., ldebugvf=.false.)
      save ifl, Lold, zm2i, linit, lgc, icb1, icb2, icount1, icount2
      save ja,jz,jva,jvz,psia,psiz,fq !,j5pm,j5mp
      save j5zzue, j5zzeu,j5zaue, j5zaeu,j5azue, j5azeu,j5aaue, j5aaeu
      double complex  zero
      parameter (zero = (0d0,0d0) )
      integer ii,ll
   
    

      integer icountmax
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


