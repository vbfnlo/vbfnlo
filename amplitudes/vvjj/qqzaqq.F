c		( this routine is for one fixed set of lep. helicities only ) 
c---------------------------------------------------------------------------
      
      subroutine qqzaqq(pbar,sign, nlo, L,bos, ps_number,
     1                  uucc,uuss,ddcc,ddss,udsc,ducs)
          use globalvars, only: ldoblha
      implicit none
c       Francisco Campanario francisco.campanario@ific.uv.es
c       Last modified: 11/2014
c       Taken from qqzzqq.F (Introduce new boxline and penlines routines) :
c	Barbara Jaeger, <baj@particle.physik.uni-karlsruhe.de>
c	Initial version:  2005 Sept. 02
c	Last modified:    2005 
C
C  QQZZQQ calculates the matrix elements**2 for electroweak
c  weak boson pair production by quark quark scattering
C
C        q1 q3    ---->   q2 q4 ZZ,   Z ---> f5-bar f6, Z ---> f7-bar f8
C
c  where the final state leptons can be 4 charged leptons (bos=22) or 
c  2 charged leptons + 2 neutrinos (bos=21).
C  Crossing related processes can be computed as well. Pauli interference terms for
c  identical fermions are neglected. In particular, only the
c  t-channel exchange of elctroweak bosons is considered. s-channel
c  production of 3 weak bosons is NOT implemented.
c
C  This code is modified to allow for virtual corrections, more precisely
C  the interference of Born with the finite part of virtual diagrams
C  for 
c
c  INPUT:  NLO = 1       return uucc = |M_born|^2 + 2Re(M_Born^* M_virt)
c          NLO = 0       return uucc = |M_born|^2   etc.
c   for steering the calculation of finite box and pentagon contributions
c   the following additional options are implemented
c          NLO = +4      set all finite box and pentagon contributions to 0
c                        but include cvirt*M_Born contribution
c          NLO = +5      true finite pentagons = 0, all other virtual contributions
c                  #included, including leading q^mu/m_V terms of 
c                        pentagons which can be related to boxes
c        remaining cases:    uucc = 2Re(M_Born^* M_virt), more precisely
c        the piece which is missing for the full NLO=1 when calling +|NLO|
c          NLO = -4      finite Born + q^mu/m_V terms of Pentagon contributions only, no cvirt
c          NLO = -5      true finite Pentagon contributions only
c
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
     1                  puz(0:4), pza(0:4)
      integer  sign(4+max_v), nlo, mu, i, j, jj, k, kk,
     1         isig, isig1, isig3,bos
      integer  ifl(4,6), js1, js3, L, Lold, is1, is3
      double complex prop21(4), prop43(4)
      double complex mat(6,-1:1,-1:1,9), matv(6,-1:1,-1:1,9)
      double complex mm(6,-1:1,-1:1), 
     1               mv12(6,-1:1,-1:1), mv34(6,-1:1,-1:1)

      double complex maa, maz, mza, mzz, mzz5, mzz6
      integer id1,id2,id3
      double complex m5e(3,2,2),m5u(3,2,2),
     1 	             m1p(6,2,2),m2p(6,2,2),m3p(6,2,2)
      double complex mv5e(3,2,2),mv5u(3,2,2),
     1 	             mv1p(6,2,2),mv2p(6,2,2),mv3p(6,2,2)
      double complex  ma(2), mz(2),
     1		      m1ae,m1ze,m2ae,m2ze,z1ae,z1ze,z2ae,z2ze,
     1	              m1au,m1zu,m2au,m2zu,z1au,z1zu,z2au,z2zu,
     1		      m3au,m3zu,m4au,m4zu,m3ae,m3ze,m4ae,m4ze,
     1		      z3ae,z3ze,z3zu,z3au,z4ae,z4ze,z4zu,z4au,
     1		      mate12,matu12,mate34,matu34

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
      double complex mva(2),mvz(2)
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
      double complex zuu(4:5),zee(4:5)
      double complex qepsCCee(2), qepsCCuu(2),qepsCCee6(2), qepsCCuu6(2)
      double complex propt1(-1:1,-1:1,4,2), propt2(-1:1,-1:1,4,2),
     1		     propt(-1:1,-1:1,5:6,2), prop_ee(3), prop_uu(3),
     1		     propbbe, propbbu,
     1		     propbbea(6,-1:1,-1:1), propbbua(6,-1:1,-1:1),
     1		     propbbez(6,-1:1,-1:1), propbbuz(6,-1:1,-1:1)
      double complex fac

      double complex contract_Tjj, dotcc, dotrc, dotqj, s1c, contract_Tjj_ck
      external contract_Tjj, dotcc, dotrc, dotqj, s1c,contract_Tjj_ck
      logical ldebug, ldebugv, ldebugvf,linit,lerror,lgc(4)
      data linit /.true./, lerror /.false./, ldebug /.false./
      data lgc /4*.false./
      
      logical vvdebug,vbdebug, vlldebug,bbdebug,pdebug
      parameter (vvdebug = .false.,vbdebug = .false.,vlldebug =.false.,
     #		 bbdebug = .false.,pdebug = .false.)
      
      integer*8 icb1, icount1, icb2, icount2
      data icb1/0/,icount1/0/,icb2/0/,icount2/0/,Lold/1/
      double precision xgc1, xgc2

      parameter (ldebugv=.false., ldebugvf=.false.)
      save ifl, Lold, zm2i, linit, lgc, icb1, icb2, icount1, icount2
      save ja,jz,jva,jvz,psia,psiz,fq !,j5pm,j5mp
      save j5zzue, j5zzeu,j5zaue, j5zaeu,j5azue, j5azeu,j5aaue, j5aaeu
      double complex  zero
      parameter (zero = (0d0,0d0) )
      integer ps_number


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
      double precision c2,c2o4pi, cvirtc
c      parameter (pi2o3=pi**2/3d0, cvirtc=pi2o3-7d0)
      parameter (c2=4d0/3d0, c2o4pi=c2/4d0/pi)
      double precision lnmusqp1p2, lnmusqp3p4
      logical lnlo, lbox, lpent, lpt
cfc
      logical debugza
      parameter (debugza=.false.)
cfc for virtuals
cfc call the new boxlines and penlines routines.
cfc add variables for wardtests
      double precision  pvirtcalc(0:3,5)
      double precision scalelowerline,scale, scaleupperline,temp

      DOUBLE COMPLEX mboxgauge(2,-1:1,0:3,0:2,12), ! check
     1   mboxgaugeborn(2,-1:1,0:3,0:2,12)!, ! check
      DOUBLE COMPLEX mpengauge(3,-1:1,0:3,0:2,12), ! check
     1   mpengaugeborn(3,-1:1,0:3,0:2,12)!, ! check

!     2     mboxborn(-1:1,-1:1,2,0:2),mbox(5,-1:1,-1:1,2,0:2)
      logical lwardtest(100),lwardtestp(100)
      integer wardIndex,div,divMax,divStep
      double complex  jez1(0:5,-1:1,4),jea1(0:5,-1:1,4),
     1		     jua1(0:5,-1:1,4),
     1               jvez1(0:5,-1:1,-1:1,4), jvea1(0:5,-1:1,-1:1,4),
     1              jvua1(0:5,-1:1,-1:1,4)
      double complex jbzaeu1(0:3,-1:1,3,2), j5zaeu1(0:3,-1:1,3,-1:1,2),jbzaeu(0:3,-1:1,3,2)
      double complex jbaaeu1(0:3,-1:1,3,2), j5aaeu1(0:3,-1:1,3,-1:1,2),jbaaeu(0:3,-1:1,3,2)
      double complex jbazue1(0:3,-1:1,3,2), j5azue1(0:3,-1:1,3,-1:1,2),jbazue(0:3,-1:1,3,2)
      double complex jbaaue1(0:3,-1:1,3,2), j5aaue1(0:3,-1:1,3,-1:1,2),jbaaue(0:3,-1:1,3,2)
      double complex ja1(0:5,-1:1,-1:1,4), jz1(0:5,-1:1,-1:1,4), 
     1               jva1(0:5,-1:1,-1:1,4), jvz1(0:5,-1:1,-1:1,4)

      logical lwardbox,lwardpent
      integer wardindexp
      logical lvirtratio
      parameter (lvirtratio=.false.)
      double precision K_contcheck 
      double complex ratio
c     lvirtratio: box- and penlines compute the 1/eps and 1/eps^2 pole 
c     terms of the virt. contributions
c     these terms are proportional to the born: mvirt(div=i)/mborn = const

      wardindex=1
      wardindexp=1
      lwardBox = .true.
      lwardPent = .true.
      divMax = 0
      divStep = 2


c      parameter (lnlo=.true.)
      lnlo = NLO.ne.0    ! include some virtual stuff if T
      lbox = NLO.eq.1 .or. NLO.eq.-4 .or. NLO.eq.5    ! call boxline if T
      lpt  = NLO.eq.1 .or. NLO.le.-5                  ! true pentagon contributions on
!      lpq  = NLO.eq.5 .or. NLO.eq.-4                  ! q^mu/m_V terms of Pentagon contributions only
      lpent= lpt !.or. lpq                             ! include pentagons if T
c
	bad_gauge = .false. ! set F at beginning of each run
c
c
c define flavors of external quarks for the 4 NC and 2 CC subprocesses
c
      if(lvirtratio) then
         divMax = 2
         divStep = 1
         lbox = .true.
         lpent = .true.
         cvirtc = 0d0
         call InitGaugeTest
      endif

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
c identify fermion line sign factors
c
      is1 = sign(1)
      is3 = sign(3)
      js1 = (3+sign(1))/2       ! 1 for sign1=-1,2 for sign1=+1
      js3 = (7+sign(3))/2       ! 3 for sign3=-1,4 for sign3=+1

c
c define the internal momenta
c
      do mu = 0,3
         do i = 1,4+n_v
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo
	 
	 if(bos.eq.212) then
           p21(mu) =   p(mu,2) - p(mu,1)
           p43(mu) =   p(mu,4) - p(mu,3)
cfc note pez should have change to puz, becuase is neutrino  Keep original notation to avoid further changes
           pez(mu) =   -(p(mu,6) - p(mu,5))
cfc note that puz, in ms2 is called qa.
cfc let-s keep the original name to avoid including
c furhter changes
           puz(mu) =   p(mu,7)  !cfc-(p(mu,8) - p(mu,7))
           pza(mu) =   pez(mu) + puz(mu)	 	 

	 elseif(bos.eq.211) then
           p21(mu) =   p(mu,2) - p(mu,1)
           p43(mu) =   p(mu,4) - p(mu,3)
           pez(mu) =   p(mu,6) - p(mu,5)
cfc note that puz, in ms2 is called qa.
cfc let-s keep the original name to avoid including
           puz(mu) =  p(mu,7)  !cfcp(mu,8) - p(mu,7)
           pza(mu) =   pez(mu) + puz(mu)
	 endif	 
	 
      enddo
      p21(4) = p21(0)**2 - p21(1)**2 - p21(2)**2 - p21(3)**2
      p43(4) = p43(0)**2 - p43(1)**2 - p43(2)**2 - p43(3)**2
      pez(4) = pez(0)**2 - pez(1)**2 - pez(2)**2 - pez(3)**2
      puz(4) = puz(0)**2 - puz(1)**2 - puz(2)**2 - puz(3)**2
      pza(4) = pza(0)**2 - pza(1)**2 - pza(2)**2 - pza(3)**2

cfc TODO
         scaleupperline = -1d0* p21(4)
         scalelowerline = -1d0* p43(4)
cfc
      if(lbox.or.lpent) then
         scaleupperline = -1d0* p21(4)
         scalelowerline = -1d0* p43(4)
      endif
	 

      if (lwarn) then
         do mu=0,3
c   changed sign of p(mu,7) since now it is an A not a l+
            dummy(mu) = p(mu,1)+p(mu,3)+p(mu,5)-p(mu,7)- 
     1                 (p(mu,2)+p(mu,4)+p(mu,6)+p(mu,8))
         enddo
         if (abs(dummy(0)).gt.1d-8 .or.
     1       abs(dummy(1)).gt.1d-8 .or.
     2       abs(dummy(2)).gt.1d-8 .or.
     3       abs(dummy(3)).gt.1d-8 ) then
            print*," momentum mismatch in qqzzqq "
            print 102," dummy = ",(dummy(mu),mu=0,3)
            print 103," sign  = ",sign
 101        format( " p(", i1, ") = ", 4(f10.3, 2x) )
 102        format(a,4f10.3)
 103        format(a,6i5)
            do i = 1,7
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
c for box-box and BV graphs we need the propagators for t-channel bosons between quark lines
c as seen from upper line these W momenta are INCOMING. They are OUTGOINg as seen from lower 
c line
      do mu = 0,3
         qee(mu) = pez(mu)+p(mu,2) - p(mu,1)    ! e+e-  emitted on upper line
         quu(mu) = puz(mu)+p(mu,2) - p(mu,1)    !mu+mu- emitted on upper line
      enddo
      qee(4) = qee(0)**2-qee(1)**2-qee(2)**2-qee(3)**2
      quu(4) = quu(0)**2-quu(1)**2-quu(2)**2-quu(3)**2

      prop_ee(1) = 1d0/qee(4)
      prop_ee(2) = 1d0/dcmplx(qee(4)-xm2(2),xmg(2))
      prop_ee(3) = 1d0/dcmplx(qee(4)-xm2(3),xmg(3))
      prop_uu(1) = 1d0/quu(4)
      prop_uu(2) = 1d0/dcmplx(quu(4)-xm2(2),xmg(2))
      prop_uu(3) = 1d0/dcmplx(quu(4)-xm2(3),xmg(3))
      
c
c get the external quark spinors (including factor sqrt(2E) )
c
      call psi0m(4,pbar(0,1),sign(1),psi)
c
c get the f-fbar currents J21^mu=jqq(mu,*,1), J43^mu=jqq(mu,*,2) 
c
      call curr6(1,psi(1,-1,2),p(0,2),psi(1,-1,1),p(0,1),jqq(0,-1,1))
      call curr6(1,psi(1,-1,4),p(0,4),psi(1,-1,3),p(0,3),jqq(0,-1,2))

      
c
c -------------------------------------------------------------------
c
c contract with vvtoww tensors to get Vertex-Vertex scattering diagrams
c
      do isig1 = -1,1,2
         do isig3 = -1,1,2
            maa = contract_Tjj(aazz(0,0,L),jqq(0,isig1,1),jqq(0,isig3,2))
            maz = contract_Tjj(azzz(0,0,L),jqq(0,isig1,1),jqq(0,isig3,2))
            mza = contract_Tjj(zazz(0,0,L),jqq(0,isig3,2),jqq(0,isig1,1))
            mzz = contract_Tjj(zzzz(0,0,L),jqq(0,isig1,1),jqq(0,isig3,2))

            do k = 1,4
               mat(k,isig1,isig3,1) = 
     1              maa*clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3)
     2            + maz*clr(ifl(1,k),1,isig1)*clr(ifl(3,k),2,isig3)
     3            + mza*clr(ifl(1,k),2,isig1)*clr(ifl(3,k),1,isig3)
     4            + mzz*clr(ifl(1,k),2,isig1)*clr(ifl(3,k),2,isig3)
cfc Virtuals ZERO. They factorize to the born
! TODO POLES
               matv(k,isig1,isig3,1) = (0d0,0d0)
            enddo
         enddo !isig3
      enddo !isig1
      
c----------------------      
      
      mzz6 = contract_Tjj(wwzz6(0,0,L),jqq(0,-1,2),jqq(0,-1,1))
      mzz5 = contract_Tjj(wwzz5(0,0,L),jqq(0,-1,1),jqq(0,-1,2))
      mat(5,-1,-1,1) = mzz5*clr(3,3,-1)**2
      mat(6,-1,-1,1) = mzz6*clr(3,3,-1)**2

cfc  Virtuals ZERO. They factorize to the born
! TODO POLES
      matv(5,-1,-1,1) = (0d0,0d0)
      matv(6,-1,-1,1) = (0d0,0d0)
c
c -------------------------------------------------------------------
c
c  prepare box diagrams: attach A/Z to external spinors
c 
c      isig = +-1   : left- and righthanded spinors coupling to A/Z"s
c
c  Notation for virtual 2-component spinors and momenta
c
c  Z -> e+e- attached to quark number i: psize(*,isig,i) with momentum fqze(mu,i)
c  Z->mu+mu- attached to quark number i: psizu(*,isig,i) with momentum fqzu(mu,i)
c  A -> e+e- attached to quark number i: psiae(*,isig,i) with momentum fqae(mu,i)
c  A->mu+mu- attached to quark number i: psiau(*,isig,i) with momentum fqau(mu,i)
c  
c  the fermion current corresponding to a quark line with the real emitted V-> l+l-
c  attached next to quark number i is stored in jlv(mu,isig,i) etc. 
c 
c  For the virtual amlitudes the notation, e.g.     jvlv(mu,isig,is,i)
c  is used for the boxline correction to a quark line with one V attached next to
c  quark #i and a free Lorentz index mu for the second attached EW boson. is=+-1
c  refers to the sign factor of this quark (vs.antiquark line). They are recalculated 
c  only if this quark line sign has not been calculated yet for this phase space point
c  (i.e. lzs = .true.)  Otherwise they are taken from saved previous calculation
 
      do i = 1,3,2
        do isig = -1,1,2
      
         call ket2c(psi(1,isig,i),.true.,p(0,i),isig,qe,ze,
     1              psize(1,isig,i),fqze(0,i))
         call ket2c(psi(1,isig,i),.true.,p(0,i),isig,qe,ae,
     1              psiae(1,isig,i),fqae(0,i))
cfc
         if(debugza)then
         call ket2c(psi(1,isig,i),.true.,p(0,i),isig,qu,zu,
     1              psizu(1,isig,i),fqzu(0,i))
         endif
cfc
         call ket2c(psi(1,isig,i),.true.,p(0,i),isig,qu,au,
     1              psiau(1,isig,i),fqau(0,i))
        

	 call bra2c(psi(1,isig,i+1),.true.,p(0,i+1),isig,qe,ze,
     1              psize(1,isig,i+1),fqze(0,i+1))
         call bra2c(psi(1,isig,i+1),.true.,p(0,i+1),isig,qe,ae,
     1              psiae(1,isig,i+1),fqae(0,i+1))
cfc
         if(debugza)then
         call bra2c(psi(1,isig,i+1),.true.,p(0,i+1),isig,qu,zu,
     1              psizu(1,isig,i+1),fqzu(0,i+1))
         endif
cfc
         call bra2c(psi(1,isig,i+1),.true.,p(0,i+1),isig,qu,au,
     1              psiau(1,isig,i+1),fqau(0,i+1))

         enddo !isig
        
         call curr6(1,psi(1,-1,i+1),p(0,i+1),
     1                 psize(1,-1,i),fqze(0,i), jez(0,-1,i)   )
         call curr6(1,psi(1,-1,i+1),p(0,i+1),
     1                 psiae(1,-1,i),fqae(0,i), jea(0,-1,i)   )
         if(debugza)then
         call curr6(1,psi(1,-1,i+1),p(0,i+1),
     1                 psizu(1,-1,i),fqzu(0,i), juz(0,-1,i)   )
         endif
         call curr6(1,psi(1,-1,i+1),p(0,i+1),
     1                 psiau(1,-1,i),fqau(0,i), jua(0,-1,i)   )

         call curr6(1,psize(1,-1,i+1),fqze(0,i+1),
     1                 psi(1,-1,i),p(0,i), jez(0,-1,i+1) )
         call curr6(1,psiae(1,-1,i+1),fqae(0,i+1),
     1                 psi(1,-1,i),p(0,i), jea(0,-1,i+1) )
         if(debugza)then
         call curr6(1,psizu(1,-1,i+1),fqzu(0,i+1),
     1                 psi(1,-1,i),p(0,i), juz(0,-1,i+1) )
         endif
         call curr6(1,psiau(1,-1,i+1),fqau(0,i+1),
     1                 psi(1,-1,i),p(0,i), jua(0,-1,i+1) )
     
   
       enddo !i


#ifdef WITH_NLO  
       if (lbox) then
         if(.not. lzs(js1)) then
            print*, "STOP in qqzaqq.F"
         ENDIF  
      endif

      if (lbox .and. lzs(js1)) then

c        if(.true.)then
         do i = 1,2 ! Two permutations of momenta
c Notation is:  quirality, psi1,psi4,1,p4, momenta of one boson ez
c    true to compute tensor integrals
c 3-i what is left as a curretn
c jez born 
c jvez virtual
c  output:   jvez, jvea,jvuz,jvua,(llamar para lineas de arriba y de abajo)
cfc           print*,"psi(1,-1,1)",psi(:,-1,1)
cfc #cfc Convert to new notation. 
cfc #c Convert to internal momenta for boxlines
cfc #c momenta of the particles
cfc #c quark 1:   p(i,1)
cfc #c quark 2:   p(i,2)
cfc #c Boson Ze momenta   qe, for ze and ae
cfc #c
cfc #c   NOTE the if statement below there are three do before
cfc #c         to optimize this. Move somewhere else
cfc #c
cfc #c To build currents to contract
cfc #c end currents to contract
cfc #c psi echange in normal notation is the entries would be:
cfc #c         psi4,psi1,p1,p4,p2,p3
cfc #c photon momenta ea
cfc #c ze: current z-> l+l-         
cfc #c ae: current a-> l+l-
cfc #c      DOUBLE COMPLEX mboxgauge(2,-1:1,-1:1,2,0:2), 
cfc #c     1     mbox(5,-1:1,-1:1,2,0:2),mboxgaugeborn(2,-1:1,-1:1,2,0:2), 
cfc #c     2     mboxborn(-1:1,-1:1,2,0:2)
cfc #c mbox(results,helicities_fermion,helicity_boson,permutation,poles)
cfc #c mbox(helicities_fermion,helicity_boson,permutation,poles)
cfc #
cfc end conversion
cfc  So the best is to search fo all the jvez, jvea output and search where
cfc  is contracted with. Then, take the polarization vectors from there and
cfc   start to build the results. Then compare.
cfc probably the cvirt should be removed.

c  output:   jvez, jvea,jvuz,jvua,(llamar para lineas de arriba y de abajo)


cfc            call boxlinec(1,psi(1,-1,1),psi(1,-1,2),p(0,1),p(0,2),
cfc     1                    .true., ze,3-i,   jez(0,-1,i), 
cfc     2                    jvez(0,-1,is1,i) )
cfc
cfc            call boxlinec(1,psi(1,-1,1),psi(1,-1,2),p(0,1),p(0,2),
cfc     1                    .true., ae,3-i,   jea(0,-1,i), 
cfc     2                    jvea(0,-1,is1,i) )

         do mu = 0,3
            pvirtcalc(mu,1) = p(mu,1)
            pvirtcalc(mu,2) = -1d0 * p(mu,2)
            pvirtcalc(mu,3) = -1d0 * qe(mu)
            pvirtcalc(mu,4) = 1d0 * (p(mu,3) - p(mu,4) - qa(mu))
            temp=pvirtcalc(mu,1)+pvirtcalc(mu,2)+pvirtcalc(mu,3)+pvirtcalc(mu,4)
c            print*,"Mom. conserv.", temp
         enddo
           scale = scaleupperline
          call BOX_VBF_CURR(!
     $      pvirtcalc(0,1), pvirtcalc(0,3), ! momenta 
     $      pvirtcalc(0,4), pvirtcalc(0,2),  ! momenta
     $      psi(1,-1,2),psi(1,-1,1),! spinors
     $      ze, jvez1(0,-1,is1,i),jez1(0,-1,i), ! currents
     $      1, ! Max helicity 
     $      1, 3-i, ! comp, number of current
     $      scale,divmax,   ! Max helicity(1 spinor 2 ext. boson)
     $      mboxgauge(1,-1,0,0,1), ! virtuals gauge
     $      mboxgaugeborn(1,-1,0,0,1),!born gauge
     $      lwardtest(wardIndex)) ! Gauge test

          wardindex=wardindex+1
cfc          stop
          if(lwardtest(wardIndex-1)) then
          call BOX_VBF_CURR(!
     $      pvirtcalc(0,1), pvirtcalc(0,3), ! momenta 
     $      pvirtcalc(0,4), pvirtcalc(0,2),  ! momenta
     $      psi(1,-1,2),psi(1,-1,1),! spinors
     $      ae, jvea1(0,-1,is1,i),jea1(0,-1,i), ! currents
     $      1, ! Max helicity 
cfc NOTE: COMP=-1
     $      -1, 3-i, ! comp, number of current
cfc ************************************** 
     $      scale,divmax,   ! Max helicity(1 spinor 2 ext. boson)
     $      mboxgauge(1,-1,0,0,2), ! virtuals gauge
     $      mboxgaugeborn(1,-1,0,0,2),!born gauge
     $      lwardtest(wardIndex)) ! Gauge test

          wardindex=wardindex+1
          else
          wardindex=wardindex+1
           jvea1=0d0
           jea1 = 0d0
          endif

cfc NOTE: There is not Z-> second ll. Photon on-shell
CFC
cfc            call boxlinec(1,psi(1,-1,1),psi(1,-1,2),p(0,1),p(0,2),
cfc     1                    .true., zu,3-i,   juz(0,-1,i), 
cfc     2                    jvuz(0,-1,is1,i) )
c
c            print*, "zu", zu
            juz=0d0
            jvuz=0d0
c	
c there are A->l+l-, but no A->vm vm~ contributions:
cfc            call boxlinec(1,psi(1,-1,1),psi(1,-1,2),p(0,1),p(0,2),
cfc     1                    .true., au,3-i,   jua(0,-1,i), 
cfc     2                    jvua(0,-1,is1,i) )

         do mu = 0,3
            pvirtcalc(mu,1) = p(mu,1)
            pvirtcalc(mu,2) = -1d0 * p(mu,2)
            pvirtcalc(mu,3) = -1d0 * qa(mu)
            pvirtcalc(mu,4) = 1d0 * (p(mu,3) - p(mu,4) - qe(mu))
            temp=pvirtcalc(mu,1)+pvirtcalc(mu,2)+pvirtcalc(mu,3)+pvirtcalc(mu,4)
c            print*,"Mom. conserv.", temp
         enddo
           scale = scaleupperline
          call BOX_VBF_CURR(!
     $      pvirtcalc(0,1), pvirtcalc(0,3), ! momenta 
     $      pvirtcalc(0,4), pvirtcalc(0,2),  ! momenta
     $      psi(1,-1,2),psi(1,-1,1),! spinors
     $      au, jvua1(0,-1,is1,i),jua1(0,-1,i), ! currents
     $      1, ! Max helicity 
     $      1, 3-i, ! comp, number of current
     $      scale,divmax,   ! Max helicity(1 spinor 2 ext. boson)
     $      mboxgauge(1,-1,0,0,3), ! virtuals gauge
     $      mboxgaugeborn(1,-1,0,0,3),!born gauge
     $      lwardtest(wardIndex)) ! Gauge test

          wardindex=wardindex+1
          enddo !i
      endif !lbox

      if (lbox) then
         if(.not. lzs(js3)) then
            print*, "STOP2 in qqzaqq.F"
        ENDIF
      endif

      if (lbox .and. lzs(js3)) then

c      if(.true.) then
         do i = 3,4	 
cfc            call boxlinec(1,psi(1,-1,3),psi(1,-1,4),p(0,3),p(0,4),
cfc     1                    .true., ze,5-i,   jez(0,-1,i), 
cfc     2                    jvez(0,-1,is3,i) )
cfc            call boxlinec(1,psi(1,-1,3),psi(1,-1,4),p(0,3),p(0,4),
cfc     1                    .true., ae,5-i,   jea(0,-1,i), 
cfc     2                    jvea(0,-1,is3,i) )


         do mu = 0,3
            pvirtcalc(mu,1) = p(mu,3)
            pvirtcalc(mu,2) = -1d0 * p(mu,4)
            pvirtcalc(mu,3) = -1d0 * qe(mu)
            pvirtcalc(mu,4) = 1d0 * (p(mu,1) - p(mu,2) - qa(mu))
            temp=pvirtcalc(mu,1)+pvirtcalc(mu,2)+pvirtcalc(mu,3)+pvirtcalc(mu,4)
c            print*,"Mom. conserv.", temp
         enddo
           scale = scalelowerline
          call BOX_VBF_CURR(!
     $      pvirtcalc(0,1), pvirtcalc(0,3), ! momenta 
     $      pvirtcalc(0,4), pvirtcalc(0,2),  ! momenta
     $      psi(1,-1,4),psi(1,-1,3),! spinors
     $      ze, jvez1(0,-1,is3,i),jez1(0,-1,i), ! currents
     $      1, ! Max helicity 
     $      1, 5-i, ! comp, number of current
     $      scale,divmax,   ! Max helicity(1 spinor 2 ext. boson)
     $      mboxgauge(1,-1,0,0,4), ! virtuals gauge
     $      mboxgaugeborn(1,-1,0,0,4),!born gauge
     $      lwardtest(wardIndex)) ! Gauge test

          wardindex=wardindex+1
          if(lwardtest(wardIndex-1)) then
          call BOX_VBF_CURR(!
     $      pvirtcalc(0,1), pvirtcalc(0,3), ! momenta 
     $      pvirtcalc(0,4), pvirtcalc(0,2),  ! momenta
     $      psi(1,-1,4),psi(1,-1,3),! spinors
     $      ae, jvea1(0,-1,is3,i),jea1(0,-1,i), ! currents
     $      1, ! Max helicity 
cfc NOTE: COMP=-1
     $      -1, 5-i, ! comp, number of current
cfc ************************************** 
     $      scale,divmax,   ! Max helicity(1 spinor 2 ext. boson)
     $      mboxgauge(1,-1,0,0,5), ! virtuals gauge
     $      mboxgaugeborn(1,-1,0,0,5),!born gauge
     $      lwardtest(wardIndex)) ! Gauge test

          wardindex=wardindex+1
          else
          wardindex=wardindex+1
          jvea1=0d0
          jea1=0d0
          endif
cfc
            juz=0d0
            jvuz=0d0
CFC
c there are A on-shel
cfc            call boxlinec(1,psi(1,-1,3),psi(1,-1,4),p(0,3),p(0,4),
cfc     1                    .true., au,5-i,   jua(0,-1,i), 
cfc     2                    jvua(0,-1,is3,i) )


         do mu = 0,3
            pvirtcalc(mu,1) = p(mu,3)
            pvirtcalc(mu,2) = -1d0 * p(mu,4)
            pvirtcalc(mu,3) = -1d0 * qa(mu)
            pvirtcalc(mu,4) = 1d0 * (p(mu,1) - p(mu,2) - qe(mu))
            temp=pvirtcalc(mu,1)+pvirtcalc(mu,2)+pvirtcalc(mu,3)+pvirtcalc(mu,4)
c            print*,"Mom. conserv.", temp
         enddo
           scale = scalelowerline
          call BOX_VBF_CURR(!
     $      pvirtcalc(0,1), pvirtcalc(0,3), ! momenta 
     $      pvirtcalc(0,4), pvirtcalc(0,2),  ! momenta
     $      psi(1,-1,4),psi(1,-1,3),! spinors
     $      au, jvua1(0,-1,is3,i),jua1(0,-1,i), ! currents
     $      1, ! Max helicity 
     $      1, 5-i, ! comp, number of current
     $      scale,divmax,   ! Max helicity(1 spinor 2 ext. boson)
     $      mboxgauge(1,-1,0,0,6), ! virtuals gauge
     $      mboxgaugeborn(1,-1,0,0,6),!born gauge
     $      lwardtest(wardIndex)) ! Gauge test

          wardindex=wardindex+1
        
         enddo !i
        jvez = jvez1
        jvea = jvea1
        jvua = jvua1
      endif !lbox
#endif

      IF(debugza) then
      do i=1,4
         if(i.lt.3) then
            isig3=is1
          else
             isig3=is3
         endif
c         print*,"k",i
         do isig1=-1,1,2
c            print*, "isig1",isig1
            do mu=0,5
             if(lbox) then
          jez(mu,isig1,i)        = jez1(mu,isig1,i)
          jea(mu,isig1,i)        = jea1(mu,isig1,i)
          jua(mu,isig1,i)        =  jua1(mu,isig1,i)
                endif

            if((abs(jvez(mu,isig1,isig3,i)/jvez1(mu,isig1,isig3,i)-1).gt.1d-3).and.(abs(jvez(mu,isig1,isig3,i)).ne.0d0)) then
            print*, " Z -> ll new result",mu,isig1,isig3,i
            print*,"jvezR",jvez(mu,isig1,isig3,i)/jvez1(mu,isig1,isig3,i)-1
            print*," jezR", jez(mu,isig1,i)/jez1(mu,isig1,i)-1
            print*," jez1", jez1(mu,isig1,i)
            print*," jez", jez(mu,isig1,i)
            print*
            print*
            endif
            if((abs(jvea(mu,isig1,isig3,i)/jvea1(mu,isig1,isig3,i)-1).gt.1d-3).and.(abs(jvea(mu,isig1,isig3,i)).ne.0d0)) then
            print*, " A -> ll new result",mu,isig1,isig3,i
            print*,"jvea2",jvea(mu,isig1,isig3,i)/jvea1(mu,isig1,isig3,i)-1
            print*," juea2", jea(mu,isig1,i)/jea1(mu,isig1,i)-1
            print*
            print*
            endif
            if((abs(jvua(mu,isig1,isig3,i)/jvua1(mu,isig1,isig3,i)-1).gt.1d-3).and.(abs(jvua(mu,isig1,isig3,i)).ne.0d0)) then
            print*, "Photon new result",mu,isig1,isig3,i
            print*,"jvua2",jvua(mu,isig1,isig3,i)/jvua1(mu,isig1,isig3,i)-1
            print*," jua21", jua(mu,isig1,i)/jua1(mu,isig1,i)-1
            print*
            print*
            endif
         enddo
        enddo
      enddo
      endif
c
c
c -------------------------------------------------------------------
c -------------------------------------------------------------------
c
c now calculate the Vertex-box diagrams; get t-channel W currents first
      do k = 1,2
         kk = 3-k
c k=1 and kk=2 is for "box correction" to upper line
c k=2 and kk=1 is for "box correction" to lower line
         if (k.eq.1) then
            zuu(4) = -dcmplx(quu(0),quu(3))
            zuu(5) = -dcmplx(quu(1),quu(2))
            zee(4) = -dcmplx(qee(0),qee(3))
            zee(5) = -dcmplx(qee(1),qee(2))
         else
            zee(4) = dcmplx(quu(0),quu(3))
            zee(5) = dcmplx(quu(1),quu(2))
            zuu(4) = dcmplx(qee(0),qee(3))
            zuu(5) = dcmplx(qee(1),qee(2))
         endif

c NC first:

	do isig = -1,1,2
	
	! VV-> e+e-: (V = A or Z)
	
	  if (k.eq.1) then
	    call contract_T2j(aaee(0,0,k,L),jqq(0,isig,kk),epsaa) !for AA 
	    call contract_T2j(azee(0,0,k,L),jqq(0,isig,kk),epsaz) !for AZ 
	    call contract_T2j(zaee(0,0,k,L),jqq(0,isig,kk),epsza) !for ZA 
	    call contract_T2j(zzee(0,0,k,L),jqq(0,isig,kk),epszz) !for ZZ 
	  else
 	    call contract_T1j(aaee(0,0,k,L),jqq(0,isig,kk),epsaa) !for AA 
	    call contract_T1j(azee(0,0,k,L),jqq(0,isig,kk),epsza) !for AZ 
	    call contract_T1j(zaee(0,0,k,L),jqq(0,isig,kk),epsaz) !for ZA 
	    call contract_T1j(zzee(0,0,k,L),jqq(0,isig,kk),epszz) !for ZZ 
	  endif
	
           do mu = 0,3 ! graphs with uncontracted A (epsaee) or Z (epszee),
	    		! coupling to up(3)- or down(4)-type quark 
               epsaee(mu,isig,3,k) =
     1              epsaa(mu)*clr(3,1,isig)+epsaz(mu)*clr(3,2,isig)  
               epsaee(mu,isig,4,k) =
     1              epsaa(mu)*clr(4,1,isig)+epsaz(mu)*clr(4,2,isig)  
               epszee(mu,isig,3,k) =
     1              epsza(mu)*clr(3,1,isig)+epszz(mu)*clr(3,2,isig)
               epszee(mu,isig,4,k) =
     1              epsza(mu)*clr(4,1,isig)+epszz(mu)*clr(4,2,isig)
            enddo

	! VV-> mu+mu-:

        if (k.eq.1) then
          call contract_T2j(aauu(0,0,k,L),jqq(0,isig,kk),epsaa) !for AA 
          call contract_T2j(azuu(0,0,k,L),jqq(0,isig,kk),epsaz) !for AZ 
          call contract_T2j(zauu(0,0,k,L),jqq(0,isig,kk),epsza) !for ZA 
          call contract_T2j(zzuu(0,0,k,L),jqq(0,isig,kk),epszz) !for ZZ 
        else
          call contract_T1j(aauu(0,0,k,L),jqq(0,isig,kk),epsaa) !for AA 
          call contract_T1j(azuu(0,0,k,L),jqq(0,isig,kk),epsza) !for AZ 
          call contract_T1j(zauu(0,0,k,L),jqq(0,isig,kk),epsaz) !for ZA 
          call contract_T1j(zzuu(0,0,k,L),jqq(0,isig,kk),epszz) !for ZZ 
        endif

            do mu = 0,3 
               epsauu(mu,isig,3,k) =
     1              epsaa(mu)*clr(3,1,isig)+epsaz(mu)*clr(3,2,isig)  
               epsauu(mu,isig,4,k) =
     1              epsaa(mu)*clr(4,1,isig)+epsaz(mu)*clr(4,2,isig)  
               epszuu(mu,isig,3,k) =
     1              epsza(mu)*clr(3,1,isig)+epszz(mu)*clr(3,2,isig)
               epszuu(mu,isig,4,k) =
     1              epsza(mu)*clr(4,1,isig)+epszz(mu)*clr(4,2,isig)
            enddo
            do mu = 4,5           ! add momentum info to the currents
               epsaee(mu,isig,3,k) = zuu(mu)
               epsaee(mu,isig,4,k) = zuu(mu)
               epszee(mu,isig,3,k) = zuu(mu)
               epszee(mu,isig,4,k) = zuu(mu)
               
	       epsauu(mu,isig,3,k) = zee(mu)
               epsauu(mu,isig,4,k) = zee(mu)
               epszuu(mu,isig,3,k) = zee(mu)
               epszuu(mu,isig,4,k) = zee(mu)
            enddo	
	
	enddo !isig

c and same for the CC processes (V attached to j43 or j21 current):

         isig = -1
	 if (k.eq.1) then
         	call contract_T2j(CCee(0,0,k,L),jqq(0,isig,kk),   epsee) 
         	call contract_T2j(CCuu(0,0,k,L),jqq(0,isig,kk),   epsuu)
         	call contract_T1j(CCee6(0,0,k,L),jqq(0,isig,kk),   epsee6) 
         	call contract_T1j(CCuu6(0,0,k,L),jqq(0,isig,kk),   epsuu6)
         else
	 	call contract_T1j(CCee(0,0,k,L),jqq(0,isig,kk),   epsee) 
         	call contract_T1j(CCuu(0,0,k,L),jqq(0,isig,kk),   epsuu)
	 	call contract_T2j(CCee6(0,0,k,L),jqq(0,isig,kk),   epsee6) 
         	call contract_T2j(CCuu6(0,0,k,L),jqq(0,isig,kk),   epsuu6)
	 endif
	 	 
         do mu = 0,3 
            epsCCee(mu,isig,k) =
     1           epsee(mu)*clr(3,3,isig) ! coupling to fermion
            epsCCuu(mu,isig,k) =
     1           epsuu(mu)*clr(3,3,isig) 
c
            epsCCee6(mu,isig,k) =
     1           epsee6(mu)*clr(3,3,isig) ! coupling to fermion
            epsCCuu6(mu,isig,k) =
     1           epsuu6(mu)*clr(3,3,isig) 
         enddo
	 
c extra terms for (q_mu q_nu)-contributions of massive boson propagators:	 
         if (k.eq.1) then 
            qepsCCee(k) = -dotrc(quu,epsee)*zm2i(3)
            qepsCCuu(k) = -dotrc(qee,epsuu)*zm2i(3)
            qepsCCee6(k) = -dotrc(quu,epsee6)*zm2i(3)
            qepsCCuu6(k) = -dotrc(qee,epsuu6)*zm2i(3)
	 else
            qepsCCee(k) = dotrc(qee,epsee)*zm2i(3)
            qepsCCuu(k) = dotrc(quu,epsuu)*zm2i(3)	        
            qepsCCee6(k) = dotrc(qee,epsee6)*zm2i(3)
            qepsCCuu6(k) = dotrc(quu,epsuu6)*zm2i(3)	        
         endif

         do mu = 4,5            ! add momentum info to the currents
            epsCCee(mu,isig,k) = zuu(mu)
            epsCCuu(mu,isig,k) = zee(mu)
c
            epsCCee6(mu,isig,k) = zee(mu)
            epsCCuu6(mu,isig,k) = zuu(mu)
         enddo 
      enddo !k
            	
      do isig = -1,1,2 	
            
      jj21au(isig) = dotcc(au,jqq(0,isig,1))
      jj21zu(isig) = dotcc(zu,jqq(0,isig,1))
      jj21ae(isig) = dotcc(ae,jqq(0,isig,1))
      jj21ze(isig) = dotcc(ze,jqq(0,isig,1))

      jj43au(isig) = dotcc(au,jqq(0,isig,2))
      jj43zu(isig) = dotcc(zu,jqq(0,isig,2))
      jj43ae(isig) = dotcc(ae,jqq(0,isig,2))
      jj43ze(isig) = dotcc(ze,jqq(0,isig,2))

      enddo !isig
c
c --------------
c
c now construct the contribution to the amplitude by current contraction 
c virtual contributions are assembled in subroutine boxline
                ! 2 bosons attached to 12 line     
c
c	NC first:
c
	do isig1 = -1,1,2
	   do isig3 = -1,1,2
	   
c  gauge term for Z propagator:
	  do k=3,4 !up/down type
	      qepszee(k) = -dotrc(quu,epszee(0,isig3,k,1))*zm2i(2) 
	      qepszuu(k) = -dotrc(qee,epszuu(0,isig3,k,1))*zm2i(2)     
	   enddo !k
	   
	   do k = 1,4 ! e+e- from upper line:
	   
	   m1aae = dotcc(jea(0,isig1,2),epsauu(0,isig3,ifl(3,k),1))
	   m1aze = dotcc(jea(0,isig1,2),epszuu(0,isig3,ifl(3,k),1))
	   m1zae = dotcc(jez(0,isig1,2),epsauu(0,isig3,ifl(3,k),1))
	   m1zze = dotcc(jez(0,isig1,2),epszuu(0,isig3,ifl(3,k),1))
	   
	   z1aze = qepszuu(ifl(3,k))*jj21ae(isig1)    
	   z1zze = qepszuu(ifl(3,k))*jj21ze(isig1)    
	   
	   mz1aae =  m1aae*clr(ifl(1,k),1,isig1)*clr(ifl(2,k),1,isig1)
	   mz1aze = (m1aze+z1aze)*clr(ifl(1,k),2,isig1)*clr(ifl(2,k),1,isig1)
	   mz1zae =  m1zae*clr(ifl(1,k),1,isig1)*clr(ifl(2,k),2,isig1)
	   mz1zze = (m1zze+z1zze)*clr(ifl(1,k),2,isig1)*clr(ifl(2,k),2,isig1)
	   	   
	   m1e = mz1aae+mz1aze+mz1zae+mz1zze
	   	   	   
	   m2aae = dotcc(jea(0,isig1,1),epsauu(0,isig3,ifl(3,k),1))
	   m2aze = dotcc(jea(0,isig1,1),epszuu(0,isig3,ifl(3,k),1))
	   m2zae = dotcc(jez(0,isig1,1),epsauu(0,isig3,ifl(3,k),1))
	   m2zze = dotcc(jez(0,isig1,1),epszuu(0,isig3,ifl(3,k),1))
	   	   
	   mz2aae =  m2aae*clr(ifl(2,k),1,isig1)*clr(ifl(1,k),1,isig1)
	   mz2aze = (m2aze-z1aze)*clr(ifl(2,k),2,isig1)*clr(ifl(1,k),1,isig1)
	   mz2zae =  m2zae*clr(ifl(2,k),1,isig1)*clr(ifl(1,k),2,isig1)
	   mz2zze = (m2zze-z1zze)*clr(ifl(2,k),2,isig1)*clr(ifl(1,k),2,isig1)
	   	   
	   m2e = mz2aae+mz2aze+mz2zae+mz2zze
	   
	   
	   ! mu+mu- from upper line:

	   m1aau = dotcc(jua(0,isig1,2),epsaee(0,isig3,ifl(3,k),1))
	   m1azu = dotcc(jua(0,isig1,2),epszee(0,isig3,ifl(3,k),1))
	   m1zau = dotcc(juz(0,isig1,2),epsaee(0,isig3,ifl(3,k),1))
	   m1zzu = dotcc(juz(0,isig1,2),epszee(0,isig3,ifl(3,k),1))
	   
	   z1azu = qepszee(ifl(3,k))*jj21au(isig1)  
	   z1zzu = qepszee(ifl(3,k))*jj21zu(isig1)  
	   
	   mz1aau =  m1aau*clr(ifl(1,k),1,isig1)*clr(ifl(2,k),1,isig1)
	   mz1azu = (m1azu+z1azu)*clr(ifl(1,k),2,isig1)*clr(ifl(2,k),1,isig1)
	   mz1zau =  m1zau*clr(ifl(1,k),1,isig1)*clr(ifl(2,k),2,isig1)
	   mz1zzu = (m1zzu+z1zzu)*clr(ifl(1,k),2,isig1)*clr(ifl(2,k),2,isig1)
	   	   
	   m1u = mz1aau+mz1azu+mz1zau+mz1zzu
	   	   
	   m2aau = dotcc(jua(0,isig1,1),epsaee(0,isig3,ifl(3,k),1))
	   m2azu = dotcc(jua(0,isig1,1),epszee(0,isig3,ifl(3,k),1))
	   m2zau = dotcc(juz(0,isig1,1),epsaee(0,isig3,ifl(3,k),1))
	   m2zzu = dotcc(juz(0,isig1,1),epszee(0,isig3,ifl(3,k),1))
	   	   
	   mz2aau =  m2aau*clr(ifl(2,k),1,isig1)*clr(ifl(1,k),1,isig1)
	   mz2azu = (m2azu-z1azu)*clr(ifl(2,k),2,isig1)*clr(ifl(1,k),1,isig1)
	   mz2zau =  m2zau*clr(ifl(2,k),1,isig1)*clr(ifl(1,k),2,isig1)
	   mz2zzu = (m2zzu-z1zzu)*clr(ifl(2,k),2,isig1)*clr(ifl(1,k),2,isig1)
	   	   
	   m2u = mz2aau+mz2azu+mz2zau+mz2zzu
	   
	   mat(k,isig1,isig3,2) = m1e+m1u+m2e+m2u
	   
	   
	   if (lbox) then
	    ! e+e- from upper line:
c this is the first permutation 	   
	   mv1aae = dotcc(jvea(0,isig1,is1,2),epsauu(0,isig3,ifl(3,k),1))
     #		   *clr(ifl(1,k),1,isig1)*clr(ifl(2,k),1,isig1)
	   mv1aze = dotcc(jvea(0,isig1,is1,2),epszuu(0,isig3,ifl(3,k),1))
     #		   *clr(ifl(1,k),2,isig1)*clr(ifl(2,k),1,isig1)
	   mv1zae = dotcc(jvez(0,isig1,is1,2),epsauu(0,isig3,ifl(3,k),1))
     #		   *clr(ifl(1,k),1,isig1)*clr(ifl(2,k),2,isig1)
	   mv1zze = dotcc(jvez(0,isig1,is1,2),epszuu(0,isig3,ifl(3,k),1))
     #		   *clr(ifl(1,k),2,isig1)*clr(ifl(2,k),2,isig1)
	  	   	   
	   mv1e = mv1aae+mv1aze+mv1zae+mv1zze

c this is the second permutation
	   	   	   
	   mv2aae = dotcc(jvea(0,isig1,is1,1),epsauu(0,isig3,ifl(3,k),1))
     #		   *clr(ifl(2,k),1,isig1)*clr(ifl(1,k),1,isig1)
	   mv2aze = dotcc(jvea(0,isig1,is1,1),epszuu(0,isig3,ifl(3,k),1))
     #		   *clr(ifl(2,k),2,isig1)*clr(ifl(1,k),1,isig1)
	   mv2zae = dotcc(jvez(0,isig1,is1,1),epsauu(0,isig3,ifl(3,k),1))
     #		   *clr(ifl(2,k),1,isig1)*clr(ifl(1,k),2,isig1)
	   mv2zze = dotcc(jvez(0,isig1,is1,1),epszuu(0,isig3,ifl(3,k),1))
     #		   *clr(ifl(2,k),2,isig1)*clr(ifl(1,k),2,isig1)
	   	   	   	   
	   mv2e = mv2aae+mv2aze+mv2zae+mv2zze
	   
	   
	   ! mu+mu- from upper line:
	   
	   mv1aau = dotcc(jvua(0,isig1,is1,2),epsaee(0,isig3,ifl(3,k),1))
     #		   *clr(ifl(1,k),1,isig1)*clr(ifl(2,k),1,isig1)
	   mv1azu = dotcc(jvua(0,isig1,is1,2),epszee(0,isig3,ifl(3,k),1))
     #		   *clr(ifl(1,k),2,isig1)*clr(ifl(2,k),1,isig1)
	   mv1zau = dotcc(jvuz(0,isig1,is1,2),epsaee(0,isig3,ifl(3,k),1))
     #		   *clr(ifl(1,k),1,isig1)*clr(ifl(2,k),2,isig1)
	   mv1zzu = dotcc(jvuz(0,isig1,is1,2),epszee(0,isig3,ifl(3,k),1))
     #		   *clr(ifl(1,k),2,isig1)*clr(ifl(2,k),2,isig1)
	   	   	   
	   mv1u = mv1aau+mv1azu+mv1zau+mv1zzu
	   	   
	   mv2aau = dotcc(jvua(0,isig1,is1,1),epsaee(0,isig3,ifl(3,k),1))
     #		   *clr(ifl(2,k),1,isig1)*clr(ifl(1,k),1,isig1)
	   mv2azu = dotcc(jvua(0,isig1,is1,1),epszee(0,isig3,ifl(3,k),1))
     #		   *clr(ifl(2,k),2,isig1)*clr(ifl(1,k),1,isig1)
	   mv2zau = dotcc(jvuz(0,isig1,is1,1),epsaee(0,isig3,ifl(3,k),1))
     #		   *clr(ifl(2,k),1,isig1)*clr(ifl(1,k),2,isig1)
	   mv2zzu = dotcc(jvuz(0,isig1,is1,1),epszee(0,isig3,ifl(3,k),1))
     #		   *clr(ifl(2,k),2,isig1)*clr(ifl(1,k),2,isig1)
	   	   	   	   
	   mv2u = mv2aau+mv2azu+mv2zau+mv2zzu
	   
	   matv(k,isig1,isig3,2) = mv1e+mv1u+mv2e+mv2u
	   
	   
	   endif !lbox
	   
	   enddo !k
	   
	   enddo !isig3
	enddo ! isig1  
c
c ------------------
c   
c	CC next:
      
      isig = -1
      k = 5                     !udsc, e+e- from upper line
      m1ae = dotcc(jea(0,-1,2),epsCCuu(0,isig,1))*clr(ifl(2,k),1,-1)   
      m1ze = dotcc(jez(0,-1,2),epsCCuu(0,isig,1))*clr(ifl(2,k),2,-1)   
      m2ae = dotcc(jea(0,-1,1),epsCCuu(0,isig,1))*clr(ifl(1,k),1,-1)   
      m2ze = dotcc(jez(0,-1,1),epsCCuu(0,isig,1))*clr(ifl(1,k),2,-1)   

      z1ae = qepsCCuu(1)*jj21ae(isig)*clr(ifl(1,k),3,isig)*clr(ifl(2,k),1,isig)
      z1ze = qepsCCuu(1)*jj21ze(isig)*clr(ifl(1,k),3,isig)*clr(ifl(2,k),2,isig)
      z2ze = qepsCCuu(1)*jj21ze(isig)*clr(ifl(2,k),3,isig)*clr(ifl(1,k),2,isig)
      z2ae = qepsCCuu(1)*jj21ae(isig)*clr(ifl(2,k),3,isig)*clr(ifl(1,k),1,isig)

     	 ! Achtung: A an u/d Kopplung unterschiedlich!! 
     
	mate12 = clr(3,3,-1)*(m2ae-z2ae+m2ze-z2ze)
     #	      +  clr(3,3,-1)*(m1ae+z1ae+m1ze+z1ze)


     	                       !udsc, mu+mu- from upper line
      m1au = dotcc(jua(0,-1,2),epsCCee(0,isig,1))*clr(ifl(2,k),1,-1)   
      m1zu = dotcc(juz(0,-1,2),epsCCee(0,isig,1))*clr(ifl(2,k),2,-1)   
      m2au = dotcc(jua(0,-1,1),epsCCee(0,isig,1))*clr(ifl(1,k),1,-1)   
      m2zu = dotcc(juz(0,-1,1),epsCCee(0,isig,1))*clr(ifl(1,k),2,-1)   

      z1au = qepsCCee(1)*jj21au(isig)*clr(ifl(1,k),3,isig)*clr(ifl(2,k),1,isig)
      z1zu = qepsCCee(1)*jj21zu(isig)*clr(ifl(1,k),3,isig)*clr(ifl(2,k),2,isig)
      z2zu = qepsCCee(1)*jj21zu(isig)*clr(ifl(2,k),3,isig)*clr(ifl(1,k),2,isig)
      z2au = qepsCCee(1)*jj21au(isig)*clr(ifl(2,k),3,isig)*clr(ifl(1,k),1,isig)


	matu12 = clr(3,3,-1)*(m2au-z2au+m2zu-z2zu)
     #	      +  clr(3,3,-1)*(m1au+z1au+m1zu+z1zu)
     
       mat(k,-1,isig,2) = mate12 + matu12
c
c-------------
c
      isig = -1
      k = 6                     !ducs, e+e- from upper line
      m1ae = dotcc(jea(0,-1,2),epsCCuu6(0,isig,1))*clr(ifl(2,k),1,-1)	    
      m1ze = dotcc(jez(0,-1,2),epsCCuu6(0,isig,1))*clr(ifl(2,k),2,-1)	      
      m2ae = dotcc(jea(0,-1,1),epsCCuu6(0,isig,1))*clr(ifl(1,k),1,-1)	    
      m2ze = dotcc(jez(0,-1,1),epsCCuu6(0,isig,1))*clr(ifl(1,k),2,-1)	    

      z1ae = qepsCCuu6(1)*jj21ae(isig)*clr(ifl(1,k),3,isig)*clr(ifl(2,k),1,isig)
      z1ze = qepsCCuu6(1)*jj21ze(isig)*clr(ifl(1,k),3,isig)*clr(ifl(2,k),2,isig)
      z2ze = qepsCCuu6(1)*jj21ze(isig)*clr(ifl(2,k),3,isig)*clr(ifl(1,k),2,isig)
      z2ae = qepsCCuu6(1)*jj21ae(isig)*clr(ifl(2,k),3,isig)*clr(ifl(1,k),1,isig)

	mate12 = clr(3,3,-1)*(m2ae-z2ae+m2ze-z2ze)
     #	      +  clr(3,3,-1)*(m1ae+z1ae+m1ze+z1ze)


     	                       !ducs, mu+mu- from upper line
      m1au = dotcc(jua(0,-1,2),epsCCee6(0,isig,1))*clr(ifl(2,k),1,-1)	   
      m1zu = dotcc(juz(0,-1,2),epsCCee6(0,isig,1))*clr(ifl(2,k),2,-1)	    
      m2au = dotcc(jua(0,-1,1),epsCCee6(0,isig,1))*clr(ifl(1,k),1,-1)	    
      m2zu = dotcc(juz(0,-1,1),epsCCee6(0,isig,1))*clr(ifl(1,k),2,-1)	    

      z1au = qepsCCee6(1)*jj21au(isig)*clr(ifl(1,k),3,isig)*clr(ifl(2,k),1,isig)
      z1zu = qepsCCee6(1)*jj21zu(isig)*clr(ifl(1,k),3,isig)*clr(ifl(2,k),2,isig)
      z2zu = qepsCCee6(1)*jj21zu(isig)*clr(ifl(2,k),3,isig)*clr(ifl(1,k),2,isig)
      z2au = qepsCCee6(1)*jj21au(isig)*clr(ifl(2,k),3,isig)*clr(ifl(1,k),1,isig)


	matu12 = clr(3,3,-1)*(m2au-z2au+m2zu-z2zu)
     #	      +  clr(3,3,-1)*(m1au+z1au+m1zu+z1zu)
     
       mat(k,-1,isig,2) = mate12 + matu12
c
c ------------------
c   
      if (lbox) then
      
      isig = -1
      k = 5                     !udsc, e+e- from upper line
      mv1ae = dotcc(jvea(0,-1,is1,2),epsCCuu(0,isig,1))*clr(ifl(2,k),1,-1)  
      mv1ze = dotcc(jvez(0,-1,is1,2),epsCCuu(0,isig,1))*clr(ifl(2,k),2,-1)  
      mv2ae = dotcc(jvea(0,-1,is1,1),epsCCuu(0,isig,1))*clr(ifl(1,k),1,-1)  
      mv2ze = dotcc(jvez(0,-1,is1,1),epsCCuu(0,isig,1))*clr(ifl(1,k),2,-1)  
     
	matve12 = clr(3,3,-1)*(mv2ae+mv2ze+mv1ae+mv1ze)


     	                       !udsc, mu+mu- from upper line
      mv1au = dotcc(jvua(0,-1,is1,2),epsCCee(0,isig,1))*clr(ifl(2,k),1,-1)  
      mv1zu = dotcc(jvuz(0,-1,is1,2),epsCCee(0,isig,1))*clr(ifl(2,k),2,-1)  
      mv2au = dotcc(jvua(0,-1,is1,1),epsCCee(0,isig,1))*clr(ifl(1,k),1,-1)  
      mv2zu = dotcc(jvuz(0,-1,is1,1),epsCCee(0,isig,1))*clr(ifl(1,k),2,-1)  

	matvu12 = clr(3,3,-1)*(mv2au+mv2zu+mv1au+mv1zu)
     
       matv(k,-1,isig,2) = matve12 + matvu12

c    
      isig = -1
      k = 6                     !udsc, e+e- from upper line
      mv1ae = dotcc(jvea(0,-1,is1,2),epsCCuu6(0,isig,1))*clr(ifl(2,k),1,-1) 
      mv1ze = dotcc(jvez(0,-1,is1,2),epsCCuu6(0,isig,1))*clr(ifl(2,k),2,-1) 
      mv2ae = dotcc(jvea(0,-1,is1,1),epsCCuu6(0,isig,1))*clr(ifl(1,k),1,-1) 
      mv2ze = dotcc(jvez(0,-1,is1,1),epsCCuu6(0,isig,1))*clr(ifl(1,k),2,-1) 
      
	matve12 = clr(3,3,-1)*(mv2ae+mv2ze+mv1ae+mv1ze)


     	                       !udsc, mu+mu- from upper line
      mv1au = dotcc(jvua(0,-1,is1,2),epsCCee6(0,isig,1))*clr(ifl(2,k),1,-1) 
      mv1zu = dotcc(jvuz(0,-1,is1,2),epsCCee6(0,isig,1))*clr(ifl(2,k),2,-1) 
      mv2au = dotcc(jvua(0,-1,is1,1),epsCCee6(0,isig,1))*clr(ifl(1,k),1,-1) 
      mv2zu = dotcc(jvuz(0,-1,is1,1),epsCCee6(0,isig,1))*clr(ifl(1,k),2,-1) 

	matvu12 = clr(3,3,-1)*(mv2au+mv2zu+mv1au+mv1zu)
     
       matv(k,-1,isig,2) = matve12 + matvu12


      endif !lbox

c ----------------------------------
      
c repeat the same for 2 bosons attached to 34 line:
c
cc	NC first:
c
	do isig1 = -1,1,2
	   do isig3 = -1,1,2
	   
c  gauge term for Z propagator:
	  do k=3,4 !up/down type
	      qepszee(k) = dotrc(qee,epszee(0,isig1,k,2))*zm2i(2) 
	      qepszuu(k) = dotrc(quu,epszuu(0,isig1,k,2))*zm2i(2) 	    
	   enddo !k
	   
	   do k = 1,4 ! e+e- from lower line:
	   
	   m3aae = dotcc(jea(0,isig3,4),epsauu(0,isig1,ifl(1,k),2))
	   m3aze = dotcc(jea(0,isig3,4),epszuu(0,isig1,ifl(1,k),2))
	   m3zae = dotcc(jez(0,isig3,4),epsauu(0,isig1,ifl(1,k),2))
	   m3zze = dotcc(jez(0,isig3,4),epszuu(0,isig1,ifl(1,k),2))
	   
	   z3aze = qepszuu(ifl(1,k))*jj43ae(isig3)    
	   z3zze = qepszuu(ifl(1,k))*jj43ze(isig3)    
	   
	   mz3aae =  m3aae*clr(ifl(3,k),1,isig3)*clr(ifl(4,k),1,isig3)
	   mz3aze = (m3aze+z3aze)*clr(ifl(3,k),2,isig3)*clr(ifl(4,k),1,isig3)
	   mz3zae =  m3zae*clr(ifl(3,k),1,isig3)*clr(ifl(4,k),2,isig3)
	   mz3zze = (m3zze+z3zze)*clr(ifl(3,k),2,isig3)*clr(ifl(4,k),2,isig3)
	   	   
	   m3e = mz3aae+mz3aze+mz3zae+mz3zze
	   	   	   
	   m4aae = dotcc(jea(0,isig3,3),epsauu(0,isig1,ifl(1,k),2))
	   m4aze = dotcc(jea(0,isig3,3),epszuu(0,isig1,ifl(1,k),2))
	   m4zae = dotcc(jez(0,isig3,3),epsauu(0,isig1,ifl(1,k),2))
	   m4zze = dotcc(jez(0,isig3,3),epszuu(0,isig1,ifl(1,k),2))
	   	   
	   mz4aae =  m4aae*clr(ifl(4,k),1,isig3)*clr(ifl(3,k),1,isig3)
	   mz4aze = (m4aze-z3aze)*clr(ifl(4,k),2,isig3)*clr(ifl(3,k),1,isig3)
	   mz4zae =  m4zae*clr(ifl(4,k),1,isig3)*clr(ifl(3,k),2,isig3)
	   mz4zze = (m4zze-z3zze)*clr(ifl(4,k),2,isig3)*clr(ifl(3,k),2,isig3)
	   	   
	   m4e = mz4aae+mz4aze+mz4zae+mz4zze
	   
	   
	   ! mu+mu- from lower line:
	   m3aau = dotcc(jua(0,isig3,4),epsaee(0,isig1,ifl(1,k),2))
	   m3azu = dotcc(jua(0,isig3,4),epszee(0,isig1,ifl(1,k),2))
	   m3zau = dotcc(juz(0,isig3,4),epsaee(0,isig1,ifl(1,k),2))
	   m3zzu = dotcc(juz(0,isig3,4),epszee(0,isig1,ifl(1,k),2))
	   
	   z3azu = qepszee(ifl(1,k))*jj43au(isig3)    
	   z3zzu = qepszee(ifl(1,k))*jj43zu(isig3)    
	   
	   mz3aau =  m3aau*clr(ifl(3,k),1,isig3)*clr(ifl(4,k),1,isig3)
	   mz3azu = (m3azu+z3azu)*clr(ifl(3,k),2,isig3)*clr(ifl(4,k),1,isig3)
	   mz3zau =  m3zau*clr(ifl(3,k),1,isig3)*clr(ifl(4,k),2,isig3)
	   mz3zzu = (m3zzu+z3zzu)*clr(ifl(3,k),2,isig3)*clr(ifl(4,k),2,isig3)
	   	   
	   m3u = mz3aau+mz3azu+mz3zau+mz3zzu
	   	   	   
	   m4aau = dotcc(jua(0,isig3,3),epsaee(0,isig1,ifl(1,k),2))
	   m4azu = dotcc(jua(0,isig3,3),epszee(0,isig1,ifl(1,k),2))
	   m4zau = dotcc(juz(0,isig3,3),epsaee(0,isig1,ifl(1,k),2))
	   m4zzu = dotcc(juz(0,isig3,3),epszee(0,isig1,ifl(1,k),2))
	   	   
	   mz4aau =  m4aau*clr(ifl(4,k),1,isig3)*clr(ifl(3,k),1,isig3)
	   mz4azu = (m4azu-z3azu)*clr(ifl(4,k),2,isig3)*clr(ifl(3,k),1,isig3)
	   mz4zau =  m4zau*clr(ifl(4,k),1,isig3)*clr(ifl(3,k),2,isig3)
	   mz4zzu = (m4zzu-z3zzu)*clr(ifl(4,k),2,isig3)*clr(ifl(3,k),2,isig3)
	   	   
	   m4u = mz4aau+mz4azu+mz4zau+mz4zzu
	  	   
	   
	   mat(k,isig1,isig3,3) = m3e+m3u+m4e+m4u
	   
	   
	   enddo !k
	   
	   enddo !isig3
	enddo ! isig1  
c
c -----------------------------------------
c
c  virtual corrections:  
c	(gauge terms for massive gauge boson propagators 
c		   not needed when contracted with boxes)
c
c
	if (lbox) then

	do isig1 = -1,1,2
	   do isig3 = -1,1,2
	   	   
	   do k = 1,4 ! e+e- from lower line:
	   
	   mv3aae = dotcc(jvea(0,isig3,is3,4),epsauu(0,isig1,ifl(1,k),2))
     #			   *clr(ifl(3,k),1,isig3)*clr(ifl(4,k),1,isig3)
	   mv3aze = dotcc(jvea(0,isig3,is3,4),epszuu(0,isig1,ifl(1,k),2))
     #			   *clr(ifl(3,k),2,isig3)*clr(ifl(4,k),1,isig3)
	   mv3zae = dotcc(jvez(0,isig3,is3,4),epsauu(0,isig1,ifl(1,k),2))
     #			   *clr(ifl(3,k),1,isig3)*clr(ifl(4,k),2,isig3)
	   mv3zze = dotcc(jvez(0,isig3,is3,4),epszuu(0,isig1,ifl(1,k),2))
     #			   *clr(ifl(3,k),2,isig3)*clr(ifl(4,k),2,isig3)
	   	   	   
	   mv3e = mv3aae+mv3aze+mv3zae+mv3zze
	   	   	   
	   mv4aae = dotcc(jvea(0,isig3,is3,3),epsauu(0,isig1,ifl(1,k),2))
     #			   *clr(ifl(4,k),1,isig3)*clr(ifl(3,k),1,isig3)
	   mv4aze = dotcc(jvea(0,isig3,is3,3),epszuu(0,isig1,ifl(1,k),2))
     #			   *clr(ifl(4,k),2,isig3)*clr(ifl(3,k),1,isig3)
	   mv4zae = dotcc(jvez(0,isig3,is3,3),epsauu(0,isig1,ifl(1,k),2))
     #			   *clr(ifl(4,k),1,isig3)*clr(ifl(3,k),2,isig3)
	   mv4zze = dotcc(jvez(0,isig3,is3,3),epszuu(0,isig1,ifl(1,k),2))
     #			   *clr(ifl(4,k),2,isig3)*clr(ifl(3,k),2,isig3)
	   	   
	   mv4e = mv4aae+mv4aze+mv4zae+mv4zze
	   
	   
	   ! mu+mu- from lower line:
	   mv3aau = dotcc(jvua(0,isig3,is3,4),epsaee(0,isig1,ifl(1,k),2))
     #			   *clr(ifl(3,k),1,isig3)*clr(ifl(4,k),1,isig3)
	   mv3azu = dotcc(jvua(0,isig3,is3,4),epszee(0,isig1,ifl(1,k),2))
     #			   *clr(ifl(3,k),2,isig3)*clr(ifl(4,k),1,isig3)
	   mv3zau = dotcc(jvuz(0,isig3,is3,4),epsaee(0,isig1,ifl(1,k),2))
     #			   *clr(ifl(3,k),1,isig3)*clr(ifl(4,k),2,isig3)
	   mv3zzu = dotcc(jvuz(0,isig3,is3,4),epszee(0,isig1,ifl(1,k),2))
     #			   *clr(ifl(3,k),2,isig3)*clr(ifl(4,k),2,isig3)
	   	   	   
	   mv3u = mv3aau+mv3azu+mv3zau+mv3zzu
	   	   	   
	   mv4aau = dotcc(jvua(0,isig3,is3,3),epsaee(0,isig1,ifl(1,k),2))
     #			   *clr(ifl(4,k),1,isig3)*clr(ifl(3,k),1,isig3)
	   mv4azu = dotcc(jvua(0,isig3,is3,3),epszee(0,isig1,ifl(1,k),2))
     #			   *clr(ifl(4,k),2,isig3)*clr(ifl(3,k),1,isig3)
	   mv4zau = dotcc(jvuz(0,isig3,is3,3),epsaee(0,isig1,ifl(1,k),2))
     #			   *clr(ifl(4,k),1,isig3)*clr(ifl(3,k),2,isig3)
	   mv4zzu = dotcc(jvuz(0,isig3,is3,3),epszee(0,isig1,ifl(1,k),2))
     #			   *clr(ifl(4,k),2,isig3)*clr(ifl(3,k),2,isig3)
	   	   	   	   
	   mv4u = mv4aau+mv4azu+mv4zau+mv4zzu
	  	   
	   
	   matv(k,isig1,isig3,3) = mv3e+mv3u+mv4e+mv4u
	   
	   
	   enddo !k
	   
	   enddo !isig3
	enddo ! isig1
	
	endif !lbox  
	

c -----------------------------------------
c
C CC:
		
      isig = -1
      k = 5                     !udsc, e+e- from lower line
      m3ae = dotcc(jea(0,-1,4),epsCCuu(0,isig,2))*clr(ifl(4,k),1,-1)	  
      m3ze = dotcc(jez(0,-1,4),epsCCuu(0,isig,2))*clr(ifl(4,k),2,-1)	    
      m4ae = dotcc(jea(0,-1,3),epsCCuu(0,isig,2))*clr(ifl(3,k),1,-1)	  
      m4ze = dotcc(jez(0,-1,3),epsCCuu(0,isig,2))*clr(ifl(3,k),2,-1)	  

      z3ae = qepsCCuu(2)*jj43ae(isig)*clr(ifl(3,k),3,isig)*clr(ifl(4,k),1,isig)
      z3ze = qepsCCuu(2)*jj43ze(isig)*clr(ifl(3,k),3,isig)*clr(ifl(4,k),2,isig)
      z4ze = qepsCCuu(2)*jj43ze(isig)*clr(ifl(4,k),3,isig)*clr(ifl(3,k),2,isig)
      z4ae = qepsCCuu(2)*jj43ae(isig)*clr(ifl(4,k),3,isig)*clr(ifl(3,k),1,isig)
     
	mate34 = clr(3,3,-1)*(m4ae-z4ae+m4ze-z4ze)
     #	      +  clr(3,3,-1)*(m3ae+z3ae+m3ze+z3ze)

	
	     	                 !udsc, mu+mu- from lower line
      m3au = dotcc(jua(0,-1,4),epsCCee(0,isig,2))*clr(ifl(4,k),1,-1)	  
      m3zu = dotcc(juz(0,-1,4),epsCCee(0,isig,2))*clr(ifl(4,k),2,-1)	  
      m4au = dotcc(jua(0,-1,3),epsCCee(0,isig,2))*clr(ifl(3,k),1,-1)	  
      m4zu = dotcc(juz(0,-1,3),epsCCee(0,isig,2))*clr(ifl(3,k),2,-1)	  

      z3au = qepsCCee(2)*jj43au(isig)*clr(ifl(3,k),3,isig)*clr(ifl(4,k),1,isig)
      z3zu = qepsCCee(2)*jj43zu(isig)*clr(ifl(3,k),3,isig)*clr(ifl(4,k),2,isig)
      z4zu = qepsCCee(2)*jj43zu(isig)*clr(ifl(4,k),3,isig)*clr(ifl(3,k),2,isig)
      z4au = qepsCCee(2)*jj43au(isig)*clr(ifl(4,k),3,isig)*clr(ifl(3,k),1,isig)


	matu34 = clr(3,3,-1)*(m4au-z4au+m4zu-z4zu)
     #	      +  clr(3,3,-1)*(m3au+z3au+m3zu+z3zu)
     
       mat(k,-1,isig,3) = mate34 + matu34


      isig = -1
      k = 6                     !udsc, e+e- from upper line
      m3ae = dotcc(jea(0,-1,4),epsCCuu6(0,isig,2))*clr(ifl(4,k),1,-1)  
      m3ze = dotcc(jez(0,-1,4),epsCCuu6(0,isig,2))*clr(ifl(4,k),2,-1)  
      m4ae = dotcc(jea(0,-1,3),epsCCuu6(0,isig,2))*clr(ifl(3,k),1,-1)  
      m4ze = dotcc(jez(0,-1,3),epsCCuu6(0,isig,2))*clr(ifl(3,k),2,-1)  

      z3ae = qepsCCuu6(2)*jj43ae(isig)*clr(ifl(3,k),3,isig)*clr(ifl(4,k),1,isig)
      z3ze = qepsCCuu6(2)*jj43ze(isig)*clr(ifl(3,k),3,isig)*clr(ifl(4,k),2,isig)
      z4ze = qepsCCuu6(2)*jj43ze(isig)*clr(ifl(4,k),3,isig)*clr(ifl(3,k),2,isig)
      z4ae = qepsCCuu6(2)*jj43ae(isig)*clr(ifl(4,k),3,isig)*clr(ifl(3,k),1,isig)

	mate34 = clr(3,3,-1)*(m4ae-z4ae+m4ze-z4ze)
     #	      +  clr(3,3,-1)*(m3ae+z3ae+m3ze+z3ze)


     	                       !udsc, mu+mu- from upper line
      m3au = dotcc(jua(0,-1,4),epsCCee6(0,isig,2))*clr(ifl(4,k),1,-1)	 
      m3zu = dotcc(juz(0,-1,4),epsCCee6(0,isig,2))*clr(ifl(4,k),2,-1)	 
      m4au = dotcc(jua(0,-1,3),epsCCee6(0,isig,2))*clr(ifl(3,k),1,-1)	 
      m4zu = dotcc(juz(0,-1,3),epsCCee6(0,isig,2))*clr(ifl(3,k),2,-1)	 

      z3au = qepsCCee6(2)*jj43au(isig)*clr(ifl(3,k),3,isig)*clr(ifl(4,k),1,isig)
      z3zu = qepsCCee6(2)*jj43zu(isig)*clr(ifl(3,k),3,isig)*clr(ifl(4,k),2,isig)
      z4zu = qepsCCee6(2)*jj43zu(isig)*clr(ifl(4,k),3,isig)*clr(ifl(3,k),2,isig)
      z4au = qepsCCee6(2)*jj43au(isig)*clr(ifl(4,k),3,isig)*clr(ifl(3,k),1,isig)


	matu34 = clr(3,3,-1)*(m4au-z4au+m4zu-z4zu)
     #	      +  clr(3,3,-1)*(m3au+z3au+m3zu+z3zu)
     
       mat(k,-1,isig,3) = mate34 + matu34
c
c ------------------
c
C virtual corrections:

      if (lbox) then
		
      isig = -1
      k = 5                     !udsc, e+e- from lower line
      mv3ae = dotcc(jvea(0,-1,is3,4),epsCCuu(0,isig,2))*clr(ifl(4,k),1,-1)	
      mv3ze = dotcc(jvez(0,-1,is3,4),epsCCuu(0,isig,2))*clr(ifl(4,k),2,-1)	  
      mv4ae = dotcc(jvea(0,-1,is3,3),epsCCuu(0,isig,2))*clr(ifl(3,k),1,-1)	
      mv4ze = dotcc(jvez(0,-1,is3,3),epsCCuu(0,isig,2))*clr(ifl(3,k),2,-1)	
     
	matve34 = clr(3,3,-1)*(mv4ae+mv4ze+mv3ae+mv3ze)

	
	     	                 !udsc, mu+mu- from lower line
      mv3au = dotcc(jvua(0,-1,is3,4),epsCCee(0,isig,2))*clr(ifl(4,k),1,-1)	  
      mv3zu = dotcc(jvuz(0,-1,is3,4),epsCCee(0,isig,2))*clr(ifl(4,k),2,-1)	  
      mv4au = dotcc(jvua(0,-1,is3,3),epsCCee(0,isig,2))*clr(ifl(3,k),1,-1)	  
      mv4zu = dotcc(jvuz(0,-1,is3,3),epsCCee(0,isig,2))*clr(ifl(3,k),2,-1)	  

	 matvu34 = clr(3,3,-1)*(mv4au+mv4zu+mv3au+mv3zu)
     
       matv(k,-1,isig,3) = matve34 + matvu34

        
      isig = -1
      k = 6                     !udsc, e+e- from upper line
      mv3ae = dotcc(jvea(0,-1,is3,4),epsCCuu6(0,isig,2))*clr(ifl(4,k),1,-1)    
      mv3ze = dotcc(jvez(0,-1,is3,4),epsCCuu6(0,isig,2))*clr(ifl(4,k),2,-1)    
      mv4ae = dotcc(jvea(0,-1,is3,3),epsCCuu6(0,isig,2))*clr(ifl(3,k),1,-1)    
      mv4ze = dotcc(jvez(0,-1,is3,3),epsCCuu6(0,isig,2))*clr(ifl(3,k),2,-1)    

	matve34 = clr(3,3,-1)*(mv4ae+mv4ze+mv3ae+mv3ze)


     	                       !udsc, mu+mu- from upper line
      mv3au = dotcc(jvua(0,-1,is3,4),epsCCee6(0,isig,2))*clr(ifl(4,k),1,-1)    
      mv3zu = dotcc(jvuz(0,-1,is3,4),epsCCee6(0,isig,2))*clr(ifl(4,k),2,-1)    
      mv4au = dotcc(jvua(0,-1,is3,3),epsCCee6(0,isig,2))*clr(ifl(3,k),1,-1)    
      mv4zu = dotcc(jvuz(0,-1,is3,3),epsCCee6(0,isig,2))*clr(ifl(3,k),2,-1)


	matvu34 = clr(3,3,-1)*(mv4au+mv4zu+mv3au+mv3zu)
     
       matv(k,-1,isig,3) = matve34 + matvu34
 
 
       endif !lbox
      
c -----------------------------------------------------------------------
c
c next come the A/Z-->ZZ currents attached to the quark lines. For the virtual 
c corrections the most effective structure is the contraction of two polarization 
c vectors with one fermion line. First build these effective polarization vectors
c from the currents zzztens(mu) and azz(mu)
c
c NOTE: the zz and azz currents are conserved. Hence there is no 
c       need to consider q^mu * q^nu/m_Z^2  terms in the Z boson propagator

       If (lbox) then
         if(.not. lzs(js1)) then
            print*, "STOP3 in qqzaqq.F"
         ENDIF
       endif

      if (lzs(js1) .or. L.ne.Lold) then

         i = 1
         do isig = -1,1,2
            call ket2c(psi(1,isig,i),.true.,p(0,i),isig,pza,azz,
     1                 psia(1,isig,is1,i),fq(0,i))
            call ket2c(psi(1,isig,i),.true.,p(0,i),isig,pza,zzztens,
     1                 psiz(1,isig,is1,i),fq(0,i))
            call bra2c(psi(1,isig,i+1),.true.,p(0,i+1),isig,pza,azz,
     1                 psia(1,isig,is1,i+1),fq(0,i+1))
            call bra2c(psi(1,isig,i+1),.true.,p(0,i+1),isig,pza,zzztens,
     1                 psiz(1,isig,is1,i+1),fq(0,i+1))
         enddo
	 
         call curr6(1,psi(1,-1,i+1),p(0,i+1),
     1              psia(1,-1,is1,i),fq(0,i), ja(0,-1,is1,i)   )
         call curr6(1,psia(1,-1,is1,i+1),fq(0,i+1),
     1              psi(1,-1,i),p(0,i), ja(0,-1,is1,i+1) )
         call curr6(1,psi(1,-1,i+1),p(0,i+1),
     1              psiz(1,-1,is1,i),fq(0,i), jz(0,-1,is1,i)   )
         call curr6(1,psiz(1,-1,is1,i+1),fq(0,i+1),
     1              psi(1,-1,i),p(0,i), jz(0,-1,is1,i+1) )
      
#ifdef WITH_NLO               
	 if (lbox) then
            do i = 1,2
cfc               call boxlinec(1,psi(1,-1,1),psi(1,-1,2),p(0,1),p(0,2),
cfc     1                       .true., azz,3-i,   ja(0,-1,is1,i), 
cfc     2                       jva(0,-1,is1,i) )
cfc               call boxlinec(1,psi(1,-1,1),psi(1,-1,2),p(0,1),p(0,2),
cfc     1                       .false., zzztens,3-i,  jz(0,-1,is1,i), 
cfc     2                       jvz(0,-1,is1,i) )
         do mu = 0,3
            pvirtcalc(mu,1) = p(mu,1)
            pvirtcalc(mu,2) = -1d0 * p(mu,2)
            pvirtcalc(mu,3) = -1d0 * (qe(mu)+qa(mu))
            pvirtcalc(mu,4) = 1d0 * (p(mu,3) - p(mu,4))
            temp=pvirtcalc(mu,1)+pvirtcalc(mu,2)+pvirtcalc(mu,3)+pvirtcalc(mu,4)
c            print*,"Mom. conserv.", temp
         enddo
           scale = scaleupperline
          call BOX_VBF_CURR(!
     $      pvirtcalc(0,1), pvirtcalc(0,3), ! momenta 
     $      pvirtcalc(0,4), pvirtcalc(0,2),  ! momenta
     $      psi(1,-1,2),psi(1,-1,1),! spinors
     $      azz, jva1(0,-1,is1,i),ja1(0,-1,is1,i), ! currents
     $      1, ! Max helicity 
     $      1, 3-i, ! comp, number of current
     $      scale,divmax,   ! Max helicity(1 spinor 2 ext. boson)
     $      mboxgauge(1,-1,0,0,1), ! virtuals gauge
     $      mboxgaugeborn(1,-1,0,0,1),!born gauge
     $      lwardtest(wardIndex)) ! Gauge test

         wardindex=wardindex+1
          if(lwardtest(wardIndex-1)) then
          call BOX_VBF_CURR(!
     $      pvirtcalc(0,1), pvirtcalc(0,3), ! momenta 
     $      pvirtcalc(0,4), pvirtcalc(0,2),  ! momenta
     $      psi(1,-1,2),psi(1,-1,1),! spinors
     $      zzztens, jvz1(0,-1,is1,i),jz1(0,-1,is1,i), ! currents
     $      1, ! Max helicity 
cfc NOTE: COMP=-1
     $      -1, 3-i, ! comp, number of current
cfc ************************************** 
     $      scale,divmax,   ! Max helicity(1 spinor 2 ext. boson)
     $      mboxgauge(1,-1,0,0,2), ! virtuals gauge
     $      mboxgaugeborn(1,-1,0,0,2),!born gauge
     $      lwardtest(wardIndex)) ! Gauge test

          wardindex=wardindex+1
          else
          wardindex=wardindex+1
          jvz1 = 0d0
          jz1 = 0d0
          endif

            enddo
         endif
#endif
	 
      endif

      if(lbox) then
         if(.not. lzs(js3)) then
            print*, "STOP4 in qqzaqq.F"
         ENDIF
       endif
      
      if (lzs(js3) .or. L.ne.Lold) then

         i = 3
         do isig = -1,1,2
            call ket2c(psi(1,isig,i),.true.,p(0,i),isig,pza,azz,
     1                 psia(1,isig,is3,i),fq(0,i))
            call ket2c(psi(1,isig,i),.true.,p(0,i),isig,pza,zzztens,
     1                 psiz(1,isig,is3,i),fq(0,i))
            call bra2c(psi(1,isig,i+1),.true.,p(0,i+1),isig,pza,azz,
     1                 psia(1,isig,is3,i+1),fq(0,i+1))
            call bra2c(psi(1,isig,i+1),.true.,p(0,i+1),isig,pza,zzztens,
     1                 psiz(1,isig,is3,i+1),fq(0,i+1))
         enddo
         call curr6(1,psi(1,-1,i+1),p(0,i+1),
     1              psia(1,-1,is3,i),fq(0,i), ja(0,-1,is3,i)   )
         call curr6(1,psia(1,-1,is3,i+1),fq(0,i+1),
     1              psi(1,-1,i),p(0,i), ja(0,-1,is3,i+1) )
         call curr6(1,psi(1,-1,i+1),p(0,i+1),
     1              psiz(1,-1,is3,i),fq(0,i), jz(0,-1,is3,i)   )
         call curr6(1,psiz(1,-1,is3,i+1),fq(0,i+1),
     1              psi(1,-1,i),p(0,i), jz(0,-1,is3,i+1) )

#ifdef WITH_NLO 
         if (lbox) then
            do i = 3,4
cfc               call boxlinec(1,psi(1,-1,3),psi(1,-1,4),p(0,3),p(0,4),
cfc     1                       .true., azz,5-i,   ja(0,-1,is3,i), 
cfc     2                       jva(0,-1,is3,i) )
cfc               call boxlinec(1,psi(1,-1,3),psi(1,-1,4),p(0,3),p(0,4),
cfc     1                       .false., zzztens,5-i,   
cfc     2                       jz(0,-1,is3,i), jvz(0,-1,is3,i) )


         do mu = 0,3
            pvirtcalc(mu,1) = p(mu,3)
            pvirtcalc(mu,2) = -1d0 * p(mu,4)
            pvirtcalc(mu,3) = -1d0 *( qe(mu)+ qa(mu))
            pvirtcalc(mu,4) = 1d0 * (p(mu,1) - p(mu,2))
            temp=pvirtcalc(mu,1)+pvirtcalc(mu,2)+pvirtcalc(mu,3)+pvirtcalc(mu,4)
c            print*,"Mom. conserv.", temp
         enddo
           scale = scalelowerline
          call BOX_VBF_CURR(!
     $      pvirtcalc(0,1), pvirtcalc(0,3), ! momenta 
     $      pvirtcalc(0,4), pvirtcalc(0,2),  ! momenta
     $      psi(1,-1,4),psi(1,-1,3),! spinors
     $      azz, jva1(0,-1,is3,i),ja1(0,-1,is3,i), ! currents
     $      1, ! Max helicity 
     $      1, 5-i, ! comp, number of current
     $      scale,divmax,   ! Max helicity(1 spinor 2 ext. boson)
     $      mboxgauge(1,-1,0,0,4), ! virtuals gauge
     $      mboxgaugeborn(1,-1,0,0,4),!born gauge
     $      lwardtest(wardIndex)) ! Gauge test

          wardindex=wardindex+1

          if(lwardtest(wardIndex-1)) then
          call BOX_VBF_CURR(!
     $      pvirtcalc(0,1), pvirtcalc(0,3), ! momenta 
     $      pvirtcalc(0,4), pvirtcalc(0,2),  ! momenta
     $      psi(1,-1,4),psi(1,-1,3),! spinors
     $      zzztens, jvz1(0,-1,is3,i),jz1(0,-1,is3,i), ! currents
     $      1, ! Max helicity 
cfc NOTE: COMP=-1
     $      -1, 5-i, ! comp, number of current
cfc ************************************** 
     $      scale,divmax,   ! Max helicity(1 spinor 2 ext. boson)
     $      mboxgauge(1,-1,0,0,5), ! virtuals gauge
     $      mboxgaugeborn(1,-1,0,0,5),!born gauge
     $      lwardtest(wardIndex)) ! Gauge test

          wardindex=wardindex+1
          else
          wardindex=wardindex+1
          jvz1 =0d0
          jz1=0d0
          endif

            enddo

      jvz = jvz1
      jva = jva1

      IF(debugza) then
      do i=1,4
         if(i.lt.3) then
            isig3=is1
          else
             isig3=is3
         endif
c         print*,"k",i
         do isig1=-1,1,2
c            print*, "isig1",isig1
            do mu=0,5
               if(lbox) then
               jz(mu,isig1,isig3,i)  =  jz1(mu,isig1,isig3,i)
               ja(mu,isig1,isig3,i)  = ja1(mu,isig1,isig3,i)
              endif

            if((abs(jvz(mu,isig1,isig3,i)/jvz1(mu,isig1,isig3,i)-1).gt.1d-3).and.(abs(jvz(mu,isig1,isig3,i)).ne.0d0)) then
            print*, " ZZZtens -> ll new result",mu,isig1,isig3,i
            print*,"jvz2",jvz(mu,isig1,isig3,i)/jvz1(mu,isig1,isig3,i)-1
            print*," jz2", jz(mu,isig1,isig3,i)/jz1(mu,isig1,isig3,i)-1
            print*
            print*
            endif
            if((abs(jva(mu,isig1,isig3,i)/jva1(mu,isig1,isig3,i)-1).gt.1d-3).and.(abs(jvea(mu,isig1,isig3,i)).ne.0d0)) then
            print*, " AAZtens -> ll new result",mu,isig1,isig3,i
            print*,"jva2",jva(mu,isig1,isig3,i)/jva1(mu,isig1,isig3,i)-1
            print*," ja2", ja(mu,isig1,isig3,i)/ja1(mu,isig1,isig3,i)-1
            print*
            print*
            endif

         enddo
        enddo
      enddo
      endif
cfc      stop

         endif
#endif
      endif
      

      do isig1 = -1,1,2
         do isig3 = -1,1,2

c---------------------------

c box correction to upper line: polarization vectors are 
c    jqq(mu,isig3,2)=j43 with momentum    p43        and
c    azz/zzztens(mu)         with momentum    pza

            ma(1) =      dotcc(ja(0,isig1,is1,2),jqq(0,isig3,2))
            ma(2) =      dotcc(ja(0,isig1,is1,1),jqq(0,isig3,2))
            mz(1) =      dotcc(jz(0,isig1,is1,2),jqq(0,isig3,2))
            mz(2) =      dotcc(jz(0,isig1,is1,1),jqq(0,isig3,2))
	    	    
            if (lbox) then
               mva(1) = dotcc(jva(0,isig1,is1,2),jqq(0,isig3,2))
               mva(2) = dotcc(jva(0,isig1,is1,1),jqq(0,isig3,2))
               mvz(1) = dotcc(jvz(0,isig1,is1,2),jqq(0,isig3,2))
               mvz(2) = dotcc(jvz(0,isig1,is1,1),jqq(0,isig3,2))
            endif
	    
	    
            do k = 1,4
               propt1(isig1,isig3,k,2) = 
     1            clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3)*prop43(1) 
     2          + clr(ifl(1,k),2,isig1)*clr(ifl(3,k),2,isig3)*prop43(2)
               propt2(isig1,isig3,k,2) = 
     1            clr(ifl(2,k),1,isig1)*clr(ifl(3,k),1,isig3)*prop43(1) 
     2          + clr(ifl(2,k),2,isig1)*clr(ifl(3,k),2,isig3)*prop43(2)
              
	       mat(k,isig1,isig3,4) = propt2(isig1,isig3,k,2) *
     1            ( ma(2)*clr(ifl(1,k),1,isig1) +
     2              mz(2)*clr(ifl(1,k),2,isig1)   )
     3    	+ propt1(isig1,isig3,k,2) *
     1            (  ma(1)*clr(ifl(2,k),1,isig1) +
     2               mz(1)*clr(ifl(2,k),2,isig1)   )
               if (lbox) matv(k,isig1,isig3,4) = propt2(isig1,isig3,k,2) *
     1            ( mva(2)*clr(ifl(1,k),1,isig1) +
     2              mvz(2)*clr(ifl(1,k),2,isig1)   )
     3    	+ propt1(isig1,isig3,k,2) *
     1            (  mva(1)*clr(ifl(2,k),1,isig1) +
     2               mvz(1)*clr(ifl(2,k),2,isig1)   )
            enddo

            
	    if (isig1.eq.-1 .and. isig3.eq.-1) then
               do k = 5,6
                  propt(-1,-1,k,2) = clr(3,3,-1)**2*prop43(3)
                  mat(k,-1,-1,4) = propt(-1,-1,k,2) *
     1       ( ma(1)*clr(ifl(2,k),1,-1) + ma(2)*clr(ifl(1,k),1,-1) 
     2       + mz(1)*clr(ifl(2,k),2,-1) + mz(2)*clr(ifl(1,k),2,-1) )
                  if (lbox) matv(k,-1,-1,4) = propt(-1,-1,k,2) *
     1       ( mva(1)*clr(ifl(2,k),1,-1) + mva(2)*clr(ifl(1,k),1,-1) 
     2       + mvz(1)*clr(ifl(2,k),2,-1) + mvz(2)*clr(ifl(1,k),2,-1) )
               enddo
	       
            else
               do k = 5,6
                  propt(isig1,isig3,k,2) = 0d0
                  mat(k,isig1,isig3,4) = 0d0
                  matv(k,isig1,isig3,4) = 0d0
               enddo
            endif
	    
 
c ------------------	    
	    
c box correction to lower line: polarization vectors are 
c    jqq(mu,isig1,1)=j21 with momentum    p21        and
c    aww/zww(mu)         with momentum    pza

            ma(1) = dotcc(ja(0,isig3,is3,4),jqq(0,isig1,1))
            ma(2) = dotcc(ja(0,isig3,is3,3),jqq(0,isig1,1))
            mz(1) = dotcc(jz(0,isig3,is3,4),jqq(0,isig1,1))
            mz(2) = dotcc(jz(0,isig3,is3,3),jqq(0,isig1,1))
            if (lbox) then
               mva(1) = dotcc(jva(0,isig3,is3,4),jqq(0,isig1,1))
               mva(2) = dotcc(jva(0,isig3,is3,3),jqq(0,isig1,1))
               mvz(1) = dotcc(jvz(0,isig3,is3,4),jqq(0,isig1,1))
               mvz(2) = dotcc(jvz(0,isig3,is3,3),jqq(0,isig1,1))
            endif

            do k = 1,4
               propt1(isig1,isig3,k,1) = 
     1            clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3)*prop21(1)
     2          + clr(ifl(1,k),2,isig1)*clr(ifl(3,k),2,isig3)*prop21(2)
               propt2(isig1,isig3,k,1) = 
     1            clr(ifl(1,k),1,isig1)*clr(ifl(4,k),1,isig3)*prop21(1)
     2          + clr(ifl(1,k),2,isig1)*clr(ifl(4,k),2,isig3)*prop21(2)
               mat(k,isig1,isig3,5) = propt1(isig1,isig3,k,1) *
     1            ( ma(1)*clr(ifl(4,k),1,isig3) +
     2              mz(1)*clr(ifl(4,k),2,isig3)   )
     3	         +  propt2(isig1,isig3,k,1) *
     1            ( ma(2)*clr(ifl(3,k),1,isig3) +
     2              mz(2)*clr(ifl(3,k),2,isig3)   )              
	       
	       if (lbox) matv(k,isig1,isig3,5) = propt1(isig1,isig3,k,1) *
     1            ( mva(1)*clr(ifl(4,k),1,isig3) +
     2              mvz(1)*clr(ifl(4,k),2,isig3)   )
     3	         +  propt2(isig1,isig3,k,1) *
     1            ( mva(2)*clr(ifl(3,k),1,isig3) +
     2              mvz(2)*clr(ifl(3,k),2,isig3)   )
            enddo


	    if (isig1.eq.-1 .and. isig3.eq.-1) then
               do k = 5,6
                  propt(-1,-1,k,1) = clr(ifl(1,k),3,-1)**2*prop21(3)
                  mat(k,-1,-1,5) = propt(-1,-1,k,1) *
     1       ( ma(1)*clr(ifl(4,k),1,-1) + ma(2)*clr(ifl(3,k),1,-1) 
     2       + mz(1)*clr(ifl(4,k),2,-1) + mz(2)*clr(ifl(3,k),2,-1) )
                  if (lbox) matv(k,-1,-1,5) = propt(-1,-1,k,1) *
     1       ( mva(1)*clr(ifl(4,k),1,-1) + mva(2)*clr(ifl(3,k),1,-1) 
     2       + mvz(1)*clr(ifl(4,k),2,-1) + mvz(2)*clr(ifl(3,k),2,-1) )
               enddo
            else
               do k = 5,6
                  propt(isig1,isig3,k,1) = 0d0
                  mat(k,isig1,isig3,5) = 0d0
                  matv(k,isig1,isig3,5) = 0d0
               enddo
            endif
         enddo
      enddo
c
c ------------------------------------------------------------------------
c
c  next do the box-box graphs with one A/Z emitted from the upper and the 
c  other from the lower line. These are possible for left- and righthanded 
c  quarks on both lines (in case of NC)
c
c
	do isig1 = -1,1,2
	   do isig3 = -1,1,2
c
c
c e+e- from upper line, mu+mu- from lower line:
c	
	mezz(1,isig1,isig3) = dotcc(jez(0,isig1,1),juz(0,isig3,3))
	meza(1,isig1,isig3) = dotcc(jez(0,isig1,1),jua(0,isig3,3))
	meaz(1,isig1,isig3) = dotcc(jea(0,isig1,1),juz(0,isig3,3))
	meaa(1,isig1,isig3) = dotcc(jea(0,isig1,1),jua(0,isig3,3))
c	
	mezz(2,isig1,isig3) = dotcc(jez(0,isig1,2),juz(0,isig3,3))
	meza(2,isig1,isig3) = dotcc(jez(0,isig1,2),jua(0,isig3,3))
	meaz(2,isig1,isig3) = dotcc(jea(0,isig1,2),juz(0,isig3,3))
	meaa(2,isig1,isig3) = dotcc(jea(0,isig1,2),jua(0,isig3,3))
c	
	mezz(3,isig1,isig3) = dotcc(jez(0,isig1,1),juz(0,isig3,4))
	meza(3,isig1,isig3) = dotcc(jez(0,isig1,1),jua(0,isig3,4))
	meaz(3,isig1,isig3) = dotcc(jea(0,isig1,1),juz(0,isig3,4))
	meaa(3,isig1,isig3) = dotcc(jea(0,isig1,1),jua(0,isig3,4))
c	
	mezz(4,isig1,isig3) = dotcc(jez(0,isig1,2),juz(0,isig3,4))
	meza(4,isig1,isig3) = dotcc(jez(0,isig1,2),jua(0,isig3,4))
	meaz(4,isig1,isig3) = dotcc(jea(0,isig1,2),juz(0,isig3,4))
	meaa(4,isig1,isig3) = dotcc(jea(0,isig1,2),jua(0,isig3,4))
	
c mu+mu- from upper line, e+e- from lower line:
c	
	muzz(1,isig1,isig3) = dotcc(juz(0,isig1,1),jez(0,isig3,3))
	muza(1,isig1,isig3) = dotcc(juz(0,isig1,1),jea(0,isig3,3))
	muaz(1,isig1,isig3) = dotcc(jua(0,isig1,1),jez(0,isig3,3))
	muaa(1,isig1,isig3) = dotcc(jua(0,isig1,1),jea(0,isig3,3))
c	
	muzz(2,isig1,isig3) = dotcc(juz(0,isig1,2),jez(0,isig3,3))
	muza(2,isig1,isig3) = dotcc(juz(0,isig1,2),jea(0,isig3,3))
	muaz(2,isig1,isig3) = dotcc(jua(0,isig1,2),jez(0,isig3,3))
	muaa(2,isig1,isig3) = dotcc(jua(0,isig1,2),jea(0,isig3,3))
c	
	muzz(3,isig1,isig3) = dotcc(juz(0,isig1,1),jez(0,isig3,4))
	muza(3,isig1,isig3) = dotcc(juz(0,isig1,1),jea(0,isig3,4))
	muaz(3,isig1,isig3) = dotcc(jua(0,isig1,1),jez(0,isig3,4))
	muaa(3,isig1,isig3) = dotcc(jua(0,isig1,1),jea(0,isig3,4))
c	
	muzz(4,isig1,isig3) = dotcc(juz(0,isig1,2),jez(0,isig3,4))
	muza(4,isig1,isig3) = dotcc(juz(0,isig1,2),jea(0,isig3,4))
	muaz(4,isig1,isig3) = dotcc(jua(0,isig1,2),jez(0,isig3,4))
	muaa(4,isig1,isig3) = dotcc(jua(0,isig1,2),jea(0,isig3,4))
	

c -------------------------------

	if (lbox) then
c e+e- from upper line, mu+mu- from lower line
c box corr. to upper line:
c	
	mvezz(1,1,isig1,isig3) = dotcc(jvez(0,isig1,is1,1),juz(0,isig3,3))
	mveza(1,1,isig1,isig3) = dotcc(jvez(0,isig1,is1,1),jua(0,isig3,3))
	mveaz(1,1,isig1,isig3) = dotcc(jvea(0,isig1,is1,1),juz(0,isig3,3))
	mveaa(1,1,isig1,isig3) = dotcc(jvea(0,isig1,is1,1),jua(0,isig3,3))
c	
	mvezz(2,1,isig1,isig3) = dotcc(jvez(0,isig1,is1,2),juz(0,isig3,3))
	mveza(2,1,isig1,isig3) = dotcc(jvez(0,isig1,is1,2),jua(0,isig3,3))
	mveaz(2,1,isig1,isig3) = dotcc(jvea(0,isig1,is1,2),juz(0,isig3,3))
	mveaa(2,1,isig1,isig3) = dotcc(jvea(0,isig1,is1,2),jua(0,isig3,3))
c	
	mvezz(3,1,isig1,isig3) = dotcc(jvez(0,isig1,is1,1),juz(0,isig3,4))
	mveza(3,1,isig1,isig3) = dotcc(jvez(0,isig1,is1,1),jua(0,isig3,4))
	mveaz(3,1,isig1,isig3) = dotcc(jvea(0,isig1,is1,1),juz(0,isig3,4))
	mveaa(3,1,isig1,isig3) = dotcc(jvea(0,isig1,is1,1),jua(0,isig3,4))
c	
	mvezz(4,1,isig1,isig3) = dotcc(jvez(0,isig1,is1,2),juz(0,isig3,4))
	mveza(4,1,isig1,isig3) = dotcc(jvez(0,isig1,is1,2),jua(0,isig3,4))
	mveaz(4,1,isig1,isig3) = dotcc(jvea(0,isig1,is1,2),juz(0,isig3,4))
	mveaa(4,1,isig1,isig3) = dotcc(jvea(0,isig1,is1,2),jua(0,isig3,4))
c	
c box corr. to lower line:
c
	mvezz(1,2,isig1,isig3) = dotcc(jez(0,isig1,1),jvuz(0,isig3,is3,3))
	mveza(1,2,isig1,isig3) = dotcc(jez(0,isig1,1),jvua(0,isig3,is3,3))
	mveaz(1,2,isig1,isig3) = dotcc(jea(0,isig1,1),jvuz(0,isig3,is3,3))
	mveaa(1,2,isig1,isig3) = dotcc(jea(0,isig1,1),jvua(0,isig3,is3,3))
c	
	mvezz(2,2,isig1,isig3) = dotcc(jez(0,isig1,2),jvuz(0,isig3,is3,3))
	mveza(2,2,isig1,isig3) = dotcc(jez(0,isig1,2),jvua(0,isig3,is3,3))
	mveaz(2,2,isig1,isig3) = dotcc(jea(0,isig1,2),jvuz(0,isig3,is3,3))
	mveaa(2,2,isig1,isig3) = dotcc(jea(0,isig1,2),jvua(0,isig3,is3,3))
c	
	mvezz(3,2,isig1,isig3) = dotcc(jez(0,isig1,1),jvuz(0,isig3,is3,4))
	mveza(3,2,isig1,isig3) = dotcc(jez(0,isig1,1),jvua(0,isig3,is3,4))
	mveaz(3,2,isig1,isig3) = dotcc(jea(0,isig1,1),jvuz(0,isig3,is3,4))
	mveaa(3,2,isig1,isig3) = dotcc(jea(0,isig1,1),jvua(0,isig3,is3,4))
c	
	mvezz(4,2,isig1,isig3) = dotcc(jez(0,isig1,2),jvuz(0,isig3,is3,4))
	mveza(4,2,isig1,isig3) = dotcc(jez(0,isig1,2),jvua(0,isig3,is3,4))
	mveaz(4,2,isig1,isig3) = dotcc(jea(0,isig1,2),jvuz(0,isig3,is3,4))
	mveaa(4,2,isig1,isig3) = dotcc(jea(0,isig1,2),jvua(0,isig3,is3,4))
	
c mu+mu- from upper line, e+e- from lower line:
c box corr. to upper line:	
c	
	mvuzz(1,1,isig1,isig3) = dotcc(jvuz(0,isig1,is1,1),jez(0,isig3,3))
	mvuza(1,1,isig1,isig3) = dotcc(jvuz(0,isig1,is1,1),jea(0,isig3,3))
	mvuaz(1,1,isig1,isig3) = dotcc(jvua(0,isig1,is1,1),jez(0,isig3,3))
	mvuaa(1,1,isig1,isig3) = dotcc(jvua(0,isig1,is1,1),jea(0,isig3,3))
c	
	mvuzz(2,1,isig1,isig3) = dotcc(jvuz(0,isig1,is1,2),jez(0,isig3,3))
	mvuza(2,1,isig1,isig3) = dotcc(jvuz(0,isig1,is1,2),jea(0,isig3,3))
	mvuaz(2,1,isig1,isig3) = dotcc(jvua(0,isig1,is1,2),jez(0,isig3,3))
	mvuaa(2,1,isig1,isig3) = dotcc(jvua(0,isig1,is1,2),jea(0,isig3,3))
c	
	mvuzz(3,1,isig1,isig3) = dotcc(jvuz(0,isig1,is1,1),jez(0,isig3,4))
	mvuza(3,1,isig1,isig3) = dotcc(jvuz(0,isig1,is1,1),jea(0,isig3,4))
	mvuaz(3,1,isig1,isig3) = dotcc(jvua(0,isig1,is1,1),jez(0,isig3,4))
	mvuaa(3,1,isig1,isig3) = dotcc(jvua(0,isig1,is1,1),jea(0,isig3,4))
c	
	mvuzz(4,1,isig1,isig3) = dotcc(jvuz(0,isig1,is1,2),jez(0,isig3,4))
	mvuza(4,1,isig1,isig3) = dotcc(jvuz(0,isig1,is1,2),jea(0,isig3,4))
	mvuaz(4,1,isig1,isig3) = dotcc(jvua(0,isig1,is1,2),jez(0,isig3,4))
	mvuaa(4,1,isig1,isig3) = dotcc(jvua(0,isig1,is1,2),jea(0,isig3,4))
	
c box corr. to lower line:	
c	
	mvuzz(1,2,isig1,isig3) = dotcc(juz(0,isig1,1),jvez(0,isig3,is3,3))
	mvuza(1,2,isig1,isig3) = dotcc(juz(0,isig1,1),jvea(0,isig3,is3,3))
	mvuaz(1,2,isig1,isig3) = dotcc(jua(0,isig1,1),jvez(0,isig3,is3,3))
	mvuaa(1,2,isig1,isig3) = dotcc(jua(0,isig1,1),jvea(0,isig3,is3,3))
c	
	mvuzz(2,2,isig1,isig3) = dotcc(juz(0,isig1,2),jvez(0,isig3,is3,3))
	mvuza(2,2,isig1,isig3) = dotcc(juz(0,isig1,2),jvea(0,isig3,is3,3))
	mvuaz(2,2,isig1,isig3) = dotcc(jua(0,isig1,2),jvez(0,isig3,is3,3))
	mvuaa(2,2,isig1,isig3) = dotcc(jua(0,isig1,2),jvea(0,isig3,is3,3))
c	
	mvuzz(3,2,isig1,isig3) = dotcc(juz(0,isig1,1),jvez(0,isig3,is3,4))
	mvuza(3,2,isig1,isig3) = dotcc(juz(0,isig1,1),jvea(0,isig3,is3,4))
	mvuaz(3,2,isig1,isig3) = dotcc(jua(0,isig1,1),jvez(0,isig3,is3,4))
	mvuaa(3,2,isig1,isig3) = dotcc(jua(0,isig1,1),jvea(0,isig3,is3,4))
c	
	mvuzz(4,2,isig1,isig3) = dotcc(juz(0,isig1,2),jvez(0,isig3,is3,4))
	mvuza(4,2,isig1,isig3) = dotcc(juz(0,isig1,2),jvea(0,isig3,is3,4))
	mvuaz(4,2,isig1,isig3) = dotcc(jua(0,isig1,2),jvez(0,isig3,is3,4))
	mvuaa(4,2,isig1,isig3) = dotcc(jua(0,isig1,2),jvea(0,isig3,is3,4))
	
	endif !lbox
	
c
c  for the q^mu*q^nu/M_V^2 terms in the gauge boson propagators we need
	
	do i = 2,3 ! type of exchanged boson (Z,W), zero for A
	
	zezz(isig1,isig3,i) = jj21ze(isig1)*jj43zu(isig3)*zm2i(i)
	zeza(isig1,isig3,i) = jj21ze(isig1)*jj43au(isig3)*zm2i(i)
	zeaz(isig1,isig3,i) = jj21ae(isig1)*jj43zu(isig3)*zm2i(i)
	zeaa(isig1,isig3,i) = jj21ae(isig1)*jj43au(isig3)*zm2i(i)

	zuzz(isig1,isig3,i) = jj21zu(isig1)*jj43ze(isig3)*zm2i(i)
	zuza(isig1,isig3,i) = jj21zu(isig1)*jj43ae(isig3)*zm2i(i)
	zuaz(isig1,isig3,i) = jj21au(isig1)*jj43ze(isig3)*zm2i(i)
	zuaa(isig1,isig3,i) = jj21au(isig1)*jj43ae(isig3)*zm2i(i)
	
	enddo !i


c -------------------------------
c
      do k = 1,4 
             
     
      	      propbbez(k,isig1,isig3) = 
     1 	     	prop_ee(2)*clr(ifl(2,k),2,isig1)*clr(ifl(3,k),2,isig3)
     
      	      propbbuz(k,isig1,isig3) = 
     1	     	prop_uu(2)*clr(ifl(2,k),2,isig1)*clr(ifl(3,k),2,isig3)
     
      	      propbbea(k,isig1,isig3) = 
     1 	     	prop_ee(1)*clr(ifl(2,k),1,isig1)*clr(ifl(3,k),1,isig3)
     
      	      propbbua(k,isig1,isig3) = 
     1	    	prop_uu(1)*clr(ifl(2,k),1,isig1)*clr(ifl(3,k),1,isig3)     

     
c for NC make use of: ifl(1,k) = ifl (2,k) and ifl(3,k) = ifl (4,k) 
c 	gauge terms zezz etc cancel in sum of NC contributions from leg1&2/3&4 
c       (not for CC)
	     
      mat(k,isig1,isig3,6) = 
     1 ((mezz(1,isig1,isig3)+mezz(2,isig1,isig3))*
     1		clr(ifl(1,k),2,isig1)*clr(ifl(3,k),2,isig3) + 
     1	(meza(1,isig1,isig3)+meza(2,isig1,isig3))*
     1		clr(ifl(1,k),2,isig1)*clr(ifl(3,k),1,isig3)+ 
     1	(meaz(1,isig1,isig3)+meaz(2,isig1,isig3))*
     1		clr(ifl(1,k),1,isig1)*clr(ifl(3,k),2,isig3)+ 
     1	(meaa(1,isig1,isig3)+meaa(2,isig1,isig3))*
     1		clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3))*
     1	      propbbea(k,isig1,isig3)
     1	+
     1 (((mezz(1,isig1,isig3)+zezz(isig1,isig3,2))+
     1	 (mezz(2,isig1,isig3)-zezz(isig1,isig3,2)))*
     1		clr(ifl(1,k),2,isig1)*clr(ifl(3,k),2,isig3) + 
     1	((meza(1,isig1,isig3)+zeza(isig1,isig3,2))+
     1	 (meza(2,isig1,isig3)-zeza(isig1,isig3,2)))*
     1		clr(ifl(1,k),2,isig1)*clr(ifl(3,k),1,isig3)+ 
     1	((meaz(1,isig1,isig3)+zeaz(isig1,isig3,2))+
     1	 (meaz(2,isig1,isig3)-zeaz(isig1,isig3,2)))*
     1		clr(ifl(1,k),1,isig1)*clr(ifl(3,k),2,isig3)+ 
     1	((meaa(1,isig1,isig3)+zeaa(isig1,isig3,2))+
     1	 (meaa(2,isig1,isig3)-zeaa(isig1,isig3,2)))*
     1		clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3))*
     1	      propbbez(k,isig1,isig3)
     1 +
     1 ((muzz(1,isig1,isig3)+muzz(2,isig1,isig3))*
     1		clr(ifl(1,k),2,isig1)*clr(ifl(3,k),2,isig3) + 
     1	(muza(1,isig1,isig3)+muza(2,isig1,isig3))*
     1		clr(ifl(1,k),2,isig1)*clr(ifl(3,k),1,isig3)+ 
     1	(muaz(1,isig1,isig3)+muaz(2,isig1,isig3))*
     1		clr(ifl(1,k),1,isig1)*clr(ifl(3,k),2,isig3)+ 
     1	(muaa(1,isig1,isig3)+muaa(2,isig1,isig3))*
     1		clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3))*
     1	      propbbua(k,isig1,isig3)
     1 +
     1 (((muzz(1,isig1,isig3)+zuzz(isig1,isig3,2))+
     1  ( muzz(2,isig1,isig3)-zuzz(isig1,isig3,2)))*
     1		clr(ifl(1,k),2,isig1)*clr(ifl(3,k),2,isig3) + 
     1	((muza(1,isig1,isig3)+zuza(isig1,isig3,2))+
     1	 (muza(2,isig1,isig3)-zuza(isig1,isig3,2)))*
     1		clr(ifl(1,k),2,isig1)*clr(ifl(3,k),1,isig3)+ 
     1	((muaz(1,isig1,isig3)+zuaz(isig1,isig3,2))+
     1	 (muaz(2,isig1,isig3)-zuaz(isig1,isig3,2)))*
     1		clr(ifl(1,k),1,isig1)*clr(ifl(3,k),2,isig3)+ 
     1	((muaa(1,isig1,isig3)+zuaa(isig1,isig3,2))+
     1	 (muaa(2,isig1,isig3)-zuaa(isig1,isig3,2)))*
     1		clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3))*
     1	      propbbuz(k,isig1,isig3)
     3 +
     3 ((mezz(3,isig1,isig3)+mezz(4,isig1,isig3))*
     3	       clr(ifl(1,k),2,isig1)*clr(ifl(4,k),2,isig3) + 
     3  (meza(3,isig1,isig3)+meza(4,isig1,isig3))*
     3	       clr(ifl(1,k),2,isig1)*clr(ifl(4,k),1,isig3)+ 
     3  (meaz(3,isig1,isig3)+meaz(4,isig1,isig3))*
     3	       clr(ifl(1,k),1,isig1)*clr(ifl(4,k),2,isig3)+ 
     3  (meaa(3,isig1,isig3)+meaa(4,isig1,isig3))*
     3	       clr(ifl(1,k),1,isig1)*clr(ifl(4,k),1,isig3))*
     3	      propbbea(k,isig1,isig3)
     3 +
     3 (((mezz(3,isig1,isig3)-zezz(isig1,isig3,2))+
     3   (mezz(4,isig1,isig3)+zezz(isig1,isig3,2)))*
     3	       clr(ifl(1,k),2,isig1)*clr(ifl(4,k),2,isig3) + 
     3  ((meza(3,isig1,isig3)-zeza(isig1,isig3,2))+
     3   (meza(4,isig1,isig3)+zeza(isig1,isig3,2)))*
     3	       clr(ifl(1,k),2,isig1)*clr(ifl(4,k),1,isig3)+ 
     3  ((meaz(3,isig1,isig3)-zeaz(isig1,isig3,2))+
     3   (meaz(4,isig1,isig3)+zeaz(isig1,isig3,2)))*
     3	       clr(ifl(1,k),1,isig1)*clr(ifl(4,k),2,isig3)+ 
     3  ((meaa(3,isig1,isig3)-zeaa(isig1,isig3,2))+
     3   (meaa(4,isig1,isig3)+zeaa(isig1,isig3,2)))*
     3	       clr(ifl(1,k),1,isig1)*clr(ifl(4,k),1,isig3))*
     3	      propbbez(k,isig1,isig3)
     3 +
     3 ((muzz(3,isig1,isig3)+muzz(4,isig1,isig3))*
     3	       clr(ifl(1,k),2,isig1)*clr(ifl(4,k),2,isig3) + 
     3  (muza(3,isig1,isig3)+muza(4,isig1,isig3))*
     3	       clr(ifl(1,k),2,isig1)*clr(ifl(4,k),1,isig3)+ 
     3  (muaz(3,isig1,isig3)+muaz(4,isig1,isig3))*
     3	       clr(ifl(1,k),1,isig1)*clr(ifl(4,k),2,isig3)+ 
     3  (muaa(3,isig1,isig3)+muaa(4,isig1,isig3))*
     3	       clr(ifl(1,k),1,isig1)*clr(ifl(4,k),1,isig3))*
     3	      propbbua(k,isig1,isig3)
     3 +
     3 (((muzz(3,isig1,isig3)-zuzz(isig1,isig3,2))+
     3   (muzz(4,isig1,isig3)+zuzz(isig1,isig3,2)))*
     3	       clr(ifl(1,k),2,isig1)*clr(ifl(4,k),2,isig3) + 
     3  ((muza(3,isig1,isig3)-zuza(isig1,isig3,2))+
     3   (muza(4,isig1,isig3)+zuza(isig1,isig3,2)))*
     3	       clr(ifl(1,k),2,isig1)*clr(ifl(4,k),1,isig3)+ 
     3  ((muaz(3,isig1,isig3)-zuaz(isig1,isig3,2))+
     3   (muaz(4,isig1,isig3)+zuaz(isig1,isig3,2)))*
     3	       clr(ifl(1,k),1,isig1)*clr(ifl(4,k),2,isig3)+ 
     3  ((muaa(3,isig1,isig3)-zuaa(isig1,isig3,2))+
     3   (muaa(4,isig1,isig3)+zuaa(isig1,isig3,2)))*
     3	       clr(ifl(1,k),1,isig1)*clr(ifl(4,k),1,isig3))*
     3	      propbbuz(k,isig1,isig3)   
      

      enddo !k

	   enddo ! isig3
	enddo ! isig1
c
c -------------------------	
c
      do k = 5,6 ! charged current
      isig1 = -1
      isig3 = -1
            	    
      propbbe =  prop_ee(3)*clr(3,3,-1)**2
      propbbu =  prop_uu(3)*clr(3,3,-1)**2
             
      mat(k,isig1,isig3,6) = 
     1 ((mezz(1,isig1,isig3)+zezz(isig1,isig3,3))*
     1		clr(ifl(1,k),2,isig1)*clr(ifl(3,k),2,isig3) + 
     1	(meza(1,isig1,isig3)+zeza(isig1,isig3,3))*
     1		clr(ifl(1,k),2,isig1)*clr(ifl(3,k),1,isig3)+ 
     1	(meaz(1,isig1,isig3)+zeaz(isig1,isig3,3))*
     1		clr(ifl(1,k),1,isig1)*clr(ifl(3,k),2,isig3)+ 
     1	(meaa(1,isig1,isig3)+zeaa(isig1,isig3,3))*
     1		clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3))*
     1		propbbe     
     2 +
     2 ((mezz(2,isig1,isig3)-zezz(isig1,isig3,3))*
     2		clr(ifl(2,k),2,isig1)*clr(ifl(3,k),2,isig3) + 
     2	(meza(2,isig1,isig3)-zeza(isig1,isig3,3))*
     2		clr(ifl(2,k),2,isig1)*clr(ifl(3,k),1,isig3)+ 
     2	(meaz(2,isig1,isig3)-zeaz(isig1,isig3,3))*
     2		clr(ifl(2,k),1,isig1)*clr(ifl(3,k),2,isig3)+ 
     2	(meaa(2,isig1,isig3)-zeaa(isig1,isig3,3))*
     2		clr(ifl(2,k),1,isig1)*clr(ifl(3,k),1,isig3))*
     2		propbbe     
     3 +
     3 ((mezz(3,isig1,isig3)-zezz(isig1,isig3,3))*
     3	       clr(ifl(1,k),2,isig1)*clr(ifl(4,k),2,isig3) + 
     3  (meza(3,isig1,isig3)-zeza(isig1,isig3,3))*
     3	       clr(ifl(1,k),2,isig1)*clr(ifl(4,k),1,isig3)+ 
     3  (meaz(3,isig1,isig3)-zeaz(isig1,isig3,3))*
     3	       clr(ifl(1,k),1,isig1)*clr(ifl(4,k),2,isig3)+ 
     3  (meaa(3,isig1,isig3)-zeaa(isig1,isig3,3))*
     3	       clr(ifl(1,k),1,isig1)*clr(ifl(4,k),1,isig3))*
     3		propbbe     
     4 +
     4 ((mezz(4,isig1,isig3)+zezz(isig1,isig3,3))*
     4	      clr(ifl(2,k),2,isig1)*clr(ifl(4,k),2,isig3) + 
     4  (meza(4,isig1,isig3)+zeza(isig1,isig3,3))*
     4	      clr(ifl(2,k),2,isig1)*clr(ifl(4,k),1,isig3)+ 
     4  (meaz(4,isig1,isig3)+zeaz(isig1,isig3,3))*
     4	      clr(ifl(2,k),1,isig1)*clr(ifl(4,k),2,isig3)+ 
     4  (meaa(4,isig1,isig3)+zeaa(isig1,isig3,3))*
     4	      clr(ifl(2,k),1,isig1)*clr(ifl(4,k),1,isig3))*
     4		propbbe 
     1	+        
     1 ((muzz(1,isig1,isig3)+zuzz(isig1,isig3,3))*
     1	       clr(ifl(1,k),2,isig1)*clr(ifl(3,k),2,isig3) + 
     1	(muza(1,isig1,isig3)+zuza(isig1,isig3,3))*
     1	       clr(ifl(1,k),2,isig1)*clr(ifl(3,k),1,isig3)+ 
     1	(muaz(1,isig1,isig3)+zuaz(isig1,isig3,3))*
     1	       clr(ifl(1,k),1,isig1)*clr(ifl(3,k),2,isig3)+ 
     1	(muaa(1,isig1,isig3)+zuaa(isig1,isig3,3))*
     1	       clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3))*
     1	       propbbu    
     2 +
     2 ((muzz(2,isig1,isig3)-zuzz(isig1,isig3,3))*
     2	       clr(ifl(2,k),2,isig1)*clr(ifl(3,k),2,isig3) + 
     2	(muza(2,isig1,isig3)-zuza(isig1,isig3,3))*
     2	       clr(ifl(2,k),2,isig1)*clr(ifl(3,k),1,isig3)+ 
     2	(muaz(2,isig1,isig3)-zuaz(isig1,isig3,3))*
     2	       clr(ifl(2,k),1,isig1)*clr(ifl(3,k),2,isig3)+ 
     2	(muaa(2,isig1,isig3)-zuaa(isig1,isig3,3))*
     2	       clr(ifl(2,k),1,isig1)*clr(ifl(3,k),1,isig3))*
     2	       propbbu     
     3 +
     3 ((muzz(3,isig1,isig3)-zuzz(isig1,isig3,3))*
     3	      clr(ifl(1,k),2,isig1)*clr(ifl(4,k),2,isig3) + 
     3  (muza(3,isig1,isig3)-zuza(isig1,isig3,3))*
     3	      clr(ifl(1,k),2,isig1)*clr(ifl(4,k),1,isig3)+ 
     3  (muaz(3,isig1,isig3)-zuaz(isig1,isig3,3))*
     3	      clr(ifl(1,k),1,isig1)*clr(ifl(4,k),2,isig3)+ 
     3  (muaa(3,isig1,isig3)-zuaa(isig1,isig3,3))*
     3	      clr(ifl(1,k),1,isig1)*clr(ifl(4,k),1,isig3))*
     3	       propbbu    
     4 +
     4 ((muzz(4,isig1,isig3)+zuzz(isig1,isig3,3))*
     4	     clr(ifl(2,k),2,isig1)*clr(ifl(4,k),2,isig3) + 
     4  (muza(4,isig1,isig3)+zuza(isig1,isig3,3))*
     4	     clr(ifl(2,k),2,isig1)*clr(ifl(4,k),1,isig3)+ 
     4  (muaz(4,isig1,isig3)+zuaz(isig1,isig3,3))*
     4	     clr(ifl(2,k),1,isig1)*clr(ifl(4,k),2,isig3)+ 
     4  (muaa(4,isig1,isig3)+zuaa(isig1,isig3,3))*
     4	      clr(ifl(2,k),1,isig1)*clr(ifl(4,k),1,isig3))*
     1		propbbu
          
      enddo !k
c
c---------------------------------------------------
c      
      if (lbox) then ! no need to include gauge terms (contraction of box with 
      		     !		momentum -> triangle)	
         do i = 1,2      ! 1 is for upper line, 2 for lower line QCD correction
            jj = 3 + 3*i ! stored in matv(...,6) and matv(...,9) respectively 
	    do isig1 = -1,1,2
	      do isig3 = -1,1,2

                do k = 1,4 
      
      matv(k,isig1,isig3,jj) = 
     1 ((mvezz(1,i,isig1,isig3)+mvezz(2,i,isig1,isig3))*
     1		clr(ifl(1,k),2,isig1)*clr(ifl(3,k),2,isig3) + 
     1	(mveza(1,i,isig1,isig3)+mveza(2,i,isig1,isig3))*
     1		clr(ifl(1,k),2,isig1)*clr(ifl(3,k),1,isig3)+ 
     1	(mveaz(1,i,isig1,isig3)+mveaz(2,i,isig1,isig3))*
     1		clr(ifl(1,k),1,isig1)*clr(ifl(3,k),2,isig3)+ 
     1	(mveaa(1,i,isig1,isig3)+mveaa(2,i,isig1,isig3))*
     1		clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3))*
     1	      propbbea(k,isig1,isig3)
     1	+
     1 (((mvezz(1,i,isig1,isig3))+(mvezz(2,i,isig1,isig3)))*
     1		clr(ifl(1,k),2,isig1)*clr(ifl(3,k),2,isig3) + 
     1	((mveza(1,i,isig1,isig3))+(mveza(2,i,isig1,isig3)))*
     1		clr(ifl(1,k),2,isig1)*clr(ifl(3,k),1,isig3)+ 
     1	((mveaz(1,i,isig1,isig3))+(mveaz(2,i,isig1,isig3)))*
     1		clr(ifl(1,k),1,isig1)*clr(ifl(3,k),2,isig3)+ 
     1	((mveaa(1,i,isig1,isig3))+(mveaa(2,i,isig1,isig3)))*
     1		clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3))*
     1	      propbbez(k,isig1,isig3)
     1 +
     1 ((mvuzz(1,i,isig1,isig3)+mvuzz(2,i,isig1,isig3))*
     1		clr(ifl(1,k),2,isig1)*clr(ifl(3,k),2,isig3) + 
     1	(mvuza(1,i,isig1,isig3)+mvuza(2,i,isig1,isig3))*
     1		clr(ifl(1,k),2,isig1)*clr(ifl(3,k),1,isig3)+ 
     1	(mvuaz(1,i,isig1,isig3)+mvuaz(2,i,isig1,isig3))*
     1		clr(ifl(1,k),1,isig1)*clr(ifl(3,k),2,isig3)+ 
     1	(mvuaa(1,i,isig1,isig3)+mvuaa(2,i,isig1,isig3))*
     1		clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3))*
     1	      propbbua(k,isig1,isig3)
     1 +
     1 (((mvuzz(1,i,isig1,isig3))+( mvuzz(2,i,isig1,isig3)))*
     1		clr(ifl(1,k),2,isig1)*clr(ifl(3,k),2,isig3) + 
     1	((mvuza(1,i,isig1,isig3))+(mvuza(2,i,isig1,isig3)))*
     1		clr(ifl(1,k),2,isig1)*clr(ifl(3,k),1,isig3)+ 
     1	((mvuaz(1,i,isig1,isig3))+(mvuaz(2,i,isig1,isig3)))*
     1		clr(ifl(1,k),1,isig1)*clr(ifl(3,k),2,isig3)+ 
     1	((mvuaa(1,i,isig1,isig3))+(mvuaa(2,i,isig1,isig3)))*
     1		clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3))*
     1	      propbbuz(k,isig1,isig3)
     3 +
     3 ((mvezz(3,i,isig1,isig3)+mvezz(4,i,isig1,isig3))*
     3	       clr(ifl(1,k),2,isig1)*clr(ifl(4,k),2,isig3) + 
     3  (mveza(3,i,isig1,isig3)+mveza(4,i,isig1,isig3))*
     3	       clr(ifl(1,k),2,isig1)*clr(ifl(4,k),1,isig3)+ 
     3  (mveaz(3,i,isig1,isig3)+mveaz(4,i,isig1,isig3))*
     3	       clr(ifl(1,k),1,isig1)*clr(ifl(4,k),2,isig3)+ 
     3  (mveaa(3,i,isig1,isig3)+mveaa(4,i,isig1,isig3))*
     3	       clr(ifl(1,k),1,isig1)*clr(ifl(4,k),1,isig3))*
     3	      propbbea(k,isig1,isig3)
     3 +
     3 (((mvezz(3,i,isig1,isig3))+(mvezz(4,i,isig1,isig3)))*
     3	       clr(ifl(1,k),2,isig1)*clr(ifl(4,k),2,isig3) + 
     3  ((mveza(3,i,isig1,isig3))+(mveza(4,i,isig1,isig3)))*
     3	       clr(ifl(1,k),2,isig1)*clr(ifl(4,k),1,isig3)+ 
     3  ((mveaz(3,i,isig1,isig3))+(mveaz(4,i,isig1,isig3)))*
     3	       clr(ifl(1,k),1,isig1)*clr(ifl(4,k),2,isig3)+ 
     3  ((mveaa(3,i,isig1,isig3))+(mveaa(4,i,isig1,isig3)))*
     3	       clr(ifl(1,k),1,isig1)*clr(ifl(4,k),1,isig3))*
     3	      propbbez(k,isig1,isig3)
     3 +
     3 ((mvuzz(3,i,isig1,isig3)+mvuzz(4,i,isig1,isig3))*
     3	       clr(ifl(1,k),2,isig1)*clr(ifl(4,k),2,isig3) + 
     3  (mvuza(3,i,isig1,isig3)+mvuza(4,i,isig1,isig3))*
     3	       clr(ifl(1,k),2,isig1)*clr(ifl(4,k),1,isig3)+ 
     3  (mvuaz(3,i,isig1,isig3)+mvuaz(4,i,isig1,isig3))*
     3	       clr(ifl(1,k),1,isig1)*clr(ifl(4,k),2,isig3)+ 
     3  (mvuaa(3,i,isig1,isig3)+mvuaa(4,i,isig1,isig3))*
     3	       clr(ifl(1,k),1,isig1)*clr(ifl(4,k),1,isig3))*
     3	      propbbua(k,isig1,isig3)
     3 +
     3 (((mvuzz(3,i,isig1,isig3))+(mvuzz(4,i,isig1,isig3)))*
     3	       clr(ifl(1,k),2,isig1)*clr(ifl(4,k),2,isig3) + 
     3  ((mvuza(3,i,isig1,isig3))+(mvuza(4,i,isig1,isig3)))*
     3	       clr(ifl(1,k),2,isig1)*clr(ifl(4,k),1,isig3)+ 
     3  ((mvuaz(3,i,isig1,isig3))+(mvuaz(4,i,isig1,isig3)))*
     3	       clr(ifl(1,k),1,isig1)*clr(ifl(4,k),2,isig3)+ 
     3  ((mvuaa(3,i,isig1,isig3))+(mvuaa(4,i,isig1,isig3)))*
     3	       clr(ifl(1,k),1,isig1)*clr(ifl(4,k),1,isig3))*
     3	      propbbuz(k,isig1,isig3)


       		enddo !k
       
       		enddo !isig3
       	    enddo !isig1

      do k = 5,6
      isig1 = -1
      isig3 = -1
      
      matv(k,isig1,isig3,jj) = 
     1 (mvezz(1,i,isig1,isig3)*
     1		clr(ifl(1,k),2,isig1)*clr(ifl(3,k),2,isig3) + 
     1	mveza(1,i,isig1,isig3)*
     1		clr(ifl(1,k),2,isig1)*clr(ifl(3,k),1,isig3)+ 
     1	mveaz(1,i,isig1,isig3)*
     1		clr(ifl(1,k),1,isig1)*clr(ifl(3,k),2,isig3)+ 
     1	mveaa(1,i,isig1,isig3)*
     1		clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3))*
     1	    propbbe
     2 +
     2 (mvezz(2,i,isig1,isig3)*
     2		clr(ifl(2,k),2,isig1)*clr(ifl(3,k),2,isig3) + 
     2	mveza(2,i,isig1,isig3)*
     2		clr(ifl(2,k),2,isig1)*clr(ifl(3,k),1,isig3)+ 
     2	mveaz(2,i,isig1,isig3)*
     2		clr(ifl(2,k),1,isig1)*clr(ifl(3,k),2,isig3)+ 
     2	mveaa(2,i,isig1,isig3)*
     2		clr(ifl(2,k),1,isig1)*clr(ifl(3,k),1,isig3))*
     2	    propbbe
     3 +
     3 (mvezz(3,i,isig1,isig3)*
     3	       clr(ifl(1,k),2,isig1)*clr(ifl(4,k),2,isig3) + 
     3  mveza(3,i,isig1,isig3)*
     3	       clr(ifl(1,k),2,isig1)*clr(ifl(4,k),1,isig3)+ 
     3  mveaz(3,i,isig1,isig3)*
     3	       clr(ifl(1,k),1,isig1)*clr(ifl(4,k),2,isig3)+ 
     3  mveaa(3,i,isig1,isig3)*
     3	       clr(ifl(1,k),1,isig1)*clr(ifl(4,k),1,isig3))*
     3	    propbbe
     4 +
     4 (mvezz(4,i,isig1,isig3)*
     4	      clr(ifl(2,k),2,isig1)*clr(ifl(4,k),2,isig3) + 
     4  mveza(4,i,isig1,isig3)*
     4	      clr(ifl(2,k),2,isig1)*clr(ifl(4,k),1,isig3)+ 
     4  mveaz(4,i,isig1,isig3)*
     4	      clr(ifl(2,k),1,isig1)*clr(ifl(4,k),2,isig3)+ 
     4  mveaa(4,i,isig1,isig3)*
     4	      clr(ifl(2,k),1,isig1)*clr(ifl(4,k),1,isig3))*
     1	    propbbe
     1	+
     1 (mvuzz(1,i,isig1,isig3)*
     1		clr(ifl(1,k),2,isig1)*clr(ifl(3,k),2,isig3) + 
     1	mvuza(1,i,isig1,isig3)*
     1		clr(ifl(1,k),2,isig1)*clr(ifl(3,k),1,isig3)+ 
     1	mvuaz(1,i,isig1,isig3)*
     1		clr(ifl(1,k),1,isig1)*clr(ifl(3,k),2,isig3)+ 
     1	mvuaa(1,i,isig1,isig3)*
     1		clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3))*
     1	    propbbu
     2 +
     2 (mvuzz(2,i,isig1,isig3)*
     2		clr(ifl(2,k),2,isig1)*clr(ifl(3,k),2,isig3) + 
     2	mvuza(2,i,isig1,isig3)*
     2		clr(ifl(2,k),2,isig1)*clr(ifl(3,k),1,isig3)+ 
     2	mvuaz(2,i,isig1,isig3)*
     2		clr(ifl(2,k),1,isig1)*clr(ifl(3,k),2,isig3)+ 
     2	mvuaa(2,i,isig1,isig3)*
     2		clr(ifl(2,k),1,isig1)*clr(ifl(3,k),1,isig3))*
     2	    propbbu
     3 +
     3 (mvuzz(3,i,isig1,isig3)*
     3	       clr(ifl(1,k),2,isig1)*clr(ifl(4,k),2,isig3) + 
     3  mvuza(3,i,isig1,isig3)*
     3	       clr(ifl(1,k),2,isig1)*clr(ifl(4,k),1,isig3)+ 
     3  mvuaz(3,i,isig1,isig3)*
     3	       clr(ifl(1,k),1,isig1)*clr(ifl(4,k),2,isig3)+ 
     3  mvuaa(3,i,isig1,isig3)*
     3	       clr(ifl(1,k),1,isig1)*clr(ifl(4,k),1,isig3))*
     3	    propbbu
     4 +
     4 (mvuzz(4,i,isig1,isig3)*
     4	      clr(ifl(2,k),2,isig1)*clr(ifl(4,k),2,isig3) + 
     4  mvuza(4,i,isig1,isig3)*
     4	      clr(ifl(2,k),2,isig1)*clr(ifl(4,k),1,isig3)+ 
     4  mvuaz(4,i,isig1,isig3)*
     4	      clr(ifl(2,k),1,isig1)*clr(ifl(4,k),2,isig3)+ 
     4  mvuaa(4,i,isig1,isig3)*
     4	      clr(ifl(2,k),1,isig1)*clr(ifl(4,k),1,isig3))*
     4	    propbbu
     
      enddo !k

	 enddo !i
       endif !lbox	 


c ------------------------------------------------------------------------------
c
c and now, finally, the pentagon contributions, i.e. two W"s emitted from the
c  same quark line
c
c prerequisites for virtual corrections:

#ifdef WITH_NLO 	
cfc       if(lpent) then
cfc         if(.not. lzs(js1)) then
cfc            print*, "PEN re-used in qqzaqq.F"
cfc         ENDIF
cfc       endif

      if (lpent .and. lzs(js1)) then     ! need new pentagon graphs for 12 line

c       if(.true.) then
c       if(.true.) then
         musq = -p21(4)
         if (lpent) then

           do mu = 0,3
            pvirtcalc(mu,1) = p(mu,1)
            pvirtcalc(mu,2) = -1d0 * p(mu,2)
            pvirtcalc(mu,3) = -1d0 * qe(mu)
            pvirtcalc(mu,4) = -1d0 * qa(mu)
            pvirtcalc(mu,5) = 1d0 * (p(mu,3) - p(mu,4))
            temp=pvirtcalc(mu,1)+pvirtcalc(mu,2)+pvirtcalc(mu,3)+pvirtcalc(mu,4)+pvirtcalc(mu,5)
cfc            print*,"Mom. conserv.", temp
cfc           print*, "pvirtcalc",  pvirtcalc(mu,5)
         enddo
           scale = scaleupperline

            do k = 1,3
cfc               call J_virtual_tri_box_pent_isig45(psi,p,1,2,sign, musq,
cfc     1              ze, zu, zero, zero, k, -5, j5zzeu(0,-1,k,is1,1))
cfc               call J_virtual_tri_box_pent_isig45(psi,p,1,2,sign, musq,
cfc     1              zu, ze, zero, zero, k, -5, j5zzue(0,-1,k,is1,1))
cfc               call J_virtual_tri_box_pent_isig45(psi,p,1,2,sign, musq,
cfc     1              zu, ae, zero, zero, k, -5, j5zaue(0,-1,k,is1,1))
             
cfc	       call J_virtual_tri_box_pent_isig45(psi,p,1,2,sign, musq,
cfc     1              ae, zu, zero, zero, k, -5, j5azeu(0,-1,k,is1,1))
          
cfc	     call J_virtual_tri_box_pent_isig45(
cfc     1	            psi,p,1,2,sign, musq,
cfc     1              ze, au, zero, zero, k, -5, j5zaeu(0,-1,k,is1,1))
cfc             call J_virtual_tri_box_pent_isig45(
cfc     1	            psi,p,1,2,sign, musq,
cfc     1  	    ae, au, zero, zero, k, -5, j5aaeu(0,-1,k,is1,1))
cfc            call J_virtual_tri_box_pent_isig45(
cfc     1	            psi,p,1,2,sign, musq,
cfc     1  	    au, ze, zero, zero, k, -5, j5azue(0,-1,k,is1,1))
cfc            call J_virtual_tri_box_pent_isig45(
cfc     1	            psi,p,1,2,sign,musq,
cfc     1  	    au, ae, zero, zero, k, -5, j5aaue(0,-1,k,is1,1))
             
cfc           print*, "pvirtcalc1",  pvirtcalc(:,5)

       call PENT_VBF_CURRENT(!
     $      pvirtcalc(0,1), pvirtcalc(0,3), ! momenta 
     $      pvirtcalc(0,4), pvirtcalc(0,5),  ! momenta
     $      pvirtcalc(0,2),  ! momenta
     $      psi(1,-1,2),psi(1,-1,1),! spinors
     $      ze, au, j5zaeu1(0,-1,k,is1,1),jbzaeu1(0,-1,k,1), ! currents
     $      1, ! Max helicity 
     $      1, k, ! comp, number of current
     $      scale,divmax,   ! Max helicity(1 spinor 2 ext. boson)
     $      mpengauge(1,-1,0,0,1), ! virtuals gauge
     $      mpengaugeborn(1,-1,0,0,1),!born gauge
     $      lwardtestp(wardIndexp)) ! Gauge test

c       print*, "lwardtestp(wardIndexp)",lwardtestp(wardIndexp)
          wardindexp=wardindexp+1

        if(lwardtestp(wardIndexp-1)) then
       call PENT_VBF_CURRENT(!
     $      pvirtcalc(0,1), pvirtcalc(0,3), ! momenta 
     $      pvirtcalc(0,4), pvirtcalc(0,5),  ! momenta
     $      pvirtcalc(0,2),  ! momenta
     $      psi(1,-1,2),psi(1,-1,1),! spinors
     $      ae, au, j5aaeu1(0,-1,k,is1,1),jbaaeu1(0,-1,k,1), ! currents
     $      1, ! Max helicity 
     $      -1, k, ! comp, number of current
     $      scale,divmax,   ! Max helicity(1 spinor 2 ext. boson)
     $      mpengauge(1,-1,0,0,1), ! virtuals gauge
     $      mpengaugeborn(1,-1,0,0,1),!born gauge
     $      lwardtestp(wardIndexp)) ! Gauge test

          wardindexp=wardindexp+1
          else
          wardindexp=wardindexp+1             
             j5aaeu1=0d0
             jbaaeu1=0d0
          endif

cfc
       call PENT_VBF_CURRENT(!
     $      pvirtcalc(0,1), pvirtcalc(0,4), ! momenta 
     $      pvirtcalc(0,3), pvirtcalc(0,5),  ! momenta
     $      pvirtcalc(0,2),  ! momenta
     $      psi(1,-1,2),psi(1,-1,1),! spinors
     $      au,ze, j5azue1(0,-1,k,is1,1),jbazue1(0,-1,k,1), ! currents
     $      1, ! Max helicity 
     $      1, k, ! comp, number of current
     $      scale,divmax,   ! Max helicity(1 spinor 2 ext. boson)
     $      mpengauge(1,-1,0,0,1), ! virtuals gauge
     $      mpengaugeborn(1,-1,0,0,1),!born gauge
     $      lwardtestp(wardIndexp)) ! Gauge test

          wardindexp=wardindexp+1

        if(lwardtestp(wardIndexp-1)) then
       call PENT_VBF_CURRENT(!
     $      pvirtcalc(0,1), pvirtcalc(0,4), ! momenta 
     $      pvirtcalc(0,3), pvirtcalc(0,5),  ! momenta
     $      pvirtcalc(0,2),  ! momenta
     $      psi(1,-1,2),psi(1,-1,1),! spinors
     $      au,ae, j5aaue1(0,-1,k,is1,1),jbaaue1(0,-1,k,1), ! currents
     $      1, ! Max helicity 
     $      -1, k, ! comp, number of current
     $      scale,divmax,   ! Max helicity(1 spinor 2 ext. boson)
     $      mpengauge(1,-1,0,0,1), ! virtuals gauge
     $      mpengaugeborn(1,-1,0,0,1),!born gauge
     $      lwardtestp(wardIndexp)) ! Gauge test

          wardindexp=wardindexp+1

         else
             wardindexp=wardindexp+1
             j5aaue1=0d0
             jbaaue1=0d0
          endif

            enddo
c check for what is used this:
            icount1= icount1+1

      j5zaeu(:,:,:,:,1) = j5zaeu1(:,:,:,:,1)
      j5aaeu(:,:,:,:,1) = j5aaeu1(:,:,:,:,1)
      j5azue(:,:,:,:,1) = j5azue1(:,:,:,:,1)
      j5aaue(:,:,:,:,1) = j5aaue1(:,:,:,:,1)
ccccccc
      IF(debugza) then
         j=1
      do i=1,3
         if(i.lt.3) then
            isig3=is1
          else
             isig3=is1
         endif
c         print*,"k",i
         do isig1=-1,1,2
c            print*, "isig1",isig1
            do mu=0,3

         if (lpent .and. lzs(js1)) then  
               jbzaeu(mu,isig1,i,j)       =    jbzaeu1(mu,isig1,i,j)
               jbaaeu(mu,isig1,i,j)       =    jbaaeu1(mu,isig1,i,j)
               jbazue(mu,isig1,i,j)       =    jbazue1(mu,isig1,i,j)
               jbaaue(mu,isig1,i,j)       =    jbaaue1(mu,isig1,i,j)
         endif

            if((abs(j5zaeu1(mu,isig1,i,isig3,j)/j5zaeu(mu,isig1,i,isig3,j)-1).gt.1d-3)
     $           .and.(abs(j5zaeu(mu,isig1,i,isig3,j)).ne.0d0)) then
            print*, " Z* A -> ll A new result",mu,isig1,isig3,i
            print*,"jvez2",j5zaeu(mu,isig1,i,isig3,j)/j5zaeu1(mu,isig1,i,isig3,j)-1
            print*,"jvez2",j5zaeu(mu,isig1,i,isig3,j)
            print*,"jvez2",j5zaeu1(mu,isig1,i,isig3,j)
            print*," juez2", jbzaeu(mu,isig1,i,j)/jbzaeu1(mu,isig1,i,j)-1
            print*
            print*
            endif

            if((abs(j5aaeu1(mu,isig1,i,isig3,j)/j5aaeu(mu,isig1,i,isig3,j)-1).gt.1d-3)
     $           .and.(abs(j5aaeu(mu,isig1,i,isig3,j)).ne.0d0)) then
            print*, "A* A -> ll A new result",mu,isig1,isig3,i
            print*,"jvez2",j5aaeu(mu,isig1,i,isig3,j)/j5aaeu1(mu,isig1,i,isig3,j)-1
            print*,"jvez2",j5aaeu(mu,isig1,i,isig3,j)
            print*,"jvez2",j5aaeu1(mu,isig1,i,isig3,j)
            print*," juez2", jbaaeu(mu,isig1,i,j)/jbaaeu1(mu,isig1,i,j)-1
            print*
            print*
            endif


            if((abs(j5azue1(mu,isig1,i,isig3,j)/j5azue(mu,isig1,i,isig3,j)-1).gt.1d-3)
     $           .and.(abs(j5azue(mu,isig1,i,isig3,j)).ne.0d0)) then
            print*, " A Z* -> A ll new result",mu,isig1,isig3,i
            print*,"jvez2",j5azue(mu,isig1,i,isig3,j)/j5azue1(mu,isig1,i,isig3,j)-1
            print*,"jvez2",j5azue(mu,isig1,i,isig3,j)
            print*,"jvez2",j5azue1(mu,isig1,i,isig3,j)
            print*," juez2", jbazue(mu,isig1,i,j)/jbazue1(mu,isig1,i,j)-1
            print*
            print*
            endif


            if((abs(j5aaue1(mu,isig1,i,isig3,j)/j5aaue(mu,isig1,i,isig3,j)-1).gt.1d-3)
     $           .and.(abs(j5aaue(mu,isig1,i,isig3,j)).ne.0d0)) then
            print*, " A A* ->  A ll new result",mu,isig1,isig3,i
            print*,"jvez2",j5aaue(mu,isig1,i,isig3,j)/j5aaue1(mu,isig1,i,isig3,j)-1
            print*,"jvez2",j5aaue(mu,isig1,i,isig3,j)
            print*,"jvez2",j5aaue1(mu,isig1,i,isig3,j)
            print*," juez2", jbaaue(mu,isig1,i,j)/jbaaue1(mu,isig1,i,j)-1
            print*
            print*
            endif


         enddo
        enddo
      enddo
      endif
c           stop         
         endif


c for bos = 21: there is no A->vm vm~ contribution:
	if (bos.eq.211.or.bos.eq.212) then
		      j5zzeu = 0d0
		      j5zzue = 0d0
		      j5azeu = 0d0
		      j5zaue = 0d0
	     do k = 1,3	
		do mu = 0,3
		   do isig = -1,1,2
		      j5zzeu(mu,isig,k,is1,1) = 0d0
		      j5zzue(mu,isig,k,is1,1) = 0d0
		      j5azeu(mu,isig,k,is1,1) = 0d0
		      j5zaue(mu,isig,k,is1,1) = 0d0
		   enddo !isig
		enddo !mu
             enddo !k
	 endif    


cfc TO DO: implement gauge test

c         lgc(js1) = .false.
         lzs(js1) = .false.
cfc         if (lpt) then
cfc	 
cfc	    icountmax = 2**(n_points(0)-9)
cfc	    if(icountmax.eq.1) icountmax = 10
cfc
c gauge check (compare contraction of pentagonline with momentum 
c	of attached current to corresponding boxline result 
c	... should give ratio ~1 ):
c            if (sign(1).eq.2 .and. sign(3).eq.1) then
cfc            if (.true.) then
cfc               call gauge_check_isig_5(j5zzeu(0,-1,1,is1,1),p43,
cfc     #	       				         1,psi,p,1,2,zet,zut,0d0)
cfc	       if(bad_gauge_isig) bad_gauge = .true.
cfc	       call gauge_check_isig_5(j5zzue(0,-1,1,is1,1),p43,
cfc     #	       				         1,psi,p,1,2,zut,zet,0d0)
cfc	       if(bad_gauge_isig) bad_gauge = .true.
cfc	      
cfc	       if (bos.eq.22) call gauge_check_isig_5(j5zaeu(0,-1,1,is1,1),
cfc     #	       				      p43,1,psi,p,1,2,zet,aut,0d0)
cfc	       if(bad_gauge_isig) bad_gauge = .true.
cfc	       call gauge_check_isig_5(j5zaue(0,-1,1,is1,1),p43,
cfc     #	       				         1,psi,p,1,2,zut,aet,0d0)
cfc	       if(bad_gauge_isig) bad_gauge = .true.
cfc	      
cfc	       call gauge_check_isig_5(j5azeu(0,-1,1,is1,1),p43,
cfc     #	       				         1,psi,p,1,2,aet,zut,0d0)
cfc	       if(bad_gauge_isig) bad_gauge = .true.
cfc	       if (bos.eq.22) call gauge_check_isig_5(j5azue(0,-1,1,is1,1),
cfc     #	       				      p43,1,psi,p,1,2,aut,zet,0d0)
cfc	       if(bad_gauge_isig) bad_gauge = .true.
cfc	      
cfc	       if (bos.eq.22) call gauge_check_isig_5(j5aaeu(0,-1,1,is1,1),
cfc     #	       				      p43,1,psi,p,1,2,aet,aut,0d0)
cfc	       if(bad_gauge_isig) bad_gauge = .true.
cfc	       if (bos.eq.22) call gauge_check_isig_5(j5aaue(0,-1,1,is1,1),
cfc     #	       				      p43,1,psi,p,1,2,aut,aet,0d0) 
cfc	       if(bad_gauge_isig) bad_gauge = .true.
cfc  	    endif
cfcc	    if (icb1.eq.0) then
cfcc		print*, "criterion chosen: mv5/mvc.gt.1d0"
cfccc		print*,"icountmax=",icountmax
cfcc	    endif !icb1	
cfc	    if(bad_gauge) then
cfc               icb1 = icb1+1
cfc               lgc(js1) = .true.
cfc            endif
cfc        
cfc           endif   ! lpt
	endif


c pentagon corr. to lower line:
cfc       if(lpent) then
cfc         if(.not. lzs(js3)) then
cfc            print*, "PEN2 re-used in qqzaqq.F"
cfc         ENDIF
cfc       endif

       if (lpent .and. lzs(js3)) then     ! need new pentagon graphs for 34 line
c     if(.true.)then
         musq = -p43(4)
         if (lpent) then
c         if(.true.) then
           do mu = 0,3
            pvirtcalc(mu,1) = p(mu,3)
            pvirtcalc(mu,2) = -1d0 * p(mu,4)
            pvirtcalc(mu,3) = -1d0 * qe(mu)
            pvirtcalc(mu,4) = -1d0 * qa(mu)
            pvirtcalc(mu,5) = 1d0 * (p(mu,1) - p(mu,2))
            temp=pvirtcalc(mu,1)+pvirtcalc(mu,2)+pvirtcalc(mu,3)+pvirtcalc(mu,4)+pvirtcalc(mu,5)
cfc            print*,"Mom. conserv.", temp
cfc           print*, "pvirtcalc",  pvirtcalc(mu,5)
         enddo
           scale = scalelowerline

            do k = 1,3
cfc               call J_virtual_tri_box_pent_isig45(psi,p,3,4,sign, musq,
cfc     1              ze, zu, zero, zero, k, -5, j5zzeu(0,-1,k,is3,2))
cfc               call J_virtual_tri_box_pent_isig45(psi,p,3,4,sign, musq,
cfc     1              zu, ze, zero, zero, k, -5, j5zzue(0,-1,k,is3,2))
cfc               call J_virtual_tri_box_pent_isig45(psi,p,3,4,sign, musq,
cfc     1              zu, ae, zero, zero, k, -5, j5zaue(0,-1,k,is3,2))
cfc	    call J_virtual_tri_box_pent_isig45(psi,p,3,4,sign, musq,
cfc     1              ae, zu, zero, zero, k, -5, j5azeu(0,-1,k,is3,2))
               
cfc            call J_virtual_tri_box_pent_isig45(
cfc     1		    psi,p,3,4,sign, musq,
cfc     1  	    ze, au, zero, zero, k, -5, j5zaeu(0,-1,k,is3,2))
cfc           call J_virtual_tri_box_pent_isig45(
cfc     1		    psi,p,3,4,sign, musq,
cfc     1        	    ae, au, zero, zero, k, -5, j5aaeu(0,-1,k,is3,2))
cfc           call J_virtual_tri_box_pent_isig45(
cfc     1		    psi,p,3,4,sign, musq,
cfc     1  	    au, ze, zero, zero, k, -5, j5azue(0,-1,k,is3,2))
cfc           call J_virtual_tri_box_pent_isig45(
cfc     1		    psi,p,3,4,sign, musq,
cfc     1  	    au, ae, zero, zero, k, -5, j5aaue(0,-1,k,is3,2))

       call PENT_VBF_CURRENT(!
     $      pvirtcalc(0,1), pvirtcalc(0,3), ! momenta 
     $      pvirtcalc(0,4), pvirtcalc(0,5),  ! momenta
     $      pvirtcalc(0,2),  ! momenta
     $      psi(1,-1,4),psi(1,-1,3),! spinors
     $      ze, au, j5zaeu1(0,-1,k,is3,2),jbzaeu1(0,-1,k,2), ! currents
     $      1, ! Max helicity 
     $      1, k, ! comp, number of current
     $      scale,divmax,   ! Max helicity(1 spinor 2 ext. boson)
     $      mpengauge(1,-1,0,0,1), ! virtuals gauge
     $      mpengaugeborn(1,-1,0,0,1),!born gauge
     $      lwardtestp(wardIndexp)) ! Gauge test
       wardindexp = wardindexp + 1
        if(lwardtestp(wardIndexp-1)) then
       call PENT_VBF_CURRENT(!
     $      pvirtcalc(0,1), pvirtcalc(0,3), ! momenta 
     $      pvirtcalc(0,4), pvirtcalc(0,5),  ! momenta
     $      pvirtcalc(0,2),  ! momenta
     $      psi(1,-1,4),psi(1,-1,3),! spinors
     $      ae, au, j5aaeu1(0,-1,k,is3,2),jbaaeu1(0,-1,k,2), ! currents
     $      1, ! Max helicity 
     $      -1, k, ! comp, number of current
     $      scale,divmax,   ! Max helicity(1 spinor 2 ext. boson)
     $      mpengauge(1,-1,0,0,1), ! virtuals gauge
     $      mpengaugeborn(1,-1,0,0,1),!born gauge
     $      lwardtestp(wardIndexp)) ! Gauge test
       wardindexp = wardindexp + 1
          else
             wardindexp=wardindexp+1             
             j5aaeu1=0d0
             jbaaeu1=0d0
          endif
       call PENT_VBF_CURRENT(!
     $      pvirtcalc(0,1), pvirtcalc(0,4), ! momenta 
     $      pvirtcalc(0,3), pvirtcalc(0,5),  ! momenta
     $      pvirtcalc(0,2),  ! momenta
     $      psi(1,-1,4),psi(1,-1,3),! spinors
     $      au,ze, j5azue1(0,-1,k,is3,2),jbazue1(0,-1,k,2), ! currents
     $      1, ! Max helicity 
     $      1, k, ! comp, number of current
     $      scale,divmax,   ! Max helicity(1 spinor 2 ext. boson)
     $      mpengauge(1,-1,0,0,1), ! virtuals gauge
     $      mpengaugeborn(1,-1,0,0,1),!born gauge
     $      lwardtestp(wardIndexp)) ! Gauge test
       wardindexp = wardindexp + 1
        if(lwardtestp(wardIndexp-1)) then
       call PENT_VBF_CURRENT(!
     $      pvirtcalc(0,1), pvirtcalc(0,4), ! momenta 
     $      pvirtcalc(0,3), pvirtcalc(0,5),  ! momenta
     $      pvirtcalc(0,2),  ! momenta
     $      psi(1,-1,4),psi(1,-1,3),! spinors
     $      au,ae, j5aaue1(0,-1,k,is3,2),jbaaue1(0,-1,k,2), ! currents
     $      1, ! Max helicity 
     $      -1, k, ! comp, number of current
     $      scale,divmax,   ! Max helicity(1 spinor 2 ext. boson)
     $      mpengauge(1,-1,0,0,1), ! virtuals gauge
     $      mpengaugeborn(1,-1,0,0,1),!born gauge
     $      lwardtestp(wardIndexp)) ! Gauge test
       wardindexp = wardindexp + 1
          else
          wardindexp=wardindexp+1             
             j5aaue1=0d0
             jbaaue1=0d0
          endif
            enddo
            icount2= icount2+1

      j5zaeu(:,:,:,:,2) = j5zaeu1(:,:,:,:,2)
      j5aaeu(:,:,:,:,2) = j5aaeu1(:,:,:,:,2)
      j5azue(:,:,:,:,2) = j5azue1(:,:,:,:,2)
      j5aaue(:,:,:,:,2) = j5aaue1(:,:,:,:,2)
ccccccc
      IF(debugza) then
         j=2
      do i=1,3
         if(i.lt.3) then
            isig3=is3
c            print*, "is3", is3, "is1",is1, "j",j
          else
             isig3=is3
         endif
c         print*,"k",i
         do isig1=-1,1,2
c            print*, "isig1",isig1
            do mu=0,3

         if (lpent .and. lzs(js3)) then  
               jbzaeu(mu,isig1,i,j)       =    jbzaeu1(mu,isig1,i,j)
               jbaaeu(mu,isig1,i,j)       =    jbaaeu1(mu,isig1,i,j)
               jbazue(mu,isig1,i,j)       =    jbazue1(mu,isig1,i,j)
               jbaaue(mu,isig1,i,j)       =    jbaaue1(mu,isig1,i,j)
         endif


            if((abs(j5zaeu(mu,isig1,i,isig3,j)/j5zaeu1(mu,isig1,i,isig3,j)-1).gt.1d-3)
     $           .and.(abs(j5zaeu(mu,isig1,i,isig3,j)).ne.0d0)) then
            print*, " Z* A -> ll A new result",mu,isig1,isig3,i
            print*,"jvez2",j5zaeu(mu,isig1,i,isig3,j)/j5zaeu1(mu,isig1,i,isig3,j)-1
            print*,"jvez2",j5zaeu(mu,isig1,i,isig3,j)
            print*,"jvez2",j5zaeu1(mu,isig1,i,isig3,j)
            print*," juez2", jbzaeu(mu,isig1,i,j)/jbzaeu1(mu,isig1,i,j)-1
            print*
            print*
            endif

            if((abs(j5aaeu(mu,isig1,i,isig3,j)/j5aaeu1(mu,isig1,i,isig3,j)-1).gt.1d-3)
     $           .and.(abs(j5aaeu(mu,isig1,i,isig3,j)).ne.0d0)) then
            print*, "A* A -> ll A new result",mu,isig1,isig3,i
            print*,"jvez2",j5aaeu(mu,isig1,i,isig3,j)/j5aaeu1(mu,isig1,i,isig3,j)-1
            print*,"jvez2",j5aaeu(mu,isig1,i,isig3,j)
            print*,"jvez2",j5aaeu1(mu,isig1,i,isig3,j)
            print*," juez2", jbaaeu(mu,isig1,i,j)/jbaaeu1(mu,isig1,i,j)-1
            print*
            print*
            endif


            if((abs(j5azue(mu,isig1,i,isig3,j)/j5azue1(mu,isig1,i,isig3,j)-1).gt.1d-3)
     $           .and.(abs(j5azue(mu,isig1,i,isig3,j)).ne.0d0)) then
            print*, " A Z* -> A ll new result",mu,isig1,isig3,i
            print*,"jvez2",j5azue(mu,isig1,i,isig3,j)/j5azue1(mu,isig1,i,isig3,j)-1
            print*,"jvez2",j5azue(mu,isig1,i,isig3,j)
            print*,"jvez2",j5azue1(mu,isig1,i,isig3,j)
            print*," juez2", jbazue(mu,isig1,i,j)/jbazue1(mu,isig1,i,j)-1
            print*
            print*
            endif


            if((abs(j5aaue(mu,isig1,i,isig3,j)/j5aaue1(mu,isig1,i,isig3,j)-1).gt.1d-3)
     $           .and.(abs(j5aaue(mu,isig1,i,isig3,j)).ne.0d0)) then
            print*, " A A* ->  A ll new result",mu,isig1,isig3,i
            print*,"jvez2",j5aaue(mu,isig1,i,isig3,j)/j5aaue1(mu,isig1,i,isig3,j)-1
            print*,"jvez2",j5aaue(mu,isig1,i,isig3,j)
            print*,"jvez2",j5aaue1(mu,isig1,i,isig3,j)
            print*," juez2", jbaaue(mu,isig1,i,j)/jbaaue1(mu,isig1,i,j)-1
            print*
            print*
            endif


         enddo
        enddo
      enddo
      endif
c      stop
         endif

c for bos=21: there is no A->vm vm~ contribution:
	  if (bos.eq.211.or.bos.eq.212) then   
		      j5zzeu = 0d0
		      j5zzue = 0d0
		      j5azeu = 0d0
		      j5zaue = 0d0
	     do k = 1,3	
		do mu = 0,3
		   do isig = -1,1,2
		      j5zzeu(mu,isig,k,is3,2) = 0d0
		      j5zzue(mu,isig,k,is3,2) = 0d0
		      j5azeu(mu,isig,k,is3,2) = 0d0
		      j5zaue(mu,isig,k,is3,2) = 0d0
		   enddo !isig
		enddo !mu
             enddo !k
	  endif   

c         lgc(js3) = .false.
         lzs(js3) = .false.
cfc         if (lpt) then
cfcc perform gauge check:
cfcc            if (sign(1).eq.2 .and. sign(3).eq.1) then
cfc            if (.true.) then
cfc               call gauge_check_isig_5(j5zzeu(0,-1,1,is3,2),p21,
cfc     #	       				         1,psi,p,3,4,zet,zut,0d0)
cfc	       if(bad_gauge_isig) bad_gauge = .true.
cfc	       call gauge_check_isig_5(j5zzue(0,-1,1,is3,2),p21,
cfc     #	       				         1,psi,p,3,4,zut,zet,0d0)
cfc	      
cfc	       if(bad_gauge_isig) bad_gauge = .true.
cfc	       if (bos.eq.22) call gauge_check_isig_5(j5zaeu(0,-1,1,is3,2),
cfc     #	       					p21,1,psi,p,3,4,zet,aut,0d0)
cfc	       if(bad_gauge_isig) bad_gauge = .true.
cfc	       call gauge_check_isig_5(j5zaue(0,-1,1,is3,2),p21,
cfc     #	       				         1,psi,p,3,4,zut,aet,0d0)
cfc	       
cfc	       if(bad_gauge_isig) bad_gauge = .true.
cfc	       call gauge_check_isig_5(j5azeu(0,-1,1,is3,2),p21,
cfc     #	       				         1,psi,p,3,4,aet,zut,0d0)
cfc	       if(bad_gauge_isig) bad_gauge = .true.
cfc	       if (bos.eq.22) call gauge_check_isig_5(j5azue(0,-1,1,is3,2),
cfc     #	       					p21,1,psi,p,3,4,aut,zet,0d0)
cfc	       
cfc	       if(bad_gauge_isig) bad_gauge = .true.
cfc	       if (bos.eq.22) call gauge_check_isig_5(j5aaeu(0,-1,1,is3,2),
cfc     #	       					p21,1,psi,p,3,4,aet,aut,0d0)
cfc	       if(bad_gauge_isig) bad_gauge = .true.
cfc	       if (bos.eq.22) call gauge_check_isig_5(j5aaue(0,-1,1,is3,2),
cfc     #	       					p21,1,psi,p,3,4,aut,aet,0d0)
cfc	       if(bad_gauge_isig) bad_gauge = .true.
cfc 
cfc            endif
cfc	    if(bad_gauge) then
cfc               icb2 = icb2+1
cfc               lgc(js3) = .true.
cfc            endif
cfc
cfc         endif ! lpt

c perform gauge test

      endif


#endif

c ------------------------------
c
c compute matrix elements for pentagon on upper line:
c
c need m5e and m5u(pos.j43=1:3,vext1=1:3,vext2=1:3), id of internal boson arbitrary 

      do isig3 = -1,1,2 
	do isig1 = -1,1,2

         call ket2c(psi(1,isig1,1),.true.,p(0,1),isig1,
     1        p43,jqq(0,isig3,2),bkjqq(1,isig1,isig3,1),dummy)
         call bra2c(psi(1,isig1,2),.true.,p(0,2),isig1,
     1        p43,jqq(0,isig3,2),bkjqq(1,isig1,isig3,2),dummy)
	
c  eps1=j43 with v1 arbitrary,eps2=v2,eps3=v3
         m5e(1,1,1) = -s1c(psiau(1,isig1,2),ae,.true.,isig1,bkjqq(1,isig1,isig3,1))
         m5e(1,1,2) = -s1c(psizu(1,isig1,2),ae,.true.,isig1,bkjqq(1,isig1,isig3,1))
         m5e(1,2,1) = -s1c(psiau(1,isig1,2),ze,.true.,isig1,bkjqq(1,isig1,isig3,1))
         m5e(1,2,2) = -s1c(psizu(1,isig1,2),ze,.true.,isig1,bkjqq(1,isig1,isig3,1))

         m5u(1,1,1) = -s1c(psiae(1,isig1,2),au,.true.,isig1,bkjqq(1,isig1,isig3,1))
         m5u(1,1,2) = -s1c(psize(1,isig1,2),au,.true.,isig1,bkjqq(1,isig1,isig3,1))
         m5u(1,2,1) = -s1c(psiae(1,isig1,2),zu,.true.,isig1,bkjqq(1,isig1,isig3,1))
         m5u(1,2,2) = -s1c(psize(1,isig1,2),zu,.true.,isig1,bkjqq(1,isig1,isig3,1))

c  eps1=v1,eps2=j43 with v2 arbitrary,eps3=v3
         m5e(2,1,1) = -s1c(psiau(1,isig1,2),jqq(0,isig3,2),.true.,isig1,psiae(1,isig1,1))
         m5e(2,1,2) = -s1c(psizu(1,isig1,2),jqq(0,isig3,2),.true.,isig1,psiae(1,isig1,1))
         m5e(2,2,1) = -s1c(psiau(1,isig1,2),jqq(0,isig3,2),.true.,isig1,psize(1,isig1,1))
         m5e(2,2,2) = -s1c(psizu(1,isig1,2),jqq(0,isig3,2),.true.,isig1,psize(1,isig1,1))

         m5u(2,1,1) = -s1c(psiae(1,isig1,2),jqq(0,isig3,2),.true.,isig1,psiau(1,isig1,1))
         m5u(2,1,2) = -s1c(psize(1,isig1,2),jqq(0,isig3,2),.true.,isig1,psiau(1,isig1,1))
         m5u(2,2,1) = -s1c(psiae(1,isig1,2),jqq(0,isig3,2),.true.,isig1,psizu(1,isig1,1))
         m5u(2,2,2) = -s1c(psize(1,isig1,2),jqq(0,isig3,2),.true.,isig1,psizu(1,isig1,1))

c  eps1=v1,eps2=v2,eps3=j43 with v3 arbitrary
         m5e(3,1,1) = -s1c(bkjqq(1,isig1,isig3,2),au,.true.,isig1,psiae(1,isig1,1))
         m5e(3,1,2) = -s1c(bkjqq(1,isig1,isig3,2),zu,.true.,isig1,psiae(1,isig1,1))
         m5e(3,2,1) = -s1c(bkjqq(1,isig1,isig3,2),au,.true.,isig1,psize(1,isig1,1))
         m5e(3,2,2) = -s1c(bkjqq(1,isig1,isig3,2),zu,.true.,isig1,psize(1,isig1,1))
         
	 m5u(3,1,1) = -s1c(bkjqq(1,isig1,isig3,2),ae,.true.,isig1,psiau(1,isig1,1))
         m5u(3,1,2) = -s1c(bkjqq(1,isig1,isig3,2),ze,.true.,isig1,psiau(1,isig1,1))
         m5u(3,2,1) = -s1c(bkjqq(1,isig1,isig3,2),ae,.true.,isig1,psizu(1,isig1,1))
         m5u(3,2,2) = -s1c(bkjqq(1,isig1,isig3,2),ze,.true.,isig1,psizu(1,isig1,1))


         do k = 1,6
	 
           do id2 = 1,2
              do id3 = 1,2	
            	 m1p(k,id2,id3) = ( m5e(1,id2,id3) + m5u(1,id2,id3))*
     #    	  	  clr(ifl(2,k),id2,isig1)*clr(ifl(2,k),id3,isig1)   
              enddo
	    enddo

	    do id1 = 1,2
	       do id3 = 1,2 
                  m2p(k,id1,id3) = ( m5e(2,id1,id3) + m5u(2,id1,id3))*
     #    		clr(ifl(1,k),id1,isig1)*clr(ifl(2,k),id3,isig1)
                enddo
	    enddo	  

	    do id1 = 1,2
	       do id2 = 1,2
	          m3p(k,id1,id2) = ( m5e(3,id1,id2) + m5u(3,id1,id2))*
     #    		clr(ifl(1,k),id1,isig1)*clr(ifl(1,k),id2,isig1)      
		enddo
	    enddo	
	     

          enddo !k =1,6
    
c -----------------	  
c virtual stuff:         
c need mv5e and mv5u(pos.j43=1:3,vext1=1:3,vext2=1:3), id of internal boson arbitrary 
	 if (lpent) then
	 
	   do j = 1,3

c  epsj=j43 with vj arbitrary,epsa=va,epsb=vb ( j = 1 or 2 or 3; a,b = A / Z )

	      mv5e(j,1,1) = -dotcc(j5aaeu(0,isig1,j,is1,1),jqq(0,isig3,2))
	      mv5e(j,1,2) = -dotcc(j5azeu(0,isig1,j,is1,1),jqq(0,isig3,2))
	      mv5e(j,2,1) = -dotcc(j5zaeu(0,isig1,j,is1,1),jqq(0,isig3,2))
	      mv5e(j,2,2) = -dotcc(j5zzeu(0,isig1,j,is1,1),jqq(0,isig3,2))

	      mv5u(j,1,1) = -dotcc(j5aaue(0,isig1,j,is1,1),jqq(0,isig3,2))
	      mv5u(j,1,2) = -dotcc(j5azue(0,isig1,j,is1,1),jqq(0,isig3,2))
	      mv5u(j,2,1) = -dotcc(j5zaue(0,isig1,j,is1,1),jqq(0,isig3,2))
	      mv5u(j,2,2) = -dotcc(j5zzue(0,isig1,j,is1,1),jqq(0,isig3,2))	   
	   
	   enddo
        
	
	 do k = 1,6
	 
           do id2 = 1,2
              do id3 = 1,2	
            	 mv1p(k,id2,id3) = ( mv5e(1,id2,id3) + mv5u(1,id2,id3))*
     #    	  	  clr(ifl(2,k),id2,isig1)*clr(ifl(2,k),id3,isig1)   
              enddo
	    enddo

	    do id1 = 1,2
	       do id3 = 1,2 
                  mv2p(k,id1,id3) = ( mv5e(2,id1,id3) + mv5u(2,id1,id3))*
     #    		clr(ifl(1,k),id1,isig1)*clr(ifl(2,k),id3,isig1)
                enddo
	    enddo	  

	    do id1 = 1,2
	       do id2 = 1,2
	          mv3p(k,id1,id2) = ( mv5e(3,id1,id2) + mv5u(3,id1,id2))*
     #    		clr(ifl(1,k),id1,isig1)*clr(ifl(1,k),id2,isig1)      
		enddo
	    enddo	
	     

          enddo !k =1,6
	 
	 endif


c -----------------	  
	  
	  do k = 1,4  
	    
	    mat(k,isig1,isig3,7) = -(m1p(k,1,1)+m1p(k,1,2)+m1p(k,2,1)+m1p(k,2,2))*
     #	        (clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3)*prop43(1)+
     #	         clr(ifl(1,k),2,isig1)*clr(ifl(3,k),2,isig3)*prop43(2))
     #		-
     #			(m2p(k,1,1)+m2p(k,1,2)+m2p(k,2,1)+m2p(k,2,2))*
     #		 (clr(ifl(3,k),1,isig3)*clr(ifl(2,k),1,isig1)*prop43(1)+
     #		  clr(ifl(3,k),2,isig3)*clr(ifl(2,k),2,isig1)*prop43(2))
     #		-
     #			(m3p(k,1,1)+m3p(k,1,2)+m3p(k,2,1)+m3p(k,2,2))*
     #		 (clr(ifl(2,k),1,isig1)*clr(ifl(3,k),1,isig3)*prop43(1)+
     #		  clr(ifl(2,k),2,isig1)*clr(ifl(3,k),2,isig3)*prop43(2))
     
     	   if (lpent) matv(k,isig1,isig3,7) =  
     #	   	       -(mv1p(k,1,1)+mv1p(k,1,2)+mv1p(k,2,1)+mv1p(k,2,2))*
     #	        (clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3)*prop43(1)+
     #	         clr(ifl(1,k),2,isig1)*clr(ifl(3,k),2,isig3)*prop43(2))
     #		-
     #			(mv2p(k,1,1)+mv2p(k,1,2)+mv2p(k,2,1)+mv2p(k,2,2))*
     #		 (clr(ifl(3,k),1,isig3)*clr(ifl(2,k),1,isig1)*prop43(1)+
     #		  clr(ifl(3,k),2,isig3)*clr(ifl(2,k),2,isig1)*prop43(2))
     #		-
     #			(mv3p(k,1,1)+mv3p(k,1,2)+mv3p(k,2,1)+mv3p(k,2,2))*
     #		 (clr(ifl(2,k),1,isig1)*clr(ifl(3,k),1,isig3)*prop43(1)+
     #		  clr(ifl(2,k),2,isig1)*clr(ifl(3,k),2,isig3)*prop43(2))

	  enddo !k
	  
	  if (isig1.eq.-1.and.isig3.eq.-1) then
            fac = clr(3,3,-1)**2*prop43(3)
	     
	    do k = 5,6  	    
	      mat(k,isig1,isig3,7) = 
     #	       -(m1p(k,1,1)+m1p(k,1,2)+m1p(k,2,1)+m1p(k,2,2)+
     #		 m2p(k,1,1)+m2p(k,1,2)+m2p(k,2,1)+m2p(k,2,2)+
     #		 m3p(k,1,1)+m3p(k,1,2)+m3p(k,2,1)+m3p(k,2,2))*fac
	      
	      if (lpent)  matv(k,isig1,isig3,7) = 
     #	       -(mv1p(k,1,1)+mv1p(k,1,2)+mv1p(k,2,1)+mv1p(k,2,2)+
     #		 mv2p(k,1,1)+mv2p(k,1,2)+mv2p(k,2,1)+mv2p(k,2,2)+
     #		 mv3p(k,1,1)+mv3p(k,1,2)+mv3p(k,2,1)+mv3p(k,2,2))*fac
	    enddo !k

	  endif !isig13


	enddo !isig1

      enddo !isig3
      
c -------------------------------------------------------
c
c lower line: mat(k,-1,-1,8) = propt(-1,-1,k,1) with prop21
 
c
c need m5e and m5u(pos.j21=1:3,vext1=1:3,vext2=1:3), vint arbitrary 

      do isig1 = -1,1,2
	do isig3 = -1,1,2
 
         call ket2c(psi(1,isig3,3),.true.,p(0,3),isig3,
     1        p21,jqq(0,isig1,1),bkjqq(1,isig3,isig1,3),dummy)
         call bra2c(psi(1,isig3,4),.true.,p(0,4),isig3,
     1        p21,jqq(0,isig1,1),bkjqq(1,isig3,isig1,4),dummy)

	
c  eps1=j21 with v1 arbitrary,eps2=v2,eps3=v3
         m5e(1,1,1) = -s1c(psiau(1,isig3,4),ae,.true.,isig3,bkjqq(1,isig3,isig1,3))
         m5e(1,1,2) = -s1c(psizu(1,isig3,4),ae,.true.,isig3,bkjqq(1,isig3,isig1,3))
         m5e(1,2,1) = -s1c(psiau(1,isig3,4),ze,.true.,isig3,bkjqq(1,isig3,isig1,3))
         m5e(1,2,2) = -s1c(psizu(1,isig3,4),ze,.true.,isig3,bkjqq(1,isig3,isig1,3))

         m5u(1,1,1) = -s1c(psiae(1,isig3,4),au,.true.,isig3,bkjqq(1,isig3,isig1,3))
         m5u(1,1,2) = -s1c(psize(1,isig3,4),au,.true.,isig3,bkjqq(1,isig3,isig1,3))
         m5u(1,2,1) = -s1c(psiae(1,isig3,4),zu,.true.,isig3,bkjqq(1,isig3,isig1,3))
         m5u(1,2,2) = -s1c(psize(1,isig3,4),zu,.true.,isig3,bkjqq(1,isig3,isig1,3))

c  eps1=v1,eps2=j21 with v2 arbitrary,eps3=v3
         m5e(2,1,1) = -s1c(psiau(1,isig3,4),jqq(0,isig1,1),.true.,isig3,psiae(1,isig3,3))
         m5e(2,1,2) = -s1c(psizu(1,isig3,4),jqq(0,isig1,1),.true.,isig3,psiae(1,isig3,3))
         m5e(2,2,1) = -s1c(psiau(1,isig3,4),jqq(0,isig1,1),.true.,isig3,psize(1,isig3,3))
         m5e(2,2,2) = -s1c(psizu(1,isig3,4),jqq(0,isig1,1),.true.,isig3,psize(1,isig3,3))

         m5u(2,1,1) = -s1c(psiae(1,isig3,4),jqq(0,isig1,1),.true.,isig3,psiau(1,isig3,3))
         m5u(2,1,2) = -s1c(psize(1,isig3,4),jqq(0,isig1,1),.true.,isig3,psiau(1,isig3,3))
         m5u(2,2,1) = -s1c(psiae(1,isig3,4),jqq(0,isig1,1),.true.,isig3,psizu(1,isig3,3))
         m5u(2,2,2) = -s1c(psize(1,isig3,4),jqq(0,isig1,1),.true.,isig3,psizu(1,isig3,3))

c  eps1=v1,eps2=v2,eps3=j21 with v3 arbitrary
         m5e(3,1,1) = -s1c(bkjqq(1,isig3,isig1,4),au,.true.,isig3,psiae(1,isig3,3))
         m5e(3,1,2) = -s1c(bkjqq(1,isig3,isig1,4),zu,.true.,isig3,psiae(1,isig3,3))
         m5e(3,2,1) = -s1c(bkjqq(1,isig3,isig1,4),au,.true.,isig3,psize(1,isig3,3))
         m5e(3,2,2) = -s1c(bkjqq(1,isig3,isig1,4),zu,.true.,isig3,psize(1,isig3,3))
         
	 m5u(3,1,1) = -s1c(bkjqq(1,isig3,isig1,4),ae,.true.,isig3,psiau(1,isig3,3))
         m5u(3,1,2) = -s1c(bkjqq(1,isig3,isig1,4),ze,.true.,isig3,psiau(1,isig3,3))
         m5u(3,2,1) = -s1c(bkjqq(1,isig3,isig1,4),ae,.true.,isig3,psizu(1,isig3,3))
         m5u(3,2,2) = -s1c(bkjqq(1,isig3,isig1,4),ze,.true.,isig3,psizu(1,isig3,3))


         do k = 1,6
	 
           do id2 = 1,2
              do id3 = 1,2	
            	 m1p(k,id2,id3) = ( m5e(1,id2,id3) + m5u(1,id2,id3))*
     #    	  	  clr(ifl(4,k),id2,isig3)*clr(ifl(4,k),id3,isig3)   
              enddo
	    enddo

	    do id1 = 1,2
	       do id3 = 1,2 
                  m2p(k,id1,id3) = ( m5e(2,id1,id3) + m5u(2,id1,id3))*
     #    		clr(ifl(3,k),id1,isig3)*clr(ifl(4,k),id3,isig3)
                enddo
	      enddo	  

	    do id1 = 1,2
	       do id2 = 1,2
	          m3p(k,id1,id2) = ( m5e(3,id1,id2) + m5u(3,id1,id2))*
     #    		clr(ifl(3,k),id1,isig3)*clr(ifl(3,k),id2,isig3)
      
		enddo
	      enddo	  

          enddo
    
c -----------------	  
         
c need mv5e and mv5u(pos.j21=1:3,vext1=1:3,vext2=1:3), id of internal boson arbitrary 
	
	 if (lpent) then
	 
	   do j = 1,3

c  epsj=j21 with vj arbitrary,epsa=va,epsb=vb ( j = 1 or 2 or 3; a,b = A / Z )

	      mv5e(j,1,1) = -dotcc(j5aaeu(0,isig3,j,is3,2),jqq(0,isig1,1))
	      mv5e(j,1,2) = -dotcc(j5azeu(0,isig3,j,is3,2),jqq(0,isig1,1))
	      mv5e(j,2,1) = -dotcc(j5zaeu(0,isig3,j,is3,2),jqq(0,isig1,1))
	      mv5e(j,2,2) = -dotcc(j5zzeu(0,isig3,j,is3,2),jqq(0,isig1,1))

	      mv5u(j,1,1) = -dotcc(j5aaue(0,isig3,j,is3,2),jqq(0,isig1,1))
	      mv5u(j,1,2) = -dotcc(j5azue(0,isig3,j,is3,2),jqq(0,isig1,1))
	      mv5u(j,2,1) = -dotcc(j5zaue(0,isig3,j,is3,2),jqq(0,isig1,1))
	      mv5u(j,2,2) = -dotcc(j5zzue(0,isig3,j,is3,2),jqq(0,isig1,1))    
	   
	   enddo
        
         do k = 1,6
	 
           do id2 = 1,2
              do id3 = 1,2	
            	 mv1p(k,id2,id3) = ( mv5e(1,id2,id3) + mv5u(1,id2,id3))*
     #    	  	  clr(ifl(4,k),id2,isig3)*clr(ifl(4,k),id3,isig3)   
              enddo
	    enddo

	    do id1 = 1,2
	       do id3 = 1,2 
                  mv2p(k,id1,id3) = ( mv5e(2,id1,id3) + mv5u(2,id1,id3))*
     #    		clr(ifl(3,k),id1,isig3)*clr(ifl(4,k),id3,isig3)
                enddo
	      enddo	  

	    do id1 = 1,2
	       do id2 = 1,2
	          mv3p(k,id1,id2) = ( mv5e(3,id1,id2) + mv5u(3,id1,id2))*
     #    		clr(ifl(3,k),id1,isig3)*clr(ifl(3,k),id2,isig3)
      
		enddo
	      enddo	  

          enddo
	  
	  endif !lpent
	
          
	  do k = 1,4  
	    
	    mat(k,isig1,isig3,8) = -(m1p(k,1,1)+m1p(k,1,2)+m1p(k,2,1)+m1p(k,2,2))*
     #	        (clr(ifl(3,k),1,isig3)*clr(ifl(1,k),1,isig1)*prop21(1)+
     #	         clr(ifl(3,k),2,isig3)*clr(ifl(1,k),2,isig1)*prop21(2))
     #		-
     #			(m2p(k,1,1)+m2p(k,1,2)+m2p(k,2,1)+m2p(k,2,2))*
     #		 (clr(ifl(1,k),1,isig1)*clr(ifl(4,k),1,isig3)*prop21(1)+
     #		  clr(ifl(1,k),2,isig1)*clr(ifl(4,k),2,isig3)*prop21(2))
     #		-
     #			(m3p(k,1,1)+m3p(k,1,2)+m3p(k,2,1)+m3p(k,2,2))*
     #		 (clr(ifl(4,k),1,isig3)*clr(ifl(1,k),1,isig1)*prop21(1)+
     #		  clr(ifl(4,k),2,isig3)*clr(ifl(1,k),2,isig1)*prop21(2))
	    
	    if (lpent)  matv(k,isig1,isig3,8) = 
     #	              - (mv1p(k,1,1)+mv1p(k,1,2)+mv1p(k,2,1)+mv1p(k,2,2))*
     #	        (clr(ifl(3,k),1,isig3)*clr(ifl(1,k),1,isig1)*prop21(1)+
     #	         clr(ifl(3,k),2,isig3)*clr(ifl(1,k),2,isig1)*prop21(2))
     #		-
     #			(mv2p(k,1,1)+mv2p(k,1,2)+mv2p(k,2,1)+mv2p(k,2,2))*
     #		 (clr(ifl(1,k),1,isig1)*clr(ifl(4,k),1,isig3)*prop21(1)+
     #		  clr(ifl(1,k),2,isig1)*clr(ifl(4,k),2,isig3)*prop21(2))
     #		-
     #			(mv3p(k,1,1)+mv3p(k,1,2)+mv3p(k,2,1)+mv3p(k,2,2))*
     #		 (clr(ifl(4,k),1,isig3)*clr(ifl(1,k),1,isig1)*prop21(1)+
     #		  clr(ifl(4,k),2,isig3)*clr(ifl(1,k),2,isig1)*prop21(2))

	  enddo !k
	  
	  if (isig1.eq.-1.and.isig3.eq.-1) then
          fac = clr(3,3,-1)**2*prop21(3)
	  
	  do k = 5,6  	    
	    mat(k,isig1,isig3,8) = 
     #	    	-(m1p(k,1,1)+m1p(k,1,2)+m1p(k,2,1)+m1p(k,2,2)+
     #		  m2p(k,1,1)+m2p(k,1,2)+m2p(k,2,1)+m2p(k,2,2)+
     #		  m3p(k,1,1)+m3p(k,1,2)+m3p(k,2,1)+m3p(k,2,2))*fac
	    if (lpent) matv(k,isig1,isig3,8) = 
     #	    	-(mv1p(k,1,1)+mv1p(k,1,2)+mv1p(k,2,1)+mv1p(k,2,2)+
     #		  mv2p(k,1,1)+mv2p(k,1,2)+mv2p(k,2,1)+mv2p(k,2,2)+
     #		  mv3p(k,1,1)+mv3p(k,1,2)+mv3p(k,2,1)+mv3p(k,2,2))*fac
	  enddo !k
	  
	  endif !isig13

	enddo !isig1

      enddo !isig3
c
c ------------------------------------------------------------------------
c
     
c sum the graphs, square them and map them onto uucc, uuss etc.

c i = 1		VV
c i = 2,3	BV
c i = 4,5	AZZ/ZZZ
c i = 6,9	BB
c i = 7,8	P


      if(lvirtratio) then!.and.lward)then
         do div = divmax,divmax
            do i= 2,9
               if (div.eq.1) then
                  K_contcheck = -3
               else
                  K_contcheck = -2
               endif

               print*,"ratio check for div = ",div
               do k = 1,4
                  do isig1 = -1,1,2
                     do isig3 = -1,1,2
                        if(i.eq.9) then
                           if( mat(k,isig1,isig3,6).eq.0d0 .and.
     $                          matv(k,isig1,isig3,i).ne.0d0) then
                              print*,"matv != 0, mat = 0"
                              print*,"k,isig1,isig3,i,matv",k,isig1,isig3,i,matv(k,isig1,isig3,i)
                           endif                        
                           if(mat(k,isig1,isig3,6).ne.0d0) then
                              print*,"k, i, sig1, sig3, mborn",k,i,isig1,isig3,mat(k,isig1,isig3,6)
                              print*,"mvirt",matv(k,isig1,isig3,i)
                              ratio=(matv(k,isig1,isig3,i))/(mat(k,isig1,isig3,6)*K_contcheck) 
                              print*,"ratio",ratio
                           endif
                        else
                           if( mat(k,isig1,isig3,i).eq.0d0 .and.
     $                          matv(k,isig1,isig3,i).ne.0d0) then
                              print*,"matv != 0, mat = 0"
                              print*,"k,isig1,isig3,i,matv",k,isig1,isig3,i,matv(k,isig1,isig3,i)
                           endif
                           if(mat(k,isig1,isig3,i).ne.0d0) then
                              print*,"k, i, sig1, sig3, mborn",k,i,isig1,isig3,mat(k,isig1,isig3,i)
                              print*,"mvirt",matv(k,isig1,isig3,i)
                              ratio = matv(k,isig1,isig3,i)/(mat(k,isig1,isig3,i)*K_contcheck) 
                              print*,"ratio",ratio

                           endif
                        endif
                     enddo
                  enddo
               enddo
               print*
               write(*,*)
               read(*,*)
            enddo               !i
         enddo
      endif


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
	
   
      if (lpt) then
         xgc1 = real(icount1+1)/real(icount1-icb1*0.5+1) ! factor 1/2 for 2 helicities
         xgc2 = real(icount2+1)/real(icount2-icb2*0.5+1)
      endif
		      

      do k = 1,6
         res(k) = 0
         resv(k) = 0
         tree(k) = 0
         do isig1 = -1,1,2
            do isig3 = -1,1,2
               mm(k,isig1,isig3) = 0
               do i = 1,8
                  mm(k,isig1,isig3) = mm(k,isig1,isig3) + 
     1                                mat(k,isig1,isig3,i)
cjp
c                  print*, "mat(",k,",",isig1,",",isig3,",",i,")",
c     1             mat(k,isig1,isig3,i), "mm", mm(k,isig1,isig3) 
               enddo
               res(k) = res(k) + dreal(mm(k,isig1,isig3))**2
     &                         + dimag(mm(k,isig1,isig3))**2
c               print*, "res", res(k)
                if (lnlo) then
                  mv12(k,isig1,isig3) = 
     1         +  matv(k,isig1,isig3,2) + matv(k,isig1,isig3,4) + 
     2            matv(k,isig1,isig3,6) + matv(k,isig1,isig3,7) 
                  mv34(k,isig1,isig3) =
     1          + matv(k,isig1,isig3,3) + matv(k,isig1,isig3,5) +
     2            matv(k,isig1,isig3,9) + matv(k,isig1,isig3,8)


                  if (lpt) then   ! Pentagon are considered for "good PS points" only
cfc                     if (lgc(js1)) then
cfc                        mv12(k,isig1,isig3)=0
cfc                     else
                        mv12(k,isig1,isig3)=mv12(k,isig1,isig3)*xgc1
cfc                     endif
                  endif
                  if (lpt) then   ! Pentagon contributes are considered for "good PS points" only
cfc                     if (lgc(js3)) then
cfc                        mv34(k,isig1,isig3)=0
cfc                     else
                        mv34(k,isig1,isig3)=mv34(k,isig1,isig3)*xgc2
cfc                     endif
                  endif

c  add Born type term and multiply by F_q = alphas*C_2/4pi
c  the factor pi^2/3+9/2 for the born term is after adding the subtraction term
c  and the counter term for the renormalization of the pdfs
                  if (nlo.gt.0) then
                     mv12(k,isig1,isig3) = als(1,1)*c2o4pi*
     1                ( mv12(k,isig1,isig3) + mm(k,isig1,isig3)*cvirt )
                     mv34(k,isig1,isig3) = als(2,1)*c2o4pi*
     1                ( mv34(k,isig1,isig3) + mm(k,isig1,isig3)*cvirt )
                  else
                     mv12(k,isig1,isig3) = 
     1                    als(1,1)*c2o4pi*mv12(k,isig1,isig3)
                     mv34(k,isig1,isig3) = 
     1                    als(2,1)*c2o4pi*mv34(k,isig1,isig3)
                  endif
                  resv(k) = resv(k) + 2*dreal(
     1                 mm(k,isig1,isig3)   *
     1              conjg( mv12(k,isig1,isig3)+mv34(k,isig1,isig3) )  )
               endif

           enddo
         enddo  
	       
         tree(k) = res(k) * 9D0 
         if (nlo.ge.0) then
            res(k) = (res(k)+resv(k))*9d0      ! 9 is the color sum factor
         else
            res(k) = resv(k)*9d0               ! 9 is the color sum factor
         endif
      enddo
 
c eliminate processes with photon virtuality below cutoff
      if ( abs(p21(4)).lt.qsqAmin .or.  
     1     abs(p43(4)).lt.qsqAmin ) then
         do k = 1,4
            res(k) = 0
            tree(k) = 0
         enddo
      endif

c TO SET AMPLITUDES TOZEO

c----------------------
      if(lbox) then
         do i = 1,wardIndex-1
            lwardBox = lwardBox.and.lwardtest(i)
         enddo                  ! i
         if(.not.lwardBox) then
             res=0d0
             jez1  = 0d0
             jea1  = 0d0
             jua1  = 0d0
             jvez1 = 0d0
             jvea1 = 0d0
             jvua1 = 0d0
             ja1   = 0d0
             jva1  = 0d0
             jz1   = 0d0
             jvz1  = 0d0
             juz=0d0
             jvuz=0d0
c            print*, "BOX_GAUGE", lwardBox
            boxcountm = boxcountm+1
            boxcountps(ps_number) = boxcountps(ps_number)+1
         else
         boxcount2m = boxcount2m+1
         boxcount2ps(ps_number) = boxcount2ps(ps_number)+1
         endif
c         pause
      endif	


      if(lpent) then
         do i = 1,wardIndexp-1
            lwardpent = lwardpent.and.lwardtestp(i)
         enddo                  ! i
         if(.not.lwardPent.or.res(1).eq.0d0) then
c            print*, "PENT_GAUGE", lwardpent
            res=0d0
            jbzaeu1 = 0d0 
            jbaaeu1 = 0d0 
            jbazue1 = 0d0 
            jbaaue1 = 0d0 
            j5zaeu1 = 0d0 
            j5aaeu1 = 0d0 
            j5azue1 = 0d0 
            j5aaue1 = 0d0 
            jbzaeu1 = 0d0 
            jbaaeu1 = 0d0 
            jbazue1 = 0d0 
            jbaaue1 = 0d0 
	    j5zzeu  = 0d0
	    j5zzue  = 0d0
	    j5azeu  = 0d0
	    j5zaue  = 0d0

            pentcountm = pentcountm+1
            pentcountps(ps_number) = pentcountps(ps_number)+1
         else
         pentcount2m = pentcount2m+1
         pentcount2ps(ps_number) = pentcount2ps(ps_number)+1
         endif
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


