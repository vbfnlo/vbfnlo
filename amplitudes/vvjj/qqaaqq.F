c---------------------------------------------------------------------------
      
      subroutine qqaaqq(pbar,sign, nlo, L,bos, ps_number,
     1                  uucc,uuss,ddcc,ddss,udsc,ducs)
          use globalvars, only: ldoblha
      implicit none

c This code is a simplification of the qqzaqq.F by Francisco Campanario.
c Francisco.Campanario@ific.uv.es.
c That in its turn is a simplification of the qqzzqq.F by
c Barbara Jaeger, baj@particle.physik.uni-karlsruhe.de.
c
c Ivan Rosario, Ivan.Rosario@ific.uv.es
c Last modified: 28/08/2020.
c
c qqaaqq calculates the matrix elements**2 for electroweak
c weak boson pair production by quark quark scattering
c
c       q1 q3    ---->   q2 q4 A A. 
c
c Only the t-channel exchange of elctroweak bosons is considered. 
c s-channel production of 3 weak bosons is NOT implemented.
c
c This code is modified to allow for virtual corrections, more precisely
c the interference of Born with the finite part of virtual diagrams
c for 
c
c INPUT:  NLO = 1       return uucc = |M_born|^2 + 2Re(M_Born^* M_virt)
c         NLO = 0       return uucc = |M_born|^2   etc.
c  for steering the calculation of finite box and pentagon contributions
c  the following additional options are implemented
c         NLO = +4      set all finite box and pentagon contributions to 0
c                       but include cvirt*M_Born contribution
c         NLO = +5      true finite pentagons = 0, all other virtual contributions
c                 #included, including leading q^mu/m_V terms of 
c                       pentagons which can be related to boxes
c       remaining cases:    uucc = 2Re(M_Born^* M_virt), more precisely
c       the piece which is missing for the full NLO=1 when calling +|NLO|
c         NLO = -4      finite Born + q^mu/m_V terms of Pentagon contributions only, no cvirt
c         NLO = -5      true finite Pentagon contributions only
c
c integer L (=1,2) is used for the subtraction code to select leptonic 
c tensors for gluon emission off the upper Line (L=1) or lower quark 
c line (L=2)

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
      double precision pbar(0:3,4+max_v), musq
      double precision uucc,uuss,ddcc,ddss,udsc,ducs,res(6),resv(6)
      double precision  tree(6)
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
     1               m1p(6,2,2),m2p(6,2,2),m3p(6,2,2)
      double complex mv5e(3,2,2),mv5u(3,2,2),
     1               mv1p(6,2,2),mv2p(6,2,2),mv3p(6,2,2)
      double complex ma(2), mz(2),
     1               m1ae,m1ze,m2ae,m2ze,z1ae,z1ze,z2ae,z2ze,
     1               m1au,m1zu,m2au,m2zu,z1au,z1zu,z2au,z2zu,
     1               m3au,m3zu,m4au,m4zu,m3ae,m3ze,m4ae,m4ze,
     1               z3ae,z3ze,z3zu,z3au,z4ae,z4ze,z4zu,z4au,
     1               mate12,matu12,mate34,matu34

      double complex mv1ae,mv1ze,mv2ae,mv2ze,mv1au,mv1zu,mv2au,mv2zu,
     1               mv3au,mv3zu,mv4au,mv4zu,mv3ae,mv3ze,mv4ae,mv4ze
      double complex matve12,matvu12,matve34,matvu34
      double complex mv1aae,mv1aze,mv1zae,mv1zze,
     1               mv2aae,mv2aze,mv2zae,mv2zze,mv1e,mv2e,
     1               mv1aau,mv1azu,mv1zau,mv1zzu,
     1               mv2aau,mv2azu,mv2zau,mv2zzu,mv1u,mv2u,
     1               mv3aae,mv3aze,mv3zae,mv3zze,
     1               mv4aae,mv4aze,mv4zae,mv4zze,mv3e,mv4e,
     1               mv3aau,mv3azu,mv3zau,mv3zzu,
     1               mv4aau,mv4azu,mv4zau,mv4zzu,mv3u,mv4u
      double complex m1aae,m1aze,m1zae,m1zze,z1aze,z1zze,
     1               m2aae,m2aze,m2zae,m2zze,m1e,m2e,
     1               m1aau,m1azu,m1zau,m1zzu,z1azu,z1zzu,
     1               m2aau,m2azu,m2zau,m2zzu,m1u,m2u,
     1               mz1aae,mz1aze,mz1zae,mz1zze,
     1               mz2aae,mz2aze,mz2zae,mz2zze,
     1               mz1aau,mz1azu,mz1zau,mz1zzu,
     1               mz2aau,mz2azu,mz2zau,mz2zzu,
     1               m3aae,m3aze,m3zae,m3zze,z3aze,z3zze,
     1               m4aae,m4aze,m4zae,m4zze,m3e,m4e,
     1               m3aau,m3azu,m3zau,m3zzu,z3azu,z3zzu,
     1               m4aau,m4azu,m4zau,m4zzu,m3u,m4u,
     1               mz3aae,mz3aze,mz3zae,mz3zze,
     1               mz4aae,mz4aze,mz4zae,mz4zze,
     1               mz3aau,mz3azu,mz3zau,mz3zzu,
     1               mz4aau,mz4azu,mz4zau,mz4zzu
      double complex mezz(4,-1:1,-1:1),meza(4,-1:1,-1:1),
     1               meaz(4,-1:1,-1:1),meaa(4,-1:1,-1:1),
     1               muzz(4,-1:1,-1:1),muza(4,-1:1,-1:1),
     1               muaz(4,-1:1,-1:1),muaa(4,-1:1,-1:1),
     1               mvezz(4,2,-1:1,-1:1),mveza(4,2,-1:1,-1:1),
     1               mveaz(4,2,-1:1,-1:1),mveaa(4,2,-1:1,-1:1),
     1               mvuzz(4,2,-1:1,-1:1),mvuza(4,2,-1:1,-1:1),
     1               mvuaz(4,2,-1:1,-1:1),mvuaa(4,2,-1:1,-1:1),
     1               zezz(-1:1,-1:1,3),zeza(-1:1,-1:1,3),
     1               zeaz(-1:1,-1:1,3),zeaa(-1:1,-1:1,3),
     1               zuzz(-1:1,-1:1,3),zuza(-1:1,-1:1,3),
     1               zuaz(-1:1,-1:1,3),zuaa(-1:1,-1:1,3)
      double complex mva(2),mvz(2)
      double complex psi(2,-1:1,4), jqq(0:5,-1:1,2) 
      double complex psize(2,-1:1,4),psiae(2,-1:1,4),
     1               psizu(2,-1:1,4),psiau(2,-1:1,4),
     1               jez(0:5,-1:1,4),jea(0:5,-1:1,4),
     1               juz(0:5,-1:1,4),jua(0:5,-1:1,4),
     1               jvez(0:5,-1:1,-1:1,4), jvea(0:5,-1:1,-1:1,4),
     1               jvuz(0:5,-1:1,-1:1,4), jvua(0:5,-1:1,-1:1,4)
      double complex bkjqq(2,-1:1,-1:1,4)
      double complex ja(0:5,-1:1,-1:1,4), jz(0:5,-1:1,-1:1,4), 
     1               jva(0:5,-1:1,-1:1,4), jvz(0:5,-1:1,-1:1,4),
     1               psia(2,-1:1,-1:1,4), psiz(2,-1:1,-1:1,4)
      double complex j5zzue(0:3,-1:1,3,-1:1,2),
     1           j5zzeu(0:3,-1:1,3,-1:1,2),
     1           j5zaue(0:3,-1:1,3,-1:1,2), j5zaeu(0:3,-1:1,3,-1:1,2),
     1           j5azue(0:3,-1:1,3,-1:1,2), j5azeu(0:3,-1:1,3,-1:1,2),
     1           j5aaue(0:3,-1:1,3,-1:1,2), j5aaeu(0:3,-1:1,3,-1:1,2)
      double complex epsee(0:3), epsuu(0:3),
     1           epsCCuu(0:5,-1:1,2),epsCCee(0:5,-1:1,2),
     1           epsee6(0:3), epsuu6(0:3),
     1           epsCCuu6(0:5,-1:1,2),epsCCee6(0:5,-1:1,2)
      double complex epsaa(0:3),epsaz(0:3),epsza(0:3),epszz(0:3),
     1           epsaee(0:5,-1:1,3:4,2), epszee(0:5,-1:1,3:4,2),
     1           epsauu(0:5,-1:1,3:4,2), epszuu(0:5,-1:1,3:4,2),
     1           qepszee(3:4),qepszuu(3:4) 
      double precision fqze(0:4,4),fqae(0:4,4),fqzu(0:4,4),fqau(0:4,4)
      double precision fq(0:4,4),qee(0:4), quu(0:4), dummy(0:4)
      double complex zm2i(2:4)
      double complex jj21ae(-1:1),jj21au(-1:1),jj43ae(-1:1),
     1            jj43au(-1:1),
     1            jj21ze(-1:1),jj21zu(-1:1),jj43ze(-1:1),jj43zu(-1:1)
      double complex zuu(4:5),zee(4:5)
      double complex qepsCCee(2), qepsCCuu(2),qepsCCee6(2), qepsCCuu6(2)
      double complex propt1(-1:1,-1:1,4,2), propt2(-1:1,-1:1,4,2),
     1           propt(-1:1,-1:1,5:6,2), prop_ee(3), prop_uu(3),
     1           propbbe, propbbu,
     1           propbbea(6,-1:1,-1:1), propbbua(6,-1:1,-1:1),
     1           propbbez(6,-1:1,-1:1), propbbuz(6,-1:1,-1:1)
      double complex fac

      double complex contract_Tjj, dotcc, dotrc, dotqj, s1c, 
     1               contract_Tjj_ck
      external contract_Tjj, dotcc, dotrc, dotqj, s1c,contract_Tjj_ck
      logical ldebug, ldebugv, ldebugvf,linit,lerror,lgc(4)
      data linit /.true./, lerror /.false./, ldebug /.false./
      data lgc /4*.false./
      
      logical vvdebug,vbdebug, vlldebug,bbdebug,pdebug
      parameter (vvdebug = .false.,vbdebug = .false.,vlldebug =.false.,
     1           bbdebug = .false.,pdebug = .false.)
      
      integer*8 icb1, icount1, icb2, icount2
      data icb1/0/,icount1/0/,icb2/0/,icount2/0/,Lold/1/
      double precision xgc1, xgc2

      parameter (ldebugv=.false., ldebugvf=.false.)
      save ifl, Lold, zm2i, linit, lgc, icb1, icb2, icount1, icount2
      save ja,jz,jva,jvz,psia,psiz,fq 
      save j5zzue, j5zzeu,j5zaue, j5zaeu,j5azue, j5azeu,j5aaue, j5aaeu
      double complex  zero
      parameter (zero = (0d0,0d0) )
      integer ps_number

c
c variables for virtual corrections
c
      logical bad_gauge,bad_gauge_sin,bad_gauge_isig
      common /vbsgauge / bad_gauge,bad_gauge_sin,bad_gauge_isig
      double precision c2,c2o4pi, cvirtc
      parameter (c2=4d0/3d0, c2o4pi=c2/4d0/pi)
      double precision lnmusqp1p2, lnmusqp3p4
      logical lnlo, lbox, lpent, lpt

      double precision  pvirtcalc(0:3,5)
      double precision scalelowerline,scale, scaleupperline,temp

      DOUBLE COMPLEX mboxgauge(2,-1:1,0:3,0:2,12),
     1   mboxgaugeborn(2,-1:1,0:3,0:2,12)
      DOUBLE COMPLEX mpengauge(3,-1:1,0:3,0:2,12), 
     1   mpengaugeborn(3,-1:1,0:3,0:2,12)

      logical lwardtest(100),lwardtestp(100)
      integer wardIndex,div,divMax,divStep
      double complex jez1(0:5,-1:1,4),jea1(0:5,-1:1,4),
     1               jua1(0:5,-1:1,4),
     1               jvez1(0:5,-1:1,-1:1,4), jvea1(0:5,-1:1,-1:1,4),
     1               jvua1(0:5,-1:1,-1:1,4)
      double complex jbzaeu1(0:3,-1:1,3,2), j5zaeu1(0:3,-1:1,3,-1:1,2),
     1               jbzaeu(0:3,-1:1,3,2)
      double complex jbaaeu1(0:3,-1:1,3,2), j5aaeu1(0:3,-1:1,3,-1:1,2),
     1               jbaaeu(0:3,-1:1,3,2)
      double complex jbazue1(0:3,-1:1,3,2), j5azue1(0:3,-1:1,3,-1:1,2),
     1               jbazue(0:3,-1:1,3,2)
      double complex jbaaue1(0:3,-1:1,3,2), j5aaue1(0:3,-1:1,3,-1:1,2),
     1               jbaaue(0:3,-1:1,3,2)
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

      wardindex = 1
      wardindexp = 1
      lwardBox = .true.
      lwardPent = .true.
      divMax = 0
      divStep = 2

      lnlo = NLO.ne.0                                  ! include some virtual stuff if T
      lbox = NLO.eq.1 .or. NLO.eq.-4 .or. NLO.eq.5     ! call boxline if T
      lpt  = NLO.eq.1 .or. NLO.le.-5                   ! true pentagon contributions on
!      lpq  = NLO.eq.5 .or. NLO.eq.-4                  ! q^mu/m_V terms of Pentagon contributions only
      lpent= lpt !.or. lpq                             ! include pentagons if T

      bad_gauge = .false. ! set F at beginning of each run

      if(lvirtratio) then
        divMax = 2
        divStep = 1
        lbox = .true.
        lpent = .true.
        cvirtc = 0d0
        call InitGaugeTest
      endif

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
         k = 6                  ! ducs
         ifl(1,k) = 4
         ifl(2,k) = 3
         ifl(3,k) = 3
         ifl(4,k) = 4
         zm2i(2) = 1/dcmplx(xm2(2),-xmg(2))
         zm2i(3) = 1/dcmplx(xm2(3),-xmg(3))
         zm2i(4) = 1/dcmplx(xm2(4),-xmg(4))
      endif
  
c  Initialization of the matrix elements to zero

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

c Initialization of the currents to zero

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

        p21(mu) = p(mu,2) - p(mu,1)
        p43(mu) = p(mu,4) - p(mu,3) 
        pez(mu) = p(mu,5)               ! Momenta of one of the photons.
        puz(mu) = p(mu,6)               ! Momenta of the other photon. 
        pza(mu) = pez(mu) + puz(mu)
      enddo

      p21(4) = p21(0)**2 - p21(1)**2 - p21(2)**2 - p21(3)**2
      p43(4) = p43(0)**2 - p43(1)**2 - p43(2)**2 - p43(3)**2
      pez(4) = pez(0)**2 - pez(1)**2 - pez(2)**2 - pez(3)**2
      puz(4) = puz(0)**2 - puz(1)**2 - puz(2)**2 - puz(3)**2
      pza(4) = pza(0)**2 - pza(1)**2 - pza(2)**2 - pza(3)**2

      scaleupperline = -1d0* p21(4)
      scalelowerline = -1d0* p43(4)

      if(lbox.or.lpent) then
        scaleupperline = -1d0* p21(4)
        scalelowerline = -1d0* p43(4)
      endif

c
c     Momentum conservation test.
c

      if (lwarn) then
         do mu=0,3
           dummy(mu) = p(mu,1)+p(mu,3)-p(mu,5)- 
     1                 (p(mu,2)+p(mu,4)+p(mu,6))
         enddo
         if (abs(dummy(0)).gt.1d-8 .or.
     1       abs(dummy(1)).gt.1d-8 .or.
     2       abs(dummy(2)).gt.1d-8 .or.
     3       abs(dummy(3)).gt.1d-8 ) then
            print*," momentum mismatch in qqaaqq "
            print 102," dummy = ",(dummy(mu),mu=0,3)
            print 103," sign  = ",sign
 101        format( " p(", i1, ") = ", 4(f10.3, 2x) )
 102        format(a,4f10.3)
 103        format(a,6i5)
            do i = 1,6
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
c as seen from upper line these W momenta are INCOMING. They are OUTGOING as seen from lower 
c line

      do mu = 0,3
        qee(mu) = pez(mu) + p(mu,2) - p(mu,1)    ! photon emitted on upper line
        quu(mu) = puz(mu) + p(mu,2) - p(mu,1)    ! photon emitted on upper line
      enddo

      qee(4) = qee(0)**2-qee(1)**2-qee(2)**2-qee(3)**2
      quu(4) = quu(0)**2-quu(1)**2-quu(2)**2-quu(3)**2

      prop_ee(1) = 1d0/qee(4)
      prop_ee(2) = 1d0/dcmplx(qee(4)-xm2(2),xmg(2))
      prop_ee(3) = 1d0/dcmplx(qee(4)-xm2(3),xmg(3))
      prop_uu(1) = 1d0/quu(4)
      prop_uu(2) = 1d0/dcmplx(quu(4)-xm2(2),xmg(2))
      prop_uu(3) = 1d0/dcmplx(quu(4)-xm2(3),xmg(3))

      ! Creates the wavefunctions of the external quarks
      call psi0m(4,pbar(0,1),sign(1),psi)

      ! Creates the quarkline currents from the external quark
      ! wavefunctions and momenta.
      call curr6(1,psi(1,-1,2),p(0,2),psi(1,-1,1),p(0,1),jqq(0,-1,1))
      call curr6(1,psi(1,-1,4),p(0,4),psi(1,-1,3),p(0,3),jqq(0,-1,2))
      
      ! Contracts the quarkline currents with the neutral currents
      ! leptonic tensors with two final state photons.
      ! The names of the leptonic tensors are the same of the 
      ! qq -> zzqq. But they don't describe our case, here we give
      ! the correspondence with the leptonic tensor name and process.
      ! aazz: a a -> a a
      ! azzz: a z -> a a
      ! zazz: z a -> a a
      ! zzzz: z z -> a a       
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
      
      ! Contracts the quarkline currents with the charged currents
      ! leptonic tensors. 
      ! wwzz5: wm wp -> a a
      ! wwzz6: wp wm -> a a
      mzz6 = contract_Tjj(wwzz6(0,0,L),jqq(0,-1,2),jqq(0,-1,1))
      mzz5 = contract_Tjj(wwzz5(0,0,L),jqq(0,-1,1),jqq(0,-1,2))
      mat(5,-1,-1,1) = mzz5*clr(3,3,-1)**2
      mat(6,-1,-1,1) = mzz6*clr(3,3,-1)**2

cfc  Virtuals ZERO. They factorize to the born
! TODO POLES
      matv(5,-1,-1,1) = (0d0,0d0)
      matv(6,-1,-1,1) = (0d0,0d0)

c------------------------------------------------------------------------
c prepare box diagrams: attach A to external spinors
c 
c     isig = +-1   : left- and righthanded spinors coupling to A
c
c Notation for virtual 2-component spinors and momenta
c
c A attached to quark number i: psiae(*,isig,i) with momentum fqae(mu,i)
c A attached to quark number i: psiau(*,isig,i) with momentum fqau(mu,i)
c 
c the fermion current corresponding to a quark line with the real 
c emitted V-> l+l-  attached next to quark number i is stored in 
c jlv(mu,isig,i) etc. 
c 
c For the virtual amlitudes the notation, e.g.     jvlv(mu,isig,is,i)
c is used for the boxline correction to a quark line with one V attached 
c next to
c quark #i and a free Lorentz index mu for the second attached EW boson. 
c is=+-1 refers to the sign factor of this quark (vs.antiquark line). 
c They are recalculated only if this quark line sign has not been 
c calculated yet for this phase space point (i.e. lzs = .true.). 
c Otherwise they are taken from saved previous calculation
c------------------------------------------------------------------------
      do i = 1,3,2
        do isig = -1,1,2
      
          call ket2c(psi(1,isig,i),.true.,p(0,i),isig,qe,ae,
     1               psiae(1,isig,i),fqae(0,i))
          call ket2c(psi(1,isig,i),.true.,p(0,i),isig,qu,au,
     1               psiau(1,isig,i),fqau(0,i))
         

          call bra2c(psi(1,isig,i+1),.true.,p(0,i+1),isig,qe,ae,
     1               psiae(1,isig,i+1),fqae(0,i+1))
          call bra2c(psi(1,isig,i+1),.true.,p(0,i+1),isig,qu,au,
     1               psiau(1,isig,i+1),fqau(0,i+1))

        enddo !isig
        
          call curr6(1,psi(1,-1,i+1),p(0,i+1),
     1                  psiae(1,-1,i),fqae(0,i), jea(0,-1,i))
          call curr6(1,psi(1,-1,i+1),p(0,i+1),
     1                  psiau(1,-1,i),fqau(0,i), jua(0,-1,i))

          call curr6(1,psiae(1,-1,i+1),fqae(0,i+1),
     1                  psi(1,-1,i),p(0,i), jea(0,-1,i+1))
          call curr6(1,psiau(1,-1,i+1),fqau(0,i+1),
     1                  psi(1,-1,i),p(0,i), jua(0,-1,i+1))
   
      enddo !i

#ifdef WITH_NLO  
       if (lbox) then
         if(.not. lzs(js1)) then
            print*, "STOP in qqaaqq.F"
         ENDIF  
      endif

      if (lbox .and. lzs(js1)) then

        do i = 1,2 ! Two permutations of momenta

          ! Box correction to the upper line for photon 1.
          do mu = 0,3
            pvirtcalc(mu,1) = p(mu,1)
            pvirtcalc(mu,2) = -1d0 * p(mu,2)
            pvirtcalc(mu,3) = -1d0 * qe(mu)
            pvirtcalc(mu,4) = 1d0 * (p(mu,3) - p(mu,4) - qa(mu))
            temp = pvirtcalc(mu,1) + pvirtcalc(mu,2)
     1             + pvirtcalc(mu,3) + pvirtcalc(mu,4)
          enddo

          scale = scaleupperline
          call BOX_VBF_CURR(!
     $      pvirtcalc(0,1), pvirtcalc(0,3),             ! momenta 
     $      pvirtcalc(0,4), pvirtcalc(0,2),             ! momenta
     $      psi(1,-1,2),psi(1,-1,1),                    ! spinors
     $      ae, jvea1(0,-1,is1,i),jea1(0,-1,i),         ! currents
     $      1,                                          ! Max helicity 
     $      1, 3-i,                                     ! comp, number of current
     $      scale, divmax,                              ! Max helicity(1 spinor 2 ext. boson)
     $      mboxgauge(1,-1,0,0,2),                      ! virtuals gauge
     $      mboxgaugeborn(1,-1,0,0,2),                  !born gauge
     $      lwardtest(wardIndex))                       ! Gauge test
    
          wardindex = wardindex + 1

          do mu = 0,3
            pvirtcalc(mu,1) = p(mu,1)
            pvirtcalc(mu,2) = -1d0 * p(mu,2)
            pvirtcalc(mu,3) = -1d0 * qa(mu)
            pvirtcalc(mu,4) = 1d0 * (p(mu,3) - p(mu,4) - qe(mu))
            temp = pvirtcalc(mu,1) + pvirtcalc(mu,2)
     1             + pvirtcalc(mu,3) + pvirtcalc(mu,4)
          enddo

          ! Box correction to the upper line for photon 2.
          scale = scaleupperline
          call BOX_VBF_CURR(!
     $      pvirtcalc(0,1), pvirtcalc(0,3),              ! momenta 
     $      pvirtcalc(0,4), pvirtcalc(0,2),              ! momenta
     $      psi(1,-1,2),psi(1,-1,1),                     ! spinors
     $      au, jvua1(0,-1,is1,i),jua1(0,-1,i),          ! currents
     $      1,                                           ! Max helicity 
     $      1, 3-i,                                      ! comp, number of current
     $      scale,divmax,                                ! Max helicity(1 spinor 2 ext. boson)
     $      mboxgauge(1,-1,0,0,3),                       ! virtuals gauge
     $      mboxgaugeborn(1,-1,0,0,3),                   ! born gauge
     $      lwardtest(wardIndex))                        ! Gauge test

          wardindex = wardindex + 1
        enddo !i
      endif !lbox

      if (lbox) then
         if(.not. lzs(js3)) then
            print*, "STOP2 in qqaaqq.F"
        ENDIF
      endif

      if (lbox .and. lzs(js3)) then
        do i = 3,4
          
          ! Box correction to the lower line for photon 1.
          do mu = 0,3
            pvirtcalc(mu,1) = p(mu,3)
            pvirtcalc(mu,2) = -1d0 * p(mu,4)
            pvirtcalc(mu,3) = -1d0 * qe(mu)
            pvirtcalc(mu,4) = 1d0 * (p(mu,1) - p(mu,2) - qa(mu))
            temp = pvirtcalc(mu,1) + pvirtcalc(mu,2)
     1             + pvirtcalc(mu,3) + pvirtcalc(mu,4)
          enddo

          scale = scalelowerline
          call BOX_VBF_CURR(!
     $      pvirtcalc(0,1), pvirtcalc(0,3),         ! momenta 
     $      pvirtcalc(0,4), pvirtcalc(0,2),         ! momenta
     $      psi(1,-1,4),psi(1,-1,3),                ! spinors
     $      ae, jvea1(0,-1,is3,i),jea1(0,-1,i),     ! currents
     $      1,                                      ! Max helicity 
     $      1, 5-i,                                 ! comp, number of current
     $      scale, divmax,                          ! Max helicity(1 spinor 2 ext. boson)
     $      mboxgauge(1,-1,0,0,5),                  ! virtuals gauge
     $      mboxgaugeborn(1,-1,0,0,5),              ! born gauge
     $      lwardtest(wardIndex))                   ! Gauge test

          wardindex = wardindex + 1


          ! Box correction to the lower line for photon 2.
          do mu = 0,3
            pvirtcalc(mu,1) = p(mu,3)
            pvirtcalc(mu,2) = -1d0 * p(mu,4)
            pvirtcalc(mu,3) = -1d0 * qa(mu)
            pvirtcalc(mu,4) = 1d0 * (p(mu,1) - p(mu,2) - qe(mu))
            temp = pvirtcalc(mu,1) + pvirtcalc(mu,2)
     1             + pvirtcalc(mu,3) + pvirtcalc(mu,4)
          enddo

          scale = scalelowerline
          call BOX_VBF_CURR(
     $      pvirtcalc(0,1), pvirtcalc(0,3),         ! momenta 
     $      pvirtcalc(0,4), pvirtcalc(0,2),         ! momenta
     $      psi(1,-1,4),psi(1,-1,3),                ! spinors
     $      au, jvua1(0,-1,is3,i),jua1(0,-1,i),     ! currents
     $      1,                                      ! Max helicity 
     $      1, 5-i,                                 ! comp, number of current
     $      scale, divmax,                          ! Max helicity(1 spinor 2 ext. boson)
     $      mboxgauge(1,-1,0,0,6),                  ! virtuals gauge
     $      mboxgaugeborn(1,-1,0,0,6),              ! born gauge
     $      lwardtest(wardIndex))                   ! Gauge test

          wardindex = wardindex + 1
        
        enddo !i

        ! Box corrected quarklines currents.
        jvea = jvea1
        jvua = jvua1
      endif !lbox
#endif

      ! now calculate the Vertex-box diagrams; 
      ! get t-channel W currents first.
      do k = 1,2

        ! k=1 and kk=2 is for "box correction" to upper line
        ! k=2 and kk=1 is for "box correction" to lower line
        kk = 3-k

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

        isig = -1
        if (k.eq.1) then
          call contract_T2j(CCee(0,0,k,L),jqq(0,isig,kk), epsee) 
          call contract_T2j(CCuu(0,0,k,L),jqq(0,isig,kk), epsuu)
          call contract_T1j(CCee6(0,0,k,L),jqq(0,isig,kk), epsee6) 
          call contract_T1j(CCuu6(0,0,k,L),jqq(0,isig,kk), epsuu6)
        else
          call contract_T1j(CCee(0,0,k,L),jqq(0,isig,kk), epsee) 
          call contract_T1j(CCuu(0,0,k,L),jqq(0,isig,kk), epsuu)
          call contract_T2j(CCee6(0,0,k,L),jqq(0,isig,kk), epsee6) 
          call contract_T2j(CCuu6(0,0,k,L),jqq(0,isig,kk), epsuu6)
        endif

        do mu = 0,3 
          epsCCee(mu,isig,k) =
     1         epsee(mu)*clr(3,3,isig) ! coupling to fermion
          epsCCuu(mu,isig,k) =
     1         epsuu(mu)*clr(3,3,isig) 
          epsCCee6(mu,isig,k) =
     1         epsee6(mu)*clr(3,3,isig) ! coupling to fermion
          epsCCuu6(mu,isig,k) =
     1         epsuu6(mu)*clr(3,3,isig) 
        enddo

        ! extra terms for (q_mu q_nu)-contributions 
        ! of massive boson propagators:	 
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

        ! add momentum info to the currents
        do mu = 4,5
          epsCCee(mu,isig,k) = zuu(mu)
          epsCCuu(mu,isig,k) = zee(mu)
          epsCCee6(mu,isig,k) = zee(mu)
          epsCCuu6(mu,isig,k) = zuu(mu)
        enddo 
      enddo !k

      do isig = -1,1,2 
        jj21au(isig) = dotcc(au,jqq(0,isig,1))
        jj21ae(isig) = dotcc(ae,jqq(0,isig,1))

        jj43au(isig) = dotcc(au,jqq(0,isig,2))
        jj43ae(isig) = dotcc(ae,jqq(0,isig,2))
      enddo !isig

      ! now construct the contribution to the amplitude by current contraction 
      ! virtual contributions are assembled in subroutine boxline
      ! 2 bosons attached to 12 line     
      isig = -1
      k = 5                     

      !udsc, photon 1 from upper line.
      m1ae = dotcc(jea(0,-1,2),epsCCuu(0,isig,1))*clr(ifl(2,k),1,-1)   
      m2ae = dotcc(jea(0,-1,1),epsCCuu(0,isig,1))*clr(ifl(1,k),1,-1)   

      z1ae = qepsCCuu(1)*jj21ae(isig)*clr(ifl(1,k),3,isig)
     1       *clr(ifl(2,k),1,isig)
      z2ae = qepsCCuu(1)*jj21ae(isig)*clr(ifl(2,k),3,isig)
     1       *clr(ifl(1,k),1,isig)

      ! Achtung: A an u/d Kopplung unterschiedlich!! 
     
      mate12 = clr(3,3,-1)*(m2ae-z2ae)
     1         + clr(3,3,-1)*(m1ae+z1ae)

      !udsc, photon 2 from upper line.
      m1au = dotcc(jua(0,-1,2),epsCCee(0,isig,1))*clr(ifl(2,k),1,-1)   
      m2au = dotcc(jua(0,-1,1),epsCCee(0,isig,1))*clr(ifl(1,k),1,-1)   

      z1au = qepsCCee(1)*jj21au(isig)*clr(ifl(1,k),3,isig)
     1       *clr(ifl(2,k),1,isig)
      z2au = qepsCCee(1)*jj21au(isig)*clr(ifl(2,k),3,isig)
     1       *clr(ifl(1,k),1,isig)


      matu12 = clr(3,3,-1)*(m2au-z2au)
     1         + clr(3,3,-1)*(m1au+z1au)
 
     
      mat(k,-1,isig,2) = mate12 + matu12

      isig = -1
      k = 6                     

      !udsc, photon 1 from upper line.
      m1ae = dotcc(jea(0,-1,2),epsCCuu6(0,isig,1))*clr(ifl(2,k),1,-1)
      m2ae = dotcc(jea(0,-1,1),epsCCuu6(0,isig,1))*clr(ifl(1,k),1,-1)

      z1ae = qepsCCuu6(1)*jj21ae(isig)*clr(ifl(1,k),3,isig)
     1       *clr(ifl(2,k),1,isig)
      z2ae = qepsCCuu6(1)*jj21ae(isig)*clr(ifl(2,k),3,isig)
     1       *clr(ifl(1,k),1,isig)

      mate12 = clr(3,3,-1)*(m2ae-z2ae)
     1         + clr(3,3,-1)*(m1ae+z1ae)


      !ducs, photon 2 from upper line
      m1au = dotcc(jua(0,-1,2),epsCCee6(0,isig,1))*clr(ifl(2,k),1,-1)
      m2au = dotcc(jua(0,-1,1),epsCCee6(0,isig,1))*clr(ifl(1,k),1,-1)

      z1au = qepsCCee6(1)*jj21au(isig)*clr(ifl(1,k),3,isig)
     1       *clr(ifl(2,k),1,isig)
      z2au = qepsCCee6(1)*jj21au(isig)*clr(ifl(2,k),3,isig)
     1       *clr(ifl(1,k),1,isig)

      matu12 = clr(3,3,-1)*(m2au-z2au)
     1         + clr(3,3,-1)*(m1au+z1au)

     
      mat(k,-1,isig,2) = mate12 + matu12

      if (lbox) then

        isig = -1
        k = 5                     

        !udsc, photon 1 from upper line
        mv1ae = dotcc(jvea(0,-1,is1,2),epsCCuu(0,isig,1))
     1          *clr(ifl(2,k),1,-1)  
        mv2ae = dotcc(jvea(0,-1,is1,1),epsCCuu(0,isig,1))
     1          *clr(ifl(1,k),1,-1)  

        matve12 = clr(3,3,-1)*(mv2ae+mv1ae)

        !udsc, photon 2 from upper line
        mv1au = dotcc(jvua(0,-1,is1,2),epsCCee(0,isig,1))
     1          *clr(ifl(2,k),1,-1)  
        mv2au = dotcc(jvua(0,-1,is1,1),epsCCee(0,isig,1))
     1          *clr(ifl(1,k),1,-1)  

        matvu12 = clr(3,3,-1)*(mv2au+mv1au)
     
        matv(k,-1,isig,2) = matve12 + matvu12

        isig = -1
        k = 6                     

        !udsc, photon 1 from upper line
        mv1ae = dotcc(jvea(0,-1,is1,2),epsCCuu6(0,isig,1))
     1          *clr(ifl(2,k),1,-1) 
        mv2ae = dotcc(jvea(0,-1,is1,1),epsCCuu6(0,isig,1))
     1          *clr(ifl(1,k),1,-1) 
      
        matve12 = clr(3,3,-1)*(mv2ae+mv1ae)

        !udsc, photon 2 from upper line
        mv1au = dotcc(jvua(0,-1,is1,2),epsCCee6(0,isig,1))
     1          *clr(ifl(2,k),1,-1) 
        mv2au = dotcc(jvua(0,-1,is1,1),epsCCee6(0,isig,1))
     1          *clr(ifl(1,k),1,-1) 

        matvu12 = clr(3,3,-1)*(mv2au+mv1au)
     
        matv(k,-1,isig,2) = matve12 + matvu12

      endif !lbox

      
      ! repeat the same for 2 bosons attached to 34 line:
      isig = -1
      k = 5                     

      !udsc, photon 1 from lower line
      m3ae = dotcc(jea(0,-1,4),epsCCuu(0,isig,2))*clr(ifl(4,k),1,-1)
      m4ae = dotcc(jea(0,-1,3),epsCCuu(0,isig,2))*clr(ifl(3,k),1,-1)

      z3ae = qepsCCuu(2)*jj43ae(isig)*clr(ifl(3,k),3,isig)
     1       *clr(ifl(4,k),1,isig)
      z4ae = qepsCCuu(2)*jj43ae(isig)*clr(ifl(4,k),3,isig)
     1       *clr(ifl(3,k),1,isig)
     
      mate34 = clr(3,3,-1)*(m4ae-z4ae)
     1         + clr(3,3,-1)*(m3ae+z3ae)

      !udsc, photon 2 from lower line
      m3au = dotcc(jua(0,-1,4),epsCCee(0,isig,2))*clr(ifl(4,k),1,-1)
      m4au = dotcc(jua(0,-1,3),epsCCee(0,isig,2))*clr(ifl(3,k),1,-1)

      z3au = qepsCCee(2)*jj43au(isig)*clr(ifl(3,k),3,isig)
     1       *clr(ifl(4,k),1,isig)
      z4au = qepsCCee(2)*jj43au(isig)*clr(ifl(4,k),3,isig)
     1       *clr(ifl(3,k),1,isig)


      matu34 = clr(3,3,-1)*(m4au-z4au)
     1         + clr(3,3,-1)*(m3au+z3au)
    
      mat(k,-1,isig,3) = mate34 + matu34

      isig = -1
      k = 6                     

      !udsc, photon 1 from upper line
      m3ae = dotcc(jea(0,-1,4),epsCCuu6(0,isig,2))*clr(ifl(4,k),1,-1)  
      m4ae = dotcc(jea(0,-1,3),epsCCuu6(0,isig,2))*clr(ifl(3,k),1,-1)  

      z3ae = qepsCCuu6(2)*jj43ae(isig)*clr(ifl(3,k),3,isig)
     1       *clr(ifl(4,k),1,isig)
      z4ae = qepsCCuu6(2)*jj43ae(isig)*clr(ifl(4,k),3,isig)
     1       *clr(ifl(3,k),1,isig)

      mate34 = clr(3,3,-1)*(m4ae-z4ae)
     1         + clr(3,3,-1)*(m3ae+z3ae)

      !udsc, photon 2 from upper line
      m3au = dotcc(jua(0,-1,4),epsCCee6(0,isig,2))*clr(ifl(4,k),1,-1)
      m4au = dotcc(jua(0,-1,3),epsCCee6(0,isig,2))*clr(ifl(3,k),1,-1)

      z3au = qepsCCee6(2)*jj43au(isig)*clr(ifl(3,k),3,isig)
     1       *clr(ifl(4,k),1,isig)
      z4au = qepsCCee6(2)*jj43au(isig)*clr(ifl(4,k),3,isig)
     1       *clr(ifl(3,k),1,isig)

      matu34 = clr(3,3,-1)*(m4au-z4au)
     1         + clr(3,3,-1)*(m3au+z3au)
    
      mat(k,-1,isig,3) = mate34 + matu34

      if (lbox) then

        isig = -1
        k = 5                     

        !udsc, photon 1 from lower line
        mv3ae = dotcc(jvea(0,-1,is3,4),epsCCuu(0,isig,2))
     1          *clr(ifl(4,k),1,-1)
        mv4ae = dotcc(jvea(0,-1,is3,3),epsCCuu(0,isig,2))
     1          *clr(ifl(3,k),1,-1)
        
        matve34 = clr(3,3,-1)*(mv4ae+mv3ae)

        !udsc, photon 2 from lower line
        mv3au = dotcc(jvua(0,-1,is3,4),epsCCee(0,isig,2))
     1          *clr(ifl(4,k),1,-1)
        mv4au = dotcc(jvua(0,-1,is3,3),epsCCee(0,isig,2))
     1          *clr(ifl(3,k),1,-1)

        matvu34 = clr(3,3,-1)*(mv4au+mv3au)
     
        matv(k,-1,isig,3) = matve34 + matvu34
        
        isig = -1
        k = 6                     

        !udsc, photon 1 from upper line
        mv3ae = dotcc(jvea(0,-1,is3,4),epsCCuu6(0,isig,2))
     1          *clr(ifl(4,k),1,-1)
        mv4ae = dotcc(jvea(0,-1,is3,3),epsCCuu6(0,isig,2))
     1          *clr(ifl(3,k),1,-1)

        matve34 = clr(3,3,-1)*(mv4ae+mv3ae)


        !udsc, photon 2 from upper line
        mv3au = dotcc(jvua(0,-1,is3,4),epsCCee6(0,isig,2))
     1          *clr(ifl(4,k),1,-1)
        mv4au = dotcc(jvua(0,-1,is3,3),epsCCee6(0,isig,2))
     1          *clr(ifl(3,k),1,-1)


        matvu34 = clr(3,3,-1)*(mv4au+mv3au)
     
        matv(k,-1,isig,3) = matve34 + matvu34
 
      endif !lbox

      ! next do the box-box graphs with one A/Z emitted from the upper and the 
      ! other from the lower line. These are possible for left- and righthanded 
      ! quarks on both lines (in case of NC)
      do isig1 = -1,1,2
        do isig3 = -1,1,2

          ! photon 1 from upper line, photon 2 from lower line:
          meaa(1,isig1,isig3) = dotcc(jea(0,isig1,1),jua(0,isig3,3))
          meaa(2,isig1,isig3) = dotcc(jea(0,isig1,2),jua(0,isig3,3))
          meaa(3,isig1,isig3) = dotcc(jea(0,isig1,1),jua(0,isig3,4))
          meaa(4,isig1,isig3) = dotcc(jea(0,isig1,2),jua(0,isig3,4))

          ! photon2 from upper line, photon 1 from lower line:
          muaa(1,isig1,isig3) = dotcc(jua(0,isig1,1),jea(0,isig3,3))
          muaa(2,isig1,isig3) = dotcc(jua(0,isig1,2),jea(0,isig3,3))
          muaa(3,isig1,isig3) = dotcc(jua(0,isig1,1),jea(0,isig3,4))
          muaa(4,isig1,isig3) = dotcc(jua(0,isig1,2),jea(0,isig3,4))

          if (lbox) then

            ! photon 1 from upper line, photon 2 from lower line

            ! box corr. to upper line:
            mveaa(1,1,isig1,isig3) = dotcc(jvea(0,isig1,is1,1),
     1                                     jua(0,isig3,3))
            mveaa(2,1,isig1,isig3) = dotcc(jvea(0,isig1,is1,2),
     1                                     jua(0,isig3,3))
            mveaa(3,1,isig1,isig3) = dotcc(jvea(0,isig1,is1,1),
     1                                     jua(0,isig3,4))
            mveaa(4,1,isig1,isig3) = dotcc(jvea(0,isig1,is1,2),
     1                                     jua(0,isig3,4))

            ! box corr. to lower line:
            mveaa(1,2,isig1,isig3) = dotcc(jea(0,isig1,1),
     1                                     jvua(0,isig3,is3,3))
            mveaa(2,2,isig1,isig3) = dotcc(jea(0,isig1,2),
     1                                     jvua(0,isig3,is3,3))
            mveaa(3,2,isig1,isig3) = dotcc(jea(0,isig1,1),
     1                                     jvua(0,isig3,is3,4))
            mveaa(4,2,isig1,isig3) = dotcc(jea(0,isig1,2),
     1                                     jvua(0,isig3,is3,4))

            ! photon 2 from upper line, photon 1 from lower line:

            ! box corr. to upper line:	
            mvuaa(1,1,isig1,isig3) = dotcc(jvua(0,isig1,is1,1),
     1                                     jea(0,isig3,3))
            mvuaa(2,1,isig1,isig3) = dotcc(jvua(0,isig1,is1,2),
     1                                     jea(0,isig3,3))
            mvuaa(3,1,isig1,isig3) = dotcc(jvua(0,isig1,is1,1),
     1                                     jea(0,isig3,4))
            mvuaa(4,1,isig1,isig3) = dotcc(jvua(0,isig1,is1,2),
     1                                     jea(0,isig3,4))

            ! box corr. to lower line:	
            mvuaa(1,2,isig1,isig3) = dotcc(jua(0,isig1,1),
     1                                     jvea(0,isig3,is3,3))
            mvuaa(2,2,isig1,isig3) = dotcc(jua(0,isig1,2),
     1                                     jvea(0,isig3,is3,3))
            mvuaa(3,2,isig1,isig3) = dotcc(jua(0,isig1,1),
     1                                     jvea(0,isig3,is3,4))
            mvuaa(4,2,isig1,isig3) = dotcc(jua(0,isig1,2),
     1                                     jvea(0,isig3,is3,4))

          endif !lbox

          !  for the q^mu*q^nu/M_V^2 terms in the gauge boson propagators we need
          do i = 2,3 

            ! type of exchanged boson (Z,W), zero for A
            zeaa(isig1,isig3,i) = jj21ae(isig1)*jj43au(isig3)*zm2i(i)
            zuaa(isig1,isig3,i) = jj21au(isig1)*jj43ae(isig3)*zm2i(i)
          
          enddo !i
          do k = 1,4 

            propbbez(k,isig1,isig3) = 
     1        prop_ee(2)*clr(ifl(2,k),2,isig1)*clr(ifl(3,k),2,isig3)
            propbbuz(k,isig1,isig3) = 
     1        prop_uu(2)*clr(ifl(2,k),2,isig1)*clr(ifl(3,k),2,isig3)
            propbbea(k,isig1,isig3) = 
     1        prop_ee(1)*clr(ifl(2,k),1,isig1)*clr(ifl(3,k),1,isig3)
            propbbua(k,isig1,isig3) = 
     1        prop_uu(1)*clr(ifl(2,k),1,isig1)*clr(ifl(3,k),1,isig3) 
     
            ! for NC make use of: ifl(1,k) = ifl (2,k) and ifl(3,k) = ifl (4,k) 
            ! gauge terms zezz etc cancel in sum of NC contributions 
            ! from leg1&2/3&4 
            ! (not for CC)
            mat(k,isig1,isig3,6) = 
     1          ((meaa(1,isig1,isig3)+meaa(2,isig1,isig3))*
     1          clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3))*
     1          propbbea(k,isig1,isig3)
     1          +
     1          (((meaa(1,isig1,isig3)+zeaa(isig1,isig3,2))+
     1          (meaa(2,isig1,isig3)-zeaa(isig1,isig3,2)))*
     1          clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3))*
     1          propbbez(k,isig1,isig3)
     1          +
     1          ((muaa(1,isig1,isig3)+muaa(2,isig1,isig3))*
     1          clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3))*
     1          propbbua(k,isig1,isig3)
     1          +
     1          (((muaa(1,isig1,isig3)+zuaa(isig1,isig3,2))+
     1           (muaa(2,isig1,isig3)-zuaa(isig1,isig3,2)))*
     1          clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3))*
     1          propbbuz(k,isig1,isig3)
     1          +
     1          ((meaa(3,isig1,isig3)+meaa(4,isig1,isig3))*
     1          clr(ifl(1,k),1,isig1)*clr(ifl(4,k),1,isig3))*
     1          propbbea(k,isig1,isig3)
     1          +
     1          (((meaa(3,isig1,isig3)-zeaa(isig1,isig3,2))+
     1          (meaa(4,isig1,isig3)+zeaa(isig1,isig3,2)))*
     1          clr(ifl(1,k),1,isig1)*clr(ifl(4,k),1,isig3))*
     1          propbbez(k,isig1,isig3)
     1          +
     1          ((muaa(3,isig1,isig3)+muaa(4,isig1,isig3))*
     1          clr(ifl(1,k),1,isig1)*clr(ifl(4,k),1,isig3))*
     1          propbbua(k,isig1,isig3)
     1          +
     1          (((muaa(3,isig1,isig3)-zuaa(isig1,isig3,2))+
     1          (muaa(4,isig1,isig3)+zuaa(isig1,isig3,2)))*
     1          clr(ifl(1,k),1,isig1)*clr(ifl(4,k),1,isig3))*
     1          propbbuz(k,isig1,isig3)   
         enddo !k
        enddo ! isig3
      enddo ! isig1

c -------------------------	

      do k = 5,6 
        ! charged current
        isig1 = -1
        isig3 = -1

        propbbe =  prop_ee(3)*clr(3,3,-1)**2
        propbbu =  prop_uu(3)*clr(3,3,-1)**2
             
        mat(k,isig1,isig3,6) = 
     1    ((meaa(1,isig1,isig3)+zeaa(isig1,isig3,3))*
     1    clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3))*
     1    propbbe     
     1    +
     1    ((meaa(2,isig1,isig3)-zeaa(isig1,isig3,3))*
     1    clr(ifl(2,k),1,isig1)*clr(ifl(3,k),1,isig3))*
     1    propbbe     
     1    +
     1    ((meaa(3,isig1,isig3)-zeaa(isig1,isig3,3))*
     1    clr(ifl(1,k),1,isig1)*clr(ifl(4,k),1,isig3))*
     1    propbbe     
     1    +
     1    ((meaa(4,isig1,isig3)+zeaa(isig1,isig3,3))*
     1    clr(ifl(2,k),1,isig1)*clr(ifl(4,k),1,isig3))*
     1    propbbe 
     1    +        
     1    ((muaa(1,isig1,isig3)+zuaa(isig1,isig3,3))*
     1    clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3))*
     1    propbbu    
     1    +
     1    ((muaa(2,isig1,isig3)-zuaa(isig1,isig3,3))*
     1    clr(ifl(2,k),1,isig1)*clr(ifl(3,k),1,isig3))*
     1    propbbu     
     1    +
     1    ((muaa(3,isig1,isig3)-zuaa(isig1,isig3,3))*
     1    clr(ifl(1,k),1,isig1)*clr(ifl(4,k),1,isig3))*
     1    propbbu    
     1    +
     1    ((muaa(4,isig1,isig3)+zuaa(isig1,isig3,3))*
     1    clr(ifl(2,k),1,isig1)*clr(ifl(4,k),1,isig3))*
     1    propbbu
      enddo !k
c
c---------------------------------------------------
c      
      if (lbox) then 
        ! no need to include gauge terms (contraction of box with 
        ! momentum -> triangle)	
        do i = 1,2         ! 1 is for upper line, 2 for lower line QCD correction
          jj = 3 + 3*i     ! stored in matv(...,6) and matv(...,9) respectively 
          do isig1 = -1,1,2
            do isig3 = -1,1,2
              do k = 1,4 
                matv(k,isig1,isig3,jj) = 
     1             ((mveaa(1,i,isig1,isig3)+mveaa(2,i,isig1,isig3))*
     1             clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3))*
     1             propbbea(k,isig1,isig3)
     1             +
     1             (((mveaa(1,i,isig1,isig3))+(mveaa(2,i,isig1,isig3)))*
     1             clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3))*
     1             propbbez(k,isig1,isig3)
     1             +
     1             ((mvuaa(1,i,isig1,isig3)+mvuaa(2,i,isig1,isig3))*
     1             clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3))*
     1             propbbua(k,isig1,isig3)
     1             +
     1             (((mvuaa(1,i,isig1,isig3))+(mvuaa(2,i,isig1,isig3)))*
     1             clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3))*
     1             propbbuz(k,isig1,isig3)
     1             +
     1             ((mveaa(3,i,isig1,isig3)+mveaa(4,i,isig1,isig3))*
     1             clr(ifl(1,k),1,isig1)*clr(ifl(4,k),1,isig3))*
     1             propbbea(k,isig1,isig3)
     1             +
     1             (((mveaa(3,i,isig1,isig3))+(mveaa(4,i,isig1,isig3)))*
     1             clr(ifl(1,k),1,isig1)*clr(ifl(4,k),1,isig3))*
     1             propbbez(k,isig1,isig3)
     1             +
     1             ((mvuaa(3,i,isig1,isig3)+mvuaa(4,i,isig1,isig3))*
     1             clr(ifl(1,k),1,isig1)*clr(ifl(4,k),1,isig3))*
     1             propbbua(k,isig1,isig3)
     1             +
     1             (((mvuaa(3,i,isig1,isig3))+(mvuaa(4,i,isig1,isig3)))*
     1             clr(ifl(1,k),1,isig1)*clr(ifl(4,k),1,isig3))*
     1             propbbuz(k,isig1,isig3)
              enddo !k
            enddo !isig3
          enddo !isig1

          do k = 5,6
            isig1 = -1
            isig3 = -1
     
            matv(k,isig1,isig3,jj) = 
     1         (mveaa(1,i,isig1,isig3)*
     1         clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3))*
     1         propbbe
     2         +
     2         (mveaa(2,i,isig1,isig3)*
     2         clr(ifl(2,k),1,isig1)*clr(ifl(3,k),1,isig3))*
     2         propbbe
     3         +
     3         (mveaa(3,i,isig1,isig3)*
     3         clr(ifl(1,k),1,isig1)*clr(ifl(4,k),1,isig3))*
     3         propbbe
     4         +
     4         (mveaa(4,i,isig1,isig3)*
     4         clr(ifl(2,k),1,isig1)*clr(ifl(4,k),1,isig3))*
     1         propbbe
     1         +
     1         (mvuaa(1,i,isig1,isig3)*
     1         clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3))*
     1         propbbu
     2         +
     2         (mvuaa(2,i,isig1,isig3)*
     2         clr(ifl(2,k),1,isig1)*clr(ifl(3,k),1,isig3))*
     2         propbbu
     3         +
     3         (mvuaa(3,i,isig1,isig3)*
     3         clr(ifl(1,k),1,isig1)*clr(ifl(4,k),1,isig3))*
     3         propbbu
     4         +
     4         (mvuaa(4,i,isig1,isig3)*
     4         clr(ifl(2,k),1,isig1)*clr(ifl(4,k),1,isig3))*
     4         propbbu

          enddo !k
        enddo !i
      endif !lbox	 

c ------------------------------------------------------------------------------
c and now, finally, the pentagon contributions, i.e. two W"s emitted from the
c  same quark line
c
c prerequisites for virtual corrections:
c ------------------------------------------------------------------------------

#ifdef WITH_NLO 	
      if (lpent .and. lzs(js1)) then     
        ! need new pentagon graphs for 12 line

        musq = -p21(4)
        if (lpent) then
          do mu = 0,3
            pvirtcalc(mu,1) = p(mu,1)
            pvirtcalc(mu,2) = -1d0 * p(mu,2)
            pvirtcalc(mu,3) = -1d0 * qe(mu)
            pvirtcalc(mu,4) = -1d0 * qa(mu)
            pvirtcalc(mu,5) = 1d0 * (p(mu,3) - p(mu,4))
            temp = pvirtcalc(mu,1) + pvirtcalc(mu,2)
     1             + pvirtcalc(mu,3) + pvirtcalc(mu,4) 
     1             + pvirtcalc(mu,5)
          enddo

          scale = scaleupperline
          do k = 1,3
            call PENT_VBF_CURRENT(
     $         pvirtcalc(0,1), pvirtcalc(0,3),                  ! momenta 
     $         pvirtcalc(0,4), pvirtcalc(0,5),                  ! momenta
     $         pvirtcalc(0,2),                                  ! momenta
     $         psi(1,-1,2),psi(1,-1,1),                         ! spinors
     $         ae, au, j5aaeu1(0,-1,k,is1,1),jbaaeu1(0,-1,k,1), ! currents
     $         1,                                               ! Max helicity 
     $         1, k,                                            ! comp, number of current
     $         scale, divmax,                                   ! Max helicity(1 spinor 2 ext. boson)
     $         mpengauge(1,-1,0,0,1),                           ! virtuals gauge
     $         mpengaugeborn(1,-1,0,0,1),                       ! born gauge
     $         lwardtestp(wardIndexp))                          ! Gauge test

            wardindexp = wardindexp+1

            call PENT_VBF_CURRENT(
     $         pvirtcalc(0,1), pvirtcalc(0,4),                  ! momenta 
     $         pvirtcalc(0,3), pvirtcalc(0,5),                  ! momenta
     $         pvirtcalc(0,2),                                  ! momenta
     $         psi(1,-1,2),psi(1,-1,1),                         ! spinors
     $         au,ae, j5aaue1(0,-1,k,is1,1),jbaaue1(0,-1,k,1),  ! currents
     $         1,                                               ! Max helicity 
     $         1, k,                                            ! comp, number of current
     $         scale, divmax,                                   ! Max helicity(1 spinor 2 ext. boson)
     $         mpengauge(1,-1,0,0,1),                           ! virtuals gauge
     $         mpengaugeborn(1,-1,0,0,1),                       ! born gauge
     $         lwardtestp(wardIndexp))                          ! Gauge test

            wardindexp = wardindexp + 1
          enddo

          ! check for what is used this:
          icount1= icount1+1

          j5aaeu(:,:,:,:,1) = j5aaeu1(:,:,:,:,1)
          j5aaue(:,:,:,:,1) = j5aaue1(:,:,:,:,1)

        endif

cfc TO DO: implement gauge test

        if (.not.lvirtratio) then
          lzs(js1) = .false.
        endif
      endif


      ! pentagon corr. to lower line:

      if (lpent .and. lzs(js3)) then     ! need new pentagon graphs for 34 line
        musq = -p43(4)
        if (lpent) then
          do mu = 0,3
            pvirtcalc(mu,1) = p(mu,3)
            pvirtcalc(mu,2) = -1d0 * p(mu,4)
            pvirtcalc(mu,3) = -1d0 * qe(mu)
            pvirtcalc(mu,4) = -1d0 * qa(mu)
            pvirtcalc(mu,5) = 1d0 * (p(mu,1) - p(mu,2))
            temp = pvirtcalc(mu,1) + pvirtcalc(mu,2) + pvirtcalc(mu,3)
     1             + pvirtcalc(mu,4) + pvirtcalc(mu,5)
          enddo
          scale = scalelowerline

          do k = 1,3
            call PENT_VBF_CURRENT(
     $         pvirtcalc(0,1), pvirtcalc(0,3),                     ! momenta 
     $         pvirtcalc(0,4), pvirtcalc(0,5),                     ! momenta
     $         pvirtcalc(0,2),                                     ! momenta
     $         psi(1,-1,4),psi(1,-1,3),                            ! spinors
     $         ae, au, j5aaeu1(0,-1,k,is3,2),jbaaeu1(0,-1,k,2),    ! currents
     $         1,                                                  ! Max helicity 
     $         1, k,                                               ! comp, number of current
     $         scale,divmax,                                       ! Max helicity(1 spinor 2 ext. boson)
     $         mpengauge(1,-1,0,0,1),                              ! virtuals gauge
     $         mpengaugeborn(1,-1,0,0,1),                          ! born gauge
     $         lwardtestp(wardIndexp))                             ! Gauge test

            wardindexp = wardindexp + 1

            call PENT_VBF_CURRENT(
     $         pvirtcalc(0,1), pvirtcalc(0,4),                     ! momenta 
     $         pvirtcalc(0,3), pvirtcalc(0,5),                     ! momenta
     $         pvirtcalc(0,2),                                     ! momenta
     $         psi(1,-1,4),psi(1,-1,3),                            ! spinors
     $         au,ae, j5aaue1(0,-1,k,is3,2),jbaaue1(0,-1,k,2),     ! currents
     $         1,                                                  ! Max helicity 
     $         1, k,                                               ! comp, number of current
     $         scale,divmax,                                       ! Max helicity(1 spinor 2 ext. boson)
     $         mpengauge(1,-1,0,0,1),                              ! virtuals gauge
     $         mpengaugeborn(1,-1,0,0,1),                          ! born gauge
     $         lwardtestp(wardIndexp))                             ! Gauge test

            wardindexp = wardindexp + 1
          enddo
          icount2= icount2+1

          j5aaeu(:,:,:,:,2) = j5aaeu1(:,:,:,:,2)
          j5aaue(:,:,:,:,2) = j5aaue1(:,:,:,:,2)

        endif

        if (.not.lvirtratio) then
          lzs(js3) = .false.
        endif
      endif
#endif

      ! compute matrix elements for pentagon on upper line:
      ! need m5e and m5u(pos.j43=1:3,vext1=1:3,vext2=1:3), 
      ! id of internal boson arbitrary 
      do isig3 = -1,1,2 
        do isig1 = -1,1,2

          call ket2c(psi(1,isig1,1),.true.,p(0,1),isig1,
     1         p43,jqq(0,isig3,2),bkjqq(1,isig1,isig3,1),dummy)
          call bra2c(psi(1,isig1,2),.true.,p(0,2),isig1,
     1         p43,jqq(0,isig3,2),bkjqq(1,isig1,isig3,2),dummy)

          !  eps1=j43 with v1 arbitrary,eps2=v2,eps3=v3
          m5e(1,1,1) = -s1c(psiau(1,isig1,2),ae,.true.,isig1,
     1                      bkjqq(1,isig1,isig3,1))
          m5u(1,1,1) = -s1c(psiae(1,isig1,2),au,.true.,isig1,
     1                      bkjqq(1,isig1,isig3,1))
          ! eps1=v1,eps2=j43 with v2 arbitrary,eps3=v3
          m5e(2,1,1) = -s1c(psiau(1,isig1,2),jqq(0,isig3,2),.true.,
     1                      isig1,psiae(1,isig1,1))
          m5u(2,1,1) = -s1c(psiae(1,isig1,2),jqq(0,isig3,2),.true.,
     1                      isig1,psiau(1,isig1,1))
          ! eps1=v1,eps2=v2,eps3=j43 with v3 arbitrary
          m5e(3,1,1) = -s1c(bkjqq(1,isig1,isig3,2),au,.true.,isig1,
     1                      psiae(1,isig1,1))
          m5u(3,1,1) = -s1c(bkjqq(1,isig1,isig3,2),ae,.true.,isig1,
     1                      psiau(1,isig1,1))

          do k = 1,6
            m1p(k,1,1) = (m5e(1,1,1) + m5u(1,1,1))*
     1            clr(ifl(2,k),1,isig1)*clr(ifl(2,k),1,isig1)   

            m2p(k,1,1) = ( m5e(2,1,1) + m5u(2,1,1))*
     1           clr(ifl(1,k),1,isig1)*clr(ifl(2,k),1,isig1)

            m3p(k,1,1) = ( m5e(3,1,1) + m5u(3,1,1))*
     1          clr(ifl(1,k),1,isig1)*clr(ifl(1,k),1,isig1)
          enddo
 
          ! virtual stuff:         
          ! need mv5e and mv5u(pos.j43 = 1:3, vext1=1:3, vext2=1:3), 
          ! id of internal boson arbitrary 
          if (lpent) then
            do j = 1,3
              ! epsj = j43 with vj arbitrary,epsa = va,epsb = vb 
              ! (j = 1 or 2 or 3; a,b = A / Z)
              mv5e(j,1,1) = -dotcc(j5aaeu(0,isig1,j,is1,1),
     1                             jqq(0,isig3,2))
              mv5u(j,1,1) = -dotcc(j5aaue(0,isig1,j,is1,1),
     1                             jqq(0,isig3,2))
            enddo
  
            do k = 1,6
              mv1p(k,1,1) = (mv5e(1,1,1) + mv5u(1,1,1))*
     1            clr(ifl(2,k),1,isig1)*clr(ifl(2,k),1,isig1)   
              mv2p(k,1,1) = (mv5e(2,1,1) + mv5u(2,1,1))*
     1          clr(ifl(1,k),1,isig1)*clr(ifl(2,k),1,isig1)
              mv3p(k,1,1) = (mv5e(3,1,1) + mv5u(3,1,1))*
     1          clr(ifl(1,k),1,isig1)*clr(ifl(1,k),1,isig1)
            enddo !k =1,6
          endif

          do k = 1,4  
            mat(k,isig1,isig3,7) = -m1p(k,1,1)*
     1        (clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3)*prop43(1)+
     1        clr(ifl(1,k),2,isig1)*clr(ifl(3,k),2,isig3)*prop43(2))
     1        -
     1        m2p(k,1,1)*
     1        (clr(ifl(3,k),1,isig3)*clr(ifl(2,k),1,isig1)*prop43(1)+
     1        clr(ifl(3,k),2,isig3)*clr(ifl(2,k),2,isig1)*prop43(2))
     1        -
     1        m3p(k,1,1)*
     1        (clr(ifl(2,k),1,isig1)*clr(ifl(3,k),1,isig3)*prop43(1)+
     1        clr(ifl(2,k),2,isig1)*clr(ifl(3,k),2,isig3)*prop43(2))
     
            if (lpent) then
              matv(k,isig1,isig3,7) =  
     1          -mv1p(k,1,1)*
     1          (clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3)*prop43(1)+
     1          clr(ifl(1,k),2,isig1)*clr(ifl(3,k),2,isig3)*prop43(2))
     1          -
     1          mv2p(k,1,1)*
     1          (clr(ifl(3,k),1,isig3)*clr(ifl(2,k),1,isig1)*prop43(1)+
     1          clr(ifl(3,k),2,isig3)*clr(ifl(2,k),2,isig1)*prop43(2))
     1          -
     1          mv3p(k,1,1)*
     1          (clr(ifl(2,k),1,isig1)*clr(ifl(3,k),1,isig3)*prop43(1)+
     1          clr(ifl(2,k),2,isig1)*clr(ifl(3,k),2,isig3)*prop43(2))
            endif
          enddo !k
  
          if (isig1.eq.-1.and.isig3.eq.-1) then
            fac = clr(3,3,-1)**2*prop43(3)
            do k = 5,6    
              mat(k,isig1,isig3,7) = 
     1          -(m1p(k,1,1) + m2p(k,1,1) + m3p(k,1,1))*fac
      
              if (lpent) then
                matv(k,isig1,isig3,7) = 
     1              -(mv1p(k,1,1) + mv2p(k,1,1) + mv3p(k,1,1))*fac
              endif
            enddo !k
          endif !isig13
        enddo !isig1
      enddo !isig3
      
      ! lower line: mat(k,-1,-1,8) = propt(-1,-1,k,1) with prop21
      ! need m5e and m5u(pos.j21=1:3,vext1=1:3,vext2=1:3), vint arbitrary 
      do isig1 = -1,1,2
        do isig3 = -1,1,2
 
          call ket2c(psi(1,isig3,3),.true.,p(0,3),isig3,
     1         p21,jqq(0,isig1,1),bkjqq(1,isig3,isig1,3),dummy)
          call bra2c(psi(1,isig3,4),.true.,p(0,4),isig3,
     1         p21,jqq(0,isig1,1),bkjqq(1,isig3,isig1,4),dummy)

          ! eps1=j21 with v1 arbitrary,eps2=v2,eps3=v3
          m5e(1,1,1) = -s1c(psiau(1,isig3,4),ae,.true.,isig3,
     1                      bkjqq(1,isig3,isig1,3))
          m5u(1,1,1) = -s1c(psiae(1,isig3,4),au,.true.,isig3,
     1                      bkjqq(1,isig3,isig1,3))
          ! eps1=v1,eps2=j21 with v2 arbitrary,eps3=v3
          m5e(2,1,1) = -s1c(psiau(1,isig3,4),jqq(0,isig1,1),.true.,
     1                      isig3,psiae(1,isig3,3))
          m5u(2,1,1) = -s1c(psiae(1,isig3,4),jqq(0,isig1,1),.true.,
     1                      isig3,psiau(1,isig3,3))
         ! eps1=v1,eps2=v2,eps3=j21 with v3 arbitrary
          m5e(3,1,1) = -s1c(bkjqq(1,isig3,isig1,4),au,.true.,isig3,
     1                      psiae(1,isig3,3))
          m5u(3,1,1) = -s1c(bkjqq(1,isig3,isig1,4),ae,.true.,isig3,
     1                      psiau(1,isig3,3))
    
          do k = 1,6
            m1p(k,1,1) = (m5e(1,1,1) + m5u(1,1,1))*
     1          clr(ifl(4,k),1,isig3)*clr(ifl(4,k),1,isig3)   
            m2p(k,1,1) = (m5e(2,1,1) + m5u(2,1,1))*
     1          clr(ifl(3,k),1,isig3)*clr(ifl(4,k),1,isig3)
            m3p(k,1,1) = (m5e(3,1,1) + m5u(3,1,1))*
     1          clr(ifl(3,k),1,isig3)*clr(ifl(3,k),1,isig3)
          enddo

          ! need mv5e and mv5u(pos.j21=1:3,vext1=1:3,vext2=1:3), 
          ! id of internal boson arbitrary 
          if (lpent) then
            do j = 1,3

              ! epsj=j21 with vj arbitrary,epsa=va,epsb=vb 
              ! (j = 1 or 2 or 3; a,b = A)
              mv5e(j,1,1) = -dotcc(j5aaeu(0,isig3,j,is3,2),
     1                             jqq(0,isig1,1))
              mv5u(j,1,1) = -dotcc(j5aaue(0,isig3,j,is3,2),
     1                             jqq(0,isig1,1))
            enddo

            do k = 1,6
              mv1p(k,1,1) = (mv5e(1,1,1) + mv5u(1,1,1))*
     1            clr(ifl(4,k),1,isig3)*clr(ifl(4,k),1,isig3)
              mv2p(k,1,1) = (mv5e(2,1,1) + mv5u(2,1,1))*
     1            clr(ifl(3,k),1,isig3)*clr(ifl(4,k),1,isig3)
              mv3p(k,1,1) = (mv5e(3,1,1) + mv5u(3,1,1))*
     1            clr(ifl(3,k),1,isig3)*clr(ifl(3,k),1,isig3)
            enddo
          endif !lpent

          do k = 1,4  
            mat(k,isig1,isig3,8) = 
     1          - m1p(k,1,1)*
     1          (clr(ifl(3,k),1,isig3)*clr(ifl(1,k),1,isig1)*prop21(1)+
     1          clr(ifl(3,k),2,isig3)*clr(ifl(1,k),2,isig1)*prop21(2))
     1          - m2p(k,1,1)*
     1          (clr(ifl(1,k),1,isig1)*clr(ifl(4,k),1,isig3)*prop21(1)+
     1          clr(ifl(1,k),2,isig1)*clr(ifl(4,k),2,isig3)*prop21(2))
     1          - m3p(k,1,1)*
     1          (clr(ifl(4,k),1,isig3)*clr(ifl(1,k),1,isig1)*prop21(1)+
     1          clr(ifl(4,k),2,isig3)*clr(ifl(1,k),2,isig1)*prop21(2))
    
            if (lpent) then
              matv(k,isig1,isig3,8) = 
     1          - mv1p(k,1,1)*
     1          (clr(ifl(3,k),1,isig3)*clr(ifl(1,k),1,isig1)*prop21(1)+
     1          clr(ifl(3,k),2,isig3)*clr(ifl(1,k),2,isig1)*prop21(2))
     1          - mv2p(k,1,1)*
     1          (clr(ifl(1,k),1,isig1)*clr(ifl(4,k),1,isig3)*prop21(1)+
     1          clr(ifl(1,k),2,isig1)*clr(ifl(4,k),2,isig3)*prop21(2))
     1          - mv3p(k,1,1)*
     1          (clr(ifl(4,k),1,isig3)*clr(ifl(1,k),1,isig1)*prop21(1)+
     1          clr(ifl(4,k),2,isig3)*clr(ifl(1,k),2,isig1)*prop21(2))
            endif
          enddo !k

          if (isig1.eq.-1.and.isig3.eq.-1) then
            fac = clr(3,3,-1)**2*prop21(3)
            do k = 5,6
              mat(k,isig1,isig3,8) = 
     1          - (m1p(k,1,1) + m2p(k,1,1) + m3p(k,1,1))*fac
              if (lpent) then
                matv(k,isig1,isig3,8) = 
     1          - (mv1p(k,1,1) + mv2p(k,1,1) + mv3p(k,1,1))*fac
              endif
            enddo !k
          endif !isig13
        enddo !isig3
      enddo !isig1

c ------------------------------------------------------------------------
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
     1                   matv(k,isig1,isig3,i).ne.0d0) then
                       print*,"matv != 0, mat = 0"
                       print*,"k,isig1,isig3,i,matv",
     1                        k,isig1,isig3,i,matv(k,isig1,isig3,i)
                    endif                        
                    if(mat(k,isig1,isig3,6).ne.0d0) then
                       print*,"k, i, sig1, sig3, mborn",
     1                        k,i,isig1,isig3,mat(k,isig1,isig3,6)
                       print*,"mvirt",matv(k,isig1,isig3,i)
                       ratio = (matv(k,isig1,isig3,i))/
     1                         (mat(k,isig1,isig3,6)*K_contcheck) 
                       print*,"ratio", ratio
                    endif
                  else
                    if( mat(k,isig1,isig3,i).eq.0d0 .and.
     1                   matv(k,isig1,isig3,i).ne.0d0) then
                       print*,"matv != 0, mat = 0"
                       print*,"k,isig1,isig3,i,matv",
     1                         k,isig1,isig3,i,matv(k,isig1,isig3,i)
                    endif
                    if(mat(k,isig1,isig3,i).ne.0d0) then
                       print*,"k, i, sig1, sig3, mborn",
     1                        k,i,isig1,isig3,mat(k,isig1,isig3,i)
                       print*,"mvirt",matv(k,isig1,isig3,i)
                       ratio = matv(k,isig1,isig3,i)/
     1                         (mat(k,isig1,isig3,i)*K_contcheck) 
                       print*,"ratio", ratio
                    endif
                  endif
                enddo
              enddo
            enddo
            print*
            write(*,*)
          enddo               !i
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
     1                            mat(k,isig1,isig3,i)
            enddo
            res(k) = res(k) + dreal(mm(k,isig1,isig3))**2
     1                      + dimag(mm(k,isig1,isig3))**2
            if (lnlo) then
              mv12(k,isig1,isig3) = 
     1            + matv(k,isig1,isig3,2) + matv(k,isig1,isig3,4) + 
     1            matv(k,isig1,isig3,6) + matv(k,isig1,isig3,7) 
              mv34(k,isig1,isig3) =
     1            + matv(k,isig1,isig3,3) + matv(k,isig1,isig3,5) +
     1            matv(k,isig1,isig3,9) + matv(k,isig1,isig3,8)

              if (lpt) then   ! Pentagon are considered for "good PS points" only
                mv12(k,isig1,isig3)=mv12(k,isig1,isig3)*xgc1
              endif
              if (lpt) then   ! Pentagon contributes are considered for "good PS points" only
                mv34(k,isig1,isig3)=mv34(k,isig1,isig3)*xgc2
              endif

              !  add Born type term and multiply by F_q = alphas*C_2/4pi
              !  the factor pi^2/3+9/2 for the born term is after adding the subtraction term
              !  and the counter term for the renormalization of the pdfs
              if (nlo.gt.0) then
                mv12(k,isig1,isig3) = als(1,1)*c2o4pi*
     1           ( mv12(k,isig1,isig3) + mm(k,isig1,isig3)*cvirt )
                mv34(k,isig1,isig3) = als(2,1)*c2o4pi*
     1           ( mv34(k,isig1,isig3) + mm(k,isig1,isig3)*cvirt )
              else
                mv12(k,isig1,isig3) = 
     1               als(1,1)*c2o4pi*mv12(k,isig1,isig3)
                mv34(k,isig1,isig3) = 
     1               als(2,1)*c2o4pi*mv34(k,isig1,isig3)
              endif
              resv(k) = resv(k) + 2*dreal(
     1             mm(k,isig1,isig3)   *
     1          conjg( mv12(k,isig1,isig3)+mv34(k,isig1,isig3) )  )
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
 
      ! eliminate processes with photon virtuality below cutoff
      if (abs(p21(4)).lt.qsqAmin .or.  
     1     abs(p43(4)).lt.qsqAmin) then
        do k = 1,4
          res(k) = 0
          tree(k) = 0
        enddo
      endif


c----------------------
c TO SET AMPLITUDES TOZEO

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
          boxcountm = boxcountm+1
          boxcountps(ps_number) = boxcountps(ps_number)+1
        else
          boxcount2m = boxcount2m+1
          boxcount2ps(ps_number) = boxcount2ps(ps_number)+1
        endif
      endif

      if(lpent) then
        do i = 1,wardIndexp-1
          lwardpent = lwardpent.and.lwardtestp(i)
        enddo                  ! i
        if(.not.lwardPent.or.res(1).eq.0d0) then
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
     1      res(blha_idsubproc(blha_cursubproc,blha_curproc))

        else if (blha_amptype(blha_curproc).eq.1) then ! loop
          lnmusqp1p2 = log(blha_scale**2/(-p21(4)))
          lnmusqp3p4 = log(blha_scale**2/(-p43(4)))

          blha_amp(1)=blha_amp(1)+
     1      blha_CF/(2d0*pi)*2d0*
     1       (als(1,1)+als(2,1))*(-1d0)*
     1       tree(blha_idsubproc(blha_cursubproc,blha_curproc))
          blha_amp(2)=blha_amp(2)+
     1      1d0/(2d0*pi)*2d0*
     1       ( als(1,1)*(-blha_gammaQuark-blha_CF*lnmusqp1p2)
     1        +als(2,1)*(-blha_gammaQuark-blha_CF*lnmusqp3p4) )*
     1       tree(blha_idsubproc(blha_cursubproc,blha_curproc))
          blha_amp(3)= blha_amp(3)+
     1      (res(blha_idsubproc(blha_cursubproc,blha_curproc))     ! tree + loop + I
     1       -tree(blha_idsubproc(blha_cursubproc,blha_curproc)))  ! -tree
     1       * blha_CF/c2                                          ! correct for Nc ?= inf
     1      + (-(als(1,1)+als(2,1))/2d0/pi*                        ! -I
     1          (blha_CF*cvirt+
     1            2d0*(blha_KQuark+blha_gammaQuark+
     1                 blha_tgammaQuark+(pi**2/6d0-1)*blha_CF) )
     1         +als(1,1)/pi*                                       ! +(mu^2/2p1.p2)^epsilon
     1         (-blha_gammaQuark*lnmusqp1p2-blha_CF/2d0*lnmusqp1p2**2)
     1         +als(2,1)/pi*                                       ! +(mu^2/2p3.p4)^epsilon
     1         (-blha_gammaQuark*lnmusqp3p4-blha_CF/2d0*lnmusqp3p4**2)
     1        )*tree(blha_idsubproc(blha_cursubproc,blha_curproc))
          blha_amp(4)= blha_amp(4)+
     1       tree(blha_idsubproc(blha_cursubproc,blha_curproc))

        else if (blha_amptype(blha_curproc).eq.2) then ! cctree
          ! 1-2
          call BLHA_cctree(1,2,
     1      -blha_CF*
     1       tree(blha_idsubproc(blha_cursubproc,blha_curproc))
     1    )
          ! 3-4
          call BLHA_cctree(3,4,
     1      -blha_CF*
     1       tree(blha_idsubproc(blha_cursubproc,blha_curproc))
     1    )

        else 
          call BLHA_amptypeerror(blha_amptype(blha_curproc),
     1                           __FILE__,__LINE__)
        endif
      endif

      Lold = L
      return
      end
