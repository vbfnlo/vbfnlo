c---------------------------------------------------------------------------
c	lepton helicities have to be summed over externally 
c		( this routine is for one fixed set of lep. helicities only )
c 
c---------------------------------------------------------------------------
      
      subroutine qqwpzqq(pbar,sign, nlo, L,bos,
     1                  ans1,ans2,ans3,ans4)
          use globalvars, only: ldoblha
      implicit none
c
c	Barbara Jaeger, <baj@particle.physik.uni-karlsruhe.de>
c	Initial version:  2006 Jan 10
c	Last modified:     
C
C  QQWZQQ calculates the matrix elements**2 for electroweak
c  weak boson pair production by quark quark scattering
C
C        q1 q3    ---->   q2 q4 W Z,   W ---> f5 f6, Z ---> f7 f8
C
c  where the final state leptons are 3 charged leptons + 1 neutrino. 
C  Crossing related processes can be computed as well. Pauli interference 
c  terms for identical fermions are neglected. In particular, only the
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
c        remaining cases:   return uucc = 2Re(M_Born^* M_virt), more precisely
c        the piece which is missing for the full NLO=1 when calling +|NLO|
c          NLO = -4      finite Born + q^mu/m_V terms of Pentagon contributions only, no cvirt
c          NLO = -5      true finite Pentagon contributions only
c
c  integer L (=1,2) is used for the subtraction code to select leptonic 
c  tensors for gluon emission off the upper Line (L=1) or lower quark 
c  line (L=2)
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"
#include "tensorw.inc"
c#include "tensorw_gaugev.inc"
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
      double precision ans1, ans2, ans3, ans4,res(4),resv(4)

      double precision tree(6)
      double precision  p(0:3,4+max_v), p21(0:4), p43(0:4), pew(0:4),
     1                  puz(0:4), pwz(0:4)
      integer  sign(4+max_v), nlo, mu, i, j, jj, k, kk,
     1         isig, isig1, isig3,bos
      integer  ifl(4,6), js1, js3, L, Lold, is1, is3
      double complex prop21(4), prop43(4)
      double complex mat(4,-1:1,-1:1,9), matv(4,-1:1,-1:1,9)
      double complex mm(4,-1:1,-1:1), 
     1               mv12(4,-1:1,-1:1), mv34(4,-1:1,-1:1)
      double complex maw,mzw,mwa,mwz
      double complex m1,m2,m3,m4,rat
      double complex mw(2),mvw(2)
      double complex mbwa(2,3:4,-1:1,-1:1),mbwz(2,3:4,-1:1,-1:1),
     1 		     mbaw(2,3:4,-1:1,-1:1),mbzw(2,3:4,-1:1,-1:1),
     1 		     mbvwa(2,3:4,-1:1,-1:1,2),mbvwz(2,3:4,-1:1,-1:1,2),
     1 		     mbvaw(2,3:4,-1:1,-1:1,2),mbvzw(2,3:4,-1:1,-1:1,2),
     1		     zwz(-1:1,-1:1,3),zwa(-1:1,-1:1,3),
     1		     zzw(-1:1,-1:1,3),zaw(-1:1,-1:1,3),
     1		     mbwv,mbvw,gbwv,gbvw,mbvwv,mbvvw
      double complex m5(3,3,3),mv5(3,3,3)
      double complex m1ae,m1ze,z1ze,mw1ae,mw1ze,m2ae,m2ze,mw2ae,mw2ze,
     1 		     m3ae,m3ze,z3ze,mw3ae,mw3ze,m4ae,m4ze,mw4ae,mw4ze,
     1		     m1we,m2we,m3we,m4we,w1we,w2we,w3we,w4we,
     2		     m1e,m2e,m3e,m4e,
     3	             m1aau,m1azu,m1zau,m1zzu,z1azu,z1zzu,
     4		     mw1aau,mw1azu,mw1zau,mw1zzu,
     3	             m2aau,m2azu,m2zau,m2zzu,
     4		     mw2aau,mw2azu,mw2zau,mw2zzu,
     3	             m3aau,m3azu,m3zau,m3zzu,z3azu,z3zzu,
     4		     mw3aau,mw3azu,mw3zau,mw3zzu,
     3	             m4aau,m4azu,m4zau,m4zzu,
     4		     mw4aau,mw4azu,mw4zau,mw4zzu,
     2		     m1u,m2u,m3u,m4u,
     2		     m1awu,m1zwu,z1awu,z1zwu,m1wau,m1wzu,
     2		     m2awu,m2zwu,m2wau,m2wzu,
     2		     m3awu,m3zwu,z3awu,z3zwu,m3wau,m3wzu,
     2		     m4awu,m4zwu,m4wau,m4wzu
      double complex mvw1ae,mvw1ze,mvw2ae,mvw2ze,
     1 		     mvw3ae,mvw3ze,mvw4ae,mvw4ze,
     2		     mv1e,mv2e,mv3e,mv4e,
     3	             mvw1aau,mvw1azu,mvw1zau,mvw1zzu,
     3	             mvw2aau,mvw2azu,mvw2zau,mvw2zzu,
     3	             mvw3aau,mvw3azu,mvw3zau,mvw3zzu,
     3	             mvw4aau,mvw4azu,mvw4zau,mvw4zzu,
     2		     mv1u,mv2u,mv3u,mv4u,
     2		     mv1wau,mv1wzu,mv2wau,mv2wzu,
     2		     mv3wau,mv3wzu,mv4wau,mv4wzu

      double complex psi(2,-1:1,4), jqq(0:5,-1:1,2) !, eps(0:5)
      double complex psiwe(2,-1:1,4),
     1 		     psizu(2,-1:1,4),psiau(2,-1:1,4),
     1               jew(0:5,-1:1,4),
     1		     juz(0:5,-1:1,4),jua(0:5,-1:1,4),
     1               jvew(0:5,-1:1,-1:1,4), 
     1               jvuz(0:5,-1:1,-1:1,4), jvua(0:5,-1:1,-1:1,4)
      double complex bkjqq(2,-1:1,-1:1,4)
      double complex jw(0:5,-1:1,-1:1,4), 
     1               jvw(0:5,-1:1,-1:1,4),
     2               psiw(2,-1:1,-1:1,4)
      double complex j5wz(0:3,3,-1:1,2), j5zw(0:3,3,-1:1,2),
     1 		     j5wa(0:3,3,-1:1,2), j5aw(0:3,3,-1:1,2)
      double complex epsaa(0:3),epsaz(0:3),epsza(0:3),epszz(0:3),epsww(0:3),
     1 		     epsaw(0:3),epszw(0:3),epswa(0:3),epswz(0:3),
     1		     epsawe(0:3,-1:1,3:4,2), epszwe(0:3,-1:1,3:4,2),
     1		     epswae(0:3,-1:1,3:4,2), epswze(0:3,-1:1,3:4,2),
     1		     epswee(0:3,-1:1,3:4,2), epswuu(0:3,-1:1,3:4,2),
     1		     epsauu(0:3,-1:1,3:4,2), epszuu(0:3,-1:1,3:4,2),
     1		     qepswee(3:4),qepszwe(3:4),qepswze(3:4),
     1		     qepszuu(3:4),qepswuu(3:4)	 
      double precision fqwe(0:4,4),fqzu(0:4,4),fqau(0:4,4)
      double precision fq(0:4,4),qee(0:4), quu(0:4), dummy(0:4)
      double complex zm2i(1:4)
      double complex jj21we(-1:1),jj43we(-1:1),jj21au(-1:1),jj43au(-1:1),
     1 		     jj21zu(-1:1),jj43zu(-1:1)

      double complex propt1(-1:1,-1:1,4,2), propt2(-1:1,-1:1,4,2),
     1		     propw(2),
     1		     prop_ee(3), prop_uu(3)
c
       double complex psiwg,jwg(2,-1:1,-1:1,3)
       integer kb,kp
c

      integer idu,idl
      double complex contract_Tjj, dotcc, dotrc, dotqj, s1c
      external contract_Tjj, dotcc, dotrc, dotqj, s1c
      logical ldebug, ldebugv, ldebugvf,linit,lerror,lgc(4)
      data linit /.true./, lerror /.false./, ldebug /.false./
      data lgc /4*.false./
      
      logical vvdebug,vbdebug, wwvdebug,bbdebug,pdebug
      parameter (vvdebug = .false.,vbdebug = .false.,wwvdebug =.false.,
     #		 bbdebug = .false.,pdebug = .false.)
      
      integer*8 icb1, icount1, icb2, icount2
      data icb1/0/,icount1/0/,icb2/0/,icount2/0/,Lold/1/
      double precision xgc1, xgc2
      double complex mvwz,mvzw,mvwa,mvaw,xvwz,xvzw,xvwa,xvaw
 
      parameter (ldebugv=.false., ldebugvf=.false.)
      save ifl, Lold, zm2i, linit, lgc, icb1, icb2, icount1, icount2
      save jw,jvw,psiw,fq !,j5pm,j5mp
      save j5wz,j5zw,j5wa,j5aw
      double complex  zero
      parameter (zero = (0d0,0d0) )

   
      integer icountmax

c      integer iterations1, iterations2, n2max1, n2max2
c      integer lflavr(5:4+max_v)
c      character*50 gridname1, gridname2
c      common /bmainpr/ iterations1, iterations2, n2max1, n2max2,
c     2                 lflavr, gridname1, gridname2
     
c
c variables for virtual corrections
c
	logical bad_gauge,bad_gauge_sin,bad_gauge_isig
	common /gauge / bad_gauge,bad_gauge_sin,bad_gauge_isig
c
      double precision c2,c2o4pi     !,pi2o3, cvirtc
c      parameter (pi2o3=pi**2/3d0, cvirtc=pi2o3-7d0)
      parameter (c2=4d0/3d0, c2o4pi=c2/4d0/pi)
      double precision lnmusqp1p2, lnmusqp3p4

      logical lnlo, lbox, lpent, lpt, lpq
c      parameter (lnlo=.true.)
      lnlo = NLO.ne.0    ! include some virtual stuff if T
      lbox = NLO.eq.1 .or. NLO.eq.-4 .or. NLO.eq.5    ! call boxline if T
      lpt  = NLO.eq.1 .or. NLO.le.-5                  ! true pentagon contributions on
      lpq  = NLO.eq.5 .or. NLO.eq.-4                  ! q^mu/m_V terms of Pentagon contributions only
      lpent= lpt .or. lpq                             ! include pentagons if T
c
	bad_gauge = .false. ! set F at beginning of each run
c
c
c-------------------------

	if (bos.ne.32) then
	    print*, "subroutine qqwpzqq contains only ME for W+ Z production"
	    print*, "but you called it with bos = ",bos
	    stop
	endif	
	
c-------------------------
c	
c
c define flavors of external quarks for the 4 subprocesses
c
      if (linit) then
         linit = .false.
         k = 1                  ! uucs
         ifl(1,k) = 3
         ifl(2,k) = 3
         ifl(3,k) = 3
         ifl(4,k) = 4
         k = 2                  ! ddcs
         ifl(1,k) = 4
         ifl(2,k) = 4
         ifl(3,k) = 3
         ifl(4,k) = 4
         k = 3                  ! udcc
         ifl(1,k) = 3
         ifl(2,k) = 4
         ifl(3,k) = 3
         ifl(4,k) = 3
         k = 4                  ! udss
         ifl(1,k) = 3
         ifl(2,k) = 4
         ifl(3,k) = 4
         ifl(4,k) = 4
c 	 
         zm2i(2) = 1/dcmplx(xm2(2),-xmg(2))
         zm2i(3) = 1/dcmplx(xm2(3),-xmg(3))
         zm2i(4) = 1/dcmplx(xm2(4),-xmg(4))	 
      endif

 701  if ((ldebug.or.ldebugv)) then
         print*," " 
         print*," enter qqwzqq "
      endif
      do k = 1,4
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
                  jvew(k,isig1,isig3,i) = 0
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
      if (ldebugv) then
         print*," fsign = ",(sign(i),i=1,4)
         print*," is1, is3 = ",is1,is3," js1, js3 = ",js1,js3
         print*," lws = ",lws," lws1 = ",lws(js1), " lws3 = ",lws(js3) 
         if (sign(1).ne.sign(2) .or. sign(3).ne.sign(4) ) then
            print*," quark sign factors ",(sign(i),i=1,4)
            print*," are not implemented yet in qqwzqq"
            stop
         endif
      endif

c
c define the internal momenta
c
      do mu = 0,3
         do i = 1,4+n_v
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo
	 
           p21(mu) =   p(mu,2) - p(mu,1)
           p43(mu) =   p(mu,4) - p(mu,3)
           pew(mu) =   p(mu,6) - p(mu,5)
           puz(mu) =   p(mu,8) - p(mu,7)
           pwz(mu) =   pew(mu) + puz(mu)	 	 
	 
      enddo
      p21(4) = p21(0)**2 - p21(1)**2 - p21(2)**2 - p21(3)**2
      p43(4) = p43(0)**2 - p43(1)**2 - p43(2)**2 - p43(3)**2
      pew(4) = pew(0)**2 - pew(1)**2 - pew(2)**2 - pew(3)**2
      puz(4) = puz(0)**2 - puz(1)**2 - puz(2)**2 - puz(3)**2
      pwz(4) = pwz(0)**2 - pwz(1)**2 - pwz(2)**2 - pwz(3)**2
	 

      if (lwarn) then
         do mu=0,3
            dummy(mu) = p(mu,1)+p(mu,3)+p(mu,5)+p(mu,7)-
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
            do i = 1,8
               write(6,101) i, p(0,i), p(1,i), p(2,i), p(3,i)
            end do
         endif 
      endif
c

c	print*,"pwz=",pwz(0)

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
         qee(mu) = pew(mu)+p(mu,2) - p(mu,1)    ! e ve  emitted on upper line
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
      do isig = -1,1,2
            maw = contract_Tjj(awwv(0,0,L),jqq(0,isig,1),jqq(0,-1,2))
            mzw = contract_Tjj(zwwv(0,0,L),jqq(0,isig,1),jqq(0,-1,2))
            mwa = contract_Tjj(wawv(0,0,L),jqq(0,isig,2),jqq(0,-1,1))
            mwz = contract_Tjj(wzwv(0,0,L),jqq(0,isig,2),jqq(0,-1,1))
            do k = 1,2
               mat(k,isig,-1,1) = 
     1              maw*clr(ifl(1,k),1,isig)*clr(ifl(3,k),3,-1)
     2            + mzw*clr(ifl(1,k),2,isig)*clr(ifl(3,k),3,-1)
               matv(k,isig,-1,1) = (0d0,0d0)
            enddo
            do k = 3,4
               mat(k,-1,isig,1) = 
     3              mwa*clr(ifl(1,k),3,-1)*clr(ifl(3,k),1,isig)
     4            + mwz*clr(ifl(1,k),3,-1)*clr(ifl(3,k),2,isig)
               matv(k,-1,isig,1) = (0d0,0d0)
            enddo
	   	   

      enddo !isig
	    
            if (vvdebug) then
	    print*,"start loop:"
	    do isig1 = -1,1,2
	    do isig3 = -1,1,2
               do k = 1,4
                  m1 = 0
                  do i = 116,151
                     m1 = m1 + tampw(i,isig1,isig3,k)
                  enddo
                  do i = 168,191
                     m1 = m1 + tampw(i,isig1,isig3,k)
                  enddo
                  do i = 200,210
                     m1 = m1 + tampw(i,isig1,isig3,k)
                  enddo
c		  print*,"fsign=",sign
                  print*," vv for k=",k, ", isig=",isig1,isig3
c		  if (m1.ne.0d0) then 
		  if (k.lt.3) then
                  print*," res:",mat(k,isig1,isig3,1),m1
                  print*," rat:",mat(k,isig1,isig3,1)/m1
		  else !for crossed processes
                  print*," res:",mat(k,isig3,isig1,1),m1
                  print*," rat:",mat(k,isig3,isig1,1)/m1
		  endif
		  print*
c		  endif
              enddo !k
              enddo !isig3
              enddo !isig1	
            endif !vvdebug
      
c -------------------------------------------------------------------
c -------------------------------------------------------------------
c
c  prepare box diagrams: attach A/Z to external spinors
c 
c      isig = +-1   : left- and righthanded spinors coupling to A/Z"s
c
c  Notation for virtual 2-component spinors and momenta
c
c  W+ -> e+ve   attached to quark number i: psiwe(*,isig,i) with momentum fqwe(mu,i)
c  Z  -> mu+mu- attached to quark number i: psizu(*,isig,i) with momentum fqzu(mu,i)
c  A  -> mu+mu- attached to quark number i: psiau(*,isig,i) with momentum fqau(mu,i)
c  
c  the fermion current corresponding to a quark line with the emitted V-> l+l- or W->l nu
c  attached next to quark number i is stored in jlv(mu,isig,i) etc. 
c 
c  For the virtual amlitudes the notation, e.g.     jvlv(mu,isig,is,i)
c  is used for the boxline correction to a quark line with one V attached next to
c  quark #i and a free Lorentz index mu for the second attached EW boson. is=+-1
c  refers to the sign factor of this quark (vs.antiquark line). They are recalculated 
c  only if this quark line sign has not been calculated yet for this phase space point
c  (i.e. lws = .true.)  Otherwise they are taken from saved previous calculation
 
      do i = 1,3,2

c W couples only to left-handed fermions:
	 isig = -1
         call ket2c(psi(1,isig,i),.true.,p(0,i),isig,qe,we,
     1              psiwe(1,isig,i),fqwe(0,i))
	 call bra2c(psi(1,isig,i+1),.true.,p(0,i+1),isig,qe,we,
     1              psiwe(1,isig,i+1),fqwe(0,i+1))
     
         call curr6(-1,psi(1,-1,i+1),p(0,i+1),
     1                 psiwe(1,-1,i),fqwe(0,i), jew(0,-1,i)   )
         call curr6(-1,psiwe(1,-1,i+1),fqwe(0,i+1),
     1                 psi(1,-1,i),p(0,i), jew(0,-1,i+1) )
     

c A/Z couple also to right-handed fermions:
        do isig = -1,1,2
      
         call ket2c(psi(1,isig,i),.true.,p(0,i),isig,qu,zu,
     1              psizu(1,isig,i),fqzu(0,i))
         call ket2c(psi(1,isig,i),.true.,p(0,i),isig,qu,au,
     1              psiau(1,isig,i),fqau(0,i))

         call bra2c(psi(1,isig,i+1),.true.,p(0,i+1),isig,qu,zu,
     1              psizu(1,isig,i+1),fqzu(0,i+1))
         call bra2c(psi(1,isig,i+1),.true.,p(0,i+1),isig,qu,au,
     1              psiau(1,isig,i+1),fqau(0,i+1))

         enddo !isig
        
         call curr6(1,psi(1,-1,i+1),p(0,i+1),
     1                 psizu(1,-1,i),fqzu(0,i), juz(0,-1,i)   )
         call curr6(1,psi(1,-1,i+1),p(0,i+1),
     1                 psiau(1,-1,i),fqau(0,i), jua(0,-1,i)   )

         call curr6(1,psizu(1,-1,i+1),fqzu(0,i+1),
     1                 psi(1,-1,i),p(0,i), juz(0,-1,i+1) )
         call curr6(1,psiau(1,-1,i+1),fqau(0,i+1),
     1                 psi(1,-1,i),p(0,i), jua(0,-1,i+1) )
     
   
       enddo !i

#ifdef WITH_NLO                    
      if (lbox .and. lws(js1)) then
         do i = 1,2
            call boxlinec(-1,psi(1,-1,1),psi(1,-1,2),p(0,1),p(0,2),
     1                    .true., we,3-i,   jew(0,-1,i), 
     2                    jvew(0,-1,is1,i) )
            call boxlinec(1,psi(1,-1,1),psi(1,-1,2),p(0,1),p(0,2),
     1                    .true., zu,3-i,   juz(0,-1,i), 
     2                    jvuz(0,-1,is1,i) )	
            call boxlinec(1,psi(1,-1,1),psi(1,-1,2),p(0,1),p(0,2),
     1                    .true., au,3-i,   jua(0,-1,i), 
     2                    jvua(0,-1,is1,i) )
	      
          enddo !i
      endif !lbox
      
      if (lbox .and. lws(js3)) then
         do i = 3,4	 
            call boxlinec(-1,psi(1,-1,3),psi(1,-1,4),p(0,3),p(0,4),
     1                    .true., we,5-i,   jew(0,-1,i), 
     2                    jvew(0,-1,is3,i) )
            call boxlinec(1,psi(1,-1,3),psi(1,-1,4),p(0,3),p(0,4),
     1                    .true., zu,5-i,   juz(0,-1,i), 
     2                    jvuz(0,-1,is3,i) )
            call boxlinec(1,psi(1,-1,3),psi(1,-1,4),p(0,3),p(0,4),
     1                    .true., au,5-i,   jua(0,-1,i), 
     2                    jvua(0,-1,is3,i) )
        
         enddo !i
	 
      endif !lbox
#endif
c
c
c -------------------------------------------------------------------
c -------------------------------------------------------------------
c
c now calculate the Vertex-box diagrams; 
      do k = 1,2
         kk = 3-k
c k=1 and kk=2 is for "box correction" to upper line
c k=2 and kk=1 is for "box correction" to lower line

	do isig = -1,1,2
	
	! VW-> e nu: (V = A or Z)
	
	  if (k.eq.1) then
	    call contract_T1j(awen(0,0,k,L),jqq(0,isig,kk),epsaw) !for AW
	    call contract_T1j(zwen(0,0,k,L),jqq(0,isig,kk),epszw) !for ZW
	    call contract_T2j(waen(0,0,k,L),jqq(0,isig,kk),epswa) !for WA
	    call contract_T2j(wzen(0,0,k,L),jqq(0,isig,kk),epswz) !for WZ
	  else
 	    call contract_T2j(awen(0,0,k,L),jqq(0,isig,kk),epsaw) !for AW
	    call contract_T2j(zwen(0,0,k,L),jqq(0,isig,kk),epszw) !for ZW	    
	    call contract_T1j(waen(0,0,k,L),jqq(0,isig,kk),epswa) !for WA
	    call contract_T1j(wzen(0,0,k,L),jqq(0,isig,kk),epswz) !for WZ
	  endif
	  
	! VV-> mu+mu-:

	  if (k.eq.1) then
	    call contract_T2j(aauu(0,0,k,L),jqq(0,isig,kk),epsaa) !for AA 
	    call contract_T2j(azuu(0,0,k,L),jqq(0,isig,kk),epsaz) !for AZ
	    call contract_T1j(zauu(0,0,k,L),jqq(0,isig,kk),epsza) !for ZA
	    call contract_T2j(zzuu(0,0,k,L),jqq(0,isig,kk),epszz) !for ZZ 
	    call contract_T1j(wwuu(0,0,k,L),jqq(0,isig,kk),epsww) !for WW 
	  else
	    call contract_T1j(aauu(0,0,k,L),jqq(0,isig,kk),epsaa) !for AA 
	    call contract_T1j(azuu(0,0,k,L),jqq(0,isig,kk),epsza) !for AZ 
	    call contract_T2j(zauu(0,0,k,L),jqq(0,isig,kk),epsaz) !for ZA 
	    call contract_T1j(zzuu(0,0,k,L),jqq(0,isig,kk),epszz) !for ZZ 
	    call contract_T1j(wwuu(0,0,k,L),jqq(0,isig,kk),epsww) !for WW 
	  endif
	
          do mu = 0,3 ! graphs with uncontracted V from VV"->mu+mu-
	  	      ! (epsvuu) or from VV->e nu (epsvwe/epswve/epswee) 
		      ! always: k=3/4 for quark type
		      
	     epsauu(mu,isig,3,k) =
     1  	 epsaa(mu)*clr(3,1,isig)+epsaz(mu)*clr(3,2,isig) 
	     epsauu(mu,isig,4,k) =
     1  	 epsaa(mu)*clr(4,1,isig)+epsaz(mu)*clr(4,2,isig)  
             epszuu(mu,isig,3,k) =
     1   	  epsza(mu)*clr(3,1,isig)+epszz(mu)*clr(3,2,isig)
             epszuu(mu,isig,4,k) =
     1   	  epsza(mu)*clr(4,1,isig)+epszz(mu)*clr(4,2,isig)
             
	     epswuu(mu,isig,3,k) = epsww(mu)*clr(3,3,isig) !k=3/4 for quark type
	     epswuu(mu,isig,4,k) = epswuu(mu,isig,3,k)
	     
             epsawe(mu,isig,3,k) = epsaw(mu)*clr(3,3,isig)  
             epszwe(mu,isig,3,k) = epszw(mu)*clr(3,3,isig)	       
             epsawe(mu,isig,4,k) = epsawe(mu,isig,3,k)  
             epszwe(mu,isig,4,k) = epszwe(mu,isig,3,k)
	     	     
             epswae(mu,isig,3,k) = epswa(mu)*clr(3,3,isig)  
             epswze(mu,isig,3,k) = epswz(mu)*clr(3,3,isig)	     
             epswae(mu,isig,4,k) = epswae(mu,isig,3,k)  
             epswze(mu,isig,4,k) = epswze(mu,isig,3,k)
	     	     	    
	    if (k.eq.2) then 
	     epswee(mu,isig,3,k) =
     1   	  epsaw(mu)*clr(3,1,isig)+epszw(mu)*clr(3,2,isig)  
	     epswee(mu,isig,4,k) =
     1   	  epsaw(mu)*clr(4,1,isig)+epszw(mu)*clr(4,2,isig)  
	    else !k=1
	     epswee(mu,isig,3,k) =
     1  	  epswa(mu)*clr(3,1,isig)+epswz(mu)*clr(3,2,isig)  
	     epswee(mu,isig,4,k) =
     1  	  epswa(mu)*clr(4,1,isig)+epswz(mu)*clr(4,2,isig)  

	     endif

          enddo
	
	enddo !isig
	 
      enddo !k
            	
      do isig = -1,1,2 	
            
      jj21au(isig) = dotcc(au,jqq(0,isig,1))
      jj21zu(isig) = dotcc(zu,jqq(0,isig,1))
      jj21we(isig) = dotcc(we,jqq(0,isig,1))

      jj43au(isig) = dotcc(au,jqq(0,isig,2))
      jj43zu(isig) = dotcc(zu,jqq(0,isig,2))
      jj43we(isig) = dotcc(we,jqq(0,isig,2))

      enddo !isig
c
c --------------
c
c now construct the contribution to the amplitude by current contraction 
c virtual contributions are assembled in subroutine boxline
                ! 2 bosons attached to 12 line     
c

	do isig1 = -1,1,2
	   do isig3 = -1,1,2
	   
c  gauge term for Z/W propagator:
	  do k=3,4 !up/down type
	      qepszwe(k) = -dotrc(quu,epszwe(0,isig3,k,1))*zm2i(2) 	      
	      qepswee(k) = -dotrc(quu,epswee(0,isig3,k,1))*zm2i(3) 
	      qepszuu(k) = -dotrc(qee,epszuu(0,isig3,k,1))*zm2i(2) 	    
	      qepswuu(k) = -dotrc(qee,epswuu(0,isig3,k,1))*zm2i(3) 	    
	   enddo !k
	   
	   do k = 1,4 ! e+ ve from upper line 
	   
	   if (isig1.eq.-1) then
	   if (k.eq.2) then  
	   	   
	     m1we = dotcc(jew(0,isig1,2),epswuu(0,isig3,ifl(3,k),1))	   
	     w1we = qepswuu(ifl(3,k))*jj21we(isig1)    
	   
	     m1e = (m1we+w1we)*clr(ifl(1,k),3,isig1)*clr(ifl(2,k),4,isig1)
	     m2e = 0d0

	   elseif (k.eq.1) then  
	   	   
	     m2we = dotcc(jew(0,isig1,1),epswuu(0,isig3,ifl(3,k),1))	   
	     w2we = qepswuu(ifl(3,k))*jj21we(isig1)    
	   
	     m2e = (m2we-w2we)*clr(ifl(1,k),3,isig1)*clr(ifl(2,k),4,isig1)
	     m1e = 0d0
	     
	   else ! k=3,4 			     
	   
	     m1ae = dotcc(jew(0,isig1,2),epsauu(0,isig3,ifl(3,k),1))
	     m1ze = dotcc(jew(0,isig1,2),epszuu(0,isig3,ifl(3,k),1))
	   
	     z1ze = qepszuu(ifl(3,k))*jj21we(isig1)    
	   
	     mw1ae =  m1ae*clr(ifl(1,k),1,isig1)*clr(ifl(2,k),3,isig1)
	     mw1ze = (m1ze+z1ze)*clr(ifl(1,k),2,isig1)*clr(ifl(2,k),3,isig1)
	   	   
	     m1e = mw1ae+mw1ze
	   	   	   
	     m2ae = dotcc(jew(0,isig1,1),epsauu(0,isig3,ifl(3,k),1))
	     m2ze = dotcc(jew(0,isig1,1),epszuu(0,isig3,ifl(3,k),1))
	   	   
	     mw2ae =  m2ae*clr(ifl(2,k),1,isig1)*clr(ifl(1,k),3,isig1)
	     mw2ze = (m2ze-z1ze)*clr(ifl(2,k),2,isig1)*clr(ifl(1,k),3,isig1)
	   	   
	     m2e = mw2ae+mw2ze
	   	   
	   endif !k
	   
	   else
	     m1e  =  0d0
	     m2e  =  0d0
	   endif !isig1
	      
	   
	   ! mu+mu- from upper line:
	   
	   if (k.lt.3.and.isig3.eq.-1) then !k=1,2
	     m1aau = dotcc(jua(0,isig1,2),epsawe(0,isig3,ifl(3,k),1))
	     m1azu = dotcc(jua(0,isig1,2),epszwe(0,isig3,ifl(3,k),1))
	     m1zau = dotcc(juz(0,isig1,2),epsawe(0,isig3,ifl(3,k),1))
	     m1zzu = dotcc(juz(0,isig1,2),epszwe(0,isig3,ifl(3,k),1))

	     z1azu = qepszwe(ifl(3,k))*jj21au(isig1)  
	     z1zzu = qepszwe(ifl(3,k))*jj21zu(isig1)  
	   
	     mw1aau =  m1aau*clr(ifl(1,k),1,isig1)*clr(ifl(2,k),1,isig1)
	     mw1azu = (m1azu+z1azu)*clr(ifl(1,k),2,isig1)*clr(ifl(2,k),1,isig1)
  	     mw1zau =  m1zau*clr(ifl(1,k),1,isig1)*clr(ifl(2,k),2,isig1)
	     mw1zzu = (m1zzu+z1zzu)*clr(ifl(1,k),2,isig1)*clr(ifl(2,k),2,isig1)
	   	   
	     m1u = mw1aau+mw1azu+mw1zau+mw1zzu
	   
	     m2aau = dotcc(jua(0,isig1,1),epsawe(0,isig3,ifl(3,k),1))
	     m2azu = dotcc(jua(0,isig1,1),epszwe(0,isig3,ifl(3,k),1))
	     m2zau = dotcc(juz(0,isig1,1),epsawe(0,isig3,ifl(3,k),1))
	     m2zzu = dotcc(juz(0,isig1,1),epszwe(0,isig3,ifl(3,k),1))
	   	   
	     mw2aau =  m2aau*clr(ifl(2,k),1,isig1)*clr(ifl(1,k),1,isig1)
	     mw2azu = (m2azu-z1azu)*clr(ifl(2,k),2,isig1)*clr(ifl(1,k),1,isig1)
	     mw2zau =  m2zau*clr(ifl(2,k),1,isig1)*clr(ifl(1,k),2,isig1)
	     mw2zzu = (m2zzu-z1zzu)*clr(ifl(2,k),2,isig1)*clr(ifl(1,k),2,isig1)
	   	   
	     m2u = mw2aau+mw2azu+mw2zau+mw2zzu
	   
	   elseif (k.gt.2) then	!k=3,4
	     m1awu = dotcc(jua(0,isig1,2),epswee(0,isig3,ifl(3,k),1))
	     m1zwu = dotcc(juz(0,isig1,2),epswee(0,isig3,ifl(3,k),1))
	   
	     z1awu = qepswee(ifl(3,k))*jj21au(isig1)  
	     z1zwu = qepswee(ifl(3,k))*jj21zu(isig1)  
	   
	     m1wau = (m1awu+z1awu)*clr(ifl(1,k),3,isig1)*clr(ifl(2,k),1,isig1)
	     m1wzu = (m1zwu+z1zwu)*clr(ifl(1,k),3,isig1)*clr(ifl(2,k),2,isig1)
	   	   
	     m1u = m1wau+m1wzu
	   	   
	     m2awu = dotcc(jua(0,isig1,1),epswee(0,isig3,ifl(3,k),1))
	     m2zwu = dotcc(juz(0,isig1,1),epswee(0,isig3,ifl(3,k),1))
	   	   
	     m2wau = (m2awu-z1awu)*clr(ifl(2,k),3,isig1)*clr(ifl(1,k),1,isig1)
	     m2wzu = (m2zwu-z1zwu)*clr(ifl(2,k),3,isig1)*clr(ifl(1,k),2,isig1)
	   	   
	     m2u = m2wau+m2wzu
	  
	   else
	     m1u = 0d0
	     m2u = 0d0	  
	   endif !k
	   
	   mat(k,isig1,isig3,2) =  m1e-m1u+m2e-m2u
	   
	   
	   if (lbox) then
	   
	    ! e+ve from upper line:
	    	   
	   if (isig1.eq.-1) then
	   if (k.eq.2) then  
	   	   
	     mv1e = dotcc(jvew(0,isig1,is1,2),epswuu(0,isig3,ifl(3,k),1))
     #		  *clr(ifl(1,k),3,isig1)*clr(ifl(2,k),4,isig1)
	     mv2e = 0d0

	   elseif (k.eq.1) then  
	   	   
	     mv2e = dotcc(jvew(0,isig1,is1,1),epswuu(0,isig3,ifl(3,k),1))
     #		  *clr(ifl(1,k),3,isig1)*clr(ifl(2,k),4,isig1)
	     mv1e = 0d0
	   
	   else ! k=3,4
	  
	     mvw1ae =  dotcc(jvew(0,isig1,is1,2),epsauu(0,isig3,ifl(3,k),1))
     #	   	    *clr(ifl(1,k),1,isig1)*clr(ifl(2,k),3,isig1)
	     mvw1ze =  dotcc(jvew(0,isig1,is1,2),epszuu(0,isig3,ifl(3,k),1))
     #	            *clr(ifl(1,k),2,isig1)*clr(ifl(2,k),3,isig1)
	   	   
	     mv1e = mvw1ae+mvw1ze
	   	   	   	   	   
	     mvw2ae =  dotcc(jvew(0,isig1,is1,1),epsauu(0,isig3,ifl(3,k),1))
     #	            *clr(ifl(2,k),1,isig1)*clr(ifl(1,k),3,isig1)
	     mvw2ze = dotcc(jvew(0,isig1,is1,1),epszuu(0,isig3,ifl(3,k),1))
     #	   	    *clr(ifl(2,k),2,isig1)*clr(ifl(1,k),3,isig1)
	   	   
	     mv2e = mvw2ae+mvw2ze
	   
	   endif !k
	   
	   else
	     mv1e  =  0d0
	     mv2e  =  0d0
	   endif !isig1
	   	   
	   ! mu+mu- from upper line:	   
	   
	   if (k.lt.3.and.isig3.eq.-1) then !k=1,2
	   
	     mvw1aau = dotcc(jvua(0,isig1,is1,2),epsawe(0,isig3,ifl(3,k),1))
     #	   	    *clr(ifl(1,k),1,isig1)*clr(ifl(2,k),1,isig1)
	     mvw1azu = dotcc(jvua(0,isig1,is1,2),epszwe(0,isig3,ifl(3,k),1))
     #	   	    *clr(ifl(1,k),2,isig1)*clr(ifl(2,k),1,isig1)
	     mvw1zau = dotcc(jvuz(0,isig1,is1,2),epsawe(0,isig3,ifl(3,k),1))
     #	   	    *clr(ifl(1,k),1,isig1)*clr(ifl(2,k),2,isig1)
	     mvw1zzu = dotcc(jvuz(0,isig1,is1,2),epszwe(0,isig3,ifl(3,k),1))
     #	   	    *clr(ifl(1,k),2,isig1)*clr(ifl(2,k),2,isig1)
	   	   
	     mv1u = mvw1aau+mvw1azu+mvw1zau+mvw1zzu
	   	   	   
	     mvw2aau = dotcc(jvua(0,isig1,is1,1),epsawe(0,isig3,ifl(3,k),1))
     #	   	    *clr(ifl(2,k),1,isig1)*clr(ifl(1,k),1,isig1)
	     mvw2azu = dotcc(jvua(0,isig1,is1,1),epszwe(0,isig3,ifl(3,k),1))
     #	   	    *clr(ifl(2,k),2,isig1)*clr(ifl(1,k),1,isig1)
	     mvw2zau = dotcc(jvuz(0,isig1,is1,1),epsawe(0,isig3,ifl(3,k),1))
     #	   	    *clr(ifl(2,k),1,isig1)*clr(ifl(1,k),2,isig1)
	     mvw2zzu = dotcc(jvuz(0,isig1,is1,1),epszwe(0,isig3,ifl(3,k),1))
     #	   	    *clr(ifl(2,k),2,isig1)*clr(ifl(1,k),2,isig1)
	   	   
	     mv2u = mvw2aau+mvw2azu+mvw2zau+mvw2zzu


	   elseif (k.gt.2) then		!k=3,4
	   
	     mv1wau = dotcc(jvua(0,isig1,is1,2),epswee(0,isig3,ifl(3,k),1))
     #	   	    *clr(ifl(1,k),3,isig1)*clr(ifl(2,k),1,isig1)
	     mv1wzu = dotcc(jvuz(0,isig1,is1,2),epswee(0,isig3,ifl(3,k),1))
     #	   	    *clr(ifl(1,k),3,isig1)*clr(ifl(2,k),2,isig1)
	   	   
	     mv1u = mv1wau+mv1wzu
	   	   
	   	   
	     mv2wau = dotcc(jvua(0,isig1,is1,1),epswee(0,isig3,ifl(3,k),1))
     #	   	    *clr(ifl(2,k),3,isig1)*clr(ifl(1,k),1,isig1)
	     mv2wzu = dotcc(jvuz(0,isig1,is1,1),epswee(0,isig3,ifl(3,k),1))
     #	   	    *clr(ifl(2,k),3,isig1)*clr(ifl(1,k),2,isig1)
	   	   
	     mv2u = mv2wau+mv2wzu
	
	  
	   else
	     mv1u = 0d0
	     mv2u = 0d0	  
	   endif !k
	   
	   matv(k,isig1,isig3,2) = mv1e-mv1u+mv2e-mv2u
	   
	   endif !lbox
	   
	   enddo !k
	   
	   enddo !isig3
	enddo ! isig1


c------------------------------------------------------------
      
c repeat the same for 2 bosons attached to 34 line:
c
c
	do isig1 = -1,1,2
	   do isig3 = -1,1,2
	   
c  gauge term for Z/W propagator:
	  do k=3,4 !up/down type
	      qepswze(k) = dotrc(qee,epswze(0,isig1,k,2))*zm2i(2) 
	      qepswee(k) = dotrc(qee,epswee(0,isig1,k,2))*zm2i(3) 
	      qepszuu(k) = dotrc(quu,epszuu(0,isig1,k,2))*zm2i(2) 	    
	      qepswuu(k) = dotrc(quu,epswuu(0,isig1,k,2))*zm2i(3) 	    
	   enddo !k
	   
	   do k = 1,4 ! e+ve from lower line:
	   
	   if(isig3.eq.-1) then
	   if (k.lt.3) then !k=1,2
	   
	     m3ae = dotcc(jew(0,isig3,4),epsauu(0,isig1,ifl(1,k),2))
	     m3ze = dotcc(jew(0,isig3,4),epszuu(0,isig1,ifl(1,k),2))
	   
	     z3ze = qepszuu(ifl(1,k))*jj43we(isig3)    
	   
	     mw3ae =  m3ae*clr(ifl(3,k),1,isig3)*clr(ifl(4,k),3,isig3)
	     mw3ze = (m3ze+z3ze)*clr(ifl(3,k),2,isig3)*clr(ifl(4,k),3,isig3)
	   	   
	     m3e = mw3ae+mw3ze
	   	   	   
	     m4ae = dotcc(jew(0,isig3,3),epsauu(0,isig1,ifl(1,k),2))
	     m4ze = dotcc(jew(0,isig3,3),epszuu(0,isig1,ifl(1,k),2))
	   	   
	     mw4ae =  m4ae*clr(ifl(4,k),1,isig3)*clr(ifl(3,k),3,isig3)
	     mw4ze = (m4ze-z3ze)*clr(ifl(4,k),2,isig3)*clr(ifl(3,k),3,isig3)
	   	   
	     m4e = mw4ae+mw4ze
	   
	   elseif (k.eq.3) then
	   
	     m3e = 0d0

	     m4we = dotcc(jew(0,isig3,3),epswuu(0,isig1,ifl(1,k),2))
	     w4we = qepswuu(ifl(1,k))*jj43we(isig3)    	   	   	   	   	   
	     m4e = (m4we-w4we)*clr(ifl(4,k),4,isig3)*clr(ifl(3,k),3,isig3)

	   elseif (k.eq.4) then
	   
	     m4e = 0d0

	     m3we = dotcc(jew(0,isig3,4),epswuu(0,isig1,ifl(1,k),2))
	     w3we = qepswuu(ifl(1,k))*jj43we(isig3)    	   	   	   	   	   
	     m3e = (m3we+w3we)*clr(ifl(4,k),3,isig3)*clr(ifl(3,k),3,isig3)
	   	   
	   endif
	   
	   else
	     m3e = 0d0
	     m4e = 0d0
	   endif !isig3  
	   
	   	   
	   ! mu+mu- from lower line:
	   
	   if (k.lt.3.and.isig3.eq.-1) then	! k =1,2		
	     m3awu = dotcc(jua(0,isig3,4),epswee(0,isig1,ifl(1,k),2))
	     m3zwu = dotcc(juz(0,isig3,4),epswee(0,isig1,ifl(1,k),2))
	   
	     z3awu = qepswee(ifl(1,k))*jj43au(isig3)  
	     z3zwu = qepswee(ifl(1,k))*jj43zu(isig3)  
	   
	     m3wau = (m3awu+z3awu)*clr(ifl(3,k),3,isig3)*clr(ifl(4,k),1,isig3)
	     m3wzu = (m3zwu+z3zwu)*clr(ifl(3,k),3,isig3)*clr(ifl(4,k),2,isig3)
	          
	     m3u = m3wau+m3wzu
	   	   
	     m4awu = dotcc(jua(0,isig3,3),epswee(0,isig1,ifl(1,k),2))
	     m4zwu = dotcc(juz(0,isig3,3),epswee(0,isig1,ifl(1,k),2))
	   	   
	     m4wau = (m4awu-z3awu)*clr(ifl(4,k),3,isig3)*clr(ifl(3,k),1,isig3)
	     m4wzu = (m4zwu-z3zwu)*clr(ifl(4,k),3,isig3)*clr(ifl(3,k),2,isig3)
	     	   
	     m4u = m4wau+m4wzu
	
	   elseif (k.gt.2) then !k=3,4
	     m3aau = dotcc(jua(0,isig3,4),epswae(0,isig1,ifl(1,k),2))
	     m3azu = dotcc(jua(0,isig3,4),epswze(0,isig1,ifl(1,k),2))
	     m3zau = dotcc(juz(0,isig3,4),epswae(0,isig1,ifl(1,k),2))
	     m3zzu = dotcc(juz(0,isig3,4),epswze(0,isig1,ifl(1,k),2))

	     z3azu = qepswze(ifl(1,k))*jj43au(isig3)  
	     z3zzu = qepswze(ifl(1,k))*jj43zu(isig3)  
	   
	     mw3aau =  m3aau*clr(ifl(3,k),1,isig3)*clr(ifl(4,k),1,isig3)
	     mw3azu = (m3azu+z3azu)*clr(ifl(3,k),2,isig3)*clr(ifl(4,k),1,isig3)
	     mw3zau =  m3zau*clr(ifl(3,k),1,isig3)*clr(ifl(4,k),2,isig3)
	     mw3zzu = (m3zzu+z3zzu)*clr(ifl(3,k),2,isig3)*clr(ifl(4,k),2,isig3)
	    	   
	     m3u = mw3aau+mw3azu+mw3zau+mw3zzu
	   
	     m4aau = dotcc(jua(0,isig3,3),epswae(0,isig1,ifl(1,k),2))
	     m4azu = dotcc(jua(0,isig3,3),epswze(0,isig1,ifl(1,k),2))
	     m4zau = dotcc(juz(0,isig3,3),epswae(0,isig1,ifl(1,k),2))
	     m4zzu = dotcc(juz(0,isig3,3),epswze(0,isig1,ifl(1,k),2))
	   	   
	     mw4aau =  m4aau*clr(ifl(4,k),1,isig3)*clr(ifl(3,k),1,isig3)
	     mw4azu = (m4azu-z3azu)*clr(ifl(4,k),2,isig3)*clr(ifl(3,k),1,isig3)
	     mw4zau =  m4zau*clr(ifl(4,k),1,isig3)*clr(ifl(3,k),2,isig3)
	     mw4zzu = (m4zzu-z3zzu)*clr(ifl(4,k),2,isig3)*clr(ifl(3,k),2,isig3)
	    	   
	     m4u = mw4aau+mw4azu+mw4zau+mw4zzu
	   
	   else
	     m3u = 0d0
	     m4u = 0d0
	   endif !k
	   
	   mat(k,isig1,isig3,3) = m3e-m3u+m4e-m4u
	   
	   enddo !k
	   
	   enddo !isig3
	enddo ! isig1


c  virtual corrections:  
c	(gauge terms for massive gauge boson propagators 
c		   not needed when contracted with boxes)
c
c
	if (lbox) then

	do isig1 = -1,1,2
	   do isig3 = -1,1,2
	   	   
	   do k = 1,4 ! e+e- from lower line:
	   	   	   
	   if(isig3.eq.-1) then
	   if (k.lt.3) then  !k=1,2	   
		   
	     mvw3ae = dotcc(jvew(0,isig3,is3,4),epsauu(0,isig1,ifl(1,k),2))
     #		     *clr(ifl(3,k),1,isig3)*clr(ifl(4,k),3,isig3)
	     mvw3ze = dotcc(jvew(0,isig3,is3,4),epszuu(0,isig1,ifl(1,k),2))
     #		     *clr(ifl(3,k),2,isig3)*clr(ifl(4,k),3,isig3)
	   	   
	     mv3e = mvw3ae+mvw3ze
	   	   	   	   	   
	     mvw4ae = dotcc(jvew(0,isig3,is3,3),epsauu(0,isig1,ifl(1,k),2))
     #		     *clr(ifl(4,k),1,isig3)*clr(ifl(3,k),3,isig3)
	     mvw4ze = dotcc(jvew(0,isig3,is3,3),epszuu(0,isig1,ifl(1,k),2))
     #		     *clr(ifl(4,k),2,isig3)*clr(ifl(3,k),3,isig3)
	   	    
	     mv4e = mvw4ae+mvw4ze
	   
	   elseif (k.eq.3) then
	   
	     mv3e = 0d0

	     mv4e = dotcc(jvew(0,isig3,is3,3),epswuu(0,isig1,ifl(1,k),2))
     #		   *clr(ifl(4,k),4,isig3)*clr(ifl(3,k),3,isig3)
	   
	   elseif (k.eq.4) then
	   
	     mv4e = 0d0

	     mv3e = dotcc(jvew(0,isig3,is3,4),epswuu(0,isig1,ifl(1,k),2))
     #	   	   *clr(ifl(4,k),3,isig3)*clr(ifl(3,k),3,isig3)
	   	   	   
	   endif !k
	   
	   else
	     mv3e = 0d0
	     mv4e = 0d0
	   endif !isig3  
	   
	   	   
	   ! mu+mu- from lower line:
	   
	   if (k.lt.3.and.isig3.eq.-1) then !k=1,2
	   
	     mv3wau = dotcc(jvua(0,isig3,is3,4),epswee(0,isig1,ifl(1,k),2))
     #		     *clr(ifl(3,k),3,isig3)*clr(ifl(4,k),1,isig3)
	     mv3wzu = dotcc(jvuz(0,isig3,is3,4),epswee(0,isig1,ifl(1,k),2))
     #	  	     *clr(ifl(3,k),3,isig3)*clr(ifl(4,k),2,isig3)
	          
	     mv3u = mv3wau+mv3wzu
	   	   	   	   
	     mv4wau = dotcc(jvua(0,isig3,is3,3),epswee(0,isig1,ifl(1,k),2))
     #		     *clr(ifl(4,k),3,isig3)*clr(ifl(3,k),1,isig3)
	     mv4wzu = dotcc(jvuz(0,isig3,is3,3),epswee(0,isig1,ifl(1,k),2))
     #		     *clr(ifl(4,k),3,isig3)*clr(ifl(3,k),2,isig3)
	   	   
	     mv4u = mv4wau+mv4wzu
	
	   elseif (k.gt.2) then !k=3,4
	   
c	     mvw3aau = dotcc(jvua(0,isig3,is3,4),epsawe(0,isig1,ifl(1,k),2))
c     #		     *clr(ifl(3,k),1,isig3)*clr(ifl(4,k),1,isig3)
	     mvw3aau = dotcc(jvua(0,isig3,is3,4),epswae(0,isig1,ifl(1,k),2))
     #		     *clr(ifl(3,k),1,isig3)*clr(ifl(4,k),1,isig3)
	     mvw3azu = dotcc(jvua(0,isig3,is3,4),epswze(0,isig1,ifl(1,k),2))
     #		     *clr(ifl(3,k),2,isig3)*clr(ifl(4,k),1,isig3)
	     mvw3zau = dotcc(jvuz(0,isig3,is3,4),epswae(0,isig1,ifl(1,k),2))
     #		     *clr(ifl(3,k),1,isig3)*clr(ifl(4,k),2,isig3)
	     mvw3zzu = dotcc(jvuz(0,isig3,is3,4),epswze(0,isig1,ifl(1,k),2))
     #		     *clr(ifl(3,k),2,isig3)*clr(ifl(4,k),2,isig3)
	   	   
	     mv3u = mvw3aau+mvw3azu+mvw3zau+mvw3zzu
	   
	     mvw4aau = dotcc(jvua(0,isig3,is3,3),epswae(0,isig1,ifl(1,k),2))
     #		     *clr(ifl(4,k),1,isig3)*clr(ifl(3,k),1,isig3)
	     mvw4azu = dotcc(jvua(0,isig3,is3,3),epswze(0,isig1,ifl(1,k),2))
     #		     *clr(ifl(4,k),2,isig3)*clr(ifl(3,k),1,isig3)
	     mvw4zau = dotcc(jvuz(0,isig3,is3,3),epswae(0,isig1,ifl(1,k),2))
     #		     *clr(ifl(4,k),1,isig3)*clr(ifl(3,k),2,isig3)
	     mvw4zzu = dotcc(jvuz(0,isig3,is3,3),epswze(0,isig1,ifl(1,k),2))
     #		     *clr(ifl(4,k),2,isig3)*clr(ifl(3,k),2,isig3)
	   	   
	     mv4u = mvw4aau+mvw4azu+mvw4zau+mvw4zzu

	   
	   else
	     mv3u = 0d0
	     mv4u = 0d0
	    endif !k
	   
	   matv(k,isig1,isig3,3) = mv3e-mv3u+mv4e-mv4u
	   
	   enddo !k
	   
	   enddo !isig3
	enddo ! isig1
	
	endif !lbox  
	

c -----------------------------------------
c---------------
c
c	check full result for VB contributions:
	
	
	if (vbdebug) then
	   do k = 1,4
	   
	   do isig1 = -1,1,2
	   do isig3 = -1,1,2
	   	
	   if (k.eq.1.or.k.eq.3) then	   
c	   if (k.eq.1) then	   
	   m1 = 0d0	   
	   do i = 12,14
	      m1 = m1 + tampw(i,isig1,isig3,k)
	   enddo   
	   do i = 42,49
	      m1 = m1 + tampw(i,isig1,isig3,k)
	   enddo   
	   do i = 82,89
	      m1 = m1 + tampw(i,isig1,isig3,k)
	   enddo   
	   do i = 100,115
	      m1 = m1 + tampw(i,isig1,isig3,k)
	   enddo   
	   do i = 152,167
	      m1 = m1 + tampw(i,isig1,isig3,k)
	   enddo   
	   do i = 192,199
	      m1 = m1 + tampw(i,isig1,isig3,k)
	   enddo   
	
	   elseif (.false.) then	   
c	   elseif (k.eq.3) then	   
	   m1 = 0d0	   
	   do i = 23,30
	      m1 = m1 + tampw(i,isig1,isig3,k)
	   enddo   
	   do i = 52,54
	      m1 = m1 + tampw(i,isig1,isig3,k)
	   enddo   
	   do i = 82,89
	      m1 = m1 + tampw(i,isig1,isig3,k)
	   enddo   
	   do i = 100,115
	      m1 = m1 + tampw(i,isig1,isig3,k)
	   enddo   
	   do i = 152,159
	      m1 = m1 + tampw(i,isig1,isig3,k)
	   enddo   
	   do i = 171,186
	      m1 = m1 + tampw(i,isig1,isig3,k)
	   enddo   
	   	   
	   elseif (k.eq.2.or.k.eq.4) then	   
	   m4 = 0d0	   
	   do i = 23,30
	      m4 = m4 + tampw(i,isig1,isig3,k)
	   enddo   
	   do i = 52,54
	      m4 = m4 + tampw(i,isig1,isig3,k)
	   enddo   
	   do i = 82,89
	      m4 = m4 + tampw(i,isig1,isig3,k)
	   enddo   
	   do i = 100,115
	      m4 = m4 + tampw(i,isig1,isig3,k)
	   enddo   
	   do i = 152,167
	      m4 = m4 + tampw(i,isig1,isig3,k)
	   enddo   
	   do i = 192,199
	      m4 = m4 + tampw(i,isig1,isig3,k)
	   enddo   
	   endif
	   
	print*," VB for k = ",k, ", isig=",isig1,isig3
	if (m1.ne.0d0.or.m2.ne.0d0.or.m3.ne.0d0.or.m4.ne.0d0) then
	   if (k.eq.1) then
	     print*,"vb for k=",k, ", isig=",isig1,isig3
	     print*,"res:",m1,(mat(k,isig1,isig3,2)+mat(k,isig1,isig3,3))
	     	rat = m1/(mat(k,isig1,isig3,2)+mat(k,isig1,isig3,3))
	     print*,"rat:",rat
	   elseif (k.eq.3) then
	     print*,"vb for k=",k, ", isig=",isig1,isig3
	     print*,"res:",m1,(mat(k,isig3,isig1,2)+mat(k,isig3,isig1,3))
	     	rat = m1/(mat(k,isig3,isig1,2)+mat(k,isig3,isig1,3))
	     print*,"rat:",rat
	   elseif (k.eq.2) then
	     print*,"vb for k=",k, ", isig=",isig1,isig3
	     print*,"res:",m4,(mat(k,isig1,isig3,2)+mat(k,isig1,isig3,3))
	     	rat = m4/(mat(k,isig1,isig3,2)+mat(k,isig1,isig3,3))
	     print*,"rat:",rat
	   elseif (k.eq.4) then
	     print*,"vb for k=",k, ", isig=",isig1,isig3
	     print*,"res:",m4,(mat(k,isig3,isig1,2)+mat(k,isig3,isig1,3))
	     	rat = m4/(mat(k,isig3,isig1,2)+mat(k,isig3,isig1,3))
	     print*,"rat:",rat
	   endif
	   if (abs(abs(rat)-1.0).gt.1e-6) read(*,*)
        print*	   
	endif
	
		enddo !isig3
		enddo !isig1
	
	   enddo !k
	endif


c -----------------------------------------------------------------------
c
c next come the W-->WV currents attached to the quark lines. For the virtual 
c corrections the most effective structure is the contraction of two polarization 
c vectors with one fermion line. First build these effective polarization vectors
c from the current wwv(mu)
c
c NOTE: the wwv currents are NOT conserved. Hence it is necessary 
c       to consider q^mu * q^nu/m_W^2  terms in the W boson propagator

      if (lws(js1) .or. L.ne.Lold) then
         i = 1
	 isig = -1 ! only left-handed q-W coupling
            call ket2c(psi(1,isig,i),.true.,p(0,i),isig,pwz,wwv,
     1                 psiw(1,isig,is1,i),fq(0,i))
            call bra2c(psi(1,isig,i+1),.true.,p(0,i+1),isig,pwz,wwv,
     1                 psiw(1,isig,is1,i+1),fq(0,i+1))
	 
         call curr6(-1,psi(1,-1,i+1),p(0,i+1),
     1              psiw(1,-1,is1,i),fq(0,i), jw(0,-1,is1,i)   )
         call curr6(-1,psiw(1,-1,is1,i+1),fq(0,i+1),
     1              psi(1,-1,i),p(0,i), jw(0,-1,is1,i+1) )
 
#ifdef WITH_NLO              
	 if (lbox) then
            do i = 1,2
               call boxlinec(-1,psi(1,-1,1),psi(1,-1,2),p(0,1),p(0,2),
     1                       .true., wwv,3-i,   jw(0,-1,is1,i), 
     2                       jvw(0,-1,is1,i) )
            enddo
         endif
#endif
	 
      endif
      
           
      if (lws(js3) .or. L.ne.Lold) then
         i = 3
         isig = -1
            call ket2c(psi(1,isig,i),.true.,p(0,i),isig,pwz,wwv,
     1                 psiw(1,isig,is3,i),fq(0,i))
            call bra2c(psi(1,isig,i+1),.true.,p(0,i+1),isig,pwz,wwv,
     1                 psiw(1,isig,is3,i+1),fq(0,i+1))

         call curr6(-1,psi(1,-1,i+1),p(0,i+1),
     1              psiw(1,-1,is3,i),fq(0,i), jw(0,-1,is3,i)   )
         call curr6(-1,psiw(1,-1,is3,i+1),fq(0,i+1),
     1              psi(1,-1,i),p(0,i), jw(0,-1,is3,i+1) )
#ifdef WITH_NLO 
         if (lbox) then
            do i = 3,4
               call boxlinec(-1,psi(1,-1,3),psi(1,-1,4),p(0,3),p(0,4),
     1                       .true., wwv,5-i,   jw(0,-1,is3,i), 
     2                       jvw(0,-1,is3,i) )
            enddo
         endif
#endif
      endif
      
c extra terms for massive W propagator:      
      psiwg = -dotrc(pwz,wwv(0))*zm2i(3)*clr(3,3,-1)
      
      do isig1 = -1,1,2
      do isig3 = -1,1,2
      do kb = 1,3 ! boson id      
        jwg(1,isig1,isig3,kb) = psiwg*dotcc(jqq(0,isig1,1),jqq(0,isig3,2))
     $			              *prop43(kb)
        jwg(2,isig1,isig3,kb) = psiwg*dotcc(jqq(0,isig1,1),jqq(0,isig3,2))
     $			              *prop21(kb)
      enddo	
      enddo
      enddo

      propw(1) = clr(3,3,-1)**2*prop21(3)
      propw(2) = clr(3,3,-1)**2*prop43(3)
	       
      do isig3 = -1,1,2
	    isig1 = -1		! only left-handed quarks couple to W
	
c---------------------------

c box correction to upper line: polarization vectors are 
c    jqq(mu,isig3,2)=j43 with momentum    p43        and
c    wwv(mu)             with momentum    pwz

            mw(1) =      dotcc(jw(0,isig1,is1,2),jqq(0,isig3,2))
            mw(2) =      dotcc(jw(0,isig1,is1,1),jqq(0,isig3,2))
	    	    
            if (lbox) then
               mvw(1) = dotcc(jvw(0,isig1,is1,2),jqq(0,isig3,2))
               mvw(2) = dotcc(jvw(0,isig1,is1,1),jqq(0,isig3,2))
            endif
	    
            do k = 3,4
               propt1(isig1,isig3,k,2) = 
     1            clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3)*prop43(1) 
     2          + clr(ifl(1,k),2,isig1)*clr(ifl(3,k),2,isig3)*prop43(2)
               propt2(isig1,isig3,k,2) = 
     1            clr(ifl(2,k),1,isig1)*clr(ifl(3,k),1,isig3)*prop43(1) 
     2          + clr(ifl(2,k),2,isig1)*clr(ifl(3,k),2,isig3)*prop43(2)
     	    enddo
               
	    
	   if (isig3.eq.-1) then ! only left-handed quarks couple to W
	    do k = 1,2
	       kk = 3-k
	       mat(k,isig1,isig3,4) = propw(2) *
     1            	mw(kk)*clr(3,3,isig1) 
     !		      - jwg(1,isig1,isig3,3)*clr(3,3,isig1)*
     !		        clr(3,3,isig3)*(3-2*k)	
	       if (lbox) matv(k,isig1,isig3,4) = propw(2) *
     1            	mvw(kk)*clr(3,3,isig1)	    	    
	    enddo !k
	   endif !isig3
	   
	    
	    do k = 3,4 
	       mat(k,isig1,isig3,4) = propt2(isig1,isig3,k,2) *
     1            mw(2)*clr(ifl(1,k),3,isig1)   
     3    	+ propt1(isig1,isig3,k,2) *
     1            mw(1)*clr(ifl(2,k),3,isig1)        
     1		+ jwg(1,isig1,isig3,1)*(clr(ifl(1,k),1,isig1)
     1			-clr(ifl(2,k),1,isig1))*clr(ifl(3,k),1,isig3)
     1		+ jwg(1,isig1,isig3,2)*(clr(ifl(1,k),2,isig1)
     1		        -clr(ifl(2,k),2,isig1))*clr(ifl(3,k),2,isig3)
	       if (lbox) matv(k,isig1,isig3,4) = propt2(isig1,isig3,k,2) *
     1            ( mvw(2)*clr(ifl(1,k),3,isig1)   )
     3    	+ propt1(isig1,isig3,k,2) *
     1            ( mvw(1)*clr(ifl(2,k),3,isig1)   )

	    enddo
	
	
	enddo !isig3	    
	
	
c ------------------	    
 
         do isig1 = -1,1,2
	    isig3 = -1		! only left-handed quarks couple to W
	    
c box correction to lower line: polarization vectors are 
c    jqq(mu,isig1,1)=j21 with momentum    p21        and
c    aww/zww(mu)         with momentum    pzz

            mw(1) = dotcc(jw(0,isig3,is3,4),jqq(0,isig1,1))
            mw(2) = dotcc(jw(0,isig3,is3,3),jqq(0,isig1,1))
	    
            if (lbox) then
               mvw(1) = dotcc(jvw(0,isig3,is3,4),jqq(0,isig1,1))
               mvw(2) = dotcc(jvw(0,isig3,is3,3),jqq(0,isig1,1))
            endif


            do k = 1,2
               propt1(isig1,isig3,k,1) = 
     1            clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3)*prop21(1) 
     2          + clr(ifl(1,k),2,isig1)*clr(ifl(3,k),2,isig3)*prop21(2)
               propt2(isig1,isig3,k,1) = 
     1            clr(ifl(1,k),1,isig1)*clr(ifl(4,k),1,isig3)*prop21(1) 
     2          + clr(ifl(1,k),2,isig1)*clr(ifl(4,k),2,isig3)*prop21(2)
     	    enddo


	    do k = 1,2 
	       mat(k,isig1,isig3,5) = propt2(isig1,isig3,k,1) *
     1             mw(2)*clr(ifl(3,k),3,isig3)   
     3    	+ propt1(isig1,isig3,k,1) *
     1             mw(1)*clr(ifl(4,k),3,isig3)   
     1		+ jwg(2,isig1,isig3,1)*(clr(ifl(3,k),1,isig3)-clr(ifl(4,k),1,isig3))*
     1		            clr(ifl(1,k),1,isig1)
     1		+ jwg(2,isig1,isig3,2)*(clr(ifl(3,k),2,isig3)-clr(ifl(4,k),2,isig3))*
     1		            clr(ifl(1,k),2,isig1)
	       if (lbox) matv(k,isig1,isig3,5) = propt2(isig1,isig3,k,1) *
     1            ( mvw(2)*clr(ifl(3,k),3,isig3)   )
     3    	+ propt1(isig1,isig3,k,1) *
     1            ( mvw(1)*clr(ifl(4,k),3,isig3)   )

	    enddo
	

	   if (isig1.eq.-1) then ! only left-handed quarks couple to W
	    do k = 3,4
	       kk = 5-k
	       mat(k,isig1,isig3,5) = propw(1) *
     1            	mw(kk)*clr(3,3,isig1)
     !		      - jwg(2,isig1,isig3,3)*clr(3,3,isig1)*clr(3,3,isig3)*(7-2*k)
	       if (lbox) matv(k,isig1,isig3,5) = propw(1) *
     1            	mvw(kk)*clr(3,3,isig1)	    	    
	    enddo !k
	   endif !isig1

	
	enddo !isig1	    
	
c---------------
c
c	check full result for W->WV contributions:
	
	
	if (wwvdebug) then
	   do k = 1,4
	   
	   do isig1 = -1,1,2
	   do isig3 = -1,1,2
	   
c	   	
	   m1 = 0d0	   
	   if (k.eq.1.or.k.eq.3) then	   
	   do i = 9,11
	      m1 = m1 + tampw(i,isig1,isig3,k)
	   enddo   
	   do i = 17,19
	      m1 = m1 + tampw(i,isig1,isig3,k)
	   enddo   
	   do i = 36,41
	      m1 = m1 + tampw(i,isig1,isig3,k)
	   enddo   
	   do i = 54,59
	      m1 = m1 + tampw(i,isig1,isig3,k)
	   enddo   
	   do i = 76,81
	      m1 = m1 + tampw(i,isig1,isig3,k)
	   enddo   
	   do i = 94,99
	      m1 = m1 + tampw(i,isig1,isig3,k)
	   enddo   
	   
	   elseif (k.eq.2.or.k.eq.4) then
	   do i = 17,22
	      m1 = m1 + tampw(i,isig1,isig3,k)
	   enddo   
	   do i = 35,40
	      m1 = m1 + tampw(i,isig1,isig3,k)
	   enddo   
	   do i = 49,51
	      m1 = m1 + tampw(i,isig1,isig3,k)
	   enddo   
	   do i = 57,59
	      m1 = m1 + tampw(i,isig1,isig3,k)
	   enddo   
	   do i = 76,81
	      m1 = m1 + tampw(i,isig1,isig3,k)
	   enddo   
	   do i = 94,99
	      m1 = m1 + tampw(i,isig1,isig3,k)
	   enddo   
	   endif
	
	   
	   print*," W->WV for k = ",k
	   print*,"wwv for k=",k, ", isig=",isig1,isig3
	   if (m1.ne.0d0) then
	   if (k.le.2) then
	     print*,"res:",m1,(mat(k,isig1,isig3,4)+mat(k,isig1,isig3,5))
	     print*,"rat:",m1/(mat(k,isig1,isig3,4)+mat(k,isig1,isig3,5))
	   elseif (k.ge.3) then
	     print*,"res:",m1,(mat(k,isig3,isig1,4)+mat(k,isig3,isig1,5))
	     print*,"rat:",m1/(mat(k,isig3,isig1,4)+mat(k,isig3,isig1,5))
	   endif
           print*	   
	   endif
	
	   enddo !isig3
	   enddo !isig1
	
	   enddo !k
	endif

c ------------------------------------------------------------------------------
c
c  next do the box-box graphs with one boson (W/Z/A) emitted from the 
c  upper and the other from the lower line. 
c
c	   
c W from upper line, V from lower line:
	
	do isig3 = -1,1,2
	   isig1 = -1 !(only left-handed fermions couple to W)
c
	do idu = 1,2 ! fermion id for external W coupling
	do idl = 3,4 ! fermion id for external V coupling
	   mbwa(idu,idl,isig1,isig3) = dotcc(jew(0,isig1,idu),
     #				     jua(0,isig3,idl))*clr(3,3,isig1)
	   mbwz(idu,idl,isig1,isig3) = dotcc(jew(0,isig1,idu),
     #				     juz(0,isig3,idl))*clr(3,3,isig1)
	enddo !idl
	enddo !idu
	
	enddo !isig3
	   
c W from lower line, V from upper line:
	
	do isig1 = -1,1,2
	   isig3 = -1
c
	do idu = 1,2 ! fermion id for external V coupling
	do idl = 3,4 ! fermion id for external W coupling
	   mbaw(idu,idl,isig1,isig3) = dotcc(jew(0,isig3,idl),
     #				     jua(0,isig1,idu))*clr(3,3,isig3)
	   mbzw(idu,idl,isig1,isig3) = dotcc(jew(0,isig3,idl),
     #				     juz(0,isig1,idu))*clr(3,3,isig3)
	enddo !idl
	enddo !idu
	
	enddo !isig1
	   
c -------------------------------

	if (lbox) then
	   
c W from upper line, V from lower line, box corr. to upper (1) or lower line (2):
	
	do isig3 = -1,1,2
	   isig1 = -1 !(only left-handed fermions couple to W)
c
	do idu = 1,2 ! fermion id for external W coupling
	do idl = 3,4 ! fermion id for external V coupling
	   mbvwa(idu,idl,isig1,isig3,1) = dotcc(jvew(0,isig1,is1,idu),
     #		                         jua(0,isig3,idl))*clr(3,3,isig1)
	   mbvwz(idu,idl,isig1,isig3,1) = dotcc(jvew(0,isig1,is1,idu),
     #				         juz(0,isig3,idl))*clr(3,3,isig1)
	   mbvwa(idu,idl,isig1,isig3,2) = dotcc(jew(0,isig1,idu),
     #				         jvua(0,isig3,is3,idl))*clr(3,3,isig1)
	   mbvwz(idu,idl,isig1,isig3,2) = dotcc(jew(0,isig1,idu),
     #				         jvuz(0,isig3,is3,idl))*clr(3,3,isig1)
	enddo !idl
	enddo !idu
	
	enddo !isig3
		   
c W from lower line, V from upper line, box corr. to upper (1) or lower line (2):
	
	do isig1 = -1,1,2
	   isig3 = -1
c
	do idu = 1,2 ! fermion id for external V coupling
	do idl = 3,4 ! fermion id for external W coupling
	   mbvaw(idu,idl,isig1,isig3,1) = dotcc(jvew(0,isig3,is3,idl),
     #				     jua(0,isig1,idu))*clr(3,3,isig3)
	   mbvzw(idu,idl,isig1,isig3,1) = dotcc(jvew(0,isig3,is3,idl),
     #				     juz(0,isig1,idu))*clr(3,3,isig3)
	   mbvaw(idu,idl,isig1,isig3,2) = dotcc(jew(0,isig3,idl),
     #				     jvua(0,isig1,is1,idu))*clr(3,3,isig3)
	   mbvzw(idu,idl,isig1,isig3,2) = dotcc(jew(0,isig3,idl),
     #				     jvuz(0,isig1,is1,idu))*clr(3,3,isig3)
	enddo !idl
	enddo !idu
	
	enddo  !isig1
		
	endif !lbox

c  for the q^mu*q^nu/M_V^2 terms in the gauge boson propagators we need
	
	do i = 1,3 ! type of exchanged boson (Z,W,A)  (Z/A contributions cancel)
	
	 do isig3 = -1,1,2
	    isig1 = -1
	   
	    zwz(isig1,isig3,i) = jj21we(isig1)*jj43zu(isig3)*
     #				     zm2i(i)*clr(3,3,isig1)
	    zwa(isig1,isig3,i) = jj21we(isig1)*jj43au(isig3)*
     #				     zm2i(i)*clr(3,3,isig1)

	 enddo  !isig3
	 do isig1 = -1,1,2
	    isig3 = -1
	   
	    zzw(isig1,isig3,i) = jj21zu(isig1)*jj43we(isig3)*
     #				     zm2i(i)*clr(3,3,isig3)
	    zaw(isig1,isig3,i) = jj21au(isig1)*jj43we(isig3)*
     #				     zm2i(i)*clr(3,3,isig3)
	
	 enddo !isig1
	 
	enddo !i

c -----------------------------
c
c now construct amplitudes:	
c
	
	do isig1 = -1,1,2
	   isig3 = -1
	
	
	do k = 1,2
			
	if (isig1.eq.-1) then  
	   mbwv = (mbwa(k,3,isig1,isig3)*clr(ifl(3,k),1,isig3)+
     #	           mbwa(k,4,isig1,isig3)*clr(ifl(4,k),1,isig3)+	 
     #	           mbwz(k,3,isig1,isig3)*clr(ifl(3,k),2,isig3)+	
     #	           mbwz(k,4,isig1,isig3)*clr(ifl(4,k),2,isig3)
     #	           )*clr(3,3,isig1)*clr(3,3,isig3)*prop_ee(3)
     	   
	   gbwv = (zwz(isig1,isig3,3)*
     #	              (clr(ifl(3,k),2,isig3)-clr(ifl(4,k),2,isig3))+
     #	           zwa(isig1,isig3,3)*
     #	              (clr(ifl(3,k),1,isig3)-clr(ifl(4,k),1,isig3))
     #	           )*clr(3,3,isig1)*clr(3,3,isig3)*prop_ee(3)
     #		    *(3-2*k)	! +/- for k=1/2	
       
	else
	   mbwv = 0d0
	   gbwv = 0d0
	endif


	mbvw = 0d0
	gbvw = 0d0
	do kb = 1,2	
	   mbvw = mbvw +(
     #	         mbaw(1,3,isig1,isig3)*clr(ifl(1,k),1,isig1)*
     #	         	clr(ifl(2,k),kb,isig1)*clr(ifl(4,k),kb,isig3)+
     #	         mbaw(2,3,isig1,isig3)*clr(ifl(2,k),1,isig1)*
     #	         	clr(ifl(1,k),kb,isig1)*clr(ifl(4,k),kb,isig3)+
     #	         mbaw(1,4,isig1,isig3)*clr(ifl(1,k),1,isig1)*
     #	         	clr(ifl(2,k),kb,isig1)*clr(ifl(3,k),kb,isig3)+
     #	         mbaw(2,4,isig1,isig3)*clr(ifl(2,k),1,isig1)*
     #	         	clr(ifl(1,k),kb,isig1)*clr(ifl(3,k),kb,isig3)+
c
     #	         mbzw(1,3,isig1,isig3)*clr(ifl(1,k),2,isig1)*
     #	         	clr(ifl(2,k),kb,isig1)*clr(ifl(4,k),kb,isig3)+
     #	         mbzw(2,3,isig1,isig3)*clr(ifl(2,k),2,isig1)*
     #	         	clr(ifl(1,k),kb,isig1)*clr(ifl(4,k),kb,isig3)+
     #	         mbzw(1,4,isig1,isig3)*clr(ifl(1,k),2,isig1)*
     #	         	clr(ifl(2,k),kb,isig1)*clr(ifl(3,k),kb,isig3)+
     #	         mbzw(2,4,isig1,isig3)*clr(ifl(2,k),2,isig1)*
     #	         	clr(ifl(1,k),kb,isig1)*clr(ifl(3,k),kb,isig3))
     #		*prop_uu(kb)
     
	enddo !kb
		
	mat(k,isig1,isig3,6) = mbwv+mbvw+gbvw+gbwv
		
	enddo !k
	enddo !isig1


c --------------
	
	
	do isig3 = -1,1,2
	   isig1 = -1
	do k = 3,4
			
	if (isig3.eq.-1) then  
	   mbvw = (mbaw(1,k,isig1,isig3)*clr(ifl(1,k),1,isig1)+
     #	           mbaw(2,k,isig1,isig3)*clr(ifl(2,k),1,isig1)+  
     #	           mbzw(1,k,isig1,isig3)*clr(ifl(1,k),2,isig1)+ 
     #	           mbzw(2,k,isig1,isig3)*clr(ifl(2,k),2,isig1)
     #	           )*clr(3,3,isig3)*clr(3,3,isig1)*prop_uu(3)
     	   
	   gbvw = (zzw(isig1,isig3,3)*
     #	              (clr(ifl(1,k),2,isig1)-clr(ifl(2,k),2,isig1))+
     #	           zaw(isig1,isig3,3)*
     #	              (clr(ifl(1,k),1,isig1)-clr(ifl(2,k),1,isig1))
     #	           )*clr(3,3,isig3)*clr(3,3,isig1)*prop_uu(3)
     #		    *(7-2*k)	! +/- for k=3/4	
       
	else
	   mbvw = 0d0
	   gbvw = 0d0
	endif

	mbwv = 0d0
	gbwv = 0d0
	do kb = 1,2	
	   mbwv = mbwv +(
     #	         mbwa(1,3,isig1,isig3)*clr(ifl(3,k),1,isig3)*
     #	               clr(ifl(2,k),kb,isig1)*clr(ifl(4,k),kb,isig3)+
     #	         mbwa(2,3,isig1,isig3)*clr(ifl(4,k),1,isig3)*
     #	               clr(ifl(1,k),kb,isig1)*clr(ifl(4,k),kb,isig3)+
     #	         mbwa(1,4,isig1,isig3)*clr(ifl(3,k),1,isig3)*
     #	               clr(ifl(2,k),kb,isig1)*clr(ifl(3,k),kb,isig3)+
     #	         mbwa(2,4,isig1,isig3)*clr(ifl(4,k),1,isig3)*
     #	         	clr(ifl(1,k),kb,isig1)*clr(ifl(3,k),kb,isig3)+
c
     #	         mbwz(1,3,isig1,isig3)*clr(ifl(3,k),2,isig3)*
     #	               clr(ifl(2,k),kb,isig1)*clr(ifl(4,k),kb,isig3)+
     #	         mbwz(2,3,isig1,isig3)*clr(ifl(4,k),2,isig3)*
     #	               clr(ifl(1,k),kb,isig1)*clr(ifl(4,k),kb,isig3)+
     #	         mbwz(1,4,isig1,isig3)*clr(ifl(3,k),2,isig3)*
     #	               clr(ifl(2,k),kb,isig1)*clr(ifl(3,k),kb,isig3)+
     #	         mbwz(2,4,isig1,isig3)*clr(ifl(4,k),2,isig3)*
     #	         	clr(ifl(1,k),kb,isig1)*clr(ifl(3,k),kb,isig3))
     #		*prop_ee(kb)
     
           gbwv = 0d0

	enddo !kb


	mat(k,isig1,isig3,6) = mbwv+mbvw+gbvw+gbwv
		
	enddo !k
	enddo !isig
		
	
c ----------------------
	
	if (lbox) then
	
	do isig1 = -1,1,2
	   isig3 = -1

           do i = 1,2      ! 1 is for upper line, 2 for lower line QCD correction
            jj = 3 + 3*i ! stored in matv(...,6) and matv(...,9) respectively 
		
	     do k = 1,2
			
	if (isig1.eq.-1) then  
	  mbvwv = (mbvwa(k,3,isig1,isig3,i)*clr(ifl(3,k),1,isig3)+
     #	           mbvwa(k,4,isig1,isig3,i)*clr(ifl(4,k),1,isig3)+	 
     #	           mbvwz(k,3,isig1,isig3,i)*clr(ifl(3,k),2,isig3)+	
     #	           mbvwz(k,4,isig1,isig3,i)*clr(ifl(4,k),2,isig3)
     #	           )*clr(3,3,isig1)*clr(3,3,isig3)*prop_ee(3)
	else
	   mbvwv = 0d0
	endif

	mbvvw = 0d0
	do kb = 1,2	
	 mbvvw = mbvvw +(
     #	         mbvaw(1,3,isig1,isig3,i)*clr(ifl(1,k),1,isig1)*
     #	         	clr(ifl(2,k),kb,isig1)*clr(ifl(4,k),kb,isig3)+
     #	         mbvaw(2,3,isig1,isig3,i)*clr(ifl(2,k),1,isig1)*
     #	         	clr(ifl(1,k),kb,isig1)*clr(ifl(4,k),kb,isig3)+
     #	         mbvaw(1,4,isig1,isig3,i)*clr(ifl(1,k),1,isig1)*
     #	         	clr(ifl(2,k),kb,isig1)*clr(ifl(3,k),kb,isig3)+
     #	         mbvaw(2,4,isig1,isig3,i)*clr(ifl(2,k),1,isig1)*
     #	         	clr(ifl(1,k),kb,isig1)*clr(ifl(3,k),kb,isig3)+
c
     #	         mbvzw(1,3,isig1,isig3,i)*clr(ifl(1,k),2,isig1)*
     #	         	clr(ifl(2,k),kb,isig1)*clr(ifl(4,k),kb,isig3)+
     #	         mbvzw(2,3,isig1,isig3,i)*clr(ifl(2,k),2,isig1)*
     #	         	clr(ifl(1,k),kb,isig1)*clr(ifl(4,k),kb,isig3)+
     #	         mbvzw(1,4,isig1,isig3,i)*clr(ifl(1,k),2,isig1)*
     #	         	clr(ifl(2,k),kb,isig1)*clr(ifl(3,k),kb,isig3)+
     #	         mbvzw(2,4,isig1,isig3,i)*clr(ifl(2,k),2,isig1)*
     #	         	clr(ifl(1,k),kb,isig1)*clr(ifl(3,k),kb,isig3))
     #		*prop_uu(kb)

	enddo !kb

	matv(k,isig1,isig3,jj) = mbvwv+mbvvw
	
	     enddo !k
	
	   enddo !i (jj)
		
	enddo !isig1

c --------------	
	
	do isig3 = -1,1,2
	   isig1 = -1
    
           do i = 1,2      ! 1 is for upper line, 2 for lower line QCD correction
            jj = 3 + 3*i ! stored in matv(...,6) and matv(...,9) respectively 
	
	    do k = 3,4
			
	if (isig3.eq.-1) then  
	  mbvvw = (mbvaw(1,k,isig1,isig3,i)*clr(ifl(1,k),1,isig1)+
     #	           mbvaw(2,k,isig1,isig3,i)*clr(ifl(2,k),1,isig1)+  
     #	           mbvzw(1,k,isig1,isig3,i)*clr(ifl(1,k),2,isig1)+ 
     #	           mbvzw(2,k,isig1,isig3,i)*clr(ifl(2,k),2,isig1)
     #	           )*clr(3,3,isig3)*clr(3,3,isig1)*prop_uu(3)
     	   
     
	else
	   mbvvw = 0d0
	endif

	mbvwv = 0d0
	do kb = 1,2	
	 mbvwv = mbvwv +(
     #	         mbvwa(1,3,isig1,isig3,i)*clr(ifl(3,k),1,isig3)*
     #	               clr(ifl(2,k),kb,isig1)*clr(ifl(4,k),kb,isig3)+
     #	         mbvwa(2,3,isig1,isig3,i)*clr(ifl(4,k),1,isig3)*
     #	               clr(ifl(1,k),kb,isig1)*clr(ifl(4,k),kb,isig3)+
     #	         mbvwa(1,4,isig1,isig3,i)*clr(ifl(3,k),1,isig3)*
     #	               clr(ifl(2,k),kb,isig1)*clr(ifl(3,k),kb,isig3)+
     #	         mbvwa(2,4,isig1,isig3,i)*clr(ifl(4,k),1,isig3)*
     #	         	clr(ifl(1,k),kb,isig1)*clr(ifl(3,k),kb,isig3)+
c
     #	         mbvwz(1,3,isig1,isig3,i)*clr(ifl(3,k),2,isig3)*
     #	               clr(ifl(2,k),kb,isig1)*clr(ifl(4,k),kb,isig3)+
     #	         mbvwz(2,3,isig1,isig3,i)*clr(ifl(4,k),2,isig3)*
     #	               clr(ifl(1,k),kb,isig1)*clr(ifl(4,k),kb,isig3)+
     #	         mbvwz(1,4,isig1,isig3,i)*clr(ifl(3,k),2,isig3)*
     #	               clr(ifl(2,k),kb,isig1)*clr(ifl(3,k),kb,isig3)+
     #	         mbvwz(2,4,isig1,isig3,i)*clr(ifl(4,k),2,isig3)*
     #	         	clr(ifl(1,k),kb,isig1)*clr(ifl(3,k),kb,isig3))
     #		*prop_ee(kb)
     
      
	enddo !kb
	
	matv(k,isig1,isig3,jj) = mbvwv+mbvvw
	

            enddo !k
	
	  enddo !i (jj)
		
	enddo !isig3


	endif !lbox
	
c -------------------------------

	if (bbdebug) then
	
	do k = 1,4
	do isig1 = -1,1,2
	do isig3 = -1,1,2
	m1 = 0d0
	
	if (k.eq.1.or.k.eq.3) then
	do i = 3,4
	   m1 = m1+tampw(i,isig1,isig3,k)
	enddo
	do i = 7,8
	   m1 = m1+tampw(i,isig1,isig3,k)
	enddo
	do i = 20,23
	   m1 = m1+tampw(i,isig1,isig3,k)
	enddo
	do i = 28,31
	   m1 = m1+tampw(i,isig1,isig3,k)
	enddo
	do i = 60,63
	   m1 = m1+tampw(i,isig1,isig3,k)
	enddo
	do i = 68,71
	   m1 = m1+tampw(i,isig1,isig3,k)
	enddo
	endif
	
	if (k.eq.2.or.k.eq.4) then
	do i = 1,4
	   m1 = m1+tampw(i,isig1,isig3,k)
	enddo
	do i = 9,12
	   m1 = m1+tampw(i,isig1,isig3,k)
	enddo
	do i = 43,44
	   m1 = m1+tampw(i,isig1,isig3,k)
	enddo
	do i = 47,48
	   m1 = m1+tampw(i,isig1,isig3,k)
	enddo
	do i = 60,63
	   m1 = m1+tampw(i,isig1,isig3,k)
	enddo
	do i = 68,71
	   m1 = m1+tampw(i,isig1,isig3,k)
	enddo
	endif

	
	print*,"comparison for BB with k=",k," isig1,3=",isig1,isig3
	if (m1.ne.0d0) then
	if (k.eq.1.or.k.eq.2) then
	  print*,"m1,mat = ",m1,mat(k,isig1,isig3,6)
	  print*,"m1/mat = ",m1/mat(k,isig1,isig3,6)
	elseif (k.eq.3.or.k.eq.4) then
	  print*,"m1,mat = ",m1,mat(k,isig3,isig1,6)
	  print*,"m1/mat = ",m1/mat(k,isig3,isig1,6)
	endif
	print*
	endif
	
	
	enddo !isig1
	enddo !isig3
	
	enddo !k
	
	endif
        if (nlo.ge.1) then
        endif
c
c ------------------------------------------------------------------------------
c
c and now, finally, the pentagon contributions, i.e. three bosons emitted 
c from the same quark line
c
c prerequisites for virtual corrections:
            if (ldebug) then	
	       if (lpt.and.icount1.eq.0) then
c		  print*, "criterion chosen: mv5/mvc.gt.0.1d0"
	       endif !icount1
            endif

#ifdef WITH_NLO 	
      if (lpent .and. lws(js1)) then     ! need new pentagon graphs for 12 line
         musq = -p21(4)
         if (NLO.eq.1) then
            do k = 1,3
               call J_virtual_tri_box_pent45(psi,p,1,2,sign, musq,
     1              we, zu, zero, zero, k, -5, j5wz(0,k,is1,1))
               call J_virtual_tri_box_pent45(psi,p,1,2,sign, musq,
     1              zu, we, zero, zero, k, -5, j5zw(0,k,is1,1))               
	       call J_virtual_tri_box_pent45(psi,p,1,2,sign, musq,
     1              we, au, zero, zero, k, -5, j5wa(0,k,is1,1))
               call J_virtual_tri_box_pent45(psi,p,1,2,sign, musq,
     1              au, we, zero, zero, k, -5, j5aw(0,k,is1,1))               
            enddo
            icount1= icount1+1
	 elseif (nlo.eq.-5) then
            do k = 1,3
               call J_virtual_tri_box_pent45(psi,p,1,2,sign, musq,
     1              wet, zut, xxew, xxuz, k, -5, j5wz(0,k,is1,1))                              
               call J_virtual_tri_box_pent45(psi,p,1,2,sign, musq,
     1              zut, wet, xxuz, xxew, k, -5, j5zw(0,k,is1,1))               
	       call J_virtual_tri_box_pent45(psi,p,1,2,sign, musq,
     1              wet, aut, xxew, xxua, k, -5, j5wa(0,k,is1,1))
               call J_virtual_tri_box_pent45(psi,p,1,2,sign, musq,
     1              aut, wet, xxua, xxew, k, -5, j5aw(0,k,is1,1))             
            enddo 
            icount1= icount1+1      
	 else
            do k = 1,3
               call J_virtual_tri_box_pent45(psi,p,1,2,sign, musq,
     1              wet, zut, xxew, xxuz, k, 5, j5wz(0,k,is1,1))
               call J_virtual_tri_box_pent45(psi,p,1,2,sign, musq,
     1              zut, wet, xxuz, xxew, k, 5, j5zw(0,k,is1,1))               
	       call J_virtual_tri_box_pent45(psi,p,1,2,sign, musq,
     1              wet, aut, xxew, xxua, k, 5, j5wa(0,k,is1,1))
               call J_virtual_tri_box_pent45(psi,p,1,2,sign, musq,
     1              aut, wet, xxua, xxew, k, 5, j5aw(0,k,is1,1))               
            enddo
         endif


         lgc(js1) = .false.
         lws(js1) = .false.
         if (lpt) then
c	 
	    icountmax = 2**(N_points(0)-9)
	    if(icountmax.eq.1) icountmax = 10
c
c gauge check (compare contraction of pentagonline with momentum 
c	of attached current to corresponding boxline result 
c	... should give ratio ~1 ):
c            if (sign(1).eq.2 .and. sign(3).eq.1) then
            if (.true.) then
               call gauge_check_5(j5wz(0,1,is1,1),p43,1,psi,p,1,2,wet,zut,0d0)
	         if(bad_gauge_sin) bad_gauge = .true.

	       call gauge_check_5(j5zw(0,1,is1,1),p43,1,psi,p,1,2,zut,wet,0d0)
	         if(bad_gauge_sin) bad_gauge = .true.
	      
	       call gauge_check_5(j5wa(0,1,is1,1),p43,1,psi,p,1,2,wet,aut,0d0)
	         if(bad_gauge_sin) bad_gauge = .true.
	       
	       call gauge_check_5(j5aw(0,1,is1,1),p43,1,psi,p,1,2,aut,wet,0d0)
	         if(bad_gauge_sin) bad_gauge = .true.
  	    endif
	    if(bad_gauge) then
               icb1 = icb1+1
               lgc(js1) = .true.
            endif

c -------------------------	   
c  check numerical stability by gauge invariance (sum of pentagons for various k 
c	should give zero after contraction with momentum: 
c
c don"t employ this gauge test but the stonger one invoked above:

	if (.false.) then
            
	    mvwz = 0
            mvzw = 0
            xvwz = 0 
            xvzw = 0
            
	    mvwa = 0
            mvaw = 0
            xvwa = 0 
            xvaw = 0
            
    	    
	    do k = 1,3	! sum up contractions of pentagons with momentum 
               mvwz = mvwz+dotrc(p43,j5wz(0,k,is1,1)) ! sum should vanish
               xvwz = xvwz+abs(dotrc(p43,j5wz(0,k,is1,1))) !normalize to "relative size" 
               mvzw = mvzw+dotrc(p43,j5zw(0,k,is1,1))
               xvzw = xvzw+abs(dotrc(p43,j5zw(0,k,is1,1)))

	       mvwa = mvwa+dotrc(p43,j5wa(0,k,is1,1))
	       xvwa = xvwa+abs(dotrc(p43,j5wa(0,k,is1,1)))
	       mvaw = mvaw+dotrc(p43,j5aw(0,k,is1,1))
	       xvaw = xvaw+abs(dotrc(p43,j5aw(0,k,is1,1)))
	
             enddo !k

	     
	    if ((abs(mvzw/xvzw).gt.3d-2 .or.abs(mvzw/xvzw).gt.3d-2).or.
     #   	(abs(mvwz/xvwz).gt.3d-2 .or.abs(mvwz/xvwz).gt.3d-2).or.
     #   	(abs(mvaw/xvaw).gt.3d-2 .or.abs(mvaw/xvaw).gt.3d-2).or.
     #   	(abs(mvwa/xvwa).gt.3d-2 .or.abs(mvwa/xvwa).gt.3d-2)) 
     #	    then
               icb1 = icb1+1
c               print*," zw upper for isig = ",isig,abs(mvzw/xvzw)
c               print*," wz upper: ",abs(mvwz/xvwz)
c               print*," aw upper: ",abs(mvaw/xvaw)
c               print*," wa upper: ",abs(mvwa/xvwa)
               lgc(js1) = .true.
            endif !abs...
	    
	    
	   endif ! omit this gauge test 
c -------------------------	 
           if(ldebug) then
	      if (mod(icount1,icountmax).eq.0 ) then
                 print*," bad gauge check fraction (upper line)",
     #				       real(icb1)/real(2*icount1), js1
              endif
           endif
         endif   ! lpt
      endif
c pentagon corr. to lower line:
      if (lpent .and. lws(js3)) then     ! need new pentagon graphs for 34 line
         musq = -p43(4)
         if (NLO.eq.1) then
            do k = 1,3
               call J_virtual_tri_box_pent45(psi,p,3,4,sign, musq,
     1              we, zu, zero, zero, k, -5, j5wz(0,k,is3,2))
               call J_virtual_tri_box_pent45(psi,p,3,4,sign, musq,
     1              zu, we, zero, zero, k, -5, j5zw(0,k,is3,2))               
	       call J_virtual_tri_box_pent45(psi,p,3,4,sign, musq,
     1              we, au, zero, zero, k, -5, j5wa(0,k,is3,2))
               call J_virtual_tri_box_pent45(psi,p,3,4,sign, musq,
     1              au, we, zero, zero, k, -5, j5aw(0,k,is3,2))               
            enddo
            icount2= icount2+1
         
	 elseif (nlo.eq.-5) then
            do k = 1,3
               call J_virtual_tri_box_pent45(psi,p,3,4,sign, musq,
     1              wet, zut, xxew, xxuz, k, -5, j5wz(0,k,is3,2))
               call J_virtual_tri_box_pent45(psi,p,3,4,sign, musq,
     1              zut, wet, xxuz, xxew, k, -5, j5zw(0,k,is3,2))               
	       call J_virtual_tri_box_pent45(psi,p,3,4,sign, musq,
     1              wet, aut, xxew, xxua, k, -5, j5wa(0,k,is3,2))
               call J_virtual_tri_box_pent45(psi,p,3,4,sign, musq,
     1              aut, wet, xxua, xxew, k, -5, j5aw(0,k,is3,2))               
            enddo
            icount2= icount2+1
        
	 else
            do k = 1,3
               call J_virtual_tri_box_pent45(psi,p,3,4,sign, musq,
     1              wet, zut, xxew, xxuz, k, 5, j5wz(0,k,is3,2))
               call J_virtual_tri_box_pent45(psi,p,3,4,sign, musq,
     1              zut, wet, xxuz, xxew, k, 5, j5zw(0,k,is3,2))               
	       call J_virtual_tri_box_pent45(psi,p,3,4,sign, musq,
     1              wet, aut, xxew, xxua, k, 5, j5wa(0,k,is3,2))
               call J_virtual_tri_box_pent45(psi,p,3,4,sign, musq,
     1              aut, wet, xxua, xxew, k, 5, j5aw(0,k,is3,2))               
            enddo
         endif
 
         lgc(js3) = .false.
         lws(js3) = .false.
         if (lpt) then
c perform gauge check:
            if (.true.) then
               call gauge_check_5(j5wz(0,1,is3,2),p21,1,psi,p,3,4,wet,zut,0d0)
	         if(bad_gauge_sin) bad_gauge = .true.
	       
	       call gauge_check_5(j5zw(0,1,is3,2),p21,1,psi,p,3,4,zut,wet,0d0)
	         if(bad_gauge_sin) bad_gauge = .true.
	      
	       call gauge_check_5(j5wa(0,1,is3,2),p21,1,psi,p,3,4,wet,aut,0d0)
	         if(bad_gauge_sin) bad_gauge = .true.
	       
	       call gauge_check_5(j5aw(0,1,is3,2),p21,1,psi,p,3,4,aut,wet,0d0)
	         if(bad_gauge_sin) bad_gauge = .true.
	       
  	    endif
	    if(bad_gauge) then
               icb2 = icb2+1
               lgc(js3) = .true.
            endif

c don"t employ this gauge test but the stonger one invoked above:
	if (.false.) then
            
	    mvwz = 0
            mvzw = 0
            xvwz = 0 
            xvzw = 0
            
	    mvwa = 0
            mvaw = 0
            xvwa = 0 
            xvaw = 0
            
    	    
	    do k = 1,3	! sum up contractions of pentagons with momentum 
               mvwz = mvwz+dotrc(p21,j5wz(0,k,is3,2)) ! sum should vanish
               xvwz = xvwz+abs(dotrc(p21,j5wz(0,k,is3,2))) !normalize to "relative size" 
               mvzw = mvzw+dotrc(p21,j5zw(0,k,is3,2))
               xvzw = xvzw+abs(dotrc(p21,j5zw(0,k,is3,2)))

	       mvwa = mvwa+dotrc(p21,j5wa(0,k,is3,2))
	       xvwa = xvwa+abs(dotrc(p21,j5wa(0,k,is3,2)))
	       mvaw = mvaw+dotrc(p21,j5aw(0,k,is3,2))
	       xvaw = xvaw+abs(dotrc(p21,j5aw(0,k,is3,2)))
	
             enddo !k
	     
	    if ((abs(mvzw/xvzw).gt.3d-2 .or.abs(mvzw/xvzw).gt.3d-2).or.
     #   	(abs(mvwz/xvwz).gt.3d-2 .or.abs(mvwz/xvwz).gt.3d-2).or.
     #   	(abs(mvaw/xvaw).gt.3d-2 .or.abs(mvaw/xvaw).gt.3d-2).or.
     #   	(abs(mvwa/xvwa).gt.3d-2 .or.abs(mvwa/xvwa).gt.3d-2)) 
     #	    then
               icb2 = icb2+1
c               print*," zw upper for isig = ",isig,abs(mvzw/xvzw)
c               print*," wz upper: ",abs(mvwz/xvwz)
c               print*," aw upper: ",abs(mvaw/xvaw)
c               print*," wa upper: ",abs(mvwa/xvwa)
               lgc(js3) = .true.
            endif !abs...
	    	    
	   endif ! omit this gauge test 

c -------------------------	   
           if(ldebug) then
              if (mod(icount2,icountmax).eq.0 ) then
                 print*," bad gauge check fraction (lower line)",
     #				       real(icb2)/real(2*icount2), js3
              endif
           endif
         endif   ! lpt
      endif
#endif
c ------------------------------
c ------------------------------

c compute matrix elements for pentagon on upper line:
c
c need m5(pos.j43=1:3,vext1=1:3,vext2=1:3), id of internal boson arbitrary 

      do isig3 = -1,1,2 
	 isig1 = -1

         call ket2c(psi(1,isig1,1),.true.,p(0,1),isig1,
     1        p43,jqq(0,isig3,2),bkjqq(1,isig1,isig3,1),dummy)
         call bra2c(psi(1,isig1,2),.true.,p(0,2),isig1,
     1        p43,jqq(0,isig3,2),bkjqq(1,isig1,isig3,2),dummy)
	

c  eps1=j43 with v1 arbitrary,eps2=v2,eps3=v3
         m5(1,3,1) = s1c(psiau(1,isig1,2),we,.true.,isig1,bkjqq(1,isig1,isig3,1))
         m5(1,3,2) = s1c(psizu(1,isig1,2),we,.true.,isig1,bkjqq(1,isig1,isig3,1))

         m5(1,1,3) = s1c(psiwe(1,isig1,2),au,.true.,isig1,bkjqq(1,isig1,isig3,1))
         m5(1,2,3) = s1c(psiwe(1,isig1,2),zu,.true.,isig1,bkjqq(1,isig1,isig3,1))

c  eps1=v1,eps2=j43 with v2 arbitrary,eps3=v3
         m5(2,3,1) = s1c(psiau(1,isig1,2),jqq(0,isig3,2),.true.,isig1,psiwe(1,isig1,1))
         m5(2,3,2) = s1c(psizu(1,isig1,2),jqq(0,isig3,2),.true.,isig1,psiwe(1,isig1,1))

         m5(2,1,3) = s1c(psiwe(1,isig1,2),jqq(0,isig3,2),.true.,isig1,psiau(1,isig1,1))
         m5(2,2,3) = s1c(psiwe(1,isig1,2),jqq(0,isig3,2),.true.,isig1,psizu(1,isig1,1))

c  eps1=v1,eps2=v2,eps3=j43 with v3 arbitrary
         m5(3,3,1) = s1c(bkjqq(1,isig1,isig3,2),au,.true.,isig1,psiwe(1,isig1,1))
         m5(3,3,2) = s1c(bkjqq(1,isig1,isig3,2),zu,.true.,isig1,psiwe(1,isig1,1))
         
	 m5(3,1,3) = s1c(bkjqq(1,isig1,isig3,2),we,.true.,isig1,psiau(1,isig1,1))
         m5(3,2,3) = s1c(bkjqq(1,isig1,isig3,2),we,.true.,isig1,psizu(1,isig1,1))


	if (lpent) then	
	   do j = 1,3	   
	      mv5(j,3,2) = dotcc(j5wz(0,j,is1,1),jqq(0,isig3,2))
	      mv5(j,2,3) = dotcc(j5zw(0,j,is1,1),jqq(0,isig3,2))
	      mv5(j,3,1) = dotcc(j5wa(0,j,is1,1),jqq(0,isig3,2))
	      mv5(j,1,3) = dotcc(j5aw(0,j,is1,1),jqq(0,isig3,2))		
	   enddo	   
	endif !lpent


	do kb = 1,2
	
	if (isig3.eq.-1) then 
	
	k = 1
	
	mat(k,isig1,isig3,7) = mat(k,isig1,isig3,7) + clr(3,3,-1)**3*prop43(3)*(
     #	       m5(2,3,kb)*clr(ifl(2,k),kb,isig1)
     #	      +m5(3,3,kb)*clr(7-ifl(1,k),kb,isig1)
     #	      +m5(3,kb,3)*clr(ifl(1,k),kb,isig1))
	
	if (lpent) matv(k,isig1,isig3,7) = matv(k,isig1,isig3,7) + 
     #		clr(3,3,-1)**3*prop43(3)*(
     #	       mv5(2,3,kb)*clr(ifl(2,k),kb,isig1)
     #	      +mv5(3,3,kb)*clr(7-ifl(1,k),kb,isig1)
     #	      +mv5(3,kb,3)*clr(ifl(1,k),kb,isig1))

	k = 2
	
	mat(k,isig1,isig3,7) = mat(k,isig1,isig3,7) + clr(3,3,-1)**3*prop43(3)*(
     #	       m5(1,3,kb)*clr(ifl(2,k),kb,isig1)
     #	      +m5(1,kb,3)*clr(7-ifl(2,k),kb,isig1)
     #	      +m5(2,kb,3)*clr(ifl(1,k),kb,isig1))
     
	if (lpent) matv(k,isig1,isig3,7) = matv(k,isig1,isig3,7) 
     #		+ clr(3,3,-1)**3*prop43(3)*(
     #	       mv5(1,3,kb)*clr(ifl(2,k),kb,isig1)
     #	      +mv5(1,kb,3)*clr(7-ifl(2,k),kb,isig1)
     #	      +mv5(2,kb,3)*clr(ifl(1,k),kb,isig1))
     
     	endif !isig3
     
        do k= 3,4
	
	do kp = 1,2
	
	mat(k,isig1,isig3,7) = mat(k,isig1,isig3,7) + 
     #			       clr(3,3,-1)*prop43(kp)*clr(ifl(3,k),kp,isig3)*(
     #	       m5(1,3,kb)*clr(ifl(2,k),kb,isig1)*clr(ifl(1,k),kp,isig1)
     #	      +m5(1,kb,3)*clr(7-ifl(2,k),kb,isig1)*clr(ifl(1,k),kp,isig1)
     #	      +m5(2,3,kb)*clr(ifl(2,k),kb,isig1)*clr(ifl(2,k),kp,isig1)
     #	      +m5(2,kb,3)*clr(ifl(1,k),kb,isig1)*clr(7-ifl(2,k),kp,isig1)
     #	      +m5(3,3,kb)*clr(7-ifl(1,k),kb,isig1)*clr(ifl(2,k),kp,isig1)
     #	      +m5(3,kb,3)*clr(ifl(1,k),kb,isig1)*clr(ifl(2,k),kp,isig1))

	
	if (lpent) matv(k,isig1,isig3,7) = matv(k,isig1,isig3,7) + 
     #		 clr(3,3,-1)*prop43(kp)*clr(ifl(3,k),kp,isig3)*(
     #	       mv5(1,3,kb)*clr(ifl(2,k),kb,isig1)*clr(ifl(1,k),kp,isig1)
     #	      +mv5(1,kb,3)*clr(7-ifl(2,k),kb,isig1)*clr(ifl(1,k),kp,isig1)
     #	      +mv5(2,3,kb)*clr(ifl(2,k),kb,isig1)*clr(ifl(2,k),kp,isig1)
     #	      +mv5(2,kb,3)*clr(ifl(1,k),kb,isig1)*clr(7-ifl(2,k),kp,isig1)
     #	      +mv5(3,3,kb)*clr(7-ifl(1,k),kb,isig1)*clr(ifl(2,k),kp,isig1)
     #	      +mv5(3,kb,3)*clr(ifl(1,k),kb,isig1)*clr(ifl(2,k),kp,isig1))


	enddo !kp

	enddo !k
	enddo !kb


      enddo !isig3
      
c ------------------------------

c same for lower line:
c
c need m5(pos.j21=1:3,vext1=1:3,vext2=1:3), id of internal boson arbitrary 


      do isig1 = -1,1,2 
	 isig3 = -1

         call ket2c(psi(1,isig3,3),.true.,p(0,3),isig3,
     1        p21,jqq(0,isig1,1),bkjqq(1,isig1,isig3,3),dummy)
         call bra2c(psi(1,isig3,4),.true.,p(0,4),isig3,
     1        p21,jqq(0,isig1,1),bkjqq(1,isig1,isig3,4),dummy)
	

c  eps1=j21 with v1 arbitrary,eps2=v2,eps3=v3
         m5(1,3,1) = s1c(psiau(1,isig3,4),we,.true.,isig3,bkjqq(1,isig1,isig3,3))
         m5(1,3,2) = s1c(psizu(1,isig3,4),we,.true.,isig3,bkjqq(1,isig1,isig3,3))

         m5(1,1,3) = s1c(psiwe(1,isig3,4),au,.true.,isig3,bkjqq(1,isig1,isig3,3))
         m5(1,2,3) = s1c(psiwe(1,isig3,4),zu,.true.,isig3,bkjqq(1,isig1,isig3,3))

c  eps1=v1,eps2=j21 with v2 arbitrary,eps3=v3
         m5(2,3,1) = s1c(psiau(1,isig3,4),jqq(0,isig1,1),.true.,isig3,psiwe(1,isig3,3))
         m5(2,3,2) = s1c(psizu(1,isig3,4),jqq(0,isig1,1),.true.,isig3,psiwe(1,isig3,3))

         m5(2,1,3) = s1c(psiwe(1,isig3,4),jqq(0,isig1,1),.true.,isig3,psiau(1,isig3,3))
         m5(2,2,3) = s1c(psiwe(1,isig3,4),jqq(0,isig1,1),.true.,isig3,psizu(1,isig3,3))

c  eps1=v1,eps2=v2,eps3=j21 with v3 arbitrary
         m5(3,3,1) = s1c(bkjqq(1,isig1,isig3,4),au,.true.,isig3,psiwe(1,isig3,3))
         m5(3,3,2) = s1c(bkjqq(1,isig1,isig3,4),zu,.true.,isig3,psiwe(1,isig3,3))
         
	 m5(3,1,3) = s1c(bkjqq(1,isig1,isig3,4),we,.true.,isig3,psiau(1,isig3,3))
         m5(3,2,3) = s1c(bkjqq(1,isig1,isig3,4),we,.true.,isig3,psizu(1,isig3,3))


	if (lpent) then	
	   do j = 1,3	   
	      mv5(j,3,2) = dotcc(j5wz(0,j,is3,2),jqq(0,isig1,1))
	      mv5(j,2,3) = dotcc(j5zw(0,j,is3,2),jqq(0,isig1,1))
	      mv5(j,3,1) = dotcc(j5wa(0,j,is3,2),jqq(0,isig1,1))
	      mv5(j,1,3) = dotcc(j5aw(0,j,is3,2),jqq(0,isig1,1)) 	    
	   enddo	   
	endif !lpent


	do kb = 1,2
	
	if (isig1.eq.-1) then 
	
	k = 3
	
	mat(k,isig1,isig3,8) = mat(k,isig1,isig3,8) + clr(3,3,-1)**3*prop21(3)*(
     #	       m5(2,3,kb)*clr(ifl(4,k),kb,isig3)
     #	      +m5(3,3,kb)*clr(7-ifl(3,k),kb,isig3)
     #	      +m5(3,kb,3)*clr(ifl(3,k),kb,isig3))
	
	if (lpent) matv(k,isig1,isig3,8) = matv(k,isig1,isig3,8) 
     #		+ clr(3,3,-1)**3*prop21(3)*(
     #	       mv5(2,3,kb)*clr(ifl(4,k),kb,isig3)
     #	      +mv5(3,3,kb)*clr(7-ifl(3,k),kb,isig3)
     #	      +mv5(3,kb,3)*clr(ifl(3,k),kb,isig3))

	k = 4
	
	mat(k,isig1,isig3,8) = mat(k,isig1,isig3,8) + clr(3,3,-1)**3*prop21(3)*(
     #	       m5(1,3,kb)*clr(ifl(4,k),kb,isig3)
     #	      +m5(1,kb,3)*clr(7-ifl(4,k),kb,isig3)
     #	      +m5(2,kb,3)*clr(ifl(3,k),kb,isig3))
     
	if (lpent)  matv(k,isig1,isig3,8) = matv(k,isig1,isig3,8) 
     #		+ clr(3,3,-1)**3*prop21(3)*(
     #	       mv5(1,3,kb)*clr(ifl(4,k),kb,isig3)
     #	      +mv5(1,kb,3)*clr(7-ifl(4,k),kb,isig3)
     #	      +mv5(2,kb,3)*clr(ifl(3,k),kb,isig3)) 
     
     	endif !isig1
     
        do k= 1,2
	
	do kp = 1,2
	
	mat(k,isig1,isig3,8) = mat(k,isig1,isig3,8) + 
     #			       clr(3,3,-1)*prop21(kp)*clr(ifl(1,k),kp,isig1)*(
     #	       m5(1,3,kb)*clr(ifl(4,k),kb,isig3)*clr(ifl(3,k),kp,isig3)
     #	      +m5(1,kb,3)*clr(7-ifl(4,k),kb,isig3)*clr(ifl(3,k),kp,isig3)
     #	      +m5(2,3,kb)*clr(ifl(4,k),kb,isig3)*clr(ifl(4,k),kp,isig3)
     #	      +m5(2,kb,3)*clr(ifl(3,k),kb,isig3)*clr(7-ifl(4,k),kp,isig3)
     #	      +m5(3,3,kb)*clr(7-ifl(3,k),kb,isig3)*clr(ifl(4,k),kp,isig3)
     #	      +m5(3,kb,3)*clr(ifl(3,k),kb,isig3)*clr(ifl(4,k),kp,isig3))

	
	if (lpent) matv(k,isig1,isig3,8) = matv(k,isig1,isig3,8) + 
     #			       clr(3,3,-1)*prop21(kp)*clr(ifl(1,k),kp,isig1)*(
     #	       mv5(1,3,kb)*clr(ifl(4,k),kb,isig3)*clr(ifl(3,k),kp,isig3)
     #	      +mv5(1,kb,3)*clr(7-ifl(4,k),kb,isig3)*clr(ifl(3,k),kp,isig3)
     #	      +mv5(2,3,kb)*clr(ifl(4,k),kb,isig3)*clr(ifl(4,k),kp,isig3)
     #	      +mv5(2,kb,3)*clr(ifl(3,k),kb,isig3)*clr(7-ifl(4,k),kp,isig3)
     #	      +mv5(3,3,kb)*clr(7-ifl(3,k),kb,isig3)*clr(ifl(4,k),kp,isig3)
     #	      +mv5(3,kb,3)*clr(ifl(3,k),kb,isig3)*clr(ifl(4,k),kp,isig3))


	enddo !kp

	enddo !k
	enddo !kb


      enddo !isig1
      
c -------------------------------------------------------------

	if (pdebug) then
	
	do k = 1,4
	
	do isig1 = -1,1,2
	do isig3 = -1,1,2
	
c	isig1 = -1
c	isig3 = -1
	
	m1 = 0d0
	m2 = 0d0

c   lower line for k=1 (upper for k=3):	
	if (k.eq.1.or.k.eq.3) then
	do i = 1,2
	   m1 = m1 + tampw(i,isig1,isig3,k)
	enddo
	do i = 5,6
	   m1 = m1 + tampw(i,isig1,isig3,k)
	enddo
        do i = 15,16
	  m1 = m1 + tampw(i,isig1,isig3,k)
        enddo
c   upper line for k=1 (lower for k=3):	
       do i = 24,27
          m1 = m1 + tampw(i,isig1,isig3,k)
       enddo	  	
       do i = 32,35
          m1 = m1 + tampw(i,isig1,isig3,k)
       enddo	  	
       do i = 50,53
          m1 = m1 + tampw(i,isig1,isig3,k)
       enddo	  	
       do i = 64,67
          m1 = m1 + tampw(i,isig1,isig3,k)
       enddo	  	
       do i = 72,75
          m1 = m1 + tampw(i,isig1,isig3,k)
       enddo	  	
       do i = 90,93
          m1 = m1 + tampw(i,isig1,isig3,k)
       enddo	  	
	
	elseif (k.eq.2.or.k.eq.4) then
	do i = 41,42 
	   m1 = m1 + tampw(i,isig1,isig3,k)
	enddo
	do i = 45,46
	   m1 = m1 + tampw(i,isig1,isig3,k)
	enddo
       do i = 55,56
	  m1 = m1 + tampw(i,isig1,isig3,k)
       enddo
       do i = 5,8
	  m1 = m1 + tampw(i,isig1,isig3,k)
       enddo
       do i = 13,16
	  m1 = m1 + tampw(i,isig1,isig3,k)
       enddo
       do i = 31,34
	  m1 = m1 + tampw(i,isig1,isig3,k)
       enddo
       do i = 64,67
	  m1 = m1 + tampw(i,isig1,isig3,k)
       enddo
       do i = 72,75
	  m1 = m1 + tampw(i,isig1,isig3,k)
       enddo
       do i = 90,93
	  m1 = m1 + tampw(i,isig1,isig3,k)
       enddo
       endif
	
	
c	if (.false.) then	
       print*,"comparison for pentagons from upper&lower line"
       print*,"with k=",k,", isig1,isig3=",isig1,isig3
       if (m1.ne.0d0) then
       if (k.le.2) then
       	print*,"m1,mat=",m1,mat(k,isig1,isig3,7)+mat(k,isig1,isig3,8)
       	print*,"m1/mat=",m1/(mat(k,isig1,isig3,7)+mat(k,isig1,isig3,8))
       	print*
       else
       	print*,"m1,mat=",m1,mat(k,isig3,isig1,7)+mat(k,isig3,isig1,8)
       	print*,"m1/mat=",m1/(mat(k,isig3,isig1,7)+mat(k,isig3,isig1,8))
       	print*
       endif
       endif
c       endif
       
	enddo !isig1
	enddo !isig3
	
	enddo !k
	
	endif      


c ------------------------------------------------------------------------------
c ------------------------------------------------------------------------------
     
c sum the graphs, square them and map them onto ansi etc.

c i = 1		VV
c i = 2,3	BV
c i = 4,5	WZZ
c i = 6,9	BB
c i = 7,8	P
   
      if (lpt) then
         xgc1 = real(icount1+1)/real(icount1-icb1*0.5+1) ! factor 1/2 for 2 helicities
         xgc2 = real(icount2+1)/real(icount2-icb2*0.5+1)
      endif
		      

      do k = 1,4
         res(k) = 0
         resv(k) = 0
         tree(k) = 0
         do isig1 = -1,1,2
            do isig3 = -1,1,2
               mm(k,isig1,isig3) = 0
               do i = 1,8
                  mm(k,isig1,isig3) = mm(k,isig1,isig3) + 
     1                                mat(k,isig1,isig3,i)
               enddo
	       
               res(k) = res(k) + dreal(mm(k,isig1,isig3))**2
     &                         + dimag(mm(k,isig1,isig3))**2
                if (lnlo) then
                  mv12(k,isig1,isig3) = 
     1         +  matv(k,isig1,isig3,2) + matv(k,isig1,isig3,4) + 
     2            matv(k,isig1,isig3,6) + matv(k,isig1,isig3,7) 
                  mv34(k,isig1,isig3) =
     1          + matv(k,isig1,isig3,3) + matv(k,isig1,isig3,5) +
     2            matv(k,isig1,isig3,9) + matv(k,isig1,isig3,8)

                  if (lpt) then   ! Pentagons are considered for "good PS points" only
                     if (lgc(js1)) then
                        mv12(k,isig1,isig3)=0
                     else
                        mv12(k,isig1,isig3)=mv12(k,isig1,isig3)*xgc1
                     endif
                  endif
                  if (lpt) then   ! Pentagon contributes are considered for "good PS points" only
                     if (lgc(js3)) then
                        mv34(k,isig1,isig3)=0
                     else
                        mv34(k,isig1,isig3)=mv34(k,isig1,isig3)*xgc2
                     endif
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
 
      if (ldebug) then
		print*,"uucs",ans1,res(1),res(1)/ans1
		print*,"ddcs",ans2,res(2),res(2)/ans2
		print*,"udcc",ans3,res(3),res(3)/ans3
		print*,"udss",ans4,res(4),res(4)/ans4
        	print*
      endif
c      if (ldebugv) 
c      read(*,"(a)") text
  
      if (lwarn) then
         resv(1) = ans1	!uucs
         resv(2) = ans2	!ddcs
         resv(3) = ans3	!udcc
         resv(4) = ans4	!udss
         do k = 1,4
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

      ans1 = res(1)	!uucs 
      ans2 = res(2)	!ddcs 
      ans3 = res(3)	!udcc 
      ans4 = res(4)	!udss 

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


